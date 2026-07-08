#ifndef IMRPHENOMD_NRT_EOS
#define IMRPHENOMD_NRT_EOS
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_spline.h>

#include <armadillo>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

#include "IMRPhenomD_NRT.h"
#include "util.h"

/**
 * @file Class that extends the IMRPhenomD_NRT waveform to sample directly on
 * equation of state (EOS) parameters.
 */

/* -------------------------------------------------------------------------- */
/*                       Interpolation Class Definition                       */
/* -------------------------------------------------------------------------- */

/**
 * @class Interpolation IMRPhenomD_NRT_EOS.h "include/gwat/IMRPhenomD_NRT_EOS.h"
 * @brief Defines a class to create an interpolated function y(x).
 */
class Interpolation {
 public:
  // Constructors and destructors
  Interpolation(gsl_interp_type* interp_type, std::vector<double> x,
                std::vector<double> y);
  Interpolation();
  ~Interpolation();

  // Function to initialize GSL splines for interpolation
  void initialize(gsl_interp_type* interp_type, std::vector<double> x,
                  std::vector<double> y);

  // Function to free GSL spline and accelerator memory
  void free();

  // Functions to evaluate points
  double yofx(double x);
  double dyofx(double x);

 protected:
  // GSL-type variables for interpolation
  gsl_interp_accel* acc;
  gsl_spline* spline;
  gsl_interp_type* type;

  // Size of data
  size_t size;
};

/* -------------------------------------------------------------------------- */
/*                     IMRPhenomD_NRT_EOS Class Definition                    */
/* -------------------------------------------------------------------------- */

/**
 * @class IMRPhenomD_NRT_EOS IMRPhenomD_NRT_EOS.h
 * "include/gwat/IMRPhenomD_NRT_EOS.h"
 * @brief Defines a class to construct a waveform with EoS parameters.
 */
template <class T>
class IMRPhenomD_NRT_EOS : public IMRPhenomD_NRT<T> {
 public:
  // Function to override IMRPhenomD_NRT construct_waveform
  virtual int construct_waveform(T* frequencies, int length,
                                 std::complex<T>* waveform,
                                 source_parameters<T>* params) override;

  // Function to calculate and update masses and tidal love numbers
  void get_m_love(source_parameters<T>* params);

 protected:
  // Function to update mass-dependent parameters
  void get_observable_params(source_parameters<T>* params);
};

/* -------------------------------------------------------------------------- */
/*                    EoS Constructor Base Class Definition                   */
/* -------------------------------------------------------------------------- */

/**
 * @struct EoSData IMRPhenomD_NRT_EOS.h "include/gwat/IMRPhenomD_NRT_EOS.h"
 * @brief Struct to store data about an EoS for the EoS-based waveform.
 */
struct EoSData {
  // Vectors
  std::vector<double> pressure;
  std::vector<double> epsilon;
  std::vector<double> nb;
  std::vector<double> cs2;

  // Central values to return
  double eps_c1;
  double eps_c2;
};

/**
 * @struct CrustEoS IMRPhenomD_NRT_EOS.h "include/gwat/IMRPhenomD_NRT_EOS.h"
 * @brief Struct to store data about a crust EoS table for the EoS-based
 * waveform.
 */
struct CrustEoS {
  // Vectors
  std::vector<double> pressure;
  std::vector<double> epsilon;
  std::vector<double> nb;
  std::vector<double> cs2;
};

/**
 * @class EOS_Constructor IMRPhenomD_NRT_EOS.h
 * "include/gwat/IMRPhenomD_NRT_EOS.h"
 * @brief Defines a base class to construct an equation of state.
 */
class EOS_Constructor {
 public:
  // Object to store built EoS information
  EoSData eos;

  // Static object to store the crust table
  static CrustEoS crust;

  // Constructor and destructor for inheritance, if necessary
  EOS_Constructor() {}
  virtual ~EOS_Constructor() {}

  // Function to load a file of crust table data
  void get_crust(std::string EoS_filepath, std::string pressure_header,
                 std::string epsilon_header, std::string nb_header);

  // Function to clear the vectors in the structs storing the EoS information
  void clear_crust();
  void clear_EOS();

  // Functions for unit conversion
  double convert_fm3_to_MeV(double x);
  double convert_nsat_to_MeV(double x);

 protected:
  // Function to create a cs2 column for the crust table using GSL spline
  // derivative calculations
  void convert_crust_to_cs2();

  // Function to convert a stored cs2(nb) EoS into pressure and epsilon using a
  // simple integration routine
  void convert_cs2_to_eos(std::size_t start_index);
};

/* -------------------------------------------------------------------------- */
/*               Bumpy EoS Constructor Derived Class Definition               */
/* -------------------------------------------------------------------------- */

/**
 * @class Bumpy_EOS_Constructor IMRPhenomD_NRT_EOS.h
 * "include/gwat/IMRPhenomD_NRT_EOS.h"
 * @brief Defines a class to construct an equation of state with a parabolic
 * bump in cs2 from given parameters.
 */
class Bumpy_EOS_Constructor : public EOS_Constructor {
 public:
  /**
   * @struct BumpyParams IMRPhenomD_NRT_EOS.h
   * "include/gwat/IMRPhenomD_NRT_EOS.h"
   * @brief Struct to store parameters describing a parabolic bump in cs2.
   */
  struct BumpyParams {
    double bump_magnitude;
    double bump_width;
    double bump_offset;
    double plat;
    double nbc;
    double n1;
    double n2;

    // Additional params
    double bump_start;
    double bump_end;
    double f1_n1;
    double max_nb;
  };

  BumpyParams eos_params;

  // Functions to store parameters
  void store_EOS_params(source_parameters<adouble>* params);
  void store_EOS_params(source_parameters<double>* params);
  virtual void get_additional_EOS_params();

  // Function to inject a parabolic bump into a crust EoS
  virtual void construct_EOS();

 protected:
  // Injection point of the bump
  std::size_t injection_index;

  // Functinos to calculate points in cs2 for the quadratic bump
  void build_cs2_one_quad_bump();
  double get_quadratic_bump_point(const double& nb, const double& f1_n1);
};

/* -------------------------------------------------------------------------- */
/*                    TOV/Tidal Integrator Class Definition                   */
/* -------------------------------------------------------------------------- */
/**
 * @class ObservablesIntegrator IMRPhenomD_NRT_EOS.h
 * "include/gwat/IMRPhenomD_NRT_EOS.h"
 * @brief Defines a class to integrate for observable values mass, radius, and
 * tidal deformability for a given EoS.
 *
 * @details Integrates in terms of enthalpy, according to 10.1086/171882. Uses
 * the classic Runge-Kutta method.
 */
class ObservablesIntegrator {
 public:
  // Function to initialize the integrator with enthalpy and central values
  ObservablesIntegrator(Interpolation& p_of_e, double ec_1, double ec_2);

  // Function to drive integration for the observable values
  void integrate_for_observables(arma::vec::fixed<3>& first_observables,
                                 arma::vec::fixed<3>& second_observables,
                                 bool& CurveIsNegative);

 protected:
  // EoS Interpolations
  Interpolation p_of_h;
  Interpolation e_of_h;

  // Central enthalpy values
  double hc_1;
  double hc_2;

  // Shifted enthalpy values
  double shift = 0.001;
  double hc_1_shift;
  double hc_2_shift;

  struct common_constants {
    double x1 = 4.0 * M_PI;
  };

  common_constants constant;

  // Functions to get p(h) and e(h)
  double calculate_dh_of_de(const double e, Interpolation& p_of_e);
  void integrate_for_eos_of_h(Interpolation& p_of_e, double ec_1, double ec_2);

  // Functions to evaluate integration for m(h) and r(h)
  void evaluate_ODE_at_point(double h, const arma::vec::fixed<3>& y,
                             arma::vec::fixed<3>& f);
  void rk4(double& t, arma::vec::fixed<3>& y, double h);

  // Functions to convert units
  double convert_dimensions(double value, std::string unit,
                            bool adimensionalize);
};

#endif