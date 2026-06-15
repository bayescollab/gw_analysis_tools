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
 *
 */

/* -------------------------------------------------------------------------- */
/*                       Interpolation Class Definition                       */
/* -------------------------------------------------------------------------- */

/* -------- Defines a class to create an interpolated function y(x). -------- */
class Interpolation {
 public:
  // Constructors and destructors
  Interpolation(gsl_interp_type* interp_type, std::vector<double> x,
                std::vector<double> y);
  Interpolation();
  ~Interpolation();

  // Initializes GSL splines for interpolation
  void initialize(gsl_interp_type* interp_type, std::vector<double> x,
                  std::vector<double> y);

  // Frees GSL spline and accelerator memory
  void free();

  // Evaluation functions
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

/* ------ Defines a class to construct a waveform with EoS parameters. ------ */
template <class T>
class IMRPhenomD_NRT_EOS : public IMRPhenomD_NRT<T> {
 public:
  // Overrides IMRPhenomD_NRT construct_waveform
  virtual int construct_waveform(T* frequencies, int length,
                                 std::complex<T>* waveform,
                                 source_parameters<T>* params) override;

  void get_m_love(source_parameters<T>* params);

 protected:
  void get_observable_params(source_parameters<T>* params);
};

/* -------------------------------------------------------------------------- */
/*                    EoS Constructor Base Class Definition                   */
/* -------------------------------------------------------------------------- */

// TODO: Format this properly.

struct EoS_data {
  // Vectors
  std::vector<double> pressure;
  std::vector<double> epsilon;
  std::vector<double> nb;
  std::vector<double> cs2;

  // Central values to return
  double eps_c1;
  double eps_c2;
};

struct crust_EoS {
  // Vectors
  std::vector<double> pressure;
  std::vector<double> epsilon;
  std::vector<double> nb;
  std::vector<double> cs2;
};

/* --------------- Defines a *base class* to construct an EoS. -------------- */
class EOS_Constructor {
 public:
  // Constructor and destructor for inheritance, if necessary
  EOS_Constructor() {}
  ~EOS_Constructor() {}

  // Object to store built EoS information
  EoS_data eos;

  // Static object to store the crust
  // For new users: Static objects will share this data across *all* instances
  // of the class. This memory is more-or-less global, but this is necessary for
  // the EoS classes as the crust table must be read and accessed *every time a
  // point is sampled*.
  static crust_EoS crust;

  void get_crust(std::string EoS_filepath, std::string pressure_header,
                 std::string epsilon_header, std::string nb_header);
  void clear_crust();
  void clear_EOS();

  // Conversion methods
  double convert_fm3_to_MeV(double x);
  double convert_nsat_to_MeV(double x);

 protected:
  void convert_crust_to_cs2();
  void convert_cs2_to_eos(std::size_t start_index);
};

/* -------------------------------------------------------------------------- */
/*               Bumpy EoS Constructor Derived Class Definition               */
/* -------------------------------------------------------------------------- */

/* ------------------ Implementation of constructor is WIP ------------------ */

/* ----- Defines a class to construct a bumpy EoS from given parameters. ---- */
class Bumpy_EOS_Constructor : public EOS_Constructor {
 public:
  struct Bumpy_Params {
    double bump_magnitude;
    double bump_width;
    double bump_offset;
    double plat;
    double nbc;
    double n1;
    double n2;

    double bump_start;
    double bump_end;
    double f1_n1;
    double max_nb;
  };

  Bumpy_Params eos_params;

  // Methods to store parameters
  void store_EOS_params(source_parameters<adouble>* params);
  void store_EOS_params(source_parameters<double>* params);
  virtual void get_additional_EOS_params();

  // Injects bump into the EOS
  virtual void construct_EOS();

 protected:
  std::size_t injection_index;

  // Methods to calculate points in cs2 for the quadratic bump
  void build_cs2_one_quad_bump();
  double get_quadratic_bump_point(const double& nb, const double& f1_n1);
};

/* -------------------------------------------------------------------------- */
/*                    TOV/Tidal Integrator Class Definition                   */
/* -------------------------------------------------------------------------- */

class ObservablesIntegrator {
 public:
  ObservablesIntegrator(Interpolation& p_of_e, double ec_1, double ec_2);

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