#ifndef IMRPHENOMD_NRT_EOS
#define IMRPHENOMD_NRT_EOS
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_spline.h>

#include <cmath>
#include <cstdlib>
#include <string>
#include <valarray>
#include <vector>

#include "IMRPhenomD_NRT.h"
#include "util.h"

using std::string;
using std::valarray;
using std::vector;

/*! \file
 */

/*! Class that extends the IMRPhenomD_NRT waveform to sample directly on
 * equation of state (EOS) parameters.
 *
 * THIS FILE IS DYSFUNCTIONAL RIGHT NOW - DO NOT USE
 */

/* -------------------------------------------------------------------------- */
/*                       Interpolation Class Definition                       */
/* -------------------------------------------------------------------------- */

/* -------- Defines a class to create an interpolated function y(x). -------- */
class Interpolation {
public:
  // Calls initialize()
  Interpolation(gsl_interp_type *type, vector<double> x, vector<double> y);
  Interpolation();
  // Calls free()
  ~Interpolation();

  // Initializes GSL splines for interpolation
  void initialize(gsl_interp_type *type, vector<double> x, vector<double> y);

  // Frees GSL spline and accelerator memory
  void free();

  // Methods to evaluate interpolated function at point x
  double yofx(double x);
  double dyofx(double x);

protected:
  // GSL-type variables for interpolation
  gsl_interp_accel *acc;
  gsl_spline *spline;
  gsl_interp_type *type;

  // Size of data to be interpolated
  size_t size;
};

/* -------------------------------------------------------------------------- */
/*                     IMRPhenomD_NRT_EOS Class Definition                    */
/* -------------------------------------------------------------------------- */

/* ------ Defines a class to construct a waveform with EoS parameters. ------ */
template <class T> class IMRPhenomD_NRT_EOS : public IMRPhenomD_NRT<T> {
public:
  // Overrides IMRPhenomD_NRT construct_waveform
  int construct_waveform(T *frequencies, int length, std::complex<T> *waveform,
                         source_parameters<T> *params) override;

  void get_m_love(source_parameters<T> &params);
};

/* -------------------------------------------------------------------------- */
/*                    EoS Constructor Base Class Definition                   */
/* -------------------------------------------------------------------------- */

/* --------------- Defines a *base class* to construct an EoS. -------------- */
class EOS_Constructor {
public:
  EOS_Constructor(string EOS_filepath);

  void get_EOS(Interpolation &p_of_e, double &central_epsilon_1,
               double &central_epsilon_2);

protected:
  struct EoS_data {
    // Vectors
    vector<double> pressure;
    vector<double> epsilon;
    vector<double> nb;
    vector<double> cs2;

    // Central values to return
    double eps_c_1;
    double eps_c_2;
  };

  EoS_data eos;

  void transpose_data_to_column_major(const vector<vector<double>> &row_major,
                                      vector<vector<double>> &column_major);

  // Conversion methods
  void convert_eos_to_cs2();
  void convert_cs2_to_eos();
  void convert_fm3_to_MeV(double &x);
  void convert_nsat_to_MeV(double &x);
};

/* -------------------------------------------------------------------------- */
/*               Bumpy EoS Constructor Derived Class Definition               */
/* -------------------------------------------------------------------------- */

/* ------------------ Implementation of constructor is WIP ------------------ */

/* ----- Defines a class to construct a bumpy EoS from given parameters. ---- */
class Bumpy_EOS_Constructor : public EOS_Constructor {
public:
  Bumpy_EOS_Constructor(string EOS_filepath) : EOS_Constructor(EOS_filepath){};

  // Injects bump into the EOS
  virtual void construct_EOS(double bump_mag, double bump_width,
                             double bump_offset, double plat, double nbc1,
                             double nbc2);

protected:
  struct Bumpy_Params {
    double bump_magnitude;
    double bump_width;
    double bump_offset;
    double plat;
    double nbc;
    double n1;
    double n2;
  };

  Bumpy_Params eos_params;

  // Methods to calculate points in cs2 for the quadratic bump
  void build_cs2_one_quad_bump();
  double get_quadratic_bump_point(const double &nb, const double &f1_n1);
};

/* -------------------------------------------------------------------------- */
/*                    TOV/Tidal Integrator Class Definition                   */
/* -------------------------------------------------------------------------- */

/* ------------------------------- WIP Section ------------------------------ */

class ObservablesIntegrator {
public:
  ObservablesIntegrator(double epsilon_c1, double epsilon_c2,
                        Interpolation eos);

  void integrate_for_observables(double &m1, double &m2, double &tidal1,
                                 double &tidal2, double step_size);

protected:
  Interpolation p_of_e;
  Interpolation p_of_h;
  Interpolation e_of_h;

  double eps_c1;
  double eps_c2;

  double h_max1;
  double h_max2;

  double calculate_dh_of_de(double epsilon);
  void integrate_for_eos_of_h();
  valarray<double> Observables_ODE(double h, const valarray<double> y);
  void rk4(double h, valarray<double> &y, double step_size);
};

#endif