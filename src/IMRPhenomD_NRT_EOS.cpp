/**
 * @file IMRPhenomD_NRT_EOS.cpp
 *
 * @brief File for the generation of equations of state with features in the
 * speed of sound squared as in arXiv:2106.03890 and related papers.
 *
 * @details Supported features: single parabolic bump ... [more to be possibly
 * added later]
 *
 * Supported crust types: SLy EOS
 *
 */

#include "IMRPhenomD_NRT_EOS.h"

#include <adolc/adouble.h>
#include <adolc/drivers/drivers.h>
#include <adolc/taping.h>
#include <gsl/gsl_randist.h>
#include <math.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "IMRPhenomD_NRT.h"
#include "io_util.h"
#include "util.h"

/* -------------------------------------------------------------------------- */
/*                         IMRPhenomD_NRT_EOS METHODS                         */
/* -------------------------------------------------------------------------- */

/* ----------------------------- Public Methods ----------------------------- */

/**
 * @brief Overrides the IMRPhenomD_NRT waveform constructor to add EoS
 * functionality
 *
 * @details Calculates the EoS using the specified parameters, then integrates
 * for the mass and tidal deformability using the central values. Updates
 * observable parameters that depend on the mass and tidal deformability, then
 * passes this information off to the base IMRPhenomD_NRT waveform constructor.
 *
 * @tparam T
 * @param frequencies
 * @param length
 * @param waveform
 * @param params
 * @return int
 */
template <class T>
int IMRPhenomD_NRT_EOS<T>::construct_waveform(T* frequencies, int length,
                                              std::complex<T>* waveform,
                                              source_parameters<T>* params) {
  get_m_love(params);

  int result = this->IMRPhenomD_NRT<T>::construct_waveform(frequencies, length,
                                                           waveform, params);

  return result;
}

/**
 * @brief Convert EOS parameters to neutron star masses and tidal
 * deformabilities.
 *
 * @details Given the equation of state (EOS) parameters (bump width, offset,
 * magnitude, and central baryon number densities), this function computes the
 * corresponding astrophysical properties (masses and tidal deformabilities)
 * for a neutron star binary.
 *
 * @param[in,out] params On output, mass1, mass2, tidal1, and tidal2 as well as
 * dependent parameters are updated.
 */
template <class T>
void IMRPhenomD_NRT_EOS<T>::get_m_love(source_parameters<T>* params) {
  // Define objects needed to store information
  std::string eos_filepath = "/opt/gw_analysis_tools/data/eos.csv";

  // Central epsilon values
  double ec1;
  double ec2;

  // Objects to store integration information
  arma::vec::fixed<3> observables1;
  arma::vec::fixed<3> observables2;

  // Constructs the bumpy EoS and stores that information in interpolated object
  Bumpy_EOS_Constructor* bumpy_eos = new Bumpy_EOS_Constructor;
  bumpy_eos->store_EOS_params(params);
  bumpy_eos->get_additional_EOS_params();
  bumpy_eos->construct_EOS();
  ec1 = bumpy_eos->eos.eps_c1;
  ec2 = bumpy_eos->eos.eps_c2;
  // Interpolated pressure in terms of epsilon
  Interpolation* p_of_e =
      new Interpolation((gsl_interp_type*)gsl_interp_steffen,
                        bumpy_eos->eos.epsilon, bumpy_eos->eos.pressure);
  delete bumpy_eos;
  bumpy_eos = 0;

  // Get observable values
  ObservablesIntegrator* integrator =
      new ObservablesIntegrator(*p_of_e, ec1, ec2);
  bool CurvatureCheck;

  integrator->integrate_for_observables(observables1, observables2,
                                        CurvatureCheck);

  delete integrator;
  integrator = 0;
  delete p_of_e;
  p_of_e = 0;

  // Check if curvature is valid
  if (CurvatureCheck) {
    std::cout << "This curvature... is bad. Sobs." << std::endl;
  }

  // TODO: (Kaitlyn) For debugging, remove when done
  /* std::cout << "Mass (1): " << observables1[1]
            << ", Mass (2): " << observables2[1]
            << ", Tidal (1): " << observables1[2]
            << ", Tidal (2): " << observables2[2] << std::endl; */

  // Store observables into params structure
  params->mass1 = observables1[1] * MSOL_SEC;
  params->mass2 = observables2[1] * MSOL_SEC;
  params->tidal1 = observables1[2];
  params->tidal2 = observables2[2];

  // TODO: (Kaitlyn) Perhaps replace this with the repopulate source parameters
  // function? The only issue with this is I only need the mass-dependent
  // parameters, the repopulate function calculates EVERYTHING.

  // Update dependent observable parameters
  get_observable_params(params);
}

/**
 * @brief Calculate mass-derived parameters for source_parameters.
 *
 * @details After integrating for the mass and tidal deformability, this
 * repopulates all of the parameters that are dependent on mass. Without this,
 * GWAT would think these values are NaN.
 *
 * @tparam T Double/adouble specification for source_parameters.
 * @param params GWAT's source_parameters object.
 */
template <class T>
void IMRPhenomD_NRT_EOS<T>::get_observable_params(
    source_parameters<T>* params) {
  params->q = params->mass2 / params->mass1;
  params->chirpmass = calculate_chirpmass(params->mass1, params->mass2);
  params->eta = calculate_eta(params->mass1, params->mass2);
  params->M = params->mass1 + params->mass2;
  params->chi_eff =
      (params->mass1 * (params->spin1z) + params->mass2 * (params->spin2z)) /
      (params->M);
  params->chi_pn =
      params->chi_eff - (38 * params->eta / 113) * (2 * params->chi_s);
  params->delta_mass = sqrt(1. - 4 * params->eta);

  if (params->sky_average) {
    params->A0 = sqrt(M_PI / 30) * params->chirpmass * params->chirpmass /
                 params->DL * pow(M_PI * params->chirpmass, -7. / 6);
  } else {
    params->A0 = sqrt(M_PI * 40. / 192.) * params->chirpmass *
                 params->chirpmass / params->DL *
                 pow(M_PI * params->chirpmass, -7. / 6);
  }
}

// Because GWAT hates me (DON'T REMOVE THIS THE COMPILER WILL FAIL TO FIND THE
// WAVEFORM CONSTRUCTOR)
template class IMRPhenomD_NRT_EOS<double>;
template class IMRPhenomD_NRT_EOS<adouble>;

/* -------------------------------------------------------------------------- */
/*                          EOS Contructor Functions                          */
/* -------------------------------------------------------------------------- */

/* ---------------------------- Public Functions ---------------------------- */

CrustEoS EOS_Constructor::crust;

/**
 * @brief Read-in and store an EoS table
 *
 * @details Assumes the EoS table is a CSV file with headers and contains
 * pressure, energy density, and baryon number density.
 *
 * @param EoS_filepath Path to file relative to execution directory
 * @param pressure_header Header for the pressure column in the CSV file
 * @param epsilon_header Header for the energy density column in the CSV file
 * @param nb_header Header for the baryon number density column in the CSV file
 */
void EOS_Constructor::get_crust(std::string EoS_filepath,
                                std::string pressure_header,
                                std::string epsilon_header,
                                std::string nb_header) {
  // Read in EOS table
  std::unordered_map<std::string, std::vector<double>> EoS_Table;
  read_csv(EoS_filepath, EoS_Table, ',', true);

  // Convert units to MeV
  std::transform(EoS_Table.at(nb_header).begin(), EoS_Table.at(nb_header).end(),
                 EoS_Table.at(nb_header).begin(),
                 std::bind(&EOS_Constructor::convert_fm3_to_MeV, this,
                           std::placeholders::_1));

  // Store table information
  crust.pressure = EoS_Table.at(pressure_header);
  crust.epsilon = EoS_Table.at(epsilon_header);
  crust.nb = EoS_Table.at(nb_header);
  convert_crust_to_cs2();
}

/**
 * @brief Function to clear the crust EoS vectors
 *
 */
void EOS_Constructor::clear_crust() {
  crust.pressure.clear();
  crust.epsilon.clear();
  crust.nb.clear();
  crust.cs2.clear();
}

/**
 * @brief Function to clear the built EoS vectors and central values
 *
 */
void EOS_Constructor::clear_EOS() {
  eos.pressure.clear();
  eos.epsilon.clear();
  eos.nb.clear();
  eos.cs2.clear();
  eos.eps_c1 = std::nan("");
  eos.eps_c2 = std::nan("");
}

/**
 * @brief Convert value from units of fm^-3 to MeV.
 *
 * @param[in] x Value in fm^-3.
 * @return Value converted to MeV.
 */
double EOS_Constructor::convert_fm3_to_MeV(double x) {
  return x = x * (197.3 * 197.3 * 197.3);
}

/**
 * @brief Convert value from units of nsat to MeV.
 *
 * @param[in] x Value in nsat.
 * @return Value converted to MeV.
 */
double EOS_Constructor::convert_nsat_to_MeV(double x) {
  double nsat = 0.16;
  x *= nsat;
  return convert_fm3_to_MeV(x);
}

/* --------------------------- Protected Functions -------------------------- */

/**
 * @brief Calculate the speed of sound squared for each entry in the provided
 * EoS table.
 *
 */
void EOS_Constructor::convert_crust_to_cs2() {
  Interpolation p_of_e((gsl_interp_type*)gsl_interp_steffen, crust.epsilon,
                       crust.pressure);

  for (double epsilon : crust.epsilon) {
    crust.cs2.push_back(p_of_e.dyofx(epsilon));
  }
}

/**
 * @brief Convert speed of sound squared values to pressure and energy density
 * values.
 *
 * @param start_index Starting index for the baryon number density
 */
void EOS_Constructor::convert_cs2_to_eos(std::size_t start_index) {
  // Get starting points for integration
  double nb = eos.nb.at(start_index - 1);
  double pressure = crust.pressure.at(start_index - 1);
  double epsilon = crust.epsilon.at(start_index - 1);

  // Copy existing crust data
  eos.pressure.insert(eos.pressure.begin(), crust.pressure.begin(),
                      crust.pressure.begin() + start_index);
  eos.epsilon.insert(eos.epsilon.begin(), crust.epsilon.begin(),
                     crust.epsilon.begin() + start_index);

  // Loop through all baryon number density values and perform integration
  for (int i = start_index; i < eos.nb.size(); i++) {
    double delta_nb = eos.nb[i] - nb;  // Calculate step size
    double delta_e =
        delta_nb * (epsilon + pressure) / nb;  // Calculate energy delta at step

    // Take steps
    nb = eos.nb[i];
    epsilon += delta_e;
    pressure += eos.cs2[i - 1] * delta_e;

    // Update pressure and energy density
    eos.pressure.push_back(pressure);
    eos.epsilon.push_back(epsilon);
  }
}

/* -------------------------------------------------------------------------- */
/*                       Bumpy EOS Constructor Functions                      */
/* -------------------------------------------------------------------------- */

/* ----------------------------- Public Methods ----------------------------- */

/**
 * @brief Converts GWAT's parameters object to a set of more useable double
 * parameters.
 *
 * @param params GWAT's source_parameters object.
 */
void Bumpy_EOS_Constructor::store_EOS_params(
    source_parameters<adouble>* params) {
  eos_params.bump_magnitude = params->bump_mag.value();
  eos_params.bump_width = convert_nsat_to_MeV(params->bump_width.value());
  eos_params.bump_offset = convert_nsat_to_MeV(params->bump_offset.value());
  eos_params.plat = params->plat.value();
  eos_params.n1 = convert_nsat_to_MeV(params->nbc1.value());
  eos_params.n2 = convert_nsat_to_MeV(params->nbc2.value());
}

/**
 * @brief Converts GWAT's parameters object to a set of more useable double
 * parameters.
 *
 * @param params GWAT's source_parameters object.
 */
void Bumpy_EOS_Constructor::store_EOS_params(
    source_parameters<double>* params) {
  eos_params.bump_magnitude = params->bump_mag;
  eos_params.bump_width = convert_nsat_to_MeV(params->bump_width);
  eos_params.bump_offset = convert_nsat_to_MeV(params->bump_offset);
  eos_params.plat = params->plat;
  eos_params.n1 = convert_nsat_to_MeV(params->nbc1);
  eos_params.n2 = convert_nsat_to_MeV(params->nbc2);
}

void Bumpy_EOS_Constructor::get_additional_EOS_params() {
  // Get additional helpful values to define the bump
  eos_params.bump_start = eos_params.bump_offset - (eos_params.bump_width / 2);
  eos_params.bump_end = eos_params.bump_offset + (eos_params.bump_width / 2);
  eos_params.max_nb = std::max(eos_params.n1, eos_params.n2);
  Interpolation cs2_of_nb((gsl_interp_type*)gsl_interp_steffen, crust.nb,
                          crust.cs2);
  eos_params.f1_n1 = cs2_of_nb.yofx(eos_params.bump_start);
}

/**
 * @brief Constructs and stores a GSL interpolation object for a bumpy EoS from
 * a given set of parameters.
 *
 */
void Bumpy_EOS_Constructor::construct_EOS() {
  // Check for valid inputs
  if (eos_params.max_nb > crust.nb.back()) {
    throw std::invalid_argument(
        "Encountered central baryon number density greater than the max value "
        "given in the EOS table. Aborting.");
  } else if (eos_params.bump_start <= convert_nsat_to_MeV(0.5)) {
    throw std::invalid_argument(
        "Encountered bump that starts below valid "
        "interpolation for SLy EoS tables. Aborting.");
  }

  // Build bump in the speed of sound
  build_cs2_one_quad_bump();
  // Convert speed of sound bump to pressure and epsilon values.
  convert_cs2_to_eos(injection_index);

  // Get central epsilon values
  Interpolation* e_of_nb = new Interpolation;

  e_of_nb->initialize((gsl_interp_type*)gsl_interp_steffen, eos.nb,
                      eos.epsilon);
  eos.eps_c1 = e_of_nb->yofx(eos_params.n1);
  eos.eps_c2 = e_of_nb->yofx(eos_params.n2);
  delete e_of_nb;
  e_of_nb = 0;
}

/* --------------------------- Protected Functions -------------------------- */

/**
 * @brief Interpolates the baryon number density and calculates the expected
 * cs^2 values for the given bump parameters.
 *
 * @details Uses the crust table points prior to the start of the bump, then
 * interpolates once the bump actually starts (step size of 100 in MeV), then
 * much less finely once the plateau starts (step size of 1000 in MeV). Calls
 * the quadratic bump calculator in the region of the bump, and simply uses
 * either crust or plateau values outside of this range.
 *
 */
void Bumpy_EOS_Constructor::build_cs2_one_quad_bump() {
  // Define parameters for interpolation in nb
  // How much to step by during the bump
  double nb_step_bump = eos_params.bump_width / 1e3;
  // How much to step by during the plateau
  double nb_step_plat = 10 * nb_step_bump;
  // Define an ending value slightly greater than the max central nb value, to
  // prevent interpolation errors
  double nb_end = eos_params.max_nb + nb_step_plat;
  double nb;

  // Simply use the crust prior to the start of the bump
  for (int i = 0; i < crust.nb.size(); i++) {
    nb = crust.nb.at(i);
    if ((nb > nb_end) || (nb >= eos_params.bump_start)) {
      injection_index = i - 1;
      nb = crust.nb.at(i - 1) + nb_step_bump;
      break;
    } else {
      eos.cs2.push_back(crust.cs2.at(i));
      eos.nb.push_back(nb);
    }
  }

  // Interpolate up to the bump start (this is a separate loop to prevent
  // creating an interpolation object when it's not necessary)
  if (nb < eos_params.bump_start) {
    // Get central epsilon values
    Interpolation* cs2_of_nb = new Interpolation;
    cs2_of_nb->initialize((gsl_interp_type*)gsl_interp_steffen, crust.nb,
                          crust.cs2);
    while ((nb < eos_params.bump_start) && (nb <= nb_end)) {
      double cs2_val = cs2_of_nb->yofx(nb);
      eos.cs2.push_back(cs2_val);
      eos.nb.push_back(nb);
      nb += nb_step_bump;
    }
    delete cs2_of_nb;
  }

  // Inject bump if maximum value has not yet been reached
  while (nb <= nb_end) {
    // Inject bump values between bump start and end
    if ((eos_params.bump_start <= nb) && (nb <= eos_params.bump_end)) {
      eos.cs2.push_back(get_quadratic_bump_point(nb, eos_params.f1_n1));
      eos.nb.push_back(nb);
      nb += nb_step_bump;
    }
    // Use plateau values after bump
    else {
      eos.cs2.push_back(eos_params.plat);
      eos.nb.push_back(nb);
      nb += nb_step_plat;
    }
  }
}

/**
 * @brief Calculates cs^2 value from point in baryon number density for given
 * EoS parameters.
 *
 * @details Defines a quadratic parabola to calculate cs2 points. The quadratic
 * is piecewise due to 2 roots resulting from whether the bump magnitude is
 * greater or less than the cs2 value at the start of the bump. For the purposes
 * of defining the function, it is useful to consider that the bump_start
 * parameter is "n1", the bump_end is "n2", and thus the cs2 value at n1 is
 * "f1_n1" and the plateau is then equivalent to the cs2 value at n2, "f3_n2."
 * Then, the quadratic is found by simply solving f1_n1 = a1 + a2 n_1 + a3
 * n_1^2 and f3_n2 = a1 + a2 n_2 + a3 n_2^2 (with the condition that the bump
 * magnitude is the peak of the parabola).
 *
 * @param nb Baryon number density at which to evaluate
 * @param f1_n1 Transition point of the bump
 * @return double
 */
double Bumpy_EOS_Constructor::get_quadratic_bump_point(const double& nb,
                                                       const double& f1_n1) {
  // Define helpful constants
  const double x1 = (nb - eos_params.bump_start);
  const double x1_2 = x1 * x1;

  const double x2 = (nb - eos_params.bump_end);
  const double x2_2 = x2 * x2;

  const double x3 = (eos_params.bump_start - eos_params.bump_end);
  const double x3_2 = x3 * x3;

  // Define unchanged terms in the function
  const double term1 = eos_params.plat * x1_2;
  const double term2 = f1_n1 * x2_2;
  const double term3 = 2 * eos_params.bump_magnitude * x1 * x2;
  const double unchanged_terms = term1 + term2 - term3;

  // Define the function roots (will change sign depending on conditions)
  const double func_root = 2 * x1 * x2 *
                           sqrt((eos_params.bump_magnitude - f1_n1) *
                                (eos_params.bump_magnitude - eos_params.plat));

  // Calculate value of piecewise function
  double return_value;
  if (f1_n1 < eos_params.bump_magnitude) {
    return_value = (1 / x3_2) * (unchanged_terms - func_root);
  } else {
    return_value = (1 / x3_2) * (unchanged_terms + func_root);
  }

  return return_value;
}

/* -------------------------------------------------------------------------- */
/*                       ObservablesIntegrator Functions                      */
/* -------------------------------------------------------------------------- */

/* ---------------------------- Public Functions ---------------------------- */

/**
 * @brief Construct a new Observables Integrator object, initialize the EoS in
 terms of enthalpy
 *
 * @details Takes in the current EoS (p(e)) and uses it to find e(h) and p(h),
 * where h is the pseudo-enthalpy. Formulation is based on doi:10.1086/171882.
 * It also takes in the central epsilon of two stars and converts it to central
 * pseudo-enthalpy.
 *
 * @param p_of_e Input equation of state
 * @param ec_1 Central epsilon of first star
 * @param ec_2 Central epsilon of second star
 */
ObservablesIntegrator::ObservablesIntegrator(Interpolation& p_of_e, double ec_1,
                                             double ec_2) {
  integrate_for_eos_of_h(p_of_e, ec_1, ec_2);
}

/**
 * @brief Integrates for the observable mass and radius based on the stored EoS.
 *
 * @details Sets some small starting radius (R_epsilon) and calculates the
 * corrected starting pseudo-enthalpy, h_c - h_epsilon (h_epsilon from
 * R_epsilon). Then, uses a classic fixed-step Runge-Kutta (rk4) routine to
 * integrate the TOV equations in terms of pseudo-enthalpy from the corrected
 * starting enthalpy to 0 (the surface of the star). Default number of steps for
 * integration is 10,000.
 *
 * @param first_observables Armadillo vector to store the first set of
 * calculated neutron star observables (Radius, Mass, Tidal).
 * @param second_observables Armadillo vector to store the second set of
 * calculated neutron star observables (Radius, Mass, Tidal).
 * @param CurveIsNegative Boolean to store if the slope of the m(epsilon) curve
 * is negative.
 */
void ObservablesIntegrator::integrate_for_observables(
    arma::vec::fixed<3>& first_observables,
    arma::vec::fixed<3>& second_observables, bool& CurveIsNegative) {
  // Define starting radius
  const double R_start = convert_dimensions(4e-4, "km", true);
  const double R_start_2 = R_start * R_start;

  // Common terms
  const double x1 = (2.0 / 3.0) * M_PI * R_start_2;
  const double x2 = 2.0 * R_start * x1;

  // Define starting enthalpy values
  double h1 = hc_1 - x1 * (e_of_h.yofx(hc_1) + 3.0 * p_of_h.yofx(hc_1));
  double h2 = hc_2 - x1 * (e_of_h.yofx(hc_2) + 3.0 * p_of_h.yofx(hc_2));

  // Define starting mass values
  const double M_start_1 = x2 * e_of_h.yofx(h1);
  const double M_start_2 = x2 * e_of_h.yofx(h2);

  // Define starting y values
  const double Y_start = 2.0;

  // Define enthalpy shift values
  double h1_shift = hc_1_shift - x1 * (e_of_h.yofx(hc_1_shift) +
                                       3.0 * p_of_h.yofx(hc_1_shift));
  double h2_shift = hc_2_shift - x1 * (e_of_h.yofx(hc_2_shift) +
                                       3.0 * p_of_h.yofx(hc_2_shift));

  // Define mass shift values
  const double M_shift_1 = x2 * e_of_h.yofx(h1_shift);
  const double M_shift_2 = x2 * e_of_h.yofx(h2_shift);

  // Define step sizes
  double N_steps = 5e4;
  double h1_step = -h1 / N_steps;
  double h2_step = -h2 / N_steps;

  // TODO: (Kaitlyn) Look into replacing rk4 routine with GSL or other
  // integration functionality.

  /* -------------------- Integrate for First Observables ------------------- */

  // Set starting values
  first_observables = {R_start, M_start_1, Y_start};

  // Run until just before h = 0
  while (h1 + h1_step > 0) {
    rk4(h1, first_observables, h1_step);
  }
  // Step exactly to h=0
  rk4(h1, first_observables, (0 - h1));

  // Calculate the dimensionless tidal love number
  double Y_1 = first_observables(2);

  // Useful definitions for faster calculations
  double C = first_observables(1) / first_observables(0);
  double C2 = C * C;
  double C3 = C * C2;
  double C4 = C2 * C2;
  double C5 = C * C4;
  double k_1 = (1.0 - 2.0 * C) * (1.0 - 2.0 * C);

  // Dimensionless tidal apsidal constant k2
  double k2_1 = (8.0 / 5.0) * k_1 * C5 * (2.0 * C * (Y_1 - 1.0) - Y_1 + 2.0) *
                (1.0 / (2.0 * C *
                            (4 * (Y_1 + 1.0) * C4 + (6.0 * Y_1 - 4.0) * C3 +
                             (26.0 - 22.0 * Y_1) * C2 +
                             3.0 * (5.0 * Y_1 - 8.0) * C - 3.0 * Y_1 + 6.0) +
                        3.0 * k_1 * (2 * C * (Y_1 - 1.0) - Y_1 + 2.0) *
                            log(1.0 - 2.0 * C)));

  // Store the dimensionless tidal love number (λ̄)
  first_observables(2) = (2.0 / 3.0) * k2_1 * (1 / C5);

  // Dimensionalize the radius
  first_observables(0) = convert_dimensions(first_observables(0), "km", false);

  /* ------------------- Integrate for Second Observables ------------------- */

  // Set starting values
  second_observables = {R_start, M_start_2, Y_start};

  // Run until just before h = 0
  while (h2 + h2_step > 0) {
    rk4(h2, second_observables, h2_step);
  }
  // Step exactly to h=0
  rk4(h2, second_observables, (0 - h2));

  // Calculate the dimensionless tidal love number
  double Y_2 = second_observables(2);

  // Useful definitions for faster calculations
  C = second_observables(1) / second_observables(0);
  C2 = C * C;
  C3 = C * C2;
  C4 = C2 * C2;
  C5 = C * C4;
  double k_2 = (1.0 - 2.0 * C) * (1.0 - 2.0 * C);

  // Dimensionless tidal apsidal constant k2
  double k2_2 = (8.0 / 5.0) * k_2 * C5 * (2.0 * C * (Y_2 - 1.0) - Y_2 + 2.0) *
                (1.0 / (2.0 * C *
                            (4 * (Y_2 + 1.0) * C4 + (6.0 * Y_2 - 4.0) * C3 +
                             (26.0 - 22.0 * Y_2) * C2 +
                             3.0 * (5.0 * Y_2 - 8.0) * C - 3.0 * Y_2 + 6.0) +
                        3.0 * k_2 * (2 * C * (Y_2 - 1.0) - Y_2 + 2.0) *
                            log(1.0 - 2.0 * C)));

  // Store the dimensionless tidal love number (λ̄)
  second_observables(2) = (2.0 / 3.0) * k2_2 * (1 / C5);

  // Dimensionalize the radius
  second_observables(0) =
      convert_dimensions(second_observables(0), "km", false);

  /* ----------------------- Integrate to Check Slope ----------------------- */

  arma::vec::fixed<3> first_shifts;
  arma::vec::fixed<3> second_shifts;

  // Set starting values
  first_shifts = {R_start, M_shift_1, Y_start};
  second_shifts = {R_start, M_shift_2, Y_start};

  // Run until just before h = 0
  while (h1_shift + h1_step > 0) {
    rk4(h1_shift, first_shifts, h1_step);
  }
  // Step exactly to h=0
  rk4(h1_shift, first_shifts, (0 - h1_shift));

  // Run until just before h = 0
  while (h2_shift + h2_step > 0) {
    rk4(h2_shift, second_shifts, h2_step);
  }
  // Step exactly to h=0
  rk4(h2_shift, second_shifts, (0 - h2_shift));

  // Get first-order differences to find the derivative
  double slope1 = (first_observables[1] - first_shifts[1]) / shift;
  double slope2 = (second_observables[1] - second_shifts[1]) / shift;

  if ((slope1 < 0) || (slope2 < 0)) {
    CurveIsNegative = true;
  } else {
    CurveIsNegative = false;
  }
}

/* --------------------------- Protected Functions -------------------------- */

/**
 * @brief Finds the derivative of enthalpy in terms of energy density for a
 * given point.
 *
 * @param e Energy density point at which to evaluate.
 * @param p_of_e GSL interpolation storing pressure as a function of energy
 * density.
 * @return double
 */
double ObservablesIntegrator::calculate_dh_of_de(const double e,
                                                 Interpolation& p_of_e) {
  double dhde = p_of_e.dyofx(e) / (e + p_of_e.yofx(e));
  return dhde;
}

/**
 * @brief Finds the equation of state (pressure/epsilon) in terms of enthalpy.
 *
 * @details Interpolates epsilon and calculates dh/de at each point to find the
 * enthalpy. Uses this to build vectors of epsilon, pressure, and enthalpy at
 * each point, then creates GSL interpolation objects for p(h) and e(h).
 *
 * @param p_of_e A GSL interpolation object storing pressure in terms of
 * enthalpy.
 * @param ec_1 Central energy density for the first star.
 * @param ec_2 Central energy density for the second star.
 */
void ObservablesIntegrator::integrate_for_eos_of_h(Interpolation& p_of_e,
                                                   double ec_1, double ec_2) {
  // Define vectors to store values for interpolation
  std::vector<double> epsilon;
  std::vector<double> pressure;
  std::vector<double> enthalpy;

  // Set starting values
  double e = 0;
  double p = 0;
  double h = 0;

  // Add starting values to vectors
  epsilon.push_back(e);
  pressure.push_back(p);
  enthalpy.push_back(h);

  // Set step value and end value for integration
  const double e_stop = std::max(ec_1, ec_2);
  double delta_e = 0.1;

  // Integrate until e_stop
  while (e < e_stop) {
    if (e + delta_e > e_stop) {
      delta_e = e_stop - e;
      // Exit if delta_e is too small to proceed
      if (delta_e < 1e-6) {
        e = e_stop;
        break;
      }
    }

    // Calculate needed step in enthalpy
    double delta_h = delta_e * calculate_dh_of_de(e + delta_e, p_of_e);

    // Take a step in the EoS
    e += delta_e;
    p = p_of_e.yofx(e);
    h += delta_h;

    // Add *adimensionalized* values to vectors
    epsilon.push_back(convert_dimensions(e, "MeV/fm^3", true));
    pressure.push_back(convert_dimensions(p, "MeV/fm^3", true));
    enthalpy.push_back(h);
  }

  // Define interpolations
  e_of_h.initialize((gsl_interp_type*)gsl_interp_steffen, enthalpy, epsilon);
  p_of_h.initialize((gsl_interp_type*)gsl_interp_steffen, enthalpy, pressure);

  // Get central enthalpy values
  Interpolation h_of_e((gsl_interp_type*)gsl_interp_steffen, epsilon, enthalpy);
  hc_1 = h_of_e.yofx(convert_dimensions(ec_1, "MeV/fm^3", true));
  hc_2 = h_of_e.yofx(convert_dimensions(ec_2, "MeV/fm^3", true));

  // Grab shifts for evaluating M(epsilon) curve
  hc_1_shift = hc_1 - shift;
  hc_2_shift = hc_2 - shift;
}

/**
 * @brief Evaluates the TOV equations at a single point of pseudo-enthalpy
 *
 * @param h Enthalpy at which to evaluate
 * @param y Equations to solve for (given in order of R(h), M(h), and Y(h))
 * @param f Object to store the results of calculating dR/dh, dM/dh, and dY/dh
 */
void ObservablesIntegrator::evaluate_ODE_at_point(double h,
                                                  const arma::vec::fixed<3>& y,
                                                  arma::vec::fixed<3>& f) {
  // Grab values of pressure and epsilon at the point
  const double p = p_of_h.yofx(h);
  const double e = e_of_h.yofx(h);
  const double dedh = e_of_h.dyofx(h);

  // Grab R and M at this point (also calculate R^2 and R^3, which is more
  // efficient to do this way than using pow())
  const double R = y(0);
  const double M = y(1);
  const double Y = y(2);

  // Common/helpful terms for derivative calculations
  const double R_inv = 1.0 / R;
  const double R_2 = R * R;
  const double R_2_pi = M_PI * R_2;
  const double R_3 = R * R_2;
  const double R_4 = R * R_3;
  const double R_3_pi = constant.x1 * R_3;

  const double M_2 = M * M;

  const double R_M_diff = (2.0 * M) - R;
  const double R_M_diff_2 = R_M_diff * R_M_diff;

  const double Y_2 = Y * Y;

  const double e_x_5 = 5 * e;

  // Calculate dR/dh, dM/dh, and dY/dh.
  double drdh = R * R_M_diff / (M + (constant.x1 * R_3 * p));
  double dmdh = constant.x1 * R_2 * e * drdh;
  double dydh =
      (-Y_2 * R_inv + Y * (4 * R_2_pi * (-e + p) + 1) / R_M_diff +
       (dedh * (1.0 / drdh)) * R_3_pi / (M + p * R_3_pi) +
       R_inv *
           (4 * M_2 + 4 * M * R * (2 * R_2_pi * (13 * p + e_x_5) - 3) +
            R_4 * constant.x1 * (p * (16 * p * R_2_pi - 9) - e_x_5) + 6 * R_2) /
           R_M_diff_2) *
      drdh;

  // Load calculated values into pass-in array
  f = {drdh, dmdh, dydh};
}

/**
 * @brief Simple classic fixed-step Runge-Kutta method for integration
 *
 * @param t Independent variable which functions are being taken with respect to
 * @param y Functions given in terms of the independent variable t
 * @param h Step size for the rk4 routine
 */
void ObservablesIntegrator::rk4(double& t, arma::vec::fixed<3>& y, double h) {
  // Calculate step size fractions
  const double h_6 = h / 6.0;
  const double h_2 = h / 2.0;

  // Define arrays to store rk4 term results in
  arma::vec::fixed<3> k1;
  arma::vec::fixed<3> k2;
  arma::vec::fixed<3> k3;
  arma::vec::fixed<3> k4;

  // Calculate rk4 terms
  evaluate_ODE_at_point(t, y, k1);
  evaluate_ODE_at_point(t + h_2, y + k1 * h_2, k2);
  evaluate_ODE_at_point(t + h_2, y + k2 * h_2, k3);
  evaluate_ODE_at_point(t + h, y + k3 * h, k4);

  // Update function values and take a step
  y += h_6 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
  t += h;
}

/**
 * @brief Converts dimensions for the integrator to solar masses (or vice
 * versa).
 *
 * @param value Double to convert.
 * @param unit Unit to convert ("km," "MeV/fm^3")
 * @param adimensionalize Boolean determining whether to convert to solar masses
 * (True) or convert from solar masses (False).
 * @return The converted double.
 */
double ObservablesIntegrator::convert_dimensions(double value, std::string unit,
                                                 bool adimensionalize) {
  // Define constants
  double eV = 1.602176634e-19;
  double MeV = 1e6 * eV;
  double fm = 1e-15;
  double km = 1e3;

  // Calculate and store the conversion factor
  double conversion_factor;

  if (unit == "km") {
    conversion_factor = km * (c * c / GWAT_G_SI) * (1 / GWAT_MSUN_SI);
  } else if (unit == "MeV/fm^3") {
    conversion_factor = (MeV / (fm * fm * fm)) *
                        sqrt(pow(GWAT_G_SI, 6) / pow(c, 16)) *
                        (GWAT_MSUN_SI * GWAT_MSUN_SI);
  } else {
    std::cerr << "Unsupported unit type passed to "
                 "ObservablesIntegrator::adimensionalize.\n";
    conversion_factor = 1;
  }

  // Store converted value
  double converted_value;

  if (adimensionalize) {
    converted_value = value * conversion_factor;
  } else {
    converted_value = value / conversion_factor;
  }

  return converted_value;
}

/* -------------------------------------------------------------------------- */
/*                        Interpolation Class Functions                       */
/* -------------------------------------------------------------------------- */

/* ---------------------------- Public Functions ---------------------------- */

// Constructors

Interpolation::Interpolation(gsl_interp_type* interp_type,
                             std::vector<double> x, std::vector<double> y) {
  initialize(interp_type, x, y);
}

Interpolation::Interpolation() {}

// Destructor

Interpolation::~Interpolation() { free(); }

// Function value calculations

/**
 * @brief Function to calculate interpolated y value for a given x using GSL
 * spline
 *
 * @param x
 * @return double
 */
double Interpolation::yofx(double x) {
  double result;
  if (x >= 0.0) {
    result = gsl_spline_eval(spline, x, acc);
  } else {
    result = 0.0;
  }
  return result;
}

/**
 * @brief Function to calculate dy/dx of intferpolated data using GSL spline
 *
 * @param x
 * @return double
 */
double Interpolation::dyofx(double x) {
  double result;
  if (x > 0.0) {
    result = gsl_spline_eval_deriv(spline, x, acc);
  } else {
    result = 0.0;
  }
  return result;
}

/**
 * @brief Function to initialize the GSL splines for interpolation
 *
 * @param interp_type
 * @param x
 * @param y
 */
void Interpolation::initialize(gsl_interp_type* interp_type,
                               std::vector<double> x, std::vector<double> y) {
  type = interp_type;
  // Size of the independent variable vector
  size = x.size();

  // Allocating GSL interpolation accelerator
  acc = gsl_interp_accel_alloc();
  if (!acc) {
    std::cerr << "Error: Failed to allocate memory for accelerator."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Allocating GSL spline of specified type and size
  spline = gsl_spline_alloc(type, size);

  // Initializing GSL spline with given data points (x, y) and size
  gsl_spline_init(spline, x.data(), y.data(), size);

  int status = gsl_spline_init(spline, x.data(), y.data(), size);
  if (status != GSL_SUCCESS) {
    std::cerr << "gsl_spline_init failed: " << gsl_strerror(status) << "\n";
    std::exit(EXIT_FAILURE);
  }
}

/**
 * @brief Function to release memory of spline and accelerator
 *
 */
void Interpolation::free() {
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
}

/* -------------------------------------------------------------------------- */