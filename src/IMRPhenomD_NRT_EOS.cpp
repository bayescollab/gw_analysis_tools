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
 * @brief Convert EOS parameters to neutron star masses and tidal
 * deformabilities.
 *
 * @details Given the equation of state (EOS) parameters (bump width, offset,
 * magnitude, and central baryon number densities), this function computes the
 * corresponding astrophysical properties (masses and tidal deformabilities)
 * for a neutron star binary.
 *
 * @param[in,out] parameters On output, mass1, mass2, tidal1, and tidal2 as well
 * as dependent parameters are updated.
 */
template <class T>
void IMRPhenomD_NRT_EOS<T>::get_m_love(gen_params_base<T>* parameters) {
  // Central epsilon values
  double ec1;
  double ec2;

  // Objects to store integration information
  Eigen::Vector3d observables1;
  Eigen::Vector3d observables2;

  // Constructs the bumpy EoS and stores that information in interpolated object
  Bumpy_EOS_Constructor* bumpy_eos = new Bumpy_EOS_Constructor;
  bumpy_eos->store_EOS_params(parameters);
  bumpy_eos->get_additional_EOS_params();
  bumpy_eos->construct_EOS();
  ec1 = bumpy_eos->eos.eps_c1;
  ec2 = bumpy_eos->eos.eps_c2;
  // Interpolated pressure in terms of epsilon
  Interpolation* p_of_e =
      new Interpolation(static_cast<const gsl_interp_type*>(gsl_interp_steffen),
                        bumpy_eos->eos.epsilon, bumpy_eos->eos.pressure);
  delete bumpy_eos;

  // Get observable values
  ObservablesIntegrator* integrator =
      new ObservablesIntegrator(*p_of_e, ec1, ec2);
  bool CurvatureCheck;

  integrator->integrate_for_observables(observables1, observables2,
                                        CurvatureCheck);

  delete integrator;
  delete p_of_e;

  /* // Check if curvature is valid
  if (CurvatureCheck) {
    std::cout << "This curvature... is bad. Sobs." << std::endl;
  } */

  // Store observables into parameters structure
  parameters->mass1 = observables1[1];
  parameters->mass2 = observables2[1];
  parameters->tidal1 = observables1[2];
  parameters->tidal2 = observables2[2];
}

// Because GWAT hates me (DON'T REMOVE THIS THE COMPILER WILL FAIL TO FIND THE
// WAVEFORM CONSTRUCTOR)
template class IMRPhenomD_NRT_EOS<double>;
template class IMRPhenomD_NRT_EOS<adouble>;

/* -------------------------------------------------------------------------- */
/*                          EOS Contructor Functions                          */
/* -------------------------------------------------------------------------- */

/* ---------------------------- Public Functions ---------------------------- */

// Static objects of the read-in crust EoS. This is to prevent unnecessary file
// re-reads and re-calculations.
CrustEoS EOS_Constructor::crust;
Interpolation EOS_Constructor::crust_cs2_of_nb;

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
  // Store interpolation object for the speed of sound squared
  crust_cs2_of_nb.initialize(
      static_cast<const gsl_interp_type*>(gsl_interp_steffen), crust.nb,
      crust.cs2);
}

/**
 * @brief Function to clear the crust EoS vectors and interpolation
 */
void EOS_Constructor::clear_crust() {
  crust.pressure.clear();
  crust.epsilon.clear();
  crust.nb.clear();
  crust.cs2.clear();
  crust_cs2_of_nb.free();
}

/**
 * @brief Function to clear the built EoS vectors and central values
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
 */
void EOS_Constructor::convert_crust_to_cs2() {
  Interpolation p_of_e(static_cast<const gsl_interp_type*>(gsl_interp_steffen),
                       crust.epsilon, crust.pressure);

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
void Bumpy_EOS_Constructor::store_EOS_params(gen_params_base<adouble>* params) {
  // Bump parameters
  eos_params.bump_magnitude = params->bump_mag.value();
  eos_params.bump_start = convert_nsat_to_MeV(params->bump_start.value());
  eos_params.bump_width = convert_nsat_to_MeV(params->bump_width.value());
  eos_params.bump_end = eos_params.bump_start + eos_params.bump_width;
  // Plateau
  eos_params.plat = params->plat.value();
  // Baryon number density parameters
  eos_params.nbc1 = convert_nsat_to_MeV(params->nbc1.value());
  eos_params.nbc2 = convert_nsat_to_MeV(params->nbc2.value());
  eos_params.max_nb = std::max(eos_params.nbc1, eos_params.nbc2);
}

/**
 * @brief Converts GWAT's parameters object to a set of more useable double
 * parameters.
 *
 * @param params GWAT's source_parameters object.
 */
void Bumpy_EOS_Constructor::store_EOS_params(gen_params_base<double>* params) {
  // Bump parameters
  eos_params.bump_magnitude = params->bump_mag;
  eos_params.bump_start = convert_nsat_to_MeV(params->bump_start);
  eos_params.bump_width = convert_nsat_to_MeV(params->bump_width);
  eos_params.bump_end = eos_params.bump_start + eos_params.bump_width;
  // Plateau
  eos_params.plat = params->plat;
  // Baryon number density parameters
  eos_params.nbc1 = convert_nsat_to_MeV(params->nbc1);
  eos_params.nbc2 = convert_nsat_to_MeV(params->nbc2);
  eos_params.max_nb = std::max(eos_params.nbc1, eos_params.nbc2);
}

/**
 * @brief Creates an interpolation of cs^2 and nb to define additional bump
 * parameters
 *
 * @details The parabolic equation used to calculate the bump in cs^2 is a
 * piecewise function, with its conditions dependent on the value in cs^2 at the
 * start of the bump, as compared to the plateau value. As such, the cs^2 at
 * bump_start for the CRUST EOS is needed to calculate bump values; this is what
 * the variable f1_n1 stores and refers to.
 */
void Bumpy_EOS_Constructor::get_additional_EOS_params() {
  // Get additional helpful value to define the bump
  eos_params.f1_n1 = crust_cs2_of_nb.yofx(eos_params.bump_start);
}

/**
 * @brief Constructs and stores a GSL interpolation object for a bumpy EoS from
 * a given set of parameters.
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

  e_of_nb->initialize(static_cast<const gsl_interp_type*>(gsl_interp_steffen),
                      eos.nb, eos.epsilon);
  eos.eps_c1 = e_of_nb->yofx(eos_params.nbc1);
  eos.eps_c2 = e_of_nb->yofx(eos_params.nbc2);
  delete e_of_nb;
}

/* --------------------------- Protected Functions --------------------------
 */

/**
 * @brief Interpolates the baryon number density and calculates the expected
 * cs^2 values for the given bump parameters.
 *
 * @details Uses the crust table points prior to the start of the bump,
 * interpolates 10 points up to the start of the bump (or maximum central baryon
 * number density, whichever is smaller), interpolates 1000 points on the bump,
 * then interpolates 100 points on the plateau (if the maximum central baryon
 * number denisty hasn't been reached yet.) Calls the quadratic bump calculator
 * in the region of the bump, and simply uses either crust or plateau values
 * outside of this range.
 */
void Bumpy_EOS_Constructor::build_cs2_one_quad_bump() {
  // Use the crust table prior to the start of the bump
  double max_table_value = std::min(eos_params.max_nb, eos_params.bump_start);
  auto nearest_index = std::distance(
      crust.nb.begin(),
      std::lower_bound(crust.nb.begin(), crust.nb.end(), max_table_value));
  // Save the index at which crust values stop being used
  injection_index = nearest_index;
  eos.nb.insert(eos.nb.end(), crust.nb.begin(),
                crust.nb.begin() + nearest_index);
  eos.cs2.insert(eos.cs2.end(), crust.cs2.begin(),
                 crust.cs2.begin() + nearest_index);

  // Interpolate up to the bump
  if (eos.nb.back() < eos_params.bump_start) {
    // Get the last saved baryon number density
    double nb = eos.nb.back();
    // Define interpolation parameters
    double max_table_interp_value =
        std::min(eos_params.bump_start, eos_params.max_nb);
    int steps_before_bump = 10;
    // The +1 is just to make sure the ending value is within range.
    double nb_step_before_bump =
        1 + (max_table_interp_value - nb) / steps_before_bump;

    // Reserve/resize vectors
    eos.nb.reserve(eos.nb.size() + steps_before_bump);
    eos.cs2.reserve(eos.cs2.size() + steps_before_bump);

    // Interpolate baryon number density and cs^2
    std::generate_n(std::back_inserter(eos.nb), steps_before_bump,
                    [&]() { return nb += nb_step_before_bump; });
    std::transform(eos.nb.end() - steps_before_bump, eos.nb.end(),
                   std::back_inserter(eos.cs2),
                   [](double nb) { return crust_cs2_of_nb.yofx(nb); });
  }
  // Calculate values during the bump
  if ((eos.nb.back() < eos_params.bump_end) ||
      (eos.nb.back() < eos_params.max_nb)) {
    // Get the last saved baryon number density
    double nb = eos.nb.back();
    // Define interpolation parameters
    double max_bump_value = std::min(eos_params.bump_end, eos_params.max_nb);
    int steps_bump = 1000;
    // The +1 is just to make sure the ending value is within range.
    double nb_bump_step = 1 + (max_bump_value - nb) / steps_bump;

    // Reserve/resize vectors
    eos.nb.reserve(eos.nb.size() + steps_bump);
    eos.cs2.reserve(eos.cs2.size() + steps_bump);

    // Interpolate baryon number density
    std::generate_n(std::back_inserter(eos.nb), steps_bump,
                    [&]() { return nb += nb_bump_step; });
    // Calculate cs2 values
    std::transform(
        eos.nb.end() - steps_bump, eos.nb.end(), std::back_inserter(eos.cs2),
        [get_quadratic_bump_point =
             &Bumpy_EOS_Constructor::get_quadratic_bump_point,
         this](double nb) { return this->get_quadratic_bump_point(nb); });
  }  // Insert plateau values once bump ends
  if (eos.nb.back() < eos_params.max_nb) {
    // Get the last saved baryon number density
    double nb = eos.nb.back();
    // Define interpolation parameters
    int steps_plat = 100;
    // The +1 is just to make sure the ending value is within range.
    double nb_plat_step = 1 + (eos_params.max_nb - nb) / steps_plat;

    // Reserve/resize vectors
    eos.nb.reserve(eos.nb.size() + steps_plat);
    eos.cs2.reserve(eos.cs2.size() + steps_plat);

    // Interpolate baryon number density
    std::generate_n(std::back_inserter(eos.nb), nb_plat_step,
                    [&]() { return nb += nb_plat_step; });
    // Fill with plateau value
    std::fill_n(std::back_inserter(eos.cs2), steps_plat, eos_params.plat);
  }
}

/**
 * @brief Calculates cs^2 value from point in baryon number density for given
 * EoS parameters.
 *
 * @details Defines a quadratic parabola to calculate cs2 points. The quadratic
 * is piecewise due to 2 roots resulting from whether the bump magnitude is
 * greater or less than the cs2 value at the start of the bump. For the
 * purposes of defining the function, it is useful to consider that the
 * bump_start parameter is "n1", the bump_end is "n2", and thus the cs2 value
 * at n1 is "f1_n1" and the plateau is then equivalent to the cs2 value at n2,
 * "f3_n2." Then, the quadratic is found by simply solving f1_n1 = a1 + a2 n_1
 * + a3 n_1^2 and f3_n2 = a1 + a2 n_2 + a3 n_2^2 (with the condition that the
 * bump magnitude is the peak of the parabola).
 *
 * @param nb Baryon number density at which to evaluate
 * @return double
 */
double Bumpy_EOS_Constructor::get_quadratic_bump_point(const double& nb) {
  // Define helpful constants
  const double dist_start = (nb - eos_params.bump_start);
  const double dist_start_2 = dist_start * dist_start;

  const double dist_end = (nb - eos_params.bump_end);
  const double dist_end_2 = dist_end * dist_end;

  const double width_2 = eos_params.bump_width * eos_params.bump_width;

  // Define constant terms in the function
  const double term1 = eos_params.plat * dist_start_2;
  const double term2 = eos_params.f1_n1 * dist_end_2;
  const double term3 = 2 * eos_params.bump_magnitude * dist_start * dist_end;
  const double constant_terms = term1 + term2 - term3;

  // Define the function roots (will change sign depending on conditions)
  const double func_root = 2 * dist_start * dist_end *
                           sqrt((eos_params.bump_magnitude - eos_params.f1_n1) *
                                (eos_params.bump_magnitude - eos_params.plat));

  // Calculate value of piecewise function
  double cs2;
  if (eos_params.f1_n1 < eos_params.bump_magnitude) {
    cs2 = (1 / width_2) * (constant_terms - func_root);
  } else {
    cs2 = (1 / width_2) * (constant_terms + func_root);
  }

  return cs2;
}

/* --------------------------------------------------------------------------
 */
/*                       ObservablesIntegrator Functions */
/* --------------------------------------------------------------------------
 */

/* ---------------------------- Public Functions ----------------------------
 */

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
 * @brief Integrates for the observable mass and radius based on the stored
 * EoS.
 *
 * @details Sets some small starting radius (R_epsilon) and calculates the
 * corrected starting pseudo-enthalpy, h_c - h_epsilon (h_epsilon from
 * R_epsilon). Then, uses a classic fixed-step Runge-Kutta (rk4) routine to
 * integrate the TOV equations in terms of pseudo-enthalpy from the corrected
 * starting enthalpy to 0 (the surface of the star). Default number of steps
 * for integration is 10,000.
 *
 * @param first_observables Eigen vector to store the first set of
 * calculated neutron star observables (Radius, Mass, Tidal).
 * @param second_observables Eigen vector to store the second set of
 * calculated neutron star observables (Radius, Mass, Tidal).
 * @param CurveIsNegative Boolean to store if the slope of the m(epsilon) curve
 * is negative.
 */
void ObservablesIntegrator::integrate_for_observables(
    Eigen::Vector3d& first_observables, Eigen::Vector3d& second_observables,
    bool& CurveIsNegative) {
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

  /* // Define enthalpy shift values
  double h1_shift = hc_1_shift - x1 * (e_of_h.yofx(hc_1_shift) +
                                       3.0 * p_of_h.yofx(hc_1_shift));
  double h2_shift = hc_2_shift - x1 * (e_of_h.yofx(hc_2_shift) +
                                       3.0 * p_of_h.yofx(hc_2_shift));

  // Define mass shift values
  const double M_shift_1 = x2 * e_of_h.yofx(h1_shift);
  const double M_shift_2 = x2 * e_of_h.yofx(h2_shift); */

  // Define step sizes
  double N_steps = 5e4;
  double h1_step = -h1 / N_steps;
  double h2_step = -h2 / N_steps;

  // TODO: (Kaitlyn) Look into replacing rk4 routine with GSL or other
  // integration functionality.

  /* -------------------- Integrate for First Observables -------------------
   */

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

  CurveIsNegative = false;
  /* Eigen::Vector3d first_shifts;
  Eigen::Vector3d second_shifts;

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
  } */
}

/* -------------------------- Protected Functions ------------------------- */

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
  e_of_h.initialize(static_cast<const gsl_interp_type*>(gsl_interp_steffen),
                    enthalpy, epsilon);
  p_of_h.initialize(static_cast<const gsl_interp_type*>(gsl_interp_steffen),
                    enthalpy, pressure);

  // Get central enthalpy values
  Interpolation h_of_e(static_cast<const gsl_interp_type*>(gsl_interp_steffen),
                       epsilon, enthalpy);
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
                                                  const Eigen::Vector3d& y,
                                                  Eigen::Vector3d& f) {
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
 * @param t Independent variable which functions are being taken with respect
 * to
 * @param y Functions given in terms of the independent variable t
 * @param h Step size for the rk4 routine
 */
void ObservablesIntegrator::rk4(double& t, Eigen::Vector3d& y, double h) {
  // Calculate step size fractions
  const double h_6 = h / 6.0;
  const double h_2 = h / 2.0;

  // Define arrays to store rk4 term results in
  Eigen::Vector3d k1;
  Eigen::Vector3d k2;
  Eigen::Vector3d k3;
  Eigen::Vector3d k4;

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
 * @param adimensionalize Boolean determining whether to convert to solar
 * masses (True) or convert from solar masses (False).
 * @return The converted double.
 */
double ObservablesIntegrator::convert_dimensions(double value, std::string unit,
                                                 bool adimensionalize) {
  // Calculate and store the conversion factor
  double conversion_factor;

  if (unit == "km") {
    conversion_factor = 1e3 * (c * c / GWAT_G_SI) * (1 / GWAT_MSUN_SI);
  } else if (unit == "MeV/fm^3") {
    conversion_factor = (GWAT_eV * 1e51) *
                        sqrt(pow(GWAT_G_SI, 6) / pow(c, 16)) *
                        (GWAT_MSUN_SI * GWAT_MSUN_SI);
  } else {
    std::cerr << "Unsupported unit type passed to "
                 "ObservablesIntegrator::adimensionalize. Returning "
                 "unconverted value. \n";
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

/* ------------------------------------------------------------------------ */
/*                       Interpolation Class Functions                      */
/* ------------------------------------------------------------------------ */

/* --------------------------- Public Functions --------------------------- */

// Constructors

Interpolation::Interpolation(const gsl_interp_type* interp_type,
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

  if (std::isnan(x)) {
    throw std::invalid_argument(
        "A NaN value was passed to an Interpolation object. Aborting to "
        "prevent a segmentation fault.");
  }

  result = gsl_spline_eval(spline, x, acc);

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

  if (std::isnan(x)) {
    throw std::invalid_argument(
        "A NaN value was passed to an Interpolation object. Aborting to "
        "prevent a segmentation fault.");
  }

  result = gsl_spline_eval_deriv(spline, x, acc);

  return result;
}

/**
 * @brief Function to initialize the GSL splines for interpolation
 *
 * @param interp_type
 * @param x
 * @param y
 */
void Interpolation::initialize(const gsl_interp_type* interp_type,
                               std::vector<double> x, std::vector<double> y) {
  type = interp_type;
  // Size of the independent variable vector
  size = x.size();

  // Allocating GSL interpolation accelerator
  acc = gsl_interp_accel_alloc();
  if (!acc) {
    std::cerr << "Failed to allocate memory for GSL accelerator.\n";
    throw std::bad_alloc();
  }

  // Allocating GSL spline of specified type and size
  spline = gsl_spline_alloc(type, size);

  // Initializing GSL spline with given data points (x, y) and size
  gsl_spline_init(spline, x.data(), y.data(), size);

  int status = gsl_spline_init(spline, x.data(), y.data(), size);
  if (status != GSL_SUCCESS) {
    std::cerr << "GSL error code: " << status << "invoked.";
    if (status == GSL_EINVAL) {
      throw std::invalid_argument(
          "Invalid data passed to create GSL spline object.");
    } else if (status == GSL_ENOMEM) {
      std::cerr << "Failed to allocate memory for GSL accelerator.\n";
      throw std::bad_alloc();
    } else {
      throw std::logic_error("GSL spline encountered a domain or range error.");
    }
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

/* ------------------------------------------------------------------------ */
/*                                End of File                               */
/* ------------------------------------------------------------------------ */