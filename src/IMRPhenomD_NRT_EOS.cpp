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
#include "IMRPhenomD_NRT.h"
#include "util.h"
#include <adolc/adouble.h>
#include <adolc/drivers/drivers.h>
#include <adolc/taping.h>
#include <cmath>
#include <complex>
#include <gsl/gsl_randist.h>
#include <iostream>
#include <math.h>
#include <stdexcept>
#include <valarray>
/* -------------------------------------------------------------------------- */
/*                         IMRPhenomD_NRT_EOS METHODS                         */
/* -------------------------------------------------------------------------- */

/* ----------------------------- Public Methods ----------------------------- */

/**
 * @brief Overrides the IMRPhenomD_NRT waveform creation to add EoS
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
int IMRPhenomD_NRT_EOS<T>::construct_waveform(T *frequencies, int length,
                                              std::complex<T> *waveform,
                                              source_parameters<T> *params) {
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
void IMRPhenomD_NRT_EOS<T>::get_m_love(source_parameters<T> &params) {
  // Define objects needed to store information
  string eos_filepath = "/opt/gw_analysis_tools/data/eos.csv";

  // Interpolated pressure in terms of epsilon
  Interpolation *p_of_e = new Interpolation;

  // Central epsilon values
  double ec1;
  double ec2;

  // Constructs the bumpy EoS and stores that information in interpolated object
  Bumpy_EOS_Constructor *eos = new Bumpy_EOS_Constructor(eos_filepath);

  eos->construct_EOS(params->bump_mag, params->bump_width, params->bump_offset,
                     params->plat, params->nbc1, params->nbc2);
  eos->get_EOS(*p_of_e, ec1, ec2);

  delete eos;

  // Get observable values
  ObservablesIntegrator *integrator =
      new ObservablesIntegrator(ec1, ec2, *p_of_e);
  integrator->integrate_for_observables(params->mass1, params.mass2,
                                        params.tidal1, params.tidal2, 0.001);
  delete integrator;

  // Update dependent observable parameters

  // TODO: I'm doing this instead of calling populate_source_params because I
  // don't want to recalculate *everything*, but maybe it'd be more clear to do
  // it anyway? At the cost of efficiency obviously...
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

/* -------------------------------------------------------------------------- */
/*                           EOS_CONSTRUCTOR METHODS                          */
/* -------------------------------------------------------------------------- */

/* ----------------------------- Public Methods ----------------------------- */

/**
 * @brief Incomplete
 *
 * @details Incomplete
 *
 * @param[in]
 * @param[out]
 */
EOS_Constructor::EOS_Constructor(string EOS_filepath) {
  // Reads the EOS file in ROW-MAJOR ORDER
  vector<vector<double>> fileRead;
  read_file<double>(EOS_filepath, fileRead, ',');

  // Convert to COLUMN-MAJOR ORDER
  vector<vector<double>> EOSvectors;
  transpose_data_to_column_major(fileRead, EOSvectors);

  // * IMPORTANT: Assumes MUSES EoS table convention.
  eos.epsilon = EOSvectors[7];
  eos.pressure = EOSvectors[8];
  eos.nb = EOSvectors[4];
  eos.cs2 = eos_to_cs2_convert(EoS_data.pressure, EoS_data.epsilon);
}

/**
 * @brief Incomplete
 *
 * @details Incomplete
 *
 * @param[in]
 * @param[out]
 */
void EOS_Constructor::get_EOS(Interpolation &p_of_e, double &central_epsilon_1,
                              double &central_epsilon_2) {
  convert_cs2_to_eos();
  p_of_e.initialize((gsl_interp_type *)gsl_interp_steffen, eos.epsilon,
                    eos.pressure);
  central_epsilon_1 = eos.eps_c_1;
  central_epsilon_2 = eos.eps_c_2;
}

/* ---------------------------- Protected Methods --------------------------- */

/**
 * @brief Convert a row-major 2D vector of doubles to column-major order.
 *
 * @details This is intended primarily for converting read-in CSV file data to
 * column-major order, since default for file-read in is row-major order.
 *
 * @param[in] row_major 2D vector of data in row-major order.
 * @param[out] column_major 2D vector of data in column-major order.
 */
void transpose_data_to_column_major(const vector<vector<double>> &row_major,
                                    vector<vector<double>> &column_major) {
  // If the row-major vector passed in is empty, exits the function
  if (row_major.empty()) {
    return;
  }

  // Grabs the number of rows and columns
  size_t num_rows = row_major.size();
  size_t num_cols = row_major[0].size();

  // Resize the column_major vector to table size in column-major order
  column_major.resize(num_cols, vector<double>(num_rows));

  // Perform the transposition
  for (size_t i = 0; i < num_rows; ++i) {
    for (size_t j = 0; j < num_cols; ++j) {
      column_major[j][i] = row_major[i][j];
    }
  }
}

/**
 * @brief Calculate the speed of sound squared from pressure and energy density.
 *
 * @details Pressure and energy density are assumed to be in the same units.
 * cs^2 is given in light speed units and calculated using gsl_interp_steffen.
 *
 * @return Vector of cs^2 values (in c).
 */
void EOS_Constructor::convert_eos_to_cs2() {
  // Initialize interpolation of epsilon and pressure
  Interpolation p_of_e((gsl_interp_type *)gsl_interp_steffen, eos.epsilon,
                       eos.pressure);

  // Calculate cs2 for each point
  for (int i = 0; i < eos.epsilon.size(); i++) {
    double cs2_val = p_of_e.dyofx(eos.epsilon[i]);
    eos.cs2.push_back(cs2_val);
  }
}

/**
 * @brief Convert speed of sound squared values to pressure and energy density
 * values.
 *
 * @details Takes in the pressure and energy density to start calculating for
 * and a list of the baryon number density and cs^2 values to evaluate for. The
 * output vectors are assumed to be empty prior to pass-in. Algorithm assumes
 * that p_base, epsilon_base, and nb_list are all given in units of MeV. All
 * inputs must start at the same baryon number density value.
 */

void EOS_Constructor::convert_cs2_to_eos() {
  // Get starting points for integration
  double nb = eos.nb[0];
  double pressure = 0;
  double epsilon = 0;

  // Clear vectors to fill with new EoS data
  eos.pressure.clear();
  eos.epsilon.clear();

  eos.pressure.push_back(pressure);
  eos.epsilon.push_back(epsilon);

  // Loop through all baryon number density values and perform integration
  for (int i = 1; i < eos.nb.size(); ++i) {
    double delta_nb = eos.nb[i] - nb; // Calculate step size
    double delta_e =
        delta_nb * (epsilon + pressure) / nb; // Calculate energy delta at step

    // Take steps
    nb = eos.nb[i];
    epsilon += delta_e;
    pressure += eos.cs2[i - 1] * delta_e;

    // Update pressure and energy density
    eos.pressure[i] = pressure;
    eos.epsilon[i] = epsilon;
  }
}

/**
 * @brief Convert value from units of fm^-3 to MeV.
 *
 * @param[in] x Value in fm^-3.
 * @return Value converted to MeV.
 */
void EOS_Constructor::convert_fm3_to_MeV(double &x) { x = x * pow(197.3, 3); }

/**
 * @brief Convert value from units of nsat to MeV.
 *
 * @param[in] x Value in nsat.
 * @return Value converted to MeV.
 */
void EOS_Constructor::convert_nsat_to_MeV(double &x) {
  double nsat = 0.16;
  x *= nsat;
  convert_fm3_to_MeV(x);
}

/* -------------------------------------------------------------------------- */
/*                        BUMPY_EOS_CONSTRUCTOR METHODS                       */
/* -------------------------------------------------------------------------- */

/* ----------------------------- Public Methods ----------------------------- */

/**
 * @brief Inject a bump into cs2 and get the corresponding pressures and energy
 * densities.
 *
 * @details Constructs cs^2 by stitching together the crust EOS, a single
 * parabolic bump, and a plateau. Modifies pressure and energy density vectors
 * according to cs^2 after bump injection.
 *
 * @param[out] pressure1 Vector of pressures for star 1 (in MeV).
 * @param[out] pressure2 Vector of pressures for star 2 (in MeV).
 * @param[out] epsilon1 Vector of energy densities for star 1 (in MeV).
 * @param[out] epsilon2 Vector of energy densities for star 2 (in MeV).
 * @param[in] params gen_params object containing EOS and bump parameters,
 * including central baryon number densities.
 */
void construct_EOS(double bump_mag, double bump_width, double bump_offset,
                   double plat, double nbc1, double nbc2) {
  // Find start of the bump
  double nb_cutoff = (bump_offset - bump_width / 2);

  // Checks for invalid input (should only be relevant if prior is bad)
  if (std::max(nbc1, nbc2) > eos.nb.back()) {
    throw std::invalid_argument(
        "Encountered central baryon number density greater than the max value "
        "given in the EOS table. Aborting.");
  } else if ((std::max(nbc1, nbc2) < nb_cutoff)) {
    throw std::invalid_argument("Encountered central baryon number density "
                                "less than the start of the bump. Aborting.");
  }

  // Copy over values with correct units
  eos_params.bump_magnitude = convert_nsat_to_MeV(bump_mag);
  eos_params.bump_width = convert_nsat_to_MeV(bump_width);
  eos_params.bump_offset = convert_nsat_to_MeV(bump_offset);
  eos_params.nbc = convert_nsat_to_MeV(std::max(nbc1, nbc2));
  eos_params.plat = plat;

  Interpolation *e_of_nb = new Interpolation;
  Interpolation *p_of_nb = new Interpolation;

  e_of_nb->initialize((gsl_interp_type *)gsl_interp_steffen, eos.nb,
                      eos.epsilon);
  p_of_nb->initialize((gsl_interp_type *)gsl_interp_steffen, eos.nb,
                      eos.pressure);

  // Find location in data where bump starts, clear data after
  auto iterator = std::upper_bound(eos.nb.begin(), eos.nb.end(), nb_cutoff);
  auto start_index = std::distance(eos.nb.begin(), iterator);

  eos.pressure.erase(eos.pressure.begin() + start_index, eos.pressure.end());
  eos.epsilon.erase(eos.epsilon.begin() + start_index, eos.epsilon.end());
  eos.nb.erase(eos.nb.begion() + start_index, eos.nb.end());

  double interpolation_steps = 0.001;

  for (double i = eos_params.nbc; i <= split_limit; i += steps) {
    // For star 1
    if (i <= nb_end1) // If statement prevents exiting limits
    {
      nb_split1.push_back(conversion_fm3_to_MeV(i));
      epsilon_split1.push_back(e_of_nb.yofx(i));
      pressure_split1.push_back(p_of_nb.yofx(i));
    }

    // For star 2
    if (i <= nb_end2) // If statement prevents exiting limits
    {
      nb_split2.push_back(conversion_fm3_to_MeV(i));
      epsilon_split2.push_back(e_of_nb.yofx(i));
      pressure_split2.push_back(p_of_nb.yofx(i));
    }
  }

  // Get bump parameters from gen_params structure
  // IMPORTANT: Bump parameters based on baryon number density values are
  // assumed to be given in units of fm^-3
  double offset = conversion_fm3_to_MeV(eos_params.bump_offset * nsat);
  double magnitude = eos_params.bump_mag;
  double width = conversion_fm3_to_MeV(eos_params.bump_width * nsat);
  double plateau = eos_params.plat;

  // Get new cs2 curve with bump
  build_cs2_one_quad_bump(nb_split1, cs2_star1, width, magnitude, offset,
                          plateau);
  build_cs2_one_quad_bump(nb_split2, cs2_star2, width, magnitude, offset,
                          plateau);

  // Define vectors to store the pressure and energy density values with the
  // bump
  vector<double> p_bump1;
  vector<double> p_bump2;
  vector<double> e_bump1;
  vector<double> e_bump2;

  // Get new bumpy EoS
  cs2_to_eos_convert(pressure_split1.front(), epsilon_split1.front(), nb_split1,
                     cs2_star1, p_bump1, e_bump1);
  cs2_to_eos_convert(pressure_split2.front(), epsilon_split2.front(), nb_split2,
                     cs2_star2, p_bump2, e_bump2);

  // Stitch new bumpy EoS onto the original crust EOS

  // Finding location in data table where prior cut-off occurred
  auto iterator = std::upper_bound(EOSvectors[4].begin(), EOSvectors[4].end(),
                                   nb_split_val); // Get iterator object
  auto split_index =
      std::distance(EOSvectors[4].begin(),
                    iterator); // Convert to index to extract values for p and e

  // Grab elements from the EOS table and add them to the pressure and epsilon
  // vectors passed in

  epsilon1.push_back(0);
  epsilon2.push_back(0);
  pressure1.push_back(0);
  pressure2.push_back(0);

  // For star 1
  pressure1.insert(pressure1.end(), EOSvectors[8].begin(),
                   EOSvectors[8].begin() + split_index);
  epsilon1.insert(epsilon1.end(), EOSvectors[7].begin(),
                  EOSvectors[7].begin() + split_index);
  // For star 2
  pressure2.insert(pressure2.end(), EOSvectors[8].begin(),
                   EOSvectors[8].begin() + split_index);
  epsilon2.insert(epsilon2.end(), EOSvectors[7].begin(),
                  EOSvectors[7].begin() + split_index);

  // Copying bump elements to the final pressure and epsilon vectors

  // For star 1
  pressure1.insert(pressure1.end(), p_bump1.begin(), p_bump1.end());
  epsilon1.insert(epsilon1.end(), e_bump1.begin(), e_bump1.end());
  // For star 2
  pressure2.insert(pressure2.end(), p_bump2.begin(), p_bump2.end());
  epsilon2.insert(epsilon2.end(), e_bump2.begin(), e_bump2.end());
}

/* ---------------------------- Protected Methods --------------------------- */
/**
 * @brief Create a speed of sound squared curve with a parabolic bump added to
 * it.
 *
 * @Details Algorithm assumes list of cs^2 values and list of nb values are the
 * same size.
 *
 * @param[in] nb_list Vector storing baryon number density values (in MeV).
 * @param[in, out] cs2_list Vector containing speed of sound squared values to
 * write over (in c).
 * @param[in] bump_width Width of the injected bump peak (in MeV).
 * @param[in] bump_magnitude Magnitude/"height" of the injected bump (in c).
 * @param[in] bump_offset Value in baryon number density where the bump peak
 * occurs (in MeV).
 * @param[in] bump_plat Plateau to set speed of sound squared to after bump
 * injection (in c).
 */
void construct_EOS::build_cs2_one_quad_bump() {
  // Initialize interpolator function to get cs2 as a function of nb
  Interpolation cs2_of_nb((gsl_interp_type *)gsl_interp_steffen, EoS_data.nb,
                          EoS_data.cs2);

  double f1_n1 =
      cs2_of_nb.yofx(n1); // Value of the crust at the transition point n1

  // Loop through baryon number density values to build cs^2
  for (int i = 0; i < EoS_data.nb.size(); ++i) {
    double nb_temp = EoS_data.nb[i];

    // Use crust EoS before transition
    if (nb_temp < n1) {
      EoS_data.cs2[i] = cs2_of_nb.yofx(nb_temp);
    }
    // Inject quadratic bump after transition
    else if (nb_temp >= n1 && nb_temp <= n2) {
      EoS_data.cs2[i] = f_quad(nb_temp, bump_width, bump_magnitude, bump_offset,
                               bump_plat, f1_n1);
    }
    // Inject plateau after parabola
    else {
      EoS_data.cs2[i] = bump_plat;
    }
  }
}

/**
 * @brief Calculate values for parabolic bump in the speed of sound squared.
 *
 * @details This function calculates the value of cs^2 at a given baryon number
 * density nb.
 *
 * @param[in] nb Baryon number density value to evaluate bump at (in MeV).
 * @param[in] bump_width Width of the injected bump peak (in MeV).
 * @param[in] bump_magnitude Magnitude/"height" of the injected bump (in c).
 * @param[in] bump_offset Value in baryon number density where the bump peak
 * occurs (in MeV).
 * @param[in] bump_plat Plateau to set speed of sound squared to after bump
 * injection (in c).
 * @param[in] f1_n1 Baryon number density value at start of bump (in MeV).
 * @return Value of cs^2 at nb (in c).
 */
void construct_EOS::f_quad(double nb, double f1_n1) {
  // Calculate the expected value of cs2
  return -0.25 *
             ((8 * pow(eos_params.bump_offset, 2) *
               (-1 + 6 * eos_params.bump_magnitude - 3 * f1_n1)) /
                  (3. * eos_params.bump_width) -
              (2 * eos_params.bump_width *
               (-1 + 6 * eos_params.bump_magnitude - 3 * f1_n1)) /
                  3. -
              4 * eos_params.bump_offset * f1_n1 -
              2 * eos_params.bump_width * f1_n1 +
              4 * eos_params.bump_offset * eos_params.plat -
              2 * eos_params.bump_width * eos_params.plat) /
             eos_params.bump_width -
         (((-4 * eos_params.bump_offset *
            (-1 + 6 * eos_params.bump_magnitude - 3 * f1_n1)) /
               (3. * eos_params.bump_width) +
           f1_n1 - eos_params.plat) *
          nb) /
             eos_params.bump_width -
         (2 * (-1 + 6 * eos_params.bump_magnitude - 3 * f1_n1) * pow(nb, 2)) /
             (3. * pow(eos_params.bump_width, 2));
}

// IMPORTANT: GWAT has issues with the difference between C++'s *normal* double
// value and ADOL-C's *adouble* value. This is to prevent conflict and make sure
// everything compiles normally. All new functions added to the
// IMRPhenomD_NRT_EOS class should be templated with "template <class T>" at the
// top (see above functions). This is not necessary for the QLIMR related
// classes.

template class IMRPhenomD_NRT_EOS<double>;
template class IMRPhenomD_NRT_EOS<adouble>;

/* -------------------------------------------------------------------------- */
/*                          TOV and λ̄ FUNCTIONALITY                          */
/* -------------------------------------------------------------------------- */

/* ----------------------------- Public Methods ----------------------------- */

ObservablesIntegrator::ObservablesIntegrator(double epsilon_c1,
                                             double epsilon_c2,
                                             Interpolation eos) {
  p_of_e = eos;
  eps_c1 = epsilon_c1;
  eps_c2 = epsilon_c2;
}

void ObservablesIntegrator::integrate_for_observables(double &m1, double &m2,
                                                      double &tidal1,
                                                      double &tidal2,
                                                      double step_size) {
  valarray<double> observables1;
  valarray<double> observables2;

  // Mass
  observables1[0] = 0;
  observables2[0] = 0;
  // Radius
  observables1[1] = 0;
  observables2[1] = 0;
  // Tidal deformability
  observables1[2] = 0;
  observables2[2] = 0;

  // TODO: There is DEFINITELY a way to condense this so they can be under the
  // same loop.

  // Integrate for first observables

  double h = h_max1;

  while (h >= 0) {
    rk4(h, observables1, step_size);
    h -= step_size;
  }

  m1 = observables1[0] * GWAT_MSUN_SI;
  tidal1 = observables1[2];

  // Integrate for second observables

  double h = h_max2;

  while (h >= 0) {
    rk4(h, observables2, step_size)
  }

  m2 = observables2[0] * GWAT_MSUN_SI;
  tidal2 = observables2[2];
}

/* ---------------------------- Protected Methods --------------------------- */

/**
 * @brief Calculates the derivative of enthalpy(epsilon) with respect to epsilon
 *
 * @details See DOI:10.1086/171882 equation (4). Enthalpy is defined by the
 * integral of dp/(epsilon(p) + p).
 *
 * @param epsilon Energy density value at which to calculate dh/de
 * @return double Enthalpy value for given epsilon
 */
double ObservablesIntegrator::calculate_dh_of_de(double epsilon) {
  double pressure = p_of_e.yofx(epsilon);
  double cs2 = p_of_e.dyofx(epsilon);

  return cs2 / (epsilon + pressure);
}

/**
 * @brief Integrates to find epsilon and pressure in terms of enthalpy
 *
 * @details Simple integrator that uses dh/de to create a table of interpolated
 * epsilon, pressure, and enthalpy values to define interpolated p_of_h and
 * e_of_h objects.
 *
 */
void ObservablesIntegrator::integrate_for_eos_of_h() {
  double delta_e = 0.01;

  double p = 0;
  double h = 0;
  double e = 0;

  vector<double> pressure;
  pressure.push_back(p);
  vector<double> enthalpy;
  enthalpy.push_back(h);
  vector<double> epsilon;
  epsilon.push_back(0);

  for (int i = 1; i <= std::max(eps_c1, eps_c2);) {
    double delta_h = delta_e * calculate_dh_of_de(e);
    p += p_of_e.yofx(e) h += delta_h;
    e += delta_e;

    pressure.push_back(p);
    enthalpy.push_back(h);
    epsilon.push_back(i);

    i += step_size;
  }

  e_of_h.initialize((gsl_interp_type *)gsl_interp_steffen, enthalpy, epsilon);
  p_of_h.initialize((gsl_interp_type *)gsl_interp_steffen, enthalpy, pressure);

  Interpolation h_of_e((gsl_interp_type *)gsl_interp_steffen, epsilon,
                       enthalpy);

  h_max1 = h_of_e(eps_c1);
  h_max2 = h_of_e(eps_c2);
}

valarray<double>
ObservablesIntegrator::Observables_ODE(double h, const valarray<double> y) {
  double e = e_of_h.yofx(h);
  double de = e_of_h.dyofx(h);
  double p = p_of_h.yofx(h);

  // Functions to be solved for (y in RK4)
  double R = y[0];
  double M = y[1];
  double Y = y[2];

  // Output of ODE equations (TOV and tidal deformability)

  valarray<double> f;

  // dM/dh
  f[0] = -(4 * M_PI * e * pow(R, 3) * (R - 2 * M)) /
         (M + 4 * M_PI * pow(R, 3) * p);
  // dR/dh
  f[1] = -(R * (R - 2 * M)) / (M + 4 * M_PI * pow(R, 3) * p);
  // dλ̄/dh (specifically, dλ̄/dR and then multiplying by dR/dh)
  f[2] =
      ((4 * pow(M, 2)) / R + (6 * pow(R, 2)) / (pow((R - 2 * M), 2)) +
       (4 * M * R) * ((10 * e + 26 * p) * (M_PI * pow(R, 2)) - 3) +
       (4 * de * M_PI * pow(R, 3)) / (M + 4 * p * M_PI * pow(R, 3)) +
       (4 * M_PI * pow(R, 4)) * (p * (16 * p * M_PI * pow(R, 2) - 9) - 5 * e) +
       ((1 + 4 * (p - e) * M_PI * pow(R, 2)) * Y) / (2 * M - R) -
       pow(Y, 2) / R) *
      f[1];

  return f;
}

void ObservablesIntegrator::rk4(double t, valarray<double> &y,
                                double step_size) {
  k1 = Observables_ODE(t, y);
  k2 = Observables_ODE(t + (step_size / 2), y + k1 * (step_size / 2));
  k3 = Observables_ODE(t + (step_size / 2), y + k2 * (step_size / 2));
  k4 = Observables_ODE(t + step_size, y + k3 * step_size);

  y += (step_size / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
}

/* -------------------------------------------------------------------------- */
/*                         Interpolation Class Methods                        */
/* -------------------------------------------------------------------------- */

/* ----------------------------- Public Methods ----------------------------- */

// Constructors

Interpolation::Interpolation(gsl_interp_type *interp_type, vector<double> x,
                             vector<double> y) {
  initialize(*interp_type, x, y);
}

Interpolation::Interpolation() {}

// Destructor

Interpolation::~Interpolation() { free(); }

// Function value calculations

// Function to calculate interpolated y value for a given x using GSL spline
double Interpolation::yofx(double x) {
  double result;
  if (x >= 0.0) {
    result = gsl_spline_eval(spline, x, acc);
  } else {
    result = 0.0;
  }
  return result;
}

// Function to calculate dy/dx of interpolated data using GSL spline
double Interpolation::dyofx(double x) {
  double result;
  if (x > 0.0) {
    result = gsl_spline_eval_deriv(spline, x, acc);
  } else {
    result = 0.0;
  }
  return result;
}

// GSL initialization and memory release functions

// Function to initialize the GSL splines for interpolation
void Interpolation::initialize(gsl_interp_type *type, vector<double> x,
                               vector<double> y) {
  type = interp_type;
  // Size of the independent variable vector
  size = x.size();

  // Allocating GSL interpolation accelerator
  acc = gsl_interp_accel_alloc();
  if (!acc) {
    std::cerr << "Error: Failed to allocate memory for accelerator."
              << std::endl;
    exit(EXIT_FAILURE);
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

// Function to release memory of spline and accelerator
void Interpolation::free() {
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
}

/* -------------------------------------------------------------------------- */