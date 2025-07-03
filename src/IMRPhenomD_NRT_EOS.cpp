#include "IMRPhenomD_NRT_EOS.h"
#include "IMRPhenomD_NRT.h"
#include <math.h>
#include <adolc/adouble.h>
#include <adolc/taping.h>
#include <adolc/drivers/drivers.h>
#include <iostream>
#include <cmath>
#include <complex>
#include "util.h"
#include "io_util.h"
#include <gsl/gsl_randist.h>
#include <algorithm>
#include <iterator>

/*! \file
 * File for the generation of equations of state with features in the speed of
 * sound squared as in arXiv:2106.03890 and related papers.
 *
 * Supported features: single parabolic bump ... [more to be possibly added later]
 *
 * Supported crust types: SLy EOS
 *
 * Also included in this file ...
 */

// TODO: cs^2 bump injection can be made more efficient. As it stands, a vector of cs^2 values is calculated for NS 1 and NS 2 separately. However, the only difference in calculation is the central baryon number density value.
// There is a lot of overlap, and the calculation can undoubedtly be made more efficient if it proves to be necessary.

//##############################################################################
// ******************* CONVERSION FROM EOS TO GW PARAMETERS *******************
//##############################################################################

/****************************************************************************************************************
 * !\brief Function to convert EOS parameters to masses and tidal deformability.
 *
 * gen_params object is passed in by reference and modified.
 *
 * Converts bump width, offset, magnitude, and central baryon number densities to masses and tidal deformabilities.
 *
 ****************************************************************************************************************/
template <class T>
void IMRPhenomD_NRT_EOS<T>::get_m_love(gen_params *params)
{
	// Inject bump into cs2 and retrieve new p(e) for star 1 and star 2

	// Define new vectors to store values
	std::vector<double> pressure1;
	std::vector<double> pressure2;
	std::vector<double> epsilon1;
	std::vector<double> epsilon2;

	// Inject bump
	inject_cs2_bump(pressure1, pressure2, epsilon1, epsilon2, params);

	// Start QLIMR routine

	// Specifies the starting radius for integration.
	double R_start = 0.0004;
	// Specifices the central epsilon to integrate at
	double single_epsilon1 = epsilon1.back();
	// Specifices the central epsilon to integrate at
	double single_epsilon2 = epsilon2.back();

	// Initialize MRLevaluator that is based on QLIMR architecture
	Second_Order MRLevaluator1;
	Second_Order MRLevaluator2;

	// Load in the starting values (epsilon and pressure from the EOS and a starting radius value for integration)
	MRLevaluator1.input(epsilon1, pressure1, R_start, single_epsilon1);
	MRLevaluator2.input(epsilon2, pressure2, R_start, single_epsilon2);

	// Interpolates the EoS in terms of the pseudo-enthalpy variable (h)
	MRLevaluator1.initialize_eos((gsl_interp_type *)gsl_interp_steffen);
	MRLevaluator2.initialize_eos((gsl_interp_type *)gsl_interp_steffen);

	// Runs the TOV integrator
	MRLevaluator1.TOV_Integrator(MRLevaluator1.params.single_epsilon, &MRLevaluator1.EoS);
	MRLevaluator2.TOV_Integrator(MRLevaluator2.params.single_epsilon, &MRLevaluator2.EoS);

	// Runs the TidalLove integrator
	MRLevaluator1.TidalLove_Integrator(&MRLevaluator1.fun);
	MRLevaluator2.TidalLove_Integrator(&MRLevaluator2.fun);

	// Note that output is *tidal deformability* (dimensionless)
	// Please see the QLIMR docs for more information: https://ce.musesframework.io/docs/modules/qlimr/contents/2_Physics_Overview.html.
	params->mass1 = MRLevaluator1.NS_M;		// This is given in solar masses
	params->mass2 = MRLevaluator2.NS_M;		// Same
	params->tidal1 = MRLevaluator1.NS_Lbar; // This has no dimensionality
	params->tidal2 = MRLevaluator2.NS_Lbar; // Same
}

/****************************************************************************************************************
 * !\brief Function to inject a bump into cs2 and get the corresponding pressures and energy densities.
 *
 * Vectors to store the pressures and energy densities are passed in by reference. Assumed to be empty prior to pass-in.
 *
 * Assumes there is a SLy crust EOS data file stored in "/data" as "eos.csv" in MUSES EOS table convention.
 * See https://ce.musesframework.io/docs/modules/qlimr/contents/5_Parameters.html for MUSES EOS table convention.
 *
 * Constructs cs^2 by stitching together the crust EOS, a single feature (parabolic bump), and a plateau.
 *
 * Modifies pressure and energy density vectors according to cs^2 after bump injection.
 *
 ****************************************************************************************************************/
template <class T>
void IMRPhenomD_NRT_EOS<T>::inject_cs2_bump(std::vector<double> &pressure1, /**< vector of pressures for star 1 */
											std::vector<double> &pressure2, /**< vector of pressures for star 2 */
											std::vector<double> &epsilon1,	/**< vector of energy densities for star 1 */
											std::vector<double> &epsilon2,	/**< vector of energy densities for star 2 */
											gen_params *params 				/**< object to store EOS parameters */)
{
	// Specifies the filepath to read the EOS data file
	string filename = "../data/eos.csv";
	// Initializes a vector to read EOS file data.
	vector<vector<double>> fileRead;

	// Reads the EOS file in **ROW-MAJOR ORDER** (this is the default for how CSV files are read in C++)
	read_file(filename, fileRead, ',');

	// Initialize vector to store EOS data in **COLUMN-MAJOR ORDER**
	vector<vector<double>> EOSvectors;

	// Transpose the row-major data to column-major
	transpose_data_to_column_major(fileRead, EOSvectors);

	// IMPORTANT: The following lines assume MUSES EOS table convention!

	// Grab the epsilon and pressure vectors from the column-major vector
	vector<double> epsilon = EOSvectors[7];
	vector<double> pressure = EOSvectors[8];
	// Grab baryon number density
	vector<double> nb = EOSvectors[4];

	// Initialize interpolator object to get energy density as a function of nb
	Interpolation e_of_nb;
	e_of_nb.initialize((gsl_interp_type *)gsl_interp_steffen, nb, epsilon);

	// Initialize interpolator object to get pressure as a function of nb
	Interpolation p_of_nb;
	p_of_nb.initialize((gsl_interp_type *)gsl_interp_steffen, nb, pressure);

	// IMPORTANT: The following assumes a SLy EOS. This would need to be updated if the EOS used is ever changed.

	// The SLy EOS table does not interpolate well for values below 0.5 nsat.
	// As such, we define a cutoff point at 0.5 nsat in baryon number density and interpolate values above (since the nb values are sparse)

	double nsat = 0.16;				  // Defining nsat in fm^-3
	double nb_split_val = 0.5 * nsat; // Defining cut-off value
	double nb_end1 = params->nbc1;	  // Getting upper limit of nb values for star 1
	double nb_end2 = params->nbc2;	  // Getting upper limit of nb values for star 2
	double steps = 0.005;			  // Defining step-size for interpolation in fm^-3 units

	// Defining vectors to store new split values for star 1 and star 2

	// Baryon number density
	std::vector<double> nb_split1;
	std::vector<double> nb_split2;
	// Pressure
	std::vector<double> pressure_split1;
	std::vector<double> pressure_split2;
	// Energy density
	std::vector<double> epsilon_split1;
	std::vector<double> epsilon_split2;

	// Get greater of the two nbc values, and use this to initialize star 1 and star 2 vectors simultaneously
	double split_limit = std::max(nb_end1, nb_end2);

	// Defining new nb vector and getting the interpolated values for epsilon and pressure
	// IMPORTANT: Baryon number density values are converted to MeV, based on assumptions of cs2_to_eos_convert

	for (double i = nb_split_val; i <= split_limit; i += steps)
	{
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

	// Get original cs2 values without bump
	vector<double> cs2_star1 = eos_to_cs2_convert(pressure_split1, epsilon_split1);
	vector<double> cs2_star2 = eos_to_cs2_convert(pressure_split2, epsilon_split2);

	// Get bump parameters from gen_params structure
	// IMPORTANT: Bump parameters based on baryon number density values are assumed to be given in units of fm^-3
	double offset = conversion_fm3_to_MeV(params->bump_offset);
	double magnitude = params->bump_mag;
	double width = conversion_fm3_to_MeV(params->bump_width);
	double plateau = params->plat;

	// Get new cs2 curve with bump
	build_cs2_one_quad_bump(nb_split1, cs2_star1, width, magnitude, offset, plateau);
	build_cs2_one_quad_bump(nb_split2, cs2_star2, width, magnitude, offset, plateau);

	// Define vectors to store the pressure and energy density values with the bump
	vector<double> p_bump1;
	vector<double> p_bump2;
	vector<double> e_bump1;
	vector<double> e_bump2;

	// Get new bumpy EoS
	cs2_to_eos_convert(pressure_split1.begin(), epsilon_split1.begin(), nb_split1, cs2_star1, p_bump1, e_bump1);
	cs2_to_eos_convert(pressure_split2.begin(), epsilon_split2.begin(), nb_split2, cs2_star2, p_bump2, e_bump2);

	// Stitch new bumpy EoS onto the original crust EOS

	// Finding location in data table where prior cut-off occurred
	auto iterator = std::upper_bound(nb.begin(), nb.end(), nb_split_val); // Get iterator object
	auto split_index = std::distance(nb.begin(), iterator);				  // Convert to index to extract values for p and e

	// Grab elements from the EOS table and add them to the pressure and epsilon vectors passed in

	for (int i; i < split_index; i++)
	{
		// Get pressure and epsilon values from EOS table
		auto p_value = pressure[i];
		auto e_value = epsilon[i];

		// Add values for star 1
		pressure1.push_back(p_value);
		epsilon1.push_back(e_value);

		// Add values for star 2
		pressure2.push_back(p_value);
		epsilon2.push_back(e_value);
	}

	// Copying bump elements to the final pressure and epsilon vectors

	// For star 1
	pressure1.insert(pressure1.end(), p_bump1.begin(), p_bump1.end());
	epsilon1.insert(epsilon1.end(), e_bump1.begin(), e_bump1.end());
	// For star 2
	pressure2.insert(pressure2.end(), p_bump2.begin(), p_bump2.end());
	epsilon2.insert(epsilon2.end(), e_bump2.begin(), e_bump2.end());
}

/****************************************************************************************************************
 * !\brief Function to convert a row-major 2D vector of doubles to column-major order
 *
 * Row-major vector and column-major vector as passed in by reference.
 *
 * This is intended primarly for converting read-in CSV file data to column-major order, 
 * since default for file-read in is row-major order.
 *
 ****************************************************************************************************************/
template <class T>
void IMRPhenomD_NRT_EOS<T>::transpose_data_to_column_major(const std::vector<std::vector<double>> &row_major, /**< 2D vector of data in row-major order */
														   std::vector<std::vector<double>> &column_major 	  /**< 2D vector to store data converted into column-major order */)
{
	// If the row-major vector passed in is empty, exits the function
	if (row_major.empty())
	{
		return;
	}

	// Grabs the number of rows and columns
	size_t num_rows = row_major.size();
	size_t num_cols = row_major[0].size();

	// Resize the column_major vector to table size in column-major order
	column_major.resize(num_cols, std::vector<double>(num_rows));

	// Perform the transposition
	for (size_t i = 0; i < num_rows; ++i)
	{
		for (size_t j = 0; j < num_cols; ++j)
		{
			column_major[j][i] = row_major[i][j];
		}
	}
}

/****************************************************************************************************************
 * !\brief Function to calculate the speed of sound squared from pressure and energy density
 *
 * Takes in double vectors of pressure and energy density, assumed to be in same units.
 *
 * Outputs a vector of doubles storing the cs^2 values in light speed units.
 *
 ****************************************************************************************************************/
template <class T>
std::vector<double> IMRPhenomD_NRT_EOS<T>::eos_to_cs2_convert(std::vector<double> pressure, /**< Vector of pressure values */
															  std::vector<double> epsilon 	/**< Vector of energy density values */)
{
	// Initialize interpolator object to get derivative of pressure as a function of energy density
	Interpolation p_of_e;
	p_of_e.initialize((gsl_interp_type *)gsl_interp_steffen, epsilon, pressure);

	// Initialize vector of doubles to store cs^2
	std::vector<double> cs2;

	// Loops through each energy density point and takes the derivative, adds to cs2 vector
	for (int i = 0; i < epsilon.size(); i++)
	{
		double cs2_val = p_of_e.dyofx(epsilon[i]);
		cs2.push_back(cs2_val);
	}

	// Returns vector of cs^2 values
	return cs2;
}

/*!\brief Function to convert value from units of fm^-3 to MeV
 *
 * Takes in double value, outputs double value.
 */
template <class T>
double IMRPhenomD_NRT_EOS<T>::conversion_fm3_to_MeV(double x)
{
	double x_new = x * pow(197.3, 3);
	return x_new;
}

/****************************************************************************************************************
 * !\brief Function to create a speed of sound squared curve with a parabolic bump added to it
 *
 * Takes in a vector of baryon number density values, a vector of cs^2 values from a crust EOS, and bump parameters.
 *
 * Takes in parameters to define the bump width, magnitude, peak location (offset), and plateau value after the bump.
 *
 * cs^2 vector is passed in by reference and loaded with the bump injection values.
 *
 * Assumes list of cs^2 values and list of nb values are the same size.
 *
 ****************************************************************************************************************/
template <class T>
void IMRPhenomD_NRT_EOS<T>::build_cs2_one_quad_bump(std::vector<double> nb_list,   /**< Vector storing baryon number density values */
													std::vector<double> &cs2_list, /**< Vector containing speed of sound squared values to write over */
													double bump_width,			   /**< Width of the injected bump peak */
													double bump_magnitude,		   /**< Magnitude/"height" of the injected bump */
													double bump_offset,			   /**< Value in baryon number density where the bump peak occurs */
													double bump_plat 			   /**< Plateau to set speed of sound squared to after bump injection */)
{
	// Initialize interpolator function to get cs2 as a function of nb
	Interpolation cs2_of_nb;
	cs2_of_nb.initialize((gsl_interp_type *)gsl_interp_steffen, nb_list, cs2_list);

	// Obtain parameters for the parabola.
	double n1 = bump_offset - bump_width / 2; // Number density at which the parabola begins
	double n2 = bump_offset + bump_width / 2; // Number density at which the parabola ends
	double f1_n1 = cs2_of_nb.yofx(n1);		  // Value of the crust at the transition point n1

	// Loop through baryon number density values to build cs^2

	for (int i = 0; i < nb_list.size(); ++i)
	{
		// Obtain the baryon number density at point i
		double nb_temp = nb_list[i];

		// Check if the baryon number density is smaller than the transition point
		if (nb_temp < n1)
		{
			cs2_list[i] = cs2_of_nb.yofx(nb_temp); // Before transition point, use the crust values
		}
		// Check if the baryon number density is in the parabolic region (between transition points)
		else if (nb_temp >= n1 && nb_temp <= n2)
		{
			cs2_list[i] = f_quad(nb_temp, bump_width, bump_magnitude, bump_offset, bump_plat, f1_n1); // Calculate the value for the parabola
		}
		// After the parabola, cs^2 is in the plateau region
		else
		{
			cs2_list[i] = bump_plat; // After parabola, use plateau value
		}
	}
}

/****************************************************************************************************************
 * !\brief Function to calculate values for parabolic bump in the speed of sound squared
 *
 * Takes in baryon number density value, bump parameters, and the starting baryon number density value for the bump.
 *
 * Parameters to define the bump are width, magnitude, peak location (offset), and plateau value after the bump.
 *
 ****************************************************************************************************************/
template <class T>
double IMRPhenomD_NRT_EOS<T>::f_quad(double nb,				/**< Baryon number density value to evaluate bump at */
									 double bump_width,		/**< Width of the injected bump peak */
									 double bump_magnitude, /**< Magnitude/"height" of the injected bump */
									 double bump_offset,	/**< Value in baryon number density where the bump peak occurs */
									 double bump_plat,		/**< Plateau to set speed of sound squared to after bump injection */
									 double f1_n1 			/**< Baryon number density value at start of bump */)
{
	// Calculate the expected value of cs2
	double cs2_val = -0.25 * ((8 * pow(bump_offset, 2) * (-1 + 6 * bump_magnitude - 3 * f1_n1)) / (3. * bump_width) - 
							  (2 * bump_width * (-1 + 6 * bump_magnitude - 3 * f1_n1)) / 3. - 
							   4 * bump_offset * f1_n1 - 2 * bump_width * f1_n1 + 4 * bump_offset * bump_plat - 2 * bump_width * bump_plat) / bump_width - 
					(((-4 * bump_offset * (-1 + 6 * bump_magnitude - 3 * f1_n1)) / (3. * bump_width) + f1_n1 - bump_plat) * nb) / bump_width -
					(2 * (-1 + 6 * bump_magnitude - 3 * f1_n1) * pow(nb, 2)) / (3. * pow(bump_width, 2));
	return cs2_val;
}

/****************************************************************************************************************
 * !\brief Function to convert speed of sound squared values to pressure and energy density values
 *
 * Takes in the pressure and energy density to start calculating for 
 * and a list of the baryon number density and cs^2 values.
 * 
 * Vectors to store the new calculated pressure and and energy density are passed in by reference. 
 * Assumed to be empty prior to pass-in.
 *
 * Algorithm assumes that p_base, epsilon_base, and nb_list are all given in units of MeV.
 * All inputs must start at the same baryon number density value.
 *
 ****************************************************************************************************************/
template <class T>
void IMRPhenomD_NRT_EOS<T>::cs2_to_eos_convert(double p_base,				 	 /**< Pressure at the start of bump (in MeV) */
											   double epsilon_base,			 	 /**< Energy density at the start of bump (in MeV) */
											   std::vector<double> nb_list,	 	 /**< List of baryon number density values (in MeV) */
											   std::vector<double> cs2_bump, 	 /**< List of speed of sound squared values (in c) */
											   std::vector<double> &p_bump,	 	 /**< Vector to store new pressure values (in MeV) */
											   std::vector<double> &epsilon_bump /**< Vector to store new energy density values (in MeV) */)
{
	// Get starting points for integration
	double p = p_base;			   // Pressure
	double epsilon = epsilon_base; // Energy density
	double nb = nb_list[0];		   // Baryon number density

	// Add starting pressure and energy density values to bump vectors
	p_bump.push_back(p);
	epsilon_bump.push_back(epsilon);

	// Loop through all baryon number density values and perform integration

	for (int i = 1; i < nb_list.size();; ++i)
	{
		double delta_nb = nb_list[i] - nb;				// Calculating the size of the current step
		double delta_e = delta_nb * (epsilon + p) / nb; // Calculating the delta in energy density at the current step

		nb = nb_list[i];				// Taking a step in number density
		epsilon += delta_e;				// Taking a step in energy density
		p += cs2_bump[i - 1] * delta_e; // Taking a step in pressure

		p_bump.push_back(p);			 // Add the pressure calculated at this step to the end of the pressure list
		epsilon_bump.push_back(epsilon); // Add the energy density calculated at this step to the end of the pressure list
	}
}

// IMPORTANT: GWAT has issues with the difference between C++'s *normal* double value and ADOL-C's *adouble* value. This is to prevent conflict and make sure everything compiles normally.
// All new functions added to the IMRPhenomD_NRT_EOS class should be templated with "template <class T>" at the top (see above functions). This is not necessary for the QLIMR related classes.
// Don't ask me how this even became an issue in the first place; I just work here.

template class IMRPhenomD_NRT_EOS<double>;
template class IMRPhenomD_NRT_EOS<adouble>;

//##############################################################################
// *********************** QLIMR TOV and λ̄ FUNCTIONALITY ***********************
//##############################################################################

// ----------------------------------------------------------------------------

Input_QLIMR::Input_QLIMR() {}

// ------------------------- Input_QLIMR: Read yaml params --------------------
void Input_QLIMR::input(vector<double> epsilon, vector<double> pressure, double R_start, double single_epsilon)
{
	params.R_start = adimensionalize(R_start, "km");

	// Check that the pressure and epsilon vectors are the same size
	if (pressure.size() == epsilon.size())
	{
		// Adimensionalize the pressure and epsilon
		for (int i = 0; i < epsilon.size(); i++)
		{
			epsilon[i] = adimensionalize(epsilon[i], "MeV/fm^3");
			pressure[i] = adimensionalize(pressure[i], "MeV/fm^3");
		}

		params.epsilon_col1 = epsilon;
		params.pressure_col2 = pressure;

		params.single_epsilon = adimensionalize(single_epsilon, "MeV/fm^3");
	}
	else
	{
		std::cout << "Error: Epsilon and pressure vectors passed to LMR structure are not the same size. Check read EOS file." << std::endl;
	}
}
// ----------------------------------------------------------------------------

// --------------- Adimensionalize conversion function ------------------------
double Input_QLIMR::adimensionalize(double value, string unit)
{
	double factor;

	if (unit == "MeV/fm^3")
	{
		factor = 1.0 / 346933.783551;
	}
	else if (unit == "km")
	{
		factor = 1.0 / 1.47663;
	}
	else if (unit == "g/cm^3")
	{
		factor = 1.0 / 6.17625e+17;
	}
	else if (unit == "MeV")
	{
		factor = 8.96162e-61;
	}
	else if (unit == "1/fm^3")
	{
		factor = 3.216297e54;
	}
	else if (unit == "-")
	{
		factor = 1.0;
	}
	else
	{
		std::cout << "no unit match to adimensionalize" << std::endl;
		exit(0);
	}

	return factor * value;
}

// --------------- Dimensionalize conversion function ------------------------
double Input_QLIMR::dimensionalize(double value, string unit)
{
	double factor;

	if (unit == "MeV/fm^3")
	{
		factor = 346933.783551;
	}
	else if (unit == "km")
	{
		factor = 1.47663;
	}
	else if (unit == "g/cm^3")
	{
		factor = 6.17625e+17;
	}
	else if (unit == "MeV")
	{
		factor = 1.0 / 8.96162e-61;
	}
	else if (unit == "1/fm^3")
	{
		factor = 1.0 / 3.216297e54;
	}
	else if (unit == "Hz")
	{
		factor = (299792 / 1.47663); // c [km/s] / lsun [km]
	}
	else if (unit == "-")
	{
		factor = 1.0;
	}
	else
	{
		std::cout << "no unit match to dimensionalize" << std::endl;
		exit(0);
	}

	return factor * value;
}

// ---------------------- Interpolation class methods -------------------------
void Interpolation::initialize(gsl_interp_type *interp_type, vector<double> x, vector<double> y)
{
	type = interp_type;
	// Size of the independent variable vector
	size = x.size();

	// Allocating GSL interpolation accelerator
	acc = gsl_interp_accel_alloc();
	if (!acc)
	{
		std::cerr << "Error: Failed to allocate memory for accelerator." << std::endl;
		exit(EXIT_FAILURE);
	}

	// Allocating GSL spline of specified type and size
	spline = gsl_spline_alloc(type, size);

	// Initializing GSL spline with given data points (x, y) and size
	gsl_spline_init(spline, x.data(), y.data(), size);
}

// Function to calculate interpolated y value for a given x using GSL spline
double Interpolation::yofx(double x)
{
	double result;
	if (x >= 0.0)
	{
		result = gsl_spline_eval(spline, x, acc);
	}
	else
	{
		result = 0.0;
	}
	return result;
}

// Function to calculate dy/dx of interpolated data using GSL spline
double Interpolation::dyofx(double x)
{
	double result;
	if (x > 0.0)
	{
		result = gsl_spline_eval_deriv(spline, x, acc);
	}
	else
	{
		result = 0.0;
	}
	return result;
}

// Method to release memory of spline and accelerator
void Interpolation::free()
{
	gsl_spline_free(spline);
	gsl_interp_accel_free(acc);
}
//-----------------------------------------------------------------------------

// ---------------------------- Defining dh/dε --------------------------------
double enthalpy_integrand(double epsilon, void *params)
{

	double p = ((Interpolation *)params)->yofx(epsilon);
	double cs2 = ((Interpolation *)params)->dyofx(epsilon);

	return cs2 / (p + epsilon);
}
// ----------------------------------------------------------------------------

// -------------------------- Default EOS constructor --------------------------
EOS::EOS() {}; // Default EOS constructor
// ----------------------------------------------------------------------------

// ----------------------- Parametric EOS constructor -------------------------
void EOS::initialize_eos(gsl_interp_type *type)
{
	// Interpolate p(ε): initialize GSL interpolation spline
	EoS.p_of_e.initialize(type, params.epsilon_col1, params.pressure_col2);

	// Compute interpolated ε(h) and p(h)
	calculate_eos_of_h(&params.epsilon_col1, type);
}
// ----------------------------------------------------------------------------

//------------------- Integration: Finding ε(h) and p(h) ----------------------
void EOS::calculate_eos_of_h(vector<double> *epsilon, gsl_interp_type *type)
{
	double delta_h;
	double error;
	double h = 0.0;

	// Integration workspace size
	size_t integration_workspace_size = 100 * (epsilon->size());

	// GSL variable type for integrating function
	gsl_function F;

	// Function to be integrated
	F.function = &enthalpy_integrand;

	// Passing object parameter which contains p(ε)
	F.params = &EoS.p_of_e;

	// Add h=0 and ε=0 according to definition of h
	EoS.h_vec.push_back(0.0);
	EoS.e_vec.push_back(0.0);
	EoS.p_vec.push_back(0.0);

	for (size_t i = 1; i < epsilon->size() - 1; i++)
	{

		// Allocate memory for GSL integration workspace
		gsl_integration_workspace *w = gsl_integration_workspace_alloc(integration_workspace_size);

		// Perform adaptive quadrature integration using GSL library
		gsl_integration_qag(&F,
							(*epsilon)[i],	   // Lower integration limit
							(*epsilon)[i + 1], // Upper integration limit
							1e-9,			   // Absolute tolerance
							1e-9,			   // Relative tolerance
							1000,			   // Maximal number of subintervals
							6,				   // Integration method (Gauss-Kronrod)
							w,				   // Work space
							&delta_h,		   // Estimated step size
							&error			   // Estimated integration error
		);

		h = h + delta_h;
		EoS.h_vec.push_back(h);
		EoS.e_vec.push_back((*epsilon)[i + 1]);
		EoS.p_vec.push_back(EoS.p_of_e.yofx((*epsilon)[i + 1]));
		gsl_integration_workspace_free(w);
	}

	EoS.h_of_e.initialize(type, EoS.e_vec, EoS.h_vec); // Interpolate h(ε)
	EoS.h_of_p.initialize(type, EoS.p_vec, EoS.h_vec); // Interpolate h(p)
	EoS.e_of_h.initialize(type, EoS.h_vec, EoS.e_vec); // Interpolate ε(h)
	EoS.p_of_h.initialize(type, EoS.h_vec, EoS.p_vec); // Interpolate p(h)
}
//-----------------------------------------------------------------------------

// TOV equations (h-formulation)
int TOV_equations(double h, const double y[], double f[], void *paraM_sol)
{

	// Recast input parameter pointer of the ODE system as an input object
	EOSinterpolation *eos = (EOSinterpolation *)paraM_sol;

	// EoS variables ε(h) and p(h) to be inserted into the ODE system
	double e = eos->e_of_h.yofx(h);
	double p = eos->p_of_h.yofx(h);

	// Functions R(h) and M(h) to be solved
	double R = y[0];
	double M = y[1];

	double x0 = 4 * M_PI * std::pow(R, 3);

	// TOV equations in h-formulation form. Here, f[0]= dR/dh and f[1]= dM/dh.
	f[0] = -R * (-2 * M + R) / (M + x0 * p);
	f[1] = -e * x0 * (-2 * M + R) / (M + p * x0);

	return GSL_SUCCESS;
}

// Initial conditions for TOVh system
TOV::Initial_conditions_TOV TOV::IC_TOV(double epsilon_c)
{

	// Values at exactly the center
	double pc = EoS.p_of_e.yofx(epsilon_c);
	double hc = EoS.h_of_e.yofx(epsilon_c);
	double Cs2 = EoS.p_of_e.dyofx(epsilon_c);

	// ################## Asymptotic solutions ###################

	double R = params.R_start;

	double M = (4.0 / 3.0) * M_PI * epsilon_c * pow(R, 3);

	// -----------------------------------------------------------

	double p = pc - (2.0 / 3.0) * M_PI * pow(R, 2) * (epsilon_c + 3.0 * pc) * (epsilon_c + pc);

	double e = epsilon_c - (2.0 / 3.0) * M_PI * pow(R, 2) * (epsilon_c + 3.0 * pc) * (epsilon_c + pc) / Cs2;

	// ####################### dM/dh, dR/dh at R=Rε ########################

	double dRdh = -(R * (R - 2.0 * M)) / (M + 4.0 * M_PI * pow(R, 3) * p);

	double dMdh = 4.0 * M_PI * e * pow(R, 2) * dRdh;

	// Filling initial conditions and initial step size into structure
	IC_tov.R_start = params.R_start;
	IC_tov.M_start = M;
	IC_tov.h_start = hc - (2.0 / 3.0) * M_PI * pow(R, 2) * (epsilon_c + 3.0 * p);
	IC_tov.h_istep = 0.1 * std::min(abs(R / dRdh), std::abs(M / dMdh));

	// cout << min(abs(R / dRdh), abs(M / dMdh)) << endl;
	// cout << hc << endl;
	// //exit(1);

	return IC_tov;
}

//  TOV Integrator using pseudo-enthalpy (h) formulation
void TOV::TOV_Integrator(double epsilon_c, EOSinterpolation *eos)
{

	// Number of dependent variables of the ODE system to be solved
	const int dim = 2;

	// Defining GSL variables: system, step, control and evolve
	gsl_odeiv2_system sys = {TOV_equations, NULL, dim, eos};
	gsl_odeiv2_step *s = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rkf45, dim);
	gsl_odeiv2_control *c = gsl_odeiv2_control_y_new(1e-10, 1e-10);
	gsl_odeiv2_evolve *e = gsl_odeiv2_evolve_alloc(dim);

	// Obtain initial conditions
	IC_TOV(epsilon_c);

	// Initial conditions array for r and m at R_start
	double y[2] = {IC_tov.R_start, IC_tov.M_start};

	// Initial value of h to begin the integration
	double h = IC_tov.h_start;

	// Initial step size for the GSL adaptative step size algorithm
	double h_istep = -1.0 * IC_tov.h_istep;

	// Final value of h to stop integration
	double h_stop = 0.0;

	// Solution vectors to store h, R, M, p, e and ν at each step in h
	vector<double> h_sol, R_sol, M_sol, p_sol, e_sol, nu_sol;

	// Add initial conditions at h_start
	h_sol.push_back(h);
	R_sol.push_back(y[0]);
	M_sol.push_back(y[1]);
	p_sol.push_back(eos->p_of_e.yofx((eos->e_of_h.yofx(h))));
	e_sol.push_back(eos->e_of_h.yofx(h));

	// Integration loop from h = h_start up to the surface when h = 0.0
	while (h > h_stop)
	{

		// Solve the ODE system using GSL with rkf45 method at each step
		int status = gsl_odeiv2_evolve_apply(e, c, s,
											 &sys, &h, h_stop, &h_istep, y);

		// Stop solving if something's wrong with the numerical integration
		if (status != GSL_SUCCESS)
			break;

		// Store h values to be used for finding ν(h) */
		h_sol.push_back(h);

		// Store radial distance (R) after each step
		R_sol.push_back(y[0]);

		// Store enclosed mass (M) after each step
		M_sol.push_back(y[1]);

		// Store energy density (ε) after each step
		e_sol.push_back(eos->e_of_h.yofx(h));

		if (e_sol.back() <= 0)
		{
			p_sol.push_back(0);
		}
		else
		{
			// Store pressure (p) after each step
			p_sol.push_back(eos->p_of_e.yofx((eos->e_of_h.yofx(h))));
		}
	}

	// // Check that integration indeed goes up to h = 0 where p = 0.
	// cout << "h = " << h  << " " << "p = " << ps.back() << " " << "e = " <<
	// es.back()  << endl; exit(1);

	// Free GSL variables of the integrator
	gsl_odeiv2_evolve_free(e);
	gsl_odeiv2_control_free(c);
	gsl_odeiv2_step_free(s);

	// Storing total mass M = m(h=0) and radius R = r(h=0)
	NS_R = y[0];
	NS_M = y[1];

	// Obtaining analytical solution for ν(h)
	for (size_t i = 0; i < h_sol.size(); i++)
	{
		nu_sol.push_back(log(1 - 2.0 * (NS_M / NS_R)) - 2 * h_sol[i]);
	}

	// Initialize GSL interpolation spline for m(r) solution
	fun.M_of_R.initialize((gsl_interp_type *)gsl_interp_steffen, R_sol, M_sol);

	// Initialize GSL interpolation spline for p(r) solution
	fun.p_of_R.initialize((gsl_interp_type *)gsl_interp_steffen, R_sol, p_sol);

	// Initialize GSL interpolation spline for ε(r) solution
	fun.e_of_R.initialize((gsl_interp_type *)gsl_interp_steffen, R_sol, e_sol);

	// Initialize GSL interpolation spline for nu(r) solution
	fun.nu_of_R.initialize((gsl_interp_type *)gsl_interp_steffen, R_sol, nu_sol);
};

// ############################################################################
// ################################ λ̄ (l=2) ###################################
// ############################################################################

// ODE for the tidal love function Y
int TidalLove_equation(double R, const double y[], double f[], void *params_ptr)
{

	// Function parameters of the ODE:
	Local_functions *fun = (Local_functions *)params_ptr;

	double e = fun->e_of_R.yofx(R);
	double p = fun->p_of_R.yofx(R);
	double M = fun->M_of_R.yofx(R);
	double de = fun->e_of_R.dyofx(R);

	double x0 = 1.0 / R;
	double x1 = 4 * M_PI;
	double x2 = std::pow(R, 3) * x1;
	double x3 = 2 * M - R;
	double x4 = std::pow(R, 2);
	double x5 = M_PI * x4;
	double x6 = 5 * e;

	double Y = y[0];

	f[0] = -std::pow(Y, 2) * x0 + Y * (4 * x5 * (-e + p) + 1) / x3 + de * x2 / (M + p * x2) + x0 * (4 * std::pow(M, 2) + 4 * M * R * (2 * x5 * (13 * p + x6) - 3) + std::pow(R, 4) * x1 * (p * (16 * p * x5 - 9) - x6) + 6 * x4) / std::pow(x3, 2);

	return GSL_SUCCESS;
}

// Initial Conditions for Y
Second_Order::Initial_conditions_Y Second_Order::IC_Y()
{

	double R = params.R_start;

	// Defining initial step size for adaptative method
	IC_y.R_istep = R / 100;

	return IC_y;
}

// Integrator for Y
void Second_Order::TidalLove_Integrator(Local_functions *fun)
{

	const int dim = 1;

	gsl_odeiv2_system sys = {TidalLove_equation, NULL, dim, fun};
	gsl_odeiv2_step *s = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rkf45, dim);
	gsl_odeiv2_control *c = gsl_odeiv2_control_y_new(1e-12, 1e-12);
	gsl_odeiv2_evolve *e = gsl_odeiv2_evolve_alloc(dim);

	// Computing initial conditions of the ODE
	IC_Y();

	double y[1] = {IC_y.Y_start};
	double R = params.R_start;
	double R_istep = IC_y.R_istep;

	// Defining step and solution vectors
	vector<double> R_sol, Y_sol;

	// Adding initial condition to the solution vector
	R_sol.push_back(R);
	Y_sol.push_back(IC_y.Y_start);

	// Solving the ODE over the specified range
	while (R < NS_R)
	{
		int status = gsl_odeiv2_evolve_apply(e, c, s, &sys, &R, NS_R, &R_istep, y);

		if (status != GSL_SUCCESS)
			break;

		R_sol.push_back(R);
		Y_sol.push_back(y[0]);
		// cout << R << " " << y[0] << endl;
	}

	gsl_odeiv2_evolve_free(e);
	gsl_odeiv2_control_free(c);
	gsl_odeiv2_step_free(s);

	// Evaluating solution Y(R) at the NS surface: NS_Y
	double NS_Y = y[0];

	// Useful definitions for faster calculations
	double C = NS_M / NS_R;
	double C2 = C * C;
	double C3 = C * C2;
	double C4 = C2 * C2;
	double C5 = C * C4;
	double k = (1.0 - 2.0 * C) * (1.0 - 2.0 * C);

	// Dimensionless tidal apsidal constant NS_k2: k2 [-]
	NS_k2 = (8.0 / 5.0) * k * C5 * (2.0 * C * (NS_Y - 1.0) - NS_Y + 2.0) *
          (1.0 / (2.0 * C *
            (4 * (NS_Y + 1.0) * C4 + (6.0 * NS_Y - 4.0) * C3 + (26.0 - 22.0 * NS_Y) * C2 + 3.0 * (5.0 * NS_Y - 8.0) * C - 3.0 * NS_Y + 6.0) -
            3.0 * k * (2 * C * (NS_Y - 1.0) - NS_Y + 2.0) * log(1.0 / (1.0 - 2.0 * C))));

	// Dimensionless tidal love number NS_Lbar: λ̄ [-]
	NS_Lbar = (2.0 / 3.0) * NS_k2 * (1 / C5);

	// Initialize GSL interpolation spline for Y(r) solution
	fun->Y_of_R.initialize((gsl_interp_type *)gsl_interp_steffen, R_sol, Y_sol);
}