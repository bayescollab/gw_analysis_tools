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

/*! \file 
 * File for the generation of equations of state with features in the speed of
 * sound squared as in arXiv:2106.03890 and related papers.
 *
 * Supported features: single Gaussian, ...[add more later]
 *
 * Supported crust types: SLy EOS
 * 
 * Also included in this file ...
 */

template<class T> 
void IMRPhenomD_NRT_EOS<T>::build_cs2_one_bump(source_parameters<T> *sp)
{
  //Will need to add some sort of array/table/pointer thing to pass back and forth the values of cs2
  
  /* A function to construct cs^2 by stitching together an EOS in the crust,
   * a single feature (Gaussian bump), and a plateau.
   *
   * The bump is described by the height/magnitude (bump_mag),
   * the width (bump_width), and the location (bump_offset).
   */

  
  std::cout<<"Calling build_cs2_one_bump"<<std::endl; 
  
}

template<class T> 
void IMRPhenomD_NRT_EOS<T>::cs2_to_eos_convert()
{
  
}

 template<class T> 
void IMRPhenomD_NRT_EOS<T>::get_m_love(gen_params* params)
{
  // A dummy function only for code infrastructure testing!!
  // If you edit what variables this takes, you'll need to update the EOS_testing.cpp script in the injections directory. 
  //TODO replace with something that is physically accurate.
  params->mass1 = params->nbc1; //As a dummy function, I just input nbc1 as mass. 
  params->mass2 = params->nbc2; //Same
  params->tidal1 = params->nbc1*100; //Gives a tidal deformability that depends on mass and has the right order of magnitude. This will be changed later. 
  params->tidal2 = params->nbc2*100; 

}
/*
 template<class T> 
std::array<double, 3> IMRPhenomD_NRT_EOS<T>::get_m_love(vector<double> epsilon, vector<double> pressure)
{
    // Specifies the starting radius for integration.
    double R_start = 0.0004;
    // Specifices the central epsilon to integrate at
    double single_epsilon = epsilon.back();

    // Initialize MRLevaluator that is based on QLIMR architecture
    Second_Order MRLevaluator;
    // Load in the starting values (epsilon and pressure rom the EOS and a starting radius value for integration)
    // Be careful to make sure the single epsilon input is FROM the epsilon vector!
    MRLevaluator.input(epsilon, pressure, R_start, single_epsilon);
    // Interpolates the EoS in terms of the pseudo-enthalpy variable (h) 
    MRLevaluator.initialize_eos((gsl_interp_type *)gsl_interp_steffen);

    // Runs the TOV integrator
    MRLevaluator.TOV_Integrator(MRLevaluator.params.single_epsilon, &MRLevaluator.EoS);
    // Runs the TidalLove integrator
    MRLevaluator.TidalLove_Integrator(&MRLevaluator.fun);

    // NOTE THAT THE TIDAL LOVE NUMBER OUTPUT IS *DIMENSIONLESS*. I need to add functionality to redimensionalize it still!
    // Please see the QLIMR docs for more information: https://ce.musesframework.io/docs/modules/qlimr/contents/2_Physics_Overview.html.

    double M = MRLevaluator.NS_M;  // This is given in solar masses!
    double R = MRLevaluator.dimensionalize(MRLevaluator.NS_R, "km");  // This is given in km!
    double L = MRLevaluator.NS_Lbar;   // This has NO DIMENSIONS!!!
    
    return {M, R, L};
}
*/



template class IMRPhenomD_NRT_EOS<double>;
template class IMRPhenomD_NRT_EOS<adouble>;

// *********************** QLIMR TOV and λ̄ FUNCTIONALITY ***********************

// ----------------------------------------------------------------------------

QLIMR_params Input_QLIMR::params;
Input_QLIMR::Input_QLIMR() {}

// ------------------------- Input_QLIMR: Read yaml params --------------------

void Input_QLIMR::input(vector<double> epsilon, vector<double> pressure, double R_start, double single_epsilon) {
    params.R_start = adimensionalize(R_start, "km");

    // Check that the pressure and epsilon vectors are the same size
    if(pressure.size() == epsilon.size())
    {
      // Adimensionalize the pressure and epsilon
      for(int i = 0; i < epsilon.size(); i++) {
        epsilon[i] = adimensionalize(epsilon[i], "MeV/fm^3");
        pressure[i] = adimensionalize(pressure[i], "MeV/fm^3");
      }

      params.epsilon_col1 = epsilon;
      params.pressure_col2 = pressure;

      params.single_epsilon = adimensionalize(single_epsilon, "MeV/fm^3");
    }
    else{
      std::cout<<"Error: Epsilon and pressure vectors passed to LMR structure are not the same size. Check read EOS file."<<std::endl;
    }
}
// ----------------------------------------------------------------------------

// --------------- Adimensionalize conversion function ------------------------
double Input_QLIMR::adimensionalize(double value, string unit) {
    double factor;
  
    if (unit == "MeV/fm^3") {
      factor = 1.0 / 346933.783551;
    } else if (unit == "km") {
      factor = 1.0 / 1.47663;
    } else if (unit == "g/cm^3") {
      factor = 1.0 / 6.17625e+17;
    } else if (unit == "MeV"){
      factor = 8.96162e-61;   
    } else if (unit == "1/fm^3"){
      factor = 3.216297e54;  
    } else if (unit == "-") {
      factor = 1.0;
    } else {
      std::cout << "no unit match to adimensionalize" << std::endl;
      exit(0);
    }
    
    return factor * value;
  }

// --------------- Dimensionalize conversion function ------------------------
double Input_QLIMR::dimensionalize(double value, string unit) {
  double factor;

  if (unit == "MeV/fm^3") {
    factor = 346933.783551;
  } else if (unit == "km") {
    factor = 1.47663;
  } else if (unit == "g/cm^3") {
    factor = 6.17625e+17;
  } else if (unit == "MeV"){
    factor = 1.0/8.96162e-61;  
  } else if (unit == "1/fm^3"){
    factor = 1.0 / 3.216297e54;  
  } else if (unit == "Hz") {
    factor = (299792/1.47663); // c [km/s] / lsun [km]
  } else if (unit == "-") { 
    factor = 1.0;
  } else {
    std::cout << "no unit match to dimensionalize" << std::endl;
    exit(0);
  }

  return factor * value;
}

// ---------------------- Interpolation class methods -------------------------
void Interpolation::initialize(gsl_interp_type *interp_type, vector<double> x,
                               vector<double> y) {
  type = interp_type;
  // Size of the independent variable vector                            
  size = x.size();    

  // Allocating GSL interpolation accelerator
  acc = gsl_interp_accel_alloc();
  if (!acc) {
    std::cerr << "Error: Failed to allocate memory for accelerator." << std::endl;
    exit(EXIT_FAILURE);
  }

  // Allocating GSL spline of specified type and size
  spline = gsl_spline_alloc(type, size);

  // Initializing GSL spline with given data points (x, y) and size
  gsl_spline_init(spline, x.data(), y.data(), size);
}

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

// Method to release memory of spline and accelerator
void Interpolation::free() {
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
}
//-----------------------------------------------------------------------------

// ---------------------------- Defining dh/dε --------------------------------
double enthalpy_integrand(double epsilon, void *params) {

  double p = ((Interpolation *)params)->yofx(epsilon);
  double cs2 = ((Interpolation *)params)->dyofx(epsilon);

  return cs2 / (p + epsilon);
}
// ----------------------------------------------------------------------------

// -------------------------- Default EOS constructor --------------------------
EOSinterpolation EOS::EoS;
EOS::EOS(){}; // Default EOS constructor
// ----------------------------------------------------------------------------

// ----------------------- Parametric EOS constructor -------------------------
void EOS::initialize_eos(gsl_interp_type *type) {

  // Interpolate p(ε): initialize GSL interpolation spline
  EoS.p_of_e.initialize(type, params.epsilon_col1, params.pressure_col2);

  // Compute interpolated ε(h) and p(h) 
  calculate_eos_of_h(&params.epsilon_col1, type);

}

// ----------------------------------------------------------------------------

//------------------- Integration: Finding ε(h) and p(h) ----------------------
void EOS::calculate_eos_of_h(vector<double> *epsilon,
                                 gsl_interp_type *type) {
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

  for (size_t i = 1; i < epsilon->size() - 1; i++) {

    // Allocate memory for GSL integration workspace
    gsl_integration_workspace *w =
    gsl_integration_workspace_alloc(integration_workspace_size);

    // Perform adaptive quadrature integration using GSL library
    gsl_integration_qag(&F,
                       (*epsilon)[i],     // Lower integration limit
                       (*epsilon)[i + 1], // Upper integration limit
                       1e-9,              // Absolute tolerance
                       1e-9,              // Relative tolerance
                       1000,              // Maximal number of subintervals
                       6,                 // Integration method (Gauss-Kronrod)
                       w,                 // Work space
                       &delta_h,          // Estimated step size
                       &error             // Estimated integration error
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

  double p = pc - (2.0 / 3.0) * M_PI * pow(R, 2) * (epsilon_c + 3.0 * pc) *
                      (epsilon_c + pc);

  double e = epsilon_c - (2.0 / 3.0) * M_PI * pow(R, 2) *
                             (epsilon_c + 3.0 * pc) * (epsilon_c + pc) / Cs2;

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
void TOV::TOV_Integrator(double epsilon_c, EOSinterpolation *eos) {

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
  while (h > h_stop) {

    // Solve the ODE system using GSL with rkf45 method at each step
    int status = gsl_odeiv2_evolve_apply(e, c, s,
                 &sys, &h, h_stop, &h_istep, y);

    // Stop solving if something's wrong with the numerical integration
    if (status != GSL_SUCCESS) break;
      
    // Store h values to be used for finding ν(h) */
    h_sol.push_back(h);

    // Store radial distance (R) after each step
    R_sol.push_back(y[0]);

    // Store enclosed mass (M) after each step
    M_sol.push_back(y[1]);

    // Store energy density (ε) after each step
    e_sol.push_back(eos->e_of_h.yofx(h));

    if(e_sol.back() <= 0) {
      p_sol.push_back(0);
    }
    else {
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
  for (size_t i = 0; i < h_sol.size(); i++) {
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

// Chemical potential as a function of radius µ(R)
double TOV::mu_of_R(double R){

  double mu_Fe;
  double nu = fun.nu_of_R.yofx(R);
  double mu_of_R;

  // See Input_QLIMR::adimensionalize in A_Input.cpp
  mu_Fe = adimensionalize(930.54, "MeV" );
  mu_of_R = mu_Fe * sqrt((1.0 - 2.0 * (NS_M/NS_R))) * exp(-nu);

  return mu_of_R;
}

// Baryon number density as a function of radius n(R)
double TOV::n_of_R(double R){

  double mu_Fe;
  double p = fun.p_of_R.yofx(R);
  double e = fun.e_of_R.yofx(R);
  double nu = fun.nu_of_R.yofx(R);
  double mu_of_R;
  double n_of_R;

  mu_Fe = adimensionalize(930.54, "MeV" );
  mu_of_R = mu_Fe * sqrt((1.0 - 2.0 * (NS_M/NS_R))) * exp(-nu);

  n_of_R = (p + e) / mu_of_R;

  return n_of_R;
}

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

  double x0 = 1.0/R;
  double x1 = 4*M_PI;
  double x2 = std::pow(R, 3)*x1;
  double x3 = 2*M - R;
  double x4 = std::pow(R, 2);
  double x5 = M_PI*x4;
  double x6 = 5*e;

  double Y = y[0];

  f[0] =  -std::pow(Y, 2)*x0 + Y*(4*x5*(-e + p) + 1)/x3 + de*x2/(M + p*x2) 
          + x0*(4*std::pow(M, 2) + 4*M*R*(2*x5*(13*p + x6) - 3) 
          + std::pow(R, 4)*x1*(p*(16*p*x5 - 9) - x6) + 6*x4)/std::pow(x3, 2);

  return GSL_SUCCESS;
}

// Initial Conditions for Y
Second_Order::Initial_conditions_Y Second_Order::IC_Y()
{

  double R = params.R_start;

  // Defining initial step size for adaptative method
  IC_y.R_istep = R/100;

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
                      (4 * (NS_Y + 1.0) * C4 + (6.0 * NS_Y - 4.0) * C3 +
                       (26.0 - 22.0 * NS_Y) * C2 +
                       3.0 * (5.0 * NS_Y - 8.0) * C - 3.0 * NS_Y + 6.0) -
                  3.0 * k * (2 * C * (NS_Y - 1.0) - NS_Y + 2.0) *
                      log(1.0 / (1.0 - 2.0 * C))));

  // Dimensionless tidal love number NS_Lbar: λ̄ [-]
  NS_Lbar = (2.0 / 3.0) * NS_k2 * (1 / C5);

  // Initialize GSL interpolation spline for Y(r) solution
  fun->Y_of_R.initialize((gsl_interp_type *)gsl_interp_steffen, R_sol, Y_sol);
  
}
