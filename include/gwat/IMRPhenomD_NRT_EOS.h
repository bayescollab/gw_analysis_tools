#ifndef IMRPHENOMD_NRT_EOS
#define IMRPHENOMD_NRT_EOS
#include "IMRPhenomD_NRT.h"
#include "util.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

using std::string;
using std::vector;

/*! \file
 */

/*! Class that extends the IMRPhenomD_NRT waveform to sample directly on
 * equation of state (EOS) parameters.
 */

template<class T>
class IMRPhenomD_NRT_EOS: public IMRPhenomD_NRT<T>
{
public:
  // Function to build bump in cs2
  virtual void build_cs2_one_bump(source_parameters<T> *sp);

  // Function to convert cs2 to p(epsilon)
  virtual void cs2_to_eos_convert(); // Empty for now, "master function" that will be loaded with C++ conversion of Jaki's code

  // Function to calculate observable variables from the EOS
  virtual void get_m_love(gen_params* params); // Empty for now, "master function" that will be loaded with QLIMR functionality

  // Function to interface with


  /* Leaving this here for now in case we need it for future development - PLEASE DELETE LATER IF NOT USED!
  
  // Inherited from IMRPhenomD_NRT which inherited from IMRPhenomD
  virtual int construct_waveform(T *frequences, int length, std::complex<T> *waveform, source)parameters<T> *params);
  */

};

// ****************************************************************************
// Structure to store input parameters
struct QLIMR_params {
    double R_start;
    double single_epsilon;
    vector<double> epsilon_col1;
    vector<double> pressure_col2;
};

class Input_QLIMR {
    public:
    // Static structure input parameters to be passed by inheritance
    static QLIMR_params params;
    
    // Default and parametric constructors
    Input_QLIMR();

    // Initialize object to read source EOS file
    void input(vector<double> epsilon, vector<double> pressure, double R_start, double single_epsilon);
    
    // Unit conversion functions
    double adimensionalize(double value, string unit);
    double dimensionalize(double value, string unit);
    
};

// ****************************************************************************
class Interpolation : public Input_QLIMR {
public:

  // GSL-type variables for interpolation
  gsl_interp_accel *acc;
  gsl_spline       *spline;
  gsl_interp_type  *type;

  // Size of data to be interpolated
  size_t size;

  // Initializing GSL spline for interpolation
  void initialize(gsl_interp_type *type, vector<double> x, vector<double> y);

  // Evaluate interpolated function at point x
  double yofx(double x);

  // Evaluate derivative of interpolated function at point x
  double dyofx(double x);

  // Free GSL spline and accelerator memory
  void free();
};
// ****************************************************************************

struct EOSinterpolation {

  // Interpolation type objects:

  Interpolation p_of_e; // Pressure as a function of energy density
  Interpolation h_of_e; // Pseudo-enthalpy as a function of energy density
  Interpolation h_of_p; // Enthalpy as a function of pressure
  Interpolation e_of_h; // Energy density as a function of pseudo-enthalpy
  Interpolation p_of_h; // Pressure as a function of pseudo-enthalpy

  // Vectors to store EoS data columns

  vector<double> e_vec;  // Energy density vector
  vector<double> p_vec;  // Pressure vector
  vector<double> h_vec;  // Enthalpy vector
};

// ****************************************************************************
class EOS : public Interpolation {
public:
  // Static structure for EoS to be used along the code 
  static EOSinterpolation EoS;

  // Default constructor
  EOS();

  // Parametric constructor: 
  void initialize_eos(gsl_interp_type *type);

  // Method function to compute EoS in terms of pseudo-enthalpy (h)
  void calculate_eos_of_h(vector<double> *epsilon, gsl_interp_type *type);
};
// ****************************************************************************

struct Local_functions {
    // Zeroth order local functions
    Interpolation M_of_R;
    Interpolation p_of_R;
    Interpolation e_of_R;
    Interpolation nu_of_R;
  
    // Second order local functions
    Interpolation Y_of_R;
  
    // ######### Integration constants ########
    // First order
    double NS_Omega;
  
  };
  
  //--------------------------------- CLASS TOV ---------------------------------
  class TOV : public EOS {
  public:
  
    // Declaring structure variable fun of type Local_functions
    Local_functions fun;
  
    // Declaring neutron star mass and radius!
    double NS_R;
    double NS_M;
  
    // Structure to hold initial conditions for TOV integrator
    struct Initial_conditions_TOV {
      double R_start;
      double M_start; 
      double h_start;
      double h_istep;
    };
  
    // Instance of Initial_conditions_TOV to store initial conditions
    Initial_conditions_TOV IC_tov;
  
    // Method to calculate initial conditions based on central energy density
    Initial_conditions_TOV IC_TOV(double epsilon_c);
  
    // TOV integrator method using pseudo-enthalpy (h) formulation
    void TOV_Integrator(double epsilon_c, EOSinterpolation *eos);
  
    // Method to calculate chemical potential as a function of radius µ(R)
    double mu_of_R(double r);
  
    // Method to calculate baryon number density as a function of radius n(R)
    double n_of_R(double r);
  
  };
  //-----------------------------------------------------------------------------

  class Second_Order : public TOV {
    public:
      double NS_k2;  
      double NS_Lbar;
      
      // ########################### Tidal Love number λ̄ #########################
    
      struct Initial_conditions_Y {
        double Y_start;
        double R_istep;
      };
    
      Initial_conditions_Y IC_y;
      Initial_conditions_Y IC_Y();
    
      // Integrator to obtain tidal love number: λ̄
      void TidalLove_Integrator(Local_functions *fun);
    
    };

#endif
