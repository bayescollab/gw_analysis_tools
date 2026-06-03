#include "gwat/util.h"
#include "gwat/fisher.h"
#include "gwat/detector_util.h"
#include "gwat/waveform_util.h"
#include "gwat/ortho_basis.h"
#include "gwat/pn_waveform_util.h"
#include "gwat/ppE_utilities.h"
#include <iostream>
#include <iomanip>

#include <fstream>

#include <chrono>

void set_arguments(double Tsignal, double deltaF, int length, double fmin,
       	double *& frequencies,
       	gen_params_base<double>*& parameters) {
  
  
  parameters->mass1 = 18.0;
	parameters->mass2 = 14.0;
  
  //Spins
	parameters->spin1[2] = .0;
	parameters->spin2[2] = .0;
	parameters->spin1[1] = .0;
	parameters->spin2[1] = .0;
	parameters->spin1[0] = .0;
	parameters->spin2[0] = .0;

	//Extrinsic
	parameters->Luminosity_Distance = 600;
	parameters->psi = .2;
  parameters->RA = 3.42;
	parameters->DEC = -.37;
	parameters->incl_angle = 2.532207345558998;
	double gps = 1187008882.4;
	parameters->gmst = gps_to_GMST_radian(gps);

	//Trivial constants
	parameters->phiRef = 2.;
	parameters->f_ref = 20.;

	//Specific flags -- probably don't modify
	parameters->equatorial_orientation = false;
	parameters->horizon_coord = false;
	parameters->shift_time = true;
	parameters->shift_phase = true;
  

  double T_merger= Tsignal*3./4.;
  
  parameters->tc=Tsignal-T_merger;
	
  
  for (int j = 0; j < length; j++) {
    frequencies[j] = fmin + j * deltaF;
  }
}
  

int main(int argc, char* argv[])
{
  std::cout<<"TESTING FISHER CALCULATE DERIVATIVES FUNCTION"<<std::endl;
  
  //variable sfrom parameters->cpp used to get length
  double fmin = 10;
  double fmax = 2048;
  double Tsignal = 4;
  double deltaF = 1./Tsignal;
  

  int length = (int)((fmax-fmin)/deltaF);
  int dimension = 11;
  
  // gen_method = injection_method
  string  gen_method = "IMRPhenomD";
  string detector = "Hanford";
  string reference_detector = "Hanford";


  gen_params_base<double> *parameters = new gen_params_base<double>();
  double                  *frequencies = new double[length];
  
  // initialize and define response_deriv arrays
  std::complex<double>    **response_deriv_original = new std::complex<double>*[dimension];
  std::complex<double>    **response_deriv_new = new std::complex<double>*[dimension];
	for (int i = 0; i < dimension; i++){
		response_deriv_original[i] = new std::complex<double>[length];
    response_deriv_new[i] = new std::complex<double>[length];
    for (int j = 0; j < length; j++) {
      response_deriv_original[i][j] = std::complex<double>(0.0,0.0);
      response_deriv_new[i][j] = std::complex<double>(0.0,0.0);
    }
	}
  

  // Running Tests
  int order = 4; // CHANGE TO TEST DIFFERENT ORDER
  parameters->sky_average = true; // CHANGE TO TEST 

  auto some_dum_thing = std::chrono::steady_clock::now();
  auto total_time = some_dum_thing - some_dum_thing; // dum way to initailize to zero

  
  float iterations = 1000;
  for (int i = 0; i < iterations; i++) {
    
    // output parameter: response_deriv_new
    set_arguments(Tsignal, deltaF, length, fmin, frequencies, parameters);

    //---------------------------------------------------------------------------------------------------------
    // INSERT NEW calculate_derivatives() FUNCTION HERE vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    //---------------------------------------------------------------------------------------------------------
    // make sure to use "response_deriv_new" as response_deriv parameter - all others use same as original
    auto start_time = std::chrono::steady_clock::now();
    // only run one or the other
    // calculate_derivatives(response_deriv_original, frequencies, length, dimension, detector, reference_detector, gen_method, parameters, order);
    calculate_derivatives_new(response_deriv_new, frequencies, length, dimension, detector, reference_detector, gen_method, parameters, order);

    //---------------------------------------------------------------------------------------------------------
    //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    //---------------------------------------------------------------------------------------------------------
    // std::cout<<"finished new\n"<<std::endl;
    auto end_time = std::chrono::steady_clock::now();
    // std::cout<<std::chrono::duration_cast<std::chrono::microseconds>(end_time-start_time).count()<<std::endl;
    total_time += (end_time - start_time);
  }
  auto avg_duration =std::chrono::duration_cast<std::chrono::microseconds>(total_time / iterations);
  std::cout<<"Average Duration: " <<avg_duration.count()<<" microseconds"<<std::endl;


  // deallocate things
  delete parameters;
  delete[] frequencies;
  for (int i = 0; i < dimension; i++){
    delete[] response_deriv_original[i];
    delete[] response_deriv_new[i];
  }
  delete[] response_deriv_original;
  delete[] response_deriv_new;

  return 0;
}
