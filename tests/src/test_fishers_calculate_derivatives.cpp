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
  int order = 2; // CHANGE TO TEST DIFFERENT ORDER
  parameters->sky_average = false; // CHANGE TO TEST 


  // Running original calculate_derivatives()
  // output parameter: response_deriv_original
  set_arguments(Tsignal, deltaF, length, fmin, frequencies, parameters);
  std::cout<<"Running original..."<<std::endl;
  calculate_derivatives(response_deriv_original, frequencies, length, dimension, detector, reference_detector, gen_method, parameters, order);
  std::cout<<"finished original\n"<<std::endl;

  // Running new calculate_derivatives()
  // output parameter: response_deriv_new
  set_arguments(Tsignal, deltaF, length, fmin, frequencies, parameters);
  std::cout<<"Running new..."<<std::endl;

  //---------------------------------------------------------------------------------------------------------
  // INSERT NEW calculate_derivatives() FUNCTION HERE vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  //---------------------------------------------------------------------------------------------------------
  // make sure to use "response_deriv_new" as response_deriv parameter - all others use same as original
  calculate_derivatives_new(response_deriv_new, frequencies, length, dimension, detector, reference_detector, gen_method, parameters, order);

  //---------------------------------------------------------------------------------------------------------
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  //---------------------------------------------------------------------------------------------------------
  std::cout<<"finished new\n"<<std::endl;

  // open files to output response_deriv content
  std::ofstream output_original("output_original_calculate_derivatives.txt");
  std::ofstream output_new("output_new_calculate_derivatives.txt");

  if (!output_original.is_open() || !output_new.is_open()) {
    std::cout<<"file failed to open sad gonna return now"<<std::endl;;
    return 1;
}

  // Compare original and new response_deriv elements and outputting to files
  bool flag = true; // true = all values of elements of arrays are same thus far
  double TOLERANCE = 1e-8; // tolerance is percentage error new to original

  std::cout<<"Comparing Elements of original and new..."<<std::endl;
  for (int i = 0; i < dimension; i++) {
    for (int j = 0; j < length; j++) {
      // compare real components
      if (abs((response_deriv_new[i][j].real() - response_deriv_original[i][j].real())/response_deriv_original[i][j].real()) > TOLERANCE) {
        flag = false;
        std::cout<<"Real component not equal at index ("<<i<<", "<<j<<"), real component percentage error: "<< 
        abs((response_deriv_new[i][j].real() - response_deriv_original[i][j].real())/response_deriv_original[i][j].real())*100.<<"%"
        <<"\nOriginal: " << response_deriv_original[i][j].real() <<"\nNew: " << response_deriv_new[i][j].real()<<std::endl;

      }
      if (abs((response_deriv_new[i][j].imag() - response_deriv_original[i][j].imag())/response_deriv_original[i][j].imag()) > TOLERANCE){
        flag = false;
        std::cout<<"Imaginary component not equal at index ("<<i<<", "<<j<<"), imag component percentage error: "<< 
        abs((response_deriv_new[i][j].imag() - response_deriv_original[i][j].imag())/response_deriv_original[i][j].imag())*100.<<"%"
        <<"\nOriginal: " << response_deriv_original[i][j].imag() <<"\nNew: " << response_deriv_new[i][j].imag()<<"\n"<<std::endl;

      }
      output_original<<(response_deriv_original[i][j].real()) << "+" << (response_deriv_original[i][j].imag())<<"j";
      output_new <<(response_deriv_new[i][j].real()) << "+" << (response_deriv_new[i][j].imag())<<"j";
      if (j < length-1) {
        output_original<<", ";
        output_new<<", ";
      }

    }
    output_original<<std::endl;
    output_new<<std::endl;
  }
  
  if (flag) {
    std::cout<<"Passed! All elements are same!"<<std::endl;
  }
  else {
    std::cout<<"Not passed! Not all elements are same!"<<std::endl;
  }

  // std::cout<<response_deriv_original[0][0]<<std::endl;
  // std::cout<<response_deriv_new[0][0]<<std::endl;
  

  output_original.close();
  output_new.close();




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
