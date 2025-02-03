#ifndef ERROR_H
#define ERROR_H
#include "util.h"
#include <string>

/*! \file
 */

using namespace std; 

void calculate_systematic_error(
    double *frequency,
	std::complex<double>* h_model,
	std::complex<double>* h_true, 
	int length,/**< if 0, standard frequency range for the detector is used*/ 
	std::string generation_method, 
	std::string* detectors, 
	std::string reference_detector,  
	double* output,/**< double [dimension][dimension]*/
	int dimension, 
	gen_params_base<double> *parameters,
	int order,/**< Order of the numerical derivative (2 or 4)**/
	//double *parameters,
	double *noise,
	int Ndetect
);

double calculate_sys_err_elements(
	double *frequency,
	std::complex<double>* h_model,
	std::complex<double>* h_true,
	std::complex<double> **dh,
	double** fisher_inverse,
	double *psd, 
	int length,
	int dimension,
	int a,
	int b
);

void calculate_statistical_error(
    double *frequency,
	//std::complex<double>* h_model,
	//std::complex<double>* h_true, 
	int length,/**< if 0, standard frequency range for the detector is used*/ 
	std::string generation_method, 
	std::string* detectors, 
	std::string reference_detector,  
	double* output,/**< double [dimension][dimension]*/
	int dimension, 
	gen_params_base<double> *parameters,
	int order,/**< Order of the numerical derivative (2 or 4)**/
	//double *parameters,
	double **noise,
	int Ndetect
);


#endif
