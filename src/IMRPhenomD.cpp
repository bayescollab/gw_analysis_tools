#include "IMRPhenomD.h"
#include "QNM_data.h"
//#include "general_parameter_structures.h"
#include "util.h"
#include <math.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <adolc/adouble.h>
#include <adolc/taping.h>
#include <adolc/drivers/drivers.h>
#include <typeinfo>
#include <omp.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
using namespace std;

#ifndef _OPENMP
#define omp ignore
#endif



//double log_64 = 1.80617997398;//log base 10...
double log_64 = 4.15888308336;
/*! \file 
 * File that includes all the low level functions that go into constructing the waveform
 *
 * Matches LALsuite -- 2019_09_25
 *
 * Two versions of QNM fitting -- one matches LALSuites' implementation, and is probably more accuate, but involves GSL fitting (ie not ADOLC friendly). The other version ensures ADOLC and numerical versions match -- for testing -- just comment out one option if testing (only fRD and fdamp)
 */
//#####################################################################################

template <class T>
void IMRPhenomD<T>::fisher_calculation_sky_averaged(double *frequency, 
			int length, 
			gen_params *parameters,
			double **amplitude_deriv, 
			double **phase_deriv, 
			double *amplitude, 
			int *amp_tapes, 
			int *phase_tapes
			)
{
	IMRPhenomD<double> modeld;
	int dimension = 7;
	//populate model
	source_parameters<double> input_params;
	//input_params = source_parameters<double>::populate_source_parameters(parameters);
	input_params.populate_source_parameters(parameters);
	//input_params.shift_time = parameters->shift_time;
	input_params.shift_time = false;
	//Need the splitting frequency	
	lambda_parameters<double> lambda, *lambda_ptr;
	modeld.assign_lambda_param(&input_params, &lambda);
	modeld.post_merger_variables(&input_params);
	input_params.f1 = 0.014/(input_params.M);
	input_params.f3 = modeld.fpeak(&input_params, &lambda);
	input_params.f1_phase = 0.018/(input_params.M);
	input_params.f2_phase = input_params.fRD/2.;
	//###########################################################################
	//populate derivative
	modeld.construct_amplitude_derivative(frequency, 
			length,
			dimension, 
			amplitude_deriv,
			&input_params,
			amp_tapes
			);
	
	//PROBLEM SECTION
	modeld.construct_phase_derivative(frequency, 
			length,
			dimension, 
			phase_deriv,
			&input_params,
			phase_tapes
			);
	modeld.construct_amplitude(frequency, length, amplitude, &input_params);
	//LOG FACTORS for A0, chirpmass, and eta)
	for (int i = 0;i <length; i++)
	{
		amplitude_deriv[0][i] = (input_params.A0)*amplitude_deriv[0][i];
		amplitude_deriv[3][i] = (input_params.chirpmass)*amplitude_deriv[3][i];
		amplitude_deriv[4][i] = (input_params.eta)*amplitude_deriv[4][i];
		phase_deriv[0][i] = (input_params.A0)*phase_deriv[0][i];
		phase_deriv[3][i] = (input_params.chirpmass)*phase_deriv[3][i];
		phase_deriv[4][i] = (input_params.eta)*phase_deriv[4][i];
	}

}

/*! \brief Construct the derivative of the amplitude for a given source evaluated by the given frequency
 *
 * Order of output: dh/d \theta : \theta \el {A0,tc, phic, chirp mass, eta, symmetric spin, antisymmetric spin}
 * 
 */
template <class T>
void IMRPhenomD<T>::construct_amplitude_derivative(double *frequencies, /**< input array of frequency*/
				int length,/**< length of the frequency array*/
				int dimension,/** < dimension of the fisher*/
				double **amplitude_derivative,/**< output array for all the derivatives double[dimension][length]*/
				source_parameters<double> *input_params,/**< Source parameters structure for the source*/
				int *tapes /**<int array of tape ids, if NULL, these will be calculated*/
				)
{
	//if(tapes==NULL)
	if(!tapes)
	{
		int tape_temp[3];
		tape_temp[0]=10;
		tape_temp[1]=11;
		tape_temp[2]=12;
		this->amplitude_tape(input_params,tape_temp);
		tapes = tape_temp;
	}
	
	double evaluate_params[dimension +1] ;
	double grad[dimension+1];
	evaluate_params[1] = input_params-> A0;
	evaluate_params[2] = input_params-> phiRef;
	evaluate_params[3] = input_params-> tc;
	evaluate_params[4] = input_params-> chirpmass;
	evaluate_params[5] = input_params-> eta;
	evaluate_params[6] = input_params-> chi_s;
	evaluate_params[7] = input_params-> chi_a;
	for (int i=0;i<length; i++)
	{
		evaluate_params[0] = frequencies[i];
		if(evaluate_params[0]>.2/input_params->M){
			for(int j=0;j<dimension;j++)
			{
				amplitude_derivative[j][i] = 0;
			} 

		}
		else{	
			if(evaluate_params[0]<input_params->f1)
			{
				gradient(tapes[0],dimension+1, evaluate_params, grad);
				for(int j=0;j<dimension;j++)
				{
					amplitude_derivative[j][i] = grad[j+1];
				} 
			}
			else if(evaluate_params[0] <input_params->f3)
			{
				gradient(tapes[1],dimension+1, evaluate_params, grad);
				for(int j=0;j<dimension;j++)
				{
					amplitude_derivative[j][i] = grad[j+1];
				} 
			}
			else
			{
				gradient(tapes[2],dimension+1, evaluate_params, grad);
				for(int j=0;j<dimension;j++)
				{
					amplitude_derivative[j][i] = grad[j+1];
				} 
			}
		}
	}
	
}
/*! \brief Construct the derivative of the phase for a given source evaluated by the given frequency
 *
 * Order of output: dh/d \theta : \theta \el {A0,tc, phic, chirp mass, eta, symmetric spin, antisymmetric spin}
 * 
 */
template <class T>
void IMRPhenomD<T>::construct_phase_derivative(double *frequencies, /**< input array of frequency*/
				int length,/**< length of the frequency array*/
				int dimension,/** < dimension of the fisher*/
				double **phase_derivative,/**< output array for all the derivatives double[dimension][length]*/
				source_parameters<double> *input_params,/**< Source parameters structure for the source*/
				int *tapes /**<int array of tape ids, if NULL, these will be calculated*/
				)
{
	//if(tapes==NULL)
	if(!tapes)
	{
		int tape_temp[3];
		tape_temp[0]=13;
		tape_temp[1]=14;
		tape_temp[2]=15;
		this->phase_tape(input_params,tape_temp);
		tapes = tape_temp;
	}
	
	double evaluate_params[dimension +1] ;
	double grad[dimension+1];
	evaluate_params[1] = input_params-> A0;
	evaluate_params[2] = input_params-> phiRef;
	evaluate_params[3] = input_params-> tc;
	evaluate_params[4] = input_params-> chirpmass;
	evaluate_params[5] = input_params-> eta;
	evaluate_params[6] = input_params-> chi_s;
	evaluate_params[7] = input_params-> chi_a;
	for (int i=0;i<length; i++)
	{
		evaluate_params[0] = frequencies[i];
		if(evaluate_params[0]>.2/input_params->M){
			for(int j=0;j<dimension;j++)
			{
				phase_derivative[j][i] = 0;
			} 

		}
		else{	
		if(evaluate_params[0]<input_params->f1_phase)
		{
			gradient(tapes[0],dimension+1, evaluate_params, grad);
			for(int j=0;j<dimension;j++)
			{
				phase_derivative[j][i] = grad[j+1];
			} 
		}
		else if(evaluate_params[0] <input_params->f2_phase)
		{
			gradient(tapes[1],dimension+1, evaluate_params, grad);
			for(int j=0;j<dimension;j++)
			{
				phase_derivative[j][i] = grad[j+1];
			} 
		}
		else
		{
			gradient(tapes[2],dimension+1, evaluate_params, grad);
			for(int j=0;j<dimension;j++)
			{
				phase_derivative[j][i] = grad[j+1];
			} 
		}
		}
	}
	
}
/*! \brief Creates the tapes for derivatives of the amplitude
 *
 * For efficiency in long runs of large sets of fishers, the tapes can be precomputed and reused
 */
template <class T>
void IMRPhenomD<T>::amplitude_tape(source_parameters<double> *input_params, /**< source parameters structure of the desired source*/
				int *tape /**<tape ids*/
				)
{
	//Find splitting frequencies between tapes
	IMRPhenomD<double> modeld;
	lambda_parameters<double> lambda, *lambda_ptr;
	modeld.assign_lambda_param(input_params, &lambda);
	modeld.post_merger_variables(input_params);
	input_params->f1 = 0.014/(input_params->M);
	input_params->f3 = modeld.fpeak(input_params, &lambda);
	double f2 = (input_params->f1 + input_params->f3)/2.;
	double freqs[3] = {input_params->f1 * .9, f2, input_params->f3 * 1.1};
	
	for (int i =0; i<3;i ++)	
	{
		trace_on(tape[i]);
		double amp_out;

		adouble freq  ;
		adouble parameters[7];
		//adouble amp ;
		IMRPhenomD<adouble> model;
		adouble new_params[7];

		freq <<= freqs[i];
		parameters[0] <<= input_params-> A0;
		parameters[1] <<= input_params-> phiRef;
		parameters[2] <<= input_params-> tc;
		parameters[3] <<= input_params-> chirpmass;
		parameters[4] <<= input_params-> eta;
		parameters[5] <<= input_params-> chi_s;
		parameters[6] <<= input_params-> chi_a;
		//parameters[3] = parameters[3]/MSOL_SEC;
		model.change_parameter_basis(parameters, new_params,input_params->sky_average);

		adouble spin1vec[3] = {0,0,new_params[3]};
		adouble spin2vec[3] = {0,0,new_params[4]};
		source_parameters<adouble> intermediate_params;
		//intermediate_params = intermediate_params.populate_source_parameters_old(new_params[0]/MSOL_SEC,
		intermediate_params.populate_source_parameters_old(new_params[0]/MSOL_SEC,
				new_params[1]/MSOL_SEC,new_params[2]/MPC_SEC,spin1vec,spin2vec,new_params[5],
				new_params[6],input_params->sky_average);
		intermediate_params.shift_time = input_params->shift_time;

		adouble freqs[1] = {freq};
		adouble amp_temp[1];
		model.construct_amplitude(freqs, 1,  amp_temp, &intermediate_params);
		amp_temp[0]>>=amp_out;
		trace_off();
	}
	
}

/*! \brief Creates the tapes for derivatives of phase
 *
 * For efficiency in long runs of large sets of fishers, the tapes can be precomputed and reused
 */
template <class T>
void IMRPhenomD<T>::phase_tape(source_parameters<double> *input_params, /**< source parameters structure of the desired source*/
				int *tape /**<tape ids*/
				)
{
	IMRPhenomD<double> modeld;
	//Find splitting frequencies between tapes
	lambda_parameters<double> lambda, *lambda_ptr;
	modeld.assign_lambda_param(input_params, &lambda);
	modeld.post_merger_variables(input_params);
	input_params->f1_phase = 0.018/(input_params->M);
	input_params->f2_phase = input_params->fRD/2.;
	double f2 = (input_params->f1_phase + input_params->f2_phase)/2.;
	double freqs[3] = {input_params->f1_phase * .9, f2, input_params->f2_phase * 1.1};
	
	for (int i =0; i<3;i ++)	
	{
		trace_on(tape[i]);
		double phase_out;

		adouble freq  ;
		adouble parameters[7];
		adouble phase ;
		IMRPhenomD<adouble> model;

		freq <<= freqs[i];
		parameters[0] <<= input_params-> A0;
		parameters[1] <<= input_params-> phiRef;
		parameters[2] <<= input_params-> tc;
		parameters[3] <<= input_params-> chirpmass;
		parameters[4] <<= input_params-> eta;
		parameters[5] <<= input_params-> chi_s;
		parameters[6] <<= input_params-> chi_a;
		//parameters[3] = parameters[3]/MSOL_SEC;
		adouble new_params[7];
		model.change_parameter_basis(parameters, new_params,input_params->sky_average);
		source_parameters<adouble> intermediate_params;
		adouble spin1vec[3] = {0,0,new_params[3]};
		adouble spin2vec[3] = {0,0,new_params[4]};
		
		//intermediate_params = intermediate_params.populate_source_parameters_old(new_params[0]/MSOL_SEC,
		intermediate_params.populate_source_parameters_old(new_params[0]/MSOL_SEC,
				new_params[1]/MSOL_SEC,new_params[2]/MPC_SEC,spin1vec,spin2vec,new_params[5],
				new_params[6], input_params->sky_average);
		intermediate_params.shift_time = input_params->shift_time;
		adouble freqs[1] = {freq};
		adouble phase_temp[1];
		model.construct_phase(freqs, 1,  phase_temp, &intermediate_params);
		phase_temp[0]>>=phase_out;
		trace_off();
	}
	
}
/*! \brief Convience method to change parameter basis between common Fisher parameters and the intrinsic parameters of IMRPhenomD
 *
 * Takes input array of old parameters and ouputs array of transformed parameters
 */
template <class T>
void IMRPhenomD<T>::change_parameter_basis(T *old_param, /**< array of old params, order {A0, tc, phic, chirpmass, eta, spin1, spin2}*/
					T *new_param, /**< output new array: order {m1,m2,DL, spin1,spin2,phic,tc}*/
					bool sky_average
					)
{
	T m1, m2, DL, spin1,spin2;
	m1 = calculate_mass1(old_param[3],old_param[4]);
	m2 = calculate_mass2(old_param[3],old_param[4]);
	//T chirpmass_sec = old_param[3]*MSOL_SEC;
	//DL =1/(old_param[0] * sqrt(30/M_PI) /(chirpmass_sec*chirpmass_sec)/
	//		pow(M_PI*chirpmass_sec,-7./6))/MPC_SEC  ;
	//#############################################
	//Sky averaged
	if(sky_average){
		DL =1/(old_param[0] * sqrt(30/M_PI) /(old_param[3]*old_param[3])/
			pow(M_PI*old_param[3],-7./6))  ;
	}
	//#############################################
	//Not sky averaged
	else{
		DL =1/(old_param[0] * sqrt(192./(40*M_PI)) /(old_param[3]*old_param[3])/
				pow(M_PI*old_param[3],-7./6))  ;
	}
	spin1 = old_param[5]+old_param[6];
	spin2 = old_param[5]-old_param[6];
	T spinvec1[3] ={0.,0.,spin1};
	T spinvec2[3] ={0.,0.,spin2};
	new_param[0] = m1;
	new_param[1] = m2;
	new_param[2] = DL;
	new_param[3] = spin1;
	new_param[4] = spin2;
	new_param[5] = old_param[1];
	new_param[6] = old_param[2];
}


/*! \brief Constructs the waveform as outlined by 
 *
 * arguments:
 * 	array of frequencies, length of that array, a complex array for the output waveform, and a source_parameters structure
 */
template <class T>
int IMRPhenomD<T>::construct_waveform(T *frequencies, /**< T array of frequencies the waveform is to be evaluated at*/
				int length, /**< integer length of the array of frequencies and the waveform*/	
				std::complex<T> *waveform,/**< complex T array for the waveform to be output*/ 
				source_parameters<T> *params /*Structure of source parameters to be initialized before computation*/
				)
{
	T M = params-> M;
	T chirpmass = params->chirpmass;
	T DL = params->DL;
	lambda_parameters<T> lambda, *lambda_ptr;
	this->assign_lambda_param(params, &lambda);

	/*Initialize the post merger quantities*/
	this->post_merger_variables(params);
	params->f1_phase = 0.018/(params->M);
	params->f2_phase = params->fRD/2.;

	params->f1 = 0.014/(params->M);
	params->f3 = this->fpeak(params, &lambda);
	
	useful_powers<T> pows;
	this->precalc_powers_PI(&pows);

	T deltas[6];
	T pn_amp_coeffs[7];
	T pn_phase_coeffs[12];

	this->assign_pn_amplitude_coeff(params, pn_amp_coeffs);
	this->assign_static_pn_phase_coeff(params, pn_phase_coeffs);	

	this->amp_connection_coeffs(params,&lambda,pn_amp_coeffs,deltas);
	this->phase_connection_coefficients(params,&lambda,pn_phase_coeffs);

	//################################################################
	//Calculate phase and coalescence time variables
	T phic, f_ref, tc, phi_shift, tc_shift;
	//If phic is unspecified - use f_ref and phiRef
	if(params->shift_phase ){
		f_ref = params->f_ref;
		precalc_powers_ins(f_ref, M, &pows);
		phi_shift = (this->build_phase(f_ref,&lambda,params,&pows,pn_phase_coeffs));
		phic = 2*params->phiRef + phi_shift;
	}
	//If phic is specified, ignore f_ref phiRef and use phic
	else{
		f_ref = 0;
		phic = params->phiRef;
	}
	
	//Assign shift: first shift so coalescence happens at t=0, then shift from there according to tc
	//This aligns more with the physical meaning of tc, but the phase is NO LONGER just
	if(params->shift_time){
		T alpha1_offset = assign_lambda_param_element(params,14);
		tc_shift = this->Dphase_mr(params->f3, params, &lambda)+(-lambda.alpha[1]+alpha1_offset)*params->M/params->eta;

		//tc_shift = this->Dphase_mr(params->f3, params, &lambda);
		//debugger_print(__FILE__,__LINE__,"SHIFTING DOES NOT MATCH LALSUITE");
	}
	else{
		tc_shift=0;
	}
	//tc = 2*M_PI*params->tc - tc_shift;
	tc = 2*M_PI*params->tc + tc_shift;

	//T A0 = sqrt(M_PI/30)*chirpmass*chirpmass/DL * pow(M_PI*chirpmass,-7./6);
	T A0 = params->A0* pow(M,7./6.);

	T f;
	std::complex<T> amp, phase;
	std::complex<T> i;
	i = std::complex<T> (0,1.);
	T fcut = .2/M; //Cutoff frequency for IMRPhenomD - all higher frequencies return 0
	for (size_t j =0; j< length; j++)
	{
		f = frequencies[j];
		if(f>fcut){
			amp = 0.0;
			waveform[j] = 0.0;
		}
		else{	
			if (f<params->f1_phase)
			{
				precalc_powers_ins(f, M, &pows);
			}
			else
			{
				pows.MFsixth= pow(M*f,1./6.);	
				pows.MF7sixth= pow_int(pows.MFsixth,7);//*pows.MFsixth*pows.MFsixth*pows.MFsixth*pows.MFsixth*pows.MFsixth*pows.MFsixth;
			}
			amp = (A0 * this->build_amp(f,&lambda,params,&pows,pn_amp_coeffs,deltas));
			phase = (this->build_phase(f,&lambda,params,&pows,pn_phase_coeffs));
			//phase +=   (T)(tc*(f-f_ref) - phic);
			//debugger_print(__FILE__,__LINE__,std::real(phase));
			phase -=   (T)(tc*(f-f_ref) + phic);
			waveform[j] = amp * std::exp(-i * phase);
		}

	}
	
	return 1;
}

/*! \brief overloaded method to evaluate the waveform for one frequency instead of an array
 */
template <class T>
std::complex<T> IMRPhenomD<T>::construct_waveform(T frequency, /**< T array of frequencies the waveform is to be evaluated at*/
				source_parameters<T> *params /*Structure of source parameters to be initialized before computation*/
				)
{
	T M = params-> M;
	T chirpmass = params->chirpmass;
	T DL = params->DL;
	lambda_parameters<T> lambda, *lambda_ptr;
	this->assign_lambda_param(params, &lambda);

	/*Initialize the post merger quantities*/
	this->post_merger_variables(params);

	params->f1_phase = 0.018/(params->M);
	params->f2_phase = params->fRD/2.;

	params->f1 = 0.014/(params->M);
	params->f3 = this->fpeak(params, &lambda);

	useful_powers<T> pows;
	this->precalc_powers_PI(&pows);

	T deltas[6];
	T pn_amp_coeffs[7];
	T pn_phase_coeffs[12];

	this->assign_pn_amplitude_coeff(params, pn_amp_coeffs);
	this->assign_static_pn_phase_coeff(params, pn_phase_coeffs);	

	this->amp_connection_coeffs(params,&lambda,pn_amp_coeffs,deltas);
	this->phase_connection_coefficients(params,&lambda,pn_phase_coeffs);

	//################################################################
	//Calculate phase and coalescence time variables
	T phic, f_ref, tc, phi_shift, tc_shift;
	//If phic is unspecified - use f_ref and phiRef
	if(params->shift_phase ){
		f_ref = params->f_ref;
		precalc_powers_ins(f_ref, M, &pows);
		phi_shift = (this->build_phase(f_ref,&lambda,params,&pows,pn_phase_coeffs));
		phic = 2*params->phiRef + phi_shift;
	}
	//If phic is specified, ignore f_ref phiRef and use phic
	else{
		f_ref = 0;
		phic = params->phiRef;
	}
	
	//Assign shift: first shift so coalescence happens at t=0, then shift from there according to tc
	//This aligns more with the physical meaning of tc, but the phase is NO LONGER just
	if(params->shift_time){
		T alpha1_offset = assign_lambda_param_element(params,14);
		tc_shift = this->Dphase_mr(params->f3, params, &lambda)+(-lambda.alpha[1]+alpha1_offset)*params->M/params->eta;
	}
	else{
		tc_shift=0;
	}
	tc = 2*M_PI*params->tc - tc_shift;
	//################################################################

	//T A0 = sqrt(M_PI/30)*chirpmass*chirpmass/DL * pow(M_PI*chirpmass,-7./6);
	T A0 = params->A0* pow(M,7./6.);

	std::complex<T> amp, phase;
	std::complex<T> i;
	i = std::complex<T> (0,1.);
	T fcut = .2/M; //Cutoff frequency for IMRPhenomD - all higher frequencies return 0
	if(frequency>fcut){
		return std::complex<T>(0.0,0.0);
	}
	else{	
	if (frequency<params->f1_phase)
	{
		this->precalc_powers_ins(frequency, M, &pows);
	}
	else
	{
		pows.MFsixth = pow(M*frequency,1./6.);	
		pows.MF7sixth= pows.MFsixth*pows.MFsixth*pows.MFsixth*pows.MFsixth*pows.MFsixth*pows.MFsixth*pows.MFsixth;
	}
	amp = (A0 * this->build_amp(frequency,&lambda,params,&pows,pn_amp_coeffs,deltas));
	phase = (this->build_phase(frequency,&lambda,params,&pows,pn_phase_coeffs));
	//CORRECTION -- 05-21-2023!! Switched to (+)
	//phase -=   (T)(tc*(frequency-f_ref) - phic);
	phase -=   (T)(tc*(frequency-f_ref) + phic);
	return amp * std::exp(-i * phase);
	}

}
/*! \brief Constructs the Amplitude as outlined by IMRPhenomD
 *
 * arguments:
 * 	array of frequencies, length of that array, T array for the output amplitude, and a source_parameters structure
 */
template <class T>
int IMRPhenomD<T>::construct_amplitude(T *frequencies, /**< T array of frequencies the waveform is to be evaulated at*/
				int length,/**< integer length of the input array of frequencies and the output array*/ 
				T *amplitude,/**< output T array for the amplitude*/ 
				source_parameters<T> *params/**< Structure of source parameters to be initilized before computation*/
				)
{
	T M = params-> M;
	T chirpmass = params->chirpmass;
	T DL = params->DL;
	lambda_parameters<T> lambda;
	this->assign_lambda_param(params, &lambda);

	/*Initialize the post merger quantities*/
	this->post_merger_variables(params);

	params->f1 = 0.014/(params->M);
	params->f3 = this->fpeak(params, &lambda);

	useful_powers<T> pows;
	this->precalc_powers_PI(&pows);

	T deltas[6];
	T pn_amp_coeffs[7];

	this->assign_pn_amplitude_coeff(params, pn_amp_coeffs);

	this->amp_connection_coeffs(params,&lambda,pn_amp_coeffs,deltas);

	//T A0 = sqrt(M_PI/30)*chirpmass*chirpmass/DL * pow(M_PI*chirpmass,-7./6);
	T A0 = params->A0 * pow(M,7./6.);

	T f;
	T fcut = .2/M; //Cutoff frequency for IMRPhenomD - all higher frequencies return 0
	for (int j =0; j< length; j++)
	{
		f = frequencies[j];
		if(f>fcut){
			amplitude[j] = 0.0;
		}
		else{	
		if (f<params->f1)
		{
			this->precalc_powers_ins_amp(f, M, &pows);
		}
		else
		{
			pows.MFsixth = pow(M*f,1./6.);	
			pows.MF7sixth= pows.MFsixth*pows.MFsixth*pows.MFsixth*pows.MFsixth*pows.MFsixth*pows.MFsixth*pows.MFsixth;
		}
		amplitude[j] = (A0 * this->build_amp(f,&lambda,params,&pows,pn_amp_coeffs,deltas));

		}
	}
	return 1;
}
/*! \brief Constructs the Phase as outlined by IMRPhenomD
 *
 * arguments:
 * 	array of frequencies, length of that array, T array for the output phase, and a source_parameters structure
 */
template <class T>
int IMRPhenomD<T>::construct_phase(T *frequencies, /**< T array of frequencies the waveform is to be evaluated at*/
				int length,/**< integer length of the input and output arrays*/ 
				T *phase,/**< output T array for the phasee*/ 
				source_parameters<T> *params/**< structure of source parameters to be calculated before computation*/
				)
{
	T M = params-> M;
	lambda_parameters<T> lambda, *lambda_ptr;
	this->assign_lambda_param(params, &lambda);

	/*Initialize the post merger quantities*/
	this->post_merger_variables(params);

	params->f1_phase = 0.018/(params->M);
	params->f2_phase = params->fRD/2.;
	params->f3 = this->fpeak(params, &lambda);

	useful_powers<T> pows;
	this->precalc_powers_PI(&pows);

	T pn_phase_coeffs[12];

	this->assign_static_pn_phase_coeff(params, pn_phase_coeffs);	

	this->phase_connection_coefficients(params,&lambda,pn_phase_coeffs);
	
	//################################################################
	//Calculate phase and coalescence time variables
	T phic, f_ref, tc, phi_shift, tc_shift;
	//If phic is unspecified - use f_ref and phiRef
	if(params->shift_phase ){
		f_ref = params->f_ref;
		precalc_powers_ins(f_ref, M, &pows);
		phi_shift = (this->build_phase(f_ref,&lambda,params,&pows,pn_phase_coeffs));
		phic = 2*params->phiRef + phi_shift;
	}
	//If phic is specified, ignore f_ref phiRef and use phic
	else{
		f_ref = 0;
		phic = params->phiRef;
	}
	
	//Assign shift: first shift so coalescence happens at t=0, then shift from there according to tc
	//This aligns more with the physical meaning of tc, but the phase is NO LONGER just
	if(params->shift_time){
		T alpha1_offset = assign_lambda_param_element(params,14);
		tc_shift = this->Dphase_mr(params->f3, params, &lambda)+(-lambda.alpha[1]+alpha1_offset)*params->M/params->eta;
	}
	else{
		tc_shift=0;
	}
	//CORRECTION -- 05-21-2023!! Switched to (+)
	//tc = 2*M_PI*params->tc - tc_shift;
	tc = 2*M_PI*params->tc + tc_shift;
	//################################################################
	
	
	T f;
	
	for (int j =0; j< length; j++)
	{
		f = frequencies[j];
		if (f<params->f1_phase)
		{
			this->precalc_powers_ins_phase(f, M, &pows);
			
		}
		phase[j] =( this->build_phase(f,&lambda,params,&pows,pn_phase_coeffs));
		//CORRECTION -- 05-21-2023!! Switched to (+)
		//phase[j] -=   (T)(tc*(f-f_ref) - phic);
		phase[j] -=   (T)(tc*(f-f_ref) + phic);
		phase[j]*=(-1);

	}
	return 1;
}
/*! \brief constructs the IMRPhenomD amplitude for frequency f
 *
 * arguments: 
 * 	numerical parameters from Khan et al lambda_parameters structure, source_parameters structure, useful_powers<T> structure, PN parameters for the inspiral portions of the waveform, and the delta parameters for the intermediate region, numerically solved for using the amp_connection_coeffs function
 */
template <class T>
T IMRPhenomD<T>::build_amp(T f, 
		lambda_parameters<T> *lambda, 
		source_parameters<T> *params, 
		useful_powers<T> *pows,
		T *amp_coeff, 
		T *deltas)
{
	//##################################################3
	//##################################################3
	//##################################################3
	//TESTING -- FOR THE LOVE OF GOD DON'T FORGET ABOUT THIS PLEASE
	//##################################################3
	//##################################################3
	//debugger_print(__FILE__,__LINE__,"DO NOT FORGET THIS");
	//T pn_amp=0;
	//T M = params->M;

	//pn_amp = amp_coeff[0];// + amp_coeff[1] * pows->PIthird * pows->MFthird + 0;
	//	 //amp_coeff[2] * pows->PI2third * pows->MF2third + 
	//	 //amp_coeff[3] * M_PI * M * f + 
	//	 //amp_coeff[4] * pows->PI4third * pows->MF4third +
	//	 //amp_coeff[5] * pows->PI5third * pows->MF5third +
	//	 //amp_coeff[6] * pows->PIsquare * pows->MFsquare ;
	//
	//T nr_amp=0;
	////nr_amp= (lambda->rho[0]) * pows->MF7third +0;
	////	(lambda->rho[1]) * pows->MF8third +
	////	(lambda->rho[2]) * pows->MFcube ;
	//
	//return 1./(pows->MF7sixth) *pn_amp+nr_amp;
	//##################################################3
	//##################################################3
	//##################################################3
	//##################################################3
	//##################################################3
	if(f<params->f1)	
	{
		//return  pow(f,-7./6)* this->amp_ins(f, params, amp_coeff, lambda, pows);
		return  1./(pows->MF7sixth) * this->amp_ins(f, params, amp_coeff, lambda, pows);
	}
	else if(f>params->f3)	
	{
		//return pow(f,-7./6)* this->amp_int(f, params,  lambda, deltas);
		return 1./pows->MF7sixth* this->amp_mr(f, params, lambda);
	}
	else 	{
		//return pow(f,-7./6)* this->amp_mr(f, params, lambda);
		return 1./pows->MF7sixth* this->amp_int(f, params,  lambda, deltas);
	}
}

/*! \brief constructs the IMRPhenomD phase for frequency f
 *
 * arguments: 
 * 	numerical parameters from Khan et al lambda_parameters structure, source_parameters structure, useful_powers structure, PN parameters for the inspiral portions of the waveform
 */
template <class T>
T IMRPhenomD<T>::build_phase(T f, 
		lambda_parameters<T> *lambda, 
		source_parameters<T> *params, 
		useful_powers<T> *pows,
		T *phase_coeff)
{
	if(f < params->f1_phase)
	{
		return this->phase_ins(f, params, phase_coeff,lambda, pows);
	}
	else if(f > params->f2_phase)
	{
		return this->phase_mr(f, params, lambda);
	}
	else {
		return this->phase_int(f,params,lambda);
	}
}





//#####################################################################################
//

/*!\brief Wrapper for the Lambda parameter assignment that handles the looping
 */
template <class T>
void IMRPhenomD<T>::assign_lambda_param(source_parameters<T> *source_param, lambda_parameters<T> *lambda)
{
	for (int i=0;i<3; i++)
		lambda->rho[i] = this->assign_lambda_param_element(source_param, i);
	lambda->v2 = assign_lambda_param_element(source_param, 3);
	for (int i=0;i<3; i++)
		lambda->gamma[i] = this->assign_lambda_param_element(source_param, i+4);
	for (int i=0;i<4; i++)
		lambda->sigma[i+1] = this->assign_lambda_param_element(source_param, i+7);
	for (int i=0;i<3; i++)
		lambda->beta[i+1] = this->assign_lambda_param_element(source_param, i+11);
	for (int i=0;i<5; i++)
		lambda->alpha[i+1] = this->assign_lambda_param_element(source_param, i+14);


	//DO NOT LEAVE IN
	//WAVEFORM SYSTEMATICS TESTING
	//lambda->sigma[4] = 0;
}

/*!\brief Calculate the lambda parameters from Khan et al for element i
 */
template <class T>
T IMRPhenomD<T>::assign_lambda_param_element(source_parameters<T> *source_param,int i)
{
	T param_list[11];
	for (int j = 0; j<11;j++)
		param_list[j] = lambda_num_params[i][j];
	
	T spin_coeff = source_param->chi_pn - 1;
	T eta = source_param-> eta;
        T parameter = param_list[0] + param_list[1]*eta + 
            (spin_coeff)*(param_list[2] + param_list[3]*eta + param_list[4]*eta*eta) + 
            pow(spin_coeff,2)*(param_list[5] + param_list[6]*eta+param_list[7]*eta*eta) +
            pow(spin_coeff,3)*(param_list[8] + param_list[9]*eta+param_list[10]*eta*eta);

        return parameter;	
}

/*!\brief Pre-calculate powers of Mf, to speed up calculations for the inspiral waveform (both amplitude and phase
 *
 * It seems the pow() function is very slow, so to speed things up, powers of Mf will be precomputed and passed to the functions within the frequency loops
 */
template <class T>
void IMRPhenomD<T>::precalc_powers_ins(T f, T M, useful_powers<T> *Mf_pows)
{
	T Mf = M*f;
	Mf_pows->MFsquare = Mf*Mf;
	Mf_pows->MFcube = Mf*Mf*Mf;
	Mf_pows->MFsixth = pow(M*f,1./6.);
	Mf_pows->MF7sixth= Mf_pows->MFsixth * M*f;
	Mf_pows->MFthird = Mf_pows->MFsixth*Mf_pows->MFsixth;
	Mf_pows->MF2third = Mf_pows->MFthird * Mf_pows->MFthird;
	Mf_pows->MF4third = Mf_pows->MFthird * Mf;
	Mf_pows->MF5third = Mf_pows->MF2third * Mf;
	Mf_pows->MF7third = Mf_pows->MF4third * Mf;
	Mf_pows->MF8third = Mf_pows->MF5third * Mf;
	Mf_pows->MFminus_5third = 1./Mf_pows->MF5third;
	
}
/*!\brief Pre-calculate powers of Mf, to speed up calculations for the inspiral amplitude
 *
 * It seems the pow() function is very slow, so to speed things up, powers of Mf will be precomputed and passed to the functions within the frequency loops
 */
template <class T>
void IMRPhenomD<T>::precalc_powers_ins_amp(T f, T M, useful_powers<T> *Mf_pows)
{
	T Mf = M*f;
	Mf_pows->MFsquare = Mf*Mf;
	Mf_pows->MFcube = Mf*Mf*Mf;
	Mf_pows->MFsixth = pow(M*f,1./6.);
	Mf_pows->MF7sixth= Mf_pows->MFsixth * M*f;
	Mf_pows->MFthird = Mf_pows->MFsixth*Mf_pows->MFsixth;
	Mf_pows->MF2third = Mf_pows->MFthird * Mf_pows->MFthird;
	Mf_pows->MF4third = Mf_pows->MFthird * Mf;
	Mf_pows->MF5third = Mf_pows->MF2third * Mf;
	Mf_pows->MF7third = Mf_pows->MF4third * Mf;
	Mf_pows->MF8third = Mf_pows->MF5third * Mf;
	
}

/*!\brief Pre-calculate powers of Mf, to speed up calculations for the inspiral phase
 *
 * It seems the pow() function is very slow, so to speed things up, powers of Mf will be precomputed and passed to the functions within the frequency loops
 */
template <class T>
void IMRPhenomD<T>::precalc_powers_ins_phase(T f, T M, useful_powers<T> *Mf_pows)
{
	T Mf = M*f;
	Mf_pows->MFsquare = Mf*Mf;
	Mf_pows->MFcube = Mf*Mf*Mf;
	Mf_pows->MFsixth = pow(M*f,1./6.);
	Mf_pows->MF7sixth= Mf_pows->MFsixth * M*f;
	Mf_pows->MFthird = Mf_pows->MFsixth*Mf_pows->MFsixth;
	Mf_pows->MF2third = Mf_pows->MFthird * Mf_pows->MFthird;
	Mf_pows->MF4third = Mf_pows->MFthird * Mf;
	Mf_pows->MF5third = Mf_pows->MF2third * Mf;
	Mf_pows->MF7third = Mf_pows->MF4third * Mf;
	Mf_pows->MFminus_5third = 1./Mf_pows->MF5third;
	
}
/*!\brief Pre-calculate powers of pi, to speed up calculations for the inspiral phase
 *
 * It seems the pow() function is very slow, so to speed things up, powers of PI will be precomputed and passed to the functions within the frequency loops
 */
template <class T>
void IMRPhenomD<T>::precalc_powers_PI( useful_powers<T> *PI_pows)
{
	PI_pows->PIsquare = M_PI*M_PI;
	PI_pows->PIcube = M_PI*M_PI*M_PI;
	PI_pows->PIthird = pow(M_PI,1./3.);
	PI_pows->PI2third = PI_pows->PIthird * PI_pows->PIthird;
	PI_pows->PI4third = PI_pows->PIthird * M_PI;
	PI_pows->PI5third = PI_pows->PI2third * M_PI;
	PI_pows->PI7third = PI_pows->PI4third * M_PI;
	PI_pows->PIminus_5third = 1./PI_pows->PI5third;
	
}
/*!
 * \brief Calculates the static PN coeffecients for the amplitude
 */
template <class T>
void IMRPhenomD<T>::assign_pn_amplitude_coeff(source_parameters<T> *source_param, 
				T*coeff)
{
	T massdelta = source_param->delta_mass;
	T eta = source_param->eta;
	T chi_a = source_param->chi_a;
	T chi_s = source_param->chi_s;
	T eta2 = eta*eta;
	T eta3 = eta2*eta;
	T chi_s2 = chi_s*chi_s;
	T chi_a2 = chi_a*chi_a;
	double pi2 = M_PI*M_PI;
	
 	coeff[0] = 1.;
    	coeff[1] = 0.;
    	coeff[2] = (-323./224 + 451*eta/168);
    	coeff[3] = (27.*massdelta*chi_a/8 + (27./8 - 11.*eta/6)*chi_s);
    	coeff[4] = (-27312085./8128512 - 1975055*eta/338688 + 
    		105271*eta2/24192 + (-81./32 + 8*eta)*chi_a2 - 
    		81*massdelta*chi_a*chi_s/16 + (-81./32 + 17*eta/8)*chi_s2);

    	coeff[5] = 1.*(-85*M_PI/64 + 85*M_PI*eta/16 + massdelta*(285197./16128 - 
    		1579*eta/4032)*chi_a + (285197./16128 - 15317*eta/672 - 
    		2227*eta2/1008)*chi_s);

    	coeff[6] = 1.*(-177520268561./8583708672 +(545384828789./5007163392 - 205*pi2/48)*eta-
    		3248849057*eta2/178827264 + 34473079*eta3/6386688 + 
    		(1614569./64512 - 1873643.*eta/16128 + 2167*eta2/42)*chi_a2 + 
    		(31*M_PI/12 - 7*M_PI*eta/3)*chi_s + (1614569./64512 - 61391*eta/1344 + 
    		57451*eta2/4032)*chi_s2 + 
    		massdelta*chi_a*(31*M_PI/12 + (1614569./32256 - 165961*eta/2688)*chi_s));
}

/*! 
 * \brief Calculates the static PN coeffecients for the phase - coeffecients 0,1,2,3,4,7
 */
template <class T>
void IMRPhenomD<T>::assign_static_pn_phase_coeff(source_parameters<T> *source_param, 
					T *coeff)
{

	T delta = source_param->delta_mass;
	T eta = source_param->eta;
	T chi_a = source_param->chi_a;
	T chi_s = source_param->chi_s;
	T eta2 = eta*eta;
	T eta3 = eta2*eta;
	T chi_s2 = chi_s*chi_s;
	T chi_a2 = chi_a*chi_a;
	double pi2 = M_PI*M_PI;
	
    	coeff[0] = 1.;
    	coeff[1] =  0.;
    	coeff[2] =  3715./756 + 55.*eta/9;
    	coeff[3] = -16.*M_PI + 113.*delta*chi_a/3 + 
    			(113./3 - 76.*eta/3)*chi_s;

    	coeff[4] = 15293365./508032 + 27145.*eta/504 + 3085.*eta2/72 + 
    			(-405./8 + 200.*eta)*chi_a2 - 
    			(405./4)*delta*chi_a*chi_s +
    			(-405./8 + 5.*eta/2)*chi_s2;

    	coeff[7] = 77096675.*M_PI/254016 + 378515.*M_PI*eta/1512 - 
    			74045.*M_PI*eta2/756 + delta*(-25150083775./3048192 + 
    			26804935.*eta/6048 - 1985.*eta2/48)*chi_a + 
    			(-25150083775./3048192 + 10566655595.*eta/762048 - 
    			1042165.*eta2/3024 + 5345.*eta3/36)*chi_s;
	coeff[8] = (38645.*M_PI/756 - 65.*M_PI*eta/9 + 
    		delta*(-732985./2268 - 140.*eta/9)*chi_a + 
    		(-732985./2268 + 24260.*eta/81 + 340.*eta2/9)*chi_s);
	coeff[9] = coeff[8];
	coeff[10] = 11583231236531./4694215680 - 6848.*gamma_E/21 -\
     		640.*pi2/3 + (-15737765635./3048192 + 2255.*pi2/12)*eta + 
     		76055.*eta2/1728 - 127825.*eta3/1296 
     		+ 2270.*delta*chi_a*M_PI/3 + 
     		(2270.*M_PI/3 - 520.*M_PI*eta)*chi_s- 6848.*( log_64)/63;
	coeff[11] = 6848./63.;
}

/*! 
 * \brief Calculates the dynamic PN phase coefficients 5,6
 *
 * f is in Hz
 */
template <class T>
void IMRPhenomD<T>::assign_nonstatic_pn_phase_coeff(source_parameters<T> *source_param, 
					T *coeff, 
					T f)
{
	T delta = source_param->delta_mass;
	T eta = source_param->eta;
	T chi_a = source_param->chi_a;
	T chi_s = source_param->chi_s;
	T eta2 = eta*eta;
	T eta3 = eta2*eta;
	T chi_s2 = chi_s*chi_s;
	T chi_a2 = chi_a*chi_a;
	double pi2 = M_PI*M_PI;
	T M = source_param->mass1 + source_param->mass2;
	T logF = log(M_PI*M*f);

	//coeff[5] = coeff[8] *(1.+logF);
	//coeff[6] = coeff[9]-6848.*(logF)/63.;
	coeff[5] = coeff[8] + logF*coeff[9];
	coeff[6] = coeff[10]-coeff[11]*(logF);
}

/*! 
 * \brief Calculates the derivative of the dynamic PN phase coefficients 5,6
 *
 * f is in Hz
 */
template <class T>
void IMRPhenomD<T>::assign_nonstatic_pn_phase_coeff_deriv(source_parameters<T> *source_param, 
					T Dcoeff[], 
					T f)
{
	T delta = source_param->delta_mass;
	T eta = source_param->eta;
	T chi_a = source_param->chi_a;
	T chi_s = source_param->chi_s;
	T eta2 = eta*eta;

	Dcoeff[0] = (1./f)* (38645.*M_PI/756 - 65*M_PI*eta/9 + 
    		delta*(-732985./2268 - 140*eta/9)*chi_a + 
    		(-732985./2268 + 24260*eta/81 + 340*eta2/9)*chi_s);
	
    	Dcoeff[1] = (-6848.*(1./f)/63);
	
}

/*! \brief Calculates the post-merger ringdown frequency and dampening frequency
 *
 * Returns in Hz - assigns fRD to var[0] and fdamp to var[1]
 */
template <class T>
void IMRPhenomD<T>::post_merger_variables(source_parameters<T> *source_param)
{
	if(!source_param->dep_postmerger){
		calc_fring(source_param);
		calc_fdamp(source_param);
	}
	else{
		_calc_fring(source_param);
		_calc_fdamp(source_param);
	}
}

//#################################################################################
//#################################################################################
//ALL taken directly from LALsuite
//#################################################################################
//#################################################################################
//Matches ADOL-C
//#################################################################################
/*! \brief Deprecated version, but ADOL-C version -- can be forced to use if dep_postmerger flag set to TRUE
 */
template <class T>
void IMRPhenomD<T>::_calc_fring( source_parameters<T> *source_param)
{
	T chi1 = source_param->spin1z;
	T chi2 = source_param->spin2z;
	T eta = source_param->eta;
	T M = source_param->M;
	T m1 = source_param->mass1;
	T m2 = source_param->mass2;
	T eta2 = eta*eta;
	T eta3 = eta2*eta;
	T eta4 = eta3*eta;
	T m12 = m1*m1;
	T m22 = m2*m2;
	T M2 = M*M;
 	T S = (chi1*m12 + chi2*m22)/M2 ;
	T S2 = S*S;
	T S3 = S2*S;
	T S4 = S3*S;

    	T S_red = S/(1.-2.*eta);

    	T a = final_spin(source_param);

    	T E_rad_ns = 0.0559745*eta +0.580951*eta2 - 
    	 0.960673*eta3 + 3.35241*eta4 ;

	T E_rad = E_rad_ns*(1.+S_red*(-0.00303023 - 2.00661*eta +7.70506*eta2)) / 
		(1+ S_red*(-0.67144 - 1.47569*eta +7.30468*eta2));
	T MWRD = (1.5251-1.1568*pow(1-a,0.1292));
	T MWdamp = ((1.5251-1.1568*pow(1.-a,0.1292))/(2.*(0.700 + 1.4187*pow(1.-a,-.4990))));
	source_param->fRD =  (1./(2*M_PI))*(MWRD)/(M*(1. - E_rad));
}
/*! \brief Deprecated version, but ADOL-C version -- can be forced to use if dep_postmerger flag set to TRUE
 */
template <class T>
void IMRPhenomD<T>::_calc_fdamp(source_parameters<T> *source_param)
{
	T chi1 = source_param->spin1z;
	T chi2 = source_param->spin2z;
	T eta = source_param->eta;
	T M = source_param->M;
	T m1 = source_param->mass1;
	T m2 = source_param->mass2;
	T eta2 = eta*eta;
	T eta3 = eta2*eta;
	T eta4 = eta3*eta;
	T m12 = m1*m1;
	T m22 = m2*m2;
	T M2 = M*M;
 	T S = (chi1*m12 + chi2*m22)/M2 ;
	T S2 = S*S;
	T S3 = S2*S;
	T S4 = S3*S;

    	T S_red = S/(1.-2.*eta);
    	T a = final_spin(source_param);

    	T E_rad_ns = 0.0559745*eta +0.580951*eta2 - 
    	 0.960673*eta3 + 3.35241*eta4 ;

	T E_rad = E_rad_ns*(1.+S_red*(-0.00303023 - 2.00661*eta +7.70506*eta2)) / 
	  (1+ S_red*(-0.67144 - 1.47569*eta +7.30468*eta2));
	T MWRD = (1.5251-1.1568*pow(1-a,0.1292));
	T MWdamp = ((1.5251-1.1568*pow(1.-a,0.1292))/(2.*(0.700 + 1.4187*pow(1.-a,-.4990))));
	source_param->fdamp = (1./(2*M_PI))*(MWdamp)/(M*(1. - E_rad));

}
//#################################################################################
//Matches LALSuite
template <>
void IMRPhenomD<double>::calc_fring( source_parameters<double> *source_param)
{
	double eta = source_param->eta;
	// Check that eta is physical.
	// Numerical round-off sometimes causes it to not be,
	// and causes nan in Erad
	// (also in final_spin, but that is checked inside of it).
	if (eta > 0.25)
	{
		eta = 0.25;
	}
	double chi1 = source_param->spin1z;
	double chi2 = source_param->spin2z;
	double finspin = final_spin(source_param);
	
	double Erad = EradRational0815(eta,chi1,chi2);
	gsl_interp_accel *acc = gsl_interp_accel_alloc();
  	gsl_spline *iFring = gsl_spline_alloc(gsl_interp_cspline, QNMData_length);
  	gsl_spline_init(iFring, QNMData_a, QNMData_fring, QNMData_length);

  	double return_val = gsl_spline_eval(iFring, finspin, acc) / (1.0 - Erad);

  	gsl_spline_free(iFring);
  	gsl_interp_accel_free(acc);
	
  	source_param->fRD = return_val/source_param->M;
}
template <>
void IMRPhenomD<double>::calc_fdamp(source_parameters<double> *source_param)
{
	double eta = source_param->eta;
	double chi1 = source_param->spin1z;
	double chi2 = source_param->spin2z;
    	double finspin = final_spin(source_param);
	double Erad = EradRational0815(eta,chi1,chi2);
	gsl_interp_accel *acc = gsl_interp_accel_alloc();
	gsl_spline *iFdamp = gsl_spline_alloc(gsl_interp_cspline, QNMData_length);
	gsl_spline_init(iFdamp, QNMData_a, QNMData_fdamp, QNMData_length);
	
	double return_val = gsl_spline_eval(iFdamp, finspin, acc) / (1.0 - Erad);
	
	gsl_spline_free(iFdamp);
	gsl_interp_accel_free(acc);
	source_param->fdamp =  return_val/source_param->M;

}
//#################################################################################
template <>
void IMRPhenomD<adouble>::calc_fring( source_parameters<adouble> *source_param)
{
	_calc_fring(source_param);
}
template <>
void IMRPhenomD<adouble>::calc_fdamp(source_parameters<adouble> *source_param)
{
	_calc_fdamp(source_param);

}

template<class T>
T IMRPhenomD<T>::final_spin(source_parameters<T> *params)
{
	return FinalSpin0815(params->eta, params->spin1z,params->spin2z);
}

/**
 * Formula to predict the final spin. Equation 3.6 arXiv:1508.07250
 * s defined around Equation 3.6.
 */
template<class T>
T IMRPhenomD<T>::FinalSpin0815_s(T eta, T s) {
	T eta2 = eta*eta;
	T eta3 = eta2*eta;
	T s2 = s*s;
	T s3 = s2*s;

	return eta*(3.4641016151377544 - 4.399247300629289*eta +
	      9.397292189321194*eta2 - 13.180949901606242*eta3 +
	      s*((1.0/eta - 0.0850917821418767 - 5.837029316602263*eta) +
	      (0.1014665242971878 - 2.0967746996832157*eta)*s +
	      (-1.3546806617824356 + 4.108962025369336*eta)*s2 +
	      (-0.8676969352555539 + 2.064046835273906*eta)*s3));
}
/**
 * Wrapper function for FinalSpin0815_s.
 */
template< class T>
T IMRPhenomD<T>::FinalSpin0815(T eta, T chi1, T chi2) {
  // Convention m1 >= m2
  // Check that eta is physical.
  // Numerical issues sometimes cause it to not be,
  // and Seta is returned as nan.
  T eta_max = 1/4;
  T Seta;
  if (eta >= eta_max)
  {
	Seta = 0;
  }
  else
  {
	Seta = sqrt(1.0 - 4.0*eta);
  }
  T m1 = 0.5 * (1.0 + Seta);
  T m2 = 0.5 * (1.0 - Seta);
  T m1s = m1*m1;
  T m2s = m2*m2;
  // s defined around Equation 3.6 arXiv:1508.07250
  T s = (m1s * chi1 + m2s * chi2);
  return FinalSpin0815_s(eta, s);
}
/**
 * Formula to predict the total radiated energy. Equation 3.7 and 3.8 arXiv:1508.07250
 * Input parameter s defined around Equation 3.7 and 3.8.
 */
template<class T>
T IMRPhenomD<T>::EradRational0815_s(T eta, T s) {
  T eta2 = eta*eta;
  T eta3 = eta2*eta;

  return (eta*(0.055974469826360077 + 0.5809510763115132*eta - 0.9606726679372312*eta2 + 3.352411249771192*eta3)*
    (1. + (-0.0030302335878845507 - 2.0066110851351073*eta + 7.7050567802399215*eta2)*s))/(1. + (-0.6714403054720589 - 1.4756929437702908*eta + 7.304676214885011*eta2)*s);
}
/**
 * Wrapper function for EradRational0815_s.
 */
template<class T>
T IMRPhenomD<T>::EradRational0815(T eta, T chi1, T chi2) {
  // Convention m1 >= m2
  T Seta = sqrt(1.0 - 4.0*eta);
  T m1 = 0.5 * (1.0 + Seta);
  T m2 = 0.5 * (1.0 - Seta);
  T m1s = m1*m1;
  T m2s = m2*m2;
  // arXiv:1508.07250
  T s = (m1s * chi1 + m2s * chi2) / (m1s + m2s);

  return EradRational0815_s(eta, s);
}

//#################################################################################
//#################################################################################




/*!\brief Solves for the peak frequency, where the waveform transitions from intermediate to merger-ringdown
 *
 * returns Hz
 */
template <class T>
T IMRPhenomD<T>::fpeak(source_parameters<T> *params, lambda_parameters<T> *lambda)
{
	T fRD = params->fRD;
	T fdamp = params->fdamp;
	T gamma3 = lambda->gamma[2];
	T gamma2 = lambda->gamma[1];
	T out;
	//return abs(fRD + fdamp*gamma3*(sqrt(1-gamma2*gamma2)-1)/gamma2);
	if(gamma2>1){
		out = fRD + (fdamp*(-1.)*gamma3)/gamma2;
		return sqrt(out*out);
	}
	else
		return sqrt( (fRD + fdamp*gamma3*(sqrt(1-gamma2*gamma2)-1)/gamma2)*(fRD + fdamp*gamma3*(sqrt(1-gamma2*gamma2)-1)/gamma2));
}


/*! \brief Calculates the scaled inspiral amplitude A/A0 for frequency f with precomputed powers of MF and PI
 *
 * return a T
 *
 * additional argument contains useful powers of MF and PI in structure userful_powers
 */
template <class T>
T IMRPhenomD<T>::amp_ins(T f, source_parameters<T> *param, T *pn_coeff, 
			lambda_parameters<T> *lambda, useful_powers<T> *pow)
{
	T pn_amp;
	T M = param->M;

	pn_amp = pn_coeff[0] + pn_coeff[1] * pow->PIthird * pow->MFthird +
		 pn_coeff[2] * pow->PI2third * pow->MF2third + 
		 pn_coeff[3] * M_PI * M * f + 
		 pn_coeff[4] * pow->PI4third * pow->MF4third +
		 pn_coeff[5] * pow->PI5third * pow->MF5third +
		 pn_coeff[6] * pow->PIsquare * pow->MFsquare ;
	
	T nr_amp;
	nr_amp= (lambda->rho[0]) * pow->MF7third +
		(lambda->rho[1]) * pow->MF8third +
		(lambda->rho[2]) * pow->MFcube ;
	
	return pn_amp+nr_amp;
}	


/*! \brief Calculates the inspiral phase for frequency f with precomputed powers of MF and PI for speed
 *
 * return a T
 *
 * extra argument of precomputed powers of MF and pi, contained in the structure useful_powers<T>
 */
template <class T>
T IMRPhenomD<T>::phase_ins(T f, source_parameters<T> *param, T *pn_coeff, 
		lambda_parameters<T> *lambda, useful_powers<T> *pow)
{
	T M = param->M;	
	T eta = param->eta;
	this->assign_nonstatic_pn_phase_coeff(param, pn_coeff, f);
	T pn_phase;
	pn_phase = pn_coeff[0] + pn_coeff[1] * pow->PIthird * pow->MFthird +
		 pn_coeff[2] * pow->PI2third * pow->MF2third + 
		 pn_coeff[3] * M_PI * M * f + 
		 pn_coeff[4] * pow->PI4third * pow->MF4third +
		 pn_coeff[5] * pow->PI5third * pow->MF5third +
		 pn_coeff[6] * pow->PIsquare * pow->MFsquare +
		 pn_coeff[7] * pow->PI7third * pow->MF7third ;
		
	if (param->PNorder != 35){
		/* Add higher order terms if needed*/
	}

	T phase_TF2 =  -M_PI/4. 
		+ 3./(128.*eta) * pow->PIminus_5third * pow->MFminus_5third * pn_phase;

	/*sigma0 and sigma1 can be reabsorbed into tc and phic*/	
	T sigma0 = 0;
        //T sigma1 =0;
        T sigma1 = lambda->sigma[1];
	T sigma2 = lambda->sigma[2];
	T sigma3 = lambda->sigma[3];
	T sigma4 = lambda->sigma[4];

	//std::cout<<pn_coeff[0]<<" "<<pn_coeff[1]<<" "<<pn_coeff[2]<<" "<<pn_coeff[3]<<" "<<pn_coeff[4]<<" "<<pn_coeff[5]<<" "<<pn_coeff[6]<<" "<<pn_coeff[7]<<" "<<lambda->sigma[1]<<" "<<lambda->sigma[2]<<" "<<lambda->sigma[3]<<" "<<lambda->sigma[4]<<" "<<std::endl;
        return phase_TF2 + (1./eta)*(sigma0 + sigma1*M*f + 
        (3./4.)*sigma2* pow->MF4third + (3./5)*sigma3* pow->MF5third + 
        (1./2.)*sigma4* pow->MFsquare);
}
/*! \brief Calculates the derivative wrt frequency for the scaled inspiral amplitude A/A0 for frequency f
 *
 * This is an analytic derivative for the smoothness condition on the amplitude connection 
 *
 * return a T
 */
template <class T>
T IMRPhenomD<T>::Damp_ins(T f, source_parameters<T> *param, T *pn_coeff, lambda_parameters<T> *lambda)
{
	T pn_amp=0;
	T M = param->M;
	for(int i =0; i<7;i++)
		pn_amp+= pn_coeff[i]*pow(M_PI*(M),i/3.)*pow(f,i/3.-1.)*(i/3.);
	
	T nr_amp=0;
	for(int i=0;i<3;i++)
		nr_amp+= (lambda->rho[i])*pow(M,(7.+i)/3.)*pow(f,(7.+i)/3.-1.)*(7.+i)/3.;
	
	return pn_amp+nr_amp;
}	

/*! \brief Calculates the derivative of the inspiral phase for frequency f
 *
 * For phase continuity and smoothness
 * return a T
 */
template <class T>
T IMRPhenomD<T>::Dphase_ins(T f, source_parameters<T> *param, T *pn_coeff, lambda_parameters<T> *lambda)
{
	T M = param->M;	
	T eta = param->eta;
	
	this->assign_nonstatic_pn_phase_coeff(param, pn_coeff, f);

	T *pnderiv = new T[3];
	this->assign_nonstatic_pn_phase_coeff_deriv(param,pnderiv,f);
	
	T pn_phase=0;
	for (int i =0; i< 8; i++)
		pn_phase += pn_coeff[i]*pow(M_PI * M * f,i/3.);
	T phase_TF2 = 	+ 3./(128*eta) *pow(M_PI * M,-5./3.)*(-5./3.)*pow(f,-5./3-1.)*pn_phase;

	T pn_phase_deriv=0;
	for (int i =0;i<8;i++)
		pn_phase_deriv += pn_coeff[i]*pow(M_PI * M ,i/3.)*(i/3.)*pow(f,i/3.- 1.);
	pn_phase_deriv+= pnderiv[0]*pow(M_PI*M*f,5./3) + pnderiv[1]*pow(M_PI*M*f,6./3);
	phase_TF2 += pn_phase_deriv * (3./(128*eta) *pow(M_PI*M*f,-5./3));

	/*sigma0 and sigma1 can be reabsorbed into tc and phic*/	
	T sigma0 = 0;
        T sigma1 =lambda->sigma[1];
        //T sigma1 =0 ;
	T sigma2 = lambda->sigma[2];
	T sigma3 = lambda->sigma[3];
	T sigma4 = lambda->sigma[4];
	
	delete [] pnderiv;
        return phase_TF2 + (1./eta)*(sigma1*M + 
        sigma2*pow(M,4./3.)*pow(f,1./3) + sigma3*pow(M,5./3.)*pow(f,2./3) + 
        sigma4*pow(M,2.)*f);
}


/*! \brief Calculates the scaled merger-ringdown amplitude A/A0 for frequency f
 *
 * return a T
 */
template <class T>
T IMRPhenomD<T>::amp_mr(T f, source_parameters<T> *param, lambda_parameters<T> *lambda)
{
	T M = param->M;
	T M2 = M*M;
	T numerator = ((lambda->gamma[0])*(lambda->gamma[2])*(param->fdamp)*(param->M)) *
				exp( -(lambda->gamma[1]) * (f - (param->fRD)) / 
				( (lambda->gamma[2])*(param->fdamp) ) );
	T denominator = M2*(( (f - (param->fRD))* (f - (param->fRD)))  
			+((lambda->gamma[2])*(param->fdamp)*(lambda->gamma[2])*(param->fdamp)) );
	return numerator/ denominator;
}

/*! \brief Calculates the merger-ringdown phase for frequency f
 *
 * return a T
 */
template <class T>
T IMRPhenomD<T>::phase_mr(T f, source_parameters<T> *param, lambda_parameters<T> *lambda)
{
	T M = param->M;
	T eta = param->eta;
	T fdamp = param->fdamp;
	T fRD = param->fRD;
	T alpha0 = lambda->alpha[0];
	T alpha1 = lambda->alpha[1];
	T alpha2 = lambda->alpha[2];
	T alpha3 = lambda->alpha[3];
	T alpha4 = lambda->alpha[4];
	T alpha5 = lambda->alpha[5];
        T Mf = M*f;
	T Mfcube = Mf*Mf*Mf;
	T Mf3fourths = sqrt(sqrt(Mfcube));
        return (1./eta)*(alpha0 +alpha1*(Mf) -alpha2*(1./(Mf)) + (4./3)*alpha3*Mf3fourths + 
        alpha4*atan((f-alpha5*fRD)/fdamp));

}

/*! \brief Calculates the derivative wrt frequency for the scaled merger-ringdown amplitude A/A0 for frequency f
 *
 * This is an analytic derivative for the smoothness condition on the amplitude connection 
 * 
 * The analytic expression was obtained from Mathematica - See the mathematica folder for code
 *
 * return a T
 */
template <class T>
T IMRPhenomD<T>::Damp_mr(T f, source_parameters<T> *param, lambda_parameters<T> *lambda)
{
	T gamma1 = lambda->gamma[0];
	T gamma2 = lambda->gamma[1];
	T gamma3 = lambda->gamma[2];
	T fdamp = param->fdamp;
	T fRD = param->fRD;
	T M = param->M;
	return -((exp(((-f + fRD)*gamma2)/
	         (fdamp*gamma3))*gamma1*
	       (pow_int(f - fRD,2)*gamma2 + 
	         2*fdamp*(f - fRD)*gamma3 + 
	         pow_int(fdamp,2)*gamma2*
	          pow(gamma3,2)))/
	     (pow(pow_int(f - fRD,2) + 
	         pow_int(fdamp,2)*pow_int(gamma3,2)
	         ,2)*M));
}

/*! \brief Calculates the derivative of the merger-ringdown phase for frequency f
 *
 * For phase continuity and smoothness
 * return a T
 */
template <class T>
T IMRPhenomD<T>::Dphase_mr(T f, source_parameters<T> *param, lambda_parameters<T> *lambda)
{
	T M = param->M;
	T eta = param->eta;
	T fdamp = param->fdamp;
	T fRD = param->fRD;
	T alpha0 = lambda->alpha[0];
	T alpha1 = lambda->alpha[1];
	T alpha2 = lambda->alpha[2];
	T alpha3 = lambda->alpha[3];
	T alpha4 = lambda->alpha[4];
	T alpha5 = lambda->alpha[5];
 	return (alpha4/(fdamp*(1. + pow(f - alpha5*fRD,2.)/pow(fdamp,2.))) + alpha2/(pow(f,2.)*M) +
         M*(alpha1 + alpha3/pow(f*M,0.25)))/eta;
}

/*! \brief Calculates the scaled intermediate range amplitude A/A0 for frequency f
 *
 * return a T
 */
template <class T>
T IMRPhenomD<T>::amp_int(T f, source_parameters<T> *param, lambda_parameters<T> *lambda, T *deltas)
{
	T M = param->M;
	T Mf = M*f;
	T Mf2 = Mf*Mf;
	T Mf3 = Mf2*Mf;
	T Mf4 = Mf3*Mf;
	return (deltas[0] +  deltas[1]*Mf +deltas[2]*Mf2 + 
        deltas[3]*Mf3 + deltas[4]*Mf4);
	
}

/*! \brief Calculates the intermediate phase for frequency f
 *
 * return a T
 */
template <class T>
T IMRPhenomD<T>::phase_int(T f, source_parameters<T> *param, lambda_parameters<T> *lambda)
{
	T beta0 = lambda->beta[0];
	T beta1 = lambda->beta[1];
	T beta2 = lambda->beta[2];
	T beta3 = lambda->beta[3];
	T M = param->M;
	T Mf = M*f;
	T Mf3 = Mf*Mf*Mf;
	T eta = param->eta;

	return (1./eta)*(beta0+ beta1*(Mf) + beta2*log(Mf) - beta3/3. *1./ Mf3);
}

/*! \brief Calculates the derivative of the intermediate phase for frequency f
 *
 * For phase continuity and smoothness
 * return a T
 */
template <class T>
T IMRPhenomD<T>::Dphase_int(T f, source_parameters<T> *param, lambda_parameters<T> *lambda)
{
	T beta0 = lambda->beta[0];
	T beta1 = lambda->beta[1];
	T beta2 = lambda->beta[2];
	T beta3 = lambda->beta[3];
	T M = param->M;
	T eta = param->eta;

	return (1./eta)*( beta1*(M) + beta2/f + beta3 *pow(M,-3.)*pow(f,-4.));
}

/*!\brief Calculates the phase connection coefficients alpha{0,1} and beta{0,1}
 *
 * Note: these coefficients are stored in the lambda parameter structure, not a separate array
 */
template <class T>
void IMRPhenomD<T>::phase_connection_coefficients(source_parameters<T> *param, 
				lambda_parameters<T> *lambda, 
				T *pn_coeffs)
{
	lambda->beta[0] = 0;
	lambda->beta[1] = 0;
	lambda->alpha[0] = 0;
	lambda->alpha[1] = 0;

	lambda->beta[1] =this->calculate_beta1(param,lambda,pn_coeffs);
	lambda->beta[0] =this->calculate_beta0(param,lambda,pn_coeffs);
	lambda->alpha[1] =this->calculate_alpha1(param,lambda);
	lambda->alpha[0] = this->calculate_alpha0(param,lambda);
}

template <class T>
T IMRPhenomD<T>::calculate_beta1(source_parameters<T> *param, lambda_parameters<T> *lambda, T *pn_coeffs)
{
	T M = param->M;
	T eta = param->eta;
	T f1 = param->f1_phase;
	T Dins = this->Dphase_ins(f1, param, pn_coeffs, lambda);
	T Dint = this->Dphase_int(f1, param, lambda);
	return (eta/M) * Dins - (eta/M)*Dint; 
}

template <class T>
T IMRPhenomD<T>::calculate_beta0(source_parameters<T> *param, lambda_parameters<T> *lambda, T *pn_coeffs)
{
	T M = param->M;
	T eta = param->eta;
	T f1 = param->f1_phase;
	useful_powers<T> powers;
	this->precalc_powers_ins(f1, M, &powers);	
	this->precalc_powers_PI(&powers);
	T ins = this->phase_ins(f1, param, pn_coeffs, lambda,&powers);
	T intval = this->phase_int(f1, param, lambda);
	return (eta) * ins - eta*intval; 
}

template <class T>
T IMRPhenomD<T>::calculate_alpha1(source_parameters<T> *param, lambda_parameters<T> *lambda)
{
	T M = param->M;
	T eta = param->eta;
	T f2 = param->f2_phase;
	T Dint = this->Dphase_int(f2, param, lambda);
	T Dmr = this->Dphase_mr(f2, param, lambda);
	return (eta/M) * Dint - eta/M*Dmr; 
}

template <class T>
T IMRPhenomD<T>::calculate_alpha0(source_parameters<T> *param, lambda_parameters<T> *lambda)
{
	T M = param->M;
	T eta = param->eta;
	T f2 = param->f2_phase;
	T intval = this->phase_int(f2, param, lambda);
	T mr = this->phase_mr(f2, param, lambda);
	return (eta) * intval - eta*mr; 
}

/*! \brief Solves for the connection coefficients to ensure the transition from inspiral to merger ringdown is continuous and smooth
 */
template <class T>
void IMRPhenomD<T>::amp_connection_coeffs(source_parameters<T> *param, 
			lambda_parameters<T> *lambda, 
			T *pn_coeffs, 
			T *coeffs)
{
	T M = param->M;
	T f1 = param->f1;
	T f3 = param->f3;
	useful_powers<T> powers;
	this->precalc_powers_ins(f1, M, &powers);	
	this->precalc_powers_PI(&powers);
	T f2 = (f1+f3)/2.;
	T v1 = this->amp_ins(f1,param,pn_coeffs,lambda, &powers);
	T v2 = lambda->v2;
	T v3 = this->amp_mr(f3,param,lambda);	
	T dd1 = this->Damp_ins(f1,param,pn_coeffs,lambda);
	T dd3 = this->Damp_mr(f3, param, lambda);
	
	coeffs[0] =this->calculate_delta_parameter_0(f1,f2,f3,v1,v2,v3,dd1,dd3,M);
	coeffs[1] =this->calculate_delta_parameter_1(f1,f2,f3,v1,v2,v3,dd1,dd3,M);
	coeffs[2] =this->calculate_delta_parameter_2(f1,f2,f3,v1,v2,v3,dd1,dd3,M);
	coeffs[3] =this->calculate_delta_parameter_3(f1,f2,f3,v1,v2,v3,dd1,dd3,M);
	coeffs[4] =this->calculate_delta_parameter_4(f1,f2,f3,v1,v2,v3,dd1,dd3,M);
}

/*!\brief Calculates the delta_0 component
 *
 * Solved in Mathematica and imported to C
 */
template <class T>
T IMRPhenomD<T>::calculate_delta_parameter_0(T f1,T f2,T f3,T v1,
				T v2,T v3,T dd1,T dd3,T M)
{
	return (-(dd3*pow(f1,2)*pow(f1 - f2,2)*f2*(f1 - f3)*(f2 - f3)*f3) + 
     dd1*f1*(f1 - f2)*f2*(f1 - f3)*pow(f2 - f3,2)*pow(f3,2) - 
     4*pow(f1,2)*pow(f2,3)*pow(f3,2)*v1 + 3*f1*pow(f2,4)*pow(f3,2)*v1 + 
     8*pow(f1,2)*pow(f2,2)*pow(f3,3)*v1 - 4*f1*pow(f2,3)*pow(f3,3)*v1 - 
     pow(f2,4)*pow(f3,3)*v1 - 4*pow(f1,2)*f2*pow(f3,4)*v1 - 
     f1*pow(f2,2)*pow(f3,4)*v1 + 2*pow(f2,3)*pow(f3,4)*v1 + 
     2*f1*f2*pow(f3,5)*v1 - pow(f2,2)*pow(f3,5)*v1 + 
     pow(f1,5)*pow(f3,2)*v2 - 3*pow(f1,4)*pow(f3,3)*v2 + 
     3*pow(f1,3)*pow(f3,4)*v2 - pow(f1,2)*pow(f3,5)*v2 + 
     pow(f1,5)*pow(f2,2)*v3 - 2*pow(f1,4)*pow(f2,3)*v3 + 
     pow(f1,3)*pow(f2,4)*v3 - 2*pow(f1,5)*f2*f3*v3 + 
     pow(f1,4)*pow(f2,2)*f3*v3 + 4*pow(f1,3)*pow(f2,3)*f3*v3 - 
     3*pow(f1,2)*pow(f2,4)*f3*v3 + 4*pow(f1,4)*f2*pow(f3,2)*v3 - 
     8*pow(f1,3)*pow(f2,2)*pow(f3,2)*v3 + 
     4*pow(f1,2)*pow(f2,3)*pow(f3,2)*v3)/
   (pow(f1 - f2,2)*pow(f1 - f3,3)*pow(f2 - f3,2));
}

/*!\brief Calculates the delta_1 component
 *
 * Solved in Mathematica and imported to C
 */
template <class T>
T IMRPhenomD<T>::calculate_delta_parameter_1(T f1,T f2,T f3,T v1,
				T v2,T v3,T dd1,T dd3,T M)
{
	return (dd3*f1*pow(f1 - f2,2)*(f1 - f3)*(f2 - f3)*(2*f2*f3 + f1*(f2 + f3)) + f3*
      (-(dd1*(f1 - f2)*(f1 - f3)*pow(f2 - f3,2)*(f2*f3 + f1*(2*f2 + f3))) - 
        2*f1*(pow(f3,4)*(v1 - v2) + f1*(2*pow(f3,3)*(-v1 + v2) - 4*pow(f2,3)*(v1 - v3) + 6*pow(f2,2)*f3*(v1 - v3)) + 3*pow(f2,4)*(v1 - v3) + pow(f1,4)*(v2 - v3) + 
           4*pow(f2,3)*f3*(-v1 + v3) + 2*pow(f1,3)*f3*(-v2 + v3))))/(pow(f1 - f2,2)*pow(f1 - f3,3)*pow(f2 - f3,2)*M);
}
/*!\brief Calculates the delta_2 component
 *
 * Solved in Mathematica and imported to C
 */
template <class T>
T IMRPhenomD<T>::calculate_delta_parameter_2(T f1,T f2,T f3,T v1,
				T v2,T v3,T dd1,T dd3,T M)
{
	return (-(dd3*pow(f1 - f2,2)*(f1 - f3)*(f2 - f3)*(pow(f1,2) + f2*f3 + 2.*f1*(f2 + f3))) + 
     dd1*(f1 - f3)*pow(f2 - f3,2)*(-(f2*f3*(2*f2 + f3)) + pow(f1,2)*(f2 + 2*f3) + f1*(-pow(f2,2) + pow(f3,2))) - 4*pow(f1,2)*pow(f2,3)*v1 + 3*f1*pow(f2,4)*v1 - 
     4*f1*pow(f2,3)*f3*v1 + 3*pow(f2,4)*f3*v1 + 12*pow(f1,2)*f2*pow(f3,2)*v1 - 4*pow(f2,3)*pow(f3,2)*v1 - 8*pow(f1,2)*pow(f3,3)*v1 + f1*pow(f3,4)*v1 + 
     pow(f3,5)*v1 + pow(f1,5)*v2 + pow(f1,4)*f3*v2 - 8*pow(f1,3)*pow(f3,2)*v2 + 8*pow(f1,2)*pow(f3,3)*v2 - f1*pow(f3,4)*v2 - pow(f3,5)*v2 - pow(f1,5)*v3 + 
     4*pow(f1,2)*pow(f2,3)*v3 - 3*f1*pow(f2,4)*v3 - pow(f1,4)*f3*v3 + 4*f1*pow(f2,3)*f3*v3 - 3*pow(f2,4)*f3*v3 + 8*pow(f1,3)*pow(f3,2)*v3 - 
     12*pow(f1,2)*f2*pow(f3,2)*v3 + 4*pow(f2,3)*pow(f3,2)*v3)/(pow(f1 - f2,2)*pow(f1 - f3,3)*pow(f2 - f3,2)*pow(M,2));
}

/*!\brief Calculates the delta_3 component
 *
 * Solved in Mathematica and imported to C
 */
template <class T>
T IMRPhenomD<T>::calculate_delta_parameter_3(T f1,T f2,T f3,T v1,
				T v2,T v3,T dd1,T dd3,T M)
{

	return (dd3*pow(f1 - f2,2)*(f1 - f3)*(f2 - f3)*(2*f1 + f2 + f3) - dd1*(f1 - f2)*(f1 - f3)*pow(f2 - f3,2)*(f1 + f2 + 2*f3) + 
     2*(pow(f3,4)*(-v1 + v2) + 2*pow(f1,2)*pow(f2 - f3,2)*(v1 - v3) + 2*pow(f2,2)*pow(f3,2)*(v1 - v3) + 2*pow(f1,3)*f3*(v2 - v3) + pow(f2,4)*(-v1 + v3) + 
        pow(f1,4)*(-v2 + v3) + 2*f1*f3*(pow(f3,2)*(v1 - v2) + pow(f2,2)*(v1 - v3) + 2*f2*f3*(-v1 + v3))))/(pow(f1 - f2,2)*pow(f1 - f3,3)*pow(f2 - f3,2)*pow(M,3));
}

/*!\brief Calculates the delta_4 component
 *
 * Solved in Mathematica and imported to C
 */
template <class T>
T IMRPhenomD<T>::calculate_delta_parameter_4(T f1,T f2,T f3,T v1,
				T v2,T v3,T dd1,T dd3,T M)
{

	return (-(dd3*pow(f1 - f2,2)*(f1 - f3)*(f2 - f3)) + dd1*(f1 - f2)*(f1 - f3)*pow(f2 - f3,2) - 3*f1*pow(f2,2)*v1 + 2*pow(f2,3)*v1 + 6*f1*f2*f3*v1 - 3*pow(f2,2)*f3*v1 - 
     3*f1*pow(f3,2)*v1 + pow(f3,3)*v1 + pow(f1,3)*v2 - 3*pow(f1,2)*f3*v2 + 3*f1*pow(f3,2)*v2 - pow(f3,3)*v2 - pow(f1,3)*v3 + 3*f1*pow(f2,2)*v3 - 2*pow(f2,3)*v3 + 
     3*pow(f1,2)*f3*v3 - 6*f1*f2*f3*v3 + 3*pow(f2,2)*f3*v3)/(pow(f1 - f2,2)*pow(f1 - f3,3)*pow(f2 - f3,2)*pow(M,4));
}
template class IMRPhenomD<double>;
template class IMRPhenomD<adouble>;
