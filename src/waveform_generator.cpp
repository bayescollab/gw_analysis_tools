#include <iostream>
#include "waveform_generator.h"
#include "TaylorT2.h"
#include "IMRPhenomD.h"
#include "IMRPhenomP.h"
#include "IMRPhenomPv3.h"
#include "ppE_IMRPhenomD.h"
#include "ppE_IMRPhenomD_NRT.h"
#include "ppE_IMRPhenomP.h"
#include "gIMRPhenomD.h"
#include "IMRPhenomD_NRT.h"
#include "EA_IMRPhenomD_NRT.h"
#include "EA_IMRPhenomD.h"
#include "IMRPhenomP_NRT.h"
#include "ppE_utilities.h"
#include "gIMRPhenomP.h"
#include "util.h"
#include <complex>
#include <time.h>
#include <gsl/gsl_randist.h>
#include <adolc/adouble.h>
#include <omp.h>
using namespace std;

/*!\file
 *File that handles the construction of the (2,2) waveform as described by IMRPhenomD by Khan et. al.
 *
 * Builds a waveform for given DETECTOR FRAME parameters
 */

template<class T>
int time_waveform(T *times, /**< double array of frequencies for the waveform to be evaluated at*/
	int length,/**<integer length of all the arrays*/
	waveform_polarizations<T> *wp,/**< [out] Output waveforms by polarization*/
	string generation_method,/**<String that corresponds to the generation method - MUST BE SPELLED EXACTLY*/
	gen_params_base<T> *parameters/**<structure containing all the source parameters*/
	)
{
	int status=1;
	bool NSflag1 = parameters->NSflag1;
	bool NSflag2 = parameters->NSflag2;


	source_parameters<T> params;

	std::string local_method = prep_source_parameters(&params, parameters,generation_method);
	if(local_method.find("Taylor")!=std::string::npos)
	{
		//std::complex<T> ci = std::complex<T>(cos(params.incl_angle),0);
		if(local_method == "TaylorT2")
		{
			TaylorT2<T> T2model;
			status = T2model.construct_waveform(times, length, wp->hplus,wp->hcross, &params);

		}
	}
	//if(check_extra_polarizations(generation_method))
	//{
	//	//TESTING MUST FIX
	//	for (int i =0;i < length; i++)
	//	{
	//		wp->hx[i] = wp->hplus[i];
	//		wp->hy[i] = wp->hplus[i];
	//		wp->hb[i] = wp->hplus[i];
	//		wp->hl[i] = wp->hplus[i];
	//	}
	//}
	cleanup_source_parameters(&params,generation_method);

	return status ;
}



/*!\brief Function to produce the plus/cross polarizations of an quasi-circular binary
 *
 * By using the structure parameter, the function is allowed to be more flexible in using different
 * method of waveform generation - not all methods use the same parameters
 *
 * This puts the responsibility on the user to pass the necessary parameters
 *
 *
 * *NEED TO OUTLINE OPTIONS FOR EACH METHOD IN DEPTH*
 *
 *
 * NEW PHASE OPTIONS for
 *
 * PHENOMD ONLY:
 *
 * If phic is assigned, the reference frequency and reference phase are IGNORED.
 *
 * If Phic is unassigned, a reference phase AND a reference frequency are looked for.If no options are found, both are set to 0.
 *
 * If tc is assigned, it is used.
 *
 * If tc is unassigned, the waveform is shifted so the merger happens at 0.
 *
 * PhenomPv2:
 *
 * PhiRef and f_ref are required, phic is not an option.
 *
 * tc, if specified, is used with the use of interpolation. If not, tc is set such that coalescence happens at t=0
 */
template<class T>
int fourier_waveform(T *frequencies, /**< double array of frequencies for the waveform to be evaluated at*/
			int length,/**<integer length of all the arrays*/
			waveform_polarizations<T> *wp,/**< [out] Output waveforms by polarization*/
			//std::complex<T> *waveform_plus, /**< [out] complex array for the output plus polarization waveform*/
			//std::complex<T> *waveform_cross, /**< [out] complex array for the output cross polarization waveform*/
			string generation_method,/**<String that corresponds to the generation method - MUST BE SPELLED EXACTLY*/
			gen_params_base<T> *parameters/**<structure containing all the source parameters*/
			)
{
	int status=1;
	bool NSflag1 = parameters->NSflag1;
	bool NSflag2 = parameters->NSflag2;


	/*Eventually, this will be where NS specific quantities are defined*/
	//if (NSflag1 || NSflag2)
	//{
	//	cout<<"NS waveforms still under develpment - BH only"<<endl;
	//	return 0;
	//}
	/* Convert all dimensionful quantities to seconds and build all needed source quantities once*/
	source_parameters<T> params;
	//params = params.populate_source_parameters(parameters);

	std::string local_method = prep_source_parameters(&params, parameters,generation_method);
	if(local_method.find("IMRPhenomD")!=std::string::npos)
	{
		std::complex<T> ci = std::complex<T>(cos(params.incl_angle),0);
		std::complex<T> si = std::complex<T>(sin(params.incl_angle),0);
		std::complex<T> s2i = std::complex<T>(sin(2*params.incl_angle),0);
		bool restrictedEval = true;
		if(local_method == "ppE_IMRPhenomD_Inspiral")
		{
			ppE_IMRPhenomD_Inspiral<T> ppemodeld;
			status = ppemodeld.construct_waveform(frequencies, length, wp->hplus, &params);

		}
		else if(local_method == "ppE_IMRPhenomD_IMR")
		{
			ppE_IMRPhenomD_IMR<T> ppemodeld;
			status = ppemodeld.construct_waveform(frequencies, length, wp->hplus, &params);

		}
		else if(local_method == "ppE_IMRPhenomD_NRT_Inspiral")
		{
			ppE_IMRPhenomD_NRT_Inspiral<T> ppemodeld;
			status = ppemodeld.construct_waveform(frequencies, length, wp->hplus, &params);

		}
		else if(local_method == "ppE_IMRPhenomD_NRT_IMR")
		{
			ppE_IMRPhenomD_NRT_IMR<T> ppemodeld;
			status = ppemodeld.construct_waveform(frequencies, length, wp->hplus, &params);

		}
		else if(local_method == "gIMRPhenomD")
		{
			gIMRPhenomD<T> gmodeld;
			status = gmodeld.construct_waveform(frequencies, length, wp->hplus, &params);

		}
		else if(local_method == "IMRPhenomD_NRT")
		  {
		    IMRPhenomD_NRT<T> modeldNRT;
		    status = modeldNRT.construct_waveform(frequencies, length, wp->hplus, &params);

		  }
		else if(local_method == "EA_IMRPhenomD_NRT")
		  {
		    restrictedEval = false;
		    EA_IMRPhenomD_NRT<T> EAmodeldNRT;
		    status = EAmodeldNRT.EA_construct_waveform(frequencies, length, wp, &params);

		  }
		else if(local_method == "EA_IMRPhenomD")
		  {
		    restrictedEval = false;
		    EA_IMRPhenomD<T> EAmodeld;
		    status = EAmodeld.EA_construct_waveform(frequencies, length, wp, &params);

		  }

		else{
			IMRPhenomD<T> modeld;
			status = modeld.construct_waveform(frequencies, length, wp->hplus, &params);
		}
		if(restrictedEval){
			//for (int i =0 ; i < length; i++){
			//	wp->hcross[i] = ci*std::complex<T>(0,-1) * wp->hplus[i];
			//	wp->hplus[i] = wp->hplus[i]* std::complex<T>(.5,0) *(std::complex<T>(1,0)+ci*ci);
			//}
			for (int i =0 ; i < length; i++){
				wp->hcross[i] = std::complex<T>(0,-1) * wp->hplus[i];
				wp->hplus[i] = wp->hplus[i];
			}
		}
		/*! Handle the inclination angle dependence*/
		/*! Iota is weird, because it's a extrinsic parameter, but one that has dependence on frequency (potentially)*/
		/*! I think it should be handled separately because of complications like LISA response, etc*/
		for (int i =0 ; i < length; i++){
			wp->hcross[i] *= ci;
			wp->hplus[i] *= std::complex<T>(.5,0) *(std::complex<T>(1,0)+ci*ci);
		}
		/*! Iota dependence for the extra four polarizations that 
		 * don't exist within GR is handled in the specific waveform 
		 * template file. For instance, for Einstein AEther, it is in 
		 * the EA_IMRPhenomD_NRT.cpp file. 
		 */
		/*
		if(wp->active_polarizations[2]){
			for (int i =0 ; i < length; i++){
			  wp->hx[i] *= s2i;
			}
		}
		if(wp->active_polarizations[3]){
			for (int i =0 ; i < length; i++){
				wp->hy[i] *= si;
			}
		}
		if(wp->active_polarizations[4]){
			for (int i =0 ; i < length; i++){
				wp->hb[i] *= si*si;
			}
		}
		if(wp->active_polarizations[5]){
			for (int i =0 ; i < length; i++){
				wp->hl[i] *= si*si;
			}	
		}*/
	}
	else if(local_method.find("IMRPhenomPv2")!=std::string::npos)
	{
		if(local_method == "ppE_IMRPhenomPv2_Inspiral")
		{
			ppE_IMRPhenomPv2_Inspiral<T> ppemodel;
			status = ppemodel.construct_waveform(frequencies, length, wp->hplus, wp->hcross,&params);

		}
		else if(local_method == "ppE_IMRPhenomPv2_IMR")
		{
			ppE_IMRPhenomPv2_IMR<T> ppemodel;
			status = ppemodel.construct_waveform(frequencies, length, wp->hplus,wp->hcross, &params);

		}
		else if(local_method == "gIMRPhenomPv2")
		{
			gIMRPhenomPv2<T> gmodel;
			status = gmodel.construct_waveform(frequencies, length, wp->hplus,wp->hcross, &params);

		}
		else if(local_method == "IMRPhenomPv2_NRT")
		  {
		    IMRPhenomPv2_NRT<T> modelNRT;
		    status = modelNRT.construct_waveform(frequencies, length, wp->hplus, wp->hcross, &params);

		  }
		else{
			IMRPhenomPv2<T> model;
			status = model.construct_waveform(frequencies, length, wp->hplus, wp->hcross, &params);
		}
		std::complex<T> tempPlus,tempCross;
		std::complex<T> c2z = std::complex<T>(cos(2.*params.zeta_polariz));
		std::complex<T> s2z = std::complex<T>(sin(2.*params.zeta_polariz));
		for (int i =0;i < length; i++)
		{
			tempPlus = wp->hplus[i];
			tempCross = wp->hcross[i];
			wp->hplus[i] = c2z*tempPlus+s2z*tempCross;
			wp->hcross[i] = c2z*tempCross-s2z*tempPlus;
		}
	}
	else if(local_method.find("IMRPhenomPv3")!=std::string::npos)
	{
		if (local_method == "IMRPhenomPv3")
		{
			IMRPhenomPv3<T> model;
			status = model.construct_waveform(frequencies, length, wp->hplus, wp->hcross, &params);
		}
	}

	//Catch all for any modifications not captured in ppE formalism like extra polarizations
	extra_modifications(generation_method, parameters,&params, wp,frequencies,length);

	//if(check_extra_polarizations(generation_method))
	//{
	//	//TESTING MUST FIX
	//	for (int i =0;i < length; i++)
	//	{
	//		wp->hx[i] = wp->hplus[i];
	//		wp->hy[i] = wp->hplus[i];
	//		wp->hb[i] = wp->hplus[i];
	//		wp->hl[i] = wp->hplus[i];
	//	}
	//}
	cleanup_source_parameters(&params,generation_method);

	return status ;
}







int fourier_waveform(double *frequencies, /**< double array of frequencies for the waveform to be evaluated at*/
			int length,/**<integer length of all the arrays*/
			double *waveform_plus_real, /**< complex array for the output waveform*/
			double *waveform_plus_imag, /**< complex array for the output waveform*/
			double *waveform_cross_real, /**< complex array for the output waveform*/
			double *waveform_cross_imag, /**< complex array for the output waveform*/
			double *waveform_x_real, /**< complex array for the output waveform*/
			double *waveform_x_imag, /**< complex array for the output waveform*/
			double *waveform_y_real, /**< complex array for the output waveform*/
			double *waveform_y_imag, /**< complex array for the output waveform*/
			double *waveform_b_real, /**< complex array for the output waveform*/
			double *waveform_b_imag, /**< complex array for the output waveform*/
			double *waveform_l_real, /**< complex array for the output waveform*/
			double *waveform_l_imag, /**< complex array for the output waveform*/
			string generation_method,/**<String that corresponds to the generation method - MUST BE SPELLED EXACTLY*/
			gen_params *parameters/**<structure containing all the source parameters*/
			)
{
	//std::complex<double> *waveform_plus = (std::complex<double> *)malloc(sizeof(std::complex<double>)*length);
	//std::complex<double> *waveform_cross = (std::complex<double> *)malloc(sizeof(std::complex<double>)*length);
	std::complex<double> *waveform_plus = new std::complex<double>[length];
	std::complex<double> *waveform_cross = new std::complex<double>[length];
	std::complex<double> *waveform_x;
	std::complex<double> *waveform_y;
	std::complex<double> *waveform_b;
	std::complex<double> *waveform_l;
	if(waveform_x_real){
		waveform_x = new std::complex<double>[length];
	}
	if(waveform_y_real){
		waveform_y = new std::complex<double>[length];
	}
	if(waveform_b_real){
		waveform_b = new std::complex<double>[length];
	}
	if(waveform_l_real){
		waveform_l = new std::complex<double>[length];
	}
	waveform_polarizations<double> *wp = new waveform_polarizations<double>;
	wp->hplus = waveform_plus;
	wp->hcross = waveform_cross;
	wp->hx = waveform_x;
	wp->hy = waveform_y;
	wp->hb = waveform_b;
	wp->hl = waveform_l;
	//int status = fourier_waveform(frequencies, length, waveform_plus,waveform_cross, generation_method, parameters);
	int status = fourier_waveform(frequencies, length, wp, generation_method, parameters);
	for(int l = 0 ; l<length;l++){
		waveform_plus_real[l] = std::real(waveform_plus[l]);
		waveform_plus_imag[l] = std::imag(waveform_plus[l]);
		waveform_cross_real[l] = std::real(waveform_cross[l]);
		waveform_cross_imag[l] = std::imag(waveform_cross[l]);
	}
	if(waveform_x_real){
		for(int l = 0 ; l<length;l++){
			waveform_x_real[l] = std::real(waveform_x[l]);
			waveform_x_imag[l] = std::imag(waveform_x[l]);
		}
	}
	if(waveform_y_real){
		for(int l = 0 ; l<length;l++){
			waveform_y_real[l] = std::real(waveform_y[l]);
			waveform_y_imag[l] = std::imag(waveform_y[l]);
		}
	}
	if(waveform_b_real){
		for(int l = 0 ; l<length;l++){
			waveform_b_real[l] = std::real(waveform_b[l]);
			waveform_b_imag[l] = std::imag(waveform_b[l]);
		}
	}
	if(waveform_l_real){
		for(int l = 0 ; l<length;l++){
			waveform_l_real[l] = std::real(waveform_l[l]);
			waveform_l_imag[l] = std::imag(waveform_l[l]);
		}
	}
	delete [] waveform_plus;
	delete [] waveform_cross;
	if(waveform_x_real){
		delete [] waveform_x;
	}
	if(waveform_y_real){
		delete [] waveform_y;
	}
	if(waveform_b_real){
		delete [] waveform_b;
	}
	if(waveform_l_real){
		delete [] waveform_l;
	}
	delete wp;
	return status;
}





/*!\brief Function to produce the (2,2) mode of an quasi-circular binary
 *
 * By using the structure parameter, the function is allowed to be more flexible in using different
 * method of waveform generation - not all methods use the same parameters
 *
 */
int fourier_waveform(double *frequencies, /**< double array of frequencies for the waveform to be evaluated at*/
			int length,/**<integer length of all the arrays*/
			std::complex<double> *waveform, /**< complex array for the output waveform*/
			string generation_method,/**<String that corresponds to the generation method - MUST BE SPELLED EXACTLY*/
			gen_params *parameters/**<structure containing all the source parameters*/
			)
{

	int status=1;
	bool NSflag1 = parameters->NSflag1;
	bool NSflag2 = parameters->NSflag2;

	/*Eventually, this will be where NS specific quantities are defined*/
	//if (NSflag1 || NSflag2)
	//{
	//	cout<<"NS waveforms still under develpment - BH only"<<endl;
	//	return 0;
	//}
	/* Convert all dimensionful quantities to seconds and build all needed source quantities once*/
	double mass1 = parameters->mass1;
	double mass2 = parameters->mass2;
	double Luminosity_Distance = parameters->Luminosity_Distance;
	double *spin1 = parameters->spin1;
	double *spin2 = parameters->spin2;
	double t_c = parameters->tc;
	source_parameters<double> params;
	//params = params.populate_source_parameters(mass1, mass2, Luminosity_Distance, spin1, spin2, phi_c,t_c);
	//params = params.populate_source_parameters(parameters);
	params.populate_source_parameters(parameters);

	params.f_ref = parameters->f_ref;
	params.phiRef = parameters->phiRef;
	params.phi = parameters->phi;
	params.theta = parameters->theta;
	params.cosmology = parameters->cosmology;
	params.shift_time = parameters->shift_time;
	params.shift_phase = parameters->shift_phase;
	params.NSflag1 = parameters->NSflag1;
	params.NSflag2 = parameters->NSflag2;
	params.dep_postmerger = parameters->dep_postmerger;
	if(generation_method == "IMRPhenomD")
	{
		IMRPhenomD<double> modeld;
		status = modeld.construct_waveform(frequencies, length, waveform, &params);
	}
	else if(generation_method == "ppE_IMRPhenomD_Inspiral")
	{
		ppE_IMRPhenomD_Inspiral<double> ppemodeld;
		params.betappe = parameters->betappe;
		params.bppe = parameters->bppe;
		params.Nmod = parameters->Nmod;
		status = ppemodeld.construct_waveform(frequencies, length, waveform, &params);
	}
	//else if(generation_method == "_dCS_IMRPhenomD")
	//{
	//	bool local_spline = false;
	//	dCS_IMRPhenomD<double> ppemodeld;
	//	params.betappe = parameters->betappe;
	//	params.Nmod = 1;
	//	int tempbppe[params.Nmod] = {-1};
	//	params.bppe = tempbppe;
	//	double temp[params.Nmod] ;
	//	for( int i = 0; i < params.Nmod; i++)
	//		temp[i] = params.betappe[i];
	//	status = ppemodeld.construct_waveform(frequencies, length, waveform, &params);
	//	for( int i = 0; i < params.Nmod; i++)
	//		parameters->betappe[i] = temp[i];
	//}
	//else if(generation_method == "EdGB_IMRPhenomD")
	//{
	//	EdGB_IMRPhenomD<double> ppemodeld;
	//	params.betappe = parameters->betappe;
	//	params.Nmod = 1;
	//	int tempbppe[params.Nmod] = {-7};
	//	params.bppe = tempbppe;
	//	double temp[params.Nmod] ;
	//	for( int i = 0; i < params.Nmod; i++)
	//		temp[i] = params.betappe[i];
	//	status = ppemodeld.construct_waveform(frequencies, length, waveform, &params);
	//	for( int i = 0; i < params.Nmod; i++)
	//		parameters->betappe[i] = temp[i];
	//}
	else if(generation_method == "ppE_IMRPhenomD_IMR")
	{
		ppE_IMRPhenomD_IMR<double> ppemodeld;
		params.betappe = parameters->betappe;
		params.bppe = parameters->bppe;
		params.Nmod = parameters->Nmod;
		status = ppemodeld.construct_waveform(frequencies, length, waveform, &params);
	}
	else if(generation_method == "IMRPhenomD_NRT")
	  {
	    IMRPhenomD_NRT<double> modeldNRT;
	    status = modeldNRT.construct_waveform(frequencies, length, waveform, &params);

	  }
	//else if(generation_method == "IMRPhenomPv2")
	//{
	//	IMRPhenomPv2<double> modeld;
	//	//Initialize Pv2 specific params
	//	status = modeld.construct_waveform(frequencies, length, waveform, &params);
	//}

	return status ;
}

int fourier_waveform(double *frequencies, /**< double array of frequencies for the waveform to be evaluated at*/
			int length,/**<integer length of all the arrays*/
			double *waveform_real, /**< complex array for the output waveform*/
			double *waveform_imag, /**< complex array for the output waveform*/
			string generation_method,/**<String that corresponds to the generation method - MUST BE SPELLED EXACTLY*/
			gen_params *parameters/**<structure containing all the source parameters*/
			)
{
	std::complex<double> *waveform = (std::complex<double> *)malloc(sizeof(std::complex<double>)*length);
	int status = fourier_waveform(frequencies, length, waveform, generation_method, parameters);
	for(int i = 0 ;i<length;i++)
	{
		waveform_real[i] = real(waveform[i]);
		waveform_imag[i] = imag(waveform[i]);
	}
	free(waveform);
	return status;
}

/*!\brief Function to produce the amplitude of the (2,2) mode of an quasi-circular binary
 *
 * By using the structure parameter, the function is allowed to be more flexible in using different
 * method of waveform generation - not all methods use the same parameters
 */
template<class T>
int fourier_amplitude(T *frequencies, /**< double array of frequencies for the waveform to be evaluated at*/
			int length,/**<integer length of all the arrays*/
			T *amplitude, /**< output array for the amplitude*/
			string generation_method,/**<String that corresponds to the generation method - MUST BE SPELLED EXACTLY*/
			gen_params_base<T> *parameters
			)
{
	int status=1;
	bool NSflag1 = parameters->NSflag1;
	bool NSflag2 = parameters->NSflag2;

	/*Eventually, this will be where NS specific quantities are defined*/
	//if (NSflag1 || NSflag2)
	//{
	//	cout<<"NS waveforms still under develpment - BH only"<<endl;
	//	return 0;
	//}
	/* Convert all dimensionful quantities to seconds and build all needed source quantities once*/
	source_parameters<T> params;
	//params = params.populate_source_parameters(parameters);
	params.populate_source_parameters(parameters);

	std::string local_method = prep_source_parameters(&params, parameters,generation_method);
	if(local_method.find("IMRPhenomD")!=std::string::npos)
	{
		std::complex<T> ci = std::complex<T>(cos(params.incl_angle),0);
		if(local_method == "ppE_IMRPhenomD_Inspiral")
		{
			ppE_IMRPhenomD_Inspiral<T> ppemodeld;
			status = ppemodeld.construct_amplitude(frequencies, length, amplitude, &params);

		}
		else if(local_method == "ppE_IMRPhenomD_IMR")
		{
			ppE_IMRPhenomD_IMR<T> ppemodeld;
			status = ppemodeld.construct_amplitude(frequencies, length, amplitude, &params);

		}
		else if(local_method == "gIMRPhenomD")
		{
			gIMRPhenomD<T> gmodeld;
			status = gmodeld.construct_amplitude(frequencies, length, amplitude, &params);

		}
		else if(local_method == "IMRPhenomD_NRT")
		  {
		    IMRPhenomD_NRT<T> modeldNRT;
		    status = modeldNRT.construct_amplitude(frequencies, length, amplitude, &params);

		  }
		else if(local_method == "EA_IMRPhenomD_NRT")
		  {
		    EA_IMRPhenomD_NRT<T> EAmodeldNRT;
		    status = EAmodeldNRT.construct_amplitude(frequencies, length, amplitude, &params);

		  }
		else if(local_method == "EA_IMRPhenomD")
		  {
		    EA_IMRPhenomD<T> EAmodeld;
		    status = EAmodeld.construct_amplitude(frequencies, length, amplitude, &params);

		  }
		else{
			IMRPhenomD<T> modeld;
			status = modeld.construct_amplitude(frequencies, length, amplitude, &params);

		}
	}
	cleanup_source_parameters(&params,generation_method);





	//params.f_ref = parameters->f_ref;
	//params.phiRef = parameters->phiRef;
	//params.cosmology = parameters->cosmology;
	//params.shift_time = parameters->shift_time;
	//params.shift_phase = parameters->shift_phase;
	//params.NSflag1 = parameters->NSflag1;
	//params.NSflag2 = parameters->NSflag2;
	//params.dep_postmerger = parameters->dep_postmerger;
	//if(generation_method == "IMRPhenomD")
	//{
	//	IMRPhenomD<T> modeld;
	//	status = modeld.construct_amplitude(frequencies, length, amplitude, &params);
	//}
	//else if(generation_method == "ppE_IMRPhenomD_Inspiral")
	//{
	//	params.bppe = parameters->bppe;
	//	params.Nmod = parameters->Nmod;
	//	params.betappe = parameters->betappe;
	//	ppE_IMRPhenomD_Inspiral<T> ppemodeld;
	//	status = ppemodeld.construct_amplitude(frequencies, length, amplitude, &params);
	//}
	//else if(generation_method == "_dCS_IMRPhenomD")
	//{
	//	dCS_IMRPhenomD<T> ppemodeld;
	//	status = ppemodeld.construct_amplitude(frequencies, length, amplitude, &params);
	//}
	//else if(generation_method == "EdGB_IMRPhenomD")
	//{
	//	EdGB_IMRPhenomD<T> ppemodeld;
	//	status = ppemodeld.construct_amplitude(frequencies, length, amplitude, &params);
	//}
	//else if(generation_method == "ppE_IMRPhenomD_IMR")
	//{
	//	params.bppe = parameters->bppe;
	//	params.Nmod = parameters->Nmod;
	//	params.betappe = parameters->betappe;
	//	ppE_IMRPhenomD_IMR<T> ppemodeld;
	//	status = ppemodeld.construct_amplitude(frequencies, length, amplitude, &params);
	//}
	//else if(generation_method == "gIMRPhenomD")
	//{
	//	gIMRPhenomD<T> gmodeld;
	//	params.delta_phi = parameters->delta_phi;
	//	params.delta_sigma = parameters->delta_sigma;
	//	params.delta_beta = parameters->delta_beta;
	//	params.delta_alpha = parameters->delta_alpha;
	//	params.phii = parameters->phii;
	//	params.sigmai = parameters->sigmai;
	//	params.betai = parameters->betai;
	//	params.alphai = parameters->alphai;
	//	params.Nmod_phi = parameters->Nmod_phi;
	//	params.Nmod_sigma = parameters->Nmod_sigma;
	//	params.Nmod_beta = parameters->Nmod_beta;
	//	params.Nmod_alpha = parameters->Nmod_alpha;
	//	status = gmodeld.construct_amplitude(frequencies, length, amplitude, &params);
	//}

	return status ;
}
template int fourier_amplitude<double>(double *, int , double * ,std::string, gen_params_base<double> *);
template int fourier_amplitude<adouble>(adouble *, int , adouble * ,std::string, gen_params_base<adouble> *);
/*!\brief Function to produce the phase of the (2,2) mode of an quasi-circular binary
 *
 * By using the structure parameter, the function is allowed to be more flexible in using different
 * method of waveform generation - not all methods use the same parameters
 */
template<class T>
int fourier_phase(T *frequencies, /**<double array of frequencies for the waveform to be evaluated at*/
			int length,/**<integer length of all the arrays*/
			T *phase, /**<output array for the phase*/
			string generation_method,/**<String that corresponds to the generation method - MUST BE SPELLED EXACTLY*/
			gen_params_base<T> *parameters
			)
{
	std::cout<<"Possible error in negative sign! Deprecated, please use the version with (+) and (x) separate"<<std::endl;
	int status=1;
	bool NSflag1 = parameters->NSflag1;
	bool NSflag2 = parameters->NSflag2;

	/*Eventually, this will be where NS specific quantities are defined*/
	//if (NSflag1 || NSflag2)
	//{
	//	cout<<"NS waveforms still under develpment - BH only"<<endl;
	//	return 0;
	//}
	/* Convert all dimensionful quantities to seconds and build all needed source quantities once*/
	source_parameters<T> params;
	//params = params.populate_source_parameters(parameters);
	std::string local_method = prep_source_parameters(&params, parameters,generation_method);
	//params.populate_source_parameters(parameters);
	//params.f_ref = parameters->f_ref;
	//params.phiRef = parameters->phiRef;
	//params.cosmology = parameters->cosmology;
	//params.shift_time = parameters->shift_time;
	//params.shift_phase = parameters->shift_phase;
	//params.NSflag1 = parameters->NSflag1;
	//params.NSflag2 = parameters->NSflag2;
	//params.dep_postmerger = parameters->dep_postmerger;

	if(local_method == "IMRPhenomD")
	{
		IMRPhenomD<T> modeld;
		status = modeld.construct_phase(frequencies, length, phase, &params);
		for(int i = 0 ; i<length; i++){
				phase[i]*= (T)(-1.);
		}
	}
	else if(local_method == "ppE_IMRPhenomD_Inspiral")
	{
		params.betappe = parameters->betappe;
		params.bppe = parameters->bppe;
		params.Nmod = parameters->Nmod;
		ppE_IMRPhenomD_Inspiral<T> ppemodeld;
		status = ppemodeld.construct_phase(frequencies, length, phase, &params);
		for(int i = 0 ; i<length; i++){
				phase[i]*= (T)(-1.);
		}
	}
	//else if(generation_method == "_dCS_IMRPhenomD")
	//{
	//	bool local_spline = false;
	//	params.betappe = parameters->betappe;
	//	params.Nmod = 1;
	//	int tempbppe[params.Nmod] = {-1};
	//	params.bppe = tempbppe;
	//	T temp[params.Nmod] ;
	//	for( int i = 0; i < params.Nmod; i++)
	//		temp[i] = params.betappe[i];
	//	dCS_IMRPhenomD<T> ppemodeld;
	//	status = ppemodeld.construct_phase(frequencies, length, phase, &params);
	//	for(int i = 0 ; i<length; i++){
	//			phase[i]*= (T)(-1.);
	//	}
	//
	//	for( int i = 0; i < params.Nmod; i++)
	//		parameters->betappe[i] = temp[i];
	//}
	//else if(generation_method == "EdGB_IMRPhenomD")
	//{
	//	params.betappe = parameters->betappe;
	//	params.Nmod = 1;
	//	int tempbppe[params.Nmod] = {-7};
	//	params.bppe = tempbppe;
	//	T temp[params.Nmod] ;
	//	for( int i = 0; i < params.Nmod; i++)
	//		temp[i] = params.betappe[i];
	//	EdGB_IMRPhenomD<T> ppemodeld;
	//	status = ppemodeld.construct_phase(frequencies, length, phase, &params);
	//	for(int i = 0 ; i<length; i++){
	//			phase[i]*= (T)(-1.);
	//	}
	//	for( int i = 0; i < params.Nmod; i++)
	//		parameters->betappe[i] = temp[i];
	//}
	else if(local_method == "ppE_IMRPhenomD_IMR")
	{
		params.betappe = parameters->betappe;
		params.bppe = parameters->bppe;
		params.Nmod = parameters->Nmod;
		ppE_IMRPhenomD_IMR<T> ppemodeld;
		status = ppemodeld.construct_phase(frequencies, length, phase, &params);
		for(int i = 0 ; i<length; i++){
				phase[i]*= (T)(-1.);
		}
	}
	else if(local_method == "gIMRPhenomD")
	{
		gIMRPhenomD<T> gmodeld;
		params.delta_phi = parameters->delta_phi;
		params.delta_sigma = parameters->delta_sigma;
		params.delta_beta = parameters->delta_beta;
		params.delta_alpha = parameters->delta_alpha;
		params.phii = parameters->phii;
		params.sigmai = parameters->sigmai;
		params.betai = parameters->betai;
		params.alphai = parameters->alphai;
		params.Nmod_phi = parameters->Nmod_phi;
		params.Nmod_sigma = parameters->Nmod_sigma;
		params.Nmod_beta = parameters->Nmod_beta;
		params.Nmod_alpha = parameters->Nmod_alpha;
		status = gmodeld.construct_phase(frequencies, length, phase, &params);
		for(int i = 0 ; i<length; i++){
				phase[i]*= (T)(-1.);
		}
	}
	else if(local_method == "IMRPhenomD_NRT")
	  {
	    IMRPhenomD_NRT<T> modeldNRT;
	    params.tidal1 = parameters->tidal1; //Is this right?!?
	    params.tidal2 = parameters->tidal2;
	    status = modeldNRT.construct_phase(frequencies, length, phase, &params);
	  }
	else if(local_method == "EA_IMRPhenomD_NRT")
	  {
	    EA_IMRPhenomD_NRT<T> EAmodeldNRT;
	    params.tidal1 = parameters->tidal1;
	    params.tidal2 = parameters->tidal2;
	    status = EAmodeldNRT.construct_phase(frequencies, length, phase, &params);
	  }
	else if(local_method == "EA_IMRPhenomD")
	  {
	    EA_IMRPhenomD<T> EAmodeld;
	    status = EAmodeld.construct_phase(frequencies, length, phase, &params);
	  }

	return status ;
}
template int fourier_phase<double>(double *, int , double * , std::string, gen_params_base<double> *);
template int fourier_phase<adouble>(adouble *, int , adouble * , std::string,gen_params_base<adouble> *);
/*!\brief Function to produce the phase of the plus and cross mode of a quasi-circular binary
 *
 * By using the structure parameter, the function is allowed to be more flexible in using different
 * method of waveform generation - not all methods use the same parameters
 */
template<class T>
int fourier_phase(T *frequencies, /**<double array of frequencies for the waveform to be evaluated at*/
			int length,/**<integer length of all the arrays*/
			T *phase_plus, /**<output array for the phase*/
			T *phase_cross, /**<output array for the phase*/
			string generation_method,/**<String that corresponds to the generation method - MUST BE SPELLED EXACTLY*/
			gen_params_base<T> *parameters
			)
{
	int status=1;
	bool NSflag1 = parameters->NSflag1;
	bool NSflag2 = parameters->NSflag2;

	/*Eventually, this will be where NS specific quantities are defined*/
	//if (NSflag1 || NSflag2)
	//{
	//	cout<<"NS waveforms still under develpment - BH only"<<endl;
	//	return 0;
	//}
	/* Convert all dimensionful quantities to seconds and build all needed source quantities once*/
	source_parameters<T> params;

	std::string local_method = prep_source_parameters(&params, parameters,generation_method);
	if(local_method.find("IMRPhenomD")!=std::string::npos)
	{
		if(local_method == "ppE_IMRPhenomD_Inspiral")
		{
			ppE_IMRPhenomD_Inspiral<T> ppemodeld;
			status = ppemodeld.construct_phase(frequencies, length, phase_plus, &params);

		}
		else if(local_method == "ppE_IMRPhenomD_IMR")
		{
			ppE_IMRPhenomD_IMR<T> ppemodeld;
			status = ppemodeld.construct_phase(frequencies, length, phase_plus, &params);

		}
		else if(local_method == "gIMRPhenomD")
		{
			gIMRPhenomD<T> gmodeld;
			status = gmodeld.construct_phase(frequencies, length, phase_plus, &params);

		}
		else if(local_method == "IMRPhenomD_NRT")
		  {
		    IMRPhenomD_NRT<T> modeldNRT;
		    status = modeldNRT.construct_phase(frequencies, length, phase_plus, &params);

		  }
		else if(local_method == "EA_IMRPhenomD_NRT")
		  {
		    EA_IMRPhenomD_NRT<T> EAmodeldNRT;
		    status = EAmodeldNRT.construct_phase(frequencies, length, phase_plus, &params);

		  }
		else if(local_method == "EA_IMRPhenomD")
		  {
		    EA_IMRPhenomD<T> EAmodeld;
		    status = EAmodeld.construct_phase(frequencies, length, phase_plus, &params);

		  }
		else{
			IMRPhenomD<T> modeld;
			status = modeld.construct_phase(frequencies, length, phase_plus, &params);
		}
		for(int i = 0 ; i<length; i++){
			phase_cross[i] = phase_plus[i]+ M_PI/2.;
		}
	}
	else if(local_method.find("IMRPhenomPv2")!=std::string::npos)
	{
		T *phase_plus_temp = new T[length];
		T *phase_cross_temp = new T[length];
		if(local_method == "ppE_IMRPhenomPv2_Inspiral")
		{
			ppE_IMRPhenomPv2_Inspiral<T> ppemodel;
			status = ppemodel.construct_phase(frequencies, length, phase_plus_temp, phase_cross_temp, &params);

		}
		else if(local_method == "ppE_IMRPhenomPv2_IMR")
		{
			ppE_IMRPhenomPv2_IMR<T> ppemodel;
			status = ppemodel.construct_phase(frequencies, length, phase_plus_temp, phase_cross_temp, &params);

		}
		else if(local_method == "gIMRPhenomPv2")
		{
			gIMRPhenomPv2<T> gmodel;
			status = gmodel.construct_phase(frequencies, length, phase_plus_temp, phase_cross_temp, &params);

		}
		else if(local_method == "IMRPhenomPv2_NRT")
		  {
		    IMRPhenomPv2_NRT<T> modelNRT;
		    status = modelNRT.construct_phase(frequencies, length, phase_plus_temp, phase_cross_temp, &params);

		  }
		else{
			IMRPhenomPv2<T> model;
			status = model.construct_phase(frequencies, length, phase_plus_temp, phase_cross_temp, &params);
		}
		unwrap_array(phase_plus_temp, phase_plus, length);
		unwrap_array(phase_cross_temp, phase_cross, length);
		delete [] phase_plus_temp;
		delete [] phase_cross_temp;
	}
	cleanup_source_parameters(&params,generation_method);




	//params = params.populate_source_parameters(mass1, mass2, Luminosity_Distance, spin1, spin2, phi_c,t_c);
	//params = params.populate_source_parameters(parameters);
	//params.populate_source_parameters(parameters);
	//params.phi = parameters->phi;
	//params.theta = parameters->theta;
	//params.incl_angle = parameters->incl_angle;
	//params.f_ref = parameters->f_ref;
	//params.phiRef = parameters->phiRef;
	//params.cosmology = parameters->cosmology;
	//params.shift_time = parameters->shift_time;
	//params.shift_phase = parameters->shift_phase;
	//params.NSflag1 = parameters->NSflag1;
	//params.NSflag2 = parameters->NSflag2;
	//params.sky_average = parameters->sky_average;
	//params.dep_postmerger = parameters->dep_postmerger;

	//if(generation_method == "IMRPhenomD")
	//{
	//	IMRPhenomD<T> modeld;
	//	status = modeld.construct_phase(frequencies, length, phase_plus, &params);
	//	for(int i = 0 ; i<length; i++){
	//		phase_cross[i] = phase_plus[i]+ M_PI/2.;
	//		//phase_plus[i]*= (T)(-1.);
	//		//phase_cross[i] = phase_plus[i]- M_PI/2.;
	//	}
	//}
	//else if(generation_method == "ppE_IMRPhenomD_Inspiral")
	//{
	//	params.betappe = parameters->betappe;
	//	params.bppe = parameters->bppe;
	//	params.Nmod = parameters->Nmod;
	//	ppE_IMRPhenomD_Inspiral<T> ppemodeld;
	//	status = ppemodeld.construct_phase(frequencies, length, phase_plus, &params);
	//	for(int i = 0 ; i<length; i++){
	//		//phase_plus[i]*= (T)(-1.);
	//		phase_cross[i] = phase_plus[i]+ M_PI/2.;
	//	}
	//}
	//else if(generation_method == "_dCS_IMRPhenomD")
	//{
	//	bool local_spline = false;
	//	params.betappe = parameters->betappe;
	//	params.Nmod = 1;
	//	int tempbppe[params.Nmod] = {-1};
	//	params.bppe = tempbppe;
	//	T temp[params.Nmod] ;
	//	for( int i = 0; i < params.Nmod; i++)
	//		temp[i] = params.betappe[i];
	//	dCS_IMRPhenomD<T> ppemodeld;
	//	status = ppemodeld.construct_phase(frequencies, length, phase_plus, &params);
	//	for(int i = 0 ; i<length; i++){
	//		//phase_plus[i]*= (T)(-1.);
	//		phase_cross[i] = phase_plus[i]+ M_PI/2.;
	//	}
	//
	//	for( int i = 0; i < params.Nmod; i++)
	//		parameters->betappe[i] = temp[i];
	//}
	//else if(generation_method == "EdGB_IMRPhenomD")
	//{
	//	params.betappe = parameters->betappe;
	//	params.Nmod = 1;
	//	int tempbppe[params.Nmod] = {-7};
	//	params.bppe = tempbppe;
	//	T temp[params.Nmod] ;
	//	for( int i = 0; i < params.Nmod; i++)
	//		temp[i] = params.betappe[i];
	//	EdGB_IMRPhenomD<T> ppemodeld;
	//	status = ppemodeld.construct_phase(frequencies, length, phase_plus, &params);
	//	for(int i = 0 ; i<length; i++){
	//		//phase_plus[i]*= (T)(-1.);
	//		phase_cross[i] = phase_plus[i]+ M_PI/2.;
	//	}

	//	for( int i = 0; i < params.Nmod; i++)
	//		parameters->betappe[i] = temp[i];
	//}
	//else if(generation_method == "ppE_IMRPhenomD_IMR")
	//{
	//	params.betappe = parameters->betappe;
	//	params.bppe = parameters->bppe;
	//	params.Nmod = parameters->Nmod;
	//	ppE_IMRPhenomD_IMR<T> ppemodeld;
	//	status = ppemodeld.construct_phase(frequencies, length, phase_plus, &params);
	//	for(int i = 0 ; i<length; i++){
	//		//phase_plus[i]*= (T)(-1.);
	//		phase_cross[i] = phase_plus[i]+ M_PI/2.;
	//	}
	//	//for(int i = 0 ; i<length; i++){
	//	//	phase_plus[i]*= (T)(-1.);
	//	//	phase_cross[i] = phase_plus[i]+ M_PI/2.;
	//	//}
	//}
	//else if(generation_method == "gIMRPhenomD")
	//{
	//	params.delta_phi = parameters->delta_phi;
	//	params.delta_sigma = parameters->delta_sigma;
	//	params.delta_beta = parameters->delta_beta;
	//	params.delta_alpha = parameters->delta_alpha;
	//	params.phii = parameters->phii;
	//	params.sigmai = parameters->sigmai;
	//	params.betai = parameters->betai;
	//	params.alphai = parameters->alphai;
	//	params.Nmod_phi = parameters->Nmod_phi;
	//	params.Nmod_sigma = parameters->Nmod_sigma;
	//	params.Nmod_beta = parameters->Nmod_beta;
	//	params.Nmod_alpha = parameters->Nmod_alpha;
	//	gIMRPhenomD<T> gmodeld;
	//	status = gmodeld.construct_phase(frequencies, length, phase_plus, &params);
	//	for(int i = 0 ; i<length; i++){
	//		//phase_plus[i]*= (T)(-1.);
	//		phase_cross[i] = phase_plus[i]+ M_PI/2.;
	//	}
	//	//for(int i = 0 ; i<length; i++){
	//	//	phase_plus[i]*= (T)(-1.);
	//	//	phase_cross[i] = phase_plus[i]+ M_PI/2.;
	//	//}
	//}
	//else if(generation_method == "gIMRPhenomPv2")
	//{
	//	gIMRPhenomPv2<T> gmodelp;
	//	if((parameters->chip +1)>DOUBLE_COMP_THRESH){
	//		params.chip = parameters->chip;
	//		params.spin1z = parameters->spin1[2];
	//		params.spin2z = parameters->spin2[2];
	//		params.phip = parameters->phip;
	//		gmodelp.PhenomPv2_Param_Transform_reduced(&params);
	//	}
	//	else {
	//		gmodelp.PhenomPv2_Param_Transform(&params);
	//	}
	//	params.delta_phi = parameters->delta_phi;
	//	params.delta_sigma = parameters->delta_sigma;
	//	params.delta_beta = parameters->delta_beta;
	//	params.delta_alpha = parameters->delta_alpha;
	//	params.phii = parameters->phii;
	//	params.sigmai = parameters->sigmai;
	//	params.betai = parameters->betai;
	//	params.alphai = parameters->alphai;
	//	params.Nmod_phi = parameters->Nmod_phi;
	//	params.Nmod_sigma = parameters->Nmod_sigma;
	//	params.Nmod_beta = parameters->Nmod_beta;
	//	params.Nmod_alpha = parameters->Nmod_alpha;
	//	T *phase_plus_temp = new T[length];
	//	T *phase_cross_temp = new T[length];
	//	status = gmodelp.construct_phase(frequencies, length, phase_plus_temp, phase_cross_temp, &params);
	//	unwrap_array(phase_plus_temp, phase_plus, length);
	//	unwrap_array(phase_cross_temp, phase_cross, length);
	//	delete [] phase_plus_temp;
	//	delete [] phase_cross_temp;
	//}
	//else if(generation_method == "IMRPhenomPv2")
	//{
	//	//std::complex<T> ci = std::complex<T>(cos(params.incl_angle),0);

	//	IMRPhenomPv2<T> modeld;
	//	//Calculate Waveform
	//	if((parameters->chip +1)>DOUBLE_COMP_THRESH){
	//		params.chip = parameters->chip;
	//		params.spin1z = parameters->spin1[2];
	//		params.spin2z = parameters->spin2[2];
	//		params.phip = parameters->phip;
	//		modeld.PhenomPv2_Param_Transform_reduced(&params);
	//	}
	//	else {
	//		modeld.PhenomPv2_Param_Transform(&params);
	//	}
	//	T *phase_plus_temp = new T[length];
	//	T *phase_cross_temp = new T[length];
	//	status = modeld.construct_phase(frequencies, length, phase_plus_temp, phase_cross_temp, &params);
	//	unwrap_array(phase_plus_temp, phase_plus, length);
	//	unwrap_array(phase_cross_temp, phase_cross, length);
	//	delete [] phase_plus_temp;
	//	delete [] phase_cross_temp;
	//	//for(int i = 0 ; i<length; i++){
	//	//	phase_plus[i]*= (T)(-1.);
	//	//	phase_cross[i]*= (T)(-1.);
	//	//}
	//}
	//else if(generation_method == "ppE_IMRPhenomPv2_Inspiral")
	//{
	//	//std::complex<T> ci = std::complex<T>(cos(params.incl_angle),0);
	//	params.betappe = parameters->betappe;
	//	params.bppe = parameters->bppe;
	//	params.Nmod = parameters->Nmod;

	//	ppE_IMRPhenomPv2_Inspiral<T> modeld;
	//	//Calculate Waveform
	//	if((parameters->chip +1)>DOUBLE_COMP_THRESH){
	//		params.chip = parameters->chip;
	//		params.spin1z = parameters->spin1[2];
	//		params.spin2z = parameters->spin2[2];
	//		params.phip = parameters->phip;
	//		modeld.PhenomPv2_Param_Transform_reduced(&params);
	//	}
	//	else {
	//		modeld.PhenomPv2_Param_Transform(&params);
	//	}
	//	T *phase_plus_temp = new T[length];
	//	T *phase_cross_temp = new T[length];
	//	status = modeld.construct_phase(frequencies, length, phase_plus_temp, phase_cross_temp, &params);
	//	unwrap_array(phase_plus_temp, phase_plus, length);
	//	unwrap_array(phase_cross_temp, phase_cross, length);
	//	delete [] phase_plus_temp;
	//	delete [] phase_cross_temp;
	//	//for(int i = 0 ; i<length; i++){
	//	//	phase_plus[i]*= (T)(-1.);
	//	//	phase_cross[i]*= (T)(-1.);
	//	//}
	//}
	//else if(generation_method == "ppE_IMRPhenomPv2_IMR")
	//{
	//	//std::complex<T> ci = std::complex<T>(cos(params.incl_angle),0);
	//	params.betappe = parameters->betappe;
	//	params.bppe = parameters->bppe;
	//	params.Nmod = parameters->Nmod;

	//	ppE_IMRPhenomPv2_IMR<T> modeld;
	//	//Calculate Waveform
	//	if((parameters->chip +1)>DOUBLE_COMP_THRESH){
	//		params.chip = parameters->chip;
	//		params.spin1z = parameters->spin1[2];
	//		params.spin2z = parameters->spin2[2];
	//		params.phip = parameters->phip;
	//		modeld.PhenomPv2_Param_Transform_reduced(&params);
	//	}
	//	else {
	//		modeld.PhenomPv2_Param_Transform(&params);
	//	}
	//	T *phase_plus_temp = new T[length];
	//	T *phase_cross_temp = new T[length];
	//	status = modeld.construct_phase(frequencies, length, phase_plus_temp, phase_cross_temp, &params);
	//	unwrap_array(phase_plus_temp, phase_plus, length);
	//	unwrap_array(phase_cross_temp, phase_cross, length);
	//	delete [] phase_plus_temp;
	//	delete [] phase_cross_temp;
	//	//for(int i = 0 ; i<length; i++){
	//	//	phase_plus[i]*= (T)(-1.);
	//	//	phase_cross[i]*= (T)(-1.);
	//	//}
	//}
	//else if(generation_method == "dCS_IMRPhenomPv2")
	//{
	//	//########################################
	//	//convert betappe for dCS (alpha**2) to the full betappe
	//	//dCS only supports one modification
	//	dCS_IMRPhenomD<T> dcs_phenomd;
	//	params.Nmod = 1;
	//	params.bppe = new int[1];
	//	params.bppe[0] = -1;
	//	params.betappe = new T[1];
	//	params.betappe[0] = parameters->betappe[0];
	//	params.betappe[0] = dcs_phenomd.dCS_phase_mod(&params);
	//	//########################################

	//	ppE_IMRPhenomPv2_Inspiral<T> model;
	//	//Initialize Pv2 specific params

	//	//########################################
	//	if((parameters->chip +1)>DOUBLE_COMP_THRESH){
	//		params.chip = parameters->chip;
	//		params.spin1z = parameters->spin1[2];
	//		params.spin2z = parameters->spin2[2];
	//		params.phip = parameters->phip;
	//		model.PhenomPv2_Param_Transform_reduced(&params);
	//	}
	//	else {
	//		model.PhenomPv2_Param_Transform(&params);
	//	}
	//	T *phase_plus_temp = new T[length];
	//	T *phase_cross_temp = new T[length];
	//	status = model.construct_phase(frequencies, length, phase_plus_temp, phase_cross_temp, &params);
	//	unwrap_array(phase_plus_temp, phase_plus, length);
	//	unwrap_array(phase_cross_temp, phase_cross, length);
	//	delete [] phase_plus_temp;
	//	delete [] phase_cross_temp;
	//	delete [] params.bppe;
	//	delete [] params.betappe;
	//}
	//else if(generation_method == "EdGB_IMRPhenomPv2")
	//{
	//	//########################################
	//	//convert betappe for dCS (alpha**2) to the full betappe
	//	//dCS only supports one modification
	//	EdGB_IMRPhenomD<T> EdGB_phenomd;
	//	params.Nmod = 1;
	//	params.bppe = new int[1];
	//	params.bppe[0] = -7;
	//	params.betappe = new T[1];
	//	params.betappe[0] = parameters->betappe[0];
	//	params.betappe[0] = EdGB_phenomd.EdGB_phase_mod(&params);
	//	//########################################

	//	ppE_IMRPhenomPv2_Inspiral<T> model;
	//	//Initialize Pv2 specific params

	//	//########################################
	//	if((parameters->chip +1)>DOUBLE_COMP_THRESH){
	//		params.chip = parameters->chip;
	//		params.spin1z = parameters->spin1[2];
	//		params.spin2z = parameters->spin2[2];
	//		params.phip = parameters->phip;
	//		model.PhenomPv2_Param_Transform_reduced(&params);
	//	}
	//	else {
	//		model.PhenomPv2_Param_Transform(&params);
	//	}
	//	T *phase_plus_temp = new T[length];
	//	T *phase_cross_temp = new T[length];
	//	status = model.construct_phase(frequencies, length, phase_plus_temp, phase_cross_temp, &params);
	//	unwrap_array(phase_plus_temp, phase_plus, length);
	//	unwrap_array(phase_cross_temp, phase_cross, length);
	//	delete [] phase_plus_temp;
	//	delete [] phase_cross_temp;
	//	delete [] params.bppe;
	//	delete [] params.betappe;
	//}

	return status ;
}



//template int fourier_waveform<double>(double *, int, std::complex<double> *,std::complex<double> *, std::string, gen_params_base<double> *);
//template int fourier_waveform<adouble>(adouble *, int, std::complex<adouble> *,std::complex<adouble> *, std::string, gen_params_base<adouble> *);
template int fourier_waveform<double>(double *, int, waveform_polarizations<double> *wp, std::string, gen_params_base<double> *);
template int fourier_waveform<adouble>(adouble *, int, waveform_polarizations<adouble> *wp, std::string, gen_params_base<adouble> *);

template int fourier_phase<double>(double *, int, double *,double *, std::string, gen_params_base<double> *);
template int fourier_phase<adouble>(adouble *, int, adouble *,adouble *, std::string, gen_params_base<adouble> *);

template int time_waveform<double>(double *, int, waveform_polarizations<double> *wp, std::string, gen_params_base<double> *);
template int time_waveform<adouble>(adouble *, int, waveform_polarizations<adouble> *wp, std::string, gen_params_base<adouble> *);

template<class T>
std::string prep_source_parameters(source_parameters<T> *out, gen_params_base<T> *in,std::string generation_method){
	std::string local_method = generation_method;
	out->populate_source_parameters(in);
	out->phi = in->phi;
	out->theta = in->theta;
	out->incl_angle = in->incl_angle;
	out->f_ref = in->f_ref;
	out->phiRef = in->phiRef;
	out->cosmology = in->cosmology;
	out->shift_time = in->shift_time;
	out->sky_average = in->sky_average;
	out->shift_phase = in->shift_phase;
	out->NSflag1 = in->NSflag1;
	out->NSflag2 = in->NSflag2;
	out->dep_postmerger = in->dep_postmerger;
	out->include_l1 = in->include_l1;
	if(generation_method.find("Pv2")!=std::string::npos){
		IMRPhenomPv2<T> model;
		if((in->chip +1)>DOUBLE_COMP_THRESH){
			out->chip = in->chip;
			out->spin1z = in->spin1[2];
			out->spin2z = in->spin2[2];
			out->phip = in->phip;
			model.PhenomPv2_Param_Transform_reduced(out);
		}
		else {
			model.PhenomPv2_Param_Transform(out);
		}
	}
	else if(generation_method.find("Pv3")!=std::string::npos)
	{
		if (in->mass1 < in->mass2)
		{
			PhenomPrecessingSpinEnforcePrimary(&(in->mass1), &(in->mass2),
			&(in->spin1[0]), &(in->spin1[1]), &(in->spin1[2]),
			&(in->spin2[0]), &(in->spin2[1]), &(in->spin2[2]));
			out->populate_source_parameters(in);
		}

		PhenomPv3_Param_Transform(out, in);
	}
	if(generation_method.find("ppE") != std::string::npos){
		out->Nmod = in->Nmod;
		out->betappe = in->betappe;
		out->bppe = in->bppe;
	}
	if(generation_method.find("gIMR") != std::string::npos){
			out->delta_phi = in->delta_phi;
			out->delta_sigma = in->delta_sigma;
			out->delta_beta = in->delta_beta;
			out->delta_alpha = in->delta_alpha;
			out->phii = in->phii;
			out->sigmai = in->sigmai;
			out->betai = in->betai;
			out->alphai = in->alphai;
			out->Nmod_phi = in->Nmod_phi;
			out->Nmod_sigma = in->Nmod_sigma;
			out->Nmod_beta = in->Nmod_beta;
			out->Nmod_alpha = in->Nmod_alpha;
	}

	if(generation_method.find("NRT") != std::string::npos){
	  //if(in->tidal_s >=0)
	  if(in->tidal_love )
	    {
	      IMRPhenomD_NRT<T> modelNRT; 
	      modelNRT.binary_love_relation(in->tidal_s, in->tidal_love_error, out);
	      in->tidal1= out->tidal1; //copying into the gen params struct
	      in->tidal2= out->tidal2;
	      
	      in->tidal_a = (in->tidal1 - in->tidal2)/2.; 

	    }
	  if((in->tidal1 < 0 || in->tidal2<0) && in->tidal_weighted >= 0) {
	    out->tidal_weighted = in->tidal_weighted;
	  }
	  else if((in->tidal1 >= 0 && in->tidal2>=0) ) {
	    out->tidal1 = in->tidal1;
	    out->tidal2 = in->tidal2;
	    //arXiv 1402.5156, Eqs (5) (tidal_weighted) and (6) (delta_tidal_weighted)
	    out->tidal_weighted = 8./13. * ( (1. + 7.*out->eta - 31.*out->eta*out->eta)*(out->tidal1 + out->tidal2)
					     + sqrt( 1. - 4.*out->eta ) * ( 1. + 9.*out->eta - 11. * out->eta*out->eta) * (out->tidal1 - out->tidal2) ) ;

	    out->delta_tidal_weighted = 1./2. * ( sqrt( 1. - 4.*out->eta ) * ( 1. - 13272./1319. * out->eta + 8944./1319. * out->eta*out->eta) *
						  (out->tidal1 + out->tidal2) + ( 1. - 15910./1319. * out->eta + 32850./1319. * out->eta*out->eta + 3380./1319.
										  * out->eta *out->eta*out->eta)*(out->tidal1-out->tidal2));
		//debugger_print(__FILE__,__LINE__,out->tidal_weighted);
		}
	  // Dissipative tidal deformability
          if((in->diss_tidal1 < 0 || in->diss_tidal2<0) && in->diss_tidal_weighted >= 0) {
	    out->diss_tidal_weighted = in->diss_tidal_weighted;
	  }
	  else if((in->diss_tidal1 >= 0 && in->diss_tidal2>=0) ) {
	    out->diss_tidal1 = in->diss_tidal1;
	    out->diss_tidal2 = in->diss_tidal2;
	    //arXiv 2306.15633, Eqs (13b) 
	    out->diss_tidal_weighted =  
               (2.0*pow(out->eta,2) - 4.0*out->eta + 1.0)*0.5*(out->diss_tidal1 + out->diss_tidal2) 
               -
               sqrt(1.0 - 4.0*out->eta)*(1.0 - 2.0*out->eta)*0.5*(out->diss_tidal1 - out->diss_tidal2) 
                ;

	  }
	}
	if(generation_method.find("EA_IMRPhenomD") != std::string::npos){
	  out->alpha_param = in->alpha_param;
	  out->EA_region1 = in->EA_region1;
	  if(in->alpha_param){
	        out->alpha1_EA = in->alpha1_EA;
		out->alpha2_EA = in->alpha2_EA;
		out->cbarw_EA = in->cbarw_EA; 
	  }
	  else{
	    if(in->EA_region1){
	      out->ca_EA = in->ca_EA;
	      out->ctheta_EA = 3.* in->ca_EA*(1. + in->delta_ctheta_EA);
	      out->cw_EA = in->cw_EA; 
	    }
	    else{
		out->ca_EA = in->ca_EA;
		out->ctheta_EA = in->ctheta_EA;
		out->cw_EA = in->cw_EA;
	    }
	  }
		out->csigma_EA = in->csigma_EA;
		if(generation_method.find("NRT") != std::string::npos)
		  {
		    EA_IMRPhenomD_NRT<T> EAmodeldNRT;
		    EAmodeldNRT.pre_calculate_EA_factors(out);
		  }
		else
		  {
		    EA_IMRPhenomD<T> EAmodeld;
		    EAmodeld.pre_calculate_EA_factors(out);
		  }
	}
	if(check_theory_support(generation_method)){
		theory_ppE_map<T> mapping;
		assign_mapping<T>(generation_method,&mapping,in);
		out->Nmod = in->Nmod;
		out->bppe = new double[mapping.Nmod];
		out->betappe = new T[out->Nmod];
		//The input beta vector might contain theory specific
		//parameters that happen at multiple PN orders --
		//save the beta output and assign at the end
		T *temp_beta = new T[mapping.Nmod];
		for(int i = 0 ; i<out->Nmod; i++){
			out->betappe[i]=in->betappe[i];
		}
		for(int i = 0 ; i<mapping.Nmod; i++){
			out->bppe[i]=mapping.bppe[i];
			temp_beta[i]=mapping.beta_fns[i](out);

		}
		delete[] out->betappe;
		out->betappe = new T[mapping.Nmod];
		out->Nmod = mapping.Nmod;
		for(int i = 0 ; i<out->Nmod; i++){
			out->betappe[i]=temp_beta[i];
		}
		delete [] temp_beta;
		local_method = mapping.ppE_method;
		deallocate_mapping(&mapping);
	}
	if(generation_method.find("TaylorT2") != std::string::npos){
		out->x0 = in->x0;
	}
	return local_method;
}
template std::string prep_source_parameters(source_parameters<adouble> *, gen_params_base<adouble> *, std::string);
template std::string prep_source_parameters(source_parameters<double> *, gen_params_base<double> *, std::string);

template<class T>
void cleanup_source_parameters(source_parameters<T> *params,std::string generation_method)
{
	if(check_theory_support(generation_method)){
		delete [] params->betappe;params->betappe=NULL;
		delete [] params->bppe;	params->bppe=NULL;
	}

}
template void cleanup_source_parameters(source_parameters<adouble> *,std::string );
template void cleanup_source_parameters(source_parameters<double> *,std::string );

template<class T>
void waveform_polarizations<T>::allocate_memory(int length)
{
	if(this->active_polarizations[0]){
		this->hplus = new std::complex<T>[length];
	}
	if(this->active_polarizations[1]){
		this->hcross = new std::complex<T>[length];
	}
	if(this->active_polarizations[2]){
		this->hx = new std::complex<T>[length];
	}
	if(this->active_polarizations[3]){
		this->hy = new std::complex<T>[length];
	}
	if(this->active_polarizations[4]){
		this->hb = new std::complex<T>[length];
	}
	if(this->active_polarizations[5]){
		this->hl = new std::complex<T>[length];
	}
	return;
}
template void waveform_polarizations<double>::allocate_memory(int length);
template void waveform_polarizations<adouble>::allocate_memory(int length);

template<class T>
void waveform_polarizations<T>::deallocate_memory()
{
	if(this->hplus){
		delete [] this->hplus;
		this->hplus = NULL;
	}
	if(this->hcross){
		delete [] this->hcross;
		this->hcross = NULL;
	}
	if(this->hx){
		delete [] this->hx;
		this->hx = NULL;
	}
	if(this->hy){
		delete [] this->hy;
		this->hy = NULL;
	}
	if(this->hb){
		delete [] this->hb;
		this->hb = NULL;
	}
	if(this->hl){
		delete [] this->hl;
		this->hl = NULL;
	}
	return;
}
template void waveform_polarizations<double>::deallocate_memory();
template void waveform_polarizations<adouble>::deallocate_memory();

bool check_extra_polarizations(std::string generation_method)
{
	if(generation_method == "polarization_test_IMRPhenomD"){
		return true;
	}
	if(generation_method.find("EA_IMRPhenomD") != std::string::npos){
		return true;
	}
	return false;
}

template<class T>
void assign_polarizations(std::string generation_method, waveform_polarizations<T> *wp)
{
	if(generation_method == "polarization_test_IMRPhenomD"){
		wp->active_polarizations[0]=true;
		wp->active_polarizations[1]=true;
		wp->active_polarizations[2]=true;
		wp->active_polarizations[3]=true;
		wp->active_polarizations[4]=true;
		wp->active_polarizations[5]=true;
	}
	else if(generation_method.find("EA_IMRPhenomD") != std::string::npos){
		wp->active_polarizations[0]=true;
		wp->active_polarizations[1]=true;
		wp->active_polarizations[2]=true;
		wp->active_polarizations[3]=true;
		wp->active_polarizations[4]=true;
		wp->active_polarizations[5]=true;
	}
	else{
		wp->active_polarizations[0]=true;
		wp->active_polarizations[1]=true;
		wp->active_polarizations[2]=false;
		wp->active_polarizations[3]=false;
		wp->active_polarizations[4]=false;
		wp->active_polarizations[5]=false;
	}
	return ;
}
template void assign_polarizations<double>(std::string generation_method, waveform_polarizations<double> *wp);
template void assign_polarizations<adouble>(std::string generation_method, waveform_polarizations<adouble> *wp);
