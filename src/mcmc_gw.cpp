#include "mcmc_gw.h"
#include "waveform_generator.h"
#include "util.h"
#include "io_util.h"
#include "detector_util.h"
#include "ppE_utilities.h"
#include "waveform_util.h"
#include "ortho_basis.h"
#include "fisher.h"
#include "mcmc_sampler.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <fftw3.h>
#include <algorithm>
#include <iostream>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <mutex>
#include <thread>
#include <condition_variable>
#include "IMRPhenomD_NRT.h" //For testing purposes only! Remove later.
#include "EA_IMRPhenomD_NRT.h" //For testing purposes only! Remove later.

/*! \file 
 * Routines for implementation in MCMC algorithms specific to GW CBC analysis
 *
 * */

/*!\brief Function to calculate the log Likelihood as defined by -1/2 (d-h|d-h) maximized over the extrinsic parameters phic and tc
 *
 * frequency array must be uniform spacing - this shouldn't be a problem when working with real data as DFT return uniform spacing
 */
double maximized_coal_log_likelihood_IMRPhenomD(double *frequencies,
				int length,
				std::complex<double> *data,
				double *noise,
				double SNR,
				double chirpmass,	/**< in solar masses*/
				double symmetric_mass_ratio, 
				double spin1,
				double spin2,
				bool NSflag,
				fftw_outline *plan
				)
{
	//Make template waveform
	gen_params params;
	params.mass1 = calculate_mass1(chirpmass, symmetric_mass_ratio);	
	params.mass2 = calculate_mass2(chirpmass, symmetric_mass_ratio);	
	params.Luminosity_Distance = 500;	
	std::complex<double> *template_strain = (std::complex<double>*)malloc(sizeof(std::complex<double>)*length);	
	params.spin1[0]= 0;
	params.spin1[1]= 0;
	params.spin1[2]= spin1;
	params.spin2[0]= 0;
	params.spin2[1]= 0;
	params.spin2[2]= spin2;
	params.phiRef=0;
	params.tc=0;
	params.NSflag1=NSflag;
	params.NSflag2=NSflag;
	//fourier_waveform(frequencies, length, template_strain, "IMRPhenomD", m1,m2,Dl,spin1vec,
	//				spin2vec);
	fourier_waveform(frequencies, length, template_strain, "IMRPhenomD", &params);
	
	//Calculate template snr and scale it to match the data snr
	//later, upgrade to non uniform spacing, cause why not
	double delta_f = frequencies[1]-frequencies[0];
	double sum = 0;
	double *integrand = (double *)malloc(sizeof(double)*length);
	for (int i =0;i< length;i++)
		integrand[i] = real(template_strain[i]*std::conj(template_strain[i]))/noise[i];
	//double integral = 4.*trapezoidal_sum_uniform(delta_f, length, integrand);
	double integral = 4.*simpsons_sum(delta_f, length, integrand);
	double normalizing_factor = SNR/sqrt(integral);
	
	//calculate the fourier transform that corresponds to maximizing over phic and tc
	//Use malloc at some point, not sure how long these arrays will be 
	//std::complex<double> g_tilde[length];
	std::complex<double> g_tilde;
	for (int i=0;i<length; i++)
	{
		g_tilde = 4.*normalizing_factor*conj(data[i]) * template_strain[i] / noise[i]; 
		plan->in[i][0] = real(g_tilde);
		plan->in[i][1] = imag(g_tilde);
	}

	double *g = (double *)malloc(sizeof(double) *length);
	
	fftw_execute(plan->p);
	
	for (int i=0;i<length; i++)
	{
		g[i] = std::abs(std::complex<double>(plan->out[i][0],plan->out[i][1])) ;
	}

	double max = *std::max_element(g, g+length)*delta_f; 
	
	free(template_strain);
	free(g);
	free(integrand);
	return -(SNR*SNR - max);
}
				
double maximized_coal_log_likelihood_IMRPhenomD(double *frequencies,
				size_t length,
				double *real_data,
				double *imag_data,
				double *noise,
				double SNR,
				double chirpmass,	/**< in solar masses*/
				double symmetric_mass_ratio, 
				double spin1,
				double spin2,
				bool NSflag)
{
	fftw_outline plan;
	allocate_FFTW_mem_forward(&plan, length);
	std::complex<double> *data = (std::complex<double> *)malloc(sizeof(std::complex<double>)*length);
	for (int i =0; i<length; i++)
		data[i] = std::complex<double>(real_data[i],imag_data[i]);
	double out =  maximized_coal_log_likelihood_IMRPhenomD(frequencies,
				length,
				data,
				noise,
				SNR,
				chirpmass,	
				symmetric_mass_ratio, 
				spin1,
				spin2,
				NSflag,
				&plan);
	deallocate_FFTW_mem(&plan);
	free(data);
	return out;
}
double maximized_coal_log_likelihood_IMRPhenomD(double *frequencies,
				size_t length,
				double *real_data,
				double *imag_data,
				double *noise,
				double SNR,
				double chirpmass,	/**< in solar masses*/
				double symmetric_mass_ratio, 
				double spin1,
				double spin2,
				bool NSflag,
				fftw_outline *plan)
{
	std::complex<double> *data = (std::complex<double> *)malloc(sizeof(std::complex<double>)*length);
	for (int i =0; i<length; i++)
		data[i] = std::complex<double>(real_data[i],imag_data[i]);

	double out =  maximized_coal_log_likelihood_IMRPhenomD(frequencies,
				length,
				data,
				noise,
				SNR,
				chirpmass,	
				symmetric_mass_ratio, 
				spin1,
				spin2,
				NSflag,
				plan);
	free(data);
	return out;
}
				
double maximized_coal_log_likelihood_IMRPhenomD_Full_Param(double *frequencies,
				int length,
				std::complex<double> *data,
				double *noise,
				double chirpmass,	/**< in solar masses*/
				double symmetric_mass_ratio, 
				double spin1,
				double spin2,
				double Luminosity_Distance,
				double theta,
				double phi,
				double iota,
				bool NSflag,
				fftw_outline *plan)
{

	//Make template waveform
	gen_params params;
	params.mass1 = calculate_mass1(chirpmass, symmetric_mass_ratio);	
	params.mass2 = calculate_mass2(chirpmass, symmetric_mass_ratio);	
	params.Luminosity_Distance = Luminosity_Distance;	
	std::complex<double> *template_strain = (std::complex<double> *)malloc(sizeof(std::complex<double>)*length);	
	params.spin1[0]= 0;
	params.spin1[1]= 0;
	params.spin1[2]= spin1;
	params.spin2[0]= 0;
	params.spin2[1]= 0;
	params.spin2[2]= spin2;
	params.phiRef=0.;
	params.tc=0.0;
	params.NSflag1=NSflag;
	params.NSflag2=NSflag;

	fourier_waveform(frequencies, length, template_strain, "IMRPhenomD", &params);
	std::complex<double> q = Q(theta,phi,iota);
	for (int i = 0;i<length;i++)
		template_strain[i] = template_strain[i]*q;
	

	////Calculate template snr and scale it to match the data snr
	////later, upgrade to non uniform spacing, cause why not
	//double delta_f = frequencies[1]-frequencies[0];
	//double sum = 0;
	//double *integrand = (double *)malloc(sizeof(double)*length);
	//for (int i =0;i< length;i++)
	//	integrand[i] = real(template_strain[i]*std::conj(template_strain[i]))/noise[i];
	////double integral = 4.*trapezoidal_sum_uniform(delta_f, length, integrand);
	//double integral = 4.*simpsons_sum(delta_f, length, integrand);
	//double HH = integral;

	//sum = 0;
	//double integral2;
	//for (int i =0;i< length;i++)
	//	integrand[i] = real(data[i]*std::conj(data[i]))/noise[i];
	////integral2 = 4*trapezoidal_sum_uniform(delta_f, length, integrand2);
	//integral2 = 4*simpsons_sum(delta_f, length, integrand);
	//double DD = integral2;

	////###################################################################
	////testing
	////###################################################################
	//sum = 0;
	//double integral3;
	//for (int i =0;i< length;i++)
	//	integrand[i] = real(template_strain[i]*std::conj(data[i]))/noise[i];
	////integral2 = 4*trapezoidal_sum_uniform(delta_f, length, integrand2);
	//integral3 = 4*simpsons_sum(delta_f, length, integrand);
	//double HD = integral3;
	////###################################################################



	////double normalizing_factor = SNR/sqrt(integral);
	//
	////calculate the fourier transform that corresponds to maximizing over phic and tc
	////Use malloc at some point, not sure how long these arrays will be 
	//std::complex<double> g_tilde;
	//for (int i=0;i<length; i++)
	//{
	//	g_tilde = 4.*conj(data[i]) * template_strain[i] / noise[i]; 
	//	plan->in[i][0] = real(g_tilde);
	//	plan->in[i][1] = imag(g_tilde);
	//}

	//double *g = (double *)malloc(sizeof(double)*length);
	//
	//fftw_execute(plan->p);
	//
	//for (int i=0;i<length; i++)
	//{
	//	g[i] = std::abs(std::complex<double>(plan->out[i][0],plan->out[i][1])) ;
	//}

	//double max = *std::max_element(g, g+length)*delta_f; 

	//free(integrand);
	//free(g);

	double out =  maximized_Log_Likelihood_aligned_spin_internal(data,
				noise,
				frequencies,
				template_strain,
				length,
				plan
				);
	free(template_strain);
	return out;
	//return -0.5*(HH- 2*HD);//TESTING
	//return -0.5*(HH- 2*max);
	//return -0.5*(DD+HH- 2*max);
}
double maximized_coal_log_likelihood_IMRPhenomD_Full_Param(double *frequencies,
				size_t length,
				double *real_data,
				double *imag_data,
				double *noise,
				double chirpmass,	/**< in solar masses*/
				double symmetric_mass_ratio, 
				double spin1,
				double spin2,
				double Luminosity_Distance,
				double theta,
				double phi,
				double iota,
				bool NSflag)
{
	fftw_outline plan;
	allocate_FFTW_mem_forward(&plan,length);
	std::complex<double> *data = (std::complex<double> *) malloc(sizeof(std::complex<double> ) *length);
	for (int i =0; i<length; i++)
		data[i] = std::complex<double>(real_data[i],imag_data[i]);
	double out =  maximized_coal_log_likelihood_IMRPhenomD_Full_Param(frequencies,
				length,
				data,
				noise,
				chirpmass,	
				symmetric_mass_ratio, 
				spin1,
				spin2,
				Luminosity_Distance,
				theta,
				phi,
				iota,
				NSflag,
				&plan);
	deallocate_FFTW_mem(&plan);

	free(data);

	return out;
}
double maximized_coal_log_likelihood_IMRPhenomD_Full_Param(double *frequencies,
				size_t length,
				double *real_data,
				double *imag_data,
				double *noise,
				double chirpmass,	/**< in solar masses*/
				double symmetric_mass_ratio, 
				double spin1,
				double spin2,
				double Luminosity_Distance,
				double theta,
				double phi,
				double iota,
				bool NSflag,
				fftw_outline *plan)
{
	std::complex<double> *data = (std::complex<double> *) malloc(sizeof(std::complex<double> ) *length);
	for (int i =0; i<length; i++)
		data[i] = std::complex<double>(real_data[i],imag_data[i]);
	double out =  maximized_coal_log_likelihood_IMRPhenomD_Full_Param(frequencies,
				length,
				data,
				noise,
				chirpmass,	
				symmetric_mass_ratio, 
				spin1,
				spin2,
				Luminosity_Distance,
				theta,
				phi,
				iota,
				NSflag,
				plan);

	free(data);

	return out;
}



/*! \brief routine to maximize over all extrinsic quantities and return the log likelihood
 *
 * IMRPhenomD -- maximizes over DL, phic, tc, \iota, \phi, \theta
 * IMRPhenomP -- maximizes over DL, phic,tc, \psi, \phi , \theta
 */
double maximized_Log_Likelihood(std::complex<double> *data, 
				double *psd,
				double *frequencies,
				size_t length,
				//gen_params *params,
				gen_params_base<double> *params,
				std::string detector,
				std::string generation_method,
				fftw_outline *plan
				)
{
	double ll = 0;
	if(	generation_method.find("IMRPhenomD")!=std::string::npos){
		std::complex<double> *response =
			(std::complex<double> *) malloc(sizeof(std::complex<double>) * length);
		fourier_detector_response_horizon(frequencies, length, response, detector, generation_method, params);
		ll = maximized_Log_Likelihood_aligned_spin_internal(data, psd, frequencies,response, length, plan);
		
		free(response);
		}
	if(	generation_method.find("IMRPhenomPv2")!=std::string::npos){
				
		//fourier_waveform(frequencies,length,hp,hc,generation_method,params);
		waveform_polarizations<double> wp;
		assign_polarizations(generation_method, &wp);
		wp.allocate_memory(length);
		fourier_waveform(frequencies,length,&wp,generation_method,params);
		ll = maximized_Log_Likelihood_unaligned_spin_internal(data,psd,frequencies,wp.hplus,wp.hcross, length, plan);
		wp.deallocate_memory();

	}
	return ll;
}

					
double maximized_coal_Log_Likelihood(double *data_real, 
				double *data_imag,
				double *psd,
				double *frequencies,
				size_t length,
				double *template_real,
				double *template_imag,
				fftw_outline *plan
				)
{
	
	std::complex<double> *data = 
			(std::complex<double> *) malloc(sizeof(std::complex<double>)*length);
	std::complex<double> *h = 
			(std::complex<double> *) malloc(sizeof(std::complex<double>)*length);
	for(int i =0; i<length; i ++){
		data[i] = std::complex<double>(data_real[i],data_imag[i]);
		h[i] = std::complex<double>(template_real[i],template_imag[i]);
	}
	double tc,phic,ll=0;
	ll = maximized_coal_Log_Likelihood_internal(data,
				psd,
				frequencies,
				h,
				length,
				plan,
				&tc,
				&phic);
	
	free(data);
	free(h);
	return ll;
}
double maximized_Log_Likelihood(double *data_real, 
				double *data_imag,
				double *psd,
				double *frequencies,
				size_t length,
				gen_params_base<double> *params,
				std::string detector,
				std::string generation_method,
				fftw_outline *plan
				)
{
	
	std::complex<double> *data = 
			(std::complex<double> *) malloc(sizeof(std::complex<double>)*length);
	for(int i =0; i<length; i ++){
		data[i] = std::complex<double>(data_real[i],data_imag[i]);
	}
	
	double ll = maximized_Log_Likelihood(data, 
				psd,
				frequencies,
				length,
				params,
				detector,
				generation_method,
				plan);
	free(data);
	return ll;
}

/*! \brief Function to maximize only over coalescence variables tc and phic, returns the maximum values used
 *
 */
double maximized_coal_Log_Likelihood(std::complex<double> *data, 
				double *psd,
				double *frequencies,
				size_t length,
				gen_params_base<double> *params,
				std::string detector,
				std::string generation_method,
				fftw_outline *plan,
				double *tc,
				double *phiRef
				)
{
	std::complex<double> *response =
		(std::complex<double> *) malloc(sizeof(std::complex<double>) * length);
	fourier_detector_response(frequencies, length, response, detector, generation_method, params);
	double ll = maximized_coal_Log_Likelihood_internal(data, psd, frequencies,
				response, length, plan, tc, phiRef);
	free(response);
	return ll;
}

double maximized_coal_Log_Likelihood_internal(std::complex<double> *data,
				double *psd,
				double *frequencies,
				std::complex<double> *detector_response,
				size_t length,
				fftw_outline *plan,
				double *tc,
				double *phiRef
				)
{

	//Calculate template snr and scale it to match the data snr
	//later, upgrade to non uniform spacing, cause why not
	double delta_f = frequencies[length/2]-frequencies[length/2-1];
	double sum = 0.;
	double *integrand = (double *)malloc(sizeof(double)*length);
	for (int i =0;i< length;i++)
		integrand[i] = real(detector_response[i]*std::conj(detector_response[i]))/psd[i];
	//double integral = 4.*trapezoidal_sum_uniform(delta_f, length, integrand);
	double integral = 4.*simpsons_sum(delta_f, length, integrand);
	double HH = integral;

	
	//calculate the fourier transform that corresponds to maximizing over phic and tc
	//Use malloc at some point, not sure how long these arrays will be 
	std::complex<double> g_tilde;

	fftw_complex *in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * length);	
	fftw_complex *out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * length);	
	for (int i=0;i<length; i++)
	{
		g_tilde = 4.*conj(data[i]) * detector_response[i] / psd[i]; 
		//g_tilde = conj(data[i]) * detector_response[i] / psd[i]; 
		//plan->in[i][0] = real(g_tilde);
		//plan->in[i][1] = imag(g_tilde);
		in[i][0] = real(g_tilde);
		in[i][1] = imag(g_tilde);
	}

	double *g = (double *)malloc(sizeof(double)*length);
	std::complex<double> *gc = (std::complex<double> *)malloc(sizeof(std::complex<double>)*length);
	
	//fftw_execute(plan->p);
	fftw_execute_dft(plan->p, in, out);
	
	for (int i=0;i<length; i++)
	{
		g[i] = std::abs(std::complex<double>(out[i][0],out[i][1])) ;
		gc[i] = std::complex<double>(out[i][0],out[i][1]) ;
		//g[i] = plan->out[i][0]*plan->out[i][0]+plan->out[i][1]*plan->out[i][1] ;
		//g[i] = out[i][0]*out[i][0]+out[i][1]*out[i][1] ;
	}
	//std::vector<int>::iterator max;
	double *max;
	max = std::max_element(g, g+length); 
	double max_val = *max*delta_f;
	//double max_val = *max;
	int max_index = std::distance(g,max);
	double Tau = 1./delta_f;
	*tc = (double)(max_index)/length * ( Tau );
	*phiRef = std::arg(gc[max_index]);
	//double max = *std::max_element(g, g+length)*delta_f; 

	free(integrand);
	free(g);
	free(gc);
	fftw_free(in);
	fftw_free(out);

	return -0.5*(HH- 2*max_val);
	//return .5*(max*max)/HH;
	//return .5*(max)/HH;
}

/*! \brief Unmarginalized log of the likelihood
 *
 */
double Log_Likelihood(std::complex<double> *data, 
				double *psd,
				double *frequencies,
				size_t length,
				gen_params_base<double> *params,
				std::string detector,
				std::string generation_method
				)
{
	double ll = 0;

	std::complex<double> *detect_response =
			(std::complex<double> *) malloc(sizeof(std::complex<double>) * length);
	fourier_detector_response(frequencies,length,detect_response,detector, generation_method,params);
	ll = Log_Likelihood_internal(data,psd,frequencies,(double*)NULL,detect_response, length, false,"SIMPSONS");

	//if(ll>0){
	//}
	free(detect_response);
	return ll;
}
/*! \brief Maximized match over coalescence variables - returns log likelihood NOT NORMALIZED for aligned spins
 *
 * Note: this function is not properly normalized for an absolute comparison. This is made for MCMC sampling, so to minimize time, constant terms like (Data|Data), which would cancel in the Metropolis-Hasting ratio, are left out for efficiency
 */
double maximized_Log_Likelihood_aligned_spin_internal(std::complex<double> *data,
				double *psd,
				double *frequencies,
				std::complex<double> *detector_response,
				size_t length,
				fftw_outline *plan
				)
{
	//Calculate template snr and scale it to match the data snr
	//later, upgrade to non uniform spacing, cause why not
	double delta_f = frequencies[1]-frequencies[0];
	double sum = 0.;
	double *integrand = (double *)malloc(sizeof(double)*length);
	for (int i =0;i< length;i++)
		integrand[i] = real(detector_response[i]*std::conj(detector_response[i]))/psd[i];
	//double integral = 4.*trapezoidal_sum_uniform(delta_f, length, integrand);
	double integral = 4.*simpsons_sum(delta_f, length, integrand);
	double HH = integral;

	
	//calculate the fourier transform that corresponds to maximizing over phic and tc
	//Use malloc at some point, not sure how long these arrays will be 
	std::complex<double> g_tilde;

	fftw_complex *in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * length);	
	fftw_complex *out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * length);	
	for (int i=0;i<length; i++)
	{
		g_tilde = 4.*conj(data[i]) * detector_response[i] / psd[i]; 
		//g_tilde = conj(data[i]) * detector_response[i] / psd[i]; 
		//plan->in[i][0] = real(g_tilde);
		//plan->in[i][1] = imag(g_tilde);
		in[i][0] = real(g_tilde);
		in[i][1] = imag(g_tilde);
	}

	double *g = (double *)malloc(sizeof(double)*length);
	
	//fftw_execute(plan->p);
	fftw_execute_dft(plan->p, in, out);
	
	for (int i=0;i<length; i++)
	{
		//g[i] = std::abs(std::complex<double>(plan->out[i][0],plan->out[i][1])) ;
		//g[i] = plan->out[i][0]*plan->out[i][0]+plan->out[i][1]*plan->out[i][1] ;
		g[i] = out[i][0]*out[i][0]+out[i][1]*out[i][1] ;
	}

	double max = *std::max_element(g, g+length)*delta_f*delta_f; 
	//double max = *std::max_element(g, g+length)*length; 

	free(integrand);
	free(g);
	fftw_free(in);
	fftw_free(out);
	//std::cout<<.5*(max)/HH<<" "<<1./length<<" "<<delta_f*delta_f<<" "<<delta_f/length<<" "<<length<<" "<<1./(delta_f*length)<<std::endl;
	return .5*(max)/HH;
}

/*! \brief log likelihood function that maximizes over extrinsic parameters tc, phic, D, and phiRef, the reference frequency - for unaligned spins 
 *
 * Ref: arXiv 1603.02444v2
 *
 * NOTE: Only works for +/x polarizations
 */
double maximized_Log_Likelihood_unaligned_spin_internal(std::complex<double> *data,
				double *psd,
				double *frequencies,
				std::complex<double> *hplus,
				std::complex<double> *hcross,
				size_t length,
				fftw_outline *plan
				)
{
	double delta_f = frequencies[1]-frequencies[0];
	double *integrand = (double *)malloc(sizeof(double)*length);
	double integral;
	//NOTE: I'm setting detector to HANFORD TEMPORARILY for testing
	//
	//THIS NEEDS TO CHANGE
	//
	//
	//
	//calculate overall template snr squared HH = <H|H> NOT NECESSARY I think
	//std::complex<double> *detector_response = 
	//	(std::complex<double> *)malloc(sizeof(std::complex<double>)*length);
	//fourier_detector_response(frequencies, length, hplus,hcross, detector_response, 0.0,0.0,"Hanford");

	//for (int i =0;i< length;i++)
	//	integrand[i] = real(detector_response[i]*std::conj(detector_response[i]))/psd[i];
	//integral = 4.*simpsons_sum(delta_f, length, integrand);
	//double HH = integral;

	//Calculate template snr for plus polarization sqrt(<H+|H+>)
	for (int i =0;i< length;i++){
		integrand[i] = real(hplus[i]*std::conj(hplus[i]))/psd[i];
	}
	integral = 4.*simpsons_sum(delta_f, length, integrand);
	double HpHproot = sqrt(integral);

	//Calculate template snr for cross polarization sqrt(<Hx|Hx>)
	for (int i =0;i< length;i++)
		integrand[i] = real(hcross[i]*std::conj(hcross[i]))/psd[i];
	integral = 4.*simpsons_sum(delta_f, length, integrand);
	double HcHcroot = sqrt(integral);
	
	//Rescale waveforms from hplus/cross to \hat{hplus/cross} 
	std::complex<double> *hpnorm = 
		(std::complex<double> *)malloc(sizeof(std::complex<double>)*length);
	std::complex<double> *hcnorm = 
		(std::complex<double> *)malloc(sizeof(std::complex<double>)*length);
	for (int i =0 ;i<length;i++)
	{
		hpnorm[i] = hplus[i]/HpHproot;
		hcnorm[i] = hcross[i]/HcHcroot;
	}

	//calculate \hat{rhoplus/cross} (just denoted rhoplus/cross) <d|hpnorm> 
	//To maximize of coalescence phase, this is an FFT (so its a vector of 
	//<d|h> at discrete tc
	double *rhoplus2 = 
		(double *)malloc(sizeof(double)*length);
	double *rhocross2 = 
		(double *)malloc(sizeof(double)*length);
	std::complex<double> *rhoplus = 
		(std::complex<double> *)malloc(sizeof(std::complex<double>)*length);
	std::complex<double> *rhocross = 
		(std::complex<double> *)malloc(sizeof(std::complex<double>)*length);
	double *gammahat = 
		(double *)malloc(sizeof(double)*length);
	std::complex<double> g_tilde;
	fftw_complex *in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * length);	
	fftw_complex *out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * length);	

	for (int i=0;i<length; i++)
	{
		g_tilde = 4.*conj(data[i]) * hpnorm[i] / psd[i]; 
		in[i][0] = real(g_tilde);
		in[i][1] = imag(g_tilde);
	}
	//fftw_execute(plan->p);
	fftw_execute_dft(plan->p, in, out);
	for (int i=0;i<length; i++)
	{
		rhoplus[i] = delta_f*(std::complex<double>(out[i][0],out[i][1]));
		//Norm of the output, squared (Re{g}^2 + Im{g}^2)
		rhoplus2[i] = delta_f*delta_f*(out[i][0]* out[i][0]+ out[i][1]* out[i][1]);
		
	}
	for (int i=0;i<length; i++)
	{
		g_tilde = 4.*conj(data[i]) * hcnorm[i] / psd[i]; 
		in[i][0] = real(g_tilde);
		in[i][1] = imag(g_tilde);
	}
	//fftw_execute(plan->p);
	fftw_execute_dft(plan->p, in, out);
	for (int i=0;i<length; i++)
	{
		rhocross[i] = delta_f*(std::complex<double>(out[i][0],out[i][1]));
		//Norm of the output, squared (Re{g}^2 + Im{g}^2)
		rhocross2[i] = delta_f*delta_f*(out[i][0]* out[i][0]+ out[i][1]* out[i][1]);
	}
	
	for (int i = 0; i <length;i++)
		gammahat[i] = real(rhoplus[i] * conj(rhocross[i]));
	
	
	for (int i =0;i< length;i++)
		integrand[i] = real(hpnorm[i]*std::conj(hcnorm[i]))/psd[i];
	integral = 4.*simpsons_sum(delta_f, length, integrand);
	double Ipc = integral;

	double *lambda = (double *)malloc(sizeof(double) * length);
	for(int i = 0; i < length; i++){
		lambda[i] =  (rhoplus2[i] + rhocross2[i] - 2*gammahat[i] * Ipc +
			sqrt( (rhoplus2[i] -rhocross2[i])*(rhoplus2[i] -rhocross2[i]) +
			4. * (Ipc*rhoplus2[i] - gammahat[i] ) * (Ipc*rhocross2[i] - gammahat[i])))
			/(1. - Ipc*Ipc);
	}
	double max = .25 * (*std::max_element(lambda, lambda+length)); 
	//double maxrhoplus2 =  (*std::max_element(rhoplus2, rhoplus2+length)); 
	//double maxrhocross2 =  (*std::max_element(rhocross2, rhocross2+length)); 
	//double maxgammahat =  (*std::max_element(gammahat, gammahat+length)); 

	free(integrand);
	free(hpnorm);
	free(hcnorm);
	free(rhoplus2);
	free(rhoplus);
	free(rhocross);
	free(rhocross2);
	free(gammahat);
	free(lambda);
	fftw_free(in);
	fftw_free(out);

	//std::cout<<max<<" "<<maxrhoplus2 <<" "<<maxrhocross2 <<" "<<maxgammahat << " "<<Ipc<<std::endl;
	//return -0.5*(HH- 2*max);
	return max;
}

/*! \brief Internal function for the unmarginalized log of the likelihood 
 *
 * .5 * ( ( h | h ) - 2 ( D | h ) )
 */
double Log_Likelihood_internal(std::complex<double> *data,
			double *psd,
			double *frequencies,
			double *weights,
			std::complex<double> *detector_response,
			int length,
			bool log10F,
			std::string integration_method
			)
{
	double delta_f = frequencies[length/2]-frequencies[length/2-1];
	double sum = 0.;
	double *integrand = (double *)malloc(sizeof(double)*length);
	for (int i =0;i< length;i++){
		integrand[i] = real(detector_response[i]*std::conj(detector_response[i]))/psd[i];
	}
	double integral = 0;
	if(integration_method=="SIMPSONS"){
		integral = 4.*simpsons_sum(delta_f, length, integrand);
	}
	else if(integration_method=="GAUSSLEG"){
		if(log10F){
			for(int i = 0 ; i<length; i++){
				integral+=weights[i]*integrand[i]*frequencies[i]*LOG10;	
			}
		}
		else{
			for(int i = 0 ; i<length; i++){
				integral+=weights[i]*integrand[i];	
			}
		}
		integral *= 4;
	}
	double HH = integral;
	integral = 0;

	for (int i =0;i< length;i++){
		integrand[i] = real(data[i]*std::conj(detector_response[i]))/psd[i];
	}
	if(integration_method=="SIMPSONS"){
		integral = 4.*simpsons_sum(delta_f, length, integrand);
	}
	else if(integration_method=="GAUSSLEG"){
		
		if(log10F){
			for(int i = 0 ; i<length; i++){
				integral+=weights[i]*integrand[i]*frequencies[i]*LOG10;	
			}
		}
		else{
			for(int i = 0 ; i<length; i++){
				integral+=weights[i]*integrand[i];	
			}
		}
		integral *= 4;
	}
	double DH = integral;

	//for (int i =0;i< length;i++)
	//	integrand[i] = real(data[i]*std::conj(data[i]))/psd[i];
	//integral = 4.*simpsons_sum(delta_f, length, integrand);
	//double DD = integral;

	free(integrand);
	//std::cout<<HH<<" "<<DH<<" "<<-0.5*(HH- 2*DH)<<std::endl;
	return -0.5*(HH- 2*DH);
	//return -2*(HH- 2*DH);
}


struct skysearch_params{
	std::complex<double> *hplus;
	std::complex<double> *hcross;
};

/*! \brief Takes in intrinsic parameters and conducts a sky search using uncorrelated MCMC
 *
 * Searches over RA, sin(DEC), PSI, INCLINATION, phiREF, tc,ln(DL)
 *
 * Only works with PhenomD
 *
 * hplus and hcross should be constructed with the value for inclination in initial_pos -- its scaled out internally
 *
 */
void SkySearch_PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW(mcmc_sampler_output *sampler_output,
	double **output,
	int dimension,
	int N_steps,
	int chain_N,
	int max_chain_N_thermo_ensemble,
	double *initial_pos,
	double *seeding_var,
	double **ensemble_initial_pos,
	double *chain_temps,
	int swp_freq,
	int t0,
	int nu,
	int max_chunk_size,
	std::string chain_distribution_scheme,
	double(*log_prior)(double *param, mcmc_data_interface *interface,void *parameters),
	int numThreads,
	bool pool,
	bool show_prog,
	int num_detectors,
	std::complex<double> **data,
	double **noise_psd,
	double **frequencies,
	int *data_length,
	double gps_time,
	std::string *detectors,
	int Nmod,
	double *bppe,
	std::complex<double> *hplus,
	std::complex<double> *hcross,
	std::string statistics_filename,
	std::string chain_filename,
	std::string likelihood_log_filename,
	std::string checkpoint_filename
	)
{
	//Create fftw plan for each detector (length of data stream may be different)
	fftw_outline *plans= (fftw_outline *)malloc(sizeof(fftw_outline)*num_detectors);
	for (int i =0;i<num_detectors;i++)
	{	
		allocate_FFTW_mem_forward(&plans[i] , data_length[i]);
	}
	void **user_parameters=NULL;
	mcmc_noise = noise_psd;	
	mcmc_frequencies = frequencies;
	mcmc_data = data;
	mcmc_data_length = data_length;
	mcmc_detectors = detectors;
	mcmc_fftw_plans = plans;
	mcmc_num_detectors = num_detectors;
	mcmc_gps_time = gps_time;
	mcmc_gmst = gps_to_GMST_radian(mcmc_gps_time);
	//std::cout<<"MCMC GMST: "<<mcmc_gmst<<std::endl;
	mcmc_Nmod = Nmod;
	mcmc_bppe = bppe;
	mcmc_log_beta = false;

	//To save time, intrinsic waveforms can be saved between detectors, if the 
	//frequencies are all the same
	//mcmc_save_waveform = true;
	//for(int i =1 ;i<mcmc_num_detectors; i++){
	//	if( mcmc_data_length[i] != mcmc_data_length[0] ||
	//		mcmc_frequencies[i][0]!= mcmc_frequencies[0][0] ||
	//		mcmc_frequencies[i][mcmc_data_length[i]-1] 
	//			!= mcmc_frequencies[0][mcmc_data_length[0]-1]){
	//		mcmc_save_waveform= false;
	//	}
	//		
	//}

	bool local_seeding ;
	if(!seeding_var){
		local_seeding = true;
	}
	else
		local_seeding = false;

	double **seeding_var_ptr = &seeding_var;
	PTMCMC_method_specific_prep("SkySearch", dimension, seeding_var_ptr, local_seeding);

	skysearch_params p;
	p.hplus = hplus;
	p.hcross = hcross;

	double cos_i_initial = (initial_pos[3]);
	double p_fac = 0.5*(1+cos_i_initial*cos_i_initial)/(exp(initial_pos[6])*MPC_SEC);
	double c_fac = (cos_i_initial)/(exp(initial_pos[6])*MPC_SEC);
	for(int i = 0 ; i<data_length[0]; i ++){
		p.hplus[i] /= p_fac;
		p.hcross[i] /= c_fac;
	}
	
	skysearch_params **P = new skysearch_params*[chain_N];
	for(int i = 0 ; i<chain_N; i ++){
		P[i]= &p;
	}
	

	PTMCMC_MH_dynamic_PT_alloc_uncorrelated(sampler_output,output, dimension, N_steps, chain_N, 
		max_chain_N_thermo_ensemble,initial_pos,seeding_var,ensemble_initial_pos, chain_temps, 
		swp_freq, t0, nu,max_chunk_size,chain_distribution_scheme,
		//log_prior,MCMC_likelihood_wrapper_SKYSEARCH, MCMC_fisher_wrapper_SKYSEARCH,(void **)P,numThreads, pool, 
		log_prior,MCMC_likelihood_wrapper_SKYSEARCH, NULL,(void **)P,numThreads, pool, 
		show_prog,statistics_filename,
		chain_filename, likelihood_log_filename,checkpoint_filename);
	
	//Deallocate fftw plans
	for (int i =0;i<num_detectors;i++)
		deallocate_FFTW_mem(&plans[i]);
	free(plans);
	if(local_seeding){ delete [] seeding_var;}
}

double MCMC_likelihood_wrapper_SKYSEARCH(double *param, mcmc_data_interface *interface,void *parameters)
{
	int dimension = interface->max_dim;
	double ll = 0 ; 
	skysearch_params *p = (skysearch_params *) parameters;
	waveform_polarizations<double> wp ;
	//assign_polarizations(generation_method,&wp);
	wp.allocate_memory(mcmc_data_length[0]);
	double p_fac = 0.5*(1+param[3]*param[3])/(exp(param[6])*MPC_SEC);
	double c_fac = (param[3])/(exp(param[6])*MPC_SEC);
	for(int i = 0 ; i<mcmc_data_length[0]; i ++){
		wp.hplus[i] = p->hplus[i]*p_fac;
		wp.hcross[i] = p->hcross[i]*c_fac;
	}
	std::complex<double> *response = new std::complex<double>[mcmc_data_length[0]];	
	double delta_t,tc;
	//double tc_ref = param[5],tc;
	double T = 1./( mcmc_frequencies[0][1]-mcmc_frequencies[0][0]);
	double tc_ref = T-param[5];
	//double llvec[3];
	for(int i = 0 ; i<mcmc_num_detectors; i++){
		delta_t = DTOA_DETECTOR(param[0],asin(param[1]),mcmc_gmst, mcmc_detectors[0],mcmc_detectors[i]);
		//tc = tc_ref + delta_t;
		tc = tc_ref - delta_t;
		tc*=2.*M_PI;
		fourier_detector_response_equatorial(mcmc_frequencies[i],mcmc_data_length[i],&wp, response,param[0],asin(param[1]),param[2],mcmc_gmst, (double *)NULL,0.,0.,0.,0.,mcmc_detectors[i]);
		for(int j = 0 ; j<mcmc_data_length[i];j++){
			//response[j]*=exp(-std::complex<double>(
			//	0,tc*(mcmc_frequencies[i][j]-20.) - param[4] ));
			response[j]*=exp(std::complex<double>(
				0,tc*(mcmc_frequencies[i][j]) - param[4] ));
		}	
		ll += Log_Likelihood_internal(mcmc_data[i], 
			mcmc_noise[i],
			mcmc_frequencies[i],
			(double*) NULL,
			response,
			(size_t) mcmc_data_length[i],
			false,
			"SIMPSONS"		
			);
		//llvec[i]=ll;
		
	}
	//debugger_print(__FILE__,__LINE__,ll);
	//std::cout<<ll<<" "<<llvec[0]<<" "<<llvec[1]-llvec[0]<<" "<<llvec[2]-llvec[1]<<std::endl;
	delete [] response;
	wp.deallocate_memory();
	return ll;
}



/*! \brief Wrapper for the MCMC_MH function, specifically for GW analysis
 *
 * Handles the details of setting up the MCMC sampler and wraps the fisher and log likelihood to conform to the format of the sampler
 *
 * *NOTE*  -- This sampler is NOT thread safe. There is global memory declared for each call to MCMC_MH_GW, so separate samplers should not be run in the same process space
 *
 * Supported parameter combinations:
 *
 * IMRPhenomD - 4 dimensions -- ln chirpmass, eta, chi1, chi2
 *
 * IMRPhenomD - 7 dimensions -- ln D_L, tc, phic, ln chirpmass, eta, chi1, chi2
 *
 * IMRPhenomD - 9 dimensions -- cos inclination, RA, DEC, ln D_L, ln chirpmass, eta, chi1, chi2, psi
 * 
 * dCS_IMRPhenomD- 8 dimensions -- cos inclination, RA, DEC, ln D_L, ln chirpmass, eta, chi1, chi2, \alpha^2 (the coupling parameter)
 *
 * dCS_IMRPhenomD_root_alpha- 8 dimensions -- cos inclination, RA, DEC, ln D_L, ln chirpmass, eta, chi1, chi2, \sqrt \alpha (in km) (the coupling parameter)
 *
 * IMRPhenomPv2 - 9 dimensions -- cos J_N, ln chirpmass, eta, |chi1|, |chi1|, theta_1, theta_2, phi_1, phi_2
 */
void PTMCMC_MH_GW(double ***output,
	int dimension,
	int N_steps,
	int chain_N,
	double *initial_pos,
	double *seeding_var,	
	double **ensemble_initial_pos,
	double *chain_temps,
	int swp_freq,
	double(*log_prior)(double *param, mcmc_data_interface *interface,void *parameters),
	int numThreads,
	bool pool,
	bool show_prog,
	int num_detectors,
	std::complex<double> **data,
	double **noise_psd,
	double **frequencies,
	int *data_length,
	double gps_time,
	std::string *detectors,
	MCMC_modification_struct *mod_struct,
	std::string generation_method,
	std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
	std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output*/
	std::string auto_corr_filename,/**< Filename to output auto correlation in some interval, if empty string, not output*/
	std::string likelihood_log_filename,
	std::string checkpoint_file/**< Filename to output data for checkpoint, if empty string, not saved*/
				)
{
	//Create fftw plan for each detector (length of data stream may be different)
	fftw_outline *plans= (fftw_outline *)malloc(sizeof(fftw_outline)*num_detectors);
	for (int i =0;i<num_detectors;i++)
	{	
		allocate_FFTW_mem_forward(&plans[i] , data_length[i]);
	}
	void **user_parameters=NULL;
	mcmc_noise = noise_psd;	
	mcmc_frequencies = frequencies;
	mcmc_data = data;
	mcmc_data_length = data_length;
	mcmc_detectors = detectors;
	mcmc_generation_method = generation_method;
	mcmc_fftw_plans = plans;
	mcmc_num_detectors = num_detectors;
	mcmc_gps_time = gps_time;
	mcmc_gmst = gps_to_GMST_radian(mcmc_gps_time);
	//mcmc_gmst = 2.456551496456678;
	//mcmc_Nmod = Nmod;
	//mcmc_bppe = bppe;
	//MCMC_modification_struct ms;
	//ms.ppE_Nmod = Nmod;
	//ms.bppe = bppe;
	//mcmc_mod_struct = &ms;
	mcmc_mod_struct=mod_struct;
	
	mcmc_log_beta = false;
	mcmc_intrinsic = false;

	//To save time, intrinsic waveforms can be saved between detectors, if the 
	//frequencies are all the same
	mcmc_save_waveform = true;
	for(int i =1 ;i<mcmc_num_detectors; i++){
		if( mcmc_data_length[i] != mcmc_data_length[0] ||
			mcmc_frequencies[i][0]!= mcmc_frequencies[0][0] ||
			mcmc_frequencies[i][mcmc_data_length[i]-1] 
				!= mcmc_frequencies[0][mcmc_data_length[0]-1]){
			mcmc_save_waveform= false;
		}
			
	}

	bool local_seeding ;
	if(!seeding_var)
		local_seeding = true;
	else
		local_seeding = false;
	double **seeding_var_ptr = &seeding_var;
	PTMCMC_method_specific_prep(generation_method, dimension, seeding_var_ptr, local_seeding);

	PTMCMC_MH(output, dimension, N_steps, chain_N, initial_pos,seeding_var,ensemble_initial_pos, chain_temps, swp_freq,
		 log_prior,MCMC_likelihood_wrapper, MCMC_fisher_wrapper,user_parameters,numThreads, pool, show_prog,statistics_filename,
		 //log_prior,MCMC_likelihood_wrapper, NULL,user_parameters,numThreads, pool, show_prog,statistics_filename,
		chain_filename, checkpoint_file);
	
	//Deallocate fftw plans
	for (int i =0;i<num_detectors;i++)
		deallocate_FFTW_mem(&plans[i]);
	free(plans);
	if(local_seeding){ delete [] seeding_var;}
}
/*! \brief Takes in an MCMC checkpoint file and continues the chain
 *
 * Obviously, the user must be sure to correctly match the dimension, number of chains, the generation_method, 
 * the prior function, the data, psds, freqs, and the detectors (number and name), and the gps_time to the 
 * previous run, otherwise the behavior of the sampler is undefined.
 *
 * numThreads and pool do not necessarily have to be the same
 */
void continue_PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW(std::string checkpoint_file_start,
	mcmc_sampler_output *sampler_output,
	double **output,
	int N_steps,
	int max_chain_N_thermo_ensemble,
	double *chain_temps,
	int swp_freq,
	int t0,
	int nu,
	int max_chunk_size,
	std::string chain_distribution_scheme,
	double(*log_prior)(double *param, mcmc_data_interface *interface,void *parameters),
	int numThreads,
	bool pool,
	bool show_prog,
	int num_detectors,
	std::complex<double> **data,
	double **noise_psd,
	double **frequencies,
	int *data_length,
	double gps_time,
	std::string *detectors,
	MCMC_modification_struct *mod_struct,
	std::string generation_method,
	std::string statistics_filename,
	std::string chain_filename,
	std::string likelihood_log_filename,
	std::string checkpoint_filename
	)
{

	std::mutex fisher_mutex;
	int dimension;
	dimension_from_checkpoint_file(checkpoint_file_start,&dimension, &dimension);
	//Create fftw plan for each detector (length of data stream may be different)
	fftw_outline *plans= (fftw_outline *)malloc(sizeof(fftw_outline)*num_detectors);
	for (int i =0;i<num_detectors;i++)
	{	
		allocate_FFTW_mem_forward(&plans[i] , data_length[i]);
	}
	mcmc_noise = noise_psd;	
	mcmc_frequencies = frequencies;
	mcmc_data = data;
	mcmc_data_length = data_length;
	mcmc_detectors = detectors;
	mcmc_generation_method = generation_method;
	mcmc_fftw_plans = plans;
	mcmc_num_detectors = num_detectors;
	mcmc_gps_time = gps_time;
	mcmc_gmst = gps_to_GMST_radian(mcmc_gps_time);
	//mcmc_Nmod = Nmod;
	//mcmc_bppe = bppe;
	mcmc_mod_struct = mod_struct;
	mcmc_log_beta = false;
	mcmc_intrinsic = false;

	//To save time, intrinsic waveforms can be saved between detectors, if the 
	//frequencies are all the same
	mcmc_save_waveform = true;
	for(int i =1 ;i<mcmc_num_detectors; i++){
		if( mcmc_data_length[i] != mcmc_data_length[0] ||
			mcmc_frequencies[i][0]!= mcmc_frequencies[0][0] ||
			mcmc_frequencies[i][mcmc_data_length[i]-1] 
				!= mcmc_frequencies[0][mcmc_data_length[0]-1]){
			mcmc_save_waveform= false;
		}
			
	}

	bool local_seeding=false;
	double *seeding_var=NULL;
	double **seeding_var_ptr = &seeding_var;
	PTMCMC_method_specific_prep(generation_method, dimension, seeding_var_ptr, local_seeding);
	//######################################################
	int T = (int)(1./(mcmc_frequencies[0][1]-mcmc_frequencies[0][0]));
	debugger_print(__FILE__,__LINE__,T);
	int burn_factor = T/4; //Take all sources to 4 seconds
	//debugger_print(__FILE__,__LINE__,burn_factor);
	std::complex<double> **burn_data = new std::complex<double>*[mcmc_num_detectors];
	double **burn_freqs = new double*[mcmc_num_detectors];
	double **burn_noise = new double*[mcmc_num_detectors];
	int *burn_lengths = new int[mcmc_num_detectors];
	fftw_outline *burn_plans= new fftw_outline[num_detectors];
	for(int j = 0; j<mcmc_num_detectors; j++){
		burn_lengths[j] = mcmc_data_length[j]/burn_factor;
		burn_data[j]= new std::complex<double>[burn_lengths[j]];
		burn_freqs[j]= new double[burn_lengths[j]];
		burn_noise[j]= new double[burn_lengths[j]];
		allocate_FFTW_mem_forward(&burn_plans[j], burn_lengths[j]);
		int ct = 0;
		for( int k = 0 ; k<mcmc_data_length[j]; k++){
			if(k%burn_factor==0 && ct<burn_lengths[j]){
				burn_data[j][ct] = mcmc_data[j][k];
				burn_freqs[j][ct] = mcmc_frequencies[j][k];
				burn_noise[j][ct] = mcmc_noise[j][k];
				ct++;
			}
		}
	}
	

	//######################################################
	//######################################################
	//Fishers sometimes need AD, but that's slow and single 
	//threaded -- use GAUSS quad with precomputed weights 
	//and abscissa 
	//double *fish_freqs = NULL;
	//double *fish_weights = NULL;
	//double **fish_psd = NULL; //Needs to be interpolated from data given
	//int fish_length = 100;
	//double flow = mcmc_frequencies[0][0];
	//double fhigh = mcmc_frequencies[0][mcmc_data_length[0]-1];
	//if(mod_struct->GAUSS_QUAD){
	//	fish_freqs = new double[fish_length];	
	//	fish_weights = new double[fish_length];	
	//	fish_psd = new double*[mcmc_num_detectors];	
	//	
	//	gauleg(log10(flow), log10(fhigh), fish_freqs, fish_weights, fish_length);
	//	for(int i = 0 ; i<fish_length; i++){
	//		fish_freqs[i]=pow(10,fish_freqs[i]);
	//	}
	//	for(int i = 0 ; i<mcmc_num_detectors; i++){
	//		fish_psd[i] = new double[fish_length];
	//		gsl_interp_accel *accel = gsl_interp_accel_alloc();
	//		gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear, mcmc_data_length[i]);
	//		gsl_spline_init(spline, mcmc_frequencies[i], mcmc_noise[i],mcmc_data_length[i]);
	//		for(int j = 0 ; j<fish_length; j++){
	//			fish_psd[i][j] = gsl_spline_eval(spline, fish_freqs[j],accel);
	//		}
	//		gsl_spline_free(spline);
	//	}
	//}
	//######################################################
	int chain_N = 0;
	int status = chain_number_from_checkpoint_file(checkpoint_file_start, &chain_N);

	debugger_print(__FILE__,__LINE__,"Number of chains: "+std::to_string(chain_N));

	MCMC_user_param **user_parameters=NULL;
	user_parameters = new MCMC_user_param*[chain_N];
	for(int i = 0 ;i<chain_N; i++){
		user_parameters[i] = new MCMC_user_param;
		
		user_parameters[i]->burn_data = burn_data;
		user_parameters[i]->burn_freqs = burn_freqs;
		user_parameters[i]->burn_noise = burn_noise;
		user_parameters[i]->burn_lengths = burn_lengths;
		user_parameters[i]->burn_plans=burn_plans;

		user_parameters[i]->mFish= &fisher_mutex;
		user_parameters[i]->GAUSS_QUAD= mod_struct->GAUSS_QUAD;
		user_parameters[i]->log10F = mod_struct->log10F;
		if(mod_struct->weights){
			user_parameters[i]->weights = mod_struct->weights;			
		}
		else{
			user_parameters[i]->weights = new double*[num_detectors];			
			for(int j = 0 ; j<num_detectors; j++){
				user_parameters[i]->weights[j]=NULL;
			}
		}

		user_parameters[i]->fisher_GAUSS_QUAD = mod_struct->fisher_GAUSS_QUAD;
		user_parameters[i]->fisher_log10F = mod_struct->fisher_log10F;
		user_parameters[i]->fisher_freq= mod_struct->fisher_freq;
		if(mod_struct->fisher_weights){
			user_parameters[i]->fisher_weights= mod_struct->fisher_weights;
		}
		else{
			user_parameters[i]->fisher_weights = new double*[num_detectors];
			for(int j = 0 ; j<num_detectors; j++){
				user_parameters[i]->fisher_weights[j]=NULL;
			}
		}	
		user_parameters[i]->fisher_PSD= mod_struct->fisher_PSD;
		user_parameters[i]->fisher_length= mod_struct->fisher_length;

		user_parameters[i]->mod_struct = mod_struct;

		//user_parameters[i]->burn_freqs = mcmc_frequencies;
		//user_parameters[i]->burn_data = mcmc_data;
		//user_parameters[i]->burn_noise = mcmc_noise;
		//user_parameters[i]->burn_lengths = mcmc_data_length;
	}
	//######################################################

	//###########################################################
	//double trial[11] = {1,.9,.2,.9,.2,3,7,4,.24,0,0};
	//mcmc_data_interface i;
	//i.min_dim = 11; 
	//i.max_dim = 11; 
	//i.chain_id  = 0; 
	//i.chain_number  = 11; 
	//double ll = MCMC_likelihood_wrapper(trial, &i ,user_parameters[0]);
	//debugger_print(__FILE__,__LINE__,ll);
	//###########################################################
	continue_PTMCMC_MH_dynamic_PT_alloc_uncorrelated(checkpoint_file_start,sampler_output,output,  N_steps,  
		max_chain_N_thermo_ensemble, chain_temps, 
		swp_freq, t0, nu, max_chunk_size,chain_distribution_scheme,
		log_prior,MCMC_likelihood_wrapper, MCMC_fisher_wrapper,(void **)user_parameters,numThreads, pool, 
		//log_prior,MCMC_likelihood_wrapper, NULL,user_parameters,numThreads, pool, 
		show_prog,statistics_filename,
		chain_filename, likelihood_log_filename,checkpoint_filename);
	
	if(!mod_struct->fisher_weights){
		for(int i = 0 ;i<chain_N;i++){
			delete[] user_parameters[i]->fisher_weights;
		}
	}
	if(!mod_struct->weights){
		for(int i = 0 ;i<chain_N;i++){
			delete[] user_parameters[i]->weights;
		}
	}
	
	//Deallocate fftw plans
	for (int i =0;i<num_detectors;i++)
		deallocate_FFTW_mem(&plans[i]);
	free(plans);
	//#################################################
	for(int i = 0 ; i<num_detectors; i++){
		delete [] burn_data[i];
		delete [] burn_freqs[i];
		delete [] burn_noise[i];
		deallocate_FFTW_mem(&burn_plans[i]);
	}
	//if(mod_struct->GAUSS_QUAD){
	//	delete [] fish_freqs;
	//	delete [] fish_weights;
	//	for(int i = 0 ;i<mcmc_num_detectors; i++){
	//		delete [] fish_psd[i];
	//	}
	//	delete [] fish_psd;
	//}
	delete [] burn_data;
	delete [] burn_lengths;
	delete [] burn_noise;
	delete [] burn_freqs;
	delete [] burn_plans;
	for(int i = 0 ; i<chain_N; i++){
		delete user_parameters[i];
	}
	delete [] user_parameters;
	//#################################################
}
/*! \brief Takes in an MCMC checkpoint file and continues the chain
 *
 * Obviously, the user must be sure to correctly match the dimension, number of chains, the generation_method, 
 * the prior function, the data, psds, freqs, and the detectors (number and name), and the gps_time to the 
 * previous run, otherwise the behavior of the sampler is undefined.
 *
 * numThreads and pool do not necessarily have to be the same
 */
void PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW(mcmc_sampler_output *sampler_output,
	double **output,
	int dimension,
	int N_steps,
	int chain_N,
	int max_chain_N_thermo_ensemble,
	double *initial_pos,
	double *seeding_var,
	double **ensemble_initial_pos,
	double *chain_temps,
	int swp_freq,
	int t0,
	int nu,
	int max_chunk_size,
	std::string chain_distribution_scheme,
	double(*log_prior)(double *param, mcmc_data_interface *interface,void *parameters),
	int numThreads,
	bool pool,
	bool show_prog,
	int num_detectors,
	std::complex<double> **data,
	double **noise_psd,
	double **frequencies,
	int *data_length,
	double gps_time,
	std::string *detectors,
	MCMC_modification_struct *mod_struct,
	std::string generation_method,
	std::string statistics_filename,
	std::string chain_filename,
	std::string likelihood_log_filename,
	std::string checkpoint_filename
	)
{
	std::mutex fisher_mutex;
	std::cout<<"Waveform Approximant: "<<generation_method<<std::endl;
	std::cout<<"Dimension: "<<dimension<<std::endl;

	//Create fftw plan for each detector (length of data stream may be different)
	fftw_outline *plans= (fftw_outline *)malloc(sizeof(fftw_outline)*num_detectors);
	for (int i =0;i<num_detectors;i++)
	{	
		allocate_FFTW_mem_forward(&plans[i] , data_length[i]);
	}
	mcmc_noise = noise_psd;	
	mcmc_init_pos = initial_pos;
	mcmc_frequencies = frequencies;
	mcmc_data = data;
	mcmc_data_length = data_length;
	mcmc_detectors = detectors;
	mcmc_generation_method = generation_method;
	mcmc_fftw_plans = plans;
	mcmc_num_detectors = num_detectors;
	mcmc_gps_time = gps_time;
	mcmc_gmst = gps_to_GMST_radian(mcmc_gps_time);
	//mcmc_Nmod = mod_struct->ppE_Nmod;
	//mcmc_bppe = mod_struct->bppe;
	mcmc_mod_struct = mod_struct;
	mcmc_log_beta = false;
	mcmc_intrinsic = false;

	//To save time, intrinsic waveforms can be saved between detectors, if the 
	//frequencies are all the same
	mcmc_save_waveform = true;
	for(int i =1 ;i<mcmc_num_detectors; i++){
		if( mcmc_data_length[i] != mcmc_data_length[0] ||
			mcmc_frequencies[i][0]!= mcmc_frequencies[0][0] ||
			mcmc_frequencies[i][mcmc_data_length[i]-1] 
				!= mcmc_frequencies[0][mcmc_data_length[0]-1]){
			mcmc_save_waveform= false;
		}
			
	}

	bool local_seeding ;
	if(!seeding_var){
		local_seeding = true;
	}
	else{
		local_seeding = false;
	}

	double **seeding_var_ptr = &seeding_var;
	PTMCMC_method_specific_prep(generation_method, dimension, seeding_var_ptr, local_seeding);
	
	//######################################################
	int T = (int)(1./(mcmc_frequencies[0][1]-mcmc_frequencies[0][0]));
	debugger_print(__FILE__,__LINE__,T);
	int burn_factor = T/4; //Take all sources to 4 seconds
	debugger_print(__FILE__,__LINE__,burn_factor);
	std::complex<double> **burn_data = new std::complex<double>*[mcmc_num_detectors];
	double **burn_freqs = new double*[mcmc_num_detectors];
	double **burn_noise = new double*[mcmc_num_detectors];
	int *burn_lengths = new int[mcmc_num_detectors];
	fftw_outline *burn_plans= new fftw_outline[num_detectors];
	for(int j = 0; j<mcmc_num_detectors; j++){
		burn_lengths[j] = mcmc_data_length[j]/burn_factor;
		burn_data[j]= new std::complex<double>[burn_lengths[j]];
		burn_freqs[j]= new double[burn_lengths[j]];
		burn_noise[j]= new double[burn_lengths[j]];
		allocate_FFTW_mem_forward(&burn_plans[j], burn_lengths[j]);
		int ct = 0;
		for( int k = 0 ; k<mcmc_data_length[j]; k++){
			if(k%burn_factor==0 && ct<burn_lengths[j]){
				burn_data[j][ct] = mcmc_data[j][k];
				burn_freqs[j][ct] = mcmc_frequencies[j][k];
				burn_noise[j][ct] = mcmc_noise[j][k];
				ct++;
			}
		}
	}
	

	//######################################################
	//######################################################
	//Fishers sometimes need AD, but that's slow and single 
	//threaded -- use GAUSS quad with precomputed weights 
	//and abscissa 
	//double *fish_freqs = NULL;
	//double *fish_weights = NULL;
	//double **fish_psd = NULL; //Needs to be interpolated from data given
	//int fish_length = 100;
	//double flow = mcmc_frequencies[0][0];
	//double fhigh = mcmc_frequencies[0][mcmc_data_length[0]-1];
	//if(mod_struct->GAUSS_QUAD){
	//	fish_freqs = new double[fish_length];	
	//	fish_weights = new double[fish_length];	
	//	fish_psd = new double*[mcmc_num_detectors];	
	//	
	//	gauleg(log10(flow), log10(fhigh), fish_freqs, fish_weights, fish_length);
	//	for(int i = 0 ; i<fish_length; i++){
	//		fish_freqs[i]=pow(10,fish_freqs[i]);
	//	}
	//	for(int i = 0 ; i<mcmc_num_detectors; i++){
	//		fish_psd[i] = new double[fish_length];
	//		gsl_interp_accel *accel = gsl_interp_accel_alloc();
	//		gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear, mcmc_data_length[i]);
	//		gsl_spline_init(spline, mcmc_frequencies[i], mcmc_noise[i],mcmc_data_length[i]);
	//		for(int j = 0 ; j<fish_length; j++){
	//			fish_psd[i][j] = gsl_spline_eval(spline, fish_freqs[j],accel);
	//		}
	//		gsl_spline_free(spline);
	//	}
	//}
	//######################################################
	MCMC_user_param **user_parameters=NULL;
	user_parameters = new MCMC_user_param*[chain_N];
	for(int i = 0 ;i<chain_N; i++){
		user_parameters[i] = new MCMC_user_param;
		
		user_parameters[i]->burn_data = burn_data;
		user_parameters[i]->burn_freqs = burn_freqs;
		user_parameters[i]->burn_noise = burn_noise;
		user_parameters[i]->burn_lengths = burn_lengths;
		user_parameters[i]->burn_plans=burn_plans;

		user_parameters[i]->mFish= &fisher_mutex;
		user_parameters[i]->GAUSS_QUAD= mod_struct->GAUSS_QUAD;
		user_parameters[i]->log10F = mod_struct->log10F;

		if(mod_struct->weights){
			user_parameters[i]->weights = mod_struct->weights;			
		}
		else{
			user_parameters[i]->weights = new double*[num_detectors];			
			for(int j = 0 ; j<num_detectors; j++){
				user_parameters[i]->weights[j]=NULL;
			}
		}

		user_parameters[i]->fisher_GAUSS_QUAD = mod_struct->fisher_GAUSS_QUAD;
		user_parameters[i]->fisher_log10F = mod_struct->fisher_log10F;
		user_parameters[i]->fisher_freq= mod_struct->fisher_freq;
		if(mod_struct->fisher_weights){
			user_parameters[i]->fisher_weights= mod_struct->fisher_weights;
		}
		else{
			user_parameters[i]->fisher_weights = new double*[num_detectors];
			for(int j = 0 ; j<num_detectors; j++){
				user_parameters[i]->fisher_weights[j]=NULL;
			}
		}	
		user_parameters[i]->fisher_PSD= mod_struct->fisher_PSD;
		user_parameters[i]->fisher_length= mod_struct->fisher_length;


		user_parameters[i]->mod_struct = mod_struct;

		//user_parameters[i]->burn_freqs = mcmc_frequencies;
		//user_parameters[i]->burn_data = mcmc_data;
		//user_parameters[i]->burn_noise = mcmc_noise;
		//user_parameters[i]->burn_lengths = mcmc_data_length;
	}
	//######################################################
	//###########################################################
	//double trial[11] = {1,.9,.2,.9,.2,3,7,4,.24,0,0};
	//mcmc_data_interface i;
	//i.min_dim = 11; 
	//i.max_dim = 11; 
	//i.chain_id  = 0; 
	//i.chain_number  = 11; 
	//double ll = MCMC_likelihood_wrapper(trial, &i ,user_parameters[0]);
	//debugger_print(__FILE__,__LINE__,ll);
	//###########################################################


	PTMCMC_MH_dynamic_PT_alloc_uncorrelated(sampler_output,output, dimension, N_steps, chain_N, 
		max_chain_N_thermo_ensemble,initial_pos,seeding_var,ensemble_initial_pos, chain_temps, 
		swp_freq, t0, nu,max_chunk_size,chain_distribution_scheme,
		log_prior,MCMC_likelihood_wrapper, MCMC_fisher_wrapper,(void**)user_parameters,numThreads, pool, 
		//log_prior,MCMC_likelihood_wrapper, NULL,(void**)user_parameters,numThreads, pool, 
		show_prog,statistics_filename,
		chain_filename, likelihood_log_filename,checkpoint_filename);
	if(!mod_struct->fisher_weights){
		for(int i = 0 ;i<chain_N;i++){
			delete[] user_parameters[i]->fisher_weights;
		}
	}
	if(!mod_struct->weights){
		for(int i = 0 ;i<chain_N;i++){
			delete[] user_parameters[i]->weights;
		}
	}
	
	//Deallocate fftw plans
	for (int i =0;i<num_detectors;i++)
		deallocate_FFTW_mem(&plans[i]);
	//#################################################
	for(int i = 0 ; i<num_detectors; i++){
		delete [] burn_data[i];
		delete [] burn_freqs[i];
		delete [] burn_noise[i];
		deallocate_FFTW_mem(&burn_plans[i]);
	}
	//if(mod_struct->GAUSS_QUAD){
	//	delete [] fish_freqs;
	//	delete [] fish_weights;
	//	for(int i = 0 ;i<mcmc_num_detectors; i++){
	//		delete [] fish_psd[i];
	//	}
	//	delete [] fish_psd;
	//}
	delete [] burn_data;
	delete [] burn_lengths;
	delete [] burn_noise;
	delete [] burn_freqs;
	delete [] burn_plans;
	for(int i = 0 ; i<chain_N; i++){
		delete user_parameters[i];
	}
	delete [] user_parameters;
	//#################################################
	free(plans);
	if(local_seeding){ delete [] seeding_var;}
}

/*! \brief Takes in an MCMC checkpoint file and continues the chain
 *
 * Obviously, the user must be sure to correctly match the dimension, number of chains, the generation_method, 
 * the prior function, the data, psds, freqs, and the detectors (number and name), and the gps_time to the 
 * previous run, otherwise the behavior of the sampler is undefined.
 *
 * numThreads and pool do not necessarily have to be the same
 */
void PTMCMC_MH_dynamic_PT_alloc_GW(double ***output,
	int dimension,
	int N_steps,
	int chain_N,
	int max_chain_N_thermo_ensemble,
	double *initial_pos,
	double *seeding_var,
	double **ensemble_initial_pos,
	double *chain_temps,
	int swp_freq,
	int t0,
	int nu,
	std::string chain_distribution_scheme,
	double(*log_prior)(double *param, mcmc_data_interface *interface,void *parameters),
	int numThreads,
	bool pool,
	bool show_prog,
	int num_detectors,
	std::complex<double> **data,
	double **noise_psd,
	double **frequencies,
	int *data_length,
	double gps_time,
	std::string *detectors,
	MCMC_modification_struct *mod_struct,
	std::string generation_method,
	std::string statistics_filename,
	std::string chain_filename,
	std::string likelihood_log_filename,
	std::string checkpoint_filename
	)
{
	//Create fftw plan for each detector (length of data stream may be different)
	fftw_outline *plans= (fftw_outline *)malloc(sizeof(fftw_outline)*num_detectors);
	for (int i =0;i<num_detectors;i++)
	{	
		allocate_FFTW_mem_forward(&plans[i] , data_length[i]);
	}
	void **user_parameters=NULL;
	mcmc_noise = noise_psd;	
	mcmc_frequencies = frequencies;
	mcmc_data = data;
	mcmc_data_length = data_length;
	mcmc_detectors = detectors;
	mcmc_generation_method = generation_method;
	mcmc_fftw_plans = plans;
	mcmc_num_detectors = num_detectors;
	mcmc_gps_time = gps_time;
	mcmc_gmst = gps_to_GMST_radian(mcmc_gps_time);
	//mcmc_Nmod = Nmod;
	//mcmc_bppe = bppe;
	mcmc_mod_struct = mod_struct;
	mcmc_log_beta = false;
	mcmc_intrinsic = false;

	//To save time, intrinsic waveforms can be saved between detectors, if the 
	//frequencies are all the same
	mcmc_save_waveform = true;
	for(int i =1 ;i<mcmc_num_detectors; i++){
		if( mcmc_data_length[i] != mcmc_data_length[0] ||
			mcmc_frequencies[i][0]!= mcmc_frequencies[0][0] ||
			mcmc_frequencies[i][mcmc_data_length[i]-1] 
				!= mcmc_frequencies[0][mcmc_data_length[0]-1]){
			mcmc_save_waveform= false;
		}
			
	}

	bool local_seeding ;
	if(!seeding_var)
		local_seeding = true;
	else
		local_seeding = false;

	double **seeding_var_ptr = &seeding_var;
	PTMCMC_method_specific_prep(generation_method, dimension, seeding_var_ptr, local_seeding);

	PTMCMC_MH_dynamic_PT_alloc(output, dimension, N_steps, chain_N, 
		max_chain_N_thermo_ensemble,initial_pos,seeding_var,ensemble_initial_pos, chain_temps, 
		swp_freq, t0, nu, chain_distribution_scheme,
		log_prior,MCMC_likelihood_wrapper, MCMC_fisher_wrapper,user_parameters,numThreads, pool, 
		//log_prior,MCMC_likelihood_wrapper, NULL,numThreads, pool, 
		show_prog,statistics_filename,
		chain_filename, checkpoint_filename);
	
	//Deallocate fftw plans
	for (int i =0;i<num_detectors;i++)
		deallocate_FFTW_mem(&plans[i]);
	free(plans);
	if(local_seeding){ delete [] seeding_var;}
}



/*! \brief Takes in an MCMC checkpoint file and continues the chain
 *
 * Obviously, the user must be sure to correctly match the dimension, number of chains, the generation_method, 
 * the prior function, the data, psds, freqs, and the detectors (number and name), and the gps_time to the 
 * previous run, otherwise the behavior of the sampler is undefined.
 *
 * numThreads and pool do not necessarily have to be the same
 */
void continue_PTMCMC_MH_GW(std::string start_checkpoint_file,
	double ***output,
	int dimension,
	int N_steps,
	int swp_freq,
	double(*log_prior)(double *param, mcmc_data_interface *interface ,void *parameters),
	int numThreads,
	bool pool,
	bool show_prog,
	int num_detectors,
	std::complex<double> **data,
	double **noise_psd,
	double **frequencies,
	int *data_length,
	double gps_time,
	std::string *detectors,
	MCMC_modification_struct *mod_struct,
	std::string generation_method,
	std::string statistics_filename,
	std::string chain_filename,
	std::string auto_corr_filename,
	std::string likelihood_log_filename,
	std::string final_checkpoint_filename,
	bool tune
	)
{
	//Create fftw plan for each detector (length of data stream may be different)
	fftw_outline *plans= (fftw_outline *)malloc(sizeof(fftw_outline)*num_detectors);
	for (int i =0;i<num_detectors;i++)
	{	
		allocate_FFTW_mem_forward(&plans[i] , data_length[i]);
	}
	void **user_parameters=NULL;
	mcmc_noise = noise_psd;	
	mcmc_frequencies = frequencies;
	mcmc_data = data;
	mcmc_data_length = data_length;
	mcmc_detectors = detectors;
	mcmc_generation_method = generation_method;
	mcmc_fftw_plans = plans;
	mcmc_num_detectors = num_detectors;
	mcmc_gps_time = gps_time;
	mcmc_gmst = gps_to_GMST_radian(mcmc_gps_time);
	//mcmc_Nmod = Nmod;
	//mcmc_bppe = bppe;
	mcmc_mod_struct = mod_struct;
	mcmc_log_beta = false;
	mcmc_intrinsic = false;

	//To save time, intrinsic waveforms can be saved between detectors, if the 
	//frequencies are all the same
	mcmc_save_waveform = true;
	for(int i =1 ;i<mcmc_num_detectors; i++){
		if( mcmc_data_length[i] != mcmc_data_length[0] ||
			mcmc_frequencies[i][0]!= mcmc_frequencies[0][0] ||
			mcmc_frequencies[i][mcmc_data_length[i]-1] 
				!= mcmc_frequencies[0][mcmc_data_length[0]-1]){
			mcmc_save_waveform= false;
		}
			
	}

	bool local_seeding=false ;

	double **seeding_var_ptr = NULL;
	PTMCMC_method_specific_prep(generation_method, dimension, seeding_var_ptr, local_seeding);

	continue_PTMCMC_MH(start_checkpoint_file,output, N_steps,swp_freq,log_prior,
			MCMC_likelihood_wrapper, MCMC_fisher_wrapper,user_parameters,numThreads, pool, 
			//MCMC_likelihood_wrapper, NULL,numThreads, pool, 
			show_prog,statistics_filename,chain_filename,
			final_checkpoint_filename,tune);
	//Deallocate fftw plans
	for (int i =0;i<num_detectors;i++)
		deallocate_FFTW_mem(&plans[i]);
	free(plans);
}

/*! \brief Unpacks MCMC parameters for method specific initiation 
 *
 * Populates seeding vector if non supplied, populates mcmc_Nmod, populates mcmc_log_beta, populates mcmc_intrinsic
 */
void PTMCMC_method_specific_prep(std::string generation_method, int dimension,double **seeding_var, bool local_seeding)
{
	int totalmod = (mcmc_mod_struct->gIMR_Nmod_phi + mcmc_mod_struct->gIMR_Nmod_sigma + mcmc_mod_struct->gIMR_Nmod_beta + mcmc_mod_struct->gIMR_Nmod_alpha  + mcmc_mod_struct->ppE_Nmod);
	//if(generation_method.find("EA") != std::string::npos){totalmod+=4;}
	if(generation_method.find("EA") != std::string::npos){
	        totalmod+=3;
	}
	if(generation_method.find("PhenomD") != std::string::npos && (dimension - totalmod) == 4)
	{
		std::cout<<"Sampling in parameters: ln chirpmass, eta, chi1, chi2";
		for(int i =0; i<totalmod; i++){
			std::cout<<", mod_"<<i;
		}
		std::cout<<std::endl;
		if(local_seeding){
			(*seeding_var) = new double[dimension];
			(*seeding_var)[0]=.5;
			(*seeding_var)[1]=.1;
			(*seeding_var)[2]=.1;
			(*seeding_var)[3]=.1;
			for(int i = 0  ;i <totalmod; i++){
				(*seeding_var)[i+4] = 1;
			}
		}
		mcmc_intrinsic=true;
	} 
	else if(generation_method.find("PhenomD_NRT") != std::string::npos &&   (dimension - totalmod) == 6)
	{
		std::cout<<"Sampling in parameters: ln chirpmass, eta, chi1, chi2, tidal1,  tidal2";
		for(int i =0; i<totalmod; i++){
			std::cout<<", mod_"<<i;
		}
		std::cout<<std::endl;
		if(local_seeding){
			(*seeding_var) = new double[dimension];
			(*seeding_var)[0]=.5;
			(*seeding_var)[1]=.1;
			(*seeding_var)[2]=.1;
			(*seeding_var)[3]=.1;
			(*seeding_var)[4]=10;
			(*seeding_var)[5]=10;
			for(int i = 0  ;i <totalmod; i++){
				(*seeding_var)[i+6] = 1;
			}
		}
		mcmc_intrinsic=true;
	} 
	else if(generation_method.find("PhenomD_NRT") != std::string::npos &&   (dimension - totalmod) == 5)
	{
		std::cout<<"Sampling in parameters: ln chirpmass, eta, chi1, chi2, tidal_s";
		for(int i =0; i<totalmod; i++){
			std::cout<<", mod_"<<i;
		}
		std::cout<<std::endl;
		if(local_seeding){
			(*seeding_var) = new double[dimension];
			(*seeding_var)[0]=.5;
			(*seeding_var)[1]=.1;
			(*seeding_var)[2]=.1;
			(*seeding_var)[3]=.1;
			(*seeding_var)[4]=10;
			for(int i = 0  ;i <totalmod; i++){
				(*seeding_var)[i+5] = 1;
			}
		}
		mcmc_intrinsic=true;
	} 
	else if(generation_method.find("PhenomPv2") != std::string::npos && (dimension - totalmod) == 8)
	{
		std::cout<<"Sampling in parameters: ln chirpmass, eta, a1, a2, tilt1, tilt2, phi1, phi2";
		for(int i =0; i<totalmod; i++){
			std::cout<<", mod_"<<i;
		}
		std::cout<<std::endl;
		if(local_seeding){
			(*seeding_var) = new double[dimension];
			(*seeding_var)[0]=.5;
			(*seeding_var)[1]=.1;
			(*seeding_var)[2]=.1;
			(*seeding_var)[3]=.1;
			(*seeding_var)[4]=.1;
			(*seeding_var)[5]=.1;
			(*seeding_var)[6]=.1;
			(*seeding_var)[7]=.1;
			for(int i = 0  ;i <totalmod; i++){
				(*seeding_var)[i+8] = 1;
			}
		}
		mcmc_intrinsic=true;
	} 
	else if(generation_method.find("PhenomD") != std::string::npos && (dimension - totalmod) == 11)
	{
		std::cout<<"Sampling in parameters: RA, sin DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, chi1, chi2"<<std::endl;
		for(int i =0; i<totalmod; i++){
			std::cout<<", mod_"<<i;
		}
		std::cout<<std::endl;
		if(local_seeding){
			(*seeding_var) = new double[dimension];
			(*seeding_var)[0]=.1;
			(*seeding_var)[1]=.1;
			(*seeding_var)[2]=.1;
			(*seeding_var)[3]=.1;
			(*seeding_var)[4]=.5;
			(*seeding_var)[5]=.1;
			(*seeding_var)[6]=.1;
			(*seeding_var)[7]=.1;
			(*seeding_var)[8]=.1;
			(*seeding_var)[9]=.1;
			(*seeding_var)[10]=.5;
			for(int i = 0  ;i <totalmod; i++){
				(*seeding_var)[i+11] = 1;
			}
		}
		mcmc_intrinsic=false;
	} 
	else if(generation_method.find("PhenomD_NRT") != std::string::npos && (dimension - totalmod) == 13)
	{
		std::cout<<"Sampling in parameters: RA, sin  DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, chi1, chi2, tidal1, tidal2"<<std::endl;
		for(int i =0; i<totalmod; i++){
			std::cout<<", mod_"<<i;
		}
		std::cout<<std::endl;
		if(local_seeding){
			(*seeding_var) = new double[dimension];
			(*seeding_var)[0]=.1;
			(*seeding_var)[1]=.1;
			(*seeding_var)[2]=.1;
			(*seeding_var)[3]=.1;
			(*seeding_var)[4]=.5;
			(*seeding_var)[5]=.1;
			(*seeding_var)[6]=.1;
			(*seeding_var)[7]=.1;
			(*seeding_var)[8]=.1;
			(*seeding_var)[9]=.1;
			(*seeding_var)[10]=.5;
			(*seeding_var)[11]=10;
			(*seeding_var)[12]=10;
			for(int i = 0  ;i <totalmod; i++){
				(*seeding_var)[i+13] = 1;
			}
		}
		mcmc_intrinsic=false;
	} 
	else if(generation_method.find("PhenomD_NRT") != std::string::npos && (dimension - totalmod) == 12)
	{
		std::cout<<"Sampling in parameters: RA, sin  DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, chi1, chi2, tidal_s"<<std::endl;
		for(int i =0; i<totalmod; i++){
			std::cout<<", mod_"<<i;
		}
		std::cout<<std::endl;
		if(local_seeding){
			(*seeding_var) = new double[dimension];
			(*seeding_var)[0]=.1;
			(*seeding_var)[1]=.1;
			(*seeding_var)[2]=.1;
			(*seeding_var)[3]=.1;
			(*seeding_var)[4]=.5;
			(*seeding_var)[5]=.1;
			(*seeding_var)[6]=.1;
			(*seeding_var)[7]=.1;
			(*seeding_var)[8]=.1;
			(*seeding_var)[9]=.1;
			(*seeding_var)[10]=.5;
			(*seeding_var)[11]=10;
			//(*seeding_var)[12]=10;
			for(int i = 0  ;i <totalmod; i++){
				(*seeding_var)[i+12] = 1;
			}
		}
		mcmc_intrinsic=false;
	} 
	else if(generation_method.find("PhenomPv2") != std::string::npos && (dimension - totalmod) == 15)
	{
		std::cout<<"Sampling in parameters: RA, sin DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, a1, a2,cos tilt1, cos tilt2, phi1, phi2"<<std::endl;
		for(int i =0; i<totalmod; i++){
			std::cout<<", mod_"<<i;
		}
		std::cout<<std::endl;
		if(local_seeding){
			(*seeding_var) = new double[dimension];
			(*seeding_var)[0]=1.;
			(*seeding_var)[1]=.5;
			(*seeding_var)[2]=1;
			(*seeding_var)[3]=1;
			(*seeding_var)[4]=1;
			(*seeding_var)[5]=.1;
			(*seeding_var)[6]=1;
			(*seeding_var)[7]=1.;
			(*seeding_var)[8]=.3;
			(*seeding_var)[9]=.5;
			(*seeding_var)[10]=.5;
			(*seeding_var)[11]=1.;
			(*seeding_var)[12]=1.;
			(*seeding_var)[13]=1.;
			(*seeding_var)[14]=1.;
			for(int i = 0  ;i <totalmod; i++){
				(*seeding_var)[i+15] = 1;
			}
		}
		mcmc_intrinsic=false;
	} 
	else{
		std::cout<<
			"Input parameters not valid, please check that input is compatible with the supported methods - dimension combinations"<<std::endl;
		std::cout<<"Using: "<<generation_method<<" "<<dimension<<std::endl;
		exit(1);
	}
	//if(dimension==4 && generation_method =="IMRPhenomD"){
	//	std::cout<<"Sampling in parameters: ln chirpmass, eta, chi1, chi2"<<std::endl;
	//	if(local_seeding){
	//		(*seeding_var) = new double[dimension];
	//		(*seeding_var)[0]=.5;
	//		(*seeding_var)[1]=.1;
	//		(*seeding_var)[2]=.1;
	//		(*seeding_var)[3]=.1;
	//	}
	//	mcmc_intrinsic=true;
	//}
	//else if(dimension==6 && generation_method =="IMRPhenomD_NRT"){
	//	std::cout<<"Sampling in parameters: ln chirpmass, eta, chi1, chi2, lambda1, lambda2"<<std::endl;
	//	if(local_seeding){
	//		(*seeding_var) = new double[dimension];
	//		(*seeding_var)[0]=.5;
	//		(*seeding_var)[1]=.1;
	//		(*seeding_var)[2]=.1;
	//		(*seeding_var)[3]=.1;
	//		(*seeding_var)[4]=.1;
	//		(*seeding_var)[5]=.1;
	//	}
	//	mcmc_intrinsic=true;
	//}
	//else if(generation_method =="SkySearch"){
	//	std::cout<<"Sampling in parameters: "<<std::endl;
	//	if(local_seeding){
	//		(*seeding_var) = new double[dimension];
	//		(*seeding_var)[0]=1;
	//		(*seeding_var)[1]=.3;
	//		(*seeding_var)[2]=1.;
	//		(*seeding_var)[3]=1.;
	//		(*seeding_var)[4]=1.;
	//		(*seeding_var)[5]=20;
	//		(*seeding_var)[6]=20;
	//	}
	//	mcmc_intrinsic=false;
	//}
	//else if(dimension==11 && generation_method =="IMRPhenomD"){
	//	std::cout<<"Sampling in parameters: RA, DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, chi1, chi2"<<std::endl;
	//	mcmc_intrinsic=false;
	//	if(local_seeding){
	//		(*seeding_var) = new double[dimension];
	//		(*seeding_var)[0]=.1;
	//		(*seeding_var)[1]=.1;
	//		(*seeding_var)[2]=.1;
	//		(*seeding_var)[3]=.1;
	//		(*seeding_var)[4]=.5;
	//		(*seeding_var)[5]=.1;
	//		(*seeding_var)[6]=.1;
	//		(*seeding_var)[7]=.1;
	//		(*seeding_var)[8]=.1;
	//		(*seeding_var)[9]=.1;
	//		(*seeding_var)[10]=.5;
	//	}
	//}
	//else if(dimension==13 && generation_method =="IMRPhenomD_NRT"){
	//	std::cout<<"Sampling in parameters: RA, DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, chi1, chi2"<<std::endl;
	//	mcmc_intrinsic=false;
	//	if(local_seeding){
	//		(*seeding_var) = new double[dimension];
	//		(*seeding_var)[0]=.1;
	//		(*seeding_var)[1]=.1;
	//		(*seeding_var)[2]=.1;
	//		(*seeding_var)[3]=.1;
	//		(*seeding_var)[4]=.5;
	//		(*seeding_var)[5]=.1;
	//		(*seeding_var)[6]=.1;
	//		(*seeding_var)[7]=.1;
	//		(*seeding_var)[8]=.1;
	//		(*seeding_var)[9]=.1;
	//		(*seeding_var)[10]=.5;
	//		(*seeding_var)[11]=20;
	//		(*seeding_var)[12]=20;
	//	}
	//}
	//else if(dimension==12 && generation_method =="dCS_IMRPhenomD"){
	//	mcmc_mod_struct->ppE_Nmod = 1;
	//	std::cout<<"Sampling in parameters: RA, DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, chi1, chi2, sqrt(alpha) (km)"<<std::endl;
	//	if(local_seeding){
	//		(*seeding_var) = new double[dimension];
	//		(*seeding_var)[0]=1.;
	//		(*seeding_var)[1]=.5;
	//		(*seeding_var)[2]=1;
	//		(*seeding_var)[3]=1;
	//		(*seeding_var)[4]=1;
	//		(*seeding_var)[5]=.1;
	//		(*seeding_var)[6]=1;
	//		(*seeding_var)[7]=1.;
	//		(*seeding_var)[8]=.3;
	//		(*seeding_var)[9]=.5;
	//		(*seeding_var)[10]=.5;
	//		(*seeding_var)[11]=1;
	//	}
	//}
	////Sampling parameter is root alpha in km
	//else if(dimension==12 && generation_method =="EdGB_IMRPhenomD"){
	//	mcmc_mod_struct->ppE_Nmod = 1;
	//	std::cout<<"Sampling in parameters: RA, DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, chi1, chi2, sqrt(alpha) (km)"<<std::endl;
	//	if(local_seeding){
	//		(*seeding_var) = new double[dimension];
	//		(*seeding_var)[0]=1.;
	//		(*seeding_var)[1]=.5;
	//		(*seeding_var)[2]=1;
	//		(*seeding_var)[3]=1;
	//		(*seeding_var)[4]=1;
	//		(*seeding_var)[5]=.1;
	//		(*seeding_var)[6]=1;
	//		(*seeding_var)[7]=1.;
	//		(*seeding_var)[8]=.3;
	//		(*seeding_var)[9]=.5;
	//		(*seeding_var)[10]=.5;
	//		(*seeding_var)[11]=1;
	//	}
	//}
	//else if(dimension==13 && 
	//	(generation_method =="EdGB_GHOv1_IMRPhenomD" || generation_method =="EdGB_GHOv2_IMRPhenomD" || generation_method =="EdGB_GHOv3_IMRPhenomD" )){
	//	mcmc_mod_struct->ppE_Nmod = 2;
	//	std::cout<<"Sampling in parameters: RA, DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, chi1, chi2, sqrt(alpha) (km), gamma (generic 1PN term)"<<std::endl;
	//	if(local_seeding){
	//		(*seeding_var) = new double[dimension];
	//		(*seeding_var)[0]=1.;
	//		(*seeding_var)[1]=.5;
	//		(*seeding_var)[2]=1;
	//		(*seeding_var)[3]=1;
	//		(*seeding_var)[4]=1;
	//		(*seeding_var)[5]=.1;
	//		(*seeding_var)[6]=1;
	//		(*seeding_var)[7]=1.;
	//		(*seeding_var)[8]=.3;
	//		(*seeding_var)[9]=.5;
	//		(*seeding_var)[10]=.5;
	//		(*seeding_var)[11]=1;
	//		(*seeding_var)[12]=1;
	//	}
	//}
	//else if(dimension==16 && generation_method =="dCS_IMRPhenomPv2"){
	//	mcmc_mod_struct->ppE_Nmod = 1;
	//	std::cout<<"Sampling in parameters: RA, DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, a1, a2,cos tilt1, cos tilt2, phi1, phi2, sqrt(alpha) (km)"<<std::endl;
	//	if(local_seeding){
	//		(*seeding_var) = new double[dimension];
	//		(*seeding_var)[0]=1.;
	//		(*seeding_var)[1]=.5;
	//		(*seeding_var)[2]=1;
	//		(*seeding_var)[3]=1;
	//		(*seeding_var)[4]=1;
	//		(*seeding_var)[5]=.1;
	//		(*seeding_var)[6]=1;
	//		(*seeding_var)[7]=1.;
	//		(*seeding_var)[8]=.3;
	//		(*seeding_var)[9]=.5;
	//		(*seeding_var)[10]=.5;
	//		(*seeding_var)[11]=1.;
	//		(*seeding_var)[12]=1.;
	//		(*seeding_var)[13]=1.;
	//		(*seeding_var)[14]=1.;
	//		(*seeding_var)[15]=1.;
	//	}
	//}
	////Sampling parameter is root alpha in km
	//else if(dimension==16 && generation_method =="EdGB_IMRPhenomPv2"){
	//	mcmc_mod_struct->ppE_Nmod = 1;
	//	std::cout<<"Sampling in parameters: RA, DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, a1, a2,cos tilt1, cos tilt2, phi1, phi2, sqrt(alpha) (km)"<<std::endl;
	//	if(local_seeding){
	//		(*seeding_var) = new double[dimension];
	//		(*seeding_var)[0]=1.;
	//		(*seeding_var)[1]=.5;
	//		(*seeding_var)[2]=1;
	//		(*seeding_var)[3]=1;
	//		(*seeding_var)[4]=1;
	//		(*seeding_var)[5]=.1;
	//		(*seeding_var)[6]=1;
	//		(*seeding_var)[7]=1.;
	//		(*seeding_var)[8]=.3;
	//		(*seeding_var)[9]=.5;
	//		(*seeding_var)[10]=.5;
	//		(*seeding_var)[11]=1.;
	//		(*seeding_var)[12]=1.;
	//		(*seeding_var)[13]=1.;
	//		(*seeding_var)[14]=1.;
	//		(*seeding_var)[15]=1.;
	//	}
	//}
	//else if(dimension==17 && 
	//	(generation_method =="EdGB_GHOv1_IMRPhenomPv2" || generation_method =="EdGB_GHOv2_IMRPhenomPv2" || generation_method =="EdGB_GHOv3_IMRPhenomPv2" )){
	//	mcmc_mod_struct->ppE_Nmod = 2;
	//	std::cout<<"Sampling in parameters: RA, DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, a1, a2,cos tilt1, cos tilt2, phi1, phi2, sqrt(alpha) (km)"<<std::endl;
	//	if(local_seeding){
	//		(*seeding_var) = new double[dimension];
	//		(*seeding_var)[0]=1.;
	//		(*seeding_var)[1]=.5;
	//		(*seeding_var)[2]=1;
	//		(*seeding_var)[3]=1;
	//		(*seeding_var)[4]=1;
	//		(*seeding_var)[5]=.1;
	//		(*seeding_var)[6]=1;
	//		(*seeding_var)[7]=1.;
	//		(*seeding_var)[8]=.3;
	//		(*seeding_var)[9]=.5;
	//		(*seeding_var)[10]=.5;
	//		(*seeding_var)[11]=1.;
	//		(*seeding_var)[12]=1.;
	//		(*seeding_var)[13]=1.;
	//		(*seeding_var)[14]=1.;
	//		(*seeding_var)[15]=1.;
	//		(*seeding_var)[16]=1.;
	//	}
	//}
	//else if(dimension==15 && generation_method =="IMRPhenomPv2"){
	//	mcmc_intrinsic=false;
	//	std::cout<<"Sampling in parameters: RA, DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, a1, a2,cos tilt1, cos tilt2, phi1, phi2"<<std::endl;
	//	if(local_seeding){
	//		(*seeding_var) = new double[dimension];
	//		(*seeding_var)[0]=1.;
	//		(*seeding_var)[1]=.5;
	//		(*seeding_var)[2]=1;
	//		(*seeding_var)[3]=1;
	//		(*seeding_var)[4]=1;
	//		(*seeding_var)[5]=.1;
	//		(*seeding_var)[6]=1;
	//		(*seeding_var)[7]=1.;
	//		(*seeding_var)[8]=.3;
	//		(*seeding_var)[9]=.5;
	//		(*seeding_var)[10]=.5;
	//		(*seeding_var)[11]=1.;
	//		(*seeding_var)[12]=1.;
	//		(*seeding_var)[13]=1.;
	//		(*seeding_var)[14]=1.;
	//	}
	//}
	//else if(dimension==8 && generation_method =="IMRPhenomPv2"){
	//	mcmc_intrinsic=true;
	//	std::cout<<"Sampling in parameters: ln chirpmass, eta, a1, a2,cos tilt1, cos tilt2, phi1, phi2"<<std::endl;
	//	if(local_seeding){
	//		(*seeding_var) = new double[dimension];
	//		(*seeding_var)[0]=1.;
	//		(*seeding_var)[1]=.3;
	//		(*seeding_var)[2]=1.;
	//		(*seeding_var)[3]=1.;
	//		(*seeding_var)[4]=1.;
	//		(*seeding_var)[5]=1.;
	//		(*seeding_var)[6]=1.;
	//		(*seeding_var)[7]=1.;
	//	}
	//}
	//else if( generation_method == "ppE_IMRPhenomPv2_Inspiral"
	//	|| generation_method == "ppE_IMRPhenomPv2_IMR"
	//	|| generation_method == "PNSeries_ppE_IMRPhenomPv2_IMR"
	//	|| generation_method == "PNSeries_ppE_IMRPhenomPv2_Inspiral"
	//	){
	//	if(dimension-mcmc_mod_struct->ppE_Nmod == 15){
	//		
	//	std::cout<<"Sampling in parameters: RA, DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, a1, a2,cos tilt1, cos tilt2, phi1, phi2";
	//		for(int i =0; i<mcmc_Nmod; i++){
	//			std::cout<<", beta"<<i;
	//		}
	//		std::cout<<endl;
	//		if(local_seeding){
	//			(*seeding_var) = new double[dimension];
	//			(*seeding_var)[0]=1.;
	//			(*seeding_var)[1]=.5;
	//			(*seeding_var)[2]=1;
	//			(*seeding_var)[3]=1;
	//			(*seeding_var)[4]=1;
	//			(*seeding_var)[5]=.1;
	//			(*seeding_var)[6]=1;
	//			(*seeding_var)[7]=1.;
	//			(*seeding_var)[8]=.3;
	//			(*seeding_var)[9]=.5;
	//			(*seeding_var)[10]=.5;
	//			(*seeding_var)[11]=1.;
	//			(*seeding_var)[12]=1.;
	//			(*seeding_var)[13]=1.;
	//			(*seeding_var)[14]=1.;
	//			for(int i =0; i< mcmc_mod_struct->ppE_Nmod;i++){
	//				(*seeding_var)[15 + i]=2;
	//			}
	//		}
	//	}
	//	if(dimension-mcmc_mod_struct->ppE_Nmod == 8){
	//		mcmc_intrinsic=true;
	//		std::cout<<"Sampling in parameters: ln chirpmass, eta, a1, a2,cos tilt1, cos tilt2, phi1, phi2";
	//		for(int i =0; i<mcmc_mod_struct->ppE_Nmod; i++){
	//			std::cout<<", beta"<<i;
	//		}
	//		std::cout<<endl;
	//		if(local_seeding){
	//			(*seeding_var) = new double[dimension];
	//			(*seeding_var)[0]=1.;
	//			(*seeding_var)[1]=.3;
	//			(*seeding_var)[2]=1.;
	//			(*seeding_var)[3]=1.;
	//			(*seeding_var)[4]=1.;
	//			(*seeding_var)[5]=1.;
	//			(*seeding_var)[6]=1.;
	//			(*seeding_var)[7]=1.;
	//			for(int i =0; i< mcmc_mod_struct->ppE_Nmod;i++){
	//				(*seeding_var)[8 + i]=2;
	//			}
	//		}
	//	}

	//}
	//else if(generation_method == "gIMRPhenomPv2" 
	//	){
	//	int totalmod = (mcmc_mod_struct->gIMR_Nmod_phi + mcmc_mod_struct->gIMR_Nmod_sigma + mcmc_mod_struct->gIMR_Nmod_beta + mcmc_mod_struct->gIMR_Nmod_alpha );
	//	if(dimension-totalmod == 15){
	//		
	//		mcmc_intrinsic=false;
	//		std::cout<<"Sampling in parameters: RA, DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, a1, a2,cos tilt1, cos tilt2, phi1, phi2";
	//		for(int i =0; i<totalmod; i++){
	//			std::cout<<", mod "<<i;
	//		}
	//		std::cout<<endl;
	//		if(local_seeding){
	//			(*seeding_var) = new double[dimension];
	//			(*seeding_var)[0]=1.;
	//			(*seeding_var)[1]=.5;
	//			(*seeding_var)[2]=1;
	//			(*seeding_var)[3]=1;
	//			(*seeding_var)[4]=1;
	//			(*seeding_var)[5]=.1;
	//			(*seeding_var)[6]=1;
	//			(*seeding_var)[7]=1.;
	//			(*seeding_var)[8]=.3;
	//			(*seeding_var)[9]=.5;
	//			(*seeding_var)[10]=.5;
	//			(*seeding_var)[11]=1.;
	//			(*seeding_var)[12]=1.;
	//			(*seeding_var)[13]=1.;
	//			(*seeding_var)[14]=1.;
	//			for(int i =0; i< totalmod;i++){
	//				(*seeding_var)[15 + i]=2;
	//			}
	//		}
	//	}
	//	else{
	//		mcmc_intrinsic=true;
	//		std::cout<<"Sampling in parameters: ln chirpmass, eta, a1, a2,cos tilt1, cos tilt2, phi1, phi2";
	//		for(int i =0; i<totalmod; i++){
	//			std::cout<<", mod "<<i;
	//		}
	//		if(local_seeding){
	//			(*seeding_var) = new double[dimension];
	//			(*seeding_var)[0]=1.;
	//			(*seeding_var)[1]=.3;
	//			(*seeding_var)[2]=1.;
	//			(*seeding_var)[3]=1.;
	//			(*seeding_var)[4]=1.;
	//			(*seeding_var)[5]=1.;
	//			(*seeding_var)[6]=1.;
	//			(*seeding_var)[7]=1.;
	//			for(int i =0; i< totalmod;i++){
	//				(*seeding_var)[8 + i]=2;
	//			}
	//		}
	//	}

	//}
	//else if(generation_method == "gIMRPhenomD" 
	//	){
	//	int totalmod = (mcmc_mod_struct->gIMR_Nmod_phi + mcmc_mod_struct->gIMR_Nmod_sigma + mcmc_mod_struct->gIMR_Nmod_beta + mcmc_mod_struct->gIMR_Nmod_alpha );
	//	if(dimension-totalmod == 4){
	//		std::cout<<"Sampling in parameters: ln chirpmass, eta, chi1, chi2";
	//		
	//		mcmc_intrinsic=true;
	//		for(int i =0; i<totalmod; i++){
	//			std::cout<<", beta"<<i;
	//		}
	//		
	//		std::cout<<endl;
	//		if(local_seeding){
	//			(*seeding_var) = new double[dimension];
	//			(*seeding_var)[0]=1.;
	//			(*seeding_var)[1]=.3;
	//			(*seeding_var)[2]=1;
	//			(*seeding_var)[3]=1;
	//			for(int i =0; i< totalmod;i++){
	//				(*seeding_var)[4 + i]=2;
	//			}
	//		}
	//	}
	//	else{
	//		std::cout<<"Sampling in parameters: RA, DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, chi1, chi2";
	//		mcmc_intrinsic=false;
	//		if(local_seeding){
	//			(*seeding_var) = new double[dimension];
	//			(*seeding_var)[0]=1.;
	//			(*seeding_var)[1]=.5;
	//			(*seeding_var)[2]=1;
	//			(*seeding_var)[3]=1;
	//			(*seeding_var)[4]=1;
	//			(*seeding_var)[5]=.1;
	//			(*seeding_var)[6]=1;
	//			(*seeding_var)[7]=1.;
	//			(*seeding_var)[8]=.3;
	//			(*seeding_var)[9]=.5;
	//			(*seeding_var)[10]=.5;
	//			for(int i =0; i< totalmod;i++){
	//				(*seeding_var)[11 + i]=2;
	//			}
	//		}
	//	}

	//}
	//else if(generation_method == "ppE_IMRPhenomD_Inspiral"
	//	|| generation_method == "ppE_IMRPhenomD_IMR"
	//	|| generation_method == "PNSeries_ppE_IMRPhenomD_IMR"
	//	|| generation_method == "PNSeries_ppE_IMRPhenomD_Inspiral"
	//	){
	//	if(dimension-mcmc_mod_struct->ppE_Nmod == 11){
	//		
	//		std::cout<<"Sampling in parameters: RA, DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, chi1, chi2";
	//		for(int i =0; i<mcmc_mod_struct->ppE_Nmod; i++){
	//			std::cout<<", beta"<<i;
	//		}
	//		std::cout<<endl;
	//		if(local_seeding){
	//			(*seeding_var) = new double[dimension];
	//			(*seeding_var)[0]=1.;
	//			(*seeding_var)[1]=.5;
	//			(*seeding_var)[2]=1;
	//			(*seeding_var)[3]=1;
	//			(*seeding_var)[4]=1;
	//			(*seeding_var)[5]=.1;
	//			(*seeding_var)[6]=1;
	//			(*seeding_var)[7]=1.;
	//			(*seeding_var)[8]=.3;
	//			(*seeding_var)[9]=.5;
	//			(*seeding_var)[10]=.5;
	//			for(int i =0; i< mcmc_mod_struct->ppE_Nmod;i++){
	//				(*seeding_var)[11 + i]=2;
	//			}
	//		}
	//	}
	//	if(dimension-mcmc_mod_struct->ppE_Nmod == 4){
	//		mcmc_intrinsic=true;
	//		std::cout<<"Sampling in parameters: ln chirpmass, eta, chi1, chi2";
	//		for(int i =0; i<mcmc_mod_struct->ppE_Nmod; i++){
	//			std::cout<<", beta"<<i;
	//		}
	//		std::cout<<endl;
	//		if(local_seeding){
	//			(*seeding_var) = new double[dimension];
	//			(*seeding_var)[0]=.1;
	//			(*seeding_var)[1]=.1;
	//			(*seeding_var)[2]=.1;
	//			(*seeding_var)[3]=1;
	//			(*seeding_var)[4]=.5;
	//			(*seeding_var)[5]=.1;
	//			(*seeding_var)[6]=.1;
	//			(*seeding_var)[7]=.1;
	//			for(int i =0; i< mcmc_mod_struct->ppE_Nmod;i++){
	//				(*seeding_var)[8 + i]=2;
	//			}
	//		}
	//	}

	//}
	//else{
	//	std::cout<<
	//		"Input parameters not valid, please check that input is compatible with the supported methods - dimension combinations"<<std::endl;
	//	exit(1);
	//}
}


void MCMC_fisher_transformations(
	double *param, 
	double **fisher, 
	int dimension,
	std::string generation_method,
	bool intrinsic,
	mcmc_data_interface *interface, 
	MCMC_modification_struct *mod_struct,
	void *parameters)
{
	if(!intrinsic){
		fisher[0][0] += 1./(4*M_PI*M_PI);//RA
		fisher[1][1] += 1./4;//sin DEC
		fisher[2][2] += 1./(4*M_PI*M_PI);//psi
		fisher[3][3] += 1./(4);//cos iota
		fisher[4][4] += 1./(4*M_PI*M_PI);//phiref
		fisher[5][5] += 1./(.01);//tc
		fisher[8][8] += 1./.25;//eta
		fisher[9][9] += 1./4;//spin1
		fisher[10][10] += 1./4;//spin2
		if(generation_method.find("PhenomPv2") != std::string::npos){
			fisher[11][11] += 1./4;//cos theta1
			fisher[12][12] += 1./4;//cos theta2
			fisher[13][13] += 1./(4*M_PI*M_PI);//phi1
			fisher[14][14] += 1./(4*M_PI*M_PI);//phi2
		}
	}
	else{
		if(generation_method.find("PhenomPv2") != std::string::npos){
			fisher[1][1] =1./(.25) ;//eta
			fisher[2][2] =1./(4);//spin1
			fisher[3][3] =1./(4);//spin2
			fisher[4][4] =1./(4);//cos theta1
			fisher[5][5] =1./(4);//cos theta2
			fisher[6][6] =1./(4*M_PI*M_PI) ;//phi1
			fisher[7][7] =1./(4*M_PI*M_PI) ;//phi2
		}
		else if (generation_method.find("PhenomD")!=std::string::npos){
			fisher[1][1] =1./(.25) ;//eta
			fisher[2][2] =1./(4) ;//spin1
			fisher[3][3] =1./(4) ;//spin2
	
		}
	}

	if(generation_method.find("dCS") != std::string::npos ||
		generation_method.find("EdGB") != std::string::npos){
		int base = dimension - mod_struct->ppE_Nmod;
		for(int i = 0 ; i<dimension; i++){
			//Transform to root alpha from alpha^2
			double factor = 4* pow(param[base], 3./4.);
			//Transform to km from sec
			factor *= 1000 / c ;
			fisher[base][i] *= factor ;
			fisher[i][base] *= factor;
		}
	}/*
	if(generation_method.find("EA") != std::string::npos){
	  //for(int i = 0 ; i <4; i++){
	    for(int i = 0 ; i <3; i++){
	      for(int j = 0 ; j<dimension; j++){
		if(i!=j){
		  fisher[dimension-1-i][dimension-1-j] = 0;
		  fisher[dimension-1-j][dimension-1-i] = 0;
		}
		/*
		  if(i==j){
		  fisher[dimension-1-i][dimension-1-j] = 1./pow(10, -6.);
		  }
	      }
	    }
	    fisher[dimension-3][dimension-3] = 1./pow(10, -2.);
	    fisher[dimension-2][dimension-2] = 1./pow(10, -4.);
	    fisher[dimension-1][dimension-1] = 1./pow(10, -2.);
	}*/
	/*if(isnan(fabs(fisher[8][8]))){
	  for(int i = 0 ; i<dimension; i++){
	    std::cout<<i<<" ";
	    for(int j = 0 ; j<dimension; j++){
	      std::cout<<std::setprecision(5)<<fisher[i][j]<<" ";
	    }
	    std::cout<<std::endl;
	  }
	  }*/
	//std::cout<<"eta nan?:"<<fisher[8][8]<<std::endl; 
	/*if(isnan(fabs(fisher[8][8])))
	  {
	    std::cout<<"eta nan problem"<<std::endl; 
	    }*/
	return;

}


void MCMC_fisher_wrapper(double *param,  double **output, mcmc_data_interface *interface,void *parameters)
{
	MCMC_user_param *user_param = (MCMC_user_param *)parameters;
	int dimension = interface->max_dim;
	double *temp_params = new double[dimension];
	//#########################################################################
	gen_params_base<double> params;
	std::string local_gen = MCMC_prep_params(param, 
		temp_params,&params, dimension, mcmc_generation_method,mcmc_mod_struct);
	//#########################################################################
	//#########################################################################
	repack_parameters(temp_params, &params, 
		"MCMC_"+mcmc_generation_method, dimension, NULL);
	//std::cout<<temp_params[11]<<" "<<temp_params[12]<<std::endl;
	//repack_parameters(mcmc_init_pos, &params, 
	//	"MCMC_"+mcmc_generation_method, dimension, NULL);
	//#########################################################################
	//#########################################################################
	//std::cout<<"INCL angle fisher: "<<params.incl_angle<<std::endl;
	for(int j =0; j<dimension; j++){
		for(int k =0; k<dimension; k++)
		{
			output[j][k] =0;
		}
	} 
	double **local_freq=mcmc_frequencies;
	double **local_noise=mcmc_noise;
	double **local_weights=user_param->weights;
	int *local_lengths= mcmc_data_length;
	std::string local_integration_method = "SIMPSONS";
	bool local_log10F = user_param->fisher_log10F;
	if(user_param->fisher_freq){
		local_freq = user_param->fisher_freq;
	}
	if(user_param->fisher_PSD){
		local_noise = user_param->fisher_PSD;
	}
	if(user_param->fisher_length){
		local_lengths = user_param->fisher_length;
	}
	if(user_param->fisher_weights){
		local_weights = user_param->fisher_weights;
	}
	if(user_param->fisher_GAUSS_QUAD){
		local_integration_method = "GAUSSLEG";
	}
	std::string local_gen_method = mcmc_generation_method;
	int local_dimension = dimension;  
	/*if(local_gen_method.find("EA") != std::string::npos)
	  {
	    local_gen_method = "IMRPhenomD_NRT";
	    local_dimension -= 3;
	    }*/
	double **temp_out = allocate_2D_array(local_dimension,local_dimension);
	for(int i =0 ; i <mcmc_num_detectors; i++){
		
		//Use AD 
		if(user_param->fisher_AD)
		{	
			std::unique_lock<std::mutex> lock{*(user_param->mFish)};
			//fisher_autodiff(mcmc_frequencies[i], mcmc_data_length[i],
			//	"MCMC_"+mcmc_generation_method, mcmc_detectors[i],mcmc_detectors[0],temp_out,dimension, 
			//	(gen_params *)(&params),  "SIMPSONS",(double *)NULL,false,mcmc_noise[i]);
			fisher_autodiff(local_freq[i], local_lengths[i],
				"MCMC_"+local_gen_method, mcmc_detectors[i],mcmc_detectors[0],temp_out,local_dimension, 
				(gen_params *)(&params),  local_integration_method,local_weights[i],true,local_noise[i]);
		}
		else{
			fisher_numerical(local_freq[i], local_lengths[i],
				"MCMC_"+local_gen_method, mcmc_detectors[i],mcmc_detectors[0],temp_out,local_dimension, 
				&params, mcmc_deriv_order, NULL, NULL, local_noise[i]);

		}
		for(int j =0; j<local_dimension; j++){
			for(int k =0; k<local_dimension; k++)
			{
				output[j][k] +=temp_out[j][k];
				//if(std::isnan(output[j][k]))
				//{
				//      std::cout<<j<<" "<<k<<" "<<temp_out[j][k]<<std::endl;
				//}
			}
		} 
	}
	//Add prior information to fisher
	//if(mcmc_generation_method.find("Pv2") && !mcmc_intrinsic){
	
	MCMC_fisher_transformations(temp_params, output,dimension,local_gen,mcmc_intrinsic,
		interface,mcmc_mod_struct, parameters);
	
	deallocate_2D_array(temp_out, local_dimension,local_dimension);
	//////////////////////////////////////////////
	//if(!interface->burn_phase)
	//{
	//	debugger_print(__FILE__,__LINE__,"Fisher MCMC");
	//	double **cov = allocate_2D_array( dimension,dimension);
	//	gsl_cholesky_matrix_invert(output, cov, dimension);
	//	for(int i = 0 ; i<dimension; i++){
	//		std::cout<<sqrt(cov[i][i])<<" ";	
	//		//for(int j = 0 ; j<dimension; j++){
	//		//	std::cout<<cov[i][j]<<" ";	
	//		//}
	//		//std::cout<<std::endl;	
	//		
	//	}
	//	std::cout<<std::endl;	
	//	deallocate_2D_array(cov, dimension,dimension);
	//}
	//////////////////////////////////////////////

	//Cleanup
	delete [] temp_params;
	if(check_mod(local_gen)){
		//if(local_gen.find("ppE") != std::string::npos ||
		//	local_gen.find("dCS")!=std::string::npos||
		//	local_gen.find("EdGB")!=std::string::npos){
		//	delete [] params.betappe;
		//}
		if(local_gen.find("ppE") != std::string::npos ||
			check_theory_support(local_gen)){
			delete [] params.betappe;
		}
		else if( local_gen.find("gIMR") != std::string::npos){
			if(mcmc_mod_struct ->gIMR_Nmod_phi !=0){
				delete [] params.delta_phi;
			}
			if(mcmc_mod_struct ->gIMR_Nmod_sigma !=0){
				delete [] params.delta_sigma;
			}
			if(mcmc_mod_struct ->gIMR_Nmod_beta !=0){
				delete [] params.delta_beta;
			}
			if(mcmc_mod_struct ->gIMR_Nmod_alpha !=0){
				delete [] params.delta_alpha;
			}

		}
	}


}


//RA, DEC, and PSI were absorbed into gen_params structure -- remove from arguments
double MCMC_likelihood_extrinsic(bool save_waveform, 
	gen_params_base<double> *parameters,
	std::string generation_method, 
	int *data_length, 
	double **frequencies, 
	std::complex<double> **data, 
	double **psd, 
	double **weights, 
	std::string integration_method, 
	bool log10F, 
	std::string *detectors, 
	int num_detectors
	)
{
	//################################################################
	//Outdated -- remove after some time 
	//
	//double *phi = new double[num_detectors];
	//double *theta = new double[num_detectors];
	////celestial_horizon_transform(RA,DEC, gps_time, detectors[0], &phi[0], &theta[0]);
	//double tc_ref, phic_ref, ll=0, delta_t;
	//double LISA_alpha0,LISA_phi0, LISA_thetal, LISA_phil;
	//double *times=NULL;
	////#######################################################################3
	//
	//std::complex<double> *response = 
	//	(std::complex<double> *)malloc(sizeof(std::complex<double>)*
	//		data_length[0]);
	//double T = 1./( frequencies[1]-frequencies[0]);
	//tc_ref = T-parameters->tc;
	//double tc = tc_ref;
	//tc*=2.*M_PI;
	//parameters->tc=0;
	//waveform_polarizations<double> wp;
	//assign_polarizations(generation_method, &wp);
	//wp.allocate_memory(data_length[0]);
	//fourier_waveform(frequencies[0], data_length[0], 
	//	&wp,generation_method, parameters);
	//fourier_detector_response_equatorial(frequencies[0], data_length[0], 
	//	&wp, response, parameters->RA, parameters->DEC, parameters->psi,
	//	parameters->gmst,times, LISA_alpha0, LISA_phi0, LISA_thetal, LISA_phil,detectors[0]);
	//for(int i = 0 ; i<data_length[0];i++){
	//	response[i]*=exp(std::complex<double>(0,tc*(frequencies[0][i])));
	//}	
	////Referecne detector first

	//ll += Log_Likelihood_internal(data[0], 
	//		psd[0],
	//		frequencies[0],
	//		weights[0],
	//		response,
	//		(size_t) data_length[0],
	//		log10F,
	//		integration_method
	//		);
	//for(int i=1; i < num_detectors; i++){
	//	delta_t = DTOA_DETECTOR(parameters->RA, parameters->DEC,parameters->gmst, detectors[0], detectors[i]);
	//	tc = tc_ref - delta_t;
	//	fourier_detector_response_equatorial(frequencies[i], data_length[i], 
	//		&wp, response, parameters->RA, parameters->DEC, 
	//		parameters->psi,parameters->gmst,times, LISA_alpha0, 
	//		LISA_phi0, LISA_thetal, LISA_phil,detectors[i]);
	//	tc*=2.*M_PI;
	//	for(int j = 0 ; j<data_length[i];j++){
	//		response[j]*=exp(std::complex<double>(
	//			0,tc*(frequencies[i][j])
	//			));
	//	}	
	//
	//	ll += Log_Likelihood_internal(data[i], 
	//		psd[i],
	//		frequencies[i],
	//		weights[i],
	//		response,
	//		(size_t) data_length[i],
	//		log10F,
	//		integration_method	
	//		);
	//}
	//wp.deallocate_memory();
	//free( response);
	//
	//delete [] phi;
	//delete [] theta;
	
	//debugger_print(__FILE__,__LINE__,ll);
	//################################################################
	//################################################################
	
	
	//#################################################################################
	//#################################################################################
	double T = 1./( frequencies[1]-frequencies[0]);
	double tc_ref = T-parameters->tc;
	double ll=0;
	std::complex<double> **responses = new std::complex<double>*[num_detectors];	
	for(int i = 0 ; i<num_detectors; i++){
		responses[i] = new std::complex<double>[data_length[i]];
	}
	parameters->tc = tc_ref;	
	create_coherent_GW_detection(detectors, num_detectors, frequencies,data_length, save_waveform, parameters, generation_method, responses);
	for(int i = 0 ;i<num_detectors; i++){
		ll += Log_Likelihood_internal(data[i],psd[i],frequencies[i],weights[i],responses[i],data_length[i], log10F,integration_method);	
	}		
	
	for(int i = 0 ; i<num_detectors; i++){
		delete [] responses[i];
	}
	delete [] responses;
	//#################################################################################
	//#################################################################################
	
	return ll;
}
/*! \brief utility to do MCMC specific transformations on the input param vector before passing to the repacking utillity
 *
 * Returns the local generation method to be used in the LL functions
 */
std::string MCMC_prep_params(double *param, double *temp_params, gen_params_base<double> *gen_params, int dimension, std::string generation_method, MCMC_modification_struct *mod_struct)
{
	if(mcmc_intrinsic) gen_params->sky_average = true;
	else gen_params->sky_average = false;
	gen_params->tidal_love = mod_struct->tidal_love;
	gen_params->tidal_love_error = mod_struct->tidal_love_error;
	gen_params->alpha_param = mod_struct->alpha_param;
	gen_params->EA_region1 = mod_struct->EA_region1;
	gen_params->EOS_plat_flag = mod_struct->EOS_plat_flag; 
	gen_params->f_ref = 20;
	gen_params->shift_time = true;
	gen_params->shift_phase = true;
	//gen_params->shift_time = false;
	//gen_params->shift_phase = false;
	gen_params->gmst = mcmc_gmst;
	gen_params->equatorial_orientation=false;
	gen_params->horizon_coord=false;

	gen_params->NSflag1 = mod_struct->NSflag1;
	gen_params->NSflag2 = mod_struct->NSflag2;
	//gen_params->NSflag1 = false;
	//gen_params->NSflag2 = false;
	for(int i = 0 ; i <dimension; i++){
		temp_params[i]=param[i];
	}
	int base = dimension;
	if(check_mod(generation_method)){
		//if(generation_method.find("ppE") != std::string::npos ||
		//	generation_method.find("dCS") !=std::string::npos||
		//	generation_method.find("EdGB") != std::string::npos){
		//	gen_params->bppe=mcmc_mod_struct->bppe;
		//	gen_params->Nmod=mcmc_mod_struct->ppE_Nmod;
		//	gen_params->betappe=new double[gen_params->Nmod];
		//}
		if(generation_method.find("ppE") != std::string::npos ||
			check_theory_support(generation_method)){
			gen_params->bppe=mod_struct->bppe;
			gen_params->Nmod=mod_struct->ppE_Nmod;
			gen_params->betappe=new double[gen_params->Nmod];
			base = dimension - mod_struct->ppE_Nmod;
			//debugger_print(__FILE__,__LINE__,mod_struct->ppE_Nmod);
		}
		else if(generation_method.find("gIMR") != std::string::npos){
			gen_params->Nmod_phi=mod_struct->gIMR_Nmod_phi;
			gen_params->phii=mod_struct->gIMR_phii;
			if(gen_params->Nmod_phi !=0){
				gen_params->delta_phi=new double[gen_params->Nmod_phi];
			}
			gen_params->Nmod_sigma=mod_struct->gIMR_Nmod_sigma;
			gen_params->sigmai=mod_struct->gIMR_sigmai;
			if(gen_params->Nmod_sigma !=0){
				gen_params->delta_sigma=new double[gen_params->Nmod_sigma];
			}
			gen_params->Nmod_beta=mod_struct->gIMR_Nmod_beta;
			gen_params->betai=mod_struct->gIMR_betai;
			if(gen_params->Nmod_beta !=0){
				gen_params->delta_beta=new double[gen_params->Nmod_beta];
			}
			gen_params->Nmod_alpha=mod_struct->gIMR_Nmod_alpha;
			gen_params->alphai=mod_struct->gIMR_alphai;
			if(gen_params->Nmod_alpha !=0){
				gen_params->delta_alpha=new double[gen_params->Nmod_alpha];
			}
			base = dimension 
				- mod_struct->gIMR_Nmod_phi 
				- mod_struct->gIMR_Nmod_sigma 
				- mod_struct->gIMR_Nmod_beta 
				- mod_struct->gIMR_Nmod_alpha; 
		}
		if((generation_method.find("dCS")!= std::string::npos ||
			generation_method.find("EdGB")!=std::string::npos)){
			//temp_params[base] = pow(temp_params[base],.25)/(c*1000);
			temp_params[base] = 
				pow_int(temp_params[base]/(c/1000.) , 4);
		}
	}
	return generation_method;
}
double MCMC_likelihood_wrapper(double *param, mcmc_data_interface *interface ,void *parameters)
{
  //double start = omp_get_wtime();
  //return 2;
  MCMC_user_param *user_param = (MCMC_user_param *)parameters;

  int dimension = interface->max_dim;
  double ll = 0;
  double *temp_params = new double[dimension];
  //#########################################################################
  gen_params_base<double> gen_params;
  std::string local_gen = MCMC_prep_params(param, 
					   temp_params,&gen_params, dimension, mcmc_generation_method,mcmc_mod_struct);
  //#########################################################################
  //#########################################################################

  //repack_non_parameters(temp_params, &gen_params, 
  //"MCMC_"+mcmc_generation_method, dimension, NULL);
  repack_parameters(temp_params, &gen_params, 
		    "MCMC_"+mcmc_generation_method, dimension, NULL);
  //#########################################################################
  //#########################################################################
  //return 1;
  std::complex<double> **local_data = mcmc_data;
  double **local_freqs = mcmc_frequencies;
  double **local_noise = mcmc_noise;
  double **local_weights = user_param->weights;
  int *local_lengths = mcmc_data_length;
  fftw_outline *local_plans = mcmc_fftw_plans;
  std::string local_integration_method="SIMPSONS";
  //if(interface->burn_phase && user_param->burn_data){
  if(false){
    local_data = user_param->burn_data;
    local_freqs = user_param->burn_freqs;
    local_noise = user_param->burn_noise;
    local_lengths = user_param->burn_lengths;
    local_plans = user_param->burn_plans;
  }
  if(user_param->GAUSS_QUAD){
    local_integration_method="GAUSSLEG";
  }
  if(mcmc_intrinsic){
    if(mcmc_generation_method.find("IMRPhenomD") != std::string::npos){
      if(!mcmc_save_waveform){
	for(int i=0; i < mcmc_num_detectors; i++){
	  gen_params.theta=0;	
	  gen_params.phi=0;	
	  gen_params.psi=0;	
	  gen_params.phiRef = 1;
	  gen_params.f_ref = 10;
	  gen_params.incl_angle=0;	
	  gen_params.tc =1;
	  std::complex<double> *response =
	    (std::complex<double> *) malloc(sizeof(std::complex<double>) * local_lengths[i]);
	  fourier_detector_response_horizon(local_freqs[i], local_lengths[i], response, mcmc_detectors[i], local_gen, &gen_params);
	  ll += maximized_Log_Likelihood_aligned_spin_internal(local_data[i], 
							       local_noise[i],
							       local_freqs[i],
							       response,
							       (size_t) local_lengths[i],
							       &local_plans[i]
							       );
	  //ll += maximized_Log_Likelihood(mcmc_data[i], 
	  //		mcmc_noise[i],
	  //		mcmc_frequencies[i],
	  //		(size_t) mcmc_data_length[i],
	  //		&gen_params,
	  //		mcmc_detectors[i],
	  //		local_gen,
	  //		&mcmc_fftw_plans[i]
	  //		);
	  free(response);
				}
      }
      else{
	gen_params.theta=0;	
	gen_params.phi=0;	
	gen_params.psi=0;	
	gen_params.phiRef = 1;
	gen_params.f_ref = 10;
	gen_params.incl_angle=0;	
	gen_params.tc =1;
	std::complex<double> *response =
	  (std::complex<double> *) malloc(sizeof(std::complex<double>) * local_lengths[0]);
	fourier_detector_response_horizon(local_freqs[0], local_lengths[0], response, mcmc_detectors[0], local_gen, &gen_params);
	//std::complex<double> *hc =
	//	(std::complex<double> *) malloc(sizeof(std::complex<double>) * mcmc_data_length[0]);
	//std::complex<double> *hp =
	//	(std::complex<double> *) malloc(sizeof(std::complex<double>) * mcmc_data_length[0]);
	//fourier_waveform(mcmc_frequencies[0], mcmc_data_length[0], hp,hc, local_gen, &gen_params);
	for(int i=0; i < mcmc_num_detectors; i++){
	  ll += maximized_Log_Likelihood_aligned_spin_internal(local_data[i], 
							       local_noise[i],
							       local_freqs[i],
							       response,
							       (size_t) local_lengths[i],
							       &local_plans[i]
							       );
	  //ll += maximized_Log_Likelihood(mcmc_data[i], 
	  //		mcmc_noise[i],
	  //		mcmc_frequencies[i],
	  //		(size_t) mcmc_data_length[i],
	  //		&gen_params,
	  //		mcmc_detectors[i],
	  //		local_gen,
	  //		&mcmc_fftw_plans[i]
	  //		);
	  //ll += maximized_Log_Likelihood_unaligned_spin_internal(mcmc_data[i], 
	  //		mcmc_noise[i],
	  //		mcmc_frequencies[i],
	  //		hp,
	  //		hc,
	  //		(size_t) mcmc_data_length[i],
	  //		&mcmc_fftw_plans[i]
	  //		);
	}
	free(response);
	//free(hp);
	//free(hc);
	
      }

    }
    else if(mcmc_generation_method.find("IMRPhenomP")!=std::string::npos){
      //if(!mcmc_save_waveform){
      if(false){
      }
      else{
	gen_params.theta=0;	
	gen_params.phi=0;	
	gen_params.psi=0;	
	gen_params.phiRef = 1;
	gen_params.f_ref = 20;
	gen_params.incl_angle=0;	
	gen_params.tc =1;
	waveform_polarizations<double> wp;
	assign_polarizations(mcmc_generation_method,&wp);
	wp.allocate_memory(local_lengths[0]);
	fourier_waveform(local_freqs[0],local_lengths[0], &wp, local_gen, &gen_params);
	for(int i=0; i < mcmc_num_detectors; i++){
	  ll += maximized_Log_Likelihood_unaligned_spin_internal(local_data[i], 
								 local_noise[i],
								 local_freqs[i],
								 wp.hplus,
								 wp.hcross,
								 (size_t) local_lengths[i],
								 &local_plans[i]
								 );
	}
	wp.deallocate_memory();
      }
      
    }
  }
  else{
    double RA = gen_params.RA;
    double DEC = gen_params.DEC;
    double PSI = gen_params.psi;
    //if(mcmc_generation_method.find("IMRPhenomD") != std::string:npos){
    
    ll =  MCMC_likelihood_extrinsic(mcmc_save_waveform, 
				    &gen_params,local_gen, local_lengths, 
				    local_freqs, local_data, local_noise, local_weights, local_integration_method, user_param->log10F,mcmc_detectors, 
				    mcmc_num_detectors);
    //ll=2;
    
    //ll = Log_Likelihood(mcmc_data[0], 
    //		mcmc_noise[0],
    //		mcmc_frequencies[0],
    //		mcmc_data_length[0],
    //		&gen_params,
    //		mcmc_detectors[0],
    //		local_gen,
    //		&mcmc_fftw_plans[0]
    //		);
    
    //}
    //else if(mcmc_generation_method.find("IMRPhenomP")!=std::string::npos){
    
    //}
  }
  //Cleanup
  delete [] temp_params;
  if(check_mod(local_gen)){
    //if( local_gen.find("ppE") != std::string::npos ||
    //	local_gen.find("dCS") != std::string::npos ||
    //	local_gen.find("EdGB") != std::string::npos){
    //	delete [] gen_params.betappe;
    //}
    if( local_gen.find("ppE") != std::string::npos ||
	check_theory_support(local_gen)){
      delete [] gen_params.betappe;
    }
    else if( local_gen.find("gIMR") != std::string::npos){
      if(mcmc_mod_struct ->gIMR_Nmod_phi !=0){
	delete [] gen_params.delta_phi;
      }
      if(mcmc_mod_struct ->gIMR_Nmod_sigma !=0){
	delete [] gen_params.delta_sigma;
      }
      if(mcmc_mod_struct ->gIMR_Nmod_beta !=0){
	delete [] gen_params.delta_beta;
      }
      if(mcmc_mod_struct ->gIMR_Nmod_alpha !=0){
	delete [] gen_params.delta_alpha;
      }
      
    }
  }
  //std::cout<<"LL time for eval: "<<(double)(omp_get_wtime() -start)<<std::endl;

  //std::cout<<"Likelihood: "<<ll<<std::endl;
  if(isnan(ll)){
    std::cout<<"NAN"<<std::endl;
    for(int i = 0 ; i<dimension; i++){
      std::cout<<param[i]<<", ";
    }
    std::cout<<std::endl;
  }
  return ll;
  
}
//######################################################################################
//######################################################################################
/*! \brief Unpacks MCMC parameters for method specific initiation (RJ version)
 *
 * Populates seeding vector if non supplied, populates mcmc_Nmod, populates mcmc_log_beta, populates mcmc_intrinsic
 */
void RJPTMCMC_method_specific_prep(std::string generation_method, int max_dim, int min_dim,double *seeding_var, bool local_seeding)
{
	if(min_dim==11 && (generation_method.find("IMRPhenomD") != std::string::npos)){
		mcmc_intrinsic=false;
		std::cout<<"Sampling in parameters: cos inclination, RA, DEC, ln DL ,ln chirpmass, eta, chi1, chi2, psi";
		for(int i =0; i<mcmc_Nmod_max; i++){
			std::cout<<", beta"<<i<<" ("<<mcmc_bppe[i]<<")";
		}
		std::cout<<endl;
		mcmc_log_beta = false;
		if(local_seeding){
			seeding_var = new double[max_dim];
			seeding_var[0]=.1;
			seeding_var[1]=.5;
			seeding_var[2]=.1;
			seeding_var[3]=.1;
			seeding_var[4]=.1;
			seeding_var[5]=.1;
			seeding_var[6]=.1;
			seeding_var[7]=.1;
			seeding_var[8]=.1;
			seeding_var[9]=.1;
			seeding_var[10]=.1;
			for(int i =0; i< mcmc_Nmod_max;i++){
				seeding_var[11 + i]=2;
			}
		}
	}
	else if(min_dim==4 && (generation_method.find("IMRPhenomD") != std::string::npos)){
		mcmc_intrinsic=true;
		std::cout<<"Sampling in parameters: cos inclination, RA, DEC, ln DL ,ln chirpmass, eta, chi1, chi2, psi";
		for(int i =0; i<mcmc_Nmod_max; i++){
			std::cout<<", beta"<<i<<" ("<<mcmc_bppe[i]<<")";
		}
		std::cout<<endl;
		mcmc_log_beta = false;
		if(local_seeding){
			seeding_var = new double[max_dim];
			seeding_var[0]=.1;
			seeding_var[1]=.5;
			seeding_var[2]=.1;
			seeding_var[3]=.1;
			for(int i =0; i< mcmc_Nmod_max;i++){
				seeding_var[4 + i]=2;
			}
		}
	}
	else if(min_dim==15 && (generation_method.find("IMRPhenomPv2") != std::string::npos)){
		mcmc_intrinsic=false;
		std::cout<<"Sampling in parameters: cos inclination, RA, DEC, ln DL ,ln chirpmass, eta, chi1, chi2, psi";
		for(int i =0; i<mcmc_Nmod_max; i++){
			std::cout<<", beta"<<i<<" ("<<mcmc_bppe[i]<<")";
		}
		std::cout<<endl;
		mcmc_log_beta = false;
		if(local_seeding){
			seeding_var = new double[max_dim];
			seeding_var[0]=.1;
			seeding_var[1]=.5;
			seeding_var[2]=.1;
			seeding_var[3]=.1;
			seeding_var[4]=.1;
			seeding_var[5]=.1;
			seeding_var[6]=.1;
			seeding_var[7]=.1;
			seeding_var[8]=.1;
			seeding_var[9]=.1;
			seeding_var[10]=.1;
			seeding_var[11]=.1;
			seeding_var[12]=.1;
			seeding_var[13]=.1;
			seeding_var[14]=.1;
			for(int i =0; i< mcmc_Nmod_max;i++){
				seeding_var[15 + i]=2;
			}
		}
	}
	else if(min_dim==8 && (generation_method.find("IMRPhenomPv2") != std::string::npos)){
		mcmc_intrinsic=true;
		std::cout<<"Sampling in parameters: cos inclination, RA, DEC, ln DL ,ln chirpmass, eta, chi1, chi2, psi";
		for(int i =0; i<mcmc_Nmod_max; i++){
			std::cout<<", beta"<<i<<" ("<<mcmc_bppe[i]<<")";
		}
		std::cout<<endl;
		mcmc_log_beta = false;
		if(local_seeding){
			seeding_var = new double[max_dim];
			seeding_var[0]=.1;
			seeding_var[1]=.5;
			seeding_var[2]=.1;
			seeding_var[3]=.1;
			seeding_var[4]=.1;
			seeding_var[5]=.1;
			seeding_var[6]=.1;
			seeding_var[7]=.1;
			for(int i =0; i< mcmc_Nmod_max;i++){
				seeding_var[8 + i]=2;
			}
		}
	}

}

/*! \brief Performs RJPTMCMC on two competeting models, waveform 1 and waveform 2. Assumes waveform 2 is a superset including waveform 1 (ie dCS Pv2 and Pv2 or NRT Pv2 and Pv2, but not dCS Pv2 and D).
 *
 */
void continue_RJPTMCMC_MH_dynamic_PT_alloc_comprehensive_2WF_GW(
	std::string checkpoint_file_start,
	mcmc_sampler_output *sampler_output,
	double **output,
	int **status,
	int *model_status,
	int nested_model_number,
	int N_steps,
	int max_chain_N_thermo_ensemble,
	double **prior_ranges,/**< Range of priors on MODIFICATIONS -- shape : double[N_mods][2] with low in 0 and high in 1**/
	double *chain_temps,
	int swp_freq,
	int t0,
	int nu,
	int max_chunk_size,
	std::string chain_distribution_scheme,
	double(*log_prior)(double *param, int *status,int model_status,mcmc_data_interface *interface,void *parameters),
	int numThreads,
	bool pool,
	bool show_prog,
	int num_detectors,
	std::complex<double> **data,
	double **noise_psd,
	double **frequencies,
	int *data_length,
	double gps_time,
	std::string *detectors,
	MCMC_modification_struct *mod_struct,
	std::string generation_method_base,
	std::string generation_method_extended,
	std::string statistics_filename,
	std::string chain_filename,
	std::string likelihood_log_filename,
	std::string checkpoint_filename
	)
{
	int max_dimension, min_dimension,chain_N;
	dimension_from_checkpoint_file(checkpoint_file_start, &min_dimension, &max_dimension);
	chain_number_from_checkpoint_file(checkpoint_file_start, &chain_N) ;
	sampler_output->RJ = true;

	std::mutex fisher_mutex;

	bool update_RJ_widths = false;
	//Create fftw plan for each detector (length of data stream may be different)
	fftw_outline *plans= (fftw_outline *)malloc(sizeof(fftw_outline)*num_detectors);
	for (int i =0;i<num_detectors;i++)
	{	
		allocate_FFTW_mem_forward(&plans[i] , data_length[i]);
	}
	mcmc_noise = noise_psd;	
	mcmc_frequencies = frequencies;
	mcmc_data = data;
	mcmc_data_length = data_length;
	mcmc_detectors = detectors;
	mcmc_generation_method_base = generation_method_base;
	mcmc_generation_method_extended = generation_method_extended;
	mcmc_fftw_plans = plans;
	mcmc_num_detectors = num_detectors;
	mcmc_gps_time = gps_time;
	mcmc_gmst = gps_to_GMST_radian(mcmc_gps_time);
	mcmc_mod_struct = mod_struct;
	mcmc_log_beta = false;
	mcmc_intrinsic = false;

	//To save time, intrinsic waveforms can be saved between detectors, if the 
	//frequencies are all the same
	mcmc_save_waveform = true;
	for(int i =1 ;i<mcmc_num_detectors; i++){
		if( mcmc_data_length[i] != mcmc_data_length[0] ||
			mcmc_frequencies[i][0]!= mcmc_frequencies[0][0] ||
			mcmc_frequencies[i][mcmc_data_length[i]-1] 
				!= mcmc_frequencies[0][mcmc_data_length[0]-1]){
			mcmc_save_waveform= false;
		}
			
	}

	bool local_seeding = false;
	double *seeding_var = NULL;

	RJPTMCMC_method_specific_prep(generation_method_extended, max_dimension, min_dimension,seeding_var, local_seeding);


	//######################################################
	int T = (int)(1./(mcmc_frequencies[0][1]-mcmc_frequencies[0][0]));
	debugger_print(__FILE__,__LINE__,T);
	int burn_factor = T/4; //Take all sources to 4 seconds
	debugger_print(__FILE__,__LINE__,burn_factor);
	std::complex<double> **burn_data = new std::complex<double>*[mcmc_num_detectors];
	double **burn_freqs = new double*[mcmc_num_detectors];
	double **burn_noise = new double*[mcmc_num_detectors];
	int *burn_lengths = new int[mcmc_num_detectors];
	fftw_outline *burn_plans= new fftw_outline[num_detectors];
	for(int j = 0; j<mcmc_num_detectors; j++){
		burn_lengths[j] = mcmc_data_length[j]/burn_factor;
		burn_data[j]= new std::complex<double>[burn_lengths[j]];
		burn_freqs[j]= new double[burn_lengths[j]];
		burn_noise[j]= new double[burn_lengths[j]];
		allocate_FFTW_mem_forward(&burn_plans[j], burn_lengths[j]);
		int ct = 0;
		for( int k = 0 ; k<mcmc_data_length[j]; k++){
			if(k%burn_factor==0 && ct<burn_lengths[j]){
				burn_data[j][ct] = mcmc_data[j][k];
				burn_freqs[j][ct] = mcmc_frequencies[j][k];
				burn_noise[j][ct] = mcmc_noise[j][k];
				ct++;
			}
		}
	}
	

	//######################################################
	//######################################################
	//Fishers sometimes need AD, but that's slow and single 
	//threaded -- use GAUSS quad with precomputed weights 
	//and abscissa 
	//double *fish_freqs = NULL;
	//double *fish_weights = NULL;
	//double **fish_psd = NULL; //Needs to be interpolated from data given
	//int fish_length = 100;
	//double flow = mcmc_frequencies[0][0];
	//double fhigh = mcmc_frequencies[0][mcmc_data_length[0]-1];
	//if(mod_struct->GAUSS_QUAD){
	//	fish_freqs = new double[fish_length];	
	//	fish_weights = new double[fish_length];	
	//	fish_psd = new double*[mcmc_num_detectors];	
	//	
	//	gauleg(log10(flow), log10(fhigh), fish_freqs, fish_weights, fish_length);
	//	for(int i = 0 ; i<fish_length; i++){
	//		fish_freqs[i]=pow(10,fish_freqs[i]);
	//	}
	//	for(int i = 0 ; i<mcmc_num_detectors; i++){
	//		fish_psd[i] = new double[fish_length];
	//		gsl_interp_accel *accel = gsl_interp_accel_alloc();
	//		gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear, mcmc_data_length[i]);
	//		gsl_spline_init(spline, mcmc_frequencies[i], mcmc_noise[i],mcmc_data_length[i]);
	//		for(int j = 0 ; j<fish_length; j++){
	//			fish_psd[i][j] = gsl_spline_eval(spline, fish_freqs[j],accel);
	//		}
	//		gsl_spline_free(spline);
	//	}
	//}
	//######################################################
	//######################################################
	const gsl_rng_type *gsl_T;
	gsl_rng **rvec = new gsl_rng*[chain_N];
	gsl_rng_env_setup();	
	gsl_T = gsl_rng_default;
	for(int i = 0 ; i<chain_N; i++){
		rvec[i] = gsl_rng_alloc(gsl_T);
	}
	
	//######################################################
	MCMC_user_param **user_parameters=NULL;
	user_parameters = new MCMC_user_param*[chain_N];
	for(int i = 0 ;i<chain_N; i++){
		user_parameters[i] = new MCMC_user_param;
		
		user_parameters[i]->burn_data = burn_data;
		user_parameters[i]->burn_freqs = burn_freqs;
		user_parameters[i]->burn_noise = burn_noise;
		user_parameters[i]->burn_lengths = burn_lengths;
		user_parameters[i]->burn_plans=burn_plans;

		user_parameters[i]->mFish= &fisher_mutex;
		user_parameters[i]->GAUSS_QUAD= mod_struct->GAUSS_QUAD;
		user_parameters[i]->log10F = mod_struct->log10F;
		if(mod_struct->GAUSS_QUAD){
			user_parameters[i]->weights = mod_struct->weights;			
		}
		else{
			user_parameters[i]->weights = new double*[num_detectors];			
			for(int j = 0 ; j<num_detectors; j++){
				user_parameters[i]->weights[j]=NULL;
			}
		}

		user_parameters[i]->fisher_GAUSS_QUAD = mod_struct->fisher_GAUSS_QUAD;
		user_parameters[i]->fisher_log10F = mod_struct->fisher_log10F;
		if(mod_struct->fisher_freq && mod_struct->fisher_weights && mod_struct->fisher_PSD){
			user_parameters[i]->fisher_freq= mod_struct->fisher_freq;
			user_parameters[i]->fisher_weights= mod_struct->fisher_weights;
			user_parameters[i]->fisher_PSD= mod_struct->fisher_PSD;
			user_parameters[i]->fisher_length= mod_struct->fisher_length;
		}

		user_parameters[i]->r = rvec[i];

		user_parameters[i]->mod_prior_ranges = prior_ranges;
		user_parameters[i]->mod_struct = mod_struct;

		//user_parameters[i]->burn_freqs = mcmc_frequencies;
		//user_parameters[i]->burn_data = mcmc_data;
		//user_parameters[i]->burn_noise = mcmc_noise;
		//user_parameters[i]->burn_lengths = mcmc_data_length;
	}
	//######################################################


	continue_RJPTMCMC_MH_dynamic_PT_alloc_comprehensive(checkpoint_file_start,sampler_output,output,status, model_status, nested_model_number,
		N_steps, 
		max_chain_N_thermo_ensemble, chain_temps, 
		swp_freq, t0, nu,max_chunk_size,chain_distribution_scheme,
		log_prior,RJMCMC_2WF_likelihood_wrapper, RJMCMC_2WF_fisher_wrapper,RJMCMC_2WF_RJ_proposal_wrapper,(void**)user_parameters,numThreads, pool, 
		show_prog,update_RJ_widths,statistics_filename,
		chain_filename, likelihood_log_filename,checkpoint_filename);
	
	//Deallocate fftw plans
	for (int i =0;i<num_detectors;i++)
		deallocate_FFTW_mem(&plans[i]);
	//#################################################
	for(int i = 0 ; i<num_detectors; i++){
		delete [] burn_data[i];
		delete [] burn_freqs[i];
		delete [] burn_noise[i];
		deallocate_FFTW_mem(&burn_plans[i]);
	}
	//if(mod_struct->GAUSS_QUAD){
	//	delete [] fish_freqs;
	//	delete [] fish_weights;
	//	for(int i = 0 ;i<mcmc_num_detectors; i++){
	//		delete [] fish_psd[i];
	//	}
	//	delete [] fish_psd;
	//}
	delete [] burn_data;
	delete [] burn_lengths;
	delete [] burn_noise;
	delete [] burn_freqs;
	delete [] burn_plans;
	for(int i = 0 ; i<chain_N; i++){
		delete user_parameters[i];
		gsl_rng_free(rvec[i]);
	}
	delete [] rvec;
	delete [] user_parameters;
	//#################################################
	free(plans);
}

/*! \brief Performs RJPTMCMC on two competeting models, waveform 1 and waveform 2. Assumes waveform 2 is a superset including waveform 1 (ie dCS Pv2 and Pv2 or NRT Pv2 and Pv2, but not dCS Pv2 and D).
 *
 */
void RJPTMCMC_MH_dynamic_PT_alloc_comprehensive_2WF_GW(
	mcmc_sampler_output *sampler_output,
	double **output,
	int **status,
	int *model_status,
	int nested_model_number,
	int max_dimension,
	int min_dimension,
	int N_steps,
	int chain_N,
	int max_chain_N_thermo_ensemble,
	double *initial_pos,
	int *initial_status,
	int initial_model_status,
	double *seeding_var,
	double **ensemble_initial_pos,
	int **ensemble_initial_status,
	int *ensemble_initial_model_status,
	double **prior_ranges,/**< Range of priors on MODIFICATIONS -- shape : double[N_mods][2] with low in 0 and high in 1**/
	double *chain_temps,
	int swp_freq,
	int t0,
	int nu,
	int max_chunk_size,
	std::string chain_distribution_scheme,
	double(*log_prior)(double *param, int *status,int model_status,mcmc_data_interface *interface,void *parameters),
	int numThreads,
	bool pool,
	bool show_prog,
	int num_detectors,
	std::complex<double> **data,
	double **noise_psd,
	double **frequencies,
	int *data_length,
	double gps_time,
	std::string *detectors,
	MCMC_modification_struct *mod_struct,
	std::string generation_method_base,
	std::string generation_method_extended,
	std::string statistics_filename,
	std::string chain_filename,
	std::string likelihood_log_filename,
	std::string checkpoint_filename
	)
{
	sampler_output->RJ = true;

	std::mutex fisher_mutex;

	bool update_RJ_widths = false;
	//Create fftw plan for each detector (length of data stream may be different)
	fftw_outline *plans= (fftw_outline *)malloc(sizeof(fftw_outline)*num_detectors);
	for (int i =0;i<num_detectors;i++)
	{	
		allocate_FFTW_mem_forward(&plans[i] , data_length[i]);
	}
	mcmc_noise = noise_psd;	
	mcmc_init_pos = initial_pos;
	mcmc_frequencies = frequencies;
	mcmc_data = data;
	mcmc_data_length = data_length;
	mcmc_detectors = detectors;
	mcmc_generation_method_base = generation_method_base;
	mcmc_generation_method_extended = generation_method_extended;
	mcmc_fftw_plans = plans;
	mcmc_num_detectors = num_detectors;
	mcmc_gps_time = gps_time;
	mcmc_gmst = gps_to_GMST_radian(mcmc_gps_time);
	mcmc_mod_struct = mod_struct;
	mcmc_log_beta = false;
	mcmc_intrinsic = false;

	//To save time, intrinsic waveforms can be saved between detectors, if the 
	//frequencies are all the same
	mcmc_save_waveform = true;
	for(int i =1 ;i<mcmc_num_detectors; i++){
		if( mcmc_data_length[i] != mcmc_data_length[0] ||
			mcmc_frequencies[i][0]!= mcmc_frequencies[0][0] ||
			mcmc_frequencies[i][mcmc_data_length[i]-1] 
				!= mcmc_frequencies[0][mcmc_data_length[0]-1]){
			mcmc_save_waveform= false;
		}
			
	}

	bool local_seeding ;
	if(!seeding_var)
		local_seeding = true;
	else
		local_seeding = false;

	RJPTMCMC_method_specific_prep(generation_method_extended, max_dimension, min_dimension,seeding_var,local_seeding);


	//######################################################
	int T = (int)(1./(mcmc_frequencies[0][1]-mcmc_frequencies[0][0]));
	debugger_print(__FILE__,__LINE__,T);
	int burn_factor = T/4; //Take all sources to 4 seconds
	debugger_print(__FILE__,__LINE__,burn_factor);
	std::complex<double> **burn_data = new std::complex<double>*[mcmc_num_detectors];
	double **burn_freqs = new double*[mcmc_num_detectors];
	double **burn_noise = new double*[mcmc_num_detectors];
	int *burn_lengths = new int[mcmc_num_detectors];
	fftw_outline *burn_plans= new fftw_outline[num_detectors];
	for(int j = 0; j<mcmc_num_detectors; j++){
		burn_lengths[j] = mcmc_data_length[j]/burn_factor;
		burn_data[j]= new std::complex<double>[burn_lengths[j]];
		burn_freqs[j]= new double[burn_lengths[j]];
		burn_noise[j]= new double[burn_lengths[j]];
		allocate_FFTW_mem_forward(&burn_plans[j], burn_lengths[j]);
		int ct = 0;
		for( int k = 0 ; k<mcmc_data_length[j]; k++){
			if(k%burn_factor==0 && ct<burn_lengths[j]){
				burn_data[j][ct] = mcmc_data[j][k];
				burn_freqs[j][ct] = mcmc_frequencies[j][k];
				burn_noise[j][ct] = mcmc_noise[j][k];
				ct++;
			}
		}
	}
	

	//######################################################
	//######################################################
	//Fishers sometimes need AD, but that's slow and single 
	//threaded -- use GAUSS quad with precomputed weights 
	//and abscissa 
	//double *fish_freqs = NULL;
	//double *fish_weights = NULL;
	//double **fish_psd = NULL; //Needs to be interpolated from data given
	//int fish_length = 100;
	//double flow = mcmc_frequencies[0][0];
	//double fhigh = mcmc_frequencies[0][mcmc_data_length[0]-1];
	//if(mod_struct->GAUSS_QUAD){
	//	fish_freqs = new double[fish_length];	
	//	fish_weights = new double[fish_length];	
	//	fish_psd = new double*[mcmc_num_detectors];	
	//	
	//	gauleg(log10(flow), log10(fhigh), fish_freqs, fish_weights, fish_length);
	//	for(int i = 0 ; i<fish_length; i++){
	//		fish_freqs[i]=pow(10,fish_freqs[i]);
	//	}
	//	for(int i = 0 ; i<mcmc_num_detectors; i++){
	//		fish_psd[i] = new double[fish_length];
	//		gsl_interp_accel *accel = gsl_interp_accel_alloc();
	//		gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear, mcmc_data_length[i]);
	//		gsl_spline_init(spline, mcmc_frequencies[i], mcmc_noise[i],mcmc_data_length[i]);
	//		for(int j = 0 ; j<fish_length; j++){
	//			fish_psd[i][j] = gsl_spline_eval(spline, fish_freqs[j],accel);
	//		}
	//		gsl_spline_free(spline);
	//	}
	//}
	//######################################################
	//######################################################
	const gsl_rng_type *gsl_T;
	gsl_rng **rvec = new gsl_rng*[chain_N];
	gsl_rng_env_setup();	
	gsl_T = gsl_rng_default;
	for(int i = 0 ; i<chain_N; i++){
		rvec[i] = gsl_rng_alloc(gsl_T);
	}
	
	//######################################################
	MCMC_user_param **user_parameters=NULL;
	user_parameters = new MCMC_user_param*[chain_N];
	for(int i = 0 ;i<chain_N; i++){
		user_parameters[i] = new MCMC_user_param;
		
		user_parameters[i]->burn_data = burn_data;
		user_parameters[i]->burn_freqs = burn_freqs;
		user_parameters[i]->burn_noise = burn_noise;
		user_parameters[i]->burn_lengths = burn_lengths;
		user_parameters[i]->burn_plans=burn_plans;

		user_parameters[i]->mFish= &fisher_mutex;
		user_parameters[i]->GAUSS_QUAD= mod_struct->GAUSS_QUAD;
		user_parameters[i]->log10F = mod_struct->log10F;
		if(mod_struct->GAUSS_QUAD){
			user_parameters[i]->weights = mod_struct->weights;			
		}
		else{
			user_parameters[i]->weights = new double*[num_detectors];			
			for(int j = 0 ; j<num_detectors; j++){
				user_parameters[i]->weights[j]=NULL;
			}
		}

		user_parameters[i]->fisher_GAUSS_QUAD = mod_struct->fisher_GAUSS_QUAD;
		user_parameters[i]->fisher_log10F = mod_struct->fisher_log10F;
		if(mod_struct->fisher_freq && mod_struct->fisher_weights && mod_struct->fisher_PSD){
			user_parameters[i]->fisher_freq= mod_struct->fisher_freq;
			user_parameters[i]->fisher_weights= mod_struct->fisher_weights;
			user_parameters[i]->fisher_PSD= mod_struct->fisher_PSD;
			user_parameters[i]->fisher_length= mod_struct->fisher_length;
		}


		user_parameters[i]->r = rvec[i];

		user_parameters[i]->mod_prior_ranges = prior_ranges;
		user_parameters[i]->mod_struct = mod_struct;

		//user_parameters[i]->burn_freqs = mcmc_frequencies;
		//user_parameters[i]->burn_data = mcmc_data;
		//user_parameters[i]->burn_noise = mcmc_noise;
		//user_parameters[i]->burn_lengths = mcmc_data_length;
	}
	//######################################################


	RJPTMCMC_MH_dynamic_PT_alloc_comprehensive(sampler_output,output,status, model_status, nested_model_number,
		max_dimension,min_dimension, N_steps, chain_N, 
		max_chain_N_thermo_ensemble,initial_pos,initial_status,initial_model_status,seeding_var, ensemble_initial_pos,ensemble_initial_status,ensemble_initial_model_status, chain_temps, 
		swp_freq, t0, nu,max_chunk_size,chain_distribution_scheme,
		log_prior,RJMCMC_2WF_likelihood_wrapper, RJMCMC_2WF_fisher_wrapper,RJMCMC_2WF_RJ_proposal_wrapper,(void**)user_parameters,numThreads, pool, 
		show_prog,update_RJ_widths,statistics_filename,
		chain_filename, likelihood_log_filename,checkpoint_filename);
	
	//Deallocate fftw plans
	for (int i =0;i<num_detectors;i++)
		deallocate_FFTW_mem(&plans[i]);
	//#################################################
	for(int i = 0 ; i<num_detectors; i++){
		delete [] burn_data[i];
		delete [] burn_freqs[i];
		delete [] burn_noise[i];
		deallocate_FFTW_mem(&burn_plans[i]);
	}
	//if(mod_struct->GAUSS_QUAD){
	//	delete [] fish_freqs;
	//	delete [] fish_weights;
	//	for(int i = 0 ;i<mcmc_num_detectors; i++){
	//		delete [] fish_psd[i];
	//	}
	//	delete [] fish_psd;
	//}
	delete [] burn_data;
	delete [] burn_lengths;
	delete [] burn_noise;
	delete [] burn_freqs;
	delete [] burn_plans;
	for(int i = 0 ; i<chain_N; i++){
		delete user_parameters[i];
		gsl_rng_free(rvec[i]);
	}
	delete [] rvec;
	delete [] user_parameters;
	//#################################################
	free(plans);
	if(local_seeding){ delete [] seeding_var;}
}

void pack_local_mod_structure(mcmc_data_interface *interface,
	double *param,
	int *status,
	std::string waveform_extended, 
	void *parameters, 
	MCMC_modification_struct *full_struct, 
	MCMC_modification_struct *local_struct )	
{
	if(waveform_extended.find("gIMR") != std::string::npos){
		int dimct = 0 ;
		int dphi_boundary = full_struct->gIMR_Nmod_phi + interface->min_dim;
		int dsigma_boundary = full_struct->gIMR_Nmod_sigma + dphi_boundary;
		int dbeta_boundary = full_struct->gIMR_Nmod_beta + dsigma_boundary;
		int dalpha_boundary = full_struct->gIMR_Nmod_alpha + dbeta_boundary;
		for(int i = 0 ; i<interface->max_dim; i++){
			if(status[i] == 1){
				dimct++;
			}
			if( i >= interface->min_dim){
				if(status[i] == 1 && i<dphi_boundary){
					local_struct->gIMR_Nmod_phi ++;
				}
				else if(status[i] == 1 && i<dsigma_boundary){
					local_struct->gIMR_Nmod_sigma ++;
				}
				else if(status[i] == 1 && i<dbeta_boundary){
					local_struct->gIMR_Nmod_beta ++;
				}
				else if(status[i] == 1 && i<dalpha_boundary){
					local_struct->gIMR_Nmod_alpha ++;
				}
			}
		}
		if(dimct != interface->min_dim){
			if(local_struct->gIMR_Nmod_phi != 0){
				local_struct->gIMR_phii = new int[local_struct->gIMR_Nmod_phi];
			}
			if(local_struct->gIMR_Nmod_sigma != 0){
				local_struct->gIMR_sigmai = new int[local_struct->gIMR_Nmod_sigma];
			}
			if(local_struct->gIMR_Nmod_beta != 0){
				local_struct->gIMR_betai = new int[local_struct->gIMR_Nmod_beta];
			}
			if(local_struct->gIMR_Nmod_alpha != 0){
				local_struct->gIMR_alphai = new int[local_struct->gIMR_Nmod_alpha];
			}
			
			int ct_phi = 0 ;
			int ct_sigma = 0 ;
			int ct_beta = 0 ;
			int ct_alpha = 0 ;
			for(int i =interface->min_dim ; i<interface->max_dim; i++){
				if(status[i] == 1){
					if(i < dphi_boundary){
						local_struct->gIMR_phii[ct_phi] = full_struct->gIMR_phii[i-dphi_boundary+full_struct->gIMR_Nmod_phi];
						ct_phi++;
					}
					else if(i < dsigma_boundary){
						local_struct->gIMR_sigmai[ct_sigma] = full_struct->gIMR_sigmai[i-dsigma_boundary+full_struct->gIMR_Nmod_sigma];
						ct_sigma++;
					}
					else if(i < dbeta_boundary){
						local_struct->gIMR_betai[ct_beta] = full_struct->gIMR_betai[i-dbeta_boundary+full_struct->gIMR_Nmod_beta];
						ct_beta++;
					}
					else if(i < dalpha_boundary){
						local_struct->gIMR_alphai[ct_alpha] = full_struct->gIMR_alphai[i-dalpha_boundary+full_struct->gIMR_Nmod_alpha];
						ct_alpha++;
					}
				}
			}
		}
	}

	return;
}


double RJMCMC_2WF_likelihood_wrapper(
	double *param, 
	int *status, 
	int model_status, 
	mcmc_data_interface *interface, 
	void *parameters)
{
	//debugger_print(__FILE__,__LINE__,mcmc_mod_struct->gIMR_phii[0]);	
	//return 2;
	MCMC_user_param *user_param = (MCMC_user_param *)parameters;

	int max_dimension = interface->max_dim;
	double ll = 0;
	double *temp_params = new double[max_dimension];
	int dimct = 0 ;
	for(int i = 0 ; i<interface->max_dim; i++){
		if(status[i] == 1){
			temp_params[dimct] = param[i];
			dimct++;
		}
	}
	bool WF2 = false;
	std::string gen_meth=mcmc_generation_method_base;
	if(dimct > interface->min_dim){ WF2 = true;gen_meth = mcmc_generation_method_extended;}
	//######################################################################
	//Pack up local mod_struct
	MCMC_modification_struct mod_struct_local;
	
	if(WF2){
		pack_local_mod_structure(interface,param,status,gen_meth, parameters , mcmc_mod_struct, &mod_struct_local);
	}
	
	//######################################################################
	//#########################################################################
	gen_params_base<double> gen_params;
	//std::string local_gen = MCMC_prep_params(param, 
	//	temp_params,&gen_params, dimct, gen_meth,&mod_struct_local);
	std::string local_gen = MCMC_prep_params(temp_params, 
		temp_params,&gen_params, dimct, gen_meth,&mod_struct_local);
	//#########################################################################
	//#########################################################################

	//repack_non_parameters(temp_params, &gen_params, 
		//"MCMC_"+mcmc_generation_method, dimension, NULL);
	repack_parameters(temp_params, &gen_params, 
		"MCMC_"+gen_meth, dimct, NULL);
	//#########################################################################
	//#########################################################################
	//return 1;
	std::complex<double> **local_data = mcmc_data;
	double **local_freqs = mcmc_frequencies;
	double **local_noise = mcmc_noise;
	double **local_weights = user_param->weights;
	int *local_lengths = mcmc_data_length;
	fftw_outline *local_plans = mcmc_fftw_plans;
	std::string local_integration_method="SIMPSONS";
	//if(interface->burn_phase && user_param->burn_data){
	if(false){
		local_data = user_param->burn_data;
		local_freqs = user_param->burn_freqs;
		local_noise = user_param->burn_noise;
		local_lengths = user_param->burn_lengths;
		local_plans = user_param->burn_plans;
	}
	if(user_param->GAUSS_QUAD){
		local_integration_method="GAUSSLEG";
	}
	if(mcmc_intrinsic){
		if(gen_meth.find("IMRPhenomD") != std::string::npos){
			if(!mcmc_save_waveform){
				for(int i=0; i < mcmc_num_detectors; i++){
					gen_params.theta=0;	
					gen_params.phi=0;	
					gen_params.psi=0;	
					gen_params.phiRef = 1;
					gen_params.f_ref = 10;
					gen_params.incl_angle=0;	
					gen_params.tc =1;
					std::complex<double> *response =
						(std::complex<double> *) malloc(sizeof(std::complex<double>) * local_lengths[i]);
					fourier_detector_response_horizon(local_freqs[i], local_lengths[i], response, mcmc_detectors[i], local_gen, &gen_params);
					ll += maximized_Log_Likelihood_aligned_spin_internal(local_data[i], 
							local_noise[i],
							local_freqs[i],
							response,
							(size_t) local_lengths[i],
							&local_plans[i]
							);
					free(response);
				}
			}
			else{
				gen_params.theta=0;	
				gen_params.phi=0;	
				gen_params.psi=0;	
				gen_params.phiRef = 1;
				gen_params.f_ref = 10;
				gen_params.incl_angle=0;	
				gen_params.tc =1;
				std::complex<double> *response =
					(std::complex<double> *) malloc(sizeof(std::complex<double>) * local_lengths[0]);
				fourier_detector_response_horizon(local_freqs[0], local_lengths[0], response, mcmc_detectors[0], local_gen, &gen_params);
				for(int i=0; i < mcmc_num_detectors; i++){
					ll += maximized_Log_Likelihood_aligned_spin_internal(local_data[i], 
							local_noise[i],
							local_freqs[i],
							response,
							(size_t) local_lengths[i],
							&local_plans[i]
							);
				}
				free(response);
			}

		}
		else if(gen_meth.find("IMRPhenomP")!=std::string::npos){
			//if(!mcmc_save_waveform){
			if(false){
			}
			else{
				gen_params.theta=0;	
				gen_params.phi=0;	
				gen_params.psi=0;	
				gen_params.phiRef = 1;
				gen_params.f_ref = 20;
				gen_params.incl_angle=0;	
				gen_params.tc =1;
				//fourier_waveform(local_freqs[0],local_lengths[0], hp,hc, local_gen, &gen_params);
				waveform_polarizations<double> wp;
				assign_polarizations(gen_meth, &wp);	
				wp.allocate_memory(local_lengths[0]);	
				fourier_waveform(local_freqs[0],local_lengths[0], &wp, local_gen, &gen_params);
				for(int i=0; i < mcmc_num_detectors; i++){
					ll += maximized_Log_Likelihood_unaligned_spin_internal(local_data[i], 
							local_noise[i],
							local_freqs[i],
							wp.hplus,
							wp.hcross,
							(size_t) local_lengths[i],
							&local_plans[i]
							);
				}
				wp.deallocate_memory();	
			}

		}
	}
	else{
		double RA = gen_params.RA;
		double DEC = gen_params.DEC;
		double PSI = gen_params.psi;
		
		ll =  MCMC_likelihood_extrinsic(mcmc_save_waveform, 
			&gen_params,local_gen, local_lengths, 
			local_freqs, local_data, local_noise,local_weights, local_integration_method, user_param->log10F, mcmc_detectors, 
			 mcmc_num_detectors);
	}
	//Cleanup
	delete [] temp_params;
	//We DO need to delete gIMR index arrays because we're making a copy, unlike in 
	//regular MCMC
	if(check_mod(local_gen)){
		if( local_gen.find("ppE") != std::string::npos ||
			local_gen.find("dCS") != std::string::npos ||
			local_gen.find("EdGB") != std::string::npos){
			delete [] gen_params.betappe;
		}
		else if( local_gen.find("gIMR") != std::string::npos){
			if(mod_struct_local.gIMR_Nmod_phi !=0){
				delete [] gen_params.delta_phi;
				delete [] gen_params.phii;
			}
			if(mod_struct_local.gIMR_Nmod_sigma !=0){
				delete [] gen_params.delta_sigma;
				delete [] gen_params.sigmai;
			}
			if(mod_struct_local.gIMR_Nmod_beta !=0){
				delete [] gen_params.delta_beta;
				delete [] gen_params.betai;
			}
			if(mod_struct_local.gIMR_Nmod_alpha !=0){
				delete [] gen_params.delta_alpha;
				delete [] gen_params.alphai;
			}

		}
	}
	return ll;
}

void RJMCMC_2WF_fisher_wrapper(
	double *param, 
	int *status, 
	int model_status, 
	double **fisher,
	mcmc_data_interface *interface, 
	void *parameters)
{


	MCMC_user_param *user_param = (MCMC_user_param *)parameters;

	int max_dimension = interface->max_dim;
	int min_dimension = interface->min_dim;
	double ll = 0;
	double *temp_params = new double[min_dimension];

	//#########################################################################
	gen_params_base<double> gen_parameters;
	std::string local_gen = MCMC_prep_params(param, 
		temp_params,&gen_parameters, min_dimension, mcmc_generation_method_base,mcmc_mod_struct);
	//#########################################################################
	//#########################################################################

	//#########################################################################
	//#########################################################################
	repack_parameters(param, &gen_parameters, 
		"MCMC_"+mcmc_generation_method_base, min_dimension, NULL);
	//#########################################################################
	//#########################################################################
	for(int j =0; j<min_dimension; j++){
		for(int k =0; k<min_dimension; k++)
		{
			fisher[j][k] =0;
		}
	} 
	double **temp_out = allocate_2D_array(min_dimension,min_dimension);
	for(int i =0 ; i <mcmc_num_detectors; i++){
		
		//Use AD 
		if(user_param->GAUSS_QUAD)
		{	
			std::unique_lock<std::mutex> lock{*(user_param->mFish)};
			//fisher_autodiff(mcmc_frequencies[i], mcmc_data_length[i],
			//	"MCMC_"+mcmc_generation_method, mcmc_detectors[i],mcmc_detectors[0],temp_out,dimension, 
			//	(gen_params *)(&params),  "SIMPSONS",(double *)NULL,false,mcmc_noise[i]);
			fisher_autodiff(user_param->fisher_freq[i], user_param->fisher_length[i],
				"MCMC_"+mcmc_generation_method_base, mcmc_detectors[i],mcmc_detectors[0],temp_out,min_dimension, 
				(gen_params *)(&gen_parameters),  "GAUSSLEG",user_param->fisher_weights[i],true,user_param->fisher_PSD[i]);
		}
		else{
			fisher_numerical(mcmc_frequencies[i], mcmc_data_length[i],
				"MCMC_"+mcmc_generation_method_base, mcmc_detectors[i],mcmc_detectors[0],temp_out,min_dimension, 
				&gen_parameters, mcmc_deriv_order, NULL, NULL, mcmc_noise[i]);

		}
		for(int j =0; j<min_dimension; j++){
			for(int k =0; k<min_dimension; k++)
			{
				fisher[j][k] +=temp_out[j][k];
				//if(std::isnan(fisher[j][k]))
				//{
				//      std::cout<<j<<" "<<k<<" "<<temp_out[j][k]<<std::endl;
				//}
			}
		} 
	}
	//Add prior information to fisher
	//if(mcmc_generation_method.find("Pv2") && !mcmc_intrinsic){
	if(!mcmc_intrinsic){
		fisher[0][0] += 1./(4*M_PI*M_PI);//RA
		fisher[1][1] += 1./4;//sin DEC
		fisher[2][2] += 1./(4*M_PI*M_PI);//psi
		fisher[3][3] += 1./(4);//cos iota
		fisher[4][4] += 1./(4*M_PI*M_PI);//phiref
		fisher[5][5] += 1./(.01);//tc
		fisher[8][8] += 1./.25;//eta
		fisher[9][9] += 1./4;//spin1
		fisher[10][10] += 1./4;//spin2
		if(mcmc_generation_method_base.find("PhenomPv2") != std::string::npos){
			fisher[11][11] += 1./4;//cos theta1
			fisher[12][12] += 1./4;//cos theta2
			fisher[13][13] += 1./(4*M_PI*M_PI);//phi1
			fisher[14][14] += 1./(4*M_PI*M_PI);//phi2
		}
	}
	else{
		if(mcmc_generation_method_base.find("PhenomPv2") != std::string::npos){
			fisher[1][1] =1./(.25) ;//eta
			fisher[2][2] =1./(4);//spin1
			fisher[3][3] =1./(4);//spin2
			fisher[4][4] =1./(4);//cos theta1
			fisher[5][5] =1./(4);//cos theta2
			fisher[6][6] =1./(4*M_PI*M_PI) ;//phi1
			fisher[7][7] =1./(4*M_PI*M_PI) ;//phi2
		}
		else if (mcmc_generation_method_base.find("PhenomD")!=std::string::npos){
			fisher[1][1] =1./(.25) ;//eta
			fisher[2][2] =1./(4) ;//spin1
			fisher[3][3] =1./(4) ;//spin2
	
		}
	}
	deallocate_2D_array(temp_out, min_dimension,min_dimension);
	//////////////////////////////////////////////
	//if(!interface->burn_phase)
	//{
	//	debugger_print(__FILE__,__LINE__,"Fisher MCMC");
	//	double **cov = allocate_2D_array( dimension,dimension);
	//	gsl_cholesky_matrix_invert(output, cov, dimension);
	//	for(int i = 0 ; i<dimension; i++){
	//		std::cout<<sqrt(cov[i][i])<<" ";	
	//		//for(int j = 0 ; j<dimension; j++){
	//		//	std::cout<<cov[i][j]<<" ";	
	//		//}
	//		//std::cout<<std::endl;	
	//		
	//	}
	//	std::cout<<std::endl;	
	//	deallocate_2D_array(cov, dimension,dimension);
	//}
	//////////////////////////////////////////////

	//Cleanup
	delete [] temp_params;
	if(check_mod(local_gen)){
		if(local_gen.find("ppE") != std::string::npos ||
			local_gen.find("dCS")!=std::string::npos||
			local_gen.find("EdGB")!=std::string::npos){
			delete [] gen_parameters.betappe;
		}
		else if( local_gen.find("gIMR") != std::string::npos){
			if(mcmc_mod_struct ->gIMR_Nmod_phi !=0){
				delete [] gen_parameters.delta_phi;
			}
			if(mcmc_mod_struct ->gIMR_Nmod_sigma !=0){
				delete [] gen_parameters.delta_sigma;
			}
			if(mcmc_mod_struct ->gIMR_Nmod_beta !=0){
				delete [] gen_parameters.delta_beta;
			}
			if(mcmc_mod_struct ->gIMR_Nmod_alpha !=0){
				delete [] gen_parameters.delta_alpha;
			}

		}
	}
	return ;
}
void RJMCMC_2WF_RJ_proposal_wrapper(
	double *current_param, 
	double *proposed_param, 
	int *current_status, 
	int *proposed_status, 
	int *current_model_status, 
	int *proposed_model_status, 
	double *MH_corrections,
	mcmc_data_interface *interface, 
	void *parameters)
{

	*(MH_corrections) = 0;
	int dimct = 0 ;
	for(int i = 0 ; i<interface->max_dim; i++){
		if(current_status[i] == 1){
			dimct++;
		}
	}


	MCMC_user_param *user_param = (MCMC_user_param *)parameters;

	int max_dimension = interface->max_dim;
	int min_dimension = interface->min_dim;
	for(int i = 0 ; i<max_dimension; i++){
		proposed_param[i] = current_param[i];
		proposed_status[i] = current_status[i];
	}
	if(dimct> min_dimension && dimct< max_dimension){
		double alpha = gsl_rng_uniform(user_param->r);
		//Add dimension
		if(  alpha> 0.5){
			int beta = (int)(gsl_rng_uniform(user_param->r)*( max_dimension - dimct));
			int ct = 0 ;
			for (int i =min_dimension  ; i<max_dimension; i++){
				if(current_status[i] == 0){
					if(ct == beta){
						proposed_status[i] = 1;
						//proposed_param[i] = gsl_ran_gaussian(user_param->r,interface->RJ_step_width);
						proposed_param[i] = gsl_rng_uniform(user_param->r) 
							* (user_param->mod_prior_ranges[i-min_dimension][1]-user_param->mod_prior_ranges[i-min_dimension][0]) 
							+ user_param->mod_prior_ranges[i-min_dimension][0];
						*(MH_corrections)-=log(1./(user_param->mod_prior_ranges[i-min_dimension][1]
							-user_param->mod_prior_ranges[i-min_dimension][0]));
						//proposed_param[i] = gsl_ran_gaussian(user_param->r,10);
						break;
					}
					ct++;
				}
			}
		}
		//Remove dimension
		else if( alpha< 0.5){
			int beta = (int)(gsl_rng_uniform(user_param->r)*( -min_dimension + dimct));
			int ct = 0 ;
			for (int i =min_dimension  ; i<max_dimension; i++){
				if(current_status[i] == 1){
					if(ct == beta){
						proposed_status[i] = 0;
						proposed_param[i] = 0;
						*(MH_corrections)+=log(1./(user_param->mod_prior_ranges[i-min_dimension][1]
							-user_param->mod_prior_ranges[i-min_dimension][0]));
						break;
					}
					ct++;
				}
			}


		}
	}
	else if(dimct ==  min_dimension){
		int beta = (int)(gsl_rng_uniform(user_param->r)*(  max_dimension - dimct));
		int ct = 0 ;
		for (int i =min_dimension  ; i<max_dimension; i++){
			if(current_status[i] == 0){
				if(ct == beta){
					proposed_status[i] = 1;
					//proposed_param[i] = gsl_ran_gaussian(user_param->r,interface->RJ_step_width);
					proposed_param[i] = gsl_rng_uniform(user_param->r) 
						* (user_param->mod_prior_ranges[i-min_dimension][1]-user_param->mod_prior_ranges[i-min_dimension][0]) 
						+ user_param->mod_prior_ranges[i-min_dimension][0];
					*(MH_corrections)-=log(1./.5/(user_param->mod_prior_ranges[i-min_dimension][1]
						-user_param->mod_prior_ranges[i-min_dimension][0]));
					//proposed_param[i] = gsl_ran_gaussian(user_param->r,10);
					break;
				}
				ct++;
			}
		}
	}
	else if(dimct ==  max_dimension){
		int beta = (int)(gsl_rng_uniform(user_param->r)*(  -min_dimension + dimct));
		int ct = 0 ;
		for (int i =min_dimension  ; i<max_dimension; i++){
			if(current_status[i] == 1){
				if(ct == beta){
					proposed_status[i] = 0;
					proposed_param[i] = 0;
					*(MH_corrections)+=log(.5*1./(user_param->mod_prior_ranges[i-min_dimension][1]
						-user_param->mod_prior_ranges[i-min_dimension][0]));
					break;
				}
				ct++;
			}
		}
	}
	//#########################################################
	//#########################################################
	//TESTING
	//if(proposed_status[max_dimension-1] == 0 ){
	//	proposed_status[max_dimension-1] = 1;
	//	proposed_param[max_dimension-1] = gsl_ran_gaussian(user_param->r,interface->RJ_step_width);
	//}
	//else{
	//	proposed_status[max_dimension-1] = 0;
	//	proposed_param[max_dimension-1] = 0;
	//}
	//#########################################################
	//#########################################################
	
	return ;
}





