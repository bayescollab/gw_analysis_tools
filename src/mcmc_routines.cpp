#include "mcmc_routines.h"
#include "waveform_generator.h"
#include "util.h"
#include "noise_util.h"
#include "waveform_util.h"
#include "fisher.h"
#include "mcmc_sampler.h"
#include <iostream>
#include <vector>
#include <complex>
#include <fftw3.h>
#include <algorithm>
#include <iostream>



//#####################################################################
//#####################################################################
//
//MEMORY LEAK IN MCMC ROUTINES -- MUST FIX
//
//#####################################################################
//#####################################################################

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
	params.phic=0;
	params.tc=0;
	params.NSflag=NSflag;
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
	initiate_likelihood_function(&plan, length);
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
	deactivate_likelihood_function(&plan);
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
	params.phic=0.;
	params.tc=0.0;
	params.NSflag=NSflag;

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

	double out =  maximized_coal_Log_Likelihood_aligned_spin_internal(data,
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
	initiate_likelihood_function(&plan,length);
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
	deactivate_likelihood_function(&plan);

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

void initiate_likelihood_function(fftw_outline *plan, int length)
{
	plan->in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * length);	
	plan->out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * length);	
	plan->p = fftw_plan_dft_1d(length, plan->in, plan->out,FFTW_FORWARD, FFTW_MEASURE);
}
void deactivate_likelihood_function(fftw_outline *plan)
{
	fftw_destroy_plan(plan->p);
	fftw_free(plan->in);
	fftw_free(plan->out);
	fftw_cleanup();
}


double maximized_coal_Log_Likelihood(std::complex<double> *data, 
				double *psd,
				double *frequencies,
				size_t length,
				gen_params *params,
				std::string detector,
				std::string generation_method,
				fftw_outline *plan
				)
{
	double ll = 0;
	if(	generation_method  == "IMRPhenomD" ||
		generation_method  == "ppE_IMRPhenomD_Inspiral" ||
		generation_method  == "ppE_IMRPhenomD_IMR" ){
		//std::cout<<params->mass1<<" "<<params->mass2<<" "<<params->spin1[2]<<" "<<params->spin2[2]<<std::endl;	
		std::complex<double> *response =
			(std::complex<double> *) malloc(sizeof(std::complex<double>) * length);
		fourier_detector_response(frequencies, length, response, detector, generation_method, params);
		ll = maximized_coal_Log_Likelihood_aligned_spin_internal(data, psd, frequencies,response, length, plan);
		
		free(response);
		}
	else if ( 	generation_method == "IMRPhenomPv2" ||
			generation_method == "ppE_IMRPhenomPv2_Inspiral" ||
			generation_method == "ppE_IMRPhenomPv2_IMR" ){
				
		std::complex<double> *hp =
				(std::complex<double> *) malloc(sizeof(std::complex<double>) * length);
		std::complex<double> *hc =
				(std::complex<double> *) malloc(sizeof(std::complex<double>) * length);
		fourier_waveform(frequencies,length,hp,hc,generation_method,params);
		ll = maximized_coal_Log_Likelihood_unaligned_spin_internal(data,psd,frequencies,hp,hc, length, plan);

		free(hp);
		free(hc);
	}
	return ll;
}
					
double maximized_coal_Log_Likelihood(double *data_real, 
				double *data_imag,
				double *psd,
				double *frequencies,
				size_t length,
				gen_params *params,
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
	
	double ll = maximized_coal_Log_Likelihood(data, 
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

/*! \brief Unmarginalized log of the likelihood
 *
 */
double Log_Likelihood(std::complex<double> *data, 
				double *psd,
				double *frequencies,
				size_t length,
				gen_params *params,
				std::string detector,
				std::string generation_method,
				fftw_outline *plan
				)
{
	double ll = 0;

	std::complex<double> *detect_response =
			(std::complex<double> *) malloc(sizeof(std::complex<double>) * length);
	fourier_detector_response(frequencies,length,detect_response,detector, generation_method,params);
	ll = Log_Likelihood_internal(data,psd,frequencies,detect_response, length, plan);

	free(detect_response);
	return ll;
}
/*! \brief Maximized match over coalescence variables - returns log likelihood NOT NORMALIZED for aligned spins
 *
 * Note: this function is not properly normalized for an absolute comparison. This is made for MCMC sampling, so to minimize time, constant terms like (Data|Data), which would cancel in the Metropolis-Hasting ratio, are left out for efficiency
 */
double maximized_coal_Log_Likelihood_aligned_spin_internal(std::complex<double> *data,
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

	free(integrand);
	free(g);
	fftw_free(in);
	fftw_free(out);
	//std::cout<<"inner products: "<<max<<" "<<HH<<std::endl;

	//return -0.5*(HH- 2*max);
	//std::cout<<"ll: "<<.5*(max)/HH<<std::endl;
	//std::cout<<"psd: "<<psd[1000]<<std::endl;
	//std::cout<<"freq: "<<frequencies[1000]<<std::endl;
	//std::cout<<"data: "<<data[1000]<<std::endl;
	//std::cout<<"SNR**2 template "<<HH<<std::endl;
	//return .5*(max*max)/HH;
	return .5*(max)/HH;
}

/*! \brief log likelihood function that maximizes over extrinsic parameters tc, phic, D, and phiRef, the reference frequency - for unaligned spins 
 *
 * Ref: arXiv 1603.02444v2
 */
double maximized_coal_Log_Likelihood_unaligned_spin_internal(std::complex<double> *data,
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
	std::complex<double> *detector_response = 
		(std::complex<double> *)malloc(sizeof(std::complex<double>)*length);
	//fourier_detector_response(frequencies, length, hplus,hcross, detector_response, 0.0,0.0,"Hanford");

	//for (int i =0;i< length;i++)
	//	integrand[i] = real(detector_response[i]*std::conj(detector_response[i]))/psd[i];
	//integral = 4.*simpsons_sum(delta_f, length, integrand);
	//double HH = integral;

	//Calculate template snr for plus polarization sqrt(<H+|H+>)
	for (int i =0;i< length;i++)
		integrand[i] = real(hplus[i]*std::conj(hplus[i]))/psd[i];
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
		rhoplus[i] = std::complex<double>(plan->out[i][0],plan->out[i][1]);
		//Norm of the output, squared (Re{g}^2 + Im{g}^2)
		rhoplus2[i] = out[i][0]* out[i][0]+ out[i][1]* out[i][1];
		
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
		rhocross[i] = std::complex<double>(plan->out[i][0],plan->out[i][1]);
		//Norm of the output, squared (Re{g}^2 + Im{g}^2)
		rhocross2[i] = out[i][0]* out[i][0]+ out[i][1]* out[i][1];
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
	double max = .25 * (*std::max_element(lambda, lambda+length))*delta_f; 

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
			std::complex<double> *detector_response,
			int length,
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
	double integral = 4.*simpsons_sum(delta_f, length, integrand);
	double HH = integral;

	for (int i =0;i< length;i++)
		integrand[i] = real(data[i]*std::conj(detector_response[i]))/psd[i];
	integral = 4.*simpsons_sum(delta_f, length, integrand);
	double DH = integral;
	//std::cout<<"inner products: "<<DH<<" "<<HH<<std::endl;

	free(integrand);

	return -0.5*(HH- 2*DH);
}

/*! \brief Wrapper for the MCMC_MH function, specifically for GW analysis
 *
 * Handles the details of setting up the MCMC sampler and wraps the fisher and log likelihood to conform to the format of the sampler
 */
void MCMC_MH_GW(double ***output,
		int dimension,
		int N_steps,
		int chain_N,
		double *initial_pos,
		double *chain_temps,
		int swp_freq,
		double(*log_prior)(double *param, int dimension),
		int num_detectors,
		std::complex<double> **data,
		double **noise_psd,
		double **frequencies,
		int *data_length,
		std::string *detectors,
		std::string generation_method,
		std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
		std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output*/
		std::string auto_corr_filename/**< Filename to output auto correlation in some interval, if empty string, not output*/
					)
{
	//Create fftw plan for each detector (length of data stream may be different)
	fftw_outline *plans= (fftw_outline *)malloc(sizeof(fftw_outline)*num_detectors);
	for (int i =0;i<num_detectors;i++)
	{	
		initiate_likelihood_function(&plans[i] , data_length[i]);
	}
	mcmc_noise = noise_psd;	
	mcmc_frequencies = frequencies;
	mcmc_data = data;
	mcmc_data_length = data_length;
	mcmc_detectors = detectors;
	mcmc_generation_method = generation_method;
	mcmc_fftw_plans = plans;
	mcmc_num_detectors = num_detectors;
	if(dimension==4 && generation_method =="IMRPhenomD"){
		std::cout<<"Sampling in parameters: chirpmass, eta, chi1, chi2"<<std::endl;
	}
	else if(dimension==7 && generation_method =="IMRPhenomD"){
		std::cout<<"Sampling in parameters: DL, tc, phic, chirpmass, eta, chi1, chi2"<<std::endl;
	}
	else{
		std::cout<<
			"Only supported options are currently 1 detector and PhenomD"<<std::endl;
		exit(1);
	}
	MCMC_MH(output, dimension, N_steps, chain_N, initial_pos, chain_temps, swp_freq,
		 log_prior,MCMC_likelihood_wrapper, MCMC_fisher_wrapper,statistics_filename,
		chain_filename,auto_corr_filename);
	//MCMC_MH(output, dimension, N_steps, chain_N, initial_pos, chain_temps, swp_freq,
	//	 log_prior,MCMC_likelihood_wrapper, NULL,statistics_filename,
	//	chain_filename,auto_corr_filename);
	
	//Deallocate fftw plans
	for (int i =0;i<num_detectors;i++)
		deactivate_likelihood_function(&plans[i]);
	free(plans);
}

void MCMC_fisher_wrapper(double *param, int dimension, double **output)
{
	//if(mcmc_num_detectors ==1 && mcmc_generation_method =="IMRPhenomD"){	
	if(dimension ==4 && mcmc_generation_method =="IMRPhenomD"){	
		//unpack parameter vector
		//double dl_prime = std::exp(param[0])/MPC_SEC;
		double dl_prime = 1000;
		double chirpmass = std::exp(param[0])/MSOL_SEC;
		double eta = param[1];
		double chi1 = param[2];
		double chi2 = param[3];
	
		//create gen_param struct
		gen_params parameters; 
		parameters.mass1 = calculate_mass1(chirpmass, eta);
		parameters.mass2 = calculate_mass2(chirpmass, eta);
		parameters.spin1[0] = 0;
		parameters.spin1[1] = 0;
		parameters.spin1[2] = chi1;
		parameters.spin2[0] = 0;
		parameters.spin2[1] = 0;
		parameters.spin2[2] = chi2;
		parameters.Luminosity_Distance = dl_prime;
		//The rest is maximized over for this option
		parameters.tc = 0;
		parameters.phic = 0;
		parameters.incl_angle = 0;
		parameters.phi=0;
		parameters.theta=0;
		parameters.NSflag = false;
		parameters.sky_average = false;
		
		for(int j =0; j<dimension; j++){
			for(int k =0; k<dimension; k++)
			{
				output[j][k] =0;
			}
		} 
		double **temp_out = allocate_2D_array(dimension,dimension);
		for (int i =0; i<mcmc_num_detectors; i++){
			fisher(mcmc_frequencies[i], mcmc_data_length[i],
				"MCMC_"+mcmc_generation_method+"_single_detect", 
				mcmc_detectors[i], temp_out, 4, &parameters, 
				NULL, NULL, mcmc_noise[i]);
			for(int j =0; j<dimension; j++){
				for(int k =0; k<dimension; k++)
				{
					output[j][k] +=temp_out[j][k];
				}
			} 
		}


		deallocate_2D_array(temp_out, dimension,dimension);
	}
	else if(dimension ==7 && mcmc_generation_method =="IMRPhenomD"){	
		//unpack parameter vector
		double dl_prime = std::exp(param[0])/MPC_SEC;
		double tc = std::exp(param[1]);
		double phic = std::exp(param[2]);
		double chirpmass = std::exp(param[3])/MSOL_SEC;
		double eta = param[4];
		double chi1 = param[5];
		double chi2 = param[6];
	
		//create gen_param struct
		gen_params parameters; 
		parameters.mass1 = calculate_mass1(chirpmass, eta);
		parameters.mass2 = calculate_mass2(chirpmass, eta);
		parameters.spin1[0] = 0;
		parameters.spin1[1] = 0;
		parameters.spin1[2] = chi1;
		parameters.spin2[0] = 0;
		parameters.spin2[1] = 0;
		parameters.spin2[2] = chi2;
		parameters.Luminosity_Distance = dl_prime;
		//The rest is maximized over for this option
		parameters.tc = tc;
		parameters.phic = phic;
		parameters.incl_angle = 0;
		parameters.phi=0;
		parameters.theta=0;
		parameters.NSflag = false;
		parameters.sky_average = false;
	
		//*NOTE* Current fisher is log \eta -- sampler is in \eta -- 
		//Leaving for now, but that should change too
		//fisher(mcmc_frequencies[0], mcmc_data_length[0],"MCMC_"+mcmc_generation_method+"_ind_spins", mcmc_detectors[0], output, 7, &parameters, NULL, NULL, mcmc_noise[0]);
		
		for(int j =0; j<dimension; j++){
			for(int k =0; k<dimension; k++)
			{
				output[j][k] =0;
			}
		} 
		double **temp_out = allocate_2D_array(dimension,dimension);
		for (int i =0; i<mcmc_num_detectors; i++){
			fisher(mcmc_frequencies[i], 
				mcmc_data_length[i],"MCMC_"+mcmc_generation_method+"_ind_spins",
				mcmc_detectors[i], output, 7, &parameters, NULL, NULL, 
				mcmc_noise[i]);
			for(int j =0; j<dimension; j++){
				for(int k =0; k<dimension; k++)
				{
					output[j][k] +=temp_out[j][k];
				}
			} 
		}


		deallocate_2D_array(temp_out, dimension,dimension);
	}
}

double MCMC_likelihood_wrapper(double *param, int dimension)
{
	double ll = 0;
	//if(mcmc_num_detectors ==1 && mcmc_generation_method =="IMRPhenomD"){	
	if(dimension ==4 && mcmc_generation_method =="IMRPhenomD"){	
	//if(false){	
		//unpack parameter vector
		//double dl_prime = std::exp(param[0])/MPC_SEC;
		double dl_prime = 1000;
		double chirpmass = std::exp(param[0])/MSOL_SEC;
		double eta = param[1];
		double chi1 = param[2];
		double chi2 = param[3];
		//std::cout<<"LL: "<<chirpmass<<" "<<eta<<" "<<chi1<<" "<<chi2<<std::endl;	
		//create gen_param struct
		gen_params parameters; 
		parameters.mass1 = calculate_mass1(chirpmass, eta);
		parameters.mass2 = calculate_mass2(chirpmass, eta);
		parameters.spin1[0] = 0;
		parameters.spin1[1] = 0;
		parameters.spin1[2] = chi1;
		parameters.spin2[0] = 0;
		parameters.spin2[1] = 0;
		parameters.spin2[2] = chi2;
		parameters.Luminosity_Distance = dl_prime;
		//The rest is maximized over for this option
		parameters.tc = 0;
		parameters.phic = 0;
		parameters.incl_angle = 0;
		parameters.phi=0;
		parameters.theta=0;
		parameters.NSflag = false;
		parameters.sky_average = false;

		//calculate log likelihood
		for(int i=0; i < mcmc_num_detectors; i++){
			ll += maximized_coal_Log_Likelihood(mcmc_data[i], 
					mcmc_noise[i],
					mcmc_frequencies[i],
					(size_t) mcmc_data_length[i],
					&parameters,
					mcmc_detectors[i],
					mcmc_generation_method,
					&mcmc_fftw_plans[i]
					);
		}
	}
	if(dimension ==7 && mcmc_generation_method =="IMRPhenomD"){	
	//if(false){	
		//unpack parameter vector
		double dl_prime = std::exp(param[0])/MPC_SEC;
		double tc = std::exp(param[1]);
		double phic = std::exp(param[2]);
		double chirpmass = std::exp(param[3])/MSOL_SEC;
		double eta = param[4];
		double chi1 = param[5];
		double chi2 = param[6];
	
		//create gen_param struct
		gen_params parameters; 
		parameters.mass1 = calculate_mass1(chirpmass, eta);
		parameters.mass2 = calculate_mass2(chirpmass, eta);
		parameters.spin1[0] = 0;
		parameters.spin1[1] = 0;
		parameters.spin1[2] = chi1;
		parameters.spin2[0] = 0;
		parameters.spin2[1] = 0;
		parameters.spin2[2] = chi2;
		parameters.Luminosity_Distance = dl_prime;
		//The rest is maximized over for this option
		parameters.tc = tc;
		parameters.phic = phic;
		parameters.incl_angle = 0;
		parameters.phi=0;
		parameters.theta=0;
		parameters.NSflag = false;
		parameters.sky_average = false;

		
		//calculate log likelihood
		for(int i=0; i < mcmc_num_detectors; i++){
			ll += Log_Likelihood(mcmc_data[i], 
					mcmc_noise[i],
					mcmc_frequencies[i],
					(size_t) mcmc_data_length[i],
					&parameters,
					mcmc_detectors[i],
					mcmc_generation_method,
					&mcmc_fftw_plans[i]
					);
		}
	}
	return ll;

	//testing detailed balance
	//return 2;
}

			

