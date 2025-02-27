#include "waveform_util.h"
#include "util.h"
#include "GWATConfig.h"
#include "ortho_basis.h"
#include "waveform_generator.h"
#include "IMRPhenomP.h"
#include "IMRPhenomD.h"
#include "gIMRPhenomD.h"
#include "gIMRPhenomP.h"
#include "IMRPhenomD_NRT.h"
#include "ppE_IMRPhenomD.h"
#include "ppE_IMRPhenomP.h"
#include "ppE_utilities.h"
#include "detector_util.h"
#include "pn_waveform_util.h"
#include <fftw3.h>
#include <algorithm>
#include <complex>
#include <vector>
#include <string>
#include <adolc/taping.h>
#include <adolc/adouble.h>
#include <adolc/drivers/drivers.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_errno.h>
/*!\file 
 * Utilities for waveforms - SNR calculation and detector response
 * 	
 * includes snr and detector response
 *
 * Includes some utilities useful for MCMC and fisher calculations, as well as time-frequency methods for detectors like LISA
 */


/*\brief Calculates the match maximized over phic  and tc 
 *
 * Taken from III A of https://arxiv.org/pdf/1807.07163.pdf
 *
 */
double match(  std::complex<double> *data1, std::complex<double> *data2, double *SN,double *frequencies,int length)
{
	std::complex<double> Gtilde ;
	double *G= new double[length];
	double delta_f = frequencies[1]-frequencies[0];

	double norms[3] ;
	norms[0] =  data_snr(frequencies,length,data1,data1,SN);
	norms[1] =  data_snr(frequencies,length,data2,data2,SN);

	//norms[2] =  data_snr(frequencies,length,data1,data2,SN);
	//return norms[2]*norms[2]/norms[1]/norms[0];
	
	fftw_complex *in, *out; 
	fftw_plan p;
	in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * length);
        out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * length);
        p = fftw_plan_dft_1d(length, in, out,FFTW_BACKWARD, FFTW_MEASURE);
        for (int i=0;i<length; i++)
        {
                Gtilde = conj(data1[i]) * data2[i] / SN[i];
                in[i][0] = real(Gtilde);
                in[i][1] = imag(Gtilde);
        }

        fftw_execute(p);

        for (int i=0;i<length; i++)
        {
                //G[i] = std::abs(std::complex<double>(out[i][0],out[i][1])) ;
                G[i] = sqrt(out[i][0]*out[i][0]+out[i][1]*out[i][1]) ;
        }
	double max = G[0];
	for(int i= 0 ; i<length; i++){
		if(G[i]>max){max = G[i];}
	}
	max*=delta_f;

        //double max = (*std::max_element(G, G+length))*delta_f;
	
	fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);
	fftw_cleanup();

	delete [] G;
	
	return 4.*max/(norms[0]*norms[1]);
}


struct gsl_snr_struct
{
	gen_params_base<double> *params;
	std::string SN;
	std::string generation_method;
	std::string detector;
	//Extra parameters for root finding
	gsl_integration_workspace *w;
	int np;
	double fmin;
	double fmax;
	double relerr;
	double T_obs;
	double T_wait;
};

double data_snr(double *frequencies, 
	int length,
	std::complex<double> *data,
	std::complex<double> *response,
	double *psd
	)	
{
	double *integrand= (double *)malloc(sizeof(double)*length);

	double delta_f = frequencies[1]-frequencies[0];

	
	for (int i = 0; i<length;i++)
		integrand[i] = 4.*real(conj(data[i])*response[i])/psd[i]; 
	double inner_prod = (simpsons_sum(delta_f,length, integrand));
	//double inner_prod = (simpsons_sum(delta_f,length, integrand));
	free(integrand);
	return sqrt(inner_prod);

}
template <class T>
void create_coherent_GW_detection(
	std::string *detectors,
	int detector_N, 
	T **frequencies, 
	int *lengths,
	bool reuse_WF,/**< If using the same exact frequencies for each detector, we can save a lot of time by only computing the waveform once*/
	gen_params_base<T> *gen_params,/**<tc in gen_params refers to the time of coalescence, relative to the initial time of the data stream, for the FIRST detector, the other detectors are shifted by the appropriate amount*/
	std::string generation_method,
	std::complex<T> **responses/**< [out] Responses for the source at each detector, same order as detectors parameter -- should be pre allocated shape [detector_N][lengths[i]] */
	)
{
	if(reuse_WF){
		create_coherent_GW_detection_reuse_WF(detectors, detector_N, frequencies[0],lengths[0], gen_params, generation_method, responses);
	}
	else{
		debugger_print(__FILE__,__LINE__,"Independent WFs for each detector currently not supported!");	
	}

}
template void create_coherent_GW_detection<double>(std::string *, int, double**, int *, bool, gen_params_base<double> *,std::string, std::complex<double> **);
template void create_coherent_GW_detection<adouble>(std::string *, int, adouble**, int *, bool, gen_params_base<adouble> *,std::string, std::complex<adouble> **);

/*Note -- this only works with equatorial coordinates*/
template <class T>
void create_coherent_GW_detection_reuse_WF(
	std::string *detectors,
	int detector_N, 
	T *frequencies, 
	int length,
	gen_params_base<T> *gen_params,/**<tc in gen_params refers to the time of coalescence, relative to the initial time of the data stream, for the FIRST detector, the other detectors are shifted by the appropriate amount*/
	std::string generation_method,
	std::complex<T> **responses/**< [out] Responses for the source at each detector, same order as detectors parameter -- should be pre allocated shape [detector_N][length] */
	)
{
	waveform_polarizations<T> wp;
	assign_polarizations(generation_method, &wp);
	wp.allocate_memory(length);
	T tc_ref = gen_params->tc, tc = 0;
	gen_params->tc=tc_ref;
	fourier_waveform(frequencies, length, &wp,generation_method, gen_params);
	T DTOA = 0;
	for (int i = 0 ; i<detector_N; i++){
		DTOA = DTOA_DETECTOR(gen_params->RA,gen_params->DEC,gen_params->gmst, detectors[0],detectors[i]);
		//tc = tc_ref-DTOA;	
		tc = -DTOA;	
		tc*=2*M_PI;
		//gen_params->tc = tc;
		fourier_detector_response_equatorial(frequencies, length,&wp, responses[i], gen_params->RA,gen_params->DEC,gen_params->psi, gen_params->gmst,(T *)NULL, gen_params->LISA_alpha0,gen_params->LISA_phi0, gen_params->theta_l, gen_params->phi_l, detectors[i]);
		for(int j = 0 ; j<length; j++){
			responses[i][j] *= exp(std::complex<T>(0,tc*frequencies[j]));
		}
	}
	wp.deallocate_memory();
	
	return;
}
template void create_coherent_GW_detection_reuse_WF<double>(std::string *, int, double*, int , gen_params_base<double> *,std::string, std::complex<double> **);
template void create_coherent_GW_detection_reuse_WF<adouble>(std::string *, int, adouble*, int , gen_params_base<adouble> *,std::string, std::complex<adouble> **);

/*! \brief Utility to calculate the snr of a fourier transformed data stream while maximizing over the coalescence parameters phic and tc
 *
 * The gen_params structure holds the parameters for the template to be used (the maximimum likelihood parameters)
 */
double data_snr_maximized_extrinsic(double *frequencies, /**< Frequencies used by data*/
	int length,/**< length of the data*/
	std::complex<double> *data,/**< input data in the fourier domain*/
	double *psd,/**< PSD for the detector that created the data*/
	std::string detector,/**< Name of the detector --See noise_util for options */
	std::string generation_method,/**< Generation method for the template -- See waveform_generation.cpp for options*/
	gen_params *param/**< gen_params structure for the template*/
	)
{
	
	/*produce noise curve*/
	double *noise = (double *)malloc(sizeof(double)*length);	
	//populate_noise(frequencies,detector, noise, length);
	for(int i = 0; i<length; i++)
		//noise[i] = noise[i]*noise[i];
		noise[i ] =psd[i];

	std::complex<double> *detector_response
			 = (std::complex<double> *)malloc(sizeof(std::complex<double>)*length);
	double *integrand
			 = (double *)malloc(sizeof(double)*length);
	/*produce the waveform*/
	fourier_detector_response(frequencies, length, detector_response, detector,generation_method, param);

	/*Calculate the template snr integrand 4*Re(h* h /S(f)) - factor of q for the plus, cross modes 
 * 	effect on the detector*/
	for (int i = 0; i<length;i++)
		integrand[i] = 4.*real(conj(detector_response[i])*detector_response[i])/noise[i]; 
	double delta_f = frequencies[1]-frequencies[0];
	double snr_template;
	snr_template = sqrt(simpsons_sum(delta_f,length, integrand));

	fftw_complex *in, *out; 
	fftw_plan p;
	in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * length);
        out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * length);
        p = fftw_plan_dft_1d(length, in, out,FFTW_FORWARD, FFTW_MEASURE);
	std::complex<double> g_tilde;
        for (int i=0;i<length; i++)
        {
                g_tilde = 4.*conj(data[i]) * detector_response[i] / noise[i];
                in[i][0] = real(g_tilde);
                in[i][1] = imag(g_tilde);
        }

        double *g = (double *)malloc(sizeof(double) * length);

        fftw_execute(p);

        for (int i=0;i<length; i++)
        {
                g[i] = std::abs(std::complex<double>(out[i][0],out[i][1])) ;
        }

        double max = *std::max_element(g, g+length)*delta_f;
	


	fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);
	fftw_cleanup();

	free(g);
	free(noise);
	free(detector_response);
	free(integrand);
	return max/(snr_template);
}
/*! \brief Light wrapper for the data_snr_maximized_extrinsic method
 *
 * Splits the data into real and imaginary, so all the arguments are C-safe
 */
double data_snr_maximized_extrinsic(double *frequencies, /**< Frequencies used by data*/
	int length,/**< length of the data*/
	double *data_real,/**< input data in the fourier domain -- real part*/
	double *data_imag,/**< input data in the fourier domain -- imaginary part*/
	double *psd,/**< PSD for the detector that created the data*/
	std::string detector,/**< Name of the detector --See noise_util for options */
	std::string generation_method,/**< Generation method for the template -- See waveform_generation.cpp for options*/
	gen_params *param/**< gen_params structure for the template*/
	)
{
	std::complex<double> *data = (std::complex<double> *)malloc(sizeof(std::complex<double>) * length);
        for (int i =0; i<length; i++)
                data[i] = std::complex<double>(data_real[i],data_imag[i]);	
	double snr;
	snr = data_snr_maximized_extrinsic(frequencies,
				length,
				data,
				psd,
				detector,
				generation_method,
				param);
	free(data);
	return snr;
}

double calculate_snr(std::string sensitivity_curve,
	std::string detector,
	std::string generation_method,
	gen_params_base<double> *params,
	double *frequencies,
	int length,
	std::string integration_method,
	double *weights,
	bool log10_freq)
{
	double *times;
	if(params->equatorial_orientation){
		transform_orientation_coords(params, generation_method,detector);
	}

	bool SA_save=params->sky_average;
	if(detector == "LISA" && !params->sky_average){
		times = new double[length];
		if(integration_method == "GAUSSLEG"){
			time_phase_corrected_autodiff(times, length, frequencies, params, generation_method, false, NULL,1);
		}
		else{
			time_phase_corrected(times, length, frequencies, params, generation_method, false,1);
		}
	}
	double snr;
	if(!params->sky_average){
		std::complex<double> *response = new std::complex<double>[length];
		fourier_detector_response(frequencies, length, response, detector, generation_method, params,times);
		snr = calculate_snr(sensitivity_curve, response, frequencies, length,integration_method, weights, log10_freq);
		delete [] response;	
	}
	else{
		if(sensitivity_curve.find("SADC") != std::string::npos && params->sky_average){
			params->sky_average=false;//We don't want the sky averaging factor typically used for terrestrial detectors, so we need to turn this off
		}
		//fourier_waveform(frequencies, length, hp,hc, generation_method, params);
		waveform_polarizations<double> wp;
		assign_polarizations(generation_method, &wp);
		wp.allocate_memory(length);
			
		fourier_waveform(frequencies, length, &wp, generation_method, params);
		snr = calculate_snr(sensitivity_curve, wp.hplus, frequencies, length,integration_method, weights, log10_freq);
		wp.deallocate_memory();

	}
	params->sky_average=SA_save;
	if(detector == "LISA"){
		if(!params->sky_average){
			delete [] times;
		}
	}
	return snr;

}
/**< \brief Routine to calculate the SNR of a template with GSL quadrature integration
 *
 * Sometimes, this is faster than the ``grid'' style integration
 *
 * Supports sky-averaged templates, but this should only be used with non-precessing waveforms
 *
 */
int calculate_snr_gsl(double *snr,/**< [out] SNR*/
	std::string sensitivity_curve,/**< Noise curve */
	std::string detector,/**<Detector to compute response -- can be empty is SA*/
	std::string generation_method,/**<Generation method */
	gen_params_base<double> *params,/**< Source Parameters*/
	double f_min,/**< Lower frequency bound*/
	double f_max,/**< Upper frequency bound*/
	double relative_error/**< Relative error threshold*/
	)
{
	if(params->equatorial_orientation){
		transform_orientation_coords(params, generation_method,detector);
	}
	int np=100;
	gsl_integration_workspace *w=gsl_integration_workspace_alloc(np) ;
	int err =  calculate_snr_gsl(snr,sensitivity_curve, detector, generation_method, params, f_min, f_max,relative_error,w,np);
	gsl_integration_workspace_free(w);
	return err;
}
/**< \brief Routine to calculate the SNR of a template with GSL quadrature integration
 *
 * Sometimes, this is faster than the ``grid'' style integration
 *
 * Supports sky-averaged templates, but this should only be used with non-precessing waveforms
 *
 */
int calculate_snr_gsl(double *snr,
	std::string sensitivity_curve,/**< Noise curve */
	std::string detector,/**<Detector to compute response -- can be empty is SA*/
	std::string generation_method,/**<Generation method */
	gen_params_base<double> *params,/**< Source Parameters*/
	double f_min,/**< Lower frequency bound*/
	double f_max,/**< Upper frequency bound*/
	double relative_error,/**< Relative error threshold*/
	gsl_integration_workspace *w, /**< User-allocated gsl_integration_workspace*/
	int np/**<Size of gsl_integration_workspace allocation*/
	)
{
	bool SA_save=params->sky_average;
	if(sensitivity_curve.find("SADC") != std::string::npos && params->sky_average){
		params->sky_average=false;	
	}
	gsl_snr_struct helper_params;
	helper_params.params = params;
	helper_params.generation_method = generation_method;
	helper_params.SN = sensitivity_curve;
	helper_params.detector = detector;
	gsl_function F;
	if(SA_save){
		F.function = [](double f, void *param){return integrand_snr_SA_subroutine(f,param);};
	}
	else{
		F.function = [](double f, void *param){return integrand_snr_subroutine(f,param);};
	}
	F.params = (void *)&helper_params;
	double result, err;
	gsl_set_error_handler_off();
	int errcode = gsl_integration_qag(&F,f_min, f_max, 0,relative_error, np, GSL_INTEG_GAUSS15,w, &result, &err);
	*snr = sqrt(result);
	params->sky_average = SA_save;
	return errcode;
}

/*! \brief Internal function to calculate the SNR integrand for sky-averaged waveforms
 *
 * NOTE: Only works for GR spin aligned 
 */
double integrand_snr_SA_subroutine(double f, void *subroutine_params)
{
	gsl_snr_struct cast_params = *(gsl_snr_struct *)subroutine_params;
	//fourier_waveform(&f, 1,&wfp, &wfc, cast_params.generation_method, cast_params.params);
	waveform_polarizations<double> wp;
	assign_polarizations(cast_params.generation_method, &wp);
	wp.allocate_memory(1);
	fourier_waveform(&f, 1,&wp, cast_params.generation_method, cast_params.params);
	double SN;
	populate_noise(&f, cast_params.SN, &SN,1);
	SN*=SN;
	double retval =  4*std::real(std::conj(wp.hplus[0])*wp.hplus[0])/SN;
	//std::cout<<retval<<" "<<SN<<" "<<wp.hplus[0]<<" "<<f<<std::endl;
	wp.deallocate_memory();
	return retval;
}
/*! \brief Internal function to calculate the SNR integrand for full waveforms
 */
double integrand_snr_subroutine(double f, void *subroutine_params)
{
	gsl_snr_struct cast_params = *(gsl_snr_struct *)subroutine_params;
	double times[2];
	if(cast_params.detector == "LISA" ){
		//double temp_f[2] = {f, f+1.e-4};
		//time_phase_corrected(times, 2, temp_f, cast_params.params, cast_params.generation_method, false);
		time_phase_corrected_autodiff(&times[0], 1, &f, cast_params.params, cast_params.generation_method, false, NULL);
	}
	double time = times[0];

	std::complex<double> response;
	fourier_detector_response(&f, 1,&response, cast_params.detector,cast_params.generation_method, cast_params.params, &time);
	double SN;
	populate_noise(&f, cast_params.SN, &SN,1);
	SN*=SN;
	return 4*std::real(std::conj(response)*response)/SN;
}

/*! \brief Caclulates the snr given a detector and waveform (complex) and frequencies
 *      
 * This function computes the un-normalized snr: \sqrt( ( H | H ) )
 */     
double calculate_snr(std::string sensitivity_curve, /**< detector name - must match the string of populate_noise precisely*/
        std::complex<double> *waveform,/**< complex waveform */
        double *frequencies,/**< double array of frequencies that the waveform is evaluated at*/
        int length,/**< length of the above two arrays*/
	std::string integration_method,
	double *weights,
	bool log10_freq)
{
        double *noise = (double *)malloc(sizeof(double)*length);
        populate_noise(frequencies,sensitivity_curve, noise,  length);
        for (int i = 0; i< length; i++){
                noise[i] = noise[i]*noise[i];
	}
	double snr  = calculate_snr_internal(noise,waveform, frequencies,length, integration_method, weights, log10_freq);

	free(noise);
        return snr;
}

double calculate_snr_internal(double *psd, 
	std::complex<double>*waveform,
	double *frequencies, 
	int length,
	std::string integration_method,
	double *weights,
	bool log10_freq)
{
        double *integrand = (double *) malloc(sizeof(double)*length);
        for (int i = 0; i<length; i++){
                integrand[i] = 4.* real(conj(waveform[i])*waveform[i]/psd[i]);
		if(log10_freq){
			integrand[i]*= (frequencies[i]*LOG10);	
		}
	}
        double integral=0;
	if(integration_method == "SIMPSONS"){
		integral = simpsons_sum((frequencies[1]-frequencies[0]), length, integrand);
	}
	else if(integration_method == "GAUSSLEG"){
        	for (int i = 0; i<length; i++){
			integral+= weights[i]*integrand[i];
		
		}
	
	}
	else{
		std::cout<<"UNSUPPORTED INTEGRATION METHOD"<<std::endl;	
	}
	free(integrand);
        return sqrt(integral);
}

// SNR calculator with Quadrature method
double calculate_snr_internal(
	double *psd,
	std::complex<double> *waveform,
	const Quadrature *QuadMethod
)
{
	int length = QuadMethod->get_length();
	// Integrand array
	double *integrand = new double [length];

	for (int i = 0; i < length; i++)
	{
		integrand[i] = real(conj(waveform[i])*waveform[i]/psd[i]);
	}

	double innerproduct = QuadMethod->integrate(integrand);
	
	// Clean-up
	delete [] integrand;

	return 2.*sqrt(innerproduct);
}

/* \brief calculates the detector response for a given waveform and detector -- polarization angle =0
 */
template<class T>
int fourier_detector_response_horizon(T *frequencies, /**<array of frequencies corresponding to waveform*/
			int length,/**< length of frequency/waveform arrays*/
			//std::complex<T> *hplus, /*<precomputed plus polarization of the waveform*/ 
			//std::complex<T> *hcross, /**<precomputed cross polarization of the waveform*/ 
			waveform_polarizations<T> *wp,
			std::complex<T> *detector_response, /**< [out] detector response*/
			T theta, /**< polar angle (rad) theta in detector frame*/
			T phi, /**< azimuthal angle (rad) phi in detector frame*/ 
			std::string detector/**< detector - list of supported detectors in noise_util*/
			)
{
	//return fourier_detector_response_horizon(frequencies, length, hplus, hcross, detector_response, theta, phi, (T)0., detector);
	return fourier_detector_response_horizon(frequencies, length, wp, detector_response, theta, phi, (T)0., detector);
	
}
/* \brief calculates the detector response for a given waveform and detector
 */
template<class T>
int fourier_detector_response_horizon(T *frequencies, /**<array of frequencies corresponding to waveform*/
			int length,/**< length of frequency/waveform arrays*/
			//std::complex<T> *hplus, /*<precomputed plus polarization of the waveform*/ 
			//std::complex<T> *hcross, /**<precomputed cross polarization of the waveform*/ 
			waveform_polarizations<T> *wp,
			std::complex<T> *detector_response, /**< [out] detector response*/
			T theta, /**< polar angle (rad) theta in detector frame*/
			T phi, /**< azimuthal angle (rad) phi in detector frame*/ 
			T psi, /**< polarization angle (rad) phi in detector frame*/ 
			std::string detector/**< detector - list of supported detectors in noise_util*/
			)
{
	int status=1;
	bool extra_polarizations = false;
	for(int i = 2 ; i<6; i++){
		if(wp->active_polarizations[i]){
			extra_polarizations=true;
		}
	}
	T fplus, fcross, Fplus, Fcross, c2psi, s2psi,Fx,Fy,Fb,Fl;
	
	if(	detector == "LIGO" || 
		detector == "Livingston" || 
		detector == "LIVINGSTON" || 
		detector == "livingston" || 
		detector == "Hanford" || 
		detector == "HANFORD" || 
		detector == "hanford" || 
		detector == "VIRGO" ||
		detector == "Virgo" ||
		detector == "virgo" ||
		detector == "CE" ||
		detector == "ET1"||
		detector == "ET2"||
		detector == "ET3"
	)
	{
		//fplus = right_interferometer_plus(theta,phi);
		//fcross = right_interferometer_cross(theta,phi);	
		//c2psi = cos(2*psi);
		//s2psi = sin(2*psi);
		//Fplus = fplus*c2psi- fcross*s2psi;
		//Fcross = fplus*s2psi+ fcross*c2psi;
		det_res_pat<T> r_pat;
		r_pat.Fplus = &Fplus;
		r_pat.Fcross = &Fcross;
		r_pat.Fx = &Fx;
		r_pat.Fy = &Fy;
		r_pat.Fb = &Fb;
		r_pat.Fl = &Fl;
		r_pat.active_polarizations = &(wp->active_polarizations[0]);
		right_interferometer(&r_pat, theta, phi,psi);
	}
	if ( 	detector == "ET1"||
		detector == "ET2"||
		detector == "ET3")
	{
		Fplus*= ET1_geometric_factor;	
		Fcross*= ET1_geometric_factor;	
		if(extra_polarizations){
			Fx*= ET1_geometric_factor;	
			Fy*= ET1_geometric_factor;	
			Fb*= ET1_geometric_factor;	
			Fl*= ET1_geometric_factor;	
		}
	}
	for (int i =0; i <length; i++)
	{
		detector_response[i] = Fplus * wp->hplus[i] 
					+ (Fcross )*wp->hcross[i];
	}	
	if(extra_polarizations){

		if(wp->active_polarizations[2]){
			for (int i =0; i <length; i++)
			{
				detector_response[i] += Fx * wp->hx[i] ;
			}	
		}
		if(wp->active_polarizations[3]){
			for (int i =0; i <length; i++)
			{
				detector_response[i] += Fy * wp->hy[i] ;
			}	
		}
		if(wp->active_polarizations[4]){
			for (int i =0; i <length; i++)
			{
				detector_response[i] += Fb * wp->hb[i] ;
			}	
		}
		if(wp->active_polarizations[5]){
			for (int i =0; i <length; i++)
			{
				detector_response[i] += Fl * wp->hl[i] ;
			}	
		}
	}
	return status;
	
}
template<class T>
int time_detector_response_horizon(T *times, /**< double array of frequencies for the waveform to be evaluated at*/
	int length,/**<integer length of all the arrays*/
	std::complex<T> *response, /**< [out] complex array for the output plus polarization waveform*/
	std::string detector,
	std::string generation_method,/**<String that corresponds to the generation method - MUST BE SPELLED EXACTLY*/
	gen_params_base<T> *parameters/**<structure containing all the source parameters*/
	)
{
	int status = 1;
	//generate waveform
	waveform_polarizations<T> wp;
	assign_polarizations(generation_method, &wp);
	wp.allocate_memory(length);

	status = time_waveform(times, 
			length,
			&wp,
			generation_method,
			parameters
			);
	//Not sure why horizon needs frequencies?? 
	status = fourier_detector_response_horizon((T *)NULL, 
			length, 
			&wp,
			response, 
			parameters->theta, 
			parameters->phi, 
			parameters->psi, 
			detector 
			) ;
	
		
	//Deallocate memory
	wp.deallocate_memory();

	return status;
}

/*!\brief Function to produce the detector response caused by impinging gravitational waves from a quasi-circular binary
 *
 * By using the structure parameter, the function is allowed to be more flexible in using different 
 * method of waveform generation - not all methods use the same parameters
 *
 * This puts the responsibility on the user to pass the necessary parameters
 *
 * Detector options include classic interferometers like LIGO/VIRGO (coming soon: ET and LISA)
 * 	
 * This is a wrapper that combines generation with response functions: if producing mulitple responses for one waveform (ie stacking Hanford, Livingston, and VIRGO), it will be considerably more efficient to calculate the waveform once, then combine each response manually 
 */
template<class T>
int fourier_detector_response_horizon(T *frequencies, /**< double array of frequencies for the waveform to be evaluated at*/
	int length,/**<integer length of all the arrays*/
	std::complex<T> *response, /**< [out] complex array for the output plus polarization waveform*/
	std::string detector,
	std::string generation_method,/**<String that corresponds to the generation method - MUST BE SPELLED EXACTLY*/
	gen_params_base<T> *parameters/**<structure containing all the source parameters*/
	)
{
	int status = 1;
	//generate waveform
	waveform_polarizations<T> wp;
	assign_polarizations(generation_method, &wp);
	wp.allocate_memory(length);

	status = fourier_waveform(frequencies, 
			length,
			&wp,
			generation_method,
			parameters
			);
	status = fourier_detector_response_horizon(frequencies, 
			length, 
			&wp,
			response, 
			parameters->theta, 
			parameters->phi, 
			parameters->psi, 
			detector 
			) ;
	
		
	//Deallocate memory
	wp.deallocate_memory();

	return status;
}

/* \brief calculates the detector response for a given waveform and detector -- using equatorial coordinates and greenwich mean sidereal time
 */
template<class T>
int fourier_detector_response_equatorial(T *frequencies, /**<array of frequencies corresponding to waveform*/
			int length,/**< length of frequency/waveform arrays*/
			//std::complex<T> *hplus, /*<precomputed plus polarization of the waveform*/ 
			//std::complex<T> *hcross, /**<precomputed cross polarization of the waveform*/ 
			waveform_polarizations<T> *wp,
			std::complex<T> *detector_response, /**< [out] detector response*/
			T ra, /**< Right Ascension in rad*/
			T dec, /**< Declination in rad*/ 
			T psi, /**< polarization angle (rad) */ 
			double gmst, /**< greenwich mean sidereal time*/ 
			T *times,
			T LISA_alpha0,
			T LISA_phi0,
			T theta_j_ecl,
			T phi_j_ecl,
			std::string detector/**< detector - list of supported detectors in noise_util*/
			)
{
	int status=1;
	//Check for extra polarizations
	bool extra_polarizations = false;
	for(int i = 2 ; i<6; i++){
		if(wp->active_polarizations[i]){
			extra_polarizations=true;
		}
	}
	//Not an elegant solution, but should work..
	T fplus ;
	T fcross ;
	T *Fplus;
	T *Fcross;
	T *Fx,*Fy,*Fb,*Fl;
	T fx,fy,fb,fl;
	det_res_pat<T> r_pat;
	r_pat.active_polarizations = &(wp->active_polarizations[0]);
	if(detector=="LISA"){
		Fplus = new T[length];
		Fcross = new T[length];
		r_pat.Fplus = Fplus;
		r_pat.Fcross = Fcross;
		if(extra_polarizations){
			if(wp->active_polarizations[2]){
				Fx = new T[length];
				r_pat.Fx = Fx;
			}
			if(wp->active_polarizations[3]){
				Fy = new T[length];
				r_pat.Fy = Fy;
			}
			if(wp->active_polarizations[4]){
				Fb = new T[length];
				r_pat.Fb = Fb;
			}
			if(wp->active_polarizations[5]){
				Fl = new T[length];
				r_pat.Fl = Fl;
			}
		}
		//detector_response_functions_equatorial(detector, ra, dec, psi, gmst,times, length,LISA_alpha0, LISA_phi0,theta_j_ecl,phi_j_ecl, Fplus, Fcross);
		detector_response_functions_equatorial(detector, ra, dec, psi, gmst,times, length,LISA_alpha0, LISA_phi0,theta_j_ecl,phi_j_ecl, &r_pat);
	}
	else{
		r_pat.Fplus = &fplus;
		r_pat.Fcross = &fcross;
		if(extra_polarizations){
			if(wp->active_polarizations[2]){
				r_pat.Fx = &fx;
			}
			if(wp->active_polarizations[3]){
				r_pat.Fy = &fy;
			}
			if(wp->active_polarizations[4]){
				r_pat.Fb = &fb;
			}
			if(wp->active_polarizations[5]){
				r_pat.Fl = &fl;
			}
		}
		detector_response_functions_equatorial(detector, ra, dec, psi, gmst,times, length,LISA_alpha0, LISA_phi0,theta_j_ecl,phi_j_ecl, &r_pat);
	}
	
	if(detector == "LISA"){
		for (int i =0; i <length; i++)
		{
			//Doppler phase shift

			//detector_response[i] = (r_pat.Fplus[i] * wp->hplus[i] 
			//			+ (r_pat.Fcross[i] )*wp->hcross[i]) * std::exp(std::complex<T>(0, - doppler_phase_shift ) );
			detector_response[i] = (r_pat.Fplus[i] * wp->hplus[i] 
						+ (r_pat.Fcross[i] )*wp->hcross[i]) ;
		}	
		if(extra_polarizations){
			if(wp->active_polarizations[2]){
				for(int i = 0 ; i<length; i++){
					detector_response[i] += r_pat.Fx[i] * wp->hx[i];
				}
			}
			if(wp->active_polarizations[3]){
				for(int i = 0 ; i<length; i++){
					detector_response[i] += r_pat.Fy[i] * wp->hy[i];
				}
			}
			if(wp->active_polarizations[4]){
				for(int i = 0 ; i<length; i++){
					detector_response[i] += r_pat.Fb[i] * wp->hb[i];
				}
			}
			if(wp->active_polarizations[5]){
				for(int i = 0 ; i<length; i++){
					detector_response[i] += r_pat.Fl[i] * wp->hl[i];
				}
			}
		}
		for (int i =0; i <length; i++)
		{
			T theta_s, phi_s, phi_t;
			ecl_from_eq((T)(M_PI/2. - dec), ra, &theta_s, &phi_s);
			phi_t = LISA_phi0 + 2. * M_PI * times[i] / T_year;
			T doppler_phase_shift = 2*M_PI * frequencies[i] * AU_SEC *sin(theta_s) *cos(phi_t - phi_s);
			detector_response[i]*= std::exp(std::complex<T>(0, - doppler_phase_shift ) );
		}
		delete [] Fplus; 
		delete [] Fcross;
		if(extra_polarizations){
			if(wp->active_polarizations[2]){
				delete [] Fx;
			}
			if(wp->active_polarizations[3]){
				delete [] Fy;
			}
			if(wp->active_polarizations[4]){
				delete [] Fb;
			}
			if(wp->active_polarizations[5]){
				delete [] Fl;
			}
		}

	}
	else{
		for (int i =0; i <length; i++)
		{
			//detector_response[i] = fplus * hplus[i] 
			//			+ (fcross )*hcross[i];
			//detector_response[i] = fplus * hplus[i] 
			//			+ (fcross )*hcross[i];
			detector_response[i] = (*(r_pat.Fplus)) * wp->hplus[i] 
						+ (*(r_pat.Fcross ))*wp->hcross[i];
		}	
		if(extra_polarizations){
			if(wp->active_polarizations[2]){
				for(int i = 0 ; i<length; i++){
					detector_response[i] += (*(r_pat.Fx)) * wp->hx[i];
				}
			}
			if(wp->active_polarizations[3]){
				for(int i = 0 ; i<length; i++){
					detector_response[i] += (*(r_pat.Fy)) * wp->hy[i];
				}
			}
			if(wp->active_polarizations[4]){
				for(int i = 0 ; i<length; i++){
					detector_response[i] += (*(r_pat.Fb)) * wp->hb[i];
				}
			}
			if(wp->active_polarizations[5]){
				for(int i = 0 ; i<length; i++){
					detector_response[i] += (*(r_pat.Fl)) * wp->hl[i];
				}
			}
		}
	}
	return status;
	
}
/*!\brief Function to produce the detector response caused by impinging gravitational waves from a quasi-circular binary for equatorial coordinates for TERRESTIAL detectors, where the earth's rotation during the signal is minor
 *
 * By using the structure parameter, the function is allowed to be more flexible in using different 
 * method of waveform generation - not all methods use the same parameters
 *
 * This puts the responsibility on the user to pass the necessary parameters
 *
 * Detector options include classic interferometers like LIGO/VIRGO (coming soon: ET and LISA)
 * 	
 * This is a wrapper that combines generation with response functions: if producing mulitple responses for one waveform (ie stacking Hanford, Livingston, and VIRGO), it will be considerably more efficient to calculate the waveform once, then combine each response manually 
 */
template<class T>
int fourier_detector_response_equatorial(T *frequencies, /**< double array of frequencies for the waveform to be evaluated at*/ 
			int length,/**<integer length of all the arrays*/
			std::complex<T> *response, /**< [out] complex array for the output plus polarization waveform*/
			std::string detector,
			std::string generation_method,/**<String that corresponds to the generation method - MUST BE SPELLED EXACTLY*/
			gen_params_base<T> *parameters/**<structure containing all the source parameters*/
			)
{
	T *times=NULL;
	return fourier_detector_response_equatorial(frequencies, length,response, detector, generation_method, parameters, times);	

}
/*!\brief Function to produce the detector response caused by impinging gravitational waves from a quasi-circular binary for equatorial coordinates
 *
 * By using the structure parameter, the function is allowed to be more flexible in using different 
 * method of waveform generation - not all methods use the same parameters
 *
 * This puts the responsibility on the user to pass the necessary parameters
 *
 * Detector options include classic interferometers like LIGO/VIRGO (coming soon: ET and LISA)
 * 	
 * This is a wrapper that combines generation with response functions: if producing mulitple responses for one waveform (ie stacking Hanford, Livingston, and VIRGO), it will be considerably more efficient to calculate the waveform once, then combine each response manually 
 */
template<class T>
int fourier_detector_response_equatorial(T *frequencies, /**< double array of frequencies for the waveform to be evaluated at*/ 
			int length,/**<integer length of all the arrays*/
			std::complex<T> *response, /**< [out] complex array for the output plus polarization waveform*/
			std::string detector,
			std::string generation_method,/**<String that corresponds to the generation method - MUST BE SPELLED EXACTLY*/
			gen_params_base<T> *parameters,/**<structure containing all the source parameters*/
			T *times
			)
{
	int status = 1;
	//generate waveform
	if(parameters->equatorial_orientation){
		transform_orientation_coords(parameters, generation_method,detector);
	}
	else{
		if(detector=="LISA"){
			std::cout<<"ERROR -- fourier_detector_response_equatorial -- LISA currently only accepts equatorial_orientation parameters"<<std::endl;
		}
	}
	waveform_polarizations<T> wp;
	assign_polarizations(generation_method, &wp);
	wp.allocate_memory(length);
	status = fourier_waveform(frequencies, 
			length,
			&wp,	
			generation_method,
			parameters
			);
	status = fourier_detector_response_equatorial(frequencies, 
			length, 
			&wp,
			response, 
			parameters->RA, 
			parameters->DEC, 
			parameters->psi, 
			parameters->gmst, 
			times,
			parameters->LISA_alpha0,
			parameters->LISA_phi0,
			parameters->theta_j_ecl,
			parameters->phi_j_ecl,
			detector 
			) ;
	
		
	//Deallocate memory
	wp.deallocate_memory();

	return status;
}

template<class T>
int time_detector_response_equatorial(T *times, /**< double array of frequencies for the waveform to be evaluated at*/ 
			int length,/**<integer length of all the arrays*/
			std::complex<T> *response, /**< [out] complex array for the output plus polarization waveform*/
			std::string detector,
			std::string generation_method,/**<String that corresponds to the generation method - MUST BE SPELLED EXACTLY*/
			gen_params_base<T> *parameters/**<structure containing all the source parameters*/
			)
{
	int status = 1;
	if(detector=="LISA"){
		debugger_print(__FILE__,__LINE__,"LISA is not supported for TD waveforms because of doppler shift issues");
		return 1;
	}
	//generate waveform
	if(parameters->equatorial_orientation){
		transform_orientation_coords(parameters, generation_method,detector);
	}
	else{
		if(detector=="LISA"){
			std::cout<<"ERROR -- fourier_detector_response_equatorial -- LISA currently only accepts equatorial_orientation parameters"<<std::endl;
		}
	}
	waveform_polarizations<T> wp;
	assign_polarizations(generation_method, &wp);
	wp.allocate_memory(length);
	status = time_waveform(times, 
			length,
			&wp,	
			generation_method,
			parameters
			);
	//Nothing changes in fourier_detector_response_equatorial for time waveforms except for 
	//LISA's doppler shift
	status = fourier_detector_response_equatorial((T *)NULL, 
			length, 
			&wp,
			response, 
			parameters->RA, 
			parameters->DEC, 
			parameters->psi, 
			parameters->gmst, 
			times,
			parameters->LISA_alpha0,
			parameters->LISA_phi0,
			parameters->theta_j_ecl,
			parameters->phi_j_ecl,
			detector 
			) ;
	
		
	//Deallocate memory
	wp.deallocate_memory();

	return status;
}
template<class T>
int time_detector_response(T *times,
	int length,
	std::complex<T> *response,
	std::string detector,
	std::string generation_method,
	gen_params_base<T> *parameters
	)
{
	int status;
	if(parameters->horizon_coord){
		status=time_detector_response_horizon(times, length, response, detector, generation_method, parameters);
	}
	else{
		status=time_detector_response_equatorial(times, length, response, detector, generation_method, parameters);

	}
	return status;

}
template int time_detector_response<double>(double *, int, std::complex<double> *, std::string, std::string, gen_params_base<double> *);
template int time_detector_response<adouble>(adouble *, int, std::complex<adouble> *, std::string, std::string, gen_params_base<adouble> *);

/*! \brief Wrapper to handle all detector_response calls -- horizon and equatorial
 *
 */
template<class T>
int fourier_detector_response(T *frequencies,
	int length,
	std::complex<T> *response,
	std::string detector,
	std::string generation_method,
	gen_params_base<T> *parameters,
	T *times)
{
	int status;
	if(parameters->horizon_coord){
		status=fourier_detector_response_horizon(frequencies, length, response, detector, generation_method, parameters);
	}
	else{
		status=fourier_detector_response_equatorial(frequencies, length, response, detector, generation_method, parameters, times);

	}
	return status;

}
template int fourier_detector_response<double>(double *, int, std::complex<double> *, std::string, std::string, gen_params_base<double> *, double *);
template int fourier_detector_response<adouble>(adouble *, int, std::complex<adouble> *, std::string, std::string, gen_params_base<adouble> *, adouble *);

/*! \brief Calculates the amplitude (magnitude) and phase (argument) of the response of a given detector
 *
 * This is for general waveforms, and will work for precessing waveforms
 *
 * Not as fast as non-precessing, but that can't be helped. MUST include plus/cross polarizations
 */
int fourier_detector_amplitude_phase(double *frequencies, 
			int length,
			double *amplitude, 
			double *phase, 
			std::string detector,
			std::string generation_method,
			gen_params *parameters
			)
{
	std::complex<double> *response = (std::complex<double>*)malloc(
			sizeof(std::complex<double> ) * length);	
	fourier_detector_response_horizon(frequencies, length,
				response, detector, generation_method, parameters);
	for(int i =0 ; i< length; i++){
		amplitude[i] = std::abs(response[i]);
		phase[i] = std::arg(response[i]);
	}

	free(response);
	return 0;
}

/*! \brief Mapping from phase to time AUTODIFFERENTIATION using ADOL-C
 *
 *
 * ***NOTE*** Currently, only IMRPhenomD is used for the time-phase relationship, as this is simpler and captures almost all of the physics. The generation_method parameter does NOTHING right now.
 *
 * This is NOT autodiff safe. ADOL-C does not support wrapping a section of code as active twice. To find the derivative of this function, you must use the hessian function of ADOL-C
 *
 * Made for use with detectors like LISA
 *
 * Uses the deprecated postmerger calculation for IMRPhenomD -- this is ADOL-C friendly
 *
 * Using https://arxiv.org/abs/1809.04799 t = (1/2PI) d phi/ d f
 *
 * correct_time is currently not supported
 *
 * For IMRPhenomPv2, the phase has to be wrapped, because arctan is taken of the waveform because of the euler rotations. This might make the numerical derivative unpredictable
 *
 * if the order is set to 2, returns (1/2PI) d^2 phi / d f^2
 *
 */
void time_phase_corrected_autodiff(double *times, int length, double *frequencies,gen_params_base<double> *params, std::string generation_method, bool correct_time, int *tapes_in,int order)
{
	std::string local_method = "IMRPhenomD";
	if(generation_method.find("ppE") != std::string::npos){
		if(generation_method.find("Inspiral") != std::string::npos){
			local_method = "ppE_IMRPhenomD_Inspiral";
		}
		else if(generation_method.find("IMR") != std::string::npos){
			local_method = "ppE_IMRPhenomD_IMR";
		}
	}
	//std::string local_method = generation_method;

	bool save_dep = params->dep_postmerger;
	params->dep_postmerger=true;

	int boundary_num = boundary_number(generation_method);
	double freq_boundaries[boundary_num];
	double grad_freqs[boundary_num];
	assign_freq_boundaries(freq_boundaries, grad_freqs, boundary_num, params, local_method);	
	int *tapes;
	gen_params_base<adouble> aparams;
	if(tapes_in){ 
		tapes = tapes_in;
	}
	else{
		tapes = new int[boundary_num];
		transform_parameters(params, &aparams);
		if(params->equatorial_orientation){
			transform_orientation_coords(&aparams,local_method,"");
		}
		for(int i = 0 ; i<boundary_num ; i++){
			tapes[i]=(i+1)*8;	
			trace_on(tapes[i]);
			adouble freq;
			freq <<= grad_freqs[i];
			adouble phasep, phasec;
			fourier_phase(&freq, 1, &phasep, &phasec, local_method, &aparams);
			double phaseout;
			phasep>>=phaseout;
			trace_off();
		}
	}
	bool eval = false;
	double freq;	
	double **temp;
	if(order==2){
		temp = new double*[1];
		temp[0]=new double[1];
	}
	for(int k = 0; k<length; k++){
		freq = frequencies[k];
		for(int n = 0; n<boundary_num ; n++){
			if(freq < freq_boundaries[n]){
				if(order==1){
					gradient(tapes[n], 1, &freq, &times[k]);
				}
				else if(order==2){	
					hessian(tapes[n], 1, &freq, temp);
					times[k]=temp[0][0];
				}
				//Mark successful derivative
				eval = true;
				//Skip the rest of the bins
				break;
			}
		}
		if(!eval){
			//std::cout<<"Not EVAL "<< freq<<std::endl;
			//New CHANGE 05-21-2023
			//times[k]=0;
			if(k !=0 && order==1){
				times[k]=times[k-1]+(-frequencies[k]*frequencies[k]+frequencies[k-1]*frequencies[k-1] );
			}
			else if(k !=0 && order==2){
				times[k]=times[k-1]+2*(-frequencies[k]+frequencies[k-1] );
			}
			else{
				source_parameters<double> s_param;
				//s_param = source_parameters<T>::populate_source_parameters(params);
				s_param.populate_source_parameters(params);
				s_param.sky_average = params->sky_average;
				s_param.f_ref = params->f_ref;
				s_param.phiRef = params->phiRef;
				s_param.cosmology=params->cosmology;
				s_param.incl_angle=params->incl_angle;
				s_param.dep_postmerger=params->dep_postmerger;
				IMRPhenomD<double> model;
				lambda_parameters<double> lambda;
				model.assign_lambda_param(&s_param,&lambda);	
				model.post_merger_variables(&s_param);
				double fRD = s_param.fRD;
				double fdamp = s_param.fdamp;
				double  fpeak = model.fpeak(&s_param , &lambda);
				if(order==1){
					times[k]=params->tc+(-frequencies[k]*frequencies[k]+fpeak*fpeak );
						
				}
				else if(order==2){
					times[k]=params->tc+2*(-frequencies[k]+fpeak );
				}
				else{
					times[k]=0;
				}
			}
				
		}
		eval = false;
	}
	if(order==2){
		delete [] temp[0];
		delete [] temp;
	}

	//divide by 2 PI
	for(int i = 0 ; i<length; i++){
		times[i]/=(2.*M_PI);
	}
	if(!tapes_in){
		delete [] tapes;
		if(check_mod(generation_method)){
			if(generation_method.find("ppE")!= std::string::npos ||
			check_theory_support(generation_method)){
				delete [] aparams.betappe;
				delete [] aparams.bppe;
			}
			else if(generation_method.find("gIMR")!= std::string::npos){
				if(aparams.Nmod_phi !=0){
					delete [] aparams.phii;
					delete [] aparams.delta_phi;
				}
				if(aparams.Nmod_sigma !=0){
					delete [] aparams.sigmai;
					delete [] aparams.delta_sigma;
				}
				if(aparams.Nmod_beta !=0){
					delete [] aparams.betai;
					delete [] aparams.delta_beta;
				}
				if(aparams.Nmod_alpha !=0){
					delete [] aparams.alphai;
					delete [] aparams.delta_alpha;
				}
			}
		}
	}
	params->dep_postmerger=save_dep;
	//params->shift_time=save_shift_time;
	
}
/*! \brief Utility to inform the fisher routine how many logical boundaries should be expected
 *
 * The automatic derivative code requires a new tape for each logical branch of the program, so each waveform_generation method needs to add the number of branches here
 */
int boundary_number(std::string method)
{
	if(method.find("IMRPhenomP") != std::string::npos || 
		method.find("IMRPhenomD")!=std::string::npos){
		if(method.find("NRT") != std::string::npos){
			return 7;	
		}
		return 5;
	}
	return -1;
}
/*! \brief Mapping from phase to time NUMERICAL
 *
 * ***NOTE*** Currently, only IMRPhenomD is used for the time-phase relationship, as this is simpler and captures almost all of the physics. The generation_method parameter does NOTHING right now.
 *
 * Made for use with detectors like LISA
 *
 * Using https://arxiv.org/abs/1809.04799 t = (1/2PI) d phi/ d f
 *
 * correct_time is currently not supported
 *
 * For IMRPhenomPv2, the phase has to be wrapped, because arctan is taken of the waveform because of the euler rotations. This might make the numerical derivative unpredictable
 *
 * Just uses a second order numerical derivative for now
 *
 * Second order is only available for 1 frequency at a time -- something's funky with this
 */
template<class T>
void time_phase_corrected(T *times, int length, T *frequencies,gen_params_base<T> *params, std::string generation_method, bool correct_time, int order)
{
	std::string local_method="IMRPhenomD";
	if(generation_method.find("ppE") != std::string::npos){
		if(generation_method.find("Inspiral") != std::string::npos){
			local_method = "ppE_IMRPhenomD_Inspiral";
		}
		else if(generation_method.find("IMR") != std::string::npos){
			local_method = "ppE_IMRPhenomD_IMR";
		}
	}

	if(params->equatorial_orientation){
		transform_orientation_coords(params,local_method,"");
	}
	if(length ==1){
		T epsilon =frequencies[0]*1e-2;
		if(order == 1){
			T temp_f[2];
			temp_f[0] = frequencies[0]-epsilon;
			temp_f[1] = frequencies[0]+epsilon;
			if(temp_f[1] > 0.2/(params->mass1 + params->mass2)/MSOL_SEC){
				times[0] = 0 ;	
				return;
			}
			T temp_deltaf = 2*epsilon;
			//T temp_deltaf = temp_f[1]-temp_f[0];
			T temp_phase_plus[2];
			T temp_phase_cross[2];
			fourier_phase(temp_f, 2, temp_phase_plus, temp_phase_cross, local_method, params);
			times[0] = (temp_phase_plus[1]-temp_phase_plus[0])/(2*M_PI*temp_deltaf);
		}
		else if(order == 2){
			T temp_f[3];
			temp_f[0] = frequencies[0]-epsilon;
			temp_f[2] = frequencies[0]+epsilon;
			temp_f[1] = frequencies[0];
			//T temp_deltaf = epsilon;
			T temp_phase_plus[3];
			T temp_phase_cross[3];
			fourier_phase(temp_f, 3, temp_phase_plus, temp_phase_cross, local_method, params);
			times[0] = (temp_phase_plus[2]-2*temp_phase_plus[1] + temp_phase_plus[0])/(2*M_PI*epsilon*epsilon);

		}
		return ;
	}
	//bool save_shift_time = params->shift_time;
	//params->shift_time = false;
	//################################################
	T *phase_plus = new T[length];
	T *phase_cross = new T[length];
	//################################################
	//################################################
	fourier_phase(frequencies, length, phase_plus, phase_cross, local_method, params);
	//################################################
	T fdamp, fRD, fpeak, deltaf;
	if(local_method.find("IMRPhenomD")!=std::string::npos){
		source_parameters<T> s_param;
		//s_param = source_parameters<T>::populate_source_parameters(params);
		s_param.populate_source_parameters(params);
		s_param.sky_average = params->sky_average;
		s_param.f_ref = params->f_ref;
		s_param.phiRef = params->phiRef;
		s_param.cosmology=params->cosmology;
		s_param.incl_angle=params->incl_angle;
		s_param.dep_postmerger=params->dep_postmerger;
		IMRPhenomD<T> model;
		lambda_parameters<T> lambda;
		model.assign_lambda_param(&s_param,&lambda);	
		model.post_merger_variables(&s_param);
		fRD = s_param.fRD;
		fdamp = s_param.fdamp;
		fpeak = model.fpeak(&s_param , &lambda);
		deltaf = frequencies[1]-frequencies[0];
	}
	else if(local_method.find("IMRPhenomPv2")!=std::string::npos){
		source_parameters<T> s_param;
		//s_param = source_parameters<T>::populate_source_parameters(params);
		s_param.populate_source_parameters(params);
		s_param.sky_average = params->sky_average;
		s_param.f_ref = params->f_ref;
		s_param.phiRef = params->phiRef;
		s_param.cosmology=params->cosmology;
		s_param.incl_angle=params->incl_angle;
		s_param.dep_postmerger=params->dep_postmerger;
		IMRPhenomPv2<T> model;
		if( ( params->chip + 1) > DOUBLE_COMP_THRESH){
			s_param.chip = params->chip;
			s_param.phip = params->phip;
		}
		else{
			s_param.spin1y = params->spin1[1];
			s_param.spin2y = params->spin2[1];
			s_param.spin1x = params->spin1[0];
			s_param.spin2x = params->spin2[0];
			model.PhenomPv2_Param_Transform(&s_param);
			
		}
		lambda_parameters<T> lambda;
		model.assign_lambda_param(&s_param,&lambda);	
		model.post_merger_variables(&s_param);
		fRD = s_param.fRD;
		fdamp = s_param.fdamp;
		fpeak = model.fpeak(&s_param , &lambda);
		deltaf = frequencies[1]-frequencies[0];
	}
	//################################################
	//Factor of 2 pi for the definition of time from frequency
	//IMRPhenomD returns (-) phase
	if(local_method.find("IMRPhenom") !=std::string::npos){
		//Currently using Nico's fix
		if(correct_time){
			T f = frequencies[0];
			bool check = true, check2=true;
			T pt1, pt2, f1,f2;
			int i = 0 ;
			//One sided, to start it off
			times[0] = (phase_plus[1]-phase_plus[0])/(2.*M_PI*deltaf);
			i++;
			while(f < .95*fRD && i<length-1)
			{
				f = frequencies[i];
				//central difference for the rest of the steps
				times[i] = (phase_plus[i+1]-phase_plus[i-1])/(4.*M_PI*deltaf);
				if(check){
					if(f>fpeak){
						pt1 = times[i];
						f1 = f;
						check=false;
					}
				}
				else{
					if(check2){
						if(f>.9*fRD){
							pt2 = times[i];
							f2 = f;
							check2=false;
						}
					}
				}
				i++;
			}	
			T f_intercept = times[i-1];
			T f_mr = f;
			T freq_slope_pm = (pt2-pt1)/(f2-f1);
			while(f<1.5*fRD && i<length-1)
			{
				f = frequencies[i];
				times[i] =f_intercept+ (f-f_mr)*freq_slope_pm;
				//Stop if observation goes past 20 years
				i++;
				if(times[i-1]>(params->tc+630720000)){
					break;
				}

			}
			T time_transition = times[i-1];
			T f_transition = f;
			while(i<length-1)
			{
				f = frequencies[i];
				times[i] = time_transition +pow_int(f-f_transition,2);
				//Stop if observation goes past 20 years
				i++;
				if(times[i-1]>(params->tc+630720000)){
					break;
				}
				
			}
			//if stopped at 20 years, fill out with frozen time
			if(i != length){
				while(i<length){
					times[i] = times[i-1];
					i++;
				}
			}
			else{
				times[length-1] = (phase_plus[length-1] - phase_plus[length-2])/(2*M_PI*deltaf);
			}
		}
		else{
			times[0] = (phase_plus[1]-phase_plus[0])/(2*M_PI*deltaf);
			if(length>2){
				for(int i = 1  ;i<length-1; i++){
					times[i] = (phase_plus[i+1]-phase_plus[i-1])/(4*M_PI*deltaf);
				}
				times[length-1] = (phase_plus[length-1]-phase_plus[length-2])/(2*M_PI*deltaf);
			}
			else{
				//only option with length 2 vector
				times[1] = times[0];
			}
		}
	}
	//################################################
	delete [] phase_plus;
	delete [] phase_cross;
	//params->shift_time = save_shift_time;
}

template<class T>
void transform_orientation_coords(gen_params_base<T> *parameters,std::string generation_method,std::string detector)
{
	if(generation_method.find("IMRPhenomP") != std::string::npos){
		T theta_s = M_PI/2. - parameters->DEC;
		T phi_s = parameters->RA;
		//N should point source to detector
		T Neq[3] = {-sin(theta_s)*cos(phi_s), -sin(theta_s)*sin(phi_s), -cos(theta_s)};
		T Leq[3] = {sin(parameters->theta_l)*cos(parameters->phi_l), 
			sin(parameters->theta_l)*sin(parameters->phi_l), 
			cos(parameters->theta_l)};
		parameters->incl_angle = acos(Neq[0]*Leq[0]+Neq[1]*Leq[1]+Neq[2]*Leq[2]);
		if(parameters->incl_angle == 0 || parameters->incl_angle == M_PI){
			std::cout<<"ERROR -- L == N || L == -N : The source frame is ill defined, and the conversion between theta_l and phi_l to theta_j and phi_j should not be trusted"<<std::endl; 
		}
		//Populate JSF here
		T JSF[3];
		if(generation_method.find("Pv2")!=std::string::npos){
			IMRPhenomPv2<T> model;
			model.PhenomPv2_JSF_from_params(parameters, JSF);
		}
		else{
			std::cout<<"ERROR -- Only Pv2 is supported right now"<<std::endl;	
		}
		T Jeq[3];
		equatorial_from_SF(JSF, parameters->theta_l, parameters->phi_l,(T) (theta_s), phi_s, parameters->incl_angle, parameters->phiRef, Jeq);
		T Jeqsph[3];
		transform_cart_sph(Jeq,Jeqsph);
		T theta_j = Jeqsph[1];
		T phi_j= Jeqsph[2];

		T iota_j;
		if(detector==""){
			terr_pol_iota_from_equat_sph(parameters->RA, parameters->DEC, theta_j, phi_j,&parameters->psi, &iota_j);
			ecl_from_eq(theta_j, phi_j, &parameters->theta_j_ecl, &parameters->phi_j_ecl);
		}
		else if(detector!="LISA"){
			terr_pol_iota_from_equat_sph(parameters->RA, parameters->DEC, theta_j, phi_j,&parameters->psi, &iota_j);
		}
		else{
			ecl_from_eq(theta_j, phi_j, &parameters->theta_j_ecl, &parameters->phi_j_ecl);
		}
	}
	else{
		T theta_s = M_PI/2. - parameters->DEC;
		T phi_s = parameters->RA;
		//N should point source to detector
		T Neq[3] = {-sin(theta_s)*cos(phi_s), -sin(theta_s)*sin(phi_s), -cos(theta_s)};
		T Leq[3] = {sin(parameters->theta_l)*cos(parameters->phi_l), 
			sin(parameters->theta_l)*sin(parameters->phi_l), 
			cos(parameters->theta_l)};
		parameters->incl_angle = acos(Neq[0]*Leq[0]+Neq[1]*Leq[1]+Neq[2]*Leq[2]);
		if(detector==""){
			terr_pol_iota_from_equat_sph(parameters->RA, parameters->DEC,parameters->theta_l, parameters->phi_l, &parameters->psi, &parameters->incl_angle);
			ecl_from_eq(parameters->theta_l, parameters->phi_l, &parameters->theta_j_ecl, &parameters->phi_j_ecl);
		}
		if(detector!="LISA"){
			terr_pol_iota_from_equat_sph(parameters->RA, parameters->DEC,parameters->theta_l, parameters->phi_l, &parameters->psi, &parameters->incl_angle);
		}
		else{
			ecl_from_eq(parameters->theta_l, parameters->phi_l, &parameters->theta_j_ecl, &parameters->phi_j_ecl);
		}
	
	}
}
template void transform_orientation_coords<double>(gen_params_base<double> *, std::string,std::string);
template void transform_orientation_coords<adouble>(gen_params_base<adouble> *, std::string,std::string);

/*! \brief For ADOL-C, assigns the frequencies at which logical breaks occur in the waveform
 *
 * ADOL-C cannot handle logically branches, so a tape must be made for each branch separately. This function returns the locations of these branches
 */
void assign_freq_boundaries(double *freq_boundaries, 
	double *intermediate_freqs, 
	int boundary_num, 
	gen_params_base<double> *input_params, 
	std::string generation_method)
{
	for(int i = 0 ; i<boundary_num; i++){
		freq_boundaries[i] = -1;
		intermediate_freqs[i] = -1;
	}
	if(input_params->equatorial_orientation){
		transform_orientation_coords(input_params,generation_method,"");
	}
	gen_params_base<adouble> internal_params;
	gen_params_base<adouble> *internal_params_ptr = &internal_params;
	transform_parameters(input_params, &internal_params_ptr);
	source_parameters<adouble> s_param;
	//s_param = source_parameters<adouble>::populate_source_parameters(&internal_params);
	s_param.populate_source_parameters(&internal_params);
	s_param.sky_average = internal_params.sky_average;
	s_param.f_ref = internal_params.f_ref;
	s_param.phiRef = internal_params.phiRef;
	s_param.cosmology=internal_params.cosmology;
	s_param.incl_angle=internal_params.incl_angle;
	s_param.shift_time = input_params->shift_time;
	if(generation_method.find("NRT") != std::string::npos){
		//Copied directly from prep_source_parameters 
		//Not ideal..
		if(input_params->tidal_love){
			/* The binary love relations are used here to compute lambda_a
	      		 * as a function of lambda_s (following equations 11-13 of
	      		 * arXiv:1903.03909). These relations were fit for Neutron stars,
	      		 * so the relevant coefficients/fit parameters are in
	      		 * IMRPhenomD_NRT.h.
	      		 */
	      		double q, Q, F;
	      		q = input_params->mass2 / input_params->mass1;
	      		Q = pow(q, 10./(3. - n_binLove));
	      		F = (1. - Q)/(1. + Q);

	      		double num = 1;
	      		double denom = 1;
	      		double q_pow[2], lambda_pow[3];
	      		q_pow[0] = q;
	      		q_pow[1] = q*q;
	      		lambda_pow[0] = pow(input_params->tidal_s, -1./5.);
	      		lambda_pow[1] = lambda_pow[0] * lambda_pow[0];
	      		lambda_pow[2] = lambda_pow[0] * lambda_pow[1];
	      		for(int i = 0; i<3; i++)
	      		  {
	      		    for(int j = 0; j<2; j++)
	      		      {
	      		        num += b_binLove[i][j]*q_pow[j]*lambda_pow[i];
	      		        denom += c_binLove[i][j]*q_pow[j]*lambda_pow[i];
	      		      }
	      		  }

	      		input_params->tidal_a = F * (num / denom) * input_params->tidal_s;


	      		s_param.tidal1 = input_params->tidal_s + input_params->tidal_a;
	      		s_param.tidal2 = input_params->tidal_s - input_params->tidal_a;
	      		
	      		//Gotta copy these over so that the next two lines run..
	      		input_params->tidal1= input_params->tidal_s + input_params->tidal_a;
	      		input_params->tidal2= input_params->tidal_s - input_params->tidal_a;

		}
		if((input_params->tidal1 < 0 || input_params->tidal2<0) && input_params->tidal_weighted >= 0) {
			s_param.tidal_weighted = input_params->tidal_weighted;
		}
		else if((input_params->tidal1 >= 0 && input_params->tidal2>=0) ) {
			s_param.tidal1 = input_params->tidal1;
			s_param.tidal2 = input_params->tidal2;
			//arXiv 1402.5156
			s_param.tidal_weighted = 8./13. * ( (1. + 7.*s_param.eta - 31.*s_param.eta*s_param.eta)*(s_param.tidal1 + s_param.tidal2) 
						+ sqrt( 1. - 4.*s_param.eta ) * ( 1. + 9.*s_param.eta - 11. * s_param.eta*s_param.eta) * (s_param.tidal1 - s_param.tidal2) ) ;
			s_param.delta_tidal_weighted = 1./2. * ( sqrt( 1. - 4.*s_param.eta ) * ( 1. - 13272./1319. * s_param.eta + 8944./1319. * s_param.eta*s_param.eta) *
						(s_param.tidal1 + s_param.tidal2) + ( 1. - 15910./1319. * s_param.eta + 32850./1319. * s_param.eta*s_param.eta + 3380./1319. 
						* s_param.eta *s_param.eta*s_param.eta)*(s_param.tidal1-s_param.tidal2));
		}
		
	}
	lambda_parameters<adouble> lambda;
	
	//IMRPhenomPv2<double> model;
	//model.PhenomPv2_inplane_spin(input_params);
	if(	(
		generation_method.find("IMRPhenomPv2") != std::string::npos
		)
		&& !input_params->sky_average){
		s_param.shift_time = false;
		double M,fRD,fpeak;
		IMRPhenomPv2<adouble> modelp;
		s_param.spin1z = internal_params.spin1[2];
		s_param.spin2z = internal_params.spin2[2];
		if((internal_params.chip + 1)>DOUBLE_COMP_THRESH){
			s_param.chip = internal_params.chip;
			s_param.phip = internal_params.phip;
			modelp.PhenomPv2_Param_Transform_reduced(&s_param);
		}
		else{
			s_param.spin1x = internal_params.spin1[0];
			s_param.spin1y = internal_params.spin1[1];
			s_param.spin2x = internal_params.spin2[0];
			s_param.spin2y = internal_params.spin2[1];
			modelp.PhenomPv2_Param_Transform(&s_param);
		}
		if(generation_method.find("gIMR") !=std::string::npos){
			gIMRPhenomPv2<adouble> gmodelp;
			gmodelp.assign_lambda_param(&s_param, &lambda);
			gmodelp.post_merger_variables(&s_param);
			M = s_param.M.value();
			fRD = s_param.fRD.value();
			fpeak = gmodelp.fpeak(&s_param, &lambda).value();
			std::cout<<fRD<<" "<<fpeak<<std::endl;
		
		}

		else{
			modelp.assign_lambda_param(&s_param, &lambda);
			modelp.post_merger_variables(&s_param);
			M = s_param.M.value();
			fRD = s_param.fRD.value();
			fpeak = modelp.fpeak(&s_param, &lambda).value();
		}
		//###########################################
		freq_boundaries[0] = .014/M;
		freq_boundaries[1] = .018/M;
		if(fRD/2. < fpeak){
			freq_boundaries[2] = fRD/2.;
			freq_boundaries[3] = fpeak;
		}
		else{
			freq_boundaries[3] = fRD/2.;
			freq_boundaries[2] = fpeak;
		}
		freq_boundaries[4] = .2/M;//End waveform
	}
	else if(
		generation_method.find("IMRPhenomD") != std::string::npos
		){

		double M, fRD, fpeak;
		if(generation_method.find("gIMR") !=std::string::npos){
			gIMRPhenomD<adouble> modeld;
			modeld.assign_lambda_param(&s_param, &lambda);
			modeld.post_merger_variables(&s_param);
			M = s_param.M.value();
			fRD = s_param.fRD.value();
			fpeak = modeld.fpeak(&s_param, &lambda).value();
			std::cout<<fRD<<" "<<fpeak<<std::endl;
		
		}
		else{
			IMRPhenomD<adouble> modeld;
			modeld.assign_lambda_param(&s_param, &lambda);
			modeld.post_merger_variables(&s_param);
			M = s_param.M.value();
			fRD = s_param.fRD.value();
			fpeak = modeld.fpeak(&s_param, &lambda).value();
		}
		//###########################################
		freq_boundaries[0] = .014/M;
		freq_boundaries[1] = .018/M;
		if(fRD/2. < fpeak){
			freq_boundaries[2] = fRD/2.;
			freq_boundaries[3] = fpeak;
		}
		else{
			freq_boundaries[3] = fRD/2.;
			freq_boundaries[2] = fpeak;
		}
		freq_boundaries[4] = .2/M;//End waveform
	}
	if(generation_method.find("NRT") != std::string::npos){
		adouble kappa_temp = (3./16.) * s_param.tidal_weighted;
		double kappa = kappa_temp.value();


		double a0 = 0.3586;
		double n1 = 3.35411203e-2;
		double n2 = 4.31460284e-5;
		double d1 = 7.54224145e-2;
		double d2 = 2.23626859e-4;
		adouble fmerger_temp = (1./(2.*s_param.M * M_PI))* a0* sqrt(s_param.mass2 / s_param.mass1) *(1.0 + n1 * kappa + n2 * kappa * kappa)/(1.0 + d1 * kappa + d2* kappa * kappa);
		double fmerger = fmerger_temp.value();
		//Sort the boundaries
		std::vector<double> temp_boundaries(boundary_num);
		double temp[2] = {fmerger, fmerger*1.2};
		for(int i =0 ;i<2; i++){
			temp_boundaries[i] = temp[i];
		}		
		for(int i =2 ;i<boundary_num; i++){
			temp_boundaries[i] = freq_boundaries[i-2];
		}		
		std::sort(temp_boundaries.begin(),temp_boundaries.end());
		for(int i =0 ;i<boundary_num; i++){
			freq_boundaries[i] = temp_boundaries[i];
		}
		
		//for (int i = 0 ; i<boundary_num; i++){
		//	std::cout<<i<<" "<<temp_boundaries[i]<<std::endl;
		//}
		
	}
	//###########################################
	intermediate_freqs[0] = freq_boundaries[0]*.9;
	for(int i = 1 ; i<boundary_num; i++){
		intermediate_freqs[i] = freq_boundaries[i-1]+(double)(freq_boundaries[i]-freq_boundaries[i-1])/2.;
	}
	if(check_mod(generation_method)){
		if(generation_method.find("ppE")!= std::string::npos ||
		check_theory_support(generation_method)){
			delete [] internal_params.betappe;
			delete [] internal_params.bppe;
		}
		else if(generation_method.find("gIMR")!= std::string::npos){
			if(internal_params.Nmod_phi !=0){
				delete [] internal_params.phii;
				delete [] internal_params.delta_phi;
			}
			if(internal_params.Nmod_sigma !=0){
				delete [] internal_params.sigmai;
				delete [] internal_params.delta_sigma;
			}
			if(internal_params.Nmod_beta !=0){
				delete [] internal_params.betai;
				delete [] internal_params.delta_beta;
			}
			if(internal_params.Nmod_alpha !=0){
				delete [] internal_params.alphai;
				delete [] internal_params.delta_alpha;
			}
		}
	}
	
	

}

/*! \brief Utility to find the integration bounds for Fisher matrices for increasing speed of Fisher evaluation
 *
 * Numerically finds the frequencies at which the Fisher should be evaluated at. 
 * 
 * Uses the bisection search algorithm for the cases where the waveform enters/leaves the band at SNR>1
 *
 * Uses a 100 pt grid search (logarithmically spaced) if the signal has SNR<1 when entering and leaving
 *
 * integrand_bounds[0] ~ frequency at which |h|/(sqrt S) ~signal_to_noise +/- tol
 *
 * integrand_bounds[1] ~ frequency at which |h|/(sqrt S) ~signal_to_noise +/- tol
 *
 * If autodiff set to true, ADOL-C routines are used. This is openmp safe, not necessarily thread safe.
 */
void integration_bounds(gen_params_base<double> *params, /**< Parameters of the waveform*/
	std::string generation_method, /*<< Generation method to use for the waveform*/
	std::string detector, /**< Detector to use for the response function*/
	std::string sensitivity_curve, /**< Sensitivity curve to use (must be one of the analytic curves in the detector_utilitiy file*/
	double fmin, /**< minimum frequency to use (specific to the detector)*/
	double fmax, /**< max frequency to use (specific to the detector)*/
	double signal_to_noise,/**< Target ratio of |h|/ sqrt(S) (typically ~0.1)*/
	double tol,/**< This is a numerical algorithm, so the tolerance must be specified*/
	double *integration_bounds,/**< [out] bounds fo the integral shape -- [2] -- (fmin,fmax)*/
	bool autodiff
	) 
{
	if(params->equatorial_orientation){
		transform_orientation_coords(params,generation_method,"");
	}
	double integration_time = 48;//Just for sensitivity curve
	bool lower=false, upper=false;
	std::complex<double> response;
	double eval_freq;
	double time;
	double psd;
	double ratio;
	
	//Check lowest frequency
	eval_freq = fmin;	
	if(autodiff){
		time_phase_corrected_autodiff(&time, 1, &eval_freq, params, 
			generation_method, false,(int *)NULL,1);
	}
	else{
		time_phase_corrected(&time, 1, &eval_freq, params, 
			generation_method, false,1);
	}
	fourier_detector_response_equatorial(&eval_freq, 1, &response, detector, 
		generation_method, params, &time);
	populate_noise(&eval_freq, sensitivity_curve, &psd, 1,integration_time);
	ratio = std::abs(response)/psd ;
	if( ratio>(signal_to_noise-tol))
	{
		integration_bounds[0] = eval_freq;
		lower = true;
	}

	//Check highest frequency
	eval_freq = fmax;
	if(autodiff){
		time_phase_corrected_autodiff(&time, 1, &eval_freq, params, 
			generation_method, false,(int *)NULL,1);
	}
	else{
		time_phase_corrected(&time, 1, &eval_freq, params, 
			generation_method, false,1);
	}
	fourier_detector_response_equatorial(&eval_freq, 1, &response, detector, 
		generation_method, params, &time);
	populate_noise(&eval_freq, sensitivity_curve, &psd, 1,integration_time);
	ratio = std::abs(response)/psd ;
	if( ratio>(signal_to_noise-tol))
	{
		integration_bounds[1] = eval_freq;
		upper = true;
	}

	//Search for bounds
	double fmin_search = fmin;	
	double fmax_search = fmax;	
	bool continue_search=true;
	if(upper && !lower)
	{
		while(continue_search)
		{
			//Bisection in log freq space
			eval_freq=sqrt(fmax_search*fmin_search);
			if(autodiff){
				time_phase_corrected_autodiff(&time, 1, &eval_freq, params, 
					generation_method, false,(int *)NULL,1);
			}
			else{
				time_phase_corrected(&time, 1, &eval_freq, params, 
					generation_method, false,1);
			}
			fourier_detector_response_equatorial(&eval_freq, 1, &response, detector, 
				generation_method, params, &time);
			populate_noise(&eval_freq, sensitivity_curve, &psd, 1,integration_time);
			ratio = std::abs(response)/psd ;
			
			//The function can be so steep the algorithm cannot determine a valid 
			//frequency with the required tolerance because of floating point error
			//check to see if the difference is at all meaningful
			if(2*(fmax_search - fmin_search)/(fmax_search+fmin_search) <1e-12){
				continue_search =false;
				integration_bounds[1]=eval_freq;

			}
			if(ratio > (signal_to_noise + tol)){
				fmax_search = eval_freq;
			}
			else if(ratio < (signal_to_noise - tol)){
				fmin_search = eval_freq;
			}
			else{
				continue_search =false;
				integration_bounds[0]=eval_freq;
			}
		}
	}
	else if(lower && !upper){
		while(continue_search)
		{
			//Bisection in log freq space
			eval_freq=sqrt(fmax_search*fmin_search);
			if(autodiff){
				time_phase_corrected_autodiff(&time, 1, &eval_freq, params, 
					generation_method, true,(int *)NULL,1);
			}
			else{
				time_phase_corrected(&time, 1, &eval_freq, params, 
					generation_method, false,1);
			}
			fourier_detector_response_equatorial(&eval_freq, 1, &response, detector, 
				generation_method, params, &time);
			populate_noise(&eval_freq, sensitivity_curve, &psd, 1,integration_time);
			ratio = std::abs(response)/psd ;

			//The function can be so steep the algorithm cannot determine a valid 
			//frequency with the required tolerance because of floating point error
			//check to see if the difference is at all meaningful
			if(2*(fmax_search - fmin_search)/(fmax_search+fmin_search) <1e-12){
				continue_search =false;
				integration_bounds[1]=eval_freq;

			}
			if(ratio > (signal_to_noise + tol)){
				fmin_search = eval_freq;
			}
			else if(ratio < (signal_to_noise - tol)){
				fmax_search = eval_freq;
			}
			else{
				continue_search =false;
				integration_bounds[1]=eval_freq;
			}
		}

	}
	//For now, just doing a sparse grid search. This can be refined later 
	//After doing grid search, it would be best to do a bisection search like above
	//with the new bounds after the grid search
	else if(!lower && !upper){
		//Grid search to see if any point is above the noise curve
		int vec_length= 100;
		std::complex<double> response_vec[vec_length];
		double time_vec[vec_length];
		double freq_vec[vec_length];
		double psd_vec[vec_length];
		double ratio_vec[vec_length];

		//Logarithmically space frequencies
		double delta_f_factor = pow(fmax/fmin,1./(vec_length-1));
		for(int i = 0 ; i<vec_length; i++){
			freq_vec[i]= pow_int(delta_f_factor,i)*fmin;
		}
		if(autodiff){
			time_phase_corrected_autodiff(time_vec, vec_length, freq_vec, params, 
				generation_method, false,(int *)NULL,1);
		}
		else{
			time_phase_corrected(time_vec, vec_length, freq_vec, params, 
				generation_method, false,1);
		}
		fourier_detector_response_equatorial(freq_vec, vec_length, response_vec, detector, 
			generation_method, params, time_vec);
		populate_noise(freq_vec, sensitivity_curve, psd_vec, vec_length,integration_time);
		for(int i = 0 ; i<vec_length; i++){
			ratio_vec[i]=std::abs(response_vec[i])/psd_vec[i];
		}
		bool search_lower = true, search_upper = true;
		for(int i = 0 ;i<vec_length; i++){
			if(ratio_vec[i] >(signal_to_noise -tol) )
			{
				fmin_search = freq_vec[i];
				search_lower=false;
				break;
			}		

		}
		for(int i = vec_length-1 ;i>0; i--){
			if(ratio_vec[i] >(signal_to_noise -tol) )
			{
				fmax_search = freq_vec[i];
				search_upper=false;
				break;
			}		

		}
		if(!search_lower && !search_upper){
			integration_bounds[0]=fmin_search;
			integration_bounds[1]=fmax_search;
		}
		else{
			integration_bounds[0]=fmax;
			integration_bounds[1]=fmax;
		}
	}
}
/*! \brief Determines the integration bounds for the log likelihood or fisher given some observation time, sampling frequency, detector, and sensitivity curve
 *
 * Sensitivity curve has to be one of the options in detector_util analytic options
 *
 * The current scheme is to use the frequency bounds determined by the SNR if the binary spends less than the integration time in band. If the merger spends more time in band than the integration time, the frequencies are determined to be (f_integration_time, f_high_band)
 *
 * Accounts for the integration time, unlike the integration_bounds routine
 *
 * returns 0 if successful, returns 1 if bounds could not be found due to roundoff, and returns 2 if entirely unsuccessful
 *
 * autodiff set to true will use ADOL-C for time-phase relationships. This is more accurate, but slower and only threadsafe with OpenMP
 */
int observation_bounds(double sampling_freq, /**< Frequency at which the detector operates*/
	double integration_time, /**< Time of observation in seconds*/
	std::string detector, /**< Detector to use for the response function*/
	std::string sensitivity_curve, /**< Sensitivity curve to use -- must match analytic choices in detector_util*/
	std::string generation_method,/**< method to use for the waveform generation*/
	gen_params_base<double> *params,/**< parameters of the source*/
	double *freq_bounds,/**< [out] Output bounds*/
	bool autodiff	
	)
{
	if(params->equatorial_orientation){
		transform_orientation_coords(params,generation_method,detector);
	}
	double fmax = sampling_freq /2.; //Nyquist limit
	double fmin= 1e-6; //DC component
	double delta_f =  1./integration_time;
	double SN_target= 0.01;
	
	double bounds_from_band[2];
	integration_bounds( params, generation_method, detector, sensitivity_curve, fmin, fmax, SN_target, .01, bounds_from_band,autodiff);
	double times[2];
	if(autodiff){
		time_phase_corrected_autodiff(times, 2, bounds_from_band, params, generation_method, false, (int *)NULL, 1);
	}
	else{
		time_phase_corrected(&times[0], 1, &bounds_from_band[0], params, generation_method, false,1);
		time_phase_corrected(&times[1], 1, &bounds_from_band[1], params, generation_method, false,1);
	}
	double T_band = -times[1]+times[0];

	//Conform the output of the band calculation to the pre-set frequency grid
	bool max_found=false, min_found=false;

	if(T_band < integration_time){
		freq_bounds[0]=bounds_from_band[0];	
		freq_bounds[1]=bounds_from_band[1];	
		return 0;
	}
	else{
		freq_bounds[1] = bounds_from_band[1];
		bool continue_search=true;
		double eval_freq, time;
		double max_search=bounds_from_band[1], min_search=bounds_from_band[0];
		double tolerance = .01*integration_time;
		
		while(continue_search)
		{
			//eval_id = (max_id_search +min_id_search)/2;
			eval_freq = (max_search +min_search)/2;
			if(autodiff){
				time_phase_corrected_autodiff(&time, 1, &eval_freq, params, 
					generation_method, false, (int *)NULL,1);
			}
			else{
				time_phase_corrected(&time, 1, &eval_freq, params, 
					generation_method, false,1);
			}
			if( 	( ( (-times[1] + time) > integration_time-tolerance )  && 
				( (-times[1] + time) < integration_time+tolerance) ) ) 
			{
				continue_search = false;
				freq_bounds[0]=eval_freq;
				return 0;
			}
			else if( fabs(max_search-min_search)<DOUBLE_COMP_THRESH){
				continue_search = false;
				freq_bounds[0] = eval_freq;
				return 1;

			}
			else if(( (-times[1] + time) < (integration_time-tolerance) )){
				max_search = eval_freq;
			}
			else if(( (-times[1] + time) > (integration_time+tolerance) )){
				min_search = eval_freq;
			}

		}
	}
	return 2;
	
}

/*! \brief ***DEPRECATED*** Determines the integration bounds for the log likelihood or fisher given some observation time, sampling frequency, detector, and sensitivity curve
 *
 * Sensitivity curve has to be one of the options in detector_util analytic options
 *
 * The current scheme is to use the frequency bounds determined by the SNR if the binary spends less than the integration time in band. If the merger spends more time in band than the integration time, the frequencies are determined to be (f_integration_time, f_high_band)
 */
void observation_bounds_discrete(double sampling_freq, /**< Frequency at which the detector operates*/
	double integration_time, /**< Time of observation in seconds*/
	std::string detector, /**< Detector to use for the response function*/
	std::string sensitivity_curve, /**< Sensitivity curve to use -- must match analytic choices in detector_util*/
	std::string generation_method,/**< method to use for the waveform generation*/
	gen_params_base<double> *params,/**< parameters of the source*/
	double *freq_bounds/**< [out] Output bounds*/
	)
{
	if(params->equatorial_orientation){
		transform_orientation_coords(params,generation_method,detector);
	}
	double fmax = sampling_freq /2.; //Nyquist limit
	//double fmin= 0; //DC component
	double fmin= 1e-6; //DC component
	double delta_f =  1./integration_time;
	int N = (fmax-fmin)/(delta_f)+1;
	double *frequencies = new double[N];
	for(int i = 0 ; i<N; i++){
		frequencies[i] = fmin + i*delta_f;	
	}
	
	bool autodiff=false;
	double bounds_from_band[2];
	integration_bounds( params, generation_method, detector, sensitivity_curve, fmin, fmax, .1, .01, bounds_from_band,autodiff);
	double times[2];
	//time_phase_corrected_autodiff(times, 2, bounds_from_band, params, generation_method, true);
	time_phase_corrected(times, 2, bounds_from_band, params, generation_method, true);
	double T_band = -times[1]+times[0];
	//std::cout<<T_band<<std::endl;

	//Conform the output of the band calculation to the pre-set frequency grid
	bool max_found=false, min_found=false;
	int max_id , min_id ;
	int i=0 ;

	while((!max_found || !min_found) && i<N){
		if(bounds_from_band[1] >= (frequencies[i]-delta_f/2. ) &&
			bounds_from_band[1] <= (frequencies[i]+delta_f/2. )){
			max_id = i;
			max_found = true;
		}
		if(bounds_from_band[0] >= (frequencies[i]-delta_f/2. ) &&
			bounds_from_band[0] <= (frequencies[i]+delta_f/2. )){
			min_id = i;
			min_found = true;
		}
		i++;
	}

	if(T_band < integration_time){
		freq_bounds[0]=frequencies[min_id];	
		freq_bounds[1]=frequencies[max_id];	
	}
	else{
		freq_bounds[1] = frequencies[max_id];
		bool continue_search=true;
		double eval_freq, time;
		int min_id_search=min_id, max_id_search=max_id, eval_id;
		double tolerance = .1*integration_time;
		
		while(continue_search)
		{
			eval_id = (max_id_search +min_id_search)/2;
			eval_freq = frequencies[eval_id];
			//time_phase_corrected_autodiff(&time, 1, &eval_freq, params, 
			//	generation_method, true);
			time_phase_corrected(&time, 1, &eval_freq, params, 
				generation_method, true);
			if( 	( ( (-times[1] + time) > integration_time-tolerance )  && 
				( (-times[1] + time) < integration_time+tolerance) )
				||
				(freq_bounds[1]-eval_freq < 100*delta_f) ) 
			{

				continue_search = false;
				freq_bounds[0]=eval_freq;
				//std::cout<<(times[1]-time)/T_year<<std::endl;
			}
			else if(( (-times[1] + time) < (integration_time-tolerance) )){
				max_id_search = eval_id;
			}
			else if(( (-times[1] + time) > (integration_time+tolerance) )){
				min_id_search = eval_id;
			}

		}
	}
	
	delete [] frequencies;
}
/*! \brief Calculates the postmerger parameters for a parameter set
 *
 * calculates fpeak, fdamp, and fRD
 */
template<class T>
void postmerger_params(gen_params_base<T>*params,
	std::string generation_method,
	T *fpeak,
	T *fdamp,
	T *fRD
	)
{
	if(generation_method.find("IMRPhenomD")!=std::string::npos){
		source_parameters<T> s_param;
		//s_param = source_parameters<T>::populate_source_parameters(params);
		s_param.populate_source_parameters(params);
		s_param.sky_average = params->sky_average;
		s_param.f_ref = params->f_ref;
		s_param.phiRef = params->phiRef;
		s_param.cosmology=params->cosmology;
		s_param.incl_angle=params->incl_angle;
		if(generation_method.find("gIMR")!=std::string::npos){
			gIMRPhenomD<T> model;
			lambda_parameters<T> lambda;
			model.assign_lambda_param(&s_param,&lambda);	
			model.post_merger_variables(&s_param);
			*fRD = s_param.fRD;
			*fdamp = s_param.fdamp;
			*fpeak = model.fpeak(&s_param , &lambda);
		}
		else{
			IMRPhenomD<T> model;
			lambda_parameters<T> lambda;
			model.assign_lambda_param(&s_param,&lambda);	
			model.post_merger_variables(&s_param);
			*fRD = s_param.fRD;
			*fdamp = s_param.fdamp;
			*fpeak = model.fpeak(&s_param , &lambda);

		}
	}
	else if(generation_method.find("IMRPhenomPv2")!=std::string::npos){
		source_parameters<T> s_param;
		//s_param = source_parameters<T>::populate_source_parameters(params);
		s_param.populate_source_parameters(params);
		s_param.sky_average = params->sky_average;
		s_param.f_ref = params->f_ref;
		s_param.phiRef = params->phiRef;
		s_param.cosmology=params->cosmology;
		s_param.incl_angle=params->incl_angle;
		s_param.chip = params->chip;
		s_param.phip = params->phip;
		if(generation_method.find("gIMR")!=std::string::npos){
			gIMRPhenomPv2<T> model;
			lambda_parameters<T> lambda;
			model.assign_lambda_param(&s_param,&lambda);	
			model.post_merger_variables(&s_param);
			*fRD = s_param.fRD;
			*fdamp = s_param.fdamp;
			*fpeak = model.fpeak(&s_param , &lambda);
		}
		else{
			IMRPhenomPv2<T> model;
			lambda_parameters<T> lambda;
			model.assign_lambda_param(&s_param,&lambda);	
			model.post_merger_variables(&s_param);
			*fRD = s_param.fRD;
			*fdamp = s_param.fdamp;
			*fpeak = model.fpeak(&s_param , &lambda);
		}
	}

}
template void postmerger_params<double>(gen_params_base<double> *,std::string, double *, double *, double*);
template void postmerger_params<adouble>(gen_params_base<adouble> *,std::string, adouble *, adouble *, adouble*);

struct Tbm_struct
{
	gen_params_base<double> *g_param;
	std::string method;
	double Tbm;
	bool ad;
};
/*! \brief Convenience function to Calculate the time before merger using numerical methods
 *
 * Tbm should be positive
 *
 * tol refers to the tolerance of the bisection search. If using numerical derivatives,
 * tolerance is not guaranteed because of errors introduced in the phase-time conversion
 *
 * Also, tolerance is associated with Tbm. Closer to merger, this can correspond to much larger errors on frequency because the relationship is so steep between time and frequency near merger.
 *
 * Uses numerical if autodiff set to false-- omp safe and thread safe
 *
 * Numerical accuracy is spotty. The second derivative isn't always reliable
 *
 * Relative time: 
 *
 * 	true -- Tbm is interpreted as the time before tc, and is shifted accordingly
 *
 * 	false -- Tbm is interpreted as some absolute time, and is not shifted. Assumed to be in the same time coordinates as tc. ie, if Tbm set to 0, freq will be the frequency at a time tc before tc
 *
 * Return values -- See util.cpp, newton_raphson_method_1d
 */
int Tbm_to_freq(gen_params_base<double> *params,/**< Generation parameters of the source*/
	std::string generation_method,/**< Generation method for the waveform*/
	double Tbm,/**< target time before merger -- in seconds or years (years if YEARS==true)*/
	double *freq,/**< Frequency at the input time before merger*/
	double tol, /**< Tolerance for the scheme*/
	bool autodiff, /**<Use autodiff routines instead of numerical derivatives*/
	int max_iterations,
	bool relative_time /**< True if shifting time by tc (time before merger) and false if not shifting time (absolute time frequency)*/
	)
{
	double scale_factor=1;
	if(params->equatorial_orientation){
		transform_orientation_coords(params,generation_method,"");
	}
	double fpeak,fRD,fdamp ;
	postmerger_params(params, generation_method, &fpeak, &fdamp, &fRD);
	
	Tbm_struct helper_struct;
	helper_struct.g_param = params;
	helper_struct.method = generation_method;
	//Shift target time or not
	if(relative_time){
		//helper_struct.Tbm = Tbm-params->tc;
		helper_struct.Tbm = Tbm+params->tc;
	}
	else{
		helper_struct.Tbm = Tbm;
	}
	helper_struct.ad = autodiff;
	//Zero PN assumes the merger time is zero
	double f0;
	f0 = f_0PN(helper_struct.Tbm,calculate_chirpmass(params->mass1,params->mass2)*MSOL_SEC);
	//std::cout<<"Initial Guess: "<<f0<<std::endl;
	//if(relative_time){
	//	f0 = f_0PN(Tbm+params->tc,calculate_chirpmass(params->mass1,params->mass2)*MSOL_SEC);
	//	//f0 = f_0PN(Tbm-params->tc,calculate_chirpmass(params->mass1,params->mass2)*MSOL_SEC);
	//}
	//else{
	//	//f0 = f_0PN(Tbm+params->tc,calculate_chirpmass(params->mass1,params->mass2)*MSOL_SEC);
	//	f0 = f_0PN(Tbm,calculate_chirpmass(params->mass1,params->mass2)*MSOL_SEC);
	//}
	int status = newton_raphson_method_1d(&Tbm_subroutine, f0, tol, max_iterations, (void *)(&helper_struct), freq);
	return status;
}

void Tbm_subroutine( double f, double *t, double *tp, void *param)
{
	Tbm_struct *local_param = (Tbm_struct *)param;
	gen_params_base<double> *gen_params = local_param->g_param;
	std::string method = local_param->method;
	if(local_param->ad){
		time_phase_corrected_autodiff(t, 1, &f, gen_params, method, false, (int *)NULL, 1);
		*t-=local_param->Tbm;
		time_phase_corrected_autodiff(tp, 1, &f, gen_params, method, false, (int *)NULL, 2);
	}
	else {
		time_phase_corrected(t, 1, &f, gen_params, method, false, 1);
		*t-=local_param->Tbm;
		//*t=local_param->Tbm-*t;
		time_phase_corrected(tp, 1, &f, gen_params, method, false,  2);
	}
	//std::cout<<f<<" "<<*t<<" "<<*tp<<std::endl;
}
/*! \brief *****OUTDATED***** Convenience function to Calculate the time before merger using numerical methods
 *
 * Tbm should be positive
 *
 * tol refers to the tolerance of the bisection search. If using numerical derivatives,
 * tolerance is not guaranteed because of errors introduced in the phase-time conversion
 *
 * Also, tolerance is associated with Tbm. Closer to merger, this can correspond to much larger errors on frequency because the relationship is so steep between time and frequency near merger.
 *
 * Uses numerical if autodiff set to false-- omp safe and thread safe
 */
void _Tbm_to_freq(gen_params_base<double> *params,/**< Generation parameters of the source*/
	std::string generation_method,/**< Generation method for the waveform*/
	double Tbm,/**< target time before merger -- in seconds or years (years if YEARS==true)*/
	double *freq,/**< Frequency at the input time before merger*/
	double tol, /**< Tolerance for the scheme*/
	bool autodiff, /**<Use autodiff routines instead of numerical derivatives*/
	bool YEARS/**<What units to use -- Years (true) or seconds(false)*/
	)
{
	double scale_factor=1;
	if(YEARS){scale_factor = T_year;}
	if(params->equatorial_orientation){
		transform_orientation_coords(params,generation_method,"");
	}
	double fpeak,fRD,fdamp ;
	postmerger_params(params, generation_method, &fpeak, &fdamp, &fRD);
	//std::cout<<"f fpeak "<<fpeak<<std::endl;	
	double time_peak=1, time_Tbm;
	if(autodiff){
		time_phase_corrected_autodiff(&time_peak, 1, &fpeak, params, 
			generation_method, false);
		time_peak/=scale_factor;
	}
	else{
		time_phase_corrected(&time_peak, 1, &fpeak, params, 
			generation_method, false);
		time_peak/=scale_factor;
	}
	//std::cout<<"Time fpeak "<<time_peak<<std::endl;	
	bool continue_search = true;
	double Tbmp = (1+tol)*Tbm;
	double Tbmm = (1-tol)*Tbm;
	double eval_freq, time;
	double fmax_search = fpeak;
	double fmin_search = fpeak;
	double T=0;
	while(T<Tbm){
		fmin_search = .1*fmin_search;
		if(autodiff){
			time_phase_corrected_autodiff(&time, 1, &fmin_search, params, 
				generation_method, false);
			time/=scale_factor;
		}
		else{
			time_phase_corrected(&time, 1, &fmin_search, params, 
				generation_method, false);
			time/=scale_factor;
		}
		//T = time_peak- time;
		T = -time_peak+ time;
		//std::cout<<"Time  "<<time<<std::endl;	
	
	}
	while(continue_search)
	{
		//Bisection in log freq space
		eval_freq=sqrt(fmax_search*fmin_search);
		if(autodiff){
			time_phase_corrected_autodiff(&time, 1, &eval_freq, params, 
				generation_method, true);
			time/=scale_factor;
		}
		else{
			time_phase_corrected(&time, 1, &eval_freq, params, 
				generation_method, false);
			time/=scale_factor;
		}
		//T = time_peak-time;	
		T = -time_peak+time;	
		//The function can be so steep the algorithm cannot determine a valid 
		//frequency with the required tolerance because of floating point error
		//check to see if the difference is at all meaningful
		if(2*(fmax_search - fmin_search)/(fabs(fmax_search)+fabs(fmin_search)) <DOUBLE_COMP_THRESH){
			std::cout<<"ROUNDING"<<std::endl;
			continue_search =false;
			
			*freq=eval_freq;

		}
		else if(T > (Tbmp )){
			fmin_search = eval_freq;
		}
		else if(T< (Tbmm)){
			fmax_search = eval_freq;
		}
		else{
			continue_search =false;
			*freq=eval_freq;
			//std::cout<<T/T_year<<std::endl;
		}
	}
}
/*! \brief ***DEPRECATED*** Utility for calculating the threshold times before merger that result in an SNR>SNR_thresh
 *
 * See arXiv 1902.00021
 *
 * Binary must merge within time T_wait
 *
 * SNR is calculated with frequencies [f(t_mer),f(t_mer-T_obs)] or [f(t_mer),0] depending on whether the binary has merged or not
 *
 * Assumes sky average -- Only supports PhenomD for now
 *
 * If no time before merger satisfies the requirements, both are set to -1
 */
void threshold_times(gen_params_base<double> *params,
	std::string generation_method, /**<Generation method to use for the waveform*/
	double T_obs, /**<Observation time -- also specifies the frequency spacing (\delta f = 1./T_obs)*/
	double T_wait, /**<Wait time -- Maximum time for binaries to coalesce */
	double f_lower,/**<Lower bound of search*/
	double f_upper,/**<upper bound of search*/
	std::string SN,/**< Noise curve name*/
	double SNR_thresh, /**< Threshold SNR */
	double *threshold_times_out,/**<[out] Output frequencies */
	double tolerance /**< Percent tolerance on SNR search*/
	)
{
	int length = (f_upper-f_lower)/T_obs;
	double deltaf = 1./T_obs;
	double *freqs = new double[length];
	for(int i = 0 ; i<length; i++){
		freqs[i]=f_lower +i*deltaf;
	}
	double *SN_curve = new double[length];
	populate_noise(freqs, SN, SN_curve, length);
	for(int i = 0 ; i<length; i++){
		SN_curve[i] = SN_curve[i]*SN_curve[i];	
	}
	threshold_times(params, generation_method, T_obs, T_wait, freqs, SN_curve, length, SNR_thresh, threshold_times_out,tolerance);
	delete [] freqs;
	delete [] SN_curve;

}

/*! \brief Utility for calculating the threshold times before merger that result in an SNR>SNR_thresh
 *
 * See arXiv 1902.00021
 *
 * Binary must merge within time T_wait
 *
 * SNR is calculated with frequencies [f(t_mer),f(t_mer-T_obs)] or [f(t_mer),0] depending on whether the binary has merged or not
 *
 * Assumes sky average -- Only supports PhenomD for now -- No angular dependence used ( only uses plus polarization -- assumes iota = psi = 0 )
 *
 * Assumes this is for multiband -- ie stellar mass BHs -- Only uses pn approximation of time frequency relation
 *
 * If no time before merger satisfies the requirements, both are set to -1
 *
 * Only supports tensor polarizations !!!
 */
void threshold_times(gen_params_base<double> *params,
	std::string generation_method, /**<Generation method to use for the waveform*/
	double T_obs, /**<Observation time -- also specifies the frequency spacing (\delta f = 1./T_obs)*/
	double T_wait, /**<Wait time -- Maximum time for binaries to coalesce */
	double *freqs,/**<Maximum frequency array*/
	double *SN,/**< Noise curve array, should be prepopulated from f_lower to f_upper with spacing 1./T_obs*/
	int length,/**< Length of maximum frequency array*/
	double SNR_thresh, /**< Threshold SNR */
	double *threshold_times_out,/**<[out] Output frequencies */
	double tolerance /**< Percent tolerance on SNR search*/
	)
{
	if(!params->sky_average){ std::cout<<"NOT sky averaged -- This is not supported by threshold_freqs"<<std::endl;}
	
	params->sky_average = false;
	double bounds[2];
	
	//Max number of iterations -- safety net
	int max_iter = 100;
	int ct = 0;
	double deltaf = freqs[1]-freqs[0];
	double chirpmass = calculate_chirpmass(params->mass1, params->mass2)*MSOL_SEC;
	std::complex<double> *hplus= new std::complex<double>[length];
	std::complex<double> *hcross= new std::complex<double>[length];
	bool not_found = true;
	int bound_id_lower = 0, bound_id_upper = length-1;//Current ids
	int bound_id_lower_prev = 0, bound_id_upper_prev = length-1;//Current ids
	double t_mer=0; //Time before merger
	double snr, snr_prev;
	int precalc_wf_id;
	//Determine if any t_mer between [0,T_wait] allows for an SNR>SNR_thresh 
	
	//Frequency T_obs before it leaves band
	//Using PN f(t) instead of local, because its faster and simpler -- maybe upgrade later
	//Shouldn't matter this far from merger
	t_mer= t_0PN(freqs[bound_id_upper], chirpmass)+T_obs;
	bound_id_lower = (f_0PN(t_mer, chirpmass)- freqs[0])/deltaf;

	//fourier_waveform(&freqs[bound_id_lower], bound_id_upper-bound_id_lower, &hplus[bound_id_lower], &hcross[bound_id_lower], generation_method,params);
	waveform_polarizations<double> wp;
	wp.hplus = &hplus[bound_id_lower];
	wp.hcross = &hcross[bound_id_lower];
	fourier_waveform(&freqs[bound_id_lower], bound_id_upper-bound_id_lower, &wp, generation_method,params);
	precalc_wf_id = bound_id_lower;
	snr = std::sqrt(1.)*calculate_snr_internal(&SN[bound_id_lower], &hplus[bound_id_lower],&freqs[bound_id_lower],bound_id_upper-bound_id_lower);
	snr_prev=snr;
	bound_id_lower_prev = bound_id_lower;
	bound_id_upper_prev = bound_id_lower;
	double t1 = t_mer, t2=t_mer;
	if(snr>SNR_thresh){not_found= false;}
	else{
		bool bound_search=true, t1_moved=true,t2_moved=true;
		while(bound_search){
			t_mer *=2.;
			bound_id_lower = (f_0PN(t_mer, chirpmass)- freqs[0])/deltaf;
			bound_id_upper = (f_0PN(t_mer-T_obs, chirpmass)- freqs[0])/deltaf;

			//fourier_waveform(&freqs[bound_id_lower], bound_id_lower_prev-bound_id_lower, &hplus[bound_id_lower], &hcross[bound_id_lower], generation_method,params);
			wp.hplus = &hplus[bound_id_lower];
			wp.hcross = &hcross[bound_id_lower];
			fourier_waveform(&freqs[bound_id_lower], bound_id_lower_prev-bound_id_lower, &wp, generation_method,params);
			precalc_wf_id = bound_id_lower;
			snr = std::sqrt(1.)*calculate_snr_internal(&SN[bound_id_lower], &hplus[bound_id_lower],&freqs[bound_id_lower],bound_id_upper-bound_id_lower);
			if(snr<snr_prev){bound_search=false;}	
		}	
		snr_prev=snr;
		t2 = t_mer;
		while(not_found &&ct < max_iter && t_mer < T_wait && (t2-t1)>T_day){
			t_mer = (t1+t2)/2.;
			bound_id_lower = ( f_0PN(t_mer ,chirpmass) - freqs[0])/deltaf;
			if(t_mer-T_obs > 0){
				bound_id_upper = ( std::min( f_0PN(  t_mer - T_obs,chirpmass) , freqs[length-1] ) - freqs[0])/deltaf;
			}
			else{
				bound_id_upper = (  freqs[length-1]  - freqs[0])/deltaf;

			}
			//Update waveform if using new frequencies
			if(bound_id_lower < precalc_wf_id){
				//fourier_waveform(&freqs[bound_id_lower], precalc_wf_id-bound_id_lower, &hplus[bound_id_lower], &hcross[bound_id_lower], generation_method,params);

				wp.hplus = &hplus[bound_id_lower];
				wp.hcross = &hcross[bound_id_lower];
				fourier_waveform(&freqs[bound_id_lower], precalc_wf_id-bound_id_lower, &wp, generation_method,params);
				precalc_wf_id = bound_id_lower;
			}
			bound_id_lower_prev= bound_id_lower;
			bound_id_upper_prev= bound_id_upper;
			snr =std::sqrt(1.)*calculate_snr_internal(&SN[bound_id_lower], &hplus[bound_id_lower],&freqs[bound_id_lower],bound_id_upper-bound_id_lower);
			ct++;
			if(snr>SNR_thresh){not_found=false;}
			else{
				if(t1_moved){
					if(snr>snr_prev){t1=t_mer;t1_moved=true; t2_moved=false;}
					else{t2=t_mer;t2_moved=true; t1_moved=false;;}
				}
				if(t2_moved){
					if(snr>snr_prev){t2=t_mer;t2_moved=true; t1_moved=false;}
					else{t1=t_mer;t1_moved=true; t2_moved=false;;}

				}
			}
			snr_prev = snr;

		}
	}	
	
	//If no SNR is larger than threshold, return 
	if(not_found){
		threshold_times_out[0] = -1;
		threshold_times_out[1] = -1;
	}
	//Find roots
	else{
		double t_save = t_mer;
		ct=0;
		bool found_lower_root=false;	
		bool found_upper_root=false;	
		t1=t_save, t2=t_save;//We know t_mer is over the threshold
		//Find a lower bound for bisection search
		while(snr>SNR_thresh){
			t1/=2.;
			bound_id_lower = ( f_0PN(t1 ,chirpmass) - freqs[0])/deltaf;
			if(t1-T_obs > 0){
				bound_id_upper = ( std::min( f_0PN(  t1 - T_obs,chirpmass) , freqs[length-1] ) - freqs[0])/deltaf;
			}
			else{
				bound_id_upper = (  freqs[length-1]  - freqs[0])/deltaf;

			}
			snr =std::sqrt(1.)*calculate_snr_internal(&SN[bound_id_lower], &hplus[bound_id_lower],&freqs[bound_id_lower],bound_id_upper-bound_id_lower);
		}
		while(!found_lower_root){
			
			t_mer = (t1+t2)/2.;
			//Find new frequency bound ids
			bound_id_lower = ( f_0PN(t_mer ,chirpmass) - freqs[0])/deltaf;
			if(t_mer-T_obs > 0){
				bound_id_upper = ( std::min( f_0PN(  t_mer - T_obs,chirpmass) , freqs[length-1] ) - freqs[0])/deltaf;
			}
			else{
				bound_id_upper = (  freqs[length-1]  - freqs[0])/deltaf;

			}
			//Update waveform if using new frequencies
			if(bound_id_lower < precalc_wf_id){
				//fourier_waveform(&freqs[bound_id_lower], precalc_wf_id-bound_id_lower, &hplus[bound_id_lower], &hcross[bound_id_lower], generation_method,params);
				wp.hplus = &hplus[bound_id_lower];
				wp.hcross = &hcross[bound_id_lower];
				fourier_waveform(&freqs[bound_id_lower], precalc_wf_id-bound_id_lower, &wp, generation_method,params);
				precalc_wf_id = bound_id_lower;
			}
			bound_id_lower_prev= bound_id_lower;
			bound_id_upper_prev= bound_id_upper;
			snr =std::sqrt(1.)*calculate_snr_internal(&SN[bound_id_lower], &hplus[bound_id_lower],&freqs[bound_id_lower],bound_id_upper-bound_id_lower);
			ct++;
			if(std::abs(snr-SNR_thresh)/SNR_thresh<tolerance ){found_lower_root=true;threshold_times_out[0]=t_mer;}
			else{
				if(snr>SNR_thresh){ t2 = t_mer;	}
				else{ t1=t_mer;}
			}
			snr_prev=snr;
			
		}
		ct=0;
		t1=t_save; t2=t_save;
		do{
			t2*=2.;
			if(t2>T_wait){t2=T_wait;}	
			bound_id_lower = ( f_0PN(t2 ,chirpmass) - freqs[0])/deltaf;
			if(t2-T_obs > 0){
				bound_id_upper = ( std::min( f_0PN(  t2 - T_obs,chirpmass) , freqs[length-1] ) - freqs[0])/deltaf;
			}
			else{
				bound_id_upper = (  freqs[length-1]  - freqs[0])/deltaf;

			}
			if(bound_id_lower < precalc_wf_id){
				//fourier_waveform(&freqs[bound_id_lower], precalc_wf_id-bound_id_lower, &hplus[bound_id_lower], &hcross[bound_id_lower], generation_method,params);
				wp.hplus = &hplus[bound_id_lower];
				wp.hcross = &hcross[bound_id_lower];
				fourier_waveform(&freqs[bound_id_lower], precalc_wf_id-bound_id_lower, &wp, generation_method,params);
				precalc_wf_id = bound_id_lower;
			}
			snr =std::sqrt(1.)*calculate_snr_internal(&SN[bound_id_lower], &hplus[bound_id_lower],&freqs[bound_id_lower],bound_id_upper-bound_id_lower);
			if(t2==T_wait && snr>SNR_thresh){ found_upper_root=true; threshold_times_out[1]=T_wait;break;}
		}while(snr>SNR_thresh  );
		while(!found_upper_root){
			
			t_mer = (t1+t2)/2.;
			//Find new frequency bound ids
			bound_id_lower = ( f_0PN(t_mer ,chirpmass) - freqs[0])/deltaf;
			if(t_mer-T_obs > 0){
				bound_id_upper = ( std::min( f_0PN(  t_mer - T_obs,chirpmass) , freqs[length-1] ) - freqs[0])/deltaf;
			}
			else{
				bound_id_upper = (  freqs[length-1]  - freqs[0])/deltaf;

			}
			//Update waveform if using new frequencies
			if(bound_id_lower < precalc_wf_id){
				//fourier_waveform(&freqs[bound_id_lower], precalc_wf_id-bound_id_lower, &hplus[bound_id_lower], &hcross[bound_id_lower], generation_method,params);
				waveform_polarizations<double> wp;
				wp.hplus = &hplus[bound_id_lower];
				wp.hcross = &hcross[bound_id_lower];
				fourier_waveform(&freqs[bound_id_lower], precalc_wf_id-bound_id_lower, &wp, generation_method,params);
				precalc_wf_id = bound_id_lower;
			}
			bound_id_lower_prev= bound_id_lower;
			bound_id_upper_prev= bound_id_upper;
			snr =std::sqrt(1.)*calculate_snr_internal(&SN[bound_id_lower], &hplus[bound_id_lower],&freqs[bound_id_lower],bound_id_upper-bound_id_lower);
			ct++;
			if(std::abs(snr-SNR_thresh)/SNR_thresh<tolerance ){found_upper_root=true;threshold_times_out[1]=t_mer;}
			else{
				if(snr>SNR_thresh){ t1 = t_mer;	}
				else{ t2=t_mer;}
			}
			snr_prev=snr;
			
		}
	}
	//Cleanup
	delete [] hplus;
	delete [] hcross;
}

/*! \brief **NOT FINISHED -- DO NOT USE** Utility for calculating the threshold times before merger that result in an SNR>SNR_thresh --GSL quad integration implementation
 *
 * NOT FINISHED -- DO NOT USE 
 *
 * See arXiv 1902.00021
 *
 * Binary must merge within time T_wait
 *
 * SNR is calculated with frequencies [f(t_mer),f(t_mer-T_obs)] or [f(t_mer),0] depending on whether the binary has merged or not
 *
 * Assumes sky average -- Only supports PhenomD for now -- No angular dependence used ( only uses plus polarization -- assumes iota = psi = 0 )
 *
 * If no time before merger satisfies the requirements, both are set to -1
 *
 * ALL temporal quantities in seconds or Hz
 *
 * Return values: 
 * 	
 * 	0 -- success
 * 	
 * 	11-- Failure: SNR was 0 in lower bound
 *
 * 	12-- Failure: SNR was 0 in upper bound
 *
 * 	13 -- partial success: Closest values output, but roundoff error prevented the routine from reaching the desired accuracy
 * 	
 * 	14 -- partial success: numerical inversion for time -> frequency had errors 
 *
 * 	15 -- partial success: numerical inversion for time -> frequency had failed, and a PN approximate had to be used
 *
 * 	16 -- max iterations reached
 */
int threshold_times_full_gsl(gen_params_base<double> *params,
	std::string generation_method, /**<Generation method to use for the waveform*/
	double T_obs, /**<Observation time -- also specifies the frequency spacing (\delta f = 1./T_obs)*/
	double T_wait, /**<Wait time -- Maximum time for binaries to coalesce */
	double fmin,/**<Maximum frequency array*/
	double fmax,/**<Maximum frequency array*/
	std::string SN,/**< Noise curve array, should be prepopulated from f_lower to f_upper with spacing 1./T_obs*/
	double SNR_thresh, /**< Threshold SNR */
	double *threshold_times_out,/**<[out] Output times */
	double tolerance, /**< Percent tolerance on SNR search*/
	gsl_integration_workspace *w,
	int np
	)
{
	int status1=0, status2=0;
	bool forced_PN_approx = false;
	bool bad_time_freq_inversion = false;
	bool round_off_error=false;	
	bool max_iterations_reached=false;	

	threshold_times_out[0]=-1;
	threshold_times_out[1]=-1;
	if(params->equatorial_orientation){
		transform_orientation_coords(params,generation_method,"");
	}
	bool save_SA = params->sky_average;
	if(!params->sky_average){ std::cout<<"NOT sky averaged -- This is not supported by threshold_freqs"<<std::endl;}
	params->sky_average = false;
	//params->sky_average = true;
	
	bool relative_time = true;
	bool autodiff = true;	
	//bool gsl_integration=false;
	bool gsl_integration=true;

	int GL_length = 500;
	double *GL_freqs;
	double *GL_w;
	if(!gsl_integration){
		GL_freqs = new double[GL_length];
		GL_w = new double[GL_length];
		params->sky_average=save_SA;
	}


	bool stellar_mass = true;	
	int max_iter = 20;//Tbm
	int max_iter_search = 200;//Box search
	int search_ct;
	double fpeak, fdamp,fRD;
	postmerger_params(params,generation_method,&fpeak,&fdamp, &fRD);
	
	if(fpeak < fmax){
		stellar_mass=false;
	}
	//stellar_mass=true;

	
	//Max number of iterations -- safety net
	double chirpmass = calculate_chirpmass(params->mass1, params->mass2)*MSOL_SEC;
	bool not_found = true;
	double f_lower = fmin, f_upper = fmax;//Current ids
	double f_lower_prev = fmin, f_upper_prev = fmax;//Current ids
	double t_mer=T_obs; //Time before merger
	double snr, snr_prev;
	double rel_err = tolerance;
	//Determine if any t_mer between [0,T_wait] allows for an SNR>SNR_thresh 
	
	//Frequency T_obs before it leaves band
	//Using PN f(t) instead of local, because its faster and simpler -- maybe upgrade later
	//Shouldn't matter this far from merger
	if(stellar_mass){
		t_mer= t_0PN(f_upper, chirpmass)+T_obs;
		f_lower = f_0PN(t_mer,chirpmass);	
	}
	else{
		f_upper = fpeak*2.;	
		f_lower = f_0PN(t_mer,chirpmass);	
		f_lower = (f_lower>fmin) ? f_lower : fmin;
		f_upper = (f_upper<fmax) ? f_upper : fmax;
	}
	

	//#######################################################
	//#######################################################
	const gsl_min_fminimizer_type *T;
	gsl_min_fminimizer *s;
	double m = T_obs;
	double a = 0.01*T_year, b=T_wait;
	gsl_snr_struct helper_params;
	helper_params.params = params;
	helper_params.T_obs = T_obs;
	helper_params.T_wait = T_wait;
	helper_params.fmax = fmax;
	helper_params.fmin = fmin;
	helper_params.w = w;
	helper_params.np = np;
	helper_params.relerr = rel_err;
	helper_params.SN =SN ;
	helper_params.generation_method =generation_method ;
	gsl_function F;
	
	//std::cout<<params->mass1<<" "<<params->mass2<<" "<<params->Luminosity_Distance<<std::endl;	
	F.function = [](double T, void *param){
		//std::cout<<"TIME: "<<T/T_year<<std::endl;
		gsl_snr_struct *param_unpacked = (gsl_snr_struct *) param;
		gen_params_base<double> *params = param_unpacked->params;
		double chirpmass = 
			calculate_chirpmass(params->mass1,params->mass2)*MSOL_SEC;

		double f_min =  f_0PN(T,chirpmass) ;
		double f_max;
		if(T-param_unpacked->T_obs > 0){
			f_max =  std::min( f_0PN(  T - param_unpacked->T_obs,chirpmass) , param_unpacked->fmax ) ;
		}
		else{
			f_max =   param_unpacked->fmax;

		}
		//std::cout<<f_min<<" "<<f_max<<" "<<(T-param_unpacked->T_obs)/T_year<<std::endl;
		double snr =snr_threshold_subroutine(f_min,f_max,param_unpacked->relerr,params, param_unpacked->generation_method,param_unpacked->SN, param_unpacked->w, param_unpacked->np); 
		//std::cout<<"SNR: "<<snr<<std::endl;
		return -1.*snr;
	};
	F.params = (void *)&helper_params;

	//double test = GSL_FN_EVAL(&F,4*T_year);

	T=gsl_min_fminimizer_brent;
	s = gsl_min_fminimizer_alloc(T);
	gsl_min_fminimizer_set(s,&F, m, a, b);


	int iter = 0, max_itergsl=20,status;
	gsl_set_error_handler_off();
	do{
		iter++;
		status = gsl_min_fminimizer_iterate(s);
		if(status == GSL_EINVAL)
			std::cout<<"NO MINIMUM"<<std::endl;
		m = gsl_min_fminimizer_x_minimum(s);
		a = gsl_min_fminimizer_x_lower(s);
		b = gsl_min_fminimizer_x_upper(s);
		status = gsl_min_test_interval(a,b,0.001,0.0);
		if(status == GSL_SUCCESS)
			std::cout<<"SUCCESS"<<std::endl;
		//std::cout<<"Parameters: "<<m<<" "<<a<<" "<<b<<std::endl;
	}while(status == GSL_CONTINUE && iter < max_itergsl);

	threshold_times_out[0]=-1;
	threshold_times_out[1]=-1;
	//std::cout<<"done "<<m/T_year<<std::endl;
	gsl_min_fminimizer_free(s);
	//#######################################################
	//#######################################################
	

	return 0;
}

/*! \brief Utility for calculating the threshold times before merger that result in an SNR>SNR_thresh --GSL quad integration implementation
 *
 * See arXiv 1902.00021
 *
 * Binary must merge within time T_wait
 *
 * SNR is calculated with frequencies [f(t_mer),f(t_mer-T_obs)] or [f(t_mer),0] depending on whether the binary has merged or not
 *
 * Assumes sky average -- Only supports PhenomD for now -- No angular dependence used ( only uses plus polarization -- assumes iota = psi = 0 )
 *
 * If no time before merger satisfies the requirements, both are set to -1
 *
 * ALL temporal quantities in seconds or Hz
 *
 * Return values: 
 * 	
 * 	0 -- success
 * 	
 * 	11-- Failure: SNR was 0 in lower bound
 *
 * 	12-- Failure: SNR was 0 in upper bound
 *
 * 	13 -- partial success: Closest values output, but roundoff error prevented the routine from reaching the desired accuracy
 * 	
 * 	14 -- partial success: numerical inversion for time -> frequency had errors 
 *
 * 	15 -- partial success: numerical inversion for time -> frequency had failed, and a PN approximate had to be used
 *
 * 	16 -- max iterations reached
 */
int threshold_times_gsl(gen_params_base<double> *params,
	std::string generation_method, /**<Generation method to use for the waveform*/
	double T_obs, /**<Observation time -- also specifies the frequency spacing (\delta f = 1./T_obs)*/
	double T_wait, /**<Wait time -- Maximum time for binaries to coalesce */
	double fmin,/**<Maximum frequency array*/
	double fmax,/**<Maximum frequency array*/
	std::string SN,/**< Noise curve array, should be prepopulated from f_lower to f_upper with spacing 1./T_obs*/
	double SNR_thresh, /**< Threshold SNR */
	double *threshold_times_out,/**<[out] Output times */
	double *T_obs_SNR,
	double tolerance, /**< Percent tolerance on SNR search*/
	gsl_integration_workspace *w,
	int np
	)
{
	int status1=0, status2=0;
	bool forced_PN_approx = false;
	bool bad_time_freq_inversion = false;
	bool round_off_error=false;	
	bool max_iterations_reached=false;	

	threshold_times_out[0]=-1;
	threshold_times_out[1]=-1;
	if(params->equatorial_orientation){
		transform_orientation_coords(params,generation_method,"");
	}
	bool save_SA = params->sky_average;
	if(!params->sky_average){ std::cout<<"NOT sky averaged -- This is not supported by threshold_freqs"<<std::endl;}
	params->sky_average = false;
	//params->sky_average = true;
	
	bool relative_time = true;
	bool autodiff = true;	
	//bool gsl_integration=false;
	bool gsl_integration=true;

	int GL_length = 500;
	double *GL_freqs;
	double *GL_w;
	if(!gsl_integration){
		GL_freqs = new double[GL_length];
		GL_w = new double[GL_length];
		params->sky_average=save_SA;
	}


	bool stellar_mass = true;	
	int max_iter = 20;//Tbm
	int max_iter_search = 200;//Box search
	int search_ct;
	double fpeak, fdamp,fRD;
	postmerger_params(params,generation_method,&fpeak,&fdamp, &fRD);
	
	if(fpeak < fmax){
		stellar_mass=false;
	}
	//stellar_mass=true;
	//if(stellar_mass){
	//	std::cout<<"Stellar Masss"<<std::endl;
	//	
	//}
	//std::cout<<"Fpeak "<<fpeak<<std::endl;

	
	//Max number of iterations -- safety net
	double chirpmass = calculate_chirpmass(params->mass1, params->mass2)*MSOL_SEC;
	bool not_found = true;
	double f_lower = fmin, f_upper = fmax;//Current ids
	double f_lower_prev = fmin, f_upper_prev = fmax;//Current ids
	double t_mer=T_obs; //Time before merger
	double snr, snr_prev;
	double rel_err = tolerance;
	//Determine if any t_mer between [0,T_wait] allows for an SNR>SNR_thresh 
	
	//Frequency T_obs before it leaves band
	//Using PN f(t) instead of local, because its faster and simpler -- maybe upgrade later
	//Shouldn't matter this far from merger
	if(stellar_mass){
		t_mer= t_0PN(f_upper, chirpmass)+T_obs;
		f_lower = f_0PN(t_mer,chirpmass);	
	}
	else{
		f_upper = fpeak*2.;	
		f_lower = f_0PN(t_mer,chirpmass);	
		f_lower = (f_lower>fmin) ? f_lower : fmin;
		f_upper = (f_upper<fmax) ? f_upper : fmax;
	}

	if(gsl_integration){
		snr = snr_threshold_subroutine(	f_lower, f_upper, rel_err,params, generation_method,SN, w,np);
	}
	else{
		gauleg(log10(f_lower), log10(f_upper), GL_freqs, GL_w, GL_length);
		for(int k =0; k<GL_length; k++){
			GL_freqs[k] = pow(10.,GL_freqs[k]);
		}
		snr = calculate_snr(SN, "LISA",generation_method, params,GL_freqs,GL_length,"GAUSSLEG",GL_w,true);
	}
	//std::cout<<"Initial SNR calc done"<<std::endl;
	//std::cout<<"Initial SNR: "<<snr<<std::endl;
	*T_obs_SNR= snr;
	snr_prev=snr;
	double t1 = t_mer, t2=t_mer;
	if(snr>SNR_thresh){not_found= false;}
	else{
		bool bound_search=true, t1_moved=true,t2_moved=true;

		//Push t_mer back till an SNR has peaked
		//First SNR is calculated T_obs after the edge of the band, 
		//which should be close to max,
		//but the SNR should increase and peak after this, 
		//so this is a one direction search
		while(bound_search){
			t_mer *=2.;
			if(stellar_mass){
				f_lower = (f_0PN(t_mer, chirpmass));
				f_upper = (f_0PN(t_mer-T_obs, chirpmass));
			}
			else{
				//Much much slower
				status1 =Tbm_to_freq(params, generation_method, t_mer, &f_lower,tolerance,autodiff,max_iter,relative_time);
				if(status1 != 0){
					forced_PN_approx = true;
					f_lower = (f_0PN(t_mer, chirpmass));

				}
				if(t_mer - T_obs > 0){
					status2= Tbm_to_freq(params, generation_method, t_mer-T_obs, &f_upper,tolerance,autodiff,max_iter,relative_time);
					if(status2 != 0){
						forced_PN_approx = true;
						f_upper = (f_0PN(t_mer-T_obs, chirpmass));

					}
					f_upper = std::min(f_upper,fmax);
				}
				else{
					f_upper = 2.*fpeak;	
				}
				if(f_lower<0)
				{
					forced_PN_approx = true;
					f_lower = (f_0PN(t_mer, chirpmass));
				}
				if(f_upper<0)
				{
					forced_PN_approx = true;
					f_upper = (f_0PN(t_mer-T_obs, chirpmass));
				}
				if(f_lower>f_upper){
					forced_PN_approx = true;
					f_lower = (f_0PN(t_mer, chirpmass));
					f_upper =  std::min( f_0PN(  t_mer - T_obs,chirpmass) , fmax );
				}
				if(status1 !=0 || status2 != 0){
					bad_time_freq_inversion = true;
				}
				status1=0;
				status2=0;
			}
			
			if(gsl_integration){
				snr = snr_threshold_subroutine(	f_lower, f_upper, rel_err,params, generation_method,SN, w,np);
			}
			else{
				gauleg(log10(f_lower), log10(f_upper), GL_freqs, GL_w, GL_length);
				for(int k =0; k<GL_length; k++){
					GL_freqs[k] = pow(10.,GL_freqs[k]);
				}
				snr = calculate_snr(SN, "LISA",generation_method, params,GL_freqs,GL_length,"GAUSSLEG",GL_w,true);
			}
			if(snr<snr_prev){bound_search=false;}	
			else if(f_lower < fmin){ 
				t_mer = t_0PN(fmin, chirpmass); 
				bound_search=false;
			}
		}	
		//Now the range is bounded, search more systematically for any times 
		//that can produce an SNR greater than the threshold
		
		snr_prev=snr;
		if(t_mer > T_wait){
			t2 = T_wait;	
		}
		else{
			t2 = t_mer;
		}
		//This loop should run until an SNR > THRESH is found, 
		//the times are pushed back too far, 
		//or the range is down to 1 day (nothing was found)
		while(not_found && t_mer < T_wait && (t2-t1)>T_day){
			t_mer = (t1+t2)/2.;
			if(stellar_mass){
				f_lower =  f_0PN(t_mer ,chirpmass) ;
				if(t_mer-T_obs > 0){
					f_upper =  std::min( f_0PN(  t_mer - T_obs,chirpmass) , fmax );
				}
				else{
					f_upper =   fmax;

				}
			}
			else{
				status1 = Tbm_to_freq(params, generation_method, t_mer, &f_lower,tolerance,autodiff,max_iter,relative_time);
				if(status1 != 0){
					forced_PN_approx = true;
					f_lower = (f_0PN(t_mer, chirpmass));

				}
				if(t_mer-T_obs > 0){
					status2=Tbm_to_freq(params, generation_method, t_mer-T_obs, &f_upper,tolerance,autodiff,max_iter,relative_time);
					if(status2 != 0){
						forced_PN_approx = true;
						f_upper = (f_0PN(t_mer-T_obs, chirpmass));

					}
					f_upper = std::min(f_upper,fmax);
				}
				else{
					f_upper =   fpeak*2;

				}
				if(f_lower<0)
				{
					forced_PN_approx = true;
					f_lower = (f_0PN(t_mer, chirpmass));
				}
				if(f_upper<0){
					forced_PN_approx = true;
					f_upper =  std::min( f_0PN(  t_mer - T_obs,chirpmass) , fmax );
				}
				if(f_lower>f_upper){
					forced_PN_approx = true;
					f_lower = (f_0PN(t_mer, chirpmass));
					f_upper =  std::min( f_0PN(  t_mer - T_obs,chirpmass) , fmax );
				}
				if(status1 !=0 || status2 != 0){
					bad_time_freq_inversion = true;
				}
				status1=0;
				status2=0;

			}
			f_lower_prev= f_lower;
			f_upper_prev= f_upper;
			if(gsl_integration){
				snr = snr_threshold_subroutine(	f_lower, f_upper, rel_err,params, generation_method,SN, w,np);
			}
			else{
				gauleg(log10(f_lower), log10(f_upper), GL_freqs, GL_w, GL_length);
				for(int k =0; k<GL_length; k++){
					GL_freqs[k] = pow(10.,GL_freqs[k]);
				}
				snr = calculate_snr(SN, "LISA",generation_method, params,GL_freqs,GL_length,"GAUSSLEG",GL_w,true);
			}
			//SNR found
			if(snr>SNR_thresh){not_found=false;}
			//Squeeze boundary
			else{
				//t1 moved last, so if SNR is larger (success), 
				//we move t1 down to accept change
				//if SNR is smaller (failure), we move t2 up
				if(t1_moved){
					if(snr>snr_prev){t1=t_mer;t1_moved=true; t2_moved=false;}
					else{t2=t_mer;t2_moved=true; t1_moved=false;}
				}
				//t2 moved last, so if SNR is larger (success), 
				//we move t2 up to accept change
				//if SNR is smaller (failure), we move t1 up
				if(t2_moved){
					if(snr>snr_prev){t2=t_mer;t2_moved=true; t1_moved=false;}
					else{t1=t_mer;t1_moved=true; t2_moved=false;}

				}
			}
			snr_prev = snr;

		}
	}	
	//If no SNR is larger than threshold, return negative
	if(not_found){
		threshold_times_out[0] = -1;
		threshold_times_out[1] = -1;
	}
	//Find roots
	else{
		double t_save = t_mer;
		bool found_lower_root=false;	
		bool found_upper_root=false;	
		double t1_prev=t_save, t2_prev=t_save;
		//We know t_mer is over the threshold 
		//(either initial snr was large enough or the boinitial search was successful, 
		//setting t_mer to a value that satisfies the threshold),
		//so start there
		t1=t_save, t2=t_save;
		//Find a extreme bound for bisection search
		//move t1 closer to merger until SNR is too low
		while(snr>SNR_thresh){
			t1/=2.;
			//std::cout<<"t1 | T_obs: "<<t1<<" | "<<T_obs<<std::endl;
			if(stellar_mass){
				f_lower =  f_0PN(t1 ,chirpmass) ;
				if(t1-T_obs > 0){
					f_upper =  std::min( f_0PN(  t1 - T_obs,chirpmass) , fmax ) ;
				}
				else{
					f_upper =   fmax;

				}
			}
			else{
				status1=Tbm_to_freq(params, generation_method, t1, &f_lower,tolerance,autodiff,max_iter,relative_time);
				
				if(status1 != 0){
					forced_PN_approx = true;
					f_lower = (f_0PN(t1, chirpmass));

				}
				if(t1-T_obs > 0){
					status2=Tbm_to_freq(params, generation_method, t1-T_obs, &f_upper,tolerance,autodiff,max_iter,relative_time);
					if(status2 != 0){
						forced_PN_approx = true;
						f_upper = (f_0PN(t1-T_obs, chirpmass));

					}
					f_upper = std::min(f_upper,fmax);
				}
				else{
					f_upper =   2*fpeak;
					//std::cout<<"Fpeak: "<<fpeak<<std::endl;
					f_upper = std::min(f_upper,fmax);
					//std::cout<<"Assigning f_upper: "<<f_upper<<std::endl;

				}
				//std::cout<<"f_lower | f_upper: "<<f_lower<<" | "<<f_upper<<std::endl;
				if(f_lower<0)
				{
					forced_PN_approx = true;
					f_lower = (f_0PN(t1, chirpmass));
				}
				if(f_upper<0){
					forced_PN_approx = true;
					f_upper =  std::min( f_0PN(  t1 - T_obs,chirpmass) , fmax );
				}
				if(f_lower>f_upper){
					forced_PN_approx = true;
					//f_lower = (f_0PN(t_mer, chirpmass));
					//f_upper =  std::min( f_0PN(  t_mer - T_obs,chirpmass) , fmax );
					f_lower = (f_0PN(t1, chirpmass));
					f_upper =  std::min( f_0PN(  t1 - T_obs,chirpmass) , fmax );
				}
				if(status1 !=0 || status2 != 0){
					bad_time_freq_inversion = true;
				}
				status1=0;
				status2=0;
			}
			//std::cout<<"Freq range: "<<f_lower<<" "<<f_upper<<std::endl;
			if(gsl_integration){
				snr = snr_threshold_subroutine(	f_lower, f_upper, rel_err,params, generation_method,SN, w,np);
			}
			else{
				gauleg(log10(f_lower), log10(f_upper), GL_freqs, GL_w, GL_length);
				for(int k =0; k<GL_length; k++){
					GL_freqs[k] = pow(10.,GL_freqs[k]);
				}
				snr = calculate_snr(SN, "LISA",generation_method, params,GL_freqs,GL_length,"GAUSSLEG",GL_w,true);
			}
			//If hit 5 day limit, exit loop with snr less than threshold, 
			//set lower bound to this time, and skip the rest of the search
			if(t1<5*T_day && snr>SNR_thresh){
				threshold_times_out[0]=t1;
				found_lower_root=true;
				break;
			}
		}
		//Search for exact time when SNR drops below threshold
		search_ct = 0;
		while(!found_lower_root && search_ct < max_iter_search){
			search_ct++;
			t_mer = (t1+t2)/2.;

			//std::cout<<"Merger time: "<<t_mer<<std::endl;
			//std::cout<<"Merger time diff: "<<t_mer-T_obs<<std::endl;
			//Find new frequency bound ids
			if(stellar_mass){
				f_lower =  f_0PN(t_mer ,chirpmass) ;
				if(t_mer-T_obs > 0){
					f_upper =  std::min( f_0PN(  t_mer - T_obs,chirpmass) , fmax);
				}
				else{
					f_upper = fmax;

				}
			}
			else{
				status1 = Tbm_to_freq(params, generation_method, t_mer, &f_lower,tolerance,autodiff,max_iter,relative_time);
				if(status1 != 0){
					forced_PN_approx = true;
					f_lower = (f_0PN(t_mer, chirpmass));

				}
				//std::cout<<"Value checking: "<<f_lower<<std::endl;
				if(t_mer-T_obs > 0){
					status2 = Tbm_to_freq(params, generation_method, t_mer-T_obs, &f_upper,tolerance,autodiff,max_iter,relative_time);
					if(status2 != 0){
						forced_PN_approx = true;
						f_upper = (f_0PN(t_mer-T_obs, chirpmass));

					}
					//std::cout<<"Value checking: "<<f_upper<<std::endl;
					f_upper = std::min(f_upper,fmax);
				}
				else{
					f_upper =   2*fpeak;
					//std::cout<<"Value checking: "<<f_upper<<std::endl;

				}
				if(f_lower<0)
				{
					forced_PN_approx = true;
					f_lower = (f_0PN(t_mer, chirpmass));
				}
				if(f_upper<0){
					forced_PN_approx = true;
					f_upper =  std::min( f_0PN(  t_mer - T_obs,chirpmass) , fmax );
				}
				if(f_lower>f_upper){
					forced_PN_approx = true;
					f_lower = (f_0PN(t_mer, chirpmass));
					f_upper =  std::min( f_0PN(  t_mer - T_obs,chirpmass) , fmax );
					//std::cout<<"Value checking (times): "<<t_mer-T_obs<<std::endl;
					//std::cout<<"Value checking: "<<f_upper<<std::endl;
				}
				if(status1 !=0 || status2 != 0){
					bad_time_freq_inversion = true;
				}
				status1=0;
				status2=0;

			}
			//std::cout<<"YIKES: Freq range: "<<f_lower<<" "<<f_upper<<std::endl;
			if(gsl_integration){
				snr = snr_threshold_subroutine(	f_lower, f_upper, rel_err,params, generation_method,SN, w,np);
			}
			else{
				gauleg(log10(f_lower), log10(f_upper), GL_freqs, GL_w, GL_length);
				for(int k =0; k<GL_length; k++){
					GL_freqs[k] = pow(10.,GL_freqs[k]);
				}
				snr = calculate_snr(SN, "LISA",generation_method, params,GL_freqs,GL_length,"GAUSSLEG",GL_w,true);
			}
			//SNR reached threshold
			if(std::abs(snr-SNR_thresh)/SNR_thresh<tolerance ){found_lower_root=true;threshold_times_out[0]=t_mer;}
			//If the bounds are too close, exit to stop rounding error
			else if( (fabs(t1-t2)*2/fabs(t1+t2)<1e-14) ){
				found_lower_root=true;
				threshold_times_out[0]=t_mer;	
				round_off_error = true;
			}
			else{
				//Sometimes, integration messes up
				if(snr==0){
					//If the SNR hits 0 for the lower root, its typically because the 
					//Waveform is sharply dropping off post merger, and this algorithm can't work
					//with such a steep function. But if this is the case, the error is on the order of 
					//seconds, which is completely negligible for most things this routine is used for
					found_lower_root=true; threshold_times_out[0]=(t1+t2 )/2.;
					
				}
				else{
					//If SNR is too high, move t2 to t_mer
					if(snr>SNR_thresh){ t2 = t_mer;	}
					//else, move t1 up 
					else{ t1=t_mer;}
				}
			}
			f_lower_prev= f_lower;
			f_upper_prev= f_upper;
			t1_prev= t1;
			t2_prev= t2;
			snr_prev=snr;
		}
		if(search_ct == max_iter_search ){
			max_iterations_reached = true;
		}
		//Do it again for lower bound
		t1=t_save; t2=t_save;
		//Initial bracketing
		//debugger_print(__FILE__,__LINE__,"Starting bracket loop");
		//debugger_print(__FILE__,__LINE__,"Parameters");
		//debugger_print(__FILE__,__LINE__,params->mass1);
		//debugger_print(__FILE__,__LINE__,params->mass2);
		//debugger_print(__FILE__,__LINE__,params->Luminosity_Distance);
		do{
			t2*=2.;
			if(t2>T_wait){t2=T_wait;}	
			if(stellar_mass){
				f_lower =  f_0PN(t2 ,chirpmass) ;
				if(t2-T_obs > 0){
					f_upper =  std::min( f_0PN(  t2 - T_obs,chirpmass) , fmax ) ;
				}
				else{
					f_upper =fmax;

				}
			}
			else{
				status1=Tbm_to_freq(params, generation_method, t2, &f_lower,tolerance,autodiff,max_iter,relative_time);
				
				if(status1 != 0){
					forced_PN_approx = true;
					f_lower = (f_0PN(t2, chirpmass));

				}
				if(t2-T_obs > 0){
					status2=Tbm_to_freq(params, generation_method, t2-T_obs, &f_upper,tolerance,autodiff,max_iter,relative_time);
					if(status2 != 0){
						forced_PN_approx = true;
						f_upper = (f_0PN(t2-T_obs, chirpmass));

					}
					f_upper = std::min(f_upper,fmax);
				}
				else{
					f_upper =   2*fpeak;

				}
				if(f_lower<0)
				{
					forced_PN_approx = true;
					f_lower = (f_0PN(t2, chirpmass));
				}
				if(f_upper<0){
					forced_PN_approx = true;
					f_upper =  std::min( f_0PN(  t2 - T_obs,chirpmass) , fmax );
				}
				if(f_lower>f_upper){
					forced_PN_approx = true;
					f_lower = (f_0PN(t_mer, chirpmass));
					f_upper =  std::min( f_0PN(  t_mer - T_obs,chirpmass) , fmax );
				}
				if(status1 !=0 || status2 != 0){
					bad_time_freq_inversion = true;
				}
				status1=0;
				status2=0;
			}
			//debugger_print(__FILE__,__LINE__,"Fmin");
			//debugger_print(__FILE__,__LINE__,f_lower);
			//debugger_print(__FILE__,__LINE__,"Fmax");
			//debugger_print(__FILE__,__LINE__,f_upper);
			
			if(gsl_integration){
				snr = snr_threshold_subroutine(	f_lower, f_upper, rel_err,params, generation_method,SN, w,np);
			}
			else{
				gauleg(log10(f_lower), log10(f_upper), GL_freqs, GL_w, GL_length);
				for(int k =0; k<GL_length; k++){
					GL_freqs[k] = pow(10.,GL_freqs[k]);
				}
				snr = calculate_snr(SN, "LISA",generation_method, params,GL_freqs,GL_length,"GAUSSLEG",GL_w,true);
			}
			if( (fabs(t2-T_wait)/T_wait < DOUBLE_COMP_THRESH) && snr>SNR_thresh){ found_upper_root=true; threshold_times_out[1]=T_wait;break;}
		}while(snr>SNR_thresh  );
		//debugger_print(__FILE__,__LINE__,"Finished bracket loop");
		//debugger_print(__FILE__,__LINE__,"Exit time");
		//debugger_print(__FILE__,__LINE__,t2/T_year);
		//debugger_print(__FILE__,__LINE__,"SNR");
		//debugger_print(__FILE__,__LINE__,snr);
		//if(found_upper_root){
		//	debugger_print(__FILE__,__LINE__,"Found upper root");
		//}
		search_ct = 0;
		double t2_save_output = t2;
		while(!found_upper_root && search_ct < max_iter_search){
			
			t_mer = (t1+t2)/2.;
			//Find new frequency bound ids
			if(stellar_mass){
				f_lower =  f_0PN(t_mer ,chirpmass) ;
				if(t_mer-T_obs > 0){
					f_upper =  std::min( f_0PN(  t_mer - T_obs,chirpmass) , fmax );
				}
				else{
					f_upper = fmax;

				}
			}
			else{
				status1=Tbm_to_freq(params, generation_method, t_mer, &f_lower,tolerance,autodiff,max_iter,relative_time);
				
				if(status1 != 0){
					forced_PN_approx = true;
					f_lower = (f_0PN(t_mer, chirpmass));

				}
				if(t_mer-T_obs > 0){
					status2=Tbm_to_freq(params, generation_method, t_mer-T_obs, &f_upper,tolerance,autodiff,max_iter,relative_time);
					if(status2 != 0){
						forced_PN_approx = true;
						f_upper = (f_0PN(t_mer-T_obs, chirpmass));

					}
					f_upper = std::min(f_upper,fmax);
				}
				else{
					f_upper =   2*fpeak;

				}
				if(f_lower<0)
				{
					forced_PN_approx = true;
					f_lower = (f_0PN(t_mer, chirpmass));
				}
				if(f_upper<0){
					forced_PN_approx = true;
					f_upper =  std::min( f_0PN(  t_mer - T_obs,chirpmass) , fmax );
				}
				if(f_lower>f_upper){
					forced_PN_approx = true;
					f_lower = (f_0PN(t_mer, chirpmass));
					f_upper =  std::min( f_0PN(  t_mer - T_obs,chirpmass) , fmax );
				}
				if(status1 !=0 || status2 != 0){
					bad_time_freq_inversion = true;
				}
				status1=0;
				status2=0;
			}
			f_lower_prev= f_lower;
			f_upper_prev= f_upper;
			if(gsl_integration){
				snr = snr_threshold_subroutine(	f_lower, f_upper, rel_err,params, generation_method,SN, w,np);
			}
			else{
				gauleg(log10(f_lower), log10(f_upper), GL_freqs, GL_w, GL_length);
				for(int k =0; k<GL_length; k++){
					GL_freqs[k] = pow(10.,GL_freqs[k]);
				}
				snr = calculate_snr(SN, "LISA",generation_method, params,GL_freqs,GL_length,"GAUSSLEG",GL_w,true);
			}
			//std::cout<<snr<<std::endl;
			if(std::abs(snr-SNR_thresh)/SNR_thresh<tolerance ){found_upper_root=true;threshold_times_out[1]=t_mer;}
			else if( (fabs(t1-t2)*2/fabs(t1+t2)<1e-14) ){
				found_upper_root=true;
				threshold_times_out[1]=t_mer;	
				round_off_error = true;
			}
			else{
				if(snr==0){
					threshold_times_out[0] = -1;
					threshold_times_out[1] = -1;
					params->sky_average = save_SA;
					return 12;
					
				}
				else{
					if(snr>SNR_thresh){ t1 = t_mer;	}
					else{ t2=t_mer;}
				}
			}
			snr_prev=snr;
			
		}
		if(search_ct == max_iter_search ){
			max_iterations_reached = true;
		}
		//std::cout<<t2_save_output/T_year<<" "<<threshold_times_out[1]/T_year<<" "<<*(T_obs_SNR)<<std::endl;
	}
	
	params->sky_average = save_SA;
	if(round_off_error){
		return 13;
	}
	if(bad_time_freq_inversion){
		return 14;
	}
	if(forced_PN_approx){
		return 15;
	}
	if(max_iterations_reached){
		return 16;
	}

	if(!gsl_integration){
		delete [] GL_freqs ;
		delete [] GL_w ;
	}
	return 0;
}
double snr_threshold_subroutine(double fmin, double fmax, double rel_err, gen_params_base<double> *params, std::string generation_method,std::string SN, gsl_integration_workspace *w, int np)
{
	//std::cout<<fmin<<" "<<fmax<<std::endl;
	gsl_snr_struct helper_params;
	helper_params.params = params;
	helper_params.generation_method = generation_method;
	helper_params.SN = SN;
	gsl_function F;
	F.function = [](double f, void *param){return integrand_snr_SA_subroutine(f,param);};
	F.params = (void *)&helper_params;
	double result, err;
	int errcode = gsl_integration_qag(&F,fmin, fmax, 0,rel_err, np, GSL_INTEG_GAUSS15,w, &result, &err);
	return sqrt(result);
}

double TaylorF2ReducedSpinChirpTime(
	const double fStart,	//< Starting GW frequency in Hertz
	const double m1,	//< Primary mass in solar units
	const double m2,	//< Secondary mass in solar units
	const double s1z,	//< Dimension-less primary aligned spin component
	const double s2z,	//< Dimension-less secondary aligned spin component
	const int PNO		//< Twice the PN phase order
)
{
    const double m = m1+m2; 
    const double eta = m1*m2/(m*m);
    const double eta2 = eta*eta;

	// Compute reduced-spin parameter (arXiv:1107.1267, Eq. (5.9))
	const double delta = (m1 - m2) / m;
	const double chi_s = (s1z + s2z) / 2.;
	const double chi_a = (s1z - s2z) / 2.;
	const double chi = chi_s * (1. - 76. * eta / 113.) + delta * chi_a;

    const double chi2 = chi*chi;
    const double sigma0 = (-12769*(-81. + 4.*eta))/(16.*(-113. + 76.*eta)*(-113. + 76.*eta));
    const double gamma0 = (565*(-146597. + 135856.*eta + 17136.*eta2))/(2268.*(-113. + 76.*eta));
 
	double v = cbrt(M_PI * m * GWAT_MTSUN_SI * fStart); // Initial dimension-less velocity
    double vk = v;  // v^k
    double tk[8];  // chirp time coefficients up to 3.5 PN
	const double Pi_p2 = M_PI*M_PI;

	// chirp time coefficients up to 3.5PN
    tk[0] = (5.*m*GWAT_MTSUN_SI)/(256.*pow_int(v,8)*eta);
    tk[1] = 0.;
    tk[2] = 2.9484126984126986 + (11*eta)/3.;
    tk[3] = (-32*M_PI)/5. + (226.*chi)/15.;
    tk[4] = 6.020630590199042 - 2*sigma0*chi2 + (5429*eta)/504. + (617*eta2)/72.;
    tk[5] = (3*gamma0*chi)/5. - (7729*M_PI)/252. + (13*M_PI*eta)/3.;
    tk[6] = -428.291776175525 + (128*Pi_p2)/3. + (6848*gamma_E)/105. + (3147553127*eta)/3.048192e6 -
       (451*Pi_p2*eta)/12. - (15211*eta2)/1728. + (25565*eta2*eta)/1296. + (6848*log(4*v))/105.;
    tk[7] = (-15419335*M_PI)/127008. - (75703*M_PI*eta)/756. + (14809*M_PI*eta2)/378.;

	// compute chirp time
    double tau = 1.;  // chirp time
    for (size_t k = 2; k <= PNO; k++) {
        vk *= v;
        tau += tk[k] * vk;
    }
    tau *= tk[0];

	return tau;
}

//###########################################################################
//template void map_extrinsic_angles<double>(gen_params_base<double> *);
//template void map_extrinsic_angles<adouble>(gen_params_base<adouble> *);

template void  time_phase_corrected<double>(double *, int, double *, gen_params_base<double> *, std::string, bool,int );
template void  time_phase_corrected<adouble>(adouble *, int, adouble *, gen_params_base<adouble> *, std::string, bool,int);

template int fourier_detector_response_horizon<double>(double *, int, waveform_polarizations<double> *,std::complex<double> *, double, double, std::string);
template int fourier_detector_response_horizon<adouble>(adouble *, int, waveform_polarizations<adouble> *,std::complex<adouble> *, adouble, adouble, std::string);
//
template int fourier_detector_response_horizon<double>(double *, int, waveform_polarizations<double> *, std::complex<double> *, double, double, double, std::string);
template int fourier_detector_response_horizon<adouble>(adouble *, int, waveform_polarizations<adouble> *, std::complex<adouble> *, adouble, adouble, adouble, std::string);
//
//
template int time_detector_response_horizon<double>(double *, int, std::complex<double> *, std::string, std::string, gen_params_base<double>*);
template int time_detector_response_horizon<adouble>(adouble *, int, std::complex<adouble> *, std::string, std::string, gen_params_base<adouble>*);
//
template int time_detector_response_equatorial<double>(double *, int , std::complex<double> *,std::string, std::string, gen_params_base<double> *);
template int time_detector_response_equatorial<adouble>(adouble *, int , std::complex<adouble> *,std::string, std::string, gen_params_base<adouble> *);
//
//
template int fourier_detector_response_equatorial<double>(double *, int , std::complex<double> *,std::string, std::string, gen_params_base<double> *);
template int fourier_detector_response_equatorial<adouble>(adouble *, int , std::complex<adouble> *,std::string, std::string, gen_params_base<adouble> *);
template int fourier_detector_response_equatorial<double>(double *, int , std::complex<double> *,std::string, std::string, gen_params_base<double> *, double *);
template int fourier_detector_response_equatorial<adouble>(adouble *, int , std::complex<adouble> *,std::string, std::string, gen_params_base<adouble> *, adouble *);
//
template int fourier_detector_response_horizon<double>(double *, int, std::complex<double> *, std::string, std::string, gen_params_base<double>*);
template int fourier_detector_response_horizon<adouble>(adouble *, int, std::complex<adouble> *, std::string, std::string, gen_params_base<adouble>*);
//
template int fourier_detector_response_equatorial<double>(double *, int, waveform_polarizations<double> *, std::complex<double> *, double, double , double, double,double *, double, double, double, double, std::string);
template int fourier_detector_response_equatorial<adouble>(adouble *, int, waveform_polarizations<adouble> *, std::complex<adouble> *, adouble, adouble , adouble, double,adouble*, adouble, adouble, adouble, adouble,  std::string);
//
