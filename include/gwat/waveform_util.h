#ifndef WAVEFORM_UTIL_H
#define WAVEFORM_UTIL_H
#include "waveform_generator.h"
#include "util.h"
#include <string>
#include <gsl/gsl_integration.h>

/*! \file 
 * Header file for waveform specific utilites
 */

double data_snr_maximized_extrinsic(double *frequencies,
	int length,
	std::complex<double> *data,
	double *psd,
	std::string detector,
	std::string generation_method,
	gen_params *param
	);
double data_snr_maximized_extrinsic(double *frequencies,
	int length,
	double *data_real,
	double *data_imag,
	double *psd,
	std::string detector,
	std::string generation_method,
	gen_params *param
	);
double calculate_snr(std::string sensitivity_curve,
        std::complex<double> *waveform,
        double *frequencies,
        int length);
double calculate_snr_internal(double *psd,
        std::complex<double> *waveform,
        double *frequencies,
        int length);
double calculate_snr(std::string sensitivity_curve,
	std::string detector,
	std::string generation_method,
        gen_params_base<double> *params,
        double *frequencies,
        int length);
double calculate_snr_gsl(std::string sensitivity_curve,
	std::string detector,
	std::string generation_method,
	gen_params_base<double> *params,
	double f_min,
	double f_max,
	double relative_error
	);
double calculate_snr_gsl(std::string sensitivity_curve,
	std::string detector,
	std::string generation_method,
	gen_params_base<double> *params,
	double f_min,
	double f_max,
	double relative_error,
	gsl_integration_workspace *w, 
	int np
	);

double integrand_snr_SA_subroutine(double f, void *subroutine_params);
double integrand_snr_subroutine(double f, void *subroutine_params);
template<class T>
int fourier_detector_response_horizon(T *frequencies, 
	int length,
	std::complex<T> *hplus, 
	std::complex<T> *hcross, 
	std::complex<T> *detector_response, 
	T theta, 
	T phi, 
	std::string detector);
template<class T>
int fourier_detector_response_horizon(T *frequencies, 
	int length,
	std::complex<T> *hplus, 
	std::complex<T> *hcross, 
	std::complex<T> *detector_response, 
	T theta, 
	T phi, 
	T psi, 
	std::string detector);
template<class T>
int fourier_detector_response_horizon(T *frequencies, 
	int length,
	std::complex<T> *response, 
	std::string detector,
	std::string generation_method,
	gen_params_base<T> *parameters
	);
template<class T>
int fourier_detector_response_equatorial(T *frequencies, 
	int length,
	std::complex<T> *hplus, 
	std::complex<T> *hcross, 
	std::complex<T> *detector_response, 
	T ra, 
	T dec, 
	T psi, 
	double gmst, 
	T *times,
	T LISA_alpha0,
	T LISA_phi0,
	T theta_j_ecl,
	T phi_j_ecl,
	std::string detector);

template<class T>
int fourier_detector_response_equatorial(T *frequencies, 
	int length,
	std::complex<T> *response, 
	std::string detector,
	std::string generation_method,
	gen_params_base<T> *parameters
	);
template<class T>
int fourier_detector_response_equatorial(T *frequencies, 
	int length,
	std::complex<T> *response, 
	std::string detector,
	std::string generation_method,
	gen_params_base<T> *parameters,
	T *times
	);

template<class T>
int fourier_detector_response(T *frequencies,
	int length,
	std::complex<T> *response,
	std::string detector,
	std::string generation_method,
	gen_params_base<T> *parameters,
	T *times=NULL);


int boundary_number(std::string method);

void time_phase_corrected_autodiff(double *times, 
	int length, 
	double *frequencies,
	gen_params_base<double> *params, 
	std::string generation_method, 
	bool correct_time,
	int *tapes_in = NULL);

template<class T>
void time_phase_corrected(T *times, 
	int length, 
	T *frequencies, 
	gen_params_base<T> *params,
	std::string generation_method,
	bool correct_time
	);
int fourier_detector_amplitude_phase(double *frequencies, 
	int length,
	double *amplitude, 
	double *phase, 
	std::string detector,
	std::string generation_method,
	gen_params *parameters
	);
template<class T>
void transform_orientation_coords(gen_params_base<T> *parameters,std::string generation_method,std::string detector);
template<class T>
void map_extrinsic_angles(gen_params_base<T> *params);

void assign_freq_boundaries(double *freq_boundaries, 
	double *intermediate_freqs, 
	int boundary_num, 
	gen_params_base<double> *input_params, 
	std::string generation_method);


void integration_bounds(gen_params_base<double> *params, 
	std::string generation_method,
	std::string detector, 
	std::string sensitivity_curve, 
	double fmin, 
	double fmax, 
	double signal_to_noise,
	double tol,
	double *integration_bounds
	) ;
void integration_interval(double sampling_freq, 
	double integration_time, 
	std::string detector, 
	std::string sensitivity_curve, 
	std::string generation_method,
	gen_params_base<double> *params,
	double *freq_bounds);
void Tbm_to_freq(gen_params_base<double> *params,
	std::string generation_method,
	double Tbm,
	double *freq,
	double tol 
	);
template<class T>
void postmerger_params(gen_params_base<T>*params,
	std::string generation_method,
	T *fpeak,
	T *fdamp,
	T *fRD
	);
void threshold_times(gen_params_base<double> *params,
	std::string generation_method,
	double T_obs, 
	double T_wait,
	double f_lower,
	double f_upper,
	std::string SN,
	double SNR_thresh, 
	double *threshold_times_out,
	double tolerance
	);
void threshold_times(gen_params_base<double> *params,
	std::string generation_method,
	double T_obs, 
	double T_wait, 
	double *freqs,
	double *SN,
	int length,
	double SNR_thresh, 
	double *threshold_times_out,
	double tolerance
	);
double integrand_threshold_subroutine(double f, void *subroutine_params);
double snr_threshold_subroutine(double fmin, double fmax, double rel_err, gen_params_base<double> *params, std::string generation_method,std::string SN, gsl_integration_workspace *w, int np);
int threshold_times_gsl(gen_params_base<double> *params,
	std::string generation_method, 
	double T_obs, 
	double T_wait, 
	double fmin,
	double fmax,
	std::string SN,
	double SNR_thresh, 
	double *threshold_times_out,
	double tolerance, 
	gsl_integration_workspace *w,
	int np
	);
#endif

