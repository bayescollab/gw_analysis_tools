#ifndef MCMC_GW_H
#define MCMC_GW_H
#include <complex>
#include <fftw3.h>
#include "util.h"
#include <iostream>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <mutex>
#include "mcmc_sampler.h"
#include "quadrature.h"
#include "adaptivelikelihoods.h"
/*! \file 
 *
 * Header file for the Graviational Wave specific MCMC routines
 */

//########################################
//MCMC global variables - facilitate wrapping 
//of likelihood function 
static double **mcmc_noise=NULL;
static double *mcmc_init_pos=NULL;
static std::complex<double> **mcmc_data=NULL;
static double **mcmc_frequencies=NULL;
static std::string *mcmc_detectors=NULL;
static std::string mcmc_generation_method="";
static std::string mcmc_generation_method_base="";
static std::string mcmc_generation_method_extended="";
static int *mcmc_data_length=NULL ;
static fftw_outline *mcmc_fftw_plans=NULL ;
static int mcmc_num_detectors=2;
static double mcmc_gps_time=0;
static double mcmc_gmst=0;
static int mcmc_max_dim;
static int mcmc_min_dim;
static int mcmc_Nmod;
static int mcmc_Nmod_max;
static double *mcmc_bppe;
static gsl_interp_accel **mcmc_accels = NULL;
static gsl_spline **mcmc_splines = NULL;
static bool mcmc_log_beta;
static bool mcmc_intrinsic;
static bool mcmc_save_waveform;
static int mcmc_deriv_order=4;
static gsl_rng **mcmc_rvec;

//########################################

struct MCMC_modification_struct
{
	int ppE_Nmod = 0; //ppE
	double *bppe = NULL; //ppE
	int gIMR_Nmod_phi = 0; //gIMR
	int *gIMR_phii = NULL; //gIMR
	int gIMR_Nmod_sigma = 0; //gIMR
	int *gIMR_sigmai = NULL; //gIMR
	int gIMR_Nmod_beta = 0; //gIMR
	int *gIMR_betai = NULL; //gIMR
	int gIMR_Nmod_alpha = 0; //gIMR
	int *gIMR_alphai = NULL; //gIMR

	bool NSflag1 =false;
	bool NSflag2 =false;

	bool tidal_love = true; 
        bool tidal_love_error = false;
        bool alpha_param = true;
        bool EA_region1 = false; 
  
	/* Whether to use Gauss-Legendre Quadrature for the LIKELIHOOD
 * 		If using GLQ, provide the weights vector for the integration */
	bool GAUSS_QUAD=false;
	bool log10F = false;
	double **weights=NULL;

	
	/*If set to anything besides NULL, this overrides the frequency/PSD arrays passed
 * 		to sampler.
 * 		These parameters are for the fishers ONLY*/
	double **fisher_freq=NULL;
	double **fisher_weights=NULL;
	double **fisher_PSD=NULL;
	bool fisher_GAUSS_QUAD=false;
	int *fisher_length =NULL;
	bool fisher_log10F = false;

	// Quadrature method to be used for Fisher and likelihood calculations
	// (unless an AdaptiveLikelihood is being used for the latter).
	Quadrature *QuadMethod = NULL;
	AdaptiveLikelihood *adaptivell = NULL;

	// Refererence frequency for waveform generation.
	// Not as trivial as you might think!
	// Default: 20 Hz.
	double f_ref = 20.;
};

static MCMC_modification_struct *mcmc_mod_struct;

struct MCMC_user_param
{
	std::complex<double> **burn_data=NULL;
	double **burn_noise=NULL;
	double **burn_freqs=NULL;
	int *burn_lengths=NULL;
	fftw_outline *burn_plans=NULL;
	std::mutex *mFish;

	bool GAUSS_QUAD=false;
	bool log10F = false;
	double **weights=NULL;

	double **fisher_freq=NULL;
	double **fisher_weights=NULL;
	double **fisher_PSD=NULL;
	bool fisher_AD=false;
	bool fisher_GAUSS_QUAD=false;
	int *fisher_length =NULL;
	bool fisher_log10F = false;

	Quadrature *QuadMethod = NULL;
		
	//RJ stuff
	double **mod_prior_ranges=NULL;
	MCMC_modification_struct *mod_struct;
	
	gsl_rng *r;
};


double maximized_coal_log_likelihood_IMRPhenomD(double *frequencies,
				int length,
				std::complex<double> *data,
				double *noise,
				double SNR,
				double chirpmass,	
				double symmetric_mass_ratio, 
				double spin1,
				double spin2,
				bool NSflag,
				fftw_outline *plan);
double maximized_coal_log_likelihood_IMRPhenomD(double *frequencies,
				size_t length,
				double *real_data,
				double *imag_data,
				double *noise,
				double SNR,
				double chirpmass,	
				double symmetric_mass_ratio, 
				double spin1,
				double spin2,
				bool NSflag);

double maximized_coal_log_likelihood_IMRPhenomD(double *frequencies,
				size_t length,
				double *real_data,
				double *imag_data,
				double *noise,
				double SNR,
				double chirpmass,	
				double symmetric_mass_ratio, 
				double spin1,
				double spin2,
				bool NSflag,
				fftw_outline *plan);
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
				fftw_outline *plan);

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
				bool NSflag);

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
				fftw_outline *plan);


double maximized_Log_Likelihood_aligned_spin_internal(std::complex<double> *data,
				double *psd,
				double *frequencies,
				std::complex<double> *detector_response,
				size_t length,
				fftw_outline *plan
				);
double Log_Likelihood(std::complex<double> *data,
				double *psd,
				double *frequencies,
				size_t length,
				gen_params_base<double> *params,
				std::string detector,
				std::string generation_method
				);

double maximized_Log_Likelihood_unaligned_spin_internal(std::complex<double> *data,
				double *psd,
				double *frequencies,
				std::complex<double> *hplus,
				std::complex<double> *hcross,
				size_t length,
				fftw_outline *plan
				);
double maximized_Log_Likelihood(std::complex<double> *data, 
				double *psd,
				double *frequencies,
				size_t length,
				gen_params_base<double> *params,
				std::string detector,
				std::string generation_method,
				fftw_outline *plan
				);
double maximized_Log_Likelihood(double *data_real, 
				double *data_imag,
				double *psd,
				double *frequencies,
				size_t length,
				//gen_params *params,
				gen_params_base<double> *params,
				std::string detector,
				std::string generation_method,
				fftw_outline *plan
				);
double maximized_coal_Log_Likelihood(double *data_real, 
				double *data_imag,
				double *psd,
				double *frequencies,
				size_t length,
				double *template_real,
				double *template_imag,
				fftw_outline *plan
				);
double maximized_coal_Log_Likelihood(std::complex<double> *data, 
				double *psd,
				double *frequencies,
				size_t length,
				gen_params_base<double> *params,
				std::string detector,
				std::string generation_method,
				fftw_outline *plan,
				double *tc,
				double *phic
				);
double maximized_coal_Log_Likelihood_internal(std::complex<double> *data,
				double *psd,
				double *frequencies,
				std::complex<double> *detector_response,
				size_t length,
				fftw_outline *plan,
				double *tc,
				double *phic
				);
double Log_Likelihood_internal(std::complex<double> *data,
			double *psd,
			double *frequencies,
			double *weights,
			std::complex<double> *detector_response,
			int length,
			bool log10F,
			std::string integration_method
			);
double Log_Likelihood_internal(std::complex<double> *data,
	double *psd,
	std::complex<double> *detector_response,
	Quadrature *QuadMethod
);

double MCMC_likelihood_wrapper_SKYSEARCH(double *param, mcmc_data_interface *interface,void *parameters);
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
	);

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
	std::string *detector,
	MCMC_modification_struct *mod_struct,
	std::string generation_method,
	std::string statistics_filename,
	std::string chain_filename,
	std::string auto_corr_filename,
	std::string likelihood_log_filename,
	std::string checkpoint_filename
	);
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
	);
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
	);
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
	std::string likelihood_log_filename,
	std::string checkpoint_filename
	);
void continue_PTMCMC_MH_GW(std::string start_checkpoint_file,
	double ***output,
	int dimension,
	int N_steps,
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
	std::string *detector,
	MCMC_modification_struct *mod_struct,
	std::string generation_method,
	std::string statistics_filename,
	std::string chain_filename,
	std::string auto_corr_filename,
	std::string likelihood_log_filename,
	std::string final_checkpoint_filename,
	bool tune
	);
void PTMCMC_method_specific_prep(std::string generation_method, int dimension,double **seeding_var_ptr, bool local_seeding);

void RJPTMCMC_method_specific_prep(std::string generation_method, int max_dim, int min_dim,double *seeding_var, bool local_seeding);

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
	int num_detectors,
	Quadrature *QuadMethod = NULL
	);

void MCMC_fisher_wrapper(double *param,  double **output, mcmc_data_interface *interface ,void *parameters) ;

std::string MCMC_prep_params(double *param, double *temp_params, gen_params_base<double> *gen_params, int dimension, std::string generation_method,MCMC_modification_struct *mod_struct);

double MCMC_likelihood_wrapper(double *param, mcmc_data_interface *interface ,void *parameters) ;

void continue_RJPTMCMC_MH_dynamic_PT_alloc_comprehensive_2WF_GW(
	std::string checkpoint_file_start,
	mcmc_sampler_output *sampler_output,
	double **output,
	int **status,
	int *model_status,
	int nested_model_number,
	int N_steps,
	int max_chain_N_thermo_ensemble,
	double **prior_ranges,
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
	);
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
	double **prior_ranges,
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
	);
void RJMCMC_2WF_RJ_proposal_wrapper(
	double *current_param, 
	double *proposed_param, 
	int *current_status, 
	int *proposed_status, 
	int *current_model_status, 
	int *proposed_model_status, 
	double *MH_corrections,
	mcmc_data_interface *interface, 
	void *user_param);
void RJMCMC_2WF_fisher_wrapper(
	double *param, 
	int *status, 
	int model_status, 
	double **fisher,
	mcmc_data_interface *interface, 
	void *user_param);
double RJMCMC_2WF_likelihood_wrapper(
	double *param, 
	int *status, 
	int model_status, 
	mcmc_data_interface *interface, 
	void *user_param);
void pack_local_mod_structure(mcmc_data_interface *interface,
	double *param,
	int *status,
	std::string waveform_extended, 
	void *parameters, 
	MCMC_modification_struct *full_struct, 
	MCMC_modification_struct *local_struct );

void MCMC_fisher_transformations(
	double *param, 
	double **fisher, 
	int dimension,
	std::string generation_method,
	bool intrinsic,
	mcmc_data_interface *interface, 
	MCMC_modification_struct *mod_struct,
	void *parameters);

#endif
