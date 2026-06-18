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
        bool EOS_plat_flag = false; 

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

	// ---- Temperature ladder settings ----
	// The maximum temperature below Inf
	double Tmax = 100.;
	// The dampening time. If left at 0, BayesShip sets it to burnIterations/8.
	double t0 = 0.;
	// Timescale of adjustments
	double nu = 100.;
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



#endif
