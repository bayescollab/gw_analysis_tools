#ifndef MCMC_GW_EXTENDED_H
#define MCMC_GW_EXTENDED_H

#include "util.h"
#include "mcmc_gw.h"
#include <bayesship/bayesshipSampler.h>
#include <bayesship/dataUtilities.h>
#include "likelihoods.h"

struct mcmcVariables 
{
	double **mcmc_noise=nullptr;
	double *mcmc_init_pos=nullptr;
	std::complex<double> **mcmc_data=nullptr;
	double **mcmc_frequencies=nullptr;
	std::string *mcmc_detectors=nullptr;
	std::string mcmc_generation_method="";
	int *mcmc_data_length = nullptr;
	fftw_outline *mcmc_fftw_plans=nullptr;
	int mcmc_num_detectors  ;
	double mcmc_gps_time = 0;
	double mcmc_gmst = 0;
	int mcmc_max_dim;
	int mcmc_min_dim;
	MCMC_modification_struct *mcmc_mod_struct=nullptr;
	bool mcmc_save_waveform=true;	
	bool mcmc_intrinsic=false;	
	int mcmc_deriv_order = 4;
	bool mcmc_log_beta = false;
	MCMC_user_param *user_parameters = nullptr;	
	double maxDim;
	bool mcmc_adaptive = false;
	GWATLikelihoods::Likelihood *adaptivell=nullptr;
	Quadrature *QuadMethod = NULL;
};
struct mcmcVariablesRJ
{
	double **mcmc_noise=nullptr;
	double *mcmc_init_pos=nullptr;
	std::complex<double> **mcmc_data=nullptr;
	double **mcmc_frequencies=nullptr;
	std::string *mcmc_detectors=nullptr;
	std::string mcmc_generation_method="";
	std::string mcmc_generation_method_extended="";
	int *mcmc_data_length = nullptr;
	fftw_outline *mcmc_fftw_plans=nullptr;
	int mcmc_num_detectors  ;
	double mcmc_gps_time = 0;
	double mcmc_gmst = 0;
	int mcmc_max_dim;
	int mcmc_min_dim;
	MCMC_modification_struct *mcmc_mod_struct=nullptr;
	bool mcmc_save_waveform=true;	
	bool mcmc_intrinsic=false;	
	int mcmc_deriv_order = 4;
	bool mcmc_log_beta = false;
	MCMC_user_param *user_parameters = nullptr;	
	int maxDim;
	int minDim;
};

void PTMCMC_method_specific_prep_v2(std::string generation_method, int dimension, bool *intrinsic,MCMC_modification_struct *mod_struct);
void RJPTMCMC_method_specific_prep_v2(std::string generation_method, int dimension, bool *intrinsic,MCMC_modification_struct *mod_struct);

std::string MCMC_prep_params_v2(double *param, double *temp_params, gen_params_base<double> *gen_params, int dimension, std::string generation_method, MCMC_modification_struct *mod_struct, bool intrinsic, double gmst);


bayesship::bayesshipSampler *  PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW_v2(
	int dimension,
	int independentSamples,
	int ensembleSize,
	int ensembleN,
	bayesship::positionInfo *initialPosition,
	bayesship::positionInfo **initialEnsemble,
	double swapProb,
	int burnIterations,
	int burnPriorIterations,
	int priorIterations,
	bool writePriorData,
	int batchSize,
	double **priorRanges,
	bayesship::probabilityFn *lp,
	int numThreads,
	bool pool,
	int num_detectors,
	std::complex<double> **data,
	double **noise_psd,
	double **frequencies,
	int *data_length,
	double gps_time,
	std::string *detectors,
	MCMC_modification_struct *mod_struct,
	std::string generation_method,
	std::string outputDir,
	std::string outputFileMoniker, 
	bool ignoreExistingCheckpoint,
	bool restrictSwapTemperatures,
	bool coldChainStorageOnly);

bayesship::bayesshipSampler *  RJPTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW_v2(
	int minDim,
	int maxDim,
	int independentSamples,
	int ensembleSize,
	int ensembleN,
	bayesship::positionInfo *initialPosition,
	bayesship::positionInfo **initialEnsemble,
	double swapProb,
	int burnIterations,
	int burnPriorIterations,
	int priorIterations,
	bool writePriorData,
	int batchSize,
	double **priorRanges,
	bayesship::probabilityFn *log_prior,
	int numThreads,
	bool pool,
	int num_detectors,
	std::complex<double> **data,
	double **noise_psd,
	double **frequencies,
	int *data_length,
	double gps_time,
	std::string *detectors,
	MCMC_modification_struct *mod_struct,
	std::string generation_method,
	std::string generation_method_extended,
	std::string outputDir,
	std::string outputFileMoniker,
	bool ignoreExistingCheckpoint,
	bool restrictSwapTemperatures,
	bool coldChainStorageOnly);

/// Runs a  Metropolis-Hastings chain from initial_params to locate
/// a good fiducial waveform for relative binning initialization.
/// Outputs detector responses:
///   fiducial_out  — response at the maximum-likelihood point found
///   test_out      — response at the last accepted point in the chain
/// Caller must pre-allocate both output arrays as
/// [num_detectors][data_length[i]]. If waveform generation fails for
/// fiducial_out, it is set to the initial parameter point. If it fails for
/// test_out, it is set to the max-likelihood point.
///
/// \param dimension      Size of parameter space.
/// \param initial_params Array of parameters to begin exploration.
/// \param log_prior      Log-prior probability function for the model.
/// \param param_scales   Typical scale for each parameter (e.g. prior width);
///  used to size finite-difference steps for the Fisher diagonal.
/// \param num_mh_steps   Number of M-H steps to be taken.
/// \param num_detectors  Number of detectors.
/// \param data           The signal at each detector whose fiducial waveform
/// will be found for.
/// \param noise_psd      The PSD for each detector.
/// \param frequencies    Frequency arrays for each detector.
/// \param data_length    Length of the signal for each detector.
/// \param gps_time       GPS time of the signal.
/// \param detectors      Names of the detectors.
/// \param generation_method Model name for the output waveforms.
/// \param mod_struct     MCMC_modification_struct
/// \param fiducial_out  [out] The fiducial signal.
/// \param test_out	     [out] The signal at the final M-H step.
/// \param test_gp_out   [out] The parameters at the final M-H step.
void find_fiducial(int dimension, double* initial_params,
                   bayesship::probabilityFn* log_prior,
                   const std::vector<double>& param_scales, int num_mh_steps,
                   int num_detectors, std::complex<double>** data,
                   double** noise_psd, double** frequencies, int* data_length,
                   double gps_time, std::string* detectors,
                   std::string generation_method,
                   MCMC_modification_struct* mod_struct,
                   std::complex<double>** fiducial_out,
                   std::complex<double>** test_out,
                   gen_params_base<double>* test_gp_out = nullptr);

#endif
