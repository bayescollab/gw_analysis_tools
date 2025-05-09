#ifndef MCMC_GW_EXTENDED_H
#define MCMC_GW_EXTENDED_H

#include "util.h"
#include "mcmc_gw.h"
#include <bayesship/bayesshipSampler.h>
#include <bayesship/dataUtilities.h>
#include "adaptivelikelihoods.h"

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
	AdaptiveLikelihood *adaptivell=nullptr;
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

//double MCMC_likelihood_wrapper_v2(bayesship::positionInfo *pos, int chainID, bayesship::bayesshipSampler *sampler ,void *userParameters);

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
	//double(*log_prior)(bayesship::positionInfo *pos, int chainID,bayesship::bayesshipSampler *sampler, void *userParameters),
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

#endif
