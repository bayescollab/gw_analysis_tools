#ifndef MCMC_GW_H
#define MCMC_GW_H
/*! \file
 *
 * Header file for the Graviational Wave specific MCMC routines
 */

#include <bayesship/bayesshipSampler.h>
#include <bayesship/dataUtilities.h>
#include <fftw3.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_spline.h>

#include <algorithm>
#include <complex>
#include <condition_variable>
#include <fstream>
#include <iostream>
#include <mutex>
#include <thread>
#include <vector>

#include "likelihoods.h"
#include "mcmc_sampler.h"
#include "parameter_space.h"
#include "quadrature.h"
#include "util.h"

static bool mcmc_intrinsic;
static bool mcmc_save_waveform;

struct MCMC_modification_struct {
  int ppE_Nmod = 0;         // ppE
  double* bppe = NULL;      // ppE
  int gIMR_Nmod_phi = 0;    // gIMR
  int* gIMR_phii = NULL;    // gIMR
  int gIMR_Nmod_sigma = 0;  // gIMR
  int* gIMR_sigmai = NULL;  // gIMR
  int gIMR_Nmod_beta = 0;   // gIMR
  int* gIMR_betai = NULL;   // gIMR
  int gIMR_Nmod_alpha = 0;  // gIMR
  int* gIMR_alphai = NULL;  // gIMR

  bool NSflag1 = false;
  bool NSflag2 = false;
  bool EOS_plat_flag = false;

  bool tidal_love = true;
  bool tidal_love_error = false;
  bool alpha_param = true;
  bool EA_region1 = false;

  /* Whether to use Gauss-Legendre Quadrature for the LIKELIHOOD
   * 		If using GLQ, provide the weights vector for the integration */
  bool GAUSS_QUAD = false;
  bool log10F = false;
  double** weights = NULL;

  /*If set to anything besides NULL, this overrides the frequency/PSD arrays
   * passed to sampler. These parameters are for the fishers ONLY*/
  double** fisher_freq = NULL;
  double** fisher_weights = NULL;
  double** fisher_PSD = NULL;
  bool fisher_GAUSS_QUAD = false;
  int* fisher_length = NULL;
  bool fisher_log10F = false;

  // Quadrature for the Fisher when param_space is not set.
  // When param_space is set, GWParameterSpace::quadrature() is used instead.
  Quadrature* QuadMethod = NULL;
  // Always set param_space and likelihood together or not at all.
  GWParameterSpace* param_space = nullptr;
  gw_likelihoods::Likelihood* likelihood = nullptr;

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

struct MCMC_user_param {
  std::complex<double>** burn_data = NULL;
  double** burn_noise = NULL;
  double** burn_freqs = NULL;
  int* burn_lengths = NULL;
  fftw_outline* burn_plans = NULL;
  std::mutex* mFish;

  bool GAUSS_QUAD = false;
  bool log10F = false;
  double** weights = NULL;

  double** fisher_freq = NULL;
  double** fisher_weights = NULL;
  double** fisher_PSD = NULL;
  bool fisher_AD = false;
  bool fisher_GAUSS_QUAD = false;
  int* fisher_length = NULL;
  bool fisher_log10F = false;

  Quadrature* QuadMethod = NULL;

  // RJ stuff
  double** mod_prior_ranges = NULL;
  MCMC_modification_struct* mod_struct;

  gsl_rng* r;
};

struct mcmcVariables {
  double** mcmc_noise = nullptr;
  double* mcmc_init_pos = nullptr;
  std::complex<double>** mcmc_data = nullptr;
  double** mcmc_frequencies = nullptr;
  std::string* mcmc_detectors = nullptr;
  std::string mcmc_generation_method = "";
  int* mcmc_data_length = nullptr;
  fftw_outline* mcmc_fftw_plans = nullptr;
  int mcmc_num_detectors;
  double mcmc_gps_time = 0;
  double mcmc_gmst = 0;
  int mcmc_max_dim;
  int mcmc_min_dim;
  MCMC_modification_struct* mcmc_mod_struct = nullptr;
  bool mcmc_save_waveform = true;
  bool mcmc_intrinsic = false;
  int mcmc_deriv_order = 4;
  bool mcmc_log_beta = false;
  MCMC_user_param* user_parameters = nullptr;
  double maxDim;
  GWParameterSpace* param_space = nullptr;
  gw_likelihoods::Likelihood* likelihood = nullptr;
  Quadrature* QuadMethod = nullptr;
};
struct mcmcVariablesRJ {
  double** mcmc_noise = nullptr;
  double* mcmc_init_pos = nullptr;
  std::complex<double>** mcmc_data = nullptr;
  double** mcmc_frequencies = nullptr;
  std::string* mcmc_detectors = nullptr;
  std::string mcmc_generation_method = "";
  std::string mcmc_generation_method_extended = "";
  int* mcmc_data_length = nullptr;
  fftw_outline* mcmc_fftw_plans = nullptr;
  int mcmc_num_detectors;
  double mcmc_gps_time = 0;
  double mcmc_gmst = 0;
  int mcmc_max_dim;
  int mcmc_min_dim;
  MCMC_modification_struct* mcmc_mod_struct = nullptr;
  bool mcmc_save_waveform = true;
  bool mcmc_intrinsic = false;
  int mcmc_deriv_order = 4;
  bool mcmc_log_beta = false;
  MCMC_user_param* user_parameters = nullptr;
  int maxDim;
  int minDim;
};

void PTMCMC_method_specific_prep(std::string generation_method, int dimension,
                                 bool* intrinsic,
                                 MCMC_modification_struct* mod_struct);
void RJPTMCMC_method_specific_prep(std::string generation_method, int dimension,
                                   bool* intrinsic,
                                   MCMC_modification_struct* mod_struct);

std::string MCMC_prep_params(double* param, double* temp_params,
                             gen_params_base<double>* gen_params, int dimension,
                             std::string generation_method,
                             MCMC_modification_struct* mod_struct,
                             bool intrinsic, double gmst);

bayesship::bayesshipSampler* PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW(
    int dimension, int independentSamples, int ensembleSize, int ensembleN,
    bayesship::positionInfo* initialPosition,
    bayesship::positionInfo** initialEnsemble, double swapProb,
    int burnIterations, int burnPriorIterations, int priorIterations,
    bool writePriorData, int batchSize, double** priorRanges,
    bayesship::probabilityFn* lp, int numThreads, bool pool, int num_detectors,
    std::complex<double>** data, double** noise_psd, double** frequencies,
    int* data_length, double gps_time, std::string* detectors,
    MCMC_modification_struct* mod_struct, std::string generation_method,
    std::string outputDir, std::string outputFileMoniker,
    bool ignoreExistingCheckpoint, bool restrictSwapTemperatures,
    bool coldChainStorageOnly);

bayesship::bayesshipSampler* RJPTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW(
    int minDim, int maxDim, int independentSamples, int ensembleSize,
    int ensembleN, bayesship::positionInfo* initialPosition,
    bayesship::positionInfo** initialEnsemble, double swapProb,
    int burnIterations, int burnPriorIterations, int priorIterations,
    bool writePriorData, int batchSize, double** priorRanges,
    bayesship::probabilityFn* log_prior, int numThreads, bool pool,
    int num_detectors, std::complex<double>** data, double** noise_psd,
    double** frequencies, int* data_length, double gps_time,
    std::string* detectors, MCMC_modification_struct* mod_struct,
    std::string generation_method, std::string generation_method_extended,
    std::string outputDir, std::string outputFileMoniker,
    bool ignoreExistingCheckpoint, bool restrictSwapTemperatures,
    bool coldChainStorageOnly);

double maximized_Log_Likelihood_aligned_spin_internal(
    std::complex<double>* data, double* psd, double* frequencies,
    std::complex<double>* detector_response, size_t length, fftw_outline* plan);

double maximized_Log_Likelihood_unaligned_spin_internal(
    std::complex<double>* data, double* psd, double* frequencies,
    std::complex<double>* hplus, std::complex<double>* hcross, size_t length,
    fftw_outline* plan);

double Log_Likelihood_internal(std::complex<double>* data, double* psd,
                               double* frequencies, double* weights,
                               std::complex<double>* detector_response,
                               int length, bool log10F,
                               std::string integration_method);

double MCMC_likelihood_extrinsic(bool save_waveform,
                                 gen_params_base<double>* parameters,
                                 std::string generation_method,
                                 int* data_length, double** frequencies,
                                 std::complex<double>** data, double** psd,
                                 double** weights,
                                 std::string integration_method, bool log10F,
                                 std::string* detectors, int num_detectors);

/// @brief Runs a Metropolis-Hastings chain from @p initial_params to locate
/// good fiducial and test waveforms for relative-binning initialization.
///
/// Detector data (frequencies, PSDs, strain) is obtained from
/// @p likelihood. The parameter-space dimension, generation
/// method, and Fisher quadrature are obtained from @p param_space.
///
/// \param param_space    Parameter-space object (dimension, to_gen_params,
/// etc.)
/// \param likelihood     Likelihood whose get_ifos() supplies detector data.
/// \param initial_params Starting MCMC parameter vector (length
/// param_space.dim()).
/// \param log_prior      Log-prior probability function for the model.
/// \param param_scales   Typical scale for each parameter (e.g. prior width).
/// \param num_mh_steps   Number of M-H steps.
/// \param fiducial_out  [out] Detector responses at the MAP point.
/// \param test_out      [out] Detector responses at the final M-H step.
/// \param test_gp_out   [out] gen_params at the final M-H step (optional).
void find_fiducial(const GWParameterSpace& param_space,
                   const gw_likelihoods::Likelihood& likelihood,
                   const double* initial_params,
                   bayesship::probabilityFn* log_prior,
                   const std::vector<double>& param_scales, int num_mh_steps,
                   std::vector<VECCPL>& fiducial_out,
                   std::vector<VECCPL>& test_out,
                   gen_params_base<double>* test_gp_out = nullptr);

#endif
