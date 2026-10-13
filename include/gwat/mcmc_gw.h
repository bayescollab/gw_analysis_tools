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

#include "gw_fisher.h"
#include "likelihoods.h"
#include "mcmc_sampler.h"
#include "parameter_map.h"
#include "quadrature.h"
#include "util.h"

/// @brief Bundle of new-API objects that together define a GW inference model.
///
/// When model != nullptr in mcmcVariables or MCMC_modification_struct, the MCMC
/// uses this path rather than the legacy generation-method decision trees.
struct GWModel {
  ParameterMap*                     param_map  = nullptr;
  bayesship::probabilityFn*         prior      = nullptr;
  gw_fishers::GWNumericalFishers*   fisher     = nullptr;
  gw_likelihoods::Likelihood*       likelihood = nullptr;

  /// Optional: used instead of `likelihood` for chains at beta==0 (the
  /// hottest PT rung, which by design explores the full prior and ignores
  /// the likelihood for acceptance -- see MCMC_likelihood_wrapper::eval()).
  /// Meant for when `likelihood` is a fast approximation (e.g. relative
  /// binning) with a validity guard that would reject almost the entire
  /// prior volume there, corrupting evidence's thermodynamic integration at
  /// that anchor point: only ensembleN chains ever sit at beta==0, so
  /// paying for an exact likelihood just there is cheap. Leave null to use
  /// `likelihood` at every beta (default, backward compatible).
  gw_likelihoods::Likelihood*       beta_zero_likelihood = nullptr;

  bool active() const { return param_map != nullptr && likelihood != nullptr; }
};

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

  Quadrature* QuadMethod = NULL;
  // Set model to use the new ParameterMap / Likelihood / Fisher path.
  GWModel* model = nullptr;

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
  GWModel* model = nullptr;
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

/// Overload when using GWModel in mod_struct
bayesship::bayesshipSampler* PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW(
    bayesship::positionInfo* initialPosition, int independentSamples,
    int ensembleSize, int ensembleN, double swapProb, int burnIterations,
    int burnPriorIterations, int priorIterations, bool writePriorData,
    int batchSize, int numThreads, bool pool, double gps_time,
    MCMC_modification_struct* mod_struct, std::string outputDir,
    std::string outputFileMoniker, bool ignoreExistingCheckpoint,
    bool restrictSwapTemperatures, bool coldChainStorageOnly) {
  if (mod_struct->model->active()) {
    const GWModel* model = mod_struct->model;
    return PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW(
        model->param_map->dim(), independentSamples, ensembleSize, ensembleN,
        initialPosition, (bayesship::positionInfo**)nullptr, swapProb,
        burnIterations, burnPriorIterations, priorIterations, writePriorData,
        batchSize, (double**)nullptr, model->prior, numThreads, pool,
        model->likelihood->detector_number(), (std::complex<double>**)nullptr,
        (double**)nullptr, (double**)nullptr, (int*)nullptr, gps_time,
        (std::string*)nullptr, mod_struct, "", outputDir, outputFileMoniker,
        ignoreExistingCheckpoint, restrictSwapTemperatures,
        coldChainStorageOnly);
  } else {
    throw std::runtime_error(
        "This is an overload of PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW "
        "meant for GWModel, but mod_struct->model is not active.");
  }
  return nullptr;
}

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
/// good fiducial and test parameter points for relative-binning initialization.
///
/// Returns the MAP parameter vector and the final M-H step parameter vector.
/// The caller is responsible for generating the corresponding waveforms.
///
/// \param model             GWModel supplying param_map, likelihood, and
/// fisher.
/// \param initial_params    Starting parameter vector (length
/// param_map->dim()).
/// \param log_prior         Log-prior probability function for the model.
/// \param num_mh_steps      Number of M-H steps.
/// \param map_params_out   [out] Parameter vector at the MAP point.
/// \param final_params_out [out] Parameter vector at the final M-H step.
///
/// M-H proposal widths are derived from the prior bounds in model.param_map:
/// sigma_i = min(Fisher, 0.5 * w_i) where w_i = hi_i - lo_i.
void find_fiducial(const GWModel& model, const VECDBL& initial_params,
                   bayesship::probabilityFn* log_prior, int num_mh_steps,
                   std::vector<double>& map_params_out,
                   std::vector<double>& final_params_out);

/// @brief Calibrates a validity_guard_tol for a
/// RelativeBinningBisectionPolarizationsLikelihood, specific to the given
/// injection/Fisher matrix, instead of hardcoding a cross-injection constant.
///
/// Runs a staged per-parameter (Gibbs) M-H chain around @p fiducial_params --
/// the same proposal scheme find_fiducial uses (sigma_i =
/// min(c_mh/sqrt(Gamma_ii), 0.5*prior_width_i), c_mh = 2.38/sqrt(dim), one
/// parameter proposed per step, Metropolis acceptance) -- in blocks of @p
/// stage_steps, each stage computing its own (independent,
/// post-burn-in-for-stage-1-only) max of rb_ll's guard residual (the same
/// diagnostic last_guard_residual() reports). Plateau is judged by comparing
/// consecutive stages' own maxes, not a cumulative running max: a cumulative
/// max can only grow or stay flat, so a single outlier step in one stage would
/// permanently floor every later stage's reading there, making the very next
/// stage trivially look plateaued regardless of whether the chain actually
/// stabilized. The function stops once consecutive stage maxes agree within @p
/// plateau_rel_tol (relative, either direction), or after @p max_stages stages,
/// whichever comes first; the returned tol is based on that plateaued per-stage
/// value, not the all-time running max across all stages. An unguarded
/// calibration chain can itself wander into a genuine phantom-mode excursion (a
/// huge, one-off residual spike) -- that is precisely the kind of point the
/// guard exists to reject, so it must stay well below tol rather than inflate
/// tol past itself. The all-time running max is tracked only to emit a stdout
/// NOTE when it significantly exceeds the plateaued value, flagging that such
/// an excursion occurred.
///
/// If a stage has zero accepted proposals, the chain has frozen: it landed
/// where every neighboring step evaluates to a non-finite RB log-likelihood,
/// so nothing can ever be accepted again. That frozen residual is real but
/// not a representative tail value, so it is discarded -- the function logs
/// a warning and returns using the last stage that still had
/// healthy mixing (or the fiducial's own residual, if even the first stage
/// froze) instead. Seeing this warning means the fiducial/test-point pair
/// used to build rb_ll's bins is not robust and the returned tol should not
/// be trusted as a tight threshold.
///
/// A plain joint (all-parameters-at-once) Gaussian proposal from the full
/// Fisher covariance does not work for this: at realistic SNR the Fisher
/// matrix is typically too ill-conditioned (the tightest- and
/// loosest-constrained directions can differ by many orders of magnitude)
/// for a joint step to mix, so it predominantly proposes points with
/// catastrophic likelihood -- the same failure find_fiducial itself avoids
/// by proposing one parameter at a time, reused here for the same reason.
///
/// rb_ll's own validity_guard_tol() is saved and restored around the
/// calibration run (temporarily set large enough to never reject, so every
/// step's residual is recorded); call rb_ll.set_validity_guard_tol() yourself
/// with the returned value once you're satisfied with it.
///
/// \param rb_ll            RB likelihood to calibrate against.
/// \param model            Supplies param_map, prior, and fisher (the same
/// GWModel bundle find_fiducial takes).
/// \param fiducial_params  Chain starting point (typically find_fiducial's
/// map_params_out).
/// \param stage_steps      Gibbs steps per stage.
/// \param max_stages       Upper bound on stages run.
/// \param plateau_rel_tol  Stop once the running max residual grows by less
/// than this fraction between consecutive stages.
/// \param safety_margin    Returned tol = safety_margin * plateaued max
/// residual.
///
/// \return Recommended validity_guard_tol.
double calibrate_validity_guard_tol(
    gw_likelihoods::RelativeBinning::RelativeBinningBisectionPolarizationsLikelihood& rb_ll,
    const GWModel& model, const VECDBL& fiducial_params,
    int stage_steps = 20000, int max_stages = 15,
    double plateau_rel_tol = 0.05, double safety_margin = 10.0);

/// @brief Generates deterministic candidate points along each parameter's
/// own Fisher axis at @p fiducial_params, for use as additional RB
/// bin-construction test points (see
/// RelativeBinningBisectionPolarizationsLikelihood's multi-test-point
/// constructor) instead of relying on a single find_fiducial draw.
///
/// find_fiducial's final_params_out is just the chain's raw position after a
/// fixed number of steps -- one draw whose excursion size from the fiducial
/// varies by an order of magnitude or more between otherwise-identical runs
/// (seen directly: the same system's RB bin count swung from 5 to 21 purely
/// from varying fiducial_mh_steps). Since RelativeBinningBisection*'s bins are
/// built from comparing the fiducial against whichever single point it is
/// handed, a single unlucky draw can leave the bins too coarse. This
/// function instead returns fiducial_params +/- n_sigma/sqrt(Gamma_ii) along
/// every parameter's own axis, deterministic and reproducible. Gamma is
/// regularized via apply_prior_regularization() first (Gamma_ii += 1/w_i^2,
/// w_i the prior width) so a weakly-constrained direction gets a sensibly
/// bounded sigma_i (<= w_i) instead of an unbounded one whose candidate
/// would otherwise just be silently dropped by the prior check below --
/// unlike find_fiducial/calibrate_validity_guard_tol, this function has no
/// separate min(sigma, 0.5*prior_width) clip to fall back on.
///
/// An earlier version of this function used the Fisher matrix's full
/// eigen-decomposition instead of per-parameter axes, to also cover
/// correlated/degenerate directions (e.g. a phiRef-tc degeneracy) that
/// per-parameter exploration can miss. That does not work for a Fisher
/// matrix this ill-conditioned (condition number ~1e13, from
/// Gamma_lnMc~2.3e16 vs Gamma_tc~900): its smallest eigenvalues are
/// dominated by numerical noise in the finite-difference off-diagonal
/// terms, not real degeneracies, giving some eigenvectors an astronomically
/// large (physically meaningless) implied sigma -- whitening by each
/// parameter's own diagonal scale first did not fix this. Per-parameter
/// axes reuse the exact diagonal-only formula find_fiducial and
/// calibrate_validity_guard_tol already rely on, which is numerically sane
/// for this same matrix, at the cost of not covering correlated directions.
///
/// A candidate is dropped (not clipped) if @p model's prior rejects it:
/// BHBPriorFn's non-rectangular mass/spin constraints mean an axis-aligned
/// Fisher step can land outside the physically allowed region even when
/// each parameter individually stays in its own prior box, and RB does not
/// need to be resolved somewhere the prior already forbids the sampler from
/// going.
///
/// \param model            Supplies param_map, prior, and fisher (the same
/// GWModel bundle find_fiducial takes).
/// \param fiducial_params  Point the Fisher matrix is evaluated at.
/// \param n_sigma          Distance along each parameter's own axis, in
/// that parameter's own stdev (1/sqrt(Gamma_ii)).
/// \return One VECDBL per kept candidate (up to 2*dim; fewer if some are
/// rejected by the prior).
std::vector<VECDBL> fisher_axis_points(const GWModel& model,
                                       const VECDBL& fiducial_params,
                                       double n_sigma = 5.0);

#endif
