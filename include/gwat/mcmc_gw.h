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
#include <string>
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

  // ---- Checkpoint / resume settings ----
  // Passed to the bayesshipSampler members of the same name; the defaults
  // match BayesShip's. See docs/checkpointing.md in BayesShip.
  // Recent steps each hot chain keeps in memory and in every checkpoint
  // (coldOnlyStorage, non-RJ only). 0 keeps the full history and saves none.
  int hotChainHistoryLength = 2000;
  // Complete checkpoint generations kept on disk.
  int checkpointGenerationsToKeep = 2;
  // On resume, the re-evaluated log likelihood/prior of every chain must match
  // the checkpoint to this relative tolerance. Negative disables the check.
  double resumeValueRelTol = 1e-10;
  // Stop cleanly after a batch once this many seconds since sampling began
  // leave no room for another. 0: no limit of its own (SLURM_JOB_END_TIME and
  // BAYESSHIP_WALL_TIME_LIMIT still apply).
  double wallTimeLimitSeconds = 0;
  // Margin kept free at the end of the wall time (BAYESSHIP_WALL_TIME_MARGIN
  // overrides it).
  double wallTimeMarginSeconds = 120;
  // Catch SIGTERM/SIGUSR1/SIGINT and stop after the current batch.
  bool handleStopSignals = true;
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

/// Settings for @ref calibrate_validity_guard_tol(). The defaults suit an
/// expensive waveform (tens of ms per call): typically a few 10^4 steps.
struct ValidityGuardCalibrationOptions {
  /// Number of chains. They share the DE history and are stepped in rounds
  /// of check_interval steps, in parallel on n_threads threads.
  int n_chains = 1;
  /// Threads the chains run on; 0 means one per chain. The result does not
  /// depend on it.
  int n_threads = 0;
  /// Hard cap on each chain's steps (burn-in included).
  int max_steps = 200000;
  /// Each chain's burn-in bounds. 0 for min_burn_steps means 800 * dim.
  int min_burn_steps = 0;
  int max_burn_steps = 50000;
  /// Post-burn-in steps, summed over the chains, required before
  /// convergence is tested.
  int min_samples = 20000;
  /// Steps between convergence tests (and progress lines).
  int check_interval = 2000;
  /// Quantile of the post-burn-in residuals the tolerance is based on.
  double quantile = 0.999;
  /// Converged once the quantiles of two halves of the post-burn-in
  /// residuals agree to this relative tolerance: the first and second half
  /// of the steps for one chain, the first and second half of the chains
  /// for several.
  double rel_tol = 0.1;
  /// Returned tol = margin * quantile.
  double margin = 5.0;
  /// Probability that a step moves one parameter; otherwise it moves along
  /// one eigen-direction of the Fisher matrix.
  double coordinate_fraction = 0.5;
  /// Probability that a step is a differential-evolution move
  /// x + gamma * (h_a - h_b) from two past positions (once the history holds
  /// 2 * dim of them); otherwise a 1-D move as below.
  double de_fraction = 0.5;
  /// Every de_thin-th position is added to the DE history.
  int de_thin = 10;
  /// Halfway through burn-in, replace the Fisher eigen-directions with those
  /// of the covariance of the burn-in positions so far.
  bool empirical_directions = true;
  /// Acceptance rate each direction's step width is adapted toward during
  /// burn-in (0.44 is optimal for a 1-D Gaussian random walk).
  double target_acceptance = 0.44;
  /// gsl RNG seed of chain 0 (0 is gsl's default seed); chain c uses
  /// seed + c * 2654435761.
  unsigned long seed = 0;
  /// If non-empty, every step of every chain is written to this CSV file.
  std::string trace_path;
};

/// What calibrate_validity_guard_tol() found.
struct ValidityGuardCalibration {
  double tol = 0;                // margin * quantile_residual
  double quantile_residual = 0;  // over all post-burn-in steps of all chains
  /// The two halves of the convergence test (steps for one chain, chains for
  /// several), at the last test.
  double quantile_first_half = 0, quantile_second_half = 0;
  double max_residual = 0;       // largest post-burn-in residual
  double fiducial_residual = 0;
  int n_chains = 1;
  int steps = 0;                 // steps per chain, burn-in included
  int burn_steps = 0;            // longest burn-in of any chain
  int samples = 0;               // post-burn-in steps, all chains
  double acceptance = 0;         // post-burn-in acceptance rate
  bool converged = false;        // false: stopped at max_steps or froze
  std::string stop_reason;
};

/// @brief Calibrates a validity_guard_tol for a
/// RelativeBinningBisectionPolarizationsLikelihood, specific to the given
/// injection/Fisher matrix, instead of hardcoding a cross-injection constant.
///
/// Runs one or more Metropolis chains on the RB posterior (RB log-likelihood
/// + the model's prior) from @p fiducial_params, records rb_ll's guard
/// residual at every step after burn-in, and returns
/// margin * (a high quantile of those residuals).
///
/// **Proposal.** Three kinds of Metropolis move:
/// - differential evolution (DE), with probability de_fraction once the
///   history holds 2*dim positions: x + gamma * (h_a - h_b) for two past
///   positions (every de_thin-th position is kept), gamma = 2.38/sqrt(2 dim)
///   (1 for a tenth of the moves), jittered by 10%. Its jumps take the shape
///   and size of the region explored so far, so they open up long, curved
///   degeneracies that fixed directions crawl along. Like DREAM(ZS), the
///   proposal depends on the chain's past.
/// - otherwise a 1-D Gaussian move along one direction, like BayesShip's
///   Fisher proposal: with probability coordinate_fraction along one
///   parameter's axis, otherwise along one eigen-direction.
///
/// The eigen-directions start as those of the Fisher matrix at the fiducial,
/// regularized with the prior widths (Gamma_ii += 1/w_i^2) and decomposed
/// after scaling it to unit diagonal (D^-1/2 Gamma D^-1/2): the raw matrix's
/// condition number (~1e13 at realistic SNR) is mostly units, and decomposing
/// it directly leaves the small eigenvalues as numerical noise. With
/// empirical_directions, halfway through burn-in they are replaced by the
/// eigenvectors of the covariance of the burn-in positions so far, which
/// follow the posterior where the finite-difference Fisher's do not. A joint
/// all-parameters Gaussian move from the Fisher covariance is not used: with
/// an ill-conditioned, finite-difference Fisher it mostly proposes
/// catastrophic points.
///
/// Each direction starts at its natural width (1/sqrt(Gamma_ii) along a
/// parameter axis; 1/sqrt(eigenvalue), or 2.4 sample standard deviations,
/// along an eigen-direction), never more than half the prior box, and its
/// width is adapted toward target_acceptance during burn-in, then frozen. A
/// direction whose starting width is wrong is shrunk by the adaptation
/// instead of wasting steps.
///
/// Tested on an IMRPhenomD injection (dim 8): the spin posterior is far
/// wider than the Fisher at the fiducial suggests (chi1 std ~0.1), and
/// one-parameter and Fisher eigen-direction moves alone needed ~10^4-10^5
/// steps per independent sample. With DE and empirical directions, 600k-step
/// chains agree on the 0.999 quantile to 5% (10% without DE), and
/// default-length runs land 0-40% below it.
///
/// **Burn-in.** The chain starts at the fiducial, near the posterior mode,
/// where residuals are smallest. Near the mode, logL_max - logL is
/// chi^2_dim/2, typically dim/2. Burn-in ends at twice the step at which the
/// log posterior first falls dim/2 below the largest value seen, bounded by
/// [min_burn_steps, max_burn_steps].
///
/// **Estimate and stopping.** The tolerance comes from a quantile of all
/// post-burn-in residuals, not a maximum: one outlier barely moves it, and an
/// unguarded calibration chain can itself wander into a phantom-mode
/// excursion (a huge one-off residual), which the guard exists to reject and
/// must not set the tolerance. Convergence: once min_samples post-burn-in
/// steps exist, every check_interval steps the quantiles of two halves are
/// compared; the calibration stops when they agree within rel_tol, or at
/// max_steps (result marked not converged).
///
/// **Several chains** (n_chains > 1). The residuals' upper tail comes from
/// rare excursions that one chain revisits only every ~10^4 steps (EFPE): a
/// single chain can spend its first few 10^4 steps where residuals are small
/// and pass its own split-half test there. Independent chains explore the
/// tail independently, and run in parallel, so the same number of samples
/// takes 1/n_chains of the time (burn-in is still paid by every chain).
/// - Chain 0 starts at the fiducial. Chain c > 0 starts at a draw from the
///   Gaussian of the Fisher at the fiducial, accepted only if its log
///   posterior is within dim of the fiducial's (each 100 failed draws halve
///   the draw's scale; after 1000 the chain starts at the fiducial).
/// - Every chain has its own RNG, proposal widths, burn-in and empirical
///   directions. The DE history is shared: a DE move draws from every
///   chain's history up to the start of the current round, plus its own
///   since, so jumps cover the region all chains have explored.
/// - The chains step in rounds of check_interval steps. After each round,
///   the new history, trace rows and messages are merged in chain order,
///   so the result depends on the seed and n_chains, not on n_threads or
///   thread timing.
/// - Convergence needs every chain past burn-in, and compares the quantile
///   of the first n_chains/2 chains' residuals with that of the rest,
///   instead of the first and second halves of the steps: two groups of
///   independent chains are unlikely to both miss the tail.
/// With n_chains = 1 the calibration is the single chain described above,
/// step for step.
///
/// If a whole check interval has zero accepted proposals in any chain, that
/// chain has frozen where every neighboring step has a non-finite RB
/// log-likelihood (a genuine RB breakdown region). The calibration stops with
/// a warning, and the result is marked not converged: the fiducial/test
/// points behind rb_ll's bins may not be robust.
///
/// rb_ll is called from several threads at once with n_threads > 1, through
/// its thread-safe log_likelihood(params, guard) overload; the waveform
/// generator and the model's prior must be thread-safe, as for BayesShip's
/// threaded sampling.
///
/// rb_ll's own validity_guard_tol() is saved and restored around the
/// calibration (temporarily set large enough to never reject, so every
/// step's residual is recorded); call rb_ll.set_validity_guard_tol() yourself
/// with the returned tol.
///
/// \param rb_ll            RB likelihood to calibrate against.
/// \param model            Supplies param_map, prior, and fisher (the same
/// GWModel bundle find_fiducial takes).
/// \param fiducial_params  Chain starting point (typically find_fiducial's
/// map_params_out).
/// \param opts             See ValidityGuardCalibrationOptions.
ValidityGuardCalibration calibrate_validity_guard_tol(
    gw_likelihoods::RelativeBinning::RelativeBinningBisectionPolarizationsLikelihood& rb_ll,
    const GWModel& model, const VECDBL& fiducial_params,
    const ValidityGuardCalibrationOptions& opts = {});

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
/// unlike find_fiducial, this function has no separate
/// min(sigma, 0.5*prior_width) clip to fall back on.
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
/// axes reuse the exact diagonal-only formula find_fiducial relies on,
/// which is numerically sane for this same matrix, at the cost of not
/// covering correlated directions. (calibrate_validity_guard_tol does use
/// eigen-directions, but only as proposal directions whose widths it adapts,
/// so a noise-dominated eigenvalue costs it nothing.)
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
