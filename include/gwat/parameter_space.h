#ifndef PARAMETER_SPACE_H
#define PARAMETER_SPACE_H

#include <memory>
#include <string>
#include <vector>

#include "gw_fisher.h"
#include "gw_prior.h"
#include "quadrature.h"
#include "util.h"

/// @brief Uniform box prior built from per-parameter [lo, hi] bounds.
///
/// Holds the bounds of each parameter in order.
class BoxPrior {
  const std::vector<DBLPAIR> bounds_;

 public:
  explicit BoxPrior(const std::vector<DBLPAIR> b) : bounds_(std::move(b)) {}

  /// @brief Compute the log-prior for a set of parameters.
  ///
  /// @param params  The parameter vector. Must match the order of @p bounds_ .
  /// @return  0 if all prior bounds are satisfied, otherwise -inf.
  double ln_prior(const double* params) const {
    for (int i = 0; i < static_cast<int>(bounds_.size()); ++i) {
      if (params[i] < bounds_[i].first || params[i] > bounds_[i].second)
        return bayesship::limitInf;
    }
    return 0.0;
  }
};

/// @brief Abstract base class defining parameter layout and Fisher support for
/// a user-defined gravitational-wave model.
///
/// Concrete subclasses define:
///   - the MCMC parameter vector layout via param_names, to_gen_params,
///   from_gen_params
///   - the generation method name for Fisher computation
///   - the prior function
///   - an optional Fisher regularization prior via fisher_prior
///
/// When a GWParameterSpace* is stored in mcmcVariables::param_space, the MCMC
/// machinery bypasses the legacy "MCMC_" generation-method decision trees for
/// the numerical Fisher matrix and parameter proposals. The likelihood is
/// evaluated separately via mcmcVariables::likelihood if set.
///
class GWParameterSpace {
 protected:
  std::vector<std::string> param_names;
  int extrinsic_params_n = 0;

 public:
  int dim() const { return static_cast<int>(param_names.size()); }
  const std::vector<std::string>& names() const { return param_names; }

  int index_of(const std::string& name) const {
    for (int i = 0; i < (int)param_names.size(); ++i)
      if (param_names[i] == name) return i;
    return -1;
  }

  int number_of_extrinsic_params() const { return extrinsic_params_n; }

  /// @brief MCMC vector → gen_params_base (replaces repack_parameters).
  ///
  /// Must populate all fields needed by the waveform generator, including
  /// non-varying quantities such as gmst, f_ref, sky_average, and any
  /// modification arrays (betappe, etc.).
  /// @param vec  Parameter vector.
  /// @param p    [out] gen_params_base.
  virtual void to_gen_params(const double* vec,
                             gen_params_base<double>& p) const = 0;

  /// @brief gen_params_base → MCMC vector (replaces unpack_parameters).
  /// @param p    gen_params_base object
  /// @param vec  [out] MCMC vector
  virtual void from_gen_params(const gen_params_base<double>& p,
                               double* vec) const = 0;

  /// Returns a pointer to the prior built by the subclass, or nullptr.
  virtual bayesship::probabilityFn* prior_fn() const { return nullptr; }

  // Optionally add a regularizing prior contribution to the Fisher diagonal.
  virtual void fisher_prior(double** fisher) const {}

  // Waveform method name.
  virtual std::string generation_method() const = 0;

  // Quadrature method to use for numerical Fisher computation.
  // Must be consistent with the integration used inside log_likelihood.
  virtual const Quadrature* quadrature() const = 0;

  // Factory: returns a heap-allocated derivative strategy for this parameter
  // space. The caller owns the returned object. The default implementation
  // (defined in gw_fisher.cpp) creates a GWParameterSingleDetectorDerivatives
  // using central differences with the given order.
  virtual std::unique_ptr<GWFisherDerivatives> fisher_derivatives(
      int order = 4) const;

  virtual ~GWParameterSpace() = default;
};

#endif  // PARAMETER_SPACE_H
