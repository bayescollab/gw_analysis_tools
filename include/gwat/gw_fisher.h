#ifndef GW_FISHER_H
#define GW_FISHER_H

/// @file
/// @brief ParameterMap-aware Fisher matrix computation.

#include <memory>
#include <vector>

#include "likelihoods.h"
#include "parameter_map.h"
#include "quadrature.h"
#include "waveform_generator_v2.h"
#include "waveform_util.h"

namespace gw_fishers {
class GWNumericalFishers {
 public:
  virtual ~GWNumericalFishers() = default;

  /// @brief Compute the Fisher matrix
  /// @param[out] fisherM  fisherM[i][j] = (∂h/∂θ_i | ∂h/∂θ_j).
  /// @param theta         MCMC parameter vector, length dim().
  virtual void compute_Fisher(double** fisherM, const double* theta) const = 0;

  /// @brief Generate polarization modes at the given frequencies for @p theta.
  /// Used by find_fiducial to extract fiducial and test waveforms.
  virtual std::vector<VECCPL> generate_modes(const double* theta,
                                              const VECDBL& freqs) const = 0;
};

class PolarizationDataFishers : public GWNumericalFishers {
 private:
  const ParameterMap& pmap_;
  const Quadrature& quad_;
  const std::vector<PolarizationData>& ifos_;
  const waveform_generator::WaveformGenerator& wf_gen_;
  const std::vector<double> sky_avg_factors_;
  int order_;
  int eta_idx_;
  static constexpr double eps_ = 1e-8;

  const double f_ref_, gmst_;
  const bool shift_time_, shift_phase_;

  std::vector<std::vector<VECCPL>> compute_derivative_in_detector(
      const double* theta, const PolarizationData& ifo) const;

 public:
  PolarizationDataFishers(const ParameterMap& pmap,
                          const std::vector<PolarizationData>& ifos,
                          const waveform_generator::WaveformGenerator& wf_gen,
                          std::vector<double> sky_avg_factors,
                          const Quadrature& quad, double f_ref, double gmst,
                          bool shift_time, bool shift_phase, int order)
      : pmap_(pmap),
        ifos_(ifos),
        wf_gen_(wf_gen),
        sky_avg_factors_(sky_avg_factors),
        quad_(quad),
        f_ref_(f_ref),
        gmst_(gmst),
        shift_time_(shift_time),
        shift_phase_(shift_phase),
        order_(order),
        eta_idx_(pmap.index_of("eta")) {};

  void compute_Fisher(double** fisherM, const double* theta) const override;
  std::vector<VECCPL> generate_modes(const double* theta,
                                     const VECDBL& freqs) const override;
};

/// @brief RB Fisher for sky-averaged polarization data.
///
/// Finite-difference derivatives are evaluated only at the bin-edge frequencies
/// from an already-initialised RelativeBinningBisectionPolarizationsLikelihood.
/// The derivative waveforms are normalised against the stored fiducial via
/// compute_waveform_ratios, and the Fisher elements are accumulated using the
/// B0/B1 template-template summary data — so no full-grid waveform evaluations
/// are needed after construction.
class RelativeBinningPolarizationsFishers : public GWNumericalFishers {
 public:
  using RBLikelihood = gw_likelihoods::RelativeBinning::
      RelativeBinningBisectionPolarizationsLikelihood;

  RelativeBinningPolarizationsFishers(
      const ParameterMap& pmap, const RBLikelihood& rb_likelihood,
      const waveform_generator::WaveformGenerator& wf_gen,
      std::vector<double> sky_avg_factors, double f_ref, double gmst,
      bool shift_time, bool shift_phase, int order = 4)
      : pmap_(pmap),
        rb_(rb_likelihood),
        wf_gen_(wf_gen),
        sky_avg_factors_(std::move(sky_avg_factors)),
        f_ref_(f_ref),
        gmst_(gmst),
        shift_time_(shift_time),
        shift_phase_(shift_phase),
        order_(order),
        eta_idx_(pmap.index_of("eta")) {}

  void compute_Fisher(double** fisherM, const double* theta) const override;
  std::vector<VECCPL> generate_modes(const double* theta,
                                     const VECDBL& freqs) const override;

 private:
  const ParameterMap& pmap_;
  const RBLikelihood& rb_;
  const waveform_generator::WaveformGenerator& wf_gen_;
  const std::vector<double> sky_avg_factors_;
  const double f_ref_, gmst_;
  const bool shift_time_, shift_phase_;
  int order_;
  int eta_idx_;
  static constexpr double eps_ = 1e-8;
};

/// @brief Add prior-based diagonal regularization to a Fisher matrix in-place.
///
/// For each sampled parameter i with prior width w_i = hi_i - lo_i (taken from
/// pmap.specs()), adds 1/w_i^2 to fisherM[i][i].  This is the diagonal of the
/// prior Hessian under a uniform prior on [lo_i, hi_i], making fisherM the
/// Hessian of -log(posterior) at the MAP point.  Parameters with zero width are
/// skipped.
///
/// @param fisherM  Square dim×dim Fisher matrix (modified in place).
/// @param pmap     ParameterMap; supplies specs(), index_of().
void apply_prior_regularization(double** fisherM, const ParameterMap& pmap);

}  // namespace gw_fishers

#endif  // GW_FISHER_H
