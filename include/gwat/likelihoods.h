#ifndef LIKELIHOODS_H
#define LIKELIHOODS_H

#include <stdexcept>
#include <vector>

#include "waveform_generator_v2.h"
#include "waveform_util.h"

namespace gw_likelihoods {

using waveform_generator::WaveformGenerator;

///@class Likelihood
///@brief Base class for likelihoods
class Likelihood {
 public:
  Likelihood() = default;
  virtual ~Likelihood() = default;

  /// @brief Compute the log likelihood -1/2[(h|h) - 2(d|h)] over all detectors
  /// @param params             Model parameters.
  virtual double log_likelihood(gen_params_base<double>* params) const = 0;

  /// @brief The detector strain data used by this likelihood.
  /// Only meaningful for strain-based likelihoods; throws for others.
  virtual const std::vector<IfoData>& get_ifos() const {
    throw std::logic_error("get_ifos() not supported for this likelihood type");
  }

  /// @brief The polarization data used by this likelihood.
  /// Returns nullptr for strain-based likelihoods.
  virtual const PolarizationData* get_polarization_data() const { return nullptr; }
};

/// @class CoherentBareLikelihood
/// @brief The most simple likelihood class, no marginalization. Uses
/// create_coherent_GW_detection_reuse_wf.
class CoherentBareLikelihood : public Likelihood {
 private:
  const Quadrature& quad;
  const std::vector<IfoData>& ifos_;
  const WaveformGenerator& waveform_generator_;
  const double f_ref_, gmst_;
  const bool equatorial_orientation_, horizon_coord_;
  const bool shift_time_, shift_phase_;

 public:
  CoherentBareLikelihood(const std::vector<IfoData>& is, const Quadrature& q,
                         const WaveformGenerator& wf_generator, double f_ref,
                         double gmst, bool equatorial_orientation,
                         bool horizon_coord, bool shift_time, bool shift_phase)
      : ifos_(is),
        quad(q),
        waveform_generator_(wf_generator),
        f_ref_(f_ref),
        gmst_(gmst),
        equatorial_orientation_(equatorial_orientation),
        horizon_coord_(horizon_coord),
        shift_time_(shift_time),
        shift_phase_(shift_phase) {};

  /// @details Compute the unmarginalized log-likelihood.
  double log_likelihood(gen_params_base<double>* params) const override;

  const std::vector<IfoData>& get_ifos() const override { return ifos_; }
};

/// @class PolarizationsLikelihood
/// @brief Sky-averaged likelihood computed directly from waveform polarizations.
///
/// The inner product is summed over all active polarization modes with per-mode
/// sky-average factors:
///   logL = sum_i factor_i * [-1/2 (h_i|h_i) + (d_i|h_i)]
///
/// For GR tensor modes, sky_avg_factors = {1/5, 1/5}. Modes are ordered as in
/// waveform_polarizations: {hplus, hcross, hx, hy, hb, hl}.
class PolarizationsLikelihood : public Likelihood {
 private:
  const PolarizationData& data_;
  const std::vector<double> sky_avg_factors_;
  const Quadrature& quad;
  const WaveformGenerator& waveform_generator_;
  const double f_ref_, gmst_;
  const bool shift_time_, shift_phase_;

 public:
  PolarizationsLikelihood(const PolarizationData& data,
                          std::vector<double> sky_avg_factors,
                          const Quadrature& q,
                          const WaveformGenerator& wf_generator,
                          double f_ref, double gmst,
                          bool shift_time, bool shift_phase)
      : data_(data),
        sky_avg_factors_(std::move(sky_avg_factors)),
        quad(q),
        waveform_generator_(wf_generator),
        f_ref_(f_ref),
        gmst_(gmst),
        shift_time_(shift_time),
        shift_phase_(shift_phase) {
    if (sky_avg_factors_.size() != data_.modes.size())
      throw std::invalid_argument(
          "sky_avg_factors size (" + std::to_string(sky_avg_factors_.size()) +
          ") must match data modes size (" +
          std::to_string(data_.modes.size()) + ")");
  }

  double log_likelihood(gen_params_base<double>* params) const override;

  const PolarizationData* get_polarization_data() const override { return &data_; }
};

namespace RelativeBinning {

// BISECTION RELATIVE BINNING
//
// Bins are placed by recursively bisecting frequency intervals until the
// absolute log-likelihood error within each bin (exact minus approximate,
// evaluated on the test waveform and summed over detectors) falls below
// epsilon.  This makes bin placement data-adaptive rather than PN-ansatz-
// driven.

///@class RelativeBinningBisectionLikelihood
///@brief Relative-binning likelihood binned with a bisection algorithm.
class RelativeBinningBisectionLikelihood : public Likelihood {
 public:
  RelativeBinningBisectionLikelihood(
      const IfoData& ifo, const VECCPL& fiducial_data, const VECCPL& test_data,
      const WaveformGenerator& waveform_generator, double epsilon, double f_ref,
      double gmst, bool equatorial_orientation, bool horizon_coord,
      bool shift_time, bool shift_phase, bool log_spacing = false);

  double log_likelihood(gen_params_base<double>* params) const override;

  VECDBL get_bin_freqs() { return bin_freqs; }
  VECINT get_bin_inds() { return bin_inds; }

  // Evaluate RB logL from waveforms already extracted at bin edges.
  // h_at_bins[det] must have exactly bin_inds.size() entries.
  double log_likelihood_at_waveform(const VECCPL& h_at_bins) const;

  const std::vector<IfoData>& get_ifos() const override { return ifos_vec_; }

 private:
  const WaveformGenerator& waveform_generator_;
  const double f_ref_, gmst_;
  const bool equatorial_orientation_, horizon_coord_;
  const bool shift_time_, shift_phase_;
  /// @struct SummaryData
  /// @brief Holds an interferometer's pre-computed summary data.
  struct SummaryData {
    VECCPL A0, A1;
    VECDBL B0, B1;
  };

  int number_of_bins;
  SummaryData ifo_summary_data_;
  IfoData ifo_;
  // Stable storage for get_ifos(); rebuilt once after setup_summary_data.
  std::vector<IfoData> ifos_vec_;

  /// @brief Indices of the bin edges corresponding to the full-resolution grid,
  /// including the final right edge.
  VECINT bin_inds;
  /// @brief The frequencies at each bin edge.
  VECDBL bin_freqs;
  /// @brief The spacing between bin edges. Size: bin_inds.size()-1
  VECDBL bin_widths;
  /// @brief The center of the bin.
  VECDBL bin_centers;

  std::pair<int, int> find_min_max_indices(const VECCPL& fiducial_data);

  /// @brief Returns the absolute logL error (exact minus approximate, summed
  /// over detectors) for the proposed bin [left_idx, right_idx).
  double bin_log_likelihood_error(int left_idx, int right_idx,
                                  const IfoData& ifo,
                                  const VECCPL& fiducial_data,
                                  const VECCPL& test_data, bool log_spacing);

  /// @brief Bisection algorithm. Populates the bin_* vectors.
  void bin_bisection(const IfoData& ifo, const VECCPL& fiducial_data,
                     const VECCPL& test_data, double epsilon, bool log_spacing);

  /// @brief Sets up @p ifos_summary_data and @p ifo_.
  void setup_summary_data(const IfoData& ifo, const VECCPL& fiducial_data,
                          bool log_spacing);

  /// @brief Computes @p r0 and @p r1 for a given signal @p h.
  void compute_waveform_ratios(VECCPL& r0, VECCPL& r1, const VECCPL& h,
                               const VECCPL& data) const;
};

/// @class RelativeBinningBisectionPolarizationsLikelihood
/// @brief Bisection relative-binning likelihood for sky-averaged polarization data.
///
/// Bin placement minimises the total log-likelihood error summed across all
/// active polarization modes. Each mode contributes independently to the
/// summary statistics (A0, A1, B0, B1) and to the log-likelihood via its
/// sky-average factor.
class RelativeBinningBisectionPolarizationsLikelihood : public Likelihood {
 public:
  RelativeBinningBisectionPolarizationsLikelihood(
      const PolarizationData& data,
      const std::vector<VECCPL>& fiducial_modes,
      const std::vector<VECCPL>& test_modes,
      const std::vector<double>& sky_avg_factors,
      const WaveformGenerator& waveform_generator,
      double epsilon, double f_ref, double gmst,
      bool shift_time, bool shift_phase,
      bool log_spacing = false);

  double log_likelihood(gen_params_base<double>* params) const override;

  /// @brief Evaluate RB logL from template modes already extracted at bin edges.
  /// h_at_bins[m] must have exactly bin_inds_.size() entries.
  double log_likelihood_at_waveform(const std::vector<VECCPL>& h_at_bins) const;

  VECDBL get_bin_freqs() const { return bin_freqs_; }
  VECINT get_bin_inds() const { return bin_inds_; }

  const PolarizationData* get_polarization_data() const override {
    return &binned_data_;
  }

  /// Per-mode template-template summary data; B0/B1 are used by the Fisher.
  struct SummaryData {
    VECCPL A0, A1;
    VECDBL B0, B1;
  };

  /// Returns the per-mode summary data (one entry per polarization mode).
  const std::vector<SummaryData>& template_summary() const {
    return summary_data_;
  }

  /// Compute per-mode RB ratios r0[m][b] and r1[m][b] from template modes
  /// at bin edges and the stored fiducial at edges.
  /// Pass derivative waveforms here to obtain normalised derivative ratios
  /// for Fisher computation.
  void compute_waveform_ratios(std::vector<VECCPL>& r0, std::vector<VECCPL>& r1,
                               const std::vector<VECCPL>& h) const;

 private:
  const WaveformGenerator& waveform_generator_;
  const double f_ref_, gmst_;
  const bool shift_time_, shift_phase_;
  const std::vector<double> sky_avg_factors_;

  int number_of_bins_;
  int number_of_modes_;
  std::vector<SummaryData> summary_data_;   // one per mode
  std::vector<VECCPL> fiducial_at_edges_;   // [mode][edge], size = bin_inds_.size() per mode
  PolarizationData binned_data_;            // freqs at bin edges; modes unused

  VECINT bin_inds_;
  VECDBL bin_freqs_;
  VECDBL bin_widths_;
  VECDBL bin_centers_;

  std::pair<int, int> find_min_max_indices(const std::vector<VECCPL>& modes);

  double bin_log_likelihood_error(int left_idx, int right_idx,
                                  const PolarizationData& data,
                                  const std::vector<VECCPL>& fiducial_modes,
                                  const std::vector<VECCPL>& test_modes,
                                  bool log_spacing);

  void bin_bisection(const PolarizationData& data,
                     const std::vector<VECCPL>& fiducial_modes,
                     const std::vector<VECCPL>& test_modes,
                     double epsilon, bool log_spacing);

  void setup_summary_data(const PolarizationData& data,
                          const std::vector<VECCPL>& fiducial_modes,
                          bool log_spacing);
};

}  // namespace RelativeBinning

}  // namespace gw_likelihoods

#endif  // LIKELIHOODS_H
