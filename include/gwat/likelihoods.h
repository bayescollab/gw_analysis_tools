#ifndef LIKELIHOODS_H
#define LIKELIHOODS_H

#include <vector>

#include "waveform_util.h"

namespace gw_likelihoods {

///@class Likelihood
///@brief Base class for likelihoods
class Likelihood {
 public:
  Likelihood() = default;
  virtual ~Likelihood() = default;

  /// @brief Compute the log likelihood -1/2[(h|h) - 2(d|h)] over all detectors
  /// @param params             Model parameters.
  /// @param generation_method  Model name.
  virtual double log_likelihood(gen_params_base<double>* params,
                                std::string generation_method) const = 0;

  /// @brief The Quadrature object used internally, if any.
  virtual const Quadrature* quadrature() const { return nullptr; }

  /// @brief The detector data used by this likelihood.
  /// For precomputed likelihoods (e.g. relative binning), returns the binned
  /// representation rather than the full-resolution detector data.
  virtual const std::vector<IfoData>& get_ifos() const = 0;
};

/// @class CoherentLikelihood
/// @brief The most simple likelihood class. Uses
/// create_coherent_GW_detection_reuse_wf.
class CoherentLikelihood : public Likelihood {
 private:
  const Quadrature& quad;
  const std::vector<IfoData>& ifos;

 public:
  CoherentLikelihood(const std::vector<IfoData>& is, const Quadrature& q)
      : ifos(is), quad(q) {};

  /// @details Compute the unmarginalized log-likelihood.
  double log_likelihood(gen_params_base<double>* params,
                        std::string generation_method) const override;

  const Quadrature* quadrature() const override { return &quad; }
  const std::vector<IfoData>& get_ifos() const override { return ifos; }
};

namespace RelativeBinning {

// BISECTION RELATIVE BINNING
//
// Bins are placed by recursively bisecting frequency intervals until the
// absolute log-likelihood error within each bin (exact minus approximate,
// evaluated on the test waveform and summed over detectors) falls below
// epsilon.  This makes bin placement data-adaptive rather than PN-ansatz-
// driven.

// Struct for each interferometer.
// Holds their binned strain and summary data.
struct SummaryData {
  // Summary data
  VECCPL A0, A1;
  VECDBL B0, B1;
};

///@class RelativeBinningBisectionLikelihood
///@brief Relative-binning likelihood binned with a bisection algorithm.
class RelativeBinningBisectionLikelihood : public Likelihood {
 public:
  RelativeBinningBisectionLikelihood(const std::vector<IfoData>& ifos_data,
                                     const std::vector<VECCPL>& fiducial_data,
                                     const std::vector<VECCPL>& test_data,
                                     const double* frequencies, double epsilon,
                                     bool log_spacing = false);

  double log_likelihood(gen_params_base<double>* params,
                        std::string generation_method) const override;

  VECDBL get_bin_freqs() { return bin_freqs; }
  VECINT get_bin_inds() { return bin_inds; }

  // Evaluate RB logL from waveforms already extracted at bin edges.
  // h_at_bins[det] must have exactly bin_inds.size() entries.
  double log_likelihood_from_waveform(const std::vector<VECCPL>& h_at_bins);

  /// @brief Returns the binned detector data (bin-edge frequencies, fiducial
  /// strain at edges, and PSD at edges). Populated during construction.
  const std::vector<IfoData>& get_ifos() const override { return ifos; }

 private:
  int number_of_bins;
  std::vector<SummaryData> ifos_summary_data;
  std::vector<IfoData> ifos;

  /// @brief Indices of the bin edges corresponding to the full-resolution grid,
  /// including the final right edge.
  VECINT bin_inds;
  /// @brief The frequencies at each bin edge.
  VECDBL bin_freqs;
  /// @brief The spacing between bin edges. Size: bin_inds.size()-1
  VECDBL bin_widths;
  /// @brief The center of the bin.
  VECDBL bin_centers;

  int find_max_index(const std::vector<VECCPL>& fiducial_data);

  /// @brief Returns the absolute logL error (exact minus approximate, summed
  /// over detectors) for the proposed bin [left_idx, right_idx).
  double bin_log_likelihood_error(int left_idx, int right_idx,
                                  const std::vector<IfoData>& ifos_data,
                                  const std::vector<VECCPL>& fiducial_data,
                                  const std::vector<VECCPL>& test_data,
                                  bool log_spacing);

  /// @brief Bisection algorithm. Populates the bin_* vectors.
  void bin_bisection(const std::vector<IfoData>& ifos_data,
                     const std::vector<VECCPL>& fiducial_data,
                     const std::vector<VECCPL>& test_data,
                     const double* frequencies, double epsilon,
                     bool log_spacing);

  /// @brief Sets up @p ifos_summary_data and @p ifos.
  void setup_summary_data(const std::vector<IfoData>& ifos_data,
                          const std::vector<VECCPL>& fiducial_data,
                          bool log_spacing);

  /// @brief Computes @p r0 and @p r1 for a given signal @p h.
  void compute_waveform_ratios(VECCPL& r0, VECCPL& r1, const VECCPL& h,
                               const VECCPL& data) const;

  /// @brief Compute the log-likelihood for detector @p i, the i-th detector in
  /// @ref ifos.
  double log_likelihood_per_detector(const VECCPL& h, size_t i) const;
};

/// Initialize a RelativeBinningBisectionLikelihood object from standard
/// C-arrays. GWAT standard for waveform structures are C-arrays so this wraps
/// them to instantiate the likelihood object, which uses C++ vectors.
RelativeBinningBisectionLikelihood RelBinBisectLL_from_c_arrays(
    double** freq, CPL** data, double** psds, CPL** fiducials, CPL** test_wfs,
    int length, int n_detect, double duration, double epsilon,
    bool log_spacing = false);

}  // namespace RelativeBinning

}  // namespace gw_likelihoods

#endif  // LIKELIHOODS_H