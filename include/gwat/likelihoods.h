#ifndef LIKELIHOODS_H
#define LIKELIHOODS_H

#include <vector>

#include "waveform_util.h"

namespace GWATLikelihoods {

using CPL = std::complex<double>;
using VECDBL = std::vector<double>;
using VECINT = std::vector<int>;
using VECCPL = std::vector<CPL>;

// Struct for an interferometer's strain and PSD at specific frequencies
struct IfoData {
  // Interferometer strain
  VECCPL strain;
  // Frequency array
  VECDBL freqs;
  // PSD
  VECDBL psd;
  // Signal duration
  double duration;
};

// Struct for each interferometer.
// Holds their binned strain and summary data.
struct SummaryData {
  VECCPL strain;
  // Summary data
  VECCPL A0, A1, B0, B1;
};

/**
 * Base class for adaptive likelihoods
 */
class Likelihood {
 public:
  Likelihood() = default;
  virtual ~Likelihood() = default;

  // Compute the log likelihood -1/2[(h|h) - 2(d|h)] over all detectors
  virtual double log_likelihood(std::string* detectors, int num_detectors,
                                gen_params_base<double>* params,
                                std::string generation_method,
                                bool reuse_WF) = 0;
};

// RELATIVE BINNING

class RelativeBinningPNansatzLikelihood : public Likelihood {
 public:
  RelativeBinningPNansatzLikelihood(
      double chi, double epsilon, const std::vector<IfoData>& ifos_data,
      const CPL* const* fiducial_data, int num_detectors,
      const double* frequencies, int data_length,
      const VECDBL gammas = VECDBL{-5. / 3., -2. / 3., 1., 5. / 3., 7. / 3.});

  double log_likelihood(std::string* detectors, int num_detectors,
                        gen_params_base<double>* params,
                        std::string generation_method, bool reuse_WF);

  VECDBL get_bin_freqs() { return bin_freqs; }

 private:
  double chi, epsilon;
  // Signal duration for tc purposes
  double duration;
  int number_of_bins;
  // Check that bins have been setup
  bool bins_are_setup = false;
  // Interferometers
  std::vector<SummaryData> ifos_fiducial_data;

  // Maximum frequency
  double max_frequency;
  // Array of frequency powers for binning criterion
  VECDBL gammas;
  // Frequency array bin edge indices
  VECINT bin_inds;
  // Frequencies at bin edges
  VECDBL bin_freqs;
  // Distances between bin indices
  VECDBL bin_sizes;
  // Distances between edge frequencies
  VECDBL bin_widths;
  // Frequency at bin centers
  VECDBL bin_centers;

  double find_max_frequency(const CPL* const* fiducial_data_in,
                            const double* frequencies, const int data_length,
                            const int num_detectors);
  void setup_bins(const double* frequencies, const int data_length);
  void setup_fiducial_data(const CPL* const* fiducial_data_in,
                           const int num_detectors);
  void compute_summary_data(const CPL* const* fiducial_data,
                            const std::vector<IfoData>& ifos_data);

  void compute_waveform_ratios(VECCPL& r0, VECCPL& r1, const CPL* h,
                               const SummaryData* fiducial);
  double log_likelihood_per_detector(const CPL* h, const SummaryData* fiducial);
};

// BISECTION RELATIVE BINNING
//
// Bins are placed by recursively bisecting frequency intervals until the
// absolute log-likelihood error within each bin (exact minus approximate,
// evaluated on the test waveform and summed over detectors) falls below
// epsilon.  This makes bin placement data-adaptive rather than PN-ansatz-
// driven.

//! \class RelativeBinningBisectionLikelihood
//! Relative-binning likelihood binned with a bisection algorithm.
class RelativeBinningBisectionLikelihood : public Likelihood {
 public:
  RelativeBinningBisectionLikelihood(const std::vector<IfoData>& ifos_data,
                                     const std::vector<VECCPL>& fiducial_data,
                                     const std::vector<VECCPL>& test_data,
                                     const double* frequencies, double epsilon);

  double log_likelihood(std::string* detectors, int num_detectors,
                        gen_params_base<double>* params,
                        std::string generation_method, bool reuse_WF);

  VECDBL get_bin_freqs() { return bin_freqs; }
  VECINT get_bin_inds() { return bin_inds; }

  // Evaluate RB logL from waveforms already extracted at bin edges.
  // h_at_bins[det] must have exactly bin_inds.size() entries.
  double log_likelihood_from_waveform(const std::vector<VECCPL>& h_at_bins);

 private:
  int number_of_bins;
  std::vector<SummaryData> ifos_summary_data;

  VECINT bin_inds;
  VECDBL bin_freqs;
  VECDBL bin_sizes;
  VECDBL bin_widths;
  VECDBL bin_centers;

  int find_max_index(const std::vector<VECCPL>& fiducial_data);

  // Returns the absolute logL error (exact minus approximate, summed over
  // detectors) for the proposed bin [left_idx, right_idx).
  double bin_log_likelihood_error(int left_idx, int right_idx,
                                  const std::vector<IfoData>& ifos_data,
                                  const std::vector<VECCPL>& fiducial_data,
                                  const std::vector<VECCPL>& test_data);

  void bin_bisection(const std::vector<IfoData>& ifos_data,
                     const std::vector<VECCPL>& fiducial_data,
                     const std::vector<VECCPL>& test_data,
                     const double* frequencies, double epsilon);

  void setup_summary_data(const std::vector<IfoData>& ifos_data,
                          const std::vector<VECCPL>& fiducial_data);

  void compute_waveform_ratios(
      VECCPL& r0, VECCPL& r1, const CPL* h, const SummaryData& fiducial);

  double log_likelihood_per_detector(const CPL* h, const SummaryData& fiducial);
};

/// Initialize a RelativeBinningBisectionLikelihood object from standard
/// C-arrays. GWAT standard for waveform structures are C-arrays so this wraps
/// them to instantiate the likelihood object, which uses C++ vectors.
RelativeBinningBisectionLikelihood RelBinBisectLL_from_c_arrays(
    double** freq, CPL** data, double** psds, CPL** fiducials, CPL** test_wfs,
    int length, int n_detect, double duration, double epsilon);

}  // namespace GWATLikelihoods

#endif  // LIKELIHOODS_H