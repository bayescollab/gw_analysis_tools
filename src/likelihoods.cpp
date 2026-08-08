#include "likelihoods.h"

#include <assert.h>
#include <waveform_util.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <set>
#include <utility>

namespace GWLikelihoods {

using std::conj;

double log_likelihood_in_detector(const VECCPL& h, const IfoData& ifo,
                                  const Quadrature& quad) {
  double hh = quad.inner_product(h, h, ifo.psd);
  double dh = quad.inner_product(ifo.strain, h, ifo.psd);
  return -0.5 * hh + dh;
}

double CoherentLikelihood::log_likelihood(gen_params_base<double>* params,
                                          std::string model) const {
  std::vector<VECCPL> responses =
      create_coherent_GW_detection_reuse_wf(ifos, params, model);

  double ll = 0.0;
  for (size_t i = 0; i < ifos.size(); i++)
    ll += log_likelihood_in_detector(responses[i], ifos[i], quad);

  return ll;
}

namespace RelativeBinning {

/* Print statement */
void RelativeBinningPrinter(std::string message) {
  std::cout << "RELATIVE BINNING:\t";
  std::cout << message << "\n";
}

// ============================================================
// RelativeBinningBisectionLikelihood
// ============================================================

/// @brief Compute the summary data of a bin with the trapezoidal
/// rule.
/// @param A0 [out]
/// @param A1 [out]
/// @param B0 [out]
/// @param B1 [out]
/// @param bin_m          The central frequency point of the bin. Use its
/// logarithm if @p log_spacing is true.
/// @param left_idx       Array index of the left bin edge.
/// @param right_idx      Array index of the right bin edge.
/// @param freqs          Full frequency array.
/// @param psd            Full PSD array.
/// @param h0             Full fiducial array.
/// @param d              Full data array.
/// @param log_spacing    Assume frequency array is log-uniform.
/// @param overall_factor Overall factor of the inner product, usually 4 ∆f.
void compute_summary_data_per_bin_trapezoid_rule(
    CPL& A0, CPL& A1, double& B0, double& B1, double bin_m, int left_idx,
    int right_idx, const VECDBL& freqs, const VECDBL& psd, const VECCPL& h0,
    const VECCPL& d, bool log_spacing, double overall_factor = 1.0) {
  A0 = A1 = 0.; B0 = B1 = 0.;

  for (int i = left_idx + 1; i < right_idx; i++) {
    double delta = (log_spacing ? std::log10(freqs[i]) : freqs[i]) - bin_m;

    CPL A_fac = conj(h0[i]) / psd[i];
    if (log_spacing) A_fac *= freqs[i];
    double B_fac = std::real(h0[i] * A_fac);
    A_fac *= d[i];
    A0 += A_fac;
    A1 += A_fac * delta;
    B0 += B_fac;
    B1 += B_fac * delta;
  }

  for (int i : {left_idx, right_idx}) {
    double delta = (log_spacing ? std::log10(freqs[i]) : freqs[i]) - bin_m;

    CPL A_fac = conj(h0[i]) / psd[i];
    if (log_spacing) A_fac *= freqs[i];
    double B_fac = std::real(h0[i] * A_fac);
    A_fac *= d[i];
    A0 += 0.5 * A_fac;
    A1 += 0.5 * A_fac * delta;
    B0 += 0.5 * B_fac;
    B1 += 0.5 * B_fac * delta;
  }

  if (overall_factor != 1.0) {
    A0 *= overall_factor;
    A1 *= overall_factor;
    B0 *= overall_factor;
    B1 *= overall_factor;
  }
}

RelativeBinningBisectionLikelihood::RelativeBinningBisectionLikelihood(
    const std::vector<IfoData>& ifos_data,
    const std::vector<VECCPL>& fiducial_data,
    const std::vector<VECCPL>& test_data, const double* frequencies,
    double epsilon, bool log_spacing) {
  std::cout << "\nRELATIVE BINNING (BISECTION) INITIALIZING\n";

  bin_bisection(ifos_data, fiducial_data, test_data, frequencies, epsilon,
                log_spacing);
  RelativeBinningPrinter(std::to_string(number_of_bins) + " bins setup");

  setup_summary_data(ifos_data, fiducial_data, log_spacing);

  double logL = 0.;
  int num_edges = (int)bin_inds.size();
  VECCPL ht(num_edges);
  for (size_t i = 0; i < ifos_summary_data.size(); i++) {
    for (int k = 0; k < num_edges; k++) ht[k] = test_data[i][bin_inds[k]];

    logL += log_likelihood_per_detector(ht, i);
  }
  RelativeBinningPrinter("logL of test data: " + std::to_string(logL));
}

int RelativeBinningBisectionLikelihood::find_max_index(
    const std::vector<VECCPL>& fiducial_data) {
  int max_idx = (int)fiducial_data.front().size() - 1;

  for (const auto& data : fiducial_data) {
    // Skip leading zeros: some waveforms (e.g. EFPE) are exactly zero at
    // f = f_start due to a floating-point boundary in the SPA frequency check.
    auto first_nz = std::find_if(data.begin(), data.end(),
                                 [](const CPL& v) { return std::norm(v) > 0.; });
    if (first_nz == data.end()) continue;  // all zeros — no cutoff from this IFO

    // Find the first zero after the waveform has become non-zero (merger cutoff).
    auto it = std::find_if(first_nz, data.end(),
                           [](const CPL& v) { return std::norm(v) == 0.; });
    if (it != data.end()) {
      int zero_idx = (int)std::distance(data.begin(), it);
      max_idx = std::min(max_idx, zero_idx - 1);
    }
  }
  // Require at least index 1 so bin_freqs always has at least two distinct edges.
  return std::max(max_idx, 1);
}

double RelativeBinningBisectionLikelihood::bin_log_likelihood_error(
    int left_idx, int right_idx, const std::vector<IfoData>& ifos_data,
    const std::vector<VECCPL>& fiducial_data,
    const std::vector<VECCPL>& test_data, bool log_spacing) {
  double exact_logL = 0., approx_logL = 0.;
  int num_det = (int)ifos_data.size();

  for (int det = 0; det < num_det; det++) {
    const IfoData& ifo = ifos_data[det];
    const VECCPL& h0 = fiducial_data[det];
    const VECCPL& ht = test_data[det];
    double weight =
        (log_spacing
             ? LOG10 * (std::log10(ifo.freqs[1]) - std::log10(ifo.freqs[0]))
             : (ifo.freqs[1] - ifo.freqs[0])) *
        4.0;
    double f_m = (log_spacing ? (std::log10(ifo.freqs[left_idx]) +
                                 std::log10(ifo.freqs[right_idx]))
                              : (ifo.freqs[left_idx] + ifo.freqs[right_idx])) *
                 0.5;

    CPL A0 = 0., A1 = 0.; double B0 = 0., B1 = 0.;

    compute_summary_data_per_bin_trapezoid_rule(
        A0, A1, B0, B1, f_m, left_idx, right_idx, ifo.freqs, ifo.psd, h0,
        ifo.strain, log_spacing, weight);
    double exact_dh = 0., exact_hh = 0.;
    for (int j = left_idx + 1; j < right_idx; j++) {
      exact_dh +=
          log_spacing
              ? real(ifo.strain[j] * conj(ht[j]) * ifo.freqs[j]) / ifo.psd[j]
              : real(ifo.strain[j] * conj(ht[j])) / ifo.psd[j];
      exact_hh += log_spacing ? std::norm(ht[j]) * ifo.freqs[j] / ifo.psd[j]
                              : std::norm(ht[j]) / ifo.psd[j];
    }
    for (int j : {left_idx, right_idx}) {
      exact_dh +=
          0.5 *
          (log_spacing
               ? real(ifo.strain[j] * conj(ht[j]) * ifo.freqs[j]) / ifo.psd[j]
               : real(ifo.strain[j] * conj(ht[j])) / ifo.psd[j]);
      exact_hh +=
          0.5 * (log_spacing ? std::norm(ht[j]) * ifo.freqs[j] / ifo.psd[j]
                             : std::norm(ht[j]) / ifo.psd[j]);
    }
    exact_dh *= weight;
    exact_hh *= weight;

    // Ratio r(f) = h_test(f) / h0(f) at bin edges
    CPL r_left =
        (h0[left_idx] != CPL(0.)) ? ht[left_idx] / h0[left_idx] : CPL(0.);
    CPL r_right =
        (h0[right_idx] != CPL(0.)) ? ht[right_idx] / h0[right_idx] : CPL(0.);
    double bin_width = log_spacing ? std::log10(ifo.freqs[right_idx]) -
                                         std::log10(ifo.freqs[left_idx])
                                   : ifo.freqs[right_idx] - ifo.freqs[left_idx];
    CPL r0 = 0.5 * (r_left + r_right);
    CPL r1 = (bin_width > 0.) ? (r_right - r_left) / bin_width : CPL(0.);

    double approx_dh = real(A0 * conj(r0) + A1 * conj(r1));
    double approx_hh =
        B0 * std::norm(r0) + 2. * B1 * real(r0 * conj(r1));

    exact_logL += -0.5 * exact_hh + exact_dh;
    approx_logL += -0.5 * approx_hh + approx_dh;
  }

  return std::abs(exact_logL - approx_logL);
}

void RelativeBinningBisectionLikelihood::bin_bisection(
    const std::vector<IfoData>& ifos_data,
    const std::vector<VECCPL>& fiducial_data,
    const std::vector<VECCPL>& test_data, const double* frequencies,
    double epsilon, bool log_spacing) {
  int max_idx = find_max_index(fiducial_data);
  RelativeBinningPrinter("Max frequency: " +
                         std::to_string(frequencies[max_idx]));

  // Find the first index where ALL detectors have non-zero fiducial amplitude.
  // Waveforms like EFPE are exactly zero at f = f_start due to a floating-point
  // boundary in the SPA evaluation; starting bins there gives NaN ratios.
  int min_idx = 0;
  for (const auto& data : fiducial_data) {
    auto it = std::find_if(data.begin(), data.end(),
                           [](const CPL& v) { return std::norm(v) > 0.; });
    if (it != data.end())
      min_idx = std::max(min_idx, (int)std::distance(data.begin(), it));
  }
  min_idx = std::min(min_idx, max_idx - 1);  // always leave room for ≥1 bin

  // Collect sorted bin boundary indices via iterative DFS bisection.
  // The set always contains min_idx (left edge of the first bin); each accepted
  // interval [left, right) contributes its right index.
  std::set<int> boundaries;
  boundaries.insert(min_idx);

  std::vector<std::pair<int, int>> stack;
  stack.push_back(std::make_pair(min_idx, max_idx));

  while (!stack.empty()) {
    int left = stack.back().first;
    int right = stack.back().second;
    stack.pop_back();

    if (right - left <= 1) {
      boundaries.insert(right);
      continue;
    }

    double err = bin_log_likelihood_error(left, right, ifos_data, fiducial_data,
                                          test_data, log_spacing);

    if (err <= epsilon) {
      boundaries.insert(right);
    } else {
      int mid = (left + right) / 2;
      stack.push_back(std::make_pair(left, mid));
      stack.push_back(std::make_pair(mid, right));
    }
  }

  for (int idx : boundaries) {
    bin_inds.push_back(idx);
    bin_freqs.push_back(frequencies[idx]);
  }

  number_of_bins = (int)bin_inds.size() - 1;
  for (int i = 1; i < (int)bin_inds.size(); i++) {
    if (log_spacing) {
      bin_widths.push_back(std::log10(bin_freqs[i]) -
                           std::log10(bin_freqs[i - 1]));
      bin_centers.push_back(
          (std::log10(bin_freqs[i]) + std::log10(bin_freqs[i - 1])) * 0.5);
    } else {
      bin_widths.push_back(bin_freqs[i] - bin_freqs[i - 1]);
      bin_centers.push_back((bin_freqs[i] + bin_freqs[i - 1]) * 0.5);
    }
  }
}

void RelativeBinningBisectionLikelihood::setup_summary_data(
    const std::vector<IfoData>& ifos_data,
    const std::vector<VECCPL>& fiducial_data, bool log_spacing) {
  int num_det = (int)ifos_data.size();
  for (int det = 0; det < num_det; det++) {
    SummaryData det_summary_data;
    const IfoData& ifo = ifos_data[det];
    IfoData new_ifo = ifo;  // copy all fields; clear the data arrays below
    new_ifo.freqs.clear();
    new_ifo.psd.clear();
    new_ifo.strain.clear();
    const VECCPL& h0 = fiducial_data[det];
    double weight =
        (log_spacing
             ? LOG10 * (std::log10(ifo.freqs[1]) - std::log10(ifo.freqs[0]))
             : (ifo.freqs[1] - ifo.freqs[0])) *
        4.0;

    for (int b = 0; b < number_of_bins; b++) {
      int start = bin_inds[b];
      int end = bin_inds[b + 1];

      new_ifo.strain.push_back(h0[start]);
      new_ifo.psd.push_back(ifo.psd[start]);
      new_ifo.freqs.push_back(ifo.freqs[start]);

      CPL A0 = 0., A1 = 0.; double B0 = 0., B1 = 0.;
      compute_summary_data_per_bin_trapezoid_rule(
          A0, A1, B0, B1, bin_centers[b], start, end, ifo.freqs, ifo.psd, h0,
          ifo.strain, log_spacing, weight);
      det_summary_data.A0.push_back(A0);
      det_summary_data.A1.push_back(A1);
      det_summary_data.B0.push_back(B0);
      det_summary_data.B1.push_back(B1);
    }
    new_ifo.strain.push_back(h0[bin_inds.back()]);
    new_ifo.psd.push_back(ifo.psd[bin_inds.back()]);
    new_ifo.freqs.push_back(ifo.freqs[bin_inds.back()]);

    ifos_summary_data.push_back(det_summary_data);
    ifos.push_back(new_ifo);
  }
}

void RelativeBinningBisectionLikelihood::compute_waveform_ratios(
    VECCPL& r0, VECCPL& r1, const VECCPL& h, const VECCPL& data) const {
  r0.clear(); r1.clear();

  CPL ratio_left = h[0] / data.front();
  CPL ratio_right;

  for (int i = 0; i < number_of_bins; i++) {
    ratio_right = h[i + 1] / data[i + 1];

    r0.push_back(0.5 * (ratio_right + ratio_left));
    r1.push_back((ratio_right - ratio_left) / bin_widths[i]);

    ratio_left = ratio_right;
  }
}

double RelativeBinningBisectionLikelihood::log_likelihood_per_detector(
    const VECCPL& h, size_t i) const {
  double d_h = 0., h_h = 0.;
  const auto fiducial = ifos_summary_data[i];
  VECCPL r0, r1;
  compute_waveform_ratios(r0, r1, h, ifos[i].strain);

  for (int b = 0; b < number_of_bins; b++) {
    d_h += real(fiducial.A0[b] * conj(r0[b]) + fiducial.A1[b] * conj(r1[b]));
    h_h += fiducial.B0[b] * std::norm(r0[b]) +
           2. * fiducial.B1[b] * real(r0[b] * conj(r1[b]));
  }

  return -0.5 * h_h + d_h;
}

double RelativeBinningBisectionLikelihood::log_likelihood_from_waveform(
    const std::vector<VECCPL>& h_at_bins) {
  int num_det = (int)ifos_summary_data.size();
  double logL = 0.;
  for (int i = 0; i < num_det; i++)
    logL += log_likelihood_per_detector(h_at_bins[i], i);
  return logL;
}

double RelativeBinningBisectionLikelihood::log_likelihood(
    gen_params_base<double>* params, std::string model) const {
  std::vector<VECCPL> responses =
      create_coherent_GW_detection_reuse_wf(ifos, params, model);

  double logL = 0.;
  for (size_t i = 0; i < ifos_summary_data.size(); i++)
    logL += log_likelihood_per_detector(responses[i], i);

  return logL;
}

RelativeBinningBisectionLikelihood RelBinBisectLL_from_c_arrays(
    double** freq, CPL** data, double** psds, CPL** fiducials, CPL** test_wfs,
    int length, int n_detect, double duration, double epsilon,
    bool log_spacing) {
  std::vector<IfoData> ifos_data;
  std::vector<VECCPL> fiducial_vecs;
  std::vector<VECCPL> test_vecs;

  for (int i = 0; i < n_detect; i++) {
    IfoData ifo;
    ifo.freqs = VECDBL(freq[i], freq[i] + length);
    ifo.psd = VECDBL(psds[i], psds[i] + length);
    ifo.strain = VECCPL(data[i], data[i] + length);
    ifos_data.push_back(ifo);

    fiducial_vecs.push_back(VECCPL(fiducials[i], fiducials[i] + length));
    test_vecs.push_back(VECCPL(test_wfs[i], test_wfs[i] + length));
  }

  return RelativeBinningBisectionLikelihood(ifos_data, fiducial_vecs, test_vecs,
                                            freq[0], epsilon, log_spacing);
}

}  // namespace RelativeBinning

}  // namespace GWLikelihoods
