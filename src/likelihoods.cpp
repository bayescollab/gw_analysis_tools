#include "likelihoods.h"

#include <assert.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <set>
#include <utility>

namespace GWATLikelihoods {

using std::conj;

/* Print statement */
void RelativeBinningPrinter(std::string message) {
  std::cout << "RELATIVE BINNING:\t";
  std::cout << message << "\n";
}

/**
 * Initialize relative binning
 */
RelativeBinningPNansatzLikelihood::RelativeBinningPNansatzLikelihood(
    double chi,                            /**< size of PN perturbation */
    double epsilon,                        /**< size of dephasing */
    const std::vector<IfoData>& ifos_data, /**< Vector of signal data */
    const CPL* const*
        fiducial_data, /**< Strain of fiducial data for each interferometer */
    int num_detectors, /**< Number of detectors */
    const double*
        frequencies,    /**< Frequency array shared among all interferometers */
    int data_length,    /**< length of each fiducial array */
    const VECDBL gammas /**< Vector of gammas for binning criterion */
) {
  std::cout << "RELATIVE BINNING INITIALIZING\n";
  this->chi = chi;
  this->epsilon = epsilon;
  this->gammas = gammas;

  std::cout << "\tUsing gamma = { ";
  for (auto& gamma : gammas) {
    std::cout << gamma << " ";
  }
  std::cout << "}\n";

  duration = ifos_data[0].duration;
  max_frequency = find_max_frequency(fiducial_data, frequencies, data_length,
                                     num_detectors);
  RelativeBinningPrinter("Max frequency: " + std::to_string(max_frequency));

  setup_bins(frequencies, data_length);
  RelativeBinningPrinter(std::to_string(number_of_bins) + " bins setup");

  setup_fiducial_data(fiducial_data, num_detectors);
  RelativeBinningPrinter("Fiducial data setup");

  compute_summary_data(fiducial_data, ifos_data);
  RelativeBinningPrinter("Summary data setup");

  // Calculate likelihood
  double logL = 0.;
  CPL* data;
  for (int i = 0; i < num_detectors; i++) {
    data = ifos_fiducial_data.at(i).strain.data();
    logL += log_likelihood_per_detector(data, &(ifos_fiducial_data.at(i)));
  }
  RelativeBinningPrinter("Fiducial likelihood: " + std::to_string(logL));
}

/**
 * Find maximum frequency at which we can evaluate the likelihood,
 * which will determine the binning.
 * Defined as the last frequency where the fiducial data is non-zero
 */
double RelativeBinningPNansatzLikelihood::find_max_frequency(
    const CPL* const* fiducial_data_in, const double* frequencies,
    const int data_length, const int num_detectors) {
  VECCPL fiducial_data;
  VECCPL::reverse_iterator ref;
  int idx;
  double maxFreq = frequencies[data_length - 1];

  for (int d = 0; d < num_detectors; d++) {
    fiducial_data =
        VECCPL(fiducial_data_in[d], fiducial_data_in[d] + data_length);

    ref = std::find_if(fiducial_data.rbegin(), fiducial_data.rend(),
                       [](CPL& dat) { return dat != CPL(0.); });
    if (ref == fiducial_data.rbegin()) {
      // If the last data point is non-zero
      continue;
    }

    idx = std::distance(ref, fiducial_data.rend()) - 1;
    maxFreq = fmin(frequencies[idx], maxFreq);
  }

  return maxFreq;
}

/**
 * Setup frequency bins based on PN ansatz.
 * See Binning Criterion in arxiv:2312.06009
 */
void RelativeBinningPNansatzLikelihood::setup_bins(
    const double* frequencies, /**< Frequency array */
    const int data_length      /**< Size of frequency array */
) {
  double min_frequency;
  double perturb_size = GWAT_TWOPI * chi;

  // Find the frequency range over all detectors
  // Assumes array is in sequential order
  min_frequency = frequencies[0];

  // Eqs. 15-16 of 2312.06009
  // To optimize computing the sum in Eq. 16 we do not solve for \Delta_\alpha
  // but rather find the f_{k*} for each k
  VECDBL freqs_gamma;
  VECINT signs_gamma;
  double freq_gamma;
  bool negative_gamma;
  for (const double& gamma : gammas) {
    negative_gamma = gamma < 0;
    freq_gamma = negative_gamma ? min_frequency : max_frequency;
    freqs_gamma.push_back(freq_gamma);
    signs_gamma.push_back(negative_gamma ? -1 : 1);
  }

  // Eq. 16, second line instead of first
  VECDBL d_phis;
  double d_phi_f, f;
  size_t k;
  for (int i = 0; i < data_length; i++) {
    f = frequencies[i];
    if (f > max_frequency) {
      break;
    }

    d_phi_f = 0.;

    for (k = 0; k < gammas.size(); k++) {
      d_phi_f += signs_gamma[k] * pow(f / freqs_gamma[k], gammas[k]);
    }

    d_phis.push_back(perturb_size * d_phi_f);
  }

  // \Delta\Psi(f) - \Delta\Psi(f_{\min})
  VECDBL d_phi_from_start;
  for (const double& d_phi : d_phis) {
    d_phi_from_start.push_back(d_phi - d_phis[0]);
  }

  // Number of bins
  int num_bins =
      static_cast<int>(std::floor(d_phi_from_start.back() / epsilon));

  // Find bin edges
  VECDBL::iterator bin_itr;
  // last_* variables to avoid starting the find_if searches from the beginning
  const double* last_base_ptr = frequencies;
  VECDBL::iterator last_itr = d_phi_from_start.begin();
  int bin_ind, last_ind = -1;
  double d_phi_lower, bin_freq;
  double d_phi_lower_den = d_phi_from_start.back() / (double)num_bins;
  // Find in integer increments of epsilon where a bin will have reached a
  // dephasing on the order of epsilon
  for (int i = 0; i < num_bins + 1; i++) {
    d_phi_lower = i * d_phi_lower_den;

    bin_itr = std::find_if(
        last_itr, d_phi_from_start.end(),
        [d_phi_lower](double& d_phi) { return d_phi >= d_phi_lower; });
    bin_ind = std::distance(d_phi_from_start.begin(), bin_itr);
    if (bin_ind == last_ind) {
      continue;
    }

    last_itr = bin_itr;
    last_ind = bin_ind;
    bin_freq = frequencies[bin_ind];
    last_base_ptr =
        std::find_if(last_base_ptr, frequencies + data_length,
                     [bin_freq](double f) { return f >= bin_freq; });
    bin_ind = std::distance(frequencies, last_base_ptr);

    bin_inds.push_back(bin_ind);
    bin_freqs.push_back(bin_freq);
  }

  // Set bin info
  number_of_bins = bin_inds.size() - 1;
  for (int i = 1; i < bin_inds.size(); i++) {
    bin_sizes.push_back(bin_inds[i] - bin_inds[i - 1]);
    bin_widths.push_back(bin_freqs[i] - bin_freqs[i - 1]);
    bin_centers.push_back((bin_freqs[i] + bin_freqs[i - 1]) / 2);
  }

  bins_are_setup = true;
}

/**
 * Store fiducial waveforms at bin edges
 * Saved within ifos_fiducial_data
 */
void RelativeBinningPNansatzLikelihood::setup_fiducial_data(
    const CPL* const*
        fiducial_data_in,   /**< Fiducial data of each interferometer */
    const int num_detectors /** Number of detectors */
) {
  assert(bins_are_setup);

  const CPL* det_wf;
  int d;

  for (d = 0; d < num_detectors; d++) {
    SummaryData ifo;
    det_wf = fiducial_data_in[d];

    for (int& ind : bin_inds) {
      ifo.strain.push_back(det_wf[ind]);
    }

    ifos_fiducial_data.push_back(ifo);
  }
}

/**
 * Compute the A0, A1, B0, and B2 summary data for each bin in each
 * interferometer. See Eq. 5 of 2312.06009
 */
void RelativeBinningPNansatzLikelihood::compute_summary_data(
    const CPL* const* fiducial_data, const std::vector<IfoData>& ifos_data) {
  const IfoData* ifo;
  const CPL* fiducial;
  SummaryData* ifo_fiducial;
  VECINT ifo_bin_ends;
  VECCPL data, h0;
  VECDBL psd, freqs;
  int start_ind, end_ind;
  CPL A_fac, B_fac, A0_b, A1_b, B0_b, B1_b;
  double delta_f, inner_product_weight;

  // Compute summary data for each interferometer
  for (int i = 0; i < ifos_data.size(); i++) {
    // The ifo holding the data
    ifo = &(ifos_data.at(i));
    // The fiducial data
    fiducial = fiducial_data[i];
    // The ifo to store the summary data
    ifo_fiducial = &ifos_fiducial_data.at(i);

    inner_product_weight = 4. / ifo->duration;

    ifo_bin_ends = VECINT(bin_inds);
    // Cover the last frequency in the bin array
    ifo_bin_ends.end() += 1;

    // Compute summary data per bin
    for (int b = 0; b < number_of_bins; b++) {
      // Get the starting and ending indices of the bin
      start_ind = ifo_bin_ends[b];
      end_ind = ifo_bin_ends[b + 1];

      // Grab the data for each bin
      data = VECCPL(ifo->strain.begin() + start_ind,
                    ifo->strain.begin() + end_ind);
      psd = VECDBL(ifo->psd.begin() + start_ind, ifo->psd.begin() + end_ind);
      h0 = VECCPL(fiducial + start_ind, fiducial + end_ind);
      freqs =
          VECDBL(ifo->freqs.begin() + start_ind, ifo->freqs.begin() + end_ind);

      // Compute the Riemann sums
      A0_b = A1_b = B0_b = B1_b = 0.;
      for (int f = 0; f < freqs.size(); f++) {
        A_fac = conj(h0[f]) / psd[f];
        B_fac = h0[f] * A_fac;
        A_fac *= data[f];
        delta_f = freqs[f] - bin_centers[b];

        A0_b += A_fac;
        A1_b += A_fac * delta_f;
        B0_b += B_fac;
        B1_b += B_fac * delta_f;
      }
      // Store the data
      ifo_fiducial->A0.push_back(inner_product_weight * A0_b);
      ifo_fiducial->A1.push_back(inner_product_weight * A1_b);
      ifo_fiducial->B0.push_back(inner_product_weight * B0_b);
      ifo_fiducial->B1.push_back(inner_product_weight * B1_b);
    }
  }
}

void RelativeBinningPNansatzLikelihood::compute_waveform_ratios(
    VECCPL& r0,   /**< [out] Array of r0 coefficients in each bin */
    VECCPL& r1,   /**< [out] Array of r1 coefficients in each bin */
    const CPL* h, /**< Template waveform evaluated at bin edges */
    const SummaryData* fiducial /**< Interferometer data */
) {
  // Ratios at left edge
  CPL ratio_left = h[0] / fiducial->strain.front();
  // Ratio at right edge
  CPL ratio_right;

  for (int i = 0; i < number_of_bins; i++) {
    ratio_right = h[i + 1] / fiducial->strain[i + 1];

    r0.push_back(0.5 * (ratio_right + ratio_left));
    r1.push_back((ratio_right - ratio_left) / bin_widths[i]);

    // For next bin
    ratio_left = ratio_right;
  }
}

double RelativeBinningPNansatzLikelihood::log_likelihood_per_detector(
    const CPL*
        h, /**< Template waveform. Must be evaluated at the bin edges only */
    const SummaryData* fiducial /**< Interferometer data */
) {
  double d_h = 0.;
  double h_h = 0.;

  // Obtain the waveform ratios
  VECCPL r0, r1;
  compute_waveform_ratios(r0, r1, h, fiducial);

  for (int b = 0; b < number_of_bins; b++) {
    // Inner product of the data with the template
    d_h += real(fiducial->A0[b] * conj(r0[b]) + fiducial->A1[b] * conj(r1[b]));
    // Inner product of the template
    h_h += real(fiducial->B0[b] * r0[b] * conj(r0[b]) +
                2. * fiducial->B1[b] * real(r0[b] * conj(r1[b])));
  }

  return -0.5 * h_h + d_h;
}

double RelativeBinningPNansatzLikelihood::log_likelihood(
    std::string* detectors,          /**< Detector names */
    int num_detectors,               /**< Number of detectors */
    gen_params_base<double>* params, /**< Template parameters */
    std::string generation_method,   /**< Template model name */
    bool reuse_WF /**< Option to obtain detector responses using the same
                     waveform */
) {
  // Set up template arrays a là MCMC_likelihood_extrinsic (mcmc_gw.cpp)
  int i;
  int* data_lengths = new int[num_detectors];
  CPL** responses = new CPL*[num_detectors];
  double** frequencies = new double*[num_detectors];

  for (int i = 0; i < num_detectors; i++) {
    data_lengths[i] = bin_freqs.size();
    responses[i] = new CPL[data_lengths[i]];
    frequencies[i] = bin_freqs.data();
  }

  // Obtain the data
  if (num_detectors == 1) {
    create_single_GW_detection(responses[0], detectors[0], frequencies[0],
                               data_lengths[0], params, generation_method);
  } else {
    create_coherent_GW_detection(detectors, num_detectors, frequencies,
                                 data_lengths, reuse_WF, params,
                                 generation_method, responses);
  }

  // Calculate likelihood
  double logL = 0.;
  for (i = 0; i < num_detectors; i++) {
    logL += log_likelihood_per_detector(responses[i], &(ifos_fiducial_data[i]));
  }

  // Clean up
  for (i = 0; i < num_detectors; i++) {
    delete[] responses[i];
  }
  delete[] frequencies;
  delete[] responses;
  delete[] data_lengths;

  return logL;
}

// ============================================================
// RelativeBinningBisectionLikelihood
// ============================================================

RelativeBinningBisectionLikelihood::RelativeBinningBisectionLikelihood(
    const std::vector<IfoData>& ifos_data,
    const std::vector<VECCPL>& fiducial_data,
    const std::vector<VECCPL>& test_data, const double* frequencies,
    double epsilon) {
  std::cout << "\nRELATIVE BINNING (BISECTION) INITIALIZING\n";

  if (bin_bisection(ifos_data, fiducial_data, test_data, frequencies,
                    epsilon)) {
    RelativeBinningPrinter(std::to_string(number_of_bins) + " bins setup");

    setup_summary_data(ifos_data, fiducial_data);

    double logL = 0.;
    int num_edges = (int)bin_inds.size();
    CPL* ht = new CPL[num_edges];
    for (size_t i = 0; i < ifos_summary_data.size(); i++) {
      for (int k = 0; k < num_edges; k++) ht[k] = test_data[i][bin_inds[k]];

      logL += log_likelihood_per_detector(ht, ifos_summary_data.at(i));
    }
    RelativeBinningPrinter("Test logL: " + std::to_string(logL));

    delete[] ht;
  }
}

int RelativeBinningBisectionLikelihood::find_max_index(
    const std::vector<VECCPL>& fiducial_data) {
  int max_idx = (int)fiducial_data.front().size() - 1;

  for (const auto& data : fiducial_data) {
    auto it = std::find_if(data.begin(), data.end(),
                           [](const CPL& v) { return std::norm(v) == 0.; });

    if (it != data.end()) {
      int zero_idx = (int)std::distance(data.begin(), it);
      max_idx = std::min(max_idx, zero_idx - 1);
    }
  }
  return max_idx;
}

double RelativeBinningBisectionLikelihood::bin_log_likelihood_error(
    int left_idx, int right_idx, const std::vector<IfoData>& ifos_data,
    const std::vector<VECCPL>& fiducial_data,
    const std::vector<VECCPL>& test_data) {
  double exact_logL = 0., approx_logL = 0.;
  int num_det = (int)ifos_data.size();

  for (int det = 0; det < num_det; det++) {
    const IfoData& ifo = ifos_data[det];
    const VECCPL& h0 = fiducial_data[det];
    const VECCPL& ht = test_data[det];
    double weight = 4. / ifo.duration;
    double f_m = (ifo.freqs[left_idx] + ifo.freqs[right_idx]) * 0.5;

    CPL A0 = 0., A1 = 0., B0 = 0., B1 = 0.;
    double exact_dh = 0., exact_hh = 0.;

    for (int j = left_idx; j < right_idx; j++) {
      double delta = ifo.freqs[j] - f_m;
      double psd = ifo.psd[j];
      CPL strain = ifo.strain[j];

      // Summary data from h0 (fiducial)
      CPL A_fac = conj(h0[j]) / psd;
      CPL B_fac = h0[j] * A_fac;
      A_fac *= strain;
      A0 += A_fac;
      A1 += A_fac * delta;
      B0 += B_fac;
      B1 += B_fac * delta;

      // Exact inner products with h_test
      exact_dh += real(strain * conj(ht[j])) / psd;
      exact_hh += std::norm(ht[j]) / psd;
    }
    A0 *= weight;
    A1 *= weight;
    B0 *= weight;
    B1 *= weight;
    exact_dh *= weight;
    exact_hh *= weight;

    // Ratio r(f) = h_test(f) / h0(f) at bin edges
    CPL r_left =
        (h0[left_idx] != CPL(0.)) ? ht[left_idx] / h0[left_idx] : CPL(0.);
    CPL r_right =
        (h0[right_idx] != CPL(0.)) ? ht[right_idx] / h0[right_idx] : CPL(0.);
    double bin_width = ifo.freqs[right_idx] - ifo.freqs[left_idx];
    CPL r0 = 0.5 * (r_left + r_right);
    CPL r1 = (bin_width > 0.) ? (r_right - r_left) / bin_width : CPL(0.);

    double approx_dh = real(A0 * conj(r0) + A1 * conj(r1));
    double approx_hh =
        real(B0) * std::norm(r0) + 2. * real(r0 * conj(r1)) * real(B1);

    exact_logL += -0.5 * exact_hh + exact_dh;
    approx_logL += -0.5 * approx_hh + approx_dh;
  }

  return std::abs(exact_logL - approx_logL);
}

int RelativeBinningBisectionLikelihood::bin_bisection(
    const std::vector<IfoData>& ifos_data,
    const std::vector<VECCPL>& fiducial_data,
    const std::vector<VECCPL>& test_data, const double* frequencies,
    double epsilon) {
  int max_idx = find_max_index(fiducial_data);
  RelativeBinningPrinter("Max frequency: " +
                         std::to_string(frequencies[max_idx]));

  // Collect sorted bin boundary indices via iterative DFS bisection.
  // The set always contains 0 (left edge of the first bin); each accepted
  // interval [left, right) contributes its right index.
  std::set<int> boundaries;
  boundaries.insert(0);

  std::vector<std::pair<int, int>> stack;
  stack.push_back(std::make_pair(0, max_idx));

  while (!stack.empty()) {
    int left = stack.back().first;
    int right = stack.back().second;
    stack.pop_back();

    if (right - left <= 1) {
      boundaries.insert(right);
      continue;
    }

    double err = bin_log_likelihood_error(left, right, ifos_data, fiducial_data,
                                          test_data);

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
    bin_sizes.push_back(bin_inds[i] - bin_inds[i - 1]);
    bin_widths.push_back(bin_freqs[i] - bin_freqs[i - 1]);
    bin_centers.push_back((bin_freqs[i] + bin_freqs[i - 1]) * 0.5);
  }

  return 1;
}

void RelativeBinningBisectionLikelihood::setup_summary_data(
    const std::vector<IfoData>& ifos_data,
    const std::vector<VECCPL>& fiducial_data) {
  int num_det = (int)ifos_data.size();
  for (int det = 0; det < num_det; det++) {
    SummaryData det_summary_data;

    const IfoData& ifo = ifos_data[det];
    const VECCPL& h0 = fiducial_data[det];
    double weight = 4. / ifo.duration;

    for (int b = 0; b < number_of_bins; b++) {
      int start = bin_inds[b];
      int end = bin_inds[b + 1];

      det_summary_data.strain.push_back(h0[start]);

      CPL A0 = 0., A1 = 0., B0 = 0., B1 = 0.;
      for (int j = start; j < end; j++) {
        double delta = ifo.freqs[j] - bin_centers[b];
        double psd = ifo.psd[j];
        CPL d = ifo.strain[j];

        CPL A_fac = conj(h0[j]) / psd;
        CPL B_fac = h0[j] * A_fac;
        A_fac *= d;

        A0 += A_fac;
        A1 += A_fac * delta;
        B0 += B_fac;
        B1 += B_fac * delta;
      }
      det_summary_data.A0.push_back(weight * A0);
      det_summary_data.A1.push_back(weight * A1);
      det_summary_data.B0.push_back(weight * B0);
      det_summary_data.B1.push_back(weight * B1);
    }
    det_summary_data.strain.push_back(h0[bin_inds.back()]);

    ifos_summary_data.push_back(det_summary_data);
  }
}

void RelativeBinningBisectionLikelihood::compute_waveform_ratios(
    VECCPL& r0, VECCPL& r1, const CPL* h, const SummaryData& fiducial) {
  r0.clear(); r1.clear();

  CPL ratio_left = h[0] / fiducial.strain.front();
  CPL ratio_right;

  for (int i = 0; i < number_of_bins; i++) {
    ratio_right = h[i + 1] / fiducial.strain[i + 1];

    r0.push_back(0.5 * (ratio_right + ratio_left));
    r1.push_back((ratio_right - ratio_left) / bin_widths[i]);

    ratio_left = ratio_right;
  }
}

double RelativeBinningBisectionLikelihood::log_likelihood_per_detector(
    const CPL* h, const SummaryData& fiducial) {
  double d_h = 0., h_h = 0.;

  VECCPL r0, r1;
  compute_waveform_ratios(r0, r1, h, fiducial);

  for (int b = 0; b < number_of_bins; b++) {
    d_h += real(fiducial.A0[b] * conj(r0[b]) + fiducial.A1[b] * conj(r1[b]));
    h_h += real(fiducial.B0[b] * r0[b] * conj(r0[b]) +
                2. * fiducial.B1[b] * real(r0[b] * conj(r1[b])));
  }

  return -0.5 * h_h + d_h;
}

double RelativeBinningBisectionLikelihood::log_likelihood_from_waveform(
    const std::vector<VECCPL>& h_at_bins) {
  int num_det = (int)ifos_summary_data.size();
  double logL = 0.;
  for (int i = 0; i < num_det; i++)
    logL +=
        log_likelihood_per_detector(h_at_bins[i].data(), ifos_summary_data[i]);
  return logL;
}

double RelativeBinningBisectionLikelihood::log_likelihood(
    std::string* detectors, int num_detectors, gen_params_base<double>* params,
    std::string generation_method, bool reuse_WF) {
  int* data_lengths = new int[num_detectors];
  CPL** responses = new CPL*[num_detectors];
  double** frequencies = new double*[num_detectors];

  for (int i = 0; i < num_detectors; i++) {
    data_lengths[i] = (int)bin_freqs.size();
    responses[i] = new CPL[data_lengths[i]];
    frequencies[i] = bin_freqs.data();
  }

  if (num_detectors == 1) {
    create_single_GW_detection(responses[0], detectors[0], frequencies[0],
                               data_lengths[0], params, generation_method);
  } else {
    create_coherent_GW_detection(detectors, num_detectors, frequencies,
                                 data_lengths, reuse_WF, params,
                                 generation_method, responses);
  }

  double logL = 0.;
  for (int i = 0; i < num_detectors; i++)
    logL += log_likelihood_per_detector(responses[i], ifos_summary_data[i]);

  for (int i = 0; i < num_detectors; i++) delete[] responses[i];
  delete[] frequencies;
  delete[] responses;
  delete[] data_lengths;

  return logL;
}

RelativeBinningBisectionLikelihood RelBinBisectLL_from_c_arrays(
    double** freq, CPL** data, double** psds, CPL** fiducials, CPL** test_wfs,
    int length, int n_detect, double duration, double epsilon) {
  std::vector<IfoData> ifos_data;
  std::vector<VECCPL> fiducial_vecs;
  std::vector<VECCPL> test_vecs;

  for (int i = 0; i < n_detect; i++) {
    IfoData ifo;
    ifo.duration = duration;
    ifo.freqs = VECDBL(freq[i], freq[i] + length);
    ifo.psd = VECDBL(psds[i], psds[i] + length);
    ifo.strain = VECCPL(data[i], data[i] + length);
    ifos_data.push_back(ifo);

    fiducial_vecs.push_back(VECCPL(fiducials[i], fiducials[i] + length));
    test_vecs.push_back(VECCPL(test_wfs[i], test_wfs[i] + length));
  }

  return RelativeBinningBisectionLikelihood(ifos_data, fiducial_vecs, test_vecs,
                                            freq[0], epsilon);
}

}  // namespace GWATLikelihoods
