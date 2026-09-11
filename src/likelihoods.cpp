#include "likelihoods.h"

#include <assert.h>
#include <waveform_util.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <set>
#include <utility>

#include "detector_util.h"

namespace gw_likelihoods {

double log_likelihood_in_detector(const VECCPL& h, const IfoData& ifo,
                                  const Quadrature& quad) {
  double hh = quad.inner_product(h, h, ifo.psd);
  double dh = quad.inner_product(ifo.strain, h, ifo.psd);
  return -0.5 * hh + dh;
}

double log_likelihood_in_detector(const CPL* h, const CPL* d, const VECDBL& psd,
                                  int length, const Quadrature& quad) {
  const VECCPL hvec(h, h + length);
  const VECCPL dvec(d, d + length);
  double hh = quad.inner_product(hvec, hvec, psd);
  double dh = quad.inner_product(dvec, hvec, psd);
  return -0.5 * hh + dh;
};

double log_likelihood_in_detector(const VECCPL& h, const VECCPL& d,
                                  const VECDBL& psd, const Quadrature& quad) {
  double hh = quad.inner_product(h, h, psd);
  double dh = quad.inner_product(d, h, psd);
  return -0.5 * hh + dh;
}

double CoherentBareLikelihood::log_likelihood(
    gen_params_base<double>* params) const {
  gen_params_base<double> local_params = *params;
  local_params.f_ref = f_ref_;
  local_params.gmst = gmst_;
  local_params.equatorial_orientation = equatorial_orientation_;
  local_params.horizon_coord = horizon_coord_;
  local_params.shift_time = shift_time_;
  local_params.shift_phase = shift_phase_;

  const IfoData& reference_ifo = ifos_[0];
  waveform_polarizations<double> wf;
  wf.allocate_memory(reference_ifo.freqs.size());
  waveform_generator_.fill_polarizations(&wf, &local_params, reference_ifo.freqs);

  double ll = 0.0;
  for (auto& ifo : ifos_) {
    int data_length = ifo.freqs.size();
    VECCPL response(data_length);
    fourier_detector_response_equatorial(
        const_cast<double*>(ifo.freqs.data()), data_length, &wf,
        response.data(), local_params.RA, local_params.DEC, local_params.psi,
        local_params.gmst, (double*)nullptr, local_params.LISA_alpha0,
        local_params.LISA_phi0, local_params.theta_l, local_params.phi_l,
        ifo.name);

    if (ifo.name != reference_ifo.name) {
      double time_shift =
          -2.0 * kPi *
          DTOA_DETECTOR(local_params.RA, local_params.DEC, local_params.gmst,
                        reference_ifo.name, ifo.name);
      for (int j = 0; j < data_length; j++)
        response[j] *= std::exp(CPL(0, time_shift * ifo.freqs[j]));
    }
    ll += log_likelihood_in_detector(response, ifo, quad);
  }
  wf.deallocate_memory();
  return ll;
}

std::vector<VECCPL> PolarizationsLikelihood::generate_modes(
    const double* theta, const VECDBL& freqs) const {
  gen_params_base<double> gp;
  pmap_.to_gen_params(theta, gp);
  gp.f_ref = f_ref_;
  gp.gmst = gmst_;
  gp.shift_time = shift_time_;
  gp.shift_phase = shift_phase_;
  return waveform_generator_.generate_polarizations(&gp, freqs);
}

double PolarizationsLikelihood::log_likelihood(
    gen_params_base<double>* params) const {
  gen_params_base<double> local_params = *params;
  local_params.f_ref = f_ref_;
  local_params.gmst = gmst_;
  local_params.shift_time = shift_time_;
  local_params.shift_phase = shift_phase_;

  auto modes = waveform_generator_.generate_polarizations(&local_params, data_.freqs);

  double ll = 0.0;
  for (int m = 0; m < static_cast<int>(modes.size()); ++m)
    ll += sky_avg_factors_[m] *
          log_likelihood_in_detector(modes[m], data_.modes[m], data_.psd, quad);
  return ll;
}

namespace RelativeBinning {

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

  duration = ifos_data[0].freqs[1] - ifos_data[0].freqs[0];
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
  double maxFreq = frequencies[data_length - 1];

  for (int d = 0; d < num_detectors; d++) {
    const CPL* data = fiducial_data_in[d];

    int idx = data_length - 1;
    while (idx >= 0 && data[idx] == CPL(0.)) {
      idx--;
    }

    if (idx >= 0) {
      maxFreq = std::min(maxFreq, frequencies[idx]);
    }
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
  CPL A_fac, B_fac, A0_b, A1_b;
  double B0_b, B1_b;
  double delta_f, inner_product_weight;

  // Compute summary data for each interferometer
  for (int i = 0; i < ifos_data.size(); i++) {
    // The ifo holding the data
    ifo = &(ifos_data.at(i));
    // The fiducial data
    fiducial = fiducial_data[i];
    // The ifo to store the summary data
    ifo_fiducial = &ifos_fiducial_data.at(i);

    inner_product_weight = 4. * (ifo->freqs[1] - ifo->freqs[0]);

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
      A0_b = A1_b = 0.; B0_b = B1_b = 0.;
      for (int f = 0; f < freqs.size(); f++) {
        A_fac = conj(h0[f]) / psd[f];
        B_fac = h0[f] * A_fac;
        A_fac *= data[f];
        delta_f = freqs[f] - bin_centers[b];

        A0_b += A_fac;
        A1_b += A_fac * delta_f;
        B0_b += std::real(B_fac);
        B1_b += std::real(B_fac) * delta_f;
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
    h_h += fiducial->B0[b] * std::norm(r0[b]) +
           2. * fiducial->B1[b] * real(r0[b] * conj(r1[b]));
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
// Shared helper
// ============================================================

/// @brief Compute the summary data of a bin with the trapezoidal rule.
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

// ============================================================
// RelativeBinningBisectionLikelihood
// ============================================================

RelativeBinningBisectionLikelihood::RelativeBinningBisectionLikelihood(
    const IfoData& ifo, const VECCPL& fiducial_data, const VECCPL& test_data,
    const WaveformGenerator& waveform_generator, double epsilon, double f_ref,
    double gmst, bool equatorial_orientation, bool horizon_coord,
    bool shift_time, bool shift_phase, bool log_spacing)
    : waveform_generator_(waveform_generator),
      f_ref_(f_ref),
      gmst_(gmst),
      equatorial_orientation_(equatorial_orientation),
      horizon_coord_(horizon_coord),
      shift_time_(shift_time),
      shift_phase_(shift_phase) {
  std::cout << "\nRELATIVE BINNING (BISECTION) INITIALIZING\n";

  bin_bisection(ifo, fiducial_data, test_data, epsilon, log_spacing);
  RelativeBinningPrinter(std::to_string(number_of_bins) + " bins setup");

  setup_summary_data(ifo, fiducial_data, log_spacing);

  int num_edges = (int)bin_inds.size();
  VECCPL ht(num_edges);
  for (int k = 0; k < num_edges; k++) ht[k] = test_data[bin_inds[k]];

  double logL = log_likelihood_at_waveform(ht);
  RelativeBinningPrinter("logL of test data: " + std::to_string(logL));
}

std::pair<int, int> RelativeBinningBisectionLikelihood::find_min_max_indices(
    const VECCPL& strain) {
  int min_idx = 0;
  int max_idx = static_cast<int>(strain.size()) - 1;

  // Find first non-zero: some waveforms (e.g. EFPE) can be exactly zero at
  // f = f_start due to a floating-point boundary in the SPA frequency check.
  auto first_nz = std::find_if(strain.begin(), strain.end(),
                               [](const CPL& v) { return std::norm(v) > 0.; });
  if (first_nz == strain.end())
    throw std::runtime_error(
        "RELATIVE BINNING: Fiducial strain entirely zero, unable to bin.");
  else
    min_idx = std::distance(strain.begin(), first_nz);

  // Find the first zero after the waveform has become non-zero (merger cutoff).
  auto it = std::find_if(first_nz, strain.end(),
                         [](const CPL& v) { return std::norm(v) == 0.; });
  if (it != strain.end()) {
    int zero_idx = std::distance(strain.begin(), it);
    max_idx = std::min(max_idx, zero_idx - 1);
  }

  // Require at least index 1 so bin_freqs always has at least two distinct
  // edges.
  return std::pair<int, int>(min_idx, std::max(max_idx, 1));
}

double RelativeBinningBisectionLikelihood::bin_log_likelihood_error(
    int left_idx, int right_idx, const IfoData& ifo, const VECCPL& h0,
    const VECCPL& ht, bool log_spacing) {
  double weight =
      (log_spacing
           ? LOG10 * (std::log10(ifo.freqs[1]) - std::log10(ifo.freqs[0]))
           : (ifo.freqs[1] - ifo.freqs[0])) *
      4.0;
  double f_m = (log_spacing ? (std::log10(ifo.freqs[left_idx]) +
                               std::log10(ifo.freqs[right_idx]))
                            : (ifo.freqs[left_idx] + ifo.freqs[right_idx])) *
               0.5;

  CPL A0 = 0., A1 = 0.;
  double B0 = 0., B1 = 0.;
  compute_summary_data_per_bin_trapezoid_rule(A0, A1, B0, B1, f_m, left_idx,
                                              right_idx, ifo.freqs, ifo.psd, h0,
                                              ifo.strain, log_spacing, weight);

  double exact_dh = 0., exact_hh = 0.;
  for (int j = left_idx + 1; j < right_idx; j++) {
    exact_dh += log_spacing ? real(ifo.strain[j] * conj(ht[j]) * ifo.freqs[j]) /
                                  ifo.psd[j]
                            : real(ifo.strain[j] * conj(ht[j])) / ifo.psd[j];
    exact_hh += log_spacing ? std::norm(ht[j]) * ifo.freqs[j] / ifo.psd[j]
                            : std::norm(ht[j]) / ifo.psd[j];
  }
  for (int j : {left_idx, right_idx}) {
    exact_dh +=
        0.5 * (log_spacing ? real(ifo.strain[j] * conj(ht[j]) * ifo.freqs[j]) /
                                 ifo.psd[j]
                           : real(ifo.strain[j] * conj(ht[j])) / ifo.psd[j]);
    exact_hh +=
        0.5 * (log_spacing ? std::norm(ht[j]) * ifo.freqs[j] / ifo.psd[j]
                           : std::norm(ht[j]) / ifo.psd[j]);
  }
  exact_dh *= weight;
  exact_hh *= weight;

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
  double approx_hh = B0 * std::norm(r0) + 2. * B1 * real(r0 * conj(r1));

  double exact_logL = -0.5 * exact_hh + exact_dh;
  double approx_logL = -0.5 * approx_hh + approx_dh;

  return std::abs(exact_logL - approx_logL);
}

void RelativeBinningBisectionLikelihood::bin_bisection(
    const IfoData& ifo, const VECCPL& test_data, const VECCPL& fiducial_data,
    double epsilon, bool log_spacing) {
  auto min_max_idx = find_min_max_indices(ifo.strain);
  int min_idx = min_max_idx.first, max_idx = min_max_idx.second;
  RelativeBinningPrinter("Max frequency: " +
                         std::to_string(ifo.freqs[max_idx]));

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

    double err = bin_log_likelihood_error(left, right, ifo, fiducial_data,
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
    bin_freqs.push_back(ifo.freqs[idx]);
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

void RelativeBinningBisectionLikelihood::setup_summary_data(const IfoData& ifo,
                                                            const VECCPL& h0,
                                                            bool log_spacing) {
  ifo_ = ifo;  // copy all fields; clear the data arrays below
  ifo_.freqs.clear();
  ifo_.psd.clear();
  ifo_.strain.clear();
  double weight =
      (log_spacing
           ? LOG10 * (std::log10(ifo.freqs[1]) - std::log10(ifo.freqs[0]))
           : (ifo.freqs[1] - ifo.freqs[0])) *
      4.0;

  for (int b = 0; b < number_of_bins; b++) {
    int start = bin_inds[b];
    int end = bin_inds[b + 1];

    ifo_.strain.push_back(h0[start]);
    ifo_.psd.push_back(ifo.psd[start]);
    ifo_.freqs.push_back(ifo.freqs[start]);

    CPL A0 = 0., A1 = 0.;
    double B0 = 0., B1 = 0.;
    compute_summary_data_per_bin_trapezoid_rule(
        A0, A1, B0, B1, bin_centers[b], start, end, ifo.freqs, ifo.psd, h0,
        ifo.strain, log_spacing, weight);
    ifo_summary_data_.A0.push_back(A0);
    ifo_summary_data_.A1.push_back(A1);
    ifo_summary_data_.B0.push_back(B0);
    ifo_summary_data_.B1.push_back(B1);
  }
  ifo_.strain.push_back(h0[bin_inds.back()]);
  ifo_.psd.push_back(ifo.psd[bin_inds.back()]);
  ifo_.freqs.push_back(ifo.freqs[bin_inds.back()]);
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

double RelativeBinningBisectionLikelihood::log_likelihood_at_waveform(
    const VECCPL& h) const {
  double d_h = 0., h_h = 0.;
  VECCPL r0, r1;
  compute_waveform_ratios(r0, r1, h, ifo_.strain);

  for (int b = 0; b < number_of_bins; b++) {
    d_h += std::real(ifo_summary_data_.A0[b] * conj(r0[b]) +
                     ifo_summary_data_.A1[b] * conj(r1[b]));
    h_h += ifo_summary_data_.B0[b] * std::norm(r0[b]) +
           2. * ifo_summary_data_.B1[b] * real(r0[b] * conj(r1[b]));
  }

  return -0.5 * h_h + d_h;
}

double RelativeBinningBisectionLikelihood::log_likelihood(
    gen_params_base<double>* params) const {
  gen_params_base<double> local_params = *params;
  local_params.f_ref = f_ref_;
  local_params.gmst = gmst_;
  local_params.equatorial_orientation = equatorial_orientation_;
  local_params.horizon_coord = horizon_coord_;
  local_params.shift_time = shift_time_;
  local_params.shift_phase = shift_phase_;

  waveform_polarizations<double> wf;
  wf.allocate_memory(ifo_.freqs.size());
  waveform_generator_.fill_polarizations(&wf, &local_params, ifo_.freqs);

  int data_length = ifo_.freqs.size();
  VECCPL response(data_length);
  fourier_detector_response_equatorial(
      const_cast<double*>(ifo_.freqs.data()), data_length, &wf, response.data(),
      local_params.RA, local_params.DEC, local_params.psi, local_params.gmst,
      (double*)nullptr, local_params.LISA_alpha0, local_params.LISA_phi0,
      local_params.theta_l, local_params.phi_l, ifo_.name);

  wf.deallocate_memory();
  return log_likelihood_at_waveform(response);
}

// ============================================================
// RelativeBinningBisectionPolarizationsLikelihood
// ============================================================

RelativeBinningBisectionPolarizationsLikelihood::
    RelativeBinningBisectionPolarizationsLikelihood(
        const ParameterMap& pmap,
        const PolarizationData& data,
        const std::vector<VECCPL>& fiducial_modes,
        const std::vector<VECCPL>& test_modes,
        const std::vector<double>& sky_avg_factors,
        const WaveformGenerator& waveform_generator,
        double epsilon, double f_ref, double gmst,
        bool shift_time, bool shift_phase, bool log_spacing)
    : pmap_(pmap),
      waveform_generator_(waveform_generator),
      f_ref_(f_ref),
      gmst_(gmst),
      shift_time_(shift_time),
      shift_phase_(shift_phase),
      sky_avg_factors_(sky_avg_factors) {
  number_of_modes_ = static_cast<int>(data.modes.size());

  if (static_cast<int>(fiducial_modes.size()) != number_of_modes_ ||
      static_cast<int>(test_modes.size()) != number_of_modes_ ||
      static_cast<int>(sky_avg_factors.size()) != number_of_modes_)
    throw std::invalid_argument(
        "RELATIVE BINNING (POLARIZATIONS): fiducial_modes, test_modes, "
        "sky_avg_factors, and data.modes must all have the same length.");

  std::cout << "\nRELATIVE BINNING (BISECTION, POLARIZATIONS) INITIALIZING\n";

  bin_bisection(data, fiducial_modes, test_modes, epsilon, log_spacing);
  RelativeBinningPrinter(std::to_string(number_of_bins_) + " bins setup");

  setup_summary_data(data, fiducial_modes, log_spacing);

  int n_edges = static_cast<int>(bin_inds_.size());
  std::vector<VECCPL> ht(number_of_modes_, VECCPL(n_edges));
  std::vector<VECCPL> d(number_of_modes_, VECCPL(n_edges));
  for (int m = 0; m < number_of_modes_; ++m)
    for (int k = 0; k < n_edges; ++k){
      ht[m][k] = test_modes[m][bin_inds_[k]];
      d[m][k] = data.modes[m][bin_inds_[k]];
    }

  double logL = log_likelihood_at_waveform(ht);
  RelativeBinningPrinter("logL of test data: " + std::to_string(logL));
  double snr = std::sqrt(2.0  * std::max(0.0, log_likelihood_at_waveform(d)));
  RelativeBinningPrinter("SNR of data: " + std::to_string(snr));
}

std::pair<int, int>
RelativeBinningBisectionPolarizationsLikelihood::find_min_max_indices(
    const std::vector<VECCPL>& modes) {
  int min_idx = 0;
  int max_idx = static_cast<int>(modes[0].size()) - 1;

  for (const auto& mode : modes) {
    auto first_nz = std::find_if(mode.begin(), mode.end(),
                                 [](const CPL& v) { return std::norm(v) > 0.; });
    if (first_nz == mode.end())
      throw std::runtime_error(
          "RELATIVE BINNING (POLARIZATIONS): A polarization mode is entirely "
          "zero, unable to bin.");
    min_idx = std::max(min_idx,
                       static_cast<int>(std::distance(mode.begin(), first_nz)));

    auto it = std::find_if(first_nz, mode.end(),
                           [](const CPL& v) { return std::norm(v) == 0.; });
    if (it != mode.end())
      max_idx = std::min(max_idx,
                         static_cast<int>(std::distance(mode.begin(), it)) - 1);
  }

  return {min_idx, std::max(max_idx, 1)};
}

double
RelativeBinningBisectionPolarizationsLikelihood::bin_log_likelihood_error(
    int left_idx, int right_idx, const PolarizationData& data,
    const std::vector<VECCPL>& fiducial_modes,
    const std::vector<VECCPL>& test_modes, bool log_spacing) {
  double weight =
      (log_spacing
           ? LOG10 * (std::log10(data.freqs[1]) - std::log10(data.freqs[0]))
           : (data.freqs[1] - data.freqs[0])) *
      4.0;
  double f_m =
      (log_spacing
           ? std::log10(data.freqs[left_idx]) + std::log10(data.freqs[right_idx])
           : data.freqs[left_idx] + data.freqs[right_idx]) *
      0.5;
  double bin_width =
      log_spacing
          ? std::log10(data.freqs[right_idx]) - std::log10(data.freqs[left_idx])
          : data.freqs[right_idx] - data.freqs[left_idx];

  double total_exact_ll = 0., total_approx_ll = 0.;

  for (int m = 0; m < number_of_modes_; ++m) {
    CPL A0 = 0., A1 = 0.;
    double B0 = 0., B1 = 0.;
    compute_summary_data_per_bin_trapezoid_rule(
        A0, A1, B0, B1, f_m, left_idx, right_idx, data.freqs, data.psd,
        fiducial_modes[m], data.modes[m], log_spacing, weight);

    double exact_dh = 0., exact_hh = 0.;
    for (int j = left_idx + 1; j < right_idx; ++j) {
      double fac = (log_spacing ? data.freqs[j] : 1.0) / data.psd[j];
      exact_dh += std::real(data.modes[m][j] * std::conj(test_modes[m][j])) * fac;
      exact_hh += std::norm(test_modes[m][j]) * fac;
    }
    for (int j : {left_idx, right_idx}) {
      double fac = 0.5 * (log_spacing ? data.freqs[j] : 1.0) / data.psd[j];
      exact_dh += std::real(data.modes[m][j] * std::conj(test_modes[m][j])) * fac;
      exact_hh += std::norm(test_modes[m][j]) * fac;
    }
    exact_dh *= weight;
    exact_hh *= weight;

    CPL r_left = (fiducial_modes[m][left_idx] != CPL(0.))
                     ? test_modes[m][left_idx] / fiducial_modes[m][left_idx]
                     : CPL(0.);
    CPL r_right = (fiducial_modes[m][right_idx] != CPL(0.))
                      ? test_modes[m][right_idx] / fiducial_modes[m][right_idx]
                      : CPL(0.);
    CPL r0 = 0.5 * (r_left + r_right);
    CPL r1 = (bin_width > 0.) ? (r_right - r_left) / bin_width : CPL(0.);

    double approx_dh = std::real(A0 * std::conj(r0) + A1 * std::conj(r1));
    double approx_hh =
        B0 * std::norm(r0) + 2. * B1 * std::real(r0 * std::conj(r1));

    total_exact_ll  += sky_avg_factors_[m] * (-0.5 * exact_hh  + exact_dh);
    total_approx_ll += sky_avg_factors_[m] * (-0.5 * approx_hh + approx_dh);
  }

  return std::abs(total_exact_ll - total_approx_ll);
}

void RelativeBinningBisectionPolarizationsLikelihood::bin_bisection(
    const PolarizationData& data, const std::vector<VECCPL>& fiducial_modes,
    const std::vector<VECCPL>& test_modes, double epsilon, bool log_spacing) {
  auto min_max_idx = find_min_max_indices(fiducial_modes);
  int min_idx = min_max_idx.first, max_idx = min_max_idx.second;
  RelativeBinningPrinter("Max frequency: " +
                         std::to_string(data.freqs[max_idx]));

  std::set<int> boundaries;
  boundaries.insert(min_idx);

  std::vector<std::pair<int, int>> stack;
  stack.push_back({min_idx, max_idx});

  while (!stack.empty()) {
    int left = stack.back().first, right = stack.back().second;
    stack.pop_back();

    if (right - left <= 1) {
      boundaries.insert(right);
      continue;
    }

    double err = bin_log_likelihood_error(left, right, data, fiducial_modes,
                                          test_modes, log_spacing);
    if (err <= epsilon) {
      boundaries.insert(right);
    } else {
      int mid = (left + right) / 2;
      stack.push_back({left, mid});
      stack.push_back({mid, right});
    }
  }

  for (int idx : boundaries) {
    bin_inds_.push_back(idx);
    bin_freqs_.push_back(data.freqs[idx]);
  }

  number_of_bins_ = static_cast<int>(bin_inds_.size()) - 1;
  for (int i = 1; i < static_cast<int>(bin_inds_.size()); ++i) {
    if (log_spacing) {
      bin_widths_.push_back(std::log10(bin_freqs_[i]) -
                            std::log10(bin_freqs_[i - 1]));
      bin_centers_.push_back(
          (std::log10(bin_freqs_[i]) + std::log10(bin_freqs_[i - 1])) * 0.5);
    } else {
      bin_widths_.push_back(bin_freqs_[i] - bin_freqs_[i - 1]);
      bin_centers_.push_back((bin_freqs_[i] + bin_freqs_[i - 1]) * 0.5);
    }
  }
}

void RelativeBinningBisectionPolarizationsLikelihood::setup_summary_data(
    const PolarizationData& data, const std::vector<VECCPL>& fiducial_modes,
    bool log_spacing) {
  double weight =
      (log_spacing
           ? LOG10 * (std::log10(data.freqs[1]) - std::log10(data.freqs[0]))
           : (data.freqs[1] - data.freqs[0])) *
      4.0;

  summary_data_.resize(number_of_modes_);
  fiducial_at_edges_.resize(number_of_modes_);

  for (int b = 0; b < number_of_bins_; ++b) {
    int start = bin_inds_[b];
    int end = bin_inds_[b + 1];

    binned_data_.freqs.push_back(data.freqs[start]);
    for (int m = 0; m < number_of_modes_; ++m)
      fiducial_at_edges_[m].push_back(fiducial_modes[m][start]);

    for (int m = 0; m < number_of_modes_; ++m) {
      CPL A0 = 0., A1 = 0.;
      double B0 = 0., B1 = 0.;
      compute_summary_data_per_bin_trapezoid_rule(
          A0, A1, B0, B1, bin_centers_[b], start, end, data.freqs, data.psd,
          fiducial_modes[m], data.modes[m], log_spacing, weight);
      summary_data_[m].A0.push_back(A0);
      summary_data_[m].A1.push_back(A1);
      summary_data_[m].B0.push_back(B0);
      summary_data_[m].B1.push_back(B1);
    }
  }

  // Add the final right edge.
  binned_data_.freqs.push_back(data.freqs[bin_inds_.back()]);
  for (int m = 0; m < number_of_modes_; ++m)
    fiducial_at_edges_[m].push_back(fiducial_modes[m][bin_inds_.back()]);
}

void RelativeBinningBisectionPolarizationsLikelihood::compute_waveform_ratios(
    std::vector<VECCPL>& r0, std::vector<VECCPL>& r1,
    const std::vector<VECCPL>& h) const {
  r0.assign(number_of_modes_, VECCPL());
  r1.assign(number_of_modes_, VECCPL());

  for (int m = 0; m < number_of_modes_; ++m) {
    CPL ratio_left = (fiducial_at_edges_[m][0] != CPL(0.))
                         ? h[m][0] / fiducial_at_edges_[m][0]
                         : CPL(0.);
    for (int b = 0; b < number_of_bins_; ++b) {
      CPL ratio_right = (fiducial_at_edges_[m][b + 1] != CPL(0.))
                            ? h[m][b + 1] / fiducial_at_edges_[m][b + 1]
                            : CPL(0.);
      r0[m].push_back(0.5 * (ratio_left + ratio_right));
      r1[m].push_back((ratio_right - ratio_left) / bin_widths_[b]);
      ratio_left = ratio_right;
    }
  }
}

double
RelativeBinningBisectionPolarizationsLikelihood::log_likelihood_at_waveform(
    const std::vector<VECCPL>& h) const {
  std::vector<VECCPL> r0, r1;
  compute_waveform_ratios(r0, r1, h);

  double ll = 0.0;
  for (int m = 0; m < number_of_modes_; ++m) {
    double d_h = 0., h_h = 0.;
    for (int b = 0; b < number_of_bins_; ++b) {
      d_h += std::real(summary_data_[m].A0[b] * std::conj(r0[m][b]) +
                       summary_data_[m].A1[b] * std::conj(r1[m][b]));
      h_h += summary_data_[m].B0[b] * std::norm(r0[m][b]) +
             2. * summary_data_[m].B1[b] *
                 std::real(r0[m][b] * std::conj(r1[m][b]));
    }
    ll += sky_avg_factors_[m] * (-0.5 * h_h + d_h);
  }
  return ll;
}

std::vector<VECCPL> RelativeBinningBisectionPolarizationsLikelihood::generate_modes(
    const double* theta, const VECDBL& freqs) const {
  gen_params_base<double> gp;
  pmap_.to_gen_params(theta, gp);
  gp.f_ref = f_ref_;
  gp.gmst = gmst_;
  gp.shift_time = shift_time_;
  gp.shift_phase = shift_phase_;
  return waveform_generator_.generate_polarizations(&gp, freqs);
}

double RelativeBinningBisectionPolarizationsLikelihood::log_likelihood(
    gen_params_base<double>* params) const {
  gen_params_base<double> local_params = *params;
  local_params.f_ref = f_ref_;
  local_params.gmst = gmst_;
  local_params.shift_time = shift_time_;
  local_params.shift_phase = shift_phase_;

  auto h_at_bins = waveform_generator_.generate_polarizations(&local_params,
                                                              binned_data_.freqs);

  if (static_cast<int>(h_at_bins.size()) != number_of_modes_)
    throw std::runtime_error(
        "RELATIVE BINNING (POLARIZATIONS): waveform generator produced " +
        std::to_string(h_at_bins.size()) + " active modes but expected " +
        std::to_string(number_of_modes_));

  return log_likelihood_at_waveform(h_at_bins);
}

}  // namespace RelativeBinning

}  // namespace gw_likelihoods
