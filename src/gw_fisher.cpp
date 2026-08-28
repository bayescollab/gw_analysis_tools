#include "gw_fisher.h"

#include <memory>

#include "waveform_util.h"

namespace gw_fishers {
std::vector<std::vector<VECCPL>>
PolarizationDataFishers::compute_derivative_in_detector(
    const double* theta, const PolarizationData& ifo) const {
  const int dim = pmap_.dim();
  const int len = static_cast<int>(ifo.freqs.size());

  auto make_wf = [&](int i, double val) -> std::vector<VECCPL> {
    std::vector<double> th(theta, theta + dim);
    th[i] = val;
    return likelihood_.generate_modes(th.data(), ifo.freqs);
  };

  const auto& factors = likelihood_.sky_avg_factors();
  const int n_modes = static_cast<int>(factors.size());
  std::vector<std::vector<VECCPL>> modes_derivs(
      dim, std::vector<VECCPL>(n_modes, VECCPL(len)));
  for (int i = 0; i < dim; i++) {
    const bool clamped =
        (i == eta_idx_ && eta_idx_ >= 0 && theta[i] > 0.25 - eps_);

    if (order_ == 4) {
      const double vpp = clamped ? theta[i] : theta[i] + 2.0 * eps_;
      const double vp = clamped ? theta[i] : theta[i] + eps_;
      const double vm = theta[i] - eps_;
      const double vmm = theta[i] - 2.0 * eps_;

      std::vector<VECCPL> h_pp = make_wf(i, vpp);
      std::vector<VECCPL> h_p = make_wf(i, vp);
      std::vector<VECCPL> h_m = make_wf(i, vm);
      std::vector<VECCPL> h_mm = make_wf(i, vmm);

      const double denom = clamped ? 6.0 * eps_ : 12.0 * eps_;
      for (int m = 0; m < n_modes; m++) {
        for (int l = 0; l < len; l++)
          modes_derivs[i][m][l] =
              (-h_pp[m][l] + 8.0 * h_p[m][l] - 8.0 * h_m[m][l] + h_mm[m][l]) /
              denom;
      }
    } else {
      const double vp = clamped ? theta[i] : theta[i] + eps_;
      const double vm = theta[i] - eps_;

      std::vector<VECCPL> h_p = make_wf(i, vp);
      std::vector<VECCPL> h_m = make_wf(i, vm);

      const double denom = clamped ? eps_ : 2.0 * eps_;
      for (int m = 0; m < n_modes; m++) {
        for (int l = 0; l < len; l++)
          modes_derivs[i][m][l] = (h_p[m][l] - h_m[m][l]) / denom;
      }
    }
  }
  return modes_derivs;
}

void PolarizationDataFishers::compute_Fisher(double** fisherM,
                                             const double* theta) const {
  const int dim = pmap_.dim();

  for (int i = 0; i < dim; i++)
    for (int j = 0; j < dim; j++) fisherM[i][j] = 0.0;

  const auto& factors = likelihood_.sky_avg_factors();
  for (const auto& ifo : ifos_) {
    auto resp_deriv = compute_derivative_in_detector(theta, ifo);
    for (int j = 0; j < dim; j++)
      for (int k = 0; k <= j; k++)
        for (int m = 0; m < static_cast<int>(factors.size()); m++)
          fisherM[j][k] +=
              factors[m] *
              quad_.inner_product(resp_deriv[j][m], resp_deriv[k][m], ifo.psd);
  }

  // Symmetrize
  for (int j = 0; j < dim; j++)
    for (int k = 0; k < j; k++) fisherM[k][j] = fisherM[j][k];
}

void RelativeBinningPolarizationsFishers::compute_Fisher(
    double** fisherM, const double* theta) const {
  const int dim = pmap_.dim();
  const auto& factors = rb_.sky_avg_factors();
  const int n_modes = static_cast<int>(factors.size());
  const VECDBL& bin_freqs = rb_.get_bin_freqs();
  const auto& summary = rb_.template_summary();
  const int n_bins = static_cast<int>(summary[0].B0.size());

  for (int i = 0; i < dim; i++)
    for (int j = 0; j < dim; j++) fisherM[i][j] = 0.0;

  // Build normalised derivative ratios g0[i][m][b] and g1[i][m][b]
  // for every sampled parameter i at every bin edge b.
  std::vector<std::vector<VECCPL>> g0(dim), g1(dim);

  auto make_modes = [&](int param, double val) -> std::vector<VECCPL> {
    std::vector<double> th(theta, theta + dim);
    th[param] = val;
    return rb_.generate_modes(th.data(), bin_freqs);
  };

  for (int i = 0; i < dim; i++) {
    const bool clamped =
        (i == eta_idx_ && eta_idx_ >= 0 && theta[i] > 0.25 - eps_);

    std::vector<VECCPL> delta_h(n_modes);
    const int n_edges = static_cast<int>(bin_freqs.size());

    if (order_ == 4) {
      const double vpp = clamped ? theta[i] : theta[i] + 2.0 * eps_;
      const double vp = clamped ? theta[i] : theta[i] + eps_;
      const double vm = theta[i] - eps_;
      const double vmm = theta[i] - 2.0 * eps_;

      auto h_pp = make_modes(i, vpp);
      auto h_p = make_modes(i, vp);
      auto h_m = make_modes(i, vm);
      auto h_mm = make_modes(i, vmm);

      const double denom = clamped ? 6.0 * eps_ : 12.0 * eps_;
      for (int m = 0; m < n_modes; m++) {
        delta_h[m].resize(n_edges);
        for (int b = 0; b < n_edges; b++)
          delta_h[m][b] =
              (-h_pp[m][b] + 8.0 * h_p[m][b] - 8.0 * h_m[m][b] + h_mm[m][b]) /
              denom;
      }
    } else {
      const double vp = clamped ? theta[i] : theta[i] + eps_;
      const double vm = theta[i] - eps_;

      auto h_p = make_modes(i, vp);
      auto h_m = make_modes(i, vm);

      const double denom = clamped ? eps_ : 2.0 * eps_;
      for (int m = 0; m < n_modes; m++) {
        delta_h[m].resize(n_edges);
        for (int b = 0; b < n_edges; b++)
          delta_h[m][b] = (h_p[m][b] - h_m[m][b]) / denom;
      }
    }

    rb_.compute_waveform_ratios(g0[i], g1[i], delta_h);
  }

  // Accumulate Fisher elements using B0/B1 from the RB summary.
  // Gamma_jk = sum_m factor_m * Re[ sum_b conj(g0_j_m_b)*g0_k_m_b * B0_m_b
  //                                      + (conj(g0_j)*g1_k +
  //                                      conj(g1_j)*g0_k)*B1 ]
  for (int j = 0; j < dim; j++) {
    for (int k = 0; k <= j; k++) {
      double gamma_jk = 0.0;
      for (int m = 0; m < n_modes; m++) {
        const VECDBL& B0 = summary[m].B0;
        const VECDBL& B1 = summary[m].B1;
        for (int b = 0; b < n_bins; b++) {
          gamma_jk += factors[m] *
                      std::real(std::conj(g0[j][m][b]) * g0[k][m][b] * B0[b] +
                                (std::conj(g0[j][m][b]) * g1[k][m][b] +
                                 std::conj(g1[j][m][b]) * g0[k][m][b]) *
                                    B1[b]);
        }
      }
      fisherM[j][k] = gamma_jk;
    }
  }

  // Symmetrize.
  for (int j = 0; j < dim; j++)
    for (int k = 0; k < j; k++) fisherM[k][j] = fisherM[j][k];
}

void apply_prior_regularization(double** fisherM, const ParameterMap& pmap) {
  for (const auto& kv : pmap.specs()) {
    if (kv.second.fixed) continue;
    double w = kv.second.hi - kv.second.lo;
    if (w == 0.0) continue;
    int i = pmap.index_of(kv.first);
    fisherM[i][i] += 1.0 / (w * w);
  }
}

}  // namespace gw_fishers
