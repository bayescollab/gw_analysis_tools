#include "gw_fisher.h"

#include <memory>

#include "parameter_space.h"
#include "waveform_util.h"

// ── Private concrete derivative class ────────────────────────────────────────

class GWParameterSingleDetectorDerivatives : public GWFisherDerivatives {
 public:
  GWParameterSingleDetectorDerivatives(const GWParameterSpace& ps, int order)
      : ps_(ps), order_(order), eta_idx_(ps.index_of("eta")) {}

  void compute(const double* theta, const IfoData& ifo,
               std::vector<VECCPL>& resp_deriv) const override;

  int dim() const override { return ps_.dim(); }
  const Quadrature& quadrature() const override { return *ps_.quadrature(); }

 private:
  const GWParameterSpace& ps_;
  int order_;
  int eta_idx_;
  static constexpr double eps_ = 1e-8;
};

void GWParameterSingleDetectorDerivatives::compute(
    const double* theta, const IfoData& ifo,
    std::vector<VECCPL>& resp_deriv) const {
  const int d = ps_.dim();
  const std::string gen_method = ps_.generation_method();
  const int len = static_cast<int>(ifo.freqs.size());
  double* freq = const_cast<double*>(ifo.freqs.data());

  resp_deriv.resize(d);

  auto make_wf = [&](int i, double val) -> VECCPL {
    std::vector<double> th(theta, theta + d);
    th[i] = val;
    gen_params_base<double> gp;
    ps_.to_gen_params(th.data(), gp);
    VECCPL h(len);
    create_single_GW_detection(h.data(), ifo.name, freq, len, &gp, gen_method);
    return h;
  };

  for (int i = 0; i < d; i++) {
    resp_deriv[i].resize(len);
    const bool clamped =
        (i == eta_idx_ && eta_idx_ >= 0 && theta[i] > 0.25 - eps_);

    if (order_ == 4) {
      const double vpp = clamped ? theta[i] : theta[i] + 2.0 * eps_;
      const double vp = clamped ? theta[i] : theta[i] + eps_;
      const double vm = theta[i] - eps_;
      const double vmm = theta[i] - 2.0 * eps_;

      VECCPL h_pp = make_wf(i, vpp);
      VECCPL h_p = make_wf(i, vp);
      VECCPL h_m = make_wf(i, vm);
      VECCPL h_mm = make_wf(i, vmm);

      const double denom = clamped ? 6.0 * eps_ : 12.0 * eps_;
      for (int l = 0; l < len; l++)
        resp_deriv[i][l] =
            (-h_pp[l] + 8.0 * h_p[l] - 8.0 * h_m[l] + h_mm[l]) / denom;
    } else {
      const double vp = clamped ? theta[i] : theta[i] + eps_;
      const double vm = theta[i] - eps_;

      VECCPL h_p = make_wf(i, vp);
      VECCPL h_m = make_wf(i, vm);

      const double denom = clamped ? eps_ : 2.0 * eps_;
      for (int l = 0; l < len; l++)
        resp_deriv[i][l] = (h_p[l] - h_m[l]) / denom;
    }
  }
}

// ── GWParameterSpace factory (declared in parameter_space.h) ─────────────────

std::unique_ptr<GWFisherDerivatives> GWParameterSpace::fisher_derivatives(int order) const {
  return std::make_unique<GWParameterSingleDetectorDerivatives>(*this, order);
}

// ── fisher_numerical ─────────────────────────────────────────────────────────

void fisher_numerical(const double* theta, const std::vector<IfoData>& ifos,
                      GWFisherDerivatives* derivatives, double** FisherM) {
  const int dim = derivatives->dim();
  const Quadrature& quad = derivatives->quadrature();

  for (int i = 0; i < dim; i++)
    for (int j = 0; j < dim; j++) FisherM[i][j] = 0.0;

  std::vector<VECCPL> resp_deriv(dim);

  for (const auto& ifo : ifos) {
    derivatives->compute(theta, ifo, resp_deriv);
    for (int j = 0; j < dim; j++)
      for (int k = 0; k <= j; k++)
        FisherM[j][k] +=
            quad.inner_product(resp_deriv[j], resp_deriv[k], ifo.psd);
  }

  for (int j = 0; j < dim; j++)
    for (int k = 0; k < j; k++) FisherM[k][j] = FisherM[j][k];
}
