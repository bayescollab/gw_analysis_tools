#include "gw_prior.h"

#include <cmath>
#include <iostream>

#include "standardPriorLibrary.h"

namespace gw_prior {

BoundsMap lnMc_eta_bounds_from_mass_bounds(double m1lo, double m1hi,
                                           double m2lo, double m2hi) {
  // Mc is monotonically increasing in each mass, so extremes are at corners.
  // Min: equal masses at the shared lower bound.
  double Mc_min = calculate_chirpmass(m2lo, m2lo);
  // Max: both masses at their upper bounds, respecting m2 <= m1.
  double Mc_max = calculate_chirpmass(m1hi, std::min(m1hi, m2hi));

  // eta = m1*m2/(m1+m2)^2; maximised at q=1 (eta=1/4), decreases toward 0.
  // Most extreme mass ratio is at (m1hi, m2lo).
  double eta_min = calculate_eta(m1hi, m2lo);
  double eta_max = 0.25;  // equal masses

  BoundsMap bounds;
  bounds["lnMc"] = {std::log(Mc_min), std::log(Mc_max)};
  bounds["eta"] = {eta_min, eta_max};
  return bounds;
}

BHBPriorFn::BHBPriorFn(const ParamSpecMap& specs, const ParameterMap& pmap,
                       BoundsMap extra_bounds) {
  // Box constraints for every sampled parameter from its ParamSpec bounds.
  for (const auto& kv : specs) {
    if (kv.second.fixed) continue;
    const std::string& name = kv.first;
    int idx = pmap.index_of(name);
    double lo = kv.second.lo, hi = kv.second.hi;
    add([idx, lo, hi](const double* v) {
      return v[idx] >= lo && v[idx] <= hi;
    });
  }

  // lnMc + eta: component-mass bounds + chirpmass-eta Jacobian.
  {
    int lnMc_idx = pmap.index_of("lnMc");
    if (lnMc_idx >= 0) {
      int eta_idx = pmap.index_of("eta");
      add_mass_constraint(extra_bounds, lnMc_idx, eta_idx);
      add_ln_prior_term([lnMc_idx, eta_idx](const double* v) {
        return std::log(chirpmass_eta_jac(std::exp(v[lnMc_idx]), v[eta_idx]));
      });
    }
  }

  // chi1, chi2: aligned-spin Jacobian.
  {
    int chi1_idx = pmap.index_of("chi1");
    if (chi1_idx >= 0) {
      add_ln_prior_term([chi1_idx](const double* v) {
        return std::log(aligned_spin_prior(v[chi1_idx]));
      });
      std::cout << "Aligned chi1 prior added.\n";
    }
    int chi2_idx = pmap.index_of("chi2");
    if (chi2_idx >= 0) {
      add_ln_prior_term([chi2_idx](const double* v) {
        return std::log(aligned_spin_prior(v[chi2_idx]));
      });
      std::cout << "Aligned chi2 prior added.\n";
    }
  }

  // lnDL: uniform-in-volume Jacobian.
  {
    int lnDL_idx = pmap.index_of("lnDL");
    if (lnDL_idx >= 0) {
      add_ln_prior_term(
          [lnDL_idx](const double* v) { return 3.0 * v[lnDL_idx]; });
    }
  }
}

void BHBPriorFn::add_mass_constraint(const BoundsMap& b, int lnMc_idx,
                                     int eta_idx) {
  auto m1_it = b.find("m1"), m2_it = b.find("m2");
  if (m1_it == b.end())
    throw std::invalid_argument("BoundsMap missing key: m1");
  if (m2_it == b.end())
    throw std::invalid_argument("BoundsMap missing key: m2");
  double m1lo = m1_it->second.first, m1hi = m1_it->second.second;
  double m2lo = m2_it->second.first, m2hi = m2_it->second.second;
  add([lnMc_idx, eta_idx, m1lo, m1hi, m2lo, m2hi](const double* v) {
    double Mc = std::exp(v[lnMc_idx]), eta = v[eta_idx];
    auto mpair = component_masses_from_Mc_eta(Mc, eta);
    double m1 = mpair.first, m2 = mpair.second;
    return m1 >= m1lo && m1 <= m1hi && m2 >= m2lo && m2 <= m2hi;
  });
}

}  // namespace gw_prior
