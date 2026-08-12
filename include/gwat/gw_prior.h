#ifndef GW_PRIOR_H
#define GW_PRIOR_H

#include <bayesship/bayesshipSampler.h>
#include <bayesship/dataUtilities.h>

#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "util.h"

/// Map from physical-space parameter name to [lo, hi] prior bounds.
using BoundsMap = std::map<std::string, DBLPAIR>;

/// @brief bayesship::probabilityFn built from a list of boolean constraints
/// and additive log-prior terms.
///
/// eval() returns bayesship::limitInf on the first violated constraint,
/// otherwise the sum of all registered log-prior terms (0.0 if none).
/// Constraints are added via add(); log-prior contributions (e.g. Jacobians)
/// are added via add_ln_prior_term().
class GWPriorFn : public bayesship::probabilityFn {
  std::vector<std::function<bool(const double*)>> constraints_;
  std::vector<std::function<double(const double*)>> terms_;

 public:
  /// Append a boolean constraint.
  void add(std::function<bool(const double*)> c) {
    constraints_.push_back(std::move(c));
  }

  /// Append an additive log-prior term (e.g. a Jacobian contribution).
  void add_ln_prior_term(std::function<double(const double*)> f) {
    terms_.push_back(std::move(f));
  }

  double eval(bayesship::positionInfo* pos, int /*chainID*/) override {
    for (const auto& c : constraints_)
      if (!c(pos->parameters)) return bayesship::limitInf;
    double lp = 0.0;
    for (const auto& f : terms_) lp += f(pos->parameters);
    return lp;
  }

  // -----------------------------------------------------------------------
  // Helpers — called from concrete subclass constructors to populate
  // a GWPriorFn with the transforms appropriate to each parameter layout.
  // -----------------------------------------------------------------------

  /// Bound the sampling-coordinate parameter at @p idx directly:
  /// lo ≤ v[idx] ≤ hi.  The key in @p b must match the sampling name.
  void add_direct_constraint(const BoundsMap& b, const std::string& name,
                             int idx) {
    auto it = b.find(name);
    if (it == b.end())
      throw std::invalid_argument("BoundsMap missing key: " + name);
    double lo = it->second.first, hi = it->second.second;
    add([idx, lo, hi](const double* v) {
      return v[idx] >= lo && v[idx] <= hi;
    });
  }

  /// Bound component masses derived from ln(Mc) and eta.  @p b must contain
  /// keys "m1" and "m2" with bounds in physical (solar-mass) units.
  void add_mass_constraint(const BoundsMap& b, int lnMc_idx, int eta_idx) {
    auto m1_it = b.find("m1"), m2_it = b.find("m2");
    if (m1_it == b.end())
      throw std::invalid_argument("BoundsMap missing key: m1");
    if (m2_it == b.end())
      throw std::invalid_argument("BoundsMap missing key: m2");
    double m1lo = m1_it->second.first, m1hi = m1_it->second.second;
    double m2lo = m2_it->second.first, m2hi = m2_it->second.second;
    add([lnMc_idx, eta_idx, m1lo, m1hi, m2lo, m2hi](const double* v) {
      double Mc = std::exp(v[lnMc_idx]), eta = v[eta_idx];
      double M = Mc * std::pow(eta, -0.6);
      double disc = std::sqrt(std::max(0.0, 1.0 - 4.0 * eta));
      double m1 = 0.5 * M * (1.0 + disc), m2 = 0.5 * M * (1.0 - disc);
      return m1 >= m1lo && m1 <= m1hi && m2 >= m2lo && m2 <= m2hi;
    });
  }
};

#endif  // GW_PRIOR_H
