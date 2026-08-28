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

#include "parameter_map.h"
#include "util.h"

namespace gw_prior {
/// Map from physical-space parameter name to [lo, hi] prior bounds.
using BoundsMap = std::map<std::string, PAIRDBL>;

/// @brief Per-parameter [lo, hi] bounds for lnMc and eta derived from
/// component-mass bounds, assuming the convention m1 >= m2 > 0.
///
/// @p m1lo   Lower bound on m1 (also minimum m2 value; must equal m2lo
///           under the shared-lower-bound convention).
/// @p m1hi   Upper bound on m1.
/// @p m2lo   Lower bound on m2 (the lighter mass).
/// @p m2hi   Upper bound on m2.
///
/// Returns a BoundsMap with keys "lnMc" and "eta".
BoundsMap lnMc_eta_bounds_from_mass_bounds(double m1lo, double m1hi,
                                            double m2lo, double m2hi);

/// @brief Formats `x` the same way std::cout would by default. Used to build
/// row cells for print_param_table() before any column widths are measured.
std::string format_value(double x);

/// @brief Prints a name column followed by right-aligned value columns,
/// sized to fit their contents, with a header and separator row. Intended
/// for GWPriorFn classes to report their sampled/fixed parameters
/// consistently during construction. No-op if @p rows is empty.
///
/// @p title         Line printed above the table.
/// @p name_header   Header for the (left-aligned) first column.
/// @p value_headers Headers for each subsequent (right-aligned) column.
/// @p rows          Row cells; row[0] is the name, row[1:] align with
///                   value_headers. Every row must have
///                   value_headers.size() + 1 cells.
void print_param_table(const std::string& title, const std::string& name_header,
                       const std::vector<std::string>& value_headers,
                       const std::vector<std::vector<std::string>>& rows);

/// @brief Uniform box prior built from per-parameter [lo, hi] bounds.
class BoxPrior {
  const std::vector<PAIRDBL> bounds_;

 public:
  explicit BoxPrior(const std::vector<PAIRDBL> b) : bounds_(std::move(b)) {}

  double ln_prior(const double* params) const {
    for (int i = 0; i < static_cast<int>(bounds_.size()); ++i) {
      if (params[i] < bounds_[i].first || params[i] > bounds_[i].second)
        return bayesship::limitInf;
    }
    return 0.0;
  }
};

/// @brief bayesship::probabilityFn built from boolean constraints and additive
/// log-prior terms. A constraint returns True if the parameter position is
/// valid, False otherwise. eval() returns limitInf on the first violated
/// constraint, otherwise the sum of all registered log-prior terms.
class GWPriorFn : public bayesship::probabilityFn {
  std::vector<std::function<bool(const double*)>>   constraints_;
  std::vector<std::function<double(const double*)>> terms_;

 public:
  void add(std::function<bool(const double*)> c) {
    constraints_.push_back(std::move(c));
  }

  void add_ln_prior_term(std::function<double(const double*)> f) {
    terms_.push_back(std::move(f));
  }

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

  double eval(bayesship::positionInfo* pos, int /*chainID*/) override {
    for (const auto& c : constraints_)
      if (!c(pos->parameters)) return bayesship::limitInf;
    double lp = 0.0;
    for (const auto& f : terms_) lp += f(pos->parameters);
    return lp;
  }

  /// @brief Print out all parameters that violate their constraints
  void check_position(bayesship::positionInfo* pos) {
    int i = 0;
    for (const auto& c : constraints_) {
      if (!c(pos->parameters))
        std::cout << "Parameter i=" << i << " failed prior with value "
                  << pos->parameters[i] << "\n";
      i++;
    }
  }
};

/// @brief GWPriorFn for black-hole binary models built from a ParamSpecMap.
///
/// Registers box constraints for all sampled parameters directly from their
/// ParamSpec bounds. Additional physics:
///   - (lnMc, eta): component-mass bounds check + chirpmass-eta Jacobian.
///     Requires "m1" and "m2" entries in @p extra_bounds.
///   - chi1, chi2 (aligned): box constraint + aligned-spin magnitude prior.
///   - a1/cosT1/phi1, a2/cosT2/phi2 (precessing): box constraints only.
///   - lnDL: uniform-volume prior (Jacobian 3*lnDL).
class BHBPriorFn : public GWPriorFn {
 public:
  /// @param specs        Full parameter specification (sampled + fixed).
  /// @param pmap         Compiled ParameterMap; used for index_of() lookups.
  /// @param extra_bounds Additional bounds not in ParamSpec — required to
  ///                     contain "m1" and "m2" whenever "lnMc" is sampled.
  BHBPriorFn(const ParamSpecMap& specs, const ParameterMap& pmap,
              BoundsMap extra_bounds = {});

  /// Bound component masses derived from ln(Mc) and eta. @p b must contain
  /// keys "m1" and "m2".
  void add_mass_constraint(const BoundsMap& b, int lnMc_idx, int eta_idx);
};

}  // namespace gw_prior

#endif  // GW_PRIOR_H
