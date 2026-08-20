#ifndef PARAMETER_MAP_H
#define PARAMETER_MAP_H

#include <functional>
#include <string>
#include <vector>

#include "util.h"

enum class ParamCategory { Intrinsic, Extrinsic };

struct ParamSpec {
  bool fixed;
  double value;            // fixed: the pinned value; sampled: unused
  double lo, hi;           // sampled: prior bounds; fixed: unused
  ParamCategory category;  // only used when !fixed

  static ParamSpec Fixed(double v) {
    return {true, v, 0., 0., ParamCategory::Intrinsic};
  }
  static ParamSpec Sampled(double lo, double hi,
                           ParamCategory cat = ParamCategory::Intrinsic) {
    return {false, 0., lo, hi, cat};
  }
};

// Ordered vector of (name, spec) pairs. Insertion order defines the
// sampling-vector layout: sampled params appear in the order they are listed,
// fixed params are omitted from the vector.
using ParamSpecMap = std::vector<std::pair<std::string, ParamSpec>>;

/// @brief Compiled parameter-space mapping built from a ParamSpecMap.
///
/// Translates between the flat MCMC sampling vector and
/// gen_params_base<double>.
///
/// Construction-time arguments:
///   specs             — ordered (name, ParamSpec) pairs defining the model
///   generation_method — waveform generation method string
///
/// Known parameter names handled by the registry and group transforms:
///   Direct 1-to-1: mass1, mass2, phiRef, RA, DEC, psi, tc, incl_angle,
///                  spin1x/y/z, spin2x/y/z, e_start, mean_anomaly_start, f_ref
///   Single reparameterized: lnDL→Luminosity_Distance, sinDEC→DEC,
///                           cosIota→incl_angle
///   Groups (coupled): lnMc+eta→mass1/mass2,
///                     chi1+chi2→aligned spins (spin[0]=spin[1]=0),
///                     a1+cosT1+phi1→spin1 cartesian,
///                     a2+cosT2+phi2→spin2 cartesian
class ParameterMap {
 public:
  explicit ParameterMap(const ParamSpecMap& specs);

  /// Map sampler vector to gen_params_base
  void to_gen_params(const double* vec, gen_params_base<double>& p) const;
  /// Inverse of to_gen_params()
  void from_gen_params(const gen_params_base<double>& p, double* vec) const;

  int dim() const { return dim_; }

  /// Index of @p name in the sampling vector, or -1 if fixed/absent.
  int index_of(const std::string& name) const;

  /// Sampled parameter names in sampling-vector order (fixed params excluded).
  const std::vector<std::string>& names() const { return sampled_names_; }

  /// Indices of extrinsic sampled parameters in the sampling vector.
  const std::vector<int>& extrinsic_indices() const { return extrinsic_idx_; }

  /// Indices of intrinsic sampled parameters in the sampling vector.
  const std::vector<int>& intrinsic_indices() const { return intrinsic_idx_; }

  /// The original ParamSpecMap passed at construction (used for prior widths).
  const ParamSpecMap& specs() const { return specs_; }

 private:
  using Setter = std::function<void(gen_params_base<double>&, double)>;
  using Getter = std::function<double(const gen_params_base<double>&)>;

  struct CompiledEntry {
    int sampling_idx;  // >=0: position in vec; -1: fixed
    double fixed_val;
    Setter set;  // null → handled by a group transform
    Getter get;  // null → handled by a group transform
  };

  struct GroupTransform {
    std::function<void(const double*, gen_params_base<double>&)> to_gen;
    std::function<void(const gen_params_base<double>&, double*)> from_gen;
  };

  ParamSpecMap specs_;
  std::vector<CompiledEntry> entries_;
  std::vector<GroupTransform> groups_;
  std::vector<std::string> sampled_names_;
  std::vector<int> extrinsic_idx_;
  std::vector<int> intrinsic_idx_;
  int dim_ = 0;
};

#endif  // PARAMETER_MAP_H
