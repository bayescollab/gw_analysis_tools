#ifndef PARAMETER_SPACES_H
#define PARAMETER_SPACES_H

/// @file
/// Concrete GWParameterSpace implementations for common waveform families.
///
/// Each class encodes one MCMC parameter layout (vector → gen_params mapping)
/// and the corresponding generation method string.  The quadrature object and
/// fixed injection parameters (f_ref, gmst) are supplied at construction time
/// so the same class can be reused across different grids.
///
/// Constructors accept a BoundsMap keyed by physical parameter names (see each
/// class's static physical_names() for the required keys).  All parameters
/// sampled directly are bounded in sampling coordinates; component masses are
/// checked in physical (solar-mass) units after the (lnMc, eta) → (m1, m2)
/// transform.
///

#include <cmath>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "parameter_space.h"
#include "util.h"

/// @brief Aligned-spin 11-parameter MCMC layout.
///
/// Vector layout:
///   [RA, sin(DEC), psi, cos(iota), phiRef, tc, ln(DL), ln(Mc), eta, χ₁z, χ₂z]
///
/// Extrinsic parameters (indices 0–6): RA, sin(DEC), psi, cos(iota), phiRef,
/// tc, ln(DL).
///
/// BoundsMap keys: see physical_names().
class AlignedSpin_ParamSpace : public GWParameterSpace {
  const Quadrature& quad_;
  const std::string model_;
  const double f_ref_;
  const double gmst_;
  std::unique_ptr<GWPriorFn> prior_fn_obj_;

 public:
  /// Physical-space parameter names that must appear as keys in the BoundsMap
  /// passed to the constructor.  All are in sampling coordinates except "m1"
  /// and "m2" (solar masses).
  static const std::vector<std::string>& physical_names() {
    static const std::vector<std::string> names = {
        "RA",   "sinDEC", "psi",   "cosIota", "phiRef", "tc",
        "lnDL", "eta",    "chi1z", "chi2z",   "m1",     "m2"};
    return names;
  }

  /// @param bounds BoundsMap with keys from physical_names().
  /// @param q      Quadrature method for Fisher inner products.
  /// @param f_ref  Reference frequency passed to the waveform generator.
  /// @param gmst   Greenwich mean sidereal time (radians) at the event epoch.
  /// @param model  Waveform model string (default "IMRPhenomD").
  AlignedSpin_ParamSpace(BoundsMap bounds, const Quadrature& q, double f_ref,
                         double gmst, std::string model = "IMRPhenomD")
      : quad_(q), model_(std::move(model)), f_ref_(f_ref), gmst_(gmst) {
    param_names = {"RA",   "sinDEC", "psi", "cosIota", "phiRef", "tc",
                   "lnDL", "lnMc",   "eta", "chi1z",   "chi2z"};
    extrinsic_params_n = 7;

    prior_fn_obj_ = std::make_unique<GWPriorFn>();
    int lnMc_idx = index_of("lnMc"), eta_idx = index_of("eta");
    for (const auto& name : {"RA", "sinDEC", "psi", "cosIota", "phiRef", "tc",
                             "lnDL", "eta", "chi1z", "chi2z"})
      prior_fn_obj_->add_direct_constraint(bounds, name, index_of(name));
    prior_fn_obj_->add_mass_constraint(bounds, lnMc_idx, eta_idx);

    prior_fn_obj_->add_ln_prior_term(
        [chi1_idx = index_of("chi1z")](const double* v) {
          return std::log(aligned_spin_prior(v[chi1_idx]));
        });
    prior_fn_obj_->add_ln_prior_term(
        [chi2_idx = index_of("chi2z")](const double* v) {
          return std::log(aligned_spin_prior(v[chi2_idx]));
        });
    prior_fn_obj_->add_ln_prior_term([lnMc_idx, eta_idx](const double* v) {
      return std::log(chirpmass_eta_jac(std::exp(v[lnMc_idx]), v[eta_idx]));
    });
    prior_fn_obj_->add_ln_prior_term(
        [lnDL_idx = index_of("lnDL")](const double* v) {
          return 3.0 * v[lnDL_idx];
        });
  }

  bayesship::probabilityFn* prior_fn() const override {
    return prior_fn_obj_.get();
  }

  void to_gen_params(const double* v,
                     gen_params_base<double>& p) const override {
    double Mc = std::exp(v[7]);
    double eta = v[8];
    double M = Mc * std::pow(eta, -0.6);
    double disc = std::sqrt(std::max(0.0, 1.0 - 4.0 * eta));
    p.mass1 = 0.5 * M * (1.0 + disc);
    p.mass2 = 0.5 * M * (1.0 - disc);
    p.spin1[0] = 0.;
    p.spin1[1] = 0.;
    p.spin1[2] = v[9];
    p.spin2[0] = 0.;
    p.spin2[1] = 0.;
    p.spin2[2] = v[10];
    p.Luminosity_Distance = std::exp(v[6]);
    p.incl_angle = std::acos(std::max(-1.0, std::min(1.0, v[3])));
    p.phiRef = v[4];
    p.f_ref = f_ref_;
    p.RA = v[0];
    p.DEC = std::asin(std::max(-1.0, std::min(1.0, v[1])));
    p.psi = v[2];
    p.tc = v[5];
    p.gmst = gmst_;
    p.equatorial_orientation = false;
    p.horizon_coord = false;
    p.shift_time = true;
    p.shift_phase = true;
  }

  void from_gen_params(const gen_params_base<double>& p,
                       double* v) const override {
    v[0] = p.RA;
    v[1] = std::sin(p.DEC);
    v[2] = p.psi;
    v[3] = std::cos(p.incl_angle);
    v[4] = p.phiRef;
    v[5] = p.tc;
    v[6] = std::log(p.Luminosity_Distance);
    v[7] = std::log(calculate_chirpmass(p.mass1, p.mass2));
    v[8] = calculate_eta(p.mass1, p.mass2);
    v[9] = p.spin1[2];
    v[10] = p.spin2[2];
  }

  std::string generation_method() const override { return model_; }
  const Quadrature* quadrature() const override { return &quad_; }
};

/// @brief Spin-precessing 15-parameter MCMC layout.
///
/// Vector layout:
///   [RA, sin(DEC), psi, cos(iota), phiRef, tc, ln(DL), ln(Mc), eta,
///    a₁, a₂, cos_tilt₁, cos_tilt₂, φ₁, φ₂]
///
/// Extrinsic parameters (indices 0–6): RA, sin(DEC), psi, cos(iota), phiRef,
/// tc, ln(DL).  Spins are parameterized by magnitude, cosine of the tilt
/// angle, and azimuthal angle; the conversion uses spin_cartesian_to_spherical.
///
/// BoundsMap keys: see physical_names().
class SpinPrecessing_ParamSpace : public GWParameterSpace {
  const Quadrature& quad_;
  const std::string model_;
  const double f_ref_;
  const double gmst_;
  std::unique_ptr<GWPriorFn> prior_fn_obj_;

 public:
  /// Physical-space parameter names that must appear as keys in the BoundsMap.
  static const std::vector<std::string>& physical_names() {
    static const std::vector<std::string> names = {
        "RA", "sinDEC", "psi",   "cosIota", "phiRef", "tc",   "lnDL", "eta",
        "a1", "a2",     "cosT1", "cosT2",   "phi1",   "phi2", "m1",   "m2"};
    return names;
  }

  /// @param bounds BoundsMap with keys from physical_names().
  /// @param q      Quadrature method for Fisher and likelihood inner products.
  /// @param f_ref  Reference frequency passed to the waveform generator.
  /// @param gmst   Greenwich mean sidereal time (radians) at the event epoch.
  /// @param model  Waveform model string (default "EFPE_circular_fast").
  SpinPrecessing_ParamSpace(BoundsMap bounds, const Quadrature& q, double f_ref,
                            double gmst,
                            std::string model = "EFPE_circular_fast")
      : quad_(q), model_(std::move(model)), f_ref_(f_ref), gmst_(gmst) {
    param_names = {"RA", "sinDEC", "psi",   "cosIota", "phiRef",
                   "tc", "lnDL",   "lnMc",  "eta",     "a1",
                   "a2", "cosT1",  "cosT2", "phi1",    "phi2"};
    extrinsic_params_n = 7;

    prior_fn_obj_ = std::make_unique<GWPriorFn>();
    int lnMc_idx = index_of("lnMc"), eta_idx = index_of("eta");
    for (const auto& name :
         {"RA", "sinDEC", "psi", "cosIota", "phiRef", "tc", "lnDL", "eta", "a1",
          "a2", "cosT1", "cosT2", "phi1", "phi2"})
      prior_fn_obj_->add_direct_constraint(bounds, name, index_of(name));
    prior_fn_obj_->add_mass_constraint(bounds, lnMc_idx, eta_idx);

    prior_fn_obj_->add_ln_prior_term([lnMc_idx, eta_idx](const double* v) {
      return std::log(chirpmass_eta_jac(std::exp(v[lnMc_idx]), v[eta_idx]));
    });
    prior_fn_obj_->add_ln_prior_term(
        [lnDL_idx = index_of("lnDL")](const double* v) {
          return 3.0 * v[lnDL_idx];
        });
  }

  bayesship::probabilityFn* prior_fn() const override {
    return prior_fn_obj_.get();
  }

  void to_gen_params(const double* v,
                     gen_params_base<double>& p) const override {
    double Mc = std::exp(v[7]);
    double eta = v[8];
    double M = Mc * std::pow(eta, -0.6);
    double disc = std::sqrt(std::max(0.0, 1.0 - 4.0 * eta));
    p.mass1 = 0.5 * M * (1.0 + disc);
    p.mass2 = 0.5 * M * (1.0 - disc);

    double a1 = v[9], a2 = v[10];
    double cosT1 = v[11], cosT2 = v[12];
    double sinT1 = std::sqrt(std::max(0.0, 1.0 - cosT1 * cosT1));
    double sinT2 = std::sqrt(std::max(0.0, 1.0 - cosT2 * cosT2));
    p.spin1[0] = a1 * sinT1 * std::cos(v[13]);
    p.spin1[1] = a1 * sinT1 * std::sin(v[13]);
    p.spin1[2] = a1 * cosT1;
    p.spin2[0] = a2 * sinT2 * std::cos(v[14]);
    p.spin2[1] = a2 * sinT2 * std::sin(v[14]);
    p.spin2[2] = a2 * cosT2;

    p.Luminosity_Distance = std::exp(v[6]);
    p.incl_angle = std::acos(std::max(-1.0, std::min(1.0, v[3])));
    p.phiRef = v[4];
    p.f_ref = f_ref_;
    p.RA = v[0];
    p.DEC = std::asin(std::max(-1.0, std::min(1.0, v[1])));
    p.psi = v[2];
    p.tc = v[5];
    p.gmst = gmst_;
    p.equatorial_orientation = false;
    p.horizon_coord = false;
    p.shift_time = true;
    p.shift_phase = true;
  }

  void from_gen_params(const gen_params_base<double>& p,
                       double* v) const override {
    v[0] = p.RA;
    v[1] = std::sin(p.DEC);
    v[2] = p.psi;
    v[3] = std::cos(p.incl_angle);
    v[4] = p.phiRef;
    v[5] = p.tc;
    v[6] = std::log(p.Luminosity_Distance);
    v[7] = std::log(calculate_chirpmass(p.mass1, p.mass2));
    v[8] = calculate_eta(p.mass1, p.mass2);
    spin_cartesian_to_spherical(p.spin1[0], p.spin1[1], p.spin1[2], v[9], v[11],
                                v[13]);
    spin_cartesian_to_spherical(p.spin2[0], p.spin2[1], p.spin2[2], v[10],
                                v[12], v[14]);
  }

  std::string generation_method() const override { return model_; }
  const Quadrature* quadrature() const override { return &quad_; }
};

/// @brief Spin-precessing, no-sky-angles 12-parameter MCMC layout.
///
/// Vector layout:
///   [cos(iota), phiRef, tc, ln(DL), ln(Mc), eta,
///    a₁, a₂, cos_tilt₁, cos_tilt₂, φ₁, φ₂]
///
/// Extrinsic parameters (indices 0–3): cos(iota), phiRef, tc, ln(DL).
///
/// BoundsMap keys: see physical_names().
class SpinPrecessing_NoSky_ParamSpace : public GWParameterSpace {
  const Quadrature& quad_;
  const std::string model_;
  const double f_ref_;
  const double gmst_;
  std::unique_ptr<GWPriorFn> prior_fn_obj_;

 public:
  /// Physical-space parameter names that must appear as keys in the BoundsMap.
  static const std::vector<std::string>& physical_names() {
    static const std::vector<std::string> names = {
        "cosIota", "phiRef", "tc",   "lnDL", "eta", "a1", "a2",
        "cosT1",   "cosT2",  "phi1", "phi2", "m1",  "m2"};
    return names;
  }

  /// @param bounds BoundsMap with keys from physical_names().
  /// @param q      Quadrature method for Fisher and likelihood inner products.
  /// @param f_ref  Reference frequency passed to the waveform generator.
  /// @param gmst   Greenwich mean sidereal time (radians) at the event epoch.
  /// @param model  Waveform model string.
  SpinPrecessing_NoSky_ParamSpace(BoundsMap bounds, const Quadrature& q,
                                  double f_ref, double gmst, std::string model)
      : quad_(q), model_(std::move(model)), f_ref_(f_ref), gmst_(gmst) {
    param_names = {"cosIota", "phiRef", "tc",    "lnDL",  "lnMc", "eta",
                   "a1",      "a2",     "cosT1", "cosT2", "phi1", "phi2"};
    extrinsic_params_n = 3;

    prior_fn_obj_ = std::make_unique<GWPriorFn>();
    int lnMc_idx = index_of("lnMc"), eta_idx = index_of("eta");
    for (const auto& name : {"cosIota", "phiRef", "tc", "lnDL", "eta", "a1",
                             "a2", "cosT1", "cosT2", "phi1", "phi2"})
      prior_fn_obj_->add_direct_constraint(bounds, name, index_of(name));
    prior_fn_obj_->add_mass_constraint(bounds, lnMc_idx, eta_idx);

    prior_fn_obj_->add_ln_prior_term([lnMc_idx, eta_idx](const double* v) {
      return std::log(chirpmass_eta_jac(std::exp(v[lnMc_idx]), v[eta_idx]));
    });
    prior_fn_obj_->add_ln_prior_term(
        [lnDL_idx = index_of("lnDL")](const double* v) {
          return 3.0 * v[lnDL_idx];
        });
  }

  bayesship::probabilityFn* prior_fn() const override {
    return prior_fn_obj_.get();
  }

  void to_gen_params(const double* v,
                     gen_params_base<double>& p) const override {
    double Mc = std::exp(v[4]);
    double eta = v[5];
    double M = Mc * std::pow(eta, -0.6);
    double disc = std::sqrt(std::max(0.0, 1.0 - 4.0 * eta));
    p.mass1 = 0.5 * M * (1.0 + disc);
    p.mass2 = 0.5 * M * (1.0 - disc);

    double a1 = v[6], a2 = v[7];
    double cosT1 = v[8], cosT2 = v[9];
    double sinT1 = std::sqrt(std::max(0.0, 1.0 - cosT1 * cosT1));
    double sinT2 = std::sqrt(std::max(0.0, 1.0 - cosT2 * cosT2));
    p.spin1[0] = a1 * sinT1 * std::cos(v[10]);
    p.spin1[1] = a1 * sinT1 * std::sin(v[10]);
    p.spin1[2] = a1 * cosT1;
    p.spin2[0] = a2 * sinT2 * std::cos(v[11]);
    p.spin2[1] = a2 * sinT2 * std::sin(v[11]);
    p.spin2[2] = a2 * cosT2;

    p.Luminosity_Distance = std::exp(v[3]);
    p.incl_angle = std::acos(std::max(-1.0, std::min(1.0, v[0])));
    p.phiRef = v[1];
    p.f_ref = f_ref_;
    p.tc = v[2];
    p.gmst = gmst_;
    p.equatorial_orientation = false;
    p.horizon_coord = false;
    p.shift_time = true;
    p.shift_phase = true;
  }

  void from_gen_params(const gen_params_base<double>& p,
                       double* v) const override {
    v[0] = std::cos(p.incl_angle);
    v[1] = p.phiRef;
    v[2] = p.tc;
    v[3] = std::log(p.Luminosity_Distance);
    v[4] = std::log(calculate_chirpmass(p.mass1, p.mass2));
    v[5] = calculate_eta(p.mass1, p.mass2);
    spin_cartesian_to_spherical(p.spin1[0], p.spin1[1], p.spin1[2], v[6], v[8],
                                v[10]);
    spin_cartesian_to_spherical(p.spin2[0], p.spin2[1], p.spin2[2], v[7], v[9],
                                v[11]);
  }

  std::string generation_method() const override { return model_; }
  const Quadrature* quadrature() const override { return &quad_; }
};

#endif  // PARAMETER_SPACES_H
