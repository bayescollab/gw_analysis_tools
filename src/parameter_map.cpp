#include "parameter_map.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <unordered_map>

#include "util.h"

// ── Static registry ─────────────────────────────────────────────────────────

namespace {

using Setter = std::function<void(gen_params_base<double>&, double)>;
using Getter = std::function<double(const gen_params_base<double>&)>;

struct DirectAccessor {
  Setter set;  // null → this name is handled entirely by a group transform
  Getter get;  // null → this name is handled entirely by a group transform
};

// Meyer's singleton: built once, thread-safe in C++11 and later.
const std::unordered_map<std::string, DirectAccessor>& registry() {
  static const std::unordered_map<std::string, DirectAccessor> reg = {
      // ── Direct 1-to-1 mappings ────────────────────────────────────────────
      {"mass1",
       {[](auto& p, double v) { p.mass1 = v; },
        [](const auto& p) { return p.mass1; }}},
      {"mass2",
       {[](auto& p, double v) { p.mass2 = v; },
        [](const auto& p) { return p.mass2; }}},
      {"phiRef",
       {[](auto& p, double v) { p.phiRef = v; },
        [](const auto& p) { return p.phiRef; }}},
      {"RA",
       {[](auto& p, double v) { p.RA = v; },
        [](const auto& p) { return p.RA; }}},
      {"DEC",
       {[](auto& p, double v) { p.DEC = v; },
        [](const auto& p) { return p.DEC; }}},
      {"psi",
       {[](auto& p, double v) { p.psi = v; },
        [](const auto& p) { return p.psi; }}},
      {"tc",
       {[](auto& p, double v) { p.tc = v; },
        [](const auto& p) { return p.tc; }}},
      {"incl_angle",
       {[](auto& p, double v) { p.incl_angle = v; },
        [](const auto& p) { return p.incl_angle; }}},
      {"spin1x",
       {[](auto& p, double v) { p.spin1[0] = v; },
        [](const auto& p) { return p.spin1[0]; }}},
      {"spin1y",
       {[](auto& p, double v) { p.spin1[1] = v; },
        [](const auto& p) { return p.spin1[1]; }}},
      {"spin1z",
       {[](auto& p, double v) { p.spin1[2] = v; },
        [](const auto& p) { return p.spin1[2]; }}},
      {"spin2x",
       {[](auto& p, double v) { p.spin2[0] = v; },
        [](const auto& p) { return p.spin2[0]; }}},
      {"spin2y",
       {[](auto& p, double v) { p.spin2[1] = v; },
        [](const auto& p) { return p.spin2[1]; }}},
      {"spin2z",
       {[](auto& p, double v) { p.spin2[2] = v; },
        [](const auto& p) { return p.spin2[2]; }}},
      {"e_start",
       {[](auto& p, double v) { p.e_start = v; },
        [](const auto& p) { return p.e_start; }}},
      {"mean_anomaly_start",
       {[](auto& p, double v) { p.mean_anomaly_start = v; },
        [](const auto& p) { return p.mean_anomaly_start; }}},
      {"f_ref",
       {[](auto& p, double v) { p.f_ref = v; },
        [](const auto& p) { return p.f_ref; }}},

      // ── Single reparameterized mappings ───────────────────────────────────
      // lnDL  →  Luminosity_Distance = exp(lnDL)
      {"lnDL",
       {[](auto& p, double v) { p.Luminosity_Distance = std::exp(v); },
        [](const auto& p) { return std::log(p.Luminosity_Distance); }}},
      // sinDEC  →  DEC = asin(clamp(sinDEC, -1, 1))
      {"sinDEC",
       {[](auto& p, double v) {
          p.DEC = std::asin(std::max(-1.0, std::min(1.0, v)));
        },
        [](const auto& p) { return std::sin(p.DEC); }}},
      // cosIota  →  incl_angle = acos(clamp(cosIota, -1, 1))
      {"cosIota",
       {[](auto& p, double v) {
          p.incl_angle = std::acos(std::max(-1.0, std::min(1.0, v)));
        },
        [](const auto& p) { return std::cos(p.incl_angle); }}},

      // ── Group-handled names: null set/get ─────────────────────────────────
      // lnMc and eta are decoded together in the lnMc+eta group transform.
      {"lnMc", {nullptr, nullptr}},
      {"eta", {nullptr, nullptr}},
      // chi1 and chi2 are decoded together in the chi1+chi2 group transform.
      {"chi1", {nullptr, nullptr}},
      {"chi2", {nullptr, nullptr}},
      // Precessing spin 1: group transform handles all three.
      {"a1", {nullptr, nullptr}},
      {"cosT1", {nullptr, nullptr}},
      {"phi1", {nullptr, nullptr}},
      // Precessing spin 2: group transform handles all three.
      {"a2", {nullptr, nullptr}},
      {"cosT2", {nullptr, nullptr}},
      {"phi2", {nullptr, nullptr}},
  };
  return reg;
}

}  // namespace

// ── Constructor───────────────────────────────────────────────────────

ParameterMap::ParameterMap(const ParamSpecMap& specs) : specs_(specs) {
  const auto& reg = registry();

  // Helper: scan specs for a named entry and return its resolved state.
  struct Resolved {
    bool in_spec;
    int idx;  // sampling index (>=0) if sampled; -1 if fixed/absent
    double fixed_val;
  };

  auto resolve = [&](const std::string& name) -> Resolved {
    for (const auto& kv : specs) {
      if (kv.first == name) {
        const ParamSpec& sp = kv.second;
        if (sp.fixed) {
          return {true, -1, sp.value};
        } else {
          // Find its position among sampled names built so far.
          for (int i = 0; i < (int)sampled_names_.size(); ++i)
            if (sampled_names_[i] == name) return {true, i, 0.0};
          // Not yet assigned — this shouldn't happen if we process in order.
          return {true, -1, 0.0};
        }
      }
    }
    return {false, -1, 0.0};
  };

  // ── First pass: assign sampling indices in insertion order ────────────────
  // We need sampled_names_ populated before resolve() can look up indices, so
  // walk specs once to number the sampled parameters.
  for (const auto& kv : specs) {
    if (!kv.second.fixed) sampled_names_.push_back(kv.first);
  }
  dim_ = static_cast<int>(sampled_names_.size());

  // Build extrinsic/intrinsic index lists.
  {
    int sidx = 0;
    for (const auto& kv : specs) {
      if (kv.second.fixed) continue;
      if (kv.second.category == ParamCategory::Extrinsic)
        extrinsic_idx_.push_back(sidx);
      else
        intrinsic_idx_.push_back(sidx);
      ++sidx;
    }
  }

  // ── Detect which groups are present ───────────────────────────────────────
  bool has_lnMc_eta = false;
  bool has_chi1_chi2 = false;
  bool has_spin1_sph = false;
  bool has_spin2_sph = false;

  for (const auto& kv : specs) {
    const std::string& n = kv.first;
    if (n == "lnMc" || n == "eta") has_lnMc_eta = true;
    if (n == "chi1" || n == "chi2") has_chi1_chi2 = true;
    if (n == "a1" || n == "cosT1" || n == "phi1") has_spin1_sph = true;
    if (n == "a2" || n == "cosT2" || n == "phi2") has_spin2_sph = true;
  }

  // ── Second pass: compile entries ───────────────────────────
  for (const auto& kv : specs) {
    const std::string& name = kv.first;
    const ParamSpec& sp = kv.second;

    // Find sampling index from sampled_names_.
    int sidx = -1;
    if (!sp.fixed) {
      for (int i = 0; i < (int)sampled_names_.size(); ++i) {
        if (sampled_names_[i] == name) {
          sidx = i;
          break;
        }
      }
    }

    // Look up the registry.
    auto it = reg.find(name);
    if (it == reg.end())
      throw std::invalid_argument("ParameterMap: unknown parameter name '" +
                                  name + "'");

    CompiledEntry ce;
    ce.sampling_idx = sidx;
    ce.fixed_val = sp.fixed ? sp.value : 0.0;
    ce.set = it->second.set;
    ce.get = it->second.get;
    entries_.push_back(ce);
  }

  // ── Register group transforms ────────────────────────────

  // Group: lnMc + eta → mass1, mass2
  if (has_lnMc_eta) {
    Resolved lnMc_r = resolve("lnMc");
    Resolved eta_r = resolve("eta");

    int lnMc_idx = lnMc_r.in_spec ? lnMc_r.idx : -1;
    double lnMc_v = lnMc_r.in_spec ? lnMc_r.fixed_val : 0.0;
    int eta_idx = eta_r.in_spec ? eta_r.idx : -1;
    double eta_v = eta_r.in_spec ? eta_r.fixed_val : 0.0;

    GroupTransform gt;
    gt.to_gen = [lnMc_idx, lnMc_v, eta_idx, eta_v](const double* vec,
                                                   gen_params_base<double>& p) {
      double lnMc = (lnMc_idx >= 0) ? vec[lnMc_idx] : lnMc_v;
      double eta = (eta_idx >= 0) ? vec[eta_idx] : eta_v;
      auto mpair = component_masses_from_Mc_eta(std::exp(lnMc), eta);
      p.mass1 = mpair.first;
      p.mass2 = mpair.second;
    };
    gt.from_gen = [lnMc_idx, eta_idx](const gen_params_base<double>& p,
                                      double* vec) {
      if (lnMc_idx >= 0)
        vec[lnMc_idx] = std::log(calculate_chirpmass(p.mass1, p.mass2));
      if (eta_idx >= 0) vec[eta_idx] = calculate_eta(p.mass1, p.mass2);
    };
    groups_.push_back(gt);
  }

  // Group: chi1 + chi2 → aligned spins (transverse components zeroed)
  if (has_chi1_chi2) {
    Resolved chi1_r = resolve("chi1");
    Resolved chi2_r = resolve("chi2");

    int chi1_idx = chi1_r.in_spec ? chi1_r.idx : -1;
    double chi1_v = chi1_r.in_spec ? chi1_r.fixed_val : 0.0;
    int chi2_idx = chi2_r.in_spec ? chi2_r.idx : -1;
    double chi2_v = chi2_r.in_spec ? chi2_r.fixed_val : 0.0;

    GroupTransform gt;
    gt.to_gen = [chi1_idx, chi1_v, chi2_idx, chi2_v](
                    const double* vec, gen_params_base<double>& p) {
      double chi1 = (chi1_idx >= 0) ? vec[chi1_idx] : chi1_v;
      double chi2 = (chi2_idx >= 0) ? vec[chi2_idx] : chi2_v;
      p.spin1[0] = 0.0;
      p.spin1[1] = 0.0;
      p.spin1[2] = chi1;
      p.spin2[0] = 0.0;
      p.spin2[1] = 0.0;
      p.spin2[2] = chi2;
    };
    gt.from_gen = [chi1_idx, chi2_idx](const gen_params_base<double>& p,
                                       double* vec) {
      if (chi1_idx >= 0) vec[chi1_idx] = p.spin1[2];
      if (chi2_idx >= 0) vec[chi2_idx] = p.spin2[2];
    };
    groups_.push_back(gt);
  }

  // Group: a1 + cosT1 + phi1 → spin1 cartesian
  if (has_spin1_sph) {
    Resolved a1_r = resolve("a1");
    Resolved cosT1_r = resolve("cosT1");
    Resolved phi1_r = resolve("phi1");

    int a1_idx = a1_r.in_spec ? a1_r.idx : -1;
    double a1_v = a1_r.in_spec ? a1_r.fixed_val : 0.0;
    int cosT1_idx = cosT1_r.in_spec ? cosT1_r.idx : -1;
    double cosT1_v = cosT1_r.in_spec ? cosT1_r.fixed_val : 0.0;
    int phi1_idx = phi1_r.in_spec ? phi1_r.idx : -1;
    double phi1_v = phi1_r.in_spec ? phi1_r.fixed_val : 0.0;

    GroupTransform gt;
    gt.to_gen = [a1_idx, a1_v, cosT1_idx, cosT1_v, phi1_idx, phi1_v](
                    const double* vec, gen_params_base<double>& p) {
      double a = (a1_idx >= 0) ? vec[a1_idx] : a1_v;
      double cosT = (cosT1_idx >= 0) ? vec[cosT1_idx] : cosT1_v;
      double phi = (phi1_idx >= 0) ? vec[phi1_idx] : phi1_v;
      spin_spherical_to_cartestion(a, cosT, phi, p.spin1[0], p.spin1[1],
                                   p.spin1[2]);
    };
    gt.from_gen = [a1_idx, cosT1_idx, phi1_idx](
                      const gen_params_base<double>& p, double* vec) {
      double a, cosT, phi;
      spin_cartesian_to_spherical(p.spin1[0], p.spin1[1], p.spin1[2], a, cosT,
                                  phi);
      if (a1_idx >= 0) vec[a1_idx] = a;
      if (cosT1_idx >= 0) vec[cosT1_idx] = cosT;
      if (phi1_idx >= 0) vec[phi1_idx] = phi;
    };
    groups_.push_back(gt);
  }

  // Group: a2 + cosT2 + phi2 → spin2 cartesian
  if (has_spin2_sph) {
    Resolved a2_r = resolve("a2");
    Resolved cosT2_r = resolve("cosT2");
    Resolved phi2_r = resolve("phi2");

    int a2_idx = a2_r.in_spec ? a2_r.idx : -1;
    double a2_v = a2_r.in_spec ? a2_r.fixed_val : 0.0;
    int cosT2_idx = cosT2_r.in_spec ? cosT2_r.idx : -1;
    double cosT2_v = cosT2_r.in_spec ? cosT2_r.fixed_val : 0.0;
    int phi2_idx = phi2_r.in_spec ? phi2_r.idx : -1;
    double phi2_v = phi2_r.in_spec ? phi2_r.fixed_val : 0.0;

    GroupTransform gt;
    gt.to_gen = [a2_idx, a2_v, cosT2_idx, cosT2_v, phi2_idx, phi2_v](
                    const double* vec, gen_params_base<double>& p) {
      double a = (a2_idx >= 0) ? vec[a2_idx] : a2_v;
      double cosT = (cosT2_idx >= 0) ? vec[cosT2_idx] : cosT2_v;
      double phi = (phi2_idx >= 0) ? vec[phi2_idx] : phi2_v;
      spin_spherical_to_cartestion(a, cosT, phi, p.spin2[0], p.spin2[1],
                                   p.spin2[2]);
    };
    gt.from_gen = [a2_idx, cosT2_idx, phi2_idx](
                      const gen_params_base<double>& p, double* vec) {
      double a, cosT, phi;
      spin_cartesian_to_spherical(p.spin2[0], p.spin2[1], p.spin2[2], a, cosT,
                                  phi);
      if (a2_idx >= 0) vec[a2_idx] = a;
      if (cosT2_idx >= 0) vec[cosT2_idx] = cosT;
      if (phi2_idx >= 0) vec[phi2_idx] = phi;
    };
    groups_.push_back(gt);
  }
}

// ── to_gen_params ───────────────────────────────

void ParameterMap::to_gen_params(const double* vec,
                                 gen_params_base<double>& p) const {
  // Phase 1: direct entries (set != null means not group-handled).
  for (const auto& e : entries_) {
    if (!e.set) continue;  // group-handled
    double val = (e.sampling_idx >= 0) ? vec[e.sampling_idx] : e.fixed_val;
    e.set(p, val);
  }

  // Phase 2: group transforms.
  for (const auto& g : groups_) g.to_gen(vec, p);
}

// ── from_gen_params ───────────────────────────

void ParameterMap::from_gen_params(const gen_params_base<double>& p,
                                   double* vec) const {
  // Direct entries that are sampled (sampling_idx >= 0) and have a getter.
  for (const auto& e : entries_) {
    if (!e.get || e.sampling_idx < 0) continue;
    vec[e.sampling_idx] = e.get(p);
  }

  // Group transforms.
  for (const auto& g : groups_) g.from_gen(p, vec);
}

// ── index_of ──────────────────────────────

int ParameterMap::index_of(const std::string& name) const {
  for (int i = 0; i < (int)sampled_names_.size(); ++i)
    if (sampled_names_[i] == name) return i;
  return -1;
}
