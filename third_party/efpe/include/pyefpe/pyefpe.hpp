#ifndef PYEFPE_CPP_PYEFPE_HPP
#define PYEFPE_CPP_PYEFPE_HPP

#include <complex>
#include <cstddef>
#include <memory>
#include <vector>

namespace pyefpe {

// Controls how the adaptive RK solution is interpolated between accepted ODE
// steps. FastHermite is cheaper and uses endpoint derivatives. DormandPrince
// uses the RK45 dense-output polynomial and is the closer validation/reference
// choice.
enum class DenseOutput {
  FastHermite,
  DormandPrince
};

// Presets are intentionally coarse-grained. The raw Parameters fields remain
// public so callers can override any individual choice after applying a preset.
enum class ParameterPreset {
  Fast,            // shortest setup/runtime; useful for exploratory runs
  FastProduction,  // Production with sua_kmax=2 and amplitude_interp_nodes=8;
                   // ~1.7x speedup over Production with <2e-4 relative logL error
  Production,      // tuned production waveform settings used by the validation suite
  Validation       // slow reference settings: direct amplitudes and strict tolerances
};

// Physical and numerical parameters for one binary system. Masses are supplied
// in solar masses and are converted internally to geometric seconds. Frequencies
// are f_22 gravitational-wave frequencies in Hz, not orbital frequencies.
struct Parameters {
  double mass1 = 0.0;      // solar masses, required
  double mass2 = 0.0;      // solar masses, required
  double e_start = 0.0;    // eccentricity at f22_start

  // Dimensionless spin vectors in the same convention as the Python pyEFPE
  // code. The constructor swaps labels internally if mass2 > mass1.
  double spin1x = 0.0;
  double spin1y = 0.0;
  double spin1z = 0.0;
  double spin2x = 0.0;
  double spin2y = 0.0;
  double spin2z = 0.0;

  // Spin-induced quadrupole parameters. For black holes these are 1.
  double q1 = 1.0;
  double q2 = 1.0;
  double distance = 100.0; // Mpc
  double inclination = 0.0;
  double f22_start = 20.0;
  double f22_end = 0.0;    // <= 0 means use ISCO
  double phi_start = 0.0;
  double mean_anomaly_start = 0.0;

  // PN truncation controls. The current implementation follows the supported
  // terms from the Python reference and defaults to the highest implemented
  // order.
  int pn_phase_order = 6;
  int pn_spin_order = 6;

  // Mode support controls. amplitude_tol is a cumulative Newtonian-power
  // truncation tolerance: smaller values keep more eccentric harmonics.
  double amplitude_tol = 1.0e-3;
  int amplitude_pmax = 100;

  // Amplitude interpolation caches slowly varying pieces of the mode amplitude
  // on Chebyshev-Lobatto nodes inside each RK segment. This is the default
  // production speedup; disabling it is useful for validation.
  bool interpolate_amplitudes = true;
  int amplitude_interp_nodes = 16;

  // Circular templates are common in recovery studies. When e_start is
  // effectively zero, this option avoids building unused eccentric harmonics and
  // uses a small inverse frequency-to-time cache before Newton refinement.
  bool circular_fast_path = true; // e_start = 0 uses only the circular quadrupole mode
  int circular_spa_iterations = 3; // Newton refinements after the circular inverse cache

  // Mode-parallel uniform-grid waveform generation. Leave this at 1 for small
  // mode counts; high-eccentricity LISA systems benefit most.
  int num_threads = 1;             // uniform-grid frequency generation; 1 is serial

  DenseOutput dense_output = DenseOutput::FastHermite;
  double dj2_tol = 1.0e-10;
  int sua_kmax = 3;

  // Dormand-Prince 5(4) tolerances for
  // [y, e2, lambda, delta_lambda, DeltaJ2, psi_p_bar, phi_z0, zeta0].
  std::vector<double> rr_sol_rtol = {1e-10, 1e-10, 1e-12, 1e-12,
                                     1e-10, 1e-12, 1e-12, 1e-12};
  std::vector<double> rr_sol_atol = {1e-12, 1e-12, 1e-8, 1e-8,
                                     1e-12, 1e-8, 1e-6, 1e-6};

  double spa_frequency_rtol = 1.0e-12;
};

// Apply one of the named presets in-place, or return a preset-modified copy.
// A common pattern is:
//
//   Parameters p = ...;
//   apply_preset(p, ParameterPreset::Production);
//   p.num_threads = 4; // optional caller override
void apply_preset(Parameters& parameters, ParameterPreset preset);
Parameters with_preset(Parameters parameters, ParameterPreset preset);

// Schwarzschild ISCO f_22 frequency for the total mass. This is a termination
// helper; callers often use min(f_ISCO, 0.01 Hz) for one-year LISA tests.
double f22_isco_frequency_hz(double mass1_solar, double mass2_solar);

// Cheap leading-order circular estimate for a start frequency that gives the
// requested duration before f22_end_hz. This is fast but not the full model.
double estimate_f22_start_for_duration_newtonian(double mass1_solar, double mass2_solar,
                                                double f22_end_hz, double duration_seconds,
                                                double min_f22_start_hz = 1.0e-5);

// More accurate start-frequency helper. It repeatedly builds validation models
// while bracketing the requested duration, so use it during setup rather than
// inside a likelihood loop.
double choose_f22_start_for_duration(Parameters parameters, double duration_seconds,
                                     double min_f22_start_hz = 1.0e-5);

// Frequency-domain polarizations on either an arbitrary grid or a uniform grid,
// depending on which generation call produced it.
struct FrequencyWaveform {
  std::vector<std::complex<double>> plus;
  std::vector<std::complex<double>> cross;
};

// Extrinsic-only transformation that can be applied after an intrinsic waveform
// has been generated. This is useful when sampling distance, coalescence time,
// phase, or polarization angle without changing the binary dynamics.
struct ExtrinsicTransform {
  double amplitude_scale = 1.0;
  double time_shift_seconds = 0.0;
  double phase_shift_radians = 0.0;
  double polarization_angle_radians = 0.0;
};

// Apply an extrinsic transformation to a uniform frequency grid. The grid is
// specified by f_i = f_min_hz + i * delta_f_hz.
void apply_extrinsic_transform_uniform(double f_min_hz, double delta_f_hz,
                                       const ExtrinsicTransform& transform,
                                       FrequencyWaveform& waveform);
FrequencyWaveform transformed_uniform(const FrequencyWaveform& waveform, double f_min_hz,
                                      double delta_f_hz,
                                      const ExtrinsicTransform& transform);

struct TimeWaveform {
  std::vector<double> plus;
  std::vector<double> cross;
};

// A Model owns all expensive precomputation for a single intrinsic parameter
// point: initial conditions, radiation-reaction segments, mode lists, and any
// interpolation caches. Build a new Model when intrinsic parameters change.
class Model {
public:
  explicit Model(Parameters parameters);
  ~Model();

  Model(Model&&) noexcept;
  Model& operator=(Model&&) noexcept;

  Model(const Model&) = delete;
  Model& operator=(const Model&) = delete;

  // Arbitrary frequency grid. If frequencies_hz is sorted, the implementation
  // traverses modes by their support intervals; otherwise it falls back to a
  // simpler frequency-by-frequency loop.
  FrequencyWaveform generate_waveform(const std::vector<double>& frequencies_hz) const;

  // Uniform-grid API intended for likelihood code. The pointer overload avoids
  // allocations when the caller owns reusable output buffers.
  FrequencyWaveform generate_waveform_uniform(double f_min_hz, double delta_f_hz,
                                              std::size_t count) const;
  void generate_waveform_uniform(double f_min_hz, double delta_f_hz, std::size_t count,
                                 std::complex<double>* plus,
                                 std::complex<double>* cross) const;

  TimeWaveform generate_time_domain_waveform(const std::vector<double>& times_s) const;

  // Diagnostics and metadata. segment_count counts accepted RK dense-output
  // intervals; mode_count counts segment-local stationary-phase mode intervals.
  double start_time() const;
  double end_time() const;
  std::size_t segment_count() const;
  std::size_t mode_count() const;
  const Parameters& parameters() const;

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

} // namespace pyefpe

#endif // PYEFPE_CPP_PYEFPE_HPP
