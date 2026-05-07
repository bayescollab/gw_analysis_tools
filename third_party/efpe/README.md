# pyEFPE C++ Port

This folder contains a standalone C++17 implementation of the core pyEFPE inspiral model from the Python package in `../pyEFPE-main` and the accompanying paper `../pyEFPE.pdf`.

The public API is intentionally small so it can be embedded in a larger C++ codebase:

```cpp
#include "pyefpe/pyefpe.hpp"

#include <algorithm>

pyefpe::Parameters p;
p.mass1 = 2.4;
p.mass2 = 1.2;
p.e_start = 0.7;
p.f22_end = std::min(0.01, pyefpe::f22_isco_frequency_hz(p.mass1, p.mass2));
p.f22_start = pyefpe::choose_f22_start_for_duration(p, 31557600.0);
pyefpe::apply_preset(p, pyefpe::ParameterPreset::Production);

pyefpe::Model model(p);
auto h_fd = model.generate_waveform({20.0, 20.25, 20.5});
auto h_uniform = model.generate_waveform_uniform(0.0, 1.0 / 31557600.0, 1000000);
auto h_td = model.generate_time_domain_waveform({model.start_time(), model.end_time()});

pyefpe::ExtrinsicTransform ext;
ext.amplitude_scale = 0.5;
ext.time_shift_seconds = 10.0;
pyefpe::apply_extrinsic_transform_uniform(0.0, 1.0 / 31557600.0, ext, h_uniform);
```

## Build

This workspace does not have CMake installed, so the port uses a simple Makefile:

```sh
cd cpp
make
make smoke
```

To use the static library from another C++ target:

```sh
c++ -std=c++17 -I/path/to/cpp/include your_code.cpp /path/to/cpp/libpyefpe.a -o your_code
```

No third-party C++ dependencies are required. The implementation uses the platform C math library for integer Bessel functions and includes local Carlson symmetric elliptic-integral routines for the elliptic functions used in the precession dynamics.

## Scope

Implemented:

- parameter handling and the same mass/spin convention as `pyEFPE.pyEFPE`
- initial-condition construction
- 3PN non-spinning and spin terms present in the Python implementation
- adaptive Dormand-Prince 5(4) radiation-reaction integration
- Newtonian eccentric Fourier amplitudes for the `l = 2`, `m = 0, +/-2` modes
- stationary-phase frequency-domain waveform generation with SUA amplitudes
- sorted-grid and uniform-grid frequency-domain generation
- cached amplitude interpolation for fast likelihood-style repeated waveform generation
- direct time-domain waveform generation

## Likelihood Integration Guide

The likelihood code should treat this library as a waveform engine. The library
does not know about detector data, noise PSDs, relative binning, priors, or MCMC
state. Those pieces should stay in the separate likelihood code.

For likelihood work, the usual workflow is:

1. Fill a `pyefpe::Parameters` object with the intrinsic and extrinsic
   parameters for the current sample.
2. Apply `pyefpe::ParameterPreset::Production`.
3. Override any production options needed by your likelihood code, especially
   `num_threads`.
4. Build a `pyefpe::Model`. This performs the expensive intrinsic setup:
   initial conditions, radiation-reaction integration, mode construction, and
   interpolation caches.
5. Call `generate_waveform_uniform(...)` on the same frequency grid as the data.
6. Evaluate the likelihood in your separate likelihood file.

The high-level pattern is:

```cpp
pyefpe::Parameters p;
p.mass1 = m1_solar;
p.mass2 = m2_solar;
p.e_start = e0;
p.spin1x = s1x;
p.spin1y = s1y;
p.spin1z = s1z;
p.spin2x = s2x;
p.spin2y = s2y;
p.spin2z = s2z;
p.inclination = inclination;
p.distance = distance_mpc;
p.f22_end = std::min(0.01, pyefpe::f22_isco_frequency_hz(p.mass1, p.mass2));
p.f22_start = pyefpe::choose_f22_start_for_duration(p, 31557600.0);

pyefpe::apply_preset(p, pyefpe::ParameterPreset::Production);
p.num_threads = 4; // optional; set after apply_preset

pyefpe::Model model(p);

std::vector<std::complex<double>> hp(count);
std::vector<std::complex<double>> hc(count);
model.generate_waveform_uniform(f_min_hz, delta_f_hz, count, hp.data(), hc.data());

// The likelihood loop then compares hp/hc to data using its own PSD, detector
// response, binning, priors, etc.
```

Important likelihood-facing choices:

- Rebuild `pyefpe::Model` when intrinsic parameters change: masses, spins,
  eccentricity, quadrupole parameters, start/end frequencies, PN orders, or
  waveform accuracy settings.
- Do not rebuild `pyefpe::Model` for purely extrinsic changes if you can avoid
  it. Distance/amplitude, coalescence time, constant phase, and polarization
  angle can be handled by `pyefpe::apply_extrinsic_transform_uniform(...)` after
  generating an intrinsic waveform.
- Prefer the pointer-filling overload of `generate_waveform_uniform(...)` inside
  likelihood code. It lets the likelihood own and reuse output buffers.
- Use `pyefpe::ParameterPreset::Production` for sampling. Use
  `pyefpe::ParameterPreset::Validation` only to create reference waveforms or
  debug accuracy; it is intentionally too slow for millions of likelihood calls.
- Choose `f22_start` outside the innermost likelihood loop when possible. The
  model-based `choose_f22_start_for_duration(...)` helper builds several
  validation models internally, so it is a setup tool, not a cheap per-call
  operation.
- If you need a cheaper start-frequency estimate inside exploratory code,
  `estimate_f22_start_for_duration_newtonian(...)` is much faster but only a
  leading-order circular estimate.

## Parameter And Flag Reference

### Presets

- `pyefpe::ParameterPreset::Fast`
  - Purpose: quick exploratory waveform generation.
  - Uses `DenseOutput::FastHermite`, cached amplitudes, the circular fast path,
    and the cheaper radiation-reaction tolerances.
  - Do not use this as the default for final LISA likelihood studies unless you
    have validated the mismatch for your target parameter region.

- `pyefpe::ParameterPreset::Production`
  - Purpose: default for likelihood waveform generation.
  - Uses `DenseOutput::DormandPrince`, cached amplitudes, the circular fast path,
    and adaptive accuracy choices tuned against the 40-system LISA validation
    suite.
  - For circular systems, it keeps only the circular quadrupole harmonic.
  - For high-eccentricity or low-mass moderate-eccentricity systems, it tightens
    the radiation-reaction tolerances and/or mode truncation to keep the
    validation mismatch below `1e-3`.
  - Set `num_threads` after applying this preset if you want threaded waveform
    generation.

- `pyefpe::ParameterPreset::Validation`
  - Purpose: slow reference waveform generation.
  - Uses `DenseOutput::DormandPrince`, direct amplitude evaluation, strict
    radiation-reaction tolerances, `amplitude_tol = 1e-6`, and disables the
    circular fast path.
  - This is the right choice for spot checks, regression tests, and generating
    "exact" reference outputs within the approximate pyEFPE model.
  - This is not meant for production likelihood sampling.

### Physical parameters

- `mass1`, `mass2`
  - Component masses in solar masses.
  - Internally converted to geometric seconds.
  - The constructor swaps labels internally if needed so the heavier object is
    treated as body 1, matching the pyEFPE equations.

- `e_start`
  - Eccentricity at `f22_start`.
  - If `e_start` is effectively zero and `circular_fast_path = true`, the model
    uses the circular fast path.

- `spin1x`, `spin1y`, `spin1z`, `spin2x`, `spin2y`, `spin2z`
  - Dimensionless spin vectors in the same convention as the Python code.
  - These are intrinsic parameters, so changing them requires rebuilding
    `pyefpe::Model`.

- `q1`, `q2`
  - Spin-induced quadrupole parameters.
  - The default value `1` corresponds to black holes.

- `distance`
  - Luminosity distance in Mpc.
  - This only changes the overall amplitude. For likelihood sampling, it can be
    cheaper to generate a reference intrinsic waveform and rescale it externally
    when only distance changes.

- `inclination`
  - Inclination angle in radians.
  - This enters the projection geometry during model construction, so changing
    it currently requires rebuilding `pyefpe::Model`.

- `f22_start`
  - Starting gravitational-wave `m=2` frequency in Hz.
  - This is not orbital frequency; internally the code uses `f_orb = f22 / 2`
    where needed.

- `f22_end`
  - Ending gravitational-wave `m=2` frequency in Hz.
  - If `f22_end <= 0`, the code uses the Schwarzschild ISCO helper.
  - For the LISA validation suite we use `min(f_ISCO, 0.01 Hz)`.

- `phi_start`
  - Initial orbital phase-like angle.
  - A constant phase shift can also be applied afterward with
    `ExtrinsicTransform::phase_shift_radians` if the intrinsic dynamics are
    unchanged.

- `mean_anomaly_start`
  - Initial mean-anomaly-like eccentric phase.
  - This affects the intrinsic eccentric waveform and therefore requires model
    reconstruction when changed.

### Numerical accuracy and speed parameters

- `dense_output`
  - `DenseOutput::FastHermite`: faster interpolation between ODE steps.
  - `DenseOutput::DormandPrince`: RK45 dense-output polynomial; better agreement
    with the Python/SciPy reference.
  - Production and Validation presets use `DormandPrince`.

- `interpolate_amplitudes`
  - If `true`, caches slowly varying amplitude pieces on interpolation nodes
    inside each ODE segment.
  - This should normally remain `true` for likelihood code.
  - Set to `false` for validation to force direct amplitude evaluation.

- `amplitude_interp_nodes`
  - Number of interpolation nodes per ODE segment.
  - Default is `16`.
  - Increasing this can improve amplitude interpolation accuracy but increases
    setup cost and cache size.

- `amplitude_tol`
  - Harmonic truncation tolerance based on cumulative Newtonian eccentric-mode
    power.
  - Smaller values keep more modes and are more accurate but slower.
  - Production changes this automatically depending on eccentricity and mass.

- `amplitude_pmax`
  - Maximum eccentric harmonic search range.
  - Larger values allow more sidebands at high eccentricity but increase setup
    and waveform-generation cost.

- `circular_fast_path`
  - If `true` and `e_start` is effectively zero, the model keeps only the
    circular quadrupole harmonic and uses a cached inverse frequency-to-time map.
  - This is very useful when recovering eccentric injections with circular
    templates, because those recovery templates may be evaluated many times.
  - Disable it only for validation/debugging.

- `circular_spa_iterations`
  - Number of Newton refinements after the circular inverse-cache guess.
  - Default production value is `3`.
  - Increasing this slightly improves stationary-time solve accuracy for
    circular templates but costs time.

- `num_threads`
  - Number of mode-parallel threads for `generate_waveform_uniform(...)`.
  - Default is `1`.
  - Useful for high-eccentricity LISA systems with many modes.
  - Often not useful for circular systems because there are too few modes to
    amortize thread overhead.
  - Set this after applying a preset, because presets clamp or reset threading
    choices.

- `rr_sol_rtol`, `rr_sol_atol`
  - Relative and absolute tolerances for the 8D radiation-reaction ODE:
    `[y, e2, lambda, delta_lambda, DeltaJ2, psi_p_bar, phi_z0, zeta0]`.
  - Production chooses fast or strict tolerance sets based on the parameter
    regime.
  - Usually you should not tune these directly unless validating a new region of
    parameter space.

- `spa_frequency_rtol`
  - Relative tolerance for the stationary-phase frequency equation
    `dPhi/dt = 2*pi*f`.
  - Tighter values improve the stationary-time solve but can cost more Newton
    iterations.

- `sua_kmax`
  - Half-width of the Shifted Uniform Asymptotics amplitude stencil.
  - Default is `3`.
  - Larger values use more shifted amplitude samples around the stationary time.

- `dj2_tol`
  - Tolerance for the initial `DeltaJ2` self-consistency iteration.
  - This is setup-only and normally does not need manual tuning.

### ExtrinsicTransform fields

- `amplitude_scale`
  - Multiplies both polarizations.
  - Useful for distance rescaling or any external amplitude factor.

- `time_shift_seconds`
  - Applies the Fourier-domain factor `exp(-2*pi*i*f*dt)`.
  - Useful for coalescence-time shifts on a fixed uniform grid.

- `phase_shift_radians`
  - Applies a constant complex phase rotation.

- `polarization_angle_radians`
  - Rotates plus/cross by twice the polarization angle.

### Uniform-grid waveform calls

- for likelihood work on a regular Fourier grid, prefer `generate_waveform_uniform(f_min_hz, delta_f_hz, count)` or the pointer-filling overload. This avoids allocating and searching a separate frequency vector.
- `pyefpe::choose_f22_start_for_duration(parameters, 31557600.0)` solves for the start frequency corresponding to a one-sidereal-year observation ending at `parameters.f22_end`, with a default floor of `1e-5 Hz`.
- `pyefpe::estimate_f22_start_for_duration_newtonian(...)` gives a cheap circular leading-order estimate. Use it inside samplers when an approximate start is enough; reserve the model-based solver for setup or validation.
- `pyefpe::apply_extrinsic_transform_uniform(...)` applies distance/amplitude scaling, a frequency-domain time shift, a constant phase shift, and a plus/cross polarization rotation to an already-generated uniform-grid waveform. Use this to avoid regenerating intrinsic waveforms when only extrinsic parameters change.
- `Parameters::interpolate_amplitudes` defaults to `true`, using 16 interpolation nodes per RK segment. Set it to `false` to force direct amplitude evaluation for validation runs.
- `Parameters::num_threads` enables optional mode-parallel uniform-grid generation. Keep it at `1` unless the waveform has many modes; high-eccentricity LISA cases benefit most.
- `Parameters::circular_fast_path` defaults to `true` and restricts `e_start = 0` models to the circular quadrupole harmonic. It also caches a per-segment inverse frequency-to-time map so circular recovery templates avoid most stationary-time Newton work.
- `Parameters::dense_output` defaults to `pyefpe::DenseOutput::FastHermite` for speed. Set it to `pyefpe::DenseOutput::DormandPrince` to use the RK45 quartic dense-output polynomial for closer agreement with the Python/SciPy reference.
- `pyefpe::apply_preset(parameters, pyefpe::ParameterPreset::Fast)` keeps the older fast settings.
- `pyefpe::apply_preset(parameters, pyefpe::ParameterPreset::Production)` uses the RK45 dense output, cached amplitudes, and an adaptive LISA-oriented accuracy setting tuned to stay below a white-noise mismatch of `1e-3` on the 40-system LISA validation suite. It uses a special circular-template branch at `e_start = 0`, fast RR tolerances for most low-eccentricity systems, strict RR tolerances for low-mass moderate-eccentricity and high-eccentricity systems, and tighter mode truncation when high eccentricity needs more harmonic support.
- `pyefpe::apply_preset(parameters, pyefpe::ParameterPreset::Validation)` uses RK45 dense output, direct amplitude evaluation, tighter RR tolerances, and `amplitude_tol = 1e-6`. Use this for reference data rather than likelihood production.

The uniform-grid APIs use

```cpp
f_i = f_min_hz + i * delta_f_hz,  i = 0, ..., count - 1.
```

Use the pointer overload when integrating into a likelihood code:

```cpp
std::vector<std::complex<double>> plus(count);
std::vector<std::complex<double>> cross(count);
model.generate_waveform_uniform(f_min_hz, delta_f_hz, count, plus.data(), cross.data());
```

The returned arrays contain complex Fourier-domain `h_plus(f_i)` and
`h_cross(f_i)`. The detector response, PSD weighting, inner product, and any
relative-binning or reduced-order likelihood machinery should be applied outside
this library.

## Examples And Diagnostic Flags

- `examples/ligo_benchmark.cpp` demonstrates a compact LIGO-like `10 + 5 Msun` eccentric/precessing waveform on an 8 s uniform Fourier grid from 20 Hz to 1024 Hz.
- `examples/lisa_benchmark.cpp` demonstrates a `m1 = 1e3 Msun`, `m2 = 1e4 Msun` eccentric/precessing waveform sampled on a one-year Fourier grid up to a 15 s cadence Nyquist frequency.
- `examples/lisa_validation_suite.cpp` runs 40 LISA-like systems with total masses from `1e3` to `1e5 Msun`, mass ratios `0.1 <= q < 1`, spin directions distributed on the sphere, and eccentricities up to `0.8`. Each system chooses `f22_start` so the model starts one sidereal year before reaching `min(f_ISCO, 0.01 Hz)`, with a floor of `f22_start >= 1e-5 Hz`. It reports the white-noise mismatch between `Validation` and `Production`; add `--write-waveforms` only when you also want validation waveform tables. Use `--force-e=0` for circular recovery-template timing, `--threads=N` for mode-parallel production waveform generation, and `--first=N --last=M` to run a subset:

```sh
cd cpp
make examples/lisa_validation_suite
examples/lisa_validation_suite validation_outputs/lisa40_new_speedups_eccentric --threads=4
examples/lisa_validation_suite validation_outputs/lisa40_new_speedups_circular --force-e=0 --threads=4
```

The validation-suite command-line flags are diagnostic flags for this example
program. They do not change the library globally.

- `OUT_DIR`
  - First non-flag argument.
  - Directory where `summary.json` is written.
  - If omitted, the default is `validation_outputs/lisa_validation_suite`.

- `--write-waveforms`
  - Writes validation waveform tables named `SYSTEM_validation.tsv`.
  - These files can be large.
  - This is for debugging and Python/C++ comparisons, not for normal likelihood
    sampling.

- `--force-e=VALUE`
  - Replaces every generated system's `e_start` by `VALUE`.
  - The suite then recomputes `f22_start` so each system still begins one year
    before `min(f_ISCO, 0.01 Hz)`, respecting the `1e-5 Hz` floor.
  - `--force-e=0` is the circular-template timing/validation mode.

- `--threads=N`
  - Sets `production_p.num_threads = N` in the example.
  - This only affects the Production waveform in the validation comparison.
  - It does not affect the Validation waveform, which remains serial and strict.
  - This is the flag to use when estimating whether threaded uniform-grid
    generation helps your target likelihood systems.

- `--first=N`
  - Runs only systems starting at 1-based index `N`.
  - Useful for rerunning a problematic case without waiting for all 40 systems.

- `--last=M`
  - Runs only through 1-based index `M`.
  - Combine with `--first=N --last=N` to run one system.

For likelihood integration, copy the library usage pattern, not the validation
loop. In particular, the validation mismatch computation is intentionally inside
the example only and should not be placed in the production likelihood loop.

The paper notes the same physical limitations as the Python model: this is an inspiral PN model and does not include merger-ringdown, higher-order modes beyond the Newtonian quadrupole content, tidal interactions, or mode asymmetries.
