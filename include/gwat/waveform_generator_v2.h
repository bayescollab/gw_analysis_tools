#ifndef WAVEFORM_GENERATOR_V2_H
#define WAVEFORM_GENERATOR_V2_H

/// @file
/// @brief Quasi-modern framework for waveform generation.
///
/// Defines abstract and derived classes for generating GWAT waveforms.
/// To be used with GWLikelihoods and GWFisherDerivatives.
/// Eventually waveform_polarizations<double> would be replaced by a better
/// container. This utilizes old paths for old models (i.e. IMRPhenomD) that
/// should be modernized. The EFPE generator is an example of what a modern
/// definition would be.

#include "util.h"
#include "waveform_generator.h"

namespace waveform_generator {
class WaveformGenerator {
 protected:
  const std::string model_;
  const int n_modes_;

 public:
  ~WaveformGenerator() = default;
  WaveformGenerator(std::string model, int n_modes)
      : model_(std::move(model)), n_modes_(n_modes) {}

  std::string model_name() const { return model_; }
  int n_modes() const { return n_modes_; }

  /// @brief Generate the polarizations.
  ///
  /// Returns the active modes in canonical order
  /// {hplus, hcross, hx, hy, hb, hl}, omitting inactive ones.
  /// For standard GR generators the result always has exactly two entries:
  /// modes[0] = hplus, modes[1] = hcross.
  virtual std::vector<VECCPL> generate_polarizations(
      const gen_params_base<double>* params,
      const VECDBL& frequencies) const = 0;

  /// @brief Fill a C-style waveform_polarizations struct.
  ///
  /// Convenience wrapper for callers that pass polarizations to legacy
  /// functions (e.g. fourier_detector_response_equatorial).
  /// wp must already have memory allocated for frequencies.size() samples.
  /// Modes are written in canonical order: modes[0]→hplus, modes[1]→hcross,
  /// modes[2]→hx, modes[3]→hy, modes[4]→hb, modes[5]→hl.
  void fill_polarizations(waveform_polarizations<double>* wp,
                          const gen_params_base<double>* params,
                          const VECDBL& frequencies) const {
    auto modes = generate_polarizations(params, frequencies);
    std::complex<double>* ptrs[6] = {wp->hplus, wp->hcross, wp->hx,
                                     wp->hy,    wp->hb,     wp->hl};
    for (int m = 0; m < static_cast<int>(modes.size()) && m < 6; ++m)
      if (ptrs[m]) std::copy(modes[m].begin(), modes[m].end(), ptrs[m]);
  }
};
}  // namespace waveform_generator

#endif  // WAVEFORM_GENERATOR_V2_H