#ifndef WAVEFORM_GENERATORS_H
#define WAVEFORM_GENERATORS_H

/// @file
/// @brief Implementations of WaveformGenerator

#include "EFPE.h"
#include "IMRPhenomD.h"
#include "waveform_generator_v2.h"

namespace waveform_generator {
class IMRPhenomD_generator : public WaveformGenerator {
 public:
  IMRPhenomD_generator() : WaveformGenerator("IMRPhenomD", 2) {}

  std::vector<VECCPL> generate_polarizations(
      const gen_params_base<double>* params,
      const VECDBL& frequencies) const override {
    int length = static_cast<int>(frequencies.size());
    source_parameters<double> converted_params;
    converted_params.populate_source_parameters(
        const_cast<gen_params_base<double>*>(params));
    converted_params.f_ref       = params->f_ref;
    converted_params.shift_time  = params->shift_time;
    converted_params.shift_phase = params->shift_phase;

    VECCPL hplus(length), hcross(length);
    IMRPhenomD<double> phenomd;
    phenomd.construct_waveform(const_cast<double*>(frequencies.data()), length,
                               hplus.data(), &converted_params);

    for (int i = 0; i < length; i++) hcross[i] = CPL(0.0, -1.0) * hplus[i];

    if (params->incl_angle != 0.0) {
      double cos_iota = std::cos(params->incl_angle);
      double hplus_factor = 0.5 * (1.0 + cos_iota * cos_iota);
      for (int i = 0; i < length; i++) {
        hplus[i] *= hplus_factor;
        hcross[i] *= cos_iota;
      }
    }
    return {std::move(hplus), std::move(hcross)};
  }
};

class EFPE_generator : public WaveformGenerator {
 private:
  const bool uniform_;

 public:
  EFPE_generator(bool uniform = false)
      : WaveformGenerator("EFPE", 2), uniform_(uniform) {}

  std::vector<VECCPL> generate_polarizations(
      const gen_params_base<double>* params,
      const VECDBL& frequencies) const override {
    return uniform_ ? efpe_fourier_waveform_uniform(frequencies, params)
                    : efpe_fourier_waveform(frequencies, params);
  }
};

class EFPEFast_generator : public WaveformGenerator {
 private:
  const bool uniform_;

 public:
  EFPEFast_generator(bool uniform = false)
      : WaveformGenerator("EFPE_fast", 2), uniform_(uniform) {}

  std::vector<VECCPL> generate_polarizations(
      const gen_params_base<double>* params,
      const VECDBL& frequencies) const override {
    return uniform_ ? efpe_fast_fourier_waveform_uniform(frequencies, params)
                    : efpe_fast_fourier_waveform(frequencies, params);
  }
};
}  // namespace waveform_generator
#endif  // WAVEFORM_GENERATORS_H
