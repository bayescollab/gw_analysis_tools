/*! \file test_efpe.cpp
 *  \brief Test that pyefpe waveform matches calling EFPE through GWAT.
 */
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

#include <pyefpe/pyefpe.hpp>

#include "EFPE.h"
#include "util.h"
#include "waveform_generator.h"

int main()
{
    // Parameters from ligo_benchmark.cpp, with f22_end=0 to match the wrapper
    // (GWAT wrapper sets f22_end=0, instructing efpe to use ISCO)
    const double f_min   = 20.0;
    const double delta_f = 1.0 / 8.0;
    const double f_max   = 1024.0;
    const int count =
        static_cast<int>(std::floor((f_max - f_min) / delta_f)) + 1;

    // -----------------------------------------------------------------------
    // Reference: call pyefpe directly, mirroring exactly what the wrapper does
    // -----------------------------------------------------------------------
    pyefpe::Parameters ep;
    ep.mass1       = 10.0;
    ep.mass2       = 5.0;
    ep.e_start     = 0.15;
    ep.spin1x      = 0.25;
    ep.spin1y      = -0.1;
    ep.spin1z      = 0.4;
    ep.spin2x      = -0.05;
    ep.spin2y      = 0.08;
    ep.spin2z      = 0.2;
    ep.inclination = 1.0;
    ep.distance    = 500.0;
    ep.f22_start   = f_min;
    ep.f22_end     = 0.0;
    pyefpe::apply_preset(ep, pyefpe::ParameterPreset::Production);

    pyefpe::Model ref_model(ep);
    pyefpe::FrequencyWaveform ref =
        ref_model.generate_waveform_uniform(f_min, delta_f,
                                            static_cast<std::size_t>(count));

    // -----------------------------------------------------------------------
    // GWAT wrapper: efpe_fourier_waveform_uniform
    // -----------------------------------------------------------------------
    std::vector<double> frequencies(count);
    for (int i = 0; i < count; i++)
        frequencies[i] = f_min + i * delta_f;

    waveform_polarizations<double> wp;
    wp.allocate_memory(count);

    gen_params_base<double> p;
    p.mass1               = 10.0;
    p.mass2               = 5.0;
    p.e_start             = 0.15;
    p.spin1[0]            = 0.25;
    p.spin1[1]            = -0.1;
    p.spin1[2]            = 0.4;
    p.spin2[0]            = -0.05;
    p.spin2[1]            = 0.08;
    p.spin2[2]            = 0.2;
    p.incl_angle          = 1.0;
    p.Luminosity_Distance = 500.0;
    p.f_ref               = f_min;

    efpe_fourier_waveform_uniform(frequencies.data(), count, &wp, &p);

    // -----------------------------------------------------------------------
    // Compare
    // -----------------------------------------------------------------------
    const double tol = 1.0e-14;
    for (int i = 0; i < count; i++) {
        double diff_plus  = std::abs(wp.hplus[i]  - ref.plus[i]);
        double diff_cross = std::abs(wp.hcross[i] - ref.cross[i]);
        
        // Unsuccesful branch
        if (diff_plus > tol || diff_cross > tol) {
            std::cerr << "mismatch at bin " << i
                      << " (f=" << frequencies[i] << " Hz)"
                      << ": |delta hplus|=" << diff_plus
                      << " |delta hcross|=" << diff_cross << "\n";
            wp.deallocate_memory();
            return 1;
        }
    }

    // Success
    wp.deallocate_memory();
    std::cout << "ok\n";
    return 0;
}
