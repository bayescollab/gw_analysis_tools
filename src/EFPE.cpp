/*! \file EFPE.cpp
 *  \brief Implementation of the EFPE frequency-domain waveform wrapper.
 *
 *  Maps GWAT source parameters to pyefpe::Parameters, calls pyefpe::Model.
 */
#include "EFPE.h"
#include <pyefpe/pyefpe.hpp>
#include <cmath>
#include <complex>
#include <vector>

namespace {

pyefpe::Parameters build_efpe_params(gen_params_base<double> *p, double f_fallback)
{
    pyefpe::Parameters ep;
    ep.mass1              = p->mass1;
    ep.mass2              = p->mass2;
    ep.spin1x             = p->spin1[0];
    ep.spin1y             = p->spin1[1];
    ep.spin1z             = p->spin1[2];
    ep.spin2x             = p->spin2[0];
    ep.spin2y             = p->spin2[1];
    ep.spin2z             = p->spin2[2];
    ep.distance           = p->Luminosity_Distance;
    ep.inclination        = p->incl_angle;
    ep.phi_start          = p->phiRef;
    ep.e_start            = p->e_start;
    ep.mean_anomaly_start = p->mean_anomaly_start;
    ep.f22_start          = (p->f_ref > 0.0) ? p->f_ref : f_fallback;
    ep.f22_end            = 0.0; // 0 instructs efpe to use ISCO
    return ep;
}

void apply_time_shift(double *frequencies, int length,
                      waveform_polarizations<double> *wp, double tc)
{
    const double two_pi_tc = 2.0 * M_PI * tc;
    for (int i = 0; i < length; i++) {
        double phase = two_pi_tc * frequencies[i];
        std::complex<double> shift(std::cos(phase), std::sin(phase));
        wp->hplus[i]  *= shift;
        wp->hcross[i] *= shift;
    }
}

} // namespace

template<>
int efpe_fourier_waveform<double>(double *frequencies, int length,
                                  waveform_polarizations<double> *wp,
                                  gen_params_base<double> *p)
{
    pyefpe::Parameters ep = build_efpe_params(p, frequencies[0]);
    pyefpe::apply_preset(ep, pyefpe::ParameterPreset::Production);
    pyefpe::Model model(ep);

    std::vector<double> freq_vec(frequencies, frequencies + length);
    pyefpe::FrequencyWaveform fw = model.generate_waveform(freq_vec);

    for (int i = 0; i < length; i++) {
        wp->hplus[i]  = fw.plus[i];
        wp->hcross[i] = fw.cross[i];
    }

    if (p->shift_time && p->tc != 0.0)
        apply_time_shift(frequencies, length, wp, p->tc);

    return 1;
}

template<>
int efpe_fourier_waveform_uniform<double>(double *frequencies, int length,
                                          waveform_polarizations<double> *wp,
                                          gen_params_base<double> *p)
{
    double f_min   = frequencies[0];
    double delta_f = frequencies[1] - frequencies[0];

    pyefpe::Parameters ep = build_efpe_params(p, f_min);
    pyefpe::apply_preset(ep, pyefpe::ParameterPreset::Production);
    pyefpe::Model model(ep);

    model.generate_waveform_uniform(f_min, delta_f, static_cast<std::size_t>(length),
                                    wp->hplus, wp->hcross);

    if (p->shift_time && p->tc != 0.0)
        apply_time_shift(frequencies, length, wp, p->tc);

    return 1;
}
