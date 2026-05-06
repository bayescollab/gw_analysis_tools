#include "EFPE.h"
#include <pyefpe/pyefpe.hpp>
#include <cmath>
#include <complex>
#include <vector>

template<>
int efpe_fourier_waveform<double>(double *frequencies, int length,
                                  waveform_polarizations<double> *wp,
                                  gen_params_base<double> *p)
{
    pyefpe::Parameters ep;
    ep.mass1               = p->mass1;
    ep.mass2               = p->mass2;
    ep.spin1x              = p->spin1[0];
    ep.spin1y              = p->spin1[1];
    ep.spin1z              = p->spin1[2];
    ep.spin2x              = p->spin2[0];
    ep.spin2y              = p->spin2[1];
    ep.spin2z              = p->spin2[2];
    ep.distance            = p->Luminosity_Distance;
    ep.inclination         = p->incl_angle;
    ep.phi_start           = p->phiRef;
    ep.e_start             = p->e_start;
    ep.mean_anomaly_start  = p->mean_anomaly_start;
    ep.f22_start           = (p->f_ref > 0.0) ? p->f_ref : frequencies[0];
    ep.f22_end             = 0.0; // 0 instructs efpe to use ISCO

    pyefpe::Model model(ep);

    std::vector<double> freq_vec(frequencies, frequencies + length);
    pyefpe::FrequencyWaveform fw = model.generate_waveform(freq_vec);

    if (p->shift_time && p->tc != 0.0) {
        const double two_pi_tc = 2.0 * M_PI * p->tc;
        for (int i = 0; i < length; i++) {
            double phase = two_pi_tc * frequencies[i];
            std::complex<double> shift(std::cos(phase), std::sin(phase));
            fw.plus[i]  *= shift;
            fw.cross[i] *= shift;
        }
    }

    for (int i = 0; i < length; i++) {
        wp->hplus[i]  = fw.plus[i];
        wp->hcross[i] = fw.cross[i];
    }

    return 1;
}
