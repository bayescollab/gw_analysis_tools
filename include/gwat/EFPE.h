#ifndef GWAT_EFPE_H
#define GWAT_EFPE_H

#include "util.h"
#include "waveform_generator.h"
#include <stdexcept>

// Default template: throws at runtime for any non-double type (e.g. ADOL-C adouble).
// EFPE does not support automatic differentiation.
template<class T>
int efpe_fourier_waveform(T * /*frequencies*/, int /*length*/,
                          waveform_polarizations<T> * /*wp*/,
                          gen_params_base<T> * /*params*/)
{
    throw std::runtime_error(
        "EFPE waveform model only supports double precision; "
        "automatic differentiation (ADOL-C) is not supported.");
    return 0;
}

template<>
int efpe_fourier_waveform<double>(double *frequencies, int length,
                                  waveform_polarizations<double> *wp,
                                  gen_params_base<double> *params);

#endif
