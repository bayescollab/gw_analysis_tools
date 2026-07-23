#ifndef GWAT_EFPE_H
#define GWAT_EFPE_H

#include "util.h"
#include "waveform_generator.h"
#include <stdexcept>

/*! \brief Frequency-domain waveform for the EFPE eccentric-precessing model.
 *
 * Evaluates the plus and cross polarizations of an eccentric, precessing binary
 * inspiral at the requested frequencies using the pyEFPE SUA. Sets tolerances with
 * pyefpe::ParameterPreset::Production.
 *
 * \param frequencies Array of gravitational-wave frequencies in Hz
 * \param length Number of elements in \p frequencies
 * \param wp Output polarizations; hplus and hcross must be pre-allocated to \p length
 * \param params Source parameters — masses in solar masses, spins dimensionless,
 *               distance in Mpc, \c e_start and \c mean_anomaly_start for eccentricity
 *
 * \return 1 on success
 *
 * \note Only double precision is supported. Instantiation with ADOL-C adouble
 *       throws std::runtime_error at runtime.
 */
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

/*! \brief Uniform-grid frequency-domain waveform for the EFPE eccentric-precessing model.
 *
 * Identical to efpe_fourier_waveform() but uses the pyEFPE uniform-grid generator,
 * which avoids per-frequency allocation and is the preferred path for 
 * where the analysis band is a fixed, evenly-spaced grid. The caller is responsible for
 * ensuring that \p frequencies is uniform; no check is performed.
 *
 * \param frequencies Uniform array of gravitational-wave frequencies in Hz;
 *                    f_min = frequencies[0], delta_f = frequencies[1] - frequencies[0]
 * \param length Number of elements in \p frequencies
 * \param wp Output polarizations; hplus and hcross must be pre-allocated to \p length
 * \param params Source parameters — same convention as efpe_fourier_waveform()
 *
 * \return 1 on success
 *
 * \note Only double precision is supported. Instantiation with ADOL-C adouble
 *       throws std::runtime_error at runtime.
 */
template<class T>
int efpe_fourier_waveform_uniform(T * /*frequencies*/, int /*length*/,
                                  waveform_polarizations<T> * /*wp*/,
                                  gen_params_base<T> * /*params*/)
{
    throw std::runtime_error(
        "EFPE waveform model only supports double precision; "
        "automatic differentiation (ADOL-C) is not supported.");
    return 0;
}

template<>
int efpe_fourier_waveform_uniform<double>(double *frequencies, int length,
                                          waveform_polarizations<double> *wp,
                                          gen_params_base<double> *params);

/*! \brief Frequency-domain waveform for the EFPE model using the FastProduction preset.
 *
 * Identical to efpe_fourier_waveform() but applies
 * pyefpe::ParameterPreset::FastProduction (sua_kmax=2, amplitude_interp_nodes=8),
 * which gives ~1.8x speedup over Production with <2e-4 relative logL error across
 * the LISA circular validation suite.
 */
template<class T>
int efpe_fast_fourier_waveform(T * /*frequencies*/, int /*length*/,
                               waveform_polarizations<T> * /*wp*/,
                               gen_params_base<T> * /*params*/)
{
    throw std::runtime_error(
        "EFPE waveform model only supports double precision; "
        "automatic differentiation (ADOL-C) is not supported.");
    return 0;
}

template<>
int efpe_fast_fourier_waveform<double>(double *frequencies, int length,
                                       waveform_polarizations<double> *wp,
                                       gen_params_base<double> *params);

/*! \brief Uniform-grid frequency-domain waveform for the EFPE model using the FastProduction preset.
 *
 * Identical to efpe_fourier_waveform_uniform() but applies
 * pyefpe::ParameterPreset::FastProduction (sua_kmax=2, amplitude_interp_nodes=8).
 */
template<class T>
int efpe_fast_fourier_waveform_uniform(T * /*frequencies*/, int /*length*/,
                                       waveform_polarizations<T> * /*wp*/,
                                       gen_params_base<T> * /*params*/)
{
    throw std::runtime_error(
        "EFPE waveform model only supports double precision; "
        "automatic differentiation (ADOL-C) is not supported.");
    return 0;
}

template<>
int efpe_fast_fourier_waveform_uniform<double>(double *frequencies, int length,
                                               waveform_polarizations<double> *wp,
                                               gen_params_base<double> *params);

#endif
