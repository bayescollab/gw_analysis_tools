#ifndef IMRPHENOMPV3_H
#define IMRPHENOMPV3_H

#include "IMRPhenomD.h"
#include "IMRPhenomP.h"
#include "util.h"
#include "IMRPhenomPv3utils.h"


/*! \file
 * Header file for IMRPhenomPv3 utilities
 */


template<class T>
class IMRPhenomPv3: public IMRPhenomPv2<T>
{
public:
    virtual int construct_waveform(
                T *frequencies, int length, 
				std::complex<T> *waveform_plus,
				std::complex<T> *waveform_cross,
				source_parameters<T> *params);

    virtual int construct_amplitude(T *frequencies,
				int length,
				T *amplitude,
				source_parameters<T> *params);

private:
    virtual void init_post_merger(source_parameters<T> *params, lambda_parameters<T> *lambda);
};

template <class T> void IMRPhenomPv3_Compute_a_b_e(
    T *alpha, T *beta, T *two_epsilon,
    const T fHz, const T pi_Msec,
    const PhenomPv3Storage<T> *params, const sysprecquant<T> *pAngles);

template <class T> void IMRPhenomPv3twist(
    std::complex<T> *hpTerm, std::complex<T> *hcTerm,
    const T alpha, const sph_harm<T> *Y2m, const T (*wignerD)[2][5]
);

template<class T> void PhenomPComputeWignerD(T d2[5], T dm2[5], T b, source_parameters<T> *params);

#endif /* IMRPHENOMPV3_H */