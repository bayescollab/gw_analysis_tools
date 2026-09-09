#ifndef PPE_IMRPHENOMD_H
#define PPE_IMRPHENOMD_H
#include "IMRPhenomD.h"
#include "util.h"

/*! \file 
 */

/*! Class that extends the IMRPhenomD waveform to include non-GR terms in the inspiral portion of
 * the phase. This is an appropriate waveform choice for generation effects, but not necessarily for
 * propagation effects
 */
 
template<class T> 
class ppE_IMRPhenomD_Inspiral: public IMRPhenomD<T>
{
public:
virtual T phase_ins(T f, source_parameters<T> *param, T *pn_coeff, 
		lambda_parameters<T> *lambda, useful_powers<T> *pow);

virtual T Dphase_ins(T f, source_parameters<T> *param, T *pn_coeff, lambda_parameters<T> *lambda);

};
/*! Class that extends the IMRPhenomD waveform to include non-GR terms in the full
 * phase. This is an appropriate waveform choice for propagation effects
 */
 
template<class T> 
class ppE_IMRPhenomD_IMR: public ppE_IMRPhenomD_Inspiral<T>
{
public:
virtual T Dphase_mr(T f, source_parameters<T> *param, lambda_parameters<T> *lambda);
virtual T phase_mr(T f, source_parameters<T> *param, lambda_parameters<T> *lambda);

virtual T phase_int(T f, source_parameters<T> *param, lambda_parameters<T> *lambda );
virtual T Dphase_int(T f, source_parameters<T> *param, lambda_parameters<T> *lambda);

};

//template<class T> 
//class dCS_IMRPhenomD_log: public ppE_IMRPhenomD_Inspiral<T>
//{
//public:
//virtual int construct_waveform(T *frequencies, int length, std::complex<T> *waveform, source_parameters<T> *params);
//virtual T dCS_phase_mod( source_parameters<T> *param);
//virtual T dCS_phase_factor(source_parameters<T> *param);
//
//virtual int construct_amplitude(T *frequencies, int length, T *amplitude, source_parameters<T> *params);
//
////virtual T construct_amplitude(T frequency,  source_parameters<T> *params);
//
//virtual int construct_phase(T *frequencies, int length, T *phase, source_parameters<T> *params);
//};
//template<class T> 
//class dCS_IMRPhenomD: public ppE_IMRPhenomD_Inspiral<T>
//{
//public:
//virtual int construct_waveform(T *frequencies, int length, std::complex<T> *waveform, source_parameters<T> *params);
//virtual T dCS_phase_mod( source_parameters<T> *param);
//virtual T dCS_phase_factor(source_parameters<T> *param);
//
//virtual int construct_amplitude(T *frequencies, int length, T *amplitude, source_parameters<T> *params);
//
////virtual T construct_amplitude(T frequency,  source_parameters<T> *params);
//
//virtual int construct_phase(T *frequencies, int length, T *phase, source_parameters<T> *params);
//};
//
//template<class T> 
//class EdGB_IMRPhenomD_log: public ppE_IMRPhenomD_Inspiral<T>
//{
//public:
//virtual int construct_waveform(T *frequencies, int length, std::complex<T> *waveform, source_parameters<T> *params);
//
//virtual T EdGB_phase_mod( source_parameters<T> *param);
//virtual T EdGB_phase_factor(source_parameters<T> *param);
//
////virtual std::complex<T> construct_waveform(T frequency,source_parameters<T> *params );
//
//virtual int construct_amplitude(T *frequencies, int length, T *amplitude, source_parameters<T> *params);
//
////virtual T construct_amplitude(T frequency,  source_parameters<T> *params);
//
//virtual int construct_phase(T *frequencies, int length, T *phase, source_parameters<T> *params);
//};
//template<class T> 
//class EdGB_IMRPhenomD: public ppE_IMRPhenomD_Inspiral<T>
//{
//public:
//virtual int construct_waveform(T *frequencies, int length, std::complex<T> *waveform, source_parameters<T> *params);
//
//virtual T EdGB_phase_mod( source_parameters<T> *param);
//virtual T EdGB_phase_factor(source_parameters<T> *param);
//
////virtual std::complex<T> construct_waveform(T frequency,source_parameters<T> *params );
//
//virtual int construct_amplitude(T *frequencies, int length, T *amplitude, source_parameters<T> *params);
//
////virtual T construct_amplitude(T frequency,  source_parameters<T> *params);
//
//virtual int construct_phase(T *frequencies, int length, T *phase, source_parameters<T> *params);
//};

#endif
