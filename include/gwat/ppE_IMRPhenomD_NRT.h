#ifndef PPE_IMRPHENOMD_NRT_H
#define PPE_IMRPHENOMD_NRT_H
#include "IMRPhenomD_NRT.h"
#include "util.h"

/*! \file 
 */

/*! Class that extends the IMRPhenomD waveform to include non-GR terms in the inspiral portion of
 * the phase. This is an appropriate waveform choice for generation effects, but not necessarily for
 * propagation effects
 */
 
template<class T> 
class ppE_IMRPhenomD_NRT_Inspiral: public IMRPhenomD_NRT<T>
{
public:
virtual T phase_ins(T f, source_parameters<T> *param, T *pn_coeff, 
		lambda_parameters<T> *lambda, useful_powers<T> *pow);
  
  virtual T phase_ins_NRT(T f, useful_powers<T> *pow, source_parameters<T> *param);

  virtual T Dphase_ins(T f, source_parameters<T> *param, T *pn_coeff, lambda_parameters<T> *lambda);

};
/*! Class that extends the IMRPhenomD_NRT waveform to include non-GR terms in the full
 * phase. This is an appropriate waveform choice for propagation effects
 */

template<class T> 
class ppE_IMRPhenomD_NRT_IMR: public ppE_IMRPhenomD_NRT_Inspiral<T>
{
public:
virtual T Dphase_mr(T f, source_parameters<T> *param, lambda_parameters<T> *lambda);
virtual T phase_mr(T f, source_parameters<T> *param, lambda_parameters<T> *lambda);

virtual T phase_int(T f, source_parameters<T> *param, lambda_parameters<T> *lambda );
virtual T Dphase_int(T f, source_parameters<T> *param, lambda_parameters<T> *lambda);

};

#endif
