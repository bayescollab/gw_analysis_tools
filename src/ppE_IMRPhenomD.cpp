#include "ppE_IMRPhenomD.h"
#include <math.h>
#include <adolc/adouble.h>
#include <adolc/taping.h>
#include <adolc/drivers/drivers.h>
#include <iostream>
#include <cmath>
#include <complex>
#include "util.h"

/*! \file 
 * File for the implementation of the ppE formalism for testing GR
 *
 * Extends the IMRPhenomD template to include non-GR phase terms
 *
 * Supported waveforms: ppE Inspiral, ppE IMR, dCS, EdGB
 */

//################################################################################
//################################################################################

/*!\brief Overloaded method for the inspiral portion of the phase
 */
template<class T>
T ppE_IMRPhenomD_Inspiral<T>::phase_ins(T f, source_parameters<T> *param, T *pn_coeff, 
		lambda_parameters<T> *lambda, useful_powers<T> *powers)
{
	IMRPhenomD<T> model;
	T PIMFcube = pow(M_PI * param->chirpmass * f, 1./3.);
	T gr_ins = model.phase_ins(f, param, pn_coeff, lambda,powers);
	T phaseout= gr_ins;
	
	for(int i = 0; i<param->Nmod; i++)
		phaseout =phaseout +  pow((PIMFcube),param->bppe[i]) * param->betappe[i];

	return phaseout;

}
template<class T>
T ppE_IMRPhenomD_Inspiral<T>::Dphase_ins(T f, source_parameters<T> *param, T *pn_coeff, 
		lambda_parameters<T> *lambda)
{
	IMRPhenomD<T> model;
	T gr_ins = model.Dphase_ins(f, param, pn_coeff,lambda);
	T phaseout= gr_ins;
	for(int i = 0; i<param->Nmod; i++)
		phaseout =phaseout +  ((T)param->bppe[i]) / 3. * pow( f ,(((T)param->bppe[i])/3.-1.))* pow((param->chirpmass *M_PI ),((T)param->bppe[i])/3.) * param->betappe[i];
	return phaseout;

}
//################################################################################
//################################################################################

template<class T>
T ppE_IMRPhenomD_IMR<T>::phase_mr(T f, source_parameters<T> *param, lambda_parameters<T> *lambda)
{
	ppE_IMRPhenomD_Inspiral<T> model;
	T gr_mr = model.phase_mr(f, param, lambda);
	T phaseout = gr_mr;
	T PIMFcube = pow(M_PI * param->chirpmass * f, 1./3.);
	for (int i = 0 ; i<param->Nmod; i++)
		phaseout+= pow((PIMFcube),param->bppe[i]) * param->betappe[i];
	return phaseout;

}
template<class T>
T ppE_IMRPhenomD_IMR<T>::Dphase_mr(T f, source_parameters<T> *param, lambda_parameters<T> *lambda)
{
	ppE_IMRPhenomD_Inspiral<T> model;
	T Dgr_mr = model.Dphase_mr(f, param, lambda);
	T phaseout = Dgr_mr;
	for (int i = 0 ; i<param->Nmod; i++)
		phaseout+= param->bppe[i] / 3. * pow( f ,((T)param->bppe[i]/3.-1))* pow((param->chirpmass *M_PI ),(T)param->bppe[i]/3.) * param->betappe[i];
	return phaseout;

}
template<class T>
T ppE_IMRPhenomD_IMR<T>::Dphase_int(T f, source_parameters<T> *param, lambda_parameters<T> *lambda)
{
	ppE_IMRPhenomD_Inspiral<T> model;
	T Dgr_int = model.Dphase_int(f, param, lambda);
	T phaseout = Dgr_int;
	for (int i = 0 ; i<param->Nmod; i++)
		phaseout+= param->bppe[i] / 3. * pow( f ,((T)param->bppe[i]/3.-1))* pow((param->chirpmass *M_PI ),(T)param->bppe[i]/3.) * param->betappe[i];
	return phaseout;

}
template<class T>
T ppE_IMRPhenomD_IMR<T>::phase_int(T f, source_parameters<T> *param, lambda_parameters<T> *lambda)
{
	ppE_IMRPhenomD_Inspiral<T> model;
	T gr_int = model.phase_int(f, param, lambda);
	T phaseout = gr_int;
	T PIMFcube = pow(M_PI * param->chirpmass * f, 1./3.);
	for (int i = 0 ; i<param->Nmod; i++)
		phaseout+= pow((PIMFcube),(T)param->bppe[i]) * param->betappe[i];
	return phaseout;

}


template class ppE_IMRPhenomD_Inspiral<double>;
template class ppE_IMRPhenomD_Inspiral<adouble>;
template class ppE_IMRPhenomD_IMR<double>;
template class ppE_IMRPhenomD_IMR<adouble>;

