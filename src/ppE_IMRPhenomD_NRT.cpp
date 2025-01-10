#include "ppE_IMRPhenomD_NRT.h"
#include <math.h>
#include <adolc/adouble.h>
#include <adolc/taping.h>
#include <adolc/drivers/drivers.h>
#include <iostream>
#include <cmath>
#include <complex>
#include "util.h"

/*!\brief Overloaded method for the inspiral portion of the phase
 */
template<class T>
T ppE_IMRPhenomD_NRT_Inspiral<T>::phase_ins(T f, source_parameters<T> *param, T *pn_coeff, 
		lambda_parameters<T> *lambda, useful_powers<T> *powers)
{
	IMRPhenomD_NRT<T> model;
	T PIMFcube = pow(M_PI * param->chirpmass * f, 1./3.);
	T gr_ins = model.phase_ins(f, param, pn_coeff, lambda,powers);
	T phaseout= gr_ins;

	//std::cout<<"Entering phase_ins function"<<__FILE__<<__LINE__<<std::endl;
	//return the model itself, no ppE corrections test
	//return phaseout;
	
	for(int i = 0; i<param->Nmod; i++)
		phaseout =phaseout +  pow((PIMFcube),param->bppe[i]) * param->betappe[i];

	return phaseout;

}
template<class T>
T ppE_IMRPhenomD_NRT_Inspiral<T>::phase_ins_NRT(T f, useful_powers<T> *powers, source_parameters<T> *param)
{
	IMRPhenomD_NRT<T> model;
	T PIMFcube = pow(M_PI * param->chirpmass * f, 1./3.);
	T gr_ins = model.phase_ins_NRT(f, powers, param);
	T phaseout= gr_ins;

	//std::cout<<"Entering phase_ins_NRT function"<<__FILE__<<__LINE__<<std::endl;
	//return the model itself, no ppE corrections test
	return phaseout;
	/*
	for(int i = 0; i<param->Nmod; i++)
		phaseout =phaseout +  pow((PIMFcube),param->bppe[i]) * param->betappe[i];

	return phaseout;
	*/
}

template<class T>
T ppE_IMRPhenomD_NRT_Inspiral<T>::Dphase_ins(T f, source_parameters<T> *param, T *pn_coeff, 
		lambda_parameters<T> *lambda)
{
	IMRPhenomD_NRT<T> model;
	T gr_ins = model.Dphase_ins(f, param, pn_coeff,lambda);
	T phaseout= gr_ins;

	//return the model itself, no ppE corrections test
	//return phaseout;
	
	for(int i = 0; i<param->Nmod; i++)
		phaseout =phaseout +  ((T)param->bppe[i]) / 3. * pow( f ,(((T)param->bppe[i])/3.-1.))* pow((param->chirpmass *M_PI ),((T)param->bppe[i])/3.) * param->betappe[i];
	return phaseout;
	

}
//################################################################################
//################################################################################

template<class T>
T ppE_IMRPhenomD_NRT_IMR<T>::phase_mr(T f, source_parameters<T> *param, lambda_parameters<T> *lambda)
{
	

	ppE_IMRPhenomD_NRT_Inspiral<T> model;
	T gr_mr = model.phase_mr(f, param, lambda);
	T phaseout = gr_mr;

	//return the model itself, no ppE corrections test
	//return phaseout;
	
	T PIMFcube = pow(M_PI * param->chirpmass * f, 1./3.);
	for (int i = 0 ; i<param->Nmod; i++)
		phaseout+= pow((PIMFcube),param->bppe[i]) * param->betappe[i];
	return phaseout;
	
}
template<class T>
T ppE_IMRPhenomD_NRT_IMR<T>::Dphase_mr(T f, source_parameters<T> *param, lambda_parameters<T> *lambda)
{
	

	ppE_IMRPhenomD_NRT_Inspiral<T> model;
	T Dgr_mr = model.Dphase_mr(f, param, lambda);
	T phaseout = Dgr_mr;

	//return the model itself, no ppE corrections test
	//return phaseout;
	
	for (int i = 0 ; i<param->Nmod; i++)
		phaseout+= param->bppe[i] / 3. * pow( f ,((T)param->bppe[i]/3.-1))* pow((param->chirpmass *M_PI ),(T)param->bppe[i]/3.) * param->betappe[i];
	return phaseout;
	
}
template<class T>
T ppE_IMRPhenomD_NRT_IMR<T>::Dphase_int(T f, source_parameters<T> *param, lambda_parameters<T> *lambda)
{
	

	ppE_IMRPhenomD_NRT_Inspiral<T> model;
	T Dgr_int = model.Dphase_int(f, param, lambda);
	T phaseout = Dgr_int;

	//return the model itself, no ppE corrections test
	//return phaseout;
	
	for (int i = 0 ; i<param->Nmod; i++)
		phaseout+= param->bppe[i] / 3. * pow( f ,((T)param->bppe[i]/3.-1))* pow((param->chirpmass *M_PI ),(T)param->bppe[i]/3.) * param->betappe[i];
	return phaseout;
	
}
template<class T>
T ppE_IMRPhenomD_NRT_IMR<T>::phase_int(T f, source_parameters<T> *param, lambda_parameters<T> *lambda)
{
	

	ppE_IMRPhenomD_NRT_Inspiral<T> model;
	T gr_int = model.phase_int(f, param, lambda);
	T phaseout = gr_int;

	//return the model itself, no ppE corrections test
	//return phaseout;
	
	T PIMFcube = pow(M_PI * param->chirpmass * f, 1./3.);
	for (int i = 0 ; i<param->Nmod; i++)
		phaseout+= pow((PIMFcube),(T)param->bppe[i]) * param->betappe[i];
	return phaseout;
	
}


template class ppE_IMRPhenomD_NRT_Inspiral<double>;
template class ppE_IMRPhenomD_NRT_Inspiral<adouble>;
template class ppE_IMRPhenomD_NRT_IMR<double>;
template class ppE_IMRPhenomD_NRT_IMR<adouble>;
