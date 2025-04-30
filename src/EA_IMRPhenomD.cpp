#include "EA_IMRPhenomD_NRT.h"
#include "IMRPhenomD_NRT.h"
#include "IMRPhenomD.h"
#include <math.h>
#include <adolc/adouble.h>
#include <adolc/taping.h>
#include <adolc/drivers/drivers.h>
#include <iostream>
#include <cmath>
#include <complex>
#include "util.h"
#include <stdio.h>

/*! \file
 * File for the construction of waveforms in Einstein AEther theory,
 * specifically for neutron stars. Tidal effects are included.
 *
 * Add relevant references.
 */

template<class T>
T EA_IMRPhenomD_NRT<T>::calculate_EA_sensitivity(int body, source_parameters<T> *p)
{
  T lambda, compact, OmRatio, s;
  /* The tidal deformability (love number), compactness, binding energy (Omega)
   * to mass ratio, and sensitivity.
   */
  //std::cout<<"body #"<<body<<std::endl; 
  if(body == 1)
    {
     lambda = p->tidal1;
    }
  else
    {
      lambda = p->tidal2;
    }

  /* Compactness computed using C-Love relation from arXiv:1903.03909,
   * equation 8, values in table 1.
   */
  T a[3] = {-919.6,330.3,-857.2};
  T b[3] = {-383.5,192.5,-811.1};

  // Useful Powers
  T lambda_pow[3];
  lambda_pow[0] = pow(lambda, - 1./5.);
  lambda_pow[1] = lambda_pow[0] * lambda_pow[0];
  lambda_pow[2] = lambda_pow[0] * lambda_pow[1];
  T K, num, denom;
  K = 0.2496;
  num = 1 + a[0]*lambda_pow[0] + a[1]*lambda_pow[1] + a[2]*lambda_pow[2];
  denom = 1 + b[0]*lambda_pow[0] + b[1]*lambda_pow[1] + b[2]*lambda_pow[2];

  //compact = p->compact1; //artifact from testing with specific values of C
  compact = K * lambda_pow[0] * (num/denom);
	//std::cout<<"Compactness: "<<compact<<std::endl;
  
  if(body == 1)
    {
      p->compact1 = compact;
    }
  else
    {
      p->compact2 = compact;
    }
  
  /* Calculation of sensitivities taken from arXiv:2104.04596v1
   * Equation 80 of that paper was inverted to get binding energy to mass ratio
   * as a function of compactness.
   * Equation 81 was used to compute sensitivity as a function of the binding
   * energy to mass ratio.
   */
  OmRatio = (-5./7.)*compact - ((18275.*p->alpha1_EA)/168168.)*pow(compact, 3.);

  //std::cout<<" compactness = "<<compact<<", OmRatio = "<<OmRatio<<std::endl;
  //std::cout<<" compactness = "<<compact<<", OmRatio = "<<OmRatio<<std::endl;
   
  T coeff1, coeff2, coeff3;
   
  coeff1 =  ((3.*p->alpha1_EA + 2.*p->alpha2_EA)/3.);
  coeff2 = ((573.*pow(p->alpha1_EA, 3.) + p->alpha1_EA*p->alpha1_EA*(67669. - 764.*p->alpha2_EA) + 96416.*p->alpha2_EA*p->alpha2_EA + 68.*p->alpha1_EA*p->alpha2_EA*(9.*p->alpha2_EA - 2632.))/(25740.*p->alpha1_EA));
  coeff3 = (1./(656370000.*p->cw_EA*p->alpha1_EA*p->alpha1_EA))*(-4.*p->alpha1_EA*p->alpha1_EA*(p->alpha1_EA + 8.)*(36773030.*p->alpha1_EA*p->alpha1_EA - 39543679.*p->alpha1_EA*p->alpha2_EA + 11403314.*p->alpha2_EA*p->alpha2_EA) + p->cw_EA*(1970100.*pow(p->alpha1_EA,5.) - 13995878400.*pow(p->alpha2_EA, 3.) - 640.*p->alpha1_EA*p->alpha2_EA*p->alpha2_EA*(-49528371. + 345040.*p->alpha2_EA) - 5.*pow(p->alpha1_EA, 4.)*(19548109. + 788040.*p->alpha2_EA) - 16.*p->alpha1_EA*p->alpha1_EA*p->alpha2_EA*(1294533212. - 29152855.*p->alpha2_EA + 212350.*p->alpha2_EA*p->alpha2_EA) + pow(p->alpha1_EA,3.)*(2699192440. - 309701434.*p->alpha2_EA + 5974000.*p->alpha2_EA*p->alpha2_EA)));
  
  s = coeff1 * (OmRatio) + coeff2 * (OmRatio*OmRatio) + coeff3 * (pow(OmRatio, 3.));
  
  /*
  //Coded s(C) directly from eqn.79 of arXiv:2104.04596v1 for comparison purposes. This gives exactly the same answer as what was used above (a good sign obviously)
  coeff1 = (5./21.)*(-3*p->alpha1_EA + 2*p->alpha2_EA);
  coeff2 = (5./(252252.*p->alpha1_EA))*(573.*pow(p->alpha1_EA, 3.) + (67669. - 746.*p->alpha2_EA)*p->alpha1_EA*p->alpha1_EA + 96416.*p->alpha2_EA*p->alpha2_EA + 68.*p->alpha1_EA*p->alpha2_EA*(-2632.+9.*p->alpha2_EA));
  coeff3 = 1./(1801079280.*p->cw_EA*p->alpha1_EA*p->alpha1_EA)*(16.*p->alpha1_EA*p->alpha1_EA*(8+p->alpha1_EA)*(36773030.*p->alpha1_EA*p->alpha1_EA - 39543679.*p->alpha1_EA*p->alpha2_EA + 11403314.*p->alpha2_EA*p->alpha2_EA) + p->cw_EA*(-1970100.*pow(p->alpha1_EA, 5.) + 13995878400.*pow(p->alpha2_EA, 3.) + 640.*p->alpha1_EA*p->alpha2_EA*p->alpha2_EA*(-49528371. + 345040.*p->alpha2_EA) + 5*pow(p->alpha1_EA, 4.)*(-19596941. + 788040.*p->alpha2_EA) + pow(p->alpha1_EA, 3.)*(-2699192440. + 440184934.*p->alpha2_EA - 5974000.*p->alpha2_EA*p->alpha2_EA)*(16.*p->alpha1_EA*p->alpha1_EA*p->alpha2_EA*(1294533212. - 29152855.*p->alpha2_EA + 212350.*p->alpha2_EA*p->alpha2_EA)))); 

  s = coeff1 * (compact) + coeff2 * (compact*compact) + coeff3 * (pow(compact, 3.));
  */
  
  //if(p->c14_EA < pow(10, -32.)){
  /* std::cout<<"c1: "<<p->c1_EA<<", c14: "<<p->c14_EA<<", c13: "<<p->c13_EA<<", cminus: "<<p->cminus_EA<<std::endl;
      
      std::cout<<"alpha1: "<<p->alpha1_EA<<std::endl;
      std::cout<<"alpha2: "<<p->alpha2_EA<<std::endl;
      std::cout<<"Coeff1: "<<coeff1<<std::endl;
      std::cout<<"Coeff2: "<<coeff2<<std::endl;
      std::cout<<"Coeff3: "<<coeff3<<std::endl;
      std::cout<<"tidal1: "<<p->tidal1<<", tidal2: "<<p->tidal2<<std::endl;
      std::cout<<"lambda: "<<lambda<<std::endl; 
      std::cout<<"lambda^(-1/5): "<<lambda_pow[0]<<std::endl; 
      std::cout<<"compact: "<<compact<<std::endl; 
      std::cout<<"OmRatio: "<<OmRatio<<std::endl;
      std::cout<<"s "<<s<<std::endl;*/
      // }
  //if(p->c14_EA < pow(10, -32.)){std::cout<<"s "<<s<<std::endl;}
  /*if(s > 1){
    std::cout<<"sensitivity = "<<s<<", coeff1 = "<<coeff1<<", coeff2 = "<<coeff2<<", coeff3 = "<<coeff3<<std::endl;
    }*/
  return s;
}

template<>
void EA_IMRPhenomD_NRT<double>::EA_check_nan(source_parameters<double> *p)
{
  if(p->EA_nan_error_message)
    {
      if(isnan(p->alpha1_EA))
	{
	  std::cout<<"WARNING: alpha1_EA is NAN"<<std::endl;
	}
      if(isnan(p->alpha2_EA))
	{
	  std::cout<<"WARNING: alpha2_EA is NAN"<<std::endl;
	}
      if(isnan(p->Z_EA))
	{
	  std::cout<<"WARNING: Z_EA is NAN"<<std::endl;
	}
      if(isnan(p->s1_EA))
	{
	  std::cout<<"WARNING: s1_EA is NAN"<<std::endl;
	}
      if(isnan(p->s2_EA))
	{
	  std::cout<<"WARNING: s1_EA is NAN"<<std::endl;
	}
      if(isnan(p->kappa3_EA))
	{
	  std::cout<<"WARNING: kappa3_EA is NAN"<<std::endl;
	}
      if(isnan(p->epsilon_x_EA))
	{
	  std::cout<<"WARNING: epsilon_x_EA is NAN"<<std::endl;
	}
    }
}


template<>
void EA_IMRPhenomD_NRT<adouble>::EA_check_nan(source_parameters<adouble> *p)
{
  if(p->EA_nan_error_message)
    {
      if(isnan(p->alpha1_EA.value()))
	{
	  std::cout<<"WARNING: alpha1_EA is NAN"<<std::endl;
	}
      if(isnan(p->alpha2_EA.value()))
	{
	  std::cout<<"WARNING: alpha2_EA is NAN"<<std::endl;
	}
      if(isnan(p->Z_EA.value()))
	{
	  std::cout<<"WARNING: Z_EA is NAN"<<std::endl;
	}
      if(isnan(p->s1_EA.value()))
	{
	  std::cout<<"WARNING: s1_EA is NAN"<<std::endl;
	}
      if(isnan(p->s2_EA.value()))
	{
	  std::cout<<"WARNING: s2_EA is NAN"<<std::endl;
	}
      if(isnan(p->kappa3_EA.value()))
	{
	  std::cout<<"WARNING: kappa3_EA is NAN"<<std::endl;
	}
      if(isnan(p->epsilon_x_EA.value()))
	{
	  std::cout<<"WARNING: epsilon_x_EA is NAN"<<std::endl;
	}
    }
}

template<class T>
void EA_IMRPhenomD_NRT<T>::pre_calculate_EA_factors(source_parameters<T> *p)
{
  p->V_x_EA = 0.;
  p->V_y_EA = 0.;
  p->V_z_EA = 0.;

  //p->alpha_param = true;

  if(p->EA_region1){
    //std::cout<<"Sampling in region 1 of parameter space"<<std::endl;
    p->ctheta_EA = 3.*p->ca_EA*(1. + p->delta_ctheta_EA); 
    //region of parameter space with ctheta = 3ca. Don't sample on ctheta
  }
  
  //If using alpha parameterization, must first define ca and ctheta in terms
  //of the alphas 
  if(p->alpha_param){
    //std::cout<<"Using alpha parameterization, computing ca & ctheta"<<std::endl;
    p->ca_EA = -p->alpha1_EA/4.;
    p->ctheta_EA = (0.75*p->alpha1_EA*p->alpha1_EA)/(8.*p->alpha2_EA + p->alpha1_EA*(p->alpha2_EA - 0.5*p->alpha1_EA - 1.));
    p->cw_EA = (1 - p->cbarw_EA)/p->cbarw_EA; 
  }
 
  //Transforming to the parameters used in arXiv:1911.10278v2 (because that is where many of these formulas come from)
  p->c1_EA = (p->cw_EA + p->csigma_EA)/2.;
  p->c2_EA = (p->ctheta_EA - p->csigma_EA)/3.;
  p->c3_EA = (p->csigma_EA - p->cw_EA)/2.;
  p->c4_EA = p->ca_EA - (p->csigma_EA + p->cw_EA)/2.;
  //std::cout<<"c1: "<<p->c1_EA<<"c2: "<<p->c2_EA<<"c3: "<<p->c3_EA<<"c4: "<<p->c4_EA<<std::endl; 
  
  //more convenient parameters
  /*p->c13_EA = p->c1_EA + p->c3_EA;
  p->cminus_EA = p->c1_EA - p->c3_EA;
  p->c14_EA = p->c1_EA + p->c4_EA;*/
  p->c13_EA = p->csigma_EA;
  p->c14_EA = p->ca_EA;
  p->cminus_EA = p->cw_EA;

  //squared speeds of the different polarizations
  p->cTsq_EA = 1./(1. - p->c13_EA);
  p->cVsq_EA = (2.*p->c1_EA - p->c13_EA*p->cminus_EA)/(2.*(1.- p->c13_EA)*p->c14_EA);
  p->cSsq_EA = ((2. - p->c14_EA)*(p->c13_EA + p->c2_EA))/((2.+3.*p->c2_EA + p->c13_EA)*(1. - p->c13_EA)*p->c14_EA);

  //speeds of the different polarizations
  p->cT_EA = sqrt(p->cTsq_EA);
  p->cV_EA = sqrt(p->cVsq_EA);
  p->cS_EA = sqrt(p->cSsq_EA);

  //Relevant combinations of parameters
  if(!p->alpha_param){
    //If not using alpha parameterization, alphas should be calculated now
    //If using alpha parameterization, don't want to overwrite them
    //std::cout<<"Not using alpha parameterization, computing alphas"<<std::endl;
    p->alpha1_EA = -8.*(p->c1_EA*p->c14_EA - p->cminus_EA*p->c13_EA)/(2.*p->c1_EA - p->cminus_EA*p->c13_EA);
    p->alpha2_EA = (1./2.)*p->alpha1_EA + ((p->c14_EA - 2.*p->c13_EA)*(3.*p->c2_EA + p->c13_EA + p->c14_EA))/((p->c2_EA + p->c13_EA)*(2. - p->c14_EA));
    p->cbarw_EA = 1./(1. + p->cw_EA); 
  }
  
  p->Z_EA = ((p->alpha1_EA - 2.*p->alpha2_EA)*(1. - p->c13_EA)) / (3.*(2.*p->c13_EA - p->c14_EA));

  p->beta1_EA = -2.* p->c13_EA / p->cV_EA;
  p->beta2_EA = (p->c14_EA - 2.* p->c13_EA)/(2.*p->c14_EA * (1 - p->c13_EA) * p->cS_EA * p->cS_EA);


  p->A1_EA = (1./p->cT_EA) + (2*p->c14_EA*p->c13_EA*p->c13_EA)/((2.*p->c1_EA - p->c13_EA*p->cminus_EA)*(2.*p->c1_EA - p->c13_EA*p->cminus_EA)*p->cV_EA) + (3.*p->c14_EA*(p->Z_EA - 1.)*(p->Z_EA - 1.))/(2.*(2. - p->c14_EA)*p->cS_EA);

  p->A2_EA = -(2.*p->c13_EA)/((2.*p->c1_EA - p->c13_EA*p->cminus_EA)*pow(p->cV_EA, 3.)) - 2.*(p->Z_EA - 1.)/((2. - p->c14_EA)*pow(p->cS_EA, 3.));

  p->A3_EA = 1./(2.*p->c14_EA* pow(p->cV_EA, 5.)) + 2./(3.*p->c14_EA * (2. - p->c14_EA)*pow(p->cS_EA, 5.));

  p->B3_EA = 1./(9.*p->c14_EA*(2. - p->c14_EA)*pow(p->cS_EA, 5.));

  p->C_EA = 4./(3.*p->c14_EA * pow(p->cV_EA, 3.)) + 4./(3.*p->c14_EA*(2. - p->c14_EA)*pow(p->cS_EA, 3.));

  p->D_EA = 1./(6.*p->c14_EA * pow(p->cV_EA, 5.));

  //Get sensitivities
  p->s1_EA = calculate_EA_sensitivity(1, p);
  p->s2_EA = calculate_EA_sensitivity(2, p);
  
  //The functions that are actually used to compute the phase
  p->S_EA = p->s1_EA*(p->mass2/p->M) + p->s2_EA*(p->mass1/p->M);
  p->kappa3_EA = p->A1_EA + p->S_EA * p->A2_EA + p->S_EA*p->S_EA * p->A3_EA;
  p->epsilon_x_EA = (((p->s1_EA - p->s2_EA)*(p->s1_EA - p->s2_EA))/(32.*p->kappa3_EA))*((21.*p->A3_EA + 90.*p->B3_EA + 5.*p->D_EA)*(p->V_x_EA*p->V_x_EA + p->V_y_EA*p->V_y_EA + p->V_z_EA*p->V_z_EA) - (3.*p->A3_EA + 90.*p->B3_EA - 5.*p->D_EA)*p->V_z_EA*p->V_z_EA + 5.*p->C_EA);
  

  EA_check_nan(p);

  //debugger_print(__FILE__,__LINE__,"EA Debugging");
  //std::cout<<"aBL "<<p->abL_EA<<std::endl;
  //std::cout<<"gb1 "<<p->gb1_EA<<std::endl;
  //std::cout<<"gX1 "<<p->gX1_EA<<std::endl;
  //std::cout<<"epsilon_x "<<p->epsilon_x_EA<<std::endl;
  //std::cout<<"S "<<p->S_EA<<std::endl;
  //std::cout<<"alpha "<<p->alpha_ppE_2T_0_EA<<std::endl;
  //std::cout<<"k3 "<<p->kappa3_EA<<std::endl;
  //std::cout<<"cT "<<p->cT_EA<<std::endl;
  //std::cout<<"cV "<<p->cV_EA<<std::endl;
  //std::cout<<"cS "<<p->cS_EA<<std::endl;

    
}
//template void pre_calculate_EA_factors(source_parameters<double> *);
//template void pre_calculate_EA_factors(source_parameters<adouble> *);
//#############################################################


template<class T>
T EA_IMRPhenomD_NRT<T>::EA_phase_ins1(T f, useful_powers<T> *powers, source_parameters<T> *p) {

  T EA_phase, GR_phase, phase_out;

  //did not include the terms that cancel (they are added later in construct_waveform) -- 2*pi*f*t_c - phi(t_c) - pi/4
  EA_phase = (3./128.) * (((1 - p->s1_EA) * (1 - p->s2_EA)) / (2 - p->c14_EA)) * (1 / p->kappa3_EA) * (pow(2, -5./3.) * powers->MFminus_5third * powers->PIminus_5third) * (1 - ((4./7.) * (1/ (pow(2, 2./3.) * powers->MF2third * powers->PI2third)) * pow(p->eta, 2./5.) * p->epsilon_x_EA));
  GR_phase = (3./256.) * (powers->MFminus_5third * powers->PIminus_5third * pow(2, -5./3.));

  phase_out = EA_phase - GR_phase;

  return phase_out;
}


template<class T>
T EA_IMRPhenomD_NRT<T>::EA_phase_ins2(T f, useful_powers<T> *powers, source_parameters<T> *p)
{
  T phase_out;
  T EA_phase, GR_phase;
  /* Here EA_phase is the leading order contribution to the l=2 mode of the
   * Einstein Aether phase (for all polarizations) and GR_phase is the leading
   * order contribution to the l=2 mode of the phase in GR. We will need to
   * subtract GR_phase off so that we are not double counting terms that were
   * already added in IMRPhenomD. Note that I have not included three terms
   * which obviously cancel.
   * These terms are 2Pi f t_c - 2 Phi(t_c) - Pi/4.
   */

  EA_phase = (3./64.)*(((1 - p->s1_EA)*(1 - p->s2_EA))/(2 - p->c14_EA))*(1/p->kappa3_EA)*(powers->MFminus_5third*powers->PIminus_5third)*(1 - (4./7.)*(1./(powers->MF2third*powers->PI2third))*pow(p->eta, 2./5.)*p->epsilon_x_EA);
  GR_phase = (3./128.)*(powers->MFminus_5third*powers->PIminus_5third);

  phase_out = EA_phase - GR_phase; //The correction due to EA theory. This correction applies equally to all polarizations.
  return phase_out;
}

template<class T>
T EA_IMRPhenomD_NRT<T>::EA_amp_ins1(T f, useful_powers<T> *powers, source_parameters<T> *p) {

  T EA_amp;

  EA_amp = (-1./4.) * sqrt(5. * M_PI / 48.) * sqrt((2. - p->c14_EA) / ((1. - p->s1_EA) * (1. - p->s2_EA))) * (1. / p->DL) * (p->s1_EA - p->s2_EA) * p->chirpmass * p->chirpmass * (1. / sqrt(p->kappa3_EA)) * pow(p->eta, 1./5.) * (1. / sqrt(powers->PIcube * powers->MFcube)) * (1. - ((1./2.) * (1. / (pow(2, 2./3.) * powers->PI2third * powers->MF2third)) * pow(p->eta, 2./5.) * p->epsilon_x_EA));

  return EA_amp;
}

template<class T>
//template<>
//double EA_IMRPhenomD_NRT<double>::EA_amp_ins2(double f, useful_powers<double> *powers, source_parameters<double> *p)
T EA_IMRPhenomD_NRT<T>::EA_amp_ins2(T f, useful_powers<T> *powers, source_parameters<T> *p)
{
  /* Here EA_amp is the leading order contribution to the l=2 mode of the
   * Einstein Aether amplitude and GR_amp is the leading order contributio
   * to the l=2 mode of the amplitude in GR. We will need to subtract 
   * GR_amp off so that we are not double counting terms that were
   * already added in IMRPhenomD.
   */
  T EA_amp, GR_amp, amp_out;

  EA_amp = (-1./2.) * sqrt(5. * M_PI / 48) * sqrt((2 - p->c14_EA) / ((1. - p->s1_EA) * (1. - p->s2_EA))) * (1. / p->DL) * p->chirpmass * p->chirpmass * (1. / sqrt(p->kappa3_EA)) * (1. / sqrt(powers->PI7third * powers->MF7third)) * (1. - ((1. / 2.) * (1. / (powers->PI2third * powers->MF2third)) * pow(p->eta, 2./5.) * p->epsilon_x_EA));
  GR_amp = (-1./2.) * sqrt(5. * M_PI / 24) * (1. / p->DL) * p->chirpmass * p->chirpmass * (1. / sqrt(powers->PI7third * powers->MF7third));

  amp_out = EA_amp - GR_amp;
  //EA_IMRPhenomD_NRT<T> model;
  //this->EA_check_nan(true, p);

  return amp_out;
}

template<class T>
int EA_IMRPhenomD_NRT<T>::EA_construct_waveform(T *frequencies, int length, waveform_polarizations<T> *waveform, source_parameters<T> *params)
{
  //this->pre_calculate_EA_factors(params);
  //pre_calculate_EA_factors is now being called in prep_source_parameters and
  //we don't want to call it twice-> don't need to call it here. Note that it is being called before we switch to a barred mass, so the sensitivities (which are calculated in that function) are being calculated with the unbarred version as they should be. And the amplitude and phase of the waveform are being calculated with the barred version of the mass, as they should be because they are calculated after we switch.

  /*TODO*/
  /*The input mass should be unbarred*/
  /*Calcualte sensitivites with unbarred quantities using C = G_N M / R^2 c^2*/
  /*Unbarred to barred */
  //T calG = (1 - params->s1_EA)*(1 - params->s2_EA) ;
  T calG = 1.;
  //std::cout<<"Before transformation:"<<std::endl; 
  //std::cout<<"M = "<<params->M<<", chirpmass = "<<params->chirpmass<<", eta = "<<params->eta<<std::endl; 
  params->M *=calG;
  params->chirpmass *=calG;
  params->delta_mass *=calG;
  params->mass1*=calG;
  params->mass2*=calG;
  // std::cout<<"After transformation:"<<std::endl; 
  //std::cout<<"M = "<<params->M<<", chirpmass = "<<params->chirpmass<<", eta = "<<params->eta<<std::endl;
  
  //std::cout<<"Used EA_construct_waveform. Print statement in line 393 of EA_IMRPhenomD_NRT.cpp"<<std::endl;
  params->NRT_phase_coeff = - (3./16.) * params->tidal_weighted * (39./(16. * params->eta));

  if(params->tidal1<=0) {
    params->oct1=1;params->quad1 = 1;
  }

  else {
    params->quad1 = this->calculate_quad_moment(params->tidal1);
    params->oct1 = this->calculate_oct_moment(params->quad1);
  }

  if(params->tidal2<=0) {
    params->oct2=1;params->quad2 = 1;
  }

  else {
    params->quad2 = this->calculate_quad_moment(params->tidal2);
    params->oct2 = this->calculate_oct_moment(params->quad2);
  }

  this->calculate_spin_coefficients_3p5(params);
  params->NRT_amp_coefficient = this->calculate_NRT_amp_coefficient(params);
  T M = params-> M;
  T chirpmass = params->chirpmass;
  T DL = params->DL;
  lambda_parameters<T> lambda, *lambda_ptr;
  this->assign_lambda_param(params, &lambda);

  /*Initialize the post merger quantities*/
  this->post_merger_variables(params);
  params->f1_phase = 0.018/(params->M);
  params->f2_phase = params->fRD/2.;

  params->f1 = 0.014/(params->M);
  params->f3 = this->fpeak(params, &lambda);

  useful_powers<T> pows;
  this->precalc_powers_PI(&pows);

  T deltas[6];
  T pn_amp_coeffs[7];
  T pn_phase_coeffs[12];

  this->assign_pn_amplitude_coeff(params, pn_amp_coeffs);
  this->assign_static_pn_phase_coeff(params, pn_phase_coeffs);

  this->amp_connection_coeffs(params,&lambda,pn_amp_coeffs,deltas);
  this->phase_connection_coefficients(params,&lambda,pn_phase_coeffs);

  //################################################################
  //Calculate phase and coalescence time variables
  T phic, f_ref, tc, phi_shift, tc_shift;
  /*Note -- to match LALSuite, this sets phiRef equal to phi_phenomD at f_ref,
   * not phenomD_NRT. Arbitrary constant, not really important
   */
  if(params->shift_phase ) {
    f_ref = params->f_ref;
    this->precalc_powers_ins(f_ref, M, &pows);
    phi_shift = (this->build_phase(f_ref,&lambda,params,&pows,pn_phase_coeffs));
    phic = 2*params->phiRef + phi_shift;
  }

  //If phic is specified, ignore f_ref phiRef and use phic
  else {
    f_ref = 0;
    phic = params->phiRef;
  }

  //Assign shift: first shift so coalescence happens at t=0, then shift from there according to tc
  /*Note -- to match LALSuite, this sets phiRef equal to phi_phenomD at f_ref, not phenomD_NRT
 * Arbitrary constant, not really important
 */
  if(params->shift_time) {
    T alpha1_offset = this->assign_lambda_param_element(params,14);
    tc_shift = this->Dphase_mr(params->f3, params, &lambda)+(-lambda.alpha[1]+alpha1_offset)*params->M/params->eta;
  }

  else {
    tc_shift=0;
  }

  tc = 2*M_PI*params->tc + tc_shift;

  T A0 = params->A0* pow(M,7./6.);

  T f;
  std::complex<T> amp, phase;
  std::complex<T> i;
  i = std::complex<T> (0,1.);
  T fcut = .2/M; //Cutoff frequency for IMRPhenomD - all higher frequencies return 0

  for (size_t j =0; j< length; j++) {

    f = frequencies[j];

    if(f>fcut){
      amp = 0.0;
	    waveform->hplus[j] = 0.0;
	    waveform->hcross[j] = 0.0;
	    waveform->hx[j] = 0.0;
	    waveform->hy[j] = 0.0;
	    waveform->hb[j] = 0.0;
	    waveform->hl[j] = 0.0;

    }

    else {
      //if (f<params->f1_phase)
      //This is always needed for NRT
      if (true) {
        this->precalc_powers_ins(f, M, &pows);
      }

      else {
        //pows.MFsixth= pow(M*f,1./6.);
	      //pows.MF7sixth= pow_int(pows.MFsixth,7);//*pows.MFsixth*pows.MFsixth*pows.MFsixth*pows.MFsixth*pows.MFsixth*pows.MFsixth;
	      //Needed for Pade function
	      //pows.MFthird = pows.MFsixth*pows.MFsixth;
	      //pows.MF2third = pows.MFthird*pows.MFthird;
      }

      amp = (A0 * this->build_amp(f,&lambda,params,&pows,pn_amp_coeffs,deltas));
      phase = (this->build_phase(f,&lambda,params,&pows,pn_phase_coeffs));
      /*Append phase_ins_NRT and amp_ins_NRT to the entire waveform*/

      {
        T phaseNRT = this->phase_ins_NRT(f,&pows,params);
        phase += phaseNRT;
        T phaseSpinNRT = this->phase_spin_NRT(f, &pows,params);
        phase += phaseSpinNRT;
        T ampNRT = (A0*this->amp_ins_NRT(f,&pows, params));
        amp +=ampNRT;
      }

      /*Append EA correction to the entire waveform*/
      {
        T EAphase2 = this->EA_phase_ins2(f, &pows, params);
        phase += EAphase2;
        T EAamp2 = (A0*this->EA_amp_ins2(f,&pows, params));
        amp += EAamp2;
      }

      /*Compute coefficients specific to the different polarizations */  
      //These are just the terms for the l=2 mode
      std::complex<T> hx2, hy2, hb2, hl2;
      
      hx2 = (params->beta1_EA /((2.*params->c1_EA - params->c13_EA*params->cminus_EA)*2*params->cV_EA))*(params->S_EA - (params->c13_EA/(1. - params->c13_EA)));
      hy2 = std::complex<T>(2.,0)*hx2*std::complex<T>(0,-1);
      hb2 = (1./(2.-params->c14_EA))*(3.*params->c14_EA*(params->Z_EA - 1.) - (2.*params->S_EA/params->cSsq_EA));
      hl2 = params->abL_EA*hb2;

      //Adding iota dependence 
      std::complex<T> ci = std::complex<T>(cos(params->incl_angle),0);
      std::complex<T> si = std::complex<T>(sin(params->incl_angle),0);
      std::complex<T> s2i = std::complex<T>(sin(2*params->incl_angle),0);

      hx2 *= s2i;
      hy2 *= si;
      hb2 *= si*si;
      hl2 *= si*si;

      //Note that the iota dependence for the plus and cross modes is handled in the waveform_generator.cpp file!
      
      /* The phase correction which depends on the speed of the
	    * different polarizations.
	    */
      T phaseTVScoeff, EAphaseT, EAphaseV, EAphaseS;
	    phaseTVScoeff = 2*M_PI*(f-f_ref)*params->DL;
	    EAphaseT = phaseTVScoeff*(1.-(1./params->cT_EA));
	    EAphaseV = phaseTVScoeff*(1.-(1./params->cV_EA));
	    EAphaseS = phaseTVScoeff*(1.-(1./params->cS_EA));

      phase -= (T)(tc*(f-f_ref) + phic);
      waveform->hplus[j] = amp*std::exp(-i * (phase - EAphaseT));
      waveform->hcross[j] = amp*std::complex<T>(0,-1)* std::exp(-i * (phase - EAphaseT));
      waveform->hx[j] = amp*hx2*std::exp(-i * (phase - EAphaseV));
      waveform->hy[j] = amp*hy2*std::exp(-i * (phase - EAphaseV));
      waveform->hb[j] = amp*hb2*std::exp(-i * (phase - EAphaseS));
      waveform->hl[j] = amp*hl2*std::exp(-i * (phase - EAphaseS));

      if (params->include_l1 == true) {

        std::complex<T> phase1, amp1;
        T f1;

        f1 = f / 2;

        //add the GR phase for l=1 mode
        phase1 = this->build_phase(f1,&lambda,params,&pows,pn_phase_coeffs);

        //calculate and add the NRT phase corrections for the l=1 mode
        T phaseNRT1 = this->phase_ins_NRT(f1,&pows,params);
        phase1 += phaseNRT1;
        T phaseSpinNRT1 = this->phase_spin_NRT(f1, &pows,params);
        phase += phaseSpinNRT1;

        //calculate and add the EA phase corrections for l=1 mode
        phase1 += this->EA_phase_ins1(f, &pows, params);

        //add the EA amp for l=1 mode (there are no NRT or GR corrections to the l=1 amp)
        amp1 = this->EA_amp_ins1(f, &pows, params);

        phase1 -= (T)(tc*(f-f_ref) + phic);

        //compute polarization coefficients (without iota dependence)
        //for hb1, it wouldn't let me multiply by 2i in the same line
        std::complex<T> hx1, hy1, hb1, hl1;
        hy1 = -1. * params->beta1_EA / ((2. * params->c1_EA) - (params->c13_EA * params->cminus_EA));
        hx1 = std::complex<T>(0,1.) * hy1;
        hb1 = 1 / ((2. - params->c14_EA) * params->cS_EA);
        hb1 *= std::complex<T>(0,-2.);
        hl1 = params->abL_EA * hb1;

	hx1 *= ci;	
	hb1 *= si;
	hl1 *= si; 
	

        waveform->hx[j] += (amp1 * hx1 * std::exp(-i * (phase1 - EAphaseV)));
        waveform->hy[j] += (amp1 * hy1 * std::exp(-i * (phase1 - EAphaseV)));
        waveform->hb[j] += (amp1 * hb1 * std::exp(-i * (phase1 - EAphaseS)));
        waveform->hl[j] += (amp1 * hl1 * std::exp(-i * (phase1 - EAphaseS)));
      }

    }

  }

  //###################################### The next part applies the taper.
  for(int i = 0; i<length; i++)
    {
      std::complex<T> Taper = std::complex<T>((T)(1.0) - this->taper(frequencies[i], length, params),0);

      waveform->hplus[i] = waveform->hplus[i] * Taper;
      waveform->hcross[i] = waveform->hcross[i] * Taper;
      waveform->hx[i] = waveform->hx[i] * Taper;
      waveform->hy[i] = waveform->hy[i] * Taper;
      waveform->hb[i] = waveform->hb[i] * Taper;
      waveform->hl[i] = waveform->hl[i] * Taper;
    }

  return 1;
}

/* Below, the construct_phase, construct_amplitude, and
 * construct_waveform functions are overloaded with the earlier
 * versions from parent classes. Because they will not have the
 * correct arguments and we had to define new functions, these
 * simply call the original and warn the user that they have not been
 * updated for Einstein AEther theory.
 */

template<class T>
int EA_IMRPhenomD_NRT<T>::construct_phase(T *frequencies, int length, T *phase, source_parameters<T> *params)
{
  IMRPhenomD<T> model;
  model.construct_phase(frequencies, length, phase, params);
  std::cout<<"WARNING: code is using IMRPhenomD_NRT version of construct_phase function. Not an updated EA version."<<std::endl;
  return 1;
}

template<class T>
int EA_IMRPhenomD_NRT<T>::construct_amplitude(T *frequencies, int length, T *amplitude, source_parameters<T> *params)
{
  IMRPhenomD<T> model;
  model.construct_amplitude(frequencies, length, amplitude, params);
  std::cout<<"WARNING: code is using IMRPhenomD_NRT version of construct_amplitude function. Not an updated EA version."<<std::endl;
  return 1;
}

template<class T>
int EA_IMRPhenomD_NRT<T>::construct_waveform(T *frequencies, int length, std::complex<T> *waveform, source_parameters<T> *params)
{
  IMRPhenomD_NRT<T> model;
  model.construct_waveform(frequencies, length, waveform, params);
  std::cout<<"WARNING: code is using IMRPhenomD_NRT version of construct_waveform function. Not EA version. Search code for call to construct_waveform function and change to EA_construct_waveform."<<std::endl;
  return 1;
}

template class EA_IMRPhenomD_NRT<double>;
template class EA_IMRPhenomD_NRT<adouble>;
