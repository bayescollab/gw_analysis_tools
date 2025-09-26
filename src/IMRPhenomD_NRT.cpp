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
#include <gsl/gsl_randist.h>

/*! \file 
 * File for the addition of tidal effects to the waveform. 
 *
 * Extends the IMRPhenomD template to include tidal effects, following the 
 * NRTidal model of arXiv:1905.06011v2. Specifically, equations 17, 18, 19, 20, 21, 24
 * and 25 were used. A Planck taper was added to the waveform following arXiv:1003.2939.
 * 
 * Also included in this file
 *  1. Binary love relations 
 *  2. Error marginalization over different equations of state
 *  3. Dissipative tidal response (see arXiv:2306.15633, equations 12b, 13b, 42)  
 */

double tidal_error(double tidal_s, double tidal_a, double q){
  /* Performing error marginalization over residual EoS 
   * dependence of the binary Love relations (following 
   * equations 15-22 of arXiv:1903.03909). The relevant 
   * coefficients/fit parameters are in IMRPhenomD_NRT.h. 
   */
  double error, error_mean, error_standardDev, sigma_L, sigma_q;
  double lambda_pow_new[3];
  double q_sq = q*q;
  
  lambda_pow_new[0] = pow(tidal_s, 1./2.);
  lambda_pow_new[1] = lambda_pow_new[0]*tidal_s;
  lambda_pow_new[2] = lambda_pow_new[1]*tidal_s; 
  //error_mean = (mu_binLove[0]*tidal_s + mu_binLove[1] + mu_binLove[2]*q*q + mu_binLove[3]*q + mu_binLove[4])/2.; //definition given in arXiv:1903.03909
  /*
    sigma_L = sigma_binLove[0]*lambda_pow_new[2] + sigma_binLove[1]*lambda_pow_new[1] + sigma_binLove[2]*tidal_s + sigma_binLove[3]*lambda_pow_new[0] + sigma_binLove[4];
    sigma_q = sigma_binLove[5]*q*q*q + sigma_binLove[6]*q*q + sigma_binLove[7]*q + sigma_binLove[8]; 
    error_standardDev = sqrt(sigma_L*sigma_L + sigma_q*sigma_q); 
    //Definition given in arXiv:1903.03909
    */

  //Using the new formulation and coefficients we found
  error_mean = bL_error[0] + bL_error[1]*q + bL_error[2]*q_sq + bL_error[3]*tidal_s + bL_error[4]*q*tidal_s + bL_error[5]*q_sq*tidal_s + bL_error[6]*q*tidal_s*tidal_s;
  sigma_q = bL_error[7] + bL_error[8]*q + bL_error[9]*q_sq + bL_error[10]*q*q_sq;
  sigma_L = bL_error[11]*lambda_pow_new[2] + bL_error[12]*lambda_pow_new[1] + bL_error[13]*tidal_s + bL_error[14]*lambda_pow_new[0];
  error_standardDev = sigma_q + sigma_L + bL_error[15]*lambda_pow_new[0]*q + bL_error[16]*q_sq*lambda_pow_new[0] + bL_error[17]*q*tidal_s;
		  
  //Random number declaration and seeding
  const gsl_rng_type *t;
  gsl_rng *r;
  
  double seed = omp_get_wtime();
  gsl_rng_env_setup();
  
  t=gsl_rng_default;
  r=gsl_rng_alloc(t);
  gsl_rng_set(r, ( seed-( (int)seed) )*pow(10, 6.) ); //seeding the random number generator with time
		  
  //Now get a random point from the normal distribution with mean (error_mean) and variance (error_variance).
  //This will be added to the binary love relation.
		  
  error = error_mean + gsl_ran_gaussian(r, error_standardDev);
  gsl_rng_free(r);
  
  tidal_a += error;
  //std::cout<<"Error flag is set correctly and tidal love error is being calculated."<<std::endl; //print for testing purposes
	       
  return tidal_a;
}

adouble tidal_error(adouble tidal_s, adouble tidal_a, adouble q)
{
  std::cout<<"Note that error marginalization over the EoS for the binary love relations is not supported for Fishers and is not being performed."<<std::endl; 
  return tidal_a;
}

template<class T>
void IMRPhenomD_NRT<T>::binary_love_relation(T tidal_s, bool tidal_love_error, source_parameters<T> *sp)
{
  /* The binary love relations are used here to compute lambda_a
   * as a function of lambda_s (following equations 11-13 of
   * arXiv:1903.03909). These relations were fit for Neutron stars,
   * so the relevant coefficients/fit parameters are in
   * IMRPhenomD_NRT.h.
   */
  T tidal_a = -1;
      
  T q, Q, F;
  q = sp->mass2 / sp->mass1;
  Q = pow(q, 10./(3. - n_binLove));
  F = (1. - Q)/(1. + Q);
      
  T num = 1;
  T denom = 1;
  T q_pow[2], lambda_pow[3];
  q_pow[0] = q;
  q_pow[1] = q*q;
  lambda_pow[0] = pow(tidal_s, -1./5.);
  lambda_pow[1] = lambda_pow[0] * lambda_pow[0];
  lambda_pow[2] = lambda_pow[0] * lambda_pow[1];
  for(int i = 0; i<3; i++)
    {
      for(int j = 0; j<2; j++)
	{
	  num += b_binLove[i][j]*q_pow[j]*lambda_pow[i];
	  denom += c_binLove[i][j]*q_pow[j]*lambda_pow[i];
	}
    }
  
  tidal_a = F * (num / denom) * tidal_s;
  
  if(tidal_love_error)
    {
      tidal_a = tidal_error(tidal_s, tidal_a, q);
    }
  /* This matches the definition in arXiv:1512.02639, but that paper used the opposite mass convention (m1 < m2) from our code.
  sp->tidal1 = tidal_s + tidal_a;
  sp->tidal2 = tidal_s - tidal_a;
  */
  // This matches the definition in arXiv:1903.03909 which is what we are following here (is consistent with all our conventions)
  sp->tidal1 = tidal_s - tidal_a;
  sp->tidal2 = tidal_s + tidal_a;
}

template<class T>
void IMRPhenomD_NRT<T>::assign_static_pn_phase_coeff(source_parameters<T> *source_param, T *coeff)
{
	IMRPhenomD<T> base_model;
	base_model.assign_static_pn_phase_coeff(source_param,coeff);

  T ssA_2PN, ssB_2PN, ss_2PN, ssA_3PN, ssB_3PN, ss_3PN, ssA_3p5PN, ssB_3p5PN, ss_3p5PN;
  T X_A, X_Asq, chi1, chi1_sq, lambda1, quad1, oct1; 
  T X_B, X_Bsq, chi2, chi2_sq, lambda2, quad2, oct2;  

  lambda1 = source_param->tidal1;
  lambda2 = source_param->tidal2;
  
  X_A = source_param->mass1 / source_param->M;    
  X_B = source_param->mass2 / source_param->M;
  chi1 = source_param->spin1z;
  chi2 = source_param->spin2z;
  
  X_Asq = X_A * X_A;
  X_Bsq = X_B * X_B; 
  chi1_sq = chi1 * chi1;
  chi2_sq = chi2 * chi2;


  quad1 = source_param->quad1;
  quad2 = source_param->quad2;
  oct1 = source_param->oct1;
  oct2 = source_param->oct2;
  
  /*Following equation 27 of NRTidal paper (arXiv:1905.06011v2)*/
  //2 PN contribution
  ssA_2PN = -50*(quad1 - 1.) * X_Asq * chi1_sq;
  ssB_2PN = -50*(quad2 - 1.) * X_Bsq * chi2_sq;
  ss_2PN = ssA_2PN + ssB_2PN;

  //3 PN contribution 
  ssA_3PN = (5/84.) * (9407 + 8218 * X_A - 2016 * X_Asq)* (quad1 - 1.) * X_Asq * chi1_sq;
  ssB_3PN = (5/84.) * (9407 + 8218 * X_B - 2016 * X_Bsq)* (quad2 - 1.) * X_Bsq * chi2_sq;
  ss_3PN = ssA_3PN + ssB_3PN;

  coeff[4]+= ss_2PN;
  coeff[10]+= ss_3PN;
}

//##############################################################################
//##############################################################################
template<class T>
T IMRPhenomD_NRT<T>::Pade(T f, source_parameters<T> *param,useful_powers<T> *powers, char deriv)
{ 
  //T x = pow((M_PI * param->M *f), 2./3.);
  
  //T xpowers[5];
  //for(int i = 0; i<5; i++)
  //  {
  //    xpowers[i] = pow(x, 1 + i/2.);
  //  }
  //
  //  Same as above, but faster because now call to pow -- just need to define powers->MF2third in 
  //  construct waveform because precalc_powers isn't called after inspiral is over
  T x = powers->MF2third * powers->PI2third;
  T x_3_2 = param->M*f*M_PI;
  T xpowers[5];
  xpowers[0] = x;
  xpowers[1] = x_3_2;
  xpowers[2] = x*x;
  xpowers[3] = x*x_3_2;
  xpowers[4] = x*x*x;
	
  
  T P_NRT = 1, P_NRTdenom = 1;
  for(int i = 0; i<5; i++)
    {
      P_NRT += (n_NRT[i] * xpowers[i]);
      if(i<3) P_NRTdenom += (d_NRT[i] * xpowers[i]); 
    }

  if(!deriv)
    {
      P_NRT = xpowers[3] * P_NRT/P_NRTdenom; //Note that this is the Pade function times x^(5/2).
      //std::cout<<P_NRT<<","<<std::endl; 
      
      return P_NRT;
    }
  else
    {
      T DP_NRT;
      T term1, term2;

      term1 = n_NRT[0] + (3./2.)*n_NRT[1]*pow(x,1./2.) + 2*n_NRT[2]*xpowers[0] + (5./2.)*n_NRT[3]*xpowers[1] + 3*n_NRT[4]*xpowers[2]; 
      term2 = d_NRT[0] + (3./2.)*d_NRT[1]*pow(x,1./2.) + 2*d_NRT[2]*xpowers[0];
  
      DP_NRT = xpowers[3]*term1/P_NRTdenom + xpowers[3]* term2*P_NRT/(P_NRTdenom*P_NRTdenom) + (5./2.)*xpowers[1]* (P_NRT/P_NRTdenom);
      //This is a derivative with respect to x. To get a derivative with respect to f, need to multiply by dx/df
      DP_NRT = DP_NRT* (2/3.)*pow((param->M*param->M* M_PI*M_PI/f), 1/3.);

      /*With the current structure of the code, this portion of this function will never be needed since
	we're not actually calculating the derivative of the phase to do matching between inspiral and 
	intermediate regions. */
      return DP_NRT;
    }
}

/* Note that this is NOT an overloaded version of phase_ins -- completely new. Only calculates the NRT phase, not the total phase_ins. It will be appended to the entire waveform later
 */
template<class T> 
T IMRPhenomD_NRT<T>::phase_ins_NRT(T f, useful_powers<T> *powers,source_parameters<T> *param)
{
  
  /*Note that this is just the NRT part now*/
  T phaseout;
  bool deriv = false; //tells the Pade function not to take a derivative

  //phaseout = - ((3./16.) * param->tidal_weighted * (39./(16. * param->eta)) * Pade(f, param, powers,deriv));
  phaseout = param->NRT_phase_coeff * Pade(f, param, powers,deriv);
  
  return phaseout;
}

/* This adds the dissipative tidal correction. There is no IMRPhenom correction--just the PN correction.
 * See 2306.15633 
 */
template<class T> 
T IMRPhenomD_NRT<T>::phase_ins_NRT_D(T f, useful_powers<T> *powers, source_parameters<T> *param)
{
  
  /*Note that this is just the NRT part now*/
  T phaseout;
  
  T piMf = M_PI * param->M * f;

  //phaseout = - ((3./16.) * param->diss_tidal_weighted * (39./(16. * param->eta)) * Pade(f, param, powers,deriv));
  phaseout = - (75./512.) * (1.0/param->eta) * param->diss_tidal_weighted * piMf * log(piMf) ;
  
  return phaseout;
}

//Old version...
/*template<class T>
T IMRPhenomD_NRT<T>::spin_spin(source_parameters<T> *param, double PNorder)
{
  //T ss, X, Xsq, chi, chisq, X_other, chi_other, lambda, quad, oct;
  T ssA_2PN, ssB_2PN, ss_2PN, ssA_3PN, ssB_3PN, ss_3PN, ssA_3p5PN, ssB_3p5PN, ss_3p5PN;
  T X_A, X_Asq, chi1, chi1_sq, lambda1, quad1, oct1; 
  T X_B, X_Bsq, chi2, chi2_sq, lambda2, quad2, oct2;  

  lambda1 = param->tidal1;
  lambda2 = param->tidal1;
  
  X_A = param->mass1 / param->M;    
  X_B = param->mass2 / param->M;
  chi1 = param->spin1z;
  chi2 = param->spin2z;
  
  X_Asq = X_A * X_A;
  X_Bsq = X_B * X_B; 
  chi1_sq = chi1 * chi1;
  chi2_sq = chi2 * chi2;
  
      /*  if(body == 1)
    {
      lambda = param->tidal1;
      X = param->mass1 / param->M;    
      X_other = param->mass2 / param->M;
      chi = param->spin1z;
      chi_other = param->spin2z; 

      //chi = sqrt(param->spin1x*param->spin1x + param->spin1y*param->spin1y + param->spin1z*param->spin1z);
      //chi_other = sqrt(param->spin2x*param->spin2x + param->spin2y*param->spin2y + param->spin2z*param->spin2z);
    }
  else
    {
      lambda = param->tidal2;
      X = param->mass2 / param->M;
      X_other = param->mass1 / param->M;
      chi = param->spin2z;
      chi_other = param->spin1z; 
      
      //chi = sqrt(param->spin2x*param->spin2x + param->spin2y*param->spin2y + param->spin2z*param->spin2z);
      //chi_other = sqrt(param->spin1x*param->spin1x + param->spin1y*param->spin1y + param->spin1z*param->spin1z);
    }
      */

    /*Numerical coefficients from tables 1 and 2 of arXiv:1608.02582*/
/*  T q0 = 0.1940;
  T q1 = 0.09163;
  T q2 = 0.04812;
  T q3 = -0.004286;
  T q4 = 0.00012450;

  T o0 = 0.003131;
  T o1 = 2.071;
  T o2 = -0.7152;
  T o3 = 0.2458;
  T o4 = -0.03309;
*/
  /*Equation 15 of arXiv:1608.02582 for quadrupolar and octupolar spin induced deformabilities*/
/*  quad1 = exp(q0 + q1*log(lambda1) + q2*log(lambda1)*log(lambda1) + q3*log(lambda1)*log(lambda1)*log(lambda1) + q4*log(lambda1)*log(lambda1)*log(lambda1)*log(lambda1));
  oct = exp(o0 + o1*log(quad1) + o2*log(quad1)*log(quad1) + o3*log(quad1)*log(quad1)*log(quad1) + o4*log(quad1)*log(quad1)*log(quad1)*log(quad1));
  
  quad2 = exp(q0 + q1*log(lambda2) + q2*log(lambda2)*log(lambda2) + q3*log(lambda2)*log(lambda2)*log(lambda2) + q4*log(lambda2)*log(lambda2)*log(lambda2)*log(lambda2));
  oct2 = exp(o0 + o1*log(quad2) + o2*log(quad2)*log(quad2) + o3*log(quad2)*log(quad2)*log(quad2) + o4*log(quad2)*log(quad2)*log(quad2)*log(quad2)); 
*/
  //quad = exp(q0) * pow(lambda, q1 + q2 * log(lambda) + q3 * pow(log(lambda),2) + q4 * pow(log(lambda),3));
  //oct = exp(o0) * pow(quad, o1 + o2 * log(quad) + o3 * pow(log(quad),2) + o4 * pow(log(quad),3));
/*  
  if(PNorder == 2.0)
    {
      //ss_2PN = 0.0;
      ssA_2PN = -50*(quad1 - 1.) * X_Asq * chi1_sq;
      ssB_2PN = -50*(quad2 - 1.) * X_Bsq * chi2_sq;
      ss_2PN = ssA_2PN + ssB_2PN;
    }
  else if(PNorder == 3.0)
    {
      //ss = 0.0;
      ssA_3PN = (5/84.) * (9407 + 8218 * X_A - 2016 * X_Asq)* (quad1 - 1.) * X_Asq * chi1_sq;
      ssB_3PN = (5/84.) * (9407 + 8218 * X_B - 2016 * X_Bsq)* (quad2 - 1.) * X_Bsq * chi2_sq;
      ss_3PN = ssA_3PN + ssB_3PN; 
    }
  else
    {
      ssA_3p5PN = 10 * ((X_Asq + (308./3.)*X_A)*chi1 + (X_Bsq - (89/3.)*X_B)*chi2 - 40* M_PI)*(quad1 - 1.)*X_Asq*chi1_sq - 440*(oct1 - 1.)*X_Asq*X_A*chi1_sq*chi1; 
      ssB_3p5PN = 10 * ((X_Bsq + (308./3.)*X_B)*chi2 + (X_Asq - (89/3.)*X_A)*chi1 - 40* M_PI)*(quad2 - 1.)*X_Bsq*chi2_sq - 440*(oct2 - 1.)*X_Bsq*X_B*chi2_sq*chi2;
      ss_3p5PN = ssA_3p5PN + ssB_3p5PN; 
    }
  
  return ss; 
}*/

template<class T>
T IMRPhenomD_NRT<T>::calculate_quad_moment(T lambda)
{
  /*Numerical coefficients from tables 1 and 2 of arXiv:1608.02582*/
  T q0 = 0.1940;
  T q1 = 0.09163;
  T q2 = 0.04812;
  //T q3 = -0.004286; //This matches LAL. Ask Nico about which one we should use
  T q3 = -0.004283; //note small discrepancy between this number in arXiv:1608.02582 vs arXiv:1905.06011v2 (last digit is different). 
  T q4 = 0.00012450;
  T quad = exp(q0 + q1*log(lambda) + q2*pow(log(lambda), 2.) + q3*pow(log(lambda), 3.) + q4*pow(log(lambda), 4.));
  return quad;

}
template<class T>
T IMRPhenomD_NRT<T>::calculate_oct_moment(T quad_moment)
{
  T o0 = 0.003131;
  T o1 = 2.071;
  T o2 = -0.7152;
  T o3 = 0.2458;
  T o4 = -0.03309;
  /*Equation 15 of arXiv:1608.02582 for quadrupolar and octupolar spin induced deformabilities. Note that this only works for lambda >= 1*/
  T oct = exp(o0 + o1*log(quad_moment) + o2*pow(log(quad_moment), 2.) + o3*pow(log(quad_moment), 3.) + o4*pow(log(quad_moment), 4.));
  return oct;
}

template<class T>
void IMRPhenomD_NRT<T>::calculate_spin_coefficients_3p5(source_parameters<T> *param)
{
  T ssA_3p5PN, ssB_3p5PN;
  T X_A, X_Asq, chi1, chi1_sq, lambda1, quad1, oct1; 
  T X_B, X_Bsq, chi2, chi2_sq, lambda2, quad2, oct2;  
  
  X_A = param->mass1 / param->M;    
  X_B = param->mass2 / param->M;
  chi1 = param->spin1z;
  chi2 = param->spin2z;
  
  X_Asq = X_A * X_A;
  X_Bsq = X_B * X_B; 
  chi1_sq = chi1 * chi1;
  chi2_sq = chi2 * chi2;
	
  quad1 = param->quad1;
  quad2 = param->quad2;
  oct1 = param->oct1;
  oct2 = param->oct2;

  ssA_3p5PN = 10 * ((X_Asq + (308./3.)*X_A)*chi1 + (X_Bsq - (89/3.)*X_B)*chi2 - 40* M_PI)*(quad1 - 1.)*X_Asq*chi1_sq - 440*(oct1 - 1.)*X_Asq*X_A*chi1_sq*chi1; 
  ssB_3p5PN = 10 * ((X_Bsq + (308./3.)*X_B)*chi2 + (X_Asq - (89/3.)*X_A)*chi1 - 40* M_PI)*(quad2 - 1.)*X_Bsq*chi2_sq - 440*(oct2 - 1.)*X_Bsq*X_B*chi2_sq*chi2;
  param->ss_3p5PN_coeff = ssA_3p5PN + ssB_3p5PN;

  return ;

}

template<class T>
T IMRPhenomD_NRT<T>::phase_spin_NRT(T f, useful_powers<T> *powers,source_parameters<T> *param)
{
  //###########################################################################
  //T ssA_2PN, ssB_2PN, ss_2PN, ssA_3PN, ssB_3PN, ss_3PN, ssA_3p5PN, ssB_3p5PN, ss_3p5PN;
  //T X_A, X_Asq, chi1, chi1_sq, lambda1, quad1, oct1; 
  //T X_B, X_Bsq, chi2, chi2_sq, lambda2, quad2, oct2;  

  //lambda1 = param->tidal1;
  //lambda2 = param->tidal2;
  //
  //X_A = param->mass1 / param->M;    
  //X_B = param->mass2 / param->M;
  //chi1 = param->spin1z;
  //chi2 = param->spin2z;
  //
  //X_Asq = X_A * X_A;
  //X_Bsq = X_B * X_B; 
  //chi1_sq = chi1 * chi1;
  //chi2_sq = chi2 * chi2;
  //      
  //quad1 = param->quad1;
  //quad2 = param->quad2;
  //oct1 = param->oct1;
  //oct2 = param->oct2;

  ///*Numerical coefficients from tables 1 and 2 of arXiv:1608.02582*/
  //T q0 = 0.1940;
  //T q1 = 0.09163;
  //T q2 = 0.04812;
  ////T q3 = -0.004286; //This matches LAL. Ask Nico about which one we should use
  //T q3 = -0.004283; //note small discrepancy between this number in arXiv:1608.02582 vs arXiv:1905.06011v2 (last digit is different). 
  //T q4 = 0.00012450;

  //T o0 = 0.003131;
  //T o1 = 2.071;
  //T o2 = -0.7152;
  //T o3 = 0.2458;
  //T o4 = -0.03309;
  ///*Equation 15 of arXiv:1608.02582 for quadrupolar and octupolar spin induced deformabilities. Note that this only works for lambda >= 1*/
  //if(lambda1<=0){oct1=1;quad1 = 1;}
  //else{
  //quad1 = exp(q0 + q1*log(lambda1) + q2*pow(log(lambda1), 2.) + q3*pow(log(lambda1), 3.) + q4*pow(log(lambda1), 4.));
  //oct1 = exp(o0 + o1*log(quad1) + o2*pow(log(quad1), 2.) + o3*pow(log(quad1), 3.) + o4*pow(log(quad1), 4.));
  //}
  //
  //if(lambda2<=0){oct2=1;quad2 = 1;}
  //else{
  //quad2 = exp(q0 + q1*log(lambda2) + q2*pow(log(lambda2), 2.) + q3*pow(log(lambda2), 3.) + q4*pow(log(lambda2), 4.));
  //oct2 = exp(o0 + o1*log(quad2) + o2*pow(log(quad2), 2.) + o3*pow(log(quad2), 3.) + o4*pow(log(quad2), 4.));
  //}

  //std::cout<<"quad1: "<<quad1<<"\t quad2: "<<quad2<<std::endl; 
  /*Following equation 27 of NRTidal paper (arXiv:1905.06011v2)*/
  //2 PN contribution
  //ssA_2PN = -50*(quad1 - 1.) * X_Asq * chi1_sq;
  //ssB_2PN = -50*(quad2 - 1.) * X_Bsq * chi2_sq;
  //ssA_2PN = -50*(quad1) * X_Asq * chi1_sq;
  //ssB_2PN = -50*(quad2) * X_Bsq * chi2_sq;
  //ss_2PN = ssA_2PN + ssB_2PN;
  //ss_2PN = 0.0;

  //3 PN contribution 
  //ssA_3PN = (5/84.) * (9407 + 8218 * X_A - 2016 * X_Asq)* (quad1 - 1.) * X_Asq * chi1_sq;
  //ssB_3PN = (5/84.) * (9407 + 8218 * X_B - 2016 * X_Bsq)* (quad2 - 1.) * X_Bsq * chi2_sq;
  //ssA_3PN = (5/84.) * (9407 + 8218 * X_A - 2016 * X_Asq)* (quad1 ) * X_Asq * chi1_sq;
  //ssB_3PN = (5/84.) * (9407 + 8218 * X_B - 2016 * X_Bsq)* (quad2 ) * X_Bsq * chi2_sq;
  //ss_3PN = ssA_3PN + ssB_3PN;
  //ss_3PN = 0.0;

  //3.5 PN contribution
  //ssA_3p5PN = 10 * ((X_Asq + (308./3.)*X_A)*chi1 + (X_Bsq - (89/3.)*X_B)*chi2 - 40* M_PI)*(quad1 - 1.)*X_Asq*chi1_sq - 440*(oct1 - 1.)*X_Asq*X_A*chi1_sq*chi1; 
  //ssB_3p5PN = 10 * ((X_Bsq + (308./3.)*X_B)*chi2 + (X_Asq - (89/3.)*X_A)*chi1 - 40* M_PI)*(quad2 - 1.)*X_Bsq*chi2_sq - 440*(oct2 - 1.)*X_Bsq*X_B*chi2_sq*chi2;
  //ss_3p5PN = ssA_3p5PN + ssB_3p5PN;
  //###########################################################################
  
  T phaseout, spin_spin;
  //T x = pow((M_PI * param->M *f), 2./3.);
  T x = powers->MF2third * powers->PI2third;
  T x_m5_2 = powers->MFminus_5third * powers->PIminus_5third;
  T x_7_2 = powers->MF7third * powers->PI7third;

  //T coeff = (3./(128.* param->eta))* pow(x, -5/2.);
  T coeff = (3./(128.* param->eta))* x_m5_2;

  //spin_spin = ss_2PN * pow(x, 2.) + ss_3PN  * pow(x, 3.) + ss_3p5PN  * pow(x, 7/2.);
  //spin_spin = ss_2PN * x*x+ ss_3PN  * x*x*x+ ss_3p5PN  * x_7_2;
  spin_spin = param->ss_3p5PN_coeff  * x_7_2;
  
  phaseout = coeff*spin_spin;
  //equation 26 of arXiv:1905.06011v2
  return phaseout; 
 
}

template<class T>
T IMRPhenomD_NRT<T>::calculate_NRT_amp_coefficient(source_parameters<T> *param)
{
   return -sqrt(5*M_PI*param->eta / 24.) * (9 * param->M * param->M / param->DL) * (3./16.) * param->tidal_weighted ;//x^(13/4) term -- overall factor

}

template<class T>
T IMRPhenomD_NRT<T>::amp_ins_NRT(T f, useful_powers<T> *powers,source_parameters<T> *param)
{
  /*IMRPhenomD<T> model;
  T gr_ins = model.amp_ins(f, param, pn_coeff, lambda, powers);
  T ampout = gr_ins;
  */

  //T x = pow((M_PI * param->M *f), 2./3.);
  T x = powers->MF2third * powers->PI2third;
  T x2 = x*x;
  T x4 = x2*x2;

  //T amp_NRT = -sqrt(5*M_PI*param->eta / 24.) * (9 * param->M * param->M / param->DL) * (3./16.) * param->tidal_weighted * pow(x, 13./4.) * (1 + (449./108)*x + (22672./9.) * pow(x, 2.89) ) / (1 + 13477.8* pow(x, 4));
  T amp_NRT = param->NRT_amp_coefficient* pow(x, 13./4.) * (1 + (449./108)*x + (22672./9.) * pow(x, 2.89) ) / (1 + 13477.8*x4);
  //has to be scaled by f^(7/6)/A0 to be consistent with the rest of GW analysis tools
  //amp_NRT = amp_NRT*powers->MF7sixth/ (param->A0*pow(param->M, 7/6.));
  //amp_NRT = amp_NRT*pow(f, 7/6.)/(param->A0);
  //amp_NRT = amp_NRT/(param->A0); 

  //ampout += amp_NRT;
  //std::cout<<"Called IMRPhenomD_NRT::amp_ins"<<std::endl; 
  return amp_NRT; 
}



template<class T>
T IMRPhenomD_NRT<T>::taper(T f, int length, source_parameters<T> *params)
{
  //std::cout<<"Called IMRPhenomD_NRT::taper"<<std::endl; 
  IMRPhenomD<T> model;
  T sigma;
  T fmerger, kappa, kappa_eff;
  kappa_eff = (3./16.) * params->tidal_weighted; //kappa effective as defined in arXiv:1804.02235 and arXiv:1905.06011v2.
  
  //kappa = 3.*(params->mass2 * pow(params->mass1, 4.)* params->tidal1/ pow(params->M,5.) + params->mass1 * pow(params->mass2, 4.)* params->tidal2/ pow(params->M,5.));
  //This is the definition for kappa_2^T given in arXiv:1804.02235. As described in that paper,
  //it is accurate enough to use kappa effective. Can test by turning this back on. No change in output.

  double a0 = 0.3586;
  double n1 = 3.35411203e-2;
  double n2 = 4.31460284e-5;
  double d1 = 7.54224145e-2;
  double d2 = 2.23626859e-4;
  /*In order to reproduce LAL's result exactly had to get these coefficients from their code
   instead of the paper because there are more digits given. */
  
  kappa = kappa_eff;
  fmerger = (1./(2.*params->M * M_PI))* a0* sqrt(params->mass2 / params->mass1) *(1.0 + n1 * kappa + n2 * kappa * kappa)/(1.0 + d1 * kappa + d2* kappa * kappa);
  //Equation 11 of arXiv:1804.02235

  T fmerger12 = fmerger*1.2;
  T z = (fmerger - fmerger12)/(f - fmerger) + (fmerger - fmerger12)/(f - fmerger12); 
  
  if(f < fmerger)
    {
      return 0.0; 
    }
  else if(fmerger < f && f < fmerger12)
    {
      sigma = 1.0/(exp(-z) + 1.0);
      return sigma; 
    }
  else if(f > fmerger12)
    {
      return 1.0; 
    }

  return -1.0; //If it gets to this point, something is wrong and the change in sign should indicate that
}


/*Just calls the phenomD construct_waveform and tacks on the phase and amplitude corrections at the end with the taper
 */
template<class T>
int IMRPhenomD_NRT<T>::construct_waveform(T *frequencies, int length, std::complex<T> *waveform, source_parameters<T> *params)
{
  params->NRT_phase_coeff = - (3./16.) * params->tidal_weighted * (39./(16. * params->eta));
  if(params->tidal1<=0){params->oct1=1;params->quad1 = 1;}
  else{
  	params->quad1 = this->calculate_quad_moment(params->tidal1);
  	params->oct1 = this->calculate_oct_moment(params->quad1);
  }
  if(params->tidal2<=0){params->oct2=1;params->quad2 = 1;}
  else{
  	params->quad2 = this->calculate_quad_moment(params->tidal2);
  	params->oct2 = this->calculate_oct_moment(params->quad2);
  }

  this->calculate_spin_coefficients_3p5(params);
  params->NRT_amp_coefficient=this->calculate_NRT_amp_coefficient(params);
  IMRPhenomD<T> model;
  //int status = model.construct_waveform(frequencies, length, waveform, params);
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

  /*This is only here so that fmerger can be printed out and then put on the plots*/
  T fmerger, kappa, kappa_eff;
  kappa_eff = (3./16.) * params->tidal_weighted; //kappa effective as defined in arXiv:1804.02235 and arXiv:1905.06011v2.
  
  //kappa = 3.*(params->mass2 * pow(params->mass1, 4.)* params->tidal1/ pow(params->M,5.) + params->mass1 * pow(params->mass2, 4.)* params->tidal2/ pow(params->M,5.));
  //This is the definition for kappa_2^T given in arXiv:1804.02235. As described in that paper,
  //it is accurate enough to use kappa effective. Can test by turning this back on. No change in output.
  
  kappa = kappa_eff;
  fmerger = (1./(2*params->M * M_PI))* 0.3586* sqrt(params->mass2 / params->mass1) *(1 + 3.354e-2 * kappa + 4.315e-5 * kappa * kappa)/(1 + 7.542e-2 * kappa + 2.236e-4* kappa * kappa);
  

  //Print statements so that we can put the transition frequencies on the plots
  //std::cout<<"f1_phase: "<<params->f1_phase<<std::endl;
  //std::cout<<"f2_phase: "<<params->f2_phase<<std::endl;
  //std::cout<<"f1_amp: "<<params->f1<<std::endl; 
  //std::cout<<"f3_amp: "<<params->f3<<std::endl;
  //std::cout<<"fmerger: "<<fmerger<<std::endl; 
  //std::cout<<"1.2*fmerger: "<<1.2*fmerger<<std::endl; 
	
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
  /*Note -- to match LALSuite, this sets phiRef equal to phi_phenomD at f_ref, not phenomD_NRT
 * Arbitrary constant, not really important
 */
  if(params->shift_phase ){
    f_ref = params->f_ref;
    this->precalc_powers_ins(f_ref, M, &pows);
    phi_shift = (this->build_phase(f_ref,&lambda,params,&pows,pn_phase_coeffs));
    phic = 2*params->phiRef + phi_shift;
  }
  //If phic is specified, ignore f_ref phiRef and use phic
  else{
    f_ref = 0;
    phic = params->phiRef;
  }
  
  //Assign shift: first shift so coalescence happens at t=0, then shift from there according to tc
  /*Note -- to match LALSuite, this sets phiRef equal to phi_phenomD at f_ref, not phenomD_NRT
 * Arbitrary constant, not really important
 */
  if(params->shift_time){
    T alpha1_offset = this->assign_lambda_param_element(params,14);
    tc_shift = this->Dphase_mr(params->f3, params, &lambda)+(-lambda.alpha[1]+alpha1_offset)*params->M/params->eta;
  }
  else{
    tc_shift=0;
  }
  tc = 2*M_PI*params->tc + tc_shift;
  
  T A0 = params->A0* pow(M,7./6.);
  
  T f;
  std::complex<T> amp, phase;
  std::complex<T> i;
  i = std::complex<T> (0,1.);
  T fcut = .2/M; //Cutoff frequency for IMRPhenomD - all higher frequencies return 0
  for (size_t j =0; j< length; j++)
    {
      f = frequencies[j];
      if(f>fcut){
	amp = 0.0;
	waveform[j] = 0.0;
      }
      else{	
	//if (f<params->f1_phase)
	//This is always needed for NRT
	if (true)
	  {
	    this->precalc_powers_ins(f, M, &pows);
	  }
	else
	  {
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
		//I don't know why this would be minus instead of plus, but it seems to get closer to LAL's result if it's minus phaseSpinNRT
		T ampNRT = (A0*this->amp_ins_NRT(f,&pows, params));
		amp +=ampNRT;
               
                // Dissipative tidal contribution to the phase 
                T phaseNRT_D = this->phase_ins_NRT_D(f,&pows,params);
                phase += phaseNRT_D;	
	}
	//phase +=   (T)(tc*(f-f_ref) - phic);
	phase -=   (T)(tc*(f-f_ref) + phic);
	waveform[j] = amp * std::exp(-i * phase);
      }
      
    }
	
  //###################################### The next part applies the taper.  
  for(int i = 0; i<length; i++)
    {
      waveform[i] = waveform[i] * std::complex<T>((T)(1.0) -  taper(frequencies[i], length, params),0);
    }
  
  return 1;
}

/*Old version*/
//template<class T> 
//T IMRPhenomD_NRT<T>::phase_ins_NRT(T f, source_parameters<T> *param, T *pn_coeff, lambda_parameters<T> *lambda, useful_powers<T> *powers)
//{
//  IMRPhenomD<T> model;
//  T gr_ins = model.phase_ins(f, param, pn_coeff, lambda, powers);
//  T phaseout = gr_ins;
//  
//  T phasefix; 
//  bool deriv = false; //tells the Pade function not to take a derivative
//
//  phasefix = - ((3./16.) * param->tidal_weighted * (39./(16. * param->eta)) * Pade(f, param, deriv));
//  phaseout = phaseout + phasefix;
//  
//  return phaseout;
//}


/* This shouldn't be needed anymore, but I'll leave it in. At the moment, it's not called at all (changed the name so its not overloaded)
 */
/*
template<class T>
T IMRPhenomD_NRT<T>::Dphase_ins_NRT(T f, source_parameters<T> *param, T *pn_coeff, lambda_parameters<T> *lambda)
{
  IMRPhenomD<T> model;
  T gr_ins = model.Dphase_ins(f, param, pn_coeff, lambda);
  T phaseout = gr_ins;
  bool deriv = true;
  
  phaseout = phaseout - ((3./16.) * param->tidal_weighted * (39./(16. * param->eta)) * Pade(f, param, deriv));
  return gr_ins; 
}
*/

/*For the amplitude, you may want to call build_amp instead of each amplitude component (ins,int,mr) separately, since LALsuite seems to append the correction
 * to the entire waveform. Up to you though. You would overload that function exactly like the others
 */

/*Still be used as overloaded function*/
 /*template<class T>
T IMRPhenomD_NRT<T>::amp_ins(T f, source_parameters<T> *param, T *pn_coeff,
                        lambda_parameters<T> *lambda, useful_powers<T> *powers)
{
  IMRPhenomD<T> model;
  T gr_ins = model.amp_ins(f, param, pn_coeff, lambda, powers);
  T ampout = gr_ins;

  T x = pow((M_PI * param->M *f), 2./3.);

  T amp_NRT = -sqrt(5*M_PI*param->eta / 24.) * (9 * param->M * param->M / param->DL) * (3./16.) * param->tidal_weighted * pow(x, 13./4.) * (1 + (449./108)*x + (22672./9.) * pow(x, 2.89) ) / (1 + 13477.8* pow(x, 4));
  //has to be scaled by f^(7/6)/A0 to be consistent with the rest of GW analysis tools
  amp_NRT = amp_NRT*powers->MF7sixth/ (param->A0*pow(param->M, 7/6.));

  ampout += amp_NRT;
  //std::cout<<"Called IMRPhenomD_NRT::amp_ins"<<std::endl; 
  return ampout; 
  }*/

/*Still be used as overloaded function*/
/*template<class T>
T IMRPhenomD_NRT<T>::Damp_ins(T f, source_parameters<T> *param, T *pn_coeff, lambda_parameters<T> *lambda)
{
  IMRPhenomD<T> model;
  T gr_ins = model.Damp_ins(f, param, pn_coeff, lambda);
  T ampout = gr_ins; 

  T x = pow((M_PI * param->M *f), 2./3.);
  
  T aNRT = (1 + (449./108.)*x + (22672./9.) * pow(x, 2.89) );
  T aNRTdenom = (1 + 13477.8* pow(x, 4));
  T coefficient =  -sqrt(5*M_PI*param->eta / 24.) * (9 * pow(param->M, 2.0) / param->DL) * (3./16.) * param->tidal_weighted;

  T amp_NRT = coefficient* pow(x,13./4.)* aNRT/aNRTdenom;

  T DaNRT = coefficient * ( (-53911.2 * pow(x, 25./4.) * aNRT /(aNRTdenom * aNRTdenom)) + (pow(x, 13./4.) * ((449./108.) + 7280.23* pow(x, 1.89))/aNRTdenom) + (13*pow(x, 9./4.) * aNRT/(4.*aNRTdenom)) );
  //This is actually the derivative of the amplitude with respect to x
  //to make it the derivative of the amplitude with respect to f, need to multiply by dx/df, done in the next line  
  DaNRT = DaNRT* (2/3.)*pow((param->M*param->M* M_PI*M_PI/f), 1/3.);

  DaNRT = DaNRT *pow(param->M * f, 7/6.)/ (param->A0*pow(param->M, 7/6.)) + amp_NRT * (7/6.)*pow(f,1/6.)/ (param->A0);
  //with appropriate scaling factors
  std::cout<<"DaNRT ="<<DaNRT<<std::endl;
  debugger_print(__FILE__,__LINE__,DaNRT);
  ampout = ampout + DaNRT; 
  return ampout; 
}
*/
 
/*Old version*/

//template<class T>
//int IMRPhenomD_NRT<T>::construct_waveform(T *frequencies, int length, std::complex<T> *waveform, source_parameters<T> *params)
//{
////TESTING 
//  T kappa_eff = (3./16.) * params->tidal_weighted; //kappa effective as defined in arXiv:1804.02235 and arXiv:1905.06011v2.
//  T kappa = kappa_eff;
//  T fmerger = (1./(2*params->M * M_PI))* 0.3586* sqrt(params->mass2 / params->mass1) *(1 + 3.354e-2 * kappa + 4.315e-5 * kappa * kappa)/(1 + 7.542e-2 * kappa + 2.236e-4* kappa * kappa); //Equation 11 of arXiv:1804.02235
////#############################
//  IMRPhenomD<T> model;
//  //int status = model.construct_waveform(frequencies, length, waveform, params);
//  T M = params-> M;
//  T chirpmass = params->chirpmass;
//  T DL = params->DL;
//  lambda_parameters<T> lambda, *lambda_ptr;
//  this->assign_lambda_param(params, &lambda);
//  
//  /*Initialize the post merger quantities*/
//  this->post_merger_variables(params);
//  params->f1_phase = 0.018/(params->M);
//  params->f2_phase = params->fRD/2.;
//  
//  params->f1 = 0.014/(params->M);
//  params->f3 = this->fpeak(params, &lambda);
//
//  /* T fmerger, kappa, kappa_eff;
//  kappa_eff = (3./16.) * params->tidal_weighted; //kappa effective as defined in arXiv:1804.02235 and arXiv:1905.06011v2.
//  
//  //kappa = 3.*(params->mass2 * pow(params->mass1, 4.)* params->tidal1/ pow(params->M,5.) + params->mass1 * pow(params->mass2, 4.)* params->tidal2/ pow(params->M,5.));
//  //This is the definition for kappa_2^T given in arXiv:1804.02235. As described in that paper,
//  //it is accurate enough to use kappa effective. Can test by turning this back on. No change in output.
//  
//  kappa = kappa_eff;
//  fmerger = (1./(2*params->M * M_PI))* 0.3586* sqrt(params->mass2 / params->mass1) *(1 + 3.354e-2 * kappa + 4.315e-5 * kappa * kappa)/(1 + 7.542e-2 * kappa + 2.236e-4* kappa * kappa);
//
//  params->f3 = fmerger;
//  */
//  //Wanted to try using the fmerger defined in the NR Tidal paper instead of fpeak. This made a huge difference...
//	
//  useful_powers<T> pows;
//  this->precalc_powers_PI(&pows);
//
//  T deltas[6];
//  T pn_amp_coeffs[7];
//  T pn_phase_coeffs[12];
//  
//  this->assign_pn_amplitude_coeff(params, pn_amp_coeffs);
//  this->assign_static_pn_phase_coeff(params, pn_phase_coeffs);	
//  
//  this->amp_connection_coeffs(params,&lambda,pn_amp_coeffs,deltas);
//  this->phase_connection_coefficients(params,&lambda,pn_phase_coeffs);
//  
//  //################################################################
//  //Calculate phase and coalescence time variables
//  T phic, f_ref, tc, phi_shift, tc_shift;
//  //If phic is unspecified - use f_ref and phiRef
//  if(params->shift_phase ){
//    f_ref = params->f_ref;
//    this->precalc_powers_ins(f_ref, M, &pows);
//    //NOTE : this is weird, but seems to agree with LIGO
//    //phi_shift = (this->build_phase(f_ref,&lambda,params,&pows,pn_phase_coeffs));
//    phi_shift = (model.build_phase(f_ref,&lambda,params,&pows,pn_phase_coeffs));
//    phic = 2*params->phiRef + phi_shift;
//  }
//  //If phic is specified, ignore f_ref phiRef and use phic
//  else{
//    f_ref = 0;
//    phic = params->phiRef;
//  }
//  
//  //Assign shift: first shift so coalescence happens at t=0, then shift from there according to tc
//  //This aligns more with the physical meaning of tc, but the phase is NO LONGER just
//  if(params->shift_time){
//    T alpha1_offset = this->assign_lambda_param_element(params,14);
//    //NOTE : this is weird, but seems to agree with LIGO
//    //tc_shift = this->Dphase_mr(params->f3, params, &lambda)+(-lambda.alpha[1]+alpha1_offset)*params->M/params->eta;
//    tc_shift = model.Dphase_mr(params->f3, params, &lambda)+(-lambda.alpha[1]+alpha1_offset)*params->M/params->eta;
//  }
//  else{
//    tc_shift=0;
//  }
//  //tc = 2*M_PI*params->tc - tc_shift;
//  tc = 2*M_PI*params->tc + tc_shift;
//  
//  //T A0 = sqrt(M_PI/30)*chirpmass*chirpmass/DL * pow(M_PI*chirpmass,-7./6);
//  T A0 = params->A0* pow(M,7./6.);
//  
//  T f;
//  std::complex<T> amp, phase;
//  std::complex<T> i;
//  i = std::complex<T> (0,1.);
//  T fcut = .2/M; //Cutoff frequency for IMRPhenomD - all higher frequencies return 0
//  for (size_t j =0; j< length; j++)
//    {
//      f = frequencies[j];
//      if(f>fcut){
//	amp = 0.0;
//	waveform[j] = 0.0;
//      }
//      else{	
//	if (f<params->f1_phase)
//	  {
//	    this->precalc_powers_ins(f, M, &pows);
//	  }
//	else
//	  {
//	    pows.MFsixth= pow(M*f,1./6.);	
//	    pows.MF7sixth= pow_int(pows.MFsixth,7);//*pows.MFsixth*pows.MFsixth*pows.MFsixth*pows.MFsixth*pows.MFsixth*pows.MFsixth;
//	  }
//	amp = (A0 * this->build_amp(f,&lambda,params,&pows,pn_amp_coeffs,deltas));
//	phase = (this->build_phase(f,&lambda,params,&pows,pn_phase_coeffs));
//	if(true){
//		T phaseNRT = this->phase_ins_NRT(f, params,pn_phase_coeffs,&lambda,&pows);
//		phase += phaseNRT;
//	}
//	//phase +=   (T)(tc*(f-f_ref) - phic);
//	phase -=   (T)(tc*(f-f_ref) + phic);
//	waveform[j] = amp * std::exp(-i * phase);
//      }
//      
//    }
//	
//  //###################################### The next part applies the taper.  
//  for(int i = 0; i<length; i++)
//    {
//      waveform[i] = waveform[i] * taper(frequencies[i], length, params);
//    }
//  
//  return 1;
//}

template class IMRPhenomD_NRT<double>;
template class IMRPhenomD_NRT<adouble>;
