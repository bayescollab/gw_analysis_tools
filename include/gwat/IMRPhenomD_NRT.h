#ifndef IMRPHENOMD_NRT_H
#define IMRPHENOMD_NRT_H
#include "IMRPhenomD.h"
#include "util.h"

/*! \file
 */

/*! Class that extends the IMRPhenomD waveform to include tidal effects in the 
 *inspiral portion of the phase.
 */

template<class T>
class IMRPhenomD_NRT: public IMRPhenomD<T>
{
public:
  virtual T Pade(T f, source_parameters<T> *param, useful_powers<T> *powers,char deriv);
  
  //virtual T phase_ins_NRT(T f, source_parameters<T> *param, T *pn_coeff, lambda_parameters<T> *lambda, useful_powers<T> *pow); 
  //This is NOT overloaded anymore -- completely different function than phase_ins
  virtual T phase_ins_NRT(T f,useful_powers<T> *powers, source_parameters<T> *param);

  //Adds the dissipative tidal correction. There is no IMRPhenom correction--just the PN correction.
  // See 2306.15633 
  virtual T phase_ins_NRT_D(T f, useful_powers<T> *powers, source_parameters<T> *param);

  //virtual T spin_spin(source_parameters<T> *param, double PNorder, int body); 
  virtual T phase_spin_NRT(T f, useful_powers<T> *powers,source_parameters<T> *param);
    
  //Probably not needed anymore, but I'll leave for now
  //virtual T Dphase_ins_NRT(T r, source_parameters<T> *param, T *pn_coeff, lambda_parameters<T> *lambda); 

  //virtual T amp_ins(T f, source_parameters<T> *param, T *pn_coeff,
  //lambda_parameters<T> *lambda,useful_powers<T> *pow);
  virtual T amp_ins_NRT(T f, useful_powers<T> *powers,source_parameters<T> *param);
//virtual T Damp_ins(T f, source_parameters<T> *param, T *pn_coeff, lambda_parameters<T> *lambda);
 
  virtual int construct_waveform(T *frequencies, int length, std::complex<T> *waveform, source_parameters<T> *params);
  
  virtual T taper(T f, int length, source_parameters<T> *params); 


  virtual void assign_static_pn_phase_coeff(source_parameters<T> *source_param, T *coeff);

  virtual T calculate_quad_moment(T lambda);
  virtual T calculate_oct_moment(T lambda);

  virtual void calculate_spin_coefficients_3p5(source_parameters<T> *param);

  virtual T calculate_NRT_amp_coefficient(source_parameters<T> *param);

  virtual void binary_love_relation(T tidal_s, bool tidal_love_error, source_parameters<T> *sp);

  //virtual T tidal_error(T tidal_s, T tidal_a, T q);

};

double tidal_error(double tidal_s, double tidal_a, double q);
adouble tidal_error(adouble tidal_s, adouble tidal_a, adouble q);

//template<class T>
//void binary_love_relation(T tidal_s, bool tidal_love_error, source_parameters<T> *sp); 

//##############################################################################

/*!Numerically calibrated coefficients of Pade approximant (PNRTidal_v2) 
 *from arXiv:1905.06011 equations 19, 20, and 21.
 */
/*const double c_NRT[4] = {3115./1248., -M_PI, 28024205./3302208., -4283.* M_PI/1092.};
const double d_NRT[3]={-15.111208, -(c_NRT[3] + (c_NRT[1] * d_NRT[0]) - 90.550822)/c_NRT[0], 8.0641096};
const double n_NRT[5]={c_NRT[0] + d_NRT[0], ((c_NRT[0] * c_NRT[1]) - c_NRT[3] - (c_NRT[1] * d_NRT[0]) + 90.550822)/c_NRT[0], c_NRT[2] + (c_NRT[0] * d_NRT[0]) + d_NRT[2], 90.550822, -60.253578};
*/
//The ones in the paper aren't as precise, so we instead copied the values used in LAL's code. (line 256-264 of LALSimNRTunedTides.c). It does make a difference
const double n_NRT[5]= {-12.615214237993088, 19.0537346970349, -21.166863146081035, 90.55082156324926, -60.25357801943598};
const double d_NRT[3]= {-15.111207827736678, 22.195327350624694, 8.064109635305156};

/*! Numerically calibrated coefficients of binary love relations from 
 * arXiv:1903.03909 equations 11, 12, and 13.
 */
const double n_binLove = 0.743; 
const double b_binLove[3][2] = {{-14.40, 14.45},{31.36,-32.25},{-22.44,20.35}};
const double c_binLove[3][2] = {{-15.25,15.37}, {37.33,-43.20}, {-29.93,35.18}};

/*! Numerically calibrated coefficients for error marginalization over residual
 * EoS dependence of binary love relations from arXiv:1903.03909 table V. These 
 * coefficients are used in equations 19-22 of that paper.
 */
const double mu_binLove[5] = {3.509*pow(10,-3.), 9.351*pow(10,-1.), -18.07, 27.56, -10.10};
const double sigma_binLove[9] = {-2.074*pow(10,-7.), 1.492*pow(10,-3.), -4.891*pow(10,-2.), 8.207*pow(10,-1.), -1.308, -63.76, 11.14, 75.25, -23.69};
const double bL_error[18] = {-12.9277, 24.3532, -11.4391, 1.14402*pow(10,-2.), -7.91351*pow(10,-3.), -3.40364*pow(10,-3.), -1.50512*pow(10, -7.), 15.4219, -46.1750, 56.3086, -25.0453, 3.49200*pow(10, -8.), -1.85544*pow(10, -4.), 6.17533*pow(10, -2.), -4.65122, 11.1415, -6.59985, -5.40616*pow(10, -2.)};

#endif
