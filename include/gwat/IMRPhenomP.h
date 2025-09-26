#ifndef IMRPHENOMP_H
#define IMRPHENOMP_H
#include "IMRPhenomD.h"
#include "util.h"

/*! \file 
 * Header file for IMRPhenomP functions
 *
 * Currently, only Pv2 is supported.
 *
 * Wrapped around IMRPhenomD
 */


//Shamelessly stolen from lalsuite
/* Macro functions to rotate the components of a vector about an axis */
#define ROTATEZ(angle, vx, vy, vz)\
tmp1 = vx*cos(angle) - vy*sin(angle);\
tmp2 = vx*sin(angle) + vy*cos(angle);\
vx = tmp1;\
vy = tmp2

#define ROTATEY(angle, vx, vy, vz)\
tmp1 = vx*cos(angle) + vz*sin(angle);\
tmp2 = - vx*sin(angle) + vz*cos(angle);\
vx = tmp1;\
vz = tmp2


typedef enum tagIMRPhenomP_version_type {
 IMRPhenomPv1_V, /**< version 1: based on IMRPhenomC */
 IMRPhenomPv2_V,  /**< version 2: based on IMRPhenomD */
 IMRPhenomPv2NRTidal_V, /**< version Pv2_NRTidal: based on IMRPhenomPv2; NRTides added before precession; can be used with both NRTidal versions defined below */
 IMRPhenomPv3_V  /**< version 3: based on IMRPhenomD and the precession angles from Katerina Chatziioannou PhysRevD.95.104004 (arxiv:1703.03967) */
} IMRPhenomP_version_type;

const double sqrt_6 = 2.44948974278317788;

template<class T>
struct alpha_coeffs
{
	T coeff1 ;
	T coeff2 ;
	T coeff3 ;
	T coeff4 ;
	T coeff5 ;
};

template<class T>
struct epsilon_coeffs
{
	T coeff1 ;
	T coeff2 ;
	T coeff3 ;
	T coeff4 ;
	T coeff5 ;
};
template<class T> 
class IMRPhenomPv2: public IMRPhenomD<T>
{
public:
virtual T alpha(T omega, T q,T chi2l, T chi2);

virtual T epsilon(T omega, T q, T chi2l, T chi2);

virtual void calculate_euler_coeffs(alpha_coeffs<T> *acoeffs, epsilon_coeffs<T> *ecoeffs, source_parameters<T> *params);

virtual T d(int l, int mp, int m,T s);

virtual void PhenomPv2_JSF_from_params(gen_params_base<T> *params, T *JSF);
virtual int construct_waveform(T *frequencies, 
				int length, 
				std::complex<T> *waveform_plus,
				std::complex<T> *waveform_cross,
				source_parameters<T> *params 
				);
virtual int construct_phase(T *frequencies, 
				int length, 
				T *phase_plus,
				T *phase_cross,
				source_parameters<T> *params 
				);

//###############################################################################
virtual T calculate_time_shift(source_parameters<T> *params, useful_powers<T> *pows, T *pn_phase_coeffs, lambda_parameters<T> *lambda);
//###############################################################################

virtual void WignerD(T d2[5], T dm2[5], useful_powers<T> *pows,source_parameters<T> *params);

virtual void calculate_twistup( T alpha, std::complex<T> *hp_factor, std::complex<T> *hc_factor, T d2[5], T dm2[5], sph_harm<T> *sph_harm);

virtual void calculate_euler_angles(T *alpha, T *epsilon, useful_powers<T> *pows, alpha_coeffs<T> *acoeffs, epsilon_coeffs<T> *ecoeffs);

virtual T PhenomPv2_inplane_spin(gen_params_base<T> *params);
virtual void PhenomPv2_Param_Transform(source_parameters<T> *params);

virtual void PhenomPv2_Param_Transform_J(source_parameters<T> *params);

virtual void PhenomPv2_Param_Transform_reduced(source_parameters<T> *params);

virtual T L2PN( T eta, useful_powers<T> *pow);

virtual T FinalSpinIMRPhenomD_all_in_plane_spin_on_larger_BH(T m1, T m2, T chi1_l, T chi2_l, T chip);
virtual T final_spin(source_parameters<T> *params);

};

#endif
