#include "IMRPhenomP.h"
#include <iostream>
#include <fstream>
#include <string>
#include <complex>
#include "IMRPhenomD.h"
#include "util.h"
#include <adolc/adouble.h>
#include <math.h>
#include <algorithm>
#include <type_traits>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
/*! \file 
 *
 * Source code for IMRPhenomP
 *
 * Verified against LALsuite 2019_09_25
 */



template<class T>
T IMRPhenomPv2<T>::alpha(T omega, T q,T chi2l, T chi2){

	T alpha;
	alpha = (-5*(338688*pow_int(1 + q,4)*(3 + 4*q) + 508032*chi2l*pow(omega,0.3333333333333333)*q*pow_int(1 + q,4)*
        (3 + 4*q) + 3024*pow(omega,0.6666666666666666)*pow_int(1 + q,2)*
        (2985 + q*(12890 + q*(15789 + 4988*q + 168*pow_int(chi2,2)*pow_int(1 + q,2)*(3 + 4*q)))) + 
       pow(omega,1.3333333333333333)*(-17660607 + 
          q*(-107348840 + 4064256*chi2l*M_PI*pow_int(1 + q,4)*(3 + 4*q) - 
             84672*pow_int(chi2l,2)*q*pow_int(1 + q,2)*(3 + 4*q)*
              (75 + q*(113 + 6*pow_int(chi2,2)*q*pow_int(1 + q,2))) + 
             q*(-271003598 + 127008*pow_int(chi2,4)*pow_int(q,2)*pow_int(1 + q,4)*(3 + 4*q) - 
                q*(327403764 + q*(181442579 + 39432548*q)) - 
                1512*pow_int(chi2,2)*pow_int(1 + q,2)*(213 + q*(1802 + q*(2909 + 956*q)))))) + 
       1008*omega*pow_int(1 + q,2)*(1344*M_PI*pow_int(1 + q,2)*(3 + 4*q) + 
          chi2l*q*(-5253 + q*(-18854 + q*(-18197 - 2972*q + 168*pow_int(chi2,2)*pow_int(1 + q,2)*(3 + 4*q)))))*
        log(omega)))/(6.5028096e7*omega*q*pow_int(1 + q,4));
	return alpha;
}


template<class T>
T IMRPhenomPv2<T>::epsilon(T omega, T q, T chi2l, T chi2)
{
	T epsilon;
	epsilon = (-5*(338688*pow_int(1 + q,4)*(3 + 4*q) + 508032*chi2l*pow(omega,0.3333333333333333)*q*pow_int(1 + q,4)*
        (3 + 4*q) + 3024*pow(omega,0.6666666666666666)*pow_int(1 + q,2)*
        (2985 + q*(12890 + q*(15789 + 4988*q))) + 
       pow(omega,1.3333333333333333)*(-17660607 + 
          q*(-107348840 + 4064256*chi2l*M_PI*pow_int(1 + q,4)*(3 + 4*q) - 
             84672*pow_int(chi2l,2)*q*pow_int(1 + q,2)*(3 + 4*q)*(75 + 113*q) - 
             q*(271003598 + q*(327403764 + q*(181442579 + 39432548*q))))) - 
       1008*omega*pow_int(1 + q,2)*(-1344*M_PI*pow_int(1 + q,2)*(3 + 4*q) + 
          chi2l*q*(5253 + q*(18854 + q*(18197 + 2972*q))))*log(omega)))/(6.5028096e7*omega*q*pow_int(1 + q,4));
	return epsilon;
}

/*! \brief Pre calculate euler angle coefficients
 *
 * Straight up stolen from LALsuite
 */
template<class T>
void IMRPhenomPv2<T>::calculate_euler_coeffs(alpha_coeffs<T> *acoeffs, epsilon_coeffs<T> *ecoeffs, source_parameters<T> *params)
{
	T m2 = params->q/(1. + params->q);
	T m1 = 1./(1. +params-> q);
	T dm = m1 - m2;
	T mtot = 1.;
	T eta = m1*m2; /* mtot = 1 */
	T eta2 = eta*eta;
	T eta3 = eta2*eta;
	T eta4 = eta3*eta;
	T mtot2 = mtot*mtot;
	T mtot4 = mtot2*mtot2;
	T mtot6 = mtot4*mtot2;
	T mtot8 = mtot6*mtot2;
	T chil2 = params->chil*params->chil;
	T chip2 = params->chip*params->chip;
	T chip4 = chip2*chip2;
	T dm2 = dm*dm;
	T dm3 = dm2*dm;
	T m2_2 = m2*m2;
	T m2_3 = m2_2*m2;
	T m2_4 = m2_3*m2;
	T m2_5 = m2_4*m2;
	T m2_6 = m2_5*m2;
	T m2_7 = m2_6*m2;
	T m2_8 = m2_7*m2;
	T chil = params->chil;
	
	acoeffs->coeff1 = (-0.18229166666666666 - (5*dm)/(64.*m2));

	acoeffs->coeff2 = ((-15*dm*m2*chil)/(128.*mtot2*eta) - (35*m2_2*chil)/(128.*mtot2*eta));
	
	acoeffs->coeff3 = (-1.7952473958333333 - (4555*dm)/(7168.*m2) -
	      (15*chip2*dm*m2_3)/(128.*mtot4*eta2) -
	      (35*chip2*m2_4)/(128.*mtot4*eta2) - (515*eta)/384. - (15*dm2*eta)/(256.*m2_2) -
	      (175*dm*eta)/(256.*m2));
	
	acoeffs->coeff4 = - (35*M_PI)/48. - (5*dm*M_PI)/(16.*m2) +
	   (5*dm2*chil)/(16.*mtot2) + (5*dm*m2*chil)/(3.*mtot2) +
	   (2545*m2_2*chil)/(1152.*mtot2) -
	   (5*chip2*dm*m2_5*chil)/(128.*mtot6*eta3) -
	   (35*chip2*m2_6*chil)/(384.*mtot6*eta3) + (2035*dm*m2*chil)/(21504.*mtot2*eta) +
	   (2995*m2_2*chil)/(9216.*mtot2*eta);
	
	acoeffs->coeff5 = (4.318908476114694 + (27895885*dm)/(2.1676032e7*m2) -
	      (15*chip4*dm*m2_7)/(512.*mtot8*eta4) -
	      (35*chip4*m2_8)/(512.*mtot8*eta4) -
	      (485*chip2*dm*m2_3)/(14336.*mtot4*eta2) +
	      (475*chip2*m2_4)/(6144.*mtot4*eta2) +
	      (15*chip2*dm2*m2_2)/(256.*mtot4*eta) + (145*chip2*dm*m2_3)/(512.*mtot4*eta) +
	      (575*chip2*m2_4)/(1536.*mtot4*eta) + (39695*eta)/86016. + (1615*dm2*eta)/(28672.*m2_2) -
	      (265*dm*eta)/(14336.*m2) + (955*eta2)/576. + (15*dm3*eta2)/(1024.*m2_3) +
	      (35*dm2*eta2)/(256.*m2_2) + (2725*dm*eta2)/(3072.*m2) - (15*dm*m2*M_PI*chil)/(16.*mtot2*eta) -
	      (35*m2_2*M_PI*chil)/(16.*mtot2*eta) + (15*chip2*dm*m2_7*chil2)/(128.*mtot8*eta4) +
	      (35*chip2*m2_8*chil2)/(128.*mtot8*eta4) +
	      (375*dm2*m2_2*chil2)/(256.*mtot4*eta) + (1815*dm*m2_3*chil2)/(256.*mtot4*eta) +
	      (1645*m2_4*chil2)/(192.*mtot4*eta));

	ecoeffs->coeff1 = (-0.18229166666666666 - (5*dm)/(64.*m2));
	
	ecoeffs->coeff2 = ((-15*dm*m2*chil)/(128.*mtot2*eta) - (35*m2_2*chil)/(128.*mtot2*eta));
	
	ecoeffs->coeff3 = (-1.7952473958333333 - (4555*dm)/(7168.*m2) - (515*eta)/384. -
	      (15*dm2*eta)/(256.*m2_2) - (175*dm*eta)/(256.*m2));
	
	ecoeffs->coeff4 = - (35*M_PI)/48. - (5*dm*M_PI)/(16.*m2) +
	   (5*dm2*chil)/(16.*mtot2) + (5*dm*m2*chil)/(3.*mtot2) +
	   (2545*m2_2*chil)/(1152.*mtot2) + (2035*dm*m2*chil)/(21504.*mtot2*eta) +
	   (2995*m2_2*chil)/(9216.*mtot2*eta);
	
	ecoeffs->coeff5 = (4.318908476114694 + (27895885*dm)/(2.1676032e7*m2) + (39695*eta)/86016. +
	      (1615*dm2*eta)/(28672.*m2_2) - (265*dm*eta)/(14336.*m2) + (955*eta2)/576. +
	      (15*dm3*eta2)/(1024.*m2_3) + (35*dm2*eta2)/(256.*m2_2) +
	      (2725*dm*eta2)/(3072.*m2) - (15*dm*m2*M_PI*chil)/(16.*mtot2*eta) - (35*m2_2*M_PI*chil)/(16.*mtot2*eta) +
	      (375*dm2*m2_2*chil2)/(256.*mtot4*eta) + (1815*dm*m2_3*chil2)/(256.*mtot4*eta) +
	      (1645*m2_4*chil2)/(192.*mtot4*eta));



}
template<class T>
T IMRPhenomPv2<T>::d(int l, int m_prime, int m,T s)
{
	T sqrt2 = sqrt(2);
	T s_sqrt = sqrt(1+s*s);
	/*cos(beta/2)*/
	T cb = (1./sqrt2)*sqrt(1 + 1./s_sqrt);
	/*sin(beta/2)*/
	T sb = (1./sqrt2) * sqrt( 1. - 1./s_sqrt);
	
	T overall_factor = sqrt(factorial(l+m) * factorial(l-m) * factorial(l +m_prime) *factorial(l-m_prime));
	
	/*Limits of the sum (determined by keeping all factorials positive)*/
	int kmax;	
	int kmax_options[3];
	if (l == std::abs(m))
		kmax = 0;
	else
	{
		if(m < m_prime) return 0;
		else
		{
			kmax_options[0] = m-m_prime;	
			kmax_options[1] = l+m;
			kmax_options[2] = l-m;
			kmax = *std::min_element(kmax_options, kmax_options+3);
		}
	}

	/*Compute the rotation matrix*/
	int k = 0;
	T sum = 0;
	while(k<=kmax)
	{
		sum=sum+ pow_int(-1.,k+m_prime-m)/ ( factorial(l+m-k)*factorial(l-m-k)*
					factorial(k)*factorial(m_prime-m+k) ) *
					pow_int(cb,2*l-2*k-m_prime+m)*pow_int(sb,2*k+m_prime-m);
		k++;
	}
	return sum;
	
	
}

/*! \brief Constructs the waveform for IMRPhenomPv2 - uses IMRPhenomD, then twists up
 *
 * arguments:
 * 	array of frequencies, length of that array, a complex array for the output waveform, and a source_parameters structure
 */
template <class T>
int IMRPhenomPv2<T>::construct_waveform(T *frequencies, /**< T array of frequencies the waveform is to be evaluated at*/
				int length, /**< integer length of the array of frequencies and the waveform*/	
				std::complex<T> *waveform_plus,/**< complex T array for the plus polariaztion waveform to be output*/ 
				std::complex<T> *waveform_cross,/**< complex T array for the cross polarization waveform to be output*/ 
				source_parameters<T> *params /*Structure of source parameters to be initialized before computation*/
				)
{

	//Initialize Spherical harmonics for polarization construction
	sph_harm<T> harmonics;
	T phiHarm = 0.;
	harmonics.Y22 = std::complex<T>(0.0,0.0);
	harmonics.Y21 = std::complex<T>(0.0,0.0);
	harmonics.Y20 = std::complex<T>(0.0,0.0);
	harmonics.Y2m1 =std::complex<T>(0.0,0.0);
	harmonics.Y2m2 =std::complex<T>(0.0,0.0);
	harmonics.Y22 = XLALSpinWeightedSphericalHarmonic(params->thetaJN,phiHarm, -2,2,2);
	harmonics.Y21 = XLALSpinWeightedSphericalHarmonic(params->thetaJN,phiHarm, -2,2,1);
	harmonics.Y20 = XLALSpinWeightedSphericalHarmonic(params->thetaJN,phiHarm, -2,2,0);
	harmonics.Y2m1 = XLALSpinWeightedSphericalHarmonic(params->thetaJN,phiHarm, -2,2,-1);
	harmonics.Y2m2 = XLALSpinWeightedSphericalHarmonic(params->thetaJN,phiHarm, -2,2,-2);
	
	T M = params-> M;
	T chirpmass = params->chirpmass;
	T DL = params->DL;
	T eta = params->eta;
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
	//T pn_phase_coeffs[10];
	T pn_phase_coeffs[12];
	

	this->assign_pn_amplitude_coeff(params, pn_amp_coeffs);
	this->assign_static_pn_phase_coeff(params, pn_phase_coeffs);	

	this->amp_connection_coeffs(params,&lambda,pn_amp_coeffs,deltas);
	this->phase_connection_coefficients(params,&lambda,pn_phase_coeffs);
	//#################################################################
	//T phic, f_ref, tc, phi_shift, tc_shift;
	//params->phiRef = params->phi_aligned;
	//f_ref = params->f_ref;
	//phic = 2*params->phiRef;
	T phic, f_ref, tc, phi_shift, tc_shift;
	f_ref = params->f_ref;
	phic = 2*params->phi_aligned;
	tc = 2*M_PI*params->tc;
	//#################################################################

	//Rescale amplitude because we aren't just using (2,2) mode anymore
	T A0 = params->A0* pow(M,7./6.) / (2. * sqrt(5. / (64.*M_PI)) );

	T q = params->mass1/params->mass2;
	params->q = q;
	T d2[5] ;
	T dm2[5];
	std::complex<T> hp_factor = std::complex<T>(0.0,0.0);
	std::complex<T> hc_factor = std::complex<T>(0.0,0.0);
	T epsilon, epsilon_offset;
	T alpha, alpha_offset;

	alpha_coeffs<T> acoeffs;
	epsilon_coeffs<T> ecoeffs;
	this->calculate_euler_coeffs(&acoeffs, &ecoeffs, params);

	T f;
	std::complex<T> amp, phase;
	std::complex<T> *amp_vec = new std::complex<T>[length];
	std::complex<T> *phase_vec =new std::complex<T>[length];
	std::complex<T> *hpfac_vec =new std::complex<T>[length];
	std::complex<T> *hcfac_vec= new std::complex<T>[length];
	
	//T amp, phase;
	std::complex<T> i;
	i = std::complex<T> (0,1.);


	//Calculate offsets at fRef
	pows.MFthird = pow(params->M* params->f_ref, 1./3.);
	pows.MF2third =pows.MFthird* pows.MFthird;
	this->calculate_euler_angles(&alpha_offset, &epsilon_offset, &pows, &acoeffs, &ecoeffs);
	T fcut = .2/M; //Cutoff frequency for IMRPhenomD - all higher frequencies return 0
	for (int j =0; j< length; j++)
	{
		f = frequencies[j];
		if(f>fcut){
			amp = 0.0;
			waveform_plus[j] = 0.0;
			waveform_cross[j] = 0.0;
			amp_vec[j] = amp;
			//################################################
			phase_vec[j] = 0;
			hpfac_vec[j] = 0;
			hcfac_vec[j] = 0;
		}
		else{	
			//if (f<params->f1_phase)
			//{
			//	this->precalc_powers_ins(f, M, &pows);
			//}
			//else{
			//	pows.MFsixth = pow(M*f,1./6 );
			//	pows.MF7sixth= pow_int(pows.MFsixth,7);
			//	pows.MFthird = pows.MFsixth * pows.MFsixth;
			//	pows.MF2third =pows.MFthird* pows.MFthird;
			//}
			this->precalc_powers_ins(f, M, &pows);
			amp = (A0 * this->build_amp(f,&lambda,params,&pows,pn_amp_coeffs,deltas));
			phase = (this->build_phase(f,&lambda,params,&pows,pn_phase_coeffs));
			//Calculate WignerD matrices -- See mathematica nb for the forms: stolen from lalsuite
			this->WignerD(d2,dm2, &pows, params);
			//Calculate Euler angles alpha and epsilon
			this->calculate_euler_angles(&alpha, &epsilon, &pows, &acoeffs, &ecoeffs);
			//Add offset to alpha
			alpha = alpha + params->alpha0 - alpha_offset;
			epsilon = epsilon - epsilon_offset;

			//Twist it up
			calculate_twistup(alpha, &hp_factor, &hc_factor, d2, dm2, &harmonics);

			phase = phase + (std::complex<T>)(2. * epsilon) ;
			//################################################
			amp_vec[j] = amp/std::complex<T>(2.,0.0);
			//################################################
			phase_vec[j] = phase;
			hpfac_vec[j] = hp_factor;
			hcfac_vec[j] = hc_factor;
			
			hp_factor = 0.;
			hc_factor = 0.;
		}

	}
	//#########################################################
	//Because this part requires the use of gsl_interpolation, and that requires
	//a specific type (adouble not supported), this must be done for each specific
	//type. The implementations are fundamentally different, so these templates must be 
	//written explicitly. No way around it.
	T t_corr_fixed;
	if(params->shift_time){
		t_corr_fixed = this->calculate_time_shift(params, &pows, pn_phase_coeffs, &lambda);
	}
	else{
		t_corr_fixed = 0;
	}
	//#########################################################
	for (int j = 0; j<length;j++){
		waveform_plus[j] = 
			amp_vec[j] * hpfac_vec[j]*std::exp(
			-i * (phase_vec[j]-std::complex<T>(tc*(frequencies[j]-f_ref) ,0)- phic+(std::complex<T>)(2*M_PI*t_corr_fixed*frequencies[j])));
		waveform_cross[j] = 
			amp_vec[j] * hcfac_vec[j]*std::exp(
			-i *(phase_vec[j]-std::complex<T>(tc*(frequencies[j]-f_ref) ,0)- phic+(std::complex<T>)(2*M_PI*t_corr_fixed)*frequencies[j]));
		
	}
	delete [] amp_vec;
	delete [] phase_vec;
	delete [] hpfac_vec;
	delete [] hcfac_vec;
	return 1;
	
}
/*! \brief Constructs the phase for IMRPhenomPv2 - uses IMRPhenomD, then twists up
 *
 * arguments:
 * 	array of frequencies, length of that array, a complex array for the output waveform, and a source_parameters structure
 */
template <class T>
int IMRPhenomPv2<T>::construct_phase(T *frequencies, /**< T array of frequencies the waveform is to be evaluated at*/
				int length, /**< integer length of the array of frequencies and the waveform*/	
				T *phase_plus,/**< complex T array for the plus polariaztion waveform to be output*/ 
				T *phase_cross,/**< complex T array for the cross polarization waveform to be output*/ 
				source_parameters<T> *params /*Structure of source parameters to be initialized before computation*/
				)
{


	//Initialize Spherical harmonics for polarization construction
	sph_harm<T> harmonics;
	T phiHarm = 0.;
	harmonics.Y22 = std::complex<T>(0.0,0.0);
	harmonics.Y21 = std::complex<T>(0.0,0.0);
	harmonics.Y20 = std::complex<T>(0.0,0.0);
	harmonics.Y2m1 =std::complex<T>(0.0,0.0);
	harmonics.Y2m2 =std::complex<T>(0.0,0.0);
	harmonics.Y22 = XLALSpinWeightedSphericalHarmonic(params->thetaJN,phiHarm, -2,2,2);
	harmonics.Y21 = XLALSpinWeightedSphericalHarmonic(params->thetaJN,phiHarm, -2,2,1);
	harmonics.Y20 = XLALSpinWeightedSphericalHarmonic(params->thetaJN,phiHarm, -2,2,0);
	harmonics.Y2m1 = XLALSpinWeightedSphericalHarmonic(params->thetaJN,phiHarm, -2,2,-1);
	harmonics.Y2m2 = XLALSpinWeightedSphericalHarmonic(params->thetaJN,phiHarm, -2,2,-2);
	
	T M = params-> M;
	T chirpmass = params->chirpmass;
	T DL = params->DL;
	T eta = params->eta;
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
	//T pn_phase_coeffs[10];
	T pn_phase_coeffs[12];
	

	this->assign_pn_amplitude_coeff(params, pn_amp_coeffs);
	this->assign_static_pn_phase_coeff(params, pn_phase_coeffs);	

	this->amp_connection_coeffs(params,&lambda,pn_amp_coeffs,deltas);
	this->phase_connection_coefficients(params,&lambda,pn_phase_coeffs);
	//#################################################################
	//T phic, f_ref, tc, phi_shift, tc_shift;
	//params->phiRef = params->phi_aligned;
	//f_ref = params->f_ref;
	//phic = 2*params->phiRef;
	T phic, f_ref, tc, phi_shift, tc_shift;
	f_ref = params->f_ref;
	phic = 2*params->phi_aligned;
	tc = 2*M_PI*params->tc;
	//#################################################################

	//Rescale amplitude because we aren't just using (2,2) mode anymore
	T A0 = params->A0* pow(M,7./6.) / (2. * sqrt(5. / (64.*M_PI)) );

	T q = params->mass1/params->mass2;
	params->q = q;
	T d2[5] ;
	T dm2[5];
	std::complex<T> hp_factor = std::complex<T>(0.0,0.0);
	std::complex<T> hc_factor = std::complex<T>(0.0,0.0);
	T epsilon, epsilon_offset;
	T alpha, alpha_offset;

	alpha_coeffs<T> acoeffs;
	epsilon_coeffs<T> ecoeffs;
	this->calculate_euler_coeffs(&acoeffs, &ecoeffs, params);

	T f;
	T amp, phase;
	T *phase_vec =new T[length];
	std::complex<T> *hpfac_vec =new std::complex<T>[length];
	std::complex<T> *hcfac_vec= new std::complex<T>[length];
	
	//T amp, phase;
	std::complex<T> i;
	i = std::complex<T> (0,1.);


	//Calculate offsets at fRef
	pows.MFthird = pow(params->M* params->f_ref, 1./3.);
	pows.MF2third =pows.MFthird* pows.MFthird;
	this->calculate_euler_angles(&alpha_offset, &epsilon_offset, &pows, &acoeffs, &ecoeffs);
	T fcut = .2/M; //Cutoff frequency for IMRPhenomD - all higher frequencies return 0
	for (int j =0; j< length; j++)
	{
		f = frequencies[j];
		if(f>fcut){
			//################################################
			phase_vec[j] = 0;
			hpfac_vec[j] = 0;
			hcfac_vec[j] = 0;
		}
		else{	
			//if (f<params->f1_phase)
			//{
			//	this->precalc_powers_ins(f, M, &pows);
			//}
			//else{
			//	pows.MFsixth = pow(M*f,1./6 );
			//	pows.MF7sixth= pow_int(pows.MFsixth,7);
			//	pows.MFthird = pows.MFsixth * pows.MFsixth;
			//	pows.MF2third =pows.MFthird* pows.MFthird;
			//}
			this->precalc_powers_ins(f, M, &pows);
			//amp = (A0 * this->build_amp(f,&lambda,params,&pows,pn_amp_coeffs,deltas));
			phase = (this->build_phase(f,&lambda,params,&pows,pn_phase_coeffs));
			//Calculate WignerD matrices -- See mathematica nb for the forms: stolen from lalsuite
			this->WignerD(d2,dm2, &pows, params);
			//Calculate Euler angles alpha and epsilon
			this->calculate_euler_angles(&alpha, &epsilon, &pows, &acoeffs, &ecoeffs);
			//Add offset to alpha
			alpha = alpha + params->alpha0 - alpha_offset;
			epsilon = epsilon - epsilon_offset;

			//Twist it up
			calculate_twistup(alpha, &hp_factor, &hc_factor, d2, dm2, &harmonics);

			phase = phase + (2. * epsilon) ;
			//################################################
			//amp_vec[j] = amp/std::complex<T>(2.,0.0);
			//################################################
			phase_vec[j] = phase;
			hpfac_vec[j] = hp_factor;
			hcfac_vec[j] = hc_factor;
			
			hp_factor = 0.;
			hc_factor = 0.;
		}

	}
	//#########################################################
	//Because this part requires the use of gsl_interpolation, and that requires
	//a specific type (adouble not supported), this must be done for each specific
	//type. The implementations are fundamentally different, so these templates must be 
	//written explicitly. No way around it.
	T t_corr_fixed;
	if(params->shift_time){
		t_corr_fixed = this->calculate_time_shift(params, &pows, pn_phase_coeffs, &lambda);
	}
	else{
		t_corr_fixed = 0;
	}
	//#########################################################
	for (int j = 0; j<length;j++){
		phase_plus[j] = 
			arg(hpfac_vec[j])-
			(phase_vec[j]-(tc*(frequencies[j]-f_ref) )- phic+(2*M_PI*t_corr_fixed*frequencies[j]));
		phase_cross[j] = 
			arg(hcfac_vec[j])-
			(phase_vec[j]-(tc*(frequencies[j]-f_ref) )- phic+(2*M_PI*t_corr_fixed*frequencies[j]));
		
	}
	delete [] phase_vec;
	delete [] hpfac_vec;
	delete [] hcfac_vec;
	return 1;
	
}
/*! \brief Shifts the time of coalescence to the desired value
 *
 * Because GSL interpolation must have double (not adouble), the two cases must behandled  separately, explicitly.
 */
template<>
double IMRPhenomPv2<double>::calculate_time_shift(source_parameters<double> *params, useful_powers<double> *pows, double *pn_phase_coeffs, lambda_parameters<double> *lambda)
{
	double M = params->M;
	//#########################################################
	//Interpolate phase to set coalescence time to 0
	gsl_interp_accel *acc_fixed = NULL;
	gsl_spline *phiI_fixed = NULL;

	int n_fixed = 10;
	
	double f_final = params->fRD;
	double freqs_fixed[n_fixed] ; 	
	double phase_fixed[n_fixed] ; 	
	double freqs_fixed_start = .8*(f_final);
	double freqs_fixed_stop = 1.2*(f_final);
	double delta_freqs_fixed  =  (freqs_fixed_stop-freqs_fixed_start)/(n_fixed-1);
	if(std::isnan(delta_freqs_fixed) || delta_freqs_fixed<= 0)
	{
		debugger_print(__FILE__,__LINE__,"Issue with phase interpolation in PhenomPv");
		debugger_print(__FILE__,__LINE__,"Delta_f, f_final: "+std::to_string(delta_freqs_fixed)+" , "+std::to_string(f_final));
		return 0;
	}
	for(int j =0; j<n_fixed; j++){
		freqs_fixed[j] = freqs_fixed_start +  j *delta_freqs_fixed;
	}
	for (int j =  0 ;j<n_fixed; j++){
		//std::complex<T> hphenom = 0.0;
		double phasing = 0  ;
		double f =  freqs_fixed[j];
		if (f<params->f1_phase)
		{
			this->precalc_powers_ins(f, M, pows);
		}
		else{
			pows->MFsixth = pow(M*f,1./6 );
			pows->MF7sixth= pows->MFsixth*pows->MFsixth*pows->MFsixth*pows->MFsixth*pows->MFsixth*pows->MFsixth*pows->MFsixth;
			pows->MFthird = pows->MFsixth * pows->MFsixth;
			pows->MF2third =pows->MFthird* pows->MFthird;
		}
		phasing = (this->build_phase(f,lambda,params,pows,pn_phase_coeffs));
		phase_fixed[j]= -phasing;
	}

	//Interpolation to  fix tc
	acc_fixed = gsl_interp_accel_alloc();
	phiI_fixed = gsl_spline_alloc(gsl_interp_cspline,n_fixed);
	gsl_spline_init(phiI_fixed, freqs_fixed, phase_fixed, n_fixed);
	//double t_corr_fixed =  gsl_spline_eval_deriv(phiI_fixed, f_final, acc_fixed)/(2.*M_PI) + params->tc;
	double t_corr_fixed =  gsl_spline_eval_deriv(phiI_fixed, f_final, acc_fixed)/(2.*M_PI) ;
	//Clean up
	gsl_spline_free(phiI_fixed);
	gsl_interp_accel_free(acc_fixed);
	//#########################################################
	return t_corr_fixed;

}
/*! \brief Shifts the time of coalescence to the desired value
 *
 * Because GSL interpolation must have double (not adouble), the two cases must behandled  separately, explicitly.
 */
template<>
adouble IMRPhenomPv2<adouble>::calculate_time_shift(source_parameters<adouble> *params, useful_powers<adouble> *pows, adouble *pn_phase_coeffs, lambda_parameters<adouble> *lambda)
{
	adouble M = params->M;
	//#########################################################
	//Interpolate phase to set coalescence time to 0
	gsl_interp_accel *acc_fixed = NULL;
	gsl_spline *phiI_fixed = NULL;

	int n_fixed = 10;
	
	adouble f_final_tmp =(adouble) params->fRD;//params->M;
	double f_final = f_final_tmp.value();
	double freqs_fixed[n_fixed] ; 	
	double phase_fixed[n_fixed] ; 	
	double freqs_fixed_start = .8*(f_final);
	double freqs_fixed_stop = 1.2*(f_final);
	double delta_freqs_fixed  =  (freqs_fixed_stop-freqs_fixed_start)/(n_fixed-1);
	for(int j =0; j<n_fixed; j++){
		freqs_fixed[j] = freqs_fixed_start +  j *delta_freqs_fixed;
	}
	for (int j =  0 ;j<n_fixed; j++){
		//std::complex<T> hphenom = 0.0;
		adouble phasing = 0  ;
		adouble f =  freqs_fixed[j];
		if (f<params->f1_phase)
		{
			this->precalc_powers_ins(f, M, pows);
		}
		else{
			pows->MFsixth = pow(M*f,1./6 );
			pows->MF7sixth= pows->MFsixth*pows->MFsixth*pows->MFsixth*pows->MFsixth*pows->MFsixth*pows->MFsixth*pows->MFsixth;
			pows->MFthird = pows->MFsixth * pows->MFsixth;
			pows->MF2third =pows->MFthird* pows->MFthird;
		}
		phasing = (this->build_phase(f,lambda,params,pows,pn_phase_coeffs));
		adouble adub_phasing = (adouble)phasing;
		phase_fixed[j]= -adub_phasing.value();
	}

	//Interpolation to  fix tc
	acc_fixed = gsl_interp_accel_alloc();
	phiI_fixed = gsl_spline_alloc(gsl_interp_cspline,n_fixed);
	gsl_spline_init(phiI_fixed, freqs_fixed, phase_fixed, n_fixed);
	//adouble t_corr_fixed =  gsl_spline_eval_deriv(phiI_fixed, f_final, acc_fixed)/(2.*M_PI) + params->tc;
	adouble t_corr_fixed =  gsl_spline_eval_deriv(phiI_fixed, f_final, acc_fixed)/(2.*M_PI) ;
	//Clean up
	gsl_spline_free(phiI_fixed);
	gsl_interp_accel_free(acc_fixed);
	//#########################################################
	return (adouble)t_corr_fixed;

}
template<class T>
void IMRPhenomPv2<T>::WignerD(T d2[5],T dm2[5], useful_powers<T> *pows, source_parameters<T> *params)
{
	T L = this->L2PN(params->eta, pows) ;//* params->M * params->M;
	T s = params->SP / ( L + params-> SL);
	T s_2 = s*s;

	T cos_beta = 1./sqrt(1.0 + s_2);
	T cos_beta_half = sqrt( (1.0 + cos_beta) / 2.0);
	T sin_beta_half = sqrt( (1.0 - cos_beta) / 2.0);
	T c2 = cos_beta_half * cos_beta_half;
	T s2 = sin_beta_half * sin_beta_half;
	T c3 = c2 * cos_beta_half;
	T s3 = s2 * sin_beta_half;
	T c4 = c3 * cos_beta_half;
	T s4 = s3 * sin_beta_half;
	
	d2[0] = s4;
	d2[1] = 2*cos_beta_half * s3;
	d2[2] = sqrt_6 * s2*c2 ;
	d2[3] = 2 * c3* sin_beta_half;
	d2[4] = c4;
	//Exploit Symmetry
	dm2[0] = d2[4];
	dm2[1] = -d2[3];
	dm2[2] = d2[2];
	dm2[3] = -d2[1];
	dm2[4] = d2[0];
	
}

template<class T>
void IMRPhenomPv2<T>::calculate_twistup( T alpha, std::complex<T> *hp_factor, std::complex<T> *hc_factor, T d2[5], T dm2[5], sph_harm<T> *sph_harm)
{
	std::complex<T> T2m;
	std::complex<T> Tm2m;
	std::complex<T> exp_a = std::exp(std::complex<T>(0,1.) *alpha);
	std::complex<T> exp_ma = std::complex<T>(1.,0.0)/exp_a;
	std::complex<T> exp_2a = exp_a*exp_a;
	std::complex<T> exp_m2a = exp_ma*exp_ma;
	std::complex<T> exp_a_vec[5] = {exp_m2a, exp_ma,std::complex<T>(1.0,0.0) , exp_a, exp_2a};
	std::complex<T> harmonics[5] = 
			{sph_harm->Y2m2, sph_harm->Y2m1, sph_harm->Y20, sph_harm->Y21, sph_harm->Y22};
	for (int m = -2; m<=2; m++)
	{
		T2m = exp_a_vec[-m+2] * dm2[m+2] * harmonics[m+2];
		Tm2m = exp_a_vec[m+2] * d2[m+2] * conj(harmonics[m+2]); 
		*hp_factor += T2m + Tm2m;
		*hc_factor += std::complex<T>(0,1.) * (T2m - Tm2m);
	}
	
}
template<class T>
void IMRPhenomPv2<T>::calculate_euler_angles(T *alpha, T *epsilon, useful_powers<T> *pows, alpha_coeffs<T> *acoeffs, epsilon_coeffs<T> *ecoeffs)
{
	T omega_cbrt = pows->MFthird * pows->PIthird ;
	T omega_cbrt2 = omega_cbrt*omega_cbrt;
	T logomega = log(omega_cbrt2*omega_cbrt);
	T omega = omega_cbrt2*omega_cbrt;
	*alpha = (acoeffs->coeff1/omega
              + acoeffs->coeff2/omega_cbrt2
              + acoeffs->coeff3/omega_cbrt
              + acoeffs->coeff4*logomega
              + acoeffs->coeff5*omega_cbrt);

	*epsilon = (ecoeffs->coeff1/omega
                + ecoeffs->coeff2/omega_cbrt2
                + ecoeffs->coeff3/omega_cbrt
                + ecoeffs->coeff4*logomega
                + ecoeffs->coeff5*omega_cbrt);

}

template<class T>
T IMRPhenomPv2<T>::PhenomPv2_inplane_spin(gen_params_base<T> *params)
{
	//Calculate spin parameters chil and chip
	T chi1_l = params->spin1[2];
	T chi2_l = params->spin2[2];

	T q = params->mass1/params->mass2;
	T chi_eff = (params->mass1*chi1_l + params->mass2*chi2_l) / (params->mass1 + params->mass2); /* Effective aligned spin */
  	T chil = (1.0+q)/q * chi_eff; /* dimensionless aligned spin of the largest BH */
	//params->chil = chil;

	T eta = calculate_eta(params->mass1, params->mass2);
	
	T m1_2 = params->mass1 * params->mass1;
	T m2_2 = params->mass2 * params->mass2;

	T m1 = q/(1+q);
	T m2 = 1./(1+q);

	T S1_perp = m1_2 * sqrt( params->spin1[1]* params->spin1[1] +
				params->spin1[0] * params->spin1[0]);
	T S2_perp = m2_2 * sqrt( params->spin2[1]* params->spin2[1] +
				params->spin2[0] * params->spin2[0]);

	T A1 = 2 + (3*params->mass2)/ ( 2 * params->mass1);
	T A2 = 2 + (3*params->mass1)/ ( 2 * params->mass2);
	T ASp1 = A1*S1_perp;
	T ASp2 = A2*S2_perp;
	T num = (ASp2>ASp1) ? ASp2 : ASp1;
	T denom = (params->mass2 > params->mass1)? A2*m2_2 : A1*m1_2;
	return num/denom;
}

/*!\brief Calculate the unit vector in the direction of the total angular momentum
 * 
 * 
 */
template<class T>
void IMRPhenomPv2<T>::PhenomPv2_JSF_from_params(gen_params_base<T> *params, T *JSF)
{
	//Calculate spin parameters chil and chip
	T chi1_l = params->spin1[2];
	T chi2_l = params->spin2[2];

	T q = params->mass1/params->mass2;
	T chi_eff = (params->mass1*chi1_l + params->mass2*chi2_l) / (params->mass1 + params->mass2); /* Effective aligned spin */
  	T chil = (1.0+q)/q * chi_eff; /* dimensionless aligned spin of the largest BH */
	//params->chil = chil;

	T eta = calculate_eta(params->mass1, params->mass2);
	
	T m1_2 = params->mass1 * params->mass1;
	T m2_2 = params->mass2 * params->mass2;

	T m1 = q/(1+q);
	T m2 = 1./(1+q);
	//params->SP = params->chip * m1*m1;
	//params->SL = chi1_l * m1 * m1 + chi2_l *m2 * m2;
	
	//Compute the rotation operations for L, J0 at fref	
	T L0 = 0.0;
	useful_powers<T> pows;
	IMRPhenomD<T> temp;
	temp.precalc_powers_PI(&pows);
	temp.precalc_powers_ins(params->f_ref, (params->mass1+params->mass2)*MSOL_SEC, &pows);
	
	L0 = pow_int(params->mass1+params->mass2,2) * this->L2PN(eta, &pows);
	
	//_sf denotes source frame - ie L_hat propto z_hat
	T J0x_sf = m1_2*params->chip *cos(params->phip)  ;	
	T J0y_sf = m1_2*params->chip *sin(params->phip)  ;	
	T J0z_sf = L0 + m1_2* chi1_l + m2_2 * chi2_l;
		
	T J0 = sqrt(J0x_sf * J0x_sf + J0y_sf * J0y_sf + J0z_sf * J0z_sf ) ;
	
	JSF[0] =J0x_sf/J0;
	JSF[1] =J0y_sf/J0;
	JSF[2] =J0z_sf/J0;
}
/*! /Brief Parameter transformation to pre-calculate needed parameters for PhenomP from source parameters
 *
 * Pretty much stolen verbatim from lalsuite
 */
template<class T>
void IMRPhenomPv2<T>::PhenomPv2_Param_Transform_reduced(source_parameters<T> *params /*< Source Parameters*/
						)
{
	//Calculate spin parameters chil and chip
	T chi1_l = params->spin1z;
	T chi2_l = params->spin2z;

	T q = params->mass1/params->mass2;
	T chi_eff = (params->mass1*chi1_l + params->mass2*chi2_l) / params->M; /* Effective aligned spin */
  	T chil = (1.0+q)/q * chi_eff; /* dimensionless aligned spin of the largest BH */
	params->chil = chil;


	
	T m1_2 = params->mass1 * params->mass1;
	T m2_2 = params->mass2 * params->mass2;

	T m1 = q/(1+q);
	T m2 = 1./(1+q);
	params->SP = params->chip * m1*m1;
	params->SL = chi1_l * m1 * m1 + chi2_l *m2 * m2;
	
	//Compute the rotation operations for L, J0 at fref	
	T L0 = 0.0;
	useful_powers<T> pows;
	IMRPhenomD<T> temp;
	temp.precalc_powers_PI(&pows);
	temp.precalc_powers_ins(params->f_ref, params->M, &pows);
	
	L0 = params->M * params-> M * this->L2PN(params->eta, &pows);
	
	//_sf denotes source frame - ie L_hat propto z_hat
	T J0x_sf = m1_2*params->chip *cos(params->phip)  ;	
	T J0y_sf = m1_2*params->chip *sin(params->phip)  ;	
	T J0z_sf = L0 + m1_2* params->spin1z + m2_2 * params->spin2z;
	T J0 = sqrt(J0x_sf * J0x_sf + J0y_sf * J0y_sf + J0z_sf * J0z_sf ) ;
	
	//thetaJ_sf is the angle between J0 and L (zhat)
	T thetaJ_sf;
	thetaJ_sf = acos(J0z_sf/J0);

	//azimuthal angle of J0 in the source frame
	T phiJ_sf;
	phiJ_sf = atan2(J0y_sf,J0x_sf); 
	params->phi_aligned = - phiJ_sf;

	//Rotation of the system s.t. the total J is pointed in zhat
	T tmp1,tmp2;
	T incl = params->incl_angle;
	T phiRef = params->phiRef;
	T Nx_sf = sin(incl) * cos(M_PI/2. - phiRef);
	T Ny_sf = sin(incl) * sin(M_PI/2. - phiRef);
	T Nz_sf = cos(incl);
	T tmp_x = Nx_sf;
	T tmp_y = Ny_sf;
	T tmp_z = Nz_sf;
	ROTATEZ(-phiJ_sf, tmp_x,tmp_y, tmp_z);
	ROTATEY(-thetaJ_sf, tmp_x,tmp_y, tmp_z);
	T kappa;
	kappa = -atan2(tmp_y,tmp_x);

	//alpha0
	tmp_x = 0.;
	tmp_y = 0.;
	tmp_z = 1.;
	
	ROTATEZ(-phiJ_sf, tmp_x,tmp_y, tmp_z);
	ROTATEY(-thetaJ_sf, tmp_x,tmp_y, tmp_z);
	ROTATEZ(kappa, tmp_x,tmp_y, tmp_z);
	
	params->alpha0 = atan2(tmp_y,tmp_x);

	tmp_x = Nx_sf;
  	tmp_y = Ny_sf;
  	tmp_z = Nz_sf;
  	ROTATEZ(-phiJ_sf, tmp_x, tmp_y, tmp_z);
  	ROTATEY(-thetaJ_sf, tmp_x, tmp_y, tmp_z);
  	ROTATEZ(kappa, tmp_x, tmp_y, tmp_z);
  	T Nx_Jf = tmp_x; // let's store those two since we will reuse them later (we don't need the y component)
  	T Nz_Jf = tmp_z;
  	params->thetaJN = acos(Nz_Jf);

	/* Finally, we need to redefine the polarizations :
	   PhenomP's polarizations are defined following Arun et al (arXiv:0810.5336)
	   i.e. projecting the metric onto the P,Q,N triad defined with P=NxJ/|NxJ| (see (2.6) in there).
	   By contrast, the triad X,Y,N used in LAL
	   ("waveframe" in the nomenclature of T1500606-v6)
	   is defined in e.g. eq (35) of this document
	   (via its components in the source frame; note we use the defautl Omega=Pi/2).
	   Both triads differ from each other by a rotation around N by an angle \zeta
	   and we need to rotate the polarizations accordingly by 2\zeta
	  */

	T Xx_sf = -cos(incl)*sin(phiRef);
  	T Xy_sf = -cos(incl)*cos(phiRef);
  	T Xz_sf = sin(incl);
  	tmp_x = Xx_sf;
  	tmp_y = Xy_sf;
  	tmp_z = Xz_sf;
  	ROTATEZ(-phiJ_sf, tmp_x, tmp_y, tmp_z);
  	ROTATEY(-thetaJ_sf, tmp_x, tmp_y, tmp_z);
  	ROTATEZ(kappa, tmp_x, tmp_y, tmp_z);
  	//now the tmp_a are the components of X in the J frame
  	//we need the polar angle of that vector in the P,Q basis of Arun et al
  	// P=NxJ/|NxJ| and since we put N in the (pos x)z half plane of the J frame
  	T PArunx_Jf = 0.;
  	T PAruny_Jf = -1.;
  	T PArunz_Jf = 0.;
  	// Q=NxP
  	T QArunx_Jf = Nz_Jf;
  	T QAruny_Jf = 0.;
  	T QArunz_Jf = -Nx_Jf;
  	T XdotPArun = tmp_x*PArunx_Jf+tmp_y*PAruny_Jf+tmp_z*PArunz_Jf;
  	T XdotQArun = tmp_x*QArunx_Jf+tmp_y*QAruny_Jf+tmp_z*QArunz_Jf;
  	params->zeta_polariz = atan2(XdotQArun , XdotPArun);
	
}
/*! /Brief Parameter transformtion to precalculate needed parameters for PhenomP from source parameters -- assumed inclination of total angular momentum J is given, not orbital angular momentum (in source frame (Lhat == zhat)
 *
 * Pretty much stolen verbatim from lalsuite
 */
template<class T>
void IMRPhenomPv2<T>::PhenomPv2_Param_Transform_J(source_parameters<T> *params /*< Source Parameters*/
						)
{
	//T phiJ_sf = params->phiJL;
	//Calculate spin parameters chil and chip
	T chi1_l = params->spin1z;
	T chi2_l = params->spin2z;
	T chip = params->chip;

	T q = params->mass1/params->mass2;
	T chi_eff = (params->mass1*chi1_l + params->mass2*chi2_l) / params->M; /* Effective aligned spin */
  	T chil = (1.0+q)/q * chi_eff; /* dimensionless aligned spin of the largest BH */
	params->chil = chil;


	
	T m1_2 = params->mass1 * params->mass1;
	T m2_2 = params->mass2 * params->mass2;
	
	//T S1_perp = m1_2 * sqrt( params->chi1_p);
	//T S2_perp = m2_2 * sqrt( params->chi2_p);

	//T A1 = 2 + (3*params->mass2)/ ( 2 * params->mass1);
	//T A2 = 2 + (3*params->mass1)/ ( 2 * params->mass2);
	//T ASp1 = A1*S1_perp;
	//T ASp2 = A2*S2_perp;
	//T num = (ASp2>ASp1) ? ASp2 : ASp1;
	//T denom = (params->mass2 > params->mass1)? A2*m2_2 : A1*m1_2;
	//params->chip = num/denom;
	T m1 = q/(1+q);
	T m2 = 1./(1+q);
	params->SP = params->chip * m1*m1;
	params->SL = chi1_l * m1 * m1 + chi2_l *m2 * m2;
	
	//Compute the rotation operations for L, J0 at fref	
	//T L0 = 0.0;
	//useful_powers<T> pows;
	//IMRPhenomD<T> temp;
	//temp.precalc_powers_PI(&pows);
	//temp.precalc_powers_ins(params->f_ref, params->M, &pows);
	//
	//L0 = params->M * params-> M * this->L2PN(params->eta, &pows);
	//
	////_sf denotes source frame - ie L_hat propto z_hat
	////T J0x_sf = m1_2 * params->spin1x + m2_2 * params->spin2x;	
	////T J0y_sf = m1_2 * params->spin1y + m2_2 * params->spin2y;	
	//
	////T J0z_sf = L0 + m1_2* params->spin1z + m2_2 * params->spin2z;
	//	
	//T J0 = sqrt(J0x_sf * J0x_sf + J0y_sf * J0y_sf + J0z_sf * J0z_sf ) ;
	
	//thetaJ_sf is the angle between J0 and L (zhat)
	T thetaJ_sf;
	//thetaJ_sf = acos(J0z_sf/J0);
	thetaJ_sf = params->thetaJL;

	//azimuthal angle of J0 in the source frame
	T phiJ_sf;
	//phiJ_sf = atan2(J0y_sf,J0x_sf); 
	phiJ_sf = params->phiJL;
	params->phi_aligned = - phiJ_sf;

	//Rotation of the system s.t. the total J is pointed in zhat
	T tmp1,tmp2;
	T incl = params->incl_angle;
	T phiRef = params->phiRef;
	T Nx_sf = sin(incl) * cos(M_PI/2. - phiRef);
	T Ny_sf = sin(incl) * sin(M_PI/2. - phiRef);
	T Nz_sf = cos(incl);
	T tmp_x = Nx_sf;
	T tmp_y = Ny_sf;
	T tmp_z = Nz_sf;
	ROTATEZ(-phiJ_sf, tmp_x,tmp_y, tmp_z);
	ROTATEY(-thetaJ_sf, tmp_x,tmp_y, tmp_z);
	T kappa;
	kappa = -atan2(tmp_y,tmp_x);

	//alpha0
	tmp_x = 0.;
	tmp_y = 0.;
	tmp_z = 1.;
	
	ROTATEZ(-phiJ_sf, tmp_x,tmp_y, tmp_z);
	ROTATEY(-thetaJ_sf, tmp_x,tmp_y, tmp_z);
	ROTATEZ(kappa, tmp_x,tmp_y, tmp_z);
	
	params->alpha0 = atan2(tmp_y,tmp_x);

	tmp_x = Nx_sf;
  	tmp_y = Ny_sf;
  	tmp_z = Nz_sf;
  	ROTATEZ(-phiJ_sf, tmp_x, tmp_y, tmp_z);
  	ROTATEY(-thetaJ_sf, tmp_x, tmp_y, tmp_z);
  	ROTATEZ(kappa, tmp_x, tmp_y, tmp_z);
  	T Nx_Jf = tmp_x; // let's store those two since we will reuse them later (we don't need the y component)
  	T Nz_Jf = tmp_z;
  	params->thetaJN = acos(Nz_Jf);

	/* Finally, we need to redefine the polarizations :
	   PhenomP's polarizations are defined following Arun et al (arXiv:0810.5336)
	   i.e. projecting the metric onto the P,Q,N triad defined with P=NxJ/|NxJ| (see (2.6) in there).
	   By contrast, the triad X,Y,N used in LAL
	   ("waveframe" in the nomenclature of T1500606-v6)
	   is defined in e.g. eq (35) of this document
	   (via its components in the source frame; note we use the defautl Omega=Pi/2).
	   Both triads differ from each other by a rotation around N by an angle \zeta
	   and we need to rotate the polarizations accordingly by 2\zeta
	  */

	T Xx_sf = -cos(incl)*sin(phiRef);
  	T Xy_sf = -cos(incl)*cos(phiRef);
  	T Xz_sf = sin(incl);
  	tmp_x = Xx_sf;
  	tmp_y = Xy_sf;
  	tmp_z = Xz_sf;
  	ROTATEZ(-phiJ_sf, tmp_x, tmp_y, tmp_z);
  	ROTATEY(-thetaJ_sf, tmp_x, tmp_y, tmp_z);
  	ROTATEZ(kappa, tmp_x, tmp_y, tmp_z);
  	//now the tmp_a are the components of X in the J frame
  	//we need the polar angle of that vector in the P,Q basis of Arun et al
  	// P=NxJ/|NxJ| and since we put N in the (pos x)z half plane of the J frame
  	T PArunx_Jf = 0.;
  	T PAruny_Jf = -1.;
  	T PArunz_Jf = 0.;
  	// Q=NxP
  	T QArunx_Jf = Nz_Jf;
  	T QAruny_Jf = 0.;
  	T QArunz_Jf = -Nx_Jf;
  	T XdotPArun = tmp_x*PArunx_Jf+tmp_y*PAruny_Jf+tmp_z*PArunz_Jf;
  	T XdotQArun = tmp_x*QArunx_Jf+tmp_y*QAruny_Jf+tmp_z*QArunz_Jf;
  	params->zeta_polariz = atan2(XdotQArun , XdotPArun);
	
}
/*! /Brief Parameter transformtion to precalculate needed parameters for PhenomP from source parameters
 *
 * Pretty much stolen verbatim from lalsuite
 */
template<class T>
void IMRPhenomPv2<T>::PhenomPv2_Param_Transform(source_parameters<T> *params /*< Source Parameters*/
						)
{
	//Calculate spin parameters chil and chip
	T chi1_l = params->spin1z;
	T chi2_l = params->spin2z;

	T q = params->mass1/params->mass2;
	T chi_eff = (params->mass1*chi1_l + params->mass2*chi2_l) / params->M; /* Effective aligned spin */
  	T chil = (1.0+q)/q * chi_eff; /* dimensionless aligned spin of the largest BH */
	params->chil = chil;


	
	T m1_2 = params->mass1 * params->mass1;
	T m2_2 = params->mass2 * params->mass2;
	
	T S1_perp = m1_2 * sqrt( params->spin1y* params->spin1y +
				params->spin1x * params->spin1x);
	T S2_perp = m2_2 * sqrt( params->spin2y* params->spin2y +
				params->spin2x * params->spin2x);

	T A1 = 2 + (3*params->mass2)/ ( 2 * params->mass1);
	T A2 = 2 + (3*params->mass1)/ ( 2 * params->mass2);
	T ASp1 = A1*S1_perp;
	T ASp2 = A2*S2_perp;
	T num = (ASp2>ASp1) ? ASp2 : ASp1;
	T denom = (params->mass2 > params->mass1)? A2*m2_2 : A1*m1_2;
	params->chip = num/denom;
	T m1 = q/(1+q);
	T m2 = 1./(1+q);
	params->SP = params->chip * m1*m1;
	params->SL = chi1_l * m1 * m1 + chi2_l *m2 * m2;
	
	//Compute the rotation operations for L, J0 at fref	
	T L0 = 0.0;
	useful_powers<T> pows;
	IMRPhenomD<T> temp;
	temp.precalc_powers_PI(&pows);
	temp.precalc_powers_ins(params->f_ref, params->M, &pows);
	
	L0 = params->M * params-> M * this->L2PN(params->eta, &pows);
	
	//_sf denotes source frame - ie L_hat propto z_hat
	T J0x_sf = m1_2 * params->spin1x + m2_2 * params->spin2x;	
	T J0y_sf = m1_2 * params->spin1y + m2_2 * params->spin2y;	
	T J0z_sf = L0 + m1_2* params->spin1z + m2_2 * params->spin2z;
		
	T J0 = sqrt(J0x_sf * J0x_sf + J0y_sf * J0y_sf + J0z_sf * J0z_sf ) ;
	
	//thetaJ_sf is the angle between J0 and L (zhat)
	T thetaJ_sf;
	thetaJ_sf = acos(J0z_sf/J0);

	//azimuthal angle of J0 in the source frame
	T phiJ_sf;
	phiJ_sf = atan2(J0y_sf,J0x_sf); 
	params->phi_aligned = - phiJ_sf;

	//Rotation of the system s.t. the total J is pointed in zhat
	T tmp1,tmp2;
	T incl = params->incl_angle;
	T phiRef = params->phiRef;
	T Nx_sf = sin(incl) * cos(M_PI/2. - phiRef);
	T Ny_sf = sin(incl) * sin(M_PI/2. - phiRef);
	T Nz_sf = cos(incl);
	T tmp_x = Nx_sf;
	T tmp_y = Ny_sf;
	T tmp_z = Nz_sf;
	ROTATEZ(-phiJ_sf, tmp_x,tmp_y, tmp_z);
	ROTATEY(-thetaJ_sf, tmp_x,tmp_y, tmp_z);
	T kappa;
	kappa = -atan2(tmp_y,tmp_x);

	//alpha0
	tmp_x = 0.;
	tmp_y = 0.;
	tmp_z = 1.;
	
	ROTATEZ(-phiJ_sf, tmp_x,tmp_y, tmp_z);
	ROTATEY(-thetaJ_sf, tmp_x,tmp_y, tmp_z);
	ROTATEZ(kappa, tmp_x,tmp_y, tmp_z);
	
	params->alpha0 = atan2(tmp_y,tmp_x);

	tmp_x = Nx_sf;
  	tmp_y = Ny_sf;
  	tmp_z = Nz_sf;
  	ROTATEZ(-phiJ_sf, tmp_x, tmp_y, tmp_z);
  	ROTATEY(-thetaJ_sf, tmp_x, tmp_y, tmp_z);
  	ROTATEZ(kappa, tmp_x, tmp_y, tmp_z);
  	T Nx_Jf = tmp_x; // let's store those two since we will reuse them later (we don't need the y component)
  	T Nz_Jf = tmp_z;
  	params->thetaJN = acos(Nz_Jf);

	/* Finally, we need to redefine the polarizations :
	   PhenomP's polarizations are defined following Arun et al (arXiv:0810.5336)
	   i.e. projecting the metric onto the P,Q,N triad defined with P=NxJ/|NxJ| (see (2.6) in there).
	   By contrast, the triad X,Y,N used in LAL
	   ("waveframe" in the nomenclature of T1500606-v6)
	   is defined in e.g. eq (35) of this document
	   (via its components in the source frame; note we use the defautl Omega=Pi/2).
	   Both triads differ from each other by a rotation around N by an angle \zeta
	   and we need to rotate the polarizations accordingly by 2\zeta
	  */

	T Xx_sf = -cos(incl)*sin(phiRef);
  	T Xy_sf = -cos(incl)*cos(phiRef);
  	T Xz_sf = sin(incl);
  	tmp_x = Xx_sf;
  	tmp_y = Xy_sf;
  	tmp_z = Xz_sf;
  	ROTATEZ(-phiJ_sf, tmp_x, tmp_y, tmp_z);
  	ROTATEY(-thetaJ_sf, tmp_x, tmp_y, tmp_z);
  	ROTATEZ(kappa, tmp_x, tmp_y, tmp_z);
  	//now the tmp_a are the components of X in the J frame
  	//we need the polar angle of that vector in the P,Q basis of Arun et al
  	// P=NxJ/|NxJ| and since we put N in the (pos x)z half plane of the J frame
  	T PArunx_Jf = 0.;
  	T PAruny_Jf = -1.;
  	T PArunz_Jf = 0.;
  	// Q=NxP
  	T QArunx_Jf = Nz_Jf;
  	T QAruny_Jf = 0.;
  	T QArunz_Jf = -Nx_Jf;
  	T XdotPArun = tmp_x*PArunx_Jf+tmp_y*PAruny_Jf+tmp_z*PArunz_Jf;
  	T XdotQArun = tmp_x*QArunx_Jf+tmp_y*QAruny_Jf+tmp_z*QArunz_Jf;
  	params->zeta_polariz = atan2(XdotQArun , XdotPArun);
	
}

template<class T>
T IMRPhenomPv2<T>::L2PN( T eta, useful_powers<T> *pow)
{
	T x = pow->MF2third * pow->PI2third;
	T x2 = x*x;
	T eta2 = eta*eta;
	return (eta*(1.0 + (1.5 + eta/6.0)*x + (3.375 - (19.0*eta)/8. - eta2/24.0)*x2)) / sqrt(x);
}

/**
 *  * Wrapper for final-spin formula based on:
 *   * - IMRPhenomD's FinalSpin0815() for aligned spins.
 *    *
 *     * We use their convention m1>m2
 *      * and put <b>all in-plane spin on the larger BH</b>.
 *       *
 *        * In the aligned limit return the FinalSpin0815 value.
 *         */
template <class T>
T IMRPhenomPv2<T>::FinalSpinIMRPhenomD_all_in_plane_spin_on_larger_BH(
  T m1,     /**< Mass of companion 1 (solar masses) */
  T m2,     /**< Mass of companion 2 (solar masses) */
  T chi1_l, /**< Aligned spin of BH 1 */
  T chi2_l, /**< Aligned spin of BH 2 */
  T chip)   /**< Dimensionless spin in the orbital plane */
{
  T M = m1+m2;
  T eta = m1*m2/(M*M);
  T af_parallel, q_factor;
  if (m1 >= m2) {
    q_factor = m1/M;
    af_parallel = this->FinalSpin0815(eta, chi1_l, chi2_l);
  }
  else {
    q_factor = m2/M;
    af_parallel = this->FinalSpin0815(eta, chi2_l, chi1_l);
  }

  T Sperp = chip * q_factor*q_factor;
  T af = copysign_internal((T)1.0, af_parallel) * sqrt(Sperp*Sperp + af_parallel*af_parallel);
  return af;
}

template<class T>
T IMRPhenomPv2<T>::final_spin(source_parameters<T> *params)
{
	return FinalSpinIMRPhenomD_all_in_plane_spin_on_larger_BH(params->mass1,params->mass2,params->spin1z, params->spin2z, params->chip);
}


template class IMRPhenomPv2<double>;
template class IMRPhenomPv2<adouble>;
