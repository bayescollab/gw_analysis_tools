#include "TaylorT2.h"
#include <complex>
#include <cmath>
#include <adolc/adouble.h>
#include <adolc/taping.h>
#include <adolc/drivers/drivers.h>

#ifndef _OPENMP
#define omp ignore
#endif



/*! \brief Constructs the waveform as outlined by 
 *
 * arguments:
 * 	array of frequencies, length of that array, a complex array for the output waveform, and a source_parameters structure
 */
template <class T>
int TaylorT2<T>::construct_waveform(T *times, /**< T array of frequencies the waveform is to be evaluated at*/
	int length, /**< integer length of the array of frequencies and the waveform*/	
	std::complex<T> *hplus,/**< complex T array for the waveform to be output*/ 
	std::complex<T> *hcross,/**< complex T array for the waveform to be output*/ 
	source_parameters<T> *params /*Structure of source parameters to be initialized before computation*/
	)
{
	T omega0 = pow(params->x0,3./2.)/params->M;
	T mu = params->mass1 * params->mass2 / params->M;
	for (int i = 0 ; i<length ; i++){
		T theta = THETA(times[i], params);
		T x = X(theta, params);
		T phi = PHI(x, params);
		T gamma = GAMMA(x,params);
		T MADM = M_ADM(gamma, params);
		//T omega = pow(x,3./2.) / MADM;
		T omega = pow(x,3./2.)/params->M ;
		T psi = phi - 2 * MADM * omega * log(omega/omega0);
		contract_H(x,psi,params, &hplus[i], &hcross[i]) ;
		hplus[i]*= 2*mu * x / params->DL;
		hcross[i]*= 2*mu * x / params->DL;
	}
	return 1;
}

template<class T>
T TaylorT2<T>::GAMMA(T x, source_parameters<T> *params)
{
	T eta = params->eta;
	T out =   1 ;
	out+= ( 1 - eta/3.) * x ;
	//out+= ( 1. - 65./12. * eta ) *x * x;
	//out+= (1 + ( -2203. / 2520. - 41. / 192. * M_PI * M_PI - 22./3. * log(r/r0))*eta + 229./36. * eta * eta + 1./81. * eta*eta*eta) * x *x *x;
	out*= x;
	return out;
}

template<class T> 
void TaylorT2<T>::contract_H(T x, T psi, source_parameters<T> *params, std::complex<T> *hplus, std::complex<T> *hcross)
{
	T ci = cos(params->incl_angle);
	T si = sin(params->incl_angle);
	T cp = cos(2.*psi);
	T sp = sin(2.*psi);

	T rootx = sqrt(x);

	*hplus = 0;
	*hplus += (-(1. + ci*ci) * cp - 1. / 96. * si*si * ( 17. + ci*ci))  ;//p = 0

	*hcross = 0;
	*hcross += -2. * ci * sp  ;//p = 0

}

template<class T>
T TaylorT2<T>::M_ADM(T gamma, source_parameters<T> *params)
{
	T eta = params->eta;
	//return params->M*(1 - eta/2. * gamma + eta / 8. * ( 7 - eta) * gamma * gamma);
	return params->M*(1 - eta/2. * gamma );
}

/*! \brief phi function from 318 1310.1528 
 */
template<class T>
T TaylorT2<T>::PHI(T x , source_parameters<T> *params)
{
	T x_2 = sqrt(x);
	T x_3_2 = x * x_2;
	T x2 = x*x;
	T x_5_2 = x_2 * x2;
	T x3 = x*x*x;
	T x_7_2 = x_2 * x3;
	T eta = params->eta;

	T out = 0;
	out += 1;
	out+= ( 3715./1008. + 55./12. * eta) * x;
        // out+= -10. * M_PI * x_3_2;
        // out += ( 15293365. / 1016064. + 27145. / 1008. * eta + 3085. / 144. *
        // eta * eta ) * x2; out += ( 38645. / 1344. - 65. / 16. * eta) * M_PI *
        // x_5_2 * log(x/params->x0); out += ( 12348611926451. / 18776862720. -
        // 160. / 3. * M_PI * M_PI  - 1712./21. * kGammaE - 856. / 21. *
        // log(16.*x) + ( -15737765635. / 12192768. + 2255. / 48. * M_PI * M_PI
        // ) * eta + 76055./6912. * eta * eta - 127825. / 5184. * eta * eta *
        // eta) * x3; out += (77096675. / 2032128. + 378515. / 12096. * eta -
        // 74045./6048. * eta * eta) * M_PI * x_7_2;
        out *= - 1./ 32. / eta / x_5_2;
	return out;	
}

/*! \brief x function from 316 1310.1528 
 */
template<class T>
T TaylorT2<T>::X(T theta , source_parameters<T> *params)
{
	T theta_8 = pow(theta, 1./8.);
	T theta_4 = theta_8*theta_8;
	T theta_3_8 = theta_4*theta_8;
	T theta_2 = theta_4*theta_4;
	T theta_5_8 = theta_2*theta_8;
	T theta_3_4 = theta_5_8*theta_8;
	T theta_7_8 = theta_3_4*theta_8;
	T eta = params->eta;
	
	T out = 0;
	out+= 1;
	out+= (743./4032. + 11./48. * eta) /theta_4;
        // out+= (-1./5.) *M_PI / theta_3_8;
        // out+= (19583./254016.  + 24401./193536. * eta + 31./288. * eta*eta )
        // / theta_2; out += (-11891./53760. + 109./1920. * eta) *M_PI /
        // theta_5_8; out += (-10052469856691./6008596070400. + 1./6. * M_PI *
        // M_PI + 107./420. * kGammaE - 107./3360. * log(theta/256.) +
        // (3147553127./780337152. - 451. / 3072. *M_PI * M_PI) * eta -
        // 15211./442368. * eta*eta + 25565./331776. * eta * eta * eta) /
        // theta_3_4; out += (-113868647./433520640. - 31821./143360. * eta +
        // 294941./3870720. * eta * eta ) * M_PI / theta_7_8;
        out *= 1./4. * 1./theta_4 ;

	//T term1 = 
	return out;
}

/*! \brief Theta function from 315 1310.1528 
 */
template<class T>
T TaylorT2<T>::THETA(T time, source_parameters<T> *params)
{
	T out = params->eta * (-1)*time / params->M/ 5;
	return out;
}



template class TaylorT2<double>;
template class TaylorT2<adouble>;
