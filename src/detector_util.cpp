#include "detector_util.h"
#include "util.h"
#include "pn_waveform_util.h"
#include "GWATConfig.h"
#include "io_util.h"
#include <fstream>
#include <iostream>
#include <string>
#include <math.h>
#include <adolc/adouble.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_integration.h>

/*! \file
 * Routines to construct noise curves for various detectors and for detector specific utilities for response functions and coordinate transformations
 */

/*! \brief Function to populate the squareroot of the noise curve for various detectors
 *
 * If frequencies are left as NULL, standard frequency spacing is applied and the frequencies are returned, in which case the frequencies argument becomes an output array
 *
 * Detector names must be spelled exactly
 *
 * Detectors include: 
 * 	
 * 	aLIGO_analytic -- analytic approximation of the advanced LIGO sensitivity curve
 *
 * 	Hanford_O1_fitted -- Fitted function to the O1 noise curve for Hanford
 * 		
 * 	AdLIGOMidHigh -- from Emanuele Berti
 *
 * 	LISA -- LISA sensitivity curve with out sky averaging for a single channel 
 * 		
 * 		-- does NOT include the factor of sqrt(3)/2, which is taken care of in the response function for LISA
 *
 * 	LISA_CONF -- LISA sensitivity curve with out sky averaging for a single channel including confusion noise from the galatic white dwarf population
 *
 * 		-- does NOT include the factor of sqrt(3)/2, which is taken care of in the response function for LISA
 *
 * 	LISA_SADC -- LISA sensitivity curve with sky averaging for a dual channel 
 *
 * 		-- does include the factor of sqrt(3)/2
 * 		
 * 		-- This DOES include the sky averaged pattern functions -- i.e. sky_average flag should NOT be used when making the waveform
 *
 * 	LISA_SADC_CONF -- LISA sensitivity curve with sky averaging for a dual channel including confusion noise from the galatic white dwarf population
 *
 * 		-- does include the factor of sqrt(3)/2
 *
 * 		-- This DOES include the sky averaged pattern functions -- i.e. sky_average flag should NOT be used when making the waveform
 *
 * 	AdLIGOAPlus -- AdLIGO A+ O5
 *
 * 	AdLIGOAPlus_smoothed -- AdLIGO A+ O5 -- without spectral lines
 *
 * 	AdLIGODesign -- Design sensitivity for LIGO O4
 *
 * 	AdLIGODesign_smoothed -- Design sensitivity for LIGO O4 -- without spectral lines
 *
 * 	CE1 -- Cosmic Explorer phase 1 -- circa 2035
 *
 * 	CE2 -- Cosmic Explorer phase 2 -- circa 2045
 *
 * 	KAGRA_opt -- KAGRA O5 05 optimistic
 *
 * 	KAGRA_pess -- KAGRA O5 05 pessimistic
 *
 * 	AdVIRGOPlus2_opt -- phase 2, O5, optimistic
 *
 * 	AdVIRGOPlus2_opt_smoothed -- phase 2, O5, optimistic -- without spectral lines
 *
 * 	AdVIRGOPlus2_pess -- phase 2, O5, pessimistic
 *
 * 	AdVIRGOPlus2_pess_smoothed -- phase 2, O5, pessimistic -- without spectral lines
 *
 * 	AdVIRGOPlus1 -- phase 1, O4 analog ( rescaled O5 curve according to BBH distance estimates in https://dcc.ligo.org/public/0161/P1900218/002/SummaryForObservers.pdf)
 *
 * 	AdVIRGOPlus1_smoothed -- phase 1, O4 analog ( rescaled O5 curve according to BBH distance estimates in https://dcc.ligo.org/public/0161/P1900218/002/SummaryForObservers.pdf) -- without spectral lines
 *
 * 	ET-D -- Einstein telescope with lines included. From ET website
 *
 * 	ET-D_smooth -- Einstein telescope without lines. From Emanuele 
 *
 * 	AdLIGOVoyager -- Voyager configuration of Advanced LIGO
 *
 */
void populate_noise(double *frequencies, /**< double array of frquencies (NULL)*/
		std::string detector, /**< String to designate the detector noise curve to be used */
		double *noise_root, /**< ouptput double array for the square root of the PSD of the noise of the specified detector*/
		int length,/**< integer length of the output and input arrays*/
		double integration_time /**< Integration time in months (only important for LISA_conf*/
		)
{
	//std::string currently_supported_dir = std::string(GWAT_ROOT_DIRECTORY)+"data/noise_data/currently_supported/";
	std::string currently_supported_dir = std::string(GWAT_SHARE_DIR)+"noise_data/";
	//std::string currently_supported_dir = std::string(GWAT_ROOT_DIRECTORY)+"data/noise_data/currently_supported_raw/";
	//std::cout<<currently_supported_dir<<std::endl;
	if(check_list(detector, (std::string *)analytic_PSD_models,analytic_PSD_models_N)){
		if(detector == "aLIGO_analytic")
		{
			if(!frequencies )
			{
				int len = 1000;
				for (int i =0; i<len;i++)
				{
					frequencies[i] = 10.+ i;
					noise_root[i] = aLIGO_analytic(frequencies[i]);
				}
			}
			
			else
			{
				for (int i =0; i<length;i++)
				{
					noise_root[i] = aLIGO_analytic(frequencies[i]);
				}
			}
		}
		else if(detector == "Hanford_O1_fitted")
		{
			if(!frequencies)
			{
				int len = 1000;
				for (int i =0; i<len;i++)
				{
					frequencies[i] = 10.+ i;
					noise_root[i] = Hanford_O1_fitted(frequencies[i]);
				}
			}
			
			else
			{
				for (int i =0; i<length;i++)
				{
					noise_root[i] = Hanford_O1_fitted(frequencies[i]);
				}
			}
		}
		else if(detector =="LISA_SADC" ){
			for(int i = 0 ; i<length; i++){
				noise_root[i] = sqrt(LISA_analytic_SADC(frequencies[i]));
			}
		}
		else if(detector == "LISA_SADC_CONF"){
			double alpha, beta, kappa, gamma, fk;
			sort_LISA_SC_coeffs(&alpha, &beta, &kappa,  &gamma, &fk, integration_time);
			for(int i = 0 ; i<length; i++){
				noise_root[i] = sqrt(LISA_analytic_SADC(frequencies[i]) +  LISA_SC(frequencies[i], alpha, beta, kappa,gamma, fk));
			}
		}
		else if(detector =="LISA" ){
			for(int i = 0 ; i<length; i++){
				noise_root[i] = sqrt(LISA_analytic(frequencies[i]));
			}
		}
		else if(detector == "LISA_CONF"){
			double alpha, beta, kappa, gamma, fk;
			sort_LISA_SC_coeffs(&alpha, &beta, &kappa,  &gamma, &fk, integration_time);
			for(int i = 0 ; i<length; i++){
				noise_root[i] = sqrt(LISA_analytic(frequencies[i]) +  LISA_SC(frequencies[i], alpha, beta, kappa,gamma, fk));
			}
		}
	}
	else if(check_list(detector, (std::string *)interp_PSD_models,interp_PSD_models_N)){
		int dat_length=0;
		std::string file;
		if(detector=="AdLIGOMidHigh"){
			//dat_length = 3000;
			dat_length = 3453;
			file = "AdLIGOMidHigh.csv";
		}
		if(detector=="_AdLIGODesign"){
			dat_length = 3000;
			file = "AdLIGODesign.csv";
		}
		else if(detector=="AdLIGODesign"){
			dat_length = 2736;
			file = "aligo_O4high.csv";
		}
		else if(detector=="AdLIGODesign_smoothed"){
			dat_length = 2736;
			file = "aligo_O4high_smoothed.csv";
		}
		else if(detector=="AdLIGOAPlus"){
			dat_length = 3000;
			file = "AplusDesign.csv";
		}
		else if(detector=="AdLIGOAPlus_smoothed"){
			dat_length = 3000;
			file = "AplusDesign_smoothed.csv";
		}
		else if(detector=="CE1"){
			dat_length = 3000;
			file = "CE1_strain.csv";
		}
		else if(detector=="CE1_smoothed"){
			dat_length = 3000;
			file = "CE1_strain_smoothed.csv";
		}
		else if(detector=="CE2"){
			dat_length = 3000;
			file = "CE2_strain.csv";
		}
		else if(detector=="CE2_smoothed"){
			dat_length = 3000;
			file = "CE2_strain_smoothed.csv";
		}
		else if(detector=="AdVIRGOPlus2_opt"){
			dat_length = 3000;
			file = "avirgo_O5high_NEW.csv";
		}
		else if(detector=="AdVIRGOPlus2_opt_smoothed"){
			dat_length = 3000;
			file = "avirgo_O5high_NEW_smoothed.csv";
		}
		else if(detector=="AdVIRGOPlus2_pess"){
			dat_length = 3000;
			file = "avirgo_O5low_NEW.csv";
		}
		else if(detector=="AdVIRGOPlus2_pess_smoothed"){
			dat_length = 3000;
			file = "avirgo_O5low_NEW_smoothed.csv";
		}
		else if(detector=="AdVIRGOPlus1"){
			dat_length = 3000;
			file = "avirgo_O4high_NEW.csv";
		}
		else if(detector=="AdVIRGOPlus1_smoothed"){
			dat_length = 3000;
			file = "avirgo_O4high_NEW_smoothed.csv";
		}
		else if(detector=="KAGRA_opt"){
			dat_length = 1000;
			file = "kagra_128Mpc.csv";
		}
		else if(detector=="KAGRA_pess"){
			dat_length = 1000;
			file = "kagra_80Mpc.csv";
		}
		else if(detector=="ET-D"){
			dat_length = 3000;
			file = "ET-0000A-18_ETDSensitivityCurveTxtFile.csv";
		}
		else if(detector=="ET-D_smoothed"){
			dat_length = 3000;
			file = "ETDXylophoneDwyer.csv";
		}
		else if(detector=="AdLIGOVoyager"){
			dat_length = 3000;
			file = "Voyager.csv";
		}
		count_lines_data_file(currently_supported_dir+file, &dat_length);
		gsl_interp_accel *accel = gsl_interp_accel_alloc();
		gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear, dat_length);
		double **data = new double*[dat_length];
		for(int i = 0 ;i<dat_length; i++){
			data[i] = new double[2];
		}
		read_file(currently_supported_dir + file,data, dat_length, 2);
		double psd[dat_length];
		double f[dat_length];
		for(int i = 0 ;i<dat_length; i++){
			psd[i]=data[i][1];	
			f[i]=data[i][0];	
		}
		gsl_spline_init(spline, f, psd, dat_length);
		for(int i = 0 ; i<length; i++){
			noise_root[i]=gsl_spline_eval(spline, frequencies[i],accel);
		}	
		gsl_spline_free(spline);
		gsl_interp_accel_free(accel);
		for(int i = 0 ; i<dat_length; i++){
			delete [] data[i];
		}
		delete [] data;
	}
	else{
		std::cout<<"Detector "<<detector<<" not supported"<<std::endl;
	}
			
	return;	
}

/*! \brief Analytic function approximating the PSD for aLIGO 
 *
 * CITE (Will?)
 */
double aLIGO_analytic(double f)
{
	double S= 3e-48;
	double fknee = 70.;
	double x = fknee/f;
	double x4 = x*x*x*x;
	return sqrt( S * (x4 + 2 + 2*x*x)/5 );
}
/*! \brief Analytic function approximating the PSD for LISA sensitivity curve -- this is S, not root S and does not sky average, treats the 2 channels separately, and the geometrical factor of \sqrt{3}/2 is included in the waveform.
 *
 * NON Sky averaged -  single channel
 *
 * arXiv:1905.08811 -- equation 14-17
 */
double LISA_analytic(double f)
{
	double L = 2.5 * pow_int(10.,9);
	double fstar = 19.09*pow_int(10.,-3);
	double twopi = 2.*M_PI;
	double POMS = LISA_POMS(f);
	double PACC = LISA_PACC(f);
	double S = 1./(  L*L) * ( POMS +2.*(1. + pow_int(  std::cos(f/fstar) ,2) )*(  PACC/pow_int( twopi * f, 4) ))*(1. + 6./10. * pow_int(f/fstar,2)) ;
	return  S;
}
/*! \brief Analytic function approximating the PSD for LISA sensitivity curve -- this is S, not root S 
 *
 * Sky averaged -  dual channel
 *
 * arXiv:1803.01944 -- equation 13
 *
 * NOTE: may need to divide by \sqrt{3}/2.. My LISA response functions already include this factor -- factor 0f 3./4. (sqrt(3)/2)**2
 *
 * Keeping the factor of sqrt(3)/2, since this will probably be used without the beam patter functions (that's where the factor of sqrt(3)/2 usually appears)
 */
double LISA_analytic_SADC(double f)
{
	double L = 2.5 * pow_int(10.,9);
	double fstar = 19.09*pow_int(10.,-3);
	double twopi = 2.*M_PI;
	double POMS = LISA_POMS(f);
	double PACC = LISA_PACC(f);
	//double S = (3./4.)*10./(3. * L*L) * ( POMS +2.*(1. + pow_int(  std::cos(f/fstar) ,2) )*(  PACC/pow_int( twopi * f, 4) ))*(1. + 6./10. * pow_int(f/fstar,2)) ;
	double S = 10./(3. * L*L) * ( POMS +2.*(1. + pow_int(  std::cos(f/fstar) ,2) )*(  PACC/pow_int( twopi * f, 4) ))*(1. + 6./10. * pow_int(f/fstar,2)) ;
	return  S;
}
/*! \breif Optical metrology noise function
 *
 * arXiv:1803.01944 -- equation 10
 */
double LISA_POMS(double f)
{
	double factor = 1.5 * pow_int(10., -11);
	return factor * factor * (1. + pow_int( (.002/f) ,4) );
}

/*! \breif Single test mass accelartion noise
 *
 * arXiv:1803.01944 -- equation 11
 */
double LISA_PACC(double f)
{
	double factor = 3.* pow_int(10.,-15);
	return factor * factor * (1. + pow_int(.0004/f,2) ) * ( 1. + pow_int(f/.008,4));

}

/*! \breif Confusion noise  -- this is S, not root  S
 *
 * arXiv:1803.01944 -- 14
 *
 * integration_time is the observation time in months -- options are 6, 12, 24, and 48
 */
double LISA_SC(double f, double alpha, double beta, double kappa, double gamma, double fk)
{
	double A = 9.* pow_int(10.,-45);
	double SC = A * pow(f,-7./3.) * exp( -pow(f,alpha) + beta * f * sin(kappa * f )) * ( 1. + tanh(gamma * (fk - f)));
	
	return  SC;
}

/*! \breif LISA confusion noise coefficients 
 *
 * arXiv:1803.01944 -- Table 1
 *
 * integration_time is the observation time in months -- options are 6, 12, 24, and 48
 *
 */
void sort_LISA_SC_coeffs(double *alpha, double *beta, double *kappa, double *gamma, double *fk, double  integration_time)
{
	if(integration_time ==6){
		*alpha = 0.133;
		*beta = 243.;
		*kappa = 482.;
		*gamma = 917.;
		*fk = .00258;
	}
	else if(integration_time ==12){
		*alpha = 0.171;
		*beta = 292.;
		*kappa = 1020.;
		*gamma = 1680.;
		*fk = .00215;
	}
	else if(integration_time ==24){
		*alpha = 0.165;
		*beta = 299.;
		*kappa = 611.;
		*gamma = 1340.;
		*fk = .00173;
	}
	else if(integration_time >=48){
		*alpha = 0.138;
		*beta = -221.;
		*kappa = 521.;
		*gamma = 1680.;
		*fk = .00113;
	}

}

/*! \brief Numerically fit PSD to the Hanford Detector's O1
 *
 * CITE (Yunes?)
 */
double Hanford_O1_fitted(double f)
{
	double avec[7]={47.8466,-92.1896,35.9273,-7.61447,0.916742,-0.0588089,0.00156345};
	double S0 = .8464;
	double x = log(f);
	return sqrt(S0) * exp( avec[0] + avec[1]*x + avec[2]*x*x +
            avec[3] * x*x*x + avec[4]*x*x*x*x + avec[5]*x*x*x*x*x + avec[6]*x*x*x*x*x*x);
}

/*! \brief Utility for the overall amplitude and phase shift for spin-aligned systems
 *
 * For spin aligned, all the extrinsic parameters have the effect of an overall amplitude modulation and phase shift
 */
std::complex<double> Q(double theta, double phi, double iota, double psi)
{
	double ct = cos(theta);
	double cp2 = cos(2.*phi);
	double sp2 = sin(2.*phi);
	double ci = cos(iota);
	double cpsi2 = cos(2*psi);
	double spsi2 = sin(2*psi);

	double Fplus = (1./2)*(1+ ct*ct)*cp2;
	double Fcross = ct * sp2;
	//std::complex<double> Q = (1+ci*ci)/2. *Fplus + std::complex<double>(0,Fcross*ci);
	std::complex<double> Q= (1+ci*ci)/2. *(Fplus*cpsi2 - Fcross*spsi2) + std::complex<double>(0,(Fcross *cpsi2 + Fplus*spsi2)*ci);
	return Q;
}

/*! \brief Utility for the overall amplitude and phase shift for spin-aligned systems
 *
 * For spin aligned, all the extrinsic parameters have the effect of an overall amplitude modulation and phase shift
 */
std::complex<double> Q(double theta, double phi, double iota)
{
	std::complex<double> Qout= Q(theta, phi, iota, 0.);
	return Qout;
}


/*! \brief Response function of a 90 deg interferometer for longitudinal (scalar) polarization
 *
 * Theta, phi, and psi are local, horizontal coordinates relative to the detector
 */
template<class T>
T right_interferometer_l(T theta, T phi, T psi)
{
	return 1;
}

/*! \brief Response function of a 90 deg interferometer for breathing (scalar) polarization
 *
 * Theta, phi, and psi are local, horizontal coordinates relative to the detector
 */
template<class T>
T right_interferometer_b(T theta, T phi, T psi)
{
	return 1;
}

/*! \brief Response function of a 90 deg interferometer for y (vector) polarization
 *
 * Theta, phi, and psi are local, horizontal coordinates relative to the detector
 */
template<class T>
T right_interferometer_y(T theta, T phi, T psi)
{
	return 1;
}

/*! \brief Response function of a 90 deg interferometer for x (vector) polarization
 *
 * Theta, phi, and psi are local, horizontal coordinates relative to the detector
 */
template<class T>
T right_interferometer_x(T theta, T phi, T psi)
{
	return 1;
}

/*! \brief Response function of a 90 deg interferometer for plus polarization
 *
 * Theta and phi are local, horizontal coordinates relative to the detector
 */
template<class T>
T right_interferometer_plus(T theta, T phi)
{
	T ct = cos(theta);
	T cp2 = cos(2.*phi);

	T Fplus = (1./2)*(1+ ct*ct)*cp2;
	return Fplus;
}
/*! \brief Response function of a 90 deg interferometer for cross polarization
 *
 * Theta and phi are local, horizontal coordinates relative to the detector
 */
template<class T>
T right_interferometer_cross(T theta, T phi)
{
	T ct = cos(theta);
	T sp2 = sin(2.*phi);
	T Fcross = ct * sp2;
	return Fcross;
}
/*! \brief Response function of a 90 deg interferometer 
 *
 * Theta and phi are local, horizontal coordinates relative to the detector
 *
 * With general psi
 */
template<class T>
void right_interferometer(det_res_pat<T> *r_pat,T theta, T phi, T psi)
{
	T Fplus = right_interferometer_plus(theta,phi);
	T Fcross = right_interferometer_cross(theta,phi);
	T c2psi = cos(2.*psi);
	T s2psi = sin(2.*psi);
	*(r_pat->Fplus) = Fplus*c2psi - Fcross*s2psi;
	*(r_pat->Fcross) = Fplus*s2psi + Fcross*c2psi;

	//###############################################
	//Stand in for TESTING MUST BE FIXED
	if(r_pat->active_polarizations[2]){
		//vector x
		*(r_pat->Fx) = right_interferometer_x(theta,phi,psi);
	}	
	if(r_pat->active_polarizations[3]){
		//vector y
		*(r_pat->Fy) = right_interferometer_y(theta,phi,psi);
	}	
	if(r_pat->active_polarizations[4]){
		//breathing
		*(r_pat->Fb) = right_interferometer_b(theta,phi,psi);
	}	
	if(r_pat->active_polarizations[5]){
		//longitudinal 
		*(r_pat->Fl) = right_interferometer_l(theta,phi,psi);
	}	
	//###############################################
}


/*! \brief Transform from celestial coordinates to local horizontal coords
 *
 * (RA,DEC) -> (altitude, azimuth)
 *
 * Need gps_time of transformation, as the horizontal coords change in time
 *
 * detector is used to specify the lat and long of the local frame
 */
template<class T>
void celestial_horizon_transform(T RA, /**< in RAD*/
		T DEC, /**< in RAD*/
		double gps_time, 
		std::string detector, 
		T *phi, /**< in RAD*/
		T *theta /**< in RAD*/
		)
{
	T LAT, LONG, azimuth_offset;
	if (detector =="Hanford" || detector == "hanford")
	{
		LAT =  H_LAT;		
		LONG =  H_LONG;		
		azimuth_offset = H_azimuth_offset;
		
	}
	else if (detector =="Livingston" || detector == "livingston")
	{
		LAT =  L_LAT;		
		LONG =  L_LONG;		
		azimuth_offset = L_azimuth_offset;

	}
	else if (detector =="Virgo" || detector == "virgo")
	{
		LAT =  V_LAT;		
		LONG =  V_LONG;		
		azimuth_offset = V_azimuth_offset;

	}
	else if (detector =="Kagra" || detector == "kagra")
	{
		LAT =  K_LAT;		
		LONG =  K_LONG;		
		azimuth_offset = K_azimuth_offset;

	}
	else if (detector =="Indigo" || detector == "indigo")
	{
		LAT =  I_LAT;		
		LONG =  I_LONG;		
		azimuth_offset = I_azimuth_offset;

	}
	else if (detector =="Cosmic Explorer" || detector == "cosmic explorer" || detector=="CE")
	{
		LAT =  CE_LAT;		
		LONG =  CE_LONG;		
		azimuth_offset = CE_azimuth_offset;

	}
	else if (detector =="Einstein Telescope 1" || detector == "einstein telescope 1" || detector=="ET1")
	{
		LAT =  ET1_LAT;		
		LONG =  ET1_LONG;		
		azimuth_offset = ET1_azimuth_offset;

	}
	else if (detector =="Einstein Telescope 2" || detector == "einstein telescope 2" || detector=="ET2")
	{
		LAT =  ET2_LAT;		
		LONG =  ET2_LONG;		
		azimuth_offset = ET2_azimuth_offset;
	}
	else if (detector =="Einstein Telescope 3" || detector == "einstein telescope 3" || detector=="ET3")
	{
		LAT =  ET2_LAT;		
		LONG =  ET2_LONG;		
		azimuth_offset = ET2_azimuth_offset;
	}
	else {
		std::cout<<"Invalid detector"<<std::endl;
		exit(1);
	}
	celestial_horizon_transform(RA,DEC, gps_time, LONG, LAT, phi, theta);
	*phi += azimuth_offset;
	if(*phi>2*M_PI) *phi -=2*M_PI;
	if(*phi<2*M_PI) *phi +=2*M_PI;
}

/*! \brief Numerical derivative of the transformation
 *
 * Planned for use in Fisher calculations, but not currently implemented anywhere
 */
void derivative_celestial_horizon_transform(double RA, /**< in RAD*/
		double DEC, /**< in RAD*/
		double gps_time, 
		std::string detector, 
		double *dphi_dRA, 
		double *dtheta_dRA, 
		double *dphi_dDEC, 
		double *dtheta_dDEC 
		)
{
	double phip, phim, thetap, thetam;
	double epsilon = 1e-5;

	celestial_horizon_transform(RA+epsilon, DEC,gps_time,detector,&phip,&thetap);
	celestial_horizon_transform(RA-epsilon, DEC,gps_time,detector,&phim,&thetam);
	*dtheta_dRA = (thetap-thetam)/(2*epsilon);
	*dphi_dRA = (phip-phim)/(2*epsilon);

	celestial_horizon_transform(RA, DEC+epsilon,gps_time,detector,&phip,&thetap);
	celestial_horizon_transform(RA, DEC-epsilon,gps_time,detector,&phim,&thetam);
	*dtheta_dDEC = (thetap-thetam)/(2*epsilon);
	*dphi_dDEC = (phip-phim)/(2*epsilon);
}

/*! \brief calculate difference in time of arrival (DTOA) for a given source location and 2 different detectors
 *
 * Full version, from LAL
 */
template <class T>
T DTOA_DETECTOR(T RA, /**< spherical polar angle for detector 1 in RAD*/
	T DEC, /**<spherical polar angle for detector 2 in RAD*/ 
	double GMST_rad,/**< Greenwich mean sidereal time of detection*/
	std::string detector1, /**< name of detector one*/
	std::string detector2 /**<name of detector two*/
	)
{
	const double *earth_centered_location1;
	const double *earth_centered_location2;
	//detector one
	if(detector1 == "Hanford" || detector1 == "hanford")
	{
		earth_centered_location1 =H_location ;	
	}
	else if(detector1 == "Livingston" || detector1 == "livingston")
	{
		earth_centered_location1 =L_location ;	
	}
	else if(detector1 == "Virgo" || detector1 == "virgo")
	{
		earth_centered_location1 =V_location ;	
	}
	else if(detector1 == "Kagra" || detector1 == "kagra")
	{
		earth_centered_location1 =K_location ;	
	}
	else if(detector1 == "Indigo" || detector1 == "indigo")
	{
		earth_centered_location1 =I_location ;	
	}
	else if(detector1 == "Cosmic Explorer" || detector1 == "cosmic explorer" || detector1=="CE")
	{
		earth_centered_location1 =CE_location ;	
	}
	else if(detector1 == "Einstein Telescope 1" || detector1 == "einstein telescope 1" || detector1=="ET1")
	{
		earth_centered_location1 =ET1_location ;	
	}
	else if(detector1 == "Einstein Telescope 2" || detector1 == "einstein telescope 2" || detector1=="ET2")
	{
		earth_centered_location1 =ET2_location ;	
	}
	else if(detector1 == "Einstein Telescope 3" || detector1 == "einstein telescope 3" || detector1=="ET3")
	{
		earth_centered_location1 =ET3_location ;	
	}
	//detector 2
	if(detector2 == "Hanford" || detector2 == "hanford")
	{
		earth_centered_location2 =H_location ;	
	}
	else if(detector2 == "Livingston" || detector2 == "livingston")
	{
		earth_centered_location2 =L_location ;	
	}
	else if(detector2 == "Virgo" || detector2 == "virgo")
	{
		earth_centered_location2 =V_location ;	
	}
	else if(detector2 == "Kagra" || detector2 == "kagra")
	{
		earth_centered_location2 =K_location ;	
	}
	else if(detector2 == "Indigo" || detector2 == "indigo")
	{
		earth_centered_location2 =I_location ;	
	}
	else if(detector2 == "Cosmic Explorer" || detector2 == "cosmic explorer" || detector2=="CE")
	{
		earth_centered_location2 =CE_location ;	
	}
	else if(detector2 == "Einstein Telescope 1" || detector2 == "einstein telescope 1" || detector2=="ET1")
	{
		earth_centered_location2 =ET1_location ;	
	}
	else if(detector2 == "Einstein Telescope 2" || detector2 == "einstein telescope 2" || detector2=="ET2")
	{
		earth_centered_location2 =ET2_location ;	
	}
	else if(detector2 == "Einstein Telescope 3" || detector2 == "einstein telescope 3" || detector2=="ET3")
	{
		earth_centered_location2 =ET3_location ;	
	}
	return DTOA_earth_centered_coord(RA,DEC,GMST_rad,earth_centered_location1,earth_centered_location2);
}
template double DTOA_DETECTOR<double>(double,double,double,std::string,std::string);
template adouble DTOA_DETECTOR<adouble>(adouble,adouble,double,std::string,std::string);
/*! \brief calculate difference in time of arrival (DTOA) for a given source location and 2 different detectors
 */
template <class T>
T DTOA_earth_centered_coord(T RA, /**< spherical polar angle for detector 1 in RAD*/
	T DEC, /**<spherical polar angle for detector 2 in RAD*/ 
	double GMST_rad,/**< Greenwich mean sidereal time of detection*/
	const double *loc1, /**< Location of the first detector in Earth centered coordinates in meters*/
	const double *loc2 /**<Location of the second detector in Earth centered coordinates in meters*/
	)
{
	T dx[3];
	dx[0] = loc1[0]-loc2[0];
	dx[1] = loc1[1]-loc2[1];
	dx[2] = loc1[2]-loc2[2];
	
	//Direction to source
	T hour_angle = GMST_rad - RA;
	T ehat[3];
	ehat[0] = cos(DEC) * cos(hour_angle);
	ehat[1] = cos(DEC) * -sin(hour_angle);
	ehat[2] = sin(DEC) ;
		
	return (dx[0]*ehat[0] + dx[1]*ehat[1] + dx[2]*ehat[2])/c;

}
template double DTOA_earth_centered_coord<double>(double,double,double,const double *,const double *);
template adouble DTOA_earth_centered_coord<adouble>(adouble,adouble,double,const double *,const double *);
/*! \brief calculate difference in time of arrival (DTOA) for a given source location and 2 different detectors
 */
double DTOA(double theta1, /**< spherical polar angle for detector 1 in RAD*/
	double theta2, /**<spherical polar angle for detector 2 in RAD*/ 
	std::string detector1, /**< name of detector one*/
	std::string detector2 /**<name of detector two*/
	)
{
	double R1, R2;
	//detector one
	if(detector1 == "Hanford" || detector1 == "hanford")
	{
		R1 = H_radius;	
	}
	else if(detector1 == "Livingston" || detector1 == "livingston")
	{
		R1 = L_radius;	
	}
	else if(detector1 == "Virgo" || detector1 == "virgo")
	{
		R1 = V_radius;	
	}
	else if(detector1 == "Kagra" || detector1 == "kagra")
	{
		R1 = K_radius;	
	}
	else if(detector1 == "Indigo" || detector1 == "indigo")
	{
		R1 = I_radius;	
	}
	else if(detector1 == "Cosmic Explorer" || detector1 == "cosmic explorer" || detector1=="CE")
	{
		R1 = CE_radius;	
	}
	else if(detector1 == "Einstein Telescope 1" || detector1 == "einstein telescope 1" || detector1=="ET1")
	{
		R1 = ET1_radius;	
	}
	else if(detector1 == "Einstein Telescope 2" || detector1 == "einstein telescope 2" || detector1=="ET2")
	{
		R1 = ET2_radius;	
	}
	else if(detector1 == "Einstein Telescope 3" || detector1 == "einstein telescope 3" || detector1=="ET3")
	{
		R1 = ET3_radius;	
	}
	//detector 2
	if(detector2 == "Hanford" || detector2 == "hanford")
	{
		R2 = H_radius;	
	}
	else if(detector2 == "Livingston" || detector2 == "livingston")
	{
		R2 = L_radius;	
	}
	else if(detector2 == "Virgo" || detector2 == "virgo")
	{
		R2 = V_radius;	
	}
	else if(detector2 == "Kagra" || detector2 == "kagra")
	{
		R2 = K_radius;	
	}
	else if(detector2 == "Indigo" || detector2 == "indigo")
	{
		R2 = I_radius;	
	}
	else if(detector2 == "Cosmic Explorer" || detector2 == "cosmic explorer" || detector2=="CE")
	{
		R2 = CE_radius;	
	}
	else if(detector2 == "Einstein Telescope 1" || detector2 == "einstein telescope 1" || detector2=="ET1")
	{
		R2 = ET1_radius;	
	}
	else if(detector2 == "Einstein Telescope 2" || detector2 == "einstein telescope 2" || detector2=="ET2")
	{
		R2 = ET2_radius;	
	}
	else if(detector2 == "Einstein Telescope 3" || detector2 == "einstein telescope 3" || detector2=="ET3")
	{
		R2 = ET3_radius;	
	}
	return (cos(theta1)*R1 - cos(theta2)*R2)/c;
}

/*! /brief Analytic approximation of the radius from the center of earth to a given location
 *
 * Just the raidus as a function of angles, modelling an oblate spheroid
 */
double radius_at_lat(double latitude, /**< latitude in degrees*/
			double elevation /**<elevation in meters*/
			)
{
	double numerator = pow(RE_equatorial*RE_equatorial * cos(latitude*M_PI/180.) ,2 ) 
			+ pow(RE_polar*RE_polar * sin(latitude*M_PI/180. ), 2);
	double denominator = pow(RE_equatorial * cos(latitude*M_PI/180.) ,2 ) 
			+ pow(RE_polar * sin(latitude*M_PI/180. ), 2);
	return sqrt(numerator/denominator) + elevation;
}

/*! \brief Calculates the response coefficients for a detector with response tensor D for a source at RA, Dec, and psi
 *
 * Taken from LALSuite
 *
 * The response tensor for each of the operational detectors is precomputed in detector_util.h, but to create a new tensor, follow the outline in Anderson et al 36 PRD 63 042003 (2001) Appendix B
 *
 * For terrestial detectors -- psi is the polarization angle from a detector at earth's center, aligned with equatorial coordinates
 */
template<class T>
void detector_response_functions_equatorial(double D[3][3],/**< Detector Response tensor (3x3)*/
	T ra,/**<Right ascension in rad*/
	T dec,/**<Declination in rad*/
	T psi,/**< polarization angle in rad*/
	double gmst,/**<Greenwich mean sidereal time (rad)*/
	det_res_pat<T> *r_pat
	)
{
	int i;
	T X[3];
	T Y[3];
	T Z[3];
	
	/* Greenwich hour angle of source (radians). */
	T gha = gmst - ra;
	/* pre-compute trig functions */
	T cosgha = cos(gha);
	T singha = sin(gha);
	T cosdec = cos(dec);
	T sindec = sin(dec);
	T cospsi = cos(psi);
	T sinpsi = sin(psi);
	
	/* Eq. (B4) of [ABCF].  Note that dec = pi/2 - theta, and gha =
	 * -phi where theta and phi are the standard spherical coordinates
	 * used in that paper. */
	X[0] = -cospsi * singha - sinpsi * cosgha * sindec;
	X[1] = -cospsi * cosgha + sinpsi * singha * sindec;
	X[2] =  sinpsi * cosdec;
	
	/* Eq. (B5) of [ABCF].  Note that dec = pi/2 - theta, and gha =
	 * -phi where theta and phi are the standard spherical coordinates
	 * used in that paper. */
	Y[0] =  sinpsi * singha - cospsi * cosgha * sindec;
	Y[1] =  sinpsi * cosgha + cospsi * singha * sindec;
	Y[2] =  cospsi * cosdec;

/* Eqns from [Blaut2012] - but converted from Blaut's theta = pi/2 - dec
 * given cos(dec) = sin(theta) and cos(theta) = sin(dec), and Blaut's phi = ra
 * - gmst, so cos(phi) = cos(gha) and sin(phi) = -sin(gha). This is consistent
 * with the convention in XLALComputeDetAMResponse.
 */
	Z[0] = -cosdec * cosgha;
	Z[1] = cosdec * singha;
	Z[2] = -sindec;
	
	/* Now compute Eq. (B7) of [ABCF] for each polarization state, i.e.,
	 * with s+=1 and sx=0 to get F+, with s+=0 and sx=1 to get Fx */
	*(r_pat->Fplus)= 0.0;
	*(r_pat->Fcross)= 0.0;
	if(r_pat->active_polarizations[2]){
		*(r_pat->Fx) = 0.0;
	}
	if(r_pat->active_polarizations[3]){
		*(r_pat->Fy) = 0.0;
	}
	if(r_pat->active_polarizations[4]){
		*(r_pat->Fb) = 0.0;
	}
	if(r_pat->active_polarizations[5]){
		*(r_pat->Fl) = 0.0;
	}
	for(i = 0; i < 3; i++) {
	        T DX = D[i][0] * X[0] + D[i][1] * X[1] + D[i][2] * X[2];
	        T DY = D[i][0] * Y[0] + D[i][1] * Y[1] + D[i][2] * Y[2];
		T DZ = D[i][0] * Z[0] + D[i][1] * Z[1] + D[i][2] * Z[2];
	        *(r_pat->Fplus ) += X[i] * DX - Y[i] * DY;
	        *(r_pat->Fcross) += X[i] * DY + Y[i] * DX;
		if(r_pat->active_polarizations[2]){
			*(r_pat->Fx) += X[i] * DZ + Z[i] * DX;
		}
		if(r_pat->active_polarizations[3]){
			*(r_pat->Fy) += Y[i] * DZ + Z[i] * DY;
		}
		if(r_pat->active_polarizations[4]){
			*(r_pat->Fb) += X[i] * DX + Y[i] * DY;
		}
		if(r_pat->active_polarizations[5]){
			*(r_pat->Fl) += Z[i] * DZ;
		}
	}

	//if(r_pat->active_polarizations[2]){
	//	//vector x
	//	*(r_pat->Fx) = 0.0;
	//	for(i = 0; i < 3; i++) {
	//		*(r_pat->Fx) += X[i] * DZ + Z[i] * DX;
	//	}
	//}	
	//if(r_pat->active_polarizations[3]){
	//	//vector y
	//	*(r_pat->Fy) = 0.0;
	//	for(i = 0; i < 3; i++) {
	//		*(r_pat->Fy) += Y[i] * DZ + Z[i] * DY;
	//	}
	//}	
	//if(r_pat->active_polarizations[4]){
	//	//breathing
	//	*(r_pat->Fb) = 0.0;
	//	for(i = 0; i < 3; i++) {
	//		*(r_pat->Fb) += X[i] * DX + Y[i] * DY;
	//	}
	//}	
	//if(r_pat->active_polarizations[5]){
	//	//longitudinal 
	//	*(r_pat->Fl) = 0.0;
	//	for(i = 0; i < 3; i++) {
	//		*(r_pat->Fl) += Z[i] * DZ;
	//	}
	//}	
	//###############################################

}
/*! \brief Wrapping of the equatorial detector response for terrestial based detectors
 *
 * For ground based detectors, the antenna pattern functions are not functions of time.
 */
template <class T>
void detector_response_functions_equatorial(std::string detector,/**< Detector */
	T ra,/**<Right ascension in rad*/
	T dec,/**<Declination in rad*/
	T psi,/**< polarization angle in rad*/
	double gmst,/**<Greenwich mean sidereal time (rad)*/
	det_res_pat<T> *r_pat
	)
{
	//FIX -- dec goes frorm pi/2 to -pi/2 while theta goes from 0 to pi
	detector_response_functions_equatorial(detector,ra,dec,psi,gmst,(T *)NULL,0,(T)0.,(T)0.,(T)0.,(T)0., r_pat);
}
/*! \brief Same as the other function, but for active and future detectors
 *
 * If terrestial detectors, it will use ra, dec, psi, and gmst
 *
 * If space detector, it will use ra, dec transformed to ecliptic coord, and LISA_alpha0, LISA_phi0 and theta_j_ecl, phi_j_ecl, which are the offsets for LISA and the spherical angles for the unit vector in the total angular momemntum in the ecliptic system
 */
template<class T>
void detector_response_functions_equatorial(std::string detector,/**< Detector */
	T ra,/**<Right ascension in rad*/
	T dec,/**<Declination in rad*/
	T psi,/**< polarization angle in rad*/
	double gmst,/**<Greenwich mean sidereal time (rad)*/
	T *times,/**<Times at which to evaluate Fplus and Fcross, in which case the Fplus and Fcross pointers are arrays*/
	int length,/**< Length of the arrays*/
	T LISA_alpha0,/**< Offset for alpha*/
	T LISA_phi0,/**<Offset for phi*/
	T theta_j_ecl,/**< Ecliptic spherical polar angle of Lhat (Jhat for precessing systems)*/
	T phi_j_ecl,/**< Ecliptic sherical azimuthal angle (Jhat  for precessing systems)*/
	det_res_pat<T> *r_pat
	)
{
	double geometric_factor=1;
	//Time dependent antenna patterns
	if(detector=="LISA" ||detector =="lisa"){
		//Polar angle theta runs 0,PI , dec runs -PI/2,PI/2
		//Ecliptic  coord
		T theta_s;
		T phi_s;
		//M_PI/2 - dec is polar angle in equatorial
		//RA is azimuthal angle in equatorial
		ecl_from_eq((T)( M_PI/2. - dec), ra, &theta_s, &phi_s);
		for(int i =0 ; i<length; i++){
			r_pat->Fplus[i]  = LISA_response_plus_time(theta_s,phi_s, theta_j_ecl,phi_j_ecl,LISA_alpha0,LISA_phi0, times[i]);
			r_pat->Fcross[i]  = LISA_response_cross_time(theta_s,phi_s, theta_j_ecl,phi_j_ecl,LISA_alpha0,LISA_phi0, times[i]);
		}
		
		if(r_pat->active_polarizations[2]){
			//vector x
		}	
		if(r_pat->active_polarizations[3]){
			//vector y
		}	
		if(r_pat->active_polarizations[4]){
			//breathing
		}	
		if(r_pat->active_polarizations[5]){
			//longitudinal 
		}	
	}
	//Time independent response functions
	else if (detector=="static") // Generic static detector
	{
		// See Eqs. (B9)-(B10) of 10.1103/PhysRevD.63.042003
		// Pre-compute common factors
		T theta = M_PI_2 - dec;
		T phi = ra - gmst;
		T cosTheta = cos(theta);
		T cosTwoPsi = cos(2*psi);
		T sinTwoPsi = sin(2*psi);
		// 1/2*(1+cos(\theta)^2)*cos(2\phi)
		T fact1 = 0.5*(1+cosTheta*cosTheta)*cos(2*phi);
		// -cos(\theta)*sin(2\phi)
		T fact2 = -cosTheta*sin(2*phi);
		// Compute the antenna patterns
		*(r_pat->Fplus) = -fact1*cosTwoPsi + fact2*sinTwoPsi;
		*(r_pat->Fcross) = fact1*sinTwoPsi + fact2*cosTwoPsi;

		// Extra polarizations not yet implemented
		if (r_pat->active_polarizations[2]\
			|| r_pat->active_polarizations[3]\
			|| r_pat->active_polarizations[4] \
			|| r_pat->active_polarizations[5])
		{
			std::cerr << "ERROR: Extra polarizations not available for this detector.\n";
			exit(1);
		}
		
	}
	else{
		double responseM[3][3];
		if(detector =="Hanford" || detector=="hanford"){
			geometric_factor = H_geometric_factor;
			for(int i =0; i<3; i++){
				for(int j =0 ;j<3; j++){
					responseM[i][j] = Hanford_D[i][j];
				}
			}
		}
		else if(detector =="Livingston" || detector=="livingston"){
			geometric_factor = L_geometric_factor;
			for(int i =0; i<3; i++){
				for(int j =0 ;j<3; j++){
					responseM[i][j] = Livingston_D[i][j];
				}
			}
		}
		else if(detector =="Virgo" || detector=="virgo"){
			geometric_factor = V_geometric_factor;
			for(int i =0; i<3; i++){
				for(int j =0 ;j<3; j++){
					responseM[i][j] = Virgo_D[i][j];
				}
			}
		}
		else if(detector =="Kagra" || detector=="kagra"){
			geometric_factor = K_geometric_factor;
			for(int i =0; i<3; i++){
				for(int j =0 ;j<3; j++){
					responseM[i][j] = Kagra_D[i][j];
				}
			}
		}
		else if(detector =="Indigo" || detector=="indigo"){
			geometric_factor = I_geometric_factor;
			for(int i =0; i<3; i++){
				for(int j =0 ;j<3; j++){
					responseM[i][j] = Indigo_D[i][j];
				}
			}
		}
		else if(detector =="Cosmic Explorer" || detector=="cosmic explorer" || detector=="CE"){
			geometric_factor = CE_geometric_factor;
			for(int i =0; i<3; i++){
				for(int j =0 ;j<3; j++){
					responseM[i][j] = CE_D[i][j];
				}
			}
		}
		else if(detector =="Einstein Telescope 1" || detector=="einstein telescope 1" || detector=="ET1"){
			geometric_factor = ET1_geometric_factor;
			for(int i =0; i<3; i++){
				for(int j =0 ;j<3; j++){
					responseM[i][j] = ET1_D[i][j];
				}
			}
		}
		else if(detector =="Einstein Telescope 2" || detector=="einstein telescope 2" || detector=="ET2"){
			geometric_factor = ET2_geometric_factor;
			for(int i =0; i<3; i++){
				for(int j =0 ;j<3; j++){
					responseM[i][j] = ET2_D[i][j];
				}
			}
		}
		else if(detector =="Einstein Telescope 3" || detector=="einstein telescope 3" || detector=="ET3"){
			geometric_factor = ET3_geometric_factor;
			for(int i =0; i<3; i++){
				for(int j =0 ;j<3; j++){
					responseM[i][j] = ET3_D[i][j];
				}
			}
		}
		else{
			std::cout<<"ERROR -- unsupported detector"<<std::endl;
			exit(1);
		}
		detector_response_functions_equatorial(responseM, ra, dec, psi, gmst, r_pat);
		(*(r_pat->Fplus))*=geometric_factor;
		(*(r_pat->Fcross))*=geometric_factor;
		if(r_pat->active_polarizations[2]){
			//vector x
			(*(r_pat->Fx))*=geometric_factor;
		}	
		if(r_pat->active_polarizations[3]){
			//vector y
			(*(r_pat->Fy))*=geometric_factor;
		}	
		if(r_pat->active_polarizations[4]){
			//breathing
			(*(r_pat->Fb))*=geometric_factor;
		}	
		if(r_pat->active_polarizations[5]){
			//longitudinal 
			(*(r_pat->Fl))*=geometric_factor;
		}	
	}
}

template<class T>
T LISA_response_plus(source_parameters<T> *params, T theta_s, T phi_s, T theta_j, T phi_j, T alpha_0, T phi_0, T f)
{
	T t = t_2PN(f, params->eta, params->chirpmass, params->spin1z, params->spin2z, params->tc);
	return LISA_response_plus_time(theta_s, phi_s, theta_j, phi_j, alpha_0, phi_0, t);
}
template<class T>
T LISA_response_cross(source_parameters<T> *params, T theta_s, T phi_s, T theta_j, T phi_j, T alpha_0, T phi_0, T f)
{
	T t = t_2PN(f, params->eta, params->chirpmass, params->spin1z, params->spin2z, params->tc);
	return LISA_response_cross_time(theta_s, phi_s, theta_j, phi_j, alpha_0, phi_0, t);
}

/*! \brief Time dependent detector response of LISA for non-precessing waveforms
 *
 * See https://arxiv.org/abs/gr-qc/0411129 or https://arxiv.org/abs/gr-qc/9703068
 *
 * All the arguments are ``barred'', using the notation in these two works. That is, they are relative to the solar system  barycenter.
 *
 * To get the second interferometer's response, evaluate with phi_j - pi/4.
 *
 * All the coordinates are in the ecliptic coordinate system
 */
template<class T>
T LISA_response_plus_time( T theta_s, T phi_s, T theta_j, T phi_j, T alpha_0, T phi_0, T t)
{
	//T_year defined in include/util.h
	T phi_t = phi_0 + 2. * M_PI * t / T_year;
	T alpha1 = alpha_0;
	T out = (0.1e1 + pow(cos(theta_s) / 0.2e1 - sqrt(0.3e1) * sin(theta_s) * cos(-phi_t + phi_s) / 0.2e1, 0.2e1)) * cos((T) (2 * alpha1) + 0.3141592654e1 / 0.6e1 - 0.2e1 * atan((sqrt(0.3e1) * cos(theta_s) + sin(theta_s) * cos(-phi_t + phi_s)) / sin(theta_s) / sin(-phi_t + phi_s) / 0.2e1)) * cos(0.2e1 * atan((-(cos(theta_j) * cos(theta_s) + sin(theta_j) * sin(theta_s) * cos(phi_j - phi_s)) * (-cos(theta_s) / 0.2e1 + sqrt(0.3e1) * sin(theta_s) * cos(-phi_t + phi_s) / 0.2e1) + cos(theta_j) / 0.2e1 - sqrt(0.3e1) * sin(theta_j) * cos(-phi_t + phi_j) / 0.2e1) / (sin(theta_j) * sin(theta_s) * sin(phi_j - phi_s) / 0.2e1 - sqrt(0.3e1) * cos(phi_t) * (cos(theta_j) * sin(theta_s) * sin(phi_s) - cos(theta_s) * sin(theta_j) * sin(phi_j)) / 0.2e1 - sqrt(0.3e1) * sin(phi_t) * (cos(theta_s) * sin(theta_j) * cos(phi_j) - cos(theta_j) * sin(theta_s) * cos(phi_s)) / 0.2e1))) / 0.2e1 - (cos(theta_s) / 0.2e1 - sqrt(0.3e1) * sin(theta_s) * cos(-phi_t + phi_s) / 0.2e1) * sin((T) (2 * alpha1) + 0.3141592654e1 / 0.6e1 - 0.2e1 * atan((sqrt(0.3e1) * cos(theta_s) + sin(theta_s) * cos(-phi_t + phi_s)) / sin(theta_s) / sin(-phi_t + phi_s) / 0.2e1)) * sin(0.2e1 * atan((-(cos(theta_j) * cos(theta_s) + sin(theta_j) * sin(theta_s) * cos(phi_j - phi_s)) * (-cos(theta_s) / 0.2e1 + sqrt(0.3e1) * sin(theta_s) * cos(-phi_t + phi_s) / 0.2e1) + cos(theta_j) / 0.2e1 - sqrt(0.3e1) * sin(theta_j) * cos(-phi_t + phi_j) / 0.2e1) / (sin(theta_j) * sin(theta_s) * sin(phi_j - phi_s) / 0.2e1 - sqrt(0.3e1) * cos(phi_t) * (cos(theta_j) * sin(theta_s) * sin(phi_s) - cos(theta_s) * sin(theta_j) * sin(phi_j)) / 0.2e1 - sqrt(0.3e1) * sin(phi_t) * (cos(theta_s) * sin(theta_j) * cos(phi_j) - cos(theta_j) * sin(theta_s) * cos(phi_s)) / 0.2e1)));
	//Factor of sqrt(3)/2 for the equilateral triangle -- see Sec III B of gr-qc/9703068
	return ROOT_THREE/2.*out;
}
template<class T>
T LISA_response_cross_time( T theta_s, T phi_s, T theta_j, T phi_j, T alpha_0, T phi_0, T t)
{
	//T_year defined in include/util.h
	T phi_t = phi_0 + 2. * M_PI * t / T_year;
	T alpha1 = alpha_0;
	T out = (0.1e1 + pow(cos(theta_s) / 0.2e1 - sqrt(0.3e1) * sin(theta_s) * cos(-phi_t + phi_s) / 0.2e1, 0.2e1)) * cos((T) (2 * alpha1) + 0.3141592654e1 / 0.6e1 - 0.2e1 * atan((sqrt(0.3e1) * cos(theta_s) + sin(theta_s) * cos(-phi_t + phi_s)) / sin(theta_s) / sin(-phi_t + phi_s) / 0.2e1)) * sin(0.2e1 * atan((-(cos(theta_j) * cos(theta_s) + sin(theta_j) * sin(theta_s) * cos(phi_j - phi_s)) * (-cos(theta_s) / 0.2e1 + sqrt(0.3e1) * sin(theta_s) * cos(-phi_t + phi_s) / 0.2e1) + cos(theta_j) / 0.2e1 - sqrt(0.3e1) * sin(theta_j) * cos(-phi_t + phi_j) / 0.2e1) / (sin(theta_j) * sin(theta_s) * sin(phi_j - phi_s) / 0.2e1 - sqrt(0.3e1) * cos(phi_t) * (cos(theta_j) * sin(theta_s) * sin(phi_s) - cos(theta_s) * sin(theta_j) * sin(phi_j)) / 0.2e1 - sqrt(0.3e1) * sin(phi_t) * (cos(theta_s) * sin(theta_j) * cos(phi_j) - cos(theta_j) * sin(theta_s) * cos(phi_s)) / 0.2e1))) / 0.2e1 + (cos(theta_s) / 0.2e1 - sqrt(0.3e1) * sin(theta_s) * cos(-phi_t + phi_s) / 0.2e1) * sin((T) (2 * alpha1) + 0.3141592654e1 / 0.6e1 - 0.2e1 * atan((sqrt(0.3e1) * cos(theta_s) + sin(theta_s) * cos(-phi_t + phi_s)) / sin(theta_s) / sin(-phi_t + phi_s) / 0.2e1)) * cos(0.2e1 * atan((-(cos(theta_j) * cos(theta_s) + sin(theta_j) * sin(theta_s) * cos(phi_j - phi_s)) * (-cos(theta_s) / 0.2e1 + sqrt(0.3e1) * sin(theta_s) * cos(-phi_t + phi_s) / 0.2e1) + cos(theta_j) / 0.2e1 - sqrt(0.3e1) * sin(theta_j) * cos(-phi_t + phi_j) / 0.2e1) / (sin(theta_j) * sin(theta_s) * sin(phi_j - phi_s) / 0.2e1 - sqrt(0.3e1) * cos(phi_t) * (cos(theta_j) * sin(theta_s) * sin(phi_s) - cos(theta_s) * sin(theta_j) * sin(phi_j)) / 0.2e1 - sqrt(0.3e1) * sin(phi_t) * (cos(theta_s) * sin(theta_j) * cos(phi_j) - cos(theta_j) * sin(theta_s) * cos(phi_s)) / 0.2e1)));
	//Factor of sqrt(3)/2 for the equilateral triangle -- see Sec III B of gr-qc/9703068 
	return ROOT_THREE/2.*out;
}
/*! \brief Utility to calculate the cumulative amplitude distribution for a single detector
 *
 * P(\omega) = \int_V \Theta(\omega'(\Omega, \psi, \iota)-\omega) d\Omega d\psi dcos\iota
 * 
 * Integrated over the volume which \omega' is larger than \omega
 *
 * Integrates using Monte Carlo integration
 *
 * Uniform sampling in \psi, cos(\iota), cos(\theta) , and \phi
 */
double p_single_detector(double omega, /**< \omega = \rho/\rho_opt**/
	int samples/**< number of monte carlo samples to use**/
	)
{
	gsl_rng_env_setup();
	const gsl_rng_type *T = gsl_rng_default;
	gsl_rng *r = gsl_rng_alloc(T);
	gsl_rng_set(r, 1323);
	double omega_prime,omega_prime_squared,Fplus, Fcross,psi, cosiota, 
		costheta, phi,theta,iota;
	double twopi = 2.*M_PI;
	double omega_squared = omega*omega;
	double sum = 0;
	bool active_polar[6] = {true, true, false, false, false, false};
	for(int i= 0 ; i<samples ; i++){
		psi = gsl_rng_uniform(r)*twopi;	
		cosiota = gsl_rng_uniform(r)*2.-1.;	
		costheta = gsl_rng_uniform(r)*2.-1.;	
		phi = gsl_rng_uniform(r)*twopi;	
		iota = acos(cosiota);
		theta = acos(costheta);
		det_res_pat<double> r_pat;
		
		r_pat.Fplus = &Fplus;
		r_pat.Fcross = &Fcross;
		r_pat.active_polarizations = &active_polar[0];
		right_interferometer(&r_pat, theta, phi, psi);
		//Do squared to avoid the expensive sqrt function
		omega_prime_squared =(pow_int( 1. +cosiota*cosiota,2)/4.* Fplus*Fplus + cosiota*cosiota*Fcross*Fcross);	
		if(omega_prime_squared>omega_squared){sum++;}
		
	}	
	gsl_rng_free(r);
	return sum/=samples;
}
/*! \brief Utility to calculate the cumulative amplitude distribution for three detectors
 *
 * P(\omega) = \int_V \Theta(\omega'(\Omega, \psi, \iota)-\omega) d\Omega d\psi dcos\iota
 * 
 * Integrated over the volume which \omega' is larger than \omega
 *
 * Integrates using Monte Carlo integration
 *
 * Uniform sampling in \psi, cos(\iota), cos(DEC) , and RA
 *
 * Sampled using any supported detectors. See documentation for supported detectors
 *
 * Rough upper limit:
 * 	HLV: 1.4
 * 	HLVK: 1.5
 * 	HLVI: 1.65
 * 	HLVKI: 1.7
 */
double p_N_detector(double omega, /**< \omega = \rho/\rho_opt**/
	int samples,/**< number of monte carlo samples to use**/
	int N_detectors, /**< Number of detectors*/
	std::string *detectors,/**<String name of detectors to use -- MUST BE SUPPORTED DETECTORS*/
	int rand_seed /**< Seed for random number generator*/
	)
{
	double gmst = 1;//Doesn't matter -- RA is uniform random anyway
	gsl_rng_env_setup();
	const gsl_rng_type *T = gsl_rng_default;
	gsl_rng *r = gsl_rng_alloc(T);
	gsl_rng_set(r, rand_seed);
	double omega_prime,omega_prime_squared,Fplus, Fcross,psi, cosiota, 
		sinDEC, RA,DEC,iota;
	double twopi = 2.*M_PI;
	double omega_squared = omega*omega;
	double sum = 0;
	det_res_pat<double> r_pat;
	r_pat.Fplus = &Fplus;	
	r_pat.Fcross = &Fcross;	
	bool active_polar[6] = {true, true, false, false, false, false};
	r_pat.active_polarizations = &active_polar[0];
	for(int i= 0 ; i<samples ; i++){
		psi = gsl_rng_uniform(r)*twopi;	
		cosiota = gsl_rng_uniform(r)*2.-1.;	
		sinDEC = gsl_rng_uniform(r)*2.-1.;	
		RA = gsl_rng_uniform(r)*twopi;	
		iota = acos(cosiota);
		DEC = asin(sinDEC);
		omega_prime_squared=0;
		for(int j = 0 ; j<N_detectors; j++){	
				
			detector_response_functions_equatorial(detectors[j],RA,DEC,psi,gmst, &r_pat);	
			//right_interferometer(&Fplus, &Fcross, theta, phi, psi);
			//Do squared to avoid the expensive sqrt function
			omega_prime_squared +=(pow_int( 1. +cosiota*cosiota,2)/4.* Fplus*Fplus + cosiota*cosiota*Fcross*Fcross);	
		}
		if(omega_prime_squared>omega_squared){sum++;}
		
	}	
	gsl_rng_free(r);
	return sum/=samples;
}
/*! \brief Utility to calculate the cumulative amplitude distribution for a single detector -- Numerical Fit 
 *
 * P(\omega) = \int_V \Theta(\omega'(\Omega, \psi, \iota)-\omega) d\Omega d\psi dcos\iota
 * 
 * Integrated over the volume which \omega' is larger than \omega
 *
 * see arXiv:1405.7016
 */
double p_single_detector_fit(double omega /**< \omega = \rho/\rho_opt**/
	)
{
	double alpha=1., a2 = 0.374222, a4 = 2.04216, a8 =-2.63948;
	double omega_term = 1-omega/alpha;
	return (a2 * pow_int( omega_term, 2) 
		+ a4* pow_int(omega_term,4) 
		+ a8* pow_int(omega_term,8)
		+ (1 - a2 - a4 - a8) * pow_int(omega_term, 10));
}
/*! \brief Utility to calculate the cumulative amplitude distribution for triple detector network -- Numerical Fit 
 *
 * P(\omega) = \int_V \Theta(\omega'(\Omega, \psi, \iota)-\omega) d\Omega d\psi dcos\iota
 * 
 * Integrated over the volume which \omega' is larger than \omega
 *
 * see arXiv:1405.7016
 *
 * Not accurate, apparently. Just use interpolation
 */
double p_triple_detector_fit(double omega /**< \omega = \rho/\rho_opt**/
	)
{
	double alpha=1.4, a2 =  1.19549, a4 =  1.61758, a8 =-4.87024;
	double omega_term = 1-omega/alpha;
	return (a2 * pow_int( omega_term, 2) 
		+ a4* pow_int(omega_term,4) 
		+ a8* pow_int(omega_term,8)
		+ (1 - a2 - a4 - a8) * pow_int(omega_term, 10));
}
/*! \brief Utility to calculate the cumulative amplitude distribution for triple detector network -- interpolated data from https://pages.jh.edu/~eberti2/research/ 
 *
 * P(\omega) = \int_V \Theta(\omega'(\Omega, \psi, \iota)-\omega) d\Omega d\psi dcos\iota
 * 
 * Integrated over the volume which \omega' is larger than \omega
 *
 * see arXiv:1405.7016
 */
double p_triple_detector_interp(double omega /**< \omega = \rho/\rho_opt**/
	)
{
	double pomega = 0;	
	int data_length = 141;
	gsl_interp_accel *accel  = gsl_interp_accel_alloc();
	gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline,data_length);

	double omegas[data_length];
	double pomega_numerical[data_length];
	double **temp =  new double*[data_length];
	for(int i=0 ; i<data_length; i ++){
		temp[i]=new double[2];
	}
	//read_file(std::string(GWAT_ROOT_DIRECTORY)+"data/detection_probability_tables/Pw_three.csv", temp, data_length,2);
	read_file( std::string(GWAT_SHARE_DIR)+"detection_probability_tables/Pw_three.csv", temp, data_length,2);
	;
	for(int i=0 ; i<data_length; i ++){
		omegas[i] = temp[i][0];	
		pomega_numerical[i] = temp[i][1];	
	}
	
	gsl_spline_init(spline, omegas, pomega_numerical,data_length);

	pomega = gsl_spline_eval(spline, omega, accel);
	gsl_spline_free(spline);
	gsl_interp_accel_free(accel);
	for(int i = 0 ; i<data_length; i++){
		delete [] temp[i];
	}
	delete[] temp;
	return pomega;
}
/*! \brief Utility to calculate the cumulative amplitude distribution for single detector network -- interpolated data from https://pages.jh.edu/~eberti2/research/ 
 *
 * P(\omega) = \int_V \Theta(\omega'(\Omega, \psi, \iota)-\omega) d\Omega d\psi dcos\iota
 * 
 * Integrated over the volume which \omega' is larger than \omega
 *
 * see arXiv:1405.7016
 */
double p_single_detector_interp(double omega /**< \omega = \rho/\rho_opt**/
	)
{
	double pomega = 0;	
	int data_length = 1001;
	gsl_interp_accel *accel  = gsl_interp_accel_alloc();
	gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline,data_length);

	double omegas[data_length];
	double pomega_numerical[data_length];
	double **temp =  new double*[data_length];
	for(int i=0 ; i<data_length; i ++){
		temp[i]=new double[2];
	}
	//read_file(std::string(GWAT_ROOT_DIRECTORY)+"data/detection_probability_tables/Pw_single.csv", temp, data_length,2);
	read_file(std::string(GWAT_SHARE_DIR)+"detection_probability_tables/Pw_single.csv", temp, data_length,2);
	for(int i=0 ; i<data_length; i ++){
		omegas[i] = temp[i][0];	
		pomega_numerical[i] = temp[i][1];	
	}
	
	gsl_spline_init(spline, omegas, pomega_numerical,data_length);

	pomega = gsl_spline_eval(spline, omega, accel);
	gsl_spline_free(spline);
	gsl_interp_accel_free(accel);
	for(int i = 0 ; i<data_length; i++){
		delete [] temp[i];
	}
	delete[] temp;
	return pomega;
}

double pdet_triple_detector_fit(double rho_thresh, double rho_opt)
{
	int np = 1e3;
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(np);
	gsl_function F;
	double result, err; 
	double abserr=0, relerr = 1e-5;
	F.function = [](double omega, void * param){ return p_triple_detector_fit(omega);};	

	int errcode = gsl_integration_qag(&F, rho_thresh/rho_opt, 1, abserr, relerr, 
		np, GSL_INTEG_GAUSS15, w, &result, &err);
	gsl_integration_workspace_free(w);
	
	return result;
}
//###########################################################
//Explicit Declarations of the temlates for double and adouble
template double LISA_response_plus_time<double>( double, double , double , double,double , double, double);
template adouble LISA_response_plus_time<adouble>( adouble, adouble , adouble , adouble,adouble , adouble ,adouble);
//
template double LISA_response_cross_time<double>( double, double , double , double,double , double, double);
template adouble LISA_response_cross_time<adouble>( adouble, adouble , adouble , adouble,adouble , adouble, adouble);
//
template double LISA_response_plus<double>( source_parameters<double> *params,double, double , double , double,double , double, double);
template adouble LISA_response_plus<adouble>( source_parameters<adouble> *params,adouble, adouble , adouble , adouble,adouble , adouble ,adouble);
//
template double LISA_response_cross<double>( source_parameters<double> *params,double, double , double , double,double , double, double);
template adouble LISA_response_cross<adouble>( source_parameters<adouble> *params,adouble, adouble , adouble , adouble,adouble , adouble, adouble);
//
template double  right_interferometer_cross<double>(double, double);
template adouble right_interferometer_cross<adouble>(adouble, adouble);
//
template double right_interferometer_plus<double>(double, double);
template adouble right_interferometer_plus<adouble>(adouble, adouble);
//
template double  right_interferometer_x<double>(double, double,double);
template adouble right_interferometer_x<adouble>(adouble, adouble,adouble);
//
template double  right_interferometer_y<double>(double, double,double);
template adouble right_interferometer_y<adouble>(adouble, adouble,adouble);
//
template double  right_interferometer_b<double>(double, double,double);
template adouble right_interferometer_b<adouble>(adouble, adouble,adouble);
//
template double  right_interferometer_l<double>(double, double,double);
template adouble right_interferometer_l<adouble>(adouble, adouble,adouble);
//
template void right_interferometer<double>(det_res_pat<double> *,double, double,double);
template void right_interferometer<adouble>(det_res_pat<adouble> *,adouble, adouble,adouble);
//
template void celestial_horizon_transform<double>(double, double, double, std::string, double *,double*);
template void celestial_horizon_transform<adouble>(adouble, adouble, double, std::string, adouble *,adouble*);
//
template void detector_response_functions_equatorial<double>(double[3][3],double, double, double, double, det_res_pat<double> *);
template void detector_response_functions_equatorial<adouble>(double[3][3],adouble, adouble, adouble, double, det_res_pat<adouble> *);
//
template void detector_response_functions_equatorial<double>(std::string, double, double, double, double, double*,int, double, double, double, double,  det_res_pat<double> *);
template void detector_response_functions_equatorial<adouble>(std::string, adouble, adouble, adouble, double, adouble*,int, adouble, adouble, adouble, adouble, det_res_pat<adouble> *);
//
template void detector_response_functions_equatorial<double>(std::string, double, double, double ,double ,det_res_pat<double> *);
template void detector_response_functions_equatorial<adouble>(std::string, adouble, adouble, adouble ,double ,det_res_pat<adouble> *);
