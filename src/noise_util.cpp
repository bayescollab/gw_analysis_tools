#include "noise_util.h"
#include "util.h"
#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

/*! \file
 * Routines to construct noise curves for various detectors*/

/*! \brief Function to populate the squareroot of the noise curve for various detectors
 *
 * If frequencies are left as NULL, standard frequency spacing is applied and the frequencies are returned, in which case the frequencies argument becomes an output array
 *
 * Detector names must be spelled exactly
 *
 * Detectors include: aLIGO_analytic, Hanford_O1_fitted
 */
void populate_noise(double *frequencies, /**< double array of frquencies (NULL)*/
		std::string detector, /**< String to designate the detector noise curve to be used */
		double *noise_root, /**< ouptput double array for the square root of the PSD of the noise of the specified detector*/
		int length/**< integer length of the output and input arrays*/
		)
{
	if(detector == "aLIGO_analytic")
	{
		if(frequencies == NULL)
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
	if(detector == "Hanford_O1_fitted")
	{
		if(frequencies == NULL)
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
			
		
}


double aLIGO_analytic(double f)
{
	double S= 3e-48;
	double fknee = 70.;
	double x = fknee/f;
	double x4 = x*x*x*x;
	return sqrt( S * (x4 + 2 + 2*x*x)/5 );
}

double Hanford_O1_fitted(double f)
{
	double avec[7]={47.8466,-92.1896,35.9273,-7.61447,0.916742,-0.0588089,0.00156345};
	double S0 = .8464;
	double x = log(f);
	return sqrt(S0) * exp( avec[0] + avec[1]*x + avec[2]*x*x +
            avec[3] * x*x*x + avec[4]*x*x*x*x + avec[5]*x*x*x*x*x + avec[6]*x*x*x*x*x*x);
}

std::complex<double> Q(double theta, double phi, double iota)
{
	double ct = cos(theta);
	double cp2 = cos(2.*phi);
	double sp2 = sin(2.*phi);
	double ci = cos(iota);

	double Fplus = (1./2)*(1+ ct*ct)*cp2;
	double Fcross = ct * sp2;
	std::complex<double> Q = (1+ci*ci)/2. *Fplus + std::complex<double>(0,Fcross*ci);
	return Q;
}


double right_interferometer_plus(double theta, double phi)
{
	double ct = cos(theta);
	double cp2 = cos(2.*phi);

	double Fplus = (1./2)*(1+ ct*ct)*cp2;
	return Fplus;
}

double right_interferometer_cross(double theta, double phi)
{
	double ct = cos(theta);
	double sp2 = sin(2.*phi);
	double Fcross = ct * sp2;
	return Fcross;
}

void celestial_horizon_transform(double RA, /**< in RAD*/
		double DEC, /**< in RAD*/
		double gps_time, 
		std::string detector, 
		double *phi, /**< in RAD*/
		double *theta /**< in RAD*/
		)
{
	double LAT, LONG, azimuth_offset;
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
	celestial_horizon_transform(RA,DEC, gps_time, LONG, LAT, phi, theta);
	*phi += azimuth_offset;
	if(*phi>2*M_PI) *phi -=2*M_PI;
}

/*! \brief calculate Difference in time of arrival (DTOA) for a given source location and 2 different detectors
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
	return (cos(theta1)*R1 - cos(theta2)*R2)/c;
}

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
