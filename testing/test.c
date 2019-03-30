#include "waveform_generator_C.h"
#include <stdio.h>
#include <math.h>

void test1(void)
{
	
	int length = 50;
	double mass1 = 200;
	double mass2 = 50;
	char *method= "IMRPhenomD";
	//string method= "ppE_IMRPhenomD_Inspiral";
	//double wpr[length];
	//double wpi[length];
	//double wcr[length];
	//double wci[length];
	double chi1 = .2;
	double chi2 = .1;
	double phic = 2.0;
	double tc = 8.0;
	double Luminosity_Distance = 800.;
	double phi = 0;
	double theta = 0;
	double incl_angle = 0;
	double f_ref = 100;
	double phiRef = 1.0;
	
	double fhigh =100;
	double flow =17;
	double df = (fhigh-flow)/(length-1);
	double *freq = (double *)malloc(sizeof(double) * length);
	double *wpr = (double *)malloc(sizeof(double) * length);
	double *test = (double *)malloc(sizeof(double) * length);
	double *wpi = (double *)malloc(sizeof(double) * length);
	double *wcr = (double *)malloc(sizeof(double) * length);
	double *wci = (double *)malloc(sizeof(double) * length);
	
	for (int i = 0; i < length; i ++)
		test[i] = i;
	for (int i = 0; i < length; i ++)
		wpi[i] = test[i];
	free(test);
	for (int i = 0; i < length; i ++)
		printf("%f, %f\n",test[i],wpi[i]);
	for (int i =  0; i<length; i++)
		freq[i] = (i+1)*df;
fourier_waveformC(freq, //Freqs                                                                                                 
                 length, //length of array                                                                                                       
                 wpr, //waveform plus real                                                                                                    
                 wpi, //waveform plus imag                                                                                                    
                 wcr, //waveform cross real                                                                                                   
                 wci, //waveform cross imag                                                                                                   
                 method, //method of waveform generation                                                                                  
                 mass1, //Mass 1 in solar masses                                                                                                 
                 mass2, //Mass2 in solar masses                                                                                                  
                 Luminosity_Distance, //Distance in Mpc                                                                                                  
                 0, //spin1x                                                                                                                  
                 0, //spin1y                                                                                                                  
                 chi1, //spin1z                                                                                                               
                 0, //spin2x                                                                                                                  
                 0, //spin2y                                                                                                                  
                 chi2, //spin2z                                                                                                               
                 0, //phic                                                                                                                    
                 tc, //tc                                                                                                                     
                 f_ref, //f_ref                                                                                                             
                 phiRef, //phiRef                                                                                                               
                 0, //ppE_beta                                                                                                                
                 0, //ppE_b                                                                                                                   
                 0, //incl_angle                                                                                                              
                 0, //theta                                                                                                                   
                 0 //phi                                                                                                                      
                 ); 
	for (int i = 0; i < length; i ++)
		printf("%f, %f\n",log(abs(wpr[i])),wpi[i]);
	free(freq);
	free(wcr);
	free(wci);
	free(wpr);
	free(wpi);
}
int main(void)
{
	test1();
	return 0;
}
