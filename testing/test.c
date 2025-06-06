#include "waveform_generator_C.h"
#include <stdio.h>
#include <math.h>

void test1(void)
{
	initiate_LumD_Z_interp_C();	
	int length = 500;
	double mass1 = 20;
	double mass2 = 10;
	//char *method= "ppE_IMRPhenomD_Inspiral";
	char *method= "EdGB_IMRPhenomD_log";
	//char *method= "IMRPhenomD";
	//string method= "ppE_IMRPhenomD_Inspiral";
	//double wpr[length];
	//double wpi[length];
	//double wcr[length];
	//double wci[length];
	double chi1 = .2;
	double chi2 = .1;
	double phic = 2.0;
	double tc = 8.0;
	double Luminosity_Distance = 200.;
	double phi = 0;
	double theta = 0;
	double incl_angle = 0;
	double f_ref = 100;
	double phiRef = 1.0;
	int mods = 1;
	double *betappe = (double *)malloc(sizeof(double)*mods);
	betappe[0] = 67.;
	//betappe[1] = 10;
	//betappe[1] = 10;
	int *bppe = (int *)malloc(sizeof(int)*mods);
	bppe[0] = -1;
	//bppe[1] = -2;
	//bppe[1] = -3;
	
	double fhigh =550;
	double flow =17;
	double df = (fhigh-flow)/(length-1);
	double *freq = (double *)malloc(sizeof(double) * length);
	double *wpr = (double *)malloc(sizeof(double) * length);
	double *wpi = (double *)malloc(sizeof(double) * length);
	double *wcr = (double *)malloc(sizeof(double) * length);
	double *wci = (double *)malloc(sizeof(double) * length);
	double *phase = (double *)malloc(sizeof(double) * length);
	double *amplitude = (double *)malloc(sizeof(double) * length);
	
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
                 betappe, //ppE_beta                                                                                                                
                 bppe, //ppE_b                                                                                                                   
		 mods, //Nmod
                 0, //incl_angle                                                                                                              
                 0, //theta                                                                                                                   
                 0 //phi                                                                                                                      
                 ); 
	//betappe[0] = 20;
fourier_amplitudeC(freq, //Freqs                                                                                                 
                 length, //length of array                                                                                                       
                 amplitude, //waveform plus real                                                                                                    
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
                 0, //incl_angle                                                                                                              
                 0, //theta                                                                                                                   
                 0 //phi                                                                                                                      
                 ); 
	//betappe[0] = 20;
fourier_phaseC(freq, //Freqs                                                                                                 
                 length, //length of array                                                                                                       
                 phase, //waveform plus real                                                                                                    
                 method, //method of waveform generation                                                                                  
                 mass1, //Mass 1 in solar masses                                                                                                 
                 mass2, //Mass2 in solar masses                                                                                                  
		Luminosity_Distance,
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
                 betappe, //ppE_beta                                                                                                                
                 bppe, //ppE_b                                                                                                                   
		 mods, //Nmod
                 0, //incl_angle                                                                                                              
                 0, //theta                                                                                                                   
                 0 //phi                                                                                                                      
                 ); 
	//for (int i = 0; i < length; i ++){
	//	//printf("%.10e, %.10e\n",log(abs(wpr[i])),log(abs(wpi[i])));
	//	printf("%.10e, %.10e\n",wpr[i],wpi[i]);
	//}
	for (int i = 0; i < length; i ++){
		printf("%.10e\n",phase[i]);
	}
	//for (int i = 0; i < length; i ++){
	//	printf("%.10e\n",amplitude[i]);
	//}
	free(freq);
	free(wcr);
	free(wci);
	free(wpr);
	free(wpi);
	free(phase);
	free(amplitude);
	free(betappe);
	free(bppe);
	free_LumD_Z_interp_C();
}
int main(void)
{
	test1();
	return 0;
}
