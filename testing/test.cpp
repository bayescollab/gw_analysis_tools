#include <iostream>
#include <fstream>
#include <time.h>
#include <complex>
#include <string>
#include "waveform_generator.h"
#include "IMRPhenomD.h"
#include "mcmc_routines.h"
#include "noise_util.h"
#include "util.h"
#include "waveform_util.h"
#include <adolc/adouble.h>
#include "fisher.h"
#include "ppE_IMRPhenomD.h"
#include "IMRPhenomP.h"
#include "waveform_generator_C.h"
#include "mcmc_sampler.h"
#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_rng.h"
#include "adolc/adouble.h"
#include "adolc/drivers/drivers.h"
#include "adolc/taping.h"
#include "limits"


using namespace std;

void test1();
void test2();
void test3();
void test4();
void test5();
void test6();
void test7();
void test8();
void test9();
void test10();
double test_ll(double *pos, int dim);
double test_lp(double *pos, int dim);
double test_lp_GW(double *pos, int dim);
void test_fisher(double *pos, int dim, double **fisher);
double log_student_t (double *x,int dim);
double log_neil_proj3 (double *x,int dim);
void fisher_neil_proj3 (double *x,int dim, double **fish);
adouble dist(adouble *pos, int dimension);

const gsl_rng_type* Y;
gsl_rng * g;


int main(){

	//gsl_rng_env_setup();
	//Y = gsl_rng_default;
	//g = gsl_rng_alloc(Y);
	test9();	
	return 0;
}
void test10()
{
	//char *ptr = (char *)malloc(sizeof(char) * 10);
	//ptr = "PhenomD";
	//std::cout<<ptr<<std::endl;
	//
	//char *ptr2;
	//ptr2 = ptr;
	//
	//std::cout<<ptr2<<std::endl;
	//free(ptr2);
	//free (ptr);
	
	std::string detectors[2];
	detectors[0] = "Hanford";
	detectors[1] = "Livingston";
	std::string *ptr1 = new std::string[2];
	ptr1[0] = detectors[0];
	ptr1[1] = detectors[1];
	std::cout<<ptr1[0]<<std::endl;
	std::cout<<ptr1[1]<<std::endl;
	std::string *new_detect = ptr1;
	//new_detect = detectors;
	std::cout<<new_detect[0]<<std::endl;
	std::cout<<new_detect[1]<<std::endl;
	delete [] ptr1;

	//int num =2;
	//int namelen = 50;
	//char **ptr = (char **)malloc(sizeof(char*) * num);
	//for (int i = 0; i<num; i++)
	//	ptr[i] = (char *)malloc(sizeof(char) * namelen);
	//ptr[0] = "PhenomD";
	//ptr[1] = "PhenomD";
	//
	//std::cout<<ptr[0]<<std::endl;
	//std::cout<<ptr[1]<<std::endl;
	//char **ptr2;
	//ptr2 = ptr;
	//
	//std::cout<<ptr2[0]<<std::endl;
	//std::cout<<ptr2[1]<<std::endl;
	//free(ptr2);
	//for (int i = 0; i<num; i++)
	//	free(ptr[i]);
	//free (ptr);
	
}
void test9()
{
	int raw_length = 28673;
	int cutoff = 600;
	//int high_cut = 11000;
	//int length = high_cut-cutoff;
	int length = raw_length-cutoff;
	int num_detectors =1;

	double **temp_data = allocate_2D_array(raw_length,2);
	double *temp_psd = (double *)malloc(sizeof(double)*raw_length);
	double *temp_freq = (double *)malloc(sizeof(double)*raw_length);
	std::string filebase = "testing/data/gw_150914_";
	read_file(filebase+"data.csv",temp_data, raw_length,2);
	read_file(filebase+"psd.csv",temp_psd);
	read_file(filebase+"freq.csv",temp_freq);

	std::complex<double> **data= (std::complex<double>**)malloc(
			sizeof(std::complex<double>*)*num_detectors);
	double **psd = (double **)malloc(sizeof(double *)*num_detectors);
	double **frequencies = (double **)malloc(sizeof(double *)*num_detectors);
	int *data_length= (int*)malloc(sizeof(int)*num_detectors);
	data_length[0] =length;

	//bool check = true;
	for (int i =0; i<num_detectors; i++){
		data[i] = (std::complex<double> *)malloc(
			sizeof(std::complex<double>)*data_length[0]);
		
		psd[i] = (double *)malloc(sizeof(double)*data_length[0]);
		frequencies[i] = (double *)malloc(sizeof(double)*data_length[0]);
		for(int j = 0; j<data_length[0]; j++){
			frequencies[i][j] = temp_freq[j+cutoff];	
			psd[i][j] = (temp_psd[j+cutoff]);	
			data[i][j] = std::complex<double>(temp_data[j+cutoff][0],temp_data[j+cutoff][1]);	
			//std::cout<<frequencies[i][j]<<std::endl;
			//std::cout<<psd[i][j]<<std::endl;
			//std::cout<<data[i][j]<<std::endl;
			//if(temp_freq[j]>400 && check){std::cout<<j<<std::endl;check=false;}
		}
	}

	deallocate_2D_array(temp_data,raw_length,2);
	free(temp_psd);
	free(temp_freq);
	//#########################################################
	//MCMC options
	int dimension = 5;
	double initial_pos[dimension]={log(400*MPC_SEC),log(30*MSOL_SEC), .24, 0,0};
	//double initial_pos[dimension]={log(200*MPC_SEC),log(20*MSOL_SEC), .15, 0,0};
	int N_steps = 5000;
	int chain_N= 20;
	double ***output;
	output = allocate_3D_array( chain_N, N_steps, dimension );
	//double *initial_pos_ptr = initial_pos;
	int swp_freq = 100;
	//double chain_temps[chain_N] ={1,2,3,10,12};
	double chain_temps[chain_N];
	double temp_step = 400./(chain_N); 
	for(int i =0; i < chain_N;  i ++)
		chain_temps[i] = 1+ temp_step * i;
	//double chain_temps[chain_N] ={1};
	
	//#########################################################
	//GW options
	std::string *detectors = new std::string[1];//(std::string*)malloc(sizeof(std::string)*50*num_detectors);
	detectors[0] = "Hanford";
	std::string generation_method = "IMRPhenomD";
	
	
	std::string autocorrfile = "testing/data/auto_corr_mcmc.csv";
	std::string chainfile = "testing/data/mcmc_output.csv";
	std::string statfilename = "testing/data/mcmc_statistics.txt";

	MCMC_MH_GW(output, dimension, N_steps, chain_N, initial_pos,chain_temps, 
			swp_freq, test_lp_GW, num_detectors, data, psd, 
			frequencies, data_length, detectors, generation_method,
			statfilename,"",autocorrfile);	
	std::cout<<"ENDED"<<std::endl;

	double **output_transform=(double **)malloc(sizeof(double*)*N_steps);
	for (int j =0; j<N_steps; j++)
		output_transform[j] = (double *)malloc(sizeof(double)*dimension);

	for(int j = 0; j<N_steps;j++){
			output_transform[j][0]=std::exp(output[0][j][0])/MPC_SEC;
			output_transform[j][1]=std::exp(output[0][j][1])/MSOL_SEC;
			output_transform[j][2]=output[0][j][2];
			output_transform[j][3]=output[0][j][3];
			output_transform[j][4]=output[0][j][4];
	}
	write_file(chainfile, output_transform, N_steps, dimension);
	//ofstream mcmc_out;
	//mcmc_out.open("testing/data/mcmc_output.csv");
	//mcmc_out.precision(15);
	////for(int i = 0;i<chain_N;i++){
	//for(int j = 0; j<N_steps;j++){
	//	//for(int k = 0; k<dimension; k++){
	//		mcmc_out<<std::exp(output[0][j][0])/MPC_SEC<<" , "<<std::exp(output[0][j][1])/MSOL_SEC<<" , "<<output[0][j][2]<<" , "<<output[0][j][3]<<" , "<<output[0][j][4]<<endl;
	//	//}
	//}
	////}
	//mcmc_out.close();

	deallocate_3D_array(output, chain_N, N_steps, dimension);
	for(int i =0; i< num_detectors; i++){
		free(data[i]);
		free(psd[i]);
		free(frequencies[i]);
	}
	for(int i =0; i< N_steps; i++){
		free(output_transform[i]);
	}
	free(output_transform);
	free(data);
	free(psd);
	free(frequencies );
	delete [] detectors;
	//free(detectors);
	free(data_length);
}
void test8()
{
	//#########################################################
	//Make trial data
	gen_params params;
	IMRPhenomD<double> modeld;
	int length = 200;
	double chirpm = 30.78;
	double eta =.21;
	params.mass1 = calculate_mass1(chirpm,eta);
	params.mass2 = calculate_mass2(chirpm,eta);
	string method= "IMRPhenomD";
	//string method= "ppE_IMRPhenomD_Inspiral";
	complex<double> waveformout[length];
	params.spin1[0] = 0;
	params.spin1[1] = 0;
	params.spin1[2] = -.0;
	params.spin2[0] = 0;
	params.spin2[1] = 0;
	params.spin2[2] = .0;
	double *spin1  = params.spin1;
	double *spin2= params.spin2;
	params.phic = .0;
	params.tc = -.0;
	params.Luminosity_Distance = 410.;
	//params.betappe = new double[1] ;
	//params.betappe[0]=-50;
	//params.bppe  =new int[1];
	//params.bppe[0] = -1;
	//params.Nmod = 1;
	params.NSflag = false;
	params.phi = 0;
	params.theta = 0;
	params.incl_angle = 0;
	params.sky_average=false;
	//params.f_ref = 30.5011;
	//params.phiRef =58.944425/2.;
	
	double fhigh =300;
	double flow =10;
	double df = (fhigh-flow)/(length-1);
	double *freq = (double *)malloc(sizeof(double) * length);

	cout<<"Freq spacing "<<df<<endl;

	for(int i=0;i<length;i++)
		freq[i]=flow+i*df;

	clock_t  start, end;
	start = clock(); 
	fourier_waveform(freq, length, waveformout,method,&params);
	end=clock();

	double noise[length];
	populate_noise(freq,"Hanford_O1_fitted", noise,length);
	for (int i =0; i<length;i++){
		noise[i] = noise[i]*noise[i];
		//std::cout<<noise[i]<<std::endl;
	}
	double snr = calculate_snr("Hanford_O1_fitted",waveformout, freq, length);
	std::cout<<"SNR of injection: "<<snr<<std::endl;
	//#########################################################
	


	//#########################################################
	//MCMC options
	int dimension = 5;
	double initial_pos[dimension]={log(params.Luminosity_Distance*MPC_SEC),log(chirpm*MSOL_SEC), eta, params.spin1[2],params.spin2[2]};
	//double initial_pos[dimension]={log(200*MPC_SEC),log(20*MSOL_SEC), .15, 0,0};
	int N_steps = 100;
	int chain_N= 10;
	double ***output;
	output = allocate_3D_array( chain_N, N_steps, dimension );
	//double *initial_pos_ptr = initial_pos;
	int swp_freq = 10;
	//double chain_temps[chain_N] ={1,2,3,10,12};
	double chain_temps[chain_N];
	double temp_step = 400./(chain_N);
	for(int i =0; i < chain_N;  i ++)
		//chain_temps[i]=1.;
		chain_temps[i] = 1+ temp_step * i;
	//double chain_temps[chain_N] ={1};
	
	//#########################################################
	//GW options
	int num_detectors =1;
	int *data_length= (int*)malloc(sizeof(int)*num_detectors);
	data_length[0] =length;
	std::complex<double> **data= (std::complex<double>**)malloc(
			sizeof(std::complex<double>*)*num_detectors);
	double **psd = (double **)malloc(sizeof(double *)*num_detectors);
	double **frequencies = (double **)malloc(sizeof(double *)*num_detectors);
	std::string *detectors = new std::string[1];//(std::string*)malloc(sizeof(std::string)*10*num_detectors);
	detectors[0] = "Hanford";
	std::string generation_method = "IMRPhenomD";
	
	for (int i =0; i<num_detectors; i++){
		data[i] = (std::complex<double> *)malloc(
			sizeof(std::complex<double>)*data_length[i]);
		
		psd[i] = (double *)malloc(sizeof(double)*data_length[i]);
		frequencies[i] = (double *)malloc(sizeof(double)*data_length[i]);
		for(int j = 0; j<data_length[i]; j++){
			frequencies[i][j] = freq[j];	
			psd[i][j] = noise[j];	
			data[i][j] = waveformout[j];	
		}
	}
	
	std::string autocorrfile = "testing/data/auto_corr_mcmc.csv";
	std::string chainfile = "testing/data/mcmc_output.csv";
	std::string statfilename = "testing/data/mcmc_statistics.txt";

	
	std::cout<<"TEST"<<std::endl;
	//fftw_outline plan;
	//initiate_likelihood_function(&plan, data_length[0]);
	//maximized_coal_Log_Likelihood(data[0], psd[0], frequencies[0], data_length[0],
	//		&params, "Hanford", "IMRPhenomD", &plan);
	//deactivate_likelihood_function(&plan);
	MCMC_MH_GW(output, dimension, N_steps, chain_N, initial_pos,chain_temps, 
			swp_freq, test_lp_GW, num_detectors, data, psd, 
			frequencies, data_length, detectors, generation_method,
			statfilename, "" ,autocorrfile);	
	std::cout<<"ENDED"<<std::endl;

	double **output_transform=(double **)malloc(sizeof(double*)*N_steps);
	for (int j =0; j<N_steps; j++)
		output_transform[j] = (double *)malloc(sizeof(double)*dimension);

	for(int j = 0; j<N_steps;j++){
			output_transform[j][0]=std::exp(output[0][j][0])/MPC_SEC;
			output_transform[j][1]=std::exp(output[0][j][1])/MSOL_SEC;
			output_transform[j][2]=output[0][j][2];
			output_transform[j][3]=output[0][j][3];
			output_transform[j][4]=output[0][j][4];
	}
	write_file(chainfile, output_transform, N_steps, dimension);
	
	//ofstream mcmc_out;
	//mcmc_out.open("testing/data/mcmc_output.csv");
	//mcmc_out.precision(15);
	////for(int i = 0;i<chain_N;i++){
	//for(int j = 0; j<N_steps;j++){
	//	//for(int k = 0; k<dimension; k++){
	//		mcmc_out<<std::exp(output[0][j][0])/MPC_SEC<<" , "<<std::exp(output[0][j][1])/MSOL_SEC<<" , "<<output[0][j][2]<<" , "<<output[0][j][3]<<" , "<<output[0][j][4]<<endl;
	//	//}
	//}
	////}
	//mcmc_out.close();

	deallocate_3D_array(output, chain_N, N_steps, dimension);
	for(int i =0; i< num_detectors; i++){
		free(data[i]);
		free(psd[i]);
		free(frequencies[i]);
	}
	for(int i =0; i< N_steps; i++){
		free(output_transform[i]);
	}
	free(output_transform);
	free(data);
	free(psd);
	free(freq);
	free(frequencies );
	//free(detectors);
	delete [] detectors;
	free(data_length);
}
void test7()
{
	int dimension = 2;
	double initial_pos[2]={6,5.};

	
	int N_steps = 50;
	int chain_N= 1;
	double ***output;
	output = allocate_3D_array( chain_N, N_steps, dimension );
	//double *initial_pos_ptr = initial_pos;
	int swp_freq = 30;
	//double chain_temps[chain_N] ={1,2,3,10,12};
	double chain_temps[chain_N];
	double temp_step = 500./(chain_N);
	for(int i =0; i < chain_N;  i ++)
		//chain_temps[i]=1.;
		chain_temps[i] = 1+ temp_step * i;
	//double chain_temps[chain_N] ={1};
	std::string autocorrfile = "testing/data/auto_corr_mcmc.csv";
	std::string chainfile = "testing/data/mcmc_output.csv";
	std::string statfilename = "testing/data/mcmc_statistics.txt";
	
	//MCMC_MH(output, dimension, N_steps, chain_N, initial_pos,chain_temps, swp_freq, test_lp, log_neil_proj3,fisher_neil_proj3 );	
	MCMC_MH(output, dimension, N_steps, chain_N, initial_pos,chain_temps, swp_freq, test_lp, log_neil_proj3,NULL,statfilename,chainfile,autocorrfile );	
	std::cout<<"ENDED"<<std::endl;

	//write_file("testing/data/mcmc_output.csv", output[0],N_steps, dimension);
	//ofstream mcmc_out;
	//mcmc_out.open("testing/data/mcmc_output.csv");
	//mcmc_out.precision(15);
	////for(int i = 0;i<chain_N;i++){
	//for(int j = 0; j<N_steps;j++){
	//	//for(int k = 0; k<dimension; k++){
	//		mcmc_out<<output[0][j][0]<<" , "<<output[0][j][1]<<endl;
	//	//}
	//}
	////}
	//mcmc_out.close();

	deallocate_3D_array(output, chain_N, N_steps, dimension);
}

void test6()
{

	gen_params params;
	IMRPhenomD<double> modeld;
	IMRPhenomD<adouble> modela;
	int length = 8;
	double chirpm = 49.78;
	double eta =.21;
	params.mass1 = calculate_mass1(chirpm,eta);
	params.mass2 = calculate_mass2(chirpm,eta);
	string method= "IMRPhenomD";
	//string method= "ppE_IMRPhenomD_Inspiral";
	double amp[length];
	double phaseout[length];
	complex<double> waveformout[length];
	params.spin1[0] = 0;
	params.spin1[1] = 0;
	params.spin1[2] = -.2;
	params.spin2[0] = 0;
	params.spin2[1] = 0;
	params.spin2[2] = .4;
	double *spin1  = params.spin1;
	double *spin2= params.spin2;
	params.phic = .0;
	params.tc = -.0;
	params.Luminosity_Distance = 410.;
	//params.betappe = new double[1] ;
	//params.betappe[0]=-100.;
	//params.bppe  =new int[1];
	//params.bppe[0] = -3;
	//params.Nmod = 1;
	params.NSflag = false;
	params.phi = 0;
	params.theta = 0;
	params.incl_angle = 0;
	params.sky_average=true;
	//params.f_ref = 100;
	//params.phiRef = 1.0;
	
	double fhigh =300;
	double flow =15;
	double df = (fhigh-flow)/(length-1);
	double *freq = (double *)malloc(sizeof(double) * length);

	cout<<"Freq spacing "<<df<<endl;

	for(int i=0;i<length;i++)
		freq[i]=flow+i*df;

	int dimension = 7;
	int dimensionmcmc = 5;

	clock_t start7,end7;

	double **output = (double **)malloc(dimension * sizeof(**output));	
	double **output2 = (double **)malloc(dimension * sizeof(**output2));	
	double **output3 = (double **)malloc(dimensionmcmc * sizeof(**output3));	

	for (int i = 0;i<dimension;i++){
		output[i] = (double *)malloc(dimension*sizeof(double));
		output2[i] = (double *)malloc(dimension*sizeof(double));
	}
	for (int i = 0;i<dimensionmcmc;i++){
	
		output3[i] = (double*)malloc(dimensionmcmc*sizeof(double));
	}
	start7 = clock();
	fisher(freq, length, "MCMC_IMRPhenomD_single_detect","Hanford_O1_fitted", 
			output3, dimensionmcmc, &params );

	end7 = clock();

	cout<<"TIMING: FISHER: "<<(double)(end7-start7)/CLOCKS_PER_SEC<<endl;

	cout.precision(5);
	for (int i = 0;i <dimensionmcmc;i++)
	{
		for (int j=0;j <dimensionmcmc; j++)
			cout<<output3[i][j]<<"   ";
		cout<<endl;
	}
	
	start7 = clock();
	fisher(freq, length, "IMRPhenomD","Hanford_O1_fitted", output, dimension, 
				&params );

	end7 = clock();
	cout<<"TIMING: FISHER: "<<(double)(end7-start7)/CLOCKS_PER_SEC<<endl;
	cout.precision(5);
	for (int i = 0;i <dimension;i++)
	{
		for (int j=0;j <dimension; j++)
			cout<<output[i][j]<<"   ";
		cout<<endl;
	}
	start7 = clock();
	fisher_autodiff(freq, length, "IMRPhenomD","Hanford_O1_fitted", output2, dimension, 
				&params );

	end7 = clock();
	cout<<"TIMING: FISHER autodiff: "<<(double)(end7-start7)/CLOCKS_PER_SEC<<endl;
	for (int i = 0;i <dimension;i++)
	{
		for (int j=0;j <dimension; j++)
			cout<<output2[i][j]<<"   ";
		cout<<endl;
	}
	std::cout<<"fractional DIFF: "<<std::endl;
	for (int i = 0;i <dimension;i++)
	{
		for (int j=0;j <dimension; j++)
			cout<<(output2[i][j]-output[i][j])/output2[i][j]<<"   ";
		cout<<endl;
	}
	for(int i =0; i<dimension; i++)
	{
		free(output[i]);
		free(output2[i]);
	}
	for(int i =0; i<dimensionmcmc; i++)
	{
		free(output3[i]);
	}
	free(output);
	free(output2);
	free(output3);
	free(freq);
}
void test5()
{

	gen_params params;
	IMRPhenomD<double> modeld;
	IMRPhenomD<adouble> modela;
	int length = 1000;
	params.mass1 = 200;
	params.mass2 = 50;
	params.spin1[0] = 0;
	params.spin1[1] = 0;
	params.spin1[2] = -.2;
	params.spin2[0] = 0;
	params.spin2[1] = 0;
	params.spin2[2] = .9;
	double *spin1  = params.spin1;
	double *spin2= params.spin2;
	params.phic = 2.0;
	params.tc = 8.0;
	params.Luminosity_Distance = 800.;
	params.betappe = new double[1] ;
	params.betappe[0]=10.;
	params.bppe  =new int[1];
	params.bppe[0] = -1;
	params.Nmod = 1;
	params.NSflag = false;
	params.phi = 0;
	params.theta = 0;
	params.incl_angle = 0;
	params.f_ref = 100;
	params.phiRef = 1.0;
	params.sky_average=true;
	
	double fhigh =100;
	double flow =17;
	double df = (fhigh-flow)/(length-1);
	double *freq = (double *)malloc(sizeof(double) * length);
	for(int i=0;i<length;i++)
		freq[i]=flow+i*df;

	double noise[length];
	populate_noise(freq,"Hanford_O1_fitted", noise,length);
	for (int i =0; i<length;i++)
		noise[i] = noise[i]*noise[i];

	int dimension =8;
	clock_t start7,end7;
	double **output = (double **)malloc(dimension * sizeof(**output));	
	for (int i = 0;i<dimension;i++)
		output[i] = (double *)malloc(dimension*sizeof(double));
	
	start7 = clock();
	fisher(freq, length, "ppE_IMRPhenomD_Inspiral","Hanford_O1_fitted", output, dimension, 
				&params );

	end7 = clock();
	cout<<"TIMING: FISHER: "<<(double)(end7-start7)/CLOCKS_PER_SEC<<endl;
	for (int i = 0;i <dimension;i++)
	{
		for (int j=0;j <dimension; j++)
			cout<<output[i][j]<<"   ";
		cout<<endl;
	}
	free(output);
	free(freq);
}
void test4()
{

	cout.precision(15);

	gen_params params;

	int length = 16000;
	double fhigh =20;
	double flow =17;
	double df = (fhigh-flow)/(length-1);
	double *freq = (double *)malloc(sizeof(double) * length);
	double *freqnew = (double *)malloc(sizeof(double) * length);

	cout<<"Freq spacing "<<df<<endl;

	for(int i=0;i<length;i++){
		freq[i]=flow+i*df;
		//freqnew[i] = freq[i]-.1;
		freqnew[i] = freq[i];
	}
	//for(int i=0;i<length;i++){
	//	cout<<freqnew[i]<<" "<<freq[i] <<endl;
	//}
	
	double chirpmass = 20;
	double eta = .2;
	params.mass1 = calculate_mass1(chirpmass,eta);
	params.mass2 = calculate_mass2(chirpmass,eta);
	//string method= "IMRPhenomD";
	//string method= "ppE_IMRPhenomD_Inspiral";
	string method= "IMRPhenomPv2";
	complex<double> *waveformout = (complex<double> *)malloc(sizeof(complex<double>) * length);
	complex<double> *waveformoutcross = (complex<double> *)malloc(sizeof(complex<double>) * length);
	params.spin1[0] = 0.01;
	params.spin1[1] = 0;
	params.spin1[2] = .1;
	params.spin2[0] = 0.01;
	params.spin2[1] = 0;
	params.spin2[2] = .3;
	double *spin1  = params.spin1;
	double *spin2= params.spin2;
	params.phic = 1.0;
	params.tc = .0;
	params.Luminosity_Distance = 100.2;
	params.betappe = new double[1] ;
	params.betappe[0]=1.;
	params.bppe  =new int[1];
	params.bppe[0] = -1;
	params.Nmod = 1;
	params.NSflag = false;
	params.phi = 1.2;
	params.theta =3.4;
	params.incl_angle = 1.3;
	//params.f_ref=10;
	//params.phiRef=0.;
	
	
	clock_t  start, end;
	start = clock(); 
	fourier_waveform(freq, length, waveformout,waveformoutcross,method,&params);
	end=clock();
	cout<<"TIMING waveform: "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;

	double *noise = (double *)malloc(sizeof(double)*length);
	populate_noise(freq,"Hanford_O1_fitted", noise,length);
	for (int i =0; i<length;i++)
		noise[i] = noise[i]*noise[i];
	
	//params.mass1 = 100;
	//params.mass2 = 5;
	fftw_outline plan;
	initiate_likelihood_function(&plan,length);
	
	int masslen = 10;
	//double chirp = calculate_chirpmass(params.mass1,params.mass2);
	//double eta = calculate_eta(params.mass1,params.mass2);
	double masses[masslen];
	for (int i =0; i <masslen;i++)
		masses[i] = (i+1.)*.1*chirpmass;
	double mass2 = params.mass2;
	cout<<"Mass2: "<<mass2<<std::endl;
	


	double *real = (double *)malloc(sizeof(double)*length);
	double *imag = (double *)malloc(sizeof(double)*length);
	for ( int i =0; i<length;i++)
	{
		real[i]=(waveformout[i]).real();
		imag[i]=(waveformout[i]).imag();
	}
	
	complex<double> *hplus_new = (complex<double> *)malloc(sizeof(complex<double>)*length);
	complex<double> *detector_response = (complex<double> *)malloc(sizeof(complex<double>)*length);
	complex<double> *hcross_new = (complex<double> *)malloc(sizeof(complex<double>)*length);
	complex<double> *hplus_old = (complex<double> *)malloc(sizeof(complex<double>)*length);
	double llnew, llold,sum;
	complex<double> q;
	for (int i =0;i<masslen;i++)
	{
		params.mass1 = calculate_mass1(masses[i],eta);
		params.mass2 = calculate_mass2(masses[i],eta);

		
		start = clock();
		llnew = maximized_coal_Log_Likelihood(waveformout, noise,freqnew,length, 
					&params,"Hanford","IMRPhenomPv2",&plan);
		//llnew = maximized_coal_Log_Likelihood(real,imag, noise,freqnew,length, 
					//&params,"Hanford","IMRPhenomD",&plan);
		end = clock();
		cout<<"logl new  TIMING: "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;
		start = clock();
		fourier_detector_response(freqnew,length, detector_response,"Hanford",
					"IMRPhenomPv2",&params);
		end = clock();
		cout<<"waveform Pv2  TIMING: "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;
		
		start = clock();
		llold = maximized_coal_log_likelihood_IMRPhenomD_Full_Param(freq, length, real,imag, 
					noise,  calculate_chirpmass(params.mass1,params.mass2), 
					calculate_eta(params.mass1,params.mass2), params.spin1[2], params.spin2[2],params.Luminosity_Distance,params.theta,params.phi,params.incl_angle, false,&plan);
		end = clock();
		
		start = clock();
		q = Q(params.theta,params.phi,params.incl_angle);
		fourier_waveform(freq,length, hplus_old,
					"IMRPhenomD",&params);
		for (int i = 0; i<length;i++)
			hplus_old[i] = q * hplus_old[i];
		end = clock();
		cout<<"waveform old  TIMING: "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;
		cout<<"logl old  TIMING: "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;
		
		//cout<<"LOGLnew: "<<llnew<<endl;
		//cout<<"LOGLold: "<<llold<<endl;
		cout<<"chirpmass: "<<masses[i]<<" new: "<<llnew<<" old: "<<llold<<" diff ll: "<<(llold-llnew)/llold<<endl;
		for (int i =0 ; i< length;i ++)
			sum += abs((detector_response[i]-hplus_old[i])/hplus_old[i]);
		cout<<"Average Diff waveform plus: "<<sum/length<<endl;
	}



	free(hplus_new);
	free(hcross_new);
	free(hplus_old);
	free(real);
	free(imag);
	free(waveformout);
	free(waveformoutcross);
	free(freq);
	free(freqnew);
	free(detector_response);
	free(noise);
	delete [] params.betappe;
	delete [] params.bppe;

	deactivate_likelihood_function(&plan);	
}
void test3()
{
	gen_params params;
	IMRPhenomPv2<double> modeld;
	IMRPhenomPv2<adouble> modela;
	int length = 900;
	params.mass1 = 200;
	params.mass2 = 50;
	string method= "IMRPhenomPv2";
	//string method= "IMRPhenomD";
	//string method= "ppE_IMRPhenomD_Inspiral";
	double amp[length];
	double phaseout[length];
	complex<double> waveformout_plus[length];
	complex<double> waveformout_cross[length];
	params.spin1[0] = .0;
	params.spin1[1] = .0;
	params.spin1[2] = -.2;
	params.spin2[0] = .0;
	params.spin2[1] = 0.0;
	params.spin2[2] = .9;
	double *spin1  = params.spin1;
	double *spin2= params.spin2;
	params.phic = 2.0;
	params.tc = 8.0;
	params.Luminosity_Distance = 800.;
	params.betappe = new double[1] ;
	params.betappe[0]=1.;
	params.bppe  =new int[1];
	params.bppe[0] = -1;
	params.Nmod = 1;
	params.NSflag = false;
	//params.phi = M_PI/3.;
	//params.theta = M_PI/3;
	params.phi = 0;
	params.theta = 0;
	params.incl_angle = 0;
	
	double freq[length];
	for(int i=0;i<length;i++)
		freq[i]=10.+i*1e-1;

	clock_t  start, end;
	start = clock(); 
	fourier_waveform(freq, length, waveformout_plus,waveformout_cross,method,&params);
	end=clock();
	cout<<"TIMING waveform: "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;
	delete [] params.betappe;
	delete [] params.bppe;
	
}
void test2()
{
	double alpha, epsilon;
	IMRPhenomPv2<double> modelP;
	alpha = modelP.alpha(2,3,1./2,.75);
	epsilon = modelP.epsilon(2,3,1./2,.75);
	cout<<alpha<<endl; 
	cout<<epsilon<<endl; 
	long factorial_num = factorial(15);
	cout<<factorial_num<<endl;

	double d = modelP.d(2,1,0,.4);
	cout<<d<<endl;
}
void test1()
{

	initiate_LumD_Z_interp();
	gen_params params;
	IMRPhenomD<double> modeld;
	IMRPhenomD<adouble> modela;
	int length = 5000;
	double chirpm = 49.78;
	double eta =.21;
	params.mass1 = calculate_mass1(chirpm,eta);
	params.mass2 = calculate_mass2(chirpm,eta);
	string method= "IMRPhenomD";
	//string method= "ppE_IMRPhenomD_Inspiral";
	double amp[length];
	double phaseout[length];
	complex<double> waveformout[length];
	params.spin1[0] = 0;
	params.spin1[1] = 0;
	params.spin1[2] = -.2;
	params.spin2[0] = 0;
	params.spin2[1] = 0;
	params.spin2[2] = .4;
	double *spin1  = params.spin1;
	double *spin2= params.spin2;
	params.phic = .0;
	params.tc = -.0;
	params.Luminosity_Distance = 410.;
	params.betappe = new double[1] ;
	params.betappe[0]=-50;
	params.bppe  =new int[1];
	params.bppe[0] = -1;
	params.Nmod = 1;
	params.NSflag = false;
	params.phi = 0;
	params.theta = 0;
	params.incl_angle = 0;
	params.sky_average=true;
	//params.f_ref = 30.5011;
	//params.phiRef =58.944425/2.;
	
	double fhigh =200;
	double flow =10;
	double df = (fhigh-flow)/(length-1);
	double *freq = (double *)malloc(sizeof(double) * length);

	cout<<"Freq spacing "<<df<<endl;

	for(int i=0;i<length;i++)
		freq[i]=flow+i*df;

	clock_t  start, end;
	start = clock(); 
	fourier_waveform(freq, length, waveformout,method,&params);
	end=clock();
	
	clock_t  start2, end2;
	start2 = clock(); 
	fourier_amplitude(freq, length, amp,method,&params);
	end2=clock();

	clock_t  start3, end3;
	start3 = clock(); 
	fourier_phase(freq, length, phaseout,method,&params);
	end3=clock();

	ofstream ampfile;
	ampfile.open("testing/data/amplitude_output.csv");
	ampfile.precision(15);
	for(int i = 0;i<length;i++)
		ampfile<<freq[i]<<','<<amp[i]<<endl;
	ampfile.close();
	
	ofstream phasefile;
	phasefile.open("testing/data/phase_output.csv");
	phasefile.precision(15);
	for(int i = 0;i<length;i++)
		phasefile<<freq[i]<<','<<phaseout[i]<<endl;
	phasefile.close();

	ofstream wavefilereal;
	wavefilereal.open("testing/data/real_waveform_output.csv");
	wavefilereal.precision(15);
	for(int i = 0;i<length;i++)
		wavefilereal<<freq[i]<<','<<real(waveformout[i])<<endl;
	wavefilereal.close();

	ofstream wavefileimag;
	wavefileimag.open("testing/data/imag_waveform_output.csv");
	wavefileimag.precision(15);
	for(int i = 0;i<length;i++)
		wavefileimag<<freq[i]<<','<<imag(waveformout[i])<<endl;
	wavefileimag.close();
	cout<<"TIMING waveform: "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;
	cout<<"TIMING amp: "<<(double)(end2-start2)/CLOCKS_PER_SEC<<endl;
	cout<<"TIMING phase: "<<(double)(end3-start3)/CLOCKS_PER_SEC<<endl;
	
	double noise[length];
	populate_noise(freq,"Hanford_O1_fitted", noise,length);
	for (int i =0; i<length;i++)
		noise[i] = noise[i]*noise[i];

	double snr = calculate_snr("Hanford_O1_fitted", waveformout, freq, length);
	cout<<"SNR: "<<snr<<endl;
	snr = data_snr_maximized_extrinsic(freq,length,waveformout,"Hanford_O1_fitted","IMRPhenomD",
			params);
	cout<<" FULL SNR: "<<snr<<endl;

	double logl;
	fftw_outline plan;
	clock_t  start4, end4;
	initiate_likelihood_function(&plan,length);
	start4 = clock();
	logl = maximized_coal_log_likelihood_IMRPhenomD(freq, length, waveformout, 
				noise, snr, calculate_chirpmass(params.mass1,params.mass2), 
				calculate_eta(params.mass1,params.mass2), spin1[2], spin2[2], false,&plan);
	end4 = clock();
	double logl2 = maximized_coal_log_likelihood_IMRPhenomD(freq, length, waveformout, 
				noise, snr, 1.2*calculate_chirpmass(params.mass1,params.mass2), 
				calculate_eta(params.mass1,params.mass2), spin1[2], spin2[2], false,&plan);
	//params.mass1=300;
	//params.mass2=100;
	//params.spin1[2] = .9;
	//params.spin2[2] = -.2;
	start = clock();
	double logl3 = maximized_coal_Log_Likelihood(waveformout, noise,freq,length, 
				&params,"Hanford","IMRPhenomD",&plan);
	end = clock();
	cout<<"logl new TIMING: "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;
	start = clock();
	double logl4 = maximized_coal_log_likelihood_IMRPhenomD_Full_Param(freq, length, waveformout, 
				noise,  calculate_chirpmass(params.mass1,params.mass2), 
				calculate_eta(params.mass1,params.mass2), params.spin1[2], params.spin2[2],params.Luminosity_Distance,params.theta,params.phi,params.incl_angle, false,&plan);
	end = clock();
	cout.precision(15);
	cout<<"logl old old TIMING: "<<(double)(end4-start4)/CLOCKS_PER_SEC<<endl;
	cout<<"logl old  TIMING: "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;
	cout<<logl<<endl;
	cout<<logl2<<endl;
	cout<<"LOGLnew: "<<logl3<<endl;
	cout<<"LOGLold: "<<logl4<<endl;
	deactivate_likelihood_function(&plan);	



	double real_data[length];
	double imag_data[length];
	double loglpy;
	for (int i = 0; i<length; i++)
	{
		real_data[i] = real(waveformout[i]);
		imag_data[i] = imag(waveformout[i]);
	}
	start4 = clock();
	loglpy = maximized_coal_log_likelihood_IMRPhenomD(freq, length, real_data, imag_data, 
				noise, snr, calculate_chirpmass(params.mass1,params.mass2), 
				calculate_eta(params.mass1,params.mass2), spin1[2], spin2[2], false);
	end4 = clock();
	cout<<"logl TIMING with setup: "<<(double)(end4-start4)/CLOCKS_PER_SEC<<endl;
	cout<<loglpy<<endl;

//###################################################################################################
	
	//method = "ppE_IMRPhenomD_Inspiral";
	//method = "dCS_IMRPhenomD_log";
	method = "EdGB_IMRPhenomD_log";
	clock_t  startppe, endppe;
	startppe = clock(); 
	fourier_waveform(freq, length, waveformout,method,&params);
	endppe=clock();
	cout<<"TIMING waveform ppE: "<<(double)(endppe-startppe)/CLOCKS_PER_SEC<<endl;
	
	startppe = clock(); 
	fourier_amplitude(freq, length, amp,method,&params);
	endppe=clock();
	cout<<"TIMING amplitude ppE: "<<(double)(endppe-startppe)/CLOCKS_PER_SEC<<endl;

	startppe = clock(); 
	fourier_phase(freq, length, phaseout,method,&params);
	endppe=clock();
	cout<<"TIMING phase ppE: "<<(double)(endppe-startppe)/CLOCKS_PER_SEC<<endl;
	params.betappe[0]=2.;

	//ofstream ampfile;
	ampfile.open("testing/data/ppeamplitude_output.csv");
	ampfile.precision(15);
	for(int i = 0;i<length;i++)
		ampfile<<freq[i]<<','<<amp[i]<<endl;
	ampfile.close();
	
	//ofstream phasefile;
	phasefile.open("testing/data/ppephase_output.csv");
	phasefile.precision(15);
	for(int i = 0;i<length;i++){
		phasefile<<freq[i]<<','<<phaseout[i]<<endl;
		//std::cout<<phaseout[i]<<std::endl;
	}
	phasefile.close();

	//ofstream wavefilereal;
	wavefilereal.open("testing/data/ppereal_waveform_output.csv");
	wavefilereal.precision(15);
	for(int i = 0;i<length;i++)
		wavefilereal<<freq[i]<<','<<real(waveformout[i])<<endl;
	wavefilereal.close();

	//ofstream wavefileimag;
	wavefileimag.open("testing/data/ppeimag_waveform_output.csv");
	wavefileimag.precision(15);
	for(int i = 0;i<length;i++)
		wavefileimag<<freq[i]<<','<<imag(waveformout[i])<<endl;
	wavefileimag.close();
	
	
	
	
	
	int dimension = 7;
	
	double parameters[dimension] = {params.mass1,params.mass2,params.Luminosity_Distance,spin1[2],spin2[2],params.phic,params.tc};
	double **amp_derivative = (double**) malloc(dimension * sizeof(**amp_derivative));
	for (int i = 0; i<dimension;i++)
		amp_derivative[i] = (double *)malloc(length * sizeof(double)); 
	double **phase_derivative = (double**) malloc(dimension * sizeof(**phase_derivative));
	for (int i = 0; i<dimension;i++)
		phase_derivative[i] = (double *)malloc(length * sizeof(double)); 
	double spin1vec[3] = {0,0,parameters[3]};
	double spin2vec[3] = {0,0,parameters[4]};
	source_parameters<double> source_params;
	source_params = source_params.populate_source_parameters_old(parameters[0],
			parameters[1],parameters[2],spin1vec,spin2vec,parameters[5],
			parameters[6],false);

	lambda_parameters<double> lambda;
	modeld.assign_lambda_param(&source_params, &lambda);
	modeld.post_merger_variables(&source_params);
	source_params.f1 = 0.014/(source_params.M);
	source_params.f3 = modeld.fpeak(&source_params, &lambda);
	source_params.f1_phase = 0.018/(source_params.M);
	source_params.f2_phase = source_params.fRD/2.;
	source_params.bppe = new int[1];
	source_params.bppe[0] =-1;
	source_params.betappe = new double[1];
	source_params.betappe[0] = 10;
	source_params.Nmod=1;

	double A0 = source_params.A0; 
	double tc = source_params.tc;
	double phic = source_params.phic;
	double chirpmass = source_params.chirpmass;
	double symm = source_params.eta;
	double chi_s = source_params.chi_s;
	double chi_a = source_params.chi_a;


	IMRPhenomD<adouble> model;
	int amptapes[3] = {10,11,12};
	int phasetapes[3] = {13,14,15};
	model.amplitude_tape(&source_params, amptapes);
	model.phase_tape(&source_params, phasetapes);



	clock_t start5,end5;
	start5=clock();
	//for (int i = 0; i<100;i++)
	modela.construct_amplitude_derivative(freq,length,dimension,amp_derivative, &source_params); 
	
	modela.construct_phase_derivative(freq,length,dimension,phase_derivative, &source_params); 
	end5=clock();
	
	cout<<"TIMING: 2 grad: "<<(double)(end5-start5)/CLOCKS_PER_SEC<<endl;
	ofstream derivA;
	derivA.open("testing/data/deriv_amp.csv");
	derivA.precision(15);
	for(int i = 0;i<length;i++)
		derivA<<freq[i]<<','<<A0*amp_derivative[0][i]<<','<<amp_derivative[1][i]<<','<<amp_derivative[2][i]<<','<<chirpmass*amp_derivative[3][i]<<','<<symm*amp_derivative[4][i]<<','<<amp_derivative[5][i]<<','<<amp_derivative[6][i]<<endl;
	derivA.close();
	ofstream derivp;
	derivp.open("testing/data/deriv_phase.csv");
	derivp.precision(15);
	for(int i = 0;i<length;i++)
		derivp<<freq[i]<<','<<A0*phase_derivative[0][i]<<','<<phase_derivative[1][i]<<','<<phase_derivative[2][i]<<','<<chirpmass*phase_derivative[3][i]<<','<<symm*phase_derivative[4][i]<<','<<phase_derivative[5][i]<<','<<phase_derivative[6][i]<<endl;
	derivp.close();

	
	clock_t start7,end7;
	double **output = (double **)malloc(dimension * sizeof(**output));	
	for (int i = 0;i<dimension;i++)
		output[i] = (double *)malloc(dimension*sizeof(double));
	
	start7 = clock();
	fisher(freq, length, "IMRPhenomD","Hanford_O1_fitted", output, dimension, 
				&params );

	end7 = clock();
	cout<<"TIMING: FISHER: "<<(double)(end7-start7)/CLOCKS_PER_SEC<<endl;
	cout.precision(4);
	for (int i = 0;i <dimension;i++)
	{
		for (int j=0;j <dimension; j++)
			cout<<output[i][j]<<"   ";
		cout<<endl;
	}

	int dimensionppe = dimension +1;
	double **outputppe = (double **)malloc(dimensionppe * sizeof(**outputppe));	
	for (int i = 0;i<dimensionppe;i++)
		outputppe[i] = (double *)malloc(dimensionppe*sizeof(double));
	start7 = clock();
	fisher(freq, length, "ppE_IMRPhenomD_IMR","Hanford_O1_fitted", outputppe, dimensionppe, 
				&params );

	end7 = clock();
	cout<<"TIMING: FISHER ppE: "<<(double)(end7-start7)/CLOCKS_PER_SEC<<endl;
	for (int i = 0;i <dimensionppe;i++)
	{
		for (int j=0;j <dimensionppe; j++)
			cout<<outputppe[i][j]<<"   ";
		cout<<endl;
	}


	double **ppeamp_derivative = (double**) malloc(dimensionppe * sizeof(**ppeamp_derivative));
	for (int i = 0; i<dimensionppe;i++)
		ppeamp_derivative[i] = (double *)malloc(length * sizeof(double)); 
	double **ppephase_derivative = (double**) malloc(dimensionppe * sizeof(**ppephase_derivative));
	for (int i = 0; i<dimensionppe;i++)
		ppephase_derivative[i] = (double *)malloc(length * sizeof(double)); 


	ppE_IMRPhenomD_Inspiral<double> modelppe;
	start5=clock();
	//for (int i = 0; i<100;i++)
	modelppe.construct_amplitude_derivative(freq,length,dimensionppe,ppeamp_derivative, &source_params); 
	
	modelppe.construct_phase_derivative(freq,length,dimensionppe,ppephase_derivative, &source_params); 
	end5=clock();
	
	cout<<"TIMING: 2 grad ppE: "<<(double)(end5-start5)/CLOCKS_PER_SEC<<endl;
	ofstream ppederivA;
	ppederivA.open("testing/data/ppederiv_amp.csv");
	ppederivA.precision(15);
	for(int i = 0;i<length;i++)
		ppederivA<<freq[i]<<','<<A0*ppeamp_derivative[0][i]<<','<<ppeamp_derivative[1][i]<<','<<ppeamp_derivative[2][i]<<','<<chirpmass*ppeamp_derivative[3][i]<<','<<symm*ppeamp_derivative[4][i]<<','<<ppeamp_derivative[5][i]<<','<<ppeamp_derivative[6][i]<<','<<ppeamp_derivative[7][i]<<endl;
	ppederivA.close();
	ofstream ppederivp;
	ppederivp.open("testing/data/ppederiv_phase.csv");
	ppederivp.precision(15);
	for(int i = 0;i<length;i++)
		ppederivp<<freq[i]<<','<<A0*ppephase_derivative[0][i]<<','<<ppephase_derivative[1][i]<<','<<ppephase_derivative[2][i]<<','<<chirpmass*ppephase_derivative[3][i]<<','<<symm*ppephase_derivative[4][i]<<','<<ppephase_derivative[5][i]<<','<<ppephase_derivative[6][i]<<','<<ppephase_derivative[7][i]<<endl;
	ppederivp.close();
	
	
	
	
	for (int i =0;i<dimension;i++)
	{
		free( amp_derivative[i]);
		free( phase_derivative[i]);
		free(output[i]);
	}
	for (int i =0;i<dimensionppe;i++)
	{
		free(outputppe[i]);
		free( ppeamp_derivative[i]);
		free( ppephase_derivative[i]);
	}
	free(amp_derivative);
	free(freq);
	free(phase_derivative);
	free(ppeamp_derivative);
	free(ppephase_derivative);
	free(output);
	free(outputppe);
	delete [] params.betappe;
	delete [] params.bppe;
	delete [] source_params.betappe;
	delete [] source_params.bppe;
	free_LumD_Z_interp();
}

void fisher_neil_proj3 (double *pos,int dimension, double **fisher)
{
	//int alpha = (int)(gsl_rng_uniform(g)*1e7);
 	adouble* x = new adouble[dimension];
 	adouble y = 1;  
 	double out =1;
 	trace_on(1);
 	for (int i =0; i< dimension; i++){
 	        x[i]<<= pos[i];
 	}
 	y =-1* log(dist(x, dimension));
 	y>>=out;
 	delete[] x;
 	trace_off();
 	hessian(1,dimension,pos,fisher);
	for (int i = 0 ; i<dimension; i++){
        	for (int j=0;j<i;j++){
        	        if (i!=j) fisher[j][i] =fisher[i][j];
        	}
	}

}
adouble dist(adouble *pos, int dimension){
        adouble x = pos[0];
        adouble y = pos[1];
        adouble exponent_1 = - pow(x,2) - pow(9 + 4*pow(x,2) + 8*y , 2);
        adouble exponent_2 = - 8*pow(x,2) - 8*pow(y - 2, 2);
        adouble out =( 16/(3 * M_PI) ) * ( exp(exponent_1) + 0.5 * exp(exponent_2) ); 
 
        return out;
}
double log_neil_proj3 (double *c,int dim)
{
	double x = c[0];
	double y = c[1];
	double prefactor = 16./(M_PI*3.);
	double pow1 = -x*x - pow((9+4*x*x +8*y),2);
	double pow2 = -8*x*x -8*pow(y-2,2);
	return log(prefactor*(std::exp(pow1) + .5*std::exp(pow2)));
	//return 2.;
}
double log_student_t (double *x,int dim){

	double  mu=1, nu=3,  sigma=1;
        double g1 = gsl_sf_gamma( (nu + 1) / 2 ) ;
        double g2 = gsl_sf_gamma( (nu/2) );
        double parenth = 1 + (1/nu) *pow( (x[0] - mu) / sigma, 2 );
        return log(g1 / (g2 * sqrt(nu * M_PI ) * sigma ) * pow(parenth,-(nu+1)/2    ));
}
void test_fisher(double *pos, int dim, double **fisher)
{
	fisher[0][0] = .5;
}
double test_ll(double *pos, int dim)
{
	//std::cout<<"LL"<<std::endl;
	//std::cout<<"Pos in LL: "<<pos[0]<<std::endl;
	return -pos[0]*pos[0]/(4.);
	//return  0;
}
double test_lp(double *pos, int dim)
{
	//return 0;
	return -pos[0]*pos[0]/(10.)- pos[1]*pos[1]/20.;
}	
double test_lp_GW(double *pos, int dim)
{
	double a = std::numeric_limits<double>::infinity();
	//Flat priors across physical regions
	if (std::exp(pos[0])/MPC_SEC<50 || std::exp(pos[0])/MPC_SEC>1000){return a;}
	if (std::exp(pos[1])/MSOL_SEC<2 || std::exp(pos[1])/MSOL_SEC>100){return a;}
	if ((pos[2])<.1 || (pos[2])>.245){return a;}
	if ((pos[3])<-.9 || (pos[3])>.9){return a;}
	if ((pos[4])<-.9 || (pos[4])>.9){return a;}
	else {return 0.;}
}
