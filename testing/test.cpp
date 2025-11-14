#include <iostream>
#include <fstream>
#include <time.h>
#include <complex>
#include <string>
#include "gwat/waveform_generator.h"
#include "gwat/io_util.h"
#include "gwat/autocorrelation.h"
#include "gwat/IMRPhenomD.h"
#include "gwat/mcmc_gw.h"
#include "gwat/mcmc_sampler_internals.h"
#include "gwat/detector_util.h"
#include "gwat/util.h"
#include "gwat/waveform_util.h"
#include <adolc/adouble.h>
#include "gwat/fisher.h"
#include "gwat/ppE_IMRPhenomD.h"
#include "gwat/IMRPhenomP.h"
#include "gwat/ppE_IMRPhenomP.h"
#include "gwat/waveform_generator_C.h"
#include "gwat/mcmc_sampler.h"
#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_rng.h"
#include "adolc/adouble.h"
#include "adolc/drivers/drivers.h"
#include "adolc/taping.h"
#include "limits"
#include "gwat/ortho_basis.h"

#include "gwat/autocorrelation_cuda.h"

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <algorithm>

using namespace std;

bool show_progress = true;
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
void test11();
void test12();
void test13();
void test14();
void test15();
void test16();
void test17();
void test18();
void test19();
void test20();
void test21();
void test22();
void test23();
void test24();
void test25();
void test26();
void test27();
void test28();
void test29();
void test30();
void test31();
void test32();
void test33();
void test34();
void test35();
void test36();
void test37();
void test38();
void test39();
void test40();
void test41();
void test42();
void test43();
void test45();
void test46();
void test47();
void test48();
void test49();
void test50();
void test51();
void test52();
void test53();
void test54();
void test55();
void test56();
void test57();
void test58();
void test59();
void test60();
void test61();
void test62();
void test63();
void test64();
void test_prob(double *prob, void *param, int d, int threadid);
double test_ll(double *pos, int dim,void *parameters);
double test_lp(double *pos, int dim,void *parameters);
double test_lp_nts(double *pos, int dim, int chain_id,void *parameters);
double test_lp_GW(double *pos, int dim, int chain_id,void *parameters);
double test_lp_GW_Pv2(double *pos, int dim, int chain_id,void *parameters);
double test_lp_GW_D(double *pos, int dim, int chain_id,void *parameters);
double test_lp_GW_Pv2_in(double *pos, int dim, int chain_id,void *parameters);
double test_lp_GW_Pv2_ppE(double *pos, int dim, int chain_id,void *parameters);
double test_lp_GW_Pv2_dCS_root_alpha(double *pos, int dim, int chain_id,void *parameters);
double test_lp_GW_7dim(double *pos, int dim, int chain_id,void *parameters);
double test_lp_GW_DFull(double *pos, int dim, int chain_id,void *parameters);
double test_lp_GW_dCS(double *pos, int dim, int chain_id,void *parameters);
double test_lp_GW_dCS_log(double *pos, int dim, int chain_id,void *parameters);
double test_lp_GW_dCS_root_alpha(double *pos, int dim , int chain_id,void *parameters);
double test_lp_GW_ppE(double *pos, int dim, int chain_id,void *parameters);
double test_lp_reduceddim(double *pos, int dim, int chain_id,void *parameters);
void test_fisher(double *pos, int dim, double **fisher,void *parameters);
double log_student_t (double *x,int dim,void *parameters);
double log_neil_proj3 (double *x,int dim,void *parameters);
double log_neil_proj3_nts (double *x,int dim, int chain_id,void *parameters);
double log_neil_proj32 (double *c,int dim,void *parameters);
void fisher_neil_proj3 (double *x,int dim, double **fish, int chain_id,void *parameters);
double chirpmass_eta_jac(double chirpmass, double eta);
double test_lp_GW_D_in(double *pos, int dim, int chain_id,void *parameters);
adouble dist(adouble *pos, int dimension);

const gsl_rng_type* Y;
gsl_rng * g;
static std::complex<double> *waveformout11=NULL;
static double *freq=NULL;
static double *psd=NULL;

int main(){

	//test38();	
	//test54();	
	//test54();	
	//test6();	
	//test45();	
	//test62();	
	test63();	
	return 0;
}

void test64()
{
	int length = 4000;
	int num_detectors =3;
	//int num_detectors =2;
	std::string *detectors = new std::string[num_detectors];//(std::string*)malloc(sizeof(std::string)*50*num_detectors);
	detectors[0] = "Hanford";
	detectors[1] = "Livingston";
	detectors[2] = "Virgo";

	double gps_time = 1135136350.6;//TESTING -- gw151226

	std::complex<double> **data= (std::complex<double>**)malloc(
			sizeof(std::complex<double>*)*num_detectors);
	double **psd = (double **)malloc(sizeof(double *)*num_detectors);
	double **frequencies = (double **)malloc(sizeof(double *)*num_detectors);
	int *data_length= (int*)malloc(sizeof(int)*num_detectors);
	data_length[0] =length;
	data_length[1] =length;
	data_length[2] =length;

	for (int i =0; i<num_detectors; i++){
		data[i] = (std::complex<double> *)malloc(
			sizeof(std::complex<double>)*data_length[i]);
		
		psd[i] = (double *)malloc(sizeof(double)*data_length[i]);
		frequencies[i] = (double *)malloc(sizeof(double)*data_length[i]);
	}
	//#########################################################
	//Make trial data
	gen_params params;
	//double RA = 5.;
	//double DEC = 1.;
	double RA = 2;
	double DEC = .4;
	double chirpm = 9.71;
	double eta =.23;
	params.mass1 = calculate_mass1(chirpm,eta);
	params.mass2 = calculate_mass2(chirpm,eta);
	complex<double> waveformout[length];
	params.spin1[0] = 0;
	params.spin1[1] = .2;
	params.spin1[2] = .43;
	params.spin2[0] = .1;
	params.spin2[1] = 0;
	params.spin2[2] = .2;
	params.phiRef = 1.0;
	double tc = 3;
	params.Luminosity_Distance = 125.;
	params.NSflag1 = false;
	params.NSflag2 = false;
	params.shift_time = true;
	params.shift_phase = true;
	params.incl_angle = .2;//M_PI/3.;
	params.sky_average=false;
	std::string injection_method = "IMRPhenomD";
	params.f_ref=20;
	params.RA = RA;
	params.gmst = gps_to_GMST_radian(gps_time);
	//params.f_ref = 30.5011;
	//params.phiRef =58.944425/2.;
	
	//############################################################
	double fhigh =2048;
	double flow =20;
	double df = (fhigh-flow)/(length-1);
	//double *freq = (double *)malloc(sizeof(double) * length);
	double freq[length];

	cout<<"Freq spacing "<<df<<endl;

	for(int i=0;i<length;i++)
		freq[i]=flow+i*df;
	//############################################################
	double noise[length];
	populate_noise(freq,"AdLIGODesign", noise,length);
	for (int i =0; i<length;i++){
		noise[i] = noise[i]*noise[i];
	}
	//############################################################

	for(int i = 0 ; i<num_detectors; i++){
		double dt = DTOA_DETECTOR(params.RA,params.DEC,params.gmst,"Hanford",detectors[i]);
		params.tc = dt+tc; 
		fourier_detector_response(freq,length, data[i],detectors[i],injection_method, &params);

		for(int j = 0; j<data_length[i]; j++){
			frequencies[i][j] = freq[j];	
			psd[i][j] = (noise[j]);	
		}
	}

	//#########################################################
	//mcmc options
	int dimension = 14;
	double initial_pos[dimension]={params.RA,sin(params.DEC),params.psi,cos(params.incl_angle),params.phiRef,params.tc,std::log(params.Luminosity_Distance),std::log(chirpm), eta,.1,.1,.9,.9,.1};
	//int dimension = 4;
	//double initial_pos[dimension]={std::log(30), .24,.1,.1};
	double *seeding_var = NULL;
	int n_steps = 5000;
	int chain_N=24 ;
	int max_thermo=24 ;
	int t0 = 1000;
	int nu = 100;
	std::string chain_alloc = "double";
	int swp_freq = 3;
	double chain_temps[chain_N];
	double c = 1.3;
	
	int Nmod = 0;
	int *bppe = NULL;
	int numThreads = 10;
	bool pool = true;
	//#########################################################
	//gw options
	//std::string generation_method = "dCS_IMRPhenomD_log";
	//std::string generation_method = "EdGB_IMRPhenomD_log";
	//std::string generation_method = "dCS_IMRPhenomD_root_alpha";
	//std::string generation_method = "IMRPhenomD";
	std::string generation_method = "IMRPhenomPv2";
	//std::string generation_method = "EdGB_IMRPhenomD_root_alpha";
	
	
	std::string chainfile = "testing/data/mcmc_output_uncorr_Pv2_in.csv";
	std::string statfilename = "testing/data/mcmc_statistics_uncorr_Pv2_in.txt";
	std::string checkfile = "testing/data/mcmc_checkpoint_uncorr_Pv2_in.csv";

	int corr_threshold = 10;
	int corr_segments = 10;
	double corr_converge_thresh = 0.1;
	double corr_target_ac = .01;

	double **output;
	output = allocate_2D_array(n_steps, dimension );
	//double ***output;
	//output = allocate_3D_array(chain_N,n_steps, dimension );
	PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW(output, dimension, n_steps, chain_N, max_thermo, initial_pos,seeding_var,chain_temps, 
			swp_freq, t0, nu, corr_threshold, corr_segments, corr_converge_thresh, corr_target_ac,chain_alloc, test_lp_GW_Pv2,numThreads, pool,show_progress,
			num_detectors, 
			data, psd,frequencies, data_length,gps_time, detectors,Nmod, bppe,
			generation_method,statfilename,chainfile, "",checkfile);	
	//std::string checkfile2="testing/data/mcmc_checkpoint_uncorr_D2.csv";
	//continue_PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW(checkfile,output,  n_steps,  max_thermo, chain_temps, 
	//		swp_freq, t0, nu, corr_threshold, corr_segments, corr_converge_thresh, corr_target_ac,chain_alloc, test_lp_GW_D,numThreads, pool,show_progress,
	//		num_detectors, 
	//		data, psd,freqs, data_length,gps_time, detectors,Nmod, bppe,
	//		generation_method,statfilename,chainfile, "",checkfile2);	
	//PTMCMC_MH_dynamic_PT_alloc_GW(output, dimension, n_steps,chain_N,max_thermo, initial_pos, seeding_var, chain_temps, swp_freq, t0,nu,"half_ensemble",test_lp_GW_D,numThreads, pool, show_progress, num_detectors, data, psd, freqs, data_length, gps_time, detectors, Nmod, bppe, generation_method, statfilename, chainfile,  "",checkfile);


	//write_file(chainfile, output[0], n_steps, dimension);

	deallocate_2D_array(output, n_steps, dimension);
	//deallocate_3D_array(output, chain_N,n_steps, dimension);
	delete [] detectors;
	free(data_length);
	//free_LOSC_data(data, psd,freqs, num_detectors, length);
	//deallocate_2D_array(psd,num_detectors, data_length[0]);
	//deallocate_2D_array(freq,num_detectors, data_length[0]);
	for(int i =0; i<num_detectors; i++)
		free(data[i]);
	free(data);
	//delete [] detector_files;

}
void test63()
{
	//Terr
	int length = 10000;
	int curves = 14;
	double *f = new double[length];
	double **psd = new double*[curves+1];
	std::string psds[curves] = {"aLIGO_analytic","Hanford_O1_fitted","AdLIGOMidHigh","AdLIGOAPlus","AdLIGODesign","CE1","CE2","KAGRA_opt","KAGRA_pess","AdVIRGOPlus1","AdVIRGOPlus2_opt","AdVIRGOPlus2_pess","AdLIGODesign_smoothed","CE1_smoothed"};
	psd[0]=new double[length];
	for(int i = 0 ; i<length; i++){
		f[i]=15. + 3000./length * i;
		psd[0][i] = f[i];
		
	}
	for(int i = 0 ; i<curves; i++){
		psd[i+1] = new double[length];
		populate_noise(f, psds[i], psd[i+1], length);
	}
	write_file("testing/data/noise_curves_terr.csv",psd,curves+1,length);
	delete [] f;
	for(int i = 0 ; i<curves+1; i++){
		delete [] psd[i];	
	}
	delete [] psd;
	//Space
	std::string psds_space[curves] = {"LISA","LISA_CONF","LISA_SADC","LISA_SADC_CONF"};
	curves = 4;
	f = new double[length];
	psd = new double*[curves+1];
	psd[0]=new double[length];
	for(int i = 0 ; i<length; i++){
		f[i]=1e-5 + 1./length * i;
		psd[0][i] = f[i];
	}
	for(int i = 0 ; i<curves; i++){
		psd[i+1] = new double[length];
		populate_noise(f, psds_space[i], psd[i+1], length);
	}
	write_file("testing/data/noise_curves_space.csv",psd,curves+1,length);
	delete [] f;
	for(int i = 0 ; i<curves; i++){
		delete [] psd[i];	
	}
	delete [] psd;
}

void test62()
{
	std::cout.precision(15);
	gen_params params;
	params.mass1 = 10e5;
	params.mass2 = 5e5;
	params.spin1[2] = .2;
	params.spin2[2] = .1;
	params.spin1[1] = .2;
	params.spin2[1] = .1;
	params.spin1[0] = .2;
	params.spin2[0] = .1;
	params.phiRef = 2.;
	params.tc = T_year;
	params.f_ref = 1e-5;
	params.NSflag1=false;
	params.NSflag2=false;
	params.shift_phase=false;
	params.shift_time=false;
	params.equatorial_orientation=true;
	params.horizon_coord=false;
	params.Luminosity_Distance  = 400;
	params.RA =1.;
	params.DEC =-1.;
	//params.psi =2.;
	//params.incl_angle =2.;
	params.phi_l =2.;
	params.theta_l =2.;
	params.gmst=1.;
	params.LISA_phi0=0;
	params.LISA_alpha0=0;
	
	
	int length = T_year;
	double deltaf = 1./(T_year);
	double *freqs = new double[length];
	//for(int i = 0 ; i<length; i++){
	//	freqs[i] = 1e-5 + i*deltaf;
	//}
	
	int lengthgl1 = 5000;
	int lengthgl2 = 5000;
	double *freqsgl1 = new double[lengthgl1];
	double *freqsgl2 = new double[lengthgl2];
	double *w1 = new double[lengthgl1];
	double *w2 = new double[lengthgl2];

	gauleg(1e-5, 1, freqsgl1, w1, lengthgl1);
	gauleg(log10(1e-5), log10(1), freqsgl2, w2, lengthgl2);
	for(int i =0; i<lengthgl2; i++){
		freqsgl2[i] = pow(10.,freqsgl2[i]);
	}

	gsl_rng_env_setup();
	const gsl_rng_type *T = gsl_rng_default;
	gsl_rng * r = gsl_rng_alloc(T);
	gsl_rng_set(r, 1);


	double snr_gsl, snr_simps, snr_gl1,snr_gl2;
	std::string noise_curve="LISA_CONF"	;
	std::string method= "IMRPhenomPv2";
	std::string detector="LISA";
	int iterations = 10;
	double **output = new double*[iterations];
	std::cout<<"Starting loop -- finished prep"<<std::endl;
	for(int i = 0 ; i<iterations; i++){
		output[i]=new double[8];
		double m1 = 10e4+1e4*gsl_rng_uniform(r);
		double m2 = 10e4+1e4*gsl_rng_uniform(r);
		if(m1>m2){
			params.mass1 = m1;
			params.mass2 = m2;
		}
		else{
			params.mass1 = m1;
			params.mass2 = m2;
		}
		params.spin1[2] = -.3+gsl_rng_uniform(r)*.5;
		params.spin2[2] = -.3+gsl_rng_uniform(r)*.5;
		params.spin1[1] = -.4 + gsl_rng_uniform(r)*.8;
		params.spin2[1] = .1;
		params.spin1[0] = .2;
		params.spin2[1] = 0;//-.4 + gsl_rng_uniform(r)*1.;
		params.Luminosity_Distance  = 50 + gsl_rng_uniform(r)*1000;
		params.RA =gsl_rng_uniform(r)*2*M_PI;
		params.DEC =-1.5 + gsl_rng_uniform(r) * M_PI;
		//params.psi =gsl_rng_uniform(r)*2*M_PI;
		//params.incl_angle =gsl_rng_uniform(r)*M_PI;
		params.phi_l =gsl_rng_uniform(r)*2*M_PI;
		params.theta_l =gsl_rng_uniform(r)*M_PI;
		params.gmst=gsl_rng_uniform(r)*2*M_PI;
		//params.chip = -1;
		transform_orientation_coords(&params,method, detector);

		double bounds[2];
		integration_interval(1, T_year, detector,noise_curve, method, &params, bounds);
		//std::cout<<bounds[0]<<" "<<bounds[1]<<std::endl;

		for(int i = 0 ; i<length; i++){
			freqs[i] = bounds[0] + i*deltaf;
		}


		clock_t start = clock();
		calculate_snr_gsl(&snr_gsl, noise_curve,detector,method, &params, bounds[0], bounds[1], 1.e-12);
		//snr_gsl*=sqrt(2.);
		double gsl_t = (double)(clock()-start)/CLOCKS_PER_SEC;
		start = clock();
		//snr_simps = calculate_snr(noise_curve, detector,method, &params,freqs,length,"SIMPSONS",(double *)NULL,false);
		double simps_t = (double)(clock()-start)/CLOCKS_PER_SEC;
		start = clock();
		gauleg(bounds[0], bounds[1], freqsgl1, w1, lengthgl1);
		snr_gl1 = calculate_snr(noise_curve, detector,method, &params,freqsgl1,lengthgl1,"GAUSSLEG",w1,false);
		double gl1_t = (double)(clock()-start)/CLOCKS_PER_SEC;
		start = clock();
		gauleg(log10(bounds[0]), log10(bounds[1]), freqsgl2, w2, lengthgl2);
		for(int i =0; i<lengthgl2; i++){
			freqsgl2[i] = pow(10.,freqsgl2[i]);
		}
		snr_gl2 = calculate_snr(noise_curve, detector,method, &params,freqsgl2,lengthgl2,"GAUSSLEG",w2,true);
		double gl2_t = (double)(clock()-start)/CLOCKS_PER_SEC;
		std::cout<<snr_gsl<<" "<<snr_simps<<" "<<snr_gl1<<" "<<snr_gl2<<std::endl;
		//std::cout<<params.mass1<<" "<<params.mass2<<" "<<params.spin1[0]<<" "<<params.spin1[1]<<" "<<params.spin1[2]<<" "<<params.spin2[0]<<" "<<params.spin2[1]<<" "<<params.spin2[2]<<std::endl;
		//std::cout<<" "<<(snr_gsl-snr_simps)/snr_gsl<<" "<<(snr_gsl-snr_gl1)/snr_gsl<<" "<<(snr_gsl-snr_gl2)/snr_gsl<<std::endl;
		//std::cout<<gsl_t<<" "<<simps_t<<" "<<gl1_t<<" "<<gl2_t<<std::endl;
		output[i][0]=snr_gsl;
		output[i][1]=snr_simps;
		output[i][2]=snr_gl1;
		output[i][3]=snr_gl2;
		output[i][4]=gsl_t;
		output[i][5]=simps_t;
		output[i][6]=gl1_t;
		output[i][7]=gl2_t;
		printProgress((double)i/iterations);
	}
	std::cout<<std::endl;
	
	write_file("testing/data/snr_comp.csv",output,iterations, 8);
	delete [] freqsgl1;
	delete [] freqsgl2;
	delete [] w1;
	delete [] w2;
	for(int i = 0 ; i<iterations; i++){
		delete [] output[i];	
	}
	delete [] output;
}

void test61()
{
	int datalength = 131075;

	std::string psd_file = "/Users/sperkins/Downloads/LOSC_data/GW170729/GWTC1_GW170729_PSDs.dat.txt";
	int num_detectors = 3, psd_length = 4016, length;
	double gps_time = 1185389807.3;//TESTING -- gw170729
	std::string *detectors = new std::string[num_detectors];//(std::string*)malloc(sizeof(std::string)*50*num_detectors);
	detectors[0] = "Hanford";
	detectors[1] = "Livingston";
	detectors[2] = "Virgo";
	std::string *detector_files = new std::string[num_detectors];
	detector_files[0] =  "/Users/sperkins/Downloads/LOSC_data/GW170729/H-H1_GWOSC_4KHZ_R1-1185389792-32.txt";
	detector_files[1] =  "/Users/sperkins/Downloads/LOSC_data/GW170729/L-L1_GWOSC_4KHZ_R1-1185389792-32.txt";
	detector_files[2] =  "/Users/sperkins/Downloads/LOSC_data/GW170729/V-V1_GWOSC_4KHZ_R1-1185389792-32.txt";
	//
	//std::string psd_file = "/home/sperkins/Downloads/LOSC_data/GW170729/GWTC1_GW170729_PSDs.dat.txt";
	//int num_detectors = 2, psd_length = 8032, length;
	//double gps_time = 1126259462.4;//TESTING -- gw150914
	//std::string *detectors = new std::string[num_detectors];//(std::string*)malloc(sizeof(std::string)*50*num_detectors);
	//detectors[0] = "Hanford";
	//detectors[1] = "Livingston";
	//std::string *detector_files = new std::string[num_detectors];
	//detector_files[0] =  "/Users/sperkins/Downloads/LOSC_data/GW150914/H-H1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[1] =  "/Users/sperkins/Downloads/LOSC_data/GW150914/L-L1_GWOSC_4KHZ_R1-1126259447-32.txt";
 	double trigger_time = gps_time;
	double **psd = allocate_2D_array(num_detectors,psd_length);
	double **freqs = allocate_2D_array(num_detectors,psd_length);
	std::complex<double> **data = (std::complex<double> **)malloc(sizeof(std::complex<double> *)*num_detectors);
	for(int i =0; i<num_detectors; i++)
		data[i] = (std::complex<double>*)malloc(sizeof(std::complex<double>)*psd_length);

	allocate_LOSC_data(detector_files, psd_file, num_detectors, psd_length, datalength, trigger_time, data, psd, freqs);

	//double **temp = allocate_2D_array(psd_length,3);
	//
	//for(int i=0 ; i<psd_length; i++){
	//	std::cout<<psd[0][i]<<data[0][i]<<std::endl;
	//	temp[i][0]=freqs[0][i];	
	//	temp[i][1]=std::real( data[0][i]/sqrt(psd[0][i]));
	//	temp[i][2]=std::imag( data[0][i]/sqrt(psd[0][i]));
	//}
	//write_file("whitened_GW150914_Hanford.csv",temp,psd_length,3);
	//deallocate_2D_array(temp,psd_length,3);

	int *data_length= (int*)malloc(sizeof(int)*num_detectors);
	data_length[0] =psd_length;
	data_length[1] =psd_length;
	data_length[2] =psd_length;

	//#########################################################
	//mcmc options
	int dimension = 7;
	double initial_pos[dimension]={std::log(30), .24,.1,.1,.9,.9,.1};
	//int dimension = 4;
	//double initial_pos[dimension]={std::log(30), .24,.1,.1};
	double *seeding_var = NULL;
	int n_steps = 10000;
	int chain_N=12 ;
	int max_thermo=12 ;
	int t0 = 1000;
	int nu = 100;
	std::string chain_alloc = "double";
	int swp_freq = 5;
	double chain_temps[chain_N];
	double c = 1.3;
	
	int Nmod = 0;
	int *bppe = NULL;
	int numThreads = 10;
	bool pool = true;
	//#########################################################
	//gw options
	//std::string generation_method = "dCS_IMRPhenomD_log";
	//std::string generation_method = "EdGB_IMRPhenomD_log";
	//std::string generation_method = "dCS_IMRPhenomD_root_alpha";
	//std::string generation_method = "IMRPhenomD";
	std::string generation_method = "IMRPhenomPv2";
	//std::string generation_method = "EdGB_IMRPhenomD_root_alpha";
	
	
	std::string chainfile = "testing/data/mcmc_output_uncorr_Pv2_in.csv";
	std::string statfilename = "testing/data/mcmc_statistics_uncorr_Pv2_in.txt";
	std::string checkfile = "testing/data/mcmc_checkpoint_uncorr_Pv2_in.csv";

	int corr_threshold = 10;
	int corr_segments = 10;
	double corr_converge_thresh = 0.05;
	double corr_target_ac = .01;

	double **output;
	output = allocate_2D_array(n_steps, dimension );
	//double ***output;
	//output = allocate_3D_array(chain_N,n_steps, dimension );
	PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW(output, dimension, n_steps, chain_N, max_thermo, initial_pos,seeding_var,chain_temps, 
			swp_freq, t0, nu, corr_threshold, corr_segments, corr_converge_thresh, corr_target_ac,chain_alloc, test_lp_GW_Pv2_in,numThreads, pool,show_progress,
			num_detectors, 
			data, psd,freqs, data_length,gps_time, detectors,Nmod, bppe,
			generation_method,statfilename,chainfile, "",checkfile);	
	//std::string checkfile2="testing/data/mcmc_checkpoint_uncorr_D2.csv";
	//continue_PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW(checkfile,output,  n_steps,  max_thermo, chain_temps, 
	//		swp_freq, t0, nu, corr_threshold, corr_segments, corr_converge_thresh, corr_target_ac,chain_alloc, test_lp_GW_D,numThreads, pool,show_progress,
	//		num_detectors, 
	//		data, psd,freqs, data_length,gps_time, detectors,Nmod, bppe,
	//		generation_method,statfilename,chainfile, "",checkfile2);	
	//PTMCMC_MH_dynamic_PT_alloc_GW(output, dimension, n_steps,chain_N,max_thermo, initial_pos, seeding_var, chain_temps, swp_freq, t0,nu,"half_ensemble",test_lp_GW_D,numThreads, pool, show_progress, num_detectors, data, psd, freqs, data_length, gps_time, detectors, Nmod, bppe, generation_method, statfilename, chainfile,  "",checkfile);


	//write_file(chainfile, output[0], n_steps, dimension);

	deallocate_2D_array(output, n_steps, dimension);
	//deallocate_3D_array(output, chain_N,n_steps, dimension);
	delete [] detectors;
	free(data_length);
	//free_LOSC_data(data, psd,freqs, num_detectors, length);
	deallocate_2D_array(psd,num_detectors, psd_length);
	deallocate_2D_array(freqs,num_detectors, psd_length);
	for(int i =0; i<num_detectors; i++)
		free(data[i]);
	free(data);
	delete [] detector_files;

}
void test60()
{
	int psd_num = 5;
	std::string psds[psd_num]={"AdLIGODesign","AdLIGOAPlus","KAGRA","CE1","CE2"};
	int length = 20000;
	double **output = new double*[psd_num+1];
	output[0]=new double[length];
	double fmin = 10;
	double fmax = 5000;
	double deltaf = (fmax-fmin)/length;
	for(int i = 0 ; i<length; i++){
		output[0][i] = fmin+i*deltaf;
	}
	for(int i = 0 ; i<psd_num; i++){
		output[i+1] = new double[length];
		populate_noise(output[0],psds[i],output[i+1],length);

	}
	write_file("testing/data/terr_psds.csv",output,  psd_num+1,length);
	for(int i = 0 ; i<psd_num+1; i++){
		delete [] output[i];
	}
	delete [] output;
}
void test59()
{
	std::string detector1("Hanford");
	std::string detector2("Livingston");
	double gps_time = 100000;
	double gmst_rad = gps_to_GMST_radian(gps_time);
	double RA = 2.8;
	double DEC = 0.8;
	double phi1,theta1,phi2,theta2;
	celestial_horizon_transform(RA,DEC,gps_time,detector1,&phi1,&theta1);	
	celestial_horizon_transform(RA,DEC,gps_time,detector2,&phi2,&theta2);	

	std::cout<<DTOA_DETECTOR(RA,DEC,gmst_rad,detector1,detector2)<<std::endl;;
	std::cout<<DTOA(theta1,theta2,detector1,detector2)<<std::endl;;
}

void test58()
{
	std::cout.precision(15);
	int length = 200;
	double **output = allocate_2D_array(length, 2);
	int N_detectors = 5;
	//std::string detectors[N_detectors] = {"Hanford","Livingston","Virgo","Kagra","Indigo"};
	std::string detectors[N_detectors] = {"Hanford","Livingston","Virgo","Indigo"};
	//std::string detectors[N_detectors] = {"Hanford","Hanford","Hanford"};
	#pragma omp parallel for
	for(int i = 0 ; i < length; i++){
		double omega = .0+.01*i;
		int samples = 5e5;
		double P=p_N_detector(omega,samples,N_detectors,detectors, i*2);
		std::cout<<"Omega: "<<omega<<"  Detection Probability: "<<P<<std::endl;
		output[i][0] = omega;
		output[i][1] = P;
	}
	write_file("testing/data/p_test.csv",output,length,2);
	deallocate_2D_array(output,length,2);

}
void test57()
{
	sampler sampler;
	sampler.max_dim = 2;
	sampler.min_dim = sampler.max_dim;
	sampler.chain_N = 128;
	sampler.chain_temps= new double[sampler.chain_N];
	sampler.N_steps= 100;
	sampler.fisher_exist=false;
	sampler.dimension=sampler.max_dim;
	sampler.num_threads=10;
	sampler.output = allocate_3D_array(sampler.chain_N, sampler.N_steps, sampler.dimension);
	sampler.param_status = allocate_3D_array_int(sampler.chain_N, sampler.N_steps, sampler.dimension);
	int ensemble_number = 15;
	for(int i = 0 ; i<sampler.chain_N; i++){
		sampler.chain_temps[i] = 1 + i%ensemble_number *2;
	}
	allocate_sampler_mem(&sampler);
	for(int i = 0 ; i<sampler.chain_N; i++){
		for(int j =0; j<sampler.dimension; j++){
			sampler.output[i][0][j]=sampler.chain_temps[i];
			sampler.param_status[i][0][j]=1;
		}
		sampler.chain_pos[i] = 0;
		sampler.current_hist_pos[i] = 0;
		for(int j =0 ; j<sampler.history_length ;j++){
			for(int k = 0 ; k<sampler.dimension; k++){

				sampler.history[i][j][k] = j;
			}
		}
		sampler.de_primed[i]=true;
	}
	write_checkpoint_file(&sampler, "testing/data/trial_checkpoint_file.csv");


	//Change Data
	ensemble_number = 13;
	for(int i = 0 ; i<sampler.chain_N; i++){
		sampler.chain_temps[i] = 1 + i%ensemble_number *3;
	}
	allocate_sampler_mem(&sampler);
	for(int i = 0 ; i<sampler.chain_N; i++){
		for(int j =0; j<sampler.dimension; j++){
			sampler.output[i][0][j]=sampler.chain_temps[i]+1;
			sampler.param_status[i][0][j]=1;
		}
		sampler.chain_pos[i] = 0;
		sampler.current_hist_pos[i] = 0;
		for(int j =0 ; j<sampler.history_length ;j++){
			for(int k = 0 ; k<sampler.dimension; k++){

				sampler.history[i][j][k] = i+1;
			}
		}
		sampler.de_primed[i]=true;
	}
	write_checkpoint_file(&sampler, "testing/data/trial_checkpoint_file_changed_data.csv");
	
	copy_base_checkpoint_properties("testing/data/trial_checkpoint_file.csv",&sampler);	
	write_checkpoint_file(&sampler, "testing/data/trial_checkpoint_file_changed_data_back.csv");

	deallocate_sampler_mem(&sampler);

	delete [] sampler.chain_temps;
	deallocate_3D_array(sampler.output, sampler.chain_N, sampler.N_steps, sampler.dimension);
	//deallocate_3D_array(sampler.param_status, sampler.chain_N, sampler.N_steps, sampler.dimension);
}
void test56()
{
	//std::string psd_file = "testing/data/GWTC1_GW150914_PSDs.dat.txt";
	//std::string psd_file = "testing/data/GWTC1_GW151226_PSDs.dat.txt";
	//std::string data_file = "testing/data/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//std::string psd_file = "testing/data/GWTC1_GW170729_PSDs.dat.txt";
	//std::string psd_file = "/home/sperkins/Downloads/LOSC_data/GW151226/GWTC1_GW151226_PSDs.dat.txt";
	//std::string psd_file = "/Users/sperkins/Downloads/LOSC_data/GW151226/GWTC1_GW151226_PSDs.dat.txt";
	//std::string psd_file = "/Users/sperkins/Downloads/LOSC_data/GW150914/GWTC1_GW150914_PSDs.dat.txt";
	std::string psd_file = "/Users/sperkins/Downloads/LOSC_data/GW150914/GWTC1_GW150914_PSDs.dat.txt";
	//std::string psd_file = "/home/sperkins/Downloads/LOSC_data/GW170729/GWTC1_GW170729_PSDs.dat.txt";
	//std::string psd_file = "/Users/sperkins/Downloads/LOSC_data/GW170729/GWTC1_GW170729_PSDs.dat.txt";
	//int rows = 8032;
	//int cols = 3;
	int datalength = 131075;
	//double **psd = allocate_2D_array(rows, cols);
	//read_LOSC_PSD_file(psd_file, psd, rows, cols);
	//double data_start_time, duration, fs;
	int num_detectors = 2, psd_length = 8032, length;
	//int num_detectors = 3, psd_length = 4016, length;
	//int num_detectors = 3, psd_length = 4016, length;
	//int num_detectors = 2, psd_length = 4016, length;
	//double gps_time = 1135136350.6;//TESTING -- gw151226
	double gps_time = 1126259462.4;//TESTING -- gw150914
	//double gps_time = 1185389807.3;//TESTING -- gw170729
	//double gps_time = 1185389807.3;//TESTING -- gw170729
	std::string *detectors = new std::string[num_detectors];//(std::string*)malloc(sizeof(std::string)*50*num_detectors);
	detectors[0] = "Hanford";
	detectors[1] = "Livingston";
	//detectors[2] = "Virgo";
	std::string *detector_files = new std::string[num_detectors];
	//detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[0] =  "/home/sperkins/Downloads/LOSC_data/GW151226/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//detector_files[1] =  "/home/sperkins/Downloads/LOSC_data/GW151226/L-L1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//detector_files[0] =  "/Users/sperkins/Downloads/LOSC_data/GW151226/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//detector_files[1] =  "/Users/sperkins/Downloads/LOSC_data/GW151226/L-L1_GWOSC_4KHZ_R1-1135136335-32.txt";
	detector_files[0] =  "/Users/sperkins/Downloads/LOSC_data/GW150914/H-H1_GWOSC_4KHZ_R1-1126259447-32.txt";
	detector_files[1] =  "/Users/sperkins/Downloads/LOSC_data/GW150914/L-L1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[0] =  "/home/sperkins/Downloads/LOSC_data/GW150914/H-H1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[1] =  "/home/sperkins/Downloads/LOSC_data/GW150914/L-L1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[0] =  "/home/sperkins/Downloads/LOSC_data/GW170729/H-H1_GWOSC_4KHZ_R1-1185389792-32.txt";
	//detector_files[1] =  "/home/sperkins/Downloads/LOSC_data/GW170729/L-L1_GWOSC_4KHZ_R1-1185389792-32.txt";
	//detector_files[2] =  "/home/sperkins/Downloads/LOSC_data/GW170729/V-V1_GWOSC_4KHZ_R1-1185389792-32.txt";
	//detector_files[0] =  "/Users/sperkins/Downloads/LOSC_data/GW170729/H-H1_GWOSC_4KHZ_R1-1185389792-32.txt";
	//detector_files[1] =  "/Users/sperkins/Downloads/LOSC_data/GW170729/L-L1_GWOSC_4KHZ_R1-1185389792-32.txt";
	//detector_files[2] =  "/Users/sperkins/Downloads/LOSC_data/GW170729/V-V1_GWOSC_4KHZ_R1-1185389792-32.txt";
	//detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1185389792-32.txt";
	//detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1185389792-32.txt";
	//detector_files[2] =  "testing/data/V-V1_GWOSC_4KHZ_R1-1185389792-32.txt";
 	//double trigger_time= 1135136350.6;
 	double trigger_time = gps_time;
	double **psd = allocate_2D_array(num_detectors,psd_length);
	double **freqs = allocate_2D_array(num_detectors,psd_length);
	std::complex<double> **data = (std::complex<double> **)malloc(sizeof(std::complex<double> *)*num_detectors);
	for(int i =0; i<num_detectors; i++)
		data[i] = (std::complex<double>*)malloc(sizeof(std::complex<double>)*psd_length);

	allocate_LOSC_data(detector_files, psd_file, num_detectors, psd_length, datalength, trigger_time, data, psd, freqs);

	//double **temp = allocate_2D_array(psd_length,3);
	//
	//for(int i=0 ; i<psd_length; i++){
	//	std::cout<<psd[0][i]<<data[0][i]<<std::endl;
	//	temp[i][0]=freqs[0][i];	
	//	temp[i][1]=std::real( data[0][i]/sqrt(psd[0][i]));
	//	temp[i][2]=std::imag( data[0][i]/sqrt(psd[0][i]));
	//}
	//write_file("whitened_GW150914_Hanford.csv",temp,psd_length,3);
	//deallocate_2D_array(temp,psd_length,3);

	int *data_length= (int*)malloc(sizeof(int)*num_detectors);
	data_length[0] =psd_length;
	data_length[1] =psd_length;
	//data_length[2] =psd_length;

	//#########################################################
	//mcmc options
	int dimension = 11;
	//double initial_pos[dimension]={2.5, sin(-.9),5.78,cos(3.1),.6,6,std::log(500),std::log(30), .24,.1,.1,.1,.1};
	//double initial_pos[dimension]={2.5, sin(-.9),.1,cos(.301),.6,6,std::log(500),std::log(30), .24,.1,.1,.1,.1};
	double initial_pos[dimension]={2.5, sin(-.9),1.78,cos(.01),.6,6,std::log(500),std::log(30), .24,.1,.1};
	//double initial_pos[dimension]={2.5, 0.05,.8,.8,.6,3,3,2, .24,.1,.1};
	double *seeding_var = NULL;
	int n_steps = 10000;
	//int chain_N=20 ;
	//int max_thermo=20 ;
	int chain_N=12 ;
	int max_thermo=12 ;
	int t0 = 1000;
	int nu = 100;
	std::string chain_alloc = "double";
	//std::string chain_alloc = "cold";
	//int swp_freq = 5;
	int swp_freq = 5;
	double chain_temps[chain_N];
	double c = 1.3;
	
	int Nmod = 0;
	int *bppe = NULL;
	int numThreads = 10;
	bool pool = true;
	//#########################################################
	//gw options
	//std::string generation_method = "dCS_IMRPhenomD_log";
	//std::string generation_method = "EdGB_IMRPhenomD_log";
	//std::string generation_method = "dCS_IMRPhenomD_root_alpha";
	//std::string generation_method = "IMRPhenomD";
	std::string generation_method = "IMRPhenomD";
	//std::string generation_method = "EdGB_IMRPhenomD_root_alpha";
	
	
	std::string chainfile = "testing/data/mcmc_output_uncorr_D.csv";
	std::string statfilename = "testing/data/mcmc_statistics_uncorr_D.txt";
	std::string checkfile = "testing/data/mcmc_checkpoint_uncorr_D.csv";

	int corr_threshold = 10;
	int corr_segments = 10;
	double corr_converge_thresh = 0.50;
	double corr_target_ac = .01;

	double **output;
	output = allocate_2D_array(n_steps, dimension );
	//double ***output;
	//output = allocate_3D_array(chain_N,n_steps, dimension );
	PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW(output, dimension, n_steps, chain_N, max_thermo, initial_pos,seeding_var,chain_temps, 
			swp_freq, t0, nu, corr_threshold, corr_segments, corr_converge_thresh, corr_target_ac,chain_alloc, test_lp_GW_D,numThreads, pool,show_progress,
			num_detectors, 
			data, psd,freqs, data_length,gps_time, detectors,Nmod, bppe,
			generation_method,statfilename,chainfile, "",checkfile);	
	//std::string checkfile2="testing/data/mcmc_checkpoint_uncorr_D2.csv";
	//continue_PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW(checkfile,output,  n_steps,  max_thermo, chain_temps, 
	//		swp_freq, t0, nu, corr_threshold, corr_segments, corr_converge_thresh, corr_target_ac,chain_alloc, test_lp_GW_D,numThreads, pool,show_progress,
	//		num_detectors, 
	//		data, psd,freqs, data_length,gps_time, detectors,Nmod, bppe,
	//		generation_method,statfilename,chainfile, "",checkfile2);	
	//PTMCMC_MH_dynamic_PT_alloc_GW(output, dimension, n_steps,chain_N,max_thermo, initial_pos, seeding_var, chain_temps, swp_freq, t0,nu,"half_ensemble",test_lp_GW_D,numThreads, pool, show_progress, num_detectors, data, psd, freqs, data_length, gps_time, detectors, Nmod, bppe, generation_method, statfilename, chainfile,  "",checkfile);


	//write_file(chainfile, output[0], n_steps, dimension);

	deallocate_2D_array(output, n_steps, dimension);
	//deallocate_3D_array(output, chain_N,n_steps, dimension);
	delete [] detectors;
	free(data_length);
	//free_LOSC_data(data, psd,freqs, num_detectors, length);
	deallocate_2D_array(psd,num_detectors, psd_length);
	deallocate_2D_array(freqs,num_detectors, psd_length);
	for(int i =0; i<num_detectors; i++)
		free(data[i]);
	free(data);
	delete [] detector_files;

}
void test55()
{
	int dimension = 2;
	double initial_pos[2]={1,-2.};
	double seeding_var[2] ={2,2} ;

	
	int N_steps = 10000;
	int chain_N= 8;
	int max_chain_N_thermo= 8;
	int t0= 1000;
	int nu= 100;
	//std::string chain_dist_method = "half_ensemble";
	std::string chain_dist_method = "double";
	//std::string chain_dist_method = "cold";
	double ***output;
	output = allocate_3D_array(chain_N, N_steps, dimension );
	//double *initial_pos_ptr = initial_pos;
	int swp_freq = 3;
	//double chain_temps[chain_N] ={1,2,3,10,12};
	double chain_temps[chain_N];
	//double chain_temps[chain_N] ={1};
	//std::string autocorrfile = "";
	//std::string chainfile = "testing/data/neil_mcmc_output_uncorr2.csv";
	std::string chainfile = "";
	std::string statfilename = "testing/data/neil_mcmc_statistics_uncorr2.txt";
	std::string checkpointfile = "testing/data/neil_mcmc_checkpoint_uncorr2.csv";
	std::string checkpointfile_start = "testing/data/neil_mcmc_checkpoint_uncorr.csv";
	
	int numThreads = 10;
	bool pool = true;

	int corr_threshold = 10;
	int corr_segments = 10;
	double corr_converge_thresh = 0.5;
	double corr_target_ac = .01;

	std::function<double(double* , int *,int, int,void *)> lp = [](double *param, int * status,int dim, int chain_id,void *parameters){ return test_lp_nts(param,dim,chain_id,parameters);};
	std::function<double(double* , int*,int, int,void *)> ll = [](double *param, int *status,int dim, int chain_id,void *parameters){ return log_neil_proj3_nts(param,dim,chain_id,parameters);};
	std::function<void(double* , int*,int, double**,int,void *)> f = [](double *param, int *status,int dim, double **fish,int chain_id,void *parameters){ fisher_neil_proj3(param,dim,fish,chain_id,parameters);};
	f=NULL;
	
	continue_PTMCMC_MH_dynamic_PT_alloc(checkpointfile_start ,output,  N_steps, max_chain_N_thermo, chain_temps, swp_freq, t0, nu,chain_dist_method,test_lp_nts, log_neil_proj3_nts,fisher_neil_proj3 ,(void **)NULL,numThreads, pool,show_progress, statfilename,"", "",checkpointfile );	
	std::cout<<"ENDED"<<std::endl;

	deallocate_3D_array( output,chain_N,N_steps, dimension);
}
void test54()
{
	int dimension = 2;
	double initial_pos[2]={5,-4.};
	double seeding_var[2] ={2,2} ;

	
	int N_steps = 510;
	int chain_N= 3;
	int max_chain_N_thermo= 3;
	int t0= 100;
	int nu= 50;
	//std::string chain_dist_method = "half_ensemble";
	std::string chain_dist_method = "double";
	//std::string chain_dist_method = "cold";
	//double **output;
	//output = allocate_2D_array(  N_steps, dimension );
	double **output;
	output = allocate_2D_array(N_steps, dimension );
	//double *initial_pos_ptr = initial_pos;
	int swp_freq = 5;
	//double chain_temps[chain_N] ={1,2,3,10,12};
	double chain_temps[chain_N];
	//double chain_temps[chain_N] ={1};
	//std::string autocorrfile = "";
	std::string chainfile = "testing/data/neil_mcmc_output_uncorr.csv";
	std::string statfilename = "testing/data/neil_mcmc_statistics_uncorr.txt";
	std::string checkpointfile = "testing/data/neil_mcmc_checkpoint_uncorr.csv";
	
	int numThreads = 10;
	bool pool = false;

	int corr_threshold = 5;
	int corr_segments = 20;
	double corr_converge_thresh = 1;
	double corr_target_ac = .01;

	std::function<double(double* , int *,int, int,void *)> lp = [](double *param, int * status,int dim, int chain_id,void *parameters){ return test_lp_nts(param,dim,chain_id,parameters);};
	std::function<double(double* , int*,int, int,void *)> ll = [](double *param, int *status,int dim, int chain_id,void *parameters){ return log_neil_proj3_nts(param,dim,chain_id,parameters);};
	std::function<void(double* , int*,int, double**,int,void *)> f = [](double *param, int *status,int dim, double **fish,int chain_id,void *parameters){ fisher_neil_proj3(param,dim,fish,chain_id,parameters);};
	//f=NULL;
	
	PTMCMC_MH_dynamic_PT_alloc_uncorrelated(output, dimension, N_steps, chain_N,max_chain_N_thermo, initial_pos,seeding_var,chain_temps, swp_freq, t0, nu,corr_threshold, corr_segments, corr_converge_thresh, corr_target_ac,chain_dist_method,test_lp_nts, log_neil_proj3_nts,NULL ,(void **)NULL,numThreads, pool,show_progress, statfilename,chainfile, "",checkpointfile );	
	//PTMCMC_MH_dynamic_PT_alloc_internal(output, dimension, N_steps, chain_N,max_chain_N_thermo, initial_pos,seeding_var,chain_temps, swp_freq, t0, nu,chain_dist_method,lp, ll,f,numThreads, pool,show_progress, statfilename,"", "",checkpointfile );	
	//PTMCMC_MH_dynamic_PT_alloc(output, dimension, N_steps, chain_N,max_chain_N_thermo, initial_pos,seeding_var,chain_temps, swp_freq, t0, nu,chain_dist_method,test_lp_nts, log_neil_proj3_nts,fisher_neil_proj3 ,(void **)NULL,numThreads, pool,show_progress, statfilename,"", "",checkpointfile );	
	std::cout<<"ENDED"<<std::endl;
	//std::cout<<"Chain temps: "<<std::endl;
	//for(int i =0; i<chain_N; i++){
	//	std::cout<<chain_temps[i]<<std::endl;
	//}


	deallocate_2D_array(output,  N_steps, dimension);
}
void test53()
{
	double T_obs = 4*T_year;
	double T_wait = 5*T_obs;
	double tolerance = 0.005;
	double times[2];
	double SNR_thresh = 8;

	gen_params params;
	params.mass1 = 36.e1;
	params.mass2 = 29.e1;
	//params.spin1[2] =-0.3;
	//params.spin2[2] = 0.8;
	params.spin1[2] = 0.0;
	params.spin2[2] = 0.0;
	params.Luminosity_Distance=4000;
	params.sky_average=true;
	params.incl_angle = 0;
	params.phi=0;params.theta=0;params.psi=0;
	params.equatorial_orientation=false;
	params.horizon_coord=true;
	params.tc = T_obs;
	params.phiRef = 1;
	params.f_ref = 20;
	params.NSflag1=false;
	params.NSflag2=false;
	params.shift_time = false;
	params.shift_phase = true;

	double f_lower = 1e-4;
	double f_upper = 1;


	gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
	params.sky_average=true;	

	std::cout<<"STARTING CALC"<<std::endl;
	clock_t start = clock();
	
	threshold_times_gsl(&params, "IMRPhenomD", T_obs, T_wait, f_lower,f_upper,"LISA_SADC_CONF",  SNR_thresh, times, tolerance,w,1000);
	std::cout<<"TIME: "<<(double)(clock() - start)/CLOCKS_PER_SEC<<std::endl;
	std::cout<<"TIMES: "<<times[0]/T_year<<" "<<times[1]/T_year<<std::endl;

	//delete [] freqs;
	//delete [] SN;
	gsl_integration_workspace_free(w);
}
void test52()
{
	int length = 1000;
	double freqs[length];
	for(int i = 0 ; i<length; i++){
		freqs[i]= 10 + i*1;
	}
	double psd[length];
	populate_noise(freqs, "AdLIGOMidHigh",psd, length);
	write_file("testing/data/AdLIGOMidHigh.csv",psd,length);
}
void test51()
{
	std::cout.precision(15);
	double rho_thresh = 8;
	double rho_opt = 9;
	double prob = p_triple_detector_interp(0.8001);
	std::cout<<fixed<<prob<<std::endl;
	prob = p_single_detector_interp(0.8001);
	std::cout<<fixed<<prob<<std::endl;
}
void test50()
{
	std::cout.precision(15);
	double prob ;
	int length =1000;
	double results[length];
	double omegas[length];
	double deltaomega = 1./(length);
	double init = deltaomega;
	#pragma omp parallel for
	for(int i = 0; i<length; i++){
		double omega = init + i * deltaomega;	
		omegas[i] = omega;
		//prob = p_single_detector(omega, 1e9);
		prob = p_single_detector_fit(omega);
		results[i]=prob;
		
	}
	for(int i = 0 ; i<length; i++){
		std::cout<<omegas[i]<<" "<<fixed<<results[i]<<std::endl;
	}
}
void test49()
{
	//std::string psd_file = "testing/data/GWTC1_GW150914_PSDs.dat.txt";
	//std::string psd_file = "testing/data/GWTC1_GW151226_PSDs.dat.txt";
	//std::string data_file = "testing/data/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//std::string psd_file = "/Users/sperkins/Downloads/LOSC_data/GW170729/GWTC1_GW170729_PSDs.dat.txt";
	std::string psd_file = "/Users/sperkins/Downloads/LOSC_data/GW150914/GWTC1_GW150914_PSDs.dat.txt";
	//int rows = 8032;
	//int cols = 3;
	int datalength = 131075;
	//double **psd = allocate_2D_array(rows, cols);
	//read_LOSC_PSD_file(psd_file, psd, rows, cols);
	//double data_start_time, duration, fs;
	int num_detectors = 2, psd_length = 8032, length;
	//int num_detectors = 3, psd_length = 4016, length;
	//double gps_time = 1135136350.6;//TESTING -- gw151226
	double gps_time = 1126259462.1;//TESTING -- gw150914
	//double gps_time = 1185389807.3;//TESTING -- gw170729
	std::string *detectors = new std::string[num_detectors];//(std::string*)malloc(sizeof(std::string)*50*num_detectors);
	detectors[0] = "Hanford";
	detectors[1] = "Livingston";
	//detectors[2] = "Virgo";
	std::string *detector_files = new std::string[num_detectors];
	//detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//detector_files[2] =  "testing/data/V-V1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1126259447-32.txt";
	detector_files[0] =  "/Users/sperkins/Downloads/LOSC_data/GW150914/H-H1_GWOSC_4KHZ_R1-1126259447-32.txt";
	detector_files[1] =  "/Users/sperkins/Downloads/LOSC_data/GW150914/L-L1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[2] =  "/home/sperkins/Downloads/LOSC_data/GW170729/V-V1_GWOSC_4KHZ_R1-1185389792-32.txt";
	//detector_files[0] =  "/Users/sperkins/Downloads/LOSC_data/GW170729/H-H1_GWOSC_4KHZ_R1-1185389792-32.txt";
	//detector_files[1] =  "/Users/sperkins/Downloads/LOSC_data/GW170729/L-L1_GWOSC_4KHZ_R1-1185389792-32.txt";
	//detector_files[2] =  "/Users/sperkins/Downloads/LOSC_data/GW170729/V-V1_GWOSC_4KHZ_R1-1185389792-32.txt";
 	//double trigger_time= 1135136350.6;
 	double trigger_time = gps_time;
	double **psd = allocate_2D_array(num_detectors,psd_length);
	double **freqs = allocate_2D_array(num_detectors,psd_length);
	std::complex<double> **data = (std::complex<double> **)malloc(sizeof(std::complex<double> *)*num_detectors);
	for(int i =0; i<num_detectors; i++)
		data[i] = (std::complex<double>*)malloc(sizeof(std::complex<double>)*psd_length);

	allocate_LOSC_data(detector_files, psd_file, num_detectors, psd_length, datalength, trigger_time, data, psd, freqs);

	int *data_length= (int*)malloc(sizeof(int)*num_detectors);
	data_length[0] =psd_length;
	data_length[1] =psd_length;
	//data_length[2] =psd_length;
	//#####################################################################	
	int dim = 11;
	int threads = 8;
	//std::string gen_meth= "IMRPhenomPv2";
	std::string gen_meth= "IMRPhenomD";
	//double init_pos[dim] = {6.,-1.,1.,.9,1.,3.,log(400),log(28),.24,0,0,.01,1.};
	double init_pos[dim] = {6.,-1.,1.,.9,1.,3.,log(400),log(28),.24,0,0};
	//double init_pos[dim] = {1.,-1.,1.,.9,1.,3.,log(400),log(28),.24,0,0,.01,1.};
	//double init_pos[dim] = {1.,-0.,1.,.9,3.,3.,log(2700),log(54),.24,.8,0,.01,1.};
	//double init_pos[dim] = {1.,-1.,1.,.9,1.,1.,log(400),log(28),.24,0,0};
	double *seeding_var = NULL;
	int swp_freq = 5;
	bool pool = true;
	//threads=1;
	bool show_prog=true;
	int Nmod= 0;
	int *bppe = NULL;
	std::string statf = "testing/data/mcmc_pv2_stat_data.txt";
	std::string chainf = "testing/data/mcmc_pv2_chain_data.csv";
	std::string checkf = "testing/data/mcmc_pv2_check_data.csv";
	std::string checkf2 = "testing/data/mcmc_pv2_check2_data.csv";
	int chain_n=8;
	int max_thermo=8;
	int steps =1e5;
	int stepsalloc =1e4;
	double temps[chain_n];
	int t0 = 1e3;
	int nu = 1e2;
	double c = 1.2;
	temps[0]=1;
	for(int i = 1 ; i<chain_n; i++){
		temps[i] = temps[i-1]*c;
	}
	double ***output=allocate_3D_array(chain_n, steps, dim);
	PTMCMC_MH_dynamic_PT_alloc_GW(output, dim, stepsalloc,chain_n,max_thermo, init_pos, seeding_var, temps, swp_freq, t0,nu,"half_ensemble",test_lp_GW_D,threads, pool, show_prog, num_detectors, data, psd, freqs, data_length, gps_time, detectors, Nmod, bppe, gen_meth, statf, chainf,  "",checkf);

	continue_PTMCMC_MH_GW(checkf,output, dim, steps,  swp_freq, test_lp_GW_D,threads, pool, show_prog, num_detectors, data, psd, freqs, data_length, gps_time, detectors, Nmod, bppe, gen_meth, statf, chainf, "", "testing/data/ll.csv",checkf2,true);
	//PTMCMC_MH_GW(output, dim, steps,chain_n, init_pos, seeding_var, temps, swp_freq, test_lp_GW_D,threads, pool, show_prog, num_detectors, data, psd, freqs, data_length, gps_time, detectors, Nmod, bppe, gen_meth, statf, chainf, "", "testing/data/ll.csv",checkf);
	
	deallocate_3D_array(output,chain_n, steps,dim);
	delete [] detectors;
	delete [] detector_files;
	deallocate_2D_array(psd, num_detectors, psd_length);
	deallocate_2D_array(freqs, num_detectors, psd_length);
	for(int i = 0 ; i<num_detectors; i++){
		free(data[i]);
	}
	free(data);
	free(data_length);
}
void test48()
{

	int lenA = 2e1;
	int lenB = 4e1;
	int lenD = 3e1;
	int lenE = 5e1;
	int **A = new int*[lenA];
	int **B = new int*[lenB];
	int **D = new int*[lenD];
	int **E = new int*[lenE];
	int *Avals = new int[lenA];
	int *Bvals = new int[lenB];
	int *Dvals = new int[lenD];
	int *Evals = new int[lenE];
	int lenC,lenG;
	for(int i = 0 ; i<lenA; i++){
		Avals[i] = 3*i;
		A[i] = &Avals[i];
	}
	for(int i = 0 ; i<lenB; i++){
		Bvals[i] = i;
		B[i] = &Bvals[i];
	}
	for(int i = 0 ; i<lenD; i++){
		Dvals[i] = 4*i;
		D[i] = &Dvals[i];
	}
	for(int i = 0 ; i<lenE; i++){
		Evals[i] = 6*i;
		E[i] = &Evals[i];
	}
	int **C = new int*[lenE];
	//int **G = new int*[lenE];
	int **G =new int*[lenE];
	clock_t start = clock();
	list_intersect_ptrs(A,lenA,B,lenB,C,&lenC);
	std::cout<<lenC<<std::endl;
	for(int i = 0 ; i<lenC ;i++){
		std::cout<<*C[i]<<std::endl;
	}
	std::cout<<std::endl;
	//G = C;
	for(int i = 0 ; i<lenC; i++){
		G[i] = C[i];	
		//std::cout<<*G[i]<<std::endl;
	}
	//std::cout<<G<<std::endl;
	//std::cout<<*G<<std::endl;
	//std::cout<<G[1]<<std::endl;
	//std::cout<<&(*(G[1]))<<std::endl;
	//std::cout<<&((*G)[1])<<std::endl;
	//std::cout<<*(G[1])<<std::endl;
	//std::cout<<(*G)[1]<<std::endl;
	lenG=lenC;
	//for(int i = 0 ; i<lenG; i++){
	//	std::cout<<G[i]<<std::endl;
	//}
	list_intersect_ptrs(G,lenG,D,lenD,C,&lenC);
	//list_intersect(D,lenD,*G,lenG,C,&lenC);
	for(int i = 0 ; i<lenC ;i++){
		std::cout<<*C[i]<<std::endl;
	}
	//std::cout<<std::endl;
	//G = C;
	for(int i = 0 ; i<lenC; i++){
		G[i] = C[i];	
		//std::cout<<G[i]<<" "<<C[i]<<std::endl;
		//std::cout<<(*G)[i]<<" "<<C[i]<<std::endl;
		//std::cout<<*G[i]<<std::endl;
	}
	//std::cout<<G<<std::endl;
	//std::cout<<*G<<std::endl;
	//std::cout<<G[1]<<std::endl;
	//std::cout<<&(*(G[1]))<<std::endl;
	//std::cout<<&((*G)[1])<<std::endl;
	//std::cout<<*(G[1])<<std::endl;
	//std::cout<<(*G)[1]<<std::endl;
	lenG=lenC;
	//for(int i = 0 ; i<lenG; i++){
	//	std::cout<<G[i]<<std::endl;
	//}
	list_intersect_ptrs(G,lenG,E,lenE,C,&lenC);
	////list_intersect(E,lenE,*G,lenG,C,&lenC);

	std::cout<<lenC<<" "<<(double)(clock()-start)/CLOCKS_PER_SEC<<std::endl;
	for(int i = 0 ; i<lenC ;i++){
		std::cout<<*C[i]<<std::endl;
	}
	////for(int i = 0 ; i<lenA ;i++){
	////	std::cout<<A[i]<<std::endl;
	////}
	//start = clock();
	//std::vector<int> v(lenB);                      // 0  0  0  0  0  0  0  0  0  0
	//std::vector<int>::iterator it;
	//it=std::set_intersection (A, A+lenA, B, B+lenB, v.begin());
	//v.resize(it-v.begin());
	//std::cout<<v.size()<<" "<<(double)(clock()-start)/CLOCKS_PER_SEC<<std::endl;
	//delete [] A;
	//delete [] B;
	//delete [] C;
}
void test47()
{
	int a[2] = {7,8};
	int b[2] = {10,11};
	int **c = new int *[2];
	int **d = new int *[2];
	c[0] = &a[0];
	c[1] = &a[1];
	d[0] = c[0];
	d[1] = c[1];
	c[0] = &b[0];	
	c[1] = &b[1];	
	c[1] = &(*d[0]);
	std::cout<<*d[0]<<std::endl;
	std::cout<<*d[1]<<std::endl;
	std::cout<<*c[0]<<std::endl;
	std::cout<<*c[1]<<std::endl;
	d[0] = c[0];
	d[1] = c[1];
	c[0] =&(* d[1]);
	c[1] = &a[0];
	std::cout<<*c[0]<<std::endl;
	std::cout<<*c[1]<<std::endl;
	

}
void test46()
{

	int lenA = 2e1;
	int lenB = 4e1;
	int lenD = 3e1;
	int lenE = 5e1;
	int *A = new int[lenA];
	int *B = new int[lenB];
	int *D = new int[lenD];
	int *E = new int[lenE];
	int lenC,lenG;
	for(int i = 0 ; i<lenA; i++){
		A[i] = 3*i;
	}
	for(int i = 0 ; i<lenB; i++){
		B[i] = i;
	}
	for(int i = 0 ; i<lenD; i++){
		D[i] = 4*i;
	}
	for(int i = 0 ; i<lenE; i++){
		E[i] = 6*i;
	}
	int **C = new int*[lenE];
	//int **G = new int*[lenE];
	int **G =new int*[lenE];
	clock_t start = clock();
	list_intersect(A,lenA,B,lenB,C,&lenC);
	for(int i = 0 ; i<lenC ;i++){
		std::cout<<*C[i]<<std::endl;
	}
	std::cout<<std::endl;
	G = C;
	//for(int i = 0 ; i<lenC; i++){
	//	G[i] = C[i];	
	//	//std::cout<<*G[i]<<std::endl;
	//}
	std::cout<<G<<std::endl;
	std::cout<<*G<<std::endl;
	std::cout<<G[1]<<std::endl;
	std::cout<<&(*(G[1]))<<std::endl;
	std::cout<<&((*G)[1])<<std::endl;
	std::cout<<*(G[1])<<std::endl;
	std::cout<<(*G)[1]<<std::endl;
	lenG=lenC;
	//for(int i = 0 ; i<lenG; i++){
	//	std::cout<<G[i]<<std::endl;
	//}
	list_intersect(*G,lenG,D,lenD,C,&lenC);
	//list_intersect(D,lenD,*G,lenG,C,&lenC);
	for(int i = 0 ; i<lenC ;i++){
		std::cout<<*C[i]<<std::endl;
	}
	std::cout<<std::endl;
	G = C;
	//for(int i = 0 ; i<lenC; i++){
	//	G[i] = C[i];	
	//	//std::cout<<G[i]<<" "<<C[i]<<std::endl;
	//	//std::cout<<(*G)[i]<<" "<<C[i]<<std::endl;
	//	//std::cout<<*G[i]<<std::endl;
	//}
	std::cout<<G<<std::endl;
	std::cout<<*G<<std::endl;
	std::cout<<G[1]<<std::endl;
	std::cout<<&(*(G[1]))<<std::endl;
	std::cout<<&((*G)[1])<<std::endl;
	std::cout<<*(G[1])<<std::endl;
	std::cout<<(*G)[1]<<std::endl;
	lenG=lenC;
	//for(int i = 0 ; i<lenG; i++){
	//	std::cout<<G[i]<<std::endl;
	//}
	list_intersect(*G,lenG,E,lenE,C,&lenC);
	//list_intersect(E,lenE,*G,lenG,C,&lenC);

	std::cout<<lenC<<" "<<(double)(clock()-start)/CLOCKS_PER_SEC<<std::endl;
	for(int i = 0 ; i<lenC ;i++){
		std::cout<<*C[i]<<std::endl;
	}
	//for(int i = 0 ; i<lenA ;i++){
	//	std::cout<<A[i]<<std::endl;
	//}
	start = clock();
	std::vector<int> v(lenB);                      // 0  0  0  0  0  0  0  0  0  0
	std::vector<int>::iterator it;
	it=std::set_intersection (A, A+lenA, B, B+lenB, v.begin());
	v.resize(it-v.begin());
	std::cout<<v.size()<<" "<<(double)(clock()-start)/CLOCKS_PER_SEC<<std::endl;
	delete [] A;
	delete [] B;
	//delete [] C;
}
void test45()
{
	int dim = 11;
	int threads = 10;
	gen_params params;
	params.mass1 = 35;
	params.mass2 = 30;
	double chirp = calculate_chirpmass(params.mass1,params.mass2);
	double eta = calculate_eta(params.mass1,params.mass2);
	std::cout<<"CHIRP: "<<chirp<<std::endl;
	std::cout<<"eta: "<<eta<<std::endl;
	//params.spin1[0]=0.5;
	//params.spin1[1]=0.7;
	params.spin1[2]=0.3;
	//params.spin2[0]=0.2;
	//params.spin2[1]=0.1;
	params.spin2[2]=0.2;
	params.phip = 0;
	params.chip = .4;
	params.Luminosity_Distance = 400;
	params.RA = 3.1;
	params.DEC = -.5;
	params.psi=1.;
	double tcref = 1;
	params.tc = tcref;
	double gps_time = 1185389807.3;//TESTING -- gw170729
	params.gmst = gps_to_GMST_radian(gps_time);
	params.shift_time = false;
	params.shift_phase = true;
	params.equatorial_orientation=false;
	params.phiRef = 2;
	params.f_ref = 20;
	params.incl_angle = 0.31;
	params.NSflag1=false;
	params.NSflag2=false;
	params.sky_average=false;
	int n_detect = 3;
	int length = 2000;
	std::string gen_meth="IMRPhenomD";
	int lengths[3] = {length,length,length};
	double fmax = 1024;
	double fmin = 10;
	double delta_f = (fmax-fmin)/(length-1);
	double **freqs= new double*[n_detect];
	double **psds= new double*[n_detect];
	std::string SN="aLIGO_analytic";
	std::string detectors[3] = {"Hanford","Livingston","Virgo"};
	std::complex<double> **data = new std::complex<double>*[n_detect];
	double phis[n_detect];
	double thetas[n_detect];
	celestial_horizon_transform(params.RA, params.DEC, gps_time, "Hanford", &phis[0], &thetas[0]);
	double snr=0;
	for(int i = 0 ; i<n_detect; i++){
		data[i] = new std::complex<double>[length];
		freqs[i] = new double[length];
		psds[i] = new double[length];
		for(int j = 0 ; j<length; j++){
			freqs[i][j] = fmin + j*delta_f;
		}	
		populate_noise(freqs[i],SN, psds[i],lengths[i]);
		for(int j = 0 ; j<length; j++){
			psds[i][j]=pow_int(psds[i][j],2);
		}
			
		celestial_horizon_transform(params.RA, params.DEC, gps_time, detectors[i], &phis[i], &thetas[i]);
		params.tc = tcref+ DTOA(thetas[0], thetas[i], "Hanford",detectors[i]); 
		//std::cout<<DTOA(thetas[0], thetas[i], "Hanford",detectors[i])<<std::endl;
		fourier_detector_response(freqs[i], lengths[i],data[i], detectors[i], gen_meth, &params, (double *)NULL);
		snr+=pow_int(calculate_snr(SN, data[i], freqs[i],lengths[i]),2);
	}
	std::cout<<"Network SNR: "<<sqrt(snr)<<std::endl;

	
	//double init_pos[dim] = {params.RA,params.DEC,params.psi,cos(params.incl_angle),params.phiRef,tcref,log(params.Luminosity_Distance),log(chirp),eta,params.spin1[2],params.spin2[2],params.chip,params.phip};
	double init_pos[dim] = {params.RA,params.DEC,params.psi,(params.incl_angle),params.phiRef,params.tc,log(params.Luminosity_Distance),log(chirp),eta,params.spin1[2],params.spin2[2]};
	double *seeding_var = NULL;
	int swp_freq = 5;
	bool pool = true;
	//threads=1;
	bool show_prog=true;
	int Nmod= 0;
	int *bppe = NULL;
	std::string statf = "testing/data/mcmc_pv2_stat.txt";
	std::string chainf = "testing/data/mcmc_pv2_chain.csv";
	std::string checkf = "testing/data/mcmc_pv2_check.csv";
	std::string checkf2 = "testing/data/mcmc_pv2_check2.csv";
	int chain_n=8;
	int max_thermo=8;
	int steps =1e3;
	int stepsalloc =1e4;
	double temps[chain_n];
	int t0 = 1e3;
	int nu = 1e2;
	double c = 1.2;
	temps[0]=1;
	for(int i = 1 ; i<chain_n; i++){
		temps[i] = temps[i-1]*c;
	}
	double ***output=allocate_3D_array(chain_n, steps, dim);
	load_temps_checkpoint_file(checkf,temps, chain_n);
	std::string log_file = "testing/data/test_log.csv";
	//PTMCMC_MH_dynamic_PT_alloc_GW(output, dim, stepsalloc,chain_n,max_thermo, init_pos, seeding_var, temps, swp_freq, t0,nu,"half_ensemble",test_lp_GW_D,threads, pool, show_prog, n_detect, data, psds, freqs, lengths, gps_time, detectors, Nmod, bppe, gen_meth, statf, chainf,  "",checkf);

	//continue_PTMCMC_MH_GW(checkf,output, dim, steps,  swp_freq, test_lp_GW_D,threads, pool, show_prog, n_detect, data, psds, freqs, lengths, gps_time, detectors, Nmod, bppe, gen_meth, statf, chainf, "", "",checkf2);
	PTMCMC_MH_GW(output, dim, steps,chain_n, init_pos, seeding_var, temps, swp_freq, test_lp_GW_D,threads, pool, show_prog, n_detect, data, psds, freqs, lengths, gps_time, detectors, Nmod, bppe, gen_meth, statf, chainf, "", log_file,checkf);
	
	deallocate_3D_array(output,chain_n, steps,dim);
	for(int i = 0 ; i<n_detect; i++){
		delete [] data[i];
		delete [] freqs[i];
		delete [] psds[i];
	}
	delete [] data;
	delete [] freqs;
	delete [] psds;
}
void test_prob(double *prob, void *param, int d, int thread_id)
{
	double *internal_param = (double *)param;
	double sigma1 = 2;
	double sigma2 = 3;
	double sigma3 = 2;
	double sigma4 = 2;
	double gauss1 = 1./(sqrt(M_PI*sigma1*sigma1)) * exp(-(internal_param[0]-2)*(internal_param[0]-2)/(2.*sigma1*sigma1));
	double gauss2 = 1./(sqrt(M_PI*sigma2*sigma2)) * exp(-internal_param[1]*internal_param[1]/(2.*sigma2*sigma2));
	double gauss3 = 1./(sqrt(M_PI*sigma3*sigma3)) * exp(-internal_param[2]*internal_param[2]/(2.*sigma3*sigma3));
	double gauss4 = 1./(sqrt(M_PI*sigma4*sigma4)) * exp(-internal_param[0]*internal_param[2]/(2.*sigma4*sigma4));
	*prob = gauss1*gauss2*gauss3*gauss4;
}
void test43()
{
	gen_params params;
	double RA=M_PI, DEC=0.0;
	double gps_time = 1185389807.3;//TESTING -- gw170729
	params.mass1 = 35.1;
	params.mass2 = 29.5;
	//params.spin1[0] = .01;
	//params.spin1[1] = 0;
	//params.spin2[0] = 0;
	//params.spin2[1] = 0;
	params.spin1[2] = .31;
	params.spin2[2] = .39;
	params.chip = .1;
	params.phip = .1;
	params.sky_average = false;
	params.NSflag1 = false;
	params.NSflag2 = false;
	params.shift_time = false;
	params.shift_phase = false;

	params.equatorial_orientation=true;
	//params.theta_l = M_PI/2. - DEC+.01;
	params.theta_l = M_PI/2.;
	params.phi_l = .01;
	//params.theta_l = 0;
	//params.phi_l = 0;
	params.RA = RA;
	params.DEC = DEC;
	//params.theta = 0;
	//params.phi = 0;
	params.Luminosity_Distance = 400;
	params.phiRef = M_PI/2.;
	params.f_ref=20;
	params.tc = 0;
	params.gmst=gps_to_GMST(gps_time);
	//std::string gen_method = "MCMC_dCS_IMRPhenomD_log_Full";
	std::string gen_method = "IMRPhenomPv2";
	
	int length = 1000;
	double freqs[length];
	//double fhigh = 1.6e-2;
	//double flow = 1.3e-2;
	double fhigh=1, flow;
	Tbm_to_freq(&params, gen_method, 5.1*T_year, &flow, .01);
	//flow = 1.5e-1;
	std::cout<<"Low frequency cutoff "<<flow<<std::endl;
	double fstep = (fhigh-flow)/length;
	for(int i =0; i<length; i++){
		freqs[i] = flow + i*fstep;
	}
	std::complex<double> response[length];
	std::complex<double> hp[length];
	std::complex<double> hc[length];
	double tc=2;
	params.tc = tc;
	transform_orientation_coords(&params, gen_method, "");
	std::cout<<"Inclination angle "<<params.incl_angle<<std::endl;
	
	std::cout<<params.incl_angle<<" "<<params.psi<<std::endl;
	fourier_waveform(freqs, length, hp,hc,  gen_method,&params);
	double snr = calculate_snr("LISA_SADC_CONF",hp, freqs, length);
	std::cout<<"SA and DC: "<<snr<<std::endl;

	params.sky_average=false;
	double times[length];
	time_phase_corrected_autodiff(times, length, freqs, &params, gen_method, false,NULL);	
	std::cout<<"Integration time: "<<(times[0]-times[length-1])/T_year<<std::endl;
	fourier_detector_response(freqs, length, response, "LISA", gen_method,&params,times);
	for(int i = 0 ; i<length; i++){
		//std::cout<<response[i]<<std::endl;
	}
	snr = calculate_snr("LISA_CONF",response, freqs, length);
	std::cout<<"Not SA, dual channel "<<sqrt(2.)*snr<<std::endl;

}
void test42()
{
	gsl_rng_env_setup();
	const gsl_rng_type *T = gsl_rng_default;
	gsl_rng * r = gsl_rng_alloc(T);
	gsl_rng_set(r, 1);
	
	int tests=2000;
	double sphvec[3];
	double cartvec[3];
	for(int i = 0 ; i<tests; i++){
		for(int j = 0 ; j<3; j++){
			cartvec[j] = gsl_rng_uniform(r);		
		}
		transform_cart_sph(cartvec,sphvec);	
		std::cout<<cartvec[0] - sphvec[0]*sin(sphvec[1]	)*cos(sphvec[2])<<" ";
		std::cout<<cartvec[1] - sphvec[0]*sin(sphvec[1]	)*sin(sphvec[2])<<" ";
		std::cout<<cartvec[2] - sphvec[0]*cos(sphvec[1]	)<<" ";
		std::cout<<std::endl;	
	}
	gsl_rng_free(r);
}
void test41()
{
	gen_params params;
	double RA=M_PI/2., DEC=.0;
	double gps_time = 1185389807.3;//TESTING -- gw170729
	params.mass1 = 80;
	params.mass2 = 40;
	//params.spin1[0] = .01;
	//params.spin1[1] = 0;
	//params.spin2[0] = 0;
	//params.spin2[1] = 0;
	params.spin1[2] = .8;
	params.spin2[2] = .3;
	params.chip = .8;
	params.phip = .1;
	params.sky_average = false;
	params.NSflag1 = false;
	params.NSflag2 = false;
	params.shift_time = false;
	params.shift_phase = false;

	params.equatorial_orientation=true;
	params.theta_l = M_PI/3.;
	params.phi_l = M_PI/3.;
	//params.theta_l = 0;
	//params.phi_l = 0;
	params.RA = RA;
	params.DEC = DEC;
	//params.theta = 0;
	//params.phi = 0;
	params.Luminosity_Distance = 400;
	params.phiRef = M_PI/2.;
	params.f_ref=20;
	params.tc = 0;
	params.gmst=gps_to_GMST(gps_time);
	//std::string gen_method = "MCMC_dCS_IMRPhenomD_log_Full";
	std::string gen_method = "IMRPhenomPv2";
	std::string detector = "Hanford";
	
	int length = 100;
	double freqs[length];
	double fhigh = 1024;
	double flow = 10;
	double fstep = (fhigh-flow)/length;
	for(int i =0; i<length; i++){
		freqs[i] = flow + i*fstep;
	}
	std::complex<double> response[length];
	std::complex<double> response2[length];
	std::complex<double> response3[length];
	double tc=2;
	params.tc = tc;

	fourier_detector_response(freqs, length, response, detector, gen_method,&params);

	params.equatorial_orientation=false;
	params.incl_angle =  0.326558168896834;
	params.psi = 1.44304830664865;
	fourier_detector_response(freqs, length, response2, detector, gen_method,&params);

	//celestial_horizon_transform(RA, DEC, gps_time, "Hanford", &params.phi, &params.theta);
	//params.horizon_coord=true;
	//params.psi = M_PI;
	//fourier_detector_response(freqs, length, response3, detector, gen_method,&params);
	double temptheta = params.theta;

	for(int i = 0 ; i<length; i++){
		std::cout<<(response[i]-response2[i])/response[i]<<std::endl;
	}
}
void test40()
{
	double phi_eq = 285.*M_PI/180.;
	double theta_eq = M_PI/2. - 20.*M_PI/180.;
	double theta_ecl;
	double phi_ecl;
	ecl_from_eq(theta_eq, phi_eq, &theta_ecl, &phi_ecl);
	std::cout<<M_PI/2. - theta_ecl<<" "<<phi_ecl<<std::endl;
}
void test39()
{
	
	int length = 100;
	double phase_in[length];
	double phase_out[length];
	double phase_corr[length];
	double deltap  = .1;
	for (int i = 0 ;i < length; i++){
		phase_in[i]= -i*deltap;
		phase_out[i]= std::atan2(std::sin(phase_in[i]), std::cos(phase_in[i]));
	}
	unwrap_array(phase_out, phase_corr, length);
	//phase_corr[0]=phase_out[0];
	//for(int i = 1 ; i<length; i++){
	//	phase_corr[i] = unwrap( phase_corr[i-1], phase_out[i]);
	//}
	for (int i = 0 ;i < length; i++){
		std::cout<<phase_in[i]<<" "<<phase_out[i]<<" "<<phase_corr[i]<<std::endl;
	}
}
void test38()
{
	//std::string psd_file = "testing/data/GWTC1_GW150914_PSDs.dat.txt";
	//std::string psd_file = "testing/data/GWTC1_GW151226_PSDs.dat.txt";
	//std::string data_file = "testing/data/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	std::string psd_file = "/Users/sperkins/Downloads/LOSC_data/GW170729/GWTC1_GW170729_PSDs.dat.txt";
	//int rows = 8032;
	//int cols = 3;
	int datalength = 131075;
	//double **psd = allocate_2D_array(rows, cols);
	//read_LOSC_PSD_file(psd_file, psd, rows, cols);
	//double data_start_time, duration, fs;
	//int num_detectors = 2, psd_length = 8032, length;
	int num_detectors = 3, psd_length = 4016, length;
	//double gps_time = 1135136350.6;//TESTING -- gw151226
	//double gps_time = 1126259462;//TESTING -- gw150914
	double gps_time = 1185389807.3;//TESTING -- gw170729
	std::string *detectors = new std::string[num_detectors];//(std::string*)malloc(sizeof(std::string)*50*num_detectors);
	detectors[0] = "Hanford";
	detectors[1] = "Livingston";
	detectors[2] = "Virgo";
	std::string *detector_files = new std::string[num_detectors];
	//detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//detector_files[2] =  "testing/data/V-V1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1126259447-32.txt";
	detector_files[0] =  "/Users/sperkins/Downloads/LOSC_data/GW170729/H-H1_GWOSC_4KHZ_R1-1185389792-32.txt";
	detector_files[1] =  "/Users/sperkins/Downloads/LOSC_data/GW170729/L-L1_GWOSC_4KHZ_R1-1185389792-32.txt";
	detector_files[2] =  "/Users/sperkins/Downloads/LOSC_data/GW170729/V-V1_GWOSC_4KHZ_R1-1185389792-32.txt";
 	//double trigger_time= 1135136350.6;
 	double trigger_time = gps_time;
	double **psd = allocate_2D_array(num_detectors,psd_length);
	double **freqs = allocate_2D_array(num_detectors,psd_length);
	std::complex<double> **data = (std::complex<double> **)malloc(sizeof(std::complex<double> *)*num_detectors);
	for(int i =0; i<num_detectors; i++)
		data[i] = (std::complex<double>*)malloc(sizeof(std::complex<double>)*psd_length);

	allocate_LOSC_data(detector_files, psd_file, num_detectors, psd_length, datalength, trigger_time, data, psd, freqs);

	int *data_length= (int*)malloc(sizeof(int)*num_detectors);
	data_length[0] =psd_length;
	data_length[1] =psd_length;
	data_length[2] =psd_length;

	//#########################################################
	//mcmc options
	int dimension = 4;
	double initial_pos[dimension]={std::log(50),.23,.01,.01};
	double *seeding_var = NULL;
	int n_steps = 50000;
	int chain_N=12 ;
	double ***output;
	output = allocate_3D_array( chain_N, n_steps, dimension );
	int swp_freq = 5;
	double chain_temps[chain_N];
	double c = 1.1;
	chain_temps[0]=1.;
	for(int i =1; i < chain_N;  i ++)
		chain_temps[i] = c*chain_temps[i-1];
	
	int Nmod = 0;
	//int *bppe = new int[Nmod];
	//bppe[0] = -1;
	int *bppe;
	int numThreads = 10;
	bool pool = true;
	//#########################################################
	//gw options
	std::string generation_method = "IMRPhenomD";
	
	std::string autocorrfile = "testing/data/auto_corr_mcmc_reduceddim.csv";
	std::string chainfile = "testing/data/mcmc_output_reduceddim.csv";
	std::string statfilename = "testing/data/mcmc_statistics_reduceddim.txt";
	std::string checkfile = "testing/data/mcmc_checkpoint_reduceddim.csv";
	//std::string autocorrfile = "testing/data/auto_corr_mcmc_EdGB.csv";
	//std::string chainfile = "testing/data/mcmc_output_EdGB.csv";
	//std::string statfilename = "testing/data/mcmc_statistics_EdGB.txt";
	//std::string checkfile = "testing/data/mcmc_checkpoint_EdGB.csv";

	PTMCMC_MH_GW(output, dimension, n_steps, chain_N, initial_pos,seeding_var,chain_temps, 
			swp_freq, test_lp_reduceddim,numThreads, pool,show_progress,
			num_detectors, 
			data, psd,freqs, data_length,gps_time, detectors,Nmod, bppe,
			generation_method,statfilename,"","", "",checkfile);	

	//double **output_transform=(double **)malloc(sizeof(double*)*n_steps);
	//for (int j =0; j<n_steps; j++)
	//	output_transform[j] = (double *)malloc(sizeof(double)*dimension);

	//for(int j = 0; j<n_steps;j++){
	//	output_transform[j][0]=output[0][j][0];
	//	output_transform[j][1]=output[0][j][1];
	//	output_transform[j][2]=output[0][j][2];
	//	output_transform[j][3]=std::exp(output[0][j][3]);
	//	output_transform[j][4]=std::exp(output[0][j][4]);
	//	output_transform[j][5]=output[0][j][5];
	//	output_transform[j][6]=output[0][j][6];
	//	output_transform[j][7]=output[0][j][7];
	//	//output_transform[j][8]=std::exp(output[0][j][8]);
	//	output_transform[j][8]=output[0][j][8];
	//}
	write_file(chainfile, output[0], n_steps, dimension);
	//std::string chain2 = "testing/data/mcmc_output_dCS_2.csv";
	//write_file(chain2, output[4], n_steps, dimension);
	//output hottest chain too
	//chainfile = "testing/data/mcmc_output_EdGB_hot.csv";
	//chainfile = "testing/data/mcmc_output_dCS_hot.csv";
	//for(int j = 0; j<n_steps;j++){
	//	output_transform[j][0]=output[chain_N-1][j][0];
	//	output_transform[j][1]=output[chain_N-1][j][1];
	//	output_transform[j][2]=output[chain_N-1][j][2];
	//	output_transform[j][3]=std::exp(output[chain_N-1][j][3]);
	//	output_transform[j][4]=std::exp(output[chain_N-1][j][4]);
	//	output_transform[j][5]=output[chain_N-1][j][5];
	//	output_transform[j][6]=output[chain_N-1][j][6];
	//	output_transform[j][7]=output[chain_N-1][j][7];
	//	//output_transform[j][8]=std::exp(output[chain_N-1][j][8]);
	//	output_transform[j][8]=output[chain_N-1][j][8];
	//}
	//write_file(chainfile, output_transform, n_steps, dimension);

	deallocate_3D_array(output, chain_N, n_steps, dimension);
	//for(int i =0; i< n_steps; i++){
	//	free(output_transform[i]);
	//}
	//free(output_transform);
	delete [] detectors;
	free(data_length);
	//free_LOSC_data(data, psd,freqs, num_detectors, length);
	deallocate_2D_array(psd,num_detectors, psd_length);
	deallocate_2D_array(freqs,num_detectors, psd_length);
	for(int i =0; i<num_detectors; i++)
		free(data[i]);
	free(data);
	delete [] detector_files;
	delete [] bppe;
	//deallocate_2D_array(psd, rows, cols);
	//free(data);
	
	//int length = 2000;
	//double x[length];
	//double y[length];
	//double xlim = 10.;
	//double xstart=0;
	//double xstep = (xlim-xstart)/length;
	//for (int i =0; i<length; i++){
	//	x[i]= i*xstep;
	//	y[i]= x[i]*x[i];
	//}
	//double sum = simpsons_sum(xstep, length, y);
	//std::cout<<sum<<std::endl;

}
void test37()
{
	gen_params params;
	int length = 16000;
	params.mass1 = 200;
	params.mass2 = 50;
	string method= "IMRPhenomPv2";
	//string method= "IMRPhenomD";
	//string method= "ppE_IMRPhenomPv2_IMR";
	//string method= "IMRPhenomD";
	//string method= "ppE_IMRPhenomD_Inspiral";
	double amp[length];
	double phaseout[length];
	complex<double> waveformout_plus[length];
	complex<double> waveformout_cross[length];
	double waveformout_plus_real[length];
	double waveformout_cross_real[length];
	double waveformout_plus_imag[length];
	double waveformout_cross_imag[length];
	complex<double> response[length];
	params.spin1[0] = .0;
	params.spin1[1] = .01;
	params.spin1[2] = -.2;
	params.spin2[0] = .0;
	params.spin2[1] = 0.1;
	params.spin2[2] = .9;
	double *spin1  = params.spin1;
	double *spin2= params.spin2;
	//params.phic = 2.0;
	//params.tc = 8.0;
	params.Luminosity_Distance = 800.;
	params.betappe = new double[1] ;
	params.betappe[0]=1.;
	params.bppe  =new int[1];
	params.bppe[0] = -1;
	params.Nmod = 1;
	params.NSflag1 = false;
	params.NSflag2 = false;
	params.sky_average = false;
	//params.phi = M_PI/3.;
	//params.theta = M_PI/3;
	params.RA = 2;
	params.DEC = 1;
	params.gmst = 10;
	params.phi = 0;
	params.theta = 0;
	params.incl_angle = 0;
	params.phiRef = 10;
	params.f_ref = 15;
	
	double freq[length];
	for(int i=0;i<length;i++)
		freq[i]=10.+i*1e-1;

	clock_t  start, end;
	start = clock(); 
	//fourier_waveform(freq, length, waveformout_plus,waveformout_cross,method,&params);
	int loops = 1000;
	for(int i = 0 ; i<loops; i++)
		fourier_waveform(freq, length, waveformout_plus_real,waveformout_plus_imag,waveformout_cross_real,waveformout_cross_imag,method,&params);
	//fourier_detector_response_equatorial(freq, length, response,"Hanford",method,&params);
	end=clock();
	//for(int i = 0 ;i<length; i++){
	//	std::cout<<waveformout_plus[i]<<" "<<waveformout_cross[i]<<std::endl;
	//}
	cout<<"TIMING waveform: "<<(double)(end-start)/(CLOCKS_PER_SEC*loops)<<endl;
	delete [] params.betappe;
	delete [] params.bppe;
	
}

//Proof of concept for PTMCMC (non-RJ) param_status array
void test36()
{
	int length =10;
	int width =12;	
	int depth =5;	
	//double **A = allocate_2D_array(length, width);
	double ***A = (double ***)malloc(sizeof(double **)*length);
	A[0] = (double**)malloc(sizeof(double*)*width);
	A[0][0] = (double*)malloc(sizeof(double)*depth);
	for(int i = 0 ;i<depth; i++){
		
		A[0][0][i]=i;
	}
	for(int i = 1 ;i<width; i++){
		A[0][i] = A[0][0];
	}
	for(int i = 0 ;i<depth; i++){
		std::cout<<A[0][0][i]<<std::endl;	
	}
	A[1]=A[0];
	A[2]=A[0];
	for(int i = 0 ;i<depth; i++){
		std::cout<<A[1][0][i]<<std::endl;	
		std::cout<<A[2][0][i]<<std::endl;	
	}
	free(A[0][0]);	
	free(A[0]);	
	free(A);
	//deallocate_2D_array(A, length,width);
}
void test35()
{
	//std::string psd_file = "/home/sperkins/Downloads/LOSC_data/GW151226/GWTC1_GW151226_PSDs.dat.txt";
	//std::string psd_file = "/Users/sperkins/Downloads/LOSC_data/GW151226/GWTC1_GW151226_PSDs.dat.txt";
	//std::string psd_file = "/Users/sperkins/Downloads/LOSC_data/GW150914/GWTC1_GW150914_PSDs.dat.txt";
	std::string psd_file = "/Users/sperkins/Downloads/LOSC_data/GW150914/GWTC1_GW150914_PSDs.dat.txt";
	//std::string psd_file = "/home/sperkins/Downloads/LOSC_data/GW170729/GWTC1_GW170729_PSDs.dat.txt";
	int datalength = 131075;
	int num_detectors = 2, psd_length = 8032, length;
	//int num_detectors = 3, psd_length = 4016, length;
	//double gps_time = 1135136350.6;//TESTING -- gw151226
	double gps_time = 1126259462.4;//TESTING -- gw150914
	//double gps_time = 1185389807.3;//TESTING -- gw170729
	std::string *detectors = new std::string[num_detectors];//(std::string*)malloc(sizeof(std::string)*50*num_detectors);
	detectors[0] = "Hanford";
	detectors[1] = "Livingston";
	//detectors[2] = "Virgo";
	std::string *detector_files = new std::string[num_detectors];
	//detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[0] =  "/home/sperkins/Downloads/LOSC_data/GW151226/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//detector_files[1] =  "/home/sperkins/Downloads/LOSC_data/GW151226/L-L1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//detector_files[0] =  "/Users/sperkins/Downloads/LOSC_data/GW151226/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//detector_files[1] =  "/Users/sperkins/Downloads/LOSC_data/GW151226/L-L1_GWOSC_4KHZ_R1-1135136335-32.txt";
	detector_files[0] =  "/Users/sperkins/Downloads/LOSC_data/GW150914/H-H1_GWOSC_4KHZ_R1-1126259447-32.txt";
	detector_files[1] =  "/Users/sperkins/Downloads/LOSC_data/GW150914/L-L1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[0] =  "/home/sperkins/Downloads/LOSC_data/GW150914/H-H1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[1] =  "/home/sperkins/Downloads/LOSC_data/GW150914/L-L1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1185389792-32.txt";
	//detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1185389792-32.txt";
	//detector_files[2] =  "testing/data/V-V1_GWOSC_4KHZ_R1-1185389792-32.txt";
	//detector_files[0] =  "/home/sperkins/Downloads/LOSC_data/GW170729/H-H1_GWOSC_4KHZ_R1-1185389792-32.txt";
	//detector_files[1] =  "/home/sperkins/Downloads/LOSC_data/GW170729/L-L1_GWOSC_4KHZ_R1-1185389792-32.txt";
	//detector_files[2] =  "/home/sperkins/Downloads/LOSC_data/GW170729/V-V1_GWOSC_4KHZ_R1-1185389792-32.txt";
 	//double trigger_time= 1135136350.6;
 	double trigger_time = gps_time;
	double **psd = allocate_2D_array(num_detectors,psd_length);
	double **freqs = allocate_2D_array(num_detectors,psd_length);
	std::complex<double> **data = (std::complex<double> **)malloc(sizeof(std::complex<double> *)*num_detectors);
	for(int i =0; i<num_detectors; i++)
		data[i] = (std::complex<double>*)malloc(sizeof(std::complex<double>)*psd_length);

	allocate_LOSC_data(detector_files, psd_file, num_detectors, psd_length, datalength, trigger_time, data, psd, freqs);


	int *data_length= (int*)malloc(sizeof(int)*num_detectors);
	data_length[0] =psd_length;
	data_length[1] =psd_length;
	//data_length[2] =psd_length;

	//#########################################################
	//mcmc options
	int dimension = 14;
	int n_steps = 1000;
	int chain_N=24 ;
	double ***output;
	output = allocate_3D_array( chain_N, n_steps, dimension );
	int swp_freq = 5;
	double chain_temps[chain_N];
	
	int Nmod = 0;
	int *bppe = NULL;
	int numThreads = 15;
	bool pool = true;
	//#########################################################
	//gw options
	std::string generation_method = "IMRPhenomPv2";
	
	
	std::string chainfile_base = "testing/data/mcmc_output_Pv2_dpt";
	std::string statfilename = "testing/data/mcmc_statistics_Pv2_dpt.txt";
	std::string checkfile = "testing/data/mcmc_checkpoint_Pv2_dpt.csv";
	std::string start_checkfile = "testing/data/mcmc_checkpoint_Pv2_dpt_alloc.csv";
	double *temps = (double *)malloc(sizeof(double)*chain_N);
	load_temps_checkpoint_file(start_checkfile,temps, chain_N);

	//data[1] = data[2];
	//psd[1] = psd[2];
	//freqs[1] = freqs[2];
	//detectors[1] ="Virgo";

	continue_PTMCMC_MH_GW(start_checkfile, output, dimension, n_steps,   
			swp_freq,  test_lp_GW_Pv2,numThreads, pool,show_progress,
			num_detectors, 
			data, psd,freqs, data_length,gps_time, detectors,Nmod, bppe,
			generation_method,statfilename,"","", "",checkfile,true);	


	std::string ind_chain_file;
	int criteria = 0 ;
	int file_ct = 0 ;
	while(criteria ==0){
		ind_chain_file = chainfile_base+std::to_string(file_ct)+".csv";
		criteria = std::remove(ind_chain_file.c_str());
		file_ct++;
	}
	file_ct = 0;
	for(int i =0; i < chain_N; i++){
		std::cout<<temps[i]<<std::endl;
		if(temps[i]==1){
			ind_chain_file = chainfile_base + to_string(file_ct) + ".csv";
			write_file(ind_chain_file , output[i], n_steps, dimension);
			file_ct++;
		}
	}

	deallocate_3D_array(output, chain_N, n_steps, dimension);
	delete [] detectors;
	free(data_length);
	deallocate_2D_array(psd,num_detectors, psd_length);
	deallocate_2D_array(freqs,num_detectors, psd_length);
	for(int i =0; i<num_detectors; i++)
		free(data[i]);
	free(data);
	delete [] detector_files;
	delete [] temps;
	
}
void test34()
{
	//std::string psd_file = "testing/data/GWTC1_GW150914_PSDs.dat.txt";
	//std::string psd_file = "testing/data/GWTC1_GW151226_PSDs.dat.txt";
	//std::string data_file = "testing/data/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//std::string psd_file = "testing/data/GWTC1_GW170729_PSDs.dat.txt";
	//std::string psd_file = "/home/sperkins/Downloads/LOSC_data/GW151226/GWTC1_GW151226_PSDs.dat.txt";
	//std::string psd_file = "/Users/sperkins/Downloads/LOSC_data/GW151226/GWTC1_GW151226_PSDs.dat.txt";
	//std::string psd_file = "/Users/sperkins/Downloads/LOSC_data/GW150914/GWTC1_GW150914_PSDs.dat.txt";
	std::string psd_file = "/home/sperkins/Downloads/LOSC_data/GW150914/GWTC1_GW150914_PSDs.dat.txt";
	//std::string psd_file = "/home/sperkins/Downloads/LOSC_data/GW170729/GWTC1_GW170729_PSDs.dat.txt";
	//std::string psd_file = "/Users/sperkins/Downloads/LOSC_data/GW170729/GWTC1_GW170729_PSDs.dat.txt";
	//int rows = 8032;
	//int cols = 3;
	int datalength = 131075;
	//double **psd = allocate_2D_array(rows, cols);
	//read_LOSC_PSD_file(psd_file, psd, rows, cols);
	//double data_start_time, duration, fs;
	int num_detectors = 2, psd_length = 8032, length;
	//int num_detectors = 3, psd_length = 4016, length;
	//int num_detectors = 3, psd_length = 4016, length;
	//int num_detectors = 2, psd_length = 4016, length;
	//double gps_time = 1135136350.6;//TESTING -- gw151226
	double gps_time = 1126259462.4;//TESTING -- gw150914
	//double gps_time = 1185389807.3;//TESTING -- gw170729
	//double gps_time = 1185389807.3;//TESTING -- gw170729
	std::string *detectors = new std::string[num_detectors];//(std::string*)malloc(sizeof(std::string)*50*num_detectors);
	detectors[0] = "Hanford";
	detectors[1] = "Livingston";
	//detectors[2] = "Virgo";
	std::string *detector_files = new std::string[num_detectors];
	//detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[0] =  "/home/sperkins/Downloads/LOSC_data/GW151226/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//detector_files[1] =  "/home/sperkins/Downloads/LOSC_data/GW151226/L-L1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//detector_files[0] =  "/Users/sperkins/Downloads/LOSC_data/GW151226/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//detector_files[1] =  "/Users/sperkins/Downloads/LOSC_data/GW151226/L-L1_GWOSC_4KHZ_R1-1135136335-32.txt";
	detector_files[0] =  "/home/sperkins/Downloads/LOSC_data/GW150914/H-H1_GWOSC_4KHZ_R1-1126259447-32.txt";
	detector_files[1] =  "/home/sperkins/Downloads/LOSC_data/GW150914/L-L1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[0] =  "/home/sperkins/Downloads/LOSC_data/GW150914/H-H1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[1] =  "/home/sperkins/Downloads/LOSC_data/GW150914/L-L1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[0] =  "/home/sperkins/Downloads/LOSC_data/GW170729/H-H1_GWOSC_4KHZ_R1-1185389792-32.txt";
	//detector_files[1] =  "/home/sperkins/Downloads/LOSC_data/GW170729/L-L1_GWOSC_4KHZ_R1-1185389792-32.txt";
	//detector_files[2] =  "/home/sperkins/Downloads/LOSC_data/GW170729/V-V1_GWOSC_4KHZ_R1-1185389792-32.txt";
	//detector_files[0] =  "/Users/sperkins/Downloads/LOSC_data/GW170729/H-H1_GWOSC_4KHZ_R1-1185389792-32.txt";
	//detector_files[1] =  "/Users/sperkins/Downloads/LOSC_data/GW170729/L-L1_GWOSC_4KHZ_R1-1185389792-32.txt";
	//detector_files[2] =  "/Users/sperkins/Downloads/LOSC_data/GW170729/V-V1_GWOSC_4KHZ_R1-1185389792-32.txt";
	//detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1185389792-32.txt";
	//detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1185389792-32.txt";
	//detector_files[2] =  "testing/data/V-V1_GWOSC_4KHZ_R1-1185389792-32.txt";
 	//double trigger_time= 1135136350.6;
 	double trigger_time = gps_time;
	double **psd = allocate_2D_array(num_detectors,psd_length);
	double **freqs = allocate_2D_array(num_detectors,psd_length);
	std::complex<double> **data = (std::complex<double> **)malloc(sizeof(std::complex<double> *)*num_detectors);
	for(int i =0; i<num_detectors; i++)
		data[i] = (std::complex<double>*)malloc(sizeof(std::complex<double>)*psd_length);

	allocate_LOSC_data(detector_files, psd_file, num_detectors, psd_length, datalength, trigger_time, data, psd, freqs);


	int *data_length= (int*)malloc(sizeof(int)*num_detectors);
	data_length[0] =psd_length;
	data_length[1] =psd_length;
	data_length[2] =psd_length;

	//#########################################################
	//mcmc options
	int dimension = 14;
	//double initial_pos[dimension]={.9, 2, 1.,std::log(500),std::log(9), .22,.4,.4,.1,.1,.1,.1,1,1};
	//double initial_pos[dimension]={.9, 5.2,-1.,std::log(2000),std::log(55), .22,.4,.4,.1,.1,.1,.1,1,1};
	//double initial_pos[dimension]={.9, 5.2,-1.,std::log(500),std::log(9), .22,.4,.4,.1,.1,.1,.1,1,1};
	double initial_pos[dimension]={.9, 5.2,-1.,std::log(300),std::log(30), .22,.4,.4,.1,.1,.1,.1,1,1};
	double *seeding_var = NULL;
	int n_steps = 30000;
	int chain_N=24 ;
	int max_thermo=12 ;
	int t0 = 10000;
	int nu = 100;
	std::string chain_alloc = "half_ensemble";
	//std::string chain_alloc = "cold";
	double ***output;
	output = allocate_3D_array( chain_N, n_steps, dimension );
	int swp_freq = 5;
	double chain_temps[chain_N];
	double c = 1.3;
	
	int Nmod = 0;
	int *bppe = NULL;
	int numThreads = 10;
	bool pool = true;
	//#########################################################
	//gw options
	//std::string generation_method = "dCS_IMRPhenomD_log";
	//std::string generation_method = "EdGB_IMRPhenomD_log";
	//std::string generation_method = "dCS_IMRPhenomD_root_alpha";
	std::string generation_method = "IMRPhenomPv2";
	//std::string generation_method = "EdGB_IMRPhenomD_root_alpha";
	
	
	std::string chainfile = "testing/data/mcmc_output_Pv2_dpt_alloc.csv";
	std::string statfilename = "testing/data/mcmc_statistics_Pv2_dpt_alloc.txt";
	std::string checkfile = "testing/data/mcmc_checkpoint_Pv2_dpt_alloc.csv";

	PTMCMC_MH_dynamic_PT_alloc_GW(output, dimension, n_steps, chain_N, max_thermo, initial_pos,seeding_var,chain_temps, 
			swp_freq, t0, nu, chain_alloc, test_lp_GW_Pv2,numThreads, pool,show_progress,
			num_detectors, 
			data, psd,freqs, data_length,gps_time, detectors,Nmod, bppe,
			generation_method,statfilename,"", "",checkfile);	


	//write_file(chainfile, output[0], n_steps, dimension);

	deallocate_3D_array(output, chain_N, n_steps, dimension);
	delete [] detectors;
	free(data_length);
	//free_LOSC_data(data, psd,freqs, num_detectors, length);
	deallocate_2D_array(psd,num_detectors, psd_length);
	deallocate_2D_array(freqs,num_detectors, psd_length);
	for(int i =0; i<num_detectors; i++)
		free(data[i]);
	free(data);
	delete [] detector_files;

}
void test33()
{
	std::string start_check = "testing/data/neil_mcmc_checkpoint1_dynamicPT.csv";
	std::string chain = "testing/data/neil_mcmc_output1_dynamicPT.csv";
	std::string check = "testing/data/neil_mcmc_checkpoint1_dynamicPT-real.csv";
	std::string stats = "testing/data/neil_mcmc_statistics1_dynamicPT-real.txt";
	int N_steps = 1000;
	int dimension = 2;
	int chain_N= 20;
	int swp_freq = 3;
	int threads = 10;
	bool show_prog = true;
	bool pool = true;
	double *temps = (double *)malloc(sizeof(double)*chain_N);
	load_temps_checkpoint_file(start_check,temps, chain_N);


	double ***output;
	output = allocate_3D_array( chain_N, N_steps, dimension );

	continue_PTMCMC_MH(start_check, output, N_steps, swp_freq, test_lp, log_neil_proj3,NULL, (void **)NULL,threads, pool, show_prog, stats, "", "", "",check,true);	
	int filecount = 1;
	std::string chainfilebase = "testing/data/neil_mcmc_output";
	for(int i =0; i<chain_N; i++){
		std::cout<<temps[i]<<std::endl;
		if(temps[i] == 1){
			write_file(chainfilebase+to_string(filecount)+"_dynamicPT.csv", output[i],N_steps, dimension);
			filecount+=1;
		}
	}
	deallocate_3D_array(output, chain_N, N_steps, dimension);
	free(temps);
}
void test32()
{
	int dimension = 2;
	double initial_pos[2]={1,-2.};
	double *seeding_var = NULL;

	
	int N_steps = 100000;
	int chain_N= 20;
	int max_chain_N_thermo= 4;
	int t0= 10000;
	int nu= 5;
	//std::string chain_dist_method = "half_ensemble";
	std::string chain_dist_method = "double";
	//std::string chain_dist_method = "cold";
	double ***output;
	output = allocate_3D_array( chain_N, N_steps, dimension );
	//double *initial_pos_ptr = initial_pos;
	int swp_freq = 3;
	//double chain_temps[chain_N] ={1,2,3,10,12};
	double chain_temps[chain_N];
	//double chain_temps[chain_N] ={1};
	//std::string autocorrfile = "";
	std::string chainfile = "testing/data/neil_mcmc_output1_dynamicPT.csv";
	std::string statfilename = "testing/data/neil_mcmc_statistics1_dynamicPT.txt";
	std::string checkpointfile = "testing/data/neil_mcmc_checkpoint1_dynamicPT.csv";
	
	int numThreads = 10;
	bool pool = false;
	
	PTMCMC_MH_dynamic_PT_alloc(output, dimension, N_steps, chain_N,max_chain_N_thermo, initial_pos,seeding_var,chain_temps, swp_freq, t0, nu,chain_dist_method,test_lp, log_neil_proj3,NULL,(void **)NULL,numThreads, pool,show_progress, statfilename,"", "",checkpointfile );	
	std::cout<<"ENDED"<<std::endl;
	std::cout<<"Chain temps: "<<std::endl;
	for(int i =0; i<chain_N; i++){
		std::cout<<chain_temps[i]<<std::endl;
	}

	//int filecount = 1;
	//std::string chainfilebase = "testing/data/neil_mcmc_output";
	//for(int i =0; i<chain_N; i++){
	//	if(chain_temps[i] == 1){
	//		write_file(chainfilebase+to_string(filecount)+"_dynamicPT.csv", output[i],N_steps, dimension);
	//		filecount+=1;
	//	}
	//}
	
	//autocorrfile = "testing/data/neil_auto_corr_mcmc2.csv";
	//chainfile = "testing/data/neil_mcmc_output2.csv";
	//statfilename = "testing/data/neil_mcmc_statistics2.txt";
	//MCMC_MH(output, dimension, N_steps, chain_N, initial_pos,chain_temps, swp_freq, test_lp, log_neil_proj32,NULL,statfilename,chainfile,autocorrfile );	
	//std::cout<<"ENDED"<<std::endl;
	


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
void test31()
{
	//std::string psd_file = "testing/data/GWTC1_GW150914_PSDs.dat.txt";
	//std::string psd_file = "testing/data/GWTC1_GW151226_PSDs.dat.txt";
	std::string psd_file = "/home/sperkins/Downloads/LOSC_data/GW151226/GWTC1_GW151226_PSDs.dat.txt";
	//std::string data_file = "testing/data/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//std::string psd_file = "testing/data/GWTC1_GW170729_PSDs.dat.txt";
	//int rows = 8032;
	//int cols = 3;
	int datalength = 131072;
	//double **psd = allocate_2D_array(rows, cols);
	//read_LOSC_PSD_file(psd_file, psd, rows, cols);
	//double data_start_time, duration, fs;
	int num_detectors = 2, psd_length = 8032, length;
	//int num_detectors = 3, psd_length = 4016, length;
	double gps_time = 1135136350.6;//TESTING -- gw151226
	//double gps_time = 1126259462;//TESTING -- gw150914
	//double gps_time = 1185389807.3;//TESTING -- gw170729
	std::string *detectors = new std::string[num_detectors];//(std::string*)malloc(sizeof(std::string)*50*num_detectors);
	detectors[0] = "Hanford";
	detectors[1] = "Livingston";
	//detectors[2] = "Virgo";
	std::string *detector_files = new std::string[num_detectors];
	//detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1135136335-32.txt";
	detector_files[0] =  "/home/sperkins/Downloads/LOSC_data/GW151226/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	detector_files[1] =  "/home/sperkins/Downloads/LOSC_data/GW151226/L-L1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//detector_files[2] =  "testing/data/V-V1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1185389792-32.txt";
	//detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1185389792-32.txt";
	//detector_files[2] =  "testing/data/V-V1_GWOSC_4KHZ_R1-1185389792-32.txt";
 	//double trigger_time= 1135136350.6;
 	double trigger_time = gps_time;
	double **psd = allocate_2D_array(num_detectors,psd_length);
	double **freqs = allocate_2D_array(num_detectors,psd_length);
	std::complex<double> **data = (std::complex<double> **)malloc(sizeof(std::complex<double> *)*num_detectors);
	for(int i =0; i<num_detectors; i++)
		data[i] = (std::complex<double>*)malloc(sizeof(std::complex<double>)*psd_length);

	allocate_LOSC_data(detector_files, psd_file, num_detectors, psd_length, datalength, trigger_time, data, psd, freqs);


	int *data_length= (int*)malloc(sizeof(int)*num_detectors);
	data_length[0] =psd_length;
	data_length[1] =psd_length;
	//data_length[2] =psd_length;

	//#########################################################
	//mcmc options
	int dimension = 15;
	//int dimension = 14;
	double initial_pos[dimension]={.9, 2, 1.,std::log(500),std::log(9), .22, .4,.4,.01,.01,.01,.01,1,1,2};
	//double initial_pos[dimension]={.9, 2, 1.,std::log(500),std::log(9), .22, .4,.4,.01,.01,.01,.01,1,1};
	//double initial_pos[dimension]={.9, 5, -1.,std::log(2000),std::log(50), .22, .4,.4,.01,.01,.01,.01,1,1,20};
	double *seeding_var = NULL;
	int n_steps = 100000;
	int n_steps_pt_alloc = 50000;
	int chain_N=20 ;
	int max_thermo_ensemble=10 ;
	double ***output;
	output = allocate_3D_array( chain_N, n_steps, dimension );
	double ***outputptalloc = allocate_3D_array( chain_N, n_steps_pt_alloc, dimension );
	int swp_freq = 5;
	double chain_temps[chain_N];
	double c = 1.3;
	chain_temps[0]=1.;
	for(int i =1; i < chain_N;  i ++)
		chain_temps[i] = c*chain_temps[i-1];
	
	int Nmod = 1;
	int *bppe = new int[Nmod];
	bppe[0] = -1;
	int numThreads = 6;
	bool pool = true;
	//#########################################################
	//gw options
	//std::string generation_method = "dCS_IMRPhenomD_log";
	//std::string generation_method = "EdGB_IMRPhenomD_log";
	//std::string generation_method = "dCS_IMRPhenomD_root_alpha";
	//std::string generation_method = "IMRPhenomPv2";
	std::string generation_method = "dCS_IMRPhenomPv2_root_alpha";
	//std::string generation_method = "EdGB_IMRPhenomD_root_alpha";
	
	
	std::string autocorrfile = "";
	//std::string chainfile_base = "testing/data/mcmc_output_dCS_Pv2_";
	//std::string statfilename = "testing/data/mcmc_statistics_dCS_Pv2.txt";
	//std::string checkfile = "testing/data/mcmc_checkpoint_dCS_Pv2.csv";
	//std::string checkfile_pt = "testing/data/mcmc_checkpoint_dCS_Pv2_dynamicPTalloc.csv";
	std::string chainfile_base = "testing/data/mcmc_output_dCS_Pv2_run2_";
	std::string statfilename = "testing/data/mcmc_statistics_dCS_Pv2_run2.txt";
	std::string checkfile = "testing/data/mcmc_checkpoint_dCS_Pv2_run2.csv";
	std::string checkfile_pt = "testing/data/mcmc_checkpoint_dCS_Pv2.csv";
	int t0 = 10000;
	int nu = 100;
	std::string chain_alloc = "half_ensemble";
	//PTMCMC_MH_dynamic_PT_alloc_GW(outputptalloc, dimension, n_steps_pt_alloc, chain_N, max_thermo_ensemble, initial_pos,seeding_var,chain_temps, 
	//		swp_freq, t0, nu, chain_alloc, test_lp_GW_Pv2_dCS_root_alpha,numThreads, pool,show_progress,
	//		num_detectors, 
	//		data, psd,freqs, data_length,gps_time, detectors,Nmod, bppe,
	//		generation_method,statfilename,"", "",checkfile_pt);	
	load_temps_checkpoint_file(checkfile_pt, chain_temps, chain_N);

	continue_PTMCMC_MH_GW(checkfile_pt,output, dimension, n_steps,  
			swp_freq, test_lp_GW_Pv2_dCS_root_alpha,numThreads, pool,show_progress,
			num_detectors, 
			data, psd,freqs, data_length,gps_time, detectors,Nmod, bppe,
			generation_method,statfilename,"",autocorrfile, "",checkfile,true);	

	std::string ind_chain_file;
	int criteria = 0 ;
	int file_ct = 0 ;
	while(criteria ==0){
		ind_chain_file = chainfile_base+std::to_string(file_ct)+".csv";
		criteria = std::remove(ind_chain_file.c_str());
		file_ct++;
	}
	file_ct = 0;
	for(int i =0; i < chain_N; i++){
		std::cout<<chain_temps[i]<<std::endl;
		if(chain_temps[i]==1){
			ind_chain_file = chainfile_base + to_string(file_ct) + ".csv";
			write_file(ind_chain_file , output[i], n_steps, dimension);
			file_ct++;
		}
	}

	deallocate_3D_array(output, chain_N, n_steps, dimension);
	deallocate_3D_array(outputptalloc, chain_N, n_steps_pt_alloc, dimension);
	delete [] detectors;
	free(data_length);
	deallocate_2D_array(psd,num_detectors, psd_length);
	deallocate_2D_array(freqs,num_detectors, psd_length);
	for(int i =0; i<num_detectors; i++)
		free(data[i]);
	free(data);
	delete [] detector_files;
	delete [] bppe;

}
void test30()
{
	//std::string psd_file = "testing/data/GWTC1_GW150914_PSDs.dat.txt";
	//std::string psd_file = "testing/data/GWTC1_GW151226_PSDs.dat.txt";
	//std::string data_file = "testing/data/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	std::string psd_file = "testing/data/GWTC1_GW170729_PSDs.dat.txt";
	//int rows = 8032;
	//int cols = 3;
	int datalength = 131075;
	//double **psd = allocate_2D_array(rows, cols);
	//read_LOSC_PSD_file(psd_file, psd, rows, cols);
	//double data_start_time, duration, fs;
	//int num_detectors = 2, psd_length = 8032, length;
	int num_detectors = 3, psd_length = 4016, length;
	//double gps_time = 1135136350.6;//TESTING -- gw151226
	//double gps_time = 1126259462;//TESTING -- gw150914
	double gps_time = 1185389807.3;//TESTING -- gw170729
	std::string *detectors = new std::string[num_detectors];//(std::string*)malloc(sizeof(std::string)*50*num_detectors);
	detectors[0] = "Hanford";
	detectors[1] = "Livingston";
	detectors[2] = "Virgo";
	std::string *detector_files = new std::string[num_detectors];
	//detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//detector_files[2] =  "testing/data/V-V1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1126259447-32.txt";
	detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1185389792-32.txt";
	detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1185389792-32.txt";
	detector_files[2] =  "testing/data/V-V1_GWOSC_4KHZ_R1-1185389792-32.txt";
 	//double trigger_time= 1135136350.6;
 	double trigger_time = gps_time;
	double **psd = allocate_2D_array(num_detectors,psd_length);
	double **freqs = allocate_2D_array(num_detectors,psd_length);
	std::complex<double> **data = (std::complex<double> **)malloc(sizeof(std::complex<double> *)*num_detectors);
	for(int i =0; i<num_detectors; i++)
		data[i] = (std::complex<double>*)malloc(sizeof(std::complex<double>)*psd_length);

	allocate_LOSC_data(detector_files, psd_file, num_detectors, psd_length, datalength, trigger_time, data, psd, freqs);


	int *data_length= (int*)malloc(sizeof(int)*num_detectors);
	data_length[0] =psd_length;
	data_length[1] =psd_length;
	data_length[2] =psd_length;

	//#########################################################
	//mcmc options
	int dimension = 15;
	double initial_pos[dimension]={.9, 2, 1.,std::log(500),std::log(9), .22, .4,.4,.01,.01,.01,.01,1,1,0};
	double *seeding_var = NULL;
	int n_steps = 50000;
	int chain_N=12 ;
	double ***output;
	output = allocate_3D_array( chain_N, n_steps, dimension );
	int swp_freq = 5;
	double chain_temps[chain_N];
	double c = 1.1;
	chain_temps[0]=1.;
	//chain_temps[1] = c;
	//chain_temps[2] = c*c;
	//chain_temps[3] = c*c*c;
	//chain_temps[4] = 1;
	//chain_temps[5] = c;
	//chain_temps[6] = c*c;
	//chain_temps[7] = c*c*c;
	//chain_temps[1] = 1.1;
	//chain_temps[2] = 1.2;
	//chain_temps[3] = 1.3;
	//chain_temps[4] = 1.4;
	//chain_temps[5] = 1.5;
	//chain_temps[6] = 1.6;
	//chain_temps[7] = 1.7;
	for(int i =1; i < chain_N;  i ++)
		chain_temps[i] = c*chain_temps[i-1];
	//double c = 1./std::sqrt(dimension);//optimum spacing of inverse temperature
	//for(int i =1; i < chain_N;  i ++){
	//	chain_temps[i] = chain_temps[i-1]/(1-chain_temps[i-1]*c);
	//	std::cout<<chain_temps[i]<<std::endl;
	//}
	
	int Nmod = 1;
	int *bppe = new int[Nmod];
	bppe[0] = -1;
	int numThreads = 10;
	bool pool = true;
	//#########################################################
	//gw options
	//std::string generation_method = "dCS_IMRPhenomD_log";
	//std::string generation_method = "EdGB_IMRPhenomD_log";
	//std::string generation_method = "dCS_IMRPhenomD_root_alpha";
	//std::string generation_method = "IMRPhenomPv2";
	std::string generation_method = "ppE_IMRPhenomPv2_Inspiral";
	//std::string generation_method = "EdGB_IMRPhenomD_root_alpha";
	
	
	//std::string autocorrfile = "";
	//std::string autocorrfile = "testing/data/auto_corr_mcmc_dCS.csv";
	//std::string chainfile = "testing/data/mcmc_output_dCS.csv";
	//std::string statfilename = "testing/data/mcmc_statistics_dCS.txt";
	//std::string checkfile = "testing/data/mcmc_checkpoint_dCS.csv";
	//std::string autocorrfile = "testing/data/auto_corr_mcmc_Pv2.csv";
	//std::string chainfile = "testing/data/mcmc_output_Pv2.csv";
	//std::string statfilename = "testing/data/mcmc_statistics_Pv2.txt";
	//std::string checkfile = "testing/data/mcmc_checkpoint_Pv2.csv";
	std::string autocorrfile = "testing/data/auto_corr_mcmc_ppE_Pv2.csv";
	std::string chainfile = "testing/data/mcmc_output_ppE_Pv2.csv";
	std::string statfilename = "testing/data/mcmc_statistics_ppE_Pv2.txt";
	std::string checkfile = "testing/data/mcmc_checkpoint_ppE_Pv2.csv";
	//std::string autocorrfile = "testing/data/auto_corr_mcmc_EdGB.csv";
	//std::string chainfile = "testing/data/mcmc_output_EdGB.csv";
	//std::string statfilename = "testing/data/mcmc_statistics_EdGB.txt";
	//std::string checkfile = "testing/data/mcmc_checkpoint_EdGB.csv";

	PTMCMC_MH_GW(output, dimension, n_steps, chain_N, initial_pos,seeding_var,chain_temps, 
			swp_freq, test_lp_GW_Pv2_ppE,numThreads, pool,show_progress,
			num_detectors, 
			data, psd,freqs, data_length,gps_time, detectors,Nmod, bppe,
			generation_method,statfilename,"",autocorrfile, "",checkfile);	

	//double **output_transform=(double **)malloc(sizeof(double*)*n_steps);
	//for (int j =0; j<n_steps; j++)
	//	output_transform[j] = (double *)malloc(sizeof(double)*dimension);

	//for(int j = 0; j<n_steps;j++){
	//	output_transform[j][0]=output[0][j][0];
	//	output_transform[j][1]=output[0][j][1];
	//	output_transform[j][2]=output[0][j][2];
	//	output_transform[j][3]=std::exp(output[0][j][3]);
	//	output_transform[j][4]=std::exp(output[0][j][4]);
	//	output_transform[j][5]=output[0][j][5];
	//	output_transform[j][6]=output[0][j][6];
	//	output_transform[j][7]=output[0][j][7];
	//	//output_transform[j][8]=std::exp(output[0][j][8]);
	//	output_transform[j][8]=output[0][j][8];
	//}
	write_file(chainfile, output[0], n_steps, dimension);
	//std::string chain2 = "testing/data/mcmc_output_dCS_2.csv";
	//write_file(chain2, output[4], n_steps, dimension);
	//output hottest chain too
	//chainfile = "testing/data/mcmc_output_EdGB_hot.csv";
	//chainfile = "testing/data/mcmc_output_dCS_hot.csv";
	//for(int j = 0; j<n_steps;j++){
	//	output_transform[j][0]=output[chain_N-1][j][0];
	//	output_transform[j][1]=output[chain_N-1][j][1];
	//	output_transform[j][2]=output[chain_N-1][j][2];
	//	output_transform[j][3]=std::exp(output[chain_N-1][j][3]);
	//	output_transform[j][4]=std::exp(output[chain_N-1][j][4]);
	//	output_transform[j][5]=output[chain_N-1][j][5];
	//	output_transform[j][6]=output[chain_N-1][j][6];
	//	output_transform[j][7]=output[chain_N-1][j][7];
	//	//output_transform[j][8]=std::exp(output[chain_N-1][j][8]);
	//	output_transform[j][8]=output[chain_N-1][j][8];
	//}
	//write_file(chainfile, output_transform, n_steps, dimension);

	deallocate_3D_array(output, chain_N, n_steps, dimension);
	//for(int i =0; i< n_steps; i++){
	//	free(output_transform[i]);
	//}
	//free(output_transform);
	delete [] detectors;
	free(data_length);
	//free_LOSC_data(data, psd,freqs, num_detectors, length);
	deallocate_2D_array(psd,num_detectors, psd_length);
	deallocate_2D_array(freqs,num_detectors, psd_length);
	for(int i =0; i<num_detectors; i++)
		free(data[i]);
	free(data);
	delete [] detector_files;
	delete [] bppe;
	//deallocate_2D_array(psd, rows, cols);
	//free(data);
	
	//int length = 2000;
	//double x[length];
	//double y[length];
	//double xlim = 10.;
	//double xstart=0;
	//double xstep = (xlim-xstart)/length;
	//for (int i =0; i<length; i++){
	//	x[i]= i*xstep;
	//	y[i]= x[i]*x[i];
	//}
	//double sum = simpsons_sum(xstep, length, y);
	//std::cout<<sum<<std::endl;

}

void test29()
{
	//std::string psd_file = "testing/data/GWTC1_GW150914_PSDs.dat.txt";
	//std::string psd_file = "testing/data/GWTC1_GW151226_PSDs.dat.txt";
	//std::string data_file = "testing/data/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//std::string psd_file = "testing/data/GWTC1_GW170729_PSDs.dat.txt";
	std::string psd_file = "/Users/sperkins/Downloads/LOSC_data/GW170729/GWTC1_GW170729_PSDs.dat.txt";
	//int rows = 8032;
	//int cols = 3;
	int datalength = 131075;
	//double **psd = allocate_2D_array(rows, cols);
	//read_LOSC_PSD_file(psd_file, psd, rows, cols);
	//double data_start_time, duration, fs;
	//int num_detectors = 2, psd_length = 8032, length;
	int num_detectors = 3, psd_length = 4016, length;
	//int num_detectors = 2, psd_length = 4016, length;
	//double gps_time = 1135136350.6;//TESTING -- gw151226
	//double gps_time = 1126259462;//TESTING -- gw150914
	double gps_time = 1185389807.3;//TESTING -- gw170729
	std::string *detectors = new std::string[num_detectors];//(std::string*)malloc(sizeof(std::string)*50*num_detectors);
	detectors[0] = "Hanford";
	detectors[1] = "Livingston";
	detectors[2] = "Virgo";
	std::string *detector_files = new std::string[num_detectors];
	//detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1185389792-32.txt";
	//detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1185389792-32.txt";
	//detector_files[2] =  "testing/data/V-V1_GWOSC_4KHZ_R1-1185389792-32.txt";

	detector_files[0] =  "/Users/sperkins/Downloads/LOSC_data/GW170729/H-H1_GWOSC_4KHZ_R1-1185389792-32.txt";
	detector_files[1] =  "/Users/sperkins/Downloads/LOSC_data/GW170729/L-L1_GWOSC_4KHZ_R1-1185389792-32.txt";
	detector_files[2] =  "/Users/sperkins/Downloads/LOSC_data/GW170729/V-V1_GWOSC_4KHZ_R1-1185389792-32.txt";
 	//double trigger_time= 1135136350.6;
 	double trigger_time = gps_time;
	double **psd = allocate_2D_array(num_detectors,psd_length);
	double **freqs = allocate_2D_array(num_detectors,psd_length);
	std::complex<double> **data = (std::complex<double> **)malloc(sizeof(std::complex<double> *)*num_detectors);
	for(int i =0; i<num_detectors; i++)
		data[i] = (std::complex<double>*)malloc(sizeof(std::complex<double>)*psd_length);

	allocate_LOSC_data(detector_files, psd_file, num_detectors, psd_length, datalength, trigger_time, data, psd, freqs);


	int *data_length= (int*)malloc(sizeof(int)*num_detectors);
	data_length[0] =psd_length;
	data_length[1] =psd_length;
	data_length[2] =psd_length;

	//#########################################################
	//mcmc options
	int dimension = 14;
	//double initial_pos[dimension]={.9, 2, 1.,std::log(500),std::log(9), .22,.4,.4,.1,.1,.1,.1,1,1};
	double initial_pos[dimension]={.9, 5.2,-1.,std::log(2000),std::log(55), .22,.4,.4,.1,.1,.1,.1,1,1};
	double *seeding_var = NULL;
	int n_steps = 10000;
	int chain_N=12 ;
	double ***output;
	output = allocate_3D_array( chain_N, n_steps, dimension );
	int swp_freq = 5;
	double chain_temps[chain_N];
	double c = 1.3;
	chain_temps[0]=1.;
	//chain_temps[1] = c;
	//chain_temps[2] = c*c;
	//chain_temps[3] = c*c*c;
	//chain_temps[4] = 1;
	//chain_temps[5] = c;
	//chain_temps[6] = c*c;
	//chain_temps[7] = c*c*c;
	//chain_temps[1] = 1.1;
	//chain_temps[2] = 1.2;
	//chain_temps[3] = 1.3;
	//chain_temps[4] = 1.4;
	//chain_temps[5] = 1.5;
	//chain_temps[6] = 1.6;
	//chain_temps[7] = 1.7;
	for(int i =1; i < chain_N;  i ++)
		chain_temps[i] = c*chain_temps[i-1];
	//double c = 1./std::sqrt(dimension);//optimum spacing of inverse temperature
	//for(int i =1; i < chain_N;  i ++){
	//	chain_temps[i] = chain_temps[i-1]/(1-chain_temps[i-1]*c);
	//	std::cout<<chain_temps[i]<<std::endl;
	//}
	
	int Nmod = 0;
	int *bppe = NULL;
	int numThreads = 10;
	bool pool = true;
	//#########################################################
	//gw options
	//std::string generation_method = "dCS_IMRPhenomD_log";
	//std::string generation_method = "EdGB_IMRPhenomD_log";
	//std::string generation_method = "dCS_IMRPhenomD_root_alpha";
	std::string generation_method = "IMRPhenomPv2";
	//std::string generation_method = "EdGB_IMRPhenomD_root_alpha";
	
	
	//std::string autocorrfile = "";
	//std::string autocorrfile = "testing/data/auto_corr_mcmc_dCS.csv";
	//std::string chainfile = "testing/data/mcmc_output_dCS.csv";
	//std::string statfilename = "testing/data/mcmc_statistics_dCS.txt";
	//std::string checkfile = "testing/data/mcmc_checkpoint_dCS.csv";
	std::string autocorrfile = "testing/data/auto_corr_mcmc_Pv2.csv";
	std::string chainfile = "testing/data/mcmc_output_Pv2.csv";
	std::string statfilename = "testing/data/mcmc_statistics_Pv2.txt";
	std::string checkfile = "testing/data/mcmc_checkpoint_Pv2.csv";
	//std::string autocorrfile = "testing/data/auto_corr_mcmc_EdGB.csv";
	//std::string chainfile = "testing/data/mcmc_output_EdGB.csv";
	//std::string statfilename = "testing/data/mcmc_statistics_EdGB.txt";
	//std::string checkfile = "testing/data/mcmc_checkpoint_EdGB.csv";

	PTMCMC_MH_GW(output, dimension, n_steps, chain_N, initial_pos,seeding_var,chain_temps, 
			swp_freq, test_lp_GW_Pv2,numThreads, pool,show_progress,
			num_detectors, 
			data, psd,freqs, data_length,gps_time, detectors,Nmod, bppe,
			generation_method,statfilename,"",autocorrfile, "",checkfile);	


	write_file(chainfile, output[0], n_steps, dimension);
	write_file("testing/data/mcmc_output_Pv2_hot.csv", output[chain_N-1], n_steps, dimension);

	deallocate_3D_array(output, chain_N, n_steps, dimension);
	delete [] detectors;
	free(data_length);
	//free_LOSC_data(data, psd,freqs, num_detectors, length);
	deallocate_2D_array(psd,num_detectors, psd_length);
	deallocate_2D_array(freqs,num_detectors, psd_length);
	for(int i =0; i<num_detectors; i++)
		free(data[i]);
	free(data);
	delete [] detector_files;
	//deallocate_2D_array(psd, rows, cols);
	//free(data);
	
	//int length = 2000;
	//double x[length];
	//double y[length];
	//double xlim = 10.;
	//double xstart=0;
	//double xstep = (xlim-xstart)/length;
	//for (int i =0; i<length; i++){
	//	x[i]= i*xstep;
	//	y[i]= x[i]*x[i];
	//}
	//double sum = simpsons_sum(xstep, length, y);
	//std::cout<<sum<<std::endl;

}
void test28()
{
	int length = 4000;
	//int num_detectors =3;
	int num_detectors =2;
	std::string *detectors = new std::string[num_detectors];//(std::string*)malloc(sizeof(std::string)*50*num_detectors);
	detectors[0] = "Hanford";
	detectors[1] = "Livingston";
	//detectors[2] = "Virgo";

	double gps_time = 1135136350.6;//TESTING -- gw151226

	std::complex<double> **data= (std::complex<double>**)malloc(
			sizeof(std::complex<double>*)*num_detectors);
	double **psd = (double **)malloc(sizeof(double *)*num_detectors);
	double **frequencies = (double **)malloc(sizeof(double *)*num_detectors);
	int *data_length= (int*)malloc(sizeof(int)*num_detectors);
	data_length[0] =length;
	data_length[1] =length;
	//data_length[2] =length;

	for (int i =0; i<num_detectors; i++){
		data[i] = (std::complex<double> *)malloc(
			sizeof(std::complex<double>)*data_length[i]);
		
		psd[i] = (double *)malloc(sizeof(double)*data_length[i]);
		frequencies[i] = (double *)malloc(sizeof(double)*data_length[i]);
	}
	//#########################################################
	//Make trial data
	gen_params params;
	//double RA = 5.;
	//double DEC = 1.;
	double RA = 2;
	double DEC = .4;
	double chirpm = 9.71;
	double eta =.2;
	params.mass1 = calculate_mass1(chirpm,eta);
	params.mass2 = calculate_mass2(chirpm,eta);
	complex<double> waveformout[length];
	params.spin1[0] = 0;
	params.spin1[1] = 0;
	params.spin1[2] = .43;
	params.spin2[0] = 0;
	params.spin2[1] = 0;
	params.spin2[2] = .0;
	params.phic = .0;
	double tc = 2;
	params.Luminosity_Distance = 125.;
	params.NSflag1 = false;
	params.NSflag2 = false;
	params.incl_angle = 0.;//M_PI/3.;
	params.sky_average=false;
	params.bppe = new int[1];
	params.betappe = new double[1];
	params.bppe[0] = -1;
	params.betappe[0] = (pow(0./(3e5),4));
	//params.betappe[0] = log(pow(100./(3e5),4));
	params.Nmod = 1;
	//std::string injection_method = "dCS_IMRPhenomD";
	std::string injection_method = "IMRPhenomD";
	//params.f_ref = 30.5011;
	//params.phiRef =58.944425/2.;
	
	//############################################################
	double fhigh =500;
	double flow =10;
	double df = (fhigh-flow)/(length-1);
	//double *freq = (double *)malloc(sizeof(double) * length);
	double freq[length];

	cout<<"Freq spacing "<<df<<endl;

	for(int i=0;i<length;i++)
		freq[i]=flow+i*df;
	//############################################################
	double noise[length];
	populate_noise(freq,"Hanford_O1_fitted", noise,length);
	for (int i =0; i<length;i++){
		noise[i] = noise[i]*noise[i];
	}
	//############################################################

	params.phi = 0;
	params.theta = 0;
	params.tc = tc; 
	celestial_horizon_transform(RA, DEC, gps_time, "Hanford", &params.phi, &params.theta);
	
	fourier_detector_response(freq,length, waveformout,"Hanford",injection_method, &params);
	double temptheta = params.theta;
	for(int j = 0; j<data_length[0]; j++){
		frequencies[0][j] = freq[j];	
		psd[0][j] = noise[j];	
		data[0][j] = waveformout[j];	
	}
	//double snr2 = pow(calculate_snr("Hanford_O1_fitted",waveformout, freq, length),2);
	double snr2 = pow(data_snr_maximized_extrinsic(freq,length,waveformout, noise,"Hanford", injection_method,&params),2);

	celestial_horizon_transform(RA, DEC, gps_time, "Livingston", &params.phi, &params.theta);
	params.tc = tc+ DTOA(temptheta, params.theta, "Hanford","Livingston"); 
	std::cout<<"Time at Livingston of injection: "<<params.tc<<std::endl;
	fourier_detector_response(freq,length, waveformout,"Livingston",injection_method, &params);

	for(int j = 0; j<data_length[0]; j++){
		frequencies[1][j] = freq[j];	
		psd[1][j] = (noise[j]);	
		data[1][j] = waveformout[j];	
	}
	//snr2+=pow( calculate_snr("Hanford_O1_fitted",waveformout, freq, length),2);
	snr2 += pow(data_snr_maximized_extrinsic(freq,length,waveformout, noise,"Livingston", injection_method,&params),2);
	//celestial_horizon_transform(RA, DEC, gps_time, "Virgo", &params.phi, &params.theta);
	//params.tc = tc+ DTOA(temptheta, params.theta, "Hanford","Virgo"); 
	//std::cout<<"Time at Virgo of injection: "<<params.tc<<std::endl;
	//fourier_detector_response(freq,length, waveformout,"Virgo",injection_method, &params);

	//for(int j = 0; j<data_length[0]; j++){
	//	frequencies[2][j] = freq[j];	
	//	psd[2][j] = (noise[j]);	
	//	data[2][j] = waveformout[j];	
	//}
	////#########################################################
	////snr2+=pow( calculate_snr("Hanford_O1_fitted",waveformout, freq, length),2);
	//snr2 += pow(data_snr_maximized_extrinsic(freq,length,waveformout, noise,"Virgo", injection_method,&params),2);
	std::cout<<"SNR of injection: "<<sqrt(snr2)<<std::endl;
	//snr = data_snr_maximized_extrinsic(freq,length, waveformout,"Hanford_O1_fitted",method,params );
	//std::cout<<"SNR of injection calculated as data: "<<snr<<std::endl;

	delete [] params.bppe;
	delete [] params.betappe;
	
	//#########################################################
	//mcmc options
	int dimension = 9;

	int n_steps = 150000;
	int chain_N= 10;
	double ***output;
	output = allocate_3D_array( chain_N, n_steps, dimension );
	int swp_freq = 3;
	
	int Nmod = 1;
	int *bppe = new int[1];
	bppe[0] = -1;
	int numThreads = 4;
	bool pool = true;
	//#########################################################
	//gw options
	//std::string generation_method = "dCS_IMRPhenomD_log";
	std::string generation_method = "dCS_IMRPhenomD_root_alpha";
	
	std::string iteration="3";
	std::string previteration="2";
	
	//std::string autocorrfile = "testing/data/auto_corr_mcmc_injection"+iteration+".csv";
	////std::string autocorrfile = "";
	//std::string chainfile = "testing/data/mcmc_output_injection"+iteration+".csv";
	//std::string statfilename = "testing/data/mcmc_statistics_injection"+iteration+".txt";
	//std::string checkfile = "testing/data/mcmc_checkpoint_injection"+iteration+".csv";
	//std::string startcheckfile = "testing/data/mcmc_checkpoint_injection"+previteration+".csv";
	std::string autocorrfile = "testing/data/auto_corr_mcmc_injection"+iteration+"_1226.csv";
	//std::string autocorrfile = "";
	std::string chainfile = "testing/data/mcmc_output_injection"+iteration+"_1226.csv";
	std::string statfilename = "testing/data/mcmc_statistics_injection"+iteration+"_1226.txt";
	std::string checkfile = "testing/data/mcmc_checkpoint_injection"+iteration+"_1226.csv";
	std::string startcheckfile = "testing/data/mcmc_checkpoint_injection"+previteration+"_1226.csv";

	continue_PTMCMC_MH_GW(startcheckfile,output, dimension, n_steps, 
			swp_freq, test_lp_GW_dCS_root_alpha,numThreads, pool,show_progress,
			num_detectors, 
			data, psd,frequencies, data_length,gps_time, detectors,Nmod, bppe,
			generation_method,statfilename,"","","",checkfile,true);	
	
	std::cout<<"ended"<<std::endl;
	double **output_transform=(double **)malloc(sizeof(double*)*n_steps);
	for (int j =0; j<n_steps; j++)
		output_transform[j] = (double *)malloc(sizeof(double)*dimension);

	for(int j = 0; j<n_steps;j++){
			output_transform[j][0]=output[0][j][0];
			output_transform[j][1]=output[0][j][1];
			output_transform[j][2]=output[0][j][2];
			output_transform[j][3]=std::exp(output[0][j][3]);
			output_transform[j][4]=std::exp(output[0][j][4]);
			output_transform[j][5]=output[0][j][5];
			output_transform[j][6]=output[0][j][6];
			output_transform[j][7]=output[0][j][7];
			//output_transform[j][8]=exp(output[0][j][8]);
			output_transform[j][8]=(output[0][j][8]);
	}
	write_file(chainfile, output_transform, n_steps, dimension);

	int segs = 10;
	double target_corr = .01;
	//write_file_auto_corr_from_data_file_accel(autocorrfile, chainfile,dimension,n_steps,segs,target_corr);
	//write_auto_corr_file_from_data_file(autocorrfile, chainfile,n_steps,dimension,segs,target_corr,numThreads);

	//output hottest chain too
	chainfile = "testing/data/mcmc_output_injection"+iteration+"_1226_hot.csv";
	for(int j = 0; j<n_steps;j++){
			output_transform[j][0]=output[chain_N-1][j][0];
			output_transform[j][1]=output[chain_N-1][j][1];
			output_transform[j][2]=output[chain_N-1][j][2];
			output_transform[j][3]=std::exp(output[chain_N-1][j][3]);
			output_transform[j][4]=std::exp(output[chain_N-1][j][4]);
			output_transform[j][5]=output[chain_N-1][j][5];
			output_transform[j][6]=output[chain_N-1][j][6];
			output_transform[j][7]=output[chain_N-1][j][7];
			output_transform[j][8]=(output[chain_N-1][j][8]);
	}
	write_file(chainfile, output_transform, n_steps, dimension);

	deallocate_3D_array(output, chain_N, n_steps, dimension);
	for(int i =0; i< num_detectors; i++){
		free(data[i]);
		free(psd[i]);
		free(frequencies[i]);
	}
	for(int i =0; i< n_steps; i++){
		free(output_transform[i]);
	}
	free(output_transform);
	free(data);
	delete [] bppe;
	free(psd);
	free(frequencies );
	delete [] detectors;
	free(data_length);
}
void test27()
{
	//std::string psd_file = "testing/data/GWTC1_GW150914_PSDs.dat.txt";
	//std::string data_file = "testing/data/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//std::string psd_file = "testing/data/GWTC1_GW170729_PSDs.dat.txt";
	std::string psd_file = "testing/data/GWTC1_GW151226_PSDs.dat.txt";
	//std::string psd_file = "/Users/sperkins/Downloads/LOSC_data/GW170608/GWTC1_GW170608_PSDs.dat.txt";
	//std::string psd_file = "/Users/sperkins/Downloads/LOSC_data/GW150914/GWTC1_GW150914_PSDs.dat.txt";
	//int rows = 8032;
	//int cols = 3;
	int datalength = 131075;
	//double **psd = allocate_2D_array(rows, cols);
	//read_LOSC_PSD_file(psd_file, psd, rows, cols);
	//double data_start_time, duration, fs;
	int num_detectors = 2, psd_length = 8032, length;
	//int num_detectors = 2, psd_length = 16064, length;
	//int num_detectors = 3, psd_length = 4016, length;
	//double gps_time = 1126259462;//TESTING -- gw150914
	double gps_time = 1135136350.6;//TESTING -- gw151226
	//double gps_time = 1185389807.3;//TESTING -- gw170729
	//double gps_time = 1180922494.5;//TESTING -- gw170608
	std::string *detectors = new std::string[num_detectors];//(std::string*)malloc(sizeof(std::string)*50*num_detectors);
	detectors[0] = "Hanford";
	detectors[1] = "Livingston";
	//detectors[2] = "Virgo";
	std::string *detector_files = new std::string[num_detectors];
	detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//detector_files[0] =  "/Users/sperkins/Downloads/LOSC_data/GW170608/H-H1_GWOSC_4KHZ_R1-1180922479-32.txt";
	//detector_files[1] =  "/Users/sperkins/Downloads/LOSC_data/GW170608/L-L1_GWOSC_4KHZ_R1-1180922479-32.txt";
	//detector_files[0] =  "/Users/sperkins/Downloads/LOSC_data/GW150914/H-H1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[1] =  "/Users/sperkins/Downloads/LOSC_data/GW150914/L-L1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[2] =  "testing/data/V-V1_GWOSC_4KHZ_R1-1185389792-32.txt";
 	//double trigger_time= 1135136350.6;
 	double trigger_time = gps_time;
	double **psd = allocate_2D_array(num_detectors,psd_length);
	double **freqs = allocate_2D_array(num_detectors,psd_length);
	std::complex<double> **data = (std::complex<double> **)malloc(sizeof(std::complex<double> *)*num_detectors);
	for(int i =0; i<num_detectors; i++)
		data[i] = (std::complex<double>*)malloc(sizeof(std::complex<double>)*psd_length);

	allocate_LOSC_data(detector_files, psd_file, num_detectors, psd_length, datalength, trigger_time, data, psd, freqs);


	int *data_length= (int*)malloc(sizeof(int)*num_detectors);
	data_length[0] =psd_length;
	data_length[1] =psd_length;
	//data_length[2] =psd_length;

	//#########################################################
	//mcmc options
	int dimension = 8;
	//double initial_pos[dimension]={.3, 2., -0.2,log(400),log(30), .24,- .0,-.0};
	//double initial_pos[dimension]={-.9, 2, -1.2,log(410),log(30), .24,-.4,.3};
	//double initial_pos[dimension]={-.0, 0, -0,log(500),log(50), .2,-.0,.0};
	//double initial_pos[dimension]={-.99, 2, -1.2,log(410),log(30.78), .24,-.4,.3};
	int n_steps = 100000;
	int chain_N= 10;
	double ***output;
	output = allocate_3D_array( chain_N, n_steps, dimension );
	int swp_freq = 5;
	
	int Nmod = 0;
	int *bppe = NULL;
	int numThreads = 4;
	bool pool = true;
	//bool pool = false;
	//#########################################################
	//gw options
	std::string generation_method = "IMRPhenomD";
	
	std::string iteration="3";	
	std::string previteration="2";	
	
	std::string autocorrfile = "testing/data/auto_corr_mcmc_DFull"+iteration+".csv";
	//std::string autocorrfile = "";
	std::string chainfile = "testing/data/mcmc_output_DFull"+iteration+".csv";
	std::string statfilename = "testing/data/mcmc_statistics_DFull"+iteration+".txt";
	std::string checkfile = "testing/data/mcmc_checkpoint_DFull"+iteration+".csv";
	std::string start_checkfile = "testing/data/mcmc_checkpoint_DFull"+previteration+".csv";
	//std::string autocorrfile = "testing/data/auto_corr_mcmc_DFull.csv";
	////std::string autocorrfile = "";
	//std::string chainfile = "testing/data/mcmc_output_DFull.csv";
	//std::string statfilename = "testing/data/mcmc_statistics_DFull.txt";
	//std::string checkfile = "testing/data/mcmc_checkpoint_DFull.csv";

	continue_PTMCMC_MH_GW(start_checkfile,output, dimension, n_steps, 
			swp_freq, test_lp_GW_DFull,numThreads, pool,show_progress,
			num_detectors, 
			data, psd,freqs, data_length,gps_time, detectors,Nmod, bppe,
			generation_method,statfilename,"","","",checkfile,true);	
	std::cout<<"ended"<<std::endl;

	double **output_transform=(double **)malloc(sizeof(double*)*n_steps);
	for (int j =0; j<n_steps; j++)
		output_transform[j] = (double *)malloc(sizeof(double)*dimension);

	for(int j = 0; j<n_steps;j++){
			output_transform[j][0]=output[0][j][0];
			output_transform[j][1]=output[0][j][1];
			output_transform[j][2]=output[0][j][2];
			output_transform[j][3]=std::exp(output[0][j][3]);
			output_transform[j][4]=std::exp(output[0][j][4]);
			output_transform[j][5]=output[0][j][5];
			output_transform[j][6]=output[0][j][6];
			output_transform[j][7]=output[0][j][7];
	}
	write_file(chainfile, output_transform, n_steps, dimension);
	//output hottest chain too
	chainfile = "testing/data/mcmc_output_DFull_hot"+iteration+".csv";
	//chainfile = "testing/data/mcmc_output_DFull_hot.csv";
	for(int j = 0; j<n_steps;j++){
			output_transform[j][0]=output[chain_N-1][j][0];
			output_transform[j][1]=output[chain_N-1][j][1];
			output_transform[j][2]=output[chain_N-1][j][2];
			output_transform[j][3]=std::exp(output[chain_N-1][j][3]);
			output_transform[j][4]=std::exp(output[chain_N-1][j][4]);
			output_transform[j][5]=output[chain_N-1][j][5];
			output_transform[j][6]=output[chain_N-1][j][6];
			output_transform[j][7]=output[chain_N-1][j][7];
	}
	write_file(chainfile, output_transform, n_steps, dimension);

	deallocate_3D_array(output, chain_N, n_steps, dimension);
	for(int i =0; i< n_steps; i++){
		free(output_transform[i]);
	}
	free(output_transform);
	delete [] detectors;
	free(data_length);
	//free_LOSC_data(data, psd,freqs, num_detectors, length);
	deallocate_2D_array(psd,num_detectors, psd_length);
	deallocate_2D_array(freqs,num_detectors, psd_length);
	for(int i =0; i<num_detectors; i++)
		free(data[i]);
	free(data);
	delete [] detector_files;
	//deallocate_2D_array(psd, rows, cols);
	//free(data);
	
	//int length = 2000;
	//double x[length];
	//double y[length];
	//double xlim = 10.;
	//double xstart=0;
	//double xstep = (xlim-xstart)/length;
	//for (int i =0; i<length; i++){
	//	x[i]= i*xstep;
	//	y[i]= x[i]*x[i];
	//}
	//double sum = simpsons_sum(xstep, length, y);
	//std::cout<<sum<<std::endl;

}

void test26()
{
	//std::string psd_file = "testing/data/GWTC1_GW150914_PSDs.dat.txt";
	std::string psd_file = "testing/data/GWTC1_GW151226_PSDs.dat.txt";
	//std::string data_file = "testing/data/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//std::string psd_file = "testing/data/GWTC1_GW170729_PSDs.dat.txt";
	int datalength = 131075;
	int num_detectors = 2, psd_length = 8032, length;
	//int num_detectors = 3, psd_length = 4016, length;
	double gps_time = 1135136350.6;//TESTING -- gw151226
	//double gps_time = 1126259462;//TESTING -- gw150914
	//double gps_time = 1185389807.3;//TESTING -- gw170729
	std::string *detectors = new std::string[num_detectors];//(std::string*)malloc(sizeof(std::string)*50*num_detectors);
	detectors[0] = "Hanford";
	detectors[1] = "Livingston";
	//detectors[2] = "Virgo";
	std::string *detector_files = new std::string[num_detectors];
	detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1185389792-32.txt";
	//detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1185389792-32.txt";
	//detector_files[2] =  "testing/data/V-V1_GWOSC_4KHZ_R1-1185389792-32.txt";
 	//double trigger_time= 1135136350.6;
 	double trigger_time = gps_time;
	double **psd = allocate_2D_array(num_detectors,psd_length);
	double **freqs = allocate_2D_array(num_detectors,psd_length);
	std::complex<double> **data = (std::complex<double> **)malloc(sizeof(std::complex<double> *)*num_detectors);
	for(int i =0; i<num_detectors; i++)
		data[i] = (std::complex<double>*)malloc(sizeof(std::complex<double>)*psd_length);

	allocate_LOSC_data(detector_files, psd_file, num_detectors, psd_length, datalength, trigger_time, data, psd, freqs);


	int *data_length= (int*)malloc(sizeof(int)*num_detectors);
	data_length[0] =psd_length;
	data_length[1] =psd_length;
	//data_length[2] =psd_length;

	//#########################################################
	//mcmc options
	int dimension = 9;
	//double initial_pos[dimension]={.0, 1, 0.,log(300),log(10), .2,- .0,-.0, log(pow(MPC_SEC,4)*pow(50000,4))};
	//double initial_pos[dimension]={.0, 1, 0.,log(300),log(10), .2,- .0,-.0, log(pow(MPC_SEC,4)*pow(5,4))};
	//double initial_pos[dimension]={.0, 1, 0.,log(300),log(10), .2,- .0,-.0, -5};
	//double initial_pos[dimension]={.0, 1, 0.,log(300),log(10), .2,- .0,-.0,4* log(50000/(3e8))};
	int n_steps = 200000;
	int chain_N = 10;
	double ***output;
	output = allocate_3D_array( chain_N, n_steps, dimension );
	int swp_freq = 5;
	int Nmod = 1;
	int *bppe = new int[Nmod];
	bppe[0] = -1;
	int numThreads = 12;
	bool pool = true;
	//#########################################################
	//gw options
	//std::string generation_method = "ppE_IMRPhenomD_Inspiral";
	std::string generation_method = "EdGB_IMRPhenomD_log";
	//std::string generation_method = "dCS_IMRPhenomD_log";
	
	
	std::string iteration = "4";
	std::string previteration = "3";
	//std::string autocorrfile = "";
	//std::string autocorrfile = "testing/data/auto_corr_mcmc_ppE.csv";
	//std::string chainfile = "testing/data/mcmc_output_ppE.csv";
	//std::string statfilename = "testing/data/mcmc_statistics_ppE.txt";
	//std::string autocorrfile = "testing/data/auto_corr_mcmc_dCS.csv";
	//std::string chainfile = "testing/data/mcmc_output_dCS.csv";
	//std::string statfilename = "testing/data/mcmc_statistics_dCS.txt";
	//std::string checkfile = "testing/data/mcmc_checkpoint_dCS.csv";
	std::string autocorrfile = "testing/data/auto_corr_mcmc_EdGB"+iteration+".csv";
	std::string chainfile = "testing/data/mcmc_output_EdGB"+iteration+".csv";
	std::string statfilename = "testing/data/mcmc_statistics_EdGB"+iteration+".txt";
	std::string checkfile = "testing/data/mcmc_checkpoint_EdGB"+iteration+".csv";
	std::string startcheckfile = "testing/data/mcmc_checkpoint_EdGB"+previteration+".csv";
	//std::string autocorrfile = "testing/data/auto_corr_mcmc_logflat_EdGB"+iteration+".csv";
	//std::string chainfile = "testing/data/mcmc_output_logflat_EdGB"+iteration+".csv";
	//std::string statfilename = "testing/data/mcmc_statistics_logflat_EdGB"+iteration+".txt";
	//std::string checkfile = "testing/data/mcmc_checkpoint_logflat_EdGB"+iteration+".csv";
	//std::string startcheckfile = "testing/data/mcmc_checkpoint_logflat_EdGB"+previteration+".csv";

	continue_PTMCMC_MH_GW(startcheckfile,output, dimension, n_steps, 
			swp_freq, test_lp_GW_dCS_log,numThreads, pool,show_progress,
			num_detectors, 
			data, psd,freqs, data_length,gps_time, detectors,Nmod, bppe,
			generation_method,statfilename,"","","",checkfile,true);	

	double **output_transform=(double **)malloc(sizeof(double*)*n_steps);
	for (int j =0; j<n_steps; j++)
		output_transform[j] = (double *)malloc(sizeof(double)*dimension);

	for(int j = 0; j<n_steps;j++){
		output_transform[j][0]=output[0][j][0];
		output_transform[j][1]=output[0][j][1];
		output_transform[j][2]=output[0][j][2];
		output_transform[j][3]=std::exp(output[0][j][3]);
		output_transform[j][4]=std::exp(output[0][j][4]);
		output_transform[j][5]=output[0][j][5];
		output_transform[j][6]=output[0][j][6];
		output_transform[j][7]=output[0][j][7];
		output_transform[j][8]=std::exp(output[0][j][8]);
		//output_transform[j][8]=output[0][j][8];
	}
	write_file(chainfile, output_transform, n_steps, dimension);
	//output hottest chain too
	chainfile = "testing/data/mcmc_output_logflat_EdGB_hot"+iteration+".csv";
	//chainfile = "testing/data/mcmc_output_ppE_hot.csv";
	//chainfile = "testing/data/mcmc_output_dCS_hot.csv";
	for(int j = 0; j<n_steps;j++){
		output_transform[j][0]=output[chain_N-1][j][0];
		output_transform[j][1]=output[chain_N-1][j][1];
		output_transform[j][2]=output[chain_N-1][j][2];
		output_transform[j][3]=std::exp(output[chain_N-1][j][3]);
		output_transform[j][4]=std::exp(output[chain_N-1][j][4]);
		output_transform[j][5]=output[chain_N-1][j][5];
		output_transform[j][6]=output[chain_N-1][j][6];
		output_transform[j][7]=output[chain_N-1][j][7];
		output_transform[j][8]=std::exp(output[chain_N-1][j][8]);
		//output_transform[j][8]=output[chain_N-1][j][8];
	}
	write_file(chainfile, output_transform, n_steps, dimension);

	//write_file(chainfile, output_transform, n_steps, dimension);

	int segs = 10;
	double target_corr = .01;

	//write_file_auto_corr_from_data_file_accel(autocorrfile, chainfile,dimension,n_steps,segs,target_corr);
	write_auto_corr_file_from_data_file(autocorrfile, chainfile,n_steps,dimension,segs,target_corr, numThreads,false);

	deallocate_3D_array(output, chain_N, n_steps, dimension);
	for(int i =0; i< n_steps; i++){
		free(output_transform[i]);
	}
	free(output_transform);
	delete [] detectors;
	free(data_length);
	//free_LOSC_data(data, psd,freqs, num_detectors, length);
	deallocate_2D_array(psd,num_detectors, psd_length);
	deallocate_2D_array(freqs,num_detectors, psd_length);
	for(int i =0; i<num_detectors; i++)
		free(data[i]);
	free(data);
	delete [] detector_files;
	delete [] bppe;
	//deallocate_2D_array(psd, rows, cols);
	//free(data);
	
	//int length = 2000;
	//double x[length];
	//double y[length];
	//double xlim = 10.;
	//double xstart=0;
	//double xstep = (xlim-xstart)/length;
	//for (int i =0; i<length; i++){
	//	x[i]= i*xstep;
	//	y[i]= x[i]*x[i];
	//}
	//double sum = simpsons_sum(xstep, length, y);
	//std::cout<<sum<<std::endl;


}

void test25()
{
	int N_steps = 1500;
	int dimension =2;
	int chain_N = 10;
	double ***output;
	output = allocate_3D_array( chain_N, N_steps, dimension );
	//double *initial_pos_ptr = initial_pos;
	int swp_freq = 3;
	//std::string autocorrfile = "";
	std::string iteration = "3";
	std::string previteration = "2";
	std::string autocorrfile = "testing/data/neil_auto_corr_mcmc"+iteration+".csv";
	std::string chainfile = "testing/data/neil_mcmc_output"+iteration+".csv";
	std::string statfilename = "testing/data/neil_mcmc_statistics"+iteration+".txt";
	std::string checkpointfile = "testing/data/neil_mcmc_checkpoint"+iteration+".csv";
	std::string start_checkfile = "testing/data/neil_mcmc_checkpoint"+previteration+".csv";
	
	int numThreads = 10;
	bool pool = true;
	
	continue_PTMCMC_MH(start_checkfile,output, N_steps, swp_freq, test_lp, log_neil_proj3,NULL,(void **)NULL,numThreads, pool,show_progress, statfilename,chainfile,autocorrfile, "",checkpointfile ,true);	
	std::cout<<"ENDED"<<std::endl;
}
void test24()
{
	//std::string data_file = "testing/data/mcmc_output_DFull.csv";
	//int length =50000;
	//int dim =8;
	std::string data_file = "testing/data/mcmc_output_7dim.csv";
	std::string acfile = "testing/data/auto_corr_mcmc_7dim_fft.csv";
	int length =40000;
	int dim =7;
	int segs = 50;
	double **data=allocate_2D_array(length, dim);
	read_file(data_file, data, length, dim);
	//int **ac=(int **)malloc(sizeof(int *)*dim);
	//for(int i =0 ; i<dim; i++){
	//	ac[i]=(int *)malloc(sizeof(int)*length);
	//}
	double accuracy = .01;
	int num_threads = 10;
	double wstart, wend;
	wstart = omp_get_wtime();
	//auto_corr_from_data(data, length, dim, ac, segs, accuracy, num_threads);
	write_auto_corr_file_from_data(acfile, data, length, dim, segs, accuracy, num_threads,false);
	wend = omp_get_wtime();
	std::cout<<"DONE: TIME: "<<wend-wstart<<std::endl;
	//for(int i =0 ; i<dim; i++){
	//	for(int j = 0; j<segs; j++){
	//		std::cout<<ac[i][j]<<std::endl;
	//	}
	//}
	deallocate_2D_array(data, length, dim);
	//for(int i =0; i<dim; i++)
	//	free(ac[i]);
	//free(ac);
}
void test23()
{

	int length = 750000;
	//int length = 40000;
	int dimension = 9;
	//int dimension = 7;
	std::string filename = "testing/data/mcmc_output_dCS.csv";
	//std::string filename = "testing/data/mcmc_output_7dim.csv";
	double *chain = (double *)malloc(sizeof(double)*length);
	double *ac = (double *)malloc(sizeof(double)*length);
	double **output = allocate_2D_array(length, dimension);
	read_file(filename, output, length, dimension);
	clock_t start, stop;
	for(int i =0; i< length; i++){
		chain[i] = output[i][0];
	}
	fftw_outline plan_forw;
	fftw_outline plan_rev;
	allocate_FFTW_mem_forward(&plan_forw, pow(2, std::ceil( std::log2(length))));
	allocate_FFTW_mem_reverse(&plan_rev, pow(2, std::ceil( std::log2(length))));
	start = clock();
	auto_correlation_spectral(chain, length, ac, &plan_forw, &plan_rev);
	stop = clock();
	std::cout<<"Spectral method timing: "<<(double)(stop-start)/CLOCKS_PER_SEC<<std::endl;
	deallocate_FFTW_mem(&plan_forw);
	deallocate_FFTW_mem(&plan_rev);
	for(int i =1; i<length; i++)
	{
		if(ac[i]<.01){
			std::cout<<i<<std::endl;
			std::cout<<ac[i]<<std::endl;
			break;
		}
	}
	start = clock();
	int lag_serial = auto_correlation_serial_old(chain, length);
	stop = clock();
	std::cout<<"old method timing: "<<(double)(stop-start)/CLOCKS_PER_SEC<<std::endl;
	std::cout<<lag_serial<<std::endl;

	deallocate_2D_array(output,length, dimension);
	free(chain);
	free(ac);


}
void test22()
{
	#ifdef CUDA_ENABLED
		std::cout<<"CUDA_ENABLED"<<std::endl;
	#endif
	//std::string outputfile = "testing/data/mcmc_output_dCS.csv";
	std::string outputfile = "testing/data/mcmc_output_injection.csv";
	//std::string outputfile = "testing/data/mcmc_output_DFull.csv";
	//std::string acfile = "testing/data/auto_corr_mcmc_dCS_GPU.csv";
	std::string acfile = "testing/data/auto_corr_mcmc_injection_GPU.csv";
	std::string acfilecpu = "testing/data/auto_corr_mcmc_dCS.csv";
	//std::string acfile = "testing/data/auto_corr_mcmc_DFull.csv";
	int dimension = 9;
	//int dimension = 8;
	//int N_steps = 750000;
	int N_steps = 3000000;
	//int N_steps = 20000;
	double **output = allocate_2D_array(N_steps,dimension);
	read_file(outputfile,output, N_steps, dimension);
	int segs =5;
	double target_corr = .1;
	//double **autocorr = allocate_2D_array(dimension,segs);
	clock_t start, end;
	double wstart, wend;
	//omp_set_num_threads(12);
	wstart = omp_get_wtime();
	//
	//write_auto_corr_file_from_data_file(acfilecpu, outputfile, segs, dimension, N_steps);	
	//wend = omp_get_wtime();
	//cout<<"TIMING cpu: "<<(double)(wend-wstart)<<endl;
	
	//start = clock();
	//write_file_auto_corr_from_data_accel(acfile, output,dimension,N_steps,segs,target_corr);
	//write_file_auto_corr_from_data_file_accel(acfile, outputfile,dimension,N_steps,segs,target_corr);
	wend = omp_get_wtime();
	//end = clock();
	//cout<<"TIMING gpu: "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;
	cout<<"TIMING GPU: "<<(double)(wend-wstart)<<endl;
	deallocate_2D_array(output, N_steps, dimension);
	//deallocate_2D_array(autocorr, dimension, segs);
	
}
void test21()
{
	//std::string psd_file = "testing/data/GWTC1_GW150914_PSDs.dat.txt";
	//std::string data_file = "testing/data/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//std::string psd_file = "testing/data/GWTC1_GW170729_PSDs.dat.txt";
	std::string psd_file = "/Users/sperkins/Downloads/LOSC_data/GW170608/GWTC1_GW170608_PSDs.dat.txt";
	//std::string psd_file = "/Users/sperkins/Downloads/LOSC_data/GW150914/GWTC1_GW150914_PSDs.dat.txt";
	//int rows = 8032;
	//int cols = 3;
	int datalength = 131075;
	//double **psd = allocate_2D_array(rows, cols);
	//read_LOSC_PSD_file(psd_file, psd, rows, cols);
	//double data_start_time, duration, fs;
	//int num_detectors = 2, psd_length = 8032, length;
	int num_detectors = 2, psd_length = 16064, length;
	//int num_detectors = 3, psd_length = 4016, length;
	//double gps_time = 1126259462;//TESTING -- gw150914
	//double gps_time = 1135136350.6;//TESTING -- gw151226
	//double gps_time = 1185389807.3;//TESTING -- gw170729
	double gps_time = 1180922494.5;//TESTING -- gw170608
	std::string *detectors = new std::string[num_detectors];//(std::string*)malloc(sizeof(std::string)*50*num_detectors);
	detectors[0] = "Hanford";
	detectors[1] = "Livingston";
	//detectors[2] = "Virgo";
	std::string *detector_files = new std::string[num_detectors];
	//detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1135136335-32.txt";
	detector_files[0] =  "/Users/sperkins/Downloads/LOSC_data/GW170608/H-H1_GWOSC_4KHZ_R1-1180922479-32.txt";
	detector_files[1] =  "/Users/sperkins/Downloads/LOSC_data/GW170608/L-L1_GWOSC_4KHZ_R1-1180922479-32.txt";
	//detector_files[0] =  "/Users/sperkins/Downloads/LOSC_data/GW150914/H-H1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[1] =  "/Users/sperkins/Downloads/LOSC_data/GW150914/L-L1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[2] =  "testing/data/V-V1_GWOSC_4KHZ_R1-1185389792-32.txt";
 	//double trigger_time= 1135136350.6;
 	double trigger_time = gps_time;
	double **psd = allocate_2D_array(num_detectors,psd_length);
	double **freqs = allocate_2D_array(num_detectors,psd_length);
	std::complex<double> **data = (std::complex<double> **)malloc(sizeof(std::complex<double> *)*num_detectors);
	for(int i =0; i<num_detectors; i++)
		data[i] = (std::complex<double>*)malloc(sizeof(std::complex<double>)*psd_length);

	allocate_LOSC_data(detector_files, psd_file, num_detectors, psd_length, datalength, trigger_time, data, psd, freqs);


	int *data_length= (int*)malloc(sizeof(int)*num_detectors);
	data_length[0] =psd_length;
	data_length[1] =psd_length;
	//data_length[2] =psd_length;

	//#########################################################
	//mcmc options
	int dimension = 9;
	double initial_pos[dimension]={0., log(10), .2, .0,- .0,.1, .1,.1,.1};

	double *seeding_var = NULL;
	int n_steps = 200;
	int chain_N= 8;
	double ***output;
	output = allocate_3D_array( chain_N, n_steps, dimension );
	int swp_freq = 10;
	double chain_temps[chain_N];
	chain_temps[0]=1.;
	double c = 1.2;
	for(int i =1; i < chain_N;  i ++)
		chain_temps[i] = c*chain_temps[i-1];
	
	int Nmod = 0;
	int *bppe = NULL;
	int numThreads = 5;
	bool pool = true;
	//#########################################################
	//gw options
	std::string generation_method = "IMRPhenomPv2";
	
	
	std::string autocorrfile = "testing/data/auto_corr_mcmc_Pv2.csv";
	//std::string autocorrfile = "";
	std::string chainfile = "testing/data/mcmc_output_Pv2.csv";
	std::string statfilename = "testing/data/mcmc_statistics_Pv2.txt";
	std::string checkpointfile = "testing/data/mcmc_checkpoint_Pv2.csv";

	PTMCMC_MH_GW(output, dimension, n_steps, chain_N, initial_pos,seeding_var,chain_temps, 
			swp_freq, test_lp_GW_DFull,numThreads, pool,show_progress,
			num_detectors, 
			data, psd,freqs, data_length,gps_time, detectors,Nmod, bppe,
			generation_method,statfilename,"",autocorrfile, "",checkpointfile);	

	double **output_transform=(double **)malloc(sizeof(double*)*n_steps);
	for (int j =0; j<n_steps; j++)
		output_transform[j] = (double *)malloc(sizeof(double)*dimension);

	for(int j = 0; j<n_steps;j++){
			output_transform[j][0]=output[0][j][0];
			output_transform[j][1]=std::exp(output[0][j][1]);
			output_transform[j][2]=output[0][j][2];
			output_transform[j][3]=output[0][j][3];
			output_transform[j][4]=output[0][j][4];
			output_transform[j][5]=output[0][j][5];
			output_transform[j][6]=output[0][j][6];
			output_transform[j][7]=output[0][j][7];
			output_transform[j][8]=output[0][j][8];
	}
	write_file(chainfile, output_transform, n_steps, dimension);
	//output hottest chain too
	chainfile = "testing/data/mcmc_output_Pv2_hot.csv";
	for(int j = 0; j<n_steps;j++){
			output_transform[j][0]=output[chain_N-1][j][0];
			output_transform[j][1]=std::exp(output[chain_N-1][j][1]);
			output_transform[j][2]=output[chain_N-1][j][2];
			output_transform[j][3]=output[chain_N-1][j][3];
			output_transform[j][4]=output[chain_N-1][j][4];
			output_transform[j][5]=output[chain_N-1][j][5];
			output_transform[j][6]=output[chain_N-1][j][6];
			output_transform[j][7]=output[chain_N-1][j][7];
			output_transform[j][8]=output[chain_N-1][j][8];
	}
	write_file(chainfile, output_transform, n_steps, dimension);

	deallocate_3D_array(output, chain_N, n_steps, dimension);
	for(int i =0; i< n_steps; i++){
		free(output_transform[i]);
	}
	free(output_transform);
	delete [] detectors;
	free(data_length);
	//free_LOSC_data(data, psd,freqs, num_detectors, length);
	deallocate_2D_array(psd,num_detectors, psd_length);
	deallocate_2D_array(freqs,num_detectors, psd_length);
	for(int i =0; i<num_detectors; i++)
		free(data[i]);
	free(data);
	delete [] detector_files;
	//deallocate_2D_array(psd, rows, cols);
	//free(data);
	
	//int length = 2000;
	//double x[length];
	//double y[length];
	//double xlim = 10.;
	//double xstart=0;
	//double xstep = (xlim-xstart)/length;
	//for (int i =0; i<length; i++){
	//	x[i]= i*xstep;
	//	y[i]= x[i]*x[i];
	//}
	//double sum = simpsons_sum(xstep, length, y);
	//std::cout<<sum<<std::endl;

}
void test20()
{
	//std::string psd_file = "testing/data/GWTC1_GW150914_PSDs.dat.txt";
	std::string psd_file = "testing/data/GWTC1_GW151226_PSDs.dat.txt";
	//std::string data_file = "testing/data/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//std::string psd_file = "testing/data/GWTC1_GW170729_PSDs.dat.txt";
	//int rows = 8032;
	//int cols = 3;
	int datalength = 131075;
	//double **psd = allocate_2D_array(rows, cols);
	//read_LOSC_PSD_file(psd_file, psd, rows, cols);
	//double data_start_time, duration, fs;
	int num_detectors = 2, psd_length = 8032, length;
	//int num_detectors = 3, psd_length = 4016, length;
	double gps_time = 1135136350.6;//TESTING -- gw151226
	//double gps_time = 1126259462;//TESTING -- gw150914
	//double gps_time = 1185389807.3;//TESTING -- gw170729
	std::string *detectors = new std::string[num_detectors];//(std::string*)malloc(sizeof(std::string)*50*num_detectors);
	detectors[0] = "Hanford";
	detectors[1] = "Livingston";
	//detectors[2] = "Virgo";
	std::string *detector_files = new std::string[num_detectors];
	detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1185389792-32.txt";
	//detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1185389792-32.txt";
	//detector_files[2] =  "testing/data/V-V1_GWOSC_4KHZ_R1-1185389792-32.txt";
 	//double trigger_time= 1135136350.6;
 	double trigger_time = gps_time;
	double **psd = allocate_2D_array(num_detectors,psd_length);
	double **freqs = allocate_2D_array(num_detectors,psd_length);
	std::complex<double> **data = (std::complex<double> **)malloc(sizeof(std::complex<double> *)*num_detectors);
	for(int i =0; i<num_detectors; i++)
		data[i] = (std::complex<double>*)malloc(sizeof(std::complex<double>)*psd_length);

	allocate_LOSC_data(detector_files, psd_file, num_detectors, psd_length, datalength, trigger_time, data, psd, freqs);


	int *data_length= (int*)malloc(sizeof(int)*num_detectors);
	data_length[0] =psd_length;
	data_length[1] =psd_length;
	//data_length[2] =psd_length;

	//#########################################################
	//mcmc options
	int dimension = 9;
	//double initial_pos[dimension]={.0, 1, 0.,log(300),log(10), .2,- .0,-.0, log(pow(MPC_SEC,4)*pow(50000,4))};
	//double initial_pos[dimension]={.0, 1, 0.,log(300),log(10), .2,- .0,-.0, log(pow(MPC_SEC,4)*pow(5,4))};
	//double initial_pos[dimension]={.0, 1, 0.,log(300),log(10), .2,- .0,-.0, -5};
	//double initial_pos[dimension]={.0, 1, 0.,log(300),log(10), .2,- .0,-.0,4* log(50000/(3e8))};
	double initial_pos[dimension]={.9, 2, 1.,log(400),log(10), .24,- .0,-.0,-40};
	//double initial_pos[dimension]={.9, 2, -1.,log(300),log(10), .24,- .0,-.0,-0};
	double *seeding_var = NULL;
	int n_steps = 15000;
	//int n_steps = 122000;
	int chain_N=10 ;
	double ***output;
	output = allocate_3D_array( chain_N, n_steps, dimension );
	int swp_freq = 5;
	double chain_temps[chain_N];
	chain_temps[0]=1.;
	double c = 1.2;
	for(int i =1; i < chain_N;  i ++)
		chain_temps[i] = c*chain_temps[i-1];
	
	int Nmod = 1;
	int *bppe = new int[Nmod];
	bppe[0] = -1;
	int numThreads = 12;
	bool pool = true;
	//#########################################################
	//gw options
	//std::string generation_method = "ppE_IMRPhenomD_Inspiral";
	std::string generation_method = "EdGB_IMRPhenomD_log";
	//std::string generation_method = "dCS_IMRPhenomD_log";
	
	
	//std::string autocorrfile = "";
	//std::string autocorrfile = "testing/data/auto_corr_mcmc_ppE.csv";
	//std::string chainfile = "testing/data/mcmc_output_ppE.csv";
	//std::string statfilename = "testing/data/mcmc_statistics_ppE.txt";
	//std::string autocorrfile = "testing/data/auto_corr_mcmc_dCS.csv";
	//std::string chainfile = "testing/data/mcmc_output_dCS.csv";
	//std::string statfilename = "testing/data/mcmc_statistics_dCS.txt";
	//std::string checkfile = "testing/data/mcmc_checkpoint_dCS.csv";
	std::string autocorrfile = "testing/data/auto_corr_mcmc_EdGB1.csv";
	std::string chainfile = "testing/data/mcmc_output_EdGB1.csv";
	std::string statfilename = "testing/data/mcmc_statistics_EdGB1.txt";
	std::string checkfile = "testing/data/mcmc_checkpoint_EdGB1.csv";
	//std::string autocorrfile = "testing/data/auto_corr_mcmc_logflat_EdGB1.csv";
	//std::string chainfile = "testing/data/mcmc_output_logflat_EdGB1.csv";
	//std::string statfilename = "testing/data/mcmc_statistics_logflat_EdGB1.txt";
	//std::string checkfile = "testing/data/mcmc_checkpoint_logflat_EdGB1.csv";
	//std::string autocorrfile = "testing/data/auto_corr_mcmc_EdGB.csv";
	//std::string chainfile = "testing/data/mcmc_output_EdGB.csv";
	//std::string statfilename = "testing/data/mcmc_statistics_EdGB.txt";
	//std::string checkfile = "testing/data/mcmc_checkpoint_EdGB.csv";

	PTMCMC_MH_GW(output, dimension, n_steps, chain_N, initial_pos,seeding_var,chain_temps, 
			swp_freq, test_lp_GW_dCS_log,numThreads, pool,show_progress,
			num_detectors, 
			data, psd,freqs, data_length,gps_time, detectors,Nmod, bppe,
			generation_method,statfilename,"","","",checkfile);	

	double **output_transform=(double **)malloc(sizeof(double*)*n_steps);
	for (int j =0; j<n_steps; j++)
		output_transform[j] = (double *)malloc(sizeof(double)*dimension);

	for(int j = 0; j<n_steps;j++){
		output_transform[j][0]=output[0][j][0];
		output_transform[j][1]=output[0][j][1];
		output_transform[j][2]=output[0][j][2];
		output_transform[j][3]=std::exp(output[0][j][3]);
		output_transform[j][4]=std::exp(output[0][j][4]);
		output_transform[j][5]=output[0][j][5];
		output_transform[j][6]=output[0][j][6];
		output_transform[j][7]=output[0][j][7];
		output_transform[j][8]=std::exp(output[0][j][8]);
		//output_transform[j][8]=output[0][j][8];
	}
	write_file(chainfile, output_transform, n_steps, dimension);
	//output hottest chain too
	chainfile = "testing/data/mcmc_output_EdGB_hot1.csv";
	//chainfile = "testing/data/mcmc_output_logflat_EdGB_hot1.csv";
	//chainfile = "testing/data/mcmc_output_EdGB_hot.csv";
	//chainfile = "testing/data/mcmc_output_ppE_hot.csv";
	//chainfile = "testing/data/mcmc_output_dCS_hot.csv";
	for(int j = 0; j<n_steps;j++){
		output_transform[j][0]=output[chain_N-1][j][0];
		output_transform[j][1]=output[chain_N-1][j][1];
		output_transform[j][2]=output[chain_N-1][j][2];
		output_transform[j][3]=std::exp(output[chain_N-1][j][3]);
		output_transform[j][4]=std::exp(output[chain_N-1][j][4]);
		output_transform[j][5]=output[chain_N-1][j][5];
		output_transform[j][6]=output[chain_N-1][j][6];
		output_transform[j][7]=output[chain_N-1][j][7];
		output_transform[j][8]=std::exp(output[chain_N-1][j][8]);
		//output_transform[j][8]=output[chain_N-1][j][8];
	}
	write_file(chainfile, output_transform, n_steps, dimension);

	write_file(chainfile, output_transform, n_steps, dimension);

	int segs = 10;
	double target_corr = .01;

	//write_file_auto_corr_from_data_file_accel(autocorrfile, chainfile,dimension,n_steps,segs,target_corr);
	write_auto_corr_file_from_data_file(autocorrfile, chainfile,n_steps,dimension,segs,target_corr, numThreads,false);

	deallocate_3D_array(output, chain_N, n_steps, dimension);
	for(int i =0; i< n_steps; i++){
		free(output_transform[i]);
	}
	free(output_transform);
	delete [] detectors;
	free(data_length);
	//free_LOSC_data(data, psd,freqs, num_detectors, length);
	deallocate_2D_array(psd,num_detectors, psd_length);
	deallocate_2D_array(freqs,num_detectors, psd_length);
	for(int i =0; i<num_detectors; i++)
		free(data[i]);
	free(data);
	delete [] detector_files;
	delete [] bppe;
	//deallocate_2D_array(psd, rows, cols);
	//free(data);
	
	//int length = 2000;
	//double x[length];
	//double y[length];
	//double xlim = 10.;
	//double xstart=0;
	//double xstep = (xlim-xstart)/length;
	//for (int i =0; i<length; i++){
	//	x[i]= i*xstep;
	//	y[i]= x[i]*x[i];
	//}
	//double sum = simpsons_sum(xstep, length, y);
	//std::cout<<sum<<std::endl;

}
void test19()
{
	gen_params params;
	double RA=1, DEC=.2;
	double gps_time = 1185389807.3;//TESTING -- gw170729
	params.mass1 = 80;
	params.mass2 = 40;
	params.spin1[0] = 0;
	params.spin1[1] = 0;
	params.spin2[0] = 0;
	params.spin2[1] = 0;
	params.spin1[2] = .1;
	params.spin2[2] = .1;
	params.Nmod = 1;
	params.betappe= new double[params.Nmod];
	params.betappe[0] = -1;
	params.sky_average = false;
	params.NSflag1 = false;
	params.NSflag2 = false;
	params.incl_angle = M_PI/3.;
	params.theta = 0;
	params.phi = 0;
	params.Luminosity_Distance = 400;
	params.phic = 0;
	params.tc = 0;
	std::string gen_method = "MCMC_dCS_IMRPhenomD_log_Full";
	std::string detector = "Hanford";
	
	int length = 4000;
	double freqs[length];
	double fhigh = 1024;
	double flow = 10;
	double fstep = (fhigh-flow)/length;
	for(int i =0; i<length; i++){
		freqs[i] = flow + i*fstep;
	}
	double noise[length];
	populate_noise(freqs,"Hanford_O1_fitted", noise,length);
	for (int i =0; i<length;i++){
		noise[i] = noise[i]*noise[i];
	}
	int dim = 9;
	double **fishout = allocate_2D_array(dim, dim);
	std::complex<double> response[length];
	std::complex<double> responsedcs[length];
	double tc=2;
	params.tc = tc;
	double phic;
	double phicmax, tcmax;
	fftw_outline plan;
	allocate_FFTW_mem_forward(&plan, length);

	//fisher(freqs, length, gen_method, detector,fishout,dim, &params,NULL, NULL, noise);
	fourier_detector_response(freqs, length, response, detector, "IMRPhenomD",&params);
	fourier_detector_response(freqs, length, responsedcs, detector, "dCS_IMRPhenomD_log",&params);


	celestial_horizon_transform(RA, DEC, gps_time, "Hanford", &params.phi, &params.theta);
	double llmax = maximized_coal_Log_Likelihood(response, noise, freqs, length, &params,
			"Hanford", "dCS_IMRPhenomD_log",&plan, &tcmax,&phicmax);
	std::cout<<tcmax<<" "<<phicmax<<std::endl;
	double temptheta = params.theta;
	celestial_horizon_transform(RA, DEC, gps_time, "Livingston", &params.phi, &params.theta);
	params.tc = tc+ DTOA(temptheta, params.theta, "Hanford","Livingston"); 
	double ll = Log_Likelihood(response, noise, freqs, length, &params,
			"Livingston", "dCS_IMRPhenomD_log", &plan);
	std::cout.precision(15);
	std::cout<<llmax<<std::endl;
	std::cout<<ll<<std::endl;
	std::cout<<ll+llmax<<std::endl;
	for(int i =0; i<length/4;i++){
	//for(int i =0; i<dim;i++){
		for(int j =0; j<dim; j++){
			//std::cout<<i<<" "<<j<<" "<<fishout[i][j]<<std::endl;
			//std::cout<<i<<" "<<j<<" "<<fishout[i][j]<<std::endl;
		}
		//std::cout<<responsedcs[i]-response[i]<<std::endl;
	}
	deallocate_2D_array(fishout, dim , dim);
	delete [] params.betappe;
	deallocate_FFTW_mem(&plan);
}
void test18()
{
	//std::string psd_file = "testing/data/GWTC1_GW150914_PSDs.dat.txt";
	std::string psd_file = "testing/data/GWTC1_GW151226_PSDs.dat.txt";
	//std::string data_file = "testing/data/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//std::string psd_file = "testing/data/GWTC1_GW170729_PSDs.dat.txt";
	//int rows = 8032;
	//int cols = 3;
	int datalength = 131075;
	//double **psd = allocate_2D_array(rows, cols);
	//read_LOSC_PSD_file(psd_file, psd, rows, cols);
	//double data_start_time, duration, fs;
	int num_detectors = 2, psd_length = 8032, length;
	//int num_detectors = 3, psd_length = 4016, length;
	double gps_time = 1135136350.6;//TESTING -- gw151226
	//double gps_time = 1126259462;//TESTING -- gw150914
	//double gps_time = 1185389807.3;//TESTING -- gw170729
	std::string *detectors = new std::string[num_detectors];//(std::string*)malloc(sizeof(std::string)*50*num_detectors);
	detectors[0] = "Hanford";
	detectors[1] = "Livingston";
	//detectors[2] = "Virgo";
	std::string *detector_files = new std::string[num_detectors];
	detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1185389792-32.txt";
	//detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1185389792-32.txt";
	//detector_files[2] =  "testing/data/V-V1_GWOSC_4KHZ_R1-1185389792-32.txt";
 	//double trigger_time= 1135136350.6;
 	double trigger_time = gps_time;
	double **psd = allocate_2D_array(num_detectors,psd_length);
	double **freqs = allocate_2D_array(num_detectors,psd_length);
	std::complex<double> **data = (std::complex<double> **)malloc(sizeof(std::complex<double> *)*num_detectors);
	for(int i =0; i<num_detectors; i++)
		data[i] = (std::complex<double>*)malloc(sizeof(std::complex<double>)*psd_length);

	allocate_LOSC_data(detector_files, psd_file, num_detectors, psd_length, datalength, trigger_time, data, psd, freqs);


	int *data_length= (int*)malloc(sizeof(int)*num_detectors);
	data_length[0] =psd_length;
	data_length[1] =psd_length;
	//data_length[2] =psd_length;

	//#########################################################
	//mcmc options
	int dimension = 9;
	//double initial_pos[dimension]={.0, 1, 0.,log(300),log(10), .2,- .0,-.0, log(pow(MPC_SEC,4)*pow(50000,4))};
	//double initial_pos[dimension]={.0, 1, 0.,log(300),log(10), .2,- .0,-.0, log(pow(MPC_SEC,4)*pow(5,4))};
	//double initial_pos[dimension]={.0, 1, 0.,log(300),log(10), .2,- .0,-.0, -5};
	//double initial_pos[dimension]={.0, 1, 0.,log(300),log(10), .2,- .0,-.0,4* log(50000/(3e8))};
	//double initial_pos[dimension]={.9, 2, 1.,log(300),log(10), .24,- .0,-.0,-40};
	//double initial_pos[dimension]={.9, 2, 1.,std::log(300),std::log(10), .24,- .0,-.0,7.7e-18};
	double initial_pos[dimension]={.9, 2, 1.,std::log(500),std::log(9), .22,- .0,.4,1};
	double *seeding_var = NULL;
	int n_steps = 100000;
	int chain_N=8 ;
	double ***output;
	output = allocate_3D_array( chain_N, n_steps, dimension );
	int swp_freq = 5;
	double chain_temps[chain_N];
	double c = 1.4;
	chain_temps[0]=1.;
	chain_temps[1] = c;
	chain_temps[2] = c*c;
	chain_temps[3] = c*c*c;
	chain_temps[4] = 1;
	chain_temps[5] = c;
	chain_temps[6] = c*c;
	chain_temps[7] = c*c*c;
	//chain_temps[1] = 1.1;
	//chain_temps[2] = 1.2;
	//chain_temps[3] = 1.3;
	//chain_temps[4] = 1.4;
	//chain_temps[5] = 1.5;
	//chain_temps[6] = 1.6;
	//chain_temps[7] = 1.7;
	//for(int i =1; i < chain_N;  i ++)
	//	chain_temps[i] = c*chain_temps[i-1];
	//double c = 1./std::sqrt(dimension);//optimum spacing of inverse temperature
	//for(int i =1; i < chain_N;  i ++){
	//	chain_temps[i] = chain_temps[i-1]/(1-chain_temps[i-1]*c);
	//	std::cout<<chain_temps[i]<<std::endl;
	//}
	
	int Nmod = 1;
	int *bppe = NULL;
	int numThreads = 4;
	bool pool = false;
	//#########################################################
	//gw options
	//std::string generation_method = "dCS_IMRPhenomD_log";
	//std::string generation_method = "EdGB_IMRPhenomD_log";
	std::string generation_method = "dCS_IMRPhenomD_root_alpha";
	//std::string generation_method = "EdGB_IMRPhenomD_root_alpha";
	
	
	//std::string autocorrfile = "";
	std::string autocorrfile = "testing/data/auto_corr_mcmc_dCS.csv";
	std::string chainfile = "testing/data/mcmc_output_dCS.csv";
	std::string statfilename = "testing/data/mcmc_statistics_dCS.txt";
	std::string checkfile = "testing/data/mcmc_checkpoint_dCS.csv";
	//std::string autocorrfile = "testing/data/auto_corr_mcmc_EdGB.csv";
	//std::string chainfile = "testing/data/mcmc_output_EdGB.csv";
	//std::string statfilename = "testing/data/mcmc_statistics_EdGB.txt";
	//std::string checkfile = "testing/data/mcmc_checkpoint_EdGB.csv";

	PTMCMC_MH_GW(output, dimension, n_steps, chain_N, initial_pos,seeding_var,chain_temps, 
			swp_freq, test_lp_GW_dCS_root_alpha,numThreads, pool,show_progress,
			num_detectors, 
			data, psd,freqs, data_length,gps_time, detectors,Nmod, bppe,
			generation_method,statfilename,"",autocorrfile, "",checkfile);	

	double **output_transform=(double **)malloc(sizeof(double*)*n_steps);
	for (int j =0; j<n_steps; j++)
		output_transform[j] = (double *)malloc(sizeof(double)*dimension);

	for(int j = 0; j<n_steps;j++){
		output_transform[j][0]=output[0][j][0];
		output_transform[j][1]=output[0][j][1];
		output_transform[j][2]=output[0][j][2];
		output_transform[j][3]=std::exp(output[0][j][3]);
		output_transform[j][4]=std::exp(output[0][j][4]);
		output_transform[j][5]=output[0][j][5];
		output_transform[j][6]=output[0][j][6];
		output_transform[j][7]=output[0][j][7];
		//output_transform[j][8]=std::exp(output[0][j][8]);
		output_transform[j][8]=output[0][j][8];
	}
	write_file(chainfile, output_transform, n_steps, dimension);
	std::string chain2 = "testing/data/mcmc_output_dCS_2.csv";
	write_file(chain2, output[4], n_steps, dimension);
	//output hottest chain too
	chainfile = "testing/data/mcmc_output_EdGB_hot.csv";
	//chainfile = "testing/data/mcmc_output_dCS_hot.csv";
	for(int j = 0; j<n_steps;j++){
		output_transform[j][0]=output[chain_N-1][j][0];
		output_transform[j][1]=output[chain_N-1][j][1];
		output_transform[j][2]=output[chain_N-1][j][2];
		output_transform[j][3]=std::exp(output[chain_N-1][j][3]);
		output_transform[j][4]=std::exp(output[chain_N-1][j][4]);
		output_transform[j][5]=output[chain_N-1][j][5];
		output_transform[j][6]=output[chain_N-1][j][6];
		output_transform[j][7]=output[chain_N-1][j][7];
		//output_transform[j][8]=std::exp(output[chain_N-1][j][8]);
		output_transform[j][8]=output[chain_N-1][j][8];
	}
	write_file(chainfile, output_transform, n_steps, dimension);

	deallocate_3D_array(output, chain_N, n_steps, dimension);
	for(int i =0; i< n_steps; i++){
		free(output_transform[i]);
	}
	free(output_transform);
	delete [] detectors;
	free(data_length);
	//free_LOSC_data(data, psd,freqs, num_detectors, length);
	deallocate_2D_array(psd,num_detectors, psd_length);
	deallocate_2D_array(freqs,num_detectors, psd_length);
	for(int i =0; i<num_detectors; i++)
		free(data[i]);
	free(data);
	delete [] detector_files;
	//deallocate_2D_array(psd, rows, cols);
	//free(data);
	
	//int length = 2000;
	//double x[length];
	//double y[length];
	//double xlim = 10.;
	//double xstart=0;
	//double xstep = (xlim-xstart)/length;
	//for (int i =0; i<length; i++){
	//	x[i]= i*xstep;
	//	y[i]= x[i]*x[i];
	//}
	//double sum = simpsons_sum(xstep, length, y);
	//std::cout<<sum<<std::endl;

}
void test17()
{
	std::cout<<Z_from_DL(1000.,"PLANCK13")<<std::endl;
	std::cout<<Z_from_DL(10000.,"PLANCK13")<<std::endl;
	std::cout<<Z_from_DL(100.,"PLANCK13")<<std::endl;
	std::cout<<DL_from_Z(10.,"PLANCK13")<<std::endl;
	std::cout<<DL_from_Z(1.,"PLANCK13")<<std::endl;
	std::cout<<DL_from_Z(.1,"PLANCK13")<<std::endl;
	//clock_t start7,end7;
	//int threads =10;
	//start7 = clock();
	//gsl_spline **z_d_spline = (gsl_spline **)malloc(sizeof(gsl_spline*)*threads);	
	//gsl_interp_accel **z_d_accel= (gsl_interp_accel **)malloc(sizeof(gsl_interp_accel *)*threads);	
	//for(int i =0; i<threads; i++){
	//	initiate_LumD_Z_interp(&z_d_accel[i], &z_d_spline[i]);
	//}
	//int iterations = 100;
	//omp_set_num_threads(threads);
	//#pragma omp parallel for 
	//for (int j =0; j<iterations; j++){
	//	gsl_spline *z_d_spline;
	//	gsl_interp_accel *z_d_accel;
	//	initiate_LumD_Z_interp(&z_d_accel, &z_d_spline);
	//	std::cout<<j<<std::endl;
	//	std::cout<<Z_from_DL_interp(j+10, z_d_accel, z_d_spline)<<std::endl;;
	//	free_LumD_Z_interp(&z_d_accel, &z_d_spline);
	//}

	////for (int i =0; i<threads;i++)
	////	free_LumD_Z_interp(&z_d_accel[i], &z_d_spline[i]);
	//end7 = clock();

	//cout<<"TIMING interp: "<<(double)(end7-start7)/CLOCKS_PER_SEC<<endl;
}

void test16()
{
	//std::string psd_file = "testing/data/GWTC1_GW150914_PSDs.dat.txt";
	//std::string data_file = "testing/data/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//std::string psd_file = "testing/data/GWTC1_GW170729_PSDs.dat.txt";
	std::string psd_file = "testing/data/GWTC1_GW151226_PSDs.dat.txt";
	//std::string psd_file = "/Users/sperkins/Downloads/LOSC_data/GW170608/GWTC1_GW170608_PSDs.dat.txt";
	//std::string psd_file = "/Users/sperkins/Downloads/LOSC_data/GW150914/GWTC1_GW150914_PSDs.dat.txt";
	//int rows = 8032;
	//int cols = 3;
	int datalength = 131075;
	//double **psd = allocate_2D_array(rows, cols);
	//read_LOSC_PSD_file(psd_file, psd, rows, cols);
	//double data_start_time, duration, fs;
	int num_detectors = 2, psd_length = 8032, length;
	//int num_detectors = 2, psd_length = 16064, length;
	//int num_detectors = 3, psd_length = 4016, length;
	//double gps_time = 1126259462;//TESTING -- gw150914
	double gps_time = 1135136350.6;//TESTING -- gw151226
	//double gps_time = 1185389807.3;//TESTING -- gw170729
	//double gps_time = 1180922494.5;//TESTING -- gw170608
	std::string *detectors = new std::string[num_detectors];//(std::string*)malloc(sizeof(std::string)*50*num_detectors);
	detectors[0] = "Hanford";
	detectors[1] = "Livingston";
	//detectors[2] = "Virgo";
	std::string *detector_files = new std::string[num_detectors];
	detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//detector_files[0] =  "/Users/sperkins/Downloads/LOSC_data/GW170608/H-H1_GWOSC_4KHZ_R1-1180922479-32.txt";
	//detector_files[1] =  "/Users/sperkins/Downloads/LOSC_data/GW170608/L-L1_GWOSC_4KHZ_R1-1180922479-32.txt";
	//detector_files[0] =  "/Users/sperkins/Downloads/LOSC_data/GW150914/H-H1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[1] =  "/Users/sperkins/Downloads/LOSC_data/GW150914/L-L1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[2] =  "testing/data/V-V1_GWOSC_4KHZ_R1-1185389792-32.txt";
 	//double trigger_time= 1135136350.6;
 	double trigger_time = gps_time;
	double **psd = allocate_2D_array(num_detectors,psd_length);
	double **freqs = allocate_2D_array(num_detectors,psd_length);
	std::complex<double> **data = (std::complex<double> **)malloc(sizeof(std::complex<double> *)*num_detectors);
	for(int i =0; i<num_detectors; i++)
		data[i] = (std::complex<double>*)malloc(sizeof(std::complex<double>)*psd_length);

	allocate_LOSC_data(detector_files, psd_file, num_detectors, psd_length, datalength, trigger_time, data, psd, freqs);


	int *data_length= (int*)malloc(sizeof(int)*num_detectors);
	data_length[0] =psd_length;
	data_length[1] =psd_length;
	//data_length[2] =psd_length;

	//#########################################################
	//mcmc options
	int dimension = 8;
	//double initial_pos[dimension]={.3, 2., -0.2,log(400),log(30), .24,- .0,-.0};
	double initial_pos[dimension]={.0, 1, 0.,log(300),log(10), .2,- .0,-.0};
	//double initial_pos[dimension]={-.9, 2, -1.2,log(410),log(30), .24,-.4,.3};
	//double initial_pos[dimension]={-.0, 0, -0,log(500),log(50), .2,-.0,.0};
	//double initial_pos[dimension]={-.99, 2, -1.2,log(410),log(30.78), .24,-.4,.3};
	double *seeding_var = NULL;
	int n_steps = 100000;
	int chain_N= 10;
	double ***output;
	output = allocate_3D_array( chain_N, n_steps, dimension );
	int swp_freq = 5;
	double chain_temps[chain_N];
	chain_temps[0]=1.;
	double c = 1.2;
	for(int i =1; i < chain_N;  i ++)
		chain_temps[i] = c*chain_temps[i-1];
	
	int Nmod = 0;
	int *bppe = NULL;
	int numThreads = 4;
	bool pool = true;
	//bool pool = false;
	//#########################################################
	//gw options
	std::string generation_method = "IMRPhenomD";
	
	
	std::string autocorrfile = "testing/data/auto_corr_mcmc_DFull1.csv";
	//std::string autocorrfile = "";
	std::string chainfile = "testing/data/mcmc_output_DFull1.csv";
	std::string statfilename = "testing/data/mcmc_statistics_DFull1.txt";
	std::string checkfile = "testing/data/mcmc_checkpoint_DFull1.csv";
	//std::string autocorrfile = "testing/data/auto_corr_mcmc_DFull.csv";
	////std::string autocorrfile = "";
	//std::string chainfile = "testing/data/mcmc_output_DFull.csv";
	//std::string statfilename = "testing/data/mcmc_statistics_DFull.txt";
	//std::string checkfile = "testing/data/mcmc_checkpoint_DFull.csv";

	PTMCMC_MH_GW(output, dimension, n_steps, chain_N, initial_pos,seeding_var,chain_temps, 
			swp_freq, test_lp_GW_DFull,numThreads, pool,show_progress,
			num_detectors, 
			data, psd,freqs, data_length,gps_time, detectors,Nmod, bppe,
			generation_method,statfilename,"",autocorrfile, "",checkfile);	
	std::cout<<"ended"<<std::endl;

	double **output_transform=(double **)malloc(sizeof(double*)*n_steps);
	for (int j =0; j<n_steps; j++)
		output_transform[j] = (double *)malloc(sizeof(double)*dimension);

	for(int j = 0; j<n_steps;j++){
			output_transform[j][0]=output[0][j][0];
			output_transform[j][1]=output[0][j][1];
			output_transform[j][2]=output[0][j][2];
			output_transform[j][3]=std::exp(output[0][j][3]);
			output_transform[j][4]=std::exp(output[0][j][4]);
			output_transform[j][5]=output[0][j][5];
			output_transform[j][6]=output[0][j][6];
			output_transform[j][7]=output[0][j][7];
	}
	write_file(chainfile, output_transform, n_steps, dimension);
	//output hottest chain too
	chainfile = "testing/data/mcmc_output_DFull_hot1.csv";
	//chainfile = "testing/data/mcmc_output_DFull_hot.csv";
	for(int j = 0; j<n_steps;j++){
			output_transform[j][0]=output[chain_N-1][j][0];
			output_transform[j][1]=output[chain_N-1][j][1];
			output_transform[j][2]=output[chain_N-1][j][2];
			output_transform[j][3]=std::exp(output[chain_N-1][j][3]);
			output_transform[j][4]=std::exp(output[chain_N-1][j][4]);
			output_transform[j][5]=output[chain_N-1][j][5];
			output_transform[j][6]=output[chain_N-1][j][6];
			output_transform[j][7]=output[chain_N-1][j][7];
	}
	write_file(chainfile, output_transform, n_steps, dimension);

	deallocate_3D_array(output, chain_N, n_steps, dimension);
	for(int i =0; i< n_steps; i++){
		free(output_transform[i]);
	}
	free(output_transform);
	delete [] detectors;
	free(data_length);
	//free_LOSC_data(data, psd,freqs, num_detectors, length);
	deallocate_2D_array(psd,num_detectors, psd_length);
	deallocate_2D_array(freqs,num_detectors, psd_length);
	for(int i =0; i<num_detectors; i++)
		free(data[i]);
	free(data);
	delete [] detector_files;
	//deallocate_2D_array(psd, rows, cols);
	//free(data);
	
	//int length = 2000;
	//double x[length];
	//double y[length];
	//double xlim = 10.;
	//double xstart=0;
	//double xstep = (xlim-xstart)/length;
	//for (int i =0; i<length; i++){
	//	x[i]= i*xstep;
	//	y[i]= x[i]*x[i];
	//}
	//double sum = simpsons_sum(xstep, length, y);
	//std::cout<<sum<<std::endl;

}
void test15()
{
	int length = 4000;
	//int num_detectors =3;
	int num_detectors =2;
	std::string *detectors = new std::string[num_detectors];//(std::string*)malloc(sizeof(std::string)*50*num_detectors);
	detectors[0] = "Hanford";
	detectors[1] = "Livingston";
	//detectors[2] = "Virgo";

	double gps_time = 1135136350.6;//TESTING -- gw151226

	std::complex<double> **data= (std::complex<double>**)malloc(
			sizeof(std::complex<double>*)*num_detectors);
	double **psd = (double **)malloc(sizeof(double *)*num_detectors);
	double **frequencies = (double **)malloc(sizeof(double *)*num_detectors);
	int *data_length= (int*)malloc(sizeof(int)*num_detectors);
	data_length[0] =length;
	data_length[1] =length;
	//data_length[2] =length;

	for (int i =0; i<num_detectors; i++){
		data[i] = (std::complex<double> *)malloc(
			sizeof(std::complex<double>)*data_length[i]);
		
		psd[i] = (double *)malloc(sizeof(double)*data_length[i]);
		frequencies[i] = (double *)malloc(sizeof(double)*data_length[i]);
	}
	//#########################################################
	//Make trial data
	gen_params params;
	//double RA = 5.;
	//double DEC = 1.;
	double RA = 2;
	double DEC = .7;
	double chirpm = 9.71;
	double eta =.2;
	params.mass1 = calculate_mass1(chirpm,eta);
	params.mass2 = calculate_mass2(chirpm,eta);
	complex<double> waveformout[length];
	params.spin1[0] = 0;
	params.spin1[1] = 0;
	params.spin1[2] = .43;
	params.spin2[0] = 0;
	params.spin2[1] = 0;
	params.spin2[2] = .0;
	params.phic = .0;
	double tc = 2;
	params.Luminosity_Distance = 125.;
	params.NSflag1 = false;
	params.NSflag2 = false;
	params.incl_angle = 0.;//M_PI/3.;
	params.sky_average=false;
	params.bppe = new int[1];
	params.betappe = new double[1];
	params.bppe[0] = -1;
	params.betappe[0] = (pow(0./(3e5),4));
	//params.betappe[0] = log(pow(100./(3e5),4));
	params.Nmod = 1;
	//std::string injection_method = "dCS_IMRPhenomD";
	std::string injection_method = "IMRPhenomD";
	//params.f_ref = 30.5011;
	//params.phiRef =58.944425/2.;
	
	//############################################################
	double fhigh =500;
	double flow =10;
	double df = (fhigh-flow)/(length-1);
	//double *freq = (double *)malloc(sizeof(double) * length);
	double freq[length];

	cout<<"Freq spacing "<<df<<endl;

	for(int i=0;i<length;i++)
		freq[i]=flow+i*df;
	//############################################################
	//############################################################
	double noise[length];
	populate_noise(freq,"Hanford_O1_fitted", noise,length);
	for (int i =0; i<length;i++){
		noise[i] = noise[i]*noise[i];
	}
	//############################################################

	params.phi = 0;
	params.theta = 0;
	params.tc = tc; 
	celestial_horizon_transform(RA, DEC, gps_time, "Hanford", &params.phi, &params.theta);
	
	fourier_detector_response(freq,length, waveformout,"Hanford",injection_method, &params);
	double temptheta = params.theta;
	for(int j = 0; j<data_length[0]; j++){
		frequencies[0][j] = freq[j];	
		psd[0][j] = noise[j];	
		data[0][j] = waveformout[j];	
	}
	//double snr2 = pow(calculate_snr("Hanford_O1_fitted",waveformout, freq, length),2);
	
	double snr2 = pow(data_snr_maximized_extrinsic(freq,length,waveformout, noise,"Hanford", injection_method,&params),2);

	celestial_horizon_transform(RA, DEC, gps_time, "Livingston", &params.phi, &params.theta);
	params.tc = tc+ DTOA(temptheta, params.theta, "Hanford","Livingston"); 
	std::cout<<"Time at Livingston of injection: "<<params.tc<<std::endl;
	fourier_detector_response(freq,length, waveformout,"Livingston",injection_method, &params);

	for(int j = 0; j<data_length[0]; j++){
		frequencies[1][j] = freq[j];	
		psd[1][j] = (noise[j]);	
		data[1][j] = waveformout[j];	
	}
	//snr2+=pow( calculate_snr("Hanford_O1_fitted",waveformout, freq, length),2);
	snr2 += pow(data_snr_maximized_extrinsic(freq,length,waveformout, noise,"Virgo", injection_method,&params),2);
	//celestial_horizon_transform(RA, DEC, gps_time, "Virgo", &params.phi, &params.theta);
	//params.tc = tc+ DTOA(temptheta, params.theta, "Hanford","Virgo"); 
	//std::cout<<"Time at Virgo of injection: "<<params.tc<<std::endl;
	//fourier_detector_response(freq,length, waveformout,"Virgo",injection_method, &params);

	//for(int j = 0; j<data_length[0]; j++){
	//	frequencies[2][j] = freq[j];	
	//	psd[2][j] = (noise[j]);	
	//	data[2][j] = waveformout[j];	
	//}
	//#########################################################
	//snr2+=pow( calculate_snr("Hanford_O1_fitted",waveformout, freq, length),2);
	snr2 += pow(data_snr_maximized_extrinsic(freq,length,waveformout, noise,"Virgo", injection_method,&params),2);
	std::cout<<"SNR of injection: "<<sqrt(snr2)<<std::endl;
	//snr = data_snr_maximized_extrinsic(freq,length, waveformout,"Hanford_O1_fitted",method,params );
	//std::cout<<"SNR of injection calculated as data: "<<snr<<std::endl;

	delete [] params.bppe;
	delete [] params.betappe;
	
	//#########################################################
	//mcmc options
	//int dimension = 9;
	int dimension = 9;
	//double initial_pos[dimension]={.3, 2., -0.2,log(400),log(40), .24,- .0,-.0};
	double initial_pos[dimension]={.9, 2, .7,log(520),log(10), .2, .5,.0, 5};
	//double initial_pos[dimension]={.9, 2, .7,log(520),log(10), .2, .5,.0};
	//double initial_pos[dimension]={-.9, 2, -1.2,log(410),log(30), .24,-.4,.3};
	//double initial_pos[dimension]={-.0, 0, -0,log(500),log(50), .2,-.0,.0};
	//double initial_pos[dimension]={-.99, 2, -1.2,log(410),log(30.78), .24,-.4,.3};
	double *seeding_var = NULL;
	int n_steps = 100000;
	int chain_N= 10;
	double ***output;
	output = allocate_3D_array( chain_N, n_steps, dimension );
	int swp_freq = 3;
	double chain_temps[chain_N];
	chain_temps[0]=1.;
	double c = 1.2;
	for(int i =1; i < chain_N;  i ++)
		chain_temps[i] = c*chain_temps[i-1];
	
	int Nmod = 1;
	int *bppe = new int[1];
	bppe[0] = -1;
	int numThreads = 4;
	bool pool = true;
	//#########################################################
	//gw options
	//std::string generation_method = "dCS_IMRPhenomD_log";
	//std::string generation_method = "IMRPhenomD";
	std::string generation_method = "dCS_IMRPhenomD_root_alpha";
	
	
	//std::string autocorrfile = "testing/data/auto_corr_mcmc_injection1_1226.csv";
	////std::string autocorrfile = "";
	//std::string chainfile = "testing/data/mcmc_output_injection1_1226.csv";
	//std::string statfilename = "testing/data/mcmc_statistics_injection1_1226.txt";
	//std::string checkfile = "testing/data/mcmc_checkpoint_injection1_1226.csv";
	//std::string autocorrfile = "testing/data/junk.csv";
	//std::string chainfile = "testing/data/mcmc_output_injection1_1226.csv";
	//std::string statfilename = "testing/data/mcmc_statistics_injection1_1226.txt";
	//std::string checkfile = "testing/data/mcmc_checkpoint_injection1_1226.csv";
	std::string autocorrfile = "testing/data/auto_corr_mcmc_injection1_1226.csv";
	std::string chainfile = "testing/data/mcmc_output_injection1_1226.csv";
	std::string statfilename = "testing/data/mcmc_statistics_injection1_1226.txt";
	std::string checkfile = "testing/data/mcmc_checkpoint_injection1_1226.csv";

	PTMCMC_MH_GW(output, dimension, n_steps, chain_N, initial_pos,seeding_var,chain_temps, 
			swp_freq, test_lp_GW_dCS_root_alpha,numThreads, pool,show_progress,
			num_detectors, 
			data, psd,frequencies, data_length,gps_time, detectors,Nmod, bppe,
			generation_method,statfilename,"","", "",checkfile);	
			//generation_method,"","","", "");	
	
	std::cout<<"ended"<<std::endl;
	double **output_transform=(double **)malloc(sizeof(double*)*n_steps);
	for (int j =0; j<n_steps; j++)
		output_transform[j] = (double *)malloc(sizeof(double)*dimension);

	for(int j = 0; j<n_steps;j++){
			output_transform[j][0]=output[0][j][0];
			output_transform[j][1]=output[0][j][1];
			output_transform[j][2]=output[0][j][2];
			output_transform[j][3]=std::exp(output[0][j][3]);
			output_transform[j][4]=std::exp(output[0][j][4]);
			output_transform[j][5]=output[0][j][5];
			output_transform[j][6]=output[0][j][6];
			output_transform[j][7]=output[0][j][7];
			//output_transform[j][8]=exp(output[0][j][8]);
			output_transform[j][8]=(output[0][j][8]);
	}
	write_file(chainfile, output_transform, n_steps, dimension);

	int segs = 10;
	double target_corr = .01;
	//write_file_auto_corr_from_data_file_accel(autocorrfile, chainfile,dimension,n_steps,segs,target_corr);
	//write_auto_corr_file_from_data_file(autocorrfile, chainfile,n_steps,dimension,segs,target_corr,numThreads);

	//output hottest chain too
	chainfile = "testing/data/mcmc_output_injection1_1226_hot.csv";
	for(int j = 0; j<n_steps;j++){
			output_transform[j][0]=output[chain_N-1][j][0];
			output_transform[j][1]=output[chain_N-1][j][1];
			output_transform[j][2]=output[chain_N-1][j][2];
			output_transform[j][3]=std::exp(output[chain_N-1][j][3]);
			output_transform[j][4]=std::exp(output[chain_N-1][j][4]);
			output_transform[j][5]=output[chain_N-1][j][5];
			output_transform[j][6]=output[chain_N-1][j][6];
			output_transform[j][7]=output[chain_N-1][j][7];
			//output_transform[j][8]=exp(output[chain_N-1][j][8]);
			output_transform[j][8]=(output[chain_N-1][j][8]);
	}
	//write_file(chainfile, output_transform, n_steps, dimension);

	deallocate_3D_array(output, chain_N, n_steps, dimension);
	for(int i =0; i< num_detectors; i++){
		free(data[i]);
		free(psd[i]);
		free(frequencies[i]);
	}
	for(int i =0; i< n_steps; i++){
		free(output_transform[i]);
	}
	free(output_transform);
	free(data);
	delete [] bppe;
	free(psd);
	free(frequencies );
	delete [] detectors;
	free(data_length);
}
void test14()
{
	//int raw_length = 28673;
	//int cutoff = 600;
	//int high_cut = 11000;
	//int length = high_cut-cutoff;
	
	int raw_length=  4096;
	//int raw_length=  8032;
	int cutoff = 50;
	int length = raw_length-cutoff;

	int num_detectors =2;
	//int num_detectors =1;
	std::string *detectors = new std::string[num_detectors];//(std::string*)malloc(sizeof(std::string)*50*num_detectors);
	detectors[0] = "Hanford";
	detectors[1] = "Livingston";

	double **temp_data = allocate_2D_array(raw_length,2);
	double *temp_psd = (double *)malloc(sizeof(double)*raw_length);
	double *temp_freq = (double *)malloc(sizeof(double)*raw_length);
	//std::string filebase = "testing/data/gw170608_";
	std::string filebase = "testing/data/gw151226_";
	//std::string filebase = "testing/data/gw150914_";
	//std::string filebase = "testing/data/gw_150914_";
	//double gps_time = 1126259462;//TESTING -- gw150914
	double gps_time = 1135136350.6;//TESTING -- gw151226
	//double gps_time = 1180922494.5; //TESTING -- gw170608

	std::complex<double> **data= (std::complex<double>**)malloc(
			sizeof(std::complex<double>*)*num_detectors);
	double **psd = (double **)malloc(sizeof(double *)*num_detectors);
	double **frequencies = (double **)malloc(sizeof(double *)*num_detectors);
	int *data_length= (int*)malloc(sizeof(int)*num_detectors);
	data_length[0] =length;
	data_length[1] =length;

	for (int i =0; i<num_detectors; i++){
		data[i] = (std::complex<double> *)malloc(
			sizeof(std::complex<double>)*data_length[i]);
		
		psd[i] = (double *)malloc(sizeof(double)*data_length[i]);
		frequencies[i] = (double *)malloc(sizeof(double)*data_length[i]);
	}
	read_file(filebase+"data_H.csv",temp_data, raw_length,2);
	read_file(filebase+"psd_H.csv",temp_psd);
	read_file(filebase+"freq_H.csv",temp_freq);
	//read_file(filebase+"data_H_new.csv",temp_data, raw_length,2);
	//read_file(filebase+"psd_H_new.csv",temp_psd);
	//read_file(filebase+"freq_H_new.csv",temp_freq);
	//read_file(filebase+"data_L.csv",temp_data, raw_length,2);
	//read_file(filebase+"psd_L.csv",temp_psd);
	//read_file(filebase+"freq_L.csv",temp_freq);
	//read_file(filebase+"data.csv",temp_data, raw_length,2);
	//read_file(filebase+"psd.csv",temp_psd);
	//read_file(filebase+"freq.csv",temp_freq);
	
	double delta_f = temp_freq[length/2]-temp_freq[length/2-1];
	//double delta_f = 1;
	//double delta_f = 1./(temp_freq[length/2]-temp_freq[length/2-1]);
	for(int j = 0; j<data_length[0]; j++){
		frequencies[0][j] = temp_freq[j+cutoff];	
		psd[0][j] = 1./delta_f * (temp_psd[j+cutoff]);	
		data[0][j] = 1./delta_f * std::complex<double>(temp_data[j+cutoff][0],temp_data[j+cutoff][1]);	
	}
	read_file(filebase+"data_L.csv",temp_data, raw_length,2);
	read_file(filebase+"psd_L.csv",temp_psd);
	read_file(filebase+"freq_L.csv",temp_freq);
	//read_file(filebase+"data_L_new.csv",temp_data, raw_length,2);
	//read_file(filebase+"psd_L_new.csv",temp_psd);
	//read_file(filebase+"freq_L_new.csv",temp_freq);
	//read_file(filebase+"data_H.csv",temp_data, raw_length,2);
	//read_file(filebase+"psd_H.csv",temp_psd);
	//read_file(filebase+"freq_H.csv",temp_freq);
	for(int j = 0; j<data_length[1]; j++){
		frequencies[1][j] = temp_freq[j+cutoff];	
		psd[1][j] = 1./delta_f *(temp_psd[j+cutoff]);	
		data[1][j] = 1./delta_f *std::complex<double>(temp_data[j+cutoff][0],temp_data[j+cutoff][1]);	
	}

	deallocate_2D_array(temp_data,raw_length,2);
	free(temp_psd);
	free(temp_freq);

	//######### SNR #####################	
	gen_params params;
	//150914
	//double RA = 1.69;
	//double DEC = -1.21;
	//double chirpm = 30.78;
	//double eta =.24;
	//params.mass1 = 38.56;
	//params.mass2 = 33.18;
	//chirpm = calculate_chirpmass(params.mass1,params.mass2);
	//params.spin1[0] = 0;
	//params.spin1[1] = 0;
	//params.spin1[2] = -.0;
	//params.spin2[0] = 0;
	//params.spin2[1] = 0;
	//params.spin2[2] = .0;
	//params.phic = .0;
	//double tc = 2;
	//params.Luminosity_Distance = 431.;
	////params.Luminosity_Distance = 800.;
	//params.NSflag = false;
	//params.incl_angle = M_PI;
	//params.sky_average=false;
	//params.phi = 0;
	//params.theta = 0;
	//params.tc = tc; 
	
	//151226
	//double RA = 3.34;
	//double DEC = -.5;
	//double chirpm = 30.78;
	//params.mass1 = 15.26;
	//params.mass2 = 8.29;
	//chirpm = calculate_chirpmass(params.mass1,params.mass2);
	//params.spin1[0] = 0;
	//params.spin1[1] = 0;
	//params.spin1[2] = .5;
	//params.spin2[0] = 0;
	//params.spin2[1] = 0;
	//params.spin2[2] = .0;
	//params.phic = .0;
	//double tc = 2;
	//params.Luminosity_Distance = 425.;
	//params.NSflag = false;
	//params.incl_angle = M_PI/3.;
	//params.sky_average=false;


	//170806
	double RA = 2.1;
	double DEC = .5;
	double chirpm = .78;
	params.mass1 = 11.9;
	params.mass2 = 8.09;
	chirpm = calculate_chirpmass(params.mass1,params.mass2);
	params.spin1[0] = 0;
	params.spin1[1] = 0;
	params.spin1[2] = .3;
	params.spin2[0] = 0;
	params.spin2[1] = 0;
	params.spin2[2] = .4;
	params.phic = .0;
	double tc = 2;
	params.Luminosity_Distance = 425.;
	params.NSflag1 = false;
	params.NSflag2 = false;
	params.incl_angle = 0;
	params.sky_average=false;

	params.phi = 0;
	params.theta = 0;
	params.tc = tc; 
	celestial_horizon_transform(RA, DEC, gps_time, detectors[0], &params.phi, &params.theta);
	double snr2 = pow(data_snr_maximized_extrinsic(frequencies[0],length,data[0], psd[0],detectors[0], "IMRPhenomD",&params),2);
	celestial_horizon_transform(RA, DEC, gps_time, detectors[1], &params.phi, &params.theta);
	snr2 += pow(data_snr_maximized_extrinsic(frequencies[1],length,data[1], psd[1],detectors[1], "IMRPhenomD",&params),2);
	std::cout<<"SNR of injection: "<<sqrt(snr2)<<std::endl;

	//#########################################################
	//mcmc options
	int dimension = 8;
	//double initial_pos[dimension]={-.9, 1.,-1.,log(400),log(30), .24,- .0,-.0};
	double initial_pos[dimension]={cos(params.incl_angle)*.9, RA,DEC,
		log(params.Luminosity_Distance),log(chirpm), .24,- .0,-.0};
	//double initial_pos[dimension]={-.0, 0.,-0.,log(400),log(10), .24,- .0,-.0};
	//double initial_pos[dimension]={log(8), .24,- .0,-.0};
	//double initial_pos[dimension]={log(200),log(20), .15, 0,0};
	double *seeding_var = NULL;
	int n_steps = 50000;
	int chain_N= 8;
	double ***output;
	output = allocate_3D_array( chain_N, n_steps, dimension );
	//double *initial_pos_ptr = initial_pos;
	int swp_freq = 10;
	//double chain_temps[chain_N] ={1,2,3,10,12};
	double chain_temps[chain_N];
	//double temp_step = 20./(chain_N);
	chain_temps[0]=1.;
	double c = 1.4;
	for(int i =1; i < chain_N;  i ++)
		chain_temps[i] = c*chain_temps[i-1];
		//chain_temps[i] = (1.+i*temp_step);
	int Nmod = 0;
	int *bppe = NULL;
	
	int numThreads = 8;
	bool pool = true;
	//#########################################################
	//gw options
	std::string generation_method = "IMRPhenomD";
	
	
	std::string autocorrfile = "testing/data/auto_corr_mcmc_DFull_Neil.csv";
	//std::string autocorrfile = "";
	std::string chainfile = "testing/data/mcmc_output_DFull_Neil.csv";
	std::string statfilename = "testing/data/mcmc_statistics_DFull_Neil.txt";
	std::string checkfile = "testing/data/mcmc_checkpoint_DFull_Neil.csv";
	PTMCMC_MH_GW(output, dimension, n_steps, chain_N, initial_pos,seeding_var,chain_temps, 
			swp_freq, test_lp_GW_DFull,numThreads, pool, show_progress, 
			num_detectors, 
			data, psd, 
			frequencies, data_length, gps_time,detectors,Nmod, bppe, generation_method,
			statfilename,"",autocorrfile, "",checkfile);	
	std::cout<<"ended"<<std::endl;

	double **output_transform=(double **)malloc(sizeof(double*)*n_steps);
	for (int j =0; j<n_steps; j++)
		output_transform[j] = (double *)malloc(sizeof(double)*dimension);

	for(int j = 0; j<n_steps;j++){
			output_transform[j][0]=output[0][j][0];
			output_transform[j][1]=output[0][j][1];
			output_transform[j][2]=output[0][j][2];
			output_transform[j][3]=std::exp(output[0][j][3]);
			output_transform[j][4]=std::exp(output[0][j][4]);
			output_transform[j][5]=output[0][j][5];
			output_transform[j][6]=output[0][j][6];
			output_transform[j][7]=output[0][j][7];
	}
	write_file(chainfile, output_transform, n_steps, dimension);
	//output hottest chain too
	chainfile = "testing/data/mcmc_output_DFull_hot_Neil.csv";
	for(int j = 0; j<n_steps;j++){
			output_transform[j][0]=output[chain_N-1][j][0];
			output_transform[j][1]=output[chain_N-1][j][1];
			output_transform[j][2]=output[chain_N-1][j][2];
			output_transform[j][3]=std::exp(output[chain_N-1][j][3]);
			output_transform[j][4]=std::exp(output[chain_N-1][j][4]);
			output_transform[j][5]=output[chain_N-1][j][5];
			output_transform[j][6]=output[chain_N-1][j][6];
			output_transform[j][7]=output[chain_N-1][j][7];
	}
	write_file(chainfile, output_transform, n_steps, dimension);

	deallocate_3D_array(output, chain_N, n_steps, dimension);
	for(int i =0; i< num_detectors; i++){
		free(data[i]);
		free(psd[i]);
		free(frequencies[i]);
	}
	for(int i =0; i< n_steps; i++){
		free(output_transform[i]);
	}
	free(output_transform);
	free(data);
	free(psd);
	free(frequencies );
	delete [] detectors;
	free(data_length);
}
void test13()
{
	double RV, RH, RL;
	RV = radius_at_lat(V_LAT, V_elevation);
	RH = radius_at_lat(H_LAT, H_elevation);
	RL = radius_at_lat(L_LAT, L_elevation);
	double dtoa = DTOA(80,80,"Hanford","Virgo");
	
	cout.precision(15);
	std::cout<<"Distance to center of the earth for Virgo: "<<RV<<std::endl;
	std::cout<<"Distance to center of the earth for Livingston: "<<RL<<std::endl;
	std::cout<<"Distance to center of the earth for Hanford: "<<RH<<std::endl;
	std::cout<<"DTOA: "<<dtoa<<std::endl;
}
void test12()
{
	int length = 4000;
	//int num_detectors =3;
	int num_detectors =2;
	std::string *detectors = new std::string[num_detectors];//(std::string*)malloc(sizeof(std::string)*50*num_detectors);
	detectors[0] = "Hanford";
	detectors[1] = "Livingston";
	//detectors[2] = "Virgo";

	double gps_time = 1135136350.6;//TESTING -- gw151226

	std::complex<double> **data= (std::complex<double>**)malloc(
			sizeof(std::complex<double>*)*num_detectors);
	double **psd = (double **)malloc(sizeof(double *)*num_detectors);
	double **frequencies = (double **)malloc(sizeof(double *)*num_detectors);
	int *data_length= (int*)malloc(sizeof(int)*num_detectors);
	data_length[0] =length;
	data_length[1] =length;
	//data_length[2] =length;

	for (int i =0; i<num_detectors; i++){
		data[i] = (std::complex<double> *)malloc(
			sizeof(std::complex<double>)*data_length[i]);
		
		psd[i] = (double *)malloc(sizeof(double)*data_length[i]);
		frequencies[i] = (double *)malloc(sizeof(double)*data_length[i]);
	}
	//#########################################################
	//Make trial data
	gen_params params;
	//double RA = 5.;
	//double DEC = 1.;
	double RA = 2;
	double DEC = .7;
	double chirpm = 9.71;
	double eta =.2;
	params.mass1 = calculate_mass1(chirpm,eta);
	params.mass2 = calculate_mass2(chirpm,eta);
	complex<double> waveformout[length];
	params.spin1[0] = 0;
	params.spin1[1] = .01;
	params.spin1[2] = .43;
	params.spin2[0] = 0;
	params.spin2[1] = 0;
	params.spin2[2] = .0;
	params.phic = .0;
	double tc = 2;
	params.phiRef=1;
	params.f_ref =20;
	params.Luminosity_Distance = 100.;
	params.NSflag1 = false;
	params.NSflag2 = false;
	params.incl_angle = 0.;//M_PI/3.;
	params.sky_average=false;
	params.bppe = new int[1];
	params.betappe = new double[1];
	params.bppe[0] = -1;
	params.betappe[0] = (pow(0./(3e5),4));
	params.psi = 0;
	params.gmst = gps_to_GMST(gps_time);
	params.RA = RA;
	params.DEC = DEC;
	//params.betappe[0] = log(pow(100./(3e5),4));
	params.Nmod = 1;
	//std::string injection_method = "dCS_IMRPhenomD";
	std::string injection_method = "IMRPhenomD";
	//std::string injection_method = "IMRPhenomPv2";
	//params.f_ref = 30.5011;
	//params.phiRef =58.944425/2.;
	
	//############################################################
	double fhigh =500;
	double flow =10;
	double df = (fhigh-flow)/(length-1);
	//double *freq = (double *)malloc(sizeof(double) * length);
	double freq[length];

	cout<<"Freq spacing "<<df<<endl;

	for(int i=0;i<length;i++)
		freq[i]=flow+i*df;
	//############################################################
	double noise[length];
	populate_noise(freq,"Hanford_O1_fitted", noise,length);
	for (int i =0; i<length;i++){
		noise[i] = noise[i]*noise[i];
	}
	//############################################################

	params.phi = 0;
	params.theta = 0;
	params.tc = tc; 
	celestial_horizon_transform(RA, DEC, gps_time, "Hanford", &params.phi, &params.theta);
	
	fourier_detector_response_equatorial(freq,length, waveformout,"Hanford",injection_method, &params);
	double temptheta = params.theta;
	for(int j = 0; j<data_length[0]; j++){
		frequencies[0][j] = freq[j];	
		psd[0][j] = noise[j];	
		data[0][j] = waveformout[j];	
	}
	//double snr2 = pow(calculate_snr("Hanford_O1_fitted",waveformout, freq, length),2);
	double snr2 = pow(data_snr_maximized_extrinsic(freq,length,waveformout, noise,"Hanford", injection_method,&params),2);

	celestial_horizon_transform(RA, DEC, gps_time, "Livingston", &params.phi, &params.theta);
	params.tc = tc+ DTOA(temptheta, params.theta, "Hanford","Livingston"); 
	std::cout<<"Time at Livingston of injection: "<<params.tc<<std::endl;
	fourier_detector_response_equatorial(freq,length, waveformout,"Livingston",injection_method, &params);

	for(int j = 0; j<data_length[0]; j++){
		frequencies[1][j] = freq[j];	
		psd[1][j] = (noise[j]);	
		data[1][j] = waveformout[j];	
	}
	//snr2+=pow( calculate_snr("Hanford_O1_fitted",waveformout, freq, length),2);
	snr2 += pow(data_snr_maximized_extrinsic(freq,length,waveformout, noise,"Livingston", injection_method,&params),2);
	//celestial_horizon_transform(RA, DEC, gps_time, "Virgo", &params.phi, &params.theta);
	//params.tc = tc+ DTOA(temptheta, params.theta, "Hanford","Virgo"); 
	//std::cout<<"Time at Virgo of injection: "<<params.tc<<std::endl;
	//fourier_detector_response(freq,length, waveformout,"Virgo",injection_method, &params);

	//for(int j = 0; j<data_length[0]; j++){
	//	frequencies[2][j] = freq[j];	
	//	psd[2][j] = (noise[j]);	
	//	data[2][j] = waveformout[j];	
	//}
	////#########################################################
	////snr2+=pow( calculate_snr("Hanford_O1_fitted",waveformout, freq, length),2);
	//snr2 += pow(data_snr_maximized_extrinsic(freq,length,waveformout, noise,"Virgo", injection_method,&params),2);
	std::cout<<"SNR of injection: "<<sqrt(snr2)<<std::endl;
	//snr = data_snr_maximized_extrinsic(freq,length, waveformout,"Hanford_O1_fitted",method,params );
	//std::cout<<"SNR of injection calculated as data: "<<snr<<std::endl;

	delete [] params.bppe;
	delete [] params.betappe;
	//#########################################################
	//mcmc options
	int dimension = 14;
	double initial_pos[dimension]={.99,2,.7,log(100),log(9.7), .2, .43,.01,.1,.1, .1,.1,1, .1};
	double *seeding_var = NULL;
	int n_steps = 150000;
	int chain_N= 12;
	double ***output;
	output = allocate_3D_array( chain_N, n_steps, dimension );
	int swp_freq = 5;
	double chain_temps[chain_N];
	chain_temps[0]=1.;
	double c = 1.1;
	for(int i =1; i < chain_N;  i ++)
		chain_temps[i] = c*chain_temps[i-1];
	
	int Nmod = 0;
	int *bppe = NULL;
	int numThreads = 10;
	bool pool = true;
	//#########################################################
	//gw options
	std::string generation_method = "IMRPhenomPv2";
	
	
	std::string autocorrfile = "testing/data/auto_corr_mcmc_Pv2.csv";
	//std::string autocorrfile = "";
	std::string chainfile = "testing/data/mcmc_output_Pv2.csv";
	std::string statfilename = "testing/data/mcmc_statistics_Pv2.txt";
	std::string checkpointfile = "testing/data/mcmc_checkpoint_Pv2.csv";

	PTMCMC_MH_GW(output, dimension, n_steps, chain_N, initial_pos,seeding_var,chain_temps, 
			swp_freq, test_lp_GW_Pv2,numThreads, pool,show_progress, num_detectors, data, psd, 
			frequencies, data_length, gps_time,detectors,Nmod, bppe, generation_method,
			statfilename,"",autocorrfile, "",checkpointfile);	
	std::cout<<"ended"<<std::endl;

	double **output_transform=(double **)malloc(sizeof(double*)*n_steps);
	for (int j =0; j<n_steps; j++)
		output_transform[j] = (double *)malloc(sizeof(double)*dimension);

	for(int j = 0; j<n_steps;j++){
			//output_transform[j][0]=std::exp(output[0][j][0])/mpc_sec;
			output_transform[j][0] = output[0][j][0];
			output_transform[j][1]=std::exp(output[0][j][1]);
			output_transform[j][2]=output[0][j][2];
			output_transform[j][3]=output[0][j][3];
			output_transform[j][4]=output[0][j][4];
			output_transform[j][5]=output[0][j][5];
			output_transform[j][6]=output[0][j][6];
	}
	write_file(chainfile, output[0], n_steps, dimension);
	//output hottest chain too
	chainfile = "testing/data/mcmc_output_Pv2_hot.csv";
	for(int j = 0; j<n_steps;j++){
			output_transform[j][0] = output[chain_N-1][j][0];
			output_transform[j][1]=std::exp(output[chain_N-1][j][1]);
			output_transform[j][2]=output[chain_N-1][j][2];
			output_transform[j][3]=output[chain_N-1][j][3];
			output_transform[j][4]=output[chain_N-1][j][4];
			output_transform[j][5]=output[chain_N-1][j][5];
			output_transform[j][6]=output[chain_N-1][j][6];
	}
	write_file(chainfile, output[chain_N-1], n_steps, dimension);

	deallocate_3D_array(output, chain_N, n_steps, dimension);
	for(int i =0; i< num_detectors; i++){
		free(data[i]);
		free(psd[i]);
		free(frequencies[i]);
	}
	for(int i =0; i< n_steps; i++){
		free(output_transform[i]);
	}
	free(output_transform);
	free(data);
	free(psd);
	free(frequencies );
	delete [] detectors;
	free(data_length);
}
void test11()
{
	int length = 4000;
	int loops = 1000;
	
	//std::complex<double> **out = (std::complex<double> **)malloc(sizeof(std::complex<double>)*loops);

	//for (int i =0 ; i< loops; i++)
	//{
	//	out[i] = (std::complex<double>*)malloc(sizeof(std::complex<double>)*length);
	//}
	double fhigh =200;
	double flow =10;
	double df = (fhigh-flow)/(length-1);
	freq = (double *)malloc(sizeof(double) * length);
	for(int i=0;i<length;i++)
		freq[i]=flow+i*df;

	//Synthetic data
	double gps_time = 1135136350.6;//TESTING -- gw151226
	gen_params params_data;
	double chirpm = 50.78;
	double eta =.24;
	double Ra = 1;
	double Dec = -1;
	params_data.mass1 = calculate_mass1(chirpm,eta);
	params_data.mass2 = calculate_mass2(chirpm,eta);
	string method= "IMRPhenomD";
	//complex<double> waveformout[length];
	params_data.spin1[0] = 0;
	params_data.spin1[1] = 0;
	params_data.spin1[2] = -.4;
	params_data.spin2[0] = 0;
	params_data.spin2[1] = 0;
	params_data.spin2[2] = .3;
	params_data.phic = .0;
	params_data.tc = 5;
	params_data.Luminosity_Distance = 400.;
	params_data.NSflag1 = false;
	params_data.NSflag2 = false;
	celestial_horizon_transform(Ra, Dec, gps_time, "Hanford", &params_data.phi, &params_data.theta);	
	//params_data.phi = 0;
	//params_data.theta = 0;
	params_data.incl_angle = M_PI/3.;
	params_data.sky_average=false;
	waveformout11 = (std::complex<double>*) malloc(sizeof(std::complex<double>)*length);
	std::complex<double> *waveformout11L = (std::complex<double>*) malloc(sizeof(std::complex<double>)*length);
	
	//fourier_detector_response(freq, length, waveformout11,method,&params_data);
	fourier_detector_response(freq, length, waveformout11,"Hanford",method,&params_data);
	double temptheta=params_data.theta;
	celestial_horizon_transform(Ra, Dec, gps_time, "Livingston", &params_data.phi, &params_data.theta);	

	params_data.tc = params_data.tc + DTOA(temptheta, params_data.theta,"Hanford", "Livingston");
	fourier_detector_response(freq, length, waveformout11L,"Livingston",method,&params_data);

	//double psd[length];
	psd = (double *) malloc(sizeof(double)*length);
	populate_noise(freq,"Hanford_O1_fitted", psd,length);
	for (int i =0; i<length;i++){
		psd[i] = psd[i]*psd[i];
	}
	//fftw_outline *plans = (fftw_outline *)malloc(sizeof(fftw_outline)*loops);	
	//for(int i =0 ; i<loops; i ++)
	//{
	//	initiate_likelihood_function(&plans[i] , length);
	//}
	fftw_outline plan;
	allocate_FFTW_mem_forward(&plan, length);
	double ll_true =maximized_Log_Likelihood(waveformout11, psd, freq,
			length, &params_data, "Livingston", "IMRPhenomD", &plan);
	std::cout<<"TRUE ll: "<<ll_true<<std::endl;	
	omp_set_num_threads(10);

	double *DLs = (double *)malloc(sizeof(double)*loops);
	double *lls1 = (double *)malloc(sizeof(double)*loops);
	double *lls2 = (double *)malloc(sizeof(double)*loops);
	double *lls3 = (double *)malloc(sizeof(double)*loops);
	double DLmin = 50;
	double DLmax = 1000;
	double DLstep = (DLmax - DLmin)/loops;
	for (int i =0; i<loops; i++)
	{
		DLs[i] = DLstep * i + DLmin;
	}
	//#pragma omp parallel 
	{
		//#pragma omp for
		for (int i =0; i<loops; i ++)
		{
			gen_params params;
			double chirpm = 50.78;
			double eta =.24;
			params.mass1 = calculate_mass1(chirpm,eta);
			params.mass2 = calculate_mass2(chirpm,eta);
			string method= "IMRPhenomD";
			params.spin1[0] = 0;
			params.spin1[1] = 0;
			params.spin1[2] = -.4;
			params.spin2[0] = 0;
			params.spin2[1] = 0;
			params.spin2[2] = .3;
			params.phic = .0;
			params.tc = 0;
			params.Luminosity_Distance = DLs[i];
			params.NSflag1 = false;
			params.NSflag2 = false;
			celestial_horizon_transform(Ra, Dec, gps_time, "Hanford",&params.phi,&params.theta);	
			params.incl_angle = 0;
			params.sky_average=false;

			double tc,phic;
			lls1[i] =maximized_coal_Log_Likelihood(waveformout11, psd, freq,
					length, &params, "Hanford", "IMRPhenomD", &plan, &tc, &phic);
			//lls1[i]=0;
			//tc=5;
			//phic = 0;

			//std::cout<<"LL 1: "<<lls1[i]<<std::endl;
			double phi1 = params.phi;
			double theta1 = params.theta;
			//std::cout<<"HANFORD ANGLES: "<<params.theta<<" "<<params.phi<<std::endl;
			celestial_horizon_transform(Ra, Dec, gps_time, "Livingston",&params.phi,&params.theta);	
			
			//std::cout<<"LIV ANGLES: "<<params.theta<<" "<<params.phi<<std::endl;
			//std::cout<<"HANFORD TIME: "<<tc<<std::endl;
			params.tc = tc + DTOA(theta1,params.theta, "Hanford", "Livingston");
			//std::cout<<"LIV TIME: "<<params.tc<<std::endl;
			//std::cout<<"LIV TIME REAL: "<<params_data.tc<<std::endl;
			params.phic = phic;
			//params.tc = params_data.tc;
			lls1[i] +=Log_Likelihood(waveformout11L, psd, freq,
					length, &params, "Livingston", "IMRPhenomD", &plan);
			//params.tc = 0;
			//params.phic =0;
			//lls1[i] +=maximized_coal_Log_Likelihood(waveformout11L, psd, freq,
			//		length, &params, "Livingston", "IMRPhenomD", &plan, &tc, &phic);
			//std::cout<<"LL 2: "<<lls1[i]<<std::endl;

			params.incl_angle = M_PI/3.;
			params.tc = 0;
			params.phic = 0;
			celestial_horizon_transform(Ra, Dec, gps_time, "Hanford",&params.phi,&params.theta);	
			lls2[i] =maximized_coal_Log_Likelihood(waveformout11, psd, freq,
					length, &params, "Hanford", "IMRPhenomD", &plan, &tc, &phic);
			//std::cout<<"HANFORD TIME: "<<tc<<std::endl;
			phi1 = params.phi;
			theta1 = params.theta;
			celestial_horizon_transform(Ra, Dec, gps_time, "Livingston",&params.phi,&params.theta);	
			
			params.tc = tc + DTOA(theta1,params.theta, "Hanford", "Livingston");
			//std::cout<<"LIV TIME: "<<params.tc<<std::endl;
			params.phic = phic;
			lls2[i] +=Log_Likelihood(waveformout11L, psd, freq,
					length, &params, "Livingston", "IMRPhenomD", &plan);

			params.incl_angle = -M_PI/4.;
			params.tc = 0;
			params.phic = 0;
			celestial_horizon_transform(Ra, Dec, gps_time, "Hanford",&params.phi,&params.theta);	
			lls3[i] =maximized_coal_Log_Likelihood(waveformout11, psd, freq,
					length, &params, "Hanford", "IMRPhenomD", &plan, &tc, &phic);
			phi1 = params.phi;
			theta1 = params.theta;
			celestial_horizon_transform(Ra, Dec, gps_time, "Livingston",&params.phi,&params.theta);	
			
			params.tc =tc+ DTOA(theta1,params.theta, "Hanford", "Livingston");
			params.phic = phic;
			lls3[i] +=Log_Likelihood(waveformout11L, psd, freq,
					length, &params, "Livingston", "IMRPhenomD", &plan);
			//if(ll != ll_true)
			//{
			//	std::cout<<ll<<std::endl;
			//}
		}
	}
	double **output = allocate_2D_array(loops, 4);
	for (int i=0;i<loops; i++){
		output[i][0] = DLs[i];
		output[i][1] = lls1[i];
		output[i][2] = lls2[i];
		output[i][3] = lls3[i];
	}
	std::string filename = "testing/data/dl_ll.csv";
	write_file(filename, output, loops, 4);
	//for(int i =1;i < loops; i++)
	//{
	//	for(int j =0; j <length; j++){
	//		if(out[i-1][j]!=out[i][j]){
	//			std::cout<<"YIKES"<<std::endl;
	//		}
	//	}
	//}

	//for(int i =0;i < loops; i++)
	//{
	//	free(out[i]);
	//}
	//free(out);
	//for(int i =0 ; i<loops; i ++)
	//{
	//	deactivate_likelihood_function(&plans[i]);
	//}
	//free(plans);
	deallocate_FFTW_mem(&plan);
	deallocate_2D_array(output, loops,2);
	free(freq);
	free(DLs);
	free(lls1);
	free(lls2);
	free(lls3);
	free(psd);
	free(waveformout11);
	free(waveformout11L);
}
void test10()
{
	//int raw_length = 28673;
	//int cutoff = 600;
	//int high_cut = 11000;
	//int length = high_cut-cutoff;
	int cutoff = 40;
	int raw_length=  4096;
	int length = raw_length-cutoff;
	//int length = raw_length-cutoff;
	int num_detectors =2;

	double **temp_data = allocate_2D_array(raw_length,2);
	double *temp_psd = (double *)malloc(sizeof(double)*raw_length);
	double *temp_freq = (double *)malloc(sizeof(double)*raw_length);
	std::string filebase = "testing/data/gw150914_";
	double gps_time = 1126259462;//TESTING -- gw150914
	//double gps_time = 1135136350.6;//TESTING -- gw151226
	//double gps_time = 1180922494.5; //TESTING -- gw170608

	std::complex<double> **data= (std::complex<double>**)malloc(
			sizeof(std::complex<double>*)*num_detectors);
	double **psd = (double **)malloc(sizeof(double *)*num_detectors);
	double **frequencies = (double **)malloc(sizeof(double *)*num_detectors);
	int *data_length= (int*)malloc(sizeof(int)*num_detectors);
	data_length[0] =length;
	data_length[1] =length;

	//bool check = true;
	for (int i =0; i<num_detectors; i++){
		data[i] = (std::complex<double> *)malloc(
			sizeof(std::complex<double>)*data_length[i]);
		
		psd[i] = (double *)malloc(sizeof(double)*data_length[i]);
		frequencies[i] = (double *)malloc(sizeof(double)*data_length[i]);
	}
	read_file(filebase+"data_H.csv",temp_data, raw_length,2);
	read_file(filebase+"psd_H.csv",temp_psd);
	read_file(filebase+"freq_H.csv",temp_freq);
	for(int j = 0; j<data_length[0]; j++){
		frequencies[0][j] = temp_freq[j+cutoff];	
		psd[0][j] = (temp_psd[j+cutoff]);	
		data[0][j] = std::complex<double>(temp_data[j+cutoff][0],temp_data[j+cutoff][1]);	
	}
	read_file(filebase+"data_L.csv",temp_data, raw_length,2);
	read_file(filebase+"psd_L.csv",temp_psd);
	read_file(filebase+"freq_L.csv",temp_freq);
	for(int j = 0; j<data_length[1]; j++){
		frequencies[1][j] = temp_freq[j+cutoff];	
		psd[1][j] = (temp_psd[j+cutoff]);	
		data[1][j] = std::complex<double>(temp_data[j+cutoff][0],temp_data[j+cutoff][1]);	
	}

	deallocate_2D_array(temp_data,raw_length,2);
	free(temp_psd);
	free(temp_freq);
	//#########################################################
	//MCMC options
	int dimension = 7;
	double initial_pos[dimension]={log(400),2,2,log(30), .24, 0,0};
	//double initial_pos[dimension]={log(200*MPC_SEC),log(20*MSOL_SEC), .15, 0,0};
	double *seeding_var = NULL;
	int N_steps = 50000;
	int chain_N= 8;
	double ***output;
	output = allocate_3D_array( chain_N, N_steps, dimension );
	//double *initial_pos_ptr = initial_pos;
	int swp_freq = 5;
	//double chain_temps[chain_N] ={1,2,3,10,12};
	double chain_temps[chain_N];
	//double temp_step = 20./(chain_N);
	chain_temps[0]=1.;
	double c = 1.2;
	for(int i =1; i < chain_N;  i ++)
		chain_temps[i] = c*chain_temps[i-1];
		//chain_temps[i] = (1.+i*temp_step);
	int Nmod = 0;
	int *bppe = NULL;
	
	int numThreads = 20;
	bool pool = false;
	//#########################################################
	//GW options
	std::string *detectors = new std::string[num_detectors];//(std::string*)malloc(sizeof(std::string)*50*num_detectors);
	detectors[0] = "Hanford";
	detectors[1] = "Livingston";
	std::string generation_method = "IMRPhenomD";
	
	
	std::string autocorrfile = "testing/data/auto_corr_mcmc_7dim.csv";
	std::string chainfile = "testing/data/mcmc_output_7dim.csv";
	std::string statfilename = "testing/data/mcmc_statistics_7dim.txt";
	std::string checkfile= "testing/data/mcmc_checkpoint_7dim.csv";

	PTMCMC_MH_GW(output, dimension, N_steps, chain_N, initial_pos,seeding_var,chain_temps, 
			swp_freq, test_lp_GW_7dim,numThreads, pool,show_progress, num_detectors, data, psd, 
			frequencies, data_length,gps_time, detectors,Nmod, bppe, generation_method,
			statfilename,"",autocorrfile, "",checkfile);	
	std::cout<<"ENDED"<<std::endl;

	double **output_transform=(double **)malloc(sizeof(double*)*N_steps);
	for (int j =0; j<N_steps; j++)
		output_transform[j] = (double *)malloc(sizeof(double)*dimension);

	for(int j = 0; j<N_steps;j++){
			output_transform[j][0]=std::exp(output[0][j][0]);
			output_transform[j][1]=output[0][j][1];
			output_transform[j][2]=output[0][j][2];
			output_transform[j][3]=std::exp(output[0][j][3]);
			output_transform[j][4]=output[0][j][4];
			output_transform[j][5]=output[0][j][5];
			output_transform[j][6]=output[0][j][6];
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
void test9()
{
	//int raw_length = 28673;
	//int cutoff = 600;
	//int high_cut = 11000;
	//int length = high_cut-cutoff;
	
	int raw_length=  4096;
	//int raw_length=  8032;
	int cutoff = 50;
	int length = raw_length-cutoff;

	int num_detectors =2;
	//int num_detectors =1;
	std::string *detectors = new std::string[num_detectors];//(std::string*)malloc(sizeof(std::string)*50*num_detectors);
	detectors[0] = "Hanford";
	detectors[1] = "Livingston";

	double **temp_data = allocate_2D_array(raw_length,2);
	double *temp_psd = (double *)malloc(sizeof(double)*raw_length);
	double *temp_freq = (double *)malloc(sizeof(double)*raw_length);
	//std::string filebase = "testing/data/gw170608_";
	//std::string filebase = "testing/data/gw151226_";
	std::string filebase = "testing/data/gw150914_";
	//std::string filebase = "testing/data/gw_150914_";
	double gps_time = 1126259462;//TESTING -- gw150914
	//double gps_time = 1135136350.6;//TESTING -- gw151226
	//double gps_time = 1180922494.5; //TESTING -- gw170608

	std::complex<double> **data= (std::complex<double>**)malloc(
			sizeof(std::complex<double>*)*num_detectors);
	double **psd = (double **)malloc(sizeof(double *)*num_detectors);
	double **frequencies = (double **)malloc(sizeof(double *)*num_detectors);
	int *data_length= (int*)malloc(sizeof(int)*num_detectors);
	data_length[0] =length;
	data_length[1] =length;

	for (int i =0; i<num_detectors; i++){
		data[i] = (std::complex<double> *)malloc(
			sizeof(std::complex<double>)*data_length[i]);
		
		psd[i] = (double *)malloc(sizeof(double)*data_length[i]);
		frequencies[i] = (double *)malloc(sizeof(double)*data_length[i]);
	}
	read_file(filebase+"data_H.csv",temp_data, raw_length,2);
	read_file(filebase+"psd_H.csv",temp_psd);
	read_file(filebase+"freq_H.csv",temp_freq);
	//read_file(filebase+"data_H_new.csv",temp_data, raw_length,2);
	//read_file(filebase+"psd_H_new.csv",temp_psd);
	//read_file(filebase+"freq_H_new.csv",temp_freq);
	//read_file(filebase+"data.csv",temp_data, raw_length,2);
	//read_file(filebase+"psd.csv",temp_psd);
	//read_file(filebase+"freq.csv",temp_freq);
	//double delta_f = temp_freq[length/2]-temp_freq[length/2-1];
	double delta_f = 1;
	for(int j = 0; j<data_length[0]; j++){
		frequencies[0][j] = temp_freq[j+cutoff];	
		psd[0][j] =  (temp_psd[j+cutoff]);	
		data[0][j] = 1./delta_f * std::complex<double>(temp_data[j+cutoff][0],temp_data[j+cutoff][1]);	
	}
	read_file(filebase+"data_L.csv",temp_data, raw_length,2);
	read_file(filebase+"psd_L.csv",temp_psd);
	read_file(filebase+"freq_L.csv",temp_freq);
	//read_file(filebase+"data_L_new.csv",temp_data, raw_length,2);
	//read_file(filebase+"psd_L_new.csv",temp_psd);
	//read_file(filebase+"freq_L_new.csv",temp_freq);
	for(int j = 0; j<data_length[1]; j++){
		frequencies[1][j] = temp_freq[j+cutoff];	
		psd[1][j] =  (temp_psd[j+cutoff]);	
		data[1][j] = 1./delta_f * std::complex<double>(temp_data[j+cutoff][0],temp_data[j+cutoff][1]);	
	}

	deallocate_2D_array(temp_data,raw_length,2);
	free(temp_psd);
	free(temp_freq);
	
	//gen_params params;
	//imrphenomd<double> modeld;
	//params.mass1 = 63.187;
	//params.mass2 = 17.331;
	//string method= "imrphenomd";
	//params.spin1[0] = 0;
	//params.spin1[1] = 0;
	//params.spin1[2] = -.16;
	//params.spin2[0] = 0;
	//params.spin2[1] = 0;
	//params.spin2[2] = .01388;
	//params.phic = .0;
	//params.tc = 0;
	//params.luminosity_distance = 100.;
	//params.nsflag = false;
	//params.phi = 0;
	//params.theta = 0;
	//params.incl_angle = 0;
	//params.sky_average=false;
	//fftw_outline plan;
	//initiate_likelihood_function(&plan, length);
	//double ll =maximized_log_likelihood(data[0], psd[0], frequencies[0],
	//		length, &params, "hanford", "imrphenomd", &plan);
	//double snr = data_snr_maximized_extrinsic(frequencies[0],length, data[0],"Hanford_O1_fitted","imrphenomd",params );
	//std::cout<<"ll "<<ll<<std::endl;
	//std::cout<<"snr "<<snr<<std::endl;
	//deactivate_likelihood_function(&plan);
	//exit(1);
	//#########################################################
	//mcmc options
	int dimension = 4;
	//double initial_pos[dimension]={log(30), .24,- .0,-.0};
	double initial_pos[dimension]={log(8), .12,- .0,-.0};
	//double initial_pos[dimension]={log(200),log(20), .15, 0,0};
	double *seeding_var = NULL;
	int n_steps = 40000;
	int chain_N= 8;
	double ***output;
	output = allocate_3D_array( chain_N, n_steps, dimension );
	//double *initial_pos_ptr = initial_pos;
	int swp_freq = 5;
	//double chain_temps[chain_N] ={1,2,3,10,12};
	double chain_temps[chain_N];
	//double temp_step = 20./(chain_N);
	chain_temps[0]=1.;
	double c = 1.2;
	for(int i =1; i < chain_N;  i ++)
		chain_temps[i] = c*chain_temps[i-1];
		//chain_temps[i] = (1.+i*temp_step);
	int Nmod = 0;
	int *bppe = NULL;
	
	int numThreads = 20;
	bool pool = true;
	//#########################################################
	//gw options
	std::string generation_method = "IMRPhenomD";
	
	
	std::string autocorrfile = "testing/data/auto_corr_mcmc.csv";
	//std::string autocorrfile = "";
	std::string chainfile = "testing/data/mcmc_output.csv";
	std::string statfilename = "testing/data/mcmc_statistics.txt";
	std::string checkfile = "testing/data/mcmc_checkpoint.csv";

	PTMCMC_MH_GW(output, dimension, n_steps, chain_N, initial_pos,seeding_var,chain_temps, 
			swp_freq, test_lp_GW,numThreads, pool,show_progress, num_detectors, data, psd, 
			frequencies, data_length,gps_time, detectors,Nmod, bppe, generation_method,
			statfilename,"",autocorrfile, "",checkfile);	
	std::cout<<"ended"<<std::endl;

	double **output_transform=(double **)malloc(sizeof(double*)*n_steps);
	for (int j =0; j<n_steps; j++)
		output_transform[j] = (double *)malloc(sizeof(double)*dimension);

	for(int j = 0; j<n_steps;j++){
			//output_transform[j][0]=std::exp(output[0][j][0])/mpc_sec;
			output_transform[j][0]=std::exp(output[0][j][0]);
			output_transform[j][1]=output[0][j][1];
			output_transform[j][2]=output[0][j][2];
			output_transform[j][3]=output[0][j][3];
	}
	write_file(chainfile, output_transform, n_steps, dimension);
	//output hottest chain too
	chainfile = "testing/data/mcmc_output_hot.csv";
	for(int j = 0; j<n_steps;j++){
			//output_transform[j][0]=std::exp(output[0][j][0])/mpc_sec;
			output_transform[j][0]=std::exp(output[chain_N-1][j][0]);
			output_transform[j][1]=output[chain_N-1][j][1];
			output_transform[j][2]=output[chain_N-1][j][2];
			output_transform[j][3]=output[chain_N-1][j][3];
	}
	write_file(chainfile, output_transform, n_steps, dimension);

	deallocate_3D_array(output, chain_N, n_steps, dimension);
	for(int i =0; i< num_detectors; i++){
		free(data[i]);
		free(psd[i]);
		free(frequencies[i]);
	}
	for(int i =0; i< n_steps; i++){
		free(output_transform[i]);
	}
	free(output_transform);
	free(data);
	free(psd);
	free(frequencies );
	delete [] detectors;
	free(data_length);
}
void test8()
{
	int length = 5000;
	int num_detectors =3;
	std::string *detectors = new std::string[num_detectors];//(std::string*)malloc(sizeof(std::string)*50*num_detectors);
	detectors[0] = "Hanford";
	detectors[1] = "Livingston";
	detectors[2] = "Virgo";

	double gps_time = 1135136350.6;//TESTING -- gw151226

	std::complex<double> **data= (std::complex<double>**)malloc(
			sizeof(std::complex<double>*)*num_detectors);
	double **psd = (double **)malloc(sizeof(double *)*num_detectors);
	double **frequencies = (double **)malloc(sizeof(double *)*num_detectors);
	int *data_length= (int*)malloc(sizeof(int)*num_detectors);
	data_length[0] =length;
	data_length[1] =length;
	data_length[2] =length;

	for (int i =0; i<num_detectors; i++){
		data[i] = (std::complex<double> *)malloc(
			sizeof(std::complex<double>)*data_length[i]);
		
		psd[i] = (double *)malloc(sizeof(double)*data_length[i]);
		frequencies[i] = (double *)malloc(sizeof(double)*data_length[i]);
	}
	//#########################################################
	//Make trial data
	gen_params params;
	std::string generation_method = "IMRPhenomD";
	double RA = 2.;
	double DEC = -1.2;
	double chirpm = 20.78;
	double eta =.249;
	params.mass1 = calculate_mass1(chirpm,eta);
	params.mass2 = calculate_mass2(chirpm,eta);
	complex<double> waveformout[length];
	params.spin1[0] = 0;
	params.spin1[1] = 0;
	params.spin1[2] = -.4;
	params.spin2[0] = 0;
	params.spin2[1] = 0;
	params.spin2[2] = .3;
	params.phic = .0;
	double tc = 1;
	params.Luminosity_Distance = 810.;
	params.NSflag1 = false;
	params.NSflag2 = false;
	params.incl_angle = 0;
	params.sky_average=false;
	//params.f_ref = 30.5011;
	//params.phiRef =58.944425/2.;
	
	//############################################################
	double fhigh =400;
	double flow =10;
	double df = (fhigh-flow)/(length-1);
	//double *freq = (double *)malloc(sizeof(double) * length);
	double freq[length];

	cout<<"Freq spacing "<<df<<endl;

	for(int i=0;i<length;i++)
		freq[i]=flow+i*df;
	//############################################################
	//############################################################
	double noise[length];
	populate_noise(freq,"Hanford_O1_fitted", noise,length);
	for (int i =0; i<length;i++){
		noise[i] = noise[i]*noise[i];
	}
	//############################################################

	params.phi = 0;
	params.theta = 0;
	params.tc = tc; 
	celestial_horizon_transform(RA, DEC, gps_time, "Hanford", &params.phi, &params.theta);
	fourier_detector_response(freq,length, waveformout,"Hanford","IMRPhenomD", &params);
	double temptheta = params.theta;
	for(int j = 0; j<data_length[0]; j++){
		frequencies[0][j] = freq[j];	
		psd[0][j] = noise[j];	
		data[0][j] = waveformout[j];	
	}
	//double snr2 = pow(calculate_snr("Hanford_O1_fitted",waveformout, freq, length),2);

	celestial_horizon_transform(RA, DEC, gps_time, "Livingston", &params.phi, &params.theta);
	params.tc = tc+ DTOA(temptheta, params.theta, "Hanford","Livingston"); 
	fourier_detector_response(freq,length, waveformout,"Livingston","IMRPhenomD", &params);
	for(int j = 0; j<data_length[0]; j++){
		frequencies[1][j] = freq[j];	
		psd[1][j] = (noise[j]);	
		data[1][j] = waveformout[j];	
	}
	//snr2+=pow( calculate_snr("Hanford_O1_fitted",waveformout, freq, length),2);
	celestial_horizon_transform(RA, DEC, gps_time, "Virgo", &params.phi, &params.theta);
	params.tc = tc+ DTOA(temptheta, params.theta, "Hanford","Virgo"); 
	fourier_detector_response(freq,length, waveformout,"Virgo","IMRPhenomD", &params);
	for(int j = 0; j<data_length[0]; j++){
		frequencies[2][j] = freq[j];	
		psd[2][j] = (noise[j]);	
		data[2][j] = waveformout[j];	
	}
	//#########################################################
	//snr2+=pow( calculate_snr("Hanford_O1_fitted",waveformout, freq, length),2);
	//std::cout<<"SNR of injection: "<<sqrt(snr2)<<std::endl;
	
	//#########################################################
	//MCMC options
	int dimension = 4;
	//double initial_pos[dimension]={log(chirpm*MSOL_SEC), eta, params.spin1[2],params.spin2[2]};
	double initial_pos[dimension]={log(10), .2, 0,0};
	//double initial_pos[dimension]={log(200),log(20), .15, 0,0};
	double *seeding_var = NULL;
	int N_steps = 15000;
	int chain_N= 8;
	double ***output;
	output = allocate_3D_array( chain_N, N_steps, dimension );
	//double *initial_pos_ptr = initial_pos;
	int swp_freq = 10;
	//double chain_temps[chain_N] ={1,2,3,10,12};
	double chain_temps[chain_N];
	//double temp_step = 20./(chain_N);
	chain_temps[0]=1.;
	double c = 1.2;
	for(int i =1; i < chain_N;  i ++)
		chain_temps[i] = c*chain_temps[i-1];
		//chain_temps[i] = (1.+i*temp_step);
	
	int Nmod = 0;
	int *bppe = NULL;
	
	std::string autocorrfile = "testing/data/auto_corr_mcmc.csv";
	std::string chainfile = "testing/data/mcmc_output.csv";
	std::string statfilename = "testing/data/mcmc_statistics.txt";
	std::string checkfile= "testing/data/mcmc_checkpoint.csv";

	
	int numThreads = 20;
	bool pool = true;

	PTMCMC_MH_GW(output, dimension, N_steps, chain_N, initial_pos,seeding_var,chain_temps, 
			swp_freq, test_lp_GW,numThreads, pool,show_progress, 
			num_detectors, 
			data, psd, 
			frequencies, data_length, gps_time,detectors,Nmod, bppe, generation_method,
			statfilename, "" ,autocorrfile, "",checkfile);	
	std::cout<<"ENDED"<<std::endl;

	double **output_transform=(double **)malloc(sizeof(double*)*N_steps);
	for (int j =0; j<N_steps; j++)
		output_transform[j] = (double *)malloc(sizeof(double)*dimension);

	for(int j = 0; j<N_steps;j++){
			output_transform[j][0]=std::exp(output[0][j][0]);
			output_transform[j][1]=output[0][j][1];
			output_transform[j][2]=output[0][j][2];
			output_transform[j][3]=output[0][j][3];
	}
	write_file(chainfile, output_transform, N_steps, dimension);
	

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
	//free(freq);
	free(frequencies );
	//free(detectors);
	delete [] detectors;
	free(data_length);
}
void test7()
{
	int dimension = 2;
	double initial_pos[2]={1,0.};
	double *seeding_var = NULL;

	
	int N_steps = 5000;
	int chain_N= 10;
	double ***output;
	output = allocate_3D_array( chain_N, N_steps, dimension );
	//double *initial_pos_ptr = initial_pos;
	int swp_freq = 3;
	//double chain_temps[chain_N] ={1,2,3,10,12};
	double chain_temps[chain_N];
	chain_temps[0] = 1.;
	double c = 1.5;
	//double temp_step = 1.2./(chain_N);
	for(int i =1; i < chain_N;  i ++)
		//chain_temps[i]=1.;
		chain_temps[i] =  chain_temps[i-1] * c;
	//double chain_temps[chain_N] ={1};
	//std::string autocorrfile = "testing/data/neil_auto_corr_mcmc.csv";
	std::string autocorrfile = "";
	std::string chainfile = "testing/data/neil_mcmc_output.csv";
	std::string statfilename = "testing/data/neil_mcmc_statistics.txt";
	std::string checkpointfile = "testing/data/neil_mcmc_checkpoint.csv";
	
	int numThreads = 1;
	bool pool = false;
	
	PTMCMC_MH(output, dimension, N_steps, chain_N, initial_pos,seeding_var,chain_temps, swp_freq, test_lp_nts, log_neil_proj3_nts,fisher_neil_proj3,(void **)NULL,numThreads, pool,show_progress, statfilename,chainfile,autocorrfile, "",checkpointfile );	
	//PTMCMC_MH(output, dimension, N_steps, chain_N, initial_pos,seeding_var,chain_temps, swp_freq, test_lp_nts, log_neil_proj3_nts,NULL,numThreads, pool,show_progress, statfilename,chainfile,autocorrfile, "",checkpointfile );	
	std::cout<<"ENDED"<<std::endl;

	//autocorrfile = "testing/data/neil_auto_corr_mcmc2.csv";
	//chainfile = "testing/data/neil_mcmc_output2.csv";
	//statfilename = "testing/data/neil_mcmc_statistics2.txt";
	//MCMC_MH(output, dimension, N_steps, chain_N, initial_pos,chain_temps, swp_freq, test_lp, log_neil_proj32,NULL,statfilename,chainfile,autocorrfile );	
	//std::cout<<"ENDED"<<std::endl;
	


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
	int length = 4000;
	double chirpm = 49.78;
	double eta =.21;
	//params.mass1 = calculate_mass1(chirpm,eta);
	//params.mass2 = calculate_mass2(chirpm,eta);
	params.mass1 = 14.9 ;
	params.mass2 = 8.4;
	string method= "IMRPhenomD";
	//string method= "ppE_IMRPhenomD_Inspiral";
	double amp[length];
	double phaseout[length];
	complex<double> waveformout[length];
	params.spin1[0] = 0;
	params.spin1[1] = 0;
	params.spin1[2] = 0.28221897810218977;
	params.spin2[0] = 0;
	params.spin2[1] = 0;
	params.spin2[2] = 0;
	double *spin1  = params.spin1;
	double *spin2= params.spin2;
	params.shift_phase=false;
	params.shift_time=false;
	params.phic = .0;
	params.tc = -.0;
	params.Luminosity_Distance = 334.;
	params.betappe = new double[1] ;
	params.betappe[0]=00.;
	params.bppe  =new int[1];
	params.bppe[0] = -3;
	params.Nmod = 1;
	params.NSflag1 = false;
	params.NSflag2 = false;
	params.phi = 0;
	params.theta = 0;
	params.incl_angle = 0;
	params.sky_average=true;
	//params.f_ref = 100;
	//params.phiRef = 1.0;
	
	double fhigh =1000;
	double flow =10;
	double df = (fhigh-flow)/(length-1);
	double *freq = (double *)malloc(sizeof(double) * length);

	cout<<"Freq spacing "<<df<<endl;

	for(int i=0;i<length;i++)
		freq[i]=flow+i*df;

	int dimension = 8;
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
	fisher_numerical(freq, length, "MCMC_IMRPhenomD_single_detect","Hanford_O1_fitted", 
			output3, dimensionmcmc, &params ,2);

	end7 = clock();

	cout<<"TIMING: FISHER: "<<(double)(end7-start7)/CLOCKS_PER_SEC<<endl;

	cout.precision(15);
	for (int i = 0;i <dimensionmcmc;i++)
	{
		for (int j=0;j <dimensionmcmc; j++){
			cout<<output3[i][j]<<"   ";
		}
		cout<<endl;
	}
	
	start7 = clock();
	fisher_numerical(freq, length, "dCS_IMRPhenomD","Hanford_O1_fitted", output, dimension, 
				&params ,2,NULL,NULL,NULL);

	end7 = clock();
	cout<<"TIMING: FISHER: "<<(double)(end7-start7)/CLOCKS_PER_SEC<<endl;
	cout.precision(15);
	for (int i = 0;i <dimension;i++)
	{
		for (int j=0;j <dimension; j++)
			cout<<output[i][j]<<"   ";
		cout<<endl;
	}
	start7 = clock();
	fisher_autodiff(freq, length, "dCS_IMRPhenomD","Hanford_O1_fitted", output2, dimension, 
				&params );

	end7 = clock();
	cout<<"TIMING: FISHER autodiff: "<<(double)(end7-start7)/CLOCKS_PER_SEC<<endl;
	for (int i = 0;i <dimension;i++)
	{
		for (int j=0;j <dimension; j++){
			cout<<output2[i][j]<<"   ";
			output2[i][j]*=(3);
		}
		cout<<endl;
	}
	std::cout<<"#########################################"<<std::endl<<"INVERSE"<<std::endl;
	double **inverse = allocate_2D_array(dimension,dimension);
	gsl_LU_matrix_invert(output2,inverse, dimension);
	for (int i = 0;i <dimension;i++)
	{
		//for (int j=0;j <dimension; j++)
		cout<<sqrt(inverse[i][i])<<"   ";
		cout<<endl;
	}
	deallocate_2D_array(inverse, dimension, dimension);
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
	params.NSflag1 = false;
	params.NSflag2 = false;
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
	fisher_numerical(freq, length, "ppE_IMRPhenomD_Inspiral","Hanford_O1_fitted", output, dimension, 
				&params ,2);

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
	params.NSflag1 = false;
	params.NSflag2 = false;
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
	allocate_FFTW_mem_forward(&plan,length);
	
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
		llnew = maximized_Log_Likelihood(waveformout, noise,freqnew,length, 
					&params,"Hanford","IMRPhenomPv2",&plan);
		//llnew = maximized_Log_Likelihood(real,imag, noise,freqnew,length, 
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

	deallocate_FFTW_mem(&plan);	
}
void test3()
{
	gen_params params;
	int length = 1000;
	params.mass1 = 200;
	params.mass2 = 50;
	string method= "IMRPhenomPv2";
	//string method= "IMRPhenomD";
	//string method= "ppE_IMRPhenomPv2_IMR";
	//string method= "IMRPhenomD";
	//string method= "ppE_IMRPhenomD_Inspiral";
	double amp[length];
	double phaseout[length];
	complex<double> waveformout_plus[length];
	complex<double> waveformout_cross[length];
	complex<double> response[length];
	params.spin1[0] = .0;
	params.spin1[1] = .01;
	params.spin1[2] = -.2;
	params.spin2[0] = .0;
	params.spin2[1] = 0.1;
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
	params.NSflag1 = false;
	params.NSflag2 = false;
	params.sky_average = false;
	//params.phi = M_PI/3.;
	//params.theta = M_PI/3;
	params.RA = 2;
	params.DEC = 1;
	params.gmst = 10;
	params.phi = 0;
	params.theta = 0;
	params.incl_angle = 0;
	params.phiRef = 10;
	params.f_ref = 15;
	
	double freq[length];
	for(int i=0;i<length;i++)
		freq[i]=10.+i*1e-1;

	clock_t  start, end;
	start = clock(); 
	//fourier_waveform(freq, length, waveformout_plus,waveformout_cross,method,&params);
	fourier_detector_response_equatorial(freq, length, response,"Hanford",method,&params);
	end=clock();
	//for(int i = 0 ;i<length; i++){
	//	std::cout<<waveformout_plus[i]<<" "<<waveformout_cross[i]<<std::endl;
	//}
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
	gsl_spline *Z_DL_spline_ptr;
	gsl_interp_accel *Z_DL_accel_ptr;
	initiate_LumD_Z_interp(&Z_DL_accel_ptr,&Z_DL_spline_ptr);
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
	params.betappe[0]=pow((  3e5/3.e5),4);
	params.bppe  =new int[1];
	params.bppe[0] = -1;
	params.Nmod = 1;
	params.NSflag1 = false;
	params.NSflag2 = false;
	params.phi = 0;
	params.theta = 0;
	params.incl_angle = 0;
	params.sky_average=true;
	params.Z_DL_accel_ptr = Z_DL_accel_ptr;
	params.Z_DL_spline_ptr = Z_DL_spline_ptr;
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
	//snr = data_snr_maximized_extrinsic(freq,length,waveformout,"Hanford_O1_fitted","IMRPhenomD",
	//		&params);
	//cout<<" FULL SNR: "<<snr<<endl;

	double logl;
	fftw_outline plan;
	clock_t  start4, end4;
	allocate_FFTW_mem_forward(&plan,length);
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
	double logl3 = maximized_Log_Likelihood(waveformout, noise,freq,length, 
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
	deallocate_FFTW_mem(&plan);	



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
	method = "dCS_IMRPhenomD";
	//method = "EdGB_IMRPhenomD_log";
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
	fisher_numerical(freq, length, "IMRPhenomD","Hanford_O1_fitted", output, dimension, 
				&params ,2);

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
	//fisher(freq, length, "ppE_IMRPhenomD_IMR","Hanford_O1_fitted", outputppe, dimensionppe, 
	//			&params );

	//end7 = clock();
	//cout<<"TIMING: FISHER ppE: "<<(double)(end7-start7)/CLOCKS_PER_SEC<<endl;
	//for (int i = 0;i <dimensionppe;i++)
	//{
	//	for (int j=0;j <dimensionppe; j++)
	//		cout<<outputppe[i][j]<<"   ";
	//	cout<<endl;
	//}


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
	free_LumD_Z_interp(&Z_DL_accel_ptr, &Z_DL_spline_ptr);
}

void fisher_neil_proj3 (double *pos,int dimension, double **fisher,int chain_id,void *parameters)
{
	
	
	int alpha = (int)(gsl_rng_uniform(g)*1e7);
	
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
	//fisher[0][0]=10;
	//fisher[1][0]=0;
	//fisher[0][1]=0;
	//fisher[1][1]=10;

}
adouble dist(adouble *pos, int dimension){
        adouble x = pos[0];
        adouble y = pos[1];
        adouble exponent_1 = - pow(x,2) - pow(9 + 4*pow(x,2) + 8*y , 2);
        adouble exponent_2 = - 8*pow(x,2) - 8*pow(y - 2, 2);
        adouble out =( 16/(3 * M_PI) ) * ( exp(exponent_1) + 0.5 * exp(exponent_2) ); 
 
        return out;
}
double log_neil_proj3_nts (double *c,int dim, int chainid,void *parameters)
{
	return log_neil_proj3(c,dim,parameters);
	//return 2;
}
double log_neil_proj3 (double *c,int dim,void *parameters)
{
	double x = c[0];
	double y = c[1];
	double prefactor = 16./(M_PI*3.);
	double pow1 = -x*x - pow((9+4*x*x +8*y),2);
	double pow2 = -8*x*x -8*pow(y-2,2);
	return log(prefactor*(std::exp(pow1) + .5*std::exp(pow2)));
	//return 2.;
}
double log_neil_proj32 (double *c,int dim,void *parameters)
{
	double x = c[0]*c[1];
	double y = c[0]/c[1];
	double prefactor = 16./(M_PI*3.);
	double pow1 = -x*x - pow((9+4*x*x +8*y),2);
	double pow2 = -8*x*x -8*pow(y-2,2);
	return log(prefactor*(std::exp(pow1) + .5*std::exp(pow2)));
	//return 2.;
}
double log_student_t (double *x,int dim,void *parameters){

	double  mu=1, nu=3,  sigma=1;
        double g1 = gsl_sf_gamma( (nu + 1) / 2 ) ;
        double g2 = gsl_sf_gamma( (nu/2) );
        double parenth = 1 + (1/nu) *pow( (x[0] - mu) / sigma, 2 );
        return log(g1 / (g2 * sqrt(nu * M_PI ) * sigma ) * pow(parenth,-(nu+1)/2    ));
}
void test_fisher(double *pos, int dim, double **fisher,void *parameters)
{
	fisher[0][0] = .5;
}
double test_ll(double *pos, int dim,void *parameters)
{
	//std::cout<<"LL"<<std::endl;
	//std::cout<<"Pos in LL: "<<pos[0]<<std::endl;
	return -pos[0]*pos[0]/(4.);
	//return  0;
}
double test_lp_nts(double *pos, int dim, int chain_id,void *parameters){
	return test_lp(pos,dim,parameters);
}
double test_lp(double *pos, int dim,void *parameters)
{
	double a = -std::numeric_limits<double>::infinity();
	if(pos[0]<-10 || pos[0]>10){return a;}
	else if(pos[1]<-10 || pos[1]>10){return a;}
	return 0;
	//return -pos[0]*pos[0]/(10.)- pos[1]*pos[1]/20.;
}	
double test_lp_GW(double *pos, int dim, int chain_id,void *parameters)
{
	double a = -std::numeric_limits<double>::infinity();
	//Flat priors across physical regions
	//if (std::exp(pos[0])/MPC_SEC<50 || std::exp(pos[0])/MPC_SEC>1000){return a;}
	if (std::exp(pos[0])<2 || std::exp(pos[0])>100){return a;}
	else if ((pos[1])<.1 || (pos[1])>.249999){return a;}
	else if ((pos[2])<-.9 || (pos[2])>.9){return a;}
	else if ((pos[3])<-.9 || (pos[3])>.9){return a;}
	//else {return 0.;}
	else {return pos[0] ;}
	//else {return log(std::exp(pos[0])/MSOL_SEC)-(std::exp(pos[0])/MSOL_SEC-30)*(std::exp(pos[0])/MSOL_SEC-30)/(2*10);}
	//else {return log(std::exp(pos[0])/MSOL_SEC)-(std::exp(pos[0])/MSOL_SEC-30)*(std::exp(pos[0])/MSOL_SEC-30)/(2*10)-(pos[1]-.24)*(pos[1]-.24)/(2*.010);}
}
double test_lp_GW_D(double *pos, int dim, int chain_id,void *parameters)
{
	double a = -std::numeric_limits<double>::infinity();
	if ((pos[0])<0 || (pos[0])>2*M_PI){return a;}//RA
	if ((pos[1])<-1 || (pos[1])>1){return a;}//sinDEC

	//if ((pos[2])<0 || (pos[2])>2*M_PI){return a;}//PSI
	if ((pos[2])<0 || (pos[2])>M_PI){return a;}//PSI
	if ((pos[3])<-1 || (pos[3])>1){return a;}//cos \iota
	//if ((pos[3])<0 || (pos[3])>1){return a;}//cos \iota

	if ((pos[4])<0 || (pos[4])>2*M_PI){return a;}//phiRef
	if ((pos[5])<0 || (pos[5])>10){return a;}//tc
	if (std::exp(pos[6])<10 || std::exp(pos[6])>10000){return a;}//DL
	if (std::exp(pos[7])<2 || std::exp(pos[7])>100 ){return a;}//chirpmass
	if ((pos[8])<.1 || (pos[8])>.249999){return a;}//eta
	if ((pos[9])<-.9 || (pos[9])>.9){return a;}//chi1 
	if ((pos[10])<-.9 || (pos[10])>.9){return a;}//chi2
	//else {return log(-sin(pos[0]))+pos[4]+pos[3];}
	//else { return pos[7]+3*pos[6]+log(abs(cos(pos[1]))) + log( abs( sin(pos[3])));}
	else { return pos[7]+3*pos[6] ;}
	//else {return pos[7]+3*pos[6];}
	//else {return pos[4]+3*pos[3] +std::log(std::abs(std::cos(pos[2])))
	//	+std::log(std::abs(std::cos(pos[10])))+std::log(std::abs(std::cos(pos[11])));}

}
double test_lp_GW_Pv2(double *pos, int dim, int chain_id,void *parameters)
{
	double a = -std::numeric_limits<double>::infinity();
	//if(sqrt(pos[9]*pos[9] + pos[10]*pos[10]+pos[11]*pos[11]) >0.95) {return a;}
	double chirp = std::exp(pos[7]);
	//double eta = pos[8];
	//double m1 = calculate_mass1(chirp,eta );
	//double m2 = calculate_mass2(chirp,eta );
	//double q =m1/m2;
	//double W = (3*q +4)/ ( 4*q*q +3*q);
	////Max values
	//double chi1l = pos[9];
	//double chi2l = pos[10];
	//double chi1p = std::sqrt(1- chi1l*chi1l);
	//double chi2p = std::sqrt(1- chi2l*chi2l);
	//double chi_thresh=chi2p ;
	//if(chi1p > W*chi2p){ chi_thresh =chi1p;}
	//if(pos[11] > chi_thresh){ return a;}

	//Flat priors across physical regions
	//if ((pos[0])<0 || (pos[0])>M_PI){return a;}
	if ((pos[0])<0 || (pos[0])>2*M_PI){return a;}//RA
	if ((pos[1])<-1 || (pos[1])>1){return a;}//sinDEC

	//if ((pos[2])<0 || (pos[2])>2*M_PI){return a;}//PSI
	if ((pos[2])<0 || (pos[2])>M_PI){return a;}//PSI
	if ((pos[3])<-1 || (pos[3])>1){return a;}//cos \iota
	//if ((pos[3])<0 || (pos[3])>1){return a;}//cos \iota

	if ((pos[4])<0 || (pos[4])>2*M_PI){return a;}//PhiRef
	if ((pos[5])<0 || (pos[5])>10){return a;}//PhiRef
	if (std::exp(pos[6])<10 || std::exp(pos[6])>10000){return a;}//DL
	if (std::exp(pos[7])<2 || std::exp(pos[7])>100 || std::isnan(pos[4])){return a;}//chirpmass
	if ((pos[8])<.1 || (pos[8])>.249999){return a;}//eta
	if ((pos[9])<.0 || (pos[9])>.9){return a;}//chi1 
	if ((pos[10])<.0 || (pos[10])>.9){return a;}//chi2
	if ((pos[11])<-1 || (pos[11])>1){return a;}//chip
	if ((pos[12])<-1 || (pos[12])>1){return a;}//phip
	if ((pos[13])<0 || (pos[13])>2*M_PI){return a;}//phip
	//else {return log(-sin(pos[0]))+pos[4]+pos[3];}
	else {return 3*pos[6]+log(chirpmass_eta_jac(chirp,pos[8]));}
	//else {return pos[7]+3*pos[6] ;}
	//else {return pos[4]+3*pos[3] +std::log(std::abs(std::cos(pos[2])))
	//	+std::log(std::abs(std::cos(pos[10])))+std::log(std::abs(std::cos(pos[11])));}
}
//double test_lp_GW_Pv2(double *pos, int dim, int chain_id)
//{
//	double a = -std::numeric_limits<double>::infinity();
//	//Flat priors across physical regions
//	//if ((pos[0])<0 || (pos[0])>M_PI){return a;}
//	if ((pos[0])<-1 || (pos[0])>1){return a;}//cos \iota
//	if ((pos[1])<0 || (pos[1])>2*M_PI){return a;}//RA
//	if ((pos[2])<-M_PI/2. || (pos[2])>M_PI/2.){return a;}//DEC
//	if (std::exp(pos[3])<10 || std::exp(pos[3])>10000){return a;}//DL
//	if (std::exp(pos[4])<2 || std::exp(pos[4])>100 || std::isnan(pos[4])){return a;}//chirpmass
//	if ((pos[5])<.1 || (pos[5])>.249999){return a;}//eta
//	if ((pos[6])<0 || (pos[6])>.9){return a;}//chi1 
//	if ((pos[7])<0 || (pos[7])>.9){return a;}//chi2
//	if ((pos[8])<0 || (pos[8])>M_PI){return a;}//theta1
//	if ((pos[9])<0 || (pos[9])>M_PI){return a;}//theta2
//	if ((pos[10])<0 || (pos[10])>2*M_PI){return a;}//phi1
//	if ((pos[11])<0 || (pos[11])>2*M_PI){return a;}//phi2
//	if ((pos[12])<0 || (pos[12])>2*M_PI){return a;}//phiRef
//	if ((pos[13])<0 || (pos[13])>2*M_PI){return a;}//polarization angle
//	//else {return log(-sin(pos[0]))+pos[4]+pos[3];}
//	//else {return pos[4]+3*pos[3];}
//	else {return pos[4]+3*pos[3] +std::log(std::abs(std::cos(pos[2])))
//		+std::log(std::abs(std::cos(pos[10])))+std::log(std::abs(std::cos(pos[11])));}
//}
double test_lp_GW_Pv2_ppE(double *pos, int dim, int chain_id,void *parameters)
{
	double a = -std::numeric_limits<double>::infinity();
	if ((pos[14])<-100 || (pos[14])>100){return a;}//ppE beta
	else{return test_lp_GW_Pv2(pos,dim,chain_id,parameters);}
}
double test_lp_GW_Pv2_dCS_root_alpha(double *pos, int dim, int chain_id,void *parameters)
{
	double a = -std::numeric_limits<double>::infinity();
	if ((pos[14])<0 || (pos[14])>100){return a;}//ppE beta
	else{return test_lp_GW_Pv2(pos,dim,chain_id,parameters);}
}
double test_lp_GW_7dim(double *pos, int dim, int chain_id,void *parameters)
{
	double a = -std::numeric_limits<double>::infinity();
	//Flat priors across physical regions
	if (std::exp(pos[0])<50 || std::exp(pos[0])>1000){return a;}
	if ((pos[1])<0 || (pos[1])>4){return a;}
	if ((pos[2])<0 || (pos[2])>2*M_PI){return a;}
	if (std::exp(pos[3])<2 || std::exp(pos[3])>100){return a;}
	if ((pos[4])<.1 || (pos[4])>.245){return a;}
	if ((pos[5])<-.9 || (pos[5])>.9){return a;}
	if ((pos[6])<-.9 || (pos[6])>.9){return a;}
	//else {return 0.;}
	//else {return log(std::exp(pos[3])*std::exp(pos[0])*std::exp(pos[0])*std::exp(pos[0]));}
	else {return pos[3]+4*pos[0];}
}
double test_lp_GW_DFull(double *pos, int dim, int chain_id,void *parameters)
{
	double a = -std::numeric_limits<double>::infinity();
	//Flat priors across physical regions
	//if ((pos[0])<0 || (pos[0])>M_PI){return a;}
	if ((pos[0])<-1 || (pos[0])>1){return a;}//cos \iota
	if ((pos[1])<0 || (pos[1])>2*M_PI){return a;}//RA
	if ((pos[2])<-M_PI/2. || (pos[2])>M_PI/2.){return a;}//DEC
	if (std::exp(pos[3])<10 || std::exp(pos[3])>10000){return a;}//DL
	if (std::exp(pos[4])<2 || std::exp(pos[4])>100){return a;}//chirpmass
	if ((pos[5])<.1 || (pos[5])>.249999){return a;}//eta
	if ((pos[6])<-.9 || (pos[6])>.9){return a;}//chi1 
	if ((pos[7])<-.9 || (pos[7])>.9){return a;}//chi2
	//else {return log(-sin(pos[0]))+pos[4]+pos[3];}
	//else {return pos[4]+3*pos[3];}
	else {return pos[4]+3*pos[3] +std::log(std::abs(std::cos(pos[2])));}
	//else {return pos[4]+1*pos[3];}
}
double test_lp_GW_dCS(double *pos, int dim , int chain_id,void *parameters)
{
	double a = -std::numeric_limits<double>::infinity();
	double alphasqmin = 0;
	double alphasqmax = pow(1.e2/(3e5),4);//>10^6 km in \sqrt{\alpha}
	
	//if((std::exp(pos[8])<alphamin) || std::exp(pos[8])>alphamax){return a;}
	//if((std::exp(pos[8])<alphamin) ){return a;}
	//if(((pos[8])<lnalphamin) ){return a;}
	//if(((pos[8])<alphasqmin )|| (pos[8]>alphasqmax  )){std::cout<<"BOUDNARY"<<std::endl;return a;}
	if(((pos[8])<alphasqmin )|| (pos[8]>alphasqmax  )){return a;}
	
	//Uniform prior on \alpha^2, not \alpha
	else { return test_lp_GW_DFull(pos,dim,chain_id,parameters)-.75*std::log(pos[8]);}
	//else { return test_lp_GW_DFull(pos,dim,chain_id);}
}
double test_lp_GW_dCS_root_alpha(double *pos, int dim , int chain_id,void *parameters)
{
	double a = -std::numeric_limits<double>::infinity();
	double alphasqmin = 0;
	double alphasqmax = 100;//>10^6 km in \sqrt{\alpha}
	
	//if((std::exp(pos[8])<alphamin) || std::exp(pos[8])>alphamax){return a;}
	//if((std::exp(pos[8])<alphamin) ){return a;}
	//if(((pos[8])<lnalphamin) ){return a;}
	//if(((pos[8])<alphasqmin )|| (pos[8]>alphasqmax  )){std::cout<<"BOUDNARY"<<std::endl;return a;}
	if(((pos[8])<alphasqmin )|| (pos[8]>alphasqmax  )){return a;}
	
	//Uniform prior on \alpha^2, not \alpha
	else { return test_lp_GW_DFull(pos,dim,chain_id,parameters);}
	//else { return test_lp_GW_DFull(pos,dim,chain_id);}
}
double test_lp_GW_dCS_log(double *pos, int dim , int chain_id,void *parameters)
{
	double a = -std::numeric_limits<double>::infinity();
	//double lnalphamin = -500;//Corresponds to a limit on root alpha of 10^-46 m
	//double lnalphamin = -100;//Corresponds to a limit on root alpha of 4 mm
	double lnalphamin = -70;//Corresponds to a limit on root alpha of 7.5 m
	double lnalphamax = -25;//Corresponds to 500km
	//if((std::exp(pos[8])<alphamin) || std::exp(pos[8])>alphamax){return a;}
	//if((std::exp(pos[8])<alphamin) ){return a;}
	//if(((pos[8])<lnalphamin) ){return a;}
	if(((pos[8])<lnalphamin)|| (pos[8]>lnalphamax) ){return a;}
	
	//Uniform prior on \alpha^.5, not \alpha
	else { return test_lp_GW_DFull(pos,dim,chain_id,parameters)+.25*pos[8];}
	//Uniform in ln \alpha^2
	//else { return test_lp_GW_DFull(pos,dim,chain_id);}
}
double test_lp_GW_ppE(double *pos, int dim , int chain_id,void *parameters)
{
	double a = -std::numeric_limits<double>::infinity();
	//if(((pos[8])<lnalphamin)|| (pos[8]>lnalphamax) ){return a;}
	//Uniform prior on \alpha^.5, not \alpha
	//else { return test_lp_GW_DFull(pos,dim,chain_id)+.25*pos[8];}
	//Uniform in ln \alpha^2
	{ return test_lp_GW_DFull(pos,dim,chain_id,parameters);}
	//else { return test_lp_GW_DFull(pos,dim,chain_id);}
}
double test_lp_reduceddim(double *pos, int dim, int chain_id,void *parameters)
{
	double a = -std::numeric_limits<double>::infinity();
	if (std::exp(pos[0])<2 || std::exp(pos[0])>100){return a;}//chirpmass
	if ((pos[1])<.1 || (pos[1])>.249999){return a;}//eta
	if ((pos[2])<-.9 || (pos[2])>.9){return a;}//chi1 
	if ((pos[3])<-.9 || (pos[3])>.9){return a;}//chi2
	else {return pos[0];}
	//else {return pos[4]+1*pos[3];}
}
double test_lp_GW_D_in(double *pos, int dim, int chain_id,void *parameters)
{
	double a = -std::numeric_limits<double>::infinity();
	double chirp = std::exp(pos[0]);
	double eta = pos[1];
	if (std::exp(pos[0])<2 || std::exp(pos[0])>100){return a;}//chirpmass
	if ((pos[1])<.1 || (pos[1])>.249999){return a;}//eta
	if ((pos[2])<-.99 || (pos[2])>.99){return a;}//chi1 
	if ((pos[3])<-.99 || (pos[3])>.99){return a;}//chi2
	else {return log(chirpmass_eta_jac(chirp,eta));}
}
double test_lp_GW_Pv2_in(double *pos, int dim, int chain_id,void *parameters)
{
	double a = -std::numeric_limits<double>::infinity();
	double chirp = std::exp(pos[0]);
	double eta = pos[1];
	if (std::exp(pos[0])<2 || std::exp(pos[0])>100){return a;}//chirpmass
	if ((pos[1])<.1 || (pos[1])>.249999){return a;}//eta
	if ((pos[2])<0 || (pos[2])>.99){return a;}//chi1 
	if ((pos[3])<0 || (pos[3])>.99){return a;}//chi2
	if ((pos[4])<-1 || (pos[4])>1){return a;}//chi1 
	if ((pos[5])<-1 || (pos[5])>1){return a;}//chi2
	if ((pos[6])<0 || (pos[6])>2*M_PI){return a;}//chi2
	else {return log(chirpmass_eta_jac(chirp,eta));}
}
//Uniform in m1 and m2, transformed to lnM and eta
double chirpmass_eta_jac(double chirpmass, double eta){
	return chirpmass*chirpmass/(sqrt(1. - 4.*eta)*pow(eta,1.2));
}
