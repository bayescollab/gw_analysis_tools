#include "mcmc_gw.h"
#include "waveform_generator.h"
#include "io_util.h"
#include "util.h"
#include <iostream>
#include <unordered_map>
#include <string>
#include <limits>

double T_mcmc_gw_tool;
/*! \file 
 *
 * Command line tool for analyzing LOSC GW data
 *
 * See mcmc_sampler.cpp documentation for more explanation on the sampler
 *
 * See mcmc_gw.cpp for GW specific explanation
 *
 * Usage:
 * 	
 * 	mcmc_gw_tool_basic /PATH/TO/PARAM/FILE
 *
 * See data/sample_config_files/mcmc_gw_tool_param_template.dat for an example parameter file
 *
 * See data/sample_init_pos.csv for an example initial position file
 *
 * NOTE: SkySearch generation method requires an initialization file with 11 dimensions. 
 * The extra parameters correspond to the injected intrinsic parameters -- exactly the format of PhenomD
 */

double **beta_priors;
double standard_log_prior_D(double *pos, int dim, int chain_id,void *parameters);
double standard_log_prior_Pv2(double *pos, int dim, int chain_id,void *parameters);
double standard_log_prior_D_intrinsic(double *pos, int dim, int chain_id,void *parameters);
double standard_log_prior_Pv2_intrinsic(double *pos, int dim, int chain_id,void *parameters);
double standard_log_prior_skysearch(double *pos, int dim, int chain_id, void *parameters);
double standard_log_prior_Pv2_ppE(double *pos, int dim, int chain_id,void *parameters);
double chirpmass_eta_jac(double m1,double m2);
int main(int argc, char *argv[])
{
	std::cout.precision(15);
	if(argc !=2){
		std::cout<<"ERROR -- A parameter file is required"<<std::endl;
		return 1;
	}
	std::string param_file(argv[1]);

	std::unordered_map<std::string, int> int_dict;
	std::unordered_map<std::string, double> dbl_dict;
	std::unordered_map<std::string, float> flt_dict;
	std::unordered_map<std::string, bool> bool_dict;
	std::unordered_map<std::string, std::string> str_dict;
	std::string input_param_file(argv[1]);
	int status = unpack_input_io_file(input_param_file, &int_dict, &str_dict, &dbl_dict, &flt_dict, &bool_dict);

	//Global
	int detector_N = int_dict["detector number"];
	std::cout<<"DETECTOR NUMBER: "<<detector_N<<std::endl;
	int samples = int_dict["samples"];
	std::cout<<"Samples: "<<samples<<std::endl;
	int max_thermo_chain_N = int_dict["max thermo chain number"];
	std::cout<<"Max thermo ensemble chain number: "<<max_thermo_chain_N<<std::endl;
	int t0 = int_dict["t0"];
	std::cout<<"t0: "<<t0<<std::endl;
	int nu = int_dict["nu"];
	std::cout<<"nu: "<<nu<<std::endl;
	
	std::string detectors[detector_N];
	std::string detector_files[detector_N];
	std::string psd_file = str_dict["PSD filepath"];
	std::cout<<"PSD file: "<<psd_file<<std::endl;
	for(int i = 0 ; i<detector_N; i++){
		detectors[i]= str_dict["detector name "+std::to_string(i)];
		detector_files[i]= str_dict["data file "+std::to_string(i)];
	}
	std::string generation_method = str_dict["generation method"];
	
	double data_time_length = dbl_dict["data length"];
	std::cout<<"data length: "<<data_time_length<<std::endl;
	int data_length=131072;
	//if((int)data_time_length == 32){
	//	data_length = 131072;
	//}	
	count_lines_LOSC_data_file(detector_files[0], &data_length);
	std::cout<<"data length (lines): "<<data_length<<std::endl;
	double gps_time = dbl_dict["gps"];
	std::cout<<"GPS time : "<<gps_time<<std::endl;
	std::string allocation_scheme = str_dict["allocation scheme"];
	std::cout<<"Allocation: "<<allocation_scheme<<std::endl;
	int swap_freq = int_dict["swap frequency"];
	std::cout<<"swap_freq: "<<swap_freq<<std::endl;
	int threads = int_dict["thread number"];
	std::cout<<"threads: "<<threads<<std::endl;
	int correlation_thresh = int_dict["correlation threshold"];
	int correlation_segs = int_dict["correlation segments"];
	double correlation_convergence_thresh = dbl_dict["correlation convergence threshold"];
	double ac_target = dbl_dict["autocorrelation target"];
	std::string output_file = str_dict["output data file"];
	std::string stat_file = str_dict["output stat file"];
	std::string check_file = str_dict["checkpoint file"];
	int dimension = int_dict["dimension"];
	
	std::string initial_position_file, initial_checkpoint_file;
	int chain_N;
	bool continue_from_checkpoint=false;
	if(str_dict.find("initial checkpoint file") == str_dict.end()){
		//Original
		initial_position_file = str_dict["initial position file"];
		chain_N = int_dict["chain number"];
		std::cout<<"Chain number: "<<chain_N<<std::endl;
	}
	else{
		//Continue	
		continue_from_checkpoint=true;
		initial_checkpoint_file = str_dict["initial checkpoint file"];
		std::cout<<"INITIAL checkpoint file: "<<initial_checkpoint_file<<std::endl;
		
	}
	

	
	int psd_length ;
	count_lines_LOSC_PSD_file(psd_file, &psd_length);
	std::cout<<"Length of PSD: "<<psd_length<<std::endl;



	double **psd = allocate_2D_array(detector_N,psd_length);
	double **freqs = allocate_2D_array(detector_N,psd_length);
	std::complex<double> **data = (std::complex<double> **)malloc(sizeof(std::complex<double> *)*detector_N);
	for(int i =0; i<detector_N; i++){
		data[i] = (std::complex<double>*)malloc(sizeof(std::complex<double>)*psd_length);
	}

	allocate_LOSC_data(detector_files, psd_file, detector_N, psd_length, data_length, gps_time, data, psd, freqs);

	std::cout<<"DATA loaded"<<std::endl;

	T_mcmc_gw_tool = 1./(freqs[0][2]-freqs[0][1]);
	double df = 1./T_mcmc_gw_tool;

	std::cout<<"Total time: "<<T_mcmc_gw_tool<<std::endl;
	int *data_lengths= (int*)malloc(sizeof(int)*detector_N);
	for(int i = 0 ; i<detector_N; i++){
		data_lengths[i] =psd_length;
	}


	//#############################################################
	//#############################################################
	
	double **whitened = allocate_2D_array(data_lengths[0],7);
	for(int i = 0 ; i<data_lengths[0]; i++){
		whitened[i][0]= freqs[0][i];
		whitened[i][1]=psd[0][i];
		whitened[i][2]=psd[1][i];
		whitened[i][3]=real(data[0][i]);
		whitened[i][4]=imag(data[0][i]);
		whitened[i][5]=real(data[1][i]);
		whitened[i][6]=imag(data[1][i]);
	}	
	write_file("data/whitened_data.csv",whitened,data_lengths[0],7);
	deallocate_2D_array(whitened, data_lengths[0],7);
	
	//#############################################################
	//#############################################################

	//########################################################################
	double **output;
	output = allocate_2D_array(samples, dimension );
	double chain_temps[chain_N];
	
	int Nmod = 0;
	int *bppe = NULL;
	std::cout<<"Generation method: "<<generation_method<<std::endl;
	if(generation_method.find("ppE") != std::string::npos){
		Nmod = int_dict["Number of modifications"];
		std::cout<<"Number of ppE modifications: "<<Nmod<<std::endl;
		std::cout<<"ppE b parmeters: "<<Nmod<<std::endl;
		bppe= new int[Nmod];
		beta_priors = new double*[Nmod];
		for(int i =0; i<Nmod ; i++){
			bppe[i] = int_dict["ppE b "+std::to_string(i)];
			beta_priors[i]= new double[2];
			std::cout<<i<<" : "<<bppe[i]<<std::endl;
			beta_priors[i][0] = dbl_dict["ppE beta "+std::to_string(i)+" minimum"];
			beta_priors[i][1] = dbl_dict["ppE beta "+std::to_string(i)+" maximum"];
		}
		
	}
	bool pool = true;
	bool show_progress = true;

	//#########################################################
	if(generation_method.find("SkySearch") != std::string::npos){
		std::complex<double> *hplus=new std::complex<double>[data_lengths[0]];
		std::complex<double> *hcross=new std::complex<double>[data_lengths[0]];

		double **initial_position = new double*[1];
		initial_position[0] = new double[11];
		read_file(initial_position_file, initial_position,1,11);

		gen_params params ;
		for(int i = 0 ; i<11; i++){
			std::cout<<initial_position[0][i]<<std::endl;
		}
		params.mass1 = calculate_mass1(exp(initial_position[0][7]),initial_position[0][8]);
		params.mass2 = calculate_mass2(exp(initial_position[0][7]),initial_position[0][8]);
		params.spin1[2]=initial_position[0][9];
		params.spin2[2]=initial_position[0][10];
		params.phiRef = 0;
		params.tc = 0;
		params.f_ref =20;
		params.NSflag1=false;
		params.NSflag2=false;
		params.shift_phase=true;
		params.shift_time=true;
		params.equatorial_orientation=false;
		params.sky_average=false;
		params.Luminosity_Distance = exp(initial_position[0][6]);
		params.incl_angle = acos(initial_position[0][3]);

		fourier_waveform(freqs[0],data_lengths[0],hplus, hcross, "IMRPhenomD",&params);
	
	
		double *seeding_var = NULL;
		SkySearch_PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW(output, dimension, samples, chain_N, 
				max_thermo_chain_N, initial_position[0],seeding_var,chain_temps, 
				swap_freq, t0, nu, correlation_thresh, correlation_segs,
				correlation_convergence_thresh , ac_target,allocation_scheme, 
				standard_log_prior_skysearch,threads, pool,show_progress,detector_N, 
				data, psd,freqs, data_lengths,gps_time, detectors,Nmod, bppe,
				hplus,hcross,stat_file,output_file, "",check_file);	
		delete []  hplus;
		delete []  hcross;
		delete [] initial_position[0]; delete [] initial_position;
	}
	else{
	
		double(*lp)(double *param, int dimension, int chain_id, void *parameters);
		if(generation_method.find("IMRPhenomD") != std::string::npos && dimension == 11){
			lp = &standard_log_prior_D;
		}
		else if(generation_method.find("IMRPhenomPv2") != std::string::npos && dimension == 14){
			lp = &standard_log_prior_Pv2;
		}
		else if(generation_method.find("ppE_IMRPhenomPv2") != std::string::npos && dimension >= 14){
			lp = &standard_log_prior_Pv2_ppE;
		}
		else if(generation_method.find("IMRPhenomPv2") != std::string::npos && dimension == 7){
			lp = &standard_log_prior_Pv2_intrinsic;
		}
		else if(generation_method.find("IMRPhenomD") != std::string::npos && dimension == 4){
			lp = &standard_log_prior_D_intrinsic;
		}
		else{
			std::cout<<"ERROR -- wrong detector/dimension combination for this tool -- Check mcmc_gw for general support"<<std::endl;
			return 1;
		}

		if(continue_from_checkpoint){
			continue_PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW(initial_checkpoint_file,output, samples,  
					max_thermo_chain_N, chain_temps, 
					swap_freq, t0, nu, correlation_thresh, correlation_segs,
					correlation_convergence_thresh , ac_target,allocation_scheme, 
					lp,threads, pool,show_progress,detector_N, 
					data, psd,freqs, data_lengths,gps_time, detectors,Nmod, bppe,
					generation_method,stat_file,output_file, "",check_file);	

		}
		else{
			double **initial_position = new double*[1];
			initial_position[0] = new double[dimension];
			read_file(initial_position_file, initial_position,1,dimension);
			double *seeding_var = NULL;
			double ***output2 = allocate_3D_array(chain_N,samples, dimension );
			double c = 1.1;
			chain_temps[0]=1;
			for(int i = 1 ; i<chain_N; i++){
				chain_temps[i]=chain_temps[i-1]*c;
			}
			PTMCMC_MH_dynamic_PT_alloc_GW(output2, dimension, (int)5*t0, chain_N, max_thermo_chain_N, 
					initial_position[0],seeding_var,chain_temps, 
					swap_freq,t0, nu, allocation_scheme,lp,threads, pool,show_progress,detector_N, 
					data, psd,freqs, data_lengths,gps_time, detectors,Nmod, bppe,
					generation_method,stat_file,"", "",check_file);	
			continue_PTMCMC_MH_GW(check_file, output2, dimension, 5*t0,  
					swap_freq, lp,threads, pool,show_progress,detector_N, 
					data, psd,freqs, data_lengths,gps_time, detectors,Nmod, bppe,
					generation_method,stat_file,"", "","",check_file,true);	
			continue_PTMCMC_MH_GW(check_file, output2, dimension, samples,  
					swap_freq, lp,threads, pool,show_progress,detector_N, 
					data, psd,freqs, data_lengths,gps_time, detectors,Nmod, bppe,
					generation_method,stat_file,output_file, "","",check_file,true);	
			deallocate_3D_array(output2,chain_N,samples,dimension);
			delete [] initial_position[0]; delete [] initial_position;
		}

	}


	//write_file(chainfile, output[0], n_steps, dimension);

	deallocate_2D_array(output, samples, dimension);
	free(data_lengths);
	//free_LOSC_data(data, psd,freqs, detector_N, length);
	deallocate_2D_array(psd,detector_N, psd_length);
	deallocate_2D_array(freqs,detector_N, psd_length);
	for(int i =0; i<detector_N; i++)
		free(data[i]);
	free(data);

	if(generation_method.find("ppE") != std::string::npos){
		delete [] bppe;	
		for(int i = 0 ; i<Nmod ; i++){
			delete [] beta_priors[i] ;
		}
		delete [] beta_priors;
	}
	
	return 0;

}
double standard_log_prior_D(double *pos, int dim, int chain_id,void *parameters)
{
	double chirp = std::exp(pos[7]);
	double eta = pos[8];
	double a = -std::numeric_limits<double>::infinity();
	if ((pos[0])<0 || (pos[0])>2*M_PI){ return a;}//RA

	if ((pos[1])<-1 || (pos[1])>1){return a;}//sinDEC
	//if ((pos[1])<-M_PI/2 || (pos[1])>M_PI/2){return a;}//sinDEC

	if ((pos[2])<0 || (pos[2])>2*M_PI){return a;}//PSI
	if ((pos[3])<-1 || (pos[3])>1){return a;}//cos \iota
	if ((pos[4])<0 || (pos[4])>2*M_PI){return a;}//phiRef
	//if ((pos[5])<0 || (pos[5])>T_mcmc_gw_tool){return a;}//tc
	if ((pos[5])<T_mcmc_gw_tool*3./4. -.1 || (pos[5])>3.*T_mcmc_gw_tool/4. + .1){return a;}//tc
	if (std::exp(pos[6])<10 || std::exp(pos[6])>1000){return a;}//DL
	if (std::exp(pos[7])<2 || std::exp(pos[7])>80 ){return a;}//chirpmass
	if ((pos[8])<.1 || (pos[8])>.249999){return a;}//eta

	//if ((pos[9])<-.95 || (pos[9])>.95){return a;}//chi1 
	//if ((pos[10])<-.95 || (pos[10])>.95){return a;}//chi2
	double chi1 = pos[9]+pos[10];	
	double chi2 = pos[9]-pos[10];	
	if ((chi1)<-.95 || (chi1)>.95){return a;}//chi1 
	if ((chi2)<-.95 || (chi2)>.95){return a;}//chi2

	else {return log(chirpmass_eta_jac(chirp,eta))+3*pos[6] ;}
	//else {return log(chirpmass_eta_jac(chirp,eta))+3*pos[6] -log(cos(asin(pos[1]))) ;}

}
double standard_log_prior_Pv2(double *pos, int dim, int chain_id,void *parameters)
{
	double a = -std::numeric_limits<double>::infinity();
	double chirp = std::exp(pos[7]);
	double eta = pos[8];

	//Flat priors across physical regions
	if ((pos[0])<0 || (pos[0])>2*M_PI){return a;}//RA
	if ((pos[1])<-1 || (pos[1])>1){return a;}//sinDEC

	if ((pos[2])<0 || (pos[2])>2*M_PI){return a;}//PSI
	if ((pos[3])<-1 || (pos[3])>1){return a;}//cos \iota

	if ((pos[4])<0 || (pos[4])>2*M_PI){return a;}//PhiRef
	if ((pos[5])<0 || (pos[5])>T_mcmc_gw_tool){return a;}//PhiRef
	if (std::exp(pos[6])<10 || std::exp(pos[6])>10000){return a;}//DL
	if (std::exp(pos[7])<2 || std::exp(pos[7])>100 ){return a;}//chirpmass
	if ((pos[8])<.1 || (pos[8])>.249999){return a;}//eta
	if ((pos[9])<0 || (pos[9])>.95){return a;}//a1 
	if ((pos[10])<0 || (pos[10])>.95){return a;}//a2
	if ((pos[11])<-1 || (pos[11])>1){return a;}//theta1
	if ((pos[12])<-1 || (pos[12])>1){return a;}//theta2
	if ((pos[13])<0 || (pos[13])>2*M_PI){return a;}//phip
	else {return log(chirpmass_eta_jac(chirp,eta))+3*pos[6];}
}
double standard_log_prior_Pv2_ppE(double *pos, int dim, int chain_id,void *parameters)
{
	double a = -std::numeric_limits<double>::infinity();
	double chirp = std::exp(pos[7]);
	double eta = pos[8];

	//Flat priors across physical regions
	if ((pos[0])<0 || (pos[0])>2*M_PI){return a;}//RA
	if ((pos[1])<-1 || (pos[1])>1){return a;}//sinDEC

	if ((pos[2])<0 || (pos[2])>2*M_PI){return a;}//PSI
	if ((pos[3])<-1 || (pos[3])>1){return a;}//cos \iota

	if ((pos[4])<0 || (pos[4])>2*M_PI){return a;}//PhiRef
	if ((pos[5])<0 || (pos[5])>T_mcmc_gw_tool){return a;}//PhiRef
	//if ((pos[5])<T_mcmc_gw_tool*3./4. -.1 || (pos[5])>3.*T_mcmc_gw_tool/4. + .1){return a;}//tc
	if (std::exp(pos[6])<10 || std::exp(pos[6])>10000){return a;}//DL
	if (std::exp(pos[7])<2 || std::exp(pos[7])>100 ){return a;}//chirpmass
	if ((pos[8])<.1 || (pos[8])>.249999){return a;}//eta
	if ((pos[9])<0 || (pos[9])>.95){return a;}//a1 
	if ((pos[10])<0 || (pos[10])>.95){return a;}//a2
	if ((pos[11])<-1 || (pos[11])>1){return a;}//theta1
	if ((pos[12])<-1 || (pos[12])>1){return a;}//theta2
	if ((pos[13])<0 || (pos[13])>2*M_PI){return a;}//phip
	for(int i = 0 ; i<dim - 14; i++){
		if( pos[14+i] <beta_priors[i][0] || pos[14+i] >beta_priors[i][1]){return a;}
	}
	return log(chirpmass_eta_jac(chirp,eta))+3*pos[6];
}
double standard_log_prior_D_intrinsic(double *pos, int dim, int chain_id,void *parameters)
{
	double a = -std::numeric_limits<double>::infinity();
	double chirp = std::exp(pos[0]);
	double eta = pos[1];
	//Flat priors across physical regions
	if (exp(pos[0])<2 || exp(pos[0])>100){return a;}//RA
	if ((pos[1])<.1 || (pos[1])>.25){return a;}//sinDEC
	if ((pos[2])<-.95 || (pos[2])>.95){return a;}//chi1 
	if ((pos[3])<-.95 || (pos[3])>.95){return a;}//chi2
	else {return log(chirpmass_eta_jac(chirp,eta)) ;}


}
double standard_log_prior_Pv2_intrinsic(double *pos, int dim, int chain_id,void *parameters)
{
	double a = -std::numeric_limits<double>::infinity();
	double chirp = std::exp(pos[0]);
	double eta = pos[1];
	//Flat priors across physical regions
	if ((exp(pos[0]))<2 || exp(pos[0])>100){return a;}//RA
	if ((pos[1])<.1 || (pos[1])>.25){return a;}//sinDEC
	if ((pos[2])<0 || (pos[2])>.95){return a;}//chi1 
	if ((pos[3])<0 || (pos[3])>.95){return a;}//chi2
	if ((pos[4])<-1 || (pos[4])>1){return a;}//chi1 
	if ((pos[5])<-1 || (pos[5])>1){return a;}//chi2
	if ((pos[6])<0 || (pos[6])>2*M_PI){return a;}//chi2
	else {return log(chirpmass_eta_jac(chirp,eta)) ;}

}

double standard_log_prior_skysearch(double *pos, int dim, int chain_id, void *parameters){

	double a = -std::numeric_limits<double>::infinity();
	double DEC = asin(pos[1]);
	if ((pos[0])<0 || (pos[0])>2*M_PI){return a;}//RA
	if ((pos[1])<-1 || (pos[1])>1){return a;}//sinDEC

	if ((pos[2])<0 || (pos[2])>M_PI){return a;}//PSI
	if ((pos[3])<-1 || (pos[3])>1){return a;}//cos \iota

	if ((pos[4])<0 || (pos[4])>2*M_PI){return a;}//PhiRef
	if ((pos[5])<0 || (pos[5])>T_mcmc_gw_tool){return a;}//PhiRef
	if (std::exp(pos[6])<10 || std::exp(pos[6])>10000){return a;}//DL
	else { return 3*pos[6];}
}

//Uniform in m1 and m2, transformed to lnM and eta
double chirpmass_eta_jac(double chirpmass, double eta){
	return chirpmass*chirpmass/(sqrt(1. - 4.*eta)*pow(eta,1.2));
}
