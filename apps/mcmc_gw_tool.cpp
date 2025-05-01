#include "mcmc_gw.h"
#include "waveform_generator.h"
#include "ppE_utilities.h"
#include "io_util.h"
#include "util.h"
#include "EA_IMRPhenomD_NRT.h"
#include <iostream>
#include <unordered_map>
#include <string>
#include <limits>
#include <eigen3/Eigen/Eigen>


double T_mcmc_gw_tool;
double T_merger; 
/*! \file 
 *
 * Command line tool for analyzing LOSC GW data
 *
 * Runs the ``uncorrelated'' MCMC sampler with a generic prior
 *
 * See mcmc_sampler.cpp documentation for more explanation on the sampler
 *
 * See mcmc_gw.cpp for GW specific explanation
 *
 * Usage:
 * 	
 * 	mcmc_gw_tool /PATH/TO/PARAM/FILE
 *
 * See data/sample_config_files/mcmc_gw_tool_param_template.dat for an example parameter file
 *
 * See data/sample_init_pos.csv for an example initial position file
 *
 * NOTE: SkySearch generation method requires an initialization file with 11 dimensions. 
 * The extra parameters correspond to the injected intrinsic parameters -- exactly the format of PhenomD
 */

double **mod_priors;
double chirpmass_prior[2];
double mass1_prior[2];
double mass2_prior[2];
double spin1_prior[2];
double spin2_prior[2];
double tidal1_prior[2];
double tidal2_prior[2];
double tidal_s_prior[2];
double EA_prior[6];
double RA_bounds[2];
double sinDEC_bounds[2];
bool tidal_love=true;
bool tidal_love_error;
bool alpha_param;
bool EOS_plat_flag=false; 
double DL_prior[2];
bool tidal_love_boundary_violation(double q,double lambda_s);
double standard_log_prior_D(double *pos, mcmc_data_interface *interface,void *parameters);
double standard_log_prior_D_NRT_EA(double *pos, mcmc_data_interface *interface, void *parameters);
double EA_current_constraints(double *pos, mcmc_data_interface *interface, void *parameters);
double standard_log_prior_D_NRT(double *pos, mcmc_data_interface *interface,void *parameters);
double standard_log_prior_D_intrinsic_NRT(double *pos, mcmc_data_interface *interface,void *parameters);
double standard_log_prior_D_intrinsic_NRT_mod(double *pos, mcmc_data_interface *interface,void *parameters);
double standard_log_prior_D_NRT_mod(double *pos, mcmc_data_interface *interface,void *parameters);
double standard_log_prior_D_mod(double *pos, mcmc_data_interface *interface,void *parameters);
double standard_log_prior_Pv2(double *pos,  mcmc_data_interface *interface,void *parameters);
double standard_log_prior_D_intrinsic(double *pos, mcmc_data_interface *interface,void *parameters);
double standard_log_prior_D_intrinsic_mod(double *pos, mcmc_data_interface *interface,void *parameters);
double standard_log_prior_Pv2_intrinsic(double *pos, mcmc_data_interface *interface,void *parameters);
double standard_log_prior_Pv2_intrinsic_mod(double *pos, mcmc_data_interface *interface,void *parameters);
double standard_log_prior_skysearch(double *pos, mcmc_data_interface *interface, void *parameters);
double standard_log_prior_Pv2_mod(double *pos, mcmc_data_interface *interface,void *parameters);
double chirpmass_eta_jac(double m1,double m2);
double standard_log_prior_D_intrinsic_Jeffreys(double *pos, mcmc_data_interface *interface,void *parameters);
double chirpmass_q_jac(double chirpmass, double q);
double aligned_spin_prior(double chi);
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
		std::cout<<"Detector stream file: "<<detectors[i]<<" : "<<detector_files[i]<<std::endl;
	}
	std::string generation_method = str_dict["generation method"];
	
	double data_time_length = dbl_dict["data length"];
	std::cout<<"data length: "<<data_time_length<<std::endl;
	int data_length=131072;
	if((int)data_time_length == 32){
		data_length = 131072;
	}	
	else if((int)data_time_length == 4096){
		data_length = 16777216;
	}	
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
	int correlation_thresh =10; 
	int correlation_segs = 10;
	double correlation_convergence_thresh =0.2;
	double ac_target = 0.01;
	std::string output_file = str_dict["output data file"];
	std::string stat_file = str_dict["output stat file"];
	std::string check_file = str_dict["checkpoint file"];
	int dimension = int_dict["dimension"];
	std::cout<<"Dimension: "<<dimension<<std::endl;
	int max_chunk_size = int_dict["Max chunk size"];
	std::cout<<"Max chunk size: "<<max_chunk_size<<std::endl;
	
	std::string initial_position_file="", initial_checkpoint_file="",initial_ensemble_position_file="";
	int chain_N;
	bool continue_from_checkpoint=false;
	if(str_dict.find("initial checkpoint file") == str_dict.end()){
		//Original
		if(str_dict.find("initial position file") != str_dict.end()){
			initial_position_file = str_dict["initial position file"];
			chain_N = int_dict["chain number"];
			std::cout<<"Chain number: "<<chain_N<<std::endl;
			std::cout<<"Initial position file: "<<initial_position_file<<std::endl;
		}
		if (str_dict.find("initial ensemble position file") != str_dict.end()){
			initial_ensemble_position_file = str_dict["initial ensemble position file"];
			chain_N = int_dict["chain number"];
			std::cout<<"Chain number: "<<chain_N<<std::endl;
			std::cout<<"Initial ensemble position file: "<<initial_ensemble_position_file<<std::endl;
		}
		if ( initial_position_file =="" && initial_ensemble_position_file ==""){
			std::cout<<"ERROR -- you need an initial checkpoint file, an initial position file, or an initial ensemble file"<<std::endl;
			exit(1);
		}
	}
	else{
		//Continue	
		continue_from_checkpoint=true;
		initial_checkpoint_file = str_dict["initial checkpoint file"];
		std::cout<<"INITIAL checkpoint file: "<<initial_checkpoint_file<<std::endl;
		chain_number_from_checkpoint_file(initial_checkpoint_file,&chain_N);
		std::cout<<"Chain number: "<<chain_N<<std::endl;
		
	}

	if(dbl_dict.find("Mass1 minimum") == dbl_dict.end()){
		mass1_prior[0]=1;
		mass1_prior[1]=80;
	}
	else{
		mass1_prior[0]=dbl_dict["Mass1 minimum"];
		mass1_prior[1]=dbl_dict["Mass1 maximum"];
	}
	if(dbl_dict.find("Mass2 minimum") == dbl_dict.end()){
		mass2_prior[0]=1;
		mass2_prior[1]=80;
	}
	else{
		mass2_prior[0]=dbl_dict["Mass2 minimum"];
		mass2_prior[1]=dbl_dict["Mass2 maximum"];
	}
	if(dbl_dict.find("Spin1 minimum") == dbl_dict.end()){
		spin1_prior[0]=-1;
		spin1_prior[1]=1;
	}
	else{
		spin1_prior[0]=dbl_dict["Spin1 minimum"];
		spin1_prior[1]=dbl_dict["Spin1 maximum"];
	}
	if(dbl_dict.find("Spin2 minimum") == dbl_dict.end()){
		spin2_prior[0]=-1;
		spin2_prior[1]=1;
	}
	else{
		spin2_prior[0]=dbl_dict["Spin2 minimum"];
		spin2_prior[1]=dbl_dict["Spin2 maximum"];
	}
	if(dbl_dict.find("Luminosity distance minimum") == dbl_dict.end()){
		DL_prior[0]=10;
		DL_prior[1]=5000;
	}
	else{
		DL_prior[0]=dbl_dict["Luminosity distance minimum"];
		DL_prior[1]=dbl_dict["Luminosity distance maximum"];
	}
	std::cout<<"Range of Mass1: "<<mass1_prior[0]<<" - "<<mass1_prior[1]<<std::endl;
	std::cout<<"Range of Mass2: "<<mass2_prior[0]<<" - "<<mass2_prior[1]<<std::endl;
	std::cout<<"Range of DL: "<<DL_prior[0]<<" - "<<DL_prior[1]<<std::endl;
	if(generation_method.find("IMRPhenomPv2") != std::string::npos ){
		std::cout<<"Range of Spin1: "<<0<<" - "<<spin1_prior[1]<<std::endl;
		std::cout<<"Range of Spin2: "<<0<<" - "<<spin2_prior[1]<<std::endl;
	}
	else{
		std::cout<<"Range of Spin1: "<<spin1_prior[0]<<" - "<<spin1_prior[1]<<std::endl;
		std::cout<<"Range of Spin2: "<<spin2_prior[0]<<" - "<<spin2_prior[1]<<std::endl;

	}

	if(dbl_dict.find("tidal_s minimum") == dbl_dict.end()){
		tidal_s_prior[0]=1;
		tidal_s_prior[2]=5000;
	}
	else{
		tidal_s_prior[0]=dbl_dict["tidal_s minimum"];
		tidal_s_prior[1]=dbl_dict["tidal_s maximum"];
	}
	if(dbl_dict.find("tidal1 minimum") == dbl_dict.end()){
		tidal1_prior[0]=1;
		tidal1_prior[1]=5000;
	}
	else{
		tidal1_prior[0]=dbl_dict["tidal1 minimum"];
		tidal1_prior[1]=dbl_dict["tidal1 maximum"];
	}
	if(dbl_dict.find("tidal2 minimum") == dbl_dict.end()){
		tidal2_prior[0]=1;
		tidal2_prior[1]=5000;
	}
	else{
		tidal2_prior[0]=dbl_dict["tidal2 minimum"];
		tidal2_prior[1]=dbl_dict["tidal2 maximum"];
	}
	tidal_love = true;
	if(bool_dict.find("Tidal love relation") != bool_dict.end())
	{
		tidal_love = bool_dict["Tidal love relation"];
	}
	tidal_love_error = false;
	if(bool_dict.find("Tidal love error marginalization") != bool_dict.end())
	{
		tidal_love_error = bool_dict["Tidal love error marginalization"];
	}
	if(generation_method.find("NRT") != std::string::npos){
		std::cout<<"Range of tidal1: "<<tidal1_prior[0]<<" - "<<tidal1_prior[1]<<std::endl;
		std::cout<<"Range of tidal2: "<<tidal2_prior[0]<<" - "<<tidal2_prior[1]<<std::endl;
		std::cout<<"Range of tidal_s: "<<tidal_s_prior[0]<<" - "<<tidal_s_prior[1]<<std::endl;
		std::cout<<"Using tidal-love relations: "<<tidal_love<<std::endl;
		std::cout<<"Using tidal-love error marginalization: "<<tidal_love_error<<std::endl;
				
	}
	bool jeff_prior = false;
	if(bool_dict.find("Jefferys prior") == bool_dict.end())
	{
		jeff_prior = bool_dict["Jeffreys prior"];
	}
	//Einstein AEther parameterization
       	alpha_param = true;
	if(bool_dict.find("alpha parameterization") != bool_dict.end())
	{
	        alpha_param = bool_dict["alpha parameterization"];
	}
	//Einstein AEther priors
	if(generation_method.find("EA") != std::string::npos){
	  std::cout<<"Using alpha parameterization: "<<alpha_param<<std::endl;
	  if(alpha_param){
	    if(dbl_dict.find("EA alpha_1 minimum") ==dbl_dict.end()){
	      EA_prior[0]= -1.;
	      EA_prior[1]= 1.; 
	    }
	    else{
	      EA_prior[0]=dbl_dict["EA alpha_1 minimum"];
	      EA_prior[1]=dbl_dict["EA alpha_1 maximum"]; 
	    }
	    if(dbl_dict.find("EA alpha_2 minimum") ==dbl_dict.end()){
	      EA_prior[2]=-.1;
	      EA_prior[3]=.1; 
	    }
	    else{
	      EA_prior[2]=dbl_dict["EA alpha_2 minimum"];
	      EA_prior[3]=dbl_dict["EA alpha_2 maximum"]; 
	    }
	    if(dbl_dict.find("EA cbar_w minimum") ==dbl_dict.end()){
	      EA_prior[4]=-1.;
	      EA_prior[5]=1.; 
	    }
	    else{
	      EA_prior[4]=dbl_dict["EA cbar_w minimum"];
	      EA_prior[5]=dbl_dict["EA cbar_w maximum"]; 
	    }
	    std::cout<<"Range of EA alpha1: "<<EA_prior[0]<<" - "<<EA_prior[1]<<std::endl;
	    std::cout<<"Range of EA alpha2: "<<EA_prior[2]<<" - "<<EA_prior[3]<<std::endl;
	    std::cout<<"Range of EA cbarw: "<<EA_prior[4]<<" - "<<EA_prior[5]<<std::endl;
	  }
	  else{
	    if(dbl_dict.find("EA c_a minimum") == dbl_dict.end()){
	      EA_prior[0]=0;
	      EA_prior[1]=pow(10,-4.);
	    }
	    else{
	      EA_prior[0]=dbl_dict["EA c_a minimum"];
	      EA_prior[1]=dbl_dict["EA c_a maximum"];
	    }
	    if(dbl_dict.find("EA c_theta minimum") == dbl_dict.end()){
	      EA_prior[2]=0;
	      EA_prior[3]=pow(10,-4.);
	    }
	    else{
	      EA_prior[2]=dbl_dict["EA c_theta minimum"];
	      EA_prior[3]=dbl_dict["EA c_theta maximum"];
	    }
	    if(dbl_dict.find("EA c_w minimum") == dbl_dict.end()){
	      EA_prior[4]=-10;
	      EA_prior[5]=10;
	    }
	    else{
	      EA_prior[4]=dbl_dict["EA c_w minimum"];
	      EA_prior[5]=dbl_dict["EA c_w maximum"];
	    }
	  
	    /*
	      if(dbl_dict.find("EA c_sigma minimum") == dbl_dict.end()){
	      EA_prior[6]=-pow(10,-15.);
	      EA_prior[7]=pow(10,-15.);
	      }
	      else{
	      EA_prior[6]=dbl_dict["EA c_sigma minimum"];
	      EA_prior[7]=dbl_dict["EA c_sigma maximum"];
	      }
	    */
	    std::cout<<"Range of EA c_a: "<<EA_prior[0]<<" - "<<EA_prior[1]<<std::endl;
	    std::cout<<"Range of EA c_theta: "<<EA_prior[2]<<" - "<<EA_prior[3]<<std::endl;
	    std::cout<<"Range of EA c_w: "<<EA_prior[4]<<" - "<<EA_prior[5]<<std::endl;
	    //std::cout<<"Range of EA c_sigma: "<<EA_prior[6]<<" - "<<EA_prior[7]<<std::endl;
	  }
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
	
	double post_merger_duration = 2;
	if( dbl_dict.find("Post merger signal duration") != dbl_dict.end()){
		post_merger_duration = dbl_dict["Post merger signal duration"];
	}
	std::cout<<"Poster merger signal duration: "<<post_merger_duration<<std::endl;
	if( dbl_dict.find("RA minimum") != dbl_dict.end() && dbl_dict.find("RA maximum") != dbl_dict.end()){
		RA_bounds[0] = dbl_dict["RA minimum"];
		RA_bounds[1] = dbl_dict["RA maximum"];
		std::cout<<"RA bounds: "<<RA_bounds[0]<<" "<<RA_bounds[1]<<std::endl;
	}
	else{
		RA_bounds[0] = 0;
		RA_bounds[1] = 2*M_PI;

	}
	if( dbl_dict.find("Sin(DEC) minimum") != dbl_dict.end() && dbl_dict.find("Sin(DEC) maximum") != dbl_dict.end()){
		sinDEC_bounds[0] = dbl_dict["Sin(DEC) minimum"];
		sinDEC_bounds[1] = dbl_dict["Sin(DEC) maximum"];
		std::cout<<"Sin(DEC) bounds: "<<sinDEC_bounds[0]<<" "<<sinDEC_bounds[1]<<std::endl;
	}
	else{
		sinDEC_bounds[0] = -1;
		sinDEC_bounds[1] = 1;

	}
	

	allocate_LOSC_data(detector_files, psd_file, detector_N, psd_length, data_length, gps_time,post_merger_duration, data, psd, freqs);
	//########################################################################
	//TEST
	//for(int i = 0 ; i<detector_N; i++){
	//	for(int j =0 ; j<psd_length; j++){
	//		psd[i][j]=unpack1[j][i+1];
	//		data[i][j]=std::complex<double>(unpack2[j][2*i+1],unpack2[j][2*i+2]);
	//		freqs[i][j]= unpack1[j][0];
	//		//std::cout<<unpack1[j][1]<<std::endl;
	//		//std::cout<<unpack2[j][4]<<std::endl;
	//		//std::cout<<unpack1[1][j]<<std::endl;
	//	}	
	//}
	//freqs[0][0]=.0001;
	//freqs[1][0]=.0001;
	//########################################################################
	std::cout<<"DATA loaded"<<std::endl;

	T_mcmc_gw_tool = 1./(freqs[0][2]-freqs[0][1]);
	T_merger = T_mcmc_gw_tool - post_merger_duration;
	double df = 1./T_mcmc_gw_tool;
	//std::cout<<df<<std::endl;
	//for(int i = 0 ; i<detector_N; i++){
	//	for(int j =0 ; j<psd_length; j++){
	//		psd[i][j]*=df;
	//	}	
	//}

	std::cout<<"Total time: "<<T_mcmc_gw_tool<<std::endl;
	int *data_lengths= (int*)malloc(sizeof(int)*detector_N);
	for(int i = 0 ; i<detector_N; i++){
		data_lengths[i] =psd_length;
	}


	//double **output_test = allocate_2D_array(psd_length, 10 );
	//
	//for(int i = 0 ; i<psd_length; i++)
	//{
	//	output_test[i][0] = freqs[0][i];	
	//	output_test[i][1] = psd[0][i];	
	//	output_test[i][2] = psd[1][i];	
	//	output_test[i][3] = psd[2][i];	
	//	output_test[i][4] =real( data[0][i]);	
	//	output_test[i][5] =imag( data[0][i]);	
	//	output_test[i][6] =real( data[1][i]);	
	//	output_test[i][7] =imag( data[1][i]);	
	//	output_test[i][8] =real( data[2][i]);	
	//	output_test[i][9] =imag( data[2][i]);	
	//}
	//write_file("data/processed_data.csv",output_test,psd_length,10);
	//deallocate_2D_array(output_test,psd_length,10);
	

	//#############################################################
	//#############################################################
	
	//double **whitened = allocate_2D_array(data_lengths[0],10);
	//for(int i = 0 ; i<data_lengths[0]; i++){
	//	whitened[i][0]= freqs[0][i];
	//	whitened[i][1]=psd[0][i];
	//	whitened[i][2]=psd[1][i];
	//	whitened[i][3]=psd[2][i];
	//	whitened[i][4]=real(data[0][i])/std::pow(psd[0][i],.5);
	//	whitened[i][5]=imag(data[0][i])/std::pow(psd[0][i],.5);
	//	whitened[i][6]=real(data[1][i])/std::pow(psd[1][i],.5);
	//	whitened[i][7]=imag(data[1][i])/std::pow(psd[1][i],.5);
	//	whitened[i][8]=real(data[2][i])/std::pow(psd[2][i],.5);
	//	whitened[i][9]=imag(data[2][i])/std::pow(psd[2][i],.5);
	//}	
	//write_file("data/whitened_data.csv",whitened,data_lengths[0],10);
	//deallocate_2D_array(whitened, data_lengths[0],10);
	
	//#############################################################
	//#############################################################

	//########################################################################
	double **output;
	output = allocate_2D_array(samples, dimension );
	double chain_temps[chain_N];
	
	int Nmod = 0;
	int gNmod_phi = 0;
	int gNmod_sigma = 0;
	int gNmod_beta = 0;
	int gNmod_alpha = 0;
	int *gIMR_phii = NULL;
	int *gIMR_sigmai = NULL;
	int *gIMR_betai = NULL;
	int *gIMR_alphai = NULL;
	double *bppe = NULL;
	std::cout<<"Generation method: "<<generation_method<<std::endl;
	if(generation_method.find("ppE") != std::string::npos || check_theory_support(generation_method)){
		Nmod = int_dict["Number of modifications"];
		std::cout<<"Number of ppE modifications: "<<Nmod<<std::endl;
		std::cout<<"ppE b parmeters: "<<Nmod<<std::endl;
		bppe= new double[Nmod];
		mod_priors = new double*[Nmod];
		for(int i =0; i<Nmod ; i++){
			bppe[i] = dbl_dict["ppE b "+std::to_string(i)];
			mod_priors[i]= new double[2];
			std::cout<<i<<" : "<<bppe[i]<<std::endl;
			mod_priors[i][0] = dbl_dict["ppE beta "+std::to_string(i)+" minimum"];
			mod_priors[i][1] = dbl_dict["ppE beta "+std::to_string(i)+" maximum"];
			std::cout<<"Min"<<" : "<<mod_priors[i][0]<<std::endl;
			std::cout<<"Max"<<" : "<<mod_priors[i][1]<<std::endl;
		}
		
	}
	//else if(generation_method.find("EA") != std::string::npos){
	//	mod_priors = new double*[4];
	//	mod_priors[0]= new double[2];
	//	mod_priors[0][0] = dbl_dict["EA c_sigma minimum"];
	//	mod_priors[0][1] = dbl_dict["EA c_sigma maximum"];
	//	mod_priors[1]= new double[2];
	//	mod_priors[1][0] = dbl_dict["EA c_theta minimum"];
	//	mod_priors[1][1] = dbl_dict["EA c_theta maximum"];
	//	mod_priors[2]= new double[2];
	//	mod_priors[2][0] = dbl_dict["EA c_omega minimum"];
	//	mod_priors[2][1] = dbl_dict["EA c_omega maximum"];
	//	mod_priors[3]= new double[2];
	//	mod_priors[3][0] = dbl_dict["EA c_a minimum"];
	//	mod_priors[3][1] = dbl_dict["EA c_a maximum"];

	//}
	//if(generation_method.find("dCS") != std::string::npos
	//|| generation_method.find("EdGB") != std::string::npos){
	//	if(generation_method.find("EdGB_GHO") != std::string::npos){
	//		Nmod = 2;
	//		std::cout<<"Number of ppE modifications: "<<Nmod<<std::endl;
	//		std::cout<<"ppE b parmeters: "<<Nmod<<std::endl;
	//		bppe= new double[Nmod];
	//		mod_priors = new double*[Nmod];
	//		bppe[0] = -7;
	//		bppe[1] = -5;
	//		mod_priors[0]= new double[2];
	//		mod_priors[1]= new double[2];
	//		mod_priors[0][0] = dbl_dict["ppE beta "+std::to_string(0)+" minimum"];
	//		mod_priors[0][1] = dbl_dict["ppE beta "+std::to_string(0)+" maximum"];
	//		mod_priors[1][0] = dbl_dict["ppE beta "+std::to_string(1)+" minimum"];
	//		mod_priors[1][1] = dbl_dict["ppE beta "+std::to_string(1)+" maximum"];
	//		std::cout<<0<<" : "<<bppe[0]<<std::endl;
	//		std::cout<<"Min"<<" : "<<mod_priors[0][0]<<std::endl;
	//		std::cout<<"Max"<<" : "<<mod_priors[0][1]<<std::endl;
	//		std::cout<<0<<" : "<<bppe[1]<<std::endl;
	//		std::cout<<"Min"<<" : "<<mod_priors[1][0]<<std::endl;
	//		std::cout<<"Max"<<" : "<<mod_priors[1][1]<<std::endl;
	//	
	//	}
	//	else{
	//		Nmod = 1;
	//		std::cout<<"Number of ppE modifications: "<<Nmod<<std::endl;
	//		std::cout<<"ppE b parmeters: "<<Nmod<<std::endl;
	//		bppe= new double[Nmod];
	//		mod_priors = new double*[Nmod];
	//		if(generation_method.find("dCS") != std::string::npos){
	//			bppe[0] = -1;
	//		}
	//		if(generation_method.find("EdGB") != std::string::npos){
	//			bppe[0] = -7;
	//		}
	//		mod_priors[0]= new double[2];
	//		std::cout<<0<<" : "<<bppe[0]<<std::endl;
	//		mod_priors[0][0] = dbl_dict["ppE beta "+std::to_string(0)+" minimum"];
	//		mod_priors[0][1] = dbl_dict["ppE beta "+std::to_string(0)+" maximum"];
	//	}
	//	
	//}
	if(generation_method.find("gIMR") != std::string::npos){
		gNmod_phi = int_dict["Number of phi modifications"];
		gNmod_sigma = int_dict["Number of sigma modifications"];
		gNmod_beta = int_dict["Number of beta modifications"];
		gNmod_alpha = int_dict["Number of alpha modifications"];
		int Nmod_tot = gNmod_phi+ gNmod_sigma+gNmod_beta+gNmod_alpha;
		std::cout<<"Number of general modifications: "<<Nmod_tot<<std::endl;
		std::cout<<"delta phi i: "<<Nmod<<std::endl;
		mod_priors = new double*[Nmod_tot];
		int ct = 0;
		if(gNmod_phi != 0){
			gIMR_phii= new int[gNmod_phi];
			for(int i =0; i<gNmod_phi ; i++){
				gIMR_phii[i] = int_dict["delta phi "+std::to_string(i)+" i"];
				mod_priors[ct]= new double[2];
				std::cout<<i<<" : "<<gIMR_phii[i]<<std::endl;
				mod_priors[ct][0] = dbl_dict["delta phi "+std::to_string(i)+" minimum"];
				mod_priors[ct][1] = dbl_dict["delta phi "+std::to_string(i)+" maximum"];
				ct++;
			}
		}
		if(gNmod_sigma != 0){
			gIMR_sigmai= new int[gNmod_sigma];
			for(int i =0; i<gNmod_sigma ; i++){
				gIMR_sigmai[i] = int_dict["delta sigma "+std::to_string(i)+" i"];
				mod_priors[ct]= new double[2];
				std::cout<<i<<" : "<<gIMR_sigmai[i]<<std::endl;
				mod_priors[ct][0] = dbl_dict["delta sigma "+std::to_string(i)+" minimum"];
				mod_priors[ct][1] = dbl_dict["delta sigma "+std::to_string(i)+" maximum"];
				ct++;
			}
		}
		if(gNmod_beta != 0){
			gIMR_betai= new int[gNmod_beta];
			for(int i =0; i<gNmod_beta ; i++){
				gIMR_betai[i] = int_dict["delta beta "+std::to_string(i)+" i"];
				mod_priors[ct]= new double[2];
				std::cout<<i<<" : "<<gIMR_betai[i]<<std::endl;
				mod_priors[ct][0] = dbl_dict["delta beta "+std::to_string(i)+" minimum"];
				mod_priors[ct][1] = dbl_dict["delta beta "+std::to_string(i)+" maximum"];
				ct++;
			}
		}
		if(gNmod_alpha != 0){
			gIMR_alphai= new int[gNmod_alpha];
			for(int i =0; i<gNmod_alpha ; i++){
				gIMR_alphai[i] = int_dict["delta alpha "+std::to_string(i)+" i"];
				mod_priors[ct]= new double[2];
				std::cout<<i<<" : "<<gIMR_alphai[i]<<std::endl;
				mod_priors[ct][0] = dbl_dict["delta alpha "+std::to_string(i)+" minimum"];
				mod_priors[ct][1] = dbl_dict["delta alpha "+std::to_string(i)+" maximum"];
				ct++;
			}
		}
		
	}
	int total_mods = Nmod+gNmod_phi+gNmod_sigma+gNmod_beta+gNmod_alpha;
	if(generation_method.find("EA") != std::string::npos){total_mods+=3;}
	bool pool = true;
	if(pool){
		debugger_print(__FILE__,__LINE__,"POOLING");
	}
	else{
		debugger_print(__FILE__,__LINE__,"NOT POOLING");
	}
	bool show_progress = true;
	MCMC_modification_struct mod_struct;
	mod_struct.tidal_love = tidal_love;
	mod_struct.tidal_love_error = tidal_love_error;
	mod_struct.alpha_param = alpha_param;
	mod_struct.EOS_plat_flag = EOS_plat_flag; 
	mod_struct.ppE_Nmod = Nmod;
	mod_struct.bppe = bppe;
	mod_struct.gIMR_Nmod_phi = gNmod_phi;
	mod_struct.gIMR_phii = gIMR_phii;
	mod_struct.gIMR_Nmod_sigma = gNmod_sigma;
	mod_struct.gIMR_sigmai = gIMR_sigmai;
	mod_struct.gIMR_Nmod_beta = gNmod_beta;
	mod_struct.gIMR_betai = gIMR_betai;
	mod_struct.gIMR_Nmod_alpha = gNmod_alpha;
	mod_struct.gIMR_alphai = gIMR_alphai;
	if(bool_dict.find("NS Flag 0") != bool_dict.end()){
		mod_struct.NSflag1 = bool_dict["NS Flag 0"];
		if(mod_struct.NSflag1){
			std::cout<<"Object 0 is a NS"<<std::endl;
		}
	}
	if(bool_dict.find("NS Flag 1") != bool_dict.end()){
		mod_struct.NSflag2 = bool_dict["NS Flag 1"];
		if(mod_struct.NSflag2){
			std::cout<<"Object 1 is a NS"<<std::endl;
		}
	}

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

		//fourier_waveform(freqs[0],data_lengths[0],hplus, hcross, "IMRPhenomD",&params);
		waveform_polarizations<double> wp;
		wp.hplus = hplus;
		wp.hcross = hcross;
		fourier_waveform(freqs[0],data_lengths[0],&wp, "IMRPhenomD",&params);
	
	
		double *seeding_var = NULL;
		double **ensemble_initial_position=NULL;
		mcmc_sampler_output sampler_output(chain_N, dimension);
		SkySearch_PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW(&sampler_output,output, dimension, samples, chain_N, 
				max_thermo_chain_N, initial_position[0],seeding_var,ensemble_initial_position,chain_temps, 
				swap_freq, t0, nu,max_chunk_size,allocation_scheme, 
				standard_log_prior_skysearch,threads, pool,show_progress,detector_N, 
				data, psd,freqs, data_lengths,gps_time, detectors,Nmod, bppe,
				hplus,hcross,stat_file,output_file, "",check_file);	
		delete []  hplus;
		delete []  hcross;
		delete [] initial_position[0]; delete [] initial_position;
	}
	else{
	  std::cout<<"dimension - total_mods="<<dimension-total_mods<<std::endl;
		double(*lp)(double *param, mcmc_data_interface *interface, void *parameters);
		if(generation_method.find("IMRPhenomD") != std::string::npos && (dimension-total_mods) == 4){
			if(total_mods == 0){
				lp = &standard_log_prior_D_intrinsic;
				if(jeff_prior){
					lp = &standard_log_prior_D_intrinsic_Jeffreys;
				}
			}
			else{
				lp = &standard_log_prior_D_intrinsic_mod;
			}
		}
		else if(generation_method.find("IMRPhenomD_NRT") != std::string::npos && ( (dimension-total_mods) == 6 || (dimension-total_mods) == 5)){
			
			if(total_mods == 0){
				lp = &standard_log_prior_D_intrinsic_NRT;
			}
			else if(generation_method.find("EA") !=std::string::npos){
				lp = NULL;
			}
			else{
				lp = &standard_log_prior_D_intrinsic_NRT_mod;
			}
		}
		else if(generation_method.find("IMRPhenomD") != std::string::npos && (dimension-total_mods) == 11){
			if(total_mods == 0){
				lp = &standard_log_prior_D;
			}
			else{
				lp = &standard_log_prior_D_mod;
			}
		}
		else if(generation_method.find("IMRPhenomD_NRT") != std::string::npos && ( (dimension-total_mods) == 13 || (dimension-total_mods) == 12)){
			if(total_mods == 0){
				lp = &standard_log_prior_D_NRT;
			}
			else if(generation_method.find("EA") !=std::string::npos){
				lp = &standard_log_prior_D_NRT_EA;
			}
			else{
				lp = &standard_log_prior_D_NRT_mod;
			}
		}
		else if(generation_method.find("IMRPhenomPv2") != std::string::npos && (dimension-total_mods) == 8){
			if(total_mods == 0){
				lp = &standard_log_prior_Pv2_intrinsic;
			}
			else{
				lp = &standard_log_prior_Pv2_intrinsic_mod;
			}
		}
		else if(generation_method.find("IMRPhenomPv2") != std::string::npos && (dimension-total_mods) == 15){
			if(total_mods == 0){
				lp = &standard_log_prior_Pv2;
			}
			else{
				lp = &standard_log_prior_Pv2_mod;
			}
		}
		else{
			std::cout<<"ERROR -- wrong detector/dimension combination for this tool -- Check mcmc_gw for general support"<<std::endl;
			return 1;
		}
		//if(generation_method.find("IMRPhenomD") != std::string::npos && dimension == 11){
		//	lp = &standard_log_prior_D;
		//}
		//else if(generation_method.find("IMRPhenomPv2") != std::string::npos && dimension == 15){
		//	lp = &standard_log_prior_Pv2;
		//}
		//else if( (generation_method.find("ppE_IMRPhenomPv2") != std::string::npos  
		//	|| generation_method.find("gIMRPhenomPv2") !=std::string::npos||
		//	generation_method.find("PNSeries_ppE_IMRPhenomPv2") !=std::string::npos||
		//	generation_method.find("dCS_IMRPhenomPv2") != std::string::npos || 
		//	generation_method.find("EdGB_IMRPhenomPv2") != std::string::npos ||
		//	generation_method.find("EdGB_GHOv1_IMRPhenomPv2") != std::string::npos|| 
		//	generation_method.find("EdGB_GHOv2_IMRPhenomPv2") != std::string::npos|| 
		//	generation_method.find("EdGB_GHOv3_IMRPhenomPv2") != std::string::npos)
		//	&& dimension >= 16 ){
		//	lp = &standard_log_prior_Pv2_mod;
		//}
		//else if( (generation_method.find("ppE_IMRPhenomD") != std::string::npos 
		//	|| generation_method.find("gIMRPhenomD") != std::string::npos
		//	|| generation_method.find("PNSeries_ppE_IMRPhenomD") != std::string::npos
		//	|| generation_method.find("dCS_IMRPhenomD") != std::string::npos 
		//	|| generation_method.find("EdGB_IMRPhenomD") != std::string::npos 
		//	|| generation_method.find("EdGB_GHOv1_IMRPhenomD") != std::string::npos
		//	|| generation_method.find("EdGB_GHOv2_IMRPhenomD") != std::string::npos
		//	|| generation_method.find("EdGB_GHOv3_IMRPhenomD") != std::string::npos)
		//	&& dimension >= 11){
		//	lp = &standard_log_prior_D_mod;
		//}
		//else if(generation_method.find("IMRPhenomPv2") != std::string::npos && dimension == 8){
		//	lp = &standard_log_prior_Pv2_intrinsic;
		//}
		//else if(generation_method.find("IMRPhenomD") != std::string::npos && dimension == 4){
		//	lp = &standard_log_prior_D_intrinsic;
		//}
		//else if((generation_method.find("ppE_IMRPhenomD") != std::string::npos 
		//	|| generation_method.find("gIMRPhenomD") != std::string::npos) 
		//	&& dimension >= 4){
		//	lp = &standard_log_prior_D_intrinsic_mod;
		//}
		//else if( (generation_method.find("ppE_IMRPhenomPv2") != std::string::npos  
		//	|| generation_method.find("gIMRPhenomPv2") !=std::string::npos)
		//	&& dimension >= 9 ){
		//	lp = &standard_log_prior_Pv2_intrinsic_mod;
		//}
		//else if((generation_method.find("IMRPhenomD_NRT") != std::string::npos) && (dimension == 6)){
		//	lp = &standard_log_prior_D_intrinsic_NRT;

		//}
		//else if((generation_method.find("IMRPhenomD_NRT") != std::string::npos) && (dimension == 13)){
		//	lp = &standard_log_prior_D_NRT;

		//}
		//else{
		//	std::cout<<"ERROR -- wrong detector/dimension combination for this tool -- Check mcmc_gw for general support"<<std::endl;
		//	return 1;
		//}

		if(continue_from_checkpoint){
			mcmc_sampler_output sampler_output(chain_N, dimension);
			continue_PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW(initial_checkpoint_file,&sampler_output,output, samples,  
					max_thermo_chain_N, chain_temps, 
					swap_freq, t0, nu,max_chunk_size,allocation_scheme, 
					lp,threads, pool,show_progress,detector_N, 
					data, psd,freqs, data_lengths,gps_time, detectors,&mod_struct,
					generation_method,stat_file,output_file, "",check_file);	

		}
		else{
			double **initial_position=NULL;
			initial_position = new double*[1];
			initial_position[0] = new double[dimension];
			double *seeding_var = NULL;
			double **ensemble_initial_position = NULL;
			if(initial_position_file != ""){
				read_file(initial_position_file, initial_position,1,dimension);
			}
			if(initial_ensemble_position_file != ""){
				ensemble_initial_position = new double*[chain_N];
				for(int i = 0 ; i<chain_N; i++){
					ensemble_initial_position[i] = new double[dimension];
				}
				read_file(initial_ensemble_position_file, ensemble_initial_position,chain_N,dimension);
				//for(int i = 0 ; i<chain_N; i++){
				//	for(int j = 0 ; j<dimension; j++){
				//		std::cout<<ensemble_initial_position[i][j]<<" ";
				//	}	
				//	std::cout<<std::endl;
				//}
			}
			std::cout<<"Running uncorrelated sampler "<<std::endl;
			mcmc_sampler_output sampler_output(chain_N, dimension);
			PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW(&sampler_output,output, dimension, samples, chain_N, 
					max_thermo_chain_N, initial_position[0],seeding_var,ensemble_initial_position,chain_temps, 
					swap_freq, t0, nu,max_chunk_size,allocation_scheme, 
					lp,threads, pool,show_progress,detector_N, 
					//data, psd,freqs, data_lengths,gps_time, detectors,Nmod, bppe,
					data, psd,freqs, data_lengths,gps_time, detectors,&mod_struct,
					generation_method,stat_file,output_file, "",check_file);	
			sampler_output.create_data_dump(true, false, output_file);
			delete [] initial_position[0]; delete [] initial_position;
			if(initial_ensemble_position_file != ""){
				for(int i = 0 ; i<chain_N; i++){
					delete [] ensemble_initial_position[i]; 
				}
				delete [] ensemble_initial_position;
			}
		}

	}



	deallocate_2D_array(output, samples, dimension);
	free(data_lengths);
	deallocate_2D_array(psd,detector_N, psd_length);
	deallocate_2D_array(freqs,detector_N, psd_length);
	for(int i =0; i<detector_N; i++)
		free(data[i]);
	free(data);

	if(gNmod_phi != 0){
		delete [] gIMR_phii;
	}
	if(gNmod_sigma != 0){
		delete [] gIMR_sigmai;
	}
	if(gNmod_beta != 0){
		delete [] gIMR_betai;
	}
	if(gNmod_alpha != 0){
		delete [] gIMR_alphai;
	}
	if(generation_method.find("ppE") != std::string::npos){
		delete [] bppe;	
		for(int i = 0 ; i<Nmod ; i++){
			delete [] mod_priors[i] ;
		}
		delete [] mod_priors;
	}
	else if(generation_method.find("EA") !=std::string::npos){
		for(int i = 0 ; i<4 ; i++){
			delete [] mod_priors[i] ;
		}
		delete [] mod_priors;

	}
	
	return 0;

}

double standard_log_prior_D_mod(double *pos, mcmc_data_interface *interface,void *parameters)
{
	int dim =  interface->max_dim;
	double a = -std::numeric_limits<double>::infinity();
	//###########
	double chirp = exp(pos[7]);
	double eta = pos[8];
	if (eta<.0 || eta>.25){return a;}//eta
	double m1 = calculate_mass1(chirp,eta );
	double m2 = calculate_mass2(chirp,eta );
	if(m1<mass1_prior[0] || m1>mass1_prior[1]){return a;}
	if(m2<mass2_prior[0] || m2>mass2_prior[1]){return a;}
	//###########
	if ((pos[0])<RA_bounds[0] || (pos[0])>RA_bounds[1]){ return a;}//RA
	if ((pos[1])<sinDEC_bounds[0] || (pos[1])>sinDEC_bounds[1]){return a;}//sinDEC
	if ((pos[2])<0 || (pos[2])>M_PI){return a;}//PSI
	if ((pos[3])<-1 || (pos[3])>1){return a;}//cos \iota
	if ((pos[4])<0 || (pos[4])>2*M_PI){return a;}//phiRef
	if( pos[5] < (T_merger - .1) || pos[5] > (T_merger + .1)) { return a; }
	if (std::exp(pos[6])<DL_prior[0] || std::exp(pos[6])>DL_prior[1]){return a;}//DL
	if ((pos[9])<spin1_prior[0] || (pos[9])>spin1_prior[1]){return a;}//chi1 
	if ((pos[10])<spin2_prior[0] || (pos[10])>spin2_prior[1]){return a;}//chi2
	for(int i = 0 ; i<dim - 11; i++){
		if( pos[11+i] <mod_priors[i][0] || pos[11+i] >mod_priors[i][1]){return a;}
	}
	return log(chirpmass_eta_jac(chirp,eta))+3*pos[6] ;

}
double EA_current_constraints(double *pos, mcmc_data_interface *interface, void *parameters)
{
  //std::cout<<"Checking EA_constraints"<<std::endl; 

  int dim =  interface->max_dim;
  double a = -std::numeric_limits<double>::infinity();

  source_parameters<double> sp;
  double lnChirpmass = pos[7];//ln(M_sol)
  double eta = pos[8];
  
  sp.mass1 = calculate_mass1(std::exp(lnChirpmass)*MSOL_SEC,eta);
  sp.mass2 = calculate_mass2(std::exp(lnChirpmass)*MSOL_SEC,eta);
  sp.M = sp.mass1 + sp.mass2;
  sp.alpha_param = alpha_param;
  
  //Setting tidal1 and tidal2 from tidal_s
  if(tidal_love)
    {
      IMRPhenomD_NRT<double> modelNRT;
      modelNRT.binary_love_relation(pos[11], tidal_love_error, &sp);
    }
  sp.csigma_EA = 0;
  
  if(tidal_love){
    if(alpha_param){
      sp.alpha1_EA = pos[12]; //alpha1
      sp.alpha2_EA = pos[13]; //alpha2
      sp.cbarw_EA = pos[14]; //cbarw
    }
    else{
      sp.ca_EA = pos[12]; //ca
      sp.ctheta_EA = pos[13]; //ctheta
      sp.cw_EA = pos[14]; //cw
    }
  }
  else{
    if(alpha_param){
      sp.alpha1_EA = pos[13]; //alpha1
      sp.alpha2_EA = pos[14]; //alpha2
      sp.cbarw_EA = pos[15]; //cbarw
    }
    else{
      sp.ca_EA = pos[13]; //ca
      sp.ctheta_EA = pos[14]; //ctheta
      sp.cw_EA = pos[15]; //cw
    }
  }
  sp.EA_nan_error_message = false;
  EA_IMRPhenomD_NRT<double> model;
  model.pre_calculate_EA_factors(&sp);

  if(sp.ca_EA < 0 || sp.ca_EA > 2.){return a;}
  /* Throws out points with ca < 0 or ca > 2 because these violate 
   * the positive energy condition for the spin-0 mode (scalar mode).   
   * See equation 40 of arXiv:gr-qc/0507059v3.
   */
  if(sp.cw_EA < (-sp.csigma_EA/(1. - sp.csigma_EA))){return a;}
  /* Throws out points with cw < -csigma/(1 - csigma) because these violate 
   * the positive energy condition for the spin-1 mode (vector mode).
   * Note that the positive energy condition for the spin-2 modes is always 
   * satisfied. See equation 40 of arXiv:gr-qc/0507059v3.
   */
  //std::cout<<"EA constraints test 1"<<std::endl; 
  if(sp.cTsq_EA < 0 || sp.cVsq_EA < 0 || sp.cSsq_EA < 0){return a;}
  /* Throws out points with speeds not greater than or equal to zero (these 
   * would produce gradient instabilities or ghosts)
   * arXiv:gr-qc/0402005 and arXiv:1108.1835
   */
  //std::cout<<"EA constraints test 2"<<std::endl; 
  
  if(isnan(sp.kappa3_EA))
    {
      //std::cout<<"kappa3:"<<sp.kappa3_EA<<std::endl; 
      return a; 
      }
  //std::cout<<"EA constraints test 3"<<std::endl; 

  if(sp.ctheta_EA < 0 || (sp.ctheta_EA + (8.*sp.ca_EA/7)) > (2./7.)){return a;}
  /*
   * Throws out points that violate Big Bang Nucleosynthesis constraints. arXiv:hep-th/0407149v3 
   */

  if(sp.s1_EA >=1 || sp.s2_EA >=1){return a;}
  
  //if(fabs(sp.alpha1_EA) > pow(10, -4.) || fabs(sp.alpha2_EA) > 4.*pow(10, -7.)){return a;}
  /* Throws out points that do not obey observational solar system constraints on 
   * alpha1 and alpha2
   * arXiv:1403.7377 and arXiv:gr-qc/0509114
   */
  //std::cout<<"EA constraints test 5"<<std::endl; 
 
  bool Cherenkov = true;
  if(Cherenkov)
  {
  	bool violate = false;
  	if(sp.cV_EA < 1)
  	  {
      		if(abs(sp.c13_EA * sp.c13_EA *(sp.c1_EA * sp.c1_EA + 2*sp.c1_EA*sp.c3_EA + sp.c3_EA * sp.c3_EA - 2*sp.c4_EA)/(2*sp.c1_EA*sp.c1_EA)) >= 7*pow(10, -32.))
		//enforcing constraint from Eq. 4.7 of arXiv:hep-ph/0505211
		{
	  	  violate = true;
		}
	      
	    }
	  if(violate){return a;}
	  //std::cout<<"EA constraints test 6"<<std::endl; 
  
	  if(sp.cS_EA < 1)
	    {
	      if(fabs((sp.c2_EA + sp.c3_EA - sp.c4_EA)/sp.c1_EA) > pow(10, -22.))
		{
		  if((sp.c3_EA - sp.c4_EA)*(sp.c3_EA - sp.c4_EA)/fabs(sp.c14_EA) >= pow(10, -30.)){return a;}
		  //enforcing constraint from Eq.4.15 of arXiv:hep-ph/0505211
		}
	      //std::cout<<"EA constraints test 7"<<std::endl; 
	      if(fabs((sp.c4_EA - sp.c2_EA - sp.c3_EA)/sp.c1_EA) >= 3*pow(10,-19.)){return a;}
	      //enforcing constraint from Eq.5.14 of arXiv:hep-ph/0505211
	
	    }
  } 

  //std::cout<<"Made it to end of EA constraints"<<std::endl;
  /*
  //Enforcing gaussian prior on alpha1 from binary pulsar and triple systems. 
  //arXiv:2104.04596
  double sigma = 1.021*pow(10, -5.); 
  double mu = -0.563*pow(10, -5.);
  double prob;
  
  prob = exp(-(1./2.)*((sp.alpha1_EA - mu)*(sp.alpha1_EA - mu))/(sigma*sigma));
  */

  sp.EA_nan_error_message = true;
  model.EA_check_nan(&sp);
  
  //return log(prob); //Use if enforcing gaussian on alpha1
  return 0; 
}
double standard_log_prior_D_NRT_EA(double *pos, mcmc_data_interface *interface, void *parameters)
{
  int dim =  interface->max_dim;
  double a = -std::numeric_limits<double>::infinity();
  if(tidal_love){
    if(pos[12]<EA_prior[0] || pos[12]>EA_prior[1]){return a;} //ca or alpha1
    if(pos[13]<EA_prior[2] || pos[13]>EA_prior[3]){return a;} //ctheta or alpha2
    if(pos[14]<EA_prior[4] || pos[14]>EA_prior[5]){return a;} //cw or cbarw
  }
  else{
    if(pos[13]<EA_prior[0] || pos[13]>EA_prior[1]){return a;} //ca or alpha1
    if(pos[14]<EA_prior[2] || pos[14]>EA_prior[3]){return a;} //ctheta or alpha2
    if(pos[15]<EA_prior[4] || pos[15]>EA_prior[5]){return a;} //cw or cbarw
  }
  
  double NS = standard_log_prior_D_NRT(pos,interface, parameters);
  if(NS == a){return a;}
  double EA_constraints =  EA_current_constraints(pos, interface, parameters);
  if (EA_constraints ==a){return a;}
  
  return EA_constraints + NS;
}
double standard_log_prior_D_NRT_mod(double *pos, mcmc_data_interface *interface,void *parameters)
{
	int dim =  interface->max_dim;
	double a = -std::numeric_limits<double>::infinity();
	for(int i = 0 ; i<dim - 13; i++){
		if( pos[13+i] <mod_priors[i][0] || pos[13+i] >mod_priors[i][1]){return a;}
	}
	return standard_log_prior_D_NRT(pos,interface, parameters);

}

bool tidal_love_boundary_violation(double q,double lambda_s)
{
	//Relation from 1903.03909v7 fit as rough threshhold for the 
	//validity of the tidal_s-tidal_a-q relation
	//Fit in log(lambda_s) - q space
	if(  q< 1.2321 - .124616*log(lambda_s)){return true;}
	
	return false;
}
double standard_log_prior_D_NRT(double *pos, mcmc_data_interface *interface,void *parameters)
{
	int dim =  interface->max_dim;
	double a = -std::numeric_limits<double>::infinity();
	double chirp = exp(pos[7]);
	double m1 = calculate_mass1(chirp,pos[8]);
	double m2 = calculate_mass2(chirp,pos[8]);
	double q = m2/m1;//<1
	if(tidal_love){
		if(pos[11]<tidal_s_prior[0] || pos[11]>tidal_s_prior[1]){return a;}
		if(tidal_love_boundary_violation(q,pos[11])){return a;}
		
	}
	else{
		if(pos[11]<tidal1_prior[0] || pos[11]>tidal1_prior[1]){return a;}
		if(pos[12]<tidal2_prior[0] || pos[12]>tidal2_prior[1]){return a;}
	}
	return standard_log_prior_D(pos,interface, parameters);

}
double standard_log_prior_D(double *pos, mcmc_data_interface *interface,void *parameters)
{
	int dim =  interface->max_dim;
	double a = -std::numeric_limits<double>::infinity();
	//###########
	double chirp = exp(pos[7]);
	double eta = pos[8];
	if (eta<.0 || eta>.25){return a;}//eta
	double m1 = calculate_mass1(chirp,eta );
	double m2 = calculate_mass2(chirp,eta );
	if(m1<mass1_prior[0] || m1>mass1_prior[1]){return a;}
	if(m2<mass2_prior[0] || m2>mass2_prior[1]){return a;}
	//###########
	if ((pos[0])<RA_bounds[0] || (pos[0])>RA_bounds[1]){ return a;}//RA
	if ((pos[1])<sinDEC_bounds[0] || (pos[1])>sinDEC_bounds[1]){return a;}//sinDEC
	if ((pos[2])<0 || (pos[2])>M_PI){return a;}//PSI
	if ((pos[3])<-1 || (pos[3])>1){return a;}//cos \iota
	if ((pos[4])<0 || (pos[4])>2*M_PI){return a;}//phiRef
	if( pos[5] < (T_merger - .1) || pos[5] > (T_merger + .1)) { return a; }
	if (std::exp(pos[6])<DL_prior[0] || std::exp(pos[6])>DL_prior[1]){return a;}//DL
	if ((pos[9])<spin1_prior[0] || (pos[9])>spin1_prior[1]){return a;}//chi1 
	if ((pos[10])<spin2_prior[0] || (pos[10])>spin2_prior[1]){return a;}//chi2
	//return log(chirpmass_eta_jac(chirp,eta))+3*pos[6] ;
	return log(aligned_spin_prior(pos[9]))+log(aligned_spin_prior(pos[10])) + log(chirpmass_eta_jac(chirp,eta))+3*pos[6] ;
	//return -.5*pow_int(pos[0]-5,2)/.01/.01-.5*pow_int(pos[1]+.9,2)/.01/.01+log(aligned_spin_prior(pos[9]))+log(aligned_spin_prior(pos[10])) + log(chirpmass_eta_jac(chirp,eta))+3*pos[6] ;

}
double standard_log_prior_Pv2(double *pos, mcmc_data_interface *interface,void *parameters)
{
	int dim = interface->max_dim;
	double a = -std::numeric_limits<double>::infinity();
	//####################
	double chirp = std::exp(pos[7]);
	double eta = pos[8];
	if (eta<.0 || eta>.25){return a;}//eta
	double m1 = calculate_mass1(chirp,eta );
	double m2 = calculate_mass2(chirp,eta );
	if(m1<mass1_prior[0] || m1>mass1_prior[1]){return a;}
	if(m2<mass2_prior[0] || m2>mass2_prior[1]){return a;}
	//####################

	if ((pos[0])<RA_bounds[0] || (pos[0])>RA_bounds[1]){ return a;}//RA
	if ((pos[1])<sinDEC_bounds[0] || (pos[1])>sinDEC_bounds[1]){return a;}//sinDEC
	if ((pos[2])<0 || (pos[2])>M_PI){return a;}//PSI
	if ((pos[3])<-1 || (pos[3])>1){return a;}//cos \iota
	if ((pos[4])<0 || (pos[4])>2*M_PI){return a;}//PhiRef
	if( pos[5] < (T_merger - .1) || pos[5] > (T_merger + .1)) { return a; }
	if (std::exp(pos[6])<DL_prior[0] || std::exp(pos[6])>DL_prior[1]){return a;}//DL
	if ((pos[9])<0 || (pos[9])>spin1_prior[1]){return a;}//a1 
	if ((pos[10])<0 || (pos[10])>spin2_prior[1]){return a;}//a2
	if ((pos[11])<-1 || (pos[11])>1){return a;}//theta1
	if ((pos[12])<-1 || (pos[12])>1){return a;}//theta2
	if ((pos[13])<0 || (pos[13])>2*M_PI){return a;}//phip
	if ((pos[14])<0 || (pos[14])>2*M_PI){return a;}//phip
	else {return log(chirpmass_eta_jac(chirp,eta))+3*pos[6];}
}
double standard_log_prior_Pv2_mod(double *pos, mcmc_data_interface *interface,void *parameters)
{
	int dim = interface->max_dim;
	double a = -std::numeric_limits<double>::infinity();
	//####################
	double chirp = std::exp(pos[7]);
	double eta = pos[8];
	if (eta<.0 || eta>.25){return a;}//eta
	double m1 = calculate_mass1(chirp,eta );
	double m2 = calculate_mass2(chirp,eta );
	if(m1<mass1_prior[0] || m1>mass1_prior[1]){return a;}
	if(m2<mass2_prior[0] || m2>mass2_prior[1]){return a;}
	//####################

	if ((pos[0])<RA_bounds[0] || (pos[0])>RA_bounds[1]){ return a;}//RA
	if ((pos[1])<sinDEC_bounds[0] || (pos[1])>sinDEC_bounds[1]){return a;}//sinDEC
	if ((pos[2])<0 || (pos[2])>M_PI){return a;}//PSI
	if ((pos[3])<-1 || (pos[3])>1){return a;}//cos \iota
	if ((pos[4])<0 || (pos[4])>2*M_PI){return a;}//PhiRef
	if( pos[5] < (T_merger - .1) || pos[5] > (T_merger + .1)) { return a; }
	if (std::exp(pos[6])<DL_prior[0] || std::exp(pos[6])>DL_prior[1]){return a;}//DL
	if ((pos[9])<0 || (pos[9])>spin1_prior[1]){return a;}//a1 
	if ((pos[10])<0 || (pos[10])>spin2_prior[1]){return a;}//a2
	if ((pos[11])<-1 || (pos[11])>1){return a;}//theta1
	if ((pos[12])<-1 || (pos[12])>1){return a;}//theta2
	if ((pos[13])<0 || (pos[13])>2*M_PI){return a;}//phip
	if ((pos[14])<0 || (pos[14])>2*M_PI){return a;}//phip
	for(int i = 0 ; i<dim - 15; i++){
		if( pos[15+i] <mod_priors[i][0] || pos[15+i] >mod_priors[i][1]){return a;}
	}
	return log(chirpmass_eta_jac(chirp,eta))+3*pos[6];
}
double standard_log_prior_D_intrinsic_NRT(double *pos, mcmc_data_interface *interface,void *parameters)
{
	int dim =  interface->max_dim;
	double a = -std::numeric_limits<double>::infinity();

	double chirp = exp(pos[0]);
	double m1 = calculate_mass1(chirp,pos[1]);
	double m2 = calculate_mass2(chirp,pos[1]);
	double q = m2/m1;//<1
	if(tidal_love){
		if(pos[4]<tidal_s_prior[0] || pos[4]>tidal_s_prior[1]){return a;}
		if(tidal_love_boundary_violation(q,pos[4])){return a;}
	}
	else{
		if(pos[4]<tidal1_prior[0] || pos[4]>tidal1_prior[1]){return a;}
		if(pos[5]<tidal2_prior[0] || pos[5]>tidal2_prior[1]){return a;}
	}
	return standard_log_prior_D_intrinsic(pos,interface, parameters);

}
double standard_log_prior_D_intrinsic_NRT_mod(double *pos, mcmc_data_interface *interface,void *parameters)
{
	int dim =  interface->max_dim;
	double a = -std::numeric_limits<double>::infinity();
	for(int i = 0 ; i<dim - 6; i++){
		if( pos[6+i] <mod_priors[i][0] || pos[6+i] >mod_priors[i][1]){return a;}
	}
	return standard_log_prior_D_intrinsic_NRT(pos,interface,parameters) ;

}
double standard_log_prior_D_intrinsic_Jeffreys(double *pos, mcmc_data_interface *interface,void *parameters)
{
	double a = -std::numeric_limits<double>::infinity();
	int dim = interface->max_dim;
	double lp = standard_log_prior_D_intrinsic(pos,interface,parameters);
	lp = 0;//Not using the numeric value of standard prior, just the range

	double **fisher	= new double*[dim];
	double *fisher1D = new double[dim*dim];
	for(int i = 0 ; i<dim; i++){
		fisher[i] = new double[dim];
	}
	MCMC_fisher_wrapper(pos, fisher, interface,parameters);
	for(int i = 0 ; i<dim; i++){
		for(int j = 0 ; j<dim ; j++){
			fisher1D[i*dim + j]  = fisher[i][j];
		}
	}
	Eigen::Map<Eigen::MatrixXd> m(fisher1D,dim,dim);
		
	double det = m.determinant();
	for(int i = 0 ; i<dim; i++){
		delete [] fisher[i];	
	}
	delete [] fisher;
	delete [] fisher1D;
	if(std::isnan(sqrt(det)) || std::isinf(sqrt(det))){return a;}
	//debugger_print(__FILE__,__LINE__,sqrt(det));
	//std::cout<<sqrt(det)<<std::endl;
	return sqrt(det);
}
double standard_log_prior_D_intrinsic(double *pos, mcmc_data_interface *interface,void *parameters)
{
	int dim = interface->max_dim;
	double a = -std::numeric_limits<double>::infinity();
	//###########
	double chirp = exp(pos[0]);
	double eta = pos[1];
	if (eta<.0 || eta>.25){return a;}//eta
	double m1 = calculate_mass1(chirp,eta );
	double m2 = calculate_mass2(chirp,eta );
	if(m1<mass1_prior[0] || m1>mass1_prior[1]){return a;}
	if(m2<mass2_prior[0] || m2>mass2_prior[1]){return a;}
	//###########
	if ((pos[2])<spin1_prior[0] || (pos[2])>spin1_prior[1]){return a;}//chi1 
	if ((pos[3])<spin2_prior[0] || (pos[3])>spin2_prior[1]){return a;}//chi2
	//else {return log(chirpmass_eta_jac(chirp,eta)) ;}
	else {return log(aligned_spin_prior(pos[2]))+log(aligned_spin_prior(pos[3])) + log(chirpmass_eta_jac(chirp,eta));} 


}
double standard_log_prior_D_intrinsic_mod(double *pos, mcmc_data_interface *interface,void *parameters)
{
	int dim = interface->max_dim;
	double a = -std::numeric_limits<double>::infinity();
	//###########
	double chirp = exp(pos[0]);
	double eta = pos[1];
	if (eta<.0 || eta>.25){return a;}//eta
	double m1 = calculate_mass1(chirp,eta );
	double m2 = calculate_mass2(chirp,eta );
	if(m1<mass1_prior[0] || m1>mass1_prior[1]){return a;}
	if(m2<mass2_prior[0] || m2>mass2_prior[1]){return a;}
	//###########
	if ((pos[2])<spin1_prior[0] || (pos[2])>spin1_prior[1]){return a;}//chi1 
	if ((pos[3])<spin2_prior[0] || (pos[3])>spin2_prior[1]){return a;}//chi2
	for(int i = 0 ; i<dim - 4; i++){
		if( pos[4+i] <mod_priors[i][0] || pos[4+i] >mod_priors[i][1]){return a;}
	}
	return log(chirpmass_eta_jac(chirp,eta)) ;


}
double standard_log_prior_Pv2_intrinsic(double *pos, mcmc_data_interface *interface ,void *parameters)
{
	int dim = interface->max_dim;
	double a = -std::numeric_limits<double>::infinity();
	//###########
	double chirp = exp(pos[0]);
	double eta = pos[1];
	if (eta<.0 || eta>.25){return a;}//eta
	double m1 = calculate_mass1(chirp,eta );
	double m2 = calculate_mass2(chirp,eta );
	if(m1<mass1_prior[0] || m1>mass1_prior[1]){return a;}
	if(m2<mass2_prior[0] || m2>mass2_prior[1]){return a;}
	//###########
	if ((pos[2])<0 || (pos[2])>spin1_prior[1]){return a;}//chi1 
	if ((pos[3])<0 || (pos[3])>spin2_prior[1]){return a;}//chi2
	if ((pos[4])<-1 || (pos[4])>1){return a;}//chi1 
	if ((pos[5])<-1 || (pos[5])>1){return a;}//chi2
	if ((pos[6])<0 || (pos[6])>2*M_PI){return a;}//chi2
	if ((pos[7])<0 || (pos[7])>2*M_PI){return a;}//chi2
	return log(chirpmass_eta_jac(chirp,eta)) ;

}
double standard_log_prior_Pv2_intrinsic_mod(double *pos, mcmc_data_interface *interface,void *parameters)
{
	int dim = interface->max_dim;
	double a = -std::numeric_limits<double>::infinity();
	//###########
	double chirp = exp(pos[0]);
	double eta = pos[1];
	if (eta<.0 || eta>.25){return a;}//eta
	double m1 = calculate_mass1(chirp,eta );
	double m2 = calculate_mass2(chirp,eta );
	if(m1<mass1_prior[0] || m1>mass1_prior[1]){return a;}
	if(m2<mass2_prior[0] || m2>mass2_prior[1]){return a;}
	//###########
	if ((pos[2])<0 || (pos[2])>spin1_prior[1]){return a;}//chi1 
	if ((pos[3])<0 || (pos[3])>spin2_prior[1]){return a;}//chi2
	if ((pos[4])<-1 || (pos[4])>1){return a;}//chi1 
	if ((pos[5])<-1 || (pos[5])>1){return a;}//chi2
	if ((pos[6])<0 || (pos[6])>2*M_PI){return a;}//chi2
	if ((pos[7])<0 || (pos[7])>2*M_PI){return a;}//chi2
	for(int i = 0 ; i<dim - 8; i++){
		if( pos[8+i] <mod_priors[i][0] || pos[8+i] >mod_priors[i][1]){return a;}
	}
	return log(chirpmass_eta_jac(chirp,eta)) ;

}

double standard_log_prior_skysearch(double *pos, mcmc_data_interface *interface, void *parameters){

	double a = -std::numeric_limits<double>::infinity();
	if ((pos[0])<RA_bounds[0] || (pos[0])>RA_bounds[1]){ return a;}//RA
	if ((pos[1])<sinDEC_bounds[0] || (pos[1])>sinDEC_bounds[1]){return a;}//sinDEC
	if ((pos[2])<0 || (pos[2])>M_PI){return a;}//PSI
	if ((pos[3])<-1 || (pos[3])>1){return a;}//cos \iota
	if ((pos[4])<0 || (pos[4])>2*M_PI){return a;}//PhiRef
	if( pos[5] < (T_merger - .1) || pos[5] > (T_merger + .1)) { return a; }
	if (std::exp(pos[6])<10 || std::exp(pos[6])>1000){return a;}//DL
	else { return 3*pos[6];}
}

//Uniform in m1 and m2, transformed to lnM and eta
double chirpmass_eta_jac(double chirpmass, double eta){
	return chirpmass*chirpmass/(sqrt(1. - 4.*eta)*pow(eta,1.2));
}

//Uniform in m1 and m2, transformed to lnM and q
double chirpmass_q_jac(double chirpmass, double q){
	return chirpmass*chirpmass/(pow(q/pow_int(q+1,2),1./5.) * q);
}

//For use with spin aligned system 
//-- mimics uniform magnitude and uniform cos \tilt, but mapped to chi
//-- Fit in mathematica 
double aligned_spin_prior(double chi){
	double a=0.0039132 , b= 3.95381;
	return a * exp(-b * abs(chi));	
}
