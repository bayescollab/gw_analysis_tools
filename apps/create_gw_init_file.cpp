#include <iostream>
#include "util.h"
#include "fisher.h"
#include "io_util.h"
#include "mcmc_gw.h"
#include "mcmc_sampler_internals.h"
#include "ppE_utilities.h"


int main(int argc, char *argv[])
{
	std::cout<<"Creating an initial ensemble file from mcmc data file"<<std::endl;
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


	int detector_N = int_dict["detector number"];
	std::cout<<"DETECTOR NUMBER: "<<detector_N<<std::endl;
	
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
	int threads = int_dict["thread number"];
	std::cout<<"threads: "<<threads<<std::endl;
	std::string output_file = str_dict["output data file"];
	std::string stat_file = str_dict["output stat file"];
	std::string check_file = str_dict["checkpoint file"];
	int dimension = int_dict["dimension"];
	std::cout<<"Dimension: "<<dimension<<std::endl;
	
	std::string initial_position_file="", initial_checkpoint_file="",initial_ensemble_position_file="";
	int chain_N;
	if(str_dict.find("initial position file") != str_dict.end()){
		initial_position_file = str_dict["initial position file"];
		chain_N = int_dict["chain number"];
		std::cout<<"Chain number: "<<chain_N<<std::endl;
		std::cout<<"Initial position file: "<<initial_position_file<<std::endl;
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

	allocate_LOSC_data(detector_files, psd_file, detector_N, psd_length, data_length, gps_time,post_merger_duration, data, psd, freqs);
	std::cout<<"DATA loaded"<<std::endl;
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
	if(has_substring(generation_method, "ppE") || check_theory_support(generation_method)){
		Nmod = int_dict["Number of modifications"];
		std::cout<<"Number of ppE modifications: "<<Nmod<<std::endl;
		std::cout<<"ppE b parmeters: "<<Nmod<<std::endl;
		bppe= new double[Nmod];
		for(int i =0; i<Nmod ; i++){
			bppe[i] = dbl_dict["ppE b "+std::to_string(i)];
			std::cout<<i<<" : "<<bppe[i]<<std::endl;
		}
		
	}
	if(has_substring(generation_method, "gIMR")){
		gNmod_phi = int_dict["Number of phi modifications"];
		gNmod_sigma = int_dict["Number of sigma modifications"];
		gNmod_beta = int_dict["Number of beta modifications"];
		gNmod_alpha = int_dict["Number of alpha modifications"];
		int Nmod_tot = gNmod_phi+ gNmod_sigma+gNmod_beta+gNmod_alpha;
		std::cout<<"Number of general modifications: "<<Nmod_tot<<std::endl;
		std::cout<<"delta phi i: "<<Nmod<<std::endl;
		int ct = 0;
		if(gNmod_phi != 0){
			gIMR_phii= new int[gNmod_phi];
			for(int i =0; i<gNmod_phi ; i++){
				gIMR_phii[i] = int_dict["delta phi "+std::to_string(i)+" i"];
				std::cout<<i<<" : "<<gIMR_phii[i]<<std::endl;
				ct++;
			}
		}
		if(gNmod_sigma != 0){
			gIMR_sigmai= new int[gNmod_sigma];
			for(int i =0; i<gNmod_sigma ; i++){
				gIMR_sigmai[i] = int_dict["delta sigma "+std::to_string(i)+" i"];
				std::cout<<i<<" : "<<gIMR_sigmai[i]<<std::endl;
				ct++;
			}
		}
		if(gNmod_beta != 0){
			gIMR_betai= new int[gNmod_beta];
			for(int i =0; i<gNmod_beta ; i++){
				gIMR_betai[i] = int_dict["delta beta "+std::to_string(i)+" i"];
				std::cout<<i<<" : "<<gIMR_betai[i]<<std::endl;
				ct++;
			}
		}
		if(gNmod_alpha != 0){
			gIMR_alphai= new int[gNmod_alpha];
			for(int i =0; i<gNmod_alpha ; i++){
				gIMR_alphai[i] = int_dict["delta alpha "+std::to_string(i)+" i"];
				std::cout<<i<<" : "<<gIMR_alphai[i]<<std::endl;
				ct++;
			}
		}
		
	}
	int total_mods = Nmod+gNmod_phi+gNmod_sigma+gNmod_beta+gNmod_alpha;
	MCMC_modification_struct mod_struct;
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
	double init_pos[dimension];
	double **fisher = new double*[dimension];
	for(int i = 0 ; i<dimension; i++){
		fisher[i] = new double[dimension];
		for(int j = 0 ; j<dimension; j++){
			fisher[i][j] = 0;
		}
	}
	for(int i = 0 ; i<dimension; i++){
		delete [] fisher[i];

	}
	if(has_substring(generation_method, "ppE")){
		delete [] bppe;	
	}
	delete [] fisher;
}
