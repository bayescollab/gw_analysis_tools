#include "mcmc_gw_extended.h"
#include "standardPriorLibrary.h"
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

#include <bayesship/bayesshipSampler.h>
#include <bayesship/dataUtilities.h>
#include <bayesship/utilities.h>


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


double EA_prior[6];
bool alpha_param; 

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
	int samples = int_dict["independent samples"];
	std::cout<<"Independent Samples: "<<samples<<std::endl;
	int ensembleSize = int_dict["ensemble size"];
	std::cout<<"ensemble size: "<<ensembleSize<<std::endl;
	int ensembleN = int_dict["ensemble number"];
	std::cout<<"ensemble number: "<<ensembleN<<std::endl;
	int chainN = ensembleN * ensembleSize;
	std::cout<<"Chain number: "<<chainN<<std::endl;
	int burnPriorIterations = int_dict["prior burn iterations"];
	std::cout<<"prior burn iterations: "<<burnPriorIterations<<std::endl;
	int priorIterations = int_dict["prior iterations"];
	std::cout<<"prior iterations: "<<priorIterations<<std::endl;
	int burnIterations = int_dict["burn iterations"];
	std::cout<<"burn iterations: "<<burnIterations<<std::endl;
	
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
	bool writePriorData = bool_dict["write prior data"];	
	std::cout<<"Write Prior Data : "<<writePriorData<<std::endl;
	
	double swapProb = dbl_dict["swap probability"];
	std::cout<<"swap probability: "<<swapProb<<std::endl;
	int threads = int_dict["thread number"];
	std::cout<<"threads: "<<threads<<std::endl;
	std::string outputDir = str_dict["output directory"];
	std::string outputMoniker = str_dict["output moniker"];
	std::cout<<"Output Directory: "<<outputDir<<std::endl;
	std::cout<<"Output moniker: "<<outputMoniker<<std::endl;
	int dimension = int_dict["dimension"];
	std::cout<<"Dimension: "<<dimension<<std::endl;
	int max_chunk_size = int_dict["Max chunk size"];
	std::cout<<"Max chunk size: "<<max_chunk_size<<std::endl;
	
	std::string initial_position_file="", initial_checkpoint_file="",initial_ensemble_position_file="";
	bool continue_from_checkpoint=false;


	priorData PD;

	if(str_dict.find("initial position file") != str_dict.end()){
		initial_position_file = str_dict["initial position file"];
		std::cout<<"Chain number: "<<chainN<<std::endl;
		std::cout<<"Initial position file: "<<initial_position_file<<std::endl;
	}
	if (str_dict.find("initial ensemble position file") != str_dict.end()){
		initial_ensemble_position_file = str_dict["initial ensemble position file"];
		std::cout<<"Chain number: "<<chainN<<std::endl;
		std::cout<<"Initial ensemble position file: "<<initial_ensemble_position_file<<std::endl;
	}
	if ( initial_position_file =="" && initial_ensemble_position_file ==""){
		std::cout<<"ERROR -- you need an initial checkpoint file, an initial position file, or an initial ensemble file"<<std::endl;
		exit(1);
	}

	if(dbl_dict.find("Mass1 minimum") == dbl_dict.end()){
		PD.mass1_prior[0]=1;
		PD.mass1_prior[1]=80;
	}
	else{
		PD.mass1_prior[0]=dbl_dict["Mass1 minimum"];
		PD.mass1_prior[1]=dbl_dict["Mass1 maximum"];
	}
	if(dbl_dict.find("Mass2 minimum") == dbl_dict.end()){
		PD.mass2_prior[0]=1;
		PD.mass2_prior[1]=80;
	}
	else{
		PD.mass2_prior[0]=dbl_dict["Mass2 minimum"];
		PD.mass2_prior[1]=dbl_dict["Mass2 maximum"];
	}
	if(dbl_dict.find("Spin1 minimum") == dbl_dict.end()){
		PD.spin1_prior[0]=-1;
		PD.spin1_prior[1]=1;
	}
	else{
		PD.spin1_prior[0]=dbl_dict["Spin1 minimum"];
		PD.spin1_prior[1]=dbl_dict["Spin1 maximum"];
	}
	if(dbl_dict.find("Spin2 minimum") == dbl_dict.end()){
		PD.spin2_prior[0]=-1;
		PD.spin2_prior[1]=1;
	}
	else{
		PD.spin2_prior[0]=dbl_dict["Spin2 minimum"];
		PD.spin2_prior[1]=dbl_dict["Spin2 maximum"];
	}
	if(dbl_dict.find("Luminosity distance minimum") == dbl_dict.end()){
		PD.DL_prior[0]=10;
		PD.DL_prior[1]=5000;
	}
	else{
		PD.DL_prior[0]=dbl_dict["Luminosity distance minimum"];
		PD.DL_prior[1]=dbl_dict["Luminosity distance maximum"];
	}
	std::cout<<"Range of Mass1: "<<PD.mass1_prior[0]<<" - "<<PD.mass1_prior[1]<<std::endl;
	std::cout<<"Range of Mass2: "<<PD.mass2_prior[0]<<" - "<<PD.mass2_prior[1]<<std::endl;
	std::cout<<"Range of DL: "<<PD.DL_prior[0]<<" - "<<PD.DL_prior[1]<<std::endl;
	if(generation_method.find("IMRPhenomPv2") != std::string::npos ){
		std::cout<<"Range of Spin1: "<<0<<" - "<<PD.spin1_prior[1]<<std::endl;
		std::cout<<"Range of Spin2: "<<0<<" - "<<PD.spin2_prior[1]<<std::endl;
	}
	else{
		std::cout<<"Range of Spin1: "<<PD.spin1_prior[0]<<" - "<<PD.spin1_prior[1]<<std::endl;
		std::cout<<"Range of Spin2: "<<PD.spin2_prior[0]<<" - "<<PD.spin2_prior[1]<<std::endl;

	}

	if(dbl_dict.find("tidal_s minimum") == dbl_dict.end()){
		PD.tidal_s_prior[0]=1;
		PD.tidal_s_prior[2]=5000;
	}
	else{
		PD.tidal_s_prior[0]=dbl_dict["tidal_s minimum"];
		PD.tidal_s_prior[1]=dbl_dict["tidal_s maximum"];
	}
	if(dbl_dict.find("tidal1 minimum") == dbl_dict.end()){
		PD.tidal1_prior[0]=1;
		PD.tidal1_prior[1]=5000;
	}
	else{
		PD.tidal1_prior[0]=dbl_dict["tidal1 minimum"];
		PD.tidal1_prior[1]=dbl_dict["tidal1 maximum"];
	}
	if(dbl_dict.find("tidal2 minimum") == dbl_dict.end()){
		PD.tidal2_prior[0]=1;
		PD.tidal2_prior[1]=5000;
	}
	else{
		PD.tidal2_prior[0]=dbl_dict["tidal2 minimum"];
		PD.tidal2_prior[1]=dbl_dict["tidal2 maximum"];
	}
	PD.tidal_love = true;
	if(bool_dict.find("Tidal love relation") != bool_dict.end())
	{
		PD.tidal_love = bool_dict["Tidal love relation"];
	}
	PD.tidal_love_error = false;
	if(bool_dict.find("Tidal love error marginalization") != bool_dict.end())
	{
		PD.tidal_love_error = bool_dict["Tidal love error marginalization"];
	}
	if(generation_method.find("NRT") != std::string::npos){
		std::cout<<"Range of tidal1: "<<PD.tidal1_prior[0]<<" - "<<PD.tidal1_prior[1]<<std::endl;
		std::cout<<"Range of tidal2: "<<PD.tidal2_prior[0]<<" - "<<PD.tidal2_prior[1]<<std::endl;
		std::cout<<"Range of tidal_s: "<<PD.tidal_s_prior[0]<<" - "<<PD.tidal_s_prior[1]<<std::endl;
		std::cout<<"Using tidal-love relations: "<<PD.tidal_love<<std::endl;
		std::cout<<"Using tidal-love error marginalization: "<<PD.tidal_love_error<<std::endl;
				
	}
	bool jeff_prior = false;
	if(bool_dict.find("Jefferys prior") == bool_dict.end())
	{
		jeff_prior = bool_dict["Jeffreys prior"];
	}
	//Einstein AEther parameterization
       	PD.alpha_param = true;
	if(bool_dict.find("alpha parameterization") != bool_dict.end())
	{
	        PD.alpha_param = bool_dict["alpha parameterization"];
	}
	//Einstein AEther priors
	if(generation_method.find("EA") != std::string::npos){
	  std::cout<<"Using alpha parameterization: "<<PD.alpha_param<<std::endl;
	  if(PD.alpha_param){
	    if(dbl_dict.find("EA alpha_1 minimum") ==dbl_dict.end()){
	      PD.EA_prior[0]= -1.;
	      PD.EA_prior[1]= 1.; 
	    }
	    else{
	      PD.EA_prior[0]=dbl_dict["EA alpha_1 minimum"];
	      PD.EA_prior[1]=dbl_dict["EA alpha_1 maximum"]; 
	    }
	    if(dbl_dict.find("EA alpha_2 minimum") ==dbl_dict.end()){
	      PD.EA_prior[2]=-.1;
	      PD.EA_prior[3]=.1; 
	    }
	    else{
	      PD.EA_prior[2]=dbl_dict["EA alpha_2 minimum"];
	      PD.EA_prior[3]=dbl_dict["EA alpha_2 maximum"]; 
	    }
	    if(dbl_dict.find("EA cbar_w minimum") ==dbl_dict.end()){
	      PD.EA_prior[4]=-1.;
	      PD.EA_prior[5]=1.; 
	    }
	    else{
	      PD.EA_prior[4]=dbl_dict["EA cbar_w minimum"];
	      PD.EA_prior[5]=dbl_dict["EA cbar_w maximum"]; 
	    }
	    std::cout<<"Range of EA alpha1: "<<PD.EA_prior[0]<<" - "<<PD.EA_prior[1]<<std::endl;
	    std::cout<<"Range of EA alpha2: "<<PD.EA_prior[2]<<" - "<<PD.EA_prior[3]<<std::endl;
	    std::cout<<"Range of EA cbarw: "<<PD.EA_prior[4]<<" - "<<PD.EA_prior[5]<<std::endl;
	  }
	  else{
	    if(dbl_dict.find("EA c_a minimum") == dbl_dict.end()){
	      PD.EA_prior[0]=0;
	      PD.EA_prior[1]=pow(10,-4.);
	    }
	    else{
	      PD.EA_prior[0]=dbl_dict["EA c_a minimum"];
	      PD.EA_prior[1]=dbl_dict["EA c_a maximum"];
	    }
	    if(dbl_dict.find("EA c_theta minimum") == dbl_dict.end()){
	      PD.EA_prior[2]=0;
	      PD.EA_prior[3]=pow(10,-4.);
		}
	}
	}
	bool restrictSwapTemperatures = false;
	if(bool_dict.find("restrict swap temperatures") != bool_dict.end()){
		restrictSwapTemperatures = bool_dict["restrict swap temperatures"];
	}

	bool ignoreExistingCheckpoint = false;
	if(bool_dict.find("ignore existing checkpoint") != bool_dict.end()){
		ignoreExistingCheckpoint = bool_dict["ignore existing checkpoint"];
	}
	bool coldChainStorageOnly = true;
	if(bool_dict.find("cold chain storage only") != bool_dict.end()){
		coldChainStorageOnly = bool_dict["cold chain storage only"];
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
		PD.RA_bounds[0] = dbl_dict["RA minimum"];
		PD.RA_bounds[1] = dbl_dict["RA maximum"];
		std::cout<<"RA bounds: "<<PD.RA_bounds[0]<<" "<<PD.RA_bounds[1]<<std::endl;
	}
	else{
		PD.RA_bounds[0] = 0;
		PD.RA_bounds[1] = 2*M_PI;

	}
	if( dbl_dict.find("Sin(DEC) minimum") != dbl_dict.end() && dbl_dict.find("Sin(DEC) maximum") != dbl_dict.end()){
		PD.sinDEC_bounds[0] = dbl_dict["Sin(DEC) minimum"];
		PD.sinDEC_bounds[1] = dbl_dict["Sin(DEC) maximum"];
		std::cout<<"Sin(DEC) bounds: "<<PD.sinDEC_bounds[0]<<" "<<PD.sinDEC_bounds[1]<<std::endl;
	}
	else{
		PD.sinDEC_bounds[0] = -1;
		PD.sinDEC_bounds[1] = 1;

	}
	

	allocate_LOSC_data(detector_files, psd_file, detector_N, psd_length, data_length, gps_time,post_merger_duration, data, psd, freqs);

	PD.T_mcmc_gw_tool = 1./(freqs[0][2]-freqs[0][1]);
	PD.T_merger = PD.T_mcmc_gw_tool - post_merger_duration;
	double df = 1./PD.T_mcmc_gw_tool;

	std::cout<<"Total time: "<<PD.T_mcmc_gw_tool<<std::endl;
	int *data_lengths= (int*)malloc(sizeof(int)*detector_N);
	for(int i = 0 ; i<detector_N; i++){
		data_lengths[i] =psd_length;
	}

	double **output;
	output = allocate_2D_array(samples, dimension );
	double chain_temps[chainN];
	
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
		PD.mod_priors = new double*[Nmod];
		for(int i =0; i<Nmod ; i++){
			bppe[i] = dbl_dict["ppE b "+std::to_string(i)];
			PD.mod_priors[i]= new double[2];
			std::cout<<i<<" : "<<bppe[i]<<std::endl;
			PD.mod_priors[i][0] = dbl_dict["ppE beta "+std::to_string(i)+" minimum"];
			PD.mod_priors[i][1] = dbl_dict["ppE beta "+std::to_string(i)+" maximum"];
			std::cout<<"Min"<<" : "<<PD.mod_priors[i][0]<<std::endl;
			std::cout<<"Max"<<" : "<<PD.mod_priors[i][1]<<std::endl;
		}
		
	}
	if(generation_method.find("gIMR") != std::string::npos){
		gNmod_phi = int_dict["Number of phi modifications"];
		gNmod_sigma = int_dict["Number of sigma modifications"];
		gNmod_beta = int_dict["Number of beta modifications"];
		gNmod_alpha = int_dict["Number of alpha modifications"];
		int Nmod_tot = gNmod_phi+ gNmod_sigma+gNmod_beta+gNmod_alpha;
		std::cout<<"Number of general modifications: "<<Nmod_tot<<std::endl;
		std::cout<<"delta phi i: "<<Nmod<<std::endl;
		PD.mod_priors = new double*[Nmod_tot];
		int ct = 0;
		if(gNmod_phi != 0){
			gIMR_phii= new int[gNmod_phi];
			for(int i =0; i<gNmod_phi ; i++){
				gIMR_phii[i] = int_dict["delta phi "+std::to_string(i)+" i"];
				PD.mod_priors[ct]= new double[2];
				std::cout<<i<<" : "<<gIMR_phii[i]<<std::endl;
				PD.mod_priors[ct][0] = dbl_dict["delta phi "+std::to_string(i)+" minimum"];
				PD.mod_priors[ct][1] = dbl_dict["delta phi "+std::to_string(i)+" maximum"];
				ct++;
			}
		}
		if(gNmod_sigma != 0){
			gIMR_sigmai= new int[gNmod_sigma];
			for(int i =0; i<gNmod_sigma ; i++){
				gIMR_sigmai[i] = int_dict["delta sigma "+std::to_string(i)+" i"];
				PD.mod_priors[ct]= new double[2];
				std::cout<<i<<" : "<<gIMR_sigmai[i]<<std::endl;
				PD.mod_priors[ct][0] = dbl_dict["delta sigma "+std::to_string(i)+" minimum"];
				PD.mod_priors[ct][1] = dbl_dict["delta sigma "+std::to_string(i)+" maximum"];
				ct++;
			}
		}
		if(gNmod_beta != 0){
			gIMR_betai= new int[gNmod_beta];
			for(int i =0; i<gNmod_beta ; i++){
				gIMR_betai[i] = int_dict["delta beta "+std::to_string(i)+" i"];
				PD.mod_priors[ct]= new double[2];
				std::cout<<i<<" : "<<gIMR_betai[i]<<std::endl;
				PD.mod_priors[ct][0] = dbl_dict["delta beta "+std::to_string(i)+" minimum"];
				PD.mod_priors[ct][1] = dbl_dict["delta beta "+std::to_string(i)+" maximum"];
				ct++;
			}
		}
		if(gNmod_alpha != 0){
			gIMR_alphai= new int[gNmod_alpha];
			for(int i =0; i<gNmod_alpha ; i++){
				gIMR_alphai[i] = int_dict["delta alpha "+std::to_string(i)+" i"];
				PD.mod_priors[ct]= new double[2];
				std::cout<<i<<" : "<<gIMR_alphai[i]<<std::endl;
				PD.mod_priors[ct][0] = dbl_dict["delta alpha "+std::to_string(i)+" minimum"];
				PD.mod_priors[ct][1] = dbl_dict["delta alpha "+std::to_string(i)+" maximum"];
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
	mod_struct.tidal_love = PD.tidal_love;
	mod_struct.tidal_love_error = PD.tidal_love_error;
	mod_struct.alpha_param = PD.alpha_param; 
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

	std::cout<<"dimension - total_mods="<<dimension-total_mods<<std::endl;
	double(*lp)(double *param, mcmc_data_interface *interface, void *parameters);
	bayesship::probabilityFn *logp;

	if(generation_method.find("IMRPhenomD") != std::string::npos && (dimension-total_mods) == 11){
		if(total_mods == 0){
			std::cout<<"Using standard all-sky IMRPhenomD prior"<<std::endl;
			logp = new logPriorStandard_D(&PD);
		}
		//else{
		//	lp = &standard_log_prior_D_mod;
		//}
	}
	else if(generation_method.find("IMRPhenomD_NRT") != std::string::npos && ( (dimension-total_mods) == 13 || (dimension-total_mods) == 12)){
		if(total_mods == 0){
			std::cout<<"Using standard all-sky IMRPhenomD/NRT prior"<<std::endl;
			logp = new logPriorStandard_D_NRT(&PD);
		}
		else if(generation_method.find("EA") !=std::string::npos){
			std::cout<<"Using standard all-sky IMRPhenomD/NRT/EA prior"<<std::endl;
			logp = new logPriorStandard_D_NRT_EA(&PD);
		}
		else if(generation_method.find("ppE") != std::string::npos || generation_method.find("ppE") != std::string::npos || check_theory_support(generation_method)){
			std::cout<<"Using standard all-sky IMRPhenomD/NRT/Mod prior"<<std::endl;
			logp = new logPriorStandard_D_NRT_mod(&PD);
		}
	}
	else{
		std::cout<<"ERROR -- wrong detector/dimension combination for this tool -- Check mcmc_gw for general support"<<std::endl;
		return 1;
	}

	double **initial_position=NULL;
	initial_position = new double*[1];
	initial_position[0] = new double[dimension];
	double *seeding_var = NULL;
	double **ensemble_initial_position = NULL;
	bayesship::positionInfo *initialPosition = nullptr;
	bayesship::positionInfo **ensembleInitialPosition = nullptr;

	if(initial_position_file != ""){
		initialPosition = new bayesship::positionInfo(dimension,false);
		read_file(initial_position_file, initial_position,1,dimension);
		for(int i = 0 ; i<dimension; i++){
			initialPosition->parameters[i] = initial_position[0][i];
		}
	}
	if(initial_ensemble_position_file != ""){
		ensembleInitialPosition = new bayesship::positionInfo*[chainN];
		for(int i = 0 ; i<chainN; i++){
			ensembleInitialPosition[i] = new bayesship::positionInfo(dimension,false);
		}
		ensemble_initial_position = new double*[chainN];
		for(int i = 0 ; i<chainN; i++){
			ensemble_initial_position[i] = new double[dimension];
		}
		read_file(initial_ensemble_position_file, ensemble_initial_position,chainN,dimension);

		for(int j = 0 ; j<chainN; j++){
			for(int i = 0 ; i<dimension; i++){
				ensembleInitialPosition[j]->parameters[i] = ensemble_initial_position[j][i];
			}
		}
	}
	//bayesship::probabilityFn *logp = new logPriorStandard_D(&PD);
	std::cout<<"Running uncorrelated sampler "<<std::endl;
	bayesship::bayesshipSampler *sampler = PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW_v2( 
			dimension, samples, ensembleSize, ensembleN, 
			 initialPosition,ensembleInitialPosition,swapProb, 
			 burnIterations, burnPriorIterations,priorIterations, writePriorData,max_chunk_size, (double **)nullptr,
			logp,threads, pool,detector_N, 
			data, psd,freqs, data_lengths,gps_time, detectors,&mod_struct,
			generation_method,outputDir, outputMoniker,ignoreExistingCheckpoint,restrictSwapTemperatures,coldChainStorageOnly);	
	delete logp;
	delete [] initial_position[0]; delete [] initial_position;
	if(initial_ensemble_position_file != ""){
		for(int i = 0 ; i<chainN; i++){
			delete [] ensemble_initial_position[i]; 
		}
		delete [] ensemble_initial_position;
	}
	if(initialPosition){
		delete initialPosition;
	}
	if(ensembleInitialPosition){
		for(int i = 0 ; i<chainN; i++){
			delete ensembleInitialPosition[i];
		}
		delete [] ensembleInitialPosition;
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
			delete [] PD.mod_priors[i] ;
		}
		delete [] PD.mod_priors;
	}
	else if(generation_method.find("EA") !=std::string::npos){
		for(int i = 0 ; i<4 ; i++){
			delete [] PD.mod_priors[i] ;
		}
		delete [] PD.mod_priors;

	}
	
	delete sampler;
	
	return 0;

}
