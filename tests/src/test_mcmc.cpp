#include <gwat/waveform_util.h>
#include <gwat/ortho_basis.h>
#include <gwat/io_util.h>
#include <gwat/mcmc_sampler.h>
#include <gwat/mcmc_gw.h>
#include <gwat/mcmc_gw_extended.h>
#include <gwat/detector_util.h>
#include <gwat/util.h>
#include <gwat/waveform_util.h>
#include <gwat/IMRPhenomP.h>
#include <gwat/fisher.h>
#include <iostream>
#include <gsl/gsl_sf_gamma.h>


double root2 = std::sqrt(2.);

void fisher_1d_test(double *c, double **fisher, mcmc_data_interface *interface, void *parameters);
int test_1d(int argc, char *argv[]);
double log_1d_test (double *c,mcmc_data_interface *interface,void *parameters);
double log_1d_test_prior (double *c,mcmc_data_interface *interface,void *parameters);
double **multi_gaussian_prior_cov;
double **multi_gaussian_like_cov;
double **multi_gaussian_like_fisher;
double **multi_gaussian_prior_fisher;
double *multi_gaussian_prior_mean;
double *multi_gaussian_like_mean;
double multi_gaussian_scale = 1.;

int stupid_test(int argc, char *argv[]);

void RT_ERROR_MSG();
int test_likelihood(int argc, char *argv[]);
int mcmc_standard_test(int argc, char *argv[]);
int mcmc_rosenbock(int argc, char *argv[]);
int mcmc_output_class(int argc, char *argv[]);
int mcmc_RJ_sin(int argc, char *argv[]);
int multiple_continue(int argc, char *argv[]);
void RJ_sin_fish(double *c,int *status,int model_status,double **fisher,mcmc_data_interface *interface, void *parameters);
//double RJ_sin_logL(double *params, mcmc_data_interface *interface, void *parameters);
double RJ_sin_logL(double *params, int *status,int model_status,mcmc_data_interface *interface, void *parameters);
void RJ_sin_proposal(double *current_params, double *proposed_params, int *current_status, int *proposed_status,int *current_model_status, int *proposed_model_status, double *MH_corrections,mcmc_data_interface *interface, void *parameters);
//double RJ_sin_prior(double *params,  mcmc_data_interface *interface, void *parameters);
double RJ_sin_prior(double *params,  int *status,int model_status,mcmc_data_interface *interface, void *parameters);
double log_test (double *c,mcmc_data_interface *interface,void *parameters);
double log_test_prior (double *c,mcmc_data_interface *interface,void *parameters);
double log_rosenbock (double *c,mcmc_data_interface *interface,void *parameters);
double log_rosenbock_prior (double *c,mcmc_data_interface *interface,void *parameters);
void fisher_test(double *c, double **fisher, mcmc_data_interface *interface, void *parameters);
double standard_log_prior_D(double *pos, mcmc_data_interface *interface,void *parameters);
double standard_log_prior_D_sa(double *pos, mcmc_data_interface *interface,void *parameters);
double standard_log_prior_P(double *pos, mcmc_data_interface *interface,void *parameters);
double standard_log_prior_P_RJ(double *pos,int *status,int model_status, mcmc_data_interface *interface,void *parameters);
double standard_log_prior_D_RJ(double *pos,int *status,int model_status, mcmc_data_interface *interface,void *parameters);
double chirpmass_eta_jac(double chirpmass, double eta);
double chirpmass_q_jac(double chirpmass, double q);
double T_mcmc_gw_tool ;
double standard_log_prior_dCS(double *pos, mcmc_data_interface *interface,void *parameters);
//void fisher_rosenbock(double *c,double **fisher,  mcmc_data_interface *interface,void *parameters);

int test_evidence(int argc, char *argv[]);
int ensemble_size(int argc, char *argv[]);
int validate_evidence(int argc, char *argv[]);
double log_prior_multi_gaussian(double *param, mcmc_data_interface *interface, void *parameters);
double log_like_multi_gaussian(double *param, mcmc_data_interface *interface, void *parameters);
void fisher_multi_gaussian(double *param, double **fisher, mcmc_data_interface *interface, void *parameters);
struct RJ_sin_param
{
	double *data;
	int N;
	gsl_rng *r;
};


#include <ptrjmcmc/PtrjmcmcSampler.h>
#include <ptrjmcmc/dataUtilities.h>

int test_ptrjmcmc_integration(int argc, char *argv[]);
double test_ptrjmcmc_integration_log_prior(ptrjmcmc::positionInfo *pos, int chainID, ptrjmcmc::PtrjmcmcSampler  *sampler, void *userParameters);


int main(int argc, char *argv[])
{
	std::cout<<"TESTING MCMC CALCULATIONS"<<std::endl;
	if(argc != 2){
		RT_ERROR_MSG();
		return 1;
	}
	
	int runtime_opt = std::stoi(argv[1]);	
	if(runtime_opt == 0){
		std::cout<<"MCMC student t"<<std::endl;
		return mcmc_standard_test(argc,argv);
	}
	else if(runtime_opt == 1){
		std::cout<<"MCMC Injected GW"<<std::endl;
		std::cout<<"Test function deprecated due to function PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW no longer being in use"<<std::endl;
		return -1;
	}
	else if(runtime_opt == 2){
		std::cout<<"MCMC data GW"<<std::endl;
		std::cout<<"Test Function depcrecated due to function PTMCMC_MH_GW no longer being in use"<<std::endl;
		return -1;
	}
	else if(runtime_opt == 3){
		std::cout<<"MCMC N-dimensional Rosenbock"<<std::endl;
		return mcmc_rosenbock(argc,argv);
	}
	else if(runtime_opt == 4){
		std::cout<<"MCMC output class testing"<<std::endl;
		return mcmc_output_class(argc,argv);
	}
	else if(runtime_opt == 5){
		std::cout<<"RJ sine wave test"<<std::endl;
		return mcmc_RJ_sin(argc,argv);
	}
	else if(runtime_opt == 6){
		std::cout<<"MCMC Injected GW RJ"<<std::endl;
		std::cout<<"Test Function deprecated due to function RJPTMCMC_MH_dynamic_PT_alloc_comprehensive_2WF_GW no longer being in use"<<std::endl;
		return -1;
	}
	else if(runtime_opt == 7){
		std::cout<<"Continue MCMC validation"<<std::endl;
		return multiple_continue(argc,argv);
	}
	else if(runtime_opt == 8){
		std::cout<<"Test Likelihood"<<std::endl;
		return test_likelihood(argc,argv);
	}
	else if(runtime_opt == 9){
		std::cout<<"Test evidence calculation"<<std::endl;
		return test_evidence(argc,argv);
	}
	else if(runtime_opt == 10){
		std::cout<<"Validate evidence calculation"<<std::endl;
		return validate_evidence(argc,argv);
	}
	else if(runtime_opt == 11){
		std::cout<<"Test chain ensemble size"<<std::endl;
		return ensemble_size(argc,argv);
	}
	else if(runtime_opt == 12){
		std::cout<<"Stupid test"<<std::endl;
		return stupid_test(argc,argv);
	}
	else if(runtime_opt == 13){
		std::cout<<"1D test"<<std::endl;
		return test_1d(argc,argv);
	}
	else if(runtime_opt == 14){
		std::cout<<"PTRJMCMC sampler testing -- GW injection"<<std::endl;
		return test_ptrjmcmc_integration(argc,argv);
	}
	else{
		RT_ERROR_MSG();
		return 1;
	}
}

int test_ptrjmcmc_integration(int argc, char *argv[])
{
	gen_params injection;
	injection.mass1 = 8.4;
	injection.mass2 = 7.3;
	double chirpmass = calculate_chirpmass(injection.mass1,injection.mass2);
	double eta = calculate_eta(injection.mass1,injection.mass2);
	double q = injection.mass2/injection.mass1;
	injection.Luminosity_Distance =1200;
	//injection.Luminosity_Distance =300;
	injection.psi = .2;
	injection.phiRef = 2.;
	injection.f_ref = 20.;
	injection.RA = 3.5;
	injection.DEC = -.8;
	injection.spin1[2] = .0;
	injection.spin2[2] = .0;
	injection.spin1[1] = .0;
	injection.spin2[1] = -.0;
	injection.spin1[0] = .00;
	injection.spin2[0] = -.00;
	injection.incl_angle = -.91;
	double gps = 1135136350.6;
	injection.gmst = gps_to_GMST_radian(gps);
	std::cout<<"GMST: "<<injection.gmst<<std::endl;
	injection.equatorial_orientation = false;
	injection.horizon_coord = false;
	injection.shift_time = true;
	injection.shift_phase = true;

	IMRPhenomPv2<double> model;
	injection.chip = model.PhenomPv2_inplane_spin(&injection);
	injection.phip = 1;


	//injection.Nmod = 1;
	//injection.bppe = new double[1];
	//injection.betappe = new double[1];
	//injection.bppe[0] = -1;
	//injection.betappe[0] = 0;
	
	int detect_number = 4;
	std::string detectors[4] = {"Hanford","Livingston","Virgo","Kagra"};
	std::string SN[4] = {"AdLIGODesign","AdLIGODesign","AdLIGODesign","KAGRA_pess"};
	std::string injection_method = "IMRPhenomD";
	double fmin = 10;
	double fmax =1024;
	T_mcmc_gw_tool= 4;
	int length ;
	double **psd;
	double **freq;
	double **weights;
	std::complex<double> **data ;
	double tc_ref=0;
	tc_ref = T_mcmc_gw_tool*(1-3./4.);
	//double tc_ref = T_mcmc_gw_tool*(3./4.);
	double deltaf = 1./T_mcmc_gw_tool;
	length = (fmax-fmin)/deltaf;
	data = new std::complex<double>*[detect_number];
	psd = new double*[detect_number];
	freq = new double*[detect_number];
	double total_snr = 0;
	double tc=0;
	injection.tc = 0;
	for(int i = 0 ; i<detect_number; i++){
		data[i]= new std::complex<double>[length];
		psd[i]= new double[length];
		freq[i]= new double[length];
		for(int j =0 ; j<length; j++){
			freq[i][j] = fmin + j*deltaf;
		}
		populate_noise(freq[i],SN[i],psd[i],length);
		for(int j =0 ; j<length; j++){
			psd[i][j] *= psd[i][j];
		}
		double deltat = DTOA_DETECTOR(injection.RA,injection.DEC,injection.gmst,detectors[0],detectors[i]);
		//injection.tc = tc_ref + deltat;
		tc = tc_ref - deltat;
		tc*=2*M_PI;
		fourier_detector_response(freq[i],length, data[i],detectors[i],injection_method, &injection, (double *)NULL);
		for(int j = 0 ; j<length; j++){
			data[i][j]*=exp(std::complex<double>(0,tc*freq[i][j]));
		}
		double snr =  calculate_snr_internal(psd[i],data[i],freq[i],length, "SIMPSONS",(double*) NULL, false);
		std::cout<<detectors[i]<<" "<<snr<<std::endl;
		total_snr += pow_int(snr,2);
	}
	injection.tc=tc_ref;
	std::cout<<"NETWORK SNR of injection: "<<sqrt(total_snr)<<std::endl;

	int data_lengths[detect_number];
	for(int i = 0 ; i<detect_number; i++){
		data_lengths[i]=length;	
	}
	//################################################################
	//################################################################
	//################################################################
	//################################################################
	int dim = 11;
	std::string recovery_method = "IMRPhenomD";
	
	double spin1sph[3];
	double spin2sph[3];
	transform_cart_sph(injection.spin1,spin1sph);
	transform_cart_sph(injection.spin2,spin2sph);

	double initial_position[dim]= {injection.RA, sin(injection.DEC),injection.psi, cos(injection.incl_angle), injection.phiRef, T_mcmc_gw_tool-tc_ref, log(injection.Luminosity_Distance),log(chirpmass), eta, injection.spin1[2],injection.spin2[2]};
	//double initial_position[dim]= {log(chirpmass), eta, injection.spin1[2],injection.spin2[2]};

	write_file("data/injections_PTRJMCMC.csv",initial_position,dim);

	ptrjmcmc::positionInfo *initialPosition = new ptrjmcmc::positionInfo(dim,false);
	for(int i = 0  ;  i<dim; i++){
		initialPosition->parameters[i] = initial_position[i];
	}

	double **priorRanges = new double*[dim];
	for(int i = 0 ; i<dim; i++){priorRanges[i] = new double[2];}

	priorRanges[0][0] = 0;
	priorRanges[0][1] = 2*M_PI;

	priorRanges[1][0] = -1;
	priorRanges[1][1] = 1;

	priorRanges[2][0] = 0;
	priorRanges[2][1] = M_PI;

	priorRanges[3][0] = -1;
	priorRanges[3][1] = 1;

	priorRanges[4][0] = 0;
	priorRanges[4][1] = 2*M_PI;

	priorRanges[5][0] = 2;
	priorRanges[5][1] = 4;

	priorRanges[6][0] = std::log(10);
	priorRanges[6][1] = std::log(5000);

	priorRanges[7][0] = std::log(1);
	priorRanges[7][1] = std::log(40);

	priorRanges[8][0] = 0.01;
	priorRanges[8][1] = .25;

	priorRanges[9][0] = -1;
	priorRanges[9][1] = 1;

	priorRanges[10][0] = -1;
	priorRanges[10][1] = 1;



	int ensembleSize = 15;
	int ensembleN = 5;
	int swapProb = .1;
	int threads = 6;
	bool pool = true;
	std::string outputDir = "data/";
	std::string outputMoniker = "PTRJMCMC_GW_injection";
	MCMC_modification_struct mod_struct;
	mod_struct.ppE_Nmod = 0;
	
	int samples = 1000;
	//int samples = 100;
	double burnIterations = 50000;
	//double burnIterations = 100;
	double burnPriorIterations = 0;
	//double burnPriorIterations = 0;
	double priorIterations =0 ;
	//double priorIterations =0 ;
	int max_chunk_size = 1e6;

	ptrjmcmc::PtrjmcmcSampler * samplerObj = PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW_v2(
		dim, samples, ensembleSize, ensembleN, initialPosition, (ptrjmcmc::positionInfo **)nullptr, swapProb, burnIterations, burnPriorIterations,priorIterations, true, max_chunk_size, &(priorRanges[0]), test_ptrjmcmc_integration_log_prior, threads, pool, detect_number, data, psd, freq, data_lengths, gps, detectors, &mod_struct, recovery_method, outputDir, outputMoniker);
	//samplerObj->data->writeStatFile(outputDir+outputMoniker+"Stat.csv");
	if(priorIterations>0){
		samplerObj->priorData->writeStatFile(outputDir+outputMoniker+"PriorStat.csv");
	}
	delete initialPosition;

	delete samplerObj;
	for(int i = 0 ; i<dim ; i++){
		delete [] priorRanges[i];
	}
	delete [] priorRanges;
	
	for(int i = 0 ; i<detect_number; i++){
		delete [] freq[i];	
		delete [] psd[i];	
		delete [] data[i];	
	}
	delete [] injection.bppe;
	delete [] injection.betappe;
	delete [] freq;
	delete [] psd;
	delete [] data;

	return 0;
}
double aligned_spin_prior(double chi){
	double a=0.0039132 , b= 3.95381;
	return a * exp(-b * abs(chi));	
}

double test_ptrjmcmc_integration_log_prior(ptrjmcmc::positionInfo *posInfo, int chainID, ptrjmcmc::PtrjmcmcSampler  *sampler, void *userParameters)
{
	int dim = sampler->maxDim;
	double pos[dim];
	for(int i = 0 ; i<dim; i++){
		pos[i] = posInfo->parameters[i];
	}
	double a = -std::numeric_limits<double>::infinity();
	//###########
	double chirp = exp(pos[7]);
	double eta = pos[8];
	if (eta<.0 || eta>.25){return a;}//eta
	double m1 = calculate_mass1(chirp,eta );
	double m2 = calculate_mass2(chirp,eta );
	if(m1<1 || m1>20){return a;}
	if(m2<1 || m2>20){return a;}
	//###########
	if ((pos[0])<0 || (pos[0])>2*M_PI){ return a;}//RA
	if ((pos[1])<-1 || (pos[1])>1){return a;}//sinDEC
	if ((pos[2])<0 || (pos[2])>M_PI){return a;}//PSI
	if ((pos[3])<-1 || (pos[3])>1){return a;}//cos \iota
	if ((pos[4])<0 || (pos[4])>2*M_PI){return a;}//phiRef
	if( pos[5] < (3 - .1) || pos[5] > (3 + .1)) { return a; }
	if (std::exp(pos[6])<10 || std::exp(pos[6])>5000){return a;}//DL
	if ((pos[9])<-.95 || (pos[9])>.95){return a;}//chi1 
	if ((pos[10])<-.95 || (pos[10])>.95){return a;}//chi2
	//return log(chirpmass_eta_jac(chirp,eta))+3*pos[6] ;
	return log(aligned_spin_prior(pos[9]))+log(aligned_spin_prior(pos[10])) + log(chirpmass_eta_jac(chirp,eta))+3*pos[6] ;
	//return -.5*pow_int(pos[0]-5,2)/.01/.01-.5*pow_int(pos[1]+.9,2)/.01/.01+log(aligned_spin_prior(pos[9]))+log(aligned_spin_prior(pos[10])) + log(chirpmass_eta_jac(chirp,eta))+3*pos[6] ;


}


//double test_ptrjmcmc_integration_log_prior_aligned(ptrjmcmc::positionInfo *posInfo, int chainID, ptrjmcmc::PtrjmcmcSampler  *sampler, void *userParameters)
//{
//	int dim = sampler->maxDim;
//	double pos[dim];
//	for(int i = 0 ; i<dim; i++){
//		pos[i] = posInfo->parameters[i];
//	}
//	
//	double a = -std::numeric_limits<double>::infinity();
//	//###########
//	double chirp = exp(pos[0]);
//	double eta = pos[1];
//	if (eta<.0 || eta>.25){return a;}//eta
//	double m1 = calculate_mass1(chirp,eta );
//	double m2 = calculate_mass2(chirp,eta );
//	if(m1<1 || m1>40){return a;}
//	if(m2<1 || m2>40){return a;}
//	//###########
//	if ((pos[2])<-.95 || (pos[2])>.95){return a;}//chi1 
//	if ((pos[3])<-.95 || (pos[3])>.95){return a;}//chi2
//	//else {return log(chirpmass_eta_jac(chirp,eta)) ;}
//	return log(aligned_spin_prior(pos[2]))+log(aligned_spin_prior(pos[3])) + log(chirpmass_eta_jac(chirp,eta));
//
//}



int test_1d(int argc, char *argv[])
{
	int dimension = 1;
	//double initial_pos[2]={90,90.};
	double initial_pos[1]={0};
	double *seeding_var = NULL;
	int N_steps = 50000;
	int chain_N= 9;
	int max_chain_N= 3;
	//double *initial_pos_ptr = initial_pos;
	int swp_freq = 5;
	//double chain_temps[chain_N] ={1,2,3,10,12};
	double chain_temps[chain_N];
	chain_temps[0] = 1.;
	double c = 1.8;
	for(int i =1; i < chain_N/2;  i ++)
		chain_temps[i] =  chain_temps[i-1] * c;
	chain_temps[chain_N/2] = 1;
	for(int i =chain_N/2+1; i < chain_N;  i ++)
		chain_temps[i] =  chain_temps[i-1] * c;
	std::string autocorrfile = "";
	std::string chainfile = "data/mcmc_output_1d.hdf5";
	std::string statfilename = "data/mcmc_statistics_1d.txt";
	std::string checkpointfile = "data/mcmc_checkpoint_1d.csv";
	//std::string LLfile = "data/mcmc_LL.csv";
	std::string LLfile = "";
	
	int numThreads = 8;
	bool pool = true;
	bool show_progress = true;
	
	//double ***output;
	//output = allocate_3D_array( chain_N, N_steps, dimension );
	//PTMCMC_MH(output, dimension, N_steps, chain_N, initial_pos,seeding_var,chain_temps, swp_freq, log_test_prior, log_test,fisher_test,(void **)NULL,numThreads, pool,show_progress, statfilename,chainfile,autocorrfile, "",checkpointfile );	
	//deallocate_3D_array(output, chain_N, N_steps, dimension);
	double **output;
	output = allocate_2D_array(  N_steps, dimension );
	int t0 = 1000;
	int nu = 100;
	std::string chain_distribution_scheme="double";
	int max_chunk_size = 10000;
	mcmc_sampler_output sampler_output(chain_N,dimension);
	PTMCMC_MH_dynamic_PT_alloc_uncorrelated(&sampler_output,output, dimension, N_steps, chain_N, max_chain_N,initial_pos,seeding_var,(double**)NULL,chain_temps, swp_freq, t0,nu, max_chunk_size,chain_distribution_scheme, log_1d_test_prior, log_1d_test,fisher_1d_test,(void **)NULL,numThreads, pool,show_progress, statfilename,chainfile, LLfile,checkpointfile );	
	//sampler_output.create_data_dump(false,false,chainfile);
	deallocate_2D_array(output, N_steps, dimension);
	//sampler_output.write_flat_thin_output( "./data/test_flat_standard",true,true);
	std::cout<<"ENDED"<<std::endl;


		
	return 0;

}
double log_1d_test (double *c,mcmc_data_interface *interface,void *parameters)
{
	return 2;
	double x = c[0];
	//return 2.;
}
double log_1d_test_prior (double *c,mcmc_data_interface *interface,void *parameters)
{
	double a = -std::numeric_limits<double>::infinity();
	if( fabs(c[0]) > 10){ return a;}
	//return -pow_int(c[0]-4,2)/4. -pow_int(c[1]+3,2)/6.;
	return 0.;
}

void fisher_1d_test(double *c, double **fisher, mcmc_data_interface *interface, void *parameters)
{
	
	fisher[0][0] = 1./10/10;
	return;
} 
//########################################################################
//########################################################################
double stupid_test_log_prior(double *param, mcmc_data_interface *interface, void *parameters)
{
	double a = -std::numeric_limits<double>::infinity();
	for(int i = 0 ; i<interface->max_dim; i++){
		if( fabs(param[i]) > 10){ return a;}
	}
	return 0;
}

double stupid_test_log_likelihood(double *param, mcmc_data_interface *interface, void *parameters)
{
	return 2;

}
void stupid_test_fisher(double *param, double **fisher, mcmc_data_interface *interface, void *parameters)
{
	double return_val=0;
	for(int i = 0 ; i<interface->max_dim; i++){
		for(int j = 0 ; j<interface->max_dim; j++){
			if(i == j){

				fisher[i][j] = 1;
			}
			else{	
				fisher[i][j] = 0;
			}
		}
	}

}

int stupid_test(int argc, char *argv[])
{
	int total_samples = 2e4;
	int dimension = 20;	
	int chain_N = 10;
	int max_chain_N_thermo_ensemble = 5;	
	int N_steps = total_samples / ( (double)chain_N / max_chain_N_thermo_ensemble);
	std::cout<<"Number of steps per chain: "<<N_steps<<std::endl;
	double *seeding_var=NULL;
	double chain_temps[chain_N];
	int swp_freq = 5;
	int t0 = 3000;
	int nu = 100;
	int max_chunk_size = 100000;	
	std::string chain_distribution_scheme = "double";
	void ** user_parameters = NULL;
	int numthreads = 10;
	bool pool = true;
	bool show_prog = true;
	std::string modifier="";
	std::string stat_file = "data/stat_stupid_test"+modifier+".txt";
	std::string chain_file = "data/output_stupid_test"+modifier+".hdf5";
	std::string likelihood_file = "data/likelihood_stupid_test"+modifier+".csv";
	std::string check_file = "data/checkpoint_stupid_test"+modifier+".csv"	;
	mcmc_sampler_output sampler_output(chain_N, dimension) ;
	double **output = allocate_2D_array(N_steps,dimension);


	//double *init_pos = &multi_gaussian_like_mean[0];
	double *init_pos = new double[dimension];
	for(int i = 0 ;i<dimension; i++){
		//init_pos[i] = -10;
		init_pos[i] = 0;
	}
	
	PTMCMC_MH_dynamic_PT_alloc_uncorrelated(&sampler_output,output, dimension, N_steps, chain_N, max_chain_N_thermo_ensemble, init_pos, seeding_var, (double**)NULL,chain_temps, swp_freq, t0, nu, max_chunk_size, chain_distribution_scheme, stupid_test_log_prior, stupid_test_log_likelihood, stupid_test_fisher, user_parameters, numthreads, pool, show_prog, stat_file, chain_file, likelihood_file, check_file);
	
	deallocate_2D_array(output, chain_N,N_steps);
	delete [] init_pos;
	return 0;
}
//########################################################################
//########################################################################
struct trans_helper
{
	int N;
	double dt;
	double * data;
	gsl_rng *r;
};
double Chebyshev_fn(int P, double *coeff, double x)
{
	double sum = 0 ;
	for (int i = 0 ; i<P; i++){
		sum += coeff[i] * std::cos( i * std::acos(x));
	}
	return sum;
}
double transdimensional_likelihood(
	double *param, 
	int *status, 
	int model_status, 
	mcmc_data_interface *interface, 
	void * parameters)
{
	//return 1;
	//Sigma first, always there
	//start at 1 so there's always at least one coeff
	int P = 1;
	for(int i = 2 ; i<interface->max_dim; i++){
		if(status[i] == 0){
			break;
		}
		P+=1;
	}
	trans_helper *h = (trans_helper *)parameters;
	double reconst_signal[h->N];
	double dn =2. / ( h->N-1); 
	for(int i = 0 ; i < h->N; i++){
		reconst_signal[i] = Chebyshev_fn(P, &param[1], -1 + i *dn );
	}
	
	double ll = 0;
	for (int i = 0 ; i<h->N; i++){
		ll -= pow_int((h->data[i] - reconst_signal[i]),2) / ( 2 * param[0]*param[0]);
	}
	ll-= (h->N / 2)*std::log(2 * M_PI * param[0]*param[0]);
	
	return ll;
}
double transdimensional_prior(
	double *param, 
	int *status, 
	int model_status, 
	mcmc_data_interface *interface, 
	void * parameters)
{
	double a = -std::numeric_limits<double>::infinity();
	double prior=1;
	if (param[0] < .01|| param[0] > 10){ return a;}
	prior*= 1. /(10-.01); 	
	for(int i = 1 ; i<interface->max_dim; i++){
		if(status[i] !=0){
			if (param[i] < -10|| param[i] > 10){ return a;}
			prior*= 1. /(20); 	
		}
	}
	return std::log(prior);
}
double transdimensional_likelihood_fixed(
	double *param, 
	mcmc_data_interface *interface, 
	void * parameters)
{
	int temp_status[interface->max_dim];
	for(int i = 0 ; i<interface->max_dim; i++){
		temp_status[i] = 1;
	}
	return transdimensional_likelihood(param, temp_status, 0, interface, parameters);

}
double transdimensional_prior_fixed(
	double *param, 
	mcmc_data_interface *interface, 
	void * parameters)
{
	int temp_status[interface->max_dim];
	for(int i = 0 ; i<interface->max_dim; i++){
		temp_status[i] = 1;
	}
	return transdimensional_prior(param, temp_status, 0, interface, parameters);
}
void transdimensional_RJprop(
	double *current_param, 
	double *prop_param, 
	int *current_status, 
	int *prop_status, 
	int *current_model_status, 
	int *prop_model_status, 
	double *MH_corrections,
	mcmc_data_interface *interface, 
	void * parameters)
{
	trans_helper *h = (trans_helper *)parameters;
	int P = 0;
	for(int i = 0 ; i<interface->max_dim; i++){
		prop_param[i] = current_param[i] ; 
		prop_status[i] = current_status[i] ; 
		P += current_status[i];
	}
	P -=1;//For sigma

	double alpha = gsl_rng_uniform(h->r);
	
	//Kill 
	if(alpha < .5 && P >1){
		prop_status[P] = 0;
		prop_param[P] = 0;

	}
	//Create 
	else if(alpha >= .5 && P <(interface->max_dim-1)){
		prop_status[P+1] = 1;
		//prop_param[P+1] = gsl_ran_gaussian(h->r,10);
		prop_param[P+1] = gsl_rng_uniform(h->r)*20-10;

	}
	return ;
}
int ensemble_size(int argc, char *argv[])
{
	std::string data_file("full_data_transdimensional_5_5_1_100.csv");
	int N = 100;
	double dt = 1;

	double *data=new double[N] ;
	read_file("data/"+data_file, data);

	//######################################################################
	//######################################################################
	int max_thermo = 10;
	int chain_N = 1000;
	int max_dim = 10;
	int min_dim = 2;
	int samples = 50000;
	int swp_freq = 5;
	int t0 = 10000;
	int nu = 100;
	int max_chunk_size = 5000;
	int nested_models=0;
	double **output = allocate_2D_array(samples, max_dim);
	int **status = allocate_2D_array_int(samples, max_dim);
	int *model_status = NULL;
	double initial_position[max_dim];	
	int initial_status[max_dim];	
	int initial_model_status =  0;	
	double seeding_var[max_dim];	
	double **ensemble_initial_pos = NULL;
	int **ensemble_initial_status = NULL;
	int *ensemble_initial_model_status = NULL;
	double chain_temps[chain_N];
	int initnum = 5;
	initial_position[0] = 1;	
	initial_status[0] = 1;	
	seeding_var[0] = 1;	
	for(int i = 1 ; i<initnum+1; i++){
		initial_position[i] = 1;	
		initial_status[i] = 1;	
		seeding_var[i] = 1;	
	}
	for(int i = initnum+1; i<max_dim; i++){
		initial_position[i] = 0;	
		initial_status[i] = 0;	
		seeding_var[i] = 1;	
	}
	

	trans_helper **helpers = new trans_helper*[chain_N];
	const gsl_rng_type *T;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	for(int i = 0 ; i<chain_N; i++){
		helpers[i] = new trans_helper;
		helpers[i]->dt = dt;
		helpers[i]->N = N;
		helpers[i]->data = data;
		helpers[i]->r = gsl_rng_alloc(T);
	}
	std::string stat_file = "data/stat_trans.txt";
	std::string chain_file = "data/output_trans.hdf5";
	std::string likelihood_file = "data/likelihood_trans.txt";
	std::string checkpoint_file = "data/checkpoint_trans.csv";
		
	mcmc_sampler_output mcmc_out(chain_N,max_dim, nested_models) ; 

	mcmc_out.RJ= true;
	RJPTMCMC_MH_dynamic_PT_alloc_comprehensive(&mcmc_out, output, status, model_status,nested_models, max_dim, min_dim , samples, chain_N, max_thermo, initial_position, initial_status, initial_model_status, seeding_var, ensemble_initial_pos, ensemble_initial_status, ensemble_initial_model_status, chain_temps, swp_freq, t0,nu, max_chunk_size, "double", transdimensional_prior, transdimensional_likelihood, NULL, transdimensional_RJprop, (void**)helpers, 10, true, true, true, stat_file, chain_file, likelihood_file, checkpoint_file);
	mcmc_out.~mcmc_sampler_output();
	
	deallocate_2D_array(output,samples,max_dim);
	deallocate_2D_array(status,samples,max_dim);
	//######################################################################
	//######################################################################
	max_dim = 5;
	min_dim = 5;
	initnum = 4;


	initial_position[0] = 1;	
	initial_status[0] = 1;	
	seeding_var[0] = 1;	
	for(int i = 1 ; i<initnum+1; i++){
		initial_position[i] = 1;	
		initial_status[i] = 1;	
		seeding_var[i] = 1;	
	}
	//######################################################################
	//######################################################################
	samples = 300;
	chain_N = 1000;
	max_chunk_size = 1e5;
	output = allocate_2D_array(samples, max_dim);
	std::string stat_file_2 = "data/stat_ensemble_size_testing_"+std::to_string(chain_N)+".txt";
	std::string chain_file_2 = "data/output_ensemble_size_testing_"+std::to_string(chain_N)+".hdf5";
	std::string likelihood_file_2 = "data/likelihood_ensemble_size_testing_"+std::to_string(chain_N)+".txt";
	std::string checkpoint_file_2 = "data/checkpoint_ensemble_size_testing_"+std::to_string(chain_N)+".txt";
	mcmc_sampler_output mcmc_out_fixed_1(chain_N,max_dim, 0) ; 
	PTMCMC_MH_dynamic_PT_alloc_uncorrelated(&mcmc_out_fixed_1, output,max_dim, samples,chain_N, max_thermo, initial_position, seeding_var, (double **)NULL, chain_temps, swp_freq, t0,nu, max_chunk_size, "double",transdimensional_prior_fixed, transdimensional_likelihood_fixed, NULL, (void **)helpers, 10, true, true, stat_file_2, chain_file_2, likelihood_file_2, checkpoint_file_2);

	deallocate_2D_array(output,samples,max_dim);
	mcmc_out_fixed_1.~mcmc_sampler_output();
	for(int i = 0 ; i<chain_N; i++){
		gsl_rng_free(helpers[i]->r);
		delete helpers[i];
	}
	delete [] helpers;
	delete [] data;
	return 0;
}

int validate_evidence(int argc, char *argv[])
{
	std::string data_file("full_data_transdimensional_5_5_1_100.csv");
	//std::string data_file("full_data_transdimensional_5_6_1_100.csv");
	int N = 100;
	double dt = 1;

	double *data=new double[N] ;
	read_file("data/"+data_file, data);


	//######################################################################
	//######################################################################
	int max_thermo = 50;
	int chain_N = 200;
	int max_dim = 10;
	int min_dim = 2;
	int samples = 1000000;
	int swp_freq = 5;
	int t0 = 10000;
	int nu = 100;
	int max_chunk_size = 500000;
	int nested_models=0;
	double **output = allocate_2D_array(samples, max_dim);
	int **status = allocate_2D_array_int(samples, max_dim);
	int *model_status = NULL;
	double initial_position[max_dim];	
	int initial_status[max_dim];	
	int initial_model_status=0;	
	double seeding_var[max_dim];	
	double **ensemble_initial_pos = NULL;
	int **ensemble_initial_status = NULL;
	int *ensemble_initial_model_status = NULL;
	double chain_temps[chain_N];
	int initnum = 5;
	initial_position[0] = 1;	
	initial_status[0] = 1;	
	seeding_var[0] = 1;	
	for(int i = 1 ; i<initnum+1; i++){
		initial_position[i] = 1;	
		initial_status[i] = 1;	
		seeding_var[i] = 1;	
	}
	for(int i = initnum+1; i<max_dim; i++){
		initial_position[i] = 0;	
		initial_status[i] = 0;	
		seeding_var[i] = 1;	
	}
	

	trans_helper **helpers = new trans_helper*[chain_N];
	const gsl_rng_type *T;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	for(int i = 0 ; i<chain_N; i++){
		helpers[i] = new trans_helper;
		helpers[i]->dt = dt;
		helpers[i]->N = N;
		helpers[i]->data = data;
		helpers[i]->r = gsl_rng_alloc(T);
	}
	std::string stat_file = "data/stat_trans.txt";
	std::string chain_file = "data/output_trans.hdf5";
	std::string likelihood_file = "data/likelihood_trans.txt";
	std::string checkpoint_file = "data/checkpoint_trans.csv";
		
	mcmc_sampler_output mcmc_out(chain_N,max_dim, nested_models) ; 

	mcmc_out.RJ= true;
	RJPTMCMC_MH_dynamic_PT_alloc_comprehensive(&mcmc_out, output, status, model_status,nested_models, max_dim, min_dim , samples, chain_N, max_thermo, initial_position, initial_status, initial_model_status, seeding_var, ensemble_initial_pos, ensemble_initial_status, ensemble_initial_model_status, chain_temps, swp_freq, t0,nu, max_chunk_size, "double", transdimensional_prior, transdimensional_likelihood, NULL, transdimensional_RJprop, (void**)helpers, 10, true, true, true, stat_file, chain_file, likelihood_file, checkpoint_file);
	mcmc_out.~mcmc_sampler_output();
	
	deallocate_2D_array(output,samples,max_dim);
	deallocate_2D_array(status,samples,max_dim);
	//######################################################################
	//######################################################################
	max_dim = 5;
	min_dim = 5;
	initnum = 4;
	samples = 6000;

	output = allocate_2D_array(samples, max_dim);

	initial_position[0] = 1;	
	initial_status[0] = 1;	
	seeding_var[0] = 1;	
	std::string stat_file_2 = "data/stat_trans_fixed_dim1.txt";
	std::string chain_file_2 = "data/output_trans_fixed_dim1.hdf5";
	std::string likelihood_file_2 = "data/likelihood_trans_fixed_dim1.txt";
	std::string checkpoint_file_2 = "data/checkpoint_trans_fixed_dim1.txt";
	for(int i = 1 ; i<initnum+1; i++){
		initial_position[i] = 1;	
		initial_status[i] = 1;	
		seeding_var[i] = 1;	
	}
	mcmc_sampler_output mcmc_out_fixed_1(chain_N,max_dim, 0) ; 
	PTMCMC_MH_dynamic_PT_alloc_uncorrelated(&mcmc_out_fixed_1, output,max_dim, samples,chain_N, max_thermo, initial_position, seeding_var, (double **)NULL, chain_temps, swp_freq, t0,nu, max_chunk_size, "double",transdimensional_prior_fixed, transdimensional_likelihood_fixed, NULL, (void **)helpers, 10, true, true, stat_file_2, chain_file_2, likelihood_file_2, checkpoint_file_2);

	deallocate_2D_array(output,samples,max_dim);
	mcmc_out_fixed_1.~mcmc_sampler_output();
	//######################################################################
	//######################################################################
	//######################################################################
	//######################################################################
	max_dim = 6;
	min_dim = 6;
	initnum = 5;
	samples = 6000;

	output = allocate_2D_array(samples, max_dim);

	initial_position[0] = 1;	
	initial_status[0] = 1;	
	seeding_var[0] = 1;	
	std::string stat_file_3 = "data/stat_trans_fixed_dim2.txt";
	std::string chain_file_3 = "data/output_trans_fixed_dim2.hdf5";
	std::string likelihood_file_3 = "data/likelihood_trans_fixed_dim2.txt";
	std::string checkpoint_file_3 = "data/checkpoint_trans_fixed_dim2.txt";
	for(int i = 1 ; i<initnum+1; i++){
		initial_position[i] = 1;	
		initial_status[i] = 1;	
		seeding_var[i] = 1;	
	}
	mcmc_sampler_output mcmc_out_fixed_2(chain_N,max_dim, 0) ; 
	PTMCMC_MH_dynamic_PT_alloc_uncorrelated(&mcmc_out_fixed_2, output,max_dim, samples,chain_N, max_thermo, initial_position, seeding_var, (double **)NULL, chain_temps, swp_freq, t0,nu, max_chunk_size, "double",transdimensional_prior_fixed, transdimensional_likelihood_fixed, NULL, (void **)helpers, 10, true, true, stat_file_3, chain_file_3, likelihood_file_3, checkpoint_file_3);

	deallocate_2D_array(output,samples,max_dim);
	mcmc_out_fixed_2.~mcmc_sampler_output();
	//######################################################################
	//######################################################################
	max_dim = 4;
	min_dim = 4;
	initnum = 3;
	samples = 6000;

	output = allocate_2D_array(samples, max_dim);

	initial_position[0] = 1;	
	initial_status[0] = 1;	
	seeding_var[0] = 1;	
	std::string stat_file_4 = "data/stat_trans_fixed_dim3.txt";
	std::string chain_file_4 = "data/output_trans_fixed_dim3.hdf5";
	std::string likelihood_file_4 = "data/likelihood_trans_fixed_dim3.txt";
	std::string checkpoint_file_4 = "data/checkpoint_trans_fixed_dim3.txt";
	for(int i = 1 ; i<initnum+1; i++){
		initial_position[i] = 1;	
		initial_status[i] = 1;	
		seeding_var[i] = 1;	
	}
	mcmc_sampler_output mcmc_out_fixed_3(chain_N,max_dim, 0) ; 
	PTMCMC_MH_dynamic_PT_alloc_uncorrelated(&mcmc_out_fixed_3, output,max_dim, samples,chain_N, max_thermo, initial_position, seeding_var, (double **)NULL, chain_temps, swp_freq, t0,nu, max_chunk_size, "double",transdimensional_prior_fixed, transdimensional_likelihood_fixed, NULL, (void **)helpers, 10, true, true, stat_file_4, chain_file_4, likelihood_file_4, checkpoint_file_4);

	deallocate_2D_array(output,samples,max_dim);
	mcmc_out_fixed_3.~mcmc_sampler_output();
	//######################################################################
	//######################################################################
	max_dim = 3;
	min_dim = 3;
	initnum = 2;
	samples = 12000;

	output = allocate_2D_array(samples, max_dim);

	initial_position[0] = 1;	
	initial_status[0] = 1;	
	seeding_var[0] = 1;	
	std::string stat_file_5 = "data/stat_trans_fixed_dim4.txt";
	std::string chain_file_5 = "data/output_trans_fixed_dim4.hdf5";
	std::string likelihood_file_5 = "data/likelihood_trans_fixed_dim4.txt";
	std::string checkpoint_file_5 = "data/checkpoint_trans_fixed_dim4.txt";
	for(int i = 1 ; i<initnum+1; i++){
		initial_position[i] = 1;	
		initial_status[i] = 1;	
		seeding_var[i] = 1;	
	}
	mcmc_sampler_output mcmc_out_fixed_4(chain_N,max_dim, 0) ; 
	PTMCMC_MH_dynamic_PT_alloc_uncorrelated(&mcmc_out_fixed_4, output,max_dim, samples,chain_N, max_thermo, initial_position, seeding_var, (double **)NULL, chain_temps, swp_freq, t0,nu, max_chunk_size, "double",transdimensional_prior_fixed, transdimensional_likelihood_fixed, NULL, (void **)helpers, 10, true, true, stat_file_5, chain_file_5, likelihood_file_5, checkpoint_file_5);

	deallocate_2D_array(output,samples,max_dim);
	mcmc_out_fixed_4.~mcmc_sampler_output();
	//######################################################################
	//######################################################################

	for(int i = 0 ; i<chain_N; i++){
		gsl_rng_free(helpers[i]->r);
		delete helpers[i];
	}
	delete [] helpers;
	delete [] data;
	return 0;
}

//########################################################################
//########################################################################
int test_evidence(int argc, char *argv[])
{
	int temps_N = 100;
	double likelihood_vals[temps_N];
	double temps[temps_N];
	double beta;
	for(int i = 0 ; i<temps_N; i++){
		likelihood_vals[i] = 10 + 90. * i/99;	
		beta = (0+(98.-i)*1.00/99);
		temps[i] = 1./beta;	
		//std::cout<<temps[i]<<" "<<1./temps[i]<<" "<<likelihood_vals[i]<<std::endl;
	}
	double evidence, error;
	thermodynamic_integration(likelihood_vals, temps, temps_N, &evidence,&error);
	std::cout<<evidence<<" "<<error<<std::endl;
	return 0;
}


int test_likelihood(int argc, char *argv[])
{
	gen_params *gp = new gen_params;

	gp->mass1 = 2;
	gp->mass2 = 1.9;
	gp->Luminosity_Distance = 90;
	gp->phiRef = 0;
	gp->f_ref = 20;
	gp->tc = 0;
	gp->spin1[2] = 0;
	gp->spin2[2] = 0;
	gp->gmst = .2;
	gp->psi = .2;
	gp->RA = .2;
	gp->DEC = .2;
	gp->incl_angle = .051;
	gp->shift_phase = true;
	gp->shift_time = true;
	gp->equatorial_orientation = false;
	gp->horizon_coord = false;
	
	//#############################################################
	//SIMPSONS
	int length = 2*2048;
	double *freqs = new double[length];
	double *psd = new double[length];
	for(int i = 0 ; i<length; i++){
		freqs[i] = 5+(2048.-5.)*i / length;
	}
	populate_noise(freqs, "aLIGO_analytic",psd,length,48);
	for(int i = 0 ; i<length; i++){
		psd[i] *=psd[i];
	}
	
	std::complex<double> *response = new std::complex<double>[length];
	fourier_detector_response(freqs, length, response, "Hanford", "IMRPhenomD",gp, (double*)NULL);
	double snr = calculate_snr_internal(psd, response,freqs, length, "SIMPSONS",(double*)NULL,false);
	std::cout<<"SNR: "<<snr<<std::endl;
	
	clock_t start = clock();
	double ll = Log_Likelihood_internal(response,psd, freqs,(double*)NULL,response,length,false,"SIMPSONS");
	double t1 = (double)(clock()-start)/CLOCKS_PER_SEC;
	std::cout<<"time: "<<t1<<std::endl;
	std::cout<<"LL: "<<ll<<std::endl;
	std::cout<<"'SNR': "<<sqrt(2*ll)<<std::endl;
	delete [] response;
	delete [] freqs;
	delete [] psd;

	//#############################################################
	//GLQ
	length = 4000;
	freqs = new double[length];
	double *weights = new double[length];
	psd = new double[length];
	gauleg(log10(5),log10(2048),freqs,weights,length);
	for(int i = 0 ; i<length; i++){
		freqs[i] = pow(10,freqs[i]);
	}
	populate_noise(freqs, "aLIGO_analytic",psd,length,48);
	for(int i = 0 ; i<length; i++){
		psd[i] *=psd[i];
	}
	
	response = new std::complex<double>[length];
	fourier_detector_response(freqs, length, response, "Hanford", "IMRPhenomD",gp, (double*)NULL);
	snr = calculate_snr_internal(psd, response,freqs, length, "GAUSSLEG",weights,true);
	std::cout<<"SNR: "<<snr<<std::endl;
	
	start = clock();
	ll = Log_Likelihood_internal(response,psd, freqs,weights,response,length,true,"GAUSSLEG");
	double t2 = (double)(clock()-start)/CLOCKS_PER_SEC;
	std::cout<<"time: "<<t2<<std::endl;
	std::cout<<"Speedup: "<<t1/t2<<std::endl;
	delete [] response;
	delete [] freqs;
	delete [] psd;

		
	std::cout<<"LL: "<<ll<<std::endl;
	std::cout<<"'SNR': "<<sqrt(2*ll)<<std::endl;
	delete gp;
}
double log_prior_multi_gaussian(double *param, mcmc_data_interface *interface, void *parameters)
{
	double a = -std::numeric_limits<double>::infinity();
	for(int i = 0 ; i<interface->max_dim; i++){
		if( fabs(param[i]) > 10){ return a;}
	}
	//double return_val=0;
	//for(int i = 0 ; i<interface->max_dim; i++){
	//	for(int j = 0 ; j<interface->max_dim; j++){
	//		return_val-= (multi_gaussian_prior_mean[i]-param[i])*(multi_gaussian_prior_mean[j]-param[j])/2*multi_gaussian_prior_fisher[i][j];
	//	}
	//}
	//return return_val/multi_gaussian_scale;
	return 0;
}

double log_like_multi_gaussian(double *param, mcmc_data_interface *interface, void *parameters)
{
	//return 2;
	double return_val1=0;
	for(int i = 0 ; i<interface->max_dim; i++){
		for(int j = 0 ; j<interface->max_dim; j++){
			return_val1-= (multi_gaussian_like_mean[i]-param[i])*(multi_gaussian_like_mean[j]-param[j])/2*multi_gaussian_like_fisher[i][j];
		}
	}
	double return_val2=0;
	for(int i = 0 ; i<interface->max_dim; i++){
		for(int j = 0 ; j<interface->max_dim; j++){
			return_val2-= (-multi_gaussian_like_mean[i]-param[i])*(-multi_gaussian_like_mean[j]-param[j])/2*multi_gaussian_like_fisher[i][j];
		}
	}
	//double long_target = -pow(param[0]-30,2)/2 ;
	//return log(exp(return_val1) + exp(return_val2))/multi_gaussian_scale;
	return (return_val1 + log(1+exp(return_val2-return_val1)))/multi_gaussian_scale;
	//return (log(exp(return_val1) + exp(return_val2)+exp(long_target)))/multi_gaussian_scale;
	//return 1;

}
void fisher_multi_gaussian(double *param, double **fisher, mcmc_data_interface *interface, void *parameters)
{
	double return_val=0;
	for(int i = 0 ; i<interface->max_dim; i++){
		for(int j = 0 ; j<interface->max_dim; j++){
			fisher[i][j] = multi_gaussian_like_fisher[i][j];
		}
	}

}

int multiple_continue(int argc, char *argv[])
{
	int total_samples = 2e4;
	int dimension = 6;	
	
	//int chain_N = 100;
	//int max_chain_N_thermo_ensemble = 5;	
	int chain_N = 100;
	int max_chain_N_thermo_ensemble = 10;	
	int N_steps = total_samples/((double)chain_N/ max_chain_N_thermo_ensemble);
	std::cout<<"Number of steps per chain: "<<N_steps<<std::endl;
	double *seeding_var=NULL;
	double chain_temps[chain_N];
	int swp_freq = 5;
	int t0 = 3000;
	int nu = 100;
	int max_chunk_size = 100000;	
	std::string chain_distribution_scheme = "double";
	void ** user_parameters = NULL;
	int numthreads = 8;
	bool pool = true;
	bool show_prog = true;
	std::string modifier="";
	std::string stat_file = "data/gaussian_stat_0"+modifier+".txt";
	std::string chain_file = "data/gaussian_output_0"+modifier+".hdf5";
	std::string likelihood_file = "data/gaussian_likelihood_0"+modifier+".csv";
	std::string check_file = "data/gaussian_checkpoint_0"+modifier+".csv"	;
	mcmc_sampler_output sampler_output(chain_N, dimension) ;
	double **output = allocate_2D_array(N_steps,dimension);

	multi_gaussian_like_cov = new double*[dimension];
	multi_gaussian_prior_cov = new double*[dimension];
	multi_gaussian_like_fisher = new double*[dimension];
	multi_gaussian_prior_fisher = new double*[dimension];
	multi_gaussian_like_mean = new double[dimension];
	multi_gaussian_prior_mean = new double[dimension];
	multi_gaussian_like_mean[0]=4;
	multi_gaussian_like_mean[1]=2;
	multi_gaussian_like_mean[2]=2;
	multi_gaussian_like_mean[3]=1;
	multi_gaussian_like_mean[4]=2;
	multi_gaussian_like_mean[5]=3;
	//multi_gaussian_like_mean[6]=2;
	//multi_gaussian_like_mean[7]=4;
	//multi_gaussian_like_mean[8]=3;
	//multi_gaussian_like_mean[9]=2;
	//multi_gaussian_like_mean[10]=3;
	//multi_gaussian_like_mean[11]=3;
	//multi_gaussian_like_mean[12]=2;
	//multi_gaussian_like_mean[13]=4;
	//multi_gaussian_like_mean[14]=5;
	
	multi_gaussian_prior_mean[0]=4;
	multi_gaussian_prior_mean[1]=-2;
	multi_gaussian_prior_mean[2]=2;
	multi_gaussian_prior_mean[3]=-1;
	multi_gaussian_prior_mean[4]=2;
	multi_gaussian_prior_mean[5]=3;
	//multi_gaussian_prior_mean[3]=1;
	//multi_gaussian_prior_mean[4]=2;
	//multi_gaussian_prior_mean[5]=3;
	//multi_gaussian_prior_mean[6]=2;
	//multi_gaussian_prior_mean[7]=4;
	//multi_gaussian_prior_mean[8]=3;
	//multi_gaussian_prior_mean[9]=2;
	for(int i = 0 ; i<dimension; i++){
		multi_gaussian_prior_cov[i]=new double[dimension];
		multi_gaussian_like_cov[i]=new double[dimension];
		multi_gaussian_like_fisher[i]=new double[dimension];
		multi_gaussian_prior_fisher[i]=new double[dimension];
		for(int j = 0 ; j<dimension; j++){
			if(j == i){
				//multi_gaussian_prior_cov[i][j] = 25*multi_gaussian_like_mean[j];
				multi_gaussian_prior_cov[i][j] = 5;
				//multi_gaussian_prior_cov[i][j] = 0;
				//multi_gaussian_prior_cov[i][j] = 1000;
				multi_gaussian_like_cov[i][j] = 1;
			}
			else{
				multi_gaussian_prior_cov[i][j] = 0;
				multi_gaussian_like_cov[i][j] = pow(-1, j+i)*.5;
				//multi_gaussian_like_cov[i][j] = 0;
			}
		}

	}
	gsl_LU_matrix_invert(multi_gaussian_prior_cov, multi_gaussian_prior_fisher, dimension);
	gsl_LU_matrix_invert(multi_gaussian_like_cov, multi_gaussian_like_fisher, dimension);
	for(int i = 0 ; i<dimension; i++){
		for(int j = 0 ; j<dimension; j++){
			std::cout<<multi_gaussian_like_fisher[i][j]<<" ";	
		}
		std::cout<<std::endl;
	
	}
	

	//double *init_pos = &multi_gaussian_like_mean[0];
	double *init_pos = new double[dimension];
	for(int i = 0 ;i<dimension; i++){
		init_pos[i] = -9;
		//init_pos[i] = -multi_gaussian_like_mean[i];
	}
	write_file("data/multi_gaussian_prior_cov.csv",multi_gaussian_prior_cov,dimension,dimension);
	write_file("data/multi_gaussian_like_cov.csv",multi_gaussian_like_cov,dimension,dimension);
	write_file("data/multi_gaussian_like_mean.csv",multi_gaussian_like_mean,dimension);
	write_file("data/multi_gaussian_prior_mean.csv",multi_gaussian_prior_mean,dimension);
	
	//PTMCMC_MH_dynamic_PT_alloc_uncorrelated(&sampler_output,output, dimension, N_steps, chain_N, max_chain_N_thermo_ensemble, init_pos, seeding_var, (double**)NULL,chain_temps, swp_freq, t0, nu, max_chunk_size, chain_distribution_scheme, log_prior_multi_gaussian, log_like_multi_gaussian, fisher_multi_gaussian, user_parameters, numthreads, pool, show_prog, stat_file, chain_file, likelihood_file, check_file);
	PTMCMC_MH_dynamic_PT_alloc_uncorrelated(&sampler_output,output, dimension, N_steps, chain_N, max_chain_N_thermo_ensemble, init_pos, seeding_var, (double**)NULL,chain_temps, swp_freq, t0, nu, max_chunk_size, chain_distribution_scheme, log_prior_multi_gaussian, log_like_multi_gaussian, NULL, user_parameters, numthreads, pool, show_prog, stat_file, chain_file, likelihood_file, check_file);
	
	for(int i = 1 ; i<10; i++){
		std::string stat_file_new = "data/gaussian_stat_"+std::to_string(i)+modifier+".txt";
		std::string chain_file_new = "data/gaussian_output_"+std::to_string(i)+modifier+".hdf5";
		std::string likelihood_file_new = "data/gaussian_likelihood_"+std::to_string(i)+modifier+".csv";
		std::string check_file_new = "data/gaussian_checkpoint_"+std::to_string(i)+modifier+".csv"	;
		std::string old_check = "data/gaussian_checkpoint_"+std::to_string(i-1)+modifier+".csv";
		
		mcmc_sampler_output s_out(chain_N,dimension);
		continue_PTMCMC_MH_dynamic_PT_alloc_uncorrelated(old_check,&s_out,output,  N_steps,  max_chain_N_thermo_ensemble,   chain_temps, swp_freq, t0, nu, max_chunk_size, chain_distribution_scheme, log_prior_multi_gaussian, log_like_multi_gaussian, fisher_multi_gaussian, user_parameters, numthreads, pool, show_prog, stat_file_new, chain_file_new, likelihood_file_new, check_file_new);
	}
	deallocate_2D_array(output, chain_N,N_steps);
	for(int i = 0 ; i<dimension; i++){
		delete [] multi_gaussian_like_cov[i];	
		delete [] multi_gaussian_prior_cov[i];	
		delete [] multi_gaussian_like_fisher[i];	
		delete [] multi_gaussian_prior_fisher[i];	
	}
	delete [] multi_gaussian_like_cov;	
	delete [] multi_gaussian_prior_cov;	
	delete [] multi_gaussian_like_fisher;	
	delete [] multi_gaussian_prior_fisher;	
	delete [] multi_gaussian_like_mean;	
	delete [] multi_gaussian_prior_mean;	
	delete [] init_pos;
	return 0;
}

double RJ_sin_tilt_range[2];
int mcmc_RJ_sin(int argc, char *argv[])
{
	RJ_sin_tilt_range[0] = 0;
	RJ_sin_tilt_range[1] = .2;
	int dim = 6;
	int N = 5000;
	//int N = 500;
	double injections[dim];
	injections[0]=200;
	//injections[0]=400;
	injections[1]=2;
	injections[2]=.1;
	injections[3]=500.;
	injections[4]=N/2.;
	injections[5]=0.05;

	write_file("data/RJ_injections_sin.csv",injections,dim);

	double *data = new double[N];
	double *data_pure = new double[N];
	gsl_rng_env_setup();
	const gsl_rng_type *T = gsl_rng_default;
	gsl_rng *rand = gsl_rng_alloc(T);
	gsl_rng_set(rand,10);
	double snr = 0;
	double alpha = injections[5];
	for(int i= 0 ; i<N; i++){
		double r = gsl_ran_gaussian(rand,1);
		double shifting_factor = (1+erf((alpha*(i-injections[4]))/root2));
		data_pure[i]=injections[0]*sin(injections[2]*(i)+injections[1])
			*exp(-.5*(injections[4]-i)*(injections[4]-i)/(injections[3]*injections[3]))/sqrt(2*injections[3]*injections[3]) *shifting_factor;
		data[i]=data_pure[i]+r;
		snr+=data[i]*data[i];
	}
	snr /= N;
	std::cout<<"SNR: "<<sqrt(snr)<<std::endl;
	
	write_file("data/RJ_sin_data_injection.csv",data,N);
	write_file("data/RJ_sin_data_pure_injection.csv",data_pure,N);
	delete [] data_pure;
	gsl_rng_free(rand);

	//int chain_N = 80;
	int chain_N = 21;
	//int max_ensemble_chain_N = 8;
	int max_ensemble_chain_N = 3;
	double initial_pos[dim];	
	double seeding_var[dim];	
	for(int i = 0 ; i<dim; i++){
		initial_pos[i] = injections[i];
		seeding_var[i] = injections[i];
	}
	int swp_freq = 5;
	std::string chain_distribution="double";
	RJ_sin_param **param = new RJ_sin_param*[chain_N];
	for(int i = 0 ; i<chain_N; i++){
		param[i]=new RJ_sin_param;
		param[i]->data = data;
		param[i]->N = N;
		param[i]->r = gsl_rng_alloc(T);
	}
	std::string stat_file = "data/stat_RJ_sin.txt";
	std::string chain_file = "data/chains_RJ_sin.hdf5";
	std::string checkpoint = "data/checkpoint_RJ_sin.csv";
	double chain_temps[chain_N];
	for(int i = 0 ; i<chain_N; i++){
		if( i %max_ensemble_chain_N == 0 ){
			chain_temps[i] = 1.;
		}
		else{
			chain_temps[i] = chain_temps[i-1]*1.5;
		}
	}
	int numthreads = 8;
	bool pool = false;
	bool show_prog = true;
	
	//############################################3
	//int N_steps = 1000;
	//mcmc_sampler_output sampler_output(chain_N,dim);
	//double **output = new double*[N_steps];
	//for(int i = 0 ; i<N_steps;i++){
	//	output[i]=new double[dim];
	//}
	//PTMCMC_MH_dynamic_PT_alloc_uncorrelated(&sampler_output,output, dim, N_steps, chain_N,max_ensemble_chain_N, initial_pos, seeding_var, chain_temps, swp_freq, t0,nu, max_chunk_size, chain_distribution, RJ_sin_prior, RJ_sin_logL, NULL, (void **)param, numthreads, pool, show_prog, stat_file, chain_file, "",checkpoint);
	//sampler_output.calc_ac_vals(true);
	//sampler_output.count_indep_samples(true);
	//sampler_output.create_data_dump(true,false,"data/chains_RJ_sin.hdf5");
	//sampler_output.create_data_dump(false,false,"data/chains_sin_full.hdf5");
	//for(int i = 0 ; i<N_steps;i++){
	//	delete [] output[i];
	//}
	//delete [] output;
	//############################################3
	int N_steps = 5*100;
	int max_dim = dim;
	int min_dim = dim-1;
	int initial_status[max_dim];
	for(int i = 0 ; i<max_dim; i++){
		initial_status[i]=1;
	}

	//##################################################	
	//###############################################
	//double ***output = new double**[chain_N];
	//int ***status = new int**[chain_N];
	//for(int i = 0 ; i<chain_N;i++){
	//	output[i]=new double*[N_steps];
	//	status[i]=new int*[N_steps];
	//	for(int j = 0 ; j<N_steps;j ++){
	//		output[i][j]=new double[dim];
	//		status[i][j]=new int[dim];
	//	}
	//}
	//RJPTMCMC_MH(output,status, max_dim,min_dim, N_steps, chain_N, initial_pos,initial_status, seeding_var, chain_temps, swp_freq,  RJ_sin_prior, RJ_sin_logL, RJ_sin_fish, RJ_sin_proposal, (void **)param, numthreads, pool, show_prog, stat_file, chain_file, "","",checkpoint);
	//for(int i = 0 ; i<chain_N;i++){
	//	for(int j = 0 ; j<N_steps;j++){
	//		delete [] output[i][j];
	//		delete [] status[i][j];
	//	}
	//	delete [] output[i];
	//	delete [] status[i];
	//}
	//delete [] output;
	//delete [] status;
	//##################################################	
	//###############################################
	int t0 = 50;
	int nu = 10;
	int max_chunksize = 500;
	bool update_RJ_width = true;
	double **output = new double*[N_steps];
	int **status = new int*[N_steps];
	int nested_model_number =3;
	int *model_status = NULL;
	int initial_model_status = 0;
	if(nested_model_number !=0){
		model_status = new int[N_steps];
	}
	for(int j = 0 ; j<N_steps;j ++){
		output[j]=new double[dim];
		status[j]=new int[dim];
	}

	mcmc_sampler_output sampler_output(chain_N,max_dim,nested_model_number);
	sampler_output.RJ = true;
	sampler_output.nested_model_number = nested_model_number;
	RJPTMCMC_MH_dynamic_PT_alloc_comprehensive_internal(&sampler_output,output,status, model_status,nested_model_number,max_dim,min_dim, N_steps, chain_N, max_ensemble_chain_N,initial_pos,initial_status, initial_model_status,seeding_var,(double**)NULL,(int**)NULL,(int*)NULL, chain_temps, swp_freq,  t0, nu, max_chunksize, chain_distribution,RJ_sin_prior, RJ_sin_logL, RJ_sin_fish, RJ_sin_proposal, (void **)param, numthreads, pool, show_prog, update_RJ_width, stat_file, chain_file, "",checkpoint);
	if(model_status){
		delete [] model_status;
	}
	for(int j = 0 ; j<N_steps;j++){
		delete [] output[j];
		delete [] status[j];
	}
	delete [] output;
	delete [] status;
	//sampler_output.~mcmc_sampler_output();

	//##################################################	
	//###############################################
	
	//Cleanup
	for(int i = 0 ; i<chain_N; i++){
		gsl_rng_free(param[i]->r);
		delete param[i];	
	}
	delete [] param;
	delete [] data;
	return 0;
}

void RJ_sin_proposal(double *current_params, double *proposed_params, int *current_status, int *proposed_status,int *current_model_status, int *proposed_model_status, double *MH_corrections,mcmc_data_interface *interface, void *parameters)
{
	int current_dim = 0;
	for(int i = 0 ; i<interface->max_dim; i++){
		current_dim+=current_status[i];
	}
	RJ_sin_param *param = (RJ_sin_param *) parameters;
	//double alpha = gsl_rng_uniform(param->r);
	double std = 1;
	double mean = .10;
	double a = 1;
	double b = 1;
	//if(alpha < .5 && current_dim == 5){
	if(current_dim == 5){
		for(int i = 0 ; i<interface->max_dim; i++){
			proposed_params[i]=current_params[i];
			proposed_status[i]=current_status[i];
		}
		//proposed_params[5] =std::fabs(gsl_ran_gaussian(param->r,std)) ;
		//proposed_params[5] =gsl_ran_gaussian(param->r,std)+mean ;
		//proposed_params[5] =gsl_ran_gamma(param->r,a,b);
		proposed_params[5] =gsl_rng_uniform(param->r)*(RJ_sin_tilt_range[1]-RJ_sin_tilt_range[0])+RJ_sin_tilt_range[0] ;
		proposed_status[5] = 1;
		//double gaussian_cor = 2*1./sqrt(2*M_PI*std*std) * exp(-0.5*pow_int(proposed_params[5],2)/(std*std));
		//double log_cor = log(2*1./sqrt(2*M_PI*std*std)) +( -0.5*pow_int(proposed_params[5],2)/(std*std));
		//double log_cor = log(1./sqrt(2*M_PI*std*std)) +( -0.5*pow_int(mean-proposed_params[5],2)/(std*std));
		//double log_cor = log(1./sqrt(2*M_PI*std*std)) +( -0.5*pow_int(mean-proposed_params[5],2)/(std*std));
		//double log_cor = log(pow(proposed_params[5],a-1.)/(gsl_sf_gamma(a)*pow(b,a) ) )-( proposed_params[5]/b);
		//double log_cor = log(gsl_ran_gamma_pdf(proposed_params[5],a,b));
		double log_cor = log(1. / (RJ_sin_tilt_range[1]-RJ_sin_tilt_range[0]));
		*(MH_corrections) = -1*(log_cor);
	}
	//else if (alpha >= .5 && current_dim ==6){
	else if (current_dim ==6){
		for(int i = 0 ; i<interface->max_dim; i++){
			proposed_params[i]=current_params[i];
			proposed_status[i]=current_status[i];
		}
		proposed_params[5] =0 ;
		proposed_status[5] = 0;
		//double log_cor = log(2*1./sqrt(2*M_PI*std*std)) +( -0.5*pow_int(proposed_params[5],2)/(std*std));
		//double log_cor = log(1./sqrt(2*M_PI*std*std)) +( -0.5*pow_int(mean-current_params[5],2)/(std*std));
		//double log_cor = log(pow(current_params[5],a-1.)/(gsl_sf_gamma(a)*pow(b,a) ) )-( current_params[5]/b);
		//double log_cor = log(gsl_ran_gamma_pdf(current_params[5],a,b));
		double log_cor = log(1. / (RJ_sin_tilt_range[1]-RJ_sin_tilt_range[0]));
		*(MH_corrections) = 1*(log_cor);

	}
	else{
		for(int i = 0 ; i<interface->max_dim; i++){
			proposed_params[i]=current_params[i];
			proposed_status[i]=current_status[i];
		}
	}
	*(proposed_model_status) = *(current_model_status);
	return ;
}

double RJ_sin_logL(double *params, int *status,int model_status,mcmc_data_interface *interface, void *parameters)
//double RJ_sin_logL(double *params, mcmc_data_interface *interface, void *parameters)
{
	return 0;
	//for(int i = 0 ; i<2; i++){
	//	std::cout<<model_status[i]<<" ";
	//}
	//std::cout<<std::endl;
	int dim = 0;
	for(int i = 0 ; i<interface->max_dim; i++){
		dim+=status[i];
	}
	if(dim == 5){
		RJ_sin_param *param = (RJ_sin_param *)parameters;
		double A = params[0];
		double phi0 = params[1];
		double f = params[2];
		double sigma = params[3];
		double N0 = params[4];
		//double t = params[3];
		double L = 0;
		for(int i= 0 ; i<param->N; i ++){
			L -= .5 * pow_int( A*sin(f*(i) +phi0)*exp(-.5*(N0-i)*(N0-i)/(sigma*sigma))/sqrt(2*sigma*sigma) - param->data[i] ,2);
		}
		return L;
	}
	else{

		RJ_sin_param *param = (RJ_sin_param *)parameters;
		double A = params[0];
		double phi0 = params[1];
		double f = params[2];
		double sigma = params[3];
		double N0 = params[4];
		double alpha = params[5];
		double L = 0;
		for(int i= 0 ; i<param->N; i ++){
			double shifting_factor = (1+erf((alpha*(i-N0))/root2));
			L -= .5 * pow_int( A*sin(f*(i) +phi0)*exp(-.5*(N0-i)*(N0-i)/(sigma*sigma))/sqrt(2*sigma*sigma)*shifting_factor - param->data[i] ,2);
		}
		return L;
	}
}
double RJ_sin_prior(double *params,  int *status,int model_status,mcmc_data_interface *interface, void *parameters)
//double RJ_sin_prior(double *params,  mcmc_data_interface *interface, void *parameters)
{
	//int dim = interface->max_dim;
	
	int dim = 0;
	for(int i = 0 ; i<interface->max_dim; i++){
		dim+=status[i];
	}
	double p = 0;
	double a = -std::numeric_limits<double>::infinity();
	RJ_sin_param *param = (RJ_sin_param *)parameters;
	if(params[0]<0 || params[0]>5000){return a;}
	//if(params[0]<0 || params[0]>10000){return a;}
	if(params[1]<0 || params[1]>2*M_PI){return a;}
	if(params[2]<0 || params[2]>.5){return a;}
	if(params[3]<0 || params[3]>1000){return a;}
	if(params[4]<0 || params[4]>param->N){return a;}
	p += log(1./(5000. - 0));
	p += log(1./(2.*M_PI - 0));
	p += log(1./(0.5 - 0));
	p += log(1./(1000 - 0));
	p += log(1./(param->N - 0));
	if(dim==6){
		if(params[5]<RJ_sin_tilt_range[0] || params[5]>RJ_sin_tilt_range[1]){return a;}
		p += log(1./(RJ_sin_tilt_range[1] -RJ_sin_tilt_range[0]));
	}
	//if(params[3]<0 || params[3]>5){return a;}
	return p;
}
void RJ_sin_fish(double *c, int *status,int model_status,double **fisher,mcmc_data_interface *interface, void *parameters)
{
	int dim = interface->min_dim;
	double epsilon = 1e-3;
	int current_dim = 0;
	for(int i = 0 ; i<interface->max_dim; i++){
		current_dim+=status[i];
	}
	double temp[current_dim];
	for (int i = 0 ; i<current_dim ; i++){
		temp[i]=c[i];	
	}
	for (int i = 0 ; i<dim ; i++){
		for (int j = 0 ; j<i ; j++){
			//temp[i] = c[i]+epsilon;
			//temp[j] = c[j]+epsilon;
			//double plusplus = logL(temp,dim,chain_id, parameters);
			//temp[i] = c[i]-epsilon;
			//temp[j] = c[j]-epsilon;
			//double minusminus = logL(temp,dim,chain_id, parameters);
			//temp[i] = c[i]+epsilon;
			//temp[j] = c[j]-epsilon;
			//double plusminus = logL(temp,dim,chain_id, parameters);
			//temp[i] = c[i]-epsilon;
			//temp[j] = c[j]+epsilon;
			//double minusplus = logL(temp,dim,chain_id, parameters);
			//temp[i] = c[i];
			//temp[j] = c[j];
			//fisher[i][j]= -(plusplus -plusminus - minusplus + minusminus)/(4.*epsilon*epsilon);
			fisher[i][j]= 0;
		}
		temp[i]= c[i]+epsilon;
		double plus = RJ_sin_logL(temp,status,model_status,interface,parameters);	
		temp[i]= c[i]-epsilon;
		double minus = RJ_sin_logL(temp,status,model_status,interface,parameters);	
		temp[i]=c[i];
		double eval = RJ_sin_logL(temp,status,model_status,interface,parameters);	
		fisher[i][i]= -(plus -2*eval +minus)/(epsilon*epsilon);
		fisher[i][i]= 1;
	}
	for (int i = 0 ; i<dim ; i++){
		for (int j = i+1 ; j<dim ; j++){
			fisher[i][j] = fisher[j][i];	
		}
	}
	
	return;
}

int mcmc_output_class(int argc, char *argv[])
{
	bool RJ = true;
	int nested_models = 0;
	if(RJ){nested_models = 1;}
	int chain_N = 3;
	int dim = 3;
	mcmc_sampler_output output(chain_N,dim,nested_models);
	output.RJ = RJ;
	double chain_temps[chain_N] = {1.,2.,1.};
	output.populate_chain_temperatures(chain_temps);
	for(int i = 0 ; i<chain_N; i++){
		std::cout<<chain_temps[i]<<std::endl;
	}
	std::cout<<"Cold chains: "<<output.cold_chain_number<<std::endl;
	std::cout<<"Cold chains id 1: "<<output.cold_chain_ids[0]<<std::endl;
	std::cout<<"Cold chains id 2: "<<output.cold_chain_ids[1]<<std::endl;
	double ***init_output = new double**[chain_N];
	double ***init_LL_LP = new double**[chain_N];
	int **init_model_status =NULL;
	int ***init_status =NULL;
	if(RJ){
		init_model_status = new int*[chain_N];
		init_status = new int**[chain_N];
	}
	
	int steps = 5;
	int init_positions[chain_N] = {5,3,5};
	for(int i = 0 ; i<chain_N; i++){
		if(RJ){
			init_model_status[i]  = new int[steps];
			init_status[i]  = new int*[steps];
		}
		init_output[i] = new double*[steps];
		init_LL_LP[i] = new double*[steps];
		for(int j =0 ; j<steps; j++){
			if(RJ){
				init_model_status[i][j] = 3;
				init_status[i][j] = new int[dim];
			}
			init_output[i][j] = new double[dim];
			init_LL_LP[i][j] = new double[2];
			if(j < init_positions[i]){
				init_LL_LP[i][j][0] = 10;
				init_LL_LP[i][j][1] = 1;
				for(int k =0 ; k<dim; k++){
					init_output[i][j][k] = i+j+k;
					if(RJ){
						init_status[i][j][k] = 1;	
					}
				}
			}
			else{
				init_LL_LP[i][j][0] = 0;
				init_LL_LP[i][j][1] = 0;
				for(int k =0 ; k<dim; k++){
					init_output[i][j][k] = 0;
					if(RJ){
						init_status[i][j][k] = 0;	
					}
				}
			}
		}
	}
	output.populate_initial_output(init_output,init_status,init_model_status, init_LL_LP,init_positions);
	std::cout<<"Initial output"<<std::endl;
	for(int i = 0 ; i<chain_N; i++){
		for(int j =0 ; j<init_positions[i]; j++){
			for(int k =0 ; k<dim; k++){
				std::cout<<output.output[i][j][k]<<" ";
			}
			std::cout<<std::endl;
		}
		std::cout<<std::endl;
	}
	output.trim_lengths[0]=2;
	output.trim_lengths[1]=2;
	output.trim_lengths[2]=2;


	std::cout<<"Creating data dump"<<std::endl;
	output.create_data_dump(false,true, "./data/testT.hdf5");
	output.create_data_dump(false,false, "./data/test.hdf5");
	//exit(1);

	double ***app_output = new double**[chain_N];
	double ***app_LL_LP = new double**[chain_N];
	int ***app_status = NULL;
	int **app_model_status = NULL;
	if(RJ){
		app_status = new int**[chain_N];
		app_model_status = new int*[chain_N];
	}
	int app_steps = 8;
	int app_positions[chain_N] = {8,5,8};
	for(int i = 0 ; i<chain_N; i++){
		if(RJ){
			app_status[i] = new int*[app_steps];
			app_model_status[i] = new int[app_steps];
		}
		app_output[i] = new double*[app_steps];
		app_LL_LP[i] = new double*[app_steps];
		for(int j =0 ; j<app_steps; j++){
			if(RJ){
				app_status[i][j] = new int[dim];
				app_model_status[i][j] = 4;
			}
			app_output[i][j] = new double[dim];
			app_LL_LP[i][j] = new double[2];
			if(j < app_positions[i]){
				app_LL_LP[i][j][0] = 10*(i+j);
				app_LL_LP[i][j][1] = 11*(i+j);
				for(int k =0 ; k<dim; k++){
					app_output[i][j][k] = 10*(i+j+k);
					if(RJ){
						app_status[i][j][k] = 1;
					}
				}
			}
			else{
				app_LL_LP[i][j][0] = 0;
				app_LL_LP[i][j][1] = 0;
				for(int k =0 ; k<dim; k++){
					app_output[i][j][k] = 0;
					if(RJ){
						app_status[i][j][k] = 0;
					}
				}
			}
		}
	}
	chain_temps[1]=5;
	output.append_to_output(app_output,app_status,app_model_status, app_LL_LP,app_positions);
	output.populate_chain_temperatures(chain_temps);
	
	std::cout<<std::endl;
	std::cout<<std::endl;
	std::cout<<"Appended output"<<std::endl;
	for(int i = 0 ; i<chain_N; i++){
		for(int j =0 ; j<output.chain_lengths[i]; j++){
			for(int k =0 ; k<dim; k++){
				std::cout<<output.output[i][j][k]<<" ";
			}
			std::cout<<std::endl;
		}
		std::cout<<std::endl;
	}

	std::cout<<std::endl;
	std::cout<<std::endl;
	std::cout<<"Appending new data to data dump"<<std::endl;
	output.append_to_data_dump("./data/testT.hdf5");
	output.append_to_data_dump("./data/test.hdf5");


	output.calc_ac_vals(true);
	for(int i = 0 ; i<output.cold_chain_number; i++){
		for(int j = 0 ; j<output.dimension; j++){
			std::cout<<output.ac_vals[i][j]<<" ";
		}
		std::cout<<std::endl;
	}
	output.count_indep_samples(false);
	std::cout<<"Indep samples "<<output.indep_samples<<std::endl;


	output.write_flat_thin_output( "./data/test_flat",false,false);
	//############################################################
	//Cleanup
	//############################################################
	if(init_model_status){
		for(int i = 0 ; i<chain_N; i ++){
			delete [] init_model_status[i];
		}
		delete [] init_model_status;

	}
	if(init_status){
		for(int i = 0 ; i<chain_N; i ++){
			for(int j = 0 ; j<steps; j++){
				delete [] init_status[i][j];
			}
			delete [] init_status[i];
		}
		delete [] init_status;

	}
	for(int i = 0 ; i<chain_N; i ++){
		for(int j = 0 ; j<steps; j++){
			delete [] init_output[i][j];
			delete [] init_LL_LP[i][j];
		}
		delete [] init_output[i];
		delete [] init_LL_LP[i];
	}
	delete [] init_output;
	delete [] init_LL_LP;
	init_output = NULL;
	init_LL_LP = NULL;
	for(int i = 0 ; i<chain_N; i ++){
		for(int j = 0 ; j<app_steps; j++){
			delete [] app_output[i][j];
			delete [] app_LL_LP[i][j];
		}
		delete [] app_output[i];
		delete [] app_LL_LP[i];
	}
	delete [] app_output;
	delete [] app_LL_LP;
	app_output = NULL;
	app_LL_LP = NULL;
	if(app_model_status){
		for(int i = 0 ; i<chain_N; i ++){
			delete [] app_model_status[i];
		}
		delete [] app_model_status;

	}
	if(app_status){
		for(int i = 0 ; i<chain_N; i ++){
			for(int j = 0 ; j<steps; j++){
				delete [] app_status[i][j];
			}
			delete [] app_status[i];
		}
		delete [] app_status;

	}
	
	return 0;
}
struct rosenbock_param
{
	double a;
	double mu;
	int n1;
	int n2;
	int n;
	double *b;
	double *bounds_high;
	double *bounds_low;
};
//https://arxiv.org/pdf/1903.09556.pdf
int mcmc_rosenbock(int argc, char *argv[])
{
	double a = 1./20.;
	//int n1 = 3;
	//int n2 = 2;
	int n1 = 5;
	int n2 = 3;
	int n = (n1-1)*n2 + 1;
	double b[n];
	for(int i = 0 ; i<n ; i++){
		b[i]=100./20.;
	}
	double mu = 1.;
	double bounds_high[n];
	double bounds_low[n];
	for(int i=0 ; i<n; i++){
		bounds_high[i]=10000;
		bounds_low[i]=-10000;
	}
	double out_param[5+n];
	out_param[0]=a;
	out_param[1]=mu;
	out_param[2]=n;
	out_param[3]=n1;
	out_param[4]=n2;
	for(int i = 0 ; i<n; i++){
		out_param[5+i]=b[i];
	}
	write_file("data/rosenbock_parameters.csv",out_param,5+n);



	double initial_pos[n];
	double seeding_var[n] ;
	for(int i = 0 ; i<n ; i++){
		initial_pos[i]=pow(-1,i)*500.;
		seeding_var[i]=10.;
	}
	int N_steps = 10000;
	int chain_N= 240;
	int max_chain_N= 80;
	//double *initial_pos_ptr = initial_pos;
	int swp_freq = 5;
	//double chain_temps[chain_N] ={1,2,3,10,12};
	double chain_temps[chain_N];
	//chain_temps[0] = 1.;
	//double c = 1.8;
	//for(int i =1; i < chain_N/2;  i ++)
	//	chain_temps[i] =  chain_temps[i-1] * c;
	//chain_temps[chain_N/2] = 1;
	//for(int i =chain_N/2+1; i < chain_N;  i ++)
	//	chain_temps[i] =  chain_temps[i-1] * c;
	std::string autocorrfile = "";
	std::string chainfile = "data/mcmc_output_RB.hdf5";
	std::string statfilename = "data/mcmc_statistics_RB.txt";
	std::string checkpointfile = "data/mcmc_checkpoint_RB.csv";
	//std::string LLfile = "data/mcmc_LL_RB.csv";
	std::string LLfile = "";
	
	int numThreads = 10;
	bool pool = true;
	bool show_progress = true;
	
	double **output;
	output = allocate_2D_array(  N_steps, n );
	int t0 = 15000;
	int nu = 100;
	double corr_threshold = 0.01;
	int corr_segments = 5;
	double corr_convergence_thresh = 0.01;
	double corr_target_ac = 0.01;
	std::string chain_distribution_scheme="double";
	int max_chunk_size = 1000000;
	
	rosenbock_param **param = new rosenbock_param*[chain_N];
	for(int i = 0 ; i<chain_N; i++){
		param[i]=new rosenbock_param;
		param[i]->a = a;
		param[i]->b = b;
		param[i]->bounds_high = bounds_high;
		param[i]->bounds_low = bounds_low;
		param[i]->mu = mu;
		param[i]->n1 = n1;
		param[i]->n2 = n2;
		param[i]->n = n;
	}
	
	mcmc_sampler_output sampler_output(chain_N,n);
	PTMCMC_MH_dynamic_PT_alloc_uncorrelated(&sampler_output,output, n, N_steps, chain_N, max_chain_N,initial_pos,seeding_var,(double**)NULL,chain_temps, swp_freq, t0,nu, max_chunk_size,chain_distribution_scheme, log_rosenbock_prior, log_rosenbock,NULL,(void **)param,numThreads, pool,show_progress, statfilename,chainfile, LLfile,checkpointfile );	
	sampler_output.calc_ac_vals(true);
	for(int i = 0 ; i<sampler_output.cold_chain_number; i++){
		std::cout<<sampler_output.max_acs[i]<<std::endl;
	}
	sampler_output.count_indep_samples(true);
	std::cout<<"Indep samples "<<sampler_output.indep_samples<<std::endl;
	sampler_output.create_data_dump(true,false, chainfile);
	sampler_output.create_data_dump(true,true,"data/mcmc_output_RB_trimmed.hdf5" );
	//sampler_output.write_flat_thin_output( "./data/test.hdf5",true,true);
	deallocate_2D_array(output, N_steps, n);
	for(int i = 0 ; i<chain_N; i++){
		delete param[i];	
	}
	delete param;	
	std::cout<<"ENDED"<<std::endl;


	return 0;
}

//void fisher_rosenbock(double *c,int dim, double **fisher,  void *parameters)
//{
//	double epsilon = 1e-3;
//	double temp[dim];
//	for (int i = 0 ; i<dim ; i++){
//		temp[i]=c[i];	
//	}
//	for (int i = 0 ; i<dim ; i++){
//		for (int j = 0 ; j<i ; j++){
//			temp[i] = c[i]+epsilon;
//			temp[j] = c[j]+epsilon;
//			double plusplus = log_rosenbock(temp,dim, parameters);
//			temp[i] = c[i]-epsilon;
//			temp[j] = c[j]-epsilon;
//			double minusminus = log_rosenbock(temp,dim, parameters);
//			temp[i] = c[i]+epsilon;
//			temp[j] = c[j]-epsilon;
//			double plusminus = log_rosenbock(temp,dim, parameters);
//			temp[i] = c[i]-epsilon;
//			temp[j] = c[j]+epsilon;
//			double minusplus = log_rosenbock(temp,dim, parameters);
//			temp[i] = c[i];
//			temp[j] = c[j];
//			fisher[i][j]= -(plusplus -plusminus - minusplus + minusminus)/(4.*epsilon*epsilon);
//			//fisher[i][j]= 0;
//		}
//		temp[i]= c[i]+epsilon;
//		double plus = log_rosenbock(temp,dim,parameters);	
//		temp[i]= c[i]-epsilon;
//		double minus = log_rosenbock(temp,dim,parameters);	
//		temp[i]=c[i];
//		double eval = log_rosenbock(temp,dim,parameters);	
//		fisher[i][i]= -(plus -2*eval +minus)/(epsilon*epsilon);
//		//fisher[i][i]= 1;
//	}
//	for (int i = 0 ; i<dim ; i++){
//		for (int j = i+1 ; j<dim ; j++){
//			fisher[i][j] = fisher[j][i];	
//		}
//	}
//	
//	return;
//} 
double log_rosenbock (double *c,mcmc_data_interface *interface,void *parameters)
{
	int dim = interface->max_dim;
	rosenbock_param *param = (rosenbock_param *)parameters;
	double LL = 0;
	LL-= param->a * pow_int(c[0] - param->mu,2)	;
	for(int j = 0 ; j<param->n2; j++){
		for(int i = 1 ; i<param->n1; i++){
			if(i == 1){
				LL -= param->b[(j)*(param->n1-1) + i]*pow_int(c[(j)*(param->n1-1) + i] - c[0]*c[0],2);
			}
			else{
				LL -= param->b[(j)*(param->n1-1) + i]*pow_int(c[(j)*(param->n1-1) + i] - c[(j)*(param->n1-1) + i-1]*c[(j)*(param->n1-1) + i-1],2);

			}
		}
	}
	return LL;
	//return 2;
}
double log_rosenbock_prior (double *c,mcmc_data_interface *interface,void *parameters)
{
	int dim =interface->max_dim;
	double a = -std::numeric_limits<double>::infinity();
	rosenbock_param *param = (rosenbock_param *)parameters;
	for(int i = 0 ; i<param->n; i++){
		if( c[i]<param->bounds_low[i] || c[i] > param->bounds_high[i]){return a;}
	}
	return 0;
}
int mcmc_standard_test(int argc, char *argv[])
{
	int dimension = 2;
	//double initial_pos[2]={90,90.};
	double initial_pos[2]={0,0.};
	double *seeding_var = NULL;
	int N_steps = 5000;
	int chain_N= 30;
	int max_chain_N= 5;
	//double *initial_pos_ptr = initial_pos;
	int swp_freq = 5;
	//double chain_temps[chain_N] ={1,2,3,10,12};
	double chain_temps[chain_N];
	chain_temps[0] = 1.;
	double c = 1.8;
	for(int i =1; i < chain_N/2;  i ++)
		chain_temps[i] =  chain_temps[i-1] * c;
	chain_temps[chain_N/2] = 1;
	for(int i =chain_N/2+1; i < chain_N;  i ++)
		chain_temps[i] =  chain_temps[i-1] * c;
	std::string autocorrfile = "";
	std::string chainfile = "data/mcmc_output.hdf5";
	std::string statfilename = "data/mcmc_statistics.txt";
	std::string checkpointfile = "data/mcmc_checkpoint.csv";
	//std::string LLfile = "data/mcmc_LL.csv";
	std::string LLfile = "";
	
	int numThreads = 8;
	bool pool = true;
	bool show_progress = true;
	
	//double ***output;
	//output = allocate_3D_array( chain_N, N_steps, dimension );
	//PTMCMC_MH(output, dimension, N_steps, chain_N, initial_pos,seeding_var,chain_temps, swp_freq, log_test_prior, log_test,fisher_test,(void **)NULL,numThreads, pool,show_progress, statfilename,chainfile,autocorrfile, "",checkpointfile );	
	//deallocate_3D_array(output, chain_N, N_steps, dimension);
	double **output;
	output = allocate_2D_array(  N_steps, dimension );
	int t0 = 1000;
	int nu = 100;
	std::string chain_distribution_scheme="double";
	int max_chunk_size = 10000;
	mcmc_sampler_output sampler_output(chain_N,dimension);
	PTMCMC_MH_dynamic_PT_alloc_uncorrelated(&sampler_output,output, dimension, N_steps, chain_N, max_chain_N,initial_pos,seeding_var,(double**)NULL,chain_temps, swp_freq, t0,nu, max_chunk_size,chain_distribution_scheme, log_test_prior, log_test,fisher_test,(void **)NULL,numThreads, pool,show_progress, statfilename,chainfile, LLfile,checkpointfile );	
	//sampler_output.create_data_dump(false,false,chainfile);
	deallocate_2D_array(output, N_steps, dimension);
	//sampler_output.write_flat_thin_output( "./data/test_flat_standard",true,true);
	std::cout<<"ENDED"<<std::endl;


		
	return 0;

}
double log_test (double *c,mcmc_data_interface *interface,void *parameters)
{
	//return 2;
	double x = c[0];
	double y = c[1];
	double prefactor = 16./(M_PI*3.);
	double pow1 = -x*x - pow_int((9+4*x*x +8*y),2);
	double pow2 = -8*x*x -8*pow_int(y-2,2);
	return log(prefactor*(std::exp(pow1) + .5*std::exp(pow2)));
	//return 2.;
}
double log_test_prior (double *c,mcmc_data_interface *interface,void *parameters)
{
	double a = -std::numeric_limits<double>::infinity();
	if( fabs(c[0]) > 100){ return a;}
	if( fabs(c[1]) > 100){ return a;}
	//return -pow_int(c[0]-4,2)/4. -pow_int(c[1]+3,2)/6.;
	return 0.;
}

double fisher11(double *c)
{
	double x = c[0];
	double y = c[1];

	return 1.6976527263135504*(-8.*exp(-8*pow(x,2) - 8*pow(-2 + y,2)) + 128.*exp(-8*pow(x,2) - 8*pow(-2 + y,2))*pow(x,2) + 
     exp(-pow(x,2) - pow(9 + 4*pow(x,2) + 8*y,2))*(-2 - 128*pow(x,2) - 16*(9 + 4*pow(x,2) + 8*y)) + 
     exp(-pow(x,2) - pow(9 + 4*pow(x,2) + 8*y,2))*pow(-2*x - 16*x*(9 + 4*pow(x,2) + 8*y),2));
}
double fisher12(double *c)
{
	double x = c[0];
	double y = c[1];
	return 1.6976527263135504*(-128*exp(-pow(x,2) - pow(9 + 4*pow(x,2) + 8*y,2))*x + 128.*exp(-8*pow(x,2) - 8*pow(-2 + y,2))*x*(-2 + y) - 
     16*exp(-pow(x,2) - pow(9 + 4*pow(x,2) + 8*y,2))*(9 + 4*pow(x,2) + 8*y)*(-2*x - 16*x*(9 + 4*pow(x,2) + 8*y)));

}
double fisher22(double *c)
{
	double x = c[0];
	double y = c[1];
	return 1.6976527263135504*(-8.*exp(-8*pow(x,2) - 8*pow(-2 + y,2)) - 128*exp(-pow(x,2) - pow(9 + 4*pow(x,2) + 8*y,2)) + 
     128.*exp(-8*pow(x,2) - 8*pow(-2 + y,2))*pow(-2 + y,2) + 256*exp(-pow(x,2) - pow(9 + 4*pow(x,2) + 8*y,2))*pow(9 + 4*pow(x,2) + 8*y,2));

}
void fisher_test(double *c, double **fisher, mcmc_data_interface *interface, void *parameters)
{
	
	fisher[0][0] = fisher11(c);
	fisher[1][1] = fisher22(c);
	fisher[1][0] = fisher12(c);
	fisher[0][1] = fisher12(c);
	return;
} 
double standard_log_prior_D_sa(double *pos, mcmc_data_interface *interface,void *parameters)
{
	double a = -std::numeric_limits<double>::infinity();
	double chirp = std::exp(pos[0]);
	//double chirp = pos[7];
	double eta = pos[1];
	//double q = pos[8];
	if ((pos[1])<.01 || (pos[1])>.25){return a;}//eta
	if (pos[0]<.01  ){return a;}//chirpmass
	double m1 = calculate_mass1(chirp,eta);
	double m2 = calculate_mass2(chirp,eta);
	if (m1<.1 || m1>50){return a;}//PSI
	if (m2<.1 || m2>50){return a;}//PSI
	if ((pos[2])<-.9 || (pos[2])>.9){return a;}//eta
	if ((pos[3])<-.9 || (pos[3])>.9){return a;}//eta
	if ((pos[4])<-10 || (pos[4])>10){return a;}//eta
	//double chi1 = pos[9]+pos[10];
	//double chi2 = pos[9]-pos[10];
	//if ((chi1)<-.95 || (chi1)>.95){return a;}//chi1 
	//if ((chi2)<-.95 || (chi2)>.95){return a;}//chi2
	//else {return log(chirpmass_eta_jac(chirp,eta))+3*pos[6] ;}
	else {return log(chirpmass_eta_jac(chirp,eta)) ;}
	//else {return 3*pos[6] ;}
	//else {return log(chirpmass_eta_jac(chirp,eta))+3*pos[6] -log(cos(asin(pos[1]))) ;}

}
	
double standard_log_prior_D(double *pos, mcmc_data_interface *interface,void *parameters)
{
	double a = -std::numeric_limits<double>::infinity();
	double chirp = std::exp(pos[7]);
	//double chirp = pos[7];
	double eta = pos[8];
	//double q = pos[8];
	if ((pos[8])<.01 || (pos[8])>.25){return a;}//eta
	if (pos[7]<.01  ){return a;}//chirpmass
	double m1 = calculate_mass1(chirp,eta);
	double m2 = calculate_mass2(chirp,eta);
	if (m1<.1 || m1>50){return a;}//PSI
	if (m2<.1 || m2>50){return a;}//PSI
	if ((pos[0])<0 || (pos[0])>2*M_PI){ return a;}//RA

	if ((pos[1])<-1 || (pos[1])>1){return a;}//sinDEC
	//if ((pos[1])<-M_PI/2 || (pos[1])>M_PI/2){return a;}//sinDEC

	if ((pos[2])<0 || (pos[2])>2*M_PI){return a;}//PSI
	if ((pos[3])<-1 || (pos[3])>1){return a;}//cos \iota
	if ((pos[4])<0 || (pos[4])>2*M_PI){return a;}//phiRef
	if ((pos[5])<T_mcmc_gw_tool*3./4.-.1 || (pos[5])>T_mcmc_gw_tool*3./4.+.1){return a;}//tc
	//if ((pos[5])<0 || (pos[5])>T_mcmc_gw_tool){return a;}//tc
	if (std::exp(pos[6])<10 || std::exp(pos[6])>5000){return a;}//DL
	//if (std::exp(pos[7])<.01 || std::exp(pos[7])>5 ){return a;}//chirpmass
	//if (pos[7]<.01 || pos[7]>14 ){return a;}//chirpmass
	//if (pos[7]<.1 || pos[7]>10 ){return a;}//chirpmass
	//if ((pos[8])<.01 || (pos[8])>.249999){return a;}//eta
	//if ((pos[8])<.01 || (pos[8])>1){return a;}//eta
	//if ((pos[8])<.01 || (pos[8])>10){return a;}//eta
	//if (pos[8]>pos[7]){return a;}//m1 greater than m2
	if ((pos[9])<-.9 || (pos[9])>.9){return a;}//eta
	if ((pos[10])<-.9 || (pos[10])>.9){return a;}//eta
	for(int i = 11 ; i<interface->max_dim; i++){

		if ((pos[i])<-10 || (pos[i])>10){return a;}//eta
	}
	//double chi1 = pos[9]+pos[10];
	//double chi2 = pos[9]-pos[10];
	//if ((chi1)<-.95 || (chi1)>.95){return a;}//chi1 
	//if ((chi2)<-.95 || (chi2)>.95){return a;}//chi2
	//else {return log(chirpmass_eta_jac(chirp,eta))+3*pos[6] ;}
	return log(chirpmass_eta_jac(chirp,eta))+3*pos[6] ;
	//else {return 3*pos[6] ;}
	//else {return log(chirpmass_eta_jac(chirp,eta))+3*pos[6] -log(cos(asin(pos[1]))) ;}

}
double standard_log_prior_dCS(double *pos, mcmc_data_interface *interface,void *parameters)
{
	double chirp = std::exp(pos[7]);
	double eta = pos[8];
	double a = -std::numeric_limits<double>::infinity();
	if ((pos[0])<0 || (pos[0])>2*M_PI){ return a;}//RA

	if ((pos[1])<-1 || (pos[1])>1){return a;}//sinDEC
	//if ((pos[1])<-M_PI/2 || (pos[1])>M_PI/2){return a;}//sinDEC

	if ((pos[2])<0 || (pos[2])>M_PI){return a;}//PSI
	if ((pos[3])<-1 || (pos[3])>1){return a;}//cos \iota
	if ((pos[4])<0 || (pos[4])>2*M_PI){return a;}//phiRef
	//if ((pos[5])<T_mcmc_gw_tool*3./4.-.1 || (pos[5])>T_mcmc_gw_tool*3./4.+.1){return a;}//tc
	if ((pos[5])<0 || (pos[5])>T_mcmc_gw_tool){return a;}//tc
	if (std::exp(pos[6])<10 || std::exp(pos[6])>500){return a;}//DL
	if (std::exp(pos[7])<.1 || std::exp(pos[7])>20 ){return a;}//chirpmass
	if ((pos[8])<.01 || (pos[8])>.249999){return a;}//eta
	if ((pos[9])<0 || (pos[9])>.95){return a;}//chi1 
	if ((pos[10])<0 || (pos[10])>.95){return a;}//chi2
	if ((pos[11])<-1 || (pos[11])>1){return a;}//chi2
	if ((pos[12])<-1 || (pos[12])>1){return a;}//chi2
	if ((pos[13])<0 || (pos[13])>2*M_PI){return a;}//chi2
	if(pos[14]<0 || pos[14] > 100) { return a;}
	else {return log(chirpmass_eta_jac(chirp,eta))+3*pos[6] ;}
	//else {return log(chirpmass_eta_jac(chirp,eta))+3*pos[6] -log(cos(asin(pos[1]))) ;}
}
double standard_log_prior_P(double *pos, mcmc_data_interface *interface,void *parameters)
{
	double chirp = std::exp(pos[7]);
	double eta = pos[8];
	double a = -std::numeric_limits<double>::infinity();
	if ((pos[8])<.01 || (pos[8])>1){return a;}//eta
	if (chirp<.01  ){return a;}//chirpmass
	double m1 = calculate_mass1(chirp,eta);
	double m2 = calculate_mass2(chirp,eta);
	if (m1<.1 || m1>5){return a;}//PSI
	if (m2<.1 || m2>5){return a;}//PSI
	if ((pos[0])<0 || (pos[0])>2*M_PI){ return a;}//RA

	if ((pos[1])<-1 || (pos[1])>1){return a;}//sinDEC
	//if ((pos[1])<-M_PI/2 || (pos[1])>M_PI/2){return a;}//sinDEC

	if ((pos[2])<0 || (pos[2])>M_PI){return a;}//PSI
	if ((pos[3])<-1 || (pos[3])>1){return a;}//cos \iota
	if ((pos[4])<0 || (pos[4])>2*M_PI){return a;}//phiRef
	//if ((pos[5])<T_mcmc_gw_tool*3./4.-.1 || (pos[5])>T_mcmc_gw_tool*3./4.+.1){return a;}//tc
	if ((pos[5])<0 || (pos[5])>T_mcmc_gw_tool){return a;}//tc
	if (std::exp(pos[6])<10 || std::exp(pos[6])>800){return a;}//DL
	//if (std::exp(pos[7])<.01 || std::exp(pos[7])>15 ){return a;}//chirpmass
	//if (pos[7]<2 || pos[7]>80 ){return a;}//chirpmass
	//if ((pos[8])<.01 || (pos[8])>.25){return a;}//eta
	if ((pos[9])<0 || (pos[9])>.95){return a;}//chi1 
	if ((pos[10])<0 || (pos[10])>.95){return a;}//chi2
	if ((pos[11])<-1 || (pos[11])>1){return a;}//chi2
	if ((pos[12])<-1 || (pos[12])>1){return a;}//chi2
	if ((pos[13])<0 || (pos[13])>2*M_PI){return a;}//chi2
	if ((pos[14])<0 || (pos[14])>2*M_PI){return a;}//chi2
	if ((pos[15])<-10 || (pos[15])>10){return a;}//chi2
	//else {return log(chirp*chirpmass_q_jac(chirp,q))+3*pos[6] ;}
	else {return log(chirpmass_eta_jac(chirp,eta))+3*pos[6] ;}
	//else {return log(chirpmass_q_jac(chirp,q))+3*pos[6] ;}
	//else {return log(chirpmass_eta_jac(chirp,eta))+3*pos[6] -log(cos(asin(pos[1]))) ;}
}
//Uniform in m1 and m2, transformed to lnM and eta
double chirpmass_eta_jac(double chirpmass, double eta){
	return chirpmass*chirpmass/(sqrt(1. - 4.*eta)*pow(eta,1.2));
}
//Uniform in m1 and m2, transformed to lnM and q
double chirpmass_q_jac(double chirpmass, double q){
	//return (pow(q/pow_int(q+1,2),1./5.) * q);
	//return chirpmass*chirpmass/(pow(q/pow_int(q+1,2),1./5.) * q);
	return chirpmass/(pow(q/pow_int(q+1,2),1./5.) * q);
}

double standard_log_prior_D_RJ(double *pos,int *status,int model_status, mcmc_data_interface *interface,void *parameters)
{
	double a = -std::numeric_limits<double>::infinity();
	double chirp = std::exp(pos[7]);
	double eta = pos[8];

	if (eta<.0 || eta>.25){return a;}//eta
	double m1 = calculate_mass1(chirp,eta );
	double m2 = calculate_mass2(chirp,eta );
	if(m1<0.01 || m1>50){return a;}
	if(m2<0.01 || m2>50){return a;}

	if ((pos[0])<0 || (pos[0])>2*M_PI){ return a;}//RA

	if ((pos[1])<-1 || (pos[1])>1){return a;}//sinDEC

	if ((pos[2])<0 || (pos[2])>M_PI){return a;}//PSI
	if ((pos[3])<-1 || (pos[3])>1){return a;}//cos \iota
	if ((pos[4])<0 || (pos[4])>2*M_PI){return a;}//phiRef
	if ((pos[5])<T_mcmc_gw_tool*3./4. -.1 || (pos[5])>3.*T_mcmc_gw_tool/4. + .1){return a;}//tc
	if (std::exp(pos[6])<10 || std::exp(pos[6])>5000){return a;}//DL
	if ((pos[9])<-.9 || (pos[9])>.95){return a;}//chi1 
	if ((pos[10])<-.9 || (pos[10])>.95){return a;}//chi2
	for(int i = interface->min_dim; i<interface->max_dim; i++){
		if(status[i] ==1 ){
			if ((pos[i])<-10 || (pos[i])>10){return a;}//gIMR
		}
	}
	return log(chirpmass_eta_jac(chirp,eta))+3*pos[6] ;
}
double standard_log_prior_P_RJ(double *pos,int *status,int model_status, mcmc_data_interface *interface,void *parameters)
{
	double a = -std::numeric_limits<double>::infinity();
	double chirp = std::exp(pos[7]);
	double eta = pos[8];

	if (eta<.0 || eta>.25){return a;}//eta
	double m1 = calculate_mass1(chirp,eta );
	double m2 = calculate_mass2(chirp,eta );
	if(m1<0.01 || m1>50){return a;}
	if(m2<0.01 || m2>50){return a;}

	if ((pos[0])<0 || (pos[0])>2*M_PI){ return a;}//RA

	if ((pos[1])<-1 || (pos[1])>1){return a;}//sinDEC

	if ((pos[2])<0 || (pos[2])>M_PI){return a;}//PSI
	if ((pos[3])<-1 || (pos[3])>1){return a;}//cos \iota
	if ((pos[4])<0 || (pos[4])>2*M_PI){return a;}//phiRef
	if ((pos[5])<T_mcmc_gw_tool*3./4. -.1 || (pos[5])>3.*T_mcmc_gw_tool/4. + .1){return a;}//tc
	if (std::exp(pos[6])<10 || std::exp(pos[6])>5000){return a;}//DL
	if ((pos[9])<0 || (pos[9])>.95){return a;}//chi1 
	if ((pos[10])<0 || (pos[10])>.95){return a;}//chi2
	if ((pos[11])<-1 || (pos[11])>1){return a;}//chi2
	if ((pos[12])<-1 || (pos[12])>1){return a;}//chi2
	if ((pos[13])<0 || (pos[13])>2*M_PI){return a;}//chi2
	if ((pos[14])<0 || (pos[14])>2*M_PI){return a;}//chi2
	for(int i = interface->min_dim; i<interface->max_dim; i++){
		if(status[i] ==1 ){
			if ((pos[i])<-10 || (pos[i])>10){return a;}//gIMR
		}
	}
	return log(chirpmass_eta_jac(chirp,eta))+3*pos[6] ;
}
void RT_ERROR_MSG()
{
	std::cout<<"ERROR -- incorrect arguments"<<std::endl;
	std::cout<<"Please supply function option:"<<std::endl;
	std::cout<<"0 --- Standard"<<std::endl;
	std::cout<<"1 --- GW injection"<<std::endl;
	std::cout<<"2 --- GW experimentation"<<std::endl;
	std::cout<<"3 --- N-dimensional Rosenbock"<<std::endl;
	std::cout<<"4 --- MCMC output class test"<<std::endl;
	std::cout<<"5 --- MCMC RJ Sine Wave test"<<std::endl;
	std::cout<<"6 --- MCMC RJ GW injection test"<<std::endl;
	std::cout<<"7 --- Multiple continue mcmc testing"<<std::endl;
	std::cout<<"8 --- Test likelihood"<<std::endl;
	std::cout<<"9 --- Test evidence calculation"<<std::endl;
	std::cout<<"10 --- Validate evidence calculation"<<std::endl;
	std::cout<<"11 --- Ensemble size testing"<<std::endl;
	std::cout<<"12 --- stupid test"<<std::endl;
	std::cout<<"13 --- 1D test"<<std::endl;
	std::cout<<"14 --- PTRJMCMC GW injection test"<<std::endl;
}
