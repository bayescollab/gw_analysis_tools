#include "mcmc_gw_extended.h"


#include "mcmc_gw.h"
#include "waveform_generator.h"
#include "util.h"
#include "io_util.h"
#include "detector_util.h"
#include "ppE_utilities.h"
#include "waveform_util.h"
#include "ortho_basis.h"
#include "fisher.h"

#include <bayesship/bayesshipSampler.h>
#include <bayesship/dataUtilities.h>
#include <bayesship/proposalFunctions.h>
#include <bayesship/utilities.h>

void MCMC_fisher_wrapper_v2(bayesship::positionInfo *pos,   double **output, void *userParameters);

void MCMC_fisher_wrapper_v3(bayesship::positionInfo *pos,   double **output, std::vector<int> ids, void *userParameters);

int invertFisherBlock(double **fisherIn, double **fisherOut, int dimIn, std::vector<int> ids);

void MCMC_fisher_wrapper_RJ_v2(bayesship::positionInfo *pos,   double **output, std::vector<int> block, void *userParameters);

void MCMC_fisher_wrapper_RJ_v3(bayesship::positionInfo *pos,   double **output, std::vector<int> block, void *userParameters);

void pack_local_mod_structure_v2(bayesship::bayesshipSampler *sampler,
	double *param,
	int *status,
	std::string waveform_extended, 
	void *parameters, 
	MCMC_modification_struct *full_struct, 
	MCMC_modification_struct *local_struct )	
{
	if(waveform_extended.find("gIMR") != std::string::npos){
		int dimct = 0 ;
		int dphi_boundary = full_struct->gIMR_Nmod_phi + sampler->minDim;
		int dsigma_boundary = full_struct->gIMR_Nmod_sigma + dphi_boundary;
		int dbeta_boundary = full_struct->gIMR_Nmod_beta + dsigma_boundary;
		int dalpha_boundary = full_struct->gIMR_Nmod_alpha + dbeta_boundary;
		for(int i = 0 ; i<sampler->maxDim; i++){
			if(status[i] == 1){
				dimct++;
			}
			if( i >= sampler->minDim){
				if(status[i] == 1 && i<dphi_boundary){
					local_struct->gIMR_Nmod_phi ++;
				}
				else if(status[i] == 1 && i<dsigma_boundary){
					local_struct->gIMR_Nmod_sigma ++;
				}
				else if(status[i] == 1 && i<dbeta_boundary){
					local_struct->gIMR_Nmod_beta ++;
				}
				else if(status[i] == 1 && i<dalpha_boundary){
					local_struct->gIMR_Nmod_alpha ++;
				}
			}
		}
		if(dimct != sampler->minDim){
			if(local_struct->gIMR_Nmod_phi != 0){
				local_struct->gIMR_phii = new int[local_struct->gIMR_Nmod_phi];
			}
			if(local_struct->gIMR_Nmod_sigma != 0){
				local_struct->gIMR_sigmai = new int[local_struct->gIMR_Nmod_sigma];
			}
			if(local_struct->gIMR_Nmod_beta != 0){
				local_struct->gIMR_betai = new int[local_struct->gIMR_Nmod_beta];
			}
			if(local_struct->gIMR_Nmod_alpha != 0){
				local_struct->gIMR_alphai = new int[local_struct->gIMR_Nmod_alpha];
			}
			
			int ct_phi = 0 ;
			int ct_sigma = 0 ;
			int ct_beta = 0 ;
			int ct_alpha = 0 ;
			for(int i =sampler->minDim ; i<sampler->maxDim; i++){
				if(status[i] == 1){
					if(i < dphi_boundary){
						local_struct->gIMR_phii[ct_phi] = full_struct->gIMR_phii[i-dphi_boundary+full_struct->gIMR_Nmod_phi];
						ct_phi++;
					}
					else if(i < dsigma_boundary){
						local_struct->gIMR_sigmai[ct_sigma] = full_struct->gIMR_sigmai[i-dsigma_boundary+full_struct->gIMR_Nmod_sigma];
						ct_sigma++;
					}
					else if(i < dbeta_boundary){
						local_struct->gIMR_betai[ct_beta] = full_struct->gIMR_betai[i-dbeta_boundary+full_struct->gIMR_Nmod_beta];
						ct_beta++;
					}
					else if(i < dalpha_boundary){
						local_struct->gIMR_alphai[ct_alpha] = full_struct->gIMR_alphai[i-dalpha_boundary+full_struct->gIMR_Nmod_alpha];
						ct_alpha++;
					}
				}
			}
		}
	}
	else if(waveform_extended.find("ppE") != std::string::npos){
		int dimct = 0 ;
		local_struct->ppE_Nmod=0;
		for(int i = 0 ; i<sampler->maxDim; i++){
			if(status[i] == 1){
				dimct++;
			}

			if( i >= sampler->minDim){
				if(status[i] == 1 ){
					local_struct->ppE_Nmod++;
				}
			}
		}
		if(dimct != sampler->minDim){
			local_struct->bppe = new double[local_struct->ppE_Nmod];
			
			int ct_ppe = 0 ;
			for(int i =sampler->minDim ; i<sampler->maxDim; i++){
				if(status[i] == 1){
					local_struct->bppe[ct_ppe] = full_struct->bppe[i-sampler->minDim];
					ct_ppe++;
				}
			}
		}
	}

	return;
}


class MCMC_likelihood_wrapper_v2_RJ: public bayesship::probabilityFn
{
public:
	bayesship::bayesshipSampler *sampler;
	mcmcVariablesRJ *mcmcVarRJ;
	virtual double eval(bayesship::positionInfo *pos, int chainID)
	{
		//return 2;
		//mcmcVariables *mcmcVar = (mcmcVariables *)userParameters;
		//MCMC_user_param *user_param = (MCMC_user_param *)userParameters;
	
		int maxDim = sampler->maxDim;
		int minDim = sampler->minDim;
		double ll = 0;
		double *temp_params = new double[maxDim];
		
		int dim = pos->countActiveDimensions();		
		
		int ct=0;
		for(int i = 0 ; i<maxDim;i++){
			if(pos->status[i]){
				temp_params[ct] = pos->parameters[i];
				//std::cout<<temp_params[ct]<<",";
				ct+=1;
			}
		}
		//std::cout<<std::endl;
		bool wfExtended=false;
		std::string gen_meth = mcmcVarRJ->mcmc_generation_method;
		if(dim > minDim){wfExtended=true;gen_meth=mcmcVarRJ->mcmc_generation_method_extended;}
		MCMC_modification_struct mod_struct_local;

		if(wfExtended){
			pack_local_mod_structure_v2(sampler,pos->parameters, pos->status, mcmcVarRJ->mcmc_generation_method_extended,mcmcVarRJ->user_parameters , mcmcVarRJ->mcmc_mod_struct, &mod_struct_local);
		}
		//#########################################################################
		gen_params_base<double> gen_params;
		std::string local_gen = MCMC_prep_params_v2(temp_params, 
			temp_params,&gen_params, dim, gen_meth,&mod_struct_local, mcmcVarRJ->mcmc_intrinsic,mcmcVarRJ->mcmc_gmst );
		//#########################################################################
		//#########################################################################
	
		//repack_non_parameters(temp_params, &gen_params, 
			//"MCMC_"+mcmc_generation_method, dimension, NULL);
		repack_parameters(temp_params, &gen_params, 
			"MCMC_"+gen_meth, dim, NULL);
		//if(gen_params.Nmod !=0){
		//	std::cout<<gen_params.Nmod<<" "<<local_gen<<std::endl;
		//	for(int i = 0 ; i<gen_params.Nmod; i++){
		//		std::cout<<" "<<gen_params.bppe[i]<<" "<<gen_params.betappe[i];
		//	}
		//	std::cout<<std::endl;
		//	
		//}
		//#########################################################################
		//#########################################################################
		//return 1;
		std::complex<double> **local_data = mcmcVarRJ->mcmc_data;
		double **local_freqs = mcmcVarRJ->mcmc_frequencies;
		double **local_noise = mcmcVarRJ->mcmc_noise;
		double **local_weights = (mcmcVarRJ->user_parameters)->weights;
		int *local_lengths = mcmcVarRJ->mcmc_data_length;
		fftw_outline *local_plans = mcmcVarRJ->mcmc_fftw_plans;
		std::string local_integration_method="SIMPSONS";
		//if(interface->burn_phase && user_param->burn_data){
		if(false){
			local_data = mcmcVarRJ->user_parameters->burn_data;
			local_freqs = mcmcVarRJ->user_parameters->burn_freqs;
			local_noise = mcmcVarRJ->user_parameters->burn_noise;
			local_lengths = mcmcVarRJ->user_parameters->burn_lengths;
			local_plans = mcmcVarRJ->user_parameters->burn_plans;
		}
		if(mcmcVarRJ->user_parameters->GAUSS_QUAD){
			local_integration_method="GAUSSLEG";
		}
		if(mcmcVarRJ->mcmc_intrinsic){
			if(local_gen.find("IMRPhenomD") != std::string::npos){
				if(!mcmcVarRJ->mcmc_save_waveform){
					for(int i=0; i < mcmcVarRJ->mcmc_num_detectors; i++){
						gen_params.theta=0;	
						gen_params.phi=0;	
						gen_params.psi=0;	
						gen_params.phiRef = 1;
						gen_params.f_ref = 10;
						gen_params.incl_angle=0;	
						gen_params.tc =1;
						std::complex<double> *response =
							(std::complex<double> *) malloc(sizeof(std::complex<double>) * local_lengths[i]);
						fourier_detector_response_horizon(local_freqs[i], local_lengths[i], response, mcmcVarRJ->mcmc_detectors[i], local_gen, &gen_params);
						ll += maximized_Log_Likelihood_aligned_spin_internal(local_data[i], 
								local_noise[i],
								local_freqs[i],
								response,
								(size_t) local_lengths[i],
								&local_plans[i]
								);
						//ll += maximized_Log_Likelihood(mcmc_data[i], 
						//		mcmc_noise[i],
						//		mcmc_frequencies[i],
						//		(size_t) mcmc_data_length[i],
						//		&gen_params,
						//		mcmc_detectors[i],
						//		local_gen,
						//		&mcmc_fftw_plans[i]
						//		);
						free(response);
					}
				}
				else{
					gen_params.theta=0;	
					gen_params.phi=0;	
					gen_params.psi=0;	
					gen_params.phiRef = 1;
					gen_params.f_ref = 10;
					gen_params.incl_angle=0;	
					gen_params.tc =1;
					std::complex<double> *response =
						(std::complex<double> *) malloc(sizeof(std::complex<double>) * local_lengths[0]);
					fourier_detector_response_horizon(local_freqs[0], local_lengths[0], response, mcmcVarRJ->mcmc_detectors[0], local_gen, &gen_params);
					//std::complex<double> *hc =
					//	(std::complex<double> *) malloc(sizeof(std::complex<double>) * mcmc_data_length[0]);
					//std::complex<double> *hp =
					//	(std::complex<double> *) malloc(sizeof(std::complex<double>) * mcmc_data_length[0]);
					//fourier_waveform(mcmc_frequencies[0], mcmc_data_length[0], hp,hc, local_gen, &gen_params);
					for(int i=0; i < mcmcVarRJ->mcmc_num_detectors; i++){
						ll += maximized_Log_Likelihood_aligned_spin_internal(local_data[i], 
								local_noise[i],
								local_freqs[i],
								response,
								(size_t) local_lengths[i],
								&local_plans[i]
								);
						//ll += maximized_Log_Likelihood(mcmc_data[i], 
						//		mcmc_noise[i],
						//		mcmc_frequencies[i],
						//		(size_t) mcmc_data_length[i],
						//		&gen_params,
						//		mcmc_detectors[i],
						//		local_gen,
						//		&mcmc_fftw_plans[i]
						//		);
						//ll += maximized_Log_Likelihood_unaligned_spin_internal(mcmc_data[i], 
						//		mcmc_noise[i],
						//		mcmc_frequencies[i],
						//		hp,
						//		hc,
						//		(size_t) mcmc_data_length[i],
						//		&mcmc_fftw_plans[i]
						//		);
					}
					free(response);
					//free(hp);
					//free(hc);
	
				}
	
			}
			else if(local_gen.find("IMRPhenomP")!=std::string::npos){
				//if(!mcmc_save_waveform){
				if(false){
				}
				else{
					gen_params.theta=0;	
					gen_params.phi=0;	
					gen_params.psi=0;	
					gen_params.phiRef = 1;
					gen_params.f_ref = 20;
					gen_params.incl_angle=0;	
					gen_params.tc =1;
					waveform_polarizations<double> wp;
					assign_polarizations(local_gen,&wp);
					wp.allocate_memory(local_lengths[0]);
					fourier_waveform(local_freqs[0],local_lengths[0], &wp, local_gen, &gen_params);
					for(int i=0; i < mcmcVarRJ->mcmc_num_detectors; i++){
						ll += maximized_Log_Likelihood_unaligned_spin_internal(local_data[i], 
								local_noise[i],
								local_freqs[i],
								wp.hplus,
								wp.hcross,
								(size_t) local_lengths[i],
								&local_plans[i]
								);
					}
					wp.deallocate_memory();
				}
	
			}
		}
		else{
			double RA = gen_params.RA;
			double DEC = gen_params.DEC;
			double PSI = gen_params.psi;
			//if(mcmc_generation_method.find("IMRPhenomD") != std::string:npos){
			
			ll =  MCMC_likelihood_extrinsic(mcmcVarRJ->mcmc_save_waveform, 
				&gen_params,local_gen, local_lengths, 
				local_freqs, local_data, local_noise, local_weights, local_integration_method, mcmcVarRJ->user_parameters->log10F,mcmcVarRJ->mcmc_detectors, 
				 mcmcVarRJ->mcmc_num_detectors);
			//ll=2;
	
			//ll = Log_Likelihood(mcmc_data[0], 
			//		mcmc_noise[0],
			//		mcmc_frequencies[0],
			//		mcmc_data_length[0],
			//		&gen_params,
			//		mcmc_detectors[0],
			//		local_gen,
			//		&mcmc_fftw_plans[0]
			//		);
	
			//}
			//else if(mcmc_generation_method.find("IMRPhenomP")!=std::string::npos){
	
			//}
		}
		//Cleanup
		delete [] temp_params;
		if(check_mod(local_gen)){
			//if( local_gen.find("ppE") != std::string::npos ||
			//	local_gen.find("dCS") != std::string::npos ||
			//	local_gen.find("EdGB") != std::string::npos){
			//	delete [] gen_params.betappe;
			//}
			if( local_gen.find("ppE") != std::string::npos ||
				check_theory_support(local_gen)){
				delete [] gen_params.betappe;
				delete [] mod_struct_local.bppe;
			}
			else if( local_gen.find("gIMR") != std::string::npos){
				if(mod_struct_local.gIMR_Nmod_phi !=0){
					delete [] gen_params.delta_phi;
					delete [] mod_struct_local.gIMR_phii;
				}
				if(mod_struct_local.gIMR_Nmod_sigma !=0){
					delete [] gen_params.delta_sigma;
					delete [] mod_struct_local.gIMR_sigmai;
				}
				if(mod_struct_local.gIMR_Nmod_beta !=0){
					delete [] gen_params.delta_beta;
					delete [] mod_struct_local.gIMR_betai;
				}
				if(mod_struct_local.gIMR_Nmod_alpha !=0){
					delete [] gen_params.delta_alpha;
					delete [] mod_struct_local.gIMR_alphai;
				}
	
			}
		}
		//debugger_print(__FILE__,__LINE__,ll);
		return ll;
	
	}
};

class ppEFisherRJVariables
{
public:
	double *bppe=nullptr;
	int fisherDim;
	double **fisher=nullptr;
	int detectN;
	double **psds=nullptr;
	double **freqs=nullptr;
	int length;
	ppEFisherRJVariables(double *bppe, double **psds, double **freqs,int detectN,int fisherDim,int length)
	{
		this->detectN = detectN;
		this->fisherDim = fisherDim;
		this->length = length;
		this->bppe = new double[fisherDim];
		for(int i = 0 ; i<fisherDim; i++){
			this->bppe[i] = bppe[i];
		}
		this->psds = new double*[detectN];
		for(int i = 0 ; i<detectN; i++){
			this->psds[i] = new double[length];
			for(int j = 0 ; j<length; j++){
				this->psds[i][j]  = psds[i][j];
			}
		}
		this->freqs = new double*[detectN];
		for(int i = 0 ; i<detectN; i++){
			this->freqs[i] = new double[length];
			for(int j = 0 ; j<length; j++){
				this->freqs[i][j]  = freqs[i][j];
			}
		}
	
	
		this->fisher = new double*[fisherDim];
		for(int i = 0 ; i<fisherDim ; i++){
			this->fisher[i] = new double[fisherDim];
			for(int j = 0 ; j<fisherDim ; j++){
				this->fisher[i][j] = 0;
			}
		}
	
	
		//Calculate the actual fisher
		double *integrand = new double[length];
		for(int l = 0 ; l< detectN;l++){
			for(int i = 0 ; i < fisherDim; i++){
				for(int j = 0 ; j <= i; j++){
					for(int k = 0 ; k<length; k++){
						integrand[k] = ( pow(freqs[l][k], bppe[i]/3.+bppe[j]/3. -7/3.)/psds[l][k] );
						integrand[k]*=(2./15.) * pow(M_PI, -4./3. + bppe[i]/3. + bppe[j]/3.);
					}
					fisher[i][j] += simpsons_sum(freqs[l][1]-freqs[l][0], length, integrand);
				}
			}	
		}
		for(int i = 0 ; i < fisherDim; i++){
			for(int j = 0 ; j <= i; j++){
				fisher[j][i] = fisher[i][j];
			}
		}	
		//std::cout<<"Fisher: "<<std::endl;
		//for(int i = 0 ; i < fisherDim; i++){
		//	for(int j = 0 ; j <fisherDim; j++){
		//		std::cout<<fisher[i][j] <<" , ";
		//	}
		//	std::cout<<std::endl;
		//}	
		
		
	
		delete [] integrand;
	};
		~ppEFisherRJVariables(){
			if(bppe){
				delete [] bppe;	
				bppe = nullptr;
			}
			if(fisher){
				for(int i = 0 ; i<fisherDim ; i++){	
					delete [] fisher[i];	
				}
				delete [] fisher;	
				fisher = nullptr;
			}
			if(psds){
				for(int i = 0 ; i<detectN ; i++){	
					delete [] psds[i];	
				}
				delete [] psds;	
				psds = nullptr;
			}
			if(freqs){
				for(int i = 0 ; i<detectN ; i++){	
					delete [] freqs[i];	
				}
				delete [] freqs;	
				freqs = nullptr;
			}
		};
	};




void MCMC_fisher_wrapper_RJ_ppE(bayesship::positionInfo *pos,   double **output, std::vector<int> block, void *userParameters)
{
	double chirp = std::exp(pos->parameters[7])*MSOL_SEC;
	double DL = std::exp(pos->parameters[6])*MPC_SEC;
	ppEFisherRJVariables *p = (ppEFisherRJVariables *) userParameters;
	for(int i = 0 ; i<p->fisherDim ; i++){
		for(int j = 0 ; j<p->fisherDim; j++){
			output[i][j] = p->fisher[i][j] * pow(chirp, 4.- 7./3. + p->bppe[i]/3. + p->bppe[j]/3.)/DL/DL;
		}
	}
	//std::cout<<"Fisher: "<<std::endl;
	//for(int i = 0 ; i < p->fisherDim; i++){
	//	for(int j = 0 ; j <p->fisherDim; j++){
	//		std::cout<<output[i][j] <<" , ";
	//	}
	//	std::cout<<std::endl;
	//}	
	return;
}

class MCMC_likelihood_wrapper_v2: public bayesship::probabilityFn
{
public:
	bayesship::bayesshipSampler *sampler;
	mcmcVariables *mcmcVar;
	virtual double eval(bayesship::positionInfo *pos, int chainID)
	{
		//return 2;
		//mcmcVariables *mcmcVar = (mcmcVariables *)userParameters;
		//MCMC_user_param *user_param = (MCMC_user_param *)userParameters;
	
		int dimension = sampler->maxDim;
		double ll = 0;
		double *temp_params = new double[dimension];
		//#########################################################################
		gen_params_base<double> gen_params;
		std::string local_gen = MCMC_prep_params_v2(pos->parameters, 
			temp_params,&gen_params, dimension, mcmcVar->mcmc_generation_method,mcmcVar->mcmc_mod_struct, mcmcVar->mcmc_intrinsic,mcmcVar->mcmc_gmst );
		//#########################################################################
		//#########################################################################
	
		//repack_non_parameters(temp_params, &gen_params, 
			//"MCMC_"+mcmc_generation_method, dimension, NULL);
		repack_parameters(temp_params, &gen_params, 
			"MCMC_"+mcmcVar->mcmc_generation_method, dimension, NULL);
		//#########################################################################
		//#########################################################################
		//return 1;
		std::complex<double> **local_data = mcmcVar->mcmc_data;
		double **local_freqs = mcmcVar->mcmc_frequencies;
		double **local_noise = mcmcVar->mcmc_noise;
		double **local_weights = (mcmcVar->user_parameters)->weights;
		int *local_lengths = mcmcVar->mcmc_data_length;
		fftw_outline *local_plans = mcmcVar->mcmc_fftw_plans;
		std::string local_integration_method="SIMPSONS";
		//if(interface->burn_phase && user_param->burn_data){
		if(false){
			local_data = mcmcVar->user_parameters->burn_data;
			local_freqs = mcmcVar->user_parameters->burn_freqs;
			local_noise = mcmcVar->user_parameters->burn_noise;
			local_lengths = mcmcVar->user_parameters->burn_lengths;
			local_plans = mcmcVar->user_parameters->burn_plans;
		}
		if(mcmcVar->user_parameters->GAUSS_QUAD){
			local_integration_method="GAUSSLEG";
		}
		if(mcmcVar->mcmc_intrinsic){
			if(mcmcVar->mcmc_generation_method.find("IMRPhenomD") != std::string::npos){
				if(!mcmcVar->mcmc_save_waveform){
					for(int i=0; i < mcmcVar->mcmc_num_detectors; i++){
						gen_params.theta=0;	
						gen_params.phi=0;	
						gen_params.psi=0;	
						gen_params.phiRef = 1;
						gen_params.f_ref = 10;
						gen_params.incl_angle=0;	
						gen_params.tc =1;
						std::complex<double> *response =
							(std::complex<double> *) malloc(sizeof(std::complex<double>) * local_lengths[i]);
						fourier_detector_response_horizon(local_freqs[i], local_lengths[i], response, mcmcVar->mcmc_detectors[i], local_gen, &gen_params);
						ll += maximized_Log_Likelihood_aligned_spin_internal(local_data[i], 
								local_noise[i],
								local_freqs[i],
								response,
								(size_t) local_lengths[i],
								&local_plans[i]
								);
						//ll += maximized_Log_Likelihood(mcmc_data[i], 
						//		mcmc_noise[i],
						//		mcmc_frequencies[i],
						//		(size_t) mcmc_data_length[i],
						//		&gen_params,
						//		mcmc_detectors[i],
						//		local_gen,
						//		&mcmc_fftw_plans[i]
						//		);
						free(response);
					}
				}
				else{
					gen_params.theta=0;	
					gen_params.phi=0;	
					gen_params.psi=0;	
					gen_params.phiRef = 1;
					gen_params.f_ref = 10;
					gen_params.incl_angle=0;	
					gen_params.tc =1;
					std::complex<double> *response =
						(std::complex<double> *) malloc(sizeof(std::complex<double>) * local_lengths[0]);
					fourier_detector_response_horizon(local_freqs[0], local_lengths[0], response, mcmcVar->mcmc_detectors[0], local_gen, &gen_params);
					//std::complex<double> *hc =
					//	(std::complex<double> *) malloc(sizeof(std::complex<double>) * mcmc_data_length[0]);
					//std::complex<double> *hp =
					//	(std::complex<double> *) malloc(sizeof(std::complex<double>) * mcmc_data_length[0]);
					//fourier_waveform(mcmc_frequencies[0], mcmc_data_length[0], hp,hc, local_gen, &gen_params);
					for(int i=0; i < mcmcVar->mcmc_num_detectors; i++){
						ll += maximized_Log_Likelihood_aligned_spin_internal(local_data[i], 
								local_noise[i],
								local_freqs[i],
								response,
								(size_t) local_lengths[i],
								&local_plans[i]
								);
						//ll += maximized_Log_Likelihood(mcmc_data[i], 
						//		mcmc_noise[i],
						//		mcmc_frequencies[i],
						//		(size_t) mcmc_data_length[i],
						//		&gen_params,
						//		mcmc_detectors[i],
						//		local_gen,
						//		&mcmc_fftw_plans[i]
						//		);
						//ll += maximized_Log_Likelihood_unaligned_spin_internal(mcmc_data[i], 
						//		mcmc_noise[i],
						//		mcmc_frequencies[i],
						//		hp,
						//		hc,
						//		(size_t) mcmc_data_length[i],
						//		&mcmc_fftw_plans[i]
						//		);
					}
					free(response);
					//free(hp);
					//free(hc);
	
				}
	
			}
			else if(mcmcVar->mcmc_generation_method.find("IMRPhenomP")!=std::string::npos){
				//if(!mcmc_save_waveform){
				if(false){
				}
				else{
					gen_params.theta=0;	
					gen_params.phi=0;	
					gen_params.psi=0;	
					gen_params.phiRef = 1;
					gen_params.f_ref = 20;
					gen_params.incl_angle=0;	
					gen_params.tc =1;
					waveform_polarizations<double> wp;
					assign_polarizations(mcmcVar->mcmc_generation_method,&wp);
					wp.allocate_memory(local_lengths[0]);
					fourier_waveform(local_freqs[0],local_lengths[0], &wp, local_gen, &gen_params);
					for(int i=0; i < mcmcVar->mcmc_num_detectors; i++){
						ll += maximized_Log_Likelihood_unaligned_spin_internal(local_data[i], 
								local_noise[i],
								local_freqs[i],
								wp.hplus,
								wp.hcross,
								(size_t) local_lengths[i],
								&local_plans[i]
								);
					}
					wp.deallocate_memory();
				}
	
			}
		}
		else if (mcmcVar->mcmc_adaptive)
		{
			ll = mcmcVar->adaptivell->log_likelihood(
				mcmcVar->mcmc_detectors, mcmcVar->mcmc_num_detectors,
				&gen_params, local_gen, mcmcVar->mcmc_save_waveform
				);
		}
		else{
			double RA = gen_params.RA;
			double DEC = gen_params.DEC;
			double PSI = gen_params.psi;
			//if(mcmc_generation_method.find("IMRPhenomD") != std::string:npos){
			
			ll =  MCMC_likelihood_extrinsic(mcmcVar->mcmc_save_waveform, 
				&gen_params,local_gen, local_lengths, 
				local_freqs, local_data, local_noise, local_weights, local_integration_method, mcmcVar->user_parameters->log10F,mcmcVar->mcmc_detectors, 
				 mcmcVar->mcmc_num_detectors, mcmcVar->QuadMethod);
			//ll=2;
	
			//ll = Log_Likelihood(mcmc_data[0], 
			//		mcmc_noise[0],
			//		mcmc_frequencies[0],
			//		mcmc_data_length[0],
			//		&gen_params,
			//		mcmc_detectors[0],
			//		local_gen,
			//		&mcmc_fftw_plans[0]
			//		);
	
			//}
			//else if(mcmc_generation_method.find("IMRPhenomP")!=std::string::npos){
	
			//}
		}
		//Cleanup
		delete [] temp_params;
		if(check_mod(local_gen)){
			//if( local_gen.find("ppE") != std::string::npos ||
			//	local_gen.find("dCS") != std::string::npos ||
			//	local_gen.find("EdGB") != std::string::npos){
			//	delete [] gen_params.betappe;
			//}
			if( local_gen.find("ppE") != std::string::npos ||
				check_theory_support(local_gen)){
				delete [] gen_params.betappe;
			}
			else if( local_gen.find("gIMR") != std::string::npos){
				if(mcmcVar->mcmc_mod_struct ->gIMR_Nmod_phi !=0){
					delete [] gen_params.delta_phi;
				}
				if(mcmcVar->mcmc_mod_struct ->gIMR_Nmod_sigma !=0){
					delete [] gen_params.delta_sigma;
				}
				if(mcmcVar->mcmc_mod_struct ->gIMR_Nmod_beta !=0){
					delete [] gen_params.delta_beta;
				}
				if(mcmcVar->mcmc_mod_struct ->gIMR_Nmod_alpha !=0){
					delete [] gen_params.delta_alpha;
				}
	
			}
		}
		if ( isnan(ll)){
			return bayesship::limitInf;
		}
		//std::cout<<ll<<std::endl;
		return ll;
	
	}
};

/*! \brief Takes in an MCMC checkpoint file and continues the chain
 *
 * Obviously, the user must be sure to correctly match the dimension, number of chains, the generation_method, 
 * the prior function, the data, psds, freqs, and the detectors (number and name), and the gps_time to the 
 * previous run, otherwise the behavior of the sampler is undefined.
 *
 * numThreads and pool do not necessarily have to be the same
 */
bayesship::bayesshipSampler *  RJPTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW_v2(
	int minDim,
	int maxDim,
	int independentSamples,
	int ensembleSize,
	int ensembleN,
	bayesship::positionInfo *initialPosition,
	bayesship::positionInfo **initialEnsemble,
	double swapProb,
	int burnIterations,
	int burnPriorIterations,
	int priorIterations,
	bool writePriorData,
	int batchSize,
	double **priorRanges,
	bayesship::probabilityFn *log_prior,
	int numThreads,
	bool pool,
	int num_detectors,
	std::complex<double> **data,
	double **noise_psd,
	double **frequencies,
	int *data_length,
	double gps_time,
	std::string *detectors,
	MCMC_modification_struct *mod_struct,
	std::string generation_method,
	std::string generation_method_extended,
	std::string outputDir,
	std::string outputFileMoniker,
	bool ignoreExistingCheckpoint,
	bool restrictSwapTemperatures,
	bool coldChainStorageOnly
	)
{
	int chainN = ensembleSize*ensembleN;
	//Create fftw plan for each detector (length of data stream may be different)
	fftw_outline *plans= (fftw_outline *)malloc(sizeof(fftw_outline)*num_detectors);
	for (int i =0;i<num_detectors;i++)
	{	
		allocate_FFTW_mem_forward(&plans[i] , data_length[i]);
	}

	std::mutex fisher_mutex;

//##########################################################
//##########################################################
	mcmcVariablesRJ mcmcVarRJ ;
	mcmcVarRJ.mcmc_noise = noise_psd;
	//mcmcVarRJ.mcmc_init_pos = initial_pos;
	mcmcVarRJ.mcmc_frequencies = frequencies;
	mcmcVarRJ.mcmc_data = data;
	mcmcVarRJ.mcmc_data_length = data_length;
	mcmcVarRJ.mcmc_detectors = detectors;
	mcmcVarRJ.mcmc_generation_method = generation_method;
	mcmcVarRJ.mcmc_generation_method_extended = generation_method_extended;
	mcmcVarRJ.mcmc_fftw_plans = plans;
	mcmcVarRJ.mcmc_num_detectors = num_detectors;
	mcmcVarRJ.mcmc_gps_time = gps_time;
	mcmcVarRJ.mcmc_gmst = gps_to_GMST_radian(gps_time);
	mcmcVarRJ.mcmc_mod_struct = mod_struct;
	mcmcVarRJ.mcmc_save_waveform = true;
	mcmcVarRJ.maxDim = maxDim;
	mcmcVarRJ.minDim = minDim;



//##########################################################
//##########################################################




	//mcmc_noise = noise_psd;	
	//mcmc_init_pos = initial_pos;
	//mcmc_frequencies = frequencies;
	//mcmc_data = data;
	//mcmc_data_length = data_length;
	//mcmc_detectors = detectors;
	//mcmc_generation_method = generation_method;
	//mcmc_fftw_plans = plans;
	//mcmc_num_detectors = num_detectors;
	//mcmc_gps_time = gps_time;
	//mcmc_gmst = gps_to_GMST_radian(mcmc_gps_time);
	////mcmc_Nmod = mod_struct->ppE_Nmod;
	////mcmc_bppe = mod_struct->bppe;
	//mcmc_mod_struct = mod_struct;
	//mcmc_log_beta = false;
	//mcmc_intrinsic = false;

	//To save time, intrinsic waveforms can be saved between detectors, if the 
	//frequencies are all the same
	mcmc_save_waveform = true;
	for(int i =1 ;i<mcmcVarRJ.mcmc_num_detectors; i++){
		if( mcmcVarRJ.mcmc_data_length[i] != mcmcVarRJ.mcmc_data_length[0] ||
			mcmcVarRJ.mcmc_frequencies[i][0]!= mcmcVarRJ.mcmc_frequencies[0][0] ||
			mcmcVarRJ.mcmc_frequencies[i][mcmcVarRJ.mcmc_data_length[i]-1] 
				!= mcmcVarRJ.mcmc_frequencies[0][mcmcVarRJ.mcmc_data_length[0]-1]){
			mcmcVarRJ.mcmc_save_waveform= false;
		}
			
	}


	//std::cout<<"Base:"<<std::endl;
	RJPTMCMC_method_specific_prep_v2(generation_method, maxDim,  &(mcmcVarRJ.mcmc_intrinsic), mcmcVarRJ.mcmc_mod_struct);
	//std::cout<<"Extending up to:"<<std::endl;
	//RJPTMCMC_method_specific_prep_v2(generation_method_extended, maxDim,  &(mcmcVarRJ.mcmc_intrinsic), mcmcVarRJ.mcmc_mod_struct);
	
	//######################################################
	int T = (int)(1./(mcmcVarRJ.mcmc_frequencies[0][1]-mcmcVarRJ.mcmc_frequencies[0][0]));
	debugger_print(__FILE__,__LINE__,T);
	int burn_factor = T/4; //Take all sources to 4 seconds
	debugger_print(__FILE__,__LINE__,burn_factor);
	std::complex<double> **burn_data = new std::complex<double>*[mcmcVarRJ.mcmc_num_detectors];
	double **burn_freqs = new double*[mcmcVarRJ.mcmc_num_detectors];
	double **burn_noise = new double*[mcmcVarRJ.mcmc_num_detectors];
	int *burn_lengths = new int[mcmcVarRJ.mcmc_num_detectors];
	fftw_outline *burn_plans= new fftw_outline[mcmcVarRJ.mcmc_num_detectors];
	for(int j = 0; j<mcmcVarRJ.mcmc_num_detectors; j++){
		burn_lengths[j] = mcmcVarRJ.mcmc_data_length[j]/burn_factor;
		burn_data[j]= new std::complex<double>[burn_lengths[j]];
		burn_freqs[j]= new double[burn_lengths[j]];
		burn_noise[j]= new double[burn_lengths[j]];
		allocate_FFTW_mem_forward(&burn_plans[j], burn_lengths[j]);
		int ct = 0;
		for( int k = 0 ; k<mcmcVarRJ.mcmc_data_length[j]; k++){
			if(k%burn_factor==0 && ct<burn_lengths[j]){
				burn_data[j][ct] = mcmcVarRJ.mcmc_data[j][k];
				burn_freqs[j][ct] = mcmcVarRJ.mcmc_frequencies[j][k];
				burn_noise[j][ct] = mcmcVarRJ.mcmc_noise[j][k];
				ct++;
			}
		}
	}
	

	MCMC_user_param **user_parameters=NULL;
	user_parameters = new MCMC_user_param*[chainN];
	for(int i = 0 ;i<chainN; i++){
		user_parameters[i] = new MCMC_user_param;
		
		user_parameters[i]->burn_data = burn_data;
		user_parameters[i]->burn_freqs = burn_freqs;
		user_parameters[i]->burn_noise = burn_noise;
		user_parameters[i]->burn_lengths = burn_lengths;
		user_parameters[i]->burn_plans=burn_plans;

		user_parameters[i]->mFish= &fisher_mutex;
		user_parameters[i]->GAUSS_QUAD= mod_struct->GAUSS_QUAD;
		user_parameters[i]->log10F = mod_struct->log10F;

		if(mod_struct->weights){
			user_parameters[i]->weights = mod_struct->weights;			
		}
		else{
			user_parameters[i]->weights = new double*[num_detectors];			
			for(int j = 0 ; j<num_detectors; j++){
				user_parameters[i]->weights[j]=NULL;
			}
		}

		user_parameters[i]->fisher_GAUSS_QUAD = mod_struct->fisher_GAUSS_QUAD;
		user_parameters[i]->fisher_log10F = mod_struct->fisher_log10F;
		user_parameters[i]->fisher_freq= mod_struct->fisher_freq;
		if(mod_struct->fisher_weights){
			user_parameters[i]->fisher_weights= mod_struct->fisher_weights;
		}
		else{
			user_parameters[i]->fisher_weights = new double*[num_detectors];
			for(int j = 0 ; j<num_detectors; j++){
				user_parameters[i]->fisher_weights[j]=NULL;
			}
		}	
		user_parameters[i]->fisher_PSD= mod_struct->fisher_PSD;
		user_parameters[i]->fisher_length= mod_struct->fisher_length;


		user_parameters[i]->mod_struct = mod_struct;

		
		//user_parameters[i]->burn_freqs = mcmc_frequencies;
		//user_parameters[i]->burn_data = mcmc_data;
		//user_parameters[i]->burn_noise = mcmc_noise;
		//user_parameters[i]->burn_lengths = mcmc_data_length;
	}
	//mcmcVarRJ.user_parameters = user_parameters;
	//######################################################
	//###########################################################
	//double trial[11] = {1,.9,.2,.9,.2,3,7,4,.24,0,0};
	//mcmc_data_interface i;
	//i.min_dim = 11; 
	//i.max_dim = 11; 
	//i.chain_id  = 0; 
	//i.chain_number  = 11; 
	//double ll = MCMC_likelihood_wrapper(trial, &i ,user_parameters[0]);
	//debugger_print(__FILE__,__LINE__,ll);
	//###########################################################
	
	MCMC_likelihood_wrapper_v2_RJ *ll = new MCMC_likelihood_wrapper_v2_RJ();
	ll->mcmcVarRJ = &mcmcVarRJ;
	bayesship::bayesshipSampler *sampler = new bayesship::bayesshipSampler(ll,log_prior );
	sampler->RJ = true;
	ll->sampler = sampler;

	sampler->iterations = independentSamples;
	sampler->batchSize = batchSize;
	sampler->outputDir = outputDir;
	sampler->outputFileMoniker = outputFileMoniker;
	sampler->burnIterations = burnIterations;
	sampler->burnPriorIterations = burnPriorIterations;
	sampler->priorIterations = priorIterations;
	sampler->writePriorData = writePriorData;
	sampler->threads = numThreads;
	sampler->threadPool = pool;
	sampler->maxDim = maxDim;
	sampler->minDim = minDim;
	sampler->ensembleSize = ensembleSize;
	sampler->ensembleN = ensembleN;
	sampler->priorRanges = priorRanges;
	sampler->initialPosition = initialPosition;
	sampler->initialPositionEnsemble = initialEnsemble;
	sampler->ignoreExistingCheckpoint = ignoreExistingCheckpoint;
	sampler->restrictSwapTemperatures = restrictSwapTemperatures;

	//Testing
	sampler->coldOnlyStorage = coldChainStorageOnly;
	//sampler->coldOnlyStorage = true;

	mcmcVariablesRJ **mcmcVarVec  = new mcmcVariablesRJ*[chainN];
	for(int i = 0 ; i<chainN; i++){
		mcmcVarVec[i] = &mcmcVarRJ;
		mcmcVarVec[i]->user_parameters = user_parameters[i];
	}
	sampler->userParameters = (void **) mcmcVarVec;

	//##########################################################
	
	int proposalFnN = 6;
	bayesship::proposal **propArray = new bayesship::proposal*[proposalFnN];
	propArray[0] = new bayesship::gaussianProposal(sampler->ensembleN*sampler->ensembleSize, sampler->maxDim, sampler);
	//propArray[1] = new bayesship::differentialEvolutionProposal(sampler);
	if(mcmcVarRJ.mcmc_intrinsic){
		propArray[1] = new bayesship::differentialEvolutionProposal(sampler);
	}
	else{
		std::vector<std::vector<int>> blocksDiff = std::vector<std::vector<int>>(3);	
		for(int i = 0 ; i<7; i++){
			blocksDiff[0].push_back(i);
		}
		for(int i = 7 ; i<sampler->minDim; i++){
			blocksDiff[1].push_back(i);
		}
		for(int i = 0 ; i<sampler->minDim; i++){
			blocksDiff[2].push_back(i);
		}
		std::vector<double> blocksProbDiff = {0.3,0.3,.4};
		propArray[1] = new bayesship::blockDifferentialEvolutionProposal(sampler, blocksDiff,blocksProbDiff);
	}

	propArray[2] = new bayesship::KDEProposal(sampler->ensembleN*sampler->ensembleSize, sampler->maxDim, sampler, false );

	propArray[3] = new bayesship::randomLayerRJProposal(sampler, .5);


	//################################################
	//std::vector<std::vector<int>> blocks = {{0,1,2,3,4,5,6,7,8,9,10}};
	//std::vector<double> blockProb = {1};
	//propArray[4] = new bayesship::blockFisherProposal(sampler->ensembleN*sampler->ensembleSize, sampler->minDim, &MCMC_fisher_wrapper_RJ_v3,   sampler->userParameters,  100,sampler,blocks, blockProb );

	if(mcmcVarRJ.mcmc_intrinsic){
		std::vector<std::vector<int>> blocks = std::vector<std::vector<int>>(1);
		for(int i = 0 ; i<sampler->minDim; i++){
			blocks[0].push_back(i);
		}
		std::vector<double> blockProb = {1};
		propArray[4] = new bayesship::blockFisherProposal(sampler->ensembleN*sampler->ensembleSize, sampler->minDim, &MCMC_fisher_wrapper_RJ_v3,   sampler->userParameters,  100,sampler,blocks, blockProb );
	}
	else{
		std::vector<std::vector<int>> blocks = std::vector<std::vector<int>>(3);
		for(int i = 0 ; i<7; i++){
			blocks[0].push_back(i);
		}
		for(int i = 7 ; i<sampler->minDim; i++){
			blocks[1].push_back(i);
		}
		for(int i = 0 ; i<sampler->minDim; i++){
			blocks[2].push_back(i);
		}
		std::vector<double> blockProb = {.3,.3,.4};
		//std::vector<std::vector<int>> blocks = {
		//				{7,8,9,10}};
		//std::vector<double> blockProb = {1};
		propArray[4] = new bayesship::blockFisherProposal(sampler->ensembleN*sampler->ensembleSize, sampler->minDim, &MCMC_fisher_wrapper_RJ_v3,   sampler->userParameters,  100,sampler,blocks, blockProb );
	}

	//################################################
	//
	std::vector<std::vector<int>> blocks2 = std::vector<std::vector<int>>(1);
	blocks2[0] = std::vector<int>(maxDim-minDim);
	for(int i = 0 ; i<maxDim-minDim; i++){
		blocks2[0][i] = minDim + i;
	}
	std::vector<double> blockProb2 = {1};
	ppEFisherRJVariables *ppEFisherObj = new ppEFisherRJVariables(mod_struct->bppe, mcmcVarRJ.mcmc_noise,mcmcVarRJ.mcmc_frequencies, mcmcVarRJ.mcmc_num_detectors, maxDim-minDim, mcmcVarRJ.mcmc_data_length[0]);
	ppEFisherRJVariables **ppEFisherObjs = new ppEFisherRJVariables*[chainN];
	for(int i = 0 ; i<chainN ;i++){
		ppEFisherObjs[i] = ppEFisherObj;
	}
	
	propArray[5] = new bayesship::blockFisherProposal(sampler->ensembleN*sampler->ensembleSize, sampler->minDim, &MCMC_fisher_wrapper_RJ_ppE,   (void**)ppEFisherObjs,  100,sampler,blocks2, blockProb2 );


	//Rough estimate of the temperatures
	double betaTemp[sampler->ensembleSize];
	betaTemp[0] = 1;
	betaTemp[sampler->ensembleSize-1] = 0;
	double deltaBeta = pow((1e2),1./sampler->ensembleSize);
	for(int i = 1 ; i<sampler->ensembleSize-1; i++){
		betaTemp[i] = betaTemp[i-1]/deltaBeta;
	}

	double **propProb = new double*[chainN];	
	for(int i = 0 ; i<chainN; i++){
		int ensemble = i / sampler->ensembleN;
		propProb[i] = new double[proposalFnN];	
		propProb[i][2] = 0.0;

		propProb[i][1] = .6 - 0.45*( betaTemp[ensemble]); //.15 to .7
		propProb[i][3] = 0.15 - .1*( betaTemp[ensemble]); //.05 to .2
		propProb[i][4] = 0.05 + .4*( betaTemp[ensemble]); //.45 to .05
		propProb[i][5] = 0.05 + .2*( betaTemp[ensemble]); //.25 to .05
		//std::cout<<"TESTING FISHER"<<std::endl;
		//propProb[i][1] = 0.0;
		//propProb[i][3] = 0.0;
		//propProb[i][4] = 0.0;
		//propProb[i][5] = 1.0;

		//propProb[i][0] = 0.05;
		//propProb[i][1] = 0.25;
		//propProb[i][3] = 0.7;
		//propProb[i][3] = 0; //.1 to .2
		//propProb[i][4] = 0;


		//propProb[i][1] = 0; //.25 to .7
		//propProb[i][3] = 0.3 + .45*( betaTemp[ensemble]); //.7 to .25

		//propProb[i][0] = 1.  - propProb[i][1]- propProb[i][2];
		
		propProb[i][0] = 0.05;
		
		double sum = 0 ;
		for(int j = 0 ; j<proposalFnN; j++){
			sum+=propProb[i][j];	
		}
		for(int j = 0 ; j<proposalFnN; j++){
			propProb[i][j]/=sum;
		}
		for(int j = 0 ; j<proposalFnN; j++){
			std::cout<<propProb[i][j]<<", ";
		}
		std::cout<<"\n";

	}
	
	//propProb[0] = 0.05;
	//propProb[1] = 0.3;
	//propProb[2] = 0.05;
	//propProb[3] = 0.6;
	//propProb[0] = 0.05;
	//propProb[1] = 0.25;
	//propProb[2] = 0.0;
	//propProb[3] = 0.7;

	//propProb[0] = 0.1;
	//propProb[1] = 0.0;
	//propProb[2] = 0.0;
	//propProb[3] = 0.9;



	bayesship::proposalData *propData = new bayesship::proposalData(chainN,proposalFnN, propArray, (double *)nullptr, propProb);
	
	
	
	

	//#######################################################################################
	sampler->proposalFns = propData;
	//#######################################################################################



	//#######################################################################################
	//#######################################################################################
	sampler->sample();
	//#######################################################################################
	//#######################################################################################
	
	
	//#######################################################################################
	//Delete data
	//#######################################################################################

	for(int i =0 ; i<chainN ; i++){
		delete [] propProb[i];
	}
	delete [] propProb;


	
	

	if(!mcmcVarRJ.mcmc_mod_struct->fisher_weights){
		for(int i = 0 ;i<chainN;i++){
			delete[] user_parameters[i]->fisher_weights;
		}
	}
	if(!mcmcVarRJ.mcmc_mod_struct->weights){
		for(int i = 0 ;i<chainN;i++){
			delete[] user_parameters[i]->weights;
		}
	}
	
	//Deallocate fftw plans
	for (int i =0;i<num_detectors;i++)
		deallocate_FFTW_mem(&plans[i]);
	//#################################################
	for(int i = 0 ; i<proposalFnN; i++){
		delete propData->proposals[i];
	}
	delete [] ppEFisherObjs;
	delete ppEFisherObj;
	delete propData;
	delete [] propArray;
	for(int i = 0 ; i<num_detectors; i++){
		delete [] burn_data[i];
		delete [] burn_freqs[i];
		delete [] burn_noise[i];
		deallocate_FFTW_mem(&burn_plans[i]);
	}
	delete [] burn_data;
	delete [] burn_lengths;
	delete [] burn_noise;
	delete [] burn_freqs;
	delete [] burn_plans;
	for(int i = 0 ; i<chainN; i++){
		delete user_parameters[i];
	}
	delete [] user_parameters;
	delete [] mcmcVarVec;
	delete ll;
	//#################################################
	free(plans);
	return sampler;
}



/*! \brief Takes in an MCMC checkpoint file and continues the chain
 *
 * Obviously, the user must be sure to correctly match the dimension, number of chains, the generation_method, 
 * the prior function, the data, psds, freqs, and the detectors (number and name), and the gps_time to the 
 * previous run, otherwise the behavior of the sampler is undefined.
 *
 * numThreads and pool do not necessarily have to be the same
 */
bayesship::bayesshipSampler *  PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW_v2(
	int dimension,
	int independentSamples,
	int ensembleSize,
	int ensembleN,
	bayesship::positionInfo *initialPosition,
	bayesship::positionInfo **initialEnsemble,
	double swapProb,
	int burnIterations,
	int burnPriorIterations,
	int priorIterations,
	bool writePriorData,
	int batchSize,
	double **priorRanges,
	bayesship::probabilityFn *log_prior,
	int numThreads,
	bool pool,
	int num_detectors,
	std::complex<double> **data,
	double **noise_psd,
	double **frequencies,
	int *data_length,
	double gps_time,
	std::string *detectors,
	MCMC_modification_struct *mod_struct,
	std::string generation_method,
	std::string outputDir,
	std::string outputFileMoniker,
	bool ignoreExistingCheckpoint,
	bool restrictSwapTemperatures,
	bool coldChainStorageOnly
	)
{
	int chainN = ensembleSize*ensembleN;
	//Create fftw plan for each detector (length of data stream may be different)
	fftw_outline *plans= (fftw_outline *)malloc(sizeof(fftw_outline)*num_detectors);
	for (int i =0;i<num_detectors;i++)
	{	
		allocate_FFTW_mem_forward(&plans[i] , data_length[i]);
	}

	std::mutex fisher_mutex;

//##########################################################
//##########################################################
	mcmcVariables mcmcVar ;
	mcmcVar.mcmc_noise = noise_psd;
	//mcmcVar.mcmc_init_pos = initial_pos;
	mcmcVar.mcmc_frequencies = frequencies;
	mcmcVar.mcmc_data = data;
	mcmcVar.mcmc_data_length = data_length;
	mcmcVar.mcmc_detectors = detectors;
	mcmcVar.mcmc_generation_method = generation_method;
	mcmcVar.mcmc_fftw_plans = plans;
	mcmcVar.mcmc_num_detectors = num_detectors;
	mcmcVar.mcmc_gps_time = gps_time;
	mcmcVar.mcmc_gmst = gps_to_GMST_radian(gps_time);
	mcmcVar.mcmc_mod_struct = mod_struct;
	mcmcVar.mcmc_save_waveform = true;
	mcmcVar.maxDim = dimension;
	mcmcVar.QuadMethod = mod_struct->QuadMethod;

	if (mod_struct->adaptivell != nullptr)
	{
		std::cout << "Sampling with adaptive likelihood\n";
		mcmcVar.adaptivell = mod_struct->adaptivell;
		mcmcVar.mcmc_adaptive = true;
	}


//##########################################################
//##########################################################




	//mcmc_noise = noise_psd;	
	//mcmc_init_pos = initial_pos;
	//mcmc_frequencies = frequencies;
	//mcmc_data = data;
	//mcmc_data_length = data_length;
	//mcmc_detectors = detectors;
	//mcmc_generation_method = generation_method;
	//mcmc_fftw_plans = plans;
	//mcmc_num_detectors = num_detectors;
	//mcmc_gps_time = gps_time;
	//mcmc_gmst = gps_to_GMST_radian(mcmc_gps_time);
	////mcmc_Nmod = mod_struct->ppE_Nmod;
	////mcmc_bppe = mod_struct->bppe;
	//mcmc_mod_struct = mod_struct;
	//mcmc_log_beta = false;
	//mcmc_intrinsic = false;

	//To save time, intrinsic waveforms can be saved between detectors, if the 
	//frequencies are all the same
	mcmc_save_waveform = true;
	for(int i =1 ;i<mcmcVar.mcmc_num_detectors; i++){
		if( mcmcVar.mcmc_data_length[i] != mcmcVar.mcmc_data_length[0] ||
			mcmcVar.mcmc_frequencies[i][0]!= mcmcVar.mcmc_frequencies[0][0] ||
			mcmcVar.mcmc_frequencies[i][mcmcVar.mcmc_data_length[i]-1] 
				!= mcmcVar.mcmc_frequencies[0][mcmcVar.mcmc_data_length[0]-1]){
			mcmcVar.mcmc_save_waveform= false;
		}
			
	}


	//PTMCMC_method_specific_prep(generation_method, dimension, seeding_var_ptr, local_seeding);
	PTMCMC_method_specific_prep_v2(generation_method, dimension, &(mcmcVar.mcmc_intrinsic), mcmcVar.mcmc_mod_struct);
	
	//######################################################
	int T = (int)(1./(mcmcVar.mcmc_frequencies[0][1]-mcmcVar.mcmc_frequencies[0][0]));
	debugger_print(__FILE__,__LINE__,T);
	int burn_factor = T/4; //Take all sources to 4 seconds
	debugger_print(__FILE__,__LINE__,burn_factor);
	std::complex<double> **burn_data = new std::complex<double>*[mcmcVar.mcmc_num_detectors];
	double **burn_freqs = new double*[mcmcVar.mcmc_num_detectors];
	double **burn_noise = new double*[mcmcVar.mcmc_num_detectors];
	int *burn_lengths = new int[mcmcVar.mcmc_num_detectors];
	fftw_outline *burn_plans= new fftw_outline[mcmcVar.mcmc_num_detectors];
	for(int j = 0; j<mcmcVar.mcmc_num_detectors; j++){
		burn_lengths[j] = mcmcVar.mcmc_data_length[j]/burn_factor;
		burn_data[j]= new std::complex<double>[burn_lengths[j]];
		burn_freqs[j]= new double[burn_lengths[j]];
		burn_noise[j]= new double[burn_lengths[j]];
		allocate_FFTW_mem_forward(&burn_plans[j], burn_lengths[j]);
		int ct = 0;
		for( int k = 0 ; k<mcmcVar.mcmc_data_length[j]; k++){
			if(k%burn_factor==0 && ct<burn_lengths[j]){
				burn_data[j][ct] = mcmcVar.mcmc_data[j][k];
				burn_freqs[j][ct] = mcmcVar.mcmc_frequencies[j][k];
				burn_noise[j][ct] = mcmcVar.mcmc_noise[j][k];
				ct++;
			}
		}
	}
	

	MCMC_user_param **user_parameters=NULL;
	user_parameters = new MCMC_user_param*[chainN];
	for(int i = 0 ;i<chainN; i++){
		user_parameters[i] = new MCMC_user_param;
		
		user_parameters[i]->burn_data = burn_data;
		user_parameters[i]->burn_freqs = burn_freqs;
		user_parameters[i]->burn_noise = burn_noise;
		user_parameters[i]->burn_lengths = burn_lengths;
		user_parameters[i]->burn_plans=burn_plans;

		user_parameters[i]->mFish= &fisher_mutex;
		user_parameters[i]->GAUSS_QUAD= mod_struct->GAUSS_QUAD;
		user_parameters[i]->log10F = mod_struct->log10F;

		if(mod_struct->weights){
			user_parameters[i]->weights = mod_struct->weights;			
		}
		else{
			user_parameters[i]->weights = new double*[num_detectors];			
			for(int j = 0 ; j<num_detectors; j++){
				user_parameters[i]->weights[j]=NULL;
			}
		}

		user_parameters[i]->fisher_GAUSS_QUAD = mod_struct->fisher_GAUSS_QUAD;
		user_parameters[i]->fisher_log10F = mod_struct->fisher_log10F;
		user_parameters[i]->fisher_freq= mod_struct->fisher_freq;
		if(mod_struct->fisher_weights){
			user_parameters[i]->fisher_weights= mod_struct->fisher_weights;
		}
		else{
			user_parameters[i]->fisher_weights = new double*[num_detectors];
			for(int j = 0 ; j<num_detectors; j++){
				user_parameters[i]->fisher_weights[j]=NULL;
			}
		}	
		user_parameters[i]->fisher_PSD= mod_struct->fisher_PSD;
		user_parameters[i]->fisher_length= mod_struct->fisher_length;
		user_parameters[i]->QuadMethod = mod_struct->QuadMethod;

		user_parameters[i]->mod_struct = mod_struct;

		
		//user_parameters[i]->burn_freqs = mcmc_frequencies;
		//user_parameters[i]->burn_data = mcmc_data;
		//user_parameters[i]->burn_noise = mcmc_noise;
		//user_parameters[i]->burn_lengths = mcmc_data_length;
	}
	//mcmcVar.user_parameters = user_parameters;
	//######################################################
	//###########################################################
	//double trial[11] = {1,.9,.2,.9,.2,3,7,4,.24,0,0};
	//mcmc_data_interface i;
	//i.min_dim = 11; 
	//i.max_dim = 11; 
	//i.chain_id  = 0; 
	//i.chain_number  = 11; 
	//double ll = MCMC_likelihood_wrapper(trial, &i ,user_parameters[0]);
	//debugger_print(__FILE__,__LINE__,ll);
	//###########################################################
	
	MCMC_likelihood_wrapper_v2 *ll = new MCMC_likelihood_wrapper_v2();
	ll->mcmcVar = &mcmcVar;
	bayesship::bayesshipSampler *sampler = new bayesship::bayesshipSampler(ll,log_prior );
	ll->sampler = sampler;
	sampler->independentSamples = independentSamples;
	sampler->batchSize = batchSize;
	sampler->outputDir = outputDir;
	sampler->outputFileMoniker = outputFileMoniker;
	sampler->burnIterations = burnIterations;
	sampler->burnPriorIterations = burnPriorIterations;
	sampler->priorIterations = priorIterations;
	sampler->writePriorData = writePriorData;
	sampler->threads = numThreads;
	sampler->threadPool = pool;
	sampler->maxDim = dimension;
	sampler->ensembleSize = ensembleSize;
	sampler->ensembleN = ensembleN;
	sampler->priorRanges = priorRanges;
	sampler->initialPosition = initialPosition;
	sampler->initialPositionEnsemble = initialEnsemble;
	sampler->ignoreExistingCheckpoint = ignoreExistingCheckpoint;
	sampler->restrictSwapTemperatures = restrictSwapTemperatures;

	//Testing
	//sampler->coldOnlyStorage = false;
	sampler->coldOnlyStorage = coldChainStorageOnly;

	mcmcVariables **mcmcVarVec  = new mcmcVariables*[chainN];
	for(int i = 0 ; i<chainN; i++){
		mcmcVarVec[i] = &mcmcVar;
		mcmcVarVec[i]->user_parameters = user_parameters[i];
	}
	sampler->userParameters = (void **) mcmcVarVec;

	//##########################################################
	
	int proposalFnN = 5;
	bayesship::proposal **propArray = new bayesship::proposal*[proposalFnN];
	propArray[0] = new bayesship::gaussianProposal(sampler->ensembleN*sampler->ensembleSize, sampler->maxDim, sampler);
	//propArray[1] = new bayesship::differentialEvolutionProposal(sampler);
	if(mcmcVar.mcmc_intrinsic){
		propArray[1] = new bayesship::differentialEvolutionProposal(sampler);
	}
	else if(generation_method.find("EA_IMRPhenomD_NRT") != std::string::npos ){
		std::vector<std::vector<int>> blocksDiff = std::vector<std::vector<int>>(4);	
		for(int i = 0 ; i<7; i++){
			blocksDiff[0].push_back(i);
		}
		for(int i = 7 ; i<sampler->maxDim; i++){
			blocksDiff[1].push_back(i);
		}
		for(int i = 0 ; i<sampler->maxDim; i++){
			blocksDiff[2].push_back(i);
		}
		blocksDiff[3].push_back(7);
		blocksDiff[3].push_back(12);
		std::vector<double> blocksProbDiff = {0.25,0.25,.25,.25};
		propArray[1] = new bayesship::blockDifferentialEvolutionProposal(sampler, blocksDiff,blocksProbDiff);
		propArray[4] = new bayesship::GMMProposal( sampler->ensembleN*sampler->ensembleSize, sampler->maxDim, sampler, blocksDiff,blocksProbDiff, 10, 10, 10, 1e-10, false, sampler->maxDim*100);

	}
	else{
		std::vector<std::vector<int>> blocksDiff = std::vector<std::vector<int>>(3);	
		for(int i = 0 ; i<7; i++){
			blocksDiff[0].push_back(i);
		}
		for(int i = 7 ; i<sampler->maxDim; i++){
			blocksDiff[1].push_back(i);
		}
		for(int i = 0 ; i<sampler->maxDim; i++){
			blocksDiff[2].push_back(i);
		}
		std::vector<double> blocksProbDiff = {0.3,0.3,.4};
		propArray[1] = new bayesship::blockDifferentialEvolutionProposal(sampler, blocksDiff,blocksProbDiff);
		propArray[4] = new bayesship::GMMProposal( sampler->ensembleN*sampler->ensembleSize, sampler->maxDim, sampler, blocksDiff,blocksProbDiff, 10, 10, 10, 1e-10, false, sampler->maxDim*100);
	}
	propArray[2] = new bayesship::KDEProposal(sampler->ensembleN*sampler->ensembleSize, sampler->maxDim, sampler, false );
	//propArray[3] = new bayesship::fisherProposal(sampler->ensembleN*sampler->ensembleSize, sampler->maxDim, &MCMC_fisher_wrapper_v2,   sampler->userParameters,  100,sampler);
	if(mcmcVar.mcmc_intrinsic){
		propArray[3] = new bayesship::fisherProposal(sampler->ensembleN*sampler->ensembleSize, sampler->maxDim, &MCMC_fisher_wrapper_v2,   sampler->userParameters,  100,sampler);
	}
	else if(generation_method.find("EA_IMRPhenomD_NRT") != std::string::npos ){
		std::vector<std::vector<int>> blocks = std::vector<std::vector<int>>(4);
		for(int i = 0 ; i<7; i++){
			blocks[0].push_back(i);
		}
		for(int i = 7 ; i<sampler->maxDim; i++){
			blocks[1].push_back(i);
		}
		for(int i = 0 ; i<sampler->maxDim; i++){
			blocks[2].push_back(i);
		}
		blocks[3].push_back(7);
		blocks[3].push_back(12);
		std::vector<double> blockProb = {.25,.25,.25,.25};
		//std::vector<std::vector<int>> blocks = {
		//				{7,8,9,10}};
		//std::vector<double> blockProb = {1};
		propArray[3] = new bayesship::blockFisherProposal(sampler->ensembleN*sampler->ensembleSize, sampler->minDim, &MCMC_fisher_wrapper_v3,   sampler->userParameters,  100,sampler,blocks, blockProb );

	}
	else{
		std::vector<std::vector<int>> blocks = std::vector<std::vector<int>>(3);
		for(int i = 0 ; i<7; i++){
			blocks[0].push_back(i);
		}
		for(int i = 7 ; i<sampler->maxDim; i++){
			blocks[1].push_back(i);
		}
		for(int i = 0 ; i<sampler->maxDim; i++){
			blocks[2].push_back(i);
		}
		std::vector<double> blockProb = {.3,.3,.4};
		//std::vector<std::vector<int>> blocks = {
		//				{7,8,9,10}};
		//std::vector<double> blockProb = {1};
		propArray[3] = new bayesship::blockFisherProposal(sampler->ensembleN*sampler->ensembleSize, sampler->minDim, &MCMC_fisher_wrapper_v3,   sampler->userParameters,  100,sampler,blocks, blockProb );
	}

	//################################################
	//std::vector<std::vector<int>> blocks = {
	//				{0,1,2,3,4,5,6},
	//				{7,8,9,10}};
	//std::vector<double> blockProb = {.5,.5};
	////std::vector<std::vector<int>> blocks = {
	////				{7,8,9,10}};
	////std::vector<double> blockProb = {1};
	//propArray[4] = new bayesship::blockFisherProposal(sampler->ensembleN*sampler->ensembleSize, sampler->minDim, &MCMC_fisher_wrapper_v3,   sampler->userParameters,  100,sampler,blocks, blockProb );

	//################################################

	//Rough estimate of the temperatures
	double betaTemp[sampler->ensembleSize];
	betaTemp[0] = 1;
	betaTemp[sampler->ensembleSize-1] = 0;
	double deltaBeta = pow((1e2),1./sampler->ensembleSize);
	for(int i = 1 ; i<sampler->ensembleSize-1; i++){
		betaTemp[i] = betaTemp[i-1]/deltaBeta;
	}

	double **propProb = new double*[chainN];	
	for(int i = 0 ; i<chainN; i++){
		int ensemble = i / sampler->ensembleN;
		propProb[i] = new double[proposalFnN];	
		propProb[i][2] = 0.0;

		//propProb[i][0] = 0.05;
		//propProb[i][1] = 0.25;
		//propProb[i][3] = 0.7;
		propProb[i][1] = .4 - 0.15*( betaTemp[ensemble]); //.25 to .7
		propProb[i][3] = 0.05 + .55*( betaTemp[ensemble]); //.7 to .15
		propProb[i][4] = 0.1 + .05*( betaTemp[ensemble]); //.7 to .15
		//propProb[i][4] = 0.1 + .2*( betaTemp[ensemble]); //.3 to .1


		//propProb[i][1] = 0; //.25 to .7
		//propProb[i][3] = 0.3 + .45*( betaTemp[ensemble]); //.7 to .25

		//propProb[i][0] = 1.- propProb[i][4] - propProb[i][3] - propProb[i][1]- propProb[i][2];
		propProb[i][0] = 0.05;
		
		double sum = 0 ;
		for(int j = 0 ; j<proposalFnN; j++){
			sum+=propProb[i][j];	
		}
		for(int j = 0 ; j<proposalFnN; j++){
			propProb[i][j]/=sum;
		}
		for(int j = 0 ; j<proposalFnN; j++){
			std::cout<<propProb[i][j]<<", ";
		}
		std::cout<<"\n";
	
		//std::cout<<propProb[i][0]<<" "<<propProb[i][1]<<" "<<propProb[i][3]<<std::endl;

	}
	
	//propProb[0] = 0.05;
	//propProb[1] = 0.3;
	//propProb[2] = 0.05;
	//propProb[3] = 0.6;
	//propProb[0] = 0.05;
	//propProb[1] = 0.25;
	//propProb[2] = 0.0;
	//propProb[3] = 0.7;

	//propProb[0] = 0.1;
	//propProb[1] = 0.0;
	//propProb[2] = 0.0;
	//propProb[3] = 0.9;



	bayesship::proposalData *propData = new bayesship::proposalData(chainN,proposalFnN, propArray, (double *)nullptr, propProb);
	
	
	
	
	//int proposalFnN = 4;
	//bayesship::proposalFn propArray[proposalFnN];
	//void *proposalFnVariables[proposalFnN];	
	//propArray[0] = bayesship::gaussianProposal;
	//propArray[1] = bayesship::differentialEvolutionProposal;
	//propArray[2] = bayesship::KDEProposal;
	//propArray[3] = bayesship::FisherProposal;

	////Rough estimate of the temperatures
	//double betaTemp[sampler->ensembleSize];
	//betaTemp[0] = 1;
	//betaTemp[sampler->ensembleSize-1] = 0;
	//double deltaBeta = pow((1e2),1./sampler->ensembleSize);
	//for(int i = 1 ; i<sampler->ensembleSize-1; i++){
	//	betaTemp[i] = betaTemp[i-1]/deltaBeta;
	//}

	//float **propProb = new float*[chainN];	
	//for(int i = 0 ; i<chainN; i++){
	//	int ensemble = i / sampler->ensembleN;
	//	//std::cout<<ensemble<<std::endl;
	//	//std::cout<<betaTemp[ensemble]<<std::endl;
	//	propProb[i] = new float[proposalFnN];	
	//	propProb[i][2] = 0.0;

	//	//propProb[i][0] = 0.05;
	//	//propProb[i][1] = 0.25;
	//	//propProb[i][3] = 0.7;
	//	propProb[i][1] = .7 - 0.45*( betaTemp[ensemble]); //.25 to .7
	//	propProb[i][3] = 0.25 + .45*( betaTemp[ensemble]); //.7 to .25

	//	propProb[i][0] = 1. - propProb[i][3] - propProb[i][1]- propProb[i][2];
	//	//std::cout<<propProb[i][0]<<" "<<propProb[i][1]<<" "<<propProb[i][3]<<std::endl;

	//}
	//
	////propProb[0] = 0.05;
	////propProb[1] = 0.3;
	////propProb[2] = 0.05;
	////propProb[3] = 0.6;
	////propProb[0] = 0.05;
	////propProb[1] = 0.25;
	////propProb[2] = 0.0;
	////propProb[3] = 0.7;

	////propProb[0] = 0.1;
	////propProb[1] = 0.0;
	////propProb[2] = 0.0;
	////propProb[3] = 0.9;


	//bayesship::gaussianProposalVariables *gpv = new bayesship::gaussianProposalVariables(sampler->ensembleN*sampler->ensembleSize, sampler->maxDim);

	////bayesship::KDEProposalVariables *kdepv = new bayesship::KDEProposalVariables(sampler->ensembleN*sampler->ensembleSize, sampler->maxDim);
	//bayesship::KDEProposalVariables kdepv(sampler->ensembleN*sampler->ensembleSize, sampler->maxDim);

	//bayesship::FisherProposalVariables *fpv = new bayesship::FisherProposalVariables(sampler->ensembleN*sampler->ensembleSize, sampler->maxDim, MCMC_fisher_wrapper_v2,   sampler->userParameters,  100);

	//proposalFnVariables[0] = (void *)gpv;
	//proposalFnVariables[1] = (void *)nullptr;
	//proposalFnVariables[2] = (void *)&kdepv;
	////proposalFnVariables[2] = (void *)nullptr;
	//proposalFnVariables[3] = (void *)fpv;

	//bayesship::proposalFnWriteCheckpoint *writeCheckpointFns = new bayesship::proposalFnWriteCheckpoint[proposalFnN];
	//writeCheckpointFns[0] = bayesship::gaussianProposalWriteCheckpoint;
	//writeCheckpointFns[1] = nullptr;
	//writeCheckpointFns[2] = nullptr;
	//writeCheckpointFns[3] = nullptr;

	//bayesship::proposalFnLoadCheckpoint *loadCheckpointFns = new bayesship::proposalFnLoadCheckpoint[proposalFnN];
	//loadCheckpointFns[0] = bayesship::gaussianProposalLoadCheckpoint;
	//loadCheckpointFns[1] = nullptr;
	//loadCheckpointFns[2] = nullptr;
	//loadCheckpointFns[3] = nullptr;



	//bayesship::proposalFnData *propData = new bayesship::proposalFnData(chainN,proposalFnN, propArray,proposalFnVariables, (float *)nullptr, propProb,writeCheckpointFns,loadCheckpointFns);

	//#######################################################################################
	sampler->proposalFns = propData;
	//##########################################################



	sampler->sample();

	for(int i =0 ; i<chainN ; i++){
		delete [] propProb[i];
	}
	delete [] propProb;


	//PTMCMC_MH_dynamic_PT_alloc_uncorrelated(sampler_output,output, dimension, N_steps, chain_N, 
	//	max_chain_N_thermo_ensemble,initial_pos,seeding_var,ensemble_initial_pos, chain_temps, 
	//	swp_freq, t0, nu,max_chunk_size,chain_distribution_scheme,
	//	log_prior,MCMC_likelihood_wrapper, MCMC_fisher_wrapper,(void**)user_parameters,numThreads, pool, 
	//	//log_prior,MCMC_likelihood_wrapper, NULL,(void**)user_parameters,numThreads, pool, 
	//	show_prog,statistics_filename,
	//	chain_filename, likelihood_log_filename,checkpoint_filename);
	
	

	if(!mcmcVar.mcmc_mod_struct->fisher_weights){
		for(int i = 0 ;i<chainN;i++){
			delete[] user_parameters[i]->fisher_weights;
		}
	}
	if(!mcmcVar.mcmc_mod_struct->weights){
		for(int i = 0 ;i<chainN;i++){
			delete[] user_parameters[i]->weights;
		}
	}
	
	//Deallocate fftw plans
	for (int i =0;i<num_detectors;i++)
		deallocate_FFTW_mem(&plans[i]);
	//#################################################
	delete propData->proposals[0];
	delete propData->proposals[1];
	delete propData->proposals[2];
	delete propData->proposals[3];
	//delete propData->proposals[4];
	delete propData;
	delete [] propArray;
	for(int i = 0 ; i<num_detectors; i++){
		delete [] burn_data[i];
		delete [] burn_freqs[i];
		delete [] burn_noise[i];
		deallocate_FFTW_mem(&burn_plans[i]);
	}
	delete [] burn_data;
	delete [] burn_lengths;
	delete [] burn_noise;
	delete [] burn_freqs;
	delete [] burn_plans;
	for(int i = 0 ; i<chainN; i++){
		delete user_parameters[i];
	}
	delete [] user_parameters;
	delete [] mcmcVarVec;
	delete ll;
	//#################################################
	free(plans);
	return sampler;
}

void RJPTMCMC_method_specific_prep_v2(std::string generation_method, int dimension, bool *intrinsic,MCMC_modification_struct *mod_struct)
{
	int totalmod = (mod_struct->gIMR_Nmod_phi + mod_struct->gIMR_Nmod_sigma + mod_struct->gIMR_Nmod_beta + mod_struct->gIMR_Nmod_alpha  + mod_struct->ppE_Nmod);
	if(generation_method.find("EA") != std::string::npos){totalmod+=3;}
	debugger_print(__FILE__,__LINE__,totalmod);
	if(generation_method.find("PhenomD") != std::string::npos && (dimension - totalmod) == 4)
	{
		std::cout<<"Sampling in parameters: ln chirpmass, eta, chi1, chi2";
		for(int i =0; i<totalmod; i++){
			std::cout<<", mod_"<<i;
		}
		std::cout<<std::endl;
		*intrinsic=true;
	} 
	else if(generation_method.find("PhenomD_NRT") != std::string::npos && (dimension - totalmod) == 6)
	{
		std::cout<<"Sampling in parameters: ln chirpmass, eta, chi1, chi2, ln tidal1,  ln tidal2";
		for(int i =0; i<totalmod; i++){
			std::cout<<", mod_"<<i;
		}
		std::cout<<std::endl;
		*intrinsic=true;
	} 
	else if((generation_method.find("PhenomPv2") != std::string::npos || generation_method.find("PhenomPv3") != std::string::npos)
		&& (dimension - totalmod) == 8)
	{
		std::cout<<"Sampling in parameters: ln chirpmass, eta, a1, a2, tilt1, tilt2, phi1, phi2";
		for(int i =0; i<totalmod; i++){
			std::cout<<", mod_"<<i;
		}
		std::cout<<std::endl;
		*intrinsic=true;
	} 
	else if(generation_method.find("PhenomD") != std::string::npos && (dimension - totalmod) == 11)
	{
		std::cout<<"Sampling in parameters: RA, sin DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, chi1, chi2"<<std::endl;
		for(int i =0; i<totalmod; i++){
			std::cout<<", mod_"<<i;
		}
		std::cout<<std::endl;
		*intrinsic=false;
	} 
	else if(generation_method.find("PhenomD_NRT") != std::string::npos && (dimension - totalmod) == 13)
	{
		std::cout<<"Sampling in parameters: RA, sin  DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, chi1, chi2, ln tidal1, ln tidal2"<<std::endl;
		for(int i =0; i<totalmod; i++){
			std::cout<<", mod_"<<i;
		}
		std::cout<<std::endl;
		*intrinsic=false;
	} 
	else if(generation_method.find("PhenomD_NRT") != std::string::npos && (dimension - totalmod) == 12)
	{
		std::cout<<"Sampling in parameters: RA, sin  DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, chi1, chi2, ln tidal_s"<<std::endl;
		for(int i =0; i<totalmod; i++){
			std::cout<<", mod_"<<i;
		}
		std::cout<<std::endl;
		mcmc_intrinsic=false;
	} 
	else if((generation_method.find("PhenomPv2") != std::string::npos || generation_method.find("PhenomPv3") != std::string::npos)
		&& (dimension - totalmod) == 15)
	{
		std::cout<<"Sampling in parameters: RA, sin DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, a1, a2,cos tilt1, cos tilt2, phi1, phi2"<<std::endl;
		for(int i =0; i<totalmod; i++){
			std::cout<<", mod_"<<i;
		}
		std::cout<<std::endl;
		*intrinsic=false;
	} 
	else{
		std::cout<<
			"Input parameters not valid, please check that input is compatible with the supported methods - dimension combinations"<<std::endl;
		exit(1);
	}
}

/*! \brief Unpacks MCMC parameters for method specific initiation 
 *
 * Populates seeding vector if non supplied, populates mcmc_Nmod, populates mcmc_log_beta, populates mcmc_intrinsic
 */
void PTMCMC_method_specific_prep_v2(std::string generation_method, int dimension, bool *intrinsic,MCMC_modification_struct *mod_struct)
{
	int totalmod = (mod_struct->gIMR_Nmod_phi + mod_struct->gIMR_Nmod_sigma + mod_struct->gIMR_Nmod_beta + mod_struct->gIMR_Nmod_alpha  + mod_struct->ppE_Nmod);
	if(generation_method.find("EA") != std::string::npos){totalmod+=3;}
	if(generation_method.find("EA") != std::string::npos){totalmod+=2;} // two dissipative tidal dissipation numbers
	debugger_print(__FILE__,__LINE__,totalmod);
	if(generation_method.find("PhenomD") != std::string::npos && (dimension - totalmod) == 4)
	{
		std::cout<<"Sampling in parameters: ln chirpmass, eta, chi1, chi2";
		for(int i =0; i<totalmod; i++){
			std::cout<<", mod_"<<i;
		}
		std::cout<<std::endl;
		*intrinsic=true;
	} 
	else if(generation_method.find("PhenomD_NRT") != std::string::npos && (dimension - totalmod) == 6)
	{
		std::cout<<"Sampling in parameters: ln chirpmass, eta, chi1, chi2, ln tidal1,  ln tidal2";
		for(int i =0; i<totalmod; i++){
			std::cout<<", mod_"<<i;
		}
		std::cout<<std::endl;
		*intrinsic=true;
	} 
	else if(generation_method.find("PhenomD_NRT") != std::string::npos && (dimension - totalmod) == 8)
	{
		std::cout<<"Sampling in parameters: ln chirpmass, eta, chi1, chi2, ln tidal1,  ln tidal2, diss_tidal1, diss_tidal2";
		for(int i =0; i<totalmod; i++){
			std::cout<<", mod_"<<i;
		}
		std::cout<<std::endl;
		*intrinsic=true;
	} 
	else if((generation_method.find("PhenomPv2") != std::string::npos || generation_method.find("PhenomPv3") != std::string::npos)
		&& (dimension - totalmod) == 8)
	{
		std::cout<<"Sampling in parameters: ln chirpmass, eta, a1, a2, tilt1, tilt2, phi1, phi2";
		for(int i =0; i<totalmod; i++){
			std::cout<<", mod_"<<i;
		}
		std::cout<<std::endl;
		*intrinsic=true;
	} 
	else if(generation_method.find("PhenomD") != std::string::npos && (dimension - totalmod) == 11)
	{
		std::cout<<"Sampling in parameters: RA, sin DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, chi1, chi2"<<std::endl;
		for(int i =0; i<totalmod; i++){
			std::cout<<", mod_"<<i;
		}
		std::cout<<std::endl;
		*intrinsic=false;
	} 
	else if(generation_method.find("PhenomD_NRT") != std::string::npos && (dimension - totalmod) == 13)
	{
		std::cout<<"Sampling in parameters: RA, sin  DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, chi1, chi2, ln tidal1, ln tidal2"<<std::endl;
		for(int i =0; i<totalmod; i++){
			std::cout<<", mod_"<<i;
		}
		std::cout<<std::endl;
		*intrinsic=false;
	} 
	else if(generation_method.find("PhenomD_NRT") != std::string::npos && (dimension - totalmod) == 15)
	{
		std::cout<<"Sampling in parameters: RA, sin  DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, chi1, chi2, ln tidal1, ln tidal2, diss_tidal1, diss_tidal2"<<std::endl;
		for(int i =0; i<totalmod; i++){
			std::cout<<", mod_"<<i;
		}
		std::cout<<std::endl;
		*intrinsic=false;
	} 
	else if(generation_method.find("PhenomD_NRT") != std::string::npos && (dimension - totalmod) == 12)
	{
		std::cout<<"Sampling in parameters: RA, sin  DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, chi1, chi2, ln tidal_s"<<std::endl;
		for(int i =0; i<totalmod; i++){
			std::cout<<", mod_"<<i;
		}
		std::cout<<std::endl;
		mcmc_intrinsic=false;
	} 
	else if(generation_method.find("PhenomD_NRT") != std::string::npos && (dimension - totalmod) == 14)
	{
		std::cout<<"Sampling in parameters: RA, sin  DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, chi1, chi2, ln tidal_s, diss_tidal_1, diss_tidal2"<<std::endl;
		for(int i =0; i<totalmod; i++){
			std::cout<<", mod_"<<i;
		}
		std::cout<<std::endl;
		mcmc_intrinsic=false;
	} 
	else if((generation_method.find("PhenomPv2") != std::string::npos || generation_method.find("PhenomPv3") != std::string::npos)
		&& (dimension - totalmod) == 15)
	{
		std::cout<<"Sampling in parameters: RA, sin DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, a1, a2,cos tilt1, cos tilt2, phi1, phi2"<<std::endl;
		for(int i =0; i<totalmod; i++){
			std::cout<<", mod_"<<i;
		}
		std::cout<<std::endl;
		*intrinsic=false;
	} 
	else{
		std::cout<<
			"Input parameters not valid, please check that input is compatible with the supported methods - dimension combinations"<<std::endl;
		exit(1);
	}
}



/*! \brief utility to do MCMC specific transformations on the input param vector before passing to the repacking utillity
 *
 * Returns the local generation method to be used in the LL functions
 */
std::string MCMC_prep_params_v2(double *param, double *temp_params, gen_params_base<double> *gen_params, int dimension, std::string generation_method, MCMC_modification_struct *mod_struct, bool intrinsic, double gmst)
{
	if(intrinsic) gen_params->sky_average = true;
	else gen_params->sky_average = false;
	gen_params->f_ref = 20;
	gen_params->shift_time = true;
	gen_params->shift_phase = true;
	//gen_params->shift_time = false;
	//gen_params->shift_phase = false;
	gen_params->gmst = gmst;
	gen_params->equatorial_orientation=false;
	gen_params->horizon_coord=false;

	gen_params->tidal_love = mod_struct->tidal_love;
	gen_params->tidal_love_error = mod_struct->tidal_love_error;
	gen_params->alpha_param = mod_struct->alpha_param;
	gen_params->EA_region1 = mod_struct->EA_region1; 

	gen_params->NSflag1 = mod_struct->NSflag1;
	gen_params->NSflag2 = mod_struct->NSflag2;
	//gen_params->NSflag1 = false;
	//gen_params->NSflag2 = false;
	for(int i = 0 ; i <dimension; i++){
		temp_params[i]=param[i];
	}
	int base = dimension;
	if(check_mod(generation_method)){
		//if(generation_method.find("ppE") != std::string::npos ||
		//	generation_method.find("dCS") !=std::string::npos||
		//	generation_method.find("EdGB") != std::string::npos){
		//	gen_params->bppe=mcmc_mod_struct->bppe;
		//	gen_params->Nmod=mcmc_mod_struct->ppE_Nmod;
		//	gen_params->betappe=new double[gen_params->Nmod];
		//}
		if(generation_method.find("ppE") != std::string::npos ||
			check_theory_support(generation_method)){
			gen_params->bppe=mod_struct->bppe;
			gen_params->Nmod=mod_struct->ppE_Nmod;
			gen_params->betappe=new double[gen_params->Nmod];
			base = dimension - mod_struct->ppE_Nmod;
		}
		else if(generation_method.find("gIMR") != std::string::npos){
			gen_params->Nmod_phi=mod_struct->gIMR_Nmod_phi;
			gen_params->phii=mod_struct->gIMR_phii;
			if(gen_params->Nmod_phi !=0){
				gen_params->delta_phi=new double[gen_params->Nmod_phi];
			}
			gen_params->Nmod_sigma=mod_struct->gIMR_Nmod_sigma;
			gen_params->sigmai=mod_struct->gIMR_sigmai;
			if(gen_params->Nmod_sigma !=0){
				gen_params->delta_sigma=new double[gen_params->Nmod_sigma];
			}
			gen_params->Nmod_beta=mod_struct->gIMR_Nmod_beta;
			gen_params->betai=mod_struct->gIMR_betai;
			if(gen_params->Nmod_beta !=0){
				gen_params->delta_beta=new double[gen_params->Nmod_beta];
			}
			gen_params->Nmod_alpha=mod_struct->gIMR_Nmod_alpha;
			gen_params->alphai=mod_struct->gIMR_alphai;
			if(gen_params->Nmod_alpha !=0){
				gen_params->delta_alpha=new double[gen_params->Nmod_alpha];
			}
			base = dimension 
				- mod_struct->gIMR_Nmod_phi 
				- mod_struct->gIMR_Nmod_sigma 
				- mod_struct->gIMR_Nmod_beta 
				- mod_struct->gIMR_Nmod_alpha; 
		}
		if((generation_method.find("dCS")!= std::string::npos ||
			generation_method.find("EdGB")!=std::string::npos)){
			//temp_params[base] = pow(temp_params[base],.25)/(c*1000);
			temp_params[base] = 
				pow_int(temp_params[base]/(c/1000.) , 4);
		}
	}
	return generation_method;
}




void MCMC_fisher_transformations_v2(
	double *param, 
	double **fisher, 
	int dimension,
	std::string generation_method,
	bool intrinsic,
	MCMC_modification_struct *mod_struct
	)
{
	if(!intrinsic){
		fisher[0][0] += 1./(4*M_PI*M_PI);//RA
		fisher[1][1] += 1./4;//sin DEC
		fisher[2][2] += 1./(4*M_PI*M_PI);//psi
		fisher[3][3] += 1./(4);//cos iota
		fisher[4][4] += 1./(4*M_PI*M_PI);//phiref
		fisher[5][5] += 1./(.01);//tc
		fisher[8][8] += 1./.25;//eta
		fisher[9][9] += 1./4;//spin1
		fisher[10][10] += 1./4;//spin2
		if(generation_method.find("PhenomPv2") != std::string::npos || generation_method.find("PhenomPv3") != std::string::npos){
			fisher[11][11] += 1./4;//cos theta1
			fisher[12][12] += 1./4;//cos theta2
			fisher[13][13] += 1./(4*M_PI*M_PI);//phi1
			fisher[14][14] += 1./(4*M_PI*M_PI);//phi2
		}
	}
	else{
		if(generation_method.find("PhenomPv2") != std::string::npos || generation_method.find("PhenomPv3") != std::string::npos){
			fisher[1][1] =1./(.25) ;//eta
			fisher[2][2] =1./(4);//spin1
			fisher[3][3] =1./(4);//spin2
			fisher[4][4] =1./(4);//cos theta1
			fisher[5][5] =1./(4);//cos theta2
			fisher[6][6] =1./(4*M_PI*M_PI) ;//phi1
			fisher[7][7] =1./(4*M_PI*M_PI) ;//phi2
		}
		else if (generation_method.find("PhenomD")!=std::string::npos){
			fisher[1][1] =1./(.25) ;//eta
			fisher[2][2] =1./(4) ;//spin1
			fisher[3][3] =1./(4) ;//spin2
	
		}
	}

	if(generation_method.find("dCS") != std::string::npos ||
		generation_method.find("EdGB") != std::string::npos){
		int base = dimension - mod_struct->ppE_Nmod;
		for(int i = 0 ; i<dimension; i++){
			//Transform to root alpha from alpha^2
			double factor = 4* pow(param[base], 3./4.);
			//Transform to km from sec
			factor *= 1000 / c ;
			fisher[base][i] *= factor ;
			fisher[i][base] *= factor;
		}
	}

	if(generation_method.find("EA") != std::string::npos){
		fisher[dimension-3][dimension-3] += 1./pow(10, -2.);
		fisher[dimension-2][dimension-2] += 1./pow(10, -4.);
		fisher[dimension-1][dimension-1] += 1./pow(10, -1.);
		//std::cout<<fisher[7][12]<<std::endl;
		//fisher[7][12] = -1.15341564e+05;
		//fisher[12][7] = -1.15341564e+05;
	}

	//if(generation_method.find("EA") != std::string::npos){
	//  //for(int i = 0 ; i <4; i++){
	//    for(int i = 0 ; i <3; i++){
	//      for(int j = 0 ; j<dimension; j++){
	//	if(i!=j){
	//	  fisher[dimension-1-i][dimension-1-j] = 0;
	//	  fisher[dimension-1-j][dimension-1-i] = 0;
	//	}
	//	/*
	//	  if(i==j){
	//	  fisher[dimension-1-i][dimension-1-j] = 1./pow(10, -6.);
	//	  }*/
	//      }
	//    }
	//    fisher[dimension-3][dimension-3] = 1./pow(10, -2.);
	//    fisher[dimension-2][dimension-2] = 1./pow(10, -4.);
	//    fisher[dimension-1][dimension-1] = 1./pow(10, -2.);
	//}
	return;

}

void MCMC_fisher_wrapper_RJ_v3(bayesship::positionInfo *pos,   double **output, std::vector<int> block, void *userParameters)
{
	mcmcVariablesRJ *mcmcVarRJ= (mcmcVariablesRJ *)userParameters;
	

//##########################################################
//##########################################################
	mcmcVariables mcmcVar ;
	mcmcVar.mcmc_noise = mcmcVarRJ->mcmc_noise;
	//mcmcVar.mcmc_init_pos = initial_pos;
	mcmcVar.mcmc_frequencies = mcmcVarRJ->mcmc_frequencies;
	mcmcVar.mcmc_data = mcmcVarRJ->mcmc_data;
	mcmcVar.mcmc_data_length = mcmcVarRJ->mcmc_data_length;
	mcmcVar.mcmc_detectors = mcmcVarRJ->mcmc_detectors;
	mcmcVar.mcmc_generation_method = mcmcVarRJ->mcmc_generation_method;
	mcmcVar.mcmc_fftw_plans = mcmcVarRJ->mcmc_fftw_plans;
	mcmcVar.mcmc_num_detectors = mcmcVarRJ->mcmc_num_detectors;
	mcmcVar.mcmc_gps_time = mcmcVarRJ->mcmc_gps_time;
	mcmcVar.mcmc_gmst = gps_to_GMST_radian(mcmcVarRJ->mcmc_gps_time);
	mcmcVar.mcmc_mod_struct = mcmcVarRJ->mcmc_mod_struct;
	mcmcVar.mcmc_save_waveform = true;
	mcmcVar.maxDim = mcmcVarRJ->minDim;
	mcmcVar.user_parameters = mcmcVarRJ->user_parameters;




	MCMC_fisher_wrapper_v3(pos,   output, block, (void*)&mcmcVar);

	return ;
}

void MCMC_fisher_wrapper_RJ_v2(bayesship::positionInfo *pos,   double **output, std::vector<int> block, void *userParameters)
{
	mcmcVariablesRJ *mcmcVarRJ= (mcmcVariablesRJ *)userParameters;
	

//##########################################################
//##########################################################
	mcmcVariables mcmcVar ;
	mcmcVar.mcmc_noise = mcmcVarRJ->mcmc_noise;
	//mcmcVar.mcmc_init_pos = initial_pos;
	mcmcVar.mcmc_frequencies = mcmcVarRJ->mcmc_frequencies;
	mcmcVar.mcmc_data = mcmcVarRJ->mcmc_data;
	mcmcVar.mcmc_data_length = mcmcVarRJ->mcmc_data_length;
	mcmcVar.mcmc_detectors = mcmcVarRJ->mcmc_detectors;
	mcmcVar.mcmc_generation_method = mcmcVarRJ->mcmc_generation_method;
	mcmcVar.mcmc_fftw_plans = mcmcVarRJ->mcmc_fftw_plans;
	mcmcVar.mcmc_num_detectors = mcmcVarRJ->mcmc_num_detectors;
	mcmcVar.mcmc_gps_time = mcmcVarRJ->mcmc_gps_time;
	mcmcVar.mcmc_gmst = gps_to_GMST_radian(mcmcVarRJ->mcmc_gps_time);
	mcmcVar.mcmc_mod_struct = mcmcVarRJ->mcmc_mod_struct;
	mcmcVar.mcmc_save_waveform = true;
	mcmcVar.maxDim = mcmcVarRJ->minDim;
	mcmcVar.user_parameters = mcmcVarRJ->user_parameters;




	MCMC_fisher_wrapper_v2(pos,   output, (void*)&mcmcVar);

	return ;
}

int invertFisherBlock(double **fisherIn, double **fisherOut, int dimIn, std::vector<int> ids)
{
	int dimOut = ids.size();
	double **covIn = allocate_2D_array(dimIn,dimIn);
	double **covOut = allocate_2D_array(dimOut,dimOut);

	int status = gsl_cholesky_matrix_invert(fisherIn, covIn, dimIn);
	
	if(status == 0){
		for(int i = 0 ; i<dimOut; i++){
			for(int j = 0 ; j<dimOut; j++){
				covOut[i][j] = covIn[ids[i]][ids[j]];
			}
		}
		status = gsl_cholesky_matrix_invert(covOut, fisherOut, dimOut);
	}

	deallocate_2D_array(covIn, dimIn, dimIn);
	deallocate_2D_array(covOut, dimOut, dimOut);

	return status;
}

void MCMC_fisher_wrapper_v3(bayesship::positionInfo *pos,   double **output, std::vector<int> ids, void *userParameters)
{
	mcmcVariables *mcmcVar= (mcmcVariables *)userParameters;


	int dimension = mcmcVar->maxDim;
	double **tempOutput= allocate_2D_array(dimension, dimension);
	double **tempCov= allocate_2D_array(dimension, dimension);
	double *temp_params = new double[dimension];
	double param[dimension];
	for(int i = 0 ; i<dimension; i++){	
		param[i] = pos->parameters[i];
	}
	//#########################################################################
	gen_params_base<double> params;
	std::string local_gen = MCMC_prep_params_v2(param, 
		temp_params,&params, dimension, mcmcVar->mcmc_generation_method,mcmcVar->mcmc_mod_struct, mcmcVar->mcmc_intrinsic,mcmcVar->mcmc_gmst );
	//#########################################################################
	//#########################################################################
	repack_parameters(temp_params, &params, 
		"MCMC_"+mcmcVar->mcmc_generation_method, dimension, NULL);
	//std::cout<<temp_params[11]<<" "<<temp_params[12]<<std::endl;
	//repack_parameters(mcmc_init_pos, &params, 
	//	"MCMC_"+mcmc_generation_method, dimension, NULL);
	//#########################################################################
	//#########################################################################
	//std::cout<<"INCL angle fisher: "<<params.incl_angle<<std::endl;
	for(int j =0; j<dimension; j++){
		for(int k =0; k<dimension; k++)
		{
			tempOutput[j][k] =0;
		}
	} 
	double **local_freq=mcmcVar->mcmc_frequencies;
	double **local_noise=mcmcVar->mcmc_noise;
	double **local_weights=mcmcVar->user_parameters->weights;
	int *local_lengths= mcmcVar->mcmc_data_length;
	std::string local_integration_method = "SIMPSONS";
	bool local_log10F = mcmcVar->user_parameters->fisher_log10F;
	if(mcmcVar->user_parameters->fisher_freq){
		local_freq = mcmcVar->user_parameters->fisher_freq;
	}
	if(mcmcVar->user_parameters->fisher_PSD){
		local_noise = mcmcVar->user_parameters->fisher_PSD;
	}
	if(mcmcVar->user_parameters->fisher_length){
		local_lengths = mcmcVar->user_parameters->fisher_length;
	}
	if(mcmcVar->user_parameters->fisher_weights){
		local_weights = mcmcVar->user_parameters->fisher_weights;
	}
	if(mcmcVar->user_parameters->fisher_GAUSS_QUAD){
		local_integration_method = "GAUSSLEG";
	}

	std::string local_gen_method = mcmcVar->mcmc_generation_method;
	int local_dimension = dimension;  
	//if(local_gen_method.find("EA") != std::string::npos && ids.size() != 2)
	//  {
	//    	local_gen_method = "IMRPhenomD_NRT";
	//    	local_dimension -= 3;
	//  }
	double **temp_out = allocate_2D_array(local_dimension,local_dimension);
	//double **temp_out = allocate_2D_array(dimension,dimension);
	for(int i =0 ; i <mcmcVar->mcmc_num_detectors; i++){
		
		//Use AD 
		if(mcmcVar->user_parameters->fisher_AD)
		{	
			std::unique_lock<std::mutex> lock{*(mcmcVar->user_parameters->mFish)};
			//fisher_autodiff(mcmcVar->mcmc_frequencies[i], mcmcVar->mcmc_data_length[i],
			//	"MCMC_"+mcmcVar->mcmc_generation_method, mcmcVar->mcmc_detectors[i],mcmcVar->mcmc_detectors[0],temp_out,dimension, 
			//	(gen_params *)(&params),  "SIMPSONS",(double *)NULL,false,mcmcVar->mcmc_noise[i]);
			fisher_autodiff(local_freq[i], local_lengths[i],
				"MCMC_"+local_gen_method, mcmcVar->mcmc_detectors[i],mcmcVar->mcmc_detectors[0],temp_out,local_dimension, 
				(gen_params *)(&params),  local_integration_method,local_weights[i],true,local_noise[i]);
		}
		else{
			fisher_numerical(local_freq[i], local_lengths[i],
				"MCMC_"+local_gen_method, mcmcVar->mcmc_detectors[i],mcmcVar->mcmc_detectors[0],temp_out,local_dimension, 
				&params, 4, NULL, NULL, local_noise[i],
				mcmcVar->user_parameters->QuadMethod);

		}
		for(int j =0; j<local_dimension; j++){
			for(int k =0; k<local_dimension; k++)
			{
				tempOutput[j][k] +=temp_out[j][k];
				//if(std::isnan(output[j][k]))
				//{
				//      std::cout<<j<<" "<<k<<" "<<temp_out[j][k]<<std::endl;
				//}
			}
		} 
	}
	//Add prior information to fisher
	//if(mcmcVar->mcmc_generation_method.find("Pv2") && !mcmcVar->mcmc_intrinsic){
	
	//if(local_gen_method.find("EA") != std::string::npos && ids.size() == 2){
	if(false){
		MCMC_fisher_transformations_v2(temp_params, tempOutput,dimension,"IMRPhenomD_NRT",mcmcVar->mcmc_intrinsic,
			mcmcVar->mcmc_mod_struct);
	    	tempOutput[dimension-3][dimension-3] = 1./pow(10, -2.);
	}
	else{
		MCMC_fisher_transformations_v2(temp_params, tempOutput,dimension,local_gen,mcmcVar->mcmc_intrinsic,
			mcmcVar->mcmc_mod_struct);
	}
	deallocate_2D_array(temp_out, local_dimension,local_dimension);

	//Try marginalizing over other parameters, otherwise just use subfisher without marginalizing
	int status = invertFisherBlock(tempOutput, output, dimension, ids);
	if(status == 1){
		for(int i = 0 ; i<ids.size(); i++){
			for(int j = 0 ; j<ids.size(); j++){
				output[i][j] = tempOutput[ids[i]][ids[j]];	
			}
		}

	}
	//if(local_gen_method.find("EA") != std::string::npos && ids.size() == 2){
	//	for(int j =0; j<ids.size(); j++){
	//		for(int k =0; k<ids.size(); k++)
	//		{
	//			std::cout<<output[j][k]<<", ";
	//		}
	//		std::cout<<std::endl;
	//	} 
	//	std::cout<<std::endl;
	//}
	//for(int i = 0 ; i<ids.size(); i++){
	//	for(int j = 0 ; j<ids.size(); j++){
	//		std::cout<<output[i][j] <<", ";	
	//	}
	//	std::cout<<std::endl;
	//}

	//if(ids.size() == dimension){
	//	for(int i = 0 ; i<dimension; i++){
	//		for(int j = 0 ; j<dimension; j++){
	//			output[i][j] = tempOutput[i][j];	
	//		}
	//	}
	//}
	//else if(
	//std::find(ids.begin(), ids.end(), 0) != ids.end() && 
	//std::find(ids.begin(), ids.end(), 1) != ids.end() && 
	//std::find(ids.begin(), ids.end(), 2) != ids.end() &&  
	//std::find(ids.begin(), ids.end(), 3) != ids.end() &&
	//std::find(ids.begin(), ids.end(), 4) != ids.end() &&
	//std::find(ids.begin(), ids.end(), 5) != ids.end() &&
	//std::find(ids.begin(), ids.end(), 6) != ids.end() 
	//){
	//	for(int i = 0 ; i<7; i++){
	//		for(int j = 0 ; j<7; j++){
	//			output[i][j] = tempOutput[i][j];	
	//		}
	//	}
	//}
	//else if(
	//std::find(ids.begin(), ids.end(), 7) != ids.end() && 
	//std::find(ids.begin(), ids.end(), 8) != ids.end() && 
	//std::find(ids.begin(), ids.end(), 9) != ids.end() && 
	//std::find(ids.begin(), ids.end(), 10) != ids.end() 
	//){
	//	for(int i = 0 ; i<ids.size(); i++){
	//		for(int j = 0 ; j<ids.size(); j++){
	//			output[i][j] = tempOutput[i+7][j+7];	
	//		}
	//	}
	//}
	//////////////////////////////////////////////
	//if(!interface->burn_phase)
	//{
	//	debugger_print(__FILE__,__LINE__,"Fisher MCMC");
	//	double **cov = allocate_2D_array( dimension,dimension);
	//	gsl_cholesky_matrix_invert(output, cov, dimension);
	//	for(int i = 0 ; i<dimension; i++){
	//		std::cout<<sqrt(cov[i][i])<<" ";	
	//		//for(int j = 0 ; j<dimension; j++){
	//		//	std::cout<<cov[i][j]<<" ";	
	//		//}
	//		//std::cout<<std::endl;	
	//		
	//	}
	//	std::cout<<std::endl;	
	//	deallocate_2D_array(cov, dimension,dimension);
	//}
	//////////////////////////////////////////////

	//Cleanup
	delete [] temp_params;
	deallocate_2D_array(tempOutput, dimension, dimension);
	deallocate_2D_array(tempCov, dimension, dimension);
	if(check_mod(local_gen)){
		//if(local_gen.find("ppE") != std::string::npos ||
		//	local_gen.find("dCS")!=std::string::npos||
		//	local_gen.find("EdGB")!=std::string::npos){
		//	delete [] params.betappe;
		//}
		if(local_gen.find("ppE") != std::string::npos ||
			check_theory_support(local_gen)){
			delete [] params.betappe;
		}
		else if( local_gen.find("gIMR") != std::string::npos){
			if(mcmcVar->mcmc_mod_struct ->gIMR_Nmod_phi !=0){
				delete [] params.delta_phi;
			}
			if(mcmcVar->mcmc_mod_struct ->gIMR_Nmod_sigma !=0){
				delete [] params.delta_sigma;
			}
			if(mcmcVar->mcmc_mod_struct ->gIMR_Nmod_beta !=0){
				delete [] params.delta_beta;
			}
			if(mcmcVar->mcmc_mod_struct ->gIMR_Nmod_alpha !=0){
				delete [] params.delta_alpha;
			}

		}
	}


}

void MCMC_fisher_wrapper_v2(bayesship::positionInfo *pos,   double **output, void *userParameters)
{
	mcmcVariables *mcmcVar= (mcmcVariables *)userParameters;

	int dimension = mcmcVar->maxDim;
	double *temp_params = new double[dimension];
	double param[dimension];
	for(int i = 0 ; i<dimension; i++){	
		param[i] = pos->parameters[i];
	}
	//#########################################################################
	gen_params_base<double> params;
	std::string local_gen = MCMC_prep_params_v2(param, 
		temp_params,&params, dimension, mcmcVar->mcmc_generation_method,mcmcVar->mcmc_mod_struct, mcmcVar->mcmc_intrinsic,mcmcVar->mcmc_gmst );
	//#########################################################################
	//#########################################################################
	repack_parameters(temp_params, &params, 
		"MCMC_"+mcmcVar->mcmc_generation_method, dimension, NULL);
	//std::cout<<temp_params[11]<<" "<<temp_params[12]<<std::endl;
	//repack_parameters(mcmc_init_pos, &params, 
	//	"MCMC_"+mcmc_generation_method, dimension, NULL);
	//#########################################################################
	//#########################################################################
	//std::cout<<"INCL angle fisher: "<<params.incl_angle<<std::endl;
	for(int j =0; j<dimension; j++){
		for(int k =0; k<dimension; k++)
		{
			output[j][k] =0;
		}
	} 
	double **local_freq=mcmcVar->mcmc_frequencies;
	double **local_noise=mcmcVar->mcmc_noise;
	double **local_weights=mcmcVar->user_parameters->weights;
	int *local_lengths= mcmcVar->mcmc_data_length;
	std::string local_integration_method = "SIMPSONS";
	bool local_log10F = mcmcVar->user_parameters->fisher_log10F;
	if(mcmcVar->user_parameters->fisher_freq){
		local_freq = mcmcVar->user_parameters->fisher_freq;
	}
	if(mcmcVar->user_parameters->fisher_PSD){
		local_noise = mcmcVar->user_parameters->fisher_PSD;
	}
	if(mcmcVar->user_parameters->fisher_length){
		local_lengths = mcmcVar->user_parameters->fisher_length;
	}
	if(mcmcVar->user_parameters->fisher_weights){
		local_weights = mcmcVar->user_parameters->fisher_weights;
	}
	if(mcmcVar->user_parameters->fisher_GAUSS_QUAD){
		local_integration_method = "GAUSSLEG";
	}
	double **temp_out = allocate_2D_array(dimension,dimension);
	for(int i =0 ; i <mcmcVar->mcmc_num_detectors; i++){
		
		//Use AD 
		if(mcmcVar->user_parameters->fisher_AD)
		{	
			std::unique_lock<std::mutex> lock{*(mcmcVar->user_parameters->mFish)};
			//fisher_autodiff(mcmcVar->mcmc_frequencies[i], mcmcVar->mcmc_data_length[i],
			//	"MCMC_"+mcmcVar->mcmc_generation_method, mcmcVar->mcmc_detectors[i],mcmcVar->mcmc_detectors[0],temp_out,dimension, 
			//	(gen_params *)(&params),  "SIMPSONS",(double *)NULL,false,mcmcVar->mcmc_noise[i]);
			fisher_autodiff(local_freq[i], local_lengths[i],
				"MCMC_"+mcmcVar->mcmc_generation_method, mcmcVar->mcmc_detectors[i],mcmcVar->mcmc_detectors[0],temp_out,dimension, 
				(gen_params *)(&params),  local_integration_method,local_weights[i],true,local_noise[i]);
		}
		else{
			fisher_numerical(local_freq[i], local_lengths[i],
				"MCMC_"+mcmcVar->mcmc_generation_method, mcmcVar->mcmc_detectors[i],mcmcVar->mcmc_detectors[0],temp_out,dimension, 
				&params, 4, NULL, NULL, local_noise[i]);

		}
		for(int j =0; j<dimension; j++){
			for(int k =0; k<dimension; k++)
			{
				output[j][k] +=temp_out[j][k];
				//if(std::isnan(output[j][k]))
				//{
				//      std::cout<<j<<" "<<k<<" "<<temp_out[j][k]<<std::endl;
				//}
			}
		} 
	}
	//Add prior information to fisher
	//if(mcmcVar->mcmc_generation_method.find("Pv2") && !mcmcVar->mcmc_intrinsic){
	
	MCMC_fisher_transformations_v2(temp_params, output,dimension,local_gen,mcmcVar->mcmc_intrinsic,
		mcmcVar->mcmc_mod_struct);
	deallocate_2D_array(temp_out, dimension,dimension);
	//////////////////////////////////////////////
	//if(!interface->burn_phase)
	//{
	//	debugger_print(__FILE__,__LINE__,"Fisher MCMC");
	//	double **cov = allocate_2D_array( dimension,dimension);
	//	gsl_cholesky_matrix_invert(output, cov, dimension);
	//	for(int i = 0 ; i<dimension; i++){
	//		std::cout<<sqrt(cov[i][i])<<" ";	
	//		//for(int j = 0 ; j<dimension; j++){
	//		//	std::cout<<cov[i][j]<<" ";	
	//		//}
	//		//std::cout<<std::endl;	
	//		
	//	}
	//	std::cout<<std::endl;	
	//	deallocate_2D_array(cov, dimension,dimension);
	//}
	//////////////////////////////////////////////

	//Cleanup
	delete [] temp_params;
	if(check_mod(local_gen)){
		//if(local_gen.find("ppE") != std::string::npos ||
		//	local_gen.find("dCS")!=std::string::npos||
		//	local_gen.find("EdGB")!=std::string::npos){
		//	delete [] params.betappe;
		//}
		if(local_gen.find("ppE") != std::string::npos ||
			check_theory_support(local_gen)){
			delete [] params.betappe;
		}
		else if( local_gen.find("gIMR") != std::string::npos){
			if(mcmcVar->mcmc_mod_struct ->gIMR_Nmod_phi !=0){
				delete [] params.delta_phi;
			}
			if(mcmcVar->mcmc_mod_struct ->gIMR_Nmod_sigma !=0){
				delete [] params.delta_sigma;
			}
			if(mcmcVar->mcmc_mod_struct ->gIMR_Nmod_beta !=0){
				delete [] params.delta_beta;
			}
			if(mcmcVar->mcmc_mod_struct ->gIMR_Nmod_alpha !=0){
				delete [] params.delta_alpha;
			}

		}
	}


}
