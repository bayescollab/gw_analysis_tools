#include "mcmc_gw.h"

#include <bayesship/bayesshipSampler.h>
#include <bayesship/dataUtilities.h>
#include <bayesship/proposalFunctions.h>
#include <bayesship/utilities.h>

#include "detector_util.h"
#include "fisher.h"
#include "gw_fisher.h"
#include "io_util.h"
#include "ortho_basis.h"
#include "ppE_utilities.h"
#include "util.h"
#include "waveform_generator.h"
#include "waveform_util.h"

void MCMC_fisher_wrapper(bayesship::positionInfo* pos, double** output,
                         void* userParameters);

void MCMC_fisher_wrapper_explicit_marginalization(bayesship::positionInfo* pos,
                                                  double** output,
                                                  std::vector<int> ids,
                                                  void* userParameters);

int invertFisherBlock(double** fisherIn, double** fisherOut, int dimIn,
                      std::vector<int> ids);

void MCMC_fisher_wrapper_RJ(bayesship::positionInfo* pos, double** output,
                            std::vector<int> block, void* userParameters);

void MCMC_fisher_GWParameterSpace_wrapper(bayesship::positionInfo* pos,
                                          double** fisherM,
                                          void* userParameters);

void MCMC_fisher_GWParameterSpace_wrapper_explicit_marginalization(
    bayesship::positionInfo* pos, double** fisherM, std::vector<int> ids,
    void* userParameters);

void pack_local_mod_structure(bayesship::bayesshipSampler* sampler,
                              double* param, int* status,
                              std::string waveform_extended, void* parameters,
                              MCMC_modification_struct* full_struct,
                              MCMC_modification_struct* local_struct) {
  if (has_substring(waveform_extended, "gIMR")) {
    int dimct = 0;
    int dphi_boundary = full_struct->gIMR_Nmod_phi + sampler->minDim;
    int dsigma_boundary = full_struct->gIMR_Nmod_sigma + dphi_boundary;
    int dbeta_boundary = full_struct->gIMR_Nmod_beta + dsigma_boundary;
    int dalpha_boundary = full_struct->gIMR_Nmod_alpha + dbeta_boundary;
    for (int i = 0; i < sampler->maxDim; i++) {
      if (status[i] == 1) {
        dimct++;
      }
      if (i >= sampler->minDim) {
        if (status[i] == 1 && i < dphi_boundary) {
          local_struct->gIMR_Nmod_phi++;
        } else if (status[i] == 1 && i < dsigma_boundary) {
          local_struct->gIMR_Nmod_sigma++;
        } else if (status[i] == 1 && i < dbeta_boundary) {
          local_struct->gIMR_Nmod_beta++;
        } else if (status[i] == 1 && i < dalpha_boundary) {
          local_struct->gIMR_Nmod_alpha++;
        }
      }
    }
    if (dimct != sampler->minDim) {
      if (local_struct->gIMR_Nmod_phi != 0) {
        local_struct->gIMR_phii = new int[local_struct->gIMR_Nmod_phi];
      }
      if (local_struct->gIMR_Nmod_sigma != 0) {
        local_struct->gIMR_sigmai = new int[local_struct->gIMR_Nmod_sigma];
      }
      if (local_struct->gIMR_Nmod_beta != 0) {
        local_struct->gIMR_betai = new int[local_struct->gIMR_Nmod_beta];
      }
      if (local_struct->gIMR_Nmod_alpha != 0) {
        local_struct->gIMR_alphai = new int[local_struct->gIMR_Nmod_alpha];
      }

      int ct_phi = 0;
      int ct_sigma = 0;
      int ct_beta = 0;
      int ct_alpha = 0;
      for (int i = sampler->minDim; i < sampler->maxDim; i++) {
        if (status[i] == 1) {
          if (i < dphi_boundary) {
            local_struct->gIMR_phii[ct_phi] =
                full_struct
                    ->gIMR_phii[i - dphi_boundary + full_struct->gIMR_Nmod_phi];
            ct_phi++;
          } else if (i < dsigma_boundary) {
            local_struct->gIMR_sigmai[ct_sigma] =
                full_struct->gIMR_sigmai[i - dsigma_boundary +
                                         full_struct->gIMR_Nmod_sigma];
            ct_sigma++;
          } else if (i < dbeta_boundary) {
            local_struct->gIMR_betai[ct_beta] =
                full_struct->gIMR_betai[i - dbeta_boundary +
                                        full_struct->gIMR_Nmod_beta];
            ct_beta++;
          } else if (i < dalpha_boundary) {
            local_struct->gIMR_alphai[ct_alpha] =
                full_struct->gIMR_alphai[i - dalpha_boundary +
                                         full_struct->gIMR_Nmod_alpha];
            ct_alpha++;
          }
        }
      }
    }
  } else if (has_substring(waveform_extended, "ppE")) {
    int dimct = 0;
    local_struct->ppE_Nmod = 0;
    for (int i = 0; i < sampler->maxDim; i++) {
      if (status[i] == 1) {
        dimct++;
      }

      if (i >= sampler->minDim) {
        if (status[i] == 1) {
          local_struct->ppE_Nmod++;
        }
      }
    }
    if (dimct != sampler->minDim) {
      local_struct->bppe = new double[local_struct->ppE_Nmod];

      int ct_ppe = 0;
      for (int i = sampler->minDim; i < sampler->maxDim; i++) {
        if (status[i] == 1) {
          local_struct->bppe[ct_ppe] = full_struct->bppe[i - sampler->minDim];
          ct_ppe++;
        }
      }
    }
  }

  return;
}

class MCMC_likelihood_wrapper_RJ : public bayesship::probabilityFn {
 public:
  bayesship::bayesshipSampler* sampler;
  mcmcVariablesRJ* mcmcVarRJ;
  virtual double eval(bayesship::positionInfo* pos, int chainID) {
    // return 2;
    // mcmcVariables *mcmcVar = (mcmcVariables *)userParameters;
    // MCMC_user_param *user_param = (MCMC_user_param *)userParameters;

    int maxDim = sampler->maxDim;
    int minDim = sampler->minDim;
    double ll = 0;
    double* temp_params = new double[maxDim];

    int dim = pos->countActiveDimensions();

    int ct = 0;
    for (int i = 0; i < maxDim; i++) {
      if (pos->status[i]) {
        temp_params[ct] = pos->parameters[i];
        // std::cout<<temp_params[ct]<<",";
        ct += 1;
      }
    }
    // std::cout<<std::endl;
    bool wfExtended = false;
    std::string gen_meth = mcmcVarRJ->mcmc_generation_method;
    if (dim > minDim) {
      wfExtended = true;
      gen_meth = mcmcVarRJ->mcmc_generation_method_extended;
    }
    MCMC_modification_struct mod_struct_local;

    if (wfExtended) {
      pack_local_mod_structure(sampler, pos->parameters, pos->status,
                               mcmcVarRJ->mcmc_generation_method_extended,
                               mcmcVarRJ->user_parameters,
                               mcmcVarRJ->mcmc_mod_struct, &mod_struct_local);
    }
    // #########################################################################
    gen_params_base<double> gen_params;
    std::string local_gen = MCMC_prep_params(
        temp_params, temp_params, &gen_params, dim, gen_meth, &mod_struct_local,
        mcmcVarRJ->mcmc_intrinsic, mcmcVarRJ->mcmc_gmst);
    // #########################################################################
    // #########################################################################

    // repack_non_parameters(temp_params, &gen_params,
    //"MCMC_"+mcmc_generation_method, dimension, NULL);
    repack_parameters(temp_params, &gen_params, "MCMC_" + gen_meth, dim);
    // if(gen_params.Nmod !=0){
    //	std::cout<<gen_params.Nmod<<" "<<local_gen<<std::endl;
    //	for(int i = 0 ; i<gen_params.Nmod; i++){
    //		std::cout<<" "<<gen_params.bppe[i]<<" "<<gen_params.betappe[i];
    //	}
    //	std::cout<<std::endl;
    //
    // }
    // #########################################################################
    // #########################################################################
    // return 1;
    std::complex<double>** local_data = mcmcVarRJ->mcmc_data;
    double** local_freqs = mcmcVarRJ->mcmc_frequencies;
    double** local_noise = mcmcVarRJ->mcmc_noise;
    double** local_weights = (mcmcVarRJ->user_parameters)->weights;
    int* local_lengths = mcmcVarRJ->mcmc_data_length;
    fftw_outline* local_plans = mcmcVarRJ->mcmc_fftw_plans;
    std::string local_integration_method = "SIMPSONS";
    // if(interface->burn_phase && user_param->burn_data){
    if (mcmcVarRJ->user_parameters->GAUSS_QUAD) {
      local_integration_method = "GAUSSLEG";
    }
    if (mcmcVarRJ->mcmc_intrinsic) {
      if (has_substring(local_gen, "IMRPhenomD")) {
        if (!mcmcVarRJ->mcmc_save_waveform) {
          for (int i = 0; i < mcmcVarRJ->mcmc_num_detectors; i++) {
            gen_params.theta = 0;
            gen_params.phi = 0;
            gen_params.psi = 0;
            gen_params.phiRef = 1;
            gen_params.f_ref = 10;
            gen_params.incl_angle = 0;
            gen_params.tc = 1;
            std::complex<double>* response = (std::complex<double>*)malloc(
                sizeof(std::complex<double>) * local_lengths[i]);
            fourier_detector_response_horizon(
                local_freqs[i], local_lengths[i], response,
                mcmcVarRJ->mcmc_detectors[i], local_gen, &gen_params);
            ll += maximized_Log_Likelihood_aligned_spin_internal(
                local_data[i], local_noise[i], local_freqs[i], response,
                (size_t)local_lengths[i], &local_plans[i]);
            free(response);
          }
        } else {
          gen_params.theta = 0;
          gen_params.phi = 0;
          gen_params.psi = 0;
          gen_params.phiRef = 1;
          gen_params.f_ref = 10;
          gen_params.incl_angle = 0;
          gen_params.tc = 1;
          std::complex<double>* response = (std::complex<double>*)malloc(
              sizeof(std::complex<double>) * local_lengths[0]);
          fourier_detector_response_horizon(
              local_freqs[0], local_lengths[0], response,
              mcmcVarRJ->mcmc_detectors[0], local_gen, &gen_params);
          // std::complex<double> *hc =
          //	(std::complex<double> *) malloc(sizeof(std::complex<double>) *
          // mcmc_data_length[0]); std::complex<double> *hp =
          //	(std::complex<double> *) malloc(sizeof(std::complex<double>) *
          // mcmc_data_length[0]); fourier_waveform(mcmc_frequencies[0],
          // mcmc_data_length[0], hp,hc, local_gen, &gen_params);
          for (int i = 0; i < mcmcVarRJ->mcmc_num_detectors; i++) {
            ll += maximized_Log_Likelihood_aligned_spin_internal(
                local_data[i], local_noise[i], local_freqs[i], response,
                (size_t)local_lengths[i], &local_plans[i]);
            // ll +=
            // maximized_Log_Likelihood_unaligned_spin_internal(mcmc_data[i],
            //		mcmc_noise[i],
            //		mcmc_frequencies[i],
            //		hp,
            //		hc,
            //		(size_t) mcmc_data_length[i],
            //		&mcmc_fftw_plans[i]
            //		);
          }
          free(response);
          // free(hp);
          // free(hc);
        }

      } else if (has_substring(local_gen, "IMRPhenomP")) {
        gen_params.theta = 0;
        gen_params.phi = 0;
        gen_params.psi = 0;
        gen_params.phiRef = 1;
        gen_params.f_ref = 20;
        gen_params.incl_angle = 0;
        gen_params.tc = 1;
        waveform_polarizations<double> wp;
        assign_polarizations(local_gen, &wp);
        wp.allocate_memory(local_lengths[0]);
        fourier_waveform(local_freqs[0], local_lengths[0], &wp, local_gen,
                         &gen_params);
        for (int i = 0; i < mcmcVarRJ->mcmc_num_detectors; i++) {
          ll += maximized_Log_Likelihood_unaligned_spin_internal(
              local_data[i], local_noise[i], local_freqs[i], wp.hplus,
              wp.hcross, (size_t)local_lengths[i], &local_plans[i]);
        }
        wp.deallocate_memory();
      }
    } else {
      double RA = gen_params.RA;
      double DEC = gen_params.DEC;
      double PSI = gen_params.psi;
      // if(mcmc_generation_method.find("IMRPhenomD") != std::string:npos){

      ll = MCMC_likelihood_extrinsic(
          mcmcVarRJ->mcmc_save_waveform, &gen_params, local_gen, local_lengths,
          local_freqs, local_data, local_noise, local_weights,
          local_integration_method, mcmcVarRJ->user_parameters->log10F,
          mcmcVarRJ->mcmc_detectors, mcmcVarRJ->mcmc_num_detectors);
      // ll=2;

      //}
      // else if(has_substring(mcmc_generation_method, "IMRPhenomP")){

      //}
    }
    // Cleanup
    delete[] temp_params;
    if (check_mod(local_gen)) {
      // if( has_substring(local_gen, "ppE") ||
      //	has_substring(local_gen, "dCS") ||
      //	has_substring(local_gen, "EdGB")){
      //	delete [] gen_params.betappe;
      // }
      if (has_substring(local_gen, "ppE") || check_theory_support(local_gen)) {
        delete[] gen_params.betappe;
        delete[] mod_struct_local.bppe;
      } else if (has_substring(local_gen, "gIMR")) {
        if (mod_struct_local.gIMR_Nmod_phi != 0) {
          delete[] gen_params.delta_phi;
          delete[] mod_struct_local.gIMR_phii;
        }
        if (mod_struct_local.gIMR_Nmod_sigma != 0) {
          delete[] gen_params.delta_sigma;
          delete[] mod_struct_local.gIMR_sigmai;
        }
        if (mod_struct_local.gIMR_Nmod_beta != 0) {
          delete[] gen_params.delta_beta;
          delete[] mod_struct_local.gIMR_betai;
        }
        if (mod_struct_local.gIMR_Nmod_alpha != 0) {
          delete[] gen_params.delta_alpha;
          delete[] mod_struct_local.gIMR_alphai;
        }
      }
    }
    // debugger_print(__FILE__,__LINE__,ll);
    return ll;
  }
};

class ppEFisherRJVariables {
 public:
  double* bppe = nullptr;
  int fisherDim;
  double** fisher = nullptr;
  int detectN;
  double** psds = nullptr;
  double** freqs = nullptr;
  int length;
  ppEFisherRJVariables(double* bppe, double** psds, double** freqs, int detectN,
                       int fisherDim, int length) {
    this->detectN = detectN;
    this->fisherDim = fisherDim;
    this->length = length;
    this->bppe = new double[fisherDim];
    for (int i = 0; i < fisherDim; i++) {
      this->bppe[i] = bppe[i];
    }
    this->psds = new double*[detectN];
    for (int i = 0; i < detectN; i++) {
      this->psds[i] = new double[length];
      for (int j = 0; j < length; j++) {
        this->psds[i][j] = psds[i][j];
      }
    }
    this->freqs = new double*[detectN];
    for (int i = 0; i < detectN; i++) {
      this->freqs[i] = new double[length];
      for (int j = 0; j < length; j++) {
        this->freqs[i][j] = freqs[i][j];
      }
    }

    this->fisher = new double*[fisherDim];
    for (int i = 0; i < fisherDim; i++) {
      this->fisher[i] = new double[fisherDim];
      for (int j = 0; j < fisherDim; j++) {
        this->fisher[i][j] = 0;
      }
    }

    // Calculate the actual fisher
    double* integrand = new double[length];
    for (int l = 0; l < detectN; l++) {
      for (int i = 0; i < fisherDim; i++) {
        for (int j = 0; j <= i; j++) {
          for (int k = 0; k < length; k++) {
            integrand[k] =
                (pow(freqs[l][k], bppe[i] / 3. + bppe[j] / 3. - 7 / 3.) /
                 psds[l][k]);
            integrand[k] *=
                (2. / 15.) * pow(M_PI, -4. / 3. + bppe[i] / 3. + bppe[j] / 3.);
          }
          fisher[i][j] +=
              simpsons_sum(freqs[l][1] - freqs[l][0], length, integrand);
        }
      }
    }
    for (int i = 0; i < fisherDim; i++) {
      for (int j = 0; j <= i; j++) {
        fisher[j][i] = fisher[i][j];
      }
    }
    // std::cout<<"Fisher: "<<std::endl;
    // for(int i = 0 ; i < fisherDim; i++){
    //	for(int j = 0 ; j <fisherDim; j++){
    //		std::cout<<fisher[i][j] <<" , ";
    //	}
    //	std::cout<<std::endl;
    // }

    delete[] integrand;
  };
  ~ppEFisherRJVariables() {
    if (bppe) {
      delete[] bppe;
      bppe = nullptr;
    }
    if (fisher) {
      for (int i = 0; i < fisherDim; i++) {
        delete[] fisher[i];
      }
      delete[] fisher;
      fisher = nullptr;
    }
    if (psds) {
      for (int i = 0; i < detectN; i++) {
        delete[] psds[i];
      }
      delete[] psds;
      psds = nullptr;
    }
    if (freqs) {
      for (int i = 0; i < detectN; i++) {
        delete[] freqs[i];
      }
      delete[] freqs;
      freqs = nullptr;
    }
  };
};

void MCMC_fisher_wrapper_RJ_ppE(bayesship::positionInfo* pos, double** output,
                                std::vector<int> block, void* userParameters) {
  double chirp = std::exp(pos->parameters[7]) * MSOL_SEC;
  double DL = std::exp(pos->parameters[6]) * MPC_SEC;
  ppEFisherRJVariables* p = (ppEFisherRJVariables*)userParameters;
  for (int i = 0; i < p->fisherDim; i++) {
    for (int j = 0; j < p->fisherDim; j++) {
      output[i][j] =
          p->fisher[i][j] *
          pow(chirp, 4. - 7. / 3. + p->bppe[i] / 3. + p->bppe[j] / 3.) / DL /
          DL;
    }
  }
  // std::cout<<"Fisher: "<<std::endl;
  // for(int i = 0 ; i < p->fisherDim; i++){
  //	for(int j = 0 ; j <p->fisherDim; j++){
  //		std::cout<<output[i][j] <<" , ";
  //	}
  //	std::cout<<std::endl;
  // }
  return;
}

class MCMC_likelihood_wrapper : public bayesship::probabilityFn {
 public:
  bayesship::bayesshipSampler* sampler;
  mcmcVariables* mcmcVar;
  virtual double eval(bayesship::positionInfo* pos, int chainID) {
    int dimension = sampler->maxDim;
    double ll = 0;
    gen_params_base<double> gen_params;

    // Run modernized code if param_space is active
    // TODO: Eventually make this THE way to compute likelihoods
    if (mcmcVar->param_space) {
      mcmcVar->param_space->to_gen_params(pos->parameters, gen_params);
      return mcmcVar->likelihood->log_likelihood(
          &gen_params, mcmcVar->param_space->generation_method());
    }

    double* temp_params = new double[dimension];
    std::string local_gen = MCMC_prep_params(
        pos->parameters, temp_params, &gen_params, dimension,
        mcmcVar->mcmc_generation_method, mcmcVar->mcmc_mod_struct,
        mcmcVar->mcmc_intrinsic, mcmcVar->mcmc_gmst);
    repack_parameters(temp_params, &gen_params,
                      "MCMC_" + mcmcVar->mcmc_generation_method, dimension);

    std::complex<double>** local_data = mcmcVar->mcmc_data;
    double** local_freqs = mcmcVar->mcmc_frequencies;
    double** local_noise = mcmcVar->mcmc_noise;
    double** local_weights = (mcmcVar->user_parameters)->weights;
    int* local_lengths = mcmcVar->mcmc_data_length;
    fftw_outline* local_plans = mcmcVar->mcmc_fftw_plans;
    std::string local_integration_method = "SIMPSONS";
    if (mcmcVar->user_parameters->GAUSS_QUAD) {
      local_integration_method = "GAUSSLEG";
    }
    if (mcmcVar->mcmc_intrinsic) {
      if (mcmcVar->mcmc_generation_method.find("IMRPhenomD") !=
          std::string::npos) {
        if (!mcmcVar->mcmc_save_waveform) {
          for (int i = 0; i < mcmcVar->mcmc_num_detectors; i++) {
            gen_params.theta = 0;
            gen_params.phi = 0;
            gen_params.psi = 0;
            gen_params.phiRef = 1;
            gen_params.f_ref = 10;
            gen_params.incl_angle = 0;
            gen_params.tc = 1;
            std::complex<double>* response = (std::complex<double>*)malloc(
                sizeof(std::complex<double>) * local_lengths[i]);
            fourier_detector_response_horizon(
                local_freqs[i], local_lengths[i], response,
                mcmcVar->mcmc_detectors[i], local_gen, &gen_params);
            ll += maximized_Log_Likelihood_aligned_spin_internal(
                local_data[i], local_noise[i], local_freqs[i], response,
                (size_t)local_lengths[i], &local_plans[i]);
            free(response);
          }
        } else {
          gen_params.theta = 0;
          gen_params.phi = 0;
          gen_params.psi = 0;
          gen_params.phiRef = 1;
          gen_params.f_ref = 10;
          gen_params.incl_angle = 0;
          gen_params.tc = 1;
          std::complex<double>* response = (std::complex<double>*)malloc(
              sizeof(std::complex<double>) * local_lengths[0]);
          fourier_detector_response_horizon(
              local_freqs[0], local_lengths[0], response,
              mcmcVar->mcmc_detectors[0], local_gen, &gen_params);
          for (int i = 0; i < mcmcVar->mcmc_num_detectors; i++) {
            ll += maximized_Log_Likelihood_aligned_spin_internal(
                local_data[i], local_noise[i], local_freqs[i], response,
                (size_t)local_lengths[i], &local_plans[i]);
          }
          free(response);
        }

      } else if (mcmcVar->mcmc_generation_method.find("IMRPhenomP") !=
                 std::string::npos) {
        gen_params.theta = 0;
        gen_params.phi = 0;
        gen_params.psi = 0;
        gen_params.phiRef = 1;
        gen_params.f_ref = 20;
        gen_params.incl_angle = 0;
        gen_params.tc = 1;
        waveform_polarizations<double> wp;
        assign_polarizations(mcmcVar->mcmc_generation_method, &wp);
        wp.allocate_memory(local_lengths[0]);
        fourier_waveform(local_freqs[0], local_lengths[0], &wp, local_gen,
                         &gen_params);
        for (int i = 0; i < mcmcVar->mcmc_num_detectors; i++) {
          ll += maximized_Log_Likelihood_unaligned_spin_internal(
              local_data[i], local_noise[i], local_freqs[i], wp.hplus,
              wp.hcross, (size_t)local_lengths[i], &local_plans[i]);
        }
        wp.deallocate_memory();
      }
    } else {
      double RA = gen_params.RA;
      double DEC = gen_params.DEC;
      double PSI = gen_params.psi;

      ll = MCMC_likelihood_extrinsic(
          mcmcVar->mcmc_save_waveform, &gen_params, local_gen, local_lengths,
          local_freqs, local_data, local_noise, local_weights,
          local_integration_method, mcmcVar->user_parameters->log10F,
          mcmcVar->mcmc_detectors, mcmcVar->mcmc_num_detectors);
    }
    // Cleanup
    delete[] temp_params;
    if (check_mod(local_gen)) {
      if (has_substring(local_gen, "ppE") || check_theory_support(local_gen)) {
        delete[] gen_params.betappe;
      } else if (has_substring(local_gen, "gIMR")) {
        if (mcmcVar->mcmc_mod_struct->gIMR_Nmod_phi != 0) {
          delete[] gen_params.delta_phi;
        }
        if (mcmcVar->mcmc_mod_struct->gIMR_Nmod_sigma != 0) {
          delete[] gen_params.delta_sigma;
        }
        if (mcmcVar->mcmc_mod_struct->gIMR_Nmod_beta != 0) {
          delete[] gen_params.delta_beta;
        }
        if (mcmcVar->mcmc_mod_struct->gIMR_Nmod_alpha != 0) {
          delete[] gen_params.delta_alpha;
        }
      }
    }
    if (isnan(ll)) {
      return bayesship::limitInf;
    }
    // std::cout<<ll<<std::endl;
    return ll;
  }
};

/*! \brief Takes in an MCMC checkpoint file and continues the chain
 *
 * Obviously, the user must be sure to correctly match the dimension, number of
 * chains, the generation_method, the prior function, the data, psds, freqs, and
 * the detectors (number and name), and the gps_time to the previous run,
 * otherwise the behavior of the sampler is undefined.
 *
 * numThreads and pool do not necessarily have to be the same
 */
bayesship::bayesshipSampler* RJPTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW(
    int minDim, int maxDim, int independentSamples, int ensembleSize,
    int ensembleN, bayesship::positionInfo* initialPosition,
    bayesship::positionInfo** initialEnsemble, double swapProb,
    int burnIterations, int burnPriorIterations, int priorIterations,
    bool writePriorData, int batchSize, double** priorRanges,
    bayesship::probabilityFn* log_prior, int numThreads, bool pool,
    int num_detectors, std::complex<double>** data, double** noise_psd,
    double** frequencies, int* data_length, double gps_time,
    std::string* detectors, MCMC_modification_struct* mod_struct,
    std::string generation_method, std::string generation_method_extended,
    std::string outputDir, std::string outputFileMoniker,
    bool ignoreExistingCheckpoint, bool restrictSwapTemperatures,
    bool coldChainStorageOnly) {
  int chainN = ensembleSize * ensembleN;
  // Create fftw plan for each detector (length of data stream may be different)
  fftw_outline* plans =
      (fftw_outline*)malloc(sizeof(fftw_outline) * num_detectors);
  for (int i = 0; i < num_detectors; i++) {
    allocate_FFTW_mem_forward(&plans[i], data_length[i]);
  }

  std::mutex fisher_mutex;

  // ##########################################################
  // ##########################################################
  mcmcVariablesRJ mcmcVarRJ;
  mcmcVarRJ.mcmc_noise = noise_psd;
  // mcmcVarRJ.mcmc_init_pos = initial_pos;
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

  // ##########################################################
  // ##########################################################

  // mcmc_noise = noise_psd;
  // mcmc_init_pos = initial_pos;
  // mcmc_frequencies = frequencies;
  // mcmc_data = data;
  // mcmc_data_length = data_length;
  // mcmc_detectors = detectors;
  // mcmc_generation_method = generation_method;
  // mcmc_fftw_plans = plans;
  // mcmc_num_detectors = num_detectors;
  // mcmc_gps_time = gps_time;
  // mcmc_gmst = gps_to_GMST_radian(mcmc_gps_time);
  ////mcmc_Nmod = mod_struct->ppE_Nmod;
  ////mcmc_bppe = mod_struct->bppe;
  // mcmc_mod_struct = mod_struct;
  // mcmc_log_beta = false;
  // mcmc_intrinsic = false;

  // To save time, intrinsic waveforms can be saved between detectors, if the
  // frequencies are all the same
  mcmc_save_waveform = true;
  for (int i = 1; i < mcmcVarRJ.mcmc_num_detectors; i++) {
    if (mcmcVarRJ.mcmc_data_length[i] != mcmcVarRJ.mcmc_data_length[0] ||
        mcmcVarRJ.mcmc_frequencies[i][0] != mcmcVarRJ.mcmc_frequencies[0][0] ||
        mcmcVarRJ.mcmc_frequencies[i][mcmcVarRJ.mcmc_data_length[i] - 1] !=
            mcmcVarRJ.mcmc_frequencies[0][mcmcVarRJ.mcmc_data_length[0] - 1]) {
      mcmcVarRJ.mcmc_save_waveform = false;
    }
  }

  // std::cout<<"Base:"<<std::endl;
  RJPTMCMC_method_specific_prep(generation_method, maxDim,
                                &(mcmcVarRJ.mcmc_intrinsic),
                                mcmcVarRJ.mcmc_mod_struct);
  // std::cout<<"Extending up to:"<<std::endl;
  // RJPTMCMC_method_specific_prep(generation_method_extended, maxDim,
  // &(mcmcVarRJ.mcmc_intrinsic), mcmcVarRJ.mcmc_mod_struct);

  // ######################################################
  int T = (int)(1. / (mcmcVarRJ.mcmc_frequencies[0][1] -
                      mcmcVarRJ.mcmc_frequencies[0][0]));
  debugger_print(__FILE__, __LINE__, T);
  int burn_factor = T / 4;  // Take all sources to 4 seconds
  debugger_print(__FILE__, __LINE__, burn_factor);
  std::complex<double>** burn_data =
      new std::complex<double>*[mcmcVarRJ.mcmc_num_detectors];
  double** burn_freqs = new double*[mcmcVarRJ.mcmc_num_detectors];
  double** burn_noise = new double*[mcmcVarRJ.mcmc_num_detectors];
  int* burn_lengths = new int[mcmcVarRJ.mcmc_num_detectors];
  fftw_outline* burn_plans = new fftw_outline[mcmcVarRJ.mcmc_num_detectors];
  for (int j = 0; j < mcmcVarRJ.mcmc_num_detectors; j++) {
    burn_lengths[j] = mcmcVarRJ.mcmc_data_length[j] / burn_factor;
    burn_data[j] = new std::complex<double>[burn_lengths[j]];
    burn_freqs[j] = new double[burn_lengths[j]];
    burn_noise[j] = new double[burn_lengths[j]];
    allocate_FFTW_mem_forward(&burn_plans[j], burn_lengths[j]);
    int ct = 0;
    for (int k = 0; k < mcmcVarRJ.mcmc_data_length[j]; k++) {
      if (k % burn_factor == 0 && ct < burn_lengths[j]) {
        burn_data[j][ct] = mcmcVarRJ.mcmc_data[j][k];
        burn_freqs[j][ct] = mcmcVarRJ.mcmc_frequencies[j][k];
        burn_noise[j][ct] = mcmcVarRJ.mcmc_noise[j][k];
        ct++;
      }
    }
  }

  MCMC_user_param** user_parameters = NULL;
  user_parameters = new MCMC_user_param*[chainN];
  for (int i = 0; i < chainN; i++) {
    user_parameters[i] = new MCMC_user_param;

    user_parameters[i]->burn_data = burn_data;
    user_parameters[i]->burn_freqs = burn_freqs;
    user_parameters[i]->burn_noise = burn_noise;
    user_parameters[i]->burn_lengths = burn_lengths;
    user_parameters[i]->burn_plans = burn_plans;

    user_parameters[i]->mFish = &fisher_mutex;
    user_parameters[i]->GAUSS_QUAD = mod_struct->GAUSS_QUAD;
    user_parameters[i]->log10F = mod_struct->log10F;

    if (mod_struct->weights) {
      user_parameters[i]->weights = mod_struct->weights;
    } else {
      user_parameters[i]->weights = new double*[num_detectors];
      for (int j = 0; j < num_detectors; j++) {
        user_parameters[i]->weights[j] = NULL;
      }
    }

    user_parameters[i]->fisher_GAUSS_QUAD = mod_struct->fisher_GAUSS_QUAD;
    user_parameters[i]->fisher_log10F = mod_struct->fisher_log10F;
    user_parameters[i]->fisher_freq = mod_struct->fisher_freq;
    if (mod_struct->fisher_weights) {
      user_parameters[i]->fisher_weights = mod_struct->fisher_weights;
    } else {
      user_parameters[i]->fisher_weights = new double*[num_detectors];
      for (int j = 0; j < num_detectors; j++) {
        user_parameters[i]->fisher_weights[j] = NULL;
      }
    }
    user_parameters[i]->fisher_PSD = mod_struct->fisher_PSD;
    user_parameters[i]->fisher_length = mod_struct->fisher_length;

    user_parameters[i]->mod_struct = mod_struct;

    // user_parameters[i]->burn_freqs = mcmc_frequencies;
    // user_parameters[i]->burn_data = mcmc_data;
    // user_parameters[i]->burn_noise = mcmc_noise;
    // user_parameters[i]->burn_lengths = mcmc_data_length;
  }
  // mcmcVarRJ.user_parameters = user_parameters;
  // ######################################################
  // ###########################################################
  // double trial[11] = {1,.9,.2,.9,.2,3,7,4,.24,0,0};
  // mcmc_data_interface i;
  // i.min_dim = 11;
  // i.max_dim = 11;
  // i.chain_id  = 0;
  // i.chain_number  = 11;
  // debugger_print(__FILE__,__LINE__,ll);
  // ###########################################################

  MCMC_likelihood_wrapper_RJ* ll = new MCMC_likelihood_wrapper_RJ();
  ll->mcmcVarRJ = &mcmcVarRJ;
  bayesship::bayesshipSampler* sampler =
      new bayesship::bayesshipSampler(ll, log_prior);
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

  sampler->coldOnlyStorage = coldChainStorageOnly;

  mcmcVariablesRJ** mcmcVarVec = new mcmcVariablesRJ*[chainN];
  for (int i = 0; i < chainN; i++) {
    mcmcVarVec[i] = &mcmcVarRJ;
    mcmcVarVec[i]->user_parameters = user_parameters[i];
  }
  sampler->userParameters = (void**)mcmcVarVec;

  int proposalFnN = 6;
  bayesship::proposal** propArray = new bayesship::proposal*[proposalFnN];
  propArray[0] = new bayesship::gaussianProposal(
      sampler->ensembleN * sampler->ensembleSize, sampler->maxDim, sampler);
  // propArray[1] = new bayesship::differentialEvolutionProposal(sampler);
  if (mcmcVarRJ.mcmc_intrinsic) {
    propArray[1] = new bayesship::differentialEvolutionProposal(sampler);
  } else {
    std::vector<std::vector<int>> blocksDiff = std::vector<std::vector<int>>(3);
    for (int i = 0; i < 7; i++) {
      blocksDiff[0].push_back(i);
    }
    for (int i = 7; i < sampler->minDim; i++) {
      blocksDiff[1].push_back(i);
    }
    for (int i = 0; i < sampler->minDim; i++) {
      blocksDiff[2].push_back(i);
    }
    std::vector<double> blocksProbDiff = {0.3, 0.3, .4};
    propArray[1] = new bayesship::blockDifferentialEvolutionProposal(
        sampler, blocksDiff, blocksProbDiff);
  }

  propArray[2] =
      new bayesship::KDEProposal(sampler->ensembleN * sampler->ensembleSize,
                                 sampler->maxDim, sampler, false);

  propArray[3] = new bayesship::randomLayerRJProposal(sampler, .5);

  if (mcmcVarRJ.mcmc_intrinsic) {
    std::vector<std::vector<int>> blocks = std::vector<std::vector<int>>(1);
    for (int i = 0; i < sampler->minDim; i++) {
      blocks[0].push_back(i);
    }
    std::vector<double> blockProb = {1};
    propArray[4] = new bayesship::blockFisherProposal(
        sampler->ensembleN * sampler->ensembleSize, sampler->minDim,
        &MCMC_fisher_wrapper_RJ, sampler->userParameters, 100, sampler, blocks,
        blockProb);
  } else {
    std::vector<std::vector<int>> blocks = std::vector<std::vector<int>>(3);
    for (int i = 0; i < 7; i++) {
      blocks[0].push_back(i);
    }
    for (int i = 7; i < sampler->minDim; i++) {
      blocks[1].push_back(i);
    }
    for (int i = 0; i < sampler->minDim; i++) {
      blocks[2].push_back(i);
    }
    std::vector<double> blockProb = {.3, .3, .4};
    // std::vector<std::vector<int>> blocks = {
    //				{7,8,9,10}};
    // std::vector<double> blockProb = {1};
    propArray[4] = new bayesship::blockFisherProposal(
        sampler->ensembleN * sampler->ensembleSize, sampler->minDim,
        &MCMC_fisher_wrapper_RJ, sampler->userParameters, 100, sampler, blocks,
        blockProb);
  }

  std::vector<std::vector<int>> blocks2 = std::vector<std::vector<int>>(1);
  blocks2[0] = std::vector<int>(maxDim - minDim);
  for (int i = 0; i < maxDim - minDim; i++) {
    blocks2[0][i] = minDim + i;
  }
  std::vector<double> blockProb2 = {1};
  ppEFisherRJVariables* ppEFisherObj = new ppEFisherRJVariables(
      mod_struct->bppe, mcmcVarRJ.mcmc_noise, mcmcVarRJ.mcmc_frequencies,
      mcmcVarRJ.mcmc_num_detectors, maxDim - minDim,
      mcmcVarRJ.mcmc_data_length[0]);
  ppEFisherRJVariables** ppEFisherObjs = new ppEFisherRJVariables*[chainN];
  for (int i = 0; i < chainN; i++) {
    ppEFisherObjs[i] = ppEFisherObj;
  }

  propArray[5] = new bayesship::blockFisherProposal(
      sampler->ensembleN * sampler->ensembleSize, sampler->minDim,
      &MCMC_fisher_wrapper_RJ_ppE, (void**)ppEFisherObjs, 100, sampler, blocks2,
      blockProb2);

  // Rough estimate of the temperatures
  double betaTemp[sampler->ensembleSize];
  betaTemp[0] = 1;
  betaTemp[sampler->ensembleSize - 1] = 0;
  double deltaBeta = pow((1e2), 1. / sampler->ensembleSize);
  for (int i = 1; i < sampler->ensembleSize - 1; i++) {
    betaTemp[i] = betaTemp[i - 1] / deltaBeta;
  }

  double** propProb = new double*[chainN];
  for (int i = 0; i < chainN; i++) {
    int ensemble = i / sampler->ensembleN;
    propProb[i] = new double[proposalFnN];
    propProb[i][2] = 0.0;

    propProb[i][1] = .6 - 0.45 * (betaTemp[ensemble]);  //.15 to .7
    propProb[i][3] = 0.15 - .1 * (betaTemp[ensemble]);  //.05 to .2
    propProb[i][4] = 0.05 + .4 * (betaTemp[ensemble]);  //.45 to .05
    propProb[i][5] = 0.05 + .2 * (betaTemp[ensemble]);  //.25 to .05
    // std::cout<<"TESTING FISHER"<<std::endl;
    // propProb[i][1] = 0.0;
    // propProb[i][3] = 0.0;
    // propProb[i][4] = 0.0;
    // propProb[i][5] = 1.0;

    // propProb[i][0] = 0.05;
    // propProb[i][1] = 0.25;
    // propProb[i][3] = 0.7;
    // propProb[i][3] = 0; //.1 to .2
    // propProb[i][4] = 0;

    // propProb[i][1] = 0; //.25 to .7
    // propProb[i][3] = 0.3 + .45*( betaTemp[ensemble]); //.7 to .25

    // propProb[i][0] = 1.  - propProb[i][1]- propProb[i][2];

    propProb[i][0] = 0.05;

    double sum = 0;
    for (int j = 0; j < proposalFnN; j++) {
      sum += propProb[i][j];
    }
    for (int j = 0; j < proposalFnN; j++) {
      propProb[i][j] /= sum;
    }
    for (int j = 0; j < proposalFnN; j++) {
      std::cout << propProb[i][j] << ", ";
    }
    std::cout << "\n";
  }

  // propProb[0] = 0.05;
  // propProb[1] = 0.3;
  // propProb[2] = 0.05;
  // propProb[3] = 0.6;
  // propProb[0] = 0.05;
  // propProb[1] = 0.25;
  // propProb[2] = 0.0;
  // propProb[3] = 0.7;

  // propProb[0] = 0.1;
  // propProb[1] = 0.0;
  // propProb[2] = 0.0;
  // propProb[3] = 0.9;

  bayesship::proposalData* propData = new bayesship::proposalData(
      chainN, proposalFnN, propArray, (double*)nullptr, propProb);

  sampler->proposalFns = propData;
  sampler->sample();

  // Delete data
  for (int i = 0; i < chainN; i++) {
    delete[] propProb[i];
  }
  delete[] propProb;

  if (!mcmcVarRJ.mcmc_mod_struct->fisher_weights) {
    for (int i = 0; i < chainN; i++) {
      delete[] user_parameters[i]->fisher_weights;
    }
  }
  if (!mcmcVarRJ.mcmc_mod_struct->weights) {
    for (int i = 0; i < chainN; i++) {
      delete[] user_parameters[i]->weights;
    }
  }

  // Deallocate fftw plans
  for (int i = 0; i < num_detectors; i++) deallocate_FFTW_mem(&plans[i]);
  for (int i = 0; i < proposalFnN; i++) {
    delete propData->proposals[i];
  }
  delete[] ppEFisherObjs;
  delete ppEFisherObj;
  delete propData;
  delete[] propArray;
  for (int i = 0; i < num_detectors; i++) {
    delete[] burn_data[i];
    delete[] burn_freqs[i];
    delete[] burn_noise[i];
    deallocate_FFTW_mem(&burn_plans[i]);
  }
  delete[] burn_data;
  delete[] burn_lengths;
  delete[] burn_noise;
  delete[] burn_freqs;
  delete[] burn_plans;
  for (int i = 0; i < chainN; i++) {
    delete user_parameters[i];
  }
  delete[] user_parameters;
  delete[] mcmcVarVec;
  delete ll;
  free(plans);
  return sampler;
}

/*! \brief Takes in an MCMC checkpoint file and continues the chain
 *
 * Obviously, the user must be sure to correctly match the dimension, number of
 * chains, the generation_method, the prior function, the data, psds, freqs, and
 * the detectors (number and name), and the gps_time to the previous run,
 * otherwise the behavior of the sampler is undefined.
 *
 * numThreads and pool do not necessarily have to be the same
 */
bayesship::bayesshipSampler* PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW(
    int dimension, int independentSamples, int ensembleSize, int ensembleN,
    bayesship::positionInfo* initialPosition,
    bayesship::positionInfo** initialEnsemble, double swapProb,
    int burnIterations, int burnPriorIterations, int priorIterations,
    bool writePriorData, int batchSize, double** priorRanges,
    bayesship::probabilityFn* log_prior, int numThreads, bool pool,
    int num_detectors, std::complex<double>** data, double** noise_psd,
    double** frequencies, int* data_length, double gps_time,
    std::string* detectors, MCMC_modification_struct* mod_struct,
    std::string generation_method, std::string outputDir,
    std::string outputFileMoniker, bool ignoreExistingCheckpoint,
    bool restrictSwapTemperatures, bool coldChainStorageOnly) {
  int chainN = ensembleSize * ensembleN;
  // Create fftw plan for each detector (length of data stream may be different)
  fftw_outline* plans =
      (fftw_outline*)malloc(sizeof(fftw_outline) * num_detectors);
  for (int i = 0; i < num_detectors; i++) {
    allocate_FFTW_mem_forward(&plans[i], data_length[i]);
  }

  std::mutex fisher_mutex;
  mcmcVariables mcmcVar;
  mcmcVar.mcmc_noise = noise_psd;
  // mcmcVar.mcmc_init_pos = initial_pos;
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

  if ((mod_struct->param_space == nullptr) !=
      (mod_struct->likelihood == nullptr)) {
    std::cerr << "ERROR: param_space and likelihood must both be set or both "
                 "be null.\n";
    std::exit(EXIT_FAILURE);
  }
  if (mod_struct->param_space != nullptr) {
    mcmcVar.param_space = mod_struct->param_space;
    mcmcVar.likelihood = mod_struct->likelihood;
  }

  // To save time, intrinsic waveforms can be saved between detectors, if the
  // frequencies are all the same
  mcmc_save_waveform = true;
  for (int i = 1; i < mcmcVar.mcmc_num_detectors; i++) {
    if (mcmcVar.mcmc_data_length[i] != mcmcVar.mcmc_data_length[0] ||
        mcmcVar.mcmc_frequencies[i][0] != mcmcVar.mcmc_frequencies[0][0] ||
        mcmcVar.mcmc_frequencies[i][mcmcVar.mcmc_data_length[i] - 1] !=
            mcmcVar.mcmc_frequencies[0][mcmcVar.mcmc_data_length[0] - 1]) {
      mcmcVar.mcmc_save_waveform = false;
    }
  }

  if (mcmcVar.param_space) {
    const auto& names = mcmcVar.param_space->names();
    std::cout << "Sampling in parameters: ";
    for (int i = 0; i < (int)names.size() - 1; i++)
      std::cout << names[i] << ", ";
    std::cout << names.back() << "\n";
    std::cout << "First " << mcmcVar.param_space->number_of_extrinsic_params()
              << " are extrinsic.\n";
  } else {
    PTMCMC_method_specific_prep(generation_method, dimension,
                                &(mcmcVar.mcmc_intrinsic),
                                mcmcVar.mcmc_mod_struct);
  }

  int T = (int)(1. / (mcmcVar.mcmc_frequencies[0][1] -
                      mcmcVar.mcmc_frequencies[0][0]));
  debugger_print(__FILE__, __LINE__, T);
  int burn_factor = T / 4;  // Take all sources to 4 seconds
  debugger_print(__FILE__, __LINE__, burn_factor);
  std::complex<double>** burn_data =
      new std::complex<double>*[mcmcVar.mcmc_num_detectors];
  double** burn_freqs = new double*[mcmcVar.mcmc_num_detectors];
  double** burn_noise = new double*[mcmcVar.mcmc_num_detectors];
  int* burn_lengths = new int[mcmcVar.mcmc_num_detectors];
  fftw_outline* burn_plans = new fftw_outline[mcmcVar.mcmc_num_detectors];
  for (int j = 0; j < mcmcVar.mcmc_num_detectors; j++) {
    burn_lengths[j] = mcmcVar.mcmc_data_length[j] / burn_factor;
    burn_data[j] = new std::complex<double>[burn_lengths[j]];
    burn_freqs[j] = new double[burn_lengths[j]];
    burn_noise[j] = new double[burn_lengths[j]];
    allocate_FFTW_mem_forward(&burn_plans[j], burn_lengths[j]);
    int ct = 0;
    for (int k = 0; k < mcmcVar.mcmc_data_length[j]; k++) {
      if (k % burn_factor == 0 && ct < burn_lengths[j]) {
        burn_data[j][ct] = mcmcVar.mcmc_data[j][k];
        burn_freqs[j][ct] = mcmcVar.mcmc_frequencies[j][k];
        burn_noise[j][ct] = mcmcVar.mcmc_noise[j][k];
        ct++;
      }
    }
  }

  MCMC_user_param** user_parameters = NULL;
  user_parameters = new MCMC_user_param*[chainN];
  for (int i = 0; i < chainN; i++) {
    user_parameters[i] = new MCMC_user_param;

    user_parameters[i]->burn_data = burn_data;
    user_parameters[i]->burn_freqs = burn_freqs;
    user_parameters[i]->burn_noise = burn_noise;
    user_parameters[i]->burn_lengths = burn_lengths;
    user_parameters[i]->burn_plans = burn_plans;

    user_parameters[i]->mFish = &fisher_mutex;
    user_parameters[i]->GAUSS_QUAD = mod_struct->GAUSS_QUAD;
    user_parameters[i]->log10F = mod_struct->log10F;

    if (mod_struct->weights) {
      user_parameters[i]->weights = mod_struct->weights;
    } else {
      user_parameters[i]->weights = new double*[num_detectors];
      for (int j = 0; j < num_detectors; j++) {
        user_parameters[i]->weights[j] = NULL;
      }
    }

    user_parameters[i]->fisher_GAUSS_QUAD = mod_struct->fisher_GAUSS_QUAD;
    user_parameters[i]->fisher_log10F = mod_struct->fisher_log10F;
    user_parameters[i]->fisher_freq = mod_struct->fisher_freq;
    if (mod_struct->fisher_weights) {
      user_parameters[i]->fisher_weights = mod_struct->fisher_weights;
    } else {
      user_parameters[i]->fisher_weights = new double*[num_detectors];
      for (int j = 0; j < num_detectors; j++) {
        user_parameters[i]->fisher_weights[j] = NULL;
      }
    }
    user_parameters[i]->fisher_PSD = mod_struct->fisher_PSD;
    user_parameters[i]->fisher_length = mod_struct->fisher_length;
    user_parameters[i]->QuadMethod = mod_struct->QuadMethod;

    user_parameters[i]->mod_struct = mod_struct;
  }

  MCMC_likelihood_wrapper* ll = new MCMC_likelihood_wrapper();
  ll->mcmcVar = &mcmcVar;
  bayesship::bayesshipSampler* sampler =
      new bayesship::bayesshipSampler(ll, log_prior);
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
  sampler->t0 = mod_struct->t0;
  sampler->nu = mod_struct->nu;
  sampler->Tmax = mod_struct->Tmax;

  sampler->coldOnlyStorage = coldChainStorageOnly;

  mcmcVariables** mcmcVarVec = new mcmcVariables*[chainN];
  for (int i = 0; i < chainN; i++) {
    mcmcVarVec[i] = &mcmcVar;
    mcmcVarVec[i]->user_parameters = user_parameters[i];
  }
  sampler->userParameters = (void**)mcmcVarVec;

  int proposalFnN = 5;
  bayesship::proposal** propArray = new bayesship::proposal*[proposalFnN];
  propArray[0] = new bayesship::gaussianProposal(
      sampler->ensembleN * sampler->ensembleSize, sampler->maxDim, sampler);
  if (mcmcVar.param_space) {
    // Attempting a more flexible block proposal setup based on
    // param_space->number_of_extrinsic_params(). Ideally these settings would
    // be more flexible and outside of PTMCMC; remains to be done.
    const int extrinsic_num = mcmcVar.param_space->number_of_extrinsic_params();
    if (extrinsic_num > 0) {
      std::vector<std::vector<int>> blocksDiff(3);
      std::vector<double> blocksProbDiff = {0.3, 0.3, .4};
      for (int i = 0; i < extrinsic_num; i++) {
        blocksDiff[0].push_back(i);
      }
      for (int i = extrinsic_num; i < sampler->maxDim; i++) {
        blocksDiff[1].push_back(i);
      }
      for (int i = 0; i < sampler->maxDim; i++) {
        blocksDiff[2].push_back(i);
      }
      propArray[1] = new bayesship::blockDifferentialEvolutionProposal(
          sampler, blocksDiff, blocksProbDiff);
      propArray[4] = new bayesship::GMMProposal(
          sampler->ensembleN * sampler->ensembleSize, sampler->maxDim, sampler,
          blocksDiff, blocksProbDiff, 10, 10, 10, 1e-10, false,
          sampler->maxDim * 100);
    } else {
      propArray[1] = new bayesship::differentialEvolutionProposal(sampler);
    }
  } else if (mcmcVar.mcmc_intrinsic) {
    propArray[1] = new bayesship::differentialEvolutionProposal(sampler);
  } else if (has_substring(generation_method, "EA_IMRPhenomD_NRT")) {
    std::vector<std::vector<int>> blocksDiff = std::vector<std::vector<int>>(4);
    for (int i = 0; i < 7; i++) {
      blocksDiff[0].push_back(i);
    }
    for (int i = 7; i < sampler->maxDim; i++) {
      blocksDiff[1].push_back(i);
    }
    for (int i = 0; i < sampler->maxDim; i++) {
      blocksDiff[2].push_back(i);
    }
    blocksDiff[3].push_back(7);
    blocksDiff[3].push_back(12);
    std::vector<double> blocksProbDiff = {0.25, 0.25, .25, .25};
    propArray[1] = new bayesship::blockDifferentialEvolutionProposal(
        sampler, blocksDiff, blocksProbDiff);
    propArray[4] = new bayesship::GMMProposal(
        sampler->ensembleN * sampler->ensembleSize, sampler->maxDim, sampler,
        blocksDiff, blocksProbDiff, 10, 10, 10, 1e-10, false,
        sampler->maxDim * 100);

  } else {
    const int extrinsic_num = 7;
    std::vector<std::vector<int>> blocksDiff(3);
    std::vector<double> blocksProbDiff = {0.3, 0.3, .4};
    for (int i = 0; i < extrinsic_num; i++) {
      blocksDiff[0].push_back(i);
    }
    for (int i = extrinsic_num; i < sampler->maxDim; i++) {
      blocksDiff[1].push_back(i);
    }
    for (int i = 0; i < sampler->maxDim; i++) {
      blocksDiff[2].push_back(i);
    }
    propArray[1] = new bayesship::blockDifferentialEvolutionProposal(
        sampler, blocksDiff, blocksProbDiff);
    propArray[4] = new bayesship::GMMProposal(
        sampler->ensembleN * sampler->ensembleSize, sampler->maxDim, sampler,
        blocksDiff, blocksProbDiff, 10, 10, 10, 1e-10, false,
        sampler->maxDim * 100);
  }
  propArray[2] =
      new bayesship::KDEProposal(sampler->ensembleN * sampler->ensembleSize,
                                 sampler->maxDim, sampler, false);

  if (mcmcVar.param_space) {
    const int extrinsic_num = mcmcVar.param_space->number_of_extrinsic_params();
    if (extrinsic_num > 0) {
      std::vector<std::vector<int>> blocks(3);
      std::vector<double> blockProb{.3, .3, .4};
      for (int i = 0; i < extrinsic_num; i++) {
        blocks[0].push_back(i);
      }
      for (int i = extrinsic_num; i < sampler->maxDim; i++) {
        blocks[1].push_back(i);
      }
      for (int i = 0; i < sampler->maxDim; i++) {
        blocks[2].push_back(i);
      }
      propArray[3] = new bayesship::blockFisherProposal(
          sampler->ensembleN * sampler->ensembleSize, sampler->minDim,
          &MCMC_fisher_GWParameterSpace_wrapper_explicit_marginalization,
          sampler->userParameters, 100, sampler, blocks, blockProb);
    } else {
      propArray[3] = new bayesship::fisherProposal(
          sampler->ensembleN * sampler->ensembleSize, sampler->maxDim,
          &MCMC_fisher_GWParameterSpace_wrapper, sampler->userParameters, 100,
          sampler);
    }

  } else if (mcmcVar.mcmc_intrinsic) {
    propArray[3] = new bayesship::fisherProposal(
        sampler->ensembleN * sampler->ensembleSize, sampler->maxDim,
        &MCMC_fisher_wrapper, sampler->userParameters, 100, sampler);
  } else if (has_substring(generation_method, "EA_IMRPhenomD_NRT")) {
    std::vector<std::vector<int>> blocks = std::vector<std::vector<int>>(4);
    for (int i = 0; i < 7; i++) {
      blocks[0].push_back(i);
    }
    for (int i = 7; i < sampler->maxDim; i++) {
      blocks[1].push_back(i);
    }
    for (int i = 0; i < sampler->maxDim; i++) {
      blocks[2].push_back(i);
    }
    blocks[3].push_back(7);
    blocks[3].push_back(12);
    std::vector<double> blockProb = {.25, .25, .25, .25};
    propArray[3] = new bayesship::blockFisherProposal(
        sampler->ensembleN * sampler->ensembleSize, sampler->minDim,
        &MCMC_fisher_wrapper_explicit_marginalization, sampler->userParameters,
        100, sampler, blocks, blockProb);

  } else {
    const int extrinsic_num = 7;
    std::vector<std::vector<int>> blocks(3);
    std::vector<double> blockProb{.3, .3, .4};

    blocks = std::vector<std::vector<int>>(3);
    for (int i = 0; i < extrinsic_num; i++) {
      blocks[0].push_back(i);
    }
      for (int i = extrinsic_num; i < sampler->maxDim; i++) {
        blocks[1].push_back(i);
      }
      for (int i = 0; i < sampler->maxDim; i++) {
        blocks[2].push_back(i);
      }

    propArray[3] = new bayesship::blockFisherProposal(
        sampler->ensembleN * sampler->ensembleSize, sampler->minDim,
        &MCMC_fisher_wrapper_explicit_marginalization, sampler->userParameters,
        100, sampler, blocks, blockProb);
  }

  // Rough estimate of the temperatures
  double betaTemp[sampler->ensembleSize];
  bayesship::geometric_beta_schedule(betaTemp, sampler->ensembleSize,
                                     sampler->Tmax);

  double** propProb = new double*[chainN];
  for (int i = 0; i < chainN; i++) {
    int ensemble = i / sampler->ensembleN;
    propProb[i] = new double[proposalFnN];
    propProb[i][2] = 0.0;

    propProb[i][1] = .4 - 0.15 * (betaTemp[ensemble]);   //.25 to .7
    propProb[i][3] = 0.05 + .55 * (betaTemp[ensemble]);  //.7 to .15
    propProb[i][4] = 0.1 + .05 * (betaTemp[ensemble]);   //.7 to .15
    propProb[i][0] = 0.05;

    double sum = 0;
    for (int j = 0; j < proposalFnN; j++) {
      sum += propProb[i][j];
    }
    for (int j = 0; j < proposalFnN; j++) {
      propProb[i][j] /= sum;
    }
  }

  bayesship::proposalData* propData = new bayesship::proposalData(
      chainN, proposalFnN, propArray, (double*)nullptr, propProb);
  sampler->proposalFns = propData;
  sampler->sample();

  for (int i = 0; i < chainN; i++) {
    delete[] propProb[i];
  }
  delete[] propProb;

  if (!mcmcVar.mcmc_mod_struct->fisher_weights) {
    for (int i = 0; i < chainN; i++) {
      delete[] user_parameters[i]->fisher_weights;
    }
  }
  if (!mcmcVar.mcmc_mod_struct->weights) {
    for (int i = 0; i < chainN; i++) {
      delete[] user_parameters[i]->weights;
    }
  }

  // Deallocate fftw plans
  for (int i = 0; i < num_detectors; i++) deallocate_FFTW_mem(&plans[i]);
  delete propData->proposals[0];
  delete propData->proposals[1];
  delete propData->proposals[2];
  delete propData->proposals[3];
  delete propData;
  delete[] propArray;
  for (int i = 0; i < num_detectors; i++) {
    delete[] burn_data[i];
    delete[] burn_freqs[i];
    delete[] burn_noise[i];
    deallocate_FFTW_mem(&burn_plans[i]);
  }
  delete[] burn_data;
  delete[] burn_lengths;
  delete[] burn_noise;
  delete[] burn_freqs;
  delete[] burn_plans;
  for (int i = 0; i < chainN; i++) {
    delete user_parameters[i];
  }
  delete[] user_parameters;
  delete[] mcmcVarVec;
  delete ll;
  free(plans);
  return sampler;
}

void RJPTMCMC_method_specific_prep(std::string generation_method, int dimension,
                                   bool* intrinsic,
                                   MCMC_modification_struct* mod_struct) {
  int totalmod = (mod_struct->gIMR_Nmod_phi + mod_struct->gIMR_Nmod_sigma +
                  mod_struct->gIMR_Nmod_beta + mod_struct->gIMR_Nmod_alpha +
                  mod_struct->ppE_Nmod);
  if (has_substring(generation_method, "EA")) {
    totalmod += 3;
  }
  debugger_print(__FILE__, __LINE__, totalmod);
  if (has_substring(generation_method, "PhenomD") &&
      (dimension - totalmod) == 4) {
    std::cout << "Sampling in parameters: ln chirpmass, eta, chi1, chi2";
    for (int i = 0; i < totalmod; i++) {
      std::cout << ", mod_" << i;
    }
    std::cout << std::endl;
    *intrinsic = true;
  } else if (has_substring(generation_method, "PhenomD_NRT") &&
             (dimension - totalmod) == 6) {
    std::cout << "Sampling in parameters: ln chirpmass, eta, chi1, chi2, ln "
                 "tidal1,  ln tidal2";
    for (int i = 0; i < totalmod; i++) {
      std::cout << ", mod_" << i;
    }
    std::cout << std::endl;
    *intrinsic = true;
  } else if ((has_substring(generation_method, "PhenomPv2") ||
              has_substring(generation_method, "PhenomPv3")) &&
             (dimension - totalmod) == 8) {
    std::cout << "Sampling in parameters: ln chirpmass, eta, a1, a2, tilt1, "
                 "tilt2, phi1, phi2";
    for (int i = 0; i < totalmod; i++) {
      std::cout << ", mod_" << i;
    }
    std::cout << std::endl;
    *intrinsic = true;
  } else if (has_substring(generation_method, "PhenomD") &&
             (dimension - totalmod) == 11) {
    std::cout << "Sampling in parameters: RA, sin DEC, psi, cos iota,phi_ref, "
                 "tc,  ln DL, ln chirpmass, eta, chi1, chi2"
              << std::endl;
    for (int i = 0; i < totalmod; i++) {
      std::cout << ", mod_" << i;
    }
    std::cout << std::endl;
    *intrinsic = false;
  } else if (has_substring(generation_method, "PhenomD_NRT")) {
    if ((dimension - totalmod) == 13) {
      std::cout
          << "Sampling in parameters: RA, sin  DEC, psi, cos iota,phi_ref, tc, "
             " ln DL, ln chirpmass, eta, chi1, chi2, ln tidal1, ln tidal2"
          << std::endl;
      for (int i = 0; i < totalmod; i++) {
        std::cout << ", mod_" << i;
      }
      std::cout << std::endl;
      *intrinsic = false;
    }
  } else if ((dimension - totalmod) == 12) {
    std::cout << "Sampling in parameters: RA, sin  DEC, psi, cos iota,phi_ref, "
                 "tc,  ln DL, ln chirpmass, eta, chi1, chi2, ln tidal_s"
              << std::endl;
    for (int i = 0; i < totalmod; i++) {
      std::cout << ", mod_" << i;
    }
    std::cout << std::endl;
    mcmc_intrinsic = false;
  } else if (has_substring(generation_method, "EOS") &&
             (dimension - totalmod) == 14) {
    std::cout
        << "Sampling in parameters: RA, sin  DEC, psi, cos iota, phi_ref, tc,  "
           "ln DL, nbc1, nbc2, chi1, chi2, bump_mag, bump_width, bump_offset"
        << std::endl;
    for (int i = 0; i < totalmod; i++) {
      std::cout << ", mod_" << i;
    }
    std::cout << std::endl;
    *intrinsic = false;
  } else if (has_substring(generation_method, "EOS") &&
             (dimension - totalmod) == 15) {
    std::cout << "Sampling in parameters: RA, sin  DEC, psi, cos iota, "
                 "phi_ref, tc,  ln DL, nbc1, nbc2, chi1, chi2, bump_mag, "
                 "bump_width, bump_offset, plat"
              << std::endl;
    for (int i = 0; i < totalmod; i++) {
      std::cout << ", mod_" << i;
    }
    std::cout << std::endl;
    *intrinsic = false;
  } else if ((has_substring(generation_method, "PhenomPv2") ||
              has_substring(generation_method, "PhenomPv3")) &&
             (dimension - totalmod) == 15) {
    std::cout
        << "Sampling in parameters: RA, sin DEC, psi, cos iota,phi_ref, tc,  "
           "ln DL, ln chirpmass, eta, a1, a2,cos tilt1, cos tilt2, phi1, phi2"
        << std::endl;
    for (int i = 0; i < totalmod; i++) {
      std::cout << ", mod_" << i;
    }
    std::cout << std::endl;
    *intrinsic = false;
  } else {
    std::cout
        << "Input parameters not valid, please check that input is compatible "
           "with the supported methods - dimension combinations"
        << std::endl;
    exit(1);
  }
}

/*! \brief Unpacks MCMC parameters for method specific initiation
 *
 * Populates seeding vector if non supplied, populates mcmc_Nmod, populates
 * mcmc_log_beta, populates mcmc_intrinsic
 */
void PTMCMC_method_specific_prep(std::string generation_method, int dimension,
                                 bool* intrinsic,
                                 MCMC_modification_struct* mod_struct) {
  int totalmod = (mod_struct->gIMR_Nmod_phi + mod_struct->gIMR_Nmod_sigma +
                  mod_struct->gIMR_Nmod_beta + mod_struct->gIMR_Nmod_alpha +
                  mod_struct->ppE_Nmod);
  if (has_substring(generation_method, "EA")) {
    totalmod += 3;
  }
  if (has_substring(generation_method, "EA")) {
    totalmod += 2;
  }  // two dissipative tidal dissipation numbers
  debugger_print(__FILE__, __LINE__, totalmod);
  const int nonModDim = dimension - totalmod;
  // EFPE branches checked before other models to avoid false matches from
  // has_substring()
  if (has_substring(generation_method, "EFPE")) {
    if (nonModDim == 15) {
      std::cout
          << "Sampling in parameters: RA, sin DEC, psi, cos iota, phi_ref, tc, "
             "ln DL, ln chirpmass, eta, a1, a2, cos tilt1, cos tilt2, phi1, "
             "phi2\n";
      *intrinsic = false;
      return;
    } else if (nonModDim == 17) {
      std::cout
          << "Sampling in parameters: RA, sin DEC, psi, cos iota, phi_ref, tc, "
             "ln DL, ln chirpmass, eta, a1, a2, cos tilt1, cos tilt2, phi1, "
             "phi2, e, mean anomaly\n";
      *intrinsic = false;
      return;
    } else if (nonModDim == 8) {
      std::cout << "Sampling in parameters: ln chirpmass, eta, a1, a2, cos "
                   "tilt1, cos tilt2, phi1, phi2\n";
      *intrinsic = true;
      return;
    } else if (nonModDim == 10) {
      std::cout << "Sampling in parameters: ln chirpmass, eta, a1, a2, cos "
                   "tilt1, cos tilt2, phi1, phi2, e, mean anomaly\n";
      *intrinsic = true;
      return;
    }
  }
  if (has_substring(generation_method, "PhenomD") && nonModDim == 4) {
    std::cout << "Sampling in parameters: ln chirpmass, eta, chi1, chi2";
    for (int i = 0; i < totalmod; i++) {
      std::cout << ", mod_" << i;
    }
    std::cout << std::endl;
    *intrinsic = true;
  } else if (has_substring(generation_method, "PhenomD_NRT") &&
             nonModDim == 6) {
    std::cout << "Sampling in parameters: ln chirpmass, eta, chi1, chi2, ln "
                 "tidal1,  ln tidal2";
    for (int i = 0; i < totalmod; i++) {
      std::cout << ", mod_" << i;
    }
    std::cout << std::endl;
    *intrinsic = true;
  } else if (has_substring(generation_method, "PhenomD_NRT") &&
             nonModDim == 8) {
    std::cout << "Sampling in parameters: ln chirpmass, eta, chi1, chi2, ln "
                 "tidal1,  ln tidal2, diss_tidal1, diss_tidal2";
    for (int i = 0; i < totalmod; i++) {
      std::cout << ", mod_" << i;
    }
    std::cout << std::endl;
    *intrinsic = true;
  } else if ((has_substring(generation_method, "PhenomPv2") ||
              has_substring(generation_method, "PhenomPv3")) &&
             nonModDim == 8) {
    std::cout << "Sampling in parameters: ln chirpmass, eta, a1, a2, tilt1, "
                 "tilt2, phi1, phi2";
    for (int i = 0; i < totalmod; i++) {
      std::cout << ", mod_" << i;
    }
    std::cout << std::endl;
    *intrinsic = true;
  } else if (has_substring(generation_method, "PhenomD") && nonModDim == 11) {
    std::cout << "Sampling in parameters: RA, sin DEC, psi, cos iota,phi_ref, "
                 "tc,  ln DL, ln chirpmass, eta, chi1, chi2"
              << std::endl;
    for (int i = 0; i < totalmod; i++) {
      std::cout << ", mod_" << i;
    }
    std::cout << std::endl;
    *intrinsic = false;
  } else if (has_substring(generation_method, "PhenomD_NRT")) {
    if (nonModDim == 13) {
      std::cout
          << "Sampling in parameters: RA, sin  DEC, psi, cos iota,phi_ref, tc, "
             " ln DL, ln chirpmass, eta, chi1, chi2, ln tidal1, ln tidal2"
          << std::endl;
      for (int i = 0; i < totalmod; i++) {
        std::cout << ", mod_" << i;
      }
      std::cout << std::endl;
      *intrinsic = false;
    } else if (nonModDim == 12) {
      std::cout
          << "Sampling in parameters: RA, sin  DEC, psi, cos iota,phi_ref, tc, "
             " ln DL, ln chirpmass, eta, chi1, chi2, ln tidal_s"
          << std::endl;
      for (int i = 0; i < totalmod; i++) {
        std::cout << ", mod_" << i;
      }
      std::cout << std::endl;
      mcmc_intrinsic = false;
    } else if (has_substring(generation_method, "EOS") && nonModDim == 14) {
      std::cout << "Sampling in parameters: RA, sin  DEC, psi, cos iota, "
                   "phi_ref, tc, ln DL, nbc1, nbc2, chi1, chi2, bump_mag, "
                   "bump_width, bump_offset"
                << std::endl;
      for (int i = 0; i < totalmod; i++) {
        std::cout << ", mod_" << i;
      }
      std::cout << std::endl;
      *intrinsic = false;
    } else if (has_substring(generation_method, "EOS") && nonModDim == 15) {
      std::cout << "Sampling in parameters: RA, sin  DEC, psi, cos iota, "
                   "phi_ref, tc, ln DL, nbc1, nbc2, chi1, chi2, bump_mag, "
                   "bump_width, bump_offset, plat"
                << std::endl;
      for (int i = 0; i < totalmod; i++) {
        std::cout << ", mod_" << i;
      }
      std::cout << std::endl;
      *intrinsic = false;
    }
  } else if (has_substring(generation_method, "PhenomPv2") ||
             has_substring(generation_method, "PhenomPv3") && nonModDim == 15) {
    std::cout
        << "Sampling in parameters: RA, sin DEC, psi, cos iota,phi_ref, tc, ln "
           "DL, ln chirpmass, eta, a1, a2,cos tilt1, cos tilt2, phi1, phi2"
        << std::endl;
    for (int i = 0; i < totalmod; i++) {
      std::cout << ", mod_" << i;
    }
    std::cout << std::endl;
    *intrinsic = false;
  } else {
    std::cout
        << "Input parameters not valid, please check that input is compatible "
           "with the supported methods - dimension combinations"
        << std::endl;
    exit(1);
  }
}

/*! \brief utility to do MCMC specific transformations on the input param vector
 * before passing to the repacking utillity
 *
 * Returns the local generation method to be used in the LL functions
 */
std::string MCMC_prep_params(double* param, double* temp_params,
                             gen_params_base<double>* gen_params, int dimension,
                             std::string generation_method,
                             MCMC_modification_struct* mod_struct,
                             bool intrinsic, double gmst) {
  if (intrinsic)
    gen_params->sky_average = true;
  else
    gen_params->sky_average = false;
  gen_params->f_ref = mod_struct->f_ref;
  gen_params->shift_time = true;
  gen_params->shift_phase = true;
  // gen_params->shift_time = false;
  // gen_params->shift_phase = false;
  gen_params->gmst = gmst;
  gen_params->equatorial_orientation = false;
  gen_params->horizon_coord = false;

  gen_params->tidal_love = mod_struct->tidal_love;
  gen_params->tidal_love_error = mod_struct->tidal_love_error;
  gen_params->alpha_param = mod_struct->alpha_param;
  gen_params->EA_region1 = mod_struct->EA_region1;

  gen_params->NSflag1 = mod_struct->NSflag1;
  gen_params->NSflag2 = mod_struct->NSflag2;
  // gen_params->NSflag1 = false;
  // gen_params->NSflag2 = false;
  for (int i = 0; i < dimension; i++) {
    temp_params[i] = param[i];
  }
  int base = dimension;
  if (check_mod(generation_method)) {
    // if(has_substring(generation_method, "ppE") ||
    //	has_substring(generation_method, "dCS")||
    //	has_substring(generation_method, "EdGB")){
    //	gen_params->bppe=mcmc_mod_struct->bppe;
    //	gen_params->Nmod=mcmc_mod_struct->ppE_Nmod;
    //	gen_params->betappe=new double[gen_params->Nmod];
    // }
    if (has_substring(generation_method, "ppE") ||
        check_theory_support(generation_method)) {
      gen_params->bppe = mod_struct->bppe;
      gen_params->Nmod = mod_struct->ppE_Nmod;
      gen_params->betappe = new double[gen_params->Nmod];
      base = dimension - mod_struct->ppE_Nmod;
    } else if (has_substring(generation_method, "gIMR")) {
      gen_params->Nmod_phi = mod_struct->gIMR_Nmod_phi;
      gen_params->phii = mod_struct->gIMR_phii;
      if (gen_params->Nmod_phi != 0) {
        gen_params->delta_phi = new double[gen_params->Nmod_phi];
      }
      gen_params->Nmod_sigma = mod_struct->gIMR_Nmod_sigma;
      gen_params->sigmai = mod_struct->gIMR_sigmai;
      if (gen_params->Nmod_sigma != 0) {
        gen_params->delta_sigma = new double[gen_params->Nmod_sigma];
      }
      gen_params->Nmod_beta = mod_struct->gIMR_Nmod_beta;
      gen_params->betai = mod_struct->gIMR_betai;
      if (gen_params->Nmod_beta != 0) {
        gen_params->delta_beta = new double[gen_params->Nmod_beta];
      }
      gen_params->Nmod_alpha = mod_struct->gIMR_Nmod_alpha;
      gen_params->alphai = mod_struct->gIMR_alphai;
      if (gen_params->Nmod_alpha != 0) {
        gen_params->delta_alpha = new double[gen_params->Nmod_alpha];
      }
      base = dimension - mod_struct->gIMR_Nmod_phi -
             mod_struct->gIMR_Nmod_sigma - mod_struct->gIMR_Nmod_beta -
             mod_struct->gIMR_Nmod_alpha;
    }
    if ((has_substring(generation_method, "dCS") ||
         has_substring(generation_method, "EdGB"))) {
      // temp_params[base] = pow(temp_params[base],.25)/(c*1000);
      temp_params[base] = pow_int(temp_params[base] / (c / 1000.), 4);
    }
  }
  return generation_method;
}

void MCMC_fisher_transformations(double* param, double** fisher, int dimension,
                                 std::string generation_method, bool intrinsic,
                                 MCMC_modification_struct* mod_struct) {
  if (!intrinsic) {
    fisher[0][0] += 1. / (4 * M_PI * M_PI);  // RA
    fisher[1][1] += 1. / 4;                  // sin DEC
    fisher[2][2] += 1. / (4 * M_PI * M_PI);  // psi
    fisher[3][3] += 1. / (4);                // cos iota
    fisher[4][4] += 1. / (4 * M_PI * M_PI);  // phiref
    fisher[5][5] += 1. / (.01);              // tc
    fisher[8][8] += 1. / .25;                // eta
    fisher[9][9] += 1. / 4;                  // spin1
    fisher[10][10] += 1. / 4;                // spin2
    if (has_substring(generation_method, "PhenomPv2") ||
        has_substring(generation_method, "PhenomPv3") ||
        has_substring(generation_method, "EFPE")) {
      fisher[11][11] += 1. / 4;                  // cos theta1
      fisher[12][12] += 1. / 4;                  // cos theta2
      fisher[13][13] += 1. / (4 * M_PI * M_PI);  // phi1
      fisher[14][14] += 1. / (4 * M_PI * M_PI);  // phi2
    }
  } else {
    if (has_substring(generation_method, "PhenomPv2") ||
        has_substring(generation_method, "PhenomPv3") ||
        has_substring(generation_method, "EFPE")) {
      fisher[1][1] = 1. / (.25);              // eta
      fisher[2][2] = 1. / (4);                // spin1
      fisher[3][3] = 1. / (4);                // spin2
      fisher[4][4] = 1. / (4);                // cos theta1
      fisher[5][5] = 1. / (4);                // cos theta2
      fisher[6][6] = 1. / (4 * M_PI * M_PI);  // phi1
      fisher[7][7] = 1. / (4 * M_PI * M_PI);  // phi2
    } else if (has_substring(generation_method, "PhenomD")) {
      fisher[1][1] = 1. / (.25);  // eta
      fisher[2][2] = 1. / (4);    // spin1
      fisher[3][3] = 1. / (4);    // spin2
    }
  }

  if (has_substring(generation_method, "dCS") ||
      has_substring(generation_method, "EdGB")) {
    int base = dimension - mod_struct->ppE_Nmod;
    for (int i = 0; i < dimension; i++) {
      // Transform to root alpha from alpha^2
      double factor = 4 * pow(param[base], 3. / 4.);
      // Transform to km from sec
      factor *= 1000 / c;
      fisher[base][i] *= factor;
      fisher[i][base] *= factor;
    }
  }

  if (has_substring(generation_method, "EA")) {
    fisher[dimension - 3][dimension - 3] += 1. / pow(10, -2.);
    fisher[dimension - 2][dimension - 2] += 1. / pow(10, -4.);
    fisher[dimension - 1][dimension - 1] += 1. / pow(10, -1.);
    // std::cout<<fisher[7][12]<<std::endl;
    // fisher[7][12] = -1.15341564e+05;
    // fisher[12][7] = -1.15341564e+05;
  }

  // if(has_substring(generation_method, "EA")){
  //   //for(int i = 0 ; i <4; i++){
  //     for(int i = 0 ; i <3; i++){
  //       for(int j = 0 ; j<dimension; j++){
  //	if(i!=j){
  //	  fisher[dimension-1-i][dimension-1-j] = 0;
  //	  fisher[dimension-1-j][dimension-1-i] = 0;
  //	}
  //	/*
  //	  if(i==j){
  //	  fisher[dimension-1-i][dimension-1-j] = 1./pow(10, -6.);
  //	  }*/
  //       }
  //     }
  //     fisher[dimension-3][dimension-3] = 1./pow(10, -2.);
  //     fisher[dimension-2][dimension-2] = 1./pow(10, -4.);
  //     fisher[dimension-1][dimension-1] = 1./pow(10, -2.);
  // }
  return;
}

void MCMC_fisher_wrapper_RJ(bayesship::positionInfo* pos, double** output,
                            std::vector<int> block, void* userParameters) {
  mcmcVariablesRJ* mcmcVarRJ = (mcmcVariablesRJ*)userParameters;

  // ##########################################################
  // ##########################################################
  mcmcVariables mcmcVar;
  mcmcVar.mcmc_noise = mcmcVarRJ->mcmc_noise;
  // mcmcVar.mcmc_init_pos = initial_pos;
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

  MCMC_fisher_wrapper_explicit_marginalization(pos, output, block,
                                               (void*)&mcmcVar);

  return;
}

int invertFisherBlock(double** fisherIn, double** fisherOut, int dimIn,
                      std::vector<int> ids) {
  int dimOut = ids.size();
  double** covIn = allocate_2D_array(dimIn, dimIn);
  double** covOut = allocate_2D_array(dimOut, dimOut);

  int status = gsl_cholesky_matrix_invert(fisherIn, covIn, dimIn);

  if (status == 0) {
    for (int i = 0; i < dimOut; i++) {
      for (int j = 0; j < dimOut; j++) {
        covOut[i][j] = covIn[ids[i]][ids[j]];
      }
    }
    status = gsl_cholesky_matrix_invert(covOut, fisherOut, dimOut);
  }

  deallocate_2D_array(covIn, dimIn, dimIn);
  deallocate_2D_array(covOut, dimOut, dimOut);

  return status;
}

void MCMC_fisher_wrapper_explicit_marginalization(bayesship::positionInfo* pos,
                                                  double** output,
                                                  std::vector<int> ids,
                                                  void* userParameters) {
  mcmcVariables* mcmcVar = (mcmcVariables*)userParameters;

  int dimension = mcmcVar->maxDim;
  double** tempOutput = allocate_2D_array(dimension, dimension);
  double* temp_params = new double[dimension];
  double param[dimension];
  for (int i = 0; i < dimension; i++) {
    param[i] = pos->parameters[i];
  }
  gen_params_base<double> params;
  std::string local_gen = MCMC_prep_params(
      param, temp_params, &params, dimension, mcmcVar->mcmc_generation_method,
      mcmcVar->mcmc_mod_struct, mcmcVar->mcmc_intrinsic, mcmcVar->mcmc_gmst);
  repack_parameters(temp_params, &params,
                    "MCMC_" + mcmcVar->mcmc_generation_method, dimension);
  for (int j = 0; j < dimension; j++) {
    for (int k = 0; k < dimension; k++) {
      tempOutput[j][k] = 0;
    }
  }
  double** local_freq = mcmcVar->mcmc_frequencies;
  double** local_noise = mcmcVar->mcmc_noise;
  double** local_weights = mcmcVar->user_parameters->weights;
  int* local_lengths = mcmcVar->mcmc_data_length;
  std::string local_integration_method = "SIMPSONS";
  bool local_log10F = mcmcVar->user_parameters->fisher_log10F;
  if (mcmcVar->user_parameters->fisher_freq) {
    local_freq = mcmcVar->user_parameters->fisher_freq;
  }
  if (mcmcVar->user_parameters->fisher_PSD) {
    local_noise = mcmcVar->user_parameters->fisher_PSD;
  }
  if (mcmcVar->user_parameters->fisher_length) {
    local_lengths = mcmcVar->user_parameters->fisher_length;
  }
  if (mcmcVar->user_parameters->fisher_weights) {
    local_weights = mcmcVar->user_parameters->fisher_weights;
  }
  if (mcmcVar->user_parameters->fisher_GAUSS_QUAD) {
    local_integration_method = "GAUSSLEG";
  }

  std::string local_gen_method = mcmcVar->mcmc_generation_method;
  int local_dimension = dimension;
  double** temp_out = allocate_2D_array(local_dimension, local_dimension);
  for (int i = 0; i < mcmcVar->mcmc_num_detectors; i++) {
    // Use AD
    if (mcmcVar->user_parameters->fisher_AD) {
      std::unique_lock<std::mutex> lock{*(mcmcVar->user_parameters->mFish)};
      fisher_autodiff(local_freq[i], local_lengths[i],
                      "MCMC_" + local_gen_method, mcmcVar->mcmc_detectors[i],
                      mcmcVar->mcmc_detectors[0], temp_out, local_dimension,
                      (gen_params*)(&params), local_integration_method,
                      local_weights[i], true, local_noise[i]);
    } else {
      fisher_numerical(local_freq[i], local_lengths[i],
                       "MCMC_" + local_gen_method, mcmcVar->mcmc_detectors[i],
                       mcmcVar->mcmc_detectors[0], temp_out, local_dimension,
                       &params, 4, nullptr, nullptr, local_noise[i],
                       mcmcVar->user_parameters->QuadMethod);
    }
    for (int j = 0; j < local_dimension; j++) {
      for (int k = 0; k < local_dimension; k++) {
        tempOutput[j][k] += temp_out[j][k];
      }
    }

    MCMC_fisher_transformations(temp_params, tempOutput, dimension, local_gen,
                                mcmcVar->mcmc_intrinsic,
                                mcmcVar->mcmc_mod_struct);
    deallocate_2D_array(temp_out, local_dimension, local_dimension);
  }

  // Try marginalizing over other parameters, otherwise just use subfisher
  // without marginalizing
  int status = invertFisherBlock(tempOutput, output, dimension, ids);
  if (status == 1) {
    for (int i = 0; i < ids.size(); i++) {
      for (int j = 0; j < ids.size(); j++) {
        output[i][j] = tempOutput[ids[i]][ids[j]];
      }
    }
  }

  // Cleanup
  delete[] temp_params;
  deallocate_2D_array(tempOutput, dimension, dimension);
  if (check_mod(local_gen)) {
    if (has_substring(local_gen, "ppE") || check_theory_support(local_gen)) {
      delete[] params.betappe;
    } else if (has_substring(local_gen, "gIMR")) {
      if (mcmcVar->mcmc_mod_struct->gIMR_Nmod_phi != 0) {
        delete[] params.delta_phi;
      }
      if (mcmcVar->mcmc_mod_struct->gIMR_Nmod_sigma != 0) {
        delete[] params.delta_sigma;
      }
      if (mcmcVar->mcmc_mod_struct->gIMR_Nmod_beta != 0) {
        delete[] params.delta_beta;
      }
      if (mcmcVar->mcmc_mod_struct->gIMR_Nmod_alpha != 0) {
        delete[] params.delta_alpha;
      }
    }
  }
}

void MCMC_fisher_wrapper(bayesship::positionInfo* pos, double** output,
                         void* userParameters) {
  mcmcVariables* mcmcVar = (mcmcVariables*)userParameters;

  int dimension = mcmcVar->maxDim;
  double* temp_params = new double[dimension];
  double param[dimension];
  for (int i = 0; i < dimension; i++) {
    param[i] = pos->parameters[i];
  }
  gen_params_base<double> params;
  std::string local_gen = MCMC_prep_params(
      param, temp_params, &params, dimension, mcmcVar->mcmc_generation_method,
      mcmcVar->mcmc_mod_struct, mcmcVar->mcmc_intrinsic, mcmcVar->mcmc_gmst);
  repack_parameters(temp_params, &params,
                    "MCMC_" + mcmcVar->mcmc_generation_method, dimension);
  for (int j = 0; j < dimension; j++) {
    for (int k = 0; k < dimension; k++) {
      output[j][k] = 0;
    }
  }
  double** local_freq = mcmcVar->mcmc_frequencies;
  double** local_noise = mcmcVar->mcmc_noise;
  double** local_weights = mcmcVar->user_parameters->weights;
  int* local_lengths = mcmcVar->mcmc_data_length;
  std::string local_integration_method = "SIMPSONS";
  bool local_log10F = mcmcVar->user_parameters->fisher_log10F;
  if (mcmcVar->user_parameters->fisher_freq) {
    local_freq = mcmcVar->user_parameters->fisher_freq;
  }
  if (mcmcVar->user_parameters->fisher_PSD) {
    local_noise = mcmcVar->user_parameters->fisher_PSD;
  }
  if (mcmcVar->user_parameters->fisher_length) {
    local_lengths = mcmcVar->user_parameters->fisher_length;
  }
  if (mcmcVar->user_parameters->fisher_weights) {
    local_weights = mcmcVar->user_parameters->fisher_weights;
  }
  if (mcmcVar->user_parameters->fisher_GAUSS_QUAD) {
    local_integration_method = "GAUSSLEG";
  }
  double** temp_out = allocate_2D_array(dimension, dimension);
  for (int i = 0; i < mcmcVar->mcmc_num_detectors; i++) {
    // Use AD
    if (mcmcVar->user_parameters->fisher_AD) {
      std::unique_lock<std::mutex> lock{*(mcmcVar->user_parameters->mFish)};
      fisher_autodiff(local_freq[i], local_lengths[i],
                      "MCMC_" + mcmcVar->mcmc_generation_method,
                      mcmcVar->mcmc_detectors[i], mcmcVar->mcmc_detectors[0],
                      temp_out, dimension, (gen_params*)(&params),
                      local_integration_method, local_weights[i], true,
                      local_noise[i]);
    } else {
      fisher_numerical(local_freq[i], local_lengths[i],
                       "MCMC_" + mcmcVar->mcmc_generation_method,
                       mcmcVar->mcmc_detectors[i], mcmcVar->mcmc_detectors[0],
                       temp_out, dimension, &params, 4, nullptr, nullptr,
                       local_noise[i]);
    }
    for (int j = 0; j < dimension; j++) {
      for (int k = 0; k < dimension; k++) {
        output[j][k] += temp_out[j][k];
      }
    }

    MCMC_fisher_transformations(temp_params, output, dimension, local_gen,
                                mcmcVar->mcmc_intrinsic,
                                mcmcVar->mcmc_mod_struct);
    deallocate_2D_array(temp_out, dimension, dimension);
  }

  // Cleanup
  delete[] temp_params;
  if (check_mod(local_gen)) {
    if (has_substring(local_gen, "ppE") || check_theory_support(local_gen)) {
      delete[] params.betappe;
    } else if (has_substring(local_gen, "gIMR")) {
      if (mcmcVar->mcmc_mod_struct->gIMR_Nmod_phi != 0) {
        delete[] params.delta_phi;
      }
      if (mcmcVar->mcmc_mod_struct->gIMR_Nmod_sigma != 0) {
        delete[] params.delta_sigma;
      }
      if (mcmcVar->mcmc_mod_struct->gIMR_Nmod_beta != 0) {
        delete[] params.delta_beta;
      }
      if (mcmcVar->mcmc_mod_struct->gIMR_Nmod_alpha != 0) {
        delete[] params.delta_alpha;
      }
    }
  }
}

/// @brief Wrapper to compute Fisher matrix for BayesShip sampler
void MCMC_fisher_GWParameterSpace_wrapper(bayesship::positionInfo* pos,
                                          double** fisherM,
                                          void* userParameters) {
  const mcmcVariables* mcmcVar = (mcmcVariables*)userParameters;
  const GWParameterSpace* param_space = mcmcVar->param_space;

  gen_params_base<double> params;
  param_space->to_gen_params(pos->parameters, params);
  auto deriv = param_space->fisher_derivatives(mcmcVar->mcmc_deriv_order);
  fisher_numerical(pos->parameters, mcmcVar->likelihood->get_ifos(),
                   deriv.get(), fisherM);
  param_space->fisher_prior(fisherM);
}

/// @brief Wrapper to compute Fisher matrix for BayesShip sampler with
/// marginalization of extrinsic parameters
void MCMC_fisher_GWParameterSpace_wrapper_explicit_marginalization(
    bayesship::positionInfo* pos, double** fisherM, std::vector<int> ids,
    void* userParameters) {
  const mcmcVariables* mcmcVar = (mcmcVariables*)userParameters;
  const GWParameterSpace* param_space = mcmcVar->param_space;

  // Full Fisher matrix calculation
  const int dim = mcmcVar->maxDim;
  double** tempMatrix = allocate_2D_array(dim, dim);
  gen_params_base<double> params;
  param_space->to_gen_params(pos->parameters, params);
  auto deriv = param_space->fisher_derivatives(mcmcVar->mcmc_deriv_order);
  fisher_numerical(pos->parameters, mcmcVar->likelihood->get_ifos(),
                   deriv.get(), tempMatrix);
  param_space->fisher_prior(tempMatrix);

  // Try marginalizing over extrinsic parameters
  int status = invertFisherBlock(tempMatrix, fisherM, dim, ids);
  // If unsuccesful, pass submatrix
  if (status == 1) {
    for (int i = 0; i < ids.size(); i++)
      for (int j = 0; j < ids.size(); j++)
        fisherM[i][j] = tempMatrix[ids[i]][ids[j]];
  }

  deallocate_2D_array(tempMatrix, dim, dim);
}

/*! \brief Maximized match over coalescence variables - returns log likelihood
 * NOT NORMALIZED for aligned spins
 *
 * Note: this function is not properly normalized for an absolute comparison.
 * This is made for MCMC sampling, so to minimize time, constant terms like
 * (Data|Data), which would cancel in the Metropolis-Hasting ratio, are left out
 * for efficiency
 */
double maximized_Log_Likelihood_aligned_spin_internal(
    std::complex<double>* data, double* psd, double* frequencies,
    std::complex<double>* detector_response, size_t length,
    fftw_outline* plan) {
  // Calculate template snr and scale it to match the data snr
  // later, upgrade to non uniform spacing, cause why not
  double delta_f = frequencies[1] - frequencies[0];
  double sum = 0.;
  double* integrand = (double*)malloc(sizeof(double) * length);
  for (int i = 0; i < length; i++)
    integrand[i] =
        real(detector_response[i] * std::conj(detector_response[i])) / psd[i];
  // double integral = 4.*trapezoidal_sum_uniform(delta_f, length, integrand);
  double integral = 4. * simpsons_sum(delta_f, length, integrand);
  double HH = integral;

  // calculate the fourier transform that corresponds to maximizing over phic
  // and tc Use malloc at some point, not sure how long these arrays will be
  std::complex<double> g_tilde;

  fftw_complex* in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * length);
  fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * length);
  for (int i = 0; i < length; i++) {
    g_tilde = 4. * conj(data[i]) * detector_response[i] / psd[i];
    in[i][0] = real(g_tilde);
    in[i][1] = imag(g_tilde);
  }

  double* g = (double*)malloc(sizeof(double) * length);

  fftw_execute_dft(plan->p, in, out);

  for (int i = 0; i < length; i++) {
    g[i] = out[i][0] * out[i][0] + out[i][1] * out[i][1];
  }

  double max = *std::max_element(g, g + length) * delta_f * delta_f;

  free(integrand);
  free(g);
  fftw_free(in);
  fftw_free(out);
  return .5 * (max) / HH;
}

/*! \brief log likelihood function that maximizes over extrinsic parameters tc,
 * phic, D, and phiRef, the reference frequency - for unaligned spins
 *
 * Ref: arXiv 1603.02444v2
 *
 * NOTE: Only works for +/x polarizations
 */
double maximized_Log_Likelihood_unaligned_spin_internal(
    std::complex<double>* data, double* psd, double* frequencies,
    std::complex<double>* hplus, std::complex<double>* hcross, size_t length,
    fftw_outline* plan) {
  double delta_f = frequencies[1] - frequencies[0];
  double* integrand = (double*)malloc(sizeof(double) * length);
  double integral;

  // Calculate template snr for plus polarization sqrt(<H+|H+>)
  for (int i = 0; i < length; i++) {
    integrand[i] = real(hplus[i] * std::conj(hplus[i])) / psd[i];
  }
  integral = 4. * simpsons_sum(delta_f, length, integrand);
  double HpHproot = sqrt(integral);

  // Calculate template snr for cross polarization sqrt(<Hx|Hx>)
  for (int i = 0; i < length; i++)
    integrand[i] = real(hcross[i] * std::conj(hcross[i])) / psd[i];
  integral = 4. * simpsons_sum(delta_f, length, integrand);
  double HcHcroot = sqrt(integral);

  // Rescale waveforms from hplus/cross to \hat{hplus/cross}
  std::complex<double>* hpnorm =
      (std::complex<double>*)malloc(sizeof(std::complex<double>) * length);
  std::complex<double>* hcnorm =
      (std::complex<double>*)malloc(sizeof(std::complex<double>) * length);
  for (int i = 0; i < length; i++) {
    hpnorm[i] = hplus[i] / HpHproot;
    hcnorm[i] = hcross[i] / HcHcroot;
  }

  // calculate \hat{rhoplus/cross} (just denoted rhoplus/cross) <d|hpnorm>
  // To maximize of coalescence phase, this is an FFT (so its a vector of
  //<d|h> at discrete tc
  double* rhoplus2 = (double*)malloc(sizeof(double) * length);
  double* rhocross2 = (double*)malloc(sizeof(double) * length);
  std::complex<double>* rhoplus =
      (std::complex<double>*)malloc(sizeof(std::complex<double>) * length);
  std::complex<double>* rhocross =
      (std::complex<double>*)malloc(sizeof(std::complex<double>) * length);
  double* gammahat = (double*)malloc(sizeof(double) * length);
  std::complex<double> g_tilde;
  fftw_complex* in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * length);
  fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * length);

  for (int i = 0; i < length; i++) {
    g_tilde = 4. * conj(data[i]) * hpnorm[i] / psd[i];
    in[i][0] = real(g_tilde);
    in[i][1] = imag(g_tilde);
  }
  // fftw_execute(plan->p);
  fftw_execute_dft(plan->p, in, out);
  for (int i = 0; i < length; i++) {
    rhoplus[i] = delta_f * (std::complex<double>(out[i][0], out[i][1]));
    // Norm of the output, squared (Re{g}^2 + Im{g}^2)
    rhoplus2[i] =
        delta_f * delta_f * (out[i][0] * out[i][0] + out[i][1] * out[i][1]);
  }
  for (int i = 0; i < length; i++) {
    g_tilde = 4. * conj(data[i]) * hcnorm[i] / psd[i];
    in[i][0] = real(g_tilde);
    in[i][1] = imag(g_tilde);
  }
  // fftw_execute(plan->p);
  fftw_execute_dft(plan->p, in, out);
  for (int i = 0; i < length; i++) {
    rhocross[i] = delta_f * (std::complex<double>(out[i][0], out[i][1]));
    // Norm of the output, squared (Re{g}^2 + Im{g}^2)
    rhocross2[i] =
        delta_f * delta_f * (out[i][0] * out[i][0] + out[i][1] * out[i][1]);
  }

  for (int i = 0; i < length; i++)
    gammahat[i] = real(rhoplus[i] * conj(rhocross[i]));

  for (int i = 0; i < length; i++)
    integrand[i] = real(hpnorm[i] * std::conj(hcnorm[i])) / psd[i];
  integral = 4. * simpsons_sum(delta_f, length, integrand);
  double Ipc = integral;

  double* lambda = (double*)malloc(sizeof(double) * length);
  for (int i = 0; i < length; i++) {
    lambda[i] =
        (rhoplus2[i] + rhocross2[i] - 2 * gammahat[i] * Ipc +
         sqrt((rhoplus2[i] - rhocross2[i]) * (rhoplus2[i] - rhocross2[i]) +
              4. * (Ipc * rhoplus2[i] - gammahat[i]) *
                  (Ipc * rhocross2[i] - gammahat[i]))) /
        (1. - Ipc * Ipc);
  }
  double max = .25 * (*std::max_element(lambda, lambda + length));

  free(integrand);
  free(hpnorm);
  free(hcnorm);
  free(rhoplus2);
  free(rhoplus);
  free(rhocross);
  free(rhocross2);
  free(gammahat);
  free(lambda);
  fftw_free(in);
  fftw_free(out);

  return max;
}

/*! \brief Internal function for the unmarginalized log of the likelihood
 *
 * .5 * ( ( h | h ) - 2 ( D | h ) )
 */
double Log_Likelihood_internal(std::complex<double>* data, double* psd,
                               double* frequencies, double* weights,
                               std::complex<double>* detector_response,
                               int length, bool log10F,
                               std::string integration_method) {
  double delta_f = frequencies[length / 2] - frequencies[length / 2 - 1];
  double sum = 0.;
  double* integrand = (double*)malloc(sizeof(double) * length);
  for (int i = 0; i < length; i++) {
    integrand[i] =
        real(detector_response[i] * std::conj(detector_response[i])) / psd[i];
  }
  double integral = 0;
  if (integration_method == "SIMPSONS") {
    integral = 4. * simpsons_sum(delta_f, length, integrand);
  } else if (integration_method == "GAUSSLEG") {
    if (log10F) {
      for (int i = 0; i < length; i++) {
        integral += weights[i] * integrand[i] * frequencies[i] * LOG10;
      }
    } else {
      for (int i = 0; i < length; i++) {
        integral += weights[i] * integrand[i];
      }
    }
    integral *= 4;
  }
  double HH = integral;
  integral = 0;

  for (int i = 0; i < length; i++) {
    integrand[i] = real(data[i] * std::conj(detector_response[i])) / psd[i];
  }
  if (integration_method == "SIMPSONS") {
    integral = 4. * simpsons_sum(delta_f, length, integrand);
  } else if (integration_method == "GAUSSLEG") {
    if (log10F) {
      for (int i = 0; i < length; i++) {
        integral += weights[i] * integrand[i] * frequencies[i] * LOG10;
      }
    } else {
      for (int i = 0; i < length; i++) {
        integral += weights[i] * integrand[i];
      }
    }
    integral *= 4;
  }
  double DH = integral;

  free(integrand);
  return -0.5 * (HH - 2 * DH);
}

struct skysearch_params {
  std::complex<double>* hplus;
  std::complex<double>* hcross;
};

double MCMC_likelihood_extrinsic(bool save_waveform,
                                 gen_params_base<double>* parameters,
                                 std::string generation_method,
                                 int* data_length, double** frequencies,
                                 std::complex<double>** data, double** psd,
                                 double** weights,
                                 std::string integration_method, bool log10F,
                                 std::string* detectors, int num_detectors) {
  double ll = 0;
  std::complex<double>** responses = new std::complex<double>*[num_detectors];
  for (int i = 0; i < num_detectors; i++) {
    responses[i] = new std::complex<double>[data_length[i]];
  }

  if (num_detectors == 1) {
    create_single_GW_detection(responses[0], detectors[0], frequencies[0],
                               data_length[0], parameters, generation_method);
  } else {
    create_coherent_GW_detection(detectors, num_detectors, frequencies,
                                 data_length, save_waveform, parameters,
                                 generation_method, responses);
  }

  for (int i = 0; i < num_detectors; i++) {
    ll += Log_Likelihood_internal(data[i], psd[i], frequencies[i], weights[i],
                                  responses[i], data_length[i], log10F,
                                  integration_method);
  }

  for (int i = 0; i < num_detectors; i++) {
    delete[] responses[i];
  }
  delete[] responses;

  return ll;
}

// ============================================================
//  find_fiducial
// ============================================================

void find_fiducial(const GWParameterSpace& param_space,
                   const gw_likelihoods::Likelihood& likelihood,
                   const double* initial_params,
                   bayesship::probabilityFn* log_prior,
                   const std::vector<double>& param_scales, int num_mh_steps,
                   std::vector<VECCPL>& fiducial_out,
                   std::vector<VECCPL>& test_out,
                   gen_params_base<double>* test_gp_out) {
  const int dimension = param_space.dim();
  const std::string gen_method = param_space.generation_method();
  const std::vector<IfoData>& ifos = likelihood.get_ifos();
  const int num_detectors = static_cast<int>(ifos.size());

  auto eval_ll = [&](const double* params) -> double {
    try {
      gen_params_base<double> gp;
      param_space.to_gen_params(params, gp);
      double ll = likelihood.log_likelihood(&gp, gen_method);
      return std::isnan(ll) ? -std::numeric_limits<double>::infinity() : ll;
    } catch (const std::exception&) {
      return -std::numeric_limits<double>::infinity();
    }
  };

  auto gen_resp = [&](const double* params) -> std::vector<VECCPL> {
    try {
      gen_params_base<double> gp;
      param_space.to_gen_params(params, gp);
      return create_coherent_GW_detection_reuse_wf(ifos, &gp, gen_method);
    } catch (const std::exception&) {
      std::vector<VECCPL> zeros(num_detectors);
      for (int d = 0; d < num_detectors; d++)
        zeros[d].assign(ifos[d].freqs.size(), CPL(0., 0.));
      return zeros;
    }
  };

  // Chain state
  std::vector<double> current(initial_params, initial_params + dimension);
  std::vector<double> proposed(dimension);
  std::vector<double> best = current;

  double current_ll = eval_ll(current.data());
  // Single positionInfo reused for all prior evaluations (avoids double-free
  // from swapping objects that have raw pointer members).
  bayesship::positionInfo pos_tmp(dimension);
  for (int i = 0; i < dimension; i++) pos_tmp.parameters[i] = current[i];
  double current_lp = log_prior->eval(&pos_tmp, 0);
  double best_ll = current_ll + current_lp;

  // Fisher diagonal for M-H proposal widths via the O2 derivative class.
  std::vector<double> sigma(dimension);
  const double c_mh = 2.38 / std::sqrt((double)dimension);
  {
    std::vector<double> fisher_storage(dimension * dimension, 0.0);
    std::vector<double*> fisher_mat(dimension);
    for (int i = 0; i < dimension; i++)
      fisher_mat[i] = fisher_storage.data() + i * dimension;
    auto deriv = param_space.fisher_derivatives(2);
    fisher_numerical(current.data(), ifos, deriv.get(), fisher_mat.data());
    std::cout << "Fiducial Fisher diagonal (MCMC params):\n";
    for (int i = 0; i < dimension; i++) {
      double gamma_ii = fisher_mat[i][i];
      double fisher_sigma =
          (gamma_ii > 0.0) ? c_mh / std::sqrt(gamma_ii) : 0.1 * param_scales[i];
      double prior_half = 0.5 * param_scales[i];
      sigma[i] = std::min(fisher_sigma, prior_half);
      std::cout << "  param " << i << ": Gamma_ii=" << gamma_ii
                << "  sigma=" << sigma[i] << "  ["
                << (sigma[i] < fisher_sigma ? "prior" : "Fisher") << "]\n";
    }
  }

  std::cout << "Fiducial M-H start logL = " << current_ll << "\n";

  // Per-parameter (Gibbs-style) M-H
  int proposals_in_prior = 0, proposals_accepted = 0;
  std::vector<int> param_in_prior(dimension, 0), param_accepted(dimension, 0);
  gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);
  auto mh_t0 = std::chrono::steady_clock::now();

  for (int step = 0; step < num_mh_steps; step++) {
    int i = step % dimension;
    double old_val = current[i];
    double new_val = old_val + gsl_ran_gaussian(rng, sigma[i]);
    pos_tmp.parameters[i] = new_val;

    double proposed_lp = log_prior->eval(&pos_tmp, 0);
    if (!std::isfinite(proposed_lp)) {
      pos_tmp.parameters[i] = old_val;
      continue;
    }
    proposals_in_prior++;
    param_in_prior[i]++;

    proposed = current;
    proposed[i] = new_val;
    double proposed_ll = eval_ll(proposed.data());
    double log_alpha = (proposed_ll + proposed_lp) - (current_ll + current_lp);

    if (log_alpha >= 0. || std::log(gsl_rng_uniform(rng)) < log_alpha) {
      current[i] = new_val;
      current_ll = proposed_ll;
      current_lp = proposed_lp;
      proposals_accepted++;
      param_accepted[i]++;
      if (current_ll + current_lp > best_ll) {
        best = current;
        best_ll = current_ll + current_lp;
      }
    } else {
      pos_tmp.parameters[i] = old_val;
    }
  }

  auto mh_t1 = std::chrono::steady_clock::now();
  double mh_ms =
      std::chrono::duration<double, std::milli>(mh_t1 - mh_t0).count();
  std::cout << "Fiducial M-H MAP logL   = " << best_ll << "\n";
  std::cout << "Fiducial M-H final logL = " << current_ll << "\n";
  std::cout << "M-H loop time: " << mh_ms << " ms\n";
  std::cout << "Proposal acceptance fraction: " << proposals_accepted << "/"
            << proposals_in_prior << " in-prior proposals accepted\n";
  std::cout << "Per-parameter acceptance (in-prior):\n";
  for (int i = 0; i < dimension; i++) {
    int tot = param_in_prior[i], acc = param_accepted[i];
    std::cout << "  param " << i << ": " << acc << "/" << tot;
    if (tot > 0) std::cout << "  (" << (100 * acc / tot) << "%)";
    std::cout << "\n";
  }

  fiducial_out = gen_resp(best.data());
  test_out = gen_resp(current.data());

  if (test_gp_out != nullptr)
    param_space.to_gen_params(current.data(), *test_gp_out);

  gsl_rng_free(rng);
}
