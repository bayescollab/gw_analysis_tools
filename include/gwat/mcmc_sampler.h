#ifndef MCMC_SAMPLER_H
#define MCMC_SAMPLER_H
#include <iostream>
#include <functional>
/*! \file 
 * Header file for mcmc_sampler
 */
#include "mcmc_sampler_internals.h"
void mcmc_step_threaded(int j);
void mcmc_swap_threaded(int i, int j);

void continue_RJPTMCMC_MH(std::string start_checkpoint_file,
	double ***output,
	int ***status,
	int N_steps,
	int swp_freq,
	double (*log_prior)(double *param, int *status, int max_dimension, int chain_id),	
	double (*log_likelihood)(double *param, int *status, int max_dimension, int chain_id),
	void (*fisher)(double *param, int *status,int max_dimension, double **fisher, int chain_id),
	void(*RJ_proposal)(double *current_param, double *proposed_param, int *current_status, int *proposed_status, int max_dimension, int chain_id, double width),
	int numThreads,
	bool pool,
	bool show_prog,
	std::string statistics_filename,
	std::string chain_filename,
	std::string auto_corr_filename,
	std::string likelihood_log_filename,
	std::string end_checkpoint_file
	);
void continue_RJPTMCMC_MH(std::string start_checkpoint_file,
	double ***output,
	int ***status,
	int N_steps,
	int swp_freq,
	double (*log_prior)(double *param, int *status, int max_dimension, int chain_id),	
	double (*log_likelihood)(double *param, int *status, int max_dimension, int chain_id),
	void (*fisher)(double *param, int *status,int max_dimension, double **fisher, int chain_id),
	void(*RJ_proposal)(double *current_param, double *proposed_param, int *current_status, int *proposed_status, int max_dimension, int chain_id),
	int numThreads,
	bool pool,
	bool show_prog,
	std::string statistics_filename,
	std::string chain_filename,
	std::string auto_corr_filename,
	std::string likelihood_log_filename,
	std::string end_checkpoint_file
	);
void continue_RJPTMCMC_MH_internal(std::string start_checkpoint_file,
	double ***output,
	int ***status,
	int N_steps,
	int swp_freq,
	std::function<double(double*,int *, int,int)> log_prior,
	std::function<double(double*,int*,int,int)> log_likelihood,
	std::function<void(double*,int*,int,double**,int)>fisher,
	std::function<void(double*,double*, int*,int*,int,int)> RJ_proposal,
	int numThreads,
	bool pool,
	bool show_prog,
	bool update_RJ_width, 
	std::string statistics_filename,
	std::string chain_filename,
	std::string auto_corr_filename,
	std::string likelihood_log_filename,
	std::string end_checkpoint_file
	);
void RJPTMCMC_MH_internal(double ***output, 
	int ***parameter_status, 	
	int max_dimension, 	
	int min_dimension, 	
	int N_steps,	
	int chain_N,	
	double *initial_pos, 	
	int *initial_status, 	
	double *seeding_var, 	
	double *chain_temps,	
	int swp_freq,	
	std::function<double(double*,int *, int,int)> log_prior,
	std::function<double(double*,int*,int,int)> log_likelihood,
	std::function<void(double*,int*,int,double**,int)>fisher,
	std::function<void(double*,double*, int*,int*,int,int)> RJ_proposal,
	int numThreads, 
	bool pool, 
	bool show_prog, 
	bool update_RJ_width, 
	std::string statistics_filename,
	std::string chain_filename,
	std::string auto_corr_filename,
	std::string likelihood_log_filename,
	std::string checkpoint_file
	);
void RJPTMCMC_MH(double ***output, 
	int ***parameter_status, 
	int max_dimension, 
	int min_dimension, 	
	int N_steps,	
	int chain_N,	
	double *initial_pos, 	
	int *initial_status, 	
	double *seeding_var, 	
	double *chain_temps,	
	int swp_freq,	
	double (*log_prior)(double *param, int *status, int max_dimension, int chain_id),	
	double (*log_likelihood)(double *param, int *status, int max_dimension, int chain_id),
	void (*fisher)(double *param, int *status,int max_dimension, double **fisher, int chain_id),
	void(*RJ_proposal)(double *current_param, double *proposed_param, int *current_status, int *proposed_status, int max_dimension, int chain_id),
	int numThreads, 
	bool pool, 
	bool show_prog, 
	std::string statistics_filename,
	std::string chain_filename,
	std::string auto_corr_filename,
	std::string likelihood_log_filename,
	std::string checkpoint_file
	);
void RJPTMCMC_MH(double ***output, 
	int ***parameter_status, 
	int max_dimension, 
	int min_dimension, 	
	int N_steps,	
	int chain_N,	
	double *initial_pos, 	
	int *initial_status, 	
	double *seeding_var, 	
	double *chain_temps,	
	int swp_freq,	
	double (*log_prior)(double *param, int *status, int max_dimension, int chain_id),	
	double (*log_likelihood)(double *param, int *status, int max_dimension, int chain_id),
	void (*fisher)(double *param, int *status,int max_dimension, double **fisher, int chain_id),
	void(*RJ_proposal)(double *current_param, double *proposed_param, int *current_status, int *proposed_status, int max_dimension, int chain_id, double gaussian_width),
	int numThreads, 
	bool pool, 
	bool show_prog, 
	std::string statistics_filename,
	std::string chain_filename,
	std::string auto_corr_filename,
	std::string likelihood_log_filename,
	std::string checkpoint_file
	);

void PTMCMC_MH_dynamic_PT_alloc_internal(double ***output, 
	int dimension, 	
	int N_steps,	
	int chain_N,	
	int max_chain_N_thermo_ensemble,	
	double *initial_pos, 	
	double *seeding_var, 	
	double *chain_temps,
	int swp_freq,	
	int t0,
	int nu,
	std::string chain_distribution_scheme, 
	std::function<double(double*,int *,int,int)> log_prior,
	std::function<double(double*,int *,int,int)> log_likelihood,
	std::function<void(double*,int*,int,double**,int)>fisher,
	int numThreads, 
	bool pool, 
	bool show_prog, 
	std::string statistics_filename,
	std::string chain_filename,
	std::string likelihood_log_filename,
	std::string checkpoint_file
	);
void PTMCMC_MH_dynamic_PT_alloc(double ***output, 
	int dimension, 	
	int N_steps,	
	int chain_N,	
	int max_chain_N_thermo_ensemble,	
	double *initial_pos, 	
	double *seeding_var, 	
	double *chain_temps,
	int swp_freq,	
	int t0,
	int nu,
	std::string chain_distribution_scheme, 
	double (*log_prior)(double *param, int dimension, int chain_id),	
	double (*log_likelihood)(double *param, int dimension, int chain_id),
	void (*fisher)(double *param, int dimension, double **fisher, int chain_id),
	int numThreads, 
	bool pool, 
	bool show_prog, 
	std::string statistics_filename,
	std::string chain_filename,
	std::string likelihood_log_filename,
	std::string checkpoint_file
	);
void PTMCMC_MH_dynamic_PT_alloc(double ***output, 
	int dimension, 	
	int N_steps,	
	int chain_N,	
	int max_chain_N_thermo_ensemble,	
	double *initial_pos, 	
	double *seeding_var, 	
	double *chain_temps,
	int swp_freq,	
	int t0,
	int nu,
	std::string chain_distribution_scheme, 
	double (*log_prior)(double *param, int dimension),	
	double (*log_likelihood)(double *param, int dimension),
	void (*fisher)(double *param, int dimension, double **fisher),
	int numThreads, 
	bool pool, 
	bool show_prog, 
	std::string statistics_filename,
	std::string chain_filename,
	std::string likelihood_log_filename,
	std::string checkpoint_file
	);
void continue_PTMCMC_MH(std::string start_checkpoint_file,
	double ***output,
	int N_steps,
	int swp_freq,
	double (*log_prior)(double *param, int dimension, int chain_id),
	double (*log_likelihood)(double *param, int dimension, int chain_id),
	void (*fisher)(double *param, int dimension, double **fisher, int chain_id),
	int numThreads,
	bool pool,
	bool show_prog,
	std::string statistics_filename,
	std::string chain_filename,
	std::string auto_corr_filename,
	std::string likelihood_log_filename,
	std::string end_checkpoint_file
	);
void continue_PTMCMC_MH(std::string start_checkpoint_file,
	double ***output,
	int N_steps,
	int swp_freq,
	double (*log_prior)(double *param, int dimension),	
	double (*log_likelihood)(double *param, int dimension),	
	void (*fisher)(double *param, int dimension, double **fisher),
	int numThreads,
	bool pool,
	bool show_prog,
	std::string statistics_filename,
	std::string chain_filename,
	std::string auto_corr_filename,
	std::string likelihood_log_filename,
	std::string end_checkpoint_file
	);
void PTMCMC_MH_loop(sampler *sampler);

void PTMCMC_MH_step_incremental(sampler *sampler, int increment);

void PTMCMC_MH(	double ***output, 
	int dimension, 	
	int N_steps,	
	int chain_N,	
	double *initial_pos, 	
	double *seeding_var,
	double *chain_temps,	
	int swp_freq,	
	double (*log_prior)(double *param, int dimension),	
	double (*log_likelihood)(double *param, int dimension),	
	void (*fisher)(double *param, int dimension, double **fisher),
	int numThreads,
	bool pool,
	bool show_prog,
	std::string statistics_filename,
	std::string chain_filename,
	std::string auto_corr_filename,
	std::string likelihood_log_filename,
	std::string checkpoint_filename
	);
void PTMCMC_MH(	double ***output, 
	int dimension, 	
	int N_steps,	
	int chain_N,	
	double *initial_pos, 	
	double *seeding_var,
	double *chain_temps,	
	int swp_freq,	
	double (*log_prior)(double *param, int dimension, int chain_id),	
	double (*log_likelihood)(double *param, int dimension, int chain_id),	
	void (*fisher)(double *param, int dimension, double **fisher, int chain_id),
	int numThreads,
	bool pool,
	bool show_prog,
	std::string statistics_filename,
	std::string chain_filename,
	std::string auto_corr_filename,
	std::string likelihood_log_filename,
	std::string checkpoint_filename
	);
void PTMCMC_MH_internal(	double ***output, 
	int dimension, 	
	int N_steps,	
	int chain_N,	
	double *initial_pos, 	
	double *seeding_var,
	double *chain_temps,	
	int swp_freq,	
	std::function<double(double*,int *,int,int)> log_prior,
	std::function<double(double*,int *,int,int)> log_likelihood,
	std::function<void(double*,int *,int,double**,int)>fisher,
	int numThreads,
	bool pool,
	bool show_prog,
	std::string statistics_filename,
	std::string chain_filename,
	std::string auto_corr_filename,
	std::string likelihood_log_filename,
	std::string checkpoint_filename
	);
void continue_PTMCMC_MH_internal(std::string start_checkpoint_file,
	double ***output,
	int N_steps,
	int swp_freq,
	std::function<double(double*,int *,int,int)> log_prior,
	std::function<double(double*,int *,int,int)> log_likelihood,
	std::function<void(double*,int *,int,double**,int)>fisher,
	int numThreads,
	bool pool,
	bool show_prog,
	std::string statistics_filename,
	std::string chain_filename,
	std::string auto_corr_filename,
	std::string likelihood_log_filename,
	std::string end_checkpoint_file
	);
#endif