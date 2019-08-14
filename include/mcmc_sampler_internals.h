#ifndef MCMC_SAMPLER_INTERNALS_H
#define MCMC_SAMPLER_INTERNALS_H
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <string>
#include <functional>
#include <limits>
#include <iomanip>

/*! \file
 * Internal functions of the generic MCMC sampler (nothing specific to GW)
 */

//typedef double (*log_prior_thread_safe)(double *param, int dimension);
//typedef double (*log_likelihood_thread_safe)(double *param, int dimension);
//typedef void (*fisher_thread_safe)(double *param, int dimension,double **fisher);

const double limit_inf = -std::numeric_limits<double>::infinity();

/*! Class storing everything that defines an instance of the sampler
 */
class sampler
{
public:
	int types_of_steps = 5;
	double **step_prob;
	double **prob_boundaries;
	double *chain_temps;
	bool *waiting;
	int *chain_pos;
	double swp_freq;
	int chain_N;
	int numThreads;
	int N_steps;
	int dimension;
	int min_dim;
	int max_dim;
	bool fisher_exist;
	bool *de_primed;
	int *priority;
	bool *ref_chain_status;

	double ***output;

	bool pool;
	int progress=0;
	bool show_progress;
	int num_threads;

	int history_length;
	int history_update;
	int *current_hist_pos;
	double ***history;
	double *current_likelihoods;


	int *check_stepsize_freq;
	double *max_target_accept_ratio;
	double *min_target_accept_ratio;
	int *gauss_last_accept_ct;
	int *gauss_last_reject_ct;
	int *de_last_accept_ct;
	int *de_last_reject_ct;
	int *fish_last_accept_ct;
	int *fish_last_reject_ct;
	double **randgauss_width;

	double ***fisher_vecs;
	double **fisher_vals;
	int *fisher_update_ct;
	int fisher_update_number;

	//log_prior lp;
	//log_likelihood ll;
	//fisher fish;
	std::function<double(double*,int, int)> lp;
	std::function<double(double*,int, int)> ll;
	std::function<void(double*,int,double**,int)> fish;
 
	gsl_rng ** rvec;

	int *nan_counter;
	int *num_gauss ;
	int *num_fish ;
	int *num_de ;
	int *num_mmala ;
	int *num_RJstep ;

	double time_elapsed_cpu;
	double time_elapsed_wall;
	double time_elapsed_cpu_ac;
	double time_elapsed_wall_ac;

	int *fish_accept_ct;
	int *fish_reject_ct;
	int *de_accept_ct;
	int *de_reject_ct;
	int *gauss_accept_ct;
	int *gauss_reject_ct;
	int *mmala_accept_ct;
	int *mmala_reject_ct;
	int *RJstep_accept_ct;
	int *RJstep_reject_ct;

	int *swap_accept_ct;
	int *swap_reject_ct;
	int *step_accept_ct;
	int *step_reject_ct;

	//Pointer for testing -- stores log likelihood and log prior for output
	//Quite a bit of memory (order chain_N * N_steps * 2), so probably not good to run every single time. Just for trouble
	//shooting.
	double ***ll_lp_output;
	bool log_ll=false;
	bool log_lp=false;



	//Parameters for dynamic PT allocation
	int *A;
	bool PT_alloc=false;

	//RJPTMCMC Parameterts
	int ***param_status;
	bool RJMCMC=false;
};


int mcmc_step(sampler *sampler, double *current_param,double *next_param, int *current_status, int *proposed_status, int chain_number);

void gaussian_step(sampler *sampler, double *current_param,double *proposed_param, int *current_status, int *proposed_status, int chain_id);

void fisher_step(sampler *sampler,double *current_param, double *proposed_param,int *current_status, int *proposed_status, int chain_index);

void update_fisher(sampler *sampler, double *current_param, int chain_index);

void mmala_step(sampler *sampler,double *current_param, double *proposed_param,int *current_status, int *proposed_status);

void diff_ev_step(sampler *sampler,double *current_param, double *proposed_param,int *current_status, int *proposed_status, int chain_id);

void RJ_step(sampler *sampler, 
	double *current_param, 
	double *proposed_param, 
	int *current_status, 
	int *proposed_status, 
	int chain_number);

void chain_swap(sampler *sampler, double ***output, int step_num,int *swp_accepted, int *swp_rejected);

int single_chain_swap(sampler *sampler, double *chain1, double *chain2,int T1_index, int T2_index);

void assign_probabilities(sampler *sampler, int chain_index);

void transfer_chain(sampler *samplerptr_dest,sampler *samplerptr_source, int id_destination, int id_source, bool transfer_output);

bool check_sampler_status(sampler *samplerptr);

void update_step_widths(sampler *samplerptr, int chain_id);

void allocate_sampler_mem(sampler *sampler);

void deallocate_sampler_mem(sampler *sampler);

void update_history(sampler *sampler, double *new_params, int chain_index);


//void auto_correlation_spectral(double *chain, int length, double *autocorr);
//
//double auto_correlation(double *arr, int length, double tolerance);
//
//double auto_correlation_serial(double *arr, int length);
//
//double auto_correlation_grid_search(double *arr, int length, int box_num=10, int final_length= 50, double target_length=.01);
//
//double auto_correlation_internal(double *arr, int length, int lag, double ave);
//
//void auto_corr_intervals(double *data, int length,double *output, int num_segments,  double accuracy);
//
//void write_auto_corr_file_from_data(std::string autocorr_filename, double **output,int intervals, int dimension, int N_steps);
//
//void write_auto_corr_file_from_data_file(std::string autocorr_filename, std::string output_file,int intervals, int dimension, int N_steps);

//void write_stat_file(sampler *sampler, std::string filename, int *accepted_steps, int *rejected_steps,int accepted_swps, int rejected_swps);
void write_stat_file(sampler *sampler, std::string filename);

void write_checkpoint_file(sampler *sampler, std::string filename);

void load_checkpoint_file(std::string check_file, sampler *sampler);

void load_temps_checkpoint_file(std::string check_file, double *temps, int chain_N);

void assign_ct_p(sampler *sampler, int step, int chain_index);
void assign_ct_m(sampler *sampler, int step, int chain_index);

void assign_initial_pos(sampler *samplerptr,double *initial_pos, double *seeding_var) ;

double PT_dynamical_timescale(int t0, int nu, int t);

void update_temperatures(sampler *samplerptr,
	int t0,
	int nu,
	int t
	);
void initiate_full_sampler(sampler *sampler_new, sampler *sampler_old, 
	int chain_N_thermo_ensemble, 
	int chain_N,
	std::string chain_allocation_scheme
	);
#endif
