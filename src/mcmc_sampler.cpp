#include "mcmc_sampler.h"
#include "autocorrelation.h"
#include "util.h"
#include "mcmc_sampler_internals.h"
#include "mcmc_io_util.h"
#include "threadPool.h"
#include "io_util.h"
#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <condition_variable>
#include <mutex>
#include <thread>
#include <vector>
#include <algorithm>
#include <utility>
#include <queue>
#include <functional>
#include <unistd.h>
#include <fstream>
#include "adolc/adolc.h"
//#include <adolc/adolc_openmp.h>

//double TMAX = 10e2;
//double TMAX = 1e8;
double TMAX = 1e10;
//double TMAX = 1e14;

#ifndef _OPENMP
#define omp ignore
#endif

#ifdef _OPENMP
#include <adolc/adolc_openmp.h>
#include <omp.h>
#endif

/*!\file 
 * Source file for the sampler foundation
 *
 * Source file for generic MCMC sampler. Sub routines that are 
 * application agnostic are housed in mcmc_sampler_internals
 */

//Random number variables
const gsl_rng_type *T;
gsl_rng * r;
sampler *samplerptr_global;


//double *testing_thread_access=NULL;
//double testing_thread_access[50000];
//int thread_access_ct = 0;

//Time for queueing threads to wait after assigning job
//Helps to even out the load
int THREADWAIT=500; //in microseconds -- should only need a few clock cycles
int QUEUE_SLEEP=int(1e2);

//######################################################################################
//######################################################################################
/*! \brief Class to facilitate the comparing of chains for priority
 *
 * 3 levels of priority: 0 (high) 1 (default) 2 (low)
 */
class Comparator
{
public:
	bool operator()(int i, int j)
	{
		return samplerptr_global->priority[i]>samplerptr_global->priority[j];	
	}
};
class Comparatorswap
{
public:
	bool operator()(int i, int j)
	{
		return false;	
	}
};
bool temp_neighborhood_check(int i, int j ){
	//std::cout<<"Comparison: "<<samplerptr->chain_temps[i]-samplerptr->chain_temps[j]<<std::endl;
	//std::cout<<"IDS: "<<i<<" "<<j<<std::endl;
	if(samplerptr_global->restrict_swapping){
		if(samplerptr_global->isolate_ensembles_cold
		&& 
		((fabs(samplerptr_global->chain_temps[i]-1) < 1e-10) ||(fabs(samplerptr_global->chain_temps[j]-1) < 1e-10) ) 
		){
			if(
			check_list(j,samplerptr_global->chain_neighborhoods_ids[i],samplerptr_global->chain_neighbors[i])
			){
				return true;
			}
			else {
				return false;
			}

		}
		else{
			if(
			(samplerptr_global->isolate_ensembles && check_list(j,samplerptr_global->chain_neighborhoods_ids[i],samplerptr_global->chain_neighbors[i])) 
			|| 
			(!samplerptr_global->isolate_ensembles && check_list(samplerptr_global->chain_temps[j],samplerptr_global->chain_neighborhoods[i],samplerptr_global->chain_neighbors[i]))
			){
				return true;
			}
			else{
				return false;
			}
		}
	}
	return true;
	//return true;
}
struct swap_struct
{
	int id1;
	int id2;
};
class ThreadPool
{
public:
	//using Task = std::function<void()>;
	
	
	explicit ThreadPool(std::size_t numThreads)
	{
		if(randomize_step){
			gsl_rng_env_setup();
			T_step = gsl_rng_default;
			r_step = gsl_rng_alloc(T_step);
		}
		start(numThreads);
	}

	~ThreadPool()
	{
		stop();
	}

	int queue_length()
	{
		int length;
		{ 
			std::unique_lock<std::mutex> lock{EventMutexStep};
			length = step_queue.size();
		}
		return length;
	}
	int queue_swap_length()
	{
		int length;
		{ 
			std::unique_lock<std::mutex> lock{EventMutexSWP_paired};
			length = swap_queue_paired.size();
		}
		return length;
	}
	int queue_swap_element(int i )
	{
		int element;
		{ 
			std::unique_lock<std::mutex> lock{EventMutexSWP_paired};
			element = swap_queue_pair.at(i);
		}
		return element;
	}
	int queue_pairs_length()
	{
		int length;
		{ 
			std::unique_lock<std::mutex> lock{EventMutexSWP_paired};
			length = swap_queue_paired.size();
		}
		return length;
	}

	void enqueue(int i)
	{
		if (randomize_step){
			std::unique_lock<std::mutex> lock{EventMutexStep};
			step_queue_vector.push_back(std::move(i));

		}
		else
		{
			std::unique_lock<std::mutex> lock{EventMutexStep};
			step_queue.emplace(std::move(i));
			//if(thread_access_ct<50000){
			//	testing_thread_access[thread_access_ct] = i;
			//	thread_access_ct ++;
			//}
			//else{
			//	write_file("data/thread_access_test.csv",testing_thread_access,50000);
			//	thread_access_ct = 0;
			//}
		}
		EventVarStep.notify_one();
		//usleep(THREADWAIT);
	}
	bool step_queue_empty(){
		if( !randomize_step && step_queue.empty()){
			return true ;
		}
		else if (randomize_step && step_queue_vector.empty()){
			return true;
		}
		else{
			return false;
		}
	}
	int step_queue_pop()
	{
		int j ; 
		if(!randomize_step){
			j = std::move(step_queue.top());
			step_queue.pop();
		}
		else{
			int alpha = (int) ( gsl_rng_uniform(r_step) * step_queue_vector.size());
			j = std::move(step_queue_vector.at(alpha));
			step_queue_vector.erase(step_queue_vector.begin()+alpha);
		}
		return j;
	}
	void swap_pop()
	{

	}
	void pair_pop()
	{

	}

	void enqueue_swap(int i)
	{
		if (randomize_swap){
			std::unique_lock<std::mutex> lock{EventMutexSWP_pre_pair};
			swap_queue_pre_pair_vector.push_back(std::move(i));

		}
		else{
			std::unique_lock<std::mutex> lock{EventMutexSWP_pre_pair};
			swap_queue_pre_pair.emplace(std::move(i));
			
		}
		EventVarSWP_pre_pair.notify_one();
		//usleep(THREADWAIT);
	}
	

	void match_swap_pairs(int i)
	{
		bool notify = false;
		{
			std::unique_lock<std::mutex> lock{EventMutexSWP_pair};

			if(swap_queue_pair.empty()){
				swap_queue_pair.emplace_back(i);
			}
			else{
				std::vector<int>::iterator ptr;
				int k=-1;
				for(ptr=swap_queue_pair.begin(); ptr!=swap_queue_pair.end(); ++ptr){
					//std::cout<<"PTR: "<<samplerptr->chain_temps[*ptr]<<std::endl;
					if(temp_neighborhood_check(i,*ptr)){
						k = *ptr;
						swap_queue_pair.erase(ptr);
						break;
					}
				}
				if(k != -1){
					swap_struct ss;
					ss.id1 = i;	
					ss.id2 = k;	
					std::unique_lock<std::mutex> lock{EventMutexSWP_paired};
					swap_queue_paired.emplace(std::move(ss));
					notify=true;
				}
				else{
					swap_queue_pair.emplace_back(i);
				}
			}
		}
		if(notify){
			EventVarSWP_paired.notify_one();
		}
		//usleep(THREADWAIT);
	}
	void flush_swap_queue()
	{
		int size = queue_swap_length();
		for(int i = 0 ; i<size; i++){
			int l;	
			{
				std::unique_lock<std::mutex> lock{EventMutexSWP_paired};
				l = swap_queue_pair.front();
				
				swap_queue_pair.erase(swap_queue_pair.begin());
			}
			enqueue_swap(l);
		
		}
	}
	
	void public_stop()
	{
		stop();
	}
private:
	bool randomize_step = true;
	bool randomize_swap = false;

	const gsl_rng_type *T_step;
	gsl_rng *r_step;

	std::vector<std::thread> mThreads;
	
	std::condition_variable EventVarStep;

	std::mutex EventMutexStep;

	bool mStopping = false;
	
	int numSwpSWPThreads= 1;
	int numSwpSWPPThreads= 1;
		
	std::condition_variable EventVarSWP_paired;

	std::condition_variable EventVarSWP_pre_pair;

	std::mutex EventMutexSWP_paired;
	std::mutex EventMutexSWP_pair;
	std::mutex EventMutexSWP_pre_pair;

	std::priority_queue<int,std::vector<int>,Comparator> step_queue;
	std::vector<int> step_queue_vector;

	std::priority_queue<int,std::vector<int>,Comparator> swap_queue_pre_pair;
	std::vector<int> swap_queue_pre_pair_vector;
	std::vector<int> swap_queue_pair;
	std::queue<swap_struct> swap_queue_paired;
	std::vector<swap_struct> swap_queue_paired_vector;

	void start(std::size_t numThreads)
	{
		for(auto i =0u; i<numThreads-numSwpSWPThreads-numSwpSWPPThreads; i++)
		{
			mThreads.emplace_back([=, this]{
				while(true)
				{
					int j;
					{
						std::unique_lock<std::mutex> lock{EventMutexStep};

						//EventVarStep.wait(lock,[=, this]{return mStopping || !step_queue.empty() ; });
						//EventVarStep.wait(lock,[=, this]{return mStopping || ( (!step_queue.empty() && !randomize) ||(!step_queue_vector.empty() && randomize)) ; });
						EventVarStep.wait(lock,[=, this]{return mStopping || !step_queue_empty() ; });
						
						
						//if (mStopping && step_queue.empty())
						if (mStopping && step_queue_empty())
							break;	
						
						j = step_queue_pop();
						//j = std::move(step_queue.top());
						//step_queue.pop();
					}
					mcmc_step_threaded(j);
					
				}
			});
		}

		//Swapping thread
		for(auto i =0u; i<numSwpSWPThreads; i++)
		{
			mThreads.emplace_back([=, this]{
				while(true)
				{
					int j, k=-1;
					{
						std::unique_lock<std::mutex> lock{EventMutexSWP_paired};

						EventVarSWP_paired.wait(lock,[=, this]{return mStopping || !(swap_queue_paired.empty()); });
						
						if (mStopping && swap_queue_paired.empty())
							break;	
						swap_struct ss = std::move(swap_queue_paired.front());
						swap_queue_paired.pop();
						j = ss.id1;
						k = ss.id2;
					}
					mcmc_swap_threaded(j,k);
					
				}
			});
		}
		//Swapping thread
		for(auto i =0u; i<numSwpSWPPThreads; i++)
		{
			mThreads.emplace_back([=, this]{
				while(true)
				{
					int j;
					{
						std::unique_lock<std::mutex> lock{EventMutexSWP_pre_pair};

						EventVarSWP_pre_pair.wait(lock,[=, this]{return mStopping || !(swap_queue_pre_pair.empty()); });
						
						if (mStopping && swap_queue_pre_pair.empty())
							break;	
						
						j = std::move(swap_queue_pre_pair.top());
						swap_queue_pre_pair.pop();
					}
						match_swap_pairs(j);
					
				}
			});
		}
	}
	void stop() noexcept
	{
		//std::cout<<std::endl;
		//std::cout<<"Stop initiated -- waiting for threads to finish"<<std::endl;
		{
			std::unique_lock<std::mutex> lock{EventMutexStep};
			std::unique_lock<std::mutex> lock2{EventMutexSWP_paired};
			std::unique_lock<std::mutex> lock3{EventMutexSWP_pre_pair};
			mStopping = true;
		}
		
		EventVarStep.notify_all();
		EventVarSWP_paired.notify_all();
		EventVarSWP_pre_pair.notify_all();
		
		for(auto &thread: mThreads)
			thread.join();
		if(randomize_step){
			gsl_rng_free(r_step);
		}
	}
};
ThreadPool *poolptr;



//######################################################################################
//######################################################################################
void dynamic_temperature_full_ensemble_internal(sampler *samplerptr, 
	int N_steps, 
	double nu, 
	int t0,
	int swp_freq, 
	bool show_prog)
{
	//Frequency to check for equilibrium
	
	double *old_temps = new double[samplerptr->chain_N];
	int ensemble_size = 1;
	std::cout<<samplerptr->chain_temps[0]<<std::endl;
	for(int i =1; i<samplerptr->chain_N; i++){
		if(fabs(samplerptr->chain_temps[i] - 1 ) < 1e-10){break;}
		std::cout<<samplerptr->chain_temps[i]<<std::endl;
		ensemble_size ++;
	}
	for(int i =0; i<samplerptr->chain_N; i++){
		old_temps[i]=samplerptr->chain_temps[i];
	}

	//Keep track of acceptance ratio in chuncks
	int *running_accept_ct = new int[samplerptr->chain_N];
	int *running_reject_ct = new int[samplerptr->chain_N];
	int *prev_reject_ct = new int[samplerptr->chain_N];
	int *prev_accept_ct = new int[samplerptr->chain_N];
	double *running_ratio = new double[samplerptr->chain_N];
	for(int i =0; i<samplerptr->chain_N; i++){
		running_accept_ct[i] = 0;
		running_reject_ct[i] = 0;
		prev_accept_ct[i] = 0;
		prev_reject_ct[i] = 0;
		running_ratio[i] = 0;
	}

	//For each loop, we walk forward till one more swap has happened, then we update temps
	//samplerptr->N_steps = swp_freq;
	std::cout<<"Dynamical PT allocation (measured by average percent change in temperature): "<<std::endl;
	
	int t = 0;
	samplerptr->show_progress = false;
	bool testing=false;
	//step equilibrium_check_freq
	//for(int i =0; i<N_steps/samplerptr->swp_freq; i++){
	int steps  = samplerptr->swp_freq;
	samplerptr->swap_rate = 0;
	//while(t<(N_steps - samplerptr->swp_freq)){
	while(t<(N_steps - steps)){
		//Now that swapping is stochastic, we need to make sure
		//a swap actually happened
		//PTMCMC_MH_step_incremental(samplerptr, samplerptr->swp_freq);	
		PTMCMC_MH_step_incremental(samplerptr, steps);	
		//t+= samplerptr->swp_freq;
		t+= steps;
		
		//#######################################
		//Trying something here to increase mixing
		if((samplerptr->linear_swapping)){
			int swp_accepted =0, swp_rejected=0;
			chain_swap(samplerptr, samplerptr->output, samplerptr->param_status,samplerptr->model_status,t, &swp_accepted, &swp_rejected);
		}
		else
		{
			ThreadPool pool(samplerptr->num_threads);
			poolptr= &pool;
			//We're gonna try swapping every chain twice to 
			//attempt to minimize the chances a chain doesn't 
			//swap at all
			for(int j = 0 ; j<2; j++){
				std::vector<int> ids;
				for(int i = 0 ; i<samplerptr->chain_N; i++){
					ids.push_back(i);	
				}
				for(int i = 0 ; i<samplerptr->chain_N; i++){
					int alpha = (int)(gsl_rng_uniform(samplerptr->rvec[0])*ids.size());
					poolptr->enqueue_swap(ids.at(alpha));
					ids.erase(ids.begin()+alpha);
				}
			}
			poolptr->flush_swap_queue();
			int pairs = pool.queue_pairs_length();
			while(pairs > 0  ){
				pairs = pool.queue_pairs_length();
				
				usleep(1e2);
			}
		}
		//#######################################
		
		//Move temperatures
		update_temperatures_full_ensemble(samplerptr, t0, nu, t);

	}
	double acc, rej;
	double average_accept[ensemble_size];
	double average_reject[ensemble_size];
	int swap_attemps[ensemble_size];
	int ensemble_number =0;
	for(int j = 0 ; j<ensemble_size ; j++){
		average_accept[j] = 0;
		average_reject[j] = 0;
		swap_attemps[j] = 0;
	}
	for (int j =0; j<samplerptr->chain_N; j++){
		//std::cout<<"TEMP "<<j<<": "<<samplerptr->chain_temps[j]<<std::endl;
		acc = samplerptr->swap_accept_ct[j];	
		rej = samplerptr->swap_reject_ct[j];	
		average_accept[j % ensemble_size] += (double)acc/(acc+rej);
		average_reject[j % ensemble_size] += (double)rej/(acc+rej);
		swap_attemps[j%ensemble_size] += acc+rej;
		//std::cout<<"Accept ratio "<<j<<": "<<(double)acc/(acc+rej)<<std::endl;
		//std::cout<<"Swap attempts "<<j<<": "<<(acc+rej)<<std::endl;
		if(j%ensemble_size == 0){ensemble_number++;}
		
	}
	for (int j =0; j<ensemble_size; j++){
		std::cout<<"TEMP "<<j<<": "<<samplerptr->chain_temps[j]<<std::endl;
		acc = average_accept[j]/ensemble_number;	
		rej = average_reject[j]/ensemble_number;	
		std::cout<<"Accept ratio "<<j<<": "<<acc<<std::endl;
		std::cout<<"Swap Attempts "<<j<<": "<<swap_attemps[j]<<std::endl;
	}
	delete [] running_accept_ct;
	delete [] running_reject_ct;
	delete [] prev_accept_ct;
	delete [] prev_reject_ct;
	delete [] running_ratio;
	delete [] old_temps;
}
/*! \brief Continue dyanmically tunes an MCMC for optimal spacing. step width, and chain number
 *
 * NOTE: nu, and t0 parameters determine the dynamics, so these are important quantities. nu is related to how many swap attempts it takes to substantially change the temperature ladder, why t0 determines the length of the total dyanimcally period. Moderate initial choices would be 10 and 1000, respectively.
 *
 * Based on arXiv:1501.05823v3
 *
 * Currently, Chain number is fixed
 *
 * max_chain_N_thermo_ensemble sets the maximium number of chains to use to in successively hotter chains to cover the likelihood surface while targeting an optimal swap acceptance target_swp_acc. 
 *
 * max_chain_N determines the total number of chains to run once thermodynamic equilibrium has been reached. This results in chains being added after the initial PT dynamics have finished according to chain_distribution_scheme.
 *
 * If no preference, set max_chain_N_thermo_ensemble = max_chain_N = numThreads = (number of cores (number of threads if hyperthreaded))-- this will most likely be the most optimal configuration. If the number of cores on the system is low, you may want to use n*numThreads for some integer n instead, depending on the system.
 *
 * chain_distribution_scheme:
 *
 * "cold": All chains are added at T=1 (untempered)
 *
 * "refine": Chains are added between the optimal temps geometrically -- this may be a good option as it will be a good approximation of the ideal distribution of chains, while keeping the initial dynamical time low 
 *
 * "double": Chains are added in order of rising temperature that mimic the distribution achieved by the earier PT dynamics
 *
 * "half_ensemble": For every cold chain added, half of the ensemble is added again. Effectively, two cold chains for every ensemble
 */
void continue_PTMCMC_MH_dynamic_PT_alloc_full_ensemble_internal(std::string checkpoint_file_start,
	double ***output, /**< [out] Output chains, shape is double[max_chain_N, N_steps,dimension]*/
	int N_steps,	/**< Number of total steps to be taken, per chain AFTER chain allocation*/
	double *chain_temps, /**<[out] Final chain temperatures used -- should be shape double[chain_N]*/
	int swp_freq,	/**< the frequency with which chains are swapped*/
	int t0,/**< Time constant of the decay of the chain dynamics  (~1000)*/
	int nu,/**< Initial amplitude of the dynamics (~100)*/
	std::function<double(double*,int* ,int ,mcmc_data_interface *,void *)> log_prior,/**< std::function for the log_prior function -- takes double *position, int dimension, int chain_id*/
	std::function<double(double*,int*,int,mcmc_data_interface *,void *)> log_likelihood,/**< std::function for the log_likelihood function -- takes double *position, int dimension, int chain_id*/
	std::function<void(double*,int*,int,double**,mcmc_data_interface *,void *)>fisher,/**< std::function for the fisher function -- takes double *position, int dimension, double **output_fisher, int chain_id*/
	void **user_parameters,/**< Void pointer to any parameters the user may need inside log_prior, log_likelihood, or fisher. Should have one pointer for each chain. If this isn't needed, use (void**) NULL**/
	int numThreads, /**< Number of threads to use (=1 is single threaded)*/
	bool pool, /**< boolean to use stochastic chain swapping (MUST have >2 threads)*/
	bool show_prog, /**< boolean whether to print out progress (for example, should be set to ``false'' if submitting to a cluster)*/
	std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
	std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output*/
	std::string checkpoint_file,/**< Filename to output data for checkpoint, if empty string, not saved*/
	bool burn_phase 
	)
{
	//std::cout<<"MEM CHECK : start continue"<<std::endl;
	clock_t start, end, acend;
	double wstart, wend, wacend;
	start = clock();
	wstart = omp_get_wtime();

	sampler samplerobj;
	sampler *samplerptr = &samplerobj;
	samplerptr_global = &samplerobj;

	//if Fisher is not provided, Fisher and MALA steps
	//aren't used
	if(fisher ==NULL){
		samplerptr->fisher_exist = false;
	}
	else 
		samplerptr->fisher_exist = true;
	
	samplerptr->log_ll = false;
	samplerptr->log_lp = false;
	
	//Construct sampler structure
	samplerptr->lp = log_prior;
	samplerptr->ll = log_likelihood;
	samplerptr->fish = fisher;
	samplerptr->swp_freq = 2.;
	samplerptr->swap_rate = 1./swp_freq;
	//For PT dynamics
	samplerptr->N_steps = N_steps;

	samplerptr->num_threads = numThreads;
	samplerptr->output =output;
	samplerptr->user_parameters=user_parameters;
	samplerptr->burn_phase = burn_phase;
	samplerptr->tune = false;

	//###############################################
	//We can use pooling, but the swap radius MUST be 1
	//The equations for the temperature dynamics are written
	//for neighboring temperatures only
	//
	//Linear swapping averages temperatures (logarithmically) 
	//at the end. Pooling averages at each step
	samplerptr->chain_radius = 1;
	//samplerptr->pool = true;
	samplerptr->pool = false;
	samplerptr->linear_swapping = true;
	//###############################################

	load_checkpoint_file(checkpoint_file_start,samplerptr);


	////###############################################
	////We can use pooling, but the swap radius MUST be 1
	////The equations for the temperature dynamics are written
	////for neighboring temperatures only
	////
	////This is being too problematic for now
	//samplerptr->chain_radius = 1;
	////samplerptr->pool = true;
	//samplerptr->pool = false;
	////###############################################

	samplerptr->numThreads = numThreads;
	samplerptr->A = new int[samplerptr->chain_N];
	for(int i =0 ; i<samplerptr->chain_N; i++){
		samplerptr->A[i]=0;
	}
	samplerptr->PT_alloc = true;
	

	//samplerptr->chain_N = ;//For allocation purposes, this needs to be the maximium number of chains
	//allocate_sampler_mem(samplerptr);
	

	for (int chain_index=0; chain_index<samplerptr->chain_N; chain_index++)
		assign_probabilities(samplerptr, chain_index);
	
	for (int j=0;j<samplerptr->chain_N;j++){
		samplerptr->current_likelihoods[j] =
			 samplerptr->ll(samplerptr->output[j][0],samplerptr->param_status[j][0],samplerptr->model_status[j][0],samplerptr->interfaces[j],samplerptr->user_parameters[j]);
		//std::cout<<samplerptr->current_likelihoods[j]<<std::endl;
		//step_accepted[j]=0;
		//step_rejected[j]=0;
	}
	
	//Set chains with temp 1 to highest priority
	if(samplerptr->prioritize_cold_chains){
		for(int i =0 ;i<samplerptr->chain_N; i++){
			if(fabs(samplerptr->chain_temps[i]-1)<DOUBLE_COMP_THRESH)
				samplerptr->priority[i] = 0;
		}
	}
	
	//#########################################################################
	//#########################################################################
	//#########################################################################
	//#########################################################################
	//#########################################################################
	
	//std::cout<<"MEM CHECK : start loop allocation"<<std::endl;
	dynamic_temperature_full_ensemble_internal(samplerptr, N_steps, nu, t0,swp_freq,  show_prog);

	//std::cout<<"MEM CHECK : start memory allocation"<<std::endl;
	//#######################################################################
	//#######################################################################
	//#######################################################################
	//#######################################################################
	
	end =clock();
	wend =omp_get_wtime();

	samplerptr->time_elapsed_cpu = (double)(end-start)/CLOCKS_PER_SEC;
	samplerptr->time_elapsed_wall = (double)(wend-wstart);

	
	acend =clock();
	wacend =omp_get_wtime();
	samplerptr->time_elapsed_cpu_ac = (double)(acend-end)/CLOCKS_PER_SEC;
	samplerptr->time_elapsed_wall_ac = (double)(wacend - wend);

	//############################################################
	//#################################################################
	//
	//Replace temperatures with averages
	//Just comment out if strictly identical ensembles are not required
	//
	//#################################################################
	
	double chains_per_ensemble = 1;
	bool ensemble_edge_found = false;
	int ensemble_number = 1 ;
	for(int i = 1 ; i < samplerptr->chain_N;i++){
		if(fabs(samplerptr->chain_temps[i] - 1)>DOUBLE_COMP_THRESH ){
			if(!ensemble_edge_found){
				chains_per_ensemble++;
			}
		}	
		else{
			ensemble_edge_found = true;
			ensemble_number++;
		}
	}
	int n_temp_groups = std::max(0, (int)(chains_per_ensemble - 2));
	std::vector<double> chain_temp_averages(n_temp_groups, 0.);
	std::vector<int> chain_N_per_group(n_temp_groups, 0);
	double temp_upper_bound = samplerptr->chain_temps[(int)(chains_per_ensemble-1)];
	int ct = 0 ;
	for(int i = 1 ; i < samplerptr->chain_N;i++){
		if(fabs(samplerptr->chain_temps[i]-1) < DOUBLE_COMP_THRESH ||
			fabs(samplerptr->chain_temps[i] - temp_upper_bound) < DOUBLE_COMP_THRESH){
			ct=0;
			continue;
		}
		else{
			//chain_temp_averages[ct] +=samplerptr->chain_temps[i];
			chain_temp_averages[ct] +=std::log(samplerptr->chain_temps[i]);
			chain_N_per_group[ct]+=1;
			ct++;

		}
	}
	for(int i = 0 ; i < n_temp_groups;i++){
		chain_temp_averages[i]/=chain_N_per_group[i];
		chain_temp_averages[i] = std::exp(chain_temp_averages[i]);
		//chain_temp_averages[i]=std::pow(chain_temp_averages[i],1./chain_N_per_group[i]);
		//debugger_print(__FILE__,__LINE__,chain_temp_averages[i]);
	}
	ct = 0;
	chain_temps[0]= 1;
	for(int i = 1 ; i<samplerptr->chain_N ; i++){
		if(fabs(samplerptr->chain_temps[i]-1) < DOUBLE_COMP_THRESH || 
			fabs(samplerptr->chain_temps[i] - temp_upper_bound) < DOUBLE_COMP_THRESH){
			ct=0;
			continue;
		}
		else{
			//samplerptr->current_likelihoods[i]*=samplerptr->chain_temps[i];
			samplerptr->chain_temps[i] = chain_temp_averages[ct];
			//samplerptr->current_likelihoods[i]/=samplerptr->chain_temps[i];
			ct++;
		}
		chain_temps[i]= samplerptr->chain_temps[i];
	}
	//###########################################################
	//###########################################################
	
	
	

	if(statistics_filename != "")
		write_stat_file(samplerptr, statistics_filename);
	write_checkpoint_file(samplerptr, checkpoint_file);

	delete [] samplerptr->A;
	deallocate_sampler_mem(samplerptr);
	free(samplerptr->chain_temps);

}
//######################################################################################
//######################################################################################
/*! \brief Continue dyanmically tunes an MCMC for optimal spacing. step width, and chain number
 *
 * NOTE: nu, and t0 parameters determine the dynamics, so these are important quantities. nu is related to how many swap attempts it takes to substantially change the temperature ladder, why t0 determines the length of the total dyanimcally period. Moderate initial choices would be 10 and 1000, respectively.
 *
 * Based on arXiv:1501.05823v3
 *
 * Currently, Chain number is fixed
 *
 * max_chain_N_thermo_ensemble sets the maximium number of chains to use to in successively hotter chains to cover the likelihood surface while targeting an optimal swap acceptance target_swp_acc. 
 *
 * max_chain_N determines the total number of chains to run once thermodynamic equilibrium has been reached. This results in chains being added after the initial PT dynamics have finished according to chain_distribution_scheme.
 *
 * If no preference, set max_chain_N_thermo_ensemble = max_chain_N = numThreads = (number of cores (number of threads if hyperthreaded))-- this will most likely be the most optimal configuration. If the number of cores on the system is low, you may want to use n*numThreads for some integer n instead, depending on the system.
 *
 * chain_distribution_scheme:
 *
 * "cold": All chains are added at T=1 (untempered)
 *
 * "refine": Chains are added between the optimal temps geometrically -- this may be a good option as it will be a good approximation of the ideal distribution of chains, while keeping the initial dynamical time low 
 *
 * "double": Chains are added in order of rising temperature that mimic the distribution achieved by the earier PT dynamics
 *
 * "half_ensemble": For every cold chain added, half of the ensemble is added again. Effectively, two cold chains for every ensemble
 */
void continue_RJPTMCMC_MH_dynamic_PT_alloc_full_ensemble_internal(std::string checkpoint_file_start,
	double ***output, /**< [out] Output chains, shape is double[max_chain_N, N_steps,dimension]*/
	int ***status,/**< [out] output parameter status array, dimensions: status[chain_N][N_steps][dimension]*/
	int **model_status,
	int nested_model_number,
	int N_steps,	/**< Number of total steps to be taken, per chain AFTER chain allocation*/
	double *chain_temps, /**<[out] Final chain temperatures used -- should be shape double[chain_N]*/
	int swp_freq,	/**< the frequency with which chains are swapped*/
	int t0,/**< Time constant of the decay of the chain dynamics  (~1000)*/
	int nu,/**< Initial amplitude of the dynamics (~100)*/
	std::function<double(double*, int*,int,mcmc_data_interface *,void *)> log_prior,/**< std::function for the log_prior function -- takes double *position, int *param_status, int dimension, int chain_id*/
	std::function<double(double*,int*,int,mcmc_data_interface *,void*)> log_likelihood,/**< std::function for the log_likelihood function -- takes double *position, int *param_status,int dimension, int chain_id*/
	std::function<void(double*,int*,int,double**,mcmc_data_interface *,void*)>fisher,/**< std::function for the fisher function -- takes double *position, int *param_status,int dimension, double **output_fisher, int chain_id*/
	std::function<void(double*,double*,int*,int*,int*,int*,double *,mcmc_data_interface *,void*)> RJ_proposal,/**< std::function for the log_likelihood function -- takes double *position, int *param_status,int dimension, int chain_id*/
	void **user_parameters,/**< Void pointer to any parameters the user may need inside log_prior, log_likelihood, or fisher. Should have one pointer for each chain. If this isn't needed, use (void**) NULL**/
	int numThreads, /**< Number of threads to use (=1 is single threaded)*/
	bool pool, /**< boolean to use stochastic chain swapping (MUST have >2 threads)*/
	bool show_prog, /**< boolean whether to print out progress (for example, should be set to ``false'' if submitting to a cluster)*/
	bool update_RJ_width, 
	std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
	std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output*/
	std::string checkpoint_file,/**< Filename to output data for checkpoint, if empty string, not saved*/
	bool burn_phase 
	)
{
	//std::cout<<"MEM CHECK : start continue"<<std::endl;
	clock_t start, end, acend;
	double wstart, wend, wacend;
	start = clock();
	wstart = omp_get_wtime();

	sampler samplerobj;
	sampler *samplerptr = &samplerobj;
	samplerptr_global = &samplerobj;

	samplerptr->update_RJ_width = update_RJ_width;

	//if Fisher is not provided, Fisher and MALA steps
	//aren't used
	if(fisher ==NULL){
		samplerptr->fisher_exist = false;
	}
	else 
		samplerptr->fisher_exist = true;
	
	samplerptr->log_ll = false;
	samplerptr->log_lp = false;


	
	//Construct sampler structure
	samplerptr->lp = log_prior;
	samplerptr->ll = log_likelihood;
	samplerptr->fish = fisher;
	samplerptr->rj = RJ_proposal;
	samplerptr->swp_freq = 2.;
	samplerptr->swap_rate = 1./swp_freq;
	//For PT dynamics
	samplerptr->N_steps = N_steps;

	samplerptr->nested_model_number = nested_model_number;
	samplerptr->model_status = model_status;
	samplerptr->RJMCMC = true;

	samplerptr->num_threads = numThreads;
	samplerptr->output =output;
	samplerptr->param_status= status;
	samplerptr->user_parameters=user_parameters;
	samplerptr->burn_phase = burn_phase;
	samplerptr->tune = false;

	//###############################################
	//We can use pooling, but the swap radius MUST be 1
	//The equations for the temperature dynamics are written
	//for neighboring temperatures only
	//
	//Linear swapping averages temperatures (logarithmically) 
	//at the end. Pooling averages at each step
	samplerptr->chain_radius = 1;
	//samplerptr->pool = true;
	samplerptr->pool = false;
	samplerptr->linear_swapping = true;
	//###############################################

	load_checkpoint_file(checkpoint_file_start,samplerptr);


	////###############################################
	////We can use pooling, but the swap radius MUST be 1
	////The equations for the temperature dynamics are written
	////for neighboring temperatures only
	////
	////This is being too problematic for now
	//samplerptr->chain_radius = 1;
	////samplerptr->pool = true;
	//samplerptr->pool = false;
	////###############################################

	samplerptr->numThreads = numThreads;
	samplerptr->A = new int[samplerptr->chain_N];
	for(int i =0 ; i<samplerptr->chain_N; i++){
		samplerptr->A[i]=0;
	}
	samplerptr->PT_alloc = true;
	

	//samplerptr->chain_N = ;//For allocation purposes, this needs to be the maximium number of chains
	//allocate_sampler_mem(samplerptr);
	

	for (int chain_index=0; chain_index<samplerptr->chain_N; chain_index++)
		assign_probabilities(samplerptr, chain_index);
	
	for (int j=0;j<samplerptr->chain_N;j++){
		samplerptr->current_likelihoods[j] =
			 samplerptr->ll(samplerptr->output[j][0],samplerptr->param_status[j][0],samplerptr->model_status[j][0],samplerptr->interfaces[j],samplerptr->user_parameters[j]);
		//std::cout<<samplerptr->current_likelihoods[j]<<std::endl;
		//step_accepted[j]=0;
		//step_rejected[j]=0;
	}
	if(samplerptr->log_ll && samplerptr->log_lp){
		for(int i = 0 ; i<samplerptr->chain_N; i++){
			samplerptr->ll_lp_output[i][0][0] = samplerptr->current_likelihoods[i];
			samplerptr->ll_lp_output[i][0][1] = samplerptr->lp(samplerptr->output[i][0],samplerptr->param_status[i][0],samplerptr->model_status[i][0],samplerptr->interfaces[i],samplerptr->user_parameters[i]);
		}
	}
	
	//Set chains with temp 1 to highest priority
	if(samplerptr->prioritize_cold_chains){
		for(int i =0 ;i<samplerptr->chain_N; i++){
			if(fabs(samplerptr->chain_temps[i]-1)<DOUBLE_COMP_THRESH)
				samplerptr->priority[i] = 0;
		}
	}
	
	//#########################################################################
	//#########################################################################
	//#########################################################################
	//#########################################################################
	//#########################################################################
	
	//std::cout<<"MEM CHECK : start loop allocation"<<std::endl;
	dynamic_temperature_full_ensemble_internal(samplerptr, N_steps, nu, t0,swp_freq,  show_prog);

	//std::cout<<"MEM CHECK : start memory allocation"<<std::endl;
	//#######################################################################
	//#######################################################################
	//#######################################################################
	//#######################################################################
	
	end =clock();
	wend =omp_get_wtime();

	samplerptr->time_elapsed_cpu = (double)(end-start)/CLOCKS_PER_SEC;
	samplerptr->time_elapsed_wall = (double)(wend-wstart);

	
	acend =clock();
	wacend =omp_get_wtime();
	samplerptr->time_elapsed_cpu_ac = (double)(acend-end)/CLOCKS_PER_SEC;
	samplerptr->time_elapsed_wall_ac = (double)(wacend - wend);

	//############################################################
	//#################################################################
	//
	//Replace temperatures with averages
	//Just comment out if strictly identical ensembles are not required
	//
	//#################################################################
	
	double chains_per_ensemble = 1;
	bool ensemble_edge_found = false;
	int ensemble_number = 1 ;
	for(int i = 1 ; i < samplerptr->chain_N;i++){
		if(fabs(samplerptr->chain_temps[i] - 1)>DOUBLE_COMP_THRESH ){
			if(!ensemble_edge_found){
				chains_per_ensemble++;
			}
		}	
		else{
			ensemble_edge_found = true;
			ensemble_number++;
		}
	}
	int n_temp_groups = std::max(0, (int)(chains_per_ensemble - 2));
	std::vector<double> chain_temp_averages(n_temp_groups, 0.);
	std::vector<int> chain_N_per_group(n_temp_groups, 0);
	double temp_upper_bound = samplerptr->chain_temps[(int)(chains_per_ensemble-1)];
	int ct = 0 ;
	for(int i = 1 ; i < samplerptr->chain_N;i++){
		if(fabs(samplerptr->chain_temps[i]-1) < DOUBLE_COMP_THRESH ||
			fabs(samplerptr->chain_temps[i] - temp_upper_bound) < DOUBLE_COMP_THRESH){
			ct=0;
			continue;
		}
		else{
			//chain_temp_averages[ct] +=samplerptr->chain_temps[i];
			chain_temp_averages[ct] +=std::log(samplerptr->chain_temps[i]);
			chain_N_per_group[ct]+=1;
			ct++;

		}
	}
	for(int i = 0 ; i < n_temp_groups;i++){
		chain_temp_averages[i]/=chain_N_per_group[i];
		chain_temp_averages[i] = std::exp(chain_temp_averages[i]);
		//chain_temp_averages[i]=std::pow(chain_temp_averages[i],1./chain_N_per_group[i]);
		debugger_print(__FILE__,__LINE__,chain_temp_averages[i]);
	}
	ct = 0;
	chain_temps[0]= 1;
	for(int i = 1 ; i<samplerptr->chain_N ; i++){
		if(fabs(samplerptr->chain_temps[i]-1) < DOUBLE_COMP_THRESH || 
			fabs(samplerptr->chain_temps[i] - temp_upper_bound) < DOUBLE_COMP_THRESH){
			ct=0;
			continue;
		}
		else{
			samplerptr->chain_temps[i] = chain_temp_averages[ct];
			ct++;
		}
		chain_temps[i]= samplerptr->chain_temps[i];
	}
	//###########################################################
	//###########################################################

	if(statistics_filename != "")
		write_stat_file(samplerptr, statistics_filename);
	write_checkpoint_file(samplerptr, checkpoint_file);

	delete [] samplerptr->A;
	free(samplerptr->chain_temps);
	deallocate_sampler_mem(samplerptr);

}
//######################################################################################
//######################################################################################





void fisher_generic(double* position,int* status,double **fish,mcmc_data_interface *interface,void *, sampler *sampler)
{
	int dim = interface->max_dim;
	int chain_id = interface->chain_number;
	int current_pos = sampler->chain_pos[chain_id];
	double **cov = allocate_2D_array(dim,dim);
	if(current_pos < 10){
		for(int i = 0 ; i<dim ; i++){
			for(int j = 0 ; j<dim; j++){
				fish[i][j]=0;
			}
			fish[i][i]=1;
		}
		return;
		
	}

	double means[dim];
	for(int i = 0 ; i<dim ; i++){
		means[i]=0;
		for(int j = 0 ; j<current_pos; j++){
			means[i]+=sampler->output[chain_id][j][i];
		}
		means[i]/=current_pos;
	}
	
	for(int i = 0 ; i<dim ; i++){
		for(int j = 0 ; j<=i; j++){
			cov[i][j]=0;
			for(int k = 0 ; k<current_pos; k++){
				cov[i][j]+=(sampler->output[chain_id][k][i] - means[i])*
					(sampler->output[chain_id][k][j]-means[j]);
			}
			cov[i][j]/=current_pos;
		}
	}	
	for(int i=0 ; i<dim; i++){
		for(int j = i ; j<dim; j++){
			cov[i][j] = cov[j][i];
		}
	}
	gsl_cholesky_matrix_invert(cov,fish,dim);
	deallocate_2D_array(cov,dim,dim);
}

//######################################################################################
//######################################################################################
/*! \brief Routine to take a checkpoint file and begin a new chain at said checkpoint
 *
 * See MCMC_MH_internal for more details of parameters (pretty much all the same)
 */
void continue_RJPTMCMC_MH_simulated_annealing_internal(sampler *sampler_in,
	std::string start_checkpoint_file,/**< File for starting checkpoint*/
	double ***output,/**< [out] output array, dimensions: output[chain_N][N_steps][dimension]*/
	int ***status,/**< [out] output parameter status array, dimensions: status[chain_N][N_steps][dimension]*/
	int **model_status,
	int nested_model_number,
	int N_steps,/**< Number of new steps to take*/
	double temp_scale_factor,
	int swp_freq,/**< frequency of swap attempts between temperatures*/
	int t0,/**< Time constant of the decay of the chain dynamics  (~1000)*/
	int nu,/**< Initial amplitude of the dynamics (~100)*/
	std::function<double(double*,int*,int,mcmc_data_interface *,void *)> log_prior,/**< std::function for the log_prior function -- takes double *position, int dimension, int chain_id*/
	std::function<double(double*,int*,int,mcmc_data_interface *,void *)> log_likelihood,/**< std::function for the log_likelihood function -- takes double *position, int dimension, int chain_id*/
	std::function<void(double*,int*,int,double**,mcmc_data_interface *,void *)>fisher,/**< std::function for the fisher function -- takes double *position, int dimension, double **output_fisher, int chain_id*/
	std::function<void(double*,double*,int*,int*,int*,int*,double *,mcmc_data_interface *,void*)> RJ_proposal,/**< std::function for the log_likelihood function -- takes double *position, int *param_status,int dimension, int chain_id*/
	void **user_parameters,/**< Void pointer to any parameters the user may need inside log_prior, log_likelihood, or fisher. Should have one pointer for each chain. If this isn't needed, use (void**) NULL**/
	int numThreads,/**<Number of threads to use*/
	bool pool,/**<Boolean for whether to use ``deterministic'' vs ``stochastic'' sampling*/
	bool show_prog,/**< Boolean for whether to show progress or not (turn off for cluster runs*/
	bool update_RJ_width, 
	std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
	std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output*/
	std::string end_checkpoint_file/**< Filename to output data for checkpoint at the end of the continued run, if empty string, not saved*/
	)
{
	clock_t start, end, acend;
	double wstart, wend, wacend;
	start = clock();
	wstart = omp_get_wtime();

	//sampler sampler;
	sampler *samplerptr = sampler_in;
	samplerptr_global = sampler_in;
	samplerptr->restrict_swapping=false;	
	samplerptr->tune=true;
	samplerptr->burn_phase = true;
	samplerptr->isolate_ensembles = false;
	samplerptr->isolate_ensembles_cold = false;
	samplerptr->update_RJ_width=update_RJ_width;
	samplerptr->rj = RJ_proposal;
	samplerptr->param_status= status;
	samplerptr->model_status = model_status;
	samplerptr->nested_model_number = nested_model_number;
	samplerptr->RJMCMC = true;
	//################################################
	//This typically isn't done, but the primary focus of annealing
	//is to get the cold chains where they need to be.
	//The hotter chains equilibrate during normal operation
	//much easier and don't need the annealing process as much
	//samplerptr->prioritize_cold_chains=true;
	samplerptr->prioritize_cold_chains=false;
	//################################################

	//if Fisher is not provided, Fisher and MALA steps
	//aren't used
	if(fisher ==NULL){
		samplerptr->fisher_exist = false;
		//samplerptr->fisher_exist = true;
		//fisher = [](double *param, int*param_status,  double **fish, mcmc_data_interface *interface, void *parameters){
		//	return fisher_generic(param,param_status,fish, interface, parameters, samplerptr);	
		//};
	}
	else {
		samplerptr->fisher_exist = true;
	}

	samplerptr->log_ll = true;
	samplerptr->log_lp = true;
	
	//Construct sampler structure
	samplerptr->lp = log_prior;
	samplerptr->ll = log_likelihood;
	samplerptr->fish = fisher;
	samplerptr->swp_freq = 2;
	samplerptr->swap_rate = 1./swp_freq;
	samplerptr->N_steps = N_steps;
	samplerptr->show_progress = show_prog;
	samplerptr->num_threads = numThreads;
	samplerptr->user_parameters=user_parameters;


	samplerptr->output = output;

	//samplerptr->pool = pool;
	samplerptr->pool = false;

	//Unpack checkpoint file -- allocates memory internally -- separate call unneccessary
	load_checkpoint_file(start_checkpoint_file, samplerptr);


	//allocate other parameters
	for (int chain_index=0; chain_index<samplerptr->chain_N; chain_index++)
		assign_probabilities(samplerptr, chain_index);
	for (int j=0;j<samplerptr->chain_N;j++){
		samplerptr->current_likelihoods[j] =
			 samplerptr->ll(samplerptr->output[j][0],samplerptr->param_status[j][0],samplerptr->model_status[j][0],samplerptr->interfaces[j],samplerptr->user_parameters[j]);
	}
	
	//Set chains with temp 1 to highest priority
	if(samplerptr->prioritize_cold_chains){
		for(int i =0 ;i<samplerptr->chain_N; i++){
			if(fabs(samplerptr->chain_temps[i]-1)<DOUBLE_COMP_THRESH)
				samplerptr->priority[i] = 0;
		}
	}
	
	//########################################################
	//########################################################
	double temp_save[samplerptr->chain_N];
	int temp_step_N = (N_steps<100)? N_steps/5: 100;
	double temp_steps[samplerptr->chain_N];
	for(int i = 0 ; i<samplerptr->chain_N; i++){
		//samplerptr->current_likelihoods[i]*=samplerptr->chain_temps[i];
		temp_save[i]=samplerptr->chain_temps[i];
		double T_max = temp_scale_factor*samplerptr->chain_temps[i];
		double T_min = samplerptr->chain_temps[i];
		temp_steps[i] = pow(T_min/T_max,1./temp_step_N);
		samplerptr->chain_temps[i]=T_max;
	}
	int steps_per_temp = N_steps/temp_step_N;
		
	for(int i = 0 ; i<temp_step_N ; i++){
		for(int j = 0 ; j<samplerptr->chain_N; j++){
			samplerptr->chain_temps[j] = temp_steps[j]*samplerptr->chain_temps[j];
			//samplerptr->current_likelihoods[j]/=samplerptr->chain_temps[j];
			//if(j == 0 ){debugger_print(__FILE__,__LINE__,samplerptr->chain_temps[j]);}
		}
			
		PTMCMC_MH_step_incremental(samplerptr, steps_per_temp);
		//for(int j = 0 ; j<samplerptr->chain_N; j++){
		//	samplerptr->current_likelihoods[j]*=samplerptr->chain_temps[j];
		//}
	}
	for(int i = 0 ; i<samplerptr->chain_N; i++){
		samplerptr->chain_temps[i]=temp_save[i];
	}
	//##############################################################
	//##############################################################
	
	end =clock();
	wend =omp_get_wtime();

	samplerptr->time_elapsed_cpu = (double)(end-start)/CLOCKS_PER_SEC;
	samplerptr->time_elapsed_wall = (double)(wend-wstart);

	if(show_prog){
		std::cout<<std::endl;
	}
	//############################################################
	
	//###########################################################
	//###########################################################
	acend =clock();
	wacend =omp_get_wtime();
	samplerptr->time_elapsed_cpu_ac = (double)(acend-end)/CLOCKS_PER_SEC;
	samplerptr->time_elapsed_wall_ac = (double)(wacend - wend);

	
	if(statistics_filename != "")
		write_stat_file(samplerptr, statistics_filename);
	
	if(chain_filename != "")
		write_output_file(chain_filename, samplerptr->N_steps, samplerptr->max_dim, samplerptr->output, samplerptr->param_status,samplerptr->model_status,samplerptr->chain_N,samplerptr->nested_model_number,samplerptr->chain_temps,false);

	if(end_checkpoint_file !=""){
		write_checkpoint_file(samplerptr, end_checkpoint_file);
	}

	//temps usually allocated by user, but for continued chains, this is done internally
	free(samplerptr->chain_temps);
}

//######################################################################################
/*! \brief Routine to take a checkpoint file and begin a new chain at said checkpoint
 *
 * See MCMC_MH_internal for more details of parameters (pretty much all the same)
 */
void continue_PTMCMC_MH_simulated_annealing_internal(sampler *sampler_in,
	std::string start_checkpoint_file,/**< File for starting checkpoint*/
	double ***output,/**< [out] output array, dimensions: output[chain_N][N_steps][dimension]*/
	int N_steps,/**< Number of new steps to take*/
	double temp_scale_factor,
	int swp_freq,/**< frequency of swap attempts between temperatures*/
	int t0,/**< Time constant of the decay of the chain dynamics  (~1000)*/
	int nu,/**< Initial amplitude of the dynamics (~100)*/
	std::function<double(double*,int*,int,mcmc_data_interface *,void *)> log_prior,/**< std::function for the log_prior function -- takes double *position, int dimension, int chain_id*/
	std::function<double(double*,int*,int,mcmc_data_interface *,void *)> log_likelihood,/**< std::function for the log_likelihood function -- takes double *position, int dimension, int chain_id*/
	std::function<void(double*,int*,int,double**,mcmc_data_interface *,void *)>fisher,/**< std::function for the fisher function -- takes double *position, int dimension, double **output_fisher, int chain_id*/
	void **user_parameters,/**< Void pointer to any parameters the user may need inside log_prior, log_likelihood, or fisher. Should have one pointer for each chain. If this isn't needed, use (void**) NULL**/
	int numThreads,/**<Number of threads to use*/
	bool pool,/**<Boolean for whether to use ``deterministic'' vs ``stochastic'' sampling*/
	bool show_prog,/**< Boolean for whether to show progress or not (turn off for cluster runs*/
	std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
	std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output*/
	std::string end_checkpoint_file/**< Filename to output data for checkpoint at the end of the continued run, if empty string, not saved*/
	)
{
	clock_t start, end, acend;
	double wstart, wend, wacend;
	start = clock();
	wstart = omp_get_wtime();

	//sampler sampler;
	sampler *samplerptr = sampler_in;
	samplerptr_global = sampler_in;
	samplerptr->restrict_swapping=false;	
	samplerptr->tune=true;
	samplerptr->burn_phase = true;
	samplerptr->isolate_ensembles = false;
	samplerptr->isolate_ensembles_cold = false;
	//################################################
	//This typically isn't done, but the primary focus of annealing
	//is to get the cold chains where they need to be.
	//The hotter chains equilibrate during normal operation
	//much easier and don't need the annealing process as much
	//samplerptr->prioritize_cold_chains=true;
	samplerptr->prioritize_cold_chains=false;
	//################################################

	//if Fisher is not provided, Fisher and MALA steps
	//aren't used
	if(fisher ==NULL){
		samplerptr->fisher_exist = false;
		//samplerptr->fisher_exist = true;
		//fisher = [](double *param, int*param_status,  double **fish, mcmc_data_interface *interface, void *parameters){
		//	return fisher_generic(param,param_status,fish, interface, parameters, samplerptr);	
		//};
	}
	else {
		samplerptr->fisher_exist = true;
	}

	samplerptr->log_ll = true;
	samplerptr->log_lp = true;
	
	//Construct sampler structure
	samplerptr->lp = log_prior;
	samplerptr->ll = log_likelihood;
	samplerptr->fish = fisher;
	samplerptr->swp_freq = 2;
	samplerptr->swap_rate = 1./swp_freq;
	samplerptr->N_steps = N_steps;
	samplerptr->show_progress = show_prog;
	samplerptr->num_threads = numThreads;
	samplerptr->user_parameters=user_parameters;


	samplerptr->output = output;
	//samplerptr->pool = pool;
	samplerptr->pool = false;

	//Unpack checkpoint file -- allocates memory internally -- separate call unneccessary
	load_checkpoint_file(start_checkpoint_file, samplerptr);


	//allocate other parameters
	for (int chain_index=0; chain_index<samplerptr->chain_N; chain_index++)
		assign_probabilities(samplerptr, chain_index);
	for (int j=0;j<samplerptr->chain_N;j++){
		samplerptr->current_likelihoods[j] =
			 samplerptr->ll(samplerptr->output[j][0],samplerptr->param_status[j][0],samplerptr->model_status[j][0],samplerptr->interfaces[j],samplerptr->user_parameters[j]);
	}
	
	//Set chains with temp 1 to highest priority
	if(samplerptr->prioritize_cold_chains){
		for(int i =0 ;i<samplerptr->chain_N; i++){
			if(fabs(samplerptr->chain_temps[i]-1)<DOUBLE_COMP_THRESH)
				samplerptr->priority[i] = 0;
		}
	}
	
	//########################################################
	//########################################################
	double temp_save[samplerptr->chain_N];
	int temp_step_N = (N_steps<100)? N_steps/5: 100;
	double temp_steps[samplerptr->chain_N];
	for(int i = 0 ; i<samplerptr->chain_N; i++){
		//samplerptr->current_likelihoods[i]*=samplerptr->chain_temps[i];
		temp_save[i]=samplerptr->chain_temps[i];
		double T_max = temp_scale_factor*samplerptr->chain_temps[i];
		double T_min = samplerptr->chain_temps[i];
		temp_steps[i] = pow(T_min/T_max,1./temp_step_N);
		samplerptr->chain_temps[i]=T_max;
	}
	int steps_per_temp = N_steps/temp_step_N;
		
	for(int i = 0 ; i<temp_step_N ; i++){
		for(int j = 0 ; j<samplerptr->chain_N; j++){
			samplerptr->chain_temps[j] = temp_steps[j]*samplerptr->chain_temps[j];
			//samplerptr->current_likelihoods[j]/=samplerptr->chain_temps[j];
			//if(j == 0 ){debugger_print(__FILE__,__LINE__,samplerptr->chain_temps[j]);}
		}
			
		PTMCMC_MH_step_incremental(samplerptr, steps_per_temp);
		//for(int j = 0 ; j<samplerptr->chain_N; j++){
		//	samplerptr->current_likelihoods[j]*=samplerptr->chain_temps[j];
		//}
	}
	for(int i = 0 ; i<samplerptr->chain_N; i++){
		samplerptr->chain_temps[i]=temp_save[i];
	}
	//##############################################################
	//##############################################################
	
	end =clock();
	wend =omp_get_wtime();

	samplerptr->time_elapsed_cpu = (double)(end-start)/CLOCKS_PER_SEC;
	samplerptr->time_elapsed_wall = (double)(wend-wstart);

	if(show_prog){
		std::cout<<std::endl;
	}
	//############################################################
	
	//###########################################################
	//###########################################################
	acend =clock();
	wacend =omp_get_wtime();
	samplerptr->time_elapsed_cpu_ac = (double)(acend-end)/CLOCKS_PER_SEC;
	samplerptr->time_elapsed_wall_ac = (double)(wacend - wend);

	
	if(statistics_filename != "")
		write_stat_file(samplerptr, statistics_filename);
	
	if(chain_filename != "")
		write_output_file(chain_filename, samplerptr->N_steps, samplerptr->max_dim, samplerptr->output, samplerptr->param_status,samplerptr->model_status,samplerptr->chain_N,samplerptr->nested_model_number,samplerptr->chain_temps,false);

	if(end_checkpoint_file !=""){
		write_checkpoint_file(samplerptr, end_checkpoint_file);
	}

	//temps usually allocated by user, but for continued chains, this is done internally
	free(samplerptr->chain_temps);
}

//######################################################################################
//######################################################################################



/*! \brief Parallel tempered, dynamic chain allocation MCMC. Reversible Jump version, which doesn't currently have any autocorrelation calculations currently.
 *
 * Runs dynamic chain allocation, then it will run the RJPTMCMC routine repeatedly 
 *
 */
void continue_RJPTMCMC_MH_dynamic_PT_alloc_comprehensive_internal(
	std::string checkpoint_file_start, 
	mcmc_sampler_output *sampler_output,
	double **output, /**< [out] Output shape is double[N_steps,dimension]*/
	int **parameter_status, /**< [out] Output shape is int[N_steps,dimension]*/
	int *model_status,
	int nested_model_number,
	int N_steps,	/**< Number of total steps to be taken, per chain AFTER chain allocation*/
	int max_chain_N_thermo_ensemble,/**< Maximum number of chains to use in the thermodynamic ensemble (may use less)*/
	double *chain_temps, /**<[out] Final chain temperatures used -- should be shape double[chain_N]*/
	int swp_freq,	/**< the frequency with which chains are swapped*/
	int t0,/**< Time constant of the decay of the chain dynamics  (~1000)*/
	int nu,/**< Initial amplitude of the dynamics (~100)*/
	int max_chunk_size,/**<Maximum number of steps to take in a single sampler run*/
	std::string chain_distribution_scheme,
	std::function<double(double*, int*,int,mcmc_data_interface *,void *)> log_prior,/**< std::function for the log_prior function -- takes double *position, int *param_status, int dimension, int chain_id*/
	std::function<double(double*,int*,int,mcmc_data_interface *,void *)> log_likelihood,/**< std::function for the log_likelihood function -- takes double *position, int *param_status,int dimension, int chain_id*/
	std::function<void(double*,int*,int,double**,mcmc_data_interface *,void *)>fisher,/**< std::function for the fisher function -- takes double *position, int *param_status,int dimension, double **output_fisher, int chain_id*/
	std::function<void(double*,double*,int*,int*,int*,int*,double *,mcmc_data_interface *, void *)> RJ_proposal,/**< std::function for the log_likelihood function -- takes double *position, int *param_status,int dimension, int chain_id*/
	void **user_parameters,/**< Void pointer to any parameters the user may need inside log_prior, log_likelihood, or fisher. Should have one pointer for each chain. If this isn't needed, use (void**) NULL**/
	int numThreads, /**< Number of threads to use (=1 is single threaded)*/
	bool pool, /**< boolean to use stochastic chain swapping (MUST have >2 threads)*/
	bool show_prog, /**< boolean whether to print out progress (for example, should be set to ``false'' if submitting to a cluster)*/
	bool update_RJ_width,
	std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
	std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output*/
	std::string likelihood_log_filename,/**< Filename to write the log_likelihood and log_prior at each step -- use empty string to skip*/
	std::string checkpoint_file/**< Filename to output data for checkpoint, if empty string, not saved*/
	)
{
	int max_dimension,min_dimension;
	dimension_from_checkpoint_file(checkpoint_file_start, &min_dimension,&max_dimension);
	int chain_N;
	chain_number_from_checkpoint_file(checkpoint_file_start, &chain_N);
	clock_t start, end, acend;
	double wstart, wend, wacend;
	start = clock();
	wstart = omp_get_wtime();
	bool internal_prog=false;

	//int dynamic_search_length = 2*t0;
	int dynamic_search_length = nu*2;
	//int dynamic_search_length = 200;
	double ***temp_output = allocate_3D_array(chain_N,dynamic_search_length, max_dimension);
	int ***temp_status = allocate_3D_array_int(chain_N,dynamic_search_length, max_dimension);
	int **temp_model_status = NULL;
	temp_model_status = allocate_2D_array_int(chain_N, dynamic_search_length);
	sampler sampler_temp;
	continue_RJPTMCMC_MH_internal(&sampler_temp,checkpoint_file_start,temp_output,
		temp_status,temp_model_status, nested_model_number, dynamic_search_length, 
		swp_freq,log_prior, log_likelihood, fisher, RJ_proposal,user_parameters,
		numThreads, pool, internal_prog, update_RJ_width, statistics_filename, 
		"",  "","",checkpoint_file,true,true);
	
	deallocate_3D_array(temp_output, chain_N, dynamic_search_length, max_dimension);
	deallocate_3D_array(temp_status, chain_N, dynamic_search_length, max_dimension);
	if(temp_model_status){
		deallocate_2D_array(temp_model_status, chain_N,dynamic_search_length);
		temp_model_status = NULL;
	}
	deallocate_sampler_mem(&sampler_temp);

	
	std::cout<<"Starting driver"<<std::endl;
 	RJPTMCMC_MH_dynamic_PT_alloc_comprehensive_internal_driver(sampler_output,
		output,
		parameter_status,
		model_status,
		nested_model_number,
		max_dimension, 	
		min_dimension, 	
		N_steps,	
		chain_N,
		max_chain_N_thermo_ensemble,
		swp_freq,	
		t0,
		nu,
		max_chunk_size,
		chain_distribution_scheme, 
		log_prior,
		log_likelihood,
		fisher,
		RJ_proposal,
		user_parameters,
		numThreads, 
		pool, 
		show_prog, 
		update_RJ_width,
		statistics_filename,
		chain_filename,
		likelihood_log_filename,
		checkpoint_file
		);
	std::cout<<std::endl;
	std::cout<<"WALL time: "<<omp_get_wtime()-wstart<<std::endl;
	return ;

}

/*! \brief Parallel tempered, dynamic chain allocation MCMC. Reversible Jump version, which doesn't currently have any autocorrelation calculations currently.
 *
 * Runs dynamic chain allocation, then it will run the RJPTMCMC routine repeatedly 
 *
 */
void RJPTMCMC_MH_dynamic_PT_alloc_comprehensive_internal(mcmc_sampler_output *sampler_output,
	double **output, /**< [out] Output shape is double[N_steps,dimension]*/
	int **parameter_status, /**< [out] Output shape is int[N_steps,dimension]*/
	int *model_status,
	int nested_model_number,
	int max_dimension, 	/**< maximum dimension of the parameter space being explored*/
	int min_dimension, 	/**< minimum dimension of the parameter space being explored*/
	int N_steps,	/**< Number of total steps to be taken, per chain AFTER chain allocation*/
	int chain_N,/**< Maximum number of chains to use */
	int max_chain_N_thermo_ensemble,/**< Maximum number of chains to use in the thermodynamic ensemble (may use less)*/
	double *initial_pos, 	/**<Initial position in parameter space - shape double[max_dimension]*/
	int *initial_status, 	/**<Initial status in parameter space - shape int[max_dimension]*/
	int initial_model_status, 	/**<Initial model status in parameter space - shape int[nested_model_number]*/
	double *seeding_var, 	/**<Variance of the normal distribution used to seed each chain higher than 0 - shape double[dimension]*/
	double **ensemble_initial_pos,
	int **ensemble_initial_status,
	int *ensemble_initial_model_status,
	double *chain_temps, /**<[out] Final chain temperatures used -- should be shape double[chain_N]*/
	int swp_freq,	/**< the frequency with which chains are swapped*/
	int t0,/**< Time constant of the decay of the chain dynamics  (~1000)*/
	int nu,/**< Initial amplitude of the dynamics (~100)*/
	int max_chunk_size,/**<Maximum number of steps to take in a single sampler run*/
	std::string chain_distribution_scheme, /*How to allocate the remaining chains once equilibrium is reached*/
	std::function<double(double*, int*,int,mcmc_data_interface *,void *)> log_prior,/**< std::function for the log_prior function -- takes double *position, int *param_status, int dimension, int chain_id*/
	std::function<double(double*,int*,int,mcmc_data_interface *,void *)> log_likelihood,/**< std::function for the log_likelihood function -- takes double *position, int *param_status,int dimension, int chain_id*/
	std::function<void(double*,int*,int,double**,mcmc_data_interface *,void *)>fisher,/**< std::function for the fisher function -- takes double *position, int *param_status,int dimension, double **output_fisher, int chain_id*/
	std::function<void(double*,double*,int*,int*,int*,int*,double *,mcmc_data_interface *, void *)> RJ_proposal,/**< std::function for the log_likelihood function -- takes double *position, int *param_status,int dimension, int chain_id*/
	void **user_parameters,/**< Void pointer to any parameters the user may need inside log_prior, log_likelihood, or fisher. Should have one pointer for each chain. If this isn't needed, use (void**) NULL**/
	int numThreads, /**< Number of threads to use (=1 is single threaded)*/
	bool pool, /**< boolean to use stochastic chain swapping (MUST have >2 threads)*/
	bool show_prog, /**< boolean whether to print out progress (for example, should be set to ``false'' if submitting to a cluster)*/
	bool update_RJ_width,
	std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
	std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output*/
	std::string likelihood_log_filename,/**< Filename to write the log_likelihood and log_prior at each step -- use empty string to skip*/
	std::string checkpoint_file/**< Filename to output data for checkpoint, if empty string, not saved*/
	)
{
	clock_t start, end, acend;
	double wstart, wend, wacend;
	start = clock();
	wstart = omp_get_wtime();
	bool internal_prog=false;
	assign_ensemble_temps(chain_temps, chain_N,max_chain_N_thermo_ensemble,TMAX);

	//int dynamic_search_length = 2*t0;
	int dynamic_search_length = nu*2;
	int search_swap_freq = 2;
	//int dynamic_search_length = 200;
	double ***temp_output = allocate_3D_array(chain_N,dynamic_search_length, max_dimension);
	int ***temp_status = allocate_3D_array_int(chain_N,dynamic_search_length, max_dimension);
	int **temp_model_status = NULL;
	temp_model_status = allocate_2D_array_int(chain_N, dynamic_search_length);
	//#####################################################################
	RJPTMCMC_MH_internal(temp_output,temp_status, temp_model_status, nested_model_number,max_dimension, min_dimension,
		dynamic_search_length, chain_N,  
		initial_pos,initial_status,initial_model_status, seeding_var,ensemble_initial_pos,ensemble_initial_status,ensemble_initial_model_status, chain_temps, search_swap_freq, 
		 log_prior, log_likelihood,fisher,RJ_proposal,user_parameters,
		numThreads, pool,internal_prog,update_RJ_width,statistics_filename,"","","",checkpoint_file);
	
	deallocate_3D_array(temp_output, chain_N, dynamic_search_length, max_dimension);
	deallocate_3D_array(temp_status, chain_N, dynamic_search_length, max_dimension);
	if(temp_model_status){
		deallocate_2D_array(temp_model_status, chain_N,dynamic_search_length);
		temp_model_status = NULL;
	}



	//#########################################################################
	std::cout<<"Entering search phase"<<std::endl;
	int status = 0;
	load_temps_checkpoint_file(checkpoint_file, chain_temps, chain_N);
	bool cumulative=true;
	bool full_explore=true;
	int coldchains = count_cold_chains(chain_temps, chain_N);

	dynamic_search_length = 2*t0;
	int temp_length = 1*N_steps;
	if( 2*dynamic_search_length>temp_length){
		temp_length = 2*dynamic_search_length;
	}



	temp_output = allocate_3D_array(chain_N,temp_length, max_dimension);
	temp_status = allocate_3D_array_int(chain_N,temp_length, max_dimension);
	temp_model_status = allocate_2D_array_int(chain_N, temp_length);
	int dynamic_ct = 0 ;
	int dynamic_temp_freq = 1;
	bool continue_dynamic_search=true;

	//#################################################
	//std::cout<<"Annealing"<<std::endl;
	//sampler sampler_ann;
	//continue_PTMCMC_MH_simulated_annealing_internal(&sampler_ann,checkpoint_file,temp_output, dynamic_search_length, 
	//	100,swp_freq,log_prior, log_likelihood, fisher, user_parameters,
	//	numThreads, pool, internal_prog, statistics_filename, 
	//	"", checkpoint_file);


	//deallocate_sampler_mem(&sampler_ann);
	//double temp_factor = TMAX;
	////int temp_factor = TMAX;
	//for(int i = 0 ; i<1; i++){
	//	if(i >=1 ){ temp_factor = 100;}
	//	std::cout<<"Annealing"<<std::endl;
	//	sampler sampler_ann;
	//	continue_RJPTMCMC_MH_simulated_annealing_internal(&sampler_ann,checkpoint_file,temp_output, temp_status, temp_model_status,nested_model_number,dynamic_search_length, 
	//		temp_factor,swp_freq,t0,nu,log_prior, log_likelihood, fisher,RJ_proposal, user_parameters,
	//		numThreads, pool, internal_prog, update_RJ_width,statistics_filename, 
	//		"", checkpoint_file);
	//	//debugger_print(__FILE__,__LINE__,"Writing out swap partners");
	//	//write_file("data/swap_partners_"+std::to_string(i)+".csv",sampler_ann.swap_partners,sampler_ann.chain_N,sampler_ann.chain_N);


	//	deallocate_sampler_mem(&sampler_ann);
	//}
	//#################################################


	while(continue_dynamic_search && dynamic_ct<2){
		if(dynamic_ct >0){search_swap_freq = swp_freq/2;}
		//#############################################
		if(dynamic_ct%dynamic_temp_freq ==0){
			std::cout<<"Temperature Relaxation"<<std::endl;
			continue_RJPTMCMC_MH_dynamic_PT_alloc_full_ensemble_internal(
				checkpoint_file,temp_output, temp_status,temp_model_status, nested_model_number,
				dynamic_search_length, chain_temps, search_swap_freq, t0, nu,
				log_prior, log_likelihood,fisher,RJ_proposal,
				user_parameters,numThreads, pool,internal_prog,update_RJ_width,"","",checkpoint_file,true);

		}

			

		sampler sampler_temp;
		std::cout<<"Exploration"<<std::endl;
		continue_RJPTMCMC_MH_internal(&sampler_temp,checkpoint_file,temp_output,
			temp_status,temp_model_status, nested_model_number, dynamic_search_length, 
			search_swap_freq,log_prior, log_likelihood, fisher, RJ_proposal,user_parameters,
			numThreads, pool, internal_prog, update_RJ_width, statistics_filename, 
			"",  "","",checkpoint_file,true,true);
		int nan_counter = 0;
		for(int i = 0  ;i<sampler_temp.chain_N;i++){
			nan_counter+= sampler_temp.nan_counter[i];
		}
		std::cout<<"NANs in Fisher: "<<nan_counter<<std::endl;
		//##############################################################
		//Reset positions and rewrite checkpoint file -- Decided against this..
		//##############################################################

		//const gsl_rng_type *T_reset;
		//gsl_rng * r_reset;
		//gsl_rng_env_setup();
		//T_reset = gsl_rng_default;
		//r_reset = gsl_rng_alloc(T_reset);
		//double tempT[chain_N];
		//load_temps_checkpoint_file(checkpoint_file, tempT, chain_N);
		//
		//sampler_temp.chain_temps = tempT;
		//
		//double mean, sigma,posterior_norm;
		//int elements;
		//int ensemble_members=sampler_temp.chain_N;//number of chains IN ensemble
		//for(int i = 1 ; i<sampler_temp.chain_N; i++){
		//	if(fabs(sampler_temp.chain_temps[i]-sampler_temp.chain_temps[0]) 
		//		< DOUBLE_COMP_THRESH){
		//		ensemble_members = i;
		//		break;
		//	}
		//}
		////Loop through ensembles
		//for(int i = 0 ; i<  ensemble_members; i++){
		//	//##############################################
		//	//top percentile method
		//	//##############################################
		//	elements = 0 ;
		//	std::vector<std::pair<double,int>> temp_arr;
		//	std::vector<double*> temp_arr_pos;
		//	std::vector<int*> temp_arr_status;
		//	for( int j = 0 ; j<sampler_temp.chain_N/ensemble_members; j++){
		//		int chain_id = i+j*ensemble_members;
		//		for(int k = 0 ; k<=sampler_temp.chain_pos[chain_id]; k++){
		//			std::pair<double,int>P = std::make_pair((sampler_temp.ll_lp_output[chain_id][k][0]/sampler_temp.chain_temps[chain_id]+sampler_temp.ll_lp_output[chain_id][k][1]), k+elements);
		//			temp_arr.push_back(P);
		//			temp_arr_pos.push_back(sampler_temp.output[chain_id][k]);
		//			temp_arr_status.push_back(sampler_temp.param_status[chain_id][k]);
		//		}
		//		elements+=sampler_temp.chain_pos[chain_id];
		//	}

		//	int selection_length = elements;
		//	if(elements > 100){
		//		selection_length = elements*.1;
		//	}
		//	//int selection_length = sampler_temp.chain_N/ensemble_members * 100;
		//	//if(selection_length > elements){selection_length = elements;}
		//	
		//	std::sort(temp_arr.begin(), temp_arr.end());
		//	for(int j = 0 ; j<sampler_temp.chain_N/ensemble_members;j++){
		//		int chain_id = i+j*ensemble_members;
		//		int pos = sampler_temp.chain_pos[chain_id];
		//		int alpha = (int)(gsl_rng_uniform(r_reset)*selection_length);
		//		int temp_id = std::get<1>(temp_arr.at(temp_arr.size()-1 - alpha ));
		//		double posterior = std::get<0>(temp_arr.at(temp_arr.size()-1 - alpha ));
		//		double *temp_pos = temp_arr_pos.at(temp_id);
		//		int *temp_status = temp_arr_status.at(temp_id);
		//		for(int l = 0 ; l<sampler_temp.max_dim; l++){
		//			sampler_temp.output[chain_id][pos][l] = temp_pos[l];
		//			sampler_temp.param_status[chain_id][pos][l] = temp_status[l];
		//		}
		//	
		//	}
		//}
		//for(int i = 0 ; i<sampler_temp.chain_N; i++){
		//	sampler_temp.de_primed[i] = false;
		//}
		//
		//write_checkpoint_file(&sampler_temp, checkpoint_file);

		//gsl_rng_free(r_reset);
		//##############################################################
		//##############################################################

		deallocate_sampler_mem(&sampler_temp);


	
		dynamic_ct++;
	}
	deallocate_3D_array(temp_output,chain_N,temp_length, max_dimension);
	deallocate_3D_array(temp_status,chain_N,temp_length, max_dimension);
	if(temp_model_status){
		deallocate_2D_array(temp_model_status, chain_N,temp_length);
		temp_model_status = NULL;
	}
	std::cout<<"Number of search iterations: "<<dynamic_ct<<std::endl;
	//#########################################################################
	

	std::cout<<"Starting driver"<<std::endl;
 	RJPTMCMC_MH_dynamic_PT_alloc_comprehensive_internal_driver(sampler_output,
		output,
		parameter_status,
		model_status,
		nested_model_number,
		max_dimension, 	
		min_dimension, 	
		N_steps,	
		chain_N,
		max_chain_N_thermo_ensemble,
		swp_freq,	
		t0,
		nu,
		max_chunk_size,
		chain_distribution_scheme, 
		log_prior,
		log_likelihood,
		fisher,
		RJ_proposal,
		user_parameters,
		numThreads, 
		pool, 
		show_prog, 
		update_RJ_width,
		statistics_filename,
		chain_filename,
		likelihood_log_filename,
		checkpoint_file
		);
	std::cout<<std::endl;
	std::cout<<"WALL time: "<<omp_get_wtime()-wstart<<std::endl;
	return ;
}
/*! \brief Parallel tempered, dynamic chain allocation MCMC with output samples with specified maximum autocorrelation. 
 *
 * Runs dynamic chain allocation until the autocorrelation lengths stabilize, then it will run the PTMCMC routine repeatedly, periodically thinning the chain to ensure the auto-correlation length is below the threshold
 *
 * Note: smallest batch size is .5*N_steps, so the target correlation length should be << less than this number
 *
 * Note: This routine only works if the requested number of samples is larger than the typical correlation length of the unthinned chains. This is because the chains are run and tested in batches the size of the requested sample number.
 *
 * Note: This method does NOT guarantee the final autocorrelation length of the chains will the be the target. It merely uses the requested autocorrelation length as a guide to thin the chains as samples are accrued and to estimate the total number of effective samples. Its best to request extra samples and thin the chains out at the end one final time, or to run multiple runs and combine the results at the end.
 */
void continue_PTMCMC_MH_dynamic_PT_alloc_uncorrelated_internal(std::string checkpoint_file_start, 
	mcmc_sampler_output *sampler_output,
	double **output, /**< [out] Output shape is double[N_steps,dimension]*/
	int N_steps,	/**< Number of total steps to be taken, per chain AFTER chain allocation*/
	int max_chain_N_thermo_ensemble,/**< Maximum number of chains to use in the thermodynamic ensemble (may use less)*/
	double *chain_temps, /**<[out] Final chain temperatures used -- should be shape double[chain_N]*/
	int swp_freq,	/**< the frequency with which chains are swapped*/
	int t0,/**< Time constant of the decay of the chain dynamics  (~1000)*/
	int nu,/**< Initial amplitude of the dynamics (~100)*/
	int max_chunk_size,/**<Maximum number of steps to take in a single sampler run*/
	std::string chain_distribution_scheme, /*How to allocate the remaining chains once equilibrium is reached*/
	std::function<double(double*,int* ,int,mcmc_data_interface *,void *)> log_prior,/**< std::function for the log_prior function -- takes double *position, int dimension, int chain_id*/
	std::function<double(double*,int*,int,mcmc_data_interface *,void *)> log_likelihood,/**< std::function for the log_likelihood function -- takes double *position, int dimension, int chain_id*/
	std::function<void(double*,int*,int,double**,mcmc_data_interface *,void *)>fisher,/**< std::function for the fisher function -- takes double *position, int dimension, double **output_fisher, int chain_id*/
	void **user_parameters,/**< Void pointer to any parameters the user may need inside log_prior, log_likelihood, or fisher. Should have one pointer for each chain. If this isn't needed, use (void**) NULL**/
	int numThreads, /**< Number of threads to use (=1 is single threaded)*/
	bool pool, /**< boolean to use stochastic chain swapping (MUST have >2 threads)*/
	bool show_prog, /**< boolean whether to print out progress (for example, should be set to ``false'' if submitting to a cluster)*/
	std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
	std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output*/
	std::string likelihood_log_filename,/**< Filename to write the log_likelihood and log_prior at each step -- use empty string to skip*/
	std::string checkpoint_file/**< Filename to output data for checkpoint, if empty string, not saved*/
	)
{
	int dimension;
	dimension_from_checkpoint_file(checkpoint_file_start, &dimension,&dimension);
	int chain_N;
	chain_number_from_checkpoint_file(checkpoint_file_start, &chain_N);
	clock_t start, end, acend;
	double wstart, wend, wacend;
	start = clock();
	wstart = omp_get_wtime();
	bool internal_prog=false;

	int dynamic_search_length = nu*2;
	double ***temp_output = allocate_3D_array(chain_N,dynamic_search_length, dimension);
	//#####################################################################
	sampler sampler_temp;
	continue_PTMCMC_MH_internal(&sampler_temp,checkpoint_file_start,temp_output, dynamic_search_length, 
		swp_freq,log_prior, log_likelihood, fisher, user_parameters,
		numThreads, pool, internal_prog, statistics_filename, 
		"",  checkpoint_file,true,true);
	deallocate_sampler_mem(&sampler_temp);
	//continue_PTMCMC_MH_internal(checkpoint_file_start, temp_output,  
	//	dynamic_search_length,  max_chain_N_thermo_ensemble, 
	//	 chain_temps, swp_freq, t0, nu,
	//	chain_distribution_scheme, log_prior, log_likelihood,fisher,user_parameters,
	//	numThreads, pool,internal_prog,true,"","",checkpoint_file,true);
	deallocate_3D_array(temp_output, chain_N, dynamic_search_length, dimension);
	
 	PTMCMC_MH_dynamic_PT_alloc_uncorrelated_internal_driver(sampler_output,
		output,
		dimension, 	
		N_steps,	
		chain_N,
		max_chain_N_thermo_ensemble,
		swp_freq,	
		t0,
		nu,
		max_chunk_size,
		chain_distribution_scheme, 
		log_prior,
		log_likelihood,
		fisher,
		user_parameters,
		numThreads, 
		pool, 
		show_prog,
		statistics_filename,
		chain_filename,
		likelihood_log_filename,
		checkpoint_file,
		false
	);
	std::cout<<"WALL time: "<<omp_get_wtime()-wstart<<std::endl;
	return ;
}
/*! \brief Parallel tempered, dynamic chain allocation MCMC with output samples with specified maximum autocorrelation. 
 *
 * Runs dynamic chain allocation until the autocorrelation lengths stabilize, then it will run the PTMCMC routine repeatedly, periodically thinning the chain to ensure the auto-correlation length is below the threshold
 *
 * Note: smallest batch size is .5*N_steps, so the target correlation length should be << less than this number
 *
 * Note: This routine only works if the requested number of samples is larger than the typical correlation length of the unthinned chains. This is because the chains are run and tested in batches the size of the requested sample number.
 *
 * Note: This method does NOT guarantee the final autocorrelation length of the chains will the be the target. It merely uses the requested autocorrelation length as a guide to thin the chains as samples are accrued and to estimate the total number of effective samples. Its best to request extra samples and thin the chains out at the end one final time, or to run multiple runs and combine the results at the end.
 */
void PTMCMC_MH_dynamic_PT_alloc_uncorrelated_internal(mcmc_sampler_output *sampler_output,
	double **output, /**< [out] Output shape is double[N_steps,dimension]*/
	int dimension, 	/**< dimension of the parameter space being explored*/
	int N_steps,	/**< Number of total steps to be taken, per chain AFTER chain allocation*/
	int chain_N,/**< Maximum number of chains to use */
	int max_chain_N_thermo_ensemble,/**< Maximum number of chains to use in the thermodynamic ensemble (may use less)*/
	double *initial_pos, 	/**<Initial position in parameter space - shape double[dimension]*/
	double *seeding_var, 	/**<Variance of the normal distribution used to seed each chain higher than 0 - shape double[dimension]*/
	double **ensemble_initial_pos,
	double *chain_temps, /**<[out] Final chain temperatures used -- should be shape double[chain_N]*/
	int swp_freq,	/**< the frequency with which chains are swapped*/
	int t0,/**< Time constant of the decay of the chain dynamics  (~1000)*/
	int nu,/**< Initial amplitude of the dynamics (~100)*/
	int max_chunk_size,/**<Maximum number of steps to take in a single sampler run*/
	std::string chain_distribution_scheme, /*How to allocate the remaining chains once equilibrium is reached*/
	std::function<double(double*,int* ,int,mcmc_data_interface *,void *)> log_prior,/**< std::function for the log_prior function -- takes double *position, int dimension, int chain_id*/
	std::function<double(double*,int*,int,mcmc_data_interface *,void *)> log_likelihood,/**< std::function for the log_likelihood function -- takes double *position, int dimension, int chain_id*/
	std::function<void(double*,int*,int,double**,mcmc_data_interface *,void *)>fisher,/**< std::function for the fisher function -- takes double *position, int dimension, double **output_fisher, int chain_id*/
	void **user_parameters,/**< Void pointer to any parameters the user may need inside log_prior, log_likelihood, or fisher. Should have one pointer for each chain. If this isn't needed, use (void**) NULL**/
	int numThreads, /**< Number of threads to use (=1 is single threaded)*/
	bool pool, /**< boolean to use stochastic chain swapping (MUST have >2 threads)*/
	bool show_prog, /**< boolean whether to print out progress (for example, should be set to ``false'' if submitting to a cluster)*/
	std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
	std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output*/
	std::string likelihood_log_filename,/**< Filename to write the log_likelihood and log_prior at each step -- use empty string to skip*/
	std::string checkpoint_file/**< Filename to output data for checkpoint, if empty string, not saved*/
	)
{
	clock_t start, end, acend;
	double wstart, wend, wacend;
	start = clock();
	wstart = omp_get_wtime();
	bool internal_prog=false;

	int search_swap_freq = 2;
	//int dynamic_search_length = 2*t0;
	int dynamic_search_length = nu*2;
	//int dynamic_search_length = 200;
	double ***temp_output = allocate_3D_array(chain_N,dynamic_search_length, dimension);
	//#####################################################################
	//PTMCMC_MH_dynamic_PT_alloc_internal(temp_output, dimension, 
	//	dynamic_search_length, chain_N, max_chain_N_thermo_ensemble, 
	//	initial_pos, seeding_var,ensemble_initial_pos, chain_temps, search_swap_freq, t0, nu,
	//	chain_distribution_scheme, log_prior, log_likelihood,fisher,user_parameters,
	//	numThreads, pool,internal_prog,true,"","",checkpoint_file,true);
	//
	//deallocate_3D_array(temp_output, chain_N, dynamic_search_length, dimension);

	if(chain_distribution_scheme == "double"){
		assign_ensemble_temps(chain_temps, chain_N, max_chain_N_thermo_ensemble,TMAX);	
	}
	else{
		debugger_print(__FILE__,__LINE__,"Ooops, only double is allowed for chain_distribution_scheme");
	}
	
	PTMCMC_MH_internal(temp_output, dimension, 
		dynamic_search_length, chain_N,  
		initial_pos, seeding_var,ensemble_initial_pos, chain_temps, search_swap_freq, 
		log_prior, log_likelihood,fisher,user_parameters,
		numThreads, pool,internal_prog,"","",checkpoint_file,true,true);
	
	deallocate_3D_array(temp_output, chain_N, dynamic_search_length, dimension);



	//#########################################################################
	std::cout<<"Entering search phase"<<std::endl;
	int status = 0;
	load_temps_checkpoint_file(checkpoint_file, chain_temps, chain_N);
	bool cumulative=true;
	bool full_explore=true;
	int coldchains = count_cold_chains(chain_temps, chain_N);

	double ave_ac;
	int corr_segments =2;
	int corr_target_ac =.01;
	double corr_converge_thresh = .2;
	int **temp_ac = allocate_2D_array_int(dimension, corr_segments);
	
	dynamic_search_length = 2*t0;
	int temp_length = 1*N_steps;
	if( 2*dynamic_search_length>temp_length){
		temp_length = 2*dynamic_search_length;
	}



	temp_output = allocate_3D_array(chain_N,temp_length, dimension);
	int dynamic_ct = 0 ;
	int dynamic_temp_freq = 1;
	bool continue_dynamic_search=true;
	double max_ac_realloc=0;

	//#################################################
	//Do the annealing process 3 times
	//Start out raising the temperature by a huge amount (x1000)
	//then make it less intense
	//#################################################
	//dynamic_search_length*=2;
	//double temp_factor = TMAX;
	//int temp_factor = TMAX;
	//for(int i = 0 ; i<1; i++){
	//	if(i >=1 ){ temp_factor = 100;}
	//	std::cout<<"Annealing"<<std::endl;
	//	sampler sampler_ann;
	//	continue_PTMCMC_MH_simulated_annealing_internal(&sampler_ann,checkpoint_file,temp_output, dynamic_search_length, 
	//		temp_factor,swp_freq,t0,nu,log_prior, log_likelihood, fisher, user_parameters,
	//		numThreads, pool, internal_prog, statistics_filename, 
	//		"", checkpoint_file);
	//	//#########################################################################
	//	//TESTING

	//	int temp_pos[sampler_ann.chain_N];
	//	for(int z = 0; z<sampler_ann.chain_N;z++){
	//		temp_pos[z] = sampler_ann.chain_pos[z];
	//		if(sampler_ann.restarted_chain[z]){
	//			temp_pos[z] = sampler_ann.N_steps;
	//		}
	//	}
	//	load_temps_checkpoint_file(checkpoint_file, chain_temps, chain_N);
	//	sampler_output->populate_chain_temperatures(chain_temps);

	//	sampler_output->populate_initial_output(temp_output, (int ***)NULL,(int **)NULL,sampler_ann.ll_lp_output,temp_pos)	;
	//	sampler_output->set_trim(0);	
	//	sampler_output->update_cold_chain_list();	
	//	sampler_output->calc_ac_vals(true);
	//	sampler_output->count_indep_samples(true);
	//	sampler_output->create_data_dump(true,false, "data/annealing_output.hdf5");
	//	//#########################################################################

	//	deallocate_sampler_mem(&sampler_ann);
	//}


	//#################################################
	//dynamic_search_length/=2;


	while(continue_dynamic_search && dynamic_ct<2){
		if(dynamic_ct >0){search_swap_freq = swp_freq/2;}
		//std::cout<<"Annealing"<<std::endl;
		//sampler sampler_ann;
		//continue_PTMCMC_MH_simulated_annealing_internal(&sampler_ann,checkpoint_file,temp_output, dynamic_search_length, 
		//	10,swp_freq,t0,nu,log_prior, log_likelihood, fisher, user_parameters,
		//	numThreads, pool, internal_prog, statistics_filename, 
		//	"", checkpoint_file);

		//deallocate_sampler_mem(&sampler_ann);
		//#############################################
		//#############################################
		if(dynamic_ct%dynamic_temp_freq ==0){
			std::cout<<"Temperature Relaxation"<<std::endl;
			continue_PTMCMC_MH_dynamic_PT_alloc_full_ensemble_internal(checkpoint_file,temp_output, 
				dynamic_search_length, chain_temps, search_swap_freq, t0, nu,
				 log_prior, log_likelihood,fisher,
				user_parameters,numThreads, pool,internal_prog,"","",checkpoint_file,true);
		}

		//#######################################
		sampler sampler_temp;
		std::cout<<"Exploration"<<std::endl;
		continue_PTMCMC_MH_internal(&sampler_temp,checkpoint_file,temp_output, dynamic_search_length, 
			search_swap_freq,log_prior, log_likelihood, fisher, user_parameters,
			numThreads, pool, internal_prog, statistics_filename, 
			"",  checkpoint_file,true,true);
		//debugger_print(__FILE__,__LINE__,"Writing out swap partners");
		//write_file("data/swap_partners_"+std::to_string(dynamic_ct)+".csv",sampler_temp.swap_partners,sampler_temp.chain_N,sampler_temp.chain_N);
		//write_file("data/swap_accepts_"+std::to_string(dynamic_ct)+".csv",sampler_temp.swap_accepts,sampler_temp.chain_N,sampler_temp.chain_N);

		//##############################################################
		//Reset positions and rewrite checkpoint file -- turning this off for now
		//##############################################################

		//const gsl_rng_type *T_reset;
		//gsl_rng * r_reset;
		//gsl_rng_env_setup();
		//T_reset = gsl_rng_default;
		//r_reset = gsl_rng_alloc(T_reset);
		//double tempT[chain_N];
		//load_temps_checkpoint_file(checkpoint_file, tempT, chain_N);
		//
		//sampler_temp.chain_temps = tempT;
		//
		//double mean, sigma,posterior_norm;
		//int elements;
		//int ensemble_members=sampler_temp.chain_N;//number of chains IN ensemble
		//for(int i = 1 ; i<sampler_temp.chain_N; i++){
		//	if(fabs(sampler_temp.chain_temps[i]-sampler_temp.chain_temps[0]) 
		//		< DOUBLE_COMP_THRESH){
		//		ensemble_members = i;
		//		break;
		//	}
		//}
		////Loop through ensembles
		//for(int i = 0 ; i<  ensemble_members; i++){
		//	//##############################################
		//	//top percentile method
		//	//##############################################
		//	elements = 0 ;
		//	std::vector<std::pair<double,int>> temp_arr;
		//	std::vector<double*> temp_arr_pos;
		//	for( int j = 0 ; j<sampler_temp.chain_N/ensemble_members; j++){
		//		int chain_id = i+j*ensemble_members;
		//		for(int k = 0 ; k<=sampler_temp.chain_pos[chain_id]; k++){
		//			std::pair<double,int>P = std::make_pair((sampler_temp.ll_lp_output[chain_id][k][0]/sampler_temp.chain_temps[chain_id]+sampler_temp.ll_lp_output[chain_id][k][1]), k+elements);
		//			temp_arr.push_back(P);
		//			temp_arr_pos.push_back(sampler_temp.output[chain_id][k]);
		//		}
		//		elements+=sampler_temp.chain_pos[chain_id];
		//	}

		//	int selection_length = elements*.1;
		//	//int selection_length = sampler_temp.chain_N/ensemble_members * 100;
		//	//if(selection_length > elements){selection_length = elements;}
		//	
		//	std::sort(temp_arr.begin(), temp_arr.end());
		//	for(int j = 0 ; j<sampler_temp.chain_N/ensemble_members;j++){
		//		int chain_id = i+j*ensemble_members;
		//		int pos = sampler_temp.chain_pos[chain_id];
		//		int alpha = (int)(gsl_rng_uniform(r_reset)*selection_length);
		//		int temp_id = std::get<1>(temp_arr.at(temp_arr.size()-1 - alpha ));
		//		double posterior = std::get<0>(temp_arr.at(temp_arr.size()-1 - alpha ));
		//		double *temp_pos = temp_arr_pos.at(temp_id);
		//		for(int l = 0 ; l<sampler_temp.dimension; l++){
		//			sampler_temp.output[chain_id][pos][l] = temp_pos[l];
		//		}
		//	
		//	}
		//	//##############################################
		//	//Averaging method -- not promising
		//	//##############################################
		//	//Loop through dimensions
		//	//for(int l = 0 ; l < sampler_temp.dimension; l++){
		//	//	mean = 0 ;
		//	//	sigma = 0 ;
		//	//	elements = 0 ;
		//	//	posterior_norm = 0 ;
		//	//	//Loop through members of ensemble
		//	//	for( int j = 0 ; j<sampler_temp.chain_N/ensemble_members; j++){
		//	//		int chain_id = i+j*ensemble_members;
		//	//		for(int k = 0 ; k<sampler_temp.chain_pos[chain_id]; k++){
		//	//			posterior_norm +=(sampler_temp.ll_lp_output[chain_id][k][0]+sampler_temp.ll_lp_output[chain_id][k][1]	);
		//	//		}
		//	//		elements+=sampler_temp.chain_pos[chain_id];
		//	//	}
		//	//	posterior_norm /= elements;
		//	//	posterior_norm=0;
		//	//	std::cout<<posterior_norm<<std::endl;
		//	//	for( int j = 0 ; j<sampler_temp.chain_N/ensemble_members; j++){
		//	//		int chain_id = i+j*ensemble_members;
		//	//		for(int k = 0 ; k<sampler_temp.chain_pos[chain_id]; k++){
		//	//			mean += exp((sampler_temp.ll_lp_output[chain_id][k][0]+sampler_temp.ll_lp_output[chain_id][k][1])-posterior_norm) * (sampler_temp.output[chain_id][k][l]);
		//	//		}
		//	//	}
		//	//	mean/=elements;
		//	//	for( int j = 0 ; j<sampler_temp.chain_N/ensemble_members; j++){
		//	//		int chain_id = i+j*ensemble_members;
		//	//		for(int k = 0 ; k<sampler_temp.chain_pos[chain_id]; k++){
		//	//			sigma += pow_int(exp((sampler_temp.ll_lp_output[chain_id][k][0]+sampler_temp.ll_lp_output[chain_id][k][1])-posterior_norm) * sampler_temp.output[chain_id][k][l]-mean,2);
		//	//		}
		//	//	}
		//	//	sigma/=elements;
		//	//	sigma = std::sqrt(sigma);
		//	//	debugger_print(__FILE__,__LINE__,"Ensemble number");
		//	//	debugger_print(__FILE__,__LINE__,i);
		//	//	debugger_print(__FILE__,__LINE__,"Dimension");
		//	//	debugger_print(__FILE__,__LINE__,l);
		//	//	debugger_print(__FILE__,__LINE__,"mean");
		//	//	debugger_print(__FILE__,__LINE__,mean);
		//	//	debugger_print(__FILE__,__LINE__,"STD");
		//	//	debugger_print(__FILE__,__LINE__,sigma);
		//	//	for(int j = 0 ; j<sampler_temp.chain_N/ensemble_members;j++){
		//	//		int chain_id = i+j*ensemble_members;
		//	//		int pos = sampler_temp.chain_pos[chain_id];
		//	//		double alpha = gsl_ran_gaussian(r_reset,sigma)+mean;
		//	//		debugger_print(__FILE__,__LINE__,"Pos");
		//	//		debugger_print(__FILE__,__LINE__,alpha);
		//	//		sampler_temp.output[chain_id][pos][l] = alpha;
		//	//	
		//	//	}
		//	//}
		//	//##############################################
		//	//##############################################
		//}
		//for(int i = 0 ; i<sampler_temp.chain_N; i++){
		//	sampler_temp.de_primed[i] = false;
		//}
		//
		//write_checkpoint_file(&sampler_temp, checkpoint_file);

		//gsl_rng_free(r_reset);
		//##############################################################
		//##############################################################

		deallocate_sampler_mem(&sampler_temp);


			
		//####################################################################################
		//####################################################################################
		/* Do analysis one chain at a time, then combine*/
		if(!full_explore){
			coldchains = count_cold_chains(chain_temps, chain_N);
			int ***full_temp_ac = allocate_3D_array_int(coldchains,dimension, corr_segments);
			double ***full_temp_output = allocate_3D_array(coldchains, dynamic_search_length,dimension);
			int ccct=0;
			for(int i = 0 ; i<chain_N;i++){
				if(fabs(chain_temps[i]-1)<DOUBLE_COMP_THRESH){
					for(int k = 0 ; k<dynamic_search_length; k++){
						for(int j = 0 ; j<dimension; j++){
							full_temp_output[ccct][k][j]=temp_output[i][k][j];
						}
					}
					ccct++;
				}
			}
			auto_corr_from_data_batch(full_temp_output, dynamic_search_length, dimension, coldchains,full_temp_ac, corr_segments, corr_target_ac, numThreads, cumulative);
			continue_dynamic_search=false;
			double ave_max_ac=0;
			max_ac_realloc=0;
			ccct=0;
			for(int i = 0 ; i<chain_N;i++){
				if(fabs(chain_temps[i]-1)<DOUBLE_COMP_THRESH){
					//auto_corr_from_data(temp_output[i], dynamic_search_length, dimension, temp_ac, corr_segments, corr_target_ac, numThreads, cumulative);
					max_ac_realloc=0;
					for(int k =0 ; k<dimension; k++){
						if(full_temp_ac[ccct][k][corr_segments-1] >max_ac_realloc){
							max_ac_realloc=full_temp_ac[ccct][k][corr_segments-1];
						}	
					}
					ave_max_ac+=max_ac_realloc;
					for(int k = 0 ; k<dimension; k++){
						double variance=0;
						variance_list(full_temp_ac[ccct][k],corr_segments,&variance);
						double mean=0;
						mean_list(full_temp_ac[ccct][k],corr_segments,&mean);
						double cv = sqrt(variance)/mean;
						
						if(cv >corr_converge_thresh){
							std::cout<<"FAILED "<<cv<<" "<<mean<<" "<<variance<<" "<<k<<std::endl;
							continue_dynamic_search=true;
						}
						//if(continue_dynamic_search){break;}
					}
					ccct++;
				}
			}

			max_ac_realloc = ave_max_ac/coldchains;
			deallocate_3D_array(full_temp_ac,coldchains,dimension, corr_segments);
			deallocate_3D_array(full_temp_output,coldchains, dynamic_search_length,dimension);
		}
	
		dynamic_ct++;
	}
	deallocate_2D_array(temp_ac, dimension, corr_segments);
	deallocate_3D_array(temp_output,chain_N,temp_length, dimension);
	std::cout<<"Number of search iterations: "<<dynamic_ct<<std::endl;
	//#########################################################################
	

	std::cout<<"Starting driver"<<std::endl;
 	PTMCMC_MH_dynamic_PT_alloc_uncorrelated_internal_driver(sampler_output,
		output,
		dimension, 	
		N_steps,	
		chain_N,
		max_chain_N_thermo_ensemble,
		swp_freq,	
		t0,
		nu,
		max_chunk_size,
		chain_distribution_scheme, 
		log_prior,
		log_likelihood,
		fisher,
		user_parameters,
		numThreads, 
		pool, 
		show_prog, 
		statistics_filename,
		chain_filename,
		likelihood_log_filename,
		checkpoint_file,
		true
		);
	std::cout<<std::endl;
	std::cout<<"WALL time: "<<omp_get_wtime()-wstart<<std::endl;
	return ;
}

/*! \brief Driver routine for the uncorrelated sampler -- trying not to repeat code
 * 
 * Assumed that the checkpoint_file has already been populated --
 *
 * It will overwrite all the file paths
 *
 */
void RJPTMCMC_MH_dynamic_PT_alloc_comprehensive_internal_driver(mcmc_sampler_output *sampler_output,
	double **output,
	int **status,
	int *model_status,
	int nested_model_number,
	int max_dimension, 	/**< dimension of the parameter space being explored*/
	int min_dimension, 	/**< dimension of the parameter space being explored*/
	int N_steps,	/**< Number of total steps to be taken, per chain AFTER chain allocation*/
	int chain_N,/**< Maximum number of chains to use */
	int max_chain_N_thermo_ensemble,/**< Maximum number of chains to use in the thermodynamic ensemble (may use less)*/
	int swp_freq,	/**< the frequency with which chains are swapped*/
	int t0,/**< Time constant of the decay of the chain dynamics  (~1000)*/
	int nu,/**< Initial amplitude of the dynamics (~100)*/
	int max_chunk_size,/**<Maximum number of steps to take in a single sampler run*/
	std::string chain_distribution_scheme, /*How to allocate the remaining chains once equilibrium is reached*/
	std::function<double(double*, int*,int,mcmc_data_interface *,void *)> log_prior,/**< std::function for the log_prior function -- takes double *position, int *param_status, int dimension, int chain_id*/
	std::function<double(double*,int*,int,mcmc_data_interface *,void *)> log_likelihood,/**< std::function for the log_likelihood function -- takes double *position, int *param_status,int dimension, int chain_id*/
	std::function<void(double*,int*,int,double**,mcmc_data_interface *,void *)>fisher,/**< std::function for the fisher function -- takes double *position, int *param_status,int dimension, double **output_fisher, int chain_id*/
	std::function<void(double*,double*,int*,int*,int*,int*,double *,mcmc_data_interface *, void *)> RJ_proposal,/**< std::function for the log_likelihood function -- takes double *position, int *param_status,int dimension, int chain_id*/
	void **user_parameters,/**< Void pointer to any parameters the user may need inside log_prior, log_likelihood, or fisher. Should have one pointer for each chain. If this isn't needed, use (void**) NULL**/
	int numThreads, /**< Number of threads to use (=1 is single threaded)*/
	bool pool, /**< boolean to use stochastic chain swapping (MUST have >2 threads)*/
	bool show_prog, /**< boolean whether to print out progress (for example, should be set to ``false'' if submitting to a cluster)*/
	bool update_RJ_width,
	std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
	std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output*/
	std::string likelihood_log_filename,/**< Filename to write the log_likelihood and log_prior at each step -- use empty string to skip*/
	std::string checkpoint_file/**< Filename to output data for checkpoint, if empty string, not saved*/
	)
{
	bool cold_chain_only_writeout = true;
	int step_status = 0;
	double chain_temps[chain_N];
	load_temps_checkpoint_file(checkpoint_file, chain_temps, chain_N);
	bool cumulative=true;
	bool internal_prog=false;
	int coldchains = count_cold_chains(chain_temps, chain_N);
	
	int dynamic_search_length = 2*t0;
	int temp_length = max_chunk_size;
	if( dynamic_search_length>temp_length){
		temp_length = dynamic_search_length;
	}



	double ***temp_output = allocate_3D_array(chain_N,temp_length, max_dimension);
	int ***temp_status = allocate_3D_array_int(chain_N,temp_length, max_dimension);
	int **temp_model_status = NULL;
	temp_model_status = allocate_2D_array_int(chain_N, temp_length);
	
	int dynamic_ct = 0 ;
	int dynamic_temp_freq = 1;
	double max_ac_realloc=0;


	coldchains = count_cold_chains(chain_temps, chain_N);
	int realloc_temps_length = 2*N_steps;//Steps before re-allocating chain temps
	int realloc_temps_thresh = realloc_temps_length;
	bool realloc = false;
	bool init = true;
	bool relax = true;
	double ac_save;
	int max_search_iterations = 3;
	int search_iterations_ct = 0;
	while(step_status<N_steps){
		if(realloc || step_status>realloc_temps_thresh){
		//if(false){
			if( t0<temp_length){
				dynamic_search_length = t0;
			}
			else{
				dynamic_search_length = temp_length;
			}
			//continue_PTMCMC_MH_dynamic_PT_alloc_internal(checkpoint_file,temp_output, 
			//	dynamic_search_length,  max_chain_N_thermo_ensemble, 
			//	 chain_temps, swp_freq, t0, nu,
			//	chain_distribution_scheme, log_prior, log_likelihood,fisher,
			//	user_parameters,numThreads, pool,internal_prog,false,"","",checkpoint_file,false);

			std::cout<<"Temperature Reallocation"<<std::endl;
			continue_RJPTMCMC_MH_dynamic_PT_alloc_full_ensemble_internal(
				checkpoint_file,temp_output, temp_status,temp_model_status, nested_model_number,
				dynamic_search_length, chain_temps, swp_freq, t0, nu,
				 log_prior, log_likelihood,fisher, RJ_proposal,
				user_parameters,numThreads, pool,internal_prog,update_RJ_width,"","",checkpoint_file,false);
			sampler sampler_temp;
			std::cout<<"Exploration"<<std::endl;
			continue_RJPTMCMC_MH_internal(
				&sampler_temp,checkpoint_file,temp_output, temp_status,temp_model_status, nested_model_number, 
				(int)(dynamic_search_length), 
				swp_freq,log_prior, log_likelihood, fisher, RJ_proposal, 
				user_parameters,numThreads, pool, internal_prog, update_RJ_width,
				statistics_filename, "",  "","",checkpoint_file,true,false);
			deallocate_sampler_mem(&sampler_temp);


			realloc_temps_thresh+=realloc_temps_length;
			realloc=false;
			//We need to recreate the data_dump_file because the chains may
			//have changed numbers in ensembles
			relax=true;
			init=true;
			step_status= 0;
		}
		double harvest_pool = true;
		sampler sampler;
		continue_RJPTMCMC_MH_internal(
			&sampler, checkpoint_file,temp_output, temp_status,temp_model_status, nested_model_number, temp_length, 
			swp_freq,log_prior, log_likelihood, fisher,RJ_proposal, user_parameters,
			numThreads, harvest_pool, internal_prog, update_RJ_width, statistics_filename, 
			"", "","",checkpoint_file,false,false);
		int nan_counter = 0;
		for(int i = 0  ;i<sampler.chain_N;i++){
			nan_counter+= sampler.nan_counter[i];
		}
		std::cout<<"NANs in Fisher: "<<nan_counter<<std::endl;

		load_temps_checkpoint_file(checkpoint_file, chain_temps, chain_N);
		sampler_output->populate_chain_temperatures(chain_temps);
		if(init){
			debugger_print(__FILE__,__LINE__,"Init structure");
			//###################################################
			//TESTING
			int temp_pos[sampler.chain_N];
			for(int z = 0; z<sampler.chain_N;z++){
				temp_pos[z] = sampler.chain_pos[z];
				if(sampler.restarted_chain[z]){
					temp_pos[z] = sampler.N_steps;
				}
			}
			sampler_output->populate_initial_output(temp_output, temp_status,temp_model_status,sampler.ll_lp_output,temp_pos)	;
			//sampler_output->populate_initial_output(temp_output, temp_status,temp_model_status,sampler.ll_lp_output,sampler.chain_pos)	;
			//###################################################
			sampler_output->set_trim(0);	
			sampler_output->update_cold_chain_list();	
			init=false;
			debugger_print(__FILE__,__LINE__,"Finished init structure");
			debugger_print(__FILE__,__LINE__,"Creating dump");
			sampler_output->create_data_dump(cold_chain_only_writeout,false, chain_filename);
			debugger_print(__FILE__,__LINE__,"Finished Creating dump");
		}
		else{
			debugger_print(__FILE__,__LINE__,"Appending structure");
			//###################################################
			//TESTING
			int temp_pos[sampler.chain_N];
			for(int z = 0; z<sampler.chain_N;z++){
				temp_pos[z] = sampler.chain_pos[z];
				if(sampler.restarted_chain[z]){
					temp_pos[z] = sampler.N_steps;
				}
			}
			sampler_output->append_to_output(temp_output,temp_status,temp_model_status,sampler.ll_lp_output,temp_pos);
			//debugger_print(__FILE__,__LINE__,"Appending FULL data dump");
			//sampler_output->append_to_output(temp_output,temp_status,temp_model_status,sampler.ll_lp_output,sampler.chain_pos);
			//###################################################
			debugger_print(__FILE__,__LINE__,"Finished appending structure");
			debugger_print(__FILE__,__LINE__,"Appending dump");
			sampler_output->append_to_data_dump(chain_filename);
			//debugger_print(__FILE__,__LINE__,"Appending FULL data dump");
			//sampler_output->append_to_data_dump("data/test_full.hdf5");
			debugger_print(__FILE__,__LINE__,"Finished appending dump");
		}
		step_status += temp_length;
		//###########################################################3
		//###########################################################3
		//###########################################################3

		double temp_temps[sampler.chain_N];
		load_temps_checkpoint_file(checkpoint_file, temp_temps, sampler.chain_N);
		double swap_targets[2] = {.1,1};
		int cold_ids[sampler.chain_N];
		int cold_chains_ct=0;
		for(int k = 0 ; k<sampler.chain_N; k++){
			if(fabs(temp_temps[k] - 1.) < DOUBLE_COMP_THRESH){
				cold_ids[cold_chains_ct]=k;
				cold_chains_ct++;
			}	
		}
		int ensemble_members=sampler.chain_N;
		if(cold_chains_ct>1){ensemble_members = cold_ids[1]-cold_ids[0];}
		int ensemble_num = ceil((double)sampler.chain_N/ensemble_members);	
		
		double averages[ensemble_members];
		double average_accept_up[ensemble_members];
		int average_accept_up_attempts[ensemble_members];
		int average_accept_down_attempts[ensemble_members];
		double average_accept_down[ensemble_members];
		double average_attempts[ensemble_members];
		int average_chain_nums[ensemble_members];
		
		for(int i = 0 ; i<ensemble_members; i++){
			averages[i]=0;
			average_attempts[i]=0;
			average_accept_up[i]=0;
			average_accept_down_attempts[i]=0;
			average_accept_up_attempts[i]=0;
			average_accept_down[i]=0;
			average_chain_nums[i]=0;
		}
		for(int k = 0 ; k<sampler.chain_N; k++){
			averages[k%ensemble_members]+=(double)sampler.swap_accept_ct[k] / 
					(sampler.swap_accept_ct[k] + sampler.swap_reject_ct[k]);
			average_attempts[k%ensemble_members]+=(double)(sampler.swap_accept_ct[k] + sampler.swap_reject_ct[k]);
			average_chain_nums[k%ensemble_members]++;
			for(int l = 0 ; l<sampler.chain_N; l++){
				if(temp_temps[l] < temp_temps[k]){
					average_accept_down[k%ensemble_members]+=sampler.swap_accepts[k][l];
					average_accept_down_attempts[k%ensemble_members]+=sampler.swap_partners[k][l];

				}
				else if(temp_temps[l] > temp_temps[k]){
					average_accept_up[k%ensemble_members]+=sampler.swap_accepts[k][l];
					average_accept_up_attempts[k%ensemble_members]+=sampler.swap_partners[k][l];
				}

			}
		}
		
		bool reset = false;
		debugger_print(__FILE__,__LINE__,"Temp, Average [acceptance, to colder chains, to hotter chains], Average attempts");
		for(int i = 0 ; i<ensemble_members; i++){
			if(average_accept_up_attempts[i] != 0){
				average_accept_up[i] /= average_accept_up_attempts[i];
			}
			else{
				average_accept_up[i] =0;
			}
			if(average_accept_down_attempts[i] != 0){
				average_accept_down[i] /= average_accept_down_attempts[i];
			}
			else{
				average_accept_down[i] =0;
			}
			averages[i]/=average_chain_nums[i];
			average_attempts[i]/=average_chain_nums[i];
			//debugger_print(__FILE__,__LINE__,std::to_string(averages[i]) + " "+std::to_string(average_attempts[i]));
			debugger_print(__FILE__,__LINE__,
				std::to_string(temp_temps[i]) + " "+
				std::to_string(averages[i]) + " "+
				std::to_string(average_accept_down[i]) + " "+
				std::to_string(average_accept_up[i]) + " "+
				std::to_string(average_attempts[i])
			);
			if(search_iterations_ct < max_search_iterations){
				if(averages[i]<swap_targets[0] || averages[i]>swap_targets[1] ){
					debugger_print(__FILE__,__LINE__,"Average swap not in range -- restarting");
					reset=true;
				}
				else if ( (average_accept_up[i] < swap_targets[0] || average_accept_up[i] > swap_targets[1] ) && (i!= 0 && i!=(ensemble_members-1))){
					debugger_print(__FILE__,__LINE__,"Average swap to hot chains not in range -- restarting");
					reset=true;
				}
				else if ( (average_accept_down[i] < swap_targets[0] || average_accept_down[i] > swap_targets[1]) && (i!= 0 && i!=(ensemble_members-1))){
					debugger_print(__FILE__,__LINE__,"Average swap to cold chains not in range -- restarting");
					reset=true;
				}
			}
		}
		if (reset){
			if(!realloc){
				search_iterations_ct++;
			}

			realloc=true;
		}
		
		//double temp_temps[sampler.chain_N];
		//load_temps_checkpoint_file(checkpoint_file, temp_temps, sampler.chain_N);
		//double swap_targets[2] = {.1,1};
		//int cold_ids[sampler.chain_N];
		//int cold_chains_ct=0;
		//for(int k = 0 ; k<sampler.chain_N; k++){
		//	if(fabs(temp_temps[k] - 1.) < DOUBLE_COMP_THRESH){
		//		cold_ids[cold_chains_ct]=k;
		//		cold_chains_ct++;
		//	}	
		//}
		//int ensemble_members=sampler.chain_N;
		//if(cold_chains_ct>1){ensemble_members = cold_ids[1]-cold_ids[0];}
		//int ensemble_num = ceil((double)sampler.chain_N/ensemble_members);	
		//
		//double averages[ensemble_members];
		//double average_attempts[ensemble_members];
		//int average_chain_nums[ensemble_members];
		//for(int i = 0 ; i<ensemble_members; i++){
		//	averages[i]=0;
		//	average_attempts[i]=0;
		//	average_chain_nums[i]=0;
		//}
		//for(int k = 0 ; k<sampler.chain_N; k++){
		//	averages[k%ensemble_members]+=(double)sampler.swap_accept_ct[k] / 
		//			(sampler.swap_accept_ct[k] + sampler.swap_reject_ct[k]);
		//	average_attempts[k%ensemble_members]+=(double)(sampler.swap_accept_ct[k] + sampler.swap_reject_ct[k]);
		//	average_chain_nums[k%ensemble_members]++;
		//}
		//debugger_print(__FILE__,__LINE__,"Swap averages:");
		//for(int i = 0 ; i<ensemble_members; i++){
		//	averages[i]/=average_chain_nums[i];
		//	average_attempts[i]/=average_chain_nums[i];
		//	debugger_print(__FILE__,__LINE__,std::to_string(averages[i]) + " "+std::to_string(average_attempts[i]));
		//	if((averages[i]<swap_targets[0] || averages[i]>swap_targets[1]) && search_iterations_ct<max_search_iterations){
		//		if(!realloc){
		//			search_iterations_ct++;
		//		}

		//		realloc=true;
		//	}
		//}
		
		//###########################################################3
		deallocate_sampler_mem(&sampler);
			
		//if(false){
		//if(temp_length < 10*max_ac_realloc){

		//	if(10*max_ac_realloc < max_chunk_size){
		//		deallocate_3D_array(temp_output, chain_N, temp_length, dimension);
		//		temp_length = 10*max_ac_realloc;
		//		temp_output = allocate_3D_array(chain_N,temp_length, dimension);
		//	}
		//	else{
		//		std::cout<<"WARNING -- hit maximum chunk size for a single sampler run"<<std::endl;
		//		std::cout<<"Independent samples per batch are projected to be "<<max_chunk_size/max_ac_realloc<<" and at least 1000 samples per AC calculation is recommended"<<std::endl;
		//		deallocate_3D_array(temp_output, chain_N, temp_length, dimension);
		//		temp_length = max_chunk_size;
		//		temp_output = allocate_3D_array(chain_N,temp_length, dimension);
		//	}
		//}
		//else if(temp_length>5000*max_ac_realloc){
		//	deallocate_3D_array(temp_output, chain_N, temp_length, dimension);
		//	temp_length = 5000*max_ac_realloc;
		//	temp_output = allocate_3D_array(chain_N,temp_length, dimension);	

		//}
		//max_ac_realloc=0;
		printProgress((double)step_status/N_steps);
	}
	
	//Maybe write new stat file format
	

	//Cleanup
	deallocate_3D_array(temp_output, chain_N, temp_length, max_dimension);
	deallocate_3D_array(temp_status, chain_N, temp_length, max_dimension);
	if(temp_model_status){
		deallocate_2D_array(temp_model_status, chain_N, temp_length);
		temp_model_status = NULL;
	}
	return;
}
/*! \brief Driver routine for the uncorrelated sampler -- trying not to repeat code
 * 
 * Assumed that the checkpoint_file has already been populated --
 *
 * It will overwrite all the file paths
 *
 */
void PTMCMC_MH_dynamic_PT_alloc_uncorrelated_internal_driver(mcmc_sampler_output *sampler_output,
	double **output,
	int dimension, 	/**< dimension of the parameter space being explored*/
	int N_steps,	/**< Number of total steps to be taken, per chain AFTER chain allocation*/
	int chain_N,/**< Maximum number of chains to use */
	int max_chain_N_thermo_ensemble,/**< Maximum number of chains to use in the thermodynamic ensemble (may use less)*/
	int swp_freq,	/**< the frequency with which chains are swapped*/
	int t0,/**< Time constant of the decay of the chain dynamics  (~1000)*/
	int nu,/**< Initial amplitude of the dynamics (~100)*/
	int max_chunk_size,/**<Maximum number of steps to take in a single sampler run*/
	std::string chain_distribution_scheme, /*How to allocate the remaining chains once equilibrium is reached*/
	std::function<double(double*,int* ,int,mcmc_data_interface *,void *)> log_prior,/**< std::function for the log_prior function -- takes double *position, int dimension, int chain_id*/
	std::function<double(double*,int*,int,mcmc_data_interface *,void *)> log_likelihood,/**< std::function for the log_likelihood function -- takes double *position, int dimension, int chain_id*/
	std::function<void(double*,int*,int,double**,mcmc_data_interface *,void *)>fisher,/**< std::function for the fisher function -- takes double *position, int dimension, double **output_fisher, int chain_id*/
	void **user_parameters,/**< Void pointer to any parameters the user may need inside log_prior, log_likelihood, or fisher. Should have one pointer for each chain. If this isn't needed, use (void**) NULL**/
	int numThreads, /**< Number of threads to use (=1 is single threaded)*/
	bool pool, /**< boolean to use stochastic chain swapping (MUST have >2 threads)*/
	bool show_prog, /**< boolean whether to print out progress (for example, should be set to ``false'' if submitting to a cluster)*/
	std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
	std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output*/
	std::string likelihood_log_filename,/**< Filename to write the log_likelihood and log_prior at each step -- use empty string to skip*/
	std::string checkpoint_file,/**< Filename to output data for checkpoint, if empty string, not saved*/
	bool continue_burn
	)
{
	bool cold_chain_only_writeout = true;
	bool calc_evidence = true;
	bool local_continue_burn = continue_burn;
	int status = 0;
	double chain_temps[chain_N];
	load_temps_checkpoint_file(checkpoint_file, chain_temps, chain_N);
	bool cumulative=true;
	bool internal_prog=false;
	int coldchains = count_cold_chains(chain_temps, chain_N);
	

	double ave_ac;
	int corr_segments = 2;
	int **temp_ac = allocate_2D_array_int(dimension, corr_segments);
	
	int dynamic_search_length = 2*t0;
	int temp_length = 1*N_steps;
	if( dynamic_search_length>temp_length){
		temp_length = dynamic_search_length;
	}



	double ***temp_output = allocate_3D_array(chain_N,temp_length, dimension);
	int dynamic_ct = 0 ;
	int dynamic_temp_freq = 1;
	double max_ac_realloc=0;


	coldchains = count_cold_chains(chain_temps, chain_N);
	int realloc_temps_length = 2*N_steps;//Steps before re-allocating chain temps
	int realloc_temps_thresh = realloc_temps_length;
	bool realloc = false;
	bool init = true;
	bool relax = true;
	double ac_save;
	int max_search_iterations = 0;
	int search_iterations_ct = 0;
	int spct=0;
	while(status<N_steps){
		if(realloc || status>realloc_temps_thresh){
		//if(false){
			if( t0<temp_length){
				dynamic_search_length = t0;
			}
			else{
				dynamic_search_length = temp_length;
			}
			//continue_PTMCMC_MH_dynamic_PT_alloc_internal(checkpoint_file,temp_output, 
			//	dynamic_search_length,  max_chain_N_thermo_ensemble, 
			//	 chain_temps, swp_freq, t0, nu,
			//	chain_distribution_scheme, log_prior, log_likelihood,fisher,
			//	user_parameters,numThreads, pool,internal_prog,false,"","",checkpoint_file,false);

			continue_PTMCMC_MH_dynamic_PT_alloc_full_ensemble_internal(checkpoint_file,temp_output, 
				dynamic_search_length, chain_temps, swp_freq, t0, nu,
				 log_prior, log_likelihood,fisher,
				user_parameters,numThreads, pool,internal_prog,"","",checkpoint_file,false);
			sampler sampler_temp;
			std::cout<<"Exploration"<<std::endl;
			continue_PTMCMC_MH_internal(&sampler_temp,checkpoint_file,temp_output, (int)(dynamic_search_length), 
				swp_freq,log_prior, log_likelihood, fisher, user_parameters,
				numThreads, pool, internal_prog, statistics_filename, 
				"",  checkpoint_file,true,true);
			deallocate_sampler_mem(&sampler_temp);


			realloc_temps_thresh+=realloc_temps_length;
			realloc=false;
			//We need to recreate the data_dump_file because the chains may
			//have changed numbers in ensembles
			relax=true;
			local_continue_burn = true;
			init=true;
		}
		double harvest_pool = pool;
		sampler sampler;

		continue_PTMCMC_MH_internal(&sampler, checkpoint_file,temp_output, temp_length, 
			swp_freq,log_prior, log_likelihood, fisher, user_parameters,
			numThreads, harvest_pool, internal_prog, statistics_filename, 
			"", checkpoint_file,false,false);
		//debugger_print(__FILE__,__LINE__,"Writing out swap partners");
		//write_file("data/swap_partners_"+std::to_string(spct)+".csv",sampler.swap_partners,sampler.chain_N,sampler.chain_N);
		//write_file("data/swap_accepts_"+std::to_string(spct)+".csv",sampler.swap_accepts,sampler.chain_N,sampler.chain_N);
		//spct++;
		

		load_temps_checkpoint_file(checkpoint_file, chain_temps, chain_N);
		sampler_output->populate_chain_temperatures(chain_temps);

		if(init){
			debugger_print(__FILE__,__LINE__,"Init structure");
			//###################################################
			//TESTING
			int temp_pos[sampler.chain_N];
			for(int z = 0; z<sampler.chain_N;z++){
				temp_pos[z] = sampler.chain_pos[z];
				if(sampler.restarted_chain[z]){
					temp_pos[z] = sampler.N_steps;
				}
			}
			sampler_output->populate_initial_output(temp_output, (int ***)NULL,(int **)NULL,sampler.ll_lp_output,temp_pos)	;
			//sampler_output->populate_initial_output(temp_output, (int ***)NULL,(int ***)NULL,sampler.ll_lp_output,sampler.chain_pos)	;
			//###################################################
			sampler_output->set_trim(0);	
			sampler_output->update_cold_chain_list();	
			init=false;
			debugger_print(__FILE__,__LINE__,"Finished init structure");
		}
		else{
			debugger_print(__FILE__,__LINE__,"Appending structure");
			//###################################################
			//TESTING
			int temp_pos[sampler.chain_N];
			for(int z = 0; z<sampler.chain_N;z++){
				temp_pos[z] = sampler.chain_pos[z];
				if(sampler.restarted_chain[z]){
					temp_pos[z] = sampler.N_steps;
				}
			}
			sampler_output->append_to_output(temp_output,(int ***)NULL,(int **)NULL,sampler.ll_lp_output,temp_pos)	;
			//sampler_output->append_to_output(temp_output,(int ***)NULL,(int ***)NULL,sampler.ll_lp_output,sampler.chain_pos)	;
			//###################################################
			debugger_print(__FILE__,__LINE__,"Finished appending structure");
		}
		sampler_output->calc_ac_vals(true);
		sampler_output->count_indep_samples(true);
		status = sampler_output->indep_samples;



		double ac_mean = 1;
		double pos_mean = 0;
		mean_list(sampler_output->max_acs, sampler_output->cold_chain_number, &ac_mean);
		double *temp_positions = new double[sampler_output->cold_chain_number];
		for(int i= 0 ; i<sampler_output->cold_chain_number; i++){
			temp_positions[i]=sampler_output->chain_lengths[sampler_output->cold_chain_ids[i]];	
		}
		mean_list(temp_positions, sampler_output->cold_chain_number, &pos_mean);
		delete [] temp_positions;
		if(relax){
			//Only considered burned in if the average (cold) chain length
			//is 500x the average ac (trimming as we go, for the ac
			debugger_print(__FILE__,__LINE__,std::string("Pos/ac: ")+std::to_string(pos_mean/ac_mean));
			if(pos_mean/ac_mean <100 && local_continue_burn){
			//if(false){
				sampler_output->set_trim(pos_mean);	
			}
			else{
				//##############################################
				//Evidence calculation
				if(calc_evidence){
					debugger_print(__FILE__,__LINE__,"Deallocating integrated likelihoods (if set)");
					sampler_output->dealloc_integrated_likelihoods();
					debugger_print(__FILE__,__LINE__,"Appending integrated likelihoods");
					sampler.chain_temps = &chain_temps[0];
					integrate_likelihood(&sampler);
					int ensemble_size= (int)(sampler.chain_N/coldchains);
					double integrated_likelihoods[ensemble_size];
					int integrated_likelihoods_terms[ensemble_size];
					combine_chain_evidence(&sampler, integrated_likelihoods, integrated_likelihoods_terms, ensemble_size);
					sampler_output->append_integrated_likelihoods(integrated_likelihoods, integrated_likelihoods_terms, ensemble_size);
					debugger_print(__FILE__,__LINE__,"Calculating evidence");
					sampler_output->calculate_evidence();
				}
				//##############################################
				
				relax=false;
				debugger_print(__FILE__,__LINE__,"Creating dump");
				//debugger_print(__FILE__,__LINE__,"Creating FULL data dump");
				//sampler_output->create_data_dump(false,false, chain_filename);
				sampler_output->create_data_dump(cold_chain_only_writeout,false, chain_filename);
				debugger_print(__FILE__,__LINE__,"Finished Creating dump");
			
				//debugger_print(__FILE__,__LINE__,"Creating FULL data dump");

				//sampler_output->create_data_dump(false,false, "data/test_full.hdf5");
			}
		}
		else{
			if(ac_mean > 3*ac_save){
				debugger_print(__FILE__,__LINE__,"Resetting trim");
				sampler_output->set_trim(pos_mean);	
				//sampler_output->create_data_dump(true,false, chain_filename);
				debugger_print(__FILE__,__LINE__,"Deallocating integrated likelihoods");
				sampler_output->dealloc_integrated_likelihoods();
			}

			//##############################################
			//Evidence calculation
			if(calc_evidence){
			debugger_print(__FILE__,__LINE__,"Appending integrated likelihoods");
			sampler.chain_temps = &chain_temps[0];
			integrate_likelihood(&sampler);
			int ensemble_size= (int)(sampler.chain_N/coldchains);
			double integrated_likelihoods[ensemble_size];
			int integrated_likelihoods_terms[ensemble_size];
			combine_chain_evidence(&sampler, integrated_likelihoods, integrated_likelihoods_terms, ensemble_size);
			sampler_output->append_integrated_likelihoods(integrated_likelihoods, integrated_likelihoods_terms, ensemble_size);
			debugger_print(__FILE__,__LINE__,"Calculating evidence");
			sampler_output->calculate_evidence();
			}
			//##############################################
			//else{
			debugger_print(__FILE__,__LINE__,"Appending dump");
			sampler_output->append_to_data_dump(chain_filename);
			//sampler_output->append_to_data_dump("data/test_full.hdf5");
			debugger_print(__FILE__,__LINE__,"Finished appending dump");
			//}

		}
		ac_save = ac_mean;
		max_ac_realloc = 0;
		mean_list(sampler_output->max_acs, sampler_output->cold_chain_number,&max_ac_realloc);
		std::cout<<"Average ac: "<<max_ac_realloc<<std::endl;
		std::cout<<"Independent samples: "<<sampler_output->indep_samples<<std::endl;
		//###########################################################3
		double temp_temps[sampler.chain_N];
		load_temps_checkpoint_file(checkpoint_file, temp_temps, sampler.chain_N);
		double swap_targets[2] = {.1,1};
		int cold_ids[sampler.chain_N];
		int cold_chains_ct=0;
		for(int k = 0 ; k<sampler.chain_N; k++){
			if(fabs(temp_temps[k] - 1.) < DOUBLE_COMP_THRESH){
				cold_ids[cold_chains_ct]=k;
				cold_chains_ct++;
			}	
		}
		int ensemble_members=sampler.chain_N;
		if(cold_chains_ct>1){ensemble_members = cold_ids[1]-cold_ids[0];}
		int ensemble_num = ceil((double)sampler.chain_N/ensemble_members);	
		
		double averages[ensemble_members];
		double average_accept_up[ensemble_members];
		int average_accept_up_attempts[ensemble_members];
		int average_accept_down_attempts[ensemble_members];
		double average_accept_down[ensemble_members];
		double average_attempts[ensemble_members];
		int average_chain_nums[ensemble_members];
		
		for(int i = 0 ; i<ensemble_members; i++){
			averages[i]=0;
			average_attempts[i]=0;
			average_accept_up[i]=0;
			average_accept_down_attempts[i]=0;
			average_accept_up_attempts[i]=0;
			average_accept_down[i]=0;
			average_chain_nums[i]=0;
		}
		for(int k = 0 ; k<sampler.chain_N; k++){
			averages[k%ensemble_members]+=(double)sampler.swap_accept_ct[k] / 
					(sampler.swap_accept_ct[k] + sampler.swap_reject_ct[k]);
			average_attempts[k%ensemble_members]+=(double)(sampler.swap_accept_ct[k] + sampler.swap_reject_ct[k]);
			average_chain_nums[k%ensemble_members]++;
			for(int l = 0 ; l<sampler.chain_N; l++){
				if(temp_temps[l] < temp_temps[k]){
					average_accept_down[k%ensemble_members]+=sampler.swap_accepts[k][l];
					average_accept_down_attempts[k%ensemble_members]+=sampler.swap_partners[k][l];

				}
				else if(temp_temps[l] > temp_temps[k]){
					average_accept_up[k%ensemble_members]+=sampler.swap_accepts[k][l];
					average_accept_up_attempts[k%ensemble_members]+=sampler.swap_partners[k][l];
				}

			}
		}
		
		bool reset = false;
		debugger_print(__FILE__,__LINE__,"Temp, Average [acceptance, to colder chains, to hotter chains], Average attempts");
		for(int i = 0 ; i<ensemble_members; i++){
			if(average_accept_up_attempts[i] != 0){
				average_accept_up[i] /= average_accept_up_attempts[i];
			}
			else{
				average_accept_up[i] =0;
			}
			if(average_accept_down_attempts[i] != 0){
				average_accept_down[i] /= average_accept_down_attempts[i];
			}
			else{
				average_accept_down[i] =0;
			}
			averages[i]/=average_chain_nums[i];
			average_attempts[i]/=average_chain_nums[i];
			//debugger_print(__FILE__,__LINE__,std::to_string(averages[i]) + " "+std::to_string(average_attempts[i]));
			debugger_print(__FILE__,__LINE__,
				std::to_string(temp_temps[i]) + " "+
				std::to_string(averages[i]) + " "+
				std::to_string(average_accept_down[i]) + " "+
				std::to_string(average_accept_up[i]) + " "+
				std::to_string(average_attempts[i])
			);
			if(search_iterations_ct < max_search_iterations){
				if(averages[i]<swap_targets[0] || averages[i]>swap_targets[1] ){
					debugger_print(__FILE__,__LINE__,"Average swap not in range -- restarting");
					reset=true;
				}
				else if ( (average_accept_up[i] < swap_targets[0] || average_accept_up[i] > swap_targets[1] ) && (i!= 0 && i!=(ensemble_members-1))){
					debugger_print(__FILE__,__LINE__,"Average swap to hot chains not in range -- restarting");
					reset=true;
				}
				else if ( (average_accept_down[i] < swap_targets[0] || average_accept_down[i] > swap_targets[1]) && (i!= 0 && i!=(ensemble_members-1))){
					debugger_print(__FILE__,__LINE__,"Average swap to cold chains not in range -- restarting");
					reset=true;
				}
			}
		}
		if (reset){
			if(!realloc){
				search_iterations_ct++;
			}

			realloc=true;
		}
		
		//###########################################################3
		deallocate_sampler_mem(&sampler);
			
		//if(false){
		if(temp_length < 10*max_ac_realloc){

			if(10*max_ac_realloc < max_chunk_size){
				deallocate_3D_array(temp_output, chain_N, temp_length, dimension);
				temp_length = 10*max_ac_realloc;
				temp_output = allocate_3D_array(chain_N,temp_length, dimension);
			}
			else{
				std::cout<<"WARNING -- hit maximum chunk size for a single sampler run"<<std::endl;
				std::cout<<"Independent samples per batch are projected to be "<<max_chunk_size/max_ac_realloc<<" and at least 1000 samples per AC calculation is recommended"<<std::endl;
				deallocate_3D_array(temp_output, chain_N, temp_length, dimension);
				temp_length = max_chunk_size;
				temp_output = allocate_3D_array(chain_N,temp_length, dimension);
			}
		}
		else if(temp_length>5000*max_ac_realloc){
			deallocate_3D_array(temp_output, chain_N, temp_length, dimension);
			temp_length = 5000*max_ac_realloc;
			temp_output = allocate_3D_array(chain_N,temp_length, dimension);	

		}
		max_ac_realloc=0;
		printProgress((double)status/N_steps);
	}
	
	//Maybe write new stat file format
	std::cout<<std::endl;
	if(calc_evidence){
		debugger_print(__FILE__,__LINE__,"Evidence: "+std::to_string(sampler_output->evidence));	
	}

	//Cleanup
	deallocate_3D_array(temp_output, chain_N, temp_length, dimension);
	deallocate_2D_array(temp_ac, dimension, corr_segments);
	return;
}
//######################################################################################
//######################################################################################
/*! \brief Routine to take a checkpoint file and begin a new chain at said checkpoint 
 *
 * See MCMC_MH_internal for more details of parameters (pretty much all the same)
 */
void continue_RJPTMCMC_MH_internal(sampler *samplerptr,
	std::string start_checkpoint_file,/**< File for starting checkpoint*/
	double ***output,/**< [out] output array, dimensions: output[chain_N][N_steps][dimension]*/
	int ***status,/**< [out] output parameter status array, dimensions: status[chain_N][N_steps][dimension]*/
	int **model_status,
	int nested_model_number,
	int N_steps,/**< Number of new steps to take*/
	int swp_freq,/**< frequency of swap attempts between temperatures*/
	std::function<double(double*, int*,int,mcmc_data_interface *,void *)> log_prior,/**< std::function for the log_prior function -- takes double *position, int *param_status, int dimension, int chain_id*/
	std::function<double(double*,int*,int,mcmc_data_interface *,void*)> log_likelihood,/**< std::function for the log_likelihood function -- takes double *position, int *param_status,int dimension, int chain_id*/
	std::function<void(double*,int*,int,double**,mcmc_data_interface *,void*)>fisher,/**< std::function for the fisher function -- takes double *position, int *param_status,int dimension, double **output_fisher, int chain_id*/
	std::function<void(double*,double*,int*,int*,int*,int*,double *,mcmc_data_interface *,void*)> RJ_proposal,/**< std::function for the log_likelihood function -- takes double *position, int *param_status,int dimension, int chain_id*/
	void **user_parameters,/**< Void pointer to any parameters the user may need inside log_prior, log_likelihood, or fisher. Should have one pointer for each chain. If this isn't needed, use (void**) NULL**/
	int numThreads,/**<Number of threads to use*/
	bool pool,/**<Boolean for whether to use ``deterministic'' vs ``stochastic'' sampling*/
	bool show_prog,/**< Boolean for whether to show progress or not (turn off for cluster runs*/
	bool update_RJ_width, 
	std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
	std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output -- if multiple cold chains, it will append each output to the other, and write out the total*/
	std::string auto_corr_filename,/**< Filename to output auto correlation in some interval, if empty string, not output*/
	std::string likelihood_log_filename,/**< Filename to write the log_likelihood and log_prior at each step -- use empty string to skip*/
	std::string end_checkpoint_file,/**< Filename to output data for checkpoint at the end of the continued run, if empty string, not saved*/
	bool tune, 	
	bool burn_phase 	
	)
{
	clock_t start, end, acend;
	double wstart, wend, wacend;
	start = clock();
	wstart = omp_get_wtime();

	samplerptr_global = samplerptr;
	samplerptr->tune = tune;
	samplerptr->burn_phase = burn_phase;

	//samplerptr->RJMCMC=true;
	samplerptr->RJMCMC = true;
	samplerptr->update_RJ_width=update_RJ_width;
	//if Fisher is not provided, Fisher and MALA steps
	//aren't used
	if(fisher ==NULL){
		samplerptr->fisher_exist = false;
	}
	else 
		samplerptr->fisher_exist = true;

	//if(likelihood_log_filename !=""){
	//	samplerptr->log_ll = true;
	//	samplerptr->log_lp = true;
	//}
	samplerptr->log_ll = true;
	samplerptr->log_lp = true;
	
	
	//Construct sampler structure
	samplerptr->lp = log_prior;
	samplerptr->ll = log_likelihood;
	samplerptr->fish = fisher;
	samplerptr->rj = RJ_proposal;
	samplerptr->swp_freq = 2;
	samplerptr->swap_rate = 1./swp_freq;
	samplerptr->N_steps = N_steps;
	samplerptr->show_progress = show_prog;
	samplerptr->num_threads = numThreads;
	samplerptr->user_parameters=user_parameters;
	
	samplerptr->output = output;
	samplerptr->param_status= status;
	samplerptr->model_status = model_status;
	samplerptr->nested_model_number = nested_model_number;
	samplerptr->pool = pool;

	//Unpack checkpoint file -- allocates memory internally -- separate call unneccessary
	load_checkpoint_file(start_checkpoint_file, samplerptr);



	//allocate other parameters
	for (int chain_index=0; chain_index<samplerptr->chain_N; chain_index++)
		assign_probabilities(samplerptr, chain_index);
	for (int j=0;j<samplerptr->chain_N;j++){
		samplerptr->current_likelihoods[j] =
			 samplerptr->ll(samplerptr->output[j][0],samplerptr->param_status[j][0],samplerptr->model_status[j][0],samplerptr->interfaces[j],samplerptr->user_parameters[j]);
	}
	if(samplerptr->log_ll && samplerptr->log_lp){
		for(int i = 0 ; i<samplerptr->chain_N; i++){
			samplerptr->ll_lp_output[i][0][0] = samplerptr->current_likelihoods[i];
			samplerptr->ll_lp_output[i][0][1] = samplerptr->lp(samplerptr->output[i][0],samplerptr->param_status[i][0],samplerptr->model_status[i][0],samplerptr->interfaces[i],samplerptr->user_parameters[i]);
		}
	}
	
	//Set chains with temp 1 to highest priority
	if(samplerptr->prioritize_cold_chains){
		for(int i =0 ;i<samplerptr->chain_N; i++){
			if(fabs(samplerptr->chain_temps[i]-1)<DOUBLE_COMP_THRESH)
				samplerptr->priority[i] = 0;
		}
	}
	
	//########################################################
	PTMCMC_MH_loop(samplerptr);	
	//##############################################################
	
	end =clock();
	wend =omp_get_wtime();

	samplerptr->time_elapsed_cpu = (double)(end-start)/CLOCKS_PER_SEC;
	samplerptr->time_elapsed_wall = (double)(wend-wstart);

	if(show_prog){
		std::cout<<std::endl;
	}
	
	//###########################################################
	//Auto-correlation
	if(auto_corr_filename != ""){
		std::cout<<"Calculating Autocorrelation: "<<std::endl;
		int segments = 50;
		double target_corr = .01;
		write_auto_corr_file_from_data(auto_corr_filename, samplerptr->output[0],samplerptr->N_steps,samplerptr->dimension,segments, target_corr, samplerptr->num_threads, false);
	}
	//###########################################################
	acend =clock();
	wacend =omp_get_wtime();
	samplerptr->time_elapsed_cpu_ac = (double)(acend-end)/CLOCKS_PER_SEC;
	samplerptr->time_elapsed_wall_ac = (double)(wacend - wend);

	//std::cout<<std::endl;
	//double accepted_percent = (double)(swp_accepted)/(swp_accepted+swp_rejected);
	//double rejected_percent = (double)(swp_rejected)/(swp_accepted+swp_rejected);
	//std::cout<<"Accepted percentage of chain swaps (all chains): "<<accepted_percent<<std::endl;
	//std::cout<<"Rejected percentage of chain swaps (all chains): "<<rejected_percent<<std::endl;
	//accepted_percent = (double)(step_accepted[0])/(step_accepted[0]+step_rejected[0]);
	//rejected_percent = (double)(step_rejected[0])/(step_accepted[0]+step_rejected[0]);
	//std::cout<<"Accepted percentage of steps (cold chain): "<<accepted_percent<<std::endl;
	//std::cout<<"Rejected percentage of steps (cold chain): "<<rejected_percent<<std::endl;
	//double nansum=0;
	//for (int i =0; i< chain_N; i++)
	//	nansum+= samplerptr->nan_counter[i];
	//std::cout<<"NANS in Fisher Calculations (all chains): "<<nansum<<std::endl;
	
	if(statistics_filename != "")
		write_stat_file(samplerptr, statistics_filename);
	
	if(chain_filename != "")
		write_output_file(chain_filename, samplerptr->N_steps, samplerptr->max_dim, samplerptr->output, samplerptr->param_status,samplerptr->model_status,samplerptr->chain_N,samplerptr->nested_model_number,samplerptr->chain_temps,true);

	if(end_checkpoint_file !=""){
		write_checkpoint_file(samplerptr, end_checkpoint_file);
	}

	//free(step_accepted);
	//free(step_rejected);
	//temps usually allocated by user, but for continued chains, this is done internally
	free(samplerptr->chain_temps);
	//deallocate_sampler_mem(samplerptr);
}
//######################################################################################
//######################################################################################
/*!\brief Generic reversable jump sampler, where the likelihood, prior, and reversable jump proposal are parameters supplied by the user
 *
 * Note: Using a min_dimension tells the sampler that there is a ``base model'', and that the dimensions from min_dim to max_dim are ``small'' corrections to that model. This helps inform some of the proposal algorithms and speeds up computation. If using discrete models with no overlap of variables (ie model A or model B), set min_dim to 0. Even if reusing certain parameters, if the extra dimensions don't describe ``small'' deviations, it's probably best to set min_dim to 0.Since the  RJ  proposal is user specified, even if there are parameters that should never be removed, it's up to the user to dictate that. Using min_dim will not affect that aspect of the  sampler. If there's a ``base-model'', the fisher function should produce a fisher matrix for the base model only. The modifications are then normally distributed around the last parameter value. Then the fisher function should take the minimum dimension instead of the maximum, like the other functions.
 *
 * Currently, no dynamic PT option, as it would be too many free parameters for the sampler to converge to a reasonable temperature distribution in a reasonable amount of time. Best use case, use the PTMCMC_MH_dynamic for the ``base'' dimension space, and use that temperature ladder.
 *
 * Base of the sampler, generic, with user supplied quantities for most of the samplers
 * properties
 * 	
 * Uses the Metropolis-Hastings method, with the option for Fisher/MALA steps if the Fisher
 * routine is supplied.
 *
 * 3 modes to use - 
 *
 * single threaded (numThreads = 1) runs single threaded
 *
 * multi-threaded ``deterministic'' (numThreads>1 ; pool = false) progresses each chain in parallel for swp_freq steps, then waits for all threads to complete before swapping temperatures in sequenctial order (j, j+1) then (j+1, j+2) etc (sequenctially)
 *
 * multi-threaded ``stochastic'' (numThreads>2 ; pool = true) progresses each chain in parallel by queueing each temperature and evaluating them in the order they were submitted. Once finished, the threads are queued to swap, where they swapped in the order they are submitted. This means the chains are swapped randomly, and the chains do NOT finish at the same time. The sampler runs until the the 0th chain reaches the step number
 *
 * Note on limits: In the prior function, if a set of parameters should be disallowed, return -std::numeric_limits<double>::infinity()  -- (this is in the <limits> file in std)
 *
 * The parameter array uses the dimensions [0,min_dim] always, and [min_dim, max_dim] in RJPTMCMC fashion
 *
 * Format for the auto_corr file (compatable with csv, dat, txt extensions): each row is a dimension of the cold chain, with the first row being the lengths used for the auto-corr calculation:
 *
 * lengths: length1 , length2 ...
 *
 * dim1: length1 , length2 ...
 *
 * .
 *
 * .
 *
 * .
 *
 *
 * Format for the chain file (compatable with csv, dat, txt extensions): each row is a step, each column a dimension:
 *
 * Step1: dim1 , dim2 , ..., max_dim, param_status1, param_status2, ...
 *
 * Step2: dim1 , dim2 , ..., max_dim, param_status1, param_status2, ...
 *
 * .
 *
 * .
 *
 * .
 *
 * Statistics_filename : should be txt extension
 *
 * checkpoint_file : This file saves the final position of all the chains, as well as other metadata, and can be loaded by the function <FUNCTION> to continue the chain from the point it left off. Not meant to be read by humans, the data order is custom to this software library. An empty string ("") means no checkpoint will be saved. For developers, the contents are:
 *
 * dimension, # of chains
 *
 * temps of chains
 *
 * Stepping widths of all chains
 *
 * Final position of all chains
 */
void RJPTMCMC_MH_internal(	double ***output, /**< [out] Output chains, shape is double[chain_N, N_steps,dimension]*/
	int ***parameter_status, /**< [out] Parameter status for each step corresponding to the output chains, shape is double[chain_N, N_steps,dimension]*/
	int **model_status,
	int nested_model_number,
	int max_dimension, 	/**< maximum dimension of the parameter space being explored -- only consideration is memory, as memory scales with dimension. Keep this reasonable, unless memory is REALLY not an issue*/
	int min_dimension, 	/**< minimum dimension of the parameter space being explored >=1*/
	int N_steps,	/**< Number of total steps to be taken, per chain*/
	int chain_N,	/**< Number of chains*/
	double *initial_pos, 	/**<Initial position in parameter space - shape double[dimension]*/
	int *initial_status, 	/**<Initial status of the parameters in the initial position in parameter space - shape int[max_dim]*/
	int initial_model_status, 	/**<Initial status of the parameters in the initial position in parameter space - shape int[max_dim]*/
	double *seeding_var, 	/**<Variance of the normal distribution used to seed each chain higher than 0 - shape double[max_dimension] -- initial seeding of zero corresponds to the dimension turned off initially*/
	double **ensemble_initial_pos,
	int **ensemble_initial_status,
	int *ensemble_initial_model_status,
	double *chain_temps,	/**<Double array of temperatures for the chains*/
	int swp_freq,	/**< the frequency with which chains are swapped*/
	std::function<double(double*, int*,int,mcmc_data_interface *,void *)> log_prior,/**< std::function for the log_prior function -- takes double *position, int *param_status, int dimension, int chain_id*/
	std::function<double(double*,int*,int,mcmc_data_interface *,void *)> log_likelihood,/**< std::function for the log_likelihood function -- takes double *position, int *param_status,int dimension, int chain_id*/
	std::function<void(double*,int*,int,double**,mcmc_data_interface *,void *)>fisher,/**< std::function for the fisher function -- takes double *position, int *param_status,int dimension, double **output_fisher, int chain_id*/
	std::function<void(double*,double*,int*,int*,int*,int*,double *,mcmc_data_interface *, void *)> RJ_proposal,/**< std::function for the log_likelihood function -- takes double *position, int *param_status,int dimension, int chain_id*/
	void **user_parameters,/**< Void pointer to any parameters the user may need inside log_prior, log_likelihood, or fisher. Should have one pointer for each chain. If this isn't needed, use (void**) NULL**/
	int numThreads, /**< Number of threads to use (=1 is single threaded)*/
	bool pool, /**< boolean to use stochastic chain swapping (MUST have >2 threads)*/
	bool show_prog, /**< boolean whether to print out progress (for example, should be set to ``false'' if submitting to a cluster)*/
	bool update_RJ_width, /**< boolean whether to print out progress (for example, should be set to ``false'' if submitting to a cluster)*/
	std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
	std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output -- if multiple cold chains, it will append each output to the other, and write out the total*/
	std::string auto_corr_filename,/**< Filename to output auto correlation in some interval, if empty string, not output*/
	std::string likelihood_log_filename,/**< Filename to write the log_likelihood and log_prior at each step -- use empty string to skip*/
	std::string checkpoint_file/**< Filename to output data for checkpoint, if empty string, not saved*/
	)
{
	//To Do:
	//	Retrofit the sampler_internal functions to accept dim-status array -- for regular PTMCMC, just submit with all 1's
	//		Includes statistics file, percentage of time spent on/off, etc
	//	retrofit PTMCMC/PTMCMC_dynamic_PT, thread pool functions
	
	//########################################################################	
	clock_t start, end, acend;
	double wstart, wend, wacend;
	start = clock();
	wstart = omp_get_wtime();

	sampler sampler_ref;
	sampler *samplerptr = &sampler_ref;
	samplerptr_global = &sampler_ref;

	//if Fisher is not provided, Fisher and MALA steps
	//aren't used
	if(fisher ==NULL){
		samplerptr->fisher_exist = false;
	}
	else 
		samplerptr->fisher_exist = true;

	samplerptr->RJMCMC=true;
	samplerptr->update_RJ_width=update_RJ_width;

	//if(likelihood_log_filename !=""){
	//	samplerptr->log_ll = true;
	//	samplerptr->log_lp = true;
	//}
	samplerptr->log_ll = true;
	samplerptr->log_lp = true;
	
	//Construct sampler structure
	samplerptr->lp = log_prior;
	samplerptr->ll = log_likelihood;
	samplerptr->rj = RJ_proposal;
	samplerptr->fish = fisher;
	samplerptr->swp_freq = 2;
	samplerptr->swap_rate = 1./swp_freq;
	samplerptr->chain_temps = chain_temps;
	samplerptr->chain_N = chain_N;
	samplerptr->N_steps = N_steps;
	samplerptr->max_dim = max_dimension;
	samplerptr->min_dim = min_dimension;
	samplerptr->show_progress = show_prog;
	samplerptr->num_threads = numThreads;


	samplerptr->output = output;
	samplerptr->param_status = parameter_status;
	samplerptr->model_status = model_status;
	samplerptr->nested_model_number = nested_model_number;
	samplerptr->pool = pool;
	samplerptr->numThreads = numThreads;
	samplerptr->user_parameters=user_parameters;

	allocate_sampler_mem(samplerptr);

	//########################################################
	//Set chains with temp 1 to highest priority
	if(samplerptr->prioritize_cold_chains){
		for(int i =0 ;i<samplerptr->chain_N; i++){
			if(fabs(samplerptr->chain_temps[i]-1)<DOUBLE_COMP_THRESH)
				samplerptr->priority[i] = 0;
		}
	}
	//########################################################

	for (int chain_index=0; chain_index<samplerptr->chain_N; chain_index++)
		assign_probabilities(samplerptr, chain_index);
	int  k=0;
	//int swp_accepted=0, swp_rejected=0;
	int *step_accepted = (int *)malloc(sizeof(int)*samplerptr->chain_N);
	int *step_rejected = (int *)malloc(sizeof(int)*samplerptr->chain_N);
	for(int j=0; j<samplerptr->chain_N; j++){
		step_accepted[j]=0;
		step_rejected[j]=0;
	}

	
	//Assign initial position to start chains
	assign_initial_pos(samplerptr, initial_pos, initial_status,initial_model_status,ensemble_initial_pos,ensemble_initial_status,ensemble_initial_model_status,seeding_var);	

		
	PTMCMC_MH_loop(samplerptr);	
	
	//############################################################
	//Write ll lp to file
	//if(samplerptr->log_ll && samplerptr->log_lp){
	//	write_file(likelihood_log_filename,samplerptr->ll_lp_output[0],samplerptr->N_steps,2);
	//	//write_file(likelihood_log_filename,samplerptr->ll_lp_output[samplerptr->chain_N-1],samplerptr->N_steps,2);
	//}
	//############################################################
	//##############################################################
	int swp_accepted=0, swp_rejected=0;
	for (int i =0;i<samplerptr->chain_N; i++)
	{
		swp_accepted+=samplerptr->swap_accept_ct[i];
		swp_rejected+=samplerptr->swap_reject_ct[i];
		step_accepted[i]+=samplerptr->step_accept_ct[i];
		step_rejected[i]+=samplerptr->step_reject_ct[i];
	}
	end =clock();
	wend =omp_get_wtime();

	samplerptr->time_elapsed_cpu = (double)(end-start)/CLOCKS_PER_SEC;
	samplerptr->time_elapsed_wall = (double)(wend-wstart);

	if(show_prog){
		std::cout<<std::endl;
	}
	
	//###########################################################
	//Auto-correlation
	if(auto_corr_filename != ""){
		std::cout<<"Calculating Autocorrelation: "<<std::endl;
		int segments = 50;
		double target_corr = .01;
		write_auto_corr_file_from_data(auto_corr_filename, samplerptr->output[0],samplerptr->N_steps,samplerptr->max_dim,segments, target_corr, samplerptr->num_threads, false);
	}
	//###########################################################
	acend =clock();
	wacend =omp_get_wtime();
	samplerptr->time_elapsed_cpu_ac = (double)(acend-end)/CLOCKS_PER_SEC;
	samplerptr->time_elapsed_wall_ac = (double)(wacend - wend);

	std::cout<<std::endl;
	double accepted_percent = (double)(swp_accepted)/(swp_accepted+swp_rejected);
	double rejected_percent = (double)(swp_rejected)/(swp_accepted+swp_rejected);
	std::cout<<"Accepted percentage of chain swaps (all chains): "<<accepted_percent<<std::endl;
	std::cout<<"Rejected percentage of chain swaps (all chains): "<<rejected_percent<<std::endl;
	accepted_percent = (double)(step_accepted[0])/(step_accepted[0]+step_rejected[0]);
	rejected_percent = (double)(step_rejected[0])/(step_accepted[0]+step_rejected[0]);
	std::cout<<"Accepted percentage of steps (cold chain): "<<accepted_percent<<std::endl;
	std::cout<<"Rejected percentage of steps (cold chain): "<<rejected_percent<<std::endl;
	double nansum=0;
	for (int i =0; i< chain_N; i++)
		nansum+= samplerptr->nan_counter[i];
	std::cout<<"NANS in Fisher Calculations (all chains): "<<nansum<<std::endl;
	
	if(statistics_filename != "")
		write_stat_file(samplerptr, statistics_filename);
	
	if(chain_filename != ""){
		write_output_file(chain_filename, samplerptr->N_steps, samplerptr->max_dim, samplerptr->output, samplerptr->param_status,samplerptr->model_status,samplerptr->chain_N,samplerptr->nested_model_number,samplerptr->chain_temps, true);
	}

	if(checkpoint_file !=""){
		write_checkpoint_file(samplerptr, checkpoint_file);
	}

	free(step_accepted);
	free(step_rejected);
	deallocate_sampler_mem(samplerptr);
		
}

//######################################################################################
//######################################################################################
/*! \brief Continue dyanmically tunes an MCMC for optimal spacing. step width, and chain number
 *
 * NOTE: nu, and t0 parameters determine the dynamics, so these are important quantities. nu is related to how many swap attempts it takes to substantially change the temperature ladder, why t0 determines the length of the total dyanimcally period. Moderate initial choices would be 10 and 1000, respectively.
 *
 * Based on arXiv:1501.05823v3
 *
 * Currently, Chain number is fixed
 *
 * max_chain_N_thermo_ensemble sets the maximium number of chains to use to in successively hotter chains to cover the likelihood surface while targeting an optimal swap acceptance target_swp_acc. 
 *
 * max_chain_N determines the total number of chains to run once thermodynamic equilibrium has been reached. This results in chains being added after the initial PT dynamics have finished according to chain_distribution_scheme.
 *
 * If no preference, set max_chain_N_thermo_ensemble = max_chain_N = numThreads = (number of cores (number of threads if hyperthreaded))-- this will most likely be the most optimal configuration. If the number of cores on the system is low, you may want to use n*numThreads for some integer n instead, depending on the system.
 *
 * chain_distribution_scheme:
 *
 * "cold": All chains are added at T=1 (untempered)
 *
 * "refine": Chains are added between the optimal temps geometrically -- this may be a good option as it will be a good approximation of the ideal distribution of chains, while keeping the initial dynamical time low 
 *
 * "double": Chains are added in order of rising temperature that mimic the distribution achieved by the earier PT dynamics
 *
 * "half_ensemble": For every cold chain added, half of the ensemble is added again. Effectively, two cold chains for every ensemble
 */
void continue_PTMCMC_MH_dynamic_PT_alloc_internal(std::string checkpoint_file_start,
	double ***output, /**< [out] Output chains, shape is double[max_chain_N, N_steps,dimension]*/
	int N_steps,	/**< Number of total steps to be taken, per chain AFTER chain allocation*/
	int max_chain_N_thermo_ensemble,/**< Maximum number of chains to use in the thermodynamic ensemble (may use less)*/
	double *chain_temps, /**<[out] Final chain temperatures used -- should be shape double[chain_N]*/
	int swp_freq,	/**< the frequency with which chains are swapped*/
	int t0,/**< Time constant of the decay of the chain dynamics  (~1000)*/
	int nu,/**< Initial amplitude of the dynamics (~100)*/
	std::string chain_distribution_scheme, /*How to allocate the remaining chains once equilibrium is reached*/
	std::function<double(double*,int* ,int,mcmc_data_interface *,void *)> log_prior,/**< std::function for the log_prior function -- takes double *position, int dimension, int chain_id*/
	std::function<double(double*,int*,int,mcmc_data_interface *,void *)> log_likelihood,/**< std::function for the log_likelihood function -- takes double *position, int dimension, int chain_id*/
	std::function<void(double*,int*,int,double**,mcmc_data_interface *,void *)>fisher,/**< std::function for the fisher function -- takes double *position, int dimension, double **output_fisher, int chain_id*/
	void **user_parameters,/**< Void pointer to any parameters the user may need inside log_prior, log_likelihood, or fisher. Should have one pointer for each chain. If this isn't needed, use (void**) NULL**/
	int numThreads, /**< Number of threads to use (=1 is single threaded)*/
	bool pool, /**< boolean to use stochastic chain swapping (MUST have >2 threads)*/
	bool show_prog, /**< boolean whether to print out progress (for example, should be set to ``false'' if submitting to a cluster)*/
	bool dynamic_chain_number,
	std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
	std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output*/
	std::string checkpoint_file,/**< Filename to output data for checkpoint, if empty string, not saved*/
	bool burn_phase 
	)
{
	//std::cout<<"MEM CHECK : start continue"<<std::endl;
	clock_t start, end, acend;
	double wstart, wend, wacend;
	start = clock();
	wstart = omp_get_wtime();

	sampler samplerobj;
	sampler *samplerptr = &samplerobj;
	samplerptr_global = &samplerobj;

	//if Fisher is not provided, Fisher and MALA steps
	//aren't used
	if(fisher ==NULL){
		samplerptr->fisher_exist = false;
	}
	else 
		samplerptr->fisher_exist = true;
	
	samplerptr->log_ll = false;
	samplerptr->log_lp = false;
	
	//Construct sampler structure
	samplerptr->lp = log_prior;
	samplerptr->ll = log_likelihood;
	samplerptr->fish = fisher;
	samplerptr->swp_freq = 2.;
	samplerptr->swap_rate = 1./swp_freq;
	//For PT dynamics
	samplerptr->N_steps = N_steps;

	samplerptr->num_threads = numThreads;
	samplerptr->output =output;
	samplerptr->user_parameters=user_parameters;
	samplerptr->burn_phase = burn_phase;
	if(burn_phase){
		samplerptr->isolate_ensembles=false;	
		samplerptr->isolate_ensembles_cold = false;
	}

	load_checkpoint_file(checkpoint_file_start,samplerptr);

	//During chain allocation, pooling isn't used
	samplerptr->pool = false;
	samplerptr->numThreads = numThreads;
	samplerptr->A = new int[samplerptr->chain_N];
	samplerptr->PT_alloc = true;
	
	for (int chain_index=0; chain_index<samplerptr->chain_N; chain_index++)
		assign_probabilities(samplerptr, chain_index);
	
	for (int j=0;j<samplerptr->chain_N;j++){
		samplerptr->current_likelihoods[j] =
			 samplerptr->ll(samplerptr->output[j][0],samplerptr->param_status[j][0],samplerptr->model_status[j][0],samplerptr->interfaces[j],samplerptr->user_parameters[j]);
	}
	
	//Set chains with temp 1 to highest priority
	if(samplerptr->prioritize_cold_chains){
		for(int i =0 ;i<samplerptr->chain_N; i++){
			if(fabs(samplerptr->chain_temps[i]-1)<DOUBLE_COMP_THRESH)
				samplerptr->priority[i] = 0;
		}
	}
	int chain_N_max = samplerptr->chain_N;
	int ct = 1  ;
	for(int i = 1 ; i<samplerptr->chain_N; i++){
		if(samplerptr->chain_temps[i]!=1){
			ct++;	
		}
		else{
			break;
		}
	}	
	samplerptr->chain_N = ct;
	
	
	//#########################################################################
	//#########################################################################
	
	dynamic_temperature_internal(samplerptr, N_steps, nu, t0,swp_freq, max_chain_N_thermo_ensemble, dynamic_chain_number,show_prog);

	//#######################################################################
	//#######################################################################
	end =clock();
	wend =omp_get_wtime();

	samplerptr->time_elapsed_cpu = (double)(end-start)/CLOCKS_PER_SEC;
	samplerptr->time_elapsed_wall = (double)(wend-wstart);

	
	acend =clock();
	wacend =omp_get_wtime();
	samplerptr->time_elapsed_cpu_ac = (double)(acend-end)/CLOCKS_PER_SEC;
	samplerptr->time_elapsed_wall_ac = (double)(wacend - wend);

	//#################################################################
	//
	//Copy sampler to new sampler and run loop, skip checkpoint nonsense
	//
	//#################################################################
	
	sampler static_sampler;
	static_sampler.A = new int[chain_N_max];
	static_sampler.chain_temps = new double[chain_N_max];
	initiate_full_sampler(&static_sampler, samplerptr, max_chain_N_thermo_ensemble,chain_N_max, chain_distribution_scheme,checkpoint_file_start);

	if(statistics_filename != "")
		write_stat_file(samplerptr, statistics_filename);
	
	if(chain_filename != "")
		write_output_file(chain_filename, samplerptr->N_steps, samplerptr->max_dim, samplerptr->output, samplerptr->param_status,samplerptr->model_status,samplerptr->chain_N,samplerptr->nested_model_number,samplerptr->chain_temps,false);

	delete [] samplerptr->A;
	//If chains were added or removed, set chain N back to max for deallocation
	samplerptr->chain_N = chain_N_max;
	//samplerptr->chain_N = max_chain_N_thermo_ensemble;
	deallocate_sampler_mem(samplerptr);
	//delete [] samplerptr->chain_temps;
	free(samplerptr->chain_temps);

	static_sampler.show_progress=show_prog;
	static_sampler.pool=pool;

	write_checkpoint_file(&static_sampler, checkpoint_file);
	for(int j = 0; j<static_sampler.chain_N; j++){
		chain_temps[j] = static_sampler.chain_temps[j];
	}
	delete [] static_sampler.A;
	deallocate_sampler_mem(&static_sampler);
	delete [] static_sampler.chain_temps;
	//std::cout<<"MEM CHECK : end continue"<<std::endl;
}
/*! \brief Dyanmically tunes an MCMC for optimal spacing. step width, and chain number
 *
 * NOTE: nu, and t0 parameters determine the dynamics, so these are important quantities. nu is related to how many swap attempts it takes to substantially change the temperature ladder, why t0 determines the length of the total dyanimcally period. Moderate initial choices would be 10 and 1000, respectively.
 *
 * Based on arXiv:1501.05823v3
 *
 * Currently, Chain number is fixed
 *
 * max_chain_N_thermo_ensemble sets the maximium number of chains to use to in successively hotter chains to cover the likelihood surface while targeting an optimal swap acceptance target_swp_acc. 
 *
 * max_chain_N determines the total number of chains to run once thermodynamic equilibrium has been reached. This results in chains being added after the initial PT dynamics have finished according to chain_distribution_scheme.
 *
 * If no preference, set max_chain_N_thermo_ensemble = max_chain_N = numThreads = (number of cores (number of threads if hyperthreaded))-- this will most likely be the most optimal configuration. If the number of cores on the system is low, you may want to use n*numThreads for some integer n instead, depending on the system.
 *
 * chain_distribution_scheme:
 *
 * "cold": All chains are added at T=1 (untempered)
 *
 * "refine": Chains are added between the optimal temps geometrically -- this may be a good option as it will be a good approximation of the ideal distribution of chains, while keeping the initial dynamical time low 
 *
 * "double": Chains are added in order of rising temperature that mimic the distribution achieved by the earier PT dynamics
 *
 * "half_ensemble": For every cold chain added, half of the ensemble is added again. Effectively, two cold chains for every ensemble
 */
void PTMCMC_MH_dynamic_PT_alloc_internal(double ***output, /**< [out] Output chains, shape is double[max_chain_N, N_steps,dimension]*/
	int dimension, 	/**< dimension of the parameter space being explored*/
	int N_steps,	/**< Number of total steps to be taken, per chain AFTER chain allocation*/
	int chain_N,/**< Maximum number of chains to use */
	int max_chain_N_thermo_ensemble,/**< Maximum number of chains to use in the thermodynamic ensemble (may use less)*/
	double *initial_pos, 	/**<Initial position in parameter space - shape double[dimension]*/
	double *seeding_var, 	/**<Variance of the normal distribution used to seed each chain higher than 0 - shape double[dimension]*/
	double **ensemble_initial_pos,
	double *chain_temps, /**<[out] Final chain temperatures used -- should be shape double[chain_N]*/
	int swp_freq,	/**< the frequency with which chains are swapped*/
	int t0,/**< Time constant of the decay of the chain dynamics  (~1000)*/
	int nu,/**< Initial amplitude of the dynamics (~100)*/
	std::string chain_distribution_scheme, /*How to allocate the remaining chains once equilibrium is reached*/
	std::function<double(double*,int* ,int,mcmc_data_interface *,void *)> log_prior,/**< std::function for the log_prior function -- takes double *position, int dimension, int chain_id*/
	std::function<double(double*,int*,int,mcmc_data_interface *,void *)> log_likelihood,/**< std::function for the log_likelihood function -- takes double *position, int dimension, int chain_id*/
	std::function<void(double*,int*,int,double**,mcmc_data_interface *,void *)>fisher,/**< std::function for the fisher function -- takes double *position, int dimension, double **output_fisher, int chain_id*/
	void **user_parameters,/**< Void pointer to any parameters the user may need inside log_prior, log_likelihood, or fisher. Should have one pointer for each chain. If this isn't needed, use (void**) NULL**/
	int numThreads, /**< Number of threads to use (=1 is single threaded)*/
	bool pool, /**< boolean to use stochastic chain swapping (MUST have >2 threads)*/
	bool show_prog, /**< boolean whether to print out progress (for example, should be set to ``false'' if submitting to a cluster)*/
	bool dynamic_chain_number,
	std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
	std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output*/
	std::string checkpoint_file,/**< Filename to output data for checkpoint, if empty string, not saved*/
	bool burn_phase
	)
{
	int initial_model_status=0;

	clock_t start, end, acend;
	double wstart, wend, wacend;
	start = clock();
	wstart = omp_get_wtime();

	sampler samplerobj;
	sampler *samplerptr = &samplerobj;
	samplerptr_global = &samplerobj;

	//if Fisher is not provided, Fisher and MALA steps
	//aren't used
	if(fisher ==NULL){
		samplerptr->fisher_exist = false;
	}
	else 
		samplerptr->fisher_exist = true;
	
	samplerptr->log_ll = false;
	samplerptr->log_lp = false;
	
	//Construct sampler structure
	samplerptr->lp = log_prior;
	samplerptr->ll = log_likelihood;
	samplerptr->fish = fisher;
	samplerptr->swp_freq = 2.;
	samplerptr->swap_rate = 1./swp_freq;
	//For PT dynamics
	samplerptr->N_steps = N_steps;
	samplerptr->burn_phase = burn_phase;
	if(burn_phase){
		samplerptr->isolate_ensembles=false;	
		samplerptr->isolate_ensembles_cold = false;
	}

	samplerptr->dimension =dimension;
	samplerptr->min_dim=dimension;
	samplerptr->max_dim= dimension;
	samplerptr->num_threads = numThreads;
	samplerptr->output =output;
	samplerptr->user_parameters=user_parameters;

	//Start out with geometrically spaced chain
	samplerptr->chain_temps =new double [max_chain_N_thermo_ensemble];
	samplerptr->chain_temps[0] = 1.;
	samplerptr->chain_N = max_chain_N_thermo_ensemble;
	assign_ensemble_temps(samplerptr->chain_temps, max_chain_N_thermo_ensemble,max_chain_N_thermo_ensemble,TMAX);
	//double temp_step = pow(TMAX/1.,1./max_chain_N_thermo_ensemble);
	//for(int i = 1; i<samplerptr->chain_N-1; i++){
	//	samplerptr->chain_temps[i] = temp_step * samplerptr->chain_temps[i-1];
	//}
	////Set last chain essentially at infinity
	//samplerptr->chain_temps[samplerptr->chain_N -1] = TMAX;


	//During chain allocation, pooling isn't used
	samplerptr->pool = false;
	samplerptr->numThreads = numThreads;
	samplerptr->A = new int[samplerptr->chain_N];
	for(int i =0 ; i<samplerptr->chain_N; i++){
		samplerptr->A[i]=0;
	}
	samplerptr->PT_alloc = true;

	allocate_sampler_mem(samplerptr);


	for (int chain_index=0; chain_index<samplerptr->chain_N; chain_index++)
		assign_probabilities(samplerptr, chain_index);
	
	int  k=0;
	//Assign initial position to start chains
	int *initial_status = new int[samplerptr->max_dim];
	for(int i = 0; i<samplerptr->max_dim; i++){
		initial_status[i]=1;	
	}
	int *ensemble_initial_model_status = NULL;
	int **ensemble_initial_status = new int*[samplerptr->chain_N];
	for(int j = 0 ; j<samplerptr->chain_N; j++){
		ensemble_initial_status[j] = new int[samplerptr->max_dim];
		for(int i = 0; i<samplerptr->max_dim; i++){
			ensemble_initial_status[j][i]=1;	
		}
	}
	assign_initial_pos(samplerptr, initial_pos,initial_status,initial_model_status,ensemble_initial_pos,ensemble_initial_status,ensemble_initial_model_status,seeding_var);	
	delete [] initial_status;
	for(int j = 0 ; j<samplerptr->chain_N; j++){
		delete [] ensemble_initial_status[j];
	}
	delete [] ensemble_initial_status;
	
	
	//#########################################################################
	//#########################################################################
	
	dynamic_temperature_internal(samplerptr, N_steps, nu, t0,swp_freq, max_chain_N_thermo_ensemble,dynamic_chain_number, show_prog);

	//#######################################################################
	//#######################################################################
	end =clock();
	wend =omp_get_wtime();

	samplerptr->time_elapsed_cpu = (double)(end-start)/CLOCKS_PER_SEC;
	samplerptr->time_elapsed_wall = (double)(wend-wstart);

	
	acend =clock();
	wacend =omp_get_wtime();
	samplerptr->time_elapsed_cpu_ac = (double)(acend-end)/CLOCKS_PER_SEC;
	samplerptr->time_elapsed_wall_ac = (double)(wacend - wend);

	//#################################################################
	//
	//Copy sampler to new sampler and run loop, skip checkpoint nonsense
	//
	//#################################################################
	
	
	sampler static_sampler;
	static_sampler.A = new int[chain_N];
	static_sampler.chain_temps = new double[chain_N];
	initiate_full_sampler(&static_sampler, samplerptr, max_chain_N_thermo_ensemble,chain_N, chain_distribution_scheme,"");

	if(statistics_filename != "")
		write_stat_file(samplerptr, statistics_filename);
	
	if(chain_filename != "")
		write_output_file(chain_filename, samplerptr->N_steps, samplerptr->max_dim, samplerptr->output, samplerptr->param_status,samplerptr->model_status,samplerptr->chain_N,samplerptr->nested_model_number,samplerptr->chain_temps,false);

	delete [] samplerptr->A;
	//If chains were added or removed, set chain N back to max for deallocation
	samplerptr->chain_N = max_chain_N_thermo_ensemble;
	deallocate_sampler_mem(samplerptr);
	delete [] samplerptr->chain_temps;

	static_sampler.show_progress=show_prog;
	static_sampler.pool=pool;

	write_checkpoint_file(&static_sampler, checkpoint_file);
	for(int j = 0; j<static_sampler.chain_N; j++){
		chain_temps[j] = static_sampler.chain_temps[j];
	}
	delete [] static_sampler.A;
	deallocate_sampler_mem(&static_sampler);
	delete [] static_sampler.chain_temps;
}

void dynamic_temperature_internal(sampler *samplerptr, int N_steps, double nu, int t0,int swp_freq, int max_chain_N_thermo_ensemble, bool dynamic_chain_number, bool show_prog)
{
	
	//NOTE: instead of dynamics, use variance over accept ratios over \nu steps
	//Average percent change in temperature 
	double ave_dynamics= 1.;
	//tolerance in percent change of temperature to determine equilibrium
	double tolerance= .001;
	int stability_ct = 0;
	int stability_tol = 10;
	//Frequency to check for equilibrium
	int equilibrium_check_freq=2*nu;
	
	double *old_temps = new double[max_chain_N_thermo_ensemble];
	int ensemble_size = 1;
	std::cout<<samplerptr->chain_temps[0]<<std::endl;
	for(int i =1; i<samplerptr->chain_N; i++){
		if(fabs(samplerptr->chain_temps[i] - 1 ) < 1e-10){break;}
		std::cout<<samplerptr->chain_temps[i]<<std::endl;
		ensemble_size ++;
	}
	for(int i =0; i<samplerptr->chain_N; i++){
		old_temps[i]=samplerptr->chain_temps[i];
	}


	//Frequency to update chain number
	int chain_pop_update_freq=5*nu;
	int pop_check_var=chain_pop_update_freq;
	//Boundaries that mark acceptable average acceptance ratios 
	//for chain swapping in a chain population
	//double chain_pop_high = .7;
	//double chain_pop_low = .4;
	double chain_pop_high = .4;
	double chain_pop_low = .3;

	//Keep track of acceptance ratio in chuncks
	int *running_accept_ct = new int[max_chain_N_thermo_ensemble];
	int *running_reject_ct = new int[max_chain_N_thermo_ensemble];
	int *prev_reject_ct = new int[max_chain_N_thermo_ensemble];
	int *prev_accept_ct = new int[max_chain_N_thermo_ensemble];
	double *running_ratio = new double[max_chain_N_thermo_ensemble];
	//bool dynamic_chain_num = true;
	//bool dynamic_chain_num = false;

	bool chain_pop_target_reached = false;
	for(int i =0; i<max_chain_N_thermo_ensemble; i++){
		running_accept_ct[i] = 0;
		running_reject_ct[i] = 0;
		prev_accept_ct[i] = 0;
		prev_reject_ct[i] = 0;
		running_ratio[i] = 0;
	}

	//For each loop, we walk forward till one more swap has happened, then we update temps
	//samplerptr->N_steps = swp_freq;
	std::cout<<"Dynamical PT allocation (measured by average percent change in temperature): "<<std::endl;
	
	int t = 0;
	samplerptr->show_progress = false;
	bool testing=false;
	int test_ct=max_chain_N_thermo_ensemble - samplerptr->chain_N;
	while( t <= (N_steps-equilibrium_check_freq))
	{
		chain_pop_target_reached = false;
		if(ave_dynamics<tolerance && chain_pop_target_reached){
			if(stability_ct > stability_tol)
				break;
			else
				stability_ct++;
		}
		else stability_ct = 0;
		//step equilibrium_check_freq
		for(int i =0; i<equilibrium_check_freq/swp_freq; i++){
			//Reset A values
			//for(int i =0 ; i<samplerptr->chain_N; i++){
			//	samplerptr->A[i]=0;
			//}
			//Now that swapping is stochastic, we need to make sure
			//a swap actually happened
			PTMCMC_MH_step_incremental(samplerptr, samplerptr->swp_freq);	
			t+= samplerptr->swp_freq;
			//t+= (int)(1./samplerptr->swap_rate)*2;
			//std::cout<<"MEM CHECK : leaving MCMC steps"<<std::endl;
			//Move temperatures
			update_temperatures(samplerptr, t0, nu, t);

			if(dynamic_chain_number){
				if(pop_check_var<t){
					pop_check_var += chain_pop_update_freq;
					for(int i =0 ;i < samplerptr->chain_N; i++){
						running_accept_ct[i] = 
							samplerptr->swap_accept_ct[i] 
							- prev_accept_ct[i];
						running_reject_ct[i] = 
							samplerptr->swap_reject_ct[i] 
							- prev_reject_ct[i];
						running_ratio[i] = ((double)running_accept_ct[i])/
							(running_accept_ct[i] 
							+ running_reject_ct[i]);	
					}
					double ave_accept = 0;
					for(int i =0; i<samplerptr->chain_N; i++){
						ave_accept+=running_ratio[i];
					}
					ave_accept/=samplerptr->chain_N;
					if(ave_accept < chain_pop_high && ave_accept> chain_pop_low){
						chain_pop_target_reached = true;
					}
					else {
						chain_pop_target_reached = false;
						if(ave_accept<chain_pop_low && samplerptr->chain_N < max_chain_N_thermo_ensemble){
						//if(testing && samplerptr->chain_N<max_chain_N_thermo_ensemble){
						//	if(test_ct!=0){
						//		test_ct--;	
						//	}
						//	else{
						//		testing = false;
						//		test_ct = 0;
						//	}
							//add chain
							//std::cout<<"MEM CHECK : entering add chain"<<std::endl;
							int min_id =1;
							int min_val =running_ratio[1];
							//Don't add chain 0 and chain chain_N-1
							for (int j =1 ;j <samplerptr->chain_N-1; j++){
								if(running_ratio[j]<min_val){
									min_id = j;
									min_val = running_ratio[j];
								}
							}	
							//std::cout<<"MEM CHECK : add chain "<<min_id<<std::endl;
							samplerptr->chain_N++;	
							for(int i = samplerptr->chain_N-2; i>=min_id; i--){
								transfer_chain(samplerptr,samplerptr, i+1, i, true);	
								running_accept_ct[i+1] = running_accept_ct[i];
								running_reject_ct[i+1] = running_reject_ct[i];
								prev_accept_ct[i+1] = prev_accept_ct[i];
								prev_reject_ct[i+1] = prev_reject_ct[i];
								old_temps[i+1] = old_temps[i];
							}
							//Add new chain between two other chains, geometrically
							samplerptr->chain_temps[min_id] = std::sqrt(samplerptr->chain_temps[min_id-1]*samplerptr->chain_temps[min_id+1]);
							old_temps[min_id] = samplerptr->chain_temps[min_id];
							//populate all the necessary chain-specific parameters
							for(int i =0 ;i<samplerptr->max_dim; i++){
								//Just use the chain's last 
								//position immediately
								//below min_id 
								samplerptr->output[min_id][0][i] = samplerptr->output[min_id-1][samplerptr->chain_pos[min_id-1]][i];
								samplerptr->param_status[min_id][0][i] = samplerptr->param_status[min_id-1][samplerptr->chain_pos[min_id-1]][i];
							}
							samplerptr->current_likelihoods[min_id] = samplerptr->ll(samplerptr->output[min_id][0],samplerptr->param_status[min_id][0],samplerptr->model_status[min_id][0],samplerptr->interfaces[min_id],samplerptr->user_parameters[min_id]);
							samplerptr->current_hist_pos[min_id] = 0;
							samplerptr->chain_pos[min_id] = 0;
							samplerptr->de_primed[min_id]=false;
							samplerptr->gauss_last_accept_ct[min_id] = 0;
							samplerptr->gauss_last_reject_ct[min_id]= 0;
							samplerptr->de_last_accept_ct[min_id] = 0;
							samplerptr->de_last_reject_ct[min_id] = 0;
							samplerptr->fish_last_accept_ct[min_id] = 0;
							samplerptr->fish_last_reject_ct[min_id] = 0;

							samplerptr->gauss_accept_ct[min_id] = 0;
							samplerptr->gauss_reject_ct[min_id] = 0;
							samplerptr->de_accept_ct[min_id] = 0;
							samplerptr->de_reject_ct[min_id] = 0;
							samplerptr->fish_accept_ct[min_id] =0;
							samplerptr->fish_reject_ct[min_id] = 0;
							samplerptr->mmala_accept_ct[min_id] = 0;
							samplerptr->mmala_reject_ct[min_id] = 0;
							assign_probabilities(samplerptr, min_id);	
							samplerptr->fisher_update_ct[min_id] = samplerptr->fisher_update_number;
							samplerptr->waiting[min_id] = true;
							samplerptr->waiting_SWP[min_id] = false;
							samplerptr->priority[min_id] =1;
							samplerptr->ref_chain_status[min_id]=true;
							samplerptr->nan_counter[min_id] = 0;
							samplerptr->num_gauss[min_id] =0;
							samplerptr->num_fish[min_id] = 0;
							samplerptr->num_de[min_id] = 0;
							samplerptr->num_mmala[min_id] = 0;
							samplerptr->swap_accept_ct[min_id] = 0;
							samplerptr->swap_reject_ct[min_id] = 0;
							samplerptr->step_accept_ct[min_id] = 0;
							samplerptr->step_reject_ct[min_id] = 0;
							samplerptr->prop_MH_factor[min_id]=0;
							if(samplerptr->log_ll){
								samplerptr->ll_lp_output[min_id][0][0] = samplerptr->current_likelihoods[min_id];
							}
							if(samplerptr->log_lp){
								samplerptr->ll_lp_output[min_id][0][1] = samplerptr->lp(samplerptr->output[min_id][0], samplerptr->param_status[min_id][0],samplerptr->model_status[min_id][0],samplerptr->interfaces[min_id],samplerptr->user_parameters[min_id]);
							}
							if(samplerptr->PT_alloc)
								samplerptr->A[min_id] = 0;
							if(samplerptr->fisher_exist){
								update_fisher(samplerptr, samplerptr->output[min_id][0], samplerptr->param_status[min_id][0],&(samplerptr->model_status[min_id][0]),min_id);	
							}
							//std::cout<<"MEM CHECK : leaving add chain"<<std::endl;
						
						}
						else if (ave_accept>chain_pop_high && samplerptr->chain_N>3){
						//else if(!testing){
						//	if(test_ct<10){
						//		test_ct++;
						//	}
						//	else{
						//		testing = true;
						//	}
						//	std::cout<<"rm "<<samplerptr->chain_N-1<<std::endl;
						//	std::cout<<"MEM CHECK : entering rm chain"<<std::endl;
							//remove chain
							int max_id =0;
							int max_val =-1;
							//Don't remove chain 0 and chain chain_N-1
							for (int j =1 ;j <samplerptr->chain_N-1; j++){
								if(running_ratio[j]>max_val){
									max_id = j;
									max_val = running_ratio[j];
								}
							}	
							for(int i = max_id; i<samplerptr->chain_N-1; i++){
								transfer_chain(samplerptr,samplerptr, i, i+1, true);	
								running_accept_ct[i] = running_accept_ct[i+1];
								running_reject_ct[i] = running_reject_ct[i+1];
								prev_accept_ct[i] = prev_accept_ct[i+1];
								prev_reject_ct[i] = prev_reject_ct[i+1];
								old_temps[i] = old_temps[i+1];
							}
							samplerptr->chain_N-=1;
							//std::cout<<"MEM CHECK : leaving rm chain"<<std::endl;
						}
						else{
							chain_pop_target_reached = true;

						}
					}
				
				}

			}
		}
		//std::cout<<"MEM CHECK : updating averages"<<std::endl;
		//Calculate average percent change in temperature
		double sum = 0;
		for (int j =0; j<samplerptr->chain_N; j++){
			sum += std::abs((samplerptr->chain_temps[j] - old_temps[j])/old_temps[j]);
			old_temps[j]=samplerptr->chain_temps[j];
		}
		ave_dynamics = sum / samplerptr->chain_N;
		if(show_prog){
			printProgress(  abs(ave_dynamics - tolerance)/(tolerance+ave_dynamics));
		}
	}
	//std::cout<<"MEM CHECK : loop finished"<<std::endl;
	double acc, rej;
	double average_accept[ensemble_size];
	double average_reject[ensemble_size];
	int swap_attemps[ensemble_size];
	int ensemble_number=0;
	for(int j = 0 ; j<ensemble_size ; j++){
		average_accept[j] = 0;
		average_reject[j] = 0;
		swap_attemps[j] = 0;
	}
	for (int j =0; j<samplerptr->chain_N; j++){
		//std::cout<<"TEMP "<<j<<": "<<samplerptr->chain_temps[j]<<std::endl;
		acc = samplerptr->swap_accept_ct[j];	
		rej = samplerptr->swap_reject_ct[j];	
		average_accept[j % ensemble_size] += (double)acc/(acc+rej);
		average_reject[j % ensemble_size] += (double)rej/(acc+rej);
		swap_attemps[j%ensemble_size] += acc+rej;
		if(j%ensemble_size == 0){ensemble_number++;}
		//std::cout<<"Accept ratio "<<j<<": "<<(double)acc/(acc+rej)<<std::endl;
		//std::cout<<"Swap attempts "<<j<<": "<<(acc+rej)<<std::endl;
		
	}
	for (int j =0; j<ensemble_size; j++){
		std::cout<<"TEMP "<<j<<": "<<samplerptr->chain_temps[j]<<std::endl;
		acc = average_accept[j]/ensemble_number;	
		rej = average_reject[j]/ensemble_number;	
		std::cout<<"Accept ratio "<<j<<": "<<acc<<std::endl;
		std::cout<<"Swap Attempts "<<j<<": "<<swap_attemps[j]<<std::endl;
	}
	delete [] running_accept_ct;
	delete [] running_reject_ct;
	delete [] prev_accept_ct;
	delete [] prev_reject_ct;
	delete [] running_ratio;
	delete [] old_temps;


}

//######################################################################################
//######################################################################################
/*!\brief Generic sampler, where the likelihood, prior are parameters supplied by the user
 *
 * Base of the sampler, generic, with user supplied quantities for most of the samplers
 * properties
 * 	
 * Uses the Metropolis-Hastings method, with the option for Fisher/MALA steps if the Fisher
 * routine is supplied.
 *
 * 3 modes to use - 
 *
 * single threaded (numThreads = 1) runs single threaded
 *
 * multi-threaded ``deterministic'' (numThreads>1 ; pool = false) progresses each chain in parallel for swp_freq steps, then waits for all threads to complete before swapping temperatures in sequenctial order (j, j+1) then (j+1, j+2) etc (sequenctially)
 *
 * multi-threaded ``stochastic'' (numThreads>2 ; pool = true) progresses each chain in parallel by queueing each temperature and evaluating them in the order they were submitted. Once finished, the threads are queued to swap, where they swapped in the order they are submitted. This means the chains are swapped randomly, and the chains do NOT finish at the same time. The sampler runs until the the 0th chain reaches the step number
 *
 * Note on limits: In the prior function, if a set of parameters should be disallowed, return -std::numeric_limits<double>::infinity()  -- (this is in the <limits> file in std)
 *
 * Format for the auto_corr file (compatable with csv, dat, txt extensions): each row is a dimension of the cold chain, with the first row being the lengths used for the auto-corr calculation:
 *
 * lengths: length1 , length2 ...
 *
 * dim1: length1 , length2 ...
 *
 * .
 *
 * .
 *
 * .
 *
 *
 * Format for the chain file (compatable with csv, dat, txt extensions): each row is a step, each column a dimension:
 *
 * Step1: dim1 , dim2 , ...
 *
 * Step2: dim1 , dim2 , ...
 *
 * .
 *
 * .
 *
 * .
 *
 * Statistics_filename : should be txt extension
 *
 * checkpoint_file : This file saves the final position of all the chains, as well as other metadata, and can be loaded by the function <FUNCTION> to continue the chain from the point it left off. Not meant to be read by humans, the data order is custom to this software library. An empty string ("") means no checkpoint will be saved. For developers, the contents are:
 *
 * dimension, # of chains
 *
 * temps of chains
 *
 * Stepping widths of all chains
 *
 * Final position of all chains
 */
void PTMCMC_MH_internal(	double ***output, /**< [out] Output chains, shape is double[chain_N, N_steps,dimension]*/
	int dimension, 	/**< dimension of the parameter space being explored*/
	int N_steps,	/**< Number of total steps to be taken, per chain*/
	int chain_N,	/**< Number of chains*/
	double *initial_pos, 	/**<Initial position in parameter space - shape double[dimension]*/
	double *seeding_var, 	/**<Variance of the normal distribution used to seed each chain higher than 0 - shape double[dimension]*/
	double **ensemble_initial_pos,
	double *chain_temps,	/**<Double array of temperatures for the chains*/
	int swp_freq,	/**< the frequency with which chains are swapped*/
	std::function<double(double*,int *,int,mcmc_data_interface *,void *)> log_prior,/**< std::function for the log_prior function -- takes double *position, int dimension, int chain_id*/
	std::function<double(double*,int*,int,mcmc_data_interface *,void *)> log_likelihood,/**< std::function for the log_likelihood function -- takes double *position, int dimension, int chain_id*/
	std::function<void(double*,int*,int,double**,mcmc_data_interface *,void *)>fisher,/**< std::function for the fisher function -- takes double *position, int dimension, double **output_fisher, int chain_id*/
	void **user_parameters,/**< Void pointer to any parameters the user may need inside log_prior, log_likelihood, or fisher. Should have one pointer for each chain. If this isn't needed, use (void**) NULL**/
	int numThreads, /**< Number of threads to use (=1 is single threaded)*/
	bool pool, /**< boolean to use stochastic chain swapping (MUST have >2 threads)*/
	bool show_prog, /**< boolean whether to print out progress (for example, should be set to ``false'' if submitting to a cluster)*/
	std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
	std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output*/
	std::string checkpoint_file,/**< Filename to output data for checkpoint, if empty string, not saved*/
	bool tune, 
	bool burn_phase
	)
{
	int initial_model_status=0;
	clock_t start, end, acend;
	double wstart, wend, wacend;
	start = clock();
	wstart = omp_get_wtime();

	sampler sampler_ref;
	sampler *samplerptr = &sampler_ref;
	samplerptr_global = &sampler_ref;
	
	samplerptr->tune = tune;
	samplerptr->burn_phase = burn_phase;
	if(burn_phase){
		samplerptr->isolate_ensembles=false;	
		samplerptr->isolate_ensembles_cold = false;
	}


	//if Fisher is not provided, Fisher and MALA steps
	//aren't used
	if(fisher ==NULL){
		samplerptr->fisher_exist = false;
	}
	else 
		samplerptr->fisher_exist = true;
	
	samplerptr->log_ll = true;
	samplerptr->log_lp = true;
	
	//Construct sampler structure
	samplerptr->lp = log_prior;
	samplerptr->ll = log_likelihood;
	samplerptr->fish = fisher;
	samplerptr->swp_freq = 2;
	samplerptr->swap_rate = 1./swp_freq;
	samplerptr->chain_temps = chain_temps;
	samplerptr->chain_N = chain_N;
	samplerptr->N_steps = N_steps;
	samplerptr->dimension =samplerptr->min_dim =samplerptr->max_dim = dimension;
	samplerptr->show_progress = show_prog;
	samplerptr->num_threads = numThreads;
	samplerptr->user_parameters=user_parameters;


	samplerptr->output = output;
	samplerptr->pool = pool;
	samplerptr->numThreads = numThreads;

	allocate_sampler_mem(samplerptr);

	//Set chains with temp 1 to highest priority
	if(samplerptr->prioritize_cold_chains){
		for(int i =0 ;i<samplerptr->chain_N; i++){
			if(samplerptr->chain_temps[i]==1)
				samplerptr->priority[i] = 0;
		}
	}

	for (int chain_index=0; chain_index<samplerptr->chain_N; chain_index++)
		assign_probabilities(samplerptr, chain_index);
	

	int  k=0;
	//int swp_accepted=0, swp_rejected=0;
	int *step_accepted = (int *)malloc(sizeof(int)*samplerptr->chain_N);
	int *step_rejected = (int *)malloc(sizeof(int)*samplerptr->chain_N);
	for(int j=0; j<samplerptr->chain_N; j++){
		step_accepted[j]=0;
		step_rejected[j]=0;
	}
	
	//Assign initial position to start chains
	int *init_status = new int[samplerptr->dimension];
	for(int i =0 ; i<samplerptr->dimension; i++)init_status[i]=1;
	int *ensemble_initial_model_status = NULL;
	int **ensemble_initial_status = new int*[samplerptr->chain_N];
	for(int j = 0 ; j<samplerptr->chain_N; j++){
		ensemble_initial_status[j] = new int[samplerptr->max_dim];
		for(int i = 0; i<samplerptr->max_dim; i++){
			ensemble_initial_status[j][i]=1;	
		}
	}
	assign_initial_pos(samplerptr, initial_pos, init_status,initial_model_status,ensemble_initial_pos,ensemble_initial_status,ensemble_initial_model_status,seeding_var);	
	delete [] init_status;
	for(int j = 0 ; j<samplerptr->chain_N; j++){
		delete [] ensemble_initial_status[j];
	}
	delete [] ensemble_initial_status;

		
	//############################################################
	//##############################################################
	PTMCMC_MH_loop(samplerptr);	
	//############################################################
	//##############################################################
	int swp_accepted=0, swp_rejected=0;
	for (int i =0;i<samplerptr->chain_N; i++)
	{
		swp_accepted+=samplerptr->swap_accept_ct[i];
		swp_rejected+=samplerptr->swap_reject_ct[i];
		step_accepted[i]+=samplerptr->step_accept_ct[i];
		step_rejected[i]+=samplerptr->step_reject_ct[i];
	}
	end =clock();
	wend =omp_get_wtime();

	samplerptr->time_elapsed_cpu = (double)(end-start)/CLOCKS_PER_SEC;
	samplerptr->time_elapsed_wall = (double)(wend-wstart);

	if(show_prog){
		std::cout<<std::endl;
	}
	
	acend =clock();
	wacend =omp_get_wtime();
	samplerptr->time_elapsed_cpu_ac = (double)(acend-end)/CLOCKS_PER_SEC;
	samplerptr->time_elapsed_wall_ac = (double)(wacend - wend);

	std::cout<<std::endl;
	double accepted_percent = (double)(swp_accepted)/(swp_accepted+swp_rejected);
	double rejected_percent = (double)(swp_rejected)/(swp_accepted+swp_rejected);
	std::cout<<"Accepted percentage of chain swaps (all chains): "<<accepted_percent<<std::endl;
	std::cout<<"Rejected percentage of chain swaps (all chains): "<<rejected_percent<<std::endl;
	accepted_percent = (double)(step_accepted[0])/(step_accepted[0]+step_rejected[0]);
	rejected_percent = (double)(step_rejected[0])/(step_accepted[0]+step_rejected[0]);
	std::cout<<"Accepted percentage of steps (cold chain): "<<accepted_percent<<std::endl;
	std::cout<<"Rejected percentage of steps (cold chain): "<<rejected_percent<<std::endl;
	double nansum=0;
	for (int i =0; i< chain_N; i++)
		nansum+= samplerptr->nan_counter[i];
	std::cout<<"NANS in Fisher Calculations (all chains): "<<nansum<<std::endl;
	
	if(statistics_filename != "")
		write_stat_file(samplerptr, statistics_filename);
	
	if(chain_filename != "")
		write_output_file(chain_filename, samplerptr->N_steps, samplerptr->max_dim, samplerptr->output, samplerptr->param_status,samplerptr->model_status,samplerptr->chain_N,samplerptr->nested_model_number,samplerptr->chain_temps,false);

	if(checkpoint_file !=""){
		write_checkpoint_file(samplerptr, checkpoint_file);
	}

	free(step_accepted);
	free(step_rejected);
	deallocate_sampler_mem(samplerptr);
}

//######################################################################################
//######################################################################################
/*! \brief Routine to take a checkpoint file and begin a new chain at said checkpoint
 *
 * See MCMC_MH_internal for more details of parameters (pretty much all the same)
 */
void continue_PTMCMC_MH_internal(sampler *sampler_in,
	std::string start_checkpoint_file,/**< File for starting checkpoint*/
	double ***output,/**< [out] output array, dimensions: output[chain_N][N_steps][dimension]*/
	int N_steps,/**< Number of new steps to take*/
	int swp_freq,/**< frequency of swap attempts between temperatures*/
	std::function<double(double*,int*,int,mcmc_data_interface *,void *)> log_prior,/**< std::function for the log_prior function -- takes double *position, int dimension, int chain_id*/
	std::function<double(double*,int*,int,mcmc_data_interface *,void *)> log_likelihood,/**< std::function for the log_likelihood function -- takes double *position, int dimension, int chain_id*/
	std::function<void(double*,int*,int,double**,mcmc_data_interface *,void *)>fisher,/**< std::function for the fisher function -- takes double *position, int dimension, double **output_fisher, int chain_id*/
	void **user_parameters,/**< Void pointer to any parameters the user may need inside log_prior, log_likelihood, or fisher. Should have one pointer for each chain. If this isn't needed, use (void**) NULL**/
	int numThreads,/**<Number of threads to use*/
	bool pool,/**<Boolean for whether to use ``deterministic'' vs ``stochastic'' sampling*/
	bool show_prog,/**< Boolean for whether to show progress or not (turn off for cluster runs*/
	std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
	std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output*/
	std::string end_checkpoint_file,/**< Filename to output data for checkpoint at the end of the continued run, if empty string, not saved*/
	bool tune, /**<Allow for dynamic tuning -- technically should be turned off when drawing final samples*/
	bool burn_phase /**<Passed to LL function for user information*/
	)
{
	//if(pool){
	//	std::cout<<"Pooling"<<std::endl;
	//}
	clock_t start, end, acend;
	double wstart, wend, wacend;
	start = clock();
	wstart = omp_get_wtime();

	//sampler sampler;
	sampler *samplerptr = sampler_in;
	samplerptr_global = sampler_in;
	
	samplerptr->tune=tune;
	samplerptr->burn_phase = burn_phase;
	if(burn_phase){
		samplerptr->isolate_ensembles=false;	
		samplerptr->isolate_ensembles_cold = false;
	}

	//if Fisher is not provided, Fisher and MALA steps
	//aren't used
	if(fisher ==NULL){
		samplerptr->fisher_exist = false;
		//samplerptr->fisher_exist = true;
		//fisher = [](double *param, int*param_status,  double **fish, mcmc_data_interface *interface, void *parameters){
		//	return fisher_generic(param,param_status,fish, interface, parameters, samplerptr);	
		//};
	}
	else 
		samplerptr->fisher_exist = true;

	samplerptr->log_ll = true;
	samplerptr->log_lp = true;
	
	//Construct sampler structure
	samplerptr->lp = log_prior;
	samplerptr->ll = log_likelihood;
	samplerptr->fish = fisher;
	samplerptr->swp_freq = 2.;
	samplerptr->swap_rate = 1./swp_freq;
	samplerptr->N_steps = N_steps;
	samplerptr->show_progress = show_prog;
	samplerptr->num_threads = numThreads;
	samplerptr->user_parameters=user_parameters;


	samplerptr->output = output;
	samplerptr->pool = pool;

	//Unpack checkpoint file -- allocates memory internally -- separate call unneccessary
	load_checkpoint_file(start_checkpoint_file, samplerptr);

	//allocate other parameters
	for (int chain_index=0; chain_index<samplerptr->chain_N; chain_index++)
		assign_probabilities(samplerptr, chain_index);
	for (int j=0;j<samplerptr->chain_N;j++){
		samplerptr->current_likelihoods[j] =
			 samplerptr->ll(samplerptr->output[j][0],samplerptr->param_status[j][0],samplerptr->model_status[j][0],samplerptr->interfaces[j],samplerptr->user_parameters[j]);
	}
	if(samplerptr->log_ll && samplerptr->log_lp){
		for(int i = 0 ; i<samplerptr->chain_N; i++){
			samplerptr->ll_lp_output[i][0][0] = samplerptr->current_likelihoods[i];
			samplerptr->ll_lp_output[i][0][1] = samplerptr->lp(samplerptr->output[i][0],samplerptr->param_status[i][0],samplerptr->model_status[i][0],samplerptr->interfaces[i],samplerptr->user_parameters[i]);
		}
	}
	
	//Set chains with temp 1 to highest priority
	if(samplerptr->prioritize_cold_chains){
		for(int i =0 ;i<samplerptr->chain_N; i++){
			if(fabs(samplerptr->chain_temps[i]-1)<DOUBLE_COMP_THRESH)
				samplerptr->priority[i] = 0;
		}
	}
	
	//########################################################
	//########################################################
	PTMCMC_MH_loop(samplerptr);	
	//########################################################
	//########################################################
	
	end =clock();
	wend =omp_get_wtime();

	samplerptr->time_elapsed_cpu = (double)(end-start)/CLOCKS_PER_SEC;
	samplerptr->time_elapsed_wall = (double)(wend-wstart);

	if(show_prog){
		std::cout<<std::endl;
	}
	
	acend =clock();
	wacend =omp_get_wtime();
	samplerptr->time_elapsed_cpu_ac = (double)(acend-end)/CLOCKS_PER_SEC;
	samplerptr->time_elapsed_wall_ac = (double)(wacend - wend);

	
	if(statistics_filename != "")
		write_stat_file(samplerptr, statistics_filename);
	
	if(chain_filename != "")
		write_output_file(chain_filename, samplerptr->N_steps, samplerptr->max_dim, samplerptr->output, samplerptr->param_status,samplerptr->model_status,samplerptr->chain_N,samplerptr->nested_model_number,samplerptr->chain_temps,false);

	if(end_checkpoint_file !=""){
		write_checkpoint_file(samplerptr, end_checkpoint_file);
	}

	//temps usually allocated by user, but for continued chains, this is done internally
	free(samplerptr->chain_temps);
}
					
//######################################################################################
//######################################################################################
/*!\brief Internal function that runs the actual loop for the sampler -- increment version
 *
 * The regular loop function runs for the entire range, this increment version will only step ``increment'' steps -- asynchronous: steps are measured by the cold chains
 *
 */
void PTMCMC_MH_step_incremental(sampler *samplerptr, int increment)
{
	//Make sure we're not going out of memory
	bool reset_ref_status=true;
	if (samplerptr->progress + increment > samplerptr->N_steps){
		increment = samplerptr->N_steps - samplerptr->progress;
		reset_ref_status = false;
	}
	//Sampler Loop - ``Deterministic'' swapping between chains
	if (!samplerptr->pool)
	{
		int k =0;
		int step_log;
		omp_set_num_threads(samplerptr->num_threads);
		#pragma omp parallel ADOLC_OPENMP
		//#pragma omp parallel 
		{
		while (k<(increment-1) ){
			#pragma omp for 
			for (int j=0; j<samplerptr->chain_N; j++)
			{
				int cutoff ;
				if( samplerptr->N_steps-samplerptr->chain_pos[j] <= samplerptr->swp_freq) 
					cutoff = samplerptr->N_steps-samplerptr->chain_pos[j]-1;	
				else cutoff = samplerptr->swp_freq;	
				if(j==0)
					step_log = cutoff;
				for (int i = 0 ; i< cutoff;i++)
				{
					int success;
					//if(!samplerptr->RJMCMC){
					//	success = mcmc_step(samplerptr, samplerptr->output[j][samplerptr->chain_pos[j]], samplerptr->output[j][samplerptr->chain_pos[j]+1],samplerptr->param_status[j][0],samplerptr->param_status[j][0],j);	
					//}
					//else
					{
						success = mcmc_step(samplerptr, samplerptr->output[j][samplerptr->chain_pos[j]], samplerptr->output[j][samplerptr->chain_pos[j]+1],samplerptr->param_status[j][samplerptr->chain_pos[j]],samplerptr->param_status[j][samplerptr->chain_pos[j]+1],&(samplerptr->model_status[j][samplerptr->chain_pos[j]]),&(samplerptr->model_status[j][samplerptr->chain_pos[j]+1]),j);	
					}
					samplerptr->chain_pos[j]+=1;
					if(success==1){samplerptr->step_accept_ct[j]+=1;}
					else{samplerptr->step_reject_ct[j]+=1;}
					if(!samplerptr->de_primed[j])
						update_history(samplerptr,samplerptr->output[j][samplerptr->chain_pos[j]],samplerptr->param_status[j][samplerptr->chain_pos[j]], j);
					else if(samplerptr->chain_pos[j]%samplerptr->history_update==0)
						update_history(samplerptr,samplerptr->output[j][samplerptr->chain_pos[j]], samplerptr->param_status[j][samplerptr->chain_pos[j]],j);
					//Log LogLikelihood and LogPrior	
					if(samplerptr->log_ll){
						samplerptr->ll_lp_output[j][samplerptr->chain_pos[j]][0]= 
							samplerptr->current_likelihoods[j];
					}
					if(samplerptr->log_lp){
						samplerptr->ll_lp_output[j][samplerptr->chain_pos[j]][1]= 
							samplerptr->lp(
							samplerptr->output[j][samplerptr->chain_pos[j]],
							samplerptr->param_status[j][samplerptr->chain_pos[j]],
							samplerptr->model_status[j][samplerptr->chain_pos[j]],
							samplerptr->interfaces[j],samplerptr->user_parameters[j]);
					}
					//Update step-widths to optimize acceptance ratio
					update_step_widths(samplerptr, j);
					if(!samplerptr->de_primed[j]) 
					{
						if ((samplerptr->chain_pos[j])>samplerptr->history_length)
						{
							samplerptr->de_primed[j]=true;
							assign_probabilities(samplerptr,j);	
						}
					}
					
				}
				//if(!samplerptr->de_primed[j]) 
				//{
				//	if ((samplerptr->chain_pos[j])>samplerptr->history_length)
				//	{
				//		samplerptr->de_primed[j]=true;
				//		assign_probabilities(samplerptr,j);	
				//	}
				//}
			}
			#pragma omp single
			{
				k+= step_log;
				samplerptr->progress+=step_log;
				int swp_accepted=0, swp_rejected=0;
				//TEST
				double alpha = gsl_rng_uniform(samplerptr->rvec[0]);
				if(alpha< samplerptr->swap_rate){
					if(samplerptr->random_swaps){
						full_random_swap(samplerptr);
					}
					else{
						chain_swap(samplerptr, samplerptr->output, samplerptr->param_status,samplerptr->model_status,k, &swp_accepted, &swp_rejected);
					}
				}
				if(samplerptr->show_progress)
					printProgress((double)samplerptr->progress/samplerptr->N_steps);	
			}
		}
		}
	}

	//POOLING  -- ``Stochastic'' swapping between chains
	else
	{
		int max_steps = increment+samplerptr->progress-1;
		ThreadPool pool(samplerptr->num_threads);
		poolptr = &pool;
		//while(samplerptr->progress<increment-1)
		while(!check_sampler_status(samplerptr))
		{
			for(int i =0; i<samplerptr->chain_N; i++)
			{
				bool waiting,waiting_SWP;
				{
					std::unique_lock<std::mutex> lock{samplerptr->queue_mutexes[i]};
					waiting = samplerptr->waiting[i];
					waiting_SWP = samplerptr->waiting_SWP[i];
				}
				//if(samplerptr->waiting[i]){
				if(waiting){
					if(samplerptr->chain_pos[i]<(max_steps))
					{
						samplerptr->waiting[i]=false;
						poolptr->enqueue(i);
					}
					//If a chain finishes before chain 0, it's wrapped around 
					//and allowed to keep stepping at low priority-- 
					//not sure if this is the best
					//method for keeping the 0th chain from finishing last or not
					//TESTING NEW METHOD
					else if( ! check_list(i,samplerptr->ref_chain_ids,samplerptr->ref_chain_num) ){
					//else if(fabs(samplerptr->chain_temps[i] -1)>DOUBLE_COMP_THRESH){

						samplerptr->waiting[i]=false;
						//TESTING
						//samplerptr->priority[i] = 2;
						int pos = samplerptr->chain_pos[i];
						for (int k =0; k<samplerptr->dimension; k++){
							samplerptr->output[i][0][k] = 
								samplerptr->output[i][pos][k] ;
							samplerptr->param_status[i][0][k] = 
								samplerptr->param_status[i][pos][k] ;
						}
						samplerptr->chain_pos[i] = 0;
						samplerptr->restarted_chain[i] = true;

						poolptr->enqueue(i);
					}
					else{
						//If 0 T chain, just wait till everything else is done
						samplerptr->waiting[i] = false;	
						samplerptr->ref_chain_status[i] = true;
					}
				}
				else if(waiting_SWP){
					samplerptr->waiting_SWP[i] = false;
					poolptr->enqueue_swap(i);

				}
				
			}
			if(samplerptr->show_progress)
				printProgress((double)samplerptr->progress/samplerptr->N_steps);
			while(poolptr->queue_swap_length() > samplerptr->num_threads 
				||  poolptr->queue_length() > samplerptr->num_threads){
				usleep(QUEUE_SLEEP);
			}
		}
		for(int i= 0 ; i<samplerptr->chain_N; i++){
			if(check_list(i,samplerptr->ref_chain_ids,samplerptr->ref_chain_num)){
				samplerptr->ref_chain_status[i]=false;
			}
			else{
				samplerptr->priority[i]=1;	
			}
			samplerptr->waiting[i]=true;	
		}
	}
}
//######################################################################################
//######################################################################################
/*!\brief Internal function that runs the actual loop for the sampler
 *
 */
void PTMCMC_MH_loop(sampler *samplerptr)
{
	int k =0;
	int cutoff ;
	//Sampler Loop - ``Deterministic'' swapping between chains
	if (!samplerptr->pool)
	{
		omp_set_num_threads(samplerptr->num_threads);
		//#pragma omp parallel 
		#pragma omp parallel ADOLC_OPENMP
		{
		while (k<samplerptr->N_steps-1){
			#pragma omp single
			{
				if( samplerptr->N_steps-k <= samplerptr->swp_freq) 
					cutoff = samplerptr->N_steps-k-1;	
				else cutoff = samplerptr->swp_freq;	
			}
			#pragma omp for 
			for (int j=0; j<samplerptr->chain_N; j++)
			{
				for (int i = 0 ; i< cutoff;i++)
				{
					int success;
					//if(!sampler->RJMCMC){
					//	success = mcmc_step(sampler, sampler->output[j][k+i], sampler->output[j][k+i+1],sampler->param_status[j][0],sampler->param_status[j][0],j);	
					//}
					//else
					{
						success = mcmc_step(samplerptr, samplerptr->output[j][k+i], samplerptr->output[j][k+i+1],samplerptr->param_status[j][k+i],samplerptr->param_status[j][k+i+1],&(samplerptr->model_status[j][k+i]),&(samplerptr->model_status[j][k+i+1]),j);	
					}
					samplerptr->chain_pos[j]+=1;
					if(success==1){samplerptr->step_accept_ct[j]+=1;}
					else{samplerptr->step_reject_ct[j]+=1;}
					if(!samplerptr->de_primed[j])
						update_history(samplerptr,samplerptr->output[j][k+i+1], samplerptr->param_status[j][k+i+1],j);
					//else if(samplerptr->chain_pos[j]%samplerptr->history_update==0)
					else if((k+i+1)%samplerptr->history_update==0)
						update_history(samplerptr,samplerptr->output[j][k+i+1],samplerptr->param_status[j][k+i+1], j);
					if(!samplerptr->de_primed[j]) 
					{
						if ((k+i+1)>samplerptr->history_length)
						{
							samplerptr->de_primed[j]=true;
							assign_probabilities(samplerptr,j);	
						}
					}
					//Log LogLikelihood and LogPrior	
					if(samplerptr->log_ll){
						samplerptr->ll_lp_output[j][k+i+1][0]= 
							samplerptr->current_likelihoods[j];
					}
					if(samplerptr->log_lp){
						samplerptr->ll_lp_output[j][k+i+1][1]= 
							samplerptr->lp(
							samplerptr->output[j][k+i+1],
							samplerptr->param_status[j][k+i+1],
							samplerptr->model_status[j][k+i+1],
							samplerptr->interfaces[j],samplerptr->user_parameters[j]);
					}
					//Update step-widths to optimize acceptance ratio
					update_step_widths(samplerptr, j);
					
				}
				//if(!samplerptr->de_primed[j]) 
				//{
				//	if ((k+cutoff)>samplerptr->history_length)
				//	{
				//		samplerptr->de_primed[j]=true;
				//		assign_probabilities(samplerptr,j);	
				//	}
				//}
			}
			#pragma omp single
			{
				k+= cutoff;
				int swp_accepted=0, swp_rejected=0;
				double alpha = gsl_rng_uniform(samplerptr->rvec[0]);
				if(alpha< samplerptr->swap_rate){
					if(samplerptr->random_swaps){
						full_random_swap(samplerptr);
					}
					else{
						chain_swap(samplerptr, samplerptr->output, samplerptr->param_status,samplerptr->model_status,k, &swp_accepted, &swp_rejected);
					}
				}

				if(samplerptr->show_progress)
					printProgress((double)k/samplerptr->N_steps);	
			}
		}
		}
	}

	//POOLING  -- ``Stochastic'' swapping between chains
	else
	{
		ThreadPool pool(samplerptr->num_threads);
		poolptr = &pool;
		//while(samplerptr->progress<samplerptr->N_steps-1)
		while(!check_sampler_status(samplerptr))
		{
			for(int i =0; i<samplerptr->chain_N; i++)
			{
				bool waiting,waiting_SWP;
				{
					std::unique_lock<std::mutex> lock{samplerptr->queue_mutexes[i]};
					waiting = samplerptr->waiting[i];
					waiting_SWP = samplerptr->waiting_SWP[i];
				}
				//if(samplerptr->waiting[i]){
				if(waiting){
					if(samplerptr->chain_pos[i]<(samplerptr->N_steps-1))
					{
						samplerptr->waiting[i]=false;
						//if(i==0) samplerptrptr->progress+=samplerptrptr->swp_freq;
						poolptr->enqueue(i);
					}
					//If a chain finishes before chain 0, it's wrapped around 
					//and allowed to keep stepping at low priority-- 
					//not sure if this is the best
					//method for keeping the 0th chain from finishing last or not
					//
					//TESTING
					else if( ! check_list(i,samplerptr->ref_chain_ids,samplerptr->ref_chain_num) ){
					//else if(fabs(samplerptr->chain_temps[i] -1)>DOUBLE_COMP_THRESH){
					//We're gonna try something here
					//else if(false){

						samplerptr->waiting[i]=false;
						//std::cout<<"Chain "<<i<<" finished-- being reset"<<std::endl;
						//TESTING
						//samplerptr->priority[i] = 2;



						//TESTING
						int pos = samplerptr->chain_pos[i];
						for (int k =0; k<samplerptr->max_dim; k++){
							samplerptr->output[i][0][k] = 
								samplerptr->output[i][pos][k] ;
							samplerptr->param_status[i][0][k] = 
								samplerptr->param_status[i][pos][k] ;
						}
						samplerptr->chain_pos[i] = 0;
						samplerptr->restarted_chain[i] = true;



						//TESTING
						//samplerptr->chain_pos[i] = samplerptr->N_steps-2;





						poolptr->enqueue(i);


					}
					else{
						//If 1 T chain, just wait till everything else is done
						samplerptr->waiting[i] = false;	
						samplerptr->ref_chain_status[i] = true;	
						//Still trying something
						//if(fabs(samplerptr->chain_temps[i] -1)<DOUBLE_COMP_THRESH){
						//	samplerptr->ref_chain_status[i] = true;	
						//}
					}
				}
				else if(waiting_SWP){
					samplerptr->waiting_SWP[i] = false;
					poolptr->enqueue_swap(i);

				}
				
			}
			if(samplerptr->show_progress)
				printProgress((double)samplerptr->progress/samplerptr->N_steps);
			//usleep(10);
			//We get in trouble when over accessing the queue -- just wait until more jobs can be run
			//This may be more inefficient short term, but it will ensure a more balanced load for 
			//the ensemble
			while(poolptr->queue_swap_length() > samplerptr->num_threads 
				||  poolptr->queue_length() > samplerptr->num_threads){
				usleep(QUEUE_SLEEP);
			}
		}
	}
}

void full_random_swap(sampler *s)
{
	//debugger_print(__FILE__,__LINE__,"Started swap");
	ThreadPool pool(s->num_threads);
	poolptr = &pool;
	std::vector<int> ids;
	int ct = 0 ;
	for(int i =0 ; i<s->chain_N; i++){
		ids.insert(ids.begin()+i,i);
		ct++;
	}
	int alpha;
	while(ids.size()>0){
		alpha=(int)(gsl_rng_uniform(s->rvec[0])*ids.size());
		poolptr->enqueue_swap(ids.at(alpha));	
		ids.erase(ids.begin()+alpha);
	}
	
	//poolptr->flush_swap_queue();
	int pairs = poolptr->queue_pairs_length();
	while(pairs > 0  ){
		pairs = poolptr->queue_pairs_length();
		usleep(1e2);
	}
	//debugger_print(__FILE__,__LINE__,"Finished swap");
}

//######################################################################################
//######################################################################################
void mcmc_step_threaded(int j)
{
	int k = samplerptr_global->chain_pos[j];
	int cutoff;
	if( samplerptr_global->N_steps-k <= samplerptr_global->swp_freq) cutoff = samplerptr_global->N_steps-k-1;	
	else cutoff = samplerptr_global->swp_freq;	
	for (int i = 0 ; i< cutoff;i++)
	{
		int success;
		//success = mcmc_step(samplerptr, samplerptr->output[j][k+i], samplerptr->output[j][k+i+1],j);	
		//if(!samplerptr->RJMCMC){
		//	success = mcmc_step(samplerptr, samplerptr->output[j][k+i], samplerptr->output[j][k+i+1],samplerptr->param_status[j][0],samplerptr->param_status[j][0],j);	
		//}
		//else
		{
			success = mcmc_step(samplerptr_global, samplerptr_global->output[j][k+i], samplerptr_global->output[j][k+i+1],samplerptr_global->param_status[j][k+i],samplerptr_global->param_status[j][k+i+1],&(samplerptr_global->model_status[j][k+i]),&(samplerptr_global->model_status[j][k+i+1]),j);	
		}
		samplerptr_global->chain_pos[j]+=1;
	
		if(success==1){samplerptr_global->step_accept_ct[j]+=1;}
		else{samplerptr_global->step_reject_ct[j]+=1;}
		if(!samplerptr_global->de_primed[j])
			update_history(samplerptr_global,samplerptr_global->output[j][k+i+1], samplerptr_global->param_status[j][k+i+1],j);
		//else if(samplerptr_global->chain_pos[j]%samplerptr_global->history_update==0)
		else if((k+i+1)%samplerptr_global->history_update==0)
			update_history(samplerptr_global,samplerptr_global->output[j][k+i+1], samplerptr_global->param_status[j][k+i+1],j);
		if(!samplerptr_global->de_primed[j]) 
		{
			if ((k+i+1)>samplerptr_global->history_length)
			{
				samplerptr_global->de_primed[j]=true;
				assign_probabilities(samplerptr_global,j);	
			}
		}
		//##############################################################
		if(samplerptr_global->log_ll){
			samplerptr_global->ll_lp_output[j][k+i+1][0]= 
				samplerptr_global->current_likelihoods[j];
		}
		if(samplerptr_global->log_lp){
			samplerptr_global->ll_lp_output[j][k+i+1][1]= 
				samplerptr_global->lp(samplerptr_global->output[j][k+i+1],
				samplerptr_global->param_status[j][k+i+1],
				samplerptr_global->model_status[j][k+i+1],
				samplerptr_global->interfaces[i],samplerptr_global->user_parameters[j]);
		}
		//##############################################################
	}
	//if(!samplerptr_global->de_primed[j]) 
	//{
	//	if ((k+cutoff)>samplerptr_global->history_length)
	//	{
	//		samplerptr_global->de_primed[j]=true;
	//		assign_probabilities(samplerptr_global,j);	
	//	}
	//}
	//samplerptr_global->chain_pos[j]+=cutoff;
	//Keep track of progress of all cold chains - track the slowest
	//Now that progress is only used for outputing progress bar, and not used
	//for stopping criteria, just track chain 0, which should be a T=1 chain anyway
	if(j==0) samplerptr_global->progress =samplerptr_global->chain_pos[j];

	//update stepsize to maximize step efficiency
	//increases in stepsizes of 10%
	update_step_widths(samplerptr_global, j);


	double alpha = gsl_rng_uniform(samplerptr_global->rvec[j]);
	if(alpha< samplerptr_global->swap_rate){
		std::unique_lock<std::mutex> lock{samplerptr_global->queue_mutexes[j]};
		//poolptr->enqueue_swap(j);
		samplerptr_global->waiting_SWP[j] = true;
	}
	else{
		std::unique_lock<std::mutex> lock{samplerptr_global->queue_mutexes[j]};
		samplerptr_global->waiting[j]=true;
	}
}

//######################################################################################
//######################################################################################
void mcmc_swap_threaded(int i, int j)
{
	//debugger_print(__FILE__,__LINE__,std::to_string(i)+" "+std::to_string(j)+" "+std::to_string(samplerptr->chain_temps[i])+" "+std::to_string(samplerptr->chain_temps[j]));
	int k = samplerptr_global->chain_pos[i];
	int l = samplerptr_global->chain_pos[j];
	int success;
	success = single_chain_swap(samplerptr_global, samplerptr_global->output[i][k], samplerptr_global->output[j][l],samplerptr_global->param_status[i][k],samplerptr_global->param_status[j][l],&(samplerptr_global->model_status[i][k]),&(samplerptr_global->model_status[j][l]),i,j);
	//success = -1;
	if(success ==1){
		samplerptr_global->swap_accept_ct[i]+=1;	
		samplerptr_global->swap_accept_ct[j]+=1;	
		//NOTE: this only works if the swap radius is 1
		if(samplerptr_global->PT_alloc){
			if(samplerptr_global->chain_temps[i]<samplerptr_global->chain_temps[j]){
				samplerptr_global->A[j] = 1;
			}
			else{
				samplerptr_global->A[i] = 1;
			}
		}
	}
	else{
		samplerptr_global->swap_reject_ct[i]+=1;	
		samplerptr_global->swap_reject_ct[j]+=1;	
		//NOTE: this only works if the swap radius is 1
		if(samplerptr_global->PT_alloc){
			if(samplerptr_global->chain_temps[i]<samplerptr_global->chain_temps[j]){
				samplerptr_global->A[j] = 0;
			}
			else{
				samplerptr_global->A[i] = 0;
			}
		}
	}
	{
		std::unique_lock<std::mutex> locki{samplerptr_global->queue_mutexes[i]};
		samplerptr_global->waiting[i]=true;
	}
	{
		std::unique_lock<std::mutex> lockj{samplerptr_global->queue_mutexes[j]};
		samplerptr_global->waiting[j]=true;
	}
}
void mcmc_ensemble_step(sampler *s, int id1, int id2)
{
	//for(int i = 0 ; i<s->chains_per_ensemble; i++)
	//{
	//	double rand_u = gsl_rng_uniform(s->rvec[id1+i,
	//	double z = ;
	//}
}
//######################################################################################
//######################################################################################
void PTMCMC_MH(	double ***output, /**< [out] Output chains, shape is double[chain_N, N_steps,dimension]*/
	int dimension, 	/**< dimension of the parameter space being explored*/
	int N_steps,	/**< Number of total steps to be taken, per chain*/
	int chain_N,	/**< Number of chains*/
	double *initial_pos, 	/**<Initial position in parameter space - shape double[dimension]*/
	double *seeding_var, 	/**<Variance of the normal distribution used to seed each chain higher than 0 - shape double[dimension]*/
	double **ensemble_initial_pos,
	double *chain_temps,	/**<Double array of temperatures for the chains*/
	int swp_freq,	/**< the frequency with which chains are swapped*/
	double (*log_prior)(double *param, mcmc_data_interface *interface,void * parameters),	/**<Funcion pointer for the log_prior*/
	double (*log_likelihood)(double *param, mcmc_data_interface *interface,void *parameters),	/**<Function pointer for the log_likelihood*/
	void (*fisher)(double *param,  double **fisher, mcmc_data_interface *interface,void *parameters),	/**<Function pointer for the fisher - if NULL, fisher steps are not used*/
	void **user_parameters,/**< Void pointer to any parameters the user may need inside log_prior, log_likelihood, or fisher. Should have one pointer for each chain. If this isn't needed, use (void**) NULL**/
	int numThreads, /**< Number of threads to use (=1 is single threaded)*/
	bool pool, /**< boolean to use stochastic chain swapping (MUST have >2 threads)*/
	bool show_prog, /**< boolean whether to print out progress (for example, should be set to ``false'' if submitting to a cluster)*/
	std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
	std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output*/
	std::string checkpoint_file/**< Filename to output data for checkpoint, if empty string, not saved*/
	)
{
	//std::function<double(double*,int,int)> lp = log_prior;
	//std::function<double(double*,int,int)> ll = log_likelihood;
	//std::function<void(double*,int,double**,int)>f = fisher;
	std::function<double(double*,int*,int ,mcmc_data_interface *,void *)> ll =NULL;
	ll = [&log_likelihood](double *param, int *param_status,int model_status,mcmc_data_interface *interf ,void *parameters){
			return log_likelihood(param, interf,parameters);};
	std::function<double(double*,int*,int ,mcmc_data_interface *,void *)> lp =NULL;
	lp = [&log_prior](double *param, int *param_status,int model_status,mcmc_data_interface *interf,void *parameters){
			return log_prior(param, interf,parameters);};
	std::function<void(double*,int*,int ,double**,mcmc_data_interface *,void *)> f =NULL;
	if(fisher){
		f = [&fisher](double *param, int *param_status, int model_status,double **fisherm, mcmc_data_interface *interf,void *parameters){
			fisher(param, fisherm,interf,parameters);};
	}
	PTMCMC_MH_internal(output, dimension, N_steps, chain_N, initial_pos, seeding_var,ensemble_initial_pos,chain_temps, swp_freq, 
			lp, ll, f,user_parameters, numThreads, pool, show_prog, 
			statistics_filename, chain_filename,  checkpoint_file,false,false);
}
//######################################################################################
//######################################################################################
void continue_PTMCMC_MH(std::string start_checkpoint_file,/**< File for starting checkpoint*/
	double ***output,/**< [out] output array, dimensions: output[chain_N][N_steps][dimension]*/
	int N_steps,/**< Number of new steps to take*/
	int swp_freq,/**< frequency of swap attempts between temperatures*/
	double (*log_prior)(double *param, mcmc_data_interface *interface,void *parameters),	/**<Funcion pointer for the log_prior*/
	double (*log_likelihood)(double *param, mcmc_data_interface *interface,void *parameters),	/**<Function pointer for the log_likelihood*/
	void (*fisher)(double *param,  double **fisher, mcmc_data_interface *,void *parameters),	/**<Function pointer for the fisher - if NULL, fisher steps are not used*/
	void **user_parameters,/**< Void pointer to any parameters the user may need inside log_prior, log_likelihood, or fisher. Should have one pointer for each chain. If this isn't needed, use (void**) NULL**/
	int numThreads,/**<Number of threads to use*/
	bool pool,/**<Boolean for whether to use ``deterministic'' vs ``stochastic'' sampling*/
	bool show_prog,/**< Boolean for whether to show progress or not (turn off for cluster runs*/
	std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
	std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output*/
	std::string end_checkpoint_file,/**< Filename to output data for checkpoint at the end of the continued run, if empty string, not saved*/
	bool tune
	)
{

	//std::function<double(double*,int,int)> lp = log_prior;
	//std::function<double(double*,int,int)> ll = log_likelihood;
	//std::function<void(double*,int,double**,int)>f = fisher;
	std::function<double(double*,int*,int ,mcmc_data_interface *,void *)> ll =NULL;
	ll = [&log_likelihood](double *param, int *param_status,int model_status,mcmc_data_interface *interf,void *parameters){
			return log_likelihood(param, interf,parameters);};
	std::function<double(double*,int*,int ,mcmc_data_interface *,void *)> lp =NULL;
	lp = [&log_prior](double *param, int *param_status,int model_status,mcmc_data_interface *interf,void *parameters){
			return log_prior(param, interf,parameters);};
	std::function<void(double*,int*, int , double**,mcmc_data_interface *,void *)> f =NULL;
	if(fisher){
		f = [&fisher](double *param, int * param_status, int model_status,double **fisherm, mcmc_data_interface *interf,void *parameters){
			fisher(param, fisherm,interf,parameters);};
	}
	sampler sampler;
	continue_PTMCMC_MH_internal(&sampler, 
			start_checkpoint_file,
			output,
			N_steps,
			swp_freq,
			lp,
			ll,
			f,
			user_parameters,
			numThreads,
			pool,
			show_prog,
			statistics_filename,
			chain_filename,
			end_checkpoint_file,
			tune,
			false);
	deallocate_sampler_mem(&sampler);
}
//######################################################################################
//######################################################################################
void continue_PTMCMC_MH_dynamic_PT_alloc(std::string checkpoint_file_start,
	double ***output, /**< [out] Output chains, shape is double[max_chain_N, N_steps,dimension]*/
	int N_steps,	/**< Number of total steps to be taken, per chain AFTER chain allocation*/
	int max_chain_N_thermo_ensemble,/**< Maximum number of chains to use in the thermodynamic ensemble (may use less)*/
	double *chain_temps, /**< Final chain temperatures used -- should be shape double[chain_N]*/
	int swp_freq,	/**< the frequency with which chains are swapped*/
	int t0,/**< Time constant of the decay of the chain dynamics  (~1000)*/
	int nu,/**< Initial amplitude of the dynamics (~100)*/
	std::string chain_distribution_scheme, /*How to allocate the remaining chains once equilibrium is reached*/
	double (*log_prior)(double *param, mcmc_data_interface *interface,void *parameters),	/**<Funcion pointer for the log_prior*/
	double (*log_likelihood)(double *param, mcmc_data_interface *interface,void *parameters),	/**<Function pointer for the log_likelihood*/
	void (*fisher)(double *param,  double **fisher, mcmc_data_interface *interface,void *parameters),	/**<Function pointer for the fisher - if NULL, fisher steps are not used*/
	void **user_parameters,/**< Void pointer to any parameters the user may need inside log_prior, log_likelihood, or fisher. Should have one pointer for each chain. If this isn't needed, use (void**) NULL**/
	int numThreads, /**< Number of threads to use (=1 is single threaded)*/
	bool pool, /**< boolean to use stochastic chain swapping (MUST have >2 threads)*/
	bool show_prog, /**< boolean whether to print out progress (for example, should be set to ``false'' if submitting to a cluster)*/
	std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
	std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output*/
	std::string checkpoint_file/**< Filename to output data for checkpoint, if empty string, not saved*/
	)
{
	std::function<double(double*,int*,int ,mcmc_data_interface *,void *)> ll =NULL;
	ll = [&log_likelihood](double *param, int *param_status,int model_status,mcmc_data_interface *interf,void *parameters){
			return log_likelihood(param, interf,parameters);};
	std::function<double(double*,int*,int ,mcmc_data_interface *,void *)> lp =NULL;
	lp = [&log_prior](double *param, int *param_status,int model_status,mcmc_data_interface *interf,void *parameters){
			return log_prior(param, interf,parameters);};
	std::function<void(double*,int*,int ,double**,mcmc_data_interface *interf,void *)> f =NULL;
	if(fisher){
		f = [&fisher](double *param, int *param_status,int model_status, double **fisherm, mcmc_data_interface *interf,void *parameters){
			fisher(param,  fisherm,interf,parameters);};
	}
	continue_PTMCMC_MH_dynamic_PT_alloc_internal(checkpoint_file_start, 
			output,
			N_steps,
			max_chain_N_thermo_ensemble,
			chain_temps,
			swp_freq,
			t0,
			nu,
			chain_distribution_scheme,
			lp,
			ll,
			f,
			user_parameters,
			numThreads,
			pool,
			show_prog,
			true,
			statistics_filename,
			chain_filename,
			checkpoint_file,
			false);

}
//######################################################################################
//######################################################################################
void PTMCMC_MH_dynamic_PT_alloc(double ***output, /**< [out] Output chains, shape is double[max_chain_N, N_steps,dimension]*/
	int dimension, 	/**< dimension of the parameter space being explored*/
	int N_steps,	/**< Number of total steps to be taken, per chain AFTER chain allocation*/
	int chain_N,/**< Maximum number of chains to use */
	int max_chain_N_thermo_ensemble,/**< Maximum number of chains to use in the thermodynamic ensemble (may use less)*/
	double *initial_pos, 	/**<Initial position in parameter space - shape double[dimension]*/
	double *seeding_var, 	/**<Variance of the normal distribution used to seed each chain higher than 0 - shape double[dimension]*/
	double **ensemble_initial_pos,
	double *chain_temps, /**< Final chain temperatures used -- should be shape double[chain_N]*/
	int swp_freq,	/**< the frequency with which chains are swapped*/
	int t0,/**< Time constant of the decay of the chain dynamics  (~1000)*/
	int nu,/**< Initial amplitude of the dynamics (~100)*/
	std::string chain_distribution_scheme, /*How to allocate the remaining chains once equilibrium is reached*/
	double (*log_prior)(double *param, mcmc_data_interface *interface,void *parameters),	/**<Funcion pointer for the log_prior*/
	double (*log_likelihood)(double *param,  mcmc_data_interface *interface,void *parameters),	/**<Function pointer for the log_likelihood*/
	void (*fisher)(double *param,  double **fisher, mcmc_data_interface *interface,void *parameters),	/**<Function pointer for the fisher - if NULL, fisher steps are not used*/
	void **user_parameters,/**< Void pointer to any parameters the user may need inside log_prior, log_likelihood, or fisher. Should have one pointer for each chain. If this isn't needed, use (void**) NULL**/
	int numThreads, /**< Number of threads to use (=1 is single threaded)*/
	bool pool, /**< boolean to use stochastic chain swapping (MUST have >2 threads)*/
	bool show_prog, /**< boolean whether to print out progress (for example, should be set to ``false'' if submitting to a cluster)*/
	std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
	std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output*/
	std::string checkpoint_file/**< Filename to output data for checkpoint, if empty string, not saved*/
	)
{
	//std::function<double(double*,int,int)> lp = log_prior;
	//std::function<double(double*,int,int)> ll = log_likelihood;
	//std::function<void(double*,int,double**,int)>f = fisher;
	std::function<double(double*,int*,int ,mcmc_data_interface *,void *)> ll =NULL;
	ll = [&log_likelihood](double *param, int *param_status,int model_status,mcmc_data_interface *interf,void *parameters){
			return log_likelihood(param, interf,parameters);};
	std::function<double(double*,int*,int ,mcmc_data_interface* ,void *)> lp =NULL;
	lp = [&log_prior](double *param, int *param_status,int model_status,mcmc_data_interface *interf,void *parameters){
			return log_prior(param, interf,parameters);};
	std::function<void(double*,int*,int ,double**,mcmc_data_interface *,void *)> f =NULL;
	if(fisher){
		f = [&fisher](double *param, int *param_status,int model_status, double **fisherm, mcmc_data_interface *interf, void *parameters){
			fisher(param,  fisherm,interf,parameters);};
	}
	PTMCMC_MH_dynamic_PT_alloc_internal(output,
			dimension,
			N_steps,
			chain_N,
			max_chain_N_thermo_ensemble,
			initial_pos,
			seeding_var,
			ensemble_initial_pos,
			chain_temps,
			swp_freq,
			t0,
			nu,
			chain_distribution_scheme,
			lp,
			ll,
			f,
			user_parameters,
			numThreads,
			pool,
			show_prog,
			true,
			statistics_filename,
			chain_filename,
			checkpoint_file,
			false);

}
//######################################################################################
//######################################################################################
void continue_PTMCMC_MH_dynamic_PT_alloc_uncorrelated(std::string checkpoint_file_start,
	mcmc_sampler_output *sampler_output,
	double **output, /**< [out] Output, shape is double[N_steps,dimension]*/
	int N_steps,	/**< Number of total steps to be taken, per chain AFTER chain allocation*/
	int max_chain_N_thermo_ensemble,/**< Maximum number of chains to use in the thermodynamic ensemble (may use less)*/
	double *chain_temps, /**< Final chain temperatures used -- should be shape double[chain_N]*/
	int swp_freq,	/**< the frequency with which chains are swapped*/
	int t0,/**< Time constant of the decay of the chain dynamics  (~1000)*/
	int nu,/**< Initial amplitude of the dynamics (~100)*/
	int max_chunk_size,/**<Maximum number of steps to take in a single sampler run*/
	std::string chain_distribution_scheme, /*How to allocate the remaining chains once equilibrium is reached*/
	double (*log_prior)(double *param, mcmc_data_interface *,void *parameters),	/**<Funcion pointer for the log_prior*/
	double (*log_likelihood)(double *param, mcmc_data_interface *,void *parameters),	/**<Function pointer for the log_likelihood*/
	void (*fisher)(double *param, double **fisher, mcmc_data_interface *,void *parameters),	/**<Function pointer for the fisher - if NULL, fisher steps are not used*/
	void **user_parameters,/**< Void pointer to any parameters the user may need inside log_prior, log_likelihood, or fisher. Should have one pointer for each chain. If this isn't needed, use (void**) NULL**/
	int numThreads, /**< Number of threads to use (=1 is single threaded)*/
	bool pool, /**< boolean to use stochastic chain swapping (MUST have >2 threads)*/
	bool show_prog, /**< boolean whether to print out progress (for example, should be set to ``false'' if submitting to a cluster)*/
	std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
	std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output*/
	std::string likelihood_log_filename,/**< Filename to write the log_likelihood and log_prior at each step -- use empty string to skip*/
	std::string checkpoint_file/**< Filename to output data for checkpoint, if empty string, not saved*/
	)
{
	//std::function<double(double*,int,int)> lp = log_prior;
	//std::function<double(double*,int,int)> ll = log_likelihood;
	//std::function<void(double*,int,double**,int)>f = fisher;
	std::function<double(double*,int*,int,mcmc_data_interface *,void *)> ll =NULL;
	ll = [&log_likelihood](double *param, int *param_status,int model_status,mcmc_data_interface *interf,void * parameters){
			return log_likelihood(param, interf,parameters);};
	std::function<double(double*,int*,int ,mcmc_data_interface *interf,void *)> lp =NULL;
	lp = [&log_prior](double *param, int *param_status,int model_status,mcmc_data_interface *interf,void * parameters){
			return log_prior(param, interf,parameters);};
	std::function<void(double*,int*,int ,double**,mcmc_data_interface *,void *)> f =NULL;
	if(fisher){
		f = [&fisher](double *param, int *param_status,int model_status, double **fisherm, mcmc_data_interface *interf,void * parameters){
			fisher(param,  fisherm,interf,parameters);};
	}
	continue_PTMCMC_MH_dynamic_PT_alloc_uncorrelated_internal(checkpoint_file_start,
			sampler_output,
			output,
			N_steps,
			max_chain_N_thermo_ensemble,
			chain_temps,
			swp_freq,
			t0,
			nu,
			max_chunk_size,
			chain_distribution_scheme,
			lp,
			ll,
			f,
			user_parameters,
			numThreads,
			pool,
			show_prog,
			statistics_filename,
			chain_filename,
			likelihood_log_filename,
			checkpoint_file);
	return ;
}
//######################################################################################
//######################################################################################
void PTMCMC_MH_dynamic_PT_alloc_uncorrelated(mcmc_sampler_output *sampler_output,
	double **output, /**< [out] Output, shape is double[N_steps,dimension]*/
	int dimension, 	/**< dimension of the parameter space being explored*/
	int N_steps,	/**< Number of total steps to be taken, per chain AFTER chain allocation*/
	int chain_N,/**< Maximum number of chains to use */
	int max_chain_N_thermo_ensemble,/**< Maximum number of chains to use in the thermodynamic ensemble (may use less)*/
	double *initial_pos, 	/**<Initial position in parameter space - shape double[dimension]*/
	double *seeding_var, 	/**<Variance of the normal distribution used to seed each chain higher than 0 - shape double[dimension]*/
	double **ensemble_initial_pos,
	double *chain_temps, /**< Final chain temperatures used -- should be shape double[chain_N]*/
	int swp_freq,	/**< the frequency with which chains are swapped*/
	int t0,/**< Time constant of the decay of the chain dynamics  (~1000)*/
	int nu,/**< Initial amplitude of the dynamics (~100)*/
	int max_chunk_size,/**<Maximum number of steps to take in a single sampler run*/
	std::string chain_distribution_scheme, /*How to allocate the remaining chains once equilibrium is reached*/
	double (*log_prior)(double *param, mcmc_data_interface *interface,void *parameters),	/**<Funcion pointer for the log_prior*/
	double (*log_likelihood)(double *param, mcmc_data_interface *interface,void *parameters),	/**<Function pointer for the log_likelihood*/
	void (*fisher)(double *param,  double **fisher, mcmc_data_interface *interface,void *parameters),	/**<Function pointer for the fisher - if NULL, fisher steps are not used*/
	void **user_parameters,/**< Void pointer to any parameters the user may need inside log_prior, log_likelihood, or fisher. Should have one pointer for each chain. If this isn't needed, use (void**) NULL**/
	int numThreads, /**< Number of threads to use (=1 is single threaded)*/
	bool pool, /**< boolean to use stochastic chain swapping (MUST have >2 threads)*/
	bool show_prog, /**< boolean whether to print out progress (for example, should be set to ``false'' if submitting to a cluster)*/
	std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
	std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output*/
	std::string likelihood_log_filename,/**< Filename to write the log_likelihood and log_prior at each step -- use empty string to skip*/
	std::string checkpoint_file/**< Filename to output data for checkpoint, if empty string, not saved*/
	)
{
	//std::function<double(double*,int,int)> lp = log_prior;
	//std::function<double(double*,int,int)> ll = log_likelihood;
	//std::function<void(double*,int,double**,int)>f = fisher;
	std::function<double(double*,int*,int,mcmc_data_interface *,void *)> ll =NULL;
	ll = [&log_likelihood](double *param, int *param_status,int model_status,mcmc_data_interface *interf,void * parameters){
			return log_likelihood(param, interf,parameters);};
	std::function<double(double*,int*,int ,mcmc_data_interface *,void *)> lp =NULL;
	lp = [&log_prior](double *param, int *param_status,int model_status,mcmc_data_interface *interf,void * parameters){
			return log_prior(param, interf,parameters);};
	std::function<void(double*,int*,int ,double**,mcmc_data_interface *,void *)> f =NULL;
	if(fisher){
		f = [&fisher](double *param, int *param_status,int model_status, double **fisherm, mcmc_data_interface *interf,void * parameters){
			fisher(param, fisherm,interf,parameters);};
	}
	PTMCMC_MH_dynamic_PT_alloc_uncorrelated_internal(sampler_output,
			output,
			dimension,
			N_steps,
			chain_N,
			max_chain_N_thermo_ensemble,
			initial_pos,
			seeding_var,
			ensemble_initial_pos,
			chain_temps,
			swp_freq,
			t0,
			nu,
			max_chunk_size,
			chain_distribution_scheme,
			lp,
			ll,
			f,
			user_parameters,
			numThreads,
			pool,
			show_prog,
			statistics_filename,
			chain_filename,
			likelihood_log_filename,
			checkpoint_file);
	return ;
}
//######################################################################################
//######################################################################################
void RJPTMCMC_MH(double ***output, /**< [out] Output chains, shape is double[chain_N, N_steps,dimension]*/
	int ***parameter_status, /**< [out] Parameter status for each step corresponding to the output chains, shape is double[chain_N, N_steps,dimension]*/
	int **model_status,
	int nested_model_number,
	int max_dimension, 	/**< maximum dimension of the parameter space being explored -- only consideration is memory, as memory scales with dimension. Keep this reasonable, unless memory is REALLY not an issue*/
	int min_dimension, 	/**< minimum dimension of the parameter space being explored >=1*/
	int N_steps,	/**< Number of total steps to be taken, per chain*/
	int chain_N,	/**< Number of chains*/
	double *initial_pos, 	/**<Initial position in parameter space - shape double[dimension]*/
	int *initial_status, 	/**<Initial status of the parameters in the initial position in parameter space - shape int[max_dim]*/
	int initial_model_status, 	/**<Initial status of the parameters in the initial position in parameter space - shape int[max_dim]*/
	double *seeding_var, 	/**<Variance of the normal distribution used to seed each chain higher than 0 - shape double[max_dimension] -- initial seeding of zero corresponds to the dimension turned off initially*/
	double **ensemble_initial_pos,
	int **ensemble_initial_status,
	int *ensemble_initial_model_status,
	double *chain_temps,	/**<Double array of temperatures for the chains*/
	int swp_freq,	/**< the frequency with which chains are swapped*/
	double (*log_prior)(double *param, int *status, int model_status,mcmc_data_interface *interface,void *parameters),	
	double (*log_likelihood)(double *param, int *status,int model_status, mcmc_data_interface *interface,void *parameters),
	void (*fisher)(double *param, int *status,int model_status, double **fisher, mcmc_data_interface *interface,void *parameters),
	void(*RJ_proposal)(double *current_param, double *proposed_param, int *current_status, int *proposed_status,int *current_model_status, int *proposed_model_status,double *MH_corrections, mcmc_data_interface *interface,void *parameters),/**< std::function for the log_likelihood function -- takes double *position, int *param_status,int dimension, int chain_id*/
	void **user_parameters,/**< Void pointer to any parameters the user may need inside log_prior, log_likelihood, or fisher. Should have one pointer for each chain. If this isn't needed, use (void**) NULL**/
	int numThreads, /**< Number of threads to use (=1 is single threaded)*/
	bool pool, /**< boolean to use stochastic chain swapping (MUST have >2 threads)*/
	bool show_prog, /**< boolean whether to print out progress (for example, should be set to ``false'' if submitting to a cluster)*/
	std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
	std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output*/
	std::string auto_corr_filename,/**< Filename to output auto correlation in some interval, if empty string, not output*/
	std::string likelihood_log_filename,/**< Filename to write the log_likelihood and log_prior at each step -- use empty string to skip*/
	std::string checkpoint_file/**< Filename to output data for checkpoint, if empty string, not saved*/
	)
{
	std::function<double(double*,int*,int,mcmc_data_interface *,void *)> ll =NULL;
	ll = [&log_likelihood](double *param, int *param_status,int model_status,mcmc_data_interface *interf,void *parameters){
			return log_likelihood(param, param_status,model_status,interf,parameters);};
	std::function<double(double*,int*,int ,mcmc_data_interface *,void *)> lp =NULL;
	lp = [&log_prior](double *param, int *param_status,int  model_status,mcmc_data_interface *interf,void *parameters){
			return log_prior(param, param_status,model_status, interf,parameters);};
	std::function<void(double*,int*,int ,double**,mcmc_data_interface *,void *)> f =NULL;
	if(fisher){
		f = [&fisher](double *param, int *param_status,int model_status, double **fisherm, mcmc_data_interface *interf,void *parameters){
			fisher(param, param_status,model_status, fisherm,interf,parameters);};
	}
	std::function<void(double*,double*, int*,int*,int*,int*,double*,mcmc_data_interface *, void *)> rj =NULL;
	rj = [&RJ_proposal](double *current_param, double *prop_param,int *current_param_status,int *prop_param_status,int *current_model_status,int *prop_model_status,double *MH_corrections,mcmc_data_interface *interf, void *parameters){
			RJ_proposal(current_param, prop_param,current_param_status, prop_param_status,current_model_status, prop_model_status,MH_corrections,interf,parameters);};
	RJPTMCMC_MH_internal(output, 
		parameter_status, 	
		model_status,
		nested_model_number,
		max_dimension, 	
		min_dimension, 	
		N_steps,	
		chain_N,	
		initial_pos, 	
		initial_status, 	
		initial_model_status, 	
		seeding_var, 	
		ensemble_initial_pos,
		ensemble_initial_status,
		ensemble_initial_model_status,
		chain_temps,	
		swp_freq,	
		lp,
		ll,
		f,
		rj,
		user_parameters,
		numThreads, 
		pool, 
		show_prog, 
		false, 
		statistics_filename,
		chain_filename,
		auto_corr_filename,
		likelihood_log_filename,
		checkpoint_file);
}
//######################################################################################
//######################################################################################
void continue_RJPTMCMC_MH(std::string start_checkpoint_file,/**< File for starting checkpoint*/
	double ***output,/**< [out] output array, dimensions: output[chain_N][N_steps][dimension]*/
	int ***status,/**< [out] output parameter status array, dimensions: status[chain_N][N_steps][dimension]*/
	int **model_status,
	int nested_model_number,
	int N_steps,/**< Number of new steps to take*/
	int swp_freq,/**< frequency of swap attempts between temperatures*/
	double (*log_prior)(double *param, int *status,int model_status, mcmc_data_interface *interface,void *parameters),	
	double (*log_likelihood)(double *param, int *status,int model_status, mcmc_data_interface *interface,void *parameters),
	void (*fisher)(double *param, int *status,int model_status, double **fisher, mcmc_data_interface *interface,void *parameters),
	void(*RJ_proposal)(double *current_param, double *proposed_param, int *current_status, int *proposed_status,int *current_model_status, int *proposed_model_status,double *MH_corrections, mcmc_data_interface *interface,void *parameters),/**< std::function for the log_likelihood function -- takes double *position, int *param_status,int dimension, int chain_id*/
	void **user_parameters,/**< Void pointer to any parameters the user may need inside log_prior, log_likelihood, or fisher. Should have one pointer for each chain. If this isn't needed, use (void**) NULL**/
	int numThreads,/**<Number of threads to use*/
	bool pool,/**<Boolean for whether to use ``deterministic'' vs ``stochastic'' sampling*/
	bool show_prog,/**< Boolean for whether to show progress or not (turn off for cluster runs*/
	std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
	std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output*/
	std::string auto_corr_filename,/**< Filename to output auto correlation in some interval, if empty string, not output*/
	std::string likelihood_log_filename,/**< Filename to write the log_likelihood and log_prior at each step -- use empty string to skip*/
	std::string end_checkpoint_file/**< Filename to output data for checkpoint at the end of the continued run, if empty string, not saved*/
	)
{

	std::function<double(double*,int*,int ,mcmc_data_interface *,void *)> ll =NULL;
	ll = [&log_likelihood](double *param, int *param_status,int model_status,mcmc_data_interface *interf,void *parameters){
			return log_likelihood(param, param_status,model_status,interf ,parameters);};
	std::function<double(double*,int*,int ,mcmc_data_interface *,void *)> lp =NULL;
	lp = [&log_prior](double *param, int *param_status,int model_status,mcmc_data_interface *interf,void *parameters){
			return log_prior(param, param_status, model_status,interf,parameters);};
	std::function<void(double*,int*,int,double**,mcmc_data_interface *,void *)> f =NULL;
	if(fisher){
		f = [&fisher](double *param, int *param_status,int model_status, double **fisherm, mcmc_data_interface *interf,void *parameters){
			fisher(param, param_status,model_status, fisherm,interf,parameters);};
	}
	std::function<void(double*,double*, int*,int*,int*,int*,double*,mcmc_data_interface *, void*)> rj =NULL;
	rj = [&RJ_proposal](double *current_param, double *prop_param,int *current_param_status,int *prop_param_status,int *current_model_status,int *prop_model_status,double *MH_corrections,mcmc_data_interface *interf,void *parameters){
			RJ_proposal(current_param, prop_param,current_param_status, prop_param_status,current_model_status, prop_model_status,MH_corrections,interf,parameters);};
	sampler samplerref;
	continue_RJPTMCMC_MH_internal(&samplerref,
		start_checkpoint_file,
		output,
		status,
		model_status,
		nested_model_number,
		N_steps,
		swp_freq,
		lp,
		ll,
		f,
		rj,
		user_parameters,
		numThreads,
		pool,
		show_prog,
		false,
		statistics_filename,
		chain_filename,
		auto_corr_filename,
		likelihood_log_filename,
		end_checkpoint_file,
		true,
		false);
	deallocate_sampler_mem(&samplerref);
}
//###################################################################################
void RJPTMCMC_MH_dynamic_PT_alloc_comprehensive(mcmc_sampler_output *sampler_output,
	double **output, 
	int **parameter_status, 
	int *model_status,
	int nested_model_number,
	int max_dimension, 	
	int min_dimension, 	
	int N_steps,	
	int chain_N,
	int max_chain_N_thermo_ensemble,
	double *initial_pos, 	
	int *initial_status, 	
	int initial_model_status, 	
	double *seeding_var, 	
	double **ensemble_initial_pos,
	int **ensemble_initial_status,
	int *ensemble_initial_model_status,
	double *chain_temps, 
	int swp_freq,	
	int t0,
	int nu,
	int max_chunk_size,
	std::string chain_distribution_scheme, 
	double (*log_prior)(double *param, int *status, int model_status,mcmc_data_interface *interface, void * parameters),	
	double (*log_likelihood)(double *param, int *status,int model_status, mcmc_data_interface *interface, void * parameters),
	void (*fisher)(double *param, int *status,int model_status, double **fisher, mcmc_data_interface *interface, void * parameters),
	void(*RJ_proposal)(double *current_param, double *proposed_param, int *current_status, int *proposed_status,int *current_model_status, int *proposed_model_status,double *MH_corrections, mcmc_data_interface *interface,  void * parameters),
	void **user_parameters,
	int numThreads, 
	bool pool, 
	bool show_prog, 
	bool update_RJ_width, 
	std::string statistics_filename,
	std::string chain_filename,
	std::string likelihood_log_filename,
	std::string checkpoint_file
	)
{
	std::function<double(double*,int*,int ,mcmc_data_interface *,void *)> ll =NULL;
	ll = [&log_likelihood](double *param, int *param_status,int model_status,mcmc_data_interface *interf,void *parameters){
			return log_likelihood(param, param_status,model_status,interf ,parameters);};
	std::function<double(double*,int*,int,mcmc_data_interface *,void *)> lp =NULL;
	lp = [&log_prior](double *param, int *param_status,int model_status,mcmc_data_interface *interf,void *parameters){
			return log_prior(param, param_status, model_status,interf,parameters);};
	std::function<void(double*,int*,int ,double**,mcmc_data_interface *,void *)> f =NULL;
	if(fisher){
		f = [&fisher](double *param, int *param_status,int model_status, double **fisherm, mcmc_data_interface *interf,void *parameters){
			fisher(param, param_status, model_status,fisherm,interf,parameters);};
	}
	std::function<void(double*,double*, int*,int*,int*,int*,double*,mcmc_data_interface *, void*)> rj =NULL;
	rj = [&RJ_proposal](double *current_param, double *prop_param,int *current_param_status,int *prop_param_status,int *current_model_status,int *prop_model_status,double *MH_corrections,mcmc_data_interface *interf,void *parameters){
			RJ_proposal(current_param, prop_param,current_param_status, prop_param_status,current_model_status, prop_model_status,MH_corrections,interf,parameters);};
	RJPTMCMC_MH_dynamic_PT_alloc_comprehensive_internal(sampler_output,
		output,
		parameter_status,
		model_status,
		nested_model_number,
		max_dimension,
		min_dimension,
		N_steps,
		chain_N,
		max_chain_N_thermo_ensemble,
		initial_pos,
		initial_status,
		initial_model_status,
		seeding_var,
		ensemble_initial_pos,
		ensemble_initial_status,
		ensemble_initial_model_status,
		chain_temps,
		swp_freq,
		t0, 
		nu, 
		max_chunk_size,
		chain_distribution_scheme,
		lp,
		ll,
		f,
		rj,
		user_parameters,
		numThreads,
		pool,
		show_prog,
		update_RJ_width,
		statistics_filename,
		chain_filename,
		likelihood_log_filename,
		checkpoint_file);
}
//#############################################################################
void continue_RJPTMCMC_MH_dynamic_PT_alloc_comprehensive(
	std::string checkpoint_file_start,
	mcmc_sampler_output *sampler_output,
	double **output, 
	int **parameter_status, 
	int *model_status,
	int nested_model_number,
	int N_steps,	
	int max_chain_N_thermo_ensemble,
	double *chain_temps, 
	int swp_freq,	
	int t0,
	int nu,
	int max_chunk_size,
	std::string chain_distribution_scheme, 
	double (*log_prior)(double *param, int *status, int model_status,mcmc_data_interface *interface, void * parameters),	
	double (*log_likelihood)(double *param, int *status,int model_status, mcmc_data_interface *interface, void * parameters),
	void (*fisher)(double *param, int *status,int model_status, double **fisher, mcmc_data_interface *interface, void * parameters),
	void(*RJ_proposal)(double *current_param, double *proposed_param, int *current_status, int *proposed_status,int *current_model_status, int *proposed_model_status,double *MH_corrections, mcmc_data_interface *interface,  void * parameters),
	void **user_parameters,
	int numThreads, 
	bool pool, 
	bool show_prog, 
	bool update_RJ_width, 
	std::string statistics_filename,
	std::string chain_filename,
	std::string likelihood_log_filename,
	std::string checkpoint_file
	)
{
	std::function<double(double*,int*,int ,mcmc_data_interface *,void *)> ll =NULL;
	ll = [&log_likelihood](double *param, int *param_status,int model_status,mcmc_data_interface *interf,void *parameters){
			return log_likelihood(param, param_status,model_status,interf ,parameters);};
	std::function<double(double*,int*,int,mcmc_data_interface *,void *)> lp =NULL;
	lp = [&log_prior](double *param, int *param_status,int model_status,mcmc_data_interface *interf,void *parameters){
			return log_prior(param, param_status, model_status,interf,parameters);};
	std::function<void(double*,int*,int ,double**,mcmc_data_interface *,void *)> f =NULL;
	if(fisher){
		f = [&fisher](double *param, int *param_status,int model_status, double **fisherm, mcmc_data_interface *interf,void *parameters){
			fisher(param, param_status, model_status,fisherm,interf,parameters);};
	}
	std::function<void(double*,double*, int*,int*,int*,int*,double*,mcmc_data_interface *, void*)> rj =NULL;
	rj = [&RJ_proposal](double *current_param, double *prop_param,int *current_param_status,int *prop_param_status,int *current_model_status,int *prop_model_status,double *MH_corrections,mcmc_data_interface *interf,void *parameters){
			RJ_proposal(current_param, prop_param,current_param_status, prop_param_status,current_model_status, prop_model_status,MH_corrections,interf,parameters);};
	continue_RJPTMCMC_MH_dynamic_PT_alloc_comprehensive_internal(
		checkpoint_file_start, 
		sampler_output,
		output,
		parameter_status,
		model_status,
		nested_model_number,
		N_steps,
		max_chain_N_thermo_ensemble,
		chain_temps,
		swp_freq,
		t0, 
		nu, 
		max_chunk_size,
		chain_distribution_scheme,
		lp,
		ll,
		f,
		rj,
		user_parameters,
		numThreads,
		pool,
		show_prog,
		update_RJ_width,
		statistics_filename,
		chain_filename,
		likelihood_log_filename,
		checkpoint_file);
}
