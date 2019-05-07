#ifndef MCMC_SAMPLER_H
#define MCMC_SAMPLER_H
#include <iostream>
void testmcmc();

void MCMC_MH(	double ***output, 
		int dimension, 	
		int N_steps,	
		int chain_N,	
		double *initial_pos, 	
		double *chain_temps,	
		int swp_freq,	
		double (*log_prior)(double *param, int dimension),	
		double (*log_likelihood)(double *param, int dimension),	
		void (*fisher)(double *param, int dimension, double **fisher),
		std::string statistics_filename,
		std::string chain_filename,
		std::string auto_corr_filename
		);
#endif
