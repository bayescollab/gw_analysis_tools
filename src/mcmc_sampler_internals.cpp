#include "mcmc_sampler_internals.h"
#include "autocorrelation.h"
#include "util.h"
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <eigen3/Eigen/Eigen>
#include <limits>
#include <iomanip>
#include <fftw3.h>
#include <stdio.h>
#ifndef _OPENMP
#define omp ignore
#endif

#ifdef _OPENMP
#include <adolc/adolc_openmp.h>
#include <omp.h>
#endif

/*! \file
 * File containing definitions for all the internal, generic mcmc subroutines
 * TODO update the acceptance and count of KDE step
 */

/*! \brief interface function between the sampler and the internal step functions
 */
int mcmc_step(sampler *sampler, double *current_param, double *next_param, int *current_status, int *next_status,int *current_model_status, int *next_model_status, int chain_number)
{
	//Random number to determine type of step
	double alpha = gsl_rng_uniform(sampler->rvec[chain_number]);
	

	double proposed_param[sampler->max_dim];
	int proposed_status[sampler->max_dim];
	int proposed_model_status= *current_model_status;

	int step, selected_dimension;	
	//Determine which step to take and calculate proposed coord.
	if (alpha<sampler->prob_boundaries[chain_number][0])
	{
		gaussian_step(sampler, current_param, proposed_param, current_status, proposed_status,current_model_status, &proposed_model_status, chain_number,&selected_dimension);
		sampler->num_gauss[chain_number]+=1;
		step = 0;
	}
	else if (alpha<sampler->prob_boundaries[chain_number][1])
	{
		diff_ev_step(sampler, current_param, proposed_param, current_status, proposed_status,current_model_status, &proposed_model_status, chain_number);
		sampler->num_de[chain_number]+=1;
		step = 1;
	}
	else if (alpha<sampler->prob_boundaries[chain_number][2])
	{
		mmala_step(sampler, current_param, proposed_param, current_status, proposed_status,current_model_status, &proposed_model_status,chain_number);
		sampler->num_mmala[chain_number]+=1;
		step= 2;
	}
	else if (alpha<sampler->prob_boundaries[chain_number][3])
	{
		fisher_step(sampler, current_param, proposed_param, current_status, proposed_status,current_model_status, &proposed_model_status, chain_number);
		sampler->num_fish[chain_number]+=1;
		step = 3;
	}
	else if (alpha<sampler->prob_boundaries[chain_number][4])
	{
		RJ_step(sampler, current_param, proposed_param, current_status, proposed_status,current_model_status, &proposed_model_status, &(sampler->prop_MH_factor[chain_number]),chain_number);
		sampler->num_RJstep[chain_number]+=1;
		step = 4;
	}
	else 
	{
		kde_proposal(sampler, current_param, proposed_param, current_status, proposed_status,current_model_status, &proposed_model_status, chain_number);
		sampler->num_kde[chain_number]+=1;
		step = 5;
	}
	
	double current_lp = sampler->lp(current_param, current_status,*current_model_status,sampler->interfaces[chain_number], sampler->user_parameters[chain_number]);
	double proposed_lp = sampler->lp(proposed_param,proposed_status, proposed_model_status,sampler->interfaces[ chain_number], sampler->user_parameters[chain_number]);
	double current_ll=0, proposed_ll=0;
	double MH_ratio;
	double power;

	//Check if proposed step is out of range by calculating prior
	//if out of range, reject step
	if(current_lp == limit_inf || proposed_lp == limit_inf){
		MH_ratio = limit_inf;
	}
	else{
		//Calculate log_likelihood and log prior
		current_ll = sampler->current_likelihoods[chain_number];
		proposed_ll = sampler->ll(proposed_param, proposed_status,proposed_model_status,sampler->interfaces[chain_number], sampler->user_parameters[chain_number]);
		//proposed_ll = (proposed_ll )/sampler->chain_temps[chain_number];
		//Calculate MH ratio
		if(std::isnan(proposed_ll)){
			MH_ratio = limit_inf;
		}
		else{
			MH_ratio = (-current_ll+proposed_ll)/sampler->chain_temps[chain_number]-current_lp + proposed_lp;
		}
	}
	//Some proposals are not symmetric
	MH_ratio += sampler->prop_MH_factor[chain_number];
	//Reset proposal factor to 0 because some proposals assume symmetry
	sampler->prop_MH_factor[chain_number]=0;

	int i;
	//Random number to determine step acceptance
	double beta = log(gsl_rng_uniform(sampler->rvec[chain_number]));
	if(MH_ratio< beta){
		for ( i=0;i<sampler->max_dim; i ++)
		{
			next_param[i] = current_param[i];
			next_status[i] = current_status[i];
		}
		*next_model_status = *current_model_status;
		
		assign_ct_m(sampler, step,chain_number,selected_dimension);


		return -1;
	}	
	else
	{
		for ( i=0;i<sampler->max_dim; i ++)
		{
			next_param[i] = proposed_param[i];
			next_status[i] = proposed_status[i];
		}
		*next_model_status = proposed_model_status;
		assign_ct_p(sampler, step, chain_number, selected_dimension);
		sampler->current_likelihoods[chain_number] = proposed_ll;
		if(step == 3 && sampler->proper_fisher){	
			iterate_fisher(sampler,chain_number);
		}

		return 1;
	}		
	
}
double log_kde_evaluate(double *pos, int *status,int dim, sampler *samplerptr, int chain_id)
{
	double val = 0;
	double *diff =new double[dim];
	
	for (int k = 0 ; k<samplerptr->history_length; k++){
		double dot = 0;
		for(int i = 0 ; i<dim ; i++){
			diff[i] = pos[i] - samplerptr->history[chain_id][k][i];
		}
		
		for(int i = 0 ; i<dim ; i++){
			for(int j = 0 ; j<dim ; j++){
				dot+= diff[i] * samplerptr->kde_fisher[chain_id][i][j] * diff[j]; 
				//dot+= diff[i] * samplerptr->kde_cov[chain_id][i][j] * diff[j]; 
			}
		}
		val+=exp(dot*(-.5));
	}
	double logval = log(val);
	logval += -.5*dim * log(2*M_PI) - .5*samplerptr->kde_cov_lndet[chain_id] - log(samplerptr->history_length);
	
	delete [] diff;
	return logval;
}

	
/*! \brief KDE step from history file
 *
 * TODO: Right now, this assumes a diagonal covariance matrix. Maybe fix this in the future
 */
void kde_proposal(sampler *sampler, /**< Sampler struct*/
		double *current_param, /**< current position in parameter space*/
		double *proposed_param, /**< [out] Proposed position in parameter space*/
		int *current_status,
		int *proposed_status,
		int *current_model_status,
		int *proposed_model_status,
		int chain_id
		)
{
	for (int i=0;i<sampler->max_dim;i++){
		if(current_status[i] == 1){
			proposed_param[i] = current_param[i];
		}
		else{
			proposed_param[i] = 0;
		}
		proposed_status[i] = current_status[i];
	}

	sampler->kde_cov_update_ct[chain_id] +=1;
	if(sampler->kde_cov_update_ct[chain_id] >= sampler->kde_cov_update_number){
		if(sampler->kde_step){
			update_kde_cov(sampler, chain_id);
		}
		sampler->kde_cov_update_ct[chain_id] = 0;
	}
	int beta = 0;
	beta = int(gsl_rng_uniform(sampler->rvec[chain_id]) * sampler->history_length);

	//double alpha = sampler->randgauss_width[chain_id][5][beta];
	if(sampler->RJMCMC){
		for (int i=0;i<sampler->max_dim;i++){
			if(current_status[i] == 1){
				proposed_param[i] = current_param[i];
			}
			else{
				proposed_param[i] = 0;
			}
			proposed_status[i] = current_status[i];
		}

	}
	
	else{
		double **temp_cov = new double*[sampler->max_dim];
		for (int i=0;i<sampler->max_dim;i++){
			temp_cov[i]  = new double[sampler->max_dim];
			for (int j=0;j<sampler->max_dim;j++){
				//temp_cov[i][j] = sampler->kde_cov[chain_id][i][j]*sampler->randgauss_width[chain_id][5][0];
				temp_cov[i][j] = sampler->kde_cov[chain_id][i][j];
			}
		}
		double **temp_prop = new double*[1];
		temp_prop[0] = new double[sampler->max_dim];
		mvn_sample(1, sampler->history[chain_id][beta],temp_cov, sampler->max_dim, sampler->rvec[chain_id],temp_prop);
		//for(int i = 0 ; i<sampler->max_dim; i++){
		//	temp_prop[0][i] = gsl_ran_gaussian(sampler->rvec[chain_id], sqrt(temp_cov[i][i])) + sampler->history[chain_id][beta][i];
		//}

		for (int i=0;i<sampler->max_dim;i++){
			delete [] temp_cov[i];
			if(current_status[i] == 1){
				//std::cout<<temp_prop[0][i]<<" ";
				proposed_param[i] = temp_prop[0][i];
			}
			else{
				proposed_param[i] = 0;
			}
			//proposed_status[i] = current_status[i];
		}
		//std::cout<<std::endl;
		delete [] temp_cov;
		delete [] temp_prop[0];
		delete [] temp_prop;
	}

	*proposed_model_status = *current_model_status;	
	//Need to calculate log( current_kde_eval / prop_kde_eval)
	double log_kde_val_c = log_kde_evaluate(current_param,current_status, sampler->max_dim, sampler,chain_id);
	double log_kde_val_p = log_kde_evaluate(proposed_param,proposed_status, sampler->max_dim, sampler,chain_id);
	sampler->prop_MH_factor[chain_id] =log_kde_val_c - log_kde_val_p ;
	//sampler->prop_MH_factor[chain_id] =log_kde_val_p - log_kde_val_c ;
	return;
}

void update_kde_cov(sampler *sampler,int chain_id)
{
	
	double **cov = sampler->kde_cov[chain_id];
	int dim = sampler->max_dim;

	double means[dim];
	for(int i = 0 ; i<dim ; i++){
		int ct = 0 ;
		means[i]=0;
		for(int j = 0 ; j<sampler->history_length; j++){
			if(sampler->history_status[chain_id][j][i] == 1){
				means[i]+=sampler->history[chain_id][j][i];
				ct++;
			}
		}
		if(ct !=0){
			means[i]/=ct;
		}
	}
	
	for(int i = 0 ; i<dim ; i++){
		for(int j = 0 ; j<=i; j++){
			sampler->kde_cov[chain_id][i][j]=0;
			int ct = 0 ;
			for(int k = 0 ; k<sampler->history_length; k++){
				if(sampler->history_status[chain_id][k][i] == 1 & sampler->history_status[chain_id][k][j] == 1){
					sampler->kde_cov[chain_id][i][j]+=(sampler->history[chain_id][k][i] - means[i])*
						(sampler->history[chain_id][k][j]-means[j]);
					ct++;
				}
			}
			if(ct != 0){
				sampler->kde_cov[chain_id][i][j]/=ct;
			}
		}
	}	
	for(int i=0 ; i<dim; i++){
		for(int j = i ; j<dim; j++){
			sampler->kde_cov[chain_id][i][j] = sampler->kde_cov[chain_id][j][i];
		}
	}
	//for(int i = 0 ; i<dim ; i++){
	//	sampler->kde_cov[chain_id][i][i]=0;
	//	int ct = 0 ;
	//	for(int k = 0 ; k<sampler->history_length; k++){
	//		if(sampler->history_status[chain_id][k][i] == 1){
	//			sampler->kde_cov[chain_id][i][i]+=(sampler->history[chain_id][k][i] - means[i])*
	//				(sampler->history[chain_id][k][i]-means[i]);
	//			ct++;
	//		}
	//	}
	//	if(ct != 0){
	//		sampler->kde_cov[chain_id][i][i]/=ct;
	//	}
	//	
	//}	
	//############################################

	for(int i=0 ; i<dim; i++){
		for(int j = 0 ; j<dim; j++){
			//if(i != j ){
			//	sampler->kde_cov[chain_id][i][j] = 0;
			//}
			//else {
			//	sampler->kde_cov[chain_id][i][j] *= pow_int(.9 * pow(sampler->history_length, -1./(sampler->max_dim + 4)),2);
			//}
			//sampler->kde_cov[chain_id][i][j] *= pow_int(.9 * pow(sampler->history_length, -1./(sampler->max_dim + 4)),2)*sampler->randgauss_width[chain_id][5][0];
			sampler->kde_cov[chain_id][i][j] *= pow_int(.9 * pow(sampler->history_length, -1./(sampler->max_dim + 4)),2);
			//std::cout<<sampler->kde_cov[chain_id][i][j]<<" ";
		}
		//std::cout<<std::endl;
	}
	//std::cout<<std::endl;
	//############################################
	//for(int i=0 ; i<dim; i++){
	//	for(int j = 0 ; j<dim; j++){
	//		sampler->kde_cov[chain_id][i][j] *= .9 * pow(sampler->history_length, -1./(sampler->max_dim + 4));
	//	}
	//}
	gsl_LU_matrix_invert(sampler->kde_cov[chain_id], sampler->kde_fisher[chain_id], dim);
	//update kde_cov_lndet;
	sampler->kde_cov_lndet[chain_id] = gsl_LU_lndet(sampler->kde_cov[chain_id],dim);
	//############################################

	//for(int i=0 ; i<dim; i++){
	//	for(int j = 0 ; j<dim; j++){
	//		std::cout<<sampler->kde_fisher[chain_id][i][j]<<" ";
	//	}
	//	std::cout<<std::endl;
	//}
	//############################################
	
	return;
}

/*! \brief Straight gaussian step
 *
 * Change this to pick one direction 
 */
void gaussian_step(sampler *sampler, /**< Sampler struct*/
		double *current_param, /**< current position in parameter space*/
		double *proposed_param, /**< [out] Proposed position in parameter space*/
		int *current_status,
		int *proposed_status,
		int *current_model_status,
		int *proposed_model_status,
		int chain_id,
		int *selected_dimension
		)
{
	int beta = 0;
	do{
		beta = gsl_rng_uniform(sampler->rvec[chain_id]) * sampler->max_dim;
	}while(current_status[beta] == 0);
	*selected_dimension = beta;

	double alpha = sampler->randgauss_width[chain_id][0][beta];
	for (int i=0;i<sampler->max_dim;i++){
		if(current_status[i] == 1){
			proposed_param[i] = current_param[i];
		}
		else{
			proposed_param[i] = 0;
		}
		proposed_status[i] = current_status[i];
	}
	
	proposed_param[beta] = gsl_ran_gaussian(sampler->rvec[chain_id], alpha)+current_param[beta];

	*proposed_model_status = *current_model_status;	

	//double alpha = sampler->randgauss_width[chain_id][0][0];
	//for (int i=0;i<sampler->max_dim;i++){
	//	if(current_status[i] == 1){
	//		proposed_param[i] = current_param[i];
	//	}
	//	else{
	//		proposed_param[i] = 0;
	//	}
	//	proposed_status[i] = current_status[i];
	//}
	//
	//proposed_param[beta] = gsl_ran_gaussian(sampler->rvec[chain_id], alpha)+current_param[beta];
	//#####################################################
	//
	//
	//for (int i=0;i<sampler->max_dim;i++){
	//	if(current_status[i] == 1){
	//		proposed_param[i] = gsl_ran_gaussian(sampler->rvec[chain_id], alpha)+current_param[i];
	//	}
	//	else{
	//		proposed_param[i] = 0;
	//	}
	//	proposed_status[i] = current_status[i];
	//}
}

/*!\brief Fisher informed gaussian step
 */
void fisher_step(sampler *sampler, /**< Sampler struct*/
		double *current_param, /**< current position in parameter space*/
		double *proposed_param, /**< [out] Proposed position in parameter space*/
		int *current_status,
		int *proposed_status,
		int *current_model_status,
		int *proposed_model_status,
		int chain_index
		)
{
	if(sampler->fisher_update_ct[chain_index]==sampler->fisher_update_number )
	{
		update_fisher(sampler, current_param, current_status,current_model_status,chain_index);	
	}
	double scaling;
	int beta;
	double gamma;
	//TODO : Need to specialize this more with RJ and min_dim == 0 case
	//look for discrete model number
	if(!sampler->RJMCMC || sampler->min_dim ==0){
		if(sampler->block_sample){
					
			gamma = gsl_rng_uniform(sampler->rvec[chain_index]);
		}
		if(sampler->block_sample && gamma<=sampler->block_sample_prob){
			int delta = (int)((sampler->block_num)*gsl_rng_uniform(sampler->rvec[chain_index]));
			int range = 0, init_id = 0;
			if (delta == 0){range = sampler->block_boundary_ids[0];}
			else{ range = sampler->block_boundary_ids[delta] - sampler->block_boundary_ids[delta-1];init_id = sampler->block_boundary_ids[delta-1];}
			beta = (int)((range)*(gsl_rng_uniform(sampler->rvec[chain_index]))) + init_id;
			
			double alpha = gsl_ran_gaussian(sampler->rvec[chain_index],
						 sampler->randgauss_width[chain_index][3][0]);

			scaling = 0.0;
			//ensure the steps aren't ridiculous
			if(abs(sampler->fisher_vals[chain_index][beta])<10){scaling = 10.;}
			//else if(abs(sampler->fisher_vals[chain_index][beta])>1000){scaling = 1000.;}
			//##########################################################33
			//TESTING -- Scaling for annealing
			else{scaling = abs(sampler->fisher_vals[chain_index][beta])/
						sampler->chain_temps[chain_index];}
			//else{scaling = abs(sampler->fisher_vals[chain_index][beta]);}
			//##########################################################33
			//Take step
			for(int i =0; i< sampler->max_dim;i++)
			{
				proposed_param[i] = current_param[i] ;
				proposed_status[i] = current_status[i];
			}
			for(int i =init_id; i< init_id+range;i++)
			{
				proposed_param[i] = current_param[i] +
					alpha/sqrt(scaling) *sampler->fisher_vecs[chain_index][beta][i];
				proposed_status[i] = current_status[i];
			}
			*proposed_model_status = *current_model_status;	

		}
		else{
			//beta determines direction to step in eigen directions
			beta = (int)((sampler->max_dim)*(gsl_rng_uniform(sampler->rvec[chain_index])));
			
			double alpha = gsl_ran_gaussian(sampler->rvec[chain_index],
						 sampler->randgauss_width[chain_index][3][0]);

			scaling = 0.0;
			//ensure the steps aren't ridiculous
			if(abs(sampler->fisher_vals[chain_index][beta])<10){scaling = 10.;}
			//else if(abs(sampler->fisher_vals[chain_index][beta])>1000){scaling = 1000.;}
			//##########################################################33
			//TESTING -- Scaling for annealing
			else{scaling = abs(sampler->fisher_vals[chain_index][beta])/
						sampler->chain_temps[chain_index];}
			//else{scaling = abs(sampler->fisher_vals[chain_index][beta]);}
			//##########################################################33
			//Take step
			for(int i =0; i< sampler->max_dim;i++)
			{
				proposed_param[i] = current_param[i] +
					alpha/sqrt(scaling) *sampler->fisher_vecs[chain_index][beta][i];
				proposed_status[i] = current_status[i];
			}
			*proposed_model_status = *current_model_status;	
		}
		
	}
	//If RJPTMCMC and there's a base model, use the fisher for the base model, and gaussian steps for the modifications
	else {
		if(sampler->block_sample){
					
			gamma = gsl_rng_uniform(sampler->rvec[chain_index]);
		}
		if(sampler->block_sample && gamma<=sampler->block_sample_prob){
			int delta = (int)((sampler->block_num)*gsl_rng_uniform(sampler->rvec[chain_index]));
			int range = 0, init_id = 0;
			if (delta == 0){range = sampler->block_boundary_ids[0];}
			else{ range = sampler->block_boundary_ids[delta] - sampler->block_boundary_ids[delta-1];init_id = sampler->block_boundary_ids[delta-1];}
			beta = (int)((range)*(gsl_rng_uniform(sampler->rvec[chain_index]))) + init_id;
			
			double alpha = gsl_ran_gaussian(sampler->rvec[chain_index],
						 sampler->randgauss_width[chain_index][3][0]);

			scaling = 0.0;
			//ensure the steps aren't ridiculous
			if(abs(sampler->fisher_vals[chain_index][beta])<10){scaling = 10.;}
			//else if(abs(sampler->fisher_vals[chain_index][beta])>10000){scaling = 1000.;}

			else{scaling = abs(sampler->fisher_vals[chain_index][beta])/
						sampler->chain_temps[chain_index];}
			//Take step
			for(int i =0; i< sampler->max_dim;i++)
			{
				proposed_param[i] = current_param[i];
				proposed_status[i] = current_status[i];
			}


			for(int i =init_id; i< range+init_id;i++)
			{
				if(i < sampler->min_dim){
					proposed_param[i] = current_param[i] +
						alpha/sqrt(scaling) *sampler->fisher_vecs[chain_index][beta][i];
					proposed_status[i] = current_status[i];
				}
				else{
					if(current_status[i] == 1){
						alpha = gsl_ran_gaussian(sampler->rvec[chain_index],
							sampler->randgauss_width[chain_index][0][i]);
						proposed_param[i] = alpha+current_param[i];
						proposed_status[i] = current_status[i];
					}
					else{
						proposed_param[i] = 0;
						proposed_status[i] = current_status[i];
					}
				}
			}
			*proposed_model_status = *current_model_status;	
		
		}
		else{
			//beta determines direction to step in eigen directions
			beta = (int)((sampler->min_dim)*(gsl_rng_uniform(sampler->rvec[chain_index])));
			
			double alpha = gsl_ran_gaussian(sampler->rvec[chain_index],
						 sampler->randgauss_width[chain_index][3][0]);

			scaling = 0.0;
			//ensure the steps aren't ridiculous
			if(abs(sampler->fisher_vals[chain_index][beta])<10){scaling = 10.;}
			//else if(abs(sampler->fisher_vals[chain_index][beta])>10000){scaling = 1000.;}

			else{scaling = abs(sampler->fisher_vals[chain_index][beta])/
						sampler->chain_temps[chain_index];}
			//Take step
			for(int i =0; i< sampler->min_dim;i++)
			{
				proposed_param[i] = current_param[i] +
					alpha/sqrt(scaling) *sampler->fisher_vecs[chain_index][beta][i];
				proposed_status[i] = current_status[i];
					//std::cout<<"Delta p "<<i<<": "<<alpha/sqrt(scaling) *sampler->fisher_vecs[chain_index][beta][i]<<std::endl;
			}
			//Generate new step for gaussian steps, using gaussian width	
			for(int i =sampler->min_dim; i< sampler->max_dim;i++)
			{
				//if(current_status[i] == 1){
				//	alpha = gsl_ran_gaussian(sampler->rvec[chain_index],
				//		sampler->randgauss_width[chain_index][0][i]);
				//	proposed_param[i] = alpha+current_param[i];
				//	proposed_status[i] = current_status[i];
				//	//std::cout<<"Delta p "<<i<<": "<<alpha<<std::endl;
				//}
				//else{
				//	proposed_param[i] = 0;
				//	proposed_status[i] = current_status[i];
				//	//std::cout<<"NO STEP "<<i<<": "<<std::endl;
				//}
				proposed_param[i] = current_param[i];
				proposed_status[i] = current_status[i];
			}
			*proposed_model_status = *current_model_status;	
		}	
	}
	if(sampler->proper_fisher){
		double lp =sampler->lp(proposed_param,proposed_status, *proposed_model_status,sampler->interfaces[chain_index], sampler->user_parameters[chain_index]) ;
		if(lp != limit_inf){
			update_fisher(sampler, proposed_param, proposed_status,proposed_model_status, chain_index);
			sampler->prop_MH_factor[chain_index] = 0;
			for(int i = 0 ; i<sampler->max_dim ; i++){
				for(int j = 0 ; j<sampler->max_dim ; j++){
					if(proposed_status[i] == 1 && proposed_status[j] == 1 &&
						current_status[i] == 1 && current_status[j] ==1){
						sampler->prop_MH_factor[chain_index] += 
							-0.5 *(current_param[i] - proposed_param[i])*
							(current_param[j] - proposed_param[j])*
							(sampler->fisher_matrix[chain_index][i][j] -
							sampler->fisher_matrix_prop[chain_index][i][j]) ;
					}
				}
			}
		}
	}
	sampler->fisher_update_ct[chain_index] += 1;

}

void iterate_fisher(sampler *samplerptr,int chain_id)
{
	for(int i = 0 ; i<samplerptr->max_dim; i++){
		for(int j = 0 ; j<samplerptr->max_dim; j++){
			samplerptr->fisher_matrix[chain_id][i][j] = samplerptr->fisher_matrix_prop[chain_id][i][j];
			samplerptr->fisher_vecs[chain_id][i][j] = samplerptr->fisher_vecs_prop[chain_id][i][j];
		}
		samplerptr->fisher_vals[chain_id][i] = samplerptr->fisher_vals_prop[chain_id][i];
	}
}


void update_fisher(sampler *sampler, double *current_param, int *param_status, int *model_status,int chain_index)
{
	int local_dim = sampler->max_dim;

	//In the case we have nested models, just calculate fisher for the base model
	if(sampler->RJMCMC && sampler->min_dim !=0) local_dim=sampler->min_dim;
	//Fisher calculation
	double **fisher=(double **)malloc(sizeof(double*)*local_dim);	
	for (int i =0; i<local_dim;i++){
		fisher[i] = (double*)malloc(sizeof(double)*local_dim);
	}
	sampler->fish(current_param, param_status,*model_status, fisher,sampler->interfaces[chain_index], sampler->user_parameters[chain_index]);

	//Convert to 1D array for Eigen
	double *oneDfisher=(double *)malloc(sizeof(double)*local_dim*local_dim);
	for (int i =0; i<local_dim;i++){
		for (int j = 0; j<local_dim; j++){
			oneDfisher[local_dim*i+j] = fisher[i][j];///
		}
		
	}
	
	//Find eigen vectors and eigen values
	Eigen::Map<Eigen::MatrixXd> m(oneDfisher,local_dim,local_dim);
 	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(m);
 	Eigen::MatrixXd eigen_vecs = eigensolver.eigenvectors();
 	Eigen::VectorXd eigen_vals = eigensolver.eigenvalues();
	
	//##############################################################
	//Numerical matrix inversion can be tricky - catch nans here and replace
	//with gaussian step just to not have the program crash because of 
	//one position in parameter space
	//
	//To see impact of nans, variable is stored in sampler->nan_counter
	//##############################################################
	int nansum = 0;
	for(int j = 0 ; j<local_dim; j++){
		for(int i =0; i<local_dim; i++)
			nansum+= std::isnan(eigen_vecs.col(j)(i));	
		nansum+= std::isnan(eigen_vals(j));
	}
	if(!nansum){
		if(!sampler->proper_fisher){
			for (int i =0; i < local_dim; i++)
			{
				for(int j = 0; j<local_dim; j++)
				{
					sampler->fisher_matrix[chain_index][i][j] = fisher[i][j];
					sampler->fisher_vecs[chain_index][i][j] = eigen_vecs.col(i)(j);
				}
				sampler->fisher_vals[chain_index][i]=eigen_vals[i];
			}
		}
		else{
			for (int i =0; i < local_dim; i++)
			{
				for(int j = 0; j<local_dim; j++)
				{
					sampler->fisher_matrix_prop[chain_index][i][j] = fisher[i][j];
					sampler->fisher_vecs_prop[chain_index][i][j] = eigen_vecs.col(i)(j);
				}
				sampler->fisher_vals_prop[chain_index][i]=eigen_vals[i];
			}
		}
		sampler->fisher_update_ct[chain_index]=0;
	}
	else{ 
		sampler->fisher_update_ct[chain_index]=sampler->fisher_update_number-1;
		sampler->nan_counter[chain_index]+=1;
		//if(sampler->nan_counter[chain_index] > 50){
		//	debugger_print(__FILE__,__LINE__,"WARNING -- NANs from Fisher inversion are over 50! You're Fisher's may be unstable!");
		//}
	}

	for (int i =0; i<local_dim;i++){
		free(fisher[i]);
	}
	free(fisher);
	free(oneDfisher);
}

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/*!!!!!!!!DON'T TRUST!!!!!!!!!!!!!!!*/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
void calc_grad(sampler *sampler,double *current_param,int *current_status, int *current_model_status,int chain_index,double *grad)
{
	double epsilon =1.e-1;
	double locationp[sampler->max_dim];
	for(int i = 0 ; i<sampler->max_dim; i++){
		locationp[i]=current_param[i];
	}
	for(int i = 0 ; i<sampler->max_dim; i++){
		locationp[i]+=epsilon; 
		double llp = sampler->ll(locationp, current_status, *current_model_status,sampler->interfaces[chain_index], sampler->user_parameters[chain_index]);
		grad[i] = (llp-sampler->current_likelihoods[chain_index])/(epsilon);
		locationp[i]= current_param[i]; 
	}

}

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/*!!!!!!!!DON'T TRUST!!!!!!!!!!!!!!!*/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/*!\brief MMALA informed step -- Currently not supported
 *
 * NOTE: This assumes the Fisher doesn't change between steps. The proposal ratio only accounts for the gradient of the log likelihood
 */
void mmala_step(sampler *sampler, /**< Sampler struct*/
		double *current_param, /**< current position in parameter space*/
		double *proposed_param, /**< [out] Proposed position in parameter space*/
		int *current_status,
		int *proposed_status,
		int *current_model_status,
		int *proposed_model_status,
		int chain_index
		)
{
	double scaling;
	int beta;
	if(!sampler->RJMCMC || sampler->min_dim ==0){
		//beta determines direction to step in eigen directions
		beta = (int)((sampler->max_dim)*(gsl_rng_uniform(sampler->rvec[chain_index])));
		
		double alpha = gsl_ran_gaussian(sampler->rvec[chain_index],
					 sampler->randgauss_width[chain_index][2][0]);

		scaling = 0.0;
		//ensure the steps aren't ridiculous
		if(abs(sampler->fisher_vals[chain_index][beta])<10){scaling = 10.;}
		//else if(abs(sampler->fisher_vals[chain_index][beta])>1000){scaling = 1000.;}
		else{scaling = abs(sampler->fisher_vals[chain_index][beta])/
					sampler->chain_temps[chain_index];}
		double grad[sampler->max_dim];
		calc_grad(sampler,current_param,current_status, current_model_status,chain_index,grad);
		double dotprod = 0;
		for(int i = 0 ; i<sampler->max_dim; i++){
			dotprod+=grad[i]*sampler->fisher_vecs[chain_index][beta][i];
		}
		//Take step
		for(int i =0; i< sampler->max_dim;i++)
		{
			proposed_param[i] = current_param[i] +
				alpha/sqrt(scaling) *sampler->fisher_vecs[chain_index][beta][i]+
				sampler->fisher_vecs[chain_index][beta][i]*dotprod/(2*scaling);
			proposed_status[i] = current_status[i];
		}
		*proposed_model_status = *current_model_status;	
		
		double grad_rev[sampler->max_dim];
		calc_grad(sampler,proposed_param,current_status, current_model_status,chain_index,grad_rev);
		double **cov = new double*[sampler->max_dim];
		for(int i = 0 ; i<sampler->max_dim; i++){
			cov[i]= new double[sampler->max_dim];
			
		}
		gsl_cholesky_matrix_invert(sampler->fisher_matrix[chain_index], cov, sampler->max_dim);
		sampler->prop_MH_factor[chain_index]= 0 ;
		double mean_forward[sampler->max_dim];
		double mean_reverse[sampler->max_dim];
		
		for(int i = 0 ; i<sampler->max_dim; i++){
			double shift_forward = 0 ;
			double shift_reverse = 0 ;
			for(int j = 0 ; j<sampler->max_dim; j++){
				shift_forward += grad[j]*cov[i][j];
				shift_reverse += grad_rev[j]*cov[i][j];
			}
			mean_forward[i] = current_param[i]+1./2. * shift_forward;
			mean_reverse[i]= proposed_param[i]+1./2. * shift_reverse;
			

		}
		for(int i = 0 ; i<sampler->max_dim; i++){
			for(int j = 0 ; j<sampler->max_dim; j++){
				sampler->prop_MH_factor[chain_index]+= 
					-.5*(mean_reverse[i] - current_param[i])*
						(mean_reverse[j]-current_param[j])*
						sampler->fisher_matrix[chain_index][i][j];
					+.5*(mean_forward[i] - proposed_param[i])*
						(mean_forward[j]-proposed_param[j])*
						sampler->fisher_matrix[chain_index][i][j];
					
			}
		}
		//sampler->prop_MH_factor[chain_index] = 0;
		if(std::isnan(sampler->prop_MH_factor[chain_index])){
			exit(1);
		}
		for(int i = 0  ; i<sampler->max_dim; i++){
			delete [] cov[i];
		}
		delete [] cov;
	}
		
}

/*!\brief differential evolution informed step
 *
 * Differential evolution uses the past history of the chain to inform the proposed step:
 *
 * Take the difference of two random, accepted previous steps and step along that with some step size,
 * determined by a gaussian
 */
void diff_ev_step(sampler *sampler, /**< Sampler struct*/
		double *current_param, /**< current position in parameter space*/
		double *proposed_param, /**< [out] Proposed position in parameter space*/
		int *current_status,
		int *proposed_status,
		int *current_model_status,
		int *proposed_model_status,
		int chain_id
		)
{
	double *eff_history_coord = new double[sampler->max_dim];
	int *eff_history_status = new int[sampler->max_dim];
	double *eff_history_coord2 = new double[sampler->max_dim];
	int *eff_history_status2 = new int[sampler->max_dim];
	//if(sampler->RJMCMC){
	double gamma = 0 ;	
	int range=sampler->max_dim, init_id=0;
	if(sampler->block_sample){
				
		gamma = gsl_rng_uniform(sampler->rvec[chain_id]);
	}
	if(sampler->block_sample && gamma<=sampler->block_sample_prob){
		int delta = (int)((sampler->block_num)*gsl_rng_uniform(sampler->rvec[chain_id]));
		if (delta == 0){range = sampler->block_boundary_ids[0];}
		else{ range = sampler->block_boundary_ids[delta] - sampler->block_boundary_ids[delta-1];init_id = sampler->block_boundary_ids[delta-1];}
	}
	if(sampler->RJMCMC && sampler->min_dim>0){

		//Pick a history member
		int i = (int)((sampler->history_length-1)
			*(gsl_rng_uniform(sampler->rvec[chain_id])));
		//Second position ID
		int j;
		do{
			j=(int)((sampler->history_length-1)*
				(gsl_rng_uniform(sampler->rvec[chain_id])));	
		}while(j==i);

		for(int k = 0 ; k <sampler->min_dim; k++){
			eff_history_coord[k] = sampler->history[chain_id][i][k];
			eff_history_coord2[k] = sampler->history[chain_id][j][k];
			eff_history_status[k] = sampler->param_status[chain_id][i][k];
			eff_history_status2[k] = sampler->param_status[chain_id][j][k];
		}	
		for(int k = sampler->min_dim ; k <sampler->max_dim; k++){
			eff_history_coord[k] = 0;
			eff_history_coord2[k] = 0;
			eff_history_status[k] = sampler->param_status[chain_id][i][k];
			eff_history_status2[k] = sampler->param_status[chain_id][j][k];
		}	
		//We'll almost definitely need to do some transformation to it

		//Since there's no way to feasibly store a substantial 
		//population of history items for every combination of dimensions:
		//I'll assume the additional dimensions have small effects on the model
		//described by dimensions 0-min_dim, so I pick an element from 
		//the history file that has the smallest deviation from the base model, 
		//but still has the additional dimension. This means that the extra dimensions
		//beyond min_dim are independent, but if the affects are indeed small, the benefit 
		//of DE should still be recovered (this is only a proposal, after all)
		//if(sampler->min_dim != 0){
		//	RJ_smooth_history(sampler, current_param,current_status, i, eff_history_coord, 
		//		eff_history_status,chain_id);
		//	RJ_smooth_history(sampler, current_param, current_status, j, eff_history_coord2, 
		//		eff_history_status2,chain_id);
		//}
		//TODO
		//For models that are composed of discrete models (ie model A or B 
		//and min_dim = 0), I just need to continue picking history members 
		//until I get two that match the correct model.If two elements that 
		//contain a given model, we just use a gaussian step and abandon DE 
		//for this step.
		//else{
		//	debugger_print(__FILE__,__LINE__,"ERROR : diff ev is broken with min_dim ==0");
		//}
		
	}
	else if(sampler->RJMCMC){
		debugger_print(__FILE__,__LINE__,"ERROR : diff ev is broken with min_dim ==0");
			
	}
	//Regular PTMCMC
	else{
		//First position ID
		int i = (int)((sampler->history_length)*(gsl_rng_uniform(sampler->rvec[chain_id])));
		//Second position ID
		int j;
		do{
			j=(int)((sampler->history_length)*(gsl_rng_uniform(sampler->rvec[chain_id])));	
		}while(j==i);

		for(int k = 0 ; k <sampler->max_dim; k++){
			eff_history_coord[k] = sampler->history[chain_id][i][k];
			eff_history_coord2[k] = sampler->history[chain_id][j][k];
			eff_history_status[k] = sampler->param_status[chain_id][0][k];
			eff_history_status2[k] = sampler->param_status[chain_id][0][k];
			//eff_history_status[k] = sampler->history_status[chain_id][i][k];
			//eff_history_status2[k] = sampler->history_status[chain_id][j][k];
		}	
			
	}

	
	double alpha = 1;
	double beta = gsl_rng_uniform(sampler->rvec[chain_id]);
	if(beta<.9){
		alpha=gsl_ran_gaussian(sampler->rvec[chain_id],sampler->randgauss_width[chain_id][1][0]);
		//alpha=gsl_ran_gaussian(sampler->rvec[chain_id],2.38/sqrt(2.*sampler->max_dim));
	}
	for (int k = 0; k<sampler->max_dim; k++)
	{
		proposed_param[k] = current_param[k];
		proposed_status[k]=current_status[k];
	}
	for (int k = init_id; k<range+init_id; k++)
	{
//		proposed_param[k] = current_param[k] + alpha*
			//(sampler->history[chain_id][i][k]-sampler->history[chain_id][j][k]);
		if(current_status[k] == 1){
			proposed_param[k] = current_param[k] + alpha *
				(eff_history_coord[k] - eff_history_coord2[k]);
		}
		else{
			proposed_param[k] = 0;
		}
		proposed_status[k]=current_status[k];
	}
	*proposed_model_status = *current_model_status;	
	delete [] eff_history_coord;
	delete [] eff_history_status;
	delete [] eff_history_coord2;
	delete [] eff_history_status2;
}

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
//TODO
/*This could be a problem -- We're not randomly selecting history elements anymore*/
/*This could actually be okay -- We're using most of the correct parameters, just filling in with almost right parameters*/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
void RJ_smooth_history(sampler *sampler, /**<Current sampler */
	double *current_param,/**<Current parameters to match*/
	int *current_param_status,/**<Current parameters to match*/
	int base_history_id, /**<Original history element*/
	double *eff_history_coord, /**<[out] Modified history coord*/
	int *eff_history_status, /**<[out] Modified History status*/
	int chain_id/**<Chain ID of the current chain*/
	)
{
	//Check which dimensions are included in the history element
	int *bad_ids = new int[sampler->max_dim-sampler->min_dim];
	int id_count = 0;
	for (int j = sampler->min_dim ; j<sampler->max_dim; j++){
		//if(current_param_status[j] != sampler->history_status[chain_id][base_history_id][j]){
		if(current_param_status[j] == 1 
			&& sampler->history_status[chain_id][base_history_id][j] != 1){

			bad_ids[id_count] = j;	
			id_count++;
		}
	}
	//Copy over all history params, even if they are zero to ensure no memory errors
	for(int j = 0 ; j<sampler->max_dim; j ++){
		if(current_param_status[j] == 1 ){
			eff_history_coord[j] = sampler->history[chain_id][base_history_id][j];
			eff_history_status[j] = 1;
		}
		else{
			eff_history_coord[j] =0;
			eff_history_status[j] =0;
		}

	}
	//Go back, and overwrite the params not part of the history params
	int min_dev_id = 0;
	double min_dev = 0;
	double current_dev = 0;
	bool at_least_one = false; 
	//For each bad_id, loop through history and find the smallest deviation
	for(int k = 0 ; k<id_count; k++){
		for(int l = 0 ; l<sampler->history_length; l++){
			if(sampler->history_status[chain_id][l][bad_ids[k]]==1){
				for(int j=0 ;j<sampler->min_dim; j++){
						current_dev+= (sampler->history[chain_id][l][j] - sampler->history[chain_id][base_history_id][j])/sampler->history[chain_id][base_history_id][j];
				}
				current_dev/= sampler->min_dim;
				if(!at_least_one){
					min_dev = current_dev;		
					min_dev_id = l;		
					at_least_one = true;
				}
				else{
					if(current_dev<min_dev){
						min_dev = current_dev;		
						min_dev_id = l;		
					}
				}
			}
			current_dev = 0 ;
		}
		if(at_least_one){
			eff_history_coord[k] = sampler->history[chain_id][min_dev_id][k];
		}
		else{
			//Gaussian step for this dimension
			double alpha = gsl_ran_gaussian(sampler->rvec[chain_id],1);
			eff_history_coord[k] = alpha+current_param[k];
		}
		at_least_one = false;
	}
	
	
	delete [] bad_ids;

}

/*! \brief RJ-proposal step for trans-dimensional MCMCs
 *
 * This extra step may seem unnecessary, I'm just adding it in in case the extra flexibility is useful in the future for preprocessing of the chain before sending it to the user's RJ_proposal
 */
void RJ_step(sampler *sampler, /**< sampler*/
	double *current_param, /**< current coordinates in parameter space*/
	double *proposed_param, /**<[out] Proposed coordinates in parameter space*/
	int *current_status, /**< Current status of parameters*/
	int *proposed_status, /**<[out] Proposed status of parameters*/
	int *current_model_status, /**< Current status of parameters*/
	int *proposed_model_status, /**<[out] Proposed status of parameters*/
	double *MH_corrections,
	int chain_number/**< chain mumber*/
	)
{
	sampler->rj(current_param, proposed_param, current_status, proposed_status, current_model_status,proposed_model_status,MH_corrections,sampler->interfaces[chain_number],sampler->user_parameters[chain_number]);
}

/*! \brief subroutine to perform chain comparison for parallel tempering
 *
 * The total output file is passed, and the chains are swapped sequentially
 *
 * This is the routine for ``Deterministic'' sampling (parallel or sequential, 
 * but not pooled)
 */
void chain_swap(sampler *sampler, /**<sampler struct*/
		double ***output, /**<output vector containing chains*/
		int ***param_status,/**< Parameter status*/
		int **model_status,/**< Parameter status*/
		int step_num,	  /**<current step number*/ 
		int *swp_accepted,
		int *swp_rejected
		)
{
	for (int i =0; i < sampler->chain_N-1; i ++)
	{
		int success = single_chain_swap(sampler, output[i][step_num],output[i+1][step_num],
			param_status[i][step_num],param_status[i+1][step_num],&(model_status[i][step_num]),&(model_status[i+1][step_num]), i, i+1);
		if(success==1){
			sampler->swap_accept_ct[i]+=1;
			sampler->swap_accept_ct[i+1]+=1;
			//PT dynamics
			if(sampler->PT_alloc)
				sampler->A[i+1] = 1;
		}
		else{
			sampler->swap_reject_ct[i]+=1;
			sampler->swap_reject_ct[i+1]+=1;
			//PT dynamics
			if(sampler->PT_alloc)
				sampler->A[i+1] = 0;
		}
	
	}
}

/*! \brief subroutine to actually swap two chains
 *
 * This is the more general subroutine, which just swaps the two chains passed to the function
 */
int single_chain_swap(sampler *sampler, /**< sampler structure*/
			double *chain1, /**< parameter position of chain that could be changed*/
			double *chain2, /**< chain that is not swapped, but provides parameters to be swapped by the other chain*/
			int *chain1_status,/**<Parameter status array for chain1*/
			int *chain2_status,/**<Parameter status array for chain2*/
			int *chain1_model_status,/**<Parameter status array for chain1*/
			int *chain2_model_status,/**<Parameter status array for chain2*/
			int T1_index,	/**<number of chain swappe in chain_temps*/
			int T2_index	/**<number of chain swapper in chain_temps*/
			)
{
	//return 0;

	sampler->swap_partners[T1_index][T2_index]++;
	sampler->swap_partners[T2_index][T1_index]++;
	
	//Unpack parameters
	double T1 = sampler->chain_temps[T1_index];
	double T2 = sampler->chain_temps[T2_index];
	//Don't swap same temperature chains
	if(T1==T2){
		return -1;
	}
	double ll1 =  sampler->current_likelihoods[T1_index];
	double ll2 =  sampler->current_likelihoods[T2_index];
	double pow = (ll1-ll2)/T2 - (ll1-ll2)/T1 ;
	double MH_ratio;
	MH_ratio = pow;
	//Averaging the two random numbers from each chains seed
	//double alpha = log( (gsl_rng_uniform(sampler->rvec[T1_index])+gsl_rng_uniform(sampler->rvec[T2_index]))/2.);
	//double alpha = (gsl_rng_uniform(sampler->rvec[T1_index])+gsl_rng_uniform(sampler->rvec[T2_index]))/2.;
	double alpha = gsl_rng_uniform(sampler->rvec[T1_index]);
	MH_ratio = exp(pow);
	if (MH_ratio<alpha)
	{
		return -1;
	}	
	else
	{
		sampler->swap_accepts[T1_index][T2_index]++;
		sampler->swap_accepts[T2_index][T1_index]++;
		double temp[sampler->max_dim];
		int tempstat[sampler->max_dim];
		for(int i =0; i < sampler->max_dim;i++)
		{
			temp[i] = chain1[i];
			tempstat[i] = chain1_status[i];
			chain1[i] = chain2[i];
			chain1_status[i] = chain2_status[i];
			chain2[i]=temp[i];
			chain2_status[i]=tempstat[i];
		}
		double templl = sampler->current_likelihoods[T1_index];
		sampler->current_likelihoods[T1_index] = 
				 sampler->current_likelihoods[T2_index];
		sampler->current_likelihoods[T2_index] = templl;
		int temp_model;
		temp_model = *chain1_model_status;
		*chain1_model_status = *chain2_model_status;
		*chain2_model_status = temp_model;
		return 1;
	}

}

/*! \brief update and initiate probabilities for each variety of step
 *
 * Type 0: Gaussian step
 *
 * Type 1: Differential Evolution step
 *
 * Type 2: MMALA step (currently not supported)
 *
 * Type 3: Fisher step
 */
void assign_probabilities(sampler *sampler, int chain_index)
{

	//non-RJPTMCMC
	if(!sampler->RJMCMC){
		sampler->step_prob[chain_index][4]=.0;
		//no fisher and de not ready
		if(!sampler->fisher_exist && !sampler->de_primed[chain_index])//Obviously must add up to 1
		//if(true)//Obviously must add up to 1
		{
			sampler->step_prob[chain_index][0]=1.;
			sampler->step_prob[chain_index][1]=0.;
			sampler->step_prob[chain_index][2]=0;
			sampler->step_prob[chain_index][3]=0;
			sampler->step_prob[chain_index][5]=0;
		}
		//fisher available, but de not yet ready
		else if (sampler->fisher_exist && !sampler->de_primed[chain_index])
		//else if(true)
		{
			//sampler->step_prob[chain_index][0]=1;
			//sampler->step_prob[chain_index][1]=0;
			//sampler->step_prob[chain_index][2]=.0;
			//sampler->step_prob[chain_index][3]=.0;

			//Tailor to temperature
			sampler->step_prob[chain_index][1]=0;
			sampler->step_prob[chain_index][2]=.0;
			sampler->step_prob[chain_index][3]=.1+.8/sampler->chain_temps[chain_index];
			sampler->step_prob[chain_index][5]=0;
			double sum = sampler->step_prob[chain_index][1]	
				+sampler->step_prob[chain_index][2]
				+sampler->step_prob[chain_index][3];
			sampler->step_prob[chain_index][0]=1-sum;


		}
		//No fisher, but de ready
		else if (!sampler->fisher_exist && sampler->de_primed[chain_index])
		{
			//sampler->step_prob[chain_index][0]=.3;
			//sampler->step_prob[chain_index][1]=.7;
			//sampler->step_prob[chain_index][2]=.0;
			//sampler->step_prob[chain_index][3]=.0;
			//sampler->step_prob[chain_index][5]=0;

			sampler->step_prob[chain_index][0]=1;
			sampler->step_prob[chain_index][1]=0;
			sampler->step_prob[chain_index][2]=.0;
			sampler->step_prob[chain_index][3]=.0;
			sampler->step_prob[chain_index][5]=0;
	

		}
		//all methods available
		else
		{
			//sampler->step_prob[chain_index][0]=.1;
			//sampler->step_prob[chain_index][1]=.4;
			//sampler->step_prob[chain_index][2]=.0;
			//sampler->step_prob[chain_index][3]=.5;
			//Tailor to temperature
			sampler->step_prob[chain_index][1]=.7-.4/sampler->chain_temps[chain_index];
			//sampler->step_prob[chain_index][1]=0;
			sampler->step_prob[chain_index][2]=.0;
			sampler->step_prob[chain_index][3]=.2+.5/sampler->chain_temps[chain_index];
			//sampler->step_prob[chain_index][5]=.1;
			sampler->step_prob[chain_index][5]=.0;
			double sum = sampler->step_prob[chain_index][1]	
				+sampler->step_prob[chain_index][2]
				+sampler->step_prob[chain_index][3]
				+sampler->step_prob[chain_index][5];
			sampler->step_prob[chain_index][0]=1-sum;
			
			//sampler->step_prob[chain_index][0]=.0;
			//sampler->step_prob[chain_index][1]=.0;
			//sampler->step_prob[chain_index][2]=.0;
			//sampler->step_prob[chain_index][3]=.0;
			//sampler->step_prob[chain_index][5]=1;
		}
	}
	else{
		//no fisher and de not ready
		if(!sampler->fisher_exist && !sampler->de_primed[chain_index])//Obviously must add up to 1
		{
			//sampler->step_prob[chain_index][1]=0.;
			//sampler->step_prob[chain_index][2]=0;
			//sampler->step_prob[chain_index][3]=0;
			//sampler->step_prob[chain_index][4]=.2;
	
			sampler->step_prob[chain_index][1]=0.;
			sampler->step_prob[chain_index][2]=0;
			sampler->step_prob[chain_index][3]=0;
			sampler->step_prob[chain_index][4]=.5-.2/sampler->chain_temps[chain_index];
			double sum = sampler->step_prob[chain_index][1]	
				+sampler->step_prob[chain_index][2]
				+sampler->step_prob[chain_index][3]
				+sampler->step_prob[chain_index][4];
			sampler->step_prob[chain_index][0]=1.-sum;
			sampler->step_prob[chain_index][5]=0;
		}
		//fisher available, but de not yet ready
		else if (sampler->fisher_exist && !sampler->de_primed[chain_index])
		{
			//sampler->step_prob[chain_index][0]=.1;
			//sampler->step_prob[chain_index][1]=0;
			//sampler->step_prob[chain_index][2]=.0;
			//sampler->step_prob[chain_index][3]=.7;
			//sampler->step_prob[chain_index][4]=.2;

			sampler->step_prob[chain_index][1]=0.;
			sampler->step_prob[chain_index][2]=0;
			sampler->step_prob[chain_index][3]=.1 +.6/sampler->chain_temps[chain_index];
			sampler->step_prob[chain_index][4]=.7-.5/sampler->chain_temps[chain_index];
			double sum = sampler->step_prob[chain_index][1]	
				+sampler->step_prob[chain_index][2]
				+sampler->step_prob[chain_index][3]
				+sampler->step_prob[chain_index][4];
			sampler->step_prob[chain_index][0]=1-sum;
			sampler->step_prob[chain_index][5]=0;

		}
		//No fisher, but de ready
		else if (!sampler->fisher_exist && sampler->de_primed[chain_index])
		{
			
			//sampler->step_prob[chain_index][0]=.1;
			//sampler->step_prob[chain_index][1]=.7;
			//sampler->step_prob[chain_index][2]=.0;
			//sampler->step_prob[chain_index][3]=.0;
			//sampler->step_prob[chain_index][4]=.2;

			//sampler->step_prob[chain_index][1]=.5+.2/sampler->chain_temps[chain_index];
			//sampler->step_prob[chain_index][2]=0;
			//sampler->step_prob[chain_index][3]=0;
			//sampler->step_prob[chain_index][4]=.4-.2/sampler->chain_temps[chain_index];
			//double sum = sampler->step_prob[chain_index][1]	
			//	+sampler->step_prob[chain_index][2]
			//	+sampler->step_prob[chain_index][3]
			//	+sampler->step_prob[chain_index][4];
			//sampler->step_prob[chain_index][0]=1-sum;
			//sampler->step_prob[chain_index][5]=0;

			sampler->step_prob[chain_index][1]=0;
			sampler->step_prob[chain_index][2]=0;
			sampler->step_prob[chain_index][3]=0;
			sampler->step_prob[chain_index][4]=.5-.2/sampler->chain_temps[chain_index];
			double sum = sampler->step_prob[chain_index][1]	
				+sampler->step_prob[chain_index][2]
				+sampler->step_prob[chain_index][3]
				+sampler->step_prob[chain_index][4];
			sampler->step_prob[chain_index][0]=1-sum;
			sampler->step_prob[chain_index][5]=0;

		}
		//all methods available
		else
		{
			//sampler->step_prob[chain_index][0]=.1;
			//sampler->step_prob[chain_index][1]=.1;
			//sampler->step_prob[chain_index][2]=.0;
			//sampler->step_prob[chain_index][3]=.6;
			//sampler->step_prob[chain_index][4]=.2;

			//sampler->step_prob[chain_index][1]=.2+.1/sampler->chain_temps[chain_index];
			//sampler->step_prob[chain_index][2]=0;
			//sampler->step_prob[chain_index][3]=.1+.4/sampler->chain_temps[chain_index];
			//sampler->step_prob[chain_index][4]=.3-.2/sampler->chain_temps[chain_index];
			//double sum = sampler->step_prob[chain_index][1]	
			//	+sampler->step_prob[chain_index][2]
			//	+sampler->step_prob[chain_index][3]
			//	+sampler->step_prob[chain_index][4];
			//sampler->step_prob[chain_index][0]=1-sum;
			//sampler->step_prob[chain_index][5]=0;

			sampler->step_prob[chain_index][1]=0;
			sampler->step_prob[chain_index][2]=0;
			sampler->step_prob[chain_index][3]=.1+.5/sampler->chain_temps[chain_index];
			sampler->step_prob[chain_index][4]=.5-.4/sampler->chain_temps[chain_index];
			double sum = sampler->step_prob[chain_index][1]	
				+sampler->step_prob[chain_index][2]
				+sampler->step_prob[chain_index][3]
				+sampler->step_prob[chain_index][4];
			sampler->step_prob[chain_index][0]=1-sum;
			sampler->step_prob[chain_index][5]=0;

		}
	}
	//Split probabilities into boundaries for if-else loop
	sampler->prob_boundaries[chain_index][0] = sampler->step_prob[chain_index][0];
	sampler->prob_boundaries[chain_index][1] = sampler->step_prob[chain_index][1]+sampler->prob_boundaries[chain_index][0];
	sampler->prob_boundaries[chain_index][2] = sampler->step_prob[chain_index][2]+sampler->prob_boundaries[chain_index][1];
	sampler->prob_boundaries[chain_index][3] = sampler->step_prob[chain_index][3]+sampler->prob_boundaries[chain_index][2];
	sampler->prob_boundaries[chain_index][4] = sampler->step_prob[chain_index][4]+sampler->prob_boundaries[chain_index][3];
	sampler->prob_boundaries[chain_index][5] = sampler->step_prob[chain_index][5]+sampler->prob_boundaries[chain_index][4];
}	

/*! \brief Copies contents of one chain to another
 *
 * Transfers id_source in samplerptr_source to id_dest samplerptr_dest 
 *
 * NOTE: This copies the VALUE, not the reference. This could be expensive, so use with caution
 *
 * id_dest is ERASED
 *
 * samplerptr_dest and samplerptr_source MUST have the same dimension, the same sampling details (like having or not having a fisher) etc
 * 	
 * samplerptr_dest must be previously allocated properly
 *
 * As output is the largest transfer by far, the transfer_output flag can be used to allow the user to handle that manually.
 */
void transfer_chain(sampler *samplerptr_dest,sampler *samplerptr_source, int id_dest, int id_source, bool transfer_output)
{	
	samplerptr_dest->user_parameters[id_dest] = samplerptr_source->user_parameters[id_source];
	//Position and output and likelihood and temp
	samplerptr_dest->chain_temps[id_dest] = 
		samplerptr_source->chain_temps[id_source];
	samplerptr_dest->chain_pos[id_dest] = samplerptr_source->chain_pos[id_source];
	if(transfer_output){
		for(int i =0 ;i <=samplerptr_dest->chain_pos[id_dest]; i++){
			for(int j =0; j<samplerptr_dest->max_dim; j++){
				samplerptr_dest->output[id_dest][i][j] = 
					samplerptr_source->output[id_source][i][j];
				//Might only need to do this for RJMCMC
				if(samplerptr_dest->RJMCMC){
					samplerptr_dest->param_status[id_dest][i][j] = 
						samplerptr_source->param_status[id_source][i][j];
				}
			}
			samplerptr_dest->model_status[id_dest][i] = 
					samplerptr_source->model_status[id_source][i];
		}
	}
	samplerptr_dest->current_likelihoods[id_dest] = samplerptr_source->current_likelihoods[id_source];
	samplerptr_dest->proper_fisher=samplerptr_source->proper_fisher;

	//Histories
	samplerptr_dest->de_primed[id_dest] 
		= samplerptr_source->de_primed[id_source];
	samplerptr_dest->current_hist_pos[id_dest] 
		= samplerptr_source->current_hist_pos[id_source];
	//If History not filled, only copy over till current pos
	if(!samplerptr_dest->de_primed[id_dest]){
		for(int i =0 ;i <=samplerptr_dest->current_hist_pos[id_dest]; i++){
			for(int j =0; j<samplerptr_dest->max_dim; j++){
				samplerptr_dest->history[id_dest][i][j] = 
					samplerptr_source->history[id_source][i][j];
			}
		}
	}
	//Else, copy over whole history
	else{
		for(int i =0 ;i <samplerptr_dest->history_length; i++){
			for(int j =0; j<samplerptr_dest->max_dim; j++){
				samplerptr_dest->history[id_dest][i][j] = 
					samplerptr_source->history[id_source][i][j];
			}
		}
	}
	//KDE
	if(samplerptr_source->kde_cov_update_ct[id_source] != samplerptr_source->kde_cov_update_number){
		for (int i =0 ; i<samplerptr_source->max_dim; i++){
			for (int j =0 ; j<samplerptr_source->max_dim; j++){
				samplerptr_dest->kde_cov[id_dest][i][j] = 
					samplerptr_source->kde_cov[id_source][i][j];
				samplerptr_dest->kde_fisher[id_dest][i][j] = 
					samplerptr_source->kde_fisher[id_source][i][j];
			}
		}
		samplerptr_dest->kde_cov_update_ct[id_dest] = samplerptr_source->kde_cov_update_ct[id_source];
		samplerptr_dest->kde_cov_lndet[id_dest] = samplerptr_source->kde_cov_lndet[id_source];
	}
	else{
		samplerptr_dest->kde_cov_update_ct[id_dest] = samplerptr_dest->kde_cov_update_number;
		if(samplerptr_dest->de_primed[id_dest] && samplerptr_dest->kde_step){
			update_kde_cov(samplerptr_dest, id_dest);	
		}
	}

	//Step parameters
	samplerptr_dest->gauss_last_accept_ct[id_dest] 
		= samplerptr_source->gauss_last_accept_ct[id_source]; 
	samplerptr_dest->gauss_last_reject_ct[id_dest] 
		= samplerptr_source->gauss_last_reject_ct[id_source]; 

	for(int i = 0 ; i<samplerptr_dest->max_dim; i++){
		samplerptr_dest->gauss_last_accept_ct_per_dim[id_dest][i]
			= samplerptr_source->gauss_last_accept_ct_per_dim[id_source][i]; 
		samplerptr_dest->gauss_last_reject_ct_per_dim[id_dest][i]
			= samplerptr_source->gauss_last_reject_ct_per_dim[id_source][i]; 

		samplerptr_dest->gauss_accept_ct_per_dim[id_dest][i]
			= samplerptr_source->gauss_accept_ct_per_dim[id_source][i]; 
		samplerptr_dest->gauss_reject_ct_per_dim[id_dest][i] 
			= samplerptr_source->gauss_reject_ct_per_dim[id_source][i]; 
	}

	samplerptr_dest->de_last_accept_ct[id_dest] 
		= samplerptr_source->de_last_accept_ct[id_source]; 
	samplerptr_dest->de_last_reject_ct[id_dest] 
		= samplerptr_source->de_last_reject_ct[id_source]; 
	samplerptr_dest->fish_last_accept_ct[id_dest] 
		= samplerptr_source->fish_last_accept_ct[id_source]; 
	samplerptr_dest->fish_last_reject_ct[id_dest] 
		= samplerptr_source->fish_last_reject_ct[id_source]; 

	samplerptr_dest->kde_last_accept_ct[id_dest] 
		= samplerptr_source->kde_last_accept_ct[id_source]; 
	samplerptr_dest->kde_last_reject_ct[id_dest] 
		= samplerptr_source->kde_last_reject_ct[id_source]; 

	samplerptr_dest->gauss_accept_ct[id_dest] 
		= samplerptr_source->gauss_accept_ct[id_source]; 
	samplerptr_dest->gauss_reject_ct[id_dest] 
		= samplerptr_source->gauss_reject_ct[id_source]; 

	samplerptr_dest->de_accept_ct[id_dest] 
		= samplerptr_source->de_accept_ct[id_source]; 
	samplerptr_dest->de_reject_ct[id_dest] 
		= samplerptr_source->de_reject_ct[id_source]; 

	samplerptr_dest->fish_accept_ct[id_dest] 
		= samplerptr_source->fish_accept_ct[id_source]; 
	samplerptr_dest->fish_reject_ct[id_dest] 
		= samplerptr_source->fish_reject_ct[id_source]; 

	samplerptr_dest->kde_accept_ct[id_dest] 
		= samplerptr_source->kde_accept_ct[id_source]; 
	samplerptr_dest->kde_reject_ct[id_dest] 
		= samplerptr_source->kde_reject_ct[id_source]; 

	samplerptr_dest->mmala_accept_ct[id_dest] 
		= samplerptr_source->mmala_accept_ct[id_source]; 
	samplerptr_dest->mmala_reject_ct[id_dest] 
		= samplerptr_source->mmala_reject_ct[id_source]; 
	for(int i =0 ;i<samplerptr_dest->types_of_steps; i++){
		samplerptr_dest->step_prob[id_dest][i] 
			= samplerptr_source->step_prob[id_source][i];
		samplerptr_dest->prob_boundaries[id_dest][i] 
			= samplerptr_source->prob_boundaries[id_source][i];
		for(int j = 0 ; j<samplerptr_dest->randgauss_width_number[id_source][i];j++){
			samplerptr_dest->randgauss_width[id_dest][i][j]
				= samplerptr_source->randgauss_width[id_source][i][j];
		}
	}
	//Fisher
	if(samplerptr_dest->fisher_exist){
		if(samplerptr_source->fisher_update_ct[id_source] != samplerptr_source->fisher_update_number){
			for (int i =0 ; i<samplerptr_source->max_dim; i++){
				for (int j =0 ; j<samplerptr_source->max_dim; j++){
					samplerptr_dest->fisher_vecs[id_dest][i][j] = 
						samplerptr_source->fisher_vecs[id_source][i][j];
					samplerptr_dest->fisher_vecs_prop[id_dest][i][j] = 
						samplerptr_source->fisher_vecs_prop[id_source][i][j];
					samplerptr_dest->fisher_matrix[id_dest][i][j]=samplerptr_source->fisher_matrix[id_source][i][j];
					samplerptr_dest->fisher_matrix_prop[id_dest][i][j]=samplerptr_source->fisher_matrix_prop[id_source][i][j];
				}
				samplerptr_dest->fisher_vals[id_dest][i] = 
					samplerptr_source->fisher_vals[id_source][i];
				samplerptr_dest->fisher_vals_prop[id_dest][i] = 
					samplerptr_source->fisher_vals_prop[id_source][i];
			}
			samplerptr_dest->fisher_update_ct[id_dest] = samplerptr_source->fisher_update_ct[id_source];
		}
		else{
			samplerptr_dest->fisher_update_ct[id_dest] = samplerptr_dest->fisher_update_number;
			update_fisher(samplerptr_dest, samplerptr_dest->output[id_dest][samplerptr_dest->chain_pos[id_dest]], samplerptr_dest->param_status[id_dest][samplerptr_dest->chain_pos[id_dest]],&(samplerptr_dest->model_status[id_dest][samplerptr_dest->chain_pos[id_dest]]),id_dest);	
			if(samplerptr_dest->proper_fisher){
				iterate_fisher(samplerptr_dest, id_dest);
			}
		}
	}



	//Sampling parameters
	samplerptr_dest->waiting[id_dest] = samplerptr_source->waiting[id_source];
	samplerptr_dest->waiting_SWP[id_dest] = samplerptr_source->waiting_SWP[id_source];
	samplerptr_dest->priority[id_dest] = samplerptr_source->priority[id_source];
	samplerptr_dest->ref_chain_status[id_dest] = samplerptr_source->ref_chain_status[id_source];
	samplerptr_dest->nan_counter[id_dest] = samplerptr_source->nan_counter[id_source];
	samplerptr_dest->num_gauss[id_dest] = samplerptr_source->num_gauss[id_source];
	samplerptr_dest->num_fish[id_dest] = samplerptr_source->num_fish[id_source];
	samplerptr_dest->num_de[id_dest] = samplerptr_source->num_de[id_source];
	samplerptr_dest->num_mmala[id_dest] = samplerptr_source->num_mmala[id_source];
	samplerptr_dest->swap_accept_ct[id_dest] = samplerptr_source->swap_accept_ct[id_source];
	samplerptr_dest->swap_reject_ct[id_dest] = samplerptr_source->swap_reject_ct[id_source];
	samplerptr_dest->step_accept_ct[id_dest] = samplerptr_source->step_accept_ct[id_source];
	samplerptr_dest->step_reject_ct[id_dest] = samplerptr_source->step_reject_ct[id_source];
	if(samplerptr_source->log_ll){
		for(int i =0 ; i<samplerptr_source->chain_pos[id_source]; i++)
			samplerptr_dest->ll_lp_output[id_dest][i][0] 
				= samplerptr_source->ll_lp_output[id_source][i][0];
	}
	if(samplerptr_source->log_lp){
		for(int i =0 ; i<samplerptr_source->chain_pos[id_source]; i++)
			samplerptr_dest->ll_lp_output[id_dest][i][1] 
				= samplerptr_source->ll_lp_output[id_source][i][1];
	}
	if(samplerptr_source->PT_alloc)
		samplerptr_dest->A[id_dest] = samplerptr_source->A[id_source];

	samplerptr_dest->prop_MH_factor[id_dest] = samplerptr_source->prop_MH_factor[id_source];
	samplerptr_dest->max_target_accept_ratio[id_dest] = samplerptr_source->max_target_accept_ratio[id_source];
	samplerptr_dest->min_target_accept_ratio[id_dest] = samplerptr_source->min_target_accept_ratio[id_source];
	samplerptr_dest->min_target_accept_ratio[id_dest] = samplerptr_source->min_target_accept_ratio[id_source];
	if(samplerptr_source->update_RJ_width){
		samplerptr_dest->RJstep_accept_ct[id_dest] = samplerptr_source->RJstep_accept_ct[id_source];
		samplerptr_dest->RJstep_reject_ct[id_dest] = samplerptr_source->RJstep_reject_ct[id_source];
	}

}

/*! \brief Checks the status of a sampler for the stochastic sampling
 *
 * Just loops through the ref_chain_status variables
 */
bool check_sampler_status(sampler *samplerptr)
{
	for(int i =0 ; i<samplerptr->chain_N; i++){
		if(!samplerptr->ref_chain_status[i])
			return false;
	}
	return true;
}
/*! \brief Updates the step widths, shooting for 20% acceptance ratios for each type of step
 */
void update_step_widths(sampler *samplerptr, int chain_id)
{
	//return;
	int j = chain_id;
	//update stepsize to maximize step efficiency
	//increases in stepsizes of 10%
	double frac, acc, rej;
	if(samplerptr->chain_pos[j]%samplerptr->check_stepsize_freq[j] == 0){
		//Gaussian
		if(samplerptr->step_prob[j][0]!= 0){
			for(int i =0 ; i<samplerptr->max_dim ; i++){
				acc = samplerptr->gauss_accept_ct_per_dim[j][i] - samplerptr->gauss_last_accept_ct_per_dim[j][i];	
				rej = samplerptr->gauss_reject_ct_per_dim[j][i] - samplerptr->gauss_last_reject_ct_per_dim[j][i];	
				frac = acc / (acc + rej);
				if(frac<samplerptr->min_target_accept_ratio[j]){
					samplerptr->randgauss_width[j][0][i] *=.9;	
				}
				else if(frac>samplerptr->max_target_accept_ratio[j]){
					samplerptr->randgauss_width[j][0][i] *=1.1;	
				}
				samplerptr->gauss_last_accept_ct_per_dim[j][i]=samplerptr->gauss_accept_ct_per_dim[j][i];
				samplerptr->gauss_last_reject_ct_per_dim[j][i]=samplerptr->gauss_reject_ct_per_dim[j][i];
			}
		}	
		//de
		if(samplerptr->step_prob[j][1]!= 0){
			acc = samplerptr->de_accept_ct[j] - samplerptr->de_last_accept_ct[j];	
			rej = samplerptr->de_reject_ct[j] - samplerptr->de_last_reject_ct[j];	
			frac = acc / (acc + rej);
			if(frac<samplerptr->min_target_accept_ratio[j]){
				samplerptr->randgauss_width[j][1][0] *=.9;	
			}
			else if(frac>samplerptr->max_target_accept_ratio[j]){
				samplerptr->randgauss_width[j][1][0] *=1.1;	
			}
			samplerptr->de_last_accept_ct[j]=samplerptr->de_accept_ct[j];
			samplerptr->de_last_reject_ct[j]=samplerptr->de_reject_ct[j];
		}	
		//fisher
		if(samplerptr->step_prob[j][3]!= 0){
			acc = samplerptr->fish_accept_ct[j] - samplerptr->fish_last_accept_ct[j];	
			rej = samplerptr->fish_reject_ct[j] - samplerptr->fish_last_reject_ct[j];	
			frac = acc / (acc + rej);
			if(frac<samplerptr->min_target_accept_ratio[j]){
				samplerptr->randgauss_width[j][3][0] *=.9;	
			}
			else if(frac>samplerptr->max_target_accept_ratio[j]){
				samplerptr->randgauss_width[j][3][0] *=1.1;	
			}
			samplerptr->fish_last_accept_ct[j]=samplerptr->fish_accept_ct[j];
			samplerptr->fish_last_reject_ct[j]=samplerptr->fish_reject_ct[j];
		}	
		if(samplerptr->step_prob[j][5]!= 0){
			acc = samplerptr->kde_accept_ct[j] - samplerptr->kde_last_accept_ct[j];	
			rej = samplerptr->kde_reject_ct[j] - samplerptr->kde_last_reject_ct[j];	
			frac = acc / (acc + rej);
			if(frac<samplerptr->min_target_accept_ratio[j]){
				samplerptr->randgauss_width[j][5][0] *=.9;	
			}
			else if(frac>samplerptr->max_target_accept_ratio[j]){
				samplerptr->randgauss_width[j][5][0] *=1.1;	
			}
			samplerptr->kde_last_accept_ct[j]=samplerptr->kde_accept_ct[j];
			samplerptr->kde_last_reject_ct[j]=samplerptr->kde_reject_ct[j];
		}	
		//RJ
		if(samplerptr->update_RJ_width){
			if(samplerptr->step_prob[j][4]!= 0){
				acc = samplerptr->RJstep_accept_ct[j] - samplerptr->RJstep_last_accept_ct[j];	
				rej = samplerptr->RJstep_reject_ct[j] - samplerptr->RJstep_last_reject_ct[j];	
				frac = acc / (acc + rej);
				if(frac<samplerptr->min_target_accept_ratio[j]){
					samplerptr->randgauss_width[j][4][0] *=.9;	
				}
				else if(frac>samplerptr->max_target_accept_ratio[j]){
					samplerptr->randgauss_width[j][4][0] *=1.1;	
				}
				samplerptr->RJstep_last_accept_ct[j]=samplerptr->RJstep_accept_ct[j];
				samplerptr->RJstep_last_reject_ct[j]=samplerptr->RJstep_reject_ct[j];
			}	
		}
		
	}
}

void allocate_sampler_mem(sampler *sampler)
{
	gsl_rng_env_setup();
	const gsl_rng_type *T= gsl_rng_default;

	int i;
	if (!sampler->user_parameters){
		sampler->local_param_allocation=true;
		sampler->user_parameters = new void*[sampler->chain_N];
	}
	sampler->step_prob = (double **)malloc(sizeof(double *) * sampler->chain_N);
	sampler->prob_boundaries = (double **)malloc(sizeof(double *) * sampler->chain_N);
	sampler->de_primed = (bool *)malloc(sizeof(bool ) * sampler->chain_N);
	sampler->waiting = (bool *)malloc(sizeof(bool ) * sampler->chain_N);
	sampler->waiting_SWP = (bool *)malloc(sizeof(bool ) * sampler->chain_N);
	sampler->restarted_chain = (bool *)malloc(sizeof(bool ) * sampler->chain_N);
	sampler->queue_mutexes = new std::mutex[sampler->chain_N];
	sampler->current_hist_pos = (int *)malloc(sizeof(int ) * sampler->chain_N);
	sampler->chain_pos = (int *)malloc(sizeof(int ) * sampler->chain_N);

	sampler->kde_accept_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->kde_reject_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->fish_accept_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->fish_reject_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->de_accept_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->de_reject_ct = (int *)malloc(sizeof(int) * sampler->chain_N);

	sampler->gauss_accept_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->gauss_reject_ct = (int *)malloc(sizeof(int) * sampler->chain_N);

	sampler->gauss_accept_ct_per_dim = new int*[sampler->chain_N];
	sampler->gauss_reject_ct_per_dim = new int*[sampler->chain_N];

	sampler->mmala_accept_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->mmala_reject_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->RJstep_accept_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->RJstep_reject_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->fisher_update_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->kde_cov_update_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->swap_accept_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->swap_reject_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->step_accept_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->step_reject_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->rvec = (gsl_rng **)malloc(sizeof(gsl_rng *) * sampler->chain_N);

	sampler->nan_counter = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->num_gauss = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->num_fish = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->num_de = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->num_mmala = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->num_RJstep = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->num_kde = (int *)malloc(sizeof(int) * sampler->chain_N);

	sampler->priority = (int *)malloc(sizeof(int) * sampler->chain_N);

	sampler->current_likelihoods = (double *)malloc(sizeof(double) * sampler->chain_N);

	sampler->check_stepsize_freq = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->max_target_accept_ratio = (double *)malloc(sizeof(double) * sampler->chain_N);
	sampler->min_target_accept_ratio = (double *)malloc(sizeof(double) * sampler->chain_N);

	sampler->gauss_last_accept_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->gauss_last_reject_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->gauss_last_accept_ct_per_dim = new int*[ sampler->chain_N];
	sampler->gauss_last_reject_ct_per_dim = new int*[ sampler->chain_N];

	sampler->kde_last_accept_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->kde_last_reject_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->fish_last_accept_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->fish_last_reject_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->de_last_accept_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->de_last_reject_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->swap_partners = new int*[sampler->chain_N];
	sampler->swap_accepts = new int*[sampler->chain_N];
	sampler->thermodynamic_integrated_likelihood = new double[sampler->chain_N];
	sampler->thermodynamic_integrated_likelihood_terms = new int[sampler->chain_N];
	if(sampler->update_RJ_width){
		sampler->RJstep_last_accept_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
		sampler->RJstep_last_reject_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	}
	//sampler->randgauss_width = allocate_2D_array(sampler->chain_N, sampler->types_of_steps); //Second dimension is types of steps
	sampler->randgauss_width = new double **[sampler->chain_N]; //Second dimension is types of steps
	sampler->randgauss_width_number = new int *[sampler->chain_N]; //Second dimension is types of steps
	sampler->prop_MH_factor = (double *)malloc(sizeof(double) *sampler->chain_N);

	//RJ parameters -- initialize status array with 1's for now, then repopulate with initial position
	if(sampler->RJMCMC){
		//sampler->param_status = allocate_3D_array_int(sampler->chain_N, sampler->N_steps,sampler->max_dim);	
		sampler->history_status = allocate_3D_array_int(sampler->chain_N, sampler->history_length,sampler->max_dim);	
	}
	else{
		//Need one copy per chain at least -- just in case of race conditions
		sampler->param_status = (int ***)malloc(sizeof(int **)*sampler->chain_N);
		sampler->history_status = (int ***)malloc(sizeof(int **)*sampler->chain_N);
	}
	sampler->model_status = (int **)malloc(sizeof(int *)*sampler->chain_N);


	sampler->ref_chain_status = (bool *) malloc(sizeof(bool)*sampler->chain_N);
	sampler->ref_chain_ids = (int *) malloc(sizeof(int)*sampler->chain_N);

	sampler->ref_chain_num = 0 ;

	//#############
	int cold_chain_ids[sampler->chain_N];
	int cold_chain_ct = 0;
	//if(sampler->pool){
	sampler->chain_neighborhoods=(double **)malloc(sizeof(double*)*sampler->chain_N);
	sampler->chain_neighbors=(int *)malloc(sizeof(int)*sampler->chain_N);
	sampler->chain_neighborhoods_ids=(int **)malloc(sizeof(int*)*sampler->chain_N);
	sampler->chain_neighbors_ids=(int *)malloc(sizeof(int)*sampler->chain_N);
	for(int i = 0  ;i<sampler->chain_N; i++){
		sampler->chain_neighborhoods[i]=NULL;	
		sampler->chain_neighborhoods_ids[i]=NULL;	
	}
	update_temp_neighborhoods(sampler);
		
	//#############
	sampler->interfaces = new mcmc_data_interface*[sampler->chain_N];
	//#############
	if(sampler->tune){
		sampler->fisher_update_number = 200;
		sampler->kde_cov_update_number = 200;
		//sampler->fisher_update_number = 1000;
		//sampler->fisher_update_number = 50;
	}
	else{
		//sampler->fisher_update_number = 50;
		//sampler->fisher_update_number = 1000;
		//sampler->fisher_update_number = 5000;
		sampler->fisher_update_number = sampler->N_steps;
		sampler->kde_cov_update_number = sampler->N_steps;
	}
	//#############
	for (i =0; i<sampler->chain_N; i++)
	{
		sampler->swap_partners[i] = new int[sampler->chain_N];
		sampler->swap_accepts[i] = new int[sampler->chain_N];
		for(int j = 0 ; j<sampler->chain_N;j++){
			sampler->swap_partners[i][j] = 0 ;
			sampler->swap_accepts[i][j] = 0 ;
		}
		//designate T=1 chains as reference chains
		if(fabs(sampler->chain_temps[i] - 1)<DOUBLE_COMP_THRESH){
			sampler->ref_chain_status[i] = false;
			sampler->ref_chain_ids[sampler->ref_chain_num]=i;
			sampler->ref_chain_num++;
				
		}
		else{
			sampler->ref_chain_status[i] = true;
		}
		sampler->interfaces[i] = new mcmc_data_interface;
		sampler->interfaces[i]->min_dim = sampler->min_dim;
		sampler->interfaces[i]->max_dim = sampler->max_dim;
		sampler->interfaces[i]->chain_number = sampler->chain_N;
		sampler->interfaces[i]->chain_id = i;
		sampler->interfaces[i]->nested_model_number = sampler->nested_model_number;
		sampler->interfaces[i]->burn_phase = sampler->burn_phase;
		//#######################
		//#######################
		//initial value set to one's
		if(sampler->RJMCMC){
		}
		else{
			//debugger_print(__FILE__,__LINE__,std::to_string(sampler->N_steps)+" "+std::to_string(sampler->max_dim));
			sampler->param_status[i] = (int **)malloc(sizeof(int *)*sampler->N_steps);
			sampler->param_status[i][0] = (int *)malloc(sizeof(int )*sampler->max_dim);
			for(int k = 0 ; k<sampler->max_dim; k++){
				sampler->param_status[i][0][k] = 1;
			}
			for(int k = 1 ; k<sampler->N_steps; k ++){
				sampler->param_status[i][k] = sampler->param_status[i][0];
			}
			sampler->history_status[i] = (int **)malloc(sizeof(int *)*sampler->history_length);
			sampler->history_status[i][0] = (int *)malloc(sizeof(int )*sampler->max_dim);
			for(int k = 0 ; k<sampler->max_dim; k++){
				sampler->history_status[i][0][k] = 1;
			}
			for(int k = 1 ; k<sampler->history_length; k ++){
				sampler->history_status[i][k] = sampler->history_status[i][0];
			}

		}
		//TODO
		sampler->model_status[i] = (int *)malloc(sizeof(int *)*sampler->N_steps);
		for(int j = 0 ; j<sampler->N_steps;j++){
			sampler->model_status[i][j] = 0 ;
		}

		sampler->step_prob[i] = (double *)malloc(sizeof(double)*sampler->types_of_steps);
		sampler->prob_boundaries[i] = (double *)malloc(sizeof(double)*sampler->types_of_steps);
		sampler->de_primed[i] = false;
		sampler->current_hist_pos[i] = 0;
		sampler->kde_accept_ct[i]=0;
		sampler->kde_reject_ct[i]=0;
		sampler->fish_accept_ct[i]=0;
		sampler->fish_reject_ct[i]=0;
		sampler->de_accept_ct[i]=0;
		sampler->de_reject_ct[i]=0;
		sampler->gauss_accept_ct[i]=0;
		sampler->gauss_reject_ct[i]=0;
		sampler->mmala_accept_ct[i]=0;
		sampler->mmala_reject_ct[i]=0;
		sampler->RJstep_accept_ct[i]=0;
		sampler->RJstep_reject_ct[i]=0;
		sampler->chain_pos[i]=0;
		sampler->waiting[i]=true;
		sampler->waiting_SWP[i]=false;
		sampler->swap_accept_ct[i]=0;
		sampler->swap_reject_ct[i]=0;
		sampler->step_accept_ct[i]=0;
		sampler->step_reject_ct[i]=0;

		sampler->nan_counter[i]=0;
		sampler->num_gauss[i]=0;
		sampler->num_fish[i]=0;
		sampler->num_de[i]=0;
		sampler->num_mmala[i]=0;
		sampler->num_RJstep[i]=0;
		sampler->num_kde[i]=0;
	
		sampler->priority[i] = 1; //Default priority
		sampler->restarted_chain[i]=false;
		sampler->thermodynamic_integrated_likelihood[i]=0;
		sampler->thermodynamic_integrated_likelihood_terms[i]=0;

		sampler->rvec[i] = gsl_rng_alloc(T);

		//Seed differently
		gsl_rng_set(sampler->rvec[i] , i+1);
	
		if(sampler->tune){
			//sampler->check_stepsize_freq[i] = 500;
			sampler->check_stepsize_freq[i] = 50;
		}
		else{
			sampler->check_stepsize_freq[i] = sampler->N_steps;
		}
		//max probability is a function of the temperature -- higher temp 
		//are allowed to step more
		//sampler->max_target_accept_ratio[i] = .90-.15/sampler->chain_temps[i];
		//sampler->min_target_accept_ratio[i] = .60;
		sampler->max_target_accept_ratio[i] = .60-.2/sampler->chain_temps[i];
		sampler->min_target_accept_ratio[i] = .2;
		sampler->gauss_last_accept_ct[i] = 0.;
		sampler->gauss_last_reject_ct[i] = 0.;

		sampler->gauss_last_accept_ct_per_dim[i] = new int[sampler->max_dim];
		sampler->gauss_last_reject_ct_per_dim[i] = new int[sampler->max_dim];

		sampler->gauss_accept_ct_per_dim[i] = new int[sampler->max_dim];
		sampler->gauss_reject_ct_per_dim[i] = new int[sampler->max_dim];
		for(int j = 0 ; j<sampler->max_dim; j++){
			sampler->gauss_last_accept_ct_per_dim[i][j]=0;
			sampler->gauss_last_reject_ct_per_dim[i][j]=0;
			sampler->gauss_accept_ct_per_dim[i][j]=0;
			sampler->gauss_reject_ct_per_dim[i][j]=0;
		}

		sampler->kde_last_accept_ct[i] = 0.;
		sampler->kde_last_reject_ct[i] = 0.;
		sampler->fish_last_accept_ct[i] = 0.;
		sampler->fish_last_reject_ct[i] = 0.;
		sampler->de_last_accept_ct[i] = 0.;
		sampler->de_last_reject_ct[i] = 0.;
		if(sampler->update_RJ_width){
			sampler->RJstep_last_accept_ct[i] = 0.;
			sampler->RJstep_last_reject_ct[i] = 0.;
		}

		//Initial width size for all chains, all steps is 1.
		sampler->randgauss_width[i] = new double*[sampler->types_of_steps];
		sampler->randgauss_width_number[i] = new int[sampler->types_of_steps];
		sampler->randgauss_width_number[i][0] = sampler->max_dim;
		sampler->randgauss_width_number[i][1] = 1;
		sampler->randgauss_width_number[i][2] = 1;
		sampler->randgauss_width_number[i][3] = 1;
		sampler->randgauss_width_number[i][4] = 1;
		sampler->randgauss_width_number[i][5] = 1;
		sampler->randgauss_width[i][0] = new double[sampler->max_dim];
		sampler->randgauss_width[i][1] = new double[1];
		sampler->randgauss_width[i][2] = new double[1];
		sampler->randgauss_width[i][3] = new double[1];
		sampler->randgauss_width[i][4] = new double[1];
		sampler->randgauss_width[i][5] = new double[1];
		
		//sampler->randgauss_width[i][0]=.01;
		for(int j = 0 ; j<sampler->max_dim; j++){
			sampler->randgauss_width[i][0][j]=.05;
		}
		//sampler->randgauss_width[i][1]=.05;
		sampler->randgauss_width[i][1][0]=1;
		sampler->randgauss_width[i][2][0]=.05;
		sampler->randgauss_width[i][3][0]=.5;
		//For RJPTMCMC, this may not be used, but it'll be available
		sampler->randgauss_width[i][4][0]=.5;
		sampler->randgauss_width[i][5][0]=.5;

		sampler->prop_MH_factor[i]=0;

		sampler->fisher_update_ct[i] = sampler->fisher_update_number;
		sampler->kde_cov_update_ct[i] = sampler->kde_cov_update_number;
	}		
	sampler->history = allocate_3D_array(sampler->chain_N, 
				sampler->history_length, sampler->max_dim);
	sampler->fisher_vecs = allocate_3D_array(sampler->chain_N, 
				sampler->max_dim, sampler->max_dim);
	sampler->fisher_vals = allocate_2D_array(sampler->chain_N, sampler->max_dim);
	sampler->fisher_vecs_prop = allocate_3D_array(sampler->chain_N, 
				sampler->max_dim, sampler->max_dim);
	sampler->fisher_vals_prop = allocate_2D_array(sampler->chain_N, sampler->max_dim);
	sampler->fisher_matrix = allocate_3D_array(sampler->chain_N, 
				sampler->max_dim, sampler->max_dim);
	sampler->fisher_matrix_prop = allocate_3D_array(sampler->chain_N, sampler->max_dim,sampler->max_dim);
	sampler->kde_cov = allocate_3D_array(sampler->chain_N, 
				sampler->max_dim, sampler->max_dim);
	sampler->kde_fisher = allocate_3D_array(sampler->chain_N, 
				sampler->max_dim, sampler->max_dim);
	sampler->kde_cov_lndet = new double[sampler->chain_N];
	

	//Trouble Shooting:
	if(sampler->log_ll || sampler->log_lp){
		sampler->ll_lp_output = allocate_3D_array(sampler->chain_N, 
			sampler->N_steps, 2);
	}
}

void deallocate_sampler_mem(sampler *sampler)
{
	int i;
	for (i =0; i<sampler->chain_N; i++)
	{
		free(sampler->prob_boundaries[i]); 
		gsl_rng_free(sampler->rvec[i]);
		delete [] sampler->swap_partners[i];
		delete [] sampler->swap_accepts[i];

	}		
	delete [] sampler->swap_partners;
	delete [] sampler->swap_accepts;
	//if(sampler->pool){
	for (i =0; i<sampler->chain_N; i++)
	{
		free(sampler->chain_neighborhoods[i]);
		free(sampler->chain_neighborhoods_ids[i]);
		free(sampler->step_prob[i]);
	}
	free(sampler->chain_neighborhoods);
	free(sampler->chain_neighbors);
	free(sampler->chain_neighborhoods_ids);
	free(sampler->chain_neighbors_ids);

	if(sampler->local_param_allocation){
		delete [] sampler->user_parameters;
	}
	for(int j =0 ; j<sampler->chain_N; j++){
		delete sampler->interfaces[j];
	}
	delete [] sampler->interfaces;
	free(sampler->step_prob); 
	free(sampler->prob_boundaries); 
	free(sampler->de_primed);
	free(sampler->current_hist_pos);
	free(sampler->kde_accept_ct);
	free(sampler->kde_reject_ct);
	free(sampler->fish_accept_ct);
	free(sampler->fish_reject_ct);
	free(sampler->de_accept_ct);
	free(sampler->de_reject_ct);
	free(sampler->gauss_accept_ct);
	free(sampler->gauss_reject_ct);
	free(sampler->mmala_accept_ct);
	free(sampler->mmala_reject_ct);
	free(sampler->RJstep_accept_ct);
	free(sampler->RJstep_reject_ct);
	free(sampler->chain_pos);
	free(sampler->waiting);
	free(sampler->waiting_SWP);
	delete [] sampler->queue_mutexes;
	free(sampler->swap_accept_ct);
	free(sampler->swap_reject_ct);
	free(sampler->step_accept_ct);
	free(sampler->step_reject_ct);

	free(sampler->num_kde);
	free(sampler->nan_counter);
	free(sampler->num_gauss);
	free(sampler->num_fish);
	free(sampler->num_de);
	free(sampler->num_mmala);
	free(sampler->num_RJstep);

	free(sampler->current_likelihoods);
	free(sampler->priority);
	free(sampler->restarted_chain);
	delete [] sampler->thermodynamic_integrated_likelihood;
	delete [] sampler->thermodynamic_integrated_likelihood_terms;

	deallocate_3D_array(sampler->history,sampler->chain_N, 
				sampler->history_length, sampler->max_dim);
	deallocate_3D_array(sampler->fisher_vecs, sampler->chain_N, sampler->max_dim, sampler->max_dim);
	deallocate_2D_array(sampler->fisher_vals, sampler->chain_N, sampler->max_dim);
	deallocate_3D_array(sampler->fisher_vecs_prop, sampler->chain_N, sampler->max_dim, sampler->max_dim);
	deallocate_2D_array(sampler->fisher_vals_prop, sampler->chain_N, sampler->max_dim);
	deallocate_3D_array(sampler->fisher_matrix, sampler->chain_N, sampler->max_dim,sampler->max_dim);
	deallocate_3D_array(sampler->fisher_matrix_prop, sampler->chain_N, sampler->max_dim,sampler->max_dim);
	deallocate_3D_array(sampler->kde_cov, sampler->chain_N, sampler->max_dim, sampler->max_dim);
	deallocate_3D_array(sampler->kde_fisher, sampler->chain_N, sampler->max_dim, sampler->max_dim);
	delete [] sampler->kde_cov_lndet;
 
	free(sampler->fisher_update_ct);
	free(sampler->kde_cov_update_ct);
	free(sampler->rvec);

	free(sampler->check_stepsize_freq);
	free(sampler->gauss_last_accept_ct);
	free(sampler->gauss_last_reject_ct);
	for(int i = 0 ; i < sampler->chain_N ; i++){
		delete [] sampler->gauss_last_reject_ct_per_dim[i];
		delete [] sampler->gauss_last_accept_ct_per_dim[i];
		delete [] sampler->gauss_reject_ct_per_dim[i];
		delete [] sampler->gauss_accept_ct_per_dim[i];
	}
	delete [] sampler->gauss_last_reject_ct_per_dim;
	delete [] sampler->gauss_last_accept_ct_per_dim;
	delete [] sampler->gauss_reject_ct_per_dim;
	delete [] sampler->gauss_accept_ct_per_dim;

	free(sampler->kde_last_accept_ct);
	free(sampler->kde_last_reject_ct);
	free(sampler->de_last_accept_ct);
	free(sampler->de_last_reject_ct);
	free(sampler->fish_last_accept_ct);
	free(sampler->fish_last_reject_ct);
	if(sampler->update_RJ_width){
		free(sampler->RJstep_last_accept_ct);
		free(sampler->RJstep_last_reject_ct);
	}
	free(sampler->max_target_accept_ratio);
	free(sampler->min_target_accept_ratio);
	//deallocate_2D_array(sampler->randgauss_width,sampler->chain_N, sampler->types_of_steps);
	for(int i = 0 ; i<sampler->chain_N; i++){
		for (int j = 0 ; j<sampler->types_of_steps; j++){
			delete [] sampler->randgauss_width[i][j];
		}
		delete [] sampler->randgauss_width[i];
		delete [] sampler->randgauss_width_number[i];
	}
	delete [] sampler->randgauss_width;
	delete [] sampler->randgauss_width_number;
	if(sampler->RJMCMC){
		deallocate_3D_array(sampler->history_status,sampler->chain_N, sampler->history_length, sampler->max_dim);
	}
	else{
		for(int j = 0; j<sampler->chain_N; j++){
			free(sampler->param_status[j][0]);
			free(sampler->param_status[j]);
		}
		free(sampler->param_status);
		for(int j = 0; j<sampler->chain_N; j++){
			free(sampler->history_status[j][0]);
			free(sampler->history_status[j]);
		}
		free(sampler->history_status);
	}
	for(int j = 0; j<sampler->chain_N; j++){
		free(sampler->model_status[j]);
	}
	free(sampler->model_status);

	free(sampler->ref_chain_status);
	free(sampler->ref_chain_ids);
	

	//Trouble shooting
	if(sampler->log_ll || sampler->log_lp){
		deallocate_3D_array(sampler->ll_lp_output,
			sampler->chain_N, sampler->N_steps, 2);
	}
	free(sampler->prop_MH_factor);
	
}

//TODO 
//shouldn't history store model numbers too?
void update_history(sampler *sampler, double *new_params, int *new_param_status, int chain_index)
{
	if(sampler->current_hist_pos[chain_index] < sampler->history_length-1)
	{
		sampler->current_hist_pos[chain_index]+=1;
	}
	else
	{
		//debugger_print(__FILE__,__LINE__,"RESTART HISTORY "+std::to_string(chain_index)+" "+std::to_string(sampler->current_hist_pos[chain_index]));
		sampler->current_hist_pos[chain_index] = 0;
	}
	for (int i =0; i < sampler->max_dim; i++)
	{
		sampler->history[chain_index][sampler->current_hist_pos[chain_index]][i] =
			new_params[i];
		sampler->history_status[chain_index][sampler->current_hist_pos[chain_index]][i] =
			new_param_status[i];
	}
		

}


void write_stat_file(sampler *sampler, 
		std::string filename
		)
{
	int rejected_swps=0, accepted_swps=0;
	int *accepted_steps=(int *)malloc(sizeof(int)*sampler->chain_N);
	int *rejected_steps=(int *)malloc(sizeof(int)*sampler->chain_N);
		
	for (int i =0;i<sampler->chain_N; i++)
	{
		accepted_swps+=sampler->swap_accept_ct[i];
		rejected_swps+=sampler->swap_reject_ct[i];
		accepted_steps[i]=sampler->step_accept_ct[i];
		rejected_steps[i]=sampler->step_reject_ct[i];
	}
	double total_swps= accepted_swps + rejected_swps;
	double accepted_swp_fraction = (double)accepted_swps/(total_swps);		
	double rejected_swp_fraction = (double)rejected_swps/(total_swps);		
	
	std::ofstream out_file;
	out_file.open(filename);	
	out_file.precision(4);
	//File variables
	int width= 80;
	int third = (int)((double)width/3.);
	int half = (int)((double)width/2.);
	int fourth = (int)((double)width/4.);
	int fifth = (int)((double)width/5.);
	int sixth = (int)((double)width/6.);
	int seventh = (int)((double)width/7.);
	int eighth = (int)((double)width/8.);
	int sixteenth = (int)((double)width/16.);
	
	//Sampler parameters
	out_file<<std::setw(width)<<std::left<<
		"Parameters of sampler: "<<std::endl;	
	out_file<<
		std::setw(eighth)<<std::left<<
		"Max/Min dimension: "<<
		std::setw(sixteenth)<<std::left<<
		sampler->max_dim<<
		std::setw(sixteenth)<<std::left<<
		sampler->min_dim<<
		std::setw(fourth)<<std::left<<
		"Length of History: "<<
		std::setw(fourth)<<std::left<<
		sampler->history_length<<
		std::endl;
	out_file<<
		std::setw(fourth)<<std::left<<
		"Number of Chains: "<<
		std::setw(fourth)<<std::left<<
		sampler->chain_N<<
		std::setw(fourth)<<std::left<<
		"Threads/Stochastic: "<<
		std::setw(fourth)<<std::left<<
		sampler->num_threads<<" / "<<sampler->pool<<
		std::endl;
	out_file<<
		std::setw(fourth)<<std::left<<
		"Steps per chain: "<<
		std::setw(fourth)<<std::left<<
		sampler->N_steps<<
		std::setw(fourth)<<std::left<<
		"Chain steps/chain swap: "<<
		std::setw(fourth)<<std::left<<
		sampler->swp_freq<<
		std::endl;
	out_file<<
		std::setw(fourth)<<std::left<<
		"CPU time - sampler (min): "<<
		std::setw(fourth)<<std::left<<
		sampler->time_elapsed_cpu/60.<<
		std::setw(fourth)<<std::left<<
		"CPU time - sampler(sec): "<<
		std::setw(fourth)<<std::left<<
		sampler->time_elapsed_cpu<<
		std::endl;
	out_file<<
		std::setw(fourth)<<std::left<<
		"Wall time - sampler (min): "<<
		std::setw(fourth)<<std::left<<
		sampler->time_elapsed_wall/60.<<
		std::setw(fourth)<<std::left<<
		"Wall time - sampler(sec): "<<
		std::setw(fourth)<<std::left<<
		sampler->time_elapsed_wall<<
		std::endl;
	out_file<<
		std::setw(fourth)<<std::left<<
		"CPU time - auto-corr (min): "<<
		std::setw(fourth)<<std::left<<
		sampler->time_elapsed_cpu_ac/60.<<
		std::setw(fourth)<<std::left<<
		"CPU time - auto-corr(sec): "<<
		std::setw(fourth)<<std::left<<
		sampler->time_elapsed_cpu_ac<<
		std::endl;
	out_file<<
		std::setw(fourth)<<std::left<<
		"Wall time - auto-corr (min): "<<
		std::setw(fourth)<<std::left<<
		sampler->time_elapsed_wall_ac/60.<<
		std::setw(fourth)<<std::left<<
		"Wall time - auto-corr(sec): "<<
		std::setw(fourth)<<std::left<<
		sampler->time_elapsed_wall_ac<<
		std::endl;
	out_file<<
		std::setw(width)<<std::left<<
		"Probabilities of steps (Gaussian, DE, MMALA, FISHER, RJ, KDE): "<<std::endl;
	out_file<<
		std::setw(fifth)<<std::left<<
		sampler->step_prob[0][0]<<
		std::setw(fifth)<<std::left<<
		sampler->step_prob[0][1]<<
		std::setw(fifth)<<std::left<<
		sampler->step_prob[0][2]<<
		std::setw(fifth)<<std::left<<
		sampler->step_prob[0][3]<<
		std::setw(fifth)<<std::left<<
		sampler->step_prob[0][4]<<
		std::setw(fifth)<<std::left<<
		sampler->step_prob[0][5]<<
		std::endl;
	out_file<<std::endl;

	out_file<<std::setw(width)<<"Chain temperature, total number of steps per chain, fraction of accepted swaps"<<std::endl;
	out_file<<
		std::setw(fourth)<<"ID" <<
		std::setw(fourth)<<"Temp" <<
		std::setw(fourth)<<"Step Number" <<
		std::setw(fourth)<<"Accepted Swaps" <<
		std::endl;
	int ts;
	int swpt;
	int swpa;
	for(int i =0; i<sampler->chain_N; i++){
		ts = sampler->fish_accept_ct[i]+sampler->fish_reject_ct[i]+
			sampler->de_accept_ct[i]+sampler->de_reject_ct[i]+
			sampler->mmala_accept_ct[i]+sampler->mmala_reject_ct[i]+
			sampler->RJstep_accept_ct[i]+sampler->RJstep_reject_ct[i]+
			sampler->gauss_accept_ct[i]+sampler->gauss_reject_ct[i]+
			sampler->kde_accept_ct[i] + sampler->kde_reject_ct[i];
		swpa = sampler->swap_accept_ct[i];
		swpt = sampler->swap_reject_ct[i] + swpa;
		out_file<<std::setw(fourth)<<i<<
			std::setw(fourth)<<sampler->chain_temps[i]<<
			std::setw(fourth)<<ts<<
			std::setw(fourth)<<(double)swpa/swpt<<
			std::endl;
	}
	out_file<<std::endl;
	double total_step_type;
	out_file<<
		std::setw(width)<<std::left<<
		"Fraction of steps per type : "<<std::endl;
	out_file<<
		std::setw(sixth)<<std::left<<
		"Chain Number"<<
		std::setw(sixth)<<std::left<<
		"Gaussian"<<
		std::setw(sixth)<<std::left<<
		"Diff. Ev."<<
		std::setw(sixth)<<std::left<<
		"MMALA"<<
		std::setw(sixth)<<std::left<<
		"Fisher"<<
		std::setw(sixth)<<std::left<<
		"RJ"<<
		std::setw(sixth)<<std::left<<
		"KDE"<<
		std::endl;
	for (int i =0; i < sampler->chain_N; i++){	
	 	total_step_type= sampler->num_gauss[i]+sampler->num_mmala[i]+
				sampler->num_de[i]+sampler->num_fish[i]+sampler->num_RJstep[i]+sampler->num_kde[i];
		out_file<<
			std::setw(sixth)<<std::left<<
			i<<
			std::setw(sixth)<<std::left<<
			(double)sampler->num_gauss[i]/total_step_type<<
			std::setw(sixth)<<std::left<<
			(double)sampler->num_de[i]/total_step_type<<
			std::setw(sixth)<<std::left<<
			(double)sampler->num_mmala[i]/total_step_type<<
			std::setw(sixth)<<std::left<<
			(double)sampler->num_fish[i]/total_step_type<<
			std::setw(sixth)<<std::left<<
			(double)sampler->num_RJstep[i]/total_step_type<<
			std::setw(sixth)<<std::left<<
			(double)sampler->num_kde[i]/total_step_type<<
			std::endl;
		
	}
	out_file<<std::endl;	
	//########################################################
	
	//out_file<<
	//	std::setw(width)<<std::left<<
	//	"Final width of Gaussian random number per step type: "<<std::endl;
	//out_file<<
	//	std::setw(sixth)<<std::left<<
	//	"Chain Number"<<
	//	std::setw(sixth)<<std::left<<
	//	"Gaussian"<<
	//	std::setw(sixth)<<std::left<<
	//	"Diff. Ev."<<
	//	std::setw(sixth)<<std::left<<
	//	"MMALA"<<
	//	std::setw(sixth)<<std::left<<
	//	"Fisher"<<
	//	std::setw(sixth)<<std::left<<
	//	"RJ"<<
	//	std::endl;
	//for (int i =0; i < sampler->chain_N; i++){	
	//	out_file<<
	//		std::setw(sixth)<<std::left<<
	//		i<<
	//		std::setw(sixth)<<std::left<<
	//		(double)sampler->randgauss_width[i][0]<<
	//		std::setw(sixth)<<std::left<<
	//		(double)sampler->randgauss_width[i][1]<<
	//		std::setw(sixth)<<std::left<<
	//		(double)sampler->randgauss_width[i][2]<<
	//		std::setw(sixth)<<std::left<<
	//		(double)sampler->randgauss_width[i][3]<<
	//		std::setw(sixth)<<std::left<<
	//		(double)sampler->randgauss_width[i][4]<<
	//		std::endl;
	//	
	//}
	//out_file<<std::endl;	

	//#######################################################

	double acc_total=0;
	double rej_total=0;

	out_file<<std::left<<"Fraction of accepted steps for each chain temp (by step): "<<std::endl;
	out_file.width(width);
	
	out_file<<
		std::left<<std::setw(seventh)<<"Chain ID"<<
		std::setw(seventh)<<std::left<<"GAUSS"<<
		std::setw(seventh)<<std::left<<"Diff Ev"<<
		std::setw(seventh)<<std::left<<"MMALA"<<
		std::setw(seventh)<<std::left<<"Fisher"<<
		std::setw(seventh)<<std::left<<"RJ"<<
		std::setw(seventh)<<std::left<<"KDE"<<
		std::setw(seventh)<<std::left<<"Total"<<
		std::endl;

	double total ;
	double acc_frac = 0;
	double rej_frac = 0;
	double gtotal ;
	double detotal ;
	double mmtotal ;
	double ftotal ;
	double kdetotal ;
	double RJtotal ;
	double gtotal_total =0;
	double detotal_total =0;
	double mmtotal_total =0;
	double ftotal_total =0;
	double kdetotal_total =0;
	double RJtotal_total =0;
	double total_total =0;
	double gacc_frac = 0;
	double deacc_frac = 0;
	double mmacc_frac = 0;
	double facc_frac = 0;
	double kdeacc_frac = 0;
	double RJacc_frac = 0;
	double gacc_frac_total = 0;
	double deacc_frac_total = 0;
	double mmacc_frac_total = 0;
	double facc_frac_total = 0;
	double kdeacc_frac_total = 0;
	double RJacc_frac_total = 0;
	double acc_frac_total = 0;
	for (int i =0; i<sampler->chain_N;i++){
		//####################################
		gtotal = sampler->gauss_accept_ct[i]+sampler->gauss_reject_ct[i];
		gacc_frac = (double)sampler->gauss_accept_ct[i]/gtotal;

		detotal = sampler->de_accept_ct[i]+sampler->de_reject_ct[i];
		deacc_frac = (double)sampler->de_accept_ct[i]/detotal;

		mmtotal = sampler->mmala_accept_ct[i]+sampler->mmala_reject_ct[i];
		mmacc_frac = (double)sampler->mmala_accept_ct[i]/mmtotal;
	
		ftotal = sampler->fish_accept_ct[i]+sampler->fish_reject_ct[i];
		facc_frac = (double)sampler->fish_accept_ct[i]/ftotal;


		kdetotal = sampler->kde_accept_ct[i]+sampler->kde_reject_ct[i];
		kdeacc_frac = (double)sampler->kde_accept_ct[i]/kdetotal;

		RJtotal = sampler->RJstep_accept_ct[i]+sampler->RJstep_reject_ct[i];
		RJacc_frac = (double)sampler->RJstep_accept_ct[i]/RJtotal;

		total = accepted_steps[i]+rejected_steps[i];
		acc_frac = (double)accepted_steps[i]/total;
		
		//####################################
		gtotal_total += gtotal;
		gacc_frac_total+=(double)sampler->gauss_accept_ct[i];

		mmtotal_total += mmtotal;
		mmacc_frac_total+=(double)sampler->mmala_accept_ct[i];

		detotal_total += detotal;
		deacc_frac_total+=(double)sampler->de_accept_ct[i];

		ftotal_total += ftotal;
		facc_frac_total+=(double)sampler->fish_accept_ct[i];

		kdetotal_total += kdetotal;
		kdeacc_frac_total+=(double)sampler->kde_accept_ct[i];

		RJtotal_total += RJtotal;
		RJacc_frac_total+=(double)sampler->RJstep_accept_ct[i];

		total_total += total;
		acc_frac_total+=(double)accepted_steps[i];
		//####################################
		out_file<<std::left<<std::setw(seventh)<<i<<
			std::left<<std::setw(seventh)<<gacc_frac<<
			std::left<<std::setw(seventh)<<deacc_frac<<
			std::left<<std::setw(seventh)<<mmacc_frac<<
			std::left<<std::setw(seventh)<<facc_frac<<
			std::left<<std::setw(seventh)<<RJacc_frac<<
			std::left<<std::setw(seventh)<<kdeacc_frac<<
			std::left<<std::setw(seventh)<<acc_frac<<
			std::endl;
	}
	out_file<<
		std::left<<std::setw(seventh)<<"TOTAL: "<<
		std::setw(seventh)<<(double)gacc_frac_total/(gtotal_total)<<
		std::setw(seventh)<<(double)deacc_frac_total/(detotal_total)<<
		std::setw(seventh)<<(double)mmacc_frac_total/(mmtotal_total)<<
		std::setw(seventh)<<(double)facc_frac_total/(ftotal_total)<<
		std::setw(seventh)<<(double)RJacc_frac_total/(RJtotal_total)<<
		std::setw(seventh)<<(double)kdeacc_frac_total/(kdetotal_total)<<
		std::setw(seventh)<<(double)acc_frac_total/(total_total)<<
		std::endl;

	out_file<<std::endl;	
	
	//Accepted rejected swaps
	out_file<<std::left<<"Number of accepted and rejected swaps for all chain temps: "
			<<std::endl;
	out_file<<std::left<<std::setw(fourth)<<"Accepted: "<<std::setw(fourth)<<accepted_swps
			<<std::setw(fourth)<<"Rejected: "<<std::setw(fourth)<< 
			rejected_swps<<std::endl;
	out_file<<std::left<<std::setw(fourth)<<"Accepted fraction: "<<std::setw(fourth)<<
			accepted_swp_fraction
			<<std::setw(fourth)<<"Rejected fraction: "<<std::setw(fourth)<< 
			rejected_swp_fraction<<std::endl;
	int nan_ct = 0;
	for(int i = 0 ; i<sampler->chain_N; i++){
		nan_ct+=sampler->nan_counter[i];
	}
	out_file<<std::left<<std::setw(half)<<"Nans in Fisher: "<<std::setw(half)<<nan_ct<<std::endl;
	out_file.close();
	free(accepted_steps);free(rejected_steps);
}
/*! \brief Routine that writes metadata and final positions of a sampler to a checkpoint file
 *
 */
void write_checkpoint_file(sampler *sampler, std::string filename)
{
	std::ofstream checkfile;
	checkfile.open(filename);
	checkfile.precision(15);
	//checkfile<<sampler->min_dim<<" , "<<sampler->max_dim<<" , "<<sampler->chain_N<<std::endl;//Dim, chain_N
	checkfile<<sampler->min_dim<<" , "<<sampler->max_dim<<" , "<<sampler->chain_N<<" , "<<sampler->nested_model_number<<std::endl;//min dim, max dim, chain_N, nested model number
	//Chain temps
	checkfile<<sampler->chain_temps[0];
	for(int i =1 ; i<sampler->chain_N; i++){
		checkfile<<" , "<<sampler->chain_temps[i];
	}
	checkfile<<std::endl;
	//step widths
	//for(int i =0 ; i<sampler->chain_N; i++){
	//	checkfile<<sampler->randgauss_width[i][0];
	//	for(int j =1 ;j<sampler->types_of_steps; j++){
	//		checkfile<<" , "<<sampler->randgauss_width[i][j];
	//	}
	//	checkfile<<std::endl;
	//}
	//step widths
	for(int i =0 ; i<sampler->chain_N; i++){
		checkfile<<sampler->randgauss_width_number[i][0];
		for(int j =1 ;j<sampler->types_of_steps; j++){
			checkfile<<" , "<<sampler->randgauss_width_number[i][j];
		}
		checkfile<<std::endl;
		for(int k = 0 ; k<sampler->types_of_steps; k++){
			checkfile<<sampler->randgauss_width[i][k][0];
			for(int j =1 ;j<sampler->randgauss_width_number[i][k]; j++){
				checkfile<<" , "<<sampler->randgauss_width[i][k][j];
			}
			checkfile<<std::endl;
		}
	}
	//final position
	for(int i =0 ; i< sampler->chain_N; i++){
		int pos = sampler->chain_pos[i];
		checkfile<<sampler->output[i][pos][0];
		for(int j =1; j<sampler->max_dim;j++){
			checkfile<<" , "<<sampler->output[i][pos][j];
		}
		for(int j =0; j<sampler->max_dim;j++){
			checkfile<<" , "<<sampler->param_status[i][pos][j];
			//std::cout<<sampler->param_status[i][pos][j]<<std::endl;
		}
		checkfile<<" , "<<sampler->model_status[i][pos];
		checkfile<<std::endl;
	}
	//history -- only adds if all chains are primed
	bool de_primed= true;
	for(int i =0 ; i<sampler->chain_N; i++){
		if(!sampler->de_primed[i])
			de_primed=false	;
	}
	if(de_primed){
		checkfile<<"true"<<std::endl;
		
		for(int i =0 ; i<sampler->chain_N; i++){
			checkfile<<sampler->history[i][0][0];
			for (int j = 1 ; j<sampler->max_dim; j++){
				checkfile<<" , "<<sampler->history[i][0][j];
			}
			for(int k=1; k<sampler->history_length; k++){
				for (int j = 0 ; j<sampler->max_dim; j++){

					checkfile<<" , "<<sampler->history[i][k][j];
				}
			}
			if(sampler->RJMCMC){
				//checkfile<<sampler->history_status[i][0][0];
				//for (int j = 1 ; j<sampler->max_dim; j++){
				//	checkfile<<" , "<<sampler->history_status[i][0][j];
				//}
				//for(int k=1; k<sampler->history_length; k++){
				//	for (int j = 0 ; j<sampler->max_dim; j++){

				//		checkfile<<" , "<<sampler->history_status[i][k][j];
				//	}
				//}
				for(int k=0; k<sampler->history_length; k++){
					for (int j = 0 ; j<sampler->max_dim; j++){

						checkfile<<" , "<<sampler->history_status[i][k][j];
					}
				}
			}
			
			checkfile<<std::endl;
		}
 		checkfile<<sampler->current_hist_pos[0];
                for(int i =1 ; i<sampler->chain_N; i++){
                        checkfile<<" , "<<sampler->current_hist_pos[i];
                }
	}
	else{
		checkfile<<"false"<<std::endl;
	}
	checkfile.close();
}


/*! \brief load temperatures from checkpoint file 
 *
 * Assumed the temps array is already allocated in memory for the correct number of chains
 * 	
 * Just a utility routine to read temperatures from checkpoint file
 *
 * It would be easy to read in the chain number and allocate memory in the function, but I prefer to leave allocation/deallocation up to the client
 */
void load_temps_checkpoint_file(std::string check_file, double *temps, int chain_N)
{
	std::fstream file_in;
	file_in.open(check_file, std::ios::in);
	std::string line;
	std::string item;
	int i;
	if(file_in){
		std::getline(file_in,line);

		//Second Row -- temps
		std::getline(file_in,line);
		std::stringstream lineStreamtemps(line);
		i=0;
		while(std::getline(lineStreamtemps, item, ',')){
			temps[i] = std::stod(item);
			i++;
		}
	}
	file_in.close();
}
/*! \brief Load dimension from checkpoint file
 *
 * If RJ -- returns min_dimension and max_dimension
 *
 * If not RJ -- returns dimension in min_dimension and max_dimension
 */
int dimension_from_checkpoint_file(std::string check_file, int *min_dimension, int *max_dimension)
{
	std::fstream file_in;
	file_in.open(check_file, std::ios::in);
	std::string line;
	std::string item;
	int i;
	if(file_in){
		//First row -- dim, chain_N
		std::getline(file_in,line);
		std::stringstream lineStream(line);
		std::getline(lineStream, item, ',');
		*min_dimension= std::stod(item);
		std::getline(lineStream, item, ',');
		*max_dimension = std::stod(item);
	}
	else{
		std::cout<<"ERROR -- checkpoint file not found"<<std::endl;	
		return 1;
	}
	file_in.close();
	return 0;
	
	
}
/*! \brief Load chain number from checkpoint file
 *
 */
int chain_number_from_checkpoint_file(std::string check_file, int *chain_N)
{
	std::fstream file_in;
	file_in.open(check_file, std::ios::in);
	std::string line;
	std::string item;
	int i;
	if(file_in){
		//First row -- dim, chain_N
		std::getline(file_in,line);
		std::stringstream lineStream(line);
		std::getline(lineStream, item, ',');//Min dim
		std::getline(lineStream, item, ',');//Max dim
		std::getline(lineStream, item, ',');//Chain_number
		*chain_N = std::stod(item);
	}
	else{
		std::cout<<"ERROR -- checkpoint file not found"<<std::endl;	
		return 1;
	}
	file_in.close();
	return 0;
	
	
}

/*! \brief load checkpoint file into sampler struct
 *
 * *NOTE* -- allocate_sampler called in function -- MUST deallocate manually
 *
 * *NOTE* -- sampler->chain_temps allocated internally -- MUST free manually
 */
void load_checkpoint_file(std::string check_file,sampler *samplerptr)
{
	std::fstream file_in;
	file_in.open(check_file, std::ios::in);
	std::string line;
	std::string item;
	int i;
	if(file_in){
		//First row -- dim, chain_N
		std::getline(file_in,line);
		std::stringstream lineStream(line);
		std::getline(lineStream, item, ',');
		samplerptr->min_dim = std::stod(item);
		std::getline(lineStream, item, ',');
		samplerptr->max_dim = std::stod(item);
		if(samplerptr->min_dim == samplerptr->max_dim){
			samplerptr->dimension=samplerptr->min_dim;
			samplerptr->RJMCMC = false;
		}
		else{
			samplerptr->RJMCMC=true;
		}
		std::getline(lineStream, item, ',');
		samplerptr->chain_N = std::stod(item);
		std::getline(lineStream, item, ',');
		samplerptr->nested_model_number = std::stod(item);

		//Second Row -- temps
		std::getline(file_in,line);
		std::stringstream lineStreamtemps(line);
		samplerptr->chain_temps = (double *)malloc(sizeof(double)*samplerptr->chain_N);
		i=0;
		while(std::getline(lineStreamtemps, item, ',')){
			samplerptr->chain_temps[i] = std::stod(item);
			i++;
		}
		
		//###################################
		//Allocate memory, now we have initial parameters
		allocate_sampler_mem(samplerptr);
		//###################################


		//third row+chain_N -- step widths
		//for(int j =0 ;j<sampler->chain_N; j++){
		//	i=0;
		//	std::getline(file_in,line);
		//	std::stringstream lineStreamwidths(line);
		//	while(std::getline(lineStreamwidths, item, ',')){
		//		sampler->randgauss_width[j][i] = std::stod(item);
		//		i++;
		//	}
		//}
		//third row+chain_N -- step widths
		for(int j =0 ;j<samplerptr->chain_N; j++){
			i=0;
			std::getline(file_in,line);
			std::stringstream lineStreamwidths(line);
			while(std::getline(lineStreamwidths, item, ',')){
				samplerptr->randgauss_width_number[j][i] = std::stoi(item);
				i++;
			}
			for(int k = 0 ; k<samplerptr->types_of_steps;k++){
				i=0;
				std::getline(file_in,line);
				std::stringstream lineStreamwidths(line);
				while(std::getline(lineStreamwidths, item, ',')){
					samplerptr->randgauss_width[j][k][i] = std::stod(item);
					i++;
				}
			}
		}

		//row 3 +chain_N  to 3+2 chain_N-- initial positions
		for(int j =0 ;j<samplerptr->chain_N; j++){
			i=0;
			std::getline(file_in,line);
			std::stringstream lineStreampos(line);
			//while(i<sampler->max_dim && std::getline(lineStreampos, item, ',')){
			while(i<samplerptr->max_dim  ){
				std::getline(lineStreampos, item, ',');
				samplerptr->output[j][0][i] = std::stod(item);
				i++;
			}
			i=0;
			//while(std::getline(lineStreampos, item, ',')){
			while(i<samplerptr->max_dim  ){
				std::getline(lineStreampos, item, ',');
				samplerptr->param_status[j][0][i] = std::stod(item);
				i++;
			}
			//if(samplerptr->nested_model_number>0){
			//	i=0;
			//	while(i<samplerptr->nested_model_number  ){
			//		std::getline(lineStreampos, item, ',');
			//		samplerptr->model_status[j][0][i] = std::stod(item);
			//		i++;
			//	}
			//}
			std::getline(lineStreampos, item, ',');
			samplerptr->model_status[j][0] = std::stoi(item);
		}
		std::getline(file_in,line);
		std::stringstream lineStreamprimed(line);
		std::getline(lineStreamprimed, item, ',');
		if(item =="true"){
			for(int j =0 ;j<samplerptr->chain_N; j++){
				std::getline(file_in,line);
				std::stringstream lineStreamhist(line);
				i=0;
				while(i < samplerptr->max_dim*samplerptr->history_length){	
					std::getline(lineStreamhist,item,','); 
					int step = i/samplerptr->max_dim ;
					int pos = i%samplerptr->max_dim;
					samplerptr->history[j][step][pos] = std::stod(item);	
					i++;
				}
				if(samplerptr->RJMCMC){
					i=0;
					while(std::getline(lineStreamhist,item,',')){	
						int step = i/samplerptr->max_dim ;
						int pos = i%samplerptr->max_dim;
						samplerptr->history_status[j][step][pos] = std::stoi(item);	
						i++;
					}
				}
			}
			std::getline(file_in,line);
			std::stringstream lineStreamhist(line);
			i=0;
			while(std::getline(lineStreamhist,item,',')){   
			        samplerptr->current_hist_pos[i] = std::stoi(item);
			
			        i++;
			}


			for(int j =0 ;j<samplerptr->chain_N; j++){
				samplerptr->de_primed[j] =true;
				if(samplerptr->kde_step){
					update_kde_cov(samplerptr, j);
				}
			}
		}
		else{
			for(int j =0 ;j<samplerptr->chain_N; j++)
				samplerptr->de_primed[j] =false;

		}
	
		//exit(1);
		
	}
	else{std::cout<<"ERROR -- File "<<check_file<<" not found"<<std::endl; exit(1);}
	file_in.close();
	//if(samplerptr->fisher_exist){
	//	//check whether or not we need to update the fisher
	//	for(int i=0 ; i<samplerptr->chain_N; i++){
	//		update_fisher(samplerptr, samplerptr->output[i][0], samplerptr->param_status[i][0],samplerptr->model_status[i][0],i);	
	//		if(samplerptr->proper_fisher){
	//			iterate_fisher(samplerptr, i);
	//		}
	//	}
	//}
	if(samplerptr->fisher_exist){
		//check whether or not we need to update the fisher
		omp_set_num_threads(samplerptr->num_threads);
		#pragma omp parallel ADOLC_OPENMP
		//#pragma omp parallel 
		{
		#pragma omp for 
		for(int i=0 ; i<samplerptr->chain_N; i++){
			update_fisher(samplerptr, samplerptr->output[i][0], samplerptr->param_status[i][0],&(samplerptr->model_status[i][0]),i);	
			if(samplerptr->fisher_update_ct[i] +1 == samplerptr->fisher_update_number){
				//If NAN in fisher, it won't set the values. 
				//Just set to junk so memory error doesn't occur
				//It should be immediately overwritten after first step
				for(int j = 0 ; j<samplerptr->max_dim; j++){
					for(int l = 0 ; l<samplerptr->max_dim ; l++){
						samplerptr->fisher_vecs[i][j][l] = 0;
						samplerptr->fisher_vecs_prop[i][j][l] = 0;
						samplerptr->fisher_matrix[i][j][l] = 0;
						samplerptr->fisher_matrix_prop[i][j][l] = 0;
					}
					samplerptr->fisher_vecs[i][j][j] = 1;
					samplerptr->fisher_vecs_prop[i][j][j] = 1;
					samplerptr->fisher_matrix[i][j][j] = 1;
					samplerptr->fisher_matrix_prop[i][j][j] = 1;
				
					samplerptr->fisher_vals[i][j] = 1;
					samplerptr->fisher_vals_prop[i][j] = 1;
				}
			}
			if(samplerptr->proper_fisher){
				iterate_fisher(samplerptr, i);
			}
		}
		}
	}

}

void assign_ct_p(sampler *samplerptr, int step, int chain_index, int gauss_dim)
{
	if(step ==0) {
		samplerptr->gauss_accept_ct[chain_index]+=1;
		samplerptr->gauss_accept_ct_per_dim[chain_index][gauss_dim]+=1;
	}
	else if(step ==1) samplerptr->de_accept_ct[chain_index]+=1;
	else if(step ==2) samplerptr->mmala_accept_ct[chain_index]+=1;
	else if(step ==3) samplerptr->fish_accept_ct[chain_index]+=1;
	else if(step ==4) samplerptr->RJstep_accept_ct[chain_index]+=1;
	else if(step ==5) samplerptr->kde_accept_ct[chain_index]+=1;
}
void assign_ct_m(sampler *samplerptr, int step, int chain_index, int gauss_dim)
{
	if(step ==0) {
		samplerptr->gauss_reject_ct[chain_index]+=1;
		samplerptr->gauss_reject_ct_per_dim[chain_index][gauss_dim]+=1;
	}
	else if(step ==1) samplerptr->de_reject_ct[chain_index]+=1;
	else if(step ==2) samplerptr->mmala_reject_ct[chain_index]+=1;
	else if(step ==3) samplerptr->fish_reject_ct[chain_index]+=1;
	else if(step ==4) samplerptr->RJstep_reject_ct[chain_index]+=1;
	else if(step ==5) samplerptr->kde_reject_ct[chain_index]+=1;
}

void assign_initial_pos(sampler *samplerptr,double *initial_pos, int *initial_status, int initial_model_status,double **ensemble_initial_pos,int **ensemble_initial_status,int *ensemble_initial_model_status,double *seeding_var) 
{
	if(ensemble_initial_pos){
		debugger_print(__FILE__,__LINE__,"Using initial ENSEMBLE position");
		if(samplerptr->nested_model_number >0){
			for(int j = 0 ; j<samplerptr->chain_N; j++){
				samplerptr->model_status[j][0] = ensemble_initial_model_status[j];
			}
		}
		bool valid_pos[samplerptr->chain_N];
		omp_set_num_threads(samplerptr->num_threads);
		#pragma omp parallel ADOLC_OPENMP
		//#pragma omp parallel 
		{
		#pragma omp for 
		for (int j=0;j<samplerptr->chain_N;j++){
			samplerptr->de_primed[j]=false;
			for (int i = 0; i<samplerptr->max_dim; i++)
			{
				samplerptr->output[j][0][i] = ensemble_initial_pos[j][i];
				samplerptr->param_status[j][0][i] = ensemble_initial_status[j][i];
				//Only doing this last loop because there is sometimes ~5 elements 
				//not initialized on the end of the output, which screw up plotting
				for(int l =1; l<samplerptr->N_steps; l++){
					samplerptr->output[j][l][i] = 0;
					samplerptr->param_status[j][l][i] = 1;
					//samplerptr->output[j][l][i] = ensemble_initial_pos[j][i];
					//samplerptr->param_status[j][l][i] = ensemble_initial_status[j][i];
				}
				
			}
			double lp = samplerptr->lp(samplerptr->output[j][0],samplerptr->param_status[j][0],samplerptr->model_status[j][0],samplerptr->interfaces[j], samplerptr->user_parameters[j]);
			if(lp == limit_inf){
				
				valid_pos[j] = false;	
				debugger_print(__FILE__,__LINE__,"Chain "+std::to_string(j)+" has invalid position");
			}

			else{
				valid_pos[j] = true;
				samplerptr->current_likelihoods[j] =
					samplerptr->ll(samplerptr->output[j][0],samplerptr->param_status[j][0],samplerptr->model_status[j][0],samplerptr->interfaces[j], samplerptr->user_parameters[j]);
			}
		}
		bool all_bad= true;
		double good_pos[samplerptr->max_dim];
		int good_status[samplerptr->max_dim];
		for(int i = 0 ; i<samplerptr->chain_N; i ++){
			if(valid_pos[i]) { 
				all_bad = false;
				for(int j = 0 ; j<samplerptr->max_dim; j++){
					good_pos[j] = samplerptr->output[i][0][j];
					good_status[j] = samplerptr->param_status[i][0][j];
				}
			}	
		}
		if(all_bad){
			debugger_print(__FILE__,__LINE__,"ALL positions are bad! Fix your initial positions!");
			exit(1);
		}
		#pragma omp for 
		for(int i = 0 ; i<samplerptr->chain_N; i ++){
			if(!valid_pos[i]){
				for(int j = 0 ; j<samplerptr->max_dim; j++){
					samplerptr->output[i][0][j] = good_pos[j];
					samplerptr->param_status[i][0][j] = good_status[j];
				}
				samplerptr->current_likelihoods[i] =
					samplerptr->ll(samplerptr->output[i][0],samplerptr->param_status[i][0],samplerptr->model_status[i][0],samplerptr->interfaces[i], samplerptr->user_parameters[i]);
			}
		}	
		}
		
	}
	else{
		if(samplerptr->nested_model_number >0){
			for(int j = 0 ; j<samplerptr->chain_N; j++){
				samplerptr->model_status[j][0] = initial_model_status;
			}
		}

		if(!seeding_var){ 
			for (int j=0;j<samplerptr->chain_N;j++){
				samplerptr->de_primed[j]=false;
				for (int i = 0; i<samplerptr->max_dim; i++)
				{
					//Only doing this last loop because there is sometimes ~5 elements 
					//not initialized on the end of the output, which screw up plotting
					for(int l =0; l<samplerptr->N_steps; l++){
						samplerptr->output[j][l][i] = initial_pos[i];
						samplerptr->param_status[j][l][i] = initial_status[i];
					}
					//std::cout<<initial_pos[i]<<" "<<initial_status[i]<<std::endl;
					
				}
				samplerptr->current_likelihoods[j] =
					 samplerptr->ll(samplerptr->output[j][0],samplerptr->param_status[j][0],samplerptr->model_status[j][0],samplerptr->interfaces[j], samplerptr->user_parameters[j]);
			}
		}
		//Seed non-zero chains normally with variance as specified
		else{
			int attempts = 0;
			int max_attempts = 10;
			double temp_pos[samplerptr->max_dim];
			int temp_status[samplerptr->max_dim];
			for (int j=0;j<samplerptr->chain_N;j++){
				samplerptr->de_primed[j]=false;
				if(j == 0){
					for (int i = 0; i<samplerptr->max_dim; i++)
					{
					//Only doing this last loop because there is sometimes ~5 elements 
					//not initialized on the end of the output, which screw up plotting
						for(int l =0; l<samplerptr->N_steps; l++){
							samplerptr->output[j][l][i] = initial_pos[i];
							samplerptr->param_status[j][l][i] = initial_status[i];
						}
					}
				}
				else{
					do{
						for(int i =0; i<samplerptr->max_dim; i++){
							if(initial_status[i] ==1){
								temp_pos[i] = gsl_ran_gaussian(samplerptr->rvec[j],seeding_var[i]) + initial_pos[i];
								temp_status[i]=1;
							}
							else{
								temp_pos[i]=0;
								temp_status[i]=0;
							}
						}
						attempts+=1;
					}while(samplerptr->lp(temp_pos, temp_status,initial_model_status,samplerptr->interfaces[j], samplerptr->user_parameters[j]) == limit_inf && attempts<max_attempts);
					attempts =0;
					if(samplerptr->lp(temp_pos, temp_status,initial_model_status,samplerptr->interfaces[j], samplerptr->user_parameters[j]) != limit_inf ){
						for(int i =0; i<samplerptr->max_dim;i++){
							for(int l =0; l<samplerptr->N_steps; l++){
								samplerptr->output[j][l][i] = temp_pos[i];
								samplerptr->param_status[j][l][i] = temp_status[i];
							}
						}
					}
					else{
						for(int i =0; i<samplerptr->max_dim;i++){
							for(int l =0; l<samplerptr->N_steps; l++){
								samplerptr->output[j][l][i] = initial_pos[i];
								samplerptr->param_status[j][l][i] = initial_status[i];
							}
						}
				
					}
				}
				
				samplerptr->current_likelihoods[j] =
					 samplerptr->ll(samplerptr->output[j][0],samplerptr->param_status[j][0],samplerptr->model_status[j][0],samplerptr->interfaces[ j], samplerptr->user_parameters[j]);
			}
		}
	}
	if(samplerptr->log_ll || samplerptr->log_lp){
		for(int i = 0 ; i<samplerptr->chain_N; i++){
			samplerptr->ll_lp_output[i][0][0] = samplerptr->current_likelihoods[i];
			samplerptr->ll_lp_output[i][0][1] = 
				samplerptr->lp(samplerptr->output[i][0],
				samplerptr->param_status[i][0],samplerptr->model_status[i][0],
				samplerptr->interfaces[i], samplerptr->user_parameters[i]);
		}
	}
	if(samplerptr->fisher_exist){
		//check whether or not we need to update the fisher
		omp_set_num_threads(samplerptr->num_threads);
		#pragma omp parallel ADOLC_OPENMP
		//#pragma omp parallel 
		{
		#pragma omp for 
		for(int i=0 ; i<samplerptr->chain_N; i++){
			update_fisher(samplerptr, samplerptr->output[i][0], samplerptr->param_status[i][0],&(samplerptr->model_status[i][0]),i);	
			if(samplerptr->fisher_update_ct[i] +1 == samplerptr->fisher_update_number){
				//If NAN in fisher, it won't set the values. 
				//Just set to junk so memory error doesn't occur
				//It should be immediately overwritten after first step
				for(int j = 0 ; j<samplerptr->max_dim; j++){
					for(int l = 0 ; l<samplerptr->max_dim ; l++){
						samplerptr->fisher_vecs[i][j][l] = 0;
						samplerptr->fisher_vecs_prop[i][j][l] = 0;
						samplerptr->fisher_matrix[i][j][l] = 0;
						samplerptr->fisher_matrix_prop[i][j][l] = 0;
					}
					samplerptr->fisher_vecs[i][j][j] = 1;
					samplerptr->fisher_vecs_prop[i][j][j] = 1;
					samplerptr->fisher_matrix[i][j][j] = 1;
					samplerptr->fisher_matrix_prop[i][j][j] = 1;
				
					samplerptr->fisher_vals[i][j] = 1;
					samplerptr->fisher_vals_prop[i][j] = 1;
				}
			}
			if(samplerptr->proper_fisher){
				iterate_fisher(samplerptr, i);
			}
		}
		}
	}
}

/*! \brief Timescale of the PT dynamics
 *
 * kappa in the the language of arXiv:1501.05823v3
 */
double PT_dynamical_timescale(int t0, /**<Timescale of the dyanmics*/
	int nu, /**<Initial amplitude (number of steps to base dynamics on)*/
	int t/**< current time*/
	)
{
	return (1./nu) * (double)(t0) / (t + t0);
}

/*! \brief updates the temperatures for a sampler such that all acceptance rates are equal
 * 
 * Follows the algorithm outlined in arXiv:1501.05823v3
 *
 * Fixed temperatures for the first and last chain
 *
 * used in MCMC_MH_dynamic_PT_alloc_internal
 *
 * For defined results, this should be used while the sampler is using non-pooling methods
 */
void update_temperatures(sampler *samplerptr,
	int t0,
	int nu,
	int t
	)
{
	//acceptance ratios -- as defined by arXiv:1501.05823v3
	//Need to transform, because I save the total acceptance count per chain
	//not the acceptance count between chains
	double *old_temps = new double[samplerptr->chain_N];

	for(int i =0 ; i<samplerptr->chain_N-1; i++){
		old_temps[i] = samplerptr->chain_temps[i];
	}
	double power;
	double kappa = PT_dynamical_timescale(t0, nu, t);
	//std::cout<<kappa<<std::endl;
	for (int i =1 ; i<samplerptr->chain_N-1; i++){
		power = kappa * (samplerptr->A[i] - samplerptr->A[i+1]);	
		samplerptr->chain_temps[i] = samplerptr->chain_temps[i-1] +
			(old_temps[i] - old_temps[i-1]) * std::exp(power);
		//Replace LL, since I store it already scaled by the old temperature
	} 
	delete [] old_temps;
}

void update_temp_neighborhoods(sampler *samplerptr){
	int cold_chain_ids[samplerptr->chain_N];
	int cold_chain_ct = 0;
	for (int i =0; i<samplerptr->chain_N; i++){
		if(fabs(samplerptr->chain_temps[i]-1)<DOUBLE_COMP_THRESH){
			cold_chain_ids[cold_chain_ct]=i;
			cold_chain_ct++;
		}
	}
	//The last ensemble won't have a cold chain on the end
	//This will only cause issues if every single chain is cold (won't happen)
	cold_chain_ids[cold_chain_ct]=samplerptr->chain_N;
	cold_chain_ct = 0;
	for(int i = 0 ; i<samplerptr->chain_N; i++){
		int neighbors_low = samplerptr->chain_radius;
		int neighbors_high = samplerptr->chain_radius;
		if(i >= cold_chain_ids[cold_chain_ct+1]){cold_chain_ct++;}
		if(i - cold_chain_ids[cold_chain_ct] < samplerptr->chain_radius){
			neighbors_low = i-cold_chain_ids[cold_chain_ct];
		}
		if(-i + cold_chain_ids[cold_chain_ct+1]-1 < samplerptr->chain_radius){
			neighbors_high = -i+cold_chain_ids[cold_chain_ct+1]-1;
		}
		samplerptr->chain_neighbors[i]=neighbors_low+neighbors_high;
		if(samplerptr->chain_neighborhoods[i]){
			free(samplerptr->chain_neighborhoods[i]);
		}
		samplerptr->chain_neighborhoods[i] = (double *)malloc(sizeof(double)*(neighbors_low+neighbors_high));
		if(samplerptr->chain_neighborhoods_ids[i]){
			free(samplerptr->chain_neighborhoods_ids[i]);
		}
		samplerptr->chain_neighborhoods_ids[i] = (int*) malloc(sizeof(int)*(neighbors_low+neighbors_high));
		for(int j = 0 ; j< neighbors_low; j++){
			samplerptr->chain_neighborhoods[i][j] = samplerptr->chain_temps[i-j-1];
			samplerptr->chain_neighborhoods_ids[i][j] = i-j-1;
		}
		for(int j = 0 ; j< neighbors_high; j++){
			samplerptr->chain_neighborhoods[i][j+neighbors_low] = samplerptr->chain_temps[i+j+1];
			samplerptr->chain_neighborhoods_ids[i][j+neighbors_low] = i+j+1;
		}
	}

}

void update_temperatures_full_ensemble(sampler *samplerptr,
	int t0,
	int nu,
	int t
	)
{
	//Either average dynamics during sampling or after sampling run
		
	if(!(samplerptr->linear_swapping)){
		double *old_temps = new double[samplerptr->chain_N];
		double max_temp = 0 ;
		int ensemble_chain_number = 0 ;
		bool search = true;
		for(int i =0 ; i<samplerptr->chain_N-1; i++){
			old_temps[i] = samplerptr->chain_temps[i];
			if(i != 0 && fabs(samplerptr->chain_temps[i] - 1) < DOUBLE_COMP_THRESH){
				max_temp = samplerptr->chain_temps[i-1];
				if(search){
					ensemble_chain_number = i;
					search=false;
				}
			}
		}
		if(max_temp <DOUBLE_COMP_THRESH ){
			max_temp = samplerptr->chain_temps[samplerptr->chain_N-1];
			ensemble_chain_number = samplerptr->chain_N;
		}
		double averaged_A[ensemble_chain_number];
		int chain_num_dt[ensemble_chain_number];
		for(int i = 0 ; i<ensemble_chain_number; i++){
			averaged_A[i]=0;
			chain_num_dt[i]=0;
		}
		for(int i = 0 ; i<samplerptr->chain_N; i++){
			averaged_A[i%ensemble_chain_number]+= samplerptr->A[i];
			chain_num_dt[i%ensemble_chain_number]+=1;
		}
		for(int i = 0 ; i<ensemble_chain_number; i++){
			if(chain_num_dt[i] != 0){
				averaged_A[i]/=chain_num_dt[i];
			}
			else{
				averaged_A[i]=0;
			}
		}
		double power;
		double kappa = PT_dynamical_timescale(t0, nu, t);
		for (int i =1 ; i<samplerptr->chain_N-1; i++){
			if( !(fabs(samplerptr->chain_temps[i] - 1) < DOUBLE_COMP_THRESH ||
				fabs(samplerptr->chain_temps[i] - max_temp) < DOUBLE_COMP_THRESH)){
				power = kappa * (averaged_A[i%ensemble_chain_number] - averaged_A[(i+1)%ensemble_chain_number]);
				samplerptr->chain_temps[i] = samplerptr->chain_temps[i-1] +
					(old_temps[i] - old_temps[i-1]) * std::exp(power);
				//Replace LL, since I store it already scaled by the old temperature
			}
		} 
		delete [] old_temps;
		update_temp_neighborhoods(samplerptr);
	}

	else{
		//######################################################
		//This worked, but we're trying something new
		//######################################################
		//acceptance ratios -- as defined by arXiv:1501.05823v3
		//Need to transform, because I save the total acceptance count per chain
		//not the acceptance count between chains
		double *old_temps = new double[samplerptr->chain_N];
		double max_temp = 0 ;
		int ensemble_chain_number = 0 ;
		bool search = true;
		for(int i =0 ; i<samplerptr->chain_N-1; i++){
			old_temps[i] = samplerptr->chain_temps[i];
			if(i != 0 && fabs(samplerptr->chain_temps[i] - 1) < DOUBLE_COMP_THRESH){
				max_temp = samplerptr->chain_temps[i-1];
				if(search){
					ensemble_chain_number = i;
					search=false;
				}
			}
		}
		if(max_temp <DOUBLE_COMP_THRESH ){
			max_temp = samplerptr->chain_temps[samplerptr->chain_N-1];
			ensemble_chain_number = samplerptr->chain_N;
		}
		double power;
		double kappa = PT_dynamical_timescale(t0, nu, t);
		//std::cout<<kappa<<std::endl;
		for (int i =1 ; i<samplerptr->chain_N-1; i++){
			if( !(fabs(samplerptr->chain_temps[i] - 1) < DOUBLE_COMP_THRESH ||
				fabs(samplerptr->chain_temps[i] - max_temp) < DOUBLE_COMP_THRESH)){
				if(i!=samplerptr->chain_N-1){
					power = kappa * (samplerptr->A[i] - samplerptr->A[i+1]);
				}
				else{
					power = kappa * (samplerptr->A[i] - samplerptr->A[i+1-ensemble_chain_number]);
					
				}
				samplerptr->chain_temps[i] = samplerptr->chain_temps[i-1] +
					(old_temps[i] - old_temps[i-1]) * std::exp(power);
				//Replace LL, since I store it already scaled by the old temperature
			}
		} 
		delete [] old_temps;
	}
}

/*!  \brief For the dynamic PT sampler, this is the function that starts the full sampler with the max number of chains
 *
 * The output file will be reused, but the positions are set back to zero (copying the current position to position zero)
 *
 * Assumes the output, chain_temps have been allocated in memory for the final number of chains chain_N and steps N_steps
 *
 * Allocates memory for the new sampler sampler_new -> it's the user's responsibility to deallocate with deallocate_sampler_mem
 *
 * checkpoint_file_start is the file used for continue_PTMCMC like functions, and can be used if the user wants to maintain the approximate positions and histories of the ensemble of chains after the dynamic temperature allocation. If the rest of the chains should just be copies of the base ensemble, pass an empty string (""). Since the chain distribution can change the number and distribution of chain temperatures (that's kinda the point..), the transfer is not perfect. Cold chains are prioritized, then the rest of the chains are filled in as close as possible.
 */
void initiate_full_sampler(sampler *sampler_new, sampler *sampler_old, /**<Dynamic sampler*/
	int chain_N_thermo_ensemble, /**< Number of chains used in the thermodynamic ensemble*/
	int chain_N,/**< Number of chains to use in the static sampler*/
	std::string chain_allocation_scheme, /**< Scheme to use to allocate any remaining chains*/
	std::string checkpoint_file_start /**< Base checkpoint file*/
	)
{
	//Check to see if more chains need to be allocated
	bool allocate_chains = true;
	bool use_checkpoint_file=false;
	if(sampler_old->chain_N == chain_N){	allocate_chains = false;}
	if(checkpoint_file_start !=""){ use_checkpoint_file = true;}
	
	//Allocate new sampler 
	sampler_new->chain_N = chain_N;
	sampler_new->ll = sampler_old->ll;
	sampler_new->lp = sampler_old->lp;
	sampler_new->fish = sampler_old->fish;
	sampler_new->N_steps = sampler_old->N_steps;
	sampler_new->dimension = sampler_old->dimension;
	sampler_new->min_dim = sampler_old->min_dim;
	sampler_new->max_dim = sampler_old->max_dim;
	sampler_new->num_threads = sampler_old->num_threads;
	sampler_new->numThreads = sampler_old->numThreads;
	//sampler_new->chain_temps = sampler_old->chain_temps;
	sampler_new->history_length = sampler_old->history_length;
	sampler_new->history_update = sampler_old->history_update;
	sampler_new->fisher_update_number = sampler_old->fisher_update_number;
	sampler_new->pool = sampler_old->pool;
	sampler_new->swp_freq = sampler_old->swp_freq;
	sampler_new->nested_model_number = sampler_old->nested_model_number;

	//This is a reference copy, not a value.
	sampler_new->output = sampler_old->output;
	sampler_new->param_status = sampler_old->param_status;
	sampler_new->model_status = sampler_old->model_status;
	//sampler_new->output = allocate_3D_array(chain_N, sampler_old->N_steps, sampler_old->dimension);

	sampler_new->fisher_exist = sampler_old->fisher_exist;
	sampler_new->progress = sampler_old->progress;

	//TODO
	//Need to populate the temps first, needed for allocation...
	for(int i = 0 ; i<sampler_old->chain_N; i++){
		sampler_new->chain_temps[i] = sampler_old->chain_temps[i];
	}
	for(int i = sampler_old->chain_N ; i<sampler_new->chain_N; i++){
		if(chain_allocation_scheme=="cold"){
			sampler_new->chain_temps[i] = sampler_old->chain_temps[0];
		}
		else if(chain_allocation_scheme=="double"){
			sampler_new->chain_temps[i] = sampler_old->chain_temps[i%sampler_old->chain_N];
		}
		else if(chain_allocation_scheme=="refine"){
			double prev_temp = sampler_old->chain_temps[i%sampler_old->chain_N];
			double next_temp = sampler_old->chain_temps[(i+1)%sampler_old->chain_N];
			sampler_new->chain_temps[i] = std::sqrt(prev_temp*next_temp);
		}
		else if(chain_allocation_scheme=="half_ensemble"){
			sampler_new->chain_temps[i] = 10;	
		}
	}
	allocate_sampler_mem(sampler_new);


	
	//Copy over pertinent sampler data: histories, step_widths, de_primed, current_hist_pos, copy current pos to first pos, current LL
	//Only for chains 0 - chain_N_thermo_ensemble
	for(int i =0; i<sampler_old->chain_N; i++){
		transfer_chain(sampler_new, sampler_old, i, i, false);
	}

	//Allocate extra chains if needed
	//Only for chains chain_N_thermo_ensemble - chain_N
	if(allocate_chains){
		if(chain_allocation_scheme =="cold"){	
			for(int i =sampler_old->chain_N; i<chain_N; i++){
				transfer_chain(sampler_new, sampler_old, i, 0,true );
			}		
		}
		if(chain_allocation_scheme =="double"){	
			for(int i =sampler_old->chain_N; i<chain_N; i++){
				transfer_chain(sampler_new, sampler_old, i, i%sampler_old->chain_N, true);
			}
		}
		//Continuously refine -- needs correction -- right now, its just adding the ``refinement'' temperatures repeatedly
		if(chain_allocation_scheme =="refine"){	
			for(int i =sampler_old->chain_N; i<chain_N; i++){
				transfer_chain(sampler_new, sampler_old, i, i%sampler_old->chain_N,true);
				double prev_temp = sampler_old->chain_temps[i%sampler_old->chain_N];
				double next_temp = sampler_old->chain_temps[(i+1)%sampler_old->chain_N];
				sampler_new->chain_temps[i] = std::sqrt(prev_temp*next_temp);
			}
		}
		if(chain_allocation_scheme =="half_ensemble"){	
			int i = sampler_old->chain_N;
			bool even_odd = true;
			int j = 2;
			while(i	< chain_N){
			//{
				transfer_chain(sampler_new, sampler_old, i, 0, true);
				i++;
				if(i==chain_N) break;
				while(j<sampler_old->chain_N){
					transfer_chain(sampler_new, sampler_old, i, j, true);
					j+=2;
					i++;
					if(i==chain_N) break;
				}
				if (even_odd){
					j = 1;
					even_odd = false;
				}
				else{
					j = 2;
					even_odd = true;
				}
			}
		}
		//Replaces the copied positions and histories with the old positions and histories
		if(use_checkpoint_file){
			copy_base_checkpoint_properties(checkpoint_file_start, sampler_new);
		}
	}


}
/*! \brief Copies positions and histories into chain ensemble, skipping the first set of temperatures
 *
 * *NOTE* -- allocate_sampler called in function -- MUST deallocate manually
 *
 * *NOTE* -- sampler->chain_temps allocated internally -- MUST free manually
 *
 * !ASSUMPTIONS! -- 
 *
 * 	The checkpoint file and new sampler must have the same total number of chains, dimension, and the histories are the same length
 *
 * 	The temperatures are in ascending order, as output in write_checkpoint
 */
void copy_base_checkpoint_properties(std::string check_file,sampler *samplerptr)
{
	std::fstream file_in;
	file_in.open(check_file, std::ios::in);
	std::string line;
	std::string item;
	int i;
	double temps_old[samplerptr->chain_N];
	double ***history_old_pos=allocate_3D_array(samplerptr->chain_N, samplerptr->history_length, samplerptr->max_dim);
	int ***history_old_status=allocate_3D_array_int(samplerptr->chain_N, samplerptr->history_length, samplerptr->max_dim);
	//double **old_gauss_width=allocate_2D_array(samplerptr->chain_N,samplerptr->types_of_steps);
	double ***old_gauss_width=new double**[samplerptr->chain_N];//,samplerptr->types_of_steps);
	int **old_gauss_width_number=new int*[samplerptr->chain_N];//,samplerptr->types_of_steps);
	double **old_initial_pos = allocate_2D_array(samplerptr->chain_N,samplerptr->max_dim);
	int **old_initial_status = allocate_2D_array_int(samplerptr->chain_N,samplerptr->max_dim);
	int **old_model_status; 
	old_model_status = allocate_2D_array_int(samplerptr->chain_N,1);
	bool no_history=false;
	if(file_in){
		//First row -- dim, chain_N
		std::getline(file_in,line);
		std::stringstream lineStream(line);
		std::getline(lineStream, item, ',');
		std::getline(lineStream, item, ',');
		std::getline(lineStream, item, ',');
		std::getline(lineStream, item, ',');
		//Second Row -- temps
		std::getline(file_in,line);
		std::stringstream lineStreamtemps(line);
		i=0;
		while(std::getline(lineStreamtemps, item, ',')){
			temps_old[i] = std::stod(item);
			i++;
		}

		//third row+chain_N -- step widths
		for(int j =0 ;j<samplerptr->chain_N; j++){
			old_gauss_width[j]= new double*[samplerptr->types_of_steps];
			old_gauss_width_number[j]= new int[samplerptr->types_of_steps];
			i=0;
			std::getline(file_in,line);
			std::stringstream lineStreamwidths(line);
			while(std::getline(lineStreamwidths, item, ',')){
				old_gauss_width_number[j][i] = std::stoi(item);
				i++;
			}
			for(int k = 0 ; k < samplerptr->types_of_steps; k ++){
				old_gauss_width[j][k]= new double[old_gauss_width_number[j][i]];
				i=0;
				std::getline(file_in,line);
				std::stringstream lineStreamwidths(line);
				while(std::getline(lineStreamwidths, item, ',')){
					old_gauss_width[j][k][i] = std::stod(item);
					i++;
				}
			}
		}

		//row 3 +chain_N  to 3+2 chain_N-- initial positions
		for(int j =0 ;j<samplerptr->chain_N; j++){
			i=0;
			std::getline(file_in,line);
			std::stringstream lineStreampos(line);
			while(i<samplerptr->max_dim  ){
				std::getline(lineStreampos, item, ',');
				old_initial_pos[j][i] = std::stod(item);
				i++;
			}
			i=0;
			while(std::getline(lineStreampos, item, ',')){
				old_initial_status[j][i] = std::stod(item);
				i++;
			}
			i=0;
			old_model_status[j][i] = std::stoi(item);
		}
		std::getline(file_in,line);
		std::stringstream lineStreamprimed(line);
		std::getline(lineStreamprimed, item, ',');
		if(item =="true"){
			for(int j =0 ;j<samplerptr->chain_N; j++){
				std::getline(file_in,line);
				std::stringstream lineStreamhist(line);
				i=0;
				while(std::getline(lineStreamhist,item,',')){	
					int step = i/samplerptr->max_dim ;
					int pos = i%samplerptr->max_dim;
					history_old_pos[j][step][pos] = std::stod(item);	
					i++;
				}
			}
		}
		else{
			no_history=true;
		}
	
		//exit(1);
		
	}
	else{std::cout<<"ERROR -- File "<<check_file<<" not found"<<std::endl; exit(1);}
	file_in.close();
	
	int old_ensemble_chain_num = 1 ;
	for(int i = 1 ; i<samplerptr->chain_N; i++){
		if(fabs(temps_old[i] -1) > DOUBLE_COMP_THRESH){
			old_ensemble_chain_num++;	
		}
		else{
			break;
		}
	}
	int new_ensemble_chain_num = 1 ;
	for(int i = 1 ; i<samplerptr->chain_N; i++){
		if(fabs(samplerptr->chain_temps[i] -1) > DOUBLE_COMP_THRESH){
			new_ensemble_chain_num++;	
		}
		else{
			break;
		}
	}
	int wrap_number = new_ensemble_chain_num - old_ensemble_chain_num;
	int new_index = new_ensemble_chain_num,old_index=old_ensemble_chain_num,cp_index=old_ensemble_chain_num ;

	bool utility_bool=true;
	//std::cout<<new_ensemble_chain_num<<" "<<old_ensemble_chain_num<<std::endl;
	if(new_ensemble_chain_num != samplerptr->chain_N 
		&& old_ensemble_chain_num != samplerptr->chain_N){

		for(int i = new_ensemble_chain_num ; i<samplerptr->chain_N; i++){
			//More chains in old ensemble
			if(wrap_number < 0){
				if(old_index%old_ensemble_chain_num==0){
					cp_index = old_index;
				}
				else if(old_index%old_ensemble_chain_num 
					>= new_ensemble_chain_num + wrap_number ){
					if(utility_bool){
						cp_index = old_index;
					}
					else{
						if(old_index-wrap_number<samplerptr->chain_N){
							cp_index = old_index-wrap_number;
						}
						else{
							cp_index = old_index- old_ensemble_chain_num-wrap_number;
						}
					}	
				}
				else{
					cp_index = old_index;
				}
				old_index++;	
		
			}
			//More chains in new ensemble
			else if(wrap_number > 0){
				if(old_index%old_ensemble_chain_num 
					>= old_ensemble_chain_num - wrap_number ){
					if(utility_bool){
						cp_index = old_index;
						utility_bool=false;
					}
					else{
						cp_index = old_index;
						utility_bool=true;
						old_index++;
					}
				}
				else{
					utility_bool=true;
					cp_index = old_index;
					old_index++;	
				}
			}	
			//Chain number didn't change, straight copy
			else{
				cp_index = old_index;
				old_index++;	
			}
			//###########################################################
			//copy old values at cp_index into chain i of new ensemble
			for(int j = 0 ; j<samplerptr->max_dim; j++){
				samplerptr->output[i][samplerptr->chain_pos[i]][j]=old_initial_pos[cp_index][j];
				samplerptr->param_status[i][samplerptr->chain_pos[i]][j]=old_initial_status[cp_index][j];
			}
			samplerptr->model_status[i][samplerptr->chain_pos[i]]=old_model_status[cp_index][0];
			//for(int j = 0 ; j<samplerptr->types_of_steps ; j++){
			//	samplerptr->randgauss_width[i][j] = old_gauss_width[cp_index][j];
			//}
			for(int j = 0 ; j<samplerptr->types_of_steps ; j++){
				for(int k = 0 ; k < samplerptr->randgauss_width_number[i][j]; k ++){
					samplerptr->randgauss_width[i][j][k] = old_gauss_width[cp_index][j][k];
				}
			}
			if(!no_history){
				for(int j =0  ;j<samplerptr->history_length; j++){
					for(int k = 0 ; k<samplerptr->max_dim; k++){
						samplerptr->history[i][j][k]=history_old_pos[cp_index][j][k];
					}
				}
			}
			//###########################################################
			
			
			new_index++;	
			//If new_index is back at T==0, bump old_index to nearest T==0
			if(new_index%new_ensemble_chain_num ==0){ 
				if(old_index%old_ensemble_chain_num !=0){
					old_index += 
						(old_ensemble_chain_num - old_index%old_ensemble_chain_num );
				}
				if(utility_bool)utility_bool = false;
				else utility_bool = true;
			}
			//If we hit the end, restart old_index
			if(old_index>samplerptr->chain_N-1){
				old_index=old_index%old_ensemble_chain_num;
			}
		}
	}

	//Cleanup
	deallocate_3D_array(history_old_pos,samplerptr->chain_N,samplerptr->history_length,samplerptr->max_dim);
	deallocate_3D_array(history_old_status,samplerptr->chain_N,samplerptr->history_length,samplerptr->max_dim);
	//deallocate_2D_array(old_gauss_width, samplerptr->chain_N,samplerptr->types_of_steps);
	for(int i = 0 ; i<samplerptr->chain_N; i ++){
		for(int j = 0  ; j<samplerptr->types_of_steps; j++){
			delete [] old_gauss_width[i][j];
		}
		delete [] old_gauss_width[i];
		delete [] old_gauss_width_number[i];
	}
	delete [] old_gauss_width_number;
	delete [] old_gauss_width;
	deallocate_2D_array(old_initial_pos,samplerptr->chain_N, samplerptr->max_dim);
	deallocate_2D_array(old_initial_status,samplerptr->chain_N, samplerptr->max_dim);
	if(old_model_status){
		deallocate_2D_array(old_model_status,samplerptr->chain_N,1);
	}

}
/*! \brief Utility to write out the parameters and status of a sampler to a file
 */
void write_output_file(std::string file, int step_num, int max_dimension, double ***output, int ***status,int **model_status, int chain_N,int nested_model_number,double *temps,bool RJ)
{
	std::ofstream out_file;
	out_file.open(file);
	out_file.precision(15);
	int coldchains=0;
	int cold_chain_ids[chain_N];
	for(int k = 0 ; k<chain_N; k++){
		if(temps[k] == 1){
			cold_chain_ids[coldchains]=k;
			coldchains++;
		}	
	}
	//Loop through chains
	for(int i =0; i<step_num; i++){
		//if RJ, write out status too
		if(RJ){
			for(int k = 0 ; k<coldchains; k++){
				for(int j=0; j<max_dimension;j++){
					out_file<<output[cold_chain_ids[k]][i][j]<<" , ";
				}
				for(int j = 0 ; j<max_dimension; j++){
					//if(j==max_dimension-1)
					if(j==max_dimension-1 && nested_model_number >0)
						out_file<<status[cold_chain_ids[k]][i][j]<<std::endl;
					else
						out_file<<status[cold_chain_ids[k]][i][j]<<" , ";
				}
				out_file<<model_status[cold_chain_ids[k]][i]<<std::endl;
			}
		}
		//Else, just parameters
		else{
			for(int k = 0 ; k<coldchains; k++){
				for(int j=0; j<max_dimension;j++){
					if(j == max_dimension -1 ){
						out_file<<output[cold_chain_ids[k]][i][j]<<std::endl;
					}
					else{
						out_file<<output[cold_chain_ids[k]][i][j]<<" , ";
					}
				}
			}
		}
	}
	out_file.close();

}

int count_cold_chains(double *temps, int chain_N)
{
	int coldchains=0;
	for(int k = 0 ; k<chain_N; k++){
		if(temps[k] == 1){
			coldchains++;
		}	
	}
	return coldchains;

}

/*! \brief Utility to write out the parameters and status of a sampler to a file
 */
void reduce_output(int step_num, int max_dimension, double ***output_old, int ***status_old,int **model_status_old,double **output_new, int **status_new,int *model_status_new,int chain_N,int nested_model_number,double *temps,bool RJ)
{
	int coldchains=0;
	int cold_chain_ids[chain_N];
	for(int k = 0 ; k<chain_N; k++){
		if(temps[k] == 1){
			cold_chain_ids[coldchains]=k;
			coldchains++;
		}	
	}
	//int final_length = step_num*coldchains;
	for(int i =0; i<step_num; i++){
		for(int j = 0 ; j<coldchains; j++){
			for(int k = 0 ; k<max_dimension; k++){
				output_new[i*coldchains + j][k] = output_old[cold_chain_ids[j]][i][k];
				if(RJ){
					status_new[i*coldchains + j][k] = 
						status_old[cold_chain_ids[j]][i][k];
				}
			}
			if(RJ){
				model_status_new[i*coldchains + j] = 
					model_status_old[cold_chain_ids[j]][i];
			}
	
		}
	}

}

//##################################################################
//Evidence Calculation
//##################################################################
void integrate_likelihood(sampler *samplerptr)
{
		
	for(int i = 0 ; i<samplerptr->chain_N; i++){
		samplerptr->thermodynamic_integrated_likelihood[i] =0;	
		samplerptr->thermodynamic_integrated_likelihood_terms[i] =0;	
	}

	for(int i = 0 ; i<samplerptr->chain_N; i++){
		int pos = samplerptr->chain_pos[i];
		if(samplerptr->restarted_chain[i]){
			pos = samplerptr->N_steps;
		}
		for(int j = 0 ; j<pos; j++){
			samplerptr->thermodynamic_integrated_likelihood[i] += samplerptr->ll_lp_output[i][j][0];
		}
		samplerptr->thermodynamic_integrated_likelihood_terms[i]+=pos;
		samplerptr->thermodynamic_integrated_likelihood[i] /= samplerptr->thermodynamic_integrated_likelihood_terms[i];
	}
	
}


void combine_chain_evidence(sampler *samplerptr, double *evidences, int *total_terms, int ensemble_size)
{
	for(int i = 0 ; i<ensemble_size; i++){
		evidences[i] = 0;		
		total_terms[i] = 0;		
	}
	for(int i = 0 ; i<samplerptr->chain_N; i++){
		evidences[ i%ensemble_size]+= samplerptr->thermodynamic_integrated_likelihood[i]*samplerptr->thermodynamic_integrated_likelihood_terms[i];
		total_terms[i%ensemble_size] += samplerptr->thermodynamic_integrated_likelihood_terms[i];
	}
	for(int i = 0 ; i<ensemble_size; i++){
		evidences[i] /= total_terms[i];		
	}
	
	return ;
}

//##################################################################
//Simple solution to integration issue
//##################################################################
struct helper_params
{
	gsl_interp_accel *a;	
	gsl_spline *s;	
};
double integration_helper(double beta, void *params)
{
	helper_params *p = (helper_params *) params;
	double likelihood =  gsl_spline_eval(p->s, beta, p->a);
	return likelihood;
}
//##################################################################

int thermodynamic_integration(double *integrated_likelihoods,double *temps,int temps_N, double *evidence, double *error)
{
	double betas[temps_N];
	double temp_integrated_likelihoods[temps_N];
	//Reverse the order because interpolation needs beta going up
	for(int i = 0 ; i<temps_N; i++){
		betas[i] = 1./temps[temps_N - i-1];
		temp_integrated_likelihoods[i] = integrated_likelihoods[temps_N - i-1];
	}
	gsl_interp_accel *acc = gsl_interp_accel_alloc();	
	gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, temps_N);
	gsl_spline_init(spline, betas, temp_integrated_likelihoods, temps_N);

	helper_params params;
	params.a = acc;
	params.s = spline;

	gsl_function F;
	F.function = &integration_helper;
	F.params = &params;

	gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

	int errorcode = gsl_integration_qags(&F, betas[0],betas[temps_N-1],0,1e-7,1000,w, evidence, error);

	gsl_integration_workspace_free(w);
	gsl_spline_free(spline);
	gsl_interp_accel_free(acc);
	return errorcode;
}

double calculate_evidence(sampler *samplerptr)
{
	
	integrate_likelihood(samplerptr);
	double integrated_likelihoods[samplerptr->chain_N];
	int total_terms[samplerptr->chain_N];
	int temps_N = (int)(samplerptr->chain_N /count_cold_chains(samplerptr->chain_temps, samplerptr->chain_N));
	combine_chain_evidence(samplerptr, integrated_likelihoods, total_terms,temps_N);
	double evidence,error;
	int errcode = thermodynamic_integration(integrated_likelihoods,samplerptr->chain_temps,temps_N, &evidence, &error);
	
	return evidence;
}

//##################################################################
//##################################################################
void assign_ensemble_temps(double *chain_temps, int chain_N,int max_chain_N_thermo_ensemble,double TMAX)
{
	//Geometric
	//double temp_step = pow(TMAX/1.,1./max_chain_N_thermo_ensemble);
	////Initial temps
	//for(int i = 0 ; i<chain_N;i++){
	//	if(i % max_chain_N_thermo_ensemble == 0 ){
	//		chain_temps[i] = 1;
	//	}
	//	else if( (i+1) % max_chain_N_thermo_ensemble  == 0 ){
	//		chain_temps[i] =TMAX;
	//	}
	//	else{
	//		chain_temps[i] = temp_step* chain_temps[i-1];	
	//	}
	//}
	//Linear in beta
	//Initial temps
	double delta_beta = (1-1./TMAX)/(max_chain_N_thermo_ensemble-1);
	for(int i = 0 ; i<chain_N;i++){
		if(i % max_chain_N_thermo_ensemble == 0 ){
			chain_temps[i] = 1;
		}
		else if( (i+1) % max_chain_N_thermo_ensemble  == 0 ){
			chain_temps[i] =TMAX;
		}
		else{
			chain_temps[i] = 1./(1 - (i%max_chain_N_thermo_ensemble)*delta_beta);	
		}
	}
	return;
}

