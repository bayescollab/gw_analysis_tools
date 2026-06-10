#include "mcmc_io_util.h"
#include "mcmc_sampler_internals.h"
#include "util.h"
#include "autocorrelation.h"
#include <iostream>
#include <vector>

#ifdef _HDF5
#include <H5Cpp.h>
#endif


//#############################################################
//mcmc_sampler_output definitions
//#############################################################
mcmc_sampler_output::mcmc_sampler_output( int chain_N, int dim, int nested_model_N)
{
	chain_number = chain_N;
	dimension = dim;
	nested_model_number = nested_model_N;
	chain_temperatures = new double[chain_number];
	chain_lengths = new int[chain_number];
	trim_lengths = new int[chain_number];
	for(int i = 0 ;i<chain_number; i++){
		trim_lengths[i]=0;
	}
	file_trim_lengths = new int[chain_number];
	for(int i = 0 ;i<chain_number; i++){
		file_trim_lengths[i]=0;
	}
};
mcmc_sampler_output::~mcmc_sampler_output()
{
	if(chain_temperatures){
		delete [] chain_temperatures;		
		chain_temperatures=NULL;
	}
	if(cold_chain_ids){
		delete [] cold_chain_ids;		
		cold_chain_ids=NULL;		
	}
	dealloc_output();
	dealloc_status();
	dealloc_model_status();
	dealloc_logL_logP();
	dealloc_integrated_likelihoods();
	if(chain_lengths){
		delete [] chain_lengths;
		chain_lengths = NULL;
	}
	if(trim_lengths){
		delete [] trim_lengths;
		trim_lengths = NULL;
	}
	if(file_trim_lengths){
		delete [] file_trim_lengths;
		file_trim_lengths = NULL;
	}
	if(ac_vals){
		for(int i = 0 ; i< cold_chain_number_ac_alloc; i++){
			delete [] ac_vals[i];	
		}
		delete [] ac_vals;	
	}
	if(max_acs){
		delete [] max_acs;
	}
	if(dump_files.size() != 0){
		for(int i = 0 ; i<dump_files.size(); i++){
			if(dump_files[i]->file_trim_lengths){
				delete [] dump_files[i]->file_trim_lengths;
				dump_files[i]->file_trim_lengths = NULL;
			}
			delete dump_files[i];
		}
	}
};
void mcmc_sampler_output::set_trim(int trim){
	for(int i = 0 ; i<chain_number; i++){
		trim_lengths[i]=trim;
	}
}
void mcmc_sampler_output::populate_chain_temperatures(double *temperatures)
{
	for(int i= 0 ; i<chain_number; i++){
		chain_temperatures[i] = temperatures[i];
	}	
	update_cold_chain_list();
}
void mcmc_sampler_output::update_cold_chain_list()
{
	int temp=0;
	int *cold_chain_ids_temp = new int[chain_number];
	for(int i = 0 ; i<chain_number; i++){
		if(fabs(chain_temperatures[i] -1)<DOUBLE_COMP_THRESH)
		{
			cold_chain_ids_temp[temp] = i;
			temp +=1;
		}
	}
	cold_chain_number = temp;
	if(cold_chain_ids)
	{
		delete [] cold_chain_ids;
		cold_chain_ids = NULL;
	}
	cold_chain_ids = new int[temp];
	for(int i = 0 ; i<temp ; i++){
		cold_chain_ids[i]=cold_chain_ids_temp[i];
	}
	delete [] cold_chain_ids_temp;
	cold_chain_ids_temp = NULL;
}
void mcmc_sampler_output::populate_initial_output(double ***new_output, int ***new_status,int **new_model_status,double ***new_logL_logP,int *chain_positions)
{
	dealloc_output();	
	dealloc_logL_logP();	
	if(RJ){
		dealloc_status();	
		dealloc_model_status();	
	}
	output = new double**[chain_number];
	logL_logP = new double**[chain_number];
	for(int i = 0 ;i<chain_number; i++){
		chain_lengths[i]=chain_positions[i];
		output[i]=new double*[chain_positions[i]];
		logL_logP[i]=new double*[chain_positions[i]];
		for(int j =0 ; j<chain_positions[i];j++){
			logL_logP[i][j] = new double[2];
			output[i][j]=new double[dimension];
			for(int k =0 ; k<dimension; k++){
				output[i][j][k]=new_output[i][j][k];
			}
			logL_logP[i][j][0] = new_logL_logP[i][j][0];
			logL_logP[i][j][1] = new_logL_logP[i][j][1];
		}
	}
	if(RJ){
		status = new int**[chain_number];
		for(int i = 0 ;i<chain_number; i++){
			status[i] = new int*[chain_positions[i]];	
			for(int j =0 ; j<chain_positions[i];j++){
				status[i][j] = new int[dimension];	
				for(int k =0 ; k<dimension; k++){
					status[i][j][k]=new_status[i][j][k];
				}
			}
		}
		model_status = new int*[chain_number];
		for(int i = 0 ;i<chain_number; i++){
			model_status[i] = new int[chain_positions[i]];	
			for(int j =0 ; j<chain_positions[i];j++){
				model_status[i][j]=new_model_status[i][j];
			}
		}
	}
}

void mcmc_sampler_output::append_to_output(double ***new_output,int ***new_status,int **new_model_status, double ***new_logL_logP, int *chain_positions)
{
	int *new_lengths= new int[chain_number];
	for(int i = 0 ; i<chain_number; i++){
		new_lengths[i]=chain_lengths[i]+chain_positions[i];
	}
	//Copy all values into new temp array
	double ***new_total_output = new double**[chain_number];
	double ***new_total_logL_logP = new double**[chain_number];
	for(int i = 0 ; i<chain_number ; i++){
		new_total_output[i] = new double*[new_lengths[i]];
		new_total_logL_logP[i] = new double*[new_lengths[i]];
		for(int j = 0 ; j<new_lengths[i]; j++){
			new_total_output[i][j] = new double[dimension];
			new_total_logL_logP[i][j] = new double[2];
			if(j <chain_lengths[i]){
				for (int k = 0 ; k<dimension ; k++){
					new_total_output[i][j][k] = output[i][j][k];
				}
				new_total_logL_logP[i][j][0] = logL_logP[i][j][0];
				new_total_logL_logP[i][j][1] = logL_logP[i][j][1];
			}
			else{
				for (int k = 0 ; k<dimension ; k++){
					new_total_output[i][j][k] 
						= new_output[i][j - chain_lengths[i]][k];
				}
				new_total_logL_logP[i][j][0] = 
					new_logL_logP[i][j - chain_lengths[i]][0];
				new_total_logL_logP[i][j][1] = 
					new_logL_logP[i][j - chain_lengths[i]][1];
			}
		}
	}
	int ***new_total_status=NULL;
	int **new_total_model_status=NULL;
	if(RJ){
		new_total_status = new int**[chain_number];
		for(int i = 0 ; i<chain_number ; i++){
			new_total_status[i] = new int*[new_lengths[i]];
			for(int j = 0 ; j<new_lengths[i]; j++){
				new_total_status[i][j] = new int[dimension];
				if(j <chain_lengths[i]){
					for (int k = 0 ; k<dimension ; k++){
						new_total_status[i][j][k] = status[i][j][k];
					}
				}
				else{
					for (int k = 0 ; k<dimension ; k++){
						new_total_status[i][j][k] 
							= new_status[i][j - chain_lengths[i]][k];
					}
				}
			}
		}
		new_total_model_status = new int*[chain_number];
		for(int i = 0 ; i<chain_number ; i++){
			new_total_model_status[i] = new int[new_lengths[i]];
			for(int j = 0 ; j<new_lengths[i]; j++){
				if(j <chain_lengths[i]){
					new_total_model_status[i][j] = model_status[i][j];
				}
				else{
					new_total_model_status[i][j] = new_model_status[i][j - chain_lengths[i]];
				}
			}
		}
	}
	//deallocate and move values into output
	dealloc_output();
	dealloc_logL_logP();
	if(RJ){
		dealloc_status();
		dealloc_model_status();
	}
	output = new double**[chain_number];
	logL_logP = new double**[chain_number];
	for(int i = 0 ; i<chain_number ; i++){
		chain_lengths[i]=new_lengths[i];
		output[i] = new double*[new_lengths[i]];
		logL_logP[i] = new double*[new_lengths[i]];
		for(int j = 0 ; j<new_lengths[i]; j++){
			output[i][j] = new double[dimension];
			logL_logP[i][j] = new double[2];
			for (int k = 0 ; k<dimension ; k++){
				output[i][j][k] = new_total_output[i][j][k];
			}
			logL_logP[i][j][0] = new_total_logL_logP[i][j][0];
			logL_logP[i][j][1] = new_total_logL_logP[i][j][1];
		}
	}
	for(int j = 0 ; j<chain_number;j ++){
		for(int i = 0 ; i<chain_lengths[j];i ++){
			delete [] new_total_output[j][i];		
			delete [] new_total_logL_logP[j][i];		
		}
		delete [] new_total_output[j];		
		delete [] new_total_logL_logP[j];		
	}
	delete [] new_total_output;		
	delete [] new_total_logL_logP;		
	new_total_output= NULL;
	new_total_logL_logP= NULL;
	delete [] new_lengths;
	new_lengths = NULL;
	if(RJ){
		status = new int**[chain_number];
		for(int i = 0 ; i<chain_number ; i++){
			status[i] = new int*[chain_lengths[i]];
			for(int j = 0 ; j<chain_lengths[i]; j++){
				status[i][j] = new int[dimension];
				for (int k = 0 ; k<dimension ; k++){
					status[i][j][k] = new_total_status[i][j][k];
				}
			}
		}
		for(int j = 0 ; j<chain_number;j ++){
			for(int i = 0 ; i<chain_lengths[j];i ++){
				delete [] new_total_status[j][i];		
			}
			delete [] new_total_status[j];		
		}
		delete [] new_total_status;		
		new_total_status= NULL;

		model_status = new int*[chain_number];
		for(int i = 0 ; i<chain_number ; i++){
			model_status[i] = new int[chain_lengths[i]];
			for(int j = 0 ; j<chain_lengths[i]; j++){
				model_status[i][j] = new_total_model_status[i][j];
			}
		}
		for(int j = 0 ; j<chain_number;j ++){
			delete [] new_total_model_status[j];		
		}
		delete [] new_total_model_status;		
		new_total_model_status= NULL;
		
	}
}
void mcmc_sampler_output::dealloc_output()
{
	if(output){
		for(int j = 0 ; j<chain_number;j ++){
			for(int i = 0 ; i<chain_lengths[j];i ++){
				delete [] output[j][i];		
			}
			delete [] output[j];		
		}
		delete [] output;		
		output= NULL;
	}
}
void mcmc_sampler_output::dealloc_status()
{
	if(status){
		for(int j = 0 ; j<chain_number;j ++){
			for(int i = 0 ; i<chain_lengths[j];i ++){
				delete [] status[j][i];		
			}
			delete [] status[j];		
		}
		delete [] status;		
		status= NULL;
	}
}
void mcmc_sampler_output::dealloc_model_status()
{
	if(model_status){
		for(int j = 0 ; j<chain_number;j ++){
			delete [] model_status[j];		
		}
		delete [] model_status;		
		model_status= NULL;
	}
}
void mcmc_sampler_output::dealloc_logL_logP()
{
	if(logL_logP){
		for(int j = 0 ; j<chain_number;j ++){
			for(int i = 0 ; i<chain_lengths[j];i ++){
				delete [] logL_logP[j][i];	
			}
			delete [] logL_logP[j];		
		}
		delete [] logL_logP;		
		logL_logP= NULL;
	}
}

void mcmc_sampler_output::calc_ac_vals_full_ensemble(bool trim)
{
	debugger_print(__FILE__,__LINE__,"Not currently working! Under development");
	if(ac_vals_full_ensemble){
		for(int i = 0 ; i<cold_chain_number_ac_alloc; i++){
			delete [] ac_vals_full_ensemble[i];
		}
		delete [] ac_vals_full_ensemble;
	}
	ac_vals_full_ensemble = new int*[cold_chain_number];
	for(int i = 0 ; i<cold_chain_number; i++){
		ac_vals_full_ensemble[i] = new int[dimension];	
	}

	int segments = 1;
	int ***temp = new int**[cold_chain_number];
	for(int i = 0 ; i<cold_chain_number; i++){
		temp[i] = new int*[dimension];
		for(int j = 0 ; j<dimension; j++){
			temp[i][j]=new int[segments];
		}
	}
	double ***temp_chains = new double**[cold_chain_number];
	int beginning_id = 0;
	for (int i = 0 ;i<cold_chain_number; i++){
		int id  = cold_chain_ids[i];
		beginning_id = 0;
		if(trim){
			beginning_id = trim_lengths[id];
		}
		temp_chains[i] = new double*[chain_lengths[id]-beginning_id];
		for(int j = beginning_id ; j<chain_lengths[id]; j++){
			temp_chains[i][j-beginning_id] = new double[dimension];
			for(int k = 0 ; k<dimension  ; k++){
				temp_chains[i][j-beginning_id][k] = output[id][j][k];
			}
		}
	}
	auto_corr_from_data_batch(temp_chains, chain_lengths[cold_chain_ids[0]]-beginning_id, dimension, cold_chain_number,temp, segments, target_correlation, threads, true);
	
	for(int i = 0 ;i<cold_chain_number; i++){
		int id  = cold_chain_ids[i];
		if(trim){
			beginning_id = trim_lengths[id];
		}
		else{
			beginning_id = 0;
		}
		for(int j = beginning_id ; j<chain_lengths[id]; j ++){
			delete [] temp_chains[i][j-beginning_id];
		}
		delete [] temp_chains[i];
	}
	delete [] temp_chains;
	temp_chains=NULL;
	for(int i = 0 ; i<cold_chain_number; i++){
		for(int j = 0 ; j<dimension; j++){
			ac_vals_full_ensemble[i][j]= temp[i][j][segments-1];
		}
	}
	for(int i = 0 ; i<cold_chain_number; i++){
		for(int j = 0  ; j<dimension; j++){
			delete [] temp[i][j];
		}
		delete [] temp[i];
	}
	delete [] temp;
	temp = NULL;	


	if(max_acs){
		delete [] max_acs;
	}
	max_acs = new int[cold_chain_number];
	for(int i = 0 ; i<cold_chain_number; i ++){
		int max_ac=0;
		for(int j = 0 ; j<dimension; j++){
			if(ac_vals_full_ensemble[i][j]>max_ac){
				max_ac = ac_vals_full_ensemble[i][j];	
			}
		}
		max_acs[i]=max_ac;
	}

}
void mcmc_sampler_output::calc_ac_vals(bool trim)
{
	if(ac_vals){
		for(int i = 0 ; i<cold_chain_number_ac_alloc; i++){
			delete [] ac_vals[i];
		}
		delete [] ac_vals;
	}
	update_cold_chain_list();	
	cold_chain_number_ac_alloc = cold_chain_number;
	
	ac_vals = new int*[cold_chain_number];
	for(int i = 0 ; i<cold_chain_number; i++){
		ac_vals[i] = new int[dimension];	
	}

	int segments = 1;
	int ***temp = new int**[cold_chain_number];
	for(int i = 0 ; i<cold_chain_number; i++){
		temp[i] = new int*[dimension];
		for(int j = 0 ; j<dimension; j++){
			temp[i][j]=new int[segments];
		}
	}
	double ***temp_chains = new double**[cold_chain_number];
	int beginning_id = 0;
	for (int i = 0 ;i<cold_chain_number; i++){
		int id  = cold_chain_ids[i];
		beginning_id = 0;
		if(trim){
			beginning_id = trim_lengths[id];
		}
		temp_chains[i] = new double*[chain_lengths[id]-beginning_id];
		for(int j = beginning_id ; j<chain_lengths[id]; j++){
			temp_chains[i][j-beginning_id] = new double[dimension];
			for(int k = 0 ; k<dimension  ; k++){
				temp_chains[i][j-beginning_id][k] = output[id][j][k];
			}
		}
	}
	auto_corr_from_data_batch(temp_chains, chain_lengths[cold_chain_ids[0]]-beginning_id, dimension, cold_chain_number,temp, segments, target_correlation, threads, true);
	
	for(int i = 0 ;i<cold_chain_number; i++){
		int id  = cold_chain_ids[i];
		if(trim){
			beginning_id = trim_lengths[id];
		}
		else{
			beginning_id = 0;
		}
		for(int j = beginning_id ; j<chain_lengths[id]; j ++){
			delete [] temp_chains[i][j-beginning_id];
		}
		delete [] temp_chains[i];
	}
	delete [] temp_chains;
	temp_chains=NULL;
	for(int i = 0 ; i<cold_chain_number; i++){
		for(int j = 0 ; j<dimension; j++){
			ac_vals[i][j]= temp[i][j][segments-1];
		}
	}
	for(int i = 0 ; i<cold_chain_number; i++){
		for(int j = 0  ; j<dimension; j++){
			delete [] temp[i][j];
		}
		delete [] temp[i];
	}
	delete [] temp;
	temp = NULL;	


	if(max_acs){
		delete [] max_acs;
	}
	max_acs = new int[cold_chain_number];
	for(int i = 0 ; i<cold_chain_number; i ++){
		int max_ac=0;
		for(int j = 0 ; j<dimension; j++){
			if(ac_vals[i][j]>max_ac){
				max_ac = ac_vals[i][j];	
			}
		}
		max_acs[i]=max_ac;
	}

}
void mcmc_sampler_output::count_indep_samples(bool trim)
{
	indep_samples = 0;
	double mean_ac=0;
	double mean_pos=0;
	for(int i = 0 ; i<cold_chain_number; i ++){
		int id = cold_chain_ids[i];
		int max_ac=1;
		for(int j = 0 ; j<dimension; j++){
			if(ac_vals[i][j]>max_ac){
				max_ac = ac_vals[i][j];	
			}
		}
		mean_ac+=max_ac;
		if(trim){
			mean_pos+=(chain_lengths[id]-trim_lengths[id]);
		}
		else{
			mean_pos+=(chain_lengths[id]);
		}
		//if(trim){
		//	indep_samples += (chain_lengths[id]-trim_lengths[id])/max_ac;	
		//}
		//else{
		//	indep_samples += (chain_lengths[id])/max_ac;
		//}
	}
	mean_ac/=cold_chain_number;
	mean_pos/=cold_chain_number;
	indep_samples= mean_pos/mean_ac;
}

//Use HDF5 if available
#ifdef _HDF5
int mcmc_sampler_output::write_flat_thin_output(std::string filename, bool use_stored_ac, bool trim)
{
	try{
		if(!use_stored_ac || !ac_vals){
			calc_ac_vals(trim);	
			count_indep_samples(trim);
		}
		double *flattened= new double[indep_samples*dimension];
		int ct = 0;
		for(int i =0  ;i<cold_chain_number; i++){
			int beginning_id = 0;
			if(trim){
				beginning_id = trim_lengths[cold_chain_ids[i]];
			}
			for(int j = beginning_id ; j<chain_lengths[cold_chain_ids[i]]; j++){
				if(j%max_acs[i] == 0 && ct<indep_samples){
					for(int k = 0 ; k<dimension; k++){
						flattened[ct*dimension + k ] 
							= output[cold_chain_ids[i]][j][k];		
					}
					ct++;
				}
			}
		}
		//#################################################################
		std::string FILE_NAME(filename+".hdf5");

		H5::H5File file(FILE_NAME,H5F_ACC_TRUNC);
		H5::Group output_group(file.createGroup("/THINNED_MCMC_OUTPUT"));
		H5::DataSpace *dataspace=NULL ;
		H5::DataSet *dataset=NULL;
		H5::DSetCreatPropList *plist=NULL;
		hsize_t chunk_dims[2] = {(hsize_t)chunk_steps,(hsize_t)dimension};
		int RANK=2;
		hsize_t dims[RANK];
		dims[0]= indep_samples;
		dims[1]= dimension;

		if(chunk_steps>dims[0]){chunk_dims[0] = dims[0];}
		else{chunk_dims[0] = chunk_steps;}

		dataspace = new H5::DataSpace(RANK,dims);
	
		plist = new H5::DSetCreatPropList;
		plist->setChunk(2,chunk_dims);
		plist->setDeflate(6);

		dataset = new H5::DataSet(
			output_group.createDataSet("THINNED FLATTENED CHAINS",
				H5::PredType::NATIVE_DOUBLE,*dataspace,*plist)
			);

		dataset->write(flattened, H5::PredType::NATIVE_DOUBLE);	
		//Cleanup
		delete dataset;
		delete dataspace;
		delete plist;
	
		//Cleanup
		output_group.close();
		delete [] flattened;

	}	
	catch( H5::FileIException error )
	{
		error.printErrorStack();
		return -1;
	}
	// catch failure caused by the DataSet operations
	catch( H5::DataSetIException error )
	{
		error.printErrorStack();
	   	return -1;
	}
	// catch failure caused by the DataSpace operations
	catch( H5::DataSpaceIException error )
	{
		error.printErrorStack();
	   	return -1;
	}
	// catch failure caused by the DataSpace operations
	catch( H5::DataTypeIException error )
	{
		error.printErrorStack();
	   	return -1;
	}
	return 0;
}
int mcmc_sampler_output::create_data_dump(bool cold_only, bool trim,std::string filename)
{
	int file_id = 0;
	bool found = false;
	if(dump_files.size() != 0){
		for(int i = 0 ; i<dump_file_names.size(); i++){
			if( filename == dump_file_names[i]){
				found = true;
				file_id = i;
			}
		}	
	}
	if(!found ){
		file_id = dump_files.size();	
		dump_file_struct *new_dump_file = new dump_file_struct;
		dump_files.push_back(new_dump_file);
		dump_files[file_id]->file_trim_lengths = new int[chain_number];
		dump_file_names.push_back(filename);
	}
	dump_files[file_id]->cold_only = cold_only;
	
	if(trim){
		dump_files[file_id]->trimmed = true;
		for(int i= 0 ; i<chain_number; i++){
			dump_files[file_id]->file_trim_lengths[i]=trim_lengths[i];
		}
	}
	else{
		dump_files[file_id]->trimmed = false;
	}
	try{
		std::string FILE_NAME(filename);
		int chains;
		int *ids=NULL;
		if(cold_only){
			ids = cold_chain_ids;
			chains = cold_chain_number;
		}
		else{
			chains = chain_number;
			ids = new int[chain_number];
			for(int i = 0  ; i<chain_number; i++){
				ids[i]=i;
			}
		}
		H5::H5File file(FILE_NAME,H5F_ACC_TRUNC);
		H5::Group output_group(file.createGroup("/MCMC_OUTPUT"));
		H5::Group output_LL_LP_group(file.createGroup("/MCMC_OUTPUT/LOGL_LOGP"));
		H5::Group status_group;
		H5::Group model_status_group;
		if(RJ){
			status_group = H5::Group(file.createGroup("/MCMC_OUTPUT/STATUS"));
			model_status_group = H5::Group(file.createGroup("/MCMC_OUTPUT/MODEL_STATUS"));
		}
		H5::Group meta_group(file.createGroup("/MCMC_METADATA"));
		double *temp_buffer=NULL;
		double *temp_ll_lp_buffer=NULL;
		int *temp_status_buffer=NULL;
		int *temp_model_status_buffer=NULL;
		H5::DataSpace *dataspace=NULL ;
		H5::DataSpace *dataspace_ll_lp=NULL ;
		H5::DataSpace *dataspace_status=NULL ;
		H5::DataSpace *dataspace_model_status=NULL ;
		H5::DataSet *dataset=NULL;
		H5::DataSet *dataset_ll_lp=NULL;
		H5::DataSet *dataset_status=NULL;
		H5::DataSet *dataset_model_status=NULL;
		H5::DSetCreatPropList *plist=NULL;
		H5::DSetCreatPropList *plist_ll_lp=NULL;
		H5::DSetCreatPropList *plist_status=NULL;
		H5::DSetCreatPropList *plist_model_status=NULL;
		hsize_t chunk_dims[2] = {(hsize_t)chunk_steps,(hsize_t)dimension};
		hsize_t chunk_dims_ll_lp[2] = {(hsize_t)chunk_steps,2};
		hsize_t chunk_dims_status[2] = {(hsize_t)chunk_steps,2};
		hsize_t chunk_dims_model_status[2] = {(hsize_t)chunk_steps,2};
		hsize_t max_dims[2] = {H5S_UNLIMITED,H5S_UNLIMITED};
		for(int i = 0 ; i<chains; i++){
			int RANK=2;
			hsize_t dims[RANK];
			hsize_t dims_ll_lp[RANK];
			hsize_t dims_status[RANK];
			//hsize_t dims_model_status[RANK];
			if(trim){
				dims[0]= chain_lengths[ids[i]]-trim_lengths[ids[i]];
				dims_ll_lp[0]= chain_lengths[ids[i]]-trim_lengths[ids[i]];
				dims_status[0]= chain_lengths[ids[i]]-trim_lengths[ids[i]];
				//dims_model_status[0]= chain_lengths[ids[i]]-trim_lengths[ids[i]];
			}
			else{
				dims[0]= chain_lengths[ids[i]];
				dims_ll_lp[0]= chain_lengths[ids[i]];
				dims_status[0]= chain_lengths[ids[i]];
				//dims_model_status[0]= chain_lengths[ids[i]];
			}
			dims[1]= dimension;
			dims_ll_lp[1]= 2;
			dims_status[1]= dimension;

			if(chunk_steps>dims[0]){chunk_dims_ll_lp[0]=dims[0];chunk_dims[0] = dims[0];}
			else{chunk_dims_ll_lp[0]=chunk_steps;chunk_dims[0] = chunk_steps;}

			dataspace = new H5::DataSpace(RANK,dims,max_dims);
			dataspace_ll_lp = new H5::DataSpace(RANK,dims_ll_lp,max_dims);
			if(RJ){
				dataspace_status = new H5::DataSpace(RANK,dims_status,max_dims);
			}
	
			plist = new H5::DSetCreatPropList;
			plist->setChunk(2,chunk_dims);
			plist->setDeflate(6);

			plist_ll_lp = new H5::DSetCreatPropList;
			plist_ll_lp->setChunk(2,chunk_dims_ll_lp);
			plist_ll_lp->setDeflate(6);

			if(RJ){
				plist_status = new H5::DSetCreatPropList;
				plist_status->setChunk(2,chunk_dims_status);
				plist_status->setDeflate(6);

			}

			dataset = new H5::DataSet(
				output_group.createDataSet("CHAIN "+std::to_string(ids[i]),
					H5::PredType::NATIVE_DOUBLE,*dataspace,*plist)
				);
			dataset_ll_lp = new H5::DataSet(
				output_LL_LP_group.createDataSet("CHAIN "+std::to_string(ids[i]),
					H5::PredType::NATIVE_DOUBLE,*dataspace_ll_lp,*plist_ll_lp)
				);
			if(RJ){
				dataset_status = new H5::DataSet(
					status_group.createDataSet("CHAIN "+std::to_string(ids[i]),
						H5::PredType::NATIVE_INT,*dataspace_status,*plist_status)
					);

			}

			temp_buffer = new double[ int(dims[0]*dims[1]) ];
			temp_ll_lp_buffer = new double[ int(dims_ll_lp[0]*dims_ll_lp[1]) ];
			int beginning_id=0;
			if(trim){ beginning_id =trim_lengths[ids[i]];}
			for(int j = 0 ; j<chain_lengths[ids[i]] - beginning_id; j++){
				for(int k = 0 ; k<dimension; k++){
					temp_buffer[j*dimension +k] = output[ids[i]][j+beginning_id][k];	
				}
				temp_ll_lp_buffer[j*2]=logL_logP[ids[i]][j+beginning_id][0];
				temp_ll_lp_buffer[j*2+1]=logL_logP[ids[i]][j+beginning_id][1];
			}
			dataset->write(temp_buffer, H5::PredType::NATIVE_DOUBLE);
			dataset_ll_lp->write(temp_ll_lp_buffer, H5::PredType::NATIVE_DOUBLE);
			if(RJ){
				temp_status_buffer = new int[ int(dims_status[0]*dims_status[1]) ];
				int beginning_id=0;
				if(trim){ beginning_id =trim_lengths[ids[i]];}
				for(int j = 0 ; j<chain_lengths[ids[i]] - beginning_id; j++){
					for(int k = 0 ; k<dimension; k++){
						temp_status_buffer[j*dimension +k] = status[ids[i]][j+beginning_id][k];	
					}
				}
				dataset_status->write(temp_status_buffer, H5::PredType::NATIVE_INT);

			}
			//Cleanup
			delete dataset;
			delete dataset_ll_lp;
			delete dataspace;
			delete dataspace_ll_lp;
			delete plist;
			delete plist_ll_lp;
			delete [] temp_buffer;
			delete [] temp_ll_lp_buffer;
			temp_buffer = NULL;
			temp_ll_lp_buffer = NULL;
			if(RJ){
				delete dataset_status;
				delete dataspace_status;
				delete plist_status;
				delete [] temp_status_buffer;
				temp_status_buffer = NULL;

			}
			
			//TODO -- this section can't be right.. 
			if(RJ){
				int RANK=2;
				hsize_t dims_model_status[RANK];
				if(trim){
					dims_model_status[0]= chain_lengths[ids[i]]-trim_lengths[ids[i]];
				}
				else{
					dims_model_status[0]= chain_lengths[ids[i]];
				}
				dims_model_status[1]= nested_model_number;

				dataspace_model_status = new H5::DataSpace(RANK,dims_model_status,max_dims);
	

				plist_model_status = new H5::DSetCreatPropList;
				plist_model_status->setChunk(2,chunk_dims_model_status);
				plist_model_status->setDeflate(6);

				dataset_model_status = new H5::DataSet(
					model_status_group.createDataSet("CHAIN "+std::to_string(ids[i]),
					H5::PredType::NATIVE_INT,*dataspace_model_status,*plist_model_status)
					);

				temp_model_status_buffer = new int[ int(dims_model_status[0]*dims_model_status[1]) ];
				int beginning_id=0;
				if(trim){ beginning_id =trim_lengths[ids[i]];}
				for(int j = 0 ; j<chain_lengths[ids[i]] - beginning_id; j++){
					temp_model_status_buffer[j*nested_model_number ] = model_status[ids[i]][j+beginning_id];	
				}
				dataset_model_status->write(temp_model_status_buffer, H5::PredType::NATIVE_INT);

				//Cleanup
				delete dataset_model_status;
				delete dataspace_model_status;
				delete plist_model_status;
				delete [] temp_model_status_buffer;
				temp_model_status_buffer = NULL;
			}
			
		}
		//#################################################
		hsize_t dimsT[1];
		dimsT[0]= chain_number;
		dataspace = new H5::DataSpace(1,dimsT);
		dataset = new H5::DataSet(
			meta_group.createDataSet("CHAIN TEMPERATURES",
				H5::PredType::NATIVE_DOUBLE,*dataspace)
			);
		dataset->write(chain_temperatures, H5::PredType::NATIVE_DOUBLE);	
		delete dataset;
		delete dataspace;
		//#################################################
		if(integrated_likelihoods){
			hsize_t dimsIL[1];
			dimsIL[0]= ensemble_size;
			dataspace = new H5::DataSpace(1,dimsIL);
			dataset = new H5::DataSet(
				meta_group.createDataSet("INTEGRATED LIKELIHOODS",
					H5::PredType::NATIVE_DOUBLE,*dataspace)
				);
			dataset->write(integrated_likelihoods, H5::PredType::NATIVE_DOUBLE);	
			delete dataset;
			delete dataspace;
		}
		//#################################################
		if(integrated_likelihoods_terms){
			hsize_t dimsILT[1];
			dimsILT[0]= ensemble_size;
			dataspace = new H5::DataSpace(1,dimsILT);
			dataset = new H5::DataSet(
				meta_group.createDataSet("INTEGRATED LIKELIHOODS TERM NUMBER",
					H5::PredType::NATIVE_INT,*dataspace)
				);
			dataset->write(integrated_likelihoods_terms, H5::PredType::NATIVE_INT);	
			delete dataset;
			delete dataspace;
		}
		//#################################################
		if(calculated_evidence){
			hsize_t dimsE[1];
			dimsE[0]= 1;
			dataspace = new H5::DataSpace(1,dimsE);
			dataset = new H5::DataSet(
				meta_group.createDataSet("EVIDENCE",
					H5::PredType::NATIVE_DOUBLE,*dataspace)
				);
			dataset->write(&evidence, H5::PredType::NATIVE_DOUBLE);	
			delete dataset;
			delete dataspace;
		}
		//#################################################
		dataspace = new H5::DataSpace(1,dimsT);
		dataset = new H5::DataSet(
			meta_group.createDataSet("SUGGESTED TRIM LENGTHS",
				H5::PredType::NATIVE_INT,*dataspace)
			);
		dataset->write(trim_lengths, H5::PredType::NATIVE_INT);	
		delete dataset;
		delete dataspace;
		//#################################################
		hsize_t dimsAC[2];
		dimsAC[0]= cold_chain_number;
		dimsAC[1]= dimension;
		dataspace = new H5::DataSpace(2,dimsAC);

		int *int_temp_buffer=NULL;
		if(ac_vals){
			int_temp_buffer = new int[cold_chain_number*dimension];
			for(int i  = 0 ; i<cold_chain_number; i++){
				for(int j = 0 ; j<dimension ; j++){
					int_temp_buffer[i*dimension +j ] = ac_vals[i][j];
				}
			}

			dataset = new H5::DataSet(
				meta_group.createDataSet("AC VALUES",
					H5::PredType::NATIVE_INT,*dataspace)
				);
			dataset->write(int_temp_buffer, H5::PredType::NATIVE_INT);	

			delete [] int_temp_buffer;
			int_temp_buffer =NULL;
			delete dataset;
			delete dataspace;
		}
		//#################################################
	
		//Cleanup
		output_LL_LP_group.close();
		output_group.close();
		if(RJ){
			status_group.close();
			model_status_group.close();
		}
		meta_group.close();
		if(!cold_only){
			delete [] ids;
			ids = NULL;
		}

	}	
	catch( H5::FileIException error )
	{
		error.printErrorStack();
		return -1;
	}
	// catch failure caused by the DataSet operations
	catch( H5::DataSetIException error )
	{
		error.printErrorStack();
	   	return -1;
	}
	// catch failure caused by the DataSpace operations
	catch( H5::DataSpaceIException error )
	{
		error.printErrorStack();
	   	return -1;
	}
	// catch failure caused by the DataSpace operations
	catch( H5::DataTypeIException error )
	{
		error.printErrorStack();
	   	return -1;
	}
	return 0;

}
int mcmc_sampler_output::append_to_data_dump( std::string filename)
{
	int file_id = 0;
	bool found=false;
	for(int i = 0 ; i<dump_file_names.size(); i++){
		if( filename == dump_file_names[i]){
			found = true;
			file_id = i;
		}
	}	
	if(!found){
		std::cout<<"ERROR -- File doesn't exist"<<std::endl;
	}
	try{
		std::string FILE_NAME(filename);
		int chains;
		int *ids=NULL;
		if(dump_files[file_id]->cold_only){
			ids = cold_chain_ids;
			chains = cold_chain_number;
		}
		else{
			chains = chain_number;
			ids = new int[chain_number];
			for(int i = 0  ; i<chain_number; i++){
				ids[i]=i;
			}
		}
		H5::H5File file(FILE_NAME,H5F_ACC_RDWR);
		H5::Group output_group(file.openGroup("/MCMC_OUTPUT"));
		H5::Group output_LL_LP_group(file.openGroup("/MCMC_OUTPUT/LOGL_LOGP"));
		H5::Group status_group;
		H5::Group model_status_group;
		if(RJ){
			status_group = H5::Group(file.openGroup("/MCMC_OUTPUT/STATUS"));
			model_status_group = H5::Group(file.openGroup("/MCMC_OUTPUT/MODEL_STATUS"));
		}
		H5::Group meta_group(file.openGroup("/MCMC_METADATA"));
		double *temp_buffer=NULL;
		double *temp_buffer_ll_lp=NULL;
		int *temp_buffer_status=NULL;
		int *temp_buffer_model_status=NULL;
		H5::DataSpace *dataspace=NULL ;
		H5::DataSpace *dataspace_ll_lp=NULL ;
		H5::DataSpace *dataspace_status=NULL ;
		H5::DataSpace *dataspace_model_status=NULL ;
		H5::DataSpace *dataspace_ext=NULL ;
		H5::DataSpace *dataspace_ext_ll_lp=NULL ;
		H5::DataSpace *dataspace_ext_status=NULL ;
		H5::DataSpace *dataspace_ext_model_status=NULL ;
		H5::DataSet *dataset=NULL;
		H5::DataSet *dataset_ll_lp=NULL;
		H5::DataSet *dataset_status=NULL;
		H5::DataSet *dataset_model_status=NULL;
		H5::DSetCreatPropList *plist=NULL;
		H5::DSetCreatPropList *plist_ll_lp=NULL;
		H5::DSetCreatPropList *plist_status=NULL;
		H5::DSetCreatPropList *plist_model_status=NULL;
		hsize_t chunk_dims[2] = {(hsize_t)chunk_steps,(hsize_t)dimension};
		hsize_t chunk_dims_ll_lp[2] = {(hsize_t)chunk_steps,2};
		hsize_t chunk_dims_status[2] = {(hsize_t)chunk_steps,(hsize_t)dimension};
		hsize_t chunk_dims_model_status[2] = {(hsize_t)chunk_steps,(hsize_t)nested_model_number};
		hsize_t max_dims[2] = {H5S_UNLIMITED,H5S_UNLIMITED};
		for(int i = 0 ; i<chains; i++){
			dataset = new H5::DataSet(output_group.openDataSet("CHAIN "+std::to_string(ids[i])));
			dataset_ll_lp = new H5::DataSet(output_LL_LP_group.openDataSet("CHAIN "+std::to_string(ids[i])));
			
			dataspace = new H5::DataSpace(dataset->getSpace());
			dataspace_ll_lp = new H5::DataSpace(dataset_ll_lp->getSpace());

			plist = new H5::DSetCreatPropList(dataset->getCreatePlist());
			plist_ll_lp = new H5::DSetCreatPropList(dataset_ll_lp->getCreatePlist());
			int RANK = dataspace->getSimpleExtentNdims();
			int RANK_ll_lp = dataspace_ll_lp->getSimpleExtentNdims();

			hsize_t base_dims[RANK];
			hsize_t base_dims_ll_lp[RANK_ll_lp];
			herr_t statusH5 = dataspace->getSimpleExtentDims(base_dims);
			statusH5 = dataspace_ll_lp->getSimpleExtentDims(base_dims_ll_lp);
			int RANK_chunked;
			int RANK_chunked_ll_lp;
			hsize_t base_chunk_dims[RANK];
			hsize_t base_chunk_dims_ll_lp[RANK_ll_lp];
			if(H5D_CHUNKED == plist->getLayout()){
				RANK_chunked= plist->getChunk(RANK,base_chunk_dims);
			}
			if(H5D_CHUNKED == plist_ll_lp->getLayout()){
				RANK_chunked_ll_lp= plist_ll_lp->getChunk(RANK_ll_lp,base_chunk_dims_ll_lp);
			}
			
			hsize_t new_size[RANK];
			hsize_t new_size_ll_lp[RANK];
			if(dump_files[file_id]->trimmed){
				new_size[0]= chain_lengths[ids[i]]-dump_files[file_id]->file_trim_lengths[ids[i]];
				new_size_ll_lp[0]= chain_lengths[ids[i]]-dump_files[file_id]->file_trim_lengths[ids[i]];
			}
			else{
				new_size[0]= chain_lengths[ids[i]];
				new_size_ll_lp[0]= chain_lengths[ids[i]];
			}
			new_size[1]= dimension;
			new_size_ll_lp[1]= 2;
			dataset->extend(new_size);
			dataset_ll_lp->extend(new_size_ll_lp);

			delete dataspace;
			delete dataspace_ll_lp;
			dataspace = new H5::DataSpace(dataset->getSpace());
			dataspace_ll_lp = new H5::DataSpace(dataset_ll_lp->getSpace());
			
			hsize_t dimext[RANK];	
			hsize_t dimext_ll_lp[RANK];	
			dimext[0]=new_size[0]-base_dims[0];
			dimext[1]=dimension;
			dimext_ll_lp[0]=new_size_ll_lp[0]-base_dims_ll_lp[0];
			dimext_ll_lp[1]=2;
			
			hsize_t offset[RANK];
			hsize_t offset_ll_lp[RANK];
			offset[0]=base_dims[0];	
			offset[1]=0;	
			offset_ll_lp[0]=base_dims_ll_lp[0];	
			offset_ll_lp[1]=0;	

			dataspace->selectHyperslab(H5S_SELECT_SET,dimext,offset);
			dataspace_ll_lp->selectHyperslab(H5S_SELECT_SET,dimext_ll_lp,offset_ll_lp);

			dataspace_ext = new H5::DataSpace(RANK, dimext,NULL);
			dataspace_ext_ll_lp = new H5::DataSpace(RANK_ll_lp, dimext_ll_lp,NULL);

			temp_buffer = new double[ dimext[0]*dimext[1] ];
			temp_buffer_ll_lp = new double[ dimext_ll_lp[0]*dimext_ll_lp[1] ];
			int beginning_id = 0 ; 
			if(dump_files[file_id]->trimmed){beginning_id = dump_files[file_id]->file_trim_lengths[ids[i]];}
			for(int j = base_dims[0] ; j<chain_lengths[ids[i]]-beginning_id; j++){
				for(int k = 0 ; k<dimension; k++){
					temp_buffer[(j-base_dims[0])*dimension +k] = output[ids[i]][j+beginning_id][k];	
				}
				temp_buffer_ll_lp[(j-base_dims_ll_lp[0])*2 ] = logL_logP[ids[i]][j+beginning_id][0];	
				temp_buffer_ll_lp[(j-base_dims_ll_lp[0])*2+1 ] = logL_logP[ids[i]][j+beginning_id][1];	
			}
			
			dataset->write(temp_buffer,H5::PredType::NATIVE_DOUBLE,*dataspace_ext, *dataspace);
			dataset_ll_lp->write(temp_buffer_ll_lp,H5::PredType::NATIVE_DOUBLE,*dataspace_ext_ll_lp, *dataspace_ll_lp);

			//Cleanup
			delete dataset;
			delete dataset_ll_lp;
			delete dataspace;
			delete dataspace_ll_lp;
			delete dataspace_ext;
			delete dataspace_ext_ll_lp;
			delete plist;
			delete plist_ll_lp;
			delete [] temp_buffer;
			delete [] temp_buffer_ll_lp;
			temp_buffer = NULL;
			temp_buffer_ll_lp = NULL;

			if(RJ){
				dataset_status = new H5::DataSet(status_group.openDataSet("CHAIN "+std::to_string(ids[i])));

				dataspace_status = new H5::DataSpace(dataset_status->getSpace());
				plist_status= new H5::DSetCreatPropList(dataset_status->getCreatePlist());
				int RANK_status = dataspace_status->getSimpleExtentNdims();
				hsize_t base_dims_status[RANK_status];
				herr_t statusH5 = dataspace_status->getSimpleExtentDims(base_dims_status);
				int RANK_chunked_status;
				hsize_t base_chunk_dims_status[RANK_status];
				if(H5D_CHUNKED == plist_status->getLayout()){
					RANK_chunked_status= plist_status->getChunk(RANK_status,base_chunk_dims_status);
				}
				
				hsize_t new_size_status[RANK];
				if(dump_files[file_id]->trimmed){
					new_size_status[0]= chain_lengths[ids[i]]-dump_files[file_id]->file_trim_lengths[ids[i]];
				}
				else{
					new_size_status[0]= chain_lengths[ids[i]];
				}
				new_size_status[1]= dimension;
				dataset_status->extend(new_size_status);

				delete dataspace_status;
				dataspace_status = new H5::DataSpace(dataset_status->getSpace());
				
				hsize_t dimext_status[RANK];	
				dimext_status[0]=new_size_status[0]-base_dims_status[0];
				dimext_status[1]=dimension;
				
				hsize_t offset_status[RANK];
				offset_status[0]=base_dims_status[0];	
				offset_status[1]=0;	

				dataspace_status->selectHyperslab(H5S_SELECT_SET,dimext_status,offset_status);

				dataspace_ext_status = new H5::DataSpace(RANK_status, dimext_status,NULL);

				temp_buffer_status = new int[ dimext_status[0]*dimext_status[1] ];
				int beginning_id = 0 ; 
				if(dump_files[file_id]->trimmed){beginning_id = dump_files[file_id]->file_trim_lengths[ids[i]];}
				for(int j = base_dims_status[0] ; j<chain_lengths[ids[i]]-beginning_id; j++){
					for(int k = 0 ; k<dimension; k++){
						temp_buffer_status[(j-base_dims_status[0])*dimension +k] = status[ids[i]][j+beginning_id][k];	
					}
				}
				
				dataset_status->write(temp_buffer_status,H5::PredType::NATIVE_INT,*dataspace_ext_status, *dataspace_status);
				//Cleanup
				delete dataset_status;
				delete dataspace_status;
				delete dataspace_ext_status;
				delete plist_status;
				delete [] temp_buffer_status;
				temp_buffer_status = NULL;

				//TODO -- this section can't be right.. 
				if(RJ ){
					dataset_model_status = new H5::DataSet(model_status_group.openDataSet("CHAIN "+std::to_string(ids[i])));

					dataspace_model_status = new H5::DataSpace(dataset_model_status->getSpace());
					plist_model_status= new H5::DSetCreatPropList(dataset_model_status->getCreatePlist());
					int RANK_model_status = dataspace_model_status->getSimpleExtentNdims();
					hsize_t base_dims_model_status[RANK_model_status];
					herr_t model_statusH5 = dataspace_model_status->getSimpleExtentDims(base_dims_model_status);
					int RANK_chunked_model_status;
					hsize_t base_chunk_dims_model_status[RANK_model_status];
					if(H5D_CHUNKED == plist_model_status->getLayout()){
						RANK_chunked_model_status= plist_model_status->getChunk(RANK_model_status,base_chunk_dims_model_status);
					}
					
					hsize_t new_size_model_status[RANK];
					if(dump_files[file_id]->trimmed){
						new_size_model_status[0]= chain_lengths[ids[i]]-dump_files[file_id]->file_trim_lengths[ids[i]];
					}
					else{
						new_size_model_status[0]= chain_lengths[ids[i]];
					}
					new_size_model_status[1]= nested_model_number;
					dataset_model_status->extend(new_size_model_status);

					delete dataspace_model_status;
					dataspace_model_status = new H5::DataSpace(dataset_model_status->getSpace());
					
					hsize_t dimext_model_status[RANK];	
					dimext_model_status[0]=new_size_model_status[0]-base_dims_model_status[0];
					dimext_model_status[1]=nested_model_number;
					
					hsize_t offset_model_status[RANK];
					offset_model_status[0]=base_dims_model_status[0];	
					offset_model_status[1]=0;	

					dataspace_model_status->selectHyperslab(H5S_SELECT_SET,dimext_model_status,offset_model_status);

					dataspace_ext_model_status = new H5::DataSpace(RANK_model_status, dimext_model_status,NULL);

					temp_buffer_model_status = new int[ dimext_model_status[0]*dimext_model_status[1] ];
					int beginning_id = 0 ; 
					if(dump_files[file_id]->trimmed){beginning_id = dump_files[file_id]->file_trim_lengths[ids[i]];}
					for(int j = base_dims_model_status[0] ; j<chain_lengths[ids[i]]-beginning_id; j++){
						temp_buffer_model_status[(j-base_dims_model_status[0])*nested_model_number ] = model_status[ids[i]][j+beginning_id];	
					}
					
					dataset_model_status->write(temp_buffer_model_status,H5::PredType::NATIVE_INT,*dataspace_ext_model_status, *dataspace_model_status);
					//Cleanup
					delete dataset_model_status;
					delete dataspace_model_status;
					delete dataspace_ext_model_status;
					delete plist_model_status;
					delete [] temp_buffer_model_status;
					temp_buffer_model_status = NULL;
				}
			}
		}



		//#####################################################
		dataset = new H5::DataSet(meta_group.openDataSet("CHAIN TEMPERATURES"));
		dataset->write(chain_temperatures, H5::PredType::NATIVE_DOUBLE);	
		delete dataset;
		//#####################################################
		if(integrated_likelihoods){
			dataset = new H5::DataSet(meta_group.openDataSet("INTEGRATED LIKELIHOODS"));
			dataset->write(integrated_likelihoods, H5::PredType::NATIVE_DOUBLE);	
			delete dataset;
		}
		//#####################################################
		if(integrated_likelihoods_terms){
			dataset = new H5::DataSet(meta_group.openDataSet("INTEGRATED LIKELIHOODS TERM NUMBER"));
			dataset->write(integrated_likelihoods_terms, H5::PredType::NATIVE_INT);	
			delete dataset;
		}
		//#####################################################
		if(calculated_evidence){
			dataset = new H5::DataSet(meta_group.openDataSet("EVIDENCE"));
			dataset->write(&evidence, H5::PredType::NATIVE_DOUBLE);	
			delete dataset;
		}
		//#####################################################
		if(!dump_files[file_id]->trimmed ){
			
			dataset = new H5::DataSet(meta_group.openDataSet("SUGGESTED TRIM LENGTHS"));
		
			dataset->write(trim_lengths, H5::PredType::NATIVE_INT);	
			delete dataset;
		}
		if(ac_vals){
			int *int_temp_buffer = new int[cold_chain_number*dimension];
			dataset = new H5::DataSet(meta_group.openDataSet("AC VALUES"));
			for(int i  = 0 ; i<cold_chain_number; i++){
				for(int j = 0 ; j<dimension ; j++){
					int_temp_buffer[i*dimension +j ] = ac_vals[i][j];
				}
			}
			dataset->write(int_temp_buffer, H5::PredType::NATIVE_INT);	
			delete [] int_temp_buffer;
			int_temp_buffer = NULL;
			delete dataset;
		}
	
		//Cleanup
		output_group.close();
		output_LL_LP_group.close();
		if(RJ){
			status_group.close();
			model_status_group.close();
		}
		meta_group.close();
		if(!dump_files[file_id]->cold_only){
			delete [] ids;
			ids = NULL;
		}

	}	
	catch( H5::FileIException error )
	{
		error.printErrorStack();
		return -1;
	}
	// catch failure caused by the DataSet operations
	catch( H5::DataSetIException error )
	{
		error.printErrorStack();
	   	return -1;
	}
	// catch failure caused by the DataSpace operations
	catch( H5::DataSpaceIException error )
	{
		error.printErrorStack();
	   	return -1;
	}
	// catch failure caused by the DataSpace operations
	catch( H5::DataTypeIException error )
	{
		error.printErrorStack();
	   	return -1;
	}
	return 0;

}
#else
int mcmc_sampler_output::create_data_dump(bool cold_only, bool trim,std::string filename)
{
	std::cout<<"ERROR -- only HDF5 is supported at the moment"<<std::endl;
	return 0;
}
int mcmc_sampler_output::append_to_data_dump( std::string filename)
{
	std::cout<<"ERROR -- only HDF5 is supported at the moment"<<std::endl;
	return 0;
}
int mcmc_sampler_output::write_flat_thin_output(std::string filename, bool use_stored_ac, bool trim)
{
	std::cout<<"ERROR -- only HDF5 is supported at the moment"<<std::endl;
	return 0;

}
#endif
void mcmc_sampler_output::append_integrated_likelihoods(double *integrated_likelihoods_new, int * integrated_likelihoods_terms_new, int ensemble_size_new)
{
	ensemble_size = ensemble_size_new;
	if(!integrated_likelihoods)
	{
		integrated_likelihoods = new double[ensemble_size];
		for(int i = 0 ; i<ensemble_size; i++){
			integrated_likelihoods[i] = 0 ;
		}
	}
	if(!integrated_likelihoods_terms)
	{	
		integrated_likelihoods_terms = new int[ensemble_size];
		for(int i = 0 ; i<ensemble_size; i++){
			integrated_likelihoods_terms[i] = 0 ;
		}
	}
	for(int i = 0; i<ensemble_size; i++){
		integrated_likelihoods[i] = integrated_likelihoods_terms[i]*integrated_likelihoods[i] + integrated_likelihoods_terms_new[i] * integrated_likelihoods_new[i];
		integrated_likelihoods_terms[i] += integrated_likelihoods_terms_new[i];
		integrated_likelihoods[i]/= integrated_likelihoods_terms[i];
	}
}
void mcmc_sampler_output::dealloc_integrated_likelihoods()
{
	calculated_evidence = false;
	evidence = 0;
	evidence_error = 0;
	if(integrated_likelihoods){
		delete [] integrated_likelihoods;
		integrated_likelihoods=NULL;
	}
	if(integrated_likelihoods_terms){
		delete [] integrated_likelihoods_terms;
		integrated_likelihoods_terms=NULL;
	}
}
void mcmc_sampler_output::calculate_evidence()
{
	//debugger_print(__FILE__,__LINE__,"Integrated likelihoods and number of terms");
	//for(int i = 0 ; i<ensemble_size; i++){
	//	std::cout<<integrated_likelihoods[i]<< " "<<integrated_likelihoods_terms[i]<<std::endl;
	//}
	calculated_evidence = true;
	int errcode = thermodynamic_integration(integrated_likelihoods, chain_temperatures, (int)(chain_number/ cold_chain_number), &evidence, &evidence_error);
	debugger_print(__FILE__,__LINE__,"Evidence: " + std::to_string(evidence));
	
}
//#############################################################
//#############################################################



