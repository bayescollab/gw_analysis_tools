#include "io_util.h"
#include "util.h"
#include <unordered_map>
#include <string>
#include <algorithm>
#include <cctype>
#include <complex>


/*! \file 
 *
 * unpacks standard parameter file for input for runtime parameters
 *
 * Returns the number of parameters unpacked, (-1 if error)
 * The format MUST be exact:
 * 	
 * 	[type] Parameter name = val
 *
 * 	[type] Parameter name = val
 *
 * The unordered_maps contain the parameter names and are mapped to the value val 
 *
 * hash symbol (#) reserved for comments
 *
 * Spelling must be exact
 *
 * Supported data_types = {int,dbl,str,flt,bool}
 *
 * If using special type, map from int to dtype separately
 *
 * Trims white space on ends
 *
 * Handles empty lines
 */
int unpack_input_io_file(std::string input_param_file, 
	std::unordered_map<std::string,int> *input_param_dict_int,
	std::unordered_map<std::string,std::string> *input_param_dict_str,
	std::unordered_map<std::string,double> *input_param_dict_dbl,
	std::unordered_map<std::string,float> *input_param_dict_flt,
	std::unordered_map<std::string,bool> *input_param_dict_bool
	)
{
	std::string line,param, dtype,value,temp;
	int dtype_assignment;
	std::fstream param_file(input_param_file, std::ios::in);	

	int ct=0;
	if(param_file){
		while(std::getline(param_file, line)){
			if(line==""){continue;}
			std::stringstream lineStream(line);
			std::getline(lineStream,temp,'[');
			if(temp.find("#") == std::string::npos){
				ct++;
				//std::cout<<line<<std::endl;
				std::getline(lineStream,dtype,']');
				std::getline(lineStream,param,'=');
				std::getline(lineStream,value);
		
				dtype = trim(dtype);
				param = trim(param);
				value = trim(value);
				dtype_assignment = find_datatype(dtype);
				if(dtype_assignment == 1){
					(*input_param_dict_int).insert({param ,stoi(value)});
				}
				else if(dtype_assignment == 2){
					(*input_param_dict_str).insert({param ,value});
				}
				else if(dtype_assignment == 3){
					(*input_param_dict_dbl).insert({param ,stod(value)});
				}
				else if(dtype_assignment == 4){
					(*input_param_dict_flt).insert({param ,stof(value)});
				}
				else if(dtype_assignment == 5){
					std::transform(value.begin(), value.end(), value.begin(),
    						[](unsigned char c){ return std::tolower(c); });
					if(value == "true"){
						(*input_param_dict_bool).insert({param ,true});
					}
					else{
						(*input_param_dict_bool).insert({param ,false});
					}
				}
				else{
					std::cout<<"ERROR -- invalid datatype"<<std::endl;
					param_file.close();
					return -1;
				}
			}
			else{
				continue;
			}
			
		}
	}
	else{
		std::cout<<"Input file not found"<<std::endl;
		param_file.close();
		return -1;
	}
	param_file.close();
	return ct;
}

/*! \brief Takes in the shorthand for primitive datatypes in c++ and returns code based on type
 *
 * Returns:
 *
 * 	0 -- Not found
 *
 * 	1 -- int
 *
 * 	2 -- string
 *
 * 	3 -- double
 *
 * 	4 -- float
 *
 * 	5 -- bool
 */
int find_datatype(std::string type)
{
	if(type == "int"){
		return 1;
	}
	else if(type == "str"){
		return 2;
	}
	else if(type == "dbl"){
		return 3;
	}
	else if(type == "flt"){
		return 4;
	}
	else if(type == "bool"){
		return 5;
	}
	return 0;

}

/* Blatantly stolen from stack exchange
 */
std::string trim(std::string str)
{	
	size_t first = str.find_first_not_of(' ');
	if (std::string::npos == first)
	{
	    return str;
	}
	size_t last = str.find_last_not_of(' ');
	return str.substr(first, (last - first + 1));
}
/*!\brief Utility to read in data
 *
 * Takes filename, and assigns to output[rows][cols]
 *
 * File must be comma separated doubles
 */
void read_file(std::string filename, /**< input filename, relative to execution directory*/
		double **output, /**<[out] array to store output, dimensions rowsXcols*/
		int rows, /**< first dimension*/
		int cols /**<second dimension*/
		)
{
	std::fstream file_in;
	file_in.open(filename, std::ios::in);
	std::string line, word;
	int i=0, j=0;
	double *temp = (double *)malloc(sizeof(double)*rows*cols);
	
	if(file_in){
		while(std::getline(file_in, line)){
			std::stringstream lineStream(line);
			std::string item;
			while(std::getline(lineStream,item, ',')){
				temp[i]=std::stod(item);	
				i+=1;	
			}	
		}	
	}
	else{std::cout<<"ERROR -- File "<<filename<<" not found"<<std::endl;exit(1);}
	for(i =0; i<rows;i++){
		for(j=0; j<cols;j++)
			output[i][j] = temp[cols*i + j];
	}
	free(temp);
}

/*!\brief Utility to read in data
 *
 * Takes filename, and assigns to output[rows][cols]
 *
 * File must be comma separated doubles
 *
 * integer version
 */
void read_file(std::string filename, /**< input filename, relative to execution directory*/
		int **output, /**<[out] array to store output, dimensions rowsXcols*/
		int rows, /**< first dimension*/
		int cols /**<second dimension*/
		)
{
	std::fstream file_in;
	file_in.open(filename, std::ios::in);
	std::string line, word;
	int i=0, j=0;
	double *temp = (double *)malloc(sizeof(double)*rows*cols);
	
	if(file_in){
		while(std::getline(file_in, line)){
			std::stringstream lineStream(line);
			std::string item;
			while(std::getline(lineStream,item, ',')){
				temp[i]=std::stoi(item);	
				i+=1;	
			}	
		}	
	}
	else{std::cout<<"ERROR -- File "<<filename<<" not found"<<std::endl;exit(1);}
	for(i =0; i<rows;i++){
		for(j=0; j<cols;j++)
			output[i][j] = temp[cols*i + j];
	}
	free(temp);
}

/*!\brief Utility to read in data
 *
 * Takes filename and delimiter of file, and assigns to ROW MAJOR 2D output vector
 *
 * File must be delimiter separated numerical entries, but can contain NaNs
 *
 * double version
 */
void read_file(std::string filename, /**< input filename, relative to execution directory*/
		std::vector<std::vector<double>>& output, /**<[out] vector to store output, dynamic dimensions to store full table*/
		char delimiter /**< input delimiter based on data file*/
	)
{
	std::fstream file_in;
	file_in.open(filename);

	if(file_in.good()) {	// Checks if the file was read in successfully with no errors
		std::string line;
		while(std::getline(file_in, line)){	// Runs until entire file is read
			std::vector<double> row;	// Row vector to store values
			std::stringstream lineStream(line);
            std::string token;
			while(std::getline(lineStream, token, delimiter))	// Reads items in line between delimiters
            {
                try {
                    double item = std::stod(token);	// Attemps to convert the read item to a double
                    row.push_back(item);
                }
                catch (...) {
                    row.push_back(std::nan(""));	// If double conversion fails, inputs NaN
                }
            }
            output.push_back(row); // Adds the row vector to the output table
        }
    }
    
	else {
		std::cout<<"ERROR -- File "<<filename<<" not found"<<std::endl;
		exit(1);
	}
}

/*!\brief Utility to read in data (single dimension vector) 
 *
 * Takes filename, and assigns to output[i*rows + cols]
 *
 * Output vector must be long enough, no check is done for the length
 *
 * File must be comma separated doubles
 */
void read_file(std::string filename, /**< input filename, relative to execution directory*/
	double *output /**<[out] output array, assumed to have the proper length of total items*/
	)
{
	std::fstream file_in;
	file_in.open(filename, std::ios::in);
	std::string line, word, temp;
	int i =0;
	if(file_in){
		while(std::getline(file_in, line)){
			std::stringstream lineStream(line);
			std::string item;
			while(std::getline(lineStream,item, ',')){
				output[i]=std::stod(item);	
				i+=1;
			}	
		}	
	}
	else{std::cout<<"ERROR -- File "<<filename<<" not found"<<std::endl;exit(1);}
}
/*!\brief Utility to read in data (single dimension vector) 
 *
 * Takes filename, and assigns to output[i*rows + cols]
 *
 * Output vector must be long enough, no check is done for the length
 *
 * File must be comma separated doubles
 *
 * Int version
 */
void read_file(std::string filename, /**< input filename, relative to execution directory*/
	int *output /**<[out] output array, assumed to have the proper length of total items*/
	)
{
	std::fstream file_in;
	file_in.open(filename, std::ios::in);
	std::string line, word, temp;
	int i =0;
	if(file_in){
		while(std::getline(file_in, line)){
			std::stringstream lineStream(line);
			std::string item;
			while(std::getline(lineStream,item, ',')){
				output[i]=std::stoi(item);	
				i+=1;
			}	
		}	
	}
	else{std::cout<<"ERROR -- File "<<filename<<" not found"<<std::endl;exit(1);}
}
/*! \brief Utility to write 2D array to file
 *
 * Grid of data, comma separated
 *
 * Grid has rows rows and cols columns
 */
void write_file(std::string filename, /**<Filename of output file, relative to execution directory*/
		double **input, /**< Input 2D array pointer array[rows][cols]*/
		int rows, /**< First dimension of array*/
		int cols /**< second dimension of array*/
		)
{
	
	std::ofstream out_file;
	out_file.open(filename);
	out_file.precision(15);
	if(out_file){
		for(int i =0; i<rows; i++){
			for(int j=0; j<cols;j++){
				if(j==cols-1)
					out_file<<input[i][j]<<std::endl;
				else
					out_file<<input[i][j]<<" , ";
			}
		}
		out_file.close();
	}
	else{
		std::cout<<"ERROR -- Could not open file"<<std::endl;
	}
}
/*! \brief Utility to write 2D array to file
 *
 * Grid of data, comma separated
 *
 * Grid has rows rows and cols columns
 *  
 * integer version
 */
void write_file(std::string filename, /**<Filename of output file, relative to execution directory*/
		int **input, /**< Input 2D array pointer array[rows][cols]*/
		int rows, /**< First dimension of array*/
		int cols /**< second dimension of array*/
		)
{
	
	std::ofstream out_file;
	out_file.open(filename);
	out_file.precision(15);
	if(out_file){
		for(int i =0; i<rows; i++){
			for(int j=0; j<cols;j++){
				if(j==cols-1)
					out_file<<input[i][j]<<std::endl;
				else
					out_file<<input[i][j]<<" , ";
			}
		}
		out_file.close();
	}
	else{
		std::cout<<"ERROR -- Could not open file"<<std::endl;
	}
}
/*! \brief Utility to write 1D array to file
 *
 * Single column of data
 */
void write_file(std::string filename, /**<Filename of output file, relative to execution directory*/
		double *input, /**< input 1D array pointer array[length]*/
		int length /**< length of array*/
		)
{
	std::ofstream out_file;
	out_file.open(filename);
	out_file.precision(15);
	if(out_file){
		for(int j =0; j<length; j++)
			out_file<<input[j]<<std::endl;
		out_file.close();
	}
	else{
		std::cout<<"ERROR -- Could not open file"<<std::endl;
	}
}
/*! \brief Utility to write 1D array to file
 *
 * Single column of data
 * 
 * integer version
 */
void write_file(std::string filename, /**<Filename of output file, relative to execution directory*/
		int *input, /**< input 1D array pointer array[length]*/
		int length /**< length of array*/
		)
{
	std::ofstream out_file;
	out_file.open(filename);
	out_file.precision(15);
	if(out_file){
		for(int j =0; j<length; j++)
			out_file<<input[j]<<std::endl;
		out_file.close();
	}
	else{
		std::cout<<"ERROR -- Could not open file"<<std::endl;
	}
}
/*! \brief Read data file from LIGO Open Science Center 
 *
 * Convenience function for cutting off the first few lines of text
 */
void read_LOSC_data_file(std::string filename, /**< input filename*/
			double *output,/**<[out] Output data*/
			double *data_start_time,/**<[out] GPS start time of the data in file*/
			double *duration,/**<[out] Duration of the signal*/
			double *fs/**<[out] Sampling frequency of the data*/
			)
{

	std::fstream file_in;
	file_in.open(filename, std::ios::in);
	std::string line, word, temp;
	int i =0;
	int j =0;
	if(file_in){
		while(std::getline(file_in, line)){
			std::stringstream lineStream(line);
			std::string item;

			//skip first three rows
			if (j>2){
				while(std::getline(lineStream,item, ',')){
					output[i]=std::stod(item);	
					i+=1;
				}	
			}
			//Extract data information from first 3 rows
			else{
				//Sampling frequency first
				if(j==1){
					std::istringstream iss(line);
					for(std::string s; iss>>s;){
						if(isdigit(s[0]))
							*fs = std::stod(s);
					}
				}
				else if (j==2){
					int k = 0;
					std::istringstream iss(line);
					for(std::string s; iss>>s;){
						if(isdigit(s[0])){
							//Time stamp
							if(k == 0){
								*data_start_time = std::stod(s);
								k++;
							}
							//Duration
							else{
								*duration = std::stod(s);
							}
						}
					}
				} 
				j+=1;
			}
		}	
	}
	else{std::cout<<"ERROR -- File "<<filename<<" not found"<<std::endl;exit(1);}
}

/*! \brief Read PSD file from LIGO Open Science Center 
 *
 * Convenience function for cutting off the first few lines of text
 */
void read_LOSC_PSD_file(std::string filename, 
			double **output,
			int rows,
			int cols
			)
{

	std::fstream file_in;
	file_in.open(filename, std::ios::in);
	std::string line;
	int i=0, j=0, k =0;
	std::vector<std::string> line_str;
	double *temp = (double *)malloc(sizeof(double)*rows*cols);
	
	if(file_in){
		while(std::getline(file_in, line)){
			std::istringstream iss(line);
			//skip the first row 
			if (k >0){	
				for(std::string s; iss>>s;){
					temp[i]=std::stod(s);	
					i+=1;	
				}
			}
			else{ k+=1;}
		}	
	}
	else{std::cout<<"ERROR -- File "<<filename<<" not found"<<std::endl;exit(1);}
	for(i =0; i<rows;i++){
		for(j=0; j<cols;j++)
			output[i][j] = temp[cols*i + j];
	}
	free(temp);

}
/*!\brief Prepare data for MCMC directly from LIGO Open Science Center
 *
 * Trims data for Tobs (determined by PSD file) 3/4*Tobs in front of trigger, and 1/4*Tobs behind or max 2 seconds on the end
 *
 * Currently, default to sampling frequency and observation time set by PSD -- cannot be customized
 *
 * Output is in order of PSD columns -- string vector of detectos MUST match order of PSD cols
 *
 * Output shapes-- 
 *
 * 		psds = [num_detectors][psd_length]
 *
 * 		data = [num_detectors][psd_length]	
 *
 * 		freqs = [num_detectors][psd_length]	
 *
 * Total observation time = 1/( freq[i] - freq[i-1]) (from PSD file)
 *
 * Sampling frequency fs = max frequency from PSD file
 *
 * ALLOCATES MEMORY -- must be freed to prevent memory leak
 */
void allocate_LOSC_data(std::string *data_files, /**< Vector of strings for each detector file from LOSC*/
			std::string psd_file, /**< String of psd file from LOSC*/
			int num_detectors,/**< Number of detectors to use*/
			int psd_length,/**< Length of the PSD file (number of rows of DATA)*/
			int data_file_length,/**< Length of the data file (number of rows of DATA)*/
			double trigger_time, /**< Time for the signal trigger (GPS)*/
			double post_merger_duration, /**< Time after trigger to end the data*/
			std::complex<double> **data,/**<[out] Output array of data for each detector*/
			double **psds,/**<[out] Output array of psds for each detector */
			double **freqs/**<[out] Output array of freqs for each detector*/
			)
{
	//Read in data from files
	double **temp_data = allocate_2D_array(num_detectors, data_file_length);
	double fs, duration, file_start;
	for(int i =0; i< num_detectors ; i++){
		read_LOSC_data_file(data_files[i],temp_data[i], &file_start, &duration, &fs);
	}	
	debugger_print(__FILE__,__LINE__,"Trigger time: "+std::to_string(trigger_time));
	debugger_print(__FILE__,__LINE__,"File initial GPS: "+std::to_string(file_start));
	debugger_print(__FILE__,__LINE__,"Data file duration: "+std::to_string(duration));
	debugger_print(__FILE__,__LINE__,"Data file sampling frequency: "+std::to_string(fs));

	//Read in frequencies and PSDs from files
	double **temp_psds = allocate_2D_array( psd_length,num_detectors+1);
	read_LOSC_PSD_file(psd_file,temp_psds,  psd_length,num_detectors+1);
	double Tobs = 1./(temp_psds[1][0] - temp_psds[0][0]);
	double df = 1./Tobs;
	double dt = 1./fs;
	int N_trimmed = Tobs*fs;
	for (int j = 0; j< psd_length; j++){
		for(int i =0; i< num_detectors ; i++){
			psds[i][j] = temp_psds[j][i+1];
			freqs[i][j] = temp_psds[j][0];
		}	
	}
		
	int N = fs*duration;
	double *times_untrimmed = (double *)malloc(sizeof(double)*N);


	for (int i =0; i < N; i++){
		times_untrimmed[i] = file_start + i*dt;
	}

	double time_start = trigger_time - (Tobs-post_merger_duration);
	double time_end = trigger_time + (post_merger_duration);
	//double time_start = trigger_time - Tobs*3./4.;
	//double time_end = trigger_time + Tobs/4.;
	
	//Shift for times longer than 8 -- no more than 2 seconds longer than merger
	//if(Tobs > 8 ){
	//	time_start = trigger_time - (Tobs - 2);
	//	time_end = trigger_time +  2;
	//}

	if(time_start < file_start || time_end > (file_start+duration) ){
		debugger_print(__FILE__,__LINE__,"ERROR -- trigger time does not fit inside data file!");
		return;
	}

	debugger_print(__FILE__,__LINE__,"Signal Duration: "+std::to_string(Tobs));
	debugger_print(__FILE__,__LINE__,"Start,End: "+std::to_string(time_start)+","+std::to_string(time_end));
	double **data_trimmed = allocate_2D_array(num_detectors, N_trimmed);
	double *times_trimmed = (double *)malloc(sizeof(double)*N_trimmed);
	double *window = (double *)malloc(sizeof(double)*N_trimmed);
	//double alpha = .4; //Standard alpha choice
	//double alpha = 2*0.4/Tobs; //Standard alpha choice
	double alpha = 2*.4/Tobs; //Standard alpha choice
	//double alpha = .1; //Standard alpha choice

	//double alpha = 0.5; //Standard alpha choice
	//double alpha = 0.5/32.; //Standard alpha choice
	tukey_window(window,N_trimmed, alpha);
	//for(int i = 0  ; i<N_trimmed; i++){
	//	window[i] = 1;
	//}
	//hann_window(window,N_trimmed);
	//Trim data to Tobs, and apply tukey windowing for fft
	int l=0 ;
	for (int i =0; i<N; i ++){
		if(times_untrimmed[i]>time_start && times_untrimmed[i]<=time_end){
			times_trimmed[l] = times_untrimmed[i];
			for (int j =0; j<num_detectors; j++){
				data_trimmed[j][l] = temp_data[j][i]*window[l];
			}	
			l++;
		}
	}
	fftw_outline plan;
	allocate_FFTW_mem_forward(&plan, N_trimmed);
	std::complex<double> **fft_data = (std::complex<double> **)
					malloc(sizeof(std::complex<double>*) * num_detectors);
	for (int i =0; i < num_detectors; i++){
		fft_data[i] = (std::complex<double>*)malloc(sizeof(std::complex<double>)*N_trimmed);
		for (int j =0; j<N_trimmed; j++){
			plan.in[j][0]=data_trimmed[i][j];
			plan.in[j][1]=0;//No imaginary part
		}
		fftw_execute(plan.p);
		for (int j =0; j<N_trimmed; j++){
			fft_data[i][j] = std::complex<double>(plan.out[j][0],plan.out[j][1]);
		}
		
	}	
	deallocate_FFTW_mem(&plan);
	double *freq_untrimmed = (double *)malloc(sizeof(double)*N_trimmed);
	for(int i =0; i<N_trimmed; i++){
		freq_untrimmed[i]=i*df;
	}
	double fmin = freqs[0][0];
	double fmax = freqs[0][psd_length-1];
	l = 0;
	//MODIFIED HERE
	for (int i =0 ; i<N_trimmed; i++){
		if(freq_untrimmed[i]>=fmin && freq_untrimmed[i]<=fmax){
			for(int j =0; j<num_detectors;j++){
	//#################################################
				//data[j][l] = fft_data[j][i];
				data[j][l] = fft_data[j][i]*dt;
	//#################################################
			}
			l++;
		}
	}

	//Deallocate temporary arrays
	free(times_trimmed);
	free(times_untrimmed);
	free(freq_untrimmed);
	free(window);
	deallocate_2D_array(temp_data,num_detectors, data_file_length);
	deallocate_2D_array(data_trimmed,num_detectors, N_trimmed);
	deallocate_2D_array(temp_psds, psd_length,num_detectors+1);
	for (int i =0; i<num_detectors; i++){
		free(fft_data[i]);
	}
	free(fft_data);
	
}

/*! /brief Free data allocated by prep_LOSC_data function
 */
void free_LOSC_data(std::complex<double> **data,
		double **psds,
		double **freqs,
		int num_detectors,
		int length
		)
{
	deallocate_2D_array(psds,num_detectors, length);
	deallocate_2D_array(freqs,num_detectors, length);
	for(int i =0; i<num_detectors; i++)
		free(data[i]);
	free(data);
}

/*! \brief Takes in data file and returns the number of data elements in file
 */
int count_lines_data_file(std::string file, int *count)
{
	std::ifstream file_input(file);
	std::string line;
	*count=0;
	if(file_input){
		while(std::getline(file_input,line)){
			*count+=1;
		}	
	}
	else{ 
		std::cout<<"ERROR -- could not open file"<<std::endl;
		return 1;	
	}
	
	return 0;
}
/*! \brief Takes in LOSC PSD file and returns the number of data elements in file
 */
int count_lines_LOSC_PSD_file(std::string file, int *count)
{
	count_lines_data_file(file,count);
	*count-=1;//Adjust for header
	return 0;
}
/*! \brief Takes in LOSC data file and returns the number of data elements in file
 */
int count_lines_LOSC_data_file(std::string file, int *count)
{
	count_lines_data_file(file,count);
	*count-=3;//Adjust for header
	return 0;
}


//###########################################################################
//HDF5
//###########################################################################
//#include <H5Cpp.h>

//class HDF5wrapper_ptr
//{
//public:
//	H5std_string DFILE_NAME="";
//	H5std_string DSET_NAME="";
//	H5::DataSet *DSET;
//	H5::H5File DFILE;
//	H5::DataSpace DSPACE;
//	int dims_outer;
//	int chunkdims_outer;
//	hsize_t *dims_inner=NULL;
//	H5::DSetCreatPropList *plist=NULL;
//	hsize_t *chunkdims=NULL;
//	bool deflate;
//	~HDF5wrapper_ptr()
//	{
//		if(!this->dims_inner){
//			delete [] this->dims_inner;
//		}
//		if(!this->chunkdims){
//			delete [] this->chunkdims;
//		}
//		if(!this->plist){
//			delete this->plist;
//		}
//		this->DFILE.close();
//
//	}
//	HDF5wrapper_ptr();
//	int allocate(int dims_outer,	
//		int *dims_inner,
//		int chunkdims_outer, 
//		int *chunkdims,
//		std::string DFILE_NAME,
//		std::string DSET_NAME,
//		bool deflate)
//	{
//		this->dims_outer = dims_outer;
//		this->dims_inner = new hsize_t[dims_outer];
//		this->chunkdims = new hsize_t[chunkdims_outer];
//		for(int i = 0 ; i<dims_outer; i++){
//			this->chunkdims[i]=chunkdims[i];
//		}
//		for(int i = 0 ; i<chunkdims_outer; i++){
//			this->dims_inner[i]=dims_inner[i];
//			this->chunkdims[i]=chunkdims[i];
//		}
//		this->deflate = deflate;
//		this->plist = new H5::DSetCreatPropList;
//		this->plist->setChunk((hsize_t)chunkdims_outer,(hsize_t *)chunkdims);
//		if(this->deflate){
//			this->plist->setDeflate(6);
//		}
//		this->DFILE_NAME = DFILE_NAME;
//		this->DFILE = H5::H5File(this->DFILE_NAME, H5F_ACC_TRUNC);
//		this->DSPACE = H5::DataSpace(dims_outer, this->dims_inner);
//		this->DSET = new H5::DataSet(this->DFILE.createDataSet(DSET_NAME, H5::PredType::NATIVE_DOUBLE,this->DSPACE, *(this->plist)));
//	}
//	int write_file(double *data)
//	{
//		this->DSET->write(data, H5::PredType::NATIVE_DOUBLE);
//		return 0;
//	}
//	int _index(int *location, int outerdim, int *dims)
//	{
//		int index = 0;
//		for(int i = 0 ; i<outerdim-1; i++){
//			index+= location[i]*dims[i+1];
//		}
//		index+=location[outerdim-1];
//		return index;
//	}
//	
//};
//
// HDF5wrapper::HDF5wrapper(int dims_outer,	
//			int *dims_inner,
//			int chunkdims_outer, 
//			int *chunkdims,
//			std::string DFILE_NAME,
//			std::string DSET_NAME,
//			bool deflate)
//{
//			this->pimpl=new HDF5wrapper_ptr();
//			//this->pimpl->allocate(dims_outer, dims_inner, chunkdims, DFILE_NAME,DSET_NAME,deflate);	
//}
//HDF5wrapper::~HDF5wrapper()
//{
//	this->pimpl->~HDF5wrapper_ptr();
//	delete this->pimpl;
//}


