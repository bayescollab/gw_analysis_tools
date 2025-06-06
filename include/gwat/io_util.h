#ifndef INPUT_IO_H
#define INPUT_IO_H
#include <string>
#include <unordered_map>
#include <complex>

/*! \file
 *
 */
int unpack_input_io_file(std::string input_param_file, 
	std::unordered_map<std::string,int> *input_param_dict_int,
	std::unordered_map<std::string,std::string> *input_param_dict_str,
	std::unordered_map<std::string,double> *input_param_dict_dbl,
	std::unordered_map<std::string,float> *input_param_dict_flt,
	std::unordered_map<std::string,bool> *input_param_dict_bool
	);
int find_datatype(std::string type);
std::string trim(std::string str);
void read_file(std::string filename,double **output, int rows, int cols );
void read_file(std::string filename,int **output, int rows, int cols );
void read_file(std::string filename, double *output );
void read_file(std::string filename, int *output );
void write_file(std::string filename, double **input, int rows, int cols);
void write_file(std::string filename, int **input, int rows, int cols);
void write_file(std::string filename, double *input, int length);
void write_file(std::string filename, int *input, int length);
void read_LOSC_data_file(std::string filename, 
			double *output,
			double *data_start_time,
			double *duration,
			double *fs);

void read_LOSC_PSD_file(std::string filename, 
			double **output,
			int rows,
			int cols);

void allocate_LOSC_data(std::string *data_files, 
			std::string psd_file, 
			int num_detectors,
			int psd_length,
			int data_file_length,
			double trigger_time,
			double post_merger_duration,
			std::complex<double> **data,
			double **psds,
			double **freqs
			);

void free_LOSC_data(std::complex<double> **data,
		double **psds,
		double **freqs,
		int num_detectors,
		int length
		);
int count_lines_data_file(std::string file, int *count);
int count_lines_LOSC_PSD_file(std::string file, int *count);
int count_lines_LOSC_data_file(std::string file, int *count);



//###########################################################################
//HDF5
//###########################################################################
//class HDF5wrapper_ptr;
//class HDF5wrapper
//{
//	public:
// 		HDF5wrapper(int dims_outer,	
//			int *dims_inner,
//			int chunkdims_outer, 
//			int *chunkdims,
//			std::string DFILE_NAME,
//			std::string DSET_NAME,
//			bool deflate);
//		~HDF5wrapper();
//		
//	private:
//		HDF5wrapper_ptr *pimpl;
//};
#endif
