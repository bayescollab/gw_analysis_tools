#include "util.h"
#include "GWATConfig.h"
#include "D_Z_Config.h"
#include <math.h>
#include <string>
#include <string.h>
#include <complex>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <adolc/adouble.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


/*! \file
 *
 * General utilities that are not necessarily specific to any part of the project at large
 */



double gsl_LU_lndet(double **matrix, int dim)
{
	double lndet = 0;
	gsl_matrix *GSLmatrix = gsl_matrix_alloc(dim, dim);
	for(int row=0; row<dim; row++){
		for(int column = 0 ;column<dim; column++){
			gsl_matrix_set(GSLmatrix, row, column, matrix[row][column]);
		}
	}
	gsl_permutation *p = gsl_permutation_alloc(dim);
    	int s;
	gsl_linalg_LU_decomp(GSLmatrix, p, &s);
	lndet = gsl_linalg_LU_lndet(GSLmatrix);
	gsl_matrix_free(GSLmatrix);
	gsl_permutation_free(p);
	return lndet;
}

int mvn_sample(
	int samples, 
	double *mean,
	double **cov, 
	int dim, 
	double **output /**< [out] Size [samples][dim]*/
	)
{
	const gsl_rng_type *T;
	gsl_rng * r ;
	gsl_rng_env_setup();
	T  = gsl_rng_default;
	r= gsl_rng_alloc(T);	
	mvn_sample(samples, mean, cov, dim , r, output);	
	gsl_rng_free(r);
}

/*! \brief Samples from a multivariate-gaussian with covariance cov and mean mean, with dimension dim and puts the output in output
 * Taken from blogpost: https://juanitorduz.github.io/multivariate_normal/
 */
int mvn_sample(
	int samples, 
	double *mean,
	double **cov, 
	int dim, 
	gsl_rng *r,
	double **output /**< [out] Size [samples][dim]*/
	)
{
	gsl_matrix *matrix = gsl_matrix_alloc(dim,dim);
	for(int i = 0 ;i<dim ; i++){
		for (int j = 0 ; j<dim ; j++){
			gsl_matrix_set(matrix, i, j , cov[i][j]);
		}
	}
	int status = gsl_linalg_cholesky_decomp1(matrix);
	double **newmat = new double*[dim];
	for(int i = 0 ; i<dim ; i++){
		newmat[i]=new double[dim];
		for(int j = 0 ; j < dim; j++){
			if(j<=i){
				newmat[i][j] = gsl_matrix_get(matrix,i,j);
			}
			else{
				newmat[i][j] = 0;
			}
		}
	}
	//double test = gsl_matrix_get(matrix, row, column);
	
	double *random_nums = new double[samples*dim];
	for (int i = 0 ; i<samples*dim ; i++){
		random_nums[i] = gsl_ran_gaussian(r, 1);
	}
	
	for(int i = 0 ; i<samples; i++){
		for(int j =0  ; j< dim ; j ++){
			output[i][j] = mean[j];
			for(int k = 0 ; k<dim ; k++){
				output[i][j]+=newmat[j][k]*random_nums[i*dim +k ];
			}
		}
		
	}

	for(int i = 0 ; i<dim ; i++){
		delete [] newmat[i];
	}
	delete [] newmat;
	delete [] random_nums;

	gsl_matrix_free(matrix);
	
	return 0;
}


/*! \brief Just a quick utility to see if an item is in a list
 *
 * Didn't want to keep rewriting this loop
 */
template<>
bool check_list(int j, int *list, int length)
{
	for(int i = 0 ; i<length; i++){
		if(j == list[i]) return true;
	}
	return false;
}
template<>
bool check_list(std::string j, std::string *list, int length)
{
	for(int i = 0 ; i<length; i++){
		if(j == list[i]) return true;
	}
	return false;
}
template<>
bool check_list(double j, double *list, int length)
{
	for(int i = 0 ; i<length; i++){
		if(fabs(j- list[i])<DOUBLE_COMP_THRESH) return true;
	}
	return false;
}
//template bool check_list<int>(int, int*, int);
//template bool check_list<double>(double ,double*, int);
//template bool check_list<std::string>(std::string ,std::string*, int);


/*! \brief Matrix multiplication
 *
 * Takes A -- [dim1,dim2], B -- [dim2,dim3], and produces C -- [dim1,dim3]
 *
 * Arrays should be pre-allocated
 */
void matrix_multiply(double **A, double **B, double **C,int dim1, int dim2,int dim3)
{
	for(int i = 0 ; i<dim1; i++){
		for(int j = 0 ; j<dim3; j++){
			C[i][j] = 0;
			for(int k = 0; k<dim2; k++){
				C[i][j]+= A[i][k]*B[k][j];
			}
		}
	}
	return;
}

template<>
void debugger_print(const char *file, const int line,std::string message)
{
	std::cout.precision(15);
	std::cout<<"DEBUGGER:["<<strip_path(std::string(file))<<","<<line<<"]:";
	std::cout<<message<<std::endl;
}
template<>
void debugger_print(const char *file, const int line,const char *message)
{
	std::cout.precision(15);
	std::cout<<"DEBUGGER:["<<strip_path(std::string(file))<<","<<line<<"]:";
	std::cout<<message<<std::endl;
}
template<>
void debugger_print(const char *file, const int line,char *message)
{
	std::cout.precision(15);
	std::cout<<"DEBUGGER:["<<strip_path(std::string(file))<<","<<line<<"]:";
	std::cout<<message<<std::endl;
}
template<>
void debugger_print(const char *file, const int line,int message)
{
	std::cout.precision(15);
	std::cout<<"DEBUGGER:["<<strip_path(std::string(file))<<","<<line<<"]:";
	std::cout<<std::to_string(message)<<std::endl;
}
template<>
void debugger_print(const char *file, const int line,double message)
{
	std::cout.precision(15);
	std::cout<<"DEBUGGER:["<<strip_path(std::string(file))<<","<<line<<"]:";
	std::cout<<std::to_string(message)<<std::endl;
}
template<>
void debugger_print(const char *file, const int line,adouble message)
{
	std::cout.precision(15);
	std::cout<<"DEBUGGER:["<<strip_path(std::string(file))<<","<<line<<"]:";
	std::cout<<std::to_string(double(message.value()))<<std::endl;
}
std::string strip_path(std::string input)
{
	char sep = '/';
	size_t i = input.rfind(sep,input.length());
	if(i != std::string::npos){
		return input.substr(i+1, input.length() -1);
	}
	return "";
}

/*! \brief Newton Raphson method for a function of one dimension
 *
 * f_on_fprime is a function pointer that should assign the value result = f(x) / f'(x) 
 * 
 * param are just an extra needed parameters
 * 	
 * return values:
 *
 * 	0 -- success
 *
 * 	1 -- Bad value for derivative
 * 	
 * 	2 -- max iterations reached
 */
int newton_raphson_method_1d(void(*f_on_fprime)(double x, double *func, double *func_prime,void *param),double initial_guess,double tolerance ,int max_iterations,void *parameters, double *solution)
{
	double resid ,x1=initial_guess,x0=initial_guess,func, func_prime;
	int i = 0 ;
	do
	{
		x0=x1;
		f_on_fprime(x0,&func,&func_prime,parameters);
		if(fabs(func_prime)>1e-12  && 
			!isnan(func) && 
			!isnan(func_prime))
		{
			x1 = x0 - func/func_prime;
			i++;
		}
		else{
			//std::cout<<"UTIL "<<x0<<" "<<func_prime<<" "<<func<<std::endl;
			*solution= x1;
			return 1;
		}
	}while( fabs(x1 - x0) >tolerance && i<max_iterations);

	if(i==max_iterations){ *solution=x0; return 2;}
	*solution = x1;	
	return 0;
}

//#######################################################################################
//Interpolate Z to DL once per import
/*! \brief Function that uses the GSL libraries to interpolate pre-calculated Z-D_L data
 *
 * Initiates the requried functions -- GSL interpolation requires allocating memory before hand
 */
void initiate_LumD_Z_interp(gsl_interp_accel **Z_DL_accel_ptr, gsl_spline **Z_DL_spline_ptr)
{
	//int npts =100000;
	int npts =10000;
	double DLvec[npts];
	double Zvec[npts];
	*Z_DL_accel_ptr = gsl_interp_accel_alloc();
	*Z_DL_spline_ptr = gsl_spline_alloc(gsl_interp_cspline,npts);
	std::fstream data_table;
	//data_table.open(std::string(GWAT_ROOT_DIRECTORY)+"/data/Cosmology_data/tabulated_LumD_Z.csv",std::ios::in);
	data_table.open(std::string("GWAT_SHARE_DIR")+"Cosmology_data/tabulated_LumD_Z.csv",std::ios::in);
	std::vector<std::string> row;
	std::string line, word, temp;
	int i =0,j=0;
	if(data_table){
		while(std::getline(data_table,line)){
			std::stringstream lineStream(line);	
			std::string item;
			while(std::getline(lineStream, item, ','))
			{
				if(i<npts){
				if(j==0){DLvec[i]=std::stod(item);}
				else if(j==1){Zvec[i]=std::stod(item);}
				j++;
				}
			}
			j = 0;
			i ++;
		}
	}
	gsl_spline_init(*Z_DL_spline_ptr, DLvec, Zvec, npts);
	data_table.close();
}

/*! \brief Frees the allocated interpolation function
 */
void free_LumD_Z_interp(gsl_interp_accel **Z_DL_accel_ptr, gsl_spline **Z_DL_spline_ptr)
{
	gsl_interp_accel_free(*Z_DL_accel_ptr);
	gsl_spline_free(*Z_DL_spline_ptr);
}
//#######################################################################################

/*! Function that returns Z from a given luminosity Distance -- only Planck15
 *
 * adouble version for ADOL-C calculations
 */
adouble Z_from_DL_interp(adouble DL,gsl_interp_accel *Z_DL_accel_ptr, gsl_spline *Z_DL_spline_ptr)
{
	adouble Z = 0;
	Z = (adouble)gsl_spline_eval(Z_DL_spline_ptr, DL.value(), Z_DL_accel_ptr);
	return Z;
}

/*! Function that returns Z from a given luminosity Distance -- only Planck15
 */
double Z_from_DL_interp(double DL,gsl_interp_accel *Z_DL_accel_ptr, gsl_spline *Z_DL_spline_ptr)
{
	double Z = 0;
	double DLtemp = DL;
	//if(DL>DLvec[npts-1]){
	//	std::cout<<"WARNING: DL exceeded limit: setting to highest value in table"<<std::endl;
	//	DLtemp=DLvec[npts-1];
	//	std::cout<<DL<<std::endl;
	//	std::cout<<DLtemp<<std::endl;
	//}
	
	Z = gsl_spline_eval(Z_DL_spline_ptr, DLtemp, Z_DL_accel_ptr);
	return Z;

}
//#######################################################################################


/*! \brief Calculates the redshift given the luminosity distance
 *
 * Based on Astropy.cosmology calculations -- see python script in the ./data folder of the project -- numerically calculated given astropy.cosmology's definitions (http://docs.astropy.org/en/stable/cosmology/) and used scipy.optimize to fit to a power series, stepping in half powers of DL. These coefficients are then output to a header file (D_Z_config.h) which are used here to calculate redshift. Custom cosmologies etc can easily be acheived by editing the python script D_Z_config.py, the c++ functions do not need modification. They use whatever data is available in the header file.
 *
 * 5 cosmological models are available (this argument must be spelled exactly, although case insensitive):
 * 
 * PLANCK15, PLANCK13, WMAP9, WMAP7, WMAP5
 */
double Z_from_DL(double DL, std::string cosmology)
{
	std::string formatted_cosmo = "";
	std::locale loc;
  	for (std::string::size_type i=0; i<cosmology.length(); ++i)
    		formatted_cosmo+=std::toupper(cosmology[i], loc);
	int cosmo_index = cosmology_lookup(formatted_cosmo);
	if (cosmo_index == -1){ std::cout<<"Invalid Cosmology"<<std::endl;return -1;}
	const double *boundaries = boundaries_D[cosmo_index];
	int interp_deg = interp_degree[cosmo_index];
	//const double (**coeffs) = &COEFF_VEC_DZ[cosmo_index][0];
	int num_seg = num_segments[cosmo_index];
	double z;
	for (int i =0; i<num_seg; i++){
		if ( DL<boundaries[i+1]){
			double *coeffs = new double [interp_deg];
			for (int j =0; j<interp_deg;j++)
				coeffs[j]=COEFF_VEC_DZ[cosmo_index][i][j];
			z =  cosmology_interpolation_function(DL,coeffs, interp_deg);
			delete[] coeffs;
			return z;
			
		}	
	}
	return -1;
	
}
/*! \brief Calculates the redshift given the luminosity distance
 * adouble version for ADOL-C implementation
 */
adouble Z_from_DL(adouble DL, std::string cosmology)
{
	std::string formatted_cosmo = "";
	std::locale loc;
  	for (std::string::size_type i=0; i<cosmology.length(); ++i)
    		formatted_cosmo+=std::toupper(cosmology[i], loc);
	int cosmo_index = cosmology_lookup(formatted_cosmo);
	if (cosmo_index == -1){ std::cout<<"Invalid Cosmology"<<std::endl;return -1;}
	const double *boundaries = boundaries_D[cosmo_index];
	int interp_deg = interp_degree[cosmo_index];
	//const double (**coeffs) = &COEFF_VEC_DZ[cosmo_index][0];
	int num_seg = num_segments[cosmo_index];
	adouble z;
	for (int i =0; i<num_seg; i++){
		if ( DL<boundaries[i+1]){
			double *coeffs = new double [interp_deg];
			for (int j =0; j<interp_deg;j++)
				coeffs[j]=COEFF_VEC_DZ[cosmo_index][i][j];
			z =  cosmology_interpolation_function(DL,coeffs, interp_deg);
			delete[] coeffs;
			return z;
			
		}	
	}
	return -1;
	
}

/*! \brief Calculates the luminosity distance given the redshift
 *
 * Based on Astropy.cosmology calculations -- see python script in the ./data folder of the project -- numerically calculated given astropy.cosmology's definitions (http://docs.astropy.org/en/stable/cosmology/) and used scipy.optimize to fit to a power series, stepping in half powers of Z. These coefficients are then output to a header file (D_Z_config.h) which are used here to calculate distance. Custom cosmologies etc can easily be acheived by editing the python script D_Z_config.py, the c++ functions do not need modification. They use whatever data is available in the header file. If the functional form of the fitting function changes, these functions DO need to change.
 *
 * 5 cosmological models are available (this argument must be spelled exactly):
 * 
 * PLANCK15, PLANCK13, WMAP9, WMAP7, WMAP5
 */
double DL_from_Z(double Z, std::string cosmology)
{
	std::string formatted_cosmo = "";
	std::locale loc;
  	for (std::string::size_type i=0; i<cosmology.length(); ++i)
    		formatted_cosmo+=std::toupper(cosmology[i], loc);
	int cosmo_index = cosmology_lookup(formatted_cosmo);
	if (cosmo_index == -1){ std::cout<<"Invalid Cosmology"<<std::endl;return -1;}
	const double *boundaries = boundaries_Z[cosmo_index];
	int interp_deg = interp_degree[cosmo_index];
	//const double (**coeffs) = &COEFF_VEC_DZ[cosmo_index][0];
	int num_seg = num_segments[cosmo_index];
	double dl;
	for (int i =0; i<num_seg; i++){
		if ( Z<boundaries[i+1]){
			//if(i == 0 ){
			//	debugger_print(__FILE__,__LINE__,"Boundaries: "+std::to_string(boundaries[i])+" "+std::to_string(boundaries[i+1]) +" "+std::to_string(Z));
			//}
			double *coeffs = new double [interp_deg];
			for (int j =0; j<interp_deg;j++)
				coeffs[j]=COEFF_VEC_ZD[cosmo_index][i][j];
			dl =  cosmology_interpolation_function(Z,coeffs, interp_deg);
			delete[] coeffs;
			return dl;
			
		}	
	}
	return -1;
}
/*! \brief Calculates the luminosity distance given the redshift
 * adouble version for ADOL-C implementation
 */
adouble DL_from_Z(adouble Z, std::string cosmology)
{
	std::string formatted_cosmo = "";
	std::locale loc;
  	for (std::string::size_type i=0; i<cosmology.length(); ++i)
    		formatted_cosmo+=std::toupper(cosmology[i], loc);
	int cosmo_index = cosmology_lookup(formatted_cosmo);
	if (cosmo_index == -1){ std::cout<<"Invalid Cosmology"<<std::endl;return -1;}
	const double *boundaries = boundaries_Z[cosmo_index];
	int interp_deg = interp_degree[cosmo_index];
	//const double (**coeffs) = &COEFF_VEC_DZ[cosmo_index][0];
	int num_seg = num_segments[cosmo_index];
	adouble dl;
	for (int i =0; i<num_seg; i++){
		if ( Z<boundaries[i+1]){
			double *coeffs = new double [interp_deg];
			for (int j =0; j<interp_deg;j++)
				coeffs[j]=COEFF_VEC_ZD[cosmo_index][i][j];
			dl =  cosmology_interpolation_function(Z,coeffs, interp_deg);
			delete[] coeffs;
			return dl;
			
		}	
	}
	return -1;
}

/*! \brief Rough threshold for when to truncate waveforms for neutron stars 
 *
 * 24 km separation
 */
template <class T>
T fcontact(T M_detector, /**< Detector frame total mass in seconds*/
	T DL, /**< Luminsosity distance in seconds*/
	std::string cosmology/**< Cosmology to use for getting the redshift from dl*/
	)
{
	T z = Z_from_DL(DL/MPC_SEC,cosmology);	
	T Rcontact = 24./3.e5;
	return (1./M_PI)*(sqrt(M_detector/(1.+z)/pow_int(Rcontact,3)));
}
template double fcontact<double>(double , double , std::string);
template adouble fcontact<adouble>(adouble , adouble , std::string);
/*! \brief Custom interpolation function used in the cosmology calculations
 *
 * Power series in half power increments of x, up to 11/2. powers of x
 *
 */
double cosmology_interpolation_function(double x,double *coeffs, int interp_degree)
{
	double sum=coeffs[0];
	double rootx = std::sqrt(x);
		
	for(int i =1; i<interp_degree;i++){
		sum+= coeffs[i]*pow_int(rootx,i);
	}
	return sum;

}
/*! \brief Custom interpolation function used in the cosmology calculations
 * adouble version for ADOL-C
 */
adouble cosmology_interpolation_function(adouble x,double *coeffs, int interp_degree)
{
	adouble sum=coeffs[0];
	adouble rootx = sqrt(x);
		
	for(int i =1; i<interp_degree;i++){
		sum+= coeffs[i]*pow_int(rootx,i);
	}
	return sum;

}

/*! \brief Helper function for mapping cosmology name to an internal index
 */
int cosmology_lookup(std::string cosmology)
{
	for (int i =0; i<num_cosmologies; i++){
		if (cosmology == std::string(cosmos[i])){
			return i;
		}
	}
	return -1;
}
/*! \brief Converts a uniformly random random number into a powerlaw distribution
 *
 * https://mathworld.wolfram.com/RandomNumber.html
 */
double powerlaw_from_uniform(double x0, double x1, double power, double uniform_random_number)
{
	return pow( (pow(x1, power+1) - pow(x0,power+1) )*uniform_random_number + pow(x0, power+1) , 1./(power+1));
} 

/*! \brief Calculates a random number from the Maxwell-Boltzmann distribution using 3 gaussian random numbers
 *
 * Maxwell-Boltzmann describes a gaussian distributed vector in 3 space, so the magnitude is the RMS of the components of that vector
 *
 * Assumed mean is at 0
 *
 * Sigma would be equivalent to sqrt(kT/m) in the real maxwell boltzmann distribution
 *
 */
double gsl_maxwell_boltzmann_distribution(double sigma, gsl_rng *r)
{
	double v1=gsl_ran_gaussian(r, sigma);
	double v2=gsl_ran_gaussian(r, sigma);
	double v3=gsl_ran_gaussian(r, sigma);
	return sqrt(v1*v1 + v2*v2 + v3*v3);
} 

template <class T>
T copysign_internal(T val, T sign)
{
	return sqrt(val*val/(sign*sign))*sign;
}
template double copysign_internal<double>(double , double);
template adouble copysign_internal<adouble>(adouble , adouble);

/*! \brief Removes a dimension from a matrix (made with fishers in mind, but this is general)
 *
 * Note: the removed_dims list is indexed from 0
 */
void rm_fisher_dim(double **input,int full_dim, double **output,  int reduced_dim, int *removed_dims)
{
	int colct=0;
	int rowct=0;
	for(int j = 0 ;j <full_dim; j++){
		if(!check_list(j, removed_dims, full_dim-reduced_dim)){
			for(int k = 0 ; k<full_dim;  k++){
				if(!check_list(k,removed_dims, full_dim-reduced_dim)){
					output[rowct][colct] =input[j][k];
					colct++;
				}
			}
			rowct++;
			colct = 0;
		}
	}
}
/*! \brief Finds union of two lists which have been pre-sorted. Returns the union, also sorted
 *
 * Allocates memory for C, needs to be deallocated
 *
 */
void vector_union(std::vector<double> A, std::vector<double> B, std::vector<double> *C )
{
	auto itA = A.begin();
	auto itB = B.begin();
	while(itA !=A.end() && itB !=B.end())
	{
		if(*itA == *itB){
			(*C).push_back(*itA);
			itA++;
			itB++;
		}
		else if(*itA < *itB){
			(*C).push_back(*itA);
			itA++;

		}
		else {
			(*C).push_back(*itB);
			itB++;

		}
	}
	while(itA != A.end()){
		(*C).push_back(*itA);
		itA++;
	}
	while(itB != B.end()){
		(*C).push_back(*itB);
		itB++;
	}
}

/*! \brief SHOULDN'T BE USED FOR DOUBLESCustom list-interesction implementation for sorted lists 
 *
 * Uses pointers for efficiency -- NOT COPIED -- used std library for that
 *
 * Takes A and B and puts the intersection in C -- C = A /\ B
 *
 * Takes the length of C and puts it in lenC
 *
 * C should have allocated memory for at least as large as the smallest list
 * 	
 * WILL NOT work in-place, A, B, and C must be separate
 */
template<class T>
void list_intersect_ptrs(T **A, int lenA,T **B, int lenB, T **C, int *lenC)
{
	int max_length = (lenA>lenB) ? lenA : lenB;
	int i_A = 0;
	int i_B = 0;
	int i_C = 0;
	while(i_A<lenA && i_B < lenB){
		if(*(A[i_A]) == *(B[i_B])){
			//std::cout<<*(A[i_A])<<std::endl;;
			C[i_C] = (A[i_A]);
			i_C++;
			i_A++;
			i_B++;
		}
		else{
			if(*(A[i_A])>*(B[i_B])){
				i_B++;
			}
			else{
				i_A++;
			}
		}
	}
	*lenC = i_C;
}
//template void list_intersect_ptrs<double>(double**,int,double**,int,double**,int *);
//template void list_intersect_ptrs<adouble>(adouble**,int,adouble**,int,adouble**,int *);
template void list_intersect_ptrs<int>(int**,int,int**,int,int**,int *);
/*! \brief Custom list-interesction implementation for sorted lists 
 *
 * Uses pointers for efficiency -- NOT COPIED -- used std library for that
 *
 * Takes A and B and puts the intersection in C -- C = A /\ B
 *
 * Takes the length of C and puts it in lenC
 *
 * C should have allocated memory for at least as large as the smallest list
 * 	
 * WILL NOT work in-place, A, B, and C must be separate
 */
template<class T>
void list_intersect(T *A, int lenA,T *B, int lenB, T **C, int *lenC)
{
	int max_length = (lenA>lenB) ? lenA : lenB;
	int i_A = 0;
	int i_B = 0;
	int i_C = 0;
	while(i_A<lenA && i_B < lenB){
		if(A[i_A] == B[i_B]){
			//std::cout<<A[i_A]<<" "<<B[i_B]<<std::endl;
			C[i_C] = &(A[i_A]);
			i_C++;
			i_A++;
			i_B++;
		}
		else{
			if(A[i_A]>B[i_B]){
				//std::cout<<"B "<<B[i_B]<<std::endl;
				i_B++;
			}
			else{
				//std::cout<<"A "<<A[i_A]<<std::endl;
				i_A++;
			}
		}
	}
	*lenC = i_C;
}
template void list_intersect<double>(double*,int,double*,int,double**,int *);
//template void list_intersect<adouble>(adouble*,int,adouble*,int,adouble**,int *);
//template void list_intersect<int>(int*,int,int*,int,int**,int *);
/*! \brief Custom list-interesction implementation for sorted lists 
 *
 * Uses pointers for efficiency -- NOT COPIED -- used std library for that
 *
 * Takes A and B and puts the intersection in C -- C = A /\ B
 *
 * Takes the length of C and puts it in lenC
 *
 * C should have allocated memory for at least as large as the smallest list
 * 	
 * WILL NOT work in-place, A, B, and C must be separate
 */
template<class T>
void list_intersect_value(T *A, int lenA,T *B, int lenB, T *C, int *lenC)
{
	int max_length = (lenA>lenB) ? lenA : lenB;
	int i_A = 0;
	int i_B = 0;
	int i_C = 0;
	while(i_A<lenA && i_B < lenB){
		if(A[i_A] == B[i_B]){
			//std::cout<<A[i_A]<<" "<<B[i_B]<<std::endl;
			C[i_C] = (A[i_A]);
			i_C++;
			i_A++;
			i_B++;
		}
		else{
			if(A[i_A]>B[i_B]){
				//std::cout<<"B "<<B[i_B]<<std::endl;
				i_B++;
			}
			else{
				//std::cout<<"A "<<A[i_A]<<std::endl;
				i_A++;
			}
		}
	}
	*lenC = i_C;
}
//template void list_intersect_value<double>(double*,int,double*,int,double*,int *);
//template void list_intersect_value<adouble>(adouble*,int,adouble*,int,adouble*,int *);
template void list_intersect_value<int>(int*,int,int*,int,int*,int *);

template<class T, class U>
void mean_list(T *list, int length, U *result)
{
	U m = 0;
	for(int i = 0 ; i<length; i++){
		m +=list[i];
	}
	m/=length;
	*result=m;
}
template void mean_list<double,double>(double *,int,double*);
template void mean_list<adouble,adouble>(adouble *,int,adouble*);
template void mean_list<int,double>(int *,int,double*);
template void mean_list<int,adouble>(int *,int,adouble*);

template<class T, class U>
void variance_list(T *list, int length,U *result)
{
	U m;
	mean_list(list,length,&m);
	U variance=0;
	for(int i = 0; i<length ;i++){
		variance += pow_int(list[i] - m, 2);
	}
	variance/=length;
	*result=m;
}
template void variance_list<double,double>(double *,int,double*);
template void variance_list<adouble,adouble>(adouble *,int,adouble*);
template void variance_list<int,double>(int *,int,double*);
template void variance_list<int,adouble>(int *,int,adouble*);


/*! \brief Just a quick utility to see if an item is in a list
 *
 * Didn't want to keep rewriting this loop
 *
 * Now returns the id of the element
 *
 * Returns -1 if not found
 */
//template<class T>
//int check_list_id(T j, T *list, int length)
//{
//	for(int i = 0 ; i<length; i++){
//		if(j == list[i]) return i;
//	}
//	return -1;
//}
//template int check_list_id<int>(int, int*, int);
//template int check_list_id<double>(double ,double*, int);
template<>
int check_list_id(int j, int *list, int length)
{
	for(int i = 0 ; i<length; i++){
		if(j == list[i]) return i;
	}
	return -1;
}
int check_list_id(double j, double *list, int length)
{
	for(int i = 0 ; i<length; i++){
		if((j - list[i] ) < DOUBLE_COMP_THRESH) return i;
	}
	return -1;
}


template<class T>
void gsl_LU_matrix_invert(T **input, T **inverse, int dim)
{
	gsl_matrix *matrix = gsl_matrix_alloc(dim, dim);
	for(int row=0; row<dim; row++){
		for(int column = 0 ;column<dim; column++){
			gsl_matrix_set(matrix, row, column, input[row][column]);
		}
	}
	gsl_permutation *p = gsl_permutation_alloc(dim);
    	int s;
	gsl_linalg_LU_decomp(matrix, p, &s);
	gsl_matrix *inv = gsl_matrix_alloc(dim, dim);
	gsl_linalg_LU_invert(matrix, p, inv);
    
	gsl_permutation_free(p);
	
	
	for(int row=0; row<dim; row++){
		for(int column = 0 ;column<dim; column++){
			inverse[row][column] = gsl_matrix_get(inv,row,column);
		}
	}
	gsl_matrix_free(matrix);
	gsl_matrix_free(inv);
	
    
}
template void gsl_LU_matrix_invert<double>(double **, double**, int );

int gsl_cholesky_matrix_invert(double **input, double **inverse, int dim)
{
	gsl_matrix *matrix = gsl_matrix_alloc(dim, dim);
	for(int row=0; row<dim; row++){
		for(int column = 0 ;column<dim; column++){
			gsl_matrix_set(matrix, row, column, input[row][column]);
		}
	}
	
	gsl_permutation *p = gsl_permutation_alloc(dim);
	gsl_matrix *inv = gsl_matrix_alloc(dim, dim);
	int status =gsl_linalg_pcholesky_decomp(matrix,p);
	status = gsl_linalg_pcholesky_invert(matrix,p,inv);
	gsl_permutation_free(p);
	
	bool failed=false;
	if(status != 0 ){
		failed = true;
	}
	else{
		for(int row=0; row<dim; row++){
			for(int column = 0 ;column<dim; column++){
				if(!std::isnan(gsl_matrix_get(inv,row,column))){
					inverse[row][column] = gsl_matrix_get(inv,row,column);
				}
				else{
					//std::cout<<"Inversion failed -- NAN in inverse matrix -- cholesky decomposition"<<std::endl;
					failed=true;
					break;
				}
			}
			if(failed){break;}
		}
	}
	gsl_matrix_free(matrix);
	gsl_matrix_free(inv);
	if(failed){return 1;}
	return 0;
	
    
}

/*! \brief Normalize the Fisher matrix before inversion to try and tame singularity issues:
 *
 * arXiv:1108.1826v2 Equation 60-61
 */
int normalized_gsl_cholesky_matrix_invert(double **input, double **inverse, int dim)
{
	int status = 1;
	double ** A = new double*[dim];
	double ** out = new double*[dim];
	for (int i = 0 ; i<dim; i++){
		A[i] = new double[dim];
		out[i] = new double[dim];
		for(int j = 0 ; j<dim ; j++){
			if((input[i][i] * input[j][j]) <0){
				status=0;
			}
			A[i][j] = input[i][j] / sqrt(input[i][i] * input[j][j]);
		}
	}
	if (status ==0){ 
		std::cout<<"Error: Product of diagonal elements in Fisher normalization are negative"<<std::endl;
		return status;
	}
	status = gsl_cholesky_matrix_invert(A, out, dim);
	
	for (int i = 0 ; i<dim; i++){
		for(int j = 0 ; j< dim; j++){
			inverse[i][j] = out[i][j]/ sqrt(input[i][i] * input[j][j]);
		}
		delete [] A[i];
		delete [] out[i];
	}
	delete [] A;
	delete [] out;
	return status;
	
}

/*! \brief map the error on RA and DEC to error on solid angle \Omega
 *
 * All quantities in rad
 *
 * From arXiv:0906.4269
 */
double std_omega(double RA, double std_RA, double std_DEC, double cov_RA_DEC)
{
	return 2.*M_PI*( std::abs(sin(RA)) * sqrt( std_RA*std_DEC - cov_RA_DEC*cov_RA_DEC));
}

/*! \brief routine to print the progress of a process to the terminal as a progress bar
 *
 * Call everytime you want the progress printed
 */
void printProgress (double percentage)
{
    	int val = (int) (percentage * 100);
    	int lpad = (int) (percentage * PBWIDTH);
    	int rpad = PBWIDTH - lpad;
    	printf ("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    	fflush (stdout);
}

/*! \brief Allocate memory for FFTW3 methods used in a lot of inner products
 * input is a locally defined structure that houses all the pertinent data
 */
void allocate_FFTW_mem_forward(fftw_outline *plan, int length)
{
	plan->in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * length);	
	plan->out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * length);	
	plan->p = fftw_plan_dft_1d(length, plan->in, plan->out,FFTW_FORWARD, FFTW_MEASURE);
}
/*! \brief Allocate memory for FFTW3 methods used in a lot of inner products --INVERSE
 * input is a locally defined structure that houses all the pertinent data
 */
void allocate_FFTW_mem_reverse(fftw_outline *plan, int length)
{
	plan->in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * length);	
	plan->out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * length);	
	plan->p = fftw_plan_dft_1d(length, plan->in, plan->out,FFTW_BACKWARD, FFTW_MEASURE);
}
/*!\brief deallocates the memory used for FFTW routines
 */
void deallocate_FFTW_mem(fftw_outline *plan)
{
	fftw_destroy_plan(plan->p);
	fftw_free(plan->in);
	fftw_free(plan->out);
	//fftw_cleanup();
}

/*! \brief Builds the structure that shuttles source parameters between functions -updated version to incorporate structure argument
 *
 * Populates the structure that is passed to all generation methods - contains all relavent source parameters 
 *
 * Template type of source parameters and gen_parameters must match
 *
 */
template <class T>
void source_parameters<T>::populate_source_parameters(
			gen_params_base<T> *param_in
			) 
{

	/* Convert all dimensionful quantities to seconds and build all needed source quantities once*/
	this->mass1 = param_in->mass1*MSOL_SEC;
	this->mass2 = param_in->mass2*MSOL_SEC;
	this->spin1x = param_in->spin1[0];
	this->spin2x = param_in->spin2[0];
	this->spin1y = param_in->spin1[1];
	this->spin2y = param_in->spin2[1];
	this->spin1z = param_in->spin1[2];
	this->spin2z = param_in->spin2[2];
	this->chi_s = (1./2)*(this->spin1z+this->spin2z);
	this->chi_a = (1./2)*(this->spin1z-this->spin2z);
	//params.chirpmass = (adouble)calculate_chirpmass((double)params.mass1.value(),(double)params.mass2.value());
	this->chirpmass = calculate_chirpmass(this->mass1,this->mass2);
	//params.eta = (adouble)calculate_eta((double)params.mass1.value(),(double)params.mass2.value());	
	this->eta = calculate_eta(this->mass1,this->mass2);	
	this->M = this->mass1 + this->mass2;
	this->chi_eff = (this->mass1*(this->spin1z)+ this->mass2*(this->spin2z))/(this->M);
	this->chi_pn = this->chi_eff - (38*this->eta/113)*(2*this->chi_s);
	this->DL = param_in->Luminosity_Distance*MPC_SEC;
	this->delta_mass = sqrt(1.-4*this->eta);
	this->phiRef = param_in->phiRef;
	this->tc = param_in->tc;
	this->sky_average = param_in->sky_average;
	this->A0 = A0_from_DL(this->chirpmass, this->DL, this->sky_average);
	//return params;
	return;
}
//template <class T>
//source_parameters<T> source_parameters<T>::populate_source_parameters(
//			gen_params_base<T> *param_in
//			) 
//{
//
//	/* Convert all dimensionful quantities to seconds and build all needed source quantities once*/
//	source_parameters<T> params;
//	params.mass1 = param_in->mass1*MSOL_SEC;
//	params.mass2 = param_in->mass2*MSOL_SEC;
//	params.spin1x = param_in->spin1[0];
//	params.spin2x = param_in->spin2[0];
//	params.spin1y = param_in->spin1[1];
//	params.spin2y = param_in->spin2[1];
//	params.spin1z = param_in->spin1[2];
//	params.spin2z = param_in->spin2[2];
//	params.chi_s = (1./2)*(params.spin1z+params.spin2z);
//	params.chi_a = (1./2)*(params.spin1z-params.spin2z);
//	//params.chirpmass = (adouble)calculate_chirpmass((double)params.mass1.value(),(double)params.mass2.value());
//	params.chirpmass = calculate_chirpmass(params.mass1,params.mass2);
//	//params.eta = (adouble)calculate_eta((double)params.mass1.value(),(double)params.mass2.value());	
//	params.eta = calculate_eta(params.mass1,params.mass2);	
//	params.M = params.mass1 + params.mass2;
//	params.chi_eff = (params.mass1*(params.spin1z)+ params.mass2*(params.spin2z))/(params.M);
//	params.chi_pn = params.chi_eff - (38*params.eta/113)*(2*params.chi_s);
//	params.DL = param_in->Luminosity_Distance*MPC_SEC;
//	params.delta_mass = sqrt(1.-4*params.eta);
//	params.phiRef = param_in->phiRef;
//	params.tc = param_in->tc;
//	params.sky_average = param_in->sky_average;
//	params.A0 = A0_from_DL(params.chirpmass, params.DL, params.sky_average);
//	return params;
//}
/*! \brief TESTING Simple utility to copy the members of param_in to param_out, for whatever types those are.
 *
 * This is sometimes required to jump from double params to adouble params
 *
 * Memory MUST be allocated for betappe and bppe (if present), and this needs to be deallocated
 *
 */
template<class T, class U>
void transform_parameters(gen_params_base<T> *param_in, gen_params_base<U> **param_out)
{
	(*param_out)->mass1 = param_in->mass1;
	(*param_out)->mass2 = param_in->mass2;
	(*param_out)->Luminosity_Distance = param_in->Luminosity_Distance;
	(*param_out)->spin1[0] = param_in->spin1[0];
	(*param_out)->spin1[1] = param_in->spin1[1];
	(*param_out)->spin1[2] = param_in->spin1[2];
	(*param_out)->spin2[0] = param_in->spin2[0];
	(*param_out)->spin2[1] = param_in->spin2[1];
	(*param_out)->spin2[2] = param_in->spin2[2];
	(*param_out)->phiRef = param_in->phiRef;
	(*param_out)->tc = param_in->tc;
	(*param_out)->Nmod = param_in->Nmod;
	(*param_out)->Nmod_phi = param_in->Nmod_phi;
	(*param_out)->Nmod_sigma = param_in->Nmod_sigma;
	(*param_out)->Nmod_beta = param_in->Nmod_beta;
	(*param_out)->Nmod_alpha = param_in->Nmod_alpha;
	if(param_in->Nmod != 0){
		(*param_out)->bppe = new double[(*param_out)->Nmod];
		(*param_out)->betappe = new U[(*param_out)->Nmod];
		for(int i = 0 ;i<param_in->Nmod; i++){
			(*param_out)->bppe[i] = param_in->bppe[i];
			(*param_out)->betappe[i] = param_in->betappe[i];
		}
	}
	if(param_in->Nmod_phi != 0){
		(*param_out)->phii = new int[(*param_out)->Nmod_phi];
		(*param_out)->delta_phi = new U[(*param_out)->Nmod_phi];
		for(int i = 0 ;i<param_in->Nmod_phi; i++){
			(*param_out)->phii[i] = param_in->phii[i];
			(*param_out)->delta_phi[i] = param_in->delta_phi[i];
		}
	}
	if(param_in->Nmod_sigma != 0){
		(*param_out)->sigmai = new int[(*param_out)->Nmod_sigma];
		(*param_out)->delta_sigma = new U[(*param_out)->Nmod_sigma];
		for(int i = 0 ;i<param_in->Nmod_sigma; i++){
			(*param_out)->sigmai[i] = param_in->sigmai[i];
			(*param_out)->delta_sigma[i] = param_in->delta_sigma[i];
		}
	}
	if(param_in->Nmod_beta != 0){
		(*param_out)->betai = new int[(*param_out)->Nmod_beta];
		(*param_out)->delta_beta = new U[(*param_out)->Nmod_beta];
		for(int i = 0 ;i<param_in->Nmod_beta; i++){
			(*param_out)->betai[i] = param_in->betai[i];
			(*param_out)->delta_beta[i] = param_in->delta_beta[i];
		}
	}
	if(param_in->Nmod_alpha != 0){
		(*param_out)->alphai = new int[(*param_out)->Nmod_alpha];
		(*param_out)->delta_alpha = new U[(*param_out)->Nmod_alpha];
		for(int i = 0 ;i<param_in->Nmod_alpha; i++){
			(*param_out)->alphai[i] = param_in->alphai[i];
			(*param_out)->delta_alpha[i] = param_in->delta_alpha[i];
		}
	}
	(*param_out)->incl_angle = param_in->incl_angle;
	(*param_out)->theta = param_in->theta;
	(*param_out)->phi = param_in->phi;
	(*param_out)->RA = param_in->RA;
	(*param_out)->DEC = param_in->DEC;
	(*param_out)->gmst = param_in->gmst;
	(*param_out)->psi = param_in->psi;
	(*param_out)->NSflag1 = param_in->NSflag1;
	(*param_out)->NSflag2 = param_in->NSflag2;
	(*param_out)->f_ref = param_in->f_ref;
	(*param_out)->thetaJN = param_in->thetaJN;
	(*param_out)->alpha0 = param_in->alpha0;
	(*param_out)->chip = param_in->chip;
	(*param_out)->phip = param_in->phip;
	(*param_out)->chi1_l = param_in->chi1_l;
	(*param_out)->chi2_l = param_in->chi2_l;
	(*param_out)->phiJL = param_in->phiJL;
	(*param_out)->thetaJL = param_in->thetaJL;
	(*param_out)->zeta_polariz = param_in->zeta_polariz;
	(*param_out)->phi_aligned = param_in->phi_aligned;
	(*param_out)->chil = param_in->chil;
	(*param_out)->sky_average = param_in->sky_average;
	(*param_out)->shift_time = param_in->shift_time;
	(*param_out)->shift_phase = param_in->shift_phase;
	(*param_out)->equatorial_orientation = param_in->equatorial_orientation;
	(*param_out)->theta_l = param_in->theta_l;
	(*param_out)->phi_l = param_in->phi_l;
	(*param_out)->LISA_alpha0 = param_in->LISA_alpha0;
	(*param_out)->LISA_phi0 = param_in->LISA_phi0;
	(*param_out)->theta_j_ecl = param_in->theta_j_ecl;
	(*param_out)->phi_j_ecl = param_in->phi_j_ecl;
	(*param_out)->precess_reduced_flag = param_in->precess_reduced_flag;
	(*param_out)->dep_postmerger = param_in->dep_postmerger;
	
}
/*! \brief Simple utility to copy the members of param_in to param_out, for whatever types those are.
 *
 * This is sometimes required to jump from double params to adouble params
 *
 * Memory MUST be allocated for betappe and bppe (if present), and this needs to be deallocated
 *
 */
template<class T, class U>
void transform_parameters(gen_params_base<T> *param_in, gen_params_base<U> *param_out)
{
	param_out->mass1 = param_in->mass1;
	param_out->mass2 = param_in->mass2;
	param_out->Luminosity_Distance = param_in->Luminosity_Distance;
	param_out->spin1[0] = param_in->spin1[0];
	param_out->spin1[1] = param_in->spin1[1];
	param_out->spin1[2] = param_in->spin1[2];
	param_out->spin2[0] = param_in->spin2[0];
	param_out->spin2[1] = param_in->spin2[1];
	param_out->spin2[2] = param_in->spin2[2];
	param_out->phiRef = param_in->phiRef;
	param_out->tc = param_in->tc;
	param_out->Nmod = param_in->Nmod;
	param_out->Nmod_phi = param_in->Nmod_phi;
	param_out->Nmod_sigma = param_in->Nmod_sigma;
	param_out->Nmod_beta = param_in->Nmod_beta;
	param_out->Nmod_alpha = param_in->Nmod_alpha;
	if(param_in->Nmod != 0){
		param_out->bppe = new double[param_out->Nmod];
		param_out->betappe = new U[param_out->Nmod];
		for(int i = 0 ;i<param_in->Nmod; i++){
			param_out->bppe[i] = param_in->bppe[i];
			param_out->betappe[i] = param_in->betappe[i];
		}
	}
	if(param_in->Nmod_phi != 0){
		param_out->phii = new int[param_out->Nmod_phi];
		param_out->delta_phi = new U[param_out->Nmod_phi];
		for(int i = 0 ;i<param_in->Nmod_phi; i++){
			param_out->phii[i] = param_in->phii[i];
			param_out->delta_phi[i] = param_in->delta_phi[i];
		}
	}
	if(param_in->Nmod_sigma != 0){
		param_out->sigmai = new int[param_out->Nmod_sigma];
		param_out->delta_sigma = new U[param_out->Nmod_sigma];
		for(int i = 0 ;i<param_in->Nmod_sigma; i++){
			param_out->sigmai[i] = param_in->sigmai[i];
			param_out->delta_sigma[i] = param_in->delta_sigma[i];
		}
	}
	if(param_in->Nmod_beta != 0){
		param_out->betai = new int[param_out->Nmod_beta];
		param_out->delta_beta = new U[param_out->Nmod_beta];
		for(int i = 0 ;i<param_in->Nmod_beta; i++){
			param_out->betai[i] = param_in->betai[i];
			param_out->delta_beta[i] = param_in->delta_beta[i];
		}
	}
	if(param_in->Nmod_alpha != 0){
		param_out->alphai = new int[param_out->Nmod_alpha];
		param_out->delta_alpha = new U[param_out->Nmod_alpha];
		for(int i = 0 ;i<param_in->Nmod_alpha; i++){
			param_out->alphai[i] = param_in->alphai[i];
			param_out->delta_alpha[i] = param_in->delta_alpha[i];
		}
	}
	param_out->incl_angle = param_in->incl_angle;
	param_out->theta = param_in->theta;
	param_out->phi = param_in->phi;
	param_out->RA = param_in->RA;
	param_out->DEC = param_in->DEC;
	param_out->gmst = param_in->gmst;
	param_out->psi = param_in->psi;
	param_out->NSflag1 = param_in->NSflag1;
	param_out->NSflag2 = param_in->NSflag2;
	param_out->f_ref = param_in->f_ref;
	param_out->thetaJN = param_in->thetaJN;
	param_out->alpha0 = param_in->alpha0;
	param_out->chip = param_in->chip;
	param_out->phip = param_in->phip;
	param_out->chi1_l = param_in->chi1_l;
	param_out->chi2_l = param_in->chi2_l;
	param_out->phiJL = param_in->phiJL;
	param_out->thetaJL = param_in->thetaJL;
	param_out->zeta_polariz = param_in->zeta_polariz;
	param_out->phi_aligned = param_in->phi_aligned;
	param_out->chil = param_in->chil;
	param_out->sky_average = param_in->sky_average;
	param_out->shift_time = param_in->shift_time;
	param_out->shift_phase = param_in->shift_phase;
	param_out->equatorial_orientation = param_in->equatorial_orientation;
	param_out->theta_l = param_in->theta_l;
	param_out->phi_l = param_in->phi_l;
	param_out->LISA_alpha0 = param_in->LISA_alpha0;
	param_out->LISA_phi0 = param_in->LISA_phi0;
	param_out->theta_j_ecl = param_in->theta_j_ecl;
	param_out->phi_j_ecl = param_in->phi_j_ecl;
	param_out->precess_reduced_flag = param_in->precess_reduced_flag;
	param_out->dep_postmerger = param_in->dep_postmerger;
	
}
template void transform_parameters<double,adouble>(gen_params_base<double> *, gen_params_base<adouble> *);
template void transform_parameters<double,adouble>(gen_params_base<double> *, gen_params_base<adouble> **);
/*! \brief Transforms between chirpmass and DL to overall amplitude factor A0
 *
 * All quantities in seconds
 */
template<class T>
T A0_from_DL(T chirpmass, T DL, bool sky_average)
{
	if (sky_average){
		return sqrt(M_PI/30)*chirpmass*chirpmass/DL * pow(M_PI*chirpmass,-7./6);
	}
	else{
		return sqrt(M_PI*40./192.)*chirpmass*chirpmass/DL * pow(M_PI*chirpmass,-7./6);
	}
}
template double A0_from_DL<double>(double , double , bool);
template adouble A0_from_DL<adouble>(adouble , adouble , bool);
/*! \brief Transforms between amplitude factor A0 and chirpmass to DL 
 *
 * All quantities in seconds
 */
template<class T>
T DL_from_A0(T chirpmass, T A0, bool sky_average)
{
	if (sky_average){
		return sqrt(M_PI/30)*chirpmass*chirpmass/A0 * pow(M_PI*chirpmass,-7./6);
	}
	else{
		return sqrt(M_PI*40./192.)*chirpmass*chirpmass/A0 * pow(M_PI*chirpmass,-7./6);
	}
}
template double DL_from_A0<double>(double , double , bool);
template adouble DL_from_A0<adouble>(adouble , adouble , bool);
/*! \brief Builds the structure that shuttles source parameters between functions- outdated in favor of structure argument 
 *
 * Populates the structure that is passed to all generation methods - contains all relavent source parameters 
 */
template <class T>
source_parameters<T> source_parameters<T>::populate_source_parameters_old(
			T mass1, /**< mass of the larger body - in Solar Masses*/ 
			T mass2, /**< mass of the smaller body - in Solar Masses*/
			T Luminosity_Distance,/**< Luminosity Distance in Mpc*/ 
			T *spin1,/** spin vector of the larger body  {sx,sy,sz}*/
			T *spin2, /** spin vector of the smaller body  {sx,sy,sz}*/
			T phi_c,/** coalescence phase*/
			T t_c ,/** coalescence time*/
			bool sky_average
			) 
{

	/* Convert all dimensionful quantities to seconds and build all needed source quantities once*/
	source_parameters params;
	params.mass1 = mass1*MSOL_SEC;
	params.mass2 = mass2*MSOL_SEC;
	params.spin1x = spin1[0];
	params.spin2x = spin2[0];
	params.spin1y = spin1[1];
	params.spin2y = spin2[1];
	params.spin1z = spin1[2];
	params.spin2z = spin2[2];
	params.chi_s = (1./2)*(params.spin1z+params.spin2z);
	params.chi_a = (1./2)*(params.spin1z-params.spin2z);
	//params.chirpmass = (adouble)calculate_chirpmass((double)params.mass1.value(),(double)params.mass2.value());
	params.chirpmass = calculate_chirpmass(params.mass1,params.mass2);
	//params.eta = (adouble)calculate_eta((double)params.mass1.value(),(double)params.mass2.value());	
	params.eta = calculate_eta(params.mass1,params.mass2);	
	params.M = params.mass1 + params.mass2;
	params.chi_eff = (params.mass1*(params.spin1z)+ params.mass2*(params.spin2z))/(params.M);
	params.chi_pn = params.chi_eff - (38*params.eta/113)*(2*params.chi_s);
	params.DL = Luminosity_Distance*MPC_SEC;
	params.delta_mass = sqrt(1.-4*params.eta);
	params.phiRef = phi_c;
	params.tc = t_c;
	//params.A0 = sqrt(M_PI/30)*params.chirpmass*params.chirpmass/params.DL * pow(M_PI*params.chirpmass,-7./6);
	if(sky_average){
		params.A0 = sqrt(M_PI/30)*params.chirpmass*params.chirpmass/params.DL * pow(M_PI*params.chirpmass,-7./6);
	}
	else{
		params.A0 = sqrt(M_PI*40./192.)*params.chirpmass*params.chirpmass/params.DL * pow(M_PI*params.chirpmass,-7./6);
	}
	return params;
}


/*! \brief Routine to transform from the equatorial coordinate system spherical polar description of the total angular momentum to the detector specific polarization angle and inclination angle (of the total angular momentum)
 *
 * For the terrestrial network, this calculates the polarization angle for a detector at the center of Earth in equatorial coordinates, which is the polarization angle used in the detector_response_equatorial functions
 *
 * Polarization angle is defined as:
 *
 *  tan \psi = ( J.z - (J.N)(z.N) ) / (N.(Jxz)) 
 *
 *  Where z is the axis of rotation of earth (equatorial z), and N is the line of sight to the source (direction of propagation is -N)
 *
 *  The inclination angle is between the direction of propagation and the total angular momentum
 */
template<class T>
void terr_pol_iota_from_equat_sph(T RA, /**< Right ascension in rad*/
	T DEC, /**< Declination in rad*/
	T thetaj, /**< spherical polar angle of the total angular momentum in equatorial coordinates*/
	T phij, /**< spherical azimuthal angle of the total angular momentum in equatorial coordinates*/
	T *pol, /**< polarization defined by tan \psi = ( J.z - (J.N)(z.N) ) / (N.(Jxz)) */
	T *iota /**< Inclination angle of the TOTAL angular momentum and the direction of propagation -N*/
	)
{
	//std::cout<<"INTERNAL: RA: "<<RA<<" DEC: "<<DEC<<" thetaj: "<<thetaj<<" phij: "<<phij<<std::endl;
	//PSI only appears as 2*PSI in detector response, so atan should be fine. Periodicity of pi
	T temp_pol = atan(cos(DEC)*1./tan(thetaj)*1./sin(phij - RA) - 1./tan(phij - RA)*sin(DEC));
	*pol=M_PI/2. - temp_pol;
	//if(*pol<0){*pol+=M_PI;}//Fix range
	//*iota = acos(cos(thetaj)*sin(DEC) + cos(DEC)*cos(phij - RA)*sin(thetaj));
	//-Neq
	*iota = acos(-(cos(thetaj)*sin(DEC)) - cos(DEC)*cos(phij - RA)*sin(thetaj));
}
template void terr_pol_iota_from_equat_sph<double>(double, double, double, double, double *, double*);
template void terr_pol_iota_from_equat_sph<adouble>(adouble, adouble, adouble, adouble, adouble *, adouble*);


/*! \brief transform spherical angles from equatorial to ecliptic 
 *
 *
 * Rotation about the vernal equinox (x-hat in both coordinate systems) by the axial tilt or obliquity
 *
 * From wikipedia
 *
 * Compared against astropy and agrees
 */
template<class T>
void ecl_from_eq(T theta_eq, /**<Equatorial spherical polar angle */
	T phi_eq, /**< Equatorial spherical azimuthal angle*/
	T *theta_ecl, /**<Ecliptic spherical polar angle*/
	T *phi_ecl/**< Ecliptic spherical Azimuthal angle*/
	)
{
	T lambda; //Longitude
	T beta; //Latitude
	T RA = phi_eq;
	T DEC = M_PI/2. - theta_eq;
	T sra = sin(RA);
	T cra = cos(RA);
	T cdec = cos(DEC);
	T sdec = sin(DEC);
	T ce = cos(AXIAL_TILT);
	T se = sin(AXIAL_TILT);
	
	lambda = atan2(sra*ce + sdec/cdec * se, cra);
	beta = asin(sdec*ce - cdec* se*sra);
	*theta_ecl =M_PI/2. - beta ;
	*phi_ecl =lambda;
	if((*phi_ecl)<0) *phi_ecl+=2*M_PI;
}
template void ecl_from_eq<double>( double , double, double *, double*);
template void ecl_from_eq<adouble>( adouble , adouble, adouble *, adouble*);
/*! \brief Utility to map source frame vectors to equatorial frame vectors
 *
 * Needs the equatorial vectors for the line of sight and the spin angular momentum
 *
 * Needs the inclination angle and the reference phase
 * 
 * Works by constructing a third vector from the cross product of the two others, in both frames. This defines a family of 3 vectors 
 * that can uniquely determine a rotation matrix from source frame to equatorial, which was done analytically in mathematica ( see the nb anglenb.nb)
 *
 * K = LxN
 *
 * A = {LSF,NSF, KSF}
 * 
 * B={LEQ,NEQ,KEQ}
 *
 * B = R.A
 *
 * R = B.A^-1
 *
 * Verified with specific cases and with dot products with vectors before and after rotation
 */
template<class T>
void equatorial_from_SF(T *SFvec,/**< Input, source frame vector, as defined by LAL (L = z_hat) in cartesian*/
	T thetal, /**< Polar angle in equatorial coordinates of the orbital angular momentum at the reference frequency*/
	T phil, /**< Azimuthal angle in equatorial coordinates of the orbital angular momentum at the reference frequency*/
	T thetas,  /**< Polar angle in equatorial coordinates of the source (not direction of propagation)*/
	T phis,/**< Azimuthal angle in equatorial coordinates of the source (not direction of propagation)*/
	T iota, /**< Inclination angle between the orbiatl angular momentum and the direction of propagation at the reference frequency*/
	T phi_ref,/**<Reference frequency for the waveform (defines the LAL frame ) ( orbital phase) at the reference frequency*/
	T *EQvec/**< [out] Out vector in equatorial coordinates in cartesian*/)
{
	T cp = cos(M_PI/2.-phi_ref);
	T sp = sin(M_PI/2.-phi_ref);
	T ci = cos(iota);
	T si = sin(iota);
	T ctl = cos(thetal);
	T stl = sin(thetal);
	T cts = cos(thetas);
	T sts = sin(thetas);
	T cps = cos(phis);
	T sps = sin(phis);
	T cpl = cos(phil);
	T spl = sin(phil);
	T Jxn = SFvec[0];
	T Jyn = SFvec[1];
	T Jzn = SFvec[2];

	//T EQvec_test[3];
	//T EQvec_old[3];
	//EQvec_old[0] =  (pow_int(cp,2)*cpl*Jzn*si*stl - ci*cpl*(cp*Jxn + Jyn*sp)*stl + sp*(cpl*Jzn*si*sp*stl - cts*Jxn*spl*stl + cps*Jyn*sts + ctl*Jxn*sps*sts) + cp*(cts*Jyn*spl*stl + cps*Jxn*sts - ctl*Jyn*sps*sts))/(si*(pow_int(cp,2) + pow_int(sp,2)));
	//EQvec_old[1] = (pow_int(cp,2)*Jzn*si*spl*stl + cp*(-(cpl*cts*Jyn*stl) - ci*Jxn*spl*stl + cps*ctl*Jyn*sts + Jxn*sps*sts) + sp*(cpl*cts*Jxn*stl - ci*Jyn*spl*stl + Jzn*si*sp*spl*stl - cps*ctl*Jxn*sts + Jyn*sps*sts))/(si*(pow_int(cp,2) + pow_int(sp,2)));
	//EQvec_old[2] = (cp*cts*Jxn + pow_int(cp,2)*ctl*Jzn*si - ci*ctl*(cp*Jxn + Jyn*sp) + cp*Jyn*(-(cps*spl) + cpl*sps)*stl*sts + sp*(cts*Jyn + ctl*Jzn*si*sp + Jxn*(cps*spl - cpl*sps)*stl*sts))/(si*(pow_int(cp,2) + pow_int(sp,2)));
	//EQvec_test[0] = (pow_int(cp,2)*cpl*Jzn*si*stl - ci*cpl*(cp*Jxn + Jyn*sp)*stl + sp*(cpl*Jzn*si*sp*stl - cts*Jxn*spl*stl + cps*Jyn*sts + ctl*Jxn*sps*sts) +cp*(cts*Jyn*spl*stl + cps*Jxn*sts - ctl*Jyn*sps*sts))/si ;
	//EQvec_test[1] = (pow_int(cp,2)*Jzn*si*spl*stl + cp*(-(cpl*cts*Jyn*stl) - ci*Jxn*spl*stl + cps*ctl*Jyn*sts + Jxn*sps*sts) + sp*(cpl*cts*Jxn*stl - ci*Jyn*spl*stl + Jzn*si*sp*spl*stl - cps*ctl*Jxn*sts + Jyn*sps*sts))/si;
	//EQvec_test[2] = (cp*cts*Jxn + pow_int(cp,2)*ctl*Jzn*si - ci*ctl*(cp*Jxn + Jyn*sp) + cp*Jyn*(-(cps*spl) + cpl*sps)*stl*sts + sp*(cts*Jyn + ctl*Jzn*si*sp + Jxn*(cps*spl - cpl*sps)*stl*sts))/si;

	EQvec[0]=(pow_int(cp,2)*cpl*Jzn*si*stl - ci*cpl*(cp*Jxn + Jyn*sp)*stl - cp*(cts*Jyn*spl*stl + cps*Jxn*sts - ctl*Jyn*sps*sts) + 
     sp*(cpl*Jzn*si*sp*stl + cts*Jxn*spl*stl - (cps*Jyn + ctl*Jxn*sps)*sts))/si;
	EQvec[1] = (pow_int(cp,2)*Jzn*si*spl*stl + sp*(-(cpl*cts*Jxn*stl) - ci*Jyn*spl*stl + Jzn*si*sp*spl*stl + cps*ctl*Jxn*sts - Jyn*sps*sts) + 
     cp*(cpl*cts*Jyn*stl - ci*Jxn*spl*stl - (cps*ctl*Jyn + Jxn*sps)*sts))/si;
	EQvec[2] = (ctl*Jzn*si - Jxn*(ci*cp*ctl + cp*cts + sp*(cps*spl - cpl*sps)*stl*sts) - Jyn*(ci*ctl*sp + cts*sp + cp*(-(cps*spl) + cpl*sps)*stl*sts))/si ;
	//std::cout<<EQvec[0]-EQvec_test[0]<<std::endl;
	//std::cout<<EQvec[1]-EQvec_test[1]<<std::endl;
	//std::cout<<EQvec[2]-EQvec_test[2]<<std::endl;
	//std::cout<<EQvec[0]-EQvec_old[0]<<std::endl;
	//std::cout<<EQvec[1]-EQvec_old[1]<<std::endl;
	//std::cout<<EQvec[2]-EQvec_old[2]<<std::endl;
}
template void equatorial_from_SF<double>(double *, double, double, double, double, double, double,double *);
template void equatorial_from_SF<adouble>(adouble *, adouble, adouble, adouble, adouble, adouble, adouble,adouble *);
/*! \brief Calculates the chirp mass from the two component masses
 *
 * The output units are whatever units the input masses are
 */
double calculate_chirpmass(double mass1, double mass2)
{
	return pow(mass1 * mass2,3./5) / pow(mass1 + mass2,1./5);
}
adouble calculate_chirpmass(adouble mass1, adouble mass2)
{
	return pow(mass1 * mass2,3./5) / pow(mass1 + mass2,1./5);
}

/*!\brief Calculates the symmetric mass ration from the two component masses
 */
double calculate_eta(double mass1, double mass2)
{
	double eta = (mass1 * mass2) / pow(mass1 + mass2 ,2);
	// Check that eta is physical.
	// Numerical round-off sometimes causes it to not be,
	// and causes nan later.
	if (eta > 0.25)
	{
		eta = 0.25;
	}

	return eta;
}
adouble calculate_eta(adouble mass1, adouble mass2)
{
	return (mass1 * mass2) / pow(mass1 + mass2 ,2);
}

/*! \brief Calculates the larger mass given a chirp mass and symmetric mass ratio
 *
 * Units of the output match the units of the input chirp mass
 */
double calculate_mass1(double chirpmass, double eta)
{
	double etapow = pow(eta,3./5);
    	return 1./2*(chirpmass / etapow + sqrt(1.-4*eta)*chirpmass / etapow);
}

adouble calculate_mass1(adouble chirpmass, adouble eta)
{
	adouble etapow = pow(eta,3./5);
    	return 1./2*(chirpmass / etapow + sqrt(1.-4*eta)*chirpmass / etapow);
}
/*! \brief Calculates the smaller mass given a chirp mass and symmetric mass ratio
 *
 * Units of the output match the units of the input chirp mass
 */
double calculate_mass2(double chirpmass, double eta)
{
	double etapow = pow(eta,3./5);
    	return 1./2*(chirpmass / etapow - sqrt(1.-4*eta)*chirpmass / etapow);
}
adouble calculate_mass2(adouble chirpmass, adouble eta)
{
	adouble etapow = pow(eta,3./5);
    	return 1./2*(chirpmass / etapow - sqrt(1.-4*eta)*chirpmass / etapow);
}
/*! \brief Calculates the larger mass given a chirp mass and asymmetric mass ratio
 *
 * Units of the output match the units of the input chirp mass
 */
double calculate_mass1_Mcq(double chirpmass, double q)
{
	return 1./(1+q) * chirpmass * pow( q/pow_int(1+q,2), -3./5.);
}

adouble calculate_mass1_Mcq(adouble chirpmass, adouble q)
{
	return 1./(1+q) * chirpmass * pow( q/pow_int(1+q,2), -3./5.);
}
/*! \brief Calculates the smaller mass given a chirp mass and asymmetric mass ratio
 *
 * Units of the output match the units of the input chirp mass
 */
double calculate_mass2_Mcq(double chirpmass, double q)
{
	return q/(1+q) * chirpmass * pow( q/pow_int(1+q,2), -3./5.);
}
adouble calculate_mass2_Mcq(adouble chirpmass, adouble q)
{
	return q/(1+q) * chirpmass * pow( q/pow_int(1+q,2), -3./5.);
}

/*! \brief Local function to calculate a factorial
 */
long factorial(long num)
{
	int prod = 1;
	int step = num;
	while (step>0)
	{
		prod *= step;
		step -=1;
	}
	return prod;
}

/*! \brief Local power function, specifically for integer powers
 *
 * Much faster than the std version, because this is only for integer powers
 */
double pow_int(double base, int power)
{
	if (power == 0) return 1.;
	double prod = 1;
	int pow = std::abs(power);
	for (int i = 0; i< pow;i++){
		prod = prod * base;
	}
	if (power>0)
		return prod;
	else
		return 1./prod;
}
adouble pow_int(adouble base, int power)
{
	if (power == 0) return 1.;
	adouble prod = 1;
	int pow = std::abs(power);
	for (int i = 0; i< pow;i++)
		prod = prod * base;
	if (power>0)
		return prod;
	else
		return 1./prod;
}

/*! \brief Fucntion that just returns the cuberoot 
 */
double cbrt_internal(double base)
{
	return cbrt(base);
}
/*! \brief Fucntion that just returns the cuberoot 
 * ADOL-C doesn't have the cbrt function (which is faster),
 * so have to use the power function
 */
adouble cbrt_internal(adouble base)
{
	return pow(base,1./3.);
}

/*! \brief Utility to malloc 2D array
 * 
 */
double** allocate_2D_array( int dim1, int dim2)
{
	double **array = (double **) malloc(sizeof(double*)*dim1);
	for (int i = 0; i<dim1; i ++)
	{
		array[i] = (double*)malloc(sizeof(double ) * dim2);
	}
	return array;
}
int** allocate_2D_array_int( int dim1, int dim2)
{
	int **array = (int **) malloc(sizeof(int*)*dim1);
	for (int i = 0; i<dim1; i ++)
	{
		array[i] = (int*)malloc(sizeof(int ) * dim2);
	}
	return array;
}

/*! \brief Utility to free malloc'd 2D array
 * 
 */
void deallocate_2D_array(double **array, int dim1, int dim2)
{
	for(int i =0; i < dim1; i++)
	{
		free(array[i]);
	}
	free(array);
}
void deallocate_2D_array(int **array, int dim1, int dim2)
{
	for(int i =0; i < dim1; i++)
	{
		free(array[i]);
	}
	free(array);
}
/*! \brief Utility to malloc 3D array
 * 
 */
double*** allocate_3D_array( int dim1, int dim2, int dim3)
{
	double ***array = (double ***) malloc(sizeof(double**)*dim1);
	for (int i = 0; i<dim1; i ++)
	{
		array[i] = (double**)malloc(sizeof(double *) * dim2);
		for (int j = 0 ; j< dim2; j++)
		{
			array[i][j] = (double *)malloc(sizeof(double) * dim3);
		}
	}
	return array;
}
int*** allocate_3D_array_int( int dim1, int dim2, int dim3)
{
	int ***array = (int ***) malloc(sizeof(int**)*dim1);
	for (int i = 0; i<dim1; i ++)
	{
		array[i] = (int**)malloc(sizeof(int *) * dim2);
		for (int j = 0 ; j< dim2; j++)
		{
			array[i][j] = (int *)malloc(sizeof(int) * dim3);
		}
	}
	return array;
}
/*! \brief Utility to free malloc'd 2D array
 * 
 */
void deallocate_3D_array(double ***array, int dim1, int dim2, int dim3)
{
	for(int i =0; i < dim1; i++)
	{
		for(int j =0 ; j<dim2; j++){
			free(array[i][j]);
		}
		free(array[i]);
	}
	free(array);
}
/*! \brief Utility to free malloc'd 2D array
 * 
 */
void deallocate_3D_array(int ***array, int dim1, int dim2, int dim3)
{
	for(int i =0; i < dim1; i++)
	{
		for(int j =0 ; j<dim2; j++){
			free(array[i][j]);
		}
		free(array[i]);
	}
	free(array);
}


void hann_window(double *window,int length){
	for (int i = 0 ; i<length; i++){
		window[i] = 0.5 * ( 1 - std::cos(2 * M_PI * i / length));
	}
}

/*! \brief Tukey window function for FFTs
 *
 * As defined by https://en.wikipedia.org/wiki/Window_function
 */
void tukey_window(double *window,
		int length,
		double alpha)
{
	double imin = (int)(alpha*(double)(length-1)/2.0);
	double imax = (int)((double)(length-1)*(1.-alpha/2.));
	for (int i =0; i<length; i++){
		if(i<imin){
			window[i] = 0.5*(1 + cos(M_PI * ( ( (double) i)/((double)(imin)) -1) ) );
			//window[i] = 0.5*(1 - cos(M_PI * ( (2. * i)/(alpha * length) ) ) );
		}
		else if(i<imax){
			window[i] = 1;
		}
		else{
			window[i] = 0.5*(1 + cos(M_PI * ( ( (double) i)/((double)(imin)) -2./alpha + 1.) ) );
			//window[i] = 0.5*(1 + cos(M_PI * ( (2. * i)/(alpha * length) - 2./alpha + 1) ) );
			//window[i] = 0.5*(1 - cos(M_PI * ( (2. * i)/(alpha * length) ) ) );
		}
	}	

}

/*! \brief Utility to transform from celestial coord RA and DEC to local horizon coord for detector response functions
 *
 * Outputs are the spherical polar angles defined by North as 0 degrees azimuth and the normal to the earth as 0 degree polar
 */
template<class T>
void celestial_horizon_transform(T RA, /**< Right acsension (rad)*/
				T DEC, /**< Declination (rad)*/
				double gps_time, /**<GPS time */
				T LONG, /**< Longitude (rad)*/
				T LAT,/**< Latitude (rad)*/
				T *phi, /**<[out] horizon azimuthal angle (rad)*/
				T *theta/**< [out] horizon polar angle (rad)*/
				)
{
	//#################################
	//NEED TRANSFORM FROM GPS TO SIDEREAL
	T GMST = gps_to_GMST(gps_time);
	//###############################
	
	//std::cout<<"GMST: "<<GMST<<std::endl;
	T LMST = GMST + (LONG*180./M_PI)/15.; //Local mean sidereal in hours
	T H = (LMST - (RA*180./M_PI)/15.)*15.*M_PI/180.;//Local hour angle in rad
	
	T alt = asin( sin(DEC) * sin(LAT) + cos(DEC) * cos(LAT) *cos(H) );//alt in rad
	T a =  acos( (sin(DEC) - sin(alt)*sin(LAT) )/ (cos(alt)*cos(LAT))) ; //azimuth in rad
	T azimuth ;
	if (sin(H)<0) azimuth = a ;
	else azimuth = 2*M_PI - a;
	*phi = azimuth ;//output in rad
	*theta = M_PI/2. - alt;//output in rad
}
template void celestial_horizon_transform<double>(double,double,double,double,double,double *,double *);
template void celestial_horizon_transform<adouble>(adouble,adouble,double,adouble,adouble,adouble *,adouble *);

/*! \brief Utility to transform from gps time to GMST
 * https://aa.usno.navy.mil/faq/docs/GAST.php
 */
template<class  T>
T gps_to_GMST(T gps_time)
{
	T J2000 = 2451545;
	T JD = gps_to_JD(gps_time);
	T JD0;
	T H;
	//if((JD - floor(JD)) >.5){ 
	//	JD0 = floor(JD)+.5;//Julian date of the previous midnight
	//	H = (JD - JD0)*24;//Hours past midnight (in hours)
	//}
	//else{
	//	JD0 = floor(JD) -1. ;
	//	H = (JD - JD0)*24;
	//}
		JD0 = floor(JD)+.5;//Julian date of the previous midnight
		H = (JD - JD0)*24;//Hours past midnight (in hours)
	//else{
	//	JD0 = floor(JD) -1. ;
	//	H = (JD - JD0)*24;
	//}
	T D0 = JD0 -J2000;
	T D = JD -J2000;
	T tau = D/ 36525.; //Centuries since J2000
	//approximation of GMST from JD (from GPST)
	T GMST_unscaled = 6.697374558 + 0.06570982441908*D0 + 1.00273790935*H + 0.000026*tau*tau;
	T hours;
	//if(std::is_same<double,T>::value){
	hours = ((int)floor(GMST_unscaled)%24);
	//}
	//else{
		//hours = trunc(floor(GMST_unscaled)/24.);
		//hours = floor(GMST_unscaled)/24.-floor(floor(GMST_unscaled)/24.);
	//}
	//T hours = trunc(trunc(floor(GMST_unscaled))%24.);
	T fraction = GMST_unscaled - floor(GMST_unscaled);
	//return (6.697374558 + 0.06570982441908*D0 + 1.00273790935*H + 0.000026*T*T);
	return hours + fraction;
}
template double gps_to_GMST<double>(double);
//template adouble gps_to_GMST<adouble>(adouble);


template<class T>
T gps_to_GMST_radian(T gps_time)
{
	T hours = gps_to_GMST(gps_time);
	//TO RAD
	return hours*M_PI/12.;
}
template double gps_to_GMST_radian<double>(double);

/*! \brief Utility to transform from gps to JD
 */
template<class T>
T gps_to_JD(T gps_time)
{
	T J2000 = 2451545;
	T J2000_GPST = 630763213.;
	return J2000 + (gps_time-J2000_GPST)/(86400.);
}
template double gps_to_JD<double>(double);
template adouble gps_to_JD<adouble>(adouble);

/*! \brief Converts decimal into hours, minutes and seconds
 *
 * I just got tired of keeping all these steps straight..
 */
void decimal_to_HMS(double decimal, int *hour, int *min, double *second)
{
 	*hour = (int)floor(decimal);	
  	*min = (int)floor((decimal-*hour)*60);
 	*second = ((decimal-*hour)*60 - *min)*60;
}

/*! \brief utility to transform a vector from cartesian to spherical (radian)
 * 	
 * order:
 *
 * cart: x, y, z
 *
 * spherical: r, polar, azimuthal
 *
 * Tested in all octants
 */
template<class T>
void transform_cart_sph(T *cartvec, T *sphvec)
{
	sphvec[0]  = sqrt(cartvec[0]*cartvec[0] +
			cartvec[2]*cartvec[2] +
			cartvec[1] * cartvec[1]) ;

	if (sphvec[0] == 0.)
	{
		sphvec[1] = 0.;
		sphvec[2] = 0.;
	}
	else
	{
		sphvec[1] = acos(cartvec[2]/sphvec[0]);
		sphvec[2] = atan2(cartvec[1], cartvec[0]);
		if(sphvec[2]<0){sphvec[2]+=2*M_PI;}
	}
}
template void transform_cart_sph<double>(double*, double*);
template void transform_cart_sph<adouble>(adouble*, adouble*);
/*! \brief utility to transform a vector from spherical (radian) to cartesian
 *
 * Range: \theta \el [0,\pi]
 *
 * \phi \el [0,2 \pi]
 *
 * order:
 *
 * cart: x, y, z
 *
 * spherical: r, polar, azimuthal
 *
 * Tested in all octants
 */
template<class T>
void transform_sph_cart(T *sphvec, T *cartvec)
{
	cartvec[0] = sphvec[0] * sin(sphvec[1]) * cos(sphvec[2]);
	cartvec[1] = sphvec[0] * sin(sphvec[1]) * sin(sphvec[2]);
	cartvec[2] = sphvec[0] * cos(sphvec[1]) ;
}
template void transform_sph_cart<double>(double*, double*);
template void transform_sph_cart<adouble>(adouble*, adouble*);

/*! \brief Unwrap angles from arctan 
 *
 * Stolen from stack exchange.. https://stackoverflow.com/questions/15634400/continous-angles-in-c-eq-unwrap-function-in-matlab
 *
 * Modified
 */
template <class T>
void unwrap_array(T *in, T *out, int len) {
    out[0] = in[0];
    T d;
    for (int i = 1; i < len; i++) {
        d = in[i] - in[i-1];
        //d = d > M_PI ? d - 2 * M_PI : (d < -M_PI ? d + 2 * M_PI : d);
        if(d>M_PI){
		d -=2*M_PI;
	}
	else if(d<-M_PI){
		d+=2.*M_PI;
	}
        out[i] = out[i-1] + d;
    }
}
template void unwrap_array<double>(double * , double *, int);
template void unwrap_array<adouble>(adouble * , adouble *, int);
//################################################################
template <class T>
std::complex<T> cpolar(T mag, T phase)
{
	return mag * std::exp(std::complex<T>(0,1) * phase);
}

/*! Shamelessly stolen from LALsuite
 *
 */
template <class T>
std::complex<T> XLALSpinWeightedSphericalHarmonic(
                                   T theta,  /**< polar angle (rad) */
                                   T phi,    /**< azimuthal angle (rad) */
                                   int s,        /**< spin weight */
                                   int l,        /**< mode number l */
                                   int m         /**< mode number m */
    )
{
  T fac = 0.0;
  std::complex<T> ans = std::complex<T>(0.0,0.0);

  /* sanity checks ... */
  //if ( l < abs(s) ) 
  //{
  //  XLALPrintError("XLAL Error - %s: Invalid mode s=%d, l=%d, m=%d - require |s| <= l\n", __func__, s, l, m );
  //  XLAL_ERROR_VAL(0, XLAL_EINVAL);
  //}
  //if ( l < abs(m) ) 
  //{
  //  XLALPrintError("XLAL Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", __func__, s, l, m );
  //  XLAL_ERROR_VAL(0, XLAL_EINVAL);
  //}
  if ( s == -2 ) 
  {
    if ( l == 2 ) 
    {
      switch ( m ) 
      {
        case -2:
          fac = sqrt( 5.0 / ( 64.0 * M_PI ) ) * ( 1.0 - cos( theta ))*( 1.0 - cos( theta ));
          break;
        case -1:
          fac = sqrt( 5.0 / ( 16.0 * M_PI ) ) * sin( theta )*( 1.0 - cos( theta ));
          break;

        case 0:
          fac = sqrt( 15.0 / ( 32.0 * M_PI ) ) * sin( theta )*sin( theta );
          break;

        case 1:
          fac = sqrt( 5.0 / ( 16.0 * M_PI ) ) * sin( theta )*( 1.0 + cos( theta ));
          break;

        case 2:
          fac = sqrt( 5.0 / ( 64.0 * M_PI ) ) * ( 1.0 + cos( theta ))*( 1.0 + cos( theta ));
          break;
        //default:
        //  XLALPrintError("XLAL Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", __func__, s, l, m );
        //  XLAL_ERROR_VAL(0, XLAL_EINVAL);
        //  break;
      } /*  switch (m) */
    }  /* l==2*/
    else if ( l == 3 ) 
    {
      switch ( m ) 
      {
        case -3:
          fac = sqrt(21.0/(2.0*M_PI))*cos(theta/2.0)*pow(sin(theta/2.0),5.0);
          break;
        case -2:
          fac = sqrt(7.0/(4.0*M_PI))*(2.0 + 3.0*cos(theta))*pow(sin(theta/2.0),4.0);
          break;
        case -1:
          fac = sqrt(35.0/(2.0*M_PI))*(sin(theta) + 4.0*sin(2.0*theta) - 3.0*sin(3.0*theta))/32.0;
          break;
        case 0:
          fac = (sqrt(105.0/(2.0*M_PI))*cos(theta)*pow(sin(theta),2.0))/4.0;
          break;
        case 1:
          fac = -sqrt(35.0/(2.0*M_PI))*(sin(theta) - 4.0*sin(2.0*theta) - 3.0*sin(3.0*theta))/32.0;
          break;

        case 2:
          fac = sqrt(7.0/M_PI)*pow(cos(theta/2.0),4.0)*(-2.0 + 3.0*cos(theta))/2.0;
          break;

        case 3:
          fac = -sqrt(21.0/(2.0*M_PI))*pow(cos(theta/2.0),5.0)*sin(theta/2.0);
          break;

        //default:
        //  XLALPrintError("XLAL Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", __func__, s, l, m );
        //  XLAL_ERROR_VAL(0, XLAL_EINVAL);
        //  break;
      }
    }   /* l==3 */
    else if ( l == 4 ) 
    {
      switch ( m ) 
      {
        case -4:
          fac = 3.0*sqrt(7.0/M_PI)*pow(cos(theta/2.0),2.0)*pow(sin(theta/2.0),6.0);
          break;
        case -3:
          fac = 3.0*sqrt(7.0/(2.0*M_PI))*cos(theta/2.0)*(1.0 + 2.0*cos(theta))*pow(sin(theta/2.0),5.0);
          break;

        case -2:
          fac = (3.0*(9.0 + 14.0*cos(theta) + 7.0*cos(2.0*theta))*pow(sin(theta/2.0),4.0))/(4.0*sqrt(M_PI));
          break;
        case -1:
          fac = (3.0*(3.0*sin(theta) + 2.0*sin(2.0*theta) + 7.0*sin(3.0*theta) - 7.0*sin(4.0*theta)))/(32.0*sqrt(2.0*M_PI));
          break;
        case 0:
          fac = (3.0*sqrt(5.0/(2.0*M_PI))*(5.0 + 7.0*cos(2.0*theta))*pow(sin(theta),2.0))/16.0;
          break;
        case 1:
          fac = (3.0*(3.0*sin(theta) - 2.0*sin(2.0*theta) + 7.0*sin(3.0*theta) + 7.0*sin(4.0*theta)))/(32.0*sqrt(2.0*M_PI));
          break;
        case 2:
          fac = (3.0*pow(cos(theta/2.0),4.0)*(9.0 - 14.0*cos(theta) + 7.0*cos(2.0*theta)))/(4.0*sqrt(M_PI));
          break;
        case 3:
          fac = -3.0*sqrt(7.0/(2.0*M_PI))*pow(cos(theta/2.0),5.0)*(-1.0 + 2.0*cos(theta))*sin(theta/2.0);
          break;
        case 4:
          fac = 3.0*sqrt(7.0/M_PI)*pow(cos(theta/2.0),6.0)*pow(sin(theta/2.0),2.0);
          break;
        //default:
        //  XLALPrintError("XLAL Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", __func__, s, l, m );
        //  XLAL_ERROR_VAL(0, XLAL_EINVAL);
        //  break;
      }
    }    /* l==4 */
    else if ( l == 5 ) 
    {
      switch ( m ) 
      {
        case -5:
          fac = sqrt(330.0/M_PI)*pow(cos(theta/2.0),3.0)*pow(sin(theta/2.0),7.0);
          break;
        case -4:
          fac = sqrt(33.0/M_PI)*pow(cos(theta/2.0),2.0)*(2.0 + 5.0*cos(theta))*pow(sin(theta/2.0),6.0);
          break;
        case -3:
          fac = (sqrt(33.0/(2.0*M_PI))*cos(theta/2.0)*(17.0 + 24.0*cos(theta) + 15.0*cos(2.0*theta))*pow(sin(theta/2.0),5.0))/4.0;
          break;
        case -2:
          fac = (sqrt(11.0/M_PI)*(32.0 + 57.0*cos(theta) + 36.0*cos(2.0*theta) + 15.0*cos(3.0*theta))*pow(sin(theta/2.0),4.0))/8.0;
          break;
        case -1:
          fac = (sqrt(77.0/M_PI)*(2.0*sin(theta) + 8.0*sin(2.0*theta) + 3.0*sin(3.0*theta) + 12.0*sin(4.0*theta) - 15.0*sin(5.0*theta)))/256.0;
          break;
        case 0:
          fac = (sqrt(1155.0/(2.0*M_PI))*(5.0*cos(theta) + 3.0*cos(3.0*theta))*pow(sin(theta),2.0))/32.0;
          break;
        case 1:
          fac = sqrt(77.0/M_PI)*(-2.0*sin(theta) + 8.0*sin(2.0*theta) - 3.0*sin(3.0*theta) + 12.0*sin(4.0*theta) + 15.0*sin(5.0*theta))/256.0;
          break;
        case 2:
          fac = sqrt(11.0/M_PI)*pow(cos(theta/2.0),4.0)*(-32.0 + 57.0*cos(theta) - 36.0*cos(2.0*theta) + 15.0*cos(3.0*theta))/8.0;
          break;
        case 3:
          fac = -sqrt(33.0/(2.0*M_PI))*pow(cos(theta/2.0),5.0)*(17.0 - 24.0*cos(theta) + 15.0*cos(2.0*theta))*sin(theta/2.0)/4.0;
          break;
        case 4:
          fac = sqrt(33.0/M_PI)*pow(cos(theta/2.0),6.0)*(-2.0 + 5.0*cos(theta))*pow(sin(theta/2.0),2.0);
          break;
        case 5:
          fac = -sqrt(330.0/M_PI)*pow(cos(theta/2.0),7.0)*pow(sin(theta/2.0),3.0);
          break;
        //default:
        //  XLALPrintError("XLAL Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", __func__, s, l, m );
        //  XLAL_ERROR_VAL(0, XLAL_EINVAL);
        //  break;
      }
    }  /* l==5 */
    else if ( l == 6 )
    {
      switch ( m )
      {
        case -6:
          fac = (3.*sqrt(715./M_PI)*pow(cos(theta/2.0),4)*pow(sin(theta/2.0),8))/2.0;
          break;
        case -5:
          fac = (sqrt(2145./M_PI)*pow(cos(theta/2.0),3)*(1. + 3.*cos(theta))*pow(sin(theta/2.0),7))/2.0;
          break;
        case -4:
          fac = (sqrt(195./(2.0*M_PI))*pow(cos(theta/2.0),2)*(35. + 44.*cos(theta) 
          + 33.*cos(2.*theta))*pow(sin(theta/2.0),6))/8.0;
          break;
        case -3:
          fac = (3.*sqrt(13./M_PI)*cos(theta/2.0)*(98. + 185.*cos(theta) + 110.*cos(2*theta) 
          + 55.*cos(3.*theta))*pow(sin(theta/2.0),5))/32.0;
          break;
        case -2:
          fac = (sqrt(13./M_PI)*(1709. + 3096.*cos(theta) + 2340.*cos(2.*theta) + 1320.*cos(3.*theta) 
          + 495.*cos(4.*theta))*pow(sin(theta/2.0),4))/256.0;
          break;
        case -1:
          fac = (sqrt(65./(2.0*M_PI))*cos(theta/2.0)*(161. + 252.*cos(theta) + 252.*cos(2.*theta) 
          + 132.*cos(3.*theta) + 99.*cos(4.*theta))*pow(sin(theta/2.0),3))/64.0;
          break;
        case 0:
          fac = (sqrt(1365./M_PI)*(35. + 60.*cos(2.*theta) + 33.*cos(4.*theta))*pow(sin(theta),2))/512.0;
          break;
        case 1:
          fac = (sqrt(65./(2.0*M_PI))*pow(cos(theta/2.0),3)*(161. - 252.*cos(theta) + 252.*cos(2.*theta) 
          - 132.*cos(3.*theta) + 99.*cos(4.*theta))*sin(theta/2.0))/64.0;
          break;
        case 2:
          fac = (sqrt(13./M_PI)*pow(cos(theta/2.0),4)*(1709. - 3096.*cos(theta) + 2340.*cos(2.*theta) 
          - 1320*cos(3*theta) + 495*cos(4*theta)))/256.0;
          break;
        case 3:
          fac = (-3.*sqrt(13./M_PI)*pow(cos(theta/2.0),5)*(-98. + 185.*cos(theta) - 110.*cos(2*theta) 
          + 55.*cos(3.*theta))*sin(theta/2.0))/32.0;
          break;
        case 4:
          fac = (sqrt(195./(2.0*M_PI))*pow(cos(theta/2.0),6)*(35. - 44.*cos(theta) 
          + 33.*cos(2*theta))*pow(sin(theta/2.0),2))/8.0;
          break;
        case 5:
          fac = -(sqrt(2145./M_PI)*pow(cos(theta/2.0),7)*(-1. + 3.*cos(theta))*pow(sin(theta/2.0),3))/2.0;
          break;
        case 6:
          fac = (3.*sqrt(715./M_PI)*pow(cos(theta/2.0),8)*pow(sin(theta/2.0),4))/2.0;
          break;
        //default:
        //  XLALPrintError("XLAL Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", __func__, s, l, m );
        //  XLAL_ERROR_VAL(0, XLAL_EINVAL);
        //  break;
      }
    } /* l==6 */
    else if ( l == 7 )
    {
      switch ( m )
      {
        case -7:
          fac = sqrt(15015./(2.0*M_PI))*pow(cos(theta/2.0),5)*pow(sin(theta/2.0),9);
          break;
        case -6:
          fac = (sqrt(2145./M_PI)*pow(cos(theta/2.0),4)*(2. + 7.*cos(theta))*pow(sin(theta/2.0),8))/2.0;
          break;
        case -5:
          fac = (sqrt(165./(2.0*M_PI))*pow(cos(theta/2.0),3)*(93. + 104.*cos(theta) 
          + 91.*cos(2.*theta))*pow(sin(theta/2.0),7))/8.0;
          break;
        case -4:
          fac = (sqrt(165./(2.0*M_PI))*pow(cos(theta/2.0),2)*(140. + 285.*cos(theta) 
          + 156.*cos(2.*theta) + 91.*cos(3.*theta))*pow(sin(theta/2.0),6))/16.0;
          break;
        case -3:
          fac = (sqrt(15./(2.0*M_PI))*cos(theta/2.0)*(3115. + 5456.*cos(theta) + 4268.*cos(2.*theta) 
          + 2288.*cos(3.*theta) + 1001.*cos(4.*theta))*pow(sin(theta/2.0),5))/128.0;
          break;
        case -2:
          fac = (sqrt(15./M_PI)*(5220. + 9810.*cos(theta) + 7920.*cos(2.*theta) + 5445.*cos(3.*theta) 
          + 2860.*cos(4.*theta) + 1001.*cos(5.*theta))*pow(sin(theta/2.0),4))/512.0;
          break;
        case -1:
          fac = (3.*sqrt(5./(2.0*M_PI))*cos(theta/2.0)*(1890. + 4130.*cos(theta) + 3080.*cos(2.*theta) 
          + 2805.*cos(3.*theta) + 1430.*cos(4.*theta) + 1001.*cos(5*theta))*pow(sin(theta/2.0),3))/512.0;
          break;
        case 0:
          fac = (3.*sqrt(35./M_PI)*cos(theta)*(109. + 132.*cos(2.*theta) 
          + 143.*cos(4.*theta))*pow(sin(theta),2))/512.0;
          break;
        case 1:
          fac = (3.*sqrt(5./(2.0*M_PI))*pow(cos(theta/2.0),3)*(-1890. + 4130.*cos(theta) - 3080.*cos(2.*theta) 
          + 2805.*cos(3.*theta) - 1430.*cos(4.*theta) + 1001.*cos(5.*theta))*sin(theta/2.0))/512.0;
          break;
        case 2:
          fac = (sqrt(15./M_PI)*pow(cos(theta/2.0),4)*(-5220. + 9810.*cos(theta) - 7920.*cos(2.*theta) 
          + 5445.*cos(3.*theta) - 2860.*cos(4.*theta) + 1001.*cos(5.*theta)))/512.0;
          break;
        case 3:
          fac = -(sqrt(15./(2.0*M_PI))*pow(cos(theta/2.0),5)*(3115. - 5456.*cos(theta) + 4268.*cos(2.*theta) 
          - 2288.*cos(3.*theta) + 1001.*cos(4.*theta))*sin(theta/2.0))/128.0;
          break;  
        case 4:
          fac = (sqrt(165./(2.0*M_PI))*pow(cos(theta/2.0),6)*(-140. + 285.*cos(theta) - 156.*cos(2*theta) 
          + 91.*cos(3.*theta))*pow(sin(theta/2.0),2))/16.0;
          break;
        case 5:
          fac = -(sqrt(165./(2.0*M_PI))*pow(cos(theta/2.0),7)*(93. - 104.*cos(theta) 
          + 91.*cos(2.*theta))*pow(sin(theta/2.0),3))/8.0;
          break;
        case 6:
          fac = (sqrt(2145./M_PI)*pow(cos(theta/2.0),8)*(-2. + 7.*cos(theta))*pow(sin(theta/2.0),4))/2.0;
          break;
        case 7:
          fac = -(sqrt(15015./(2.0*M_PI))*pow(cos(theta/2.0),9)*pow(sin(theta/2.0),5));
          break;
        //default:
        //  XLALPrintError("XLAL Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", __func__, s, l, m );
        //  XLAL_ERROR_VAL(0, XLAL_EINVAL);
        //  break;
      }
    } /* l==7 */
    else if ( l == 8 )
    {
      switch ( m )
      {
        case -8:
          fac = sqrt(34034./M_PI)*pow(cos(theta/2.0),6)*pow(sin(theta/2.0),10);
          break;
        case -7:
          fac = sqrt(17017./(2.0*M_PI))*pow(cos(theta/2.0),5)*(1. + 4.*cos(theta))*pow(sin(theta/2.0),9);
          break;
        case -6:
          fac = sqrt(255255./M_PI)*pow(cos(theta/2.0),4)*(1. + 2.*cos(theta))
          *sin(M_PI/4.0 - theta/2.0)*sin(M_PI/4.0 + theta/2.0)*pow(sin(theta/2.0),8);
          break;
        case -5:
          fac = (sqrt(12155./(2.0*M_PI))*pow(cos(theta/2.0),3)*(19. + 42.*cos(theta) 
          + 21.*cos(2.*theta) + 14.*cos(3.*theta))*pow(sin(theta/2.0),7))/8.0;
          break;
        case -4:
          fac = (sqrt(935./(2.0*M_PI))*pow(cos(theta/2.0),2)*(265. + 442.*cos(theta) + 364.*cos(2.*theta) 
          + 182.*cos(3.*theta) + 91.*cos(4.*theta))*pow(sin(theta/2.0),6))/32.0;
          break;
        case -3:
          fac = (sqrt(561./(2.0*M_PI))*cos(theta/2.0)*(869. + 1660.*cos(theta) + 1300.*cos(2.*theta) 
          + 910.*cos(3.*theta) + 455.*cos(4.*theta) + 182.*cos(5.*theta))*pow(sin(theta/2.0),5))/128.0;
          break;
        case -2:
          fac = (sqrt(17./M_PI)*(7626. + 14454.*cos(theta) + 12375.*cos(2.*theta) + 9295.*cos(3.*theta) 
          + 6006.*cos(4.*theta) + 3003.*cos(5.*theta) + 1001.*cos(6.*theta))*pow(sin(theta/2.0),4))/512.0;
          break;
        case -1:
          fac = (sqrt(595./(2.0*M_PI))*cos(theta/2.0)*(798. + 1386.*cos(theta) + 1386.*cos(2.*theta) 
          + 1001.*cos(3.*theta) + 858.*cos(4.*theta) + 429.*cos(5.*theta) + 286.*cos(6.*theta))*pow(sin(theta/2.0),3))/512.0;
          break;
        case 0:
          fac = (3.*sqrt(595./M_PI)*(210. + 385.*cos(2.*theta) + 286.*cos(4.*theta) 
          + 143.*cos(6.*theta))*pow(sin(theta),2))/4096.0;
          break;
        case 1:
          fac = (sqrt(595./(2.0*M_PI))*pow(cos(theta/2.0),3)*(798. - 1386.*cos(theta) + 1386.*cos(2.*theta) 
          - 1001.*cos(3.*theta) + 858.*cos(4.*theta) - 429.*cos(5.*theta) + 286.*cos(6.*theta))*sin(theta/2.0))/512.0;
          break;
        case 2:
          fac = (sqrt(17./M_PI)*pow(cos(theta/2.0),4)*(7626. - 14454.*cos(theta) + 12375.*cos(2.*theta) 
          - 9295.*cos(3.*theta) + 6006.*cos(4.*theta) - 3003.*cos(5.*theta) + 1001.*cos(6.*theta)))/512.0;
          break;
        case 3:
          fac = -(sqrt(561./(2.0*M_PI))*pow(cos(theta/2.0),5)*(-869. + 1660.*cos(theta) - 1300.*cos(2.*theta) 
          + 910.*cos(3.*theta) - 455.*cos(4.*theta) + 182.*cos(5.*theta))*sin(theta/2.0))/128.0;
          break;
        case 4:
          fac = (sqrt(935./(2.0*M_PI))*pow(cos(theta/2.0),6)*(265. - 442.*cos(theta) + 364.*cos(2.*theta) 
          - 182.*cos(3.*theta) + 91.*cos(4.*theta))*pow(sin(theta/2.0),2))/32.0;
          break;
        case 5:
          fac = -(sqrt(12155./(2.0*M_PI))*pow(cos(theta/2.0),7)*(-19. + 42.*cos(theta) - 21.*cos(2.*theta) 
          + 14.*cos(3.*theta))*pow(sin(theta/2.0),3))/8.0;
          break;
        case 6:
          fac = sqrt(255255./M_PI)*pow(cos(theta/2.0),8)*(-1. + 2.*cos(theta))*sin(M_PI/4.0 - theta/2.0)
          *sin(M_PI/4.0 + theta/2.0)*pow(sin(theta/2.0),4);
          break;
        case 7:
          fac = -(sqrt(17017./(2.0*M_PI))*pow(cos(theta/2.0),9)*(-1. + 4.*cos(theta))*pow(sin(theta/2.0),5));
          break;
        case 8:
          fac = sqrt(34034./M_PI)*pow(cos(theta/2.0),10)*pow(sin(theta/2.0),6);
          break;
        //default:
        //  XLALPrintError("XLAL Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", __func__, s, l, m );
        //  XLAL_ERROR_VAL(0, XLAL_EINVAL);
        //  break;
      }
    } /* l==8 */
    //else 
    //{
    //  XLALPrintError("XLAL Error - %s: Unsupported mode l=%d (only l in [2,8] implemented)\n", __func__, l);
    //  XLAL_ERROR_VAL(0, XLAL_EINVAL);
    //}
  }
  //else 
  //{
  //  XLALPrintError("XLAL Error - %s: Unsupported mode s=%d (only s=-2 implemented)\n", __func__, s);
  //  XLAL_ERROR_VAL(0, XLAL_EINVAL);
  //}
  if (m)
    ans = cpolar((T)(1.0), (T)(m*phi)) * fac;
  else
    ans = fac;
  return ans;
}

template std::complex<double> XLALSpinWeightedSphericalHarmonic<double>(double,double,int,int,int);
template std::complex<adouble> XLALSpinWeightedSphericalHarmonic<adouble>(adouble,adouble,int,int,int);
template std::complex<double> cpolar<double>(double, double);
template std::complex<adouble> cpolar<adouble>(adouble,adouble);
template class source_parameters<double>;
template class source_parameters<adouble>;
