#ifndef UTIL_H
#define UTIL_H
//#include "general_parameter_structures.h"
#include <string>
#include <complex>
#include "adolc/adouble.h"
#include <fftw3.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
//#define adouble double
/*! \file
 *General utilities (functions and structures) independent of modelling method
 */

/*! Euler number*/
const double gamma_E = 0.5772156649015328606065120900824024310421;
/*!Speed of light m/s*/
const double c = 299792458.;
/*!Gravitational constant in m**3/(s**2 SolMass)*/
//const double G =6.67430e-11*(1.98855e30);
const double G =6.67430e-11*(1.988409902147041637325262574352366540e30);
/*! G/c**3 seconds per solar mass*/
//const double MSOL_SEC =492549095.e-14;
//const double MSOL_SEC =492549095.e-14;
const double MSOL_SEC =4.925491025543575903411922162094833998e-6 ;
/*!consts.kpc.to('m')*1000/c Mpc in sec*/
//const double MPC_SEC = 3085677581.e13/c;
const double MPC_SEC = 3.085677581491367278913937957796471611e22/c;
/*!1 year in seconds -- ie seconds/year */
const double T_year = 31557600.;
/*!1 day in seconds*/
const double T_day = 24*3600.;
/*! Earth's Axial tilt in radian*/
const double AXIAL_TILT=0.409092627749;
//const double AXIAL_TILT=0.4075;
/*! 1AU in seconds*/
const double AU_SEC = 499.005;

/*! Planck's constant (unbarred) in ev sec*/
const double h_planck = 4.135667696e-15;

/*! Sqrt(3) comes up often, precompute it*/
const double ROOT_THREE = std::sqrt(3.);

const double LOG10=std::log(10.);

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

const double DOUBLE_COMP_THRESH = 1e-10;
//GSL versions of the constants, but G seems off..
//const double gamma_E = M_EULER;
//const double c = GSL_CONST_MKSA_SPEED_OF_LIGHT;
//const double G =GSL_CONST_MKSA_GRAVITATIONAL_CONSTANT;
//const double MSOL_SEC = GSL_CONST_MKSA_SOLAR_MASS*(GSL_CONST_MKSA_GRAVITATIONAL_CONSTANT/(c*c*c));
//const double MPC_SEC = GSL_CONST_MKSA_PARSEC*1e6/c;

struct fftw_outline
{
	fftw_complex *in, *out;
	fftw_plan p;
};



/*!\struct
 *
 * Container for set of spin weighted spherical harmonics (all spins are -2)
 */
template <class T>
struct sph_harm
{
	std::complex<T> Y22;
	std::complex<T> Y21;
	std::complex<T> Y20;
	std::complex<T> Y2m1;
	std::complex<T> Y2m2;
};

/*!\struct
 * \brief Structure for interfacing with the libraries
 *
 * Structure to interface with the libraries - Units are in solar masses and mpc
 *  contains the generation parameters including source parameters and theory parameters
 *
 *  *NOTE* not all the members of this structure need to be assigned for usage. In fact, some are reduntant. It's up to the user to determine what fields require an assignment. (Sorry)
 */

template<class T>
class gen_params_base
{
public:
	T x0 = 1;
	std::string cosmology="PLANCK15";
	/*!mass of the larger body in Solar Masses*/
	T mass1;
	/*!mass of the smaller body in Solar Masses*/
	T mass2;
	/*!Luminosity distance to the source*/
	T Luminosity_Distance;
	/*!Spin vector of the larger mass [Sx,Sy,Sz]*/
	T spin1[3];
	/*!Spin vector of the smaller mass [Sx,Sy,Sz]*/
	T spin2[3];
	/*!coalescence time of the binary*/
	T tc=0;
	/*! tidal deformability of the larger component*/
        T tidal1=-1;
        /*! tidal deformability of the smaller component*/
        T tidal2=-1;
	/*! symmetric tidal deformability*/
	T tidal_s=-1;
	/*! antisymmetric tidal deformability*/
	T tidal_a=-1;
	T tidal_weighted=-1;
        /*! Boolean flag indicating binary love relations should be used*/  
	bool tidal_love = true;
        /*! Boolean flag indicating error marginalization over residual EoS 
        * dependence in binary love relations should be performed. 
        */
        bool tidal_love_error = false; 
	T delta_tidal_weighted=-1;

	//Polarization angle
	T psi =0 ;
	/*!*angle between angular momentum and the total momentum */
	T incl_angle;

	/*! boolean flag indicating equatorial orientation coordinates should be used*/
	bool equatorial_orientation=false;
	/*! Equatorial Spherical angles for the orbital angular momentum*/
	T theta_l;
	/*! Equatorial Spherical angles for the orbital angular momentum*/
	T phi_l;

	/*! Boolean flag indicating local, horizon coordinates should be used*/
	bool horizon_coord=false;
	/*! Polar angle in detector-centered coordinates*/
	T theta;
	/*! azimuthal angle in detector-centered coordinates*/
	T phi;
	/*! Equatorial coordinates of source RA*/
	T RA;
	/*! Equatorial coordinates of source DEC*/
	T DEC;
	/*! Greenwich Mean Sidereal time (for detector orientation - start of data -- IN RADIANS NOT HOURS*/
	double gmst;
	/*! BOOL flag for early termination of NS binaries*/
	//bool NSflag;
	bool NSflag1=false;
	bool NSflag2=false;

  	//########## bumpy EOS parameters #######
        T bump_mag; 	/*! Magnitude of feature */
        T bump_width;    /*! Width of feature */
        T bump_offset;  /*! Location of feature */
        T plat; /*! Constant value (or plateau)  of cs^2 after feature*/
        T nbc1; /*! Central baryon number density of star 1 -- to replace mass in sampling*/
        T nbc2; /*! Central baryon number density of star 2 -- to replace mass in sampling*/
  
	/*! Flag to force the use of deprecated postmerger calculations -- ADOLC friendly*/
	bool dep_postmerger = false;
	/*! Reference frequency for PhenomPv2*/
	T f_ref=0;

	/*! Shift time detemines if times are shifted so coalescence is more accurately*/
	bool shift_time = true;
	/*! Shift time detemines if phic or phiRef is used*/
	bool shift_phase = true;
	T phiRef=0;

	bool sky_average=false;
	//###################################################
	//Either define all these parameters for Pv2, or define
	//all the source frame parameters above
	/*!thetaJ -- optional domain is [0,M_PI]  */
	T thetaJN = -10;

	T alpha0 = 0;

	T chip = -1;

	//Azimuthal angle of chip in plane
	T phip = -1;

	bool precess_reduced_flag=false;

	//LISA SPECIFIC OPTIONS
	T  LISA_alpha0=0;
	T  LISA_phi0=0;
	/*! Polar angle in ecliptic coordinates*/
	//T  LISA_thetal;
	//T  LISA_phil;
	/*! Azimuthal angle in ecliptic coordinates for the total angular momentum -- internal, and should not be specified*/
	T  theta_j_ecl;
	T  phi_j_ecl;

	//gIMR quantities
	int Nmod_beta=0;
	int Nmod_alpha=0;
	int Nmod_sigma=0;
	int Nmod_phi=0;

	int *betai=NULL;
	int *alphai=NULL;
	int *sigmai=NULL;
	int *phii=NULL;

	T *delta_beta=NULL;
	T *delta_alpha=NULL;
	T *delta_sigma=NULL;
	T *delta_phi=NULL;
	/*!ppE b parameter (power of the frequency) - vector for multiple modifications*/
	double *bppe=NULL;
	double *appe=NULL;
	T *alphappe=NULL;
	/*!ppE coefficient for the phase modification - vector for multiple modifications*/
	T *betappe=NULL;
	/*!Number of phase modificatinos*/
	int Nmod=0;
	//###################################################
	//T chi1_p = 0;
	//T chi2_p = 0;
	T chi1_l = 0;
	T chi2_l = 0;
	T phiJL = 0 ;
	T thetaJL = 0 ;
	T zeta_polariz =0;
	T phi_aligned = 0;

	T chil = 0;
	//###################################################
	
	int PNorder = 35;

	//coupling constants in one parameterization for Einstein Aether
	T ca_EA;
	T ctheta_EA;
	T cw_EA;
	T csigma_EA;
        bool EA_region1 = false; //region of parameter space with ctheta = 3ca(1+delta_ctheta)
        T delta_ctheta_EA; 
  	//###################################################
        //coupling constants in another parameterization for Einstein Aether
        T alpha1_EA;
        T alpha2_EA;
        T cbarw_EA;
        bool alpha_param = true; 
	//###################################################
        bool include_l1=false;
	gsl_spline *Z_DL_spline_ptr=NULL;

	gsl_interp_accel *Z_DL_accel_ptr=NULL;

	void print_properties()
	{
		std::cout<<"Source Properties: "<<std::endl;
		std::cout<<"mass1 | mass2 | DL: "<<this->mass1<<" | "<<this->mass2<<" | "<<this->Luminosity_Distance<<std::endl;
		std::cout<<"spin1 | spin2: "<<"["<<this->spin1[0]<<", "<< this->spin1[1]<<", "<<this->spin1[2]<< "]"<<" | "<<"["<<this->spin1[0]<<", "<< this->spin1[1]<<", "<<this->spin1[2]<< "]"<<std::endl;
		std::cout<<"psi | iota | RA | DEC | gmst: "<<this->psi<<" | "<<this->incl_angle<<" | "<<this->RA<<" | "<<this->DEC<<" | "<<this->gmst<<std::endl;
		std::cout<<"f_ref | shift_time | shift_phase: " <<this->f_ref<<" | "<<this->shift_time<<" | "<<this->shift_phase<<std::endl;
		std::cout<<"tc | phiRef: " <<this->tc<<" | "<<this->phiRef<<std::endl;
		std::cout<<"horizon_coord | equatorial_orientation | sky_average | dep_postmerger: " <<this->horizon_coord<<" | "<<this->equatorial_orientation<<" | "<<this->sky_average<< " | "<<this->dep_postmerger<<std::endl;
		std::cout<<"Cosmology: "<<this->cosmology<<std::endl;
		if(this->Nmod > 0 ){
			std::cout<<"Nmod: "<<this->Nmod<<std::endl;
			std::cout<<"Betas: ";
			for(int i = 0 ; i<this->Nmod; i++){
				std::cout<<this->betappe[i]<<" | ";

			}
			std::cout<<std::endl;
			std::cout<<"bppe: ";
			for(int i = 0 ; i<this->Nmod; i++){
				std::cout<<this->bppe[i]<<" | ";

			}
			std::cout<<std::endl;
		}
		if(this->Nmod_phi > 0 ){
			std::cout<<"Nmod_phi: "<<this->Nmod_phi<<std::endl;
			std::cout<<"delta phi: ";
			for(int i = 0 ; i<this->Nmod_phi; i++){
				std::cout<<this->delta_phi[i]<<" | ";

			}
			std::cout<<std::endl;
			std::cout<<"phi i: ";
			for(int i = 0 ; i<this->Nmod_phi; i++){
				std::cout<<this->phii[i]<<" | ";

			}
			std::cout<<std::endl;
		}
		if(this->Nmod_sigma > 0 ){
			std::cout<<"Nmod_sigma: "<<this->Nmod_sigma<<std::endl;
			std::cout<<"delta sigma: ";
			for(int i = 0 ; i<this->Nmod_sigma; i++){
				std::cout<<this->delta_sigma[i]<<" | ";

			}
			std::cout<<std::endl;
			std::cout<<"sigma i: ";
			for(int i = 0 ; i<this->Nmod_sigma; i++){
				std::cout<<this->sigmai[i]<<" | ";

			}
			std::cout<<std::endl;
		}
		if(this->Nmod_beta > 0 ){
			std::cout<<"Nmod_beta: "<<this->Nmod_beta<<std::endl;
			std::cout<<"delta beta: ";
			for(int i = 0 ; i<this->Nmod_beta; i++){
				std::cout<<this->delta_beta[i]<<" | ";

			}
			std::cout<<std::endl;
			std::cout<<"beta i: ";
			for(int i = 0 ; i<this->Nmod_beta; i++){
				std::cout<<this->betai[i]<<" | ";

			}
			std::cout<<std::endl;
		}
		if(this->Nmod_alpha > 0 ){
			std::cout<<"Nmod_alpha: "<<this->Nmod_alpha<<std::endl;
			std::cout<<"delta alpha: ";
			for(int i = 0 ; i<this->Nmod_alpha; i++){
				std::cout<<this->delta_alpha[i]<<" | ";

			}
			std::cout<<std::endl;
			std::cout<<"alpha i: ";
			for(int i = 0 ; i<this->Nmod_alpha; i++){
				std::cout<<this->alphai[i]<<" | ";

			}
			std::cout<<std::endl;
		}

	}

};

/*! \brief convience wrapper for the gen_params_base class
 *
 * If using the code in the intended way, this is all the user should ever have to use. Just allows the user to drop the template parameter
 *
 * Also implemented for backwards compatibility with previous versions of the code
 */
class gen_params:public gen_params_base<double>
{

};
/*!\brief To speed up calculations within the for loops, we pre-calculate reoccuring powers of M*F and Pi, since the pow() function is prohibatively slow
 *
 * Powers of PI are initialized once, and powers of MF need to be calculated once per for loop (if in the inspiral portion).
 *
 * use the functions precalc_powers_ins_amp, precalc_powers_ins_phase, precalc_powers_pi to initialize
 */
template <class T>
struct useful_powers
{
	T MFthird;
	T MFsixth;
	T MF7sixth;
	T MF2third;
	T MF4third;
	T MF5third;
	T MFsquare;
	T MF7third;
	T MF8third;
	T MFcube;
	T MFminus_5third;
	T MF3fourth;
	double PIsquare;
	double PIcube;
	double PIthird;
	double PI2third;
	double PI4third;
	double PI5third;
	double PI7third;
	double PIminus_5third;
};

/*!\struct
 * \brief For internal data transfers
 *
 * Structure to facililate parameter tranfers - All dimensionful quantities are in seconds
 */
template <class T>
struct source_parameters
{

	/*! mass of the larger component*/
	T mass1;
	/*! mass of the smaller component*/
	T mass2;
	/*! Total mass*/
	T M;
	/*Mass ratio*/
	T q;
	/*! z-Spin component of the larger body*/
	T spin1z;
	/*! z-Spin component of the smaller body*/
	T spin2z;
	/*! x-Spin component of the larger body*/
	T spin1x;
	/*! x-Spin component of the smaller body*/
	T spin2x;
	/*! y-Spin component of the larger body*/
	T spin1y;
	/*! y-Spin component of the smaller body*/
	T spin2y;
	/*!Chirp mass of the binary*/
	T chirpmass;
	/*!Symmetric mass ratio*/
	T eta;
	/*!Symmetric spin combination*/
	T chi_s;
	/*!Antisymmetric spin combination*/
	T chi_a;
	/*!Effective spin */
	T chi_eff;
	/*! PN spin */
	T chi_pn;
	/*!Luminoisity Distance*/
	T DL;
	/*! Delta mass comibination*/
	T delta_mass;
	/*!Ringdown frequency after merger*/
	T fRD;
	/*!Dampening frequency after merger*/
	T fdamp;
	/*! Transition Frequency 1 for the amplitude*/
	T f1;
	/*! Transition Frequency 2 for the amplitude*/
	T f3;
	/*! Transition frequency 1 for the phase*/
	T f1_phase;
	/*! Transition frequency 2 for the phase*/
	T f2_phase;
	/*! Coalescence time*/
	T tc;
	/*overall amplitude factor*/
	T A0;
	/*! Shift time detemines if phic or phiRef is used*/
	bool shift_phase = true;
	/*! Flag to force the use of deprecated postmerger calculations -- ADOLC friendly*/
	bool dep_postmerger = false;

	bool NSflag1;
	bool NSflag2;
	//############################################
	//PhenomP
	T s ;

	T chil;

	T chip;

	//Azimuthal angle of chip in plane
	T phip = -1;

	T f_ref=0;

	T phi_aligned;

	T incl_angle;

	T phiRef;

	T alpha0;

	T thetaJN;

	T zeta_polariz;

	T x0 = 1;

	int PNorder = 35;
	//#####################################
	T chi1_p = 0;
	T chi2_p = 0;
	T chi1_l = 0;
	T chi2_l = 0;
	T phiJL = 0 ;
	T thetaJL = -1 ;
	//######## NRT parameters ###############
	/*! tidal deformability of the larger component*/
        T tidal1=-1;
        /*! tidal deformability of the smaller component*/
        T tidal2=-1;
	/*! mass-weighted tidal deformability*/
	T tidal_weighted=-1;
	T delta_tidal_weighted=-1;
	T NRT_phase_coeff ; // Save for faster waveform calcs
	T quad1;
	T quad2;
	T oct1;
	T oct2;
	T ss_3p5PN_coeff;
	T NRT_amp_coefficient;

  	//########## bumpy EOS parameters #######
        T bump_mag; 	/*! Magnitude of feature */
        T bump_width;    /*! Width of feature */
        T bump_offset;  /*! Location of feature */
        T plat; /*! Constant value (or plateau)  of cs^2 after feature*/
        T nbc1; /*! Central baryon number density of star 1 -- to replace mass in sampling*/
        T nbc2; /*! Central baryon number density of star 2 -- to replace mass in sampling*/

	//######### ppE parameters ##############
	/*Beta factor for ppE formalism*/
	T *betappe;
	/*b power for ppE formalism*/
	double *bppe;
	double *appe;
	T *alphappe;

	/*! Number of modifications to phase*/
	int Nmod;

	//Spherical polar angles for the sky location relative to the detector in question
	T phi;

	T theta;

	T SP;

	T SL;

	bool sky_average;

	/*! Boolean -- shift time to 0 before shifting to tc or not*/
	bool shift_time = true;

	gsl_spline *Z_DL_spline_ptr=NULL;

	gsl_interp_accel *Z_DL_accel_ptr=NULL;

	std::string cosmology;
	//gIMR quantities
	int Nmod_beta=0;
	int Nmod_alpha=0;
	int Nmod_sigma=0;
	int Nmod_phi=0;
	int *betai;
	int *alphai;
	int *sigmai;
	int *phii;
	T *delta_beta;
	T *delta_alpha;
	T *delta_sigma;
	T *delta_phi;

  //Einstein-Aether stuff
  T kappa3_EA;
  T epsilon_x_EA;
  //coupling constants in one parameterization
  T ca_EA;
  T ctheta_EA;
  T cw_EA;
  T csigma_EA;
  //coupling constants in the other parameterization
  T c1_EA;
  T c2_EA;
  T c3_EA;
  T c4_EA;
  T c13_EA;
  T c14_EA;
  T cminus_EA;
  //Speeds of the different polarizations
  T cTsq_EA;
  T cVsq_EA;
  T cSsq_EA;
  T cT_EA;
  T cV_EA;
  T cS_EA;
  
  //Necessary functions of the c's
  T alpha1_EA;
  T alpha2_EA;
  T cbarw_EA;
  T beta1_EA;

  T beta2_EA;
  T Z_EA;
  T A1_EA;
  T A2_EA;
  T A3_EA;
  T B3_EA;
  T C_EA;
  T D_EA;
  //Center of mass velocity of binary (normally we'll just set this to zero)
  T V_x_EA;
  T V_y_EA;
  T V_z_EA;
  //Sensitivities for NSs
  T s1_EA;
  T s2_EA;
  T S_EA;

  //Compactness
  T compact1 = 0;
  T compact2 = 0;

  T alpha_ppE_2T_0_EA = 0 ;
  T gb1_EA = 0 ;
  T abL_EA = 0 ;
  T gX1_EA = 0 ;

  bool include_l1=false; //include l = 1 mode or not
  bool EA_nan_error_message = false; //check for nans --default must be false!
  bool alpha_param = true; //use alpha parameterization of EA theory
  bool EA_region1 = false; //region of parameter space with ctheta = 3ca(1+deltactheta)
  T delta_ctheta_EA = 0;
//static source_parameters<T> populate_source_parameters(gen_params_base<T> *param_in);
void populate_source_parameters(gen_params_base<T> *param_in);
static source_parameters<T> populate_source_parameters_old(
			T mass1,
			T mass2,
			T Luminosity_Distance,
			T *spin1,
			T *spin2,
			T phi_c,
			T t_c,
			bool sky_average) ;
};
double gsl_LU_lndet(double **matrix, int dim);
int mvn_sample(int samples, double *mean, double **cov, int dim, double **output );
int mvn_sample(int samples, double *mean, double **cov, int dim, gsl_rng *r,double **output );
void vector_union(std::vector<double> A, std::vector<double> B, std::vector<double> *C );
void matrix_multiply(double **A, double **B, double **C,int dim1, int dim2,int dim3);
template<class T>
void debugger_print(const char *file, const int line, T message);
std::string strip_path(std::string input);
template<class T, class U>
void transform_parameters(gen_params_base<T> *param_in, gen_params_base<U> *param_out);
template<class T, class U>
void transform_parameters(gen_params_base<T> *param_in, gen_params_base<U> **param_out);

int newton_raphson_method_1d(void(*f_on_fprime)(double x, double *func, double *func_prime,void *param),double initial_guess,double tolerance ,int max_iterations,void *parameters,double *solution);
template<class T>
T A0_from_DL(T chirpmass, T DL, bool sky_average);

template<class T>
T DL_from_A0(T chirpmass, T A0, bool sky_average);

void initiate_LumD_Z_interp(gsl_interp_accel **Z_DL_accel_ptr, gsl_spline **Z_DL_spline_ptr);
void free_LumD_Z_interp(gsl_interp_accel **Z_DL_accel_ptr, gsl_spline **Z_DL_spline_ptr);
adouble Z_from_DL_interp(adouble DL,gsl_interp_accel *Z_DL_accel_ptr, gsl_spline *Z_DL_spline_ptr);
double Z_from_DL_interp(double DL,gsl_interp_accel *Z_DL_accel_ptr, gsl_spline *Z_DL_spline_ptr);

double Z_from_DL(double DL, std::string cosmology);
double DL_from_Z(double Z, std::string cosmology);
double cosmology_interpolation_function(double x, double *coeffs, int interp_degree);
int cosmology_lookup(std::string cosmology);

template <class T>
T fcontact(T M_detector, T DL, std::string cosmology);

double powerlaw_from_uniform(double x0, double x1, double power, double uniform_random_number);
double gsl_maxwell_boltzmann_distribution(double sigma, gsl_rng *r);
template<class T>
void list_intersect_ptrs(T **A, int lenA,T **B, int lenB, T **C, int *lenC);
template<class T>
void list_intersect(T *A, int lenA,T *B,int lenB, T **C,int *lenC);
template<class T>
void list_intersect_value(T *A, int lenA,T *B,int lenB, T *C,int *lenC);
template<class T, class U>
void variance_list(T *list, int length,U *result);
template<class T,class U>
void mean_list(T *list, int length, U *result);

template<class T>
bool check_list(T j, T *list, int length);
template<class T>
int check_list_id(T j, T *list, int length);
template<class T>
T copysign_internal(T val, T  sign);
void rm_fisher_dim(double **input,int full_dim, double **output,  int reduced_dim, int *removed_dims);

template<class T>
void gsl_LU_matrix_invert(T **input, T **inverse, int dim);

int gsl_cholesky_matrix_invert(double **input, double **inverse, int dim);

int normalized_gsl_cholesky_matrix_invert(double **input, double **inverse, int dim);

adouble Z_from_DL(adouble DL, std::string cosmology);
adouble DL_from_Z(adouble Z, std::string cosmology);
adouble cosmology_interpolation_function(adouble x, double *coeffs, int interp_degree);

double std_omega(double RA, double std_RA, double std_DEC, double cov_RA_DEC);
void printProgress (double percentage);

void allocate_FFTW_mem_forward(fftw_outline *plan,int length);
void allocate_FFTW_mem_reverse(fftw_outline *plan,int length);
void deallocate_FFTW_mem(fftw_outline *plan);

double** allocate_2D_array( int dim1, int dim2);
int** allocate_2D_array_int( int dim1, int dim2);
void deallocate_2D_array(double **array, int dim1, int dim2);
void deallocate_2D_array(int **array, int dim1, int dim2);
double*** allocate_3D_array( int dim1, int dim2, int dim3);
int*** allocate_3D_array_int( int dim1, int dim2, int dim3);

void deallocate_3D_array(double ***array, int dim1, int dim2, int dim3);
void deallocate_3D_array(int ***array, int dim1, int dim2, int dim3);

void hann_window(double *window,int length);
void tukey_window(double *window,
		int length,
		double alpha);


template<class T>
void terr_pol_iota_from_equat_sph(T RA, T DEC, T thetaj, T phij, T *pol, T *iota);

template<class T>
void ecl_from_eq(T theta_eq, T phi_eq, T *theta_ecl, T *phi_ecl);

template<class T>
void equatorial_from_SF(T *SFvec,T thetal, T phil, T thetas, T phis, T iota, T phi_ref,T *EQvec);

double calculate_eta(double mass1, double mass2);
adouble calculate_eta(adouble mass1, adouble mass2);

double calculate_chirpmass(double mass1, double mass2);
adouble calculate_chirpmass(adouble mass1, adouble mass2);

double calculate_mass1(double chirpmass, double eta);
adouble calculate_mass1(adouble chirpmass, adouble eta);

double calculate_mass2(double chirpmass, double eta);
adouble calculate_mass2(adouble chirpmass, adouble eta);

double calculate_mass1_Mcq(double chirpmass, double q);
adouble calculate_mass1_Mcq(adouble chirpmass, adouble q);

double calculate_mass2_Mcq(double chirpmass, double q);
adouble calculate_mass2_Mcq(adouble chirpmass, adouble q);

template<class T>
void celestial_horizon_transform(T RA, T DEC, double gps_time, T LONG, T LAT,
				T *phi, T *theta);

template<class T>
T gps_to_GMST(T gps_time);
template<class T>
T gps_to_GMST_radian(T gps_time);
template<class T>
T gps_to_JD(T gps_time);

void decimal_to_HMS(double decimal, int *hour, int *min, double *second);

template<class T>
void transform_cart_sph(T *cartvec, T *sphvec);

template<class T>
void transform_sph_cart(T *sphvec, T *cartvec);

template<class T>
void unwrap_array(T *in, T *out, int len) ;
/*!\brief Trapezoidal sum rule to approximate discrete integral - Uniform spacing
 *
 * This version is faster than the general version, as it has half the function calls
 *
 * Something may be wrong with this function - had an overall offset for real data that was
 * fixed by using the simpsons rule - not sure if this was because of a boost in accuracy or
 * because something is off with the trapezoidal sum
 */
template <class T>
T trapezoidal_sum_uniform(double delta_x, int length, T *integrand);
template <class T>
T trapezoidal_sum_uniform(double delta_x, int length, T *integrand)
{
	T sum=0;
	for (int i =1;i<length-1;i++)
		sum += 2*integrand[i];
	sum += integrand[0]+integrand[length-1];
	return sum*delta_x/2;
}

/*!\brief Trapezoidal sum rule to approximate discrete integral - Non-Uniform spacing
 *
 * This version is slower than the uniform version, but will handle non-uniform spacing
 */
template <class T>
T trapezoidal_sum(double *delta_x, int length, T *integrand);
template <class T>
T trapezoidal_sum(double *delta_x, int length, T *integrand)
{
	T sum=0;
	for (int i =1;i<length-1;i++)
		sum += (integrand[i]+integrand[i-1])*(delta_x[i]-delta_x[i-1])/2;
	return sum;
}

/*!\brief Simpsons sum rule to approximate discrete integral - Uniform spacing
 *
 * More accurate than the trapezoidal rule, but must be uniform
 */
template <class T>
T simpsons_sum(double delta_x, int length, T *integrand);
template <class T>
T simpsons_sum(double delta_x, int length, T *integrand)
{
	T sum=0;
	for (int i =1;i<length-1;i++)
	{
		if( i % 2 == 0)
			sum += 2*integrand[i];
		else
			sum+=4*integrand[i];
	}
	sum += integrand[0]+integrand[length-1];
	return sum*delta_x/3;
}
long factorial(long num);
double pow_int(double base, int power);
adouble pow_int(adouble base, int power);

template <class T>
std::complex<T> cpolar(T mag, T phase);

template <class T>
std::complex<T> XLALSpinWeightedSphericalHarmonic(
                                   T theta,  /**< polar angle (rad) */
                                   T phi,    /**< azimuthal angle (rad) */
                                   int s,        /**< spin weight */
                                   int l,        /**< mode number l */
                                   int m         /**< mode number m */
    );
double cbrt_internal(double base);
adouble cbrt_internal(adouble base);
#endif
