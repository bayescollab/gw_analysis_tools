#ifndef IMRPHENOMD_H
#define IMRPHENOMD_H
#include <math.h>
#include <gsl/gsl_const_mksa.h>
#include <complex>
#include "util.h"

///*! \file 
// * Header file for utilities
// */
///*! \struct
// * Structure to facilitate IMRPhenomD parameter transfers
// */



template <class T>
struct lambda_parameters
{

	T rho[4];
	T v2;
	T gamma[4];
	T sigma[5];
	T beta[5];
	T alpha[7];
};


//##########################################################################
//function declarations 
template <class T>
class IMRPhenomD
{
public:

virtual ~IMRPhenomD() = default;

virtual void fisher_calculation_sky_averaged(double *frequency, 
			int length, 
			//double *parameters,
			gen_params *parameters,
			double **amplitude_deriv, 
			double **phase_deriv, 
			double *amplitude, 
			int *amp_tapes, 
			int *phase_tapes
			);
virtual void change_parameter_basis(T *old_param,
					T *new_param,
					bool sky_average
					);
virtual void construct_amplitude_derivative(double *frequencies, 
				int length,
				int dimension, 
				double **amplitude_derivative,
				source_parameters<double> *input_params,
				int *tapes=NULL
				);
virtual void construct_phase_derivative(double *frequencies, 
				int length,
				int dimension, 
				double **phase_derivative,
				source_parameters<double> *input_params,
				int *tapes=NULL
				);
virtual void amplitude_tape(source_parameters<double> *input_params, int *tape);

virtual void phase_tape(source_parameters<double> *input_params, int *tape);


virtual int construct_waveform(T *frequencies, int length, std::complex<T> *waveform, source_parameters<T> *params);

virtual std::complex<T> construct_waveform(T frequency, 				
			source_parameters<T> *params );

virtual int construct_amplitude(T *frequencies, int length, T *amplitude, source_parameters<T> *params);

//virtual T construct_amplitude(T frequency,  source_parameters<T> *params);

virtual int construct_phase(T *frequencies, int length, T *phase, source_parameters<T> *params);

virtual T build_amp(T f, 
		lambda_parameters<T> *lambda, 
		source_parameters<T> *params, 
		useful_powers<T> *pows,
		T *amp_coeff, 
		T *deltas);

virtual T build_phase(T f, 
		lambda_parameters<T> *lambda, 
		source_parameters<T> *params, 
		useful_powers<T> *pows,
		T *phase_coeff);

virtual T assign_lambda_param_element(source_parameters<T> *source_param,int i);

virtual void assign_lambda_param(source_parameters<T> *source_param, lambda_parameters<T> *lambda);

virtual void precalc_powers_ins(T f, T M, useful_powers<T> *Mf_pows);

virtual void precalc_powers_PI( useful_powers<T> *PI_pows);

virtual void precalc_powers_ins_phase(T f, T M, useful_powers<T> *Mf_pows);

virtual void precalc_powers_ins_amp(T f, T M, useful_powers<T> *Mf_pows);

virtual void assign_pn_amplitude_coeff(source_parameters<T> *source_param, T *coeff);

virtual void assign_static_pn_phase_coeff(source_parameters<T> *source_param, T *coeff);

virtual void assign_nonstatic_pn_phase_coeff(source_parameters<T> *source_param, 
					T *coeff, 
					T f);

virtual void assign_nonstatic_pn_phase_coeff_deriv(source_parameters<T> *source_param, 
					T *Dcoeff, 
					T f);
virtual void post_merger_variables(source_parameters<T> *source_param);

virtual void calc_fring( source_parameters<T> *source_params);
virtual void calc_fdamp( source_parameters<T> *source_params);
virtual void _calc_fring( source_parameters<T> *source_params);
virtual void _calc_fdamp( source_parameters<T> *source_params);
virtual T final_spin(source_parameters<T> *params);
virtual T FinalSpin0815_s(T eta, T s);
virtual T FinalSpin0815(T eta, T chi1,T chi2);
virtual T EradRational0815_s(T eta, T s);
virtual T EradRational0815(T eta, T chi1,T chi2);

virtual T fpeak(source_parameters<T> *params, lambda_parameters<T> *lambda);

virtual T amp_ins(T f, source_parameters<T> *param, T *pn_coeff, 
			lambda_parameters<T> *lambda,useful_powers<T> *pow);

virtual T Damp_ins(T f, source_parameters<T> *param, T *pn_coeff, lambda_parameters<T> *lambda);

virtual T phase_ins(T f, source_parameters<T> *param, T *pn_coeff, 
		lambda_parameters<T> *lambda, useful_powers<T> *pow);

virtual T Dphase_ins(T f, source_parameters<T> *param, T *pn_coeff, lambda_parameters<T> *lambda);

virtual T amp_mr(T f, source_parameters<T> *param, lambda_parameters<T> *lambda);

virtual T phase_mr(T f, source_parameters<T> *param, lambda_parameters<T> *lambda);

virtual T Damp_mr(T f, source_parameters<T> *param, lambda_parameters<T> *lambda);

virtual T Dphase_mr(T f, source_parameters<T> *param, lambda_parameters<T> *lambda);

virtual T amp_int(T f, source_parameters<T> *param, lambda_parameters<T> *lambda, T *deltas);

virtual T phase_int(T f, source_parameters<T> *param, lambda_parameters<T> *lambda );

virtual T Dphase_int(T f, source_parameters<T> *param, lambda_parameters<T> *lambda);

virtual void phase_connection_coefficients(source_parameters<T> *param, 
				lambda_parameters<T> *lambda, 
				T *pn_coeffs);

virtual T calculate_beta1(source_parameters<T> *param, lambda_parameters<T> *lambda, T *pn_coeffs);

virtual T calculate_beta0(source_parameters<T> *param, lambda_parameters<T> *lambda, T *pn_coeffs);

virtual T calculate_alpha1(source_parameters<T> *param, lambda_parameters<T> *lambda);

virtual T calculate_alpha0(source_parameters<T> *param, lambda_parameters<T> *lambda);

virtual void amp_connection_coeffs(source_parameters<T> *param, 
			lambda_parameters<T> *lambda, 
			T *pn_coeffs, 
			T *coeffs);

virtual T calculate_delta_parameter_0(T f1,T f2,T f3,T v1,
				T v2,T v3,T dd1,T dd3,T M);
virtual T calculate_delta_parameter_1(T f1,T f2,T f3,T v1,
				T v2,T v3,T dd1,T dd3,T M);
virtual T calculate_delta_parameter_2(T f1,T f2,T f3,T v1,
				T v2,T v3,T dd1,T dd3,T M);
virtual T calculate_delta_parameter_3(T f1,T f2,T f3,T v1,
				T v2,T v3,T dd1,T dd3,T M);
virtual T calculate_delta_parameter_4(T f1,T f2,T f3,T v1,
				T v2,T v3,T dd1,T dd3,T M);
};
//###########################################################################

/*!Numerically calibrated parameters from arXiv:1508.07253
 * see the table in the data directory for labeled version from lalsuite
 */
const double lambda_num_params[19][11]=
//Rho1
{{ 3931.8979897196696 ,  -17395.758706812805 ,  3132.375545898835 , 343965.86092361377 ,  -1.2162565819981997e6 ,  -70698.00600428853 ,
    1.383907177859705e6,  -3.9662761890979446e6 ,  -60017.52423652596 ,  803515.1181825735 , 
   -2.091710365941658e6 }, 
//Rho2
{ -40105.47653771657 ,  112253.0169706701 ,  23561.696065836168 ,  -3.476180699403351e6 , 
   1.137593670849482e7,  754313.1127166454 ,  -1.308476044625268e7 ,  3.6444584853928134e7 , 
   596226.612472288 ,  -7.4277901143564405e6 ,  1.8928977514040343e7 }, 
//Rho3
{ 83208.35471266537 ,  -191237.7264145924 , 
   -210916.2454782992 ,  8.71797508352568e6 ,  -2.6914942420669552e7 ,  -1.9889806527362722e6 ,  3.0888029960154563e7 , 
   -8.390870279256162e7 ,  -1.4535031953446497e6 ,  1.7063528990822166e7 ,  -4.2748659731120914e7 }, 
//v2
{ 0.8149838730507785 , 
   2.5747553517454658,  1.1610198035496786 ,  -2.3627771785551537 ,  6.771038707057573 ,  0.7570782938606834 ,  -2.7256896890432474 ,
    7.1140380397149965 ,  0.1766934149293479 ,  -0.7978690983168183 ,  2.1162391502005153 }, 
//Gamma1
{ 0.006927402739328343 , 
   0.03020474290328911 , 0.006308024337706171 ,  - 0.12074130661131138 , 0.26271598905781324 , 0.0034151773647198794 , 
   - 0.10779338611188374 , 0.27098966966891747,  0.0007374185938559283,  - 0.02749621038376281 , 
   0.0733150789135702 }, 
//Gamma2
{ 1.010344404799477,  0.0008993122007234548 , 0.283949116804459 ,  -4.049752962958005 , 
   13.207828172665366, 0.10396278486805426 ,  -7.025059158961947 ,  24.784892370130475, 0.03093202475605892 ,  -2.6924023896851663 ,
    9.609374464684983 }, 
//Gamma3
{ 1.3081615607036106 ,  -0.005537729694807678 ,  -0.06782917938621007 ,  -0.6689834970767117 , 
   3.403147966134083 ,  -0.05296577374411866 ,  -0.9923793203111362,  4.820681208409587,  -0.006134139870393713 , 
   -0.38429253308696365 ,  1.7561754421985984}, 
//Sigma1
{ 2096.551999295543 ,  1463.7493168261553,  1312.5493286098522 , 
   18307.330017082117 ,  -43534.1440746107 ,  -833.2889543511114 ,  32047.31997183187 ,  -108609.45037520859,  452.25136398112204 ,
    8353.439546391714 ,  -44531.3250037322 }, 
//Sigma2
{ -10114.056472621156 ,  -44631.01109458185 ,  -6541.308761668722 , 
   -266959.23419307504 ,  686328.3229317984 ,  3405.6372187679685 ,  -437507.7208209015 ,  1.6318171307344697e6 , 
   -7462.648563007646 ,  -114585.25177153319 ,  674402.4689098676 }, 
//Sigma3
{ 22933.658273436497 ,  230960.00814979506 , 
   14961.083974183695 ,  1.1940181342318142e6,  -3.1042239693052764e6 ,  -3038.166617199259 ,  1.8720322849093592e6, 
   -7.309145012085539e6 ,  42738.22871475411 ,  467502.018616601 ,  -3.064853498512499e6 }, 
//Sigma4
{ -14621.71522218357 , 
   -377812.8579387104 ,  -9608.682631509726 ,  -1.7108925257214056e6 ,  4.332924601416521e6,  -22366.683262266528 , 
   -2.5019716386377467e6 ,  1.0274495902259542e7 ,  -85360.30079034246 ,  -570025.3441737515 , 
   4.396844346849777e6}, 
//Beta1
{ 97.89747327985583 ,  -42.659730877489224 ,  153.48421037904913 ,  -1417.0620760768954 , 
   2752.8614143665027 ,  138.7406469558649 ,  -1433.6585075135881 ,  2857.7418952430758 ,  41.025109467376126 ,  -423.680737974639 , 
   850.3594335657173 }, 
//Beta2
{ -3.282701958759534 ,  -9.051384468245866 ,  -12.415449742258042 ,  55.4716447709787 , 
   -106.05109938966335 ,  -11.953044553690658 , 76.80704618365418 ,  -155.33172948098394 ,  -3.4129261592393263 ,  25.572377569952536 , 
   -54.408036707740465 }, 
//Beta3
{ -0.000025156429818799565 ,  0.000019750256942201327 , 
   -0.000018370671469295915 ,  0.000021886317041311973 ,  0.00008250240316860033 ,  7.157371250566708e-6  , 
   -0.000055780000112270685 ,  0.00019142082884072178 ,  5.447166261464217e-6 , 
   -0.00003220610095021982,  0.00007974016714984341 }, 
//Alpha1
{ 43.31514709695348 , 638.6332679188081 ,  -32.85768747216059 ,
    2415.8938269370315 ,  -5766.875169379177 ,  -61.85459307173841 ,  2953.967762459948 ,  -8986.29057591497 , 
   -21.571435779762044 , 981.2158224673428 ,  -3239.5664895930286 }, 
//Alpha2
{ -0.07020209449091723 ,  -0.16269798450687084 , 
   -0.1872514685185499 , 1.138313650449945 ,  -2.8334196304430046 ,  -0.17137955686840617 ,  1.7197549338119527 , 
   -4.539717148261272 ,  -0.049983437357548705 ,  0.6062072055948309 ,  -1.682769616644546 }, 
//Alpha3
{ 9.5988072383479 , 
   -397.05438595557433 ,  16.202126189517813 ,  -1574.8286986717037 ,  3600.3410843831093 ,  27.092429659075467 ,  -1786.482357315139 ,
    5152.919378666511 ,  11.175710130033895 ,  -577.7999423177481 ,  1808.730762932043}, 
//Alpha4
{ -0.02989487384493607 , 
   1.4022106448583738 ,  -0.07356049468633846 ,  0.8337006542278661 ,  0.2240008282397391 ,  -0.055202870001177226 , 
   0.5667186343606578 ,  0.7186931973380503 ,  -0.015507437354325743 ,  0.15750322779277187 , 
   0.21076815715176228 }, 
//Alpha5
{ 0.9974408278363099 ,  -0.007884449714907203 ,  -0.059046901195591035 ,  1.3958712396764088 , 
   -4.516631601676276 ,  -0.05585343136869692 , 1.7516580039343603 ,  -5.990208965347804 ,  -0.017945336522161195 , 
   0.5965097794825992 ,  -2.0608879367971804 }};
//From paper
//const double lambda_num_params[19][11]=
//{{ 3931.9 ,  -17395.8 ,  3132.38 ,  343966. ,  -1.21626e6 ,  -70698. ,
//    1.38391e6 ,  -3.96628e6 ,  -60017.5 ,  803515. , 
//   -2.09171e6 }, 
//{ -40105.5 ,  112253. ,  23561.7 ,  -3.47618e6 , 
//   1.1375900000000002e7 ,  754313. ,  -1.30848e7 ,  3.64446e7 , 
//   596227. ,  -7.42779e6 ,  1.8929e7 }, 
//{ 83208.4 ,  -191238. , 
//   -210916. ,  8.71798e6 ,  -2.69149e7 ,  -1.98898e6 ,  3.0888e7 , 
//   -8.39087e7 ,  -1.4535e6 ,  1.70635e7 ,  -4.27487e7 }, 
//{ 0.814984 , 
//   2.57476 ,  1.16102 ,  -2.36278 ,  6.77104 ,  0.757078 ,  -2.72569 ,
//    7.11404 ,  0.176693 ,  -0.797869 ,  2.11624 }, 
//{ 0.0069274 , 
//   0.0302047 ,  0.00630802 ,  -0.120741 ,  0.262716 ,  0.00341518 , 
//   -0.107793 ,  0.27099 ,  0.000737419 ,  -0.0274962 , 
//   0.0733151 }, 
//{ 1.01034 ,  0.000899312 ,  0.283949 ,  -4.04975 , 
//   13.2078 ,  0.103963 ,  -7.02506 ,  24.7849 ,  0.030932 ,  -2.6924 ,
//    9.60937 }, 
//{ 1.30816 ,  -0.00553773 ,  -0.0678292 ,  -0.668983 , 
//   3.40315 ,  -0.0529658 ,  -0.992379 ,  4.82068 ,  -0.00613414 , 
//   -0.384293 ,  1.75618 }, 
//{ 2096.55 ,  1463.75 ,  1312.55 , 
//   18307.3 ,  -43534.1 ,  -833.289 ,  32047.3 ,  -108609. ,  452.251 ,
//    8353.44 ,  -44531.3 }, 
//{ -10114.1 ,  -44631. ,  -6541.31 , 
//   -266959. ,  686328. ,  3405.64 ,  -437508. ,  1.63182e6 , 
//   -7462.65 ,  -114585. ,  674402. }, 
//{ 22933.7 ,  230960. , 
//   14961.1 ,  1.19402e6 ,  -3.10422e6 ,  -3038.17 ,  1.87203e6 , 
//   -7.30915e6 ,  42738.2 ,  467502. ,  -3.06485e6 }, 
//{ -14621.7 , 
//   -377813. ,  -9608.68 ,  -1.71089e6 ,  4.33292e6 ,  -22366.7 , 
//   -2.50197e6 ,  1.02745e7 ,  -85360.3 ,  -570025. , 
//   4.39684e6 }, 
//{ 97.8975 ,  -42.6597 ,  153.484 ,  -1417.06 , 
//   2752.86 ,  138.741 ,  -1433.66 ,  2857.74 ,  41.0251 ,  -423.681 , 
//   850.359 }, 
//{ -3.2827 ,  -9.05138 ,  -12.4154 ,  55.4716 , 
//   -106.051 ,  -11.953 ,  76.807 ,  -155.332 ,  -3.41293 ,  25.5724 , 
//   -54.408 }, 
//{ -0.000025156400000000002 ,  0.000019750300000000003 , 
//   -0.0000183707 ,  0.0000218863 ,  0.0000825024 ,  7.15737e-6 , 
//   -0.000055780000000000005 ,  0.000191421 ,  5.44717e-6 , 
//   -0.0000322061 ,  0.0000797402 }, 
//{ 43.3151 ,  638.633 ,  -32.8577 ,
//    2415.89 ,  -5766.88 ,  -61.8546 ,  2953.97 ,  -8986.29 , 
//   -21.5714 ,  981.216 ,  -3239.57 }, 
//{ -0.0702021 ,  -0.162698 , 
//   -0.187251 ,  1.13831 ,  -2.83342 ,  -0.17138 ,  1.71975 , 
//   -4.53972 ,  -0.0499834 ,  0.606207 ,  -1.68277 }, 
//{ 9.59881 , 
//   -397.054 ,  16.2021 ,  -1574.83 ,  3600.34 ,  27.0924 ,  -1786.48 ,
//    5152.92 ,  11.1757 ,  -577.8 ,  1808.73 }, 
//{ -0.0298949 , 
//   1.40221 ,  -0.0735605 ,  0.833701 ,  0.224001 ,  -0.0552029 , 
//   0.566719 ,  0.718693 ,  -0.0155074 ,  0.157503 , 
//   0.210768 }, 
//{ 0.997441 ,  -0.00788445 ,  -0.0590469 ,  1.39587 , 
//   -4.51663 ,  -0.0558534 ,  1.75166 ,  -5.99021 ,  -0.0179453 , 
//   0.59651 ,  -2.06089 }};
#endif
