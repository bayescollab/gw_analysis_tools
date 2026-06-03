#ifndef FISHER_H
#define FISHER_H
#include "util.h"
#include "quadrature.h"
#include <string>

/*! \file 
 */
using namespace std;

struct gsl_subroutine
{
	string detector;
	string reference_detector;
	string generation_method;
	string sensitivity_curve;
	gen_params *gen_params_in;
	int dim;
	int id1;
	int id2;
	int *waveform_tapes;
	int *time_tapes;
	int *phase_tapes;
	double *freq_boundaries;
	double *grad_freqs;
	int boundary_num;
	bool *log_factors;
	
};

void tape_waveform_gsl_subroutine(gsl_subroutine * params_packed);
void tape_time_gsl_subroutine(gsl_subroutine * params_packed);
void tape_phase_gsl_subroutine(gsl_subroutine * params_packed);
void prep_gsl_subroutine(gsl_subroutine *params_packed);

void fisher_numerical(double *frequency, 
	int length,
	string generation_method, 
	string detector, 
	string reference_detector, 
	double **output,
	int dimension, 
	//double *parameters,
	gen_params_base<double> *parameters,
	int order,
	int *amp_tapes = NULL,
	int *phase_tapes = NULL,
	double *noise = NULL,
	Quadrature *quadMethod = NULL
	);
void calculate_fisher_elements(double *frequency, 
	int length, 
	int dimension, 
	std::complex<double> **response_deriv, 
	double **output, 
	double *psd, 
	std::string integration_method, 
	double *weights, 
	bool log10_f);

void calculate_fisher_elements(
	double **output,
	std::complex<double> **response_deriv,
	double *psd,
	int dimension,
	Quadrature *quadMethod
);

void calculate_fisher_elements_batch(double *frequency, 
	int length, 
	int base_dimension, 
	int full_dimension, 
	std::complex<double> **response_deriv, 
	double **output,
	double *psd, 
	std::string integration_method, 
	double *weights, 
	bool log10_f);


void calculate_derivatives(std::complex<double>  **response_deriv, 
       	double *frequencies,
       	int length, 
       	int dimension, 
       	string detector, 
       	string reference_detector, 
       	string  gen_method,
       	gen_params_base<double> *parameters, 
	int order);

void fisher_autodiff(double *frequency, 
	int length,
	string generation_method, 
	string detector, 
	string reference_detector, 
	double **output,
	int dimension, 
	gen_params *parameters,
	std::string integration_method="SIMPSONS",
	double *weights = NULL,
	bool log10_f=false,
	double *noise = NULL,
	int *amp_tapes = NULL,
	int *phase_tapes = NULL
	);
void fisher_autodiff_gq_internal(double *frequency, 
	int length,
	string generation_method, 
	string detector, 
	double **output,
	int dimension, 
	//double *parameters,
	gen_params *parameters,
	double *weights,
	bool log_freq,
	int *amp_tapes = NULL,
	int *phase_tapes = NULL,
	double *noise = NULL
	);

void fisher_autodiff_batch_mod(double *frequency, 
	int length,
	string generation_method, 
	string detector, 
	string reference_detector, 
	double **output,
	int base_dimension, 
	int full_dimension, 
	gen_params *parameters,
	std::string integration_method="SIMPSONS",
	double *weights = NULL,
	bool log10_f=false,
	double *noise = NULL,
	int *amp_tapes = NULL,
	int *phase_tapes = NULL
	);


void PhenomP_fisher(double *frequency,
	int length,
	gen_params *parameters,
	std::complex<double> **waveform_derivative,
	int *amp_tapes,
	int *phase_tapes,
	std::string detector
	);

void construct_waveform_derivative(double *frequency, 
	int length,
	int dimension,
	std::complex<double> **waveform_deriv,
	source_parameters<double> *input_params,
	int *waveform_tapes
	);

void calculate_derivatives_autodiff(double *frequency,
	int length,
	int dimension,
	std::string generation_method,
	gen_params *parameters,
	std::complex<double> **waveform_deriv,
	int *waveform_tapes,/*<< Waveform tapes -- length=6*/
	std::string detector,
	bool autodiff_time_deriv,
	std::string reference_detector
	);
void time_phase_corrected_derivative_autodiff_full_hess(double **dt, 
	int length, 
	double *frequencies,
	gen_params_base<double> *params, 
	std::string generation_method, 
	int dimension, 
	bool correct_time);

void time_phase_corrected_derivative_autodiff_numerical(double **dt, int length, double *frequencies,gen_params_base<double> *params, std::string generation_method, int dimension, bool correct_time);

std::string local_generation_method(std::string generation_method);


void prep_fisher_calculation(double *parameters, bool *, double *,double*, int,gen_params_base<double> *input_params, std::string generation_method, int dim);

void detect_adjust_parameters( double *freq_boundaries,
	double *grad_freqs, 
	int *boundary_num,
	gen_params_base<double> *input_params, 
	std::string generation_method, 
	std::string detector, 
	int dim);

void unpack_parameters(double *parameters, 
	gen_params_base<double> *input_params, 
	std::string generation_method, 
	int dimension, 
	bool *log_factors);

template<class T>
void repack_parameters(T *avec_parameters, gen_params_base<T> *a_params, std::string generation_method, int dim);

template<class T>
void repack_non_parameter_options(gen_params_base<T> *waveform_params, gen_params_base<double> *input_params, std::string gen_method);

template<class T>
void deallocate_non_param_options(gen_params_base<T> *waveform_params, gen_params_base<double> *input_params, std::string gen_method);
double calculate_integrand_autodiff_gsl_subroutine(double frequency, void *params_in);
void fisher_autodiff_gsl_integration(double *frequency_bounds, 
	string generation_method, 
	string sensitivity_curve, 
	string detector, 
	string reference_detector, 
	double **output,
	double **error,
	int dimension, 
	gen_params *parameters,
	double abserr,
	double relerr
	);
void fisher_autodiff_gsl_integration(double *frequency_bounds, 
	string generation_method, 
	string sensitivity_curve, 
	string detector, 
	string reference_detector, 
	double **output,
	double **error,
	int dimension, 
	gen_params *parameters,
	double abserr,
	double relerr,
	std::string error_log,
	bool logerr
	);
void fisher_autodiff_gsl_integration_batch_mod(double *frequency_bounds, 
	string generation_method, 
	string sensitivity_curve, 
	string detector, 
	string reference_detector, 
	double **output,
	double **error,
	int base_dimension, 
	int full_dimension, 
	gen_params *parameters,
	double abserr,
	double relerr
	);
void fisher_autodiff_gsl_integration_batch_mod(double *frequency_bounds, 
	string generation_method, 
	string sensitivity_curve, 
	string detector, 
	string reference_detector, 
	double **output,
	double **error,
	int base_dimension, 
	int full_dimension, 
	gen_params *parameters,
	double abserr,
	double relerr,
	std::string error_log,
	bool logerr
	);
void ppE_theory_fisher_transformation(std::string original_method, 
	std::string new_method,
	int dimension, 
	gen_params_base<double> *param,
	double **old_fisher, 
	double **new_fisher);
void ppE_theory_transformation_jac(
	std::string original_method, 
	std::string new_method,
	int dimension, 
	gen_params_base<double> *param, 
	double **jac);
void ppE_theory_fisher_transformation_calculate_derivatives(
	std::string original_method, 
	std::string new_method,
	int dimension,
	int base_dim, 
	gen_params_base<double> *param, 
	double **derivatives);
void ppE_theory_covariance_transformation(std::string original_method, 
	std::string new_method,
	int dimension, 
	gen_params_base<double> *param,
	double **old_cov, 
	double **new_cov);


// new calcualtederivatives stuff ----------------------------------------------------------------------------------

void calculate_derivatives_new(std::complex<double> **response_deriv,
							   double *frequencies,
							   int length,
							   int dimension,
							   string detector,
							   string reference_detector,
							   string gen_method,
							   gen_params_base<double> *parameters,
							   int order);

/**
 * @class CentralDifferenceShfits
 * @brief Base class for central difference shifts.
 *
 */
class CentralDifferenceShifts
{
protected:
	int dimension;
	double *parameters_vec;
	bool *log_factors;

public:
	const double epsilon = 1e-8;

    CentralDifferenceShifts(int dim, std::string local_gen_method, gen_params_base<double> *parameters);
		virtual ~CentralDifferenceShifts();

	int get_dimension() { return dimension; }
	double *get_parameters_vec() { return parameters_vec; }
	bool *get_log_factors() { return log_factors; }

	/**
	 * @brief Shift the value of a parameter in the parameters vector to take its derivative.
	 * @param i	Index of the parameter to shift.
	 * @param local_gen_method	Name of the model.
	 */
	virtual void calculate_parameters_shift(int i, std::string local_gen_method) = 0;
	/**
	 * @brief Unshift the value of parameter in parameters vector
	 * @param i Index of the parameter to unshift
	 */
	virtual void calculate_parameters_unshift(int i) = 0;
};

/**
 * @class CentralDifferenceShfits_O2
 * @brief Order 2 class for central difference shifts
 */
class CentralDifferenceShifts_O2 : public CentralDifferenceShifts
{
protected:
	double *param_p;
	double *param_m;

public:
	CentralDifferenceShifts_O2(int dim, std::string gen_method, gen_params_base<double> *parameters);
	~CentralDifferenceShifts_O2() override;

	double *get_param_p() { return param_p; }
	double *get_param_m() { return param_m; }

	/**
	 * @copydoc CentralDifferenceShifts::calculate_parameters_shift
	 * @note order 2 version
	 */
	void calculate_parameters_shift(int i, std::string local_gen_method) override;
	/**
	 * @copydoc CentralDifferenceShifts::calculate_parameters_unshift
	 * @note order 2 version
	 */
	void calculate_parameters_unshift(int i) override;
};

/**
 * @class CentralDifferenceShfits_O4
 * @brief Order 4 class for central difference shifts
 */
class CentralDifferenceShifts_O4 : public CentralDifferenceShifts_O2
{
protected:
	double *param_pp;
	double *param_mm;

public:
	CentralDifferenceShifts_O4(int dim, std::string gen_method, gen_params_base<double> *parameters);
	~CentralDifferenceShifts_O4() override;

	double *get_param_pp() { return param_pp; }
	double *get_param_mm() { return param_mm; }

	/**
	 * @copydoc CentralDifferenceShifts::calculate_parameters_shift
	 * @note order 4 version
	 */
	void calculate_parameters_shift(int i, std::string local_gen_method) override;
	/**
	 * @copydoc CentralDifferenceShifts::calculate_parameters_unshift
	 * @note order 4 version
	 */
	void calculate_parameters_unshift(int i) override;
};

/**
 * @class FiniteFisherDerivatives
 * @brief Base class for the finite-difference methods for computing Fisher derivatives.
 *
 */
class FiniteFisherDerivatives
{
protected:
	CentralDifferenceShifts *centralDifferenceShifts; // if its general class CentralDifferenceShifts then can't call child methods
	const double epsilon = 1e-8;
	int length;
	std::string local_gen_method;

public:
	FiniteFisherDerivatives(int len, int dim, int order, std::string local_gen_method, std::string gen_method, gen_params_base<double> *parameters);
	virtual ~FiniteFisherDerivatives();

	void calculate_derivatives(
		std::complex<double> **response_deriv,
		double *frequencies,
		int length,
		int dimension,
		string detector,
		string reference_detector,
		string gen_method,
		gen_params_base<double> *parameters);

	/**
	 * @brief Calculate the response quantities needed for the Fisher derivatives.
	 * @param waveform_params
	 * @param gen_method
	 * @param frequencies
	 * @param local_gen_method
	 * @param reference_detector
	 * @param detector
	 */
	virtual void calculate_response(gen_params &waveform_params, string gen_method, double *frequencies,
									std::string local_gen_method, string reference_detector, string detector) = 0;

	/**
	 * @brief Calculate the response derivative w.r.t. a single parameter.
	 * @param[in] i					The parameter being varied.
	 * @param[in] local_gen_method	Name of the model.
	 * @param[out] response_deriv	The row-array to be filled with the derivative at each frequency point.
	 */
	virtual void calculate_response_derivative(int i, std::string &local_gen_method,
											   std::complex<double> **response_deriv) = 0;
};

/**
 * @class AmplitudePhaseDerivatives
 * @brief intermediate base class for sky-average based calculations, child of FiniteFisherDerivatives
 */
class AmplitudePhaseDerivatives : public FiniteFisherDerivatives
{
protected:
	double *amplitude;

public:
	AmplitudePhaseDerivatives(int len, int dim, int order, double *frequencies,
							  std::string local_gen_method, std::string gen_method, gen_params_base<double> *parameters);
	~AmplitudePhaseDerivatives() override;

	/**
	 * @copydoc FiniteFisherDerivatives::calculate_response
	 * @note parameters reference_detector and detector are unused in this implementation (required by FiniteFisherDerivatives)
	 */
	virtual void calculate_response(gen_params &, string, double *, std::string, string, string) = 0;
	virtual void calculate_response_derivative(int i, std::string &local_gen_method,
											   std::complex<double> **response_deriv) = 0;
};

/**
 * @class AmplitudePhaseDerivatives_O2
 * @brief Order 2 class for AmplitudePhaseDerivatives
 */
class AmplitudePhaseDerivatives_O2 : public AmplitudePhaseDerivatives
{
protected:
	double *amplitude_plus;
	double *phase_plus;
	double *phase_plusc;
	double *amplitude_minus;
	double *phase_minus;
	double *phase_minusc;

public:
	AmplitudePhaseDerivatives_O2(int len, int dim, int order, double *frequencies,
								 std::string local_gen_method, std::string gen_method, gen_params_base<double> *parameters);
	~AmplitudePhaseDerivatives_O2() override;

	void calculate_response(gen_params &waveform_params, string gen_method, double *frequencies,
							std::string local_gen_method, string reference_detector, string detector) override;
	void calculate_response_derivative(int i, std::string &local_gen_method,
									   std::complex<double> **response_deriv) override;
};

/**
 * @class AmplitudePhaseDerivatives_O4
 * @brief Order 4 class for AmplitudePhaseDerivatives
 */
class AmplitudePhaseDerivatives_O4:public AmplitudePhaseDerivatives_O2 {
  protected:
    double *amplitude_plus_plus; 
		double *amplitude_minus_minus;
		double *phase_plus_plus;
		double *phase_minus_minus;
		double *phase_plus_plusc;
		double *phase_minus_minusc;

  public:
    AmplitudePhaseDerivatives_O4(int len, int dim, int order, double* frequencies, 
			std::string local_gen_method, std::string gen_method, gen_params_base<double>* parameters);
		~AmplitudePhaseDerivatives_O4() override;

    void calculate_response(gen_params& waveform_params, string gen_method, double* frequencies, 
			std::string local_gen_method, string reference_detector, string detector) override;
    
		void calculate_response_derivative(int i, std::string& local_gen_method, 
			std::complex<double>  **response_deriv) override;
	};

/**
 * @class FullResponseDerivatives_O2
 * @brief Order 2 class for general case of FiniteFisherDerivatives
 */
class FullResponseDerivatives_O2:public FiniteFisherDerivatives {
  protected:
		std::complex<double> *response_plus;
		std::complex<double> *response_minus;

  public:
    FullResponseDerivatives_O2(int len, int dim, int order, std::string local_gen_method, std::string gen_method, gen_params_base<double> *parameters);
		~FullResponseDerivatives_O2() override;

		void calculate_response(gen_params& waveform_params, string gen_method, double* frequencies, 
			std::string local_gen_method, string reference_detector, string detector) override;
		void calculate_response_derivative(int i, std::string& local_gen_method, 
			std::complex<double>** response_deriv) override;
	};

/**
 * @class FullResponseDerivatives_O4
 * @brief Order 4 class for general case of FiniteFisherDerivatives
 */
class FullResponseDerivatives_O4:public FullResponseDerivatives_O2 {
  protected:
		std::complex<double> *response_plus_plus;
		std::complex<double> *response_minus_minus;

  public:
    FullResponseDerivatives_O4(int len, int dim, int order, std::string local_gen_method, std::string gen_method, gen_params_base<double> *parameters);
		~FullResponseDerivatives_O4() override;
    

		void calculate_response(gen_params& waveform_params, string gen_method, double* frequencies, 
			std::string local_gen_method, string reference_detector, string detector) override;
		void calculate_response_derivative(int i, std::string& local_gen_method,
			 std::complex<double>** response_deriv) override;
	};


//------------------------------------------------------------------------------------------------------------------
#endif
