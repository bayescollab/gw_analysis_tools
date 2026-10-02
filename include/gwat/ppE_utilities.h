#ifndef PPE_UTILITIES_H
#define PPE_UTILITIES_H
#include <functional>
#include "util.h"
#include "waveform_util.h"

//typedef T (*beta_fn)(source_parameters<T> *);
template <class T>
//using beta_fn = T (*)(source_parameters<T> *);
using beta_fn = std::function<T(source_parameters<T> *)>;



template<class T>
struct theory_ppE_map{
	double *bppe=NULL;
	int Nmod;
	bool extra_polarizations=false;
	//T (*beta_fn)(source_parameters<T> *);
	beta_fn<T> *beta_fns=NULL;
	beta_fn<T> **beta_fns_ptrs=NULL;//Only necessary for certain theories
	std::string ppE_method;
};


int check_num_polar(std::string generation_method);
bool check_mod(std::string generation_method);
bool check_theory_support(std::string generation_method);

template<class T>
void assign_mapping(std::string generation_method,theory_ppE_map<T> *mapping,gen_params_base<T> *params_in);
template<class T>
void deallocate_mapping(theory_ppE_map<T> *mapping);

template<class T>
T PNSeries_beta(int term,source_parameters<T> *param);

template<class T>
T ppEAlt_beta(int term,source_parameters<T> *param);

template<class T>
T dCS_beta(source_parameters<T> *param);
template<class T>
T dCS_phase_factor(source_parameters<T> *param);

//#####################################################
template<class T>
T EdGB_HO_0PN_beta( source_parameters<T> *param);

template<class T>
T EdGB_HO_1PN_beta( source_parameters<T> *param);

template<class T>
T EdGB_HO_2PN_beta( source_parameters<T> *param);
//#####################################################

template<class T>
T EdGB_beta( source_parameters<T> *param);

template<class T>
T EdGB_phase_factor( source_parameters<T> *param);

template<class T>
T EdGB_GHO_betav1( source_parameters<T> *param);

template<class T>
T EdGB_GHO_betav2( source_parameters<T> *param);

template<class T>
T EdGB_GHO_betav3( source_parameters<T> *param);

template<class T>
T ExtraDimension_beta( source_parameters<T> *param);
	
template<class T>
T BHEvaporation_beta( source_parameters<T> *param);

template<class T>
T TVG_beta( source_parameters<T> *param);

template<class T>
T DipRad_beta( source_parameters<T> *param);

template<class T>
T NonComm_beta( source_parameters<T> *param);

/*template<class T>
T EA_fully_restricted_phase0(source_parameters<T> *p);
template<class T>
T EA_fully_restricted_phase1(source_parameters<T> *p);
*/
template<class T>
T ModDispersion_beta( source_parameters<T> *param);
int dispersion_lookup(double alpha);
template <class T>
T cosmology_interpolation_function_MD(T x, double *coeffs, int interp_degree);
template <class T>
T DL_from_Z_MD(T Z, double alpha);

//template<class T>
//void pre_calculate_EA_factors(source_parameters<T> *p);


#endif
