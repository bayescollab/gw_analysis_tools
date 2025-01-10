#include "ppE_utilities.h"
#include "util.h"
#include "waveform_generator.h"
#include "D_Z_Config_modified_dispersion.h"
/*! \brief Utilities for ppE parameterizations of gravitational waves
 *
 * Supported theories:
 *
 * dCS -- input beta is \alpha^2 in seconds^4
 * 
 * EdGB -- input beta is \alpha^2 in seconds^4
 * 
 * EdGB_GHOv1 -- EdGB Generic Higher Order version 1-- input beta is \alpha^2 in seconds^4 and \gamma, which is dimensionless
 *
 * EdGB_GHOv2 -- EdGB Generic Higher Order version 2-- input beta is \alpha^2 in seconds^4 and \gamma, which is dimensionless
 * 
 * EdGB_GHOv3 -- EdGB Generic Higher Order version 3-- input beta is \alpha^2 in seconds^4 and \gamma, which is dimensionless
 *
 * PNSeries_ppE_*_Inspiral -- creates a series of generic modifications that take the form betappe[0]* ( u^(bppe[0]/3) + betappe[1] * u^(bppe[1]/3) . . .) -- for inspiral only
 *
 * PNSeries_ppE_*_IMR -- Creates a ppE waveform, but with the expansion parameter of (pi M f)^1/3 instead of (pi Mc f)^1/3 -- IMR
 *
 * PNSeries_ppE_*_Inspiral -- Creates a ppE waveform, but with the expansion parameter of (pi M f)^1/3 instead of (pi Mc f)^1/3 -- Inspiral
 *
 * ppEAlt_*_IMR -- creates a series of generic modifications that take the form betappe[0]* ( u^(bppe[0]/3) + betappe[1] * u^(bppe[1]/3) . . .) for full IMR
 *
 * ExtraDimension -- input beta l^2 in seconds^2
 * 
 * BHEvaporation -- input dm / dt (dimensionless)
 *
 * TVG -- Time varying G -- input beta is \dot{G}_{z} in Hz
 *
 * DipRad -- Generic Dipole Radiation -- input beta is \delta \dot{E} (dimensionless)
 *
 * NonComm -- Noncommutatitve gravity -- input beta is \Lambda^2 in units of Planck energies (dimensionless), where \sqrt(\Lambda) defines the energy scale of noncommutativity 
 *
 * ModDispersion -- Modified Dispersion -- input beta is A_alpha (eV^(2-alpha))
 *
 * EA_fully_restricted_v1 -- Einstein-Aether neglecting other polarizations, higher harmonics, and amplitude corrections  -- input betas are c1, c2, c3, c4, tc_T (time of coalescence of for tensor mode)
 */

template<class T>
void extra_modifications(std::string generation_method,gen_params_base<T> *gp, source_parameters<T> *p, waveform_polarizations<T> *wp,T *freqs, int length)
{
  /*	if(generation_method.find("EA_fully_restricted_v1") != std::string::npos){
		//debugger_print(__FILE__,__LINE__,"POST");
		source_parameters<T> temp_sp;
		temp_sp.populate_source_parameters(gp);
		temp_sp.phiRef = gp->phiRef;
		temp_sp.f_ref = gp->f_ref;
		temp_sp.shift_phase = gp->shift_phase;
		temp_sp.shift_time = gp->shift_time;
		temp_sp.incl_angle = gp->incl_angle;
		temp_sp.betappe = gp->betappe;
		pre_calculate_EA_factors(&temp_sp);
		return EA_fully_restricted_v1_additional_modifications(&temp_sp,wp,freqs,length);
	}
  */
	return ;
}
template void extra_modifications(std::string, gen_params_base<double> * gp,source_parameters<double> *, waveform_polarizations<double> *,double *, int );
template void extra_modifications(std::string, gen_params_base<adouble> * gp,source_parameters<adouble> *, waveform_polarizations<adouble> *, adouble *, int);


bool check_mod(std::string generation_method)
{
	if(generation_method.find("ppE") != std::string::npos)
	{
		return true;
	}
	if(generation_method.find("gIMR") != std::string::npos)
	{
		return true;
	}
	if(generation_method.find("EA") != std::string::npos)
	{
		return true;
	}
	if(check_theory_support(generation_method)){
		return true;
	}
	return false;
}
bool check_theory_support(std::string generation_method)
{
	if(generation_method.find("dCS")!=std::string::npos){
		return true;
	}
	if(generation_method.find("EdGB")!=std::string::npos){
		return true;
	}
	if(generation_method.find("EdGB_HO")!=std::string::npos){
		return true;
	}
	if(generation_method.find("EdGB_HO_LO")!=std::string::npos){
		return true;
	}
	if(generation_method.find("EdGB_GHO")!=std::string::npos){
		return true;
	}
	if(generation_method.find("ExtraDimension")!=std::string::npos){
		return true;
	}
	if(generation_method.find("BHEvaporation")!=std::string::npos){
		return true;
	}
	if(generation_method.find("TVG")!=std::string::npos){
		return true;
	}
	if(generation_method.find("DipRad")!=std::string::npos){
		return true;
	}
	if(generation_method.find("NonComm")!=std::string::npos){
		return true;
	}
	if(generation_method.find("ModDispersion")!=std::string::npos){
		return true;
	}
	//if(generation_method.find("EA_fully_restricted")!=std::string::npos){
	//if(generation_method.find("EA_IMRPhenomD_NRT")!=std::string::npos){
	//        return true;
	//}
	if(generation_method.find("PNSeries_ppE")!=std::string::npos){
		return true;
	}
	if(generation_method.find("ppEAlt")!=std::string::npos){
		return true;
	}

	if(generation_method.find("polarization_test")!=std::string::npos){
		return true;
	}
	return false;
} 
template<class T>
void deallocate_mapping(theory_ppE_map<T> *mapping){
	if(mapping->bppe){
		delete [] mapping->bppe;
	}
	if(mapping->beta_fns){
		delete [] mapping->beta_fns;
	}
	//Only used for a select few theories 
	if(mapping->beta_fns_ptrs){
		for(int i = 0 ;i < mapping->Nmod; i++){
			delete mapping->beta_fns_ptrs[i];
		}
		delete [] mapping->beta_fns_ptrs;
	}
	//if(mapping->add_mod_fn){
	//	delete mapping->add_mod_fn;
	//}
}
template void deallocate_mapping(theory_ppE_map<adouble> *mapping);
template void deallocate_mapping(theory_ppE_map<double> *mapping);

template<class T>
void assign_mapping(std::string generation_method,theory_ppE_map<T> *mapping, gen_params_base<T> *params_in)
{
	mapping->extra_polarizations = check_extra_polarizations(generation_method);
	bool ins=false;
	if(generation_method.find("dCS")!= std::string::npos){
		mapping->Nmod = 1;
		mapping->bppe = new double[1];
		mapping->beta_fns = new beta_fn<T>[mapping->Nmod]; 
		mapping->bppe[0] = -1;
		mapping->beta_fns[0] = [](source_parameters<T> *p){return dCS_beta(p);} ;
		ins = true;
	}
	else if(generation_method.find("EdGB")!= std::string::npos){
		if(generation_method.find("EdGB_HO")!= std::string::npos){
			mapping->Nmod = 3;
			mapping->bppe = new double[3];
			mapping->bppe[0] =-7;
			mapping->bppe[1] =-5;
			mapping->bppe[2] =-3;
			mapping->beta_fns = new beta_fn<T>[mapping->Nmod]; 
			mapping->beta_fns[0] = [](source_parameters<T> *p){return EdGB_HO_0PN_beta(p);} ;
			mapping->beta_fns[1] = [](source_parameters<T> *p){return EdGB_HO_1PN_beta(p);} ;
			mapping->beta_fns[2] = [](source_parameters<T> *p){return EdGB_HO_2PN_beta(p);} ;
		}
		if(generation_method.find("EdGB_HO_LO")!= std::string::npos){
			mapping->Nmod = 1;
			mapping->bppe = new double[1];
			mapping->bppe[0] =-7;
			mapping->beta_fns = new beta_fn<T>[mapping->Nmod]; 
			mapping->beta_fns[0] = [](source_parameters<T> *p){return EdGB_HO_0PN_beta(p);} ;
		}
		else if(generation_method.find("EdGB_GHOv1")!= std::string::npos){
			mapping->Nmod = 2;
			mapping->bppe = new double[2];
			mapping->bppe[0] =-7;
			mapping->bppe[1] =-5;
			mapping->beta_fns = new beta_fn<T>[mapping->Nmod]; 
			mapping->beta_fns[0] = [](source_parameters<T> *p){return EdGB_beta(p);} ;
			mapping->beta_fns[1] = [](source_parameters<T> *p){return EdGB_GHO_betav1(p);} ;
		}
		else if(generation_method.find("EdGB_GHOv2")!= std::string::npos){
			mapping->Nmod = 2;
			mapping->bppe = new double[2];
			mapping->bppe[0] =-7;
			mapping->bppe[1] =-5;
			mapping->beta_fns = new beta_fn<T>[mapping->Nmod]; 
			mapping->beta_fns[0] = [](source_parameters<T> *p){return EdGB_beta(p);} ;
			mapping->beta_fns[1] = [](source_parameters<T> *p){return EdGB_GHO_betav2(p);} ;
		}
		else if(generation_method.find("EdGB_GHOv3")!= std::string::npos){
			mapping->Nmod = 2;
			mapping->bppe = new double[2];
			mapping->bppe[0] =-7;
			mapping->bppe[1] =-5;
			mapping->beta_fns = new beta_fn<T>[mapping->Nmod]; 
			mapping->beta_fns[0] = [](source_parameters<T> *p){return EdGB_beta(p);} ;
			mapping->beta_fns[1] = [](source_parameters<T> *p){return EdGB_GHO_betav3(p);} ;
		}
		else{
			mapping->Nmod = 1;
			mapping->bppe = new double[1];
			mapping->bppe[0] =-7;
			mapping->beta_fns = new beta_fn<T>[mapping->Nmod]; 
			mapping->beta_fns[0] = [](source_parameters<T> *p){return EdGB_beta(p);} ;
		}
		ins = true;
	}
	else if(generation_method.find("ExtraDimension")!= std::string::npos){
		mapping->Nmod = 1;
		mapping->bppe = new double[1];
		mapping->bppe[0] =-13;
		mapping->beta_fns = new beta_fn<T>[mapping->Nmod]; 
		mapping->beta_fns[0] = [](source_parameters<T> *p){return ExtraDimension_beta(p);} ;
		ins = true;
	}
	//else if(generation_method.find("EA_fully_restricted_v1")!= std::string::npos){
	/*else if(generation_method.find("EA_IMRPhenomD_NRT")!= std::string::npos){
		//Need to pre-calculate EA variables here -- Needs to actually be in waveform_generator file..
		//pre_calculate_EA_factors(params_in);
		mapping->Nmod = 2;
		mapping->bppe = new double[2];
		mapping->bppe[0] =-7;
		mapping->bppe[1] =-5;
		mapping->beta_fns = new beta_fn<T>[mapping->Nmod]; 
		mapping->beta_fns[0] = [](source_parameters<T> *p){return EA_fully_restricted_phase0(p);} ;
		mapping->beta_fns[1] = [](source_parameters<T> *p){return EA_fully_restricted_phase1(p);} ;
		ins = true;
		}*/
	else if(generation_method.find("BHEvaporation")!= std::string::npos){
		mapping->Nmod = 1;
		mapping->bppe = new double[1];
		mapping->bppe[0] =-13;
		mapping->beta_fns = new beta_fn<T>[mapping->Nmod]; 
		mapping->beta_fns[0] = [](source_parameters<T> *p){return BHEvaporation_beta(p);} ;
		ins = true;
	}
	else if(generation_method.find("TVG")!= std::string::npos){
		mapping->Nmod = 1;
		mapping->bppe = new double[1];
		mapping->bppe[0] =-13;
		mapping->beta_fns = new beta_fn<T>[mapping->Nmod]; 
		mapping->beta_fns[0] = [](source_parameters<T> *p){return TVG_beta(p);} ;
		ins = true;
	}
	else if(generation_method.find("DipRad")!= std::string::npos){
		mapping->Nmod = 1;
		mapping->bppe = new double[1];
		mapping->bppe[0] =-7;
		mapping->beta_fns = new beta_fn<T>[mapping->Nmod]; 
		mapping->beta_fns[0] = [](source_parameters<T> *p){return DipRad_beta(p);} ;
		ins = true;
			
	}
	else if(generation_method.find("NonComm")!= std::string::npos){
		mapping->Nmod = 1;
		mapping->bppe = new double[1];
		mapping->bppe[0] =-1;
		mapping->beta_fns = new beta_fn<T>[mapping->Nmod]; 
		mapping->beta_fns[0] = [](source_parameters<T> *p){return NonComm_beta(p);} ;
		ins = true;
			
	}
	else if(generation_method.find("ModDispersion")!= std::string::npos){
		mapping->Nmod = 1;
		mapping->bppe = new double[1];
		mapping->bppe[0] =params_in->bppe[0];
		mapping->beta_fns = new beta_fn<T>[mapping->Nmod]; 
		mapping->beta_fns[0] = [](source_parameters<T> *p){return ModDispersion_beta(p);} ;
		ins = false;
			
	}
	else if(generation_method.find("PNSeries")!= std::string::npos){
		mapping->Nmod = params_in->Nmod;
		mapping->bppe = new double[params_in->Nmod];
		for(int i = 0 ; i<params_in->Nmod; i++){
			mapping->bppe[i] =params_in->bppe[i];
		}
		mapping->beta_fns = new beta_fn<T>[mapping->Nmod]; 
		mapping->beta_fns_ptrs = new beta_fn<T>*[mapping->Nmod]; 
		for(int i = 0 ; i<params_in->Nmod; i++){
			auto lam = new std::function<T(source_parameters<T> *)>( [i](source_parameters<T> *p){ return PNSeries_beta(i,p);} );
			mapping->beta_fns_ptrs[i] = lam;
			mapping->beta_fns[i] = *(lam); 
		}
		if(generation_method.find("Inspiral")!= std::string::npos){
			ins = true;
		}
		else{
			ins = false;
		}
			
	}
	else if(generation_method.find("ppEAlt")!= std::string::npos){
		mapping->Nmod = params_in->Nmod;
		mapping->bppe = new double[params_in->Nmod];
		for(int i = 0 ; i<params_in->Nmod; i++){
			mapping->bppe[i] =params_in->bppe[i];
		}
		mapping->beta_fns = new beta_fn<T>[mapping->Nmod]; 
		mapping->beta_fns_ptrs = new beta_fn<T>*[mapping->Nmod]; 
		for(int i = 0 ; i<params_in->Nmod; i++){
			auto lam = new std::function<T(source_parameters<T> *)>( [i](source_parameters<T> *p){ return ppEAlt_beta(i,p);} );
			mapping->beta_fns_ptrs[i] = lam;
			mapping->beta_fns[i] = *(lam); 
		}
		if(generation_method.find("Inspiral")!= std::string::npos){
			ins = true;
		}
		else{
			ins = false;
		}
			
	}
	else if(generation_method.find("polarization_test")!= std::string::npos){
		mapping->Nmod = params_in->Nmod;
		ins = true;
	}
	if(generation_method.find("Pv2")!=std::string::npos){
		if(ins){
			mapping->ppE_method = "ppE_IMRPhenomPv2_Inspiral";
		}
		else{
			mapping->ppE_method = "ppE_IMRPhenomPv2_IMR";
		}
	}
	else if(generation_method.find("PhenomD")!=std::string::npos){
		if(ins && generation_method.find("NRT") == std::string::npos){
			mapping->ppE_method = "ppE_IMRPhenomD_Inspiral";
		}
		else if(!ins && generation_method.find("NRT") == std::string::npos){
			mapping->ppE_method = "ppE_IMRPhenomD_IMR";
		}
		else if(ins && generation_method.find("NRT") != std::string::npos){
			mapping->ppE_method = "ppE_IMRPhenomD_NRT_Inspiral";
		}
		else if(!ins && generation_method.find("NRT") != std::string::npos){
			mapping->ppE_method = "ppE_IMRPhenomD_NRT_IMR";
		}

	}
	return ;
}
template void assign_mapping(std::string,theory_ppE_map<double>*,gen_params_base<double>*);
template void assign_mapping(std::string,theory_ppE_map<adouble>*,gen_params_base<adouble>*);


/*template<class T>
void EA_fully_restricted_v1_additional_modifications(source_parameters<T> *param, waveform_polarizations<T> *wp, T *freqs, int length)
{
	std::complex<T> *hall = new std::complex<T>[length];//Waveform common to all polarizations

	T ci = cos(param->incl_angle);
	T si = sin(param->incl_angle);
	for(int i = 0 ; i<length; i++){
		hall[i] = wp->hplus[i] / std::complex<T>(-1*(1+ci*ci),0) ;	
	}
	//T dtV = param->tc*(1./param->cT_EA-1./param->cV_EA)/(1-1./param->cT_EA);
	//T dtS = param->tc*(1./param->cT_EA-1./param->cS_EA)/(1-1./param->cT_EA);
	T dtV = (1.- param->cT_EA/param->cV_EA);
	T dtS = (1.- param->cT_EA/param->cS_EA);
	std::complex<T> shift_V ;
	std::complex<T> shift_S ;
	T u2m2=0;
	for(int i = 0 ; i<length; i++){
		u2m2 = pow(M_PI * param->chirpmass * freqs[i],-2./3.);
		//std::complex<T> time_shift = ( exp(std::complex<T>(0,2*M_PI*freqs[i]*dtV)+exp(std::complex<T>(0,2*M_PI*freqs[i]*dtS))));
		//wp->hplus[i]*=time_shift;
		//wp->hcross[i]*=time_shift;
		shift_V = exp(std::complex<T>(0,-2*M_PI*freqs[i] *param->tc* dtV));
		shift_S = exp(std::complex<T>(0,-2*M_PI*freqs[i] *param->tc* dtS));

		wp->hplus[i] *= (std::complex<T>(1 + u2m2 * param->alpha_ppE_2T_0_EA,0) + shift_V + shift_S);
		wp->hcross[i] *= (std::complex<T>(1 + u2m2 * param->alpha_ppE_2T_0_EA,0) + shift_V + shift_S);

		wp->hb[i] = hall[i] * std::complex<T>(0.5,0) * u2m2* param->alpha_ppE_2T_0_EA* param->gb1_EA * (std::complex<T>(1,0) - param->abL_EA)*shift_S*si*si;
		wp->hx[i] = hall[i] *u2m2 * param->alpha_ppE_2T_0_EA*param->gX1_EA * si * shift_V * ci;
		wp->hy[i] = hall[i] *u2m2 * param->alpha_ppE_2T_0_EA*param->gX1_EA * si * shift_V * std::complex<T>(0,1);
	} 
	delete [] hall;
	return ;
}
template void EA_fully_restricted_v1_additional_modifications(source_parameters<double> *param, waveform_polarizations<double> *wp, double *, int);
template void EA_fully_restricted_v1_additional_modifications(source_parameters<adouble> *param, waveform_polarizations<adouble> *wp, adouble *, int);
*/

/* \brief ppEAlt conversion
 *
 * The basis is M f and not \mathcal{M} f, so there's a conversion to account for using ppE waveforms
 */
template<class T>
T ppEAlt_beta(int term,source_parameters<T> *param)
{
	T out = 0;
	T chirpmass = calculate_chirpmass(param->mass1,param->mass2);
	T total_m = param->mass1 + param->mass2;
	
	if(term == 0 ){
		out = param->betappe[0] * pow(total_m/chirpmass,param->bppe[0]/3.);
	}
	else {
		out = param->betappe[term]* pow(total_m/chirpmass,param->bppe[term]/3.);
	}
	return out;
}
template double ppEAlt_beta(int ,source_parameters<double> *);
template adouble ppEAlt_beta(int , source_parameters<adouble> *);

/* \brief PNSeries conversion
 *
 * The basis is M f and not \mathcal{M} f, so there's a conversion to account for using ppE waveforms
 */
template<class T>
T PNSeries_beta(int term,source_parameters<T> *param)
{
	T out = 0;
	T chirpmass = calculate_chirpmass(param->mass1,param->mass2);
	T total_m = param->mass1 + param->mass2;
	
	if(term == 0 ){
		out = param->betappe[0] * pow(total_m/chirpmass,param->bppe[0]/3.);
	}
	else {
		out = param->betappe[0]* param->betappe[term]* pow(total_m/chirpmass,param->bppe[term]/3.);
	}
	return out;
}
template double PNSeries_beta(int ,source_parameters<double> *);
template adouble PNSeries_beta(int , source_parameters<adouble> *);
/*
template<class T>
T EA_fully_restricted_phase0(source_parameters<T> *p)
{
	pre_calculate_EA_factors(p);
	T out = 0;
	out = -3./224. * 1./p->kappa3_EA * pow(p->eta,2./5.) * p->epsilon_x_EA;
	
	//Minus one comes from sign choice from paper
	return -1*out;
}
template double EA_fully_restricted_phase0(source_parameters<double> *);
template adouble EA_fully_restricted_phase0(source_parameters<adouble> *);

template<class T>
T EA_fully_restricted_phase1(source_parameters<T> *p)
{
	pre_calculate_EA_factors(p);
	T out = 0;
	out = -3./128. * ( -2./3. * ( p->s1_EA+p->s2_EA) - 1./2. * (p->c1_EA + p->c4_EA) + ( p->kappa3_EA - 1.)) ;
	//Minus one comes from sign choice from paper
	return -1*out;
}
template double EA_fully_restricted_phase1(source_parameters<double> *);
template adouble EA_fully_restricted_phase1(source_parameters<adouble> *);
*/
template<class T>
T dCS_beta(source_parameters<T> *param)
{
 	T M = param->M;	
	T DL = param->DL;
	T Z= Z_from_DL(DL/MPC_SEC,param->cosmology);
	T unredshiftedM = M/(1.+Z);
	T phase_mod = (param->betappe[0]);
	T out =  16.*M_PI*phase_mod/(pow_int(unredshiftedM,4)) *dCS_phase_factor(param);
	return out;
}
template double dCS_beta(source_parameters<double> *);
template adouble dCS_beta(source_parameters<adouble> *);

template<class T>
T dCS_phase_factor(source_parameters<T> *param)
{
	T g=0;
 	T M = param->M;	
 	T chirpmass = param->chirpmass;	
 	T eta = param->eta;	
	T coeff1 = -5./8192.;
	T coeff2 = 15075./114688.;
	T m1 = calculate_mass1(chirpmass,eta);
	T m2 = calculate_mass2(chirpmass,eta);
	T m = m1+m2;
	T chi1 = param->chi_s+param->chi_a;
	T chi2 = param->chi_s-param->chi_a;
	T s1temp = 2.+2.*pow_int(chi1,4) - 2.*sqrt((1.-chi1*chi1)) - chi1*chi1 * ((3. - 2.*sqrt(1.-chi1*chi1)));
	T s2temp = 2.+2.*pow_int(chi2,4) - 2.*sqrt((1.-chi2*chi2)) - chi2*chi2 * ((3. - 2.*sqrt(1.-chi2*chi2)));
	
	T s1  ;
	T s2  ;
	if( fabs(chi1) <DOUBLE_COMP_THRESH){
		s1 = 0;		
	}
	else{
		s1  = s1temp/(2.*chi1*chi1*chi1);
	}
	if( fabs(chi2) <DOUBLE_COMP_THRESH){
		s2 = 0;		
	}
	else{
		s2  = s2temp/(2.*chi2*chi2*chi2);
	}
	//Neutron stars don't source scalar charge
	if(param->NSflag1){s1 =0;}	
	if(param->NSflag2){s2 =0;}	
	g+=coeff1/(pow(eta,14./5.)) * pow((m1*s2 - m2 * s1),2.)/(m*m);
	g+=coeff2/(pow(eta,14./5.)) * (m2*m2* chi1*chi1 - 350./201. * m1*m2*chi1*chi2 + m1*m1 * chi2*chi2)/(m*m);

	return g;
}
template double dCS_phase_factor(source_parameters<double> *);
template adouble dCS_phase_factor(source_parameters<adouble> *);

//###############################################################
template<class T>
T EdGB_HO_0PN_beta( source_parameters<T> *param)
{
 	T M = param->M;	
	T DL = param->DL;
	T Z= Z_from_DL(DL/MPC_SEC,param->cosmology);
	T unredshiftedM = M/(1.+Z);
	T alphaSq = 16.*M_PI*param->betappe[0];
	//T alphaSq = param->betappe[0];
	//return 16.*M_PI*phase_mod/(pow_int(unredshiftedM,4)) * EdGB_phase_factor(param);
	T out =  -5.*alphaSq/pow_int(unredshiftedM,4) / 7168. / pow(param->eta,18./5.) *( 4 *param->eta -1)  ;
	return out;
} 

template<class T>
T EdGB_HO_1PN_beta( source_parameters<T> *param)
{
 	T M = param->M;	
	T DL = param->DL;
	T Z= Z_from_DL(DL/MPC_SEC,param->cosmology);
	T unredshiftedM = M/(1.+Z);
	T alphaSq = 16.*M_PI*param->betappe[0];
	//T alphaSq = param->betappe[0];
	//return 16.*M_PI*phase_mod/(pow_int(unredshiftedM,4)) * EdGB_phase_factor(param);
	T out =  -5.*alphaSq/pow_int(unredshiftedM,4) / 688128. / pow_int(param->eta,4) *( 685. - 3916*param->eta + 2016* param->eta * param->eta)  ;
	return out;
} 

template<class T>
T EdGB_HO_2PN_beta( source_parameters<T> *param)
{
 	T M = param->M;	
	T DL = param->DL;
	T Z= Z_from_DL(DL/MPC_SEC,param->cosmology);
	T unredshiftedM = M/(1.+Z);
	T alphaSq = 16.*M_PI*param->betappe[0];
	//T alphaSq = param->betappe[0];
	//return 16.*M_PI*phase_mod/(pow_int(unredshiftedM,4)) * EdGB_phase_factor(param);
	T out =  5.*alphaSq/pow_int(unredshiftedM,4) / 387072. / pow(param->eta,22./5.) *pow_int( 1- 2. * param->eta,2)*(995. + 952.*param->eta)  ;
	return out;
} 

//###############################################################
template<class T>
T EdGB_beta( source_parameters<T> *param)
{
 	T M = param->M;	
	T DL = param->DL;
	T Z= Z_from_DL(DL/MPC_SEC,param->cosmology);
	T unredshiftedM = M/(1.+Z);
	T phase_mod = param->betappe[0];
	//return 16.*M_PI*phase_mod/(pow_int(unredshiftedM,4)) * EdGB_phase_factor(param);
	T out =  16.*M_PI*phase_mod/(pow_int(unredshiftedM,4)) * EdGB_phase_factor(param);
	return out;
} 
template adouble EdGB_beta(source_parameters<adouble> *);
template double EdGB_beta(source_parameters<double> *);


template<class T>
T EdGB_phase_factor( source_parameters<T> *param)
{
 	T M = param->M;	
 	T chirpmass = param->chirpmass;	
 	T eta = param->eta;	
	T m1 = calculate_mass1(chirpmass, eta);
	T m2 = calculate_mass2(chirpmass, eta);
	T chi1 = param->chi_s + param->chi_a;
	T chi2 = param->chi_s - param->chi_a;
	T temp1 = 2.*(sqrt(1.-chi1*chi1) - 1. + chi1*chi1);
	T temp2 = 2.*(sqrt(1.-chi2*chi2) - 1. + chi2*chi2);
	T s1;
        T s2;
	if( fabs(chi1) <DOUBLE_COMP_THRESH){
		s1 = 0;		
	}
	else{
		s1  = temp1/(chi1*chi1);
	}
	if( fabs(chi2) <DOUBLE_COMP_THRESH){
		s2 = 0;		
	}
	else{
		s2  = temp2/(chi2*chi2);
	}
	if(param->NSflag1){debugger_print(__FILE__,__LINE__,"NS 1");s1 =0;}	
	if(param->NSflag2){s2 =0;}	
	return (-5./7168.)* pow_int((m1*m1 * s2 - m2*m2 * s1),2) / (pow_int(M,4) * pow(eta,(18./5)));
} 
template adouble EdGB_phase_factor(source_parameters<adouble> *);
template double EdGB_phase_factor(source_parameters<double> *);

template<class T>
T EdGB_GHO_betav1( source_parameters<T> *param)
{
 	T M = param->M;	
	T DL = param->DL;
	T Z= Z_from_DL(DL/MPC_SEC,param->cosmology);
	T unredshiftedM = M/(1.+Z);
	T phase_mod = param->betappe[0];
	T generic_mod = param->betappe[1];
	return 16.*M_PI*phase_mod/(pow_int(unredshiftedM,4)) * generic_mod;
} 
template adouble EdGB_GHO_betav1(source_parameters<adouble> *);
template double EdGB_GHO_betav1(source_parameters<double> *);

template<class T>
T EdGB_GHO_betav2( source_parameters<T> *param)
{
 	T M = param->M;	
	T DL = param->DL;
	T Z= Z_from_DL(DL/MPC_SEC,param->cosmology);
	T unredshiftedM = M/(1.+Z);
	T phase_mod = param->betappe[0];
	T generic_mod = param->betappe[1];
	T phase_factor = EdGB_phase_factor(param);
	return 16.*M_PI*phase_mod/(pow_int(unredshiftedM,4)) *phase_factor* generic_mod;
} 
template adouble EdGB_GHO_betav2(source_parameters<adouble> *);
template double EdGB_GHO_betav2(source_parameters<double> *);

template<class T>
T EdGB_GHO_betav3( source_parameters<T> *param)
{
 	T M = param->M;	
	T DL = param->DL;
	T Z= Z_from_DL(DL/MPC_SEC,param->cosmology);
	T unredshiftedM = M/(1.+Z);
	T phase_mod = param->betappe[0];
	T generic_mod = param->betappe[1];
	T phase_factor = EdGB_phase_factor(param);
	T phase_GR_1PN = (3715./756. + 55. * param->eta / 9.);
	return 16.*M_PI*phase_mod/(pow_int(unredshiftedM,4)) * generic_mod* phase_factor*phase_GR_1PN;
} 
template adouble EdGB_GHO_betav3(source_parameters<adouble> *);
template double EdGB_GHO_betav3(source_parameters<double> *);

template<class T>
T ExtraDimension_beta( source_parameters<T> *param)
{
	T ED_length_sq = param->betappe[0]; //Length in seconds
	T DL = param->DL;
	T Z= Z_from_DL(DL/MPC_SEC,param->cosmology);
	T ten_mircometer = 10.e-6 / c; //10 micrometers in seconds
	T m1dot = -2.8e-7 * pow_int(MSOL_SEC * (1+Z)/ param->mass1,2)* ED_length_sq /pow_int( ten_mircometer, 2) *MSOL_SEC/T_year;
	T m2dot = -2.8e-7 * pow_int(MSOL_SEC* (1+Z)/ param->mass2,2)* ED_length_sq/pow_int( ten_mircometer, 2) *MSOL_SEC/T_year;
	T beta = (m1dot + m2dot) * (25. / 851968.) * 
		( ( 3. - 26.*param->eta + 34. * param->eta*param->eta) / 
		( pow(param->eta, 2./5.) * ( 1- 2*param->eta)) );

	return beta;
} 
template adouble ExtraDimension_beta(source_parameters<adouble> *);
template double ExtraDimension_beta(source_parameters<double> *);
	
template<class T>
T BHEvaporation_beta( source_parameters<T> *param)
{
	T mdot = param->betappe[0]; //evaporation rate -- dimensionless
	T beta = (mdot) * (25. / 851968.) * 
		( ( 3. - 26.*param->eta + 34. * param->eta*param->eta) / 
		( pow(param->eta, 2./5.) * ( 1- 2*param->eta)) );

	return beta;
} 
template adouble BHEvaporation_beta(source_parameters<adouble> *);
template double BHEvaporation_beta(source_parameters<double> *);

template<class T>
T TVG_beta( source_parameters<T> *param)
{
	T DL = param->DL;
	T Z= Z_from_DL(DL/MPC_SEC,param->cosmology);
	T Gdot = param->betappe[0];
	T beta = ( -25. / 65526. ) * ( Gdot * param->chirpmass/(1+Z)) ;
	return beta;
} 
template adouble TVG_beta(source_parameters<adouble> *);
template double TVG_beta(source_parameters<double> *);

template<class T>
T DipRad_beta( source_parameters<T> *param)
{
	T deltaEdot =param->betappe[0];
	T beta = ( -3. / 224. ) * pow(param->eta,2./5.) * deltaEdot;
	return beta;
} 
template adouble DipRad_beta(source_parameters<adouble> *);
template double DipRad_beta(source_parameters<double> *);

template<class T>
T NonComm_beta( source_parameters<T> *param)
{
	T Lambda_sq = param->betappe[0];
	T beta = ( -75. / 256. ) * pow(param->eta,-4./5.) * ( 2.*param->eta - 1.) * Lambda_sq;
	return beta;
} 
template adouble NonComm_beta(source_parameters<adouble> *);
template double NonComm_beta(source_parameters<double> *);

//#############################################################
/* \brief Modified dispersion beta
 *
 * https://arxiv.org/abs/1110.2720
 *
 * b = (alpha*3 - 3)
 *
 * Supported alpha and corresponding b and PN order:
 *
 * 	0 	-3	1
 *
 * 	0.5 	-1.5	1.75
 *
 * 	1.0 	0 	2.5
 *
 * 	1.5 	1.5 	3.25
 *
 * 	2.0 	3	4
 *
 * 	2.5 	4.5	4.75
 *
 * 	3.0 	6	5.5
 *
 * 	3.5 	7.5	6.25
 *
 * 	4.0 	9	7
 */
template<class T>
T ModDispersion_beta( source_parameters<T> *param)
{
	double alpha = (param->bppe[0]+3.)/3.;
	T Z= Z_from_DL(param->DL/MPC_SEC,param->cosmology);
	T Dalpha = DL_from_Z_MD(Z,alpha)*MPC_SEC;
	//T lambdaA = h_planck* pow(param->betappe[0],1./(alpha-2));	
	//T beta = ( pow(M_PI,2.-alpha) / (1.-alpha) ) 
	//	* ( Dalpha / pow(param->betappe[0],2.-alpha) ) 
	//	* ( pow(param->chirpmass,1.-alpha) / pow(1+Z, 1.-alpha) ) ;
	T beta = ( pow(M_PI,2.-alpha) / (1.-alpha) ) 
		* ( Dalpha * param->betappe[0]/pow(h_planck, 2.-alpha) ) 
		* ( pow(param->chirpmass,1.-alpha) / pow(1+Z, 1.-alpha) ) ;
	return beta;
} 
template adouble ModDispersion_beta(source_parameters<adouble> *);
template double ModDispersion_beta(source_parameters<double> *);

/*! \brief Calculates the modified distance given the redshift for modified dispersion relationships
 *
 * Based on Astropy.cosmology calculations -- see python script in the ./data folder of the project -- numerically calculated given astropy.cosmology's definitions (http://docs.astropy.org/en/stable/cosmology/) and used scipy.optimize to fit to a power series, stepping in half powers of Z. These coefficients are then output to a header file (D_Z_config_modified_dispersion.h) which are used here to calculate distance. Custom cosmologies etc can easily be acheived by editing the python script D_Z_config.py, the c++ functions do not need modification. They use whatever data is available in the header file. If the functional form of the fitting function changes, these functions DO need to change.
 *
 * Only Planck15 is implemented.
 *
 * Supported alphas include 0,.5,1,1.5,2,2.5,3,3.5,4
 *
 * See https://arxiv.org/abs/1110.2720 for full definitions
 * 
 */
template <class T>
T DL_from_Z_MD(T Z, double alpha)
{
	int dispersion_index = dispersion_lookup(alpha);
	if (dispersion_index == -1){ std::cout<<"Invalid Dispersion power"<<std::endl;return -1;}
	const double *boundaries = MD_boundaries_Z[dispersion_index];
	int interp_deg = interp_MD_degree[dispersion_index];
	int num_seg = num_MD_segments[dispersion_index];
	T dl;
	for (int i =0; i<num_seg; i++){
		if ( Z<boundaries[i+1]){
			//if(i == 0 ){
			//	debugger_print(__FILE__,__LINE__,"Boundaries: "+std::to_string(boundaries[i])+" "+std::to_string(boundaries[i+1]) +" "+std::to_string(Z));
			//}
			double *coeffs = new double [interp_deg];
			for (int j =0; j<interp_deg;j++)
				coeffs[j]=MD_COEFF_VEC_ZD[dispersion_index][i][j];
			dl =  cosmology_interpolation_function_MD(Z,coeffs, interp_deg);
			delete[] coeffs;
			return dl;
			
		}	
	}
	return -1;
}
template adouble DL_from_Z_MD(adouble, double);
template double DL_from_Z_MD(double, double);

template <class T>
T cosmology_interpolation_function_MD(T x, double *coeffs, int interp_degree)
{
	T result=0;
	for(int i = 0 ; i<interp_degree; i++){
		result+= coeffs[i] * pow(x,-3.5 + i*0.5);
	}
	return result;

}
template adouble cosmology_interpolation_function_MD(adouble, double*,int);
template double cosmology_interpolation_function_MD(double, double*,int);
/*! \brief Helper function for mapping cosmology name to an internal index
 */
int dispersion_lookup(double alpha)
{
	for (int i =0; i<num_MD_alphas; i++){
		if (fabs(alpha- MD_alphas[i]) < DOUBLE_COMP_THRESH){
			return i;
		}
	}
	return -1;
}

/*
template<class T>
void pre_calculate_EA_factors(source_parameters<T> *p)
{
  p->ca_EA = p->betappe[0];
  p->ctheta_EA = p->betappe[1];
  p->cw_EA = p->betappe[2];
  p->csigma_EA = p->betappe[3];
  p->s1_EA = 2e-5;
  p->s2_EA = 1e-5;
  p->V_x_EA = 0;
  p->V_y_EA = 0;
  p->V_z_EA = 0;

  //Transforming to the parameters used in arXiv:1911.10278v2 (because that is where these formulas come from)
  p->c1_EA = (p->cw_EA + p->csigma_EA)/2.;
  p->c2_EA = (p->ctheta_EA - p->csigma_EA)/3.;
  p->c3_EA = (p->csigma_EA - p->cw_EA)/2.;
  p->c4_EA = p->ca_EA - (p->csigma_EA + p->cw_EA)/2.; 
  
  //more convenient parameters
  p->c13_EA = p->c1_EA + p->c3_EA;
  p->cminus_EA = p->c1_EA - p->c3_EA;
  p->c14_EA = p->c1_EA + p->c4_EA;

  //squared speeds of the different polarizations
  p->cTsq_EA = 1./(1. - p->c13_EA);
  p->cVsq_EA = (2.*p->c1_EA - p->c13_EA*p->cminus_EA)/(2.*(1.- p->c13_EA)*p->c14_EA);
  p->cSsq_EA = ((2. - p->c14_EA)*(p->c13_EA + p->c2_EA))/((2.+3.*p->c2_EA + p->c13_EA)*(1. - p->c13_EA)*p->c14_EA);

  //speeds of the different polarizations
  p->cT_EA = sqrt(p->cTsq_EA);
  p->cV_EA = sqrt(p->cVsq_EA);
  p->cS_EA = sqrt(p->cSsq_EA);

  //Relevant combinations of parameters
  p->alpha1_EA = -8.*(p->c1_EA*p->c14_EA - p->cminus_EA*p->c13_EA)/(2.*p->c1_EA - p->cminus_EA*p->c13_EA);
  p->alpha2_EA = (1./2.)*p->alpha1_EA + ((p->c14_EA - 2.*p->c13_EA)*(3.*p->c2_EA + p->c13_EA + p->c14_EA))/((p->c2_EA + p->c13_EA)*(2. - p->c14_EA));
  p->beta1_EA = -2.* p->c13_EA / p->cV_EA; 
  p->beta2_EA = (p->c14_EA - 2.* p->c13_EA)/(2.*p->c14_EA * (1 - p->c13_EA) * p->cS_EA * p->cS_EA); 
  
  p->Z_EA = ((p->alpha1_EA - 2.*p->alpha2_EA)*(1. - p->c13_EA)) / (3.*(2.*p->c13_EA - p->c14_EA));
  
  p->A1_EA = (1./p->cT_EA) + (2*p->c14_EA*p->c13_EA*p->c13_EA)/((2.*p->c1_EA - p->c13_EA*p->cminus_EA)*(2.*p->c1_EA - p->c13_EA*p->cminus_EA)*p->cV_EA) + (3.*p->c14_EA*(p->Z_EA - 1.)*(p->Z_EA - 1.))/(2.*(2. - p->c14_EA)*p->cS_EA);
  
  p->A2_EA = -(2.*p->c13_EA)/((2.*p->c1_EA - p->c13_EA*p->cminus_EA)*pow(p->cV_EA, 3.)) - 2.*(p->Z_EA - 1.)/((2. - p->c14_EA)*pow(p->cS_EA, 3.));
  
  p->A3_EA = 1./(2.*p->c14_EA* pow(p->cV_EA, 5.)) + 2./(3.*p->c14_EA * (2. - p->c14_EA)*pow(p->cS_EA, 5.));
 
  p->B3_EA = 1./(9.*p->c14_EA*(2. - p->c14_EA)*pow(p->cS_EA, 5.));
  
  p->C_EA = 4./(3.*p->c14_EA * pow(p->cV_EA, 3.)) + 4./(3.*p->c14_EA*(2. - p->c14_EA)*pow(p->cS_EA, 3.));
  
  p->D_EA = 1./(6.*p->c14_EA * pow(p->cV_EA, 5.));

  /*
  //Sensitivities
  p->compact1 = 1.; //FIX
  p->compact2 = 1.; //FIX
  
  p->O_m1_EA = (-5./7.)*p->compact1 - ((18275.*p->alpha1_EA)/168168.)*pow(p->compact1, 3.);
  p->O_m2_EA = (-5./7.)*p->compact2 - ((18275.*p->alpha1_EA)/168168.)*pow(p->compact2, 3.);
  
  p->s1_EA = ((3.*p->alpha1_EA + 2.*p->alpha2_EA)/3.) * (p->O_m1_EA)
    +((573.*pow(p->alpha1_EA, 3.) + p->alpha1_EA*p->alpha1_EA*(67669. - 764.*p->alpha2_EA) + 96416.*p->alpha2_EA*p->alpha2_EA + 68.*p->alpha1_EA*p->alpha2_EA*(9.*p->alpha2_EA - 2632.))/(25740.*p->alpha1_EA)) * (p->O_m1_EA*p->O_m1_EA)
    + (1./(656370000.*p->cw_EA*p->alpha1_EA*p->alpha1_EA))*(-4.*p->alpha1_EA*p->alpha1_EA*(p->alpha1_EA + 8.)*(36773030.*p->alpha1_EA*p->alpha1_EA - 39543679.*p->alpha1_EA*p->alpha2_EA + 11403314.*p->alpha2_EA*p->alpha2_EA) + p->cw_EA*(1970100.*pow(p->alpha1_EA,5.) - 13995878400.*pow(p->alpha2_EA, 3.) - 640.*p->alpha1_EA*p->alpha2_EA*p->alpha2_EA*(-49528371. + 345040.*p->alpha2_EA) - 5.*pow(p->alpha1_EA, 4.)*(19548109. + 788040.*p->alpha2_EA) - 16.*p->alpha1_EA*p->alpha1_EA*p->alpha2_EA*(1294533212. - 29152855.*p->alpha2_EA + 212350.*p->alpha2_EA*p->alpha2_EA) + pow(p->alpha1_EA,3.)*(2699192440. - 309701434.*p->alpha2_EA + 5974000.*p->alpha2_EA*p->alpha2_EA))) * (pow(p->O_m1_EA, 3.));

    p->s2_EA = ((3.*p->alpha1_EA + 2.*p->alpha2_EA)/3.) * (p->O_m2_EA)
    +((573.*pow(p->alpha1_EA, 3.) + p->alpha1_EA*p->alpha1_EA*(67669. - 764.*p->alpha2_EA) + 96416.*p->alpha2_EA*p->alpha2_EA + 68.*p->alpha1_EA*p->alpha2_EA*(9.*p->alpha2_EA - 2632.))/(25740.*p->alpha1_EA)) * (p->O_m2_EA*p->O_m2_EA)
    + (1./(656370000.*p->cw_EA*p->alpha1_EA*p->alpha1_EA))*(-4.*p->alpha1_EA*p->alpha1_EA*(p->alpha1_EA + 8.)*(36773030.*p->alpha1_EA*p->alpha1_EA - 39543679.*p->alpha1_EA*p->alpha2_EA + 11403314.*p->alpha2_EA*p->alpha2_EA) + p->cw_EA*(1970100.*pow(p->alpha1_EA,5.) - 13995878400.*pow(p->alpha2_EA, 3.) - 640.*p->alpha1_EA*p->alpha2_EA*p->alpha2_EA*(-49528371. + 345040.*p->alpha2_EA) - 5.*pow(p->alpha1_EA, 4.)*(19548109. + 788040.*p->alpha2_EA) - 16.*p->alpha1_EA*p->alpha1_EA*p->alpha2_EA*(1294533212. - 29152855.*p->alpha2_EA + 212350.*p->alpha2_EA*p->alpha2_EA) + pow(p->alpha1_EA,3.)*(2699192440. - 309701434.*p->alpha2_EA + 5974000.*p->alpha2_EA*p->alpha2_EA))) * (pow(p->O_m2_EA, 3.));
  */
/*
  //The functions that are actually used to compute the phase
  p->S_EA = p->s1_EA*(p->mass2/p->M) + p->s2_EA*(p->mass1/p->M); 
  p->kappa3_EA = p->A1_EA + p->S_EA * p->A2_EA + p->S_EA*p->S_EA * p->A3_EA;
  p->epsilon_x_EA = (((p->s1_EA - p->s2_EA)*(p->s1_EA - p->s2_EA))/(32.*p->kappa3_EA))*((21.*p->A3_EA + 90.*p->B3_EA + 5.*p->D_EA)*(p->V_x_EA*p->V_x_EA + p->V_y_EA*p->V_y_EA + p->V_z_EA*p->V_z_EA) - (3.*p->A3_EA + 90.*p->B3_EA - 5.*p->D_EA)*p->V_z_EA*p->V_z_EA + 5.*p->C_EA);

  //Functions necessary for corrections to the amplitude
  p->alpha_ppE_2T_0_EA = -(1./2.)*(1./sqrt(p->kappa3_EA)) * pow(p->eta, 2./5.) * p->epsilon_x_EA;
  p->abL_EA = 1. + 2*p->beta2_EA; 
  p->gb1_EA = (2./(2. - p->c14_EA)) * (-3. *p->c14_EA * (p->Z_EA - 1) * p->cS_EA * p->cS_EA + 2.*p->S_EA)/(p->cS_EA * p->cS_EA);
  p->gX1_EA = - (p->beta1_EA)/(2*p->c1_EA - p->c13_EA*p->cminus_EA) * (1./p->cV_EA) * (p->S_EA - p->c13_EA/(1 - p->c13_EA)); 
  //debugger_print(__FILE__,__LINE__,"EA Debugging");
  //std::cout<<"aBL "<<p->abL_EA<<std::endl;
  //std::cout<<"gb1 "<<p->gb1_EA<<std::endl;
  //std::cout<<"gX1 "<<p->gX1_EA<<std::endl;
  //std::cout<<"epsilon_x "<<p->epsilon_x_EA<<std::endl;
  //std::cout<<"S "<<p->S_EA<<std::endl;
  //std::cout<<"alpha "<<p->alpha_ppE_2T_0_EA<<std::endl;
  //std::cout<<"k3 "<<p->kappa3_EA<<std::endl;
  //std::cout<<"cT "<<p->cT_EA<<std::endl;
  //std::cout<<"cV "<<p->cV_EA<<std::endl;
  //std::cout<<"cS "<<p->cS_EA<<std::endl;
  
}
template void pre_calculate_EA_factors(source_parameters<double> *);
template void pre_calculate_EA_factors(source_parameters<adouble> *);
*/
//#############################################################
