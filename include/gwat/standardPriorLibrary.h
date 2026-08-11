#ifndef STANDARDPRIORLIBRARY_H
#define STANDARDPRIORLIBRARY_H

#include "util.h"
#include <bayesship/bayesshipSampler.h>
#include <bayesship/dataUtilities.h>

struct priorData {
  double** mod_priors;
  double mass1_prior[2];
  double mass2_prior[2];
  double nbc1_prior[2];
  double nbc2_prior[2];
  double spin1_prior[2];
  double spin2_prior[2];
  double tidal1_prior[2];
  double tidal2_prior[2];
  double tidal_s_prior[2];
  double diss_tidal1_prior[2];
  double diss_tidal2_prior[2];
  double EA_prior[6];
  double RA_bounds[2];
  double sinDEC_bounds[2];
  bool tidal_love = true;
  bool tidal_love_error;
  double EOS_prior[8];
  bool EOS_plat_flag;
  bool alpha_param;
  double DL_prior[2];
  double T_mcmc_gw_tool;
  double a1_prior[2];
  double a2_prior[2];
  double ctheta1_prior[2];
  double ctheta2_prior[2];
  double phi1_prior[2];
  double phi2_prior[2];
  double ciota_bounds[2]{-1., 1.};
  double PSI_bounds[2]{0, M_PI};
  double phi_ref_bounds[2]{0, 2 * M_PI};
  double T_merger;
};

bool tidal_love_boundary_violation(double q,double lambda_s);

double EA_current_constraints(bayesship::positionInfo *position, int dim, priorData *PD);

//Uniform in m1 and m2, transformed to lnM and eta
double chirpmass_eta_jac(double chirpmass, double eta);

//Uniform in m1 and m2, transformed to lnM and q
double chirpmass_q_jac(double chirpmass, double q);

//For use with spin aligned system 
//-- mimics uniform magnitude and uniform cos \tilt, but mapped to chi
//-- Fit in mathematica 
double aligned_spin_prior(double chi);


class logPriorStandard_D: public bayesship::probabilityFn
{
public:
	priorData *PD=nullptr;
	logPriorStandard_D(priorData *PD){this->PD = PD;};
	virtual double eval(bayesship::positionInfo *position, int chainID);
};

class logPriorStandard_P: public bayesship::probabilityFn
{
public:
	priorData *PD=nullptr;
	logPriorStandard_P(priorData *PD){this->PD = PD;};
	virtual double eval(bayesship::positionInfo *position, int chainID);
};

class logPriorStandard_D_NRT: public logPriorStandard_D
{
public:
	logPriorStandard_D_NRT(priorData *PD): logPriorStandard_D(PD){this->PD = PD;};
	virtual double eval(bayesship::positionInfo *position, int chainID);
};

class logPriorStandard_D_NRT_EA: public logPriorStandard_D_NRT
{
public:
	logPriorStandard_D_NRT_EA(priorData *PD): logPriorStandard_D_NRT(PD){this->PD = PD;};
	virtual double eval(bayesship::positionInfo *position, int chainID);
};

class logPriorStandard_D_mod: public logPriorStandard_D
{
public:
	logPriorStandard_D_mod(priorData *PD): logPriorStandard_D(PD){this->PD = PD;};
	virtual double eval(bayesship::positionInfo *position, int chainID);
};

class logPriorStandard_P_mod: public logPriorStandard_P
{
public:
	logPriorStandard_P_mod(priorData *PD): logPriorStandard_P(PD){this->PD = PD;};
	virtual double eval(bayesship::positionInfo *position, int chainID);
};

class logPriorStandard_D_NRT_mod: public logPriorStandard_D_NRT
{
public:
	logPriorStandard_D_NRT_mod(priorData *PD): logPriorStandard_D_NRT(PD){this->PD = PD;};
	virtual double eval(bayesship::positionInfo *position, int chainID);
};

class logPriorStandard_P_NRT: public logPriorStandard_P
{
public:
	logPriorStandard_P_NRT(priorData *PD): logPriorStandard_P(PD){this->PD = PD;};
	virtual double eval(bayesship::positionInfo *position, int chainID);
};

class logPriorStandard_P_NRT_mod: public logPriorStandard_P_NRT
{
public:
	logPriorStandard_P_NRT_mod(priorData *PD): logPriorStandard_P_NRT(PD){this->PD = PD;};
	virtual double eval(bayesship::positionInfo *position, int chainID);
};

class logPriorStandard_D_NRT_EOS: public bayesship::probabilityFn
{
public:
	priorData *PD=nullptr;
	logPriorStandard_D_NRT_EOS(priorData *PD){this->PD = PD;};
	virtual double eval(bayesship::positionInfo *position, int chainID);
};

#endif
