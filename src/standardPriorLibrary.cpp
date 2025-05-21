#include "standardPriorLibrary.h"
#include "EA_IMRPhenomD_NRT.h"
#include <bayesship/bayesshipSampler.h>
#include <bayesship/dataUtilities.h>
#include "waveform_generator.h"
#include "ppE_utilities.h"



double chirpmass_eta_jac(double chirpmass, double eta){
  double epsilon = 1e-12;
  double delta = sqrt(1.-4.*eta);
  if(eta >.25-epsilon)
    {
      delta = sqrt(1.-4.*(eta - epsilon));
    }
	
  return chirpmass*chirpmass/(delta*pow(eta,1.2));
}

double chirpmass_q_jac(double chirpmass, double q){
	return chirpmass*chirpmass/(pow(q/pow_int(q+1,2),1./5.) * q);
}

double aligned_spin_prior(double chi){
	double a=0.0039132 , b= 3.95381;
	return a * exp(-b * abs(chi));
}

bool tidal_love_boundary_violation(double q,double lambda_s)
{
	//Relation from 1903.03909v7 fit as rough threshhold for the
	//validity of the tidal_s-tidal_a-q relation
	//Fit in log(lambda_s) - q space
	if(  q< 1.2321 - .124616*log(lambda_s)){return true;}

	return false;
}

double EA_current_constraints(bayesship::positionInfo *position,  priorData *PD)
{
  //std::cout<<"Checking EA_constraints"<<std::endl;

	double * pos = position->parameters;
	int dim = position->dimension;
	double a = -std::numeric_limits<double>::infinity();

  source_parameters<double> sp;
  double lnChirpmass = pos[7];//ln(M_sol)
  double eta = pos[8];

  sp.mass1 = calculate_mass1(std::exp(lnChirpmass)*MSOL_SEC,eta);
  sp.mass2 = calculate_mass2(std::exp(lnChirpmass)*MSOL_SEC,eta);
  sp.M = sp.mass1 + sp.mass2;
  sp.alpha_param = PD->alpha_param;

  //Setting tidal1 and tidal2 from tidal_s
  if(PD->tidal_love)
    {
      IMRPhenomD_NRT<double> modelNRT;
      modelNRT.binary_love_relation(pos[11], PD->tidal_love_error, &sp);
    }
  sp.csigma_EA = 0;

  if(PD->tidal_love){
    if(PD->alpha_param){
      sp.alpha1_EA = pos[12]; //alpha1
      sp.alpha2_EA = pos[13]; //alpha2
      sp.cbarw_EA = pos[14]; //cbarw
    }
    else{
      sp.ca_EA = pos[12]; //ca
      sp.ctheta_EA = pos[13]; //ctheta
      sp.cw_EA = pos[14]; //cw
    }
  }
  else{
    if(PD->alpha_param){
      sp.alpha1_EA = pos[13]; //alpha1
      sp.alpha2_EA = pos[14]; //alpha2
      sp.cbarw_EA = pos[15]; //cbarw
    }
    else{
      sp.ca_EA = pos[13]; //ca
      sp.ctheta_EA = pos[14]; //ctheta
      sp.cw_EA = pos[15]; //cw
    }
  }
  sp.EA_nan_error_message = false;
  EA_IMRPhenomD_NRT<double> model;
  model.pre_calculate_EA_factors(&sp);

  if(sp.ca_EA < 0 || sp.ca_EA > 2.){return a;}
  /* Throws out points with ca < 0 or ca > 2 because these violate
   * the positive energy condition for the spin-0 mode (scalar mode).
   * See equation 40 of arXiv:gr-qc/0507059v3.
   */
  if(sp.cw_EA < (-sp.csigma_EA/(1. - sp.csigma_EA))){return a;}
  /* Throws out points with cw < -csigma/(1 - csigma) because these violate
   * the positive energy condition for the spin-1 mode (vector mode).
   * Note that the positive energy condition for the spin-2 modes is always
   * satisfied. See equation 40 of arXiv:gr-qc/0507059v3.
   */
  //std::cout<<"EA constraints test 1"<<std::endl;
  if(sp.cTsq_EA < 0 || sp.cVsq_EA < 0 || sp.cSsq_EA < 0){return a;}
  /* Throws out points with speeds not greater than or equal to zero (these
   * would produce gradient instabilities or ghosts)
   * arXiv:gr-qc/0402005 and arXiv:1108.1835
   */
  //std::cout<<"EA constraints test 2"<<std::endl;

  if(isnan(sp.kappa3_EA))
    {
      //std::cout<<"kappa3:"<<sp.kappa3_EA<<std::endl;
      return a;
      }
  //std::cout<<"EA constraints test 3"<<std::endl;

  if(sp.ctheta_EA < 0 || (sp.ctheta_EA + (8.*sp.ca_EA/7)) > (2./7.)){return a;}
  /*
   * Throws out points that violate Big Bang Nucleosynthesis constraints. arXiv:hep-th/0407149v3
   */

  if(sp.s1_EA >=1 || sp.s2_EA >=1){return a;}

  //if(fabs(sp.alpha1_EA) > pow(10, -4.) || fabs(sp.alpha2_EA) > 4.*pow(10, -7.)){return a;}
  /* Throws out points that do not obey observational solar system constraints on
   * alpha1 and alpha2
   * arXiv:1403.7377 and arXiv:gr-qc/0509114
   */
  //std::cout<<"EA constraints test 5"<<std::endl;

  bool Cherenkov = true;
  if(Cherenkov)
  {
  	bool violate = false;
  	if(sp.cV_EA < 1)
  	  {
      		if(abs(sp.c13_EA * sp.c13_EA *(sp.c1_EA * sp.c1_EA + 2*sp.c1_EA*sp.c3_EA + sp.c3_EA * sp.c3_EA - 2*sp.c4_EA)/(2*sp.c1_EA*sp.c1_EA)) >= 7*pow(10, -32.))
		//enforcing constraint from Eq. 4.7 of arXiv:hep-ph/0505211
		{
	  	  violate = true;
		}

	    }
	  if(violate){return a;}
	  //std::cout<<"EA constraints test 6"<<std::endl;

	  if(sp.cS_EA < 1)
	    {
	      if(fabs((sp.c2_EA + sp.c3_EA - sp.c4_EA)/sp.c1_EA) > pow(10, -22.))
		{
		  if((sp.c3_EA - sp.c4_EA)*(sp.c3_EA - sp.c4_EA)/fabs(sp.c14_EA) >= pow(10, -30.)){return a;}
		  //enforcing constraint from Eq.4.15 of arXiv:hep-ph/0505211
		}
	      //std::cout<<"EA constraints test 7"<<std::endl;
	      if(fabs((sp.c4_EA - sp.c2_EA - sp.c3_EA)/sp.c1_EA) >= 3*pow(10,-19.)){return a;}
	      //enforcing constraint from Eq.5.14 of arXiv:hep-ph/0505211

	    }
  }

  //std::cout<<"Made it to end of EA constraints"<<std::endl;
  /*
  //Enforcing gaussian prior on alpha1 from binary pulsar and triple systems.
  //arXiv:2104.04596
  double sigma = 1.021*pow(10, -5.);
  double mu = -0.563*pow(10, -5.);
  double prob;

  prob = exp(-(1./2.)*((sp.alpha1_EA - mu)*(sp.alpha1_EA - mu))/(sigma*sigma));
  */

  sp.EA_nan_error_message = true;
  model.EA_check_nan(&sp);

  //return log(prob); //Use if enforcing gaussian on alpha1
  return 0;
}

double logPriorStandard_D_mod::eval(bayesship::positionInfo *position, int chainID)
{
	int dim =  position->dimension;
	double *pos =  position->parameters;
	double a = -std::numeric_limits<double>::infinity();
	int initial_nongr_id = 11;
	for (int i = initial_nongr_id ; i<dim; i++){
		if(pos[i]<PD->mod_priors[i-initial_nongr_id][0] || pos[i]>PD->mod_priors[i-initial_nongr_id][1]){return a;}
	}

	double PhenomD = logPriorStandard_D::eval(position,chainID);
	return PhenomD;
}

double logPriorStandard_D_NRT_mod::eval(bayesship::positionInfo *position, int chainID)
{
	int dim =  position->dimension;
	double *pos =  position->parameters;
	double a = -std::numeric_limits<double>::infinity();
	int initial_nongr_id = 12;
	if(! PD->tidal_love){
		initial_nongr_id = 13;
	}
	for (int i = initial_nongr_id ; i<dim; i++){
		if(pos[i]<PD->mod_priors[i-initial_nongr_id][0] || pos[i]>PD->mod_priors[i-initial_nongr_id][1]){return a;}
	}

	double NS = logPriorStandard_D_NRT::eval(position,chainID);
	//std::cout<<NS<<std::endl;
	return  NS;
}

double logPriorStandard_D_NRT_EA::eval(bayesship::positionInfo *position, int chainID)
{
	int dim =  position->dimension;
	double *pos =  position->parameters;
	double a = -std::numeric_limits<double>::infinity();
	if(PD->tidal_love){
	  if(pos[12]<PD->EA_prior[0] || pos[12]>PD->EA_prior[1]){return a;} //ca or alpha1
	  if(pos[13]<PD->EA_prior[2] || pos[13]>PD->EA_prior[3]){return a;} //ctheta or alpha2
	  if(pos[14]<PD->EA_prior[4] || pos[14]>PD->EA_prior[5]){return a;} //cw or cbarw
	}
	else{
	  if(pos[13]<PD->EA_prior[0] || pos[13]>PD->EA_prior[1]){return a;} //ca or alpha1
	  if(pos[14]<PD->EA_prior[2] || pos[14]>PD->EA_prior[3]){return a;} //ctheta or alpha2
	  if(pos[15]<PD->EA_prior[4] || pos[15]>PD->EA_prior[5]){return a;} //cw or cbarw
	}

	double NS = logPriorStandard_D_NRT::eval(position,chainID);
	if(NS == a){return a;}
	double EA_constraints =  EA_current_constraints(position, PD);
	if (EA_constraints ==a){return a;}

	return EA_constraints + NS;
}

double logPriorStandard_D_NRT::eval(bayesship::positionInfo *position, int chainID)
{
	int dim =  position->dimension;
	double *pos = position->parameters;
	double a = -std::numeric_limits<double>::infinity();
	double chirp = exp(pos[7]);
	double m1 = calculate_mass1(chirp,pos[8]);
	double m2 = calculate_mass2(chirp,pos[8]);
	double q = m2/m1;//<1
	double factor = 0;
	if(PD->tidal_love){
		if(exp(pos[11])<PD->tidal_s_prior[0] || exp(pos[11])>PD->tidal_s_prior[1]){return a;}
		if(tidal_love_boundary_violation(q,exp(pos[11]))){return a;}
		factor += pos[11];

	}
	else{
		if(exp(pos[11])<PD->tidal1_prior[0] || exp(pos[11])>PD->tidal1_prior[1]){return a;}
		if(exp(pos[12])<PD->tidal2_prior[0] || exp(pos[12])>PD->tidal2_prior[1]){return a;}
		factor += pos[11];
		factor += pos[12];
	}
	return logPriorStandard_D::eval(position, chainID) +factor ;

}

double logPriorStandard_D::eval(bayesship::positionInfo *position, int chainID)
{
	int dim =  position->dimension;

	double a = -std::numeric_limits<double>::infinity();
	double *pos = position->parameters;
	//###########
	double chirp = exp(pos[7]);
	double eta = pos[8];
	if (eta<.0 || eta>.25){return a;}//eta
	double m1 = calculate_mass1(chirp,eta );
	double m2 = calculate_mass2(chirp,eta );
	if(m1<PD->mass1_prior[0] || m1>PD->mass1_prior[1]){return a;}
	if(m2<PD->mass2_prior[0] || m2>PD->mass2_prior[1]){return a;}
	//###########
	if ((pos[0])<PD->RA_bounds[0] || (pos[0])>PD->RA_bounds[1]){ return a;}//RA
	if ((pos[1])<PD->sinDEC_bounds[0] || (pos[1])>PD->sinDEC_bounds[1]){return a;}//sinDEC
	if ((pos[2])<0 || (pos[2])>M_PI){return a;}//PSI
	if ((pos[3])<-1 || (pos[3])>1){return a;}//cos \iota
	if ((pos[4])<0 || (pos[4])>2*M_PI){return a;}//phiRef
	if( pos[5] < (PD->T_merger - .1) || pos[5] > (PD->T_merger + .1)) { return a; }
	if (std::exp(pos[6])<PD->DL_prior[0] || std::exp(pos[6])>PD->DL_prior[1]){return a;}//DL
	if ((pos[9])<PD->spin1_prior[0] || (pos[9])>PD->spin1_prior[1]){return a;}//chi1
	if ((pos[10])<PD->spin2_prior[0] || (pos[10])>PD->spin2_prior[1]){return a;}//chi2
	//return log(chirpmass_eta_jac(chirp,eta))+3*pos[6] ;
	return log(aligned_spin_prior(pos[9]))+log(aligned_spin_prior(pos[10])) + log(chirpmass_eta_jac(chirp,eta))+3*pos[6] ;

}

double logPriorStandard_D_NRT_EOS::eval(bayesship::positionInfo *position, int chainID)
{
	int dim =  position->dimension;

	double a = -std::numeric_limits<double>::infinity();
	double *pos = position->parameters;
	//###########
	//double chirp = exp(pos[7]);
	//double eta = pos[8];
	//if (eta<.0 || eta>.25){return a;}//eta
	//double m1 = calculate_mass1(chirp,eta );
	//double m2 = calculate_mass2(chirp,eta );
	if(pos[7]<PD->nbc1_prior[0] || pos[7]>PD->nbc1_prior[1]){return a;}
	if(pos[8]<PD->nbc2_prior[0] || pos[8]>PD->nbc2_prior[1]){return a;}
	//###########
	if ((pos[0])<PD->RA_bounds[0] || (pos[0])>PD->RA_bounds[1]){ return a;}//RA
	if ((pos[1])<PD->sinDEC_bounds[0] || (pos[1])>PD->sinDEC_bounds[1]){return a;}//sinDEC
	if ((pos[2])<0 || (pos[2])>M_PI){return a;}//PSI
	if ((pos[3])<-1 || (pos[3])>1){return a;}//cos \iota
	if ((pos[4])<0 || (pos[4])>2*M_PI){return a;}//phiRef
	if( pos[5] < (PD->T_merger - .1) || pos[5] > (PD->T_merger + .1)) { return a; }
	if (std::exp(pos[6])<PD->DL_prior[0] || std::exp(pos[6])>PD->DL_prior[1]){return a;}//DL
	if ((pos[9])<PD->spin1_prior[0] || (pos[9])>PD->spin1_prior[1]){return a;}//chi1
	if ((pos[10])<PD->spin2_prior[0] || (pos[10])>PD->spin2_prior[1]){return a;}//chi2
	//###########
	if ((pos[11])<PD->EOS_prior[0] || (pos[11])>PD->EOS_prior[1]){return a;}//bump_mag
	if ((pos[12])<PD->EOS_prior[2] || (pos[12])>PD->EOS_prior[3]){return a;}//bump_width
	if ((pos[13])<PD->EOS_prior[4] || (pos[13])>PD->EOS_prior[5]){return a;}//bump_offset
	if(PD->EOS_plat_flag){
	  if ((pos[14])<PD->EOS_prior[6] || (pos[14])>PD->EOS_prior[7]){return a;}//plat
	}
	//###########
	//return log(chirpmass_eta_jac(chirp,eta))+3*pos[6] ;
	return log(aligned_spin_prior(pos[9]))+log(aligned_spin_prior(pos[10])) +3*pos[6] ;

}
