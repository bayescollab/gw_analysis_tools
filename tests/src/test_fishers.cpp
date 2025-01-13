#include "gwat/util.h"
#include "gwat/fisher.h"
#include "gwat/detector_util.h"
#include "gwat/waveform_util.h"
#include "gwat/ortho_basis.h"
#include "gwat/pn_waveform_util.h"
#include "gwat/ppE_utilities.h"
#include <iostream>
#include <iomanip>



int AD_v_N(int argc, char *argv[]);
int network_fishers(int argc, char *argv[]);
int dCS_EdGB(int argc, char *argv[]);
int test_jac_transform(int argc, char *argv[]);
int test_MCMC_fisher(int argc, char *argv[]);
int test_LISA_fisher(int argc, char *argv[]);
int test_EA_fisher(int argc, char *argv[]);
void RT_ERROR_MSG();

int main(int argc, char *argv[])
{
	std::cout<<"TESTING FISHER CALCULATIONS"<<std::endl;
		
	if(argc != 2){
		RT_ERROR_MSG();
		return 1;
	}
	
	int runtime_opt = stoi(argv[1]);	
	if(runtime_opt == 0){
		return AD_v_N(argc,argv);
	}
	if(runtime_opt == 1){
		return network_fishers(argc,argv);
	}
	else if(runtime_opt == 2){
		return dCS_EdGB(argc,argv);
	}
	else if(runtime_opt == 3){
		return test_jac_transform(argc,argv);
	}
	else if(runtime_opt == 4){
		return test_MCMC_fisher(argc,argv);
	}
	else if(runtime_opt == 5){
		return test_LISA_fisher(argc,argv);
	}
	else if(runtime_opt == 6){
		return test_EA_fisher(argc,argv);
	}
	else{
		RT_ERROR_MSG();
		return 1;
	}
}
int test_EA_fisher(int argc, char *argv[])
{
	std::cout.precision(15);
	//Create injection structure
	gen_params params;	

	params.mass1 = 1.9;
	//params.mass1 = 1.4;
	//params.mass1 = (1.4 - pow(10, -1)); 
	params.mass2 = 1.4;

  params.spin1[2] = -.03;
	params.spin2[2] = .03 ;
	params.Luminosity_Distance = 30;
	params.incl_angle = 3*M_PI/4;

	params.NSflag1 = true;
	params.NSflag2 =true;

	params.phiRef = .0;
	params.RA = 1.;
	params.DEC = 0.6;
	params.f_ref = 20;
	double chirpmass = calculate_chirpmass(params.mass1,params.mass2);//*MSOL_SEC;
	double eta = calculate_eta(params.mass1,params.mass2);//*MSOL_SEC;

	params.horizon_coord = false;
	params.shift_time=false;
	params.shift_phase=false;
	params.dep_postmerger=true;
	params.sky_average=false;
	params.tidal_love=true;
	params.tidal_s=200;
	
	//params.tidal_love=false;
	//params.tidal1=200;
	//params.tidal2=100;
	
	params.psi = 1.;
	params.gmst = 2.;
	params.sky_average = false;


	//#######################################
	//EA parameters
	//#######################################
	params.ca_EA = 5.01510850153863e-07; 
	params.ctheta_EA = 1.50544346457309e-06; 
	params.cw_EA = 5.11795863509178; 
	//params.ca_EA = 1e-7;
	//params.ctheta_EA = 2e-7;
	//params.cw_EA = 2e-7;
	//params.csigma_EA = 1e-30;
	//#######################################
	//#######################################
	
	std::cout<<"m1: "<<params.mass1<<", m2: "<<params.mass2<<std::endl;
	std::cout<<"chirpmass: "<<chirpmass<<", eta: "<<eta<<std::endl; 

	//#######################################
	//#######################################
	//Detector characteristics
	double fmin = 5;
	double fmax = 2048;
	double T = 16;
	params.tc = 3.*T/4.;

	int length = 2000;
	double *frequency = new double[length];
	int Ndetect = 3;
	double **psd = new double*[Ndetect];
	//Noise curves to use
	std::string SN[3] = {"AdLIGOMidHigh","AdLIGOMidHigh","AdLIGOMidHigh"}; //"AdVIRGOPlus1"};
	
	//Calculate freq/weight array (using gauss-legendre quadrature)
	bool AD = true;
	bool GL = false;
	double *weights = new double[length];
	if(AD && GL){
	//if(false){
		gauleg(log10(fmin), log10(fmax),frequency,weights,length);
		for(int i = 0 ; i<length; i++){
			frequency[i] = pow(10,frequency[i]);	
		}
	}
	else{
		double deltaF = (fmax-fmin)/length;	
		for(int i = 0 ; i<length; i++){
			frequency[i] = fmin + deltaF*i;
		}
	}
	

	for(int i = 0 ; i<Ndetect; i++){
		psd[i]= new double[length];
		populate_noise(frequency, SN[i],psd[i], length, 48);
		for(int j = 0 ; j<length; j++){
			psd[i][j]*=psd[i][j];	
		}
	}
	std::cout<<"frequency[10]:"<<frequency[10]<<std::endl; 
	
	//Set detector names
	std::string detectors[3] = {"Hanford","Livingston","Virgo"};
		
	//###############################################
	//Dimension and model
	//###############################################
	//int dim = 15;
	//std::string method = "EA_IMRPhenomD_NRT";
	int dim = 11; 
	std::string method = "IMRPhenomD"; 
	//int dim = 12;
	//std::string method = "IMRPhenomD_NRT"; 
	
	//###############################################
	//Allocate output
	//###############################################
	double **output_AD = allocate_2D_array(dim,dim);
	double **output_AD_temp = allocate_2D_array(dim,dim);
	double **COV_AD = allocate_2D_array(dim,dim);
	for(int i = 0 ; i<dim; i++){
		for(int j = 0 ; j<dim; j++){
			output_AD[i][j]= 0;
			output_AD_temp[i][j]= 0;
		}
	}


	double snr; 
	double total_snr = 0;

	//###############################################
	//Calculate Fishers
	//###############################################
	for(int i = 0 ;i < Ndetect; i++){
		if(AD){
			if(GL){
				total_snr += pow_int( calculate_snr(SN[i],detectors[i],method, &params, frequency, length, "GAUSSLEG", weights, true), 2);
				fisher_autodiff(frequency, length, method, detectors[i],detectors[0], output_AD_temp, dim, &params, "GAUSSLEG",weights,true, psd[i],NULL,NULL);
				//debugger_print(__FILE__,__LINE__, total_snr);

			}
			else{
				total_snr += pow_int( calculate_snr(SN[i],detectors[i],method, &params, frequency, length, "SIMPSONS", weights, false), 2);
				fisher_autodiff(frequency, length, method, detectors[i],detectors[0], output_AD_temp, dim, &params, "SIMPSONS",weights,false, psd[i],NULL,NULL);
				//debugger_print(__FILE__,__LINE__, total_snr);
			}
		}
		else{
		        total_snr += pow_int( calculate_snr(SN[i],detectors[i],method, &params, frequency, length, "SIMPSONS", weights, false), 2);
			//debugger_print(__FILE__,__LINE__,total_snr);
			fisher_numerical(frequency, length, method, detectors[i],detectors[0], output_AD_temp, dim, &params, 2,NULL,NULL, psd[i]);
		}
		for(int k = 0 ; k<dim; k++){
			//std::cout<<i<<": "<<std::endl;
			for(int j = 0 ; j<dim; j++){
				output_AD[k][j]+= output_AD_temp[k][j];
				//std::cout<<std::setprecision(5)<<output_AD[i][j]<<" ";
			}
			//std::cout<<std::endl;
		}
	}
	std::cout<<"Total SNR: "<<sqrt(total_snr)<<std::endl;

	//####################################
	//Add prior
	//double sigma_a = 1e-8;
	//double sigma_theta = 1e-4;
	//double sigma_omega = 1e-5;
	//double sigma_sigma = 1e-15;

	//output_AD[dim-4][dim-4]+= 1./sigma_a/sigma_a;
	//output_AD[dim-3][dim-3]+= 1./sigma_theta/sigma_theta;
	//output_AD[dim-2][dim-2]+= 1./sigma_omega/sigma_omega;
	//output_AD[dim-1][dim-1]+= 1./sigma_sigma/sigma_sigma;
	//####################################

	//####Removing correlations between EA coupling constants for simplicity
	/*for(int i = 0 ; i <3; i++){
	  for(int j = 0 ; j<dim; j++){
	    if(i!=j){
	      output_AD[dim-1-i][dim-1-j] = 0;
	      output_AD[dim-1-j][dim-1-i] = 0;
	    }
	  }
	  for(int j = 0 ; j<dim; j++){
	    if(i!=j){
	      output_AD[dim-1-j][dim-1-i] = 0; 
	    }
	  }
	  }*/

	//###############################################
	//Invert and print -- uncomment to see full cov or fisher
	//###############################################
	
	std::cout<<"SNR: "<<sqrt(output_AD[6][6])<<std::endl;
	//std::cout<<"AD Fisher:"<<std::endl;
	for(int i = 0 ; i<dim; i++){
		std::cout<<i<<" ";
		for(int j = 0 ; j<dim; j++){
		  std::cout<<std::setprecision(5)<<output_AD[i][j]<<" ";
		}
		std::cout<<std::endl;
	}

	gsl_LU_matrix_invert(output_AD,COV_AD,dim);
	//gsl_cholesky_matrix_invert(output_AD,COV_AD,dim);
	std::cout<<"COV AD:"<<std::endl;
	for(int i = 0 ; i<dim; i++){
	  std::cout<<i<<" ";
	  for(int j = 0 ; j<dim; j++){
	    std::cout<<COV_AD[i][j]<<" ";
	  }
	  std::cout<<std::endl;
	}
	std::cout<<"Standard Deviations:"<<std::endl;
	for(int i = 0 ; i<dim; i++){
		std::cout<<i<<": "<<std::sqrt(COV_AD[i][i])<<"\n";
	}


	//###############################################
	//Prints cov.Fisher, which should be the identity (since cov = fisher^-1)
	//###############################################
	//double **identity_full = allocate_2D_array(dim,dim);
	//matrix_multiply(output_AD,COV_AD,identity_full,dim,dim,dim);
	//std::cout<<"IDENTITY: "<<std::endl;
	//for(int i = 0 ; i<dim; i++){
	//	for(int j = 0 ; j<dim; j++){
	//		std::cout<<identity_full[i][j]<<" ";
	//	}
	//	std::cout<<std::endl;
	//}
	//deallocate_2D_array(identity_full, dim,dim);




	deallocate_2D_array(output_AD,dim,dim);
	deallocate_2D_array(output_AD_temp,dim,dim);
	deallocate_2D_array(COV_AD,dim,dim);
	
	delete [] frequency;
	for(int i = 0 ; i<Ndetect; i++){
		delete [] psd[i];
	}
	delete [] psd;
	delete [] weights;
	return 0;
}
int test_LISA_fisher(int argc, char *argv[])
{
	std::cout.precision(15);
	gen_params params;	
	params.mass1 = 36.9e5;
	params.mass2 = 35.93e5;
	params.spin1[2] = .8;
	params.spin2[2] = .8 ;
	params.chip = .07;
	params.phip = 1.0;
	params.Luminosity_Distance = 48000;
	//params.Luminosity_Distance = DL_from_Z(10,"PLANCK15");
	//params.Luminosity_Distance = DL_from_Z(0.0019,"PLANCK15");
	//params.incl_angle = acos(.75);
	params.incl_angle = 0;
	params.phi_l = M_PI/4;
	params.theta_l = M_PI/2;

	params.NSflag1 = false;
	params.NSflag2 =false;

	params.phiRef = .0;
	params.RA = .234;
	params.DEC = asin(-.45);
	params.f_ref = 1;
	double chirpmass = calculate_chirpmass(params.mass1,params.mass2)*MSOL_SEC;
	std::cout<<"chirpmass: "<<chirpmass/MSOL_SEC<<std::endl;
	//params.spin1[2] = .38;

	params.horizon_coord = false;
	params.shift_time=false;
	params.shift_phase=false;
	params.dep_postmerger=true;
	//params.equatorial_orientation=true;
	params.equatorial_orientation=false;
	
	params.gmst = 2.;
	//params.sky_average = false;
	params.sky_average = true;
	params.Nmod = 1;
	params.betappe = new double[1];
	params.bppe = new double[1];
	params.bppe[0]=-13;
	params.betappe[0]=0;


	double fmin = 1e-5;
	double fmax = 1e-3;
	//double fmin = .006508;
	//double fmax = .0067506;
	//double fmin = .01208;
	//double fmax = 1.00;
	double T =(t_0PN(fmin,chirpmass)- t_0PN(fmax,chirpmass));
	std::cout<<"TIME: "<<T/T_year<<std::endl;

	params.tc = 3.*T/4.;
	int length = 5000;
	double *frequency = new double[length];
	double *weights = new double[length];
	double *psd = new double[length];
	gauleg(log10(fmin),log10(fmax), frequency, weights, length);
	for(int j = 0 ; j<length; j++){
		frequency[j] = pow(10,frequency[j]);	
	}
	string SN = "LISA_SADC_CONF";
	populate_noise(frequency, SN,psd, length, 48);
	for(int j = 0 ; j<length; j++){
		psd[j]*=psd[j];	
	}


	int dim = 8;
	double **output =NULL;
	output= allocate_2D_array(dim,dim);
	double **output_temp = allocate_2D_array(dim,dim);
	double **COV = allocate_2D_array(dim,dim);
	for(int i = 0 ; i<dim; i++){
		for(int j = 0 ; j<dim; j++){
			output[i][j]= 0;
			output_temp[i][j]= 0;
		}
	}
	//std::string method = "ppE_IMRPhenomPv2_Inspiral";
	std::string method = "ppE_IMRPhenomD_IMR";
	std::string detector = "LISA";
	fisher_autodiff(frequency, length, method, detector,detector, output, dim, &params, "GAUSSLEG",weights,true, psd,NULL,NULL);
	//fisher_autodiff(frequency, length, method, detector,detector, output, dim, &params, "GAUSSLEG",weights,true, psd,NULL,NULL);
	std::cout<<"SNR: "<<sqrt(output[6][6])<<std::endl;

	gsl_LU_matrix_invert(output,COV,dim);
	//for(int i = 0 ; i<dim; i++){
	//	for(int j = 0 ; j<dim; j++){
	//		std::cout<<output[i][j]<<" ";
	//	}
	//	std::cout<<std::endl;
	//}
	//for(int i = 0 ; i<dim; i++){
	//	for(int j = 0 ; j<dim; j++){
	//		std::cout<<COV[i][j]<<" ";
	//	}
	//	std::cout<<std::endl;
	//}
	for(int i = 0 ; i<dim; i++){
		std::cout<<sqrt(COV[i][i])<<std::endl;
	}

	deallocate_2D_array(output,dim,dim);
	deallocate_2D_array(output_temp,dim,dim);
	deallocate_2D_array(COV,dim,dim);
	
	delete [] frequency;
	delete [] weights;
	delete [] psd;
	return 0;


}
int test_MCMC_fisher(int argc, char *argv[])
{
	std::cout.precision(15);
	gen_params params;	
	params.mass1 = 20.9;
	params.mass2 = 1.93;
	params.spin1[2] = .0* (params.mass1+params.mass2)/params.mass1;
	params.spin2[2] = 0. ;
	params.chip = .01;
	params.phip = 1.0;
	//params.Luminosity_Distance = 290;
	params.Luminosity_Distance = 80;
	//params.incl_angle = acos(.75);
	params.incl_angle = M_PI/2;

	params.NSflag1 = true;
	params.NSflag2 =true;
	params.tidal_s = 10;
	params.tidal_love = true;

	params.phiRef = .0;
	params.RA = .234;
	params.DEC = asin(-.45);
	params.f_ref = 20;
	double chirpmass = calculate_chirpmass(params.mass1,params.mass2)*MSOL_SEC;
	std::cout<<"chirpmass: "<<chirpmass/MSOL_SEC<<std::endl;
	//params.spin1[2] = .38;

	params.horizon_coord = false;
	params.shift_time=false;
	params.shift_phase=false;
	params.dep_postmerger=true;
	params.sky_average=false;
	params.equatorial_orientation=false;
	
	params.psi = 1.;
	params.gmst = 2.;
	params.sky_average = false;
	//params.sky_average=true;
	params.Nmod = 0;


	double fmin = 5;
	double fmax = 2048;
	double T = 16;

	params.tc = 3.*T/4.;
	int length = (1024-20)*T;
	double *frequency = new double[length];
	int Ndetect = 3;
	//int Ndetect = 2;
	double **psd = new double*[Ndetect];
	//std::string SN[3] = {"AdLIGODesign_smoothed","AdLIGODesign_smoothed","AdLIGODesign_smoothed"};
	std::string SN[3] = {"AdLIGOMidHigh","AdLIGOMidHigh","AdVIRGOPlus1"};
	//std::string SN[3] = {"CE1","AdLIGOMidHigh","AdVIRGOPlus1"};
	double deltaf = 1./T;
	for(int i = 0 ; i<length; i++){
		frequency[i] = 20+i*deltaf;	
	}
	for(int i = 0 ; i<Ndetect; i++){
		psd[i]= new double[length];
		populate_noise(frequency, SN[i],psd[i], length, 48);
		//populate_noise(frequency, "LISA_CONF",psd, length, 12);
		for(int j = 0 ; j<length; j++){
			psd[i][j]*=psd[i][j];	
		}
	}


	int dim = 12;
	//int dim = 5;
	double **output =NULL;
	output= allocate_2D_array(dim,dim);
	double **output_temp = allocate_2D_array(dim,dim);
	double **COV = allocate_2D_array(dim,dim);
	for(int i = 0 ; i<dim; i++){
		for(int j = 0 ; j<dim; j++){
			output[i][j]= 0;
			output_temp[i][j]= 0;
		}
	}
	std::string method = "MCMC_IMRPhenomD_NRT";
	std::string detectors[3] = {"Hanford","Livingston","Virgo"};
	for(int i = 0 ;i < Ndetect; i++){
		//fisher_autodiff(frequency, length, method, detectors[i],detectors[0], output_temp, dim, &params, "SIMPSONS",(double*)NULL,false, psd[i],NULL,NULL);
		fisher_numerical(frequency, length, method, detectors[i],detectors[0], output_temp, dim, &params, 2,NULL,NULL, psd[i]);
		for(int k = 0 ; k<dim; k++){
			for(int j = 0 ; j<dim; j++){
				output[k][j]+= output_temp[k][j];
			}
		}
	}
	std::cout<<"SNR: "<<sqrt(output[6][6])<<std::endl;

	gsl_LU_matrix_invert(output,COV,dim);
	for(int i = 0 ; i<dim; i++){
		for(int j = 0 ; j<dim; j++){
			std::cout<<output[i][j]<<" ";
		}
		std::cout<<std::endl;
	}
	for(int i = 0 ; i<dim; i++){
		for(int j = 0 ; j<dim; j++){
			std::cout<<COV[i][j]<<" ";
		}
		std::cout<<std::endl;
	}
	for(int i = 0 ; i<dim; i++){
		std::cout<<sqrt(COV[i][i])<<std::endl;
	}

	deallocate_2D_array(output,dim,dim);
	deallocate_2D_array(output_temp,dim,dim);
	deallocate_2D_array(COV,dim,dim);
	
	delete [] frequency;
	for(int i = 0 ; i<Ndetect; i++){
		delete [] psd[i];
	}
	delete [] psd;
	return 0;


}
int test_jac_transform(int argc, char *argv[])
{
	std::cout.precision(15);
	gen_params params;	
	params.mass1 = 33.9;
	params.mass2 = 29.3;
	params.spin1[2] = .18* (params.mass1+params.mass2)/params.mass1;
	params.spin2[2] = 0.1 ;
	params.chip = .01;
	params.phip = 1.0;
	params.Luminosity_Distance = 300;
	params.incl_angle = M_PI/10.;

	params.NSflag1 = false;
	params.NSflag2 =false;

	params.phiRef = .0;
	params.RA = 1.;
	params.DEC = 0.6;
	params.f_ref = 20;
	double chirpmass = calculate_chirpmass(params.mass1,params.mass2)*MSOL_SEC;
	std::cout<<"chirpmass: "<<chirpmass/MSOL_SEC<<std::endl;
	//params.spin1[2] = .38;

	params.horizon_coord = false;
	params.shift_time=false;
	params.shift_phase=false;
	params.dep_postmerger=true;
	params.sky_average=false;
	
	params.psi = 1.;
	params.gmst = 2.;
	params.sky_average = false;
	params.Nmod = 1;
	double mod_value = 1;
	params.betappe = new double[1];
	params.bppe = new double[1];
	//params.betappe[0] =0;
	params.betappe[0] =mod_value;
	params.bppe[0] = -7;
	//params.bppe[0] = -1.5;
	//params.bppe[0] = -3;
	

	double fmin = 5;
	double fmax = 2048;
	double T = 32;

	params.tc = 3.*T/4.;
	int length = 1000;
	double *frequency = new double[length];
	int Ndetect = 3;
	//int Ndetect = 2;
	double **psd = new double*[Ndetect];
	//std::string SN[3] = {"AdLIGODesign_smoothed","AdLIGODesign_smoothed","AdLIGODesign_smoothed"};
	std::string SN[3] = {"AdLIGOMidHigh","AdLIGOMidHigh","AdVIRGOPlus1"};
	//std::string SN[3] = {"CE1","AdLIGOMidHigh","AdVIRGOPlus1"};
	
	double *weights = new double[length];
	gauleg(log10(fmin), log10(fmax),frequency,weights,length);
	for(int i = 0 ; i<length; i++){
		frequency[i] = pow(10,frequency[i]);	
	}
	for(int i = 0 ; i<Ndetect; i++){
		psd[i]= new double[length];
		populate_noise(frequency, SN[i],psd[i], length, 48);
		//populate_noise(frequency, "LISA_CONF",psd, length, 12);
		for(int j = 0 ; j<length; j++){
			psd[i][j]*=psd[i][j];	
		}
	}

	std::string method = "DipRad_IMRPhenomPv2";
	//std::string method = "ExtraDimension_IMRPhenomPv2";
	//std::string method = "TVG_IMRPhenomPv2";
	//std::string method = "ModDispersion_IMRPhenomPv2";
	std::string detectors[3] = {"Hanford","Livingston","Virgo"};

	//#########################
	source_parameters<double> sparam;
	std::string local_method = prep_source_parameters(&sparam,&params,method);
	double beta_value =  sparam.betappe[0];
	
	std::cout<<beta_value<<std::endl;;
	std::cout<<local_method<<std::endl;;
	cleanup_source_parameters(&sparam, method);
	//#########################


		
	int dim = 14;
	int dimDSA = 8;

	double **output_AD = allocate_2D_array(dim,dim);
	double **output_AD_temp = allocate_2D_array(dim,dim);
	double **output_AD2 = allocate_2D_array(dim,dim);
	double **output_AD2_T = allocate_2D_array(dim,dim);
	double **output_AD2_temp = allocate_2D_array(dim,dim);
	double **output_ADSA = allocate_2D_array(dimDSA,dimDSA);
	double **output_ADSA_temp = allocate_2D_array(dimDSA,dimDSA);
	double **output_ADSA2 = allocate_2D_array(dimDSA,dimDSA);
	double **output_ADSA2_T = allocate_2D_array(dimDSA,dimDSA);
	double **output_ADSA2_temp = allocate_2D_array(dimDSA,dimDSA);
	double **COV_AD = allocate_2D_array(dim,dim);
	double **COV_AD2 = allocate_2D_array(dim,dim);
	double **COV_AD2_T = allocate_2D_array(dim,dim);
	double **COV_ADSA = allocate_2D_array(dimDSA,dimDSA);
	double **COV_ADSA2 = allocate_2D_array(dimDSA,dimDSA);
	for(int i = 0 ; i<dim; i++){
		for(int j = 0 ; j<dim; j++){
			output_AD[i][j]= 0;
			output_AD_temp[i][j]= 0;
			output_AD2[i][j]= 0;
			output_AD2_T[i][j]= 0;
			output_AD2_temp[i][j]= 0;
		}
	}
	for(int i = 0 ; i<dimDSA; i++){
		for(int j = 0 ; j<dimDSA; j++){
			output_ADSA[i][j]= 0;
			output_ADSA_temp[i][j]= 0;
			output_ADSA2[i][j]= 0;
			output_ADSA2_temp[i][j]= 0;
		}
	}



	double snr; 

	for(int i = 0 ;i < Ndetect; i++){
		fisher_autodiff(frequency, length, method, detectors[i],detectors[0], output_AD_temp, dim, &params, "GAUSSLEG",weights,true, psd[i],NULL,NULL);
		for(int k = 0 ; k<dim; k++){
			for(int j = 0 ; j<dim; j++){
				output_AD[k][j]+= output_AD_temp[k][j];
			}
		}
	}
	std::cout<<"SNR: "<<sqrt(output_AD[6][6])<<std::endl;

	gsl_LU_matrix_invert(output_AD,COV_AD,dim);


	double **identity_full = allocate_2D_array(dim,dim);
	matrix_multiply(output_AD,COV_AD,identity_full,dim,dim,dim);
	std::cout<<"IDENTITY: "<<std::endl;
	for(int i = 0 ; i<dim; i++){
		for(int j = 0 ; j<dim; j++){
			std::cout<<identity_full[i][j]<<" ";
		}
		std::cout<<std::endl;
	}
	
	deallocate_2D_array(identity_full, dim,dim);



	std::cout<<"ppE fishers: "<<std::endl;
	method = "ppE_IMRPhenomPv2_Inspiral";
	params.betappe[0]=beta_value;
	for(int i = 0 ;i < Ndetect; i++){
		fisher_autodiff(frequency, length, method, detectors[i],detectors[0], output_AD2_temp, dim, &params, "GAUSSLEG",weights,true, psd[i],NULL,NULL);
		for(int k = 0 ; k<dim; k++){
			for(int j = 0 ; j<dim; j++){
				output_AD2[k][j]+= output_AD2_temp[k][j];
			}
		}
	}
	params.betappe[0]=mod_value;
	ppE_theory_fisher_transformation("ppE_IMRPhenomPv2_Inspiral","DipRad_IMRPhenomPv2",dim, &params, output_AD2,output_AD2_T);
	std::cout<<"Fisher comparison: i j frac_diff F1 F2 F3"<<std::endl;	
	for(int i = 0 ; i<dim; i++){
		for(int j = 0 ; j<dim; j++){
			std::cout<<i<<" "<<j<<" "<<(output_AD[i][j]-output_AD2_T[i][j])*2./
				(fabs(output_AD[i][j])+fabs(output_AD2_T[i][j]))
				//<<std::endl;
				<<" "<<output_AD[i][j]<<" "<<output_AD2_T[i][j]<<" "<<output_AD2[i][j]<<std::endl;;
		}
	}
	gsl_LU_matrix_invert(output_AD2,COV_AD2,dim);
	ppE_theory_covariance_transformation("ppE_IMRPhenomPv2_Inspiral","DipRad_IMRPhenomPv2",dim, &params, COV_AD2,COV_AD2_T);
	double cov_ppe[dim-1];
	double cov_dcs[dim-1];
	std::cout<<"Cov comparison: dcs ppe dCS-cov-comp ppE-cov-comp"<<std::endl;
	for(int i = 0 ; i<dim-1; i++){
		cov_ppe[i] = COV_AD2[i][dim-1] / sqrt( COV_AD2[i][i]*COV_AD2[dim-1][dim-1]);
		cov_dcs[i] = COV_AD[i][dim-1] / sqrt( COV_AD[i][i]*COV_AD[dim-1][dim-1]);
		std::cout<<cov_dcs[i]<<" "<<cov_ppe[i]<<" "<<COV_AD2[i][dim-1]<<" "<<COV_AD[i][dim-1]<<std::endl;
	}
	
	std::cout<<"Covariance difference: i j frac_diff Cov1 Cov2 ppE"<<std::endl;
	for(int i = 0 ; i<dim; i++){
		for(int j = 0 ; j<dim; j++){
			std::cout<<i<<" "<<j<<" "<<(COV_AD[i][j]-COV_AD2_T[i][j])*2./
				(fabs(COV_AD[i][j])+fabs(COV_AD2_T[i][j]))
				//<<std::endl;
				<<" "<<COV_AD[i][j]<<" "<<COV_AD2_T[i][j]<<" "<<COV_AD2[i][j]<<std::endl;;
		}
	}



	double **identity_full2 = allocate_2D_array(dim,dim);
	matrix_multiply(output_AD2_T,COV_AD2_T,identity_full2,dim,dim,dim);
	std::cout<<"IDENTITY: "<<std::endl;
	for(int i = 0 ; i<dim; i++){
		for(int j = 0 ; j<dim; j++){
			std::cout<<identity_full2[i][j]<<" ";
		}
		std::cout<<std::endl;
	}
	
	deallocate_2D_array(identity_full2, dim,dim);

	
	std::cout<<"AD-2:"<<std::endl;

	params.sky_average = true;
	params.incl_angle = 0;
	method = "DipRad_IMRPhenomD";
	for(int i = 0 ;i < Ndetect; i++){
		fisher_autodiff(frequency, length, method, detectors[i],detectors[i], output_ADSA_temp, dimDSA, &params, "GAUSSLEG",weights,true, psd[i],NULL,NULL);
		for(int k = 0 ; k<dimDSA; k++){
			for(int j = 0 ; j<dimDSA; j++){
				output_ADSA[k][j]+= output_ADSA_temp[k][j];
			}
		}
	}
	
	method = "ppE_IMRPhenomD_Inspiral";
	params.betappe[0]=beta_value;
	for(int i = 0 ;i < Ndetect; i++){
		fisher_autodiff(frequency, length, method, detectors[i],detectors[i], output_ADSA2_temp, dimDSA, &params, "GAUSSLEG",weights,true, psd[i],NULL,NULL);
		for(int k = 0 ; k<dimDSA; k++){
			for(int j = 0 ; j<dimDSA; j++){
				output_ADSA2[k][j]+= output_ADSA2_temp[k][j];
			}
		}
	}
	params.betappe[0]=mod_value;
	ppE_theory_fisher_transformation("ppE_IMRPhenomD_Inspiral","DipRad_IMRPhenomD",dimDSA, &params, output_ADSA2,output_ADSA2_T);
	
	std::cout<<"Sky averaged fisher diff: i j fracdiff F1 F2 F3"<<std::endl;
	for(int i = 0 ; i<dimDSA; i++){
		for(int j = 0 ; j<dimDSA; j++){
			std::cout<<i<<" "<<j<<" "<<(output_ADSA[i][j]-output_ADSA2_T[i][j])*2./
				(fabs(output_ADSA[i][j])+fabs(output_ADSA2_T[i][j]))
				//<<std::endl;
				<<" "<<output_ADSA[i][j]<<" "<<output_ADSA2_T[i][j]<<" "<<output_ADSA2[i][j]<<std::endl;;
		}
	}

	deallocate_2D_array(output_AD,dim,dim);
	deallocate_2D_array(output_AD_temp,dim,dim);
	deallocate_2D_array(COV_AD,dim,dim);
	deallocate_2D_array(COV_AD2,dim,dim);
	deallocate_2D_array(COV_AD2_T,dim,dim);
	deallocate_2D_array(COV_ADSA,dimDSA,dimDSA);
	deallocate_2D_array(COV_ADSA2,dimDSA,dimDSA);
	deallocate_2D_array(output_ADSA,dimDSA,dimDSA);
	deallocate_2D_array(output_ADSA_temp,dimDSA,dimDSA);
	deallocate_2D_array(output_ADSA2,dimDSA,dimDSA);
	deallocate_2D_array(output_ADSA2_T,dimDSA,dimDSA);
	deallocate_2D_array(output_ADSA2_temp,dimDSA,dimDSA);
	deallocate_2D_array(output_AD2,dim,dim);
	deallocate_2D_array(output_AD2_T,dim,dim);
	deallocate_2D_array(output_AD2_temp,dim,dim);
	
	delete [] frequency;
	for(int i = 0 ; i<Ndetect; i++){
		delete [] psd[i];
	}
	delete [] params.betappe;
	delete [] params.bppe;
	//delete [] params.phii;
	//delete [] params.delta_phi;
	delete [] psd;
	delete [] weights;
	return 0;

}
int network_fishers(int argc, char *argv[])
{
	std::cout.precision(5);
	gen_params params;	
	params.spin1[2] = .1;
	params.spin2[2] = -.1;
	params.chip = .3;
	params.phip = 1.0;
	params.Luminosity_Distance = 600;
	params.phiRef = .0;
	params.RA = 1.;
	params.DEC = -0.0;
	params.f_ref = 20;

	params.NSflag1 = false;
	params.NSflag2 = false;
	params.horizon_coord = false;
	params.shift_time=false;
	params.shift_phase=false;
	params.dep_postmerger=true;
	params.sky_average=false;
	
	params.mass1 = 15;
	params.mass2 = 8;
	double chirpmass = calculate_chirpmass(params.mass1,params.mass2)*MSOL_SEC;
	params.psi = 1.;
	params.gmst = 2.;

	//double fmin = 3e-2;
	//double fmax = 1e-1;
	//double T = T_year/2;

	double fmin = 5;
	double fmax = 2048;
	//double fmin = f_0PN(4*T_year,chirpmass);
	//double fmax = 1;
	double T = 32;
	//double T = 4*T_year;

	params.tc = 3.*T/4.;
	//params.tc = T;
	//params.tc = 0;

	//int length = T*(fmax-fmin);
	int length = 1000;
	double *frequency = new double[length];
	int Ndetect = 2;
	double **psd = new double*[Ndetect];
	std::string SN[3] = {"AdLIGODesign_smoothed","AdLIGODesign_smoothed","AdLIGODesign_smoothed"};
	
	double *weights = new double[length];
	gauleg(log10(fmin), log10(fmax),frequency,weights,length);
	for(int i = 0 ; i<length; i++){
		frequency[i] = pow(10,frequency[i]);	
	}
	for(int i = 0 ; i<Ndetect; i++){
		psd[i]= new double[length];
		populate_noise(frequency, SN[i],psd[i], length, 48);
		//populate_noise(frequency, "LISA_CONF",psd, length, 12);
		for(int j = 0 ; j<length; j++){
			psd[i][j]*=psd[i][j];	
		}
	}

	std::string method = "PNSeries_ppE_IMRPhenomPv2_Inspiral";
	//transform_orientation_coords(&params, method, detector);

	params.Nmod = 2;
	params.betappe = new double[params.Nmod];
	params.bppe = new double[params.Nmod];
	params.betappe[0] = 3;
	params.betappe[1] = 2;
	//params.bppe[0] = -1.;
	params.bppe[0] = -7;
	params.bppe[1] = -5;
	//params.bppe[0] = 1;
	//params.Nmod_phi = 1;
	//params.delta_phi = new double[1];
	//params.phii = new int[1];
	//params.delta_phi[0] = 0;
	//params.bppe[0] = -1.;
	//params.phii[0] = 4;
	//params.bppe[0] = 1;

	//std::string detectors[4] = {"CE","Hanford","Livingston","Virgo"};
	std::string detectors[3] = {"Hanford","Livingston","Virgo"};
		
	int dim = 15;

	double **jac_spins = allocate_2D_array(dim,dim);
	for (int i = 0 ;i<dim; i++){
		for(int j =0 ;j<dim; j++){
			if(i == 9 and j ==10){
				jac_spins[i][j] = .5;
			}
			else if(i == 10 and j ==9){
				jac_spins[i][j] = .5;
			}
			else if(i == 10 and j ==10){
				jac_spins[i][j] = -.5;
			}
			else if(i == 9 and j ==9){
				jac_spins[i][j] =.5;
			}
			else if(i != j ){
				jac_spins[i][j] =0;
			}
			else {
				jac_spins[i][j] =1;
			}
		}
	}

	double **output_AD = allocate_2D_array(dim,dim);
	double **output_AD_temp = allocate_2D_array(dim,dim);
	double **output_AD3 = allocate_2D_array(dim,dim);
	double **output_AD3_temp = allocate_2D_array(dim,dim);
	double **COV_AD = allocate_2D_array(dim,dim);
	double **COV_AD3 = allocate_2D_array(dim,dim);
	for(int i = 0 ; i<dim; i++){
		for(int j = 0 ; j<dim; j++){
			output_AD[i][j]= 0;
			output_AD_temp[i][j]= 0;
			output_AD3[i][j]= 0;
			output_AD3_temp[i][j]= 0;
		}
	}

	params.equatorial_orientation = false;
	params.sky_average = false;
	params.incl_angle = M_PI-.01;
	params.incl_angle = M_PI-.01;


	for(int i = 0 ;i < Ndetect; i++){
		fisher_autodiff(frequency, length, method, detectors[i],detectors[i], output_AD_temp, dim, &params, "GAUSSLEG",weights,true, psd[i],NULL,NULL);
		for(int k = 0 ; k<dim; k++){
			for(int j = 0 ; j<dim; j++){
				output_AD[k][j]+= output_AD_temp[k][j];
			}
		}
	}
	//matrix_multiply(output_AD, jac_spins,output_AD_temp,dim,dim,dim);
	//matrix_multiply(jac_spins,output_AD_temp, output_AD,dim,dim,dim);
	
	
	std::cout<<"AD:"<<std::endl;
	for(int i = 0 ; i<dim; i++){
		std::cout<<i<<" ";
		for(int j = 0 ; j<dim; j++){
			std::cout<<output_AD[i][j]<<" ";
		}
		std::cout<<std::endl;
	}

	gsl_LU_matrix_invert(output_AD,COV_AD,dim);
	std::cout<<"COV AD:"<<std::endl;
	for(int i = 0 ; i<dim; i++){
		std::cout<<i<<" ";
		for(int j = 0 ; j<dim; j++){
			std::cout<<COV_AD[i][j]<<" ";
		}
		std::cout<<std::endl;
	}
	std::cout<<"Variances:"<<std::endl;
	for(int i = 0 ; i<dim; i++){
		std::cout<<i<<" "<<1.64*sqrt(COV_AD[i][i])<<std::endl;
	}
	std::cout<<std::endl;

	for(int i = 0 ;i < Ndetect; i++){
		fisher_autodiff(frequency, length, method, detectors[i],detectors[0], output_AD3_temp, dim, &params,"GAUSSLEG",weights,true, psd[i],NULL,NULL);
		for(int k = 0 ; k<dim; k++){
			for(int j = 0 ; j<dim; j++){
				output_AD3[k][j]+= output_AD3_temp[k][j];
			}
		}
	}
	matrix_multiply(output_AD3, jac_spins,output_AD3_temp,dim,dim,dim);
	matrix_multiply(jac_spins,output_AD3_temp, output_AD3,dim,dim,dim);
	
	
	std::cout<<"AD-Network:"<<std::endl;
	for(int i = 0 ; i<dim; i++){
		std::cout<<i<<" ";
		for(int j = 0 ; j<dim; j++){
			std::cout<<output_AD3[i][j]<<" ";
		}
		std::cout<<std::endl;
	}
	gsl_LU_matrix_invert(output_AD3,COV_AD3,dim);
	std::cout<<"COV AD - Network:"<<std::endl;
	for(int i = 0 ; i<dim; i++){
		std::cout<<i<<" ";
		for(int j = 0 ; j<dim; j++){
			std::cout<<COV_AD3[i][j]<<" ";
		}
		std::cout<<std::endl;
	}
	std::cout<<"Variances:"<<std::endl;
	for(int i = 0 ; i<dim; i++){
		std::cout<<i<<" "<<1.64*sqrt(COV_AD3[i][i])<<std::endl;
	}
	std::cout<<std::endl;


	std::cout<<"FRACTIONAL DIFF (Netowrk-AD)*2/(Network+AD):"<<std::endl;
	for(int i = 0 ; i<dim; i++){
		std::cout<<i<<" ";
		for(int j = 0 ; j<dim; j++){
			std::cout<<(output_AD3[i][j] - output_AD[i][j])*2./(fabs(output_AD3[i][j]) + fabs(output_AD[i][j]))<<" ";
		}
		std::cout<<std::endl;
	}


	deallocate_2D_array(output_AD,dim,dim);
	deallocate_2D_array(output_AD_temp,dim,dim);
	deallocate_2D_array(COV_AD,dim,dim);
	deallocate_2D_array(COV_AD3,dim,dim);
	deallocate_2D_array(output_AD3,dim,dim);
	deallocate_2D_array(output_AD3_temp,dim,dim);
	deallocate_2D_array(jac_spins,dim,dim);
	
	delete [] frequency;
	for(int i = 0 ; i<Ndetect; i++){
		delete [] psd[i];
	}
	delete [] psd;
	return 0;
}
int dCS_EdGB(int argc, char *argv[])
{
	
	std::cout.precision(15);
	gen_params params;	
	//params.mass1 = 31.5;
	//params.mass2 = 10.3;
	//params.spin1[2] = .38;
	//params.spin2[2] = (.21*(params.mass1+params.mass2)-params.mass1* params.spin1[2])/params.mass2 ;
	//params.chip = .29;
	//params.phip = 1.0;
	//params.Luminosity_Distance = 730;
	//params.incl_angle = .76;
	
	params.mass1 = 5.9;
	params.mass2 = 1.4;
	//params.spin1[2] = .08* (params.mass1+params.mass2)/params.mass1;
	params.spin1[2] = .3;
	params.spin2[2] = .03 ;
	//params.spin1[1] = 0.52*sqrt(1-.51*.51)*sin(2.79);
	//params.spin2[1] = 0.43*sqrt(1-.2*.2)*sin(2.99);
	//params.spin1[0] = 0.52*sqrt(1-.51*.51)*cos(2.79);
	//params.spin2[0] = 0.43*sqrt(1-.2*.2)*cos(2.99);
	params.chip = .2;
	params.phip = 1.0;
	params.Luminosity_Distance = 30;
	params.incl_angle = 3*M_PI/4;

	//params.mass1 = 5.0;
	//params.mass2 = 1.4;
	//params.spin1[2] = .801;
	//params.spin2[2] = -.001;//(.21*(params.mass1+params.mass2)-params.mass1* params.spin1[2])/params.mass2 ;
	//params.chip = .001;
	//params.phip = 1.0;
	//params.Luminosity_Distance = 70;
	//params.incl_angle = M_PI/3.;
	

	params.NSflag1 = false;
	params.NSflag2 =true;

	params.phiRef = .0;
	params.RA = 1.;
	params.DEC = 0.6;
	params.f_ref = 20;
	double chirpmass = calculate_chirpmass(params.mass1,params.mass2)*MSOL_SEC;
	double eta = calculate_eta(params.mass1,params.mass2)*MSOL_SEC;
	//params.spin1[2] = .38;

	params.horizon_coord = false;
	params.shift_time=false;
	params.shift_phase=false;
	params.dep_postmerger=true;
	params.sky_average=false;
	
	params.psi = 1.;
	params.gmst = 2.;
	params.sky_average = false;
	//params.Nmod_phi =1;
	//params.phii = new int[1];
	//params.phii[0]=4;
	//params.delta_phi=new double[1];
	//params.delta_phi[0]=0;
	
	
	
	//params.Nmod = 2;
	//params.betappe = new double[2];
	//params.bppe = new double[2];
	////params.betappe[0] = pow_int(8*1000./(c),4);
	//params.betappe[0] =1e-30;
	//params.betappe[1] =0;
	////params.betappe[0] =1;
	//params.bppe[0] = -7;
	//params.bppe[1] = -5;
	
	params.Nmod = 1;
	params.betappe = new double[1];
	params.bppe = new double[1];
	//params.betappe[0] = pow_int(8*1000./(c),4);
	params.betappe[0] =1e-30;
	//params.betappe[0] =1;
	params.bppe[0] = -7;
	

	double fmin = 5;
	double fmax = 2048;
	double T = 16;

	params.tc = 3.*T/4.;
	int length = 1000;
	double *frequency = new double[length];
	int Ndetect = 3;
	//int Ndetect = 2;
	double **psd = new double*[Ndetect];
	//std::string SN[3] = {"AdLIGODesign_smoothed","AdLIGODesign_smoothed","AdLIGODesign_smoothed"};
	std::string SN[3] = {"AdLIGOMidHigh","AdLIGOMidHigh","AdVIRGOPlus1"};
	//std::string SN[3] = {"CE1","AdLIGOMidHigh","AdVIRGOPlus1"};
	
	double *weights = new double[length];
	gauleg(log10(fmin), log10(fmax),frequency,weights,length);
	for(int i = 0 ; i<length; i++){
		frequency[i] = pow(10,frequency[i]);	
	}
	for(int i = 0 ; i<Ndetect; i++){
		psd[i]= new double[length];
		populate_noise(frequency, SN[i],psd[i], length, 48);
		//populate_noise(frequency, "LISA_CONF",psd, length, 12);
		for(int j = 0 ; j<length; j++){
			psd[i][j]*=psd[i][j];	
		}
	}

	//std::string method = "dCS_IMRPhenomPv2";
	//std::string method = "ppE_IMRPhenomPv2_IMR";
	//std::string method = "gIMRPhenomPv2";
	//std::string method = "EdGB_IMRPhenomPv2";
	std::string method = "EdGB_HO_IMRPhenomPv2";
	//std::string method = "EdGB_GHOv1_IMRPhenomPv2";
	

	//source_parameters<double> sp;
	//std::string temp_meth = prep_source_parameters(&sp, &params, "EdGB_IMRPhenomPv2");
	//double phase_factor = EdGB_phase_factor(&sp);
	//double fisco = pow(6.,-3./2.) * pow(eta,3./5.)/(M_PI * chirpmass);
	//std::cout<<phase_factor * pow( M_PI * chirpmass * fisco,-2./3.)<<std::endl;
	//exit(0);


	std::string detectors[3] = {"Hanford","Livingston","Virgo"};
		
	int dim = 14;
	int dimD = 12;
	int dimDSA = 8;

	//int dim = 15;
	//int dimD = 13;
	//int dimDSA = 9;

	double **jac_spins = allocate_2D_array(dim,dim);
	for (int i = 0 ;i<dim; i++){
		for(int j =0 ;j<dim; j++){
			if(i == 9 and j ==10){
				jac_spins[i][j] = .5;
			}
			else if(i == 10 and j ==9){
				jac_spins[i][j] = .5;
			}
			else if(i == 10 and j ==10){
				jac_spins[i][j] = -.5;
			}
			else if(i == 9 and j ==9){
				jac_spins[i][j] =.5;
			}
			else if(i != j ){
				jac_spins[i][j] =0;
			}
			else {
				jac_spins[i][j] =1;
			}
		}
	}

	double **output_AD = allocate_2D_array(dim,dim);
	double **output_AD_temp = allocate_2D_array(dim,dim);
	double **output_AD2 = allocate_2D_array(dimD,dimD);
	double **output_ADSA = allocate_2D_array(dimDSA,dimDSA);
	double **output_ADSA_temp = allocate_2D_array(dimDSA,dimDSA);
	double **output_AD2_temp = allocate_2D_array(dimD,dimD);
	double **COV_AD = allocate_2D_array(dim,dim);
	double **COV_AD2 = allocate_2D_array(dim,dim);
	double **COV_ADSA = allocate_2D_array(dimDSA,dimDSA);
	for(int i = 0 ; i<dim; i++){
		for(int j = 0 ; j<dim; j++){
			output_AD[i][j]= 0;
			output_AD_temp[i][j]= 0;
		}
	}
	for(int i = 0 ; i<dimD; i++){
		for(int j = 0 ; j<dimD; j++){
			output_AD2[i][j]= 0;
			output_AD2_temp[i][j]= 0;
		}
	}
	for(int i = 0 ; i<dimDSA; i++){
		for(int j = 0 ; j<dimDSA; j++){
			output_ADSA[i][j]= 0;
			output_ADSA_temp[i][j]= 0;
		}
	}



	double snr; 

	for(int i = 0 ;i < Ndetect; i++){
		fisher_autodiff(frequency, length, method, detectors[i],detectors[0], output_AD_temp, dim, &params, "GAUSSLEG",weights,true, psd[i],NULL,NULL);
		for(int k = 0 ; k<dim; k++){
			for(int j = 0 ; j<dim; j++){
				output_AD[k][j]+= output_AD_temp[k][j];
			}
		}
	}
	matrix_multiply(output_AD, jac_spins,output_AD_temp,dim,dim,dim);
	matrix_multiply(jac_spins,output_AD_temp, output_AD,dim,dim,dim);
	std::cout<<"SNR: "<<sqrt(output_AD[6][6])<<std::endl;
	std::cout<<"AD:"<<std::endl;
	//for(int i = 0 ; i<dim; i++){
	//	std::cout<<i<<" ";
	//	for(int j = 0 ; j<dim; j++){
	//		std::cout<<output_AD[i][j]<<" ";
	//	}
	//	std::cout<<std::endl;
	//}

	gsl_LU_matrix_invert(output_AD,COV_AD,dim);
	//gsl_cholesky_matrix_invert(output_AD,COV_AD,dim);
	std::cout<<"COV AD:"<<std::endl;
	//for(int i = 0 ; i<dim; i++){
	//	std::cout<<i<<" ";
	//	for(int j = 0 ; j<dim; j++){
	//		std::cout<<COV_AD[i][j]<<" ";
	//	}
	//	std::cout<<std::endl;
	//}
	std::cout<<"Variances (90%):"<<std::endl;
	for(int i = 0 ; i<dim; i++){
		std::cout<<i<<" "<<1.64*sqrt(COV_AD[i][i])<<std::endl;
	}
	std::cout<<std::endl;
	std::cout<<"Variances root delta alpha (90%):"<<std::endl;
	if(params.NSflag2){
		std::cout<<"Cutoff:"<<params.mass1*.5*1.5<<std::endl;
	}
	else{
		std::cout<<"Cutoff:"<<params.mass2*.5*1.5<<std::endl;
	}
	if(method.find("GHO") == std::string::npos){
		std::cout<<"(delta alpha^2)^(1/4) (KM) (90%): "<<1.64*pow(COV_AD[dim-1][dim-1],1./8.)*3.e5<<std::endl;
		std::cout<<"(delta alpha^2)^(1/4) (KM): "<<pow(COV_AD[dim-1][dim-1],1./8.)*3.e5<<std::endl;
	}
	else{
		std::cout<<"(delta alpha^2)^(1/4) (KM) (90%): "<<1.64*pow(COV_AD[dim-2][dim-2],1./8.)*3.e5<<std::endl;
		std::cout<<"(delta alpha^2)^(1/4) (KM): "<<pow(COV_AD[dim-2][dim-2],1./8.)*3.e5<<std::endl;
	}
	std::cout<<std::endl;


	double **identity_full = allocate_2D_array(dim,dim);
	matrix_multiply(output_AD,COV_AD,identity_full,dim,dim,dim);
	//std::cout<<"IDENTITY: "<<std::endl;
	//for(int i = 0 ; i<dim; i++){
	//	for(int j = 0 ; j<dim; j++){
	//		std::cout<<identity_full[i][j]<<" ";
	//	}
	//	std::cout<<std::endl;
	//}
	
	deallocate_2D_array(identity_full, dim,dim);




	double **sub_AD_F = allocate_2D_array(dim-4,dim-4);
	int ids[4] = {0,1,2,3};
	rm_fisher_dim(output_AD,dim, sub_AD_F,dim-4,ids);
	gsl_LU_matrix_invert(sub_AD_F,COV_AD,dim-4);
	std::cout<<"SUB COV AD:"<<std::endl;
	//for(int i = 0 ; i<dim-4; i++){
	//	std::cout<<i<<" ";
	//	for(int j = 0 ; j<dim-4; j++){
	//		std::cout<<sub_AD_F[i][j]<<" ";
	//	}
	//	std::cout<<std::endl;
	//}
	std::cout<<"Variances (90%):"<<std::endl;
	for(int i = 0 ; i<dim-4; i++){
		std::cout<<i<<" "<<1.64*sqrt(COV_AD[i][i])<<std::endl;
	}
	std::cout<<std::endl;
	std::cout<<"Variances root delta alpha (90%):"<<std::endl;
	if(params.NSflag2){
		std::cout<<"Cutoff:"<<params.mass1*.5*1.5<<std::endl;
	}
	else{
		std::cout<<"Cutoff:"<<params.mass2*.5*1.5<<std::endl;
	}
	if(method.find("GHO") == std::string::npos){
		std::cout<<"(delta alpha^2)^(1/4) (KM) (90%): "<<1.64*pow(COV_AD[dim-5][dim-5],1./8.)*3.e5<<std::endl;
		std::cout<<"(delta alpha^2)^(1/4) (KM): "<<pow(COV_AD[dim-5][dim-5],1./8.)*3.e5<<std::endl;
	}
	else{
		std::cout<<"(delta alpha^2)^(1/4) (KM) (90%): "<<1.64*pow(COV_AD[dim-5][dim-5],1./8.)*3.e5<<std::endl;
		std::cout<<"(delta alpha^2)^(1/4) (KM): "<<pow(COV_AD[dim-5][dim-5],1./8.)*3.e5<<std::endl;
	}
	std::cout<<std::endl;
	deallocate_2D_array(sub_AD_F,dim-4,dim-4);






	//method = "dCS_IMRPhenomD";
	//method = "ppE_IMRPhenomD_IMR";
	//method = "gIMRPhenomD";
	//method = "EdGB_IMRPhenomD";
	method = "EdGB_HO_IMRPhenomD";
	//method = "EdGB_GHOv1_IMRPhenomD";
	for(int i = 0 ;i < Ndetect; i++){
		fisher_autodiff(frequency, length, method, detectors[i],detectors[0], output_AD2_temp, dimD, &params, "GAUSSLEG",weights,true, psd[i],NULL,NULL);
		for(int k = 0 ; k<dimD; k++){
			for(int j = 0 ; j<dimD; j++){
				output_AD2[k][j]+= output_AD2_temp[k][j];
			}
		}
	}
	//matrix_multiply(output_AD2, jac_spins,output_AD2_temp,dimD,dimD,dimD);
	//matrix_multiply(jac_spins,output_AD2_temp, output_AD2,dimD,dimD,dimD);
	
	
	std::cout<<"AD-D:"<<std::endl;
	//for(int i = 0 ; i<dimD; i++){
	//	std::cout<<i<<" ";
	//	for(int j = 0 ; j<dimD; j++){
	//		std::cout<<output_AD2[i][j]<<" ";
	//	}
	//	std::cout<<std::endl;
	//}
	gsl_LU_matrix_invert(output_AD2,COV_AD2,dimD);
	std::cout<<"COV AD - D:"<<std::endl;
	std::cout<<"Variances (90%):"<<std::endl;
	for(int i = 0 ; i<dimD; i++){
		std::cout<<i<<" "<<1.64*sqrt(COV_AD2[i][i])<<std::endl;
	}
	std::cout<<std::endl;
	std::cout<<"Variances root delta alpha (90%):"<<std::endl;
	if(params.NSflag2){
		std::cout<<"Cutoff:"<<params.mass1*.5*1.5<<std::endl;
	}
	else{
		std::cout<<"Cutoff:"<<params.mass2*.5*1.5<<std::endl;
	}
	if(method.find("GHO") == std::string::npos){
		std::cout<<"(delta alpha^2)^(1/4) (KM): "<<1.64*pow(COV_AD2[dimD-1][dimD-1],1./8.)*3.e5<<std::endl;
	}
	else{
		std::cout<<"(delta alpha^2)^(1/4) (KM): "<<1.64*pow(COV_AD2[dimD-2][dimD-2],1./8.)*3.e5<<std::endl;
	}
	std::cout<<std::endl;

	params.sky_average = true;
	params.incl_angle = 0;
	for(int i = 0 ;i < Ndetect; i++){
		fisher_autodiff(frequency, length, method, detectors[i],detectors[i], output_ADSA_temp, dimDSA, &params, "GAUSSLEG",weights,true, psd[i],NULL,NULL);
		for(int k = 0 ; k<dimDSA; k++){
			for(int j = 0 ; j<dimDSA; j++){
				output_ADSA[k][j]+= output_ADSA_temp[k][j];
			}
		}
	}
	
	
	std::cout<<"AD-DSA:"<<std::endl;
	//for(int i = 0 ; i<dimDSA; i++){
	//	std::cout<<i<<" ";
	//	for(int j = 0 ; j<dimDSA; j++){
	//		std::cout<<output_ADSA[i][j]<<" ";
	//	}
	//	std::cout<<std::endl;
	//}
	gsl_LU_matrix_invert(output_ADSA,COV_ADSA,dimDSA);
	double **identity = allocate_2D_array(dimDSA,dimDSA);
	matrix_multiply(output_ADSA,COV_ADSA,identity,dimDSA,dimDSA,dimDSA);
	//std::cout<<"IDENTITY: "<<std::endl;
	//for(int i = 0 ; i<dimDSA; i++){
	//	for(int j = 0 ; j<dimDSA; j++){
	//		std::cout<<identity[i][j]<<" ";
	//	}
	//	std::cout<<std::endl;
	//}
	
	deallocate_2D_array(identity, dimDSA,dimDSA);
	std::cout<<"Variances (90%):"<<std::endl;
	for(int i = 0 ; i<dimDSA; i++){
		std::cout<<i<<" "<<1.64*sqrt(COV_ADSA[i][i])<<std::endl;
	}
	std::cout<<std::endl;
	std::cout<<"Variances root delta alpha (90%):"<<std::endl;
	if(params.NSflag2){
		std::cout<<"Cutoff:"<<params.mass1*.5*1.5<<std::endl;
	}
	else{
		std::cout<<"Cutoff:"<<params.mass2*.5*1.5<<std::endl;
	}
	if(method.find("GHO") == std::string::npos){
		std::cout<<"(delta alpha^2)^(1/4) (KM): "<<1.64*pow(COV_ADSA[dimDSA-1][dimDSA-1],1./8.)*3.e5<<std::endl;
	}
	else{
		std::cout<<"(delta alpha^2)^(1/4) (KM): "<<1.64*pow(COV_ADSA[dimDSA-2][dimDSA-2],1./8.)*3.e5<<std::endl;

	}
	std::cout<<std::endl;
	std::cout<<"SNR (SA): "<<sqrt(output_ADSA[0][0])<<std::endl;

	deallocate_2D_array(output_AD,dim,dim);
	deallocate_2D_array(output_AD_temp,dim,dim);
	deallocate_2D_array(COV_AD,dim,dim);
	deallocate_2D_array(COV_AD2,dim,dim);
	deallocate_2D_array(COV_ADSA,dimDSA,dimDSA);
	deallocate_2D_array(output_ADSA,dimDSA,dimDSA);
	deallocate_2D_array(output_ADSA_temp,dimDSA,dimDSA);
	deallocate_2D_array(output_AD2,dimD,dimD);
	deallocate_2D_array(output_AD2_temp,dimD,dimD);
	deallocate_2D_array(jac_spins,dim,dim);
	
	delete [] frequency;
	for(int i = 0 ; i<Ndetect; i++){
		delete [] psd[i];
	}
	delete [] params.betappe;
	delete [] params.bppe;
	//delete [] params.phii;
	//delete [] params.delta_phi;
	delete [] psd;
	delete [] weights;
	return 0;
}
int AD_v_N(int argc, char *argv[])
{
	
	std::cout.precision(5);
	gen_params params;	
	params.spin1[2] = .1;
	params.spin2[2] = -.1;
	params.chip = .03;
	params.phip = 1.0;
	params.Luminosity_Distance = 1000;
	params.phiRef = .0;
	params.RA = 2.;
	params.DEC = -0.9;
	params.f_ref = 20;

	params.NSflag1 = false;
	params.NSflag2 = false;
	params.horizon_coord = false;
	params.shift_time=false;
	params.shift_phase=false;
	params.dep_postmerger=true;
	params.sky_average=false;
	
	params.mass1 = 36;
	params.mass2 = 29;
	double chirpmass = calculate_chirpmass(params.mass1,params.mass2)*MSOL_SEC;
	//params.theta_l = 1.1;
	//params.phi_l = 2.1;
	params.psi = 1.;
	params.gmst = 2.;

	//double fmin = 3e-2;
	//double fmax = 1e-1;
	//double T = T_year/2;

	double fmin = 5;
	double fmax = 2048;
	//double fmin = f_0PN(4*T_year,chirpmass);
	//double fmax = 1;
	double T = 32;
	//double T = 4*T_year;

	params.tc = 3.*T/4.;
	//params.tc = T;
	//params.tc = 0;

	//int length = T*(fmax-fmin);
	int length = 1000;
	double *frequency = new double[length];
	//int Ndetect = 4;
	int Ndetect = 1;
	double **psd = new double*[Ndetect];
	//std::string SN[4] = {"CE2_smoothed","AdLIGODesign_smoothed","AdLIGODesign_smoothed","AdLIGODesign_smoothed"};
	std::string SN[3] = {"AdLIGODesign_smoothed","AdLIGODesign_smoothed","AdLIGODesign_smoothed"};
	//std::string SN[4] = {"LISA_SADC_CONF","AdLIGODesign_smoothed","AdLIGODesign_smoothed","AdLIGODesign_smoothed"};
	//std::string SN[3] = {"Hanford_O1_fitted","Hanford_O1_fitted","Hanford_O1_fitted"};
	
	double *weights = new double[length];
	gauleg(log10(fmin), log10(fmax),frequency,weights,length);
	for(int i = 0 ; i<length; i++){
		frequency[i] = pow(10,frequency[i]);	
	}
	//for(int i = 0 ; i<length; i++){
	//	frequency[i]=fmin + (double)i /T;
	//}
	for(int i = 0 ; i<Ndetect; i++){
		psd[i]= new double[length];
		populate_noise(frequency, SN[i],psd[i], length, 48);
		//populate_noise(frequency, "LISA_CONF",psd, length, 12);
		for(int j = 0 ; j<length; j++){
			psd[i][j]*=psd[i][j];	
		}
	}

	//std::string detector = "LISA";
	std::string detector = "Hanford";
	//std::string method = "IMRPhenomD";
	//std::string method = "ppE_IMRPhenomPv2_Inspiral";
	//std::string method = "ppE_IMRPhenomD_Inspiral";
	//std::string method = "gIMRPhenomD";
	std::string method = "gIMRPhenomPv2";
	//transform_orientation_coords(&params, method, detector);

	//params.Nmod = 1;
	//params.betappe = new double[1];
	//params.bppe = new int[1];
	//params.betappe[0] = 0;
	//params.bppe[0] = -1.;
	//params.bppe[0] = -1;
	//params.bppe[0] = 1;
	//
	params.Nmod_phi = 1;
	params.delta_phi = new double[1];
	params.phii = new int[1];
	params.delta_phi[0] = 0;
	params.phii[0] = 2;
	std::cout<<"Phase power: "<<params.phii[0]<<std::endl;

	//params.Nmod_alpha = 1;
	//params.delta_alpha = new double[1];
	//params.alphai = new int[1];
	//params.delta_alpha[0] = 0;
	//params.alphai[0] = 2;

	//std::string detectors[4] = {"CE","Hanford","Livingston","Virgo"};
	std::string detectors[3] = {"Hanford","Livingston","Virgo"};
	//std::string detectors[4] = {"LISA","Hanford","Livingston","Virgo"};
		
	int dim = 13;
	int dimD = 12;
	int dimDSA = 8;

	double **jac_spins = allocate_2D_array(dim,dim);
	for (int i = 0 ;i<dim; i++){
		for(int j =0 ;j<dim; j++){
			if(i == 9 and j ==10){
				jac_spins[i][j] = .5;
			}
			else if(i == 10 and j ==9){
				jac_spins[i][j] = .5;
			}
			else if(i == 10 and j ==10){
				jac_spins[i][j] = -.5;
			}
			else if(i == 9 and j ==9){
				jac_spins[i][j] =.5;
			}
			else if(i != j ){
				jac_spins[i][j] =0;
			}
			else {
				jac_spins[i][j] =1;
			}
		}
	}

	double **output_N = allocate_2D_array(dim,dim);
	double **output_N_temp = allocate_2D_array(dim,dim);
	double **output_AD = allocate_2D_array(dim,dim);
	double **output_AD_temp = allocate_2D_array(dim,dim);
	double **output_AD3 = allocate_2D_array(dimD,dimD);
	double **output_ADSA = allocate_2D_array(dimDSA,dimDSA);
	double **output_ADSA_temp = allocate_2D_array(dimDSA,dimDSA);
	double **output_AD3_temp = allocate_2D_array(dimD,dimD);
	double **COV_AD = allocate_2D_array(dim,dim);
	double **COV_AD3 = allocate_2D_array(dim,dim);
	double **COV_ADSA = allocate_2D_array(dimDSA,dimDSA);
	double **output_AD2 = allocate_2D_array(dim,dim);
	for(int i = 0 ; i<dim; i++){
		for(int j = 0 ; j<dim; j++){
			output_AD[i][j]= 0;
			output_AD_temp[i][j]= 0;
			output_N[i][j]= 0;
		}
	}
	for(int i = 0 ; i<dimD; i++){
		for(int j = 0 ; j<dimD; j++){
			output_AD3[i][j]= 0;
			output_AD3_temp[i][j]= 0;
		}
	}
	for(int i = 0 ; i<dimDSA; i++){
		for(int j = 0 ; j<dimDSA; j++){
			output_ADSA[i][j]= 0;
			output_ADSA_temp[i][j]= 0;
		}
	}

	//params.equatorial_orientation = true;
	params.equatorial_orientation = false;
	params.sky_average = true;
	//params.theta_l = 0.12;
	//params.phi_l = 0.1;
	//params.incl_angle  = 0;
	double snr;
	//snr = calculate_snr(SN[0],"LISA",method, &params, frequency, length, "GAUSSLEG",weights,true);
	//std::cout<<snr<<std::endl;
	//double SNR_TARGET = 100;
	//params.Luminosity_Distance = snr/SNR_TARGET*params.Luminosity_Distance;
	params.sky_average = false;
	params.incl_angle = .33;



	//for(int i = 0 ;i < Ndetect; i++){
	//	fisher_numerical(frequency, length, method, detectors[i], output_N_temp, dim, &params, 2, NULL,NULL, psd[i]);
	//	for(int k = 0 ; k<dim; k++){
	//		for(int j = 0 ; j<dim; j++){
	//			output_N[k][j]+= output_N_temp[k][j];
	//		}
	//	}
	//}
	
	//std::cout<<"NUMERICAL:"<<std::endl;
	//for(int i = 0 ; i<dim; i++){
	//	std::cout<<i<<" ";
	//	for(int j = 0 ; j<dim; j++){
	//		std::cout<<output_N[i][j]<<" ";
	//	}
	//	std::cout<<std::endl;
	//}

	for(int i = 0 ;i < Ndetect; i++){
		fisher_autodiff(frequency, length, method, detectors[i],detectors[0], output_AD_temp, dim, &params, "GAUSSLEG",weights,true, psd[i],NULL,NULL);
		for(int k = 0 ; k<dim; k++){
			for(int j = 0 ; j<dim; j++){
				output_AD[k][j]+= output_AD_temp[k][j];
			}
		}
	}
	//matrix_multiply(output_AD, jac_spins,output_AD_temp,dim,dim,dim);
	//matrix_multiply(jac_spins,output_AD_temp, output_AD,dim,dim,dim);
	std::cout<<"SNR: "<<sqrt(output_AD[6][6])<<std::endl;
	snr = calculate_snr(SN[0],"CE",method, &params, frequency, length, "GAUSSLEG",weights,true);
	std::cout<<"SNR: "<<snr<<std::endl;
	//snr = calculate_snr(SN[1],"Hanford",method, &params, frequency, length, "SIMPSONS",NULL,false);
	//std::cout<<"SNR (Hanford): "<<snr<<std::endl;
	
	
	//std::cout<<"AD:"<<std::endl;
	//for(int i = 0 ; i<dim; i++){
	//	std::cout<<i<<" ";
	//	for(int j = 0 ; j<dim; j++){
	//		std::cout<<output_AD[i][j]<<" ";
	//	}
	//	std::cout<<std::endl;
	//}

	gsl_LU_matrix_invert(output_AD,COV_AD,dim);
	std::cout<<"COV AD:"<<std::endl;
	//for(int i = 0 ; i<dim; i++){
	//	std::cout<<i<<" ";
	//	for(int j = 0 ; j<dim; j++){
	//		std::cout<<COV_AD[i][j]<<" ";
	//	}
	//	std::cout<<std::endl;
	//}
	std::cout<<"Variances:"<<std::endl;
	for(int i = 0 ; i<dim; i++){
		std::cout<<i<<" "<<1.64*sqrt(COV_AD[i][i])<<std::endl;
	}
	std::cout<<std::endl;

	//method = "ppE_IMRPhenomD_Inspiral";
	method = "gIMRPhenomD";
	for(int i = 0 ;i < Ndetect; i++){
		fisher_autodiff(frequency, length, method, detectors[i],detectors[0], output_AD3_temp, dimD, &params, "GAUSSLEG",weights,true, psd[i],NULL,NULL);
		for(int k = 0 ; k<dimD; k++){
			for(int j = 0 ; j<dimD; j++){
				output_AD3[k][j]+= output_AD3_temp[k][j];
			}
		}
	}
	//matrix_multiply(output_AD3, jac_spins,output_AD3_temp,dimD,dimD,dimD);
	//matrix_multiply(jac_spins,output_AD3_temp, output_AD3,dimD,dimD,dimD);
	
	
	//std::cout<<"AD-D:"<<std::endl;
	//for(int i = 0 ; i<dimD; i++){
	//	std::cout<<i<<" ";
	//	for(int j = 0 ; j<dimD; j++){
	//		std::cout<<output_AD3[i][j]<<" ";
	//	}
	//	std::cout<<std::endl;
	//}
	gsl_LU_matrix_invert(output_AD3,COV_AD3,dimD);
	std::cout<<"COV AD - D:"<<std::endl;
	//for(int i = 0 ; i<dimD; i++){
	//	std::cout<<i<<" ";
	//	for(int j = 0 ; j<dimD; j++){
	//		std::cout<<COV_AD3[i][j]<<" ";
	//	}
	//	std::cout<<std::endl;
	//}
	std::cout<<"Variances:"<<std::endl;
	for(int i = 0 ; i<dimD; i++){
		std::cout<<i<<" "<<1.64*sqrt(COV_AD3[i][i])<<std::endl;
	}
	std::cout<<std::endl;

	//method = "ppE_IMRPhenomD_Inspiral";
	//method = "IMRPhenomD";
	//dimDSA = 7;
	//params.sky_average = true;
	//params.incl_angle = 0;
	//for(int i = 0 ;i < Ndetect; i++){
	//	fisher_autodiff(frequency, length, method, detectors[i],detectors[i], output_ADSA_temp, dimDSA, &params, "GAUSSLEG",weights,true, psd[i],NULL,NULL);
	//	for(int k = 0 ; k<dimDSA; k++){
	//		for(int j = 0 ; j<dimDSA; j++){
	//			output_ADSA[k][j]+= output_ADSA_temp[k][j];
	//		}
	//	}
	//}
	//
	//
	//std::cout<<"AD-DSA:"<<std::endl;
	////for(int i = 0 ; i<dimDSA; i++){
	////	std::cout<<i<<" ";
	////	for(int j = 0 ; j<dimDSA; j++){
	////		std::cout<<output_ADSA[i][j]<<" ";
	////	}
	////	std::cout<<std::endl;
	////}
	//gsl_LU_matrix_invert(output_ADSA,COV_ADSA,dimDSA);
	//std::cout<<"COV AD - DSA:"<<std::endl;
	////for(int i = 0 ; i<dimDSA; i++){
	////	std::cout<<i<<" ";
	////	for(int j = 0 ; j<dimDSA; j++){
	////		std::cout<<COV_ADSA[i][j]<<" ";
	////	}
	////	std::cout<<std::endl;
	////}
	//std::cout<<"Variances:"<<std::endl;
	//for(int i = 0 ; i<dimDSA; i++){
	//	std::cout<<i<<" "<<1.64*sqrt(COV_ADSA[i][i])<<std::endl;
	//}
	//std::cout<<std::endl;
	//std::cout<<"SNR (SA): "<<sqrt(output_ADSA[0][0])<<std::endl;
	//snr = calculate_snr(SN[0],"CE",method, &params, frequency, length, "GAUSSLEG",weights,true);
	//std::cout<<"SNR (SA -- CE): "<<snr<<std::endl;
	//snr = calculate_snr(SN[1],"Hanford",method, &params, frequency, length, "SIMPSONS",NULL,false);
	//std::cout<<"SNR (SA -- Hanford): "<<snr<<std::endl;



	//##########################################################################
	//TESTING
	//If you use 13 parameters (phip and phiRef), these two parameters are exactly,
	//linearly correlated. Simply use one. 
	//##########################################################################
	params.phiRef-=.5;
	params.phip+=.5;
	//fisher_autodiff(frequency, length, method, detectors[i],detectors[0], output_AD_temp, dim, &params, "GAUSSLEG",weights,true, psd[i],NULL,NULL);
	method = "gIMRPhenomPv2";
	fisher_autodiff(frequency, length, method, detectors[0],detectors[0], output_AD2, dim, &params, "GAUSSLEG",weights,true, psd[0],NULL,NULL);
	
	std::cout<<"FD AD / AD2:"<<std::endl;
	for(int i = 0 ; i<dim; i++){
		//std::cout<<i<<" ";
		for(int j = 0 ; j<dim; j++){
			if((output_AD2[i][j] - output_AD[i][j])*2./(output_AD2[i][j] + output_AD[i][j]) > 1e-10){ 
				std::cout<<i<<" "<<j<<" "<<output_AD2[i][j]<<" "<<output_AD[i][j]<<std::endl;
			}
			//std::cout<<(output_AD2[i][j] - output_AD[i][j])*2./(output_AD2[i][j] + output_AD[i][j])<<" ";
		}
		std::cout<<std::endl;
	}
	//##########################################################################
	//##########################################################################
	//std::cout<<"FRACTIONAL DIFF (N-AD)*2/(N+AD):"<<std::endl;
	//for(int i = 0 ; i<dim; i++){
	//	std::cout<<i<<" ";
	//	for(int j = 0 ; j<dim; j++){
	//		std::cout<<(output_N[i][j] - output_AD[i][j])*2./(fabs(output_N[i][j]) + fabs(output_AD[i][j]))<<" ";
	//	}
	//	std::cout<<std::endl;
	//}

	//std::cout<<"DIFF (N-AD):"<<std::endl;
	//for(int i = 0 ; i<dim; i++){
	//	std::cout<<i<<" ";
	//	for(int j = 0 ; j<dim; j++){
	//		std::cout<<(output_N[i][j] - output_AD[i][j])<<" ";
	//	}
	//	std::cout<<std::endl;
	//}

	deallocate_2D_array(output_AD,dim,dim);
	deallocate_2D_array(output_AD_temp,dim,dim);
	deallocate_2D_array(COV_AD,dim,dim);
	deallocate_2D_array(COV_AD3,dim,dim);
	deallocate_2D_array(COV_ADSA,dimDSA,dimDSA);
	deallocate_2D_array(output_AD2,dim,dim);
	deallocate_2D_array(output_ADSA,dimDSA,dimDSA);
	deallocate_2D_array(output_ADSA_temp,dimDSA,dimDSA);
	deallocate_2D_array(output_AD3,dimD,dimD);
	deallocate_2D_array(output_AD3_temp,dimD,dimD);
	deallocate_2D_array(output_N,dim,dim);
	deallocate_2D_array(output_N_temp,dim,dim);
	deallocate_2D_array(jac_spins,dim,dim);
	
	delete [] frequency;
	for(int i = 0 ; i<Ndetect; i++){
		delete [] psd[i];
	}
	//delete [] params.betappe;
	//delete [] params.bppe;
	delete [] params.phii;
	delete [] params.delta_phi;
	delete [] psd;
	delete [] weights;
	return 0;
}
void RT_ERROR_MSG()
{
	std::cout<<"ERROR -- incorrect arguments"<<std::endl;
	std::cout<<"Please supply function option:"<<std::endl;
	std::cout<<"0 --- Compare AD to Numerical"<<std::endl;
	std::cout<<"1 --- Network Fishers"<<std::endl;
	std::cout<<"2 --- dCS or EdGB"<<std::endl;
	std::cout<<"3 --- Check ppE-theory jac transform"<<std::endl;
	std::cout<<"4 --- Test MCMC fisher"<<std::endl;
	std::cout<<"5 --- Test LISA fisher"<<std::endl;
	std::cout<<"6 --- Test EA fisher"<<std::endl;
}


