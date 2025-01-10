#include <math.h>
#include "gwat/util.h"
#include "gwat/io_util.h"
#include "gwat/fisher.h"
#include "gwat/detector_util.h"
#include "gwat/waveform_util.h"
#include "gwat/ortho_basis.h"
#include "gwat/pn_waveform_util.h"
#include "gwat/ppE_utilities.h"
#include "gwat/error.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

int wf_consistency_test(int argc, char *argv[]);
int test_error(int argc, char *argv[]); 
void RT_ERROR_MSG();

int main(int argc, char *argv[]){

	if(argc != 2){
		RT_ERROR_MSG();
		return 1;
	}
	
	int runtime_opt = stoi(argv[1]);	
	if(runtime_opt == 0){
			return wf_consistency_test(argc,argv);
	}
	if(runtime_opt == 1){
			return test_error(argc,argv);
	}
	else{
		RT_ERROR_MSG();
		return 1;
	}
}

int wf_consistency_test(int argc, char *argv[]){
	std::cout<<"WAVEFORM CONSISTENCY TEST"<<std::endl;
	gen_params params;	
	params.spin1[1] = .0;
	params.spin2[1] = .0;
	params.spin1[0] = .0;
	params.spin2[0] = .0;
	params.f_ref = 20;
	params.NSflag1 = true;
	params.NSflag2 = true;
	params.horizon_coord = false;
	params.shift_time=true;
	params.shift_phase=true;
	
	params.tc = 6;
	params.equatorial_orientation = false;
	params.psi = 1.;
	params.gmst=3;
	params.tidal_love = true;
	params.mass1 = 1.44;
	params.mass2 = 1.29399;
	params.spin1[2] = .003;
	params.spin2[2] = -.002 ;
	params.Luminosity_Distance = 40;
	params.incl_angle = 2.532207345558998;
	params.phiRef = 2.;
	params.RA = 3.42;
	params.DEC = -.37;
	params.tidal_s=242;
	source_parameters<double> sp ;


	//int iterations = 1;
	int samples = 2*8032;
	double **output = allocate_2D_array(samples, 7);
		

	double FMIN = 5;
	double FMAX = 2608.32; //Set to match 1.2*fmerger for this system.
	/* IF YOU CHANGE THE PARAMETERS, NEED TO UPDATE THIS NUMBER. This is
	 * calculated in the taper function of IMRPhenomD_NRT and it's easier to
	 * print it from there than copying all that code over here.
	 */
	//double FMAX = 4096;
	//double FMAX = 100;
	double deltaf = (FMAX-FMIN)/samples;


	double *freqs= new double[samples];
	for (int i = 0 ; i<samples; i++){
		freqs[i] = FMIN + deltaf*i;
	}

	const gsl_rng_type *T;
	gsl_rng *r ;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	double *testPSD=new double[samples];
	populate_noise(freqs, "AdLIGOMidHigh", testPSD, samples, 48);
	for(int i = 0 ; i<samples; i++){
		testPSD[i] *= testPSD[i];
	}

	//for (int i = 0 ; i<iterations; i++){
	  
	  //params.mass1 = gsl_rng_uniform(r) +1;
	  //params.mass2 = gsl_rng_uniform(r) +1;
	  //if(params.mass2>params.mass1){
	    //double temp = params.mass2;
	    //params.mass2 = params.mass1;
	    //params.mass1 = temp;
	  //}

	  //params.spin1[2] = gsl_rng_uniform(r)*.05 -.025;
	  //params.spin2[2] = gsl_rng_uniform(r)*.05 -.025;

	  //params.tidal_s = gsl_rng_uniform(r)*1000+5; 
	  //params.tidal1 = gsl_rng_uniform(r)*100+5;
	  //params.tidal2 = gsl_rng_uniform(r)*100+5;

	  
	  double beta = .01;
	  //int b = 11.;
	  int b = -5.;

	  std::string setdir = "data/error/";
	  if(beta == 0){
	    setdir += "betazero/"; 
	  }
	  else{
	    setdir += "b"+std::to_string(abs(b))+"/"; 
	  }
	
	  params.Nmod = 1;
	  params.bppe = new double[1];
	  params.bppe[0] = b;
	  params.betappe = new double[1];
	  params.betappe[0] = beta;
		
	  std::complex<double> *response1 =  new std::complex<double>[samples];
	  std::complex<double> *response2 =  new std::complex<double>[samples];

	  //std::string method1 = "IMRPhenomD";
	  //std::string method2 = "ppE_IMRPhenomD_IMR";
	  std::string method1 = "IMRPhenomD_NRT";
	  //std::string method2 = "ppE_IMRPhenomD_NRT_Inspiral";
	  std::string method2 = "ppE_IMRPhenomD_NRT_IMR";
	  fourier_detector_response(freqs, samples, response1, "Hanford", method1, &params, (double *) NULL);
	  fourier_detector_response(freqs, samples, response2, "Hanford", method2, &params, (double *) NULL);

	  double matchresponse = match(response1, response2, testPSD, freqs, samples);

	  std::cout<<"Match: "<<matchresponse<<std::endl;
	  
	  double *phase_1 = new double[samples];
	  double *phase_2 = new double[samples];
	  double *phase_1_unwrap = new double[samples];
	  double *phase_2_unwrap = new double[samples];
	  for(int i = 0 ; i<samples ; i++){
	    phase_1[i]= std::atan2(std::imag(response1[i]),std::real(response1[i]));
	    phase_2[i]= std::atan2(std::imag(response2[i]),std::real(response2[i]));
	  }
	  unwrap_array(phase_1, phase_1_unwrap,samples);
	  unwrap_array(phase_2, phase_2_unwrap,samples);
	
	  for(int j = 0 ; j < samples ; j++){
	    output[j][0] = std::real(response1[j]);
	    output[j][1] = std::imag(response1[j]);
	    output[j][2] = std::real(response2[j]);
	    output[j][3] = std::imag(response2[j]);
	    output[j][4] = phase_1_unwrap[j];
	    output[j][5] = phase_2_unwrap[j];
	    output[j][6] = freqs[j]; 
	  }
	  //write_file("data/error/waveform_COMP_"+std::to_string(i)+".csv", output, samples, 6);
	  write_file(setdir+method1+"_comp_"+method2+".csv", output, samples, 7);
	  //write_file("data/error/b"+std::to_string(b)+".csv", output, samples, 6);

	  
	  delete [] response1;
	  delete [] response2;
	  delete [] phase_1;
	  delete [] phase_2;
	  delete [] phase_1_unwrap;
	  delete [] phase_2_unwrap;
	  //	}
	delete [] testPSD;
	gsl_rng_free(r);
	delete [] freqs;
	deallocate_2D_array(output, samples, 4);
	delete [] params.betappe;
	delete [] params.bppe;
	return 0;
}

int test_error(int argc, char *argv[])
{
	/*
	The params come from test_LISA_fishers in test_fishers.cpp
	Except I have substituted the detector to LIGO
	*/


    //std::cout.precision(15);
	//Create injection structure
	gen_params params;	
	//params.mass1 = 1.9 *MSOL_SEC;
	params.mass1 = 1.44;
	params.mass2 = 1.29399;
	params.spin1[2] = .003;
	params.spin2[2] = -.002 ;
	params.Luminosity_Distance = 40;
	params.incl_angle = 2.532207345558998;

	params.NSflag1 = true;
	params.NSflag2 =true;
	//params.NSflag1 = false;
	//params.NSflag2 =false;

	params.phiRef = 2.;
	params.RA = 3.42;
	params.DEC = -.37;
	params.f_ref = 20;
	double chirpmass = calculate_chirpmass(params.mass1,params.mass2);//*MSOL_SEC;
	double eta = calculate_eta(params.mass1,params.mass2);//*MSOL_SEC;

	params.horizon_coord = false;
	params.shift_time=false;
	params.shift_phase=false;
	params.dep_postmerger=true;
	params.sky_average=false;
	params.tidal_love=true;
	params.tidal_s=242;
	//params.tidal_s=2420;
	//params.tidal1 = 0;
	//params.tidal2 = 0; 
	

	double beta = 0.;
	//int b = 5.;
	int b = -7.;
	std::cout<<"Beta: "<<beta<<std::endl;
	std::cout<<"b: "<<b<<std::endl;

	
	
	params.Nmod = 1;
	params.bppe = new double[1];
	params.bppe[0] = b;
	params.betappe = new double[1];
	params.betappe[0] = beta;

	params.alphappe = new double[1];
	params.appe = new double[1];
	params.alphappe[0] = 0.;
	params.appe[0] = -2.;

	
	
	//params.tidal_love=false;
	//params.tidal1=200;
	//params.tidal2=100;
	params.psi = 2.;
	double gps = 1187008882.4;
	params.gmst = gps_to_GMST_radian(gps);
	params.sky_average = false;

	

	//double fmin = 5;
	//double fmax = 2048;

	double fmin = 10;
	double fmax = 2048;
	
	//double **psd = new double*[Ndetect];
	//double fmin = .006508;
	//double fmax = .0067506;
	//double fmin = .01208;
	//double fmax = 1.00;
	//double T =(t_0PN(fmin,chirpmass)- t_0PN(fmax,chirpmass));
	double T = 16;
	//std::cout<<"TIME: "<<T/T_year<<std::endl;

	double Tsignal = 4;
	//double Tsignal = 128; 
	double deltaF = 1./Tsignal;
	//Merger time -- 3/4 of total signal length. If using >4 seconds, change to Tsignal - 2
	double T_merger= Tsignal*3./4.;
	int length = (int)((fmax-fmin)/deltaF);
	params.tc=Tsignal-T_merger;
	double *frequency = new double[length];
	int Ndetect = 3;
	double **psd = new double*[Ndetect];
	
	bool AD = false;
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
		//double deltaF = (fmax-fmin)/length;	
		for(int i = 0 ; i<length; i++){
			frequency[i] = fmin + deltaF*i;
		}
	}

	std::cout<<"Freq populated"<<std::endl;

	std::string SN[3] = {"AdLIGOMidHigh","AdLIGOMidHigh","AdLIGOMidHigh"}; //"AdVIRGOPlus1"};
	for(int i = 0 ; i<Ndetect; i++){
		psd[i]= new double[length];
		populate_noise(frequency, SN[i],psd[i], length, 48);
		for(int j = 0 ; j<length; j++){
			psd[i][j]*=psd[i][j];	
		}
	}
	
	std::cout<<"frequency[10]:"<<frequency[10]<<std::endl;

	int dim = 11;	
	double* output = new double[dim];

	std::string method = "IMRPhenomD_NRT";
	std::string true_method = "IMRPhenomD";

	//std::string true_method = "ppE_IMRPhenomD_NRT_Inspiral";

	std::string detectors[3] = {"Hanford","Livingston","Virgo"};
	//std::string detector = "Hanford";

	double **output_AD = allocate_2D_array(dim,dim);
	double **output_AD_temp = allocate_2D_array(dim,dim);
	double **COV_AD = allocate_2D_array(dim,dim);
	for(int i = 0 ; i<dim; i++){
		output[i] = 0;
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

	//Get Waveforms
	
	std::complex<double> hpg[length];
	std::complex<double> hcg[length];
	std::complex<double> hpppE[length];
	std::complex<double> hcppE[length];

	std::complex<double> h_true[length];
	std::complex<double> h_model[length];
	waveform_polarizations<double> wp;

	wp.hplus = hpg;	
	wp.hcross = hcg;	
	fourier_waveform(frequency, length, &wp, method,&params);
	wp.hplus = hpppE;	
	wp.hcross = hcppE;	
	fourier_waveform(frequency, length, &wp, true_method,&params);
 


	std::cout<<"Computed fourier_waveform, file: "<<__FILE__<<" line: "<<__LINE__<<std::endl; 

	fourier_detector_response(frequency, length, h_true, detectors[0], true_method, &params);
	fourier_detector_response(frequency, length, h_model, detectors[0], method, &params);
	std::cout<<"Computed detector response, file: "<<__FILE__<<" line: "<<__LINE__<<std::endl; 

	//exit(1); 
	double* output_sys = new double[dim];
	double* output_stat = new double[dim];

	

	//Generate a plot of systematic error vs SNR, based on a range of luminosity distances

	std::ofstream output_file("data/error/sys_err_data.csv");
	std::ofstream stat_output_file("data/error/stat_err_data.csv");
	std::ofstream waveform_output_file("data/error/waveform_data.csv");
	std::ofstream bgr_waveform_output_file("data/error/bgr_waveform_data.csv");
	std::ofstream noise_curve_file("data/error/noise_curve.csv");


	
	int no_of_DL_steps = 20;
	int DL_step_size = 10.5;
	int DLeval = 40;

	double total_snr_temp = 0;

	double *phase_EA = new double[length];
	double *phase_GR = new double[length];
	  double *phase_EA_unwrap = new double[length];
	  double *phase_GR_unwrap = new double[length];
	  for(int i = 0 ; i<length ; i++){
	    phase_EA[i]= std::atan2(std::imag(h_true[i]),std::real(h_true[i]));
	    phase_GR[i]= std::atan2(std::imag(h_model[i]),std::real(h_model[i]));
	  }
	  unwrap_array(phase_EA, phase_EA_unwrap,length);
	  unwrap_array(phase_GR, phase_GR_unwrap,length);
	

	for(int i = 0; i < length; i++){
		waveform_output_file<<frequency[i]<<","<<std::real(hcg[i])<<","<<std::imag(hcg[i])<<","<<
		  //std::arg(hcg[i])
		phase_GR_unwrap[i]
		<<std::endl;
		bgr_waveform_output_file<<frequency[i]<<","<<std::real(hcppE[i])<<","<<std::imag(hcppE[i])<<","<<
		  //std::arg(hcppE[i])
		phase_EA_unwrap[i]
		<<std::endl;
		noise_curve_file<<frequency[i]<<","<<psd[1][i]<<std::endl;
	}
	exit(1);

	
	for(int a = 0; a < no_of_DL_steps; a++){
	params.Luminosity_Distance = DLeval;
	total_snr_temp = 0;

	for(int i = 0 ; i<dim; i++){
		output_sys[i] = 0;
		output_stat[i] = 0;
	}
	for(int i = 0; i < Ndetect; i++){
	total_snr_temp += pow_int( calculate_snr(SN[i],detectors[i],method, &params, frequency, length, "SIMPSONS", weights, false), 2);
	}

	wp.hplus = hpg;	
	wp.hcross = hcg;	
	fourier_waveform(frequency, length, &wp, method,&params);
	wp.hplus = hpppE;	
	wp.hcross = hcppE;	
	fourier_waveform(frequency, length, &wp, true_method,&params);

	

	calculate_systematic_error(frequency, hcg, hcppE, length, method, detectors, detectors[0], output_sys, dim, &params, 2, psd[1]);
	calculate_statistical_error(frequency, length, method, detectors, detectors[0], output_stat, dim, &params, 2, psd);

	output_file<<sqrt(total_snr_temp)<< ",";
	stat_output_file<<sqrt(total_snr_temp)<< ",";
	for(int i = 0; i < dim-1; i++){
		output_file<<output_sys[i]<<",";
		stat_output_file<<output_stat[i]<<",";
	}
	output_file<<output_sys[dim-1]<<std::endl;
	stat_output_file<<output_stat[dim-1]<<std::endl;

	DLeval+= DL_step_size;

	}
	output_file.close();
	std::vector<std::string> param_info = {"RA", "DEC", "psi", "phiRef", "tc", "iota_L", "ln DL", "ln chirpmass", "eta", "chi1", "chi2"};

	

	for(int i = 0; i < dim; i++){
	std::cout<<__LINE__<<" ? "<< param_info[i]<<" - Systematic Error "<<i<<": "<<output_sys[i]<<std::endl;
	}
	

	calculate_statistical_error(frequency, length, method, detectors, detectors[0], output_stat, dim, &params, 2, psd);
	//std::cout<<"SNR: "<<sqrt(output[0])<<std::endl;

	for(int i = 0; i < dim; i++){
	std::cout<<param_info[i]<<" - Statistical Error "<<i<<": "<<output_stat[i]<<std::endl;
	}

	for(int i = 0; i < dim; i++){
	std::cout<<param_info[i]<<" - Systematic Error / Statistical Error "<<i<<": "<<output_sys[i]/output_stat[i]<<std::endl;
	}

	

	delete [] frequency;
	delete [] weights;
	delete [] psd;
	
	return 0;
}

void RT_ERROR_MSG()
{
	std::cout<<"ERROR -- incorrect arguments"<<std::endl;
	std::cout<<"Please supply function option:"<<std::endl;
	std::cout<<"0 --- compare waveform models"<<std::endl;
	std::cout<<"1 --- test error functions"<<std::endl;
}
