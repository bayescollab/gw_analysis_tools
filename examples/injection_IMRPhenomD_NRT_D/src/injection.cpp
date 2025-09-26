#include <iostream>
#include <gwat/waveform_util.h>
#include <gwat/waveform_generator.h>
#include <gwat/util.h>
#include <gwat/detector_util.h>
#include <gwat/io_util.h>
#include <gwat/mcmc_gw_extended.h>
#include <gwat/mcmc_sampler.h>
#include <gwat/standardPriorLibrary.h>
#include <bayesship/dataUtilities.h>
#include <bayesship/utilities.h>
#include <bayesship/bayesshipSampler.h>


int main(int argc, char *argv[])
{
	//######################################
	//######################################
	//######################################
	//Injection Parameters
	//######################################
	//######################################
	//######################################

	gen_params injection;

	//Masses in M_SOL, plus useful other parameters
	injection.mass1 = 1.44;
	injection.mass2 = 1.29399; //These two together give a chirpmass of 1.188 (like in GW170817)
	//injection.mass1 = 1.81;
	//injection.mass2 = 1.514;
	double chirpmass = calculate_chirpmass(injection.mass1,injection.mass2);
	double eta = calculate_eta(injection.mass1,injection.mass2);
	double q = injection.mass2/injection.mass1;

	//Spins
	injection.spin1[2] = .003;
	injection.spin2[2] = -.002;
	injection.spin1[1] = .0;
	injection.spin2[1] = .0;
	injection.spin1[0] = .0;
	injection.spin2[0] = .0;

	//Extrinsic
	injection.Luminosity_Distance = 63;
	injection.psi = .2;
        injection.RA = 3.42;
	injection.DEC = -.37;
	injection.incl_angle = 2.532207345558998;
	double gps = 1187008882.4;
	injection.gmst = gps_to_GMST_radian(gps);

	//Tidal parameters -- either use tidal love or not
	//injection.tidal1 = 100;	
	//injection.tidal2 = 100;	
	//injection.tidal_s = .5*(injection.tidal1 + injection.tidal2);
	injection.tidal_love = true;
	injection.tidal_s = 242;
	
        injection.diss_tidal1 = 100;	
	injection.diss_tidal2 = 100;	
	
	//Trivial constants
	injection.phiRef = 2.;
	injection.f_ref = 20.;

	//Specific flags -- probably don't modify
	injection.equatorial_orientation = false;
	injection.horizon_coord = false;
	injection.shift_time = true;
	injection.shift_phase = true;

	//PhenomD_NRT injection
	std::string injection_method = "IMRPhenomD_NRT";
	
	//Detector properties
	int detect_number = 3;
	std::string detectors[3] = {"Hanford","Livingston","Virgo"};
	std::string SN[3] = {"AdLIGOMidHigh","AdLIGOMidHigh","AdVIRGOPlus1"};
	//Frequency array properties
	double fmin = 10;
	double fmax =2048;
	//Length of total signal -- leave at 4 seconds for now
	double Tsignal = 4;
	double deltaF = 1./Tsignal;
	//Merger time -- 3/4 of total signal length. If using >4 seconds, change to Tsignal - 2
	double T_merger= Tsignal*3./4.;
	int length = (int)((fmax-fmin)/deltaF);
	
	//Input PSD, freq, and data arrays -- these will actually be passed into the MCMC
	std::complex<double> **data = new std::complex<double>*[detect_number];
	double **psd = new double*[detect_number];
	double **freq = new double*[detect_number];

	//Set merger time of injection
	injection.tc=Tsignal-T_merger;
	
	//Populate freq, psd, and allocate array for data
	for(int i = 0 ; i<detect_number; i++){
		data[i]= new std::complex<double>[length];
		psd[i]= new double[length];
		freq[i]= new double[length];
		//Populate freq array
		for(int j =0 ; j<length; j++){
			freq[i][j] = fmin + j*deltaF;
		}
		//Populate PSD array
		populate_noise(freq[i],SN[i],psd[i],length);
		for(int j =0 ; j<length; j++){
			psd[i][j] *= psd[i][j];
		}
	}

	//Set data_length array 
	int data_lengths[detect_number];
	for(int i = 0 ; i<detect_number; i++){
		data_lengths[i]=length;	
	}

	//Create coherent gravitational wave detector between all detectors
	create_coherent_GW_detection(detectors,detect_number, freq, data_lengths,  true, &injection, injection_method, data);
	
	//Calculate the SNR and print
	double total_snr=0;
	for(int i = 0 ; i<detect_number; i++){
		total_snr += pow_int( calculate_snr_internal(psd[i],data[i],freq[i],length, "SIMPSONS",(double*) NULL, false), 2);
	}
	std::cout<<"NETWORK SNR of injection: "<<sqrt(total_snr)<<std::endl;
	//###################################################################
	//###################################################################
	//MCMC Settings
	//###################################################################
	//###################################################################

	int dim = 14;
	std::string recovery_method = "IMRPhenomD_NRT";
	int ensembleSize = 20;
	int ensembleN = 5;
	
	//Initial Guess  -- also the exact right answer.. Perks of injections!!
	//Needs to match dimension/model we're using -- e.g. if using IMRPhenomD_NRT with tidal love, 14 dimensions. If not, 15.
	double initial_position[dim]; 
	initial_position[0]= injection.RA; 
	initial_position[1]= sin(injection.DEC);
	initial_position[2]= injection.psi; 
	initial_position[3]= cos(injection.incl_angle); 
	initial_position[4]= injection.phiRef;
	initial_position[5]= T_merger;
	initial_position[6]= log(injection.Luminosity_Distance);
	initial_position[7]= log(chirpmass);
	initial_position[8]= eta;
	initial_position[9]= injection.spin1[2];
	initial_position[10]= injection.spin2[2];
	initial_position[11]= log(injection.tidal_s);
	initial_position[12]= injection.diss_tidal1;
	initial_position[13]= injection.diss_tidal2;
	  
	//Write it out to file
	write_file("data/injections.csv",initial_position,dim);

	bayesship::positionInfo initialPosition(dim, false);
	for(int i = 0 ; i<dim ; i++){
		initialPosition.parameters[i] = initial_position[i];
	}	


	//MCMC parameters
	//Frequency of swaps -- every N steps attempt a swap
	int swapProb = .5;
	//Number of threads to use
	//int threads = 64;
	int threads = 16;
	//Whether to pool or not, leave at true
	bool pool = true;

	//Additional modifications
	MCMC_modification_struct mod_struct;
	mod_struct.ppE_Nmod = 0;
	//double b = -1;
	//mod_struct.bppe =&b ;

	//Use tidal love or not
	bool tidal_love = true;
	mod_struct.tidal_love = tidal_love;
        bool tidal_love_error = false;
	mod_struct.tidal_love_error = tidal_love_error; 

	priorData PD;
	PD.mass1_prior[0] =.5 ;
	PD.mass1_prior[1] = 3;
	PD.mass2_prior[0] = .5;
	PD.mass2_prior[1] = 3;
	PD.spin1_prior[0] = -.01;
	PD.spin1_prior[1] = .01;
	PD.spin2_prior[0] = -.01;
	PD.spin2_prior[1] = .01;
	PD.DL_prior[0] = 1;
	PD.DL_prior[1] = 500;
	PD.tidal_s_prior[0] = 1;
	PD.tidal_s_prior[1] = 5000;
	PD.diss_tidal1_prior[0] = 0;
	PD.diss_tidal1_prior[1] = 1000;
	PD.diss_tidal2_prior[0] = 0;
	PD.diss_tidal2_prior[1] = 1000;
	PD.RA_bounds[0] = 0;
	PD.RA_bounds[1] = 2*M_PI;
	PD.sinDEC_bounds[0] = -1;
	PD.sinDEC_bounds[1] = 1;
	PD.T_merger = T_merger;
	PD.tidal_love = tidal_love;
	PD.tidal_love_error = tidal_love_error;

	//Initialize standard prior -- see in gw_analysis_tools/src/standardPriorLibrary.cpp
	bayesship::probabilityFn *log_prior = new logPriorStandard_D_NRT(&PD);;

	//Samples, burn in parameters 	
	int independentSamples = 10000;
	int burnIterations = 50000;
	int burnPriorIterations = 0;
	int priorIterations = 0;
	bool writePriorData = true;
	int batchsize = 1e4;
	std::string outputDir("data/");
	std::string outputMoniker("EA_injection");
	bool ignoreExistingCheckpoint = false;
	bool  restrictSwapTemperature= false;
	bool  coldChainStorageOnly= true;



	/* Test initial point for prior compliance */
	double lp = log_prior->eval(&initialPosition, 0);
	std::cout<<"Prior of initial point: "<<lp<<std::endl;	

	bayesship::bayesshipSampler *sampler = PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW_v2(dim, independentSamples, ensembleSize, ensembleN, &initialPosition, (bayesship::positionInfo**)nullptr,
		swapProb, burnIterations, burnPriorIterations, priorIterations, writePriorData, batchsize, (double **) nullptr, log_prior, threads, pool, detect_number, data, psd, freq, data_lengths, 
		gps, detectors, &mod_struct, recovery_method, outputDir, outputMoniker, ignoreExistingCheckpoint, restrictSwapTemperature, coldChainStorageOnly);	
	
	//Output structures

	delete log_prior;
	
	for(int i = 0 ; i<detect_number; i++){
		delete [] freq[i];	
		delete [] psd[i];	
		delete [] data[i];	
	}
	delete [] injection.bppe;
	delete [] injection.betappe;
	delete [] freq;
	delete [] psd;
	delete [] data;

	return 0;
}
