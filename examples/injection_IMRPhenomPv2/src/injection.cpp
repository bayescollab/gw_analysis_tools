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

void spin_spherical_coord(double spinx, double spiny, double spinz, double *add_mag, double *add_ctheta, double *add_phi );

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
	injection.mass1 = 8.0;
	injection.mass2 = 12.0; //These two together give a chirpmass of 1.188 (like in GW170817)

	double chirpmass = calculate_chirpmass(injection.mass1,injection.mass2);
	double eta = calculate_eta(injection.mass1,injection.mass2);
	double q = injection.mass2/injection.mass1;

	//Spins
	injection.spin1[2] = .1; //spin_1z
	injection.spin2[2] = -.15; //spin_2z
	injection.spin1[1] = 0.0; //spin_1y
	injection.spin2[1] = 0.0;//spin_2y
	injection.spin1[0] = 0.9;//spin_1x
	injection.spin2[0] = 0.0;//spin_2x

	//Convert spins to spherical coordinates
	double a1,a2,ctheta1,ctheta2,phi1,phi2;
	spin_spherical_coord(injection.spin1[0],injection.spin1[1], injection.spin1[2], &a1, &ctheta1, &phi1 );
	spin_spherical_coord(injection.spin2[0],injection.spin2[1], injection.spin2[2], &a2, &ctheta2, &phi2 );
	std::cout<<"spherical coordinates for spin1 (a1,ctheta1,phi1): "<<a1<<" ,"<<ctheta1<<" ,"<<phi1<<std::endl;
	std::cout<<"spherical coordinates for spin2 (a2,ctheta2,phi2): "<<a2<<" ,"<<ctheta2<<" ,"<<phi2<<std::endl;
	//Extrinsic
	injection.Luminosity_Distance = 451.189 ; //std::atof(argv[3]);// Mpc
	injection.psi = 2.659;
	injection.RA = 1.375;
	injection.DEC = -1.2108;
	injection.incl_angle = M_PI/4.0; // default is 0.51
	double gps = 1126259642.413;
	injection.gmst = gps_to_GMST_radian(gps);
	std::cout<<"GMST: "<<injection.gmst<<std::endl;
	//Trivial constants
	injection.phiRef = 1.3;
	injection.f_ref = 50.; // the reference frequency at which the spins are defined.

	//Specific flags -- probably don't modify
	injection.equatorial_orientation = false;
	injection.horizon_coord = false;
	injection.shift_time = true;
	injection.shift_phase = true;

	//Set injection string -- for injections that are read in, this will be ignored.
	std::string injection_method = "IMRPhenomPv2";

	//Detector properties
	int detect_number = 4;
	std::string detectors[detect_number] = {"Hanford","Livingston","Virgo","Kagra"};
	std::string SN[detect_number] = {"AdLIGODesign","AdLIGODesign","AdLIGODesign","KAGRA_pess"};
	//Frequency array properties
	//Length of total signal -- 
	double Tsignal = 16.; //16 or 4
	double deltaF = 1./Tsignal;
	double fmin = 10.0;
	double fmax =2048.0+deltaF; //1024.0+deltaF or 2049.0
	//Merger time -- 3/4 of total signal length. If using >4 seconds, change to Tsignal - 2. RC:where does this come from?
	double T_merger=Tsignal-2;
	int length = (int)((fmax-fmin)/deltaF);

	//Input PSD, freq, and data arrays -- these will actually be passed into the MCMC
	std::complex<double> **data = new std::complex<double>*[detect_number];
	double **psd = new double*[detect_number];
	double **freq = new double*[detect_number];

	//Set merger time of injection
	injection.tc=Tsignal-T_merger;
	std::cout<<"tc: "<<injection.tc<<std::endl;
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

	int out_rows = detect_number+1; //Number of cols in output array
	double **data_out = new double*[out_rows]; //Declare output array
	for(int i = 0 ; i<out_rows; i++){
		data_out[i] = new double[length];
	}
	/*
	// *******************************
	//		Write out data
	// *******************************
	for(int i = 0 ; i<length; i++){
		data_out[0][i] = freq[0][i] ; //Write first freq set, since they're identical
		int ct=0;
		for(int j = 1 ; j<out_rows; j+=1){ //change detector_number to out_rows
			data_out[j][i] = psd[ct][i];
			ct+=1;
		}
	}
	write_file("data/detectors_PSD_updated_Tsignal32.csv",data_out, out_rows, length);
	for(int i = 0 ; i<out_rows; i++){
		delete [] data_out[i];
	}
	delete [] data_out;
*/
	//Create coherent gravitational wave detector between all detectors
	create_coherent_GW_detection(detectors,detect_number, freq, data_lengths, true, &injection, injection_method, data);

	//Calculate the SNR and print
	double total_snr=0;
	for(int i = 0 ; i<detect_number; i++){
		total_snr += pow_int( calculate_snr_internal(psd[i],data[i],freq[i],length, "SIMPSONS",(double*) NULL, false), 2);
	}
	std::cout<<"NETWORK SNR of injection: "<<sqrt(total_snr)<<std::endl;

	// #####################################################################################################################
		//WRITE DATA OUT FOR TROUBLESHOOTING
	// #####################################################################################################################
		int data_out_rows = detect_number*2+1; //Number of cols in output array
		double **data_out_new = new double*[data_out_rows]; //Declare output array
		for(int i = 0 ; i<data_out_rows; i++){
			data_out_new[i] = new double[length];
		}
		//Write out data
		for(int i = 0 ; i<length; i++){
			data_out_new[0][i] = freq[0][i] ; //Write first freq set, since they're identical
			int ct = 0 ;
			for(int j = 1 ; j<data_out_rows; j+=2){ //change detector_number to out_rows
				data_out_new[j][i] = std::real(data[ct][i]);
				data_out_new[j+1][i] = std::imag(data[ct][i]);
				ct+=1;
			}
		}
		//Write out freq, real(response_1), imag(reasponse_1) ...
		write_file("data/responses_debugging.csv",data_out_new, data_out_rows, length);

	//###################################################################
	//###################################################################
	//MCMC Settings
	//###################################################################
	//###################################################################

	int dim = 15;
	std::string recovery_method = "IMRPhenomPv2";
	// std::string recovery_method = "IMRPhenomD";
	int ensembleSize = 20; //temperature ladder -- ensembleSize gives the number of temperature values.
	int ensembleN = 6; // number of chains per temperature

	//Initial Guess  -- also the exact right answer.. Perks of injections!!
	//Needs to match dimension/model we're using -- e.g. if using IMRPhenomD_NRT with tidal love, 12 dimensions. If not, 13.
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
	initial_position[9]= a1;
	initial_position[10]= a2;
	initial_position[11]= ctheta1;
	initial_position[12]= ctheta2;
	initial_position[13]= phi1;
	initial_position[14]= phi2;
	//Write it out to file
	write_file("data/injections.csv",initial_position,dim);
	std::cout<<"Written to file "<<std::endl;
	bayesship::positionInfo initialPosition(dim, false);
	for(int i = 0 ; i<dim ; i++){
		initialPosition.parameters[i] = initial_position[i];
	}

	//MCMC parameters
	//Frequency of swaps -- every N steps attempt a swap
	int swapProb = .5;
	//Number of threads to use
	int threads = 4;
	//Whether to pool or not, leave at true
	bool pool = true;

	//Additional modifications
	MCMC_modification_struct mod_struct;
	mod_struct.ppE_Nmod = 0;
	//double b = -5; // b=-7 is -1PN, b=-5 is 0PN and so on.
	//mod_struct.bppe =&b ;

	priorData PD;
	PD.mass1_prior[0] =5.0 ;
	PD.mass1_prior[1] = 50.0;
	PD.mass2_prior[0] = 5.0;
	PD.mass2_prior[1] = 50.0;
	PD.a1_prior[0]=0.0;
	PD.a1_prior[1]=1.0;
	PD.a2_prior[0]=0.0;
	PD.a2_prior[1]=1.0;	
	PD.ctheta1_prior[0]=-1.0;
	PD.ctheta1_prior[1]=1.0;
	PD.ctheta2_prior[0]=-1.0;
	PD.ctheta2_prior[1]=1.0;
	PD.phi1_prior[0]=0;
	PD.phi1_prior[1]=2*M_PI;
	PD.phi2_prior[0]=0;
	PD.phi2_prior[1]=2*M_PI;
	PD.DL_prior[0] = 50;
	PD.DL_prior[1] = 800;
	PD.RA_bounds[0] = 0;
	PD.RA_bounds[1] = 2*M_PI;
	PD.sinDEC_bounds[0] = -1;
	PD.sinDEC_bounds[1] = 1;
	PD.T_merger = T_merger;
	
	//allocate memory for ppE priors
	//PD.mod_priors = new double*[1]; //The argument here should be [Nmod]
	//PD.mod_priors[0] = new double[2]; //If there are additional params, we can replace this with a for loop that assigns lower and upper
																		// bounds for each ppE mod.
	// //Assign prior range one ppE modification at a time
	//PD.mod_priors[0][0]=-3./128;
	//PD.mod_priors[0][1]=3./128;

	std::cout<<"Priors set"<<std::endl;

	//Initialize standard prior -- see in gw_analysis_tools/src/standardPriorLibrary.cpp
	bayesship::probabilityFn *log_prior = new logPriorStandard_P(&PD);;
	std::cout<<"Prior library defined"<<std::endl;

	//Samples, burn in parameters
	int independentSamples = 2000;
	int burnIterations = 5000;
	int burnPriorIterations = 0;
	int priorIterations = 0;
	bool writePriorData = true;
	int batchsize = 1e3;
	std::string outputDir("data/");
	std::string outputMoniker("IMRPhenomPv2_injection_IMRPhenomPv2_recovery");
	bool ignoreExistingCheckpoint = false;
	bool  restrictSwapTemperature= false;
	bool  coldChainStorageOnly= true;


	/* Test initial point for prior compliance */
	double lp = log_prior->eval(&initialPosition, 0);
	std::cout<<"Prior of initial point: "<<lp<<std::endl;

	bayesship::bayesshipSampler *sampler = PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW_v2(dim, independentSamples,
		ensembleSize, ensembleN, &initialPosition, (bayesship::positionInfo**)nullptr,
		swapProb, burnIterations, burnPriorIterations, priorIterations, writePriorData, batchsize, (double **) nullptr,
		log_prior, threads, pool, detect_number, data, psd, freq, data_lengths,
		gps, detectors, &mod_struct, recovery_method, outputDir, outputMoniker, ignoreExistingCheckpoint,
		restrictSwapTemperature, coldChainStorageOnly);

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
	//delete [] PD.mod_priors[0];
	//delete [] PD.mod_priors;

	return 0;
}

void spin_spherical_coord(double spinx, double spiny, double spinz, double *add_mag, double *add_ctheta, double *add_phi ){
	*add_mag=sqrt(pow(spinx,2)+pow(spiny,2)+pow(spinz,2));
	if(spiny==0.0 ){
		*add_phi = 0.0;
	}
	else{
		*add_phi = atan(spiny/spinx);
	}
	*add_ctheta= spinz/ *add_mag;
}

