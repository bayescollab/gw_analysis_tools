#include <iostream>
#include <gwat/fisher.h>
#include <gwat/util.h>
#include <gwat/io_util.h>
#include <gwat/ortho_basis.h>
#include <gwat/detector_util.h>
#include <gwat/IMRPhenomP.h>

/*! \file 
 *
 * Command line tool for calculating a Fisher Information Matrix (FIM) for a GW event
 *
 * Usage:
 * 	
 * 	fisher_tool /PATH/TO/PARAM/FILE
 *
 * See data/sample_config_files/fisher_param_template.dat for an example parameter file
 *
 */
int help();
int main(int argc, char *argv[])
{
	//std::cout<<"Fisher Information Matrix Tool"<<std::endl;
	std::cout.precision(15);
	if(argc !=2){
		std::cout<<"ERROR -- A parameter file is required"<<std::endl;
		return 1;
	}
	if(std::string(argv[1]) == "help"){
		return help();
	}
	std::string param_file(argv[1]);

	std::unordered_map<std::string, int> int_dict;
	std::unordered_map<std::string, double> dbl_dict;
	std::unordered_map<std::string, float> flt_dict;
	std::unordered_map<std::string, bool> bool_dict;
	std::unordered_map<std::string, std::string> str_dict;
	std::string input_param_file(argv[1]);
	int status = unpack_input_io_file(input_param_file, &int_dict, &str_dict, &dbl_dict, &flt_dict, &bool_dict);

	int dim = int_dict["Dimension"];
	std::string gen_method = str_dict["Generation Method"];
	int N_d = int_dict["Number of Detectors"];
	std::string detectors[N_d];
	std::string PSDs[N_d];
	//std::cout<<"Using Detectors|Noise curves: ";
	for( int i = 0 ; i<N_d ; i++){
		detectors[i] = str_dict["Detector Name "+std::to_string(i)];
		PSDs[i] = str_dict["PSD "+std::to_string(i)];
		//std::cout<<detectors[i]<<"|"<<PSDs[i]<<" ";
	}
	//std::cout<<std::endl;
	int Nmod = 0; double *bppe = NULL; double *betappe = NULL;
	if(has_substring(gen_method, "ppE")){
		Nmod = int_dict["Number of Modifications"];	
		bppe = new double[Nmod];
		betappe = new double[Nmod];
		for(int i = 0 ; i<Nmod ; i++){
			bppe[i]= dbl_dict["b ppE "+std::to_string(i)];	
			betappe[i]= dbl_dict["beta ppE "+std::to_string(i)];	
		}
	}

	gen_params params;
	params.gmst = dbl_dict["GMST"];
	params.Nmod = Nmod;
	params.bppe = bppe;
	params.betappe = betappe;
	params.mass1 = dbl_dict["Mass 1"];
	params.mass2 = dbl_dict["Mass 2"];
	params.Luminosity_Distance = dbl_dict["Distance"];
	params.spin1[0] = dbl_dict["SPIN1X"];
	params.spin1[1] = dbl_dict["SPIN1Y"];
	params.spin1[2] = dbl_dict["SPIN1Z"];
	params.spin2[0] = dbl_dict["SPIN2X"];
	params.spin2[1] = dbl_dict["SPIN2Y"];
	params.spin2[2] = dbl_dict["SPIN2Z"];
	params.f_ref = 20;
	params.phiRef = dbl_dict["reference phase"];
	params.tc = dbl_dict["t"];
	params.RA = dbl_dict["RA"];
	params.DEC = dbl_dict["DEC"];
	params.psi = dbl_dict["PSI"];
	params.incl_angle = dbl_dict["IOTA"];
	params.shift_time = false;
	params.shift_phase = false;
	params.horizon_coord = false;
	params.NSflag1 = false;
	params.NSflag2 = false;
	params.sky_average = false;
	params.equatorial_orientation = false;

	IMRPhenomPv2<double> model;
	params.chip = model.PhenomPv2_inplane_spin(&params);
	params.phip = 0;

	int length = 1000;
	double **psds = allocate_2D_array(N_d,length);
	double *freqs =new  double[length];
	double *weights =new  double[length];
	double **fisher_temp = allocate_2D_array(dim,dim);
	double **fisher = allocate_2D_array(dim,dim);
	for(int j = 0 ; j<dim; j++){
		for(int k = 0 ; k<dim; k++){
			fisher[j][k] = 0;
		}
	}
	
	double fmin = 10;
	double fmax = 2048;
	
	gauleg(log10(fmin), log10(fmax),freqs, weights,length);	
	for(int i =0 ; i<length ; i++){
		freqs[i]=pow(10,freqs[i]);
	}
	for(int i = 0 ; i<N_d; i++){
		populate_noise(freqs, PSDs[i],psds[i],length,48);
		for(int j = 0 ; j<length; j++){
			psds[i][j]*=psds[i][j];
		}
	}
	for(int i = 0 ; i<N_d ; i++){
		fisher_autodiff(freqs, length, gen_method, detectors[i],detectors[0],fisher_temp, dim, &params, "GAUSSLEG",weights,true, psds[i],NULL,NULL);
		for(int j = 0 ; j<dim; j++){
			for(int k = 0 ; k<dim; k++){
				fisher[j][k] += fisher_temp[j][k];
			}
		}
	}
	
	write_file(str_dict["Output File Path"], fisher, dim,dim);

	if(has_substring(gen_method, "ppE")){
		delete [] bppe;
		delete [] betappe;
	}
	delete [] freqs;
	delete [] weights;
	deallocate_2D_array(psds, N_d,length);
	deallocate_2D_array(fisher, dim,dim);
	deallocate_2D_array(fisher_temp, dim,dim);
	return 0;
}
int help()
{
	std::cout<<"Fisher Information Matrix Command Line Tool:"<<std::endl;
	std::cout<<"Use the template in data/sample_config_files "<<std::endl;
	std::cout<<"to model your parameter file"<<std::endl;
	std::cout<<std::endl;
	std::cout<<"Waveform options: "<<std::endl;
	std::cout<<"IMRPhenomPv2: "<<std::endl;
	std::cout<<"12 dimension: RA, DEC, PSI, IOTA, phiRef, tc, ln(Luminosity Distance), "<<std::endl;
	std::cout<<"ln (Chirpmass), eta, aligned spin 1, aligned spin 2, In-plane effective spin"<<std::endl;
	std::cout<<std::endl;
	std::cout<<"IMRPhenomD: "<<std::endl;
	std::cout<<"11 dimension: RA, DEC, PSI, IOTA, phiRef, tc, ln(Luminosity Distance), "<<std::endl;
	std::cout<<"ln (Chirpmass), eta, aligned spin 1, aligned spin 2"<<std::endl;
	std::cout<<std::endl;
	std::cout<<"ppE:"<<std::endl;
	std::cout<<"ppE_(base)_IMR: Modification in full waveform"<<std::endl;
	std::cout<<"ppE_(base)_Inspiral: Modification in inspiral waveform"<<std::endl;
	std::cout<<"Dimension: base dimension + Number of modifications"<<std::endl;
	std::cout<<std::endl;
	std::cout<<"Detector and noise curve options -- See src/detector_util.cpp for latest detector options"<<std::endl;
	std::cout<<"Current detector options: Hanford, Livingston, Virgo, Indigo, Kagra, CE"<<std::endl;
	return 0;
}
