#include "gwat/util.h"
#include "gwat/fisher.h"
#include "gwat/detector_util.h"
#include <iostream>



int AD_v_N(int argc, char *argv[]);
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
	else{
		RT_ERROR_MSG();
		return 1;
	}
}
int AD_v_N(int argc, char *argv[])
{
	std::cout.precision(5);
	gen_params params;	
	params.spin1[2] = .3;
	params.spin2[2] = .3;
	params.chip = .7;
	params.phip = 0.1;
	params.Luminosity_Distance = 100;
	params.phic = 1;
	params.RA = 2.;
	params.DEC = -1.1;
	params.f_ref = 1e-5;
	params.NSflag1 = false;
	params.NSflag2 = false;
	params.horizon_coord = false;
	params.shift_time=false;
	params.shift_phase=false;
	
	//params.mass1 = 5e6;
	//params.mass2 = 5e6;
	//params.f_ref = 20;
	//params.psi = .1;
	//double fmin = 20;
	//double fmax = 2048;
	//params.incl_angle = .1;
	//params.tc = 10;
	//params.equatorial_orientation = false;
	//double T = 32;
	params.mass1 = 9e1;
	params.mass2 = 4e1;
	params.theta_l = 1;
	params.phi_l = 2;
	params.tc = T_year;
	params.equatorial_orientation = true;
	double fmin = 3e-2;
	double fmax = 1e-1;
	double T = T_year/2;

	int length = T*(fmax-fmin);
	double *frequency = new double[length];
	double *psd = new double[length];
	
	for(int i = 0 ; i<length; i++){
		frequency[i]=fmin + (double)i /T;
	}
	//populate_noise(frequency, "aLIGO_analytic",psd, length, 48);
	populate_noise(frequency, "LISA_CONF",psd, length, 12);
	for(int i = 0 ; i<length; i++){
		psd[i]*=psd[i];	
	}

	std::string detector = "LISA";
	std::string method = "IMRPhenomPv2";
		
	int dim = 13;
	double **output_N = allocate_2D_array(dim,dim);
	double **output_AD = allocate_2D_array(dim,dim);
	double **COV_AD = allocate_2D_array(dim,dim);
	double **output_AD2 = allocate_2D_array(dim,dim);

	fisher_numerical(frequency, length, method, detector, output_N, dim, &params, 2, NULL,NULL, psd);
	
	std::cout<<"NUMERICAL:"<<std::endl;
	for(int i = 0 ; i<dim; i++){
		std::cout<<i<<" ";
		for(int j = 0 ; j<dim; j++){
			std::cout<<output_N[i][j]<<" ";
		}
		std::cout<<std::endl;
	}

	fisher_autodiff(frequency, length, method, detector, output_AD, dim, &params, "SIMPSONS",NULL,false, psd,NULL,NULL);
	
	gsl_LU_matrix_invert(output_AD,COV_AD,dim);
	
	std::cout<<"AD:"<<std::endl;
	for(int i = 0 ; i<dim; i++){
		std::cout<<i<<" ";
		for(int j = 0 ; j<dim; j++){
			std::cout<<output_AD[i][j]<<" ";
		}
		std::cout<<std::endl;
	}
	std::cout<<"COV AD:"<<std::endl;
	for(int i = 0 ; i<dim; i++){
		std::cout<<i<<" ";
		for(int j = 0 ; j<dim; j++){
			std::cout<<COV_AD[i][j]<<" ";
		}
		std::cout<<std::endl;
	}

	//params.phiRef-=1;
	//params.phip+=1;
	//fisher_autodiff(frequency, length, method, detector, output_AD2, dim, &params, "SIMPSONS",NULL,false, psd,NULL,NULL);
	//
	//std::cout<<"FD AD / AD2:"<<std::endl;
	//for(int i = 0 ; i<dim; i++){
	//	std::cout<<i<<" ";
	//	for(int j = 0 ; j<dim; j++){
	//		std::cout<<(output_AD2[i][j] - output_AD[i][j])*2./(output_AD2[i][j] + output_AD[i][j])<<" ";
	//	}
	//	std::cout<<std::endl;
	//}
	std::cout<<"FRACTIONAL DIFF (N-AD)*2/(N+AD):"<<std::endl;
	for(int i = 0 ; i<dim; i++){

		for(int j = 0 ; j<dim; j++){
			std::cout<<(output_N[i][j] - output_AD[i][j])*2./(output_N[i][j] + output_AD[i][j])<<" ";
		}
		std::cout<<std::endl;
	}

	deallocate_2D_array(output_AD,dim,dim);
	deallocate_2D_array(COV_AD,dim,dim);
	deallocate_2D_array(output_AD2,dim,dim);
	deallocate_2D_array(output_N,dim,dim);
	
	delete [] frequency;
	delete [] psd;
	return 0;
}
void RT_ERROR_MSG()
{
	std::cout<<"ERROR -- incorrect arguments"<<std::endl;
	std::cout<<"Please supply function option:"<<std::endl;
	std::cout<<"0 --- Compare AD to Numerical"<<std::endl;
}

