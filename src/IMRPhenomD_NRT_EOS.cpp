#include "IMRPhenomD_NRT_EOS.h"
#include "IMRPhenomD_NRT.h"
#include <math.h>
#include <adolc/adouble.h>
#include <adolc/taping.h>
#include <adolc/drivers/drivers.h>
#include <iostream>
#include <cmath>
#include <complex>
#include "util.h"
#include <gsl/gsl_randist.h> 

/*! \file 
 * File for the generation of equations of state with features in the speed of
 * sound squared as in arXiv:2106.03890 and related papers.
 *
 * Supported features: single Gaussian, ...[add more later]
 *
 * Supported crust types: SLy EOS
 * 
 * Also included in this file ...
 */

template<class T> 
void IMRPhenomD_NRT_EOS<T>::build_cs2_one_bump(source_parameters<T> *sp)
{
  //Will need to add some sort of array/table/pointer thing to pass back and forth the values of cs2
  
  /* A function to construct cs^2 by stitching together an EOS in the crust,
   * a single feature (Gaussian bump), and a plateau.
   *
   * The bump is described by the height/magnitude (bump_mag),
   * the width (bump_width), and the location (bump_offset).
   */

  
  std::cout<<"Calling build_cs2_one_bump"<<std::endl; 
  
}


template class IMRPhenomD_NRT_EOS<double>;
template class IMRPhenomD_NRT_EOS<adouble>; 
