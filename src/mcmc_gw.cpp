#include "mcmc_gw.h"
#include "waveform_generator.h"
#include "util.h"
#include "io_util.h"
#include "detector_util.h"
#include "ppE_utilities.h"
#include "waveform_util.h"
#include "ortho_basis.h"
#include "fisher.h"
#include "mcmc_sampler.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <fftw3.h>
#include <algorithm>
#include <iostream>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <mutex>
#include <thread>
#include <condition_variable>
#include "IMRPhenomD_NRT.h" //For testing purposes only! Remove later.
#include "EA_IMRPhenomD_NRT.h" //For testing purposes only! Remove later.

/*! \file 
 * Routines for implementation in MCMC algorithms specific to GW CBC analysis
 *
 * */





//######################################################################################
//######################################################################################
