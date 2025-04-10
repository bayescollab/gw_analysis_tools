#ifndef IMRPHENOMD_NRT_EOS
#define IMRPHENOMD_NRT_EOS
#include "IMRPhenomD_NRT.h"
#include "util.h"

/*! \file
 */

/*! Class that extends the IMRPhenomD_NRT waveform to sample directly on
 * equation of state (EOS) parameters.
 */

template<class T>
class IMRPhenomD_NRT_EOS: public IMRPhenomD_NRT<T>
{
public:
  virtual void build_cs2_one_bump(source_parameters<T> *sp);

};

#endif
