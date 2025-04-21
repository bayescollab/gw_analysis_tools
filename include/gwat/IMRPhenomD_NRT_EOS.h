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
  // Function to build bump in cs2
  virtual void build_cs2_one_bump(source_parameters<T> *sp);

  // Function to convert cs2 to p(epsilon)
  virtual void cs2_to_eos_convert(); // Empty for now, "master function" that will be loaded with C++ conversion of Jaki's code

  // Function to calculate observable variables from the EOS
  virtual void get_m_love(); // Empty for now, "master function" that will be loaded with QLIMR functionality

  /* Leaving this here for now in case we need it for future development - PLEASE DELETE LATER IF NOT USED!
  
  // Inherited from IMRPhenomD_NRT which inherited from IMRPhenomD
  virtual int construct_waveform(T *frequences, int length, std::complex<T> *waveform, source)parameters<T> *params);
  */

};

#endif
