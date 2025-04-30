#ifndef IMRPHENOMD_NRT_EA_H
#define IMRPHENOMD_NRT_EA_H
#include "IMRPhenomD.h"
#include "IMRPhenomD_NRT.h"
#include "util.h"
#include "waveform_generator.h"

template<class T>
class EA_IMRPhenomD_NRT: public IMRPhenomD_NRT<T>
{
public:
  virtual T calculate_EA_sensitivity(int body, source_parameters<T> *p);

  virtual void EA_check_nan(source_parameters<T> *p);

  virtual void pre_calculate_EA_factors(source_parameters<T> *p);

  virtual T EA_phase_ins1(T f, useful_powers<T> *powers, source_parameters<T> *p);

  virtual T EA_phase_ins2(T f, useful_powers<T> *powers, source_parameters<T> *p);

  virtual T EA_amp_ins1(T f, useful_powers<T> *powers, source_parameters<T> *p);

  virtual T EA_amp_ins2(T f, useful_powers<T> *powers, source_parameters<T> *p);
  //virtual double EA_amp_ins2(double f, useful_powers<double> *powers, source_parameters<double> *p);
  //virtual adouble EA_amp_ins2(adouble f, useful_powers<adouble> *powers, source_parameters<adouble> *p);

  virtual int EA_construct_waveform(T *frequencies, int length, waveform_polarizations<T> *waveform, source_parameters<T> *params);

  virtual int construct_phase(T *frequencies, int length, T *phase, source_parameters<T> *params);

  virtual int construct_amplitude(T *frequencies, int length, T *amplitude, source_parameters<T> *params);

  virtual int construct_waveform(T *frequencies, int length, std::complex<T> *waveform, source_parameters<T> *params);
};

#endif
