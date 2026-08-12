#ifndef GW_FISHER_H
#define GW_FISHER_H

/// @file
/// @brief GWParameterSpace-aware Fisher matrix computation.
///
/// A more modern approach to Fisher setup and calculation. This frameowrk
/// should eventually replace the old FIsher implementation.

#include <vector>

#include "quadrature.h"
#include "waveform_util.h"

/// @brief Abstract base for computing GW detector-response derivatives.
///
/// Concrete subclasses bind the waveform generation strategy and
/// finite-difference order at construction time, replacing the decision trees
/// inside the legacy calculate_derivatives function.
class GWFisherDerivatives {
 public:
  virtual ~GWFisherDerivatives() = default;

  /// @brief Compute ∂h_detector/∂θ_i for all parameters i.
  /// @param theta        MCMC parameter vector, length dim().
  /// @param ifo          Detector data; only ifo.name and ifo.freqs are used.
  /// @param[out] resp_deriv  resp_deriv[i] = ∂h/∂θ_i at each frequency.
  virtual void compute(const double* theta, const IfoData& ifo,
                       std::vector<VECCPL>& resp_deriv) const = 0;

  virtual int dim() const = 0;
  virtual const Quadrature& quadrature() const = 0;
};

/// @brief Compute the Fisher information matrix summed over all detectors.
///
/// @param theta        MCMC parameter vector, length derivatives->dim().
/// @param ifos         Detector data.
/// @param derivatives  Derivative strategy; also provides dim(), quadrature(),
///                     and the optional apply_fisher_prior().
/// @param[out] FisherM  Pre-allocated dim×dim matrix; zeroed on entry.
void fisher_numerical(const double* theta, const std::vector<IfoData>& ifos,
                      GWFisherDerivatives* derivatives, double** FisherM);

#endif  // GW_FISHER_H
