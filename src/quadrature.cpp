#include "quadrature.h"
#include "util.h"
#include <iostream>
#include <algorithm>


// Error message for when Simpson's rule is not ideal
const std::string SIMPSONS_ERANGE_MESSAGE {"Not enough points for Simpson's rule. Result may be inaccurate."};

SimpsonsQuad::SimpsonsQuad(
    int length,     //< Length of integrand
    double delta    //< Spacing of integrand
)
{
  del = 4.0 * delta / 3.;
  this->length = length;

  // Check if there are not enough samples
  if (length < 3) {
    std::cerr << SIMPSONS_ERANGE_MESSAGE << "\n";
  }

    if ((length % 2) == 0)
    {
        // If length is even, set up trapezoidal rule for the last interval
        evenFlag = true;
        SimpsonsEnd = length-2;
        TrapDel =
            2.0 * delta;  // 4(inner product factor) * 0.5(trapezoid factor)
    }
    else
    {
        // Otherwise the whole integrand can be integrated with Simpson's
        SimpsonsEnd = length-1;
    }
}

double SimpsonsQuad::inner_product(const std::complex<double>* h1,
                                   const std::complex<double>* h2,
                                   const double* psd) const {
  // Inner sum: odd indices get weight 4, even indices get weight 2.
  double integral = 0.;
  for (int i = 1; i < SimpsonsEnd; i++)
    integral +=
        (i % 2 == 1 ? 4. : 2.) * std::real(h1[i] * std::conj(h2[i])) / psd[i];

  // Add endpoints and multiply by overall factor
  for (int i : {0, SimpsonsEnd})
    integral += std::real(h1[i] * std::conj(h2[i])) / psd[i];
  integral *= del;

  if (evenFlag) {
    // Integrate final interval with trapezoidal rule.
    for (int i : {SimpsonsEnd, SimpsonsEnd + 1})
      integral += TrapDel * std::real(h1[i] * std::conj(h2[i])) / psd[i];
  }

    return integral;
}

//! Create a SimpsonsQuad given the endpoints of the integral and number of 
//! points desired. Return the points in x as well.
SimpsonsQuad CreateSimpsonsQuad(
    std::vector<double> &xPoints, //< [out] Vector of generated points in x
    double a,   //< Start of integral
    double b,   //< End of integral
    int num     //< Number of points to be sampled
)
{
    // Calculate the spacing
    double delta = (b-a)/(num-1);
    // Hold the current x value to be stored
    double currentx = a;
    // Store the first point
    xPoints.push_back(currentx);

    while (currentx < b)
    {
        currentx += delta;
        xPoints.push_back(currentx);
    }

    return SimpsonsQuad(num, delta);
}

SimpsonsLogQuad::SimpsonsLogQuad(
    int length,         //< Length of integrand
    double logdelta,    //< Logarithmic spacing of integrand
    const double *xArray      //< Array of points in the integrand
)
// Hand SimpsonsQuad the spacing as ln(10)*Δ(log x)
 : SimpsonsQuad(length, LOG10*logdelta)
{
    // Store the xArray points
    // Copy points to this->xArray
    for (int i=0; i < length; i++)
    {
        this->xArray.push_back(xArray[i]);
    }

    // Shrink the capacity if possible
    this->xArray.shrink_to_fit();
}

// Use SimpsonsQuad method by weighting the PSD array
double SimpsonsLogQuad::inner_product(const std::complex<double>* h1,
                                      const std::complex<double>* h2,
                                      const double* psd) const {
  std::vector<double> weighted_psd;

  for (int i = 0; i < length; i++) {
    weighted_psd.push_back(psd[i] / xArray.at(i));
  }

  return SimpsonsQuad::inner_product(h1, h2, weighted_psd.data());
}

//! Create a SimpsonsLogQuad given the endpoints of the integral and number of 
//! points desired. Return the points in x as well.
SimpsonsLogQuad CreateSimpsonsLogQuad(
    std::vector<double> &xPoints,
    double a,
    double b,
    int num
)
{
    // Calculate the spacing in log-space
    double logDelta = (std::log10(b) - std::log10(a)) / (num-1);
    std::cout << "log-spacing = " << logDelta << ".\n";

    // Hold the current x value
    double x = a;
    // Store the initial log(x) value
    double loga = std::log10(x);
    // Hold the current log(x) value
    double logx = loga;
    // Store the point
    xPoints.push_back(x);

    for (int i = 1; i < num; i++)
    {
        logx = loga + i*logDelta;
        x = std::pow(10, logx);
        xPoints.push_back(x);
    }
    
    return SimpsonsLogQuad(num, logDelta, xPoints.data());
}
