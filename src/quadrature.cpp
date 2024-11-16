#include "quadrature.h"
#include "util.h"
#include <iostream>


SimpsonsQuad::SimpsonsQuad(
    int length,     //< Length of integrand
    double delta    //< Spacing of integrand
)
{
    del = delta/3.;
    this->length = length;

    // Check if there are not enough samples
    if (length < 3)
    {
        std::cerr << "Not enough points for Simpson's rule. "
         << "Result may be inaccurate.\n";
    }

    if ((length % 2) == 0)
    {
        // If length is even, set up trapezoidal rule for the last interval
        evenFlag = true;
        SimpsonsEnd = length-2;
        TrapDel = 0.5*delta;
    }
    else
    {
        // Otherwise the whole integrand can be integrated with Simpson's
        SimpsonsEnd = length-1;
    }
}

double SimpsonsQuad::integrate(const double *integrand)
{
    // Integrand array index
    int i = 1;
    // Result of sum
    double integral = 4.*integrand[i++]; // First term of inner sum

    // Inner sum, 2*f(k) + 4*f(k+1) for even k
    while (i < SimpsonsEnd)
    {
        integral += 2.*integrand[i++];
        integral += 4.*integrand[i++];
    }

    // Add in endpoints
    integral += integrand[0] + integrand[SimpsonsEnd];
    // Multiply by overall factor
    integral *= del;

    if (evenFlag)
    {
        // Integrate final interval with trapezoidal rule.
        // Inaccurate compared to Simpson's,
        // but should be enough for small-enough spacing.
        integral += TrapDel*(integrand[SimpsonsEnd] + integrand[SimpsonsEnd+1]);
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

double SimpsonsLogQuad::integrate(const double *integrand)
{
    // Integrand weighted by the x points
    std::vector<double> wintegrand;

    for (int i=0; i < length; i++)
    {
        wintegrand.push_back(xArray.at(i)*integrand[i]);
    }

    return SimpsonsQuad::integrate(wintegrand.data());
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

    // Hold the current x value
    double x = a;
    // Hold the current log(x) value
    double logx = std::log10(x);
    // Store the point
    xPoints.push_back(x);

    while (x < b)
    {
        logx += logDelta;
        x = std::pow(10, logx);
        xPoints.push_back(x);
    }
    
    return SimpsonsLogQuad(num, logDelta, xPoints.data());
}
