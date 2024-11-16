#ifndef QUADRATURE_H
#define QUADRATURE_H


#include <vector>


//! \class Quadrature
//! \brief Class to evaluate integrals with established spacing and weights.
class Quadrature
{
protected:
    // Length of integrand
    int length = 0;

public:
    virtual ~Quadrature() = default;

    virtual double integrate(const double *integrand) = 0;
    virtual int get_length() {return length;}
};

//! \brief Simpson's rule for uniformly-spaced integrals.
//!
//! Quadrature with the classic extended 3-point rule
//! (see, e.g., Numerical Recipes, extended Simpsons rule).
//! For even lengths, the trapezoidal rule is used at the last interval.
class SimpsonsQuad : public Quadrature
{
private:
    double del;         //< Overall factor h/3
    bool evenFlag = false;  //< Track if the length is even
    int SimpsonsEnd;    //< End index of integrand array for Simpson's
    double TrapDel;     //< Overall factor for trapezoid h/2

public:
    SimpsonsQuad(int length, double delta);

    virtual double integrate(const double *integrand);
};

SimpsonsQuad CreateSimpsonsQuad(
    std::vector<double> &xPoints,
    double a,
    double b,
    int num
);

//! \brief Simpson's rule for logarithmically uniformly-spaced integrals.
//!
//! Quadrature with the classic extended 3-point rule
//! (see, e.g., Numerical Recipes, extended Simpsons rule).
//! For even lengths, the trapezoidal rule is used at the last interval.
//! 
class SimpsonsLogQuad : public SimpsonsQuad
{
private:
    std::vector<double> xArray;
    // double *xArray;     //< Array of points in the integrand

public:
    SimpsonsLogQuad(int length, double logdelta, const double *xArray);
    double integrate(const double *integrand) override;
};

SimpsonsLogQuad CreateSimpsonsLogQuad(
    std::vector<double> &xPoints,
    double a,
    double b,
    int num
);


#endif // QUADRATURE_H