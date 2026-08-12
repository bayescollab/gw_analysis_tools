#ifndef QUADRATURE_H
#define QUADRATURE_H

#include <complex>
#include <string>
#include <vector>

#include "util.h"

// ----- UTILITIES -----

int unrefine_uniform_grid_with_trapezoid_quad(
    int N_initial, double f_min, double f_max, std::vector<std::string>& psds,
    int N_detectors, std::string* detectors, gen_params_base<double>* params,
    std::string model, double tol, bool log_spacing, int max_iterations = 10);

// ----- QUADRATURE CLASSES -----

//! \class Quadrature
//! \brief Class to evaluate integrals with established spacing and weights.
class Quadrature {
 protected:
  // Length of integrand
  int length = 0;

 public:
  virtual ~Quadrature() = default;

  /// @brief Compute the inner product 4•∫ h1(f)•h2*(f) / psd(f) df.
  /// @param h1
  /// @param h2
  /// @param psd
  virtual double inner_product(const VECCPL& h1, const VECCPL& h2,
                               const VECDBL& psd) const = 0;
  virtual int get_length() const { return length; }
};

//! \brief Simpson's rule for uniformly-spaced integrals.
//!
//! Quadrature with the classic extended 3-point rule
//! (see, e.g., Numerical Recipes, extended Simpsons rule).
//! For even lengths, the trapezoidal rule is used at the last interval.
class SimpsonsQuad : public Quadrature {
 private:
  double del;             //< Overall factor h/3
  bool evenFlag = false;  //< Track if the length is even
  int SimpsonsEnd;        //< End index of integrand array for Simpson's
  double TrapDel;         //< Overall factor for trapezoid h/2

 public:
  SimpsonsQuad(int length, double delta);

  virtual double inner_product(const VECCPL& h1, const VECCPL& h2,
                               const VECDBL& psd) const;
};

SimpsonsQuad CreateSimpsonsQuad(std::vector<double>& xPoints, double a,
                                double b, int num);

//! \brief Simpson's rule for logarithmically uniformly-spaced integrals.
//!
//! Quadrature with the classic extended 3-point rule
//! (see, e.g., Numerical Recipes, extended Simpsons rule).
//! For even lengths, the trapezoidal rule is used at the last interval.
//!
class SimpsonsLogQuad : public SimpsonsQuad {
 private:
  std::vector<double> xArray;
  // double *xArray;     //< Array of points in the integrand

 public:
  SimpsonsLogQuad(int length, double logdelta, const double* xArray);
  double inner_product(const VECCPL& h1, const VECCPL& h2,
                       const VECDBL& psd) const override;
};

SimpsonsLogQuad CreateSimpsonsLogQuad(std::vector<double>& xPoints, double a,
                                      double b, int num);

//! \brief Trapezoidal rule for uniformly-spaced integrals.
class TrapezoidQuad : public Quadrature {
 private:
  double weight;  //< Overall factor 4*h

 public:
  TrapezoidQuad(int length, double delta);

  virtual double inner_product(const VECCPL& h1, const VECCPL& h2,
                               const VECDBL& psd) const;
};

//! \brief Trapezoidal rule for logarithmically uniformly-spaced integrals.
class TrapezoidLogQuad : public TrapezoidQuad {
 private:
  /// @brief Array of values of the dependent variable.
  std::vector<double> xArray;

 public:
  TrapezoidLogQuad(int length, double logdelta, const double* xArray);
  double inner_product(const VECCPL& h1, const VECCPL& h2,
                       const VECDBL& psd) const override;
};

#endif  // QUADRATURE_H