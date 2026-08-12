#include "quadrature.h"

#include <algorithm>
#include <iostream>
#include <memory>

#include "detector_util.h"
#include "waveform_util.h"

// ----- UTILITIES -----

int unrefine_uniform_grid_with_trapezoid_quad(
    const int N_initial, const double f_min, const double f_max,
    std::vector<std::string>& psds, const int N_detectors,
    std::string* detectors, gen_params_base<double>* params, std::string model,
    const double tol, const bool log_spacing, const int max_iterations) {
  std::cout << "Halving from N_initial=" << N_initial
            << " to find N_converged (tol=" << tol << ") ...\n";

  const double fmin = log_spacing ? std::log10(f_min) : f_min;
  const double range = (log_spacing ? std::log10(f_max) : f_max) - fmin;
  std::vector<std::vector<double>> nvs(N_detectors);
  std::vector<std::complex<double>*> data(N_detectors, nullptr);
  std::vector<double*> freqs(N_detectors, nullptr);
  std::vector<int> lengths(N_detectors, 0);

  auto eval_hh = [&](int N) -> double {
    const double dx = range / (N - 1);
    std::vector<double> xv(N);
    for (int j = 0; j < N; j++) {
      if (log_spacing)
        xv[j] = std::pow(10.0, fmin + j * dx);
      else
        xv[j] = fmin + j * dx;
    }
    xv[N - 1] = f_max;

    for (int d = 0; d < N_detectors; d += 1) {
      nvs[d].assign(N, 0.0);
      data[d] = new std::complex<double>[N];
      freqs[d] = xv.data();
      lengths[d] = N;

      populate_noise(xv.data(), psds[d], nvs[d].data(), N);
      for (int j = 0; j < N; j++) nvs[d][j] *= nvs[d][j];
    }

    if (N_detectors == 1) {
      create_single_GW_detection(data[0], detectors[0], xv.data(), N, params,
                                 model);
    } else {
      create_coherent_GW_detection(detectors, N_detectors, freqs.data(),
                                   lengths.data(), true, params, model,
                                   data.data());
    }

    std::unique_ptr<Quadrature> q;
    if (log_spacing)
      q = std::make_unique<TrapezoidLogQuad>(N, dx, xv.data());
    else
      q = std::make_unique<TrapezoidQuad>(N, dx);

    double total_hh = 0.0;
    for (int i = 0; i < N_detectors; i++) {
      VECCPL di(data[i], data[i] + N);
      total_hh += q->inner_product(di, di, nvs[i]);
    }

    for (int d = 0; d < N_detectors; d += 1) delete[] data[d];

    return total_hh;
  };

  const double hh_ref = eval_hh(N_initial);
  int N_converged = N_initial;
  int N_halve = N_initial;
  int counter = 0;
  while (counter <= max_iterations && N_halve > 3) {
    int N_half = (N_halve - 1) / 2 + 1;
    if (N_half < 3) N_half = 3;
    double rel = std::abs(eval_hh(N_half) - hh_ref) / std::abs(hh_ref);
    if (rel < tol) {
      N_converged = N_half;
      counter++;
    } else
      break;
    N_halve = N_half;
    if (N_halve <= 3) break;
  }
  std::cout << "N_converged = " << N_converged << "  (" << N_initial << " / "
            << N_converged << " = "
            << static_cast<double>(N_initial) / N_converged << "x reduction)\n";

  return N_converged;
}

// ----- QUADRATURE CLASSES -----
// Error message for when Simpson's rule is not ideal
const std::string SIMPSONS_ERANGE_MESSAGE{
    "Not enough points for Simpson's rule. Result may be inaccurate."};

SimpsonsQuad::SimpsonsQuad(int length,   //< Length of integrand
                           double delta  //< Spacing of integrand
) {
  del = 4.0 * delta / 3.;
  this->length = length;

  // Check if there are not enough samples
  if (length < 3) {
    std::cerr << SIMPSONS_ERANGE_MESSAGE << "\n";
  }

  if ((length % 2) == 0) {
    // If length is even, set up trapezoidal rule for the last interval
    evenFlag = true;
    SimpsonsEnd = length - 2;
    TrapDel = 2.0 * delta;  // 4(inner product factor) * 0.5(trapezoid factor)
  } else {
    // Otherwise the whole integrand can be integrated with Simpson's
    SimpsonsEnd = length - 1;
  }
}

double SimpsonsQuad::inner_product(const VECCPL& h1, const VECCPL& h2,
                                   const VECDBL& psd) const {
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
    std::vector<double>& xPoints,  //< [out] Vector of generated points in x
    double a,                      //< Start of integral
    double b,                      //< End of integral
    int num                        //< Number of points to be sampled
) {
  // Calculate the spacing
  double delta = (b - a) / (num - 1);
  // Hold the current x value to be stored
  double currentx = a;
  // Store the first point
  xPoints.push_back(currentx);

  while (currentx < b) {
    currentx += delta;
    xPoints.push_back(currentx);
  }

  return SimpsonsQuad(num, delta);
}

SimpsonsLogQuad::SimpsonsLogQuad(
    int length,           //< Length of integrand
    double logdelta,      //< Logarithmic spacing of integrand
    const double* xArray  //< Array of points in the integrand
    )
    // Hand SimpsonsQuad the spacing as ln(10)*Δ(log x)
    : SimpsonsQuad(length, LOG10 * logdelta), xArray(xArray, xArray + length) {}

// Use SimpsonsQuad method by weighting the PSD array
double SimpsonsLogQuad::inner_product(const VECCPL& h1, const VECCPL& h2,
                                      const VECDBL& psd) const {
  VECDBL weighted_psd(length);
  for (int i = 0; i < length; i++)
    weighted_psd[i] = psd[i] / xArray[i];
  return SimpsonsQuad::inner_product(h1, h2, weighted_psd);
}

//! Create a SimpsonsLogQuad given the endpoints of the integral and number of
//! points desired. Return the points in x as well.
SimpsonsLogQuad CreateSimpsonsLogQuad(std::vector<double>& xPoints, double a,
                                      double b, int num) {
  // Calculate the spacing in log-space
  double logDelta = (std::log10(b) - std::log10(a)) / (num - 1);
  std::cout << "log-spacing = " << logDelta << ".\n";

  // Hold the current x value
  double x = a;
  // Store the initial log(x) value
  double loga = std::log10(x);
  // Hold the current log(x) value
  double logx = loga;
  // Store the point
  xPoints.push_back(x);

  for (int i = 1; i < num; i++) {
    logx = loga + i * logDelta;
    x = std::pow(10, logx);
    xPoints.push_back(x);
  }

  return SimpsonsLogQuad(num, logDelta, xPoints.data());
}

TrapezoidQuad::TrapezoidQuad(int length,   //< Length of integrand
                             double delta  //< Spacing of integrand
) {
  weight = 4.0 * delta;
  this->length = length;
}

double TrapezoidQuad::inner_product(const VECCPL& h1, const VECCPL& h2,
                                    const VECDBL& psd) const {
  // Inner sum
  double integral = 0.;
  for (int i = 1; i < length - 1; i++)
    integral += std::real(h1[i] * std::conj(h2[i])) / psd[i];

  // Add endpoints and multiply by overall factor
  for (int i : {0, length - 1})
    integral += 0.5 * std::real(h1[i] * std::conj(h2[i])) / psd[i];
  integral *= weight;

  return integral;
}

TrapezoidLogQuad::TrapezoidLogQuad(
    int length,           //< Length of integrand
    double logdelta,      //< Logarithmic spacing of integrand
    const double* xArray  //< Array of points in the integrand
    )
    // Hand TrapezoidQuad the spacing as ln(10)*Δ(log x)
    : TrapezoidQuad(length, LOG10 * logdelta),
      xArray(xArray, xArray + length) {}

// Use TrapezoidQuad method by weighting the PSD array
double TrapezoidLogQuad::inner_product(const VECCPL& h1, const VECCPL& h2,
                                       const VECDBL& psd) const {
  VECDBL weighted_psd(length);
  for (int i = 0; i < length; i++)
    weighted_psd[i] = psd[i] / xArray[i];
  return TrapezoidQuad::inner_product(h1, h2, weighted_psd);
}
