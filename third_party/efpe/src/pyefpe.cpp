#include "pyefpe/pyefpe.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <thread>
#include <tuple>
#include <utility>

namespace pyefpe {
namespace {

// Implementation map:
//
// 1. Convert physical input parameters into the dimensionless state used by the
//    pyEFPE radiation-reaction equations.
// 2. Integrate that state with RK45 and store dense-output segments.
// 3. For each segment, decide which eccentric Fourier harmonics have support.
// 4. For each requested frequency, solve the stationary-phase condition
//    d Phi_mode / dt = 2 pi f and evaluate the SUA amplitude there.
//
// The code is intentionally self-contained: no Boost/GSL dependency is required,
// and the public API only exposes Parameters, Model, and waveform containers.

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kEulerGamma = 0.577215664901532860606512090082402431;

// Geometric-unit conversions. Internally masses and distances are expressed as
// times in seconds, so the phasing formulas can be used directly.
constexpr double kTSunSeconds = 4.92549094831e-6;
constexpr double kMpcSeconds = 1.02927125054339e14;
const double kY2Prefactor = 1.0 / std::sqrt(0.8 * kPi);

using C = std::complex<double>;

// Radiation-reaction state vector:
//   [0] y              PN velocity-like variable
//   [1] e^2            squared eccentricity
//   [2] lambda         mean longitude-like orbital phase
//   [3] delta_lambda   pericenter-related phase difference
//   [4] DeltaJ^2       precession invariant correction
//   [5] psi_p_bar      normalized precession phase
//   [6] phi_z0         Euler angle integration constant
//   [7] zeta0          Euler angle integration constant
using State = std::array<double, 8>;

struct Vec3 {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
};

Vec3 operator+(const Vec3& a, const Vec3& b) {
  return {a.x + b.x, a.y + b.y, a.z + b.z};
}

Vec3 operator*(double a, const Vec3& v) {
  return {a * v.x, a * v.y, a * v.z};
}

double dot(const Vec3& a, const Vec3& b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec3 cross(const Vec3& a, const Vec3& b) {
  return {a.y * b.z - a.z * b.y,
          a.z * b.x - a.x * b.z,
          a.x * b.y - a.y * b.x};
}

double norm(const Vec3& v) {
  return std::sqrt(dot(v, v));
}

Vec3 normalized(Vec3 v) {
  const double n = norm(v);
  if (!(n > 0.0)) {
    throw std::runtime_error("cannot normalize zero vector");
  }
  return (1.0 / n) * v;
}

double clamp(double x, double lo, double hi) {
  return std::max(lo, std::min(hi, x));
}

double safe_sqrt(double x) {
  return std::sqrt(std::max(0.0, x));
}

double sqr(double x) {
  return x * x;
}

struct MassParams {
  double M = 0.0;
  double mu1 = 0.0;
  double mu2 = 0.0;
  double nu = 0.0;
  double dmu = 0.0;
};

// Dimensionless mass combinations used throughout the PN formulas. At this
// point m1 and m2 are already in seconds, not solar masses.
MassParams mass_params_from_m1_m2(double m1, double m2) {
  const double M = m1 + m2;
  const double mu1 = m1 / M;
  const double mu2 = m2 / M;
  return {M, mu1, mu2, mu1 * mu2, mu1 - mu2};
}

struct SphericalTrig {
  double sin_theta = 0.0;
  double cos_theta = 0.0;
  double sin_phi = 0.0;
  double cos_phi = 0.0;
};

SphericalTrig spherical_trig_from_vec(Vec3 v) {
  v = normalized(v);
  const double cth = clamp(v.z, -1.0, 1.0);
  const double ph = std::atan2(v.y, v.x);
  return {safe_sqrt(1.0 - cth * cth), cth, std::sin(ph), std::cos(ph)};
}

std::pair<double, double> spherical_coords_from_vec(Vec3 v) {
  v = normalized(v);
  return {clamp(v.z, -1.0, 1.0), std::atan2(v.y, v.x)};
}

Vec3 matvec(const std::array<std::array<double, 3>, 3>& m, const Vec3& v) {
  return {m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
          m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
          m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z};
}

class Poly {
public:
  Poly() = default;
  explicit Poly(std::vector<double> coefficients) : coefficients_(std::move(coefficients)) {}
  Poly(std::initializer_list<double> coefficients) : coefficients_(coefficients) {}

  double operator()(double x) const {
    if (coefficients_.empty()) {
      return 0.0;
    }
    // Horner evaluation keeps the many PN polynomials compact and stable.
    double result = coefficients_.back();
    for (std::size_t i = coefficients_.size() - 1; i-- > 0;) {
      result = coefficients_[i] + result * x;
    }
    return result;
  }

private:
  std::vector<double> coefficients_;
};

std::vector<double> vec(std::initializer_list<double> values) {
  return std::vector<double>(values);
}

std::vector<double> scale(double a, std::vector<double> v) {
  for (double& x : v) {
    x *= a;
  }
  return v;
}

std::vector<double> add(std::vector<double> a, const std::vector<double>& b) {
  if (a.size() < b.size()) {
    a.resize(b.size(), 0.0);
  }
  for (std::size_t i = 0; i < b.size(); ++i) {
    a[i] += b[i];
  }
  return a;
}

double bessel_j(int n, double x) {
  if (n >= 0) {
    return ::jn(n, x);
  }
  const int an = -n;
  const double sign = (an % 2 == 0) ? 1.0 : -1.0;
  return sign * ::jn(an, x);
}

// Carlson symmetric elliptic integrals. These small routines avoid a hard
// dependency on Boost/GSL while covering the parameter range used by pyEFPE.
double carlson_rc(double x, double y) {
  if (x < 0.0 || y == 0.0) {
    throw std::domain_error("invalid arguments to Carlson RC");
  }
  constexpr double errtol = 0.0012;
  constexpr double c1 = 0.3;
  constexpr double c2 = 1.0 / 7.0;
  constexpr double c3 = 0.375;
  constexpr double c4 = 9.0 / 22.0;

  double xt = x;
  double yt = y;
  double w = 1.0;
  if (yt < 0.0) {
    xt = x - y;
    yt = -y;
    w = std::sqrt(x) / std::sqrt(xt);
  }

  double ave = 0.0;
  double s = 0.0;
  for (int iter = 0; iter < 100; ++iter) {
    const double sx = std::sqrt(xt);
    const double sy = std::sqrt(yt);
    const double lambda = 2.0 * sx * sy + yt;
    xt = 0.25 * (xt + lambda);
    yt = 0.25 * (yt + lambda);
    ave = (xt + yt + yt) / 3.0;
    s = (yt - ave) / ave;
    if (std::abs(s) < errtol) {
      const double poly = 1.0 + s * s * (c1 + s * (c2 + s * (c3 + s * c4)));
      return w * poly / std::sqrt(ave);
    }
  }
  throw std::runtime_error("Carlson RC did not converge");
}

double carlson_rf(double x, double y, double z) {
  if (x < 0.0 || y < 0.0 || z < 0.0) {
    throw std::domain_error("invalid arguments to Carlson RF");
  }
  constexpr double errtol = 0.0025;
  constexpr double c1 = 1.0 / 24.0;
  constexpr double c2 = 3.0 / 44.0;
  constexpr double c3 = 1.0 / 14.0;

  double xt = x;
  double yt = y;
  double zt = z;
  double ave = 0.0;
  double delx = 0.0;
  double dely = 0.0;
  double delz = 0.0;
  for (int iter = 0; iter < 100; ++iter) {
    const double sx = std::sqrt(xt);
    const double sy = std::sqrt(yt);
    const double sz = std::sqrt(zt);
    const double lambda = sx * (sy + sz) + sy * sz;
    xt = 0.25 * (xt + lambda);
    yt = 0.25 * (yt + lambda);
    zt = 0.25 * (zt + lambda);
    ave = (xt + yt + zt) / 3.0;
    delx = (ave - xt) / ave;
    dely = (ave - yt) / ave;
    delz = (ave - zt) / ave;
    if (std::max({std::abs(delx), std::abs(dely), std::abs(delz)}) < errtol) {
      const double e2 = delx * dely - delz * delz;
      const double e3 = delx * dely * delz;
      return (1.0 + (c1 * e2 - 0.1 - c2 * e3) * e2 + c3 * e3) / std::sqrt(ave);
    }
  }
  throw std::runtime_error("Carlson RF did not converge");
}

double carlson_rd(double x, double y, double z) {
  if (x < 0.0 || y < 0.0 || z <= 0.0) {
    throw std::domain_error("invalid arguments to Carlson RD");
  }
  constexpr double errtol = 0.0015;
  constexpr double c1 = 3.0 / 14.0;
  constexpr double c2 = 1.0 / 6.0;
  constexpr double c3 = 9.0 / 22.0;
  constexpr double c4 = 3.0 / 26.0;
  constexpr double c5 = 0.25 * c3;
  constexpr double c6 = 1.5 * c4;

  double xt = x;
  double yt = y;
  double zt = z;
  double sum = 0.0;
  double fac = 1.0;
  double ave = 0.0;
  double delx = 0.0;
  double dely = 0.0;
  double delz = 0.0;
  for (int iter = 0; iter < 100; ++iter) {
    const double sx = std::sqrt(xt);
    const double sy = std::sqrt(yt);
    const double sz = std::sqrt(zt);
    const double lambda = sx * (sy + sz) + sy * sz;
    sum += fac / (sz * (zt + lambda));
    fac *= 0.25;
    xt = 0.25 * (xt + lambda);
    yt = 0.25 * (yt + lambda);
    zt = 0.25 * (zt + lambda);
    ave = 0.2 * (xt + yt + 3.0 * zt);
    delx = (ave - xt) / ave;
    dely = (ave - yt) / ave;
    delz = (ave - zt) / ave;
    if (std::max({std::abs(delx), std::abs(dely), std::abs(delz)}) < errtol) {
      const double ea = delx * dely;
      const double eb = delz * delz;
      const double ec = ea - eb;
      const double ed = ea - 6.0 * eb;
      const double ee = ed + ec + ec;
      const double poly =
          1.0 + ed * (-c1 + c5 * ed - c6 * delz * ee) +
          delz * (c2 * ee + delz * (-c3 * ec + delz * c4 * ea));
      return 3.0 * sum + fac * poly / (ave * std::sqrt(ave));
    }
  }
  throw std::runtime_error("Carlson RD did not converge");
}

double carlson_rj(double x, double y, double z, double p) {
  if (x < 0.0 || y < 0.0 || z < 0.0 || p == 0.0) {
    throw std::domain_error("invalid arguments to Carlson RJ");
  }
  constexpr double errtol = 0.0015;
  constexpr double c1 = 3.0 / 14.0;
  constexpr double c2 = 1.0 / 3.0;
  constexpr double c3 = 3.0 / 22.0;
  constexpr double c4 = 3.0 / 26.0;
  constexpr double c5 = 0.75 * c3;
  constexpr double c6 = 1.5 * c4;
  constexpr double c7 = 0.5 * c2;
  constexpr double c8 = 2.0 * c3;

  double xt = x;
  double yt = y;
  double zt = z;
  double pt = p;
  double sum = 0.0;
  double fac = 1.0;
  double ave = 0.0;
  double delx = 0.0;
  double dely = 0.0;
  double delz = 0.0;
  double delp = 0.0;

  for (int iter = 0; iter < 100; ++iter) {
    const double sx = std::sqrt(xt);
    const double sy = std::sqrt(yt);
    const double sz = std::sqrt(zt);
    const double lambda = sx * (sy + sz) + sy * sz;
    double alpha = pt * (sx + sy + sz) + sx * sy * sz;
    alpha *= alpha;
    const double beta = pt * sqr(pt + lambda);
    sum += fac * carlson_rc(alpha, beta);
    fac *= 0.25;
    xt = 0.25 * (xt + lambda);
    yt = 0.25 * (yt + lambda);
    zt = 0.25 * (zt + lambda);
    pt = 0.25 * (pt + lambda);
    ave = 0.2 * (xt + yt + zt + 2.0 * pt);
    delx = (ave - xt) / ave;
    dely = (ave - yt) / ave;
    delz = (ave - zt) / ave;
    delp = (ave - pt) / ave;
    if (std::max({std::abs(delx), std::abs(dely), std::abs(delz), std::abs(delp)}) < errtol) {
      const double ea = delx * (dely + delz) + dely * delz;
      const double eb = delx * dely * delz;
      const double ec = delp * delp;
      const double ed = ea - 3.0 * ec;
      const double ee = eb + 2.0 * delp * (ea - ec);
      const double poly =
          1.0 + ed * (-c1 + c5 * ed - c6 * ee) +
          eb * (c7 + delp * (-c8 + delp * c4)) +
          delp * ea * (c2 - delp * c3) - c2 * delp * ec;
      return 3.0 * sum + fac * poly / (ave * std::sqrt(ave));
    }
  }
  throw std::runtime_error("Carlson RJ did not converge");
}

double ellipk(double m) {
  m = clamp(m, 0.0, 1.0 - 1.0e-14);
  return carlson_rf(0.0, 1.0 - m, 1.0);
}

double ellipe(double m) {
  m = clamp(m, 0.0, 1.0 - 1.0e-14);
  return carlson_rf(0.0, 1.0 - m, 1.0) - (m / 3.0) * carlson_rd(0.0, 1.0 - m, 1.0);
}

double ellipf_inc(double phi, double m) {
  m = clamp(m, 0.0, 1.0 - 1.0e-14);
  const double s = std::sin(phi);
  const double c2 = sqr(std::cos(phi));
  return s * carlson_rf(c2, 1.0 - m * s * s, 1.0);
}

double ellipe_inc(double phi, double m) {
  m = clamp(m, 0.0, 1.0 - 1.0e-14);
  const double s = std::sin(phi);
  const double s2 = s * s;
  const double c2 = sqr(std::cos(phi));
  return s * carlson_rf(c2, 1.0 - m * s2, 1.0) -
         (m * s2 * s / 3.0) * carlson_rd(c2, 1.0 - m * s2, 1.0);
}

double ellippi(double n, double m) {
  m = clamp(m, 0.0, 1.0 - 1.0e-14);
  return ellipk(m) + (n / 3.0) * carlson_rj(0.0, 1.0 - m, 1.0, 1.0 - n);
}

double ellippi_inc(double n, double phi, double m) {
  m = clamp(m, 0.0, 1.0 - 1.0e-14);
  const double s = std::sin(phi);
  const double s2 = s * s;
  const double ns2 = n * s2;
  const double c2 = sqr(std::cos(phi));
  return s * carlson_rf(c2, 1.0 - m * s2, 1.0) +
         (ns2 * s / 3.0) * carlson_rj(c2, 1.0 - m * s2, 1.0, 1.0 - ns2);
}

struct Jacobi {
  double sn = 0.0;
  double cn = 1.0;
  double dn = 1.0;
  double am = 0.0;
};

Jacobi ellipj_from_u(double u, double m) {
  m = clamp(m, 0.0, 1.0 - 1.0e-14);
  double phi = clamp(u, -0.5 * kPi, 0.5 * kPi);
  for (int i = 0; i < 16; ++i) {
    const double s = std::sin(phi);
    const double delta = u - ellipf_inc(phi, m);
    phi += delta * std::sqrt(std::max(1.0e-30, 1.0 - m * s * s));
    phi = clamp(phi, -0.5 * kPi, 0.5 * kPi);
    if (std::abs(delta) < 1.0e-13 * std::max(1.0, std::abs(u))) {
      break;
    }
  }
  const double sn = std::sin(phi);
  const double cn = std::cos(phi);
  return {sn, cn, safe_sqrt(1.0 - m * sn * sn), phi};
}

std::pair<std::array<C, 5>, std::array<C, 5>> compute_necessary_Wigner_small_d2(double cth) {
  cth = clamp(cth, -1.0, 1.0);
  const double cth2 = cth * cth;
  const double sth = safe_sqrt(1.0 - cth2);
  std::array<C, 5> d2_mp2 = {
      C{sqr(0.5 * (1.0 - cth)), 0.0},
      C{0.5 * sth * (1.0 - cth), 0.0},
      C{std::sqrt(0.375) * sth * sth, 0.0},
      C{0.5 * sth * (1.0 + cth), 0.0},
      C{sqr(0.5 * (1.0 + cth)), 0.0}};

  const double d2_10 = -std::sqrt(1.5) * sth * cth;
  std::array<C, 5> d2_mp0 = {
      d2_mp2[2],
      C{-d2_10, 0.0},
      C{0.5 * (3.0 * cth2 - 1.0), 0.0},
      C{d2_10, 0.0},
      d2_mp2[2]};
  return {d2_mp2, d2_mp0};
}

std::pair<std::array<C, 5>, std::array<C, 5>> compute_necessary_Wigner_D2(
    double phi, double cos_theta, double zeta) {
  auto [D2_mp2, D2_mp0] = compute_necessary_Wigner_small_d2(cos_theta);
  const C exp_miph = std::polar(1.0, -phi);
  const std::array<C, 2> exp_mimpphi = {exp_miph, exp_miph * exp_miph};

  D2_mp2[0] *= std::conj(exp_mimpphi[1]);
  D2_mp2[1] *= std::conj(exp_mimpphi[0]);
  D2_mp2[3] *= exp_mimpphi[0];
  D2_mp2[4] *= exp_mimpphi[1];
  const C exp_m2iz = std::polar(1.0, -2.0 * zeta);
  for (C& v : D2_mp2) {
    v *= exp_m2iz;
  }

  D2_mp0[3] *= exp_mimpphi[0];
  D2_mp0[4] *= exp_mimpphi[1];
  D2_mp0[0] = std::conj(D2_mp0[4]);
  D2_mp0[1] = -std::conj(D2_mp0[3]);
  return {D2_mp2, D2_mp0};
}

std::array<C, 5> compute_m2_Y2(double cos_theta, double phi) {
  auto d = compute_necessary_Wigner_D2(phi, cos_theta, 0.0).first;
  for (C& v : d) {
    v = kY2Prefactor * std::conj(v);
  }
  return d;
}

double N20_Newtonian(int p, double e2) {
  if (p == 0) {
    return 0.0;
  }
  return std::sqrt(2.0 / 3.0) * bessel_j(p, static_cast<double>(p) * safe_sqrt(e2));
}

double N22_Newtonian(int p, double e2) {
  const double e = safe_sqrt(e2);
  const double sq1me2 = safe_sqrt(1.0 - e2);
  const double Jpm2_fact = 0.5 * (sq1me2 + 1.0 - 0.5 * e2);
  const double Jpp2_fact = -0.0625 * e2 * e2 / Jpm2_fact;
  const int k = p + 2;
  const double ke = static_cast<double>(k) * e;
  return static_cast<double>(k) *
         (-sq1me2 * bessel_j(k, ke) +
          0.5 * e * (bessel_j(k + 1, ke) - bessel_j(k - 1, ke)) +
          Jpm2_fact * bessel_j(k - 2, ke) + Jpp2_fact * bessel_j(k + 2, ke));
}

struct ModeOrder {
  int m_abs = 0;
  int p = 0;
};

// Decide which Newtonian eccentric harmonics to keep in a segment. The Python
// implementation builds many (m, p) contributions; here we rank candidates by
// their Newtonian power and stop once the cumulative omitted power is below tol.
// This is the central mode-count/accuracy tradeoff for eccentric waveforms.
std::vector<ModeOrder> Newtonian_orders_needed(double e2, int pmax, double tol) {
  struct Candidate {
    ModeOrder mode;
    double weight = 0.0;
  };
  std::vector<Candidate> candidates;
  candidates.reserve(static_cast<std::size_t>(3 * pmax + 1));
  for (int p = 1; p <= pmax; ++p) {
    const double n = N20_Newtonian(p, e2);
    candidates.push_back({{0, p}, n * n});
  }
  for (int p = -pmax; p <= pmax; ++p) {
    const double n = N22_Newtonian(p, e2);
    candidates.push_back({{2, p}, n * n});
  }
  const double norm = (4.0 / 3.0) * (-1.0 + 4.0 / safe_sqrt(1.0 - e2));
  for (auto& c : candidates) {
    c.weight /= norm;
  }
  std::sort(candidates.begin(), candidates.end(), [](const Candidate& a, const Candidate& b) {
    return a.weight > b.weight;
  });

  std::vector<ModeOrder> result;
  double cumulative = 0.0;
  for (const auto& c : candidates) {
    result.push_back(c.mode);
    cumulative += c.weight;
    if (cumulative > 1.0 - tol) {
      break;
    }
  }
  return result;
}

std::vector<C> solve_linear_complex(std::vector<std::vector<C>> A, std::vector<C> b) {
  const std::size_t n = b.size();
  for (std::size_t i = 0; i < n; ++i) {
    std::size_t pivot = i;
    double pivot_abs = std::abs(A[i][i]);
    for (std::size_t r = i + 1; r < n; ++r) {
      const double v = std::abs(A[r][i]);
      if (v > pivot_abs) {
        pivot_abs = v;
        pivot = r;
      }
    }
    if (pivot_abs == 0.0) {
      throw std::runtime_error("singular complex linear system");
    }
    if (pivot != i) {
      std::swap(A[i], A[pivot]);
      std::swap(b[i], b[pivot]);
    }
    const C diag = A[i][i];
    for (std::size_t c = i; c < n; ++c) {
      A[i][c] /= diag;
    }
    b[i] /= diag;
    for (std::size_t r = 0; r < n; ++r) {
      if (r == i) {
        continue;
      }
      const C factor = A[r][i];
      for (std::size_t c = i; c < n; ++c) {
        A[r][c] -= factor * A[i][c];
      }
      b[r] -= factor * b[i];
    }
  }
  return b;
}

// SUA replaces the slowly varying amplitude by a short symmetric stencil around
// the stationary time. These coefficients are independent of the binary and are
// found by solving the small linear system implied by the SUA matching
// conditions for the chosen kmax.
std::vector<C> compute_ak_SUA(int kmax) {
  std::vector<std::vector<C>> M(static_cast<std::size_t>(kmax + 1),
                                std::vector<C>(static_cast<std::size_t>(kmax + 1), C{0.0, 0.0}));
  std::vector<C> b(static_cast<std::size_t>(kmax + 1), C{0.5, 0.0});
  M[0][0] = C{0.5, 0.0};
  for (int k = 1; k <= kmax; ++k) {
    M[0][static_cast<std::size_t>(k)] = C{1.0, 0.0};
  }
  for (int q = 1; q <= kmax; ++q) {
    for (int k = 1; k <= kmax; ++k) {
      C value{1.0, 0.0};
      for (int r = 1; r <= q; ++r) {
        value *= C{0.0, static_cast<double>(k * k)} / static_cast<double>(2 * r - 1);
      }
      M[static_cast<std::size_t>(q)][static_cast<std::size_t>(k)] = value;
    }
  }
  return solve_linear_complex(M, b);
}

// Stores the PN flux/phasing polynomials used by the averaged radiation
// reaction equations. The member names mirror the Python implementation and the
// paper notation: "a" and "b" are the y/e^2 evolution pieces, "k" is the
// pericenter-precession piece, and suffixes NS/SO/SS/SOSS label non-spinning,
// spin-orbit, spin-spin, and mixed spin terms. Keeping them as Poly objects
// avoids re-parsing or rebuilding coefficient vectors during every RHS call.
class PNDerivatives {
public:
  PNDerivatives(double m1, double m2, double chi_eff, double s2_1, double s2_2,
                double q1, double q2, int pn_phase_order, int pn_spin_order);

  std::array<double, 4> Dy_De2_Dl_Ddl(double y, double e2, double dchi, double dchi2,
                                      double sperp2) const;

private:
  int pn_phase_order_ = 6;
  int pn_spin_order_ = 6;
  int pn_max_order_ = 6;
  double nu_ = 0.25;

  Poly p_a0NS;
  Poly p_a2NS;
  Poly p_a3NS;
  Poly p_a4NS;
  Poly sqrt_a4NS;
  Poly p_a5NS;
  Poly p_a6NS;
  Poly sqrt_a6NS;
  Poly log_a6NS;
  Poly p_b0NS;
  Poly p_b2NS;
  Poly p_b3NS;
  Poly p_b4NS;
  Poly sqrt_b4NS;
  Poly p_b5NS;
  Poly p_b6NS;
  Poly sqrt_b6NS;
  double one_m_sqrt_b6NS = 0.0;
  Poly log_b6NS;
  Poly chi_p_a3SO;
  Poly dch_p_a3SO;
  Poly chi_p_a5SO;
  Poly chi_e2sqrt_a5SO;
  Poly dch_p_a5SO;
  Poly dch_e2sqrt_a5SO;
  Poly chi_p_b3SO;
  Poly dch_p_b3SO;
  Poly chi_p_b5SO;
  Poly chi_sqrt_b5SO;
  Poly dch_p_b5SO;
  Poly dch_sqrt_b5SO;
  Poly const_a4SS;
  Poly sperp2_a4SS;
  Poly chidch_a4SS;
  Poly dch2_a4SS;
  Poly const_p_a6SOSS;
  Poly chi2_sqrt_a6SS;
  Poly dch_p_a6SOSS;
  Poly chidch_sqrt_a6SS;
  Poly dch2_p_a6SS;
  Poly dch2_sqrt_a6SS;
  Poly const_b4SS;
  Poly sperp2_b4SS;
  Poly chidch_b4SS;
  Poly dch2_b4SS;
  Poly const_p_b6SOSS;
  Poly chi2_sqrt_b6SS;
  Poly dch_p_b6SOSS;
  Poly chidch_sqrt_b6SS;
  Poly dch2_p_b6SS;
  Poly dch2_sqrt_b6SS;
  double k0NS = 3.0;
  Poly p_k2NS;
  Poly p_k4NS;
  Poly sqrt_k4NS;
  double chi_k1SO = 0.0;
  double dch_k1SO = 0.0;
  Poly chi_p_k3SO;
  Poly dch_p_k3SO;
  double const_k2SS = 0.0;
  double sperp2_k2SS = 0.0;
  double chidch_k2SS = 0.0;
  double dch2_k2SS = 0.0;
  Poly chi2_p_k4SS;
  Poly chidch_p_k4SS;
  Poly dch2_p_k4SS;
};

PNDerivatives::PNDerivatives(double m1, double m2, double chi_eff, double s2_1, double s2_2,
                             double q1, double q2, int pn_phase_order, int pn_spin_order) {
  // A requested order of -1 follows the Python convention: use all implemented
  // terms. The constructor then precomputes all mass-ratio and quadrupole
  // dependent coefficients once per Model.
  pn_phase_order_ = (pn_phase_order == -1) ? 6 : pn_phase_order;
  pn_spin_order_ = (pn_spin_order == -1) ? 6 : pn_spin_order;
  pn_max_order_ = std::max(pn_phase_order_, pn_spin_order_);

  const auto mp = mass_params_from_m1_m2(m1, m2);
  nu_ = mp.nu;
  const double nu = mp.nu;
  const double dmu = mp.dmu;
  const double nu2 = nu * nu;
  const double nu3 = nu2 * nu;
  const double pi2 = kPi * kPi;

  const double dqS = q1 + q2 - 2.0;
  const double dqA = q1 - q2;
  const double dqAdmu = dqA * dmu;
  const double s2iS = s2_1 + s2_2;
  const double s2iA = s2_1 - s2_2;
  const double chi_eff2 = chi_eff * chi_eff;

  const std::vector<double> c_phiy =
      vec({1.0, 97.0 / 32.0, 49.0 / 128.0, -49.0 / 18432.0,
           -109.0 / 147456.0, -2567.0 / 58982400.0});
  const std::vector<double> c_phie =
      vec({1.0, 5969.0 / 3940.0, 24217.0 / 189120.0, 623.0 / 4538880.0,
           -96811.0 / 363110400.0, -5971.0 / 4357324800.0});
  const std::vector<double> c_psiy =
      vec({1.0, -207671.0 / 8318.0, -8382869.0 / 266176.0,
           -8437609.0 / 4791168.0, 10075915.0 / 306634752.0,
           -38077159.0 / 15331737600.0});
  const std::vector<double> c_zetay =
      vec({1.0, 113002.0 / 11907.0, 6035543.0 / 762048.0,
           253177.0 / 571536.0, -850489.0 / 877879296.0,
           -1888651.0 / 10973491200.0});
  const std::vector<double> c_psie =
      vec({1.0, -9904271.0 / 891056.0, -101704075.0 / 10692672.0,
           -217413779.0 / 513248256.0, 35703577.0 / 6843310080.0,
           -3311197679.0 / 9854366515200.0});
  const std::vector<double> c_zetae =
      vec({1.0, 11228233.0 / 2440576.0, 37095275.0 / 14643456.0,
           151238443.0 / 1405771776.0, -118111.0 / 611205120.0,
           -407523451.0 / 26990818099300.0});

  const auto c_kappay = add(
      add(scale(244.0 * std::log(2.0),
                vec({0.0, 1.0, -18881.0 / 1098.0, 6159821.0 / 39528.0,
                     -16811095.0 / 19764.0, 446132351.0 / 123525.0})),
          scale(-243.0 * std::log(3.0),
                vec({0.0, 1.0, -39.0 / 4.0, 2735.0 / 64.0, 25959.0 / 512.0,
                     -638032239.0 / 409600.0}))),
      add(scale(-48828125.0 * std::log(5.0) / 5184.0,
                vec({0.0, 0.0, 0.0, 1.0, -83.0 / 8.0, 12637.0 / 256.0})),
          scale(-4747561509943.0 * std::log(7.0) / 33177600.0,
                vec({0.0, 0.0, 0.0, 0.0, 0.0, 1.0}))));

  const auto c_kappae = add(
      add(scale(6536.0 * std::log(2.0),
                vec({1.0, -22314.0 / 817.0, 7170067.0 / 19608.0,
                     -10943033.0 / 4128.0, 230370959.0 / 15480.0,
                     -866124466133.0 / 8823600.0})),
          scale(-6561.0 * std::log(3.0),
                vec({1.0, -49.0 / 4.0, 4369.0 / 64.0, 214449.0 / 512.0,
                     -623830739.0 / 81920.0, 76513915569.0 / 1638400.0}))),
      add(scale(-48828125.0 * std::log(5.0) / 64.0,
                vec({0.0, 0.0, 1.0, -293.0 / 24.0, 159007.0 / 2304.0,
                     -6631171.0 / 27648.0})),
          scale(-4747561509943.0 * std::log(7.0) / 245760.0,
                vec({0.0, 0.0, 0.0, 0.0, 1.0, -259.0 / 20.0}))));

  const std::vector<double> c_thyc =
      vec({1.0, 21263.0 / 3008.0, 52387.0 / 12032.0,
           253973.0 / 1732608.0, -82103.0 / 13860864.0});
  const std::vector<double> c_thyd =
      vec({1.0, 1897.0 / 592.0, -461.0 / 2368.0, -42581.0 / 340992.0,
           -3803.0 / 1363968.0});
  const std::vector<double> c_thec =
      vec({1.0, 377077.0 / 92444.0, 7978379.0 / 4437312.0,
           5258749.0 / 106495488.0});
  const std::vector<double> c_thed =
      vec({1.0, 37477.0 / 19748.0, 95561.0 / 947904.0,
           -631523.0 / 22749696.0});

  p_a0NS = Poly({32.0 / 5.0, 28.0 / 5.0});
  p_a2NS = Poly({-1486.0 / 105.0 - (88.0 / 5.0) * nu,
                 12296.0 / 105.0 - (5258.0 / 45.0) * nu,
                 3007.0 / 84.0 - (244.0 / 9.0) * nu});
  p_a3NS = Poly(scale((128.0 / 5.0) * kPi, c_phiy));
  p_a4NS = Poly({34103.0 / 2835.0 + (13661.0 / 315.0) * nu + (944.0 / 45.0) * nu2,
                 -489191.0 / 1890.0 - (209729.0 / 630.0) * nu + (147443.0 / 270.0) * nu2,
                 2098919.0 / 7560.0 - (2928257.0 / 2520.0) * nu + (34679.0 / 45.0) * nu2,
                 53881.0 / 2520.0 - (7357.0 / 90.0) * nu + (9392.0 / 135.0) * nu2});
  sqrt_a4NS = Poly({16.0 - (32.0 / 5.0) * nu,
                    266.0 - (532.0 / 5.0) * nu,
                    -859.0 / 2.0 + (859.0 / 5.0) * nu,
                    -65.0 + 26.0 * nu});
  p_a5NS = Poly(add(scale(-(4159.0 / 105.0) * kPi, c_psiy),
                    scale(-(756.0 / 5.0) * kPi * nu, c_zetay)));
  p_a6NS = Poly(add(
      vec({16447322263.0 / 21829500.0 - (54784.0 / 525.0) * kEulerGamma +
               (512.0 / 15.0) * pi2 + (-(56198689.0 / 34020.0) + (902.0 / 15.0) * pi2) * nu +
               (541.0 / 140.0) * nu2 - (1121.0 / 81.0) * nu3,
           33232226053.0 / 10914750.0 - (392048.0 / 525.0) * kEulerGamma +
               (3664.0 / 15.0) * pi2 + (-(588778.0 / 1701.0) + (2747.0 / 40.0) * pi2) * nu -
               (846121.0 / 1260.0) * nu2 - (392945.0 / 324.0) * nu3,
           -227539553251.0 / 58212000.0 - (93304.0 / 175.0) * kEulerGamma +
               (872.0 / 5.0) * pi2 + ((124929721.0 / 12960.0) - (41287.0 / 960.0) * pi2) * nu +
               (148514441.0 / 30240.0) * nu2 - (2198212.0 / 405.0) * nu3,
           -300856627.0 / 67375.0 - (4922.0 / 175.0) * kEulerGamma +
               (46.0 / 5.0) * pi2 + ((1588607.0 / 432.0) - (369.0 / 80.0) * pi2) * nu +
               (12594313.0 / 3780.0) * nu2 - (44338.0 / 15.0) * nu3,
           -243511057.0 / 887040.0 + (4179523.0 / 15120.0) * nu +
               (83701.0 / 3780.0) * nu2 - (1876.0 / 15.0) * nu3,
           0.0}),
      scale(1284.0 / 175.0, c_kappay)));
  sqrt_a6NS = Poly({-616471.0 / 1575.0 + ((9874.0 / 315.0) - (41.0 / 30.0) * pi2) * nu +
                        (632.0 / 15.0) * nu2,
                    2385427.0 / 1050.0 + (-(274234.0 / 45.0) + (4223.0 / 240.0) * pi2) * nu +
                        (70946.0 / 45.0) * nu2,
                    8364697.0 / 4200.0 + ((1900517.0 / 630.0) - (32267.0 / 960.0) * pi2) * nu -
                        (47443.0 / 90.0) * nu2,
                    -167385119.0 / 25200.0 + ((4272491.0 / 504.0) - (123.0 / 160.0) * pi2) * nu -
                        (43607.0 / 18.0) * nu2,
                    -65279.0 / 168.0 + (510361.0 / 1260.0) * nu - (5623.0 / 45.0) * nu2});
  log_a6NS = Poly({54784.0 / 525.0, 392048.0 / 525.0, 93304.0 / 175.0,
                   4922.0 / 175.0});

  p_b0NS = Poly({608.0 / 15.0, 242.0 / 15.0});
  p_b2NS = Poly({-1878.0 / 35.0 - (8168.0 / 45.0) * nu,
                 59834.0 / 105.0 - (7753.0 / 15.0) * nu,
                 13929.0 / 140.0 - (3328.0 / 45.0) * nu});
  p_b3NS = Poly(scale((788.0 / 3.0) * kPi, c_phie));
  p_b4NS = Poly({-949877.0 / 945.0 + (18763.0 / 21.0) * nu + (1504.0 / 5.0) * nu2,
                 -3082783.0 / 1260.0 - (988423.0 / 420.0) * nu + (64433.0 / 20.0) * nu2,
                 23289859.0 / 7560.0 - (13018711.0 / 2520.0) * nu + (127411.0 / 45.0) * nu2,
                 420727.0 / 1680.0 - (362071.0 / 1260.0) * nu + (1642.0 / 9.0) * nu2});
  sqrt_b4NS = Poly({2672.0 / 3.0 - (5344.0 / 15.0) * nu,
                    2321.0 - (4642.0 / 5.0) * nu,
                    565.0 / 3.0 - (226.0 / 3.0) * nu});
  p_b5NS = Poly(add(scale(-(55691.0 / 105.0) * kPi, c_psie),
                    scale(-(610144.0 / 315.0) * kPi * nu, c_zetae)));
  p_b6NS = Poly(add(
      vec({61669369961.0 / 4365900.0 - (2633056.0 / 1575.0) * kEulerGamma +
               (24608.0 / 45.0) * pi2 + ((50099023.0 / 56700.0) + (779.0 / 5.0) * pi2) * nu -
               (4088921.0 / 1260.0) * nu2 - (61001.0 / 243.0) * nu3,
           66319591307.0 / 21829500.0 - (9525568.0 / 1575.0) * kEulerGamma +
               (89024.0 / 45.0) * pi2 + ((28141879.0 / 450.0) - (139031.0 / 480.0) * pi2) * nu -
               (21283907.0 / 1512.0) * nu2 - (86910509.0 / 9720.0) * nu3,
           -1149383987023.0 / 58212000.0 - (4588588.0 / 1575.0) * kEulerGamma +
               (42884.0 / 45.0) * pi2 + ((11499615139.0 / 453600.0) - (271871.0 / 960.0) * pi2) * nu +
               (61093675.0 / 2016.0) * nu2 - (2223241.0 / 90.0) * nu3,
           40262284807.0 / 4312000.0 - (20437.0 / 175.0) * kEulerGamma +
               (191.0 / 5.0) * pi2 + (-(5028323.0 / 280.0) - (6519.0 / 320.0) * pi2) * nu +
               (24757667.0 / 1260.0) * nu2 - (11792069.0 / 1215.0) * nu3,
           302322169.0 / 887040.0 - (1921387.0 / 5040.0) * nu +
               (41179.0 / 108.0) * nu2 - (386792.0 / 1215.0) * nu3,
           0.0}),
      scale(428.0 / 1575.0, c_kappae)));
  sqrt_b6NS = Poly({-22713049.0 / 7875.0 + (-(11053982.0 / 945.0) + (8323.0 / 90.0) * pi2) * nu +
                        (108664.0 / 45.0) * nu2,
                    178791374.0 / 7875.0 + (-(38295557.0 / 630.0) + (94177.0 / 480.0) * pi2) * nu +
                        (681989.0 / 45.0) * nu2,
                    5321445613.0 / 189000.0 + (-(26478311.0 / 756.0) + (2501.0 / 1440.0) * pi2) * nu +
                        (450212.0 / 45.0) * nu2,
                    186961.0 / 168.0 - (289691.0 / 252.0) * nu + (3197.0 / 9.0) * nu2});
  one_m_sqrt_b6NS = 1460336.0 / 23625.0;
  log_b6NS = Poly({2633056.0 / 1575.0, 9525568.0 / 1575.0,
                   4588588.0 / 1575.0, 20437.0 / 175.0});

  chi_p_a3SO = Poly(scale(chi_eff, vec({-752.0 / 15.0, -138.0, -611.0 / 30.0})));
  dch_p_a3SO = Poly(scale(dmu, vec({-152.0 / 15.0, -154.0 / 15.0, 17.0 / 30.0})));
  chi_p_a5SO = Poly(scale(chi_eff, vec({-5861.0 / 45.0 + (4004.0 / 15.0) * nu,
                                        -968539.0 / 630.0 + (259643.0 / 135.0) * nu,
                                        -4856917.0 / 2520.0 + (943721.0 / 540.0) * nu,
                                        -64903.0 / 560.0 + (5081.0 / 45.0) * nu})));
  chi_e2sqrt_a5SO = Poly(scale(chi_eff, vec({-1416.0 / 5.0 + (1652.0 / 15.0) * nu,
                                             2469.0 / 5.0 - (5761.0 / 30.0) * nu,
                                             222.0 / 5.0 - (259.0 / 15.0) * nu})));
  dch_p_a5SO = Poly(scale(dmu, vec({-21611.0 / 315.0 + (632.0 / 15.0) * nu,
                                    -55415.0 / 126.0 + (36239.0 / 135.0) * nu,
                                    -72631.0 / 360.0 + (12151.0 / 108.0) * nu,
                                    909.0 / 560.0 - (143.0 / 45.0) * nu})));
  dch_e2sqrt_a5SO = Poly(scale(dmu, vec({-472.0 / 5.0 + (236.0 / 15.0) * nu,
                                         823.0 / 5.0 - (823.0 / 30.0) * nu,
                                         74.0 / 5.0 - (37.0 / 15.0) * nu})));

  const auto chi_p_a6SO_arr = scale(-(3008.0 / 15.0) * kPi * chi_eff, c_thyc);
  const auto dch_p_a6SO_arr = scale(-(592.0 / 15.0) * kPi * dmu, c_thyd);

  chi_p_b3SO = Poly(scale(chi_eff, vec({-3272.0 / 9.0, -26263.0 / 45.0, -812.0 / 15.0})));
  dch_p_b3SO = Poly(scale(dmu, vec({-3328.0 / 45.0, -1993.0 / 45.0, 23.0 / 15.0})));
  chi_p_b5SO = Poly(scale(chi_eff, vec({-13103.0 / 35.0 + (289208.0 / 135.0) * nu,
                                        -548929.0 / 63.0 + (61355.0 / 6.0) * nu,
                                        -6215453.0 / 840.0 + (1725437.0 / 270.0) * nu,
                                        -87873.0 / 280.0 + (13177.0 / 45.0) * nu})));
  chi_sqrt_b5SO = Poly(scale(chi_eff, vec({-1184.0 + (4144.0 / 9.0) * nu,
                                           -13854.0 / 5.0 + (16163.0 / 15.0) * nu,
                                           -626.0 / 5.0 + (2191.0 / 45.0) * nu})));
  dch_p_b5SO = Poly(scale(dmu, vec({-32857.0 / 105.0 + (52916.0 / 135.0) * nu,
                                    -1396159.0 / 630.0 + (126833.0 / 90.0) * nu,
                                    -203999.0 / 280.0 + (56368.0 / 135.0) * nu,
                                    5681.0 / 1120.0 - (376.0 / 45.0) * nu})));
  dch_sqrt_b5SO = Poly(scale(dmu, vec({-1184.0 / 3.0 + (592.0 / 9.0) * nu,
                                       -4618.0 / 5.0 + (2309.0 / 15.0) * nu,
                                       -626.0 / 15.0 + (313.0 / 45.0) * nu})));
  const auto chi_p_b6SO_arr = scale(-(92444.0 / 45.0) * kPi * chi_eff, c_thec);
  const auto dch_p_b6SO_arr = scale(-(19748.0 / 45.0) * kPi * dmu, c_thed);

  const auto c_s2iS_a4SS = scale(s2iS, vec({8.0 / 5.0 - 8.0 * dqS,
                                            24.0 / 5.0 - (108.0 / 5.0) * dqS,
                                            3.0 / 5.0 - (63.0 / 20.0) * dqS}));
  const auto c_s2iA_a4SS = scale(dqA * s2iA, vec({-8.0, -108.0 / 5.0, -63.0 / 20.0}));
  const auto c_chi2_a4SS = scale(chi_eff2, vec({156.0 / 5.0 + 12.0 * dqS,
                                                84.0 + (162.0 / 5.0) * dqS,
                                                123.0 / 10.0 + (189.0 / 40.0) * dqS}));
  const_a4SS = Poly(add(add(c_s2iS_a4SS, c_s2iA_a4SS), c_chi2_a4SS));
  sperp2_a4SS = Poly({-84.0 / 5.0, -228.0 / 5.0, -33.0 / 5.0});
  chidch_a4SS = Poly(scale(chi_eff * dqA, vec({24.0, 324.0 / 5.0, 189.0 / 20.0})));
  dch2_a4SS = Poly({-2.0 / 5.0 + 12.0 * dqS,
                    -6.0 / 5.0 + (162.0 / 5.0) * dqS,
                    -3.0 / 20.0 + (189.0 / 40.0) * dqS});

  const auto chi2_p_a6SS_arr = scale(chi_eff2, vec({
      30596.0 / 105.0 + (2539.0 / 105.0) * dqS + (443.0 / 30.0) * dqAdmu + (-(688.0 / 5.0) - (172.0 / 5.0) * dqS) * nu,
      115078.0 / 45.0 + (21317.0 / 60.0) * dqS + (3253.0 / 60.0) * dqAdmu + (-(3962.0 / 3.0) - (1981.0 / 6.0) * dqS) * nu,
      4476649.0 / 2520.0 + (133703.0 / 420.0) * dqS + (481.0 / 48.0) * dqAdmu + (-(53267.0 / 45.0) - (53267.0 / 180.0) * dqS) * nu,
      17019.0 / 140.0 + (29831.0 / 1120.0) * dqS + (29.0 / 160.0) * dqAdmu + (-(1343.0 / 15.0) - (1343.0 / 60.0) * dqS) * nu,
      0.0}));
  chi2_sqrt_a6SS = Poly(scale(chi_eff2, vec({
      -(244.0 / 15.0) - (52.0 / 15.0) * dqS - (4.0 / 15.0) * dqAdmu + (16.0 / 5.0 + (4.0 / 5.0) * dqS) * nu,
      6283.0 / 30.0 + (1339.0 / 30.0) * dqS + (103.0 / 30.0) * dqAdmu + (-(206.0 / 5.0) - (103.0 / 10.0) * dqS) * nu,
      -(48007.0 / 120.0) - (10231.0 / 120.0) * dqS - (787.0 / 120.0) * dqAdmu + (787.0 / 10.0 + (787.0 / 40.0) * dqS) * nu,
      -(183.0 / 20.0) - (39.0 / 20.0) * dqS - (3.0 / 20.0) * dqAdmu + (9.0 / 5.0 + (9.0 / 20.0) * dqS) * nu})));
  const auto chidch_p_a6SS_arr = scale(chi_eff, vec({
      (3134.0 / 15.0 + (443.0 / 15.0) * dqS) * dmu + (5078.0 / 105.0 - (344.0 / 5.0) * nu) * dqA,
      (30421.0 / 45.0 + (3253.0 / 30.0) * dqS) * dmu + (21317.0 / 30.0 - (1981.0 / 3.0) * nu) * dqA,
      (-(111.0 / 5.0) + (481.0 / 24.0) * dqS) * dmu + (133703.0 / 210.0 - (53267.0 / 90.0) * nu) * dqA,
      (-(149.0 / 40.0) + (29.0 / 80.0) * dqS) * dmu + (29831.0 / 560.0 - (1343.0 / 30.0) * nu) * dqA,
      0.0}));
  chidch_sqrt_a6SS = Poly(scale(chi_eff, vec({
      (-(104.0 / 15.0) - (8.0 / 15.0) * dqS) * dmu + (-(104.0 / 15.0) + (8.0 / 5.0) * nu) * dqA,
      (1339.0 / 15.0 + (103.0 / 15.0) * dqS) * dmu + (1339.0 / 15.0 - (103.0 / 5.0) * nu) * dqA,
      (-(10231.0 / 60.0) - (787.0 / 60.0) * dqS) * dmu + (-(10231.0 / 60.0) + (787.0 / 20.0) * nu) * dqA,
      (-(39.0 / 10.0) - (3.0 / 10.0) * dqS) * dmu + (-(39.0 / 10.0) + (9.0 / 10.0) * nu) * dqA})));
  dch2_p_a6SS = Poly({39.0 / 5.0 + (2539.0 / 105.0) * dqS + (443.0 / 30.0) * dqAdmu + (-(1163.0 / 15.0) - (172.0 / 5.0) * dqS) * nu,
                      659.0 / 15.0 + (21317.0 / 60.0) * dqS + (3253.0 / 60.0) * dqAdmu + (-(2399.0 / 15.0) - (1981.0 / 6.0) * dqS) * nu,
                      1769.0 / 90.0 + (133703.0 / 420.0) * dqS + (481.0 / 48.0) * dqAdmu + (2021.0 / 72.0 - (53267.0 / 180.0) * dqS) * nu,
                      19.0 / 10.0 + (29831.0 / 1120.0) * dqS + (29.0 / 160.0) * dqAdmu + (-(3.0 / 10.0) - (1343.0 / 60.0) * dqS) * nu});
  dch2_sqrt_a6SS = Poly({-(4.0 / 15.0) - (52.0 / 15.0) * dqS - (4.0 / 15.0) * dqAdmu + (32.0 / 15.0 + (4.0 / 5.0) * dqS) * nu,
                         103.0 / 30.0 + (1339.0 / 30.0) * dqS + (103.0 / 30.0) * dqAdmu + (-(412.0 / 15.0) - (103.0 / 10.0) * dqS) * nu,
                         -(787.0 / 120.0) - (10231.0 / 120.0) * dqS - (787.0 / 120.0) * dqAdmu + (787.0 / 15.0 + (787.0 / 40.0) * dqS) * nu,
                         -(3.0 / 20.0) - (39.0 / 20.0) * dqS - (3.0 / 20.0) * dqAdmu + (6.0 / 5.0 + (9.0 / 20.0) * dqS) * nu});

  const auto c_s2iS_b4SS = scale(s2iS, vec({-4.0 / 3.0,
                                            34.0 / 3.0 - (938.0 / 15.0) * dqS,
                                            49.0 / 2.0 - (595.0 / 6.0) * dqS,
                                            9.0 / 4.0 - (37.0 / 4.0) * dqS}));
  const auto c_s2iA_b4SS = scale(dqA * s2iA, vec({0.0, -938.0 / 15.0, -595.0 / 6.0, -37.0 / 4.0}));
  const auto c_chi2_b4SS = scale(chi_eff2, vec({2.0 / 3.0,
                                                3667.0 / 15.0 + (469.0 / 5.0) * dqS,
                                                4613.0 / 12.0 + (595.0 / 4.0) * dqS,
                                                287.0 / 8.0 + (111.0 / 8.0) * dqS}));
  const_b4SS = Poly(add(add(c_s2iS_b4SS, c_s2iA_b4SS), c_chi2_b4SS));
  sperp2_b4SS = Poly({2.0 / 3.0, -1961.0 / 15.0, -2527.0 / 12.0, -157.0 / 8.0});
  chidch_b4SS = Poly(scale(chi_eff * dqA, vec({0.0, 938.0 / 5.0, 595.0 / 2.0, 111.0 / 4.0})));
  dch2_b4SS = Poly({2.0 / 3.0, 1.0 / 3.0 + (469.0 / 5.0) * dqS,
                    -13.0 / 4.0 + (595.0 / 4.0) * dqS,
                    -3.0 / 8.0 + (111.0 / 8.0) * dqS});

  const auto chi2_p_b6SS_arr = scale(chi_eff2, vec({
      1468414.0 / 945.0 + (2852.0 / 105.0) * dqS + (3461.0 / 30.0) * dqAdmu + (-(57844.0 / 45.0) - (14461.0 / 45.0) * dqS) * nu,
      47715853.0 / 3780.0 + (1464091.0 / 840.0) * dqS + (11007.0 / 40.0) * dqAdmu + (-(21865.0 / 3.0) - (21865.0 / 12.0) * dqS) * nu,
      4255831.0 / 504.0 + (166844.0 / 105.0) * dqS + (2941.0 / 48.0) * dqAdmu + (-(222533.0 / 45.0) - (222533.0 / 180.0) * dqS) * nu,
      414027.0 / 1120.0 + (365363.0 / 4480.0) * dqS + (511.0 / 640.0) * dqAdmu + (-(1287.0 / 5.0) - (1287.0 / 20.0) * dqS) * nu}));
  chi2_sqrt_b6SS = Poly(scale(chi_eff2, vec({
      49532.0 / 45.0 + (10556.0 / 45.0) * dqS + (812.0 / 45.0) * dqAdmu + (-(3248.0 / 15.0) - (812.0 / 15.0) * dqS) * nu,
      140117.0 / 60.0 + (29861.0 / 60.0) * dqS + (2297.0 / 60.0) * dqAdmu + (-(2297.0 / 5.0) - (2297.0 / 20.0) * dqS) * nu,
      3721.0 / 180.0 + (793.0 / 180.0) * dqS + (61.0 / 180.0) * dqAdmu + (-(61.0 / 15.0) - (61.0 / 60.0) * dqS) * nu})));
  const auto chidch_p_b6SS_arr = scale(chi_eff, vec({
      (176426.0 / 135.0 + (3461.0 / 15.0) * dqS) * dmu + (5704.0 / 105.0 - (28922.0 / 45.0) * nu) * dqA,
      (387212.0 / 135.0 + (11007.0 / 20.0) * dqS) * dmu + (1464091.0 / 420.0 - (21865.0 / 6.0) * nu) * dqA,
      (2562.0 / 5.0 + (2941.0 / 24.0) * dqS) * dmu + (333688.0 / 105.0 - (222533.0 / 90.0) * nu) * dqA,
      (-(33.0 / 32.0) + (511.0 / 320.0) * dqS) * dmu + (365363.0 / 2240.0 - (1287.0 / 10.0) * nu) * dqA}));
  chidch_sqrt_b6SS = Poly(scale(chi_eff, vec({
      (21112.0 / 45.0 + (1624.0 / 45.0) * dqS) * dmu + (21112.0 / 45.0 - (1624.0 / 15.0) * nu) * dqA,
      (29861.0 / 30.0 + (2297.0 / 30.0) * dqS) * dmu + (29861.0 / 30.0 - (2297.0 / 10.0) * nu) * dqA,
      (793.0 / 90.0 + (61.0 / 90.0) * dqS) * dmu + (793.0 / 90.0 - (61.0 / 30.0) * nu) * dqA})));
  dch2_p_b6SS = Poly({8887.0 / 135.0 + (2852.0 / 105.0) * dqS + (3461.0 / 30.0) * dqAdmu + (-(13127.0 / 27.0) - (14461.0 / 45.0) * dqS) * nu,
                      161077.0 / 540.0 + (1464091.0 / 840.0) * dqS + (11007.0 / 40.0) * dqAdmu + (-(185723.0 / 270.0) - (21865.0 / 12.0) * dqS) * nu,
                      14827.0 / 90.0 + (166844.0 / 105.0) * dqS + (2941.0 / 48.0) * dqAdmu + (-(45373.0 / 360.0) - (222533.0 / 180.0) * dqS) * nu,
                      283.0 / 32.0 + (365363.0 / 4480.0) * dqS + (511.0 / 640.0) * dqAdmu + (-(117.0 / 20.0) - (1287.0 / 20.0) * dqS) * nu});
  dch2_sqrt_b6SS = Poly({812.0 / 45.0 + (10556.0 / 45.0) * dqS + (812.0 / 45.0) * dqAdmu + (-(6496.0 / 45.0) - (812.0 / 15.0) * dqS) * nu,
                         2297.0 / 60.0 + (29861.0 / 60.0) * dqS + (2297.0 / 60.0) * dqAdmu + (-(4594.0 / 15.0) - (2297.0 / 20.0) * dqS) * nu,
                         61.0 / 180.0 + (793.0 / 180.0) * dqS + (61.0 / 180.0) * dqAdmu + (-(122.0 / 45.0) - (61.0 / 60.0) * dqS) * nu});

  const_p_a6SOSS = Poly(add(chi_p_a6SO_arr, chi2_p_a6SS_arr));
  const_p_b6SOSS = Poly(add(chi_p_b6SO_arr, chi2_p_b6SS_arr));
  dch_p_a6SOSS = Poly(add(dch_p_a6SO_arr, chidch_p_a6SS_arr));
  dch_p_b6SOSS = Poly(add(dch_p_b6SO_arr, chidch_p_b6SS_arr));

  k0NS = 3.0;
  p_k2NS = Poly({27.0 / 2.0 - 7.0 * nu, 51.0 / 4.0 - (13.0 / 2.0) * nu});
  p_k4NS = Poly({105.0 / 2.0 + (-(625.0 / 4.0) + (123.0 / 32.0) * pi2) * nu + 7.0 * nu2,
                 573.0 / 4.0 + (-(357.0 / 2.0) + (123.0 / 128.0) * pi2) * nu + 40.0 * nu2,
                 39.0 / 2.0 - (55.0 / 4.0) * nu + (65.0 / 8.0) * nu2});
  sqrt_k4NS = Poly({15.0 - 6.0 * nu, 30.0 - 12.0 * nu});
  chi_k1SO = -(7.0 / 2.0) * chi_eff;
  dch_k1SO = -0.5 * dmu;
  chi_p_k3SO = Poly(scale(chi_eff, vec({-26.0 + 8.0 * nu,
                                        -(105.0 / 4.0) + (49.0 / 4.0) * nu})));
  dch_p_k3SO = Poly(scale(dmu, vec({-8.0 + 0.5 * nu,
                                    -(15.0 / 4.0) + (7.0 / 4.0) * nu})));
  const_k2SS = -(3.0 / 8.0) * dqS * s2iS -
               (3.0 / 8.0) * dqA * s2iA +
               (3.0 / 2.0 + (9.0 / 16.0) * dqS) * chi_eff2;
  sperp2_k2SS = -3.0 / 4.0;
  chidch_k2SS = (9.0 / 8.0) * dqA * chi_eff;
  dch2_k2SS = (9.0 / 16.0) * dqS;
  chi2_p_k4SS = Poly(scale(chi_eff2, vec({
      181.0 / 8.0 + (33.0 / 8.0) * dqS + (3.0 / 4.0) * dqAdmu + (-(5.0 / 2.0) - (5.0 / 8.0) * dqS) * nu,
      369.0 / 16.0 + (75.0 / 16.0) * dqS + (3.0 / 16.0) * dqAdmu + (-(29.0 / 4.0) - (29.0 / 16.0) * dqS) * nu})));
  chidch_p_k4SS = Poly(scale(chi_eff, vec({
      (43.0 / 4.0 + (3.0 / 2.0) * dqS) * dmu + (33.0 / 4.0 - (5.0 / 4.0) * nu) * dqA,
      (21.0 / 8.0 + (3.0 / 8.0) * dqS) * dmu + (75.0 / 8.0 - (29.0 / 8.0) * nu) * dqA})));
  dch2_p_k4SS = Poly({1.0 / 8.0 + (33.0 / 8.0) * dqS + (3.0 / 4.0) * dqAdmu + (-(7.0 / 2.0) - (5.0 / 8.0) * dqS) * nu,
                      -(3.0 / 16.0) + (75.0 / 16.0) * dqS + (3.0 / 16.0) * dqAdmu - (29.0 / 16.0) * dqS * nu});
}

std::array<double, 4> PNDerivatives::Dy_De2_Dl_Ddl(double y, double e2, double dchi,
                                                   double dchi2, double sperp2) const {
  // Return the four averaged derivatives that drive the inspiral:
  //   dy/dt, d(e^2)/dt, d lambda/dt, d delta_lambda/dt.
  // The long body below is mostly polynomial bookkeeping from the PN expansion.
  // The branch structure mirrors PN order: higher-order pieces are accumulated
  // first and then multiplied by powers of y.
  e2 = clamp(e2, 0.0, 1.0 - 1.0e-12);
  const double sqrt1me2 = safe_sqrt(1.0 - e2);
  const double one_m_sqrt = 1.0 - sqrt1me2;
  const double sqrt_a = (1.0 - sqrt1me2) / sqrt1me2;
  const double e2sqrt = e2 / sqrt1me2;
  const double log_fact = std::log((1.0 + sqrt1me2) / (8.0 * y * sqrt1me2 * (1.0 - e2)));
  const double y2 = y * y;
  const double y3 = y2 * y;
  const double y8 = std::pow(y, 8);

  double Dy = 0.0;
  double De2 = 0.0;
  double k = 0.0;

  if (pn_max_order_ >= 6) {
    if (pn_phase_order_ >= 6) {
      Dy += p_a6NS(e2) + sqrt_a * sqrt_a6NS(e2) + log_fact * log_a6NS(e2);
      De2 += e2 * (p_b6NS(e2) + sqrt1me2 * sqrt_b6NS(e2) + log_fact * log_b6NS(e2)) +
             one_m_sqrt * one_m_sqrt_b6NS;
      k += p_k4NS(e2) + sqrt1me2 * sqrt_k4NS(e2);
    }
    if (pn_spin_order_ >= 6) {
      Dy += const_p_a6SOSS(e2) + dchi * dch_p_a6SOSS(e2) +
            sqrt_a * chi2_sqrt_a6SS(e2) +
            dchi * sqrt_a * chidch_sqrt_a6SS(e2) +
            dchi2 * (dch2_p_a6SS(e2) + sqrt_a * dch2_sqrt_a6SS(e2));
      De2 += e2 * (const_p_b6SOSS(e2) + dchi * dch_p_b6SOSS(e2) +
                   sqrt1me2 * chi2_sqrt_b6SS(e2) +
                   dchi * sqrt1me2 * chidch_sqrt_b6SS(e2) +
                   dchi2 * (dch2_p_b6SS(e2) + sqrt1me2 * dch2_sqrt_b6SS(e2)));
      k += chi2_p_k4SS(e2) + dchi * chidch_p_k4SS(e2) + dchi2 * dch2_p_k4SS(e2);
    }
    Dy *= y;
    De2 *= y;
    k *= y;
  }
  if (pn_max_order_ >= 5) {
    if (pn_phase_order_ >= 5) {
      Dy += p_a5NS(e2);
      De2 += e2 * p_b5NS(e2);
    }
    if (pn_spin_order_ >= 5) {
      Dy += chi_p_a5SO(e2) + e2sqrt * chi_e2sqrt_a5SO(e2) +
            dchi * (dch_p_a5SO(e2) + e2sqrt * dch_e2sqrt_a5SO(e2));
      De2 += e2 * (chi_p_b5SO(e2) + sqrt1me2 * chi_sqrt_b5SO(e2) +
                    dchi * (dch_p_b5SO(e2) + sqrt1me2 * dch_sqrt_b5SO(e2)));
      k += chi_p_k3SO(e2) + dchi * dch_p_k3SO(e2);
    }
    Dy *= y;
    De2 *= y;
    k *= y;
  }
  if (pn_max_order_ >= 4) {
    if (pn_phase_order_ >= 4) {
      Dy += p_a4NS(e2) + sqrt_a * sqrt_a4NS(e2);
      De2 += e2 * (p_b4NS(e2) + sqrt1me2 * sqrt_b4NS(e2));
      k += p_k2NS(e2);
    }
    if (pn_spin_order_ >= 4) {
      Dy += const_a4SS(e2) + sperp2 * sperp2_a4SS(e2) +
            dchi * chidch_a4SS(e2) + dchi2 * dch2_a4SS(e2);
      De2 += const_b4SS(e2) + sperp2 * sperp2_b4SS(e2) +
             dchi * chidch_b4SS(e2) + dchi2 * dch2_b4SS(e2);
      k += const_k2SS + sperp2 * sperp2_k2SS +
           dchi * chidch_k2SS + dchi2 * dch2_k2SS;
    }
    Dy *= y;
    De2 *= y;
    k *= y;
  }
  if (pn_max_order_ >= 3) {
    if (pn_phase_order_ >= 3) {
      Dy += p_a3NS(e2);
      De2 += e2 * p_b3NS(e2);
    }
    if (pn_spin_order_ >= 3) {
      Dy += chi_p_a3SO(e2) + dchi * dch_p_a3SO(e2);
      De2 += e2 * (chi_p_b3SO(e2) + dchi * dch_p_b3SO(e2));
      k += chi_k1SO + dchi * dch_k1SO;
    }
    Dy *= y;
    De2 *= y;
    k *= y;
  }
  if (pn_phase_order_ >= 2) {
    Dy = y2 * (Dy + p_a2NS(e2));
    De2 = y2 * (De2 + e2 * p_b2NS(e2));
    k = y2 * (k + k0NS);
  }

  Dy = y * y8 * nu_ * (p_a0NS(e2) + Dy);
  De2 = -y8 * nu_ * (e2 * p_b0NS(e2) + De2);
  const double Dl = y3;
  const double Ddl = k * Dl / (1.0 + k);
  return {Dy, De2, Dl, Ddl};
}

struct BasicPrecQuantities {
  double m = 0.0;
  double dchi_av = 0.0;
  double dchi_diff = 0.0;
  double chi_eff = 0.0;
  double J = 0.0;
  double L = 0.0;
  double sqY3mYm = 0.0;
  double Pp = 0.0;
  double Pm = 0.0;
  double dmudchiav_m_dchi0 = 0.0;
  double dmudchi_diff = 0.0;
  double PI_fact_p = 0.0;
  double PI_fact_m = 0.0;
  double PI_arg_p = 0.0;
  double PI_arg_m = 0.0;
};

std::pair<double, double> precession_average_factors_betasigma_small_m(double m) {
  const double m_factor_dchi_prec_avg =
      -m * (1.0 + m * (-1.0 + m * (71.0 / 384.0))) /
      (16.0 + m * (-24.0 + m * ((59.0 / 6.0) - m * (11.0 / 12.0))));
  const double m2 = m * m;
  const double m_factor_sigma =
      m2 / (1024.0 + m * (-1024.0 + m * (96.0 - m2 * (133.0 / 32.0) * (1.0 + m))));
  return {m_factor_dchi_prec_avg, m_factor_sigma};
}

std::pair<double, double> precession_average_factors_betasigma_large_m(double m) {
  const double E_m = ellipe(m);
  const double K_m = ellipk(m);
  const double E_K_m = E_m / K_m;
  const double m_factor_dchi_prec_avg = (E_K_m - 1.0 + 0.5 * m) / m;
  const double m_factor_sigma =
      ((1.0 / 3.0) + m * ((-1.0 / 3.0) + m / 8.0) +
       E_K_m * ((2.0 / 3.0) * (m - 2.0) + E_K_m)) /
      (m * m);
  return {m_factor_dchi_prec_avg, m_factor_sigma};
}

std::pair<double, double> precession_average_factors_betasigma(double m) {
  if (m < 0.3) {
    return precession_average_factors_betasigma_small_m(m);
  }
  return precession_average_factors_betasigma_large_m(m);
}

BasicPrecQuantities basic_prec_quantities(double y, double DJ2, double m1, double m2,
                                          double sz_1, double sz_2, double sp2_1,
                                          double sp2_2) {
  // Algebraic precession quantities at a single y and DeltaJ^2. This is shared
  // by the averaged derivatives and by Euler-angle reconstruction.
  const auto mp = mass_params_from_m1_m2(m1, m2);
  const double mu1 = mp.mu1;
  const double mu2 = mp.mu2;
  const double nu = mp.nu;
  const double dmu = mp.dmu;

  const double chi_eff = sz_1 + sz_2;
  const double dchi0 = sz_1 - sz_2;
  const double dmu2 = dmu * dmu;
  const double j2a = dmu2 / y + (chi_eff * chi_eff) * y;
  const double dChi = dmu * dchi0;
  const double sp2_tot = sp2_1 + sp2_2 + DJ2;
  const double bp = y * sp2_tot;
  const double dp = y * dmu2 * (4.0 * sp2_1 * sp2_2 - DJ2 * DJ2);

  const double pal = (j2a - 2.0 * dChi + bp) / 3.0;
  const double pal2 = pal * pal;
  const double pperp = bp * dChi - dmu * (DJ2 * dmu + (sp2_1 - sp2_2) * y * chi_eff);
  const double pal_pal2pperp = pal * (pal2 + pperp);
  const double discy6 = (pperp * pperp) * (9.0 * pal2 + 8.0 * pperp) / 27.0 +
                        dp * (pal_pal2pperp - 0.25 * dp);
  const double argG_3 =
      discy6 >= 0.0 ? std::atan2(std::sqrt(discy6), pal_pal2pperp - 0.5 * dp) / 3.0 : 0.0;
  const double sargG_3 = std::sin(argG_3);
  const double cargG_3 = std::cos(argG_3);
  const double cargG_3_m_pi_6 = 0.5 * (std::sqrt(3.0) * cargG_3 + sargG_3);
  const double py2 = 3.0 * pal2 + 2.0 * pperp;
  const double sqpy2 = std::sqrt(std::abs(py2));

  double dmudchiav_m_dchi0 = 0.0;
  if ((py2 > 0.0) && (pal > 0.0)) {
    dmudchiav_m_dchi0 =
        (pal2 * sargG_3 * sargG_3 - (2.0 / 3.0) * pperp * cargG_3 * cargG_3) /
        (pal + std::pow(3.0, -0.5) * sqpy2 * cargG_3);
  } else {
    dmudchiav_m_dchi0 = pal - std::pow(3.0, -0.5) * sqpy2 * cargG_3;
  }
  const double dmudchi_diff = sqpy2 * sargG_3;

  double dchi_av = 0.0;
  double dchi_diff = 0.0;
  if (dmu > 1.0e-12) {
    dchi_av = dchi0 + dmudchiav_m_dchi0 / dmu;
    dchi_diff = dmudchi_diff / dmu;
  } else {
    const double s2_tot = sp2_tot + chi_eff * chi_eff;
    if (s2_tot != 0.0) {
      dchi_av = chi_eff * (sp2_1 - sp2_2 + dchi0 * chi_eff) / s2_tot;
      dchi_diff =
          std::sqrt(std::max(0.0, sp2_tot * (4.0 * (sp2_2 * sz_1 * sz_1 +
                                                   sp2_1 * (sz_2 * sz_2 + sp2_2)) -
                                             DJ2 * (4.0 * sz_1 * sz_2 + DJ2)))) /
          s2_tot;
    }
  }

  double m = sargG_3 / cargG_3_m_pi_6;
  m = clamp(m, 0.0, 1.0 - 1.0e-14);
  const double sqY3mYm = std::sqrt(std::max(0.0, 2.0 * sqpy2 * cargG_3_m_pi_6 / y));
  const double L = nu / y;
  const double Sperp2_1 = (mu1 * mu1) * sp2_1;
  const double Sperp2_2 = (mu2 * mu2) * sp2_2;
  const double J0Lh = L + 0.5 * (chi_eff + dChi);
  const double J0Lh2 = J0Lh * J0Lh;
  const double Jperp2 = Sperp2_1 + Sperp2_2 + nu * DJ2;
  const double J = std::sqrt(std::max(0.0, J0Lh2 + Jperp2));
  const double xJ = 0.5 * Jperp2 / J0Lh2;
  const bool small_x = std::abs(xJ) < 1.0e-6;
  const double dJ_small_x = J0Lh * xJ * (1.0 - 0.5 * xJ * (1.0 - xJ));
  const double J_p_J0Lh = (small_x && J0Lh < 0.0) ? -dJ_small_x : J + J0Lh;
  const double J_m_J0Lh = (small_x && J0Lh > 0.0) ? dJ_small_x : J - J0Lh;

  const double muSz = 2.0 * (mu1 * mu1 * sz_1 + mu2 * mu2 * sz_2);
  const double dmudSp2 = dmu * (Sperp2_1 - Sperp2_2);
  const double Np = J_p_J0Lh * (J_p_J0Lh - muSz) - dmudSp2;
  const double Nm = J_m_J0Lh * (J_m_J0Lh + muSz) - dmudSp2;
  const double min_2J1pmcthL = 1.0e-10;
  double Bp_m_Cp = 2.0 * J_p_J0Lh + dmudchiav_m_dchi0 - dmudchi_diff;
  Bp_m_Cp = std::max(Bp_m_Cp, min_2J1pmcthL);
  double Bm_p_Cm = 2.0 * J_m_J0Lh - dmudchiav_m_dchi0 - dmudchi_diff;
  Bm_p_Cm = std::max(Bm_p_Cm, min_2J1pmcthL);
  const double Bm_m_Cm = Bm_p_Cm + 2.0 * dmudchi_diff;
  const double PI_fact_p = Np / Bp_m_Cp;
  const double PI_fact_m = Nm / Bm_m_Cm;
  const double PI_arg_p = -2.0 * dmudchi_diff / Bp_m_Cp;
  const double PI_arg_m = 2.0 * dmudchi_diff / Bm_m_Cm;
  const double Pp = PI_fact_p * ellippi(PI_arg_p, m);
  const double Pm = PI_fact_m * ellippi(PI_arg_m, m);

  return {m, dchi_av, dchi_diff, chi_eff, J, L, sqY3mYm, Pp, Pm,
          dmudchiav_m_dchi0, dmudchi_diff, PI_fact_p, PI_fact_m, PI_arg_p, PI_arg_m};
}

struct PrecAvgForDv {
  double J = 0.0;
  double L = 0.0;
  double chi_eff = 0.0;
  double K_m = 0.0;
  double sqY3mYm = 0.0;
  double Pp = 0.0;
  double Pm = 0.0;
  double dchi_prec_avg = 0.0;
  double dmudchi_prec_avg_m_dchi0 = 0.0;
  double dchi2_prec_avg = 0.0;
  double sperp2_prec_avg = 0.0;
};

PrecAvgForDv prec_avg_quantities_for_Dv(double y, double DJ2, double m1, double m2,
                                        double sz_1, double sz_2, double sp2_1,
                                        double sp2_2) {
  // Convert the instantaneous precession solution into the averaged quantities
  // that appear in the radiation-reaction ODE.
  const auto b = basic_prec_quantities(y, DJ2, m1, m2, sz_1, sz_2, sp2_1, sp2_2);
  const auto factors = precession_average_factors_betasigma(b.m);
  const double m_factor_dchi_prec_avg = factors.first;
  const double m_factor_sigma = factors.second;
  const double dchi_prec_avg = b.dchi_av - 2.0 * b.dchi_diff * m_factor_dchi_prec_avg;
  const double dmudchi_prec_avg_m_dchi0 =
      b.dmudchiav_m_dchi0 - 2.0 * b.dmudchi_diff * m_factor_dchi_prec_avg;
  const double dchi2_prec_avg =
      dchi_prec_avg * dchi_prec_avg +
      b.dchi_diff * b.dchi_diff * (0.5 - 4.0 * m_factor_sigma);
  const double sperp2_prec_avg = sp2_1 + sp2_2 + DJ2 - dmudchi_prec_avg_m_dchi0 / y;
  return {b.J, b.L, b.chi_eff, ellipk(b.m), b.sqY3mYm, b.Pp, b.Pm,
          dchi_prec_avg, dmudchi_prec_avg_m_dchi0, dchi2_prec_avg, sperp2_prec_avg};
}

struct EulerAngles {
  double dphiz = 0.0;
  double dzeta = 0.0;
  double costhL = 1.0;
};

EulerAngles precession_Euler_angles(double bpsip, double y, double DJ2, double m1, double m2,
                                    double sz_1, double sz_2, double sp2_1, double sp2_2) {
  // Recover Euler-angle offsets from the precession phase. These angles rotate
  // the co-precessing quadrupole modes into the inertial observer frame.
  const auto mp = mass_params_from_m1_m2(m1, m2);
  const auto b = basic_prec_quantities(y, DJ2, m1, m2, sz_1, sz_2, sp2_1, sp2_2);
  const double dzeta_E_fact = y * b.sqY3mYm / (3.0 * (1.0 - y * b.chi_eff));
  const double E_m = ellipe(b.m);
  const double K_m = ellipk(b.m);
  double hbpsip_pi_2 = std::fmod((2.0 / kPi) * bpsip + 1.0, 2.0);
  if (hbpsip_pi_2 < 0.0) {
    hbpsip_pi_2 += 2.0;
  }
  hbpsip_pi_2 -= 1.0;
  const double hpsip = K_m * hbpsip_pi_2;
  const auto jac = ellipj_from_u(hpsip, b.m);
  const double E_m_inc = dzeta_E_fact * (ellipe_inc(jac.am, b.m) - hbpsip_pi_2 * E_m);
  const double nusqY3mYm = mp.nu * b.sqY3mYm;
  double ellipPI_p_inc = 0.0;
  double ellipPI_m_inc = 0.0;
  if (nusqY3mYm != 0.0) {
    ellipPI_p_inc =
        (b.PI_fact_p * ellippi_inc(b.PI_arg_p, jac.am, b.m) - hbpsip_pi_2 * b.Pp) /
        nusqY3mYm;
    ellipPI_m_inc =
        (b.PI_fact_m * ellippi_inc(b.PI_arg_m, jac.am, b.m) - hbpsip_pi_2 * b.Pm) /
        nusqY3mYm;
  }
  const double dphiz = ellipPI_p_inc + ellipPI_m_inc;
  const double dzeta = E_m_inc + ellipPI_p_inc - ellipPI_m_inc;
  const double dchi = b.dchi_av - b.dchi_diff * (1.0 - 2.0 * jac.sn * jac.sn);
  const double costhL = (b.L + 0.5 * (b.chi_eff + mp.dmu * dchi)) / b.J;
  return {dphiz, dzeta, clamp(costhL, -1.0, 1.0)};
}

State derivatives_prec_avg(const State& v, const PNDerivatives& PN, double m1, double m2,
                           double sz_1, double sz_2, double sp2_1, double sp2_2) {
  // Right-hand side of the 8D radiation-reaction system. The first four entries
  // come from PNDerivatives; the last four evolve the precession invariant and
  // Euler-angle integration constants.
  const double y = v[0];
  const double e2 = clamp(v[1], 0.0, 1.0 - 1.0e-12);
  const double DJ2 = v[4];
  const auto mp = mass_params_from_m1_m2(m1, m2);
  const auto q = prec_avg_quantities_for_Dv(y, DJ2, m1, m2, sz_1, sz_2, sp2_1, sp2_2);
  const auto ders = PN.Dy_De2_Dl_Ddl(y, e2, q.dchi_prec_avg, q.dchi2_prec_avg,
                                     q.sperp2_prec_avg);
  const double DDJ2 = -ders[0] * q.dmudchi_prec_avg_m_dchi0 / (y * y);
  const double y6 = std::pow(y, 6);
  const double A = 1.0 - y * q.chi_eff;
  const double Dbpsip = 0.375 * kPi * A * y6 * q.sqY3mYm / q.K_m;
  const double P_pref = 0.75 * A / (mp.nu * q.K_m);
  const double Dphiz0 = y6 * (0.5 * q.J + P_pref * (q.Pp + q.Pm));
  const double Dzeta0 =
      y6 * (-0.25 * (2.0 * q.L + q.chi_eff + mp.dmu * q.dchi_prec_avg) -
            (1.5 / mp.nu) * (q.L + mp.nu * q.chi_eff) * A +
            P_pref * (q.Pp - q.Pm));
  const double D_fact = std::pow(1.0 - e2, 1.5) / mp.M;
  return {D_fact * ders[0], D_fact * ders[1], D_fact * ders[2], D_fact * ders[3],
          D_fact * DDJ2, D_fact * Dbpsip, D_fact * Dphiz0, D_fact * Dzeta0};
}

double polyval(const std::vector<double>& coefficients, double x) {
  double result = 0.0;
  for (std::size_t i = coefficients.size(); i-- > 0;) {
    result = result * x + coefficients[i];
  }
  return result;
}

double F_tLO_series_at_0(double x) {
  static const std::vector<double> coefs = vec({
      1.000000000000000, -0.1511627906976744, 0.2656836084021005,
      0.007463780007501875, 0.08800790590714085, 0.03153077124184580,
      0.04185392210761341, 0.02761371124737777, 0.02642686119635599,
      0.02145414169943131, 0.01926563742403364, 0.01677217450181226,
      0.01503898450331388, 0.01347003088516929, 0.01220348405026613,
      0.01110346195121149, 0.01016749330223870, 0.009352447459389354,
      0.008642114232972108, 0.008016709107501086, 0.007463457161231162,
      0.006970902964332086, 0.006530253812462707, 0.006134111124121483,
      0.005776451398258840, 0.005452229768143720, 0.005157232928828976,
      0.004887901117379935, 0.004641215172072290, 0.004414596725230578,
      0.004205833180162198, 0.004013015784562471, 0.003834490347048246,
      0.003668816871424952, 0.003514736646109113, 0.003371145103099870,
      0.003237069368178693, 0.003111649567641214, 0.002994123194217794,
      0.002883811964918853, 0.002780110725151045, 0.002682478039498533,
      0.002590428180410570, 0.002503524280447905, 0.002421372457429333,
      0.002343616756405161, 0.002269934780185545, 0.002200033902499887,
      0.002133647975961232, 0.002070534461714599, 0.002010471919657125});
  return polyval(coefs, x);
}

double F_tLO_series_at_1(double x) {
  static const std::vector<double> cns = vec({
      0.40941176470588236, 0.02286320645905421, 0.008142925951557094,
      0.004155878512401501, 0.0024847568765827364, 0.001633290511588348,
      0.0011455395465032605, 0.000842577094447224, 0.0006427195446513564,
      0.0005045715814217113, 0.00040544477831148924, 0.0003321116545345162,
      0.0002764636962467965, 0.00023331895675436905, 0.0001992473575067662,
      0.00017190919726595898, 0.00014966654751355843, 0.000131346469484909,
      0.00011609207361015848, 0.00010326617525262309, 9.238740896587563e-05,
      8.308691874613337e-05, 7.50784086790562e-05, 6.813705814151314e-05,
      6.20844345948999e-05, 5.677753690123865e-05, 5.210072977646673e-05,
      4.7959732146031274e-05, 4.427708468062032e-05, 4.0988696117937945e-05,
      3.80411855904995e-05, 3.538981870077757e-05, 3.2996890965966614e-05,
      3.083045152872919e-05, 2.886328796018351e-05, 2.707211306399751e-05,
      2.543690918050699e-05, 2.3940396192787112e-05, 2.256759736012561e-05,
      2.1305483020854193e-05, 2.014267666038337e-05, 1.906921121897629e-05,
      1.8076326095558655e-05, 1.7156297290322297e-05, 1.6302294667351972e-05,
      1.550826151746742e-05, 1.4768812541410176e-05, 1.407914711455732e-05});
  const double u = 1.0 - x;
  return ((48.0 / 19.0) * std::pow(425.0 / 304.0, 1181.0 / 2299.0)) *
         (1.0 - 1.4555165803216864 * std::sqrt(u) + u * polyval(cns, u)) *
         std::pow(x, -24.0 / 19.0) * std::pow(1.0 + (121.0 / 304.0) * x, -3480.0 / 2299.0);
}

double F_tLO_series(double x) {
  return x < 0.4 ? F_tLO_series_at_0(x) : F_tLO_series_at_1(x);
}

double tLO_func(double y, double e2, double m1, double m2) {
  const auto mp = mass_params_from_m1_m2(m1, m2);
  e2 = clamp(e2, 1.0e-14, 1.0 - 1.0e-12);
  return -(5.0 / 256.0) * (mp.M / mp.nu) * std::pow(y, -8.0) *
         std::pow(1.0 - e2, -0.5) * F_tLO_series(e2);
}

// Small vector helpers used by the explicit RK45 implementation below. We keep
// them local because State has fixed length and does not need a general linear
// algebra dependency.
State add_scaled(const State& y, const std::vector<std::pair<double, State>>& terms, double h) {
  State out = y;
  for (const auto& [a, k] : terms) {
    for (std::size_t i = 0; i < out.size(); ++i) {
      out[i] += h * a * k[i];
    }
  }
  return out;
}

double rms_norm_divided(const State& x, const std::array<double, 8>& scale) {
  double accum = 0.0;
  for (std::size_t i = 0; i < x.size(); ++i) {
    const double value = x[i] / scale[i];
    accum += value * value;
  }
  return std::sqrt(accum / static_cast<double>(x.size()));
}

template <typename Func>
double select_initial_step_rk45(Func&& f, const State& y0, const State& f0,
                                const std::vector<double>& rtol,
                                const std::vector<double>& atol) {
  std::array<double, 8> scale{};
  for (std::size_t i = 0; i < y0.size(); ++i) {
    scale[i] = atol[i] + std::abs(y0[i]) * rtol[i];
  }
  const double d0 = rms_norm_divided(y0, scale);
  const double d1 = rms_norm_divided(f0, scale);
  double h0 = (d0 < 1.0e-5 || d1 < 1.0e-5) ? 1.0e-6 : 0.01 * d0 / d1;
  State y1 = y0;
  for (std::size_t i = 0; i < y1.size(); ++i) {
    y1[i] += h0 * f0[i];
  }
  const State f1 = f(y1);
  State df{};
  for (std::size_t i = 0; i < df.size(); ++i) {
    df[i] = f1[i] - f0[i];
  }
  const double d2 = rms_norm_divided(df, scale) / h0;
  const double h1 =
      (d1 <= 1.0e-15 && d2 <= 1.0e-15)
          ? std::max(1.0e-6, h0 * 1.0e-3)
          : std::pow(0.01 / std::max(d1, d2), 1.0 / 5.0);
  return std::min(100.0 * h0, h1);
}

// Dormand-Prince dense-output coefficients. The adaptive solver stores these so
// later waveform calls can evaluate y(t), y'(t), and y''(t) inside each accepted
// ODE step without reintegrating.
std::array<State, 4> rk45_dense_coefficients(const State& k1, const State& k2,
                                             const State& k3, const State& k4,
                                             const State& k5, const State& k6,
                                             const State& k7) {
  static constexpr std::array<std::array<double, 4>, 7> P = {{
      {1.0, -8048581381.0 / 2820520608.0, 8663915743.0 / 2820520608.0,
       -12715105075.0 / 11282082432.0},
      {0.0, 0.0, 0.0, 0.0},
      {0.0, 131558114200.0 / 32700410799.0, -68118460800.0 / 10900136933.0,
       87487479700.0 / 32700410799.0},
      {0.0, -1754552775.0 / 470086768.0, 14199869525.0 / 1410260304.0,
       -10690763975.0 / 1880347072.0},
      {0.0, 127303824393.0 / 49829197408.0, -318862633887.0 / 49829197408.0,
       701980252875.0 / 199316789632.0},
      {0.0, -282668133.0 / 205662961.0, 2019193451.0 / 616988883.0,
       -1453857185.0 / 822651844.0},
      {0.0, 40617522.0 / 29380423.0, -110615467.0 / 29380423.0,
       69997945.0 / 29380423.0},
  }};
  const std::array<const State*, 7> K = {{&k1, &k2, &k3, &k4, &k5, &k6, &k7}};
  std::array<State, 4> q{};
  for (std::size_t stage = 0; stage < K.size(); ++stage) {
    for (std::size_t col = 0; col < q.size(); ++col) {
      for (std::size_t i = 0; i < q[col].size(); ++i) {
        q[col][i] += P[stage][col] * (*K[stage])[i];
      }
    }
  }
  return q;
}

struct Segment {
  // One accepted RK step. t0/t1 are shifted after integration so that the final
  // physical time agrees with the leading-order coalescence-time convention used
  // by the Python reference.
  double t0 = 0.0;
  double t1 = 0.0;
  double h = 0.0;
  double dense_h = 0.0;
  State y0{};
  State y1{};
  State f0{};
  State f1{};
  DenseOutput dense_output = DenseOutput::FastHermite;
  std::array<State, 4> dense_q{};

  // Evaluate the stored dense solution. Dormand-Prince uses the RK polynomial;
  // FastHermite uses endpoint values and derivatives only.
  State eval(double t) const {
    const double h_eval = dense_h != 0.0 ? dense_h : h;
    const double s = clamp((t - t0) / h_eval, 0.0, 1.0);
    if (dense_output == DenseOutput::DormandPrince) {
      const double s2 = s * s;
      const double s3 = s2 * s;
      const double s4 = s3 * s;
      State out = y0;
      for (std::size_t i = 0; i < out.size(); ++i) {
        out[i] += h_eval * (dense_q[0][i] * s + dense_q[1][i] * s2 +
                            dense_q[2][i] * s3 + dense_q[3][i] * s4);
      }
      return out;
    }
    const double s2 = s * s;
    const double s3 = s2 * s;
    const double h00 = 2.0 * s3 - 3.0 * s2 + 1.0;
    const double h10 = s3 - 2.0 * s2 + s;
    const double h01 = -2.0 * s3 + 3.0 * s2;
    const double h11 = s3 - s2;
    State out{};
    for (std::size_t i = 0; i < out.size(); ++i) {
      out[i] = h00 * y0[i] + h10 * h * f0[i] + h01 * y1[i] + h11 * h * f1[i];
    }
    return out;
  }

  // First derivative of the dense solution. This is needed for the
  // stationary-phase frequency d Phi / dt.
  State deriv(double t) const {
    const double h_eval = dense_h != 0.0 ? dense_h : h;
    const double s = clamp((t - t0) / h_eval, 0.0, 1.0);
    if (dense_output == DenseOutput::DormandPrince) {
      const double s2 = s * s;
      State out{};
      for (std::size_t i = 0; i < out.size(); ++i) {
        out[i] = dense_q[0][i] + 2.0 * dense_q[1][i] * s +
                 3.0 * dense_q[2][i] * s2 + 4.0 * dense_q[3][i] * s2 * s;
      }
      return out;
    }
    const double s2 = s * s;
    const double dh00 = 6.0 * s2 - 6.0 * s;
    const double dh10 = 3.0 * s2 - 4.0 * s + 1.0;
    const double dh01 = -6.0 * s2 + 6.0 * s;
    const double dh11 = 3.0 * s2 - 2.0 * s;
    State out{};
    for (std::size_t i = 0; i < out.size(); ++i) {
      out[i] = (dh00 * y0[i] + dh10 * h * f0[i] + dh01 * y1[i] + dh11 * h * f1[i]) / h;
    }
    return out;
  }

  // Second derivative of the dense solution. Its inverse square root gives the
  // SPA timescale T_spa.
  State second_deriv(double t) const {
    const double h_eval = dense_h != 0.0 ? dense_h : h;
    const double s = clamp((t - t0) / h_eval, 0.0, 1.0);
    if (dense_output == DenseOutput::DormandPrince) {
      const double s2 = s * s;
      State out{};
      for (std::size_t i = 0; i < out.size(); ++i) {
        out[i] =
            (2.0 * dense_q[1][i] + 6.0 * dense_q[2][i] * s +
             12.0 * dense_q[3][i] * s2) /
            h_eval;
      }
      return out;
    }
    const double d2h00 = 12.0 * s - 6.0;
    const double d2h10 = 6.0 * s - 4.0;
    const double d2h01 = -12.0 * s + 6.0;
    const double d2h11 = 6.0 * s - 2.0;
    State out{};
    for (std::size_t i = 0; i < out.size(); ++i) {
      out[i] =
          (d2h00 * y0[i] + d2h10 * h * f0[i] + d2h01 * y1[i] + d2h11 * h * f1[i]) /
          (h * h);
    }
    return out;
  }
};

Segment make_segment(double t0, double t1, double dense_h, const State& y0, const State& y1,
                     const State& f0, const State& f1, DenseOutput dense_output,
                     const std::array<State, 4>& dense_q = {}) {
  Segment segment;
  segment.t0 = t0;
  segment.t1 = t1;
  segment.h = t1 - t0;
  segment.dense_h = dense_output == DenseOutput::DormandPrince ? dense_h : segment.h;
  segment.y0 = y0;
  segment.y1 = y1;
  segment.f0 = f0;
  segment.f1 = f1;
  segment.dense_output = dense_output;
  segment.dense_q = dense_q;
  return segment;
}

double phase_value(const Segment& s, double t, int n, int m) {
  const State y = s.eval(t);
  return n * y[2] + (m - n) * y[3];
}

double phase_deriv(const Segment& s, double t, int n, int m) {
  const State dy = s.deriv(t);
  return n * dy[2] + (m - n) * dy[3];
}

double phase_second_deriv(const Segment& s, double t, int n, int m) {
  const State ddy = s.second_deriv(t);
  return n * ddy[2] + (m - n) * ddy[3];
}

struct InitialConditions {
  State v_ini{};
  double yf = 0.0;
  double sz_1 = 0.0;
  double sz_2 = 0.0;
  double chi_eff = 0.0;
  double sp2_1 = 0.0;
  double sp2_2 = 0.0;
  double s2_1 = 0.0;
  double s2_2 = 0.0;
  double cos_theta_JN = 1.0;
  double phi_JN = 0.0;
};

InitialConditions initial_conditions_for_RR_eqs(double f0_orb, double ff_orb, bool has_ff_orb,
                                                double e0, double m1, double m2,
                                                Vec3 spin0_1, Vec3 spin0_2,
                                                double inclination, double phi0,
                                                double phi_e0, double DJ2_tol) {
  // The public API uses f_22, but the ODE variables use orbital quantities.
  // At this point the caller has already divided f_22 by two.
  const auto mp = mass_params_from_m1_m2(m1, m2);
  const double e20 = e0 * e0;
  const double y0 = std::pow(2.0 * kPi * mp.M * f0_orb, 1.0 / 3.0) / safe_sqrt(1.0 - e20);
  const double yf_ISCO = std::pow(6.0, -0.5);
  double yf = yf_ISCO;
  if (has_ff_orb) {
    yf = std::pow(2.0 * kPi * mp.M * ff_orb, 1.0 / 3.0);
    yf = std::min(yf, yf_ISCO);
  }

  const double l0 = phi0;
  const double dl0 = phi0 - phi_e0;
  const Vec3 s0_1 = mp.mu1 * spin0_1;
  const Vec3 s0_2 = mp.mu2 * spin0_2;
  const double sz_1 = s0_1.z;
  const double sz_2 = s0_2.z;
  const double sp2_1 = s0_1.x * s0_1.x + s0_1.y * s0_1.y;
  const double sp2_2 = s0_2.x * s0_2.x + s0_2.y * s0_2.y;
  const double s2_1 = sp2_1 + sz_1 * sz_1;
  const double s2_2 = sp2_2 + sz_2 * sz_2;
  const double dchi0 = sz_1 - sz_2;

  // Rotate from the input frame into a frame aligned with the initial total
  // angular momentum J. The projection factors later use this J-frame line of
  // sight.
  const double cinc = std::cos(inclination);
  const double sinc = std::sin(inclination);
  const std::array<std::array<double, 3>, 3> R_spin = {
      std::array<double, 3>{cinc, 0.0, sinc},
      std::array<double, 3>{0.0, 1.0, 0.0},
      std::array<double, 3>{-sinc, 0.0, cinc}};
  const Vec3 s0p_1 = matvec(R_spin, s0_1);
  const Vec3 s0p_2 = matvec(R_spin, s0_2);
  const double L0 = mp.nu / y0;
  const Vec3 L0_hat{sinc, 0.0, cinc};
  const Vec3 L0_vec = L0 * L0_hat;
  const Vec3 J_vec = L0_vec + mp.mu1 * s0p_1 + mp.mu2 * s0p_2;

  const auto Jtrig = spherical_trig_from_vec(J_vec);
  const std::array<std::array<double, 3>, 3> R_J_to_z = {
      std::array<double, 3>{Jtrig.cos_theta * Jtrig.cos_phi,
                            Jtrig.cos_theta * Jtrig.sin_phi, -Jtrig.sin_theta},
      std::array<double, 3>{-Jtrig.sin_phi, Jtrig.cos_phi, 0.0},
      std::array<double, 3>{Jtrig.sin_theta * Jtrig.cos_phi,
                            Jtrig.sin_theta * Jtrig.sin_phi, Jtrig.cos_theta}};
  const Vec3 L0_hat_J_aligned = matvec(R_J_to_z, L0_hat);
  const auto Ltrig = spherical_trig_from_vec(L0_hat_J_aligned);
  const std::array<std::array<double, 3>, 3> R_L0p_to_xz = {
      std::array<double, 3>{Ltrig.cos_phi, Ltrig.sin_phi, 0.0},
      std::array<double, 3>{-Ltrig.sin_phi, Ltrig.cos_phi, 0.0},
      std::array<double, 3>{0.0, 0.0, 1.0}};
  const Vec3 k_hat{0.0, 0.0, -1.0};
  const Vec3 k_hat_Jframe = matvec(R_L0p_to_xz, matvec(R_J_to_z, k_hat));
  const auto [cos_theta_JN, phi_JN] = spherical_coords_from_vec(k_hat_Jframe);

  double DJ20 = 2.0 * (s0_1.x * s0_2.x + s0_1.y * s0_2.y);
  const double Vol = dot(cross(L0_hat, s0p_1), s0p_2);
  const double DJ20_OG = DJ20;
  double m0 = 0.0;
  double psip0 = 0.0;
  // DeltaJ^2 is implicit in the precession solution. Iterate the elliptic
  // integral correction until the initial invariant is self-consistent.
  for (int iJ = 0; iJ < 10; ++iJ) {
    const auto b = basic_prec_quantities(y0, DJ20, m1, m2, sz_1, sz_2, sp2_1, sp2_2);
    m0 = b.m;
    const double cos2am0 =
        b.dchi_diff != 0.0 ? clamp((b.dchi_av - dchi0) / b.dchi_diff, -1.0, 1.0) : 1.0;
    double am0 = 0.5 * std::acos(cos2am0);
    psip0 = ellipf_inc(am0, m0);
    if (Vol < 0.0) {
      psip0 = -psip0;
      am0 = -am0;
    }
    const double dDJ2 =
        (4.0 / 3.0) * mp.nu * b.sqY3mYm * ((32.0 + 28.0 * e20) / 5.0) *
        (y0 * y0 / (1.0 - y0 * b.chi_eff)) *
        (ellipe_inc(am0, m0) - ellipe(m0) * psip0 / ellipk(m0));
    const double DJ20_new = DJ20_OG - dDJ2;
    if (std::abs(DJ20_new - DJ20) < DJ2_tol) {
      DJ20 = DJ20_new;
      break;
    }
    DJ20 = DJ20_new;
  }

  const double bpsip0 = 0.5 * kPi * psip0 / ellipk(m0);
  const auto euler0 = precession_Euler_angles(bpsip0, y0, DJ20, m1, m2,
                                              sz_1, sz_2, sp2_1, sp2_2);
  const double phiz00 = -euler0.dphiz;
  const double zeta00 = -euler0.dzeta;

  return {{y0, e20, l0, dl0, DJ20, bpsip0, phiz00, zeta00},
          yf, sz_1, sz_2, sz_1 + sz_2, sp2_1, sp2_2, s2_1, s2_2,
          cos_theta_JN, phi_JN};
}

struct IntegrationResult {
  std::vector<Segment> segments;
  double physical_end_time = 0.0;
};

IntegrationResult solve_ivp_RR_eqs_t(const State& v_ini, const PNDerivatives& PN,
                                     double yf, double m1, double m2, double sz_1,
                                     double sz_2, double sp2_1, double sp2_2,
                                     const std::vector<double>& rtol_in,
                                     const std::vector<double>& atol_in,
                                     DenseOutput dense_output) {
  // Integrate forward in the leading-order time coordinate until y reaches the
  // requested final value. The accepted segments are later traversed many times
  // during waveform generation, so we store dense-output data with each step.
  std::vector<double> rtol = rtol_in;
  std::vector<double> atol = atol_in;
  rtol.resize(8, rtol.empty() ? 1e-10 : rtol.back());
  atol.resize(8, atol.empty() ? 1e-12 : atol.back());

  auto f = [&](const State& y) {
    return derivatives_prec_avg(y, PN, m1, m2, sz_1, sz_2, sp2_1, sp2_2);
  };

  double t = tLO_func(v_ini[0], v_ini[1], m1, m2);
  State y = v_ini;
  State k1 = f(y);
  double h = select_initial_step_rk45(f, y, k1, rtol, atol);
  std::vector<Segment> segments;
  segments.reserve(4096);
  bool reached_final = false;
  double physical_final_t = 0.0;
  State physical_final_state{};
  bool previous_rejected = false;

  constexpr int max_steps = 200000;
  for (int step = 0; step < max_steps; ++step) {
    if (reached_final || y[0] >= yf) {
      break;
    }
    const State k2 = f(add_scaled(y, {{1.0 / 5.0, k1}}, h));
    const State k3 = f(add_scaled(y, {{3.0 / 40.0, k1}, {9.0 / 40.0, k2}}, h));
    const State k4 = f(add_scaled(y, {{44.0 / 45.0, k1}, {-56.0 / 15.0, k2},
                                      {32.0 / 9.0, k3}}, h));
    const State k5 = f(add_scaled(y, {{19372.0 / 6561.0, k1}, {-25360.0 / 2187.0, k2},
                                      {64448.0 / 6561.0, k3}, {-212.0 / 729.0, k4}}, h));
    const State k6 = f(add_scaled(y, {{9017.0 / 3168.0, k1}, {-355.0 / 33.0, k2},
                                      {46732.0 / 5247.0, k3}, {49.0 / 176.0, k4},
                                      {-5103.0 / 18656.0, k5}}, h));
    State y5 = add_scaled(y, {{35.0 / 384.0, k1}, {500.0 / 1113.0, k3},
                              {125.0 / 192.0, k4}, {-2187.0 / 6784.0, k5},
                              {11.0 / 84.0, k6}}, h);
    const State k7 = f(y5);
    State y4 = add_scaled(y, {{5179.0 / 57600.0, k1}, {7571.0 / 16695.0, k3},
                              {393.0 / 640.0, k4}, {-92097.0 / 339200.0, k5},
                              {187.0 / 2100.0, k6}, {1.0 / 40.0, k7}}, h);
    const auto dense_q =
        dense_output == DenseOutput::DormandPrince
            ? rk45_dense_coefficients(k1, k2, k3, k4, k5, k6, k7)
            : std::array<State, 4>{};

    // Standard embedded RK error estimate: y5 is the 5th-order solution, y4 is
    // the companion 4th-order solution. The RMS norm is scaled componentwise by
    // atol + rtol * max(|old|, |new|).
    double err2 = 0.0;
    for (std::size_t i = 0; i < 8; ++i) {
      const double sc = atol[i] + rtol[i] * std::max(std::abs(y[i]), std::abs(y5[i]));
      const double e = (y5[i] - y4[i]) / sc;
      err2 += e * e;
    }
    const double err = std::sqrt(err2 / 8.0);
    if (err <= 1.0 || h <= 1.0e-12) {
      double factor =
          err == 0.0 ? 10.0 : std::min(10.0, 0.9 * std::pow(err, -0.2));
      if (previous_rejected) {
        factor = std::min(1.0, factor);
      }
      if (y5[0] >= yf) {
        const Segment trial = make_segment(t, t + h, h, y, y5, k1, k7, dense_output, dense_q);
        double lo = 0.0;
        double hi = 1.0;
        // The final accepted RK step may cross yf. Bisect inside the dense
        // polynomial so physical_end_time corresponds to exactly yf.
        for (int it = 0; it < 60; ++it) {
          const double mid = 0.5 * (lo + hi);
          if (trial.eval(t + mid * h)[0] >= yf) {
            hi = mid;
          } else {
            lo = mid;
          }
        }
        physical_final_t = t + hi * h;
        physical_final_state = trial.eval(physical_final_t);
        physical_final_state[0] = yf;
        segments.push_back(trial);
        reached_final = true;
        break;
      }
      segments.push_back(make_segment(t, t + h, h, y, y5, k1, k7, dense_output, dense_q));
      t += h;
      y = y5;
      k1 = k7;
      h *= factor;
      previous_rejected = false;
    } else {
      h *= std::max(0.2, 0.9 * std::pow(err, -0.2));
      previous_rejected = true;
    }
  }

  if (segments.empty()) {
    throw std::runtime_error("radiation-reaction integration produced no segments");
  }
  if (!reached_final) {
    throw std::runtime_error("radiation-reaction integration did not reach final y");
  }
  const double target_final =
      tLO_func(physical_final_state[0], physical_final_state[1], m1, m2);
  const double shift = target_final - physical_final_t;
  // Match the Python convention for the absolute time origin: the final state
  // lands at the leading-order time associated with that final (y, e^2).
  for (auto& s : segments) {
    s.t0 += shift;
    s.t1 += shift;
  }
  return {std::move(segments), physical_final_t + shift};
}

double compute_h0_pref(double m1, double m2, double dL) {
  return 4.0 * std::sqrt(kPi / 5.0) * (m1 * m2 / (m1 + m2)) / dL;
}

struct ModeData {
  // A "mode" here is local to a single RK segment. The same harmonic on two
  // adjacent segments appears as two ModeData entries because each segment has
  // its own dense polynomial and frequency support [w0, w1].
  int segment = 0;
  int m = 0;
  int p = 0;
  int n = 0;
  double w0 = 0.0;
  double w1 = 0.0;
  // Interpolation of the Newtonian eccentric amplitude N_{2m,p}(e^2).
  std::vector<double> n_interp;
  // Circular fast path: coarse monotone map omega -> t used as Newton's first
  // guess for stationary-time solves.
  std::vector<double> inverse_omega;
  std::vector<double> inverse_time;
};

struct ApcSegmentCache {
  // Projection/precession amplitudes sampled on the segment interpolation nodes.
  // m = -2 is obtained by conjugation symmetry, so only m=0 and m=+2 are stored.
  std::vector<std::array<C, 2>> m0;
  std::vector<std::array<C, 2>> m2;
};

std::array<C, 2> zero_amp() {
  return {C{0.0, 0.0}, C{0.0, 0.0}};
}

} // namespace

double f22_isco_frequency_hz(double mass1_solar, double mass2_solar) {
  // For this inspiral-only model we use the Schwarzschild ISCO as a simple
  // termination scale. y = sqrt(M/r), so y_ISCO = 1/sqrt(6), and f_22 = Omega/pi.
  const double total_mass_seconds = kTSunSeconds * (mass1_solar + mass2_solar);
  const double y_isco = std::pow(6.0, -0.5);
  return (y_isco * y_isco * y_isco) / (kPi * total_mass_seconds);
}

double estimate_f22_start_for_duration_newtonian(double mass1_solar, double mass2_solar,
                                                double f22_end_hz, double duration_seconds,
                                                double min_f22_start_hz) {
  if (!(mass1_solar > 0.0) || !(mass2_solar > 0.0)) {
    throw std::invalid_argument("masses must be positive");
  }
  if (!(f22_end_hz > 0.0) || !(duration_seconds >= 0.0) || !(min_f22_start_hz > 0.0)) {
    throw std::invalid_argument("frequencies and duration must be positive");
  }
  const double total_mass = mass1_solar + mass2_solar;
  const double eta = mass1_solar * mass2_solar / (total_mass * total_mass);
  const double chirp_mass_seconds = kTSunSeconds * total_mass * std::pow(eta, 3.0 / 5.0);
  // Leading-order circular chirp time scales as f^{-8/3}; invert that relation
  // to get a cheap start-frequency estimate.
  const double end_term = std::pow(f22_end_hz, -8.0 / 3.0);
  const double duration_term =
      (256.0 / 5.0) * std::pow(kPi * chirp_mass_seconds, 5.0 / 3.0) * duration_seconds;
  return std::max(min_f22_start_hz, std::pow(end_term + duration_term, -3.0 / 8.0));
}

double choose_f22_start_for_duration(Parameters parameters, double duration_seconds,
                                     double min_f22_start_hz) {
  if (!(duration_seconds >= 0.0) || !(min_f22_start_hz > 0.0)) {
    throw std::invalid_argument("duration and minimum start frequency must be positive");
  }
  if (!(parameters.mass1 > 0.0) || !(parameters.mass2 > 0.0)) {
    throw std::invalid_argument("mass1 and mass2 must be positive");
  }
  if (!(parameters.f22_end > 0.0)) {
    parameters.f22_end = f22_isco_frequency_hz(parameters.mass1, parameters.mass2);
  }
  const double f_end = parameters.f22_end;
  const double upper = 0.95 * f_end;
  double guess = estimate_f22_start_for_duration_newtonian(
      parameters.mass1, parameters.mass2, f_end, duration_seconds, min_f22_start_hz);
  guess = clamp(guess, min_f22_start_hz, upper);

  // The model duration is not exactly the leading-order circular duration,
  // especially for eccentric/precessing systems. Use the cheap estimate only as
  // a bracket seed, then evaluate full validation models during bisection.
  auto duration_for = [&](double f22_start_hz) {
    Parameters p = with_preset(parameters, ParameterPreset::Validation);
    p.f22_start = f22_start_hz;
    try {
      Model model(p);
      return model.end_time() - model.start_time();
    } catch (const std::runtime_error&) {
      return 0.0;
    }
  };

  double lo = guess;
  double hi = guess;
  double duration = duration_for(guess);
  if (duration > duration_seconds) {
    while (duration > duration_seconds && hi < upper) {
      lo = hi;
      hi = std::min(upper, 1.35 * hi);
      duration = duration_for(hi);
    }
  } else {
    while (duration < duration_seconds && lo > min_f22_start_hz) {
      hi = lo;
      lo = std::max(min_f22_start_hz, lo / 1.35);
      duration = duration_for(lo);
    }
  }

  for (int i = 0; i < 24; ++i) {
    const double mid = 0.5 * (lo + hi);
    const double mid_duration = duration_for(mid);
    if (mid_duration > duration_seconds) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  return 0.5 * (lo + hi);
}

void apply_extrinsic_transform_uniform(double f_min_hz, double delta_f_hz,
                                       const ExtrinsicTransform& transform,
                                       FrequencyWaveform& waveform) {
  if (!(delta_f_hz > 0.0)) {
    throw std::invalid_argument("delta_f_hz must be positive");
  }
  if (waveform.plus.size() != waveform.cross.size()) {
    throw std::invalid_argument("plus and cross arrays must have the same size");
  }
  const double c2p = std::cos(2.0 * transform.polarization_angle_radians);
  const double s2p = std::sin(2.0 * transform.polarization_angle_radians);
  for (std::size_t i = 0; i < waveform.plus.size(); ++i) {
    const double f = f_min_hz + delta_f_hz * static_cast<double>(i);
    // A time shift in the Fourier domain is exp(-2 pi i f dt). The constant
    // phase and distance/amplitude scale multiply both polarizations.
    const C phase = transform.amplitude_scale *
                    std::polar(1.0, transform.phase_shift_radians -
                                         2.0 * kPi * f * transform.time_shift_seconds);
    const C hp = waveform.plus[i];
    const C hc = waveform.cross[i];
    waveform.plus[i] = phase * (c2p * hp + s2p * hc);
    waveform.cross[i] = phase * (-s2p * hp + c2p * hc);
  }
}

FrequencyWaveform transformed_uniform(const FrequencyWaveform& waveform, double f_min_hz,
                                      double delta_f_hz,
                                      const ExtrinsicTransform& transform) {
  FrequencyWaveform out = waveform;
  apply_extrinsic_transform_uniform(f_min_hz, delta_f_hz, transform, out);
  return out;
}

namespace {

constexpr double kCircularEccentricity = 1.0e-14;
constexpr double kLowMassSolar = 3000.0;
constexpr double kHighEccentricity = 0.30;
constexpr double kVeryHighEccentricity = 0.60;

// Two tolerance sets are used throughout the presets. "Fast" matches the
// cheaper production ODE integration; "Strict" is the validation setting and is
// also used in production for regimes where the 40-system suite showed that the
// ODE error matters.
const std::vector<double> kFastRtol = {1e-10, 1e-10, 1e-12, 1e-12,
                                       1e-10, 1e-12, 1e-12, 1e-12};
const std::vector<double> kFastAtol = {1e-12, 1e-12, 1e-8, 1e-8,
                                       1e-12, 1e-8, 1e-6, 1e-6};
const std::vector<double> kStrictRtol = {1e-12, 1e-12, 1e-14, 1e-14,
                                         1e-12, 1e-14, 1e-14, 1e-14};
const std::vector<double> kStrictAtol = {1e-14, 1e-14, 1e-10, 1e-10,
                                         1e-14, 1e-10, 1e-8, 1e-8};

bool is_effectively_circular(double e_start) {
  return e_start <= kCircularEccentricity;
}

double total_mass_solar(const Parameters& p) {
  return p.mass1 + p.mass2;
}

void set_fast_rr_tolerances(Parameters& p) {
  p.rr_sol_rtol = kFastRtol;
  p.rr_sol_atol = kFastAtol;
}

void set_strict_rr_tolerances(Parameters& p) {
  p.rr_sol_rtol = kStrictRtol;
  p.rr_sol_atol = kStrictAtol;
}

void set_common_preset_options(Parameters& p, DenseOutput dense_output, bool interpolate) {
  p.dense_output = dense_output;
  p.interpolate_amplitudes = interpolate;
  p.amplitude_interp_nodes = 16;
}

void allow_fast_waveform_options(Parameters& p) {
  p.circular_fast_path = true;
  p.circular_spa_iterations = 3;
  p.num_threads = std::max(1, p.num_threads);
  p.spa_frequency_rtol = 1.0e-12;
}

void force_reference_waveform_options(Parameters& p) {
  p.circular_fast_path = false;
  p.circular_spa_iterations = 8;
  p.num_threads = 1;
  p.spa_frequency_rtol = 1.0e-14;
}

bool needs_strict_production_rr(const Parameters& p) {
  // These thresholds are empirical but intentionally visible: they are the
  // regimes where the production waveform needs validation-like ODE tolerances
  // to stay below the 1e-3 white-noise mismatch target.
  if (is_effectively_circular(p.e_start)) {
    return false;
  }
  if (total_mass_solar(p) < kLowMassSolar && p.e_start >= 0.10) {
    return true;
  }
  return p.e_start >= kHighEccentricity;
}

bool needs_tight_high_e_mode_cutoff(const Parameters& p) {
  // High eccentricity spreads power over many p-harmonics. The tighter cutoff
  // keeps more of that support while leaving easier systems on the faster
  // 2e-4 truncation.
  const double mass = total_mass_solar(p);
  const bool low_mass_mid_eccentricity = mass < 5000.0 && p.e_start < 0.45;
  return p.e_start >= kVeryHighEccentricity || low_mass_mid_eccentricity;
}

void set_production_accuracy(Parameters& p) {
  // This is the production compromise: cached amplitudes and fast waveform
  // options are always on, but the harmonic cutoff and RR tolerances adapt to
  // eccentricity and mass. Validation uses a separate preset rather than these
  // branches.
  if (is_effectively_circular(p.e_start)) {
    p.amplitude_tol = 1.0e-3;
    p.amplitude_pmax = 4;
  } else if (total_mass_solar(p) < kLowMassSolar && p.e_start >= 0.10) {
    p.amplitude_tol = 1.0e-4;
    p.amplitude_pmax = 200;
  } else if (p.e_start >= kHighEccentricity) {
    p.amplitude_tol = needs_tight_high_e_mode_cutoff(p) ? 1.0e-4 : 2.0e-4;
    p.amplitude_pmax = 200;
  } else {
    p.amplitude_tol = 1.0e-4;
    p.amplitude_pmax = 160;
  }

  if (needs_strict_production_rr(p)) {
    set_strict_rr_tolerances(p);
  } else {
    set_fast_rr_tolerances(p);
  }
}

} // namespace

void apply_preset(Parameters& p, ParameterPreset preset) {
  switch (preset) {
  case ParameterPreset::Fast:
    set_common_preset_options(p, DenseOutput::FastHermite, true);
    allow_fast_waveform_options(p);
    p.amplitude_tol = 1.0e-3;
    p.amplitude_pmax = 100;
    set_fast_rr_tolerances(p);
    break;
  case ParameterPreset::Production:
    set_common_preset_options(p, DenseOutput::DormandPrince, true);
    allow_fast_waveform_options(p);
    set_production_accuracy(p);
    break;
  case ParameterPreset::Validation:
    set_common_preset_options(p, DenseOutput::DormandPrince, false);
    force_reference_waveform_options(p);
    p.amplitude_tol = 1.0e-6;
    p.amplitude_pmax = 200;
    set_strict_rr_tolerances(p);
    break;
  }
}

Parameters with_preset(Parameters parameters, ParameterPreset preset) {
  apply_preset(parameters, preset);
  return parameters;
}

struct Model::Impl {
  explicit Impl(Parameters input) : params(std::move(input)) {
    if (!(params.mass1 > 0.0) || !(params.mass2 > 0.0)) {
      throw std::invalid_argument("mass1 and mass2 must be positive");
    }

    // Convert external solar-mass inputs to geometric seconds once. Everything
    // after this point follows the units of the Python equations.
    m1 = kTSunSeconds * params.mass1;
    m2 = kTSunSeconds * params.mass2;
    Vec3 spin0_1{params.spin1x, params.spin1y, params.spin1z};
    Vec3 spin0_2{params.spin2x, params.spin2y, params.spin2z};
    double phi0 = params.phi_start;
    double phi_e0 = params.mean_anomaly_start;
    q1 = params.q1;
    q2 = params.q2;

    // The pyEFPE formulas assume body 1 is the heavier object. Preserve the
    // physical binary by swapping spins/quadrupoles along with the masses and
    // shifting phases by pi.
    if (m2 > m1) {
      std::swap(m1, m2);
      std::swap(q1, q2);
      std::swap(spin0_1, spin0_2);
      phi0 += kPi;
      phi_e0 += kPi;
    }

    const double f0_orb = 0.5 * params.f22_start;
    const bool has_ff_orb = params.f22_end > 0.0;
    const double ff_orb = has_ff_orb ? 0.5 * params.f22_end : 0.0;
    const double dL = kMpcSeconds * params.distance;

    // Initial conditions are the only place where user-frame geometry enters
    // the ODE state. The rest of the model evolves scalar PN/precession
    // variables and reconstructs sky-frame polarizations later.
    const auto ic = initial_conditions_for_RR_eqs(
        f0_orb, ff_orb, has_ff_orb, params.e_start, m1, m2, spin0_1, spin0_2,
        params.inclination, phi0, phi_e0, params.dj2_tol);
    v_ini = ic.v_ini;
    yf = ic.yf;
    sz_1 = ic.sz_1;
    sz_2 = ic.sz_2;
    chi_eff = ic.chi_eff;
    sp2_1 = ic.sp2_1;
    sp2_2 = ic.sp2_2;
    s2_1 = ic.s2_1;
    s2_2 = ic.s2_2;
    cos_theta_JN = ic.cos_theta_JN;
    phi_JN = ic.phi_JN;

    pn = std::make_unique<PNDerivatives>(m1, m2, chi_eff, s2_1, s2_2, q1, q2,
                                         params.pn_phase_order, params.pn_spin_order);

    // Precompute the projection from the co-precessing tensor modes onto plus
    // and cross for the fixed line of sight. This avoids repeating spherical
    // harmonic bookkeeping at every frequency sample.
    const auto y2mp = compute_m2_Y2(cos_theta_JN, phi_JN);
    std::array<C, 5> y2mp_mod{};
    for (std::size_t i = 0; i < 5; ++i) {
      y2mp_mod[i] = std::conj(y2mp[4 - i]);
    }
    y2mp_mod[1] = -y2mp_mod[1];
    y2mp_mod[3] = -y2mp_mod[3];
    for (std::size_t i = 0; i < 5; ++i) {
      const C Ap = 0.5 * (y2mp[i] + y2mp_mod[i]);
      const C Ac = C{0.0, -0.5} * (y2mp[i] - y2mp_mod[i]);
      Apc_proj[i] = {Ap, Ac};
    }

    h0_pref = compute_h0_pref(m1, m2, dL);
    // The expensive setup work happens here: integrate the radiation reaction
    // once, then build all caches needed for fast repeated waveform evaluation.
    const auto integration =
        solve_ivp_RR_eqs_t(v_ini, *pn, yf, m1, m2, sz_1, sz_2, sp2_1, sp2_2,
                           params.rr_sol_rtol, params.rr_sol_atol,
                           params.dense_output);
    segments = integration.segments;
    physical_end_time = integration.physical_end_time;
    ak_SUA = compute_ak_SUA(params.sua_kmax);
    initialize_interpolation_nodes();
    // The order matters: amplitude caches depend on interpolation nodes, modes
    // depend on the ODE segments, and the circular inverse cache depends on the
    // final segment-local mode list.
    build_apc_cache();
    build_modes();
    build_circular_inverse_cache();
    build_mode_amplitude_cache();
  }

  void initialize_interpolation_nodes() {
    const int n = std::max(2, params.amplitude_interp_nodes);
    interp_nodes.resize(static_cast<std::size_t>(n));
    interp_weights.assign(static_cast<std::size_t>(n), 1.0);
    // Chebyshev-Lobatto nodes on [0, 1]. They include both segment endpoints,
    // which is helpful because stationary times often sit close to support
    // boundaries.
    for (int i = 0; i < n; ++i) {
      interp_nodes[static_cast<std::size_t>(i)] =
          0.5 * (1.0 - std::cos(kPi * static_cast<double>(i) / static_cast<double>(n - 1)));
    }
    // Barycentric weights for polynomial interpolation on those fixed nodes.
    // n is small (default 16), so the O(n^2) setup is negligible.
    for (int i = 0; i < n; ++i) {
      double w = 1.0;
      const double xi = interp_nodes[static_cast<std::size_t>(i)];
      for (int j = 0; j < n; ++j) {
        if (i == j) {
          continue;
        }
        w *= xi - interp_nodes[static_cast<std::size_t>(j)];
      }
      interp_weights[static_cast<std::size_t>(i)] = 1.0 / w;
    }
  }

  C barycentric_complex(const std::vector<std::array<C, 2>>& values, double x,
                        std::size_t pol) const {
    // Stable Lagrange interpolation. If x is exactly one of the nodes, return
    // the stored value to avoid the removable singularity.
    C numerator{0.0, 0.0};
    double denominator = 0.0;
    for (std::size_t i = 0; i < interp_nodes.size(); ++i) {
      const double dx = x - interp_nodes[i];
      if (std::abs(dx) < 1.0e-14) {
        return values[i][pol];
      }
      const double w = interp_weights[i] / dx;
      numerator += w * values[i][pol];
      denominator += w;
    }
    return numerator / denominator;
  }

  double barycentric_real(const std::vector<double>& values, double x) const {
    // Same barycentric formula for the scalar Newtonian eccentric amplitude.
    double numerator = 0.0;
    double denominator = 0.0;
    for (std::size_t i = 0; i < interp_nodes.size(); ++i) {
      const double dx = x - interp_nodes[i];
      if (std::abs(dx) < 1.0e-14) {
        return values[i];
      }
      const double w = interp_weights[i] / dx;
      numerator += w * values[i];
      denominator += w;
    }
    return numerator / denominator;
  }

  bool use_circular_fast_path() const {
    // Treat only exactly/near-zero eccentricity as circular. Small nonzero e is
    // still allowed to generate eccentric sidebands.
    return params.circular_fast_path && is_effectively_circular(params.e_start);
  }

  double linear_time_guess(const Segment& seg, const ModeData& mode, double omega) const {
    // Cheap first guess from assuming omega(t) is linear over one RK segment.
    return seg.t0 + (omega - mode.w0) * (seg.t1 - seg.t0) / (mode.w1 - mode.w0);
  }

  double circular_inverse_time_guess(const Segment& seg, const ModeData& mode, double omega) const {
    // In circular templates each segment has only one relevant harmonic, so a
    // tiny cached inverse map gives a much better Newton starting point than the
    // generic linear guess.
    if (!use_circular_fast_path() || mode.inverse_omega.size() < 2) {
      return linear_time_guess(seg, mode, omega);
    }

    const auto upper = std::lower_bound(mode.inverse_omega.begin(), mode.inverse_omega.end(), omega);
    if (upper == mode.inverse_omega.begin()) {
      return mode.inverse_time.front();
    }
    if (upper == mode.inverse_omega.end()) {
      return mode.inverse_time.back();
    }

    const std::size_t j = static_cast<std::size_t>(std::distance(mode.inverse_omega.begin(), upper));
    const double w0 = mode.inverse_omega[j - 1];
    const double w1 = mode.inverse_omega[j];
    if (!(w1 > w0)) {
      return linear_time_guess(seg, mode, omega);
    }
    const double x = (omega - w0) / (w1 - w0);
    return mode.inverse_time[j - 1] + x * (mode.inverse_time[j] - mode.inverse_time[j - 1]);
  }

  double initial_time_guess(const Segment& seg, const ModeData& mode, double omega,
                            double previous_guess) const {
    // For sorted frequency grids the stationary time moves monotonically. Reuse
    // the previous solve as a warm start, except in the circular path where the
    // inverse cache is consistently better.
    if (use_circular_fast_path()) {
      return circular_inverse_time_guess(seg, mode, omega);
    }
    if (previous_guess > seg.t0 && previous_guess < seg.t1) {
      return previous_guess;
    }
    return linear_time_guess(seg, mode, omega);
  }

  std::vector<ModeOrder> mode_orders_for_segment(double e2, int pmax) const {
    // The circular waveform has no eccentric sidebands, so keeping only
    // (m_abs=2, p=0) saves both setup and frequency-domain traversal.
    if (use_circular_fast_path()) {
      return {{2, 0}};
    }
    return Newtonian_orders_needed(e2, pmax, params.amplitude_tol);
  }

  void build_modes() {
    modes.clear();
    modes_by_segment.assign(segments.size(), {});
    // Eccentricity decays along the inspiral, so the previous segment's support
    // bounds the harmonic search for the next segment.
    int pmax_ceiling = params.amplitude_pmax;
    for (std::size_t iseg = 0; iseg < segments.size(); ++iseg) {
      const double e2 = std::max(segments[iseg].y0[1], 0.0);
      const std::vector<ModeOrder> needed = mode_orders_for_segment(e2, pmax_ceiling);
      int next_pmax_ceiling = 1;
      for (const auto& mo : needed) {
        // pyEFPE labels Newtonian amplitudes by (m_abs, p). The stationary
        // phase uses n = |p + m_abs| and assigns the sign to m.
        const int raw_n = mo.p + mo.m_abs;
        const int m = raw_n >= 0 ? mo.m_abs : -mo.m_abs;
        const int n = std::abs(raw_n);
        const double w0 = phase_deriv(segments[iseg], segments[iseg].t0, n, m);
        const double w1 = phase_deriv(segments[iseg], segments[iseg].t1, n, m);
        // A mode contributes only where its instantaneous angular frequency is
        // positive and increasing over the segment. Otherwise there is no
        // well-behaved stationary point for positive Fourier frequencies.
        if (w1 > w0 && w1 > 0.0) {
          const int idx = static_cast<int>(modes.size());
          ModeData mode;
          mode.segment = static_cast<int>(iseg);
          mode.m = m;
          mode.p = mo.p;
          mode.n = n;
          mode.w0 = w0;
          mode.w1 = w1;
          modes.push_back(std::move(mode));
          modes_by_segment[iseg].push_back(idx);
        }
        next_pmax_ceiling = std::max(next_pmax_ceiling, 1 + std::abs(mo.p));
      }
      pmax_ceiling = next_pmax_ceiling;
    }
  }

  void build_apc_cache() {
    apc_cache.clear();
    apc_cache.resize(segments.size());
    if (!params.interpolate_amplitudes) {
      return;
    }
    // A_pc contains the precession/projection part of the amplitude. It changes
    // slowly compared with the carrier phase, so interpolation is accurate and
    // avoids recomputing Wigner rotations at every frequency bin.
    for (std::size_t iseg = 0; iseg < segments.size(); ++iseg) {
      auto& cache = apc_cache[iseg];
      cache.m0.resize(interp_nodes.size());
      cache.m2.resize(interp_nodes.size());
      const Segment& seg = segments[iseg];
      for (std::size_t inode = 0; inode < interp_nodes.size(); ++inode) {
        const double t = seg.t0 + seg.h * interp_nodes[inode];
        const State state = seg.eval(t);
        cache.m0[inode] = compute_Apc_prec_from_state_exact(state, 0);
        cache.m2[inode] = compute_Apc_prec_from_state_exact(state, 2);
      }
    }
  }

  void build_mode_amplitude_cache() {
    if (!params.interpolate_amplitudes) {
      return;
    }
    // N_{2m,p}(e^2) is the Newtonian eccentric Fourier coefficient. It depends
    // only on e^2 and the mode labels, so cache it separately from A_pc.
    for (auto& mode : modes) {
      mode.n_interp.resize(interp_nodes.size());
      const Segment& seg = segments[static_cast<std::size_t>(mode.segment)];
      for (std::size_t inode = 0; inode < interp_nodes.size(); ++inode) {
        const State state = seg.eval(seg.t0 + seg.h * interp_nodes[inode]);
        mode.n_interp[inode] = compute_N2m_from_state_exact(state, mode);
      }
    }
  }

  void build_circular_inverse_cache() {
    if (!use_circular_fast_path()) {
      return;
    }
    // Eight nodes are enough because this cache is only a first guess; Newton
    // refinement still enforces the stationary-phase equation.
    constexpr int nnode = 8;
    for (auto& mode : modes) {
      mode.inverse_omega.resize(nnode);
      mode.inverse_time.resize(nnode);
      const Segment& seg = segments[static_cast<std::size_t>(mode.segment)];
      for (int inode = 0; inode < nnode; ++inode) {
        const double x = static_cast<double>(inode) / static_cast<double>(nnode - 1);
        const double t = seg.t0 + x * seg.h;
        mode.inverse_time[static_cast<std::size_t>(inode)] = t;
        mode.inverse_omega[static_cast<std::size_t>(inode)] = phase_deriv(seg, t, mode.n, mode.m);
      }
    }
  }

  std::size_t find_segment(double t) const {
    if (t <= segments.front().t0) {
      return 0;
    }
    if (t >= segments.back().t1) {
      return segments.size() - 1;
    }
    std::size_t lo = 0;
    std::size_t hi = segments.size();
    while (lo + 1 < hi) {
      const std::size_t mid = (lo + hi) / 2;
      if (segments[mid].t0 <= t) {
        lo = mid;
      } else {
        hi = mid;
      }
    }
    return lo;
  }

  State eval_state(double t) const {
    return segments[find_segment(t)].eval(t);
  }

  std::pair<std::array<C, 5>, std::array<C, 5>> compute_D2_mp_from_state(const State& yv) const {
    // Reconstruct the Wigner-D pieces that rotate the co-precessing l=2 modes
    // into the observer frame. Non-precessing systems take the cheap branch with
    // zero Euler angles.
    double y = yv[0];
    double e2 = std::max(yv[1], 0.0);
    double DJ2 = yv[4];
    double bpsip = yv[5];
    double phiz = yv[6];
    double zeta = yv[7];
    double costhL = 1.0;
    if (sp2_1 != 0.0 || sp2_2 != 0.0) {
      const auto euler = precession_Euler_angles(bpsip, y, DJ2, m1, m2,
                                                 sz_1, sz_2, sp2_1, sp2_2);
      phiz += euler.dphiz;
      zeta += euler.dzeta;
      costhL = euler.costhL;
    } else {
      phiz = 0.0;
      zeta = 0.0;
    }
    auto D = compute_necessary_Wigner_D2(phiz, costhL, zeta);
    // The Newtonian quadrupole amplitude contributes the common (1-e^2) y^2
    // factor here; the eccentric harmonic coefficient N_{2m,p} is applied
    // later because it is mode-dependent.
    const double omega_factor = (1.0 - e2) * y * y;
    for (C& v : D.first) {
      v *= omega_factor;
    }
    for (C& v : D.second) {
      v *= omega_factor;
    }
    return D;
  }

  std::pair<std::array<C, 5>, std::array<C, 5>> compute_D2_mp_exact(double t) const {
    return compute_D2_mp_from_state(eval_state(t));
  }

  std::array<C, 2> project(const std::array<C, 5>& D) const {
    // Contract the five m' components with the line-of-sight projection factors
    // computed in the constructor. The result is [A_plus, A_cross].
    std::array<C, 2> out = zero_amp();
    for (std::size_t i = 0; i < 5; ++i) {
      out[0] += D[i] * Apc_proj[i][0];
      out[1] += D[i] * Apc_proj[i][1];
    }
    return out;
  }

  std::array<C, 2> compute_Apc_prec(double t, int m) const {
    // Production runs interpolate the expensive projection/precession amplitude;
    // validation runs call the exact reconstruction every time.
    if (!params.interpolate_amplitudes) {
      return compute_Apc_prec_from_state_exact(eval_state(t), m);
    }
    const std::size_t iseg = find_segment(t);
    const Segment& seg = segments[iseg];
    const auto& cache = apc_cache[iseg];
    const double x = clamp((t - seg.t0) / seg.h, 0.0, 1.0);
    std::array<C, 2> A = zero_amp();
    if (m == 0) {
      A = {barycentric_complex(cache.m0, x, 0), barycentric_complex(cache.m0, x, 1)};
    } else if (std::abs(m) == 2) {
      A = {barycentric_complex(cache.m2, x, 0), barycentric_complex(cache.m2, x, 1)};
      if (m == -2) {
        A[0] = std::conj(A[0]);
        A[1] = std::conj(A[1]);
      }
    }
    return A;
  }

  std::array<C, 2> compute_Apc_prec_from_state_exact(const State& state, int m) const {
    // Exact here means "no amplitude interpolation." The underlying equations
    // are still the same approximate PN/pyEFPE model.
    const auto D = compute_D2_mp_from_state(state);
    std::array<C, 2> A = zero_amp();
    if (m == 0) {
      A = project(D.second);
    } else if (std::abs(m) == 2) {
      A = project(D.first);
      if (m == -2) {
        A[0] = std::conj(A[0]);
        A[1] = std::conj(A[1]);
      }
    }
    return A;
  }

  double compute_N2m(double t, const ModeData& mode) const {
    // Use the segment-local N_{2m,p} cache only if the time lies inside that
    // segment. SUA stencil points can fall outside the segment near boundaries.
    if (!params.interpolate_amplitudes || mode.n_interp.empty()) {
      return compute_N2m_from_state_exact(eval_state(t), mode);
    }
    const Segment& seg = segments[static_cast<std::size_t>(mode.segment)];
    if (t < seg.t0 || t > seg.t1) {
      return compute_N2m_from_state_exact(eval_state(t), mode);
    }
    return barycentric_real(mode.n_interp, clamp((t - seg.t0) / seg.h, 0.0, 1.0));
  }

  double compute_N2m_from_state_exact(const State& state, const ModeData& mode) const {
    const double e2 = std::max(state[1], 0.0);
    if (mode.m == 0) {
      return N20_Newtonian(mode.p, e2);
    }
    if (std::abs(mode.m) == 2) {
      return N22_Newtonian(mode.p, e2);
    }
    return 0.0;
  }

  std::array<C, 2> compute_amplitudes(double t, const ModeData& mode) const {
    // Fast path: interpolate both A_pc and N on the same segment nodes and
    // multiply them. Fallbacks preserve correctness at segment/SUA boundaries.
    if (params.interpolate_amplitudes && !mode.n_interp.empty()) {
      const Segment& seg = segments[static_cast<std::size_t>(mode.segment)];
      if (t >= seg.t0 && t <= seg.t1) {
        const double x = clamp((t - seg.t0) / seg.h, 0.0, 1.0);
        const auto& cache = apc_cache[static_cast<std::size_t>(mode.segment)];
        std::array<C, 2> A = zero_amp();
        if (mode.m == 0) {
          A = {barycentric_complex(cache.m0, x, 0), barycentric_complex(cache.m0, x, 1)};
        } else if (std::abs(mode.m) == 2) {
          A = {barycentric_complex(cache.m2, x, 0), barycentric_complex(cache.m2, x, 1)};
          if (mode.m == -2) {
            A[0] = std::conj(A[0]);
            A[1] = std::conj(A[1]);
          }
        }
        const double N = barycentric_real(mode.n_interp, x);
        A[0] *= N;
        A[1] *= N;
        return A;
      }
    }
    std::array<C, 2> A = compute_Apc_prec(t, mode.m);
    const double N = compute_N2m(t, mode);
    A[0] *= N;
    A[1] *= N;
    return A;
  }

  std::array<C, 2> SUA_amplitudes(double t, const ModeData& mode, double T_SPA) const {
    // Shifted Uniform Asymptotics samples the slow amplitude around t by
    // multiples of T_spa. Near the beginning/end of the waveform the stencil
    // would leave the available ODE support, so we safely fall back to the
    // unshifted SPA amplitude.
    const bool in_range =
        (t - params.sua_kmax * T_SPA >= start_time()) &&
        (t + params.sua_kmax * T_SPA <= end_time());
    if (!in_range) {
      return compute_amplitudes(t, mode);
    }
    std::array<C, 2> out = zero_amp();
    for (int k = -params.sua_kmax; k <= params.sua_kmax; ++k) {
      const auto A = compute_amplitudes(t + static_cast<double>(k) * T_SPA, mode);
      const C coeff = ak_SUA[static_cast<std::size_t>(std::abs(k))];
      out[0] += coeff * A[0];
      out[1] += coeff * A[1];
    }
    return out;
  }

  bool stationary_time_for_mode(const ModeData& mode, double omega, double& t_spa,
                                double& phase, double& T_spa) const {
    // Reference stationary-time solve. It starts from a linear guess and uses a
    // bracketed Newton iteration. The bracket makes the solve robust even when
    // the local dense polynomial is not perfectly linear in frequency.
    if (omega < mode.w0 || omega > mode.w1) {
      return false;
    }
    const Segment& seg = segments[static_cast<std::size_t>(mode.segment)];
    double lo = seg.t0;
    double hi = seg.t1;
    double t = lo + (omega - mode.w0) * (hi - lo) / (mode.w1 - mode.w0);
    for (int i = 0; i < 12; ++i) {
      const double w = phase_deriv(seg, t, mode.n, mode.m);
      const double dw = phase_second_deriv(seg, t, mode.n, mode.m);
      const double residual = w - omega;
      if (std::abs(residual) < params.spa_frequency_rtol * std::abs(omega)) {
        break;
      }
      if (residual > 0.0) {
        hi = t;
      } else {
        lo = t;
      }
      if (dw > 0.0 && std::isfinite(dw)) {
        const double candidate = t - residual / dw;
        t = (candidate > lo && candidate < hi) ? candidate : 0.5 * (lo + hi);
      } else {
        t = 0.5 * (lo + hi);
      }
    }
    const double dd = phase_second_deriv(seg, t, mode.n, mode.m);
    if (!(dd > 0.0) || !std::isfinite(dd)) {
      return false;
    }
    t_spa = t;
    phase = phase_value(seg, t, mode.n, mode.m);
    T_spa = 1.0 / std::sqrt(dd);
    return true;
  }

  bool stationary_time_for_mode_with_guess(const ModeData& mode, double omega, double& t_guess,
                                           double& phase, double& T_spa) const {
    // Production solve for sorted/uniform grids. t_guess is updated in-place so
    // the next nearby frequency starts almost at the solution.
    if (omega < mode.w0 || omega > mode.w1) {
      return false;
    }
    const Segment& seg = segments[static_cast<std::size_t>(mode.segment)];
    double lo = seg.t0;
    double hi = seg.t1;
    double t = initial_time_guess(seg, mode, omega, t_guess);
    const int max_iterations =
        use_circular_fast_path() ? std::max(1, params.circular_spa_iterations) : 8;
    for (int i = 0; i < max_iterations; ++i) {
      const double w = phase_deriv(seg, t, mode.n, mode.m);
      const double dw = phase_second_deriv(seg, t, mode.n, mode.m);
      const double residual = w - omega;
      if (std::abs(residual) < params.spa_frequency_rtol * std::abs(omega)) {
        break;
      }
      if (residual > 0.0) {
        hi = t;
      } else {
        lo = t;
      }
      // Newton gives fast convergence when the local curvature is sane; the
      // bisection fallback keeps the iterate inside the support interval.
      if (dw > 0.0 && std::isfinite(dw)) {
        const double candidate = t - residual / dw;
        t = (candidate > lo && candidate < hi) ? candidate : 0.5 * (lo + hi);
      } else {
        t = 0.5 * (lo + hi);
      }
    }
    const double dd = phase_second_deriv(seg, t, mode.n, mode.m);
    if (!(dd > 0.0) || !std::isfinite(dd)) {
      return false;
    }
    t_guess = t;
    phase = phase_value(seg, t, mode.n, mode.m);
    T_spa = 1.0 / std::sqrt(dd);
    return true;
  }

  bool is_sorted_frequency_grid(const std::vector<double>& frequencies_hz) const {
    // Sorted grids let us traverse by mode support intervals rather than trying
    // every mode at every frequency.
    for (std::size_t i = 1; i < frequencies_hz.size(); ++i) {
      if (frequencies_hz[i] < frequencies_hz[i - 1]) {
        return false;
      }
    }
    return true;
  }

  FrequencyWaveform generate_waveform(const std::vector<double>& frequencies_hz) const {
    FrequencyWaveform out;
    out.plus.assign(frequencies_hz.size(), C{0.0, 0.0});
    out.cross.assign(frequencies_hz.size(), C{0.0, 0.0});
    const double pref = std::sqrt(2.0 * kPi) * h0_pref;
    if (is_sorted_frequency_grid(frequencies_hz)) {
      // Mode-major traversal: find the frequency bins inside each mode's
      // support [w0, w1], then walk those bins with warm-started stationary-time
      // solves. This is much faster for likelihood-style sorted grids.
      for (const auto& mode : modes) {
        const double f0 = mode.w0 / (2.0 * kPi);
        const double f1 = mode.w1 / (2.0 * kPi);
        auto it0 = std::lower_bound(frequencies_hz.begin(), frequencies_hz.end(), f0);
        auto it1 = std::upper_bound(frequencies_hz.begin(), frequencies_hz.end(), f1);
        if (it0 == it1) {
          continue;
        }
        const std::size_t i0 = static_cast<std::size_t>(std::distance(frequencies_hz.begin(), it0));
        const std::size_t i1 = static_cast<std::size_t>(std::distance(frequencies_hz.begin(), it1));
        const Segment& seg = segments[static_cast<std::size_t>(mode.segment)];
        const double first_omega = 2.0 * kPi * frequencies_hz[i0];
        double t_guess = initial_time_guess(seg, mode, first_omega, seg.t0);
        for (std::size_t i = i0; i < i1; ++i) {
          const double omega = 2.0 * kPi * frequencies_hz[i];
          double phase = 0.0;
          double T_spa = 0.0;
          if (!stationary_time_for_mode_with_guess(mode, omega, t_guess, phase, T_spa)) {
            continue;
          }
          const double psi = omega * t_guess - phase - 0.25 * kPi;
          const C spa = T_spa * std::polar(1.0, psi);
          const auto A = SUA_amplitudes(t_guess, mode, T_spa);
          out.plus[i] += pref * std::conj(spa * A[0]);
          out.cross[i] += pref * std::conj(spa * A[1]);
        }
      }
      return out;
    }

    // Generic fallback for unsorted user grids. It is simple and robust, but it
    // cannot reuse stationary-time guesses across frequencies.
    for (std::size_t i = 0; i < frequencies_hz.size(); ++i) {
      const double omega = 2.0 * kPi * frequencies_hz[i];
      C hp{0.0, 0.0};
      C hc{0.0, 0.0};
      for (const auto& mode : modes) {
        double t_spa = 0.0;
        double phase = 0.0;
        double T_spa = 0.0;
        if (!stationary_time_for_mode(mode, omega, t_spa, phase, T_spa)) {
          continue;
        }
        const double psi = omega * t_spa - phase - 0.25 * kPi;
        const C spa = T_spa * std::polar(1.0, psi);
        const auto A = SUA_amplitudes(t_spa, mode, T_spa);
        hp += spa * A[0];
        hc += spa * A[1];
      }
      out.plus[i] = pref * std::conj(hp);
      out.cross[i] = pref * std::conj(hc);
    }
    return out;
  }

  struct BinRange {
    std::size_t begin = 0;
    std::size_t end = 0;

    bool empty() const {
      return begin >= end;
    }
  };

  BinRange uniform_bin_range_for_mode(const ModeData& mode, double f_min_hz,
                                      double delta_f_hz, std::size_t count) const {
    // Convert a mode's continuous frequency support into half-open uniform-grid
    // indices [begin, end). The tiny padding prevents roundoff from dropping an
    // endpoint bin that mathematically lies inside support.
    if (count == 0) {
      return {};
    }
    constexpr double edge_padding = 1.0e-12;
    const double inv_df = 1.0 / delta_f_hz;
    const double f0 = mode.w0 / (2.0 * kPi);
    const double f1 = mode.w1 / (2.0 * kPi);
    const long long count_ll = static_cast<long long>(count);
    const long long begin_ll =
        static_cast<long long>(std::ceil((f0 - f_min_hz) * inv_df - edge_padding));
    const long long end_ll =
        static_cast<long long>(std::floor((f1 - f_min_hz) * inv_df + edge_padding)) + 1;

    if (end_ll <= 0 || begin_ll >= count_ll) {
      return {};
    }
    return {static_cast<std::size_t>(std::max<long long>(0, begin_ll)),
            static_cast<std::size_t>(std::min<long long>(count_ll, end_ll))};
  }

  FrequencyWaveform generate_waveform_uniform(double f_min_hz, double delta_f_hz,
                                              std::size_t count) const {
    FrequencyWaveform out;
    out.plus.resize(count);
    out.cross.resize(count);
    generate_waveform_uniform(f_min_hz, delta_f_hz, count, out.plus.data(), out.cross.data());
    return out;
  }

  void accumulate_waveform_uniform_modes(double f_min_hz, double delta_f_hz, std::size_t count,
                                         C* plus, C* cross, std::size_t mode_begin,
                                         std::size_t mode_end) const {
    // Accumulate a subset of modes into caller-owned arrays. In serial this is
    // the whole waveform; in threaded mode each worker writes into a private
    // tile buffer and the caller reduces those buffers after joining.
    const double pref = std::sqrt(2.0 * kPi) * h0_pref;
    for (std::size_t imode = mode_begin; imode < mode_end; ++imode) {
      const auto& mode = modes[imode];
      const BinRange bins = uniform_bin_range_for_mode(mode, f_min_hz, delta_f_hz, count);
      if (bins.empty()) {
        continue;
      }
      const Segment& seg = segments[static_cast<std::size_t>(mode.segment)];
      const double omega0 = 2.0 * kPi * (f_min_hz + delta_f_hz * static_cast<double>(bins.begin));
      double t_guess = initial_time_guess(seg, mode, omega0, seg.t0);
      for (std::size_t i = bins.begin; i < bins.end; ++i) {
        const double omega = 2.0 * kPi * (f_min_hz + delta_f_hz * static_cast<double>(i));
        double phase = 0.0;
        double T_spa = 0.0;
        if (!stationary_time_for_mode_with_guess(mode, omega, t_guess, phase, T_spa)) {
          continue;
        }
        const double psi = omega * t_guess - phase - 0.25 * kPi;
        const C spa = T_spa * std::polar(1.0, psi);
        const auto A = SUA_amplitudes(t_guess, mode, T_spa);
        plus[i] += pref * std::conj(spa * A[0]);
        cross[i] += pref * std::conj(spa * A[1]);
      }
    }
  }

  void generate_waveform_uniform(double f_min_hz, double delta_f_hz, std::size_t count,
                                 C* plus, C* cross) const {
    if (!(delta_f_hz > 0.0)) {
      throw std::invalid_argument("delta_f_hz must be positive");
    }
    if ((count > 0) && (plus == nullptr || cross == nullptr)) {
      throw std::invalid_argument("output pointers must not be null");
    }
    if (count == 0) {
      return;
    }
    std::fill(plus, plus + count, C{0.0, 0.0});
    std::fill(cross, cross + count, C{0.0, 0.0});
    const int requested_threads = std::max(1, params.num_threads);
    const std::size_t nthreads =
        std::min<std::size_t>(static_cast<std::size_t>(requested_threads), modes.size());
    if (nthreads <= 1 || count == 0 || modes.size() < 2) {
      // Serial path is also the right choice for circular waveforms: there are
      // too few modes for threading to pay for itself.
      accumulate_waveform_uniform_modes(f_min_hz, delta_f_hz, count, plus, cross, 0,
                                        modes.size());
      return;
    }

    constexpr std::size_t tile_size = 32768;
    // Do not allocate a full waveform per thread. Tiling keeps memory bounded
    // while still letting threads split the expensive mode work.
    std::vector<FrequencyWaveform> partials(nthreads);
    for (auto& partial : partials) {
      partial.plus.resize(tile_size);
      partial.cross.resize(tile_size);
    }
    for (std::size_t tile0 = 0; tile0 < count; tile0 += tile_size) {
      const std::size_t tile_count = std::min(tile_size, count - tile0);
      const double tile_f_min = f_min_hz + delta_f_hz * static_cast<double>(tile0);
      std::vector<std::thread> workers;
      workers.reserve(nthreads);
      for (std::size_t ithread = 0; ithread < nthreads; ++ithread) {
        // Split by mode, not by frequency. Each mode has compact frequency
        // support, so the worker can skip bins outside its modes cheaply.
        const std::size_t begin = ithread * modes.size() / nthreads;
        const std::size_t end = (ithread + 1) * modes.size() / nthreads;
        std::fill(partials[ithread].plus.begin(),
                  partials[ithread].plus.begin() + static_cast<std::ptrdiff_t>(tile_count),
                  C{0.0, 0.0});
        std::fill(partials[ithread].cross.begin(),
                  partials[ithread].cross.begin() + static_cast<std::ptrdiff_t>(tile_count),
                  C{0.0, 0.0});
        workers.emplace_back([&, ithread, begin, end, tile_f_min, tile_count] {
          accumulate_waveform_uniform_modes(tile_f_min, delta_f_hz, tile_count,
                                            partials[ithread].plus.data(),
                                            partials[ithread].cross.data(), begin, end);
        });
      }
      for (auto& worker : workers) {
        worker.join();
      }
      for (const auto& partial : partials) {
        // The physical waveform is a linear sum over modes, so reduction is just
        // complex addition of the per-thread tile buffers.
        for (std::size_t i = 0; i < tile_count; ++i) {
          plus[tile0 + i] += partial.plus[i];
          cross[tile0 + i] += partial.cross[i];
        }
      }
    }
  }

  TimeWaveform generate_time_domain_waveform(const std::vector<double>& times_s) const {
    // Time-domain generation is mainly a diagnostic/convenience path. It sums
    // the real oscillatory modes directly on the segment that contains each
    // requested time, rather than using the stationary-phase approximation.
    TimeWaveform out;
    out.plus.assign(times_s.size(), 0.0);
    out.cross.assign(times_s.size(), 0.0);
    for (std::size_t i = 0; i < times_s.size(); ++i) {
      const double t = times_s[i];
      if (t < start_time() || t > end_time()) {
        continue;
      }
      const std::size_t iseg = find_segment(t);
      const Segment& seg = segments[iseg];
      double hp = 0.0;
      double hc = 0.0;
      for (int midx : modes_by_segment[iseg]) {
        const auto& mode = modes[static_cast<std::size_t>(midx)];
        const auto A = compute_amplitudes(t, mode);
        const double phi = phase_value(seg, t, mode.n, mode.m);
        // Add twice the real part of A exp(-i phi), matching the positive and
        // negative frequency conjugate contributions.
        const double c = std::cos(phi);
        const double s = std::sin(phi);
        hp += 2.0 * (A[0].real() * c + A[0].imag() * s);
        hc += 2.0 * (A[1].real() * c + A[1].imag() * s);
      }
      out.plus[i] = h0_pref * hp;
      out.cross[i] = h0_pref * hc;
    }
    return out;
  }

  double start_time() const { return segments.front().t0; }
  double end_time() const { return physical_end_time; }

  Parameters params;

  // Internal masses are in seconds and ordered with m1 >= m2. q1/q2 and spins
  // follow the same swap so the physical system is unchanged.
  double m1 = 0.0;
  double m2 = 0.0;
  double q1 = 1.0;
  double q2 = 1.0;

  // Initial and terminal ODE state information.
  State v_ini{};
  double yf = 0.0;

  // Spin/precession invariants reused by derivatives and Euler-angle recovery.
  double sz_1 = 0.0;
  double sz_2 = 0.0;
  double chi_eff = 0.0;
  double sp2_1 = 0.0;
  double sp2_2 = 0.0;
  double s2_1 = 0.0;
  double s2_2 = 0.0;
  double cos_theta_JN = 1.0;
  double phi_JN = 0.0;

  // Overall strain scale and precomputed line-of-sight projection factors.
  double h0_pref = 0.0;
  std::unique_ptr<PNDerivatives> pn;
  std::array<std::array<C, 2>, 5> Apc_proj{};

  // Dense ODE support and segment-local mode support. modes_by_segment is used
  // by the time-domain path; frequency-domain paths mostly traverse modes.
  std::vector<Segment> segments;
  double physical_end_time = 0.0;
  std::vector<ModeData> modes;
  std::vector<std::vector<int>> modes_by_segment;

  // SUA and interpolation caches.
  std::vector<C> ak_SUA;
  std::vector<double> interp_nodes;
  std::vector<double> interp_weights;
  std::vector<ApcSegmentCache> apc_cache;
};

Model::Model(Parameters parameters) : impl_(std::make_unique<Impl>(std::move(parameters))) {}

Model::~Model() = default;

Model::Model(Model&&) noexcept = default;

Model& Model::operator=(Model&&) noexcept = default;

FrequencyWaveform Model::generate_waveform(const std::vector<double>& frequencies_hz) const {
  return impl_->generate_waveform(frequencies_hz);
}

FrequencyWaveform Model::generate_waveform_uniform(double f_min_hz, double delta_f_hz,
                                                   std::size_t count) const {
  return impl_->generate_waveform_uniform(f_min_hz, delta_f_hz, count);
}

void Model::generate_waveform_uniform(double f_min_hz, double delta_f_hz, std::size_t count,
                                      std::complex<double>* plus,
                                      std::complex<double>* cross) const {
  impl_->generate_waveform_uniform(f_min_hz, delta_f_hz, count, plus, cross);
}

TimeWaveform Model::generate_time_domain_waveform(const std::vector<double>& times_s) const {
  return impl_->generate_time_domain_waveform(times_s);
}

double Model::start_time() const {
  return impl_->start_time();
}

double Model::end_time() const {
  return impl_->end_time();
}

std::size_t Model::segment_count() const {
  return impl_->segments.size();
}

std::size_t Model::mode_count() const {
  return impl_->modes.size();
}

const Parameters& Model::parameters() const {
  return impl_->params;
}

} // namespace pyefpe
