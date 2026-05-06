#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "pyefpe/pyefpe.hpp"

namespace {

void check_model(const pyefpe::Parameters& p) {
  pyefpe::Model waveform(p);

  if (!(waveform.start_time() < waveform.end_time())) {
    throw std::runtime_error("invalid waveform time range");
  }
  if (waveform.segment_count() == 0 || waveform.mode_count() == 0) {
    throw std::runtime_error("model did not create integration segments or modes");
  }

  std::vector<double> freqs;
  for (double f = 30.0; f < 192.0; f += 2.0) {
    freqs.push_back(f);
  }
  const auto h = waveform.generate_waveform(freqs);
  if (h.plus.size() != freqs.size() || h.cross.size() != freqs.size()) {
    throw std::runtime_error("frequency-domain waveform has the wrong size");
  }

  double norm = 0.0;
  for (std::size_t i = 0; i < freqs.size(); ++i) {
    if (!std::isfinite(h.plus[i].real()) || !std::isfinite(h.plus[i].imag()) ||
        !std::isfinite(h.cross[i].real()) || !std::isfinite(h.cross[i].imag())) {
      throw std::runtime_error("frequency-domain waveform contains non-finite values");
    }
    norm += std::norm(h.plus[i]) + std::norm(h.cross[i]);
  }
  if (!(norm > 0.0)) {
    throw std::runtime_error("frequency-domain waveform is identically zero");
  }

  std::vector<double> times;
  const double t0 = waveform.start_time();
  const double t1 = waveform.end_time();
  for (int i = 0; i < 1024; ++i) {
    times.push_back(t0 + (t1 - t0) * static_cast<double>(i) / 1023.0);
  }
  const auto ht = waveform.generate_time_domain_waveform(times);
  if (ht.plus.size() != times.size() || ht.cross.size() != times.size()) {
    throw std::runtime_error("time-domain waveform has the wrong size");
  }
}

void check_helpers() {
  const double f_isco = pyefpe::f22_isco_frequency_hz(1000.0, 10000.0);
  if (!(f_isco > 0.01)) {
    throw std::runtime_error("unexpected ISCO frequency");
  }
  const double f_start =
      pyefpe::estimate_f22_start_for_duration_newtonian(1000.0, 10000.0, 0.01, 31557600.0);
  if (!(f_start >= 1.0e-5 && f_start < 0.01)) {
    throw std::runtime_error("invalid Newtonian start-frequency estimate");
  }

  pyefpe::FrequencyWaveform h;
  h.plus = {{1.0, 2.0}, {3.0, 4.0}};
  h.cross = {{-2.0, 1.0}, {0.5, -0.25}};
  pyefpe::ExtrinsicTransform transform;
  transform.amplitude_scale = 2.0;
  pyefpe::apply_extrinsic_transform_uniform(0.0, 1.0, transform, h);
  if (std::abs(h.plus[0] - std::complex<double>{2.0, 4.0}) > 1.0e-14 ||
      std::abs(h.cross[1] - std::complex<double>{1.0, -0.5}) > 1.0e-14) {
    throw std::runtime_error("extrinsic amplitude transform failed");
  }
}

} // namespace

int main() {
  check_helpers();

  pyefpe::Parameters p;
  p.mass1 = 8.0;
  p.mass2 = 6.0;
  p.e_start = 0.2;
  p.spin1z = 0.1;
  p.spin2z = -0.05;
  p.f22_start = 20.0;
  p.f22_end = 256.0;
  p.amplitude_tol = 1.0e-2;

  check_model(p);

  p.dense_output = pyefpe::DenseOutput::DormandPrince;
  check_model(p);

  std::cout << "ok\n";
}
