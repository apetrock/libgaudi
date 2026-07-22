#ifndef __GAUDI_TEST_MEDIAL_BENCH_TESTS_HPP__
#define __GAUDI_TEST_MEDIAL_BENCH_TESTS_HPP__

// Timing/robustness bench for the two medial-axis backends.
// Easy to tear down: delete this file and the single include line in
// tests/gaudi_tests.cpp.
//
// Method: seed a random Darboux-fit generator from 12 real bunny fits, jitter
// each coefficient by a fraction of its observed magnitude, and feed the
// jittered Qs to both backends with a mix of inward / outward / random mesh
// normals. Reports per-call timing, acceptance, and "convex interior" validity
// (D(center) < 0 and largest-|eigenvalue| of W is positive — the medial-ridge
// signature). Jittered Qs naturally span "seed in trust region" (D(0)<0) and
// "seed outside / degenerate" (D(0)>=0), so this exercises both the happy path
// and the trust-region return / rejection logic.

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

#include "gaudi/albers/darboux_cyclide.hpp"
#include "gaudi/albers/darboux_medial_geometry.hpp"
#include "gaudi/albers/medial_shape_search.hpp"
#include "gaudi/duchamp/medial_search_helpers.hpp"
#include "gaudi/test/test.hpp"

namespace gaudi {
namespace test {

namespace medial_bench_detail {

// 12 real bunny darboux fits (probe seeds): A,B,C,D,E,F,G,H,I,J,lambda,mu,nu,
// kappa. Used to seed the random-fit generator so jittered samples stay in a
// realistic region of Q-space rather than uniform noise.
inline const std::vector<albers::vec14> &seed_fits() {
  static const std::vector<albers::vec14> seeds = {
    albers::vec14(0.442477, 0.341304, 0.932373, -0.077291, 0.131518, 0.177204, -0.154419, -0.388669, 0.184669, -0.049614, -0.210913, 0.247227, 0.526376, -0.255675),
    albers::vec14(0.410120, 0.394230, 0.951950, -0.125180, 0.087223, 0.147610, -0.248520, -0.310760, 0.097911, -0.089221, -0.203810, 0.354020, 0.455170, -0.099132),
    albers::vec14(0.420940, 0.322640, 0.950180, -0.126830, 0.094400, 0.129160, -0.187220, -0.386180, 0.104270, -0.063973, -0.225350, 0.313080, 0.569100, -0.139940),
    albers::vec14(0.443610, 0.314100, 0.958650, -0.143890, 0.072473, 0.101640, -0.204270, -0.379200, 0.056541, -0.076218, -0.231530, 0.332440, 0.586090, -0.055712),
    albers::vec14(0.429560, 0.333180, 0.958350, -0.131680, 0.088652, 0.125780, -0.201600, -0.376310, 0.090129, -0.074587, -0.226270, 0.325480, 0.561300, -0.113410),
    albers::vec14(1.507800, 0.966500, 0.850530, -0.207190, 0.113820, 0.496200, -0.111810, 0.316080, -0.064034, -0.045334, 0.112930, -0.706960, 0.332950, 0.382360),
    albers::vec14(0.418720, 0.389060, 0.940320, -0.083723, 0.116490, 0.165270, -0.213420, -0.336630, 0.164070, -0.074723, -0.199500, 0.314770, 0.476280, -0.132850),
    albers::vec14(1.018700, 0.956790, 0.572690, -0.318900, 0.162460, 0.705720, 0.081715, 0.137790, -0.113220, -0.057185, 0.198140, -0.770430, 0.611190, 0.541770),
    albers::vec14(0.956450, 0.947310, 0.667550, -0.278900, 0.152040, 0.690250, 0.053077, 0.176700, -0.162210, -0.051104, 0.121460, -0.607550, 0.458430, 0.593020),
    albers::vec14(0.434440, -1.044300, 2.260700, -0.117590, 0.446830, -0.708740, 0.014087, 0.446330, 0.112450, -0.012427, 0.393670, 0.466000, -0.319460, 1.809000),
    albers::vec14(0.449050, 0.353830, 0.944160, -0.069251, 0.095648, 0.110940, -0.168090, -0.394550, 0.121030, -0.061465, -0.179580, 0.213140, 0.468170, -0.147270),
    albers::vec14(0.446360, 0.383120, 0.928250, -0.059742, 0.106620, 0.121740, -0.180050, -0.384870, 0.151070, -0.055336, -0.174260, 0.216350, 0.436720, -0.171530),
  };
  return seeds;
}

// Per-coefficient jitter scale = median |value| across the seed set, floored
// to a small epsilon so near-zero coefficients still get nonzero jitter.
inline albers::vec14 jitter_scale() {
  const auto &s = seed_fits();
  albers::vec14 sc = albers::vec14::Zero();
  for (int k = 0; k < 14; ++k) {
    std::vector<real> col;
    col.reserve(s.size());
    for (const auto &q : s) {
      col.push_back(std::abs(q[k]));
    }
    std::sort(col.begin(), col.end());
    real med = col[col.size() / 2];
    sc[k] = std::max(med, real(1e-3));
  }
  return sc;
}

// Generate N jittered Qs. Each is a random seed fit + Gaussian noise scaled by
// jitter_frac * per-coefficient median magnitude. jitter_frac ~ 0.3 keeps most
// samples geometrically plausible; larger values push into degenerate regime.
inline std::vector<albers::vec14> random_fits(size_t n, uint32_t seed,
                                              real jitter_frac) {
  std::mt19937 rng(seed);
  std::normal_distribution<real> nd(real(0.0), real(1.0));
  std::uniform_int_distribution<size_t> pick(0, seed_fits().size() - 1);
  const albers::vec14 sc = jitter_scale();
  std::vector<albers::vec14> out;
  out.reserve(n);
  for (size_t i = 0; i < n; ++i) {
    albers::vec14 Q = seed_fits()[pick(rng)];
    for (int k = 0; k < 14; ++k) {
      Q[k] += jitter_frac * sc[k] * nd(rng);
    }
    out.push_back(Q);
  }
  return out;
}

// Random mesh normals: 1/3 inward (-g0, good seed direction), 1/3 outward
// (+g0, bad — points away from interior), 1/3 random. This mixes "seed points
// the search can use" with "seed points outside/degenerate" so the bench
// exercises both the happy path and the trust-region return / rejection.
inline std::vector<vec3> random_normals(const std::vector<albers::vec14> &Qs,
                                        uint32_t seed) {
  std::mt19937 rng(seed);
  std::normal_distribution<real> nd(real(0.0), real(1.0));
  std::uniform_int_distribution<int> mode(0, 2);
  std::vector<vec3> out;
  out.reserve(Qs.size());
  for (const auto &Q : Qs) {
    vec3 g0 = albers::darboux_grad(Q, vec3::Zero());
    if (!g0.allFinite() || g0.norm() < real(1e-12)) {
      g0 = vec3(real(0.0), real(0.0), real(1.0));
    }
    g0.normalize();
    vec3 n;
    const int m = mode(rng);
    if (m == 0) {
      n = -g0; // inward
    } else if (m == 1) {
      n = g0; // outward
    } else {
      vec3 r(nd(rng), nd(rng), nd(rng));
      n = r.normalized(); // random
    }
    out.push_back(n);
  }
  return out;
}

struct bench_stats {
  size_t n = 0;
  int accepted = 0;
  int convex_interior = 0; // accepted AND D(center)<0 AND k1>0 (medial-ridge signature)
  int seed_in_trust = 0;   // D(0) < 0 (foot inside the trust region)
  double total_us = 0.0;
  double mean_travel = 0.0;
  real max_travel = 0.0;
};

// Classify a medial result: convex-interior iff D(center)<0 and the
// largest-magnitude eigenvalue of W is positive (ball-like / convex level set,
// as opposed to the concave far-field artifact).
inline bool convex_interior_ok(const albers::vec14 &Q, const vec3 &center_local) {
  const real D = albers::eval_darboux(Q, center_local);
  if (!(D < real(0.0))) {
    return false;
  }
  const mat3 W = albers::shape_operator_at(Q, center_local);
  if (!W.allFinite()) {
    return false;
  }
  Eigen::SelfAdjointEigenSolver<mat3> es(W);
  if (es.info() != Eigen::Success) {
    return false;
  }
  std::array<real, 3> ev = {es.eigenvalues()[0], es.eigenvalues()[1],
                            es.eigenvalues()[2]};
  std::sort(ev.begin(), ev.end(),
            [](real a, real b) { return std::abs(a) > std::abs(b); });
  return ev[0] > real(0.0);
}

inline bench_stats run_legacy(const std::vector<albers::vec14> &Qs,
                              const std::vector<vec3> &Ns, real max_travel,
                              int max_iters, real tol, real max_newton_step) {
  bench_stats s;
  s.n = Qs.size();
  real travel_sum = real(0.0);
  int travel_n = 0;
  auto t0 = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < Qs.size(); ++i) {
    if (albers::eval_darboux(Qs[i], vec3::Zero()) < real(0.0)) {
      ++s.seed_in_trust;
    }
    const auto r = duchamp::search_medial_legacy_ridge(
        Qs[i], vec3::Zero(), Ns[i], max_travel, max_iters, tol,
        max_newton_step, real(-2.0) // lax normal alignment — bench both paths
    );
    if (r.accepted) {
      ++s.accepted;
      const vec3 xl = r.point; // foot = origin => center_local == center_world
      if (convex_interior_ok(Qs[i], xl)) {
        ++s.convex_interior;
      }
      const real travel = xl.norm();
      travel_sum += travel;
      ++travel_n;
      s.max_travel = std::max(s.max_travel, travel);
    }
  }
  auto t1 = std::chrono::high_resolution_clock::now();
  s.total_us = std::chrono::duration<double, std::micro>(t1 - t0).count();
  s.mean_travel = travel_n ? double(travel_sum) / double(travel_n) : 0.0;
  return s;
}

inline bench_stats run_shape(const std::vector<albers::vec14> &Qs,
                             const std::vector<vec3> &Ns, real max_travel,
                             const albers::medial_shape_search_params &sp) {
  bench_stats s;
  s.n = Qs.size();
  real travel_sum = real(0.0);
  int travel_n = 0;
  auto t0 = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < Qs.size(); ++i) {
    if (albers::eval_darboux(Qs[i], vec3::Zero()) < real(0.0)) {
      ++s.seed_in_trust;
    }
    const auto r =
        duchamp::search_medial_shape_energy(Qs[i], vec3::Zero(), Ns[i],
                                            max_travel, sp);
    if (r.accepted) {
      ++s.accepted;
      const vec3 xl = r.point;
      if (convex_interior_ok(Qs[i], xl)) {
        ++s.convex_interior;
      }
      const real travel = xl.norm();
      travel_sum += travel;
      ++travel_n;
      s.max_travel = std::max(s.max_travel, travel);
    }
  }
  auto t1 = std::chrono::high_resolution_clock::now();
  s.total_us = std::chrono::duration<double, std::micro>(t1 - t0).count();
  s.mean_travel = travel_n ? double(travel_sum) / double(travel_n) : 0.0;
  return s;
}

inline void report(const char *name, const bench_stats &s) {
  const double per = s.n ? s.total_us / double(s.n) : 0.0;
  const int pct_acc = s.n ? int(100.0 * double(s.accepted) / double(s.n)) : 0;
  const int pct_conv = s.n ? int(100.0 * double(s.convex_interior) / double(s.n)) : 0;
  std::cerr << "  " << std::left << std::setw(22) << name << " n=" << s.n
            << " seed_in_trust=" << s.seed_in_trust
            << " accepted=" << s.accepted << " (" << pct_acc << "%)"
            << " convex_interior=" << s.convex_interior << " (" << pct_conv
            << "%)"
            << " mean_travel=" << std::fixed << std::setprecision(3)
            << s.mean_travel << " max_travel=" << s.max_travel
            << " | total=" << std::setprecision(1) << s.total_us << "us"
            << " per=" << per << "us\n";
}

} // namespace medial_bench_detail

GAUDI_TEST(medial_backend_bench) {
  namespace d = medial_bench_detail;
  const real max_travel = 5.0;     // ~bunny medial scale
  const int max_iters = 100;
  const real tol = 1e-8;
  const real max_newton_step = 5.0;
  const size_t N = 4000;
  const uint32_t seed = 0xC0FFEEu;
  const real jitter = 0.3;

  auto Qs = d::random_fits(N, seed, jitter);
  auto Ns = d::random_normals(Qs, seed + 1u);
  albers::medial_shape_search_params sp;
  sp.max_iters = 24;

  std::cerr << "medial backend bench: N=" << N << " jitter=" << jitter
            << " max_travel=" << max_travel
            << " (mix of inward/outward/random normals)\n";
  auto legacy = d::run_legacy(Qs, Ns, max_travel, max_iters, tol, max_newton_step);
  auto shape = d::run_shape(Qs, Ns, max_travel, sp);
  d::report("LegacyRidge", legacy);
  d::report("ShapeEnergy(hybrid)", shape);
  GAUDI_ASSERT(legacy.n == shape.n);
}

// Ridge smoke test: verify an accepted medial point is actually a local
// minimum of E = 1/(eps + ||W||^2) along the search ray, by probing E at
// t +/- delta. A true medial ridge of the (interior-truncated) energy must
// satisfy E(t) <= E(t-delta) and E(t) <= E(t+delta). Also checks D(center)<0
// (interior) and |E'(t)| small. Run on the real seed fits with inward normals
// so the trust region is clean.
GAUDI_TEST(medial_ridge_local_min_smoke) {
  namespace d = medial_bench_detail;
  const real max_travel = 5.0;
  albers::medial_shape_search_params sp;
  sp.max_iters = 24;

  int checked = 0;
  int local_min_ok = 0;
  int interior_ok = 0;
  int deriv_ok = 0;
  for (const auto &Q : d::seed_fits()) {
    vec3 g0 = albers::darboux_grad(Q, vec3::Zero());
    if (!g0.allFinite() || g0.norm() < real(1e-12)) {
      continue;
    }
    // Outward surface normal (g0 points outward for these convex fits); the
    // search's aligned_inward_ray_dir expects the outward normal so it can flip
    // the march to be opposite it (inward). Passing -g0 here would invert that
    // flip and march outward.
    const vec3 N = g0.normalized();
    const auto r =
        duchamp::search_medial_shape_energy(Q, vec3::Zero(), N, max_travel, sp);
    const real travel = r.point.norm();
    if (!r.accepted || travel <= real(1e-6)) {
      continue;
    }
    ++checked;

    // Reconstruct the ray and the energy at t, t-delta, t+delta.
    const vec3 dir = albers::aligned_inward_ray_dir(Q, N);
    if (dir.squaredNorm() < real(1e-24)) {
      continue;
    }
    const albers::darboux_geometry_bundle geom(Q);
    const albers::ray_line_bundle line(vec3::Zero(), dir.normalized());
    // t along the ray = center_local . dir_hat (foot at origin).
    const real t = r.point.dot(dir.normalized());
    const real delta = std::max(real(1e-3), real(0.05) * travel);

    auto E_at = [&](real tt) -> real {
      real E = 0.0, Ep = 0.0, Epp = 0.0;
      albers::medial_generated::eval_medial_energy_at_t(geom, line, tt, E, Ep,
                                                         Epp, real(1e-12));
      return E;
    };
    const real E_m = E_at(t - delta);
    const real E_0 = E_at(t);
    const real E_p = E_at(t + delta);

    if (E_0 <= E_m + real(1e-9) * (std::abs(E_0) + std::abs(E_m) + real(1.0)) &&
        E_0 <= E_p + real(1e-9) * (std::abs(E_0) + std::abs(E_p) + real(1.0))) {
      ++local_min_ok;
    }
    if (albers::eval_darboux(Q, r.point) < real(0.0)) {
      ++interior_ok;
    }
    // |E'(t)| should be small at the accepted minimum.
    real E = 0.0, Ep = 0.0, Epp = 0.0;
    albers::medial_generated::eval_medial_energy_at_t(geom, line, t, E, Ep, Epp,
                                                       real(1e-12));
    if (std::abs(Ep) < real(1e-3) * std::max(real(1.0), std::abs(E_0))) {
      ++deriv_ok;
    }
  }

  std::cerr.copyfmt(std::ios(nullptr)); // reset any leaked std::fixed/precision
  std::cerr << "medial ridge local-min smoke: checked=" << checked
            << " local_min=" << local_min_ok << " interior=" << interior_ok
            << " deriv_small=" << deriv_ok << "\n";
  // We must have checked at least a handful of the seed fits, and every
  // checked accepted point must be a local min of E and lie in the interior.
  GAUDI_ASSERT(checked >= 6);
  GAUDI_ASSERT(local_min_ok == checked);
  GAUDI_ASSERT(interior_ok == checked);
}

} // namespace test
} // namespace gaudi

#endif // __GAUDI_TEST_MEDIAL_BENCH_TESTS_HPP__
