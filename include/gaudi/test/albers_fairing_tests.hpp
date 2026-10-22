#ifndef GAUDI_ALBERS_FAIRING_TESTS_HPP
#define GAUDI_ALBERS_FAIRING_TESTS_HPP

#include "gaudi/hepworth/block/albers_fairing.hpp"
#include "gaudi/test/test.hpp"

#include <cmath>
#include <vector>

namespace gaudi {
namespace test {
namespace {

using hepworth::block::albers_fit_sample;
using hepworth::block::albers_quadric_value;
using hepworth::block::fit_albers_cyclide;
using hepworth::block::fit_albers_quadric;
using hepworth::block::fit_albers_sphere;
using hepworth::block::land_albers_cyclide;
using hepworth::block::land_albers_quadric;
using hepworth::block::land_albers_sphere;
using hepworth::block::land_albers_torus;

vec3 sphere_point(real theta, real phi) {
  const real st = std::sin(theta);
  return vec3(st * std::cos(phi), st * std::sin(phi), std::cos(theta));
}

albers_fit_sample on_sphere(const vec3 &foot, const vec3 &p) {
  albers_fit_sample s;
  s.dp = p - foot;
  s.N = p.normalized();
  s.w = 1;
  return s;
}

std::vector<albers_fit_sample> sphere_cloud(const vec3 &foot) {
  std::vector<albers_fit_sample> samples;
  for (int i = 0; i < 6; ++i) {
    const real theta = real(0.35) + real(2.4) * real(i) / real(5);
    for (int j = 0; j < 8; ++j) {
      const real phi = real(2) * real(M_PI) * real(j) / real(8);
      const vec3 p = sphere_point(theta, phi);
      if ((p - foot).norm() < real(1e-3))
        continue;
      samples.push_back(on_sphere(foot, p));
    }
  }
  return samples;
}

} // namespace

/// Gradient walk from outside the unit sphere lands on it.
GAUDI_TEST(albers_quadric_line_hits_the_sphere) {
  albers::vec10 Q = albers::vec10::Zero();
  Q[0] = Q[1] = Q[2] = 1;
  Q[9] = -1;
  const vec3 d(0, 0, 2);
  vec3 hit;
  GAUDI_EXPECT(land_albers_quadric(Q, d, hit));
  GAUDI_EXPECT((hit - vec3(0, 0, 1)).norm() < real(1e-8));
  GAUDI_EXPECT(std::abs(albers_quadric_value(Q, hit)) < real(1e-8));
}

/// Calder's constrained sphere, in the foot frame. A sample already on the
/// sphere does not move.
GAUDI_TEST(albers_sphere_fit_recovers_the_osculating_sphere) {
  const vec3 foot = vec3::UnitZ();
  const std::vector<albers_fit_sample> samples = sphere_cloud(foot);
  vec3 center;
  real radius = 0;
  GAUDI_EXPECT(fit_albers_sphere(samples, real(0.5), vec3::UnitZ(), center, radius));
  GAUDI_EXPECT((center - vec3(0, 0, -1)).norm() < real(0.05));
  GAUDI_EXPECT(std::abs(radius - real(1)) < real(0.05));

  const albers_fit_sample spoke = on_sphere(foot, sphere_point(0.6, 0.4));
  vec3 hit;
  GAUDI_EXPECT(land_albers_sphere(center, radius, spoke.dp, hit));
  GAUDI_EXPECT((hit - spoke.dp).norm() < real(0.05));
}

/// Same cloud, through the 10-coefficient quadric. The spoke stays put.
GAUDI_TEST(albers_quadric_fit_holds_a_sphere_sample) {
  const vec3 foot = vec3::UnitZ();
  const std::vector<albers_fit_sample> samples = sphere_cloud(foot);
  albers::vec10 Q;
  GAUDI_EXPECT(fit_albers_quadric(samples, real(0.5), vec3::UnitZ(), Q));
  GAUDI_EXPECT(std::abs(albers_quadric_value(Q, vec3::Zero())) < real(1e-6));
  const albers_fit_sample spoke = on_sphere(foot, sphere_point(0.8, 1.1));
  GAUDI_EXPECT(std::abs(albers_quadric_value(Q, spoke.dp)) < real(1e-3));
  vec3 hit;
  GAUDI_EXPECT(land_albers_quadric(Q, spoke.dp, hit));
  GAUDI_EXPECT((hit - spoke.dp).norm() < real(1e-2));
}

/// Closest point on the tube, not a ray. A point on the tube stays.
GAUDI_TEST(albers_torus_closest_point_is_a_fixed_point) {
  albers::osculating_torus tor;
  tor.center = vec3(0, 0, -2);
  tor.axis = vec3::UnitX();
  tor.R = 1;
  tor.r = 1;
  tor.sign = 1;
  tor.valid = true;
  const vec3 on(0, 2, -2);
  GAUDI_EXPECT(std::abs(albers::torus_sdf(on, tor)) < real(1e-8));
  vec3 hit;
  GAUDI_EXPECT(land_albers_torus(tor, on, hit));
  GAUDI_EXPECT((hit - on).norm() < real(1e-6));

  const vec3 off = on + vec3(0, real(0.1), 0);
  GAUDI_EXPECT(land_albers_torus(tor, off, hit));
  GAUDI_EXPECT((hit - on).norm() < real(1e-6));
}

/// Normal-constrained cyclide of a torus. A sample on that torus stays on it.
GAUDI_TEST(albers_cyclide_fit_holds_a_torus_sample) {
  const real R = real(1.25);
  const real r = real(0.35);
  auto point = [&](real u, real v) {
    return vec3((R + r * std::cos(v)) * std::cos(u),
                (R + r * std::cos(v)) * std::sin(u), r * std::sin(v));
  };
  auto normal = [&](real u, real v) {
    const vec3 radial(std::cos(u), std::sin(u), 0);
    return (point(u, v) - R * radial).normalized();
  };
  const vec3 foot = point(0.3, 0.4);
  std::vector<albers_fit_sample> samples;
  for (int i = 0; i < 12; ++i) {
    const real u = real(2) * real(M_PI) * real(i) / real(12);
    for (int j = 0; j < 8; ++j) {
      const real v = real(2) * real(M_PI) * real(j) / real(8);
      const vec3 p = point(u, v);
      if ((p - foot).norm() < real(1e-3))
        continue;
      albers_fit_sample s;
      s.dp = p - foot;
      s.N = normal(u, v);
      s.w = 1;
      samples.push_back(s);
    }
  }
  albers::vec14 Q;
  GAUDI_EXPECT(fit_albers_cyclide(samples, r, normal(0.3, 0.4), Q));
  const vec3 spoke = point(1.0, 1.2) - foot;
  GAUDI_EXPECT(std::abs(albers::eval_darboux(Q, spoke)) < real(1e-2));
  vec3 hit;
  GAUDI_EXPECT(land_albers_cyclide(Q, spoke, hit));
  GAUDI_EXPECT((hit - spoke).norm() < real(0.05));
}

} // namespace test
} // namespace gaudi

#endif
