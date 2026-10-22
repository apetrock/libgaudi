#ifndef GAUDI_HOPF_SWEEP_TESTS_HPP
#define GAUDI_HOPF_SWEEP_TESTS_HPP

#include "gaudi/hepworth/block/hopf_sweep.hpp"
#include "gaudi/hepworth/block/shell_constraints.hpp"
#include "gaudi/test/test.hpp"

#include <cmath>

namespace gaudi {
namespace test {
namespace {

using hepworth::block::hopf_circle_at;
using hepworth::block::hopf_make_table;
using hepworth::block::hopf_modulus;
using hepworth::block::hopf_profile;
using hepworth::block::hopf_fiber_axes;
using hepworth::block::hopf_profile_make;
using hepworth::block::hopf_table;
using hepworth::block::jacobi_sn_cn_dn;
using hepworth::block::hopf_project_fiber;
using hepworth::block::hopf_project_spin;
using hepworth::block::hopf_wave_phase;
using hepworth::block::hopf_wave_sample;
using hepworth::block::hopf_wave_slope;
using hepworth::block::k_hopf_spin_max;
using hepworth::block::hopf_elastic_kappa;
using hepworth::block::hopf_pinkall_height;
using hepworth::block::hopf_cylinder_land;
using hepworth::block::hopf_aligned_phase;
using hepworth::block::hopf_phase_shift;
using hepworth::block::hopf_direction_curvature;
using hepworth::block::hopf_normal_angle;
using hepworth::block::hopf_revolution_axes;
using hepworth::block::hopf_section_curvature;
using hepworth::block::hopf_curvature_move;
using hepworth::block::hopf_curvature_normal_offset;
using hepworth::block::hopf_curvature_plus_sin;
using hepworth::block::hopf_wave_move;
using hepworth::block::hopf_wave_scale;
using hepworth::block::hopf_section_phase;
using hepworth::block::hopf_section_sample;
using hepworth::block::hopf_ring_node;
using hepworth::block::hopf_ring_phase;
using hepworth::block::hopf_radial_land;
using hepworth::block::hopf_tube_curvature;
using hepworth::block::k_hopf_sphere_curvature;
using hepworth::block::land_hopf;

hopf_profile circle_profile(real kappa) {
  return hopf_profile_make(vec3::UnitX(), vec3::UnitY(), vec3::UnitZ(), kappa,
                           0);
}

} // namespace

GAUDI_TEST(hopf_jacobi_matches_cos_and_sech) {
  real sn = 0, cn = 0, dn = 0;
  jacobi_sn_cn_dn(0.7, 0, sn, cn, dn);
  GAUDI_EXPECT(std::abs(cn - std::cos(0.7)) < real(1e-12));
  GAUDI_EXPECT(std::abs(sn - std::sin(0.7)) < real(1e-12));

  jacobi_sn_cn_dn(0.4, 1, sn, cn, dn);
  GAUDI_EXPECT(std::abs(cn - real(1) / std::cosh(0.4)) < real(1e-12));

  const real u = 0.3;
  const real m = 0.5;
  jacobi_sn_cn_dn(u, m, sn, cn, dn);
  GAUDI_EXPECT(std::abs(sn * sn + cn * cn - real(1)) < real(1e-9));
  GAUDI_EXPECT(std::abs(dn * dn - (real(1) - m * sn * sn)) < real(1e-9));

  real sn2 = 0, cn2 = 0, dn2 = 0;
  const real du = real(1e-5);
  jacobi_sn_cn_dn(u + du, m, sn2, cn2, dn2);
  const real dcn = (cn2 - cn) / du;
  GAUDI_EXPECT(std::abs(dcn + sn * dn) < real(1e-5));
}

GAUDI_TEST(hopf_twist_zero_is_the_circle) {
  GAUDI_EXPECT(hopf_modulus(0) == real(0));
  GAUDI_EXPECT(hopf_modulus(-0.4) == real(0));
  GAUDI_EXPECT(std::abs(hopf_modulus(0.3) - real(0.3)) < real(1e-15));
  GAUDI_EXPECT(hopf_modulus(4) < real(1));

  const hopf_profile hp = circle_profile(2);
  GAUDI_EXPECT(hp.ok);
  GAUDI_EXPECT(hp.m == real(0));

  const real s = 0.05;
  const real axial = 0.4;
  real x = 0, y = 0;
  hopf_circle_at(hp.kappa0, s, x, y);
  const real series = s * s * hp.kappa0 / real(2);
  GAUDI_EXPECT(std::abs(y - series) < real(1e-4));

  const vec3 dp = axial * hp.axis + x * hp.tube + y * hp.N;
  vec3 hit = vec3::Zero();
  GAUDI_EXPECT(land_hopf(hp, dp, hit));
  GAUDI_EXPECT((hit - dp).norm() < real(1e-8));
}

GAUDI_TEST(hopf_circle_pulls_a_radial_miss_onto_the_section) {
  const hopf_profile hp = circle_profile(2);
  const real axial = 0.4;
  const vec3 dp = axial * hp.axis + real(0.2) * hp.N;
  vec3 hit = vec3::Zero();
  GAUDI_EXPECT(land_hopf(hp, dp, hit));
  GAUDI_EXPECT((hit - axial * hp.axis).norm() < real(1e-8));
  GAUDI_EXPECT(std::abs(hit.dot(hp.tube)) < real(1e-8));
}

GAUDI_TEST(hopf_negative_curvature_circle_is_a_fixed_point) {
  const hopf_profile hp = circle_profile(-2);
  const real s = 0.12;
  real x = 0, y = 0;
  hopf_circle_at(hp.kappa0, s, x, y);
  GAUDI_EXPECT(y < real(0));
  const vec3 dp = real(0.25) * hp.axis + x * hp.tube + y * hp.N;
  vec3 hit = vec3::Zero();
  GAUDI_EXPECT(land_hopf(hp, dp, hit));
  GAUDI_EXPECT((hit - dp).norm() < real(1e-8));
}

GAUDI_TEST(hopf_spin_turns_the_landing_about_the_normal) {
  const vec3 hit(0.3, 0, 0.1);
  const vec3 N = vec3::UnitZ();
  const vec3 t(0.3, 0, 0);
  const real wz = 0.4;
  const vec3 expect =
      std::cos(wz) * t + std::sin(wz) * N.cross(t) + vec3(0, 0, 0.1);
  GAUDI_EXPECT((hopf_project_spin(hit, N, wz) - expect).norm() < real(1e-8));

  const vec3 flip = hopf_project_spin(hit, N, k_hopf_spin_max + real(0.1));
  GAUDI_EXPECT((flip - hit).norm() < real(1e-12));
}

GAUDI_TEST(hopf_elliptic_sample_is_a_fixed_point) {
  const hopf_profile hp =
      hopf_profile_make(vec3::UnitX(), vec3::UnitY(), vec3::UnitZ(), 1.5, 0.4);
  GAUDI_EXPECT(hp.ok);
  GAUDI_EXPECT(hp.table.ok);
  const int i = hopf_table::N + 10;
  const vec3 dp = real(0.2) * hp.axis + hp.table.x[i] * hp.tube +
                  hp.table.y[i] * hp.N;
  vec3 hit = vec3::Zero();
  GAUDI_EXPECT(land_hopf(hp, dp, hit));
  GAUDI_EXPECT((hit - dp).norm() < real(1e-6));

  real sn = 0, cn = 0, dn = 0;
  jacobi_sn_cn_dn(0, hp.m, sn, cn, dn);
  GAUDI_EXPECT(std::abs(cn - real(1)) < real(1e-12));
}

GAUDI_TEST(hopf_wave_phase_recovers_a_shifted_crest) {
  const real kappa = 1.5;
  const real m = 0.4;
  const real phi = 0.35;
  real sn = 0, cn = 1, dn = 1;
  jacobi_sn_cn_dn(phi, m, sn, cn, dn);
  const real amp = kappa / cn;
  const real slope = -amp * std::abs(amp) * sn * dn;
  const real found = hopf_wave_phase(kappa, slope, m);
  GAUDI_EXPECT(std::abs(found - phi) < real(0.05));
  const real found_fast = hopf_wave_phase(kappa, real(2) * slope, m, 2);
  GAUDI_EXPECT(std::abs(found_fast - phi) < real(0.05));
  GAUDI_EXPECT(std::abs(hopf_wave_phase(kappa, 0, m)) < real(1e-12));

  hopf_wave_sample samples[2];
  samples[0].sigma = 0.2;
  samples[0].alpha = 0;
  samples[0].kappa = real(2) + real(0.4) * real(0.1);
  samples[1].sigma = -0.2;
  samples[1].alpha = 0;
  samples[1].kappa = real(2) - real(0.4) * real(0.1);
  GAUDI_EXPECT(std::abs(hopf_wave_slope(2, 0.5, samples, 2) - real(0.4)) <
               real(1e-8));
}

GAUDI_TEST(hopf_phased_profile_keeps_the_foot_curvature) {
  const real kappa = 1.5;
  const real phase = 0.35;
  const hopf_profile hp = hopf_profile_make(
      vec3::UnitX(), vec3::UnitY(), vec3::UnitZ(), kappa, 0.4, phase);
  GAUDI_EXPECT(hp.ok);
  GAUDI_EXPECT(std::abs(hp.phase - phase) < real(1e-12));
  real sn = 0, cn = 1, dn = 1;
  jacobi_sn_cn_dn(phase, hp.m, sn, cn, dn);
  const real foot = hp.table.kappa0 * cn;
  GAUDI_EXPECT(std::abs(foot - kappa) < real(1e-8));
}

GAUDI_TEST(hopf_fiber_rolls_the_section_about_the_ruling) {
  const vec3 axis = vec3::UnitX();
  const vec3 hit(0.4, 0.2, 0.1);
  const real wz = 0.3;
  const real c = std::cos(wz);
  const real s = std::sin(wz);
  const vec3 expect = c * hit + s * axis.cross(hit) + (real(1) - c) * hit.x() * axis;
  GAUDI_EXPECT((hopf_project_fiber(hit, axis, wz) - expect).norm() < real(1e-8));
  const vec3 along = hopf_project_fiber(vec3(0.4, 0, 0), axis, wz);
  GAUDI_EXPECT((along - vec3(0.4, 0, 0)).norm() < real(1e-8));
  const vec3 flip = hopf_project_fiber(hit, axis, k_hopf_spin_max + real(0.1));
  GAUDI_EXPECT((flip - hit).norm() < real(1e-12));
}

GAUDI_TEST(hopf_fiber_angle_rotates_the_ruling_about_the_normal) {
  const vec3 axis0 = vec3::UnitX();
  const vec3 tube0 = vec3::UnitY();
  vec3 axis = vec3::Zero();
  vec3 tube = vec3::Zero();
  hopf_fiber_axes(axis0, tube0, 0, axis, tube);
  GAUDI_EXPECT((axis - axis0).norm() < real(1e-12));
  GAUDI_EXPECT((tube - tube0).norm() < real(1e-12));
  hopf_fiber_axes(axis0, tube0, real(0.5) * M_PI, axis, tube);
  GAUDI_EXPECT((axis - tube0).norm() < real(1e-8));
  GAUDI_EXPECT((tube + axis0).norm() < real(1e-8));
}

GAUDI_TEST(hopf_pinkall_quadric_matches_the_jet) {
  const real kappa = 2;
  const real u = real(0.3);
  const real v = real(-0.2);
  const real h = hopf_pinkall_height(kappa, u, v);
  const real expect = kappa * u * u + k_hopf_sphere_curvature * u * v;
  GAUDI_EXPECT(std::abs(h - expect) < real(1e-12));
  GAUDI_EXPECT(std::abs(hopf_pinkall_height(kappa, 0, 0)) < real(1e-12));
  const real H = real(-1.5);
  GAUDI_EXPECT(std::abs(hopf_elastic_kappa(H, real(0.2), 6, 0, 0) - H) <
               real(1e-8));
  const real shifted = hopf_elastic_kappa(H, real(0.2), 6, 0, real(0.05));
  GAUDI_EXPECT(std::abs(shifted - H) > real(1e-3));
  const vec3 spoke = u * vec3::UnitX() + v * vec3::UnitY() + h * vec3::UnitZ();
  const vec3 target =
      u * vec3::UnitX() + v * vec3::UnitY() +
      hopf_pinkall_height(kappa, spoke.x(), spoke.y()) * vec3::UnitZ();
  GAUDI_EXPECT((target - spoke).norm() < real(1e-12));
}

GAUDI_TEST(hopf_radial_land_shrinks_the_tube_circle) {
  const real k_hi = 3;
  const real kappa = real(4) / real(3);
  const real lam = hopf_tube_curvature(kappa);
  GAUDI_EXPECT(std::abs(std::abs(lam) - k_hi) < real(1e-8));
  real uo = 0;
  real ho = 0;
  GAUDI_EXPECT(hopf_radial_land(k_hi, kappa, 0, 0, uo, ho));
  GAUDI_EXPECT(std::abs(uo) < real(1e-8));
  GAUDI_EXPECT(std::abs(ho) < real(1e-6));

  const real grown = real(1.9);
  const real rho = real(1) / std::abs(hopf_tube_curvature(grown));
  const real ch = real(1) / k_hi;
  GAUDI_EXPECT(hopf_radial_land(k_hi, grown, 0, 0, uo, ho));
  GAUDI_EXPECT(std::abs(uo) < real(1e-8));
  GAUDI_EXPECT(std::abs(ho - (ch - rho)) < real(1e-6));
  GAUDI_EXPECT(hopf_radial_land(k_hi, grown, 0, ch - rho, uo, ho));
  GAUDI_EXPECT(std::abs(ho - (ch - rho)) < real(1e-6));

  real x = 0;
  real y = 0;
  hopf_circle_at(k_hi, real(0.15), x, y);
  GAUDI_EXPECT(hopf_radial_land(k_hi, grown, x, y, uo, ho));
  const real du = uo;
  const real dh = ho - ch;
  GAUDI_EXPECT(std::abs(std::sqrt(du * du + dh * dh) - rho) < real(1e-6));
}

GAUDI_TEST(hopf_ring_phase_reads_the_transport_angle) {
  const vec3 N = vec3::UnitZ();
  const vec3 V = vec3::UnitX();
  hopf_ring_node nodes[6];
  for (int k = 0; k < 6; ++k) {
    const real ang = real(k) * M_PI / real(3);
    nodes[k].t =
        real(0.2) * (std::cos(ang) * vec3::UnitX() + std::sin(ang) * vec3::UnitY());
    nodes[k].spoke_z = nodes[k].t.dot(V);
  }
  for (int k = 0; k < 6; ++k)
    nodes[k].link_z = nodes[(k + 1) % 6].spoke_z - nodes[k].spoke_z;

  GAUDI_EXPECT(std::abs(hopf_ring_phase(N, V, nodes, 6)) < real(1e-5));

  const real deg = real(30) * M_PI / real(180);
  const vec3 L = std::cos(deg) * vec3::UnitX() + std::sin(deg) * vec3::UnitY();
  GAUDI_EXPECT(std::abs(hopf_ring_phase(N, L, nodes, 6) - deg) < real(1e-4));

  hopf_ring_node flipped[6];
  for (int k = 0; k < 6; ++k)
    flipped[k] = nodes[k];
  for (int k = 0; k < 6; k += 2)
    flipped[k].spoke_z += M_PI;
  for (int k = 0; k < 6; ++k)
    flipped[k].link_z = flipped[(k + 1) % 6].spoke_z - flipped[k].spoke_z;
  GAUDI_EXPECT(std::abs(hopf_ring_phase(N, L, flipped, 6) - deg) < real(1e-4));

  for (int k = 0; k < 6; ++k) {
    nodes[k].spoke_z = 0;
    nodes[k].link_z = 0;
  }
  GAUDI_EXPECT(std::abs(hopf_ring_phase(N, L, nodes, 6)) < real(1e-8));
}

GAUDI_TEST(hopf_cylinder_land_keeps_the_measured_center) {
  const real k_hi = 3;
  real uo = 0;
  real ho = 0;
  GAUDI_EXPECT(hopf_cylinder_land(k_hi, 1, 0, 0, uo, ho));
  GAUDI_EXPECT(std::abs(uo) < real(1e-8));
  GAUDI_EXPECT(std::abs(ho) < real(1e-8));

  // A looser circle moves the center out. The foot stays on the near side.
  GAUDI_EXPECT(hopf_cylinder_land(k_hi, real(0.5), 0, 0, uo, ho));
  GAUDI_EXPECT(std::abs(uo) < real(1e-8));
  GAUDI_EXPECT(std::abs(ho) < real(1e-8));

  GAUDI_EXPECT(hopf_cylinder_land(k_hi, real(1.5), 0, 0, uo, ho));
  GAUDI_EXPECT(std::abs(uo) < real(1e-8));
  GAUDI_EXPECT(std::abs(ho) < real(1e-8));

  GAUDI_EXPECT(hopf_cylinder_land(k_hi, 1, real(0.1), 0, uo, ho));
  GAUDI_EXPECT(ho > real(0));

  GAUDI_EXPECT(hopf_cylinder_land(-k_hi, 1, real(0.1), 0, uo, ho));
  GAUDI_EXPECT(ho < real(0));

  GAUDI_EXPECT(hopf_cylinder_land(k_hi, -1, real(0.1), 0, uo, ho));
  GAUDI_EXPECT(ho < real(0));
}

GAUDI_TEST(hopf_aligned_phase_follows_the_rotated_ruling) {
  const real minor = real(0.4);
  const real major = real(2);
  GAUDI_EXPECT(std::abs(hopf_aligned_phase(minor, major, 0, 1, 1, 0) - minor) <
               real(1e-5));
  const real turned =
      hopf_aligned_phase(minor, major, real(M_PI) / 2, 2, 1, 0);
  const real expect = std::atan2(std::sin(real(-4)), std::cos(real(-4)));
  GAUDI_EXPECT(std::abs(turned - expect) < real(1e-5));
}

GAUDI_TEST(hopf_section_curvature_follows_the_rotated_tube) {
  GAUDI_EXPECT(std::abs(hopf_section_curvature(3, 1, 0) - 3) < real(1e-8));
  GAUDI_EXPECT(std::abs(hopf_section_curvature(3, 1, real(M_PI) / 2) - 1) <
               real(1e-8));
  GAUDI_EXPECT(std::abs(hopf_section_curvature(3, 1, real(M_PI) / 4) - 2) <
               real(1e-8));
}

GAUDI_TEST(hopf_wave_move_is_only_the_radius_change) {
  real du = 1;
  real dh = 1;
  GAUDI_EXPECT(hopf_wave_move(3, 1, real(0.2), real(0.1), du, dh));
  GAUDI_EXPECT(std::abs(du) < real(1e-8));
  GAUDI_EXPECT(std::abs(dh) < real(1e-8));

  GAUDI_EXPECT(hopf_wave_move(3, real(1.5), 0, 0, du, dh));
  GAUDI_EXPECT(std::abs(du) < real(1e-8));
  GAUDI_EXPECT(std::abs(dh) < real(1e-8));
}

GAUDI_TEST(hopf_wave_scale_swings_from_minus_R_to_R) {
  GAUDI_EXPECT(std::abs(hopf_wave_scale(1, real(M_PI) / 2) - 1) < real(1e-8));
  GAUDI_EXPECT(std::abs(hopf_wave_scale(1, -real(M_PI) / 2) + 1) < real(1e-8));
}

GAUDI_TEST(hopf_curvature_target_adds_sin_x) {
  const real k0 = 3;
  GAUDI_EXPECT(std::abs(hopf_curvature_plus_sin(k0, 0) - k0) < real(1e-8));
  GAUDI_EXPECT(std::abs(hopf_curvature_plus_sin(k0, real(M_PI) / 2) - (k0 + 1)) <
               real(1e-8));
  GAUDI_EXPECT(std::abs(hopf_curvature_plus_sin(k0, -real(M_PI) / 2) - (k0 - 1)) <
               real(1e-8));

  real du = 1;
  real dh = 1;
  GAUDI_EXPECT(hopf_curvature_move(k0, k0, real(0.2), real(0.1), du, dh));
  GAUDI_EXPECT(std::abs(du) < real(1e-8));
  GAUDI_EXPECT(std::abs(dh) < real(1e-8));

  GAUDI_EXPECT(hopf_curvature_move(k0, k0 + 1, 0, 0, du, dh));
  GAUDI_EXPECT(std::abs(du) < real(1e-8));
  GAUDI_EXPECT(std::abs(dh) < real(1e-8));

  GAUDI_EXPECT(hopf_curvature_move(k0, k0 + 1, real(0.2), 0, du, dh));
  GAUDI_EXPECT(du * du + dh * dh > real(1e-8));
}

GAUDI_TEST(hopf_curvature_offset_adds_sin_x_to_mean_curvature) {
  const real area = real(0.01);
  const real sum_k = 8;
  const real beta =
      hopf_curvature_normal_offset(area, sum_k, real(M_PI) / 2);
  GAUDI_EXPECT(std::abs(beta - real(0.005)) < real(1e-12));
  GAUDI_EXPECT(std::abs(sum_k * beta - real(4) * area) < real(1e-12));
  const real neg =
      hopf_curvature_normal_offset(area, sum_k, -real(M_PI) / 2);
  GAUDI_EXPECT(std::abs(neg + real(0.005)) < real(1e-12));
  GAUDI_EXPECT(std::abs(hopf_curvature_normal_offset(area, sum_k, 0)) <
               real(1e-12));
}

GAUDI_TEST(hopf_revolution_tube_is_the_minor_direction) {
  vec3 major = vec3::Zero();
  vec3 tube = vec3::Zero();
  GAUDI_EXPECT(hopf_revolution_axes(vec3(real(1.33), 0, 0), vec3::UnitX(), major,
                                    tube));
  GAUDI_EXPECT((major - vec3::UnitY()).norm() < real(1e-6));
  GAUDI_EXPECT((tube - vec3::UnitZ()).norm() < real(1e-6));
  const real ang = hopf_normal_angle(vec3::UnitX(), major);
  GAUDI_EXPECT(std::abs(std::abs(ang) - real(M_PI) / 2) < real(1e-5));
}

GAUDI_TEST(hopf_section_phase_reads_the_osculating_plane) {
  const vec3 N = vec3::UnitZ();
  const vec3 fiber = vec3::UnitX();
  const vec3 tube = vec3::UnitY();
  // Z × (x, y, 0) = +X, so this fiber is already the major tangent.
  const vec3 point(0, -1, 0);
  const real beta = real(0.4);
  const quat turned(Eigen::AngleAxis<real>(beta, fiber));
  const real expect = beta;
  hopf_section_sample samples[5];
  const real angs[5] = {real(-2.7), real(-0.4), real(0.3), real(1.5), real(2.8)};
  for (int k = 0; k < 5; ++k) {
    samples[k].t = std::cos(angs[k]) * fiber + std::sin(angs[k]) * tube;
    samples[k].q = turned;
  }
  samples[3].q = quat(Eigen::AngleAxis<real>(M_PI, N));
  const real phase = hopf_section_phase(point, N, fiber, tube, samples, 5);
  GAUDI_EXPECT(std::abs(phase - expect) < real(1e-5));
  const real flipped =
      hopf_section_phase(point, N, -fiber, -tube, samples, 5);
  GAUDI_EXPECT(std::abs(flipped - expect) < real(1e-5));
  GAUDI_EXPECT(std::abs(hopf_phase_shift(phase, beta) - (expect - beta)) <
               real(1e-5));
  GAUDI_EXPECT(std::abs(hopf_phase_shift(flipped, beta) - (expect - beta)) <
               real(1e-5));
}

GAUDI_TEST(hyperboloid_x_axis_uses_k_plus_sin_x) {
  hepworth::block::shape_w W;
  W.k11 = 2;
  W.k22 = -1;
  W.ok = true;
  const vec3 T = vec3::UnitX();
  const vec3 B = vec3::UnitY();
  const vec3 N = vec3::UnitZ();
  const real x = real(M_PI) / 2;
  const hepworth::block::shape_w Wh =
      hepworth::block::shape_with_x_axis_curvature(W, T, B, N, x);
  const real hx = hepworth::block::monge_fairing_height(
      hepworth::block::monge_shape::hyperboloid, Wh, 1, 0, nullptr, T, B, N);
  const real hy = hepworth::block::monge_fairing_height(
      hepworth::block::monge_shape::hyperboloid, Wh, 0, 1, nullptr, T, B, N);
  GAUDI_EXPECT(std::abs(hx - real(1.5)) < real(1e-8));
  GAUDI_EXPECT(std::abs(hy + real(0.5)) < real(1e-8));
}

GAUDI_TEST(hopf_line_parabola_is_half_k_cos_phase) {
  hepworth::block::shape_w W;
  W.k11 = 2;
  W.k22 = -1;
  W.ok = true;
  const vec3 point(1, 0, 0);
  const vec3 N = vec3::UnitZ();
  vec3 major = vec3::Zero();
  vec3 tube = vec3::Zero();
  GAUDI_ASSERT(hepworth::block::hopf_revolution_axes(point, N, major, tube));
  const vec3 T = tube;
  const vec3 B = major;
  const real phase = hepworth::block::hopf_normal_angle(N, major);
  const hepworth::block::shape_w Wh =
      hepworth::block::shape_with_hopf_line_curvature(W, T, B, N, point);
  const real hx = hepworth::block::monge_fairing_height(
      hepworth::block::monge_shape::quadric, Wh, 1, 0, nullptr, T, B, N);
  const real hy = hepworth::block::monge_fairing_height(
      hepworth::block::monge_shape::quadric, Wh, 0, 1, nullptr, T, B, N);
  GAUDI_EXPECT(std::abs(hx - real(0.5) * real(2) * std::cos(phase)) < real(1e-8));
  GAUDI_EXPECT(std::abs(hy) < real(1e-8));
  const hepworth::block::shape_w W2 =
      hepworth::block::shape_with_hopf_line_curvature(W, T, B, N, point, 2);
  const real hx2 = hepworth::block::monge_fairing_height(
      hepworth::block::monge_shape::quadric, W2, 1, 0, nullptr, T, B, N);
  GAUDI_EXPECT(std::abs(hx2 - real(0.5) * real(2) * std::cos(real(2) * phase)) <
               real(1e-8));
  const hepworth::block::hopf_geodesic turned =
      hepworth::block::hopf_phase_geodesic(W, T, B, N, point, 1,
                                           real(0.5) * M_PI);
  GAUDI_ASSERT(turned.ok);
  GAUDI_EXPECT(std::abs(turned.tube.dot(tube)) < real(1e-8));
  GAUDI_EXPECT(std::abs(std::abs(turned.tube.dot(major)) - real(1)) < real(1e-8));
}

GAUDI_TEST(fiber_skew_splits_mean_curvature_by_the_cross_ratio) {
  hepworth::block::shape_w W;
  W.k1 = 2;
  W.k2 = -1;
  W.e1x = 1;
  W.e1y = 0;
  W.e2x = 0;
  W.e2y = 1;
  W.ok = true;
  const real H = real(0.5) * (W.k1 + W.k2);
  const vec3 T = vec3::UnitX();
  const vec3 B = vec3::UnitY();
  const vec3 N = vec3::UnitZ();
  const real skew = real(1.01);
  real k_u = 0;
  real k_v = 0;
  hepworth::block::cross_ratio_curvatures(H, skew, k_u, k_v);
  GAUDI_EXPECT(std::abs(k_u / k_v - skew) < real(1e-8));
  GAUDI_EXPECT(std::abs(real(0.5) * (k_u + k_v) - H) < real(1e-8));
  real ku = 0;
  real kv = 0;
  hepworth::block::hyperboloid_curvature_pair(k_u, k_v, ku, kv);
  vec3 u = vec3::Zero();
  vec3 v = vec3::Zero();
  GAUDI_ASSERT(hepworth::block::principal_axes_twisted(W, T, B, N, 0, u, v));
  const hepworth::block::shape_w saddle_line =
      hepworth::block::shape_along_tangent(T, B, N, u, ku);
  const hepworth::block::shape_w saddle =
      hepworth::block::shape_with_direction_curvature(saddle_line, T, B, N, v, kv);
  const real hx = hepworth::block::monge_fairing_height(
      hepworth::block::monge_shape::quadric, saddle, 1, 0, nullptr, T, B, N);
  const real hy = hepworth::block::monge_fairing_height(
      hepworth::block::monge_shape::quadric, saddle, 0, 1, nullptr, T, B, N);
  GAUDI_EXPECT(std::abs(hx - real(0.5) * ku) < real(1e-8));
  GAUDI_EXPECT(std::abs(hy - real(0.5) * kv) < real(1e-8));
  GAUDI_EXPECT(ku * kv < real(0));
  vec3 u_t = vec3::Zero();
  vec3 v_t = vec3::Zero();
  GAUDI_ASSERT(hepworth::block::principal_axes_twisted(W, T, B, N,
                                                      real(0.5) * M_PI, u_t, v_t));
  const hepworth::block::shape_w line_t =
      hepworth::block::shape_along_tangent(T, B, N, u_t, ku);
  const hepworth::block::shape_w turned =
      hepworth::block::shape_with_direction_curvature(line_t, T, B, N, v_t, kv);
  const real along_x = hepworth::block::curvature_along(turned, T, B, N, T);
  const real along_y = hepworth::block::curvature_along(turned, T, B, N, B);
  GAUDI_EXPECT(std::abs(along_x - kv) < real(1e-8));
  GAUDI_EXPECT(std::abs(along_y - ku) < real(1e-8));
}

GAUDI_TEST(conformal_hopf_lands_on_the_circle) {
  hepworth::block::shape_w W;
  W.k11 = 2;
  W.k22 = -1;
  W.ok = true;
  const vec3 point(1, 0, 0);
  const vec3 N = vec3::UnitZ();
  const vec3 T = vec3::UnitX();
  const vec3 B = vec3::UnitY();
  const hepworth::block::hopf_geodesic g =
      hepworth::block::hopf_phase_geodesic(W, T, B, N, point);
  GAUDI_ASSERT(g.ok);
  vec3 on_tube = vec3::Zero();
  GAUDI_ASSERT(hepworth::block::conformal_hopf_target(g, N, g.tube, on_tube));
  real uo = 0;
  real ho = 0;
  GAUDI_ASSERT(hepworth::block::hopf_curvature_land(g.kappa, 1, 0, uo, ho));
  GAUDI_EXPECT((on_tube - (uo * g.tube + ho * N)).norm() < real(1e-8));
  vec3 on_major = vec3::Zero();
  GAUDI_ASSERT(hepworth::block::conformal_hopf_target(g, N, g.major, on_major));
  GAUDI_EXPECT((on_major - g.major).norm() < real(1e-8));
}

} // namespace test
} // namespace gaudi

#endif
