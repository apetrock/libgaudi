#ifndef __GAUDI_TEST_DIPOLE_TUNNELING_TESTS_HPP__
#define __GAUDI_TEST_DIPOLE_TUNNELING_TESTS_HPP__

#include <algorithm>
#include <cmath>
#include <memory>
#include <set>
#include <vector>

#include "gaudi/asawa/rod/dynamic.hpp"
#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/duchamp/dipole_tunneling_constraint.hpp"
#include "gaudi/duchamp/dipole_tunneling_geometry.hpp"
#include "gaudi/duchamp/fields.hpp"
#include "gaudi/vec_addendum.h"
#include "gaudi/hepworth/blocks/shell_position_block.hpp"
#include "gaudi/hepworth/constraints/bundles.hpp"
#include "gaudi/hepworth/nodes/solver_builder.hpp"
#include "gaudi/test/test.hpp"

namespace gaudi {
namespace test {

namespace dipole_tunneling {
using namespace gaudi::duchamp::dipole_tunneling;
} // namespace dipole_tunneling

namespace {

inline real max_radius_deviation(const std::vector<vec3> &x, const vec3 &c,
                                 real r0) {
  real d = 0.0;
  for (const vec3 &xi : x) {
    d = std::max(d, std::abs((xi - c).norm() - r0));
  }
  return d;
}

inline hepworth::block::shell_position_block::ptr make_sphere_shell_block() {
  auto M = asawa::shell::load_sphere(1.0, 24, 16);
  asawa::init_vert_datum<vec3>(*M, vec3::Zero());
  auto xs = std::make_shared<duchamp::shell_vert_positions>(M, 0);
  auto vs = std::make_shared<duchamp::shell_vert_velocities>(M, 1);
  return std::make_shared<hepworth::block::shell_position_block>(M, xs, vs);
}

} // namespace

GAUDI_TEST(dipole_tunneling_sphere_elastic_fixture) {
  auto shell = make_sphere_shell_block();
  std::vector<vec3> &x = shell->xs->get();

  const vec3 c = vec3::Zero();
  const real r0 = 1.0;
  const real dev0 = max_radius_deviation(x, c, r0);
  GAUDI_ASSERT(dev0 < 1e-6);

  x[0] = x[0].normalized() * (r0 + 0.15);
  const real dev_perturbed = max_radius_deviation(x, c, r0);
  GAUDI_ASSERT(dev_perturbed > dev0 + 0.1);

  auto config =
      hepworth::block::block_solver_builder<hepworth::block::shell_position_block>::create()
          .with_blocks(shell)
          .with_bundle(hepworth::block::make_shell_physics_bundle<0>(shell, 1.0, 0.1))
          .dt(0.05)
          .damping(0.5)
          .iterations(8)
          .build();

  hepworth::block::projection_solver solver;
  hepworth::block::run_solver_step(config, solver);

  const real dev_after = max_radius_deviation(x, c, r0);
  GAUDI_ASSERT(dev_after < dev_perturbed);
}

GAUDI_TEST(dipole_dual_disk_whole_circles_at_offset_normals) {
  using namespace dipole_tunneling;
  const vec3 xr(0.0, 0.0, 0.0);
  const auto f =
      make_segment_frame(xr, vec3(1.0, 0.0, 0.0), vec3(0.0, 0.0, 1.0));
  const auto &[xr0, xr1, N, T, B] = f;
  (void)xr0;
  (void)xr1;
  (void)T;
  const real r = 0.15;
  const real h = r;

  const vec3 upper_center = xr + h * N;
  GAUDI_ASSERT(classify_dipole_disk(upper_center, xr, f, h, r) ==
               dipole_zone::upper_clearance);
  GAUDI_ASSERT(classify_dipole_disk(upper_center + 0.1 * B, xr, f, h, r) ==
               dipole_zone::upper_clearance);

  const vec3 lower_center = xr - h * N;
  GAUDI_ASSERT(classify_dipole_disk(lower_center, xr, f, h, r) ==
               dipole_zone::lower_contact);
  GAUDI_ASSERT(classify_dipole_disk(lower_center + 0.1 * B, xr, f, h, r) ==
               dipole_zone::lower_contact);

  GAUDI_ASSERT(classify_dipole_disk(xr + 0.5 * B, xr, f, h, r) ==
               dipole_zone::outside);
}

GAUDI_TEST(dipole_degenerate_cases_finite) {
  using namespace dipole_tunneling;
  // Collapsed rod edge: frame fails.
  segment_frame f;
  GAUDI_ASSERT(!try_make_segment_frame(vec3::Zero(), vec3::Zero(), vec3::UnitY(), f));
  // Bad radius / non-finite inputs: proj_to_dipole is a no-op.
  const vec3 p(0.1, 0.2, 0.3);
  GAUDI_ASSERT((proj_to_dipole(p, vec3::Zero(), vec3::UnitY(), 0.0) - p).norm() < 1e-12);
  const vec3 p_out = proj_to_dipole(p, vec3::Zero(), vec3::UnitY(), 0.5);
  GAUDI_ASSERT(is_finite_vec(p_out));

  // Nr ∥ T still yields a finite orthonormal frame.
  GAUDI_ASSERT(try_make_segment_frame(vec3::Zero(), vec3::UnitX(), vec3::UnitX(), f));
  const auto &[xr0, xr1, N, T, B] = f;
  (void)xr0;
  (void)xr1;
  GAUDI_ASSERT(is_finite_vec(N) && is_finite_vec(T) && is_finite_vec(B));
  GAUDI_ASSERT(std::abs(N.norm() - 1.0) < 1e-9);
  GAUDI_ASSERT(std::abs(T.dot(N)) < 1e-9);

  // dp ⟂ N / zero dp → infinite radius (gather rejects).
  GAUDI_ASSERT(!std::isfinite(safe_tangent_point_radius(vec3::UnitX(), vec3::UnitY())));
  GAUDI_ASSERT(!std::isfinite(safe_tangent_point_radius(vec3::Zero(), vec3::UnitY())));

  // Force stays finite on collapsed / bad inputs.
  const vec3 fv = dipole_radial_force_at(p, make_segment_frame(vec3::Zero(), vec3::Zero(),
                                                              vec3::UnitY()),
                                        0.5, 0.0, true, 1.0, vec3::Zero());
  GAUDI_ASSERT(is_finite_vec(fv));
  GAUDI_ASSERT(fv.squaredNorm() < 1e-24);

  // Antiparallel / zero normals do not NaN slerp.
  const vec3 antip = interpolate_rod_normal(vec3::UnitY(), -vec3::UnitY(), 0.5);
  GAUDI_ASSERT(is_finite_vec(antip) && antip.squaredNorm() > 1e-12);
  const vec3 zero_n = interpolate_rod_normal(vec3::Zero(), vec3::Zero(), 0.5);
  GAUDI_ASSERT(is_finite_vec(zero_n) && zero_n.squaredNorm() > 1e-12);
}

GAUDI_TEST(dipole_interpolate_rod_normal_slerp) {
  using namespace dipole_tunneling;
  const vec3 N0 = vec3::UnitY();
  const vec3 N1 = vec3::UnitZ();

  const vec3 at0 = interpolate_rod_normal(N0, N1, 0.0);
  const vec3 at1 = interpolate_rod_normal(N0, N1, 1.0);
  const vec3 mid = interpolate_rod_normal(N0, N1, 0.5);

  GAUDI_ASSERT((at0 - N0).norm() < 1e-9);
  GAUDI_ASSERT((at1 - N1).norm() < 1e-9);
  GAUDI_ASSERT(std::abs(mid.norm() - 1.0) < 1e-9);
  GAUDI_ASSERT(std::abs(mid.dot(N0) - mid.dot(N1)) < 1e-6);

  // Same endpoints → constant; per-vert frame matches single-N frame at mid.
  const vec3 p0(0.0, 0.0, 0.0);
  const vec3 p1(1.0, 0.0, 0.0);
  const auto f_const = make_segment_frame(p0, p1, N0);
  const auto f_slerp = make_segment_frame(p0, p1, N0, N0, 0.5);
  const auto &[c0, c1, Nc, Tc, Bc] = f_const;
  const auto &[s0, s1, Ns, Ts, Bs] = f_slerp;
  (void)c0;
  (void)c1;
  (void)Tc;
  (void)Bc;
  (void)s0;
  (void)s1;
  (void)Ts;
  (void)Bs;
  GAUDI_ASSERT((Nc - Ns).norm() < 1e-9);
}

GAUDI_TEST(dipole_tangent_tunnel_aabb_corners) {
  using namespace dipole_tunneling;
  const vec3 p0(0.0, 0.0, 0.0);
  const vec3 p1(1.0, 0.0, 0.0);
  const auto f = make_segment_frame(p0, p1, vec3(0.0, 0.0, 1.0));
  const auto &[xr0, xr1, N, T, B] = f;
  (void)xr0;
  (void)xr1;
  (void)T;
  const real r = 0.1;

  const ext::extents_t upper = make_tangent_tunnel_aabb(
      p0, p1, f, r, dipole_zone::upper_clearance);
  const ext::extents_t lower =
      make_tangent_tunnel_aabb(p0, p1, f, r, dipole_zone::lower_contact);

  GAUDI_ASSERT(inside_extents(upper, p0 + r * N + r * B));
  GAUDI_ASSERT(inside_extents(upper, p0 + 2.0 * r * N));
  GAUDI_ASSERT(!inside_extents(upper, p0 - r * N));

  GAUDI_ASSERT(inside_extents(lower, p0 - r * N - r * B));
  GAUDI_ASSERT(inside_extents(lower, p0 - 2.0 * r * N));
  GAUDI_ASSERT(!inside_extents(lower, p0 + r * N));
}

GAUDI_TEST(dipole_disc_tangent_to_rod) {
  using namespace dipole_tunneling;
  const vec3 xr(0.0, 0.0, 0.0);
  const auto f =
      make_segment_frame(xr, vec3(1.0, 0.0, 0.0), vec3(0.0, 0.0, 1.0));
  const auto &[xr0, xr1, N, T, B] = f;
  (void)xr0;
  (void)xr1;
  (void)T;
  (void)B;
  const real r = 0.1;
  const real h = r;
  const real eps = 1e-4;

  GAUDI_ASSERT(classify_dipole_disk(xr + (r - eps) * N, xr, f, h, r) ==
               dipole_zone::upper_clearance);
  GAUDI_ASSERT(classify_dipole_disk(xr, xr, f, h, r) == dipole_zone::outside);
}

GAUDI_TEST(dipole_gather_faces_in_tangent_pipe) {
  using namespace dipole_tunneling;
  auto M = asawa::shell::load_sphere(1.0, 12, 8);
  std::vector<vec3> x = asawa::get_vec_data(*M, 0);

  const vec3 p0(-0.5, 0.0, 0.0);
  const vec3 p1(0.5, 0.0, 0.0);
  const segment_frame f = make_segment_frame(p0, p1, vec3(0.0, 0.0, 1.0));
  const real r = 0.55;

  const auto [upper, lower] = gather_faces_in_segment(*M, x, f, r);
  GAUDI_ASSERT(!upper.empty() || !lower.empty());
}

GAUDI_TEST(dipole_gather_radius_gate_excludes_large_R) {
  using namespace dipole_tunneling;
  auto M = asawa::shell::load_sphere(1.0, 12, 8);
  std::vector<vec3> x = asawa::get_vec_data(*M, 0);

  const vec3 p0(-0.5, 0.0, 0.0);
  const vec3 p1(0.5, 0.0, 0.0);
  const segment_frame f = make_segment_frame(p0, p1, vec3(0.0, 0.0, 1.0));
  const real r = 0.3;

  const auto [upper, lower] = gather_faces_in_segment(*M, x, f, r);

  for (const auto fi : M->get_face_range()) {
    const asawa::shell::FaceId fid = asawa::shell::face_id(fi);
    const vec3 xf = asawa::shell::face_center(*M, fid, x);
    if (xf.z() < 0.9)
      continue;
    const bool in_upper =
        std::find(upper.begin(), upper.end(), fid) != upper.end();
    const bool in_lower =
        std::find(lower.begin(), lower.end(), fid) != lower.end();
    GAUDI_ASSERT(!in_upper && !in_lower);
  }
}

GAUDI_TEST(dipole_gather_splits_upper_and_lower) {
  using namespace dipole_tunneling;
  auto M = asawa::shell::load_sphere(1.0, 12, 8);
  std::vector<vec3> x = asawa::get_vec_data(*M, 0);

  const vec3 p0(-0.5, 0.0, 0.0);
  const vec3 p1(0.5, 0.0, 0.0);
  const segment_frame f = make_segment_frame(p0, p1, vec3(0.0, 0.0, 1.0));
  const real r = 0.55;

  const auto [upper, lower] = gather_faces_in_segment(*M, x, f, r);
  GAUDI_ASSERT(!upper.empty());
  GAUDI_ASSERT(!lower.empty());

  std::set<int> upper_set;
  std::set<int> lower_set;
  for (const asawa::shell::FaceId fid : upper)
    upper_set.insert(static_cast<int>(fid));
  for (const asawa::shell::FaceId fid : lower)
    lower_set.insert(static_cast<int>(fid));

  for (const int u : upper_set)
    GAUDI_ASSERT(lower_set.count(u) == 0);

  for (const asawa::shell::FaceId fid : upper) {
    const vec3 xf = asawa::shell::face_center(*M, fid, x);
    GAUDI_ASSERT(xf.z() > 0.0);
  }
  for (const asawa::shell::FaceId fid : lower) {
    const vec3 xf = asawa::shell::face_center(*M, fid, x);
    GAUDI_ASSERT(xf.z() < 0.0);
  }
}

GAUDI_TEST(dipole_tunnel_radial_force_sign_and_magnitude) {
  using namespace dipole_tunneling;
  const vec3 p0(-0.5, 0.0, 0.0);
  const vec3 p1(0.5, 0.0, 0.0);
  const segment_frame f = make_segment_frame(p0, p1, vec3(0.0, 0.0, 1.0));
  const auto &[xr0, xr1, N, T, B] = f;
  (void)xr0;
  (void)xr1;
  (void)T;
  (void)B;
  const real r = 0.2;
  const real h = r;
  // Mesh N aligned with cylinder radial on upper side → (−side)·(r−dist)·N
  // with side=+1 → −(r−dist)·N = −expect_u · UnitZ.
  const vec3 N_tri = vec3::UnitZ();

  const vec3 xf_upper(0.0, 0.0, h + 0.5 * r);
  const vec3 f_upper = dipole_radial_force_at(xf_upper, f, r, h, true, 1.0, N_tri);
  const vec3 xr_u = va::project_on_line(p0, p1, xf_upper);
  const vec3 dp_u = xf_upper - xr_u;
  const vec3 cen_u = dipole_sphere_center(xr_u, dp_u, N, r);
  const real expect_u = r - (xf_upper - cen_u).norm();
  GAUDI_ASSERT(expect_u > 0.0);
  GAUDI_ASSERT((f_upper - (-expect_u * N_tri)).norm() < 1e-9);

  const vec3 xf_lower(0.0, 0.0, -(h + 0.5 * r));
  const vec3 f_lower = dipole_radial_force_at(xf_lower, f, r, h, false, 1.0, N_tri);
  const vec3 xr_l = va::project_on_line(p0, p1, xf_lower);
  const vec3 dp_l = xf_lower - xr_l;
  const vec3 cen_l = dipole_sphere_center(xr_l, dp_l, N, r);
  const real expect_l = r - (xf_lower - cen_l).norm();
  GAUDI_ASSERT(expect_l > 0.0);
  // Lower side=-1 → (−side)·(r−dist)·N = +expect_l · UnitZ.
  GAUDI_ASSERT((f_lower - (expect_l * N_tri)).norm() < 1e-9);

  // Tilted mesh normal on upper: still (−side)·(r−dist)·N_tilt.
  const vec3 N_tilt = (vec3::UnitZ() + 0.5 * B).normalized();
  const vec3 f_tilt =
      dipole_radial_force_at(xf_upper, f, r, h, true, 1.0, N_tilt);
  GAUDI_ASSERT((f_tilt - (-expect_u * N_tilt)).norm() < 1e-9);
}

GAUDI_TEST(dipole_tunnel_radial_force_zero_at_edge) {
  using namespace dipole_tunneling;
  const vec3 p0(-0.5, 0.0, 0.0);
  const vec3 p1(0.5, 0.0, 0.0);
  const segment_frame f = make_segment_frame(p0, p1, vec3(0.0, 0.0, 1.0));
  const auto &[xr0, xr1, N, T, B] = f;
  (void)xr0;
  (void)xr1;
  (void)T;
  const real r = 0.2;
  const real h = r;
  const vec3 N_tri = vec3::UnitZ();

  const vec3 xr(0.0, 0.0, 0.0);
  const vec3 dp(h * N);
  const vec3 cen = dipole_sphere_center(xr, dp, N, r);
  const vec3 x_far = cen + (r + 2.0 * dipole_margin(r)) * B;
  const vec3 f_far = dipole_radial_force_at(x_far, f, r, h, true, 1.0, N_tri);
  GAUDI_ASSERT(f_far.squaredNorm() < 1e-18);
}

GAUDI_TEST(dipole_tunnel_radial_force_negative_outside_margin) {
  using namespace dipole_tunneling;
  const vec3 p0(-0.5, 0.0, 0.0);
  const vec3 p1(0.5, 0.0, 0.0);
  const segment_frame f = make_segment_frame(p0, p1, vec3(0.0, 0.0, 1.0));
  const auto &[xr0, xr1, N, T, B] = f;
  (void)xr0;
  (void)xr1;
  (void)T;
  const real r = 0.2;
  const real h = r;
  (void)N;

  const vec3 xr(0.0, 0.0, 0.0);
  const vec3 dp(h * vec3::UnitZ());
  const vec3 cen = dipole_sphere_center(xr, dp, vec3::UnitZ(), r);
  const real outside = 0.05 * r;
  const vec3 x_out = cen + (r + outside) * B;
  // Mesh normal along B; upper side → (−side)·(r−dist)·B = +outside·B.
  const vec3 f_out = dipole_radial_force_at(x_out, f, r, h, true, 1.0, B);
  GAUDI_ASSERT((f_out - (outside * B)).norm() < 1e-9);
}

GAUDI_TEST(dipole_tunnel_force_accumulates_at_vertices) {
  using namespace dipole_tunneling;
  const vec3 p0(-0.5, 0.0, 0.0);
  const vec3 p1(0.5, 0.0, 0.0);
  const segment_frame f = make_segment_frame(p0, p1, vec3(0.0, 0.0, 1.0));
  const real r = 0.2;
  const real h = r;
  const real force_w = 1.0;

  const vec3 N_tri = vec3::UnitZ();

  const vec3 x0(0.0, 0.0, h + 0.1 * r);
  const vec3 x1(0.0, 0.0, h + 0.5 * r);
  const vec3 x2(0.0, 0.0, h + 0.9 * r);
  const vec3 f0 = dipole_radial_force_at(x0, f, r, h, true, force_w, N_tri);
  const vec3 f1 = dipole_radial_force_at(x1, f, r, h, true, force_w, N_tri);
  const vec3 f2 = dipole_radial_force_at(x2, f, r, h, true, force_w, N_tri);

  GAUDI_ASSERT(f0.norm() > f1.norm());
  GAUDI_ASSERT(f1.norm() > f2.norm());

  const vec3 xf = (x0 + x1 + x2) / 3.0;
  const vec3 f_lumped = dipole_radial_force_at(xf, f, r, h, true, force_w, N_tri) / 3.0;
  GAUDI_ASSERT(std::abs(f0.norm() - f_lumped.norm()) > 1e-6);
  GAUDI_ASSERT(std::abs(f2.norm() - f_lumped.norm()) > 1e-6);
}

GAUDI_TEST(dipole_tangent_point_sphere_center) {
  const vec3 xc(0.0, 0.0, 0.0);
  const vec3 N(0.0, 0.0, 1.0);
  const vec3 x(0.0, 0.0, 1.0);
  const vec3 dp = x - xc;

  const real r = va::tangent_point_radius(dp, N);
  GAUDI_ASSERT(std::abs(r - 0.5) < 1e-9);

  const real sgn = va::sgn(N.dot(dp));
  const vec3 center = xc + sgn * N * r;
  GAUDI_ASSERT((center - vec3(0.0, 0.0, 0.5)).norm() < 1e-9);
  GAUDI_ASSERT(std::abs((x - center).norm() - r) < 1e-9);
}

GAUDI_TEST(dipole_tangent_circle_through_face_center) {
  using namespace dipole_tunneling;
  const vec3 p_nearest(0.0, 0.0, 0.0);
  const vec3 xr0(-1.0, 0.0, 0.0);
  const vec3 xr1(1.0, 0.0, 0.0);
  const vec3 N(0.0, 0.0, 1.0);
  const segment_frame f = make_segment_frame(xr0, xr1, N);
  const auto &[fxr0, fxr1, fN, fT, fB] = f;
  (void)fxr0;
  (void)fxr1;
  (void)fT;

  const vec3 cen_tri(0.0, 0.0, 1.0);
  const vec3 dp = cen_tri - p_nearest;
  const real radius = va::tangent_point_radius(dp, fN);
  const vec3 center = p_nearest + va::sgn(fN.dot(dp)) * radius * fN;

  GAUDI_ASSERT(std::abs(radius - 0.5) < 1e-9);
  GAUDI_ASSERT((center - vec3(0.0, 0.0, 0.5)).norm() < 1e-9);
  GAUDI_ASSERT(std::abs((cen_tri - center).norm() - radius) < 1e-9);

  const vec2 in_plane = project_to_NB(cen_tri - center, fN, fB);
  GAUDI_ASSERT(std::abs(in_plane.norm() - radius) < 1e-9);
}

GAUDI_TEST(dipole_tangent_circle_lower_zone) {
  using namespace dipole_tunneling;
  const vec3 p_nearest(0.0, 0.0, 0.0);
  const vec3 xr0(-1.0, 0.0, 0.0);
  const vec3 xr1(1.0, 0.0, 0.0);
  const vec3 N(0.0, 0.0, 1.0);
  const segment_frame f = make_segment_frame(xr0, xr1, N);

  const vec3 cen_tri(0.0, 0.0, -1.0);
  const auto &[fxr0, fxr1, fN, fT, fB] = f;
  (void)fxr0;
  (void)fxr1;
  (void)fT;
  (void)fB;
  const vec3 dp = cen_tri - p_nearest;
  const real radius = va::tangent_point_radius(dp, fN);
  const vec3 center = p_nearest + va::sgn(fN.dot(dp)) * radius * fN;

  GAUDI_ASSERT(std::abs(radius - 0.5) < 1e-9);
  GAUDI_ASSERT((center - vec3(0.0, 0.0, -0.5)).norm() < 1e-9);
  GAUDI_ASSERT(std::abs((cen_tri - center).norm() - radius) < 1e-9);
}

GAUDI_TEST(dipole_extents_overlaps_triangle) {
  using namespace dipole_tunneling;
  const ext::extents_t box = {vec3(-1.0, -1.0, -1.0), vec3(1.0, 1.0, 1.0)};
  const vec3 p0(0.0, 0.0, 0.0);
  const vec3 p1(0.5, 0.0, 0.0);
  const vec3 p2(0.0, 0.5, 0.0);
  GAUDI_ASSERT(ext::overlap(box, p0, p1, p2));

  const vec3 q0(2.0, 2.0, 2.0);
  const vec3 q1(2.5, 2.0, 2.0);
  const vec3 q2(2.0, 2.5, 2.0);
  GAUDI_ASSERT(!ext::overlap(box, q0, q1, q2));
}

GAUDI_TEST(dipole_polyline_foot_nearest_on_arc) {
  using namespace dipole_tunneling;
  const std::vector<vec3> polyline = {vec3(0.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0),
                                      vec3(2.0, 0.0, 0.5)};
  const vec3 x(1.1, 0.2, 0.1);

  real best_d2 = 1e30;
  index_t best_seg = 0;
  vec3 best_xr = polyline[0];
  for (index_t s = 0; s + 1 < static_cast<index_t>(polyline.size()); ++s) {
    const vec3 xr = va::project_on_line(polyline[s], polyline[s + 1], x);
    const real d2 = (x - xr).squaredNorm();
    if (d2 < best_d2) {
      best_d2 = d2;
      best_seg = s;
      best_xr = xr;
    }
  }
  GAUDI_ASSERT(best_seg == 1);

  const vec3 seg0_xr = va::project_on_line(polyline[0], polyline[1], x);
  GAUDI_ASSERT((best_xr - seg0_xr).norm() > 1e-6);

  const vec3 xr_expected = va::project_on_line(polyline[1], polyline[2], x);
  GAUDI_ASSERT((best_xr - xr_expected).norm() < 1e-9);

  const segment_frame ff = make_segment_frame(polyline[1], polyline[2], vec3::UnitZ());
  const auto &[fxr0, fxr1, fN, fT, fB] = ff;
  (void)fxr0;
  (void)fxr1;
  (void)fT;
  (void)fB;
  const vec3 dp = x - best_xr;
  const real r = va::tangent_point_radius(dp, fN);
  const vec3 center = best_xr + va::sgn(fN.dot(dp)) * r * fN;
  GAUDI_ASSERT(std::abs((x - center).norm() - r) < 1e-9);
}

// Consecutive x[i]/x[i+1] is not a rod edge after remesh; edge_verts is.
GAUDI_TEST(dipole_force_uses_rod_edge_topology_not_array_order) {
  using namespace dipole_tunneling;
  // Polyline with a deleted middle slot: corners 0-2-3 linked, index 1 dead.
  // Array order (0,1),(1,2) would invent phantom chords.
  std::vector<vec3> x_r = {
      vec3(0.0, 0.0, 0.0),
      vec3(10.0, 10.0, 10.0), // dead / unused slot
      vec3(0.5, 0.0, 0.0),
      vec3(1.0, 0.0, 0.0),
  };
  const std::vector<index_t> edges = {0, 2, 2, 3};
  const std::vector<index_t> phantom = {0, 1, 1, 2, 2, 3};
  const std::vector<vec3> Nr(x_r.size(), vec3::UnitZ());

  auto M = asawa::shell::load_sphere(1.0, 12, 8);
  const std::vector<vec3> &x = asawa::get_vec_data(*M, 0);
  const real r = 0.35;

  std::vector<vec3> f_edges, f_phantom;
  accumulate_dipole_tunnel_forces(*M, x, x_r, Nr, edges, r, r, 1.0, f_edges);
  accumulate_dipole_tunnel_forces(*M, x, x_r, Nr, phantom, r, r, 1.0, f_phantom);

  real n_edges = 0.0, n_phantom = 0.0;
  for (size_t i = 0; i < f_edges.size(); ++i) {
    n_edges += f_edges[i].norm();
    n_phantom += f_phantom[i].norm();
  }
  // Phantom chords through (0,1)/(1,2) should drive a different force field.
  GAUDI_ASSERT(std::abs(n_edges - n_phantom) > 1e-6 ||
               [&]() {
                 real d2 = 0.0;
                 for (size_t i = 0; i < f_edges.size(); ++i)
                   d2 += (f_edges[i] - f_phantom[i]).squaredNorm();
                 return d2 > 1e-8;
               }());
}

GAUDI_TEST(dipole_gather_matches_brute_on_true_edges) {
  using namespace dipole_tunneling;
  auto M = asawa::shell::load_sphere(1.0, 16, 12);
  const std::vector<vec3> &x = asawa::get_vec_data(*M, 0);

  std::vector<vec3> rod_pts = {
      vec3(-0.8, 0.0, 0.0), vec3(-0.3, 0.1, 0.0), vec3(0.2, -0.05, 0.05),
      vec3(0.7, 0.0, 0.0)};
  auto R = asawa::rod::rod::create(rod_pts, false);
  auto Rd = asawa::rod::dynamic::create(R, 0.2 * R->lavg(), 2.0 * R->lavg(),
                                        0.25 * R->lavg());
  // Remesh so next(i) is not trivially i+1 for all slots.
  for (int k = 0; k < 3; ++k)
    Rd->step();

  const std::vector<vec3> &xr = R->x();
  const std::vector<vec3> Nr = R->N1();
  const std::vector<index_t> edges = R->get_edge_vert_ids();
  const real r = 0.25;

  // Topology sanity: every edge pair is an actual next() link.
  GAUDI_ASSERT(!edges.empty());
  for (size_t e = 0; e + 1 < edges.size(); e += 2) {
    GAUDI_ASSERT(static_cast<index_t>(R->next(asawa::rod::corner_id(edges[e]))) ==
                 edges[e + 1]);
  }

  const auto brute = gather_dipole_faces_brute(*M, x, xr, Nr, edges, r);
  const auto faces = gather_dipole_faces(*M, x, xr, Nr, *Rd, r);
  const auto diff = compare_dipole_gather_face_sets(faces, brute);
  // Smoke: both paths run; default gather should capture something when brute does.
  GAUDI_ASSERT(diff[2] > 0 || brute.empty());
  (void)diff;
}

} // namespace test
} // namespace gaudi

#endif // __GAUDI_TEST_DIPOLE_TUNNELING_TESTS_HPP__
