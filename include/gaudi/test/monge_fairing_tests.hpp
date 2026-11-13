#ifndef GAUDI_MONGE_FAIRING_TESTS_HPP
#define GAUDI_MONGE_FAIRING_TESTS_HPP

#include "gaudi/common.h"
#include "gaudi/hepworth/block/osculating_dupin_ring.hpp"
#include "gaudi/hepworth/block/osculating_torus_ring.hpp"
#include "gaudi/hepworth/block/ring_variational_fit.hpp"
#include "gaudi/hepworth/block/shell_constraints.hpp"
#include "gaudi/test/test.hpp"

#include <cmath>
#include <vector>

namespace gaudi {
namespace test {
namespace {

using hepworth::block::hopf_curvature_land;
using hepworth::block::monge_cylinder_target;
using hepworth::block::monge_fairing;
using hepworth::block::monge_fairing_height;
using hepworth::block::monge_shape;
using hepworth::block::objective;
using hepworth::block::select_monge_height;
using hepworth::block::select_monge_spoke;
using hepworth::block::project_willmore_curvature;
using hepworth::block::shape_quad;
using hepworth::block::shape_w;

const vec3 kT = vec3::UnitX();
const vec3 kB = vec3::UnitY();
const vec3 kN = vec3::UnitZ();

shape_w diag_w(real k11, real k22) {
  shape_w W;
  W.k11 = k11;
  W.k12 = 0;
  W.k22 = k22;
  W.ok = true;
  return W;
}

real height(monge_shape shape, const shape_w &W, real tx, real ty,
            const shape_w *Wj = nullptr) {
  return monge_fairing_height(shape, W, tx, ty, Wj, kT, kB, kN);
}

} // namespace

GAUDI_TEST(monge_quadric_matches_jet_plus_cubic) {
  const shape_w Wi = diag_w(1.0, 0.5);
  shape_w Wj = diag_w(1.4, 0.2);
  Wj.k12 = 0.1;
  const real tx = 0.2;
  const real ty = 0.3;
  const real h2 = real(0.5) * shape_quad(Wi, tx, ty);
  const real h3 =
      real(1.0 / 6.0) * (shape_quad(Wj, tx, ty) - shape_quad(Wi, tx, ty));
  const real h = height(monge_shape::quadric, Wi, tx, ty, &Wj);
  GAUDI_EXPECT(std::abs(h - (h2 + h3)) < real(1e-12));
}

GAUDI_TEST(monge_plane_heights_are_zero) {
  const shape_w W = diag_w(0, 0);
  for (const monge_shape shape :
       {monge_shape::quadric, monge_shape::hyperboloid,
        monge_shape::cylinder_minor, monge_shape::cylinder_major,
        monge_shape::sphere, monge_shape::osculating_torus}) {
    const real h = height(shape, W, 0.2, -0.1);
    GAUDI_EXPECT(std::isfinite(h));
    GAUDI_EXPECT(std::abs(h) < real(1e-8));
  }
}

GAUDI_TEST(monge_cylinder_projects_onto_the_chosen_tube) {
  const shape_w W = diag_w(2.0, 0.5);
  const real k_minor = 2.0;
  const real k_major = 0.5;
  vec3 foot;
  GAUDI_EXPECT(monge_cylinder_target(monge_shape::cylinder_minor, W, kT, kB, kN,
                                    vec3::Zero(), foot));
  GAUDI_EXPECT(foot.norm() < real(1e-12));
  GAUDI_EXPECT(monge_cylinder_target(monge_shape::cylinder_major, W, kT, kB, kN,
                                    vec3::Zero(), foot));
  GAUDI_EXPECT(foot.norm() < real(1e-12));

  const vec3 ruling(0, 0.35, 0);
  vec3 kept;
  GAUDI_EXPECT(monge_cylinder_target(monge_shape::cylinder_minor, W, kT, kB, kN,
                                    ruling, kept));
  GAUDI_EXPECT((kept - ruling).norm() < real(1e-12));

  const real u = 0.15;
  const real v = 0.2;
  const vec3 spoke = u * kT + v * kB;
  vec3 minor_hit;
  vec3 major_hit;
  GAUDI_EXPECT(monge_cylinder_target(monge_shape::cylinder_minor, W, kT, kB, kN,
                                    spoke, minor_hit));
  GAUDI_EXPECT(monge_cylinder_target(monge_shape::cylinder_major, W, kT, kB, kN,
                                    spoke, major_hit));
  real uo = 0;
  real ho = 0;
  GAUDI_EXPECT(hopf_curvature_land(k_minor, u, 0, uo, ho));
  const vec3 minor_expect = v * kB + uo * kT + ho * kN;
  GAUDI_EXPECT((minor_hit - minor_expect).norm() < real(1e-12));
  GAUDI_EXPECT(hopf_curvature_land(k_major, v, 0, uo, ho));
  const vec3 major_expect = u * kT + uo * kB + ho * kN;
  GAUDI_EXPECT((major_hit - major_expect).norm() < real(1e-12));
  GAUDI_EXPECT((minor_hit - major_hit).norm() > real(1e-3));

  const real rho = real(1) / k_minor;
  const real uu = real(0.1);
  const real hh = rho - std::sqrt(rho * rho - uu * uu);
  const vec3 on_circle = uu * kT + real(0.3) * kB + hh * kN;
  vec3 stayed;
  GAUDI_EXPECT(monge_cylinder_target(monge_shape::cylinder_minor, W, kT, kB, kN,
                                    on_circle, stayed));
  GAUDI_EXPECT((stayed - on_circle).norm() < real(1e-12));
}

GAUDI_TEST(monge_hyperboloid_flips_same_sign_curvature) {
  const shape_w W = diag_w(2.0, 0.5);
  const real tx = 0.2;
  const real ty = 0.4;
  const real h = height(monge_shape::hyperboloid, W, tx, ty);
  const real expect = tx * tx - real(0.25) * ty * ty;
  GAUDI_EXPECT(std::abs(h - expect) < real(1e-12));
}

GAUDI_TEST(monge_sphere_is_the_umbilic_jet) {
  const shape_w W = diag_w(2.0, 0.5);
  const real tx = 0.2;
  const real ty = 0.4;
  const real H = real(0.5) * (2.0 + 0.5);
  const real expect = real(0.5) * H * (tx * tx + ty * ty);
  const real h = height(monge_shape::sphere, W, tx, ty);
  GAUDI_EXPECT(std::abs(h - expect) < real(1e-12));
}

GAUDI_TEST(monge_torus_matches_sphere_jet_through_second_order) {
  const shape_w W = diag_w(1.0, 1.0);
  const real rho = 0.05;
  const real h = height(monge_shape::osculating_torus, W, rho, 0.0);
  const real jet = real(0.5) * rho * rho;
  const real quartic = (rho * rho * rho * rho) / real(8);
  GAUDI_EXPECT(std::abs(h - jet) < real(2) * quartic);
  GAUDI_EXPECT(std::abs(h) > real(0.5) * jet);
}

GAUDI_TEST(monge_selection_picks_the_closer_graph) {
  const shape_w W = diag_w(2.0, 0.5);
  const real tx = 0.2;
  const real ty = 0.4;
  const real rho = real(0.5);
  const real uu = real(0.1);
  const real hh = rho - std::sqrt(rho * rho - uu * uu);
  const vec3 on_minor = uu * kT + ty * kB + hh * kN;
  const std::vector<objective> cyl_and_sphere{
      monge_fairing{monge_shape::cylinder_minor},
      monge_fairing{monge_shape::sphere}};
  vec3 picked;
  GAUDI_EXPECT(select_monge_spoke(cyl_and_sphere, on_minor, W, nullptr, kT, kB,
                                 kN, picked));
  GAUDI_EXPECT((picked - on_minor).norm() < real(1e-12));

  const vec3 on_quad = tx * kT + ty * kB +
                       height(monge_shape::quadric, W, tx, ty) * kN;
  const std::vector<objective> with_quad{
      monge_fairing{monge_shape::quadric},
      monge_fairing{monge_shape::cylinder_minor}};
  vec3 quad_wins;
  GAUDI_EXPECT(select_monge_spoke(with_quad, on_quad, W, nullptr, kT, kB, kN,
                                 quad_wins));
  GAUDI_EXPECT((quad_wins - on_quad).norm() < real(1e-12));
}

/// Opposite principal curvatures are the inner equator. Same sign is the outer.
GAUDI_TEST(monge_torus_sheet_follows_curvature_sign) {
  using hepworth::block::monge_torus_of;
  using hepworth::block::monge_torus_section_height;

  const real txs[] = {0.15, -0.12, 0.0, 0.08, -0.05};
  const real tys[] = {0.0, 0.04, 0.14, -0.1, -0.12};

  shape_w inner_w = diag_w(1.0, real(-0.5));
  const auto inner_tor = monge_torus_of(inner_w, kT, kB, kN);
  GAUDI_EXPECT(inner_tor.valid);
  GAUDI_EXPECT(inner_tor.R > inner_tor.r);
  for (int i = 0; i < 5; ++i) {
    const real h = height(monge_shape::osculating_torus, inner_w, txs[i], tys[i]);
    const real hi =
        monge_torus_section_height(inner_tor, txs[i], tys[i], kT, kB, kN, true);
    GAUDI_EXPECT(std::isfinite(hi));
    GAUDI_EXPECT(std::abs(h - hi) < real(1e-12));
  }
  const real apart =
      monge_torus_section_height(inner_tor, 0.0, 0.14, kT, kB, kN, true) -
      monge_torus_section_height(inner_tor, 0.0, 0.14, kT, kB, kN, false);
  GAUDI_EXPECT(std::abs(apart) > real(1e-4));

  shape_w outer_w = diag_w(1.0, real(0.25));
  const auto outer_tor = monge_torus_of(outer_w, kT, kB, kN);
  GAUDI_EXPECT(outer_tor.valid);
  for (int i = 0; i < 5; ++i) {
    const real h = height(monge_shape::osculating_torus, outer_w, txs[i], tys[i]);
    const real ho =
        monge_torus_section_height(outer_tor, txs[i], tys[i], kT, kB, kN, false);
    GAUDI_EXPECT(std::isfinite(ho));
    GAUDI_EXPECT(std::abs(h - ho) < real(1e-12));
  }
}

/// Inner sheet: center at 1/|κ_min| = R−r, revolution axis n×e_lo.
/// A point off the tube lands on the spine ray.
GAUDI_TEST(spine_torus_seats_the_inner_sheet_and_lands_on_the_ray) {
  using hepworth::block::choose_minor_tube_landing;
  using hepworth::block::land_albers_torus;
  using hepworth::block::minor_tube_fit;
  using hepworth::block::monge_torus_of;

  const shape_w inner_w = diag_w(1.0, real(-0.5));
  const auto tor = monge_torus_of(inner_w, kT, kB, kN);
  GAUDI_EXPECT(tor.valid);
  const real inv_lo = real(1) / real(0.5);
  GAUDI_EXPECT(std::abs((tor.center).norm() - inv_lo) < real(1e-8));
  GAUDI_EXPECT(std::abs((tor.R - tor.r) - inv_lo) < real(1e-8));
  GAUDI_EXPECT(std::abs(gaudi::albers::torus_sdf(vec3::Zero(), tor)) < real(1e-8));
  // Smaller curvature is along B, so the revolution axis is ±N×B.
  GAUDI_EXPECT(std::abs(std::abs(tor.axis.dot(kN.cross(kB))) - real(1)) < real(1e-8));

  const vec3 off = real(0.1) * kB + real(0.04) * kN;
  vec3 hit;
  GAUDI_EXPECT(land_albers_torus(tor, off, hit));
  GAUDI_EXPECT(std::abs(gaudi::albers::torus_sdf(hit, tor)) < real(1e-6));
  const vec3 g = gaudi::albers::torus_sdf_grad(off, tor);
  const vec3 step = hit - off;
  GAUDI_EXPECT(step.norm() > real(1e-8));
  GAUDI_EXPECT(std::abs(std::abs(step.normalized().dot(g.normalized())) - real(1)) <
               real(1e-6));

  // Vanishing major κ → large-R torus (Points as Tori), not a refused circle.
  const auto cylinder = monge_torus_of(diag_w(1.0, 0.0), kT, kB, kN);
  GAUDI_EXPECT(cylinder.valid);
  GAUDI_EXPECT(cylinder.tube);
  GAUDI_EXPECT(cylinder.R > real(1e7));
  GAUDI_EXPECT(std::abs(cylinder.r - real(1)) < real(1e-8));
  const vec3 on_line =
      cylinder.spine + real(0.2) * cylinder.ruling + cylinder.r * cylinder.radial;
  const auto line_only =
      choose_minor_tube_landing(cylinder, std::vector<vec3>{on_line});
  GAUDI_EXPECT(line_only.fit == minor_tube_fit::torus);
  GAUDI_EXPECT(line_only.ok.size() == 1 && line_only.ok[0]);
  GAUDI_EXPECT((line_only.hit[0] - on_line).norm() < real(1e-6));
}

/// Points on the minor-radius cylinder are closer to that line than to the
/// circle. Points on the circle-spine tube stay with the torus.
GAUDI_TEST(minor_tube_landing_picks_the_closer_spine) {
  using hepworth::block::choose_minor_tube_landing;
  using hepworth::block::minor_tube_fit;
  using hepworth::block::monge_torus_of;

  const auto tor = monge_torus_of(diag_w(1.0, real(0.25)), kT, kB, kN);
  GAUDI_EXPECT(tor.valid);
  GAUDI_EXPECT(tor.tube);
  GAUDI_EXPECT(tor.R > tor.r);

  std::vector<vec3> on_cyl;
  for (const real t : {real(-0.35), real(0.2), real(0.4)})
    on_cyl.push_back(tor.spine + t * tor.ruling + tor.r * tor.radial);
  const auto cyl_wins = choose_minor_tube_landing(tor, on_cyl);
  GAUDI_EXPECT(cyl_wins.fit == minor_tube_fit::cylinder);
  for (size_t i = 0; i < on_cyl.size(); ++i) {
    GAUDI_EXPECT(cyl_wins.ok[i]);
    GAUDI_EXPECT((cyl_wins.hit[i] - on_cyl[i]).norm() < real(1e-8));
  }

  const real sheet =
      ((-tor.spine).dot(tor.radial) >= real(0)) ? real(1) : real(-1);
  std::vector<vec3> on_tor;
  for (const real a : {real(-0.18), real(0.12), real(0.22)}) {
    const vec3 rad = std::cos(a) * tor.radial + std::sin(a) * tor.ruling;
    const vec3 ps = tor.center + tor.R * rad;
    on_tor.push_back(ps + sheet * tor.r * rad);
  }
  const auto tor_wins = choose_minor_tube_landing(tor, on_tor);
  GAUDI_EXPECT(tor_wins.fit == minor_tube_fit::torus);
  for (size_t i = 0; i < on_tor.size(); ++i) {
    GAUDI_EXPECT(tor_wins.ok[i]);
    GAUDI_EXPECT((tor_wins.hit[i] - on_tor[i]).norm() < real(1e-6));
  }
}

/// A large major radius stays a circle measured from the nearby spine, and
/// that landing agrees with the minor-radius line.
GAUDI_TEST(large_major_radius_lands_with_the_nearby_spine) {
  using hepworth::block::land_albers_torus;
  using hepworth::block::monge_torus_of;

  const auto tor = monge_torus_of(diag_w(1.0, real(1e-8)), kT, kB, kN);
  GAUDI_EXPECT(tor.valid);
  GAUDI_EXPECT(tor.tube);
  GAUDI_EXPECT(tor.R > real(1e7));

  const vec3 p = real(0.05) * kT + real(0.04) * kB + real(0.01) * kN;
  vec3 ps_circle;
  vec3 ps_line;
  vec3 on_circle;
  vec3 on_line;
  GAUDI_EXPECT(gaudi::albers::nearest_circle_spine(tor, p, ps_circle));
  GAUDI_EXPECT(gaudi::albers::nearest_line_spine(tor, p, ps_line));
  GAUDI_EXPECT(gaudi::albers::land_from_spine(ps_circle, tor.r, p, on_circle));
  GAUDI_EXPECT(gaudi::albers::land_from_spine(ps_line, tor.r, p, on_line));
  GAUDI_EXPECT((on_circle - on_line).norm() < real(1e-4));
  GAUDI_EXPECT((on_circle - p).norm() < real(0.05));
  GAUDI_EXPECT(std::abs((on_circle - ps_circle).norm() - tor.r) < real(1e-8));

  vec3 hit;
  GAUDI_EXPECT(land_albers_torus(tor, p, hit));
  GAUDI_EXPECT((hit - on_circle).norm() < real(1e-8));
}

/// The minor spine is the cylinder_minor center on both sheets. A far major
/// radius still keeps the torus (large-R limit), not a forced cylinder drop.
GAUDI_TEST(minor_spine_is_the_cylinder_center_and_a_far_major_keeps_torus) {
  using hepworth::block::choose_minor_tube_landing;
  using hepworth::block::minor_tube_fit;
  using hepworth::block::monge_cylinder_target;
  using hepworth::block::monge_shape;
  using hepworth::block::monge_torus_of;

  const shape_w outer_w = diag_w(1.0, real(0.25));
  const shape_w inner_w = diag_w(1.0, real(-0.5));
  const auto outer = monge_torus_of(outer_w, kT, kB, kN);
  const auto inner = monge_torus_of(inner_w, kT, kB, kN);
  GAUDI_EXPECT(outer.valid && inner.valid);
  GAUDI_EXPECT((outer.spine - kN).norm() < real(1e-8));
  GAUDI_EXPECT((inner.spine - kN).norm() < real(1e-8));
  GAUDI_EXPECT((outer.center - real(4) * kN).norm() < real(1e-8));
  GAUDI_EXPECT((inner.center + real(2) * kN).norm() < real(1e-8));

  const vec3 spoke = real(0.2) * kT + real(0.04) * kB + real(0.03) * kN;
  vec3 monge;
  vec3 ps;
  vec3 lined;
  GAUDI_EXPECT(monge_cylinder_target(monge_shape::cylinder_minor, outer_w, kT, kB,
                                    kN, spoke, monge));
  GAUDI_EXPECT(gaudi::albers::nearest_line_spine(outer, spoke, ps));
  GAUDI_EXPECT(gaudi::albers::land_from_spine(ps, outer.r, spoke, lined));
  GAUDI_EXPECT((lined - monge).norm() < real(1e-8));

  const auto huge = monge_torus_of(diag_w(1.0, real(1e-8)), kT, kB, kN);
  GAUDI_EXPECT(huge.valid);
  GAUDI_EXPECT(huge.R > real(1e7));
  const auto kept =
      choose_minor_tube_landing(huge, std::vector<vec3>{spoke});
  GAUDI_EXPECT(kept.fit == minor_tube_fit::torus);
}

/// Matching κ: landings seat. Kind filters isolate cylinder / sphere /
/// ellipsoid when several primitives share a zero residual. Flat major κ
/// keeps a large-R torus (*Points as Tori*).
GAUDI_TEST(spine_shapes_keep_the_outer_torus_and_a_flat_major) {
  using hepworth::block::choose_spine_shapes;
  using hepworth::block::spine_curvature;
  using hepworth::block::spine_kind;
  using hepworth::block::spine_set;

  const real R = 1;
  const real r = real(0.33);
  const real du = real(0.13);
  spine_curvature outer;
  outer.tube = vec3::UnitZ();
  outer.ruling = vec3::UnitY();
  outer.normal = vec3::UnitX();
  outer.k_hi = -real(1) / r;
  outer.k_lo = -real(1) / (R + r);
  const vec3 foot((R + r), 0, 0);
  const vec3 nbr((R + r) * std::cos(du), (R + r) * std::sin(du), 0);
  const real theta = real(0.4);
  const vec3 tube_pt = r * (std::cos(theta) - real(1)) * outer.normal +
                       r * std::sin(theta) * outer.tube;
  // Tube sample seats the osculating centerpoint; major sample is the probe.
  spine_set torus_only;
  torus_only.torus = true;
  const auto landed = choose_spine_shapes(
      outer, std::vector<vec3>{tube_pt, nbr - foot}, torus_only);
  GAUDI_EXPECT(landed.fit == spine_kind::torus);
  GAUDI_ASSERT(landed.ok.size() == 2 && landed.ok[1]);
  const vec3 hit = foot + landed.hit[1];
  const real rad = std::hypot(hit[0], hit[1]);
  const real off = std::abs(std::hypot(rad - R, hit[2]) - r);
  GAUDI_EXPECT(off < real(1e-4));

  spine_curvature flat = outer;
  flat.k_lo = 0;
  const vec3 along = real(0.08) * outer.ruling + real(0.02) * outer.tube;
  const auto flat_land = choose_spine_shapes(flat, std::vector<vec3>{along});
  GAUDI_EXPECT(flat_land.fit == spine_kind::torus);

  const auto foot_pt =
      choose_spine_shapes(outer, std::vector<vec3>{vec3::Zero()});
  GAUDI_EXPECT(foot_pt.ok.size() == 1 && foot_pt.ok[0]);
  GAUDI_EXPECT(foot_pt.hit[0].norm() < real(1e-6));

  spine_set cylinder_only;
  cylinder_only.cylinder = true;
  const auto on_cyl =
      choose_spine_shapes(outer, std::vector<vec3>{tube_pt}, cylinder_only);
  GAUDI_EXPECT(on_cyl.fit == spine_kind::cylinder);
  GAUDI_EXPECT(on_cyl.ok.size() == 1 && on_cyl.ok[0]);
  GAUDI_EXPECT((on_cyl.hit[0] - tube_pt).norm() < real(1e-6));

  const vec3 sphere_pt =
      -r * outer.normal +
      (r / std::sqrt(real(2))) * (outer.ruling + outer.normal);
  spine_set sphere_only;
  sphere_only.sphere = true;
  const auto on_sphere =
      choose_spine_shapes(outer, std::vector<vec3>{sphere_pt}, sphere_only);
  GAUDI_EXPECT(on_sphere.fit == spine_kind::sphere);
  GAUDI_EXPECT(on_sphere.ok.size() == 1 && on_sphere.ok[0]);
  GAUDI_EXPECT((on_sphere.hit[0] - sphere_pt).norm() < real(1e-6));

  const real semi = std::sqrt(r * (R + r));
  const vec3 ell_pt = -r * outer.normal + semi * outer.ruling;
  spine_set ell_only;
  ell_only.ellipsoid = true;
  const auto on_ell =
      choose_spine_shapes(outer, std::vector<vec3>{tube_pt, ell_pt}, ell_only);
  GAUDI_EXPECT(on_ell.fit == spine_kind::ellipsoid);
  GAUDI_ASSERT(on_ell.ok.size() == 2 && on_ell.ok[1]);
  GAUDI_EXPECT((on_ell.hit[1] - ell_pt).norm() < real(1e-6));
}

/// Osculating centerpoint at r0, goal radius r1, new ruling from κ: the
/// centerpoint is kept, the axis is the live ruling, and the spoke is pushed
/// to distance r1 (|Δ| ≈ |r1 − r0|).
GAUDI_TEST(spine_shapes_keep_centerpoint_replace_axis_grow_tube) {
  using hepworth::block::choose_spine_shapes;
  using hepworth::block::spine_curvature;
  using hepworth::block::spine_kind;
  using hepworth::block::spine_set;

  const real r0 = real(0.5);
  const real r1 = real(0.8);
  const real theta = real(0.35);
  spine_curvature geom;
  geom.tube = vec3::UnitZ();
  geom.ruling = vec3::UnitY();
  geom.normal = vec3::UnitX();
  geom.k_hi = -real(1) / r0;
  geom.k_lo = 0;
  const vec3 spoke = r0 * (std::cos(theta) - real(1)) * geom.normal +
                     r0 * std::sin(theta) * geom.tube;
  spine_set cylinder_only;
  cylinder_only.cylinder = true;
  // Rest C from canal fit on the r0 tube; goal κ asks for r1 on project.
  spine_curvature goal = geom;
  goal.k_hi = -real(1) / r1;
  const std::vector<vec3> fit{spoke};
  const auto land =
      choose_spine_shapes(goal, std::vector<vec3>{spoke}, cylinder_only, fit);
  GAUDI_EXPECT(land.fit == spine_kind::cylinder);
  GAUDI_ASSERT(land.ok.size() == 1 && land.ok[0]);
  const vec3 hit = land.hit[0];
  const vec3 canal = -r0 * geom.normal;
  GAUDI_EXPECT(std::abs((hit - canal).norm() - r1) < real(1e-5));
  GAUDI_EXPECT((hit - spoke).norm() > real(0.2));
  GAUDI_EXPECT(std::abs((hit - spoke).norm() - std::abs(r1 - r0)) < real(1e-4));
  GAUDI_EXPECT(std::abs(hit.dot(geom.ruling)) < real(1e-6));

  // Half bite: canal radial moves halfway from r0 toward r1.
  const auto half = choose_spine_shapes(goal, std::vector<vec3>{spoke},
                                        cylinder_only, fit, real(0.5));
  GAUDI_ASSERT(half.ok.size() == 1 && half.ok[0]);
  const real r_half = (half.hit[0] - canal).norm();
  GAUDI_EXPECT(std::abs(r_half - real(0.5) * (r0 + r1)) < real(1e-5));
}

GAUDI_TEST(willmore_curvature_projection_holds_a_sphere) {
  const real H = 2;
  const real area = 0.05;
  const real sum_w = 6;
  const real Z =
      project_willmore_curvature(H, H, H * H, area, sum_w, sum_w * H);
  GAUDI_EXPECT(std::abs(Z - H) < real(1e-8));
}

GAUDI_TEST(willmore_curvature_projection_leaves_a_cylinder) {
  const real H = -2;
  const real area = 0.05;
  const real sum_w = 6;
  const real Z = project_willmore_curvature(H, H, 0, area, sum_w, sum_w * H);
  GAUDI_EXPECT(std::abs(Z) < real(1e-8));
}

GAUDI_TEST(willmore_curvature_projection_prefers_the_umbilic) {
  const real H = 2;
  const real Kcurv = 1;
  const real area = 0.05;
  const real sum_w = 6;
  const real Z =
      project_willmore_curvature(H, H, Kcurv, area, sum_w, sum_w * H);
  GAUDI_EXPECT(std::abs(Z - real(1)) < real(1e-8));
}

/// A cotangent term much larger than H³ has one real root far from H.
/// Taking it is what overflows the next jet.
GAUDI_TEST(willmore_curvature_projection_rejects_a_runaway_root) {
  const real H = 1;
  const real area = real(1e-8);
  const real sum_w = 6;
  const real Z = project_willmore_curvature(H, H, 0, area, sum_w, real(0));
  GAUDI_EXPECT(std::abs(Z - H) < real(1e-8));
}

/// A spoke longer than the principal radius is not a graph landing.
GAUDI_TEST(monge_graph_leaves_a_spoke_outside_the_chart) {
  const shape_w W = diag_w(2.0, 0.5);
  const vec3 far = real(5) * kT;
  const std::vector<objective> hyperboloid{
      monge_fairing{monge_shape::hyperboloid}};
  vec3 hit = vec3::Ones();
  GAUDI_EXPECT(!select_monge_spoke(hyperboloid, far, W, nullptr, kT, kB, kN, hit));
}

/// The stable height lies on the exact osculating torus, on both equators.
GAUDI_TEST(osculating_torus_height_lies_on_the_torus) {
  using hepworth::block::osculating_torus_height;
  for (const real k_lo : {real(0.5), real(-0.5), real(0.05)}) {
    const real k_hi = 2;
    const mat3 W = k_hi * kT * kT.transpose() + k_lo * kB * kB.transpose();
    const albers::osculating_torus tor =
        albers::osculating_torus_from_shape_operator(vec3::Zero(), kN, W);
    GAUDI_ASSERT(tor.valid);
    real worst = 0;
    for (real u = -0.3; u <= 0.3; u += 0.1) {
      for (real v = -0.3; v <= 0.3; v += 0.1) {
        real h = 0;
        GAUDI_ASSERT(osculating_torus_height(k_hi, k_lo, u, v, h));
        const vec3 p = u * kT + v * kB + h * kN;
        worst = std::max(worst, std::abs(albers::torus_sdf(p, tor)));
      }
    }
    GAUDI_EXPECT(worst < real(1e-10));
  }
}

/// Flat and cylindrical limits are exact, with no floor on either curvature.
GAUDI_TEST(osculating_torus_height_flat_limits) {
  using hepworth::block::circle_sag;
  using hepworth::block::osculating_torus_height;
  real h = 1;
  GAUDI_ASSERT(osculating_torus_height(0, 0, 0.2, 0.3, h));
  GAUDI_EXPECT(h == real(0));
  real cyl = 0;
  GAUDI_ASSERT(circle_sag(1.5, 0.2, cyl));
  GAUDI_ASSERT(osculating_torus_height(1.5, 1e-14, 0.2, 0.3, h));
  GAUDI_EXPECT(std::abs(h - cyl) < real(1e-14));
  // Tiny curvatures keep full relative precision: h ≈ ½ k t².
  GAUDI_ASSERT(osculating_torus_height(1e-9, 1e-10, 0.2, 0.3, h));
  const real jet = real(0.5) * (real(1e-9) * 0.04 + real(1e-10) * 0.09);
  GAUDI_EXPECT(std::abs(h - jet) < real(1e-12) * jet);
  GAUDI_EXPECT(!osculating_torus_height(10, 0, 0.2, 0, h));
}

/// Ring points on an offset, tilted copy of the store's torus: the pose fit
/// recovers the offset and tilt, and every point lands on itself.
GAUDI_TEST(osculating_torus_ring_recovers_the_pose) {
  using namespace hepworth::block;
  const torus_ring_shape shape{.k_hi = 2, .k_lo = -0.6, .hx = 1, .hy = 0, .lx = 0, .ly = 1};
  const torus_ring_pose truth{.offset = 0.01, .tilt_t = 0.05, .tilt_b = -0.03};
  const torus_ring_frame f = tilt_frame(kT, kB, kN, truth);
  GAUDI_ASSERT(f.ok);
  const vec3 eh = shape.hx * f.T + shape.hy * f.B;
  const vec3 el = shape.lx * f.T + shape.ly * f.B;
  std::vector<torus_ring_point> pts;
  for (int ring = 1; ring <= 2; ++ring) {
    for (int k = 0; k < 6 * ring; ++k) {
      const real a = real(2 * M_PI) * k / (6 * ring) + real(0.1) * ring;
      const real u = real(0.1) * ring * std::cos(a);
      const real v = real(0.1) * ring * std::sin(a);
      real h = 0;
      GAUDI_ASSERT(osculating_torus_height(shape.k_hi, shape.k_lo, u, v, h));
      pts.push_back({.d = truth.offset * f.N + u * eh + v * el + h * f.N,
                     .s = ring == 1 ? real(1) : real(0.5)});
    }
  }
  std::vector<char> in;
  const osculating_torus_ring_config cfg{.two_ring = 0.5, .iterations = 8, .prior = 1e-14};
  const torus_ring_fit fit = fit_torus_ring_pose(shape, kT, kB, kN, pts, in, cfg);
  GAUDI_ASSERT(fit.ok);
  GAUDI_EXPECT(std::abs(fit.pose.offset - truth.offset) < real(1e-6));
  GAUDI_EXPECT(std::abs(fit.pose.tilt_t - truth.tilt_t) < real(1e-6));
  GAUDI_EXPECT(std::abs(fit.pose.tilt_b - truth.tilt_b) < real(1e-6));
  GAUDI_EXPECT(fit.cost < real(1e-14));
  const torus_ring_frame g = tilt_frame(kT, kB, kN, fit.pose);
  real worst = 0;
  for (const torus_ring_point &pt : pts) {
    vec3 land;
    real gap = 0;
    GAUDI_ASSERT(land_on_posed_torus(shape, g, fit.pose.offset, pt.d, land, gap));
    worst = std::max(worst, (land - pt.d).norm());
  }
  GAUDI_EXPECT(worst < real(1e-6));
}

/// The landing is the closest point: a point off the torus by a known normal
/// distance lands back at its foot, and `gap` is that distance.
GAUDI_TEST(osculating_torus_landing_is_the_closest_point) {
  using namespace hepworth::block;
  const torus_ring_shape shape{.k_hi = 2, .k_lo = 0.5, .hx = 1, .hy = 0, .lx = 0, .ly = 1};
  const torus_ring_frame f = tilt_frame(kT, kB, kN, {});
  const real u = 0.3, v = -0.2;
  real h = 0, hu = 0, hv = 0;
  GAUDI_ASSERT(osculating_torus_height_grad(shape.k_hi, shape.k_lo, u, v, h, hu, hv));
  real hp = 0, hm = 0;
  GAUDI_ASSERT(osculating_torus_height(shape.k_hi, shape.k_lo, u + 1e-6, v, hp));
  GAUDI_ASSERT(osculating_torus_height(shape.k_hi, shape.k_lo, u - 1e-6, v, hm));
  GAUDI_EXPECT(std::abs((hp - hm) / 2e-6 - hu) < real(1e-6));
  GAUDI_ASSERT(osculating_torus_height(shape.k_hi, shape.k_lo, u, v + 1e-6, hp));
  GAUDI_ASSERT(osculating_torus_height(shape.k_hi, shape.k_lo, u, v - 1e-6, hm));
  GAUDI_EXPECT(std::abs((hp - hm) / 2e-6 - hv) < real(1e-6));
  const vec3 foot = u * kT + v * kB + h * kN;
  const vec3 normal = vec3(-hu, -hv, 1).normalized();
  const real off = 0.02;
  vec3 land;
  real gap = 0;
  GAUDI_ASSERT(land_on_posed_torus(shape, f, 0, foot + off * normal, land, gap, /*closest=*/true));
  GAUDI_EXPECT((land - foot).norm() < real(1e-7));
  GAUDI_EXPECT(std::abs(gap - off) < real(1e-7));
}

/// A sphere carried by `S_b` lands on another sphere through the foot, tangent
/// there, so every landing is equidistant from a centre on the normal and the
/// cyclide tensor is umbilic everywhere at one over that radius.
GAUDI_TEST(osculating_dupin_sphere_stays_a_sphere) {
  using namespace hepworth::block;
  const real k = 1.5;
  const dupin_ring_shape shape{.torus = {.k_hi = k, .k_lo = k}, .b_hi = 0.2, .b_lo = -0.1};
  std::vector<vec3> pts;
  for (real u = -0.3; u <= 0.31; u += 0.15)
    for (real v = -0.3; v <= 0.31; v += 0.15) {
      vec3 land;
      real gap = 0;
      GAUDI_ASSERT(dupin_land_local(shape, vec3(u, v, 0), land, gap));
      pts.push_back(land);
    }
  // A sphere through the origin tangent to z = 0 has its centre at c e_z with
  // c = |p|² / (2 p_z) for any of its points.
  const vec3 far = pts.front();
  const real c = far.squaredNorm() / (real(2) * far[2]);
  real worst = 0;
  for (const vec3 &x : pts)
    worst = std::max(worst, std::abs((x - c * kN).norm() - std::abs(c)));
  GAUDI_EXPECT(worst < real(1e-9));
  for (const vec3 &x : pts) {
    mat3 W;
    vec3 n;
    GAUDI_ASSERT(dupin_tensor_local(shape, x, real(1e-4), W, n));
    vec3 t = n.unitOrthogonal();
    const vec3 bt = n.cross(t);
    GAUDI_EXPECT(std::abs(t.dot(W * t) - real(1) / c) < real(1e-5));
    GAUDI_EXPECT(std::abs(bt.dot(W * bt) - real(1) / c) < real(1e-5));
    GAUDI_EXPECT(std::abs(t.dot(W * bt)) < real(1e-5));
  }
}

/// To first order in `b` the cyclide gains exactly the two cross cubic terms
/// over its torus, and no ridge terms.
GAUDI_TEST(osculating_dupin_adds_the_cross_cubic_terms) {
  using namespace hepworth::block;
  const real kh = 3, kl = 0.5, bh = 2e-3, bl = -1.5e-3;
  const dupin_ring_shape shape{.torus = {.k_hi = kh, .k_lo = kl}, .b_hi = bh, .b_lo = bl};
  real worst = 0, scale = 0;
  for (real u = -0.05; u <= 0.051; u += 0.025)
    for (real v = -0.05; v <= 0.051; v += 0.025) {
      // A cyclide point from a torus point, then its height over the torus.
      real h = 0, hu = 0, hv = 0;
      GAUDI_ASSERT(osculating_torus_height_grad(kh, kl, u, v, h, hu, hv));
      const vec3 x = special_conformal(shape.b(), vec3(u, v, h));
      real ht = 0;
      GAUDI_ASSERT(osculating_torus_height(kh, kl, x[0], x[1], ht));
      const real cubic = (kl - kh) * bh * x[0] * x[1] * x[1] + (kh - kl) * bl * x[0] * x[0] * x[1];
      worst = std::max(worst, std::abs((x[2] - ht) - cubic));
      scale = std::max(scale, std::abs(cubic));
    }
  GAUDI_EXPECT(scale > real(1e-8));
  GAUDI_EXPECT(worst < real(0.03) * scale);
}

/// Ring points and star edge tensors sampled on a posed cyclide: starting from
/// its plain torus and a slightly wrong store, the shape fit recovers both
/// curvatures, `b` and the pose, and every point lands on itself.
GAUDI_TEST(osculating_dupin_ring_recovers_shape_and_pose) {
  using namespace hepworth::block;
  const dupin_ring_shape truth{.torus = {.k_hi = 2.5, .k_lo = 0.4}, .b_hi = 0.15, .b_lo = -0.2};
  const torus_ring_pose pose{.offset = 0.004, .tilt_t = 0.03, .tilt_b = -0.02};
  const torus_ring_frame f = tilt_frame(kT, kB, kN, pose);
  GAUDI_ASSERT(f.ok);
  auto on_cyclide = [&](real u, real v) {
    real h = 0, hu = 0, hv = 0;
    osculating_torus_height_grad(truth.torus.k_hi, truth.torus.k_lo, u, v, h, hu, hv);
    const vec3 l = special_conformal(truth.b(), vec3(u, v, h));
    return pose.offset * f.N + l[0] * f.T + l[1] * f.B + l[2] * f.N;
  };
  std::vector<torus_ring_point> pts;
  std::vector<vec3> ring1;
  for (int ring = 1; ring <= 2; ++ring)
    for (int k = 0; k < 6 * ring; ++k) {
      const real a = real(2 * M_PI) * k / (6 * ring) + real(0.2) * ring;
      const vec3 d = on_cyclide(real(0.12) * ring * std::cos(a), real(0.12) * ring * std::sin(a));
      pts.push_back({.d = d, .s = ring == 1 ? real(1) : real(0.5)});
      if (ring == 1)
        ring1.push_back(d);
    }
  // Star edges from the foot, which sits off the cyclide by the pose offset.
  std::vector<dupin_ring_edge> edges;
  for (const vec3 &d : ring1) {
    dupin_ring_edge e{.lo = vec3::Zero(), .hi = d, .mass = 0.01};
    GAUDI_ASSERT(dupin_edge_reading(truth, f, pose.offset, e, real(1e-5), e.store));
    edges.push_back(e);
  }
  dupin_ring_shape start = truth;
  start.torus.k_hi = 2.3;
  start.torus.k_lo = 0.5;
  start.b_hi = 0;
  start.b_lo = 0;
  std::vector<char> in;
  // Exact recovery checks the fit itself, so the ridge on b is off.
  const osculating_dupin_ring_config cfg{.iterations = 20, .prior = 1e-14, .b_ridge = 0};
  const dupin_ring_fit fit = fit_dupin_ring(start, kT, kB, kN, pts, edges, in, cfg, true);
  GAUDI_ASSERT(fit.ok);
  GAUDI_EXPECT(std::abs(fit.shape.torus.k_hi - truth.torus.k_hi) < real(1e-4));
  GAUDI_EXPECT(std::abs(fit.shape.torus.k_lo - truth.torus.k_lo) < real(1e-4));
  GAUDI_EXPECT(std::abs(fit.shape.b_hi - truth.b_hi) < real(1e-3));
  GAUDI_EXPECT(std::abs(fit.shape.b_lo - truth.b_lo) < real(1e-3));
  GAUDI_EXPECT(std::abs(fit.pose.offset - pose.offset) < real(1e-6));
  GAUDI_EXPECT(std::abs(fit.pose.tilt_t - pose.tilt_t) < real(1e-5));
  GAUDI_EXPECT(std::abs(fit.pose.tilt_b - pose.tilt_b) < real(1e-5));
  const torus_ring_frame g = tilt_frame(kT, kB, kN, fit.pose);
  real worst = 0;
  for (const torus_ring_point &pt : pts) {
    vec3 land;
    real gap = 0;
    GAUDI_ASSERT(land_on_posed_dupin(fit.shape, g, fit.pose.offset, pt.d, land, gap));
    worst = std::max(worst, (land - pt.d).norm());
  }
  GAUDI_EXPECT(worst < real(1e-6));
}

/// The sympy-generated gradient of a ring point's gap matches central
/// differences of `land_on_posed_dupin` in all seven fit parameters.
GAUDI_TEST(osculating_dupin_gap_gradient_matches_differences) {
  using namespace hepworth::block;
  const real a = 0.7;
  const vec3 T = vec3(std::cos(0.2), std::sin(0.2), 0), B = vec3(-std::sin(0.2), std::cos(0.2), 0), N = kN;
  const std::array<real, 7> x0{0.006, 0.04, -0.025, 2.4, -0.7, 0.11, -0.16};
  auto unpack = [&](const std::array<real, 7> &x, torus_ring_pose &pose, dupin_ring_shape &shape) {
    pose = {.offset = x[0], .tilt_t = x[1], .tilt_b = x[2]};
    shape = {.torus = {.k_hi = x[3], .k_lo = x[4], .hx = std::cos(a), .hy = std::sin(a),
                       .lx = -std::sin(a), .ly = std::cos(a)},
             .b_hi = x[5], .b_lo = x[6]};
  };
  auto gap_at = [&](const std::array<real, 7> &x, const vec3 &d) {
    torus_ring_pose pose;
    dupin_ring_shape shape;
    unpack(x, pose, shape);
    vec3 land;
    real gap = 0;
    GAUDI_ASSERT(land_on_posed_dupin(shape, tilt_frame(T, B, N, pose), pose.offset, d, land, gap));
    return gap;
  };
  real worst = 0;
  for (const vec3 &d : {vec3(0.1, 0.05, 0.02), vec3(-0.12, 0.08, -0.01), vec3(0.03, -0.14, 0.03)}) {
    torus_ring_pose pose;
    dupin_ring_shape shape;
    unpack(x0, pose, shape);
    real gap = 0;
    Eigen::Matrix<real, 7, 1> grad;
    dupin_gap_and_gradient(shape, T, B, N, pose, d, gap, grad);
    GAUDI_EXPECT(std::abs(gap - gap_at(x0, d)) < real(1e-12));
    for (int k = 0; k < 7; ++k) {
      std::array<real, 7> xp = x0, xm = x0;
      xp[k] += 1e-6;
      xm[k] -= 1e-6;
      const real fd = (gap_at(xp, d) - gap_at(xm, d)) / 2e-6;
      worst = std::max(worst, std::abs(fd - grad[k]));
    }
  }
  GAUDI_EXPECT(worst < real(1e-7));
}

/// A later store reading whose tube role swapped and whose eigenvectors came
/// back flipped is matched to the frame-start axes, so the reference
/// curvatures and `b` keep their meaning through the step.
GAUDI_TEST(osculating_dupin_axes_follow_the_reference) {
  using namespace hepworth::block;
  const real c = std::cos(0.3), s = std::sin(0.3);
  const torus_ring_shape ref{.k_hi = 2.0, .k_lo = -1.8, .hx = c, .hy = s, .lx = -s, .ly = c};
  // The saddle grew so |k_lo| > |k_hi|: the reading calls -2.1 the tube and
  // returns both directions negated.
  monge_axes later;
  later.k_hi = -2.1;
  later.k_lo = 2.0;
  later.hx = s;
  later.hy = -c;
  later.lx = -c;
  later.ly = -s;
  const monge_axes ax = align_axes_to_reference(later, ref);
  GAUDI_EXPECT(ax.k_hi == real(2.0));
  GAUDI_EXPECT(ax.k_lo == real(-2.1));
  GAUDI_EXPECT(std::abs(ax.hx - c) < real(1e-12) && std::abs(ax.hy - s) < real(1e-12));
  GAUDI_EXPECT(std::abs(ax.lx + s) < real(1e-12) && std::abs(ax.ly - c) < real(1e-12));
}

/// The one-ring cone: a flat hexagon of area A under an apex at height a has
/// volume A a / 3, and the gradient matches central differences.
GAUDI_TEST(ring_cone_volume_and_gradient) {
  using namespace hepworth::block;
  std::vector<vec3> d;
  for (int k = 0; k < 6; ++k) {
    const real t = real(2 * M_PI) * k / 6;
    d.push_back(vec3(std::cos(t), std::sin(t), 0) * real(0.1));
  }
  const real area = real(3) * std::sqrt(real(3)) / real(2) * real(0.01);
  const vec3 apex(0, 0, -0.25);
  GAUDI_EXPECT(std::abs(std::abs(fan_cone_volume(d, 6, apex)) - area * real(0.25) / real(3)) < real(1e-14));
  d[2] += vec3(0.01, -0.02, 0.015);
  std::vector<vec3> g;
  fan_cone_gradient(d, 6, apex, g);
  real worst = 0;
  for (size_t j = 0; j < 6; ++j)
    for (int c = 0; c < 3; ++c) {
      std::vector<vec3> dp = d, dm = d;
      dp[j][c] += 1e-6;
      dm[j][c] -= 1e-6;
      const real fd = (fan_cone_volume(dp, 6, apex) - fan_cone_volume(dm, 6, apex)) / 2e-6;
      worst = std::max(worst, std::abs(fd - g[j][c]));
    }
  GAUDI_EXPECT(worst < real(1e-10));
}

/// Hard mode restores the cone volume of a shrunken ring exactly. Soft mode at
/// weight 1 restores half of it. The 2-ring entries past `n` are untouched.
GAUDI_TEST(ring_cone_volume_hard_and_soft) {
  using namespace hepworth::block;
  std::vector<vec3> d;
  for (int k = 0; k < 6; ++k) {
    const real t = real(2 * M_PI) * k / 6;
    d.push_back(vec3(std::cos(t), std::sin(t), real(0.02) * std::cos(2 * t)) * real(0.1));
  }
  d.push_back(vec3(0.3, 0.1, 0.05));
  const std::vector<real> scales{1, 1, 1, 1, 1, 1, 0.5};
  const vec3 apex(0, 0, -0.25);
  const real v0 = fan_cone_volume(d, 6, apex);
  std::vector<vec3> shrunk = d;
  for (size_t j = 0; j < 6; ++j)
    shrunk[j] = real(0.9) * d[j] + vec3(0, 0, -0.004);
  const real v1 = fan_cone_volume(shrunk, 6, apex);
  GAUDI_ASSERT(std::abs(v1 - v0) > real(1e-2) * std::abs(v0));

  std::vector<vec3> hard = shrunk;
  hold_fan_cone_volume(hard, 6, scales, apex, v0, {.mode = ring_volume_mode::hard});
  GAUDI_EXPECT(std::abs(fan_cone_volume(hard, 6, apex) - v0) < real(1e-12) * std::abs(v0));
  GAUDI_EXPECT(hard[6] == shrunk[6]);

  // With the foot free too, the ring and foot move together: the cone over the
  // ring at delta + t with its foot at delta has the volume back.
  std::vector<vec3> with_foot = shrunk;
  const vec3 delta = hold_fan_cone_volume(with_foot, 6, scales, apex, v0, {.mode = ring_volume_mode::hard}, 1.0);
  GAUDI_EXPECT(delta.norm() > real(0));
  std::vector<vec3> moved(6);
  for (size_t j = 0; j < 6; ++j)
    moved[j] = delta + with_foot[j];
  GAUDI_EXPECT(std::abs(fan_cone_volume(moved, 6, apex, delta) - v0) < real(1e-12) * std::abs(v0));
  // A shrunken patch moves away from the apex below it.
  GAUDI_EXPECT(delta[2] > real(0));

  std::vector<vec3> soft = shrunk;
  hold_fan_cone_volume(soft, 6, scales, apex, v0, {.mode = ring_volume_mode::soft, .weight = 1});
  GAUDI_EXPECT(std::abs(fan_cone_volume(soft, 6, apex) - real(0.5) * (v0 + v1)) < real(1e-12) * std::abs(v0));
}

/// A flat foot puts the apex at the cap, on the side fixed at the first call.
GAUDI_TEST(ring_spine_point_caps_the_flat_limit) {
  using namespace hepworth::block;
  GAUDI_EXPECT((ring_spine_point(4.0, 0.0, kN, 1, 2.0) - real(0.25) * kN).norm() < real(1e-15));
  GAUDI_EXPECT((ring_spine_point(0.0, 0.0, kN, -1, 2.0) + real(2.0) * kN).norm() < real(1e-15));
  GAUDI_EXPECT((ring_spine_point(1e-9, 0.01, kN, 1, 2.0) - real(2.01) * kN).norm() < real(1e-12));
}

/// `land = false` only means "volume only" when a volume mode is on. With the
/// mode off the row still lands, so a leftover `land = false` cannot silently
/// switch the shape projection off.
GAUDI_TEST(ring_volume_land_needs_a_volume_mode) {
  using namespace hepworth::block;
  GAUDI_EXPECT((ring_volume_config{.mode = ring_volume_mode::off, .land = false}).lands());
  GAUDI_EXPECT(!(ring_volume_config{.mode = ring_volume_mode::hard, .land = false}).lands());
  GAUDI_EXPECT((ring_volume_config{.mode = ring_volume_mode::soft}).lands());
}

/// The closed-form frame-start shape: from an exact cyclide's star of edge
/// tensors, the slopes of the store across the ring give back `b` to first
/// order, on a convex patch and on a saddle.
GAUDI_TEST(osculating_dupin_direct_shape_recovers_b) {
  using namespace hepworth::block;
  for (const real kl : {real(0.4), real(-0.9)}) {
    const dupin_ring_shape truth{.torus = {.k_hi = 2.5, .k_lo = kl}, .b_hi = 0.12, .b_lo = -0.15};
    const torus_ring_frame f = tilt_frame(kT, kB, kN, {});
    std::vector<dupin_ring_edge> edges;
    for (int k = 0; k < 6; ++k) {
      const real a = real(2 * M_PI) * k / 6 + real(0.2);
      real h = 0, hu = 0, hv = 0;
      const real u = real(0.08) * std::cos(a), v = real(0.08) * std::sin(a);
      osculating_torus_height_grad(truth.torus.k_hi, truth.torus.k_lo, u, v, h, hu, hv);
      const vec3 d = special_conformal(truth.b(), vec3(u, v, h));
      dupin_ring_edge e{.lo = vec3::Zero(), .hi = d, .mass = 0.01};
      GAUDI_ASSERT(dupin_edge_reading(truth, f, 0, e, real(1e-5), e.store));
      edges.push_back(e);
    }
    dupin_ring_shape store = truth;
    store.b_hi = 0;
    store.b_lo = 0;
    const dupin_ring_shape got = direct_dupin_shape(store, kT, kB, kN, edges);
    GAUDI_EXPECT(std::abs(got.b_hi - truth.b_hi) < real(0.15) * std::abs(truth.b_hi));
    GAUDI_EXPECT(std::abs(got.b_lo - truth.b_lo) < real(0.15) * std::abs(truth.b_lo));
  }
}

/// The ridge on b: on a near-umbilic star with a cross slope that no Dupin
/// cyclide can carry there, the closed form keeps `b` small instead of
/// dividing the slope by the vanishing curvature gap.
GAUDI_TEST(osculating_dupin_ridge_tames_an_umbilic) {
  using namespace hepworth::block;
  const dupin_ring_shape store{.torus = {.k_hi = 2.0, .k_lo = 1.99}};
  std::vector<dupin_ring_edge> edges;
  for (int k = 0; k < 6; ++k) {
    const real a = real(2 * M_PI) * k / 6;
    const vec3 d(real(0.1) * std::cos(a), real(0.1) * std::sin(a), real(0.01));
    // An isotropic tensor tilted by a cross slope of 1 along u.
    const real m = real(0.5) * d[0];
    dupin_ring_edge e{.lo = vec3::Zero(), .hi = d, .mass = 0.01};
    e.store = vec3(real(2.0), real(0), real(1.99) + m);
    edges.push_back(e);
  }
  const dupin_ring_shape loose = direct_dupin_shape(store, kT, kB, kN, edges, real(0), real(1e9));
  const dupin_ring_shape tamed = direct_dupin_shape(store, kT, kB, kN, edges);
  const real reach = 0.05;
  // Without the ridge the slope is divided by the 0.01 gap; with it, b
  // collapses toward the plain torus, as a cyclide must at an umbilic.
  GAUDI_EXPECT(std::abs(loose.b_hi) > real(20) * std::abs(tamed.b_hi));
  GAUDI_EXPECT(std::abs(tamed.b_hi) * reach < real(0.01));
}

namespace {
/// A torus-graph model in the chart (T, B, N) = (x, y, z), principal axes on
/// x and y, for the variational fit tests.
struct test_torus_model {
  hepworth::block::ring_model_point land(const hepworth::block::ring_fit_params &x, const vec3 &d) const {
    using namespace hepworth::block;
    ring_model_point e;
    const torus_ring_frame fr = tilt_frame(kT, kB, kN, {.offset = x[0], .tilt_t = x[1], .tilt_b = x[2]});
    if (!fr.ok)
      return e;
    const vec3 p = d - x[0] * fr.N;
    const real u = p.dot(fr.T), v = p.dot(fr.B);
    real h = 0, hu = 0, hv = 0;
    if (!osculating_torus_height_grad(x[3], x[4], u, v, h, hu, hv))
      return e;
    e.gap = p.dot(fr.N) - h;
    e.land = x[0] * fr.N + u * fr.T + v * fr.B + h * fr.N;
    e.normal = (fr.N - hu * fr.T - hv * fr.B).normalized();
    e.ok = true;
    return e;
  }
  void foot(const hepworth::block::ring_fit_params &x, vec3 &F, vec3 &NF) const {
    const auto fr = hepworth::block::tilt_frame(kT, kB, kN, {.offset = x[0], .tilt_t = x[1], .tilt_b = x[2]});
    NF = fr.N;
    F = x[0] * NF;
  }
};

std::vector<vec3> torus_ring(real kh, real kl, real reach, int n) {
  std::vector<vec3> d;
  for (int k = 0; k < n; ++k) {
    const real a = real(2 * M_PI) * k / n + real(0.1);
    const real u = reach * std::cos(a), v = reach * std::sin(a);
    real h = 0;
    hepworth::block::osculating_torus_height(kh, kl, u, v, h);
    d.push_back(vec3(u, v, h));
  }
  return d;
}
} // namespace

/// The medial prism: on a sphere every curvature-sphere centre coincides, so
/// it is the cone to the centre exactly; on a cylinder the centres lie on the
/// axis and it is the wedge, about half the fan area times the radius.
GAUDI_TEST(ring_medial_prism_is_cone_on_sphere_and_wedge_on_cylinder) {
  using namespace hepworth::block;
  const test_torus_model model;
  auto land = [&](const ring_fit_params &x, const vec3 &d) { return model.land(x, d); };
  auto foot = [&](const ring_fit_params &x, vec3 &F, vec3 &NF) { model.foot(x, F, NF); };
  {
    const real k = 4;
    const std::vector<vec3> d = torus_ring(k, k, 0.06, 6);
    ring_fit_params x;
    x << 0, 0, 0, k, k;
    real v = 0, a = 0;
    GAUDI_ASSERT(ring_medial_volume(x, d, 6, land, foot, 10.0, v, a));
    std::vector<vec3> L(6);
    for (int j = 0; j < 6; ++j)
      L[j] = model.land(x, d[j]).land;
    // Same tetrahedra, listed foot-first here and apex-first there: an odd
    // permutation, so equal and opposite.
    GAUDI_EXPECT(std::abs(v + fan_cone_volume(L, 6, vec3(0, 0, 1 / k))) < real(1e-12));
  }
  {
    const real k = 4;
    const std::vector<vec3> d = torus_ring(k, 0, 0.05, 6);
    ring_fit_params x;
    x << 0, 0, 0, k, 0;
    real v = 0, a = 0;
    GAUDI_ASSERT(ring_medial_volume(x, d, 6, land, foot, 10.0, v, a));
    GAUDI_EXPECT(std::abs(std::abs(v) - real(0.5) * a / k) < real(0.02) * a / k);
  }
}

/// The variational fit on an exact torus ring: with its own volume as the
/// reference it stays on the torus; asked for ten percent more volume with
/// the volume term dominant, it meets it by fattening the tube.
GAUDI_TEST(ring_variational_fit_meets_volume_by_shape) {
  using namespace hepworth::block;
  const test_torus_model model;
  auto land = [&](const ring_fit_params &x, const vec3 &d) { return model.land(x, d); };
  auto foot = [&](const ring_fit_params &x, vec3 &F, vec3 &NF) { model.foot(x, F, NF); };
  const real kh = 4, kl = 1;
  const std::vector<vec3> d = torus_ring(kh, kl, 0.06, 6);
  const std::vector<real> s(6, 1);
  ring_fit_params start;
  start << 0, 0, 0, kh, kl;
  real v0 = 0, a0 = 0;
  GAUDI_ASSERT(ring_medial_volume(start, d, 6, land, foot, 10.0, v0, a0));

  const ring_variational_config cfg{.w_store = 10, .w_volume = 1, .prior = 1e-12, .iterations = 10};
  const ring_variational_fit same = fit_ring_variational(start, kh, kl, d, s, 6, v0, cfg, land, foot);
  GAUDI_ASSERT(same.ok);
  GAUDI_EXPECT((same.theta - start).norm() < real(1e-6));

  const ring_variational_config hard{.w_store = 1, .w_volume = 1e6, .prior = 1e-12, .iterations = 20};
  const ring_variational_fit more =
      fit_ring_variational(start, kh, kl, d, s, 6, real(1.1) * v0, hard, land, foot);
  GAUDI_ASSERT(more.ok);
  GAUDI_EXPECT(std::abs(more.volume - real(1.1) * v0) < real(1e-3) * std::abs(v0));
  GAUDI_EXPECT(more.theta[3] < kh);
}

} // namespace test
} // namespace gaudi

#endif
