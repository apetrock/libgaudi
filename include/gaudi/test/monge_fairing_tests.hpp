#ifndef GAUDI_MONGE_FAIRING_TESTS_HPP
#define GAUDI_MONGE_FAIRING_TESTS_HPP

#include "gaudi/common.h"
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

} // namespace test
} // namespace gaudi

#endif
