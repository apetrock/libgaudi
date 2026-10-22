#ifndef GAUDI_SHAPE_BLOCK_TESTS_HPP
#define GAUDI_SHAPE_BLOCK_TESTS_HPP

#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/dynamic.hpp"
#include "gaudi/hepworth/block/face_element.hpp"
#include "gaudi/hepworth/block/shape_projection.hpp"
#include "gaudi/hepworth/block/shell_constraints_init.hpp"
#include "gaudi/hepworth/block/solver.hpp"
#include "gaudi/hepworth/block/variational_audit.hpp"
#include "gaudi/hepworth/blocks/shell_cr_block.hpp"
#include "gaudi/hepworth/blocks/shell_shape_block.hpp"
#include "gaudi/test/test.hpp"

#include <Eigen/Geometry>

#include <algorithm>
#include <cmath>
#include <random>
#include <utility>

namespace gaudi {
namespace test {

struct analytic_torus {
  real k_tube = 0;
  real k_major = 0;
  vec3 tube = vec3::Zero();
  vec3 major = vec3::Zero();
};

inline analytic_torus torus_at(real R, real r, const vec3 &c) {
  const real u = std::atan2(c.y(), c.x());
  const real rho = std::hypot(c.x(), c.y());
  const real v = std::atan2(c.z(), rho - R);
  analytic_torus a;
  a.k_tube = real(1) / r;
  a.k_major = std::cos(v) / (R + r * std::cos(v));
  a.major = vec3(-std::sin(u), std::cos(u), 0).normalized();
  a.tube = vec3(-std::cos(u) * std::sin(v), -std::sin(u) * std::sin(v),
                std::cos(v))
               .normalized();
  return a;
}


inline vec3 project_to_plane(vec3 dir, const vec3 &N) {
  dir -= N * dir.dot(N);
  return dir.normalized();
}

inline real axis_degrees(const vec3 &axis, const vec3 &dir) {
  const real c = std::min(real(1), std::abs(axis.normalized().dot(dir)));
  return std::acos(c) * real(180) / std::acos(real(-1));
}

inline real frame_curvature(const hepworth::block::projected_face_frame &fr,
                            const vec3 &dir) {
  const real c1 = dir.dot(fr.e1);
  const real c2 = dir.dot(fr.e2);
  return fr.k1 * c1 * c1 + fr.k2 * c2 * c2;
}

GAUDI_TEST(frame_angle_turns_the_axes_and_keeps_the_radii) {
  const real half = std::acos(real(-1)) * real(0.5);
  const vec3 turned =
      hepworth::block::rotate_frame_tangent(vec3::UnitX(), vec3::UnitZ(), half);
  GAUDI_EXPECT((turned - vec3::UnitY()).norm() < real(1e-8));
  const vec3 held =
      hepworth::block::rotate_frame_tangent(vec3::UnitX(), vec3::UnitZ(), real(0));
  GAUDI_EXPECT((held - vec3::UnitX()).norm() < real(1e-12));
  const vec3 coeffs(real(3), real(0), real(1));
  const hepworth::block::shape_w W =
      hepworth::block::shape_w_from_coeffs(coeffs);
  GAUDI_EXPECT(W.ok);
  GAUDI_EXPECT(std::abs(W.k1 - real(3)) < real(1e-12));
  GAUDI_EXPECT(std::abs(W.k2 - real(1)) < real(1e-12));
  real gp = real(1);
  real gq = real(0);
  hepworth::block::rotate_edge_pair(gp, gq, half);
  GAUDI_EXPECT(std::abs(gp) < real(1e-8));
  GAUDI_EXPECT(std::abs(gq - real(1)) < real(1e-8));
  real gp0 = real(1);
  real gq0 = real(0);
  hepworth::block::rotate_edge_pair(gp0, gq0, real(0));
  GAUDI_EXPECT(std::abs(gp0 - real(1)) < real(1e-12));
  GAUDI_EXPECT(std::abs(gq0) < real(1e-12));
}

GAUDI_TEST(stretch_axis_scale_follows_curvature_along_the_push) {
  hepworth::block::shape_w W;
  W.k11 = real(4);
  W.k12 = real(0);
  W.k22 = real(0.25);
  W.ok = true;
  const real held = hepworth::block::stretch_axis_scale(
      W, real(1), real(0), real(0.2), real(0));
  GAUDI_EXPECT(std::abs(held - real(0.2)) < real(1e-12));
  const real flat_power = hepworth::block::stretch_axis_scale(
      W, real(0), real(0), real(0.2), real(0));
  GAUDI_EXPECT(std::abs(flat_power - real(0.2)) < real(1e-12));
  const real tight = hepworth::block::stretch_axis_scale(
      W, real(1), real(0), real(0.2), real(-1));
  GAUDI_EXPECT(std::abs(tight - real(0.2) / real(4)) < real(1e-12));
  const real loose = hepworth::block::stretch_axis_scale(
      W, real(0), real(1), real(0.2), real(-1));
  GAUDI_EXPECT(std::abs(loose - real(0.2) / real(0.25)) < real(1e-12));
  const real s = std::sqrt(real(0.5));
  const real tilt = hepworth::block::stretch_axis_scale(
      W, s, s, real(0.2), real(-1));
  GAUDI_EXPECT(std::abs(tilt - real(0.2) / real(2.125)) < real(1e-12));
  const real negative = hepworth::block::stretch_axis_scale(
      W, real(-1), real(0), real(0.2), real(-1));
  GAUDI_EXPECT(std::abs(negative - tight) < real(1e-12));
}

GAUDI_TEST(directed_stretch_lengthens_the_chart_axis_only) {
  const vec3 axis = vec3::UnitX();
  const vec3 along = hepworth::block::scale_spoke(vec3(1, 0, 0), axis, real(1));
  const vec3 across = hepworth::block::scale_spoke(vec3(0, 1, 0), axis, real(1));
  const vec3 flipped = hepworth::block::scale_spoke(vec3(1, 0, 0), -axis, real(1));
  const vec3 both = hepworth::block::scale_spoke(
      hepworth::block::scale_spoke(vec3(2, 4, 1), vec3::UnitX(), real(0.5)),
      vec3::UnitY(), real(-0.25));
  GAUDI_EXPECT((along - vec3(2, 0, 0)).norm() < real(1e-12));
  GAUDI_EXPECT((across - vec3(0, 1, 0)).norm() < real(1e-12));
  GAUDI_EXPECT((flipped - along).norm() < real(1e-12));
  GAUDI_EXPECT((both - vec3(3, 3, 1)).norm() < real(1e-12));
  vec3 u, v;
  hepworth::block::stretch_frame(vec3::UnitX(), vec3::UnitY(),
                                 real(0.5) * std::acos(real(-1)), u, v);
  const vec3 turned = hepworth::block::scale_spoke(
      hepworth::block::scale_spoke(vec3(2, 4, 1), u, real(0.5)), v,
      real(-0.25));
  GAUDI_EXPECT((turned - vec3(1.5, 6, 1)).norm() < real(1e-12));

  const std::vector<vec3> rest_x{vec3::Zero(), vec3(2, 0, 0), vec3(0, 4, 0)};
  std::vector<vec3> live_x{vec3::Zero(), vec3(3, 0, 0), vec3(0, 3, 0)};
  std::vector<vec3> mass(3, vec3::Ones());
  std::vector<vec3> vel(3, vec3::Zero());
  std::vector<vec3> force(3, vec3::Zero());
  auto pos = hepworth::vec3_block::create(mass, live_x, vel, force);
  std::vector<vec3> shape_x{vec3(real(3), real(0), real(1)), vec3::Zero(),
                            vec3::Zero()};
  auto shape = hepworth::vec3_block::create(mass, shape_x, vel, force);
  std::vector<real> cu{real(1), real(0)};
  std::vector<real> cv{real(0), real(0)};
  std::vector<real> cm{real(1), real(1)};
  auto cr = hepworth::block::cr_coeff_block::create(cu, cv, cm, 1);
  std::vector<hepworth::sim_block::ptr> blocks{pos, cr, shape};
  hepworth::block::fan_local fan;
  fan.vids = {0, 1, 2};
  fan.center = 0;
  hepworth::block::cr_edge edge;
  edge.v0 = 0;
  edge.v1 = 1;
  fan.ledges.push_back(edge);
  fan.dof_of = {0};
  fan.n_incident = 1;
  hepworth::block::cr_tri tri;
  tri.v = {0, 1, 2};
  tri.e = {0, -1, -1};
  fan.tris.push_back(tri);
  fan.tangent_e = {0};
  std::vector<index_t> ids{0, 1, 2};
  std::vector<real> weights{real(1), real(1)};
  std::vector<hepworth::block::vertex_chart> charts(1);
  charts[0].ok = true;
  charts[0].T = vec3::UnitX();
  charts[0].B = vec3::UnitY();
  charts[0].N = vec3::UnitZ();
  charts[0].h = vec3(real(3), real(0), real(1));
  auto row = hepworth::block::directed_stretch::create(
      ids, weights, rest_x, fan, 1, real(0), real(0.5), real(-0.25), real(0),
      real(1), blocks, &charts);
  vecX q;
  for (const auto &b : blocks)
    b->map_to_x(q);
  index_t id0 = 0;
  std::vector<trip> trips;
  row->fill_A(id0, trips);
  GAUDI_EXPECT(id0 == static_cast<index_t>(6));
  for (const trip &t : trips) {
    GAUDI_EXPECT(t.col() >= 0);
    GAUDI_EXPECT(t.col() < q.size());
  }
  vecX p = vecX::Zero(id0);
  row->project(q, p);
  GAUDI_EXPECT(p.allFinite());
  GAUDI_EXPECT((row->rest()[0] - vec3(3, 0, 0)).norm() < real(1e-8));
  GAUDI_EXPECT((row->rest()[1] - vec3(0, 3, 0)).norm() < real(1e-8));
  GAUDI_EXPECT((vec3(p[0], p[1], p[2]) - vec3(3, 0, 0)).norm() < real(1e-6));
  GAUDI_EXPECT((vec3(p[3], p[4], p[5]) - vec3(0, 3, 0)).norm() < real(1e-6));

  const real half = real(0.5) * std::acos(real(-1));
  auto angled = hepworth::block::directed_stretch::create(
      ids, weights, rest_x, fan, 1, half, real(0.5), real(-0.25), real(0),
      real(1), blocks, &charts);
  index_t id1 = 0;
  std::vector<trip> trips1;
  angled->fill_A(id1, trips1);
  vecX p1 = vecX::Zero(id1);
  angled->project(q, p1);
  GAUDI_EXPECT((angled->rest()[0] - vec3(real(1.5), 0, 0)).norm() < real(1e-6));
  GAUDI_EXPECT((angled->rest()[1] - vec3(0, 6, 0)).norm() < real(1e-6));
  GAUDI_EXPECT((vec3(p1[0], p1[1], p1[2]) - vec3(real(1.5), 0, 0)).norm() < real(1e-6));
  GAUDI_EXPECT((vec3(p1[3], p1[4], p1[5]) - vec3(0, 6, 0)).norm() < real(1e-6));
}

GAUDI_TEST(map_ring_log_skew_holds_product_and_adds_gap) {
  const hepworth::block::shape_w pair =
      hepworth::block::shape_w_from_coeffs(vec3(4, 0, 1));
  const hepworth::block::shape_w shifted = hepworth::block::map_ring_log_skew(
      pair, {}, std::log(real(4)));
  GAUDI_ASSERT(shifted.ok);
  GAUDI_EXPECT(std::abs(shifted.k1 - real(8)) < real(1e-6));
  GAUDI_EXPECT(std::abs(shifted.k2 - real(0.5)) < real(1e-6));
  GAUDI_EXPECT(std::abs(shifted.k1 * shifted.k2 - real(4)) < real(1e-6));

  const hepworth::block::shape_w saddle =
      hepworth::block::shape_from_principals(real(4), real(-1), real(1), real(0),
                                             real(0), real(1));
  const hepworth::block::shape_w saddle_shift =
      hepworth::block::map_ring_log_skew(saddle, {}, std::log(real(4)));
  GAUDI_ASSERT(saddle_shift.ok);
  GAUDI_EXPECT(std::abs(saddle_shift.k1 - real(8)) < real(1e-6));
  GAUDI_EXPECT(std::abs(saddle_shift.k2 + real(0.5)) < real(1e-6));
  GAUDI_EXPECT(std::abs(saddle_shift.k1 * saddle_shift.k2 + real(4)) <
               real(1e-6));

  const hepworth::block::shape_w identity =
      hepworth::block::map_ring_log_skew(pair, {}, real(0));
  GAUDI_ASSERT(identity.ok);
  GAUDI_EXPECT(std::abs(identity.k1 - pair.k1) < real(1e-10));
  GAUDI_EXPECT(std::abs(identity.k2 - pair.k2) < real(1e-10));
}

GAUDI_TEST(absolute_ratio_pair_sets_kappa_ratio_and_holds_product) {
  real o1 = 0;
  real o2 = 0;
  GAUDI_ASSERT(hepworth::block::absolute_ratio_pair(real(2), real(2), real(4),
                                                    o1, o2));
  GAUDI_EXPECT(std::abs(o1 - real(4)) < real(1e-6));
  GAUDI_EXPECT(std::abs(o2 - real(1)) < real(1e-6));
  GAUDI_EXPECT(std::abs(o1 * o2 - real(4)) < real(1e-6));
  GAUDI_ASSERT(hepworth::block::absolute_ratio_pair(real(4), real(1), real(4),
                                                    o1, o2));
  GAUDI_EXPECT(std::abs(o1 - real(4)) < real(1e-12));
  GAUDI_EXPECT(std::abs(o2 - real(1)) < real(1e-12));
  // `ratio` 1 equalizes |κ| and holds the product.
  GAUDI_ASSERT(hepworth::block::absolute_ratio_pair(real(4), real(1), real(1),
                                                    o1, o2));
  GAUDI_EXPECT(std::abs(std::abs(o1 / o2) - real(1)) < real(1e-6));
  GAUDI_EXPECT(std::abs(o1 * o2 - real(4)) < real(1e-6));
}

GAUDI_TEST(cross_ratio_spoke_map_keeps_scale_and_sets_gaps) {
  hepworth::block::linear_map2 M;
  M.m00 = real(2);
  M.m01 = real(0.3);
  M.m10 = real(-0.1);
  M.m11 = real(2.4);
  M.ok = true;
  const hepworth::block::linear_map_parts before =
      hepworth::block::split_linear_map(M);
  const hepworth::block::linear_map2 O =
      hepworth::block::cross_ratio_spoke_map(M, before.omega + real(0.2),
                                             before.h00 + real(0.1));
  GAUDI_ASSERT(O.ok);
  const hepworth::block::linear_map_parts after =
      hepworth::block::split_linear_map(O);
  GAUDI_EXPECT(std::abs(after.sigma - before.sigma) < real(1e-10));
  GAUDI_EXPECT(std::abs(after.omega - (before.omega + real(0.2))) < real(1e-10));
  GAUDI_EXPECT(std::abs(after.h00 - (before.h00 + real(0.1))) < real(1e-10));
}

namespace axis_ratio_detail {

/// Eight spokes on the unit circle in the xy plane, slightly off the plane,
/// with a chart whose k_hi tangent is x. `uneven` weighs one spoke more,
/// which makes the weighted ring anisotropic (the rows weigh spoke k by
/// `w_k^2`), so a free skew leans toward it.
struct ring {
  std::vector<vec3> rest;
  std::vector<real> weights;
  std::vector<index_t> ids;
  std::vector<hepworth::block::vertex_chart> charts{1};
  explicit ring(bool uneven = false) {
    rest.push_back(vec3::Zero());
    for (int k = 0; k < 8; ++k) {
      const real a = real(2) * real(M_PI) * k / 8 + real(0.1);
      rest.push_back(vec3(std::cos(a), std::sin(a), real(0.05) * std::sin(real(3) * a)));
    }
    weights.assign(8, real(1));
    if (uneven)
      weights[2] = real(1.5);
    for (size_t i = 0; i < rest.size(); ++i)
      ids.push_back(static_cast<index_t>(i));
    charts[0].ok = true;
    charts[0].T = vec3::UnitX();
    charts[0].B = vec3::UnitY();
    charts[0].N = vec3::UnitZ();
    charts[0].h = vec3(real(3), real(0), real(1));
  }
  /// The rest ring's tangent part mapped by `U = [[a, b], [b, c]]` in the
  /// frame of `u = (cos th, sin th)`, then turned rigidly.
  std::vector<vec3> mapped(real th, real a, real b, real c, const mat3 &R) const {
    const vec3 u(std::cos(th), std::sin(th), 0), v(-std::sin(th), std::cos(th), 0);
    std::vector<vec3> out;
    for (const vec3 &d : rest) {
      const real x = d.dot(u), y = d.dot(v);
      out.push_back(R * ((a * x + b * y) * u + (b * x + c * y) * v + d.z() * vec3::UnitZ()));
    }
    return out;
  }
  /// The row's projected spokes for live positions `live`.
  std::vector<vec3> project(hepworth::block::axis_ratio &row, std::vector<vec3> live) const {
    std::vector<vec3> mass(live.size(), vec3::Ones()), vel(live.size(), vec3::Zero()),
        force(live.size(), vec3::Zero());
    auto pos = hepworth::vec3_block::create(mass, live, vel, force);
    row._blocks = {pos};
    vecX q;
    pos->map_to_x(q);
    index_t id0 = 0;
    std::vector<trip> trips;
    row.fill_A(id0, trips);
    vecX p = vecX::Zero(id0);
    row.project(q, p);
    return row.rest();
  }
  /// Spoke rms ratio along `th` over across it.
  real ratio(const std::vector<vec3> &spokes, real th) const {
    const vec3 u(std::cos(th), std::sin(th), 0), v(-std::sin(th), std::cos(th), 0);
    real nu = 0, nv = 0;
    for (size_t k = 0; k < spokes.size(); ++k) {
      nu += weights[k] * std::pow(spokes[k].dot(u), 2);
      nv += weights[k] * std::pow(spokes[k].dot(v), 2);
    }
    return std::sqrt(nu / nv);
  }
  /// The stretch ratio an absolute row asks for along `th`: `ratio` times
  /// the rest ring's weighted rms across over along.
  real absolute_t(real ratio, real th) const {
    return ratio / this->ratio(spokes_of_rest(), th);
  }
  std::vector<vec3> spokes_of_rest() const {
    return std::vector<vec3>(rest.begin() + 1, rest.end());
  }
  /// A row with the area free unless a band is given.
  hepworth::block::axis_ratio::ptr row(
      real ratio, bool absolute, real angle,
      hepworth::block::area_band area = hepworth::block::area_band::free()) const {
    return hepworth::block::axis_ratio::create(ids, weights, rest, ratio, absolute, angle,
                                               real(1), std::vector<hepworth::sim_block::ptr>{},
                                               &charts, absolute ? "axis_ratio" : "shear",
                                               area);
  }
  /// Weighted spoke area of `spokes` in the xy plane over the rest ring's.
  real area_ratio(const std::vector<vec3> &spokes) const {
    mat2 M = mat2::Zero(), M0 = mat2::Zero();
    for (size_t k = 0; k < spokes.size(); ++k) {
      const vec2 d(spokes[k].x(), spokes[k].y()), d0(rest[k + 1].x(), rest[k + 1].y());
      M += weights[k] * d * d.transpose();
      M0 += weights[k] * d0 * d0.transpose();
    }
    return std::sqrt(M.determinant() / M0.determinant());
  }
};

inline std::vector<vec3> spokes_of(const std::vector<vec3> &x) {
  std::vector<vec3> out;
  for (size_t i = 1; i < x.size(); ++i)
    out.push_back(x[i] - x[0]);
  return out;
}

} // namespace axis_ratio_detail

GAUDI_TEST(axis_ratio_locks_the_ratio_and_leaves_skew_area_and_rotation) {
  axis_ratio_detail::ring r;
  auto row = r.row(real(2), true, real(0));
  // From the round ring it reaches the locked ratio along k_hi (x).
  const std::vector<vec3> got = r.project(*row, r.rest);
  std::cout << "  axis_ratio from round: " << r.ratio(got, 0) << " (target 2)" << std::endl;
  GAUDI_EXPECT(std::abs(r.ratio(got, 0) - real(2)) < real(0.02));
  // Any ring already at the ratio, whatever its skew, area and rotation, is
  // left where it is.
  const mat3 R = Eigen::AngleAxis<real>(real(0.4), vec3(1, 2, 2).normalized()).matrix();
  for (const real b : {real(0), real(0.3), real(-0.5)})
    for (const real c : {real(0.7), real(1.3)}) {
      const std::vector<vec3> live = r.mapped(0, r.absolute_t(real(2), 0) * c, b * c, c, R);
      const std::vector<vec3> held = r.project(*row, live);
      const std::vector<vec3> want = axis_ratio_detail::spokes_of(live);
      real err = 0;
      for (size_t k = 0; k < want.size(); ++k)
        err = std::max(err, (held[k] - want[k]).norm());
      GAUDI_EXPECT(err < real(1e-8));
    }
}

GAUDI_TEST(shear_drives_its_axis_and_leaves_the_principal_ratio) {
  axis_ratio_detail::ring r;
  const real g = real(0.25), th = real(0.25) * real(M_PI);
  auto row = r.row((1 + g) / (1 - g), false, th);
  const std::vector<vec3> got = r.project(*row, r.rest);
  std::cout << "  shear along 45: " << r.ratio(got, th) << " (target " << (1 + g) / (1 - g)
            << "), principal " << r.ratio(got, 0) << std::endl;
  GAUDI_EXPECT(std::abs(r.ratio(got, th) - (1 + g) / (1 - g)) < real(0.02));
  GAUDI_EXPECT(std::abs(r.ratio(got, 0) - real(1)) < real(0.02));
  // At 45 degrees the principal ratio is the shear axis's skew, which it
  // leaves: an anisotropic ring already sheared is held as it is.
  const mat3 R = Eigen::AngleAxis<real>(real(-0.7), vec3(0, 1, 3).normalized()).matrix();
  for (const real b : {real(0.2), real(-0.35)}) {
    const real c = real(1.1);
    const std::vector<vec3> live = r.mapped(th, (1 + g) / (1 - g) * c, b * c, c, R);
    const std::vector<vec3> held = r.project(*row, live);
    const std::vector<vec3> want = axis_ratio_detail::spokes_of(live);
    real err = 0;
    for (size_t k = 0; k < want.size(); ++k)
      err = std::max(err, (held[k] - want[k]).norm());
    GAUDI_EXPECT(err < real(1e-8));
  }
}

/// The default band keeps the ring's area: from round the ratio is reached
/// by stretching along and shrinking across, not by shrinking across alone;
/// a ring at the ratio and the rest area is held.
GAUDI_TEST(axis_ratio_keeps_area_by_default) {
  axis_ratio_detail::ring r;
  auto kept = r.row(real(1.5), true, real(0), {1, 1});
  auto free = r.row(real(1.5), true, real(0));
  const std::vector<vec3> a = r.project(*kept, r.rest), b = r.project(*free, r.rest);
  std::cout << "  ratio 1.5 from round: kept area " << r.area_ratio(a) << " ratio " << r.ratio(a, 0)
            << ", free area " << r.area_ratio(b) << " ratio " << r.ratio(b, 0) << std::endl;
  GAUDI_EXPECT(std::abs(r.area_ratio(a) - real(1)) < real(1e-3));
  GAUDI_EXPECT(std::abs(r.ratio(a, 0) - real(1.5)) < real(0.02));
  GAUDI_EXPECT(r.area_ratio(b) < real(0.95));
  // At the ratio and the rest area (any skew, any rotation): held.
  const mat3 R = Eigen::AngleAxis<real>(real(0.5), vec3(1, 1, 0).normalized()).matrix();
  const real t = r.absolute_t(real(1.5), 0), bb = real(0.2);
  const real c = std::sqrt((real(1) + bb * bb) / t);
  const std::vector<vec3> live = r.mapped(0, t * c, bb, c, R);
  const std::vector<vec3> held = r.project(*kept, live);
  const std::vector<vec3> want = axis_ratio_detail::spokes_of(live);
  real err = 0;
  for (size_t k = 0; k < want.size(); ++k)
    err = std::max(err, (held[k] - want[k]).norm());
  GAUDI_EXPECT(err < real(1e-6));
}

/// The face shear's premise and its row. With a triangle's edge normals
/// held, an in-plane map `F` of its vertices takes its CR normal derivative
/// `sum n_i grad(phi_i)^T` to that times `F^-1` (in the face basis); the
/// stretch from `face_shear_stretch` carries `S` to the sheared `S*`
/// exactly, keeps area, and is the identity at zero shear; the row holds a
/// triangle already at `R U` of its rest and passes the variational audit.
GAUDI_TEST(face_shear_moves_the_cr_tensor_and_holds_its_triangle) {
  using m2 = Eigen::Matrix<real, 2, 2>;
  const std::array<vec3, 3> X{vec3(0, 0, 0), vec3(1.2, 0.1, 0), vec3(0.3, 0.9, 0)};
  const std::array<vec3, 3> nE{vec3(0.2, -0.1, 1).normalized(), vec3(-0.3, 0.25, 1).normalized(),
                               vec3(0.05, 0.4, 1).normalized()};
  // dN in the face basis (u, v) = (x, y): CR basis of local edge i (vertex i
  // to i + 1), as the Dupin rows build it.
  auto dN = [&](const std::array<vec3, 3> &P) {
    const vec3 nf = (P[1] - P[0]).cross(P[2] - P[0]);
    const real A = real(0.5) * nf.norm();
    const vec3 n = nf.normalized();
    mat3 D = mat3::Zero();
    for (int i = 0; i < 3; ++i)
      D += nE[static_cast<size_t>(i)] *
           (-n.cross(P[static_cast<size_t>((i + 1) % 3)] - P[static_cast<size_t>(i)]) / A).transpose();
    m2 M;
    M << D(0, 0), D(0, 1), D(1, 0), D(1, 1);
    return M;
  };
  m2 Fm;
  Fm << 1.15, 0.2, -0.05, 0.9;
  Fm /= std::sqrt(Fm.determinant());
  std::array<vec3, 3> Y;
  for (int i = 0; i < 3; ++i) {
    const Eigen::Matrix<real, 2, 1> p = Fm * Eigen::Matrix<real, 2, 1>(X[i].x(), X[i].y());
    Y[static_cast<size_t>(i)] = vec3(p[0], p[1], 0);
  }
  const real cr_gap = (dN(Y) - dN(X) * Fm.inverse()).norm() / dN(X).norm();
  // The stretch algebra on a symmetric tensor.
  m2 S;
  S << 2.5, 0.4, 0.4, 0.8;
  m2 U, F, Sstar;
  GAUDI_ASSERT(hepworth::block::face_shear_stretch(S, real(0.2), real(0.4), U, real(1e-2), &F, &Sstar));
  const real alg_gap = (S * F.inverse() - Sstar).norm() / Sstar.norm();
  m2 U0;
  GAUDI_ASSERT(hepworth::block::face_shear_stretch(S, real(0), real(0.4), U0));
  std::cout << "  face shear: CR relation gap " << cr_gap << ", S F^-1 vs S* " << alg_gap
            << ", det U " << U.determinant() << ", |U0 - I| " << (U0 - m2::Identity()).norm()
            << std::endl;
  GAUDI_EXPECT(cr_gap < real(1e-12));
  GAUDI_EXPECT(alg_gap < real(1e-12));
  GAUDI_EXPECT(std::abs(U.determinant() - real(1)) < real(1e-12));
  GAUDI_EXPECT((U - U.transpose()).norm() < real(1e-12));
  GAUDI_EXPECT((U0 - m2::Identity()).norm() < real(1e-12));
  // The row: a triangle at R U of its rest is held; a perturbed one passes the audit.
  std::vector<vec3> rest{X[0], X[1], X[2]};
  const mat3 R = Eigen::AngleAxis<real>(real(0.6), vec3(1, -1, 2).normalized()).matrix();
  // U acts in the row's face basis (u along the first edge).
  const vec3 fu = (X[1] - X[0]).normalized(), fv = vec3::UnitZ().cross(fu);
  std::vector<vec3> live(3);
  for (int i = 0; i < 3; ++i) {
    const vec3 e = X[i] - X[0];
    const Eigen::Matrix<real, 2, 1> p = U * Eigen::Matrix<real, 2, 1>(e.dot(fu), e.dot(fv));
    live[static_cast<size_t>(i)] = R * (p[0] * fu + p[1] * fv) + vec3(0.3, -0.2, 0.5);
  }
  std::vector<vec3> mass(3, vec3::Ones()), vel(3, vec3::Zero()), force(3, vec3::Zero());
  auto pos = hepworth::vec3_block::create(mass, live, vel, force);
  auto row = hepworth::block::face_stretch_target::create(
      std::array<index_t, 3>{0, 1, 2}, rest, U, real(1),
      std::vector<hepworth::sim_block::ptr>{pos}, "face_shear");
  vecX q;
  pos->map_to_x(q);
  index_t id0 = 0;
  std::vector<trip> trips;
  row->fill_A(id0, trips);
  vecX Aq = vecX::Zero(id0), p = vecX::Zero(id0);
  for (const trip &t : trips)
    Aq[t.row()] += t.value() * q[t.col()];
  row->project(q, p);
  GAUDI_EXPECT((Aq - p).norm() < real(1e-12));
  for (vec3 &x : live)
    x += vec3(0.03, -0.02, 0.04) * real(&x - &live[0] + 1);
  auto pos2 = hepworth::vec3_block::create(mass, live, vel, force);
  auto row2 = hepworth::block::face_stretch_target::create(
      std::array<index_t, 3>{0, 1, 2}, rest, U, real(1),
      std::vector<hepworth::sim_block::ptr>{pos2}, "face_shear");
  vecX q2;
  pos2->map_to_x(q2);
  for (const auto &a : hepworth::block::audit_variational({row2}, q2))
    GAUDI_EXPECT(a.max < real(1e-4));
}

GAUDI_TEST(axis_ratio_rows_are_variational) {
  axis_ratio_detail::ring r(true);
  std::vector<vec3> live = r.mapped(real(0.3), real(1.4), real(0.2), real(0.9),
                                    Eigen::AngleAxis<real>(real(0.2), vec3::UnitY()).matrix());
  std::mt19937 rng(4);
  std::uniform_real_distribution<real> pick(real(-0.05), real(0.05));
  for (vec3 &p : live)
    p += vec3(pick(rng), pick(rng), pick(rng));
  std::vector<vec3> mass(live.size(), vec3::Ones()), vel(live.size(), vec3::Zero()),
      force(live.size(), vec3::Zero());
  auto pos = hepworth::vec3_block::create(mass, live, vel, force);
  const std::vector<hepworth::sim_block::ptr> blocks{pos};
  std::vector<hepworth::projection_constraint::ptr> rows{
      hepworth::block::axis_ratio::create(r.ids, r.weights, r.rest, real(2), true, real(0),
                                          real(1), blocks, &r.charts, "axis_ratio"),
      hepworth::block::axis_ratio::create(r.ids, r.weights, r.rest, real(1.6), false,
                                          real(0.5), real(1), blocks, &r.charts, "shear"),
      hepworth::block::axis_ratio::create(r.ids, r.weights, r.rest, real(1.6), false,
                                          real(0.5), real(1), blocks, &r.charts, "shear_free",
                                          hepworth::block::area_band::free())};
  vecX q;
  pos->map_to_x(q);
  for (const auto &a : hepworth::block::audit_variational(rows, q)) {
    std::cout << "  " << a.name << " variational mismatch max " << a.max << std::endl;
    GAUDI_EXPECT(a.max < real(1e-4));
  }
}

GAUDI_TEST(face_operator_matches_the_analytic_torus) {
  const real R = 1.0;
  const real r = 0.33;
  asawa::shell::shell::ptr M = asawa::shell::load_torus(R, r, 48, 24);
  auto cr = std::make_shared<hepworth::block::shell_cr_block>(M);
  hepworth::block::solver_context ctx;
  cr->prepare(ctx);
  std::vector<hepworth::sim_block::ptr> emitted;
  cr->emit_dof_blocks(emitted);
  GAUDI_ASSERT(cr->face_dof());
  std::vector<vec3> &x = asawa::get_vec_data(*M, 0);
  const std::vector<vec3> x0 = x;
  std::vector<vec3> vel(x.size(), vec3::Zero());
  std::vector<vec3> force(x.size(), vec3::Zero());
  std::vector<vec3> mass(x.size(), vec3::Ones());
  auto pos = hepworth::vec3_block::create(mass, x, vel, force);
  std::vector<hepworth::cache_compute::ptr> computes;
  hepworth::block::init_shell_caches(*M, *cr, pos, computes);
  std::vector<hepworth::projection_constraint::ptr> constraints;
  std::vector<hepworth::sim_block::ptr> blocks{pos, cr->face_dof()};
  hepworth::block::init_area(*M, constraints, x, real(1), blocks);
  hepworth::block::init_face_shape_pin(*M, *cr, pos, constraints, real(1));
  hepworth::block::init_bending(*M, constraints, x, real(1), blocks);
  std::vector<hepworth::block::objective> spines{
      hepworth::block::spine_fairing{hepworth::block::spine_kind::cylinder}};
  hepworth::block::init_face_fair(*M, *cr, pos, constraints, real(0.4), spines);
  hepworth::block::face_shape_pin::ptr pin;
  for (const hepworth::projection_constraint::ptr &constraint : constraints) {
    pin = std::dynamic_pointer_cast<hepworth::block::face_shape_pin>(constraint);
    if (pin)
      break;
  }
  GAUDI_ASSERT(pin);
  hepworth::block::projection_solver solver;
  solver.set_constraints(constraints);
  solver.set_computes(computes);
  solver.step(blocks, real(0.1), real(0.1), 1, &ctx);
  vecX q_view;
  pos->map_to_x(q_view);
  pin->record(q_view);

  const real sign = hepworth::block::k_cr_shape_sign;
  int n = 0;
  real max_tube_deg = 0;
  real max_major_deg = 0;
  real max_tight_deg = 0;
  real max_tube_k = 0;
  real max_major_k = 0;
  real max_move = 0;
  for (int i = 0; i < static_cast<int>(x.size()); ++i)
    max_move = std::max(
        max_move,
        (x[static_cast<size_t>(i)] - x0[static_cast<size_t>(i)]).norm());
  for (const hepworth::block::projected_face_frame &fr : pin->_projected) {
    if (!fr.ok)
      continue;
    const analytic_torus exact = torus_at(R, r, fr.origin);
    const vec3 tube = project_to_plane(exact.tube, fr.N);
    const vec3 major = project_to_plane(exact.major, fr.N);
    const bool tight_is_e1 = std::abs(fr.k1) >= std::abs(fr.k2);
    const vec3 tight = tight_is_e1 ? fr.e1 : fr.e2;
    const real tube_deg =
        std::min(axis_degrees(fr.e1, tube), axis_degrees(fr.e2, tube));
    const real major_deg =
        std::min(axis_degrees(fr.e1, major), axis_degrees(fr.e2, major));
    const real tube_k = std::abs(frame_curvature(fr, tube) - sign * exact.k_tube);
    const real major_k =
        std::abs(frame_curvature(fr, major) - sign * exact.k_major);
    max_tube_deg = std::max(max_tube_deg, tube_deg);
    max_major_deg = std::max(max_major_deg, major_deg);
    max_tight_deg = std::max(max_tight_deg, axis_degrees(tight, tube));
    max_tube_k = std::max(max_tube_k, tube_k);
    max_major_k = std::max(max_major_k, major_k);
    ++n;
  }
  GAUDI_EXPECT(max_move < real(1e-8));
  GAUDI_EXPECT(n == static_cast<int>(pin->_projected.size()));
  GAUDI_EXPECT(n == 48 * 24 * 2);
  GAUDI_EXPECT(max_tube_deg < real(3));
  GAUDI_EXPECT(max_major_deg < real(3));
  GAUDI_EXPECT(max_tight_deg < real(3));
  GAUDI_EXPECT(max_tube_k < real(0.05));
  GAUDI_EXPECT(max_major_k < real(0.10));
}

GAUDI_TEST(diamond_polar_chart_aligns_with_the_edge) {
  const vec3 p0(0, 0, 0);
  const vec3 p1(1, 0, 0);
  const vec3 a(real(0.5), real(0.5), 0);
  const vec3 b(real(0.5), real(-0.5), 0);
  const std::array<int, 3> t0{{0, 1, 2}};
  const std::array<int, 3> t1{{0, 3, 1}};
  vec3 N, T, B;
  GAUDI_ASSERT(hepworth::block::diamond_polar_chart(p0, p1, a, b, t0, t1, 0, 1, 2,
                                                    3, N, T, B));
  GAUDI_EXPECT(std::abs(N.z()) > real(0.9));
  GAUDI_EXPECT(std::abs(std::abs(T.dot(vec3::UnitX())) - real(1)) < real(1e-6));
  GAUDI_EXPECT(std::abs(T.dot(N)) < real(1e-9));
  GAUDI_EXPECT(std::abs(B.dot(N)) < real(1e-9));
}

GAUDI_TEST(blend_shape_on_diamond_keeps_shared_cylinder) {
  const vec3 p0(0, 0, 0);
  const vec3 p1(1, 0, 0);
  const vec3 a(real(0.5), real(0.5), 0);
  const vec3 b(real(0.5), real(-0.5), 0);
  const std::array<int, 3> t0{{0, 1, 2}};
  const std::array<int, 3> t1{{0, 3, 1}};
  const vec3 N(0, 0, 1);
  const vec3 Th(1, 0, 0);
  const vec3 h(real(0), real(0), real(-2));
  vec3 hout, Tout;
  GAUDI_ASSERT(hepworth::block::blend_shape_on_diamond(
      p0, p1, a, b, t0, t1, 0, 1, 2, 3, real(0.5), N, N, h, h, &Th, &Th, hout,
      Tout));
  GAUDI_EXPECT((hout - h).norm() < real(1e-8));
  GAUDI_EXPECT(std::abs(std::abs(Tout.dot(Th)) - real(1)) < real(1e-9));
}

GAUDI_TEST(blend_shape_on_diamond_keeps_aligned_off_edge_frames) {
  // Neighbors share a chart 45° off the edge. Destination T must keep that
  // angle — not snap to the edge from a bogus covariance polar.
  const vec3 p0(0, 0, 0);
  const vec3 p1(1, 0, 0);
  const vec3 a(real(0.5), real(0.5), 0);
  const vec3 b(real(0.5), real(-0.5), 0);
  const std::array<int, 3> t0{{0, 1, 2}};
  const std::array<int, 3> t1{{0, 3, 1}};
  const vec3 N(0, 0, 1);
  const vec3 Th = vec3(1, 1, 0).normalized();
  const vec3 h(real(-2), real(0), real(0));
  vec3 hout, Tout;
  GAUDI_ASSERT(hepworth::block::blend_shape_on_diamond(
      p0, p1, a, b, t0, t1, 0, 1, 2, 3, real(0.5), N, N, h, h, &Th, &Th, hout,
      Tout));
  GAUDI_EXPECT((hout - h).norm() < real(1e-8));
  GAUDI_EXPECT(std::abs(std::abs(Tout.dot(Th)) - real(1)) < real(1e-6));
  GAUDI_EXPECT(std::abs(Tout.dot(vec3::UnitX())) < real(0.9));
}

GAUDI_TEST(blend_shape_on_edge_keeps_a_shared_cylinder) {
  // Same chart H at both ends of an edge: world blend must return that H,
  // not a garbage lerp of mismatched bases.
  const vec3 x0(0, 0, 0);
  const vec3 x1(1, 0, 0);
  const vec3 N(0, 0, 1);
  const vec3 T(1, 0, 0);
  const vec3 h(real(0), real(0), real(-2));
  vec3 hout, Tout;
  GAUDI_ASSERT(hepworth::block::blend_shape_on_edge(x0, x1, real(0.5), N, N, h, h,
                                                   &T, &T, hout, Tout));
  GAUDI_EXPECT((hout - h).norm() < real(1e-9));
  GAUDI_EXPECT(std::abs(std::abs(Tout.dot(T)) - real(1)) < real(1e-9));
}

GAUDI_TEST(blend_shape_on_edge_survives_a_chart_twist) {
  // Endpoint charts differ by 90°; lerping coeffs would mix k11/k22. World
  // blend keeps the principal values.
  const vec3 x0(0, 0, 0);
  const vec3 x1(1, 0, 0);
  const vec3 N(0, 0, 1);
  const vec3 T0(1, 0, 0);
  const vec3 T1(0, 1, 0);
  const vec3 h0(real(-2), real(0), real(0)); // κ along T0
  const vec3 h1(real(0), real(0), real(-2)); // same κ along T1 = former B0
  vec3 hout, Tout;
  GAUDI_ASSERT(hepworth::block::blend_shape_on_edge(
      x0, x1, real(0.5), N, N, h0, h1, &T0, &T1, hout, Tout));
  real k1 = 0, k2 = 0, ex = 0, ey = 0;
  hepworth::block::shape_sym2(hout[0], hout[1], hout[2], k1, k2, ex, ey);
  GAUDI_EXPECT(std::abs(k1 + real(2)) < real(1e-6) ||
               std::abs(k2 + real(2)) < real(1e-6));
  GAUDI_EXPECT(std::abs(k1) < real(1e-6) || std::abs(k2) < real(1e-6));
}

GAUDI_TEST(store_solved_tangents_reseed_replaces_carried_u) {
  // After remesh, prepare loads carried u. Reseed must overwrite those slots
  // so H/chart_T and the edge field stay one basis.
  auto M = asawa::shell::load_cube();
  GAUDI_ASSERT(M);
  asawa::shell::triangulate(*M);
  hepworth::block::shell_cr_block cr(M);
  hepworth::block::solver_context ctx;
  cr.prepare(ctx);
  GAUDI_ASSERT(cr.topology.nE > 0);
  const int nE = cr.topology.nE;
  cr.u.assign(static_cast<size_t>(2 * nE), real(0));
  for (int e = 0; e < nE; ++e) {
    cr.u[static_cast<size_t>(e)] = real(1);
    cr.u[static_cast<size_t>(e + nE)] = real(0);
  }
  cr.guidance_tangents = false;
  std::vector<real> coeff(static_cast<size_t>(2 * nE), real(0));
  for (int e = 0; e < nE; ++e) {
    coeff[static_cast<size_t>(e)] = real(0);
    coeff[static_cast<size_t>(e + nE)] = real(1);
  }
  cr.store_solved_tangents(coeff, real(0), /*replace_all=*/false);
  GAUDI_EXPECT(std::abs(cr.u[0] - real(1)) < real(1e-12));
  cr.store_solved_tangents(coeff, real(0), /*replace_all=*/true);
  GAUDI_EXPECT(std::abs(cr.u[0]) < real(1e-12));
  GAUDI_EXPECT(std::abs(cr.u[static_cast<size_t>(nE)] - real(1)) < real(1e-12));
}

GAUDI_TEST(sync_charts_from_transports_H_on_a_chart_turn) {
  auto M = asawa::shell::load_cube();
  GAUDI_ASSERT(M);
  asawa::shell::triangulate(*M);
  hepworth::block::shell_shape_block shape(M);
  const int n = M->vert_count();
  GAUDI_ASSERT(n > 0);
  shape.descriptor()->fit_sides(n);
  shape.chart_T()[0] = vec3::UnitX();
  shape.shape()[0] = vec3(real(-2), real(0), real(0));
  std::vector<hepworth::block::vertex_chart> charts(static_cast<size_t>(n));
  charts[0].ok = true;
  charts[0].T = vec3::UnitY();
  charts[0].B = -vec3::UnitX();
  charts[0].N = vec3::UnitZ();
  charts[0].h = vec3(real(0), real(0), real(-2));
  shape.sync_charts_from(charts);
  GAUDI_EXPECT(std::abs(std::abs(shape.chart_T()[0].dot(vec3::UnitY())) - real(1)) <
               real(1e-9));
  // κ stayed on the world hoop: in the new chart that is k22.
  GAUDI_EXPECT(std::abs(shape.shape()[0][2] + real(2)) < real(1e-6));
  GAUDI_EXPECT(std::abs(shape.shape()[0][0]) < real(1e-6));
}

GAUDI_TEST(face_fair_spoke_lands_on_the_projection) {
  std::vector<vec3> x(4), vel(4, vec3::Zero()), force(4, vec3::Zero());
  std::vector<vec3> mass(4, vec3::Ones());
  x[0] = vec3::Zero();
  x[1] = vec3(2, 0, 0);
  x[2] = vec3(real(1.3), real(0.5), 0);
  x[3] = vec3(1, real(-0.4), 0);
  std::vector<vec3> sh{vec3(real(-2), real(0), real(-0.5))};
  std::vector<vec3> sv(1, vec3::Zero()), sf(1, vec3::Zero()), sm(1, vec3::Ones());
  auto pos = hepworth::vec3_block::create(mass, x, vel, force);
  auto shape = hepworth::vec3_block::create(sm, sh, sv, sf);
  hepworth::block::face_patch patch;
  patch.face = 0;
  patch.vid = {0, 1, 2, 3, -1, -1};
  patch.nvert = 4;
  patch.loc[0] = {0, 1};
  patch.loc[1] = {1, 2};
  patch.loc[2] = {2, 0};
  patch.loc[3] = {0, 3};
  patch.loc[4] = {3, 1};
  patch.nedge = 5;
  std::vector<hepworth::block::objective> objectives{
      hepworth::block::spine_fairing{hepworth::block::spine_kind::cylinder}};
  auto fair = hepworth::block::face_fair::create(
      patch, real(1), objectives,
      std::vector<hepworth::sim_block::ptr>{pos, shape});
  vecX q;
  pos->map_to_x(q);
  shape->map_to_x(q);
  index_t id0 = 0;
  std::vector<trip> trips;
  fair->fill_A(id0, trips);
  GAUDI_EXPECT(id0 == 15);
  vecX p = vecX::Zero(id0);
  fair->project(q, p);
  vec3 N, T, B;
  GAUDI_ASSERT(hepworth::block::face_tri_frame(x[0], x[1], x[2], N, T, B));
  const vec3 center = (x[0] + x[1] + x[2]) / real(3);
  const vec3 face_edge(p[0], p[1], p[2]);
  GAUDI_EXPECT((face_edge - (x[1] - x[0])).norm() < real(1e-6));
  const auto axes = hepworth::block::monge_axes_of(hepworth::block::shape_with_axes(
      hepworth::block::shape_w_from_coeffs(sh[0])));
  hepworth::block::spine_curvature sc;
  sc.tube = axes.hx * T + axes.hy * B;
  sc.ruling = axes.lx * T + axes.ly * B;
  sc.normal = N;
  sc.k_hi = axes.k_hi;
  sc.k_lo = axes.k_lo;
  const auto land = hepworth::block::choose_spine_shapes(
      sc, std::vector<vec3>{x[3] - center},
      hepworth::block::spine_set_of(objectives));
  GAUDI_ASSERT(land.fit == hepworth::block::spine_kind::cylinder);
  GAUDI_ASSERT(!land.ok.empty() && land.ok[0]);
  const vec3 target = (center + land.hit[0]) - x[0];
  const vec3 spoke(p[9], p[10], p[11]);
  GAUDI_EXPECT((spoke - target).norm() < real(1e-6));
}

GAUDI_TEST(face_fair_holds_a_point_on_the_graph) {
  const vec3 x0 = vec3::Zero();
  const vec3 x1(1, 0, 0);
  const vec3 x2(real(0.2), real(0.8), 0);
  const vec3 center = (x0 + x1 + x2) / real(3);
  const real tx = real(0.35);
  const vec3 x3 = center + vec3(tx, 0, real(0.5) * tx * tx);
  std::vector<vec3> x{x0, x1, x2, x3};
  std::vector<vec3> vel(4, vec3::Zero()), force(4, vec3::Zero());
  std::vector<vec3> mass(4, vec3::Ones());
  std::vector<vec3> sh{vec3(real(1), real(0), real(0))};
  std::vector<vec3> sv(1, vec3::Zero()), sf(1, vec3::Zero()), sm(1, vec3::Ones());
  auto pos = hepworth::vec3_block::create(mass, x, vel, force);
  auto shape = hepworth::vec3_block::create(sm, sh, sv, sf);
  hepworth::block::face_patch patch;
  patch.face = 0;
  patch.vid = {0, 1, 2, 3, -1, -1};
  patch.nvert = 4;
  patch.loc[0] = {0, 1};
  patch.loc[1] = {1, 2};
  patch.loc[2] = {2, 0};
  patch.loc[3] = {0, 3};
  patch.loc[4] = {3, 1};
  patch.nedge = 5;
  std::vector<hepworth::block::objective> objectives{
      hepworth::block::monge_fairing{hepworth::block::monge_shape::quadric}};
  auto fair = hepworth::block::face_fair::create(
      patch, real(1), objectives,
      std::vector<hepworth::sim_block::ptr>{pos, shape});
  vecX q;
  pos->map_to_x(q);
  shape->map_to_x(q);
  index_t id0 = 0;
  std::vector<trip> trips;
  fair->fill_A(id0, trips);
  vecX p = vecX::Zero(id0);
  fair->project(q, p);
  const vec3 face_edge(p[0], p[1], p[2]);
  const vec3 outer(p[9], p[10], p[11]);
  GAUDI_EXPECT((face_edge - (x1 - x0)).norm() < real(1e-8));
  GAUDI_EXPECT((outer - (x3 - x0)).norm() < real(1e-6));
}

GAUDI_TEST(svd_smear_keeps_a_flipped_tangent) {
  const vec3 N = vec3::UnitZ();
  const std::vector<vec3> samples{vec3(1, 0, 0), vec3(-1, 0, 0), vec3(0.2, 0, 0)};
  vec3 T;
  GAUDI_ASSERT(hepworth::block::svd_tangent_fit(samples, N, T));
  GAUDI_EXPECT(std::abs(std::abs(T.dot(vec3::UnitX())) - real(1)) < real(1e-6));
  GAUDI_EXPECT(std::abs(T.dot(N)) < real(1e-6));
}

GAUDI_TEST(albers_ring_fair_reads_only_positions) {
  asawa::shell::shell::ptr M = asawa::shell::load_torus(1.0, 0.33, 8, 6);
  std::vector<vec3> &x = asawa::get_vec_data(*M, 0);
  const int nv = static_cast<int>(x.size());
  std::vector<vec3> mass(static_cast<size_t>(nv), vec3::Ones());
  std::vector<vec3> vel(static_cast<size_t>(nv), vec3::Zero());
  std::vector<vec3> force(static_cast<size_t>(nv), vec3::Zero());
  auto pos = hepworth::vec3_block::create(mass, x, vel, force);
  std::vector<char> live(static_cast<size_t>(nv), 1);
  std::vector<hepworth::sim_block::ptr> blocks{pos};
  std::vector<hepworth::projection_constraint::ptr> constraints;
  const std::vector<hepworth::block::objective> shapes{
      hepworth::block::albers_fairing{hepworth::block::albers_shape::sphere}};
  hepworth::block::init_albers_fair(*M, constraints, x, live, real(1), blocks,
                                    shapes);
  GAUDI_ASSERT(!constraints.empty());
  vecX q;
  pos->map_to_x(q);
  index_t id0 = 0;
  std::vector<trip> trips;
  bool named = false;
  for (const auto &c : constraints) {
    GAUDI_EXPECT(c->name() == "albers_ring_fair");
    named = true;
    c->fill_A(id0, trips);
  }
  GAUDI_EXPECT(named);
  const index_t begin = pos->get_offset_idx(0);
  const index_t end = begin + 3 * nv;
  for (const trip &t : trips) {
    GAUDI_EXPECT(t.col() >= begin);
    GAUDI_EXPECT(t.col() < end);
  }
  vecX p = vecX::Zero(id0);
  for (const auto &c : constraints)
    c->project(q, p);
  GAUDI_EXPECT(p.allFinite());
  GAUDI_EXPECT(p.norm() > real(0));
}

GAUDI_TEST(ring_scalar_gradient_recovers_a_linear_field) {
  const real s = real(1);
  const std::vector<hepworth::block::ring_scalar_sample> samples{
      {real(7), real(2), real(0), real(1)},
      {real(21), real(0), real(4), real(1)},
      {real(-5), real(-2), real(0), real(2)},
      {real(100), real(1), real(0), real(0)},
  };
  const hepworth::block::ring_scalar_gradient g =
      hepworth::block::fit_ring_scalar_gradient(s, samples);
  GAUDI_EXPECT(g.ok);
  GAUDI_EXPECT(std::abs(g.gx - real(3)) < real(1e-10));
  GAUDI_EXPECT(std::abs(g.gy - real(5)) < real(1e-10));
  GAUDI_EXPECT(std::abs(g.div - real(14)) < real(1e-10));

  const std::vector<hepworth::block::ring_scalar_sample> flat{
      {s, real(1), real(0), real(1)},
      {s, real(0), real(1), real(2)},
  };
  const hepworth::block::ring_scalar_gradient z =
      hepworth::block::fit_ring_scalar_gradient(s, flat);
  GAUDI_EXPECT(z.ok);
  GAUDI_EXPECT(std::abs(z.gx) < real(1e-12));
  GAUDI_EXPECT(std::abs(z.gy) < real(1e-12));
  GAUDI_EXPECT(std::abs(z.div) < real(1e-12));

  const vec3 diag =
      hepworth::block::principal_coeffs(real(3), real(1), real(1), real(0),
                                        real(0), real(1));
  GAUDI_EXPECT((diag - vec3(real(3), real(0), real(1))).norm() < real(1e-12));
  const hepworth::block::shape_w back =
      hepworth::block::shape_w_from_coeffs(diag);
  GAUDI_EXPECT(back.ok);
  GAUDI_EXPECT(std::abs(back.k1 - real(3)) < real(1e-10));
  GAUDI_EXPECT(std::abs(back.k2 - real(1)) < real(1e-10));
}

GAUDI_TEST(invert_principal_frame_keeps_the_polar_tangent) {
  const vec3 T0 = vec3::UnitX();
  const vec3 B0 = vec3::UnitY();
  const vec3 N = vec3::UnitZ();
  const vec3 k0_tan = T0;
  const vec3 k1_tan = B0;
  vec3 T, B, h;
  GAUDI_EXPECT(hepworth::block::invert_principal_frame(
      real(4), real(1), k0_tan, k1_tan, real(1), real(0), real(0), real(1), N,
      T, B, h));
  GAUDI_EXPECT((T - T0).norm() < real(1e-10));
  GAUDI_EXPECT((B - B0).norm() < real(1e-10));
  GAUDI_EXPECT((h - vec3(real(4), real(0), real(1))).norm() < real(1e-10));
}

GAUDI_TEST(polar_smear_keeps_a_shared_edge_axis) {
  const vec3 N = vec3::UnitZ();
  const std::vector<vec3> samples{vec3(1, 0, 0), vec3(-1, 0, 0), vec3(0.2, 0, 0)};
  vec3 T;
  GAUDI_ASSERT(hepworth::block::polar_tangent(samples, N, T));
  GAUDI_EXPECT(std::abs(std::abs(T.dot(vec3::UnitX())) - real(1)) < real(1e-6));
  GAUDI_EXPECT(std::abs(T.dot(N)) < real(1e-6));
}

/// Edge/vertex frame gates against the analytic torus (land + remesh).
struct torus_tangent_report {
  int n_ok = 0;
  int n_vert = 0;
  int n_edge = 0;
  int n_edge_empty = 0;
  int n_ring_pairs = 0;
  int n_e1_swap = 0; ///< neighbor |e1·e1| near 0 (≈90° label flip)
  int n_e2_swap = 0;
  real min_e1_abs_dot = 1;
  real min_e2_abs_dot = 1;
  real max_e1_neighbor_deg = 0;
  real max_e2_neighbor_deg = 0;
  real max_edge_ring_deg = 0;
  real max_tube_deg = 0;
  real max_major_deg = 0;
  real max_tight_tube_deg = 0; ///< larger-|κ| axis vs analytic tube
  real max_tube_k = 0;
  real max_major_k = 0;
  real max_tube_radius_err = 0; ///< |1/|κ_tube| − r|
  real max_major_radius_err = 0;
};

/// Same axes the draw paints: `e1` is red (larger eig), `e2` green.
struct chart_principals {
  bool ok = false;
  vec3 e1 = vec3::Zero();
  vec3 e2 = vec3::Zero();
  real k1 = 0;
  real k2 = 0;
};

inline chart_principals extract_chart_principals(
    const hepworth::block::vertex_chart &c) {
  chart_principals out;
  const hepworth::block::shape_w W =
      hepworth::block::shape_with_axes(hepworth::block::shape_w_from_coeffs(c.h));
  if (!W.ok || !c.T.allFinite() || !c.B.allFinite() || !c.N.allFinite())
    return out;
  vec3 N = c.N.normalized();
  vec3 e1 = W.e1x * c.T + W.e1y * c.B;
  vec3 e2 = W.e2x * c.T + W.e2y * c.B;
  e1 -= N * e1.dot(N);
  e2 -= N * e2.dot(N);
  if (!(e1.squaredNorm() > real(1e-16)) || !(e2.squaredNorm() > real(1e-16)))
    return out;
  out.e1 = e1.normalized();
  out.e2 = e2.normalized();
  out.k1 = W.k1;
  out.k2 = W.k2;
  out.ok = true;
  return out;
}

/// Neighbor same-slot axes. A 90° label swap makes |e1_i·e1_j| ≈ 0 — that
/// must fail. Do not re-pick by |κ|; that hides the draw failure.
inline void measure_one_ring_labeled_principal_coherence(
    const asawa::shell::shell &M,
    const std::vector<hepworth::block::vertex_chart> &charts,
    torus_tangent_report &rep, real swap_abs_dot = real(0.35)) {
  const int nV = M.vert_count();
  for (int i = 0; i < nV; ++i) {
    if (static_cast<size_t>(i) >= charts.size() || !charts[static_cast<size_t>(i)].ok)
      continue;
    const chart_principals ai =
        extract_chart_principals(charts[static_cast<size_t>(i)]);
    if (!ai.ok)
      continue;
    const auto ring = M.get_one_ring(asawa::shell::vert_id(i));
    for (const asawa::shell::VertId nj : ring) {
      const int j = static_cast<int>(nj);
      if (j <= i || j >= nV)
        continue;
      if (static_cast<size_t>(j) >= charts.size() || !charts[static_cast<size_t>(j)].ok)
        continue;
      const chart_principals aj =
          extract_chart_principals(charts[static_cast<size_t>(j)]);
      if (!aj.ok)
        continue;
      const real d1 = std::abs(ai.e1.dot(aj.e1));
      const real d2 = std::abs(ai.e2.dot(aj.e2));
      ++rep.n_ring_pairs;
      rep.min_e1_abs_dot = std::min(rep.min_e1_abs_dot, d1);
      rep.min_e2_abs_dot = std::min(rep.min_e2_abs_dot, d2);
      rep.max_e1_neighbor_deg =
          std::max(rep.max_e1_neighbor_deg, axis_degrees(ai.e1, aj.e1));
      rep.max_e2_neighbor_deg =
          std::max(rep.max_e2_neighbor_deg, axis_degrees(ai.e2, aj.e2));
      if (d1 < swap_abs_dot)
        ++rep.n_e1_swap;
      if (d2 < swap_abs_dot)
        ++rep.n_e2_swap;
    }
  }
}

inline void measure_one_ring_edge_tangent_coherence(
    asawa::shell::shell &M, const std::vector<vec3> &x,
    const std::vector<real> &coeff,
    const std::vector<hepworth::block::vertex_chart> &charts,
    torus_tangent_report &rep) {
  const std::vector<vec3> edge_tan = hepworth::block::extract_tan(M, x, coeff);
  std::vector<int> slot_map;
  const int nE = kusama::build_compact_edge_dof_map(M, slot_map);
  rep.n_edge = nE;
  for (int e = 0; e < nE; ++e) {
    if (e >= static_cast<int>(edge_tan.size()) ||
        !(edge_tan[static_cast<size_t>(e)].squaredNorm() > real(1e-16)))
      ++rep.n_edge_empty;
  }
  const int nV = M.vert_count();
  for (int i = 0; i < nV; ++i) {
    if (static_cast<size_t>(i) >= charts.size() || !charts[static_cast<size_t>(i)].ok)
      continue;
    const asawa::shell::VertId vi = asawa::shell::vert_id(i);
    if (static_cast<int>(M.vbegin(vi)) < 0)
      continue;
    vec3 N = charts[static_cast<size_t>(i)].N;
    if (!N.allFinite() || !(N.squaredNorm() > real(1e-16)))
      continue;
    N.normalize();
    vec3 polar = charts[static_cast<size_t>(i)].T;
    polar -= N * polar.dot(N);
    if (!(polar.squaredNorm() > real(1e-16)))
      continue;
    polar.normalize();
    M.const_for_each_vertex(
        vi, [&](asawa::shell::CornerId ci, const asawa::shell::shell &) {
          const int slot = static_cast<int>(ci) / 2;
          const int e =
              slot >= 0 && slot < static_cast<int>(slot_map.size())
                  ? slot_map[static_cast<size_t>(slot)]
                  : -1;
          if (e < 0 || e >= static_cast<int>(edge_tan.size()))
            return;
          vec3 t = edge_tan[static_cast<size_t>(e)];
          t -= N * t.dot(N);
          if (!(t.squaredNorm() > real(1e-16)))
            return;
          t.normalize();
          rep.max_edge_ring_deg =
              std::max(rep.max_edge_ring_deg, axis_degrees(polar, t));
        });
  }
}

/// Chord midpoints from remesh leave the surface; snap back for analytic checks.
inline void snap_vertices_to_torus(std::vector<vec3> &x, real R, real r) {
  for (vec3 &p : x) {
    const real rho = std::hypot(p.x(), p.y());
    if (!(rho > real(1e-12)))
      continue;
    const real u = std::atan2(p.y(), p.x());
    const vec3 center(R * std::cos(u), R * std::sin(u), real(0));
    vec3 d = p - center;
    if (!(d.squaredNorm() > real(1e-16)))
      continue;
    p = center + r * d.normalized();
  }
}

inline void measure_charts_vs_analytic_torus(
    const asawa::shell::shell &M, const std::vector<vec3> &x,
    const std::vector<hepworth::block::vertex_chart> &charts, real R, real r,
    torus_tangent_report &rep) {
  const real sign = hepworth::block::k_cr_shape_sign;
  const int nV = M.vert_count();
  rep.n_vert = nV;
  for (int i = 0; i < nV; ++i) {
    if (i >= static_cast<int>(x.size()) || static_cast<size_t>(i) >= charts.size())
      continue;
    const hepworth::block::vertex_chart &c = charts[static_cast<size_t>(i)];
    if (!c.ok || !c.h.allFinite())
      continue;
    const chart_principals ax = extract_chart_principals(c);
    if (!ax.ok)
      continue;
    ++rep.n_ok;
    const analytic_torus exact = torus_at(R, r, x[static_cast<size_t>(i)]);
    const vec3 tube = project_to_plane(exact.tube, c.N);
    const vec3 major = project_to_plane(exact.major, c.N);
    const real tube_deg =
        std::min(axis_degrees(ax.e1, tube), axis_degrees(ax.e2, tube));
    const real major_deg =
        std::min(axis_degrees(ax.e1, major), axis_degrees(ax.e2, major));
    rep.max_tube_deg = std::max(rep.max_tube_deg, tube_deg);
    rep.max_major_deg = std::max(rep.max_major_deg, major_deg);
    const vec3 tight =
        std::abs(ax.k1) >= std::abs(ax.k2) ? ax.e1 : ax.e2;
    rep.max_tight_tube_deg =
        std::max(rep.max_tight_tube_deg, axis_degrees(tight, tube));
    const real k_tube = ax.k1 * std::pow(ax.e1.dot(tube), 2) +
                        ax.k2 * std::pow(ax.e2.dot(tube), 2);
    const real k_major = ax.k1 * std::pow(ax.e1.dot(major), 2) +
                         ax.k2 * std::pow(ax.e2.dot(major), 2);
    rep.max_tube_k =
        std::max(rep.max_tube_k, std::abs(k_tube - sign * exact.k_tube));
    rep.max_major_k =
        std::max(rep.max_major_k, std::abs(k_major - sign * exact.k_major));
    // Absolute 1/κ blows up when |κ| is small (outer equator). Gate relative
    // radius error only when the analytic |κ| is well away from zero.
    if (std::abs(k_tube) > real(1e-8) && std::abs(exact.k_tube) > real(0.5)) {
      const real r_tube = real(1) / std::abs(k_tube);
      rep.max_tube_radius_err =
          std::max(rep.max_tube_radius_err, std::abs(r_tube - r) / r);
    }
    if (std::abs(k_major) > real(1e-8) && std::abs(exact.k_major) > real(0.25)) {
      const real r_maj = real(1) / std::abs(k_major);
      const real r_exact = real(1) / std::abs(exact.k_major);
      rep.max_major_radius_err = std::max(
          rep.max_major_radius_err, std::abs(r_maj - r_exact) / r_exact);
    }
  }
}

inline bool land_guidance_and_charts(
    asawa::shell::shell::ptr M, hepworth::block::shell_cr_block &cr,
    std::vector<hepworth::block::vertex_chart> &charts) {
  if (!M)
    return false;
  hepworth::block::solver_context ctx;
  cr.prepare(ctx);
  if (cr.topology.nE <= 0 || static_cast<int>(cr.u.size()) < 2 * cr.topology.nE)
    return false;
  cr.commit_edge_tangents();
  const std::vector<vec3> &x = asawa::const_get_vec_data(*M, 0);
  const int nV = M->vert_count();
  std::vector<vec3> H(static_cast<size_t>(nV), vec3::Zero());
  std::vector<vec3> tang(static_cast<size_t>(nV), vec3::Zero());
  std::vector<vec3> axes(static_cast<size_t>(nV), vec3::Zero());
  std::vector<char> wrote(static_cast<size_t>(nV), 0);
  charts.clear();
  if (!hepworth::block::vertex_weingarten(*M, x, cr.u, H, tang, axes, wrote,
                                          &charts))
    return false;
  cr.vertex_frames = charts;
  return true;
}

inline torus_tangent_report report_torus_tangent_gates(
    asawa::shell::shell &M, hepworth::block::shell_cr_block &cr, real R,
    real r) {
  torus_tangent_report rep;
  const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);
  measure_one_ring_edge_tangent_coherence(M, x, cr.u, cr.vertex_frames, rep);
  measure_one_ring_labeled_principal_coherence(M, cr.vertex_frames, rep);
  measure_charts_vs_analytic_torus(M, x, cr.vertex_frames, R, r, rep);
  return rep;
}

/// Shared pass criteria: every vert ok, labeled neighbors aligned, analytic
/// angle + κ/radius. A 90° draw swap trips `n_e1_swap` / `min_e1_abs_dot`.
inline void expect_torus_frames_good(const torus_tangent_report &rep,
                                     real max_neighbor_deg, real max_axis_deg,
                                     real max_k_tube, real max_k_major,
                                     real max_r_tube, real max_r_major) {
  GAUDI_EXPECT(rep.n_vert > 0);
  GAUDI_EXPECT(rep.n_ok == rep.n_vert);
  GAUDI_EXPECT(rep.n_edge_empty == 0);
  GAUDI_EXPECT(rep.n_ring_pairs > 0);
  GAUDI_EXPECT(rep.n_e1_swap == 0);
  GAUDI_EXPECT(rep.n_e2_swap == 0);
  GAUDI_EXPECT(rep.min_e1_abs_dot > real(0.90));
  GAUDI_EXPECT(rep.min_e2_abs_dot > real(0.90));
  GAUDI_EXPECT(rep.max_e1_neighbor_deg < max_neighbor_deg);
  GAUDI_EXPECT(rep.max_e2_neighbor_deg < max_neighbor_deg);
  GAUDI_EXPECT(rep.max_tube_deg < max_axis_deg);
  GAUDI_EXPECT(rep.max_major_deg < max_axis_deg);
  GAUDI_EXPECT(rep.max_tight_tube_deg < max_axis_deg);
  GAUDI_EXPECT(rep.max_tube_k < max_k_tube);
  GAUDI_EXPECT(rep.max_major_k < max_k_major);
  GAUDI_EXPECT(rep.max_tube_radius_err < max_r_tube);
  GAUDI_EXPECT(rep.max_major_radius_err < max_r_major);
}

inline bool rebuild_charts_from_u(asawa::shell::shell &M,
                                  hepworth::block::shell_cr_block &cr) {
  const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);
  const int nV = M.vert_count();
  std::vector<vec3> H(static_cast<size_t>(nV), vec3::Zero());
  std::vector<vec3> tang(static_cast<size_t>(nV), vec3::Zero());
  std::vector<vec3> axes(static_cast<size_t>(nV), vec3::Zero());
  std::vector<char> wrote(static_cast<size_t>(nV), 0);
  std::vector<hepworth::block::vertex_chart> charts;
  if (!hepworth::block::vertex_weingarten(M, x, cr.u, H, tang, axes, wrote,
                                          &charts))
    return false;
  cr.vertex_frames = std::move(charts);
  return true;
}

/// Gate 1: full torus land — labeled e1/e2 neighbors + analytic angle/κ/radius.
GAUDI_TEST(torus_landed_tangents_coherent_vs_analytic) {
  const real R = 1.0;
  const real r = 0.33;
  asawa::shell::shell::ptr M = asawa::shell::load_torus(R, r, 48, 24);
  auto cr = std::make_shared<hepworth::block::shell_cr_block>(M);
  std::vector<hepworth::block::vertex_chart> charts;
  GAUDI_ASSERT(land_guidance_and_charts(M, *cr, charts));
  const torus_tangent_report rep = report_torus_tangent_gates(*M, *cr, R, r);
  GAUDI_EXPECT(rep.n_vert == 48 * 24);
  // Relative radius thresholds (fraction of analytic radius).
  expect_torus_frames_good(rep, /*neighbor*/ real(20), /*axis*/ real(8),
                           /*k*/ real(0.15), real(0.20),
                           /*radius*/ real(0.08), real(0.25));
}

/// Gate 2a: after subdiv + snap, carried tangents keep labeled e1/e2 coherent
/// and analytic angles. A 90° draw swap fails here (|e1·e1|≈0).
GAUDI_TEST(torus_remeshed_carried_axes_coherent_vs_analytic) {
  const real R = 1.0;
  const real r = 0.33;
  asawa::shell::shell::ptr M = asawa::shell::load_torus(R, r, 48, 24);
  auto cr = std::make_shared<hepworth::block::shell_cr_block>(M);
  std::vector<hepworth::block::vertex_chart> charts;
  GAUDI_ASSERT(land_guidance_and_charts(M, *cr, charts));

  const std::vector<vec3> &x0 = asawa::const_get_vec_data(*M, 0);
  const real l0 = asawa::shell::avg_length(*M, x0);
  const int nE0 = M->edge_count();
  auto dyn = asawa::shell::dynamic::create(M, real(0.05) * l0, real(0.55) * l0,
                                           real(0.25) * l0);
  dyn->subdivide_edges();
  GAUDI_EXPECT(M->edge_count() > nE0);
  GAUDI_ASSERT(M->verts_are_dense_packed());

  std::vector<vec3> &x = asawa::get_vec_data(*M, 0);
  snap_vertices_to_torus(x, R, r);
  hepworth::block::solver_context ctx;
  cr->prepare(ctx);
  GAUDI_ASSERT(rebuild_charts_from_u(*M, *cr));
  const torus_tangent_report rep = report_torus_tangent_gates(*M, *cr, R, r);
  GAUDI_EXPECT(rep.n_ok == rep.n_vert);
  GAUDI_EXPECT(rep.n_edge_empty == 0);
  GAUDI_EXPECT(rep.n_e1_swap == 0);
  GAUDI_EXPECT(rep.n_e2_swap == 0);
  GAUDI_EXPECT(rep.min_e1_abs_dot > real(0.90));
  GAUDI_EXPECT(rep.min_e2_abs_dot > real(0.90));
  GAUDI_EXPECT(rep.max_e1_neighbor_deg < real(25));
  GAUDI_EXPECT(rep.max_e2_neighbor_deg < real(25));
  GAUDI_EXPECT(rep.max_tube_deg < real(10));
  GAUDI_EXPECT(rep.max_major_deg < real(10));
  GAUDI_EXPECT(rep.max_tight_tube_deg < real(10));
}

/// Gate 2b: remeshed + snap + full Dirichlet rebake. Labeled axes and analytic
/// angles must stay land-quality. κ on irregular post-subdiv rings is weaker
/// than the regular grid (face Weingarten); still gated, but looser than land.
GAUDI_TEST(torus_remeshed_rebake_coherent_vs_analytic) {
  const real R = 1.0;
  const real r = 0.33;
  asawa::shell::shell::ptr M = asawa::shell::load_torus(R, r, 48, 24);
  auto cr = std::make_shared<hepworth::block::shell_cr_block>(M);
  std::vector<hepworth::block::vertex_chart> charts;
  GAUDI_ASSERT(land_guidance_and_charts(M, *cr, charts));

  const std::vector<vec3> &x0 = asawa::const_get_vec_data(*M, 0);
  const real l0 = asawa::shell::avg_length(*M, x0);
  auto dyn = asawa::shell::dynamic::create(M, real(0.05) * l0, real(0.55) * l0,
                                           real(0.25) * l0);
  dyn->subdivide_edges();
  GAUDI_ASSERT(M->verts_are_dense_packed());
  std::vector<vec3> &x = asawa::get_vec_data(*M, 0);
  snap_vertices_to_torus(x, R, r);

  GAUDI_ASSERT(hepworth::block::seed_cr_from_dirichlet(*M, *cr, /*reseed=*/true,
                                                       real(100), real(0)));
  const torus_tangent_report rep = report_torus_tangent_gates(*M, *cr, R, r);
  GAUDI_EXPECT(rep.n_ok == rep.n_vert);
  GAUDI_EXPECT(rep.n_edge_empty == 0);
  GAUDI_EXPECT(rep.n_e1_swap == 0);
  GAUDI_EXPECT(rep.n_e2_swap == 0);
  GAUDI_EXPECT(rep.min_e1_abs_dot > real(0.90));
  GAUDI_EXPECT(rep.min_e2_abs_dot > real(0.90));
  GAUDI_EXPECT(rep.max_e1_neighbor_deg < real(20));
  GAUDI_EXPECT(rep.max_e2_neighbor_deg < real(20));
  GAUDI_EXPECT(rep.max_tube_deg < real(8));
  GAUDI_EXPECT(rep.max_major_deg < real(8));
  GAUDI_EXPECT(rep.max_tight_tube_deg < real(8));
  // Face operators on irregular rings: document gap vs land (tube_k≈0.08).
  GAUDI_EXPECT(rep.max_tube_k < real(2.5));
  GAUDI_EXPECT(rep.max_major_k < real(1.5));
  GAUDI_EXPECT(rep.max_tube_radius_err < real(0.50));
  GAUDI_EXPECT(rep.max_major_radius_err < real(0.60));
}

/// Bucketed spine drive under anisotropy `ratio=0.7`: tube `r_tgt` vs osculating
/// `r_geom`, torus major `R_tgt`, and spoke hit gap. Crown (analytic κ_major~0)
/// still carries a tiny measured major κ, so ratio fires and R_tgt / r_tgt
/// swing hard while local hit gaps stay small (large canal ≈ planar on ring).
GAUDI_TEST(spine_anisotropy_ratio_metrics_on_torus) {
  const real R = 1.0;
  const real r = 0.33;
  const real ratio = 0.7;
  asawa::shell::shell::ptr M = asawa::shell::load_torus(R, r, 48, 24);
  auto cr = std::make_shared<hepworth::block::shell_cr_block>(M);
  std::vector<hepworth::block::vertex_chart> charts;
  GAUDI_ASSERT(land_guidance_and_charts(M, *cr, charts));
  const std::vector<vec3> &x = asawa::const_get_vec_data(*M, 0);

  struct bucket {
    const char *name = "";
    int n = 0;
    int n_ratio = 0;
    int n_ratio_skip = 0;
    real sum_r_gap0 = 0;
    real sum_r_gap1 = 0;
    real sum_R0 = 0;
    real sum_R1 = 0;
    real sum_dR = 0;
    real sum_dr = 0;
    real sum_hit0 = 0;
    real sum_hit1 = 0;
    real max_hit1 = 0;
    real max_dR = 0;
    real max_dr = 0;
  };
  bucket outer{"outer_equator"};
  bucket crown{"crown"};
  bucket inner{"inner_equator"};
  bucket mid{"mid_latitude"};
  bucket all{"all"};

  auto pick = [&](real k_maj) -> bucket & {
    const real a = std::abs(k_maj);
    if (a < real(0.15))
      return crown;
    if (k_maj > real(0.45))
      return outer;
    if (k_maj < real(-0.45))
      return inner;
    return mid;
  };

  auto accumulate = [](bucket &b, bool ratio_ok, const hepworth::block::spine_ring_metric &m0,
                       const hepworth::block::spine_ring_metric &m1) {
    ++b.n;
    if (ratio_ok)
      ++b.n_ratio;
    else
      ++b.n_ratio_skip;
    if (!m0.ok || !m1.ok)
      return;
    b.sum_r_gap0 += m0.r_gap;
    b.sum_r_gap1 += m1.r_gap;
    b.sum_R0 += m0.R_tgt;
    b.sum_R1 += m1.R_tgt;
    const real dR = std::abs(m1.R_tgt - m0.R_tgt);
    const real dr = std::abs(m1.r_tgt - m0.r_tgt);
    b.sum_dR += dR;
    b.sum_dr += dr;
    b.max_dR = std::max(b.max_dR, dR);
    b.max_dr = std::max(b.max_dr, dr);
    b.sum_hit0 += m0.mean_hit_gap;
    b.sum_hit1 += m1.mean_hit_gap;
    b.max_hit1 = std::max(b.max_hit1, m1.max_hit_gap);
  };

  auto dump = [](const bucket &b) {
    if (b.n == 0)
      return;
    const real inv = real(1) / real(b.n);
    console_logger::info << "spine_ratio0.7 " << b.name << " n=" << b.n
                         << " ratio_ok=" << b.n_ratio
                         << " ratio_skip=" << b.n_ratio_skip
                         << " mean_r_gap0=" << (b.sum_r_gap0 * inv)
                         << " mean_r_gap1=" << (b.sum_r_gap1 * inv)
                         << " mean_R0=" << (b.sum_R0 * inv)
                         << " mean_R1=" << (b.sum_R1 * inv)
                         << " mean_|dR|=" << (b.sum_dR * inv)
                         << " max_|dR|=" << b.max_dR
                         << " mean_|dr_tube|=" << (b.sum_dr * inv)
                         << " max_|dr_tube|=" << b.max_dr
                         << " mean_hit0=" << (b.sum_hit0 * inv)
                         << " mean_hit1=" << (b.sum_hit1 * inv)
                         << " max_hit1=" << b.max_hit1 << std::endl;
  };

  for (int vi = 0; vi < M->vert_count(); ++vi) {
    hepworth::block::shape_w W;
    vec3 T, B, N;
    if (!hepworth::block::vertex_chart_at(&charts, vi, W, T, B, N))
      continue;
    const hepworth::block::shape_w axes = hepworth::block::shape_with_axes(W);
    if (!axes.ok)
      continue;
    const hepworth::block::monge_axes ax0 = hepworth::block::monge_axes_of(axes);
    hepworth::block::spine_curvature cr0;
    cr0.tube = ax0.hx * T + ax0.hy * B;
    cr0.ruling = ax0.lx * T + ax0.ly * B;
    cr0.normal = N;
    cr0.k_hi = ax0.k_hi;
    cr0.k_lo = ax0.k_lo;

    std::vector<vec3> project_spokes;
    const std::vector<asawa::shell::VertId> ring =
        M->get_one_ring(asawa::shell::vert_id(vi));
    for (const asawa::shell::VertId nj : ring) {
      const int j = static_cast<int>(nj);
      if (j < 0 || j == vi)
        continue;
      project_spokes.push_back(x[static_cast<size_t>(j)] -
                               x[static_cast<size_t>(vi)]);
    }
    if (project_spokes.size() < 3)
      continue;

    // Same as frame_fair: 2-ring canal fit, 1-ring project.
    std::vector<vec3> fit_spokes;
    for (const index_t j : hepworth::block::vertex_two_ring_ids(*M, vi)) {
      if (j < 0 || j >= static_cast<int>(x.size()))
        continue;
      fit_spokes.push_back(x[static_cast<size_t>(j)] -
                           x[static_cast<size_t>(vi)]);
    }
    const hepworth::block::spine_set kinds =
        hepworth::block::all_spine_kinds();
    const hepworth::block::spine_ring_metric m0 =
        hepworth::block::measure_spine_ring(cr0, project_spokes, kinds,
                                           fit_spokes);

    real k1 = axes.k1;
    real k2 = axes.k2;
    const bool ratio_ok =
        hepworth::block::absolute_ratio_pair(axes.k1, axes.k2, ratio, k1, k2);
    hepworth::block::shape_w Wan = axes;
    if (ratio_ok)
      Wan = hepworth::block::shape_from_principals(k1, k2, axes.e1x, axes.e1y,
                                                   axes.e2x, axes.e2y);
    const hepworth::block::monge_axes ax1 =
        hepworth::block::monge_axes_of(Wan.ok ? Wan : axes);
    hepworth::block::spine_curvature cr1 = cr0;
    cr1.tube = ax1.hx * T + ax1.hy * B;
    cr1.ruling = ax1.lx * T + ax1.ly * B;
    cr1.k_hi = ax1.k_hi;
    cr1.k_lo = ax1.k_lo;
    const hepworth::block::spine_ring_metric m1 =
        hepworth::block::measure_spine_ring(cr1, project_spokes, kinds,
                                           fit_spokes);

    const analytic_torus exact = torus_at(R, r, x[static_cast<size_t>(vi)]);
    accumulate(all, ratio_ok, m0, m1);
    accumulate(pick(exact.k_major), ratio_ok, m0, m1);
  }

  dump(all);
  dump(outer);
  dump(crown);
  dump(inner);
  dump(mid);

  GAUDI_EXPECT(all.n > 0);
  GAUDI_EXPECT(outer.n > 0);
  GAUDI_EXPECT(crown.n > 0);
  // Outer equator: both κ live ⇒ ratio applies; major canal shrinks hard.
  GAUDI_EXPECT(outer.n_ratio == outer.n);
  GAUDI_EXPECT(outer.max_dR > real(0.3));
  GAUDI_EXPECT(outer.sum_dR / real(outer.n) > real(0.2));
  GAUDI_EXPECT(outer.sum_dr / real(outer.n) > real(0.1));
  // Held centerpoint + larger minor tube ⇒ O(|Δr|) spoke push.
  GAUDI_EXPECT(outer.max_hit1 > real(0.05));
  GAUDI_EXPECT(outer.sum_hit1 / real(outer.n) > real(0.05));
  GAUDI_EXPECT(crown.n_ratio == crown.n);
  GAUDI_EXPECT(crown.max_dR > real(10));
  GAUDI_EXPECT(crown.sum_dr / real(crown.n) > real(1));
  GAUDI_EXPECT(crown.max_hit1 > real(0.5));
}

} // namespace test
} // namespace gaudi

#endif
