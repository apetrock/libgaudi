#ifndef GAUDI_TEST_LIE_SPHERE_TESTS_HPP
#define GAUDI_TEST_LIE_SPHERE_TESTS_HPP

#include "gaudi/albers/lie_dupin.hpp"
#include "gaudi/albers/lie_sphere.hpp"
#include "gaudi/test/test.hpp"

#include <Eigen/Geometry>
#include <cmath>
#include <random>
#include <vector>

namespace gaudi {
namespace test {
namespace lie_sphere_detail {

using namespace albers;

/// Point and outward normal of the canonical torus at spine angle `th` and
/// tube angle `ph`.
inline void torus_contact(real R, real r, real th, real ph, vec3 &p, vec3 &n) {
  n = vec3(std::cos(ph) * std::cos(th), std::cos(ph) * std::sin(th),
           std::sin(ph));
  p = vec3(R * std::cos(th), R * std::sin(th), 0) + r * n;
}

/// Pushes the contact element `(p, n)` through `A` and decodes it.
inline bool push_contact(const mat6 &A, vec3 &p, vec3 &n) {
  const lie_contact_element L = lie_contact(p, n);
  return lie_contact_decode(A * L.point, A * L.plane, p, n);
}

inline real sigma_distance(const lie_dupin &a, const lie_dupin &b) {
  return (a.sigma - b.sigma).cwiseAbs().maxCoeff();
}

/// Distance up to the sign that swaps the two curvature families.
inline real unsigned_sigma_distance(const lie_dupin &a, const lie_dupin &b) {
  return std::min((a.sigma - b.sigma).cwiseAbs().maxCoeff(),
                  (a.sigma + b.sigma).cwiseAbs().maxCoeff());
}

/// A Möbius-and-rigid map that turns the canonical torus into a generic Dupin
/// cyclide: inversion in a sphere off the torus, then a rigid motion.
inline mat6 cyclide_map() {
  const mat3 R =
      Eigen::AngleAxis<real>(0.7, vec3(1, 2, -0.5).normalized()).toRotationMatrix();
  return lie_rigid(R, vec3(0.3, -0.2, 0.5)) *
         lie_reflection(lie_inversion_sphere(vec3(0.4, 2.2, 1.1), 1.5));
}

} // namespace lie_sphere_detail

GAUDI_TEST(lie_sphere_vectors_are_null_and_decode) {
  using namespace albers;
  const vec3 c(0.3, -1.2, 2.0);
  const vec6 s = lie_sphere(c, -0.7);
  GAUDI_EXPECT(std::abs(lie_dot(s, s)) < 1e-12);
  const lie_sphere_decoded ds = lie_decode(3.5 * s);
  GAUDI_EXPECT(ds.type == lie_sphere_decoded::kind::sphere);
  GAUDI_EXPECT((ds.c - c).norm() < 1e-12);
  GAUDI_EXPECT(std::abs(ds.r + 0.7) < 1e-12);

  const vec3 n = vec3(1, -2, 0.5).normalized();
  const vec6 t = lie_plane(n, 0.8);
  GAUDI_EXPECT(std::abs(lie_dot(t, t)) < 1e-12);
  const lie_sphere_decoded dt = lie_decode(-2.0 * t);
  GAUDI_EXPECT(dt.type == lie_sphere_decoded::kind::plane);
  GAUDI_EXPECT((dt.n - n).norm() < 1e-12);
  GAUDI_EXPECT(std::abs(dt.h - 0.8) < 1e-12);

  GAUDI_EXPECT(std::abs(lie_dot(lie_infinity(), lie_infinity())) < 1e-15);
  GAUDI_EXPECT(std::abs(lie_dot(lie_infinity(), t)) < 1e-15);
  GAUDI_EXPECT(lie_decode(lie_infinity()).type ==
               lie_sphere_decoded::kind::infinity);
}

GAUDI_TEST(lie_sphere_contact_is_orthogonality) {
  using namespace albers;
  // <s, t> = (|r1 - r2|^2 - |c1 - c2|^2) / 2.
  const vec3 c1(0.1, 0.2, 0.3), c2(1.0, -0.5, 0.7);
  const real r1 = 0.4, r2 = -0.9;
  const real expect = 0.5 * ((r1 - r2) * (r1 - r2) - (c1 - c2).squaredNorm());
  GAUDI_EXPECT(std::abs(lie_dot(lie_sphere(c1, r1), lie_sphere(c2, r2)) - expect) <
               1e-12);

  // Internally tangent, same orientation.
  const vec3 d = vec3(1, 1, -1).normalized();
  GAUDI_EXPECT(std::abs(lie_dot(lie_sphere(c1, 2.0), lie_sphere(c1 + 0.5 * d, 1.5))) <
               1e-12);
  // Point on a sphere, point on a plane, sphere touching a plane at
  // x = c - r n.
  const vec3 n = vec3(0.2, -0.3, 0.9).normalized();
  GAUDI_EXPECT(std::abs(lie_dot(lie_point(c1 + 0.8 * d), lie_sphere(c1, -0.8))) < 1e-12);
  GAUDI_EXPECT(std::abs(lie_dot(lie_point(c2), lie_plane(n, n.dot(c2)))) < 1e-12);
  GAUDI_EXPECT(std::abs(lie_dot(lie_sphere(c2 + 0.6 * n, 0.6), lie_plane(n, n.dot(c2)))) <
               1e-12);
  // The opposite orientation does not touch.
  GAUDI_EXPECT(std::abs(lie_dot(lie_sphere(c2 + 0.6 * n, -0.6), lie_plane(n, n.dot(c2)))) >
               0.1);
}

GAUDI_TEST(lie_sphere_curvature_spheres_are_the_contact_pencil) {
  using namespace albers;
  const vec3 p(0.4, -0.1, 1.3);
  const vec3 n = vec3(-0.3, 0.8, 0.5).normalized();
  const vec6 s = lie_curvature_sphere(p, n, 2.5);
  const lie_sphere_decoded ds = lie_decode(s);
  GAUDI_EXPECT(ds.type == lie_sphere_decoded::kind::sphere);
  GAUDI_EXPECT((ds.c - (p + n / 2.5)).norm() < 1e-12);
  GAUDI_EXPECT(std::abs(ds.r - 1.0 / 2.5) < 1e-12);
  // k = 0 is the tangent plane.
  const lie_sphere_decoded d0 = lie_decode(lie_curvature_sphere(p, n, 0.0));
  GAUDI_EXPECT(d0.type == lie_sphere_decoded::kind::plane);
  GAUDI_EXPECT((d0.n - n).norm() < 1e-12);

  // Every member is null and touches every other: the pencil is a null line.
  const vec6 a = lie_curvature_sphere(p, n, -1.7);
  GAUDI_EXPECT(std::abs(lie_dot(a, a)) < 1e-12);
  GAUDI_EXPECT(std::abs(lie_dot(a, s)) < 1e-12);

  // Any two curvature spheres give the contact element back.
  vec3 q, m;
  GAUDI_EXPECT(lie_contact_decode(a, s, q, m));
  GAUDI_EXPECT((q - p).norm() < 1e-10);
  GAUDI_EXPECT((m - n).norm() < 1e-10);
}

GAUDI_TEST(lie_sphere_transforms_are_in_the_group) {
  using namespace albers;
  const mat3 R =
      Eigen::AngleAxis<real>(1.1, vec3(0.3, -1, 0.4).normalized()).toRotationMatrix();
  const std::vector<mat6> A = {
      lie_rotation(R),
      lie_translation(vec3(0.5, -1.5, 2.0)),
      lie_rigid(R, vec3(-0.2, 0.1, 0.9)),
      lie_scale(1.7),
      lie_normal_offset(0.35),
      lie_orientation_flip(),
      lie_reflection(lie_mirror_plane(vec3(0, 0.6, 0.8), 0.25)),
      lie_reflection(lie_inversion_sphere(vec3(1, 2, 3), 0.75)),
  };
  for (const mat6 &a : A) {
    // Round-off grows with the entries (an inversion far from the origin has
    // entries in the hundreds).
    const real tol = 1e-14 * std::max(real(1), a.squaredNorm());
    GAUDI_EXPECT(lie_group_defect(a) < tol);
    GAUDI_EXPECT((lie_inverse(a) * a - mat6::Identity()).cwiseAbs().maxCoeff() < tol);
  }

  vec15 w;
  for (int k = 0; k < 15; ++k)
    w[k] = 0.1 * std::sin(1.0 + 2.3 * k);
  const mat6 X = lie_algebra(w);
  GAUDI_EXPECT(lie_algebra_defect(X) < 1e-15);
  GAUDI_EXPECT((lie_algebra_coords(X) - w).norm() < 1e-15);
  GAUDI_EXPECT(lie_group_defect(lie_exp(X)) < 1e-12);
  GAUDI_EXPECT(lie_group_defect(lie_cayley(X)) < 1e-12);
  GAUDI_EXPECT((lie_log(lie_exp(X)) - X).cwiseAbs().maxCoeff() < 1e-10);
  // Rotations are positive under the invariant form, boosts negative.
  GAUDI_EXPECT(lie_algebra_form(lie_algebra_generator(2, 3)) > 0);
  GAUDI_EXPECT(lie_algebra_form(lie_algebra_generator(0, 2)) < 0);
}

GAUDI_TEST(lie_sphere_transforms_act_as_named) {
  using namespace albers;
  const vec3 c(0.3, 0.4, -0.5);
  const real r = 0.6;
  const mat3 R = Eigen::AngleAxis<real>(0.9, vec3::UnitY()).toRotationMatrix();
  const vec3 t(1, -2, 0.5);

  lie_sphere_decoded d = lie_decode(lie_rigid(R, t) * lie_sphere(c, r));
  GAUDI_EXPECT((d.c - (R * c + t)).norm() < 1e-12);
  GAUDI_EXPECT(std::abs(d.r - r) < 1e-12);

  d = lie_decode(lie_scale(2.5) * lie_sphere(c, r));
  GAUDI_EXPECT((d.c - 2.5 * c).norm() < 1e-12);
  GAUDI_EXPECT(std::abs(d.r - 2.5 * r) < 1e-12);

  d = lie_decode(lie_normal_offset(0.2) * lie_sphere(c, r));
  GAUDI_EXPECT((d.c - c).norm() < 1e-12);
  GAUDI_EXPECT(std::abs(d.r - (r - 0.2)) < 1e-12);

  const vec3 n = vec3(0, 0.6, 0.8);
  d = lie_decode(lie_normal_offset(0.3) * lie_plane(n, 1.0));
  GAUDI_EXPECT(d.type == lie_sphere_decoded::kind::plane);
  GAUDI_EXPECT(std::abs(d.h - 1.3) < 1e-12);

  // Mirror in the plane n.x = h.
  const vec3 x(0.7, -0.4, 1.9);
  d = lie_decode(lie_reflection(lie_mirror_plane(n, 0.25)) * lie_point(x));
  GAUDI_EXPECT((d.c - (x - 2.0 * (n.dot(x) - 0.25) * n)).norm() < 1e-12);

  // Inversion in the sphere (c, rho): x -> c + rho^2 (x - c) / |x - c|^2.
  const vec3 ic(0.1, 0.2, -0.3);
  const real rho = 0.8;
  d = lie_decode(lie_reflection(lie_inversion_sphere(ic, rho)) * lie_point(x));
  const vec3 xi = ic + rho * rho * (x - ic) / (x - ic).squaredNorm();
  GAUDI_EXPECT((d.c - xi).norm() < 1e-12);
  GAUDI_EXPECT(std::abs(d.r) < 1e-12);

  // A rigid motion moves contact elements: (p, n) -> (R p + t, R n).
  vec3 p = x, m = n;
  GAUDI_EXPECT(lie_sphere_detail::push_contact(lie_rigid(R, t), p, m));
  GAUDI_EXPECT((p - (R * x + t)).norm() < 1e-12);
  GAUDI_EXPECT((m - R * n).norm() < 1e-12);
  // The offset slides it along its own normal.
  p = x;
  m = n;
  GAUDI_EXPECT(lie_sphere_detail::push_contact(lie_normal_offset(0.4), p, m));
  GAUDI_EXPECT((p - (x + 0.4 * n)).norm() < 1e-12);
  GAUDI_EXPECT((m - n).norm() < 1e-12);
}

GAUDI_TEST(lie_dupin_torus_is_a_cyclide_coordinate) {
  using namespace albers;
  const lie_dupin d = lie_dupin_torus(1.3, 0.4);
  GAUDI_EXPECT(lie_dupin_defect(d) < 1e-12);
  GAUDI_EXPECT(lie_group_defect(d.sigma) < 1e-12);
  GAUDI_EXPECT(lie_dupin_signature(d, 1) == Eigen::Vector2i(2, 1));
  GAUDI_EXPECT(lie_dupin_signature(d, -1) == Eigen::Vector2i(2, 1));
  // A spindle torus (r > R) is one too.
  const lie_dupin spindle = lie_dupin_torus(0.3, 0.5);
  GAUDI_EXPECT(lie_dupin_defect(spindle) < 1e-12);
  GAUDI_EXPECT(lie_dupin_signature(spindle, 1) == Eigen::Vector2i(2, 1));
}

GAUDI_TEST(lie_dupin_torus_contains_its_contact_elements) {
  using namespace albers;
  const real R = 1.3, r = 0.4;
  const lie_dupin d = lie_dupin_torus(R, r);
  real on = 0, off_point = 1e9, off_normal = 1e9;
  for (int i = 0; i < 12; ++i)
    for (int j = 0; j < 9; ++j) {
      vec3 p, n;
      lie_sphere_detail::torus_contact(R, r, 0.5 * i, 0.7 * j + 0.1, p, n);
      on = std::max(on, lie_dupin_contact_residual(d, p, n));
      off_point = std::min(off_point, lie_dupin_contact_residual(d, p + 0.05 * n, n));
      const vec3 tilt = (n + 0.1 * n.cross(vec3(0.3, 0.1, 1)).normalized()).normalized();
      off_normal = std::min(off_normal, lie_dupin_contact_residual(d, p, tilt));
    }
  GAUDI_EXPECT(on < 1e-12);
  GAUDI_EXPECT(off_point > 1e-3);
  GAUDI_EXPECT(off_normal > 1e-3);
  // Reversing the normal leaves the surface (orientation is part of it).
  vec3 p, n;
  lie_sphere_detail::torus_contact(R, r, 0.3, 0.9, p, n);
  GAUDI_EXPECT(lie_dupin_contact_residual(d, p, -n) > 1e-2);
  GAUDI_EXPECT(lie_dupin_contact_residual(lie_dupin_transform(lie_orientation_flip(), d),
                                          p, -n) < 1e-12);
}

GAUDI_TEST(lie_dupin_torus_curvatures_are_its_principal_curvatures) {
  using namespace albers;
  const real R = 1.3, r = 0.4;
  const lie_dupin d = lie_dupin_torus(R, r);
  for (real ph : {0.0, 0.6, 1.4, 2.5, M_PI}) {
    vec3 p, n;
    lie_sphere_detail::torus_contact(R, r, 0.8, ph, p, n);
    const lie_dupin_curvature_pair k = lie_dupin_curvatures(d, p, n);
    GAUDI_EXPECT(std::abs(k.k_u + 1.0 / r) < 1e-10);
    GAUDI_EXPECT(std::abs(k.k_v + std::cos(ph) / (R + r * std::cos(ph))) < 1e-10);
    GAUDI_EXPECT(k.res_u < 1e-12);
    GAUDI_EXPECT(k.res_v < 1e-12);
  }
  // On top of the tube the parallel is a flat circle's tangent: k_v = 0, with
  // no special case.
  vec3 p, n;
  lie_sphere_detail::torus_contact(R, r, 0.2, 0.5 * M_PI, p, n);
  GAUDI_EXPECT(std::abs(lie_dupin_curvatures(d, p, n).k_v) < 1e-10);
}

GAUDI_TEST(lie_dupin_advects_by_conjugation) {
  using namespace albers;
  const real R = 1.1, r = 0.3;
  const vec3 c(0.5, -0.2, 0.8);
  const vec3 axis = vec3(0.3, 0.4, 1).normalized();
  const lie_dupin moved = lie_dupin_torus(c, axis, R, r);
  GAUDI_EXPECT(lie_dupin_defect(moved) < 1e-12);

  // Pushing the canonical torus's contact elements through the same motion
  // lands them on the moved torus.
  const mat3 Rz = Eigen::Quaternion<real>::FromTwoVectors(vec3::UnitZ(), axis)
                      .toRotationMatrix();
  real worst = 0;
  for (int i = 0; i < 8; ++i) {
    vec3 p, n;
    lie_sphere_detail::torus_contact(R, r, 0.8 * i, 0.5 + 0.6 * i, p, n);
    worst = std::max(worst, lie_dupin_contact_residual(moved, Rz * p + c, Rz * n));
  }
  GAUDI_EXPECT(worst < 1e-12);

  // A Laguerre offset by d is the torus with tube r + d.
  const lie_dupin grown =
      lie_dupin_transform(lie_normal_offset(0.15), lie_dupin_torus(R, r));
  GAUDI_EXPECT(lie_sphere_detail::sigma_distance(grown, lie_dupin_torus(R, r + 0.15)) <
               1e-12);
  // So is a uniform scale.
  const lie_dupin scaled = lie_dupin_transform(lie_scale(2.0), lie_dupin_torus(R, r));
  GAUDI_EXPECT(lie_sphere_detail::sigma_distance(scaled, lie_dupin_torus(2 * R, 2 * r)) <
               1e-12);
}

GAUDI_TEST(lie_dupin_mobius_image_is_a_cyclide) {
  using namespace albers;
  const real R = 1.3, r = 0.4;
  const mat6 A = lie_sphere_detail::cyclide_map();
  const lie_dupin cyc = lie_dupin_transform(A, lie_dupin_torus(R, r));
  GAUDI_EXPECT(lie_dupin_defect(cyc) < 1e-10);
  GAUDI_EXPECT(lie_dupin_signature(cyc, 1) == Eigen::Vector2i(2, 1));

  real worst = 0;
  real k_spread = 0;
  real k_first = 0;
  for (int i = 0; i < 10; ++i) {
    vec3 p, n;
    lie_sphere_detail::torus_contact(R, r, 0.6 * i, 0.4 + 0.5 * i, p, n);
    GAUDI_EXPECT(lie_sphere_detail::push_contact(A, p, n));
    worst = std::max(worst, lie_dupin_contact_residual(cyc, p, n));
    const lie_dupin_curvature_pair k = lie_dupin_curvatures(cyc, p, n);
    GAUDI_EXPECT(k.res_u < 1e-10);
    GAUDI_EXPECT(k.res_v < 1e-10);
    if (i == 0)
      k_first = k.k_u;
    k_spread = std::max(k_spread, std::abs(k.k_u - k_first));
  }
  GAUDI_EXPECT(worst < 1e-10);
  // Not a torus: the tube curvature varies around the surface.
  GAUDI_EXPECT(k_spread > 1e-2);
}

GAUDI_TEST(lie_dupin_fit_recovers_torus_and_cyclide) {
  using namespace albers;
  const real R = 1.3, r = 0.4;
  const mat6 A = lie_sphere_detail::cyclide_map();
  const lie_dupin torus = lie_dupin_torus(R, r);
  const lie_dupin cyc = lie_dupin_transform(A, torus);

  for (int which = 0; which < 2; ++which) {
    const lie_dupin &truth = which == 0 ? torus : cyc;
    lie_dupin_fit fit;
    for (int i = 0; i < 12; ++i) {
      vec3 p, n;
      lie_sphere_detail::torus_contact(R, r, 0.53 * i, 0.31 + 0.77 * i, p, n);
      if (which == 1)
        GAUDI_EXPECT(lie_sphere_detail::push_contact(A, p, n));
      fit.accumulate(1.0, p, n);
    }
    const lie_dupin_fit::result res = fit.solve();
    GAUDI_EXPECT(res.ok);
    GAUDI_EXPECT(res.gap < 1e-8);
    GAUDI_EXPECT(lie_sphere_detail::unsigned_sigma_distance(res.cyclide, truth) < 1e-6);
  }

  // A local patch (a one-ring sized neighbourhood) still pins the torus.
  lie_dupin_fit local;
  for (int i = 0; i < 7; ++i) {
    vec3 p, n;
    lie_sphere_detail::torus_contact(R, r, 0.4 + 0.08 * std::cos(0.9 * i),
                                     1.0 + 0.08 * std::sin(0.9 * i), p, n);
    local.accumulate(1.0, p, n);
  }
  const lie_dupin_fit::result lres = local.solve(1e-2);
  GAUDI_EXPECT(lie_sphere_detail::unsigned_sigma_distance(lres.cyclide, torus) < 1e-4);

  // Three contact elements leave a pencil: no unique cyclide.
  lie_dupin_fit three;
  for (int i = 0; i < 3; ++i) {
    vec3 p, n;
    lie_sphere_detail::torus_contact(R, r, 1.1 * i, 0.4 + 1.3 * i, p, n);
    three.accumulate(1.0, p, n);
  }
  GAUDI_EXPECT(!three.solve().ok);
}

GAUDI_TEST(lie_dupin_normalize_projects_onto_cyclides) {
  using namespace albers;
  const lie_dupin truth = lie_dupin_torus(1.3, 0.4);
  mat6 noise;
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j)
      noise(i, j) = 1e-3 * std::sin(1.0 + 7.0 * i + 3.0 * j);
  const lie_dupin d = lie_dupin_normalize(truth.sigma + noise);
  GAUDI_EXPECT(lie_dupin_defect(d) < 1e-12);
  GAUDI_EXPECT(lie_dupin_signature(d, 1) == Eigen::Vector2i(2, 1));
  GAUDI_EXPECT(lie_sphere_detail::sigma_distance(d, truth) < 1e-2);
  // A cyclide is a fixed point.
  GAUDI_EXPECT(lie_sphere_detail::sigma_distance(lie_dupin_normalize(truth.sigma), truth) <
               1e-14);
}

GAUDI_TEST(lie_dupin_transvection_log_and_geodesic) {
  using namespace albers;
  const lie_dupin a = lie_dupin_torus(1.3, 0.4);
  vec15 w;
  for (int k = 0; k < 15; ++k)
    w[k] = 0.15 * std::cos(0.5 + 1.7 * k);
  const lie_dupin b = lie_dupin_transform(lie_exp(lie_algebra(w)), a);

  const mat6 g = lie_dupin_transvection(a, b);
  GAUDI_EXPECT(lie_group_defect(g) < 1e-10);
  GAUDI_EXPECT(lie_sphere_detail::sigma_distance(lie_dupin_transform(g, a), b) < 1e-10);

  const mat6 xi = lie_dupin_log(a, b);
  GAUDI_EXPECT(lie_algebra_defect(xi) < 1e-10);
  // A tangent vector anticommutes with the base point.
  GAUDI_EXPECT((a.sigma * xi + xi * a.sigma).cwiseAbs().maxCoeff() < 1e-10);

  GAUDI_EXPECT(lie_sphere_detail::sigma_distance(lie_dupin_geodesic(a, b, 0), a) < 1e-12);
  GAUDI_EXPECT(lie_sphere_detail::sigma_distance(lie_dupin_geodesic(a, b, 1), b) < 1e-10);
  const lie_dupin mid = lie_dupin_geodesic(a, b, 0.5);
  GAUDI_EXPECT(lie_dupin_defect(mid) < 1e-10);
  // The midpoint is equidistant: the geodesic symmetry at mid swaps a and b.
  const lie_dupin swapped = lie_dupin_transform(mid.sigma, a);
  GAUDI_EXPECT(lie_sphere_detail::sigma_distance(swapped, b) < 1e-9);
}

GAUDI_TEST(lie_dupin_holonomy_fixes_the_base_cyclide) {
  using namespace albers;
  const lie_dupin a = lie_dupin_torus(1.3, 0.4);
  const lie_dupin b = lie_dupin_torus(vec3(0.1, 0, 0), vec3(0, 0.1, 1).normalized(), 1.25, 0.42);
  const lie_dupin c = lie_dupin_torus(vec3(0, -0.1, 0.05), vec3(0.1, 0, 1).normalized(), 1.35, 0.38);
  const mat6 H = lie_holonomy(lie_dupin_transvection(a, b), lie_dupin_transvection(b, c),
                              lie_dupin_transvection(c, a));
  GAUDI_EXPECT(lie_group_defect(H) < 1e-10);
  GAUDI_EXPECT(lie_sphere_detail::sigma_distance(lie_dupin_transform(H, a), a) < 1e-10);
  // A generic triangle has holonomy; three equal cyclides have none.
  GAUDI_EXPECT((H - mat6::Identity()).cwiseAbs().maxCoeff() > 1e-6);
  const mat6 H0 = lie_holonomy(lie_dupin_transvection(a, a), lie_dupin_transvection(a, a),
                               lie_dupin_transvection(a, a));
  GAUDI_EXPECT((H0 - mat6::Identity()).cwiseAbs().maxCoeff() < 1e-12);
}

GAUDI_TEST(lie_dupin_tangent_space_is_nine_dimensional) {
  using namespace albers;
  const lie_dupin a = lie_dupin_transform(lie_sphere_detail::cyclide_map(),
                                          lie_dupin_torus(1.3, 0.4));
  const std::vector<mat6> basis = lie_dupin_tangent_basis(a);
  GAUDI_EXPECT(basis.size() == 9);
  for (const mat6 &X : basis) {
    GAUDI_EXPECT(lie_algebra_defect(X) < 1e-12);
    GAUDI_EXPECT((a.sigma * X + X * a.sigma).cwiseAbs().maxCoeff() < 1e-10);
  }
}

GAUDI_TEST(lie_dupin_correct_moves_minimally_onto_new_data) {
  using namespace albers;
  const real R = 1.3, r = 0.4;
  const lie_dupin truth = lie_dupin_torus(R, r);
  std::vector<vec3> ps, ns;
  for (int i = 0; i < 10; ++i) {
    vec3 p, n;
    lie_sphere_detail::torus_contact(R, r, 0.61 * i, 0.27 + 0.83 * i, p, n);
    ps.push_back(p);
    ns.push_back(n);
  }
  const std::vector<real> w(ps.size(), 1.0);

  // Start from a nearby cyclide, as an advected state would be.
  vec15 v;
  for (int k = 0; k < 15; ++k)
    v[k] = 0.02 * std::sin(2.0 + 1.3 * k);
  const lie_dupin start = lie_dupin_transform(lie_exp(lie_algebra(v)), truth);
  GAUDI_EXPECT(lie_sphere_detail::sigma_distance(start, truth) > 1e-3);

  const lie_dupin fixed = lie_dupin_correct(start, ps, ns, w, 1e-12, 6);
  GAUDI_EXPECT(lie_dupin_defect(fixed) < 1e-10);
  GAUDI_EXPECT(lie_sphere_detail::sigma_distance(fixed, truth) < 1e-7);
  real worst = 0;
  for (size_t k = 0; k < ps.size(); ++k)
    worst = std::max(worst, lie_dupin_contact_residual(fixed, ps[k], ns[k]));
  GAUDI_EXPECT(worst < 1e-8);

  // A heavy prior holds the advected state.
  const lie_dupin held = lie_dupin_correct(start, ps, ns, w, 1e6, 2);
  GAUDI_EXPECT(lie_sphere_detail::sigma_distance(held, start) < 1e-5);
}

} // namespace test
} // namespace gaudi

#endif
