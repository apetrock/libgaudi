#ifndef GAUDI_ALBERS_LIE_SPHERE_HPP
#define GAUDI_ALBERS_LIE_SPHERE_HPP

#include "gaudi/common.h"

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <algorithm>
#include <cmath>

/// Lie sphere geometry: oriented spheres, planes and points of R^3 as null
/// vectors of R^{4,2}, and the Lie sphere group O(4,2) that acts on them.
///
/// The inner product is
///
///     <x, y> = -x0 y0 + x1 y1 + x2 y2 + x3 y3 + x4 y4 - x5 y5
///
/// (`lie_metric()`, written eta). An oriented sphere with centre c and signed
/// radius r is
///
///     s = [ (1 + |c|^2 - r^2) / 2,  (1 - |c|^2 + r^2) / 2,  c,  r ],
///
/// a point is a sphere of radius 0 and an oriented plane n.x = h is
/// [h, -h, n, 1]. All are null (<s, s> = 0), and two of them are in oriented
/// contact exactly when <s, t> = 0. Coordinates are projective: s and a s are
/// the same sphere.
///
/// Orientation: a sphere of radius r touches the oriented plane with unit
/// normal n at x when its centre is x + r n. So r > 0 means the normal points
/// at the centre, and the curvature sphere of normal curvature k at a surface
/// point (p, n) has centre p + n / k and radius 1 / k (k > 0 when the surface
/// bends towards n, the sign of the monge sag `h = k t^2 / 2`).
///
/// A contact element (a point with a unit normal) is the null 2-plane spanned
/// by its point sphere P and its tangent plane T. Every sphere tangent at that
/// point lies in it, with no division: the sphere of curvature k is
/// `k P + T` (`lie_curvature_sphere`), k = 0 is the tangent plane and
/// k = infinity the point. A surface's two curvature spheres are two points
/// on this line.
///
/// Transformations are 6x6 matrices A with A^T eta A = eta. They act on
/// sphere vectors by multiplication and on Dupin cyclides by conjugation
/// (`lie_dupin.hpp`). The Euclidean motions, scalings, Laguerre normal
/// offsets and Möbius inversions are all in the group. A transport `g` between
/// neighbouring shape coordinates is an element of this group. The Lie algebra
/// so(4,2) (X^T eta + eta X = 0) is 15-dimensional.
///
/// Reference: T. Cecil, *Lie Sphere Geometry* (2nd ed., 2008), chapters 1-2.

namespace gaudi {
namespace albers {

typedef Eigen::Matrix<real, 15, 1> vec15;

/// The R^{4,2} metric eta = diag(-1, 1, 1, 1, 1, -1).
inline const mat6 &lie_metric() {
  static const mat6 eta = vec6(-1, 1, 1, 1, 1, -1).asDiagonal();
  return eta;
}

inline real lie_dot(const vec6 &a, const vec6 &b) {
  return -a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3] + a[4] * b[4] -
         a[5] * b[5];
}

/// Oriented sphere with centre `c` and signed radius `r` (see the orientation
/// note above). `r = 0` is the point `c`.
inline vec6 lie_sphere(const vec3 &c, real r) {
  const real c2 = c.squaredNorm();
  vec6 s;
  s << (real(1) + c2 - r * r) / real(2), (real(1) - c2 + r * r) / real(2),
      c[0], c[1], c[2], r;
  return s;
}

inline vec6 lie_point(const vec3 &x) { return lie_sphere(x, real(0)); }

/// Oriented plane `n.x = h` with unit normal `n`.
inline vec6 lie_plane(const vec3 &n, real h) {
  vec6 s;
  s << h, -h, n[0], n[1], n[2], real(1);
  return s;
}

/// The point at infinity, which every plane contains.
inline vec6 lie_infinity() {
  vec6 s;
  s << real(1), real(-1), 0, 0, 0, 0;
  return s;
}

/// Sphere of normal curvature `k` at the surface point `(p, n)`: `k P + T`.
/// Centre `p + n / k`, radius `1 / k`; the tangent plane at `k = 0`.
inline vec6 lie_curvature_sphere(const vec3 &p, const vec3 &n, real k) {
  return k * lie_point(p) + lie_plane(n, n.dot(p));
}

struct lie_sphere_decoded {
  enum class kind { sphere, plane, infinity };
  kind type = kind::infinity;
  /// Sphere: centre and signed radius (a point when the radius is 0).
  vec3 c = vec3::Zero();
  real r = 0;
  /// Plane: unit normal and height, `n.x = h`.
  vec3 n = vec3::Zero();
  real h = 0;
};

/// Reads a sphere, plane or the point at infinity back from a (projective)
/// sphere vector. `tol` is relative to the vector's size.
inline lie_sphere_decoded lie_decode(const vec6 &s, real tol = 1e-12) {
  lie_sphere_decoded out;
  const real scale = std::max(s.norm(), real(1e-300));
  const real w = s[0] + s[1];
  if (std::abs(w) > tol * scale) {
    out.type = lie_sphere_decoded::kind::sphere;
    out.c = s.segment<3>(2) / w;
    out.r = s[5] / w;
  } else if (std::abs(s[5]) > tol * scale) {
    out.type = lie_sphere_decoded::kind::plane;
    out.n = s.segment<3>(2) / s[5];
    out.h = s[0] / s[5];
  }
  return out;
}

/// A contact element: the null 2-plane spanned by `point` and `plane`.
struct lie_contact_element {
  vec6 point;
  vec6 plane;
};

inline lie_contact_element lie_contact(const vec3 &p, const vec3 &n) {
  return {lie_point(p), lie_plane(n, n.dot(p))};
}

/// Recovers the point and normal of the contact element spanned by `a` and
/// `b` (any basis of it, e.g. its two curvature spheres). False when the pair
/// does not span a proper contact element.
inline bool lie_contact_decode(const vec6 &a, const vec6 &b, vec3 &p,
                               vec3 &n, real tol = 1e-12) {
  // The point sphere is the member with radius coordinate 0, the plane the
  // member with x0 + x1 = 0.
  const vec6 point = b[5] * a - a[5] * b;
  const vec6 plane = (b[0] + b[1]) * a - (a[0] + a[1]) * b;
  const lie_sphere_decoded dp = lie_decode(point, tol);
  const lie_sphere_decoded dt = lie_decode(plane, tol);
  if (dp.type != lie_sphere_decoded::kind::sphere ||
      dt.type != lie_sphere_decoded::kind::plane)
    return false;
  p = dp.c;
  n = dt.n;
  return std::abs(n.norm() - real(1)) < real(1e-6);
}

/// Inverse of a group element: `eta A^T eta`.
inline mat6 lie_inverse(const mat6 &A) {
  const mat6 &eta = lie_metric();
  return eta * A.transpose() * eta;
}

/// Largest entry of `A^T eta A - eta`; zero for a group element.
inline real lie_group_defect(const mat6 &A) {
  const mat6 &eta = lie_metric();
  return (A.transpose() * eta * A - eta).cwiseAbs().maxCoeff();
}

/// Largest entry of `X^T eta + eta X`; zero for a Lie algebra element.
inline real lie_algebra_defect(const mat6 &X) {
  const mat6 &eta = lie_metric();
  return (X.transpose() * eta + eta * X).cwiseAbs().maxCoeff();
}

/// Rotation `x -> R x` of R^3.
inline mat6 lie_rotation(const mat3 &R) {
  mat6 A = mat6::Identity();
  A.block<3, 3>(2, 2) = R;
  return A;
}

/// Translation `x -> x + t` of R^3. A sphere keeps its radius.
inline mat6 lie_translation(const vec3 &t) {
  const real t2 = real(0.5) * t.squaredNorm();
  mat6 A = mat6::Identity();
  A(0, 0) += t2;
  A(0, 1) = t2;
  A(1, 0) = -t2;
  A(1, 1) -= t2;
  A.block<1, 3>(0, 2) = t.transpose();
  A.block<1, 3>(1, 2) = -t.transpose();
  A.block<3, 1>(2, 0) = t;
  A.block<3, 1>(2, 1) = t;
  return A;
}

/// Rigid motion `x -> R x + t`.
inline mat6 lie_rigid(const mat3 &R, const vec3 &t) {
  return lie_translation(t) * lie_rotation(R);
}

/// Uniform scaling `x -> lambda x` (`lambda > 0`); radii scale with it.
inline mat6 lie_scale(real lambda) {
  const real a = real(0.5) * (real(1) / lambda + lambda);
  const real b = real(0.5) * (real(1) / lambda - lambda);
  mat6 A = mat6::Identity();
  A(0, 0) = a;
  A(0, 1) = b;
  A(1, 0) = b;
  A(1, 1) = a;
  return A;
}

/// Laguerre normal offset: moves every contact element `d` along its own
/// oriented normal, `(p, n) -> (p + d n, n)`. Sphere radii change by `-d`,
/// planes `n.x = h` move to `n.x = h + d`, a point becomes a sphere.
inline mat6 lie_normal_offset(real d) {
  const real d2 = real(0.5) * d * d;
  mat6 A = mat6::Identity();
  A(0, 0) -= d2;
  A(0, 1) = -d2;
  A(1, 0) = d2;
  A(1, 1) += d2;
  A(0, 5) = d;
  A(1, 5) = -d;
  A(5, 0) = -d;
  A(5, 1) = -d;
  return A;
}

/// Reverses every orientation (radius and normal signs).
inline mat6 lie_orientation_flip() {
  return vec6(1, 1, 1, 1, 1, -1).asDiagonal();
}

/// Reflection in a non-null vector `v`: `x -> x - 2 <x, v> / <v, v> v`.
inline mat6 lie_reflection(const vec6 &v) {
  return mat6::Identity() -
         (real(2) / lie_dot(v, v)) * v * (lie_metric() * v).transpose();
}

/// Unoriented mirror plane `n.x = h` (unit `n`); `lie_reflection` of it is
/// the Euclidean reflection in that plane.
inline vec6 lie_mirror_plane(const vec3 &n, real h) {
  vec6 s;
  s << h, -h, n[0], n[1], n[2], real(0);
  return s;
}

/// Unoriented sphere `(c, rho)` with `<v, v> = 1`; `lie_reflection` of it is
/// the Möbius inversion in that sphere.
inline vec6 lie_inversion_sphere(const vec3 &c, real rho) {
  const real c2 = c.squaredNorm();
  vec6 s;
  s << (real(1) + c2 - rho * rho) / real(2),
      (real(1) - c2 + rho * rho) / real(2), c[0], c[1], c[2], real(0);
  return s / rho;
}

/// Generator of so(4,2) in the plane `(i, j)`, `i < j`:
/// `x -> <x, e_j> e_i - <x, e_i> e_j`. A rotation when both axes have the
/// same sign under eta, a boost otherwise.
inline mat6 lie_algebra_generator(int i, int j) {
  const mat6 &eta = lie_metric();
  mat6 X = mat6::Zero();
  X(i, j) = eta(j, j);
  X(j, i) = -eta(i, i);
  return X;
}

/// The element of so(4,2) with coordinates `w` in the generators
/// `(0,1), (0,2), ..., (4,5)`.
inline mat6 lie_algebra(const vec15 &w) {
  mat6 X = mat6::Zero();
  int k = 0;
  for (int i = 0; i < 6; ++i)
    for (int j = i + 1; j < 6; ++j)
      X += w[k++] * lie_algebra_generator(i, j);
  return X;
}

/// Coordinates of an algebra element in `lie_algebra`'s generators.
inline vec15 lie_algebra_coords(const mat6 &X) {
  const mat6 &eta = lie_metric();
  vec15 w;
  int k = 0;
  for (int i = 0; i < 6; ++i)
    for (int j = i + 1; j < 6; ++j)
      w[k++] = X(i, j) * eta(j, j);
  return w;
}

/// Cayley transform `(I - X/2)^-1 (I + X/2)`: an exact group element for any
/// algebra element `X`, and `exp(X)` to second order. The cheap update step.
inline mat6 lie_cayley(const mat6 &X) {
  const mat6 I = mat6::Identity();
  return (I - real(0.5) * X).partialPivLu().solve(I + real(0.5) * X);
}

inline mat6 lie_exp(const mat6 &X) { return X.exp(); }

/// Principal logarithm; defined when `A` has no eigenvalues on the negative
/// real axis.
inline mat6 lie_log(const mat6 &A) { return A.log(); }

/// `-tr(X^2) / 2`, the invariant (Killing) form, up to scale. Indefinite:
/// positive on rotations, negative on boosts.
inline real lie_algebra_form(const mat6 &X) {
  return real(-0.5) * (X * X).trace();
}

} // namespace albers
} // namespace gaudi

#endif
