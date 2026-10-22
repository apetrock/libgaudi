#ifndef GAUDI_HEPWORTH_SHAPE_OPERATOR_DATUM_HPP
#define GAUDI_HEPWORTH_SHAPE_OPERATOR_DATUM_HPP

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <vector>

#include "gaudi/asawa/datums.hpp"
#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/define_create_func.h"
#include "gaudi/hepworth/block/cr_shape.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

/// Chao, Pinkall, Sanan, Schröder 2010 polar factor. Same determinant fix as
/// `arap::project`: `R = UVᵀ`, and `U.col(2)` flips when `det R < 0`.
inline bool polar_rotation(const mat3 &M, mat3 &R) {
  if (!M.allFinite() || !(M.squaredNorm() > real(1e-16)))
    return false;
  Eigen::JacobiSVD<mat3> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
  mat3 U = svd.matrixU();
  mat3 V = svd.matrixV();
  R = U * V.transpose();
  if (R.determinant() < real(0)) {
    U.col(2) *= real(-1);
    R = U * V.transpose();
  }
  return R.allFinite() && R.determinant() > real(0);
}

/// Turn `T` toward `N × T` by `theta` radians. Zero leaves `T` as it is.
inline vec3 rotate_frame_tangent(const vec3 &T, const vec3 &N, real theta) {
  if (!T.allFinite() || !N.allFinite() || !std::isfinite(theta) ||
      !(std::abs(theta) > real(0)) || !(N.squaredNorm() > real(1e-16)))
    return T;
  const vec3 n = N.normalized();
  vec3 t = T - n * T.dot(n);
  if (!(t.squaredNorm() > real(1e-16)) || !t.allFinite())
    return T;
  t.normalize();
  vec3 b = n.cross(t);
  if (!(b.squaredNorm() > real(1e-16)) || !b.allFinite())
    return T;
  b.normalize();
  vec3 turned = std::cos(theta) * t + std::sin(theta) * b;
  if (!turned.allFinite() || !(turned.squaredNorm() > real(1e-16)))
    return T;
  return turned.normalized();
}

/// Tangent chart of the draw. `T` is `dir` projected into the plane of `N`,
/// and `B` is `N × T`.
inline bool surface_chart(const vec3 &dir, const vec3 &N, vec3 &T, vec3 &B) {
  if (!N.allFinite() || !(N.squaredNorm() > real(1e-16)) || !dir.allFinite())
    return false;
  const vec3 Nn = N.normalized();
  vec3 t = dir - Nn * dir.dot(Nn);
  if (!(t.squaredNorm() > real(1e-16)) || !t.allFinite())
    return false;
  t.normalize();
  vec3 b = Nn.cross(t);
  if (!(b.squaredNorm() > real(1e-16)) || !b.allFinite())
    return false;
  b.normalize();
  t = b.cross(Nn);
  if (!t.allFinite() || !b.allFinite())
    return false;
  T = t;
  B = b;
  return true;
}

/// Orthonormal frame. Column 0 is `dir` projected into the plane of `N`.
inline mat3 tangent_frame(const vec3 &dir, const vec3 &N) {
  const vec3 Nn = N.normalized();
  vec3 t = dir - Nn * dir.dot(Nn);
  t.normalize();
  vec3 b = Nn.cross(t);
  b.normalize();
  t = b.cross(Nn);
  mat3 F;
  F.col(0) = t;
  F.col(1) = b;
  F.col(2) = Nn;
  return F;
}

/// Polar factor of one-ring edge frames. `t` and `-t` build the same frame,
/// so a flipped guidance vector does not cancel the tangent.
inline bool polar_tangent(const std::vector<vec3> &samples, const vec3 &N,
                          vec3 &T) {
  if (!N.allFinite() || !(N.squaredNorm() > real(1e-16)))
    return false;
  const vec3 Nn = N.normalized();
  mat3 M = mat3::Zero();
  vec3 seed = vec3::Zero();
  int used = 0;
  for (vec3 s : samples) {
    if (!s.allFinite())
      continue;
    s -= Nn * s.dot(Nn);
    const real w = s.norm();
    if (!(w > real(1e-8)))
      continue;
    s /= w;
    if (seed.squaredNorm() > real(1e-16) && seed.dot(s) < real(0))
      s = -s;
    seed += w * s;
    M += w * tangent_frame(s, Nn);
    ++used;
  }
  if (used == 0)
    return false;
  mat3 R;
  if (!polar_rotation(M, R))
    return false;
  T = R.col(0);
  T -= Nn * T.dot(Nn);
  if (!(T.squaredNorm() > real(1e-16)) || !T.allFinite())
    return false;
  T.normalize();
  if (seed.dot(T) < real(0))
    T = -T;
  return true;
}

/// Leading singular vector of tangent samples. `t` and `-t` are the same
/// sample, so a sum would cancel them.
inline bool svd_tangent_fit(const std::vector<vec3> &samples, const vec3 &N,
                            vec3 &T) {
  if (!N.allFinite() || !(N.squaredNorm() > real(1e-16)))
    return false;
  const vec3 Nn = N.normalized();
  mat3 C = mat3::Zero();
  int used = 0;
  for (vec3 s : samples) {
    if (!s.allFinite())
      continue;
    s -= Nn * s.dot(Nn);
    if (!(s.squaredNorm() > real(1e-16)))
      continue;
    C += s * s.transpose();
    ++used;
  }
  if (used == 0 || !C.allFinite())
    return false;
  Eigen::JacobiSVD<mat3> svd(C, Eigen::ComputeFullU);
  T = svd.matrixU().col(0);
  T -= Nn * T.dot(Nn);
  if (!(T.squaredNorm() > real(1e-16)) || !T.allFinite())
    return false;
  T.normalize();
  return true;
}

inline void shape_sym2(real a, real b, real c, real &k1, real &k2, real &x1,
                       real &y1) {
  const real tr = a + c;
  const real det = a * c - b * b;
  const real disc = std::sqrt(std::max(real(0), tr * tr * real(0.25) - det));
  k1 = real(0.5) * tr + disc;
  k2 = real(0.5) * tr - disc;
  real x = b;
  real y = k1 - a;
  const real n2 = x * x + y * y;
  if (n2 < real(1e-24)) {
    x1 = real(1);
    y1 = real(0);
  } else {
    const real inv = real(1) / std::sqrt(n2);
    x1 = x * inv;
    y1 = y * inv;
  }
}

/// World-space shape operator from chart coeffs. Same embedding that
/// `transport_shape` reads when it maps one chart into another.
inline mat3 shape_world_from_coeffs(const vec3 &h, const vec3 &T, const vec3 &B) {
  return h[0] * (T * T.transpose()) +
         h[1] * (T * B.transpose() + B * T.transpose()) +
         h[2] * (B * B.transpose());
}

/// Chart coeffs of a world operator. Matches `transport_shape` output when
/// `S` is the world image of a neighbor chart.
inline vec3 shape_coeffs_from_world(const mat3 &S, const vec3 &T, const vec3 &B) {
  const real k11 = T.dot(S * T);
  const real k22 = B.dot(S * B);
  const real k12 = real(0.5) * (T.dot(S * B) + B.dot(S * T));
  return vec3(k11, k12, k22);
}

/// Interior-edge diamond from a mesh corner. Same verts the projection pin uses.
inline bool diamond_from_corner(const asawa::shell::shell &M,
                                asawa::shell::CornerId c0, int &v0, int &v1,
                                int &o0, int &o1, std::array<int, 3> &t0,
                                std::array<int, 3> &t1) {
  const asawa::shell::CornerId c1 = M.other(c0);
  v0 = static_cast<int>(M.vert(c0));
  v1 = static_cast<int>(M.vert(c1));
  const asawa::shell::FaceId f0 = M.face(c0);
  const asawa::shell::FaceId f1 = M.face(c1);
  if (static_cast<int>(f0) < 0 || static_cast<int>(f1) < 0)
    return false;
  if (M.fsize(f0) != 3 || M.fsize(f1) != 3)
    return false;
  const std::array<asawa::shell::VertId, 3> a = M.get_tri(f0);
  const std::array<asawa::shell::VertId, 3> b = M.get_tri(f1);
  o0 = -1;
  o1 = -1;
  for (int i = 0; i < 3; ++i) {
    t0[static_cast<size_t>(i)] = static_cast<int>(a[static_cast<size_t>(i)]);
    t1[static_cast<size_t>(i)] = static_cast<int>(b[static_cast<size_t>(i)]);
    if (t0[static_cast<size_t>(i)] != v0 && t0[static_cast<size_t>(i)] != v1)
      o0 = t0[static_cast<size_t>(i)];
    if (t1[static_cast<size_t>(i)] != v0 && t1[static_cast<size_t>(i)] != v1)
      o1 = t1[static_cast<size_t>(i)];
  }
  return v0 >= 0 && v1 >= 0 && o0 >= 0 && o1 >= 0;
}

/// Static diamond chart: edge-parallel `diamond_frame` (gp=1, gq=0).
/// Polar of `Σ spoke spokeᵀ` is wrong here — that covariance is symmetric, so
/// its polar is near identity and the chart snaps to a world axis.
inline bool diamond_polar_chart(const vec3 &p0, const vec3 &p1, const vec3 &a,
                                const vec3 &b, const std::array<int, 3> &t0,
                                const std::array<int, 3> &t1, int v0, int v1,
                                int o0, int o1, vec3 &N, vec3 &T, vec3 &B) {
  return diamond_frame(p0, p1, a, b, t0, t1, v0, v1, o0, o1, real(1), real(0),
                       N, T, B);
}

/// Blend endpoint operators on an interior diamond. Endpoint charts use stored
/// polar tangents when present. Destination `T` is the sign-aligned average of
/// those hints, projected into the diamond tangent plane — so aligned neighbors
/// stay aligned, including when their shared angle is not the edge. World
/// blend carries chart twists; fallback chart is edge-parallel `diamond_frame`.
inline bool blend_shape_on_diamond(const vec3 &p0, const vec3 &p1, const vec3 &a,
                                   const vec3 &b, const std::array<int, 3> &t0,
                                   const std::array<int, 3> &t1, int v0, int v1,
                                   int o0, int o1, real s, const vec3 &N0,
                                   const vec3 &N1, const vec3 &h0, const vec3 &h1,
                                   const vec3 *T0_hint, const vec3 *T1_hint,
                                   vec3 &h_out, vec3 &T_out) {
  h_out = vec3::Zero();
  T_out = vec3::Zero();
  if (!p0.allFinite() || !p1.allFinite() || !a.allFinite() || !b.allFinite() ||
      !N0.allFinite() || !N1.allFinite() || !h0.allFinite() || !h1.allFinite() ||
      !std::isfinite(s))
    return false;
  vec3 Nd;
  if (!diamond_normal(p0, p1, a, b, t0, t1, v0, v1, o0, o1, Nd))
    return false;
  const vec3 edge = p1 - p0;
  vec3 T0, B0, T1, B1;
  const vec3 hint0 =
      (T0_hint && T0_hint->allFinite() && T0_hint->squaredNorm() > real(1e-16))
          ? *T0_hint
          : edge;
  const vec3 hint1 =
      (T1_hint && T1_hint->allFinite() && T1_hint->squaredNorm() > real(1e-16))
          ? *T1_hint
          : edge;
  if (!surface_chart(hint0, N0, T0, B0) || !surface_chart(hint1, N1, T1, B1))
    return false;
  if (T1.dot(T0) < real(0)) {
    T1 = -T1;
    B1 = -B1;
  }
  vec3 Td, Bd;
  const vec3 Tm = (real(1) - s) * T0 + s * T1;
  if (!surface_chart(Tm, Nd, Td, Bd) &&
      !diamond_frame(p0, p1, a, b, t0, t1, v0, v1, o0, o1, real(1), real(0), Nd,
                     Td, Bd))
    return false;
  const mat3 S0 = shape_world_from_coeffs(h0, T0, B0);
  const mat3 S1 = shape_world_from_coeffs(h1, T1, B1);
  if (!S0.allFinite() || !S1.allFinite())
    return false;
  const mat3 S = (real(1) - s) * S0 + s * S1;
  h_out = shape_coeffs_from_world(S, Td, Bd);
  T_out = Td;
  return h_out.allFinite() && T_out.allFinite();
}

/// Fallback when the edge is not an interior diamond (boundary).
inline bool blend_shape_on_edge(const vec3 &x0, const vec3 &x1, real s,
                                const vec3 &N0, const vec3 &N1, const vec3 &h0,
                                const vec3 &h1, const vec3 *T0_hint,
                                const vec3 *T1_hint, vec3 &h_out, vec3 &T_out) {
  h_out = vec3::Zero();
  T_out = vec3::Zero();
  if (!x0.allFinite() || !x1.allFinite() || !N0.allFinite() || !N1.allFinite() ||
      !h0.allFinite() || !h1.allFinite() || !std::isfinite(s))
    return false;
  const vec3 edge = x1 - x0;
  if (!(edge.squaredNorm() > real(1e-16)))
    return false;
  vec3 T0, B0, T1, B1;
  const vec3 hint0 =
      (T0_hint && T0_hint->allFinite() && T0_hint->squaredNorm() > real(1e-16))
          ? *T0_hint
          : edge;
  const vec3 hint1 =
      (T1_hint && T1_hint->allFinite() && T1_hint->squaredNorm() > real(1e-16))
          ? *T1_hint
          : edge;
  if (!surface_chart(hint0, N0, T0, B0) || !surface_chart(hint1, N1, T1, B1))
    return false;
  if (T1.dot(T0) < real(0)) {
    T1 = -T1;
    B1 = -B1;
  }
  vec3 Nn = (real(1) - s) * N0 + s * N1;
  if (!(Nn.squaredNorm() > real(1e-16)) || !Nn.allFinite())
    return false;
  Nn.normalize();
  vec3 Tm = (real(1) - s) * T0 + s * T1;
  vec3 Tn, Bn;
  if (!surface_chart(Tm, Nn, Tn, Bn) && !surface_chart(edge, Nn, Tn, Bn))
    return false;
  const mat3 S0 = shape_world_from_coeffs(h0, T0, B0);
  const mat3 S1 = shape_world_from_coeffs(h1, T1, B1);
  if (!S0.allFinite() || !S1.allFinite())
    return false;
  const mat3 S = (real(1) - s) * S0 + s * S1;
  h_out = shape_coeffs_from_world(S, Tn, Bn);
  T_out = Tn;
  return h_out.allFinite() && T_out.allFinite();
}

/// Vertex Weingarten sample. `data()` is `(k11, k12, k22)` in `chart_T`.
/// Remesh averages world operators (same map as `transport_shape`), not raw
/// coeffs, and writes the destination chart into `chart_T`.
struct shape_operator_datum : public asawa::datum_t<vec3> {
  DEFINE_CREATE_FUNC(shape_operator_datum)

  shape_operator_datum(asawa::prim_type type, const std::vector<vec3> &data)
      : asawa::datum_t<vec3>(type, data) {
    fit_sides(static_cast<int>(data.size()));
  }

  std::vector<vec3> rest;
  /// Chart tangent for `data()` / `rest`. Remesh and seed keep this aligned
  /// with the projection polar frames.
  std::vector<vec3> chart_T;
  std::vector<real> mark;

  void fit_sides(int n) {
    auto fit3 = [n](std::vector<vec3> &v) {
      if (static_cast<int>(v.size()) < n)
        v.resize(static_cast<size_t>(n), vec3::Zero());
      else if (static_cast<int>(v.size()) > n)
        v.resize(static_cast<size_t>(n));
    };
    fit3(rest);
    fit3(chart_T);
    if (static_cast<int>(mark.size()) < n)
      mark.resize(static_cast<size_t>(n), real(0));
    else if (static_cast<int>(mark.size()) > n)
      mark.resize(static_cast<size_t>(n));
  }

  void subdivide(const asawa::shell::shell &M, const index_t &i0,
                 const index_t &i1, const index_t &i2, const index_t &i3,
                 const real &s, const index_t &source_corner) override {
    (void)i1;
    (void)i2;
    (void)i3;
    const asawa::shell::CornerId c0 = asawa::shell::corner_id(source_corner);
    int v0 = -1, v1 = -1, o0 = -1, o1 = -1;
    std::array<int, 3> t0{{0, 0, 0}};
    std::array<int, 3> t1{{0, 0, 0}};
    const bool have_diamond = diamond_from_corner(M, c0, v0, v1, o0, o1, t0, t1);
    if (!have_diamond) {
      v0 = static_cast<int>(M.vert(c0));
      v1 = static_cast<int>(M.vert(M.other(c0)));
    }
    const int n = static_cast<int>(data().size());
    fit_sides(n);
    if (i0 < 0 || i0 >= n || v0 < 0 || v1 < 0 || v0 >= n || v1 >= n)
      return;
    const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);
    const int xmax = have_diamond ? std::max(o0, std::max(o1, std::max(v0, v1)))
                                  : std::max(v0, v1);
    if (static_cast<int>(x.size()) <= xmax)
      return;
    const vec3 &x0 = x[static_cast<size_t>(v0)];
    const vec3 &x1 = x[static_cast<size_t>(v1)];
    const vec3 N0 = asawa::shell::vert_normal(M, asawa::shell::vert_id(v0), x);
    const vec3 N1 = asawa::shell::vert_normal(M, asawa::shell::vert_id(v1), x);
    const vec3 *T0 = chart_T[static_cast<size_t>(v0)].squaredNorm() > real(1e-16)
                         ? &chart_T[static_cast<size_t>(v0)]
                         : nullptr;
    const vec3 *T1 = chart_T[static_cast<size_t>(v1)].squaredNorm() > real(1e-16)
                         ? &chart_T[static_cast<size_t>(v1)]
                         : nullptr;
    auto blend_h = [&](const vec3 &ha, const vec3 &hb, vec3 &hout,
                       vec3 &Tout) -> bool {
      if (have_diamond) {
        const vec3 &oa = x[static_cast<size_t>(o0)];
        const vec3 &ob = x[static_cast<size_t>(o1)];
        if (blend_shape_on_diamond(x0, x1, oa, ob, t0, t1, v0, v1, o0, o1, s, N0,
                                   N1, ha, hb, T0, T1, hout, Tout))
          return true;
      }
      return blend_shape_on_edge(x0, x1, s, N0, N1, ha, hb, T0, T1, hout, Tout);
    };
    vec3 h_new, T_new;
    if (blend_h(data()[static_cast<size_t>(v0)], data()[static_cast<size_t>(v1)],
                h_new, T_new)) {
      data()[static_cast<size_t>(i0)] = h_new;
      chart_T[static_cast<size_t>(i0)] = T_new;
    } else {
      data()[static_cast<size_t>(i0)] =
          va::mix(s, data()[static_cast<size_t>(v0)],
                  data()[static_cast<size_t>(v1)]);
      chart_T[static_cast<size_t>(i0)] =
          va::mix(s, chart_T[static_cast<size_t>(v0)],
                  chart_T[static_cast<size_t>(v1)]);
    }
    vec3 r_new, Tr;
    if (blend_h(rest[static_cast<size_t>(v0)], rest[static_cast<size_t>(v1)],
                r_new, Tr))
      rest[static_cast<size_t>(i0)] = r_new;
    else
      rest[static_cast<size_t>(i0)] = va::mix(
          s, rest[static_cast<size_t>(v0)], rest[static_cast<size_t>(v1)]);
    mark[static_cast<size_t>(i0)] =
        va::mix(s, mark[static_cast<size_t>(v0)], mark[static_cast<size_t>(v1)]);
  }

  void collapse(const asawa::shell::shell &M, const index_t &i) override {
    const asawa::shell::CornerId c0 = asawa::shell::corner_id(i);
    int v0 = -1, v1 = -1, o0 = -1, o1 = -1;
    std::array<int, 3> t0{{0, 0, 0}};
    std::array<int, 3> t1{{0, 0, 0}};
    const bool have_diamond = diamond_from_corner(M, c0, v0, v1, o0, o1, t0, t1);
    if (!have_diamond) {
      v0 = static_cast<int>(M.vert(c0));
      v1 = static_cast<int>(M.vert(M.other(c0)));
    }
    const int n = static_cast<int>(data().size());
    fit_sides(n);
    if (v0 < 0 || v1 < 0 || v0 >= n || v1 >= n)
      return;
    const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);
    const int xmax = have_diamond ? std::max(o0, std::max(o1, std::max(v0, v1)))
                                  : std::max(v0, v1);
    if (static_cast<int>(x.size()) <= xmax) {
      asawa::datum_t<vec3>::collapse(M, i);
      return;
    }
    const vec3 &x0 = x[static_cast<size_t>(v0)];
    const vec3 &x1 = x[static_cast<size_t>(v1)];
    const vec3 N0 = asawa::shell::vert_normal(M, asawa::shell::vert_id(v0), x);
    const vec3 N1 = asawa::shell::vert_normal(M, asawa::shell::vert_id(v1), x);
    const vec3 *T0 = chart_T[static_cast<size_t>(v0)].squaredNorm() > real(1e-16)
                         ? &chart_T[static_cast<size_t>(v0)]
                         : nullptr;
    const vec3 *T1 = chart_T[static_cast<size_t>(v1)].squaredNorm() > real(1e-16)
                         ? &chart_T[static_cast<size_t>(v1)]
                         : nullptr;
    auto blend_h = [&](const vec3 &ha, const vec3 &hb, vec3 &hout,
                       vec3 &Tout) -> bool {
      if (have_diamond) {
        const vec3 &oa = x[static_cast<size_t>(o0)];
        const vec3 &ob = x[static_cast<size_t>(o1)];
        if (blend_shape_on_diamond(x0, x1, oa, ob, t0, t1, v0, v1, o0, o1,
                                   real(0.5), N0, N1, ha, hb, T0, T1, hout,
                                   Tout))
          return true;
      }
      return blend_shape_on_edge(x0, x1, real(0.5), N0, N1, ha, hb, T0, T1, hout,
                                 Tout);
    };
    vec3 h_new, T_new;
    if (blend_h(data()[static_cast<size_t>(v0)], data()[static_cast<size_t>(v1)],
                h_new, T_new)) {
      data()[static_cast<size_t>(v0)] = h_new;
      chart_T[static_cast<size_t>(v0)] = T_new;
    } else {
      asawa::datum_t<vec3>::collapse(M, i);
    }
    vec3 r_new, Tr;
    if (blend_h(rest[static_cast<size_t>(v0)], rest[static_cast<size_t>(v1)],
                r_new, Tr))
      rest[static_cast<size_t>(v0)] = r_new;
    else
      rest[static_cast<size_t>(v0)] =
          va::mix(real(0.5), rest[static_cast<size_t>(v0)],
                  rest[static_cast<size_t>(v1)]);
    mark[static_cast<size_t>(v0)] =
        va::mix(real(0.5), mark[static_cast<size_t>(v0)],
                mark[static_cast<size_t>(v1)]);
  }

  void merge(const asawa::shell::shell &M, const index_t &i0, const index_t &i1,
             const index_t &i2, const index_t &i3, const index_t vA0,
             const index_t vA1, const index_t vB0, const index_t vB1) override {
    const int n = static_cast<int>(data().size());
    fit_sides(n);
    const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);
    auto blend_pair = [&](index_t dst, index_t a, index_t b, real s) {
      if (dst < 0 || dst >= n || a < 0 || b < 0 || a >= n || b >= n)
        return;
      if (static_cast<int>(x.size()) <= std::max(a, b)) {
        data()[static_cast<size_t>(dst)] = va::mix(
            s, data()[static_cast<size_t>(a)], data()[static_cast<size_t>(b)]);
        rest[static_cast<size_t>(dst)] = va::mix(
            s, rest[static_cast<size_t>(a)], rest[static_cast<size_t>(b)]);
        chart_T[static_cast<size_t>(dst)] = va::mix(
            s, chart_T[static_cast<size_t>(a)], chart_T[static_cast<size_t>(b)]);
        mark[static_cast<size_t>(dst)] = va::mix(
            s, mark[static_cast<size_t>(a)], mark[static_cast<size_t>(b)]);
        return;
      }
      const vec3 &xa = x[static_cast<size_t>(a)];
      const vec3 &xb = x[static_cast<size_t>(b)];
      const vec3 Na = asawa::shell::vert_normal(M, asawa::shell::vert_id(a), x);
      const vec3 Nb = asawa::shell::vert_normal(M, asawa::shell::vert_id(b), x);
      const vec3 *Ta = chart_T[static_cast<size_t>(a)].squaredNorm() > real(1e-16)
                           ? &chart_T[static_cast<size_t>(a)]
                           : nullptr;
      const vec3 *Tb = chart_T[static_cast<size_t>(b)].squaredNorm() > real(1e-16)
                           ? &chart_T[static_cast<size_t>(b)]
                           : nullptr;
      vec3 h_new, T_new;
      if (blend_shape_on_edge(xa, xb, s, Na, Nb, data()[static_cast<size_t>(a)],
                              data()[static_cast<size_t>(b)], Ta, Tb, h_new,
                              T_new)) {
        data()[static_cast<size_t>(dst)] = h_new;
        chart_T[static_cast<size_t>(dst)] = T_new;
      } else {
        data()[static_cast<size_t>(dst)] = va::mix(
            s, data()[static_cast<size_t>(a)], data()[static_cast<size_t>(b)]);
        chart_T[static_cast<size_t>(dst)] = va::mix(
            s, chart_T[static_cast<size_t>(a)], chart_T[static_cast<size_t>(b)]);
      }
      vec3 r_new, Tr;
      if (blend_shape_on_edge(xa, xb, s, Na, Nb, rest[static_cast<size_t>(a)],
                              rest[static_cast<size_t>(b)], Ta, Tb, r_new, Tr))
        rest[static_cast<size_t>(dst)] = r_new;
      else
        rest[static_cast<size_t>(dst)] = va::mix(
            s, rest[static_cast<size_t>(a)], rest[static_cast<size_t>(b)]);
      mark[static_cast<size_t>(dst)] = va::mix(
          s, mark[static_cast<size_t>(a)], mark[static_cast<size_t>(b)]);
    };
    blend_pair(i0, vA0, vB0, real(0.5));
    blend_pair(i1, vA1, vB1, real(0.5));
    blend_pair(i2, vA0, vB0, real(1));
    blend_pair(i3, vA1, vB1, real(1));
  }

  void do_alloc(const size_t &sz) override {
    asawa::datum_t<vec3>::do_alloc(sz);
    fit_sides(static_cast<int>(data().size()));
  }

  void do_resize(const size_t &sz) override {
    asawa::datum_t<vec3>::do_resize(sz);
    fit_sides(static_cast<int>(sz));
  }

  void do_permute(const std::vector<index_t> &permute) override {
    asawa::datum_t<vec3>::do_permute(permute);
    permute_vec(rest, permute);
    permute_vec(chart_T, permute);
    if (mark.empty())
      return;
    std::vector<real> next = mark;
    const size_t n = std::min(mark.size(), permute.size());
    for (size_t i = 0; i < n; ++i) {
      if (permute[i] >= 0 &&
          static_cast<size_t>(permute[i]) < mark.size())
        next[i] = mark[static_cast<size_t>(permute[i])];
    }
    mark = std::move(next);
  }

private:
  static void permute_vec(std::vector<vec3> &v,
                          const std::vector<index_t> &permute) {
    if (v.empty())
      return;
    std::vector<vec3> next = v;
    const size_t n = std::min(v.size(), permute.size());
    for (size_t i = 0; i < n; ++i) {
      if (permute[i] >= 0 && static_cast<size_t>(permute[i]) < v.size())
        next[i] = v[static_cast<size_t>(permute[i])];
    }
    v = std::move(next);
  }
};

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif
