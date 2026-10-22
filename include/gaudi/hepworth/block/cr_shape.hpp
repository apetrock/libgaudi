#ifndef GAUDI_HEPWORTH_CR_SHAPE_HPP
#define GAUDI_HEPWORTH_CR_SHAPE_HPP

/// Shape operator from the Crouzeix–Raviart derivative of the embedding.
/// Stein, Wardetzky, Jacobson, Grinspun 2020, eq. (9) and Appendix A.2
/// (`reference/stein2020-vector-dirichlet-energy.pdf`).
/// The vector Dirichlet energy is the tangential half of that derivative.
/// This is the normal half. Tables are vertex and edge indices. No shell.

#include "gaudi/common.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace gaudi {
namespace hepworth {
namespace block {

/// One triangle. Vertex indices into the position array, edge-dof indices.
struct cr_tri {
  std::array<int, 3> v{-1, -1, -1};
  std::array<int, 3> e{-1, -1, -1};
  /// Face indices of the 1-ring, including this face. Samples for the curvature jet.
  std::vector<int> stencil;
};

/// Incident face of one vertex: three vertex ids and three edge-dof ids.
struct fan_face {
  std::array<int, 3> v{-1, -1, -1};
  std::array<int, 3> e{-1, -1, -1};
};

/// Edge dof plus its canonical endpoints. `dof` indexes the coefficient block.
struct fan_edge {
  int dof = -1;
  int v0 = -1;
  int v1 = -1;
};

/// One Stein edge. Endpoints are vert(c) and vert(other(c)). The parallel
/// axis is sign(v0 <= v1) * (x[v1] - x[v0]), lesser vertex id toward greater.
struct cr_edge {
  int v0 = -1;
  int v1 = -1;
  std::vector<int> faces;
  /// Neighbor edge dof and the face index they share.
  std::vector<std::pair<int, int>> nbr;
};

/// Per-face operator. `S` maps the tangent plane into itself. `area` is |f|.
struct cr_face_operator {
  mat3 S = mat3::Zero();
  real area = 0;
};

/// One vertex, after the polar of its incident edge tangents. `h` is
/// `(k11, k12, k22)` of the area-averaged operator in `(T, B)`.
struct vertex_chart {
  vec3 T = vec3::Zero();
  vec3 B = vec3::Zero();
  vec3 N = vec3::Zero();
  vec3 h = vec3::Zero();
  bool ok = false;
};

/// Appendix A.2 dotted with N. The fairing jet stores −dN, so the paper's
/// operator is multiplied by −1 once here.
constexpr real k_cr_shape_sign = real(-1);

inline vec3 cr_edge_parallel(const std::vector<vec3> &x, const cr_edge &e) {
  if (e.v0 < 0 || e.v1 < 0 || e.v0 >= static_cast<int>(x.size()) ||
      e.v1 >= static_cast<int>(x.size()))
    return vec3::UnitX();
  vec3 t = x[static_cast<size_t>(e.v1)] - x[static_cast<size_t>(e.v0)];
  const real len = t.norm();
  if (!(len > real(1e-20)))
    return vec3::UnitX();
  const real s = e.v0 <= e.v1 ? real(1) : real(-1);
  return s * t / len;
}

inline vec3 cr_edge_perp(const vec3 &e_par, const vec3 &n_unit) {
  vec3 e = n_unit.cross(e_par);
  const real lp = e.norm();
  if (!(lp > real(1e-20)))
    return vec3::UnitY();
  return e / lp;
}

inline vec3 cr_face_cross(const std::vector<vec3> &x, const cr_tri &tri) {
  if (tri.v[0] < 0 || tri.v[1] < 0 || tri.v[2] < 0)
    return vec3::Zero();
  if (tri.v[0] >= static_cast<int>(x.size()) ||
      tri.v[1] >= static_cast<int>(x.size()) ||
      tri.v[2] >= static_cast<int>(x.size()))
    return vec3::Zero();
  const vec3 &x0 = x[static_cast<size_t>(tri.v[0])];
  const vec3 &x1 = x[static_cast<size_t>(tri.v[1])];
  const vec3 &x2 = x[static_cast<size_t>(tri.v[2])];
  return (x1 - x0).cross(x2 - x0);
}

/// Stein eq. (7) mass is A_e / 3, A_e the sum of the incident face areas.
/// G = M^{-1} D, then each face accumulates the two outer products in A.2.
///
/// Parallel and perpendicular coordinate gradients, one pair per edge.
struct cr_edge_gradient {
  std::vector<vec3> par;
  std::vector<vec3> perp;
};

/// `G = M^{-1} D(u)` for the three coordinate functions. An interior edge
/// accumulates both incident faces.
inline cr_edge_gradient
cr_coordinate_gradients(const std::vector<vec3> &x, const std::vector<cr_tri> &tris,
                        const std::vector<cr_edge> &edges) {
  cr_edge_gradient out;
  const int nF = static_cast<int>(tris.size());
  const int nE = static_cast<int>(edges.size());
  out.par.assign(static_cast<size_t>(std::max(nE, 0)), vec3::Zero());
  out.perp.assign(static_cast<size_t>(std::max(nE, 0)), vec3::Zero());
  if (nF <= 0 || nE <= 0)
    return out;

  struct acc {
    vec3 par = vec3::Zero();
    vec3 perp = vec3::Zero();
    real area = 0;
  };
  std::vector<acc> D(static_cast<size_t>(nE));

  auto edge_len = [&](const cr_tri &tri, int i) -> real {
    const int a = tri.v[static_cast<size_t>(i)];
    const int b = tri.v[static_cast<size_t>((i + 1) % 3)];
    if (a < 0 || b < 0 || a >= static_cast<int>(x.size()) ||
        b >= static_cast<int>(x.size()))
      return 0;
    return (x[static_cast<size_t>(b)] - x[static_cast<size_t>(a)]).norm();
  };

  for (int f = 0; f < nF; ++f) {
    const cr_tri &tri = tris[static_cast<size_t>(f)];
    const vec3 cr = cr_face_cross(x, tri);
    const real area = real(0.5) * cr.norm();
    if (!(area > real(1e-20)) || !cr.allFinite())
      continue;
    real len[3];
    for (int i = 0; i < 3; ++i)
      len[i] = edge_len(tri, i);
    if (!(len[0] > real(1e-12)) || !(len[1] > real(1e-12)) ||
        !(len[2] > real(1e-12)))
      continue;
    for (int i = 0; i < 3; ++i) {
      const int e = tri.e[static_cast<size_t>(i)];
      if (e < 0 || e >= nE || !(len[i] > real(1e-20)))
        continue;
      const int tail = tri.v[static_cast<size_t>(i)];
      const int head = tri.v[static_cast<size_t>((i + 1) % 3)];
      const int opp = tri.v[static_cast<size_t>((i + 2) % 3)];
      if (tail < 0 || head < 0 || opp < 0 ||
          tail >= static_cast<int>(x.size()) ||
          head >= static_cast<int>(x.size()) ||
          opp >= static_cast<int>(x.size()))
        continue;
      const vec3 e_par = cr_edge_parallel(x, edges[static_cast<size_t>(e)]);
      const vec3 local = x[static_cast<size_t>(head)] - x[static_cast<size_t>(tail)];
      const real si = local.dot(e_par) >= real(0) ? real(1) : real(-1);
      const vec3 ui = x[static_cast<size_t>(tail)];
      const vec3 uj = x[static_cast<size_t>(head)];
      const vec3 uk = x[static_cast<size_t>(opp)];
      const real lenj = len[(i + 1) % 3];
      const real lenk = len[(i + 2) % 3];
      const real li = len[i];
      // cr_embedding_derivative.cpp: ejk and eki are squared lengths.
      const real ejk = lenj * lenj;
      const real eki = lenk * lenk;
      acc &d = D[static_cast<size_t>(e)];
      d.area += area;
      d.par += si * area / (real(3) * li) * (uj - ui);
      d.perp += -si / (real(12) * li) *
                    (li * li * (uj + ui) + (eki - ejk) * (uj - ui)) +
                si * li / real(6) * uk;
    }
  }

  for (int e = 0; e < nE; ++e) {
    const acc &d = D[static_cast<size_t>(e)];
    const real M = d.area / real(3);
    if (!(M > real(1e-20)))
      continue;
    out.par[static_cast<size_t>(e)] = d.par / M;
    out.perp[static_cast<size_t>(e)] = d.perp / M;
  }
  return out;
}

inline std::vector<cr_face_operator>
cr_face_operators(const std::vector<vec3> &x, const std::vector<cr_tri> &tris,
                  const std::vector<cr_edge> &edges,
                  const std::vector<vec3> *Gpar = nullptr,
                  const std::vector<vec3> *Gperp = nullptr) {
  const int nF = static_cast<int>(tris.size());
  const int nE = static_cast<int>(edges.size());
  std::vector<cr_face_operator> out(static_cast<size_t>(std::max(nF, 0)));
  if (nF <= 0 || nE <= 0)
    return out;

  for (int f = 0; f < nF; ++f) {
    const cr_tri &tri = tris[static_cast<size_t>(f)];
    const vec3 cr = cr_face_cross(x, tri);
    const real area = real(0.5) * cr.norm();
    if (!(area > real(1e-20)) || !cr.allFinite())
      continue;
    bool open = false;
    for (int i = 0; i < 3; ++i) {
      const int a = tri.v[static_cast<size_t>(i)];
      const int b = tri.v[static_cast<size_t>((i + 1) % 3)];
      if (a < 0 || b < 0 || a >= static_cast<int>(x.size()) ||
          b >= static_cast<int>(x.size()) ||
          !((x[static_cast<size_t>(b)] - x[static_cast<size_t>(a)]).norm() >
            real(1e-12)))
        open = true;
    }
    if (open)
      continue;
    out[static_cast<size_t>(f)].area = area;
  }

  cr_edge_gradient G;
  const bool cached = Gpar && Gperp &&
                      static_cast<int>(Gpar->size()) >= nE &&
                      static_cast<int>(Gperp->size()) >= nE;
  if (!cached)
    G = cr_coordinate_gradients(x, tris, edges);
  const std::vector<vec3> &par = cached ? *Gpar : G.par;
  const std::vector<vec3> &perp = cached ? *Gperp : G.perp;

  for (int f = 0; f < nF; ++f) {
    cr_face_operator &face = out[static_cast<size_t>(f)];
    if (!(face.area > real(1e-20)))
      continue;
    const cr_tri &tri = tris[static_cast<size_t>(f)];
    const vec3 cr = cr_face_cross(x, tri);
    if (!(cr.squaredNorm() > real(1e-24)))
      continue;
    const vec3 N = cr.normalized();
    mat3 S = mat3::Zero();
    for (int i = 0; i < 3; ++i) {
      const int e = tri.e[static_cast<size_t>(i)];
      if (e < 0 || e >= nE)
        continue;
      const int tail = tri.v[static_cast<size_t>(i)];
      const int head = tri.v[static_cast<size_t>((i + 1) % 3)];
      if (tail < 0 || head < 0 || tail >= static_cast<int>(x.size()) ||
          head >= static_cast<int>(x.size()))
        continue;
      vec3 vpar = x[static_cast<size_t>(head)] - x[static_cast<size_t>(tail)];
      const real li = vpar.norm();
      if (!(li > real(1e-20)))
        continue;
      vpar /= li;
      const vec3 e_par = cr_edge_parallel(x, edges[static_cast<size_t>(e)]);
      const real si = vpar.dot(e_par) >= real(0) ? real(1) : real(-1);
      vec3 vperp = N.cross(vpar);
      if (!(vperp.squaredNorm() > real(1e-24)))
        continue;
      vperp.normalize();
      const real scale = si * li / face.area;
      const real c_par = scale * N.dot(par[static_cast<size_t>(e)]);
      const real c_perp = scale * N.dot(perp[static_cast<size_t>(e)]);
      if (!std::isfinite(c_par) || !std::isfinite(c_perp))
        continue;
      S += c_par * vperp * vpar.transpose();
      S += c_perp * vperp * vperp.transpose();
    }
    if (S.allFinite())
      face.S = k_cr_shape_sign * S;
  }
  return out;
}

/// One face's `(k11, k12, k22)` from cached edge gradients. `T` is the first
/// edge, in the face plane. `B` is `N × T`. `v[0], v[1], v[2]` is the winding.
inline bool face_coeff_from_gradients(const std::vector<vec3> &x,
                                      const std::array<int, 3> &v,
                                      const std::array<int, 3> &edge,
                                      const std::vector<cr_edge> &edges,
                                      const std::vector<vec3> &Gpar,
                                      const std::vector<vec3> &Gperp, vec3 &coeff,
                                      vec3 &N, vec3 &T, vec3 &B) {
  for (int i = 0; i < 3; ++i) {
    if (v[static_cast<size_t>(i)] < 0 ||
        v[static_cast<size_t>(i)] >= static_cast<int>(x.size()))
      return false;
  }
  const vec3 &x0 = x[static_cast<size_t>(v[0])];
  const vec3 &x1 = x[static_cast<size_t>(v[1])];
  const vec3 &x2 = x[static_cast<size_t>(v[2])];
  const vec3 crv = (x1 - x0).cross(x2 - x0);
  const real area = real(0.5) * crv.norm();
  if (!(area > real(1e-20)) || !crv.allFinite())
    return false;
  N = crv.normalized();
  vec3 t = x1 - x0;
  t -= N * t.dot(N);
  if (!(t.squaredNorm() > real(1e-16)) || !t.allFinite())
    return false;
  t.normalize();
  vec3 b = N.cross(t);
  if (!(b.squaredNorm() > real(1e-16)) || !b.allFinite())
    return false;
  b.normalize();
  T = t;
  B = b;
  const vec3 xv[3] = {x0, x1, x2};
  mat3 S = mat3::Zero();
  const int nE = static_cast<int>(edges.size());
  for (int i = 0; i < 3; ++i) {
    const int e = edge[static_cast<size_t>(i)];
    if (e < 0 || e >= nE || e >= static_cast<int>(Gpar.size()) ||
        e >= static_cast<int>(Gperp.size()))
      continue;
    vec3 vpar = xv[(i + 1) % 3] - xv[i];
    const real li = vpar.norm();
    if (!(li > real(1e-20)))
      continue;
    vpar /= li;
    const vec3 e_par = cr_edge_parallel(x, edges[static_cast<size_t>(e)]);
    const real si = vpar.dot(e_par) >= real(0) ? real(1) : real(-1);
    vec3 vperp = N.cross(vpar);
    if (!(vperp.squaredNorm() > real(1e-24)))
      continue;
    vperp.normalize();
    const real scale = si * li / area;
    const real c_par = scale * N.dot(Gpar[static_cast<size_t>(e)]);
    const real c_perp = scale * N.dot(Gperp[static_cast<size_t>(e)]);
    if (!std::isfinite(c_par) || !std::isfinite(c_perp))
      continue;
    S += c_par * vperp * vpar.transpose();
    S += c_perp * vperp * vperp.transpose();
  }
  S = k_cr_shape_sign * S;
  if (!S.allFinite())
    return false;
  coeff = vec3(T.dot(S * T), T.dot(S * B), B.dot(S * B));
  return coeff.allFinite();
}

/// Area-weighted average of the incident face operators.
inline mat3 cr_vertex_operator(const std::vector<cr_face_operator> &faces,
                               const std::vector<int> &ids) {
  mat3 S = mat3::Zero();
  real w = 0;
  for (const int f : ids) {
    if (f < 0 || f >= static_cast<int>(faces.size()))
      continue;
    const cr_face_operator &face = faces[static_cast<size_t>(f)];
    if (!(face.area > real(0)) || !face.S.allFinite())
      continue;
    S += face.area * face.S;
    w += face.area;
  }
  if (!(w > real(0)))
    return mat3::Zero();
  return S / w;
}

/// Angle-weighted fan normal. `vertex` is the index shared by `faces`.
inline bool cr_fan_normal(const std::vector<vec3> &x,
                          const std::vector<cr_tri> &tris,
                          const std::vector<int> &faces, int vertex, vec3 &N) {
  N = vec3::Zero();
  if (vertex < 0 || vertex >= static_cast<int>(x.size()))
    return false;
  const vec3 &xv = x[static_cast<size_t>(vertex)];
  for (const int f : faces) {
    if (f < 0 || f >= static_cast<int>(tris.size()))
      continue;
    const cr_tri &tri = tris[static_cast<size_t>(f)];
    int k = -1;
    for (int i = 0; i < 3; ++i) {
      if (tri.v[static_cast<size_t>(i)] == vertex)
        k = i;
    }
    if (k < 0)
      continue;
    const vec3 cr = cr_face_cross(x, tri);
    const int ip = tri.v[static_cast<size_t>((k + 1) % 3)];
    const int in = tri.v[static_cast<size_t>((k + 2) % 3)];
    if (ip < 0 || in < 0 || ip >= static_cast<int>(x.size()) ||
        in >= static_cast<int>(x.size()))
      continue;
    vec3 e0 = x[static_cast<size_t>(ip)] - xv;
    vec3 e1 = x[static_cast<size_t>(in)] - xv;
    if (e0.squaredNorm() < real(1e-24) || e1.squaredNorm() < real(1e-24))
      continue;
    e0.normalize();
    e1.normalize();
    const real ede = std::clamp(e0.dot(e1), real(0), real(1) - real(1e-8));
    N += std::acos(ede) * cr;
  }
  if (!N.allFinite() || N.squaredNorm() < real(1e-24))
    return false;
  N.normalize();
  return true;
}

/// Red axis. `coeff` is parallel on [0, nE) and perpendicular on [nE, 2 nE).
/// `lock`, when set, flips T so it stays in the same half-plane.
inline bool cr_dirichlet_tangent(const std::vector<vec3> &x,
                                 const std::vector<cr_edge> &edges,
                                 const std::vector<real> &coeff,
                                 const std::vector<int> &edge_ids, const vec3 &N,
                                 vec3 &T, const vec3 *lock) {
  T = vec3::Zero();
  if (!N.allFinite() || N.squaredNorm() < real(1e-12))
    return false;
  const int nE = static_cast<int>(edges.size());
  if (static_cast<int>(coeff.size()) < 2 * nE)
    return false;
  const vec3 Nn = N.normalized();
  for (const int e : edge_ids) {
    if (e < 0 || e >= nE)
      continue;
    const cr_edge &ed = edges[static_cast<size_t>(e)];
    const vec3 e_par = cr_edge_parallel(x, ed);
    const vec3 W = coeff[static_cast<size_t>(e)] * e_par +
                   coeff[static_cast<size_t>(e + nE)] * cr_edge_perp(e_par, Nn);
    T += W - Nn * W.dot(Nn);
  }
  T -= T.dot(Nn) * Nn;
  if (!T.allFinite() || T.squaredNorm() < real(1e-16))
    return false;
  if (lock && lock->allFinite() && lock->squaredNorm() > real(1e-16) &&
      T.dot(*lock) < real(0))
    T = -T;
  T.normalize();
  return true;
}

/// Tangent chart of `(gp, gq)` in the plane of `N`. Parallel is
/// sign(v0 <= v1) * (p1 - p0).
inline bool diamond_chart(const vec3 &N, const vec3 &p0, const vec3 &p1, int v0,
                          int v1, real gp, real gq, vec3 &T, vec3 &B) {
  if (!N.allFinite() || !(N.squaredNorm() > real(1e-24)))
    return false;
  if (!std::isfinite(gp) || !std::isfinite(gq))
    return false;
  const vec3 Nn = N.normalized();
  vec3 e = p1 - p0;
  const real el = e.norm();
  if (!(el > real(1e-20)) || !e.allFinite())
    return false;
  const real s = v0 <= v1 ? real(1) : real(-1);
  const vec3 e_par = s * e / el;
  vec3 e_perp = Nn.cross(e_par);
  if (!(e_perp.squaredNorm() > real(1e-24)) || !e_perp.allFinite())
    return false;
  e_perp.normalize();
  const vec3 dir = gp * e_par + gq * e_perp;
  vec3 t = dir - Nn * dir.dot(Nn);
  if (!(t.squaredNorm() > real(1e-16)) || !t.allFinite())
    return false;
  t.normalize();
  vec3 b = Nn.cross(t);
  if (!(b.squaredNorm() > real(1e-16)) || !b.allFinite())
    return false;
  b.normalize();
  T = b.cross(Nn);
  B = b;
  return T.allFinite() && B.allFinite();
}

/// Area-weighted face normal of the two triangles. The cross uses the mesh
/// winding, so a sharp diamond still has one plane.
inline bool diamond_normal(const vec3 &p0, const vec3 &p1, const vec3 &a,
                           const vec3 &b, const std::array<int, 3> &t0,
                           const std::array<int, 3> &t1, int v0, int v1, int o0,
                           int o1, vec3 &N) {
  auto at = [&](int g) -> vec3 {
    if (g == v0)
      return p0;
    if (g == v1)
      return p1;
    if (g == o0)
      return a;
    if (g == o1)
      return b;
    return vec3::Constant(std::numeric_limits<real>::quiet_NaN());
  };
  vec3 nsum = vec3::Zero();
  const std::array<int, 3> *ts[2] = {&t0, &t1};
  for (int f = 0; f < 2; ++f) {
    const vec3 x0 = at((*ts[f])[0]);
    const vec3 x1 = at((*ts[f])[1]);
    const vec3 x2 = at((*ts[f])[2]);
    if (!x0.allFinite() || !x1.allFinite() || !x2.allFinite())
      return false;
    const vec3 c = (x1 - x0).cross(x2 - x0);
    if (!c.allFinite() || !(c.squaredNorm() > real(1e-24)))
      return false;
    nsum += c;
  }
  if (!nsum.allFinite() || !(nsum.squaredNorm() > real(1e-24)))
    return false;
  N = nsum.normalized();
  return true;
}

/// `(T, B, N)` of the diamond chart. The coefficients are not rebuilt here.
inline bool diamond_frame(const vec3 &p0, const vec3 &p1, const vec3 &a,
                          const vec3 &b, const std::array<int, 3> &t0,
                          const std::array<int, 3> &t1, int v0, int v1, int o0,
                          int o1, real gp, real gq, vec3 &N, vec3 &T, vec3 &B) {
  if (!diamond_normal(p0, p1, a, b, t0, t1, v0, v1, o0, o1, N))
    return false;
  return diamond_chart(N, p0, p1, v0, v1, gp, gq, T, B);
}

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif
