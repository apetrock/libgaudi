#ifndef GAUDI_KUSAMA_PRINCIPAL_CURVATURE_HPP
#define GAUDI_KUSAMA_PRINCIPAL_CURVATURE_HPP

/// Per-vertex principal curvatures and their position derivatives.
/// Eigensatz 2010, eq. (3.13) and Appendix A.1 (eq. A.1-A.16)
/// (`reference/eigensatz2010-constraint-based-surface-processing.pdf`).
/// The region B is the barycentric area of the vertex (the paper's smallest
/// scale), so only spokes contribute, each with half its length.

#include "gaudi/common.h"

#include <Eigen/Dense>

#include <array>
#include <cmath>
#include <map>
#include <utility>
#include <vector>

namespace gaudi {
namespace kusama {

using tri3 = std::array<int, 3>;

/// Edge `v -> j`. `f0` holds the directed edge `v -> j`, `f1` holds `j -> v`.
/// Either is -1 on a boundary edge.
struct principal_curvature_spoke {
  int j = -1;
  int f0 = -1;
  int f1 = -1;
};

struct principal_curvature_mesh {
  int n = 0;
  std::vector<tri3> tris;
  std::vector<std::vector<principal_curvature_spoke>> spokes;
  std::vector<char> boundary;

  static principal_curvature_mesh build(int n, const std::vector<tri3> &tris) {
    principal_curvature_mesh m;
    m.n = n;
    m.tris = tris;
    m.spokes.assign(static_cast<size_t>(n), {});
    m.boundary.assign(static_cast<size_t>(n), 0);
    std::map<std::pair<int, int>, int> directed;
    for (int f = 0; f < static_cast<int>(tris.size()); ++f)
      for (int k = 0; k < 3; ++k)
        directed[{tris[f][k], tris[f][(k + 1) % 3]}] = f;
    std::map<std::pair<int, int>, bool> seen;
    for (const auto &kv : directed) {
      const int a = kv.first.first;
      const int b = kv.first.second;
      for (const auto &ab : {std::make_pair(a, b), std::make_pair(b, a)}) {
        if (seen[ab])
          continue;
        seen[ab] = true;
        principal_curvature_spoke s;
        s.j = ab.second;
        auto f0 = directed.find(ab);
        auto f1 = directed.find({ab.second, ab.first});
        s.f0 = f0 == directed.end() ? -1 : f0->second;
        s.f1 = f1 == directed.end() ? -1 : f1->second;
        if (s.f0 < 0 || s.f1 < 0)
          m.boundary[static_cast<size_t>(ab.first)] = 1;
        m.spokes[static_cast<size_t>(ab.first)].push_back(s);
      }
    }
    return m;
  }
};

inline vec3 principal_curvature_face_cross(const std::vector<vec3> &x, const tri3 &t) {
  return (x[t[1]] - x[t[0]]).cross(x[t[2]] - x[t[0]]);
}

/// Vertex of `t` that is neither `a` nor `b`.
inline int principal_curvature_third(const tri3 &t, int a, int b) {
  for (int k = 0; k < 3; ++k)
    if (t[k] != a && t[k] != b)
      return t[k];
  return -1;
}

/// Signed dihedral angle across `a -> b`, positive on a convex crease.
/// `f0` holds `a -> b`, `f1` holds `b -> a`.
inline real principal_curvature_dihedral(const std::vector<vec3> &x,
                               const principal_curvature_mesh &m, int a, int b, int f0,
                               int f1) {
  const vec3 n0 = principal_curvature_face_cross(x, m.tris[f0]).normalized();
  const vec3 n1 = principal_curvature_face_cross(x, m.tris[f1]).normalized();
  const vec3 e = (x[b] - x[a]).normalized();
  return std::atan2(n0.cross(n1).dot(e), n0.dot(n1));
}

/// Gradient of `principal_curvature_dihedral` on the hinge `a -> b` with wings `c`
/// (in `f0`) and `d` (in `f1`), eq. (A.12-A.15). A wing moving along its
/// face normal turns that face about the edge by `1 / h`; the edge ends
/// take the lever-arm split, so the four gradients sum to zero.
struct principal_curvature_dihedral_grad {
  vec3 ga, gb, gc, gd;
};

inline principal_curvature_dihedral_grad
principal_curvature_dihedral_gradient(const vec3 &a, const vec3 &b, const vec3 &c,
                            const vec3 &d) {
  const vec3 e = b - a;
  const real l = e.norm();
  const vec3 N0 = e.cross(c - a);
  const vec3 N1 = (a - b).cross(d - b);
  principal_curvature_dihedral_grad g;
  g.gc = -(l / N0.squaredNorm()) * N0;
  g.gd = -(l / N1.squaredNorm()) * N1;
  const real t3 = (c - a).dot(e) / (l * l);
  const real t4 = (d - a).dot(e) / (l * l);
  g.ga = -(real(1) - t3) * g.gc - (real(1) - t4) * g.gd;
  g.gb = -t3 * g.gc - t4 * g.gd;
  return g;
}

/// `d A_f / d x_k` for corner `k` of `t`: `n̂ × (x_{k+2} - x_{k+1}) / 2`,
/// eq. (A.1).
inline vec3 principal_curvature_area_gradient(const std::vector<vec3> &x, const tri3 &t,
                                    int k) {
  const vec3 n = principal_curvature_face_cross(x, t).normalized();
  return real(0.5) * n.cross(x[t[(k + 2) % 3]] - x[t[(k + 1) % 3]]);
}

/// Barycentric area: a third of every incident face. Each incident face holds
/// exactly one edge leaving `v`, so it is counted once through `f0`.
inline real principal_curvature_vertex_area(const std::vector<vec3> &x,
                                  const principal_curvature_mesh &m, int v) {
  real A = 0;
  for (const principal_curvature_spoke &s : m.spokes[v])
    if (s.f0 >= 0)
      A += principal_curvature_face_cross(x, m.tris[s.f0]).norm() / real(6);
  return A;
}

inline vec3 principal_curvature_vertex_normal(const std::vector<vec3> &x,
                                    const principal_curvature_mesh &m, int v) {
  vec3 N = vec3::Zero();
  for (const principal_curvature_spoke &s : m.spokes[v])
    if (s.f0 >= 0)
      N += principal_curvature_face_cross(x, m.tris[s.f0]);
  return N.normalized();
}

/// `K(v) = 1/|B| sum_e beta(e) |e ∩ B| ebar ebarT`, eq. (3.13).
inline mat3 principal_curvature_tensor(const std::vector<vec3> &x,
                             const principal_curvature_mesh &m, int v) {
  mat3 K = mat3::Zero();
  const real A = principal_curvature_vertex_area(x, m, v);
  if (!(A > real(0)))
    return K;
  for (const principal_curvature_spoke &s : m.spokes[v]) {
    if (s.f0 < 0 || s.f1 < 0)
      continue;
    const vec3 e = x[s.j] - x[v];
    const real l = e.norm();
    if (!(l > real(0)))
      continue;
    const vec3 eb = e / l;
    const real beta = principal_curvature_dihedral(x, m, v, s.j, s.f0, s.f1);
    K += beta * (real(0.5) * l) * (eb * eb.transpose());
  }
  return K / A;
}

/// Vertices whose positions enter `K(v)`: `v` and its one-ring.
inline std::vector<int> principal_curvature_stencil(const principal_curvature_mesh &m, int v) {
  std::vector<int> ids{v};
  for (const principal_curvature_spoke &s : m.spokes[v])
    ids.push_back(s.j);
  return ids;
}

/// `dK / dx` for every stencil coordinate, eq. (A.7-A.16):
/// `dK = dKhat / |B| - Khat d|B| / |B|^2`, with `dKhat` from the dihedral,
/// half-length, and `ebar ebarT` derivatives of each spoke. Entry
/// `3 s + c` is coordinate `c` of `stencil[s]`.
inline bool principal_curvature_tensor_jacobian(const std::vector<vec3> &x,
                                      const principal_curvature_mesh &m, int v,
                                      const std::vector<int> &stencil,
                                      std::vector<mat3> &dK) {
  const int ns = static_cast<int>(stencil.size());
  dK.assign(static_cast<size_t>(3 * ns), mat3::Zero());
  std::map<int, int> local;
  for (int s = 0; s < ns; ++s)
    local[stencil[s]] = s;
  const real A = principal_curvature_vertex_area(x, m, v);
  if (!(A > real(0)))
    return false;
  mat3 Khat = mat3::Zero();
  std::vector<vec3> dA(static_cast<size_t>(ns), vec3::Zero());
  auto add = [&](int vert, const vec3 &g, const mat3 &M) {
    const auto it = local.find(vert);
    if (it == local.end())
      return false;
    for (int c = 0; c < 3; ++c)
      dK[3 * it->second + c] += g[c] * M;
    return true;
  };
  for (const principal_curvature_spoke &s : m.spokes[v]) {
    if (s.f0 >= 0) {
      const tri3 &t = m.tris[s.f0];
      for (int k = 0; k < 3; ++k) {
        const auto it = local.find(t[k]);
        if (it == local.end())
          return false;
        dA[it->second] += principal_curvature_area_gradient(x, t, k) / real(3);
      }
    }
    if (s.f0 < 0 || s.f1 < 0)
      continue;
    const vec3 e = x[s.j] - x[v];
    const real l = e.norm();
    if (!(l > real(0)))
      continue;
    const vec3 eb = e / l;
    const mat3 P = eb * eb.transpose();
    const real half = real(0.5) * l;
    const int c = principal_curvature_third(m.tris[s.f0], v, s.j);
    const int d = principal_curvature_third(m.tris[s.f1], v, s.j);
    const real beta = principal_curvature_dihedral(x, m, v, s.j, s.f0, s.f1);
    Khat += beta * half * P;
    const principal_curvature_dihedral_grad g =
        principal_curvature_dihedral_gradient(x[v], x[s.j], x[c], x[d]);
    if (!add(v, g.ga, half * P) || !add(s.j, g.gb, half * P) ||
        !add(c, g.gc, half * P) || !add(d, g.gd, half * P))
      return false;
    // d(l/2) = ebar . dx_j / 2, and d(ebar ebarT) along axis u of x_j is
    // (w ebarT + ebar wT) / l with w = (I - P) u; both flip sign on x_v.
    const int sj = local[s.j];
    const int sv = local[v];
    for (int ax = 0; ax < 3; ++ax) {
      vec3 u = vec3::Zero();
      u[ax] = real(1);
      const vec3 w = u - eb * eb[ax];
      const mat3 dP = (w * eb.transpose() + eb * w.transpose()) / l;
      const mat3 dterm = beta * (real(0.5) * eb[ax] * P + half * dP);
      dK[3 * sj + ax] += dterm;
      dK[3 * sv + ax] -= dterm;
    }
  }
  for (int s = 0; s < ns; ++s)
    for (int ax = 0; ax < 3; ++ax)
      dK[3 * s + ax] = dK[3 * s + ax] / A - Khat * (dA[s][ax] / (A * A));
  for (const mat3 &M : dK)
    if (!M.allFinite())
      return false;
  return true;
}

/// Same entries as `principal_curvature_tensor_jacobian`, by central differences of
/// `K`. `x` is scratch: entries are perturbed and restored.
inline bool principal_curvature_tensor_jacobian_fd(std::vector<vec3> &x,
                                         const principal_curvature_mesh &m, int v,
                                         const std::vector<int> &stencil,
                                         real eps, std::vector<mat3> &dK) {
  dK.assign(3 * stencil.size(), mat3::Zero());
  for (size_t s = 0; s < stencil.size(); ++s)
    for (int c = 0; c < 3; ++c) {
      const real x0 = x[stencil[s]][c];
      x[stencil[s]][c] = x0 + eps;
      const mat3 Kp = principal_curvature_tensor(x, m, v);
      x[stencil[s]][c] = x0 - eps;
      const mat3 Km = principal_curvature_tensor(x, m, v);
      x[stencil[s]][c] = x0;
      dK[3 * s + c] = (Kp - Km) / (real(2) * eps);
    }
  return true;
}

/// The paper's derivative cases, Appendix A.1.2.
enum class principal_curvature_case { distinct, cylindric, umbilic };

/// Eigenvalues of the 3x3 `K`: the one nearest zero is discarded, the larger
/// remaining signed one is `k1`, the smaller `k2` (text after eq. 3.13).
struct principal_curvatures {
  bool ok = false;
  real k1 = 0;
  real k2 = 0;
  vec3 u1 = vec3::Zero();
  vec3 u2 = vec3::Zero();
  principal_curvature_case kind = principal_curvature_case::distinct;
  /// In the cylindric case: which of `k1` / `k2` sits next to the discarded
  /// eigenvalue (0 or 1). The other one is the distinct eigenvalue.
  int near_zero = -1;
};

/// `gap_tol` scales with the largest |eigenvalue| and decides when two
/// eigenvalues count as equal (cylindric / umbilic).
inline principal_curvatures principal_curvature_at(const std::vector<vec3> &x,
                                                  const principal_curvature_mesh &m,
                                                  int v,
                                                  real gap_tol = real(1e-6)) {
  principal_curvatures out;
  if (m.boundary[v] || m.spokes[v].size() < 3)
    return out;
  const mat3 K = principal_curvature_tensor(x, m, v);
  if (!K.allFinite())
    return out;
  Eigen::SelfAdjointEigenSolver<mat3> es(real(0.5) * (K + K.transpose()));
  const vec3 lam = es.eigenvalues();
  int drop = 0;
  for (int i = 1; i < 3; ++i)
    if (std::abs(lam[i]) < std::abs(lam[drop]))
      drop = i;
  int keep[2], n = 0;
  for (int i = 0; i < 3; ++i)
    if (i != drop)
      keep[n++] = i;
  // eigenvalues are ascending, so keep[1] holds the larger one.
  out.k1 = lam[keep[1]];
  out.k2 = lam[keep[0]];
  out.u1 = es.eigenvectors().col(keep[1]);
  out.u2 = es.eigenvectors().col(keep[0]);
  const real scale = lam.cwiseAbs().maxCoeff() + real(1e-12);
  const real tol = gap_tol * scale;
  if (std::abs(out.k1 - out.k2) <= tol) {
    out.kind = principal_curvature_case::umbilic;
  } else {
    const real g1 = std::abs(out.k1 - lam[drop]);
    const real g2 = std::abs(out.k2 - lam[drop]);
    if (std::min(g1, g2) <= tol) {
      out.kind = principal_curvature_case::cylindric;
      out.near_zero = g1 < g2 ? 0 : 1;
    }
  }
  out.ok = std::isfinite(out.k1) && std::isfinite(out.k2);
  return out;
}

/// `d kappa / dx` from `dK`, eq. (A.4-A.6). Distinct: `u_jT dK u_j`.
/// Cylindric: the distinct eigenvalue by (A.4), the one next to the
/// discarded eigenvalue by `2 dH - d kappa_distinct`, `dH = tr(dK) / 2`
/// (A.5). Umbilic: both `dH` (A.6). Rows: kappa1, kappa2.
inline bool principal_curvature_jacobian(const principal_curvatures &pr,
                                         const std::vector<mat3> &dK,
                                         Eigen::MatrixXd &J) {
  if (!pr.ok)
    return false;
  const int nc = static_cast<int>(dK.size());
  J = Eigen::MatrixXd::Zero(2, nc);
  for (int i = 0; i < nc; ++i) {
    const mat3 D = real(0.5) * (dK[i] + dK[i].transpose());
    const real dH = real(0.5) * D.trace();
    const real d1 = pr.u1.dot(D * pr.u1);
    const real d2 = pr.u2.dot(D * pr.u2);
    switch (pr.kind) {
    case principal_curvature_case::distinct:
      J(0, i) = d1;
      J(1, i) = d2;
      break;
    case principal_curvature_case::cylindric:
      if (pr.near_zero == 1) {
        J(0, i) = d1;
        J(1, i) = real(2) * dH - d1;
      } else {
        J(1, i) = d2;
        J(0, i) = real(2) * dH - d2;
      }
      break;
    case principal_curvature_case::umbilic:
      J(0, i) = dH;
      J(1, i) = dH;
      break;
    }
  }
  return J.allFinite();
}

} // namespace kusama
} // namespace gaudi

#endif
