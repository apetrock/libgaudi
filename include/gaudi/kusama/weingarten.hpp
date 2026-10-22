#ifndef GAUDI_KUSAMA_WEINGARTEN_HPP
#define GAUDI_KUSAMA_WEINGARTEN_HPP

/// Weingarten (shape-operator) curvature solve on the CR edge element.
///
/// The curvature tensor lives on edges: `h_e = (k11, k12, k22)` in the edge's
/// unfolded frame `(t_e, N x t_e)`. Positions couple to it through the mean
/// of the edge's two faces' live CR read-outs (Stein et al. 2020, App. A.2,
/// with the CR derivative recomputed from each face's stencil), the metric
/// is held by Stein's vector Dirichlet energy on the intrinsic edge frames,
/// and the tensor field is faired by the CR spin-2 element (Stein's
/// construction applied to symmetric tensors). Derivations and experiments:
/// `scripts/sympy/cr_poisson/MATH.md` sections 11 to 14.

#include "gaudi/asawa/shell/face_quality.hpp"
#include "gaudi/common.h"
#include "gaudi/kusama/stencil_table.hpp"
#include "gaudi/kusama/vector_dirichlet.hpp"

#include "cr_curvature_edge_generated.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace gaudi {
namespace kusama {
namespace weingarten {

using matS = Eigen::SparseMatrix<real>;
using trip = Eigen::Triplet<real>;

// ------------------------------------------------------------ CR derivative

/// One edge's CR derivative of the embedding, `g_par = sum par_w x` and
/// `g_perp = sum perp_w x`, frozen at the rest metric. Slots: the edge's
/// lower and higher vertex ids, then the two opposite vertices (-1 absent).
/// `g_par` runs from the lower id to the higher (the canonical orientation).
/// Faces below this quality (`face_quality`) are degenerate: their CR
/// read-outs, Jacobians and cotangent stiffness divide by vanishing areas
/// and heights.
constexpr real k_degenerate_quality = real(1e-4);

using asawa::shell::face_quality;

struct cr_edge_rows {
  std::array<int, 4> v{-1, -1, -1, -1};
  std::array<real, 4> par{0, 0, 0, 0};
  std::array<real, 4> perp{0, 0, 0, 0};
  real mass = 0; // A_e / 3
};

inline std::vector<cr_edge_rows> build_cr_derivative(const stencil_table &t,
                                                     const std::vector<vec3> &x) {
  std::vector<cr_edge_rows> rows(static_cast<size_t>(t.n_edges));
  std::vector<real> area_sum(static_cast<size_t>(t.n_edges), real(0));
  auto slot = [](cr_edge_rows &r, int vid) {
    for (int k = 0; k < 4; ++k) {
      if (r.v[static_cast<size_t>(k)] == vid)
        return k;
      if (r.v[static_cast<size_t>(k)] < 0) {
        r.v[static_cast<size_t>(k)] = vid;
        return k;
      }
    }
    return -1;
  };
  for (int f = 0; f < t.n_faces; ++f) {
    const auto &id = t.ids[static_cast<size_t>(f)];
    const vec3 &a = x[static_cast<size_t>(id[0])];
    const vec3 &b = x[static_cast<size_t>(id[1])];
    const vec3 &c = x[static_cast<size_t>(id[2])];
    const real A = real(0.5) * (b - a).cross(c - a).norm();
    if (!(A > real(1e-30)))
      continue;
    for (int i = 0; i < 3; ++i) {
      const int e = t.face_edges[static_cast<size_t>(f)][static_cast<size_t>(i)];
      if (e < 0)
        continue;
      const int tail = id[static_cast<size_t>(i)];
      const int head = id[static_cast<size_t>((i + 1) % 3)];
      const int opp = id[static_cast<size_t>((i + 2) % 3)];
      const vec3 &xt = x[static_cast<size_t>(tail)];
      const vec3 &xh = x[static_cast<size_t>(head)];
      const vec3 &xo = x[static_cast<size_t>(opp)];
      const real si = tail < head ? real(1) : real(-1);
      const real l = (xh - xt).norm();
      if (!(l > real(1e-30)))
        continue;
      const real ejk = (xo - xh).squaredNorm();
      const real eki = (xt - xo).squaredNorm();
      cr_edge_rows &r = rows[static_cast<size_t>(e)];
      if (r.v[0] < 0) {
        r.v[0] = std::min(tail, head);
        r.v[1] = std::max(tail, head);
      }
      area_sum[static_cast<size_t>(e)] += A;
      const real pa = si * A / (real(3) * l);
      r.par[static_cast<size_t>(slot(r, head))] += pa;
      r.par[static_cast<size_t>(slot(r, tail))] -= pa;
      const real k = -si / (real(12) * l);
      r.perp[static_cast<size_t>(slot(r, tail))] += k * (l * l - (eki - ejk));
      r.perp[static_cast<size_t>(slot(r, head))] += k * (l * l + (eki - ejk));
      r.perp[static_cast<size_t>(slot(r, opp))] += si * l / real(6);
    }
  }
  for (int e = 0; e < t.n_edges; ++e) {
    cr_edge_rows &r = rows[static_cast<size_t>(e)];
    r.mass = area_sum[static_cast<size_t>(e)] / real(3);
    if (!(r.mass > real(0)))
      continue;
    for (int k = 0; k < 4; ++k) {
      r.par[static_cast<size_t>(k)] /= r.mass;
      r.perp[static_cast<size_t>(k)] /= r.mass;
    }
  }
  return rows;
}

inline void edge_tangents(const cr_edge_rows &r, const std::vector<vec3> &x, vec3 &gp,
                          vec3 &gq) {
  gp = vec3::Zero();
  gq = vec3::Zero();
  for (int k = 0; k < 4; ++k) {
    const int v = r.v[static_cast<size_t>(k)];
    if (v < 0)
      continue;
    gp += r.par[static_cast<size_t>(k)] * x[static_cast<size_t>(v)];
    gq += r.perp[static_cast<size_t>(k)] * x[static_cast<size_t>(v)];
  }
}

// ------------------------------------------------------------ edge readings

/// Mean of each edge's faces' live read-outs in the edge frame, minus the
/// target `h`, and its Jacobian. Rows `3 e + c`, scaled by `sqrt(mass_e)` so
/// the energy is `1/2 |r|^2 = 1/2 sum_e A_e / 3 |mean_e - h_e|^2`. With
/// `h == nullptr` the readings themselves (unscaled) are returned in `r`.
inline void edge_coupling(const stencil_table &t, const std::vector<vec3> &x,
                          const std::vector<real> &mass, const std::vector<vec3> *h,
                          Eigen::VectorXd &r, std::vector<trip> *J, index_t row0 = 0,
                          real weight = real(1)) {
  const int nE = t.n_edges;
  r = Eigen::VectorXd::Zero(3 * nE);
  std::vector<int> count(static_cast<size_t>(nE), 0);
  for (int f = 0; f < t.n_faces; ++f)
    for (int i = 0; i < 3; ++i)
      ++count[static_cast<size_t>(t.face_edges[static_cast<size_t>(f)][static_cast<size_t>(i)])];
  auto scale = [&](int e) {
    const real m = h ? std::sqrt(weight * std::max(mass[static_cast<size_t>(e)], real(0)))
                     : real(1);
    return m / real(std::max(count[static_cast<size_t>(e)], 1));
  };
  for (int f = 0; f < t.n_faces; ++f) {
    vec3 hf;
    const std::array<vec3, 6> P = stencil_positions(t, x, f, hf);
    Eigen::Matrix<real, 9, 1> R;
    Eigen::Matrix<real, 9, 18> JR;
    hepworth::block::cr_generated::cr_curvature_edge_readings(P, hf, R, JR);
    if (!R.allFinite() || !JR.allFinite())
      continue;
    const auto &id = t.ids[static_cast<size_t>(f)];
    for (int i = 0; i < 3; ++i) {
      const int e = t.face_edges[static_cast<size_t>(f)][static_cast<size_t>(i)];
      const real s = scale(e);
      r.segment<3>(3 * e) += s * R.segment<3>(3 * i);
      if (!J)
        continue;
      for (int k = 0; k < 6; ++k) {
        if (k >= 3 && hf[k - 3] == real(0))
          continue;
        const int v = id[static_cast<size_t>(k)];
        for (int c = 0; c < 3; ++c)
          for (int ax = 0; ax < 3; ++ax) {
            const real val = s * JR(3 * i + c, 3 * k + ax);
            if (val != real(0))
              J->emplace_back(row0 + 3 * e + c, 3 * v + ax, val);
          }
      }
    }
  }
  if (h)
    for (int e = 0; e < nE; ++e)
      r.segment<3>(3 * e) -= scale(e) * real(std::max(count[static_cast<size_t>(e)], 1)) *
                             (*h)[static_cast<size_t>(e)];
}

/// The edge store's natural value at `x`: the mean of each edge's faces'
/// readings in its frame (the coupling is idle there).
inline std::vector<vec3> edge_readings(const stencil_table &t, const std::vector<vec3> &x) {
  Eigen::VectorXd r;
  const std::vector<real> unit(static_cast<size_t>(t.n_edges), real(1));
  edge_coupling(t, x, unit, nullptr, r, nullptr);
  std::vector<vec3> out(static_cast<size_t>(t.n_edges));
  for (int e = 0; e < t.n_edges; ++e)
    out[static_cast<size_t>(e)] = r.segment<3>(3 * e);
  return out;
}

// ------------------------------------------------------- intrinsic memory

/// Intrinsic frame of the CR frame `(g_par, g_perp)` in its own basis
/// `t = g_par / |g_par|`, `b = n x t`: `U = [[|g_par|, g_par.g_perp / |g_par|],
/// [0, |g_par x g_perp| / |g_par|]]`, returned as `(U00, U01, U11)` with
/// `dU / d(g_par, g_perp)` (3 x 6). Rotation-invariant; an isometric bend
/// changes it only through the `|g_perp|` shortening.
inline vec3 intrinsic_frame(const vec3 &gp, const vec3 &gq,
                            Eigen::Matrix<real, 3, 6> *dU = nullptr) {
  const real a = gp.norm();
  const vec3 w = gp.cross(gq);
  const real wn = w.norm();
  if (!(a > real(1e-30)) || !(wn > real(1e-30))) {
    if (dU)
      dU->setZero();
    return vec3(a, real(0), real(0));
  }
  const real pq = gp.dot(gq);
  if (dU) {
    const vec3 wh = w / wn;
    dU->setZero();
    dU->block<1, 3>(0, 0) = (gp / a).transpose();
    dU->block<1, 3>(1, 0) = (gq / a - pq * gp / (a * a * a)).transpose();
    dU->block<1, 3>(1, 3) = (gp / a).transpose();
    dU->block<1, 3>(2, 0) = (gq.cross(wh) / a - wn * gp / (a * a * a)).transpose();
    dU->block<1, 3>(2, 3) = (wh.cross(gp) / a).transpose();
  }
  return vec3(a, pq / a, wn / a);
}

/// Rest data of the memory: the CR derivative, the rest intrinsic frames,
/// and per face `K_f^{1/2}` of Stein's vector Dirichlet element (6 x 6 over
/// the face's three edges, `comp * 3 + edge`), plus the edge mass term
/// `mu M` for the element's null space.
struct memory_rest {
  std::vector<cr_edge_rows> D;
  std::vector<vec3> U0;
  std::vector<Eigen::Matrix<real, 6, 6>> K_half;
  real mu = 0;
};

inline memory_rest build_memory_rest(const stencil_table &t, const std::vector<vec3> &x,
                                     real mu_rel = real(1e-3)) {
  memory_rest m;
  m.D = build_cr_derivative(t, x);
  m.U0.resize(static_cast<size_t>(t.n_edges));
  for (int e = 0; e < t.n_edges; ++e) {
    vec3 gp, gq;
    edge_tangents(m.D[static_cast<size_t>(e)], x, gp, gq);
    m.U0[static_cast<size_t>(e)] = intrinsic_frame(gp, gq);
  }
  m.K_half.resize(static_cast<size_t>(t.n_faces));
  real trK = 0, trM = 0;
  for (int f = 0; f < t.n_faces; ++f) {
    const auto &id = t.ids[static_cast<size_t>(f)];
    // A degenerate face's cotangents are unbounded: no stiffness, and kept
    // out of `trK`, which scales every edge's mass row through `mu`.
    if (!(face_quality(x[static_cast<size_t>(id[0])], x[static_cast<size_t>(id[1])],
                       x[static_cast<size_t>(id[2])]) >= k_degenerate_quality)) {
      m.K_half[static_cast<size_t>(f)].setZero();
      continue;
    }
    std::array<int, 3> sign{};
    for (int i = 0; i < 3; ++i)
      sign[static_cast<size_t>(i)] =
          id[static_cast<size_t>(i)] <= id[static_cast<size_t>((i + 1) % 3)] ? 1 : -1;
    const Eigen::Matrix<real, 6, 6> K = laplace::vector_dirichlet::face_element(
        x[static_cast<size_t>(id[0])], x[static_cast<size_t>(id[1])],
        x[static_cast<size_t>(id[2])], sign);
    trK += K.trace();
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<real, 6, 6>> es(K);
    const Eigen::Matrix<real, 6, 1> ev = es.eigenvalues().cwiseMax(real(0)).cwiseSqrt();
    m.K_half[static_cast<size_t>(f)] = es.eigenvectors() * ev.asDiagonal() *
                                      es.eigenvectors().transpose();
  }
  for (const cr_edge_rows &r : m.D)
    trM += real(2) * r.mass;
  m.mu = trM > real(0) ? mu_rel * trK / trM : real(0);
  return m;
}

/// Per-edge intrinsic-frame change `c_e = U(D x) - U0` and `dc_e / dx` per
/// slot vertex (3 x 3 blocks: component x coordinate).
struct memory_edge_eval {
  vec3 c = vec3::Zero();
  std::array<Eigen::Matrix<real, 3, 3>, 4> dc{};
};

inline memory_edge_eval memory_edge(const memory_rest &m, int e, const std::vector<vec3> &x) {
  memory_edge_eval out;
  const cr_edge_rows &r = m.D[static_cast<size_t>(e)];
  vec3 gp, gq;
  edge_tangents(r, x, gp, gq);
  Eigen::Matrix<real, 3, 6> dU;
  out.c = intrinsic_frame(gp, gq, &dU) - m.U0[static_cast<size_t>(e)];
  for (int k = 0; k < 4; ++k)
    out.dc[static_cast<size_t>(k)] = r.par[static_cast<size_t>(k)] * dU.block<3, 3>(0, 0) +
                                     r.perp[static_cast<size_t>(k)] * dU.block<3, 3>(0, 3);
  return out;
}

/// Memory residual `r = [K_f^{1/2} c_f per face and column; sqrt(mu M_e) c_e]`
/// (energy `1/2 |r|^2 = 1/2 sum_j c_j^T (L_dir + mu M) c_j`). Columns of `U`:
/// column 0 is `(U00, U10 = 0)`, column 1 is `(U01, U11)`; the zero component
/// is constant and contributes nothing, so it is skipped.
inline void memory_residual(const stencil_table &t, const memory_rest &m,
                            const std::vector<vec3> &x, real weight, Eigen::VectorXd &r,
                            std::vector<trip> *J, index_t row0 = 0) {
  const real sw = std::sqrt(weight);
  std::vector<memory_edge_eval> ev(static_cast<size_t>(t.n_edges));
  for (int e = 0; e < t.n_edges; ++e)
    ev[static_cast<size_t>(e)] = memory_edge(m, e, x);
  // Rows: per face 12 (two columns of 6), per edge 3 (mass rows on U00, U01, U11).
  r = Eigen::VectorXd::Zero(12 * t.n_faces + 3 * t.n_edges);
  // c components per edge per column: column 0 = (U00, 0), column 1 = (U01, U11).
  auto comp = [](int col, int part) { return col == 0 ? (part == 0 ? 0 : -1) : (part == 0 ? 1 : 2); };
  for (int f = 0; f < t.n_faces; ++f) {
    const auto &K = m.K_half[static_cast<size_t>(f)];
    for (int col = 0; col < 2; ++col) {
      Eigen::Matrix<real, 6, 1> c = Eigen::Matrix<real, 6, 1>::Zero();
      for (int part = 0; part < 2; ++part) {
        const int u = comp(col, part);
        if (u < 0)
          continue;
        for (int i = 0; i < 3; ++i) {
          const int e = t.face_edges[static_cast<size_t>(f)][static_cast<size_t>(i)];
          c[3 * part + i] = ev[static_cast<size_t>(e)].c[u];
        }
      }
      const index_t base = 12 * f + 6 * col;
      r.segment<6>(base) = sw * (K * c);
      if (!J)
        continue;
      for (int part = 0; part < 2; ++part) {
        const int u = comp(col, part);
        if (u < 0)
          continue;
        for (int i = 0; i < 3; ++i) {
          const int e = t.face_edges[static_cast<size_t>(f)][static_cast<size_t>(i)];
          const cr_edge_rows &rows = m.D[static_cast<size_t>(e)];
          for (int k = 0; k < 4; ++k) {
            const int v = rows.v[static_cast<size_t>(k)];
            if (v < 0)
              continue;
            const Eigen::Matrix<real, 1, 3> g =
                ev[static_cast<size_t>(e)].dc[static_cast<size_t>(k)].row(u);
            for (int q = 0; q < 6; ++q) {
              const real kq = sw * K(q, 3 * part + i);
              if (kq == real(0))
                continue;
              for (int ax = 0; ax < 3; ++ax)
                if (g[ax] != real(0))
                  J->emplace_back(row0 + base + q, 3 * v + ax, kq * g[ax]);
            }
          }
        }
      }
    }
  }
  const index_t mrow = 12 * t.n_faces;
  for (int e = 0; e < t.n_edges; ++e) {
    const real s = sw * std::sqrt(m.mu * std::max(m.D[static_cast<size_t>(e)].mass, real(0)));
    r.segment<3>(mrow + 3 * e) = s * ev[static_cast<size_t>(e)].c;
    if (!J)
      continue;
    const cr_edge_rows &rows = m.D[static_cast<size_t>(e)];
    for (int k = 0; k < 4; ++k) {
      const int v = rows.v[static_cast<size_t>(k)];
      if (v < 0)
        continue;
      for (int u = 0; u < 3; ++u)
        for (int ax = 0; ax < 3; ++ax) {
          const real val = s * ev[static_cast<size_t>(e)].dc[static_cast<size_t>(k)](u, ax);
          if (val != real(0))
            J->emplace_back(row0 + mrow + 3 * e + u, 3 * v + ax, val);
        }
    }
  }
}

// ------------------------------------------------------ spin-2 smoothness

/// `(k11, k12, k22)` in a frame rotated by `phi` (cos `c`, sin `s`) -> the same
/// tensor in the unrotated frame.
inline Eigen::Matrix<real, 3, 3> tensor_map(real c, real s) {
  Eigen::Matrix<real, 3, 3> M;
  M << c * c, -2 * c * s, s * s, c * s, c * c - s * s, -c * s, s * s, 2 * c * s, c * c;
  return M;
}

/// `q(h, h') = z . z' + 1/2 tau tau'`, `tau = k11 + k22`, `z = ((k11 - k22) / 2, k12)`.
inline Eigen::Matrix<real, 3, 3> tensor_inner() {
  Eigen::Matrix<real, 3, 3> Q;
  Q << real(0.75), 0, real(0.25), 0, 1, 0, real(0.25), 0, real(0.75);
  return Q;
}

/// CR scalar stiffness of one triangle (3 x 3), indexed by local edge
/// (4 x the P1 stiffness of the opposite vertices).
inline Eigen::Matrix<real, 3, 3> cr_face_stiffness(const vec3 &X0, const vec3 &X1,
                                                   const vec3 &X2) {
  const vec3 X[3] = {X0, X1, X2};
  real cot[3];
  for (int k = 0; k < 3; ++k) {
    const vec3 u = X[(k + 1) % 3] - X[k], v = X[(k + 2) % 3] - X[k];
    cot[k] = u.dot(v) / u.cross(v).norm();
  }
  Eigen::Matrix<real, 3, 3> K1 = Eigen::Matrix<real, 3, 3>::Zero();
  for (int k = 0; k < 3; ++k) {
    const int i = (k + 1) % 3, j = (k + 2) % 3;
    K1(i, j) = K1(j, i) = real(-0.5) * cot[k];
  }
  for (int i = 0; i < 3; ++i)
    K1(i, i) = -K1.row(i).sum();
  Eigen::Matrix<real, 3, 3> K;
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b)
      K(a, b) = real(4) * K1((a + 2) % 3, (b + 2) % 3);
  return K;
}

/// CR spin-2 element of one triangle (9 x 9 over its three edge tensors,
/// `3 edge + comp`): `K_ab M_a^T Q M_b`, `M_a` the map of edge `a`'s frame
/// into the face frame. Energy `h^T E h` approximates
/// `int |grad z|^2 + 1/2 |grad tau|^2` (second order; MATH.md section 14).
inline Eigen::Matrix<real, 9, 9> spin2_face_element(const vec3 &X0, const vec3 &X1,
                                                    const vec3 &X2) {
  const vec3 X[3] = {X0, X1, X2};
  const Eigen::Matrix<real, 3, 3> K = cr_face_stiffness(X0, X1, X2);
  const vec3 N = (X1 - X0).cross(X2 - X0).normalized();
  const vec3 T = (X1 - X0).normalized();
  const vec3 B = N.cross(T);
  std::array<Eigen::Matrix<real, 3, 3>, 3> M;
  for (int a = 0; a < 3; ++a) {
    const vec3 v = (X[(a + 1) % 3] - X[a]).normalized();
    M[static_cast<size_t>(a)] = tensor_map(v.dot(T), v.dot(B));
  }
  const Eigen::Matrix<real, 3, 3> Q = tensor_inner();
  Eigen::Matrix<real, 9, 9> E;
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b)
      E.block<3, 3>(3 * a, 3 * b) =
          K(a, b) * M[static_cast<size_t>(a)].transpose() * Q * M[static_cast<size_t>(b)];
  return E;
}

/// Assembled spin-2 energy matrix (3 nE x 3 nE), `h^T L h`.
inline matS build_spin2(const stencil_table &t, const std::vector<vec3> &x) {
  std::vector<trip> trips;
  trips.reserve(static_cast<size_t>(81 * t.n_faces));
  for (int f = 0; f < t.n_faces; ++f) {
    const auto &id = t.ids[static_cast<size_t>(f)];
    const Eigen::Matrix<real, 9, 9> E = spin2_face_element(
        x[static_cast<size_t>(id[0])], x[static_cast<size_t>(id[1])],
        x[static_cast<size_t>(id[2])]);
    if (!E.allFinite())
      continue;
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b) {
        const int ea = t.face_edges[static_cast<size_t>(f)][static_cast<size_t>(a)];
        const int eb = t.face_edges[static_cast<size_t>(f)][static_cast<size_t>(b)];
        for (int r = 0; r < 3; ++r)
          for (int c = 0; c < 3; ++c)
            trips.emplace_back(3 * ea + r, 3 * eb + c, E(3 * a + r, 3 * b + c));
      }
  }
  matS L(3 * t.n_edges, 3 * t.n_edges);
  L.setFromTriplets(trips.begin(), trips.end());
  return L;
}

// ------------------------------------------------------- row factor

/// Calibrates the curvature rows against the cotan `bend` row: at equal
/// weights their mean position stiffness matches (measured 9.0 on the demo
/// torus, `weingarten_coupling_matches_bend_at_equal_weight`). Both are
/// dimensionless, so it holds across mesh size and resolution.
constexpr real k_coupling = real(9);

/// Faces of quality below `q_min` (`n_faces`, `worst_face`), and what their
/// stencils reach: `face` marks them and every face sharing an edge with one
/// (a face's CR read-out uses the wing vertices across its three edges), and
/// `edge` the edges of those faces.
struct degenerate_set {
  std::vector<char> face, edge;
  /// Vertices of the degenerate faces themselves.
  std::vector<char> vert;
  int n_faces = 0, n_edges = 0, worst_face = -1;
  real worst = 1;
};

inline degenerate_set find_degenerate(const stencil_table &t, const std::vector<vec3> &x,
                                      real q_min = k_degenerate_quality) {
  degenerate_set d;
  const size_t nF = static_cast<size_t>(std::max(t.n_faces, 0));
  const size_t nE = static_cast<size_t>(std::max(t.n_edges, 0));
  d.face.assign(nF, 0);
  d.edge.assign(nE, 0);
  d.vert.assign(static_cast<size_t>(std::max(t.n_verts, 0)), 0);
  std::vector<char> core(nF, 0), core_edge(nE, 0);
  for (int f = 0; f < t.n_faces; ++f) {
    const auto &id = t.ids[static_cast<size_t>(f)];
    const real q = face_quality(x[static_cast<size_t>(id[0])], x[static_cast<size_t>(id[1])],
                                x[static_cast<size_t>(id[2])]);
    if (q < d.worst) {
      d.worst = q;
      d.worst_face = f;
    }
    if (!(q >= q_min)) {
      core[static_cast<size_t>(f)] = 1;
      ++d.n_faces;
      for (int i = 0; i < 3; ++i)
        if (id[static_cast<size_t>(i)] >= 0 && id[static_cast<size_t>(i)] < t.n_verts)
          d.vert[static_cast<size_t>(id[static_cast<size_t>(i)])] = 1;
      for (int i = 0; i < 3; ++i) {
        const int e = t.face_edges[static_cast<size_t>(f)][static_cast<size_t>(i)];
        if (e >= 0 && static_cast<size_t>(e) < nE)
          core_edge[static_cast<size_t>(e)] = 1;
      }
    }
  }
  if (d.n_faces == 0)
    return d;
  for (int f = 0; f < t.n_faces; ++f) {
    bool touched = core[static_cast<size_t>(f)];
    for (int i = 0; i < 3 && !touched; ++i) {
      const int e = t.face_edges[static_cast<size_t>(f)][static_cast<size_t>(i)];
      touched = e >= 0 && static_cast<size_t>(e) < nE && core_edge[static_cast<size_t>(e)];
    }
    if (!touched)
      continue;
    d.face[static_cast<size_t>(f)] = 1;
    for (int i = 0; i < 3; ++i) {
      const int e = t.face_edges[static_cast<size_t>(f)][static_cast<size_t>(i)];
      if (e >= 0 && static_cast<size_t>(e) < nE && !d.edge[static_cast<size_t>(e)]) {
        d.edge[static_cast<size_t>(e)] = 1;
        ++d.n_edges;
      }
    }
  }
  return d;
}

/// The curvature rows' factor per edge: `k_coupling A_e`, `A_e` the edge's
/// dual area (its incident face areas over 3). A curvature times `A_e` is an
/// integrated curvature, the quantity `bend` weighs (`L x`).
inline std::vector<real> coupling_factor(const stencil_table &t, const std::vector<vec3> &x) {
  std::vector<real> c(static_cast<size_t>(t.n_edges), real(0));
  for (int f = 0; f < t.n_faces; ++f) {
    const auto &id = t.ids[static_cast<size_t>(f)];
    const real A = real(0.5) * (x[static_cast<size_t>(id[1])] - x[static_cast<size_t>(id[0])])
                                   .cross(x[static_cast<size_t>(id[2])] - x[static_cast<size_t>(id[0])])
                                   .norm();
    for (int i = 0; i < 3; ++i)
      c[static_cast<size_t>(t.face_edges[f][i])] += k_coupling * A / real(3);
  }
  return c;
}

// -------------------------------------------------------- edge to vertex

/// Incident edges of each vertex with their weights (`A_e / 3`) and
/// canonical endpoints (lower id first). The local polar of edge to vertex.
struct vertex_stars {
  struct spoke {
    int edge = -1;
    int lo = -1;
    int hi = -1;
    real w = 0;
  };
  std::vector<std::vector<spoke>> at;
};

inline vertex_stars build_vertex_stars(const stencil_table &t, const std::vector<vec3> &x) {
  vertex_stars s;
  s.at.resize(static_cast<size_t>(t.n_verts));
  std::vector<real> mass(static_cast<size_t>(t.n_edges), real(0));
  std::vector<std::array<int, 2>> ends(static_cast<size_t>(t.n_edges), {-1, -1});
  for (int f = 0; f < t.n_faces; ++f) {
    const auto &id = t.ids[static_cast<size_t>(f)];
    const real A = real(0.5) * (x[static_cast<size_t>(id[1])] - x[static_cast<size_t>(id[0])])
                                   .cross(x[static_cast<size_t>(id[2])] - x[static_cast<size_t>(id[0])])
                                   .norm();
    for (int i = 0; i < 3; ++i) {
      const int e = t.face_edges[static_cast<size_t>(f)][static_cast<size_t>(i)];
      const int a = id[static_cast<size_t>(i)], b = id[static_cast<size_t>((i + 1) % 3)];
      mass[static_cast<size_t>(e)] += A / real(3);
      ends[static_cast<size_t>(e)] = {std::min(a, b), std::max(a, b)};
    }
  }
  for (int e = 0; e < t.n_edges; ++e) {
    const auto [lo, hi] = ends[static_cast<size_t>(e)];
    if (lo < 0)
      continue;
    const vertex_stars::spoke sp{e, lo, hi, mass[static_cast<size_t>(e)]};
    s.at[static_cast<size_t>(lo)].push_back(sp);
    s.at[static_cast<size_t>(hi)].push_back(sp);
  }
  return s;
}

/// `(k11, k12, k22)` at vertex `v` in its chart `(T, B, N)`: the weighted mean
/// of the incident edge tensors, each placed in the vertex plane by its
/// projected edge direction. `pos(i)` and `store(e)` read positions and the
/// edge store (a vector or a live `q`).
template <typename Pos, typename Store>
vec3 vertex_tensor(const vertex_stars &s, int v, Pos &&pos, Store &&store, const vec3 &T,
                   const vec3 &B, const vec3 &N) {
  mat3 S = mat3::Zero();
  real wsum = 0;
  for (const vertex_stars::spoke &sp : s.at[static_cast<size_t>(v)]) {
    vec3 t = pos(sp.hi) - pos(sp.lo);
    t -= N * t.dot(N);
    const real l = t.norm();
    const vec3 h = store(sp.edge);
    if (!(l > real(1e-20)) || !h.allFinite() || !(sp.w > real(0)))
      continue;
    t /= l;
    const vec3 b = N.cross(t);
    S += sp.w * (h[0] * t * t.transpose() + h[1] * (t * b.transpose() + b * t.transpose()) +
                 h[2] * b * b.transpose());
    wsum += sp.w;
  }
  if (!(wsum > real(0)))
    return vec3::Constant(std::numeric_limits<real>::quiet_NaN());
  S /= wsum;
  return vec3(T.dot(S * T), T.dot(S * B), B.dot(S * B));
}

// ------------------------------------------------------------------ solver

struct solver_report {
  real energy0 = 0;
  real energy1 = 0;
  real step = 0;
  int halvings = 0;
  bool accepted = false;
  bool refactored = false;
};

/// Uncoupled Weingarten solve over positions: edge tensor targets `h`, the
/// intrinsic vector Dirichlet memory, frozen-factor Gauss–Newton with
/// Armijo backtracking, and the restart protocol (`solve`: converge, re-rest,
/// repeat). `w_memory` is relative: the memory Hessian is scaled so its
/// trace matches `D^T M D`, as in the Python study.
class weingarten_solver {
public:
  real w_memory = real(1e-3);
  real mu_rel = real(1e-3);
  real eps = real(1e-10);
  int backtracks = 12;

  void reset(const stencil_table &t, const std::vector<vec3> &x_rest) {
    _t = t;
    reset_rest(x_rest);
  }

  /// Restart: metric, read-out weights and memory rest move to `x`.
  void reset_rest(const std::vector<vec3> &x) {
    _m = build_memory_rest(_t, x, mu_rel);
    _mass.assign(static_cast<size_t>(_t.n_edges), real(0));
    for (int e = 0; e < _t.n_edges; ++e)
      _mass[static_cast<size_t>(e)] = _m.D[static_cast<size_t>(e)].mass;
    real trL = 0;
    for (const cr_edge_rows &r : _m.D)
      for (int k = 0; k < 4; ++k)
        trL += real(3) * r.mass *
               (r.par[static_cast<size_t>(k)] * r.par[static_cast<size_t>(k)] +
                r.perp[static_cast<size_t>(k)] * r.perp[static_cast<size_t>(k)]);
    Eigen::VectorXd rm;
    std::vector<trip> jt;
    memory_residual(_t, _m, x, real(1), rm, &jt);
    real trH = 0;
    for (const trip &q : jt)
      trH += q.value() * q.value();
    _scale = trH > real(0) ? trL / trH : real(1);
    _factored = false;
  }

  const stencil_table &table() const { return _t; }
  std::vector<vec3> readout(const std::vector<vec3> &x) const { return edge_readings(_t, x); }

  real energy(const std::vector<vec3> &x, const std::vector<vec3> &h) const {
    Eigen::VectorXd r;
    residual(x, h, r, nullptr);
    return real(0.5) * r.squaredNorm();
  }

  solver_report iterate(std::vector<vec3> &x, const std::vector<vec3> &h) {
    solver_report rep;
    Eigen::VectorXd r;
    matS J;
    linearize(x, h, r, J);
    rep.energy0 = real(0.5) * r.squaredNorm();
    rep.energy1 = rep.energy0;
    const Eigen::VectorXd g = J.transpose() * r;
    if (!(g.norm() > real(0)))
      return rep;
    if (!_factored)
      factor(J);
    for (int attempt = 0; attempt < 2; ++attempt) {
      const Eigen::VectorXd d = -_ldlt.solve(g);
      const real slope = g.dot(d);
      real a = 1;
      for (int k = 0; k < backtracks && slope < real(0) && d.allFinite(); ++k, a *= real(0.5)) {
        std::vector<vec3> xt = x;
        for (size_t v = 0; v < x.size(); ++v)
          xt[v] += a * d.segment<3>(3 * static_cast<index_t>(v));
        const real e1 = energy(xt, h);
        if (std::isfinite(e1) && e1 <= rep.energy0 + real(1e-4) * a * slope) {
          rep.step = a * d.cwiseAbs().maxCoeff();
          rep.energy1 = e1;
          rep.halvings = k;
          rep.accepted = true;
          x = std::move(xt);
          return rep;
        }
      }
      if (attempt == 0) {
        factor(J);
        rep.refactored = true;
      }
    }
    return rep;
  }

  /// Restart protocol: `restarts` times, converge on the frozen system
  /// (`max_its` iterates or a step below `tol`), then re-rest at the result.
  void solve(std::vector<vec3> &x, const std::vector<vec3> &h, int restarts = 5,
             int max_its = 30, real tol = real(1e-10)) {
    for (int k = 0; k < restarts; ++k) {
      reset_rest(x);
      for (int it = 0; it < max_its; ++it) {
        const solver_report rep = iterate(x, h);
        if (!rep.accepted || rep.step < tol)
          break;
      }
    }
  }

private:
  void residual(const std::vector<vec3> &x, const std::vector<vec3> &h, Eigen::VectorXd &r,
                std::vector<trip> *jt) const {
    Eigen::VectorXd rc, rm;
    edge_coupling(_t, x, _mass, &h, rc, jt, 0);
    memory_residual(_t, _m, x, w_memory * _scale, rm, jt, static_cast<index_t>(rc.size()));
    r.resize(rc.size() + rm.size());
    r << rc, rm;
  }

  void linearize(const std::vector<vec3> &x, const std::vector<vec3> &h, Eigen::VectorXd &r,
                 matS &J) const {
    std::vector<trip> jt;
    residual(x, h, r, &jt);
    J.resize(r.size(), 3 * _t.n_verts);
    J.setFromTriplets(jt.begin(), jt.end());
  }

  void factor(const matS &J) {
    matS H = J.transpose() * J;
    matS I(H.rows(), H.cols());
    I.setIdentity();
    H += eps * I;
    _ldlt.compute(H);
    _factored = _ldlt.info() == Eigen::Success;
  }

  stencil_table _t;
  memory_rest _m;
  std::vector<real> _mass;
  real _scale = 1;
  bool _factored = false;
  Eigen::SimplicialLDLT<matS> _ldlt;
};

} // namespace weingarten
} // namespace kusama
} // namespace gaudi

#endif
