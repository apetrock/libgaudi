#ifndef GAUDI_HEPWORTH_BLOCK_EXPONENT_COUPLING_HPP
#define GAUDI_HEPWORTH_BLOCK_EXPONENT_COUPLING_HPP

/// The exponent curvature solver (scripts/sympy/cr_poisson/MATH.md section 16;
/// prototype `exponent_solver.py`).
///
/// Each face carries an exponent `xi_f = [omega_f]x + P_f L_f P_f^T`: its
/// displacement gradient, a rotation `omega` (3) and an in-plane log-stretch
/// `L` (3, `(l11, l12, l22)` in the face frame `T, B`). One constant system per
/// restart at the frame-start positions:
///  - curvature rows `sqrt(W) ([J_theta C, J_ell E] xi - r)`, `r = s_g - s(X)`:
///    Stein's face read-out in the discrete fundamental forms (the generated
///    rebuild `cr_ff_rebuild` composed with the generated stencil Jacobian
///    `cr_curvature_stencil`); `C` takes `dtheta_e = a_e . (omega_- - omega_+)`,
///    `E` takes `dlog l_e = t^T L t` (the stencil face owning the edge);
///  - compatibility `xi_f d_e = xi_g d_e` on every interior edge, stiff: the
///    field is the gradient of a displacement;
///  - priors `kb |omega|^2` (removes the rigid null space) and, elastic,
///    `km |L_cur + dL|^2` against a persistent rest.
/// Each iterate: `xi` by one backsolve, `F_f = R(omega_f) U(L_f)` per face
/// (rotations leave the solve through the exponential), the advected frames
/// `sum_f F_f g_f / M_e`, and positions by the CR Dirichlet operator
/// `D^T M D` in increment form, with a line search on
/// `|r|_W^2 + km tau sum A_f / mean A |L_cur|^2`. The factor (`factor`) is
/// reused across steps in each face's body frame (`factor_cache`).

#include "gaudi/hepworth/block/block_constraint.hpp"
#include "gaudi/hepworth/block/shape_det_conserve.hpp"
#include "gaudi/kusama/stencil_table.hpp"
#include "gaudi/kusama/weingarten.hpp"

#include "cr_curvature_stencil_generated.hpp"
#include "cr_ff_rebuild_generated.hpp"

#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <utility>
#include <vector>

namespace gaudi {
namespace hepworth {
namespace block {
namespace exponent {

using matS = Eigen::SparseMatrix<real>;
using trip = Eigen::Triplet<real>;
using mat2 = Eigen::Matrix<real, 2, 2>;
using mat32 = Eigen::Matrix<real, 3, 2>;

/// The curvature rows: `tensor` the face read-out `s_f` (3 per face),
/// `principal` its eigenvalues (Eigensatz 2010 eq. 3.12, 2 per face).
enum class residual { tensor, principal };

/// Material and solver knobs. `kb` the bending prior (tiny: the rigid null
/// space), `km` the growth stiffness (0 free; with `elastic`, against the
/// rest), `w_compat` the compatibility weight (relative), `eps` the gather's
/// regularization, `line_search` the step halvings tried, `max_step` a cap on
/// `|xi|` per face.
struct config {
  real kb = 1e-6;
  real km = 0.0;
  real w_compat = 1e8;
  real eps = 1e-10;
  int line_search = 12;
  real max_step = 0.5;
  /// Exponent steps per solver step (`exponent_coupling`): later projections
  /// only hold the reached frames. The solve converges in a few steps.
  int steps_per_frame = 4;
  /// Mass on the exponents, `mass tau |xi|^2` in `H` (relative to the
  /// curvature rows): shorter, safer steps through a reused factor.
  real mass = 0.0;
  /// Solver steps a factor is reused for (`factor_cache`); 1 builds one per
  /// step, so every tangent-dependent row is reseeded from the live mesh. A
  /// reused factor keeps its rows and `H` together (frozen at its build,
  /// turned per face); a stalled step refactors at once.
  int refresh_every = 1;
  /// The curvature rows (`residual`).
  residual kind = residual::tensor;
  /// Eigensatz's metric terms (eq. 3.11) on the stretch, alongside `km`:
  /// `w_conformal` its trace-free part (angles), `w_area` its trace (log
  /// area). Measured like `km` (against the rest, or damping per step).
  real w_conformal = 0.0;
  real w_area = 0.0;
  /// `principal`: within this fraction of `max |kappa|` of an umbilic both
  /// rows take the mean curvature's derivative (Eigensatz A.6).
  real umbilic = 1e-3;
};

/// The stretch stiffness on `(l11, l12, l22)`: `km` on all of it, the
/// conformal rows on its trace-free part, the areal row on its trace.
inline Eigen::Matrix<real, 3, 3> stretch_stiffness(const config &cfg) {
  const vec3 c1(real(0.5), 0, real(-0.5)), c2(0, 1, 0), ct(1, 0, 1);
  return cfg.km * Eigen::Matrix<real, 3, 3>::Identity() +
         cfg.w_conformal * (c1 * c1.transpose() + c2 * c2.transpose()) + cfg.w_area * ct * ct.transpose();
}

inline bool has_stretch(const config &cfg) {
  return cfg.km > real(0) || cfg.w_conformal > real(0) || cfg.w_area > real(0);
}

/// The principal curvatures of `s = (k11, k12, k22)`, larger first, as
/// `(k1, k2, 0)`.
inline vec3 principal_values(const vec3 &s) {
  const real m = real(0.5) * (s[0] + s[2]);
  const real d = std::sqrt(real(0.25) * (s[0] - s[2]) * (s[0] - s[2]) + s[1] * s[1]);
  return vec3(m + d, m - d, 0);
}

/// `d (k1, k2, 0) / d s` (Eigensatz A.4); both `dH` (A.6) within `umbilic`.
inline Eigen::Matrix<real, 3, 3> principal_rows(const vec3 &s, real umbilic) {
  const real m = real(0.5) * (s[0] + s[2]);
  const real h = real(0.5) * (s[0] - s[2]);
  const real d = std::sqrt(h * h + s[1] * s[1]);
  Eigen::Matrix<real, 3, 3> K = Eigen::Matrix<real, 3, 3>::Zero();
  K.row(0) << real(0.5), 0, real(0.5);
  K.row(1) = K.row(0);
  if (d > umbilic * (std::abs(m) + d)) {
    const Eigen::Matrix<real, 1, 3> dd(real(0.5) * h / d, s[1] / d, real(-0.5) * h / d);
    K.row(0) += dd;
    K.row(1) -= dd;
  }
  return K;
}

/// A face's residual against its goal: `goal - s` (tensor) or
/// `goal - principal_values(s)`.
inline vec3 face_residual(residual kind, const vec3 &goal, const vec3 &s) {
  return kind == residual::principal ? vec3(goal - principal_values(s)) : vec3(goal - s);
}

/// Interior edges with their + face (the face whose local edge runs from the
/// lower vertex id to the higher), and per face the global columns of its 9
/// lengths (own 3, then each wing's edges from tail and head to the opposite
/// vertex) and 3 dihedrals (interior index, -1 on a boundary).
struct topology {
  struct inner_edge {
    int e = -1, fp = -1, ip = -1, fm = -1, im = -1;
  };
  std::vector<inner_edge> inner;
  std::vector<std::array<int, 2>> ends;
  std::vector<std::array<int, 9>> lcol;
  std::vector<std::array<int, 3>> tcol;
  std::vector<std::vector<std::pair<int, int>>> faces_of_edge;  // (face, local index)
};

inline topology build_topology(const kusama::stencil_table &t) {
  topology tp;
  const int nF = t.n_faces, nE = t.n_edges;
  tp.ends.assign(static_cast<size_t>(nE), {-1, -1});
  tp.faces_of_edge.assign(static_cast<size_t>(nE), {});
  std::map<std::pair<int, int>, int> index;
  for (int f = 0; f < nF; ++f)
    for (int i = 0; i < 3; ++i) {
      const int e = t.face_edges[static_cast<size_t>(f)][static_cast<size_t>(i)];
      const int a = t.ids[static_cast<size_t>(f)][static_cast<size_t>(i)];
      const int b = t.ids[static_cast<size_t>(f)][static_cast<size_t>((i + 1) % 3)];
      tp.ends[static_cast<size_t>(e)] = {std::min(a, b), std::max(a, b)};
      tp.faces_of_edge[static_cast<size_t>(e)].push_back({f, i});
      index[{std::min(a, b), std::max(a, b)}] = e;
    }
  std::vector<int> inner_of(static_cast<size_t>(nE), -1);
  for (int e = 0; e < nE; ++e) {
    const auto &fs = tp.faces_of_edge[static_cast<size_t>(e)];
    if (fs.size() != 2)
      continue;
    topology::inner_edge ie;
    ie.e = e;
    for (const auto &[f, i] : fs) {
      const int a = t.ids[static_cast<size_t>(f)][static_cast<size_t>(i)];
      const int b = t.ids[static_cast<size_t>(f)][static_cast<size_t>((i + 1) % 3)];
      if (a < b) {
        ie.fp = f;
        ie.ip = i;
      } else {
        ie.fm = f;
        ie.im = i;
      }
    }
    if (ie.fp < 0 || ie.fm < 0)
      continue;
    inner_of[static_cast<size_t>(e)] = static_cast<int>(tp.inner.size());
    tp.inner.push_back(ie);
  }
  tp.lcol.assign(static_cast<size_t>(nF), {});
  tp.tcol.assign(static_cast<size_t>(nF), {-1, -1, -1});
  for (int f = 0; f < nF; ++f) {
    auto &lc = tp.lcol[static_cast<size_t>(f)];
    lc.fill(-1);
    for (int i = 0; i < 3; ++i) {
      const int e = t.face_edges[static_cast<size_t>(f)][static_cast<size_t>(i)];
      lc[static_cast<size_t>(i)] = e;
      if (t.boundary(f, i))
        continue;
      const int o = t.ids[static_cast<size_t>(f)][static_cast<size_t>(3 + i)];
      const int a = t.ids[static_cast<size_t>(f)][static_cast<size_t>(i)];
      const int b = t.ids[static_cast<size_t>(f)][static_cast<size_t>((i + 1) % 3)];
      lc[static_cast<size_t>(3 + 2 * i)] = index.at({std::min(a, o), std::max(a, o)});
      lc[static_cast<size_t>(4 + 2 * i)] = index.at({std::min(b, o), std::max(b, o)});
      tp.tcol[static_cast<size_t>(f)][static_cast<size_t>(i)] = inner_of[static_cast<size_t>(e)];
    }
  }
  return tp;
}

/// Face frame `(T, B, N)` (T along the first local edge) and area.
inline void face_frame(const vec3 &a, const vec3 &b, const vec3 &c, vec3 &T, vec3 &B, vec3 &N,
                       real &A) {
  const vec3 n = (b - a).cross(c - a);
  A = real(0.5) * n.norm();
  N = n.normalized();
  T = (b - a).normalized();
  B = N.cross(T);
}

/// The fundamental forms of a face's live stencil `P` (boundary slots per
/// `h`): `l` its 9 lengths, `th` its 3 dihedrals (sync convention).
inline void face_forms(const std::array<vec3, 6> &P, const vec3 &h, std::array<real, 9> &l,
                       vec3 &th) {
  l.fill(real(1));
  th.setZero();
  const vec3 Nf = (P[1] - P[0]).cross(P[2] - P[0]).normalized();
  for (int i = 0; i < 3; ++i)
    l[static_cast<size_t>(i)] = (P[static_cast<size_t>((i + 1) % 3)] - P[static_cast<size_t>(i)]).norm();
  for (int i = 0; i < 3; ++i) {
    if (!(h[i] > real(0)))
      continue;
    const vec3 &Pt = P[static_cast<size_t>(i)], &Ph = P[static_cast<size_t>((i + 1) % 3)],
               &Po = P[static_cast<size_t>(3 + i)];
    l[static_cast<size_t>(3 + 2 * i)] = (Po - Pt).norm();
    l[static_cast<size_t>(4 + 2 * i)] = (Po - Ph).norm();
    const vec3 ax = (Ph - Pt).normalized();
    const vec3 Ng = (Po - Ph).cross(Pt - Ph).normalized();
    const real phi = std::atan2(Nf.cross(Ng).dot(ax), Nf.dot(Ng));
    th[i] = std::remainder(phi - real(M_PI), real(2) * real(M_PI));
  }
}

/// `d s / d (log l, th)` (3 x 12) of a face from its live stencil, and its
/// read-out `s` (Stein, App. A.2), through the rebuild from the forms.
inline Eigen::Matrix<real, 3, 12> face_forms_jacobian(const std::array<vec3, 6> &Plive,
                                                      const vec3 &h, vec3 &s) {
  std::array<real, 9> l;
  vec3 th;
  face_forms(Plive, h, l, th);
  std::array<vec3, 6> P;
  Eigen::Matrix<real, 18, 12> dP;
  cr_generated::cr_ff_rebuild(l, th, P, dP);
  for (int i = 0; i < 3; ++i)
    if (!(h[i] > real(0)))
      P[static_cast<size_t>(3 + i)] = P[static_cast<size_t>(i)] + P[static_cast<size_t>((i + 1) % 3)] -
                                      P[static_cast<size_t>((i + 2) % 3)];
  Eigen::Matrix<real, 3, 18> Jp;
  cr_generated::cr_curvature_stencil(P, h, s, Jp);
  for (int i = 0; i < 3; ++i)
    if (!(h[i] > real(0)))
      Jp.block<3, 3>(0, 9 + 3 * i).setZero();
  Eigen::Matrix<real, 3, 12> J = Jp * dP;
  for (int j = 0; j < 9; ++j)
    J.col(j) *= l[static_cast<size_t>(j)];
  return J;
}

/// The live read-out of every face: `(k11, k12, k22)` in its frame.
inline std::vector<vec3> face_readouts(const kusama::stencil_table &t, const std::vector<vec3> &x) {
  std::vector<vec3> out(static_cast<size_t>(t.n_faces));
  Eigen::Matrix<real, 3, 18> J;
  for (int f = 0; f < t.n_faces; ++f) {
    vec3 h;
    const std::array<vec3, 6> P = kusama::stencil_positions(t, x, f, h);
    cr_generated::cr_curvature_stencil(P, h, out[static_cast<size_t>(f)], J);
  }
  return out;
}

/// Per face, the log of the left stretch rest -> current in the current face
/// frame, `(l11, l12, l22)`.
inline std::vector<vec3> current_log_stretch(const kusama::stencil_table &t, const std::vector<vec3> &x,
                                             const std::vector<std::array<vec3, 3>> &rest) {
  std::vector<vec3> out(static_cast<size_t>(t.n_faces), vec3::Zero());
  for (int f = 0; f < t.n_faces && f < static_cast<int>(rest.size()); ++f) {
    const auto &id = t.ids[static_cast<size_t>(f)];
    vec3 T, B, N, Tr, Br, Nr;
    real A, Ar;
    const vec3 &a = x[static_cast<size_t>(id[0])], &b = x[static_cast<size_t>(id[1])],
               &c = x[static_cast<size_t>(id[2])];
    const auto &R = rest[static_cast<size_t>(f)];
    face_frame(a, b, c, T, B, N, A);
    face_frame(R[0], R[1], R[2], Tr, Br, Nr, Ar);
    mat2 Am, Bm;
    Am << Tr.dot(R[1] - R[0]), Tr.dot(R[2] - R[0]), Br.dot(R[1] - R[0]), Br.dot(R[2] - R[0]);
    Bm << T.dot(b - a), T.dot(c - a), B.dot(b - a), B.dot(c - a);
    const mat2 Fd = Bm * Am.inverse();
    Eigen::SelfAdjointEigenSolver<mat2> es(Fd * Fd.transpose());
    const Eigen::Matrix<real, 2, 1> w = es.eigenvalues().cwiseMax(real(1e-30));
    const mat2 Lg = es.eigenvectors() * (real(0.5) * w.array().log()).matrix().asDiagonal() *
                    es.eigenvectors().transpose();
    out[static_cast<size_t>(f)] = vec3(Lg(0, 0), Lg(0, 1), Lg(1, 1));
  }
  return out;
}

inline Eigen::Matrix<real, 3, 3> skew(const vec3 &w) {
  Eigen::Matrix<real, 3, 3> K;
  K << 0, -w[2], w[1], w[2], 0, -w[0], -w[1], w[0], 0;
  return K;
}

/// `exp(L)` of a symmetric 2x2 `L = (l11, l12, l22)`, closed form.
inline mat2 sym_exp2(const vec3 &L) {
  const real m = real(0.5) * (L[0] + L[2]);
  const real d = std::sqrt(real(0.25) * (L[0] - L[2]) * (L[0] - L[2]) + L[1] * L[1]);
  mat2 D;
  D << real(0.5) * (L[0] - L[2]), L[1], L[1], -real(0.5) * (L[0] - L[2]);
  const real c = std::cosh(d), sh = d > real(1e-12) ? std::sinh(d) / d : real(1);
  return std::exp(m) * (c * mat2::Identity() + sh * D);
}

/// The per-face deformation `R(omega) U(L)`: a unit-quaternion rotation and
/// the stretch `P exp(L) P^T + N N^T` in the face plane.
inline Eigen::Matrix<real, 3, 3> face_deformation(const vec3 &omega, const vec3 &L, const mat32 &P,
                                                  const vec3 &N) {
  const real a = omega.norm();
  const Eigen::Quaternion<real> q =
      a > real(1e-300) ? Eigen::Quaternion<real>(Eigen::AngleAxis<real>(a, omega / a))
                       : Eigen::Quaternion<real>::Identity();
  const Eigen::Matrix<real, 3, 3> U = P * sym_exp2(L) * P.transpose() + N * N.transpose();
  return q.toRotationMatrix() * U;
}

/// `sum_f F_f (face f's part of each edge's CR frame at `x`) / M_e`: rows
/// `[g_par (nE); g_perp (nE)]`.
inline std::vector<vec3> advected_frames(const kusama::stencil_table &t, const std::vector<vec3> &x,
                                         const std::vector<Eigen::Matrix<real, 3, 3>> &F,
                                         const std::vector<real> &mass) {
  const int nE = t.n_edges;
  std::vector<vec3> G(static_cast<size_t>(2 * nE), vec3::Zero());
  for (int f = 0; f < t.n_faces; ++f) {
    const auto &id = t.ids[static_cast<size_t>(f)];
    const std::array<vec3, 3> Y = {x[static_cast<size_t>(id[0])], x[static_cast<size_t>(id[1])],
                                   x[static_cast<size_t>(id[2])]};
    const real A = real(0.5) * (Y[1] - Y[0]).cross(Y[2] - Y[0]).norm();
    std::array<real, 3> L;
    for (int i = 0; i < 3; ++i)
      L[static_cast<size_t>(i)] = (Y[static_cast<size_t>((i + 1) % 3)] - Y[static_cast<size_t>(i)]).norm();
    for (int i = 0; i < 3; ++i) {
      const int e = t.face_edges[static_cast<size_t>(f)][static_cast<size_t>(i)];
      const int tl = i, hd = (i + 1) % 3, op = (i + 2) % 3;
      const real si = id[static_cast<size_t>(tl)] < id[static_cast<size_t>(hd)] ? real(1) : real(-1);
      const real li = L[static_cast<size_t>(i)];
      const real ejk = L[static_cast<size_t>((i + 1) % 3)] * L[static_cast<size_t>((i + 1) % 3)];
      const real eki = L[static_cast<size_t>((i + 2) % 3)] * L[static_cast<size_t>((i + 2) % 3)];
      const vec3 gp = si * A / (real(3) * li) * (Y[static_cast<size_t>(hd)] - Y[static_cast<size_t>(tl)]);
      const real k = -si / (real(12) * li);
      const vec3 gq = k * (li * li - (eki - ejk)) * Y[static_cast<size_t>(tl)] +
                      k * (li * li + (eki - ejk)) * Y[static_cast<size_t>(hd)] +
                      si * li / real(6) * Y[static_cast<size_t>(op)];
      G[static_cast<size_t>(e)] += F[static_cast<size_t>(f)] * gp;
      G[static_cast<size_t>(nE + e)] += F[static_cast<size_t>(f)] * gq;
    }
  }
  for (int e = 0; e < nE; ++e) {
    const real m = mass[static_cast<size_t>(e)] > real(0) ? mass[static_cast<size_t>(e)] : real(1);
    G[static_cast<size_t>(e)] /= m;
    G[static_cast<size_t>(nE + e)] /= m;
  }
  return G;
}

/// The factored half of a step: curvature rows `Jc`, compatibility and
/// priors at the positions it is built at, `H` factored once (LDLT). Steps
/// are solved in each face's body frame (`basis`, `[T B N]` at the build),
/// where `H` is a metric core (MATH.md section 16), so one factor serves
/// later positions of the same mesh: `system` turns the rotations back.
class factor {
public:
  factor(const kusama::stencil_table &t, const topology &tp, const std::vector<vec3> &x, const config &cfg)
      : _cfg(cfg), _ids(t.ids) {
    const int nF = t.n_faces, nI = static_cast<int>(tp.inner.size());
    _nF = nF;
    // frames, areas
    _P.resize(static_cast<size_t>(nF));
    _N.resize(static_cast<size_t>(nF));
    _Q.resize(static_cast<size_t>(nF));
    _W.resize(3 * nF);
    for (int f = 0; f < nF; ++f) {
      const auto &id = t.ids[static_cast<size_t>(f)];
      vec3 T, B, N;
      real A;
      face_frame(x[static_cast<size_t>(id[0])], x[static_cast<size_t>(id[1])], x[static_cast<size_t>(id[2])], T, B,
                 N, A);
      _P[static_cast<size_t>(f)].col(0) = T;
      _P[static_cast<size_t>(f)].col(1) = B;
      _N[static_cast<size_t>(f)] = N;
      _Q[static_cast<size_t>(f)] << T, B, N;
      _W.segment<3>(3 * f).setConstant(A);
    }
    std::vector<vec3> axis(static_cast<size_t>(nI));
    for (int k = 0; k < nI; ++k) {
      const auto &ie = tp.inner[static_cast<size_t>(k)];
      axis[static_cast<size_t>(k)] = (x[static_cast<size_t>(tp.ends[static_cast<size_t>(ie.e)][1])] -
                                      x[static_cast<size_t>(tp.ends[static_cast<size_t>(ie.e)][0])])
                                         .normalized();
    }
    // curvature rows Jc (3 nF x 6 nF). E: dlog l_e = t^T L_g t with g the
    // stencil face owning e (f for its own edges, the wing face for the wing
    // edges); the compatibility rows make the faces of an edge agree, and the
    // single owner halves the width of H.
    std::vector<trip> jt;
    for (int f = 0; f < nF; ++f) {
      vec3 h, s;
      const std::array<vec3, 6> Pl = kusama::stencil_positions(t, x, f, h);
      const Eigen::Matrix<real, 3, 12> Jf = face_forms_jacobian(Pl, h, s);
      for (int j = 0; j < 9; ++j) {
        const int e = tp.lcol[static_cast<size_t>(f)][static_cast<size_t>(j)];
        if (e < 0)
          continue;
        int g = f;
        if (j >= 3) {
          const auto &ie =
              tp.inner[static_cast<size_t>(tp.tcol[static_cast<size_t>(f)][static_cast<size_t>((j - 3) / 2)])];
          g = ie.fp == f ? ie.fm : ie.fp;
        }
        const vec3 d = (x[static_cast<size_t>(tp.ends[static_cast<size_t>(e)][1])] -
                        x[static_cast<size_t>(tp.ends[static_cast<size_t>(e)][0])])
                           .normalized();
        const Eigen::Matrix<real, 2, 1> tv = _P[static_cast<size_t>(g)].transpose() * d;
        const real cf[3] = {tv[0] * tv[0], real(2) * tv[0] * tv[1], tv[1] * tv[1]};
        for (int b = 0; b < 3; ++b)
          for (int r = 0; r < 3; ++r)
            jt.emplace_back(3 * f + r, 3 * nF + 3 * g + b, Jf(r, j) * cf[b]);
      }
      for (int i = 0; i < 3; ++i) {
        const int k = tp.tcol[static_cast<size_t>(f)][static_cast<size_t>(i)];
        if (k < 0)
          continue;
        const auto &ie = tp.inner[static_cast<size_t>(k)];
        const vec3 &a = axis[static_cast<size_t>(k)];
        for (int r = 0; r < 3; ++r)
          for (int c = 0; c < 3; ++c) {
            jt.emplace_back(3 * f + r, 3 * ie.fm + c, Jf(r, 9 + i) * a[c]);
            jt.emplace_back(3 * f + r, 3 * ie.fp + c, -Jf(r, 9 + i) * a[c]);
          }
      }
    }
    _Jc.resize(3 * nF, 6 * nF);
    _Jc.setFromTriplets(jt.begin(), jt.end());
    if (cfg.kind == residual::principal) {
      const std::vector<vec3> s0 = face_readouts(t, x);
      std::vector<trip> kt;
      for (int f = 0; f < nF; ++f) {
        const Eigen::Matrix<real, 3, 3> K = principal_rows(s0[static_cast<size_t>(f)], cfg.umbilic);
        for (int a = 0; a < 3; ++a)
          for (int b = 0; b < 3; ++b)
            if (K(a, b) != real(0))
              kt.emplace_back(3 * f + a, 3 * f + b, K(a, b));
      }
      matS Km(3 * nF, 3 * nF);
      Km.setFromTriplets(kt.begin(), kt.end());
      _Jc = matS(Km * _Jc);
    }
    // compatibility (3 nI x 6 nF)
    std::vector<trip> ct;
    for (int k = 0; k < nI; ++k) {
      const auto &ie = tp.inner[static_cast<size_t>(k)];
      const vec3 d = x[static_cast<size_t>(tp.ends[static_cast<size_t>(ie.e)][1])] -
                     x[static_cast<size_t>(tp.ends[static_cast<size_t>(ie.e)][0])];
      const Eigen::Matrix<real, 3, 3> Dx = -skew(d);
      for (const auto &[f, sgn] : {std::pair<int, real>{ie.fp, real(1)}, std::pair<int, real>{ie.fm, real(-1)}}) {
        for (int r = 0; r < 3; ++r)
          for (int c = 0; c < 3; ++c)
            if (Dx(r, c) != real(0))
              ct.emplace_back(3 * k + r, 3 * f + c, sgn * Dx(r, c));
        const mat32 &P = _P[static_cast<size_t>(f)];
        const std::array<mat2, 3> Eb = {(mat2() << 1, 0, 0, 0).finished(), (mat2() << 0, 1, 1, 0).finished(),
                                        (mat2() << 0, 0, 0, 1).finished()};
        for (int b = 0; b < 3; ++b) {
          const vec3 v = P * Eb[static_cast<size_t>(b)] * P.transpose() * d;
          for (int r = 0; r < 3; ++r)
            ct.emplace_back(3 * k + r, 3 * nF + 3 * f + b, sgn * v[r]);
        }
      }
    }
    matS Comp(3 * nI, 6 * nF);
    Comp.setFromTriplets(ct.begin(), ct.end());
    const matS JWJ = matS(_Jc.transpose() * _W.asDiagonal() * _Jc);
    real dsum = 0;
    for (int i = 0; i < 3 * nF; ++i)
      dsum += JWJ.coeff(i, i);
    _tt = dsum / std::max(3 * nF, 1);
    const matS CC = matS(Comp.transpose() * Comp);
    real csum = 0;
    for (int i = 0; i < CC.rows(); ++i)
      csum += CC.coeff(i, i);
    const real sc = csum / std::max<real>(CC.rows(), 1);
    // The stretch prior weighs faces by `A_f / mean A`, as the merit does.
    const real Am = _W.sum() / std::max<real>(_W.size(), 1);
    _Ms = stretch_stiffness(cfg);
    std::vector<trip> pt;
    for (int i = 0; i < 3 * nF; ++i) {
      pt.emplace_back(i, i, (cfg.kb + cfg.mass) * _tt);
      pt.emplace_back(3 * nF + i, 3 * nF + i, (cfg.kb + cfg.mass) * _tt);
    }
    for (int f = 0; f < nF; ++f)
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
          if (_Ms(a, b) != real(0))
            pt.emplace_back(3 * nF + 3 * f + a, 3 * nF + 3 * f + b, _Ms(a, b) * _W[3 * f] / Am * _tt);
    matS Pr(6 * nF, 6 * nF);
    Pr.setFromTriplets(pt.begin(), pt.end());
    const matS H = JWJ + (cfg.w_compat * _tt / std::max(sc, real(1e-30))) * CC + Pr;
    _H.compute(H);
  }

  bool ok() const { return _H.info() == Eigen::Success; }
  /// Same faces (and so the same stencils) as `t`.
  bool matches(const kusama::stencil_table &t) const { return t.ids == _ids; }
  real tau() const { return _tt; }
  const vecX &weights() const { return _W; }
  /// `[T B N]` of each face at the build.
  const std::vector<Eigen::Matrix<real, 3, 3>> &basis() const { return _Q; }

  /// The body-frame step for the residual `r` (3 nF); elastic, `L_cur` the
  /// live log-stretch against the rest (the prior's gradient).
  vecX solve(const vecX &r, const std::vector<vec3> *L_cur = nullptr) const {
    vecX b = _Jc.transpose() * (_W.cwiseProduct(r));
    if (L_cur && has_stretch(_cfg)) {
      const real Am = _W.sum() / std::max<real>(_W.size(), 1);
      for (int f = 0; f < _nF && f < static_cast<int>(L_cur->size()); ++f)
        b.segment<3>(3 * _nF + 3 * f) -= _tt * _W[3 * f] / Am * (_Ms * (*L_cur)[static_cast<size_t>(f)]);
    }
    return _H.solve(b);
  }

private:
  config _cfg;
  decltype(kusama::stencil_table::ids) _ids;
  int _nF = 0;
  std::vector<mat32> _P;
  std::vector<vec3> _N;
  std::vector<Eigen::Matrix<real, 3, 3>> _Q;
  vecX _W;
  matS _Jc;
  Eigen::Matrix<real, 3, 3> _Ms = Eigen::Matrix<real, 3, 3>::Zero();
  real _tt = 1;
  Eigen::SimplicialLDLT<matS> _H;
};

/// One step's system at `x`: the CR derivative and its Dirichlet gather
/// `D^T M D` (cheap, always at `x`), the face frames at `x`, and a `factor`
/// (built here unless one of the same mesh is passed in).
class system {
public:
  system(const kusama::stencil_table &t, const topology &tp, const std::vector<vec3> &x, const config &cfg,
         std::shared_ptr<const factor> F = nullptr)
      : _t(&t), _tp(&tp), _cfg(cfg) {
    const int nF = t.n_faces, nE = t.n_edges;
    _nF = nF;
    _D = kusama::weingarten::build_cr_derivative(t, x);
    _mass.resize(static_cast<size_t>(nE));
    for (int e = 0; e < nE; ++e)
      _mass[static_cast<size_t>(e)] = _D[static_cast<size_t>(e)].mass;
    // gather: D^T M D + eps I
    {
      std::vector<trip> tr;
      for (int e = 0; e < nE; ++e) {
        const auto &r = _D[static_cast<size_t>(e)];
        for (int a = 0; a < 4; ++a)
          for (int b = 0; b < 4; ++b) {
            if (r.v[static_cast<size_t>(a)] < 0 || r.v[static_cast<size_t>(b)] < 0)
              continue;
            const real v = r.mass * (r.par[static_cast<size_t>(a)] * r.par[static_cast<size_t>(b)] +
                                     r.perp[static_cast<size_t>(a)] * r.perp[static_cast<size_t>(b)]);
            tr.emplace_back(r.v[static_cast<size_t>(a)], r.v[static_cast<size_t>(b)], v);
          }
      }
      for (int v = 0; v < t.n_verts; ++v)
        tr.emplace_back(v, v, cfg.eps);
      matS L(t.n_verts, t.n_verts);
      L.setFromTriplets(tr.begin(), tr.end());
      _gather.compute(L);
    }
    set_frames(x);
    use(F && F->matches(t) ? std::move(F) : std::make_shared<const factor>(t, tp, x, cfg), false);
  }

  /// A fresh factor and face frames at `x` (a restart after a stalled step);
  /// the gather stays at the step start, where the rows were filled.
  void refactor(const std::vector<vec3> &x) {
    set_frames(x);
    use(std::make_shared<const factor>(*_t, *_tp, x, _cfg), true);
  }

  bool ok() const { return _F->ok() && _gather.info() == Eigen::Success; }
  /// The factor was built at this system's positions (or by `refactor`).
  bool fresh() const { return _fresh; }
  const std::shared_ptr<const factor> &factor_ptr() const { return _F; }
  real tau() const { return _F->tau(); }
  const std::vector<real> &mass() const { return _mass; }
  const std::vector<kusama::weingarten::cr_edge_rows> &derivative() const { return _D; }
  const vecX &weights() const { return _F->weights(); }

  /// `xi` for the residual `r` (3 nF): the factor's body-frame step with
  /// `omega` turned to the world, capped per face.
  vecX solve(const vecX &r, const std::vector<vec3> *L_cur = nullptr) const {
    vecX xi = _F->solve(r, L_cur);
    for (int f = 0; f < _nF; ++f)
      xi.segment<3>(3 * f) = _R[static_cast<size_t>(f)] * vec3(xi.segment<3>(3 * f));
    for (int f = 0; f < 2 * _nF; ++f) {
      const real n = xi.segment<3>(3 * f).norm();
      if (n > _cfg.max_step)
        xi.segment<3>(3 * f) *= _cfg.max_step / n;
    }
    return xi;
  }

  /// The advected frames of `x` after the step `a xi`.
  std::vector<vec3> frames(const std::vector<vec3> &x, const vecX &xi, real a) const {
    std::vector<Eigen::Matrix<real, 3, 3>> F(static_cast<size_t>(_nF));
    for (int f = 0; f < _nF; ++f)
      F[static_cast<size_t>(f)] = face_deformation(a * vec3(xi.segment<3>(3 * f)),
                                                   a * vec3(xi.segment<3>(3 * _nF + 3 * f)),
                                                   _P[static_cast<size_t>(f)], _N[static_cast<size_t>(f)]);
    return advected_frames(*_t, x, F, _mass);
  }

  /// The CR frames of `x` itself (no step).
  std::vector<vec3> frames(const std::vector<vec3> &x) const {
    return advected_frames(*_t, x, std::vector<Eigen::Matrix<real, 3, 3>>(static_cast<size_t>(_nF),
                                                                          Eigen::Matrix<real, 3, 3>::Identity()),
                           _mass);
  }

  /// `D x` with this system's derivative (the rows `exponent_coupling` fills).
  std::vector<vec3> derivative_of(const std::vector<vec3> &x) const {
    const int nE = _t->n_edges;
    std::vector<vec3> out(static_cast<size_t>(2 * nE), vec3::Zero());
    for (int e = 0; e < nE; ++e) {
      const auto &r = _D[static_cast<size_t>(e)];
      for (int k = 0; k < 4; ++k) {
        const int v = r.v[static_cast<size_t>(k)];
        if (v < 0)
          continue;
        out[static_cast<size_t>(e)] += r.par[static_cast<size_t>(k)] * x[static_cast<size_t>(v)];
        out[static_cast<size_t>(nE + e)] += r.perp[static_cast<size_t>(k)] * x[static_cast<size_t>(v)];
      }
    }
    return out;
  }

  /// Positions whose CR frames fit `G` (the Dirichlet gather), in increment
  /// form `x + (D^T M D)^-1 D^T M (G - G(x))`: the operator is this system's,
  /// `x` may have moved since, and `G = G(x)` returns `x` exactly.
  std::vector<vec3> gather(const std::vector<vec3> &G, const std::vector<vec3> &x) const {
    const int n = _t->n_verts, nE = _t->n_edges;
    const std::vector<vec3> Gx = frames(x);
    std::vector<vec3> out(static_cast<size_t>(n));
    for (int c = 0; c < 3; ++c) {
      vecX b = vecX::Zero(n);
      for (int e = 0; e < nE; ++e) {
        const auto &r = _D[static_cast<size_t>(e)];
        const real dp = G[static_cast<size_t>(e)][c] - Gx[static_cast<size_t>(e)][c];
        const real dq = G[static_cast<size_t>(nE + e)][c] - Gx[static_cast<size_t>(nE + e)][c];
        for (int k = 0; k < 4; ++k) {
          const int v = r.v[static_cast<size_t>(k)];
          if (v < 0)
            continue;
          b[v] += r.mass * (r.par[static_cast<size_t>(k)] * dp + r.perp[static_cast<size_t>(k)] * dq);
        }
      }
      const vecX s = _gather.solve(b);
      for (int v = 0; v < n; ++v)
        out[static_cast<size_t>(v)][c] = x[static_cast<size_t>(v)][c] + s[v];
    }
    return out;
  }

private:
  void set_frames(const std::vector<vec3> &x) {
    _P.resize(static_cast<size_t>(_nF));
    _N.resize(static_cast<size_t>(_nF));
    _Q.resize(static_cast<size_t>(_nF));
    for (int f = 0; f < _nF; ++f) {
      const auto &id = _t->ids[static_cast<size_t>(f)];
      vec3 T, B, N;
      real A;
      face_frame(x[static_cast<size_t>(id[0])], x[static_cast<size_t>(id[1])], x[static_cast<size_t>(id[2])], T, B,
                 N, A);
      _P[static_cast<size_t>(f)].col(0) = T;
      _P[static_cast<size_t>(f)].col(1) = B;
      _N[static_cast<size_t>(f)] = N;
      _Q[static_cast<size_t>(f)] << T, B, N;
    }
  }

  void use(std::shared_ptr<const factor> F, bool fresh) {
    _F = std::move(F);
    _fresh = fresh || _F->basis() == _Q;
    _R.resize(static_cast<size_t>(_nF));
    for (int f = 0; f < _nF; ++f)
      _R[static_cast<size_t>(f)] = _Q[static_cast<size_t>(f)] * _F->basis()[static_cast<size_t>(f)].transpose();
  }

  const kusama::stencil_table *_t = nullptr;
  const topology *_tp = nullptr;
  config _cfg;
  int _nF = 0;
  std::vector<kusama::weingarten::cr_edge_rows> _D;
  std::vector<real> _mass;
  std::vector<mat32> _P;
  std::vector<vec3> _N;
  std::vector<Eigen::Matrix<real, 3, 3>> _Q, _R;
  std::shared_ptr<const factor> _F;
  bool _fresh = true;
  Eigen::SimplicialLDLT<matS> _gather;
};

/// The merit `|r|_W^2 (+ tau sum A_f / mean A L_cur^T M L_cur)` and the
/// residual (`cfg->kind`; tensor without `cfg`), `M = stretch_stiffness`
/// against `rest`.
inline real merit(const kusama::stencil_table &t, const std::vector<vec3> &x, const std::vector<vec3> &goal,
                  const vecX &W, vecX &r, const config *cfg = nullptr, real tau = 0,
                  const std::vector<std::array<vec3, 3>> *rest = nullptr) {
  const residual kind = cfg ? cfg->kind : residual::tensor;
  const std::vector<vec3> s = face_readouts(t, x);
  r.resize(3 * t.n_faces);
  real e = 0;
  for (int f = 0; f < t.n_faces; ++f) {
    r.segment<3>(3 * f) = face_residual(kind, goal[static_cast<size_t>(f)], s[static_cast<size_t>(f)]);
    e += W[3 * f] * r.segment<3>(3 * f).squaredNorm();
  }
  if (cfg && rest && has_stretch(*cfg)) {
    const Eigen::Matrix<real, 3, 3> M = stretch_stiffness(*cfg);
    const std::vector<vec3> L = current_log_stretch(t, x, *rest);
    const real Am = W.sum() / std::max<real>(W.size(), 1);
    for (int f = 0; f < t.n_faces; ++f)
      e += tau * W[3 * f] / Am * L[static_cast<size_t>(f)].dot(M * L[static_cast<size_t>(f)]);
  }
  return std::sqrt(e);
}

/// Line search along `xi` from `x`: the largest `a` in 1, 1/2, ... whose
/// gathered positions lower the merit below `m0`. Returns `a` (0 when none),
/// the accepted positions in `x_out` and their merit in `m_out`.
inline real line_search(const kusama::stencil_table &t, const system &sys, const std::vector<vec3> &x,
                        const vecX &xi, const std::vector<vec3> &goal, real m0, const config &cfg,
                        std::vector<vec3> &x_out, real &m_out,
                        const std::vector<std::array<vec3, 3>> *rest = nullptr) {
  real a = 1;
  vecX r;
  for (int k = 0; k < std::max(cfg.line_search, 1); ++k, a *= real(0.5)) {
    std::vector<vec3> xn = sys.gather(sys.frames(x, xi, a), x);
    const real mn = merit(t, xn, goal, sys.weights(), r, &cfg, sys.tau(), rest);
    if (std::isfinite(mn) && mn < m0) {
      x_out = std::move(xn);
      m_out = mn;
      return a;
    }
  }
  return 0;
}

struct solve_report {
  int iterations = 0;
  std::vector<real> merit;
};

/// The standalone solver: restarts (a fresh system at the current positions)
/// every `restart` iterates, after an iterate that cuts the merit by less
/// than 10x (a frozen factor contracts linearly, a fresh one quadratically),
/// and on a stall; one backsolve and a line search per iterate,
/// `goal` the per-face target read-outs. `rest` makes the material elastic
/// against those triangles (with `cfg.km`).
inline solve_report solve(const kusama::stencil_table &t, std::vector<vec3> &x, const std::vector<vec3> &goal,
                          const config &cfg, int iterations = 120, int restart = 30, real tol = real(1e-13),
                          const std::vector<std::array<vec3, 3>> *rest = nullptr) {
  const topology tp = build_topology(t);
  solve_report rep;
  std::unique_ptr<system> sys;
  vecX r;
  real m = 0;
  bool slow = false;
  for (int it = 0; it < iterations; ++it) {
    bool fresh = false;
    if (it % restart == 0 || !sys || slow) {
      sys = std::make_unique<system>(t, tp, x, cfg);
      fresh = true;
      m = merit(t, x, goal, sys->weights(), r, &cfg, sys->tau(), rest);
      if (rep.merit.empty())
        rep.merit.push_back(m);
    }
    std::vector<vec3> xn;
    real mn = 0;
    std::vector<vec3> Lc = rest ? current_log_stretch(t, x, *rest) : std::vector<vec3>{};
    real a = line_search(t, *sys, x, sys->solve(r, rest ? &Lc : nullptr), goal, m, cfg, xn, mn, rest);
    if (a == real(0) && !fresh) {
      sys = std::make_unique<system>(t, tp, x, cfg);
      m = merit(t, x, goal, sys->weights(), r, &cfg, sys->tau(), rest);
      a = line_search(t, *sys, x, sys->solve(r, rest ? &Lc : nullptr), goal, m, cfg, xn, mn, rest);
    }
    if (a == real(0))
      break;
    x = std::move(xn);
    slow = mn > real(0.1) * m;
    m = mn;
    merit(t, x, goal, sys->weights(), r, &cfg, sys->tau(), rest);
    rep.merit.push_back(m);
    ++rep.iterations;
    if (m < tol)
      break;
  }
  return rep;
}

/// Each face's goal from the edge store: the mean of its three edge tensors
/// (edge frame `(t_e, N_f x t_e)`, `t_e` from the lower vertex id) turned
/// into the face frame.
inline std::vector<vec3> store_face_goals(const kusama::stencil_table &t, const std::vector<vec3> &x,
                                          const std::function<vec3(int)> &store) {
  std::vector<vec3> g(static_cast<size_t>(t.n_faces), vec3::Zero());
  for (int f = 0; f < t.n_faces; ++f) {
    const auto &id = t.ids[static_cast<size_t>(f)];
    vec3 T, B, N;
    real A;
    face_frame(x[static_cast<size_t>(id[0])], x[static_cast<size_t>(id[1])], x[static_cast<size_t>(id[2])], T, B, N, A);
    for (int i = 0; i < 3; ++i) {
      const int a = id[static_cast<size_t>(i)], b = id[static_cast<size_t>((i + 1) % 3)];
      const vec3 te = (x[static_cast<size_t>(std::max(a, b))] - x[static_cast<size_t>(std::min(a, b))]).normalized();
      g[static_cast<size_t>(f)] +=
          turn_tensor(store(t.face_edges[static_cast<size_t>(f)][static_cast<size_t>(i)]), te.dot(T), te.dot(B)) /
          real(3);
    }
  }
  return g;
}

/// Calibrates the coupling's gather rows against the cotan `bend` row, as
/// `kusama::weingarten::k_coupling` does the `weingarten` rows: at equal
/// weights their mean position stiffness matches (unscaled 0.029 of `bend`
/// on both test tori, `exponent_coupling_matches_bend_at_equal_weight`), so
/// one `w_coupling` means the same under either solver.
constexpr real k_gather = real(34.5);

/// A factor carried between solver steps (owned by the caller, one per
/// mesh): reused while the faces match and it is younger than
/// `config::refresh_every` steps.
struct factor_cache {
  std::shared_ptr<const factor> F;
  int age = 0;
  /// A refactor after a stall stalled too (a large-residual balance: Gauss-
  /// Newton has no descent left); stalls wait for the age refresh.
  bool futile = false;
};

} // namespace exponent

/// The exponent solver as a projection row on positions (`blocks[0]`), the
/// curvature targets read from the edge store (`blocks[1]`, read only). Rows:
/// the CR derivative `w sqrt(k_gather M_e) D` at the frame start (constant: with the
/// other rows the global step is the Dirichlet gather). `project` takes one
/// exponent step from the live positions (residual, backsolve, per-face
/// exponentials, line search on the merit) and writes `D x + G_step - G(x)`,
/// so a zero step pulls nowhere; up to `steps_per_frame` steps per change of
/// the goal, the last target held in between. Targets in increment form, so
/// the row is idle at the frame start:
/// `s_g,f = s_f(x0) + goal_f(store) - goal_f(store at x0)`; with
/// `residual::principal` the same on the eigenvalues of each (Eigensatz:
/// the store's change of principal curvature, directions free). `cache`
/// carries the factor between frames.
class exponent_coupling : public block_constraint {
public:
  DEFINE_CREATE_FUNC(exponent_coupling)

  exponent_coupling(const kusama::stencil_table *table, const std::vector<vec3> &x0,
                    const std::vector<vec3> &h0, real w, std::vector<sim_block::ptr> blocks,
                    exponent::config cfg, std::shared_ptr<const std::vector<std::array<vec3, 3>>> rest = nullptr,
                    std::shared_ptr<exponent::factor_cache> cache = nullptr)
      : block_constraint({}, w, std::move(blocks)), _t(table), _cfg(cfg), _rest(std::move(rest)),
        _cache(std::move(cache)) {
    _tp = exponent::build_topology(*_t);
    std::shared_ptr<const exponent::factor> F;
    if (_cache && _cache->F && _cache->F->matches(*_t) && _cache->age < std::max(_cfg.refresh_every, 1))
      F = _cache->F;
    _sys = std::make_unique<exponent::system>(*_t, _tp, x0, _cfg, F);
    if (_cache) {
      const bool reused = _sys->factor_ptr() == _cache->F;
      _cache->age = reused ? _cache->age + 1 : 1;
      _cache->futile = reused && _cache->futile;
      _cache->F = _sys->factor_ptr();
    }
    _s0 = exponent::face_readouts(*_t, x0);
    _g0 = exponent::store_face_goals(*_t, x0, [&](int e) { return h0[static_cast<size_t>(e)]; });
    if (_cfg.kind == exponent::residual::principal)
      for (size_t f = 0; f < _s0.size(); ++f) {
        _s0[f] = exponent::principal_values(_s0[f]);
        _g0[f] = exponent::principal_values(_g0[f]);
      }
    _x0 = x0;
    _rest_initialized = true;
  }

  std::string name() override { return "exponent_coupling"; }

  void fill_A(index_t &id0, std::vector<trip> &triplets) override {
    _id0 = id0;
    const auto &D = _sys->derivative();
    for (size_t e = 0; e < D.size(); ++e) {
      const real s = _w * std::sqrt(exponent::k_gather * std::max(D[e].mass, real(0)));
      for (int j = 0; j < 2; ++j)
        for (int c = 0; c < 3; ++c) {
          for (int k = 0; k < 4; ++k) {
            const int v = D[e].v[static_cast<size_t>(k)];
            if (v < 0)
              continue;
            const real coef = j == 0 ? D[e].par[static_cast<size_t>(k)] : D[e].perp[static_cast<size_t>(k)];
            if (coef != real(0))
              triplets.push_back(trip(id0, _blocks[0]->get_offset_idx(v) + c, s * coef));
          }
          ++id0;
        }
    }
  }

  void project(const vecX &q, vecX &p) override {
    const int n = _t->n_verts, nE = _t->n_edges;
    std::vector<vec3> x(static_cast<size_t>(n));
    for (int v = 0; v < n; ++v)
      x[static_cast<size_t>(v)] = _blocks[0]->get_vec3(v, q);
    std::vector<vec3> goal = _s0;
    if (_blocks.size() > 1 && _blocks[1]) {
      const std::vector<vec3> g =
          exponent::store_face_goals(*_t, _x0, [&](int e) { return vec3(_blocks[1]->get_vec3(e, q)); });
      const bool principal = _cfg.kind == exponent::residual::principal;
      for (size_t f = 0; f < goal.size(); ++f)
        goal[f] += (principal ? exponent::principal_values(g[f]) : g[f]) - _g0[f];
    }
    // The store settles over the solver's iterates: a goal that moved re-arms
    // the step budget; a settled one holds the reached frames.
    real moved = 0, scale = 0;
    for (size_t f = 0; f < goal.size(); ++f) {
      moved = std::max(moved, (goal[f] - (_goal_last.empty() ? _s0[f] : _goal_last[f])).cwiseAbs().maxCoeff());
      scale = std::max(scale, goal[f].cwiseAbs().maxCoeff());
    }
    if (moved > real(1e-6) * std::max(scale, real(1e-12))) {
      _goal_last = goal;
      _steps = 0;
      _done = false;
    }
    if (_done || _steps >= _cfg.steps_per_frame) {
      // Keep pulling toward the last target; once no step improves, x is it.
      write(_held.empty() ? _sys->derivative_of(x) : _held, p);
      return;
    }
    ++_steps;
    const std::vector<vec3> Lc = _rest ? exponent::current_log_stretch(*_t, x, *_rest) : std::vector<vec3>{};
    vecX r, xi;
    real m0 = 0, mn = 0, a = 0;
    std::vector<vec3> xn;
    for (int attempt = 0; attempt < 2; ++attempt) {
      m0 = exponent::merit(*_t, x, goal, _sys->weights(), r, &_cfg, _sys->tau(), _rest.get());
      xi = _sys->solve(r, _rest ? &Lc : nullptr);
      a = exponent::line_search(*_t, *_sys, x, xi, goal, m0, _cfg, xn, mn, _rest.get());
      if (a > real(0) || _sys->fresh()) {
        if (a == real(0) && attempt > 0 && _cache)
          _cache->futile = true;
        break;
      }
      if (_cache && _cache->futile)
        break;
      // A reused factor that cannot descend: restart at the live positions.
      _sys->refactor(x);
      if (_cache) {
        _cache->F = _sys->factor_ptr();
        _cache->age = 1;
      }
    }
    _last_merit = a > real(0) ? mn : m0;
    _done = a == real(0);
    (void)nE;
    // Rows `D x = D x_k + G_step - G(x_k)`: no pull at a = 0 wherever x is.
    _held = _sys->derivative_of(x);
    if (a > real(0)) {
      const std::vector<vec3> Gs = _sys->frames(x, xi, a), Gx = _sys->frames(x);
      for (size_t i = 0; i < _held.size(); ++i)
        _held[i] += Gs[i] - Gx[i];
    }
    write(_held, p);
  }

  real last_merit() const { return _last_merit; }
  int steps() const { return _steps; }

private:
  void write(const std::vector<vec3> &G, vecX &p) const {
    const int nE = _t->n_edges;
    const auto &D = _sys->derivative();
    index_t row = _id0;
    for (int e = 0; e < nE; ++e) {
      const real s = _w * std::sqrt(exponent::k_gather * std::max(D[static_cast<size_t>(e)].mass, real(0)));
      for (int j = 0; j < 2; ++j)
        for (int c = 0; c < 3; ++c)
          p[row++] = s * G[static_cast<size_t>(j * nE + e)][c];
    }
  }

  const kusama::stencil_table *_t = nullptr;
  exponent::config _cfg;
  std::shared_ptr<const std::vector<std::array<vec3, 3>>> _rest;
  exponent::topology _tp;
  std::unique_ptr<exponent::system> _sys;
  std::shared_ptr<exponent::factor_cache> _cache;
  std::vector<vec3> _s0, _g0, _x0, _goal_last, _held;
  real _last_merit = 0;
  int _steps = 0;
  bool _done = false;
};
} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif
