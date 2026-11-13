#ifndef GAUDI_HEPWORTH_BLOCK_KITE_GROWTH_HPP
#define GAUDI_HEPWORTH_BLOCK_KITE_GROWTH_HPP

/// Moebius kites as an edge growth (scripts/sympy/cr_poisson/MATH.md
/// 19.20). A Moebius map is locally a conformal scale `e^u`; on a mesh its
/// discrete form scales each edge by its endpoints' factors,
/// `l_ij -> e^{(u_i + u_j) / 2} l_ij`, so a triangle whose vertices carry
/// different `u` changes shape (the kite at the triangle level). The
/// per-frame vertex log-scale `u` comes from `source`:
///
/// - `mean_curvature`: `u = alpha (-H - level)`, each element's Moebius map
///   anchored on its mean-curvature sphere: convex grows, concave shrinks,
///   saddles (`H ~ 0`) stay (outward normals make convex `H` negative, hence
///   the sign); `anisotropy`: `u = alpha (|z| - level)`. With `balanced` the
///   mean is removed (redistributes scale, no net growth).
/// - `tangent`: the Dirichlet edge tangents (sign-aligned to the previous
///   frame's) turned by `angle` about the normal (0 along them, trapezoids;
///   45 degrees, kites), times `alpha |z|`, as a growth gradient; one
///   Poisson solve `argmin sum_f A_f |grad u - g_f|^2` keeps its integrable
///   part, mean zero (a direction field has no level of its own).
///
/// Per edge an `edge_strain` row with goal length `e^{(u_i + u_j) / 2}`
/// times its frame-start length.
///
/// Rows are rebuilt each frame against the frame-start mesh, so the growth
/// accumulates plastically. Positions only.

#include "gaudi/hepworth/block/shell_constraints.hpp"
#include "gaudi/kusama/stencil_table.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <vector>

#include <Eigen/SparseCholesky>

namespace gaudi {
namespace hepworth {
namespace block {

/// Where the log-scale comes from: the store's mean curvature or anisotropy
/// `|z|` directly, or the Dirichlet edge tangents through a Poisson solve.
enum class kite_source { tangent, mean_curvature, anisotropy };

/// `alpha`: log-scale per unit curvature per frame (`mean_curvature`,
/// `anisotropy`), or growth gradient per unit `|z|` per frame (`tangent`);
/// `angle`: turn of the tangents (radians, `tangent` only); `max_u`: cap on
/// the per-frame `|u|`; `source`; `level`: the curvature at which nothing
/// changes (scalar sources); `balanced`: remove the mean (scalar sources;
/// the tangent source is always balanced).
struct kite_target {
  real alpha = 0.02;
  real angle = real(M_PI) / 4;
  real max_u = 0.1;
  kite_source source = kite_source::mean_curvature;
  real level = 0;
  bool balanced = false;
};

namespace kite_detail {

inline vec3 face_normal(const kusama::stencil_table &t, const std::vector<vec3> &x, int f, real *area = nullptr) {
  const auto &id = t.ids[static_cast<size_t>(f)];
  const vec3 a = x[static_cast<size_t>(id[0])], b = x[static_cast<size_t>(id[1])], c = x[static_cast<size_t>(id[2])];
  vec3 n = (b - a).cross(c - a);
  const real l = n.norm();
  if (area)
    *area = real(0.5) * l;
  return l > real(0) ? vec3(n / l) : vec3::Zero();
}

/// Gradients of the three hat functions on face `f` (zero when degenerate).
inline std::array<vec3, 3> hat_gradients(const kusama::stencil_table &t, const std::vector<vec3> &x, int f,
                                         real &A, vec3 &N) {
  const auto &id = t.ids[static_cast<size_t>(f)];
  N = face_normal(t, x, f, &A);
  std::array<vec3, 3> gr{vec3::Zero(), vec3::Zero(), vec3::Zero()};
  if (!(A > real(0)))
    return gr;
  for (int k = 0; k < 3; ++k) {
    const vec3 p = x[static_cast<size_t>(id[static_cast<size_t>((k + 1) % 3)])];
    const vec3 q = x[static_cast<size_t>(id[static_cast<size_t>((k + 2) % 3)])];
    gr[static_cast<size_t>(k)] = N.cross(q - p) / (real(2) * A);
  }
  return gr;
}

inline std::array<int, 2> edge_ends(const kusama::stencil_table &t, int e) {
  for (int s = 0; s < 2; ++s) {
    const int f = t.edge_faces[static_cast<size_t>(e)][static_cast<size_t>(s)];
    const int i = t.edge_local[static_cast<size_t>(e)][static_cast<size_t>(s)];
    if (f < 0 || i < 0)
      continue;
    const auto &id = t.ids[static_cast<size_t>(f)];
    return {id[static_cast<size_t>(i)], id[static_cast<size_t>((i + 1) % 3)]};
  }
  return {-1, -1};
}

inline real store_z(const vec3 &h) { return std::hypot(real(0.5) * (h[0] - h[2]), h[1]); }

/// Per vertex, the mean over its edges of the store's convexity `-H`
/// (positive where both curvatures bend away from the outward normal) or
/// of its anisotropy `|z|`.
inline std::vector<real> store_vertex_scalar(const kusama::stencil_table &t, const std::vector<vec3> &store,
                                             bool aniso) {
  const size_t n = static_cast<size_t>(t.n_verts), nE = static_cast<size_t>(t.n_edges);
  std::vector<real> s(n, real(0)), c(n, real(0));
  for (size_t e = 0; e < nE && e < store.size(); ++e) {
    const std::array<int, 2> v = edge_ends(t, static_cast<int>(e));
    if (v[0] < 0)
      continue;
    const real val = aniso ? store_z(store[e]) : -real(0.5) * (store[e][0] + store[e][2]);
    for (int k = 0; k < 2; ++k) {
      s[static_cast<size_t>(v[static_cast<size_t>(k)])] += val;
      c[static_cast<size_t>(v[static_cast<size_t>(k)])] += 1;
    }
  }
  for (size_t i = 0; i < n; ++i)
    s[i] = c[i] > 0 ? s[i] / c[i] : real(0);
  return s;
}

/// Per vertex, its share of the area (a third of each incident face).
inline std::vector<real> vertex_areas(const kusama::stencil_table &t, const std::vector<vec3> &x) {
  std::vector<real> a(static_cast<size_t>(t.n_verts), real(0));
  for (int f = 0; f < t.n_faces; ++f) {
    real A = 0;
    face_normal(t, x, f, &A);
    for (int k = 0; k < 3; ++k)
      a[static_cast<size_t>(t.ids[static_cast<size_t>(f)][static_cast<size_t>(k)])] += A / real(3);
  }
  return a;
}

} // namespace kite_detail

/// The vertex log-scale whose gradient best matches the per-edge vectors
/// `g` (tangential parts used), area weighted, mean zero.
inline std::vector<real> kite_poisson(const kusama::stencil_table &t, const std::vector<vec3> &x,
                                      const std::vector<vec3> &g) {
  const int n = t.n_verts;
  std::vector<real> u(static_cast<size_t>(std::max(n, 0)), real(0));
  if (n <= 0 || static_cast<int>(x.size()) != n || static_cast<int>(g.size()) != t.n_edges)
    return u;
  std::vector<Eigen::Triplet<real>> K;
  vecX b = vecX::Zero(n);
  real atot = 0;
  std::vector<real> vA(static_cast<size_t>(n), real(0));
  for (int f = 0; f < t.n_faces; ++f) {
    const auto &id = t.ids[static_cast<size_t>(f)];
    real A = 0;
    vec3 N;
    const std::array<vec3, 3> gr = kite_detail::hat_gradients(t, x, f, A, N);
    if (!(A > real(0)))
      continue;
    vec3 gf = vec3::Zero();
    for (int i = 0; i < 3; ++i)
      gf += g[static_cast<size_t>(t.face_edges[static_cast<size_t>(f)][static_cast<size_t>(i)])];
    gf /= real(3);
    gf -= N * N.dot(gf);
    for (int k = 0; k < 3; ++k) {
      const int vk = id[static_cast<size_t>(k)];
      b[vk] += A * gr[static_cast<size_t>(k)].dot(gf);
      vA[static_cast<size_t>(vk)] += A / real(3);
      for (int l = 0; l < 3; ++l)
        K.emplace_back(vk, id[static_cast<size_t>(l)], A * gr[static_cast<size_t>(k)].dot(gr[static_cast<size_t>(l)]));
    }
    atot += A;
  }
  Eigen::SparseMatrix<real> M(n, n);
  M.setFromTriplets(K.begin(), K.end());
  real diag = 0;
  for (int i = 0; i < n; ++i)
    diag += M.coeff(i, i);
  Eigen::SparseMatrix<real> I(n, n);
  I.setIdentity();
  M += I * (real(1e-9) * diag / std::max(n, 1));
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<real>> ldlt(M);
  if (ldlt.info() != Eigen::Success)
    return u;
  const vecX s = ldlt.solve(b);
  real mean = 0;
  for (int i = 0; i < n; ++i)
    mean += vA[static_cast<size_t>(i)] * s[i];
  mean /= std::max(atot, real(1e-300));
  for (int i = 0; i < n; ++i)
    u[static_cast<size_t>(i)] = s[i] - mean;
  return u;
}

/// The tangent source's growth gradient per edge: `dirs` sign-aligned to
/// `carried` (globally, so the field's smoothness is kept), turned by
/// `angle` about the edge normal, times `alpha |z|`.
inline std::vector<vec3> kite_gradient(const kusama::stencil_table &t, const std::vector<vec3> &x,
                                       std::vector<vec3> dirs, const std::vector<vec3> &carried,
                                       const std::vector<vec3> &store, const kite_target &spec) {
  const size_t nE = static_cast<size_t>(t.n_edges);
  std::vector<vec3> g(nE, vec3::Zero());
  if (store.size() != nE || dirs.size() != nE)
    return g;
  if (carried.size() == nE) {
    real s = 0;
    for (size_t e = 0; e < nE; ++e)
      s += dirs[e].dot(carried[e]);
    if (s < real(0))
      for (vec3 &d : dirs)
        d = -d;
  }
  const real ca = std::cos(spec.angle), sa = std::sin(spec.angle);
  for (size_t e = 0; e < nE; ++e) {
    vec3 N = vec3::Zero();
    for (int s = 0; s < 2; ++s) {
      const int f = t.edge_faces[e][static_cast<size_t>(s)];
      if (f >= 0)
        N += kite_detail::face_normal(t, x, f);
    }
    if (!(N.squaredNorm() > real(0)))
      continue;
    N.normalize();
    vec3 d = dirs[e] - N * N.dot(dirs[e]);
    const real l = d.norm();
    if (!(l > real(1e-12)))
      continue;
    d /= l;
    g[e] = spec.alpha * kite_detail::store_z(store[e]) * (ca * d + sa * N.cross(d));
  }
  return g;
}

/// The per-frame log-scale, capped at `max_u`.
inline std::vector<real> kite_log_scale(const kusama::stencil_table &t, const std::vector<vec3> &x,
                                        const std::vector<vec3> &dirs, const std::vector<vec3> &carried,
                                        const std::vector<vec3> &store, const kite_target &spec) {
  std::vector<real> u;
  if (spec.source == kite_source::tangent) {
    u = kite_poisson(t, x, kite_gradient(t, x, dirs, carried, store, spec));
  } else {
    if (store.size() != static_cast<size_t>(t.n_edges) || static_cast<int>(x.size()) != t.n_verts)
      return std::vector<real>(static_cast<size_t>(std::max(t.n_verts, 0)), real(0));
    u = kite_detail::store_vertex_scalar(t, store, spec.source == kite_source::anisotropy);
    for (real &v : u)
      v = spec.alpha * (v - spec.level);
    if (spec.balanced) {
      const std::vector<real> a = kite_detail::vertex_areas(t, x);
      real su = 0, sa = 0;
      for (size_t i = 0; i < u.size(); ++i) {
        su += a[i] * u[i];
        sa += a[i];
      }
      const real mean = su / std::max(sa, real(1e-300));
      for (real &v : u)
        v -= mean;
    }
  }
  if (spec.max_u > real(0))
    for (real &v : u)
      v = std::clamp(v, -spec.max_u, spec.max_u);
  return u;
}

/// One `edge_strain` row per edge with goal length `e^{(u_i + u_j) / 2}`
/// times its frame-start length. Returns the log-scale used (empty when
/// nothing was added).
inline std::vector<real> init_kite_growth(const kusama::stencil_table &t, const std::vector<vec3> &x,
                                          const std::vector<vec3> &dirs, const std::vector<vec3> &carried,
                                          const std::vector<vec3> &store, const kite_target &spec, real w,
                                          std::vector<projection_constraint::ptr> &constraints,
                                          std::vector<sim_block::ptr> blocks) {
  if (!(w > real(0)) || static_cast<int>(x.size()) != t.n_verts)
    return {};
  const std::vector<real> u = kite_log_scale(t, x, dirs, carried, store, spec);
  for (int e = 0; e < t.n_edges; ++e) {
    const std::array<int, 2> v = kite_detail::edge_ends(t, e);
    if (v[0] < 0)
      continue;
    const real l0 = (x[static_cast<size_t>(v[1])] - x[static_cast<size_t>(v[0])]).norm();
    if (!(l0 > real(0)))
      continue;
    auto row = edge_strain::create(std::vector<index_t>{v[0], v[1]}, x, w, blocks);
    row->set_goal_length(l0 * std::exp(real(0.5) * (u[static_cast<size_t>(v[0])] + u[static_cast<size_t>(v[1])])));
    constraints.push_back(row);
  }
  return u;
}

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif
