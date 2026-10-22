#ifndef GAUDI_DUCHAMP_TP_SOBOLEV_JACOBI_HPP
#define GAUDI_DUCHAMP_TP_SOBOLEV_JACOBI_HPP

/// Sobolev metric A = α_B B + α_0 B_0. Smoothing is Eigen CG (Appendix B.3.1).
/// The coarsest system is a dense saddle factorization.
/// B / B_0 follow Crane edge–edge assembly (§5.2, Eq. 21 and the low-order
/// term): vertex fields lift to edges, then the adjoint scatters back.
/// B lifts each coordinate by the unit tangent, (Δu/ℓ) T, and scatters T·r.
/// Leaf K is the two edge midpoints, the sample in AfFullProduct and
/// AfFullProductLow. The paper trapezoid {0, 1} stays on the energy Ê.
/// A far node stays one tangent-point lump. B0: K_IJ = (k+k_swap)·ℓ_I·ℓ_J
/// with k_2^{p+4}.
/// B and B0 both use r_I = Σ w (φ_I − φ_J) (Appendix B.3.4 Eq. 31).
/// A far node substitutes Σ (ℓ_J φ_J) / L for φ_J (Appendix B.1). The pyramid
/// stores sums, and w already contains L, so the partner is (w/L) Σ (ℓ φ).
/// B0 lifts with the edge average. Kernel lengths are
/// geometric ℓ. Coarse levels rediscretize A on that level's rod.
/// Writeup: docs/sobolev_B_tangent_filter.md.
/// The Sobolev right-hand side is dE for the vertex tangent-point energy.
/// Each visit is dual_i * m_j * |n_i|^α / |r|^β. bh_eps ≤ 0 walks every
/// vertex; bh_eps > 0 lumps a far cluster as one mass at its center of mass.

#include "gaudi/asawa/rod/dynamic.hpp"
#include "gaudi/asawa/rod/residual_forest.hpp"
#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/calder/rod_integrators.hpp"
#include "gaudi/calder/tree_code.hpp"
#include "gaudi/common.h"
#include "gaudi/duchamp/arnoldi.hpp"
#include "gaudi/duchamp/demo_log.hpp"
#include "gaudi/duchamp/modules/tp_config.hpp"

#include <iostream>

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/LU>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

#include <algorithm>
#include <cmath>
#include <functional>
#include <memory>
#include <vector>

namespace gaudi {
namespace duchamp {

namespace sobolev_detail {

inline std::vector<real> sample_masses(const asawa::rod::rod &rod) {
  const size_t n = rod.corner_count();
  const auto &L = rod.l0();
  std::vector<real> m(n, real(0));
  for (size_t i = 0; i < n; ++i) {
    const auto ci = asawa::rod::corner_id(static_cast<int>(i));
    const auto cn = rod.next(ci);
    const auto cp = rod.prev(ci);
    if (cn >= asawa::rod::corner_id(0) &&
        static_cast<size_t>(static_cast<int>(ci)) < L.size())
      m[i] += real(0.5) * std::max(L[static_cast<size_t>(static_cast<int>(ci))],
                                   real(0));
    if (cp >= asawa::rod::corner_id(0)) {
      const size_t ip = static_cast<size_t>(static_cast<int>(cp));
      if (ip < L.size())
        m[i] += real(0.5) * std::max(L[ip], real(0));
    }
    if (!(m[i] > real(0)))
      m[i] = real(1);
  }
  return m;
}

/// Parameters along an edge, t=0 at the start vertex and t=1 at the end.
/// Paper trapezoid: the two endpoints, each with weight 1/2.
struct edge_quad {
  std::vector<real> t{real(0), real(1)};
  std::vector<real> w{real(0.5), real(0.5)};
};

inline const edge_quad &trapezoid_edge_quad() {
  static const edge_quad q;
  return q;
}

/// Packed edges in get_edge_vert_ids() order, the tree's primitive order.
struct edge_table {
  size_t n = 0;      ///< edge count
  size_t n_vert = 0; ///< corner count
  std::vector<int> ep0;
  std::vector<int> ep1;
  std::vector<real> ell; ///< geometric length ℓ_I
  std::vector<vec3> T;   ///< unit chord
  std::vector<vec3> mid; ///< opening-test query, the edge midpoint
};

inline std::vector<index_t> edge_ordinals(size_t n) {
  std::vector<index_t> ids(n);
  for (size_t i = 0; i < n; ++i)
    ids[i] = static_cast<index_t>(i);
  return ids;
}

inline edge_table build_edges(asawa::rod::rod &rod) {
  edge_table t;
  t.n_vert = rod.corner_count();
  const std::vector<index_t> ev = rod.get_edge_vert_ids();
  t.n = ev.size() / 2;
  t.ep0.resize(t.n);
  t.ep1.resize(t.n);
  t.ell.resize(t.n);
  t.T.resize(t.n);
  t.mid.resize(t.n);
  const auto &x = rod.x();
  for (size_t e = 0; e < t.n; ++e) {
    const int i0 = static_cast<int>(ev[2 * e]);
    const int i1 = static_cast<int>(ev[2 * e + 1]);
    t.ep0[e] = i0;
    t.ep1[e] = i1;
    const vec3 d = x[static_cast<size_t>(i1)] - x[static_cast<size_t>(i0)];
    const real len = d.norm();
    t.ell[e] = std::max(len, real(1e-12));
    t.T[e] = (len > real(1e-12)) ? vec3(d / len) : vec3::Zero();
    t.mid[e] = real(0.5) * (x[static_cast<size_t>(i0)] + x[static_cast<size_t>(i1)]);
  }
  return t;
}

inline std::vector<vec3> edge_avg_E(const edge_table &E,
                                    const std::vector<vec3> &u) {
  std::vector<vec3> Eu(E.n, vec3::Zero());
  for (size_t i = 0; i < E.n; ++i)
    Eu[i] = real(0.5) * (u[static_cast<size_t>(E.ep0[i])] +
                         u[static_cast<size_t>(E.ep1[i])]);
  return Eu;
}

/// Per-coordinate slope (u1 − u0) / ℓ. §5.2.1 multiplies this by T before the
/// residual and dots with T on the scatter.
inline std::vector<vec3> edge_diff_D(const edge_table &E,
                                     const std::vector<vec3> &u) {
  std::vector<vec3> Du(E.n, vec3::Zero());
  for (size_t i = 0; i < E.n; ++i)
    Du[i] = (u[static_cast<size_t>(E.ep1[i])] -
             u[static_cast<size_t>(E.ep0[i])]) /
            E.ell[i];
  return Du;
}

inline void scatter_ET(const edge_table &E, const std::vector<vec3> &r_edge,
                       std::vector<vec3> &out) {
  out.assign(E.n_vert, vec3::Zero());
  for (size_t i = 0; i < E.n; ++i) {
    const vec3 h = real(0.5) * r_edge[i];
    out[static_cast<size_t>(E.ep0[i])] += h;
    out[static_cast<size_t>(E.ep1[i])] += h;
  }
}

/// Paper tangent-point kernel k_α^β = |T×dp|^α / (|dp|^β + ε^β).
/// B0 uses k_2^{2σ+5} (Eq. 32): with B's distance power p=2σ+1 ⇒ α=2, β=p+4.
inline real paper_tp_k(const vec3 &dp, const vec3 &T, real eps, real p) {
  if (T.squaredNorm() < real(1e-24) || dp.squaredNorm() < real(1e-24))
    return real(0);
  vec3 Ti = T.normalized();
  const real cross_n = Ti.cross(dp).norm();
  const real r = dp.norm();
  const real alpha = real(2);
  const real beta = p + real(4);
  const real denom =
      std::pow(r, beta) + std::pow(std::max(eps, real(1e-12)), beta);
  if (!(denom > real(0)))
    return real(0);
  return std::pow(cross_n, alpha) / denom;
}

/// One direction. B0 is k_2^{p+4}. B is k_0^p = 1 / (|r|^p + ε^p).
inline real kernel_k(const vec3 &dp, const vec3 &T, real eps, real p,
                     bool use_tp_k) {
  if (dp.squaredNorm() < real(1e-24))
    return real(0);
  if (use_tp_k)
    return paper_tp_k(dp, T, eps, p);
  const real r = dp.norm();
  const real denom =
      std::pow(r, p) + std::pow(std::max(eps, real(1e-12)), p);
  if (!(denom > real(0)))
    return real(0);
  return real(1) / denom;
}

inline vec3 quad_point(const vec3 &a, const vec3 &b, real t) {
  return (real(1) - t) * a + t * b;
}

/// Far lump: one tangent-point, both directions, times ℓ_I L.
inline real directed_pair_weight(const vec3 &pi, const vec3 &pj, const vec3 &Ti,
                                 const vec3 &Tj, real li, real lj, real eps,
                                 real p, bool use_tp_k) {
  if (!(li > real(0)) || !(lj > real(0)))
    return real(0);
  const real k = kernel_k(pj - pi, Ti, eps, p, use_tp_k) +
                 kernel_k(pi - pj, Tj, eps, p, use_tp_k);
  return li * lj * k;
}

/// Leaf pair. k̂ = Σ_a Σ_b w_a w_b k, then K = (k̂(I→J) + k̂(J→I)) ℓ_I ℓ_J.
inline real quad_pair_weight(const vec3 &a0, const vec3 &a1, const vec3 &Ta,
                             const vec3 &b0, const vec3 &b1, const vec3 &Tb,
                             real li, real lj, real eps, real p, bool use_tp_k,
                             const edge_quad &q) {
  if (!(li > real(0)) || !(lj > real(0)) || q.t.empty())
    return real(0);
  real khat = real(0);
  const size_t nq = q.t.size();
  for (size_t ia = 0; ia < nq; ++ia) {
    const real wa = ia < q.w.size() ? q.w[ia] : real(0);
    const vec3 pa = quad_point(a0, a1, q.t[ia]);
    for (size_t ib = 0; ib < nq; ++ib) {
      const real wb = ib < q.w.size() ? q.w[ib] : real(0);
      const vec3 pb = quad_point(b0, b1, q.t[ib]);
      const real w = wa * wb;
      khat += w * kernel_k(pb - pa, Ta, eps, p, use_tp_k);
      khat += w * kernel_k(pa - pb, Tb, eps, p, use_tp_k);
    }
  }
  return li * lj * khat;
}

inline bool leaf_edges_adjacent_or_self(const edge_table &E, size_t ii, int j0,
                                         int j1) {
  if (j0 < 0 || j1 < 0)
    return false;
  const int a0 = E.ep0[ii];
  const int a1 = E.ep1[ii];
  return a0 == j0 || a0 == j1 || a1 == j0 || a1 == j1;
}

inline bool leaf_edges_same(const edge_table &E, size_t ii, int j0, int j1) {
  if (j0 < 0 || j1 < 0)
    return false;
  const int a0 = E.ep0[ii];
  const int a1 = E.ep1[ii];
  return (a0 == j0 && a1 == j1) || (a0 == j1 && a1 == j0);
}

/// Leaf K is the two edge midpoints. Far uses the lumped point. `quad` is the
/// energy trapezoid and is not read here.
inline real sample_pair_weight(const edge_table &E, size_t ii, const vec3 &pi,
                               const vec3 &pj, const vec3 &Tj, real lj, int j0,
                               int j1, const std::vector<vec3> &x,
                               calder::Node_Type node_type, real eps, real p,
                               bool use_tp_k, const edge_quad &quad,
                               bool drop_adjacent = true) {
  (void)quad;
  const real li = E.ell[ii];
  if (node_type == calder::LEAF) {
    if (leaf_edges_same(E, ii, j0, j1))
      return real(0);
    if (drop_adjacent && leaf_edges_adjacent_or_self(E, ii, j0, j1))
      return real(0);
    if (j0 < 0 || j1 < 0)
      return real(0);
    const size_t b0 = static_cast<size_t>(j0);
    const size_t b1 = static_cast<size_t>(j1);
    if (b0 >= x.size() || b1 >= x.size())
      return real(0);
    const vec3 mid_j = real(0.5) * (x[b0] + x[b1]);
    return directed_pair_weight(E.mid[ii], mid_j, E.T[ii], Tj, li, lj, eps, p,
                                use_tp_k);
  }
  return directed_pair_weight(pi, pj, E.T[ii], Tj, li, lj, eps, p, use_tp_k);
}

/// Energy kernel k_2^{p+2}. paper_tp_k stays k_2^{p+4} for B0.
struct energy_dk {
  real k = 0;
  vec3 dk_dp = vec3::Zero();
  vec3 dk_dT = vec3::Zero();
};

inline energy_dk energy_kernel_deriv(const vec3 &dp, const vec3 &T, real eps,
                                     real p) {
  energy_dk out;
  if (T.squaredNorm() < real(1e-24) || dp.squaredNorm() < real(1e-24))
    return out;
  const vec3 Ti = T.normalized();
  const real r = dp.norm();
  const vec3 v = Ti.cross(dp);
  const real c2 = v.squaredNorm();
  const real beta = p + real(2);
  const real denom =
      std::pow(r, beta) + std::pow(std::max(eps, real(1e-12)), beta);
  if (!(denom > real(0)) || !std::isfinite(denom))
    return out;
  out.k = c2 / denom;
  const vec3 dN_dp = real(2) * v.cross(Ti);
  const vec3 dN_dT = real(2) * dp.cross(v);
  const vec3 dD_dp = beta * std::pow(r, beta - real(2)) * dp;
  const real invD2 = real(1) / (denom * denom);
  out.dk_dp = (dN_dp * denom - c2 * dD_dp) * invD2;
  out.dk_dT = dN_dT / denom;
  if (!std::isfinite(out.k) || !out.dk_dp.allFinite() || !out.dk_dT.allFinite())
    return energy_dk{};
  return out;
}

/// Full symmetric leaf term (k̂(I→J)+k̂(J→I)) ℓ_I ℓ_J and its endpoint gradient.
struct quad_term {
  real e = 0;
  vec3 ga0 = vec3::Zero();
  vec3 ga1 = vec3::Zero();
  vec3 gb0 = vec3::Zero();
  vec3 gb1 = vec3::Zero();
};

inline quad_term energy_quad_term(const vec3 &a0, const vec3 &a1, const vec3 &Ta,
                                  real ella, const vec3 &b0, const vec3 &b1,
                                  const vec3 &Tb, real ellb, real eps, real p,
                                  const edge_quad &q) {
  quad_term t;
  if (!(ella > real(0)) || !(ellb > real(0)) || q.t.empty())
    return t;
  real S = real(0);
  vec3 s0 = vec3::Zero();
  vec3 s1 = vec3::Zero();
  vec3 s2 = vec3::Zero();
  vec3 s3 = vec3::Zero();
  const size_t nq = q.t.size();
  for (size_t ia = 0; ia < nq; ++ia) {
    const real wa = ia < q.w.size() ? q.w[ia] : real(0);
    const real ta = q.t[ia];
    const vec3 pa = quad_point(a0, a1, ta);
    for (size_t ib = 0; ib < nq; ++ib) {
      const real wb = ib < q.w.size() ? q.w[ib] : real(0);
      const real w = wa * wb;
      if (!(w > real(0)))
        continue;
      const real tb = q.t[ib];
      const vec3 pb = quad_point(b0, b1, tb);
      const vec3 dp = pb - pa;
      const energy_dk d1 = energy_kernel_deriv(dp, Ta, eps, p);
      const energy_dk d2 = energy_kernel_deriv(-dp, Tb, eps, p);
      S += w * (d1.k + d2.k);

      const vec3 dp1 = w * d1.dk_dp;
      const vec3 dp2 = w * d2.dk_dp;
      s0 -= (real(1) - ta) * dp1;
      s1 -= ta * dp1;
      s2 += (real(1) - tb) * dp1;
      s3 += tb * dp1;
      s0 += (real(1) - ta) * dp2;
      s1 += ta * dp2;
      s2 -= (real(1) - tb) * dp2;
      s3 -= tb * dp2;

      const vec3 c1 = (d1.dk_dT - Ta * Ta.dot(d1.dk_dT)) / ella;
      const vec3 c2 = (d2.dk_dT - Tb * Tb.dot(d2.dk_dT)) / ellb;
      const vec3 wc1 = w * c1;
      const vec3 wc2 = w * c2;
      s0 -= wc1;
      s1 += wc1;
      s2 -= wc2;
      s3 += wc2;
    }
  }
  t.e = ella * ellb * S;
  const vec3 gla = Ta * (ellb * S);
  const vec3 glb = Tb * (ella * S);
  const real lab = ella * ellb;
  t.ga0 = lab * s0 - gla;
  t.ga1 = lab * s1 + gla;
  t.gb0 = lab * s2 - glb;
  t.gb1 = lab * s3 + glb;
  return t;
}

inline real energy_quad_weight(const vec3 &a0, const vec3 &a1, const vec3 &Ta,
                               real ella, const vec3 &b0, const vec3 &b1,
                               const vec3 &Tb, real ellb, real eps, real p,
                               const edge_quad &q) {
  return energy_quad_term(a0, a1, Ta, ella, b0, b1, Tb, ellb, eps, p, q).e;
}

/// Far lump at the query midpoint. Cluster center, tangent, and length are fixed.
struct far_term {
  real e = 0;
  vec3 ga0 = vec3::Zero();
  vec3 ga1 = vec3::Zero();
};

inline far_term energy_far_term(const vec3 &Ta, real ella, const vec3 &pi,
                                const vec3 &pj, const vec3 &Tj, real ellc,
                                real eps, real p) {
  far_term t;
  if (!(ella > real(0)) || !(ellc > real(0)))
    return t;
  const energy_dk d1 = energy_kernel_deriv(pj - pi, Ta, eps, p);
  const energy_dk d2 = energy_kernel_deriv(pi - pj, Tj, eps, p);
  const real S = d1.k + d2.k;
  t.e = ella * ellc * S;
  const vec3 gla = Ta * (ellc * S);
  const vec3 dpi = ella * ellc * real(0.5) * (d2.dk_dp - d1.dk_dp);
  const vec3 chord = (d1.dk_dT - Ta * Ta.dot(d1.dk_dT)) / ella;
  const vec3 gt = ella * ellc * chord;
  t.ga0 = dpi - gla - gt;
  t.ga1 = dpi + gla + gt;
  return t;
}

/// One accepted Barnes–Hut sample. `j` is the leaf id or the internal node id
/// the walker passed to the kernel. `w` does not depend on φ.
/// `align` is 1 for both B and B0 (Eq. 31). Kept so a stored row can still
/// record a coefficient on φ_J.
struct cov_weight_hit {
  index_t j = 0;
  real w = 0;
  bool far = false;
  real align = 1;
};

/// Per-query-edge samples for one frozen pose. `sum` is Σ w, the Jacobi row.
struct cov_weight_rows {
  std::vector<std::vector<cov_weight_hit>> row;
  std::vector<real> sum;
};

inline void bind_edge_lengths(const edge_table &E, std::vector<real> &ell_b,
                              std::vector<real> &ep0_b,
                              std::vector<real> &ep1_b) {
  const size_t n = E.n;
  ell_b.resize(n);
  ep0_b.resize(n);
  ep1_b.resize(n);
  for (size_t i = 0; i < n; ++i) {
    ell_b[i] = E.ell[i];
    ep0_b[i] = real(E.ep0[i]);
    ep1_b[i] = real(E.ep1[i]);
  }
}

/// Walk once and store w. Later matvecs pyramid φ and multiply.
inline cov_weight_rows
build_cov_weights(asawa::rod::rod &rod, const edge_table &E, real p,
                  real l0_abs, real bh_eps, real cone_eps,
                  bool use_tp_k, const arp::edge_cov_tree<2> &cov,
                  const edge_quad &quad = trapezoid_edge_quad()) {
  const size_t n = E.n;
  cov_weight_rows rows;
  rows.row.resize(n);
  rows.sum.assign(n, real(0));
  if (n == 0)
    return rows;

  std::vector<real> ell_b, ep0_b, ep1_b;
  bind_edge_lengths(E, ell_b, ep0_b, ep1_b);
  const std::vector<vec3> T_b = E.T;
  const std::vector<vec3> &x = rod.x();
  const real eps = std::max(l0_abs, real(1e-12));
  const std::vector<index_t> ids = edge_ordinals(n);

  auto bind = [&ell_b, &T_b, &ep0_b, &ep1_b](const std::vector<index_t> &edge_ids,
                                             calder::cov_rod_sum &sum) {
    sum.bind(calder::scalar_datum::create(edge_ids, ell_b));
    sum.bind(calder::vec3_datum::create(edge_ids, T_b));
    sum.bind(calder::scalar_datum::create(edge_ids, ep0_b));
    sum.bind(calder::scalar_datum::create(edge_ids, ep1_b));
  };
  auto kernel = [&E, &rows, &x, &quad, p, eps, use_tp_k](
                    const index_t i, const index_t j, const vec3 &pi,
                    const vec3 &pj, const std::vector<calder::datum::ptr> &data,
                    calder::Node_Type node_type,
                    const arp::edge_cov_tree<2> &) -> vec3 {
    const size_t ii = static_cast<size_t>(i);
    if (ii >= E.n)
      return vec3::Zero();
    const real lj = calder::get_data<real>(node_type, j, 0, data);
    vec3 Tj = calder::get_data<vec3>(node_type, j, 1, data);
    const int j0 = static_cast<int>(calder::get_data<real>(node_type, j, 2, data));
    const int j1 = static_cast<int>(calder::get_data<real>(node_type, j, 3, data));
    if (!(lj > real(0)))
      return vec3::Zero();
    const real w = sample_pair_weight(E, ii, pi, pj, Tj, lj, j0, j1, x,
                                       node_type, eps, p, use_tp_k,
                                       quad);
    if (!(w > real(0)) || !std::isfinite(w))
      return vec3::Zero();
    // Scalar kernel. ⟨T_I, T_J⟩ enters because B's edge unknown is (Δu/ℓ) T
    // (§5.2.1), not because the weight is scaled by the tangent dot.
    const real align = real(1);
    rows.row[ii].push_back(
        cov_weight_hit{j, w, node_type == calder::BRANCH, align});
    rows.sum[ii] += w;
    return vec3::Zero();
  };
  calder::integrate_over_cov_tree<vec3>(cov, ids, E.mid, bind, kernel, kernel,
                                        bh_eps, 1, cone_eps);
  return rows;
}

/// Sums of a frozen field on one level's covariance tree. Built once.
struct pyr_real {
  std::vector<real> leaf;
  std::vector<real> node;
};

struct pyr_vec3 {
  std::vector<vec3> leaf;
  std::vector<vec3> node;
};

struct pyr_mat3 {
  std::vector<mat3> leaf;
  std::vector<mat3> node;
};

/// Edge table plus the length pyramid. Both are constant for a frozen pose.
struct level_frozen {
  edge_table E;
  pyr_real ell;
  bool on = false;
};

template <typename Datum, typename T>
inline void pyramid_field(const arp::edge_cov_tree<2> &cov,
                          const std::vector<index_t> &ids,
                          const std::vector<T> &values, std::vector<T> &leaf,
                          std::vector<T> &node) {
  auto datum = Datum::create(ids, values);
  datum->pyramid(cov);
  leaf = datum->sorted_leaf_data();
  node = datum->node_data();
}

inline pyr_real pyramid_real(const arp::edge_cov_tree<2> &cov,
                             const std::vector<real> &values) {
  pyr_real out;
  pyramid_field<calder::scalar_datum>(cov, edge_ordinals(values.size()), values,
                                      out.leaf, out.node);
  return out;
}

inline pyr_vec3 pyramid_vec3(const arp::edge_cov_tree<2> &cov,
                             const std::vector<vec3> &values) {
  pyr_vec3 out;
  pyramid_field<calder::vec3_datum>(cov, edge_ordinals(values.size()), values,
                                    out.leaf, out.node);
  return out;
}

inline pyr_mat3 pyramid_mat3(const arp::edge_cov_tree<2> &cov,
                             const std::vector<mat3> &values) {
  pyr_mat3 out;
  pyramid_field<calder::mat3_datum>(cov, edge_ordinals(values.size()), values,
                                    out.leaf, out.node);
  return out;
}

inline const edge_table &level_edges(asawa::rod::rod &rod, const level_frozen *frozen,
                                   edge_table &local) {
  if (frozen && frozen->on)
    return frozen->E;
  local = build_edges(rod);
  return local;
}

inline const pyr_real &level_lengths(const arp::edge_cov_tree<2> &cov,
                                     const edge_table &E,
                                     const level_frozen *frozen,
                                     pyr_real &local) {
  if (frozen && frozen->on)
    return frozen->ell;
  local = pyramid_real(cov, E.ell);
  return local;
}

inline level_frozen freeze_level(asawa::rod::rod &rod,
                                 const arp::edge_cov_tree<2> &cov) {
  level_frozen f;
  f.E = build_edges(rod);
  f.ell = pyramid_real(cov, f.E.ell);
  f.on = true;
  return f;
}

/// r_I = Σ w (φ_I − align φ_J). align is 1.
/// A leaf uses φ_J. A far node uses Σ (ℓ φ) / L. `ell` is the cached length
/// pyramid. `mom` is Σ (ℓ φ) for this multiply.
inline std::vector<vec3>
apply_cov_weights(const edge_table &E, const std::vector<vec3> &phi,
                  const cov_weight_rows &rows, const pyr_real &ell,
                  const pyr_vec3 &mom) {
  const size_t n = E.n;
  std::vector<vec3> r(n, vec3::Zero());
  if (n == 0)
    return r;
  for (size_t i = 0; i < n && i < rows.row.size(); ++i) {
    const vec3 phi_i = i < phi.size() ? phi[i] : vec3::Zero();
    vec3 acc = vec3::Zero();
    for (const cov_weight_hit &hit : rows.row[i]) {
      const size_t j = static_cast<size_t>(hit.j);
      real L = real(0);
      vec3 M = vec3::Zero();
      if (hit.far) {
        if (j < ell.node.size())
          L = ell.node[j];
        if (j < mom.node.size())
          M = mom.node[j];
      } else {
        if (j < ell.leaf.size())
          L = ell.leaf[j];
        if (j < mom.leaf.size())
          M = mom.leaf[j];
      }
      if (!(L > real(0)))
        continue;
      acc += hit.w * phi_i - hit.align * (hit.w / L) * M;
    }
    r[i] = acc;
  }
  return r;
}

/// Edge residual r_I = Σ_J w_IJ (φ_I − φ_J). Opening query is the midpoint.
inline std::vector<vec3>
edge_K_diff(asawa::rod::rod &rod, const edge_table &E,
            const std::vector<vec3> &phi, real p, real l0_abs, real bh_eps,
            bool use_tp_k,
            const arp::edge_cov_tree<2> *cov = nullptr, real cone_eps = -1,
            const cov_weight_rows *weights = nullptr,
            const edge_quad &quad = trapezoid_edge_quad(),
            bool drop_adjacent = true) {
  const size_t n = E.n;
  std::vector<vec3> r(n, vec3::Zero());
  if (n == 0)
    return r;
  if (weights && cov) {
    std::vector<vec3> phi_b = phi;
    if (phi_b.size() < n)
      phi_b.resize(n, vec3::Zero());
    std::vector<vec3> moment(n, vec3::Zero());
    for (size_t i = 0; i < n; ++i)
      moment[i] = E.ell[i] * phi_b[i];
    const pyr_vec3 mom = pyramid_vec3(*cov, moment);
    const pyr_real ell = pyramid_real(*cov, E.ell);
    return apply_cov_weights(E, phi_b, *weights, ell, mom);
  }

  std::vector<real> ell_b, ep0_b, ep1_b;
  bind_edge_lengths(E, ell_b, ep0_b, ep1_b);
  const std::vector<vec3> T_b = E.T;
  std::vector<vec3> phi_b = phi;
  if (phi_b.size() < n)
    phi_b.resize(n, vec3::Zero());
  // Pyramid sums. Slot 2 is ℓφ so a far node is Σ(ℓ φ), and lj is L.
  std::vector<vec3> moment_b(n, vec3::Zero());
  for (size_t i = 0; i < n; ++i)
    moment_b[i] = ell_b[i] * phi_b[i];
  const std::vector<vec3> &x = rod.x();
  const real eps = std::max(l0_abs, real(1e-12));
  const std::vector<index_t> ids = edge_ordinals(n);

  auto bind = [&ell_b, &T_b, &moment_b, &ep0_b, &ep1_b](
                  const std::vector<index_t> &edge_ids, calder::cov_rod_sum &sum) {
    sum.bind(calder::scalar_datum::create(edge_ids, ell_b));
    sum.bind(calder::vec3_datum::create(edge_ids, T_b));
    sum.bind(calder::vec3_datum::create(edge_ids, moment_b));
    sum.bind(calder::scalar_datum::create(edge_ids, ep0_b));
    sum.bind(calder::scalar_datum::create(edge_ids, ep1_b));
  };
  auto kernel = [&E, &phi_b, &x, &quad, p, eps, use_tp_k, drop_adjacent](
                    const index_t i, const index_t j, const vec3 &pi,
                    const vec3 &pj, const std::vector<calder::datum::ptr> &data,
                    calder::Node_Type node_type,
                    const arp::edge_cov_tree<2> &) -> vec3 {
        const size_t ii = static_cast<size_t>(i);
        if (ii >= E.n)
          return vec3::Zero();

        const real lj = calder::get_data<real>(node_type, j, 0, data);
        vec3 Tj = calder::get_data<vec3>(node_type, j, 1, data);
        const vec3 moment_j = calder::get_data<vec3>(node_type, j, 2, data);
        const int j0 =
            static_cast<int>(calder::get_data<real>(node_type, j, 3, data));
        const int j1 =
            static_cast<int>(calder::get_data<real>(node_type, j, 4, data));
        if (!(lj > real(0)))
          return vec3::Zero();
        const real w = sample_pair_weight(E, ii, pi, pj, Tj, lj, j0, j1, x,
                                           node_type, eps, p,
                                           use_tp_k, quad, drop_adjacent);
        if (!(w > real(0)) || !std::isfinite(w))
          return vec3::Zero();
        const real align = real(1);
        // Leaf: moment/ℓ = φ. Far: moment/L = Σ(ℓ φ)/L.
        return w * phi_b[ii] - align * (w / lj) * moment_j;
      };
  arp::edge_cov_tree<2>::ptr owned;
  const arp::edge_cov_tree<2> *tree = cov;
  if (!tree) {
    owned = calder::build_rod_cov_tree(rod);
    tree = owned.get();
  }
  const real cone = calder::cov_cone_eps(bh_eps, cone_eps);
  return calder::integrate_over_cov_tree<vec3>(*tree, ids, E.mid, bind, kernel,
                                               kernel, bh_eps, 1, cone);
}

inline std::vector<real>
edge_K_row_sum(asawa::rod::rod &rod, const edge_table &E, real p, real l0_abs,
               real bh_eps, bool use_tp_k,
               const arp::edge_cov_tree<2> *cov = nullptr, real cone_eps = -1,
               const cov_weight_rows *weights = nullptr,
               const edge_quad &quad = trapezoid_edge_quad()) {
  const size_t n = E.n;
  std::vector<real> S(n, real(0));
  if (n == 0)
    return S;
  if (weights) {
    for (size_t i = 0; i < n && i < weights->sum.size(); ++i)
      S[i] = weights->sum[i];
    return S;
  }

  std::vector<real> ell_b, ep0_b, ep1_b;
  bind_edge_lengths(E, ell_b, ep0_b, ep1_b);
  const std::vector<vec3> T_b = E.T;
  const std::vector<vec3> &x = rod.x();
  const real eps = std::max(l0_abs, real(1e-12));
  const std::vector<index_t> ids = edge_ordinals(n);

  auto bind = [&ell_b, &T_b, &ep0_b, &ep1_b](const std::vector<index_t> &edge_ids,
                                             calder::cov_rod_sum &sum) {
    sum.bind(calder::scalar_datum::create(edge_ids, ell_b));
    sum.bind(calder::vec3_datum::create(edge_ids, T_b));
    sum.bind(calder::scalar_datum::create(edge_ids, ep0_b));
    sum.bind(calder::scalar_datum::create(edge_ids, ep1_b));
  };
  auto kernel = [&E, &x, &quad, p, eps, use_tp_k](
                    const index_t i, const index_t j, const vec3 &pi,
                    const vec3 &pj, const std::vector<calder::datum::ptr> &data,
                    calder::Node_Type node_type,
                    const arp::edge_cov_tree<2> &) -> vec3 {
        const size_t ii = static_cast<size_t>(i);
        if (ii >= E.n)
          return vec3::Zero();
        const real lj = calder::get_data<real>(node_type, j, 0, data);
        vec3 Tj = calder::get_data<vec3>(node_type, j, 1, data);
        const int j0 =
            static_cast<int>(calder::get_data<real>(node_type, j, 2, data));
        const int j1 =
            static_cast<int>(calder::get_data<real>(node_type, j, 3, data));
        if (!(lj > real(0)))
          return vec3::Zero();
        const real w = sample_pair_weight(E, ii, pi, pj, Tj, lj, j0, j1, x,
                                           node_type, eps, p,
                                           use_tp_k, quad);
        if (!(w > real(0)) || !std::isfinite(w))
          return vec3::Zero();
        return vec3(w, 0, 0);
      };
  arp::edge_cov_tree<2>::ptr owned;
  const arp::edge_cov_tree<2> *tree = cov;
  if (!tree) {
    owned = calder::build_rod_cov_tree(rod);
    tree = owned.get();
  }
  const real cone = calder::cov_cone_eps(bh_eps, cone_eps);
  const std::vector<vec3> packed = calder::integrate_over_cov_tree<vec3>(
      *tree, ids, E.mid, bind, kernel, kernel, bh_eps, 1, cone);
  for (size_t i = 0; i < n; ++i)
    S[i] = packed[i][0];
  return S;
}

inline void matvec_B0(asawa::rod::rod &rod, const std::vector<vec3> &u, real p,
                      real l0_abs, real bh_eps, std::vector<vec3> &out,
                      const arp::edge_cov_tree<2> *cov = nullptr,
                      real cone_eps = -1,
                      const cov_weight_rows *weights = nullptr,
                      bool drop_adjacent = true,
                      const level_frozen *frozen = nullptr) {
  edge_table local_E;
  const edge_table &E = level_edges(rod, frozen, local_E);
  const std::vector<vec3> Eu = edge_avg_E(E, u);
  if (weights && cov) {
    std::vector<vec3> moment(E.n, vec3::Zero());
    for (size_t i = 0; i < E.n; ++i)
      moment[i] = E.ell[i] * Eu[i];
    const pyr_vec3 mom = pyramid_vec3(*cov, moment);
    pyr_real local_ell;
    const pyr_real &ell = level_lengths(*cov, E, frozen, local_ell);
    const std::vector<vec3> r =
        apply_cov_weights(E, Eu, *weights, ell, mom);
    scatter_ET(E, r, out);
    return;
  }
  const std::vector<vec3> r =
      edge_K_diff(rod, E, Eu, p, l0_abs, bh_eps, /*use_tp_k=*/true,
                  cov, cone_eps, weights, trapezoid_edge_quad(), drop_adjacent);
  scatter_ET(E, r, out);
}

inline void matvec_B(asawa::rod::rod &rod, const std::vector<vec3> &u, real p,
                     real l0_abs, real bh_eps, std::vector<vec3> &out,
                     const arp::edge_cov_tree<2> *cov = nullptr,
                     real cone_eps = -1,
                     const cov_weight_rows *weights = nullptr,
                     bool drop_adjacent = true,
                     const level_frozen *frozen = nullptr) {
  edge_table local_E;
  const edge_table &E = level_edges(rod, frozen, local_E);
  const std::vector<vec3> slope = edge_diff_D(E, u);
  out.assign(E.n_vert, vec3::Zero());
  // §5.2.1–5.2.2: each coordinate is its own scalar, lifted by T. Eq. 19
  // copies that scalar operator onto x, y, and z. One mat3 pyramid carries
  // all three moments; the length pyramid is the cached one.
  if (weights && cov) {
    std::vector<mat3> moment(E.n, mat3::Zero());
    for (size_t i = 0; i < E.n; ++i) {
      for (int c = 0; c < 3; ++c)
        moment[i].col(c) = slope[i][c] * E.ell[i] * E.T[i];
    }
    const pyr_mat3 mom = pyramid_mat3(*cov, moment);
    pyr_real local_ell;
    const pyr_real &ell = level_lengths(*cov, E, frozen, local_ell);
    for (size_t i = 0; i < E.n && i < weights->row.size(); ++i) {
      for (int c = 0; c < 3; ++c) {
        const vec3 phi = slope[i][c] * E.T[i];
        vec3 acc = vec3::Zero();
        for (const cov_weight_hit &hit : weights->row[i]) {
          const size_t j = static_cast<size_t>(hit.j);
          real L = real(0);
          vec3 M = vec3::Zero();
          if (hit.far) {
            if (j < ell.node.size())
              L = ell.node[j];
            if (j < mom.node.size())
              M = mom.node[j].col(c);
          } else {
            if (j < ell.leaf.size())
              L = ell.leaf[j];
            if (j < mom.leaf.size())
              M = mom.leaf[j].col(c);
          }
          if (!(L > real(0)))
            continue;
          acc += hit.w * phi - hit.align * (hit.w / L) * M;
        }
        const real s = E.T[i].dot(acc) / E.ell[i];
        out[static_cast<size_t>(E.ep0[i])][c] -= s;
        out[static_cast<size_t>(E.ep1[i])][c] += s;
      }
    }
    return;
  }
  for (int c = 0; c < 3; ++c) {
    std::vector<vec3> phi(E.n, vec3::Zero());
    for (size_t i = 0; i < E.n; ++i)
      phi[i] = slope[i][c] * E.T[i];
    const std::vector<vec3> r =
        edge_K_diff(rod, E, phi, p, l0_abs, bh_eps, /*use_tp_k=*/false, cov,
                    cone_eps, weights, trapezoid_edge_quad(), drop_adjacent);
    for (size_t i = 0; i < E.n; ++i) {
      const real s = E.T[i].dot(r[i]) / E.ell[i];
      out[static_cast<size_t>(E.ep0[i])][c] -= s;
      out[static_cast<size_t>(E.ep1[i])][c] += s;
    }
  }
}

inline void matvec_A(asawa::rod::rod &rod, real alpha_B, real alpha_0, real p,
                     real l0_abs, real bh_eps, const std::vector<vec3> &u,
                     std::vector<vec3> &out,
                     const arp::edge_cov_tree<2> *cov = nullptr,
                     real cone_eps = -1, const cov_weight_rows *wB = nullptr,
                     const cov_weight_rows *wB0 = nullptr,
                     bool drop_adjacent = true,
                     const level_frozen *frozen = nullptr) {
  const size_t n = u.size();
  out.assign(n, vec3::Zero());

  if (alpha_0 != real(0)) {
    std::vector<vec3> B0u;
    matvec_B0(rod, u, p, l0_abs, bh_eps, B0u, cov, cone_eps, wB0,
              drop_adjacent, frozen);
    for (size_t i = 0; i < n; ++i)
      out[i] += alpha_0 * B0u[i];
  }
  if (alpha_B != real(0)) {
    std::vector<vec3> Bu;
    matvec_B(rod, u, p, l0_abs, bh_eps, Bu, cov, cone_eps, wB, drop_adjacent,
             frozen);
    for (size_t i = 0; i < n; ++i)
      out[i] += alpha_B * Bu[i];
  }
}

/// All-pairs leaf kernel. Same midpoint weight as a leaves-only tree walk,
/// without building or walking a covariance tree.
inline std::vector<vec3>
edge_K_diff_dense(const edge_table &E, const std::vector<vec3> &phi,
                  const std::vector<vec3> &x, real p, real l0_abs,
                  bool use_tp_k, bool drop_adjacent = true) {
  const size_t n = E.n;
  std::vector<vec3> r(n, vec3::Zero());
  if (n == 0)
    return r;
  const real eps = std::max(l0_abs, real(1e-12));
  const edge_quad quad = trapezoid_edge_quad();
  for (size_t i = 0; i < n; ++i) {
    const vec3 phi_i = i < phi.size() ? phi[i] : vec3::Zero();
    vec3 acc = vec3::Zero();
    for (size_t j = 0; j < n; ++j) {
      const real w = sample_pair_weight(
          E, i, E.mid[i], E.mid[j], E.T[j], E.ell[j], E.ep0[j], E.ep1[j], x,
          calder::LEAF, eps, p, use_tp_k, quad, drop_adjacent);
      if (!(w > real(0)) || !std::isfinite(w))
        continue;
      const vec3 phi_j = j < phi.size() ? phi[j] : vec3::Zero();
      acc += w * (phi_i - phi_j);
    }
    r[i] = acc;
  }
  return r;
}

inline void matvec_B0_dense(asawa::rod::rod &rod, const std::vector<vec3> &u,
                            real p, real l0_abs, std::vector<vec3> &out,
                            bool drop_adjacent = true) {
  const edge_table E = build_edges(rod);
  const std::vector<vec3> Eu = edge_avg_E(E, u);
  const std::vector<vec3> r = edge_K_diff_dense(
      E, Eu, rod.x(), p, l0_abs, /*use_tp_k=*/true, drop_adjacent);
  scatter_ET(E, r, out);
}

inline void matvec_B_dense(asawa::rod::rod &rod, const std::vector<vec3> &u,
                           real p, real l0_abs, std::vector<vec3> &out,
                           bool drop_adjacent = true) {
  const edge_table E = build_edges(rod);
  const std::vector<vec3> slope = edge_diff_D(E, u);
  out.assign(E.n_vert, vec3::Zero());
  for (int c = 0; c < 3; ++c) {
    std::vector<vec3> phi(E.n, vec3::Zero());
    for (size_t i = 0; i < E.n; ++i)
      phi[i] = slope[i][c] * E.T[i];
    const std::vector<vec3> r = edge_K_diff_dense(
        E, phi, rod.x(), p, l0_abs, /*use_tp_k=*/false, drop_adjacent);
    for (size_t i = 0; i < E.n; ++i) {
      const real s = E.T[i].dot(r[i]) / E.ell[i];
      out[static_cast<size_t>(E.ep0[i])][c] -= s;
      out[static_cast<size_t>(E.ep1[i])][c] += s;
    }
  }
}

inline void matvec_A_dense(asawa::rod::rod &rod, real alpha_B, real alpha_0,
                           real p, real l0_abs, const std::vector<vec3> &u,
                           std::vector<vec3> &out, bool drop_adjacent = true) {
  const size_t n = u.size();
  out.assign(n, vec3::Zero());
  if (alpha_0 != real(0)) {
    std::vector<vec3> B0u;
    matvec_B0_dense(rod, u, p, l0_abs, B0u, drop_adjacent);
    for (size_t i = 0; i < n; ++i)
      out[i] += alpha_0 * B0u[i];
  }
  if (alpha_B != real(0)) {
    std::vector<vec3> Bu;
    matvec_B_dense(rod, u, p, l0_abs, Bu, drop_adjacent);
    for (size_t i = 0; i < n; ++i)
      out[i] += alpha_B * Bu[i];
  }
}

inline std::vector<real> jacobi_diag(
    asawa::rod::rod &rod, real alpha_B, real alpha_0, real p, real l0_abs,
    real bh_eps, const arp::edge_cov_tree<2> *cov = nullptr,
    real cone_eps = -1, const cov_weight_rows *wB = nullptr,
    const cov_weight_rows *wB0 = nullptr) {
  const size_t n = rod.corner_count();
  std::vector<real> d(n, real(0));

  const edge_table E = build_edges(rod);

  if (alpha_0 != real(0)) {
    const auto S0 =
        edge_K_row_sum(rod, E, p, l0_abs, bh_eps, /*use_tp_k=*/true,
                       cov, cone_eps, wB0);
    for (size_t i = 0; i < E.n; ++i) {
      // E and Eᵀ each contribute 1/2, so the vertex diagonal is a quarter of the
      // edge row sum (Eq. 31).
      const real lump = real(0.25) * alpha_0 * S0[i];
      d[static_cast<size_t>(E.ep0[i])] += lump;
      d[static_cast<size_t>(E.ep1[i])] += lump;
    }
  }
  if (alpha_B != real(0)) {
    const auto SB = edge_K_row_sum(rod, E, p, l0_abs, bh_eps,
                                   /*use_tp_k=*/false, cov, cone_eps, wB);
    for (size_t i = 0; i < E.n; ++i) {
      const real inv_l2 = real(1) / (E.ell[i] * E.ell[i]);
      const real lump = alpha_B * SB[i] * inv_l2;
      d[static_cast<size_t>(E.ep0[i])] += lump;
      d[static_cast<size_t>(E.ep1[i])] += lump;
    }
  }
  for (real &x : d)
    x = std::max(x, real(1e-12));
  return d;
}

struct paper_tp_eval {
  real energy = 0;
  std::vector<vec3> grad;
};

/// One vertex of the reference tangent-point energy. `T` is the unit sum of
/// the oriented incident edge tangents. `dual` is the average incident length.
struct tpe_vertex {
  vec3 T = vec3::Zero();
  real dual = 0;
  int prev = -1;
  int next = -1;
  int deg = 0;
};

inline std::vector<tpe_vertex> tpe_vertices(const asawa::rod::rod &rod) {
  const auto &x = rod.x();
  const size_t n = x.size();
  std::vector<tpe_vertex> V(n);
  for (size_t i = 0; i < n; ++i) {
    const int ip = static_cast<int>(rod.prev(asawa::rod::corner_id(static_cast<int>(i))));
    const int in = static_cast<int>(rod.next(asawa::rod::corner_id(static_cast<int>(i))));
    tpe_vertex v;
    v.prev = ip;
    v.next = in;
    vec3 sum = vec3::Zero();
    real len = 0;
    if (ip >= 0 && static_cast<size_t>(ip) < n) {
      const vec3 e = x[i] - x[static_cast<size_t>(ip)];
      const real L = e.norm();
      if (L > real(1e-12)) {
        sum += e / L;
        len += L;
        ++v.deg;
      }
    }
    if (in >= 0 && static_cast<size_t>(in) < n && in != ip) {
      const vec3 e = x[static_cast<size_t>(in)] - x[i];
      const real L = e.norm();
      if (L > real(1e-12)) {
        sum += e / L;
        len += L;
        ++v.deg;
      }
    }
    if (v.deg > 0) {
      v.dual = len / real(v.deg);
      if (sum.squaredNorm() > real(1e-24))
        v.T = sum.normalized();
    }
    V[i] = v;
  }
  return V;
}

/// tpe_Kf_pts. `d` is x_i − x_j. Tangent is unit. No ε.
inline real tpe_kf(const vec3 &d, const vec3 &T, real alpha, real beta) {
  if (d.squaredNorm() < real(1e-24) || T.squaredNorm() < real(1e-24))
    return real(0);
  const vec3 n = d - T * T.dot(d);
  const real rn = n.norm();
  const real r = d.norm();
  if (!(rn > real(0)) || !(r > real(0)))
    return real(0);
  return std::pow(rn, alpha) / std::pow(r, beta);
}

struct energy_scatter {
  real e = 0;
  int n = 0;
  int v[4] = {-1, -1, -1, -1};
  vec3 g0 = vec3::Zero();
  vec3 g1 = vec3::Zero();
  vec3 g2 = vec3::Zero();
  vec3 g3 = vec3::Zero();
};

using tpe_mat3 = Eigen::Matrix<real, 3, 3>;

/// Unit edge tangent (head − tail) / ℓ, differentiated wrt the tail (0) or head (1).
inline tpe_mat3 tpe_dT_edge(const vec3 &tail, const vec3 &head, int wrt_head) {
  const vec3 e = head - tail;
  const real L = e.norm();
  tpe_mat3 Z = tpe_mat3::Zero();
  if (!(L > real(1e-12)))
    return Z;
  const vec3 t = e / L;
  tpe_mat3 dT = (tpe_mat3::Identity() - t * t.transpose()) / L;
  if (!wrt_head)
    dT = -dT;
  return dT;
}

/// Derivative of the vertex tangent at `i` with respect to corner `wrt`.
inline tpe_mat3 tpe_dT(const std::vector<vec3> &x, const tpe_vertex &v, int i,
                       int wrt) {
  tpe_mat3 dsum = tpe_mat3::Zero();
  vec3 sum = vec3::Zero();
  if (v.prev >= 0) {
    const vec3 e = x[static_cast<size_t>(i)] - x[static_cast<size_t>(v.prev)];
    if (e.norm() > real(1e-12)) {
      sum += e.normalized();
      if (wrt == i)
        dsum += tpe_dT_edge(x[static_cast<size_t>(v.prev)],
                            x[static_cast<size_t>(i)], 1);
      if (wrt == v.prev)
        dsum += tpe_dT_edge(x[static_cast<size_t>(v.prev)],
                            x[static_cast<size_t>(i)], 0);
    }
  }
  if (v.next >= 0 && v.next != v.prev) {
    const vec3 e = x[static_cast<size_t>(v.next)] - x[static_cast<size_t>(i)];
    if (e.norm() > real(1e-12)) {
      sum += e.normalized();
      if (wrt == v.next)
        dsum += tpe_dT_edge(x[static_cast<size_t>(i)],
                            x[static_cast<size_t>(v.next)], 1);
      if (wrt == i)
        dsum += tpe_dT_edge(x[static_cast<size_t>(i)],
                            x[static_cast<size_t>(v.next)], 0);
    }
  }
  const real ns = sum.norm();
  if (!(ns > real(1e-12)))
    return tpe_mat3::Zero();
  const vec3 u = sum / ns;
  const tpe_mat3 P = (tpe_mat3::Identity() - u * u.transpose()) / ns;
  return P * dsum;
}

/// Derivative of dual length at `i` with respect to corner `wrt`.
inline vec3 tpe_d_dual(const std::vector<vec3> &x, const tpe_vertex &v, int i,
                       int wrt) {
  if (v.deg <= 0)
    return vec3::Zero();
  vec3 g = vec3::Zero();
  if (v.prev >= 0) {
    const vec3 e = x[static_cast<size_t>(i)] - x[static_cast<size_t>(v.prev)];
    const real L = e.norm();
    if (L > real(1e-12)) {
      const vec3 t = e / L;
      if (wrt == i)
        g += t;
      if (wrt == v.prev)
        g -= t;
    }
  }
  if (v.next >= 0 && v.next != v.prev) {
    const vec3 e = x[static_cast<size_t>(v.next)] - x[static_cast<size_t>(i)];
    const real L = e.norm();
    if (L > real(1e-12)) {
      const vec3 t = e / L;
      if (wrt == i)
        g -= t;
      if (wrt == v.next)
        g += t;
    }
  }
  return g / real(v.deg);
}

/// One ordered pair E(i, j) = dual_i * m_j * k(T_i, x_i − p_j).
/// `jcorner < 0` freezes the partner (a Barnes–Hut cluster).
inline void tpe_accum_pair(const std::vector<vec3> &x,
                           const std::vector<tpe_vertex> &V, int i,
                           int jcorner, const vec3 &pj, real mj, real alpha,
                           real beta, real &energy, std::vector<vec3> &grad) {
  if (i < 0 || static_cast<size_t>(i) >= V.size())
    return;
  const tpe_vertex &vi = V[static_cast<size_t>(i)];
  if (vi.deg <= 0 || vi.deg > 2 || vi.T.squaredNorm() < real(1e-24))
    return;
  if (!(vi.dual > real(0)) || !(mj > real(0)))
    return;
  const vec3 d = x[static_cast<size_t>(i)] - pj;
  const real r = d.norm();
  const vec3 n = d - vi.T * vi.T.dot(d);
  const real rn = n.norm();
  if (!(r > real(1e-12)) || !(rn > real(1e-14)))
    return;
  const real A = std::pow(rn, alpha);
  const real B = std::pow(r, beta);
  const real k = A / B;
  const real li = vi.dual;
  energy += li * mj * k;

  const real rn_am2 = std::pow(rn, alpha - real(2));
  const real r_bm2 = std::pow(r, beta - real(2));
  const vec3 dA_dd = alpha * rn_am2 * n;
  const vec3 dB_dd = beta * r_bm2 * d;
  const vec3 dk_dd = (dA_dd * B - A * dB_dd) / (B * B);
  const vec3 gT = (-alpha * rn_am2 * d.dot(vi.T) * n) / B;

  int touch[6] = {i, vi.prev, vi.next, -1, -1, -1};
  int ntouch = 3;
  if (jcorner >= 0) {
    const tpe_vertex &vj = V[static_cast<size_t>(jcorner)];
    touch[3] = jcorner;
    touch[4] = vj.prev;
    touch[5] = vj.next;
    ntouch = 6;
  }
  for (int t = 0; t < ntouch; ++t) {
    const int w = touch[t];
    if (w < 0 || static_cast<size_t>(w) >= grad.size())
      continue;
    bool seen = false;
    for (int s = 0; s < t; ++s)
      if (touch[s] == w)
        seen = true;
    if (seen)
      continue;
    vec3 dk = tpe_dT(x, vi, i, w).transpose() * gT;
    if (w == i)
      dk += dk_dd;
    if (jcorner >= 0 && w == jcorner)
      dk -= dk_dd;
    const vec3 dli = tpe_d_dual(x, vi, i, w);
    vec3 dlj = vec3::Zero();
    if (jcorner >= 0)
      dlj = tpe_d_dual(x, V[static_cast<size_t>(jcorner)], jcorner, w);
    const vec3 ge = dli * (mj * k) + dlj * (li * k) + (li * mj) * dk;
    if (ge.allFinite())
      grad[static_cast<size_t>(w)] += ge;
  }
}

/// Vertex tangent-point energy. The visitor is Calder's Barnes–Hut walk on
/// vertex positions. bh_eps ≤ 0 never opens, so every partner is a vertex.
/// A far node is one partner at the dual-length center of mass, with mass
/// equal to the sum of dual lengths. The tangent in the kernel is the query
/// vertex tangent, the same pair the leaf evaluates.
inline paper_tp_eval eval_paper_tp(asawa::rod::rod &rod, real alpha, real beta,
                                   real bh_eps) {
  paper_tp_eval out;
  const std::vector<vec3> &x = rod.x();
  const size_t n = x.size();
  out.grad.assign(n, vec3::Zero());
  if (n < 2)
    return out;
  const std::vector<tpe_vertex> V = tpe_vertices(rod);

  std::vector<int> live;
  live.reserve(n);
  for (size_t i = 0; i < n; ++i)
    if (V[i].deg > 0)
      live.push_back(static_cast<int>(i));
  const size_t m = live.size();
  if (m < 2)
    return out;

  std::vector<vec3> pos(m);
  std::vector<real> mass(m);
  std::vector<vec3> mx(m), mT(m);
  std::vector<index_t> ids(m), adj(m);
  for (size_t k = 0; k < m; ++k) {
    const int c = live[k];
    pos[k] = x[static_cast<size_t>(c)];
    mass[k] = V[static_cast<size_t>(c)].dual;
    mx[k] = mass[k] * pos[k];
    mT[k] = mass[k] * V[static_cast<size_t>(c)].T;
    ids[k] = static_cast<index_t>(k);
    adj[k] = static_cast<index_t>(k);
  }

  arp::T1::ptr tree = arp::T1::create(adj, pos, 12);
  if (!tree || tree->internal_nodes_.empty())
    return out;

  calder::fast_summation<arp::T1> sum(*tree);
  sum.bind(calder::scalar_datum::create(ids, mass));
  sum.bind(calder::vec3_datum::create(ids, mx));
  sum.bind(calder::vec3_datum::create(ids, mT));
  for (auto &d : sum.__data)
    d->pyramid(*tree);

  const auto &data = sum.__data;
  real energy = 0;
#pragma omp parallel
  {
    std::vector<vec3> lg(n, vec3::Zero());
    real le = 0;
#pragma omp for
    for (int k = 0; k < static_cast<int>(m); ++k) {
      const int ic = live[static_cast<size_t>(k)];
      calder::traverse_bh_opening(
          *tree, static_cast<index_t>(k), pos[static_cast<size_t>(k)], bh_eps,
          [&](index_t node_id, index_t, const vec3 &) {
            const real M =
                calder::get_data<real>(calder::BRANCH, node_id, 0, data);
            if (!(M > real(0)))
              return;
            const vec3 COM =
                calder::get_data<vec3>(calder::BRANCH, node_id, 1, data) / M;
            tpe_accum_pair(x, V, ic, -1, COM, M, alpha, beta, le, lg);
          },
          [&](index_t leaf_id, index_t, index_t, const vec3 &) {
            const index_t orig = tree->get_index(leaf_id);
            if (orig < 0 || static_cast<size_t>(orig) >= m)
              return;
            const int jc = live[static_cast<size_t>(orig)];
            if (jc == ic)
              return;
            tpe_accum_pair(x, V, ic, jc, x[static_cast<size_t>(jc)],
                           V[static_cast<size_t>(jc)].dual, alpha, beta, le,
                           lg);
          });
    }
#pragma omp critical
    {
      energy += le;
      for (size_t i = 0; i < n; ++i)
        out.grad[i] += lg[i];
    }
  }
  out.energy = energy;
  return out;
}

} // namespace sobolev_detail

/// P = I - Bᵀ (BBᵀ)⁻¹ B. BBᵀ is factored once with Eigen::SparseLU.
/// Reference: NullSpaceProjector in repulsive-curves.
struct sobolev_projector {
  struct term {
    int i = -1;
    vec3 c = vec3::Zero();
  };
  std::vector<std::vector<term>> rows;
  /// Ridge on each row of BBᵀ. 0 is an exact nullspace row.
  std::vector<real> ridge;
  /// Ridge for the rows that survived factor(), aligned with B.
  std::vector<real> ridge_kept;
  Eigen::SparseMatrix<double> B;
  Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
      bbt_lu;
  bool on = false;

  void clear() {
    rows.clear();
    ridge.clear();
    ridge_kept.clear();
    B.resize(0, 0);
    on = false;
  }

  void factor() {
    on = false;
    B.resize(0, 0);
    const int m = static_cast<int>(rows.size());
    if (m == 0)
      return;
    int n = 0;
    for (const auto &row : rows)
      for (const term &t : row)
        n = std::max(n, t.i + 1);
    if (n <= 0)
      return;
    std::vector<Eigen::Triplet<double>> trips;
    trips.reserve(static_cast<size_t>(m) * 6);
    std::vector<int> kept;
    kept.reserve(static_cast<size_t>(m));
    for (int a = 0; a < m; ++a) {
      const size_t before = trips.size();
      for (const term &t : rows[static_cast<size_t>(a)]) {
        if (t.i < 0 || t.i >= n)
          continue;
        for (int c = 0; c < 3; ++c) {
          if (t.c[c] == real(0))
            continue;
          trips.emplace_back(static_cast<int>(kept.size()), 3 * t.i + c,
                             static_cast<double>(t.c[c]));
        }
      }
      if (trips.size() > before)
        kept.push_back(a);
      else
        trips.resize(before);
    }
    if (trips.empty())
      return;
    const int m_kept = static_cast<int>(kept.size());
    B.resize(m_kept, 3 * n);
    B.setFromTriplets(trips.begin(), trips.end());
    const Eigen::SparseMatrix<double> BBT = B * B.transpose();
    bbt_lu.compute(BBT);
    on = bbt_lu.info() == Eigen::Success;
  }

  void apply(std::vector<vec3> &u) const {
    if (!on || u.empty() || B.rows() == 0)
      return;
    const int nB = static_cast<int>(B.cols() / 3);
    const int n = static_cast<int>(
        std::min(u.size(), static_cast<size_t>(std::max(nB, 0))));
    if (n <= 0)
      return;
    Eigen::VectorXd v(3 * n);
    for (int i = 0; i < n; ++i) {
      v(3 * i) = u[static_cast<size_t>(i)][0];
      v(3 * i + 1) = u[static_cast<size_t>(i)][1];
      v(3 * i + 2) = u[static_cast<size_t>(i)][2];
    }
    const Eigen::VectorXd Bv = Eigen::VectorXd(B * v);
    const Eigen::VectorXd z = bbt_lu.solve(Bv);
    const Eigen::VectorXd corr = B.transpose() * z;
    for (int i = 0; i < n; ++i) {
      u[static_cast<size_t>(i)][0] -= corr(3 * i);
      u[static_cast<size_t>(i)][1] -= corr(3 * i + 1);
      u[static_cast<size_t>(i)][2] -= corr(3 * i + 2);
    }
  }
};

/// Least-squares block Cᵀ W C. Each stored coefficient already includes w,
/// and `phi` already includes w, so the products are w².
struct sobolev_gram {
  struct term {
    int i = -1;
    vec3 c = vec3::Zero();
  };
  struct row {
    real phi = 0;
    std::vector<term> terms;
  };
  std::vector<row> rows;

  void clear() { rows.clear(); }
  bool empty() const { return rows.empty(); }

  void apply(const std::vector<vec3> &u, std::vector<vec3> &out) const {
    if (rows.empty() || u.empty())
      return;
    if (out.size() < u.size())
      out.resize(u.size(), vec3::Zero());
    for (const row &r : rows) {
      real s = 0;
      for (const term &t : r.terms) {
        if (t.i < 0 || static_cast<size_t>(t.i) >= u.size())
          continue;
        s += t.c.dot(u[static_cast<size_t>(t.i)]);
      }
      for (const term &t : r.terms) {
        if (t.i < 0 || static_cast<size_t>(t.i) >= out.size())
          continue;
        out[static_cast<size_t>(t.i)] += s * t.c;
      }
    }
  }

  /// Cᵀ W Φ, the gradient of ½ ‖W¹ᐟ² Φ‖².
  void add_force(std::vector<vec3> &out) const {
    for (const row &r : rows) {
      if (r.phi == real(0))
        continue;
      for (const term &t : r.terms) {
        if (t.i < 0 || static_cast<size_t>(t.i) >= out.size())
          continue;
        out[static_cast<size_t>(t.i)] += r.phi * t.c;
      }
    }
  }

  void add_diag(std::vector<real> &d) const {
    for (const row &r : rows) {
      for (const term &t : r.terms) {
        if (t.i < 0 || static_cast<size_t>(t.i) >= d.size())
          continue;
        d[static_cast<size_t>(t.i)] += t.c.squaredNorm();
      }
    }
  }
};

/// Per-level hard projector and projective Gram. Null build leaves both off.
struct sobolev_mg_hook {
  void (*build)(asawa::rod::rod &rod, sobolev_projector &P, sobolev_gram &G,
                void *ctx) = nullptr;
  void *ctx = nullptr;
};

/// Jacobi smooth Ag≈f for n_iters (in-place on g).
inline void sobolev_jacobi_smooth(
    asawa::rod::rod &rod, const std::vector<vec3> &f,
    const std::vector<real> &diag, const sobolev_tangent_point_force_config &cfg,
    real l0_abs, int n_iters, std::vector<vec3> &g,
    const arp::edge_cov_tree<2> *cov = nullptr,
    const sobolev_detail::cov_weight_rows *wB = nullptr,
    const sobolev_detail::cov_weight_rows *wB0 = nullptr) {
  const size_t n = f.size();
  if (n == 0 || n_iters <= 0)
    return;
  const real omega = real(0.85);
  const real bh = cfg.bh_eps;
  for (int it = 0; it < n_iters; ++it) {
    std::vector<vec3> Ag;
    sobolev_detail::matvec_A(rod, cfg.alpha_B, cfg.alpha_0, cfg.p1, l0_abs, bh,
                             g, Ag, cov, cfg.cone_eps, wB, wB0);
    for (size_t i = 0; i < n; ++i)
      g[i] += omega * (f[i] - Ag[i]) / diag[i];
  }
}

/// Weighted Jacobi on Ag≈b. `apply_A` is the same projected matvec CG uses.
template <typename ApplyA>
inline void sobolev_jacobi_relax(const std::vector<vec3> &b,
                                 const std::vector<real> &diag, int n_iters,
                                 std::vector<vec3> &x, ApplyA &&apply_A) {
  const size_t n = b.size();
  if (n == 0 || n_iters <= 0)
    return;
  if (x.size() != n)
    x.assign(n, vec3::Zero());
  const real omega = real(0.85);
  for (int it = 0; it < n_iters; ++it) {
    std::vector<vec3> Ax;
    apply_A(x, Ax);
    for (size_t i = 0; i < n; ++i) {
      real d = real(1);
      if (i < diag.size() && std::isfinite(diag[i]) && diag[i] > real(0))
        d = diag[i];
      const vec3 a = i < Ax.size() ? Ax[i] : vec3::Zero();
      x[i] += omega * (b[i] - a) / d;
    }
  }
}

class sobolev_matvec;

} // namespace duchamp
} // namespace gaudi

namespace Eigen {
namespace internal {
template <>
struct traits<gaudi::duchamp::sobolev_matvec>
    : public Eigen::internal::traits<Eigen::SparseMatrix<double>> {};
} // namespace internal
} // namespace Eigen

namespace gaudi {
namespace duchamp {

/// Matrix-free A for Eigen::ConjugateGradient. Unknowns are packed x,y,z.
class sobolev_matvec : public Eigen::EigenBase<sobolev_matvec> {
public:
  typedef double Scalar;
  typedef double RealScalar;
  typedef int StorageIndex;
  enum {
    ColsAtCompileTime = Eigen::Dynamic,
    MaxColsAtCompileTime = Eigen::Dynamic,
    IsRowMajor = false
  };

  Eigen::Index n = 0;
  std::function<void(const Eigen::VectorXd &, Eigen::VectorXd &)> apply;

  Eigen::Index rows() const { return n; }
  Eigen::Index cols() const { return n; }

  template <typename Rhs>
  Eigen::Product<sobolev_matvec, Rhs, Eigen::AliasFreeProduct>
  operator*(const Eigen::MatrixBase<Rhs> &x) const {
    return Eigen::Product<sobolev_matvec, Rhs, Eigen::AliasFreeProduct>(
        *this, x.derived());
  }
};

} // namespace duchamp
} // namespace gaudi

namespace Eigen {
namespace internal {
template <typename Rhs>
struct generic_product_impl<gaudi::duchamp::sobolev_matvec, Rhs, SparseShape,
                            DenseShape, GemvProduct>
    : generic_product_impl_base<
          gaudi::duchamp::sobolev_matvec, Rhs,
          generic_product_impl<gaudi::duchamp::sobolev_matvec, Rhs>> {
  typedef typename Product<gaudi::duchamp::sobolev_matvec, Rhs>::Scalar Scalar;
  template <typename Dest>
  static void scaleAndAddTo(Dest &dst, const gaudi::duchamp::sobolev_matvec &lhs,
                            const Rhs &rhs, const Scalar &alpha) {
    Eigen::VectorXd Ax;
    lhs.apply(rhs, Ax);
    dst += alpha * Ax;
  }
};
} // namespace internal
} // namespace Eigen

namespace gaudi {
namespace duchamp {

/// CG smooth of Ag≈b for n_iters, starting from x (Appendix B.3.1).
/// Identity preconditioner, matching EigenCG in the reference hierarchy.
/// tol > 0 stops on the relative residual; otherwise Eigen keeps its default.
template <typename ApplyA>
inline void sobolev_cg_smooth(const std::vector<vec3> &b, int n_iters,
                              std::vector<vec3> &x, ApplyA &&apply_A,
                              real tol = 0) {
  const size_t n = b.size();
  if (n == 0 || n_iters <= 0)
    return;
  if (x.size() != n)
    x.assign(n, vec3::Zero());

  sobolev_matvec A;
  A.n = static_cast<Eigen::Index>(3 * n);
  A.apply = [&](const Eigen::VectorXd &v, Eigen::VectorXd &out) {
    std::vector<vec3> u(n), Au;
    for (size_t i = 0; i < n; ++i)
      u[i] = vec3(v(3 * static_cast<Eigen::Index>(i)),
                  v(3 * static_cast<Eigen::Index>(i) + 1),
                  v(3 * static_cast<Eigen::Index>(i) + 2));
    apply_A(u, Au);
    out.resize(A.n);
    for (size_t i = 0; i < n; ++i) {
      const vec3 a = (i < Au.size()) ? Au[i] : vec3::Zero();
      out(3 * static_cast<Eigen::Index>(i)) = a[0];
      out(3 * static_cast<Eigen::Index>(i) + 1) = a[1];
      out(3 * static_cast<Eigen::Index>(i) + 2) = a[2];
    }
  };

  Eigen::VectorXd rhs(A.n), guess(A.n);
  for (size_t i = 0; i < n; ++i) {
    rhs(3 * static_cast<Eigen::Index>(i)) = b[i][0];
    rhs(3 * static_cast<Eigen::Index>(i) + 1) = b[i][1];
    rhs(3 * static_cast<Eigen::Index>(i) + 2) = b[i][2];
    guess(3 * static_cast<Eigen::Index>(i)) = x[i][0];
    guess(3 * static_cast<Eigen::Index>(i) + 1) = x[i][1];
    guess(3 * static_cast<Eigen::Index>(i) + 2) = x[i][2];
  }

  Eigen::ConjugateGradient<sobolev_matvec, Eigen::Lower | Eigen::Upper,
                           Eigen::IdentityPreconditioner>
      cg;
  cg.compute(A);
  // One conjugate sequence for the whole budget, as in VCycleSolve.
  // Restarting every step drops the previous search direction.
  cg.setMaxIterations(n_iters);
  if (tol > real(0))
    cg.setTolerance(static_cast<double>(tol));
  const Eigen::VectorXd sol = cg.solveWithGuess(rhs, guess);
  if (!sol.allFinite())
    return;
  for (size_t i = 0; i < n; ++i)
    x[i] = vec3(sol(3 * static_cast<Eigen::Index>(i)),
                sol(3 * static_cast<Eigen::Index>(i) + 1),
                sol(3 * static_cast<Eigen::Index>(i) + 2));
}

/// Largest Ritz value of `apply_A` from a degree-6 Arnoldi run.
/// The seed is the current residual, which already lives in the projected
/// range. A zero residual falls back to an alternating vector so the
/// estimate is not the constant kernel.
template <typename ApplyA>
inline real sobolev_arnoldi_lambda(const std::vector<vec3> &b,
                                   const std::vector<vec3> &x,
                                   ApplyA &&apply_A) {
  const size_t n = b.size();
  if (n == 0)
    return 0;
  std::vector<vec3> v(n), Av;
  apply_A(x, Av);
  real r2 = 0;
  for (size_t i = 0; i < n; ++i) {
    const vec3 a = (i < Av.size()) ? Av[i] : vec3::Zero();
    v[i] = (i < b.size() ? b[i] : vec3::Zero()) - a;
    r2 += v[i].squaredNorm();
  }
  if (!(r2 > real(0))) {
    for (size_t i = 0; i < n; ++i)
      v[i] = (i & 1) ? vec3(1, 0, 0) : vec3(-1, 0, 0);
  }
  auto matvec = [&](const std::vector<vec3> &u) {
    std::vector<vec3> out;
    apply_A(u, out);
    if (out.size() != u.size())
      out.assign(u.size(), vec3::Zero());
    return out;
  };
  const int k = static_cast<int>(std::min<size_t>(6, n));
  const std::vector<ritz_pair> pairs = arnoldi(matvec, v, k);
  if (pairs.empty())
    return 0;
  const real lambda = pairs.front().eigenvalue;
  return (lambda > real(0) && std::isfinite(lambda)) ? lambda : real(0);
}

/// Chebyshev semi-iteration on [λ_min, λ_max]. Identity preconditioner.
/// Templates for the Solution of Linear Systems, Barrett et al.
/// `lambda_max` is the largest Ritz value; the interval is padded
/// up by 1.1 so a slightly low estimate does not make the polynomial grow.
template <typename ApplyA>
inline void sobolev_chebyshev_smooth(const std::vector<vec3> &b, int n_iters,
                                     real lambda_max, real lower,
                                     std::vector<vec3> &x, ApplyA &&apply_A) {
  const size_t n = b.size();
  if (n == 0 || n_iters <= 0 || !(lambda_max > real(0)))
    return;
  if (x.size() != n)
    x.assign(n, vec3::Zero());
  real frac = lower;
  if (!(frac > real(0)) || frac >= real(1))
    frac = real(0.3);
  const real lmax = real(1.1) * lambda_max;
  const real lmin = frac * lambda_max;
  const real d = real(0.5) * (lmax + lmin);
  const real c = real(0.5) * (lmax - lmin);
  if (!(d > real(0)) || !(c > real(0)))
    return;

  std::vector<vec3> y = x;
  std::vector<vec3> Ax, r(n), p(n), Ap;
  apply_A(y, Ax);
  for (size_t i = 0; i < n; ++i) {
    const vec3 a = (i < Ax.size()) ? Ax[i] : vec3::Zero();
    r[i] = b[i] - a;
  }
  real alpha = 0;
  for (int it = 0; it < n_iters; ++it) {
    if (it == 0) {
      p = r;
      alpha = real(1) / d;
    } else {
      const real beta = (c * alpha / real(2)) * (c * alpha / real(2));
      const real denom = d - beta / alpha;
      if (!(std::abs(denom) > real(0)) || !std::isfinite(denom))
        return;
      alpha = real(1) / denom;
      for (size_t i = 0; i < n; ++i)
        p[i] = r[i] + beta * p[i];
    }
    apply_A(p, Ap);
    bool finite = std::isfinite(alpha);
    for (size_t i = 0; i < n && finite; ++i) {
      const vec3 a = (i < Ap.size()) ? Ap[i] : vec3::Zero();
      finite = std::isfinite(a[0]) && std::isfinite(a[1]) && std::isfinite(a[2]);
    }
    if (!finite)
      return;
    for (size_t i = 0; i < n; ++i) {
      const vec3 a = (i < Ap.size()) ? Ap[i] : vec3::Zero();
      y[i] += alpha * p[i];
      r[i] -= alpha * a;
    }
  }
  for (const vec3 &a : y) {
    if (!std::isfinite(a[0]) || !std::isfinite(a[1]) || !std::isfinite(a[2]))
      return;
  }
  x.swap(y);
}

/// Smoother on every level above the coarse saddle. `solver.type` selects
/// Jacobi, Eigen CG, or Chebyshev for `solver.iter` steps. The coarsest rod
/// does not call this; it factors the KKT system. `lambda_max`, when set,
/// caches the power estimate for this level's frozen operator.
template <typename ApplyA>
inline void sobolev_relax(const sobolev_tangent_point_force_config &cfg,
                          const std::vector<vec3> &b,
                          const std::vector<real> &diag, std::vector<vec3> &x,
                          ApplyA &&apply_A, real *lambda_max = nullptr) {
  const int n = std::max(1, cfg.solver.iter);
  if (cfg.solver.type == sobolev_solver::cg) {
    sobolev_cg_smooth(b, n, x, apply_A);
    return;
  }
  if (cfg.solver.type == sobolev_solver::chebyshev) {
    real lam = (lambda_max && *lambda_max > real(0)) ? *lambda_max : real(0);
    if (!(lam > real(0))) {
      lam = sobolev_arnoldi_lambda(b, x, apply_A);
      if (lambda_max)
        *lambda_max = lam;
    }
    sobolev_chebyshev_smooth(b, n, lam, cfg.solver.chebyshev_lower, x, apply_A);
    return;
  }
  sobolev_jacobi_relax(b, diag, n, x, apply_A);
}

/// |projected| is roundoff next to |raw|. The constrained load is zero.
inline bool rhs_is_roundoff(const std::vector<vec3> &raw,
                            const std::vector<vec3> &projected) {
  real p2 = 0;
  real r2 = 0;
  const size_t n = std::max(raw.size(), projected.size());
  for (size_t i = 0; i < n; ++i) {
    if (i < projected.size())
      p2 += projected[i].squaredNorm();
    if (i < raw.size())
      r2 += raw[i].squaredNorm();
  }
  if (!(p2 > real(0)))
    return true;
  if (!(r2 > real(0)))
    return false;
  return p2 <= r2 * real(1e-24);
}

/// Solve Ag ≈ f by CG; returns g (same space as f).
inline std::vector<vec3>
sobolev_jacobi_solve(asawa::rod::rod &rod, const std::vector<vec3> &f,
                     const sobolev_tangent_point_force_config &cfg,
                     real l0_abs, arp::edge_cov_tree<2>::ptr cov = {},
                     const sobolev_mg_hook &hook = {}) {
  const size_t n = f.size();
  std::vector<vec3> g(n, vec3::Zero());
  if (n == 0 || cfg.solver.iter <= 0)
    return f;

  if (!cov)
    cov = calder::build_rod_cov_tree(rod);
  const auto E = sobolev_detail::build_edges(rod);
  sobolev_detail::cov_weight_rows wB, wB0;
  if (cfg.alpha_B != real(0))
    wB = sobolev_detail::build_cov_weights(rod, E, cfg.p1, l0_abs, cfg.bh_eps,
                                           cfg.cone_eps, /*use_tp_k=*/false,
                                           *cov);
  if (cfg.alpha_0 != real(0))
    wB0 = sobolev_detail::build_cov_weights(rod, E, cfg.p1, l0_abs, cfg.bh_eps,
                                            cfg.cone_eps, /*use_tp_k=*/true,
                                            *cov);
  const sobolev_detail::cov_weight_rows *wBp =
      cfg.alpha_B != real(0) ? &wB : nullptr;
  const sobolev_detail::cov_weight_rows *wB0p =
      cfg.alpha_0 != real(0) ? &wB0 : nullptr;
  sobolev_projector P;
  sobolev_gram G;
  if (hook.build)
    hook.build(rod, P, G, hook.ctx);
  auto apply = [&](const std::vector<vec3> &u, std::vector<vec3> &Au) {
    std::vector<vec3> up = u;
    if (P.on)
      P.apply(up);
    sobolev_detail::matvec_A(rod, cfg.alpha_B, cfg.alpha_0, cfg.p1, l0_abs,
                             cfg.bh_eps, up, Au, cov.get(), cfg.cone_eps, wBp,
                             wB0p);
    G.apply(up, Au);
    if (P.on)
      P.apply(Au);
  };
  std::vector<vec3> rhs = f;
  if (P.on)
    P.apply(rhs);
  // A load in the row space of C arrives as roundoff. CG treats that as a
  // right-hand side and walks the kernel. The constrained solution is 0.
  if (rhs_is_roundoff(f, rhs))
    return g;
  std::vector<real> diag;
  if (cfg.solver.type == sobolev_solver::jacobi) {
    diag = sobolev_detail::jacobi_diag(rod, cfg.alpha_B, cfg.alpha_0, cfg.p1,
                                       l0_abs, cfg.bh_eps, cov.get(),
                                       cfg.cone_eps, wBp, wB0p);
    G.add_diag(diag);
  }
  real lambda_max = 0;
  sobolev_relax(cfg, rhs, diag, g, apply, &lambda_max);
  if (P.on)
    P.apply(g);
  return g;
}

namespace sobolev_mg_detail {

struct level_ops {
  asawa::rod::rod::ptr rod;
  // Covariance tree of this level's rod. A coarse level builds its own from
  // that rod. Storing the fine tree here would cluster the wrong edges and
  // smooth out the scale the residual forest is meant to keep.
  arp::edge_cov_tree<2>::ptr cov;
  // Pair weights of this level's rod. Built once for the frozen pose.
  sobolev_detail::cov_weight_rows wB;
  sobolev_detail::cov_weight_rows wB0;
  // Edge table and length pyramid. φ is the only field re-summed per multiply.
  sobolev_detail::level_frozen frozen;
  // Null-space projector for this level's rod. Identity when empty.
  sobolev_projector P;
  // Projective Gram Cᵀ W C for this level's rod. Empty when unused.
  sobolev_gram G;
  // Jacobi diagonal of this level's A. Built once for the frozen pose.
  std::vector<real> diag;
  // Top Rayleigh quotient of this level's projected matvec. Filled on the
  // first Chebyshev smooth and reused for the rest of the solve.
  real lambda_max = 0;
  // n < dense_cg_threshold: all-pairs matvec, no covariance tree.
  bool dense = false;
};

inline void apply_A_level(std::vector<level_ops> &ops,
                          const sobolev_tangent_point_force_config &cfg,
                          real l0_abs, int ell, const std::vector<vec3> &u,
                          std::vector<vec3> &Au, bool project = true) {
  auto &L = ops[static_cast<size_t>(ell)];
  std::vector<vec3> up = u;
  if (project && L.P.on)
    L.P.apply(up);
  if (L.dense) {
    sobolev_detail::matvec_A_dense(*L.rod, cfg.alpha_B, cfg.alpha_0, cfg.p1,
                                   l0_abs, up, Au);
  } else {
    const sobolev_detail::cov_weight_rows *wB =
        cfg.alpha_B != real(0) ? &L.wB : nullptr;
    const sobolev_detail::cov_weight_rows *wB0 =
        cfg.alpha_0 != real(0) ? &L.wB0 : nullptr;
    sobolev_detail::matvec_A(*L.rod, cfg.alpha_B, cfg.alpha_0, cfg.p1, l0_abs,
                             cfg.bh_eps, up, Au, L.cov.get(), cfg.cone_eps, wB,
                             wB0, true, L.frozen.on ? &L.frozen : nullptr);
  }
  L.G.apply(up, Au);
  if (project && L.P.on)
    L.P.apply(Au);
}

inline std::vector<vec3> level_rhs(std::vector<level_ops> &ops, int ell,
                                 const std::vector<vec3> &b) {
  std::vector<vec3> rhs = b;
  if (ops[static_cast<size_t>(ell)].P.on)
    ops[static_cast<size_t>(ell)].P.apply(rhs);
  return rhs;
}

inline real field_peak(const std::vector<vec3> &v) {
  real p = 0.0;
  for (const vec3 &a : v) {
    const real n = a.norm();
    if (std::isfinite(n))
      p = std::max(p, n);
  }
  return p;
}

/// Peak of `rhs - Ax`. A zero `Ax` is the residual of a zero guess, which is `|rhs|`.
inline real residual_peak(const std::vector<vec3> &rhs,
                          const std::vector<vec3> &Ax) {
  real p = 0.0;
  const size_t n = std::max(rhs.size(), Ax.size());
  for (size_t i = 0; i < n; ++i) {
    const vec3 a = i < rhs.size() ? rhs[i] : vec3::Zero();
    const vec3 b = i < Ax.size() ? Ax[i] : vec3::Zero();
    const real m = (a - b).norm();
    if (std::isfinite(m))
      p = std::max(p, m);
  }
  return p;
}

struct level_note {
  int pass = -1; ///< -1 is the FMG guess, otherwise the V-cycle index
  int ell = 0;
  size_t n = 0;
  real b_max = 0; ///< level load, after the constraint projection
  real r_in = 0;  ///< residual peak before the pre-smooth
  real x_pre = 0; ///< field peak after the pre-smooth, before adding J e
  real r_max = 0; ///< residual peak after the pre-smooth
  real ec_max = 0;
  real ef_max = 0;
  real x_add = 0;
  real x_post = 0;
};

inline std::vector<level_note> &level_notes() {
  static thread_local std::vector<level_note> notes;
  return notes;
}

inline int &mg_pass() {
  static thread_local int pass = -1;
  return pass;
}

inline void note_level(int ell, size_t n, real b_max, real r_in, real x_pre,
                       real r_max, real ec_max, real ef_max, real x_add,
                       real x_post) {
  level_notes().push_back(level_note{mg_pass(), ell, n, b_max, r_in, x_pre,
                                     r_max, ec_max, ef_max, x_add, x_post});
}

inline void print_level_notes(real g_max) {
  bool hot = g_max > real(1);
  for (const level_note &s : level_notes()) {
    if (s.x_pre > real(1) || s.ef_max > real(1) || s.x_add > real(1) ||
        s.x_post > real(1) || s.ec_max > real(1))
      hot = true;
  }
  if (!hot)
    return;
  for (const level_note &s : level_notes()) {
    demo_log() << "[sobolev level] pass=" << s.pass << " ell=" << s.ell
              << " n=" << s.n << " |b|_max=" << s.b_max
              << " |r|_in=" << s.r_in << " |x|_pre=" << s.x_pre
              << " |r|_max=" << s.r_max << " |e_c|_max=" << s.ec_max
              << " |e|_max=" << s.ef_max << " |x|_add=" << s.x_add
              << " |x|_post=" << s.x_post << std::endl;
  }
}

inline void coarsest_direct(std::vector<level_ops> &ops,
                            const sobolev_tangent_point_force_config &cfg,
                            real l0_abs, int ell, const std::vector<vec3> &rhs,
                            std::vector<vec3> &x) {
  const int nv = static_cast<int>(rhs.size());
  x.assign(static_cast<size_t>(std::max(nv, 0)), vec3::Zero());
  if (nv <= 0)
    return;
  auto &L = ops[static_cast<size_t>(ell)];
  const int n3 = 3 * nv;
  const bool constrained =
      L.P.on && L.P.B.rows() > 0 && L.P.B.cols() == n3;
  const int m = constrained ? static_cast<int>(L.P.B.rows()) : 0;
  const int N = n3 + m;
  Eigen::MatrixXd S = Eigen::MatrixXd::Zero(N, N);
  std::vector<vec3> u(static_cast<size_t>(nv), vec3::Zero()), Au;
  for (int j = 0; j < n3; ++j) {
    std::fill(u.begin(), u.end(), vec3::Zero());
    u[static_cast<size_t>(j / 3)][j % 3] = real(1);
    apply_A_level(ops, cfg, l0_abs, ell, u, Au, /*project=*/false);
    for (int i = 0; i < nv; ++i) {
      const vec3 a = (static_cast<size_t>(i) < Au.size())
                         ? Au[static_cast<size_t>(i)]
                         : vec3::Zero();
      S(3 * i, j) = a[0];
      S(3 * i + 1, j) = a[1];
      S(3 * i + 2, j) = a[2];
    }
  }
  if (!L.G.empty()) {
    std::vector<vec3> basis(static_cast<size_t>(nv), vec3::Zero());
    std::vector<vec3> Gb(static_cast<size_t>(nv), vec3::Zero());
    for (int j = 0; j < n3; ++j) {
      std::fill(basis.begin(), basis.end(), vec3::Zero());
      std::fill(Gb.begin(), Gb.end(), vec3::Zero());
      basis[static_cast<size_t>(j / 3)][j % 3] = real(1);
      L.G.apply(basis, Gb);
      for (int i = 0; i < nv; ++i) {
        const vec3 a = Gb[static_cast<size_t>(i)];
        S(3 * i, j) += a[0];
        S(3 * i + 1, j) += a[1];
        S(3 * i + 2, j) += a[2];
      }
    }
  }
  if (m > 0) {
    const Eigen::MatrixXd C = Eigen::MatrixXd(L.P.B);
    S.block(n3, 0, m, n3) = C;
    S.block(0, n3, n3, m) = C.transpose();
  }
  Eigen::VectorXd b = Eigen::VectorXd::Zero(N);
  for (int i = 0; i < nv; ++i) {
    b(3 * i) = rhs[static_cast<size_t>(i)][0];
    b(3 * i + 1) = rhs[static_cast<size_t>(i)][1];
    b(3 * i + 2) = rhs[static_cast<size_t>(i)][2];
  }
  const Eigen::VectorXd sol = S.partialPivLu().solve(b);
  for (int i = 0; i < nv; ++i)
    x[static_cast<size_t>(i)] =
        vec3(sol(3 * i), sol(3 * i + 1), sol(3 * i + 2));
}

inline void vcycle(asawa::rod::residual_forest &F, std::vector<level_ops> &ops,
                   const sobolev_tangent_point_force_config &cfg, real l0_abs,
                   int ell, int n_levels, const std::vector<vec3> &b,
                   std::vector<vec3> &x) {
  const std::vector<vec3> rhs = level_rhs(ops, ell, b);
  if (rhs_is_roundoff(b, rhs)) {
    x.assign(rhs.size(), vec3::Zero());
    note_level(ell, x.size(), 0, 0, 0, 0, 0, 0, 0, 0);
    return;
  }
  const real b_max = field_peak(rhs);
  // Coarsest rod: one dense KKT factorization. Finer levels project and Jacobi.
  if (ell + 1 >= n_levels) {
    coarsest_direct(ops, cfg, l0_abs, ell, rhs, x);
    if (ops[static_cast<size_t>(ell)].P.on)
      ops[static_cast<size_t>(ell)].P.apply(x);
    auto apply = [&](const std::vector<vec3> &u, std::vector<vec3> &Au) {
      apply_A_level(ops, cfg, l0_abs, ell, u, Au);
    };
    std::vector<vec3> Ax;
    apply(x, Ax);
    const real x_peak = field_peak(x);
    note_level(ell, x.size(), b_max, b_max, 0, residual_peak(rhs, Ax), 0,
               x_peak, x_peak, x_peak);
    return;
  }
  auto apply = [&](const std::vector<vec3> &u, std::vector<vec3> &Au) {
    apply_A_level(ops, cfg, l0_abs, ell, u, Au);
  };
  real r_in = b_max;
  if (field_peak(x) > real(0)) {
    std::vector<vec3> Ax_in;
    apply(x, Ax_in);
    r_in = residual_peak(rhs, Ax_in);
  }
  sobolev_relax(cfg, rhs, ops[static_cast<size_t>(ell)].diag, x, apply,
                &ops[static_cast<size_t>(ell)].lambda_max);
  if (ops[static_cast<size_t>(ell)].P.on)
    ops[static_cast<size_t>(ell)].P.apply(x);
  const real x_pre = field_peak(x);

  std::vector<vec3> Ax;
  apply(x, Ax);
  std::vector<vec3> r(rhs.size());
  for (size_t i = 0; i < rhs.size(); ++i)
    r[i] = rhs[i] - (i < Ax.size() ? Ax[i] : vec3::Zero());
  const real r_max = field_peak(r);

  // Reference restrictWithTranspose / prolong: project on the way in and out.
  auto &Pf = ops[static_cast<size_t>(ell)].P;
  auto &Pc = ops[static_cast<size_t>(ell + 1)].P;
  if (Pf.on)
    Pf.apply(r);
  std::vector<vec3> b_c = F.restrict_JT(ell, r);
  if (Pc.on)
    Pc.apply(b_c);
  std::vector<vec3> e_c(b_c.size(), vec3::Zero());
  vcycle(F, ops, cfg, l0_abs, ell + 1, n_levels, b_c, e_c);
  if (Pc.on)
    Pc.apply(e_c);
  std::vector<vec3> e_f = F.prolong_J(ell, e_c);
  if (Pf.on)
    Pf.apply(e_f);
  for (size_t i = 0; i < x.size() && i < e_f.size(); ++i)
    x[i] += e_f[i];
  const real x_add = field_peak(x);

  sobolev_relax(cfg, rhs, ops[static_cast<size_t>(ell)].diag, x, apply,
                &ops[static_cast<size_t>(ell)].lambda_max);
  if (ops[static_cast<size_t>(ell)].P.on)
    ops[static_cast<size_t>(ell)].P.apply(x);
  note_level(ell, x.size(), b_max, r_in, x_pre, r_max, field_peak(e_c),
             field_peak(e_f), x_add, field_peak(x));
}

/// Appendix B.3.1 initial guess: restrict b to the coarsest rod, solve,
/// prolong, and smooth on the way up.
inline void fmg_guess(asawa::rod::residual_forest &F,
                      std::vector<level_ops> &ops,
                      const sobolev_tangent_point_force_config &cfg,
                      real l0_abs, int ell, int n_levels,
                      const std::vector<vec3> &b, std::vector<vec3> &x) {
  const std::vector<vec3> rhs = level_rhs(ops, ell, b);
  if (ell + 1 >= n_levels) {
    x.assign(rhs.size(), vec3::Zero());
    vcycle(F, ops, cfg, l0_abs, ell, n_levels, rhs, x);
    return;
  }
  auto &Pf = ops[static_cast<size_t>(ell)].P;
  auto &Pc = ops[static_cast<size_t>(ell + 1)].P;
  std::vector<vec3> b_in = rhs;
  if (Pf.on)
    Pf.apply(b_in);
  std::vector<vec3> b_c = F.restrict_JT(ell, b_in);
  if (Pc.on)
    Pc.apply(b_c);
  std::vector<vec3> e_c(b_c.size(), vec3::Zero());
  fmg_guess(F, ops, cfg, l0_abs, ell + 1, n_levels, b_c, e_c);
  if (Pc.on)
    Pc.apply(e_c);
  x = F.prolong_J(ell, e_c);
  if (Pf.on)
    Pf.apply(x);
  if (x.size() < rhs.size())
    x.resize(rhs.size(), vec3::Zero());
  else if (x.size() > rhs.size())
    x.resize(rhs.size());
  const real b_max = field_peak(rhs);
  const real x_pre = field_peak(x);
  auto apply = [&](const std::vector<vec3> &u, std::vector<vec3> &Au) {
    apply_A_level(ops, cfg, l0_abs, ell, u, Au);
  };
  real r_in = b_max;
  if (x_pre > real(0)) {
    std::vector<vec3> Ax_in;
    apply(x, Ax_in);
    r_in = residual_peak(rhs, Ax_in);
  }
  sobolev_relax(cfg, rhs, ops[static_cast<size_t>(ell)].diag, x, apply,
                &ops[static_cast<size_t>(ell)].lambda_max);
  if (ops[static_cast<size_t>(ell)].P.on)
    ops[static_cast<size_t>(ell)].P.apply(x);
  std::vector<vec3> Ax;
  apply(x, Ax);
  note_level(ell, x.size(), b_max, r_in, x_pre, residual_peak(rhs, Ax), 0,
             x_pre, x_pre, field_peak(x));
}

} // namespace sobolev_mg_detail

/// Covariance tree for one residual-forest level, built once for this pose.
/// Level 0 may reuse `fine` because that pointer is already this rod's tree.
/// Every coarser level builds from its own rod.
inline arp::edge_cov_tree<2>::ptr
sobolev_level_cov(asawa::rod::rod &level_rod, int ell,
                  const arp::edge_cov_tree<2>::ptr &fine) {
  if (ell == 0 && fine)
    return fine;
  return calder::build_rod_cov_tree(level_rod);
}

/// Solve Ag ≈ f by geometric V-cycle (Crane App. B.3) on residual_forest.
inline std::vector<vec3>
sobolev_mg_solve(asawa::rod::rod &rod, const std::vector<vec3> &f,
                 const sobolev_tangent_point_force_config &cfg, real l0_abs,
                 arp::edge_cov_tree<2>::ptr fine = {},
                 const sobolev_mg_hook &hook = {}) {
  const size_t n = f.size();
  if (n == 0)
    return f;
  if (cfg.num_v_cycles <= 0 || cfg.max_levels == 1)
    return sobolev_jacobi_solve(rod, f, cfg, l0_abs, std::move(fine), hook);

  // Non-owning wrap: forest is local to this solve.
  asawa::rod::rod::ptr Rptr(&rod, [](asawa::rod::rod *) {});
  auto F = asawa::rod::residual_forest::create(Rptr);
  if (!F || !F->valid() || F->n_levels() < 2)
    return sobolev_jacobi_solve(rod, f, cfg, l0_abs, std::move(fine), hook);

  int n_levels = F->n_levels();
  if (cfg.max_levels > 0)
    n_levels = std::min(n_levels, cfg.max_levels);
  if (n_levels < 2)
    return sobolev_jacobi_solve(rod, f, cfg, l0_abs, std::move(fine), hook);
  // Extend past max_levels until the coarse rod is under dense_cg_threshold,
  // then stop on that rod. It is solved by Eigen CG with the dense matvec.
  if (cfg.dense_cg_threshold > 0) {
    const auto below = static_cast<size_t>(cfg.dense_cg_threshold);
    while (n_levels < F->n_levels() &&
           F->level_rod(n_levels - 1)->corner_count() >= below)
      ++n_levels;
    for (int ell = 1; ell < n_levels; ++ell) {
      if (F->level_rod(ell)->corner_count() < below) {
        n_levels = ell + 1;
        break;
      }
    }
  }

  // Sync coarse geometry from fine survivors; freeze for this solve.
  for (int ell = 0; ell + 1 < n_levels; ++ell)
    F->restrict_up(ell);

  std::vector<sobolev_mg_detail::level_ops> ops(
      static_cast<size_t>(n_levels));
  for (int ell = 0; ell < n_levels; ++ell) {
    auto &op = ops[static_cast<size_t>(ell)];
    op.rod = F->level_rod(ell);
    op.dense = ell > 0 && cfg.dense_cg_threshold > 0 &&
               op.rod->corner_count() <
                   static_cast<size_t>(cfg.dense_cg_threshold);
    if (hook.build)
      hook.build(*op.rod, op.P, op.G, hook.ctx);
    if (op.dense)
      continue;
    op.cov = sobolev_level_cov(*op.rod, ell, fine);
    op.frozen = sobolev_detail::freeze_level(*op.rod, *op.cov);
    const auto &E = op.frozen.E;
    if (cfg.alpha_B != real(0))
      op.wB = sobolev_detail::build_cov_weights(
          *op.rod, E, cfg.p1, l0_abs, cfg.bh_eps, cfg.cone_eps,
          /*use_tp_k=*/false, *op.cov);
    if (cfg.alpha_0 != real(0))
      op.wB0 = sobolev_detail::build_cov_weights(
          *op.rod, E, cfg.p1, l0_abs, cfg.bh_eps, cfg.cone_eps,
          /*use_tp_k=*/true, *op.cov);
    const sobolev_detail::cov_weight_rows *wB =
        cfg.alpha_B != real(0) ? &op.wB : nullptr;
    const sobolev_detail::cov_weight_rows *wB0 =
        cfg.alpha_0 != real(0) ? &op.wB0 : nullptr;
    op.diag = sobolev_detail::jacobi_diag(
        *op.rod, cfg.alpha_B, cfg.alpha_0, cfg.p1, l0_abs, cfg.bh_eps,
        op.cov.get(), cfg.cone_eps, wB, wB0);
    op.G.add_diag(op.diag);
  }

  sobolev_mg_detail::level_notes().clear();
  sobolev_mg_detail::mg_pass() = -1;
  std::vector<vec3> g(n, vec3::Zero());
  sobolev_mg_detail::fmg_guess(*F, ops, cfg, l0_abs, 0, n_levels, f, g);
  const int cycles = std::max(1, cfg.num_v_cycles);
  for (int c = 0; c < cycles; ++c) {
    sobolev_mg_detail::mg_pass() = c;
    sobolev_mg_detail::vcycle(*F, ops, cfg, l0_abs, 0, n_levels, f, g);
  }
  sobolev_mg_detail::print_level_notes(sobolev_mg_detail::field_peak(g));
  sobolev_mg_detail::level_notes().clear();

  return g;
}

/// dÊ, before A^{-1}. The tree is the one the leaf evaluation just built.
struct sobolev_tp_rhs {
  std::vector<vec3> f;
  real l0 = 0.0;
  arp::edge_cov_tree<2>::ptr cov;
};

inline sobolev_tp_rhs
sobolev_tp_energy_grad(asawa::rod::rod &rod, const asawa::rod::dynamic &dynamic,
                       const sobolev_tangent_point_force_config &cfg) {
  sobolev_tp_rhs out;
  out.l0 = cfg.l0 * dynamic._Cc;
  out.cov = calder::build_rod_cov_tree(rod);
  auto ev = sobolev_detail::eval_paper_tp(rod, cfg.tpe_alpha, cfg.tpe_beta,
                                          cfg.bh_eps);
  out.f = std::move(ev.grad);
  if (cfg.w != real(1.0)) {
    for (auto &a : out.f)
      a *= cfg.w;
  }
  auto add_binormal = [&](real weight, const std::vector<vec3> &force) {
    if (weight == real(0))
      return;
    const auto &x = rod.x();
    for (size_t i = 0; i < out.f.size() && i < force.size(); ++i) {
      const auto ci = asawa::rod::corner_id(static_cast<int>(i));
      const auto cp = rod.prev(ci);
      const auto cn = rod.next(ci);
      if (cp < asawa::rod::corner_id(0) || cn < asawa::rod::corner_id(0))
        continue;
      const size_t ip = static_cast<size_t>(static_cast<int>(cp));
      const size_t in = static_cast<size_t>(static_cast<int>(cn));
      if (ip >= x.size() || in >= x.size())
        continue;
      const vec3 tang = x[in] - x[ip];
      const vec3 curv = x[ip] - real(2) * x[i] + x[in];
      vec3 B = tang.cross(curv);
      const real b2 = B.squaredNorm();
      if (!(b2 > real(1e-24)))
        continue;
      B /= std::sqrt(b2);
      out.f[i] += weight * force[i].dot(B) * B;
    }
  };
  if (cfg.coulomb_w != real(0))
    add_binormal(cfg.coulomb_w,
                 calder::coulomb_force(rod, rod.x(), out.l0, cfg.coulomb_p));
  if (cfg.vortex_w != real(0))
    add_binormal(cfg.vortex_w,
                 calder::vortex_force(rod, rod.x(), rod.l0(), out.l0,
                                      cfg.vortex_p, cfg.vortex_q));
  return out;
}

/// dÊ then Ag = dÊ. The fine tree is reused only on level 0. Each coarser
/// forest rod gets its own covariance tree for the diagonal and every matvec.
inline std::vector<vec3>
sobolev_tp_gradient(asawa::rod::rod &rod, const asawa::rod::dynamic &dynamic,
                    const sobolev_tangent_point_force_config &cfg) {
  sobolev_tp_rhs rhs = sobolev_tp_energy_grad(rod, dynamic, cfg);
  return sobolev_mg_solve(rod, rhs.f, cfg, rhs.l0, std::move(rhs.cov));
}

inline real sobolev_quadratic_B(asawa::rod::rod &rod, real p, real l0_abs,
                                real bh_eps, const std::vector<vec3> &u,
                                const std::vector<vec3> &v) {
  std::vector<vec3> Bv;
  sobolev_detail::matvec_B(rod, v, p, l0_abs, bh_eps, Bv);
  real s = 0;
  for (size_t i = 0; i < u.size(); ++i)
    s += u[i].dot(Bv[i]);
  return s;
}

inline real sobolev_quadratic_B0(asawa::rod::rod &rod, real p, real l0_abs,
                                 real bh_eps, const std::vector<vec3> &u,
                                 const std::vector<vec3> &v) {
  std::vector<vec3> B0v;
  sobolev_detail::matvec_B0(rod, v, p, l0_abs, bh_eps, B0v);
  real s = 0;
  for (size_t i = 0; i < u.size(); ++i)
    s += u[i].dot(B0v[i]);
  return s;
}

} // namespace duchamp
} // namespace gaudi

#endif // GAUDI_DUCHAMP_TP_SOBOLEV_JACOBI_HPP
