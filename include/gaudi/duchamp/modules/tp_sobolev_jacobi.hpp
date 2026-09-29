#ifndef GAUDI_DUCHAMP_TP_SOBOLEV_JACOBI_HPP
#define GAUDI_DUCHAMP_TP_SOBOLEV_JACOBI_HPP

/// Sobolev metric A = M + α_B B + α_0 B_0 and Jacobi Ag=f.
/// B / B_0 follow Crane edge–edge assembly (Eq. 31–32): vertex fields lift to
/// edges via E (average) or D (1/ℓ difference), then Eᵀ / Dᵀ back to vertices.
/// Leaf K is a quadrature over parameters on each edge. The list is the paper
/// trapezoid {0, 1} with weights {1/2, 1/2}. A far node stays one tangent-point
/// lump. B0: K_IJ = (k̂+k̂_swap)·ℓ_I·ℓ_J with k_2^{p+4}.
/// B: the same with k_0^p, and the edge residual is φ_I − ⟨T_I, T_J⟩ φ_J
/// (Eq. 11, §5.2.1–5.2.2). B0 stays the plain average φ_I − φ_J.
/// Writeup: docs/sobolev_B_tangent_filter.md.
/// The Sobolev right-hand side is dÊ for Ê = k_2^{p+2} on that same quadrature.

#include "gaudi/asawa/rod/residual_forest.hpp"
#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/calder/rod_integrators.hpp"
#include "gaudi/calder/tree_code.hpp"
#include "gaudi/common.h"
#include "gaudi/duchamp/modules/tp_config.hpp"

#include <algorithm>
#include <cmath>
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

/// Per-coordinate slope (u1 - u0) / ℓ. ⟨T_I, T_J⟩ is applied in the B residual.
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

inline void scatter_DT(const edge_table &E, const std::vector<vec3> &r_edge,
                       std::vector<vec3> &out) {
  out.assign(E.n_vert, vec3::Zero());
  for (size_t i = 0; i < E.n; ++i) {
    const vec3 s = r_edge[i] / E.ell[i];
    out[static_cast<size_t>(E.ep0[i])] -= s;
    out[static_cast<size_t>(E.ep1[i])] += s;
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

/// Leaf uses the quadrature points on both edges. Far uses the lumped point.
inline real sample_pair_weight(const edge_table &E, size_t ii, const vec3 &pi,
                               const vec3 &pj, const vec3 &Tj, real lj, int j0,
                               int j1, const std::vector<vec3> &x,
                               calder::Node_Type node_type, real mass_scale,
                               real eps, real p, bool use_tp_k,
                               const edge_quad &quad) {
  const real li = mass_scale * E.ell[ii];
  if (node_type == calder::LEAF) {
    if (leaf_edges_adjacent_or_self(E, ii, j0, j1))
      return real(0);
    if (j0 < 0 || j1 < 0)
      return real(0);
    const size_t a0 = static_cast<size_t>(E.ep0[ii]);
    const size_t a1 = static_cast<size_t>(E.ep1[ii]);
    const size_t b0 = static_cast<size_t>(j0);
    const size_t b1 = static_cast<size_t>(j1);
    if (a0 >= x.size() || a1 >= x.size() || b0 >= x.size() || b1 >= x.size())
      return real(0);
    return quad_pair_weight(x[a0], x[a1], E.T[ii], x[b0], x[b1], Tj, li, lj,
                            eps, p, use_tp_k, quad);
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
/// `align` is ⟨T_I, T_J⟩ for B and 1 for B0.
struct cov_weight_hit {
  index_t j = 0;
  real w = 0;
  bool far = false;
  real align = 1;
};

inline real tangent_align(const vec3 &Ti, const vec3 &Tj) {
  const real ni = Ti.norm();
  const real nj = Tj.norm();
  if (!(ni > real(1e-12)) || !(nj > real(1e-12)))
    return real(0);
  return Ti.dot(Tj) / (ni * nj);
}

/// Per-query-edge samples for one frozen pose. `sum` is Σ w, the Jacobi row.
struct cov_weight_rows {
  std::vector<std::vector<cov_weight_hit>> row;
  std::vector<real> sum;
};

inline void bind_edge_lengths(const edge_table &E, real mass_scale,
                              std::vector<real> &ell_b, std::vector<real> &ep0_b,
                              std::vector<real> &ep1_b) {
  const size_t n = E.n;
  ell_b.resize(n);
  ep0_b.resize(n);
  ep1_b.resize(n);
  for (size_t i = 0; i < n; ++i) {
    ell_b[i] = mass_scale * E.ell[i];
    ep0_b[i] = real(E.ep0[i]);
    ep1_b[i] = real(E.ep1[i]);
  }
}

/// Walk once and store w. Later matvecs pyramid φ and multiply.
inline cov_weight_rows
build_cov_weights(asawa::rod::rod &rod, const edge_table &E, real p,
                  real l0_abs, real mass_scale, real bh_eps, real cone_eps,
                  bool use_tp_k, const arp::edge_cov_tree<2> &cov,
                  const edge_quad &quad = trapezoid_edge_quad()) {
  const size_t n = E.n;
  cov_weight_rows rows;
  rows.row.resize(n);
  rows.sum.assign(n, real(0));
  if (n == 0)
    return rows;

  std::vector<real> ell_b, ep0_b, ep1_b;
  bind_edge_lengths(E, mass_scale, ell_b, ep0_b, ep1_b);
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
  auto kernel = [&E, &rows, &x, &quad, p, eps, use_tp_k, mass_scale](
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
                                       node_type, mass_scale, eps, p, use_tp_k,
                                       quad);
    if (!(w > real(0)) || !std::isfinite(w))
      return vec3::Zero();
    const real align = use_tp_k ? real(1) : tangent_align(E.T[ii], Tj);
    rows.row[ii].push_back(
        cov_weight_hit{j, w, node_type == calder::BRANCH, align});
    rows.sum[ii] += w;
    return vec3::Zero();
  };
  calder::integrate_over_cov_tree<vec3>(cov, ids, E.mid, bind, kernel, kernel,
                                        bh_eps, 1, cone_eps);
  return rows;
}

/// r_I = Σ w (φ_I − align φ_J). Far φ_J is the pyramided node sum.
/// B stores align = ⟨T_I, T_J⟩. B0 stores align = 1.
inline std::vector<vec3>
apply_cov_weights(const edge_table &E, const std::vector<vec3> &phi,
                  const cov_weight_rows &rows,
                  const arp::edge_cov_tree<2> &cov) {
  const size_t n = E.n;
  std::vector<vec3> r(n, vec3::Zero());
  if (n == 0)
    return r;
  std::vector<vec3> phi_b = phi;
  if (phi_b.size() < n)
    phi_b.resize(n, vec3::Zero());
  // datum holds references, so the index vector has to outlive the pyramid.
  const std::vector<index_t> edge_ids = edge_ordinals(n);
  auto phi_datum = calder::vec3_datum::create(edge_ids, phi_b);
  phi_datum->pyramid(cov);
  const std::vector<vec3> &leaf = phi_datum->sorted_leaf_data();
  const std::vector<vec3> &node = phi_datum->node_data();

  for (size_t i = 0; i < n && i < rows.row.size(); ++i) {
    vec3 acc = vec3::Zero();
    for (const cov_weight_hit &hit : rows.row[i]) {
      const vec3 phij =
          hit.far ? (static_cast<size_t>(hit.j) < node.size() ? node[hit.j]
                                                              : vec3::Zero())
                  : (static_cast<size_t>(hit.j) < leaf.size() ? leaf[hit.j]
                                                              : vec3::Zero());
      acc += hit.w * (phi_b[i] - hit.align * phij);
    }
    r[i] = acc;
  }
  return r;
}

/// Edge residual r_I = Σ_J w_IJ (φ_I − align φ_J). Opening query is the midpoint.
/// align is ⟨T_I, T_J⟩ for B and 1 for B0.
inline std::vector<vec3>
edge_K_diff(asawa::rod::rod &rod, const edge_table &E,
            const std::vector<vec3> &phi, real p, real l0_abs, real mass_scale,
            real bh_eps, bool use_tp_k,
            const arp::edge_cov_tree<2> *cov = nullptr, real cone_eps = -1,
            const cov_weight_rows *weights = nullptr,
            const edge_quad &quad = trapezoid_edge_quad()) {
  const size_t n = E.n;
  std::vector<vec3> r(n, vec3::Zero());
  if (n == 0)
    return r;
  if (weights && cov)
    return apply_cov_weights(E, phi, *weights, *cov);

  std::vector<real> ell_b, ep0_b, ep1_b;
  bind_edge_lengths(E, mass_scale, ell_b, ep0_b, ep1_b);
  const std::vector<vec3> T_b = E.T;
  std::vector<vec3> phi_b = phi;
  if (phi_b.size() < n)
    phi_b.resize(n, vec3::Zero());
  const std::vector<vec3> &x = rod.x();
  const real eps = std::max(l0_abs, real(1e-12));
  const std::vector<index_t> ids = edge_ordinals(n);

  auto bind = [&ell_b, &T_b, &phi_b, &ep0_b, &ep1_b](
                  const std::vector<index_t> &edge_ids, calder::cov_rod_sum &sum) {
    sum.bind(calder::scalar_datum::create(edge_ids, ell_b));
    sum.bind(calder::vec3_datum::create(edge_ids, T_b));
    sum.bind(calder::vec3_datum::create(edge_ids, phi_b));
    sum.bind(calder::scalar_datum::create(edge_ids, ep0_b));
    sum.bind(calder::scalar_datum::create(edge_ids, ep1_b));
  };
  auto kernel = [&E, &phi, &x, &quad, p, eps, use_tp_k, mass_scale](
                    const index_t i, const index_t j, const vec3 &pi,
                    const vec3 &pj, const std::vector<calder::datum::ptr> &data,
                    calder::Node_Type node_type,
                    const arp::edge_cov_tree<2> &) -> vec3 {
        const size_t ii = static_cast<size_t>(i);
        if (ii >= E.n)
          return vec3::Zero();

        const real lj = calder::get_data<real>(node_type, j, 0, data);
        vec3 Tj = calder::get_data<vec3>(node_type, j, 1, data);
        const vec3 phij = calder::get_data<vec3>(node_type, j, 2, data);
        const int j0 =
            static_cast<int>(calder::get_data<real>(node_type, j, 3, data));
        const int j1 =
            static_cast<int>(calder::get_data<real>(node_type, j, 4, data));
        if (!(lj > real(0)))
          return vec3::Zero();
        const real w = sample_pair_weight(E, ii, pi, pj, Tj, lj, j0, j1, x,
                                           node_type, mass_scale, eps, p,
                                           use_tp_k, quad);
        if (!(w > real(0)) || !std::isfinite(w))
          return vec3::Zero();
        const real align = use_tp_k ? real(1) : tangent_align(E.T[ii], Tj);
        return w * (phi[ii] - align * phij);
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
               real mass_scale, real bh_eps, bool use_tp_k,
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
  bind_edge_lengths(E, mass_scale, ell_b, ep0_b, ep1_b);
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
  auto kernel = [&E, &x, &quad, p, eps, use_tp_k, mass_scale](
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
                                           node_type, mass_scale, eps, p,
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

inline void matvec_B0(asawa::rod::rod &rod, const std::vector<vec3> &u,
                      const std::vector<real> & /*m_pov*/, real p, real l0_abs,
                      real mass_scale, real bh_eps, std::vector<vec3> &out,
                      const arp::edge_cov_tree<2> *cov = nullptr,
                      real cone_eps = -1,
                      const cov_weight_rows *weights = nullptr) {
  const edge_table E = build_edges(rod);
  const std::vector<vec3> Eu = edge_avg_E(E, u);
  const std::vector<vec3> r =
      edge_K_diff(rod, E, Eu, p, l0_abs, mass_scale, bh_eps, /*use_tp_k=*/true,
                  cov, cone_eps, weights);
  scatter_ET(E, r, out);
}

inline void matvec_B(asawa::rod::rod &rod, const std::vector<vec3> &u,
                     const std::vector<real> & /*m_pov*/, real p, real l0_abs,
                     real mass_scale, real bh_eps, std::vector<vec3> &out,
                     const arp::edge_cov_tree<2> *cov = nullptr,
                     real cone_eps = -1,
                     const cov_weight_rows *weights = nullptr) {
  const edge_table E = build_edges(rod);
  const std::vector<vec3> Du = edge_diff_D(E, u);
  const std::vector<vec3> r =
      edge_K_diff(rod, E, Du, p, l0_abs, mass_scale, bh_eps, /*use_tp_k=*/false,
                  cov, cone_eps, weights);
  scatter_DT(E, r, out);
}

inline void matvec_A(asawa::rod::rod &rod, const std::vector<real> &m_pov,
                     real alpha_B, real alpha_0, real p, real l0_abs,
                     real mass_scale, real bh_eps, const std::vector<vec3> &u,
                     std::vector<vec3> &out,
                     const arp::edge_cov_tree<2> *cov = nullptr,
                     real cone_eps = -1, const cov_weight_rows *wB = nullptr,
                     const cov_weight_rows *wB0 = nullptr) {
  const size_t n = u.size();
  out.assign(n, vec3::Zero());
  for (size_t i = 0; i < n; ++i)
    out[i] = m_pov[i] * u[i];

  if (alpha_0 != real(0)) {
    std::vector<vec3> B0u;
    matvec_B0(rod, u, m_pov, p, l0_abs, mass_scale, bh_eps, B0u, cov, cone_eps,
              wB0);
    for (size_t i = 0; i < n; ++i)
      out[i] += alpha_0 * B0u[i];
  }
  if (alpha_B != real(0)) {
    std::vector<vec3> Bu;
    matvec_B(rod, u, m_pov, p, l0_abs, mass_scale, bh_eps, Bu, cov, cone_eps,
             wB);
    for (size_t i = 0; i < n; ++i)
      out[i] += alpha_B * Bu[i];
  }
}

inline std::vector<real> jacobi_diag(
    asawa::rod::rod &rod, const std::vector<real> &m_pov, real alpha_B,
    real alpha_0, real p, real l0_abs, real mass_scale, real bh_eps,
    const arp::edge_cov_tree<2> *cov = nullptr, real cone_eps = -1,
    const cov_weight_rows *wB = nullptr, const cov_weight_rows *wB0 = nullptr) {
  const size_t n = m_pov.size();
  std::vector<real> d = m_pov;
  for (real &x : d)
    x = std::max(x, real(1e-12));

  const edge_table E = build_edges(rod);

  if (alpha_0 != real(0)) {
    const auto S0 =
        edge_K_row_sum(rod, E, p, l0_abs, mass_scale, bh_eps, /*use_tp_k=*/true,
                       cov, cone_eps, wB0);
    for (size_t i = 0; i < E.n; ++i) {
      const real lump = real(0.5) * alpha_0 * S0[i];
      d[static_cast<size_t>(E.ep0[i])] += lump;
      d[static_cast<size_t>(E.ep1[i])] += lump;
    }
  }
  if (alpha_B != real(0)) {
    const auto SB = edge_K_row_sum(rod, E, p, l0_abs, mass_scale, bh_eps,
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

struct energy_scatter {
  real e = 0;
  int n = 0;
  int v[4] = {-1, -1, -1, -1};
  vec3 g0 = vec3::Zero();
  vec3 g1 = vec3::Zero();
  vec3 g2 = vec3::Zero();
  vec3 g3 = vec3::Zero();
};

/// Ê and dÊ on the covariance walk. Leaf visits add half the symmetric trapezoid
/// term (each unordered pair is seen twice). A far visit adds the lumped kernel
/// and scatters only onto the query edge. Lengths are geometric; cfg.mass stays
/// in the metric.
inline paper_tp_eval eval_paper_tp(asawa::rod::rod &rod, real l0_abs, real p,
                                   real bh_eps, real cone_eps,
                                   const arp::edge_cov_tree<2> *cov = nullptr,
                                   const edge_quad &quad = trapezoid_edge_quad()) {
  paper_tp_eval out;
  const edge_table E = build_edges(rod);
  out.grad.assign(E.n_vert, vec3::Zero());
  if (E.n == 0)
    return out;

  const real eps = std::max(l0_abs, real(1e-12));
  const real mass_scale = real(1);
  std::vector<real> ell_b, ep0_b, ep1_b;
  bind_edge_lengths(E, mass_scale, ell_b, ep0_b, ep1_b);
  const std::vector<vec3> T_b = E.T;
  const std::vector<vec3> &x = rod.x();
  const std::vector<index_t> ids = edge_ordinals(E.n);
  std::vector<std::vector<energy_scatter>> hits(E.n);

  auto bind = [&ell_b, &T_b, &ep0_b, &ep1_b](const std::vector<index_t> &edge_ids,
                                             calder::cov_rod_sum &sum) {
    sum.bind(calder::scalar_datum::create(edge_ids, ell_b));
    sum.bind(calder::vec3_datum::create(edge_ids, T_b));
    sum.bind(calder::scalar_datum::create(edge_ids, ep0_b));
    sum.bind(calder::scalar_datum::create(edge_ids, ep1_b));
  };
  auto kernel = [&](const index_t i, const index_t j, const vec3 &pi,
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
    energy_scatter hit;
    if (node_type == calder::LEAF) {
      if (leaf_edges_adjacent_or_self(E, ii, j0, j1))
        return vec3::Zero();
      if (j0 < 0 || j1 < 0)
        return vec3::Zero();
      const size_t a0 = static_cast<size_t>(E.ep0[ii]);
      const size_t a1 = static_cast<size_t>(E.ep1[ii]);
      const size_t b0 = static_cast<size_t>(j0);
      const size_t b1 = static_cast<size_t>(j1);
      if (a0 >= x.size() || a1 >= x.size() || b0 >= x.size() || b1 >= x.size())
        return vec3::Zero();
      const quad_term term = energy_quad_term(
          x[a0], x[a1], E.T[ii], E.ell[ii], x[b0], x[b1], Tj, lj, eps, p, quad);
      if (!(term.e > real(0)) || !std::isfinite(term.e))
        return vec3::Zero();
      const real half = real(0.5);
      hit.e = half * term.e;
      hit.n = 4;
      hit.v[0] = E.ep0[ii];
      hit.v[1] = E.ep1[ii];
      hit.v[2] = j0;
      hit.v[3] = j1;
      hit.g0 = half * term.ga0;
      hit.g1 = half * term.ga1;
      hit.g2 = half * term.gb0;
      hit.g3 = half * term.gb1;
    } else {
      const far_term term = energy_far_term(E.T[ii], E.ell[ii], pi, pj, Tj, lj,
                                            eps, p);
      if (!(term.e > real(0)) || !std::isfinite(term.e))
        return vec3::Zero();
      hit.e = term.e;
      hit.n = 2;
      hit.v[0] = E.ep0[ii];
      hit.v[1] = E.ep1[ii];
      hit.g0 = term.ga0;
      hit.g1 = term.ga1;
    }
    hits[ii].push_back(hit);
    return vec3::Zero();
  };

  arp::edge_cov_tree<2>::ptr owned;
  const arp::edge_cov_tree<2> *tree = cov;
  if (!tree) {
    owned = calder::build_rod_cov_tree(rod);
    tree = owned.get();
  }
  if (!tree)
    return out;
  calder::integrate_over_cov_tree<vec3>(*tree, ids, E.mid, bind, kernel, kernel,
                                        bh_eps, 1, cone_eps);

  auto add = [&](int vi, const vec3 &g) {
    if (vi < 0)
      return;
    const size_t s = static_cast<size_t>(vi);
    if (s < out.grad.size() && g.allFinite())
      out.grad[s] += g;
  };
  for (const auto &row : hits) {
    for (const energy_scatter &hit : row) {
      out.energy += hit.e;
      add(hit.v[0], hit.g0);
      add(hit.v[1], hit.g1);
      if (hit.n > 2) {
        add(hit.v[2], hit.g2);
        add(hit.v[3], hit.g3);
      }
    }
  }
  return out;
}

} // namespace sobolev_detail

/// Jacobi smooth Ag≈f for n_iters (in-place on g).
inline void sobolev_jacobi_smooth(
    asawa::rod::rod &rod, const std::vector<vec3> &f,
    const std::vector<real> &m_pov, const std::vector<real> &diag,
    const sobolev_tangent_point_force_config &cfg, real l0_abs, int n_iters,
    std::vector<vec3> &g, const arp::edge_cov_tree<2> *cov = nullptr,
    const sobolev_detail::cov_weight_rows *wB = nullptr,
    const sobolev_detail::cov_weight_rows *wB0 = nullptr) {
  const size_t n = f.size();
  if (n == 0 || n_iters <= 0)
    return;
  const real omega = real(0.85);
  const real bh = cfg.bh_eps;
  for (int it = 0; it < n_iters; ++it) {
    std::vector<vec3> Ag;
    sobolev_detail::matvec_A(rod, m_pov, cfg.alpha_B, cfg.alpha_0, cfg.p,
                             l0_abs, cfg.mass, bh, g, Ag, cov, cfg.cone_eps, wB,
                             wB0);
    for (size_t i = 0; i < n; ++i)
      g[i] += omega * (f[i] - Ag[i]) / diag[i];
  }
}

/// Solve Ag ≈ f by Jacobi; returns g (same space as f).
inline std::vector<vec3>
sobolev_jacobi_solve(asawa::rod::rod &rod, const std::vector<vec3> &f,
                     const sobolev_tangent_point_force_config &cfg,
                     real l0_abs,
                     arp::edge_cov_tree<2>::ptr cov = {}) {
  const size_t n = f.size();
  std::vector<vec3> g(n, vec3::Zero());
  if (n == 0 || cfg.jacobi_iters <= 0)
    return f;

  if (!cov)
    cov = calder::build_rod_cov_tree(rod);
  auto m_pov = sobolev_detail::sample_masses(rod);
  for (real &m : m_pov)
    m *= cfg.mass;
  const auto E = sobolev_detail::build_edges(rod);
  sobolev_detail::cov_weight_rows wB, wB0;
  if (cfg.alpha_B != real(0))
    wB = sobolev_detail::build_cov_weights(
        rod, E, cfg.p, l0_abs, cfg.mass, cfg.bh_eps, cfg.cone_eps,
        /*use_tp_k=*/false, *cov);
  if (cfg.alpha_0 != real(0))
    wB0 = sobolev_detail::build_cov_weights(
        rod, E, cfg.p, l0_abs, cfg.mass, cfg.bh_eps, cfg.cone_eps,
        /*use_tp_k=*/true, *cov);
  const auto diag = sobolev_detail::jacobi_diag(
      rod, m_pov, cfg.alpha_B, cfg.alpha_0, cfg.p, l0_abs, cfg.mass, cfg.bh_eps,
      cov.get(), cfg.cone_eps, cfg.alpha_B != real(0) ? &wB : nullptr,
      cfg.alpha_0 != real(0) ? &wB0 : nullptr);
  sobolev_jacobi_smooth(rod, f, m_pov, diag, cfg, l0_abs, cfg.jacobi_iters, g,
                        cov.get(), cfg.alpha_B != real(0) ? &wB : nullptr,
                        cfg.alpha_0 != real(0) ? &wB0 : nullptr);
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
  std::vector<real> m_pov;
  std::vector<real> diag;
};

inline void vcycle(asawa::rod::residual_forest &F, std::vector<level_ops> &ops,
                   const sobolev_tangent_point_force_config &cfg, real l0_abs,
                   int ell, int n_levels, const std::vector<vec3> &b,
                   std::vector<vec3> &x) {
  auto &L = ops[static_cast<size_t>(ell)];
  const int smooth = std::max(1, cfg.jacobi_iters);
  const sobolev_detail::cov_weight_rows *wB =
      cfg.alpha_B != real(0) ? &L.wB : nullptr;
  const sobolev_detail::cov_weight_rows *wB0 =
      cfg.alpha_0 != real(0) ? &L.wB0 : nullptr;
  sobolev_jacobi_smooth(*L.rod, b, L.m_pov, L.diag, cfg, l0_abs, smooth, x,
                        L.cov.get(), wB, wB0);

  if (ell + 1 >= n_levels) {
    // Extra coarsest iters
    sobolev_jacobi_smooth(*L.rod, b, L.m_pov, L.diag, cfg, l0_abs, smooth, x,
                          L.cov.get(), wB, wB0);
    return;
  }

  std::vector<vec3> Ax;
  sobolev_detail::matvec_A(*L.rod, L.m_pov, cfg.alpha_B, cfg.alpha_0, cfg.p,
                           l0_abs, cfg.mass, cfg.bh_eps, x, Ax, L.cov.get(),
                           cfg.cone_eps, wB, wB0);
  std::vector<vec3> r(b.size());
  for (size_t i = 0; i < b.size(); ++i)
    r[i] = b[i] - Ax[i];

  std::vector<vec3> b_c = F.restrict_JT(ell, r);
  std::vector<vec3> e_c(b_c.size(), vec3::Zero());
  vcycle(F, ops, cfg, l0_abs, ell + 1, n_levels, b_c, e_c);

  std::vector<vec3> e_f = F.prolong_J(ell, e_c);
  for (size_t i = 0; i < x.size() && i < e_f.size(); ++i)
    x[i] += e_f[i];

  sobolev_jacobi_smooth(*L.rod, b, L.m_pov, L.diag, cfg, l0_abs, smooth, x,
                        L.cov.get(), wB, wB0);
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
                 arp::edge_cov_tree<2>::ptr fine = {}) {
  const size_t n = f.size();
  if (n == 0)
    return f;
  if (cfg.num_v_cycles <= 0 || cfg.max_levels == 1)
    return sobolev_jacobi_solve(rod, f, cfg, l0_abs, std::move(fine));

  // Non-owning wrap: forest is local to this solve.
  asawa::rod::rod::ptr Rptr(&rod, [](asawa::rod::rod *) {});
  auto F = asawa::rod::residual_forest::create(Rptr);
  if (!F || !F->valid() || F->n_levels() < 2)
    return sobolev_jacobi_solve(rod, f, cfg, l0_abs, std::move(fine));

  int n_levels = F->n_levels();
  if (cfg.max_levels > 0)
    n_levels = std::min(n_levels, cfg.max_levels);
  if (n_levels < 2)
    return sobolev_jacobi_solve(rod, f, cfg, l0_abs, std::move(fine));

  // Sync coarse geometry from fine survivors; freeze for this solve.
  for (int ell = 0; ell + 1 < n_levels; ++ell)
    F->restrict_up(ell);

  std::vector<sobolev_mg_detail::level_ops> ops(
      static_cast<size_t>(n_levels));
  for (int ell = 0; ell < n_levels; ++ell) {
    auto &op = ops[static_cast<size_t>(ell)];
    op.rod = F->level_rod(ell);
    op.cov = sobolev_level_cov(*op.rod, ell, fine);
    op.m_pov = sobolev_detail::sample_masses(*op.rod);
    for (real &m : op.m_pov)
      m *= cfg.mass;
    const auto E = sobolev_detail::build_edges(*op.rod);
    if (cfg.alpha_B != real(0))
      op.wB = sobolev_detail::build_cov_weights(
          *op.rod, E, cfg.p, l0_abs, cfg.mass, cfg.bh_eps, cfg.cone_eps,
          /*use_tp_k=*/false, *op.cov);
    if (cfg.alpha_0 != real(0))
      op.wB0 = sobolev_detail::build_cov_weights(
          *op.rod, E, cfg.p, l0_abs, cfg.mass, cfg.bh_eps, cfg.cone_eps,
          /*use_tp_k=*/true, *op.cov);
    op.diag = sobolev_detail::jacobi_diag(
        *op.rod, op.m_pov, cfg.alpha_B, cfg.alpha_0, cfg.p, l0_abs, cfg.mass,
        cfg.bh_eps, op.cov.get(), cfg.cone_eps,
        cfg.alpha_B != real(0) ? &op.wB : nullptr,
        cfg.alpha_0 != real(0) ? &op.wB0 : nullptr);
  }

  std::vector<vec3> g(n, vec3::Zero());
  const int cycles = std::max(1, cfg.num_v_cycles);
  for (int c = 0; c < cycles; ++c)
    sobolev_mg_detail::vcycle(*F, ops, cfg, l0_abs, 0, n_levels, f, g);

  return g;
}

/// dÊ then Ag = dÊ. The fine tree is reused only on level 0. Each coarser
/// forest rod gets its own covariance tree for the diagonal and every matvec.
inline std::vector<vec3>
sobolev_tp_gradient(asawa::rod::rod &rod, const asawa::rod::dynamic &dynamic,
                    const sobolev_tangent_point_force_config &cfg) {
  const real l0 = cfg.l0 * dynamic._Cc;
  auto cov = calder::build_rod_cov_tree(rod);
  auto ev = sobolev_detail::eval_paper_tp(rod, l0, cfg.p, cfg.bh_eps,
                                          cfg.cone_eps, cov.get());
  std::vector<vec3> f = std::move(ev.grad);
  if (cfg.w != real(1.0)) {
    for (auto &a : f)
      a *= cfg.w;
  }
  return sobolev_mg_solve(rod, f, cfg, l0, std::move(cov));
}

inline real sobolev_quadratic_B(asawa::rod::rod &rod,
                                const std::vector<real> &m_pov, real p,
                                real l0_abs, real mass_scale, real bh_eps,
                                const std::vector<vec3> &u,
                                const std::vector<vec3> &v) {
  std::vector<vec3> Bv;
  sobolev_detail::matvec_B(rod, v, m_pov, p, l0_abs, mass_scale, bh_eps, Bv);
  real s = 0;
  for (size_t i = 0; i < u.size(); ++i)
    s += u[i].dot(Bv[i]);
  return s;
}

inline real sobolev_quadratic_B0(asawa::rod::rod &rod,
                                 const std::vector<real> &m_pov, real p,
                                 real l0_abs, real mass_scale, real bh_eps,
                                 const std::vector<vec3> &u,
                                 const std::vector<vec3> &v) {
  std::vector<vec3> B0v;
  sobolev_detail::matvec_B0(rod, v, m_pov, p, l0_abs, mass_scale, bh_eps, B0v);
  real s = 0;
  for (size_t i = 0; i < u.size(); ++i)
    s += u[i].dot(B0v[i]);
  return s;
}

} // namespace duchamp
} // namespace gaudi

#endif // GAUDI_DUCHAMP_TP_SOBOLEV_JACOBI_HPP
