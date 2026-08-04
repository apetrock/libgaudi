#ifndef GAUDI_DUCHAMP_MODULES_CFL_NEIGHBORHOOD_HPP
#define GAUDI_DUCHAMP_MODULES_CFL_NEIGHBORHOOD_HPP

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "gaudi/arp/hash_tree.hpp"
#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/common.h"
#include "gaudi/geometry_types.hpp"

namespace gaudi {
namespace duchamp {

/// Aggregation of per-site CFL votes.
enum class cfl_agg_mode { hard_min, softmin };

/// Which term produced the limiting scale.
enum class cfl_limit_term { none, edge, env };

/// Displacement-envelope CFL / soft self-intersection guard:
/// per-edge truncation + swept-AABB overlap + relative gap clamp.
struct cfl_neighborhood_config {
  real cfl_frac = 0.4;
  /// hard_min avoids softmin tau·log(n) collapse under many gated votes.
  cfl_agg_mode agg = cfl_agg_mode::hard_min;
  real softmin_tau = 0.0; ///< ≤0 ⇒ 0.05·hardmin(dti) when softmin is used
  /// Drop sites with |u| < gate_frac · max_disp before aggregate.
  real gate_frac = 0.1;
  real u_abs_floor = 0.0;
  /// Enable swept-AABB displacement-envelope (broad-phase overlap).
  bool use_envelope = true;
  real h_ref = 0.0; ///< predicted u = h_ref² f
  real dt_max = 0.0;
  real dt_min = 1.0e-8;
  real eps = 1.0e-16;
  /// Tube clearance: |R|_eff = |R| - radius; also inflates swept boxes.
  real radius = 0.0;
  /// Skip contributions with |R| below this × lavg (self / coincident).
  real self_skip_frac = 1.0e-3;
  /// Skip same-strand corners within this many hops (envelope is for
  /// cross-strand / spatial contact; along-strand is edge truncate).
  int chain_skip_hops = 8;
};

struct cfl_neighborhood_result {
  real dt = 0.0;
  real dt_pre_clamp = 0.0;
  real max_disp = 0.0;
  real peak_f = 0.0;
  real p50_f = 0.0;
  real lavg = 0.0;
  real u_over_lavg = 0.0;
  real s_edge_min = 0.0;
  real s_env_min = 0.0;
  real s_used_min = 0.0;
  real min_dti = 0.0;
  real softmin_tau = 0.0;
  real closing_max = 0.0;
  real gap_min = 0.0; ///< smallest |R|_eff among closing env samples
  /// Limiting envelope pair (narrow-phase), for stall diagnosis.
  int env_i = -1;
  int env_j = -1;
  real env_R = 0.0;
  real env_reff = 0.0;
  real env_closing = 0.0;
  real env_s = 0.0;
  int n_gated = 0;
  int n_total = 0;
  int n_edge_limited = 0;
  int n_env_limited = 0;
  int n_env_overlap = 0;
  int n_env_closing = 0;
  int n_env_hits = 0; ///< broad-phase overlap hits that reached narrow phase
  cfl_limit_term limiting = cfl_limit_term::none;
  bool used_peak_scaled = false;
  bool used_envelope = false;
  bool floored = false; ///< dt_pre_clamp < dt_min before clamp
};

inline const char *cfl_limit_term_name(cfl_limit_term t) {
  switch (t) {
  case cfl_limit_term::edge:
    return "edge";
  case cfl_limit_term::env:
    return "env";
  default:
    return "none";
  }
}

inline real cfl_softmin(const std::vector<real> &vals, real tau) {
  if (vals.empty())
    return std::numeric_limits<real>::infinity();
  if (!(tau > 0.0) || vals.size() == 1) {
    real m = vals[0];
    for (real v : vals)
      m = std::min(m, v);
    return m;
  }
  real m = vals[0];
  for (real v : vals)
    m = std::min(m, v);
  real acc = 0.0;
  for (real v : vals)
    acc += std::exp(-(v - m) / tau);
  return m - tau * std::log(acc);
}

inline real cfl_hardmin(const std::vector<real> &vals) {
  if (vals.empty())
    return std::numeric_limits<real>::infinity();
  real m = vals[0];
  for (real v : vals)
    m = std::min(m, v);
  return m;
}

/// Relative closing scale: s = cfl_frac · |R|_eff / max(−U·R̂, 0).
/// |R|_eff = max(|R| - radius, 0). Same geometry for edges and envelope.
inline real cfl_rel_scale(const vec3 &R, const vec3 &U, real cfl_frac, real eps,
                          real radius = 0.0, real *closing_out = nullptr,
                          real *reff_out = nullptr) {
  const real ell = R.norm();
  if (!(ell > eps)) {
    if (closing_out)
      *closing_out = 0.0;
    if (reff_out)
      *reff_out = 0.0;
    return std::numeric_limits<real>::infinity();
  }
  const vec3 rhat = R / ell;
  const real closing = -(U.dot(rhat));
  if (closing_out)
    *closing_out = std::max(closing, real(0.0));
  const real reff = ell - std::max(radius, real(0.0));
  if (reff_out)
    *reff_out = reff;
  if (!(closing > eps))
    return std::numeric_limits<real>::infinity();
  // Already penetrating / contacting: not a predicted-contact CFL vote
  // (s = cfl_frac·eps/closing crushed dti to dt_min and stalled the sim).
  // Collision constraints own the overlap; envelope only limits approach.
  if (!(reff > eps))
    return std::numeric_limits<real>::infinity();
  return cfl_frac * reff / closing;
}

inline real cfl_edge_scale(const vec3 &R, const vec3 &U, real cfl_frac,
                           real eps, real *closing_out = nullptr) {
  return cfl_rel_scale(R, U, cfl_frac, eps, /*radius=*/0.0, closing_out,
                       nullptr);
}

inline real cfl_scale_to_dt(real s, real h_ref, real eps) {
  if (!(s > 0.0) || !(h_ref > 0.0))
    return 0.0;
  if (!std::isfinite(s))
    return std::numeric_limits<real>::infinity();
  return h_ref * std::sqrt(std::max(s, eps));
}

inline bool cfl_passes_disp_gate(real u_norm, real max_disp,
                                 const cfl_neighborhood_config &cfg) {
  if (u_norm < cfg.u_abs_floor)
    return false;
  if (cfg.gate_frac > 0.0 && max_disp > cfg.eps &&
      u_norm < cfg.gate_frac * max_disp)
    return false;
  return true;
}

inline void cfl_force_stats(const std::vector<vec3> &forces, real &peak_f,
                            real &p50_f) {
  peak_f = 0.0;
  p50_f = 0.0;
  if (forces.empty())
    return;
  std::vector<real> norms;
  norms.reserve(forces.size());
  for (const vec3 &f : forces) {
    const real n = f.norm();
    if (!std::isfinite(n))
      continue;
    peak_f = std::max(peak_f, n);
    norms.push_back(n);
  }
  if (norms.empty())
    return;
  std::nth_element(norms.begin(), norms.begin() + norms.size() / 2,
                   norms.end());
  p50_f = norms[norms.size() / 2];
}

struct cfl_envelope_pass {
  std::vector<real> s;
  real gap_min = std::numeric_limits<real>::infinity();
  int n_overlap = 0;
  int n_closing = 0;
  int n_hits = 0;
  int lim_i = -1;
  int lim_j = -1;
  real lim_R = 0.0;
  real lim_reff = 0.0;
  real lim_closing = 0.0;
  real lim_s = std::numeric_limits<real>::infinity();
};

/// Per-corner strand id + arc index along that strand (for chain skips).
inline void cfl_rod_strand_arcs(const asawa::rod::rod &rod,
                                std::vector<int> &strand, std::vector<int> &arc,
                                std::vector<int> &strand_len) {
  const size_t n = rod.x().size();
  strand.assign(n, -1);
  arc.assign(n, -1);
  strand_len.clear();
  std::vector<char> visited(n, 0);
  int sid = 0;
  for (size_t seed = 0; seed < n; ++seed) {
    if (visited[seed])
      continue;
    const asawa::rod::CornerId cseed =
        asawa::rod::corner_id(static_cast<int>(seed));
    if (rod.next(cseed) < asawa::rod::corner_id(0) &&
        rod.prev(cseed) < asawa::rod::corner_id(0)) {
      visited[seed] = 1;
      strand[seed] = sid;
      arc[seed] = 0;
      strand_len.push_back(1);
      ++sid;
      continue;
    }

    asawa::rod::CornerId start = cseed;
    asawa::rod::CornerId s = start;
    for (;;) {
      const asawa::rod::CornerId p = rod.prev(s);
      if (p < asawa::rod::corner_id(0) || p == start)
        break;
      s = p;
    }

    int a = 0;
    asawa::rod::CornerId i = s;
    do {
      const size_t ii = static_cast<size_t>(static_cast<int>(i));
      visited[ii] = 1;
      strand[ii] = sid;
      arc[ii] = a++;
      const asawa::rod::CornerId j = rod.next(i);
      if (j < asawa::rod::corner_id(0))
        break;
      i = j;
    } while (i != s && !visited[static_cast<size_t>(static_cast<int>(i))]);
    strand_len.push_back(a);
    ++sid;
  }
}

inline bool cfl_same_strand_near(int i, int k, int hops,
                                 const std::vector<int> &strand,
                                 const std::vector<int> &arc,
                                 const std::vector<int> &strand_len) {
  if (i < 0 || k < 0 || hops < 0)
    return false;
  if (static_cast<size_t>(i) >= strand.size() ||
      static_cast<size_t>(k) >= strand.size())
    return false;
  if (strand[i] < 0 || strand[i] != strand[k])
    return false;
  const int len = strand_len[static_cast<size_t>(strand[i])];
  int d = std::abs(arc[i] - arc[k]);
  if (len > 0)
    d = std::min(d, len - d); // closed loops
  return d <= hops;
}

/// Swept AABB of point x under displacement u, inflated by radius.
inline ext::extents_t cfl_swept_box(const vec3 &x, const vec3 &u, real radius) {
  ext::extents_t e = ext::init();
  e = ext::expand(e, x);
  e = ext::expand(e, x + u);
  if (radius > 0.0)
    e = ext::inflate(e, radius);
  return e;
}

/// BVH whose leaf extents are swept boxes (Morton order on box centers).
struct cfl_swept_bvh {
  std::vector<index_t> indices; ///< sorted leaf → original site
  std::vector<arp::radix_tree_node> internal_nodes;
  std::vector<arp::radix_tree_node> leaf_nodes;
  arp::TreeResult<ext::extents_t> bvh;
};

inline cfl_swept_bvh
cfl_build_swept_bvh(const std::vector<ext::extents_t> &boxes) {
  cfl_swept_bvh out;
  const size_t n = boxes.size();
  if (n == 0)
    return out;

  std::vector<vec3> centers(n);
  for (size_t i = 0; i < n; ++i)
    centers[i] = real(0.5) * (boxes[i][0] + boxes[i][1]);

  auto [hashes, indices, internal_nodes, leaf_nodes] =
      arp::make_hash_tree(centers);
  (void)hashes;
  if (indices.empty() || internal_nodes.empty())
    return out;

  std::vector<ext::extents_t> leaf_sorted(n);
  for (size_t i = 0; i < n; ++i)
    leaf_sorted[i] = boxes[static_cast<size_t>(indices[i])];

  auto internal = arp::build_pyramid(
      leaf_sorted, internal_nodes, leaf_nodes,
      [](const ext::extents_t &a, const ext::extents_t &b) {
        return ext::expand(b, a);
      },
      ext::init());

  out.indices = std::move(indices);
  out.internal_nodes = std::move(internal_nodes);
  out.leaf_nodes = std::move(leaf_nodes);
  out.bvh.leaf = std::move(leaf_sorted);
  out.bvh.internal = std::move(internal);
  return out;
}

/// Collect original site ids whose swept boxes overlap Q.
template <typename OnHit>
inline void cfl_query_box_overlaps(const cfl_swept_bvh &tree,
                                   const ext::extents_t &Q, OnHit &&on_hit) {
  if (tree.internal_nodes.empty() || tree.bvh.leaf.empty())
    return;
  arp::traverse_bfs(
      tree.internal_nodes, tree.leaf_nodes,
      [](index_t, const arp::radix_tree_node &) {},
      [&](index_t, index_t leaf_id, const arp::radix_tree_node &) {
        if (leaf_id < 0 || static_cast<size_t>(leaf_id) >= tree.bvh.leaf.size())
          return;
        if (!ext::overlap(Q, tree.bvh.leaf[static_cast<size_t>(leaf_id)]))
          return;
        on_hit(tree.indices[static_cast<size_t>(leaf_id)]);
      },
      [&](index_t node_id, const arp::radix_tree_node &) -> bool {
        if (node_id < 0 ||
            static_cast<size_t>(node_id) >= tree.bvh.internal.size())
          return false;
        return ext::overlap(Q, tree.bvh.internal[static_cast<size_t>(node_id)]);
      });
}

/// Point swept-AABB envelope: build BVH once, overlap-query each site.
/// `skip_pair(i, j)` drops topology-local neighbors edge truncate owns.
template <typename SkipPair>
inline cfl_envelope_pass
cfl_point_envelope_scales(const std::vector<vec3> &x, const std::vector<vec3> &u,
                          const cfl_neighborhood_config &cfg, real lavg,
                          SkipPair &&skip_pair) {
  cfl_envelope_pass pass;
  const size_t n = x.size();
  pass.s.assign(n, std::numeric_limits<real>::infinity());
  if (n == 0 || u.size() < n)
    return pass;

  const real skip = std::max(cfg.self_skip_frac * lavg, cfg.eps);
  const real rad = std::max(cfg.radius, real(0.0));

  std::vector<ext::extents_t> boxes(n);
  for (size_t i = 0; i < n; ++i)
    boxes[i] = cfl_swept_box(x[i], u[i], rad);

  const cfl_swept_bvh tree = cfl_build_swept_bvh(boxes);
  if (tree.internal_nodes.empty())
    return pass;

#pragma omp parallel for
  for (int ii = 0; ii < static_cast<int>(n); ++ii) {
    const size_t i = static_cast<size_t>(ii);
    real s_i = std::numeric_limits<real>::infinity();
    int best_j = -1;
    real best_R = 0.0, best_reff = 0.0, best_closing = 0.0;
    cfl_query_box_overlaps(tree, boxes[i], [&](index_t j_raw) {
      if (j_raw < 0)
        return;
      const size_t j = static_cast<size_t>(j_raw);
      if (j >= n || j == i)
        return;
      if (skip_pair(static_cast<int>(i), static_cast<int>(j)))
        return;

      const vec3 R = x[j] - x[i];
      const real g0 = R.norm();
      if (!(g0 > skip))
        return;

#pragma omp atomic
      pass.n_hits += 1;

      const vec3 U = u[j] - u[i];
      real closing = 0.0;
      real reff = 0.0;
      const real s = cfl_rel_scale(R, U, cfg.cfl_frac, cfg.eps, cfg.radius,
                                   &closing, &reff);
      if (closing > cfg.eps) {
#pragma omp atomic
        pass.n_closing += 1;
        if (!(reff > cfg.eps)) {
#pragma omp atomic
          pass.n_overlap += 1;
        }
        if (reff < pass.gap_min)
          pass.gap_min = reff;
      }
      if (std::isfinite(s) && s < s_i) {
        s_i = s;
        best_j = static_cast<int>(j);
        best_R = g0;
        best_reff = reff;
        best_closing = closing;
      }
    });
    pass.s[i] = s_i;
    if (std::isfinite(s_i) && s_i < pass.lim_s) {
#pragma omp critical(cfl_env_lim)
      {
        if (s_i < pass.lim_s) {
          pass.lim_s = s_i;
          pass.lim_i = static_cast<int>(i);
          pass.lim_j = best_j;
          pass.lim_R = best_R;
          pass.lim_reff = best_reff;
          pass.lim_closing = best_closing;
        }
      }
    }
  }

  if (!std::isfinite(pass.gap_min))
    pass.gap_min = 0.0;
  return pass;
}

/// Rod envelope: corner swept boxes + same-strand hop skip.
inline cfl_envelope_pass
cfl_point_envelope_scales(const asawa::rod::rod &rod, const std::vector<vec3> &u,
                          const cfl_neighborhood_config &cfg, real lavg) {
  std::vector<int> strand, arc, strand_len;
  cfl_rod_strand_arcs(rod, strand, arc, strand_len);
  const int hops = std::max(cfg.chain_skip_hops, 1);
  return cfl_point_envelope_scales(
      rod.x(), u, cfg, lavg,
      [&](int i, int j) {
        return cfl_same_strand_near(i, j, hops, strand, arc, strand_len);
      });
}

/// Shell envelope: vertex swept boxes; skip is a no-op (self already dropped).
inline cfl_envelope_pass
cfl_point_envelope_scales(const asawa::shell::shell &M,
                          const std::vector<vec3> &u,
                          const cfl_neighborhood_config &cfg, real lavg) {
  const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);
  return cfl_point_envelope_scales(x, u, cfg, lavg,
                                   [](int, int) { return false; });
}

struct cfl_safe_scale_result {
  real alpha = 1.0;
  real alpha_peak = 1.0;
  real alpha_geom = 1.0;
  cfl_neighborhood_result detail;
};

/// Rod CFL from proposed per-site displacement u (edge + envelope).
inline cfl_neighborhood_result
cfl_estimate_rod_from_u(const asawa::rod::rod &rod,
                        const std::vector<vec3> &u,
                        const cfl_neighborhood_config &cfg) {
  cfl_neighborhood_result out;
  out.used_envelope = cfg.use_envelope;
  const auto &x = rod.x();
  const size_t n = x.size();
  if (n == 0 || u.size() != n || !(cfg.h_ref > 0.0)) {
    out.dt = cfg.dt_max > 0.0 ? cfg.dt_max : 0.0;
    out.dt_pre_clamp = out.dt;
    return out;
  }

  for (const vec3 &q : x) {
    if (!q.array().isFinite().all()) {
      out.dt = cfg.dt_max > 0.0 ? cfg.dt_max : 0.0;
      out.dt_pre_clamp = out.dt;
      return out;
    }
  }

  const real inv_h2 = real(1.0) / std::max(cfg.h_ref * cfg.h_ref, cfg.eps);
  std::vector<vec3> pseudo_f(n);
  real max_disp = 0.0;
  for (size_t i = 0; i < n; ++i) {
    if (!u[i].array().isFinite().all()) {
      out.dt = cfg.dt_max > 0.0 ? cfg.dt_max : 0.0;
      out.dt_pre_clamp = out.dt;
      return out;
    }
    pseudo_f[i] = u[i] * inv_h2;
    max_disp = std::max(max_disp, u[i].norm());
  }
  cfl_force_stats(pseudo_f, out.peak_f, out.p50_f);
  out.max_disp = max_disp;

  real lavg_est = 0.0;
  int lavg_n = 0;
  for (size_t i = 0; i < n; ++i) {
    const asawa::rod::CornerId ci = asawa::rod::corner_id(static_cast<int>(i));
    const asawa::rod::CornerId cn = rod.next(ci);
    if (cn < asawa::rod::corner_id(0))
      continue;
    const size_t j = static_cast<size_t>(static_cast<int>(cn));
    if (j >= n)
      continue;
    lavg_est += (x[j] - x[i]).norm();
    lavg_n += 1;
  }
  if (lavg_n > 0)
    lavg_est /= real(lavg_n);
  lavg_est = std::max(lavg_est, real(1e-6));
  out.lavg = lavg_est;
  out.u_over_lavg = max_disp / lavg_est;

  cfl_envelope_pass env_pass;
  if (cfg.use_envelope)
    env_pass = cfl_point_envelope_scales(rod, u, cfg, lavg_est);
  out.n_env_overlap = env_pass.n_overlap;
  out.n_env_closing = env_pass.n_closing;
  out.n_env_hits = env_pass.n_hits;
  out.gap_min = env_pass.gap_min;
  out.env_i = env_pass.lim_i;
  out.env_j = env_pass.lim_j;
  out.env_R = env_pass.lim_R;
  out.env_reff = env_pass.lim_reff;
  out.env_closing = env_pass.lim_closing;
  out.env_s = std::isfinite(env_pass.lim_s) ? env_pass.lim_s : 0.0;

  std::vector<real> dts;
  dts.reserve(n);
  real s_used_min = std::numeric_limits<real>::infinity();
  real s_edge_min = std::numeric_limits<real>::infinity();
  real s_env_min = std::numeric_limits<real>::infinity();
  real closing_max = 0.0;
  cfl_limit_term limiting = cfl_limit_term::none;

  for (size_t i = 0; i < n; ++i) {
    const asawa::rod::CornerId ci = asawa::rod::corner_id(static_cast<int>(i));
    const asawa::rod::CornerId cn = rod.next(ci);
    if (cn < asawa::rod::corner_id(0))
      continue;

    out.n_total += 1;
    const size_t j = static_cast<size_t>(static_cast<int>(cn));
    if (j >= n)
      continue;

    const real u_gate = std::max(u[i].norm(), u[j].norm());
    if (!cfl_passes_disp_gate(u_gate, max_disp, cfg))
      continue;

    const vec3 R = x[j] - x[i];
    const vec3 U = u[j] - u[i];
    const real ell = R.norm();
    if (!(ell > cfg.eps))
      continue;

    real closing = 0.0;
    const real s_edge =
        cfl_edge_scale(R, U, cfg.cfl_frac, cfg.eps, &closing);
    closing_max = std::max(closing_max, closing);
    if (std::isfinite(s_edge))
      s_edge_min = std::min(s_edge_min, s_edge);

    real s_e = std::numeric_limits<real>::infinity();
    if (cfg.use_envelope && env_pass.s.size() == n) {
      s_e = std::min(env_pass.s[i], env_pass.s[j]);
      if (std::isfinite(s_e))
        s_env_min = std::min(s_env_min, s_e);
    }

    real s = s_edge;
    cfl_limit_term term = cfl_limit_term::edge;
    if (cfg.use_envelope && s_e < s) {
      s = s_e;
      term = cfl_limit_term::env;
    }

    if (!std::isfinite(s) || !(s > 0.0))
      continue;

    const real dti = cfl_scale_to_dt(s, cfg.h_ref, cfg.eps);
    if (!std::isfinite(dti) || !(dti > 0.0))
      continue;

    dts.push_back(dti);
    if (s < s_used_min) {
      s_used_min = s;
      limiting = term;
    }
    if (term == cfl_limit_term::edge)
      out.n_edge_limited += 1;
    else if (term == cfl_limit_term::env)
      out.n_env_limited += 1;
    out.n_gated += 1;
  }

  out.s_edge_min = std::isfinite(s_edge_min) ? s_edge_min : 0.0;
  out.s_env_min = std::isfinite(s_env_min) ? s_env_min : 0.0;
  out.s_used_min = std::isfinite(s_used_min) ? s_used_min : 0.0;
  out.closing_max = closing_max;
  out.limiting = limiting;
  out.min_dti = dts.empty() ? 0.0 : cfl_hardmin(dts);

  real dt;
  if (dts.empty()) {
    dt = cfg.dt_max > 0.0 ? cfg.dt_max : 0.0;
    out.softmin_tau = 0.0;
  } else if (cfg.agg == cfl_agg_mode::hard_min) {
    dt = cfl_hardmin(dts);
    out.softmin_tau = 0.0;
  } else {
    // Softmin undershoots via -τ log n when τ is large; keep τ on the
    // scale of the hard min so it cannot collapse to dt_min.
    const real hard = cfl_hardmin(dts);
    real tau = cfg.softmin_tau;
    if (!(tau > 0.0))
      tau = 0.05 * hard;
    out.softmin_tau = tau;
    dt = cfl_softmin(dts, tau);
    if (dt < hard)
      dt = hard;
  }

  out.dt_pre_clamp = dt;
  out.floored = (cfg.dt_min > 0.0 && std::isfinite(dt) && dt < cfg.dt_min);
  if (cfg.dt_max > 0.0)
    dt = std::min(dt, cfg.dt_max);
  if (cfg.dt_min > 0.0 && std::isfinite(dt))
    dt = std::max(dt, cfg.dt_min);
  if (!std::isfinite(dt) || !(dt > 0.0))
    dt = cfg.dt_max > 0.0 ? cfg.dt_max : cfg.dt_min;

  out.dt = dt;
  return out;
}

/// Peak + geometric scale for a proposed displacement field δ.
/// δ_safe = α·δ; equivalent substep dti ≈ α·h_ref.
inline cfl_safe_scale_result
cfl_safe_scale_rod(const asawa::rod::rod &rod, const std::vector<vec3> &delta,
                   const cfl_neighborhood_config &cfg, real dx_max = 0.0) {
  cfl_safe_scale_result out;
  out.detail = cfl_estimate_rod_from_u(rod, delta, cfg);

  if (cfg.h_ref > 0.0 && out.detail.dt_pre_clamp > 0.0)
    out.alpha_geom =
        std::min(real(1.0), out.detail.dt_pre_clamp / cfg.h_ref);
  else
    out.alpha_geom = 1.0;

  if (dx_max > 0.0 && out.detail.max_disp > cfg.eps)
    out.alpha_peak =
        std::min(real(1.0), dx_max / out.detail.max_disp);
  else
    out.alpha_peak = 1.0;

  out.alpha = std::min(out.alpha_peak, out.alpha_geom);
  return out;
}

/// Rod CFL: per-edge truncation + swept-AABB relative envelope.
inline cfl_neighborhood_result
cfl_estimate_rod_dt(const asawa::rod::rod &rod, const std::vector<vec3> &forces,
                    const cfl_neighborhood_config &cfg) {
  const auto &x = rod.x();
  const size_t n = x.size();
  if (n == 0 || forces.size() < n || !(cfg.h_ref > 0.0)) {
    cfl_neighborhood_result out;
    out.dt = cfg.dt_max > 0.0 ? cfg.dt_max : 0.0;
    out.dt_pre_clamp = out.dt;
    return out;
  }

  std::vector<vec3> u(n);
  for (size_t i = 0; i < n; ++i)
    u[i] = (cfg.h_ref * cfg.h_ref) * forces[i];
  return cfl_estimate_rod_from_u(rod, u, cfg);
}

} // namespace duchamp
} // namespace gaudi

#endif
