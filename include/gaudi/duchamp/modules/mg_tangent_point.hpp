#ifndef GAUDI_DUCHAMP_MG_TANGENT_POINT_HPP
#define GAUDI_DUCHAMP_MG_TANGENT_POINT_HPP

/// Hierarchical TP on rod edge BVH (no topology forest).
/// Ascent: extract → make_bvh → Calder points TP → Kabsch R_sim → T←R T → prolong_x.
/// Descent: prolong_x; R_sim from before/after prolong for T when needed.

#include "gaudi/arp/arp.h"
#include "gaudi/arp/procrustes_pyramid.hpp"
#include "gaudi/arp/tree_view.hpp"
#include "gaudi/asawa/rod/dynamic.hpp"
#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/calder/point_tangent_point.hpp"
#include "gaudi/common.h"
#include "gaudi/duchamp/modules/tp_config.hpp"
#include "gaudi/duchamp/modules/tp_energies.hpp"
#include "gaudi/duchamp/modules/tp_line_search.hpp"

#include <Eigen/Eigenvalues>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>

namespace gaudi {
namespace duchamp {

struct mg_tp_config {
  bool enable = false;
  /// Coarsest simulated depth = max(0, max_L - max_sim_offset).
  int max_sim_offset = 4;
  int num_v_cycles = 1;
  /// GD iters per level on the extract set.
  int coarse_iters = 4;
  real coarse_step_frac = 0.25;
  real coarse_l0 = 0.5;
  real coarse_p = 3.0;
  real bh_eps = 0.25;
  bool refine_fine = true;
  /// Legacy (unused by new driver); kept so old cfg still compiles.
  real max_leaves = 32.0;
};

inline vec3 cov_principal_dir(const mat3 &C) {
  Eigen::SelfAdjointEigenSolver<mat3> es(C);
  if (es.info() != Eigen::Success)
    return vec3::UnitZ();
  return es.eigenvectors().col(2).normalized();
}

inline real mean_nn_spacing(const std::vector<vec3> &c) {
  if (c.size() < 2)
    return real(1);
  real sum = 0;
  size_t cnt = 0;
  for (size_t i = 0; i < c.size(); ++i) {
    real best = real(1e30);
    for (size_t j = 0; j < c.size(); ++j) {
      if (i == j)
        continue;
      best = std::min(best, (c[i] - c[j]).norm());
    }
    if (best < real(1e29)) {
      sum += best;
      ++cnt;
    }
  }
  return cnt ? sum / real(cnt) : real(1);
}

/// One level: Calder points TP + R_sim tangent transport. Updates x, T in place.
inline mat3 level_tp_step(std::vector<vec3> &x, std::vector<vec3> &T,
                          const std::vector<real> &w, const mg_tp_config &cfg) {
  if (x.size() < 2)
    return mat3::Identity();
  const std::vector<vec3> x0 = x;
  const real spacing = mean_nn_spacing(x);
  const real eps = std::max(cfg.coarse_l0 * spacing, real(1e-8));
  const real step0 = std::max(cfg.coarse_step_frac * spacing, real(1e-8));

  for (int it = 0; it < cfg.coarse_iters; ++it) {
    std::vector<vec3> g = calder::tangent_point_gradient_points_cauchy(
        x, T, w, eps, cfg.coarse_p, cfg.bh_eps);
    real gmax = 0;
    for (const auto &gi : g)
      gmax = std::max(gmax, gi.norm());
    if (!(gmax > real(1e-16)))
      break;
    const real scale = step0 / gmax;
    // g from Calder is ∇E (same as rod cauchy / tp_grad_energy); descend.
    for (size_t i = 0; i < x.size(); ++i)
      x[i] += -scale * g[i];
  }

  const mat3 R_sim = calder::kabsch_sim_rotation(x0, x);
  for (size_t i = 0; i < T.size(); ++i) {
    T[i] = R_sim * T[i];
    if (T[i].squaredNorm() > real(1e-24))
      T[i].normalize();
  }
  return R_sim;
}

struct mg_rod_pyramid {
  arp::T2::ptr edge_tree;
  arp::procrustes_pyramid P;
  std::vector<index_t> edge_vert_ids;
  std::vector<vec3> leaf_T; ///< sorted-leaf order tangents (edge N2)

  static mg_rod_pyramid build(asawa::rod::rod &rod) {
    mg_rod_pyramid out;
    const std::vector<vec3> &x = rod.x();
    out.edge_vert_ids = rod.get_edge_vert_ids();
    const size_t nE = out.edge_vert_ids.size() / 2;
    std::vector<vec3> mid(nE);
    for (size_t e = 0; e < nE; ++e) {
      const index_t v0 = out.edge_vert_ids[2 * e];
      const index_t v1 = out.edge_vert_ids[2 * e + 1];
      mid[e] = real(0.5) * (x[static_cast<size_t>(v0)] +
                            x[static_cast<size_t>(v1)]);
    }
    out.edge_tree = arp::T2::create(out.edge_vert_ids, x, 12);
    out.P.build(arp::make_tree_view(*out.edge_tree), mid);

    // Edge tangents from rod N2c, placed in sorted-leaf order.
    std::vector<vec3> Tc = rod.N2c();
    out.leaf_T.assign(out.P.leaf_count(), vec3::UnitZ());
    const arp::tree_view &view = out.P.view();
    for (size_t i = 0; i < out.P.leaf_count(); ++i) {
      const index_t orig = view.index(static_cast<index_t>(i));
      if (static_cast<size_t>(orig) < Tc.size()) {
        out.leaf_T[i] = Tc[static_cast<size_t>(orig)];
        if (out.leaf_T[i].squaredNorm() > real(1e-24))
          out.leaf_T[i].normalize();
      }
    }
    return out;
  }

  void sync_from_rod_positions(asawa::rod::rod &rod) {
    const std::vector<vec3> &x = rod.x();
    const size_t nE = edge_vert_ids.size() / 2;
    std::vector<vec3> mid(nE);
    for (size_t e = 0; e < nE; ++e) {
      const index_t v0 = edge_vert_ids[2 * e];
      const index_t v1 = edge_vert_ids[2 * e + 1];
      mid[e] = real(0.5) * (x[static_cast<size_t>(v0)] +
                            x[static_cast<size_t>(v1)]);
    }
    P.set_leaf_positions_original(mid);
    P.restrict_from_leaves();
  }

  /// Extract at depth-from-leaves \p lvl (lvl≥1). lvl==0 → all leaves.
  void extract_level(index_t lvl, std::vector<index_t> &nodes,
                     std::vector<vec3> &x, std::vector<vec3> &T,
                     std::vector<real> &w) const {
    nodes.clear();
    x.clear();
    T.clear();
    w.clear();
    if (lvl <= 0) {
      const auto &live = P.leaf_live_sorted();
      nodes.resize(live.size(), arp::UNULL);
      x = live;
      T = leaf_T;
      w.assign(live.size(), real(1));
      return;
    }
    nodes = P.extract_at_depth(lvl);
    const auto counts = P.subtree_leaf_counts();
    x.reserve(nodes.size());
    T.reserve(nodes.size());
    w.reserve(nodes.size());
    for (index_t n : nodes) {
      x.push_back(P.avg_live()[static_cast<size_t>(n)]);
      // Average unsigned tangents in subtree via structure tensor.
      mat3 C = mat3::Zero();
      real mass = 0;
      P.for_each_subtree_leaf(n, [&](index_t lid) {
        const vec3 &t = leaf_T[static_cast<size_t>(lid)];
        C += t * t.transpose();
        mass += real(1);
      });
      vec3 t = cov_principal_dir(C);
      // Sign vs first leaf tangent.
      bool signed_ok = false;
      P.for_each_subtree_leaf(n, [&](index_t lid) {
        if (signed_ok)
          return;
        if (t.dot(leaf_T[static_cast<size_t>(lid)]) < real(0))
          t = -t;
        signed_ok = true;
      });
      T.push_back(t);
      w.push_back(std::max(counts[static_cast<size_t>(n)], real(1)));
    }
  }

  /// Apply level centers x1: set poses, prolong_x; transport leaf_T by R_sim
  /// under each node (or all leaves if lvl==0).
  void apply_level_update(index_t lvl, const std::vector<index_t> &nodes,
                          const std::vector<vec3> &x0,
                          const std::vector<vec3> &x1, const mat3 &R_sim) {
    if (lvl <= 0) {
      // Write leaf positions in original order then restrict.
      std::vector<vec3> orig(x1.size(), vec3::Zero());
      const arp::tree_view &view = P.view();
      for (size_t i = 0; i < x1.size(); ++i) {
        const index_t o = view.index(static_cast<index_t>(i));
        if (static_cast<size_t>(o) < orig.size())
          orig[static_cast<size_t>(o)] = x1[i];
        else {
          // x1 is sorted-leaf sized; orig must cover max original index.
        }
      }
      // Resize orig to max original index.
      index_t omax = 0;
      for (size_t i = 0; i < x1.size(); ++i)
        omax = std::max(omax, view.index(static_cast<index_t>(i)));
      orig.assign(static_cast<size_t>(omax) + 1, vec3::Zero());
      for (size_t i = 0; i < x1.size(); ++i)
        orig[static_cast<size_t>(view.index(static_cast<index_t>(i)))] = x1[i];
      P.set_leaf_positions_original(orig);
      for (size_t i = 0; i < leaf_T.size(); ++i) {
        leaf_T[i] = R_sim * leaf_T[i];
        if (leaf_T[i].squaredNorm() > real(1e-24))
          leaf_T[i].normalize();
      }
      P.restrict_from_leaves();
      return;
    }
    for (size_t i = 0; i < nodes.size(); ++i) {
      const index_t n = nodes[i];
      const mat3 R0 = P.R_live()[static_cast<size_t>(n)];
      P.set_node_pose(n, x1[i], R0);
      P.prolong_from(n);
      P.for_each_subtree_leaf(n, [&](index_t lid) {
        leaf_T[static_cast<size_t>(lid)] =
            R_sim * leaf_T[static_cast<size_t>(lid)];
        if (leaf_T[static_cast<size_t>(lid)].squaredNorm() > real(1e-24))
          leaf_T[static_cast<size_t>(lid)].normalize();
      });
    }
    P.encode_from_leaves();
  }

  std::vector<vec3>
  leaf_mid_delta_to_vert_disp(const std::vector<vec3> &mid0) const {
    const auto mid1 = P.leaf_positions_original();
    index_t vmax = 0;
    for (index_t v : edge_vert_ids)
      vmax = std::max(vmax, v);
    std::vector<vec3> dvert(static_cast<size_t>(vmax) + 1, vec3::Zero());
    const size_t nE = std::min(mid0.size(), mid1.size());
    for (size_t e = 0; e < nE; ++e) {
      const vec3 dm = mid1[e] - mid0[e];
      const index_t v0 = edge_vert_ids[2 * e];
      const index_t v1 = edge_vert_ids[2 * e + 1];
      dvert[static_cast<size_t>(v0)] += real(0.5) * dm;
      dvert[static_cast<size_t>(v1)] += real(0.5) * dm;
    }
    return dvert;
  }
};

/// Run \p cfg.num_v_cycles of ascent (sim+restrict) then descent (prolong_x).
inline void mg_vcycle(mg_rod_pyramid &pyr, const mg_tp_config &cfg) {
  const index_t max_L = pyr.P.max_depth_from_leaves();
  const index_t max_sim_L = static_cast<index_t>(
      std::max(0, static_cast<int>(max_L) - cfg.max_sim_offset));

  for (int cycle = 0; cycle < cfg.num_v_cycles; ++cycle) {
    // Ascent: j = 0 .. max_sim_L
    for (index_t j = 0; j <= max_sim_L; ++j) {
      std::vector<index_t> nodes;
      std::vector<vec3> x, T;
      std::vector<real> w;
      pyr.extract_level(j, nodes, x, T, w);
      if (x.size() < 2)
        continue;
      const std::vector<vec3> x0 = x;
      const mat3 R_sim = level_tp_step(x, T, w, cfg);
      pyr.apply_level_update(j, nodes, x0, x, R_sim);
    }

    // Descent: prolong_x from coarse toward fine using current hierarchy.
    // Re-encode then prolong from root so children see updated poses.
    pyr.P.encode_from_leaves();
    for (index_t j = max_sim_L; j > 0; --j) {
      std::vector<index_t> nodes = pyr.P.extract_at_depth(j);
      for (index_t n : nodes) {
        // Snap children before prolong for R_sim on the prolong map.
        std::vector<vec3> child_x0;
        std::vector<index_t> child_ids;
        const arp::tree_view &view = pyr.P.view();
        (void)view;
        pyr.P.for_each_subtree_leaf(n, [&](index_t lid) {
          child_ids.push_back(lid);
          child_x0.push_back(pyr.P.leaf_live_sorted()[static_cast<size_t>(lid)]);
        });
        pyr.P.prolong_from(n);
        std::vector<vec3> child_x1;
        child_x1.reserve(child_x0.size());
        for (index_t lid : child_ids)
          child_x1.push_back(pyr.P.leaf_live_sorted()[static_cast<size_t>(lid)]);
        const mat3 R_pro =
            calder::kabsch_sim_rotation(child_x0, child_x1);
        for (index_t lid : child_ids) {
          pyr.leaf_T[static_cast<size_t>(lid)] =
              R_pro * pyr.leaf_T[static_cast<size_t>(lid)];
          if (pyr.leaf_T[static_cast<size_t>(lid)].squaredNorm() > real(1e-24))
            pyr.leaf_T[static_cast<size_t>(lid)].normalize();
        }
      }
      pyr.P.encode_from_leaves();
    }
  }
}

inline tp_curvature_propose
tp_propose_mg_step(asawa::rod::rod &rod, const asawa::rod::dynamic &dynamic,
                   const tangent_point_solver_config &cfg,
                   const mg_tp_config &mg) {
  using clock = std::chrono::steady_clock;
  const auto t0 = clock::now();
  tp_curvature_propose out;

  if (!mg.enable) {
    out = tp_propose_curvature_step(rod, dynamic, cfg);
    out.propose_ms =
        real(std::chrono::duration<double, std::milli>(clock::now() - t0)
                 .count());
    return out;
  }

  mg_rod_pyramid pyr = mg_rod_pyramid::build(rod);
  pyr.sync_from_rod_positions(rod);
  const std::vector<vec3> mid0 = pyr.P.leaf_positions_original();
  const index_t max_L = pyr.P.max_depth_from_leaves();
  const index_t max_sim_L = static_cast<index_t>(
      std::max(0, static_cast<int>(max_L) - mg.max_sim_offset));

  mg_vcycle(pyr, mg);

  std::vector<vec3> dvert = pyr.leaf_mid_delta_to_vert_disp(mid0);
  if (dvert.size() < rod.x().size())
    dvert.resize(rod.x().size(), vec3::Zero());

  {
    real peak = 0, rms = 0;
    size_t nnz = 0;
    for (const auto &d : dvert) {
      const real n = d.norm();
      peak = std::max(peak, n);
      rms += n * n;
      if (n > real(1e-12))
        ++nnz;
    }
    if (!dvert.empty())
      rms = std::sqrt(rms / real(dvert.size()));
    std::cout << "[mg_tp] Lmax=" << max_L << " sim_L=" << max_sim_L
              << " cycles=" << mg.num_v_cycles
              << " refine_fine=" << (mg.refine_fine ? 1 : 0)
              << " |δ|_peak=" << peak << " rms=" << rms << " nnz=" << nnz
              << std::endl;
  }

  if (mg.refine_fine) {
    std::vector<vec3> &x = rod.x();
    const std::vector<vec3> x0 = x;
    for (size_t i = 0; i < x.size(); ++i)
      x[i] += dvert[i];
    rod._update_frames();
    tp_curvature_propose fine = tp_propose_curvature_step(rod, dynamic, cfg);
    for (size_t i = 0; i < x.size(); ++i)
      x[i] = x0[i];
    rod._update_frames();
    out = fine;
    const size_t n = std::min(dvert.size(), fine.displacement.size());
    for (size_t i = 0; i < n; ++i)
      out.displacement[i] = dvert[i] + fine.displacement[i];
    out.direction = out.displacement;
    tp_peak_normalize(out.direction);
  } else {
    out.displacement = std::move(dvert);
    out.direction = out.displacement;
    tp_peak_normalize(out.direction);
    const real E0 = tp_energy(rod, dynamic, cfg);
    out.ls.E0 = E0;
    out.ls.E_tp0 = E0;
    out.ls.E_accepted = E0;
    out.ls.alpha = real(1);
  }

  out.propose_ms =
      real(std::chrono::duration<double, std::milli>(clock::now() - t0).count());
  return out;
}

} // namespace duchamp
} // namespace gaudi

#endif // GAUDI_DUCHAMP_MG_TANGENT_POINT_HPP
