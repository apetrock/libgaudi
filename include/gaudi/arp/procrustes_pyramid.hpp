#ifndef GAUDI_ARP_PROCRUSTES_PYRAMID_HPP
#define GAUDI_ARP_PROCRUSTES_PYRAMID_HPP

#include "gaudi/arp/tree_view.hpp"
#include "gaudi/common.h"

#include <Eigen/SVD>
#include <algorithm>
#include <vector>

namespace gaudi {
namespace arp {

/// Kabsch SO(3) from cross-covariance M = Σ (x−c)(x0−c0)ᵀ. Det +1.
inline mat3 procrustes_rotation(const mat3 &M) {
  if (M.squaredNorm() < real(1e-24))
    return mat3::Identity();
  Eigen::JacobiSVD<mat3> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
  mat3 U = svd.matrixU();
  mat3 V = svd.matrixV();
  mat3 R = U * V.transpose();
  if (R.determinant() < real(0.0)) {
    U.col(2) *= real(-1.0);
    R = U * V.transpose();
  }
  return R;
}

/// Uncentered position moments for build_pyramid (mass, Σx, Σxxᵀ).
struct pos_moments {
  real m = real(0);
  vec3 sum = vec3::Zero();
  mat3 sum_xx = mat3::Zero();

  static pos_moments from_point(const vec3 &x) {
    pos_moments q;
    q.m = real(1);
    q.sum = x;
    q.sum_xx = x * x.transpose();
    return q;
  }

  pos_moments &operator+=(const pos_moments &o) {
    m += o.m;
    sum += o.sum;
    sum_xx += o.sum_xx;
    return *this;
  }

  friend pos_moments operator+(pos_moments a, const pos_moments &b) {
    a += b;
    return a;
  }

  vec3 mean() const {
    if (!(m > real(0)))
      return vec3::Zero();
    return sum / m;
  }

  /// Centered covariance (Steiner): Σ(x−c)(x−c)ᵀ / m.
  mat3 centered_cov() const {
    if (!(m > real(0)))
      return mat3::Zero();
    const vec3 c = mean();
    return sum_xx / m - c * c.transpose();
  }
};

inline pos_moments moments_identity() { return pos_moments{}; }

/// Procrustes residual pyramid on an ARP tree_view.
/// Positions only: δ = Rᵀ(x − c), x = c + R δ. Tangents deferred (TP phase).
class procrustes_pyramid {
public:
  void build(const tree_view &view, const std::vector<vec3> &leaf_x_original) {
    view_ = view;
    const size_t nL = view_.leaf_count();
    const size_t nN = view_.internal_count();

    leaf_rest_.assign(nL, vec3::Zero());
    leaf_live_.assign(nL, vec3::Zero());
    leaf_delta_.assign(nL, vec3::Zero());
    for (size_t i = 0; i < nL; ++i) {
      const index_t orig = view_.index(static_cast<index_t>(i));
      leaf_rest_[i] = leaf_x_original[static_cast<size_t>(orig)];
      leaf_live_[i] = leaf_rest_[i];
    }

    avg_rest_.assign(nN, vec3::Zero());
    avg_live_.assign(nN, vec3::Zero());
    cov_R_.assign(nN, mat3::Zero());
    node_delta_.assign(nN, vec3::Zero());
    R_live_.assign(nN, mat3::Identity());
    R_rel_.assign(nN, mat3::Identity());

    if (nL == 0)
      return;

    rebuild_avgs_from_leaves(/*use_live=*/false);
    avg_live_ = avg_rest_;
    for (size_t i = 0; i < nN; ++i) {
      R_live_[i] = mat3::Identity();
      R_rel_[i] = mat3::Identity();
    }

    encode_deltas();
  }

  index_t root() const {
    for (size_t i = 0; i < view_.internal_count(); ++i) {
      if (view_.internal_parent(static_cast<index_t>(i)) == UNULL)
        return static_cast<index_t>(i);
    }
    return UNULL;
  }

  const tree_view &view() const { return view_; }
  size_t leaf_count() const { return leaf_live_.size(); }
  size_t node_count() const { return avg_live_.size(); }

  const std::vector<vec3> &avg_rest() const { return avg_rest_; }
  const std::vector<vec3> &avg_live() const { return avg_live_; }
  const std::vector<mat3> &cov_R() const { return cov_R_; }
  const std::vector<vec3> &leaf_delta() const { return leaf_delta_; }
  const std::vector<vec3> &node_delta() const { return node_delta_; }
  const std::vector<mat3> &R_live() const { return R_live_; }

  /// Leaf positions in original (pre-sort) order.
  std::vector<vec3> leaf_positions_original() const {
    std::vector<vec3> out(leaf_live_.size(), vec3::Zero());
    for (size_t i = 0; i < leaf_live_.size(); ++i) {
      const index_t orig = view_.index(static_cast<index_t>(i));
      out[static_cast<size_t>(orig)] = leaf_live_[i];
    }
    return out;
  }

  const std::vector<vec3> &leaf_rest_sorted() const { return leaf_rest_; }
  const std::vector<vec3> &leaf_live_sorted() const { return leaf_live_; }

  /// Write leaf positions in original (pre-sort) order into live state.
  void set_leaf_positions_original(const std::vector<vec3> &x_original) {
    for (size_t i = 0; i < leaf_live_.size(); ++i) {
      const index_t orig = view_.index(static_cast<index_t>(i));
      leaf_live_[i] = x_original[static_cast<size_t>(orig)];
    }
  }

  void set_node_pose(index_t node, const vec3 &c, const mat3 &R) {
    if (node < 0 || static_cast<size_t>(node) >= avg_live_.size())
      return;
    avg_live_[static_cast<size_t>(node)] = c;
    R_live_[static_cast<size_t>(node)] = R;
  }

  /// Fine → coarse: moments from live leaves, Kabsch \(R\) vs rest, encode δ/R_rel.
  void restrict_from_leaves() {
    rebuild_avgs_from_leaves(/*use_live=*/true);
    compute_kabsch_R_live();
    encode_deltas();
    encode_R_rel();
  }

  /// Subtree leaf counts (internal nodes), via ones-pyramid.
  std::vector<real> subtree_leaf_counts() const {
    std::vector<real> ones(leaf_live_.size(), real(1));
    return build_pyramid<real, real>(
        view_, ones,
        [](const real &a, const real &) { return a; },
        [](const real &a, const real &b) { return a + b; }, real(0));
  }

  /// Partition cut: emit a node when its subtree leaf count ≤ \p max_leaves
  /// (larger \p max_leaves → coarser cut). Root always covered.
  std::vector<index_t> extract_cut(real max_leaves) const {
    std::vector<index_t> out;
    if (view_.internal_count() == 0)
      return out;
    const std::vector<real> counts = subtree_leaf_counts();
    const real thresh = std::max(max_leaves, real(1));
    const index_t r = root();
    if (r == UNULL)
      return out;
    std::vector<index_t> stack;
    stack.push_back(r);
    while (!stack.empty()) {
      const index_t n = stack.back();
      stack.pop_back();
      const real cnt = counts[static_cast<size_t>(n)];
      if (cnt <= thresh) {
        out.push_back(n);
        continue;
      }
      bool pushed = false;
      visit_children(n, [&](index_t /*p*/, bool is_leaf, index_t id) {
        if (is_leaf) {
          // Leaf under an oversize node: still emit the parent once.
          // (Handled by falling through if no internal children.)
        } else {
          stack.push_back(id);
          pushed = true;
        }
      });
      if (!pushed)
        out.push_back(n);
    }
    return out;
  }

  /// Height of each internal node = 1 + max child height (leaf height = 0).
  std::vector<index_t> depth_from_leaves() const {
    std::vector<index_t> d(view_.internal_count(), 0);
    // Post-order: process nodes with larger ids first is not guaranteed.
    // Iterate until stable (tree height ≤ 64).
    for (int pass = 0; pass < 64; ++pass) {
      bool changed = false;
      for (size_t n = 0; n < view_.internal_count(); ++n) {
        index_t best = 0;
        visit_children(static_cast<index_t>(n),
                       [&](index_t /*p*/, bool is_leaf, index_t id) {
                         const index_t h =
                             is_leaf ? 0 : d[static_cast<size_t>(id)];
                         best = std::max(best, static_cast<index_t>(h + 1));
                       });
        if (best != d[n]) {
          d[n] = best;
          changed = true;
        }
      }
      if (!changed)
        break;
    }
    return d;
  }

  index_t max_depth_from_leaves() const {
    const auto d = depth_from_leaves();
    index_t m = 0;
    for (index_t v : d)
      m = std::max(m, v);
    return m;
  }

  /// Nodes whose depth-from-leaves equals \p lvl (lvl≥1). Empty if none.
  std::vector<index_t> extract_at_depth(index_t lvl) const {
    std::vector<index_t> out;
    if (lvl < 1)
      return out;
    const auto d = depth_from_leaves();
    for (size_t n = 0; n < d.size(); ++n) {
      if (d[n] == lvl)
        out.push_back(static_cast<index_t>(n));
    }
    return out;
  }

  template <typename Fn>
  void for_each_subtree_leaf(index_t node, Fn &&fn) const {
    if (node < 0 || static_cast<size_t>(node) >= view_.internal_count())
      return;
    std::vector<index_t> stack;
    stack.push_back(node);
    while (!stack.empty()) {
      const index_t n = stack.back();
      stack.pop_back();
      visit_children(n, [&](index_t /*p*/, bool is_leaf, index_t id) {
        if (is_leaf)
          fn(id);
        else
          stack.push_back(id);
      });
    }
  }

  /// Top-down: x = c + R δ for descendants of \p node.
  void prolong_from(index_t node) {
    if (node < 0 || static_cast<size_t>(node) >= avg_live_.size())
      return;
    std::vector<index_t> stack;
    stack.push_back(node);
    while (!stack.empty()) {
      const index_t n = stack.back();
      stack.pop_back();
      const vec3 &c = avg_live_[static_cast<size_t>(n)];
      const mat3 &R = R_live_[static_cast<size_t>(n)];
      visit_children(n, [&](index_t /*parent*/, bool is_leaf, index_t id) {
        if (is_leaf) {
          leaf_live_[static_cast<size_t>(id)] =
              c + R * leaf_delta_[static_cast<size_t>(id)];
        } else {
          avg_live_[static_cast<size_t>(id)] =
              c + R * node_delta_[static_cast<size_t>(id)];
          R_live_[static_cast<size_t>(id)] =
              R * R_rel_[static_cast<size_t>(id)];
          stack.push_back(id);
        }
      });
    }
  }

  /// Recompute δ and R_rel from current leaf_live / R_live (avgs refreshed).
  void encode_from_leaves() {
    rebuild_avgs_from_leaves(/*use_live=*/true);
    encode_deltas();
    encode_R_rel();
  }

  void prolong_identity_from_root() {
    const index_t r = root();
    if (r == UNULL)
      return;
    prolong_from(r);
  }

private:
  tree_view view_;
  std::vector<vec3> leaf_rest_;
  std::vector<vec3> leaf_live_;
  std::vector<vec3> leaf_delta_;
  std::vector<vec3> avg_rest_;
  std::vector<vec3> avg_live_;
  std::vector<mat3> cov_R_;
  std::vector<vec3> node_delta_;
  std::vector<mat3> R_live_;
  /// Rotation of node frame relative to parent (root: Identity; absolute in R_live).
  std::vector<mat3> R_rel_;

  void rebuild_avgs_from_leaves(bool use_live) {
    const std::vector<vec3> &src = use_live ? leaf_live_ : leaf_rest_;
    if (view_.internal_count() == 0)
      return;

    auto moments = build_pyramid<vec3, pos_moments>(
        view_, src,
        [](const vec3 &x, const pos_moments &) {
          return pos_moments::from_point(x);
        },
        [](const pos_moments &a, const pos_moments &b) { return a + b; },
        moments_identity());

    auto &avg = use_live ? avg_live_ : avg_rest_;
    avg.resize(moments.size());
    cov_R_.resize(moments.size());
    for (size_t i = 0; i < moments.size(); ++i) {
      avg[i] = moments[i].mean();
      cov_R_[i] = moments[i].centered_cov();
    }
  }

  void compute_kabsch_R_live() {
    const size_t nN = avg_live_.size();
    std::vector<mat3> sum_lr(nN, mat3::Zero());
    std::vector<real> mass(nN, real(0));
    for (size_t i = 0; i < leaf_live_.size(); ++i) {
      const vec3 &xl = leaf_live_[i];
      const vec3 &xr = leaf_rest_[i];
      index_t parent = view_.leaf_parent(static_cast<index_t>(i));
      index_t depth = 0;
      while (depth < 64 && parent != UNULL) {
        sum_lr[static_cast<size_t>(parent)] += xl * xr.transpose();
        mass[static_cast<size_t>(parent)] += real(1);
        parent = view_.internal_parent(parent);
        ++depth;
      }
    }
    for (size_t n = 0; n < nN; ++n) {
      if (!(mass[n] > real(0))) {
        R_live_[n] = mat3::Identity();
        continue;
      }
      // Σ(xl−cl)(xr−cr)ᵀ = Σ xl xrᵀ − m cl crᵀ
      const mat3 M =
          sum_lr[n] - mass[n] * (avg_live_[n] * avg_rest_[n].transpose());
      R_live_[n] = procrustes_rotation(M);
    }
  }

  void encode_R_rel() {
    for (size_t n = 0; n < R_live_.size(); ++n) {
      const index_t p = view_.internal_parent(static_cast<index_t>(n));
      if (p == UNULL) {
        R_rel_[n] = mat3::Identity();
        continue;
      }
      R_rel_[n] = R_live_[static_cast<size_t>(p)].transpose() * R_live_[n];
    }
  }

  void encode_deltas() {
    for (size_t i = 0; i < leaf_live_.size(); ++i) {
      const index_t p = view_.leaf_parent(static_cast<index_t>(i));
      if (p == UNULL) {
        leaf_delta_[i] = vec3::Zero();
        continue;
      }
      const mat3 &Rp = R_live_[static_cast<size_t>(p)];
      const vec3 &cp = avg_live_[static_cast<size_t>(p)];
      leaf_delta_[i] = Rp.transpose() * (leaf_live_[i] - cp);
    }
    for (size_t n = 0; n < avg_live_.size(); ++n) {
      const index_t p = view_.internal_parent(static_cast<index_t>(n));
      if (p == UNULL) {
        node_delta_[n] = vec3::Zero();
        continue;
      }
      const mat3 &Rp = R_live_[static_cast<size_t>(p)];
      const vec3 &cp = avg_live_[static_cast<size_t>(p)];
      node_delta_[n] = Rp.transpose() * (avg_live_[n] - cp);
    }
  }

  template <typename Fn>
  void visit_children(index_t node, Fn &&fn) const {
    const auto &cnode = view_.internal_nodes()[static_cast<size_t>(node)];
    if (cnode.split == UNULL)
      return;
    if (both_leaves(cnode)) {
      fn(node, true, cnode.start);
      fn(node, true, cnode.end);
    } else if (left_leaf(cnode)) {
      fn(node, true, cnode.start);
      fn(node, false, cnode.split + 1);
    } else if (right_leaf(cnode)) {
      fn(node, true, cnode.end);
      fn(node, false, cnode.split);
    } else {
      fn(node, false, cnode.split);
      fn(node, false, cnode.split + 1);
    }
  }
};

} // namespace arp
} // namespace gaudi

#endif
