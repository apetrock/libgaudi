#ifndef __GAUDI_ARP_EDGE_COV_TREE__
#define __GAUDI_ARP_EDGE_COV_TREE__

#include "gaudi/arp/hash_tree.hpp"
#include "gaudi/common.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <tuple>
#include <utility>
#include <vector>

namespace gaudi {
namespace arp {

// Sign convention for an axis defined only up to ±. The largest-magnitude
// component is made positive so two runs do not swap children.
inline void orient_largest_positive(vec3 &v) {
  int k = 0;
  for (int i = 1; i < 3; ++i)
    if (std::abs(v[i]) > std::abs(v[k]))
      k = i;
  if (v[k] < 0)
    v = -v;
}

inline void orient_largest_positive(vec6 &v) {
  int k = 0;
  for (int i = 1; i < 6; ++i)
    if (std::abs(v[i]) > std::abs(v[k]))
      k = i;
  if (v[k] < 0)
    v = -v;
}

inline vec3 leading_axis(const mat3 &M) {
  Eigen::SelfAdjointEigenSolver<mat3> es(M);
  vec3 axis = vec3::UnitX();
  if (es.info() == Eigen::Success)
    axis = es.eigenvectors().col(2);
  if (axis.squaredNorm() < real(1e-24))
    axis = vec3::UnitX();
  axis.normalize();
  orient_largest_positive(axis);
  return axis;
}

template <class Simplex> inline vec3 tuple_mean(const Simplex &s) {
  using S = std::decay_t<Simplex>;
  constexpr int kN = static_cast<int>(std::tuple_size<S>::value);
  vec3 c = vec3::Zero();
  for (int k = 0; k < kN; ++k)
    c += s[k];
  return c / real(kN);
}

// Direction is whatever the caller returns for this simplex. Unit length so
// the tangent tensor lives on the sphere.
template <class DirectionFn, class Simplex>
inline vec3 unit_direction(DirectionFn &direction, index_t i, const Simplex &s) {
  vec3 T = direction(i, s);
  const real len = T.norm();
  return len > real(1e-12) ? vec3(T / len) : vec3::UnitX();
}

template <class TupleView, class DirectionFn>
inline std::pair<vec3, vec3> sample_center_tangent(const TupleView &tuples,
                                                   DirectionFn &direction,
                                                   index_t i) {
  const auto s = tuples[static_cast<size_t>(i)];
  return {tuple_mean(s), unit_direction(direction, i, s)};
}

// One outer product of the 6D sample. build_cov divides by the count.
struct cov_acc {
  mat6 C = mat6::Zero();
  index_t count = 0;
};

inline void append_cov(cov_acc &acc, const vec6 &u) {
  acc.C += u * u.transpose();
  ++acc.count;
}

inline mat6 build_cov(const cov_acc &acc) {
  if (acc.count == 0)
    return mat6::Zero();
  return acc.C / real(acc.count);
}

struct cov_axis {
  vec6 axis = vec6::Zero();
  bool degenerate = true;
};

inline cov_axis solve(const mat6 &cov) {
  Eigen::SelfAdjointEigenSolver<mat6> es(cov);
  cov_axis out;
  if (es.info() != Eigen::Success)
    return out;
  vec6 axis = es.eigenvectors().col(5);
  if (es.eigenvalues()(5) < real(1e-18) || axis.squaredNorm() < real(1e-24))
    return out;
  axis.normalize();
  orient_largest_positive(axis);
  out.axis = axis;
  out.degenerate = false;
  return out;
}

// Angle from axis to T after the split's hemisphere flip. Result is in [0, π/2].
inline real tangent_cone_angle(const vec3 &axis, vec3 T) {
  if (T.dot(axis) < real(0))
    T = -T;
  const real c = std::min(real(1), std::max(real(-1), T.dot(axis)));
  return std::acos(c);
}

// Euclidean radius of a unit tangent at angle theta from the axis.
inline real cone_radius(real theta) {
  return real(2) * std::sin(real(0.5) * theta);
}

// Leaf measure the far kernel scales by: edge length, triangle area, or 1.
template <size_t kN>
inline real simplex_measure(const std::array<vec3, kN> &s) {
  if constexpr (kN == 1)
    return real(1);
  else if constexpr (kN == 2)
    return (s[1] - s[0]).norm();
  else {
    static_assert(kN == 3, "simplex measure is a point, an edge, or a triangle");
    return real(0.5) * (s[1] - s[0]).cross(s[2] - s[0]).norm();
  }
}

// Half-angle of each internal node about its own μ_T. A merged cone would
// grow a different axis. Each leaf is measured against every ancestor.
inline std::vector<real>
make_cones(const std::vector<vec3> &leaf_T, const std::vector<vec3> &node_T,
           const std::vector<radix_tree_node> &internal_nodes,
           const std::vector<radix_tree_node> &leaf_nodes) {
  const index_t n = static_cast<index_t>(leaf_T.size());
  if (n <= 1)
    return {};

  std::vector<real> theta(static_cast<size_t>(n - 1), real(0));
  for (index_t i = 0; i < n; ++i) {
    const vec3 &T = leaf_T[static_cast<size_t>(i)];
    index_t parent = leaf_nodes[static_cast<size_t>(i)].parent;
    index_t depth = 0;
    while (depth < 64 && parent != UNULL) {
      const real ang =
          tangent_cone_angle(node_T[static_cast<size_t>(parent)], T);
      real &slot = theta[static_cast<size_t>(parent)];
      if (ang > slot)
        slot = ang;
      parent = internal_nodes[static_cast<size_t>(parent)].parent;
      ++depth;
    }
  }
  return theta;
}

// Breadth-first half-space build. One plane per bin: the leading eigenvector
// of the centered, scale-balanced covariance of (dx, dT). Node ids follow the
// same range labeling as build_aabb_radix (internal end is inclusive).
template <class TupleView, class DirectionFn>
inline std::tuple<std::vector<index_t>, std::vector<radix_tree_node>,
                  std::vector<radix_tree_node>, std::vector<vec3>,
                  std::vector<vec3>, std::vector<vec3>, std::vector<vec3>>
build_edge_cov(const TupleView &tuples, DirectionFn &&direction_fn) {
  DirectionFn &direction = direction_fn;
  const index_t n = static_cast<index_t>(tuples.size());
  std::vector<index_t> order(static_cast<size_t>(n));
  std::iota(order.begin(), order.end(), 0);

  std::vector<radix_tree_node> leaf(static_cast<size_t>(n));
  for (index_t i = 0; i < n; ++i) {
    leaf[static_cast<size_t>(i)].start = i;
    leaf[static_cast<size_t>(i)].end = i + 1;
    leaf[static_cast<size_t>(i)].split = UNULL;
    leaf[static_cast<size_t>(i)].parent = UNULL;
  }

  auto fill_leaves = [&](std::vector<vec3> &mu_c, std::vector<vec3> &mu_T) {
    mu_c.assign(static_cast<size_t>(n), vec3::Zero());
    mu_T.assign(static_cast<size_t>(n), vec3::UnitX());
    for (index_t i = 0; i < n; ++i) {
      auto [c, T] = sample_center_tangent(tuples, direction, order[static_cast<size_t>(i)]);
      mu_c[static_cast<size_t>(i)] = c;
      orient_largest_positive(T);
      mu_T[static_cast<size_t>(i)] = T;
    }
  };

  if (n <= 1) {
    std::vector<vec3> mu_c_leaf;
    std::vector<vec3> mu_T_leaf;
    fill_leaves(mu_c_leaf, mu_T_leaf);
    return {std::move(order), {}, std::move(leaf), {}, {}, std::move(mu_c_leaf),
            std::move(mu_T_leaf)};
  }

  std::vector<radix_tree_node> internal(static_cast<size_t>(n - 1));
  std::vector<vec3> mu_c(static_cast<size_t>(n - 1), vec3::Zero());
  std::vector<vec3> mu_T(static_cast<size_t>(n - 1), vec3::UnitX());

  struct Frame {
    index_t s, e, id, parent;
  };
  std::vector<Frame> level;
  std::vector<Frame> next;
  level.push_back({0, n - 1, 0, UNULL});

  constexpr real kREps = real(1e-12);

  while (!level.empty()) {
    next.clear();
    for (const Frame &f : level) {
      const index_t count = f.e - f.s + 1;
      const real inv_count = real(1) / real(count);

      mat3 M = mat3::Zero();
      vec3 mean_c = vec3::Zero();
      for (index_t i = f.s; i <= f.e; ++i) {
        auto [c, T] = sample_center_tangent(tuples, direction, order[static_cast<size_t>(i)]);
        M += T * T.transpose();
        mean_c += c;
      }
      M *= inv_count;
      mean_c *= inv_count;
      const vec3 axis = leading_axis(M);
      mu_c[static_cast<size_t>(f.id)] = mean_c;
      mu_T[static_cast<size_t>(f.id)] = axis;

      real r2 = 0;
      for (index_t i = f.s; i <= f.e; ++i) {
        auto [c, T] = sample_center_tangent(tuples, direction, order[static_cast<size_t>(i)]);
        (void)T;
        r2 += (c - mean_c).squaredNorm();
      }
      const real R = std::sqrt(r2 * inv_count);
      const real invR = R > kREps ? real(1) / R : real(1);

      auto state_of = [&](index_t prim) {
        auto [c, T] = sample_center_tangent(tuples, direction, prim);
        if (T.dot(axis) < real(0))
          T = -T;
        vec6 u;
        u << (c - mean_c) * invR, (T - axis);
        return u;
      };

      cov_acc acc;
      std::vector<vec6> states(static_cast<size_t>(count));
      for (index_t k = 0; k < count; ++k) {
        states[static_cast<size_t>(k)] =
            state_of(order[static_cast<size_t>(f.s + k)]);
        append_cov(acc, states[static_cast<size_t>(k)]);
      }
      const mat6 cov = build_cov(acc);
      const cov_axis solved = solve(cov);
      const bool degenerate = solved.degenerate;
      const vec6 E = solved.axis;

      index_t n_left = 0;
      if (!degenerate) {
        for (index_t k = 0; k < count; ++k)
          if (states[static_cast<size_t>(k)].dot(E) <= real(0))
            ++n_left;
      }
      const bool empty_side = degenerate || n_left == 0 || n_left == count;

      index_t m = UNULL;
      if (!empty_side) {
        std::stable_partition(order.begin() + f.s, order.begin() + f.e + 1,
                              [&](index_t prim) {
                                return state_of(prim).dot(E) <= real(0);
                              });
        m = f.s + n_left - 1;
      } else {
        int slot = 0;
        if (!degenerate) {
          real best = cov(0, 0);
          for (int k = 1; k < 6; ++k)
            if (cov(k, k) > best) {
              best = cov(k, k);
              slot = k;
            }
        }
        std::vector<vec6> by_prim(static_cast<size_t>(n));
        for (index_t k = 0; k < count; ++k)
          by_prim[static_cast<size_t>(order[static_cast<size_t>(f.s + k)])] =
              states[static_cast<size_t>(k)];
        const index_t mid = f.s + (f.e - f.s) / 2;
        std::nth_element(
            order.begin() + f.s, order.begin() + mid, order.begin() + f.e + 1,
            [&](index_t a, index_t b) {
              return by_prim[static_cast<size_t>(a)][slot] <
                     by_prim[static_cast<size_t>(b)][slot];
            });
        m = mid;
      }

      internal[static_cast<size_t>(f.id)].start = f.s;
      internal[static_cast<size_t>(f.id)].end = f.e;
      internal[static_cast<size_t>(f.id)].split = m;
      internal[static_cast<size_t>(f.id)].parent = f.parent;

      if (m == f.s)
        leaf[static_cast<size_t>(f.s)].parent = f.id;
      else
        next.push_back({f.s, m, m, f.id});

      if (m + 1 == f.e)
        leaf[static_cast<size_t>(f.e)].parent = f.id;
      else
        next.push_back({m + 1, f.e, m + 1, f.id});
    }
    level.swap(next);
  }

  std::vector<vec3> mu_c_leaf;
  std::vector<vec3> mu_T_leaf;
  fill_leaves(mu_c_leaf, mu_T_leaf);
  return {std::move(order), std::move(internal), std::move(leaf),
          std::move(mu_c), std::move(mu_T), std::move(mu_c_leaf),
          std::move(mu_T_leaf)};
}

// Binary tree over simplices. Each node stores the mean center, the tangent
// axis, a position box, and the cone about that axis. The walk that opens on
// r_x and r_T is separate.
template <int N> class edge_cov_tree {
public:
  static constexpr int kSimplexN = N;

  std::vector<index_t> indices_;
  std::vector<radix_tree_node> internal_nodes_;
  std::vector<radix_tree_node> leaf_nodes_;
  std::vector<vec3> mu_c_;
  std::vector<vec3> mu_T_;
  std::vector<vec3> leaf_mu_c_;
  std::vector<vec3> leaf_mu_T_;
  TreeResult<ext::extents_t> bvh_;
  std::vector<real> theta_;
  std::vector<real> mass_;
  std::vector<real> leaf_mass_;
  std::vector<vec3> data_;
  std::vector<index_t> adjacency_;

  using ptr = std::shared_ptr<edge_cov_tree<N>>;

  static int build_count() { return builds_; }

  edge_cov_tree() { ++builds_; }

  static vec3 geometric_direction(index_t, const std::array<vec3, N> &s) {
    if constexpr (N == 1)
      return vec3::UnitX();
    else if constexpr (N == 2)
      return s[1] - s[0];
    else {
      static_assert(N == 3, "edge_cov_tree stride is 1, 2, or 3");
      return (s[1] - s[0]).cross(s[2] - s[0]);
    }
  }

  template <class TupleView, class DirectionFn>
  static ptr create(const TupleView &tuples, DirectionFn &&direction) {
    auto tree = std::make_shared<edge_cov_tree<N>>();
    tree->update(tuples, std::forward<DirectionFn>(direction));
    return tree;
  }

  static ptr create(const std::vector<index_t> &adjacency,
                    const std::vector<vec3> &vertices) {
    return std::make_shared<edge_cov_tree<N>>(vertices, adjacency);
  }

  edge_cov_tree(const std::vector<vec3> &vertices,
                const std::vector<index_t> &adjacency) {
    ++builds_;
    update(vertices, adjacency);
  }

  index_t get_index(size_t i) const { return indices_[i]; }

  const vec3 &mu_c(index_t internal_id) const {
    return mu_c_[static_cast<size_t>(internal_id)];
  }
  const vec3 &mu_T(index_t internal_id) const {
    return mu_T_[static_cast<size_t>(internal_id)];
  }

  // Distance from μ_c to the farthest corner of the position box.
  real r_x(index_t internal_id) const {
    return ext::dist_from_farthest(bvh_.internal[static_cast<size_t>(internal_id)],
                                   mu_c_[static_cast<size_t>(internal_id)]);
  }
  real leaf_r_x(index_t leaf_id) const {
    return ext::dist_from_farthest(bvh_.leaf[static_cast<size_t>(leaf_id)],
                                   leaf_mu_c_[static_cast<size_t>(leaf_id)]);
  }

  real theta(index_t internal_id) const {
    return theta_[static_cast<size_t>(internal_id)];
  }
  // Chord length of the cone: ||T - μ_T|| for a unit tangent at the rim.
  real r_T(index_t internal_id) const { return cone_radius(theta(internal_id)); }

  real mass(index_t internal_id) const {
    return mass_[static_cast<size_t>(internal_id)];
  }
  real leaf_mass(index_t leaf_id) const {
    return leaf_mass_[static_cast<size_t>(leaf_id)];
  }

  template <class TupleView, class DirectionFn>
  void update(const TupleView &tuples, DirectionFn &&direction) {
    adopt(build_edge_cov(tuples, std::forward<DirectionFn>(direction)));
    if (indices_.empty()) {
      bvh_ = {};
      theta_.clear();
      mass_.clear();
      leaf_mass_.clear();
      return;
    }
    auto permuted = permuted_simplex(tuples, indices_);
    bvh_ = make_bvh(permuted, internal_nodes_, leaf_nodes_);
    theta_ = make_cones(leaf_mu_T_, mu_T_, internal_nodes_, leaf_nodes_);
    leaf_mass_.resize(permuted.size());
    for (size_t i = 0; i < permuted.size(); ++i)
      leaf_mass_[i] = simplex_measure(permuted[i]);
    mass_ = build_pyramid(
        leaf_mass_, internal_nodes_, leaf_nodes_,
        [](real a, real b) { return a + b; }, real(0));
  }

  void update(const std::vector<vec3> &vertices,
              const std::vector<index_t> &adjacency) {
    data_ = vertices;
    adjacency_ = adjacency;
    const index_t n = static_cast<index_t>(adjacency_.size() / N);
    std::vector<index_t> identity(static_cast<size_t>(n));
    std::iota(identity.begin(), identity.end(), 0);
    permuted_simplex_view<N, std::vector<vec3>, std::vector<index_t>> view(
        data_, adjacency_, identity);
    update(view, geometric_direction);
  }

private:
  static int builds_;

  template <class Built> void adopt(Built &&built) {
    indices_ = std::move(std::get<0>(built));
    internal_nodes_ = std::move(std::get<1>(built));
    leaf_nodes_ = std::move(std::get<2>(built));
    mu_c_ = std::move(std::get<3>(built));
    mu_T_ = std::move(std::get<4>(built));
    leaf_mu_c_ = std::move(std::get<5>(built));
    leaf_mu_T_ = std::move(std::get<6>(built));
  }
};

template <int N> int edge_cov_tree<N>::builds_ = 0;

// Point query against the curve. The opening test runs on every internal node,
// including one whose child is a leaf. A passing cluster is one far sample.
// Otherwise a leaf child is evaluated and an internal child is pushed.
// box_eps ≤ 0 or cone_eps ≤ 0 never opens.
template <int N, class OnFar, class OnLeaf>
void traverse_cov_opening(const edge_cov_tree<N> &tree, index_t query_id,
                          const vec3 &query, real box_eps, real cone_eps,
                          OnFar &&on_far, OnLeaf &&on_leaf) {
  const auto &internal_nodes = tree.internal_nodes_;
  if (internal_nodes.empty())
    return;

  std::vector<index_t> stack;
  stack.reserve(128);
  stack.push_back(0);

  while (!stack.empty()) {
    const index_t node_id = stack.back();
    stack.pop_back();
    const auto &node = internal_nodes[static_cast<size_t>(node_id)];
    if (node.split == UNULL)
      continue;

    const real dist = (query - tree.mu_c(node_id)).norm();
    if (box_eps > real(0) && cone_eps > real(0) &&
        tree.r_x(node_id) < dist * box_eps &&
        tree.r_T(node_id) < cone_eps) {
      on_far(node_id, query_id, query);
      continue;
    }

    if (left_leaf(node)) {
      const index_t leaf_id = node.start;
      on_leaf(leaf_id, tree.get_index(static_cast<size_t>(leaf_id)), query_id,
              query);
    } else {
      stack.push_back(node.split);
    }

    if (right_leaf(node)) {
      const index_t leaf_id = node.end;
      on_leaf(leaf_id, tree.get_index(static_cast<size_t>(leaf_id)), query_id,
              query);
    } else {
      stack.push_back(node.split + 1);
    }
  }
}

} // namespace arp
} // namespace gaudi

#endif
