#ifndef __GAUDI_TEST_EDGE_COV_TREE_TESTS_HPP__
#define __GAUDI_TEST_EDGE_COV_TREE_TESTS_HPP__

#include "gaudi/arp/edge_cov_tree.hpp"
#include "gaudi/test/test.hpp"

#include <array>
#include <cmath>
#include <set>
#include <tuple>
#include <vector>

namespace gaudi {
namespace test {

inline arp::edge_cov_tree<2>::ptr make_edge_tree(const std::vector<vec3> &verts,
                                                 const std::vector<index_t> &adj) {
  return arp::edge_cov_tree<2>::create(adj, verts);
}

inline bool box_contains(const ext::extents_t &outer, const ext::extents_t &inner) {
  for (int k = 0; k < 3; ++k) {
    if (inner[0][k] < outer[0][k] - real(1e-9))
      return false;
    if (inner[1][k] > outer[1][k] + real(1e-9))
      return false;
  }
  return true;
}

template <int N>
inline void check_edge_cov_topology(const arp::edge_cov_tree<N> &tree) {
  const index_t n = static_cast<index_t>(tree.indices_.size());
  GAUDI_ASSERT(static_cast<index_t>(tree.leaf_nodes_.size()) == n);
  GAUDI_ASSERT(static_cast<index_t>(tree.leaf_mu_c_.size()) == n);
  GAUDI_ASSERT(static_cast<index_t>(tree.leaf_mu_T_.size()) == n);
  std::set<index_t> seen(tree.indices_.begin(), tree.indices_.end());
  GAUDI_ASSERT(static_cast<index_t>(seen.size()) == n);
  for (index_t i = 0; i < n; ++i)
    GAUDI_ASSERT(seen.count(i) == 1);

  GAUDI_ASSERT(static_cast<index_t>(tree.bvh_.leaf.size()) == n);
  GAUDI_ASSERT(static_cast<index_t>(tree.leaf_mass_.size()) == n);
  if (n <= 1)
    return;

  GAUDI_ASSERT(static_cast<index_t>(tree.internal_nodes_.size()) == n - 1);
  GAUDI_ASSERT(static_cast<index_t>(tree.mu_c_.size()) == n - 1);
  GAUDI_ASSERT(static_cast<index_t>(tree.mu_T_.size()) == n - 1);
  GAUDI_ASSERT(static_cast<index_t>(tree.bvh_.internal.size()) == n - 1);
  GAUDI_ASSERT(static_cast<index_t>(tree.theta_.size()) == n - 1);
  GAUDI_ASSERT(static_cast<index_t>(tree.mass_.size()) == n - 1);
  GAUDI_ASSERT(tree.internal_nodes_[0].parent == arp::UNULL);

  for (index_t id = 0; id < n - 1; ++id) {
    const auto &node = tree.internal_nodes_[static_cast<size_t>(id)];
    GAUDI_ASSERT(node.start != arp::UNULL);
    GAUDI_ASSERT(node.split >= node.start);
    GAUDI_ASSERT(node.split < node.end);
    if (arp::left_leaf(node))
      GAUDI_ASSERT(tree.leaf_nodes_[static_cast<size_t>(node.start)].parent == id);
    else
      GAUDI_ASSERT(tree.internal_nodes_[static_cast<size_t>(node.split)].parent ==
                   id);
    if (arp::right_leaf(node))
      GAUDI_ASSERT(tree.leaf_nodes_[static_cast<size_t>(node.end)].parent == id);
    else
      GAUDI_ASSERT(
          tree.internal_nodes_[static_cast<size_t>(node.split + 1)].parent == id);

    const ext::extents_t &parent = tree.bvh_.internal[static_cast<size_t>(id)];
    const ext::extents_t &left =
        arp::left_leaf(node) ? tree.bvh_.leaf[static_cast<size_t>(node.start)]
                             : tree.bvh_.internal[static_cast<size_t>(node.split)];
    const ext::extents_t &right =
        arp::right_leaf(node) ? tree.bvh_.leaf[static_cast<size_t>(node.end)]
                              : tree.bvh_.internal[static_cast<size_t>(node.split + 1)];
    GAUDI_ASSERT(box_contains(parent, left));
    GAUDI_ASSERT(box_contains(parent, right));

    const real left_mass =
        arp::left_leaf(node) ? tree.leaf_mass(node.start) : tree.mass(node.split);
    const real right_mass = arp::right_leaf(node)
                                ? tree.leaf_mass(node.end)
                                : tree.mass(node.split + 1);
    GAUDI_ASSERT(std::abs(tree.mass(id) - (left_mass + right_mass)) < real(1e-8));
  }
}

GAUDI_TEST(edge_cov_permutation_and_parents) {
  std::vector<vec3> v;
  std::vector<index_t> adj;
  for (int i = 0; i < 12; ++i) {
    const real x = real(0.15) * i;
    const vec3 a(x, real(0.01) * (i % 3), real(0.02) * ((i / 3) % 2));
    const vec3 b = a + vec3(real(0.08), real(0.01) * (i % 2), real(0.03));
    adj.push_back(static_cast<index_t>(v.size()));
    v.push_back(a);
    adj.push_back(static_cast<index_t>(v.size()));
    v.push_back(b);
  }
  auto tree = make_edge_tree(v, adj);
  check_edge_cov_topology(*tree);
}

GAUDI_TEST(edge_cov_crossed_bundles_split) {
  std::vector<vec3> v;
  std::vector<index_t> adj;
  auto add = [&](const vec3 &a, const vec3 &b) {
    adj.push_back(static_cast<index_t>(v.size()));
    v.push_back(a);
    adj.push_back(static_cast<index_t>(v.size()));
    v.push_back(b);
  };
  for (int i = 0; i < 6; ++i) {
    const real y = real(0.05) * i;
    add(vec3(0, y, 0), vec3(real(0.2), y, 0));
  }
  for (int i = 0; i < 6; ++i) {
    const real x = real(5) + real(0.05) * i;
    add(vec3(x, 0, 0), vec3(x, real(0.2), 0));
  }
  auto tree = make_edge_tree(v, adj);
  check_edge_cov_topology(*tree);

  const auto &root = tree->internal_nodes_[0];
  int x_side[2] = {0, 0};
  int y_side[2] = {0, 0};
  for (index_t i = root.start; i <= root.end; ++i) {
    const index_t prim = tree->indices_[static_cast<size_t>(i)];
    const bool left = i <= root.split;
    if (prim < 6)
      x_side[left ? 0 : 1]++;
    else
      y_side[left ? 0 : 1]++;
  }
  const bool separated = (x_side[0] == 6 && y_side[1] == 6) ||
                         (x_side[1] == 6 && y_side[0] == 6);
  GAUDI_ASSERT(separated);
  GAUDI_ASSERT(tree->theta(0) > real(0.7));
}

GAUDI_TEST(edge_cov_mixed_sign_tangent_axis) {
  std::vector<vec3> v;
  std::vector<index_t> adj;
  vec3 mean = vec3::Zero();
  const int n = 8;
  for (int i = 0; i < n; ++i) {
    const vec3 a(real(i), 0, 0);
    const vec3 b = a + vec3(real(0.4), 0, 0);
    const bool flip = (i % 2) == 1;
    const vec3 p0 = flip ? b : a;
    const vec3 p1 = flip ? a : b;
    adj.push_back(static_cast<index_t>(v.size()));
    v.push_back(p0);
    adj.push_back(static_cast<index_t>(v.size()));
    v.push_back(p1);
    mean += real(0.5) * (a + b);
  }
  mean /= real(n);

  auto tree = make_edge_tree(v, adj);
  check_edge_cov_topology(*tree);
  const vec3 &axis = tree->mu_T(0);
  GAUDI_ASSERT(std::abs(axis.dot(vec3::UnitX())) > real(0.99));
  GAUDI_ASSERT(axis[0] > real(0));
  GAUDI_ASSERT((tree->mu_c(0) - mean).norm() < real(1e-9));
  GAUDI_ASSERT(tree->theta(0) < real(1e-4));
  GAUDI_ASSERT(tree->r_T(0) < real(1e-4));
}

GAUDI_TEST(edge_cov_long_edge_radius) {
  std::vector<vec3> v = {vec3(0, 0, 0), vec3(4, 0, 0)};
  std::vector<index_t> adj = {0, 1};
  auto tree = make_edge_tree(v, adj);
  check_edge_cov_topology(*tree);
  GAUDI_ASSERT(std::abs(tree->leaf_r_x(0) - real(2)) < real(1e-9));
  GAUDI_ASSERT(std::abs(tree->leaf_mass(0) - real(4)) < real(1e-9));
}

GAUDI_TEST(edge_cov_opening) {
  auto count_walk = [](const arp::edge_cov_tree<2> &tree, const vec3 &query,
                       real box, real cone = -1) {
    if (cone < real(0))
      cone = box;
    int far = 0;
    int leaves = 0;
    bool root_far = false;
    arp::traverse_cov_opening(
        tree, 0, query, box, cone,
        [&](index_t node_id, index_t, const vec3 &) {
          ++far;
          if (node_id == 0)
            root_far = true;
        },
        [&](index_t, index_t, index_t, const vec3 &) { ++leaves; });
    return std::tuple<int, int, bool>(far, leaves, root_far);
  };

  std::vector<vec3> parallel_v;
  std::vector<index_t> parallel_adj;
  for (int i = 0; i < 8; ++i) {
    const vec3 a(real(i), 0, 0);
    const vec3 b = a + vec3(real(0.4), 0, 0);
    const bool flip = (i % 2) == 1;
    parallel_adj.push_back(static_cast<index_t>(parallel_v.size()));
    parallel_v.push_back(flip ? b : a);
    parallel_adj.push_back(static_cast<index_t>(parallel_v.size()));
    parallel_v.push_back(flip ? a : b);
  }
  auto parallel = make_edge_tree(parallel_v, parallel_adj);
  check_edge_cov_topology(*parallel);
  const vec3 far_query = parallel->mu_c(0) + vec3(0, real(1000), 0);
  const auto [far_n, leaf_n, root_far] =
      count_walk(*parallel, far_query, real(0.05));
  GAUDI_ASSERT(root_far);
  GAUDI_ASSERT(far_n == 1);
  GAUDI_ASSERT(leaf_n == 0);
  const auto [near_far, near_leaves, near_root] =
      count_walk(*parallel, parallel->mu_c(0), real(0.05));
  (void)near_far;
  (void)near_leaves;
  GAUDI_ASSERT(!near_root);

  std::vector<vec3> crossed_v;
  std::vector<index_t> crossed_adj;
  auto add = [&](const vec3 &a, const vec3 &b) {
    crossed_adj.push_back(static_cast<index_t>(crossed_v.size()));
    crossed_v.push_back(a);
    crossed_adj.push_back(static_cast<index_t>(crossed_v.size()));
    crossed_v.push_back(b);
  };
  for (int i = 0; i < 6; ++i) {
    const real y = real(0.05) * i;
    add(vec3(0, y, 0), vec3(real(0.2), y, 0));
  }
  for (int i = 0; i < 6; ++i) {
    const real x = real(5) + real(0.05) * i;
    add(vec3(x, 0, 0), vec3(x, real(0.2), 0));
  }
  auto crossed = make_edge_tree(crossed_v, crossed_adj);
  const vec3 crossed_query = crossed->mu_c(0) + vec3(0, real(1000), 0);
  const auto [cross_far, cross_leaves, cross_root] =
      count_walk(*crossed, crossed_query, real(0.2));
  GAUDI_ASSERT(!cross_root);
  GAUDI_ASSERT(cross_far >= 1);
  GAUDI_ASSERT(cross_leaves < 12);
}

GAUDI_TEST(edge_cov_triangle_view_normal) {
  std::vector<vec3> flat;
  auto add = [&](const vec3 &a, const vec3 &b, const vec3 &c) {
    flat.push_back(a);
    flat.push_back(b);
    flat.push_back(c);
  };
  for (int i = 0; i < 4; ++i) {
    const vec3 o(real(i), 0, 0);
    const vec3 p1 = o + vec3(real(0.4), 0, 0);
    const vec3 p2 = o + vec3(0, real(0.4), 0);
    if (i % 2 == 0)
      add(o, p1, p2);
    else
      add(o, p2, p1);
  }
  for (int i = 0; i < 4; ++i) {
    const vec3 o(real(i), real(6), 0);
    add(o, o + vec3(real(0.4), 0, 0), o + vec3(0, real(0.4), 0));
  }

  simplex_view<3, std::vector<vec3>> view(flat);
  auto tree = arp::edge_cov_tree<3>::create(
      view, [](index_t, const std::array<vec3, 3> &s) {
        return (s[1] - s[0]).cross(s[2] - s[0]);
      });
  check_edge_cov_topology(*tree);

  const auto &root = tree->internal_nodes_[0];
  int side_a[2] = {0, 0};
  int side_b[2] = {0, 0};
  for (index_t i = root.start; i <= root.end; ++i) {
    const index_t prim = tree->indices_[static_cast<size_t>(i)];
    const bool left = i <= root.split;
    if (prim < 4)
      side_a[left ? 0 : 1]++;
    else
      side_b[left ? 0 : 1]++;
  }
  const bool separated = (side_a[0] == 4 && side_b[1] == 4) ||
                         (side_a[1] == 4 && side_b[0] == 4);
  GAUDI_ASSERT(separated);

  const bool a_on_left = side_a[0] == 4;
  const index_t child =
      a_on_left ? root.split : (arp::right_leaf(root) ? root.end : root.split + 1);
  GAUDI_ASSERT(!arp::left_leaf(root) || !a_on_left);
  if (a_on_left)
    GAUDI_ASSERT(!arp::left_leaf(root));
  else
    GAUDI_ASSERT(!arp::right_leaf(root));
  const vec3 &axis = tree->mu_T(child);
  GAUDI_ASSERT(std::abs(axis.dot(vec3::UnitZ())) > real(0.99));
  GAUDI_ASSERT(axis[2] > real(0));
}

GAUDI_TEST(edge_cov_point_view_captured_normals) {
  std::vector<vec3> points;
  std::vector<vec3> normals;
  for (int i = 0; i < 6; ++i) {
    points.push_back(vec3(real(0.1) * i, 0, 0));
    normals.push_back((i % 2 == 0) ? vec3(vec3::UnitX()) : vec3(-vec3::UnitX()));
  }
  for (int i = 0; i < 6; ++i) {
    points.push_back(vec3(real(5), real(0.1) * i, 0));
    normals.push_back(vec3::UnitY());
  }

  simplex_view<1, std::vector<vec3>> view(points);
  auto tree = arp::edge_cov_tree<1>::create(
      view, [&](index_t i, const std::array<vec3, 1> &) { return normals[i]; });
  check_edge_cov_topology(*tree);

  const auto &root = tree->internal_nodes_[0];
  int side_a[2] = {0, 0};
  int side_b[2] = {0, 0};
  for (index_t i = root.start; i <= root.end; ++i) {
    const index_t prim = tree->indices_[static_cast<size_t>(i)];
    const bool left = i <= root.split;
    if (prim < 6)
      side_a[left ? 0 : 1]++;
    else
      side_b[left ? 0 : 1]++;
  }
  const bool separated = (side_a[0] == 6 && side_b[1] == 6) ||
                         (side_a[1] == 6 && side_b[0] == 6);
  GAUDI_ASSERT(separated);

  const bool a_on_left = side_a[0] == 6;
  if (a_on_left)
    GAUDI_ASSERT(!arp::left_leaf(root));
  else
    GAUDI_ASSERT(!arp::right_leaf(root));
  const index_t child = a_on_left ? root.split : root.split + 1;
  const vec3 &axis = tree->mu_T(child);
  GAUDI_ASSERT(std::abs(axis.dot(vec3::UnitX())) > real(0.99));
  GAUDI_ASSERT(axis[0] > real(0));
}

// A far admissible cluster is one sample, including a node whose children are leaves.
GAUDI_TEST(edge_cov_walk_far_lumps) {
  auto count_walk = [](const arp::edge_cov_tree<2> &tree, const vec3 &query,
                       real box, real cone) {
    int far = 0;
    int leaves = 0;
    arp::traverse_cov_opening(
        tree, 0, query, box, cone,
        [&](index_t, index_t, const vec3 &) { ++far; },
        [&](index_t, index_t, index_t, const vec3 &) { ++leaves; });
    return std::pair<int, int>(far, leaves);
  };

  std::vector<vec3> two_v = {vec3(0, 0, 0), vec3(0.4, 0, 0), vec3(0.5, 0, 0),
                             vec3(0.9, 0, 0)};
  std::vector<index_t> two_adj = {0, 1, 2, 3};
  auto two = make_edge_tree(two_v, two_adj);
  const auto [two_far, two_leaves] =
      count_walk(*two, vec3(0, real(1000), 0), real(0.25), real(0.25));
  GAUDI_ASSERT(two_far == 1);
  GAUDI_ASSERT(two_leaves == 0);

  std::vector<vec3> bundle_v;
  std::vector<index_t> bundle_adj;
  for (int i = 0; i < 32; ++i) {
    const vec3 a(real(i) * real(0.4), 0, 0);
    bundle_adj.push_back(static_cast<index_t>(bundle_v.size()));
    bundle_v.push_back(a);
    bundle_adj.push_back(static_cast<index_t>(bundle_v.size()));
    bundle_v.push_back(a + vec3(real(0.4), 0, 0));
  }
  auto bundle = make_edge_tree(bundle_v, bundle_adj);
  const vec3 far_q = bundle->mu_c(0) + vec3(0, real(1000), 0);
  const auto [far_n, far_leaves] =
      count_walk(*bundle, far_q, real(0.25), real(0.25));
  GAUDI_ASSERT(far_n == 1);
  GAUDI_ASSERT(far_leaves == 0);
  const auto [closed_n, closed_leaves] =
      count_walk(*bundle, far_q, real(0), real(0));
  GAUDI_ASSERT(closed_n == 0);
  GAUDI_ASSERT(closed_leaves == 32);

  std::vector<vec3> split_v;
  std::vector<index_t> split_adj;
  auto add = [&](const vec3 &a, const vec3 &b) {
    split_adj.push_back(static_cast<index_t>(split_v.size()));
    split_v.push_back(a);
    split_adj.push_back(static_cast<index_t>(split_v.size()));
    split_v.push_back(b);
  };
  for (int i = 0; i < 6; ++i) {
    const real y = real(0.05) * i;
    add(vec3(0, y, 0), vec3(real(0.2), y, 0));
  }
  for (int i = 0; i < 6; ++i) {
    const real x = real(5) + real(0.05) * i;
    add(vec3(x, 0, 0), vec3(x, real(0.2), 0));
  }
  auto split = make_edge_tree(split_v, split_adj);
  const vec3 on_first = vec3(real(0.1), real(0.125), 0);
  const auto [near_far, near_leaves] =
      count_walk(*split, on_first, real(0.25), real(1));
  GAUDI_ASSERT(near_far >= 1);
  GAUDI_ASSERT(near_leaves > 0);
  GAUDI_ASSERT(near_far + near_leaves < 12);
}

} // namespace test
} // namespace gaudi

#endif
