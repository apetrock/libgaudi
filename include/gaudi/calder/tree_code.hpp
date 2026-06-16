#ifndef __M2TREE_CODE__
#define __M2TREE_CODE__

#include "gaudi/arp/datums.hpp"
#include "gaudi/common.h"
#include "gaudi/geometry_types.hpp"
#include "gaudi/arp/hash_tree.hpp"
#include "gaudi/arp/tree_view.hpp"
#include "gaudi/logger.hpp"
#include <vector>
#include "gaudi/geometry_logger.hpp"

namespace gaudi {
namespace calder {

template <typename TREE> void test_extents(const TREE &tree) {
  for (const auto &ext : tree.bvh_.internal) {
    vec4 c(0.5, 0.5, 0.1, 1.0);
    geometry_logger::ext(ext[0], ext[1], c);
  }
}

template <typename TREE>
void test_pyramid(const TREE &tree, const std::vector<index_t> &q_indices,
                  const std::vector<vec3> &q, const std::vector<real> &q_weights) {
  std::vector<vec3> wq(q);
  for (size_t i = 0; i < q.size(); i++)
    wq[i] *= q_weights[i];

  datum_t<vec3>::ptr x_datum = datum_t<vec3>::create(q_indices, wq);
  x_datum->pyramid(tree);
  for (size_t i = 0; i < tree.internal_nodes_.size(); i++) {
    const ext::extents_t &ext = tree.bvh_.internal[i];
    vec3 cen = 0.5 * (ext[0] + ext[1]);
    vec3 N = x_datum->node_data()[i];
    vec4 c(0.0, 0.5, 0.8, 1.0);
    geometry_logger::line(cen, cen + N, c);
  }
}

template <typename TREE>
void test_pyramid_scalar(const TREE &tree, const std::vector<index_t> &q_indices,
                         const std::vector<real> &q) {
  datum_t<real>::ptr x_datum = datum_t<real>::create(q_indices, q);
  x_datum->pyramid(tree);
  for (size_t i = 0; i < tree.internal_nodes_.size(); i++) {
    const ext::extents_t &ext = tree.bvh_.internal[i];
    real w = x_datum->node_data()[i];
    (void)w;
    (void)ext;
  }
}

// Shared Barnes-Hut opening traversal for a single query point.
//
// Walks the BVH using the BH opening test: sc < dc * eps, where sc is the
// equivalent sphere radius of the node's bounding box and dc is the distance
// from the query to the box center. "Far" nodes (test passes) call on_far;
// "near" paths descend until they reach leaf nodes, which call on_leaf.
//
// on_far:  void(index_t node_id, index_t query_id, const vec3 &query)
// on_leaf: void(index_t leaf_id, index_t orig_id, index_t query_id, const vec3 &query)
//
// The walk is sourced entirely through arp::tree_view (canonical radix node
// arrays + bounding boxes), so it is identical for any backend (Morton BVH or
// the modernized legacy AABB) that exposes that representation.
template <typename TREE, typename OnFar, typename OnLeaf>
void traverse_bh_opening(const TREE &tree, index_t query_id,
                         const vec3 &query, real eps,
                         OnFar &&on_far, OnLeaf &&on_leaf) {
  const arp::tree_view view = arp::make_tree_view(tree);
  const auto &internal_nodes = view.internal_nodes();

  if (internal_nodes.empty())
    return;

  std::vector<int> stack;
  stack.reserve(128);
  stack.push_back(0);

  while (!stack.empty()) {
    int node_id = stack.back();
    stack.pop_back();
    const auto &node = internal_nodes[node_id];
    if (node.split == arp::UNULL)
      continue;

    // Leaf children are single primitives -- always process them exactly.
    // Only apply the opening test when both children are internal subtrees.
    if (arp::left_leaf(node) || arp::right_leaf(node)) {
      if (arp::left_leaf(node)) {
        index_t leaf_id = node.start;
        index_t orig_id = view.get_index(leaf_id);
        on_leaf(leaf_id, orig_id, query_id, query);
      } else {
        stack.push_back(node.split);
      }

      if (arp::right_leaf(node)) {
        index_t leaf_id = node.end;
        index_t orig_id = view.get_index(leaf_id);
        on_leaf(leaf_id, orig_id, query_id, query);
      } else {
        stack.push_back(node.split + 1);
      }
    } else {
      const ext::extents_t &e = view.internal_bbox(node_id);
      vec3 de = e[1] - e[0];
      real V = de[0] * de[1] * de[2];
      real sc = pow(0.75 * V / M_PI, 1.0 / 3.0);
      vec3 pj = 0.5 * (e[0] + e[1]);
      real dc = va::dist(query, pj);

      if (sc < dc * eps) {
        on_far(node_id, query_id, query);
      } else {
        stack.push_back(node.split);
        stack.push_back(node.split + 1);
      }
    }
  }
}

// Visualize the Barnes-Hut traversal decisions as bounding boxes.
// "Far" internal nodes are drawn in far_color; leaf nodes reached by
// the "near" path are drawn in near_color.
template <typename TREE>
void log_bvh_barnes_hut(const TREE &tree, const std::vector<vec3> &queries,
                        real eps,
                        const vec4 &far_color = vec4(0.5, 0.5, 0.1, 1.0),
                        const vec4 &near_color = vec4(0.1, 0.8, 0.2, 1.0)) {
  const arp::tree_view view = arp::make_tree_view(tree);
  for (int qi = 0; qi < static_cast<int>(queries.size()); qi++) {
    traverse_bh_opening(
        tree, qi, queries[qi], eps,
        [&](index_t node_id, index_t, const vec3 & query) {
          const ext::extents_t &e = view.internal_bbox(node_id);
          vec3 cen = 0.5 * (e[0] + e[1]);
          geometry_logger::line(cen, query, far_color);
          geometry_logger::ext(e[0], e[1], far_color);
        },
        [&](index_t leaf_id, index_t, index_t, const vec3 &) {
          const ext::extents_t &e = view.leaf_bbox(leaf_id);
          vec3 cen = 0.5 * (e[0] + e[1]);
          geometry_logger::ext(e[0], e[1], near_color);
          geometry_logger::point(queries[qi], near_color);
        });
  }
}

enum Node_Type {
  LEAF,
  BRANCH,
};

template <typename T>
T get_data(Node_Type node_type, index_t j, index_t data_id,
           const std::vector<calder::datum::ptr> &data) {
  const typename calder::datum_t<T>::ptr F_datum =
      static_pointer_cast<typename calder::datum_t<T>>(data[data_id]);
  if (node_type == LEAF) {
    return F_datum->sorted_leaf_data()[j];
  } else {
    return F_datum->node_data()[j];
  }
}

template <typename TREE> class fast_summation {
public:
  typedef TREE Tree;
  using Node_Type = calder::Node_Type;

  fast_summation(const TREE &tree) : __tree(tree) {}

  void bind(const datum::ptr &x) { __data.push_back(x); }

  template <typename T>
  void bind(const std::vector<index_t> &indices, const std::vector<T> &x) {
    __data.push_back(datum_t<T>::create(indices, x));
  }

  template <typename Q>
  using ComputeFcn =
      std::function<Q(const index_t &, const index_t &, const vec3 &,
                      const std::vector<datum::ptr> &, Node_Type,
                      const Tree &)>;

  template <typename Q>
  std::vector<Q>
  calc(const std::vector<vec3> &pov, ComputeFcn<Q> leafComputeFcn,
       ComputeFcn<Q> nodeComputeFcn, real eps = 0.5, bool debug = false) {
    for (size_t i = 0; i < __data.size(); i++) {
      __data[i]->pyramid(__tree);
    }

    vector<Q> u(pov.size(), z::zero<Q>());

    if (__tree.internal_nodes_.empty())
      return u;

#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(pov.size()); i++) {
      traverse_bh_opening(
          __tree, i, pov[i], eps,
          [&](index_t node_id, index_t qi, const vec3 &pi) {
            u[qi] += nodeComputeFcn(qi, node_id, pi, __data, BRANCH, __tree);
          },
          [&](index_t leaf_id, index_t, index_t qi, const vec3 &pi) {
            u[qi] += leafComputeFcn(qi, leaf_id, pi, __data, LEAF, __tree);
          });
    }

    return u;
  }
  std::vector<datum::ptr> __data;
  const TREE &__tree;
};
} // namespace calder
} // namespace gaudi
#endif
