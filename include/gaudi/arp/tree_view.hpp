#ifndef __GAUDI_ARP_TREE_VIEW__
#define __GAUDI_ARP_TREE_VIEW__

#include "gaudi/arp/hash_tree.hpp"
#include <vector>

namespace gaudi {
namespace arp {

// ---------------------------------------------------------------------------
// tree_view: backend-agnostic window onto a binary acceleration structure.
// ---------------------------------------------------------------------------
// Both the Morton/radix BVH (bvh_tree / bvh_tree_t) and the (to-be) modernized
// legacy AABB tree expose the SAME canonical representation:
//
//   indices_         sorted-leaf id  -> original primitive id (permutation)
//   internal_nodes_  radix_tree_node[] with parent/child links  (N-1 nodes)
//   leaf_nodes_      radix_tree_node[] with parent links        (N   leaves)
//   bvh_             per-node extents (leaf[] + internal[]) for bounding boxes
//
// Consumers in calder (pyramid / datums / Barnes-Hut traversal / integrators)
// talk to this view instead of reaching into a specific tree's members, so the
// backend can be swapped behind USE_HASH without touching the FMM code.
//
// A binary tree over N singleton leaves has exactly N leaves and N-1 internal
// nodes, so internal-node ids and leaf ids each index a contiguous space that
// lines up 1:1 with the datum arrays.
//
// The bbox (bvh_) reference is optional: pyramid construction only needs the
// parent-walk surface, so a view can be built from the three arrays alone.
struct tree_view {
  const std::vector<index_t> *indices_ = nullptr;
  const std::vector<radix_tree_node> *internal_ = nullptr;
  const std::vector<radix_tree_node> *leaf_ = nullptr;
  const TreeResult<ext::extents_t> *bvh_ = nullptr;

  tree_view() = default;
  tree_view(const std::vector<index_t> &indices,
            const std::vector<radix_tree_node> &internal_nodes,
            const std::vector<radix_tree_node> &leaf_nodes,
            const TreeResult<ext::extents_t> *bvh = nullptr)
      : indices_(&indices), internal_(&internal_nodes), leaf_(&leaf_nodes),
        bvh_(bvh) {}

  // --- structural ---------------------------------------------------------
  size_t leaf_count() const { return leaf_->size(); }
  size_t internal_count() const { return internal_->size(); }
  bool empty() const { return leaf_->empty(); }

  // Parent links (UNULL at the root). These drive the bottom-up pyramid.
  index_t leaf_parent(index_t leaf_id) const {
    return (*leaf_)[leaf_id].parent;
  }
  index_t internal_parent(index_t node_id) const {
    return (*internal_)[node_id].parent;
  }

  // sorted-leaf id -> original primitive id.
  index_t index(index_t sorted_leaf_id) const {
    return (*indices_)[sorted_leaf_id];
  }
  index_t get_index(index_t sorted_leaf_id) const {
    return (*indices_)[sorted_leaf_id];
  }

  // --- raw array access (for traversal adapters / integrators) ------------
  const std::vector<index_t> &indices() const { return *indices_; }
  const std::vector<radix_tree_node> &internal_nodes() const {
    return *internal_;
  }
  const std::vector<radix_tree_node> &leaf_nodes() const { return *leaf_; }

  // --- bounding boxes (require the optional bvh_ reference) ---------------
  bool has_bbox() const { return bvh_ != nullptr; }
  const ext::extents_t &leaf_bbox(index_t leaf_id) const {
    return bvh_->leaf[leaf_id];
  }
  const ext::extents_t &internal_bbox(index_t node_id) const {
    return bvh_->internal[node_id];
  }
};

// Build a tree_view from any tree exposing the canonical members. Works for
// both arp::bvh_tree<N> and arp::bvh_tree_t<SimplexType>.
template <typename TreeT>
inline tree_view make_tree_view(const TreeT &tree) {
  return tree_view(tree.indices_, tree.internal_nodes_, tree.leaf_nodes_,
                   &tree.bvh_);
}

// ---------------------------------------------------------------------------
// Generic build_pyramid over a tree_view.
// ---------------------------------------------------------------------------
// Identical bottom-up map/reduce as the array-based build_pyramid in
// hash_tree.hpp, but driven by the view's parent links so it is independent of
// the backing tree type. See the array overload for the monoid requirements.
//
//   Q0 = leaf input type, Q1 = node output type
//   leaf_map(const Q0&, const Q1&) -> Q1
//   node_reduce(const Q1&, const Q1&) -> Q1
template <typename Q0, typename Q1>
inline std::vector<Q1> build_pyramid(const tree_view &view,
                                     const std::vector<Q0> &leaf_data,
                                     auto &&leaf_map, auto &&node_reduce,
                                     const Q1 &identity) {
  std::vector<Q1> internal_reduce(view.internal_count(), identity);
  for (size_t i = 0; i < leaf_data.size(); i++) {
    Q1 mapped = leaf_map(leaf_data[i], identity);
    index_t depth = 0;
    index_t parent = view.leaf_parent(static_cast<index_t>(i));
    while (depth < 64 && parent != UNULL) {
      internal_reduce[parent] = node_reduce(mapped, internal_reduce[parent]);
      parent = view.internal_parent(parent);
      depth++;
    }
  }
  return internal_reduce;
}

// Convenience: single-type pyramid where Q0 == Q1 and leaf_map == reduce.
template <TypeArray TTYPE>
inline TTYPE build_pyramid(const tree_view &view, const TTYPE &data,
                           auto &&reduce_func,
                           const typename TTYPE::value_type &default_val) {
  using O = typename TTYPE::value_type;
  return build_pyramid<O, O>(
      view, data,
      [&](const O &a, const O &b) -> O { return reduce_func(a, b); },
      [&](const O &a, const O &b) -> O { return reduce_func(a, b); },
      default_val);
}

} // namespace arp
} // namespace gaudi

#endif // __GAUDI_ARP_TREE_VIEW__
