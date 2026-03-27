#ifndef __GAUDI_ARP_HASH_TREE__
#define __GAUDI_ARP_HASH_TREE__

#include "gaudi/arp/morton.hpp"
#include "gaudi/arp/pairwise_tests.hpp"
#include "gaudi/arp/simplex_set.hpp"
#include "gaudi/common.h"
#include "gaudi/console_logger.hpp"
#include "gaudi/geometry_logger.hpp"
#include "gaudi/geometry_types.hpp"
#include <algorithm>
#include <functional>
#include <iostream>
#include <optional>
#include <queue>
#include <stack>
#include <tuple>
#include <vector>

#ifdef _MSC_VER
#include <intrin.h>
#endif

namespace gaudi {
namespace arp {

const index_t UNULL = -1;

struct radix_tree_node {
  index_t start = UNULL;
  index_t end = UNULL;
  index_t split =
      UNULL; // split is the index of the last element in the left child
  index_t parent = UNULL;
};

// Generic result type for tree operations returning leaf and internal data
template <typename T> struct TreeResult {
  std::vector<T> leaf;
  std::vector<T> internal;
};

// Convenient type alias for node results
using NodeResult = TreeResult<radix_tree_node>;

// Common leading zero calculation for tree building.
template <typename MortonT>
inline int clz(index_t i, index_t j, const std::vector<MortonT> &hash) {
  if (j < 0)
    return -1;
  if (j > static_cast<index_t>(hash.size()) - 1)
    return -1;
  MortonT diff = xor_morton(hash[i], hash[j]);
  if (is_zero(diff))
    return static_cast<int>(MortonT::total_bits);
  return clz_morton(diff);
}

// Find split point in radix tree
template <typename MortonT>
inline index_t find_split(index_t start, index_t end,
                          const std::vector<MortonT> &hash) {
  int common_prefix_dist = clz(start, end, hash);
  index_t split = start;

  index_t step = end - start;
  while (step > 1) {
    step = (step + 1) >> 1; // exponential decrease
    index_t new_split = split + step;
    if (new_split < end) {
      int new_prefix_dist = clz(start, new_split, hash);
      if (new_prefix_dist > common_prefix_dist) {
        split = new_split;
      }
    }
  }

  return split;
}

// Find range for radix tree construction
template <typename MortonT>
inline std::pair<index_t, index_t>
find_range(index_t i, const std::vector<MortonT> &hash) {
  index_t N = static_cast<index_t>(hash.size());

  int dir = (clz(i, i + 1, hash) - clz(i, i - 1, hash)) > 0 ? 1 : -1;
  int sig_min = clz(i, i - dir, hash);

  index_t lmax = 2;
  while (clz(i, i + lmax * dir, hash) > sig_min) {
    lmax *= 2;
  }

  index_t l = 0;
  index_t t = lmax;

  while (t >= 1) {
    t = t >> 1;
    if (clz(i, i + (l + t) * dir, hash) > sig_min) {
      l += t;
    }
  }
  index_t j = i + l * dir;

  return dir < 0 ? std::make_pair(j, i) : std::make_pair(i, j);
}

// Build radix tree from sorted indices and hashes
template <typename MortonT>
inline NodeResult build_tree(const std::vector<MortonT> &hash) {
  std::vector<radix_tree_node> internal_nodes(hash.size() - 1);
  std::vector<radix_tree_node> leaf_nodes(hash.size());

  // Initialize leaf nodes
  for (index_t i = 0; i < static_cast<index_t>(hash.size()); i++) {
    leaf_nodes[i].start = i;
    leaf_nodes[i].end = i + 1;
    leaf_nodes[i].split = UNULL;
    leaf_nodes[i].parent = UNULL;
  }

  // Build internal nodes
  for (index_t i = 0; i < static_cast<index_t>(hash.size()) - 1; i++) {
    auto range = find_range(i, hash);
    index_t split = find_split(range.first, range.second, hash);

    internal_nodes[i].start = range.first;
    internal_nodes[i].end = range.second;
    internal_nodes[i].split = split;

    if (split == range.first) {
      leaf_nodes[range.first].parent = i;
    } else {
      internal_nodes[split].parent = i;
    }

    if (split + 1 == range.second) {
      leaf_nodes[range.second].parent = i;
    } else {
      internal_nodes[split + 1].parent = i;
    }
  }

  return {std::move(leaf_nodes), std::move(internal_nodes)};
}

// Test tree construction and traversal
template <typename MortonT>
inline void test_tree(const std::vector<radix_tree_node> &internal_nodes,
                      const std::vector<radix_tree_node> &leaf_nodes,
                      const std::vector<index_t> &ids,
                      const std::vector<MortonT> &hash) {
  std::vector<bool> visited(ids.size(), false);
  std::stack<index_t> stack;
  stack.push(0);

  while (!stack.empty()) {
    index_t i = stack.top();
    stack.pop();

    if (internal_nodes[i].split + 0 == internal_nodes[i].start ||
        internal_nodes[i].split + 1 == internal_nodes[i].end) {
      index_t j0 = internal_nodes[i].start;
      index_t j1 = internal_nodes[i].end;
      visited[j0] = true;
      visited[j1] = true;
      assert(leaf_nodes[j0].parent == i);
      assert(leaf_nodes[j1].parent == i);
    } else {
      index_t j0 = internal_nodes[i].split;
      index_t j1 = internal_nodes[i].split + 1;
      // Assert parents
      assert(internal_nodes[j0].parent == i);
      assert(internal_nodes[j1].parent == i);

      // Assert ranges
      assert(internal_nodes[j0].start >= internal_nodes[i].start);
      assert(internal_nodes[j0].end <= internal_nodes[i].end);
      assert(internal_nodes[j1].start >= internal_nodes[i].start);
      assert(internal_nodes[j1].end <= internal_nodes[i].end);

      stack.push(j0);
      stack.push(j1);
    }
  }

  // Check if all leaves were visited
  bool visited_all = true;
  for (bool v : visited) {
    if (!v) {
      visited_all = false;
      break;
    }
  }
  assert(visited_all);
}

// Unit test for tree construction
inline void unit_test_tree() {
  std::vector<morton_t> hash = {
      morton_t::from_uint64(0), morton_t::from_uint64(1), morton_t::from_uint64(2),
      morton_t::from_uint64(3), morton_t::from_uint64(4), morton_t::from_uint64(5),
      morton_t::from_uint64(6), morton_t::from_uint64(7), morton_t::from_uint64(8),
      morton_t::from_uint64(9), morton_t::from_uint64(10), morton_t::from_uint64(11)};
  std::vector<index_t> ids = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
  auto [leaf_nodes, internal_nodes] = build_tree(hash);
  test_tree(internal_nodes, leaf_nodes, ids, hash);
}

inline bool has_leaves(const radix_tree_node &node) {
  return node.split == node.start || node.split + 1 == node.end;
}

inline bool both_leaves(const radix_tree_node &node) {
  return node.split == node.start && node.split + 1 == node.end;
}

inline bool left_leaf(const radix_tree_node &node) {
  return node.split == node.start;
}

inline bool right_leaf(const radix_tree_node &node) {
  return node.split + 1 == node.end;
}

template <typename C>
auto next(C &c) -> decltype(auto) {
  if constexpr (requires { c.top(); })
    return c.top();
  else
    return c.front();
}

template <typename Container, typename FuncNode, typename FuncLeaf,
          typename FuncContinue>
inline void traverse_impl(const std::vector<radix_tree_node> &internal_nodes,
                          const std::vector<radix_tree_node> &leaf_nodes,
                          FuncNode &&process_node, FuncLeaf &&process_leaf,
                          FuncContinue &&should_continue) {
  if (internal_nodes.empty())
    return;

  Container container;
  container.push(0);

  std::vector<bool> visited(internal_nodes.size(), false);

  while (!container.empty()) {
    index_t cid = next(container);
    container.pop();


    const auto &cnode = internal_nodes[cid];
    if (cnode.split == UNULL)
      continue;

    process_node(cid, cnode);

    if (has_leaves(cnode)) {
      if (both_leaves(cnode)) {
        process_leaf(cid, cnode.start, leaf_nodes[cnode.start]);
        process_leaf(cid, cnode.end, leaf_nodes[cnode.end]);
      } else if (left_leaf(cnode)) {
        process_leaf(cid, cnode.start, leaf_nodes[cnode.start]);
        container.push(cnode.split + 1);
      } else if (right_leaf(cnode)) {
        process_leaf(cid, cnode.end, leaf_nodes[cnode.end]);
        container.push(cnode.split);
      }
    } else if(should_continue(cid, cnode)){
      container.push(cnode.split);
      container.push(cnode.split + 1);
    }
  }
}

template <typename FuncNode, typename FuncLeaf, typename FuncContinue>
inline void traverse_dfs(const std::vector<radix_tree_node> &internal_nodes,
                         const std::vector<radix_tree_node> &leaf_nodes,
                         FuncNode &&process_node, FuncLeaf &&process_leaf,
                         FuncContinue &&should_continue) {
  traverse_impl<std::stack<index_t>>(
      internal_nodes, leaf_nodes, std::forward<FuncNode>(process_node),
      std::forward<FuncLeaf>(process_leaf),
      std::forward<FuncContinue>(should_continue));
}

template <typename FuncNode, typename FuncLeaf, typename FuncContinue>
inline void traverse_bfs(const std::vector<radix_tree_node> &internal_nodes,
                         const std::vector<radix_tree_node> &leaf_nodes,
                         FuncNode &&process_node, FuncLeaf &&process_leaf,
                         FuncContinue &&should_continue) {
  traverse_impl<std::queue<index_t>>(
      internal_nodes, leaf_nodes, std::forward<FuncNode>(process_node),
      std::forward<FuncLeaf>(process_leaf),
      std::forward<FuncContinue>(should_continue));
}

// Best-first traversal using a priority queue.
// compute_priority(child_id, child_node) -> real:
//   negative = prune, otherwise priority (lower = explored first)
template <typename FuncNode, typename FuncLeaf, typename FuncPriority>
inline void traverse_best(const std::vector<radix_tree_node> &internal_nodes,
                          const std::vector<radix_tree_node> &leaf_nodes,
                          FuncNode &&process_node, FuncLeaf &&process_leaf,
                          FuncPriority &&compute_priority) {
  if (internal_nodes.empty())
    return;

  using pq_entry = std::pair<real, index_t>;
  std::priority_queue<pq_entry, std::vector<pq_entry>, std::greater<pq_entry>>
      pq;
  pq.push({0.0, 0});

  std::vector<bool> visited(internal_nodes.size(), false);

  while (!pq.empty()) {
    auto [pri, cid] = pq.top();
    pq.pop();

    if (visited[cid])
      continue;
    visited[cid] = true;

    const auto &cnode = internal_nodes[cid];
    if (cnode.split == UNULL)
      continue;

    process_node(cid, cnode);

    auto push_child = [&](index_t child_id) {
      real p = compute_priority(child_id, internal_nodes[child_id]);
      if (p >= 0.0)
        pq.push({p, child_id});
    };

    if (has_leaves(cnode)) {
      if (both_leaves(cnode)) {
        process_leaf(cid, cnode.start, leaf_nodes[cnode.start]);
        process_leaf(cid, cnode.end, leaf_nodes[cnode.end]);
      } else if (left_leaf(cnode)) {
        process_leaf(cid, cnode.start, leaf_nodes[cnode.start]);
        push_child(cnode.split + 1);
      } else if (right_leaf(cnode)) {
        process_leaf(cid, cnode.end, leaf_nodes[cnode.end]);
        push_child(cnode.split);
      }
    } else {
      push_child(cnode.split);
      push_child(cnode.split + 1);
    }
  }
}

// Build pyramid using map-reduce pattern
// if input type is != output type, then you  need to map
// it first.
// takes an input vector of data, maps it to an output type using a map
// function, then walks the data up the tree and reduces it using a reduce
// function
template <TypeArray TTYPE>
inline TTYPE build_pyramid(const TTYPE &data,
                           const std::vector<radix_tree_node> &internal_nodes,
                           const std::vector<radix_tree_node> &leaf_nodes,
                           auto &&reduce_func,
                           const typename TTYPE::value_type &default_val) {

  using O = typename TTYPE::value_type;
  auto thread_safe_reduce = [&](const typename TTYPE::value_type &a,
                                const typename TTYPE::value_type &b) {
    // lock here
    return reduce_func(a, b);
    // unlock here
  };
  std::vector<O> internalReduce(data.size() - 1, default_val);
  for (int i = 0; i < data.size(); i++) {
    const O &datai = data[i];
    index_t j = 0;
    index_t parent = leaf_nodes[i].parent;
    while (j < 64 && parent != UNULL) {
      const O &dataj = internalReduce[parent];
      internalReduce[parent] = thread_safe_reduce(datai, dataj);
      parent = internal_nodes[parent].parent;
      j++;
    }
  }

  return internalReduce;
}

// Build complete hash tree with pyramid (Morton-key typed).
template <typename MortonT = morton_t, Vec3View TTYPE>
inline std::tuple<std::vector<MortonT>, std::vector<index_t>,
                  std::vector<radix_tree_node>, std::vector<radix_tree_node>>
make_hash_tree(const TTYPE &data) {
  auto [hashes, indices] = make_hash_3d_t<MortonT>(data);
  if (hashes.empty()) {
    return {{}, {}, {}, {}};
  }

  auto [leaf_nodes, internal_nodes] = build_tree(hashes);
  return std::make_tuple(hashes, indices, internal_nodes, leaf_nodes);
}

// make_bvh for flat Vec3View (legacy)
template <int N, Vec3View TTYPE>
inline TreeResult<ext::extents_t>
make_bvh(const TTYPE &data, const std::vector<radix_tree_node> &internal_nodes,
         const std::vector<radix_tree_node> &leaf_nodes) {

  if (data.size() % N != 0) {
    throw std::runtime_error("Data size must be a multiple of " +
                             std::to_string(N));
  }

  const auto exts = calc_extents<N>(data);
  const auto default_val = ext::init();
  auto reduce_function = [&](const ext::extents_t &a, const ext::extents_t &b) {
    ext::extents_t c = ext::expand(b, a);
    return ext::expand(b, a);
  };

  auto internalReduce = build_pyramid(exts, internal_nodes, leaf_nodes,
                                      reduce_function, default_val);

  return {
      exts,
      internalReduce,
  };
}

// make_bvh for SimplexView (type-based)
// Extracts stride from the SimplexView type automatically
template <SimplexView STYPE>
inline TreeResult<ext::extents_t>
make_bvh(const STYPE &data, const std::vector<radix_tree_node> &internal_nodes,
         const std::vector<radix_tree_node> &leaf_nodes) {

  const auto exts = calc_extents(data);
  const auto default_val = ext::init();
  auto reduce_function = [&](const ext::extents_t &a, const ext::extents_t &b) {
    return ext::expand(b, a);
  };

  auto internalReduce = build_pyramid(exts, internal_nodes, leaf_nodes,
                                      reduce_function, default_val);

  return {
      exts,
      internalReduce,
  };
}

// make_points for flat Vec3View (legacy)
template <int N, Vec3View TTYPE>
inline TreeResult<vec3>
make_points(const TTYPE &data,
            const std::vector<radix_tree_node> &internal_nodes,
            const std::vector<radix_tree_node> &leaf_nodes) {
  const auto default_val = MassPoint{0.0f, vec3::Zero()};
  const auto coms = calc_com<N>(data); // Specify template parameter N=1

  auto reduce_function = [&](const MassPoint &a,
                             const MassPoint &b) -> MassPoint {
    const real &ma = std::get<0>(a);
    const real &mb = std::get<0>(b);
    const vec3 &pa = ma * std::get<1>(a);
    const vec3 &pb = std::get<1>(b); // accumulate point sum
    return MassPoint{ma + mb, pa + pb};
  };

  // Build pyramid with map-reduce
  auto internalReduce = build_pyramid(coms, internal_nodes, leaf_nodes,
                                      reduce_function, default_val);

  std::vector<vec3> leaf_points(coms.size());
  for (int i = 0; i < coms.size(); i++) {
    leaf_points[i] = std::get<1>(coms[i]); // average point
  }

  std::vector<vec3> internal_points(internalReduce.size());
  for (int i = 0; i < internalReduce.size(); i++) {
    const vec3 &p = std::get<1>(internalReduce[i]);
    const real &m = std::get<0>(internalReduce[i]);
    internal_points[i] = p / m; // average point
  }
  return {
      leaf_points,    // leaf points ordered by indices
      internal_points // internal points
  };
}

// make_points for SimplexView (type-based)
// Extracts stride from the SimplexView type automatically
template <SimplexView STYPE>
inline TreeResult<vec3>
make_points(const STYPE &data,
            const std::vector<radix_tree_node> &internal_nodes,
            const std::vector<radix_tree_node> &leaf_nodes) {
  const auto default_val = MassPoint{0.0f, vec3::Zero()};
  const auto coms = calc_com(data);

  auto reduce_function = [&](const MassPoint &a,
                             const MassPoint &b) -> MassPoint {
    const real &ma = std::get<0>(a);
    const real &mb = std::get<0>(b);
    const vec3 &pa = ma * std::get<1>(a);
    const vec3 &pb = std::get<1>(b); // accumulate point sum
    return MassPoint{ma + mb, pa + pb};
  };

  // Build pyramid with map-reduce
  auto internalReduce = build_pyramid(coms, internal_nodes, leaf_nodes,
                                      reduce_function, default_val);

  std::vector<vec3> leaf_points(coms.size());
  for (size_t i = 0; i < coms.size(); i++) {
    leaf_points[i] = std::get<1>(coms[i]); // average point
  }

  std::vector<vec3> internal_points(internalReduce.size());
  for (size_t i = 0; i < internalReduce.size(); i++) {
    const vec3 &p = std::get<1>(internalReduce[i]);
    const real &m = std::get<0>(internalReduce[i]);
    internal_points[i] = p / m; // average point
  }
  return {
      leaf_points,    // leaf points ordered by indices
      internal_points // internal points
  };
}

// log_hierarchy for SimplexView types (tuple-based)
template <SimplexView STYPE>
void log_hierarchy(const STYPE &data,
                   const std::vector<radix_tree_node> &internal_nodes,
                   const std::vector<radix_tree_node> &leaf_nodes) {
  constexpr size_t N = simplex_stride_v<STYPE>;
  
  if (internal_nodes.empty() || leaf_nodes.empty()) {
    console_logger::debug << "No internal or leaf nodes to log" << std::endl;
    return;
  }
  // Compute centerpoints for each stride group for visualization
  const auto [leaf_points, internal_points] =
      make_points(data, internal_nodes, leaf_nodes);
  // draw the simplex for logging (SimplexView provides tuples)
  if constexpr (N > 1) {
    for (size_t i = 0; i < data.size(); i++) {
      auto simplex = data[i];
      for (size_t j = 0; j < N; j++) {
        const size_t j1 = (j + 1) % N;
        geometry_logger::line(simplex[j], simplex[j1], vec4(0.2, 0.2, 0.2, 0.5f));
      }
    }
  }
  // draw centerpoints
  for (size_t i = 0; i < leaf_points.size(); i++) {
    geometry_logger::point(leaf_points[i], vec4(0.0, 1.0, 0.0, 1.0f));
  }
  // Continue with original logic using make_points
  traverse_bfs(
      internal_nodes, leaf_nodes,
      [&](index_t node_id, const radix_tree_node &node) {
        if (node.parent != UNULL) {
          const auto &ppoint = internal_points[node.parent];
          const auto &cpoint = internal_points[node_id];
          geometry_logger::line(ppoint, cpoint, vec4(0.0, 1.0, 0.0, 1.0f));
        }
      },
      [&](index_t parent_id, index_t leaf_id, const radix_tree_node &node) {
        if (node.parent != UNULL) {
          const auto &ppoint = internal_points[node.parent];
          const auto &cpoint = leaf_points[leaf_id];
          geometry_logger::line(ppoint, cpoint, vec4(0.0, 0.0, 1.0, 1.0f));
          geometry_logger::point(ppoint, vec4(1.0, 0.0, 1.0, 1.0f));
        }
      },
      [&](index_t, const radix_tree_node &) { return true; });
}

// log_bvh for SimplexView types (tuple-based)
template <SimplexView STYPE>
void log_bvh(const STYPE &data,
             const std::vector<radix_tree_node> &internal_nodes,
             const std::vector<radix_tree_node> &leaf_nodes) {
  constexpr size_t N = simplex_stride_v<STYPE>;

  const auto bvh_result = make_bvh(data, internal_nodes, leaf_nodes);
  const auto &leaf_bvh = bvh_result.leaf;
  const auto &internal_bvh = bvh_result.internal;
  
  // Draw all vertices from simplices
  for (size_t i = 0; i < data.size(); i++) {
    auto simplex = data[i];
    for (size_t j = 0; j < N; j++) {
      geometry_logger::point(simplex[j], vec4(1.0, 0.0, 0.0, 1.0));
    }
  }
  for (size_t i = 0; i < internal_bvh.size(); i++) {
    const ext::extents_t &ext = internal_bvh[i];
    geometry_logger::ext(ext[0], ext[1], vec4(0.0, 1.0, 0.0, 0.5));
  }
}

template <size_t N>
inline ext::extents_t calc_simplex_extents(const std::array<vec3, N> &simplex) {
  auto out = ext::init();
  for (size_t i = 0; i < N; ++i) {
    out = ext::expand(out, simplex[i]);
  }
  return out;
}

// BFS contracting radius: find single closest element
template <SimplexView PTYPE, SimplexView STYPE>
index_t get_nearest(const typename PTYPE::value_type &prim, const STYPE &data,
                    const std::vector<radix_tree_node> &internal_nodes,
                    const std::vector<radix_tree_node> &leaf_nodes,
                    const TreeResult<ext::extents_t> &bvh_result,
                    auto &&testAB) {

  constexpr size_t Nprim = simplex_stride_v<PTYPE>;

  index_t idMin = -1;
  real mMin = std::numeric_limits<real>::max();
  ext::extents_t ext_t = calc_simplex_extents<Nprim>(prim);
  ext_t = ext::inflate(ext_t, mMin);

  traverse_bfs(
      internal_nodes, leaf_nodes,
      [&](index_t, const radix_tree_node &) {},
      [&](index_t, index_t leaf_id, const radix_tree_node &) {
        auto datum = data[leaf_id];
        real dist = testAB(prim, datum);
        if (dist < mMin) {
          mMin = dist;
          idMin = data.get_index(leaf_id);
        }
        ext_t = ext::inflate(calc_simplex_extents<Nprim>(prim), mMin);
      },
      [&](index_t node_id, const radix_tree_node &) -> bool {
        const ext::extents_t &ext_s = bvh_result.internal[node_id];
        ext_t = ext::inflate(calc_simplex_extents<Nprim>(prim), mMin);
        return ext::overlap(ext_t, ext_s);
      });

  return idMin;
}

// Best-first contracting radius: find single closest element
template <SimplexView PTYPE, SimplexView STYPE>
index_t get_nearest_best(const typename PTYPE::value_type &prim,
                         const STYPE &data,
                         const std::vector<radix_tree_node> &internal_nodes,
                         const std::vector<radix_tree_node> &leaf_nodes,
                         const TreeResult<ext::extents_t> &bvh_result,
                         auto &&testAB) {

  constexpr size_t Nprim = simplex_stride_v<PTYPE>;

  index_t idMin = -1;
  real mMin = std::numeric_limits<real>::max();
  vec3 query_center = calc_simplex_extents<Nprim>(prim)[0];
  for (size_t i = 0; i < Nprim; ++i)
    query_center = (query_center + prim[i]) / 2.0;

  traverse_best(
      internal_nodes, leaf_nodes,
      [&](index_t, const radix_tree_node &) {},
      [&](index_t, index_t leaf_id, const radix_tree_node &) {
        auto datum = data[leaf_id];
        real dist = testAB(prim, datum);
        if (dist < mMin) {
          mMin = dist;
          idMin = data.get_index(leaf_id);
        }
      },
      [&](index_t child_id, const radix_tree_node &) -> real {
        const ext::extents_t &ext_s = bvh_result.internal[child_id];
        real d = ext::min_distance_to_aabb(query_center, ext_s);
        return (d < mMin) ? d : -1.0;
      });

  return idMin;
}

// BFS fixed radius: collect all elements within tolerance
template <SimplexView PTYPE, SimplexView STYPE>
std::vector<index_t>
get_neighbors(const typename PTYPE::value_type &prim, const STYPE &data,
              const std::vector<radix_tree_node> &internal_nodes,
              const std::vector<radix_tree_node> &leaf_nodes,
              const TreeResult<ext::extents_t> &bvh_result, real tol,
              auto &&testAB) {

  constexpr size_t Nprim = simplex_stride_v<PTYPE>;

  ext::extents_t ext_t = calc_simplex_extents<Nprim>(prim);
  ext_t = ext::inflate(ext_t, tol);
  std::vector<index_t> collisions;

  traverse_bfs(
      internal_nodes, leaf_nodes,
      [&](index_t, const radix_tree_node &) {},
      [&](index_t, index_t leaf_id, const radix_tree_node &) {
        auto datum = data[leaf_id];
        real dist = testAB(prim, datum);
        if (dist < tol)
          collisions.push_back(data.get_index(leaf_id));
      },
      [&](index_t node_id, const radix_tree_node &) -> bool {
        const ext::extents_t &ext_s = bvh_result.internal[node_id];
        return ext::overlap(ext_t, ext_s);
      });

  return collisions;
}

// make_hash for flat Vec3View (legacy)
template <int N, Vec3View TTYPE, typename MortonT = morton_t>
inline std::tuple<std::vector<MortonT>, std::vector<index_t>,
                  std::vector<radix_tree_node>, std::vector<radix_tree_node>>
make_hash(const TTYPE &data) {
  const auto mass_points = calc_com<N>(data);
  std::vector<vec3> averaged;
  averaged.reserve(mass_points.size());
  for (const auto &[mass, com] : mass_points) {
    averaged.push_back(com);
  }
  return make_hash_tree<MortonT>(averaged);
}

// make_hash for SimplexView (type-based)
template <typename MortonT = morton_t, SimplexView STYPE>
inline std::tuple<std::vector<MortonT>, std::vector<index_t>,
                  std::vector<radix_tree_node>, std::vector<radix_tree_node>>
make_hash(const STYPE &data) {
  const auto mass_points = calc_com(data);
  std::vector<vec3> averaged;
  averaged.reserve(mass_points.size());
  for (const auto &[mass, com] : mass_points) {
    averaged.push_back(com);
  }
  return make_hash_tree<MortonT>(averaged);
}

template <int N, typename MortonT = morton_t>
inline std::tuple<std::vector<MortonT>, std::vector<index_t>,
                  std::vector<radix_tree_node>, std::vector<radix_tree_node>>
make_hash(const simplex_set<N> &set) {
  return make_hash_tree<MortonT>(set.centroids());
}

template <int N>
inline simplex_set<N> make_simplex_set(const std::vector<vec3> &vertices,
                                       const std::vector<index_t> &adjacency) {
  return simplex_set<N>(vertices, adjacency);
}

// New bvh_tree templated on SimplexType
// This version extracts N from the SimplexType automatically
template <SimplexView SimplexType, typename MortonT = morton_t>
class bvh_tree_t {
public:
  static constexpr size_t N = simplex_stride_v<SimplexType>;
  
  // Member variables
  std::vector<index_t> indices_;
  std::vector<radix_tree_node> internal_nodes_;
  std::vector<radix_tree_node> leaf_nodes_;
  TreeResult<ext::extents_t> bvh_;
  std::vector<vec3> data_;
  std::vector<index_t> adjacency_;
  std::vector<MassPoint> coms_;
  std::vector<MortonT> hashes_;

  // Type aliases
  using ptr = std::shared_ptr<bvh_tree_t<SimplexType, MortonT>>;
  using simplex_view_type = SimplexType;
  using permutation_index_type = std::vector<index_t>;
  using permuted_view_type = permuted_simplex_view<N, std::vector<vec3>, std::vector<index_t>>;

  // Optional views - constructed after data is available
  std::optional<permuted_view_type> permuted_data_view_;

  static ptr create(const std::vector<index_t> &adjacency,
                    const std::vector<vec3> &vertices, int lvl = 8) {
    return std::make_shared<bvh_tree_t<SimplexType, MortonT>>(vertices, adjacency);
  }

  bvh_tree_t(const std::vector<vec3> &data,
             const std::vector<index_t> &adjacency) {
    update(data, adjacency);
  }

  index_t get_index(size_t i) const { return indices_[i]; }

  void update(const std::vector<vec3> &data,
              const std::vector<index_t> &adjacency) {
    data_ = data;
    adjacency_ = adjacency;
    
    // Create initial simplex view to compute hash
    permuted_simplex_view<N, std::vector<vec3>, std::vector<index_t>> 
        initial_view(data_, adjacency_, std::vector<index_t>{});
    
    // We need an unpermuted view first to compute the hash
    // Create identity permutation
    std::vector<index_t> identity(adjacency_.size() / N);
    for (size_t i = 0; i < identity.size(); ++i) {
      identity[i] = static_cast<index_t>(i);
    }
    
    // Create unpermuted simplex view
    permuted_simplex_view<N, std::vector<vec3>, std::vector<index_t>>
        unpermuted_view(data_, adjacency_, identity);
    
    auto [hashes, indices, internal_nodes, leaf_nodes] = make_hash<MortonT>(unpermuted_view);
    
    // Copy to members FIRST so the view can reference stable storage
    indices_ = std::move(indices);
    internal_nodes_ = std::move(internal_nodes);
    leaf_nodes_ = std::move(leaf_nodes);
    hashes_ = std::move(hashes);

    // Now create the permuted view pointing to member storage
    permuted_data_view_.emplace(data_, adjacency_, indices_);

    bvh_ = make_bvh(*permuted_data_view_, internal_nodes_, leaf_nodes_);
    coms_ = calc_com(*permuted_data_view_);
  }

  std::array<index_t, N> get_tuple_ids(const index_t &i) {
    return permuted_data_view_->get_tuple_ids(i);
  }

  // Get center of mass for a leaf node
  vec3 get_com(index_t i) const {
    return std::get<1>(coms_[i]);
  }

  template <SimplexView QueryType>
  auto dispatch_test() const {
    constexpr size_t Nq = simplex_stride_v<QueryType>;
    if constexpr (Nq == 1 && N == 1)
      return [](const auto &q, const auto &d) { return test_point_point_tuple(q, d); };
    else if constexpr (Nq == 1 && N == 2)
      return [](const auto &q, const auto &d) { return test_point_line_tuple(q, d); };
    else if constexpr (Nq == 1 && N == 3)
      return [](const auto &q, const auto &d) { return test_point_tri_tuple(q, d); };
    else if constexpr (Nq == 2 && N == 1)
      return [](const auto &q, const auto &d) { return test_point_line_tuple(d, q); };
    else if constexpr (Nq == 2 && N == 2)
      return [](const auto &q, const auto &d) { return test_line_line_tuple(q, d); };
    else if constexpr (Nq == 2 && N == 3)
      return [](const auto &q, const auto &d) { return test_line_tri_tuple(q, d); };
    else if constexpr (Nq == 3 && N == 1)
      return [](const auto &q, const auto &d) { return test_point_tri_tuple(d, q); };
    else if constexpr (Nq == 3 && N == 2)
      return [](const auto &q, const auto &d) { return test_line_tri_tuple(d, q); };
    else if constexpr (Nq == 3 && N == 3)
      return [](const auto &q, const auto &d) { return test_tri_tri_tuple(q, d); };
  }

  template <SimplexView QueryType>
  index_t find_nearest(const typename QueryType::value_type &query) {
    return arp::get_nearest<QueryType, permuted_view_type>(
        query, *permuted_data_view_, internal_nodes_, leaf_nodes_, bvh_,
        dispatch_test<QueryType>());
  }

  template <SimplexView QueryType>
  std::vector<index_t> find_neighbors(const typename QueryType::value_type &query, real tol) {
    return arp::get_neighbors<QueryType, permuted_view_type>(
        query, *permuted_data_view_, internal_nodes_, leaf_nodes_, bvh_, tol,
        dispatch_test<QueryType>());
  }
};

// Legacy bvh_tree templated on N - now uses permuted_simplex_view internally
template <int N, typename MortonT = morton_t>
class bvh_tree {
  public:
    // Member variables first
    std::vector<index_t> indices_;
    std::vector<radix_tree_node> internal_nodes_;
    std::vector<radix_tree_node> leaf_nodes_;
    TreeResult<ext::extents_t> bvh_;
    std::vector<vec3> data_;
    std::vector<index_t> adjacency_;
    std::vector<MassPoint> coms_;
    std::vector<MortonT> hashes_;

    // Type aliases - now using permuted_simplex_view
    using ptr = std::shared_ptr<bvh_tree<N, MortonT>>;
    using permutation_index_type = std::vector<index_t>;
    using permuted_view = permuted_simplex_view<N, std::vector<vec3>, std::vector<index_t>>;
    using value_type = typename permuted_view::value_type; // std::array<vec3, N>

    // Member variables for views - optional because views are immutable after construction
    std::optional<permuted_view> permuted_data_view_;
    
    static ptr create(const std::vector<index_t> &adjacency,
                      const std::vector<vec3> &vertices, int lvl = 8) {
      return std::make_shared<bvh_tree<N, MortonT>>(vertices, adjacency);
    }

    static ptr create(const simplex_set<N> &set, int lvl = 8) {
      return std::make_shared<bvh_tree<N, MortonT>>(set);
    }

    bvh_tree(const std::vector<vec3> &data,
             const std::vector<index_t> &adjacency) {
      update(data, adjacency);
    }

    bvh_tree(const simplex_set<N> &set) { update(set); }

    index_t get_index(size_t i) const { return indices_[i]; }
    
    void update(const std::vector<vec3> &data,
                const std::vector<index_t> &adjacency) {
      data_ = data;
      adjacency_ = adjacency;
      
      // Create identity permutation for initial hashing
      std::vector<index_t> identity(adjacency_.size() / N);
      for (size_t i = 0; i < identity.size(); ++i) {
        identity[i] = static_cast<index_t>(i);
      }
      
      // Create initial simplex view to compute hash
      permuted_simplex_view<N, std::vector<vec3>, std::vector<index_t>>
          initial_view(data_, adjacency_, identity);
      
      auto [hashes, indices, internal_nodes, leaf_nodes] = make_hash<MortonT>(initial_view);
      
      // Copy to members FIRST so the view can reference stable storage
      indices_ = std::move(indices);
      internal_nodes_ = std::move(internal_nodes);
      leaf_nodes_ = std::move(leaf_nodes);
      hashes_ = std::move(hashes);

      // Now create the permuted view pointing to member storage
      permuted_data_view_.emplace(data_, adjacency_, indices_);

      bvh_ = make_bvh(*permuted_data_view_, internal_nodes_, leaf_nodes_);
      coms_ = calc_com(*permuted_data_view_);
    }

    void update(const simplex_set<N> &set) {
      data_ = set.vertices();
      adjacency_ = set.adjacency();
      
      // Create identity permutation for initial hashing
      std::vector<index_t> identity(adjacency_.size() / N);
      for (size_t i = 0; i < identity.size(); ++i) {
        identity[i] = static_cast<index_t>(i);
      }
      
      // Create initial simplex view to compute hash
      permuted_simplex_view<N, std::vector<vec3>, std::vector<index_t>>
          initial_view(data_, adjacency_, identity);
      
      auto [hashes, indices, internal_nodes, leaf_nodes] = make_hash<MortonT>(initial_view);
      
      // Copy to members FIRST so the view can reference stable storage
      indices_ = std::move(indices);
      internal_nodes_ = std::move(internal_nodes);
      leaf_nodes_ = std::move(leaf_nodes);
      hashes_ = std::move(hashes);

      // Now create the permuted view pointing to member storage
      permuted_data_view_.emplace(data_, adjacency_, indices_);

      bvh_ = make_bvh(*permuted_data_view_, internal_nodes_, leaf_nodes_);
      coms_ = calc_com(*permuted_data_view_);
    }

    std::array<index_t, N> get_tuple_ids(const index_t &i) {
      return permuted_data_view_->get_tuple_ids(i);
    }

    // Get the simplex at index i (returns std::array<vec3, N>)
    value_type get_simplex(index_t i) const {
      return (*permuted_data_view_)[i];
    }

    // Get center of mass for a leaf node
    vec3 get_com(index_t i) const {
      return std::get<1>(coms_[i]);
    }

    template <int Nq>
    auto dispatch_test() const {
      if constexpr (Nq == 1 && N == 1)
        return [](const auto &q, const auto &d) { return test_point_point_tuple(q, d); };
      else if constexpr (Nq == 1 && N == 2)
        return [](const auto &q, const auto &d) { return test_point_line_tuple(q, d); };
      else if constexpr (Nq == 1 && N == 3)
        return [](const auto &q, const auto &d) { return test_point_tri_tuple(q, d); };
      else if constexpr (Nq == 2 && N == 1)
        return [](const auto &q, const auto &d) { return test_point_line_tuple(d, q); };
      else if constexpr (Nq == 2 && N == 2)
        return [](const auto &q, const auto &d) { return test_line_line_tuple(q, d); };
      else if constexpr (Nq == 2 && N == 3)
        return [](const auto &q, const auto &d) { return test_line_tri_tuple(q, d); };
      else if constexpr (Nq == 3 && N == 1)
        return [](const auto &q, const auto &d) { return test_point_tri_tuple(d, q); };
      else if constexpr (Nq == 3 && N == 2)
        return [](const auto &q, const auto &d) { return test_line_tri_tuple(d, q); };
      else if constexpr (Nq == 3 && N == 3)
        return [](const auto &q, const auto &d) { return test_tri_tri_tuple(q, d); };
    }

    template <int Nq>
    static auto make_query(const auto &view) {
      if constexpr (Nq == 1) return std::array<vec3, 1>{view[0]};
      else if constexpr (Nq == 2) return std::array<vec3, 2>{view[0], view[1]};
      else if constexpr (Nq == 3) return std::array<vec3, 3>{view[0], view[1], view[2]};
    }

    template <Vec3View PTYPE>
    index_t find_nearest(const PTYPE &query) {
      constexpr int Nq = view_stride_v<PTYPE>;
      return arp::get_nearest<Singulus<Nq>, permuted_view>(
          make_query<Nq>(query), *permuted_data_view_,
          internal_nodes_, leaf_nodes_, bvh_, dispatch_test<Nq>());
    }

    template <Vec3View PTYPE>
    std::vector<index_t> find_neighbors(const PTYPE &query, real tol) {
      constexpr int Nq = view_stride_v<PTYPE>;
      return arp::get_neighbors<Singulus<Nq>, permuted_view>(
          make_query<Nq>(query), *permuted_data_view_,
          internal_nodes_, leaf_nodes_, bvh_, tol, dispatch_test<Nq>());
    }

    // Legacy wrapper: returns vector with idMin as last element
    template <Vec3View PTYPE>
    std::vector<index_t> get_nearest(const PTYPE &query, real tol) {
      std::vector<index_t> result;
      if (tol > 999.9) {
        result.push_back(find_nearest<PTYPE>(query));
      } else {
        result = find_neighbors<PTYPE>(query, tol);
        result.push_back(-1);
      }
      return result;
    }
  };

template <int N>
using BVH_T = bvh_tree<N, morton_t>;

// Explicit instantiation declarations - controlled by CMake option
#if defined(GAUDI_USE_EXPLICIT_INSTANTIATIONS) &&                              \
    GAUDI_USE_EXPLICIT_INSTANTIATIONS

extern template std::vector<MassPoint> calc_com<1>(const std::vector<vec3> &);
extern template std::vector<MassPoint> calc_com<2>(const std::vector<vec3> &);
extern template std::vector<MassPoint> calc_com<3>(const std::vector<vec3> &);

extern template std::vector<ext::extents_t>
calc_extents<1>(const std::vector<vec3> &);
extern template std::vector<ext::extents_t>
calc_extents<2>(const std::vector<vec3> &);
extern template std::vector<ext::extents_t>
calc_extents<3>(const std::vector<vec3> &);

#endif

} // namespace arp

} // namespace gaudi

#endif // __GAUDI_ARP_HASH_TREE__