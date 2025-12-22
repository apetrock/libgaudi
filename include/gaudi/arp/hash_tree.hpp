#ifndef __GAUDI_ARP_HASH_TREE__
#define __GAUDI_ARP_HASH_TREE__

#include "gaudi/arp/morton.hpp"
#include "gaudi/arp/pairwise_tests.hpp"
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

// Common leading zero calculation for tree building
inline int clz(index_t i, index_t j, const std::vector<uint32_t> &hash) {
  if (j < 0)
    return -1;
  if (j > static_cast<index_t>(hash.size()) - 1)
    return -1;
  uint32_t code_i = hash[i];
  uint32_t code_j = hash[j];

// Use __builtin_clz for GCC/Clang, or implement fallback
#if defined(__GNUC__) || defined(__clang__)
  return __builtin_clz(code_i ^ code_j);
#else
  // Fallback implementation
  uint32_t diff = code_i ^ code_j;
  if (diff == 0)
    return 32;
  int leading_zeros = 0;
  while ((diff & 0x80000000) == 0) {
    diff <<= 1;
    leading_zeros++;
  }
  return leading_zeros;
#endif
}

// Find split point in radix tree
inline index_t find_split(index_t start, index_t end,
                          const std::vector<uint32_t> &hash) {
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
inline std::pair<index_t, index_t>
find_range(index_t i, const std::vector<uint32_t> &hash) {
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
inline NodeResult build_tree(const std::vector<uint32_t> &hash) {
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
inline void test_tree(const std::vector<radix_tree_node> &internal_nodes,
                      const std::vector<radix_tree_node> &leaf_nodes,
                      const std::vector<index_t> &ids,
                      const std::vector<uint32_t> &hash) {
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
  std::vector<uint32_t> hash = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
  std::vector<index_t> ids = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
  auto [leaf_nodes, internal_nodes] = build_tree(hash);
  test_tree(internal_nodes, leaf_nodes, ids, hash);
}

template <typename F>
inline void traverse_dfs(const std::vector<radix_tree_node> &internal_nodes,
                         const std::vector<radix_tree_node> &leaf_nodes,
                         F &&fcn) {
  if (internal_nodes.empty())
    return;

  std::stack<index_t> stack;
  stack.push(0);

  while (!stack.empty()) {
    index_t cid = stack.top();
    stack.pop();

    const auto &cnode = internal_nodes[cid];
    if (cnode.split == UNULL)
      continue;

    // Handle leaf cases f
    else if (cnode.split + 0 == cnode.start || cnode.split + 1 == cnode.end) {
      // process node that has leaves, but process as node
      fcn(cid, -1, cnode);
      if (cnode.split + 0 == cnode.start && cnode.split + 1 == cnode.end) {
        // Both children are leaves
        fcn(cid, cnode.start, leaf_nodes[cnode.start]);
        fcn(cid, cnode.end, leaf_nodes[cnode.end]);
      } else if (cnode.split + 0 == cnode.start) {
        // Left child is leaf
        fcn(cid, cnode.start, leaf_nodes[cnode.start]);
        stack.push(cnode.split + 1);
      } else if (cnode.split + 1 == cnode.end) {
        // Right child is leaf
        fcn(cid, cnode.end, leaf_nodes[cnode.end]);
        stack.push(cnode.split);
      }
    } else if (fcn(cid, -1, cnode)) { // <- KEY: conditional traversal
      // Continue traversing if callback returns true
      const auto &node_split = internal_nodes[cnode.split];
      const auto &node_split1 = internal_nodes[cnode.split + 1];
      stack.push(cnode.split);
      stack.push(cnode.split + 1);
    }
    // If fcn returns false, we stop propagating down this branch
  }
}

template <typename F>
inline void traverse_bfs(const std::vector<radix_tree_node> &internal_nodes,
                         const std::vector<radix_tree_node> &leaf_nodes,
                         F &&fcn) {
  if (internal_nodes.empty())
    return;

  std::queue<index_t> queue;
  queue.push(0);

  while (!queue.empty()) {
    index_t cid = queue.front();
    queue.pop();

    const auto &cnode = internal_nodes[cid];
    if (cnode.split == UNULL)
      continue;

    // Handle leaf cases
    else if (cnode.split + 0 == cnode.start || cnode.split + 1 == cnode.end) {
      // process node that has leaves, but process as node
      fcn(cid, -1, cnode);
      if (cnode.split + 0 == cnode.start && cnode.split + 1 == cnode.end) {
        // Both children are leaves
        fcn(cid, cnode.start, leaf_nodes[cnode.start]);
        fcn(cid, cnode.end, leaf_nodes[cnode.end]);
      } else if (cnode.split + 0 == cnode.start) {
        // Left child is leaf
        fcn(cid, cnode.start, leaf_nodes[cnode.start]);
        queue.push(cnode.split + 1);
      } else if (cnode.split + 1 == cnode.end) {
        // Right child is leaf
        fcn(cid, cnode.end, leaf_nodes[cnode.end]);
        queue.push(cnode.split);
      }
    } else if (fcn(cid, -1, cnode)) { // <- KEY: conditional traversal
      // Continue traversing if callback returns true
      const auto &node_split = internal_nodes[cnode.split];
      const auto &node_split1 = internal_nodes[cnode.split + 1];
      queue.push(cnode.split);
      queue.push(cnode.split + 1);
    }
    // If fcn returns false, we stop propagating down this branch
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
    while (j < 32 && parent != UNULL) {
      const O &dataj = internalReduce[parent];
      internalReduce[parent] = thread_safe_reduce(datai, dataj);
      parent = internal_nodes[parent].parent;
      j++;
    }
  }

  return internalReduce;
}

// Build complete hash tree with pyramid
template <Vec3View TTYPE>
inline std::tuple<std::vector<uint32_t>, std::vector<index_t>,
                  std::vector<radix_tree_node>, std::vector<radix_tree_node>>
make_hash_tree(const TTYPE &data) {
  auto [hashes, indices] = make_hash_3d(data);
  if (hashes.empty()) {
    return {{}, {}, {}, {}};
  }

  auto [leaf_nodes, internal_nodes] = build_tree(hashes);
  return std::make_tuple(hashes, indices, internal_nodes, leaf_nodes);
}

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

template <int N, Vec3View TTYPE>
void log_hierarchy(const TTYPE &data,
                   const std::vector<radix_tree_node> &internal_nodes,
                   const std::vector<radix_tree_node> &leaf_nodes) {
  if (internal_nodes.empty() || leaf_nodes.empty()) {
    console_logger::debug << "No internal or leaf nodes to log" << std::endl;
    return;
  }
  // Compute centerpoints for each stride group for visualization
  const auto [leaf_points, internal_points] =
      make_points<N>(data, internal_nodes, leaf_nodes);
  // draw the simplex for logging
  if (N > 1) {
    for (size_t i = 0; i < data.size(); i += N) {
      for (int j = 0; j < N; j++) {
        const int j0 = i + j;
        const int j1 = i + (j + 1) % N;
        geometry_logger::line(data[j0], data[j1], vec4(0.2, 0.2, 0.2, 0.5f));
      }
    }
  }
  // draw centerpoints
  for (int i = 0; i < leaf_points.size(); i++) {
    geometry_logger::point(leaf_points[i], vec4(0.0, 1.0, 0.0, 1.0f));
  }
  // Continue with original logic using make_points
  traverse_bfs(
      internal_nodes, leaf_nodes,
      [&](index_t node_id, index_t leaf_id, const radix_tree_node &node) {
        if (leaf_id != -1 && node.parent != UNULL) {
          const auto &ppoint = internal_points[node.parent];
          const auto &cpoint = leaf_points[leaf_id];
          geometry_logger::line(ppoint, cpoint, vec4(0.0, 0.0, 1.0, 1.0f));
          geometry_logger::point(ppoint, vec4(1.0, 0.0, 1.0, 1.0f));
        } else if (node.parent != UNULL) {
          const auto &ppoint = internal_points[node.parent];
          const auto &cpoint = internal_points[node_id];
          geometry_logger::line(ppoint, cpoint, vec4(0.0, 1.0, 0.0, 1.0f));
        }
        return true;
      });
}

template <int N, Vec3View TTYPE>
void log_bvh(const TTYPE &data,
             const std::vector<radix_tree_node> &internal_nodes,
             const std::vector<radix_tree_node> &leaf_nodes) {

  const auto bvh_result = make_bvh<N>(data, internal_nodes, leaf_nodes);
  const auto &leaf_bvh = bvh_result.leaf;
  const auto &internal_bvh = bvh_result.internal;
  for (int i = 0; i < data.size(); i++) {
    const vec3 &point = data[i];
    geometry_logger::point(point, vec4(1.0, 0.0, 0.0, 1.0));
  }
  for (int i = 0; i < internal_bvh.size(); i++) {
    const ext::extents_t &ext = internal_bvh[i];
    geometry_logger::ext(ext[0], ext[1], vec4(0.0, 1.0, 0.0, 0.5));
  }
}

// SLICE = permuted_slice<vec3, NT> || slice<vec3, NT> || array<vec3, NT>

// user should'nt know that this is a slice its just an array
template <int N, Vec3View PTYPE,
          Vec3View TTYPE> // T=test, S=set... DOH! T could equal tree...
std::vector<index_t>
getNearest(const PTYPE &prim, const TTYPE &data,
           const std::vector<radix_tree_node> &internal_nodes,
           const std::vector<radix_tree_node> &leaf_nodes,
           const TreeResult<ext::extents_t> &bvh_result, real tol,
           auto &&testAB) {

  bool contracting_rad = tol > 999.9;
  ext::extents_t ext_t = ext::calc_extents(prim);
  ext_t = ext::inflate(ext_t, tol);
  const vec3 cen_t = ext::center(ext_t);
  index_t idMin = -1;
  real mMin = std::numeric_limits<real>::max();
  std::vector<index_t> collisions;

  using T = typename TTYPE::value_type;
  traverse_bfs(
      internal_nodes, leaf_nodes,
      [&](index_t node_id, index_t leaf_id, const radix_tree_node &node) {
        if (leaf_id != -1 && node.parent != UNULL) {
          // if contracting_rad
          const slice<N, TTYPE> datum(data, leaf_id);
          const ext::extents_t &ext_s = ext::calc_extents(datum);
          geometry_logger::ext(ext_s[0], ext_s[1], vec4(0.0, 1.0, 0.0, 0.5));
          // geometry_logger::ext(ext_t[0], ext_t[1], vec4(1.0, 0.0, 0.0, 0.5));
          if (contracting_rad) {
            real dist = testAB(prim, datum);
            if (dist < mMin) {
              mMin = dist;
              idMin = data.get_index(leaf_id);
            }
            ext_t = ext::inflate(ext::calc_extents(prim), mMin);
          } else if (ext::overlap(ext_t, ext_s)) {
            real dist = testAB(prim, datum);
            if (dist < tol) {
              collisions.push_back(data.get_index(leaf_id));
            }
          }

          return false;
        } else if (node.parent != UNULL) {
          const ext::extents_t &ext_s = bvh_result.internal[node_id];
          if (contracting_rad) {
            const real dist = ext::dist_from_center(ext_s, cen_t);
            tol = std::min(dist, tol);
            ext_t = ext::inflate(ext::calc_extents(prim), tol);
            return ext::overlap(ext_t, ext_s);
          }

          geometry_logger::ext(ext_t[0], ext_t[1], vec4(1.0, 0.0, 0.0, 0.5));
          // reset the extents to the original
          return ext::overlap(ext_t, ext_s);
        }
        return true;
      });

  if (contracting_rad) {
    collisions.push_back(idMin); // mintol always in the back
  } else {
    collisions.push_back(-1);
  }
  return collisions;
};

template <int N, Vec3View TTYPE>
inline std::tuple<std::vector<uint32_t>, std::vector<index_t>,
                  std::vector<radix_tree_node>, std::vector<radix_tree_node>>
make_hash(const TTYPE &data) {
  const auto mass_points = calc_com<N>(data);
  // Extract just the center of mass vectors from MassPoint tuples
  std::vector<vec3> averaged;
  averaged.reserve(mass_points.size());
  for (const auto &[mass, com] : mass_points) {
    averaged.push_back(com);
  }
  return make_hash_tree(averaged);
}


template <int N>
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
    std::vector<uint32_t> hashes_;

    // Type aliases using the member variables
    using ptr = std::shared_ptr<bvh_tree<N>>;
    using view = adjacency_view<std::vector<vec3>, std::vector<index_t>>;
    using permutation_index_type = std::vector<index_t>;
    using permuted_view = permuted_adjacency_view<N, std::vector<vec3>, std::vector<index_t>>;

    // Member variables for views - optional because views are immutable after construction
    std::optional<view> data_view_;
    std::optional<permuted_view> permuted_data_view_;
    
    static ptr create(const std::vector<index_t> &adjacency,
                      const std::vector<vec3> &vertices, int lvl = 8) {
      return std::make_shared<bvh_tree<N>>(vertices, adjacency);
    }

    bvh_tree(const std::vector<vec3> &data,
             const std::vector<index_t> &adjacency) {
      update(data, adjacency);
    }

    index_t get_index(size_t i) const { return indices_[i]; }
    
    void update(const std::vector<vec3> &data,
                const std::vector<index_t> &adjacency) {
      data_ = data;
      adjacency_ = adjacency;
      data_view_.emplace(data_, adjacency_);
      auto [hashes, indices, internal_nodes, leaf_nodes] =
          make_hash<N>(*data_view_);
      permuted_data_view_.emplace(data_, adjacency_, indices);

      indices_ = indices;
      internal_nodes_ = internal_nodes;
      leaf_nodes_ = leaf_nodes;
      bvh_ = make_bvh<N>(*permuted_data_view_, internal_nodes_, leaf_nodes_);
      coms_ = calc_com<N>(*permuted_data_view_);
      hashes_ = hashes;
    }

    std::array<index_t, N> get_tuple_ids(const index_t & i){
      return permuted_data_view_->get_tuple_ids(i);
    }

    // Get center of mass for a leaf node
    vec3 get_com(index_t i) const {
      return std::get<1>(coms_[i]);
    }

    // Templated get_nearest dispatches to appropriate test function
    // based on query stride (deduced from PTYPE) and N (data primitive stride)
    template <Vec3View PTYPE>
    std::vector<index_t> get_nearest(const PTYPE &query, real tol) {
      constexpr int Nquery = view_stride_v<PTYPE>;
      if constexpr (Nquery == 1 && N == 1) {
        // point query against point data
        return getNearest<N, PTYPE, permuted_view>(
            query, *permuted_data_view_, internal_nodes_, leaf_nodes_, bvh_, tol,
            [](const PTYPE &q, const slice<N, permuted_view> &d) {
              return test_point_point(q, d);
            });
      } else if constexpr (Nquery == 1 && N == 2) {
        // point query against line data
        return getNearest<N, PTYPE, permuted_view>(
            query, *permuted_data_view_, internal_nodes_, leaf_nodes_, bvh_, tol,
            [](const PTYPE &q, const slice<N, permuted_view> &d) {
              return test_point_line(q, d);
            });
      } else if constexpr (Nquery == 1 && N == 3) {
        // point query against triangle data
        return getNearest<N, PTYPE, permuted_view>(
            query, *permuted_data_view_, internal_nodes_, leaf_nodes_, bvh_, tol,
            [](const PTYPE &q, const slice<N, permuted_view> &d) {
              return test_point_tri(q, d);
            });
      } else if constexpr (Nquery == 2 && N == 1) {
        // line query against point data - use point-line with reversed args
        return getNearest<N, PTYPE, permuted_view>(
            query, *permuted_data_view_, internal_nodes_, leaf_nodes_, bvh_, tol,
            [](const PTYPE &q, const slice<N, permuted_view> &d) {
              return test_point_line(d, q);
            });
      } else if constexpr (Nquery == 2 && N == 2) {
        // line query against line data
        return getNearest<N, PTYPE, permuted_view>(
            query, *permuted_data_view_, internal_nodes_, leaf_nodes_, bvh_, tol,
            [](const PTYPE &q, const slice<N, permuted_view> &d) {
              return test_line_line(q, d);
            });
      } else if constexpr (Nquery == 2 && N == 3) {
        // line query against triangle data
        return getNearest<N, PTYPE, permuted_view>(
            query, *permuted_data_view_, internal_nodes_, leaf_nodes_, bvh_, tol,
            [](const PTYPE &q, const slice<N, permuted_view> &d) {
              return test_line_tri(q, d);
            });
      } else if constexpr (Nquery == 3 && N == 1) {
        // triangle query against point data - use point-tri with reversed args
        return getNearest<N, PTYPE, permuted_view>(
            query, *permuted_data_view_, internal_nodes_, leaf_nodes_, bvh_, tol,
            [](const PTYPE &q, const slice<N, permuted_view> &d) {
              return test_point_tri(d, q);
            });
      } else if constexpr (Nquery == 3 && N == 2) {
        // triangle query against line data - use line-tri with reversed args
        return getNearest<N, PTYPE, permuted_view>(
            query, *permuted_data_view_, internal_nodes_, leaf_nodes_, bvh_, tol,
            [](const PTYPE &q, const slice<N, permuted_view> &d) {
              return test_line_tri(d, q);
            });
      } else if constexpr (Nquery == 3 && N == 3) {
        // triangle query against triangle data
        return getNearest<N, PTYPE, permuted_view>(
            query, *permuted_data_view_, internal_nodes_, leaf_nodes_, bvh_, tol,
            [](const PTYPE &q, const slice<N, permuted_view> &d) {
              return test_tri_tri(q, d);
            });
      }
      return {};
    }
  };

// Explicit instantiation declarations - controlled by CMake option
#if defined(GAUDI_USE_EXPLICIT_INSTANTIATIONS) &&                              \
    GAUDI_USE_EXPLICIT_INSTANTIATIONS
extern template TreeResult<vec3>
make_points<1>(const std::vector<vec3> &, const std::vector<index_t> &,
               const std::vector<uint32_t> &,
               const std::vector<radix_tree_node> &,
               const std::vector<radix_tree_node> &);
extern template TreeResult<vec3>
make_points<2>(const std::vector<vec3> &, const std::vector<index_t> &,
               const std::vector<uint32_t> &,
               const std::vector<radix_tree_node> &,
               const std::vector<radix_tree_node> &);
extern template TreeResult<vec3>
make_points<3>(const std::vector<vec3> &, const std::vector<index_t> &,
               const std::vector<uint32_t> &,
               const std::vector<radix_tree_node> &,
               const std::vector<radix_tree_node> &);

extern template TreeResult<ext::extents_t>
make_bvh<1>(const std::vector<vec3> &, const std::vector<index_t> &,
            const std::vector<uint32_t> &, const std::vector<radix_tree_node> &,
            const std::vector<radix_tree_node> &);
extern template TreeResult<ext::extents_t>
make_bvh<2>(const std::vector<vec3> &, const std::vector<index_t> &,
            const std::vector<uint32_t> &, const std::vector<radix_tree_node> &,
            const std::vector<radix_tree_node> &);
extern template TreeResult<ext::extents_t>
make_bvh<3>(const std::vector<vec3> &, const std::vector<index_t> &,
            const std::vector<uint32_t> &, const std::vector<radix_tree_node> &,
            const std::vector<radix_tree_node> &);

extern template std::vector<MassPoint> calc_com<1>(const std::vector<vec3> &);
extern template std::vector<MassPoint> calc_com<2>(const std::vector<vec3> &);
extern template std::vector<MassPoint> calc_com<3>(const std::vector<vec3> &);

extern template std::vector<ext::extents_t>
calc_extents<1>(const std::vector<vec3> &);
extern template std::vector<ext::extents_t>
calc_extents<2>(const std::vector<vec3> &);
extern template std::vector<ext::extents_t>
calc_extents<3>(const std::vector<vec3> &);

extern template std::vector<index_t> getNearest<2>(
    const std::array<vec3, 1> &, const std::vector<vec3> &,
    const std::vector<index_t> &, const std::vector<radix_tree_node> &,
    const std::vector<radix_tree_node> &, const TreeResult<ext::extents_t> &,
    real,
    std::function<real(const std::vector<vec3> &, const near_array<1> &)>);

extern template std::vector<index_t> getNearest<2>(
    const std::array<vec3, 2> &, const std::vector<vec3> &,
    const std::vector<index_t> &, const std::vector<radix_tree_node> &,
    const std::vector<radix_tree_node> &, const TreeResult<ext::extents_t> &,
    real,
    std::function<real(const std::vector<vec3> &, const near_array<1> &)>);

#endif

} // namespace arp

} // namespace gaudi

#endif // __GAUDI_ARP_HASH_TREE__