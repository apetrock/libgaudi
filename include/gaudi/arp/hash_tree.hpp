#ifndef __GAUDI_ARP_HASH_TREE__
#define __GAUDI_ARP_HASH_TREE__

#include "gaudi/arp/morton.hpp"
#include "gaudi/common.h"
#include "gaudi/console_logger.hpp"
#include "gaudi/geometry_logger.hpp"
#include "gaudi/geometry_types.hpp"
#include <algorithm>
#include <functional>
#include <iostream>
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
      console_logger::debug << "Leaf node at " << range.first << " with parent "
                            << i << std::endl;
      leaf_nodes[range.first].parent = i;
    } else {
      console_logger::debug << "Internal node at " << split << " with parent "
                            << i << std::endl;
      internal_nodes[split].parent = i;
    }

    if (split + 1 == range.second) {
      console_logger::debug << "Leaf node at " << range.second
                            << " with parent " << i << std::endl;
      leaf_nodes[range.second].parent = i;
    } else {
      console_logger::debug << "Internal node at " << split + 1
                            << " with parent " << i << std::endl;
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
inline void traverse_bfs(const std::vector<radix_tree_node> &internal_nodes,
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
      console_logger::debug << "Traversing internal node " << cid << std::endl;
      // inspect children and verify parents of children equal to cnode
      const auto &node_split = internal_nodes[cnode.split];
      console_logger::debug << "cnode.split=" << cnode.split
                            << ", in[split]={start:" << node_split.start
                            << ", end:" << node_split.end
                            << ", split:" << node_split.split
                            << ", parent:" << node_split.parent
                            << "}, parent=" << node_split.parent << std::endl;
      const auto &node_split1 = internal_nodes[cnode.split + 1];
      console_logger::debug << "cnode.split+1=" << cnode.split + 1
                            << ", in[split+1]={start:" << node_split1.start
                            << ", end:" << node_split1.end
                            << ", split:" << node_split1.split
                            << ", parent:" << node_split1.parent
                            << "}, parent=" << node_split1.parent << std::endl;
      stack.push(cnode.split);
      stack.push(cnode.split + 1);
    }
    // If fcn returns false, we stop propagating down this branch
  }
}

// Build pyramid using map-reduce pattern
// I input type, O output type
// takes an input vector of data, maps it to an output type using a map
// function, then walks the data up the tree and reduces it using a reduce
// function
template <typename O>
inline std::vector<O> build_pyramid(
    const std::vector<O> data, const std::vector<index_t> &data_indices,
    const std::vector<radix_tree_node> &internal_nodes,
    const std::vector<radix_tree_node> &leaf_nodes,
    std::function<O(const O &, const O &)> reduce_func, const O &default_val) {

  auto thread_safe_reduce = [&](const O &a, const O &b) {
    // lock here
    return reduce_func(a, b);
    // unlock here
  };
  std::vector<O> internalReduce(data.size() - 1, default_val);
  for (int i = 0; i < data.size(); i++) {
    const O &datai = data[data_indices[i]];
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
inline std::tuple<std::vector<uint32_t>, std::vector<index_t>,
                  std::vector<radix_tree_node>, std::vector<radix_tree_node>>
make_hash_tree(const std::vector<vec3> &data) {
  auto [hashes, indices] = make_hash_3d(data);
  if (hashes.empty()) {
    return {{}, {}, {}, {}};
  }

  auto [leaf_nodes, internal_nodes] = build_tree(hashes);
  return std::make_tuple(hashes, indices, internal_nodes, leaf_nodes);
}

template <int N>
inline TreeResult<ext::extents_t>
make_bvh(const std::vector<vec3> &data, const std::vector<index_t> &indices,
         const std::vector<uint32_t> &hashes,
         const std::vector<radix_tree_node> &internal_nodes,
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

  console_logger::debug << "Building BVH with " << indices.size() << " nodes"
                        << std::endl;
  auto internalReduce = build_pyramid<ext::extents_t>(
      exts, indices, internal_nodes, leaf_nodes, reduce_function, default_val);
  console_logger::debug << "Built BVH with " << exts.size() << " leaf nodes, "
                        << internalReduce.size() << " internal nodes"
                        << std::endl;

  return {
      exts,
      internalReduce,
  };
}

template <int N>
inline TreeResult<vec3>
make_points(const std::vector<vec3> &data, const std::vector<index_t> &indices,
            const std::vector<uint32_t> &hashes,
            const std::vector<radix_tree_node> &internal_nodes,
            const std::vector<radix_tree_node> &leaf_nodes) {
  using ITYPE = vec3;
  using OTYPE = MassPoint; // accumulate type
  const auto default_val = OTYPE{0.0f, vec3::Zero()};
  const auto coms = calc_com<N>(data); // Specify template parameter N=1

  auto reduce_function = [&](const OTYPE &a, const OTYPE &b) -> OTYPE {
    const real &ma = std::get<0>(a);
    const real &mb = std::get<0>(b);
    const vec3 &pa = ma * std::get<1>(a);
    const vec3 &pb = std::get<1>(b); // accumulate point sum
    return OTYPE{ma + mb, pa + pb};
  };

  // Build pyramid with map-reduce
  auto internalReduce = build_pyramid<OTYPE>(
      coms, indices, internal_nodes, leaf_nodes, reduce_function, default_val);

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

template <int N = 1>
void log_hierarchy(const std::vector<vec3> &data,
                   const std::vector<index_t> &indices,
                   const std::vector<uint32_t> &hashes,
                   const std::vector<radix_tree_node> &internal_nodes,
                   const std::vector<radix_tree_node> &leaf_nodes) {
  if (internal_nodes.empty() || leaf_nodes.empty()) {
    console_logger::debug << "No internal or leaf nodes to log" << std::endl;
    return;
  }
  // Compute centerpoints for each stride group for visualization
  const auto [leaf_points, internal_points] =
      make_points<N>(data, indices, hashes, internal_nodes, leaf_nodes);
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
          const auto &cpoint = leaf_points[indices[leaf_id]];
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

template <int N>
void log_bvh(const std::vector<vec3> &data, const std::vector<index_t> &indices,
             const std::vector<uint32_t> &hashes,
             const std::vector<radix_tree_node> &internal_nodes,
             const std::vector<radix_tree_node> &leaf_nodes) {

  const auto bvh_result =
      make_bvh<N>(data, indices, hashes, internal_nodes, leaf_nodes);
  const auto &leaf_bvh = bvh_result.leaf;
  const auto &internal_bvh = bvh_result.internal;
  for (int i = 0; i < data.size(); i++) {
    const vec3 &point = data[i];
    geometry_logger::point(point, vec4(1.0, 0.0, 0.0, 1.0));
  }
  for (int i = 0; i < internal_bvh.size(); i++) {
    const ext::extents_t &ext = internal_bvh[i];
    console_logger::debug << "Internal BVH extents: " << ext[0].transpose()
                          << ", " << ext[1].transpose() << std::endl;
    geometry_logger::ext(ext[0], ext[1], vec4(0.0, 1.0, 0.0, 0.5));
  }
}

/*
#if 1
template <int ST, int SS> // T=test, S=set... DOH! T could equal tree...
std::vector<index_t>
getNearest(index_t &idT, const std::vector<index_t> &t_inds,
           const vector<vec3> &t_verts, //
           const std::vector<vec3> &data,
            const std::vector<index_t> &indices,
            const std::vector<uint32_t> &hashes,
            const std::vector<radix_tree_node> &internal_nodes,
            const std::vector<radix_tree_node> &leaf_nodes
              const aabb_tree<SS> &s_tree, real tol,
           std::function<real(const index_t &idT, //
                              const std::vector<index_t> &t_inds,
                              const vector<vec3> &t_verts, //
                              const index_t &idS,          //
                              const std::vector<index_t> &s_inds,
                              const std::vector<vec3> &s_verts)>
               testAB)
{

  // TIMER function//TIMER(__FUNCTION__);
  typedef aabb_tree<SS> tree_type;
  typedef typename tree_type::node Node;

  bool expanding_rad = tol > 999.9;

  index_t idMin = -1;
  real dmin = std::numeric_limits<real>::max();

  const Node &root = s_tree.nodes[0];
  std::stack<int> cstack;
  cstack.push(0);
  bool hit = false;
  // T tol = 0.05;
  ext::extents_t extT = calc_extents<ST>(idT, t_inds, t_verts);
  extT = ext::inflate(extT, tol);
  std::vector<index_t> collisions;
  while (cstack.size() > 0)
  {
    int cId = cstack.top();
    cstack.pop();
    const Node &cnode = s_tree.nodes[cId];

    if (expanding_rad && cnode.size > 0)
    {
      real d = ext::distance(extT, cnode.half.cen);
      tol = std::min(tol, d);
      tol = std::max(tol, cnode.half.mag);
      extT = calc_extents<ST>(idT, t_inds, t_verts);
      extT = ext::inflate(extT, tol);
#if 0
  if (idT == 17260) {

    std::cout << "d/tol: " << d << " " << tol << std::endl;
    vec3 cT = 0.5 * (extT[0] + extT[1]);
    geometry_logger::line(cT, cnode.half.cen, vec4(0.0, 1.0, 0.0, 1.0));
    geometry_logger::ext(extT[0], extT[1], vec4(1.0, 0.0, 0.0, 1.0));
  }
#endif
    }

    if (cnode.children[0] == -1 && cnode.children[1] == -1)
    {

      for (int k = cnode.begin; k < cnode.begin + cnode.size; k++)
      {

        const index_t &idS = s_tree.permutation[k];

        ext::extents_t extS =
            calc_extents<SS>(idS, s_tree.indices(), s_tree.verts());

        if (!ext::overlap(extT, extS))
        {

          continue;
        }
#if 0
    if (idT == 15978 && 0) {
      std::cout << "idS: " << idS << std::endl;
      vec3 cT = 0.5 * (extT[0] + extT[1]);
      vec3 cS = 0.5 * (extS[0] + extS[1]);
      geometry_logger::ext(extS[0], extS[1], vec4(0.0, 1.0, 0.0, 1.0));
      geometry_logger::ext(extT[0], extT[1], vec4(1.0, 0.0, 0.0, 1.0));
      geometry_logger::line(cT, cS, vec4(1.0, 1.0, 0.0, 1.0));
    }
#endif
        real dist = testAB(idT, t_inds, t_verts, //
                           idS, s_tree.indices(), s_tree.verts());

        if (dist < dmin)
        {
          dmin = dist;
          idMin = idS;
        }
        if (dist < tol && !expanding_rad)
        {
          collisions.push_back(idS);
        }
      }
    }

    index_t itx = cnode.half.intersect(extT);

    if (itx <= 0 && cnode.children[0] > 0)
    {
      cstack.push(cnode.children[0]);
    }

    if (itx >= 0 && cnode.children[1] > 0)
    {
      cstack.push(cnode.children[1]);
    }
  }

  if (dmin < tol)
  {
    collisions.push_back(idMin);
  }
  else
  {
    collisions.push_back(-1);
#if 0
if (expanding_rad) {
  std::cout << "no collisions:" << idT << std::endl;
  // exit(0);
}
#endif
  }

  return collisions;
};
#endif
*/
template <int N>
inline std::tuple<std::vector<uint32_t>, std::vector<index_t>,
                  std::vector<radix_tree_node>, std::vector<radix_tree_node>>
make_hash_N(const std::vector<vec3> &data) {
  const auto mass_points = calc_com<N>(data);
  // Extract just the center of mass vectors from MassPoint tuples
  std::vector<vec3> averaged;
  averaged.reserve(mass_points.size());
  for (const auto &[mass, com] : mass_points) {
    averaged.push_back(com);
  }
  return make_hash_tree(averaged);
}

} // namespace arp

} // namespace gaudi

#endif // __GAUDI_ARP_HASH_TREE__