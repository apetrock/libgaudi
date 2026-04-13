#ifndef __GAUDI_PYRAMID_DATUM_TESTS_HPP__
#define __GAUDI_PYRAMID_DATUM_TESTS_HPP__

#include "gaudi/arp/hash_tree.hpp"
#include "gaudi/arp/datums.hpp"
#include "gaudi/common.h"
#include "gaudi/test/test.hpp"
#include <cmath>
#include <numeric>
#include <vector>

namespace gaudi {
namespace test {

// Build a small radix tree from a regular 4x4x4 grid (64 points).
// Returns {hashes, indices, internal_nodes, leaf_nodes, sorted_centers}.
struct SmallTreeFixture {
  std::vector<vec3> centers;
  std::vector<arp::morton_t> hashes;
  std::vector<index_t> indices;
  std::vector<arp::radix_tree_node> internal_nodes;
  std::vector<arp::radix_tree_node> leaf_nodes;
  std::vector<vec3> sorted_centers;

  SmallTreeFixture() {
    for (int x = 0; x < 4; x++)
      for (int y = 0; y < 4; y++)
        for (int z = 0; z < 4; z++)
          centers.emplace_back(real(x), real(y), real(z));

    auto [h, idx] = arp::make_hash_3d_t<arp::morton_t>(centers);
    hashes = std::move(h);
    indices = std::move(idx);

    auto [lnodes, inodes] = arp::build_tree(hashes);
    leaf_nodes = std::move(lnodes);
    internal_nodes = std::move(inodes);

    sorted_centers.resize(centers.size());
    for (size_t i = 0; i < indices.size(); i++)
      sorted_centers[i] = centers[indices[i]];
  }
};

GAUDI_TEST(pyramid_additive_scalar_root_is_sum) {
  SmallTreeFixture fix;

  std::vector<real> weights(fix.sorted_centers.size());
  for (size_t i = 0; i < weights.size(); i++)
    weights[i] = real(i + 1);

  real expected_sum = 0.0;
  for (auto w : weights) expected_sum += w;

  auto pyramid = arp::build_pyramid(
      weights, fix.internal_nodes, fix.leaf_nodes,
      [](const real &a, const real &b) -> real { return a + b; },
      real(0));

  GAUDI_ASSERT(pyramid.size() == weights.size() - 1);

  index_t root = 0;
  for (size_t i = 0; i < fix.internal_nodes.size(); i++) {
    if (fix.internal_nodes[i].parent == arp::UNULL) {
      root = static_cast<index_t>(i);
      break;
    }
  }
  GAUDI_EXPECT(std::abs(pyramid[root] - expected_sum) < 1e-10);
}

GAUDI_TEST(pyramid_extent_union_root_is_global_bbox) {
  SmallTreeFixture fix;

  auto pyramid = arp::build_pyramid<vec3, ext::extents_t>(
      fix.sorted_centers, fix.internal_nodes, fix.leaf_nodes,
      [](const vec3 &pt, const ext::extents_t &e) -> ext::extents_t {
        return ext::expand(e, pt);
      },
      [](const ext::extents_t &a, const ext::extents_t &b) -> ext::extents_t {
        return ext::expand(b, a);
      },
      ext::init());

  GAUDI_ASSERT(pyramid.size() == fix.sorted_centers.size() - 1);

  vec3 expected_min(0, 0, 0);
  vec3 expected_max(3, 3, 3);

  index_t root = 0;
  for (size_t i = 0; i < fix.internal_nodes.size(); i++) {
    if (fix.internal_nodes[i].parent == arp::UNULL) {
      root = static_cast<index_t>(i);
      break;
    }
  }

  const auto &root_ext = pyramid[root];
  for (int k = 0; k < 3; k++) {
    GAUDI_EXPECT(std::abs(root_ext[0][k] - expected_min[k]) < 1e-10);
    GAUDI_EXPECT(std::abs(root_ext[1][k] - expected_max[k]) < 1e-10);
  }
}

GAUDI_TEST(pyramid_edge_frame_heterogeneous_root_is_sum_eeT) {
  SmallTreeFixture fix;

  std::vector<vec3> edges(fix.sorted_centers.size());
  mat3 expected = mat3::Zero();
  for (size_t i = 0; i < edges.size(); i++) {
    edges[i] = vec3(real(i % 3), real((i + 1) % 3), real((i + 2) % 3));
    expected += edges[i] * edges[i].transpose();
  }

  auto pyramid = arp::build_pyramid<vec3, mat3>(
      edges, fix.internal_nodes, fix.leaf_nodes,
      [](const vec3 &e, const mat3 &F) -> mat3 {
        return F + e * e.transpose();
      },
      [](const mat3 &a, const mat3 &b) -> mat3 { return a + b; },
      mat3::Zero());

  GAUDI_ASSERT(pyramid.size() == edges.size() - 1);

  index_t root = 0;
  for (size_t i = 0; i < fix.internal_nodes.size(); i++) {
    if (fix.internal_nodes[i].parent == arp::UNULL) {
      root = static_cast<index_t>(i);
      break;
    }
  }

  const mat3 &root_val = pyramid[root];
  for (int r = 0; r < 3; r++)
    for (int c = 0; c < 3; c++)
      GAUDI_EXPECT(std::abs(root_val(r, c) - expected(r, c)) < 1e-10);
}

GAUDI_TEST(pyramid_extents_datum_matches_make_bvh) {
  SmallTreeFixture fix;

  auto bvh_result = arp::make_bvh<1>(fix.sorted_centers,
                                     fix.internal_nodes, fix.leaf_nodes);

  std::vector<index_t> identity_adj(fix.centers.size());
  std::iota(identity_adj.begin(), identity_adj.end(), 0);

  calder::extents_datum ext_datum(identity_adj, fix.centers);
  ext_datum.do_pyramid(fix.indices, fix.internal_nodes, fix.leaf_nodes);

  GAUDI_ASSERT(ext_datum.node_data().size() == bvh_result.internal.size());

  index_t root = 0;
  for (size_t i = 0; i < fix.internal_nodes.size(); i++) {
    if (fix.internal_nodes[i].parent == arp::UNULL) {
      root = static_cast<index_t>(i);
      break;
    }
  }

  const auto &datum_root = ext_datum.node_data()[root];
  const auto &bvh_root = bvh_result.internal[root];
  for (int k = 0; k < 3; k++) {
    GAUDI_EXPECT(std::abs(datum_root[0][k] - bvh_root[0][k]) < 1e-10);
    GAUDI_EXPECT(std::abs(datum_root[1][k] - bvh_root[1][k]) < 1e-10);
  }
}

} // namespace test
} // namespace gaudi

#endif // __GAUDI_PYRAMID_DATUM_TESTS_HPP__
