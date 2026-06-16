#ifndef __GAUDI_BVH_TESTS_HPP__
#define __GAUDI_BVH_TESTS_HPP__

#include "gaudi/arp/brute_force.hpp"
#include "gaudi/arp/hash_tree.hpp"
#include "gaudi/arp/simplex_set.hpp"
#include "gaudi/asawa/faceloader.hpp"
#include "gaudi/asawa/primitive_objects.hpp"
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/common.h"
#include "gaudi/paths.hpp"
#include "gaudi/test/test.hpp"
#include "gaudi/vec_addendum.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <random>
#include <set>
#include <string>
#include <vector>

namespace gaudi {
namespace test {

struct MeshData {
  std::vector<vec3> vertices;
  std::vector<std::vector<int>> faces;
  std::vector<index_t> point_adj;
  std::vector<index_t> edge_adj;
  std::vector<index_t> tri_adj;
  asawa::shell::shell::ptr shell;
};

// Procedural sphere mesh -- self-contained, no asset file dependency.
inline MeshData load_sphere_mesh() {
  MeshData mesh;
  asawa::make_sphere(mesh.vertices, mesh.faces);
  GAUDI_ASSERT(!mesh.vertices.empty());
  GAUDI_ASSERT(!mesh.faces.empty());

  std::vector<index_t> corners_next;
  std::vector<index_t> corners_vert;
  std::vector<index_t> corners_face;
  asawa::assemble_table(mesh.vertices, mesh.faces, corners_next, corners_vert,
                        corners_face);
  mesh.shell =
      asawa::shell::shell::create(corners_next, corners_vert, corners_face);
  GAUDI_ASSERT(mesh.shell != nullptr);

  {
    auto vr = mesh.shell->get_vert_range();
    mesh.point_adj.assign(vr.begin(), vr.end());
  }
  mesh.edge_adj = mesh.shell->get_edge_vert_ids();
  mesh.tri_adj = mesh.shell->get_face_vert_ids(true);
  GAUDI_ASSERT(!mesh.point_adj.empty());
  GAUDI_ASSERT(!mesh.edge_adj.empty());
  GAUDI_ASSERT(!mesh.tri_adj.empty());
  GAUDI_ASSERT(mesh.tri_adj.size() % 3 == 0);
  return mesh;
}

inline std::vector<vec3> random_points_in_bounds(
    const std::vector<vec3> &vertices, size_t count) {
  vec3 min_val = vertices.front();
  vec3 max_val = vertices.front();
  for (const auto &v : vertices) {
    min_val = va::min(min_val, v);
    max_val = va::max(max_val, v);
  }

  std::mt19937 rng(1337);
  std::uniform_real_distribution<real> dx(min_val[0], max_val[0]);
  std::uniform_real_distribution<real> dy(min_val[1], max_val[1]);
  std::uniform_real_distribution<real> dz(min_val[2], max_val[2]);

  std::vector<vec3> points;
  points.reserve(count);
  for (size_t i = 0; i < count; ++i) {
    points.emplace_back(dx(rng), dy(rng), dz(rng));
  }
  return points;
}

GAUDI_TEST(bvh_small_grid_traverse) {
  // 4x4x4 = 64 points -- small enough to trace every iteration
  std::vector<vec3> verts;
  std::vector<index_t> adj;
  for (int x = 0; x < 4; x++)
    for (int y = 0; y < 4; y++)
      for (int z = 0; z < 4; z++) {
        adj.push_back(static_cast<index_t>(verts.size()));
        verts.emplace_back(real(x), real(y), real(z));
      }

  arp::simplex_set<1> pts(verts, adj);
  auto bvh = arp::bvh_tree<1>::create(pts);

  // BFS traverse with should_continue=true -- must visit all nodes/leaves
  size_t node_count = 0, leaf_count = 0;
  arp::traverse_bfs(
      bvh->internal_nodes_, bvh->leaf_nodes_,
      [&](index_t, const arp::radix_tree_node &) { node_count++; },
      [&](index_t, index_t, const arp::radix_tree_node &) { leaf_count++; },
      [&](index_t, const arp::radix_tree_node &) { return true; });

  GAUDI_ASSERT(node_count == bvh->internal_nodes_.size());
  GAUDI_ASSERT(leaf_count == bvh->leaf_nodes_.size());
}

GAUDI_TEST(bvh_point_to_point) {
  auto mesh = load_sphere_mesh();
  arp::simplex_set<1> points(mesh.vertices, mesh.point_adj);
  auto bvh = arp::bvh_tree<1>::create(points);

  adjacency_view<std::vector<vec3>, std::vector<index_t>> point_view(
      mesh.vertices, mesh.point_adj);

  const real tol = 1000.0;
  auto queries = random_points_in_bounds(mesh.vertices, 10);
  for (const auto &q : queries) {
    std::array<vec3, 1> query = {q};
    auto brute_res =
        arp::brute_force_nearest_auto<1, 1>(query, point_view, tol);
    auto bvh_res = bvh->get_nearest(query, tol);

    index_t bvh_id = bvh_res.empty() ? -1 : bvh_res.back();
    index_t brute_id = brute_res.empty() ? -1 : brute_res.back();
    GAUDI_ASSERT(bvh_id >= 0);
    GAUDI_ASSERT(brute_id >= 0);

    const auto bvh_dist = arp::test_point_point(query, slice<1, decltype(point_view)>(point_view, bvh_id));
    const auto brute_dist = arp::test_point_point(query, slice<1, decltype(point_view)>(point_view, brute_id));
    GAUDI_EXPECT(std::abs(bvh_dist - brute_dist) < 1e-6);
  }
}

GAUDI_TEST(bvh_point_to_edge) {
  auto mesh = load_sphere_mesh();
  arp::simplex_set<2> edges(mesh.vertices, mesh.edge_adj);
  auto bvh = arp::bvh_tree<2>::create(edges);

  adjacency_view<std::vector<vec3>, std::vector<index_t>> edge_view(
      mesh.vertices, mesh.edge_adj);

  const real tol = 1000.0;
  auto queries = random_points_in_bounds(mesh.vertices, 10);
  for (const auto &q : queries) {
    std::array<vec3, 1> query = {q};
    auto bvh_res = bvh->get_nearest(query, tol);
    auto brute_res =
        arp::brute_force_nearest_auto<1, 2>(query, edge_view, tol);

    index_t bvh_id = bvh_res.empty() ? -1 : bvh_res.back();
    index_t brute_id = brute_res.empty() ? -1 : brute_res.back();
    GAUDI_ASSERT(bvh_id >= 0);
    GAUDI_ASSERT(brute_id >= 0);

    const auto bvh_dist = arp::test_point_line(query, slice<2, decltype(edge_view)>(edge_view, bvh_id));
    const auto brute_dist = arp::test_point_line(query, slice<2, decltype(edge_view)>(edge_view, brute_id));
    GAUDI_EXPECT(std::abs(bvh_dist - brute_dist) < 1e-6);
  }
}

GAUDI_TEST(bvh_point_to_tri) {
  auto mesh = load_sphere_mesh();
  arp::simplex_set<3> tris(mesh.vertices, mesh.tri_adj);
  auto bvh = arp::bvh_tree<3>::create(tris);

  adjacency_view<std::vector<vec3>, std::vector<index_t>> tri_view(
      mesh.vertices, mesh.tri_adj);

  const real tol = 1000.0;
  auto queries = random_points_in_bounds(mesh.vertices, 10);
  for (const auto &q : queries) {
    std::array<vec3, 1> query = {q};
    auto bvh_res = bvh->get_nearest(query, tol);
    auto brute_res =
        arp::brute_force_nearest_auto<1, 3>(query, tri_view, tol);

    index_t bvh_id = bvh_res.empty() ? -1 : bvh_res.back();
    index_t brute_id = brute_res.empty() ? -1 : brute_res.back();
    GAUDI_ASSERT(bvh_id >= 0);
    GAUDI_ASSERT(brute_id >= 0);

    const auto bvh_dist = arp::test_point_tri(query, slice<3, decltype(tri_view)>(tri_view, bvh_id));
    const auto brute_dist = arp::test_point_tri(query, slice<3, decltype(tri_view)>(tri_view, brute_id));
    GAUDI_EXPECT(std::abs(bvh_dist - brute_dist) < 1e-6);
  }
}

// ============================================================================
// New tests using bvh_tree_t<SimplexType> (type-based templating)
// ============================================================================

GAUDI_TEST(bvh_tree_t_point_to_edge) {
  auto mesh = load_sphere_mesh();

  // Use the new bvh_tree_t templated on SimplexType
  using edge_view = permuted_simplex_view<2, std::vector<vec3>, std::vector<index_t>>;
  arp::bvh_tree_t<edge_view> bvh(mesh.vertices, mesh.edge_adj);

  // Verify construction
  GAUDI_ASSERT(bvh.indices_.size() > 0);
  GAUDI_ASSERT(bvh.internal_nodes_.size() > 0);
  GAUDI_ASSERT(bvh.leaf_nodes_.size() > 0);

  // Verify stride extraction
  GAUDI_EXPECT(arp::bvh_tree_t<edge_view>::N == 2);
}

GAUDI_TEST(bvh_tree_t_point_to_tri) {
  auto mesh = load_sphere_mesh();

  // Use the new bvh_tree_t templated on SimplexType
  using tri_view = permuted_simplex_view<3, std::vector<vec3>, std::vector<index_t>>;
  arp::bvh_tree_t<tri_view> bvh(mesh.vertices, mesh.tri_adj);

  // Verify construction
  GAUDI_ASSERT(bvh.indices_.size() > 0);
  GAUDI_ASSERT(bvh.internal_nodes_.size() > 0);
  GAUDI_ASSERT(bvh.leaf_nodes_.size() > 0);

  // Verify stride extraction
  GAUDI_EXPECT(arp::bvh_tree_t<tri_view>::N == 3);
}


// ---------------------------------------------------------------------------
// bvh_tree_t brute-force vs BVH ID agreement for each simplex type (1,2,3)
//
// We generate query points at a known delta from a specific simplex so the
// ground-truth nearest ID is deterministic (no tie-breaking ambiguity).
// Both BVH and brute force must return the same simplex index.
// ---------------------------------------------------------------------------

// Brute-force O(n) nearest scan over a permuted_simplex_view.
// Returns {original_simplex_index, distance}.
template <size_t N>
std::pair<index_t, real> brute_nearest_psv(
    const std::array<vec3, 1> &query,
    const permuted_simplex_view<N, std::vector<vec3>, std::vector<index_t>> &view) {
  index_t best_id = -1;
  real best_dist = std::numeric_limits<real>::max();
  for (size_t i = 0; i < view.size(); ++i) {
    auto simplex = view[i];
    real d;
    if constexpr (N == 1)
      d = arp::test_point_point_tuple(query, simplex);
    else if constexpr (N == 2)
      d = arp::test_point_line_tuple(query, simplex);
    else if constexpr (N == 3)
      d = arp::test_point_tri_tuple(query, simplex);
    if (d < best_dist) {
      best_dist = d;
      best_id = view.get_index(i);
    }
  }
  return {best_id, best_dist};
}

// Generate a point at distance |delta| from a triangle's centroid along
// its face normal.  Returns {query_point, expected_simplex_index}.
template <size_t N>
std::pair<vec3, index_t> offset_from_tri(
    size_t walk_idx, real delta,
    const permuted_simplex_view<N, std::vector<vec3>, std::vector<index_t>> &view) {
  static_assert(N == 3);
  auto tri = view[walk_idx];
  vec3 centroid = (tri[0] + tri[1] + tri[2]) / 3.0;
  vec3 e1 = tri[1] - tri[0];
  vec3 e2 = tri[2] - tri[0];
  vec3 n = e1.cross(e2).normalized();
  return {centroid + delta * n, view.get_index(walk_idx)};
}

// Return the edge midpoint itself — expected distance is 0.
// Computing a proper edge normal requires adjacent triangles, so we
// just query a point that lies exactly on the edge.
template <size_t N>
std::pair<vec3, index_t> point_on_edge(
    size_t walk_idx,
    const permuted_simplex_view<N, std::vector<vec3>, std::vector<index_t>> &view) {
  static_assert(N == 2);
  auto edge = view[walk_idx];
  vec3 mid = (edge[0] + edge[1]) * 0.5;
  return {mid, view.get_index(walk_idx)};
}

// Generate a point at distance |delta| from a vertex along +Z (or +X if
// degenerate).
template <size_t N>
std::pair<vec3, index_t> offset_from_point(
    size_t walk_idx, real delta,
    const permuted_simplex_view<N, std::vector<vec3>, std::vector<index_t>> &view) {
  static_assert(N == 1);
  auto pt = view[walk_idx];
  return {pt[0] + delta * vec3::UnitZ(), view.get_index(walk_idx)};
}

GAUDI_TEST(bvh_tree_t_brute_vs_bvh_point) {
  auto mesh = load_sphere_mesh();

  using pt_view = permuted_simplex_view<1, std::vector<vec3>, std::vector<index_t>>;
  arp::bvh_tree_t<pt_view> bvh(mesh.vertices, mesh.point_adj);
  const auto &view = *bvh.permuted_data_view_;

  // Small delta so the generated query is closest to the source simplex
  // and not to any neighbor.
  const real delta = 1e-4;
  std::mt19937 rng(42);
  std::uniform_int_distribution<size_t> pick(0, view.size() - 1);

  for (int trial = 0; trial < 20; ++trial) {
    auto [q, expected_id] = offset_from_point<1>(pick(rng), delta, view);
    std::array<vec3, 1> query = {q};

    index_t bvh_id = bvh.find_nearest<Singulus<1>>(query);
    auto [brute_id, brute_dist] = brute_nearest_psv<1>(query, view);

    GAUDI_EXPECT(bvh_id == expected_id);
    GAUDI_EXPECT(brute_id == expected_id);
    GAUDI_EXPECT(std::abs(brute_dist - delta) < 1e-6);
  }
}

GAUDI_TEST(bvh_tree_t_brute_vs_bvh_edge) {
  auto mesh = load_sphere_mesh();

  using edge_view = permuted_simplex_view<2, std::vector<vec3>, std::vector<index_t>>;
  arp::bvh_tree_t<edge_view> bvh(mesh.vertices, mesh.edge_adj);
  const auto &view = *bvh.permuted_data_view_;

  std::mt19937 rng(42);
  std::uniform_int_distribution<size_t> pick(0, view.size() - 1);

  for (int trial = 0; trial < 20; ++trial) {
    auto [q, expected_id] = point_on_edge<2>(pick(rng), view);
    std::array<vec3, 1> query = {q};

    index_t bvh_id = bvh.find_nearest<Singulus<1>>(query);
    auto [brute_id, brute_dist] = brute_nearest_psv<2>(query, view);

    GAUDI_EXPECT(bvh_id == expected_id);
    GAUDI_EXPECT(brute_id == expected_id);
    GAUDI_EXPECT(brute_dist < 1e-10);
  }
}

GAUDI_TEST(bvh_tree_t_brute_vs_bvh_tri) {
  auto mesh = load_sphere_mesh();

  using tri_view = permuted_simplex_view<3, std::vector<vec3>, std::vector<index_t>>;
  arp::bvh_tree_t<tri_view> bvh(mesh.vertices, mesh.tri_adj);
  const auto &view = *bvh.permuted_data_view_;

  const real delta = 1e-4;
  std::mt19937 rng(42);
  std::uniform_int_distribution<size_t> pick(0, view.size() - 1);

  for (int trial = 0; trial < 20; ++trial) {
    auto [q, expected_id] = offset_from_tri<3>(pick(rng), delta, view);
    std::array<vec3, 1> query = {q};

    index_t bvh_id = bvh.find_nearest<Singulus<1>>(query);
    auto [brute_id, brute_dist] = brute_nearest_psv<3>(query, view);

    GAUDI_EXPECT(bvh_id == expected_id);
    GAUDI_EXPECT(brute_id == expected_id);
    GAUDI_EXPECT(std::abs(brute_dist - delta) < 1e-6);
  }
}

} // namespace test
} // namespace gaudi

#endif // __GAUDI_BVH_TESTS_HPP__
