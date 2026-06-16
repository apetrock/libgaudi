#ifndef __GAUDI_TREE_BACKEND_VALIDATION_TESTS_HPP__
#define __GAUDI_TREE_BACKEND_VALIDATION_TESTS_HPP__

// ---------------------------------------------------------------------------
// 3-way nearest-query equivalence: legacy aabb_tree vs Morton bvh_tree vs
// brute force, over the bunny (hard), the procedural sphere (smooth), and the
// procedural torus (non-convex genus-1).
//
// This is the USE_HASH equivalence check: both backends are instantiated
// directly here (regardless of the arp.h toggle), fed the SAME fixed-seed
// random query set, and compared against an O(n) brute-force scan for point,
// edge, and triangle simplices.
//
// Assertions (tie-tolerant):
//   primary  : each backend's nearest *distance* matches brute force within eps
//              (tie-proof: thin bunny features admit multiple equidistant ids).
//   secondary: id matches brute force OR the distances agree within eps
//              (catches a genuinely wrong pick, tolerates geometric ties).
// ---------------------------------------------------------------------------

#include "gaudi/arp/aabb.hpp"      // legacy backend (always available here)
#include "gaudi/arp/brute_force.hpp"
#include "gaudi/arp/hash_tree.hpp" // Morton backend
#include "gaudi/arp/pairwise_tests.hpp"
#include "gaudi/arp/simplex_set.hpp"
#include "gaudi/asawa/objloader_refactor.hpp"
#include "gaudi/asawa/primitive_objects.hpp"
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/common.h"
#include "gaudi/paths.hpp"
#include "gaudi/test/test.hpp"
#include "gaudi/vec_addendum.h"
#include <cmath>
#include <random>
#include <string>
#include <vector>

namespace gaudi {
namespace test {

struct BackendValMesh {
  std::vector<vec3> vertices;
  std::vector<index_t> point_adj;
  std::vector<index_t> edge_adj;
  std::vector<index_t> tri_adj;
};

inline BackendValMesh make_backend_val_mesh(
    const std::vector<vec3> &verts,
    const std::vector<std::vector<int>> &faces) {
  BackendValMesh m;
  m.vertices = verts;

  std::vector<index_t> corners_next, corners_vert, corners_face;
  asawa::assemble_table(verts, faces, corners_next, corners_vert, corners_face);
  auto shell =
      asawa::shell::shell::create(corners_next, corners_vert, corners_face);
  GAUDI_ASSERT(shell != nullptr);

  {
    auto vr = shell->get_vert_range();
    m.point_adj.assign(vr.begin(), vr.end());
  }
  m.edge_adj = shell->get_edge_vert_ids();
  m.tri_adj = shell->get_face_vert_ids(true);

  GAUDI_ASSERT(!m.point_adj.empty());
  GAUDI_ASSERT(!m.edge_adj.empty());
  GAUDI_ASSERT(!m.tri_adj.empty());
  GAUDI_ASSERT(m.tri_adj.size() % 3 == 0);
  return m;
}

inline BackendValMesh load_val_sphere() {
  std::vector<vec3> v;
  std::vector<std::vector<int>> f;
  asawa::make_sphere(v, f);
  return make_backend_val_mesh(v, f);
}

inline BackendValMesh load_val_torus() {
  std::vector<vec3> v;
  std::vector<std::vector<int>> f;
  asawa::make_torus(v, f);
  return make_backend_val_mesh(v, f);
}

inline BackendValMesh load_val_bunny() {
  std::vector<vec3> v;
  std::vector<std::vector<int>> f;
  const std::string file = gaudi::resolve_path("assets/bunny.obj").string();
  asawa::loadObjfile(file, v, f);
  GAUDI_ASSERT(!v.empty());
  GAUDI_ASSERT(!f.empty());
  return make_backend_val_mesh(v, f);
}

inline std::vector<vec3> backend_val_queries(const std::vector<vec3> &verts,
                                             size_t count, unsigned seed) {
  vec3 lo = verts.front();
  vec3 hi = verts.front();
  for (const auto &v : verts) {
    lo = va::min(lo, v);
    hi = va::max(hi, v);
  }
  std::mt19937 rng(seed);
  std::uniform_real_distribution<real> dx(lo[0], hi[0]);
  std::uniform_real_distribution<real> dy(lo[1], hi[1]);
  std::uniform_real_distribution<real> dz(lo[2], hi[2]);

  std::vector<vec3> q;
  q.reserve(count);
  for (size_t i = 0; i < count; ++i)
    q.emplace_back(dx(rng), dy(rng), dz(rng));
  return q;
}

template <int N, typename ViewT>
real backend_val_dist(const std::array<vec3, 1> &q, ViewT &view, index_t id) {
  if constexpr (N == 1)
    return arp::test_point_point(q, slice<1, ViewT>(view, id));
  else if constexpr (N == 2)
    return arp::test_point_line(q, slice<2, ViewT>(view, id));
  else
    return arp::test_point_tri(q, slice<3, ViewT>(view, id));
}

// Runs the 3-way comparison for one simplex stride N over one adjacency set.
template <int N>
void run_backend_validation(BackendValMesh &mesh, std::vector<index_t> &adj,
                            const std::vector<vec3> &queries) {
  arp::simplex_set<N> set(mesh.vertices, adj);
  auto morton = arp::bvh_tree<N>::create(set);
  auto legacy = arp::aabb_tree<N>::create(set);

  adjacency_view<std::vector<vec3>, std::vector<index_t>> view(mesh.vertices,
                                                              adj);
  using ViewT = decltype(view);

  const real tol = 1000.0; // contracting-radius / single-nearest mode
  const real eps = 1e-6;

  for (const auto &qp : queries) {
    std::array<vec3, 1> query = {qp};

    auto m_res = morton->get_nearest(query, tol);
    auto l_res = legacy->get_nearest(query, tol);
    auto b_res = arp::brute_force_nearest_auto<1, N>(query, view, tol);

    index_t mid = m_res.empty() ? -1 : m_res.back();
    index_t lid = l_res.empty() ? -1 : l_res.back();
    index_t bid = b_res.empty() ? -1 : b_res.back();

    GAUDI_ASSERT(mid >= 0);
    GAUDI_ASSERT(lid >= 0);
    GAUDI_ASSERT(bid >= 0);

    const real md = backend_val_dist<N, ViewT>(query, view, mid);
    const real ld = backend_val_dist<N, ViewT>(query, view, lid);
    const real bd = backend_val_dist<N, ViewT>(query, view, bid);

    // Primary, tie-proof invariant: same minimal distance as brute force.
    GAUDI_EXPECT(std::abs(md - bd) < eps);
    GAUDI_EXPECT(std::abs(ld - bd) < eps);

    // Secondary: id agreement, tolerating exact geometric ties.
    GAUDI_EXPECT(mid == bid || std::abs(md - bd) < eps);
    GAUDI_EXPECT(lid == bid || std::abs(ld - bd) < eps);
  }
}

inline void run_all_backend_validation(BackendValMesh &mesh, size_t n_queries,
                                       unsigned seed) {
  auto queries = backend_val_queries(mesh.vertices, n_queries, seed);
  run_backend_validation<1>(mesh, mesh.point_adj, queries);
  run_backend_validation<2>(mesh, mesh.edge_adj, queries);
  run_backend_validation<3>(mesh, mesh.tri_adj, queries);
}

GAUDI_TEST(tree_backend_3way_sphere) {
  auto mesh = load_val_sphere();
  run_all_backend_validation(mesh, 40, 0x5EEDu);
}

GAUDI_TEST(tree_backend_3way_torus) {
  auto mesh = load_val_torus();
  run_all_backend_validation(mesh, 40, 0xB16B00B5u);
}

GAUDI_TEST(tree_backend_3way_bunny) {
  auto mesh = load_val_bunny();
  run_all_backend_validation(mesh, 24, 0x1337u);
}

} // namespace test
} // namespace gaudi

#endif // __GAUDI_TREE_BACKEND_VALIDATION_TESTS_HPP__
