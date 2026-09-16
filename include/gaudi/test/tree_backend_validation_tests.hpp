#ifndef __GAUDI_TREE_BACKEND_VALIDATION_TESTS_HPP__
#define __GAUDI_TREE_BACKEND_VALIDATION_TESTS_HPP__

// Nearest-query equivalence: Morton bvh_tree vs legacy aabb_tree vs
// brute_force_tree, validated against O(n) brute_force reference.

#include "gaudi/arp/aabb.hpp"
#include "gaudi/arp/brute_force_tree.hpp"
#include "gaudi/arp/hash_tree.hpp"
#include "gaudi/arp/simplex_set.hpp"
#include "gaudi/asawa/objloader_refactor.hpp"
#include "gaudi/asawa/primitive_objects.hpp"
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/common.h"
#include "gaudi/paths.hpp"
#include "gaudi/test/arp_query_harness.hpp"
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
  const std::string file =
      gaudi::resolve_path("assets/models/obj/bunny.obj").string();
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

template <int N, typename BackendFactory>
void run_backend_nearest_suite(BackendValMesh &mesh, std::vector<index_t> &adj,
                               const std::vector<vec3> &queries,
                               BackendFactory &&make_tree) {
  arp::simplex_set<N> set(mesh.vertices, adj);
  auto tree = make_tree(set);
  run_nearest_equivalence_mesh<N, decltype(*tree)>(*tree, mesh.vertices, adj,
                                                   queries);
}

template <int N>
void run_all_backends_nearest(BackendValMesh &mesh, std::vector<index_t> &adj,
                              const std::vector<vec3> &queries) {
  run_backend_nearest_suite<N>(
      mesh, adj, queries, [](const arp::simplex_set<N> &set) {
        return arp::bvh_tree<N>::create(set);
      });
  run_backend_nearest_suite<N>(
      mesh, adj, queries, [](const arp::simplex_set<N> &set) {
        return arp::aabb_tree<N>::create(set);
      });
  run_backend_nearest_suite<N>(
      mesh, adj, queries, [](const arp::simplex_set<N> &set) {
        return arp::brute_force_tree<N>::create(set);
      });
}

inline void run_all_backend_validation(BackendValMesh &mesh, size_t n_queries,
                                       unsigned seed) {
  auto queries = backend_val_queries(mesh.vertices, n_queries, seed);
  run_all_backends_nearest<1>(mesh, mesh.point_adj, queries);
  run_all_backends_nearest<2>(mesh, mesh.edge_adj, queries);
  run_all_backends_nearest<3>(mesh, mesh.tri_adj, queries);
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

GAUDI_TEST(tree_backend_brute_force_tree_nearest) {
  auto mesh = load_val_sphere();
  auto queries = backend_val_queries(mesh.vertices, 20, 42);
  arp::simplex_set<3> set(mesh.vertices, mesh.tri_adj);
  auto tree = arp::brute_force_tree<3>::create(set);
  run_nearest_equivalence_mesh<3>(*tree, mesh.vertices, mesh.tri_adj, queries);
}

} // namespace test
} // namespace gaudi

#endif // __GAUDI_TREE_BACKEND_VALIDATION_TESTS_HPP__
