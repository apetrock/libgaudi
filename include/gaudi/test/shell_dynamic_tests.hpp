#ifndef __GAUDI_SHELL_DYNAMIC_TESTS_HPP__
#define __GAUDI_SHELL_DYNAMIC_TESTS_HPP__

#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/dynamic.hpp"
#include "gaudi/test/test.hpp"
#include "gaudi/vec_addendum.h"
#include <limits>
#include <vector>

namespace gaudi {
namespace test {

namespace shell_dynamic_test_detail {

inline index_t brute_force_nearest_tri(
    const std::vector<index_t> &face_verts_m, const std::vector<vec3> &x_m,
    const std::vector<index_t> &face_map_m, index_t q_vert,
    const std::vector<vec3> &x_t, real tol) {
  const int nf = static_cast<int>(face_verts_m.size() / 3);
  const std::vector<index_t> t_inds = {q_vert};
  index_t best = -1;
  real best_d = std::numeric_limits<real>::max();
  for (int fi = 0; fi < nf; ++fi) {
    real d = arp::pnt_tri_min(0, t_inds, x_t, fi, face_verts_m, x_m);
    if (d < tol && d < best_d) {
      best_d = d;
      best = face_map_m[fi];
    }
  }
  return best;
}

} // namespace shell_dynamic_test_detail

GAUDI_TEST(shell_dynamic_pnt_tri_collisions) {
  using shell_dynamic_test_detail::brute_force_nearest_tri;

  auto M = asawa::shell::load_tet();
  GAUDI_ASSERT(M != nullptr);

  auto x_datum = static_pointer_cast<asawa::vec3_datum>(M->get_datum(0));
  std::vector<vec3> x_mesh = x_datum->data();

  std::vector<index_t> face_verts_m = M->get_face_vert_ids(true);
  std::vector<index_t> face_map_m = M->get_face_map(true);
  GAUDI_ASSERT(face_verts_m.size() == 12);
  GAUDI_ASSERT(face_map_m.size() == 4);

  vec3 v0 = x_mesh[face_verts_m[0]];
  vec3 v1 = x_mesh[face_verts_m[1]];
  vec3 v2 = x_mesh[face_verts_m[2]];
  vec3 c = (v0 + v1 + v2) / 3.0;
  vec3 n = (v1 - v0).cross(v2 - v0);
  GAUDI_ASSERT(n.norm() > 1e-10);
  n.normalize();
  vec3 query = c + 0.05 * n;

  std::vector<vec3> x_t = x_mesh;
  const index_t q_vert = 0;
  x_t[q_vert] = query;

  auto dyn = asawa::shell::dynamic::create(M, 0.25, 2.0, 0.05);
  dyn->update_trees();

  const real tol = 1.0;
  std::vector<index_t> verts_t = {q_vert};
  std::vector<index_t> verts_map_t = {q_vert};
  auto hits =
      dyn->get_pnt_tri_collisions(verts_t, verts_map_t, x_t, *M, tol);

  GAUDI_ASSERT(hits.size() == 1);
  index_t expected_face =
      brute_force_nearest_tri(face_verts_m, x_mesh, face_map_m, q_vert, x_t,
                              tol);
  GAUDI_EXPECT(expected_face >= 0);
  GAUDI_EXPECT(hits[0][1] == expected_face);
}

} // namespace test
} // namespace gaudi

#endif // __GAUDI_SHELL_DYNAMIC_TESTS_HPP__
