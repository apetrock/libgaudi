#ifndef GAUDI_TEST_SHELL_TRAVERSAL_TESTS_HPP
#define GAUDI_TEST_SHELL_TRAVERSAL_TESTS_HPP

#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/test/test.hpp"

#include <vector>

namespace gaudi {
namespace test {

namespace shell_traversal_detail {

using namespace asawa::shell;

inline std::vector<shell::ptr> meshes() {
  return {load_tet(), load_torus(1.0, 0.33, 12, 8)};
}

/// Outward direction at `p` on the `load_torus` torus (z axis, radius `R`).
inline vec3 torus_out(const vec3 &p, real R) {
  vec3 core(p[0], p[1], 0);
  core = R * core.normalized();
  return (p - core).normalized();
}

} // namespace shell_traversal_detail

GAUDI_TEST(shell_traversal_rotation_identities) {
  using namespace asawa::shell;
  for (const auto &M : shell_traversal_detail::meshes()) {
    for (int i = 0; i < static_cast<int>(M->corner_count()); ++i) {
      const CornerId c = corner_id(i);
      if (M->next(c) < 0)
        continue;
      GAUDI_EXPECT(M->vnext(M->vprev(c)) == c);
      GAUDI_EXPECT(M->vprev(M->vnext(c)) == c);
      GAUDI_EXPECT(M->vert(M->vnext(c)) == M->vert(c));
      GAUDI_EXPECT(M->vert(M->vprev(c)) == M->vert(c));
      GAUDI_EXPECT(M->vert(M->other(c)) == M->vert(M->next(c)));
      GAUDI_EXPECT(M->prev(c) == M->other(M->vnext(c)));
    }
  }
}

GAUDI_TEST(shell_traversal_rings_and_ends) {
  using namespace asawa::shell;
  for (const auto &M : shell_traversal_detail::meshes()) {
    for (VertId v : M->get_vert_range()) {
      std::vector<CornerId> walked;
      M->for_each_vertex(v, [&](CornerId c, shell &) { walked.push_back(c); });
      GAUDI_ASSERT(!walked.empty());
      GAUDI_EXPECT(walked.size() == M->vsize(v));
      GAUDI_EXPECT(walked.front() == M->vbegin(v));
      GAUDI_EXPECT(walked.back() == M->vend(v));
      GAUDI_EXPECT(M->vnext(M->vend(v)) == M->vbegin(v));
      const std::vector<VertId> ring = M->get_one_ring(v);
      const std::vector<CornerId> link = M->link_corners(v);
      GAUDI_ASSERT(ring.size() == walked.size() && link.size() == walked.size());
      const size_t n = walked.size();
      for (size_t k = 0; k < n; ++k) {
        GAUDI_EXPECT(M->vert(walked[k]) == v);
        GAUDI_EXPECT(k == 0 || M->vnext(walked[k - 1]) == walked[k]);
        GAUDI_EXPECT(ring[k] == M->vert(M->other(walked[k])));
        GAUDI_EXPECT(link[k] == M->next(walked[k]));
        GAUDI_EXPECT(M->vert(link[k]) == ring[k]);
        GAUDI_EXPECT(M->vert(M->other(link[k])) == ring[(k + 1) % n]);
      }
    }
    for (FaceId f : M->get_face_range(false)) {
      std::vector<CornerId> walked;
      M->for_each_face(f, [&](CornerId c, shell &) { walked.push_back(c); });
      GAUDI_ASSERT(!walked.empty());
      GAUDI_EXPECT(walked.front() == M->fbegin(f));
      GAUDI_EXPECT(walked.back() == M->fend(f));
      GAUDI_EXPECT(static_cast<int>(walked.size()) == face_len(*M, M->fbegin(f)));
      for (CornerId c : walked)
        GAUDI_EXPECT(M->face(c) == f);
    }
  }
}

/// The loaders wind faces counterclockwise about the outward normal, so
/// `vnext` turns counterclockwise about it too: face(c) is (v, w_k, w_{k+1}).
GAUDI_TEST(shell_traversal_orientation_is_outward_ccw) {
  using namespace asawa::shell;
  const real R = 1.0;
  auto M = load_torus(R, 0.33, 12, 8);
  const std::vector<vec3> &x = asawa::const_get_vec_data(*M, 0);
  int bad = 0;
  for (int i = 0; i < static_cast<int>(M->corner_count()); ++i) {
    const CornerId c = corner_id(i);
    if (M->next(c) < 0)
      continue;
    const vec3 &p = x[M->vert(c)];
    const vec3 &a = x[M->vert(M->other(c))];
    const vec3 &b = x[M->vert(M->other(M->vnext(c)))];
    if (!((a - p).cross(b - p).dot(shell_traversal_detail::torus_out(p, R)) > 0))
      ++bad;
  }
  GAUDI_EXPECT(bad == 0);
}

} // namespace test
} // namespace gaudi

#endif
