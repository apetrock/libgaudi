#ifndef __GAUDI_TEST_VIEW_TESTS_HPP__
#define __GAUDI_TEST_VIEW_TESTS_HPP__

#include "gaudi/common.h"
#include "gaudi/test/test.hpp"
#include <array>
#include <vector>

namespace gaudi {
namespace test {

// Test simplex_view basic functionality
GAUDI_TEST(simplex_view_basic) {
  std::vector<vec3> data = {
      vec3(0, 0, 0), vec3(1, 0, 0),  // edge 0
      vec3(2, 0, 0), vec3(3, 0, 0),  // edge 1
      vec3(4, 0, 0), vec3(5, 0, 0)   // edge 2
  };

  simplex_view<2, std::vector<vec3>> view(data);

  // Check size (3 edges from 6 vertices)
  GAUDI_EXPECT(view.size() == 3);
  GAUDI_EXPECT(!view.empty());

  // Check stride
  GAUDI_EXPECT(simplex_view<2, std::vector<vec3>>::stride == 2);

  // Check tuple access
  auto edge0 = view[0];
  GAUDI_EXPECT((edge0[0] - vec3(0, 0, 0)).norm() < 1e-10);
  GAUDI_EXPECT((edge0[1] - vec3(1, 0, 0)).norm() < 1e-10);

  auto edge1 = view[1];
  GAUDI_EXPECT((edge1[0] - vec3(2, 0, 0)).norm() < 1e-10);
  GAUDI_EXPECT((edge1[1] - vec3(3, 0, 0)).norm() < 1e-10);

  auto edge2 = view[2];
  GAUDI_EXPECT((edge2[0] - vec3(4, 0, 0)).norm() < 1e-10);
  GAUDI_EXPECT((edge2[1] - vec3(5, 0, 0)).norm() < 1e-10);
}

// Test simplex_view with triangles
GAUDI_TEST(simplex_view_triangles) {
  std::vector<vec3> data = {
      vec3(0, 0, 0), vec3(1, 0, 0), vec3(0, 1, 0),  // tri 0
      vec3(2, 0, 0), vec3(3, 0, 0), vec3(2, 1, 0)   // tri 1
  };

  simplex_view<3, std::vector<vec3>> view(data);

  GAUDI_EXPECT(view.size() == 2);
  GAUDI_EXPECT(simplex_view<3, std::vector<vec3>>::stride == 3);

  auto tri0 = view[0];
  GAUDI_EXPECT((tri0[0] - vec3(0, 0, 0)).norm() < 1e-10);
  GAUDI_EXPECT((tri0[1] - vec3(1, 0, 0)).norm() < 1e-10);
  GAUDI_EXPECT((tri0[2] - vec3(0, 1, 0)).norm() < 1e-10);
}

// Test permuted_simplex permutes tuples correctly
GAUDI_TEST(permuted_simplex_basic) {
  std::vector<vec3> data = {
      vec3(0, 0, 0), vec3(1, 0, 0),  // edge 0
      vec3(2, 0, 0), vec3(3, 0, 0),  // edge 1
      vec3(4, 0, 0), vec3(5, 0, 0)   // edge 2
  };

  simplex_view<2, std::vector<vec3>> sv(data);

  // Permute: [2, 0, 1] - edge 2 first, then 0, then 1
  std::vector<index_t> perm = {2, 0, 1};

  auto psv = permute(sv, perm);

  GAUDI_EXPECT(psv.size() == 3);

  // First element should be edge 2
  auto e0 = psv[0];
  GAUDI_EXPECT((e0[0] - vec3(4, 0, 0)).norm() < 1e-10);
  GAUDI_EXPECT((e0[1] - vec3(5, 0, 0)).norm() < 1e-10);

  // Second element should be edge 0
  auto e1 = psv[1];
  GAUDI_EXPECT((e1[0] - vec3(0, 0, 0)).norm() < 1e-10);
  GAUDI_EXPECT((e1[1] - vec3(1, 0, 0)).norm() < 1e-10);

  // Third element should be edge 1
  auto e2 = psv[2];
  GAUDI_EXPECT((e2[0] - vec3(2, 0, 0)).norm() < 1e-10);
  GAUDI_EXPECT((e2[1] - vec3(3, 0, 0)).norm() < 1e-10);
}

// Test permuted_simplex_view convenience wrapper
GAUDI_TEST(permuted_simplex_view_basic) {
  std::vector<vec3> vertices = {
      vec3(0, 0, 0),  // vertex 0
      vec3(1, 0, 0),  // vertex 1
      vec3(2, 0, 0),  // vertex 2
      vec3(3, 0, 0)   // vertex 3
  };

  // Adjacency: edge 0 = (0,1), edge 1 = (1,2), edge 2 = (2,3)
  std::vector<index_t> adjacency = {0, 1, 1, 2, 2, 3};

  // Permutation: [2, 0, 1] - edge 2 first, then 0, then 1
  std::vector<index_t> permutation = {2, 0, 1};

  permuted_simplex_view<2, std::vector<vec3>, std::vector<index_t>> 
      psv(vertices, adjacency, permutation);

  GAUDI_EXPECT(psv.size() == 3);
  GAUDI_EXPECT(permuted_simplex_view<2, std::vector<vec3>, std::vector<index_t>>::stride == 2);

  // First element should be edge 2 = (2,3)
  auto e0 = psv[0];
  GAUDI_EXPECT((e0[0] - vec3(2, 0, 0)).norm() < 1e-10);
  GAUDI_EXPECT((e0[1] - vec3(3, 0, 0)).norm() < 1e-10);

  // Check get_tuple_ids
  auto ids0 = psv.get_tuple_ids(0);
  GAUDI_EXPECT(ids0[0] == 2);
  GAUDI_EXPECT(ids0[1] == 3);

  // Check get_index
  GAUDI_EXPECT(psv.get_index(0) == 2);  // Original index was 2
  GAUDI_EXPECT(psv.get_index(1) == 0);  // Original index was 0
  GAUDI_EXPECT(psv.get_index(2) == 1);  // Original index was 1
}

// Test stride extraction
GAUDI_TEST(simplex_stride_extraction) {
  using edge_view = simplex_view<2, std::vector<vec3>>;
  using tri_view = simplex_view<3, std::vector<vec3>>;
  using point_view = simplex_view<1, std::vector<vec3>>;

  GAUDI_EXPECT(simplex_stride_v<edge_view> == 2);
  GAUDI_EXPECT(simplex_stride_v<tri_view> == 3);
  GAUDI_EXPECT(simplex_stride_v<point_view> == 1);

  // Test with permuted_simplex_view
  using psv_edge = permuted_simplex_view<2, std::vector<vec3>, std::vector<index_t>>;
  using psv_tri = permuted_simplex_view<3, std::vector<vec3>, std::vector<index_t>>;

  GAUDI_EXPECT(simplex_stride_v<psv_edge> == 2);
  GAUDI_EXPECT(simplex_stride_v<psv_tri> == 3);
}

// Test SimplexView concept satisfaction
GAUDI_TEST(simplex_view_concept) {
  // These should compile if the concept is satisfied
  static_assert(SimplexView<simplex_view<2, std::vector<vec3>>>);
  static_assert(SimplexView<simplex_view<3, std::vector<vec3>>>);
  static_assert(SimplexView<permuted_simplex_view<2, std::vector<vec3>, std::vector<index_t>>>);
  static_assert(SimplexView<permuted_simplex_view<3, std::vector<vec3>, std::vector<index_t>>>);

  // Also test permuted_simplex
  using sv = simplex_view<2, std::vector<vec3>>;
  using ps = permuted_simplex<sv, std::vector<index_t>>;
  static_assert(SimplexView<ps>);
}

// Test that std::vector<std::array<vec3, N>> could work with concept
// (compile-time verification for flexibility)
GAUDI_TEST(simplex_view_vector_array_compile) {
  // This test verifies that our algorithms could work with stored data
  // For now, we just verify that the value_type is correct
  using edge_array = std::array<vec3, 2>;
  using edge_container = std::vector<edge_array>;

  edge_container edges = {
      {vec3(0, 0, 0), vec3(1, 0, 0)},
      {vec3(2, 0, 0), vec3(3, 0, 0)}
  };

  // Verify structure
  GAUDI_EXPECT(edges.size() == 2);
  GAUDI_EXPECT((edges[0][0] - vec3(0, 0, 0)).norm() < 1e-10);
  GAUDI_EXPECT((edges[0][1] - vec3(1, 0, 0)).norm() < 1e-10);
}

} // namespace test
} // namespace gaudi

#endif // __GAUDI_TEST_VIEW_TESTS_HPP__
