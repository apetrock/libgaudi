#ifndef __GAUDI_TEST_SIMPLEX_ALGORITHM_TESTS_HPP__
#define __GAUDI_TEST_SIMPLEX_ALGORITHM_TESTS_HPP__

#include "gaudi/arp/hash_tree.hpp"
#include "gaudi/arp/morton.hpp"
#include "gaudi/common.h"
#include "gaudi/test/test.hpp"
#include <cmath>
#include <vector>

namespace gaudi {
namespace test {

// Test calc_com with SimplexView (edges)
GAUDI_TEST(calc_com_simplex_view_edges) {
  std::vector<vec3> vertices = {
      vec3(0, 0, 0),  // 0
      vec3(2, 0, 0),  // 1
      vec3(0, 0, 0),  // 2
      vec3(0, 4, 0)   // 3
  };

  // Adjacency: edge 0 = (0,1), edge 1 = (2,3)
  std::vector<index_t> adjacency = {0, 1, 2, 3};
  std::vector<index_t> identity = {0, 1};

  permuted_simplex_view<2, std::vector<vec3>, std::vector<index_t>> 
      psv(vertices, adjacency, identity);

  auto coms = arp::calc_com(psv);

  GAUDI_EXPECT(coms.size() == 2);

  // Edge 0: midpoint at (1, 0, 0)
  auto [mass0, com0] = coms[0];
  GAUDI_EXPECT((com0 - vec3(1, 0, 0)).norm() < 1e-10);

  // Edge 1: midpoint at (0, 2, 0)
  auto [mass1, com1] = coms[1];
  GAUDI_EXPECT((com1 - vec3(0, 2, 0)).norm() < 1e-10);
}

// Test calc_com with SimplexView (triangles)
GAUDI_TEST(calc_com_simplex_view_triangles) {
  std::vector<vec3> vertices = {
      vec3(0, 0, 0),  // 0
      vec3(3, 0, 0),  // 1
      vec3(0, 3, 0)   // 2
  };

  std::vector<index_t> adjacency = {0, 1, 2};
  std::vector<index_t> identity = {0};

  permuted_simplex_view<3, std::vector<vec3>, std::vector<index_t>> 
      psv(vertices, adjacency, identity);

  auto coms = arp::calc_com(psv);

  GAUDI_EXPECT(coms.size() == 1);

  // Triangle centroid at (1, 1, 0)
  auto [mass, com] = coms[0];
  GAUDI_EXPECT((com - vec3(1, 1, 0)).norm() < 1e-10);

  // Mass = 0.5 * |cross(e1, e2)| = 0.5 * |(3,0,0) x (0,3,0)| = 0.5 * 9 = 4.5
  GAUDI_EXPECT(std::abs(mass - 4.5) < 1e-10);
}

// Test calc_extents with SimplexView (edges)
GAUDI_TEST(calc_extents_simplex_view_edges) {
  std::vector<vec3> vertices = {
      vec3(0, 0, 0),  // 0
      vec3(2, 1, 0),  // 1
      vec3(5, 5, 5),  // 2
      vec3(7, 7, 7)   // 3
  };

  std::vector<index_t> adjacency = {0, 1, 2, 3};
  std::vector<index_t> identity = {0, 1};

  permuted_simplex_view<2, std::vector<vec3>, std::vector<index_t>> 
      psv(vertices, adjacency, identity);

  auto extents = arp::calc_extents(psv);

  GAUDI_EXPECT(extents.size() == 2);

  // Edge 0: bounds [0,0,0] to [2,1,0]
  const auto &ext0 = extents[0];
  GAUDI_EXPECT((ext0[0] - vec3(0, 0, 0)).norm() < 1e-10);
  GAUDI_EXPECT((ext0[1] - vec3(2, 1, 0)).norm() < 1e-10);

  // Edge 1: bounds [5,5,5] to [7,7,7]
  const auto &ext1 = extents[1];
  GAUDI_EXPECT((ext1[0] - vec3(5, 5, 5)).norm() < 1e-10);
  GAUDI_EXPECT((ext1[1] - vec3(7, 7, 7)).norm() < 1e-10);
}

// Test calc_extents with SimplexView (triangles)
GAUDI_TEST(calc_extents_simplex_view_triangles) {
  std::vector<vec3> vertices = {
      vec3(0, 0, 0),  // 0
      vec3(3, 0, 0),  // 1
      vec3(1, 2, 1)   // 2
  };

  std::vector<index_t> adjacency = {0, 1, 2};
  std::vector<index_t> identity = {0};

  permuted_simplex_view<3, std::vector<vec3>, std::vector<index_t>> 
      psv(vertices, adjacency, identity);

  auto extents = arp::calc_extents(psv);

  GAUDI_EXPECT(extents.size() == 1);

  // Triangle bounds [0,0,0] to [3,2,1]
  const auto &ext = extents[0];
  GAUDI_EXPECT((ext[0] - vec3(0, 0, 0)).norm() < 1e-10);
  GAUDI_EXPECT((ext[1] - vec3(3, 2, 1)).norm() < 1e-10);
}

// Test map with SimplexView
GAUDI_TEST(map_simplex_view) {
  std::vector<vec3> vertices = {
      vec3(0, 0, 0), vec3(2, 0, 0),  // edge 0
      vec3(0, 0, 0), vec3(0, 4, 0)   // edge 1
  };

  std::vector<index_t> adjacency = {0, 1, 2, 3};
  std::vector<index_t> identity = {0, 1};

  permuted_simplex_view<2, std::vector<vec3>, std::vector<index_t>> 
      psv(vertices, adjacency, identity);

  // Map to compute edge lengths
  auto lengths = arp::map<real>(psv, 
      [](const std::array<vec3, 2> &edge, real default_val) {
        return (edge[1] - edge[0]).norm();
      }, 0.0);

  GAUDI_EXPECT(lengths.size() == 2);
  GAUDI_EXPECT(std::abs(lengths[0] - 2.0) < 1e-10);  // edge 0 length
  GAUDI_EXPECT(std::abs(lengths[1] - 4.0) < 1e-10);  // edge 1 length
}

// Test that SimplexView-based calc_com matches legacy calc_com
GAUDI_TEST(calc_com_simplex_view_matches_legacy) {
  std::vector<vec3> flat_data = {
      vec3(0, 0, 0), vec3(2, 0, 0),  // edge 0
      vec3(0, 0, 0), vec3(0, 4, 0),  // edge 1
      vec3(1, 1, 1), vec3(3, 3, 3)   // edge 2
  };

  // Legacy approach
  auto legacy_coms = arp::calc_com<2>(flat_data);

  // SimplexView approach
  simplex_view<2, std::vector<vec3>> sv(flat_data);
  auto simplex_coms = arp::calc_com(sv);

  GAUDI_EXPECT(legacy_coms.size() == simplex_coms.size());

  for (size_t i = 0; i < legacy_coms.size(); ++i) {
    auto [legacy_mass, legacy_com] = legacy_coms[i];
    auto [simplex_mass, simplex_com] = simplex_coms[i];

    GAUDI_EXPECT(std::abs(legacy_mass - simplex_mass) < 1e-10);
    GAUDI_EXPECT((legacy_com - simplex_com).norm() < 1e-10);
  }
}

// ---------------------------------------------------------------------------
// Mesh-topology simplex view tests
//
// Minimal mesh: 4 shared vertices, 2 triangles, 6 edges.
// Exercises permuted_simplex_view with separate adjacency lists on the
// same vertex buffer — the pattern used by shell.hpp / BVH construction.
// ---------------------------------------------------------------------------

// Helper: compare two simplices (arrays of vec3) for approximate equality,
// independent of vertex ordering within the simplex.
template <size_t N>
bool simplex_eq(const std::array<vec3, N> &a, const std::array<vec3, N> &b,
                double eps = 1e-10) {
  std::array<bool, N> matched{};
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < N; ++j) {
      if (!matched[j] && (a[i] - b[j]).norm() < eps) {
        matched[j] = true;
        break;
      }
    }
  }
  for (size_t i = 0; i < N; ++i)
    if (!matched[i])
      return false;
  return true;
}

GAUDI_TEST(mesh_topology_triangle_view) {
  std::vector<vec3> verts = {
      vec3(0, 0, 0), // 0
      vec3(1, 0, 0), // 1
      vec3(0.5, 1, 0), // 2
      vec3(1.5, 1, 0), // 3
  };
  std::vector<index_t> adj_tris = {0, 1, 2, 0, 2, 3};
  std::vector<index_t> identity = {0, 1};

  permuted_simplex_view<3, std::vector<vec3>, std::vector<index_t>>
      tri_view(verts, adj_tris, identity);

  GAUDI_EXPECT(tri_view.size() == 2);

  auto t0 = tri_view[0];
  GAUDI_EXPECT((t0[0] - verts[0]).norm() < 1e-10);
  GAUDI_EXPECT((t0[1] - verts[1]).norm() < 1e-10);
  GAUDI_EXPECT((t0[2] - verts[2]).norm() < 1e-10);

  auto t1 = tri_view[1];
  GAUDI_EXPECT((t1[0] - verts[0]).norm() < 1e-10);
  GAUDI_EXPECT((t1[1] - verts[2]).norm() < 1e-10);
  GAUDI_EXPECT((t1[2] - verts[3]).norm() < 1e-10);
}

GAUDI_TEST(mesh_topology_edge_view) {
  std::vector<vec3> verts = {
      vec3(0, 0, 0),
      vec3(1, 0, 0),
      vec3(0.5, 1, 0),
      vec3(1.5, 1, 0),
  };
  // 6 edges of the two-triangle mesh
  std::vector<index_t> adj_edges = {0, 1, 1, 2, 2, 0, 0, 2, 2, 3, 3, 0};
  std::vector<index_t> identity = {0, 1, 2, 3, 4, 5};

  permuted_simplex_view<2, std::vector<vec3>, std::vector<index_t>>
      edge_view(verts, adj_edges, identity);

  GAUDI_EXPECT(edge_view.size() == 6);

  // Spot-check a few edges
  auto e0 = edge_view[0]; // (0,1)
  GAUDI_EXPECT((e0[0] - verts[0]).norm() < 1e-10);
  GAUDI_EXPECT((e0[1] - verts[1]).norm() < 1e-10);

  auto e5 = edge_view[5]; // (3,0)
  GAUDI_EXPECT((e5[0] - verts[3]).norm() < 1e-10);
  GAUDI_EXPECT((e5[1] - verts[0]).norm() < 1e-10);
}

GAUDI_TEST(mesh_topology_shared_verts_different_adjacency) {
  // Same vertex buffer, two different adjacency lists (edges vs triangles)
  std::vector<vec3> verts = {
      vec3(0, 0, 0),
      vec3(1, 0, 0),
      vec3(0.5, 1, 0),
      vec3(1.5, 1, 0),
  };
  std::vector<index_t> adj_tris = {0, 1, 2, 0, 2, 3};
  std::vector<index_t> adj_edges = {0, 1, 1, 2, 2, 0, 0, 2, 2, 3, 3, 0};
  std::vector<index_t> tri_id = {0, 1};
  std::vector<index_t> edge_id = {0, 1, 2, 3, 4, 5};

  permuted_simplex_view<3, std::vector<vec3>, std::vector<index_t>>
      tri_view(verts, adj_tris, tri_id);
  permuted_simplex_view<2, std::vector<vec3>, std::vector<index_t>>
      edge_view(verts, adj_edges, edge_id);

  // Triangle centroids
  auto tri_coms = arp::calc_com(tri_view);
  GAUDI_EXPECT(tri_coms.size() == 2);
  vec3 c0 = (verts[0] + verts[1] + verts[2]) / 3.0;
  vec3 c1 = (verts[0] + verts[2] + verts[3]) / 3.0;
  GAUDI_EXPECT((std::get<1>(tri_coms[0]) - c0).norm() < 1e-10);
  GAUDI_EXPECT((std::get<1>(tri_coms[1]) - c1).norm() < 1e-10);

  // Edge midpoints
  auto edge_coms = arp::calc_com(edge_view);
  GAUDI_EXPECT(edge_coms.size() == 6);
  vec3 m0 = (verts[0] + verts[1]) * 0.5;
  GAUDI_EXPECT((std::get<1>(edge_coms[0]) - m0).norm() < 1e-10);
}

GAUDI_TEST(mesh_topology_morton_permutation_preserves_simplices) {
  std::vector<vec3> verts = {
      vec3(0, 0, 0),
      vec3(1, 0, 0),
      vec3(0.5, 1, 0),
      vec3(1.5, 1, 0),
  };
  std::vector<index_t> adj_edges = {0, 1, 1, 2, 2, 0, 0, 2, 2, 3, 3, 0};
  std::vector<index_t> identity = {0, 1, 2, 3, 4, 5};

  // Unpermuted view
  permuted_simplex_view<2, std::vector<vec3>, std::vector<index_t>>
      unpermuted(verts, adj_edges, identity);

  // Build Morton-sorted permutation from edge midpoint COMs
  auto [hashes, perm, internal_nodes, leaf_nodes] = arp::make_hash(unpermuted);

  // Permuted view
  permuted_simplex_view<2, std::vector<vec3>, std::vector<index_t>>
      permuted_view(verts, adj_edges, perm);

  GAUDI_EXPECT(permuted_view.size() == unpermuted.size());

  // O(n^2) check: every simplex in the unpermuted view must appear
  // exactly once in the permuted view (same set, possibly different order).
  size_t n = unpermuted.size();
  std::vector<bool> matched(n, false);
  for (size_t i = 0; i < n; ++i) {
    auto orig = unpermuted[i];
    bool found = false;
    for (size_t j = 0; j < n; ++j) {
      if (!matched[j] && simplex_eq<2>(orig, permuted_view[j])) {
        matched[j] = true;
        found = true;
        break;
      }
    }
    GAUDI_EXPECT(found);
  }
  // Verify no unmatched (no duplicates in permuted view)
  for (size_t i = 0; i < n; ++i) {
    GAUDI_EXPECT(matched[i]);
  }
}

GAUDI_TEST(mesh_topology_triangle_permutation_preserves_simplices) {
  std::vector<vec3> verts = {
      vec3(0, 0, 0),
      vec3(1, 0, 0),
      vec3(0.5, 1, 0),
      vec3(1.5, 1, 0),
  };
  std::vector<index_t> adj_tris = {0, 1, 2, 0, 2, 3};
  std::vector<index_t> identity = {0, 1};

  permuted_simplex_view<3, std::vector<vec3>, std::vector<index_t>>
      unpermuted(verts, adj_tris, identity);

  auto [hashes, perm, internal_nodes, leaf_nodes] = arp::make_hash(unpermuted);

  permuted_simplex_view<3, std::vector<vec3>, std::vector<index_t>>
      permuted_view(verts, adj_tris, perm);

  GAUDI_EXPECT(permuted_view.size() == unpermuted.size());

  size_t n = unpermuted.size();
  std::vector<bool> matched(n, false);
  for (size_t i = 0; i < n; ++i) {
    auto orig = unpermuted[i];
    bool found = false;
    for (size_t j = 0; j < n; ++j) {
      if (!matched[j] && simplex_eq<3>(orig, permuted_view[j])) {
        matched[j] = true;
        found = true;
        break;
      }
    }
    GAUDI_EXPECT(found);
  }
  for (size_t i = 0; i < n; ++i) {
    GAUDI_EXPECT(matched[i]);
  }
}

GAUDI_TEST(mesh_topology_vertex_degenerate_case) {
  std::vector<vec3> verts = {
      vec3(0, 0, 0),
      vec3(1, 0, 0),
      vec3(0, 1, 0),
      vec3(0, 0, 1),
  };
  // N=1: adjacency is just ascending indices
  std::vector<index_t> adj_pts = {0, 1, 2, 3};
  std::vector<index_t> identity = {0, 1, 2, 3};

  permuted_simplex_view<1, std::vector<vec3>, std::vector<index_t>>
      pt_view(verts, adj_pts, identity);

  GAUDI_EXPECT(pt_view.size() == 4);
  for (size_t i = 0; i < 4; ++i) {
    auto pt = pt_view[i];
    GAUDI_EXPECT((pt[0] - verts[i]).norm() < 1e-10);
  }
}

} // namespace test
} // namespace gaudi

#endif // __GAUDI_TEST_SIMPLEX_ALGORITHM_TESTS_HPP__
