
#ifndef __ASAWA_PRIM_OBS__
#define __ASAWA_PRIM_OBS__

#include <cassert>
#include <cmath>
#include <vector>

#include "gaudi/vec_addendum.h"

namespace gaudi {
namespace asawa {
typedef double real;
typedef int index_t;
typedef Eigen::Matrix<real, 3, 1> vec3;
typedef Eigen::Matrix<real, 4, 1> vec4;

void make_cube(std::vector<vec3> &vertices,
               std::vector<std::vector<int>> &faces) {

  /*
     1------5
   /       / \
  /       /   \
  0------4     7
  \       \   /
   \       \ /
     2------6

      z
     /
    o---x
     \
      y
  */
  vertices = {
      vec3(0, 0, 0), // 0
      vec3(0, 0, 1), // 1
      vec3(0, 1, 0), // 2
      vec3(0, 1, 1), // 3
      vec3(1, 0, 0), // 4
      vec3(1, 0, 1), // 5
      vec3(1, 1, 0), // 6
      vec3(1, 1, 1)  // 7
  };

  faces = {
      {1, 0, 2, 3}, //
      {0, 4, 6, 2}, //
      {4, 5, 7, 6}, //
      {1, 5, 4, 0}, //
      {5, 1, 3, 7}, //
      {7, 3, 2, 6},
  };
}

void make_tet(std::vector<vec3> &vertices,
              std::vector<std::vector<int>> &faces) {
  real isq2 = 1.0 / sqrt(2.0);
  //     0
  //    / \
  //   1 - 2
  //  / \ / \
  // 0 - 3 - 0

  vertices = {
      vec3(1, 0, -isq2),  // 0
      vec3(-1, 0, -isq2), // 1
      vec3(0, 1, isq2),   // 2
      vec3(0, -1, isq2),  // 3
  };

  faces = {
      {0, 1, 2}, //
      {1, 0, 3}, //
      {2, 1, 3},
      {3, 0, 2} //
  };
}

// Simple radial (UV) sphere centered at the origin.  Three parameters:
//   radius, u_segments (longitude slices), v_segments (latitude stacks).
// Triangulated throughout (poles fanned, interior bands split on a diagonal)
// so it is a pure triangle mesh like the other test primitives.  Winding is
// consistent (outward), so assemble_table builds a closed, oriented manifold.
// No subdivision -- an icosphere is left for the future.
inline void make_sphere(std::vector<vec3> &vertices,
                        std::vector<std::vector<int>> &faces,
                        real radius = 1.0, int u_segments = 24,
                        int v_segments = 16) {
  vertices.clear();
  faces.clear();
  if (u_segments < 3)
    u_segments = 3;
  if (v_segments < 2)
    v_segments = 2;

  const int top = 0;
  vertices.push_back(vec3(0.0, radius, 0.0));

  // Interior latitude rings r = 1 .. v_segments-1 (poles excluded).
  for (int r = 1; r < v_segments; ++r) {
    const real phi = M_PI * real(r) / real(v_segments); // 0 -> pi
    const real y = radius * std::cos(phi);
    const real ring_radius = radius * std::sin(phi);
    for (int s = 0; s < u_segments; ++s) {
      const real theta = 2.0 * M_PI * real(s) / real(u_segments);
      vertices.push_back(
          vec3(ring_radius * std::cos(theta), y, ring_radius * std::sin(theta)));
    }
  }

  const int bottom = int(vertices.size());
  vertices.push_back(vec3(0.0, -radius, 0.0));

  auto ring = [u_segments](int r, int s) {
    return 1 + (r - 1) * u_segments + (s % u_segments);
  };

  // Winding is CCW seen from outside so face_cross / vert_normal point outward.

  // Top cap.
  for (int s = 0; s < u_segments; ++s)
    faces.push_back({top, ring(1, s + 1), ring(1, s)});

  // Interior bands (two triangles per quad, split on the a-c diagonal).
  for (int r = 1; r < v_segments - 1; ++r)
    for (int s = 0; s < u_segments; ++s) {
      const int a = ring(r, s), b = ring(r + 1, s);
      const int c = ring(r + 1, s + 1), d = ring(r, s + 1);
      faces.push_back({a, c, b});
      faces.push_back({a, d, c});
    }

  // Bottom cap.
  const int last = v_segments - 1;
  for (int s = 0; s < u_segments; ++s)
    faces.push_back({bottom, ring(last, s), ring(last, s + 1)});
}

// Torus centered at the origin, axis = +Z.  Parameters mirror the offset-torus
// fixture in test/darboux_cyclide_torus_fixture.hpp: major_radius (ring),
// minor_radius (tube), and segment counts around each.  Genus-1, closed and
// consistently wound -- a good non-convex stress case for proximity queries.
inline void make_torus(std::vector<vec3> &vertices,
                       std::vector<std::vector<int>> &faces,
                       real major_radius = 1.0, real minor_radius = 0.33,
                       int major_segments = 48, int minor_segments = 24) {
  vertices.clear();
  faces.clear();
  if (major_segments < 3)
    major_segments = 3;
  if (minor_segments < 3)
    minor_segments = 3;

  for (int i = 0; i < major_segments; ++i) {
    const real u = 2.0 * M_PI * real(i) / real(major_segments);
    const real cu = std::cos(u), su = std::sin(u);
    for (int j = 0; j < minor_segments; ++j) {
      const real v = 2.0 * M_PI * real(j) / real(minor_segments);
      const real rho = major_radius + minor_radius * std::cos(v);
      vertices.push_back(
          vec3(rho * cu, rho * su, minor_radius * std::sin(v)));
    }
  }

  auto vid = [major_segments, minor_segments](int i, int j) {
    return ((i % major_segments + major_segments) % major_segments) *
               minor_segments +
           ((j % minor_segments + minor_segments) % minor_segments);
  };

  for (int i = 0; i < major_segments; ++i) {
    for (int j = 0; j < minor_segments; ++j) {
      const int i1 = (i + 1) % major_segments;
      const int j1 = (j + 1) % minor_segments;
      faces.push_back({vid(i, j), vid(i1, j), vid(i1, j1)});
      faces.push_back({vid(i, j), vid(i1, j1), vid(i, j1)});
    }
  }
}

} // namespace asawa
} // namespace gaudi
#endif