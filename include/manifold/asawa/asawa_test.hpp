#include "m2.hpp"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "manifold/vec_addendum.h"

#include "GaudiGraphics/geometry_logger.h"

#include "dynamic_surface.hpp"
#include "m2_refactor.hpp"
#include "primitive_operations.hpp"

#include <vector>
#include <zlib.h>

#ifndef __M2REFACTOR_TEST__
#define __M2REFACTOR_TEST__
namespace asawa {
typedef double real;
typedef int index_t;
typedef Eigen::Matrix<real, 3, 1> vec3;
typedef Eigen::Matrix<real, 4, 1> vec4;

void debug_manifold(manifold &M, const std::vector<vec3> verts) {
  for (int i = 0; i < M.__corners_next.size(); i += 2) {
    if (M.__corners_next[i] < 0)
      continue;
    int i0 = i;
    int i1 = M.other(i0);
    int v0 = M.vert(i0);
    int v1 = M.vert(i1);
    gg::geometry_logger::line(verts[v0], verts[v1], vec4(0.5, 0.5, 0.5, 1.0));
  }
}

void test() {
  std::vector<vec3> vertices;
  std::vector<std::vector<int>> faces;
  make_cube(vertices, faces);
  std::vector<index_t> corners_next, corners_prev, corners_vert, corners_face;
  assemble_table(vertices, faces, corners_next, corners_prev, corners_vert,
                 corners_face);

  manifold *M = new manifold(corners_next, corners_vert, corners_face);

  triangulate(*M);
  subdivide_edges(*M, vertices);
  // collapse_edges(*M, vertices);
  gather_edges_parallel(*M, vertices, 1.0);
  // merge_face(*M, 3, M->other(3));
  debug_manifold(*M, vertices);
  delete M;
}

} // namespace asawa
#endif