
#ifndef __ASAWA_PRIM_OBS__
#define __ASAWA_PRIM_OBS__

#include <cassert>

#include "manifold/vec_addendum.h"

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

} // namespace asawa

#endif