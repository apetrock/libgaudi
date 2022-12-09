
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
  vertices = {
      vec3(1, 0, -isq2),  // 0
      vec3(-1, 0, -isq2), // 1
      vec3(0, 1, isq2),   // 2
      vec3(0, -1, isq2),  // 3
  };

  faces = {
      {1, 0, 2, 3}, //
      {0, 4, 6, 2}, //
      {4, 5, 7, 6}, //
      {1, 5, 4, 0}, //
  };
}

void assemble_table(const std::vector<vec3> &vertices,
                    const std::vector<std::vector<index_t>> &faces,
                    std::vector<index_t> &corners_next,
                    std::vector<index_t> &corners_prev,
                    std::vector<index_t> &corners_vert,
                    std::vector<index_t> &corners_face) {
  std::vector<std::vector<index_t>> face_corners;
  std::vector<index_t> flat_corners;

  std::vector<std::vector<index_t>> faces_on_vertices(vertices.size());

  index_t idx = 0;
  for (index_t i = 0; i < faces.size(); i++) {
    const std::vector<index_t> &face = faces[i];
    std::vector<index_t> corners(face);

    for (index_t j = 0; j < face.size(); j++) {
      faces_on_vertices[face[j]].push_back(i);
      index_t corner_idx = idx++;
      corners[j] = corner_idx;
      flat_corners.push_back(idx);
    }
    face_corners.push_back(corners);
  }
  /*
  for (const auto &vert : faces_on_vertices) {
    for (const auto &face : vert) {
      std::cout << face << " ";
    }
    std::cout << std::endl;
  }
  */

  std::vector<index_t> interaction_list;
  for (index_t idx_vert = 0; idx_vert < faces_on_vertices.size(); idx_vert++) {
    const std::vector<index_t> &faces_on_vert = faces_on_vertices[idx_vert];
    for (index_t i = 0; i < faces_on_vert.size(); i++) {
      for (index_t j = i + 1; j < faces_on_vert.size(); j++) {
        interaction_list.push_back(faces_on_vert[i]);
        interaction_list.push_back(faces_on_vert[j]);
      }
    }
  }

  std::vector<index_t> corner_index(flat_corners);
  std::vector<bool> corner_allocated(flat_corners.size(), false);

  std::vector<index_t> corners(flat_corners.size(), -1);
  corners_next = std::vector<index_t>(flat_corners.size(), -1);
  corners_prev = std::vector<index_t>(flat_corners.size(), -1);
  corners_vert = std::vector<index_t>(flat_corners.size(), -1);
  corners_face = std::vector<index_t>(flat_corners.size(), -1);

  for (index_t k = 0; k < interaction_list.size(); k++) {
    std::cout << interaction_list[k] << " ";
  }
  std::cout << std::endl;

  int ii = 0;
  for (index_t k = 0; k < interaction_list.size(); k += 2) {
    std::cout << k + 0 << " " << interaction_list.size() << std::endl;
    index_t fi = interaction_list[k + 0];
    index_t fj = interaction_list[k + 1];

    for (int i = 0; i < faces[fi].size(); i++) {
      for (int j = 0; j < faces[fj].size(); j++) {
        index_t i0 = i;
        index_t i1 = (i + 1) % faces[fi].size();

        index_t vi0 = faces[fi][i0];
        index_t vi1 = faces[fi][i1];
        index_t fi0 = face_corners[fi][i0];
        index_t fi1 = face_corners[fi][i1];

        index_t j0 = j;
        index_t j1 = (j + 1) % face_corners[fj].size();

        index_t vj0 = faces[fj][j0];
        index_t vj1 = faces[fj][j1];
        index_t fj0 = face_corners[fj][j0];
        index_t fj1 = face_corners[fj][j1];

        if (corner_allocated[fi0])
          continue;
        if (corner_allocated[fj0])
          continue;

        if (vi0 == vj1 && vj0 == vi1) {
          corner_index[fi0] = 2 * ii + 0;
          corner_index[fj0] = 2 * ii + 1;
          corner_allocated[fi0] = true;
          corner_allocated[fj0] = true;
          corners[2 * ii + 0] = fi0;
          corners[2 * ii + 1] = fj0;
          corners_next[2 * ii + 0] = fi1;
          corners_next[2 * ii + 1] = fj1;
          corners_face[2 * ii + 0] = fi;
          corners_face[2 * ii + 1] = fj;
          corners_vert[2 * ii + 0] = vi0;
          corners_vert[2 * ii + 1] = vj0;
          ii++;
        }
      }
    }
  }

  std::vector<index_t> face_start(faces.size());
  std::vector<index_t> vert_start(vertices.size());

  for (int i = 0; i < corners.size(); i++) {
    corners[i] = corner_index[corners[i]];
  }

  for (int i = 0; i < corners_next.size(); i++) {
    corners_next[i] = corner_index[corners_next[i]];
  }

  for (int i = 0; i < corners_face.size(); i++) {
    face_start[corners_face[i]] = corners[i];
  }

  for (int i = 0; i < corners_vert.size(); i++) {
    vert_start[corners_vert[i]] = corners[i];
  }

  for (int i = 0; i < face_start.size(); i++) {
    int j_end = corners[face_start[i]];
    int j0 = corners_next[j_end];
    bool it = true;
    while (it) {
      it = j0 != j_end;
      corners_prev[corners_next[j0]] = j0;
      j0 = corners_next[j0];
    }
  }
}
} // namespace asawa

#endif