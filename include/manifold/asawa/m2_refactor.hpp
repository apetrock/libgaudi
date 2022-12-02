#include "manifold/vec_addendum.h"

#include "GaudiGraphics/geometry_logger.h"

#include <vector>
#include <zlib.h>

#ifndef __ASAWA_MANIFOLD_SURFACE__
#define __ASAWA_MANIFOLD_SURFACE__
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

class manifold {
public:
  manifold(const std::vector<index_t> &corners_next,
           const std::vector<index_t> &corners_vert,
           const std::vector<index_t> &corners_face)
      : __corners_next(corners_next), __corners_vert(corners_vert),
        __corners_face(corners_face) {
    this->update_head();
    this->update_prev();
  }

  void update_head() {
    auto update = [](std::vector<index_t> &p_beg,
                     const std::vector<index_t> &c_ptr) {
      int size = *std::max_element(c_ptr.begin(), c_ptr.end());

      p_beg = std::vector<index_t>(size + 1, -1);

      for (int i = 0; i < c_ptr.size(); i++) {
        p_beg[c_ptr[i]] = i;
      }
    };

    update(__face_begin, __corners_face);
    update(__vert_begin, __corners_vert);
  }

  void update_prev() {
    // can't use the iterator because the iterator needs prev.
    __corners_prev = std::vector<index_t>(__corners_next.size(), -1);
    for (int i = 0; i < __face_begin.size(); i++) {
      int j_end = __face_begin[i];
      int j0 = __corners_next[j_end];
      bool it = true;
      while (it) {
        it = j0 != j_end;
        __corners_prev[__corners_next[j0]] = j0;
        j0 = __corners_next[j0];
      }
    }
  }
  /*
  void pack() {
    // TODO: safer pack, iterating from mRecycle[i] to mRecycle[i+1]
    // if (mFaceRecycle.size() > 0) {
    face_array tFaces;

    int j = 0;
    for (int i = 0; i < mFaces.size(); i++) {
      if (mHasFace[i]) {
        tFaces.push_back(mFaces[i]);
        tFaces.back()->position_in_set() = j;
        j++;
      }
    }
    mHasFace = std::vector<bool>(tFaces.size(), true);
    swap(mFaces, tFaces);
  }
  */
  size_t vert_count() { return __vert_begin.size(); }
  size_t face_count() { return __face_begin.size(); }
  size_t corner_count() { return __corners_next.size(); }

  index_t other(int id) { return 2 * (id / 2) + (id + 1) % 2; };
  index_t vnext(index_t id) { return this->next(this->other(id)); }
  index_t vprev(index_t id) { return this->other(this->prev(id)); }

  index_t next(index_t id) { return __corners_next[id]; }
  index_t prev(index_t id) { return __corners_prev[id]; }
  index_t vert(index_t id) { return __corners_vert[id]; }
  index_t face(index_t id) { return __corners_face[id]; }

  index_t fbegin(index_t id) { return __face_begin[id]; }
  index_t fend(index_t id) { return prev(__face_begin[id]); }
  index_t vbegin(index_t id) { return __vert_begin[id]; }
  index_t vend(index_t id) { return vprev(__vert_begin[id]); }

  void set_vbegin(index_t id, index_t c) { __vert_begin[id] = c; }
  void set_fbegin(index_t id, index_t c) { __face_begin[id] = c; }

  void set_next(index_t id, index_t c) { __corners_next[id] = c; }
  void set_prev(index_t id, index_t c) { __corners_prev[id] = c; }

  void set_vert(index_t id, index_t v) {
    __corners_vert[id] = v;
    set_vbegin(v, id);
  }
  void set_face(index_t id, index_t f) {
    __corners_face[id] = f;
    set_fbegin(f, id);
  }

  void link(index_t c0, index_t c1) {
    set_next(c0, c1);
    set_prev(c1, c0);
  }

  index_t insert_vertex() {
    // returns id of new vertex
    __vert_begin.push_back(-1);
    return __vert_begin.size() - 1;
  }

  index_t insert_face() {
    // returns id of new vertex
    __face_begin.push_back(-1);
    return __face_begin.size() - 1;
  }

  index_t insert_edge_pair() {
    // returns id of new vertex
    __corners_next.push_back(-1);
    __corners_next.push_back(-1);
    __corners_prev.push_back(-1);
    __corners_prev.push_back(-1);
    __corners_vert.push_back(-1);
    __corners_vert.push_back(-1);
    __corners_face.push_back(-1);
    __corners_face.push_back(-1);

    return __corners_next.size() - 2;
  }

  void remove_vertex(index_t i) { __vert_begin[i] = -1; }

  void remove_face(index_t i) { __face_begin[i] = -1; }

  void remove_edge_pair(index_t i0) {
    // returns id of new vertex
    index_t i1 = other(i0);
    __corners_next[i0] = -1;
    __corners_next[i1] = -1;

    __corners_prev[i0] = -1;
    __corners_prev[i1] = -1;

    __corners_vert[i0] = -1;
    __corners_vert[i1] = -1;

    __corners_face[i0] = -1;
    __corners_face[i1] = -1;
  }

  void for_each_face(index_t i,
                     std::function<void(index_t cid, manifold &m)> func) {
    int j0 = this->fbegin(i);
    std::cout << "i: " << i << " j0 " << j0 << std::endl;
    int j_end = this->fend(i);
    bool it = true;
    while (it) {
      it = j0 != j_end;
      func(j0, *this);
      j0 = this->next(j0);
    }
  }

  void for_each_vertex(index_t i,
                       std::function<void(index_t cid, manifold &m)> func) {
    int j0 = this->vbegin(i);
    int j_end = this->vend(i);
    bool it = true;
    while (it) {
      it = j0 != j_end;
      func(j0, *this);
      j0 = this->vnext(j0);
    }
  }

  void fupdate(index_t f) {
    std::cout << " f" << f << ": ";
    for_each_face(f, [f](index_t cid, manifold &m) {
      std::cout << " " << cid;
      m.set_face(cid, f);
    });
    std::cout << std::endl;
  }

  void vupdate(index_t v) {
    for_each_vertex(v, [v](index_t cid, manifold &m) { m.set_vert(cid, v); });
  }

  std::vector<index_t> __corners_next;
  std::vector<index_t> __corners_prev;
  std::vector<index_t> __corners_vert;
  std::vector<index_t> __corners_face;

  std::vector<index_t> __vert_begin;
  std::vector<index_t> __face_begin;
};

} // namespace asawa
#endif