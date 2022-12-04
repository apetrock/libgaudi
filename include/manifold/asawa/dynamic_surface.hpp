#include "m2.hpp"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "manifold/vec_addendum.h"

#include "GaudiGraphics/geometry_logger.h"

#include "m2_refactor.hpp"
#include "primitive_operations.hpp"

#include <vector>
#include <zlib.h>

#ifndef __ASAWA_DYNAMIC_SURFACE__
#define __ASAWA_DYNAMIC_SURFACE__
namespace asawa {
typedef double real;
typedef int index_t;
typedef Eigen::Matrix<real, 3, 1> vec3;
typedef Eigen::Matrix<real, 4, 1> vec4;

void op_edges(manifold &M,              //
              std::vector<vec3> &verts, //
              std::vector<index_t> &edges_to_op,
              std::function<index_t(index_t cid, manifold &m)> func) {

  for (auto d : M.get_data()) {
    if (d->type() != VERTEX)
      continue;
    d->alloc(edges_to_op.size());
  }

  for (index_t i = 0; i < edges_to_op.size(); i++) {
    for (auto d : M.get_data()) {
      if (d->type() != VERTEX)
        continue;
      d->calc(i, M, edges_to_op[i]);
    }
  }

  std::vector<index_t> post_edges(edges_to_op.size(), -1);
  for (index_t i = 0; i < edges_to_op.size(); i++) {
    index_t ic0 = edges_to_op[i];
    std::cout << ic0 << std::endl;
    post_edges[i] = func(ic0, M);
  }
  // verts.resize(M.vert_count());
  for (auto d : M.get_data()) {
    if (d->type() != VERTEX)
      continue;
    d->resize(M.vert_count());
  }
  for (index_t i = 0; i < post_edges.size(); i++) {
    for (auto d : M.get_data()) {
      if (d->type() != VERTEX)
        continue;
      d->map(M.vert(post_edges[i]), i);
    }
  }
}

void subdivide_edges(manifold &M, std::vector<vec3> &verts) {
  std::vector<index_t> edges_to_divide;
  edges_to_divide.push_back(3);
  edges_to_divide.push_back(5);
  edges_to_divide.push_back(8);
  edges_to_divide.push_back(13);
  edges_to_divide.push_back(21);
  op_edges(M, verts, edges_to_divide,
           [](index_t cid, manifold &m) { return subdivide_edge(m, cid); });
}

void collapse_edges(manifold &M, std::vector<vec3> &verts) {
  std::vector<index_t> edges_to_divide;
  edges_to_divide.push_back(4);
  // edges_to_divide.push_back(6);
  // edges_to_divide.push_back(9);
  //  edges_to_divide.push_back(14);
  edges_to_divide.push_back(22);
  op_edges(M, verts, edges_to_divide,
           [](index_t cid, manifold &m) { return collapse_edge(m, cid); });
}

void gather_edges_parallel(manifold &M, std::vector<vec3> &verts,
                           const real &eps) {
  std::vector<bool> _flags(M.face_count(), false);
  std::vector<bool> edges(M.corner_count() / 2, false);
  omp_lock_t writelock;
  omp_init_lock(&writelock);
#pragma omp parallel for shared(_flags, edges)
  for (int i = 0; i < M.corner_count(); i += 2) {
    int tid = omp_get_thread_num();
    index_t c0 = i;
    index_t c1 = M.other(i);
    index_t f0 = M.face(c0);
    index_t f1 = M.face(c1);
    omp_set_lock(&writelock);
    index_t ff0 = _flags[f0];
    index_t ff1 = _flags[f1];
    omp_unset_lock(&writelock);

    omp_set_lock(&writelock);
    if (!ff0 && !ff1) {
      _flags[f0] = true;
      _flags[f1] = true;
      edges[i / 2] = true;
    }
    omp_unset_lock(&writelock);
  }

  std::cout << "res: ";
  for (auto e : edges)
    std::cout << e << " ";
  std::cout << std::endl;
}

class dynamic_surface {
public:
  dynamic_surface() {}
};
} // namespace asawa
#endif