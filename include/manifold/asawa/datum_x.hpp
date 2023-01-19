#include <cassert>
#include <cstddef>
#include <cxxabi.h>

#include <execinfo.h>
#include <iostream>
#include <memory.h>
#include <ostream>
#include <stdio.h>
#include <vector>
#include <zlib.h>

#include "datums.hpp"
#include "manifold/geometry_types.hpp"
#include "manifold/vec_addendum.h"

#ifndef __ASAWA_X_DATUM__
#define __ASAWA_X_DATUM__

namespace asawa {

using real = double;
using vec3 = Eigen::Matrix<real, 3, 1>;
using vec4 = Eigen::Matrix<real, 4, 1>;
using real_datum = datum_t<real>;
using vec3_datum = datum_t<vec3>;

// break this out into new file at some point

real cotan(manifold &M, index_t ci, const std::vector<vec3> &x) {

  vec3 xp = x[M.vert(M.prev(ci))];
  vec3 x0 = x[M.vert(ci)];
  vec3 xn = x[M.vert(M.next(ci))];
  return va::cotan(x0, xp, xn);
}

vec3 face_normal(manifold &M, index_t fi, const std::vector<vec3> &x) {
  vec3 N = vec3::Zero();
  M.for_each_face_tri(
      fi, [&N, &x](index_t c0, index_t c1, index_t c2, asawa::manifold &M) {
        vec3 x0 = x[M.vert(c0)];
        vec3 x1 = x[M.vert(c1)];
        vec3 x2 = x[M.vert(c2)];
        N += (x1 - x0).cross(x2 - x0);
      });
  return N.normalized();
}

vec3 vert_normal(manifold &M, index_t vi, const std::vector<vec3> &x) {
  vec3 N = vec3::Zero();
  M.for_each_vertex(vi, [&N, &x](index_t ci, asawa::manifold &M) {
    N += face_normal(M, M.face(ci), x);
  });
  return N.normalized();
}

vec3 edge_normal(manifold &M, index_t c0, const std::vector<vec3> &x) {
  index_t c1 = M.other(c0);
  vec3 N = face_normal(M, M.face(c0), x);
  N += face_normal(M, M.face(c1), x);
  return N.normalized();
}

std::vector<vec3> face_normals(manifold &M, const std::vector<vec3> &x) {
  auto range = M.get_face_range();
  std::vector<vec3> Ns(range.size());
  int i = 0;
  for (auto vi : range) {
    Ns[i++] = face_normal(M, i, x);
  }
  return Ns;
}

std::vector<vec3> vertex_normals(manifold &M, const std::vector<vec3> &x) {
  auto range = M.get_vert_range();
  std::vector<vec3> Ns(range.size());
  int i = 0;
  for (auto i : range) {
    Ns[i++] = vert_normal(M, i, x);
  }
  return Ns;
}

vec3 edge_center(manifold &M, index_t c0, const std::vector<vec3> &x) {
  index_t c1 = M.other(c0);
  vec3 x0 = x[M.vert(c0)];
  vec3 x1 = x[M.vert(c1)];
  return 0.5 * (x0 + x1);
}
#if 0
template <typename SPACE>
std::vector<typename SPACE::coordinate_type>
gradient(const std::vector<typename SPACE::real> &u,
         typename surf<SPACE>::surf_ptr s) {
  M2_TYPEDEFS;
  ///////////////
  // gradient
  ///////////////
  face_array &faces = s->get_faces();
  edge_array &edges = s->get_edges();

  coordinate_array gradU(faces.size(), coordinate_type::Zero());

  for (auto e : edges) {
    coordinate_type c0 = ci::get_coordinate<SPACE>(e->v1());
    coordinate_type c1 = ci::get_coordinate<SPACE>(e->v2());
    real u0 = u[e->v1()->prev()->vertex()->position_in_set()];
    real u1 = u[e->v2()->prev()->vertex()->position_in_set()];

    real A0 = asawa::ci::area<SPACE>(e->v1()->face());
    real A1 = asawa::ci::area<SPACE>(e->v2()->face());

    coordinate_type N0 = asawa::ci::normal<SPACE>(e->v1()->face());
    coordinate_type N1 = asawa::ci::normal<SPACE>(e->v2()->face());
    coordinate_type dp = c1 - c0;

    coordinate_type M0 = dp.cross(N0);
    coordinate_type M1 = dp.cross(N1);

    gradU[e->v1()->face()->position_in_set()] -= 0.5 * M0 * u0 / A0;
    gradU[e->v2()->face()->position_in_set()] += 0.5 * M1 * u1 / A1;
  }
  return gradU;
}

template <typename SPACE>
std::vector<typename SPACE::real>
divergence(const std::vector<typename SPACE::coordinate_type> &gradU,
           typename surf<SPACE>::surf_ptr s) {
  M2_TYPEDEFS;
  ///////////////
  // divergence
  ///////////////

  edge_array &edges = s->get_edges();
  std::vector<real> divu(s->get_vertices().size(), 0.0);

  for (auto e : edges) {
    coordinate_type c0 = ci::get_coordinate<SPACE>(e->v1());
    coordinate_type c1 = ci::get_coordinate<SPACE>(e->v2());

    coordinate_type g0 = gradU[e->v1()->face()->position_in_set()];
    coordinate_type g1 = gradU[e->v2()->face()->position_in_set()];
    coordinate_type dp = c1 - c0;

    real cot0 = ci::cotan<SPACE>(e->v1()->prev());
    real cot1 = ci::cotan<SPACE>(e->v2()->prev());

    real l = 0.5 * (cot0 * dp.dot(g0) + cot1 * dp.dot(g1));

    divu[e->v1()->vertex()->position_in_set()] += l;
    divu[e->v2()->vertex()->position_in_set()] -= l;
  }

  return divu;
}
#endif

} // namespace asawa
#endif