#include <cassert>
#include <cmath>
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
  return va::abs_cotan(x0, xp, xn);
  // return va::cotan(x0, xp, xn);
}

vec3 edge_dir(manifold &M, index_t c0, const std::vector<vec3> &x) {
  index_t c1 = M.other(c0);
  return x[M.vert(c0)] - x[M.vert(M.next(c0))];
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

vec3 edge_center(manifold &M, index_t c0, const std::vector<vec3> &x) {
  index_t c1 = M.other(c0);
  vec3 x0 = x[M.vert(c0)];
  vec3 x1 = x[M.vert(c1)];
  return 0.5 * (x0 + x1);
}

vec3 face_center(manifold &M, index_t fi, const std::vector<vec3> &x) {
  vec3 c = vec3::Zero();
  int N = 0;
  M.for_each_face(fi, [&c, &x, &N](index_t c0, asawa::manifold &M) {
    c += x[M.vert(c0)];
    N++;
  });
  return c / real(N);
}

real face_area(manifold &M, index_t fi, const std::vector<vec3> &x) {
  real a = 0.0;
  M.for_each_face_tri(
      fi, [&a, &x](index_t c0, index_t c1, index_t c2, asawa::manifold &M) {
        vec3 x0 = x[M.vert(c0)];
        vec3 x1 = x[M.vert(c1)];
        vec3 x2 = x[M.vert(c2)];
        a += (x1 - x0).cross(x2 - x0).norm();
      });
  return a;
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

std::vector<real> face_areas(manifold &M, const std::vector<vec3> &x) {
  auto range = M.get_face_range();
  std::vector<real> A(range.size());
  int i = 0;
  for (auto vi : range) {
    A[i++] = face_area(M, i, x);
  }
  return A;
}
#if 1

std::vector<vec3> gradient(manifold &M, const std::vector<real> &u,
                           const std::vector<vec3> &x) {

  ///////////////
  // gradient
  ///////////////

  std::vector<index_t> edges = M.get_edge_range();

  std::vector<vec3> gradU(M.face_count(), vec3::Zero());

  for (int i = 0; i < edges.size(); i++) {
    index_t c0 = edges[i];
    index_t c1 = M.other(c0);
    vec3 x0 = x[M.vert(c0)];
    vec3 x1 = x[M.vert(c1)];
    real u0 = u[M.vert(M.prev(c0))];
    real u1 = u[M.vert(M.prev(c1))];

    real A0 = asawa::face_area(M, M.face(c0), x);
    real A1 = asawa::face_area(M, M.face(c1), x);
    real iA0 = A0 < 1e-6 ? 0.0 : 1.0 / A0;
    real iA1 = A1 < 1e-6 ? 0.0 : 1.0 / A1;

    vec3 N0 = asawa::face_normal(M, M.face(c0), x);
    vec3 N1 = asawa::face_normal(M, M.face(c1), x);

    vec3 dp0 = edge_dir(M, c0, x);
    vec3 dp1 = edge_dir(M, c1, x);

    // real sgn = va::sgn(N0, N1, dp);
    vec3 M0 = dp0.cross(N0);
    vec3 M1 = dp1.cross(N1);
#if 0
    if (M.face(c0) == 100) {
      vec3 e = edge_center(M, c0, x);
      std::cout << " A0: " << A0 << std::endl;
      gg::geometry_logger::line(e, e + M0, vec4(0.8, 1.0, 0.35, 1.0));
      gg::geometry_logger::line(e, e + N0, vec4(0.2, 1.0, 0.65, 1.0));
    }
    if (M.face(c1) == 100) {
      vec3 e = edge_center(M, c1, x);
      std::cout << " A1: " << A1 << std::endl;
      gg::geometry_logger::line(e, e + M1, vec4(0.8, 1.0, 0.35, 1.0));
      gg::geometry_logger::line(e, e + N1, vec4(0.2, 1.0, 0.65, 1.0));
    }
#endif
    gradU[M.face(c0)] += 0.5 * M0 * u0 * iA0;
    gradU[M.face(c1)] += 0.5 * M1 * u1 * iA1;
  }
  return gradU;
}

std::vector<real> divergence(manifold &M, const std::vector<vec3> &g,
                             const std::vector<vec3> &x) {

  ///////////////
  // divergence
  ///////////////

  std::vector<index_t> edges = M.get_edge_range();

  std::vector<real> divu(M.vert_count(), 0.0);

  for (int i = 0; i < edges.size(); i++) {
    index_t c0 = edges[i];
    index_t c1 = M.other(c0);
    vec3 v0 = x[M.vert(c0)];
    vec3 v1 = x[M.vert(c1)];
    vec3 g0 = g[M.face(c0)];
    vec3 g1 = g[M.face(c1)];

    vec3 dp = v1 - v0;
    vec3 dp0 = edge_dir(M, c0, x);
    real s = va::sgn(dp0.dot(dp));

    real cot0 = cotan(M, c0, x);
    real cot1 = cotan(M, c1, x);

    real l = 0.5 * (cot0 * dp.dot(g0) + cot1 * dp.dot(g1));
    assert(!isnan(l));
    divu[M.vert(c0)] -= s * l;
    divu[M.vert(c1)] += s * l;
  }

  return divu;
}
#endif

} // namespace asawa
#endif