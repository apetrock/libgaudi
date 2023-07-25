#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cxxabi.h>

#include <execinfo.h>
#include <iostream>
#include <map>
#include <memory.h>
#include <ostream>
#include <set>
#include <stack>
#include <stdio.h>
#include <vector>
#include <zlib.h>

#include "gaudi/geometry_types.hpp"
#include "gaudi/vec_addendum.h"

#include "gaudi/common.h"
#include "shell.hpp"

#ifndef __ASAWA_X_DATUM__
#define __ASAWA_X_DATUM__
namespace gaudi {
namespace asawa {

void center(std::vector<vec3> &coords) {
  real accum = 0.0;
  vec3 min = coords[0];
  vec3 max = coords[0];

  for (auto &c : coords) {
    min = va::min(c, min);
    max = va::max(c, max);
  }

  vec3 dl = (max - min);
  real maxl = dl[0];
  maxl = maxl > dl[1] ? maxl : dl[1];
  maxl = maxl > dl[2] ? maxl : dl[2];
  real s = 2.0 / maxl;
  vec3 cen = (0.5 * (max + min) - min) * s;
  std::cout << " scale: " << s << std::endl;
  cen -= min;
  for (auto &c : coords) {
    c -= min;
    c = s * c;
    c -= cen;
  }
}

std::array<vec3, 2> extents(std::vector<vec3> &coords) {
  real accum = 0.0;
  vec3 min = coords[0];
  vec3 max = coords[0];

  for (auto &c : coords) {
    min = va::min(c, min);
    max = va::max(c, max);
  }
  return {min, max};
}

namespace shell {
/*TODO: these could all be namespaced...*/
real cotan(const shell &M, index_t ci, const std::vector<vec3> &x) {

  vec3 xp = x[M.vert(M.prev(ci))];
  vec3 x0 = x[M.vert(ci)];
  vec3 xn = x[M.vert(M.next(ci))];

  return va::abs_cotan(x0, xp, xn);
  // return va::cotan(x0, xp, xn);
}

real angle(const shell &M, index_t ci, const std::vector<vec3> &x) {

  vec3 xp = x[M.vert(M.prev(ci))];
  vec3 x0 = x[M.vert(ci)];
  vec3 xn = x[M.vert(M.next(ci))];
  vec3 e0 = (xp - x0).normalized();
  vec3 e1 = (xn - x0).normalized();
  real ede = e0.dot(e1);
  ede = va::clamp(ede, 0.0, 1.0 - 1e-8);
  ede = acos(ede);

  return ede;
}

vec3 face_cross(const shell &M, index_t fi, const std::vector<vec3> &x) {
  vec3 X = vec3::Zero();
  M.const_for_each_face_tri(
      fi, [&X, &x](index_t c0, index_t c1, index_t c2, const shell &M) {
        vec3 x0 = x[M.vert(c0)];
        vec3 x1 = x[M.vert(c1)];
        vec3 x2 = x[M.vert(c2)];
        X += (x1 - x0).cross(x2 - x0);
      });
  return X;
}

vec3 face_normal(const shell &M, index_t fi, const std::vector<vec3> &x) {
  vec3 N = face_cross(M, fi, x);
  return N.normalized();
}

real face_area(const shell &M, index_t fi, const std::vector<vec3> &x) {
  real a = 0.0;
  vec3 N = face_cross(M, fi, x);
  return N.norm();
}

vec3 face_center(const shell &M, index_t fi, const std::vector<vec3> &x) {
  vec3 c = vec3::Zero();
  int N = 0;
  M.const_for_each_face(fi, [&c, &x, &N](index_t c0, const shell &M) {
    c += x[M.vert(c0)];
    N++;
  });
  return c / real(N);
}
vec3 face_interp(std::array<real, 3> s, const shell &M, index_t fi,
                 const std::vector<vec3> &x) {
  index_t c0 = M.fbegin(fi);
  index_t c1 = M.next(c0);
  index_t c2 = M.next(c1);
  vec3 x0 = x[M.vert(c0)];
  vec3 x1 = x[M.vert(c1)];
  vec3 x2 = x[M.vert(c2)];
  return s[0] * x0 + s[1] * x1 + s[2] * x2;
}
vec3 face_pnt(vec3 pt, const shell &M, index_t fi, const std::vector<vec3> &x) {
  index_t c0 = M.fbegin(fi);
  index_t c1 = M.next(c0);
  index_t c2 = M.next(c1);
  vec3 x0 = x[M.vert(c0)];
  vec3 x1 = x[M.vert(c1)];
  vec3 x2 = x[M.vert(c2)];
  std::array<real, 4> dist = va::closest_point({x0, x1, x2}, pt);
  return face_interp({dist[1], dist[2], dist[3]}, M, fi, x);
}

vec3 vert_normal(const shell &M, index_t vi, const std::vector<vec3> &x) {
  vec3 N = vec3::Zero();
  M.const_for_each_vertex(vi, [&N, &x](index_t ci, const shell &M) {
    real ede = angle(M, ci, x);
    N += ede * face_cross(M, M.face(ci), x);
  });
  return N.normalized();
}

real vert_area(const shell &M, index_t vi, const std::vector<vec3> &x) {
  real A = 0.0;
  M.const_for_each_vertex(vi, [&A, &x](index_t ci, const shell &M) {
    A += face_area(M, M.face(ci), x);
  });
  return A / 3.0;
}

real vert_cotan_weight(const shell &M, index_t vi, const std::vector<vec3> &x) {
  real w = 0.0;
  M.const_for_each_vertex(vi, [&w, &x](index_t ci, const shell &M) {
    index_t c0p = M.prev(ci);
    index_t c1p = M.prev(M.other(ci));
    w += cotan(M, c0p, x) + cotan(M, c1p, x);
  });
  return w;
}

std::vector<real> vert_cotan_weights(const shell &M, index_t vi,
                                     const std::vector<vec3> &x) {
  std::vector<real> w;
  M.const_for_each_vertex(vi, [&w, &x](index_t ci, const shell &M) {
    index_t c0p = M.prev(ci);
    index_t c1p = M.prev(M.other(ci));
    w.push_back(cotan(M, c0p, x) + cotan(M, c1p, x));
  });
  return w;
}

std::vector<real> vert_angle_weights(const shell &M, index_t vi,
                                     const std::vector<vec3> &x) {
  std::vector<real> w;
  M.const_for_each_vertex(vi, [&w, &x](index_t ci, const shell &M) {
    real thet = angle(M, ci, x);
    real A = face_area(M, M.face(ci), x);
    w.push_back(thet * A);
  });
  return w;
}

std::vector<real> vert_unitary_weights(const shell &M, index_t vi,
                                       const std::vector<vec3> &x) {
  std::vector<real> w;
  M.const_for_each_vertex(
      vi, [&w, &x](index_t ci, const shell &M) { w.push_back(1.0); });
  return w;
}

vec3 edge_tangent(const shell &M, index_t c0, const std::vector<vec3> &x) {
  index_t c1 = M.other(c0);
  return x[M.vert(c0)] - x[M.vert(M.next(c0))];
}

vec3 edge_normal(const shell &M, index_t c0, const std::vector<vec3> &x) {
  index_t c1 = M.other(c0);
  vec3 N = face_normal(M, M.face(c0), x);
  N += face_normal(M, M.face(c1), x);
  return N.normalized();
}

vec3 edge_vert(const shell &M, index_t c0, const real &s,
               const std::vector<vec3> &x) {
  index_t c1 = M.other(c0);
  vec3 x0 = x[M.vert(c0)];
  vec3 x1 = x[M.vert(c1)];
  return va::mix(s, x0, x1);
}

vec3 edge_center(const shell &M, index_t c0, const std::vector<vec3> &x) {
  return edge_vert(M, c0, 0.5, x);
}

real edge_length(const shell &M, index_t c0, const std::vector<vec3> &x) {
  index_t c1 = M.other(c0);
  std::cout << c0 << " " << c1 << " " << M.vert(c0) << " " << M.vert(c1) << " "
            << x.size() << std::endl;
  vec3 x0 = x[M.vert(c0)];
  vec3 x1 = x[M.vert(c1)];
  return (x1 - x0).norm();
}

std::vector<vec3> face_normals(const shell &M, const std::vector<vec3> &x) {
  auto range = M.get_face_range();
  std::vector<vec3> Ns(M.face_count(), vec3::Zero());
  int i = 0;
  for (auto vi : range) {
    Ns[vi] = face_normal(M, vi, x);
  }
  return Ns;
}

std::vector<vec3> vertex_normals(const shell &M, const std::vector<vec3> &x) {
  auto range = M.get_vert_range();
  std::vector<vec3> Ns(M.vert_count(), vec3::Zero());
  int i = 0;
  for (auto vi : range) {
    Ns[vi] = vert_normal(M, vi, x);
  }
  return Ns;
}

std::vector<vec3> vertex_areas(const shell &M, const std::vector<vec3> &x) {
  auto range = M.get_vert_range();
  std::vector<vec3> Ns(M.vert_count(), vec3::Zero());
  int i = 0;
  for (auto vi : range) {
    real area = vert_area(M, vi, x);
    Ns[vi] = vec3(area, area, area);
  }
  return Ns;
}

std::vector<real> edge_cotan_weights(const shell &M,
                                     const std::vector<vec3> &x) {
  auto range = M.get_edge_range();
  std::vector<real> ws(M.corner_count() / 2, 0.0);
  int i = 0;
  for (auto ci : range) {
    index_t c0p = M.prev(ci);
    index_t c1p = M.prev(M.other(ci));
    real ct = cotan(M, c0p, x) + cotan(M, c1p, x);
    ws[ci / 2] = ct;
  }
  return ws;
}
/*
std::vector<real> align_edges(shell &M, const std::vector<vec3> &x) {
  auto range = M.get_edge_range();
  std::vector<real> ws(range.size());
  int i = 0;
  for (auto ci : range) {
    vec3 dx = edge_tangent(M, ci, x);

    vec3 xf0 = x[M.vert(ci)];
    vec3 xf1 = x[M.vert(M.next(ci))];
    if (dx.dot(xf0 - xf1) < 0)
      M.flip_edge(ci);
  }
  return ws;
}
*/

std::vector<real> edge_lengths(const shell &M, const std::vector<vec3> &x) {
  auto range = M.get_edge_range();
  std::vector<real> l(range.size());
  int i = 0;
  for (auto ci : range) {
    l[i++] = edge_length(M, ci, x);
  }
  return l;
}

std::vector<vec3> edge_centers(const shell &M, const std::vector<vec3> &x) {
  auto range = M.get_edge_range();
  std::vector<vec3> cens(range.size());
  int i = 0;
  for (auto ci : range) {
    cens[i++] = edge_center(M, ci, x);
  }
  return cens;
}

std::vector<vec3> edge_normals(const shell &M, const std::vector<vec3> &x) {
  auto range = M.get_edge_range();
  std::vector<vec3> Ns(range.size());
  int i = 0;
  for (auto ci : range) {
    Ns[i++] = edge_normal(M, ci, x);
  }
  return Ns;
}

std::vector<real> edge_areas(const shell &M, const std::vector<vec3> &x) {
  auto range = M.get_edge_range();
  std::vector<real> ws(range.size());
  int i = 0;
  for (auto ci : range) {
    int i0 = ci;
    int i1 = M.other(i0);
    int f0 = M.face(i0);
    int f1 = M.face(i1);
    ws[i++] = (face_area(M, f0, x) + face_area(M, f1, x)) / 3.0;
  }
  return ws;
}

std::vector<vec3> edge_tangents(const shell &M, const std::vector<vec3> &x) {
  auto range = M.get_edge_range();
  std::vector<vec3> dirs(range.size());
  int i = 0;
  for (auto ci : range) {
    dirs[i++] = edge_tangent(M, ci, x);
  }
  return dirs;
}

std::vector<vec3> face_centers(const shell &M, const std::vector<vec3> &x) {
  auto range = M.get_face_range();
  std::vector<vec3> xc(M.face_count(), vec3::Zero());
  int i = 0;
  for (auto fi : range) {
    xc[fi] = face_center(M, fi, x);
  }
  return xc;
}

std::vector<real> face_areas(const shell &M, const std::vector<vec3> &x) {
  auto range = M.get_face_range();
  std::vector<real> A(M.face_count(), 0.0);
  int i = 0;
  for (auto fi : range) {
    A[i++] = face_area(M, fi, x);
  }
  return A;
}

template <typename TYPE>
std::vector<TYPE> vert_to_face(const shell &M, const std::vector<TYPE> &x) {
  auto range = M.get_face_range();
  std::vector<TYPE> vals(range.size());
  for (auto fi : range) {
    TYPE c = z::zero<TYPE>();
    M.const_for_each_face(fi, [&c, &x](index_t c0, const shell &M) {
      c += 0.33333 * x[M.vert(c0)];
    });
    vals[fi] = c;
  }
  return vals;
}

template <typename TYPE>
std::vector<TYPE> face_to_vert(const shell &M, const std::vector<TYPE> &x) {
  auto range = M.get_vert_range();
  std::vector<TYPE> vals(range.size());
  int i = 0;
  for (auto vi : range) {
    TYPE c = z::zero<TYPE>();
    M.const_for_each_vertex(vi, [&c, &x](index_t c0, const shell &M) {
      c += 0.33333 * x[M.face(c0)];
    });
    vals[vi] = c;
  }
  return vals;
}

real surface_area(const shell &M, const std::vector<vec3> &x) {
  auto range = M.get_face_range();
  real A;
  int i = 0;
  for (auto vi : range) {
    A += face_area(M, i, x);
  }
  return A;
}

real avg_length(const shell &M, const std::vector<vec3> &coords) {
  real accum = 0.0;
  for (int i = 0; i < M.__corners_next.size(); i += 2) {
    if (M.__corners_next[i] < 0)
      continue;
    int i0 = i;
    int i1 = M.other(i0);
    int v0 = M.vert(i0);
    int v1 = M.vert(i1);
    accum += (coords[M.vert(i0)] - coords[M.vert(i1)]).norm();
  }
  return 0.5 * accum / real(M.corner_count());
}

ext::extents_t ext(const shell &M, const std::vector<vec3> &x) {
  auto range = M.get_vert_range();
  double inf = std::numeric_limits<double>::max();
  ext::extents_t ext = {vec3(inf, inf, inf), vec3(-inf, -inf, -inf)};
  int i = 0;
  for (auto vi : range) {
    ext = ext::expand(ext, x[vi]);
  }
  return ext;
}

std::vector<vec3> gradient(shell &M, const std::vector<real> &u,
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

    real A0 = face_area(M, M.face(c0), x);
    real A1 = face_area(M, M.face(c1), x);
    real iA0 = A0 < 1e-6 ? 0.0 : 1.0 / A0;
    real iA1 = A1 < 1e-6 ? 0.0 : 1.0 / A1;

    vec3 N0 = face_normal(M, M.face(c0), x);
    vec3 N1 = face_normal(M, M.face(c1), x);

    vec3 dp0 = edge_tangent(M, c0, x);
    vec3 dp1 = edge_tangent(M, c1, x);

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

std::vector<real> divergence(shell &M, const std::vector<vec3> &g,
                             const std::vector<vec3> &x) {

  ///////////////
  // divergence
  ///////////////

  std::vector<index_t> edges = M.get_edge_range();

  std::vector<real> divu(M.vert_count(), 0.0);

  for (int i = 0; i < edges.size(); i++) {
    index_t c0 = edges[i];
    index_t c1 = M.other(c0);
    index_t c0p = M.prev(c0);
    index_t c1p = M.prev(c1);
    vec3 v0 = x[M.vert(c0)];
    vec3 v1 = x[M.vert(c1)];
    vec3 g0 = g[M.face(c0)];
    vec3 g1 = g[M.face(c1)];

    vec3 dp = v1 - v0;
    vec3 dp0 = edge_tangent(M, c0, x);
    real s = va::sgn(dp0.dot(dp));

    real cot0 = cotan(M, c0p, x);
    real cot1 = cotan(M, c1p, x);

    real l = 0.5 * (cot0 * dp.dot(g0) + cot1 * dp.dot(g1));
    assert(!isnan(l));
    divu[M.vert(c0)] -= s * l;
    divu[M.vert(c1)] += s * l;
  }

  return divu;
}

} // namespace shell
} // namespace asawa
} // namespace gaudi
#endif