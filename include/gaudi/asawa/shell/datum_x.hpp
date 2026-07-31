#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>

#include <iostream>
#include <limits>
#include <map>
#include <memory.h>
#include <ostream>
#include <set>
#include <stack>
#include <stdio.h>
#include <vector>

#include "gaudi/geometry_types.hpp"
#include "gaudi/vec_addendum.h"
#include "gaudi/logger.hpp"
#include "gaudi/common.h"
#include "shell.hpp"
#include "gaudi/geometry_logger.hpp"

#include <Eigen/Dense>
#include <Eigen/SVD>

#ifndef __ASAWA_X_DATUM__
#define __ASAWA_X_DATUM__
namespace gaudi {
namespace asawa {

void center(std::vector<vec3> &coords, real scale = 2.0) {
  real accum = 0.0;
  vec3 min = coords[0];
  vec3 max = coords[0];

  for (auto &c : coords) {
    min = va::min(c, min);
    max = va::max(c, max);
  }

  vec3 cen = 0.5 * (max + min);
  for (auto &c : coords) {
    c -= cen;
  }

  vec3 dl = (max - min);
  real maxl = dl[0];
  maxl = maxl > dl[1] ? maxl : dl[1];
  maxl = maxl > dl[2] ? maxl : dl[2];

  for (auto &c : coords) {
    c = scale * c / maxl;
  }

  std::cout << " scale: " << scale / maxl << std::endl;
}

std::array<vec3, 2> extents(const std::vector<vec3> &coords) {
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
real cotan(const shell &M, CornerId ci, const std::vector<vec3> &x) {

  vec3 xp = x[M.vert(M.prev(ci))];
  vec3 x0 = x[M.vert(ci)];
  vec3 xn = x[M.vert(M.next(ci))];

  return va::abs_cotan(x0, xp, xn);
  // return va::cotan(x0, xp, xn);
}

/// Signed half-angle cot weights for the **weak cotan Laplacian** (Pinkall–Polthier).
/// Degenerate corners contribute **0**; |cot| is capped (see `va::cotan_robust`).
/// `cotan()` above stays nonnegative for constraints / legacy weights.
real weak_cotan(const shell &M, CornerId ci, const std::vector<vec3> &x) {
  vec3 xp = x[M.vert(M.prev(ci))];
  vec3 x0 = x[M.vert(ci)];
  vec3 xn = x[M.vert(M.next(ci))];
  return va::cotan_robust(x0, xp, xn);
}

real angle(const shell &M, CornerId ci, const std::vector<vec3> &x) {

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

vec3 face_cross(const shell &M, FaceId fi, const std::vector<vec3> &x) {
  vec3 X = vec3::Zero();
  M.const_for_each_face_tri(
      fi, [&X, &x](CornerId c0, CornerId c1, CornerId c2, const shell &M) {
        vec3 x0 = x[M.vert(c0)];
        vec3 x1 = x[M.vert(c1)];
        vec3 x2 = x[M.vert(c2)];
        X += (x1 - x0).cross(x2 - x0);
      });
  return X;
}

vec3 face_normal(const shell &M, FaceId fi, const std::vector<vec3> &x) {
  vec3 N = face_cross(M, fi, x);
  return N.normalized();
}

real face_area(const shell &M, FaceId fi, const std::vector<vec3> &x) {
  real a = 0.0;
  vec3 N = face_cross(M, fi, x);
  return 0.5 * N.norm();
}

vec3 face_center(const shell &M, FaceId fi, const std::vector<vec3> &x) {
  vec3 c = vec3::Zero();
  int N = 0;
  M.const_for_each_face(fi, [&c, &x, &N](CornerId c0, const shell &M) {
    c += x[M.vert(c0)];
    N++;
  });
  return c / real(N);
}
vec3 face_interp(std::array<real, 3> s, const shell &M, FaceId fi,
                 const std::vector<vec3> &x) {
  CornerId c0 = M.fbegin(fi);
  CornerId c1 = M.next(c0);
  CornerId c2 = M.next(c1);
  vec3 x0 = x[M.vert(c0)];
  vec3 x1 = x[M.vert(c1)];
  vec3 x2 = x[M.vert(c2)];
  return s[0] * x0 + s[1] * x1 + s[2] * x2;
}
vec3 face_pnt(vec3 pt, const shell &M, FaceId fi, const std::vector<vec3> &x) {
  CornerId c0 = M.fbegin(fi);
  CornerId c1 = M.next(c0);
  CornerId c2 = M.next(c1);
  vec3 x0 = x[M.vert(c0)];
  vec3 x1 = x[M.vert(c1)];
  vec3 x2 = x[M.vert(c2)];
  std::array<real, 4> dist = va::closest_point({x0, x1, x2}, pt);
  return face_interp({dist[1], dist[2], dist[3]}, M, fi, x);
}

vec3 vert_normal(const shell &M, VertId vi, const std::vector<vec3> &x) {
  vec3 N = vec3::Zero();
  M.const_for_each_vertex(vi, [&N, &x](CornerId ci, const shell &M) {
    real ede = angle(M, ci, x);
    N += ede * face_cross(M, M.face(ci), x);
  });
  return N.normalized();
}

real vert_area(const shell &M, VertId vi, const std::vector<vec3> &x) {
  real A = 0.0;
  M.const_for_each_vertex(vi, [&A, &x](CornerId ci, const shell &M) {
    A += face_area(M, M.face(ci), x);
  });
  return A / 3.0;
}

/// Barycentric-dual mass for an undirected edge: half the sum of endpoint
/// \ref vert_area values (one third of adjacent face areas per vertex).
inline real edge_barycentric_dual_mass(const shell &M, CornerId c,
                                       const std::vector<vec3> &x) {
  return 0.5 * (vert_area(M, M.vert(c), x) + vert_area(M, M.vert(M.next(c)), x));
}

real vert_cotan_weight(const shell &M, VertId vi, const std::vector<vec3> &x) {
  real w = 0.0;
  M.const_for_each_vertex(vi, [&w, &x](CornerId ci, const shell &M) {
    CornerId c0p = M.prev(ci);
    CornerId c1p = M.prev(M.other(ci));
    w += cotan(M, c0p, x) + cotan(M, c1p, x);
  });
  return w;
}

std::vector<real> vert_cotan_weights(const shell &M, VertId vi,
                                     const std::vector<vec3> &x) {
  std::vector<real> w;
  M.const_for_each_vertex(vi, [&w, &x](CornerId ci, const shell &M) {
    CornerId c0p = M.prev(ci);
    CornerId c1p = M.prev(M.other(ci));
    w.push_back(cotan(M, c0p, x) + cotan(M, c1p, x));
  });
  return w;
}

std::vector<real> vert_angle_weights(const shell &M, VertId vi,
                                     const std::vector<vec3> &x) {
  std::vector<real> w;
  M.const_for_each_vertex(vi, [&w, &x](CornerId ci, const shell &M) {
    real thet = angle(M, ci, x);
    real A = face_area(M, M.face(ci), x);
    w.push_back(thet * A);
  });
  return w;
}

std::vector<real> vert_unitary_weights(const shell &M, VertId vi,
                                       const std::vector<vec3> &x) {
  std::vector<real> w;
  M.const_for_each_vertex(
      vi, [&w, &x](CornerId ci, const shell &M) { w.push_back(1.0); });
  return w;
}

vec3 edge_tangent(const shell &M, CornerId c0, const std::vector<vec3> &x) {
  return x[M.vert(c0)] - x[M.vert(M.next(c0))];
}

vec3 g_edge_tangent(const shell &M, CornerId c0, const std::vector<vec3> &x) {
  real s = c0 % 2 == 0 ? 1.0 : -1.0;
  vec3 tan = x[M.vert(c0)] - x[M.vert(M.next(c0))];
  return s * tan;
}

vec3 edge_normal(const shell &M, CornerId c0, const std::vector<vec3> &x) {
  CornerId c1 = M.other(c0);
  vec3 N = face_normal(M, M.face(c0), x);
  N += face_normal(M, M.face(c1), x);
  return N.normalized();
}

vec3 edge_vert(const shell &M, CornerId c0, const real &s,
               const std::vector<vec3> &x) {
  CornerId c1 = M.other(c0);
  vec3 x0 = x[M.vert(c0)];
  vec3 x1 = x[M.vert(c1)];
  return va::mix(s, x0, x1);
}

vec3 edge_center(const shell &M, CornerId c0, const std::vector<vec3> &x) {
  return edge_vert(M, c0, 0.5, x);
}

real edge_length(const shell &M, CornerId c0, const std::vector<vec3> &x) {
  CornerId c1 = M.other(c0);
  vec3 x0 = x[M.vert(c0)];
  vec3 x1 = x[M.vert(c1)];
  return (x1 - x0).norm();
}

std::vector<vec3> face_normals(const shell &M, const std::vector<vec3> &x) {
  auto range = M.get_face_range();
  std::vector<vec3> Ns(M.face_count(), vec3::Zero());
  for (auto vi : range) {
    Ns[vi] = face_normal(M, face_id(vi), x);
  }
  return Ns;
}

std::vector<vec3> vertex_normals(const shell &M, const std::vector<vec3> &x) {
  auto range = M.get_vert_range();
  std::vector<vec3> Ns(M.vert_count(), vec3::Zero());
  for (auto vi : range) {
    Ns[vi] = vert_normal(M, vert_id(vi), x);
  }
  return Ns;
}

std::vector<real> vertex_areas(const shell &M, const std::vector<vec3> &x) {
  auto range = M.get_vert_range();
  std::vector<real> Ns(M.vert_count(), 0.0);
  for (auto vi : range) {
    real area = vert_area(M, vert_id(vi), x);
    Ns[vi] = area;
  }
  return Ns;
}

std::vector<vec3> vertex_areas_3(const shell &M, const std::vector<vec3> &x) {
  auto range = M.get_vert_range();
  std::vector<vec3> Ns(M.vert_count(), vec3::Zero());
  for (auto vi : range) {
    real area = vert_area(M, vert_id(vi), x);
    Ns[vi] = vec3(area, area, area);
  }
  return Ns;
}

std::vector<real> edge_cotan_weights(const shell &M,
                                     const std::vector<vec3> &x) {
  auto range = M.get_edge_range();
  std::vector<real> ws(M.edge_count(), 0.0);
  for (auto ci : range) {
    CornerId cid = corner_id(ci);
    CornerId c0p = M.prev(cid);
    CornerId c1p = M.prev(M.other(cid));
    real ct = cotan(M, c0p, x) + cotan(M, c1p, x);
    ws[cid / 2] = ct;
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
  std::vector<real> l(M.edge_count());
  for (auto ci : range) {
    CornerId cid = corner_id(ci);
    l[cid / 2] = edge_length(M, cid, x);
  }
  return l;
}

std::vector<vec3> edge_centers(const shell &M, const std::vector<vec3> &x) {
  auto range = M.get_edge_range();
  std::vector<vec3> cens(M.edge_count(), vec3::Zero());
  for (auto ci : range) {
    CornerId cid = corner_id(ci);
    cens[cid / 2] = edge_center(M, cid, x);
  }
  return cens;
}

std::vector<vec3> edge_normals(const shell &M, const std::vector<vec3> &x) {
  auto range = M.get_edge_range();
  std::vector<vec3> Ns(M.edge_count(), vec3::Zero());
  for (auto ci : range) {
    CornerId cid = corner_id(ci);
    Ns[cid / 2] = edge_normal(M, cid, x);
  }
  return Ns;
}

std::vector<real> edge_areas(const shell &M, const std::vector<vec3> &x) {
  // well this is wrong...
  auto range = M.get_edge_range();
  std::vector<real> ws(M.edge_count(), 0.0);
  for (auto ci : range) {
    CornerId i0 = corner_id(ci);
    CornerId i1 = M.other(i0);
    FaceId f0 = M.face(i0);
    FaceId f1 = M.face(i1);
    ws[i0 / 2] =
        (face_area(M, f0, x) + face_area(M, f1, x)) / 3.0;
  }
  return ws;
}

std::vector<vec3> edge_tangents(const shell &M, const std::vector<vec3> &x) {
  auto range = M.get_edge_range();
  std::vector<vec3> dirs(M.edge_count(), vec3::Zero());
  for (auto ci : range) {
    CornerId cid = corner_id(ci);
    dirs[cid / 2] = edge_tangent(M, cid, x);
  }
  return dirs;
}

std::vector<vec3> face_centers(const shell &M, const std::vector<vec3> &x) {
  auto range = M.get_face_range();
  std::vector<vec3> xc(M.face_count(), vec3::Zero());
  for (auto fi : range) {
    xc[fi] = face_center(M, face_id(fi), x);
  }
  return xc;
}

std::vector<real> face_areas(const shell &M, const std::vector<vec3> &x) {
  auto range = M.get_face_range();
  std::vector<real> A(M.face_count(), 0.0);
  for (auto fi : range) {
    A[fi] = face_area(M, face_id(fi), x);
  }
  return A;
}

template <typename TYPE>
std::vector<TYPE> expand_from_vert_range(const shell &M,
                                         const std::vector<TYPE> &x) {
  auto v_range = M.get_vert_range();
  std::vector<TYPE> x_exp(M.vert_count());
  int i = 0;
  for (auto vi : v_range) {
    x_exp[vi] = x[i++];
  }
  return x_exp;
}

template <typename TYPE, typename ID>
std::vector<TYPE> compress_to_range(const std::vector<ID> &element_range,
                                    const std::vector<TYPE> &v) {
  assert(element_range.size() == v.size());
  std::vector<TYPE> v_comp(element_range.size());
  int i = 0;
  for (auto ei : element_range) {
    v_comp[i++] = v[ei];
  }
  return v_comp;
}

template <typename TYPE>
std::vector<TYPE> compress_to_vert_range(const shell &M,
                                         const std::vector<TYPE> &v) {
  return compress_to_range<TYPE>(M.get_vert_range(), v);
}

template <typename TYPE>
std::vector<TYPE> compress_to_face_range(const shell &M,
                                         const std::vector<TYPE> &v) {
  return compress_to_range<TYPE>(M.get_face_range(), v);
}

template <typename TYPE>
std::vector<TYPE> vert_to_face(const shell &M, const std::vector<vec3> &x,
                               const std::vector<TYPE> &v) {
  auto range = M.get_face_range();
  std::vector<TYPE> vals(M.face_count());
  for (auto fi : range) {
    TYPE c = z::zero<TYPE>();
    M.const_for_each_face(face_id(fi), [&c, &v](CornerId c0, const shell &M) {
      c += 0.33333 * v[M.vert(c0)];
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
    M.const_for_each_vertex(vert_id(vi), [&c, &x](CornerId c0, const shell &M) {
      c += 0.33333 * x[M.face(c0)];
    });
    vals[vi] = c;
  }
  return vals;
}

template <typename TYPE>
std::vector<TYPE> edge_to_face(const shell &M, const std::vector<vec3> &x,
                               const std::vector<TYPE> &v) {
  auto range = M.get_face_range();
  std::vector<TYPE> vals(M.face_count());
  for (auto fi : range) {
    TYPE c = z::zero<TYPE>();
    real A = 0.0;
    FaceId fii = face_id(fi);
    vec3 cen = face_center(M, fii, x);
    M.const_for_each_face(fii, [&](CornerId c0, const shell &M) {
      CornerId c1 = M.other(c0);
      vec3 x0 = x[M.vert(c0)];
      vec3 x1 = x[M.vert(c1)];
      vec3 d0c = x0 - cen;
      vec3 d1c = x1 - cen;
      real Ai = 0.5 * d0c.cross(d1c).norm();
      A += Ai;
      c += Ai * v[c0 / 2];
    });
    vals[fi] = c / A;
  }
  return vals;
}

real surface_area(const shell &M, const std::vector<vec3> &x) {
  auto range = M.get_face_range();
  real A = 0.0;
  for (auto vi : range) {
    A += face_area(M, face_id(vi), x);
  }
  return A;
}

/*
real avg_length(const shell &M, const std::vector<vec3> &x) {
  real accum = 0.0;
  auto range = M.get_edge_range();
  int count = 0;
  for (auto ci : range) {
    accum += edge_length(M, ci, x);
    count++;
  }

  return accum / real(count);
}
*/
real avg_length(const shell &M, const std::vector<vec3> &coords) {
  real accum = 0.0;
  for (int i = 0; i < static_cast<int>(M.__corners_next.size()); i += 2) {
    if (M.__corners_next[static_cast<size_t>(i)] < 0)
      continue;
    CornerId i0 = corner_id(i);
    CornerId i1 = M.other(i0);
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

std::vector<vec3> circulation(shell &M, const std::vector<real> &u,
                              const std::vector<vec3> &x) {

  ///////////////
  // circulation
  ///////////////

  auto edges = M.get_edge_range();

  std::vector<vec3> circU(M.face_count(), vec3::Zero());

  for (int i = 0; i < static_cast<int>(edges.size()); i++) {
    CornerId c0 = corner_id(edges[static_cast<size_t>(i)]);
    CornerId c1 = M.other(c0);
    vec3 x0 = x[M.vert(c0)];
    vec3 x1 = x[M.vert(c1)];

    real A0 = face_area(M, M.face(c0), x);
    real A1 = face_area(M, M.face(c1), x);
    real iA0 = A0 < 1e-6 ? 0.0 : 1.0 / A0;
    real iA1 = A1 < 1e-6 ? 0.0 : 1.0 / A1;
    vec3 e0 = edge_tangent(M, c0, x);
    vec3 e1 = edge_tangent(M, c1, x);

#if 0
    real u0 = u[M.vert(c0)];
    real u1 = u[M.vert(c1)];
    real ui = 0.5 * (u0 + u1);
    circU[M.face(c0)] += 0.5 * e0 * ui * iA0;
    circU[M.face(c1)] += 0.5 * e1 * ui * iA1;
#else
    real u0 = u[M.vert(M.prev(c0))];
    real u1 = u[M.vert(M.prev(c1))];
    circU[M.face(c0)] += 0.5 * e0 * u0 * iA0;
    circU[M.face(c1)] += 0.5 * e1 * u1 * iA1;
#endif
  }

  return circU;
}

std::vector<vec3> gradient(shell &M, const std::vector<real> &u,
                           const std::vector<vec3> &x) {

  ///////////////
  // gradient
  ///////////////

  auto edges = M.get_edge_range();

  std::vector<vec3> gradU(M.face_count(), vec3::Zero());

  for (int i = 0; i < static_cast<int>(edges.size()); i++) {
    CornerId c0 = corner_id(edges[static_cast<size_t>(i)]);
    CornerId c1 = M.other(c0);
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

    vec3 M0 = dp0.cross(N0);
    vec3 M1 = dp1.cross(N1);
#if 0
    if (M.face(c0) == 100) {
      vec3 e = edge_center(M, c0, x);
      std::cout << " A0: " << A0 << std::endl;
      geometry_logger::line(e, e + M0, vec4(0.8, 1.0, 0.35, 1.0));
      geometry_logger::line(e, e + N0, vec4(0.2, 1.0, 0.65, 1.0));
    }
    if (M.face(c1) == 100) {
      vec3 e = edge_center(M, c1, x);
      std::cout << " A1: " << A1 << std::endl;
      geometry_logger::line(e, e + M1, vec4(0.8, 1.0, 0.35, 1.0));
      geometry_logger::line(e, e + N1, vec4(0.2, 1.0, 0.65, 1.0));
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

  auto edges = M.get_edge_range();

  std::vector<real> divu(M.vert_count(), 0.0);

  for (int i = 0; i < static_cast<int>(edges.size()); i++) {
    CornerId c0 = corner_id(edges[static_cast<size_t>(i)]);
    CornerId c1 = M.other(c0);
    CornerId c0p = M.prev(c0);
    CornerId c1p = M.prev(c1);
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

/// Stencil for per-face curvature fitting (quadric height field over neighbor face centers).
enum class face_curvature_stencil { one_ring, butterfly, two_ring };

struct face_curvature_frame {
  vec3 n = vec3::UnitZ();
  vec3 t_min = vec3::UnitX();
  vec3 t_max = vec3::UnitY();
  real k_min = 0;
  real k_max = 0;
};

/// Monge / height jet in foot frame: z ≈ a x² + b x y + c y² over (u,v), axis n.
struct face_height_jet {
  vec3 foot = vec3::Zero();
  vec3 u = vec3::UnitX();
  vec3 v = vec3::UnitY();
  vec3 n = vec3::UnitZ();
  real a = 0.0;
  real b = 0.0;
  real c = 0.0;
  bool valid = false;

  /// ∇F for F = dp·n − a s² − b s t − c t², with s=dp·u, t=dp·v (foot-relative dp).
  vec3 grad_at_rel(const vec3 &dp) const {
    const real s = dp.dot(u);
    const real t = dp.dot(v);
    return n - (2.0 * a * s + b * t) * u - (b * s + 2.0 * c * t) * v;
  }
};

inline void face_tangent_basis_from_normal(const vec3 &n_in, vec3 *u_out,
                                           vec3 *v_out) {
  vec3 n = n_in.normalized();
  vec3 a = std::abs(n[0]) > 0.5 ? vec3(0, 1, 0) : vec3(1, 0, 0);
  *u_out = n.cross(a);
  if (u_out->norm() < 1e-12) {
    *u_out = vec3(0, 0, 1).cross(n);
  }
  *u_out = u_out->normalized();
  *v_out = n.cross(*u_out).normalized();
}

inline std::vector<FaceId>
face_curvature_stencil_faces(const shell &M, FaceId f,
                             face_curvature_stencil stencil) {
  if (stencil == face_curvature_stencil::one_ring)
    return M.face_one_ring_face_ids(f);
  if (stencil == face_curvature_stencil::two_ring)
    return M.face_two_ring_face_ids(f);
  CornerId c0 = M.fbegin(f);
  return M.butterfly_face_ids(c0);
}

/// Fit Monge height jet at face center; also fills principal frame when possible.
inline face_height_jet face_height_jet_fit(
    const shell &M, const std::vector<vec3> &x, FaceId f,
    face_curvature_stencil stencil = face_curvature_stencil::one_ring) {
  face_height_jet jet;
  if (M.fbegin(f) < 0 || M.fsize(f) != 3)
    return jet;

  const std::vector<FaceId> stencil_faces =
      face_curvature_stencil_faces(M, f, stencil);
  jet.foot = face_center(M, f, x);
  jet.n = face_normal(M, f, x);
  face_tangent_basis_from_normal(jet.n, &jet.u, &jet.v);

  const int m = static_cast<int>(stencil_faces.size());
  if (m < 3) {
    return jet;
  }

  Eigen::MatrixXd A(m, 3);
  Eigen::VectorXd bz(m);
  for (int i = 0; i < m; ++i) {
    FaceId fi = stencil_faces[static_cast<size_t>(i)];
    vec3 p = face_center(M, fi, x);
    vec3 d = p - jet.foot;
    real xi = d.dot(jet.u);
    real yi = d.dot(jet.v);
    real zi = d.dot(jet.n);
    A(i, 0) = xi * xi;
    A(i, 1) = xi * yi;
    A(i, 2) = yi * yi;
    bz(i) = zi;
  }

  Eigen::JacobiSVD<Eigen::MatrixXd> svd(
      A, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::Vector3d coef = svd.solve(bz);
  jet.a = coef(0);
  jet.b = coef(1);
  jet.c = coef(2);
  jet.valid = std::isfinite(jet.a) && std::isfinite(jet.b) &&
              std::isfinite(jet.c);
  return jet;
}

/// Height jet at a vertex (foot = vertex, n = given normal).
inline face_height_jet vertex_height_jet_fit(
    const shell &M, const std::vector<vec3> &x, VertId vi, const vec3 &n_in,
    face_curvature_stencil stencil = face_curvature_stencil::two_ring) {
  face_height_jet jet;
  const int vid = static_cast<int>(vi);
  if (vid < 0 || M.vbegin(vi) < 0)
    return jet;

  std::set<FaceId> seed;
  M.const_for_each_vertex(vi, [&](CornerId c, const shell &Mm) {
    FaceId fi = Mm.face(c);
    if (static_cast<int>(fi) >= 0)
      seed.insert(fi);
  });
  if (seed.empty())
    return jet;

  std::set<FaceId> stencil_faces = seed;
  if (stencil == face_curvature_stencil::two_ring ||
      stencil == face_curvature_stencil::butterfly) {
    for (FaceId f0 : seed) {
      for (FaceId f1 : face_curvature_stencil_faces(M, f0, stencil)) {
        if (static_cast<int>(f1) >= 0)
          stencil_faces.insert(f1);
      }
    }
  } else {
    for (FaceId f0 : seed) {
      for (FaceId f1 : M.face_one_ring_face_ids(f0)) {
        if (static_cast<int>(f1) >= 0)
          stencil_faces.insert(f1);
      }
    }
  }

  jet.foot = x[static_cast<size_t>(vid)];
  jet.n = (n_in.norm() > 1e-12) ? n_in.normalized() : vec3::UnitZ();
  face_tangent_basis_from_normal(jet.n, &jet.u, &jet.v);

  const int m = static_cast<int>(stencil_faces.size());
  if (m < 3)
    return jet;

  Eigen::MatrixXd A(m, 3);
  Eigen::VectorXd bz(m);
  int row = 0;
  for (FaceId fi : stencil_faces) {
    vec3 p = face_center(M, fi, x);
    vec3 d = p - jet.foot;
    real xi = d.dot(jet.u);
    real yi = d.dot(jet.v);
    real zi = d.dot(jet.n);
    A(row, 0) = xi * xi;
    A(row, 1) = xi * yi;
    A(row, 2) = yi * yi;
    bz(row) = zi;
    ++row;
  }

  Eigen::JacobiSVD<Eigen::MatrixXd> svd(
      A, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::Vector3d coef = svd.solve(bz);
  jet.a = coef(0);
  jet.b = coef(1);
  jet.c = coef(2);
  jet.valid = std::isfinite(jet.a) && std::isfinite(jet.b) &&
              std::isfinite(jet.c);
  return jet;
}

/// Discrete principal frame via least-squares quadratic height field
/// \(z \approx a x^2 + b x y + c y^2\) on face barycenters in the stencil
/// (Rusinkiewicz-style jet on a triangle mesh; see Rusinkiewicz, 3DPVT 2004).
inline face_curvature_frame face_curvature_frame_fit(
    const shell &M, const std::vector<vec3> &x, FaceId f,
    face_curvature_stencil stencil = face_curvature_stencil::one_ring) {
  face_curvature_frame out;
  face_height_jet jet = face_height_jet_fit(M, x, f, stencil);
  out.n = jet.n;
  out.t_min = jet.u;
  out.t_max = jet.v;
  if (!jet.valid)
    return out;

  Eigen::Matrix2d H;
  H(0, 0) = 2.0 * jet.a;
  H(0, 1) = jet.b;
  H(1, 0) = jet.b;
  H(1, 1) = 2.0 * jet.c;

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es(H);
  if (es.info() != Eigen::Success) {
    return out;
  }

  Eigen::Vector2d ev0 = es.eigenvectors().col(0);
  Eigen::Vector2d ev1 = es.eigenvectors().col(1);
  out.k_min = es.eigenvalues()[0];
  out.k_max = es.eigenvalues()[1];
  out.t_min = (ev0[0] * jet.u + ev0[1] * jet.v).normalized();
  out.t_max = (ev1[0] * jet.u + ev1[1] * jet.v).normalized();
  return out;
}

} // namespace shell
} // namespace asawa
} // namespace gaudi
#endif