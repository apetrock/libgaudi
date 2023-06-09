#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cxxabi.h>

#include <deque>
#include <execinfo.h>
#include <functional>
#include <iostream>
#include <iterator>
#include <map>

#include <list>

#include <math.h>
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

#include "datum_x.hpp"
#ifndef __ASAWA_WALK__
#define __ASAWA_WALK__
namespace gaudi {
namespace asawa {
namespace shell {

bool LineLineIntersect(vec2 p0, vec2 p1, vec2 q0, vec2 q1, real &intersection,
                       bool considerCollinearOverlapAsIntersect = false) {
  intersection = -1.0;

  auto cross = [](vec2 x, vec2 y) { return x[0] * y[1] - x[1] * y[0]; };
  auto eps = [](real x) { return std::abs(x) < 1e-8; };
  vec2 r = p1 - p0;
  vec2 s = q1 - q0;
  vec2 qp = q0 - p0;
  real rxs = cross(r, s);
  real qpxr = cross(qp, r);

  // If r x s = 0 and (q - p) x r = 0, then the two lines are collinear.
  if (eps(rxs) && eps(qpxr)) {
    // 1. If either  0 <= (q - p) * r <= r * r or 0 <= (p - q) * s <= * s
    // then the two lines are overlapping,
    if (considerCollinearOverlapAsIntersect)
      if ((0 <= qp.dot(r) && qp.dot(r) <= r.dot(r)) ||
          (0 <= -qp.dot(s) && -qp.dot(s) <= s.dot(s)))
        return true;

    // 2. If neither 0 <= (q - p) * r = r * r nor 0 <= (p - q) * s <= s * s
    // then the two lines are collinear but disjoint.
    // No need to implement this expression, as it follows from the expression
    // above.
    return false;
  }

  // 3. If r x s = 0 and (q - p) x r != 0, then the two lines are parallel and
  // non-intersecting.
  // td::cout << "rxs: " << rxs << std::endl;

  if (eps(rxs) && !eps(qpxr))
    return false;

  // t = (q - p) x s / (r x s)
  real t = cross(qp, s) / rxs;

  // u = (q - p) x r / (r x s)

  real u = cross(qp, r) / rxs;

  // 4. If r x s != 0 and 0 <= t <= 1 and 0 <= u <= 1
  // the two line segments meet at the point p + t r = q + u s.

  if (!rxs < 1e-12 && (0 <= t && t <= 1) && (0 <= u && u <= 1)) {
    // We can calculate the intersection point using either t or u.
    intersection = u;
    // An intersection was found.
    return true;
  }

  // 5. Otherwise, the two line segments are not parallel but do not
  // intersect.
  return false;
}

auto proj = [](const vec3 &x, const vec3 &xc, const vec3 &N0, const vec3 &N1) {
  return vec2((x - xc).dot(N0), (x - xc).dot(N1));
};

auto unproj = [](const vec2 &x, const vec3 &xc, const vec3 &N0,
                 const vec3 &N1) { return xc + x[0] * N0 + x[1] * N1; };

std::vector<vec2> project2(const std::vector<vec3> x, const vec3 &xc,
                           const vec3 T, const vec3 B, const vec3 N) {
  // B *= sg;
  // gg::geometry_logger::line(xc, xc + 0.02 * T, vec4(0.5, 0.0, 0.0, 1.0));
  // gg::geometry_logger::line(xc, xc + 0.02 * B, vec4(0.0, 0.5, 0.0, 1.0));
  // gg::geometry_logger::line(xc, xc + 0.02 * N, vec4(0.0, 0.0, 0.5, 1.0));
  std::vector<vec2> out(x.size());
  for (int i = 0; i < x.size(); i++) {
    out[i] = proj(x[i], xc, T, B);
  }
  return out;
}

bool has(const std::set<index_t> vis, index_t c) {
  return vis.find(c) != vis.end();
}
std::vector<vec2> project_set(const shell &M, const std::vector<vec3> &x,
                              const std::vector<index_t> &C, const vec3 dir,
                              const index_t &c0) {

  vec3 xc = x[M.vert(c0)];
  std::vector<vec3> x_col(C.size());
  for (int k = 0; k < C.size(); k++) {
    x_col[k] = x[M.vert(C[k])];
  }

  mat3 D = mat3::Zero();
  std::set<index_t> visited;
  M.const_for_each_vertex(M.vert(c0), [&](index_t ci, const shell &M) {
    index_t cn0 = M.next(ci);
    M.const_for_each_vertex(M.vert(cn0), [&](index_t ci, const shell &M) {
      index_t cn1 = M.next(ci);
      if (!has(visited, cn1)) {
        visited.insert(cn1);
        vec3 xn = x[M.vert(cn1)];
        vec3 dx = xn - xc;
        // gg::geometry_logger::line(xc, xc + 0.01 * dx,
        // vec4(0.0, 1.0, 1.0, 1.0));
        D += dx * dx.transpose();
      }
    });
  });

  vec3 x0 = x[M.vert(c0)];
  vec3 x1 = x0 + dir;
  x_col.push_back(x0);
  x_col.push_back(x1);

  Eigen::JacobiSVD<mat3> svd(D, Eigen::ComputeFullU);
  mat3 U = svd.matrixU();
  vec3 s = svd.singularValues();

  vec3 T = U.col(0);
  vec3 B = U.col(1);
  vec3 N = U.col(2);

  vec3 xN = x0 + 0.005 * N;
  // gg::geometry_logger::line(xN, xN + 0.01 * T, vec4(0.5, 0.0, 0.0, 1.0));
  // gg::geometry_logger::line(xN, xN + 0.01 * B, vec4(0.0, 0.5, 0.0, 1.0));
  // gg::geometry_logger::line(xN, xN + 0.01 * N, vec4(0.0, 0.0, 0.5, 1.0));
  std::vector<vec2> x_col2 = project2(x_col, xc, T, B, N);
  int Ni = x_col2.size() - 2;
  for (int i = 0; i < Ni; i++) {
    int i0 = i;
    int i1 = (i + 1) % Ni;
    vec3 x20 = unproj(x_col2[i0], xc, B, T);
    vec3 x21 = unproj(x_col2[i1], xc, B, T);
    // gg::geometry_logger::line(x20 + 0.01 * N, x21 + 0.01 * N,
    //                          vec4(0.5, 0.5, 1.0, 1.0));
  }
  return project2(x_col, xc, T, B, N);
}

std::vector<vec2> project_dihedral(const shell &M, const std::vector<vec3> &x,
                                   const index_t &corner, const vec3 &dir,
                                   const real &s) {
  index_t c0 = corner;
  index_t c1 = M.other(c0);
  index_t c0p = M.prev(c0);
  index_t c1p = M.prev(c1);

  real cot0 = cotan(M, c0, x);
  real cot1 = cotan(M, c1, x);

  vec3 x0 = x[M.vert(c0)];
  vec3 x0p = x[M.vert(c0p)];
  vec3 x1 = x[M.vert(c1)];
  vec3 x1p = x[M.vert(c1p)];
  vec3 p0 = va::mix(s, x0, x1);
  vec3 xc = 0.25 * (x0 + x1 + x0p + x1p);
  vec3 T = (x1 - x0).normalized();
  vec3 Nf0 = (x0 - x0p).cross(T).normalized();
  vec3 Nf1 = (x1p - x1).cross(T).normalized();
  vec3 Dn = dir.normalized().cross(T).normalized();
  vec3 N = Nf1;
  real cs = Nf0.dot(Nf1);
  real sn = (Nf0.cross(Nf1)).dot(T);
  real theta = atan2(sn, cs);

  if ((x0p - p0).dot(dir) > (x1p - p0).dot(dir)) {
    N = Nf0;
    theta *= -1.0;
    vec3 fc = 1.0 / 3.0 * (x0 + x1 + x0p);
    // gg::geometry_logger::line(fc, fc + 0.01 * N, vec4(1.0, 0.0, 0.0, 1.0));
  } else {
    vec3 fc = 1.0 / 3.0 * (x0 + x1 + x1p);
    // gg::geometry_logger::line(fc, fc + 0.01 * N, vec4(1.0, 0.0, 0.0, 1.0));
  }

  mat3 R = Eigen::AngleAxis<real>(theta, T).toRotationMatrix();
  vec3 Rdir = R * dir;
  vec3 p1 = p0 + Rdir;
  // p1 = R2 * (p1 - xc) + xc;

  // gg::geometry_logger::line(p0, p0 + 0.02 * Rdir, vec4(1.0, 1.0, 0.0, 1.0));
  vec3 B = (T.cross(N)).normalized();

  return project2({x0p, x0, x1p, x1, p0, p1}, xc, T, B, N);
}

using EdgeFcn =
    std::function<bool(const shell &, const std::vector<vec3> &,
                       const index_t &, const real &, const real &, vec3 &)>;

bool walk(const shell &M, const std::vector<vec3> &x, const index_t &start,
          const vec3 &dir0, const real &s0, const int max_iter = 1000,
          const real &tol = 1e-2, EdgeFcn fcn = nullptr) {

  auto log_seg = [&](index_t c00, index_t c01, index_t c10, index_t c11,
                     real s0, real s1, real n_o, vec4 col) {
    index_t v00 = M.vert(c00);
    index_t v01 = M.vert(c01);
    index_t v10 = M.vert(c10);
    index_t v11 = M.vert(c11);

    vec3 x0 = va::mix(s0, x[v00], x[v01]);
    vec3 x1 = va::mix(s1, x[v10], x[v11]);
    vec3 N0 = edge_normal(M, c00, x);
    vec3 N1 = edge_normal(M, c10, x);

    gg::geometry_logger::line(x0 + n_o * N0, x1 + n_o * N1, col);
  };

  auto get_one_ring = [&](index_t c) {
    std::vector<index_t> range;
    M.const_for_each_vertex(M.vert(c), [&range](int ci, const shell &M) {
      range.push_back(M.next(ci));
    });
    return range;
  };

  index_t corner = start;
  vec3 dir = dir0;
  real s = s0;
  real accumulated_length = 0.0;
  int it = 0;

  while (corner > -1 && it < max_iter) {

    index_t c0 = corner;
    index_t c1 = M.other(c0);
    index_t c0p = M.prev(c0);
    index_t c1p = M.prev(c1);

    vec3 x0 = x[M.vert(c0)];
    vec3 x1 = x[M.vert(c1)];
    vec3 xc;

    std::vector<index_t> C;
    std::vector<vec2> X2;

    if (fcn && !fcn(M, x, corner, s, accumulated_length, dir)) {
      break;
    }

    if (s < tol) {
      C = get_one_ring(c0);
      X2 = project_set(M, x, C, dir, c0);

      s = 0.0;
    } else if (s > 1.0 - tol) {
      C = get_one_ring(c1);
      X2 = project_set(M, x, C, dir, c1);
      s = 1.0;
    } else {
      C = {c0p, M.next(c1), c1p, M.next(c0)};
      X2 = project_dihedral(M, x, c0, dir, s);
    }

    int N = C.size();
    vec2 p0 = X2[N + 0];
    vec2 p1 = X2[N + 1];

    // gg::geometry_logger::line(unproj(p0, xc, B, T), unproj(p1, xc, B, T),
    //                           vec4(0.0, 1.0, 1.0, 1.0));
    corner = -1;

    for (int i = 0; i < N; i++) {
      int i0 = i;
      int i1 = (i + 1) % N;
      real si, yc;
      vec2 xi0 = X2[i0];
      vec2 xi1 = X2[i1];

      index_t ci0 = C[i0];
      index_t ci1 = C[i1];

      bool itx = LineLineIntersect(p0, p1, xi0, xi1, si, true);
      if (itx) {

        vec3 d1 = va::mix(si, x[M.vert(ci0)], x[M.vert(ci1)]);
        vec3 d0 = va::mix(s, x0, x1);
        accumulated_length += (d1 - d0).norm();
        dir = (d1 - d0).normalized();

        if (si < tol || si > 1.0 - tol) {
          si = std::round(si);
        }

        log_seg(c0, c1, ci0, ci1, s, si, 0.005, vec4(1.0, 0.0, 1.0, 1.0));
        corner = ci0;
        s = si;
      }
    }
    it++;
  }
  return it < max_iter;
}

real dist_from_line(const shell &M, const std::vector<vec3> &x,
                    const index_t &c0, const vec3 &pt) {
  index_t c1 = M.other(c0);
  vec3 x0 = x[M.vert(c0)];
  vec3 x1 = x[M.vert(c1)];
  return va::distance_from_line(x0, x1, pt);
}

bool contained_in_edge(const shell &M, const std::vector<vec3> &x,
                       const index_t &c0, const vec3 &pt) {

  index_t c1 = M.other(c0);
  vec3 x0 = x[M.vert(c0)];
  vec3 x1 = x[M.vert(c1)];
  return va::contained_in_line(x0, x1, pt);
}

real get_s(const shell &M, const std::vector<vec3> &x, const index_t &c0,
           const vec3 &pt) {
  index_t c1 = M.other(c0);
  vec3 x0 = x[M.vert(c0)];
  vec3 x1 = x[M.vert(c1)];
  real d0 = (pt - x0).norm();
  real d1 = (x1 - x0).norm();
  return d0 / d1;
}

index_t find(const shell &M, const std::vector<vec3> &x, const index_t &c0,
             const vec3 &pt,
             std::function<bool(const shell &M, const std::vector<vec3> &x,
                                const index_t &ci, const vec3 &pt, index_t &cn,
                                real &lmin)>
                 func) {

  std::vector<index_t> stack;
  std::set<index_t> v_visited;
  std::set<index_t> c_visited;

  stack.push_back(M.vert(c0));
  stack.push_back(M.vert(M.other(c0)));

  auto cmp = [&](const index_t &a, const index_t &b) {
    real da = (pt - x[a]).norm();
    real db = (pt - x[b]).norm();
    return da > db;
  };

  std::make_heap(stack.begin(), stack.end(), cmp);

  bool searching = true;
  index_t cn = -1;
  real lmin = 99999;
  index_t it = 0;
  std::cout << " " << M.vert(c0) << " " << M.vert(M.other(c0)) << std::endl;
  while (stack.size() > 0 && searching && it < 200) {

    std::pop_heap(stack.begin(), stack.end(), cmp);
    index_t v = stack.back();
    stack.pop_back();
    std::cout << " -" << v << " " << (pt - x[v]).norm() << std::endl;
    M.const_for_each_vertex(v, [&](int ci, const shell &M) {
      index_t vn = M.vert(M.next(ci));

      if (!has(v_visited, vn) && lmin > 1e-10) {
        v_visited.insert(vn);
        stack.push_back(vn);
        std::push_heap(stack.begin(), stack.end(), cmp);
      }
      if (!has(c_visited, ci / 2)) {
        it++;
        c_visited.insert(ci / 2);
        searching = func(M, x, ci, pt, cn, lmin);
      }
    });
  }

  return cn;
}

index_t find_edge(const shell &M, const std::vector<vec3> &x, const index_t &c0,
                  const vec3 &pt) {

  index_t cmin =
      find(M, x, c0, pt,
           [](const shell &M, const std::vector<vec3> &x, const index_t &ci,
              const vec3 &pt, index_t &cn, real &lmin) {
             real l = dist_from_line(M, x, ci, pt);
             bool in_edge = contained_in_edge(M, x, ci, pt);

             if (in_edge && l < lmin) {
               cn = ci;
               lmin = l;
               if (lmin < 1e12)
                 return false;
             }
             return true;
           });

  return cmin;
}

index_t find_point(const shell &M, const std::vector<vec3> &x,
                   const index_t &c0, const vec3 &pt) {

  index_t cmin =
      find(M, x, c0, pt,
           [](const shell &M, const std::vector<vec3> &x, const index_t &ci,
              const vec3 &pt, index_t &cn, real &lmin) {
             real l = (pt - x[M.vert(ci)]).norm();

             if (l < lmin) {

               cn = ci;
               lmin = l;
               if (lmin < 1e12)
                 return false;
             }
             return true;
           });

  return cmin;
}

std::vector<vec3> get_x(shell &M, const std::vector<vec3> &x,
                        std::vector<index_t> &E_in, std::vector<real> &S_in) {

  std::vector<vec3> x_out(E_in.size());
  for (int i = 0; i < E_in.size(); i++) {
    x_out[i] = edge_vert(M, E_in[i], S_in[i], x);
  }
  return x_out;
}

void log_adj_edges(shell &M, const std::vector<vec3> &x, index_t v) {
  M.const_for_each_vertex(v, [&](index_t c0, const shell &M) {
    index_t v0 = M.vert(c0);
    index_t v1 = M.vert(M.other(c0));

    vec3 x0 = x[v0];
    vec3 x1 = x[v1];
    vec3 N0 = vert_normal(M, v0, x);
    vec3 N1 = vert_normal(M, v1, x);
    real n_o = 0.0025;
    vec4 col = vec4(0.75, 0.25, 0.5, 1.0);
    gg::geometry_logger::line(x0 + n_o * N0, x1 + n_o * N1, col);
  });
}

bool edge_between(shell &M, index_t v0, index_t v1) {
  bool has_edge = false;
  M.const_for_each_vertex(v0, [&](index_t c0, const shell &M) {
    index_t vn = M.vert(M.next(c0));
    has_edge = has_edge || (vn == v1);
  });
  return has_edge;
};

void walk_new_seg(shell &M, const std::vector<vec3> &x,            //
                  index_t c0, index_t c1, std::vector<index_t> &E, //
                  std::vector<real> &S, const real &tol = 1e-2) {
  vec3 gdir = (x[M.vert(c1)] - x[M.vert(c0)]).normalized();
  asawa::shell::walk(
      M, x, c0, gdir, 0.0, 20, tol,
      [&](const asawa::shell::shell &M, const std::vector<vec3> &x,
          const index_t &ci, const real &s, const real &accumulated_length,
          vec3 &dir) {
        dir = gdir;
        real si = s;

        if (M.vert(ci) == M.vert(c0) || M.vert(M.other(ci)) == M.vert(c0))
          return true;
        if (M.vert(ci) == M.vert(c1) || M.vert(M.other(ci)) == M.vert(c1))
          return false;

        S.push_back(s);
        E.push_back(ci);

        return true;
      });
}

inline void log_edge(shell &M, const std::vector<vec3> &x, index_t c0, real n_o,
                     vec4 col) {
  index_t c1 = M.other(c0);

  index_t v0 = M.vert(c0);
  index_t v1 = M.vert(c1);

  vec3 x0 = x[v0];
  vec3 x1 = x[v1];
  vec3 N0 = edge_normal(M, c0, x);
  gg::geometry_logger::line(x0 + n_o * N0, x1 + n_o * N0, col);
};

inline void log_seg_v(shell &M, const std::vector<vec3> &x, //
                      index_t v0, index_t v1, real n_o, vec4 col) {
  vec3 x0 = x[v0];
  vec3 x1 = x[v1];
  vec3 N0 = vert_normal(M, v0, x);
  vec3 N1 = vert_normal(M, v1, x);
  gg::geometry_logger::line(x0 + n_o * N0, x1 + n_o * N1, col);
};

inline void log_seg_e(shell &M, const std::vector<vec3> &x, //
                      const index_t &c0, const index_t &c1, //
                      const real &s0, const real &s1,       //
                      real n_o, vec4 col) {
  vec3 x0 = edge_vert(M, c0, s0, x);
  vec3 x1 = edge_vert(M, c1, s1, x);
  vec3 N0 = edge_normal(M, c0, x);
  vec3 N1 = edge_normal(M, c1, x);

  gg::geometry_logger::line(x0 + n_o * N0, x1 + n_o * N1, col);
};

std::vector<index_t> stitch_walk(shell &M, const std::vector<vec3> &x,
                                 const std::vector<index_t> &walk_0) {
  std::vector<index_t> walk;
  std::vector<bool> has_walk(M.vert_count(), false);

  for (int i = 0; i < walk_0.size() - 1; i++) {
    index_t v0 = M.vert(walk_0[i + 0]);
    has_walk[v0] = true;
  }

  for (int i = 0; i < walk_0.size() - 1; i++) {
    index_t v0 = M.vert(walk_0[i + 0]);
    index_t v1 = M.vert(walk_0[i + 1]);
    log_seg_v(M, x, v0, v1, 0.0045, vec4(1.0, 1.0, 0.0, 1.0));

    asawa::index_t cm = -1;

    M.const_for_each_vertex(v0, [&](index_t c0, const shell &M) {
      index_t vn = M.vert(M.next(c0));
      if (vn == v1) {
        cm = c0;
      }
    });

    vec3 dv0 = (x[v1] - x[v0]).normalized();
    std::cout << v0 << " " << v1 << std::endl;
    int k = 0;

    if (cm == -1) {
      std::cout << " gah! " << std::endl;
      log_seg_v(M, x, v0, v1, 0.0025, vec4(1.0, 0.0, 0.0, 1.0));
    }

    while (cm == -1 && k < 10) {
      index_t co = -1;
      index_t vo = -1;
      real max = -1e12;
      M.const_for_each_vertex(v0, [&](index_t c0, const shell &M) {
        index_t vn = M.vert(M.next(c0));
        vec3 dvn = (x[vn] - x[v0]).normalized();
        real d = dvn.dot(dv0);
        if (d > max) {
          max = d;
          co = c0;
          vo = vn;
          std::cout << " d: " << d << " " << co << " " << vn << std::endl;
        }
        if (vn == v1) {
          cm = c0;
          max = 2.0;
        }
      });
      v0 = vo;
      walk.push_back(co);
      std::cout << co << " " << cm << " " << v0 << std::endl;
      k++;
    }
    std::cout << " cm: " << cm << std::endl;
    walk.push_back(cm);
  }

  return  walk;
}

inline bool stol(real s, real tol) { return s < tol || s > 1.0 - tol; };

index_t get_seg_range(index_t t0, shell &M, const std::vector<vec3> &x,
                      const std::vector<vec3> &X, std::vector<index_t> &E,
                      std::vector<real> &S, std::vector<index_t> &sub_verts,
                      real tol) {
  std::set<index_t> visited;
  bool testing = true;
  index_t i = t0;
  index_t t1 = t0;

  while (testing && i < E.size()) {
    real s = S[i];
    index_t c0 = E[i];
    index_t c1 = M.other(c0);

    vec3 xi = X[i];
    vec3 x0 = x[M.vert(c0)];
    vec3 x1 = x[M.vert(c1)];
    bool in_line = va::contained_in_line(x0, x1, xi);

    if (!in_line) {
      c0 = find_edge(M, x, c0, xi);
      c1 = M.other(c0);
      E[i] = c0;
    }

    s = get_s(M, x, c0, xi);
    S[i] = s;

    if (has(visited, c0 / 2) || stol(s, tol)) {
      std::cout << " set_verts: " << i << ": " << sub_verts[i] << " " << c0
                << std::endl;

      if (s < 0.5)
        sub_verts[i] = c0;
      else
        sub_verts[i] = c1;
      if (has(visited, c0 / 2))
        log_edge(M, x, c0, 0.0025, vec4(1.0, 0.5, 0.0, 1.0));
      else
        log_edge(M, x, c0, 0.0025, vec4(0.5, 0.0, 1.0, 1.0));
      testing = false;
      t1 = i;
    }
    if (i == E.size() - 1) {
      t1 = E.size();
      testing = false;
    }

    visited.insert(c0 / 2);
    i++;
  }

  return t1;
}

std::vector<index_t> subdivide_seg(shell &M, const std::vector<vec3> &x, //
                                   std::vector<index_t> &Ei,
                                   std::vector<real> &Si) {

  std::vector<std::array<index_t, 1>> l_verts_a =
      subdivide_op(M, Ei, Si, x,
                   [](index_t i,  //
                      index_t cs, //
                      index_t vs, //
                      index_t fs, //
                      const std::vector<index_t> &edges, shell &m) -> corner1 {
                     return {subdivide_edge(m, edges[i],    //
                                            vs + i,         //
                                            cs + 6 * i + 0, //
                                            cs + 6 * i + 2, //
                                            cs + 6 * i + 4, //
                                            fs + 2 * i + 0, //
                                            fs + 2 * i + 1  //
                                            )};
                   });

  // std::cout << "l_verts: " << l_verts.size() << std::endl;

  for (int i = 0; i < l_verts_a.size() - 1; i++) {
    // flip segments if there isn't an edge between them
    // log_seg_v(i, 0.01, vec4(0.0, 1.0, 1.0, 1.0), l_verts);
  }
  std::vector<index_t> l_verts;
  std::transform(l_verts_a.begin(), l_verts_a.end(),
                 std::back_inserter(l_verts),
                 [&](const std::array<index_t, 1> &x) { return x[0]; });

  return l_verts;
}

std::vector<index_t> crack_edges(shell &M, std::vector<index_t> &E,
                                 std::vector<real> &S, const real &tol = 1e-2) {

  using comp_great = shell_data_comp<vec3, std::greater<real>>;
  const std::vector<vec3> &x = get_vec_data(M, 0);

  auto lam_log_seg_vv = [&](index_t i, real n_o, vec4 col,
                            const std::vector<index_t> &l_verts) {
    if (l_verts.size() < 2)
      return;
    log_seg_v(M, x, M.vert(l_verts[i + 0]), M.vert(l_verts[i + 1]), n_o, col);
  };

  std::vector<vec3> X = get_x(M, x, E, S);
  std::vector<index_t> sub_verts(E.size(), -1);

  index_t t0 = 0;
  index_t t1 = t0;
  int k = 0;

  // log_edge(M, x, E[t0], 0.0025, vec4(1.0, 0.0, 0.0, 1.0));

  while (t0 < E.size()) {
    std::cout << "k: " << k << std::endl;
    k++;

    while (stol(S[t0], tol) && t0 < E.size()) {
      real s = S[t0];
      index_t c0 = E[t0];
      index_t c1 = M.other(c0);
      s = round(s);
      S[t0] = s;
      if (s < 0.5)
        sub_verts[t0] = c0;
      else
        sub_verts[t0] = c1;
      t0++;
    }

    index_t i = t0;
    t1 = get_seg_range(t0, M, x, X, E, S, sub_verts, tol);
    // std::cout << "t0/t1: " << t0 << " " << t1 << std::endl;
    if (t0 == t1) {
      t0 = t1;
      continue;
    }

    std::vector<real> Si(S.begin() + t0, S.begin() + t1);
    std::vector<index_t> Ei(E.begin() + t0, E.begin() + t1);

    std::vector<index_t> l_verts = subdivide_seg(M, x, Ei, Si);
    std::cout << "t: " << t0 << " " << l_verts.size() + t0 << " "
              << sub_verts.size() << " " << l_verts.size() << std::endl;
    std::copy(l_verts.begin(), l_verts.end(), sub_verts.begin() + t0);

    // if (t1 < E.size() - 1)
    //   log_edge(M, x, E[t1], 0.0025, vec4(0.0, 0.0, 1.0, 1.0));

    t0 = t1;
  }

  t0 = 0;
  t1 = t0;
  for (int i = 0; i < sub_verts.size(); i++) {
    std::cout << sub_verts[i] << " ";
  }
  std::cout << std::endl;

  for (int i = 0; i < sub_verts.size() - 1; i++) {
    index_t c0 = sub_verts[i];
    std::cout << i << " " << c0 << std::endl;
    if ((X[i] - x[M.vert(c0)]).norm() > 1e-8)
      sub_verts[i] = find_point(M, x, c0, X[i]);
  }

  std::vector<index_t> sub_verts_2;
  std::vector<vec3> X2;
  while (t0 < sub_verts.size()) {
    std::cout << t1 << " " << sub_verts.size() << std::endl;

    if (t1 > sub_verts.size() - 2)
      break;

    bool eb = edge_between(M, M.vert(sub_verts[t1]), M.vert(sub_verts[t1 + 1]));
    if (!eb) {
      std::copy(sub_verts.begin() + t0, sub_verts.begin() + t1 + 1,
                std::back_inserter(sub_verts_2));
      std::copy(X.begin() + t0, X.begin() + t1 + 1, std::back_inserter(X2));
      std::vector<index_t> Ei;
      std::vector<real> Si;

      walk_new_seg(M, x, sub_verts[t1], sub_verts[t1 + 1], Ei, Si);
      std::vector<vec3> Xi = get_x(M, x, Ei, Si);
      std::vector<index_t> l_verts = subdivide_seg(M, x, Ei, Si);

      std::copy(l_verts.begin(), l_verts.end(),
                std::back_inserter(sub_verts_2));

      std::copy(Xi.begin(), Xi.end(), std::back_inserter(X2));

      t0 = t1 + 1;
      t1 = t0;
    } else
      t1++;
  }

  // fix them again...
  for (int i = 0; i < sub_verts_2.size() - 1; i++) {
    index_t c0 = sub_verts_2[i];
    if ((X2[i] - x[M.vert(c0)]).norm() > 1e-8)
      sub_verts_2[i] = find_point(M, x, c0, X2[i]);
  }

#if 1
  for (int i = 0; i < sub_verts_2.size() - 1; i++) {
    lam_log_seg_vv(i, 0.03, vec4(0.0, 1.0, 0.0, 1.0), sub_verts_2);
    // log_adj_edges(M, x, M.vert(sub_verts_2[i]));
  }
#endif
  return sub_verts_2; // stitch_walk(M, x, sub_verts_2);
}

} // namespace shell

} // namespace asawa
} // namespace gaudi
#endif