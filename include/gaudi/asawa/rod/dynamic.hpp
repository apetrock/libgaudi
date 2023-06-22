
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "gaudi/arp/arp.h"

#include "gaudi/common.h"
#include "gaudi/vec_addendum.h"

#include "GaudiGraphics/geometry_logger.h"

#include "rod.hpp"

// #include "subdivide.hpp"

#include <array>
#include <set>

#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>
#include <zlib.h>

#ifndef __ASAWA_DYNAMIC_ROD__
#define __ASAWA_DYNAMIC_ROD__

namespace gaudi {
namespace asawa {
namespace rod {
// this is ugly, but we have to do it this way, with two lists because
// the callback on the data
real line_line(const index_t &idT, //
               const std::vector<index_t> &t_inds,
               const vector<vec3> &t_x, //
               const index_t &idS,      //
               const std::vector<index_t> &s_inds, const vector<vec3> &s_x) {

  if (t_inds[2 * idT + 0] == s_inds[2 * idS + 0])
    return std::numeric_limits<real>::infinity();
  if (t_inds[2 * idT + 1] == s_inds[2 * idS + 1])
    return std::numeric_limits<real>::infinity();

  if (t_inds[2 * idT + 1] == s_inds[2 * idS + 0])
    return std::numeric_limits<real>::infinity();
  if (t_inds[2 * idT + 0] == s_inds[2 * idS + 1])
    return std::numeric_limits<real>::infinity();
  //  std::cout << t_inds[2 * idT + 0] << " " << t_inds[2 * idT + 1] << "-"
  //            << t_inds[2 * idS + 0] << " " << t_inds[2 * idS + 1] <<
  //            std::endl;
  const vec3 &xA0 = t_x[t_inds[2 * idT + 0]];
  const vec3 &xA1 = t_x[t_inds[2 * idT + 1]];

  const vec3 &xB0 = s_x[s_inds[2 * idS + 0]];
  const vec3 &xB1 = s_x[s_inds[2 * idS + 1]];

  std::array<real, 3> d = va::distance_Segment_Segment(xA0, xA1, xB0, xB1);
  real s = d[1];
  real t = d[2];

  vec3 xA = va::mix(s, xA0, xA1);
  vec3 xB = va::mix(t, xB0, xB1);
  vec3 dA = (xA1 - xA0).normalized();
  vec3 dB = (xB1 - xB0).normalized();

  vec3 xAB = (xB - xA).normalized();
  // gg::geometry_logger::line(xA, xB, vec4(1.0, 1.0, 1.0, 1.0));
  /*
    real s = d[1];
    real t = d[2];
    vec3 xA = (1.0 - s) * xA0 + s * xA1;
    vec3 xB = (1.0 - t) * xB0 + t * xB1;

    vec3 dAB = (xA - xB).normalized();
  */

  if (abs(dA.dot(xAB)) > 0.95)
    return std::numeric_limits<real>::infinity();
  if (abs(dB.dot(xAB)) > 0.95)
    return std::numeric_limits<real>::infinity();

  return d[0];
};

class dynamic {
public:
  typedef std::shared_ptr<dynamic> ptr;

  static ptr create(rod::ptr R, real Cc, real Cs, real Cm) {
    ptr Rd = std::make_shared<dynamic>(R, Cc, Cs, Cm);
    return Rd;
  }

  dynamic(rod::ptr R, real Cc, real Cs, real Cm)
      : _Cc(Cc), _Cs(Cs), _Cm(Cm), __R(R) {}

  template <typename T>
  void set(const index_t &cnew, const T &xnew, std::vector<T> &x) {
    if (cnew >= x.size()) {
      x.push_back(xnew);
    } else {
      x[cnew] = xnew;
    }
  }

  template <typename T> T avg(index_t c0, index_t c1, const std::vector<T> &x) {
    T x0 = x[c0];
    T x1 = x[c1];
    T xnew = 0.5 * (x0 + x1);
    return xnew;
  }

  void interp(index_t c0, index_t c1, index_t cnew, std::vector<real> &x) {
    real xnew = avg<real>(c0, c1, x);
    set<real>(cnew, xnew, x);
  }

  void interp(index_t c0, index_t c1, index_t cnew, std::vector<vec3> &x) {
    vec3 xnew = avg<vec3>(c0, c1, x);
    set<vec3>(cnew, xnew, x);
  }

  void interp(index_t cp, index_t c0, index_t c1, index_t c2, index_t cnew,
              std::vector<vec3> &x) {

    if (cp < 0 || c2 < 0) {
      vec3 xnew = avg<vec3>(c0, c1, x);
      set<vec3>(cnew, xnew, x);
    } else {
      vec3 xp = x[cp];
      vec3 x0 = x[c0];
      vec3 x1 = x[c1];
      vec3 x2 = x[c2];
      vec3 xnew = va::catmull_rom(xp, x0, x1, x2, 0.5);
      set<vec3>(cnew, xnew, x);
    }
  }

  void interp(index_t c0, index_t c1, index_t cnew, std::vector<vec4> &x) {
    vec4 xnew = avg<vec4>(c0, c1, x);
    set<vec4>(cnew, xnew, x);
  }

  void interp(index_t c0, index_t c1, index_t cnew, std::vector<quat> &x) {
    quat x0 = x[c0];
    quat x1 = x[c1];
    quat xnew = va::slerp(x0, x1, 0.5);
    set<quat>(cnew, xnew, x);
  }

  void interp_l(index_t c0, index_t c1, index_t cnew, std::vector<real> &z) {

    std::cout << c0 << " " << c1 << " " << cnew << std::endl;
    std::cout << __R->__x.size() << " " << __R->__corners_next.size()
              << std::endl;
    vec3 x0 = __R->__x[c0];
    vec3 x1 = __R->__x[c1];
    vec3 xn = __R->__x[cnew];
    real l10 = (x1 - x0).norm();
    real ln0 = (xn - x0).norm();
    real l1n = (x1 - xn).norm();
    real lnew = ln0 + l1n;

    real l0 = z[c0];
    // real r = 0.0;
    // real l = va::mix(r, l0, lnew);
    //  std::cout << " interp: " << (ln0 + l1n) / l10 << " " << ln0 / l10 << " "
    //            << l1n / l10 << std::endl;

    // set<real>(c0, ln0 / l10 * l, z);
    // set<real>(cnew, l1n / l10 * l, z);
    set<real>(c0, ln0 / l10 * l0, z);
    set<real>(cnew, l1n / l10 * l0, z);
  }

  void split_edge(index_t c00) {
    index_t c10 = __R->next(c00);
    index_t c11 = __R->next(c10);
    index_t c01 = __R->prev(c00);
    index_t cnew = __R->insert_edge();
    std::cout << c01 << " " << c00 << " " << c10 << " " << c11 << " " << cnew
              << " " << std::endl;
    __R->link(c00, cnew);
    __R->link(cnew, c10);
    std::cout << "interp0" << std::endl;
    std::cout << __R->__x.size() << std::endl;
    std::cout << __R->__v.size() << std::endl;
    std::cout << __R->__u.size() << std::endl;
    std::cout << __R->__o.size() << std::endl;
    std::cout << __R->__M.size() << std::endl;
    std::cout << __R->__J.size() << std::endl;
    std::cout << __R->__l0.size() << std::endl;

    interp(c01, c00, c10, c11, cnew, __R->__x);
    // interp(c00, c10, cnew, __R->__x);
    std::cout << "interp1" << std::endl;
    interp(c00, c10, cnew, __R->__v);
    interp(c00, c10, cnew, __R->__u);
    interp(c00, c10, cnew, __R->__o);
    interp(c00, c10, cnew, __R->__M);
    interp(c00, c10, cnew, __R->__J);
    std::cout << "interp2" << std::endl;
    interp_l(c00, c10, cnew, __R->__l0);
  }

  void collapse_edge(index_t c00) {
    index_t c10 = __R->next(c00);
    index_t c11 = __R->next(c10);
    index_t c01 = __R->prev(c00);
    if (c10 < -1)
      return;

    __R->link(c01, c10);
    __R->set_next(c00, -1);
    __R->set_prev(c00, -1);

    interp(c01, c10, c10, __R->__x);
    // interp(c00, c10, cnew, __R->__x);

    interp(c01, c10, c10, __R->__v);
    interp(c01, c10, c10, __R->__u);
    interp(c01, c10, c10, __R->__o);
    interp(c01, c10, c10, __R->__M);
    interp(c01, c10, c10, __R->__J);
    interp_l(c01, c10, c10, __R->__l0);
  }

#if 1
  vector<std::array<index_t, 4>>
  get_collisions(const std::vector<index_t> &edge_verts_B, //
                 const std::vector<vec3> &x_B) {
    rod &R = *__R;

    std::vector<vec3> &x_A = R.__x;
    std::vector<index_t> verts_A = R.get_vert_range();
    std::vector<index_t> edge_verts_A = R.get_edge_vert_ids();
    std::vector<index_t> edge_map_A = R.get_edge_map();

    edge_tree = arp::aabb_tree<2>::create(edge_verts_A, x_A, 16);
    // calder::test_extents(*edge_tree, edge_verts, x);
    // edge_tree->debug();
    real tol = 1.0 * R._r;

    std::cout << " 0 " << std::endl;

    std::vector<std::array<index_t, 4>> collected(edge_verts_B.size() / 2);
    for (int k = 0; k < edge_verts_B.size(); k += 2) {
      index_t i = k / 2;

      std::vector<index_t> collisions =
          arp::getNearest<2, 2>(i, edge_verts_B, x_B, //
                                *edge_tree,           //
                                tol, &line_line);

      for (index_t j : collisions) {
        if (j > -1) {
          collected[i] = {edge_verts_B[2 * i + 0], edge_verts_B[2 * i + 1],
                          edge_verts_A[2 * j + 0], edge_verts_A[2 * j + 1]};
        } else {
          collected[i] = {-1, -1, -1, -1};
        }
      }
    }

    return collected;
  }
#endif

#if 1
  vector<std::array<index_t, 4>> get_internal_collisions() {
    rod &R = *__R;
    std::vector<vec3> &x = R.__x;

    std::vector<index_t> verts = R.get_vert_range();
    std::vector<index_t> edge_verts = R.get_edge_vert_ids();

    return get_collisions(edge_verts, x);
  }

#endif
  void step() {

    for (int i = 0; i < __R->corner_count(); i++) {
      index_t j = __R->next(i);
      if (j < 0)
        continue;

      std::cout << i << " " << j << std::endl;
      vec3 q0 = __R->__x[i];
      vec3 q1 = __R->__x[j];
      real l = (q1 - q0).norm();
      if (l > _Cs) {
        std::cout << " split " << std::endl;
        split_edge(i);
      }

      if (l < _Cc) {
        std::cout << " collapse " << std::endl;
        collapse_edge(i);
      }
    }
    __R->update_mass();
  }
  rod::ptr __R;
  arp::aabb_tree<2>::ptr edge_tree;
  real _Cc, _Cs, _Cm; // collapse, stretch, bridge
  // arp::aabb_tree<2>::ptr edge_tree;
};
} // namespace rod
} // namespace asawa
} // namespace gaudi

#endif