
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "gaudi/arp/arp.h"

#include "gaudi/common.h"
#include "gaudi/vec_addendum.h"

#include "GaudiGraphics/geometry_logger.h"

#include "rod.hpp"

//#include "subdivide.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <type_traits>
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
  /*
    real s = d[1];
    real t = d[2];
    vec3 xA = (1.0 - s) * xA0 + s * xA1;
    vec3 xB = (1.0 - t) * xB0 + t * xB1;
    vec3 dA = (xA1 - xA0).normalized();
    vec3 dB = (xB1 - xB0).normalized();
    vec3 dAB = (xA - xB).normalized();
    if (abs(xA.dot(dAB)) < 0.1)
      return std::numeric_limits<real>::infinity();
    if (abs(xB.dot(dAB)) < 0.1)
      return std::numeric_limits<real>::infinity();
  */
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
    quat xnew = x0.slerp(0.5, x1);
    set<quat>(cnew, xnew, x);
  }

  void interp_l(index_t c0, index_t c1, index_t cnew, std::vector<real> &z) {
    vec3 x0 = __R->__x[c0];
    vec3 x1 = __R->__x[c1];
    vec3 xn = __R->__x[cnew];
    real l10 = (x1 - x0).norm();
    real ln0 = (xn - x0).norm();
    real l1n = (x1 - xn).norm();

    real l = z[c0];
    // std::cout << " interp: " << (ln0 + l1n) / l10 << " " << ln0 / l10 << " "
    //           << l1n / l10 << std::endl;
    /// set<real>(c0, ln0 / l10, z);
    // set<real>(cnew, l1n / l10, z);
    set<real>(c0, ln0 / l10 * l, z);
    set<real>(cnew, l1n / l10 * l, z);
  }

  void split_edge(index_t c00) {
    index_t c10 = __R->next(c00);
    index_t c11 = __R->next(c10);
    index_t c01 = __R->prev(c00);
    index_t cnew = __R->insert_edge();
    __R->link(c00, cnew);
    __R->link(cnew, c10);

    interp(c01, c00, c10, c11, cnew, __R->__x);
    // interp(c00, c10, cnew, __R->__x);

    interp(c00, c10, cnew, __R->__v);
    interp(c00, c10, cnew, __R->__u);
    interp(c00, c10, cnew, __R->__o);
    interp(c00, c10, cnew, __R->__M);
    interp(c00, c10, cnew, __R->__J);
    interp_l(c00, c10, cnew, __R->__l0);
  }

#if 1
  vector<std::array<index_t, 4>> get_collisions() {
    rod &R = *__R;
    std::vector<vec3> &x = R.__x;

    std::vector<index_t> verts = __R->get_vert_range();
    std::vector<index_t> edge_verts = __R->get_edge_vert_ids();
    std::vector<index_t> edge_map = __R->get_edge_map();

    edge_tree = arp::aabb_tree<2>::create(edge_verts, x, 16);
    // calder::test_extents(*edge_tree, edge_verts, x);
    // edge_tree->debug();
    real tol = 0.95 * R._r;

    std::vector<std::array<index_t, 4>> collected(verts.size());
    //#pragma omp parallel for
    for (int k = 0; k < edge_verts.size(); k += 2) {
      index_t i = k / 2;
      std::vector<index_t> collisions =
          arp::getNearest<2, 2>(i, edge_verts, x, //
                                *edge_tree,       //
                                tol, &line_line);
      for (index_t j : collisions) {
        if (j > 0) {
          collected[i] = {edge_verts[2.0 * i + 0], edge_verts[2.0 * i + 1],
                          edge_verts[2.0 * j + 0], edge_verts[2.0 * j + 1]};
        } else
          collected[i] = {-1, -1, -1, -1};
      }
    }

    return collected;
  }
#endif
  void step() {

    for (int i = 0; i < __R->corner_count(); i++) {
      index_t j = __R->next(i);
      vec3 q0 = __R->__x[i];
      vec3 q1 = __R->__x[j];
      real l = (q1 - q0).norm();
      if (l > _Cs)
        split_edge(i);
    }
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