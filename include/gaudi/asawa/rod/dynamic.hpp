
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
real point_line_min(const index_t &idT, //
                    const std::vector<index_t> &t_inds,
                    const vector<vec3> &t_x, //
                    const index_t &idS,      //
                    const std::vector<index_t> &s_inds,
                    const vector<vec3> &s_x) {

  if (t_inds[idT] == s_inds[2 * idS + 0])
    return std::numeric_limits<real>::infinity();
  if (t_inds[idT] == s_inds[2 * idS + 1])
    return std::numeric_limits<real>::infinity();

  const vec3 &xA = t_x[t_inds[idT]];
  const vec3 &xB0 = s_x[s_inds[2 * idS + 0]];
  const vec3 &xB1 = s_x[s_inds[2 * idS + 1]];
  vec3 xB = 0.5 * (xB0 + xB1);
  vec3 dAB = (xA - xB).normalized();
  vec3 dBB = (xB1 - xB0).normalized();
  real cosAB = dAB.dot(dBB);

  if (abs(cosAB) > 0.1)
    return std::numeric_limits<real>::infinity();

  return va::distance_from_line(xB0, xB1, xA);
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

  void interp_l(index_t c0, index_t c1, index_t cnew, std::vector<real> &x) {
    real l = x[c0];
    x[c0] = 0.5 * l;
    set<real>(cnew, 0.5 * l, x);
  }

  void split_edge(index_t c0) {
    index_t c1 = __R->next(c0);
    index_t cnew = __R->insert_edge();
    __R->link(c0, cnew);
    __R->link(cnew, c1);
    interp(c0, c1, cnew, __R->__M);
    interp(c0, c1, cnew, __R->__J);
    interp(c0, c1, cnew, __R->__x);
    interp(c0, c1, cnew, __R->__v);
    interp(c0, c1, cnew, __R->__u);
    interp(c0, c1, cnew, __R->__o);
    interp_l(c0, c1, cnew, __R->__l0);
  }

#if 1
  vector<std::array<index_t, 2>> get_collisions() {
    rod &R = *__R;
    std::vector<vec3> &x = R.__x;

    std::vector<index_t> verts = __R->get_vert_range();
    std::vector<index_t> edge_verts = __R->get_edge_vert_ids();
    std::vector<index_t> edge_map = __R->get_edge_map();

    edge_tree = arp::aabb_tree<2>::create(edge_verts, x, 16);
    // calder::test_extents(*edge_tree, edge_verts, x);
    // edge_tree->debug();
    real tol = R._r;

    std::vector<std::array<index_t, 2>> collected(verts.size());
    //#pragma omp parallel for
    for (int i = 0; i < verts.size(); i++) {
      index_t j0 = arp::getNearest<1, 2>(i, verts, x, //
                                         *edge_tree,  //
                                         tol, &point_line_min);
      if (j0 > 0) {
        vec3 xA = x[i];
        vec3 xB0 = x[edge_verts[2.0 * j0 + 0]];
        vec3 xB1 = x[edge_verts[2.0 * j0 + 1]];
        vec3 dAB = xB1 - xB0;

        real s = (xA - xB0).dot(dAB) / dAB.dot(dAB);
        vec3 xAB = xB0 + s * dAB;
        gg::geometry_logger::line(xA, xAB, vec4(1.0, 0.0, 1.0, 1.0));
      }
      collected[i] = {i, j0};
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