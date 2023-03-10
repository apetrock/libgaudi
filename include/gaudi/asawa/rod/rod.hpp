

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "gaudi/arp/arp.h"

#include "gaudi/common.h"
#include "gaudi/vec_addendum.h"

#include "GaudiGraphics/geometry_logger.h"

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

#ifndef __ASAWA_ROD__
#define __ASAWA_ROD__

namespace gaudi {
namespace asawa {
namespace rod {

class rod {
public:
  typedef std::shared_ptr<rod> ptr;

  static ptr create() {
    ptr R = std::make_shared<rod>();
    return R;
  }
  static ptr create(const std::vector<vec3> &verts, bool loop = true) {
    ptr R = std::make_shared<rod>();
    R->insert_strand(verts, loop);
    return R;
  }

  rod() {}

  void insert_strand(const std::vector<vec3> &verts, bool loop) {

    __x.insert(__x.end(), verts.begin(), verts.end());

    int vN = verts.size();
    int cbegin = __corners_prev.size();

    for (int i = 0; i < vN; i++) {
      if (i == 0) {
        if (loop) {
          __corners_prev.push_back(cbegin + vN - 1);
          __corners_next.push_back(cbegin + i + 1);
        } else {
          __corners_prev.push_back(-1);
          __corners_next.push_back(cbegin + 1);
        }
      } else if (i == vN - 1) {
        if (loop) {
          __corners_prev.push_back(cbegin + i - 1);
          __corners_next.push_back(cbegin + 0);
        } else {
          __corners_prev.push_back(cbegin + i - 1);
          __corners_next.push_back(-1);
        }
      } else {
        __corners_prev.push_back(cbegin + i - 1);
        __corners_next.push_back(cbegin + i + 1);
      }
    }
    _init_params();
  }

  void _init_params() {
    __l0.resize(__corners_next.size());
    __v.resize(__corners_next.size());
    __M.resize(__corners_next.size());
    __J.resize(__corners_next.size());
    __o.resize(__corners_next.size());
    _lmax = 0.0;
    for (int i = 0; i < __corners_next.size(); i++) {
      if (__corners_next[i] == -1)
        continue;
      index_t j = __corners_next[i];
      vec3 q0 = __x[i];
      vec3 q1 = __x[j];
      real l = (q1 - q0).norm();
      _lmax += l;
      __l0[i] = l;
      __v[i] = vec3::Zero();
      real rho = 1.0;

      real M = M_PI * _r * _r * l;
      real J = l * rho * M * _r * _r;

      __M[i][0] = M;
      __M[i][1] = M;
      __M[i][2] = M;

      __J[i][0] = 0.25 * J;
      __J[i][1] = 0.25 * J;
      __J[i][2] = 0.5 * J;
      __J[i][3] = 0.0;
      /*
      __J[i][0] = 0.0;
      __J[i][1] = 0.0;
      __J[i][2] = 0.0;
      __J[i][3] = 0.0;
*/
    }
    _lmax /= real(__corners_next.size());
    _update_frames();
  }

  quat _calc_frame(const index_t &ip, const index_t &i, const index_t &in) {
    vec3 cp = __x[ip];
    vec3 c0 = __x[i];
    vec3 cn = __x[in];

    vec3 dn0 = cn - c0;
    vec3 d0p = c0 - cp;

    vec3 d0 = dn0.cross(d0p).normalized();
    vec3 d1 = dn0.cross(d0).normalized();
    vec3 d2 = d0.cross(d1).normalized();

    mat3 F;
    F.col(0) = d0;
    F.col(1) = d1;
    F.col(2) = d2;

    quat q(F);
    q.normalize();

    return q;
  }

  void _update_frames() {
    __u.resize(__corners_next.size());
    for (int i = 0; i < __corners_next.size(); i++) {

      if (next(i) == -1)
        __u[i] = _calc_frame(prev(prev(i)), prev(i), i);
      else if (prev(i) == -1)
        __u[i] = _calc_frame(i, next(i), next(next(i)));
      else
        __u[i] = _calc_frame(prev(i), i, next(i));

      index_t j = __corners_next[i];
      vec3 q0 = __x[i];
      vec3 q1 = __x[j];
      vec3 dq = q1 - q0;
      dq.normalize();
      //__o[i] = quat(0.0, dq[0], dq[1], dq[2]);
      __o[i] = quat(0.0, 0.0, 0.0, 0.0);
    }
  }
  real lavg() {
    real lbar = 0.0;
    for (int i = 0; i < corner_count(); i++) {
      if (next(i) == -1)
        continue;
      index_t j = next(i);
      vec3 q0 = __x[i];
      vec3 q1 = __x[j];
      real l = (q1 - q0).norm();
      lbar += l;
    }
    lbar /= real(corner_count());
    return lbar;
  }

  index_t next(index_t i) const { return i < 0 ? -1 : __corners_next[i]; }
  index_t prev(index_t i) const { return i < 0 ? -1 : __corners_prev[i]; }

  void set_next(index_t id, index_t c) { __corners_next[id] = c; }
  void set_prev(index_t id, index_t c) { __corners_prev[id] = c; }
  void link(index_t c0, index_t c1) {
    set_next(c0, c1);
    set_prev(c1, c0);
  }

  size_t corner_count() const { return __corners_next.size(); }
  index_t insert_edge() {
    // returns id of new vertex
    __corners_next.push_back(-1);
    __corners_prev.push_back(-1);
    return __corners_next.size() - 1;
  }

  std::vector<vec3> &corner_verts() { return __x; }
  const std::vector<vec3> &corner_verts() const { return __x; }

  std::vector<index_t> get_edge_vert_ids() {
    std::vector<index_t> range;
    range.reserve(corner_count());
    // replace this with some c++isms
    for (int i = 0; i < corner_count(); i += 2) {
      if (__corners_next[i] < 0)
        continue;
      range.push_back(i);
      range.push_back(__corners_next[i]);
    }
    return range;
  }

  std::vector<index_t> get_range_map(const std::vector<index_t> indices,
                                     int stride) {
    std::vector<index_t> map;
    map.reserve(indices.size() / stride);
    // replace this with some c++isms
    for (int i = 0; i < indices.size(); i += stride) {
      if (indices[i] < 0)
        continue;
      map.push_back(i);
    }
    return map;
  }

  std::vector<index_t> get_edge_map() {
    return get_range_map(__corners_next, 2);
  }

  std::vector<index_t> get_vert_range() const {
    std::vector<index_t> range;
    range.reserve(corner_count());
    // replace this with some c++isms
    for (int i = 0; i < corner_count(); i++) {
      range.push_back(i);
    }

    return range;
  }

  void debug() {
    for (int i = 0; i < __corners_next.size(); i++) {
      if (__corners_next[i] == -1)
        continue;
      vec3 c1 = __x[__corners_next[i]];
      vec3 c0 = __x[i];
      vec3 v = __v[i];
      // std::cout << i << " " << __l0[i] << " " << (c1 - c0).norm() <<
      // std::endl;
      gg::geometry_logger::point(c0, vec4(0.0, 1.0, 0.0, 1.0));
      gg::geometry_logger::line(c0, c1, vec4(0.0, 0.7, 1.0, 1.0));
      gg::geometry_logger::line(c0, c0 + v, vec4(0.8, 0.3, 0.0, 1.0));
      quat u = __u[i];
#if 0
      vec3 d0 = u * vec3(1, 0, 0);
      vec3 d1 = u * vec3(0, 1, 0);
      vec3 d2 = u * vec3(0, 0, 1);

      gg::geometry_logger::line(c0, c0 + 0.03 * d0, vec4(1.0, 0.0, 0.0, 1.0));
      gg::geometry_logger::line(c0, c0 + 0.03 * d1, vec4(0.0, 1.0, 0.0, 1.0));
      gg::geometry_logger::line(c0, c0 + 0.03 * d2, vec4(0.0, 0.0, 1.0, 1.0));
#endif
    }
  }

  real _r = 0.2;
  real _lmax = 0.0;
  std::vector<index_t> __corners_next;
  std::vector<index_t> __corners_prev;

  std::vector<vec3> __M;
  std::vector<vec4> __J;

  std::vector<vec3> __x;
  std::vector<vec3> __v; // velocity;
  std::vector<quat> __u;
  std::vector<quat> __o; // omega

  std::vector<real> __l0;
};
} // namespace rod
} // namespace asawa
} // namespace gaudi
#endif