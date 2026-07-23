

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "Eigen/src/Geometry/AngleAxis.h"

#include "gaudi/common.h"
#include "gaudi/vec_addendum.h"
#include "gaudi/logger.hpp"

#include "gaudi/asawa/rod/rod_id.hpp"

// #include "subdivide.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <type_traits>
#include <vector>
#include "gaudi/geometry_logger.hpp"

#ifndef __ASAWA_ROD__
#define __ASAWA_ROD__

namespace gaudi {
namespace asawa {
namespace rod {

using consec_t = std::array<CornerId, 3>;

class rod {
public:
  typedef std::shared_ptr<rod> ptr;

  static ptr create() {
    ptr R = std::make_shared<rod>();
    return R;
  }
  static ptr create(const std::vector<vec3> &verts, bool loop = true) {
    ptr R = std::make_shared<rod>();
    if (loop) {
      R->append_loop(verts);
    } else {
      R->insert_strand(verts, false);
    }
    return R;
  }

  rod() {}

  void append_loop(const std::vector<vec3> &verts) { insert_strand(verts, true); }

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

  real length(CornerId i) const {
    if (__corners_next[i] == -1)
      return 0.0;
    index_t j = __corners_next[i];
    vec3 q0 = __x[i];
    vec3 q1 = __x[j];
    real l = (q1 - q0).norm();
    return l;
  }

  vec3 dir(CornerId i) const {
    auto idx = consec(i);
    vec3 cp = __x[idx[0]];
    vec3 c0 = __x[idx[1]];
    vec3 cn = __x[idx[2]];
    return cn - cp;
  }

  void update_lengths() {
    __l0.resize(__corners_next.size(), 0.0);
    __v.resize(__corners_next.size(), vec3::Zero());
    __M.resize(__corners_next.size(), vec3::Zero());
    __J.resize(__corners_next.size(), vec4::Zero());
    __o.resize(__corners_next.size(), quat::Identity());
    _lmax = 0.0;
    for (int i = 0; i < __corners_next.size(); i++) {
      if (__corners_next[i] == -1)
        continue;
      index_t j = __corners_next[i];
      vec3 q0 = __x[i];
      vec3 q1 = __x[j];
      real l = (q1 - q0).norm();
      __l0[i] = l;
    }
  }

  void _init_params() {
    __l0.resize(__corners_next.size(), 0.0);
    __v.resize(__corners_next.size(), vec3::Zero());
    __M.resize(__corners_next.size(), vec3::Zero());
    __J.resize(__corners_next.size(), vec4::Zero());
    __o.resize(__corners_next.size(), quat::Identity());
    _lmax = 0.0;
    init_t();

    update_lengths();
    update_mass();
    _update_frames();
  }

  void init_t() {

    __t.resize(__corners_next.size(), 0.0);
    real total_length = 0.0;
    for (int i = 0; i < __corners_next.size(); i++) {
      if (__corners_next[i] == -1)
        continue;
      index_t j = __corners_next[i];
      vec3 q0 = __x[i];
      vec3 q1 = __x[j];
      __t[i] = total_length;
      total_length += (q1 - q0).norm();
    }

    for (int i = 0; i < __t.size(); i++) {
      __t[i] /= total_length;
    }
  }

  void update_mass() {

    _lmax = 0.0;
    for (int i = 0; i < __corners_next.size(); i++) {
      CornerId ci = corner_id(i);
      real mass = 0.0;
      real rho = 1.0;
      
      if (prev(ci) == corner_id(-1)) {
        if (next(ci) != corner_id(-1)) {
          real l_right = __l0[i];
          mass = 0.5 * l_right;
        }
      } else if (next(ci) == corner_id(-1)) {
        real l_left = __l0[prev(ci)];
        mass = 0.5 * l_left;
      } else {
        real l_left = __l0[prev(ci)];
        real l_right = __l0[i];
        mass = 0.5 * (l_left + l_right);
      }
      
      real M = M_PI * _r * _r * mass;
      // Soler et al.: J'_n = l ρ diag(J1,J2,J3) with circular
      // J1=J2=πr⁴/4, J3=πr⁴/2. Stored as Eigen quat coeffs (x,y,z,w) with w=0.
      const real r2 = _r * _r;
      const real r4 = r2 * r2;
      const real J1 = M_PI * r4 / 4.0;
      const real J3 = M_PI * r4 / 2.0;
      const real lrho = mass * rho;

      __M[i][0] = M;
      __M[i][1] = M;
      __M[i][2] = M;

      __J[i][0] = lrho * J1;
      __J[i][1] = lrho * J1;
      __J[i][2] = lrho * J3;
      __J[i][3] = 0.0;
      
      if (next(ci) != corner_id(-1)) {
        _lmax += __l0[i];
      }
    }
    
    _lmax /= real(__corners_next.size());
  }

  real get_total_volume() {
    real v = 0.0;
    for (auto m : __M)
      v += m[0];
    return v;
  }

  consec_t consec(CornerId i) const {
    CornerId in = next(i);

    CornerId ip = prev(i);
    if (ip == corner_id(-1)) {
      return {i, i, in};
    } else if (in == corner_id(-1)) {
      return {ip, i, i};
    }

    return {ip, i, in};
  };

  quat get_rotation(CornerId i) {
    auto idx = consec(i);
    vec3 cp = __x[idx[0]];
    vec3 c0 = __x[idx[1]];
    vec3 cn = __x[idx[2]];

    vec3 t0 = c0 - cp, t1 = cn - c0;
    t0.normalize();
    t1.normalize();

    return quat::FromTwoVectors(t0, t1);
  }

  mat3 _get_frenet_mat(CornerId i) {
    auto idx = consec(i);

    vec3 cp = __x[idx[0]];
    vec3 c0 = __x[idx[1]];
    vec3 cn = __x[idx[2]];

    vec3 dn0 = cn - c0;
    vec3 d0p = c0 - cp;
    vec3 dd0 = dn0 - d0p;

    // vec3 T = (dn0 + d0p).normalized();          // tan
    // vec3 B = dn0.cross(d0p).normalized();       // norm
    // vec3 N = T.cross(B).normalized();           // binorm
    vec3 T = (dn0 + d0p).normalized();          // tan
    vec3 B = (dn0 - d0p).cross(T).normalized(); // norm
    vec3 N = T.cross(B).normalized();           // binorm

    mat3 F;
    F.col(0) = B;
    F.col(1) = N;
    F.col(2) = T;
    return F;
  }

  quat _get_frenet(CornerId i) { return quat(_get_frenet_mat(i)).normalized(); }

  quat _calc_frame(CornerId i) {
    auto idx = consec(i);
    vec3 cp = __x[idx[0]];
    vec3 c0 = __x[idx[1]];
    vec3 cn = __x[idx[2]];

    vec3 d0 = cp - cn;

    quat q(Eigen::AngleAxisd(2.0 * M_PI, d0));

    q.normalize();

    return q;
  }

  void _update_frames(const std::vector<vec3> &N_vec) {
    __u.clear();
    __u.resize(__corners_next.size(), quat::Identity());
    for (int i = 0; i < __corners_next.size(); i++) {
      if (__corners_next[i] == -1)
        continue;
      auto idx = consec(corner_id(i));

      vec3 cp = __x[idx[0]];
      vec3 c0 = __x[idx[1]];
      vec3 cn = __x[idx[2]];

      vec3 dn0 = cn - c0;
      vec3 d0p = c0 - cp;
      vec3 dd0 = dn0 - d0p;

      vec3 T = (dn0 + d0p).normalized();         // tan
      vec3 B = (N_vec[i]).cross(T).normalized(); // norm
      vec3 N = N_vec[i];
      mat3 F;
      F.col(0) = B;
      F.col(1) = N;
      F.col(2) = T;

      __u[i] = quat(F).normalized();

      if (__u[i].coeffs().hasNaN()) {
        std::cout << __PRETTY_FUNCTION__ << i << std::endl;
        std::cout << dn0 << " " << d0p << " " << dd0 << std::endl;
        std::cout << "F: " << F << std::endl;
        exit(0);
      }
    }

    for (int i = 0; i < __corners_next.size(); i++) {
      __o[i] = quat(0.0, 0.0, 0.0, 0.0);
    }

    // fix_frame();
  }

  void _update_frames() {
    __u.resize(__corners_next.size());
    std::vector<vec3> N_vec(__corners_next.size());
    vec3 N0 = _get_frenet_mat(corner_id(0)).col(1);
    for (int i = 0; i < __corners_next.size(); i++) {
      if (__corners_next[i] == -1)
        continue;
      auto idx = consec(corner_id(i));

      vec3 cp = __x[idx[0]];
      vec3 c0 = __x[idx[1]];
      vec3 cn = __x[idx[2]];

      vec3 dn0 = cn - c0;
      vec3 d0p = c0 - cp;
      vec3 dd0 = dn0 - d0p;
      quat qi = Eigen::Quaterniond::FromTwoVectors(d0p, dn0);
      if (qi.coeffs().hasNaN()) {
        std::cout << __PRETTY_FUNCTION__ << i << std::endl;
        std::cout << dn0 << " " << d0p << " " << dd0 << std::endl;
        std::cout << qi << std::endl;
        exit(0);
      }
      N_vec[i] = N0;
      N0 = qi * N0;
    }
    _update_frames(N_vec);
  }

  real lavg() {
    real lbar = 0.0;
    for (int i = 0; i < corner_count(); i++) {
      CornerId ci = corner_id(i);
      if (next(ci) == corner_id(-1))
        continue;
      CornerId j = next(ci);
      vec3 q0 = __x[ci];
      vec3 q1 = __x[j];
      real l = (q1 - q0).norm();
      lbar += l;
    }
    lbar /= real(corner_count());
    return lbar;
  }

  CornerId next(CornerId i) const {
    return i < 0 ? corner_id(-1) : corner_id(__corners_next[i]);
  }
  CornerId prev(CornerId i) const {
    return i < 0 ? corner_id(-1) : corner_id(__corners_prev[i]);
  }

  void set_next(CornerId id, CornerId c) { __corners_next[id] = c; }
  void set_prev(CornerId id, CornerId c) { __corners_prev[id] = c; }
  void link(CornerId c0, CornerId c1) {
    set_next(c0, c1);
    set_prev(c1, c0);
  }

  std::vector<CornerId> get_endpoints() const {
    std::vector<CornerId> endpoints;
    endpoints.reserve(corner_count());
    
    for (int i = 0; i < corner_count(); i++) {
      CornerId ci = corner_id(i);
      if (next(ci) == corner_id(-1) || prev(ci) == corner_id(-1)) {
        endpoints.push_back(ci);
      }
    }
    return endpoints;
  }

  size_t corner_count() const { return __corners_next.size(); }

  CornerId insert_edge() {
    __corners_next.push_back(-1);
    __corners_prev.push_back(-1);
    return corner_id(static_cast<int>(__corners_next.size() - 1));
  }

  std::vector<vec3> &corner_verts() { return __x; }
  const std::vector<vec3> &corner_verts() const { return __x; }

  std::vector<CornerId> get_ordered_verts() {
    std::vector<CornerId> verts;
    verts.reserve(corner_count());
    std::vector<bool> visited(corner_count(), false);

    for (int i = 0; i < corner_count(); i++) {
      CornerId ci = corner_id(i);
      if (next(ci) < 0 || visited[i])
        continue;
      CornerId j = ci;
      while (j > -1 && !visited[j]) {
        verts.push_back(j);
        visited[j] = true;
        j = next(j);
      }
    }
    return verts;
  }

  std::vector<index_t> get_range_map(const std::vector<index_t> indices,
                                     int stride) {
    std::vector<index_t> map;
    map.reserve(indices.size() / stride);
    for (int i = 0; i < indices.size(); i += stride) {
      if (indices[i] < 0)
        continue;
      map.push_back(i);
    }
    return map;
  }

  std::vector<index_t> get_edge_map() {
    return get_range_map(__corners_next, 1);
  }

  std::vector<index_t> get_vert_map() {
    return get_range_map(__corners_next, 1);
  }

  std::vector<CornerId> get_vert_range() const {
    std::vector<CornerId> range;
    range.reserve(corner_count());
    for (int i = 0; i < corner_count(); i++) {
      if (__corners_next[i] < 0)
        continue;
      range.push_back(corner_id(i));
    }

    return range;
  }

  template <typename T>
  std::vector<T> spread(index_t N, const std::vector<T> &v0) const {
    std::vector<index_t> v1(v0.size() * N);
    for (int i = 0; i < v0.size(); ++i) {
      for (int k = 0; k < N; ++k) {
        v1[i * N + k] = v0[i];
      }
    }
    return v1;
  }

  std::vector<index_t> get_edge_vert_ids() {
    std::vector<index_t> range;
    range.reserve(corner_count());
    for (int i = 0; i < corner_count(); i++) {
      CornerId ci = corner_id(i);
      if (next(ci) < 0)
        continue;
      range.push_back(i);
      range.push_back(__corners_next[i]);
    }
    return range;
  }
  // Const-safe accessors
  const std::vector<vec3> &x() const { return __x; }
  const std::vector<vec3> &v() const { return __v; }
  const std::vector<quat> &u() const { return __u; }
  const std::vector<quat> &o() const { return __o; }
  const std::vector<vec3> &M() const { return __M; }
  const std::vector<vec4> &J() const { return __J; }
  const std::vector<real> &l0() const { return __l0; }
  const std::vector<real> &t() const { return __t; }

  // Mutable accessors
  std::vector<vec3> &x() { return __x; }
  std::vector<vec3> &v() { return __v; }
  std::vector<quat> &u() { return __u; }
  std::vector<quat> &o() { return __o; }
  std::vector<vec3> &M() { return __M; }
  std::vector<vec4> &J() { return __J; }
  std::vector<real> &l0() { return __l0; }

  std::vector<real> &t() { return __t; }

  std::vector<vec3> xc() {
    std::vector<vec3> xc(__x);
    std::vector<CornerId> rverts = get_vert_range();
    for (CornerId i : rverts) {
      auto cons = this->consec(i);
      if (cons[2] < 0)
        continue;
      xc[i] = 0.5 * (__x[cons[1]] + __x[cons[2]]);
    }
    return xc;
  }

  std::vector<vec3> vert_infill(const std::vector<vec3> &x, index_t n_interp,
                                bool normalize = true) const {
    std::vector<CornerId> rverts = get_vert_range();
    std::vector<vec3> xc(n_interp * x.size());
    for (CornerId i : rverts) {
      auto cons = this->consec(i);
      int ii = cons[1];
      for (int k = 0; k < n_interp; k++) {
        if (cons[2] < 0)
          xc[n_interp * ii + k] = x[cons[1]];

        real t = real(k) / real(n_interp);
        xc[n_interp * ii + k] = va::mix(t, x[cons[1]], x[cons[2]]);
      }
    }
    return xc;
  }

  std::vector<vec3> x_infill(index_t n_interp) const {
    return vert_infill(__x, n_interp, false);
  }

  std::vector<vec3> vert_val(const std::vector<vec3> &x,
                             bool normalize = true) {
    std::vector<vec3> xp(x);
    std::vector<CornerId> rverts = get_vert_range();
    for (CornerId i : rverts) {
      auto cons = this->consec(i);
      if (cons[0] < 0)
        continue;
      if (cons[2] < 0)
        continue;
      xp[i] = 0.5 * (x[cons[0]] + x[cons[2]]);
      xp[i] = normalize ? xp[i].normalized() : xp[i];
    }
    return xp;
  }

  std::vector<vec3> dirs() {
    std::vector<vec3> d(__x.size());
    for (int i = 0; i < d.size(); i++) {
      d[i] = this->dir(corner_id(i));
    }
    return d;
  }

  std::vector<vec3> Ni(const vec3 &N0) {
    std::vector<vec3> N(__x);
    for (int i = 0; i < N.size(); i++) {
      N[i] = __u[i] * N0;
    }
    return N;
  }

  std::vector<vec3> N0c() { return Ni(vec3(1.0, 0.0, 0.0)); }
  std::vector<vec3> N1c() { return Ni(vec3(0.0, 1.0, 0.0)); }
  std::vector<vec3> N2c() { return Ni(vec3(0.0, 0.0, 1.0)); }
  std::vector<vec3> N0() { return vert_val(N0c()); }
  std::vector<vec3> N1() { return vert_val(N1c()); }
  std::vector<vec3> N2() { return vert_val(N2c()); }

  void pack() {
    std::vector<index_t> cmap(__corners_next.size(), -1);

    std::vector<vec3> M;
    std::vector<vec4> J;

    std::vector<vec3> x;
    std::vector<vec3> v; // velocity;
    std::vector<quat> u;
    std::vector<quat> o; // omega

    std::vector<real> l0;
    std::vector<real> t;

    int ii = 0;
    for (int i = 0; i < __corners_next.size(); i++) {
      if (__corners_next[i] < 0 && __corners_prev[i] == -1) {
        continue;
      }
      cmap[i] = ii++;
      M.push_back(__M[i]);
      J.push_back(__J[i]);
      x.push_back(__x[i]);
      v.push_back(__v[i]);
      u.push_back(__u[i]);
      o.push_back(__o[i]);
      t.push_back(__t[i]);
      l0.push_back(__l0[i]);
    }
    std::vector<index_t> new_next;
    std::vector<index_t> new_prev;
    for (int i = 0; i < cmap.size(); i++) {
      if (cmap[i] < 0)
        continue;
      if (__corners_next[i] > -1)
        new_next.push_back(cmap[__corners_next[i]]);
      else
        new_next.push_back(-1);

      if (__corners_prev[i] > -1)
        new_prev.push_back(cmap[__corners_prev[i]]);
      else
        new_prev.push_back(-1);
    }

    __corners_next = new_next;
    __corners_prev = new_prev;
    __M = M;
    __J = J;
    __x = x;
    __v = v;
    __u = u;
    __o = o;
    __l0 = l0;
    __t = t;
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
      //geometry_logger::point(c0, vec4(0.0, 1.0, 0.0, 1.0));
      geometry_logger::line(c0, c1, vec4(0.0, 0.7, 1.0, 1.0));
      //geometry_logger::line(c0, c0 + v, vec4(0.8, 0.3, 0.0, 1.0));
      quat u = __u[i];
#if 0
      vec3 d0 = u * vec3(1, 0, 0);
      vec3 d1 = u * vec3(0, 1, 0);
      vec3 d2 = u * vec3(0, 0, 1);

      geometry_logger::line(c0, c0 + 0.03 * d0, vec4(1.0, 0.0, 0.0, 1.0));
      geometry_logger::line(c0, c0 + 0.03 * d1, vec4(0.0, 1.0, 0.0, 1.0));
      geometry_logger::line(c0, c0 + 0.03 * d2, vec4(0.0, 0.0, 1.0, 1.0));
#endif
    }
  }
  void set_radius(real r) { _r = r; }
  real _r = 0.05;
  real _lmax = 0.0;

  std::vector<index_t> __corners_next;
  std::vector<index_t> __corners_prev;

  std::vector<real> __t;
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
