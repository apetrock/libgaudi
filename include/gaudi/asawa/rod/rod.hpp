

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "Eigen/src/Geometry/AngleAxis.h"

#include "gaudi/common.h"
#include "gaudi/vec_addendum.h"
#include "gaudi/logger.hpp"

#include "gaudi/asawa/datums.hpp"
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

  index_t insert_datum(datum_ptr d) {
    const index_t id = __channels.insert(d);
    if (d)
      d->slot_id = id;
    return id;
  }
  datum_ptr &get_datum(index_t i) { return __channels.get(i); }
  const datum_ptr &const_get_datum(index_t i) const {
    return __channels.get(i);
  }
  std::vector<datum_ptr> &get_data() { return __channels.pool(); }
  const std::vector<datum_ptr> &get_data() const { return __channels.pool(); }
  GeometryChannels &channels() { return __channels; }
  const GeometryChannels &channels() const { return __channels; }

  void append_loop(const std::vector<vec3> &verts) { insert_strand(verts, true); }

  void insert_strand(const std::vector<vec3> &verts, bool loop);

  real length(CornerId i) const {
    if (__corners_next[i] == -1)
      return 0.0;
    index_t j = __corners_next[i];
    vec3 q0 = x()[i];
    vec3 q1 = x()[j];
    real l = (q1 - q0).norm();
    return l;
  }

  vec3 dir(CornerId i) const {
    auto idx = consec(i);
    vec3 cp = x()[idx[0]];
    vec3 c0 = x()[idx[1]];
    vec3 cn = x()[idx[2]];
    return cn - cp;
  }

  void update_lengths() {
    auto &L = l0();
    if (L.size() < corner_count())
      L.resize(corner_count(), 0.0);
    _lmax = 0.0;
    for (int i = 0; i < static_cast<int>(__corners_next.size()); i++) {
      if (__corners_next[i] == -1)
        continue;
      index_t j = __corners_next[i];
      vec3 q0 = x()[i];
      vec3 q1 = x()[j];
      L[i] = (q1 - q0).norm();
    }
  }

  void _init_params();

  void init_t() {
    auto &T = t();
    if (T.size() < corner_count())
      T.resize(corner_count(), 0.0);
    real total_length = 0.0;
    for (int i = 0; i < static_cast<int>(__corners_next.size()); i++) {
      if (__corners_next[i] == -1)
        continue;
      index_t j = __corners_next[i];
      vec3 q0 = x()[i];
      vec3 q1 = x()[j];
      T[i] = total_length;
      total_length += (q1 - q0).norm();
    }
    if (total_length > 0.0) {
      for (int i = 0; i < static_cast<int>(T.size()); i++)
        T[i] /= total_length;
    }
  }

  void update_mass() {
    _lmax = 0.0;
    auto &Mass = M();
    auto &Inertia = J();
    auto &L = l0();
    for (int i = 0; i < static_cast<int>(__corners_next.size()); i++) {
      CornerId ci = corner_id(i);
      real mass = 0.0;
      real rho = 1.0;

      if (prev(ci) == corner_id(-1)) {
        if (next(ci) != corner_id(-1)) {
          real l_right = L[i];
          mass = 0.5 * l_right;
        }
      } else if (next(ci) == corner_id(-1)) {
        real l_left = L[prev(ci)];
        mass = 0.5 * l_left;
      } else {
        real l_left = L[prev(ci)];
        real l_right = L[i];
        mass = 0.5 * (l_left + l_right);
      }

      real m = M_PI * _r * _r * mass;
      const real r2 = _r * _r;
      const real r4 = r2 * r2;
      const real J1 = M_PI * r4 / 4.0;
      const real J3 = M_PI * r4 / 2.0;
      const real lrho = mass * rho;

      Mass[i][0] = m;
      Mass[i][1] = m;
      Mass[i][2] = m;

      Inertia[i][0] = lrho * J1;
      Inertia[i][1] = lrho * J1;
      Inertia[i][2] = lrho * J3;
      Inertia[i][3] = 0.0;

      if (next(ci) != corner_id(-1)) {
        _lmax += L[i];
      }
    }

    _lmax /= real(__corners_next.size());
  }

  real get_total_volume() {
    real v = 0.0;
    for (auto m : M())
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
    vec3 cp = x()[idx[0]];
    vec3 c0 = x()[idx[1]];
    vec3 cn = x()[idx[2]];

    vec3 t0 = c0 - cp, t1 = cn - c0;
    t0.normalize();
    t1.normalize();

    return quat::FromTwoVectors(t0, t1);
  }

  mat3 _get_frenet_mat(CornerId i) {
    auto idx = consec(i);

    vec3 cp = x()[idx[0]];
    vec3 c0 = x()[idx[1]];
    vec3 cn = x()[idx[2]];

    vec3 dn0 = cn - c0;
    vec3 d0p = c0 - cp;

    vec3 T = (dn0 + d0p);
    if (T.norm() < 1e-14) {
      T = (std::abs(dn0.norm()) > 1e-14) ? dn0 : d0p;
    }
    T.normalize();

    vec3 B = (dn0 - d0p).cross(T);
    if (B.norm() < 1e-10) {
      // Nearly straight: pick any normal orthogonal to T.
      vec3 arb = (std::abs(T[0]) < 0.9) ? vec3::UnitX() : vec3::UnitY();
      B = T.cross(arb);
    }
    B.normalize();
    vec3 N = T.cross(B).normalized();

    mat3 F;
    F.col(0) = B;
    F.col(1) = N;
    F.col(2) = T;
    return F;
  }

  quat _get_frenet(CornerId i) { return quat(_get_frenet_mat(i)).normalized(); }

  quat _calc_frame(CornerId i) {
    auto idx = consec(i);
    vec3 cp = x()[idx[0]];
    vec3 c0 = x()[idx[1]];
    vec3 cn = x()[idx[2]];

    vec3 d0 = cp - cn;

    quat q(Eigen::AngleAxisd(2.0 * M_PI, d0));

    q.normalize();

    return q;
  }

  void _update_frames(const std::vector<vec3> &N_vec) {
    u().clear();
    u().resize(__corners_next.size(), quat::Identity());
    for (int i = 0; i < static_cast<int>(__corners_next.size()); i++) {
      if (__corners_next[i] == -1)
        continue;
      auto idx = consec(corner_id(i));

      vec3 cp = x()[idx[0]];
      vec3 c0 = x()[idx[1]];
      vec3 cn = x()[idx[2]];

      vec3 dn0 = cn - c0;
      vec3 d0p = c0 - cp;

      vec3 T = (dn0 + d0p);
      if (T.norm() < 1e-14) {
        T = (dn0.norm() > 1e-14) ? dn0 : d0p;
      }
      T.normalize();

      vec3 N = N_vec[i];
      if (!N.array().isFinite().all() || N.norm() < 1e-10) {
        N = _get_frenet_mat(corner_id(i)).col(1);
      }
      N = (N - N.dot(T) * T);
      if (N.norm() < 1e-10) {
        N = _get_frenet_mat(corner_id(i)).col(1);
      }
      N.normalize();
      vec3 B = N.cross(T).normalized();

      mat3 F;
      F.col(0) = B;
      F.col(1) = N;
      F.col(2) = T;

      u()[i] = quat(F).normalized();
      if (u()[i].coeffs().hasNaN()) {
        u()[i] = _get_frenet(corner_id(i));
      }
    }

    o().assign(__corners_next.size(), quat(0.0, 0.0, 0.0, 0.0));
  }

  /// Bishop / parallel-transport frames along each connected component.
  /// Walks `next`/`prev` (not corner index order) so multi-strand rods work.
  void _update_frames() {
    const int n = static_cast<int>(__corners_next.size());
    std::vector<vec3> N_vec(static_cast<size_t>(n), vec3::Zero());
    std::vector<char> visited(static_cast<size_t>(n), 0);

    for (int seed = 0; seed < n; ++seed) {
      if (__corners_next[seed] < 0 || visited[static_cast<size_t>(seed)]) {
        continue;
      }

      // Prefer an open endpoint; otherwise start at seed (closed loop).
      CornerId start = corner_id(seed);
      CornerId s = start;
      for (;;) {
        const CornerId p = prev(s);
        if (p < corner_id(0) || p == start) {
          break;
        }
        s = p;
      }

      vec3 N0 = _get_frenet_mat(s).col(1);
      if (!N0.array().isFinite().all() || N0.norm() < 1e-10) {
        auto idx = consec(s);
        vec3 T = (x()[idx[2]] - x()[idx[0]]);
        if (T.norm() < 1e-14) {
          T = vec3::UnitX();
        } else {
          T.normalize();
        }
        vec3 arb = (std::abs(T[0]) < 0.9) ? vec3::UnitX() : vec3::UnitY();
        N0 = T.cross(arb).normalized();
      } else {
        N0.normalize();
      }

      CornerId i = s;
      do {
        const int ii = static_cast<int>(i);
        visited[static_cast<size_t>(ii)] = 1;
        N_vec[static_cast<size_t>(ii)] = N0;

        const CornerId j = next(i);
        if (j < corner_id(0)) {
          break;
        }

        auto idx = consec(i);
        vec3 d0p = x()[idx[1]] - x()[idx[0]];
        vec3 dn0 = x()[idx[2]] - x()[idx[1]];
        if (d0p.norm() > 1e-14 && dn0.norm() > 1e-14) {
          const quat qi =
              quat::FromTwoVectors(d0p.normalized(), dn0.normalized());
          if (!qi.coeffs().hasNaN()) {
            N0 = qi * N0;
          }
        }
        i = j;
      } while (i != s && !visited[static_cast<size_t>(i)]);
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
      vec3 q0 = x()[ci];
      vec3 q1 = x()[j];
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

  CornerId insert_edge();

  std::vector<vec3> &corner_verts() { return x(); }
  const std::vector<vec3> &corner_verts() const { return x(); }

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
  // Const-safe accessors (datum pool). Lazy-install conventional channels
  // until a proper factory owns init of positions / frames / etc.
  void ensure_default_datums();
  const std::vector<vec3> &x() const;
  const std::vector<vec3> &v() const;
  const std::vector<quat> &u() const;
  const std::vector<quat> &o() const;
  const std::vector<vec3> &M() const;
  const std::vector<vec4> &J() const;
  const std::vector<real> &l0() const;
  const std::vector<real> &t() const;
  const std::vector<real> &rest_l() const;

  std::vector<vec3> &x();
  std::vector<vec3> &v();
  std::vector<quat> &u();
  std::vector<quat> &o();
  std::vector<vec3> &M();
  std::vector<vec4> &J();
  std::vector<real> &l0();
  std::vector<real> &t();
  std::vector<real> &rest_l();

  std::vector<vec3> xc() {
    std::vector<vec3> xc(x());
    std::vector<CornerId> rverts = get_vert_range();
    for (CornerId i : rverts) {
      auto cons = this->consec(i);
      if (cons[2] < 0)
        continue;
      xc[i] = 0.5 * (x()[cons[1]] + x()[cons[2]]);
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
    return vert_infill(x(), n_interp, false);
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
    std::vector<vec3> d(x().size());
    for (int i = 0; i < d.size(); i++) {
      d[i] = this->dir(corner_id(i));
    }
    return d;
  }

  std::vector<vec3> Ni(const vec3 &N0) {
    std::vector<vec3> N(x());
    for (int i = 0; i < N.size(); i++) {
      N[i] = u()[i] * N0;
    }
    return N;
  }

  std::vector<vec3> N0c() { return Ni(vec3(1.0, 0.0, 0.0)); }
  std::vector<vec3> N1c() { return Ni(vec3(0.0, 1.0, 0.0)); }
  std::vector<vec3> N2c() { return Ni(vec3(0.0, 0.0, 1.0)); }
  std::vector<vec3> N0() { return vert_val(N0c()); }
  std::vector<vec3> N1() { return vert_val(N1c()); }
  std::vector<vec3> N2() { return vert_val(N2c()); }

  /// Compact live corners. \c old_to_new[old] = new index, or -1 if dropped.
  void pack() {
    std::vector<index_t> old_to_new(__corners_next.size(), -1);
    int ii = 0;
    for (int i = 0; i < static_cast<int>(__corners_next.size()); i++) {
      if (__corners_next[i] < 0 && __corners_prev[i] == -1)
        continue;
      old_to_new[i] = ii++;
    }
    // permute[new] = old — datum::permute expects this orientation.
    std::vector<index_t> permute(static_cast<size_t>(ii), -1);
    for (int i = 0; i < static_cast<int>(old_to_new.size()); i++) {
      if (old_to_new[i] >= 0)
        permute[static_cast<size_t>(old_to_new[i])] = i;
    }

    std::vector<index_t> new_next;
    std::vector<index_t> new_prev;
    new_next.reserve(static_cast<size_t>(ii));
    new_prev.reserve(static_cast<size_t>(ii));
    for (int i = 0; i < static_cast<int>(old_to_new.size()); i++) {
      if (old_to_new[i] < 0)
        continue;
      if (__corners_next[i] > -1)
        new_next.push_back(old_to_new[__corners_next[i]]);
      else
        new_next.push_back(-1);
      if (__corners_prev[i] > -1)
        new_prev.push_back(old_to_new[__corners_prev[i]]);
      else
        new_prev.push_back(-1);
    }

    __corners_next = new_next;
    __corners_prev = new_prev;
    for (auto &d : __channels.pool()) {
      d->permute(permute);
      d->resize(static_cast<size_t>(ii));
    }
  }

  void debug() {
    for (int i = 0; i < __corners_next.size(); i++) {
      if (__corners_next[i] == -1)
        continue;
      vec3 c1 = x()[__corners_next[i]];
      vec3 c0 = x()[i];
      vec3 vi = v()[i];
      (void)vi;
      // std::cout << i << " " << l0()[i] << " " << (c1 - c0).norm() <<
      // std::endl;
      //geometry_logger::point(c0, vec4(0.0, 1.0, 0.0, 1.0));
      geometry_logger::line(c0, c1, vec4(0.0, 0.7, 1.0, 1.0));
      //geometry_logger::line(c0, c0 + v, vec4(0.8, 0.3, 0.0, 1.0));
      quat ui = u()[i];
#if 0
      vec3 d0 = ui * vec3(1, 0, 0);
      vec3 d1 = ui * vec3(0, 1, 0);
      vec3 d2 = ui * vec3(0, 0, 1);

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
  GeometryChannels __channels;
};
} // namespace rod
} // namespace asawa
} // namespace gaudi

#include "gaudi/asawa/rod/rod_datums.hpp"

namespace gaudi {
namespace asawa {
namespace rod {

inline CornerId rod::insert_edge() {
  __corners_next.push_back(-1);
  __corners_prev.push_back(-1);
  // Pool may still be empty (lazy via x()/…); only grow existing channels.
  for (auto &d : __channels.pool())
    d->alloc(1);
  return corner_id(static_cast<int>(__corners_next.size() - 1));
}

inline void rod::insert_strand(const std::vector<vec3> &verts, bool loop) {
  int vN = static_cast<int>(verts.size());
  int cbegin = static_cast<int>(__corners_prev.size());

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

  if (!__channels.empty()) {
    for (auto &d : __channels.pool())
      d->alloc(static_cast<size_t>(vN));
  }

  auto &X = x(); // lazy-installs default pool if needed
  for (int i = 0; i < vN; i++)
    X[cbegin + i] = verts[static_cast<size_t>(i)];

  _init_params();
}

inline void rod::_init_params() {
  ensure_default_datums();
  _lmax = 0.0;
  init_t();
  update_lengths();
  update_mass();
  _update_frames();
}

inline void rod::ensure_default_datums() {
  if (__channels.empty())
    install_default_datums(*this);
}

inline std::vector<vec3> &rod::x() {
  ensure_default_datums();
  return rod_vec3(*this, datum_id::x);
}
inline const std::vector<vec3> &rod::x() const {
  const_cast<rod *>(this)->ensure_default_datums();
  return rod_vec3(*this, datum_id::x);
}
inline std::vector<vec3> &rod::v() {
  ensure_default_datums();
  return rod_vec3(*this, datum_id::v);
}
inline const std::vector<vec3> &rod::v() const {
  const_cast<rod *>(this)->ensure_default_datums();
  return rod_vec3(*this, datum_id::v);
}
inline std::vector<quat> &rod::u() {
  ensure_default_datums();
  return rod_quat(*this, datum_id::u);
}
inline const std::vector<quat> &rod::u() const {
  const_cast<rod *>(this)->ensure_default_datums();
  return rod_quat(*this, datum_id::u);
}
inline std::vector<quat> &rod::o() {
  ensure_default_datums();
  return rod_quat(*this, datum_id::o);
}
inline const std::vector<quat> &rod::o() const {
  const_cast<rod *>(this)->ensure_default_datums();
  return rod_quat(*this, datum_id::o);
}
inline std::vector<vec3> &rod::M() {
  ensure_default_datums();
  return rod_vec3(*this, datum_id::M);
}
inline const std::vector<vec3> &rod::M() const {
  const_cast<rod *>(this)->ensure_default_datums();
  return rod_vec3(*this, datum_id::M);
}
inline std::vector<vec4> &rod::J() {
  ensure_default_datums();
  return rod_vec4(*this, datum_id::J);
}
inline const std::vector<vec4> &rod::J() const {
  const_cast<rod *>(this)->ensure_default_datums();
  return rod_vec4(*this, datum_id::J);
}
inline std::vector<real> &rod::l0() {
  ensure_default_datums();
  return rod_real(*this, datum_id::l0);
}
inline const std::vector<real> &rod::l0() const {
  const_cast<rod *>(this)->ensure_default_datums();
  return rod_real(*this, datum_id::l0);
}
inline std::vector<real> &rod::t() {
  ensure_default_datums();
  return rod_real(*this, datum_id::t);
}
inline const std::vector<real> &rod::t() const {
  const_cast<rod *>(this)->ensure_default_datums();
  return rod_real(*this, datum_id::t);
}
inline std::vector<real> &rod::rest_l() {
  ensure_default_datums();
  return rod_real(*this, datum_id::rest_l);
}
inline const std::vector<real> &rod::rest_l() const {
  const_cast<rod *>(this)->ensure_default_datums();
  return rod_real(*this, datum_id::rest_l);
}

/// Pack rod topology together with an explicit channel pool (usually R.channels()).
inline void pack(rod &R, GeometryChannels &channels) {
  (void)channels;
  R.pack();
}

} // namespace rod
} // namespace asawa
} // namespace gaudi

#endif
