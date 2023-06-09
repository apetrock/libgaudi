
#ifndef __HEP_SIM_BLOCKS__
#define __HEP_SIM_BLOCKS__

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cxxabi.h>
#include <execinfo.h>
#include <iostream>
#include <memory.h>
#include <numeric>
#include <ostream>
#include <stdio.h>
#include <type_traits>
#include <vector>
#include <zlib.h>

#include "gaudi/common.h"

namespace gaudi {
namespace hepworth {

template <class T> T get_elem(index_t ii, const vecX &q, index_t offset){};

template <> vec3 get_elem<vec3>(index_t ii, const vecX &q, index_t offset) {
  return q.block(offset + 3 * ii, 0, 3, 1);
}

template <> quat get_elem<quat>(index_t ii, const vecX &q, index_t offset) {
  return quat(q.block(offset + 4 * ii, 0, 4, 1).data());
}

class sim_block {
public:
  typedef std::shared_ptr<sim_block> ptr;
  static ptr create() { return std::make_shared<sim_block>(); }

  sim_block() {}
  virtual ~sim_block() {}

  // vecX x = to(x_);
  // vecX q = concat(x0, u0);
  // from(u_, u);
  vec3 get_vec3(index_t ii, const vecX &q) {
    return get_elem<vec3>(ii, q, _offset);
  }
  quat get_quat(index_t ii, const vecX &q) {
    return get_elem<quat>(ii, q, _offset);
  }

  template <int N, typename T> void map_to_x(std::vector<T> &x_, vecX &q) {
    _offset = q.size();
    std::cout << __PRETTY_FUNCTION__ << ", offset: " << _offset << std::endl;
    vecX x = to<N, T>(x_);
    q = concat(q, x);
  }

  template <class T> void map_from_x(const vecX &q, std::vector<T> &x) {
    for (index_t ii = 0; ii < x.size(); ii++)
      x[ii] = get_elem<T>(ii, q, _offset);
  }

  virtual index_t get_offset_idx(index_t ii) { return -1; };
  virtual void map_to_x(vecX &q){};
  virtual void map_mass(vecX &q){};

  virtual void map_from_x(const vecX &q, const real &h, const real &damp){};

  virtual void integrate_inertia(const real &h) {}
  virtual void update_inertia(const real &h) {}
  index_t _offset;
};

class vec3_block : public sim_block {
public:
  typedef std::shared_ptr<vec3_block> ptr;
  static ptr create(std::vector<vec3> &M, std::vector<vec3> &x,
                    std::vector<vec3> &v, std::vector<vec3> &f) {
    return std::make_shared<vec3_block>(M, x, v, f);
  }

  vec3_block(std::vector<vec3> &M, std::vector<vec3> &x, std::vector<vec3> &v,
             std::vector<vec3> &f)
      : _M(M), _x(x), _v(v), _f(f) {}
  virtual ~vec3_block() {}

  virtual void map_to_x(vecX &q) { sim_block::map_to_x<3, vec3>(_x, q); }
  virtual void map_mass(vecX &q) { sim_block::map_to_x<3, vec3>(_M, q); }

  virtual void map_from_x(const vecX &q, const real &h, const real &damp) {
    std::vector<vec3> x1(_x.size());
    std::vector<vec3> v(_x.size());
    sim_block::map_from_x(q, x1);

    for (int i = 0; i < _x.size(); i++)
      v[i] = (1.0 - damp) / h * (x1[i] - _x[i]);

    _x = x1;
    _v = v;
  }

  virtual void integrate_inertia(const real &h) {
    for (size_t i = 0; i < _v.size(); i++) {
      _x[i] += h * _v[i] + h * h * _f[i];
    }
  }

  virtual index_t get_offset_idx(index_t ii) { return _offset + 3 * ii; };

  std::vector<vec3> &_x;
  std::vector<vec3> &_v;
  std::vector<vec3> _f;
  std::vector<vec3> &_M;
};

class quat_block : public sim_block {
public:
  typedef std::shared_ptr<quat_block> ptr;
  static ptr create(std::vector<vec4> &J, std::vector<quat> &u,
                    std::vector<quat> &o) {
    return std::make_shared<quat_block>(J, u, o);
  }

  quat_block(std::vector<vec4> &J, std::vector<quat> &u, std::vector<quat> &o)
      : _J(J), _u(u), _o(o) {}
  virtual ~quat_block() {}

  virtual void map_to_x(vecX &q) { sim_block::map_to_x<4, quat>(_u, q); }
  virtual void map_mass(vecX &q) { sim_block::map_to_x<4, vec4>(_J, q); }

  virtual void map_from_x(const vecX &q, const real &h, const real &damp) {
    sim_block::map_from_x(q, _u);

    std::vector<quat> u1(_u.size());
    std::vector<quat> o(_o.size());
    sim_block::map_from_x(q, u1);

    for (int i = 0; i < _u.size(); i++) {
      o[i] = _u[i].conjugate() * u1[i];
      o[i].coeffs() *= 2.0 * (1.0 - damp) / h;
    }

    _u = u1;
    _o = o;
  }

  virtual void integrate_inertia(const real &h) {
    for (size_t i = 0; i < _o.size(); i++) {
      quat o = _o[i];
      quat u = _u[i];

      // real J = J;
      quat sO = o; // torque terms + h / J
      quat su = u;
      su.coeffs() += 0.5 * h * (u * sO).coeffs();
      _u[i] = su;
      _u[i].normalize();
    }
  }

  virtual index_t get_offset_idx(index_t ii) { return _offset + 4 * ii; };
  std::vector<quat> &_u;
  std::vector<quat> &_o;
  std::vector<vec4> &_J;
};

} // namespace hepworth
} // namespace gaudi
#endif