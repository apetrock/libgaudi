#ifndef __HEP_SIM_BLOCKS__
#define __HEP_SIM_BLOCKS__

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cxxabi.h>

#include <iostream>
#include <memory>
#include <numeric>
#include <ostream>
#include <stdio.h>
#include <type_traits>
#include <vector>

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
    vecX x = to<N, T>(x_);
    q = concat(q, x);
    if (q.hasNaN()) {
      std::cout << "NAN:" << __PRETTY_FUNCTION__ << std::endl;
    }
  }

  template <class T> void map_from_x(const vecX &q, std::vector<T> &x) {
    for (index_t ii = 0; ii < x.size(); ii++)
      x[ii] = get_elem<T>(ii, q, _offset);
  }

  virtual index_t get_offset_idx(index_t ii) const { return -1; };
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
    sim_block::map_from_x(q, x1);

    // Paper / Bouaziz: v^{t+1} = (x^{t+1} - x^t) / h  (not vs the prediction s).
    const std::vector<vec3> &x0 = _x_prev.size() == _x.size() ? _x_prev : _x;
    for (int i = 0; i < _x.size(); i++)
      _v[i] = (1.0 - damp) / h * (x1[i] - x0[i]);

    _x = x1;
  }

  virtual void integrate_inertia(const real &h) {
    _x_prev = _x;
    for (size_t i = 0; i < _v.size(); i++) {
      // s = x + h v + h² a  (f stored as acceleration, matching existing force path)
      _x[i] += h * _v[i] + h * h * _f[i];
    }
  }

  virtual index_t get_offset_idx(index_t ii) const { return _offset + 3 * ii; };

  std::vector<vec3> &_x;
  std::vector<vec3> &_v;
  std::vector<vec3> &_f;
  std::vector<vec3> &_M;
  std::vector<vec3> _x_prev;
};

class quat_block : public sim_block {
public:
  typedef std::shared_ptr<quat_block> ptr;
  static ptr create(std::vector<vec4> &J, std::vector<quat> &u,
                    std::vector<quat> &o) {
    return std::make_shared<quat_block>(J, u, o);
  }

  static ptr create(std::vector<vec4> &J, std::vector<quat> &u, std::vector<quat> &o,
                    std::vector<vec3> &torques) {
    return std::make_shared<quat_block>(J, u, o, torques);
  }

  quat_block(std::vector<vec4> &J, std::vector<quat> &u, std::vector<quat> &o)
      : _J(J), _u(u), _o(o), _torques(nullptr) {}

  quat_block(std::vector<vec4> &J, std::vector<quat> &u, std::vector<quat> &o,
             std::vector<vec3> &torques)
      : _J(J), _u(u), _o(o), _torques(&torques) {}
  virtual ~quat_block() {}

  virtual void map_to_x(vecX &q) { sim_block::map_to_x<4, quat>(_u, q); }
  virtual void map_mass(vecX &q) { sim_block::map_to_x<4, vec4>(_J, q); }

  virtual void map_from_x(const vecX &q, const real &h, const real &damp) {
    std::vector<quat> u1(_u.size());
    sim_block::map_from_x(q, u1);

    // Soler et al. Alg.1 line 12: ω^{t+1} = (2/h) Im(ū^t ◦ u^{t+1})
    // Must use u^t saved before integrate — not the prediction s_u, and not u^{t+1}.
    const std::vector<quat> &u0 = _u_prev.size() == _u.size() ? _u_prev : _u;
    for (int i = 0; i < static_cast<int>(_u.size()); i++) {
      quat dq = u0[i].conjugate() * u1[i];
      if (dq.w() < 0.0)
        dq.coeffs() *= -1.0; // short arc
      dq.coeffs() *= 2.0 * (1.0 - damp) / h;
      _o[i] = dq;
    }
    _u = u1;
  }

  virtual void integrate_inertia(const real &h) {
    // Soler et al. Alg.1 lines 3–4:
    //   s_ω = ω + h J⁻¹[τ − ω × (Jω)]
    //   s_u = u + (h/2)(u ◦ s_ω)
    // τ is world-space; convert to body for Newton–Euler.
    _u_prev = _u;
    for (size_t i = 0; i < _o.size(); i++) {
      vec3 omega(_o[i].x(), _o[i].y(), _o[i].z());
      const vec4 &Ji = _J[i];
      const vec3 Jw(Ji[0] * omega[0], Ji[1] * omega[1], Ji[2] * omega[2]);

      vec3 tau_body = vec3::Zero();
      if (_torques && i < _torques->size())
        tau_body = _u[i].inverse() * (*_torques)[i];

      const vec3 torque_net = tau_body - omega.cross(Jw);
      for (int k = 0; k < 3; ++k) {
        const real Jk = std::max(std::abs(Ji[k]), real(1e-12));
        omega[k] += h * torque_net[k] / Jk;
      }

      const quat sO(0.0, omega[0], omega[1], omega[2]);
      quat su = _u[i];
      su.coeffs() += 0.5 * h * (_u[i] * sO).coeffs();
      _u[i] = su;
      _u[i].normalize();
    }
  }

  virtual index_t get_offset_idx(index_t ii) const { return _offset + 4 * ii; };
  std::vector<quat> &_u;
  std::vector<quat> &_o;
  std::vector<vec4> &_J;
  std::vector<vec3> *_torques = nullptr;
  std::vector<quat> _u_prev;
};

} // namespace hepworth
} // namespace gaudi
#endif