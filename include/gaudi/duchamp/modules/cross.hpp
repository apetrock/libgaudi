
#ifndef __DUCHAMP_VORTEX_MODULE__
#define __DUCHAMP_VORTEX_MODULE__

#include "gaudi/asawa/datums.hpp"
#include "gaudi/calder/shell_integrators.hpp"

#include "gaudi/common.h"
#include "gaudi/logger.hpp"
#include "module_base_shell.hpp"
#include <array>
#include <vector>

namespace gaudi {

// probably doesn't need to be a module, could go straight in test class...
/////////////////////////////////////////////////////////////

namespace duchamp {

class cross : public module_base_shell {
public:
  DEFINE_CREATE_FUNC(cross)
  cross(asawa::shell::shell::ptr M) : module_base_shell(M) {
    const std::vector<vec3> &x = asawa::const_get_vec_data(*M, 0);
    _eps = asawa::shell::avg_length(*M, x);
  };

  virtual ~cross(){};

  //     i
  // D  / \ A
  //   /   \
  // l  ---  j
  //   \   /
  // B  \ / C
  //     k

  real _f(asawa::shell::shell &M, index_t c0, const std::vector<vec3> &x) {

    index_t c1 = M.other(c0);
    index_t ci = M.prev(c0);
    index_t cj = M.next(c0);
    index_t ck = M.prev(c1);
    index_t cl = M.next(c1);
    vec3 da = asawa::shell::edge_tangent(M, cj, x);
    vec3 db = asawa::shell::edge_tangent(M, cl, x);
    vec3 dc = asawa::shell::edge_tangent(M, ck, x);
    vec3 dd = asawa::shell::edge_tangent(M, ci, x);

    real la = da.norm();
    real lb = db.norm();
    real lc = dc.norm();
    real ld = dd.norm();
    real cross = la * lb / (lc * ld);
    return cross;
  }

  std::array<vec3, 4> _dlnf(asawa::shell::shell &M, index_t c0,
                            const std::vector<vec3> &x) {
    index_t c1 = M.other(c0);
    index_t ci = M.prev(c0);
    index_t cj = M.next(c0);
    index_t ck = M.prev(c1);
    index_t cl = M.next(c1);
    vec3 da = asawa::shell::edge_tangent(M, cj, x);
    vec3 db = asawa::shell::edge_tangent(M, cl, x);
    vec3 dc = asawa::shell::edge_tangent(M, ck, x);
    vec3 dd = asawa::shell::edge_tangent(M, ci, x);
    real la = da.norm();
    real lb = db.norm();
    real lc = dc.norm();
    real ld = dd.norm();
    vec3 ga = da / la / la;
    vec3 gb = db / lb / lb;
    vec3 gc = dc / lc / lc;
    vec3 gd = dd / ld / ld;
    vec3 gi = ga - gd;
    vec3 gj = ga - gc;
    vec3 gk = gb - gc;
    vec3 gl = gb - gd;

    return {gi, gj, gk, gl};
  }

  std::array<vec3, 4> _df(asawa::shell::shell &M, index_t c0,
                          const std::vector<vec3> &x, bool invert = false) {

    real f = _f(M, c0, x);
    if (invert)
      f = -1.0 / f;
    std::array<vec3, 4> dlnf = _dlnf(M, c0, x);
    return {f * dlnf[0], //
            f * dlnf[1], //
            f * dlnf[2], //
            f * dlnf[3]};
  }

  std::vector<real> calc_cross() {
    asawa::shell::shell &M = *_M;

    std::vector<index_t> edge_range = M.get_edge_range();
    std::vector<vec3> &x = asawa::get_vec_data(M, 0);
    std::vector<real> X(M.corner_count() / 2, 0.0);

    for (auto c0 : edge_range) {
      X[c0 / 2] = _f(M, c0, x);
    }

    return X;
  }

  std::vector<vec3> calc_cross_grad() {
    asawa::shell::shell &M = *_M;

    std::vector<index_t> edge_range = M.get_edge_range();
    std::vector<vec3> &x = asawa::get_vec_data(M, 0);
    std::vector<vec3> dX(M.vert_count(), vec3::Zero());
    for (auto c0 : edge_range) {
      index_t c1 = M.other(c0);
      index_t vi = M.vert(M.prev(c0));
      index_t vj = M.vert(M.next(c0));
      index_t vk = M.vert(M.prev(c1));
      index_t vl = M.vert(M.next(c1));
      std::array<vec3, 4> df = _df(M, c0, x);
      dX[vi] += df[0];
      dX[vj] += df[1];
      dX[vk] += df[2];
      dX[vl] += df[3];
    }

    return dX;
  }

  std::vector<real> calc_objective() {
    asawa::shell::shell &M = *_M;

    std::vector<index_t> edge_range = M.get_edge_range();
    std::vector<vec3> &x = asawa::get_vec_data(M, 0);
    std::vector<real> X(M.corner_count() / 2, 0.0);
    for (auto c0 : edge_range) {
      real f = _f(M, c0, x);
      real dobj0 = 2.0 * (f - _lambda);
      real dobj1 = 2.0 * (1.0 / f - _lambda);
      bool invert = abs(dobj0) > abs(dobj1);
      real dobj = invert ? dobj1 : dobj0;
      X[c0 / 2] = pow(dobj, 2.0);
    }
    return X;
  }

  std::vector<vec3> calc_obj_grad() {
    asawa::shell::shell &M = *_M;

    std::vector<index_t> edge_range = M.get_edge_range();
    std::vector<vec3> &x = asawa::get_vec_data(M, 0);
    std::vector<vec3> dX(M.vert_count(), vec3::Zero());
    for (auto c0 : edge_range) {

      index_t c1 = M.other(c0);
      index_t vi = M.vert(M.prev(c0));
      index_t vj = M.vert(M.next(c0));
      index_t vk = M.vert(M.prev(c1));
      index_t vl = M.vert(M.next(c1));

      real f = _f(M, c0, x);
      real dobj0 = 2.0 * (f - _lambda);
      real dobj1 = 2.0 * (1.0 / f - _lambda);
      bool invert = abs(dobj0) > abs(dobj1);
      // invert = false;

      real dobj = invert ? dobj1 : dobj0;
      std::array<vec3, 4> df = _df(M, c0, x, invert);
      dX[vi] += dobj * df[0];
      dX[vj] += dobj * df[1];
      dX[vk] += dobj * df[2];
      dX[vl] += dobj * df[3];
    }
    return dX;
  }

  std::vector<vec3> get_smoothd_cross_grad() {
    asawa::shell::shell &M = *_M;
    std::vector<vec3> &x = asawa::get_vec_data(M, 0);

    std::vector<real> X = calc_objective();
    std::vector<vec3> X_grad = calc_obj_grad();
    std::vector<real> X_f = asawa::shell::edge_to_face<real>(M, x, X);
    std::vector<vec3> X_grad_f = asawa::shell::vert_to_face<vec3>(M, x, X_grad);
    real exp = 3.0;
    std::vector<vec3> s_grad_s = calder::gradient_scalar(M, x, X_f, _eps, exp);
    std::vector<vec3> s_grad_g =
        calder::smoothed_gradient(M, x, X_grad_f, _eps, exp);

    std::vector<vec3> dX(M.vert_count(), vec3::Zero());
    for (int i = 0; i < dX.size(); i++) {
      vec3 gs = s_grad_s[i];
      vec3 gg = s_grad_g[i];
      logger::line(x[i], x[i] + 1e-7 * gs, vec4(0.7, 0.2, 0.1, 0.5));
      logger::line(x[i], x[i] + 1e-7 * gg, vec4(0.2, 0.1, 0.7, 0.5));

      dX[i] = s_grad_s[i] + s_grad_g[i];
      // dX[i] = s_grad_s[i];
    }

    return dX;
  }

  virtual void step(real h) {

    asawa::shell::shell &M = *_M;
    std::vector<vec3> &x = asawa::get_vec_data(M, 0);
    std::vector<vec3> X = get_smoothd_cross_grad();
    for (int i = 0; i < X.size(); i++) {
      x[i] += 1e-7 * h * X[i];
    }
  }

  real _lambda = 0.95;
  real _eps;
};

} // namespace duchamp
} // namespace gaudi

#endif