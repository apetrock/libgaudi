
#ifndef __DUCHAMP_MORPHING_MODULE__
#define __DUCHAMP_MORPHING_MODULE__

#include "gaudi/asawa/datums.hpp"
#include "gaudi/calder/integrators.hpp"

#include "gaudi/common.h"
#include "gaudi/logger.hpp"
#include "module_base_shell.hpp"

#include "module_base.hpp"
#include <algorithm>
#include <array>
#include <vector>

namespace gaudi {

// probably doesn't need to be a module, could go straight in test class...
/////////////////////////////////////////////////////////////

namespace duchamp {

class mighty_morphin : public module_base_shell {
public:
  DEFINE_CREATE_FUNC(mighty_morphin)
  mighty_morphin(const asawa::shell::shell::ptr &M,         //
                 const std::vector<vec3> &x,                //
                 const std::vector<index_t> &face_vert_ids, //
                 real s, vec3 offset)
      : module_base_shell(M) {

    _face_vert_ids = face_vert_ids;
    _x = x;

    for (int i = 0; i < _x.size(); i++) {
      _x[i] = s * _x[i] - offset;
    }

    _face_tree = arp::T3::create(_face_vert_ids, _x, 24);
  };

  virtual ~mighty_morphin(){};

  virtual std::vector<real> calc_dist(const std::vector<vec3> &x) {
    std::vector<index_t> vert_ids(x.size(), -1);
    for (int i = 0; i < x.size(); i++) {
      vert_ids[i] = i;
    }
    std::vector<real> dists(x.size(), 0.0);
    for (int i = 0; i < x.size(); i++) {

      std::vector<index_t> nearest = arp::getNearest<1, 3>(
          i, vert_ids, x, //
          *_face_tree,    //
          std::numeric_limits<real>::infinity(), &arp::pnt_tri_min);

      index_t j = nearest.back();
      if (j < 0)
        continue;
      vec3 xi = x[i];

      vec3 x0 = _x[_face_vert_ids[3 * j + 0]];
      vec3 x1 = _x[_face_vert_ids[3 * j + 1]];
      vec3 x2 = _x[_face_vert_ids[3 * j + 2]];

      std::array<real, 4> cp = va::closest_point({x0, x1, x2}, xi);

      vec3 xT = cp[1] * x0 + cp[2] * x1 + cp[3] * x2;
      // if (i == 0)
      dists[i] = cp[0];
      // logger::line(xi, xT, vec4(0.0, 0.5, 1.0, 1.0));
    }
    return dists;
  }

  virtual void step(real h) {
    std::vector<vec3> &x = asawa::get_vec_data(*_M, 0);

    std::vector<vec3> xf = asawa::shell::face_centers(*_M, x);

    std::vector<real> winding = calder::fast_winding(_face_tree, xf, 1.0);
    // std::vector<vec3> Nv = asawa::shell::vertex_normals(*_M, x);

    real l0 = asawa::shell::avg_length(*_M, x);
    std::vector<vec3> Nf = asawa::shell::face_normals(*_M, x);
    std::vector<real> dists = calder::fast_dist(_face_tree, xf, 1.0);

    std::vector<vec3> f(xf.size(), vec3::Zero());
    std::cout << "winding: " << winding[0] << std::endl;

    // std::vector<real> dists = calc_dist(x);

    for (int i = 0; i < xf.size(); i++) {
      vec3 xi = xf[i];
      vec3 fi = dists[i] * Nf[i].normalized() * (winding[i] - 0.5);
      logger::line(xi, xi + 0.1 * fi, vec4(0.0, 1.0, 0.5, 1.0));
      f[i] += fi;
    }

    _f = calder::mls_avg(*_M, f, x, 1.0 * l0, 2.0);

    for (int i = 0; i < x.size(); i++) {
      logger::line(x[i], x[i] + 0.1 * _f[i], vec4(1.0, 0.0, 0.5, 1.0));
    }
  }

  virtual void debug_target() {
    for (int i = 0; i < _face_vert_ids.size(); i += 3) {
      vec3 x0 = _x[_face_vert_ids[i + 0]];
      vec3 x1 = _x[_face_vert_ids[i + 1]];
      vec3 x2 = _x[_face_vert_ids[i + 2]];
      logger::line(x0, x1, vec4(0.5, 0.0, 1.0, 1.0));
      logger::line(x1, x2, vec4(0.5, 0.0, 1.0, 1.0));
      logger::line(x2, x0, vec4(0.5, 0.0, 1.0, 1.0));
    }
  }

  virtual const std::vector<vec3> &forces() { return _f; }
  std::vector<index_t> _face_vert_ids;
  std::vector<vec3> _x;
  std::vector<vec3> _f;
  arp::T3::ptr _face_tree;
  real t = 0.0;
};

} // namespace duchamp
} // namespace gaudi

#endif