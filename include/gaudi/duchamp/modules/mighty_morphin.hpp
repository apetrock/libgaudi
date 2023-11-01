
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
#include <limits>
#include <vector>

namespace gaudi {

// probably doesn't need to be a module, could go straight in test class...
/////////////////////////////////////////////////////////////

namespace duchamp {

struct triangle_set {
  DEFINE_CREATE_FUNC(triangle_set)
  triangle_set(const std::vector<vec3> &x_in,
               const std::vector<index_t> &fid_in, real scale, vec3 offset)
      : x(x_in), face_vert_ids(fid_in) {

    for (int i = 0; i < x.size(); i++) {
      x[i] = scale * x[i] - offset;
    }

    face_tree = arp::T3::create(face_vert_ids, x, 24);
  }
  std::vector<vec3> x;
  std::vector<index_t> face_vert_ids;
  arp::T3::ptr face_tree;
};

class mighty_morphin : public module_base_shell {
public:
  DEFINE_CREATE_FUNC(mighty_morphin)
  mighty_morphin(const asawa::shell::shell::ptr &M) : module_base_shell(M){};

  virtual ~mighty_morphin(){};

  virtual void add_geometry(const std::vector<vec3> &x,                //
                            const std::vector<index_t> &face_vert_ids, //
                            real s, vec3 offset) {

    _sets.push_back(triangle_set::create(x, face_vert_ids, s, offset));
  };

  virtual std::vector<real> calc_dist(const std::vector<vec3> &x,
                                      triangle_set &tri_set) {

    std::vector<index_t> vert_ids(x.size(), -1);
    for (int i = 0; i < x.size(); i++) {
      vert_ids[i] = i;
    }
    std::vector<real> dists(x.size(), 0.0);
    for (int i = 0; i < x.size(); i++) {

      std::vector<index_t> nearest = arp::getNearest<1, 3>(
          i, vert_ids, x,     //
          *tri_set.face_tree, //
          std::numeric_limits<real>::infinity(), &arp::pnt_tri_min);

      index_t j = nearest.back();
      if (j < 0)
        continue;
      vec3 xi = x[i];

      vec3 x0 = tri_set.x[tri_set.face_vert_ids[3 * j + 0]];
      vec3 x1 = tri_set.x[tri_set.face_vert_ids[3 * j + 1]];
      vec3 x2 = tri_set.x[tri_set.face_vert_ids[3 * j + 2]];

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

    // std::vector<vec3> Nv = asawa::shell::vertex_normals(*_M, x);
    real l0 = asawa::shell::avg_length(*_M, x);
    std::vector<vec3> Nv = asawa::shell::vertex_normals(*_M, x);

    std::vector<real> dists(x.size(), 100.0);
    for (int k = 0; k < _sets.size(); k++) {
      triangle_set &t = *_sets[k];

      std::vector<real> dists_i = calder::fast_dist(t.face_tree, x, 1.0);
      std::vector<real> winding = calder::fast_winding(t.face_tree, x, 1.0);
      std::vector<real> max_dists_i =
          calder::fast_view(t.face_tree, x, Nv, 1.0);
      for (int i = 0; i < x.size(); i++) {
        real di = va::sgn(0.5 - winding[i]) * dists_i[i];
        if (di < 0.0)
          di *= pow(max_dists_i[i], 1.0);
        dists[i] = std::min(dists[i], di);
      }
    }

    _f = std::vector<vec3>(x.size(), vec3::Zero());
    auto maxElement =
        std::max_element(dists.begin(), dists.end(),
                         [](const real &v1, real &v2) { return v1 < v2; });
    real ih_max = *maxElement;
    ih_max = std::max(ih_max, 0.1 * l0);
    for (int i = 0; i < x.size(); i++) {
      vec3 xi = x[i];
      vec3 fi = -dists[i] * Nv[i].normalized();
      // logger::line(xi, xi + 0.1 * fi, vec4(0.0, 1.0, 0.5, 1.0));
      _f[i] += fi / ih_max;
    }
    /*
    for (int i = 0; i < x.size(); i++) {
      logger::line(x[i], x[i] + 0.1 * _f[i], vec4(1.0, 0.0, 0.5, 1.0));
    }
*/
  }

  virtual void debug_target() {

    for (int k = 0; k < _sets.size(); k++) {
      triangle_set &t = *_sets[k];
      for (int i = 0; i < t.face_vert_ids.size(); i += 3) {
        vec3 x0 = t.x[t.face_vert_ids[i + 0]];
        vec3 x1 = t.x[t.face_vert_ids[i + 1]];
        vec3 x2 = t.x[t.face_vert_ids[i + 2]];
        logger::line(x0, x1, vec4(0.5, 0.0, 1.0, 1.0));
        logger::line(x1, x2, vec4(0.5, 0.0, 1.0, 1.0));
        logger::line(x2, x0, vec4(0.5, 0.0, 1.0, 1.0));
      }
    }
  }

  virtual const std::vector<vec3> &forces() { return _f; }

  std::vector<vec3> _f;
  std::vector<triangle_set::ptr> _sets;
  arp::T3::ptr _face_tree;
  real t = 0.0;
};

} // namespace duchamp
} // namespace gaudi

#endif