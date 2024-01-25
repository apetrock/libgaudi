#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "gaudi/vec_addendum.h"

#include "GaudiGraphics/geometry_logger.h"

#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/dynamic.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/asawa/shell/shell.hpp"

#include "gaudi/calder/integrators.hpp"

#include "gaudi/asawa/primitive_objects.hpp"
#include "gaudi/common.h"
#include "gaudi/define_create_func.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <iterator>
#include <memory>
#include <set>
#include <vector>
#include <zlib.h>

#ifndef __SDF_FUNCTIONS__
#define __SDF_FUNCTIONS__

namespace gaudi {
namespace duchamp {
using namespace asawa;

// these will be moved to arp, I think
// power smooth min (k=8)
real smin_pow(real a, real b, real k) {
  a = pow(a, k);
  b = pow(b, k);
  return pow((a * b) / (a + b), 1.0 / k);
}

float smin_cubic(real a, real b, real k) {
  real h = std::max(k - abs(a - b), 0.0) / k;
  return std::min(a, b) - h * h * h * k * (1.0 / 6.0);
}

class sdf_base {
public:
  DEFINE_CREATE_FUNC(sdf_base)
  sdf_base() {}
  virtual ~sdf_base() {}
  virtual std::vector<real> distance(const std::vector<vec3> &x) const = 0;
  virtual std::vector<vec3> grad_distance(const std::vector<vec3> &x) const = 0;
  virtual void debug_target() const {}; // optional
};

class sdf_sphere : public sdf_base {
public:
  DEFINE_CREATE_FUNC(sdf_sphere)
  sdf_sphere(const vec3 &c, real r) : _c(c), _r(r) {}
  virtual ~sdf_sphere() {}

  virtual std::vector<real> distance(const std::vector<vec3> &x) const {
    std::vector<real> dists(x.size(), 0.0);
    for (int i = 0; i < x.size(); i++) {
      dists[i] = (x[i] - _c).norm() - _r;
    }
    return dists;
  }

  virtual std::vector<vec3> grad_distance(const std::vector<vec3> &x) const {
    std::vector<vec3> grads(x.size(), vec3(0.0, 0.0, 0.0));
    for (int i = 0; i < x.size(); i++) {
      grads[i] = (x[i] - _c).normalized();
    }
    return grads;
  };

  vec3 _c;
  real _r;
};

class sdf_multi_sphere : public sdf_base {
public:
  DEFINE_CREATE_FUNC(sdf_multi_sphere)
  sdf_multi_sphere(const std::vector<vec3> &c, real r) : _c(c), _r(r) {}
  virtual ~sdf_multi_sphere() {}

  virtual std::vector<real> distance(const std::vector<vec3> &x) const {

    std::vector<real> dists(x.size(), 999.0);
    for (int i = 0; i < x.size(); i++) {
      real d = 999.0;
      int jmin = -1;
      for (int j = 0; j < _c.size(); j++) {
        real dj = (x[i] - _c[j]).norm();
        if (dj < d) {
          d = dj;
          jmin = j;
        }
      }
      for (int j = 0; j < _c.size(); j++) {
        dists[i] = (x[i] - _c[jmin]).norm() - _r;
      }
    }
    return dists;
  }

  virtual std::vector<vec3> grad_distance(const std::vector<vec3> &x) const {
    std::vector<vec3> grads(x.size(), vec3(0.0, 0.0, 0.0));

    for (int i = 0; i < x.size(); i++) {
      // find index of closest point
      real d = 999.0;
      int jmin = -1;
      for (int j = 0; j < _c.size(); j++) {
        real dj = (x[i] - _c[j]).norm();
        if (dj < d) {
          d = dj;
          jmin = j;
        }
      }
      if (jmin > -1) {
        grads[i] = (x[i] - _c[jmin]).normalized();
      }
    }
    return grads;
  };

  std::vector<vec3> _c;
  real _r;
};

class sdf_cylinder : public sdf_base {
public:
  DEFINE_CREATE_FUNC(sdf_cylinder)
  sdf_cylinder(const vec3 &x0, const vec3 &x1, real r)
      : _x0(x0), _x1(x1), _r(r) {}
  virtual ~sdf_cylinder() {}

  virtual std::vector<real> distance(const std::vector<vec3> &x) const {
    std::vector<real> dists(x.size(), 0.0);
    for (int i = 0; i < x.size(); i++) {
      dists[i] = va::distance_from_line(_x0, _x1, x[i]) - _r;
    }
    return dists;
  }

  virtual std::vector<vec3> grad_distance(const std::vector<vec3> &x) const {};

  vec3 _x0, _x1;
  real _r;
};

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

class geometry_sdf : public sdf_base {
public:
  DEFINE_CREATE_FUNC(geometry_sdf)
  geometry_sdf(){};

  virtual ~geometry_sdf(){};

  virtual void add_geometry(const std::vector<vec3> &x,                //
                            const std::vector<index_t> &face_vert_ids, //
                            real s, vec3 offset) {

    _sets.push_back(triangle_set::create(x, face_vert_ids, s, offset));
  };

  virtual std::vector<real> distance(const std::vector<vec3> &x) const {

    std::vector<real> dists(x.size(), 999.0);
    for (int k = 0; k < _sets.size(); k++) {
      triangle_set &t = *_sets[k];

      std::vector<real> dists_i = calder::fast_dist(t.face_tree, x, 1.0);
      std::vector<real> winding = calder::fast_winding(t.face_tree, x, 1.0);
      for (int i = 0; i < x.size(); i++) {
        real di = va::sgn(0.5 - winding[i]) * dists_i[i];
        dists[i] = std::min(dists[i], di);
        // std::cout << "winding: " << winding[i] << " " << dists[i] <<
        // std::endl;
      }
    }
    return dists;
  }

  virtual std::vector<vec3> grad_distance(const std::vector<vec3> &x) const {
    std::vector<real> dists(x.size(), 100.0);
    std::vector<vec3> grads(x.size(), vec3(0.0, 0.0, 0.0));
    for (int k = 0; k < _sets.size(); k++) {

      triangle_set &t = *_sets[k];

      std::vector<real> dists_i = calder::fast_dist(t.face_tree, x, 1.0);
      std::vector<vec3> grads_i =
          calder::fast_dist_gradient(t.face_tree, x, 1.0);

      std::vector<real> winding = calder::fast_winding(t.face_tree, x, 1.0);

      for (int i = 0; i < x.size(); i++) {
        // logger::line(x[i], x[i] + 0.1 * grads_i[i], vec4(0.5,
        // 0.0, 1.0, 1.0));
        if (dists_i[i] < dists[i]) {
          dists[i] = dists_i[i];
          grads[i] = grads_i[i];
        }
      }
    }

    // debug_target();
    return grads;
  };

  virtual void debug_target() const {

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
  std::vector<triangle_set::ptr> _sets;
  arp::T3::ptr _face_tree;
  real t = 0.0;
};
} // namespace duchamp
} // namespace gaudi
#endif