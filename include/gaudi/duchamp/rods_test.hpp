#include "Eigen/src/Geometry/Scaling.h"
#include "gaudi/common.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "gaudi/geometry_types.hpp"
#include "gaudi/vec_addendum.h"

#include "GaudiGraphics/geometry_logger.h"

#include "gaudi/arp/arp.h"

#include "gaudi/calder/integrators.hpp"

#include "gaudi/bontecou/laplacian.hpp"

#include "gaudi/asawa/dynamic_surface.hpp"
#include "gaudi/asawa/shell.hpp"

#include "gaudi/asawa/rod.hpp"

#include "gaudi/asawa/datum_x.hpp"

#include "gaudi/asawa/primitive_objects.hpp"
#include "gaudi/asawa/shell_operations.hpp"

#include "gaudi/asawa/asset_loader.hpp"

#include "gaudi/calder/tree_code.hpp"

#include <array>

#include <math.h>
#include <random>

#include <cmath>
#include <memory>
#include <vector>
#include <zlib.h>

#ifndef __M2REFACTOR_TEST__
#define __M2REFACTOR_TEST__

namespace gaudi {
namespace duchamp {

using namespace asawa;

void center(std::vector<vec3> &coords) {
  real accum = 0.0;
  vec3 min = coords[0];
  vec3 max = coords[0];

  for (auto &c : coords) {
    min = va::min(c, min);
    max = va::max(c, max);
  }

  vec3 dl = (max - min);
  real maxl = dl[0];
  maxl = maxl > dl[1] ? maxl : dl[1];
  maxl = maxl > dl[2] ? maxl : dl[2];
  real s = 2.0 / maxl;
  vec3 cen = (0.5 * (max + min) - min) * s;
  std::cout << " scale: " << s << std::endl;
  cen -= min;
  for (auto &c : coords) {
    c -= min;
    c = s * c;
    c -= cen;
  }
}

class rods_test {
public:
  typedef std::shared_ptr<rods_test> ptr;

  static ptr create() { return std::make_shared<rods_test>(); }

  rods_test() {
    //__M = load_cube();

    __M = load_bunny();
    load_loop_rod();

    triangulate(*__M);
    for (int i = 0; i < __M->face_count(); i++) {
      if (__M->fbegin(i) > 0) {
        assert(__M->fsize(i) == 3);
      }
    }

    vec3_datum::ptr x_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(0));

    std::vector<vec3> &x = x_datum->data();
    center(x);
  };

  void load_loop_rod() {
    std::uniform_real_distribution<real> dist(0.5, 1.0);
    std::mt19937_64 re;
    double a_random_double = unif(re);
    vec3 p0(dist(re), dist(re), dist(re));
    vec3 p1(dist(re), dist(re), dist(re));
    real r0 = p0.norm();
    real r1 = p1.norm();

    p0.normalize();
    p1.normalize();
    vec3 f2 = p0.cross(p1).normalized();
    vec3 f1 = p0.cross(f2).normalized();
    vec3 f0 = f1.cross(f2).normalized();

    std::cout << p0 << ::std::endl;
    std::cout << ::std::endl;
    std::cout << p1 << ::std::endl;

    std::vector<vec3> points;
    int N = 1024;
    for (int i = 0; i < N; i++) {
      real thet = 2.0 * M_PI * real(i) / real(N);
      std::cout << cos(thet) << " " << sin(thet) << std::endl;
      vec3 pi = r0 * cos(thet) * f0 + r1 * sin(thet) * f1;
      points.push_back(pi);
    }
    __R = asawa::rod::create(points);
  }

  void step(int frame) {

    _frame = frame;

    vec3_datum::ptr x_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(0));
    std::vector<vec3> &x = x_datum->data();
    __R->step();
    __R->debug();
  }
  int _frame;
  shell::ptr __M;
  rod::ptr __R;
};

} // namespace duchamp
} // namespace gaudi
#endif