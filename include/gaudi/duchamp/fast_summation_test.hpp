#include "Eigen/src/Geometry/Scaling.h"
#include "gaudi/common.h"

#include <Eigen/Dense>

#include "gaudi/asawa/datums.hpp"
#include "gaudi/geometry_types.hpp"
#include "gaudi/geometry_logger.hpp"

#include "gaudi/calder/integrators.hpp"

#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/asawa/shell/shell_id.hpp"
#include "gaudi/asawa/shell/shell.hpp"

#include <array>

#include <math.h>
#include <random>

#include <cmath>
#include <memory>
#include <vector>

#ifndef __M2REFACTOR_TEST__
#define __M2REFACTOR_TEST__

namespace gaudi {
namespace duchamp {

using namespace asawa;

class fast_summation_test {
public:
  typedef std::shared_ptr<fast_summation_test> ptr;

  static ptr create() { return std::make_shared<fast_summation_test>(); }

  fast_summation_test() {
    //__M = load_cube();
    //__M = shell::load_messer();
    __M = shell::load_bunny();

    shell::triangulate(*__M);
    for (int i = 0; i < __M->face_count(); i++) {
      const auto fi = asawa::shell::face_id(i);
      if (__M->fbegin(fi) > asawa::shell::corner_id(0)) {
        assert(__M->fsize(fi) == 3);
      }
    }

    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    asawa::center(x, 2.0);
  };

  std::vector<vec3> createPoints(int N) {
    auto randNormalVec = [](real mean, real std) {
      auto randomFunc =
          [distribution_ = std::normal_distribution<double>(mean, std),
           random_engine_ = std::mt19937{std::random_device{}()}]() mutable {
            return vec3(distribution_(random_engine_),
                        distribution_(random_engine_),
                        distribution_(random_engine_));
            ;
          };
      return randomFunc;
    };

    std::vector<vec3> points;
    std::generate_n(std::back_inserter(points), N, randNormalVec(0, 0.5));
    return points;
  }

  std::vector<vec3> get_random_points(const int &N, const ext::extents_t &ext_t,
                                      shell::shell &M,
                                      const std::vector<vec3> &x) {
    ext::extents_t ext_m = asawa::shell::ext(M, x);
    std::uniform_real_distribution<real> dist(0.0, 1.0);
    std::mt19937_64 re;

    int i = 0;
    auto scaled_rand_vec = [dist, re, ext_t, ext_m, &i]() mutable {
      vec3 p(dist(re), dist(re), dist(re));
      if (i++ == 0) {
        p = vec3(0.57, 0.63, 0.0);
      }
      vec3 dt = ext_t[1] - ext_t[0];
      vec3 dm = ext_m[1] - ext_m[0];

      p = Eigen::Scaling(dt) * p + ext_t[0];
      p = Eigen::Scaling(dm) * p;
      p += ext_m[0];
      return p;
    };

    std::vector<vec3> tpoints;
    std::generate_n(std::back_inserter(tpoints), N, scaled_rand_vec);

    return tpoints;
  }

  void test_fast_winding(shell::shell &M, const std::vector<vec3> &x, real l0) {

    real zp = 0.5 + 0.5 * sin(M_PI * real(_frame) / 100.0);
    // real zp = 1.0;

    real zs = 0.01;
    std::vector<vec3> pov = get_random_points(
        10000, {vec3(0.0, 0.0, zp - zs), vec3(1.0, 1.0, zp + zs)}, M, x);
    std::vector<real> u = calder::fast_winding(M, x, pov, l0);

    for (int i = 0; i < pov.size(); i++) {
      geometry_logger::line(pov[i], pov[i] + 1e-4 * vec3(1.0, 0.0, 0.0),
                                geometry_logger::sdf4(1.0 * u[i]));
    }
  }

  void step(int frame) {

    _frame = frame;

    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    real l0 = 1.0 * asawa::shell::avg_length(*__M, x);

    test_fast_winding(*__M, x, l0);
  }
  int _frame;
  shell::shell::ptr __M;
};

} // namespace duchamp
} // namespace gaudi
#endif