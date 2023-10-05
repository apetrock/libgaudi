#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "Eigen/src/Geometry/AngleAxis.h"
#include "gaudi/vec_addendum.h"

#include "GaudiGraphics/geometry_logger.h"

#include "gaudi/asawa/rod/dynamic.hpp"
#include "gaudi/asawa/rod/rod.hpp"

#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/dynamic.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/asawa/shell/shell.hpp"

#include "gaudi/asawa/primitive_objects.hpp"
#include "gaudi/common.h"

#include "gaudi/bontecou/laplacian.hpp"

#include "modules/module_base.hpp"
#include "modules/reaction_diffusion.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <iterator>
#include <math.h>
#include <memory>
#include <set>
#include <vector>
#include <zlib.h>

#ifndef __VORTEX_STUDY__
#define __VORTEX_STUDY__

namespace gaudi {
namespace duchamp {
using namespace asawa;
class reaction_diffusion_test {
public:
  typedef std::shared_ptr<reaction_diffusion_test> ptr;

  static ptr create() { return std::make_shared<reaction_diffusion_test>(); }

  reaction_diffusion_test() { load_shell(); };

  void load_shell() {
    //__M = load_cube();
    __M = shell::load_bunny();
    //__M = shell::load_crab();

    shell::triangulate(*__M);

    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    asawa::center(x, 2.0);

    /////////
    // dynamic surface
    /////////
    real l0 = asawa::shell::avg_length(*__M, x);
    real C = 0.85;
    // real C = 2.0;
    __surf = shell::dynamic::create(__M, C * l0, 2.5 * C * l0, C * l0);
    real da = 5.00e-4, db = 0.4 * da;
    reaction_diffusion::ptr rx =
        reaction_diffusion::create(__M, 0.025, 0.0535, da, db);
    _rx = std::dynamic_pointer_cast<module_base>(rx);
    init_normals();
  }

  void init_normals() {
    _iN = gaudi::asawa::init_vert_datum<vec3>(*__M, vec3::Zero());
    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    std::vector<vec3> &N = asawa::get_vec_data(*__M, _iN);
    std::vector<vec3> Nv = asawa::shell::vertex_normals(*__M, x);

    for (int i = 0; i < N.size(); i++) {
      N[i] = Nv[i];
    }
  }

  std::vector<vec4> get_mesh_colors() {
    std::vector<vec4> colors(__M->vert_count(), vec4(1.0, 0.0, 0.0, 1.0));
    std::vector<real> &rxa =
        std::dynamic_pointer_cast<reaction_diffusion>(_rx)->get_rxa();
    std::vector<real> &rxb =
        std::dynamic_pointer_cast<reaction_diffusion>(_rx)->get_rxb();

    vec4 col_a(0.75, 0.35, 0.0, 1.0);
    vec4 col_b(0.0, 2.0, 1.5, 1.0);

    for (int k = 0; k < __M->vert_count(); k++) {
      colors[k] = rxa[k] * col_a + 1.5 * rxb[k] * col_b;
      // if (k % 2 == 0)
      //   colors[k] = vec4(0.0, 0.0, 1.0, 1.0);
    }
    return colors;
  }

  void test_circulation() {
    std::vector<real> &rxa =
        std::dynamic_pointer_cast<reaction_diffusion>(_rx)->get_rxa();
    std::vector<real> &rxb =
        std::dynamic_pointer_cast<reaction_diffusion>(_rx)->get_rxb();

    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    std::vector<real> u = rxa;
    for (int i = 0; i < u.size(); i++) {
      u[i] = 1.0 * rxa[i] - 1.0 * rxb[i];
    }

    std::vector<vec3> circulation = asawa::shell::circulation(*__M, u, x);

    for (int i = 0; i < circulation.size(); i++) {
      vec3 c = asawa::shell::face_center(*__M, i, x);
      gg::geometry_logger::line(c - 0.001 * circulation[i],
                                c + 0.002 * circulation[i],
                                vec4(1.0, 0.0, 0.75, 1.0));
    }
  }

  vec3 alerp(const vec3 &a, const vec3 &b, double t) {
    // Ensure the input vectors are normalized
    vec3 aN = a.normalized();
    vec3 bN = b.normalized();
    real adb = aN.dot(bN);
    real angle = std::acos(adb);
    // std::cout << "adb: " << adb << std::endl;
    // std::cout << "angle: " << angle << std::endl;
    // std::cout << aN.transpose() << " " << bN.transpose() << std::endl;
    if (adb < 1e-8) {
      return va::mix(t, a, b).normalized();
    }

    vec3 axis = aN.cross(bN).normalized();
    // std::cout << "axis: " << axis.transpose() << std::endl;

    real interpolatedAngle = angle * t;
    Eigen::AngleAxisd R(interpolatedAngle, axis);
    return R * a;
  }

  void translate_normal(int frame) {
    std::vector<real> &rxa =
        std::dynamic_pointer_cast<reaction_diffusion>(_rx)->get_rxa();
    std::vector<real> &rxb =
        std::dynamic_pointer_cast<reaction_diffusion>(_rx)->get_rxb();
    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    std::vector<real> u = rxa;

    std::vector<vec3> &N = asawa::get_vec_data(*__M, _iN);
    std::vector<vec3> Nv = asawa::shell::vertex_normals(*__M, x);
    index_t i = 0;

    for (auto &n : N) {
      // n = alerp(n, Nv[i++], 0.05);
      n = va::mix(0.01, n, Nv[i++]);
      n.normalize();
    }

    for (int i = 0; i < u.size(); i++) {
      u[i] = -0.90 * rxa[i] + 2.5 * rxb[i];
    }

    for (int i = 0; i < x.size(); i++) {
      // real C = 2.0e-3 * pow(sin(M_PI * real(frame) / 400.0), 5.0);
      real C = 5.0e-4 * u[i];
      x[i] += C * N[i];
    }
  }

  void step_dynamics(int frame) {
    _rx->step(_h);
    translate_normal(frame);
    //  test_circulation();
  }

  void step(int frame) {
    std::cout << "frame: " << frame << std::endl;
    std::cout << "  -surface" << std::endl;
    std::cout << "    -corners: " << __M->corner_count() << std::endl;
    std::cout << "    -verts: " << __M->vert_count() << std::endl;
    std::cout << "    -faces: " << __M->face_count() << std::endl;
    for (int k = 0; k < 5; k++) {
      step_dynamics(frame);
    }

    for (int k = 0; k < 1; k++) {
      __surf->step(true);
    }
    if (frame > 2400)
      exit(0);
    // step_sdf(frame);
  }
  // std::map<index_t, index_t> _rod_adjacent_edges;
  real _h = 8.0e-1;

  index_t _iN;

  module_base::ptr _rx;

  shell::shell::ptr __M;
  shell::dynamic::ptr __surf;
};

} // namespace duchamp
} // namespace gaudi
#endif