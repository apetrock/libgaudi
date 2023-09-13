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
    init_rx();
  }

  index_t _init_rx() {
    std::vector<real> rx(__M->vert_count(), 0.0);
    datum_t<real>::ptr adata = datum_t<real>::create(prim_type::VERTEX, rx);
    return __M->insert_datum(adata);
  }

  void init_rx() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0f, 1.0f);

    _irxa = _init_rx();
    _irxb = _init_rx();
    std::vector<real> &rxa = asawa::get_real_data(*__M, _irxa);
    std::vector<real> &rxb = asawa::get_real_data(*__M, _irxb);
    for (int i = 0; i < rxa.size(); i++) {
      real ta = dis(gen);
      real tb = dis(gen);
      // rxb[i] = 1.0;
      rxa[i] = 1.0;
      rxb[i] = 0.0;

      if (tb > 0.95) {
        rxb[i] = 1.0;
      }
    }
  }

  std::vector<vec4> get_mesh_colors() {
    std::vector<vec4> colors(__M->vert_count(), vec4(1.0, 0.0, 0.0, 1.0));
    std::vector<real> &rxa = asawa::get_real_data(*__M, _irxa);
    std::vector<real> &rxb = asawa::get_real_data(*__M, _irxb);
    vec4 col_a(0.75, 0.35, 0.0, 1.0);
    vec4 col_b(0.0, 2.0, 1.5, 1.0);

    for (int k = 0; k < __M->vert_count(); k++) {
      colors[k] = rxa[k] * col_a + 1.5 * rxb[k] * col_b;
      // if (k % 2 == 0)
      //   colors[k] = vec4(0.0, 0.0, 1.0, 1.0);
    }
    return colors;
  }

  void _diffuse(const std::vector<vec3> &x, std::vector<real> &f, real dt) {

    std::cout << x.size() << std::endl;
    bontecou::laplacian L(__M, x);
    std::vector<real> f_comp =
        asawa::shell::compress_to_vert_range<real>(*__M, f);

    std::vector<real> d = L.diffuse2(f_comp, dt);
    std::vector<real> d_exp =
        asawa::shell::expand_from_vert_range<real>(*__M, d);

    for (int k = 0; k < f.size(); k++) {
      f[k] = d_exp[k];
    }
  }

  void test_circulation() {
    std::vector<real> &rxa = asawa::get_real_data(*__M, _irxa);
    std::vector<real> &rxb = asawa::get_real_data(*__M, _irxb);
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

  void translate_normal(int frame) {
    std::vector<real> &rxa = asawa::get_real_data(*__M, _irxa);
    std::vector<real> &rxb = asawa::get_real_data(*__M, _irxb);
    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    std::vector<real> u = rxa;
    for (int i = 0; i < u.size(); i++) {
      u[i] = -0.90 * rxa[i] + 2.5 * rxb[i];
    }

    std::vector<vec3> Nv = asawa::shell::vertex_normals(*__M, x);
    for (int i = 0; i < x.size(); i++) {
      // real C = 2.0e-3 * pow(sin(M_PI * real(frame) / 400.0), 5.0);
      real C = 5.0e-4 * u[i];
      x[i] += C * Nv[i];
    }
  }

  void rx(real h) {
    std::vector<real> &rxa = asawa::get_real_data(*__M, _irxa);
    std::vector<real> &rxb = asawa::get_real_data(*__M, _irxb);
    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    // real f = 0.0438;
    // real k = 0.0598;
    real f = 0.025;
    real k = 0.0535;
    // real f = 0.0141;
    // real k = 0.0462;
    //   real f = 0.0485;
    //   real k = 0.0612;
    //    real f = 0.0475;
    //    real k = 0.0605;

    real da = 5.00e-4;
    real db = 0.4 * da;

    for (int i = 0; i < rxa.size(); i++) {
      real u0 = rxa[i];
      real v0 = rxb[i];
      vec2 un0 = vec2(u0, v0);
      auto [rxai, rxbi] = bontecou::grey_scott_2({rxa[i], rxb[i]}, f, k, h);
      rxa[i] = rxai, rxb[i] = rxbi;
    }

    _diffuse(x, rxa, h * da);
    _diffuse(x, rxb, h * db);
  }

  void step_dynamics(int frame) {
    rx(_h);
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
  index_t _il0;
  index_t _irxa, _irxb;
  shell::shell::ptr __M;
  shell::dynamic::ptr __surf;
};

} // namespace duchamp
} // namespace gaudi
#endif