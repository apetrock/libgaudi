#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "Eigen/src/Geometry/AngleAxis.h"
#include "gaudi/calder/shell_integrators.hpp"
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
#include "gaudi/define_create_func.h"

#include "gaudi/bontecou/laplacian.hpp"
#include "gaudi/duchamp/modules/vortex.hpp"
#include "gaudi/calder/least_squares_fit.hpp"
#include "modules/module_base.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <iterator>
#include <memory>
#include <set>
#include <vector>
#include <zlib.h>

#ifndef __TEST_TEMPLATE__
#define __TEST_TEMPLATE__

namespace gaudi
{
  namespace duchamp
  {
    using namespace asawa;

    class vortex_test
    {
    public:
      typedef std::shared_ptr<vortex_test> ptr;

      static ptr create() { return std::make_shared<vortex_test>(); }

      vortex_test() { load_shell(); };

      void load_shell()
      {
        //__M = load_cube();
        __M = shell::load_bunny();
        //__M = shell::load_crab();

        shell::triangulate(*__M);

        std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
        asawa::center(x, 2.0);

        real l0 = asawa::shell::avg_length(*__M, x);

        /////////
        // dynamic surface
        /////////
        // real C = 0.5;
        real C = 1.5;
        __surf = shell::dynamic::create(__M, C * l0, 2.5 * C * l0, 0.25 * C * l0);

        _vortex_edge = std::dynamic_pointer_cast<module_base>( //
            vortex_edge::create(__M, 1e-5));

        _vortex_vert = std::dynamic_pointer_cast<module_base>( //
            vortex_vert::create(__M, 1e-5, 5e-2, 1e-4));

        for (int i = 0; i < 5; i++)
        {
          __surf->step();
        }

        _eps = asawa::shell::avg_length(*__M, x);
      }

      std::vector<vec3> darboux_forces(asawa::shell::shell &shell, real C)
      {

        const std::vector<vec3> &x = asawa::const_get_vec_data(shell, 0);
        std::vector<vec3> N_v = asawa::shell::vertex_normals(shell, x);
        std::vector<real> sdf =
            calder::darboux_cyclide_sdf(shell, x, N_v, C * _eps, 3.0);

        std::vector<vec3> f(x.size(), vec3::Zero());
        auto range = shell.get_vert_range();

        for (auto i : range)
        {
          if (isnan(sdf[i]))
            continue;
          vec3 Ni = N_v[i];
          f[i] = -sdf[i] * N_v[i];
        }
        return f;
      }

      void step_dynamics(int frame)
      {

        asawa::shell::shell &M = *__M;
        std::vector<vec3> &x = asawa::get_vec_data(M, 0);
#if 1
        std::vector<vec3> fd = darboux_forces(M, 1.0);
#endif
#if 1
        std::dynamic_pointer_cast<vortex_edge>(_vortex_edge)->step(_h);
        std::vector<vec3> w0 =
            std::dynamic_pointer_cast<vortex_edge>(_vortex_edge)->get_vorticity();
        std::cout << "calder" << std::endl;
        std::vector<vec3> f0 = calder::vortex_force(M, x, w0, 6.0 * _eps, 3.0);
#endif
#if 0
    real mollifier1 = 4.0 * _eps;
    std::dynamic_pointer_cast<vortex_vert>(_vortex_vert)->step(_h);
    std::vector<vec3> w1 =
        std::dynamic_pointer_cast<vortex_vert>(_vortex_vert)->get_vorticity();
    std::vector<vec3> f1 = calder::vortex_force(M, x, w1, mollifier1, 3.0);
#endif
        // std::vector<vec3> w = calc_edge_vorticity(M, 1e-4 * _h);

#if 1
        // std::vector<vec3> fc = asawa::shell::vert_to_face<vec3>(M, f);
        //  std::vector<vec3> fs = calder::mls_avg<vec3>(M, fc, x, 4.0 * _eps, 3.0);

        for (int i = 0; i < x.size(); i++)
        {
          // x[i] += _h * (f0[i] + f1[i]);
           x[i] += _h * f0[i] + 1e-1 * _h *fd[i];
          //x[i] += 1e-1 * fd[i];
        }
#endif
        // asawa::center(x, 2.0);
      }

      void smoothMesh(real C, int N)
      {

        vec3_datum::ptr x_datum =
            static_pointer_cast<vec3_datum>(__M->get_datum(0));
        std::vector<vec3> &x = x_datum->data();
        bontecou::laplacian3 M(__M, x, true);
        M.init();
        real cc = C / 100.0;
        for (int k = 0; k < N; k++)
        {
          std::cout << "." << std::flush;
          x = M.smooth(x, C - cc, C + cc);
          int i = 0;
          for (auto xi : x)
          {
            if (!std::isfinite(xi.dot(xi)))
            {
              std::cout << xi.transpose() << std::endl;
              __M->vprintv(i);
              i++;
            }
          }
        }
        // x_datum->data() = x;
        std::cout << "done!" << std::endl;
      }

      void step(int frame)
      {
        std::cout << "frame: " << frame << std::endl;
        std::cout << "  -surface" << std::endl;
        std::cout << "    -corners: " << __M->corner_count() << std::endl;
        std::cout << "    -verts: " << __M->vert_count() << std::endl;
        std::cout << "    -faces: " << __M->face_count() << std::endl;

        for (int k = 0; k < 2; k++)
        {
          __surf->step(true);
        }
        step_dynamics(frame);
        smoothMesh(0.02, 10);
        // step_sdf(frame);
      }
      // std::map<index_t, index_t> _rod_adjacent_edges;
      real _h = 0.1;
      real _eps;
      index_t _ig0, _ig1;

      module_base::ptr _vortex_vert;
      module_base::ptr _vortex_edge;

      shell::shell::ptr __M;
      shell::dynamic::ptr __surf;
    };

  } // namespace duchamp
} // namespace gaudi
#endif