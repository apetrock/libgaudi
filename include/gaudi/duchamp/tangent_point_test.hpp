#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "gaudi/vec_addendum.h"

#include "GaudiGraphics/geometry_logger.h"

#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/dynamic.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/asawa/shell/shell.hpp"

#include "gaudi/asawa/primitive_objects.hpp"
#include "gaudi/common.h"
#include "gaudi/define_create_func.h"

#include "gaudi/bontecou/laplacian.hpp"

#include "modules/module_base.hpp"
#include "modules/tangent_point.hpp"
#include "modules/tunnel.hpp"
#include "utils/sdf.hpp"
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

#ifndef __TEST_TEMPLATE__
#define __TEST_TEMPLATE__

namespace gaudi {
namespace duchamp {
using namespace asawa;

class tangent_point_test {
public:
  typedef std::shared_ptr<tangent_point_test> ptr;

  static ptr create() { return std::make_shared<tangent_point_test>(); }

  tangent_point_test() { load_shell(); };

  void load_shell() {
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
    __surf = shell::dynamic::create(__M, C * l0, 2.5 * C * l0, 0.5 * C * l0);

    _tangent_point =
        std::dynamic_pointer_cast<module_base>(tangent_point::create(__M));
    _tunnel = std::dynamic_pointer_cast<module_base>(tunnel::create(__M));

    for (int i = 0; i < 5; i++) {
      __surf->step();
    }

    _eps = asawa::shell::avg_length(*__M, x);
  }

  virtual void snap() {
    asawa::shell::shell &M = *__M;
    std::vector<vec3> &x = asawa::get_vec_data(M, 0);
    std::vector<index_t> edge_indices = M.get_edge_range();
    // random integer from 0 to size
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, edge_indices.size() - 1);
    index_t cA = -1;
    while (cA < 0) {
      index_t cAi = dis(gen);
      if (M.next(cAi) > 0) {
        cA = cAi;
      };
    }
    vec3 cen_A = asawa::shell::edge_center(M, cA, x);

    index_t cB = -1;
    real max = 0;
    for (int i = 0; i < edge_indices.size(); i++) {
      index_t cBi = edge_indices[i];
      vec3 cen_B = asawa::shell::edge_center(M, cBi, x);
      real d = (cen_B - cen_A).norm();
      if (d > max) {
        max = d;
        cB = cBi;
      }
    }

    index_t cA0 = cA;
    index_t cA1 = M.other(cA0);

    index_t cA0p = M.prev(cA0);
    index_t cA1p = M.prev(cA1);
    index_t cA0pp = M.prev(cA0p);
    index_t cA1pp = M.prev(cA1p);

    index_t cB0 = cB;
    index_t cB1 = M.other(cB0);

    index_t cB0p = M.prev(cB0);
    index_t cB1p = M.prev(cB1);
    index_t cB0pp = M.prev(cB0p);
    index_t cB1pp = M.prev(cB1p);

    std::cout << " ====== " << std::endl;
    std::cout << cA0 << " " << cA1 << " " << cA0p << " " << cA1p << " " << cA0pp
              << " " << cA1pp << std::endl;
    std::cout << cB0 << " " << cB1 << " " << cB0p << " " << cB1p << " " << cB0pp
              << " " << cB1pp << std::endl;
    vec3 xA0 = x[M.vert(cA0)];
    vec3 xA1 = x[M.vert(cA1)];
    vec3 xA0p = x[M.vert(cA0p)];
    vec3 xA1p = x[M.vert(cA1p)];

    vec3 xB0 = x[M.vert(cB0)];
    vec3 xB1 = x[M.vert(cB1)];
    vec3 xB0p = x[M.vert(cB0p)];
    vec3 xB1p = x[M.vert(cB1p)];

    x[M.vert(cA0)] = 0.5 * (xA0 + xB0);
    x[M.vert(cA1)] = 0.5 * (xA1 + xB1);
    x[M.vert(cA0p)] = 0.5 * (xA0p + xB0p);
    x[M.vert(cA1p)] = 0.5 * (xA1p + xB1p);
    x[M.vert(cB0)] = 0.5 * (xA0 + xB0);
    x[M.vert(cB1)] = 0.5 * (xA1 + xB1);
    x[M.vert(cB0p)] = 0.5 * (xA0p + xB0p);
    x[M.vert(cB1p)] = 0.5 * (xA1p + xB1p);
    asawa::shell::merge_edge(M, cA0, cB0);
    asawa::shell::merge_edge(M, cA0p, cB0p);
    asawa::shell::merge_edge(M, cA1p, cB1p);
    asawa::shell::merge_edge(M, cA0pp, cB0pp);
    asawa::shell::merge_edge(M, cA1pp, cB1pp);
  }

  void drill_sdf_tunnel(const std::vector<sdf_base::ptr> &sdfs) {
    asawa::shell::shell &M = *__M;

    for (int k = 0; k < 200; k++) {
      std::vector<vec3> &x = asawa::get_vec_data(M, 0);
      std::vector<vec3> xf = asawa::shell::face_centers(M, x);
      std::vector<vec3> Nf = asawa::shell::face_normals(M, x);

      std::vector<vec3> dxf(xf.size(), vec3::Zero());
      std::vector<real> df(xf.size(), 0.0);

      for (int j = 0; j < sdfs.size(); j++) {
        for (int i = 0; i < xf.size(); i++) {

          real d = sdfs[j]->distance(xf[i]);

          df[i] = std::min(df[i], d);
          // d = d < 0.0 ? d : 0.0;
        }
      }
      for (int i = 0; i < df.size(); i++) {
        dxf[i] = df[i] * Nf[i];
      }

      std::vector<vec3> dx = calder::mls_avg(M, dxf, x, 1.0 * _eps, 2.0);

      for (int i = 0; i < x.size(); i++) {
        x[i] += 2.0 * dx[i];
      }

      __surf->step();
    }
  }

  sdf_cylinder::ptr rand_cylinder() {
    asawa::shell::shell &M = *__M;
    std::vector<vec3> &x = asawa::get_vec_data(M, 0);
    std::vector<vec3> Nv = asawa::shell::vertex_normals(M, x);

    std::vector<index_t> vert_indices = M.get_vert_range();
    // random integer from 0 to size
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, vert_indices.size() - 1);
    index_t vA = -1;
    while (vA < 0) {
      index_t vAi = dis(gen);
      if (M.vbegin(vAi) > 0) {
        vA = vAi;
      };
    }

    index_t vB = -1;
    real max = 0;
    vec3 xA = x[vA];
    for (int i = 0; i < vert_indices.size(); i++) {
      index_t vBi = vert_indices[i];

      if (vBi == vA)
        continue;

      vec3 xB = x[vBi];
      vec3 dx = xB - xA;
      vec3 Na = Nv[vA];
      vec3 Nb = Nv[vBi];
      vec3 dxNa = Na * Na.transpose() * dx;
      vec3 dxNb = Nb * Nb.transpose() * dx;
      real d = 0.5 * (dxNa + dxNb).norm();

      if (d > max) {
        max = d;
        vB = vBi;
      }
    }
    vec3 xB = x[vB];
    logger::line(xA, xB, vec4(0.0, 1.0, 0.0, 1.0));

    real r = 4.0 * _eps;
    return sdf_cylinder::create(xA, xB, r);
  }

  void make_sdf_tunnel(int n = 4) {
    std::vector<sdf_base::ptr> sdfs;
    for (int i = 0; i < 1; i++)
      sdfs.push_back(rand_cylinder());
    drill_sdf_tunnel(sdfs);
  }

  void make_tunnel() {
    asawa::shell::shell &M = *__M;

    std::vector<vec3> &x = asawa::get_vec_data(M, 0);
    std::vector<vec3> Nv = asawa::shell::vertex_normals(M, x);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, x.size() - 1);
    index_t ii = dis(gen);

    std::cout << ii << " " << x.size() << " " << Nv.size() << std::endl;

    std::vector<vec3> start_points = {x[ii]};
    std::vector<vec3> headings = {-Nv[ii]};
    std::cout << "start_points: " << start_points[0].transpose() << std::endl;
    std::cout << "headings: " << headings[0].transpose() << std::endl;
    std::vector<std::vector<vec3>> paths =
        std::dynamic_pointer_cast<tunnel>(_tunnel)->drill(start_points,
                                                          headings);
  }

  void step_dynamics(int frame) {

    asawa::shell::shell &M = *__M;
    std::vector<vec3> &x = asawa::get_vec_data(M, 0);

#if 1
    real mollifier = 8.0 * _eps;
    std::dynamic_pointer_cast<tangent_point>(_tangent_point)->step(_h);
    // std::vector<vec3> X =
    //     std::dynamic_pointer_cast<tangent_point>(_tangent_point)
    //         ->get_smoothd_cross_grad();
#endif
    // asawa::center(x, 2.0);
  }

  void step(int frame) {
    std::cout << "frame: " << frame << std::endl;
    std::cout << "  -surface" << std::endl;
    std::cout << "    -corners: " << __M->corner_count() << std::endl;
    std::cout << "    -verts: " << __M->vert_count() << std::endl;
    std::cout << "    -faces: " << __M->face_count() << std::endl;
    // make_tunnel();
    if ((frame) % 200 == 1) {
      // std::cout << "  -snap" << frame << std::endl;
      make_sdf_tunnel();
      // snap();
    }

    for (int k = 0; k < 1; k++) {
      __surf->step(false, false);
    }

    for (int k = 0; k < 1; k++) {
      step_dynamics(frame);
    }

    if (frame > 1600)
      exit(0);
    // step_sdf(frame);
  }
  // std::map<index_t, index_t> _rod_adjacent_edges;
  real _h = 0.05;
  real _eps;
  index_t _ig0, _ig1;

  module_base::ptr _tangent_point;
  module_base::ptr _tunnel;

  shell::shell::ptr __M;
  shell::dynamic::ptr __surf;
};

} // namespace duchamp
} // namespace gaudi
#endif