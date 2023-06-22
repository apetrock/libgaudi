#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "gaudi/vec_addendum.h"

#include "GaudiGraphics/geometry_logger.h"

#include "gaudi/asawa/rod/dynamic.hpp"
#include "gaudi/asawa/rod/rod.hpp"

#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/dynamic.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/asawa/shell/shell.hpp"

#include "gaudi/asawa/shell/walk.hpp"

#include "gaudi/bontecou/laplacian.hpp"
#include "gaudi/hepworth/block/generic_constraints.hpp"
#include "gaudi/hepworth/block/generic_constraints_init.hpp"

#include "gaudi/hepworth/block/rod_constraints.hpp"
#include "gaudi/hepworth/block/rod_constraints_init.hpp"

#include "gaudi/hepworth/block/shell_constraints.hpp"
#include "gaudi/hepworth/block/shell_constraints_init.hpp"

#include "gaudi/hepworth/block/coupling_collisions_init.hpp"
#include "gaudi/hepworth/block/sim_block.hpp"
#include "gaudi/hepworth/block/solver.hpp"

#include "gaudi/calder/integrators.hpp"

#include "gaudi/asawa/primitive_objects.hpp"
#include "gaudi/common.h"

#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>
#include <zlib.h>

#ifndef __M2REFACTOR_TEST__
#define __M2REFACTOR_TEST__

namespace gaudi {
namespace duchamp {

using namespace asawa;

class shell_normal_constraints {
public:
  typedef std::shared_ptr<shell_normal_constraints> ptr;

  static ptr create() { return std::make_shared<shell_normal_constraints>(); }

  shell_normal_constraints() {
    //__M = load_cube();
    __M = shell::load_bunny();
    //__M = shell::load_crab();

    shell::triangulate(*__M);
    for (int i = 0; i < __M->face_count(); i++) {
      if (__M->fbegin(i) > 0) {
        assert(__M->fsize(i) == 3);
      }
    }

    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    asawa::center(x);

    /////////
    // dynamic surface
    /////////

    real l0 = 0.75 * asawa::shell::avg_length(*__M, x);
    __surf = shell::dynamic::create(__M, 0.5 * l0, 3.0 * l0, 1.0 * l0);

    /////////
    // constraints setup
    /////////
    std::vector<real> lr = asawa::shell::edge_lengths(*__M, x);

    datum_t<real>::ptr ldata = datum_t<real>::create(prim_type::EDGE, lr);
    _il0 = __M->insert_datum(ldata);
    walk_and_subdivide();
    // load_loop_rod();

    __surf->step(true);
  };

  void init_rod_shell_weld(
      asawa::rod::rod &R, asawa::rod::dynamic &rod_d, asawa::shell::shell &M,
      asawa::shell::dynamic &shell_d,
      std::vector<hepworth::projection_constraint::ptr> &constraints,
      const real &w, std::vector<hepworth::sim_block::ptr> blocks) {

    const std::vector<vec3> &x0 = R.__x;
    std::vector<vec3> &x1 = asawa::get_vec_data(M, 0);
    std::vector<index_t> edge_verts_0 = R.get_edge_vert_ids();
    std::vector<index_t> edge_map_0 = R.get_edge_map();
    real eps = 1e-4;
    vector<std::array<index_t, 2>> rs_collisions =
        shell_d.get_edge_edge_collisions(edge_verts_0, edge_map_0, x0, eps);
    for (auto &c : rs_collisions) {
      if (c[0] > -1) {
        index_t cr0 = c[0];
        index_t cr1 = R.next(cr0);
        index_t cs0 = c[1];
        index_t cs1 = M.other(cs0);

        if (cr1 < 0 || cr1 < 0 || cs0 < 0 || cs1 < 0)
          continue;
        index_t vs0 = M.vert(cs0);
        index_t vs1 = M.vert(cs1);

        if (vs0 < 0 || vs1 < 0)
          continue;

        vec3 xA0 = x0[cr0];
        vec3 xA1 = x0[cr1];
        vec3 xB0 = x1[vs0];
        vec3 xB1 = x1[vs1];

        constraints.push_back(hepworth::block::edge_edge_weld::create(
            std::vector<index_t>({cr0, cr1, vs0, vs1}), w, blocks));
      }
    }

    std::vector<index_t> edge_verts_1 = M.get_edge_vert_ids();
    std::vector<index_t> edge_map_1 = M.get_edge_map();

    vector<std::array<index_t, 4>> sr_collisions =
        rod_d.get_collisions(edge_verts_1, x1);
    for (auto &c : sr_collisions) {
      if (c[0] > -1) {
        index_t vs0 = c[0];
        index_t vs1 = c[1];
        index_t vr0 = c[2];
        index_t vr1 = c[3];

        if (vr1 < 0 || vr1 < 0 || vs0 < 0 || vs1 < 0)
          continue;
        vec3 xA0 = x0[vr0];
        vec3 xA1 = x0[vr1];
        vec3 xB0 = x1[vs0];
        vec3 xB1 = x1[vs1];
        std::array<real, 3> d =
            va::distance_Segment_Segment(xA0, xA1, xB0, xB1);
        if (d[0] > eps)
          continue;

        constraints.push_back(hepworth::block::edge_edge_weld::create(
            std::vector<index_t>({vr0, vr1, vs0, vs1}), w, blocks));
      }
    }
  }

  void load_loop_rod() {
    std::uniform_real_distribution<real> dist(0.5, 1.0);
    std::mt19937_64 re;
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
    int N = 256;
    for (int i = 0; i < N; i++) {
      real thet = 2.0 * M_PI * real(i) / real(N);
      std::cout << cos(thet) << " " << sin(thet) << std::endl;
      vec3 pi = r0 * cos(thet) * f0 + r1 * sin(thet) * f1;
      points.push_back(pi);
    }
    __R = rod::rod::create(points);

    real lavg = __R->lavg();
    __Rd = rod::dynamic::create(__R, 0.25 * lavg, 2.5 * lavg, 0.25 * lavg);
  }

  void walk_and_subdivide() {
    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    std::vector<vec3> &v = asawa::get_vec_data(*__M, 1);
    index_t test = 1000;

    vec3 N = asawa::shell::edge_normal(*__M, test, x);
    vec3 T = asawa::shell::edge_tangent(*__M, test, x).normalized();
    vec3 B = N.cross(T).normalized();

    real thet = 1e-2 * 2.0 * M_PI * 0.23;

    vec3 dir = std::cos(thet) * T + std::sin(thet) * B;
    std::vector<index_t> corners;
    std::vector<real> S;
    std::vector<vec3> points;
    std::vector<vec3> normals;

    asawa::shell::walk(*__M, x, test, dir, 0.5, 1000, 1e-2,
                       [&](const asawa::shell::shell &M,
                           const std::vector<vec3> &x, const index_t &corner,
                           const real &s, const real &accumulated_length,
                           vec3 &dir) {
                         asawa::shell::index_t ci = corner;
                         S.push_back(s);
                         corners.push_back(corner);
                         points.push_back(asawa::shell::edge_vert(M, ci, s, x));
                         vec3 Ni = asawa::shell::edge_normal(*__M, corner, x);
                         normals.push_back(Ni);
                         return true;
                       });

    __R = rod::rod::create(points, false);
    __R->_update_frames(normals);
    real lavg = __R->lavg();
    __R->_r = 0.015;
    __Rd = rod::dynamic::create(__R, 0.25 * lavg, 2.5 * lavg, 0.25 * lavg);
  }
  std::vector<vec3> compute_tangent_point_gradient() {
    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    std::vector<vec3> N = asawa::shell::vertex_normals(*__M, x);
    std::vector<vec3> g =
        calder::tangent_point_grad(*__R, x, N, 0.5 * __surf->_Cm);

    std::vector<vec3> Nr = calder::mls_edge_interp(*__R, x, __surf->_Cm);
    std::vector<vec3> dp = calder::mls_dp(*__R, x, __surf->_Cm);

    for (int i = 0; i < Nr.size(); i++) {
      // std::cout << g[i] << std::endl;
      // g[i] = compute_grad(dp[i], Nr[i], 2.0);
      g[i] *= 1e-6;
      gg::geometry_logger::line(x[i], x[i] + 0.01 * g[i],
                                vec4(0.6, 0.0, 0.8, 1.0));
    }
    return g;
  }

  void init_normal_constraints(
      std::vector<hepworth::projection_constraint::ptr> &constraints,
      std::vector<hepworth::sim_block::ptr> blocks) {
    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    std::vector<vec3> x_f = asawa::shell::face_centers(*__M, x);

    std::vector<vec3> Nr = calder::mls_edge_interp(*__R, x_f, __surf->_Cm);

    std::vector<real> d = calder::mls_dist(*__R, x_f, __surf->_Cm);

    real w = 0.0;
    for (int i = 0; i < d.size(); i++) {
      vec3 xi = x_f[i];
      w = std::max(w, 1.0 / pow(d[i], 2.0));
    }

#if 1
    auto range = __M->get_face_range();
    int i = 0;
    for (auto fi : range) {
      if (__M->fsize(fi) != 3)
        continue;
      if (asawa::shell::face_area(*__M, fi, x) < 1e-8)
        continue;

      vec3 xi = x_f[i];
      vec3 Ni = Nr[i].normalized();
      real di = 1.0 / w / pow(d[i], 2.0);

      gg::geometry_logger::line(xi, xi + 0.1 * di * Ni,
                                vec4(0.1, 0.7, 0.2, 1.0));

      auto tri = __M->get_tri(fi);

      constraints.push_back(
          hepworth::block::normal::create(std::vector<index_t>({
                                              tri[0],
                                              tri[1],
                                              tri[2],
                                          }),
                                          x, Ni, 0.1 * di, blocks));
      i++;
    }
#endif
  }

  void step_dynamics(int frame) {

    hepworth::block::projection_solver solver;

    std::vector<hepworth::projection_constraint::ptr> constraints;

    std::vector<real> &l0 = asawa::get_real_data(*__M, _il0);
    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    std::vector<vec3> &v = asawa::get_vec_data(*__M, 1);

    std::vector<vec3> M = asawa::shell::vertex_areas(*__M, x);
    std::vector<real> li = asawa::shell::edge_lengths(*__M, x);

    std::vector<real> &lr = __R->__l0;
    // for (auto &l : lr)
    //   l *= 1.01;

    int i = 0;
    auto range = __M->get_edge_range();
    for (auto c0 : range) {
      if (std::isnan(l0[c0 / 2]))
        continue;
      l0[c0 / 2] = asawa::shell::edge_length(*__M, c0, x);
    }

    std::vector<vec3> Ns = asawa::shell::vertex_normals(*__M, x);
    // std::vector<vec3> f(x.size(), vec3::Zero());
    std::vector<vec3> f = compute_tangent_point_gradient();

    hepworth::vec3_block::ptr Xs = hepworth::vec3_block::create(M, x, v, f);

    hepworth::vec3_block::ptr Xr =
        hepworth::vec3_block::create(__R->__M, __R->__x, __R->__v, f);
    hepworth::quat_block::ptr Ur =
        hepworth::quat_block::create(__R->__J, __R->__u, __R->__o);

    std::vector<hepworth::sim_block::ptr> blocks;
    // blocks.push_back(Xs);
    blocks.push_back(Xr);
    blocks.push_back(Ur);

    std::cout << "l[0]: " << l0[0] << std::endl;

    // hepworth::block::init_edge_growth(*__M, constraints, x, l0, 0.99, 1.0e-2,
    //                                  {Xs});

    // hepworth::block::init_edge_strain(*__M, constraints, x, l0, 1e-2, {Xs});
    // hepworth::block::init_bending(*__M, constraints, x, 1.0e-2, {Xs});
    // hepworth::block::init_laplacian(*__M, constraints, x, 0, 1.0e-0, {Xs});
    // hepworth::block::init_area(*__M, constraints, x, 1.50e-1, {Xs});
    //   init_normal_constraints(constraints, {Xs});

    hepworth::block::init_pinned(*__R, constraints, __R->__x, 0.1, {Xr});
    hepworth::block::init_stretch_shear(*__R, constraints, lr, 0.1, {Xr, Ur});
    hepworth::block::init_bend_twist(*__R, constraints, 0.1, {Ur});

    real eps = 3.0 * __surf->_Cm;
    /*
    hepworth::block::init_pnt_tri_collisions(
        *__M, *__surf, constraints, x, 1.0 * eps, 0.5 * eps, 1.0, {Xs, Xs});
    init_rod_shell_weld(*__R, *__Rd, //
                        *__M, *__surf, constraints, 100.0, {Xr, Xs});
    */
    solver.set_constraints(constraints);

    solver.step(blocks, _h, 20);
  }

  void step(int frame) {
    std::cout << "frame: " << frame << std::endl;
    step_dynamics(frame);
    __surf->step(true);
    __Rd->step();
  }

  real _h = 0.1;
  index_t _il0;
  index_t _iwalk;
  std::vector<index_t> _walk;
  shell::shell::ptr __M;
  shell::dynamic::ptr __surf;
  rod::rod::ptr __R;
  rod::dynamic::ptr __Rd;
};

} // namespace duchamp
} // namespace gaudi
#endif