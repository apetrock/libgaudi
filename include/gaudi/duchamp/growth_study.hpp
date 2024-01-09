#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "gaudi/vec_addendum.h"

#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/dynamic.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/asawa/shell/shell.hpp"

#include "gaudi/asawa/shell/walk.hpp"

#include "gaudi/bontecou/laplacian.hpp"

#include "gaudi/hepworth/block/generic_constraints.hpp"
#include "gaudi/hepworth/block/generic_constraints_init.hpp"

#include "gaudi/hepworth/block/shell_constraints.hpp"
#include "gaudi/hepworth/block/shell_constraints_init.hpp"
#include "gaudi/hepworth/block/sim_block.hpp"
#include "gaudi/hepworth/block/solver.hpp"

#include "gaudi/asawa/primitive_objects.hpp"
#include "gaudi/calder/quadric_fit.hpp"
#include "gaudi/calder/tangent_point_integrators.hpp"
#include "gaudi/common.h"
#include "gaudi/logger.hpp"

#include "modules/reaction_diffusion.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <math.h>
#include <memory>
#include <vector>
#include <zlib.h>

#ifndef __M2REFACTOR_TEST__
#define __M2REFACTOR_TEST__

namespace gaudi {
namespace duchamp {

using namespace asawa;

class growth_study {
public:
  typedef std::shared_ptr<growth_study> ptr;

  static ptr create() { return std::make_shared<growth_study>(); }

  growth_study() {
    //__M = load_cube();
    __M = shell::load_bunny();

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

    real l0 = 1.5 * asawa::shell::avg_length(*__M, x);
    _eps = l0;
    __surf = shell::dynamic::create(__M, 1.0 * l0, 2.5 * l0, 1.0 * l0);

    __surf->step(true);
    __surf->step(true);
    __surf->step(true);
    // real f = 0.075, k = 0.0615;
    real f0 = 0.031, k0 = 0.0585;
    // real f = 0.04, k = 0.065;
    //  real f = 0.049, k = 0.0629;
    //  real f = 0.025, k = 0.0535;

    real da0 = 3.0e-4, db0 = 0.5 * da0;
    reaction_diffusion::ptr rx0 =
        reaction_diffusion::create(__M, f0, k0, da0, db0);
    _rx0 = std::dynamic_pointer_cast<module_base>(rx0);
    /*
    real f1 = 0.04, k1 = 0.065;
    real da1 = 2.00e-4, db1 = 0.5 * da1;
    reaction_diffusion::ptr rx1 =
        reaction_diffusion::create(__M, f1, k1, da1, db1);
    _rx1 = std::dynamic_pointer_cast<module_base>(rx1);
*/
    init_origin();
  };

  void init_origin() {
    const std::vector<vec3> &x = asawa::const_get_vec_data(*__M, 0);

    std::array<vec3, 2> ext = asawa::extents(x);
    _origin = vec3(                    //
        0.5 * (ext[0][0] + ext[1][0]), //
        ext[0][1],                     //
        0.5 * (ext[0][2] + ext[1][2]));
  }

  std::vector<vec4> get_mesh_colors() {
    std::vector<vec4> colors(__M->vert_count(), vec4(1.0, 0.0, 0.0, 1.0));
    std::vector<real> &rx0a =
        std::dynamic_pointer_cast<reaction_diffusion>(_rx0)->get_rxa();
    std::vector<real> &rx0b =
        std::dynamic_pointer_cast<reaction_diffusion>(_rx0)->get_rxb();
    /*
        std::vector<real> &rx1a =
            std::dynamic_pointer_cast<reaction_diffusion>(_rx1)->get_rxa();
        std::vector<real> &rx1b =
            std::dynamic_pointer_cast<reaction_diffusion>(_rx1)->get_rxb();

    */
    vec4 col_a(1.0, 0.0, 1.0, 1.0);
    vec4 col_b(0.0, 1.0, 1.0, 1.0);
    vec4 col_c(1.0, 0.0, 1.0, 1.0);
    for (int k = 0; k < __M->vert_count(); k++) {
      vec4 col0 = 1.0 * rx0a[k] * col_a + //
                  3.0 * rx0b[k] * col_b;
      colors[k] = col0;
      // if (k % 2 == 0)
      //   colors[k] = vec4(0.0, 0.0, 1.0, 1.0);
    }
    return colors;
  }
  vec3 get_origin() { return _origin; }

  index_t get_new_origin(asawa::shell::shell &shell) {
    const std::vector<vec3> &x = asawa::const_get_vec_data(shell, 0);
    index_t imin = 0.0;
    real min = std::numeric_limits<real>::max();
    for (int i = 0; i < x.size(); i++) {
      real d = (x[i] - _origin).norm();
      if (d < min) {
        min = d;
        imin = i;
      }
    }
    _origin = x[imin];
    return imin;
  }

  std::vector<real> vertex_geodesic_weight(asawa::shell::shell &shell) {
    // weights on edges
    const std::vector<vec3> &x = asawa::const_get_vec_data(shell, 0);

    index_t imin = get_new_origin(shell);

    std::vector<real> f(shell.vert_count(), 0.0);
    f[imin] = 1.0;

    bontecou::laplacian L(__M, x);
    std::vector<real> d = L.heatDist(f, 0.2);
    return d;
  }

  std::vector<real> medial_axis(asawa::shell::shell &shell, real w = 1.0) {

    const std::vector<vec3> &x = asawa::const_get_vec_data(shell, 0);
    std::vector<vec3> N_v = asawa::shell::vertex_normals(shell, x);

    std::vector<vec3> x_e = asawa::shell::edge_centers(shell, x);
    std::vector<real> w_e = asawa::shell::edge_areas(shell, x);
    std::vector<vec3> N_f = asawa::shell::face_normals(shell, x);
    std::vector<vec3> N_e =
        calder::mls_avg<vec3>(shell, N_f, x_e, 2.0 * _eps, 2.0);
    std::vector<gaudi::calder::vec7> M =
        calder::cylinder(shell, x_e, N_e, 4.0 * _eps, 3.0);

    std::vector<real> g_edge(__M->edge_count(), 0.0);
    auto range = shell.get_edge_range();

    for (auto c0 : range) {
      int i = shell.vert(c0);
      int j = shell.vert(shell.other(c0));

      vec3 xi = x[i];
      vec3 xj = x[j];
      vec3 xc = 0.5 * (xi + xj);
      vec3 e = (xj - xi).normalized();
      vec3 cen = M[c0 / 2].segment(0, 3);
      vec3 f0 = M[c0 / 2].segment(3, 3);
      vec3 dx = xj - xi;
      real e0 = pow(e.dot(f0), 2.0);
      g_edge[c0 / 2] = e0 * w;
      // if (c0 / 2 == 1000)
      //   logger::line(xc, cen, vec4(1.0, 0.0, 0.0, 1.0));
    }
    return g_edge;
  }

  std::vector<real> covariant_growth_weight(asawa::shell::shell &shell,
                                            vec2 w = vec2(1.1, 0.1)) {

    const std::vector<vec3> &x = asawa::const_get_vec_data(shell, 0);

    std::vector<vec3> x_e = asawa::shell::edge_centers(shell, x);
    std::vector<vec3> N_f = asawa::shell::face_normals(shell, x);
    std::vector<vec3> N_e =
        calder::mls_avg<vec3>(shell, N_f, x_e, 2.0 * _eps, 2.0);

    // std::vector<vec3> N_e = asawa::shell::edge_normals(shell, x);

    std::vector<mat3> F =
        calder::quadric_curvature(shell, x_e, N_e, 8.0 * _eps, 3.0);

    // std::vector<mat3> F = calder::covariant_frame(shell, x_e, 4.0 *
    // _eps, 3.0);

    std::vector<real> g_edge(__M->edge_count(), 0.0);
    auto range = shell.get_edge_range();

    for (auto c0 : range) {
      int i = shell.vert(c0);
      int j = shell.vert(shell.other(c0));

      vec3 xi = x[i];
      vec3 xj = x[j];
      vec3 xc = 0.5 * (xi + xj);

      mat3 Fi = F[c0 / 2];
      vec3 f0 = Fi.col(0);
      vec3 f1 = Fi.col(1);
      real cc = 0.01;
      // logger::line(xc - cc * f0, xc + cc * f0, vec4(1.0, 0.0, 0.0, 1.0));
      // logger::line(xc - cc * f1, xc + cc * f1, vec4(0.0, 0.0, 1.0, 1.0));

      real l0 = pow(f0.norm(), 1.0);
      real l1 = pow(f1.norm(), 1.0);

      f0.normalize();
      f1.normalize();

      real lt = l0 + l1;

      if (lt < 1.0e-8) {
        continue;
      }
      vec3 dx = xj - xi;
      vec3 e = dx.normalized();

      real t = 2.0 * (l0 / (l0 + l1) - 0.5);
      real C = pow(t, 1.0);

      real e0 = pow(abs(e.dot(f0)), 2.0);
      real e1 = pow(abs(e.dot(f1)), 2.0);

      // real g = va::mix(C, w[1] * e0, w[0] * e1);
      real g = va::mix(C, w[0] * e0 / l0, w[1] * e1 / l1);

      // logger::line(xc - g * dx, xc + g * dx, vec4(1.0, 0.0, 0.0, 1.0));
      //   real g = C * l1 * pow(e.dot(f1), 2.0);
      g_edge[c0 / 2] = g;
    }
    return g_edge;
  }

  std::vector<real> edge_anisotropic_weight(asawa::shell::shell &shell,
                                            const std::vector<real> &d,
                                            const std::vector<mat3> &F,

                                            const vec2 &w) {
    // anisotropic weights on edges

    const std::vector<vec3> &x = asawa::const_get_vec_data(shell, 0);

    index_t imin = get_new_origin(shell);

    std::vector<real> f(shell.vert_count(), 0.0);
    f[imin] = 1.0;

    auto range = shell.get_edge_range();
    std::vector<real> g_edge(__M->edge_count(), 0.0);

    for (auto c0 : range) {
      int vi = shell.vert(c0);
      int vj = shell.vert(shell.other(c0));
      int fi = shell.face(c0);
      int fj = shell.face(shell.other(c0));

      vec3 xi = x[vi];
      vec3 xj = x[vj];
      vec3 xc = 0.5 * (xi + xj);

      real di = d[vi];
      real dj = d[vj];
      real dij = 0.5 * (di + dj);

      mat3 Fi = F[fi];
      mat3 Fj = F[fj];
      mat3 Fij = 0.5 * (Fi + Fj);
      vec3 T = Fij.col(0);
      vec3 B = Fij.col(1);
      vec3 N = Fij.col(2);
      T = T.normalized();
      N = N.normalized();
      B = T.cross(N).normalized();

      vec3 e = (x[vj] - x[vi]).normalized();
      real gT = pow(e.dot(T), 2.0);
      real gB = pow(e.dot(B), 2.0);

#if 0
      real d1 = pow(dij, 2.0);
      real d0 = 1.0 - d1;
      vec2 wp = vec2(gT, gB).array() * w.array();
      real C = 0.01;
      logger::line(xc - C * wp[0] * T, xc + C * wp[0] * T,
                   vec4(1.0, 0.0, 0.0, 1.0));
      logger::line(xc - C * wp[1] * B, xc + C * wp[1] * B,
                   vec4(0.0, 0.0, 1.0, 1.0));
#endif
      g_edge[c0 / 2] = 1.0 * vec2(gT, gB).dot(w);
    }
    return g_edge;
  }

  std::vector<real> edge_rx_weights(asawa::shell::shell &shell) {
    const std::vector<vec3> &x = asawa::const_get_vec_data(shell, 0);
    auto range = shell.get_edge_range();
    std::vector<real> g_edge(__M->edge_count(), 0.0);

    std::vector<real> &rxa =
        std::dynamic_pointer_cast<reaction_diffusion>(_rx0)->get_rxa();
    std::vector<real> &rxb =
        std::dynamic_pointer_cast<reaction_diffusion>(_rx0)->get_rxb();

    for (auto c0 : range) {
      int i = shell.vert(c0);
      int j = shell.vert(shell.other(c0));

      real ra = rxa[i] + rxa[j];
      real rb = rxb[i] + rxb[j];
      real dra = rxa[i] - rxa[j];
      real drb = rxb[i] - rxb[j];
      real dgrx = 3.0 * abs(drb) - 1.0 * abs(dra);
      // real grx = (4.0 * rb - 1.0 * ra);

      g_edge[c0 / 2] = dgrx;
    }
    return g_edge;
  }

  std::vector<real> growth_weights(asawa::shell::shell &shell) {

    const std::vector<vec3> &x = asawa::const_get_vec_data(shell, 0);

    index_t imin = get_new_origin(shell);

    std::vector<real> f(shell.vert_count(), 0.0);
    f[imin] = 1.0;

    bontecou::laplacian L(__M, x);

    std::vector<real> d = L.heatDist(f, 0.2);
    std::vector<mat3> F = L.heatFrame(f, 0.2);

    real golden = (1.0 + sqrt(5.0)) / 2.0;

    //////weights on edges
    std::cout << "edge_anisotropic_weight" << std::endl;
    std::vector<real> w_aniso =
        edge_anisotropic_weight(shell, d, F, vec2(0.25, 1.5));
    std::cout << "covariant_growth_weight" << std::endl;
    std::vector<real> w_covariant =
        covariant_growth_weight(shell, vec2(0.75, 0.5));
    std::cout << "rx weights" << std::endl;
    std::vector<real> g_rx = edge_rx_weights(shell);
    ////////

    auto range = shell.get_edge_range();
    std::vector<real> g_edge(__M->edge_count(), 0.0);

    for (auto c0 : range) {
      int i = shell.vert(c0);
      int j = shell.vert(shell.other(c0));
      real di = d[i];
      real dj = d[j];
      real dij = 0.5 * (di + dj);

      real gc = w_covariant[c0 / 2];
      real ga = w_aniso[c0 / 2];
      real grx = g_rx[c0 / 2];

      real g = pow(dij, 1.5) * gc;
      vec3 xi = x[i];
      vec3 xj = x[j];
      vec3 xc = 0.5 * (xi + xj);
      vec3 dx = xj - xi;
#if 0
      real C = 0.25;
      if (g < 0.0) {
        logger::line(xc - C * g * dx, xc + C * g * dx,
                     vec4(1.0, 0.0, 0.0, 1.0));
      } else {
        logger::line(xc - C * g * dx, xc + C * g * dx,
                     vec4(0.0, 1.0, 0.0, 1.0));
      }
#endif
      g_edge[c0 / 2] = 1.0 + 0.01 * g;
    }
    return g_edge;
  }

  std::vector<vec3> covariant_forces(asawa::shell::shell &shell, vec3 w) {
    std::vector<real> g_geodesic = vertex_geodesic_weight(shell);
    const std::vector<vec3> &x = asawa::const_get_vec_data(shell, 0);
    std::vector<vec3> N_v = asawa::shell::vertex_normals(shell, x);
    std::vector<vec3> N_f = asawa::shell::face_normals(shell, x);
    // std::vector<vec3> N_vs =
    //     calder::mls_avg<vec3>(shell, N_f, x, 2.0 * _eps, 2.0);
    // std::vector<mat3> C =
    //    calder::taubin_curvature(shell, x, N_v, 2.0 * _eps, 3.0);
    std::vector<mat3> C =
        calder::quadric_curvature(shell, x, N_v, 6.0 * _eps, 3.0);

    // std::vector<mat3> C = calder::covariant_frame(shell, x, 2.0 * _eps, 3.0);

    std::vector<vec3> f(x.size(), vec3::Zero());
    auto range = shell.get_vert_range();

    for (auto i : range) {
      vec3 Ni = N_v[i];
      // vec3 Ni_s = N_vs[i];

      mat3 Ci = C[i];
      if (Ci.norm() < 1.0e-8)
        continue;
      vec3 c0 = Ci.col(0);
      vec3 c1 = Ci.col(1);
      vec3 c2 = Ci.col(2);

      vec3 s = vec3(c0.norm(), c1.norm(), c2.norm());
      c0.normalize();
      c1.normalize();
      c2.normalize();
      real st = s[0] + s[1] + s[2];

      real t0 = s[1] / s[0];
      real t1 = s[0] + s[1];

      real t2 = 1.0 + pow(t0, 4.0) - exp(pow(t0, 8.0));

      vec3 S = 10.0 * (t1 * t2) * Ni;

      if (S.norm() > 12.0 * _eps)
        S = 12.0 * _eps * S.normalized();
      if (S.hasNaN())
        S = vec3::Zero();

      f[i] = g_geodesic[i] * S;
      // f[i] = 0.01 * Ni / (1.0 - s[1] / s[0]);
      //  f[i] = (w[0] * s[0] / st + w[1] * s[1] / st) * Ni;
      //   f[i] = vec3::Zero();
      // logger::line(x[i], x[i] + 0.1 * f[i], vec4(1.0, 0.0, 1.0, 1.0));

      // logger::line(x[i], x[i] + 0.1 * s[0] * c0, vec4(1.0, 0.0, 0.0, 1.0));
      // logger::line(x[i], x[i] + 0.1 * s[1] * c1, vec4(0.0, 1.0, 0.0, 1.0));
      //  logger::line(x[i], x[i] + 1.0 * s[2] * c2, vec4(0.0, 0.0, 1.0, 1.0));
    }
    std::vector<vec3> f_f = asawa::shell::vert_to_face<vec3>(shell, x, f);
    std::vector<vec3> f_s = calder::mls_avg<vec3>(shell, f_f, x, 4.0 * _eps);
    return f_s;
  }

  std::vector<vec3> cylinder_forces(asawa::shell::shell &shell,
                                    vec2 w = vec2(1.0, 1.0)) {

    const std::vector<vec3> &x = asawa::const_get_vec_data(shell, 0);
    index_t imin = get_new_origin(shell);
    std::vector<real> fh(shell.vert_count(), 0.0);
    fh[imin] = 1.0;
    bontecou::laplacian L(__M, x);
    std::vector<real> d = L.heatDist(fh, 0.2);
    std::vector<mat3> F_f = L.heatFrame(fh, 0.2);
    std::vector<mat3> F_v = asawa::shell::face_to_vert<mat3>(shell, F_f);
    std::vector<vec3> N_v = asawa::shell::vertex_normals(shell, x);
    std::vector<gaudi::calder::vec7> M =
        calder::cylinder(shell, x, N_v, 8.0 * _eps, 3.0);

    std::vector<vec3> f(x.size(), vec3::Zero());
    auto range = shell.get_vert_range();

    for (auto i : range) {
      if (M[i].hasNaN())
        continue;
      vec3 f1 = F_v[i].col(0).normalized();
      vec3 cen = M[i].segment(0, 3);
      vec3 f0 = M[i].segment(3, 3);
      f0 = va::sgn(f0, f1) * f0;
      real r = M[i][6];

      vec3 dp = x[i] - cen;
      vec3 dpN = dp.normalized();
      real dpd = dp.norm();

      real dr = dpd - r;
      dr = dpd - r;

      real di = d[i];
      // logger::line(x[i], x[i] + 0.1 * di * f0, vec4(0.0, 1.0, 0.0, 1.0));
      f[i] = di * (w[0] * f0 - w[1] * dr * N_v[i]);
      // f[i] = d[i] * f0;
    }
    return f;
  }

  std::vector<vec3> bulk_force(asawa::shell::shell &shell) {
    std::vector<real> g_geodesic = vertex_geodesic_weight(shell);
    const std::vector<vec3> &x = asawa::const_get_vec_data(shell, 0);
    std::vector<vec3> N_v = asawa::shell::vertex_normals(shell, x);

    std::vector<vec3> f(x.size(), vec3::Zero());
    auto range = shell.get_vert_range();
    vec3 sun = 100.0 * vec3(0.5, 1.0, 0.0);

    for (auto i : range) {
      real d = g_geodesic[i];
      vec3 xi = x[i];
      vec3 dx_s = sun - xi;
      vec3 dx_o = xi - _origin;
      dx_o[1] *= 0.1;
      dx_o.normalize();
      dx_s.normalize();
      f[i] = d * (dx_s + dx_o);
    }
    return f;
  }

  std::vector<vec3> rx_forces(asawa::shell::shell &shell) {
    std::vector<real> g_geodesic = vertex_geodesic_weight(shell);
    auto range = shell.get_vert_range();

    const std::vector<vec3> &x = asawa::const_get_vec_data(shell, 0);

    std::vector<real> &rxa =
        std::dynamic_pointer_cast<reaction_diffusion>(_rx0)->get_rxa();
    std::vector<real> &rxb =
        std::dynamic_pointer_cast<reaction_diffusion>(_rx0)->get_rxb();
    std::vector<vec3> N = asawa::shell::vertex_normals(shell, x);

    for (auto i : range) {
      real d = g_geodesic[i];
      real d2 = pow(d, 3.0);
      real d12 = pow(d, 1.0);
      real ra = rxa[i];
      real rb = rxb[i];
      vec3 n = N[i];
      vec3 f = d2 * (rb - 0.15 * ra) * n;
      // vec3 f = d * (6.0 * rb - 0.5 * ra) * n;
      N[i] = 8.0 * f;
    }
    return N;
  }

  void calc_kf(const std::vector<real> &t) {
    real f0 = 0.034, k0 = 0.059; // fingerprints
    // real f1 = 0.032, k1 = 0.058; // fingerprints

    // real f1 = 0.025, k1 = 0.057;
    //  real f1 = 0.0531, k1 = 0.0626;

    // real f1 = 0.030, k1 = 0.060; // turing patterns
    // real f1 = 0.031, k1 = 0.06  ; // turing patterns
    // real f1 = 0.0571, k1 = 0.063; // mazy
    real f1 = 0.0257, k1 = 0.0555; // mazy

    //  real f1 = 0.034, k1 = 0.0618; // spots and worms
    _k = std::vector<real>(t.size(), 0.0);
    _f = std::vector<real>(t.size(), 0.0);
    for (int i = 0; i < t.size(); i++) {
      _k[i] = va::mix(t[i], k0, k1);
      _f[i] = va::mix(t[i], f0, f1);
    }
  }

  void step_rx(int frame) {

    const std::vector<real> d = vertex_geodesic_weight(*__M);
    calc_kf(d);

    int N = frame == 1 ? 10 : 10;
    for (int i = 0; i < N; i++) {
      std::dynamic_pointer_cast<reaction_diffusion>(_rx0)->step_anisotropic(
          16.0 * _h, _f, _k);
      //_rx0->step(20.0 * _h);
      //_rx1->step(20.0 * _h);
    }
  }

  void step_dynamics(int frame) {
    hepworth::block::projection_solver solver;

    std::vector<hepworth::projection_constraint::ptr> constraints;

    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    std::vector<vec3> &v = asawa::get_vec_data(*__M, 1);

    std::vector<vec3> M = asawa::shell::vertex_areas_3(*__M, x);
    std::vector<real> li = asawa::shell::edge_lengths(*__M, x);

    std::vector<vec3> Ns = asawa::shell::vertex_normals(*__M, x);
    std::vector<vec3> f(x.size(), vec3::Zero());

    std::vector<vec3> f0 = covariant_forces(*__M, vec3(1.0, 1.0, 1.0));
    std::vector<vec3> f1 = rx_forces(*__M);
    std::vector<vec3> f2 = cylinder_forces(*__M);
    // std::vector<vec3> f3 = bulk_force(*__M);

    // add all fi to f
    for (int i = 0; i < x.size(); i++) {
      f[i] += 0.1 * f0[i];
      f[i] += 0.2 * f1[i];
      f[i] += 0.25 * f2[i];
      // f[i] += 0.25 * f2[i];
      // f[i] += 0.1 * f3[i];
    }

    std::vector<real> g = growth_weights(*__M);

    for (int i = 0; i < li.size(); i++) {
      // std::cout << g[i] << " " << 1.0 / g[i] << std::endl;
      li[i] = g[i] * li[i];
    }

    hepworth::vec3_block::ptr X = hepworth::vec3_block::create(M, x, v, f);

    hepworth::block::init_edge_strain(*__M, constraints, x, li, 1.0e-1, {X});
    //   hepworth::block::init_edge_strain(*__M, constraints, x, 1.0e-2, {X});
    hepworth::block::init_bending(*__M, constraints, x, 9.5e-1, {X});
    hepworth::block::init_edge_willmore(*__M, constraints, 3.0e-1, {X});

    // hepworth::block::init_triangle_strain(*__M, constraints, x, 1.0e-1, {X});
    hepworth::block::init_area(*__M, constraints, x, 2.0e-1, {X}, false);
    real eps = 3.0 * __surf->_Cm;

    hepworth::block::init_pnt_tri_collisions(*__M, *__surf, constraints, x,
                                             3.0 * eps, 1.0 * eps, 1.0, {X, X});

    solver.set_constraints(constraints);

    std::vector<hepworth::sim_block::ptr> blocks = {X};
    solver.step(blocks, _h, 0.5, 30);
  }

  void step_f() {
    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    std::vector<vec3> f = covariant_forces(*__M, vec3(1.0, 1.0, 1.0));
    std::vector<vec3> f2 = cylinder_forces(*__M, vec2(0.5, 1.0));
    //  f = bulk_force(*__M);
    //    f = cylinder_forces(*__M);

    for (int i = 0; i < x.size(); i++) {
      x[i] += _h * (0.1 * f2[i] + 0.05 * f[i]);
    }

    smoothMesh(0.01, 200);
  }

  void smoothMesh(real C, int N) {

    vec3_datum::ptr x_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(0));
    std::vector<vec3> &x = x_datum->data();
    bontecou::laplacian3 M(__M, x, true);
    M.init();
    real cc = C / 100.0;
    for (int k = 0; k < N; k++) {
      std::cout << "." << std::flush;
      x = M.smooth(x, C - cc, C + cc);
      int i = 0;
      for (auto xi : x) {
        if (!std::isfinite(xi.dot(xi))) {
          std::cout << xi.transpose() << std::endl;
          __M->vprintv(i);
          i++;
        }
      }
    }
    // x_datum->data() = x;
    std::cout << "done!" << std::endl;
  }

  void step(int frame) {
    std::cout << "frame: " << frame << std::endl;
    step_rx(frame);
    step_dynamics(frame);
    // step_f();
    __surf->step(false);
  }

  module_base::ptr _rx0;
  // module_base::ptr _rx1;

  vec3 _origin;
  real _h = 0.1;
  real _eps = 0.1;
  std::vector<real> _f;
  std::vector<real> _k;

  shell::shell::ptr __M;
  shell::dynamic::ptr __surf;
};

} // namespace duchamp
} // namespace gaudi
#endif