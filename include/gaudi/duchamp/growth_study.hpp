#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "gaudi/vec_addendum.h"

#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/dynamic.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/asawa/shell/triangulation_debug.hpp"

#include "gaudi/asawa/shell/walk.hpp"
#include "gaudi/asawa/datums.hpp"

#include "gaudi/kusama/laplacian.hpp"
#include "gaudi/kusama/laplacian_anisotropic.hpp"
#include "gaudi/duchamp/modules/rx_diffuse.hpp"

#include "gaudi/hepworth/block/generic_constraints.hpp"
#include "gaudi/hepworth/block/generic_constraints_init.hpp"

#include "gaudi/hepworth/block/shell_constraints.hpp"
#include "gaudi/hepworth/block/shell_constraints_init.hpp"
#include "gaudi/hepworth/block/sim_block.hpp"
#include "gaudi/hepworth/block/solver.hpp"

#include "gaudi/asawa/primitive_objects.hpp"
#include "gaudi/calder/least_squares_fit.hpp"
#include "gaudi/calder/tangent_point_integrators.hpp"
#include "gaudi/common.h"
#include "gaudi/logger.hpp"

#include "modules/rx/ginzburg_landau.hpp"
#include "modules/rx/grey_scott.hpp"
#include "modules/rx_colormap.hpp"
#include "modules/rx/swift_hohenberg.hpp"
#include <Eigen/Sparse>
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <math.h>
#include <memory>
#include <limits>
#include <optional>
#include <vector>
#include "gaudi/geometry_logger.hpp"

#ifndef __M2REFACTOR_TEST__
#define __M2REFACTOR_TEST__

namespace gaudi {
namespace duchamp {

using namespace asawa;

/// Selects which reaction–diffusion module backs \ref growth_study.
/// Use `growth_study::create(growth_rx_model::grey_scott)` etc.
enum class growth_rx_model {
  grey_scott,
  swift_hohenberg,
  ginzburg_landau,
};

class growth_study {
public:
  typedef std::shared_ptr<growth_study> ptr;

  /// Default \p grey_scott matches pre–multi-model behavior (implicit GS substeps).
  /// Swift–Hohenberg / CGLE use explicit reaction substeps and need smaller per-step `h`
  /// (see \ref step_rx).
  static ptr create(growth_rx_model model = growth_rx_model::grey_scott) {
    return std::make_shared<growth_study>(model);
  }

  static bool should_trace_frame(index_t frame) {
    return frame < 3 || frame % 60 == 0;
  }

  /// When enabled, \ref step_dynamics logs the first non-finite value in key
  /// arrays (before sanitization) and a triangle/edge quality summary. Use
  /// `frame_stride` to probe every N frames only.
  void set_nan_probe(bool on, int frame_stride = 1) {
    _nan_probe = on;
    _nan_probe_stride = std::max(1, frame_stride);
  }

  explicit growth_study(growth_rx_model model = growth_rx_model::grey_scott)
      : _rx_model(model) {
    std::cerr << "[growth_study_core] ctor start" << std::endl;
    //__M = shell::load_cube();
    __M = shell::load_bunny();
    std::cerr << "[growth_study_core] mesh loaded" << std::endl;

    {
      const std::vector<vec3> &xpre = asawa::const_get_vec_data(*__M, 0);
      asawa::shell::dump_pre_triangulation_report(*__M, std::cerr, &xpre);
    }
    shell::triangulate(*__M);
    std::cerr << "[growth_study_core] triangulated" << std::endl;
    for (int i = 0; i < __M->face_count(); i++) {
      asawa::shell::FaceId fi = asawa::shell::face_id(i);
      if (__M->fbegin(fi) >= 0) {
        assert(__M->fsize(fi) == 3);
      }
    }

    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    asawa::center(x);

    /////////
    // dynamic surface
    /////////

    real l0 = 4.0 * asawa::shell::avg_length(*__M, x);
    _eps = l0;
    __surf = shell::dynamic::create(__M, 1.0 * l0, 2.5 * l0, 0.5 * l0);
    std::cerr << "[growth_study_core] shell dynamic created" << std::endl;

    std::cout << "[growth_study_core] stepping the surface for warmup: 0" << std::endl;
    __surf->step(true);
    std::cout << "[growth_study_core] stepping the surface for warmup: 1" << std::endl;
    __surf->step(true);
    std::cout << "[growth_study_core] stepping the surface for warmup: 2" << std::endl;
    __surf->step(true);
    std::cout << "[growth_study_core] stepping the surface for warmup: 3" << std::endl;
    std::cerr << "[growth_study_core] warmup steps done" << std::endl;
    // real f = 0.075, k = 0.0615;
    real f0 = 0.031, k0 = 0.0585;
    // real f = 0.04, k = 0.065;
    //  real f = 0.049, k = 0.0629;
    //  real f = 0.025, k = 0.0535;

    real da0 = 3.0e-4, db0 = 0.5 * da0;
    switch (_rx_model) {
    case growth_rx_model::grey_scott: {
      rx::grey_scott::ptr rx0 =
          rx::grey_scott::create(__M, f0, k0, da0, db0);
      _rx0 = std::dynamic_pointer_cast<module_base>(rx0);
      std::cerr << "[growth_study_core] reaction diffusion (Grey–Scott) created"
                  << std::endl;
      break;
    }
    case growth_rx_model::swift_hohenberg: {
      rx::swift_hohenberg::ptr sh = rx::swift_hohenberg::create(__M, 0.04, 1.0);
      _rx0 = std::dynamic_pointer_cast<module_base>(sh);
      std::cerr << "[growth_study_core] Swift–Hohenberg created" << std::endl;
      break;
    }
    case growth_rx_model::ginzburg_landau: {
      rx::ginzburg_landau::ptr gl = rx::ginzburg_landau::create(__M, 1.2, 1.0);
      _rx0 = std::dynamic_pointer_cast<module_base>(gl);
      std::cerr << "[growth_study_core] Ginzburg–Landau created" << std::endl;
      break;
    }
    }
    /*
    real f1 = 0.04, k1 = 0.065;
    real da1 = 2.00e-4, db1 = 0.5 * da1;
    rx::grey_scott::ptr rx1 =
        rx::grey_scott::create(__M, f1, k1, da1, db1);
    _rx1 = std::dynamic_pointer_cast<module_base>(rx1);
*/
    init_origin();
    std::cerr << "[growth_study_core] ctor complete" << std::endl;
  };

  void init_origin() {
    const std::vector<vec3> &x = asawa::const_get_vec_data(*__M, 0);

    std::array<vec3, 2> ext = asawa::extents(x);
    _origin = vec3(                    //
        0.5 * (ext[0][0] + ext[1][0]), //
        ext[0][1],                     //
        0.5 * (ext[0][2] + ext[1][2]));
  }

  /// Optional min–max color map \f$t\in[0,1]\to\f$ RGBA (defaults to cool–warm).
  void set_mesh_color_gradient(std::function<vec4(real)> g) {
    _mesh_color_gradient = std::move(g);
  }

  std::vector<vec4> get_mesh_colors() {
    auto default_grad = [](real t) {
      return vec4(t, 0.15 + 0.85 * (1.0 - t), 1.0 - 0.4 * t, 1.0);
    };
    std::function<vec4(real)> grad =
        _mesh_color_gradient ? *_mesh_color_gradient : default_grad;

    switch (_rx_model) {
    case growth_rx_model::grey_scott: {
      std::vector<vec4> colors(__M->vert_count(), vec4(1.0, 0.0, 0.0, 1.0));
      std::vector<real> &rx0a =
          std::dynamic_pointer_cast<rx::grey_scott>(_rx0)->get_rxa();
      std::vector<real> &rx0b =
          std::dynamic_pointer_cast<rx::grey_scott>(_rx0)->get_rxb();
      vec4 col_a(1.0, 0.0, 1.0, 1.0);
      vec4 col_b(0.0, 1.0, 1.0, 1.0);
      for (int k = 0; k < __M->vert_count(); k++) {
        colors[k] = 1.0 * rx0a[k] * col_a + 3.0 * rx0b[k] * col_b;
      }
      return colors;
    }
    case growth_rx_model::swift_hohenberg: {
      std::vector<real> u =
          std::dynamic_pointer_cast<rx::swift_hohenberg>(_rx0)->get_u();
      return field_colors_minmax(u, grad);
    }
    case growth_rx_model::ginzburg_landau: {
      auto gl = std::dynamic_pointer_cast<rx::ginzburg_landau>(_rx0);
      std::vector<real> mag = gl->amplitude_abs();
      return field_colors_minmax(mag, grad);
    }
    }
    return std::vector<vec4>(__M->vert_count(), vec4(1, 0, 0, 1));
  }
  vec3 get_origin() { return _origin; }

  index_t get_new_origin(asawa::shell::shell &shell) {
    //oof, we have acceleration structures we should use them.
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

    kusama::laplacian L(__M, x);
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
      asawa::shell::CornerId cid = asawa::shell::corner_id(c0);
      int i = shell.vert(cid);
      int j = shell.vert(shell.other(cid));

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
      //   geometry_logger::line(xc, cen, vec4(1.0, 0.0, 0.0, 1.0));
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
      asawa::shell::CornerId cid = asawa::shell::corner_id(c0);
      int i = shell.vert(cid);
      int j = shell.vert(shell.other(cid));

      vec3 xi = x[i];
      vec3 xj = x[j];
      vec3 xc = 0.5 * (xi + xj);

      mat3 Fi = F[c0 / 2];
      vec3 f0 = Fi.col(0);
      vec3 f1 = Fi.col(1);
      real cc = 0.01;
      // geometry_logger::line(xc - cc * f0, xc + cc * f0, vec4(1.0, 0.0, 0.0, 1.0));
      // geometry_logger::line(xc - cc * f1, xc + cc * f1, vec4(0.0, 0.0, 1.0, 1.0));

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

      // geometry_logger::line(xc - g * dx, xc + g * dx, vec4(1.0, 0.0, 0.0, 1.0));
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
      asawa::shell::CornerId cid = asawa::shell::corner_id(c0);
      int vi = shell.vert(cid);
      int vj = shell.vert(shell.other(cid));
      int fi = shell.face(cid);
      int fj = shell.face(shell.other(cid));

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
      geometry_logger::line(xc - C * wp[0] * T, xc + C * wp[0] * T,
                   vec4(1.0, 0.0, 0.0, 1.0));
      geometry_logger::line(xc - C * wp[1] * B, xc + C * wp[1] * B,
                   vec4(0.0, 0.0, 1.0, 1.0));
#endif
      g_edge[c0 / 2] = 1.0 * vec2(gT, gB).dot(w);
    }
    return g_edge;
  }

  std::vector<real> edge_rx_weights(asawa::shell::shell &shell) {
    auto range = shell.get_edge_range();
    std::vector<real> g_edge(__M->edge_count(), 0.0);

    switch (_rx_model) {
    case growth_rx_model::grey_scott: {
      std::vector<real> &rxa =
          std::dynamic_pointer_cast<rx::grey_scott>(_rx0)->get_rxa();
      std::vector<real> &rxb =
          std::dynamic_pointer_cast<rx::grey_scott>(_rx0)->get_rxb();
      for (auto c0 : range) {
        asawa::shell::CornerId cid = asawa::shell::corner_id(c0);
        int i = shell.vert(cid);
        int j = shell.vert(shell.other(cid));
        real dra = rxa[i] - rxa[j];
        real drb = rxb[i] - rxb[j];
        real dgrx = 3.0 * abs(drb) - 1.0 * abs(dra);
        g_edge[c0 / 2] = dgrx;
      }
      break;
    }
    case growth_rx_model::swift_hohenberg: {
      std::vector<real> u =
          std::dynamic_pointer_cast<rx::swift_hohenberg>(_rx0)->get_u();
      for (auto c0 : range) {
        asawa::shell::CornerId cid = asawa::shell::corner_id(c0);
        int i = shell.vert(cid);
        int j = shell.vert(shell.other(cid));
        g_edge[c0 / 2] = std::abs(u[i] - u[j]);
      }
      break;
    }
    case growth_rx_model::ginzburg_landau: {
      auto gl = std::dynamic_pointer_cast<rx::ginzburg_landau>(_rx0);
      const std::vector<real> &ru = gl->get_u();
      const std::vector<real> &rv = gl->get_v();
      for (auto c0 : range) {
        asawa::shell::CornerId cid = asawa::shell::corner_id(c0);
        int i = shell.vert(cid);
        int j = shell.vert(shell.other(cid));
        real ai = std::sqrt(ru[i] * ru[i] + rv[i] * rv[i]);
        real aj = std::sqrt(ru[j] * ru[j] + rv[j] * rv[j]);
        g_edge[c0 / 2] = std::abs(ai - aj);
      }
      break;
    }
    }
    return g_edge;
  }

  std::vector<real> growth_weights(asawa::shell::shell &shell) {

    const std::vector<vec3> &x = asawa::const_get_vec_data(shell, 0);

    index_t imin = get_new_origin(shell);

    std::vector<real> f(shell.vert_count(), 0.0);
    f[imin] = 1.0;

    kusama::laplacian L(__M, x);

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
      asawa::shell::CornerId cid = asawa::shell::corner_id(c0);
      int i = shell.vert(cid);
      int j = shell.vert(shell.other(cid));
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
        geometry_logger::line(xc - C * g * dx, xc + C * g * dx,
                     vec4(1.0, 0.0, 0.0, 1.0));
      } else {
        geometry_logger::line(xc - C * g * dx, xc + C * g * dx,
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
      // geometry_logger::line(x[i], x[i] + 0.1 * f[i], vec4(1.0, 0.0, 1.0, 1.0));

      // geometry_logger::line(x[i], x[i] + 0.1 * s[0] * c0, vec4(1.0, 0.0, 0.0, 1.0));
      // geometry_logger::line(x[i], x[i] + 0.1 * s[1] * c1, vec4(0.0, 1.0, 0.0, 1.0));
      //  geometry_logger::line(x[i], x[i] + 1.0 * s[2] * c2, vec4(0.0, 0.0, 1.0, 1.0));
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
    kusama::laplacian L(__M, x);
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
      // geometry_logger::line(x[i], x[i] + 0.1 * di * f0, vec4(0.0, 1.0, 0.0, 1.0));
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
    std::vector<vec3> N = asawa::shell::vertex_normals(shell, x);

    switch (_rx_model) {
    case growth_rx_model::grey_scott: {
      std::vector<real> &rxa =
          std::dynamic_pointer_cast<rx::grey_scott>(_rx0)->get_rxa();
      std::vector<real> &rxb =
          std::dynamic_pointer_cast<rx::grey_scott>(_rx0)->get_rxb();
      for (auto i : range) {
        real d = g_geodesic[i];
        real d2 = pow(d, 3.0);
        real ra = rxa[i];
        real rb = rxb[i];
        vec3 n = N[i];
        vec3 f = d2 * (rb - 0.15 * ra) * n;
        N[i] = 8.0 * f;
      }
      break;
    }
    case growth_rx_model::swift_hohenberg: {
      std::vector<real> u =
          std::dynamic_pointer_cast<rx::swift_hohenberg>(_rx0)->get_u();
      for (auto i : range) {
        real d = g_geodesic[i];
        real d2 = pow(d, 3.0);
        vec3 n = N[i];
        real ui = std::isfinite(u[i]) ? u[i] : 0.0;
        vec3 f = d2 * ui * n;
        N[i] = 8.0 * f;
      }
      break;
    }
    case growth_rx_model::ginzburg_landau: {
      auto gl = std::dynamic_pointer_cast<rx::ginzburg_landau>(_rx0);
      const std::vector<real> &ru = gl->get_u();
      const std::vector<real> &rv = gl->get_v();
      for (auto i : range) {
        real d = g_geodesic[i];
        real d2 = pow(d, 3.0);
        real rui = std::isfinite(ru[i]) ? ru[i] : 0.0;
        real rvi = std::isfinite(rv[i]) ? rv[i] : 0.0;
        real amp = std::sqrt(rui * rui + rvi * rvi);
        vec3 n = N[i];
        vec3 f = d2 * (amp - 0.15 * rui) * n;
        N[i] = 8.0 * f;
      }
      break;
    }
    }
    return N;
  }

  void calc_sh_params(const std::vector<real> &t) {
    real e0 = 0.02, e1 = 0.07;
    real g0 = 0.8, g1 = 1.6;
    real l0 = 2.0e-4, l1 = 7.0e-4;
    _eps_sh.resize(t.size());
    _g_sh.resize(t.size());
    _lam_sh.resize(t.size());
    for (int i = 0; i < static_cast<int>(t.size()); i++) {
      _eps_sh[i] = va::mix(t[i], e0, e1);
      _g_sh[i] = va::mix(t[i], g0, g1);
      _lam_sh[i] = va::mix(t[i], l0, l1);
    }
  }

  void calc_gl_params(const std::vector<real> &t) {
    real a0 = 0.8, a1 = 1.8;
    real b0 = 0.6, b1 = 1.4;
    _alpha_gl.resize(t.size());
    _beta_gl.resize(t.size());
    for (int i = 0; i < static_cast<int>(t.size()); i++) {
      _alpha_gl[i] = va::mix(t[i], a0, a1);
      _beta_gl[i] = va::mix(t[i], b0, b1);
    }
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
    if (should_trace_frame(frame)) {
      std::cerr << "[growth_study_core] step_rx begin frame=" << frame
                << std::endl;
    }

    const std::vector<vec3> &x_pos = asawa::const_get_vec_data(*__M, 0);
    _C_aniso = kusama::build_curvature_aligned_laplacian(
        *__M, x_pos, _aniso_vd_lambda, _aniso_sigma_u, _aniso_sigma_v,
        asawa::shell::face_curvature_stencil::one_ring, _aniso_kind);

    if (_rx_model == growth_rx_model::grey_scott) {
      std::dynamic_pointer_cast<rx::grey_scott>(_rx0)->set_diffuse_smooth(
          std::make_optional<rx_smooth_fn>([this](std::vector<real> &f, real dt) {
            const std::vector<vec3> &xp = asawa::const_get_vec_data(*__M, 0);
            rx_diffuse_scalar_implicit_cotan(__M, xp, f, dt, &_C_aniso);
          }));
    } else if (_rx_model == growth_rx_model::swift_hohenberg) {
      std::dynamic_pointer_cast<rx::swift_hohenberg>(_rx0)->set_diffuse_smooth(
          std::make_optional<rx_smooth_fn>([this](std::vector<real> &f, real dt) {
            const std::vector<vec3> &xp = asawa::const_get_vec_data(*__M, 0);
            rx_diffuse_scalar_implicit_cotan(__M, xp, f, dt, &_C_aniso);
          }));
    } else {
      std::dynamic_pointer_cast<rx::ginzburg_landau>(_rx0)->set_dispersive_stiffness(
          &_C_aniso);
    }

    const std::vector<real> d = vertex_geodesic_weight(*__M);
    calc_kf(d);
    calc_sh_params(d);
    calc_gl_params(d);

    // Grey--Scott: implicit reaction + implicit diffusion per substep.
    // Swift--Hohenberg / GL: one macro step per frame (reaction + implicit
    // diffusion) after operator-split refactor; h chosen from old 40x substeps.
    const int N_gs = 10;
    const real h_gs = 16.0 * _h;
    const int N_sh = 1;
    const real h_sh = 0.01; // was 40 * 2.5e-4
    const int N_gl = 1;
    const real h_gl = 0.006; // was 40 * 1.5e-4

    switch (_rx_model) {
    case growth_rx_model::grey_scott:
      for (int i = 0; i < N_gs; i++) {
        std::dynamic_pointer_cast<rx::grey_scott>(_rx0)->step_anisotropic(
            h_gs, _f, _k, nullptr);
      }
      break;
    case growth_rx_model::swift_hohenberg:
      for (int i = 0; i < N_sh; i++) {
        std::dynamic_pointer_cast<rx::swift_hohenberg>(_rx0)->step(
            h_sh, _eps_sh, _g_sh, _lam_sh, nullptr);
      }
      break;
    case growth_rx_model::ginzburg_landau:
      for (int i = 0; i < N_gl; i++) {
        std::dynamic_pointer_cast<rx::ginzburg_landau>(_rx0)->step(
            h_gl, _alpha_gl, _beta_gl, nullptr);
      }
      break;
    }
    if (should_trace_frame(frame)) {
      std::cerr << "[growth_study_core] step_rx end" << std::endl;
    }
  }

  void step_dynamics(int frame) {
    if (should_trace_frame(frame)) {
      std::cerr << "[growth_study_core] step_dynamics begin frame=" << frame
                << std::endl;
    }
    hepworth::block::projection_solver solver;

    std::vector<hepworth::projection_constraint::ptr> constraints;

    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    std::vector<vec3> &v = asawa::get_vec_data(*__M, 1);

    if (_nan_probe && (frame % _nan_probe_stride == 0)) {
      nan_probe_vec3("x(raw)", x);
      nan_probe_vec3("v(raw)", v);
    }

    // Break NaN/Inf cascades into hepworth: integrate_inertia uses f → s in
    // projection_solver (b = M*s + A^T*p); bending normals need finite q.
    for (index_t ii = 0; ii < static_cast<index_t>(x.size()); ++ii) {
      if (!x[ii].allFinite())
        x[ii] = vec3::Zero();
      if (!v[ii].allFinite())
        v[ii] = vec3::Zero();
    }

    std::vector<vec3> M = asawa::shell::vertex_areas_3(*__M, x);
    const real m_min = 1e-14;
    for (auto &mi : M) {
      if (!mi.allFinite() || mi[0] < m_min)
        mi = vec3(m_min, m_min, m_min);
    }
    std::vector<real> li = asawa::shell::edge_lengths(*__M, x);
    if (_nan_probe && (frame % _nan_probe_stride == 0)) {
      nan_probe_scalar("li_edge_lengths(raw)", li);
      nan_probe_tris_edges(frame, x, li);
    }

    std::vector<vec3> Ns = asawa::shell::vertex_normals(*__M, x);
    std::vector<vec3> f(x.size(), vec3::Zero());

    std::vector<vec3> f0 = covariant_forces(*__M, vec3(1.0, 1.0, 1.0));
    std::vector<vec3> f1 = rx_forces(*__M);
    std::vector<vec3> f2 = cylinder_forces(*__M);
    if (_nan_probe && (frame % _nan_probe_stride == 0)) {
      nan_probe_vec3("f0_covariant", f0);
      nan_probe_vec3("f1_rx", f1);
      nan_probe_vec3("f2_cylinder", f2);
    }
    if (should_trace_frame(frame)) {
      std::cerr << "[growth_study_core] forces ready" << std::endl;
    }
    // std::vector<vec3> f3 = bulk_force(*__M);

    // add all fi to f
    for (int i = 0; i < x.size(); i++) {
      f[i] += 0.1 * f0[i];
      f[i] += 0.2 * f1[i];
      f[i] += 0.25 * f2[i];
      // f[i] += 0.25 * f2[i];
      // f[i] += 0.1 * f3[i];
    }

    const real f_cap = 1e4 * std::max(_eps, real(1e-6));
    for (auto &fi : f) {
      if (!fi.allFinite()) {
        fi = vec3::Zero();
        continue;
      }
      real fn = fi.norm();
      if (fn > f_cap && fn > 0.0)
        fi *= f_cap / fn;
    }

    std::vector<real> g = growth_weights(*__M);
    if (_nan_probe && (frame % _nan_probe_stride == 0))
      nan_probe_scalar("g_growth_weight(raw)", g);
    if (should_trace_frame(frame)) {
      std::cerr << "[growth_study_core] growth weights ready" << std::endl;
    }

    for (real &gi : g) {
      if (!std::isfinite(gi) || gi <= 0.0)
        gi = 1.0;
      gi = std::max(0.2, std::min(5.0, gi));
    }
    for (int i = 0; i < li.size(); i++) {
      // std::cout << g[i] << " " << 1.0 / g[i] << std::endl;
      li[i] = g[i] * li[i];
      if (!std::isfinite(li[i]) || li[i] <= 0.0)
        li[i] = 1e-8;
    }

    hepworth::vec3_block::ptr X = hepworth::vec3_block::create(M, x, v, f);
    if (should_trace_frame(frame)) {
      std::cerr << "[growth_study_core] sim block created" << std::endl;
    }

    hepworth::block::init_edge_strain(*__M, constraints, x, li, 1.0e-1, {X});
    //   hepworth::block::init_edge_strain(*__M, constraints, x, 1.0e-2, {X});
    hepworth::block::init_bending(*__M, constraints, x, 9.5e-1, {X});
    hepworth::block::init_edge_willmore(*__M, constraints, 3.0e-1, {X});

    // hepworth::block::init_triangle_strain(*__M, constraints, x, 1.0e-1, {X});
    hepworth::block::init_area(*__M, constraints, x, 2.0e-1, {X}, false);
    real eps = 3.0 * __surf->_Cm;

    hepworth::block::init_pnt_tri_collisions(*__M, *__surf, constraints, x,
                                             3.0 * eps, 1.0 * eps, 1.0, {X, X});
    if (should_trace_frame(frame)) {
      std::cerr << "[growth_study_core] constraints ready count="
                << constraints.size() << std::endl;
    }

    solver.set_constraints(constraints);

    std::vector<hepworth::sim_block::ptr> blocks = {X};
    solver.step(blocks, _h, 0.5, 30);
    if (should_trace_frame(frame)) {
      std::cerr << "[growth_study_core] step_dynamics end" << std::endl;
    }
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
    kusama::laplacian3 M(__M, x, true);
    M.init();
    real cc = C / 100.0;
    for (int k = 0; k < N; k++) {
      std::cout << "." << std::flush;
      x = M.smooth(x, C - cc, C + cc);
      int i = 0;
      for (auto xi : x) {
        if (!std::isfinite(xi.dot(xi))) {
          std::cout << xi.transpose() << std::endl;
          __M->vprintv(asawa::shell::vert_id(i));
          i++;
        }
      }
    }
    // x_datum->data() = x;
    std::cout << "done!" << std::endl;
  }

  void step(int frame) {
    if (should_trace_frame(frame)) {
      std::cerr << "[growth_study_core] step begin frame=" << frame
                << std::endl;
    }
    step_rx(frame);
    step_dynamics(frame);
    // step_f();
    __surf->step(false);
    if (should_trace_frame(frame)) {
      std::cerr << "[growth_study_core] shell dynamic step end" << std::endl;
    }
  }

  module_base::ptr _rx0;
  growth_rx_model _rx_model;
  std::optional<std::function<vec4(real)>> _mesh_color_gradient;

  vec3 _origin;
  real _h = 0.1;
  real _eps = 0.1;
  std::vector<real> _f;
  std::vector<real> _k;
  std::vector<real> _eps_sh, _g_sh, _lam_sh;
  std::vector<real> _alpha_gl, _beta_gl;

  /// Curvature-aligned anisotropic cotan stiffness; rebuilt each @ref step_rx.
  kusama::laplacian::sparmat _C_aniso;
  real _aniso_vd_lambda = 3.0;
  real _aniso_sigma_u = 1.0;
  real _aniso_sigma_v = 0.001;
  kusama::anisotropic_laplacian_kind _aniso_kind =
      kusama::anisotropic_laplacian_kind::fem_d;

  shell::shell::ptr __M;
  shell::dynamic::ptr __surf;

private:
  bool _nan_probe = false;
  int _nan_probe_stride = 1;

  void nan_probe_vec3(const char *tag, const std::vector<vec3> &a) const {
    if (!_nan_probe)
      return;
    for (index_t i = 0; i < static_cast<index_t>(a.size()); ++i) {
      if (!a[i].allFinite()) {
        std::cerr << "[growth_study nan_probe] " << tag << " first bad vert "
                  << i << " -> " << a[i].transpose() << std::endl;
        return;
      }
    }
  }

  void nan_probe_scalar(const char *tag, const std::vector<real> &a) const {
    if (!_nan_probe)
      return;
    for (index_t i = 0; i < static_cast<index_t>(a.size()); ++i) {
      if (!std::isfinite(a[i])) {
        std::cerr << "[growth_study nan_probe] " << tag << " first bad idx "
                  << i << " -> " << a[i] << std::endl;
        return;
      }
    }
  }

  /// Degenerate / NaN faces and near-zero edges (bad triangulation / inverted
  /// elements show up here before hepworth sees them).
  void nan_probe_tris_edges(int frame, const std::vector<vec3> &x,
                            const std::vector<real> &li) const {
    if (!_nan_probe || (frame % _nan_probe_stride) != 0)
      return;

    std::vector<real> fa = asawa::shell::face_areas(*__M, x);
    index_t badf = 0;
    real amin = std::numeric_limits<real>::infinity();
    for (index_t fi = 0; fi < static_cast<index_t>(fa.size()); ++fi) {
      if (!std::isfinite(fa[fi]) || fa[fi] <= 1e-20) {
        if (badf < 12)
          std::cerr << "[growth_study nan_probe] face " << fi
                    << " area=" << fa[fi] << std::endl;
        ++badf;
      } else {
        amin = std::min(amin, fa[fi]);
      }
    }
    std::cerr << "[growth_study nan_probe] frame=" << frame
              << " faces=" << fa.size() << " bad_face_area=" << badf
              << " min_pos_face_area="
              << (std::isfinite(amin) ? amin : std::numeric_limits<real>::quiet_NaN())
              << std::endl;

    index_t bade = 0;
    real lmin = std::numeric_limits<real>::infinity();
    for (index_t e = 0; e < static_cast<index_t>(li.size()); ++e) {
      if (!std::isfinite(li[e]) || li[e] <= 1e-20) {
        if (bade < 12)
          std::cerr << "[growth_study nan_probe] edge " << e
                    << " len=" << li[e] << std::endl;
        ++bade;
      } else {
        lmin = std::min(lmin, li[e]);
      }
    }
    std::cerr << "[growth_study nan_probe] edges=" << li.size()
              << " bad_edge_len=" << bade << " min_pos_edge_len="
              << (std::isfinite(lmin) ? lmin : std::numeric_limits<real>::quiet_NaN())
              << std::endl;
  }
};

} // namespace duchamp
} // namespace gaudi
#endif