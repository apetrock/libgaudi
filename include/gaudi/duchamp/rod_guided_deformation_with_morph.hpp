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

#include "gaudi/asawa/shell/walk.hpp"

#include "gaudi/bontecou/laplacian.hpp"
#include "gaudi/calder/least_squares_fit.hpp"
#include "gaudi/calder/tangent_point_integrators.hpp"

#include "gaudi/asawa/primitive_objects.hpp"
#include "gaudi/common.h"

#include "modules/ccd.hpp"
#include "modules/knotted_surface.hpp"
#include "modules/mighty_morphin.hpp"
#include "utils/point_things.hpp"
#include "utils/string_things.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <memory>
#include <set>
#include <vector>
#include <zlib.h>

#ifndef __ROD_GUIDED_WITH_MORPH__
#define __ROD_GUIDED_WITH_MORPH__

namespace gaudi {
namespace duchamp {

int ROT = 600;

struct walk_config {
  walk_config(index_t i0 = 0, index_t N_steps = 4000, real thet = 0.0, //
              bool rotate = false, vec2 cr = vec2(5.0, -8.0),
              bool align = false, vec4 ca = vec4(1.0, 0.5, 0.0, 0.0))
      : i0(i0), N_steps(N_steps), thet(thet), rotate(rotate), align(align) {
    this->ca = ca;
    this->cr = cr;
  };
  index_t i0 = 0;
  index_t N_steps = 4000;
  real thet = 0.0;
  bool align = false;
  vec4 ca = vec4(5.0, -8.0, 0.0, 0.0);
  bool rotate = false;
  vec2 cr = vec2(1.0, 0.5);
};

const int N_walk_configs = 6;

int iw0 = 2;

int iw1 = (iw0 + 1) % N_walk_configs;
int iwp = (iw0 + N_walk_configs - 1) % N_walk_configs;
int Nw = 24000;
walk_config _wc[N_walk_configs] = {
    walk_config(100, Nw, 0.021, false, vec2(0.0, 0.0), false,
                vec4(0.0, 0.0, 0.0, 0.0)),
    walk_config(100, Nw, 0.55, false, vec2(0.0, 0.0), false,
                vec4(0.0, 0.0, 0.0, 0.0)), // jennifer_0
    walk_config(100, Nw, 0.23, true, vec2(3.0, -5.0), false,
                vec4(0.0, 0.6, 0.0, 0.0)),
    walk_config(100, Nw, 0.6, true, vec2(6.0, -8.0), true,
                vec4(1.0, 0.1, 0.1, 0.25)),
    walk_config(100, Nw, 0.265, false, vec2(0.0, 0.0), true,
                vec4(1.0, 0.2, 0.3, 0.00)),
    walk_config(100, Nw, 2.5, false, vec2(0.0, 0.0), true,
                vec4(0.5, 0.0, 0.65, 0.22))};
// walk_config(0, 4000, 0.0, false, vec2(5.0, -8.0), false,
// vec4(1.0, 0.5, 0.0, 0.0)),
// walk_config(0, 4000, 0.0, false, vec2(5.0, -8.0), false,
// vec4(1.0, 0.5, 0.0, 0.0))};

inline vec3 hsv_mix(real t, vec3 a, vec3 b) {

  vec3 ah = va::rgb_to_hsv(a);
  vec3 bh = va::rgb_to_hsv(b);

  real ahr = 2.0 * M_PI * ah[0] / 360.0;
  real bhr = 2.0 * M_PI * bh[0] / 360.0;
  real hr = va::polar_mix(t, ahr, bhr);

  real h = 360.0 * hr / (2.0 * M_PI);
  real s = va::mix(t, ah[1], bh[1]);
  real v = va::mix(t, ah[2], bh[2]);

  std::cout << va::hsv_to_rgb(vec3(h, s, v)).transpose() << std::endl;

  return va::hsv_to_rgb(vec3(h, s, v));
}

inline vec3 lab_mix(real t, vec3 a, vec3 b) {
  // this wrong, replace with good code
  vec3 aLab = va::rgb_to_lab(a);
  vec3 bLab = va::rgb_to_lab(b);
  vec3 lerp_lab = va::mix(t, aLab, bLab);
  return va::lab_to_rgb(lerp_lab);
}

inline std::array<vec3, 2> _get_colors(index_t frame) {
  //
  const vec3 colors[6][2] = {
      {vec3(0.74855, 0.006349, 0.617022), vec3(0.251846, 1.0849, 0.291323)},
      {vec3(0.54710, 0.774614, 0.0197271), vec3(0.621809, 0.0433273, 0.993413)},
      {vec3(0.36217, 0.105595, 0.858886), vec3(0.789299, 0.0025459, 0.881501)},
      {vec3(1.03731, 0.44225, 0.124936), vec3(0.113354, 0.377527, 0.904618)},
      {vec3(0.19887, 0.263064, 0.919388), vec3(0.848492, 0.57854, 0.091643)},
      {vec3(0.02987, 0.76066, 0.489041), vec3(0.848492, 0.57854, 0.091643)}}; //

  const vec3 *cP = colors[iwp];
  const vec3 *c0 = colors[iw0];
  const vec3 *cN = colors[iw1];
  std::array<vec3, 2> out = {vec3::Zero(), vec3::Zero()};
  if (frame < 600) {
    real t = real(frame) / 600.0;
    // out[0] = hsv_mix(t, cP[0], c0[0]);
    // out[1] = hsv_mix(t, cP[1], c0[1]);
    out[0] = lab_mix(t, cP[0], c0[0]);
    out[1] = lab_mix(t, cP[1], c0[1]);

  } else if (frame >= 600 && frame < 2400) {
    out = std::array<vec3, 2>{c0[0], c0[1]};
  } else if (frame >= 2400) {
    real t = real(frame - 2400) / 600.0;
    out[0] = va::mix(t, c0[0], cN[0]);
    out[1] = va::mix(t, c0[1], cN[1]);
  }
  out[0] = colors[iw0][0];
  out[1] = colors[iw0][1];

  return out;
};

using namespace asawa;

///////////////////////////////////////
class rod_guided_deformation_with_morph {
public:
  typedef std::shared_ptr<rod_guided_deformation_with_morph> ptr;

  static ptr create() {
    return std::make_shared<rod_guided_deformation_with_morph>();
  }

  rod_guided_deformation_with_morph() {
    //__M = load_cube();
    //__M = shell::load_bunny();
    __M = shell::load_obj("assets/dennis.obj");
    //__M = shell::load_obj("assets/jennifer_0.obj");
    //__M = shell::load_obj("assets/joel_0.obj");
    //__M = shell::load_obj("assets/hand.obj");
    //__M = shell::load_obj("assets/washington.obj");
    //__M = shell::load_obj("assets/messer.obj");
    

    //__M = shell::load_crab();

    shell::triangulate(*__M);
    for (int i = 0; i < __M->face_count(); i++) {
      if (__M->fbegin(i) > 0) {
        assert(__M->fsize(i) == 3);
      }
    }

    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    asawa::center(x, 2.0);

    /////////
    // dynamic surface
    /////////
    real l0 = asawa::shell::avg_length(*__M, x);
    // real C = 0.5;
    real C = 2.0;
    C = 1.0;
    __surf = shell::dynamic::create(__M, C * l0, 2.5 * C * l0, 0.75 * C * l0);

    /////////////////////
    // Rod
    /////////////////////
    walk_config wc0 = _wc[iw0];
    std::vector<vec3> x_w =                             //
        silly_walk(*__M, wc0.thet, wc0.i0, wc0.N_steps, //
                   wc0.rotate, wc0.cr,                  //
                   wc0.align, wc0.ca, _eps);

    walk_config wc1 = _wc[iw1];
    _target =                                           //
        silly_walk(*__M, wc1.thet, wc1.i0, wc1.N_steps, //
                   wc1.rotate, wc1.cr,                  //
                   wc1.align, wc1.ca, _eps);

    __R = rod::rod::create(x_w, false);
    //__R->_update_frames(normals);

    real lavg = 1.5 * l0;
    __R->_r = 2.5 * lavg;
    __Rd = rod::dynamic::create(__R, 0.75 * lavg, 2.25 * lavg, 0.25 * lavg);
    _eps = l0;
    #if 0
    std::cout << "stepping: " << std::endl;
    for (int i = 0; i < 1; i++) {
      __surf->step();
      __Rd->step();
    }
#endif
    _knotted_surface = knotted_surface_module::create(__M, __surf, __R, __Rd);
    if (iw0 == 2) {
      _knotted_surface->add_angle_constraint(vec3(0.0, 0.0, 1.0), 0.15 * M_PI,
                                             0.1);
      _knotted_surface->add_angle_constraint(vec3(0.0, 1.0, 0.0), 0.10 * M_PI,
                                             0.1);
    } else if (iw0 == 3) {
      _knotted_surface->add_angle_constraint(vec3(0.0, 0.0, 1.0), 0.025 * M_PI,
                                             0.1);
      _knotted_surface->add_angle_constraint(vec3(0.0, 1.0, 0.0), 0.037 * M_PI,
                                             0.1);
    } else if (iw0 == 5) {
      _knotted_surface->add_angle_constraint(vec3(0.0, 0.0, 1.0), 0.05 * M_PI,
                                             0.1);
    }

    std::vector<index_t> face_vert_ids = __M->get_face_vert_ids();
    _morphin = mighty_morphin::create(__M);
    _morphin->add_geometry(x, face_vert_ids, 1.0, vec3(0.0, 0.0, 0.0));
    _ccd = continuous_collision_detection::create(__M, __surf, 0.1 * l0);
  };
  
  std::array<vec3, 2> get_colors(int frame) { return _get_colors(frame); }

  std::vector<vec3> rod_target(const std::vector<vec3> &target) {
    const std::vector<vec3> &x_r = __R->x();
    const std::vector<real> &t = __R->t();

    std::vector<vec3> x_t(x_r.size(), vec3::Zero());

    real Nr = real(x_r.size());
    real Nt = real(target.size());
    rod::rod::ptr trod = rod::rod::create(target, false);

    for (int i = 0; i < x_r.size(); i++) {
      real ti = t[i];
      real it = ti * real(Nt);
      index_t i0 = index_t(floor(it));
      index_t i1 = index_t(i0 + 1.0);
      real t = it - real(i0);
      x_t[i] = va::mix(0.5 * t, target[i0], target[i1]);
    }
    return x_t;
  }

  std::vector<vec3> tween_rod(const std::vector<vec3> &target) {

    const std::vector<vec3> &x_r = __R->x();
    const std::vector<real> &t = __R->t();

    std::vector<vec3> x_t = rod_target(target);
#if 0
    for (int i = 0; i < target.size() - 1; i++) {
      logger::line(target[i], target[i + 1], vec4(0.0, 1.0, 1.0, 1.0));
    }
#endif
    std::vector<vec3> dx(x_t.size(), vec3::Zero());
    for (int i = 0; i < x_t.size(); i++) {
      // logger::line(x_t[i], x_r[i], vec4(0.0, 1.0, 1.0, 1.0));
      dx[i] = (x_t[i] - x_r[i]);
    }
    return dx;
  }
  /*
    void test_quadric_sdf() {
      const std::vector<vec3> &x_s = asawa::get_vec_data(*__M, 0);
      const std::vector<vec3> &x_r = __R->x();

      std::vector<vec3> x = createPoints(200000, 0.5);
      real eps = 0.5 * _knotted_surface->get_eps();
      std::vector<vec3> Nr = _knotted_surface->get_rod_normals(*__R, *__M, eps);

      std::vector<real> sdf = calder::quadric_sdf(*__R, Nr, x, eps);

      for (int i = 0; i < x.size(); i++) {
        gg::geometry_logger::line(x[i], x[i] + 1e-4 * vec3(1.0, 1.0, 1.0),
                                  1000.0 * gg::sdf4(sdf[i]));
      }
    }
  */
  std::vector<vec3> calc_quadric_normal_flow(real cN, real cQ, real cS,
                                             real t) {
    asawa::shell::shell &M = *__M;
    asawa::rod::rod &R = *__R;
    asawa::shell::dynamic &Md = *__surf;
    asawa::rod::dynamic &Rd = *__Rd;

    std::vector<vec3> x_s = asawa::get_vec_data(M, 0);
    std::vector<vec3> x_s_f = asawa::shell::face_centers(M, x_s);
    std::vector<vec3> N_s_f = asawa::shell::face_normals(*__M, x_s);
    const std::vector<vec3> &x_r = R.x();

    std::vector<vec3> Nr = _knotted_surface->get_rod_normals(R, M, cN * _eps);

    std::vector<real> Q =
        calder::quadric_sdf(R, Nr, x_s_f, N_s_f, cQ * _eps, 3.0);

#if 1
    int i = 0;
    for (vec3 &N : N_s_f) {
      N *= -Q[i];
      i++;
    }
#endif
    std::vector<vec3> Nss =
        calder::mls_avg<vec3>(*__M, N_s_f, x_s, cS * _eps, 2.0);
    return Nss;
  }

  std::vector<vec3> calc_quadric_grad(real cN, real cQ, real cS) {
    asawa::shell::shell &M = *__M;
    asawa::rod::rod &R = *__R;
    asawa::shell::dynamic &Md = *__surf;
    asawa::rod::dynamic &Rd = *__Rd;

    std::vector<vec3> x_s = asawa::get_vec_data(M, 0);
    std::vector<vec3> x_s_f = asawa::shell::face_centers(M, x_s);
    std::vector<vec3> N_s_f = asawa::shell::face_normals(*__M, x_s);
    const std::vector<vec3> &x_r = R.x();

    std::vector<vec3> Nr = _knotted_surface->get_rod_normals(R, M, cN * _eps);

    std::vector<vec3> Q =
        calder::quadric_grad(R, Nr, x_s_f, N_s_f, cQ * _eps, 3.0);
    std::vector<vec3> Nss = calder::mls_avg<vec3>(*__M, Q, x_s, cS * _eps, 2.0);
#if 0
    int i = 0;
    for (vec3 &N : Nss) {
      gg::geometry_logger::line(x_s[i], x_s[i] + 1.0 * N,
                                vec4(1.0, 1.0, 0.0, 1.0));
      i++;
    }
#endif
    return Nss;
  }

  std::vector<vec3> calc_cyclide_normal_flow(real cN, real cQ) {
    asawa::shell::shell &M = *__M;
    asawa::rod::rod &R = *__R;
    asawa::shell::dynamic &Md = *__surf;
    asawa::rod::dynamic &Rd = *__Rd;

    std::vector<vec3> x_s = asawa::get_vec_data(M, 0);
    std::vector<vec3> N_s_v = asawa::shell::vertex_normals(*__M, x_s);
    std::vector<vec3> x_s_f = asawa::shell::face_centers(M, x_s);
    std::vector<vec3> N_s_f = asawa::shell::face_normals(*__M, x_s);
    const std::vector<vec3> &x_r = R.x();

    std::vector<vec3> Nr = _knotted_surface->get_rod_normals(R, M, cN * _eps);

    //    std::vector<real> Q =
    //        calder::darboux_cyclide_sdf(R, Nr, x_s_f, N_s_f, cQ * _eps, 3.0);

    std::vector<real> Q =
        calder::darboux_cyclide_sdf_quadric_centers(R, Nr, x_s, N_s_v, cQ * _eps, 3.0);

#if 1
    int i = 0;
    for (vec3 &N : N_s_v) {
      N *= -Q[i];
      // gg::geometry_logger::line(x_s[i], x_s[i] + 1.0e0 * N,
      //                           vec4(0.6, 0.0, 0.8, 1.0));

      i++;
    }
#endif

    return N_s_v;
  }

#if 1
  std::vector<vec3> calc_tangent_point_gradient_from_shell() {

    real eps = _knotted_surface->get_eps();
    std::vector<vec3> &x_s = asawa::get_vec_data(*__M, 0);
    std::vector<vec3> &x = __R->x();
    std::vector<real> l = __R->l0();
    std::vector<vec3> T = __R->N2c();
    std::vector<vec3> xc = __R->xc();

    std::vector<vec3> g0 =
        calder::tangent_point_gradient(*__R, x_s, l, T, 2.0 * eps, 6.0);
    for (int i = 0; i < g0.size(); i++) {
      // gg::geometry_logger::line(x[i], x[i] + 1.0e-7 * g0[i],
      //                           vec4(0.6, 0.0, 0.8, 1.0));
      g0[i] *= -1.0;
    }
    return g0;
  }
#endif

#if 1
  std::vector<vec3> calc_tangent_point_gradient() {

    real eps = _knotted_surface->get_eps();
    std::vector<vec3> &x = __R->x();
    std::vector<real> l = __R->l0();
    std::vector<vec3> T = __R->N2c();
    std::vector<vec3> xc = __R->xc();

    std::vector<vec3> g0 =
        calder::tangent_point_gradient(*__R, x, l, T, 2.0 * eps, 6.0);
    for (int i = 0; i < g0.size(); i++) {
      // gg::geometry_logger::line(x[i], x[i] + 1.0e-7 * g0[i],
      //                           vec4(0.6, 0.0, 0.8, 1.0));
      g0[i] *= -1.0;
    }
    return g0;
  }
#endif

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

  void morph_rod_step(real h, real t) {

    real t_s = va::smoothstep<real>(t);
    std::vector<vec3> &x_r = __R->x();
    std::vector<vec3> &x_s = asawa::get_vec_data(*__M, 0);
    __R->update_lengths();
    std::vector<real> &l_r = __R->l0();

    real h_step = 0.0;
    std::vector<vec3> tween_force = tween_rod(_target);
    std::vector<vec3> x_t = rod_target(_target);

    std::vector<vec3> x_tt(x_r);
    for (int i = 0; i < x_r.size(); i++) {
      x_tt[i] = va::mix(0.1 * t, x_r[i], x_t[i]);
    }

    hepworth::block::projection_solver solver;
    std::vector<hepworth::projection_constraint::ptr> constraints;
    std::vector<vec3> f(__R->__v.size(), vec3::Zero());
    std::vector<vec3> fr = calc_tangent_point_gradient();
    for (int i = 0; i < x_r.size(); i++) {
      f[i] += (1.0 - t_s) * 2e-6 * fr[i];
    }

    // capture before
    std::vector<vec3> x_r_0 = x_r;
    hepworth::vec3_block::ptr Xr =
        hepworth::vec3_block::create(__R->__M, __R->__x, __R->__v, f);
    hepworth::quat_block::ptr Ur =
        hepworth::quat_block::create(__R->__J, __R->__u, __R->__o);

    hepworth::block::init_smooth(*__R, constraints, 0.2, {Xr});
    hepworth::block::init_stretch_shear(*__R, constraints, l_r, 1e-3, {Xr, Ur});
    hepworth::block::init_bend_twist(*__R, constraints, 1e-3, {Ur});
    // hepworth::block::init_collisions(*__R, *__Rd, constraints, 1.0, {Xr,
    // Xr},
    //                                  1.0);
    hepworth::block::init_pinned(*__R, x_tt, constraints, x_r, 6e-2, {Xr});

    solver.set_constraints(constraints);
    std::vector<hepworth::sim_block::ptr> blocks = {Xr, Ur};
    solver.step(blocks, h);
    // difference

    std::vector<vec3> dx_r = x_r;
    for (int i = 0; i < x_r.size(); i++) {
      dx_r[i] = x_r[i] - x_r_0[i];
    }
    // splat the tween onto the surface
    std::vector<vec3> dxr_smooth =
        calder::mls_avg<vec3>(*__R, dx_r, x_s, 3.0 * _eps, 2.0);
    for (int i = 0; i < dxr_smooth.size(); i++) {
      // gg::geometry_logger::line(x_s[i], x_s[i] + dxr_smooth[i],
      //                           vec4(0.6, 0.0, 0.8, 1.0));
      x_s[i] += (1.0 - t_s) * dxr_smooth[i];
    }
#if 1
    int its = 2;
    _ccd->set_kmin(0.05);
    for (int k = 0; k < its; k++) {

      _ccd->init_step(0.1);
      // std::vector<vec3> quad_force_0 = calc_quadric_grad(1.0, 0.5, 2.0);
      //_ccd->add_shell_force(quad_force_0, (1.0 - t));
      // std::vector<vec3> quad_force_1 =
      //    calc_quadric_normal_flow(1.0, 1.0, 2.0, 0.5);
      // std::vector<vec3> tp_force = calc_tangent_point_gradient_from_shell();
      //_ccd->add_shell_force(quad_force_1, (1.0 - t_s) / _h);
      //_ccd->add_shell_force(tp_force, (1.0 - t_s) / _h);
      //_ccd->add_shell_force(tp_force, (1.0 - t_s) * 1e-6);
      // std::vector<vec3> quad_force_2 =
      // calc_cyclide_normal_flow(1.0, 6.0, 2.0);
      //_ccd->add_shell_force(quad_force_2, (1.0 - t_s) / _h);

      if (t > 0.5) {
        _morphin->step(h);
        std::vector<vec3> morph_force = _morphin->forces();
        real t1 = 2.0 * (t_s - 0.5);
        _ccd->add_shell_force(morph_force, t1);
      }
      _ccd->step(5e-3);
    }
#endif
    smoothMesh(0.01, 200);
  }

  void morph_step(real h, real t) {
#if 1
    int its = 1;
    _ccd->set_kmin(0.25);
    for (int k = 0; k < its; k++) {
      _morphin->step(h);
      std::vector<vec3> morph_force = _morphin->forces();
      _ccd->init_step(0.1);
      _ccd->add_shell_force(morph_force, 1.0);
      _ccd->step(1e-3);
    }
#endif
  }

  void step(int frame) {
    std::cout << "frame: " << frame << std::endl;
    std::cout << "  -surface" << std::endl;
    std::cout << "    -corners: " << __M->corner_count() << std::endl;
    std::cout << "    -verts: " << __M->vert_count() << std::endl;
    std::cout << "    -faces: " << __M->face_count() << std::endl;
    std::cout << "  -curve" << std::endl;
    std::cout << "    -verts: " << __R->x().size() << std::endl;

    auto calc_frame_t = [](index_t frame, index_t start_frame,
                           index_t end_frame) {
      return real(frame - start_frame) / real(end_frame - start_frame);
    };

    frame = frame;
    // walk(__surf->_Cc);
    if (frame < 3000) {

      if (frame >= 600 && frame < 1800) {

        _knotted_surface->init_step(_h);
        real t = calc_frame_t(frame, 600, 1800);
        real offset = 1.0;
        offset = min(1.0 + 15.0 * t, 20.0);
        //_knotted_surface->set_helicity_constraint(true);
        //_knotted_surface->set_helicity_weight(1.0e-1);
        _knotted_surface->set_rod_offset(offset);
        _knotted_surface->set_willmore_weight(5.0e-1);
        _knotted_surface->set_area_weight(5.0e-2);
        _knotted_surface->set_shell_strain_weight(1.0e-1);
        _knotted_surface->set_shell_bending_weight(1.0e-1);

        _knotted_surface->set_rod_pin_weight(1.0e-3);
        _knotted_surface->set_rod_strain_weight(1.0e-1);
        _knotted_surface->set_rod_bending_weight(1.0e-1);
        //_knotted_surface->set_rod_offset(10.0);
        // calc_torus_gradient();
        _knotted_surface->add_rod_force(calc_tangent_point_gradient(), 3.0e-7);
        //_knotted_surface->add_shell_force(calc_quadric_grad(0.5, 2.0, 4.0), 4.0);
        /*
        _knotted_surface->add_shell_force(
            calc_quadric_normal_flow(1.0, 6.0, 2.0, 0.5), 0.5 / _h);
        */

        _knotted_surface->add_shell_force(calc_cyclide_normal_flow(0.5, 1.0),
                                          1.0 / _h);

        _knotted_surface->step(_h);

        for (int k = 0; k < 1; k++) {
          __surf->step(true);
          __Rd->step();
        }
      }
    }

    if (frame > 2400)
      exit(0);

    // step_sdf(frame);
  }
  // std::map<index_t, index_t> _rod_adjacent_edges;
  std::vector<vec3> _target;
  knotted_surface_module::ptr _knotted_surface;
  mighty_morphin::ptr _morphin;
  continuous_collision_detection::ptr _ccd;

  real _h = 0.025;
  real _eps = 0.1;
  shell::shell::ptr __M;
  shell::dynamic::ptr __surf;
  rod::rod::ptr __R;
  rod::dynamic::ptr __Rd;
};

} // namespace duchamp
} // namespace gaudi
#endif