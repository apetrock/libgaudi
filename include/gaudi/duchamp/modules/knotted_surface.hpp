
#ifndef __DUCHAMP_ROD_CONTROL_MODULE__
#define __DUCHAMP_ROD_CONTROL_MODULE__

#include "Eigen/src/Geometry/AngleAxis.h"
#include "gaudi/asawa/datums.hpp"

#include "gaudi/hepworth/block/generic_constraints.hpp"
#include "gaudi/hepworth/block/generic_constraints_init.hpp"

#include "gaudi/hepworth/block/rod_constraints.hpp"
#include "gaudi/hepworth/block/rod_constraints_init.hpp"

#include "gaudi/hepworth/block/shell_constraints.hpp"
#include "gaudi/hepworth/block/shell_constraints_init.hpp"

#include "gaudi/calder/rod_integrators.hpp"
#include "gaudi/calder/shell_integrators.hpp"
#include "gaudi/hepworth/block/coupling_collisions_init.hpp"
#include "gaudi/hepworth/block/sim_block.hpp"
#include "gaudi/hepworth/block/solver.hpp"

#include "gaudi/common.h"
#include "gaudi/logger.hpp"
#include "module_base_shell.hpp"
#include <algorithm>
#include <array>
#include <vector>

namespace gaudi {

namespace duchamp {

class knotted_surface_module : public module_base {
public:
  DEFINE_CREATE_FUNC(knotted_surface_module)
  knotted_surface_module(asawa::shell::shell::ptr &M,
                         asawa::shell::dynamic::ptr &D, asawa::rod::rod::ptr &R,
                         asawa::rod::dynamic::ptr &Rd)
      : __M(M), __surf(D), __R(R), __Rd(Rd) {

    D->set_flip_pred([&](asawa::shell::shell &M, const index_t &c0) {
      index_t c1 = M.other(c0);
      return _adjacent.find(c0) == _adjacent.end() &&
             _adjacent.find(c1) == _adjacent.end();
    });

    const std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    _eps = asawa::shell::avg_length(*__M, x);
  };

  std::vector<vec3> get_rod_normals(asawa::rod::rod &R, asawa::shell::shell &M,
                                    real eps) {
    std::vector<vec3> &xr = __R->__x;
    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    std::vector<vec3> Nf = asawa::shell::face_normals(M, x);
    std::vector<vec3> Nr = calder::mls_avg<vec3>(M, Nf, xr, eps, 2.0);
    for (int i = 0; i < Nr.size(); i++) {
      vec3 T = R.dir(i);
      vec3 B = Nr[i].cross(T);
      Nr[i] = T.cross(B);
    }
    return Nr;
  }

  std::vector<real> calc_conservative_collisions(asawa::rod::rod &R,
                                                 asawa::shell::shell &M,
                                                 real eps, int N_spread = 4) {
    std::vector<vec3> x_s = asawa::get_vec_data(M, 0);

    // TODO:, these assume that everything is tightly packed, this is wrong
    // assumption if knots get more complicated
    std::vector<vec3> xr = R.x_infill(N_spread);
    std::vector<index_t> rverts(xr.size());
    for (int i = 0; i < rverts.size(); i++) {
      rverts[i] = i;
    }
    vector<std::array<index_t, 2>> nearest =
        __surf->get_pnt_tri_collisions(rverts, rverts, xr, M, eps);

    std::vector<real> dist0(xr.size(), 0.0);

    std::cout << "x.size() " << R.x().size() << std::endl;

    std::cout << "xr.size() " << xr.size() << std::endl;
    std::cout << "rverts.size() " << rverts.size() << std::endl;
    real lavg = R.lavg();
    for (auto &c : nearest) {
      auto consec = R.consec(c[0] / N_spread);
      if (consec[2] < 0)
        continue;

      if (c[1] < 0) {
        dist0[c[0]] = 4.0 * lavg;
        continue;
      }
#if 0
      index_t ivr = c[0];
      index_t ifs = c[1];

      vec3 xri = xr[ivr];
      vec3 xf = asawa::shell::face_pnt(xri, M, ifs, x_s);
      gg::geometry_logger::line(xri, xf, vec4(1.0, 0.0, 1.0, 1.0));
#endif
    }

    std::vector<real> dist(R.x().size(), 0.0);
    for (int i = 0; i < dist.size(); i++) {
      for (int k = 0; k < N_spread; k++) {
        dist[i] = std::max(dist0[i * N_spread + k], dist[i]);
      }
    }

    return dist;
  }

  std::vector<real> calc_dist_1(asawa::rod::rod &R, asawa::shell::shell &M,
                                real eps, int N_spread = 4) {

    const std::vector<vec3> &xr = R.x();
    std::vector<real> dist = calc_conservative_collisions(R, M, eps, N_spread);

    auto assign = [](index_t ip, index_t im, const std::vector<vec3> &x,
                     std::vector<real> &dist) {
      dist[ip] = std::min(dist[ip], dist[im] + (x[ip] - x[im]).norm());
    };
    // sweep forward

    std::vector<index_t> rverts = R.get_ordered_verts();

    for (int k = 0; k < 2; k++) {
      // this should only need one it, why not working?

      for (int i = 0; i < rverts.size(); i++) {
        index_t ip = R.next(rverts[i]);
        index_t im = i;
        if (ip < 0)
          continue;
        assign(ip, im, xr, dist);
      }

      for (int i = rverts.size() - 1; i > -1; i--) {
        index_t ip = R.next(rverts[i]);
        index_t im = i;
        if (ip < 0)
          continue;
        assign(im, ip, xr, dist);
      }
    }

    for (int k = 0; k < 8; k++) {
      int i0 = k % 2;
      for (int i = i0; i < rverts.size(); i += 2) {
        auto cons = R.consec(rverts[i]);
        if (cons[2] < 0)
          continue;
        index_t ip = cons[2];
        index_t i0 = cons[1];
        index_t im = cons[0];

        dist[i0] = 0.25 * dist[im] + 0.5 * dist[i0] + 0.25 * dist[ip];
      }
    }
#if 0
    std::vector<vec3> Nr = get_rod_normals(*__R, *__M, 4.0 * eps);

    for (int i : rverts) {
      auto cons = R.consec(i);
      index_t i0 = cons[1];
      if (cons[2] < 0)
        continue;

      gg::geometry_logger::line(xr[i0],
                                xr[i0] + 1.0 * dist[i0] * Nr[i0].normalized(),
                                vec4(0.0, 1.0, 1.0, 1.0));
    }
#endif
    return dist;
  }

  std::vector<real> calc_rod_dist_grad(asawa::rod::rod &R,
                                       asawa::shell::shell &M, real eps,
                                       int N_spread = 4) {
    std::vector<real> dist = calc_dist_1(R, M, eps, N_spread);
    std::vector<real> g_d(dist.size(), 0.0);
    real lavg = R.lavg();
    std::vector<vec3> xr = R.x();
    for (int i = 0; i < dist.size(); i++) {
      auto idx = R.consec(i);
      index_t im = idx[0];
      index_t ip = idx[2];
      vec3 xr1 = xr[ip];
      vec3 xr0 = xr[im];
      real di = (xr1 - xr0).norm();
      di = std::max(di, lavg);
      real ddi = dist[ip] - dist[im];
      ddi = va::sgn(ddi) * std::min(abs(ddi), lavg);
      g_d[i] = ddi / di;
      // g_d[i] = va::sgn(g_d[i]) * std::min(abs(g_d[i]), 2.0);
    }
#if 0

    std::vector<index_t> rverts = R.get_ordered_verts();
    std::vector<vec3> Nr = get_rod_normals(R, *__M, 4.0 * eps);

    for (int i : rverts) {
      auto cons = R.consec(i);
      index_t i0 = cons[1];
      if (cons[2] < 0)
        continue;

      gg::geometry_logger::line(xr[i0],
                                xr[i0] + 0.1 * g_d[i0] * Nr[i0].normalized(),
                                vec4(0.0, 1.0, 1.0, 1.0));
    }
#endif
    // g_d = __R->vert_avg(g_d);
    return g_d;
  }

  void init_rod_shell_weld(
      asawa::rod::rod &R, asawa::rod::dynamic &rod_d, asawa::shell::shell &M,
      asawa::shell::dynamic &shell_d, const std::vector<vec3> &fm,
      std::vector<hepworth::projection_constraint::ptr> &constraints,
      const real &wr, const real &ws, real eps,
      std::vector<hepworth::sim_block::ptr> blocks) {

    const std::vector<vec3> &x0 = R.x();
    std::vector<vec3> x1 = asawa::get_vec_data(M, 0);

    std::vector<index_t> edge_verts_R = R.get_edge_vert_ids();
    std::vector<index_t> edge_map_R = R.get_edge_map();

    std::vector<vec3> Nr = get_rod_normals(R, M, eps);
    std::vector<vec3> Nr0 = R.N1();

#if 1
    std::vector<index_t> edges = M.get_edge_range();
    std::vector<index_t> edge_verts_M = M.get_edge_vert_ids();
    std::vector<index_t> edge_map_M = M.get_edge_map();

    // auto g_d = calc_rod_dist_grad(R, M, 2.0 * eps);

    auto g_d = calc_rod_dist_grad(R, M, 0.5 * eps, 4);
    _willmore_mask = std::vector<real>(x1.size(), 0.0);
    vector<std::array<index_t, 4>> sr_collisions =
        rod_d.get_collisions(edge_verts_M, x1, 1.0 * eps);

    for (auto &c : sr_collisions) {
      if (c[0] < 0)
        continue;
      if (c[1] < 0)
        continue;
      index_t vs0 = c[0];
      index_t vs1 = c[1];
      index_t vr0 = c[2];
      index_t vr1 = c[3];
      // real gr = 3e1 * eps;
      // if (abs(g_d[vr0]) > gr || abs(g_d[vr1]) > gr)
      //   continue;

      if (vr1 < 0 || vr1 < 0 || vs0 < 0 || vs1 < 0)
        continue;

      vec3 xr0 = x0[vr0];
      vec3 xr1 = x0[vr1];
      vec3 xs0 = x1[vs0];
      vec3 xs1 = x1[vs1];

      if ((xr1 - xr0).norm() < 1e-8)
        continue;

      if ((xs1 - xs0).norm() < 1e-8)
        continue;

      std::array<real, 3> d = va::distance_Segment_Segment(xr0, xr1, xs0, xs1);
      real g_di = va::mix(d[1], g_d[vr0], g_d[vr1]);

      if (abs(g_di) > 1e-3) {
        _willmore_mask[vs0] = 1.0;
        _willmore_mask[vs1] = 1.0;
      }

      vec3 xr = va::mix(d[1], xr0, xr1);
      vec3 xs = va::mix(d[2], xs0, xs1);
      vec3 dr = xr1 - xr0;
      vec3 dx = xr - xs;
      vec3 Ns0 = asawa::shell::vert_normal(M, vs0, x1);
      vec3 Ns1 = asawa::shell::vert_normal(M, vs1, x1);

      vec3 Ns = va::mix(d[2], Ns0, Ns1);

      real is_perp = pow(Ns.dot(dx.normalized()), 2.0);
      // if (8.0 * d[0] * (1.0 - is_perp) < eps) {
      if ((is_perp > 0.75 && d[0] < 1.0 * eps)) {

        // gg::geometry_logger::line(xs0, xs1, vec4(0.0, 0.0, 1.0, 1.0));

        real lr = (xr1 - xr0).norm();
        vec3 Nri0 = Nr[vr0].normalized();
        vec3 Nri1 = Nr[vr1].normalized();
        vec3 Nri = va::mix(d[1], Nri0, Nri1);

        vec3 Nr0i0 = Nr0[vr0].normalized();
        vec3 Nr0i1 = Nr0[vr1].normalized();
        vec3 Nr0i = va::mix(d[1], Nr0i0, Nr0i1);
        /*
        real lbound = 1e-8;
        real ubound = 1.0 - lbound;
        if (d[2] < lbound)
          continue;
        if (d[2] > ubound)
          continue;
    */
        // gg::geometry_logger::line(xs, xs + 1.0 * g_di * dr,
        //                           vec4(0.5, 0.5, 1.0, 1.0));
        // gg::geometry_logger::line(xs, xs + 0.1 * g_di * Nri,
        //                          vec4(0.5, 0.5, 1.0, 1.0));
        _adjacent.insert(M.find_edge_from_verts(vs0, vs1));
        // gg::geometry_logger::line(xs0, xs1, vec4(0.0, 0.0, 1.0, 1.0));

        hepworth::block::edge_edge_weld::ptr constraint =
            hepworth::block::edge_edge_weld::create(
                std::vector<index_t>({vr0, vr1, vs0, vs1}), wr, ws, blocks);
        constraint->set_slide(d[1] + 8.0 * g_di);

        // constraint->set_rotate_to(0.001, Nri, Nr0i);
        constraints.push_back(constraint);
      }
    }
#endif
  }

  std::vector<real> get_dist_rod(const std::vector<index_t> &vert_ids,
                                 const std::vector<vec3> &x,
                                 const asawa::rod::rod &R,
                                 const asawa::rod::dynamic &Rd) {
    std::vector<vec3> &xr = __R->x();
    std::vector<real> df(x.size(), 0.0);
    vector<std::array<index_t, 3>> sr_collisions =
        __Rd->get_vert_collisions(vert_ids, x, 99999.9);
    for (auto &c : sr_collisions) {
      index_t vs = c[0];
      index_t vr0 = c[1];
      index_t vr1 = c[2];

      if (vr1 < 0 || vr1 < 0 || vs < 0)
        continue;
      vec3 xr0 = xr[vr0];
      vec3 xr1 = xr[vr1];
      vec3 xs = x[vs];
      df[vs] = va::distance_from_line(xr0, xr1, xs);
    }
    return df;
  }
#if 0
  void init_weighted_willmore(
      asawa::shell::shell &M, asawa::rod::rod &R,
      std::vector<hepworth::projection_constraint::ptr> &constraints,
      const real &w_min, const real &w_max,
      std::vector<hepworth::sim_block::ptr> blocks) {

    std::vector<vec3> &xv = asawa::get_vec_data(*__M, 0);
    std::vector<vec3> xf = asawa::shell::face_centers(*__M, xv);
    std::vector<vec3> xe = asawa::shell::edge_centers(*__M, xv);

    std::vector<vec3> Ne = asawa::shell::edge_normals(*__M, xv);
    real eps = 0.5 * __surf->_Cc;

    std::vector<real> dist = calc_dist_1(R, M, eps, 4);
    std::transform(dist.begin(), dist.end(), dist.begin(),
                          [](real x) { return x > eps ? 1.0 : 0.0; });

    std::vector<real> df = calder::mls_avg<real>(R, dist, xe, eps, 2.0);
    // scale by w_min/w_max x = x +
    std::transform(df.begin(), df.end(), df.begin(),
                   [w_min, w_max](real x) {
                     return w_min + (w_max - w_min) * x;
                   });

#if 0
    for (int i = 0; i < df.size(); i++) {
      gg::geometry_logger::line(xe[i], xe[i] + 0.1 * df[i] * Ne[i],
                                vec4(0.0, 0.5, 1.0, 1.0));
    }
#endif

    hepworth::block::init_edge_willmore(M, constraints, df, blocks);
  }
#endif
  void init_weighted_area(
      const asawa::shell::shell &M,
      std::vector<hepworth::projection_constraint::ptr> &constraints,
      const real &w, std::vector<hepworth::sim_block::ptr> blocks) {

    std::vector<vec3> &xv = asawa::get_vec_data(*__M, 0);
    std::vector<vec3> xf = asawa::shell::face_centers(*__M, xv);
    std::vector<vec3> Nf = asawa::shell::face_normals(*__M, xv);
    std::vector<vec3> &xr = __R->x();
    std::vector<index_t> verts_F = M.get_face_range();

    real eps = 2.0 * __surf->_Cc;
    std::vector<real> df = get_dist_rod(verts_F, xf, *__R, *__Rd);

#if 0
    for (int i = 0; i < df.size(); i++) {
      gg::geometry_logger::line(xf[verts_F[i]],
                                xf[verts_F[i]] + 0.1 * df[verts_F[i]] * Nf[i],
                                vec4(0.0, 0.5, 1.0, 1.0));
    }
#endif
    std::transform(df.begin(), df.end(), df.begin(),
                   [w, eps](double x) { return max(x - 2.5 * eps, 0.0); });
    auto [min_it, max_it] = std::minmax_element(df.begin(), df.end());
    real vmin = *min_it;
    real vmax = *max_it;
    std::transform(df.begin(), df.end(), df.begin(), [w, vmin, vmax](real x) {
      return w * (x - vmin) / (vmax - vmin);
    });

    hepworth::block::init_area(M, constraints, xv, df, blocks, true);
    // return g;
  }

  void assert_nan(index_t k) {
    std::cout << k << std::endl;
    for (int i = 0; i < __R->__u.size(); ++i) {
      if (__R->__u[i].coeffs().hasNaN()) {
        std::cout << "nan at " << i << std::endl;
        exit(0);
      }
    }
  }

  void init_step(real h) {
    std::vector<vec3> &xs = asawa::get_vec_data(*__M, 0);
    std::vector<vec3> &xr = __R->x();
    _pr = std::vector<vec3>(0, vec3::Zero());
    _fs = std::vector<vec3>(xs.size(), vec3::Zero());
    _fr = std::vector<vec3>(xr.size(), vec3::Zero());
  }

  void set_rod_pin(const std::vector<vec3> &x, const real &w = 1e-1) {
    _pin_rod = true;
    _pr = x;
    _config.w_rod_pin = w;
  }

  void add_rod_force(const std::vector<vec3> &fr, const real &h = 1.0) {
    for (int i = 0; i < fr.size(); i++) {
      _fr[i] += h * fr[i];
    }
  }

  void add_shell_force(const std::vector<vec3> &fs, const real &h = 1.0) {
    for (int i = 0; i < fs.size(); i++) {
      _fs[i] += h * fs[i];
    }
  }

  void step(real h) {
    real eps = 0.5 * _eps;

    hepworth::block::projection_solver solver;

    std::vector<vec3> &xs = asawa::get_vec_data(*__M, 0);
    std::vector<vec3> &v = asawa::get_vec_data(*__M, 1);

    std::vector<vec3> M = asawa::shell::vertex_areas_3(*__M, xs);

    std::vector<real> &lr = __R->l0();
    std::vector<vec3> &xr = __R->x();

    __R->update_lengths();

    std::vector<real> li = asawa::shell::edge_lengths(*__M, xs);

    std::vector<vec3> fs(xs.size(), vec3::Zero());
    std::vector<vec3> fr(xr.size(), vec3::Zero());

    std::vector<hepworth::projection_constraint::ptr> constraints;
    //      std::vector<vec3> f = compute_ribbon_charge();

    // fs = calc_ribbon_sdf();

    hepworth::vec3_block::ptr Xs = hepworth::vec3_block::create(M, xs, v, _fs);

    hepworth::vec3_block::ptr Xr =
        hepworth::vec3_block::create(__R->__M, __R->__x, __R->__v, _fr);
    hepworth::quat_block::ptr Ur =
        hepworth::quat_block::create(__R->__J, __R->__u, __R->__o);

    std::vector<hepworth::sim_block::ptr> blocks;
    blocks.push_back(Xs);
    blocks.push_back(Xr);
    blocks.push_back(Ur);

    std::cout << "init weld" << std::endl;
#if 1

    init_rod_shell_weld(*__R, *__Rd, //
                        *__M, *__surf, fs,
                        constraints,          //
                        _config.w_rod_weld,   //
                        _config.w_shell_weld, //
                        4.0 * eps, {Xr, Xs});
#endif

#if 1
    if (_helicity_constraint)
      for (int i = 0; i < lr.size(); i++) {
        lr[i] *= 1.01;
      }
#endif
    std::cout << "main constraints" << std::endl;
    hepworth::block::init_stretch_shear(*__R, constraints, lr,
                                        _config.w_rod_strain, {Xr, Ur});
    hepworth::block::init_bend_twist(*__R, constraints, _config.w_rod_bending,
                                     {Ur});
    if (_helicity_constraint) {
      hepworth::block::init_helicity(*__R, constraints, _config.w_helicity,
                                     {Xr});
    }

    for (int i = 0; i < _angle_constraints.size(); i++) {
      real theta = _angle_constraints[i].theta;
      real weight = _angle_constraints[i].weight;
      vec3 axis = _angle_constraints[i].axis;
      hepworth::block::init_angle(*__R, constraints, axis, theta, weight, {Ur});
    }

    // we could add a list of targets for the rod to hit
    if (_pin_rod) {
      if (_pr.size() == 0) {
        hepworth::block::init_pinned(*__R, constraints, xr, _config.w_rod_pin,
                                     {Xr});
      } else {
        hepworth::block::init_pinned(*__R, _pr, constraints, xr,
                                     _config.w_rod_pin, {Xr});
      }
    }

    std::cout << "main collisions" << std::endl;

    if (_repel_rods) {
      hepworth::block::init_collisions(*__R, *__Rd, constraints, 1.0, {Xr, Xr},
                                       _config.rod_offset);
    }
#if 1
    std::cout << "shell strain" << std::endl;

    hepworth::block::init_triangle_strain(*__M, constraints, xs,
                                          _config.w_shell_strain, {Xs});
#endif
#if 1
    std::cout << "shell bending" << std::endl;
    hepworth::block::init_bending(*__M, constraints, xs,
                                  _config.w_shell_bending, {Xs});
#endif
    std::cout << "shell willmore" << std::endl;

#if 1
    hepworth::block::init_edge_willmore(*__M, constraints, _config.w_willmore,
                                        blocks);
#endif
#if 1
    std::cout << "shell area" << std::endl;
    init_weighted_area(*__M, constraints, _config.w_area, blocks);
#endif
    if (_shell_collisions) {
      std::cout << "init pnt trie collisions" << std::endl;
      hepworth::block::init_pnt_tri_collisions(
          *__M, *__surf, constraints, xs, 0.5 * eps, 0.5 * eps, 1.0, {Xs, Xs});
    }

    solver.set_constraints(constraints);

    // copy initial state
    std::vector<vec3> xs0 = xs;
    std::vector<vec3> xr0 = xr;

    solver.step(blocks, h, 0.5, 10);

    // final state is mix of initial and final
    real t = 0.5; // make this a parameter
#if 1
    for (int i = 0; i < xr.size(); i++)
      xr[i] = va::mix(t, xr0[i], xr[i]);
    for (int i = 0; i < xs.size(); i++)
      xs[i] = va::mix(t, xs0[i], xs[i]);
#endif
    _frame++;
  }

  void set_helicity_constraint(bool b) { _helicity_constraint = b; }
  void set_helicity_weight(real w) { _config.w_helicity = w; }
  void set_willmore_weight(real w) { _config.w_willmore = w; }
  void set_area_weight(real w) { _config.w_area = w; }
  void set_shell_strain_weight(real w) { _config.w_shell_strain = w; }
  void set_shell_bending_weight(real w) { _config.w_shell_bending = w; }
  void set_rod_strain_weight(real w) { _config.w_rod_strain = w; }
  void set_rod_bending_weight(real w) { _config.w_rod_bending = w; }

  void set_rod_pin_weight(const real &w) { _config.w_rod_pin = w; }
  void set_rod_weld_weight(real w) { _config.w_rod_weld = w; }
  void set_shell_weld_weight(real w) { _config.w_shell_weld = w; }

  void clear_angle_constraints() { _angle_constraints.clear(); }
  void add_angle_constraint(vec3 axis, real theta, real w) {
    _angle_constraints.push_back(angle_constraint{axis, theta, w});
  }

  void set_repel_rods(bool b) { _repel_rods = b; }
  void set_rod_offset(real o) { _config.rod_offset = o; }
  void set_shell_collisions(bool b) { _shell_collisions = b; }
  void set_pin_rod(bool b) { _pin_rod = b; }

  real get_eps() { return _eps; }

  // std::map<index_t, index_t> _rod_adjacent_edges;
  std::vector<vec3> _fr;
  std::vector<vec3> _fs;

  std::vector<vec3> _pr;

  index_t _frame = 0;
  real _eps = 0.5;

  std::set<index_t> _adjacent;

  asawa::shell::shell::ptr __M;
  asawa::shell::dynamic::ptr __surf;
  asawa::rod::rod::ptr __R;
  asawa::rod::dynamic::ptr __Rd;
  std::vector<real> _willmore_mask;

  struct angle_constraint {
    vec3 axis = vec3(0.0, 0.0, 1.0);
    real theta = 0.1 * M_PI;
    real weight = 1.0;
  };

  std::vector<angle_constraint> _angle_constraints;

  bool _repel_rods = true;
  bool _shell_collisions = true;
  bool _pin_rod = true;
  bool _helicity_constraint = false;
  struct {
    real w_helicity = 1.0e-1;
    real w_willmore = 8e-1;
    real w_area = 2e-2;
    real w_shell_strain = 1.0e-2;
    real w_shell_bending = 1.0e-1;
    real w_rod_strain = 1.0e-1;
    real w_rod_bending = 1.0e-1;
    real w_rod_weld = 1.0;
    real w_shell_weld = 1.0;
    real w_rod_pin = 1.0e-2;
    real rod_offset = 1.0;
  } _config;
};

} // namespace duchamp
} // namespace gaudi

#endif