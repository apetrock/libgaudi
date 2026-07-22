
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

#include "gaudi/calder/least_squares_fit.hpp"
#include "gaudi/calder/rod_integrators.hpp"
#include "gaudi/calder/shell_integrators.hpp"
#include "gaudi/duchamp/dipole_tunneling_constraint.hpp"
#include "gaudi/kusama/cyclide_jet_smooth.hpp"
#include "gaudi/duchamp/fields.hpp"
#include "gaudi/hepworth/block/coupling_collisions_init.hpp"
#include "gaudi/hepworth/block/sim_block.hpp"
#include "gaudi/hepworth/block/solver.hpp"
#include "gaudi/hepworth/block/solver_composition.hpp"
#include "gaudi/hepworth/blocks/rod_position_block.hpp"
#include "gaudi/hepworth/blocks/rod_quaternion_block.hpp"
#include "gaudi/hepworth/blocks/shell_position_block.hpp"
#include "gaudi/hepworth/constraints/bundles.hpp"
#include "gaudi/hepworth/nodes/solver_builder.hpp"

#include "gaudi/common.h"
#include "gaudi/logger.hpp"
#include "module_base_shell.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <vector>
#include "gaudi/geometry_logger.hpp"

namespace gaudi
{

  namespace duchamp
  {

    class knotted_surface_module : public module_base
    {
    public:
      DEFINE_CREATE_FUNC(knotted_surface_module)
      knotted_surface_module(asawa::shell::shell::ptr &M,
                             asawa::shell::dynamic::ptr &D, asawa::rod::rod::ptr &R,
                             asawa::rod::dynamic::ptr &Rd)
          : __M(M), __surf(D), __R(R), __Rd(Rd)
      {

        D->set_flip_pred([&](asawa::shell::shell &M, asawa::shell::CornerId c0)
                         {
      asawa::shell::CornerId c1 = M.other(c0);
      return _adjacent.find(c0) == _adjacent.end() &&
             _adjacent.find(c1) == _adjacent.end(); });

        const std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
        _eps = asawa::shell::avg_length(*__M, x);

        _shell_xs = std::make_shared<shell_vert_positions>(__M, 0);
        _shell_vs =
            std::make_shared<shell_vert_velocities>(__M, __surf->__vdatum_id);
        _shell = std::make_shared<hepworth::block::shell_position_block>(
            __M, _shell_xs, _shell_vs);
        _rod = std::make_shared<hepworth::block::rod_position_block>(__R, __Rd);
        _rod_quat =
            std::make_shared<hepworth::block::rod_quaternion_block>(__R, __Rd);

        _shell->with_force([this]() { return _fs; });
        _rod->with_force([this]() { return _fr; });

        _config_solver =
            hepworth::block::block_solver_builder<
                hepworth::block::shell_position_block,
                hepworth::block::rod_position_block,
                hepworth::block::rod_quaternion_block>::create()
                .with_blocks(_shell, _rod, _rod_quat)
                .with_presolve([this](hepworth::block::solver_context &) {
                  __R->update_lengths();
                  if (_helicity_constraint) {
                    std::vector<real> &lr = __R->l0();
                    for (real &l : lr)
                      l *= 1.02;
                  }
                  // Stash once per step for dipole constraints / forces.
                  refresh_dipole_cache();
                })

                .with_recompute([this](hepworth::block::solver_context &ctx) {
                  const real eps = 0.5 * _eps;
                  const std::vector<vec3> fm;
                  init_rod_shell_weld(*__R, *__Rd, *__M, *__surf, fm,
                                      ctx.constraints, _config.w_rod_weld,
                                      _config.w_shell_weld, 4.0 * eps,
                                      hepworth::block::select_blocks<1, 0>(ctx));
                })
                .with_recompute([this](hepworth::block::solver_context &ctx) {
                  if (_config.w_tunnel_orientation <= 0.0 &&
                      _config.w_dipole_weld <= 0.0 &&
                      _config.w_dipole_weld_rod <= 0.0)
                    return;
                  if (_dipole_Nr.size() != __R->x().size())
                    return;
                  dipole_tunneling::init_dipole_clearance(
                      *__M, *__Rd, ctx.constraints, _shell_xs->get(), __R->x(),
                      _dipole_Nr, _tunnel_r, _config.w_tunnel_orientation,
                      _config.w_dipole_weld, _config.w_dipole_weld_rod,
                      hepworth::block::select_blocks<0, 1>(ctx));
                })
                .with_recompute([this](hepworth::block::solver_context &ctx) {
                  if (_config.w_rod_strain <= 0.0)
                    return;
                  hepworth::block::init_stretch_shear(
                      *__R, ctx.constraints, __R->l0(), _config.w_rod_strain,
                      hepworth::block::select_blocks<1, 2>(ctx));
                })
                .with_recompute([this](hepworth::block::solver_context &ctx) {
                  if (_config.w_rod_straight <= 0.0)
                    return;
                  hepworth::block::init_straight(
                      *__R, ctx.constraints, _config.w_rod_straight,
                      hepworth::block::select_blocks<2>(ctx));
                })
                .with_recompute([this](hepworth::block::solver_context &ctx) {
                  if (_config.w_rod_bending <= 0.0)
                    return;
                  hepworth::block::init_bend_twist(
                      *__R, ctx.constraints, _config.w_rod_bending,
                      hepworth::block::select_blocks<2>(ctx));
                })
                .with_recompute([this](hepworth::block::solver_context &ctx) {
                  if (!_helicity_constraint || _config.w_helicity <= 0.0)
                    return;
                  hepworth::block::init_helicity(
                      *__R, ctx.constraints, _config.w_helicity,
                      hepworth::block::select_blocks<1>(ctx));
                })
                .with_recompute([this](hepworth::block::solver_context &ctx) {
                  auto blocks = hepworth::block::select_blocks<2>(ctx);
                  for (const auto &ac : _angle_constraints) {
                    hepworth::block::init_angle(*__R, ctx.constraints, ac.axis,
                                                ac.theta, ac.weight, blocks);
                  }
                })
                .with_recompute([this](hepworth::block::solver_context &ctx) {
                  if (!_pin_rod)
                    return;
                  auto blocks = hepworth::block::select_blocks<1>(ctx);
                  if (_pr.empty()) {
                    hepworth::block::init_pinned(*__R, ctx.constraints,
                                                 __R->x(), _config.w_rod_pin,
                                                 blocks);
                  } else {
                    hepworth::block::init_pinned(*__R, _pr, ctx.constraints,
                                                 __R->x(), _config.w_rod_pin,
                                                 blocks);
                  }
                })
                .with_recompute([this](hepworth::block::solver_context &ctx) {
                  if (!_repel_rods)
                    return;
                  auto blocks = hepworth::block::select_blocks<1>(ctx);
                  // chain_sep=3: skip AB↔BC and AB↔CD (short segs + growing
                  // rod_offset otherwise inserts cyan mid-edge links along chain)
                  hepworth::block::init_collisions(
                      *__R, *__Rd, ctx.constraints, 1.0, {blocks[0], blocks[0]},
                      _config.rod_offset, /*chain_sep=*/3, /*geom_margin=*/0.1,
                      /*r_override=*/_config.rod_offset * __R->_r);
                })
                .with_recompute([this](hepworth::block::solver_context &ctx) {
                  if (_config.w_shell_strain <= 0.0)
                    return;
                  hepworth::block::init_triangle_strain(
                      *__M, ctx.constraints, _shell_xs->get(),
                      _config.w_shell_strain,
                      hepworth::block::select_blocks<0>(ctx));
                })
                .with_recompute([this](hepworth::block::solver_context &ctx) {
                  if (_config.w_shell_bending <= 0.0)
                    return;
                  hepworth::block::init_bending(
                      *__M, ctx.constraints, _shell_xs->get(),
                      _config.w_shell_bending,
                      hepworth::block::select_blocks<0>(ctx));
                })
                .with_recompute([this](hepworth::block::solver_context &ctx) {
                  if (_config.w_willmore <= 0.0)
                    return;
                  const auto atten = force_attenuation_verts();
                  const auto we =
                      edge_weights_from_vert_attenuation(*__M, atten,
                                                         _config.w_willmore);
                  hepworth::block::init_edge_willmore(
                      *__M, ctx.constraints, we,
                      hepworth::block::select_blocks<0>(ctx));
                })
                .with_recompute([this](hepworth::block::solver_context &ctx) {
                  if (_config.w_area <= 0.0)
                    return;
                  init_weighted_area(*__M, ctx.constraints, _config.w_area,
                                     hepworth::block::select_blocks<0>(ctx));
                })
                .with_recompute([this](hepworth::block::solver_context &ctx) {
                  if (!_shell_collisions)
                    return;
                  const real eps = 0.5 * _eps;
                  auto blocks = hepworth::block::select_blocks<0>(ctx);
                  hepworth::block::init_pnt_tri_collisions(
                      *__M, *__surf, ctx.constraints, _shell_xs->get(),
                      0.5 * eps, 0.5 * eps, 1.0, {blocks[0], blocks[0]});
                })
                .dt(0.05)
                .damping(0.5)
                .iterations(10)
                .build();
      };

      std::vector<vec3> get_rod_normals(asawa::rod::rod &R, asawa::shell::shell &M,
                                        real eps)
      {
        return R.N1();
      }

      std::vector<vec3> get_rod_normals_from_surface(asawa::rod::rod &R, asawa::shell::shell &M,
                                                     real eps)
      {
        std::vector<vec3> &xr = __R->__x;
        std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
        std::vector<vec3> Nf = asawa::shell::face_normals(M, x);
        std::vector<vec3> Nr = calder::mls_avg<vec3>(M, Nf, xr, eps, 2.0);
        for (int i = 0; i < Nr.size(); i++)
        {
          vec3 T = R.dir(asawa::rod::corner_id(i));
          vec3 B = Nr[i].cross(T);
          Nr[i] = T.cross(B);
          Nr[i].normalize();
        }
        return Nr;
      }

      std::vector<real> calc_conservative_collisions(asawa::rod::rod &R,
                                                     asawa::shell::shell &M,
                                                     real eps, int N_spread = 4)
      {
        std::vector<vec3> x_s = asawa::get_vec_data(M, 0);

        // TODO:, these assume that everything is tightly packed, this is wrong
        // assumption if knots get more complicated
        std::vector<vec3> xr = R.x_infill(N_spread);
        std::vector<index_t> rverts(xr.size());
        for (int i = 0; i < rverts.size(); i++)
        {
          rverts[i] = i;
        }
        vector<std::array<index_t, 2>> nearest =
            __surf->get_pnt_tri_collisions(rverts, rverts, xr, M, eps);

        std::vector<real> dist0(xr.size(), 0.0);

        real lavg = R.lavg();
        for (auto &c : nearest)
        {
          auto consec = R.consec(asawa::rod::corner_id(c[0] / N_spread));
          if (consec[2] < 0)
            continue;

          if (c[1] < 0)
          {
            dist0[c[0]] = 4.0 * lavg;
            continue;
          }
#if 0
      index_t ivr = c[0];
      index_t ifs = c[1];

      vec3 xri = xr[ivr];
      vec3 xf = asawa::shell::face_pnt(xri, M, ifs, x_s);
      geometry_logger::line(xri, xf, vec4(1.0, 0.0, 1.0, 1.0));
#endif
        }

        std::vector<real> dist(R.x().size(), 0.0);
        for (int i = 0; i < dist.size(); i++)
        {
          for (int k = 0; k < N_spread; k++)
          {
            dist[i] = std::max(dist0[i * N_spread + k], dist[i]);
          }
        }

        return dist;
      }

      std::vector<real> calc_dist_1(asawa::rod::rod &R, asawa::shell::shell &M,
                                    real eps, int N_spread = 4)
      {

        const std::vector<vec3> &xr = R.x();
        std::vector<real> dist = calc_conservative_collisions(R, M, eps, N_spread);

        auto assign = [](index_t ip, index_t im, const std::vector<vec3> &x,
                         std::vector<real> &dist)
        {
          dist[ip] = std::min(dist[ip], dist[im] + (x[ip] - x[im]).norm());
        };
        // sweep forward

        auto rverts = R.get_ordered_verts();

        for (int k = 0; k < 2; k++)
        {
          for (int i = 0; i < rverts.size(); i++)
          {
            auto ip = R.next(rverts[i]);
            index_t im = i;
            if (ip < 0)
              continue;
            assign(ip, im, xr, dist);
          }

          for (int i = rverts.size() - 1; i > -1; i--)
          {
            auto ip = R.next(rverts[i]);
            index_t im = i;
            if (ip < 0)
              continue;
            assign(im, ip, xr, dist);
          }
        }

        for (int k = 0; k < 8; k++)
        {
          int i0 = k % 2;
          for (int i = i0; i < rverts.size(); i += 2)
          {
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

      geometry_logger::line(xr[i0],
                                xr[i0] + 1.0 * dist[i0] * Nr[i0].normalized(),
                                vec4(0.0, 1.0, 1.0, 1.0));
    }
#endif
        return dist;
      }

      std::vector<real> calc_rod_dist_grad(asawa::rod::rod &R,
                                           asawa::shell::shell &M, real eps,
                                           int N_spread = 4)
      {
        std::vector<real> dist = calc_dist_1(R, M, eps, N_spread);
        std::vector<real> g_d(dist.size(), 0.0);
        real lavg = R.lavg();
        std::vector<vec3> xr = R.x();
        for (int i = 0; i < dist.size(); i++)
        {
          auto idx = R.consec(asawa::rod::corner_id(i));
          auto im = idx[0];
          auto ip = idx[2];
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

      geometry_logger::line(xr[i0],
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
          std::vector<hepworth::sim_block::ptr> blocks)
      {

        const std::vector<vec3> &x0 = R.x();
        std::vector<vec3> x1 = asawa::get_vec_data(M, 0);

        std::vector<index_t> edge_verts_R = R.get_edge_vert_ids();
        std::vector<index_t> edge_map_R = R.get_edge_map();

#if 1
        std::vector<index_t> edge_verts_M = M.get_edge_vert_ids();
        std::vector<index_t> edge_map_M = M.get_edge_map();

        // auto g_d = calc_rod_dist_grad(R, M, 2.0 * eps);

        auto g_d = calc_rod_dist_grad(R, M, 0.5 * eps, 4);
        _willmore_mask = std::vector<real>(x1.size(), 0.0);
        vector<std::array<index_t, 4>> sr_collisions =
            rod_d.get_collisions(edge_verts_M, x1, 2.0 * eps);

        for (auto &c : sr_collisions)
        {
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

          if (abs(g_di) > 1e-3)
          {
            _willmore_mask[vs0] = 1.0;
            _willmore_mask[vs1] = 1.0;
          }

          vec3 xr = va::mix(d[1], xr0, xr1);
          vec3 xs = va::mix(d[2], xs0, xs1);
          vec3 dr = xr1 - xr0;
          vec3 dx = xr - xs;
          vec3 Ns0 =
              asawa::shell::vert_normal(M, asawa::shell::vert_id(vs0), x1);
          vec3 Ns1 =
              asawa::shell::vert_normal(M, asawa::shell::vert_id(vs1), x1);

          vec3 Ns = va::mix(d[2], Ns0, Ns1);

          real is_perp = pow(Ns.dot(dx.normalized()), 2.0);
          // d[0] is squared distance; keep weld while within capture radius.
          const real max_sep = 2.0 * eps;
          if ((is_perp > 0.75 && d[0] < max_sep * max_sep))
          {

            // geometry_logger::line(xs0, xs1, vec4(0.0, 0.0, 1.0, 1.0));

            real lr = (xr1 - xr0).norm();

            // geometry_logger::line(xs, xs + 1.0 * g_di * dr,
            //                           vec4(0.5, 0.5, 1.0, 1.0));
            // geometry_logger::line(xs, xs + 0.1 * g_di * Nri,
            //                          vec4(0.5, 0.5, 1.0, 1.0));
            _adjacent.insert(
                M.find_edge_from_verts(asawa::shell::vert_id(vs0),
                                       asawa::shell::vert_id(vs1))
);
            // geometry_logger::line(xs0, xs1, vec4(0.0, 0.0, 1.0, 1.0));

            hepworth::block::edge_edge_weld::ptr constraint =
                hepworth::block::edge_edge_weld::create(
                    std::vector<index_t>({vr0, vr1, vs0, vs1}), wr, ws, blocks);

            constraints.push_back(constraint);
          }
        }
#endif
      }

      void
      init_torus_flow_constraint_0(
          asawa::rod::rod &R, asawa::rod::dynamic &rod_d, asawa::shell::shell &M,
          asawa::shell::dynamic &shell_d,
          std::vector<hepworth::projection_constraint::ptr> &constraints,
          const real &w, real eps, std::vector<hepworth::sim_block::ptr> blocks)
      {

        const std::vector<vec3> &x0 = R.x();
        std::vector<vec3> x1 = asawa::get_vec_data(M, 0);

        std::vector<index_t> edge_verts_R = R.get_edge_vert_ids();
        std::vector<index_t> edge_map_R = R.get_edge_map();

        std::vector<vec3> Nr = get_rod_normals_from_surface(R, M, eps);
        std::vector<vec3> Nr0 = R.N1();
        std::vector<index_t> edge_verts_M = M.get_edge_vert_ids();
        std::vector<index_t> edge_map_M = M.get_edge_map();

        vector<std::array<index_t, 4>> sr_collisions =
            rod_d.get_collisions(edge_verts_M, x1, 2.0 * eps);

        for (auto &c : sr_collisions)
        {
          if (c[0] < 0)
            continue;
          if (c[1] < 0)
            continue;
          index_t vs0 = c[0];
          index_t vs1 = c[1];
          index_t vr0 = c[2];
          index_t vr1 = c[3];

          asawa::shell::CornerId c0 =
              M.find_edge_from_verts(asawa::shell::vert_id(vs0),
                                     asawa::shell::vert_id(vs1));
          asawa::shell::CornerId c1 = M.other(c0);
          index_t vs00 = M.vert(c0);
          index_t vs10 = M.vert(c1);
          index_t vs01 = M.vert(M.prev(c0));
          index_t vs11 = M.vert(M.prev(c1));

          vec3 xr0 = x0[vr0];
          vec3 xr1 = x0[vr1];
          vec3 xs0 = x1[vs00];
          vec3 xs1 = x1[vs10];
          std::array<real, 3> d = va::distance_Segment_Segment(xr0, xr1, xs0, xs1);

          vec3 xr = va::mix(d[1], xr0, xr1);
          vec3 xs = va::mix(d[2], xs0, xs1);
          vec3 dr = xr1 - xr0;
          vec3 dx = xr - xs;

          vec3 Nri0 = Nr[vr0].normalized();
          vec3 Nri1 = Nr[vr1].normalized();
          vec3 Nri = va::mix(d[1], Nri0, Nri1);
          vec3 Tr = (xr1 - xr0).normalized();
          vec3 Br = Nri.cross(Tr);

          // project xs into null plane of rod, Nri/Br
          vec3 dxsr = xs - xr;
          vec2 xsp = vec2(dxsr.dot(Nri), dxsr.dot(Br));

          real Ndx = Nri.dot(dx);
          vec3 cen = xr - va::sgn(Ndx) * eps * Nri;
          real d_torus = eps - (xs - cen).norm();

          // geometry_logger::line(xs, cen, vec4(1.0, 0.0, 0.5, 1.0));

          if (d_torus > 0.0)
          {

            vec3 Ns = asawa::shell::edge_normal(M, c0, x1);
            real d = va::sgn(Ndx) * d_torus;
            // std::cout << d  << " " << d_torus << std::endl;

            if (d < 0.0)
              geometry_logger::line(xs, xs + 1.0 * d * Ns, vec4(0.0, 1.0, 0.5, 1.0));
            if (d > 0.0)
              geometry_logger::line(xs, xs + 1.0 * d * Ns, vec4(0.5, 0.0, 1.0, 1.0));

            hepworth::block::edge_normal_flow::ptr constraint =
                hepworth::block::edge_normal_flow::create(
                    std::vector<index_t>({vs00, vs01, vs10, vs11}), x1, d, w, blocks);
            constraints.push_back(constraint);
          }
        }
      }

      real min_dist(vec2 pt, const std::vector<vec2> &shape)
      {
        real d = std::numeric_limits<real>::max();
        int N = shape.size();
        for (int i = 0; i < shape.size(); i++)
        {
          vec2 p0 = shape[i];
          vec2 p1 = shape[(i + 1) % N];
          real di = va::dist_from_segment(p0, p1, pt);
          d = std::min(d, di);
        }
        return d;
      }

      // these will tell you if you are inside, then which direction to push
      std::array<real, 2> vec_cross_section(const vec2 &pt, const std::vector<vec2> &shape)
      {
        real wn = va::winding_number(pt, shape);
        real d = min_dist(pt, shape);
        d = wn > 0.0 ? -d : d;
        return {wn, 1.0 * d};
      }

      std::array<real, 2> circ_cross_section(const vec2 &pt, const real &r, const vec2 &c)
      {
        real d = (pt - c).norm() - r;
        real wn = d < 0.0 ? 1.0 : -1.0;
        return {wn, d};
      }

      void
      init_torus_flow_constraint(
          asawa::rod::rod &R, asawa::rod::dynamic &rod_d, asawa::shell::shell &M,
          asawa::shell::dynamic &shell_d,
          std::vector<hepworth::projection_constraint::ptr> &constraints,
          const real &w, real eps, std::vector<hepworth::sim_block::ptr> blocks)
      {

        const std::vector<vec3> &x0 = R.x();
        std::vector<vec3> x1 = asawa::get_vec_data(M, 0);

        std::vector<index_t> edge_verts_R = R.get_edge_vert_ids();
        std::vector<index_t> edge_map_R = R.get_edge_map();

        std::vector<vec3> Nr = get_rod_normals_from_surface(R, M, eps);
        std::vector<vec3> Nr0 = R.N1();
        std::vector<index_t> edge_verts_M = M.get_edge_vert_ids();
        std::vector<index_t> edge_map_M = M.get_edge_map();

        vector<std::array<index_t, 4>> sr_collisions =
            rod_d.get_collisions(edge_verts_M, x1, 2.0 * eps);

        for (auto &c : sr_collisions)
        {
          if (c[0] < 0)
            continue;
          if (c[1] < 0)
            continue;
          index_t vs0 = c[0];
          index_t vs1 = c[1];
          index_t vr0 = c[2];
          index_t vr1 = c[3];

          asawa::shell::CornerId c0 =
              M.find_edge_from_verts(asawa::shell::vert_id(vs0),
                                     asawa::shell::vert_id(vs1));
          asawa::shell::CornerId c1 = M.other(c0);
          index_t vs00 = M.vert(c0);
          index_t vs10 = M.vert(c1);
          index_t vs01 = M.vert(M.prev(c0));
          index_t vs11 = M.vert(M.prev(c1));

          vec3 xr0 = x0[vr0];
          vec3 xr1 = x0[vr1];
          vec3 xs0 = x1[vs00];
          vec3 xs1 = x1[vs10];
          std::array<real, 3> d = va::distance_Segment_Segment(xr0, xr1, xs0, xs1);

          vec3 xr = va::mix(d[1], xr0, xr1);
          vec3 xs = va::mix(d[2], xs0, xs1);
          vec3 dr = xr1 - xr0;

          vec3 Nri0 = Nr[vr0].normalized();
          vec3 Nri1 = Nr[vr1].normalized();
          vec3 Nri = va::mix(d[1], Nri0, Nri1);
          vec3 Tr = (xr1 - xr0).normalized();
          vec3 Br = Nri.cross(Tr);

          // project xs into null plane of rod, Nri/Br
          vec3 dxsr = xs - xr;
          vec2 xsp = vec2(dxsr.dot(Nri), dxsr.dot(Br));

          real Ndx = Nri.dot(xr - xs);
          vec2 cen_U(0.0, eps);
          vec2 cen_L(0.0, -eps);
          vec3 dx = xs - xr;
          vec2 dx2 = vec2(dx.dot(Br), dx.dot(Nri));

          vec2 shape_U[] = {vec2(0, 0), vec2(1, 1), vec2(0, 2), vec2(-1, 1)};
          //vec2 shape_L[] = {vec2(0, 0), vec2(-0.5, -0.2), vec2(-0.5, -1.0), vec2(0.5, -1.0), vec2(0.5, -0.2)};
          //vec2 shape_L[] = {vec2(0, 0), vec2(-0.5, -0.2), vec2(1.0, -2.0), vec2(2.0, -2.0),  vec2(0.5, -0.2)};
          vec2 shape_L[] = {vec2(0.25, 0), vec2(-0.25, 0.0), vec2(-0.125, -1.0), vec2(-0.25, -2.0),  vec2(0.25, -2.0),  vec2(0.125, -1.0)};


          for (vec2 &p : shape_U)
            p *= eps;
          for (vec2 &p : shape_L)
            p *= eps;

          auto [wn_U, d_U] = vec_cross_section(dx2, std::vector<vec2>(shape_U, shape_U + 4));
          // auto [wn_U, d_U] = circ_cross_section(dx2, eps, cen_U);

          auto [wn_L, d_L] = vec_cross_section(dx2, std::vector<vec2>(shape_L, shape_L + 6));
          // auto [wn_L, d_L] = circ_cross_section(dx2, eps, cen_L);

          if (wn_U > 0.0 || wn_L > 0.0)
          {
            real d = wn_U > 0 ? d_U : -d_L;

            // std::cout << d  << " " << d_torus << std::endl;
            vec3 Ns = asawa::shell::edge_normal(M, c0, x1);
            /*
            if (d < 0.0)
              geometry_logger::line(xs, xs + 1.0 * d * Ns, vec4(0.0, 1.0, 0.5, 1.0));
            if (d > 0.0)
              geometry_logger::line(xs, xs + 1.0 * d * Ns, vec4(0.5, 0.0, 1.0, 1.0));
            */
            hepworth::block::edge_normal_flow::ptr constraint =
                hepworth::block::edge_normal_flow::create(
                    std::vector<index_t>({vs00, vs01, vs10, vs11}), x1, d, w, blocks);
            constraints.push_back(constraint);
          }
        }
      }

      std::vector<real> get_dist_rod(const std::vector<index_t> &vert_ids,
                                     const std::vector<vec3> &x,
                                     const asawa::rod::rod &R,
                                     const asawa::rod::dynamic &Rd)
      {
        std::vector<vec3> &xr = __R->x();
        std::vector<real> df(x.size(), 0.0);
        vector<std::array<index_t, 3>> sr_collisions =
            __Rd->get_vert_collisions(vert_ids, x, 99999.9);
        for (auto &c : sr_collisions)
        {
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
      geometry_logger::line(xe[i], xe[i] + 0.1 * df[i] * Ne[i],
                                vec4(0.0, 0.5, 1.0, 1.0));
    }
#endif

    hepworth::block::init_edge_willmore(M, constraints, df, blocks);
  }
#endif
      // Attenuate smoothing where tunnel force is strong:
      // a[v] = 1 - ||f[v]|| / f_max  (1 away from force, 0 at peak).
      std::vector<real> force_attenuation_verts() const {
        std::vector<real> a(_fs.size(), 1.0);
        if (_fs.empty())
          return a;
        real f_max = 0.0;
        for (const vec3 &f : _fs)
          f_max = std::max(f_max, f.norm());
        if (!(f_max > 1e-18))
          return a;
        for (size_t i = 0; i < _fs.size(); ++i)
          a[i] = 1.0 - std::min(1.0, _fs[i].norm() / f_max);
        return a;
      }

      static std::vector<real>
      edge_weights_from_vert_attenuation(const asawa::shell::shell &M,
                                         const std::vector<real> &atten,
                                         real w0) {
        std::vector<real> we(M.edge_count(), w0);
        for (asawa::shell::CornerId c0 : M.get_edge_range()) {
          const index_t i = M.vert(c0);
          const index_t j = M.vert(M.other(c0));
          const real ai =
              (static_cast<size_t>(i) < atten.size()) ? atten[i] : 1.0;
          const real aj =
              (static_cast<size_t>(j) < atten.size()) ? atten[j] : 1.0;
          we[static_cast<index_t>(c0) / 2] = w0 * 0.5 * (ai + aj);
        }
        return we;
      }

      static std::vector<real>
      face_weights_from_vert_attenuation(const asawa::shell::shell &M,
                                         const std::vector<real> &atten,
                                         real w0) {
        std::vector<real> wf;
        wf.reserve(M.face_count());
        for (asawa::shell::FaceId fi : M.get_face_range()) {
          const auto tri = M.get_tri(fi);
          real a = 0.0;
          for (asawa::shell::VertId vi : tri) {
            const index_t i = static_cast<index_t>(vi);
            a += (static_cast<size_t>(i) < atten.size()) ? atten[i] : 1.0;
          }
          wf.push_back(w0 * a / 3.0);
        }
        return wf;
      }

      void init_weighted_area(
          const asawa::shell::shell &M,
          std::vector<hepworth::projection_constraint::ptr> &constraints,
          const real &w, std::vector<hepworth::sim_block::ptr> blocks)
      {

        std::vector<vec3> &xv = asawa::get_vec_data(*__M, 0);
        std::vector<vec3> xf = asawa::shell::face_centers(*__M, xv);
        std::vector<vec3> Nf = asawa::shell::face_normals(*__M, xv);
        std::vector<vec3> &xr = __R->x();
        auto faces_typed = M.get_face_range();
        std::vector<index_t> verts_F(faces_typed.begin(), faces_typed.end());

        real eps = 2.0 * __surf->_Cc;
        std::vector<real> df = get_dist_rod(verts_F, xf, *__R, *__Rd);

#if 0
    for (int i = 0; i < df.size(); i++) {
      geometry_logger::line(xf[verts_F[i]],
                                xf[verts_F[i]] + 0.1 * df[verts_F[i]] * Nf[i],
                                vec4(0.0, 0.5, 1.0, 1.0));
    }
#endif
        std::transform(df.begin(), df.end(), df.begin(),
                       [w, eps](double x)
                       { return max(x - 2.5 * eps, 0.0); });
        auto [min_it, max_it] = std::minmax_element(df.begin(), df.end());
        real vmin = *min_it;
        real vmax = *max_it;
        const real denom = std::max(vmax - vmin, 1e-18);
        std::transform(df.begin(), df.end(), df.begin(),
                       [w, vmin, denom](real x)
                       { return w * (x - vmin) / denom; });

        // Also back off where tunnel force is fighting the area pull.
        const auto atten = force_attenuation_verts();
        const auto face_atten =
            face_weights_from_vert_attenuation(M, atten, /*w0=*/1.0);
        const size_t n = std::min(df.size(), face_atten.size());
        for (size_t i = 0; i < n; ++i)
          df[i] *= face_atten[i];

        hepworth::block::init_area(M, constraints, xv, df, blocks,
                                     hepworth::block::area_mode::zero);
        // return g;
      }

      void assert_nan(index_t k)
      {
        (void)k;
        for (int i = 0; i < __R->__u.size(); ++i)
        {
          if (__R->__u[i].coeffs().hasNaN())
          {
            std::cerr << "nan at " << i << std::endl;
            exit(0);
          }
        }
      }

      void init_step(real h)
      {
        std::vector<vec3> &xs = asawa::get_vec_data(*__M, 0);
        std::vector<vec3> &xr = __R->x();
        _pr = std::vector<vec3>(0, vec3::Zero());
        _fs = std::vector<vec3>(xs.size(), vec3::Zero());
        _fr = std::vector<vec3>(xr.size(), vec3::Zero());
      }

      void set_rod_pin(const std::vector<vec3> &x, const real &w = 1e-1)
      {
        _pin_rod = true;
        _pr = x;
        _config.w_rod_pin = w;
      }

      void add_rod_force(const std::vector<vec3> &fr, const real &h = 1.0)
      {
        for (int i = 0; i < fr.size(); i++)
        {
          _fr[i] += h * fr[i];
        }
      }

      void add_shell_force(const std::vector<vec3> &fs, const real &h = 1.0)
      {
        for (int i = 0; i < fs.size(); i++)
        {
          _fs[i] += h * fs[i];
        }
      }

      void step(real h)
      {
        _config_solver.dt = h;
        _config_solver.damping = 0.5;
        _config_solver.iterations = 10;

        std::vector<vec3> &xs = asawa::get_vec_data(*__M, 0);
        std::vector<vec3> &xr = __R->x();
        const std::vector<vec3> xs0 = xs;
        const std::vector<vec3> xr0 = xr;

        // Shell forces: compute Nr once, then each force gates on its own weight.
        if (_config.w_tunnel_force > 0.0 || _config.w_darboux_force > 0.0) {
          refresh_dipole_cache();

          // Tunnel: f = (w/h²)*dx along mesh N onto dipole cylinder.
          if (_config.w_tunnel_force > 0.0 && _dipole_Nr.size() == xr.size()) {
            std::vector<vec3> tunnel;
            dipole_tunneling::accumulate_dipole_tunnel_forces(
                *__M, xs, xr, _dipole_Nr, __R->get_edge_vert_ids(), _tunnel_r,
                _tunnel_r, /*force_w unused=*/1.0, tunnel);
            if (tunnel.size() == _fs.size()) {
              const real inv_h2 = _config.w_tunnel_force / (h * h);
              for (size_t i = 0; i < _fs.size(); ++i)
                _fs[i] += inv_h2 * tunnel[i];
            }
          }

          // Rod LS-Darboux: first-order SDF pull onto D=0 along mesh N.
          // Homogeneous fit → ||Q||=1, so raw D is gauge junk; use D/|∇D|.
          // (Eigenvalue of A measures fit quality, not geometric scale.)
          if (_config.w_darboux_force > 0.0 && _dipole_Nr.size() == xr.size() &&
              _fs.size() == xs.size()) {
            const std::vector<vec3> Ns =
                asawa::shell::vertex_normals(*__M, xs);
            const real l0 = _config.darboux_fit_scale * _eps;
            auto Q = calder::darboux_cyclide_tangent_plane(
                *__R, _dipole_Nr, xs, Ns, l0, _config.darboux_fit_p0,
                _config.darboux_foot_normal_w, _config.darboux_fit_p1);
            // Same mesh jet-smooth as medial-axis: kills BH / eigen C0 noise.
            if (_config.darboux_smooth)
              Q = kusama::cyclide_jet_smooth(*__M, xs, Q,
                                             _config.darboux_smooth_params);
            const real inv_h2 = _config.w_darboux_force / (h * h);
            const real max_travel = 4.0 * _eps;
            for (size_t i = 0; i < xs.size(); ++i) {
              if (i >= Q.size())
                continue;
              const real D0 = albers::eval_darboux(Q[i], vec3::Zero());
              const vec3 g0 = albers::darboux_grad(Q[i], vec3::Zero());
              const vec3 &Ni = Ns[i];
              const real g_n = g0.norm();
              if (!std::isfinite(D0) || !g0.allFinite() || !(g_n > 1e-12) ||
                  !Ni.allFinite() || !(Ni.squaredNorm() > 1e-24))
                continue;

              // Signed distance ≈ D/|∇D|; push along mesh N toward D=0.
              const real dist = D0 / g_n;
              vec3 dx = -dist * Ni.normalized();
              if (!dx.allFinite())
                continue;
              const real dn = dx.norm();
              if (dn > max_travel)
                dx *= max_travel / dn;

              if (dx.squaredNorm() > 1e-18)
                geometry_logger::line(xs[i], xs[i] + dx,
                                      vec4(0.2, 0.9, 0.4, 1.0));
              _fs[i] = inv_h2 * dx;
            }
          }
        }

        hepworth::block::run_solver_step(_config_solver, _solver);

        const real t = 0.5;
        for (size_t i = 0; i < xr.size(); i++)
          xr[i] = va::mix(t, xr0[i], xr[i]);
        for (size_t i = 0; i < xs.size(); i++)
          xs[i] = va::mix(t, xs0[i], xs[i]);
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
      void set_rod_straight_weight(real w) { _config.w_rod_straight = w; }

      void set_rod_pin_weight(const real &w) { _config.w_rod_pin = w; }
      void set_rod_weld_weight(real w) { _config.w_rod_weld = w; }
      void set_shell_weld_weight(real w) { _config.w_shell_weld = w; }
      void set_tunnel_orientation_weight(real w) {
        _config.w_tunnel_orientation = w;
      }
      void set_dipole_weld_weight(real w) { _config.w_dipole_weld = w; }
      void set_dipole_weld_rod_weight(real w) { _config.w_dipole_weld_rod = w; }
      void set_tunnel_force_weight(real w) { _config.w_tunnel_force = w; }
      void set_darboux_force_weight(real w) { _config.w_darboux_force = w; }
      void set_darboux_smooth(bool b) { _config.darboux_smooth = b; }

      void clear_angle_constraints() { _angle_constraints.clear(); }
      void add_angle_constraint(vec3 axis, real theta, real w)
      {
        _angle_constraints.push_back(angle_constraint{axis, theta, w});
      }

      void set_repel_rods(bool b) { _repel_rods = b; }
      void set_rod_offset(real o) { _config.rod_offset = o; }
      void set_dipole_radius(real r) { _config.dipole_radius = r; }
      void set_shell_collisions(bool b) { _shell_collisions = b; }
      void set_pin_rod(bool b) { _pin_rod = b; }

      real get_eps() { return _eps; }

      // Dipole tunnel / weld / force cylinder radius.
      // If dipole_radius > 0: use it as an absolute radius (independent of
      // rod_offset / __R->_r). Else legacy: rod_offset * __R->_r.
      // Always capped by a few mesh edge lengths.
      real tunnel_r() const
      {
        const real r = (_config.dipole_radius > 0.0)
                           ? _config.dipole_radius
                           : _config.rod_offset * __R->_r;
        return std::min(r, 16.0 * _eps);
      }

      void refresh_dipole_cache()
      {
        _tunnel_r = tunnel_r();
        _dipole_Nr = get_rod_normals_from_surface(*__R, *__M, 4.0 * _eps);
      }

      // std::map<index_t, index_t> _rod_adjacent_edges;
      std::vector<vec3> _fr;
      std::vector<vec3> _fs;
      std::vector<vec3> _dipole_Nr;
      real _tunnel_r = 0.0;

      std::vector<vec3> _pr;

      index_t _frame = 0;
      real _eps = 0.5;

      std::set<index_t> _adjacent;

      asawa::shell::shell::ptr __M;
      asawa::shell::dynamic::ptr __surf;
      asawa::rod::rod::ptr __R;
      asawa::rod::dynamic::ptr __Rd;
      std::vector<real> _willmore_mask;

      std::shared_ptr<shell_vert_positions> _shell_xs;
      std::shared_ptr<shell_vert_velocities> _shell_vs;
      hepworth::block::shell_position_block::ptr _shell;
      hepworth::block::rod_position_block::ptr _rod;
      hepworth::block::rod_quaternion_block::ptr _rod_quat;
      hepworth::block::block_solver_config<hepworth::block::shell_position_block,
                                           hepworth::block::rod_position_block,
                                           hepworth::block::rod_quaternion_block>
          _config_solver;
      hepworth::block::projection_solver _solver;

      struct angle_constraint
      {
        vec3 axis = vec3(0.0, 0.0, 1.0);
        real theta = 0.1 * M_PI;
        real weight = 1.0;
      };

      std::vector<angle_constraint> _angle_constraints;

      bool _repel_rods = true;
      bool _shell_collisions = true;
      bool _pin_rod = true;
      bool _helicity_constraint = false;
      struct
      {
        real w_helicity = 1.0e-1;
        real w_willmore = 5e-1;
        real w_area = 1e-2;
        real w_shell_strain = 1.0e-2;
        real w_shell_bending = 2.0e-1;
        real w_rod_straight = 1.0e-2;
        real w_rod_strain = 1.0e-1;
        real w_rod_bending = 4.0e-1;
        real w_rod_weld = 1.0;
        real w_shell_weld = 1.0;
        real w_tunnel_orientation = 0.01; // triangle_dipole_tunneling (face N)
        real w_dipole_weld = 0.1;       // dipole_weld shell side
        real w_dipole_weld_rod = 0.01;  // dipole_weld rod side
        real w_tunnel_force = 0.1;       // f = (w/h²)*dx; mesh-N ray → dipole cyl
        real w_darboux_force = 0.0;      // f = (w/h²)*(d0-d); rod LS-Darboux zero set
        real darboux_fit_scale = 0.1;    // fit length = scale * _eps
        real darboux_fit_p0 = 3.0;       // κ_inv_dist power
        real darboux_fit_p1 = 8.0;       // sin(φ) radial-gate power
        // Soft tip-in of mesh N at query: w_foot = foot_normal_w * Σ w_MLS.
        // Encourages ∇D ∥ N (and D≈0) at the vert — pipe tangent to the shell.
        real darboux_foot_normal_w = 10.0;
        bool darboux_smooth = true;      // kusama jet-smooth Q on mesh (as medial)
        kusama::cyclide_jet_smooth_params darboux_smooth_params{
            .wi = 0.5,               // weaker anchor → more neighbor agreement
            .alpha_G = 2.0,          // gradient (C1) match across edges
            .alpha_H = 0.1,          // Hessian match
            .use_cotan_weights = true,
            .sweeps = 2,             // Jacobi iterations
        };
        real w_rod_pin = 1.0e-2;
        real rod_offset = 1.0;           // rod–rod collision scale on __R->_r
        // Absolute dipole cylinder radius for weld/tunnel. <=0 → legacy
        // rod_offset * __R->_r.
        real dipole_radius = 0.0;
      } _config;
    };

  } // namespace duchamp
} // namespace gaudi

#endif
