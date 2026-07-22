#pragma once

#include <cmath>
#include <memory>
#include <set>
#include <vector>

#include "gaudi/asawa/rod/dynamic.hpp"
#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/asawa/shell/dynamic.hpp"
#include "gaudi/duchamp/dipole_tunneling_constraint.hpp"
#include "gaudi/duchamp/dipole_tunneling_geometry.hpp"
#include "gaudi/duchamp/fields.hpp"
#include "gaudi/geometry_logger.hpp"
#include "gaudi/hepworth/blocks/rod_position_block.hpp"
#include "gaudi/hepworth/blocks/shell_position_block.hpp"
#include "gaudi/hepworth/constraints/bundles.hpp"
#include "gaudi/hepworth/nodes/solver_builder.hpp"

namespace gaudi {
namespace duchamp {

// Visual harness: unit sphere + elastic stretch/bend, with a dipole rod
// (3D line + constant normal) for disk / tunnel constraint development.
class dipole_tunneling_demo {
public:
  using ptr = std::shared_ptr<dipole_tunneling_demo>;

  static ptr create() { return std::make_shared<dipole_tunneling_demo>(); }

  dipole_tunneling_demo() { reset_scene(); }

  void reset_scene() {
    __M = asawa::shell::load_sphere(1.0, 32, 20);
    asawa::init_vert_datum<vec3>(*__M, vec3::Zero());

    auto xs = std::make_shared<shell_vert_positions>(__M, 0);
    auto vs = std::make_shared<shell_vert_velocities>(__M, 1);
    _shell = std::make_shared<hepworth::block::shell_position_block>(__M, xs, vs);
    _shell->with_force([this]() {
      if (!_use_tunnel_force)
        return std::vector<vec3>(__M->vert_count(), vec3::Zero());
      return calc_tunnel_radial_forces();
    });

    const real l0 = 1.0 * asawa::shell::avg_length(*__M, xs->get());
    __surf = asawa::shell::dynamic::create(__M,  l0, 1.25 * l0, 0.5 * l0);

    // Line through the sphere; material-frame N1() drives dipole Nr (slerp along
    // segments). _dipole_N is debug-only broadcast reference.
    _dipole_N = vec3(0.18, 0.88, 0.42).normalized();
    const std::vector<vec3> rod_pts = {
        vec3(-1.35, 0.12, 0.08),
        vec3(-0.45, 0.04, -0.02),
        vec3(0.35, -0.06, 0.05),
        vec3(1.25, -0.10, 0.12),
    };
    __R = asawa::rod::rod::create(rod_pts, false);
    const real rod_lavg = __R->lavg();
    __Rd = asawa::rod::dynamic::create(__R, 0.25 * rod_lavg, 2.5 * rod_lavg,
                                       0.25 * rod_lavg);
    _rod = std::make_shared<hepworth::block::rod_position_block>(__R, __Rd);

    _disk_r = 0.16;
    _disk_h = _disk_r;
    _stretch_w = 1.0;
    _bend_w = 0.1;
    _laplacian_w = 0.1;
    _area_w = 0.1;
    _rod_pin_w = 1.0;
    _tunnel_w = 1.5;        // triangle_dipole_tunneling
    _weld_w = 1.5;          // dipole_weld shell side
    _weld_rod_w = 0.01;     // dipole_weld rod side (independent)
    _use_tunnel_force = false;
    _tunnel_force_w = 0.1;
    _solver_dt = 0.05;
    _solver_iters = 16;

    auto Nr_fn = [this]() { return __R->N1(); };

    // Pins allocate/anchor rod DOFs in the joint matrix so coupling can hang
    // off the rod block and update rod positions during relaxation.
    _config =
        hepworth::block::block_solver_builder<
            hepworth::block::shell_position_block,
            hepworth::block::rod_position_block>::create()
            .with_blocks(_shell, _rod)
            .with_bundle(hepworth::block::make_shell_physics_bundle<0>(
                _shell, _stretch_w, _bend_w))
            .with_bundle(hepworth::block::make_shell_laplacian_bundle<0>(
                _shell, _laplacian_w, hepworth::block::laplacian_mode::null,
                hepworth::block::laplacian_stencil::cotan))
            .with_bundle(hepworth::block::make_shell_area_bundle<0>(
                _shell, _area_w, hepworth::block::area_mode::zero))
            .with_bundle(hepworth::block::make_rod_pin_bundle<1>(_rod, _rod_pin_w))
            .with_recompute(duchamp::make_dipole_clearance_recompute<0, 1>(
                _shell, _rod, Nr_fn, [this]() { return _disk_r; },
                [this]() { return _tunnel_w; }, [this]() { return _weld_w; },
                [this]() { return _weld_rod_w; }))
            .dt(_solver_dt)
            .damping(0.2)
            .iterations(_solver_iters)
            .build();
  }

  void step(int frame) {
    (void)frame;
    run_elastic_step();
    draw_dipole_debug();
    gather_and_draw_captured_faces();
    draw_tunnel_forces();
    __surf->step(true);
  }

  asawa::shell::shell::ptr __M;
  asawa::shell::dynamic::ptr __surf;
  asawa::rod::rod::ptr __R;
  asawa::rod::dynamic::ptr __Rd;

  vec3 dipole_normal() const { return _dipole_N; }
  real disk_offset() const { return _disk_h; }
  real disk_radius() const { return _disk_r; }

private:
  std::vector<vec3> calc_tunnel_radial_forces() const {
    std::vector<vec3> disp(__M->vert_count(), vec3::Zero());
    dipole_tunneling::accumulate_dipole_tunnel_forces(
        *__M, _shell->xs->get(), __R->x(), __R->N1(), __R->get_edge_vert_ids(),
        _disk_r, _disk_h, /*force_w unused=*/1.0, disp);
    const real inv_h2 = _tunnel_force_w / (_solver_dt * _solver_dt);
    for (vec3 &f : disp)
      f *= inv_h2;
    return disp;
  }

  void run_elastic_step() {
    hepworth::block::run_solver_step(_config, _solver);
  }

  void draw_disk_ring(const vec3 &cen, const vec3 &B, const vec3 &N, real radius,
                      const vec4 &color, int segments = 32) const {
    for (int i = 0; i < segments; ++i) {
      const real a0 = 2.0 * M_PI * real(i) / real(segments);
      const real a1 = 2.0 * M_PI * real(i + 1) / real(segments);
      const vec3 p0 = cen + radius * (std::cos(a0) * B + std::sin(a0) * N);
      const vec3 p1 = cen + radius * (std::cos(a1) * B + std::sin(a1) * N);
      geometry_logger::line(p0, p1, color, 0.005);
    }
  }

  void draw_face_edges(asawa::shell::FaceId fi, const vec4 &color) const {
    const std::vector<vec3> &x = _shell->xs->get();
    std::vector<vec3> pts;
    __M->const_for_each_face(fi, [&](asawa::shell::CornerId c0, const asawa::shell::shell &) {
      pts.push_back(x[__M->vert(c0)]);
    });
    if (pts.size() < 2)
      return;
    for (size_t i = 0; i < pts.size(); ++i) {
      geometry_logger::line(pts[i], pts[(i + 1) % pts.size()], color, 0.003);
    }
  }

  void draw_tunnel_forces() const {
    if (!_use_tunnel_force)
      return;
    const std::vector<vec3> &x = _shell->xs->get();
    const std::vector<vec3> forces = calc_tunnel_radial_forces();
    const real scale = 0.002;
    const vec4 upper_col(1.0, 0.55, 0.2, 1.0);
    const vec4 lower_col(0.35, 0.65, 1.0, 1.0);

    for (size_t vi = 0; vi < forces.size(); ++vi) {
      const vec3 &f = forces[vi];
      if (f.squaredNorm() <= 1e-18)
        continue;
      const vec4 &col = f.dot(_dipole_N) >= 0.0 ? upper_col : lower_col;
      geometry_logger::line(x[vi], x[vi] + scale * f, col, 0.002);
    }
  }

  void draw_dipole_debug() {
    const std::vector<vec3> &rod_pts = __R->x();
    const std::vector<vec3> Nr = __R->N1();
    const vec4 upper_col(1.0, 0.0, 0.0, 1.0);
    const vec4 lower_col(0.0, 0.2, 1.0, 1.0);
    const vec4 upper_box_col(1.0, 0.35, 0.35, 1.0);
    const vec4 lower_box_col(0.35, 0.45, 1.0, 1.0);

    const std::vector<index_t> debug_edges = __R->get_edge_vert_ids();
    for (size_t e = 0; e + 1 < debug_edges.size(); e += 2) {
      const index_t s0 = debug_edges[e];
      const index_t s1 = debug_edges[e + 1];
      dipole_tunneling::segment_frame f;
      if (!dipole_tunneling::try_make_segment_frame(
              rod_pts[s0], rod_pts[s1],
              dipole_tunneling::interpolate_rod_normal(Nr[s0], Nr[s1], 0.5), f))
        continue;
      const auto &[xr0, xr1, N, T, B] = f;
      (void)xr0;
      (void)xr1;
      (void)N;
      (void)T;
      (void)B;

      const auto &[upper, lower] =
          dipole_tunneling::make_dual_tunnel_aabbs(f, _disk_r);
      geometry_logger::ext(upper[0], upper[1], upper_box_col);
      geometry_logger::ext(lower[0], lower[1], lower_box_col);

      const int samples = 5;
      // Spikes stick out of each tube circle along ±N (N lies in the disk
      // plane, so rod→disk lines read as diameters and get lost in the rings).
      const real spike = 1.5 * _disk_r;
      const real thick = 0.012;
      for (int k = 0; k <= samples; ++k) {
        const real t = real(k) / real(samples);
        const vec3 xr = va::mix(t, rod_pts[s0], rod_pts[s1]);
        const vec3 Ni =
            dipole_tunneling::interpolate_rod_normal(Nr[s0], Nr[s1], t);
        dipole_tunneling::segment_frame ft;
        if (!dipole_tunneling::try_make_segment_frame(rod_pts[s0], rod_pts[s1],
                                                     Ni, ft))
          continue;
        const auto &[fxr0, fxr1, fN, fT, fB] = ft;
        (void)fxr0;
        (void)fxr1;
        (void)fT;
        const vec3 cen_u = xr + _disk_h * fN;
        const vec3 cen_l = xr - _disk_h * fN;
        draw_disk_ring(cen_u, fB, fN, _disk_r, upper_col);
        draw_disk_ring(cen_l, fB, fN, _disk_r, lower_col);
        geometry_logger::line(cen_u, cen_u + spike * fN, upper_col, thick);
        geometry_logger::line(cen_l, cen_l - spike * fN, lower_col, thick);
      }
    }
  }



  void gather_and_draw_captured_faces() {
    const std::vector<vec3> &rod_pts = __R->x();
    const std::vector<vec3> Nr = __R->N1();
    const std::vector<vec3> &x = _shell->xs->get();
    const std::vector<vec3> face_centers = asawa::shell::face_centers(*__M, x);

    const vec4 circle_col(0.85, 0.35, 0.35, 1.0);
    const vec4 ring_col(0.85, 0.85, 0.85, 1.0);
    const vec4 upper_col(1.0, 0.85, 0.35, 1.0);
    const vec4 lower_col(0.35, 0.65, 1.0, 1.0);

    const std::vector<index_t> edges = __R->get_edge_vert_ids();
    std::vector<std::tuple<asawa::shell::FaceId, index_t, index_t, bool>>
        total_captured;
    for (size_t e = 0; e + 1 < edges.size(); e += 2) {
      const index_t v0 = edges[e];
      const index_t v1 = edges[e + 1];
      dipole_tunneling::segment_frame f;
      if (!dipole_tunneling::try_make_segment_frame(
              rod_pts[v0], rod_pts[v1],
              dipole_tunneling::interpolate_rod_normal(Nr[v0], Nr[v1], 0.5), f))
        continue;
      const auto [upper, lower] =
          dipole_tunneling::gather_faces_in_segment(*__M, x, f, _disk_r);
      for (const asawa::shell::FaceId fid : upper)
        total_captured.emplace_back(fid, v0, v1, true);
      for (const asawa::shell::FaceId fid : lower)
        total_captured.emplace_back(fid, v0, v1, false);
    }

    std::set<int> drawn_ring_faces;
    for (const auto &[face_id, v0, v1, above] : total_captured) {
      dipole_tunneling::segment_frame f;
      if (!dipole_tunneling::try_make_segment_frame(
              rod_pts[v0], rod_pts[v1],
              dipole_tunneling::interpolate_rod_normal(Nr[v0], Nr[v1], 0.5), f))
        continue;
      const auto &[xr0, xr1, N, T, B] = f;
      (void)T;
      const vec3 xf = face_centers[static_cast<size_t>(face_id)];
      const vec3 xr = dipole_tunneling::safe_project_on_line(xr0, xr1, xf);
      const vec3 dp = xf - xr;
      const real R = dipole_tunneling::safe_tangent_point_radius(dp, N);
      if (!std::isfinite(R))
        continue;
      const real sg = va::sgn(N.dot(dp));
      const vec3 cen = xr + ((sg == 0.0) ? 1.0 : sg) * R * N;
      draw_disk_ring(cen, B, N, R, circle_col);
      draw_face_edges(face_id, above ? upper_col : lower_col);
      //geometry_logger::line(xr, xf, vec4(0.5, 0.5, 0.5, 1.0));

      const std::vector<asawa::shell::FaceId> ring =
          __M->face_one_ring_face_ids(face_id);
      for (const asawa::shell::FaceId &rf : ring) {
        const int rid = static_cast<int>(rf);
        if (drawn_ring_faces.insert(rid).second)
          draw_face_edges(rf, ring_col);
      }
    }
  }

  hepworth::block::shell_position_block::ptr _shell;
  hepworth::block::rod_position_block::ptr _rod;
  hepworth::block::block_solver_config<hepworth::block::shell_position_block,
                                       hepworth::block::rod_position_block>
      _config;
  hepworth::block::projection_solver _solver;
  vec3 _dipole_N = vec3::UnitY();
  real _disk_h = 0.16;
  real _disk_r = 0.16;
  real _stretch_w = 1.0;
  real _bend_w = 0.12;
  real _laplacian_w = 0.1;
  real _area_w = 0.0;
  real _rod_pin_w = 0.0;
  real _tunnel_w = 0.0;
  real _weld_w = 0.0;
  real _weld_rod_w = 0.0;
  bool _use_tunnel_force = false;
  real _tunnel_force_w = 0.0;
  real _solver_dt = 0.05;
  int _solver_iters = 6;
};

} // namespace duchamp
} // namespace gaudi
