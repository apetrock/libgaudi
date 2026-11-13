#include <algorithm>
#include <cstdio>
#include <array>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <string>

#include "gaudi/duchamp/monge_fairing_demo.hpp"
#include "gaudi/kusama/principal_curvature_solver.hpp"

// Runs the monge fairing demo headless with the curvature solver and prints,
// per frame: mean / max |dx|, the fitted torus ring and tube radii, the
// principal-curvature roughness (mean |H_i - H_j| over edges), the angular
// momentum of the vertex velocities, the vertex count, and wall time.
//
//   principal_curvature_probe [frames] [label] [w] [pd iterations]
//                   [remesh] [nodrive] [refresh=K] [backtracks=K] [areafree] [areaband=lo:hi] [warea=W]
//                   [solver=contact_net|off] [cntan=vertex|edge] [cnframes] [nmass=C]
//                   [wdcircle=W] [dcmode=closest|slide|even|graded|hyperbolic|parabolic|loxodromic] [dcfam=0|1|2]
//                   [dcamt=A] [dcangle=DEG] [dcref=x,y,z] [wamap=W] [amapinv=A] [wconf=W] [wtension=W]
//                   [dupin=W] [dupinreach=R] [dupindeg=D] [wgradtl=W] [wotorus=W] [sarap=W] [sarapfree]
//                   [wkite=W] [kalpha=A] [kangle=DEG] [kmaxu=U] [ksrc=tangent|mean|aniso] [klevel=L] [kbalanced]
//                   [wcmc=W] [cmch0=H] [cmcvol] [wtwist=W] [tau0=T] [twangle=A] [wmin=W] [growth=x[:y]] [gpower=P] [wgrowth=W] [wdet=W]
//                   [accel=anderson|chebyshev|anderson_chebyshev|none] [freshstop=0|1] [maxascent=K]
//                   [dump=K] [res=ring:tube] [perturb=P] [wcenter=W] [recenter]
//                   [wstretch=W] [wind=ring:tube[:hand]] [sdir=angle|material|stripes] [msmooth=L] [mfollow=F]
//                   [stripes=P] [tilt=G] [wtilt=W] [tiltaxis=k0|k1] [tiltmode=rigid|spin|twist|bend]
//                   [tiltsign=gradient|propagate] [aratio=W] [bend=W] [damp=D] [memory=W]
//                   [qcwarn=Q] [noqc] [ascentrefresh] [ascenttol=T] [audit] [modes=K] [tau=T] [lmax=L]
//
// `tilt=` sets the position tilt gain (`w_tilt` 1 unless `wtilt=`). Each
// frame also prints the mean log ratio and anisotropy q of the live edge
// readings, elliptic and saddle edges apart, and the mean turn of their
// principal axes since the previous frame. `label` is printed only; `w` is
// the coupling weight.

namespace gaudi {

/// Standalone solver timing on the 48 x 24 torus with the x1.2 targets, for
/// each `analytic` / `prefactor` pair: mean per-phase seconds per iterate and
/// iterations to converge.
int principal_curvature_bench() {
  std::vector<vec3> x0;
  std::vector<std::vector<int>> faces;
  asawa::make_torus(x0, faces, 1.0, 0.33, 48, 24);
  std::vector<kusama::tri3> tris;
  for (const auto &f : faces)
    tris.push_back({f[0], f[1], f[2]});
  std::vector<vec2> kh(x0.size());
  for (size_t i = 0; i < x0.size(); ++i) {
    const real phi = 2.0 * M_PI * real(static_cast<int>(i) % 24) / 24.0;
    kh[i] = vec2(1.0 / 0.396, std::cos(phi) / (1.2 + 0.396 * std::cos(phi)));
  }
  for (const bool analytic : {false, true})
    for (const bool prefactor : {false, true}) {
      kusama::principal_curvature_solver s;
      s.analytic = analytic;
      s.prefactor = prefactor;
      s.reset(static_cast<int>(x0.size()), tris, x0);
      std::vector<vec3> x = x0;
      kusama::principal_curvature_timing sum;
      int its = 0, stalls = 0;
      const real e0 = s.energy(x, kh);
      real e = e0;
      for (; its < 400 && stalls < 4; ++its) {
        const kusama::principal_curvature_report rep = s.iterate(x, kh);
        sum.jacobian += s.last_timing.jacobian;
        sum.factor += s.last_timing.factor;
        sum.solve += s.last_timing.solve;
        sum.energy += s.last_timing.energy;
        sum.total += s.last_timing.total;
        if (!rep.accepted) {
          ++stalls;
          continue;
        }
        stalls = 0;
        e = rep.energy1;
        if (rep.energy0 - rep.energy1 <= 1e-12 * rep.energy0)
          break;
      }
      const double n = std::max(its, 1);
      std::cout << "analytic " << analytic << " prefactor " << prefactor
                << "  its " << its << "  energy " << e0 << " -> " << e
                << "  per-iterate ms: jacobian " << 1e3 * sum.jacobian / n
                << " factor " << 1e3 * sum.factor / n << " solve "
                << 1e3 * sum.solve / n << " energy " << 1e3 * sum.energy / n
                << " total " << 1e3 * sum.total / n << "  wall ms "
                << 1e3 * sum.total << "\n";
    }
  return 0;
}

int principal_curvature_probe_main(int argc, char **argv) {
  if (argc > 1 && std::string(argv[1]) == "bench")
    return principal_curvature_bench();
  const int frames = argc > 1 ? std::atoi(argv[1]) : 10;
  const std::string target = argc > 2 ? argv[2] : "delta";
  gaudi::duchamp::monge_fairing_demo_config cfg{};
  auto &eig = cfg.project.curvature;
  if (argc > 3)
    eig.w_coupling = std::atof(argv[3]);
  if (argc > 4)
    cfg.step.iterations = std::atoi(argv[4]);
  bool remesh = false;
  int dump_every = 0;
  for (int a = 5; a < argc; ++a) {
    const std::string flag = argv[a];
    if (flag == "remesh") {
      remesh = true;
    } else if (flag == "nodrive") {
      cfg.shape.w_growth = 0;
      cfg.shape.w_det = 0;
      cfg.shape.w_dupin = 0;
      cfg.shape.w_grad_traceless = 0;
      cfg.shape.w_shape_arap = 0;
      std::cout << "shape drives off\n";
    } else if (flag.rfind("refresh=", 0) == 0) {
      cfg.step.gn_refresh_every = std::atoi(flag.c_str() + 8);
    } else if (flag.rfind("areaband=", 0) == 0) {
      // lo:hi on the triangle area rows.
      const std::string v = flag.substr(9);
      const size_t c = v.find(':');
      auto &band = cfg.position.area.band;
      band.lo = std::atof(v.substr(0, c).c_str());
      band.hi = c == std::string::npos ? band.lo : std::atof(v.substr(c + 1).c_str());
    } else if (flag.rfind("warea=", 0) == 0) {
      cfg.position.area.w = std::atof(flag.c_str() + 6);
    } else if (flag.rfind("growth=", 0) == 0) {
      // x[:y] along the tight direction (or the first entry's angle).
      const std::string v = flag.substr(7);
      const size_t c = v.find(':');
      hepworth::block::curvature_growth_target g =
          cfg.shape.growth.empty() ? hepworth::block::curvature_growth_target{} : cfg.shape.growth[0];
      g.x = std::atof(v.substr(0, c).c_str());
      g.y = c == std::string::npos ? real(0) : std::atof(v.substr(c + 1).c_str());
      cfg.shape.growth = {g};
      if (!(cfg.shape.w_growth > 0))
        cfg.shape.w_growth = 1;
    } else if (flag.rfind("gpower=", 0) == 0) {
      // Curvature feedback on the growth entry (`curvature_growth_target::power`).
      if (cfg.shape.growth.empty())
        cfg.shape.growth = {hepworth::block::curvature_growth_target{}};
      for (auto &g : cfg.shape.growth)
        g.power = std::atof(flag.c_str() + 7);
    } else if (flag.rfind("wgrowth=", 0) == 0) {
      cfg.shape.w_growth = std::atof(flag.c_str() + 8);
    } else if (flag.rfind("wdet=", 0) == 0) {
      cfg.shape.w_det = std::atof(flag.c_str() + 5);
    } else if (flag == "areafree") {
      // The axis-ratio rows' area free, as before the band.
      for (auto &e : cfg.position.axis_ratio)
        e.area = hepworth::block::area_band::free();
    } else if (flag.rfind("tiltsign=", 0) == 0) {
      cfg.position.tilt.at(0).sign = std::string(flag.c_str() + 9) == "propagate"
                                         ? hepworth::block::tilt_sign::propagate
                                         : hepworth::block::tilt_sign::gradient;
    } else if (flag.rfind("tilt=", 0) == 0) {
      cfg.position.tilt.at(0).gain = std::atof(flag.c_str() + 5);
      if (!(cfg.position.w_tilt > 0))
        cfg.position.w_tilt = 1;
    } else if (flag.rfind("wtilt=", 0) == 0) {
      cfg.position.w_tilt = std::atof(flag.c_str() + 6);
    } else if (flag.rfind("tiltaxis=", 0) == 0) {
      cfg.position.tilt.at(0).axis = std::string(flag.c_str() + 9) == "k0"
                                         ? hepworth::block::principal_axis::k0
                                         : hepworth::block::principal_axis::k1;
    } else if (flag.rfind("tiltmode=", 0) == 0) {
      const std::string m = flag.c_str() + 9;
      cfg.position.tilt.at(0).mode = m == "twist"  ? hepworth::block::tilt_mode::twist
                                     : m == "bend" ? hepworth::block::tilt_mode::bend
                                     : m == "spin" ? hepworth::block::tilt_mode::spin
                                                   : hepworth::block::tilt_mode::rigid;
    } else if (flag.rfind("sarap=", 0) == 0) {
      cfg.shape.w_shape_arap = std::atof(flag.c_str() + 6);
    } else if (flag == "sarapfree") {
      cfg.shape.shape_arap_mean = false;
    } else if (flag.rfind("wind=", 0) == 0) {
      // ring:tube[:handedness] on the first stretch entry.
      auto &w = cfg.position.stretch.at(0).winding;
      const std::string v = flag.substr(5);
      const size_t a = v.find(':'), b = a == std::string::npos ? a : v.find(':', a + 1);
      w.ring = std::atoi(v.substr(0, a).c_str());
      w.tube = a == std::string::npos ? 0 : std::atoi(v.substr(a + 1, b - a - 1).c_str());
      w.handedness = b == std::string::npos ? 1 : std::atoi(v.substr(b + 1).c_str());
      if (!(cfg.position.w_stretch > 0))
        cfg.position.w_stretch = 0.5;
    } else if (flag.rfind("wstretch=", 0) == 0) {
      cfg.position.w_stretch = std::atof(flag.c_str() + 9);
    } else if (flag.rfind("dupinreach=", 0) == 0) {
      cfg.shape.dupin_reach = std::atof(flag.c_str() + 11);
    } else if (flag.rfind("dupindeg=", 0) == 0) {
      cfg.shape.dupin_degree = std::atoi(flag.c_str() + 9);
    } else if (flag.rfind("wgradtl=", 0) == 0) {
      cfg.shape.w_grad_traceless = std::atof(flag.c_str() + 8);
    } else if (flag.rfind("wotorus=", 0) == 0) {
      cfg.project.w_osculating_torus = std::atof(flag.c_str() + 8);
    } else if (flag.rfind("wcenter=", 0) == 0) {
      cfg.position.w_center = std::atof(flag.c_str() + 8);
    } else if (flag == "recenter") {
      cfg.position.recenter = true;
    } else if (flag.rfind("solver=", 0) == 0) {
      const std::string v = flag.substr(7);
      using cs = gaudi::duchamp::curvature_solver;
      eig.kind = v == "off" ? cs::off : cs::contact_net;
    } else if (flag.rfind("accel=", 0) == 0) {
      const std::string a = flag.substr(6);
      using sa = gaudi::hepworth::block::step_accel;
      cfg.step.accel = a == "none"                 ? sa::none
                       : a == "chebyshev"          ? sa::chebyshev
                       : a == "anderson_chebyshev" ? sa::anderson_chebyshev
                                                   : sa::anderson;
    } else if (flag.rfind("cntan=", 0) == 0) {
      eig.contact_tangency = flag.substr(6) == "edge" ? gaudi::hepworth::block::contact_net_gn::tangency::edge
                                                      : gaudi::hepworth::block::contact_net_gn::tangency::vertex;
    } else if (flag == "cnframes") {
      eig.contact_frames = true;
    } else if (flag.rfind("nmass=", 0) == 0) {
      eig.contact_normal_mass = std::atof(flag.c_str() + 6);
    } else if (flag.rfind("wdcircle=", 0) == 0) {
      cfg.position.w_dupin_circle = std::atof(flag.c_str() + 9);
    } else if (flag.rfind("dcmode=", 0) == 0) {
      using dm = hepworth::block::dupin_circle_spec::mode;
      const std::string s = flag.substr(7);
      dm m = dm::closest;
      if (s == "slide")
        m = dm::slide;
      else if (s == "even")
        m = dm::even;
      else if (s == "graded")
        m = dm::graded;
      else if (s == "hyperbolic")
        m = dm::hyperbolic;
      else if (s == "parabolic")
        m = dm::parabolic;
      else if (s == "loxodromic")
        m = dm::loxodromic;
      for (auto &d : cfg.position.dupin_circle)
        d.kind = m;
    } else if (flag.rfind("dcfam=", 0) == 0) {
      for (auto &d : cfg.position.dupin_circle)
        d.family = std::atoi(flag.c_str() + 6);
    } else if (flag.rfind("dcamt=", 0) == 0) {
      // The active mode's amount: delta, g, lambda, beta or psi.
      const real a = std::atof(flag.c_str() + 6);
      using dm = hepworth::block::dupin_circle_spec::mode;
      for (auto &d : cfg.position.dupin_circle) {
        if (d.kind == dm::slide)
          d.delta = a;
        else if (d.kind == dm::graded)
          d.g = a;
        else if (d.kind == dm::hyperbolic)
          d.lambda = a;
        else if (d.kind == dm::parabolic)
          d.beta = a;
        else if (d.kind == dm::loxodromic)
          d.psi = a;
      }
    } else if (flag.rfind("dcangle=", 0) == 0) {
      for (auto &d : cfg.position.dupin_circle)
        d.max_angle = std::atof(flag.c_str() + 8) * real(M_PI / 180.0);
    } else if (flag.rfind("dcref=", 0) == 0) {
      real a = 0, b = 0, c = 0;
      std::sscanf(flag.c_str() + 6, "%lf,%lf,%lf", &a, &b, &c);
      cfg.position.dupin_circle_ref = vec3(a, b, c);
    } else if (flag.rfind("wamap=", 0) == 0) {
      cfg.position.w_amap = std::atof(flag.c_str() + 6);
    } else if (flag.rfind("amapinv=", 0) == 0) {
      cfg.position.amap_inversion = std::atof(flag.c_str() + 8);
    } else if (flag.rfind("wconf=", 0) == 0) {
      cfg.position.w_conformal = std::atof(flag.c_str() + 6);
    } else if (flag.rfind("wtension=", 0) == 0) {
      cfg.position.w_tension = std::atof(flag.c_str() + 9);
    } else if (flag.rfind("wkite=", 0) == 0) {
      cfg.position.w_kite = std::atof(flag.c_str() + 6);
    } else if (flag.rfind("kalpha=", 0) == 0) {
      cfg.position.kite.alpha = std::atof(flag.c_str() + 7);
    } else if (flag.rfind("kangle=", 0) == 0) {
      cfg.position.kite.angle = std::atof(flag.c_str() + 7) * real(M_PI) / real(180);
    } else if (flag.rfind("kmaxu=", 0) == 0) {
      cfg.position.kite.max_u = std::atof(flag.c_str() + 6);
    } else if (flag.rfind("ksrc=", 0) == 0) {
      const std::string v = flag.substr(5);
      cfg.position.kite.source = v == "mean"    ? hepworth::block::kite_source::mean_curvature
                                 : v == "aniso" ? hepworth::block::kite_source::anisotropy
                                                : hepworth::block::kite_source::tangent;
    } else if (flag.rfind("klevel=", 0) == 0) {
      cfg.position.kite.level = std::atof(flag.c_str() + 7);
    } else if (flag == "kbalanced") {
      cfg.position.kite.balanced = true;
    } else if (flag.rfind("wcmc=", 0) == 0) {
      cfg.shape.w_tension = std::atof(flag.c_str() + 5);
    } else if (flag.rfind("cmch0=", 0) == 0) {
      cfg.shape.tension.kind = hepworth::block::tension_target::mode::constant;
      cfg.shape.tension.h0 = std::atof(flag.c_str() + 6);
    } else if (flag == "cmcvol") {
      cfg.shape.tension.kind = hepworth::block::tension_target::mode::volume;
    } else if (flag.rfind("wtwist=", 0) == 0) {
      cfg.shape.w_twist = std::atof(flag.c_str() + 7);
    } else if (flag.rfind("tau0=", 0) == 0) {
      cfg.shape.tau0 = std::atof(flag.c_str() + 5);
    } else if (flag.rfind("twangle=", 0) == 0) {
      cfg.shape.twist_angle = std::atof(flag.c_str() + 8);
    } else if (flag.rfind("freshstop=", 0) == 0) {
      cfg.step.gn_stop_on_fresh_ascent = std::atoi(flag.c_str() + 10) != 0;
    } else if (flag.rfind("maxascent=", 0) == 0) {
      cfg.step.gn_max_ascent_refreshes = std::atoi(flag.c_str() + 10);
    } else if (flag.rfind("wmin=", 0) == 0) {
      cfg.shape.w_minimizer = std::atof(flag.c_str() + 5);
    } else if (flag.rfind("dump=", 0) == 0) {
      dump_every = std::atoi(flag.c_str() + 5);
    } else if (flag.rfind("res=", 0) == 0) {
      // ring:tube vertex counts of the built torus.
      const std::string v = flag.substr(4);
      const size_t c = v.find(':');
      cfg.torus_ring = std::atoi(v.substr(0, c).c_str());
      if (c != std::string::npos)
        cfg.torus_tube = std::atoi(v.substr(c + 1).c_str());
    } else if (flag.rfind("perturb=", 0) == 0) {
      cfg.perturb = std::atof(flag.c_str() + 8);
    } else if (flag.rfind("aratio=", 0) == 0) {
      cfg.position.w_axis_ratio = std::atof(flag.c_str() + 7);
    } else if (flag.rfind("bend=", 0) == 0) {
      cfg.position.w_bend = std::atof(flag.c_str() + 5);
    } else if (flag.rfind("damp=", 0) == 0) {
      cfg.step.damping = std::atof(flag.c_str() + 5);
    } else if (flag.rfind("memory=", 0) == 0) {
      cfg.position.w_memory = std::atof(flag.c_str() + 7);
    } else if (flag.rfind("backtracks=", 0) == 0) {
      cfg.step.gn_backtracks = std::atoi(flag.c_str() + 11);
    } else if (flag.rfind("qcwarn=", 0) == 0) {
      cfg.qc.q_warn = std::atof(flag.c_str() + 7);
    } else if (flag == "noqc") {
      cfg.qc.enabled = false;
    } else if (flag.rfind("ascenttol=", 0) == 0) {
      cfg.step.gn_ascent_tol = std::atof(flag.c_str() + 10);
    } else if (flag.rfind("dupin=", 0) == 0) {
      cfg.shape.w_dupin = std::atof(flag.c_str() + 6);
    } else if (flag.rfind("modes=", 0) == 0) {
      cfg.step.preconditioner = hepworth::block::step_preconditioner::modal;
      cfg.step.modal.modes = std::atoi(flag.c_str() + 6);
    } else if (flag.rfind("tau=", 0) == 0) {
      cfg.step.modal.tau = std::atof(flag.c_str() + 4);
    } else if (flag.rfind("lmax=", 0) == 0) {
      cfg.step.modal.lambda_max = std::atof(flag.c_str() + 5);
    } else if (flag == "audit") {
      cfg.step.audit_variational = true;
    } else if (flag == "ascentrefresh") {
      cfg.step.gn_refresh_on_ascent = true;
    } else if (flag.rfind("sdir=", 0) == 0) {
      const std::string v = flag.substr(5);
      for (auto &s : cfg.position.stretch)
        s.direction = v == "material"  ? hepworth::block::direction_source::material
                      : v == "stripes" ? hepworth::block::direction_source::stripes
                                       : hepworth::block::direction_source::angle;
    } else if (flag.rfind("msmooth=", 0) == 0) {
      cfg.position.material_smooth = std::atof(flag.c_str() + 8);
    } else if (flag.rfind("mfollow=", 0) == 0) {
      cfg.position.material_follow = std::atof(flag.c_str() + 8);
    } else if (flag.rfind("stripes=", 0) == 0) {
      cfg.position.stripes = std::atof(flag.c_str() + 8);
    }
  }
  std::cout << "target " << target << "  w_coupling " << eig.w_coupling << "  backtracks " << cfg.step.gn_backtracks << "  refresh "
            << cfg.step.gn_refresh_every << "  pd iterations "
            << cfg.step.iterations << "  remesh " << remesh << "\n";

  auto demo = gaudi::duchamp::monge_fairing_demo::create(cfg);
  std::vector<vec3> h_prev;
  std::vector<real> ring_k_prev;
  for (int f = 0; f < frames; ++f) {
    const std::vector<vec3> before =
        gaudi::asawa::const_get_vec_data(*demo->__M, 0);
    const auto t0 = std::chrono::steady_clock::now();
    demo->step(f, remesh);
    const double wall =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t0)
            .count();
    const std::vector<vec3> &x = gaudi::asawa::const_get_vec_data(*demo->__M, 0);
    real dmax = 0, dmean = 0;
    if (before.size() == x.size()) {
      for (size_t i = 0; i < x.size(); ++i) {
        const real d = (x[i] - before[i]).norm();
        dmax = std::max(dmax, d);
        dmean += d;
      }
      dmean /= real(x.size());
    }
    // Spin: angular momentum of the vertex velocities about the centroid.
    const std::vector<vec3> &vel =
        gaudi::asawa::const_get_vec_data(*demo->__M, 1);
    vec3 cx = vec3::Zero();
    for (const vec3 &p : x)
      cx += p;
    cx /= real(x.size());
    vec3 Lmom = vec3::Zero();
    for (size_t i = 0; i < x.size() && i < vel.size(); ++i)
      Lmom += (x[i] - cx).cross(vel[i]);

    real rho_min = 1e30, rho_max = 0;
    for (const vec3 &p : x) {
      const real rho = std::hypot(p[0], p[1]);
      rho_min = std::min(rho_min, rho);
      rho_max = std::max(rho_max, rho);
    }
    const real Rc = 0.5 * (rho_min + rho_max);
    real tube = 0;
    for (const vec3 &p : x)
      tube += std::hypot(std::hypot(p[0], p[1]) - Rc, p[2]);
    tube /= real(x.size());
    // Total area, and the part on the inner (saddle) side of the ring.
    // Area-weighted centroid too (|center|: drift from the origin).
    real area = 0, area_in = 0;
    vec3 center = vec3::Zero();
    for (auto fi : demo->__M->get_face_range()) {
      const asawa::shell::FaceId fid = asawa::shell::face_id(fi);
      if (demo->__M->fsize(fid) != 3)
        continue;
      const real a = asawa::shell::face_area(*demo->__M, fid, x);
      area += a;
      const auto t = demo->__M->get_tri(fid);
      const vec3 c = (x[t[0]] + x[t[1]] + x[t[2]]) / real(3);
      center += a * c;
      if (std::hypot(c[0], c[1]) < Rc)
        area_in += a;
    }
    if (area > 0)
      center /= area;
    // Symmetry breaking: spread of the mean tube radius over 24 sectors
    // around the ring, relative to the mean (0 on an axisymmetric torus).
    real asym = 0;
    {
      std::array<real, 24> sum{}, cnt{};
      for (const vec3 &p : x) {
        const real phi = std::atan2(p[1], p[0]) + real(M_PI);
        const int b = std::min(23, static_cast<int>(phi / (2 * M_PI) * 24));
        sum[static_cast<size_t>(b)] += std::hypot(std::hypot(p[0], p[1]) - Rc, p[2]);
        cnt[static_cast<size_t>(b)] += 1;
      }
      real m = 0, m2 = 0;
      int nb = 0;
      for (int b = 0; b < 24; ++b)
        if (cnt[static_cast<size_t>(b)] > 0) {
          const real t = sum[static_cast<size_t>(b)] / cnt[static_cast<size_t>(b)];
          m += t;
          m2 += t * t;
          ++nb;
        }
      if (nb > 0) {
        m /= nb;
        asym = std::sqrt(std::max(m2 / nb - m * m, real(0))) / std::max(m, real(1e-12));
      }
    }

    std::vector<kusama::tri3> tris;
    for (auto fi : demo->__M->get_face_range()) {
      const asawa::shell::FaceId fid = asawa::shell::face_id(fi);
      if (demo->__M->fsize(fid) != 3)
        continue;
      const auto t = demo->__M->get_tri(fid);
      tris.push_back({static_cast<int>(t[0]), static_cast<int>(t[1]), static_cast<int>(t[2])});
    }
    // Mode spectrum of the normal offset from the fitted torus: per vertex
    // `d = |p - ring point| - tube` at ring angle `phi` and tube angle
    // `theta`, area-weighted Fourier amplitudes `2 |A(m, n)|` for m (ring)
    // 0..12 and n (tube) -6..6, relative to the tube radius. A buckle is a few
    // low modes that grow; noise is a flat spread at high m and n.
    struct mode_amp {
      int m = 0, n = 0;
      real a = 0;
    };
    std::vector<mode_amp> modes;
    real na_rms = 0, axi_best = 0;
    int axi_n = 0;
    if (tube > real(0) && Rc > tube) {
      std::vector<real> va(x.size(), real(0));
      for (const auto &t : tris) {
        const real a = real(0.5) * (x[static_cast<size_t>(t[1])] - x[static_cast<size_t>(t[0])])
                                       .cross(x[static_cast<size_t>(t[2])] - x[static_cast<size_t>(t[0])])
                                       .norm();
        for (int k = 0; k < 3; ++k)
          va[static_cast<size_t>(t[k])] += a / real(3);
      }
      std::vector<real> ph(x.size()), th(x.size()), d(x.size());
      real wsum = 0;
      for (size_t i = 0; i < x.size(); ++i) {
        const vec3 &p = x[i];
        const real rho = std::hypot(p[0], p[1]);
        ph[i] = std::atan2(p[1], p[0]);
        th[i] = std::atan2(p[2], rho - Rc);
        d[i] = std::hypot(rho - Rc, p[2]) - tube;
        wsum += va[i];
      }
      for (int m = 0; m <= 12; ++m)
        for (int n = -6; n <= 6; ++n) {
          if (m == 0 && n <= 0)
            continue;
          real re = 0, im = 0;
          for (size_t i = 0; i < x.size(); ++i) {
            const real arg = m * ph[i] + n * th[i];
            re += va[i] * d[i] * std::cos(arg);
            im -= va[i] * d[i] * std::sin(arg);
          }
          const real amp = real(2) * std::hypot(re, im) / std::max(wsum, real(1e-30)) / tube;
          if (m == 0) {
            if (n >= 2 && amp > axi_best) {
              axi_best = amp;
              axi_n = n;
            }
            continue;
          }
          na_rms += real(0.5) * amp * amp;
          modes.push_back({m, n, amp});
        }
      na_rms = std::sqrt(na_rms);
      std::sort(modes.begin(), modes.end(),
                [](const mode_amp &a, const mode_amp &b) { return a.a > b.a; });
      if (modes.size() > 3)
        modes.resize(3);
    }
    if (dump_every > 0 && ((f + 1) % dump_every == 0 || f + 1 == frames)) {
      char path[64];
      std::snprintf(path, sizeof(path), "/tmp/probe_f%02d.obj", f);
      if (FILE *o = std::fopen(path, "w")) {
        for (const vec3 &p : x)
          std::fprintf(o, "v %.9g %.9g %.9g\n", p[0], p[1], p[2]);
        for (const auto &t : tris)
          std::fprintf(o, "f %d %d %d\n", t[0] + 1, t[1] + 1, t[2] + 1);
        std::fclose(o);
      }
    }

    real rough = 0;
    int ne = 0;
    {
      gaudi::kusama::principal_curvature_solver s;
      s.reset(static_cast<int>(x.size()), tris, x);
      const std::vector<vec2> k = s.principal(x);
      for (const auto &t : tris)
        for (int a = 0; a < 3; ++a) {
          const int i = t[a], j = t[(a + 1) % 3];
          if (i < j && k[i].allFinite() && k[j].allFinite()) {
            rough += std::abs((k[i][0] + k[i][1]) - (k[j][0] + k[j][1]));
            ++ne;
          }
        }
      rough /= std::max(ne, 1);
    }
    // Live edge readings on the solved shape: mean log ratio log(|k_hi| /
    // |k_lo|) (capped at log 1000) and anisotropy q in [0, 1] on elliptic
    // (k0 k1 > 0) and saddle edges, and the mean turn of the principal axes
    // since the previous frame.
    real rho_ell = 0, rho_sad = 0, an_ell = 0, an_sad = 0, turn = 0, turn2 = 0, skew = 0;
    real kdrift = std::numeric_limits<real>::quiet_NaN();
    int n_ell = 0, n_sad = 0, nt = 0;
    const auto &tb = demo->cr()->table;
    if (tb.n_verts == static_cast<int>(x.size())) {
      const std::vector<vec3> h = kusama::weingarten::edge_readings(tb, x);
      {
        // RMS of the per-face, per-axis curvature gradients along their own
        // lines (the Dupin rows' residual on the live read-out).
        std::vector<vec3> hm(h.size(), vec3::Ones()), hz(h.size(), vec3::Zero());
        std::vector<vec3> hh = h;
        auto blk = hepworth::vec3_block::create(hm, hh, hz, hz);
        auto row = hepworth::block::edge_shape_fair::create(
            &tb, x, hepworth::block::edge_fair_rows::dupin, real(1),
            std::vector<hepworth::sim_block::ptr>{blk}, hepworth::block::goal_zero(), real(1), &h);
        index_t id0 = 0;
        std::vector<trip> trips;
        row->fill_A(id0, trips);
        vecX qh;
        blk->map_to_x(qh);
        vecX r = vecX::Zero(id0);
        for (const trip &t : trips)
          r[t.row()] += t.value() * qh[t.col()];
        skew = id0 > 0 ? r.norm() / std::sqrt(real(id0)) : real(0);
      }
      {
        // Each ring's sum K dA on the live read-out, against the last frame's.
        const std::vector<real> c = kusama::weingarten::coupling_factor(tb, x);
        const auto rings = hepworth::block::vertex_edge_rings(tb);
        std::vector<real> ring_k(rings.size(), real(0));
        for (size_t v = 0; v < rings.size(); ++v)
          for (const index_t e : rings[v])
            ring_k[v] += c[static_cast<size_t>(e)] *
                         (h[static_cast<size_t>(e)][0] * h[static_cast<size_t>(e)][2] -
                          h[static_cast<size_t>(e)][1] * h[static_cast<size_t>(e)][1]);
        if (ring_k_prev.size() == ring_k.size()) {
          real num = 0, den = 0;
          for (size_t v = 0; v < ring_k.size(); ++v) {
            num += std::abs(ring_k[v] - ring_k_prev[v]);
            den += std::abs(ring_k_prev[v]);
          }
          kdrift = den > real(0) ? num / den : real(0);
        }
        ring_k_prev = ring_k;
      }
      for (size_t e = 0; e < h.size(); ++e) {
        const auto p = hepworth::block::principals_of(h[e]);
        real k0, k1;
        hepworth::block::signed_principals(h[e], k0, k1);
        const real rho = hepworth::block::curvature_measure(
            hepworth::block::curvature_metric::log_ratio, k0, k1);
        const real an = hepworth::block::curvature_measure(
            hepworth::block::curvature_metric::aniso, k0, k1);
        if (k0 * k1 > 0) {
          rho_ell += rho;
          an_ell += an;
          ++n_ell;
        } else {
          rho_sad += rho;
          an_sad += an;
          ++n_sad;
        }
        if (h_prev.size() == h.size()) {
          const auto q = hepworth::block::principals_of(h_prev[e]);
          if (p.ok && q.ok && std::abs(q.k_hi) > real(1.2) * std::abs(q.k_lo)) {
            const real t = std::remainder(std::atan2(p.hy, p.hx) - std::atan2(q.hy, q.hx), real(M_PI));
            turn += t;
            turn2 += t * t;
            ++nt;
          }
        }
      }
      h_prev = h;
    }
    std::cout << "frame " << f << "  |dx| mean " << dmean << " max " << dmax
              << "  R " << Rc << " tube " << tube << "  area " << area << " (inner " << area_in
              << ")  |center| " << center.norm() << "  rough " << rough
              << "  log r ell " << rho_ell / std::max(n_ell, 1) << " sad "
              << rho_sad / std::max(n_sad, 1) << "  aniso ell " << an_ell / std::max(n_ell, 1)
              << " sad " << an_sad / std::max(n_sad, 1) << " (" << n_ell << "/" << n_sad << ")"
              << "  skew " << skew << "  kdrift " << kdrift << "  turn " << turn / std::max(nt, 1) << " sd "
              << std::sqrt(std::max(turn2 / std::max(nt, 1) -
                                        std::pow(turn / std::max(nt, 1), 2),
                                    real(0)))
              << "  asym " << asym << "  dcres " << demo->dupin_circle_residual()
              << "  kiteu " << demo->kite_u_range()
              << "  defl " << demo->material_deflection()
              << "  stripes " << demo->last_stripes().lambda << " skip "
              << demo->last_stripes().skipped
              << "  |L| " << Lmom.norm() << "  verts " << x.size()
              << "  wall " << wall << " s"
              << "  solve its " << demo->last_solve().its << " E " << demo->last_solve().E_first
              << " -> " << demo->last_solve().E_last << " |dq| " << demo->last_solve().dq_last
              << "  lambda_min " << demo->last_solve().lambda_min << " neg "
              << demo->last_solve().modal_negative << " accepted "
              << demo->last_solve().modal_accepted << "/" << demo->last_solve().modal_tried
              << "  refactors " << demo->last_solve().refreshes << "  accel " << demo->last_solve().accel_taken
              << "/" << demo->last_solve().accel_tried << "\n";
    std::cout << "  modes: non-axisymmetric rms " << na_rms << "  top";
    for (const mode_amp &md : modes)
      std::cout << " (" << md.m << "," << md.n << ") " << md.a;
    std::cout << "  axisymmetric n=" << axi_n << " " << axi_best << "\n";
    if (remesh && cfg.qc.enabled) {
      const asawa::shell::qc_report &qc = demo->__surf->last_qc();
      std::cout << "  qc after flip f " << qc.flip.faces << " v " << qc.flip.val_lo << "/"
                << qc.flip.val_hi << "  fixed c/f/d/t " << qc.collapsed << "/" << qc.flipped << "/"
                << qc.dissolved << "/" << qc.retriangulated << "  left f " << qc.final.faces
                << " hard " << qc.final.hard << " v " << qc.final.val_lo << "/" << qc.final.val_hi
                << " worst " << qc.final.worst_q << "\n";
    }
  }
  return 0;
}

} // namespace gaudi

int main(int argc, char **argv) {
  return gaudi::principal_curvature_probe_main(argc, argv);
}
