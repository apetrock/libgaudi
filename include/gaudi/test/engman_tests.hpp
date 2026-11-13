#ifndef GAUDI_ENGMAN_TESTS_HPP
#define GAUDI_ENGMAN_TESTS_HPP

/// The Engman solver (`engman/solver.hpp`, MATH.md section 17). Included
/// after weingarten_tests.hpp, whose helpers (torus_tris, cube_tris,
/// rigid_error, run_hepworth_frames) it uses.

#include "gaudi/engman/constraint.hpp"
#include "gaudi/engman/energies.hpp"
#include "gaudi/engman/solver.hpp"
#include "gaudi/engman/store_metric.hpp"
#include "gaudi/engman/scene.hpp"
#include "gaudi/test/test.hpp"
#include "gaudi/test/weingarten_tests.hpp"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <vector>

namespace gaudi {
namespace test {
namespace {

/// One scene: positions and the edge store on `tris`, unit vertex masses,
/// the memory and a store pin to `target` as energies.
struct engman_scene {
  kusama::stencil_table t;
  std::vector<vec3> x, v, f, m, hs, hv, hf, hm;
  hepworth::sim_block::ptr pos, store;
  std::vector<hepworth::sim_block::ptr> blocks;

  engman_scene(const std::vector<vec3> &x0, const std::vector<std::array<int, 3>> &tris)
      : t(kusama::build_stencil_table(asawa::shell::triangle_walker(static_cast<int>(x0.size()), tris))),
        x(x0) {
    const int n = static_cast<int>(x.size());
    v.assign(static_cast<size_t>(n), vec3::Zero());
    f = v;
    m.assign(static_cast<size_t>(n), vec3::Ones());
    pos = hepworth::vec3_block::create(m, x, v, f);
    store = hepworth::vec3_block::create(hm, hs, hv, hf);
    blocks = {pos, store};
    reset_store();
  }
  /// The frame start: the store reads the live mesh, velocities dropped.
  void reset_store() {
    hs = wg::edge_readings(t, x);
    hv.assign(hs.size(), vec3::Zero());
    hf.assign(hs.size(), vec3::Zero());
    hm.assign(hs.size(), vec3::Constant(real(1e-8)));
    for (vec3 &vi : v)
      vi.setZero();
  }
  std::vector<hepworth::projection_constraint::ptr> energies(real w_m, real w_pin,
                                                             const std::vector<vec3> &target) {
    return {hb::dirichlet_memory_gn::create(&t, x, w_m, std::vector<hepworth::sim_block::ptr>{pos}),
            hb::edge_shape_fair::create(&t, x, hb::edge_fair_rows::pin, w_pin,
                                        std::vector<hepworth::sim_block::ptr>{store}, hb::goal_fixed(target))};
  }
  real size() const {
    real s = 0;
    for (const vec3 &p : x)
      s = std::max(s, p.cwiseAbs().maxCoeff());
    return s;
  }
};

inline real max_diff(const std::vector<vec3> &a, const std::vector<vec3> &b) {
  real d = 0;
  for (size_t i = 0; i < a.size(); ++i)
    d = std::max(d, (a[i] - b[i]).cwiseAbs().maxCoeff());
  return d;
}

/// Max |C| of the coupling relative to the size of its terms `c |R|`.
inline real coupling_feasibility(engman_scene &sc) {
  const engman::constraint::ptr c = engman::coupling::create(&sc.t, sc.x, sc.blocks);
  vecX q, C;
  for (auto &b : sc.blocks)
    b->map_to_x(q);
  c->evaluate(q, C, nullptr);
  real scale = 0;
  const std::vector<vec3> R = wg::edge_readings(sc.t, sc.x);
  const std::vector<real> cf = wg::coupling_factor(sc.t, sc.x);
  for (size_t e = 0; e < R.size(); ++e)
    scale = std::max(scale, cf[e] * R[e].cwiseAbs().maxCoeff());
  return C.cwiseAbs().maxCoeff() / std::max(scale, real(1e-30));
}

/// With the duals off and `rho = w^2`, Engman minimizes the same energy as
/// `projection_solver` with a `w`-weighted `weingarten_gn` row: one frame on
/// the torus tube lands on the same positions and store.
GAUDI_TEST(engman_penalty_parity) {
  std::vector<vec3> x0, xt;
  std::vector<std::array<int, 3>> tris;
  torus_tris(2.0, 0.7, 24, 12, x0, tris);
  torus_tris(2.0, 0.9, 24, 12, xt, tris);
  const real w = 10, w_m = 1, w_pin = 1e2;
  engman_scene a(x0, tris), b(x0, tris);
  const std::vector<vec3> target = wg::edge_readings(a.t, xt);
  {
    hb::projection_solver ps;
    ps.gn_refresh_every = 1;
    ps.gn_backtracks = 8;
    std::vector<hepworth::projection_constraint::ptr> cs = a.energies(w_m, w_pin, target);
    cs.push_back(hb::weingarten_gn::create(&a.t, a.x, w, a.blocks));
    ps.set_constraints(cs);
    ps.step(a.blocks, real(0.1), real(1), 200);
  }
  engman::solver es;
  es.cfg.rho = w * w;
  es.cfg.dual_updates = false;
  es.set_energies(b.energies(w_m, w_pin, target));
  es.set_constraints({engman::coupling::create(&b.t, b.x, b.blocks)});
  es.step(b.blocks, real(0.1), real(1), 200);
  const real dx = max_diff(a.x, b.x) / a.size();
  real hs = 0;
  for (const vec3 &h : a.hs)
    hs = std::max(hs, h.cwiseAbs().maxCoeff());
  const real dh = max_diff(a.hs, b.hs) / hs;
  console_logger::info << "engman penalty parity: positions " << dx << ", store " << dh << " (engman "
                       << es.last.iterations << " iterates, " << es.last.refreshes << " factors)" << std::endl;
  GAUDI_EXPECT(dx < real(1e-6));
  GAUDI_EXPECT(dh < real(1e-6));
}

/// With the duals on the coupling is met (the store is the mesh's reading),
/// and the result does not depend on `rho`.
GAUDI_TEST(engman_coupling_is_met) {
  std::vector<vec3> x0, xt;
  std::vector<std::array<int, 3>> tris;
  torus_tris(2.0, 0.7, 24, 12, x0, tris);
  torus_tris(2.0, 0.9, 24, 12, xt, tris);
  std::vector<std::vector<vec3>> xs;
  std::vector<real> feas;
  for (const real rho_rel : {real(1e2), real(1e4), real(1e6)}) {
    engman_scene sc(x0, tris);
    const std::vector<vec3> target = wg::edge_readings(sc.t, xt);
    engman::solver es;
    es.cfg.rho_rel = rho_rel;
    es.set_energies(sc.energies(1, 1e2, target));
    es.set_constraints({engman::coupling::create(&sc.t, sc.x, sc.blocks)});
    es.step(sc.blocks, real(0.1), real(1), 200);
    feas.push_back(coupling_feasibility(sc));
    console_logger::info << "engman coupling rho_rel " << rho_rel << ": |C| " << es.last.C_first << " -> "
                         << es.last.C_last << ", relative " << feas.back() << ", iterates "
                         << es.last.iterations << ", factors " << es.last.refreshes << std::endl;
    xs.push_back(sc.x);
  }
  const real size = engman_scene(x0, tris).size();
  const real spread = std::max(max_diff(xs[0], xs[1]), max_diff(xs[1], xs[2])) / size;
  console_logger::info << "engman coupling: rho spread " << spread << std::endl;
  for (real fz : feas)
    GAUDI_EXPECT(fz < real(1e-10));
  GAUDI_EXPECT(spread < real(1e-7));
}

/// Constraints are met together, not traded: holding a ring of vertices
/// exactly leaves the coupling met; as penalties the two pull against each
/// other and neither holds.
GAUDI_TEST(engman_constraints_are_met_together) {
  std::vector<vec3> x0, xt;
  std::vector<std::array<int, 3>> tris;
  torus_tris(2.0, 0.7, 24, 12, x0, tris);
  torus_tris(2.0, 0.9, 24, 12, xt, tris);
  std::vector<index_t> ids;
  std::vector<vec3> goals;
  for (int i = 0; i < 24; ++i) {
    ids.push_back(i * 12);
    goals.push_back(x0[static_cast<size_t>(i * 12)]);
  }
  auto run = [&](bool hard, real &pin_err, real &feas) {
    engman_scene sc(x0, tris);
    const std::vector<vec3> target = wg::edge_readings(sc.t, xt);
    engman::solver es;
    es.cfg.dual_updates = hard;
    if (!hard)
      es.cfg.rho = 1e2;
    es.set_energies(sc.energies(1, 1e2, target));
    es.set_constraints({engman::coupling::create(&sc.t, sc.x, sc.blocks),
                        engman::position_pin::create(ids, goals, sc.pos)});
    es.step(sc.blocks, real(0.1), real(1), 200);
    pin_err = 0;
    for (size_t k = 0; k < ids.size(); ++k)
      pin_err = std::max(pin_err, (sc.x[static_cast<size_t>(ids[k])] - goals[k]).norm());
    feas = coupling_feasibility(sc);
  };
  real pin_hard, feas_hard, pin_soft, feas_soft;
  run(true, pin_hard, feas_hard);
  run(false, pin_soft, feas_soft);
  console_logger::info << "engman two constraints: hard pin " << pin_hard << " coupling " << feas_hard
                       << "; as penalties pin " << pin_soft << " coupling " << feas_soft << std::endl;
  GAUDI_EXPECT(pin_hard < real(1e-12) && feas_hard < real(1e-10));
  GAUDI_EXPECT(pin_soft > real(1e-6) && feas_soft > real(1e-6));
}

/// Realizable targets over frames (each frame resets the store to the
/// reading and the memory rest to the frame start), against the
/// `projection_solver` runs of the same scenes.
GAUDI_TEST(engman_realizable_targets_over_frames) {
  struct scene {
    const char *name;
    std::vector<vec3> x0, xt;
    std::vector<std::array<int, 3>> tris;
    real w_c, w_pin;
  };
  std::vector<scene> rows(2);
  rows[0].name = "torus tube 0.7 -> 0.9";
  torus_tris(2.0, 0.7, 24, 12, rows[0].x0, rows[0].tris);
  torus_tris(2.0, 0.9, 24, 12, rows[0].xt, rows[0].tris);
  rows[0].w_c = 1e3;
  rows[0].w_pin = 1e4;
  rows[1].name = "cube -> sphere";
  cube_tris(6, rows[1].x0, rows[1].tris);
  {
    const real R = std::sqrt(real(24) / (4 * M_PI));
    for (const vec3 &p : rows[1].x0)
      rows[1].xt.push_back(R * p.normalized());
  }
  rows[1].w_c = 1e4;
  rows[1].w_pin = 1e5;
  for (const scene &rw : rows) {
    const real e0 = rigid_error(rw.x0, rw.xt);
    const real ep = run_hepworth_frames(rw.x0, rw.xt, rw.tris, 20, rw.w_c, real(1), rw.w_pin);
    engman_scene sc(rw.x0, rw.tris);
    const std::vector<vec3> target = wg::edge_readings(sc.t, rw.xt);
    engman::solver es;
    real worst_feas = 0;
    for (int k = 0; k < 20; ++k) {
      sc.reset_store();
      es.set_energies(sc.energies(1, rw.w_pin, target));
      es.set_constraints({engman::coupling::create(&sc.t, sc.x, sc.blocks)});
      es.step(sc.blocks, real(0.1), real(1), 30);
      const real fz = coupling_feasibility(sc);
      if (fz > real(1e-8))
        console_logger::info << "  frame " << k << ": feasibility " << fz << ", iterates " << es.last.iterations
                             << ", restorations " << es.last.restorations << ", stopped " << es.last.stopped
                             << std::endl;
      worst_feas = std::max(worst_feas, fz);
    }
    const real ee = rigid_error(sc.x, rw.xt);
    console_logger::info << "engman " << rw.name << " over 20 frames: rigid err " << e0 << " -> engman " << ee
                         << ", projection " << ep << "; worst coupling " << worst_feas << std::endl;
    GAUDI_EXPECT(ee <= std::max(ep, real(1e-10) * e0));
    GAUDI_EXPECT(worst_feas < real(1e-8));
  }
}

/// `|A h|` of `edge_shape_fair` rows with a zero goal over the field `h`.
inline real fair_rows_norm(const kusama::stencil_table &t, const std::vector<vec3> &x, hb::edge_fair_rows rows,
                           const std::vector<vec3> &h) {
  std::vector<vec3> hs = h, hm(h.size(), vec3::Ones()), hv(h.size(), vec3::Zero()), hf = hv;
  auto blk = hepworth::vec3_block::create(hm, hs, hv, hf);
  auto c = hb::edge_shape_fair::create(&t, x, rows, real(1), std::vector<hepworth::sim_block::ptr>{blk},
                                       hb::goal_zero());
  vecX q;
  blk->map_to_x(q);
  index_t id0 = 0;
  std::vector<trip> trips;
  c->fill_A(id0, trips);
  wg::matS A(id0, q.size());
  A.setFromTriplets(trips.begin(), trips.end());
  return (A * q).norm();
}

/// The spin-2 element splits exactly into its traceless part and the trace
/// part (`1/2 K` on `tau`); the row energies add the same way; the
/// traceless rows vanish on a constant planar field and on any umbilic
/// field (`z = 0`), whose varying mean curvature only the trace rows see.
GAUDI_TEST(engman_spin2_trace_traceless_split) {
  std::mt19937 rng(5);
  std::normal_distribution<real> nd(0, 1);
  const vec3 tt(1, 0, 1);
  real e_el = 0;
  for (int k = 0; k < 20; ++k) {
    const vec3 X0(nd(rng), nd(rng), nd(rng)), X1(nd(rng), nd(rng), nd(rng)), X2(nd(rng), nd(rng), nd(rng));
    const Eigen::Matrix<real, 9, 9> E = wg::spin2_face_element(X0, X1, X2);
    const Eigen::Matrix<real, 9, 9> Ez = wg::spin2_traceless_face_element(X0, X1, X2);
    const Eigen::Matrix<real, 3, 3> K = wg::cr_face_stiffness(X0, X1, X2);
    Eigen::Matrix<real, 9, 9> Et;
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b)
        Et.block<3, 3>(3 * a, 3 * b) = real(0.5) * K(a, b) * tt * tt.transpose();
    e_el = std::max(e_el, (E - Ez - Et).cwiseAbs().maxCoeff() / E.cwiseAbs().maxCoeff());
  }
  store_fixture fx;
  const real full = fair_rows_norm(fx.t, fx.x, hb::edge_fair_rows::dirichlet, fx.hs);
  const real tr = fair_rows_norm(fx.t, fx.x, hb::edge_fair_rows::dirichlet_trace, fx.hs);
  const real tl = fair_rows_norm(fx.t, fx.x, hb::edge_fair_rows::dirichlet_traceless, fx.hs);
  const real e_rows = std::abs(full * full - tr * tr - tl * tl) / (full * full);
  std::vector<vec3> xg;
  std::vector<std::array<int, 3>> tris;
  unit_grid(9, xg, tris);
  const kusama::stencil_table tg =
      kusama::build_stencil_table(asawa::shell::triangle_walker(static_cast<int>(xg.size()), tris));
  const std::vector<vec3> hc = sample_edges_planar(tg, xg, [](const vec3 &) {
    mat3 S;
    S << 1.0, 0.3, 0, 0.3, -0.5, 0, 0, 0, 0;
    return S;
  });
  std::vector<vec3> hr(hc.size());
  for (vec3 &v : hr)
    v = vec3(nd(rng), nd(rng), nd(rng));
  const real planar = fair_rows_norm(tg, xg, hb::edge_fair_rows::dirichlet_traceless, hc) /
                      fair_rows_norm(tg, xg, hb::edge_fair_rows::dirichlet_traceless, hr);
  std::vector<vec3> hu(fx.hs.size());
  for (vec3 &v : hu) {
    const real k = real(1) + real(0.3) * nd(rng);
    v = vec3(k, 0, k);
  }
  const real umb_tl = fair_rows_norm(fx.t, fx.x, hb::edge_fair_rows::dirichlet_traceless, hu);
  const real umb_tr = fair_rows_norm(fx.t, fx.x, hb::edge_fair_rows::dirichlet_trace, hu);
  console_logger::info << "spin-2 split: element " << e_el << ", rows (full^2 - trace^2 - traceless^2) / full^2 "
                       << e_rows << ", traceless on a constant planar field / random " << planar
                       << ", umbilic field traceless " << umb_tl << " trace " << umb_tr << std::endl;
  GAUDI_EXPECT(e_el < real(1e-13));
  GAUDI_EXPECT(e_rows < real(1e-10));
  GAUDI_EXPECT(planar < real(1e-8));
  GAUDI_EXPECT(umb_tl < real(1e-10) * umb_tr && umb_tr > real(0));
}

/// The store metric (`store_metric`) shapes the steps and not the solution:
/// with it the coupled solve lands on the same positions and store, and the
/// first step's store change is smoother (lower spin-2 energy per norm).
GAUDI_TEST(engman_store_metric_keeps_the_solution) {
  std::vector<vec3> x0, xt;
  std::vector<std::array<int, 3>> tris;
  torus_tris(2.0, 0.7, 24, 12, x0, tris);
  torus_tris(2.0, 0.9, 24, 12, xt, tris);
  std::mt19937 rng(9);
  std::normal_distribution<real> nd(0, 1);
  auto run = [&](real lambda, int iterations, std::vector<vec3> &x_out, std::vector<vec3> &dS) {
    engman_scene sc(x0, tris);
    std::vector<vec3> target = wg::edge_readings(sc.t, xt);
    for (vec3 &g : target)
      g += real(0.05) * vec3(nd(rng), nd(rng), nd(rng));
    const std::vector<vec3> S0 = sc.hs;
    engman::solver es;
    std::vector<hepworth::projection_constraint::ptr> en = sc.energies(1, 1e2, target);
    for (auto &e : engman::store_metric(&sc.t, sc.x, sc.store, lambda, lambda))
      en.push_back(e);
    es.set_energies(en);
    es.set_constraints({engman::coupling::create(&sc.t, sc.x, sc.blocks)});
    es.step(sc.blocks, real(0.1), real(1), iterations);
    x_out = sc.x;
    dS.resize(S0.size());
    for (size_t e = 0; e < S0.size(); ++e)
      dS[e] = sc.hs[e] - S0[e];
    return sc.size();
  };
  std::vector<vec3> xa, xb, da, db;
  rng.seed(9);
  const real size = run(0, 200, xa, da);
  rng.seed(9);
  run(1e2, 200, xb, db);
  const real same = max_diff(xa, xb) / size;
  rng.seed(9);
  run(0, 1, xa, da);
  rng.seed(9);
  run(1e2, 1, xb, db);
  engman_scene sc(x0, tris);
  auto rough = [&](const std::vector<vec3> &d) {
    real n = 0;
    for (const vec3 &v : d)
      n += v.squaredNorm();
    return fair_rows_norm(sc.t, sc.x, hb::edge_fair_rows::dirichlet, d) / std::sqrt(std::max(n, real(1e-300)));
  };
  const real ra = rough(da), rb = rough(db);
  console_logger::info << "engman store metric: solution change " << same << ", first-step roughness " << ra
                       << " -> " << rb << std::endl;
  GAUDI_EXPECT(same < real(1e-8));
  GAUDI_EXPECT(rb < ra);
}

/// The Schur reduction onto positions lands where the full KKT does: with
/// only the coupling (positions-only SPD), with an unpaired position pin
/// (reduced quasi-definite), and `automatic` falls back to full when a
/// Dupin row couples store edges.
GAUDI_TEST(engman_schur_matches_full) {
  std::vector<vec3> x0, xt;
  std::vector<std::array<int, 3>> tris;
  torus_tris(2.0, 0.7, 24, 12, x0, tris);
  torus_tris(2.0, 0.9, 24, 12, xt, tris);
  std::vector<index_t> ids;
  std::vector<vec3> goals;
  for (int i = 0; i < 24; ++i) {
    ids.push_back(i * 12);
    goals.push_back(x0[static_cast<size_t>(i * 12)]);
  }
  struct out {
    std::vector<vec3> x, h;
    bool schur;
  };
  auto run = [&](engman::kkt_mode mode, bool pin, bool dupin) {
    engman_scene sc(x0, tris);
    const std::vector<vec3> target = wg::edge_readings(sc.t, xt);
    std::vector<hepworth::projection_constraint::ptr> en = sc.energies(1, 1e2, target);
    if (dupin) {
      const std::vector<vec3> h0 = sc.hs;
      en.push_back(hb::edge_shape_fair::create(&sc.t, sc.x, hb::edge_fair_rows::dupin, real(1),
                                               std::vector<hepworth::sim_block::ptr>{sc.store}, hb::goal_zero(),
                                               real(1), &h0));
    }
    std::vector<engman::constraint::ptr> cs{engman::coupling::create(&sc.t, sc.x, sc.blocks)};
    if (pin)
      cs.push_back(engman::position_pin::create(ids, goals, sc.pos));
    engman::solver es;
    es.cfg.kkt = mode;
    es.set_energies(en);
    es.set_constraints(cs);
    es.step(sc.blocks, real(0.1), real(1), 200);
    return out{sc.x, sc.hs, es.last.schur};
  };
  const real size = engman_scene(x0, tris).size();
  for (const bool pin : {false, true}) {
    const out f = run(engman::kkt_mode::full, pin, false), s = run(engman::kkt_mode::schur, pin, false);
    const real dx = max_diff(f.x, s.x) / size, dh = max_diff(f.h, s.h);
    console_logger::info << "engman schur vs full" << (pin ? " (with a position pin)" : "") << ": positions "
                         << dx << ", store " << dh << ", schur used " << s.schur << std::endl;
    GAUDI_EXPECT(s.schur && !f.schur);
    GAUDI_EXPECT(dx < real(1e-10) && dh < real(1e-8));
  }
  const out a = run(engman::kkt_mode::automatic, false, true), f = run(engman::kkt_mode::full, false, true);
  console_logger::info << "engman automatic with a Dupin row: schur used " << a.schur << ", vs full "
                       << max_diff(a.x, f.x) / size << std::endl;
  GAUDI_EXPECT(!a.schur);
  GAUDI_EXPECT(max_diff(a.x, f.x) / size < real(1e-12));
}

/// Factor cost, reported (not asserted): the projection solver's matrix
/// (`M / h^2 + A^T A` with the `weingarten_gn` rows), the full KKT matrix,
/// and the Schur reduction, on the demo torus and a finer one.
GAUDI_TEST(engman_kkt_cost) {
  for (const auto [nu, nv] : {std::pair{48, 24}, std::pair{128, 64}}) {
    std::vector<vec3> x0, xt;
    std::vector<std::array<int, 3>> tris;
    torus_tris(1.0, 0.33, nu, nv, x0, tris);
    torus_tris(1.0, 0.4, nu, nv, xt, tris);
    engman_scene sp(x0, tris);
    const std::vector<vec3> target = wg::edge_readings(sp.t, xt);
    double t_proj = 0;
    index_t nnz_proj = 0, n_proj = 0;
    {
      std::vector<hepworth::projection_constraint::ptr> cs = sp.energies(1, 1e2, target);
      cs.push_back(hb::weingarten_gn::create(&sp.t, sp.x, real(10), sp.blocks));
      vecX q, M;
      for (auto &b : sp.blocks) {
        b->map_to_x(q);
        b->map_mass(M);
      }
      std::vector<trip> tr;
      index_t id0 = 0;
      for (auto &c : cs) {
        if (c->gauss_newton())
          c->fill_J(q, id0, tr);
        else
          c->fill_A(id0, tr);
      }
      wg::matS A(id0, q.size());
      A.setFromTriplets(tr.begin(), tr.end());
      std::vector<trip> mt;
      for (index_t i = 0; i < q.size(); ++i)
        mt.emplace_back(i, i, std::max(M[i], real(1e-8)) / real(0.01));
      wg::matS Mh(q.size(), q.size());
      Mh.setFromTriplets(mt.begin(), mt.end());
      const wg::matS H = Mh + wg::matS(A.transpose() * A);
      const auto t0 = std::chrono::steady_clock::now();
      Eigen::SimplicialLDLT<wg::matS> L(H);
      t_proj = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      nnz_proj = L.matrixL().nestedExpression().nonZeros();
      n_proj = H.rows();
    }
    auto engman_cost = [&](engman::kkt_mode mode, double &t, index_t &nnz, index_t &n) {
      engman_scene sc(x0, tris);
      engman::solver es;
      es.cfg.kkt = mode;
      es.cfg.restore = 0;
      es.set_energies(sc.energies(1, 1e2, target));
      es.set_constraints({engman::coupling::create(&sc.t, sc.x, sc.blocks)});
      es.step(sc.blocks, real(0.1), real(1), 1);
      t = es.last.factor_seconds / std::max(es.last.refreshes, 1);
      nnz = es.last.factor_nonzeros;
      n = es.last.factored_size;
    };
    double t_full = 0, t_schur = 0;
    index_t nnz_full = 0, nnz_schur = 0, n_full = 0, n_schur = 0;
    engman_cost(engman::kkt_mode::full, t_full, nnz_full, n_full);
    engman_cost(engman::kkt_mode::schur, t_schur, nnz_schur, n_schur);
    console_logger::info << "engman kkt cost, torus " << nu << " x " << nv << " (" << x0.size()
                         << " verts): projection n " << n_proj << " L " << nnz_proj << " " << t_proj
                         << " s; full KKT n " << n_full << " L " << nnz_full << " " << t_full << " s; schur n "
                         << n_schur << " L " << nnz_schur << " " << t_schur << " s" << std::endl;
    GAUDI_EXPECT(nnz_schur > 0 && nnz_full > 0);
  }
}

// ------------------------------------------------- shape-space twist pieces

/// The edge store of a planar field on the unit grid (normal +z): trace
/// `tau`, traceless `z0 e^{i omega x}` in world axes, each edge's tensor in
/// its frame `(t_e, z x t_e)`.
inline std::vector<vec3> planar_twisted_store(const kusama::stencil_table &t, const std::vector<vec3> &x,
                                              const std::vector<std::array<int, 3>> &tris, real tau,
                                              std::complex<real> z0, real omega) {
  std::vector<vec3> h(static_cast<size_t>(t.n_edges), vec3::Zero());
  for (int f = 0; f < t.n_faces; ++f)
    for (int i = 0; i < 3; ++i) {
      const int e = t.face_edges[static_cast<size_t>(f)][static_cast<size_t>(i)];
      const int a = tris[static_cast<size_t>(f)][static_cast<size_t>(i)];
      const int b = tris[static_cast<size_t>(f)][static_cast<size_t>((i + 1) % 3)];
      const vec3 m = real(0.5) * (x[static_cast<size_t>(a)] + x[static_cast<size_t>(b)]);
      const std::complex<real> z = z0 * std::exp(std::complex<real>(0, omega * m[0]));
      mat3 W = mat3::Zero();
      W(0, 0) = tau / 2 + z.real();
      W(1, 1) = tau / 2 - z.real();
      W(0, 1) = W(1, 0) = z.imag();
      const vec3 tv = (x[static_cast<size_t>(std::max(a, b))] - x[static_cast<size_t>(std::min(a, b))]).normalized();
      const vec3 bv = vec3::UnitZ().cross(tv);
      h[static_cast<size_t>(e)] = vec3(tv.dot(W * tv), tv.dot(W * bv), bv.dot(W * bv));
    }
  return h;
}

inline real face_element_energy(const kusama::stencil_table &t, const std::vector<vec3> &x,
                                const std::vector<vec3> &h, real omega) {
  real E = 0;
  for (int f = 0; f < t.n_faces; ++f) {
    const auto &id = t.ids[static_cast<size_t>(f)];
    const Eigen::Matrix<real, 9, 9> Ef = engman::twisted_spin2_face_element(
        x[static_cast<size_t>(id[0])], x[static_cast<size_t>(id[1])], x[static_cast<size_t>(id[2])],
        vec3::UnitX(), omega);
    Eigen::Matrix<real, 9, 1> hf;
    for (int i = 0; i < 3; ++i)
      hf.segment<3>(3 * i) = h[static_cast<size_t>(t.face_edges[static_cast<size_t>(f)][static_cast<size_t>(i)])];
    E += hf.dot(Ef * hf);
  }
  return E;
}

/// The twisted traceless spin-2 element is the untwisted one at `omega = 0`
/// and vanishes on a field turning at `omega` along its axis (to round-off:
/// the demodulated field is constant on each face); the readout converges to
/// that rate.
GAUDI_TEST(engman_twisted_spin2_and_readout) {
  {
    const vec3 X0(0.1, 0, 0.02), X1(1, 0.2, -0.1), X2(0.3, 0.9, 0.05);
    const real d = (engman::twisted_spin2_face_element(X0, X1, X2, vec3(0.3, 1, 0), 0) -
                    wg::spin2_traceless_face_element(X0, X1, X2))
                       .cwiseAbs()
                       .maxCoeff();
    GAUDI_EXPECT(d < real(1e-13));
  }
  const real omega = 6;
  const std::complex<real> z0(0.8, 0.3);
  real prev = std::numeric_limits<real>::infinity();
  for (const int n : {11, 21, 41}) {
    std::vector<vec3> x;
    std::vector<std::array<int, 3>> tris;
    unit_grid(n, x, tris);
    const kusama::stencil_table t =
        kusama::build_stencil_table(asawa::shell::triangle_walker(static_cast<int>(x.size()), tris));
    const std::vector<vec3> h = planar_twisted_store(t, x, tris, real(0.7), z0, omega);
    const real e_tw = face_element_energy(t, x, h, omega);
    const real e_0 = face_element_energy(t, x, h, 0);
    const real e_neg = face_element_energy(t, x, h, -omega);
    const engman::twist_reading rd =
        engman::read_twist(t, x, h, [](const vec3 &, const vec3 &) { return vec3(vec3::UnitX()); });
    console_logger::info << "twisted spin-2, grid " << n << ": E(omega) / E(0) " << e_tw / e_0
                         << ", E(-omega) / E(0) " << e_neg / e_0 << ", readout " << rd.mean << " of " << omega
                         << std::endl;
    // Exact up to round-off: the CR element is exact on the demodulated field.
    GAUDI_EXPECT(std::abs(e_tw) < real(1e-10) * e_0);
    GAUDI_EXPECT(e_neg > real(2) * e_0);
    GAUDI_EXPECT(std::abs(rd.mean - omega) < prev);
    GAUDI_EXPECT(std::abs(rd.mean - omega) < real(0.02) * omega);
    prev = std::abs(rd.mean - omega);
  }
}

/// The twist-rate energy: its gradient matches central differences, scaling
/// the store barely changes its residual (only the regularization sees the
/// scale), and on a field turning at the target rate the residual is the
/// readout's discretization error plus the regularization bias.
GAUDI_TEST(engman_twist_rate_energy) {
  std::vector<vec3> x;
  std::vector<std::array<int, 3>> tris;
  unit_grid(21, x, tris);
  const kusama::stencil_table t =
      kusama::build_stencil_table(asawa::shell::triangle_walker(static_cast<int>(x.size()), tris));
  const real omega = 6;
  const std::vector<vec3> h = planar_twisted_store(t, x, tris, real(0.7), std::complex<real>(0.8, 0.3), omega);
  std::vector<vec3> hm(h.size(), vec3::Ones()), hv(h.size(), vec3::Zero()), hf = hv, hs = h;
  auto blk = hepworth::vec3_block::create(hm, hs, hv, hf);
  vecX q;
  blk->map_to_x(q);
  const engman::axis_field ax = [](const vec3 &, const vec3 &) { return vec3(vec3::UnitX()); };
  auto e_on = engman::twist_rate::create(t, x, h, ax, omega, 1, blk);
  auto e_off = engman::twist_rate::create(t, x, h, ax, real(0), 1, blk);
  vecX q2 = 2 * q;
  const real on = e_on->residual_sq(q), off = e_off->residual_sq(q), off2 = e_off->residual_sq(q2);
  // The gradient (J^T r) against central differences of |r|^2.
  std::mt19937 rng(5);
  std::normal_distribution<real> nd(0, 1);
  vecX dq(q.size());
  for (index_t i = 0; i < q.size(); ++i)
    dq[i] = nd(rng);
  const real eps = 1e-6;
  const real fd = (e_off->residual_sq(q + eps * dq) - e_off->residual_sq(q - eps * dq)) / (2 * eps);
  vecX g = vecX::Zero(q.size());
  e_off->add_gradient(q, g);
  const real an = 2 * g.dot(dq);
  console_logger::info << "twist rate energy: E(target) / E(0) " << on / off << ", scaled store " << off2 / off
                       << ", d|r|^2 analytic " << an << " fd " << fd << ", rows " << e_on->face_count() << std::endl;
  GAUDI_EXPECT(on < real(1e-3) * off);
  // Scale-free up to the 10% regularization of |z|^2 (a 1% bias at the start).
  GAUDI_EXPECT(std::abs(off2 / off - 1) < real(0.03));
  GAUDI_EXPECT(std::abs(an - fd) < real(1e-6) * std::abs(fd));
}

/// Face area: the Jacobian matches central differences, rigid motions do not
/// change it, and as a hard constraint it is met after one step.
GAUDI_TEST(engman_face_area_constraint) {
  std::vector<vec3> x0;
  std::vector<std::array<int, 3>> tris;
  torus_tris(2.0, 0.7, 16, 8, x0, tris);
  engman_scene sc(x0, tris);
  const std::vector<real> a0 = engman::face_area::areas(tris, x0);
  std::vector<real> grown = a0;
  for (real &a : grown)
    a *= real(1.21);
  const engman::constraint::ptr c = engman::face_area::create(tris, grown, sc.pos);
  vecX q;
  for (auto &b : sc.blocks)
    b->map_to_x(q);
  vecX C;
  std::vector<trip> jt;
  c->evaluate(q, C, &jt);
  wg::matS J(c->rows(), q.size());
  J.setFromTriplets(jt.begin(), jt.end());
  std::mt19937 rng(3);
  std::normal_distribution<real> nd(0, 1);
  vecX dq = vecX::Zero(q.size());
  for (index_t i = 0; i < 3 * static_cast<index_t>(x0.size()); ++i)
    dq[i] = nd(rng);
  const real eps = 1e-6;
  vecX Cp, Cm;
  c->evaluate(q + eps * dq, Cp, nullptr);
  c->evaluate(q - eps * dq, Cm, nullptr);
  const real fd = ((Cp - Cm) / (2 * eps) - J * dq).cwiseAbs().maxCoeff() / (J * dq).cwiseAbs().maxCoeff();
  vecX dr = vecX::Zero(q.size());
  for (size_t i = 0; i < x0.size(); ++i)
    dr.segment<3>(3 * static_cast<index_t>(i)) = vec3(0, 0, 1).cross(x0[i]) + vec3(0.3, -0.2, 0.1);
  const real rigid = (J * dr).cwiseAbs().maxCoeff();

  engman::solver es;
  es.set_energies({hb::dirichlet_memory_gn::create(&sc.t, sc.x, 1, std::vector<hepworth::sim_block::ptr>{sc.pos})});
  es.set_constraints({c});
  es.step(sc.blocks, real(0.1), real(1), 60);
  console_logger::info << "engman face area: J vs FD " << fd << ", rigid " << rigid << ", |C| " << es.last.C_first
                       << " -> " << es.last.C_last << " in " << es.last.iterations << " iterates" << std::endl;
  GAUDI_EXPECT(fd < real(1e-6));
  GAUDI_EXPECT(rigid < real(1e-10));
  GAUDI_EXPECT(es.last.C_last < real(1e-9));
}

/// A constraint as a soft energy approaches the hard one as `1 / w^2`.
GAUDI_TEST(engman_as_energy_approaches_hard) {
  std::vector<vec3> x0;
  std::vector<std::array<int, 3>> tris;
  torus_tris(2.0, 0.7, 16, 8, x0, tris);
  std::vector<index_t> ids;
  std::vector<vec3> goals;
  for (int i = 0; i < 16; ++i) {
    ids.push_back(i * 8);
    goals.push_back(x0[static_cast<size_t>(i * 8)] + vec3(0, 0, 0.1));
  }
  auto pin_error = [&](real w) {
    engman_scene sc(x0, tris);
    engman::solver es;
    es.set_energies({hb::dirichlet_memory_gn::create(&sc.t, sc.x, 1, std::vector<hepworth::sim_block::ptr>{sc.pos})});
    const engman::constraint::ptr pin = engman::position_pin::create(ids, goals, sc.pos);
    if (w > 0) {
      std::vector<hepworth::projection_constraint::ptr> e = {
          hb::dirichlet_memory_gn::create(&sc.t, sc.x, 1, std::vector<hepworth::sim_block::ptr>{sc.pos}),
          engman::as_energy::create(pin, w)};
      es.set_energies(e);
    } else {
      es.set_constraints({pin});
    }
    es.step(sc.blocks, real(0.1), real(1), 60);
    real err = 0;
    for (size_t k = 0; k < ids.size(); ++k)
      err = std::max(err, (sc.x[static_cast<size_t>(ids[k])] - goals[k]).norm());
    return err;
  };
  const real hard = pin_error(0), s2 = pin_error(1e2), s3 = pin_error(1e3);
  console_logger::info << "engman as_energy: pin error hard " << hard << ", soft w 1e2 " << s2 << ", 1e3 " << s3
                       << " (ratio " << s2 / s3 << ")" << std::endl;
  GAUDI_EXPECT(hard < real(1e-10));
  GAUDI_EXPECT(s2 / s3 > real(30));
}

/// A twist-rate scene runs: the coupling holds every frame, nothing goes
/// non-finite, and the mesh's measured twist moves toward the target with the
/// anisotropy kept.
GAUDI_TEST(engman_twist_scene_runs) {
  engman::scene_config cfg = engman::preset_twist_rate();
  cfg.mesh.ring = 32;
  cfg.mesh.tube = 12;
  cfg.log = false;
  engman::scene sc(cfg);
  const real start = sc.reading.mean, z_start = sc.reading.mean_z;
  const real target = engman::find<engman::terms::twist_rate>(cfg.terms)->omega;
  real worst = 0;
  bool finite = true;
  for (int k = 0; k < 6; ++k) {
    sc.step();
    worst = std::max(worst, sc.coupling_rel);
    for (const vec3 &p : sc.x)
      finite = finite && p.allFinite();
  }
  console_logger::info << "engman twist scene: mesh twist " << start << " -> " << sc.reading.mean << " (target "
                       << target << "), |z| " << sc.reading.mean_z << ", worst coupling " << worst
                       << ", last its " << sc.last.iterations << std::endl;
  GAUDI_EXPECT(finite);
  GAUDI_EXPECT(worst < real(1e-6));
  // At h = 0.02 the twist builds over frames: 0.2 of the target by frame 6.
  GAUDI_EXPECT(sc.reading.mean > real(0.1) * target);
  GAUDI_EXPECT(sc.reading.mean_z > real(0.8) * z_start);
}

/// The torsion target touches only the geodesic torsion along its axis: zero
/// on a field with that torsion, still zero after adding any mean curvature
/// or any anisotropy along the axis, and nonzero once the torsion changes.
GAUDI_TEST(engman_torsion_target_only_torsion) {
  std::vector<vec3> x;
  std::vector<std::array<int, 3>> tris;
  unit_grid(11, x, tris);
  const kusama::stencil_table t =
      kusama::build_stencil_table(asawa::shell::triangle_walker(static_cast<int>(x.size()), tris));
  const real tau0 = 0.7;
  const vec3 a = vec3(1, 2, 0).normalized();
  const real phi = std::atan2(a[1], a[0]);
  // z with W(a, b) = tau0: tau0 e^{i (2 phi + pi/2)}; along the axis, z ~ e^{2 i phi} adds no torsion.
  auto field = [&](real tau, real H, real along) {
    return planar_twisted_store(t, x, tris, 2 * H,
                                std::polar(tau, real(2) * phi + real(M_PI / 2)) + std::polar(along, real(2) * phi),
                                real(0));
  };
  const engman::axis_field ax = [a](const vec3 &, const vec3 &) { return a; };
  auto energy = [&](const std::vector<vec3> &h) {
    std::vector<vec3> hm(h.size(), vec3::Ones()), hv(h.size(), vec3::Zero()), hf = hv, hs = h;
    auto blk = hepworth::vec3_block::create(hm, hs, hv, hf);
    vecX q;
    blk->map_to_x(q);
    auto e = engman::torsion_target::create(t, x, ax, tau0, 1, blk);
    index_t id0 = 0;
    std::vector<trip> trips;
    e->fill_A(id0, trips);
    wg::matS A(id0, q.size());
    A.setFromTriplets(trips.begin(), trips.end());
    vecX p = vecX::Zero(id0);
    e->project(q, p);
    return (A * q - p).squaredNorm();
  };
  const real e_on = energy(field(tau0, 0, 0)), e_H = energy(field(tau0, 1.3, 0)),
             e_along = energy(field(tau0, 0, 0.9)), e_off = energy(field(0.2, 0, 0));
  const engman::torsion_reading tr = engman::read_torsion(t, x, field(tau0, 1.3, 0.9), ax);
  console_logger::info << "torsion target: on " << e_on << ", + mean curvature " << e_H << ", + anisotropy along "
                       << e_along << ", other torsion " << e_off << ", readout " << tr.mean << std::endl;
  GAUDI_EXPECT(e_on < real(1e-24));
  GAUDI_EXPECT(e_H < real(1e-24));
  GAUDI_EXPECT(e_along < real(1e-24));
  GAUDI_EXPECT(e_off > real(1e-3));
  GAUDI_EXPECT(std::abs(tr.mean - tau0) < real(1e-12));
}

/// The torsion scene twists the mesh itself: geodesic torsion of the
/// read-out grows with the sign of the target, the coupling holds.
GAUDI_TEST(engman_torsion_scene_twists) {
  for (const real sgn : {real(1), real(-1)}) {
    engman::scene_config cfg = engman::preset_torsion();
    cfg.mesh.ring = 32;
    cfg.mesh.tube = 12;
    cfg.log = false;
    engman::terms::torsion *tq = engman::find<engman::terms::torsion>(cfg.terms);
    tq->tau *= sgn;
    engman::scene sc(cfg);
    const real start = sc.torsion.mean;
    real worst = 0;
    for (int k = 0; k < 10; ++k) {
      sc.step();
      worst = std::max(worst, sc.coupling_rel);
    }
    console_logger::info << "engman torsion scene tau " << tq->tau << ": torsion " << start << " -> "
                         << sc.torsion.mean << " (sin2d " << sc.torsion.relative << "), displacement/r "
                         << sc.displacement() << ", worst coupling " << worst << ", last its " << sc.last.iterations
                         << std::endl;
    GAUDI_EXPECT(worst < real(1e-6));
    // Grows from the noise (an instability): about 0.05 by frame 10 on the demo torus.
    GAUDI_EXPECT(sgn * sc.torsion.mean > sgn * start + real(0.01));
  }
}

/// Exploration, run only with `GAUDI_SWEEP` set: one scene over 30 frames.
/// `GAUDI_PRESET` names the preset; `GAUDI_TAU`, `GAUDI_TORSION_W`, `GAUDI_BEND_W`,
/// `GAUDI_OMEGA`, `GAUDI_REFRESH`, `GAUDI_ITS`, `GAUDI_H`, `GAUDI_RHO_REL`,
/// `GAUDI_FRAMES` override knobs.
GAUDI_TEST(engman_sweep) {
  if (!std::getenv("GAUDI_SWEEP"))
    return;
  auto env = [](const char *k) { return std::getenv(k); };
  engman::scene_config cfg = engman::preset(env("GAUDI_PRESET") ? env("GAUDI_PRESET") : "torsion");
  cfg.log = env("GAUDI_LOG") != nullptr;
  if (auto *tq = engman::find<engman::terms::torsion>(cfg.terms)) {
    if (env("GAUDI_TAU"))
      tq->tau = real(std::atof(env("GAUDI_TAU")));
    if (env("GAUDI_TORSION_W"))
      tq->w = real(std::atof(env("GAUDI_TORSION_W")));
  }
  if (auto *bd = engman::find<engman::terms::bending>(cfg.terms))
    if (env("GAUDI_BEND_W"))
      bd->w = real(std::atof(env("GAUDI_BEND_W")));
  if (auto *tw = engman::find<engman::terms::twist_rate>(cfg.terms))
    if (env("GAUDI_OMEGA"))
      tw->omega = real(std::atof(env("GAUDI_OMEGA")));
  if (env("GAUDI_REFRESH"))
    cfg.step.solver.refresh_every = std::atoi(env("GAUDI_REFRESH"));
  if (env("GAUDI_ITS"))
    cfg.step.iterations = std::atoi(env("GAUDI_ITS"));
  if (env("GAUDI_H"))
    cfg.step.h = real(std::atof(env("GAUDI_H")));
  if (env("GAUDI_RHO_REL"))
    cfg.step.solver.rho_rel = real(std::atof(env("GAUDI_RHO_REL")));
  if (env("GAUDI_REMESH"))
    cfg.remesh.enabled = std::atoi(env("GAUDI_REMESH")) != 0;
  const int frames = env("GAUDI_FRAMES") ? std::atoi(env("GAUDI_FRAMES")) : 30;
  engman::scene sc(cfg);
  std::ostringstream os;
  const auto clock0 = std::chrono::steady_clock::now();
  for (int k = 0; k < frames; ++k) {
    sc.step();
    if (k % 6 == 5)
      os << " [" << k + 1 << "] twist " << sc.reading.mean << " torsion " << sc.torsion.mean << " sin2d "
         << sc.torsion.relative << " coupling " << sc.coupling_rel << " disp/r " << sc.displacement() << " verts "
         << sc.x.size();
  }
  std::ostringstream terms;
  for (const engman::term &t : cfg.terms)
    terms << engman::name(t) << ";";
  console_logger::info << "sweep " << terms.str() << os.str() << " |z| " << sc.reading.mean_z << " its "
                       << sc.last.iterations << (sc.last.stopped ? " stopped" : "") << " s/frame "
                       << std::chrono::duration<double>(std::chrono::steady_clock::now() - clock0).count() / frames
                       << std::endl;
}

/// Exploration, run only with `GAUDI_REMESH_PROBE` set: the remesher alone on
/// the startup torus, vertex count and edge lengths per pass.
GAUDI_TEST(engman_remesh_probe) {
  if (!std::getenv("GAUDI_REMESH_PROBE"))
    return;
  engman::scene_config cfg = engman::preset_torsion();
  auto M = asawa::shell::load_torus(cfg.mesh.R, cfg.mesh.r, cfg.mesh.ring, cfg.mesh.tube);
  asawa::shell::triangulate(*M);
  auto stats = [&](const char *tag) {
    const std::vector<vec3> &x = asawa::get_vec_data(*M, 0);
    real mn = 1e9, mx = 0, sum = 0;
    int n = 0;
    for (auto c : M->get_edge_range()) {
      const real l = asawa::shell::edge_length(*M, c, x);
      mn = std::min(mn, l);
      mx = std::max(mx, l);
      sum += l;
      ++n;
    }
    console_logger::info << tag << ": verts " << M->vert_count() << " faces " << M->face_count() << " edges " << n
                         << " length min " << mn << " mean " << sum / std::max(n, 1) << " max " << mx << std::endl;
    return sum / std::max(n, 1);
  };
  const real l0 = stats("start");
  console_logger::info << "avg_length " << asawa::shell::avg_length(*M, asawa::get_vec_data(*M, 0)) << std::endl;
  auto surf = asawa::shell::dynamic::create(M, cfg.remesh.collapse * l0, cfg.remesh.split * l0, cfg.remesh.bridge * l0);
  for (int k = 0; k < 4; ++k) {
    surf->step(true);
    stats("pass");
  }
}

} // namespace
} // namespace test
} // namespace gaudi

#endif
