#ifndef GAUDI_ENGMAN_SCENE_HPP
#define GAUDI_ENGMAN_SCENE_HPP

/// A torus on the Engman solver, built from a list of terms
/// (`engman/terms.hpp`). The scene owns the mesh, the edge store, the rest
/// positions and the solver; everything that shapes the motion is a term in
/// `scene_config::terms`, so trying a combination is editing one list:
///
///     scene_config cfg{.terms = {terms::coupling{},
///                                terms::memory{},
///                                terms::torsion{.tau = 1},
///                                terms::face_area{.m = mode::soft}}};
///
/// Each frame restarts: the store reads the live mesh, the rest relaxes
/// toward it (`rest_config`), every term builds, one Engman step runs, and the
/// surface is remeshed (`remesh_config`, the `asawa::shell::dynamic` remesher
/// set up as in `monge_fairing_demo`). Positions, velocities, rest and startup
/// positions are vertex datums, so splits interpolate them and collapses merge
/// them.

#include "gaudi/asawa/datums.hpp"
#include "gaudi/asawa/primitive_objects.hpp"
#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/asawa/shell/dynamic.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/asawa/shell/stencil_walker.hpp"
#include "gaudi/engman/constraint.hpp"
#include "gaudi/engman/energies.hpp"
#include "gaudi/engman/solver.hpp"
#include "gaudi/engman/terms.hpp"
#include "gaudi/hepworth/block/sim_block.hpp"
#include "gaudi/kusama/stencil_table.hpp"
#include "gaudi/kusama/weingarten.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <string>
#include <vector>

namespace gaudi {
namespace engman {

/// The torus (axis +Z) and its startup noise along the normals, in mean edge
/// lengths (symmetry breaking).
struct mesh_config {
  real R = 1.0;
  real r = 0.4;
  int ring = 48;
  int tube = 16;
  real noise = 0.02;
  unsigned seed = 1;
};

struct step_config {
  /// Frame time step. 0.05 lets a strong twist move the shape too far in one
  /// frame for the coupling to converge; 0.02 keeps it exact.
  real h = 0.02;
  /// Velocity kept per frame is `1 - damping`.
  real damping = 0.5;
  int iterations = 20;
  solver_config solver{};
};

/// Rest positions for the metric terms (`face_area`, `conformal`): they relax
/// toward the live shape by `relax` per frame. 1: rest at the frame start, a
/// viscous sheet. 0: the startup shape, an elastic one.
struct rest_config {
  real relax = 1.0;
};

/// The remesher (`asawa::shell::dynamic`, as in `monge_fairing_demo`):
/// collapse edges shorter than `collapse l0`, split edges longer than
/// `split l0`, bridge within `bridge l0`, `l0` the startup mean edge length.
/// The band must contain the startup edges or every pass splits (or
/// collapses) most of them: the torus's edges span about 0.5 to 1.6 `l0`.
/// `passes` remesher steps per frame after the solve.
struct remesh_config {
  bool enabled = true;
  int every = 1;
  int passes = 1;
  real collapse = 0.5;
  real split = 1.6;
  real bridge = 0.25;
  asawa::shell::qc_config qc{};
};

/// The axis the twist and torsion readouts measure along.
struct readout_config {
  axis_spec axis{};
};

enum class overlay_kind { plain, twist_rate, anisotropy, mean_curvature };

struct draw_config {
  overlay_kind overlay = overlay_kind::twist_rate;
  /// Vertex frames as in `monge_fairing_demo`: red along k0 (the larger
  /// principal curvature), green along k1, blue the normal.
  bool frames = true;
  /// Scale red and green by `|k| / mean |k|` instead of a fixed length.
  bool frames_curvature = false;
  /// Frame half-length in mean edge lengths.
  real frame_length = 0.5;
  /// Principal direction (larger |kappa|) per face, as a short black line.
  bool principal_lines = false;
  real line_scale = 0.6;
  /// Color range of `twist_rate` (rad per length) and `anisotropy` (|z|).
  real twist_range = 2.0;
  real anisotropy_range = 2.0;
};

struct scene_config {
  mesh_config mesh{};
  step_config step{};
  rest_config rest{};
  remesh_config remesh{};
  std::vector<term> terms{};
  readout_config readout{};
  draw_config draw{};
  /// Log one line per frame (`[engman] ...`).
  bool log = true;
};

// ----------------------------------------------------------------- presets

/// A torsion target along the ring with plate bending at its own weight:
/// only the twist is asked for, so curvature changes and the surface can
/// buckle to carry it (the default).
inline scene_config preset_torsion() {
  return {.terms = {terms::coupling{}, terms::memory{}, terms::bending{.w = 0.1},
                    terms::torsion{.w = 10, .tau = 4}}};
}

/// The same along the (1, 1) winding, near the Hopf (Villarceau) fibers.
inline scene_config preset_torsion_hopf() {
  const axis_spec hopf{.kind = axis_kind::winding};
  return {.terms = {terms::coupling{}, terms::memory{}, terms::bending{.w = 0.1},
                    terms::torsion{.w = 10, .tau = 4, .axis = hopf}},
          .readout = {.axis = hopf}};
}

/// Torsion with the metric held: hard face areas against the startup shape
/// and a shear modulus, so the twist has to be bought by bending.
inline scene_config preset_torsion_elastic() {
  return {.rest = {.relax = 0.0},
          .terms = {terms::coupling{}, terms::memory{}, terms::bending{.w = 0.1},
                    terms::torsion{.w = 10, .tau = 4}, terms::face_area{}, terms::conformal{.mu = 1}}};
}

/// A twist-rate target around the ring on the store, plate bending kept.
inline scene_config preset_twist_rate() {
  return {.terms = {terms::coupling{}, terms::memory{}, terms::bending{.w = 1e-3}, terms::twist_rate{.w = 10}}};
}

/// Willmore-type relaxation of a noisy torus, no twist.
inline scene_config preset_relax() {
  return {.terms = {terms::coupling{}, terms::memory{}, terms::bending{.w = 1}}};
}

/// A preset by name (`torsion`, `torsion_hopf`, `torsion_elastic`,
/// `twist_rate`, `relax`); unknown names give `torsion`.
inline scene_config preset(const std::string &name) {
  if (name == "torsion_hopf")
    return preset_torsion_hopf();
  if (name == "torsion_elastic")
    return preset_torsion_elastic();
  if (name == "twist_rate")
    return preset_twist_rate();
  if (name == "relax")
    return preset_relax();
  return preset_torsion();
}

// ------------------------------------------------------------------- scene

class scene {
public:
  using ptr = std::shared_ptr<scene>;
  static ptr create(const scene_config &cfg = {}) { return std::make_shared<scene>(cfg); }

  explicit scene(const scene_config &cfg = {}) : cfg(cfg) { reset(); }

  scene_config cfg;

  // Live state (read by the demo for drawing).
  std::vector<vec3> x;
  std::vector<std::array<int, 3>> tris;
  kusama::stencil_table table;
  std::vector<vec3> hs; ///< the edge store
  /// The twist of the mesh (its own CR read-out, not the store).
  twist_reading reading;
  /// Geodesic torsion of the mesh along the twist axis.
  torsion_reading torsion;
  /// Max coupling residual `|R(x) - S|` over edges relative to the max `|R(x)|`.
  real coupling_rel = 0;
  int frame = 0;
  solver_report last;
  real mean_edge = 1;

  void reset() {
    M = asawa::shell::load_torus(cfg.mesh.R, cfg.mesh.r, std::max(cfg.mesh.ring, 3),
                                 std::max(cfg.mesh.tube, 3));
    asawa::shell::triangulate(*M);
    std::vector<vec3> &x0 = asawa::get_vec_data(*M, 0);
    _slot_v = asawa::init_vert_datum<vec3>(*M, vec3::Zero());
    _slot_xbar = asawa::init_vert_datum<vec3>(*M, vec3::Zero());
    _slot_start = asawa::init_vert_datum<vec3>(*M, vec3::Zero());
    // The true mean edge length (`asawa::shell::avg_length` returns about a
    // quarter of it on this torus, which would put every edge above `split`).
    {
      real sum = 0;
      int n = 0;
      for (auto c : M->get_edge_range()) {
        sum += asawa::shell::edge_length(*M, c, x0);
        ++n;
      }
      _l0 = n > 0 ? sum / n : real(1);
    }
    _surf = asawa::shell::dynamic::create(M, cfg.remesh.collapse * _l0, cfg.remesh.split * _l0,
                                          cfg.remesh.bridge * _l0);
    _surf->set_qc(cfg.remesh.qc);
    rebuild();
    add_noise(cfg.mesh.noise * mean_edge, cfg.mesh.seed);
    xbar = x;
    x_start = x;
    to_shell();
    rebuild();
    frame = 0;
    measure();
    if (cfg.log) {
      std::cout << "[engman] terms:";
      for (const term &t : cfg.terms)
        std::cout << " " << name(t) << ";";
      std::cout << std::endl;
    }
  }

  /// The shell: positions in datum 0, velocities, rest and startup positions
  /// in their own vertex datums.
  asawa::shell::shell::ptr M;

  /// Max displacement from the startup shape over the tube radius.
  real displacement() const {
    real d = 0;
    for (size_t i = 0; i < x.size() && i < x_start.size(); ++i)
      d = std::max(d, (x[i] - x_start[i]).norm());
    return d / std::max(cfg.mesh.r, real(1e-12));
  }

  /// Twist and coupling residual of the live mesh.
  void measure() {
    const std::vector<vec3> R = kusama::weingarten::edge_readings(table, x);
    reading = read_twist(table, x, R, axis());
    torsion = read_torsion(table, x, R, axis());
    real num = 0, den = 0;
    for (size_t e = 0; e < R.size(); ++e) {
      if (!R[e].allFinite() || !hs[e].allFinite())
        continue;
      num = std::max(num, (R[e] - hs[e]).cwiseAbs().maxCoeff());
      den = std::max(den, R[e].cwiseAbs().maxCoeff());
    }
    coupling_rel = den > real(0) ? num / den : real(0);
  }

  /// The readout axis (`cfg.readout.axis`).
  axis_field axis() const { return make_axis(cfg.readout.axis, cfg.mesh.r); }

  void step() {
    const auto clock0 = std::chrono::steady_clock::now();
    // Frame start: the store reads the live mesh (the coupling is idle), the
    // rest relaxes toward it.
    const real eta = std::clamp(cfg.rest.relax, real(0), real(1));
    for (size_t i = 0; i < x.size(); ++i)
      xbar[i] += eta * (x[i] - xbar[i]);
    hs = kusama::weingarten::edge_readings(table, x);
    for (vec3 &h : hs)
      if (!h.allFinite())
        h.setZero();
    std::fill(hv.begin(), hv.end(), vec3::Zero());
    std::fill(hf.begin(), hf.end(), vec3::Zero());

    std::vector<hepworth::projection_constraint::ptr> energies;
    std::vector<constraint::ptr> constraints;
    build(energies, constraints);

    _solver.cfg = cfg.step.solver;
    _solver.set_energies(energies);
    _solver.set_constraints(constraints);
    _solver.step(blocks, cfg.step.h, cfg.step.damping, cfg.step.iterations);
    last = _solver.last;
    // Readouts on the solved mesh, before the remesh (the curvature survives
    // the remesh with the geometry; the next frame re-reads it).
    measure();
    to_shell();
    if (cfg.remesh.enabled && cfg.remesh.every > 0 && (frame + 1) % cfg.remesh.every == 0) {
      const auto r0 = std::chrono::steady_clock::now();
      for (int k = 0; k < cfg.remesh.passes; ++k)
        _surf->step(true);
      const auto r1 = std::chrono::steady_clock::now();
      rebuild();
      if (cfg.log)
        std::cout << "[engman] remesh " << std::chrono::duration<double>(r1 - r0).count() << " s, rebuild "
                  << std::chrono::duration<double>(std::chrono::steady_clock::now() - r1).count()
                  << " s, verts " << x.size() << std::endl;
    }
    if (cfg.log)
      std::cout << "[engman] frame " << frame << " verts " << x.size() << " s "
                << std::chrono::duration<double>(std::chrono::steady_clock::now() - clock0).count()
                << " torsion " << torsion.mean << " (sin2d " << torsion.relative << ") twist rate "
                << reading.mean << " coupling " << coupling_rel << " |z| " << reading.mean_z << " its "
                << last.iterations << " restore " << last.restorations << " |C| " << last.C_first << " -> "
                << last.C_last << " E " << last.E_first << " -> " << last.E_last
                << (last.stopped ? " (stopped)" : "") << std::endl;
    ++frame;
  }

  /// Per-vertex scalar of `cfg.draw.overlay`, mapped to a color.
  std::vector<vec3> colors() const {
    const size_t n = x.size();
    std::vector<vec3> c(n, vec3(0.62, 0.66, 0.72));
    if (cfg.draw.overlay == overlay_kind::plain)
      return c;
    std::vector<real> face_val(tris.size(), std::numeric_limits<real>::quiet_NaN());
    if (cfg.draw.overlay == overlay_kind::twist_rate) {
      face_val = reading.omega;
    } else {
      for (int fi = 0; fi < table.n_faces; ++fi) {
        vec3 hsum = vec3::Zero();
        for (int i = 0; i < 3; ++i)
          hsum +=
              hs[static_cast<size_t>(table.face_edges[static_cast<size_t>(fi)][static_cast<size_t>(i)])] / 3;
        // Edge frames differ per edge, but |z| and the trace do not depend on the frame.
        real val = 0;
        if (cfg.draw.overlay == overlay_kind::anisotropy) {
          real zs = 0;
          for (int i = 0; i < 3; ++i) {
            const vec3 &h =
                hs[static_cast<size_t>(table.face_edges[static_cast<size_t>(fi)][static_cast<size_t>(i)])];
            zs += std::hypot((h[0] - h[2]) / 2, h[1]) / 3;
          }
          val = zs;
        } else {
          val = (hsum[0] + hsum[2]) / 2;
        }
        face_val[static_cast<size_t>(fi)] = val;
      }
    }
    std::vector<real> acc(n, 0), wsum(n, 0);
    for (size_t fi = 0; fi < tris.size(); ++fi) {
      if (!std::isfinite(face_val[fi]))
        continue;
      for (int k = 0; k < 3; ++k) {
        acc[static_cast<size_t>(tris[fi][static_cast<size_t>(k)])] += face_val[fi];
        wsum[static_cast<size_t>(tris[fi][static_cast<size_t>(k)])] += 1;
      }
    }
    const real range = cfg.draw.overlay == overlay_kind::twist_rate
                           ? cfg.draw.twist_range
                           : (cfg.draw.overlay == overlay_kind::anisotropy ? cfg.draw.anisotropy_range
                                                                           : real(2) / cfg.mesh.r);
    for (size_t i = 0; i < n; ++i) {
      if (!(wsum[i] > 0))
        continue;
      const real s = std::clamp(acc[i] / wsum[i] / std::max(range, real(1e-12)), real(-1), real(1));
      // Diverging: blue (negative), grey (zero), red (positive).
      const vec3 grey(0.62, 0.66, 0.72), pos_c(0.85, 0.25, 0.2), neg_c(0.2, 0.4, 0.85);
      c[i] = s >= 0 ? vec3((1 - s) * grey + s * pos_c) : vec3((1 + s) * grey - s * neg_c);
    }
    return c;
  }

  /// One segment per face along its principal direction of larger |kappa|.
  /// One frame per vertex: `(position, k0 axis, k1 axis, normal, k0, k1)`.
  struct vertex_frame {
    vec3 x, d0, d1, n;
    real k0 = 0, k1 = 0;
  };
  void vertex_frames(std::vector<vertex_frame> &out) const {
    out.clear();
    std::vector<vec3> N(x.size(), vec3::Zero());
    for (const auto &t : tris) {
      const vec3 n = (x[static_cast<size_t>(t[1])] - x[static_cast<size_t>(t[0])])
                         .cross(x[static_cast<size_t>(t[2])] - x[static_cast<size_t>(t[0])]);
      for (int k = 0; k < 3; ++k)
        N[static_cast<size_t>(t[static_cast<size_t>(k)])] += n;
    }
    for (int v = 0; v < static_cast<int>(x.size()); ++v) {
      vec3 n = N[static_cast<size_t>(v)];
      if (!(n.norm() > real(0)))
        continue;
      n.normalize();
      vec3 T = n.unitOrthogonal();
      const vec3 B = n.cross(T);
      const vec3 S = kusama::weingarten::vertex_tensor(
          _stars, v, [&](int i) { return x[static_cast<size_t>(i)]; },
          [&](int e) { return hs[static_cast<size_t>(e)]; }, T, B, n);
      if (!S.allFinite())
        continue;
      Eigen::Matrix<real, 2, 2> M;
      M << S[0], S[1], S[1], S[2];
      Eigen::SelfAdjointEigenSolver<Eigen::Matrix<real, 2, 2>> es(M);
      // Eigen sorts ascending: k0 (larger) is index 1.
      out.push_back({.x = x[static_cast<size_t>(v)],
                     .d0 = es.eigenvectors()(0, 1) * T + es.eigenvectors()(1, 1) * B,
                     .d1 = es.eigenvectors()(0, 0) * T + es.eigenvectors()(1, 0) * B,
                     .n = n,
                     .k0 = es.eigenvalues()[1],
                     .k1 = es.eigenvalues()[0]});
    }
  }

  void principal_segments(std::vector<std::array<vec3, 2>> &out) const {
    namespace wg = kusama::weingarten;
    out.clear();
    const real len = cfg.draw.line_scale * mean_edge;
    for (int fi = 0; fi < table.n_faces; ++fi) {
      const auto &id = table.ids[static_cast<size_t>(fi)];
      const vec3 X[3] = {x[static_cast<size_t>(id[0])], x[static_cast<size_t>(id[1])],
                         x[static_cast<size_t>(id[2])]};
      vec3 T, B, N;
      detail::face_frame(X[0], X[1], X[2], T, B, N);
      if (!N.allFinite())
        continue;
      vec3 hf = vec3::Zero();
      for (int i = 0; i < 3; ++i) {
        const int e = table.face_edges[static_cast<size_t>(fi)][static_cast<size_t>(i)];
        const vec3 d = (X[(i + 1) % 3] - X[i]).normalized();
        hf += wg::tensor_map(d.dot(T), d.dot(B)) * hs[static_cast<size_t>(e)] / real(3);
      }
      Eigen::Matrix<real, 2, 2> S;
      S << hf[0], hf[1], hf[1], hf[2];
      Eigen::SelfAdjointEigenSolver<Eigen::Matrix<real, 2, 2>> es(S);
      const int k = std::abs(es.eigenvalues()[0]) > std::abs(es.eigenvalues()[1]) ? 0 : 1;
      const vec3 dir = es.eigenvectors()(0, k) * T + es.eigenvectors()(1, k) * B;
      const vec3 c = (X[0] + X[1] + X[2]) / real(3);
      out.push_back({c - real(0.5) * len * dir, c + real(0.5) * len * dir});
    }
  }

private:
  asawa::shell::dynamic::ptr _surf;
  index_t _slot_v = -1, _slot_xbar = -1, _slot_start = -1;
  real _l0 = 1;

  /// Copies the live state into the shell's datums (before a remesh).
  void to_shell() {
    asawa::get_vec_data(*M, 0) = x;
    asawa::get_vec_data(*M, _slot_v) = v;
    asawa::get_vec_data(*M, _slot_xbar) = xbar;
    asawa::get_vec_data(*M, _slot_start) = x_start;
  }

  /// Rebuilds everything topology-dependent from the shell: the state
  /// vectors, the stencil table and triangles, masses, the store, blocks.
  void rebuild() {
    x = asawa::get_vec_data(*M, 0);
    const size_t n = x.size();
    auto take = [&](index_t slot, std::vector<vec3> &out, const std::vector<vec3> &fallback) {
      const std::vector<vec3> &d = asawa::get_vec_data(*M, slot);
      out = d.size() == n ? d : fallback;
    };
    take(_slot_v, v, std::vector<vec3>(n, vec3::Zero()));
    take(_slot_xbar, xbar, x);
    take(_slot_start, x_start, x);
    table = kusama::build_stencil_table(asawa::shell::shell_walker(*M));
    tris.assign(static_cast<size_t>(table.n_faces), {0, 0, 0});
    for (int fi = 0; fi < table.n_faces; ++fi)
      for (int k = 0; k < 3; ++k)
        tris[static_cast<size_t>(fi)][static_cast<size_t>(k)] =
            table.ids[static_cast<size_t>(fi)][static_cast<size_t>(k)];
    mean_edge = mean_edge_length();
    _stars = kusama::weingarten::build_vertex_stars(table, x);
    f.assign(n, vec3::Zero());
    m.assign(n, vec3::Ones());
    {
      const std::vector<real> A = vertex_areas();
      real mean = 0;
      for (real a : A)
        mean += a;
      mean /= std::max<real>(real(n), 1);
      for (size_t i = 0; i < n; ++i)
        m[i] = vec3::Constant(mean > real(0) ? std::max(A[i] / mean, real(1e-3)) : real(1));
    }
    hs = kusama::weingarten::edge_readings(table, x);
    for (vec3 &h : hs)
      if (!h.allFinite())
        h.setZero();
    hv.assign(hs.size(), vec3::Zero());
    hf = hv;
    hm.assign(hs.size(), vec3::Constant(real(1e-8)));
    pos = hepworth::vec3_block::create(m, x, v, f);
    store = hepworth::vec3_block::create(hm, hs, hv, hf);
    blocks = {pos, store};
  }

  std::vector<vec3> v, f, m, hv, hf, hm;
  std::vector<vec3> xbar, x_start;
  kusama::weingarten::vertex_stars _stars;
  hepworth::sim_block::ptr pos, store;
  std::vector<hepworth::sim_block::ptr> blocks;
  solver _solver;

  void build(std::vector<hepworth::projection_constraint::ptr> &energies,
             std::vector<constraint::ptr> &constraints) {
    build_context c{.table = table,
                    .tris = tris,
                    .x = x,
                    .xbar = xbar,
                    .hs = hs,
                    .pos = pos,
                    .store = store,
                    .blocks = blocks,
                    .frame = frame,
                    .tube_r = cfg.mesh.r,
                    .energies = energies,
                    .constraints = constraints};
    for (const term &t : cfg.terms)
      engman::build(t, c);
  }

  real mean_edge_length() const {
    real s = 0;
    int n = 0;
    for (const auto &t : tris)
      for (int k = 0; k < 3; ++k) {
        s += (x[static_cast<size_t>(t[static_cast<size_t>((k + 1) % 3)])] -
              x[static_cast<size_t>(t[static_cast<size_t>(k)])])
                 .norm();
        ++n;
      }
    return n > 0 ? s / n : real(1);
  }

  std::vector<real> vertex_areas() const {
    std::vector<real> A(x.size(), real(0));
    for (const auto &t : tris) {
      const real a = detail::face_area(x[static_cast<size_t>(t[0])], x[static_cast<size_t>(t[1])],
                                       x[static_cast<size_t>(t[2])]);
      for (int k = 0; k < 3; ++k)
        A[static_cast<size_t>(t[static_cast<size_t>(k)])] += a / real(3);
    }
    return A;
  }

  void add_noise(real amplitude, unsigned seed) {
    if (!(amplitude > real(0)))
      return;
    std::vector<vec3> N(x.size(), vec3::Zero());
    for (const auto &t : tris) {
      const vec3 n = (x[static_cast<size_t>(t[1])] - x[static_cast<size_t>(t[0])])
                         .cross(x[static_cast<size_t>(t[2])] - x[static_cast<size_t>(t[0])]);
      for (int k = 0; k < 3; ++k)
        N[static_cast<size_t>(t[static_cast<size_t>(k)])] += n;
    }
    std::mt19937 rng(seed);
    std::uniform_real_distribution<real> u(real(-1), real(1));
    for (size_t i = 0; i < x.size(); ++i)
      if (N[i].norm() > real(0))
        x[i] += amplitude * u(rng) * N[i].normalized();
  }
};

} // namespace engman
} // namespace gaudi

#endif
