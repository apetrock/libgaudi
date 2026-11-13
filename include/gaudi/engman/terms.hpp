#ifndef GAUDI_ENGMAN_TERMS_HPP
#define GAUDI_ENGMAN_TERMS_HPP

/// Composable terms for an Engman scene (`engman/scene.hpp`). A scene is a
/// mesh plus a list of terms; each term is a small struct with its own knobs
/// and a `build` that adds its energies or constraints for one frame. Terms
/// can repeat (two `torsion` terms on different axes) and appear in any order.
///
/// Constraint terms carry a `mode` (`off | soft | hard`): hard is a KKT row
/// with a multiplier, soft the same residual as a weighted energy
/// (`as_energy`). Energies are always soft.
///
/// To add a term: write a struct with `name()` and `build(build_context &)`
/// in `namespace terms`, then add it to the `term` variant at the bottom.

#include "gaudi/engman/constraint.hpp"
#include "gaudi/engman/energies.hpp"
#include "gaudi/engman/store_metric.hpp"
#include "gaudi/hepworth/block/shell_constraints.hpp"
#include "gaudi/hepworth/block/sim_block.hpp"
#include "gaudi/hepworth/block/weingarten_gn.hpp"
#include "gaudi/kusama/stencil_table.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <variant>
#include <vector>

namespace gaudi {
namespace engman {

/// How a constraint term enters the solve.
enum class mode {
  off,  ///< not built
  soft, ///< a weighted Gauss-Newton energy (`as_energy`), a compromise
  hard, ///< a KKT row with a multiplier, met
};

inline const char *mode_name(mode m) { return m == mode::off ? "off" : (m == mode::soft ? "soft" : "hard"); }

// --------------------------------------------------------------- axis

/// A tangent direction field on a torus about +Z. `toroidal` runs around the
/// ring, `poloidal` around the tube, `winding` along the `(ring, tube)`
/// torus-knot direction `ring rho e_phi + tube r e_theta` (`(1, 1)` is close
/// to the Hopf / Villarceau family).
enum class axis_kind { toroidal, poloidal, winding };

struct axis_spec {
  axis_kind kind = axis_kind::toroidal;
  int ring = 1;
  int tube = 1;
};

/// The field of `a` on a torus of tube radius `r`.
inline axis_field make_axis(const axis_spec &a, real r) {
  const axis_kind kind = a.kind;
  const real p = real(a.ring), q = real(a.tube);
  return [=](const vec3 &c, const vec3 &N) -> vec3 {
    vec3 ephi = vec3::UnitZ().cross(c);
    const real rho = ephi.norm();
    if (!(rho > real(0)))
      return vec3::Zero();
    ephi /= rho;
    const vec3 etheta = N.cross(ephi);
    if (kind == axis_kind::toroidal)
      return ephi;
    if (kind == axis_kind::poloidal)
      return etheta;
    return p * rho * ephi + q * r * etheta;
  };
}

/// `start + rate * frame`, clamped to `|limit|`.
inline real ramp(real start, real rate, real limit, int frame) {
  return std::clamp(start + rate * real(frame), -std::abs(limit), std::abs(limit));
}

// ------------------------------------------------------------ context

/// What a term sees when it builds: the frame-start geometry and store, the
/// rest positions (`xbar`, see `rest_config`), the blocks, and the lists to
/// add to.
struct build_context {
  const kusama::stencil_table &table;
  const std::vector<std::array<int, 3>> &tris;
  const std::vector<vec3> &x;
  const std::vector<vec3> &xbar;
  const std::vector<vec3> &hs;
  hepworth::sim_block::ptr pos, store;
  const std::vector<hepworth::sim_block::ptr> &blocks;
  int frame = 0;
  /// The torus tube radius (for `winding` axes).
  real tube_r = 1;
  std::vector<hepworth::projection_constraint::ptr> &energies;
  std::vector<constraint::ptr> &constraints;

  void add(const constraint::ptr &c, mode m, real w_soft) {
    if (m == mode::hard)
      constraints.push_back(c);
    else if (m == mode::soft)
      energies.push_back(as_energy::create(c, w_soft));
  }
  axis_field axis(const axis_spec &a) const { return make_axis(a, tube_r); }
};

// -------------------------------------------------------------- terms

namespace terms {

/// The position-store coupling `R(x) = S`: the store is the mesh's shape
/// operator. Soft at `w_soft` is the `projection_solver`-style penalty.
struct coupling {
  mode m = mode::hard;
  real w_soft = 1e2;
  std::string name() const { return std::string("coupling ") + mode_name(m); }
  void build(build_context &c) const {
    if (m != mode::off)
      c.add(engman::coupling::create(&c.table, c.x, c.blocks), m, w_soft);
  }
};

/// The intrinsic vector Dirichlet memory on positions, rest at the frame
/// start: a proximal term that charges rough metric change.
struct memory {
  real w = 1;
  std::string name() const { return "memory"; }
  void build(build_context &c) const {
    if (w > 0)
      c.energies.push_back(hepworth::block::dirichlet_memory_gn::create(
          &c.table, c.x, w, std::vector<hepworth::sim_block::ptr>{c.pos}));
  }
};

/// The store's split spin-2 metric (`store_metric`): shapes steps, leaves
/// the solution.
struct store_metric {
  real lambda_H = 0;
  real lambda_z = 0;
  std::string name() const { return "store_metric"; }
  void build(build_context &c) const {
    for (auto &e : engman::store_metric(&c.table, c.x, c.store, lambda_H, lambda_z))
      c.energies.push_back(e);
  }
};

/// Bending on the store: `plate` `w int |W|^2` (Willmore-type), `helfrich`
/// `w int (2H - c0)^2`.
struct bending {
  bending_mode mode = bending_mode::plate;
  real w = 1;
  real c0 = 0;
  std::string name() const { return mode == bending_mode::helfrich ? "helfrich" : "plate"; }
  void build(build_context &c) const {
    if (mode != bending_mode::off && w > 0)
      c.energies.push_back(store_bending::create(c.table, c.x, mode, w, c0, c.store));
  }
};

/// A torsion target along `axis` (`torsion_target`): only the geodesic
/// torsion `W(a, N x a)` is asked for, the principal curvatures and their
/// directions are free, so curvature can change and the surface buckle while
/// it twists. The sign of `tau` is the handedness; `d_tau` ramps it per frame
/// up to `tau_max`. Pair it with `bending` (its own weight) to set what the
/// curvature change costs.
struct torsion {
  real w = 1;
  real tau = 1;
  real d_tau = 0;
  real tau_max = 10;
  axis_spec axis{};
  std::string name() const { return "torsion"; }
  void build(build_context &c) const {
    const real t = ramp(tau, d_tau, tau_max, c.frame);
    if (w > 0)
      c.energies.push_back(torsion_target::create(c.table, c.x, c.axis(axis), t, w, c.store));
  }
};

/// A twist-rate target on the store (`twist_rate`): the principal field turns
/// at `omega` along `axis`, scale-free in `|z|`. Ramped per frame.
struct twist_rate {
  real w = 1;
  real omega = 1;
  real d_omega = 0;
  real omega_max = 2;
  axis_spec axis{};
  std::string name() const { return "twist_rate"; }
  void build(build_context &c) const {
    if (w > 0)
      c.energies.push_back(engman::twist_rate::create(c.table, c.x, c.hs, c.axis(axis),
                                                      ramp(omega, d_omega, omega_max, c.frame), w, c.store));
  }
};

/// The twisted spin-2 Dirichlet energy (`twisted_spin2`). Quadratic, but also
/// lowered by shrinking `|z|`: pair it with `anisotropy_hold`.
struct twisted_dirichlet {
  real w = 1;
  real omega = 1;
  real d_omega = 0;
  real omega_max = 2;
  axis_spec axis{};
  std::string name() const { return "twisted_dirichlet"; }
  void build(build_context &c) const {
    if (w > 0)
      c.energies.push_back(twisted_spin2::create(c.table, c.x, c.axis(axis),
                                                 ramp(omega, d_omega, omega_max, c.frame), w, c.store));
  }
};

/// Holds each edge's anisotropy `|z|` at its frame-start value
/// (`deviator_hold`), its phase free.
struct anisotropy_hold {
  real w = 1;
  std::string name() const { return "anisotropy_hold"; }
  void build(build_context &c) const {
    if (w > 0)
      c.energies.push_back(deviator_hold::create(c.table, c.x, deviator_hold::magnitudes(c.hs), w, c.store));
  }
};

/// Per-face area against the rest positions (`face_area`).
struct face_area {
  mode m = mode::hard;
  real w_soft = 1e2;
  std::string name() const { return std::string("face_area ") + mode_name(m); }
  void build(build_context &c) const {
    if (m != mode::off)
      c.add(engman::face_area::create(c.tris, engman::face_area::areas(c.tris, c.xbar), c.pos), m, w_soft);
  }
};

/// The shear modulus: traceless strain against the rest positions
/// (`triangle_conformal`).
struct conformal {
  real mu = 1;
  std::string name() const { return "conformal"; }
  void build(build_context &c) const {
    if (!(mu > 0))
      return;
    const std::vector<hepworth::sim_block::ptr> pos_only{c.pos};
    for (const auto &t : c.tris)
      c.energies.push_back(hepworth::block::triangle_conformal::create(std::vector<index_t>{t[0], t[1], t[2]},
                                                                       c.xbar, std::sqrt(mu), pos_only));
  }
};

} // namespace terms

/// Every term. Add new ones here.
using term = std::variant<terms::coupling, terms::memory, terms::store_metric, terms::bending, terms::torsion,
                          terms::twist_rate, terms::twisted_dirichlet, terms::anisotropy_hold,
                          terms::face_area, terms::conformal>;

inline void build(const term &t, build_context &c) {
  std::visit([&](const auto &v) { v.build(c); }, t);
}
inline std::string name(const term &t) {
  return std::visit([](const auto &v) { return v.name(); }, t);
}

/// The first term of type `T` in `ts`, or null: for editing a preset,
/// `find<terms::torsion>(cfg.terms)->tau = 2`.
template <class T>
T *find(std::vector<term> &ts) {
  for (term &t : ts)
    if (T *p = std::get_if<T>(&t))
      return p;
  return nullptr;
}

} // namespace engman
} // namespace gaudi

#endif
