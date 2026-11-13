#pragma once

#include <array>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>

#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/common.h"
#include "ring_variational_fit.hpp"
#include "shell_constraints.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

/// Height of a circle of signed curvature `k` over its tangent line, at the
/// tangent-line offset `t`. Exact, with no subtraction, so `k = 0` is the line.
/// False past the circle's silhouette (`|k t| > 1`).
inline bool circle_sag(real k, real t, real &h) {
  const real kt = k * t;
  const real under = real(1) - kt * kt;
  if (!(under >= real(0)) || !std::isfinite(under))
    return false;
  h = k * t * t / (real(1) + std::sqrt(under));
  return std::isfinite(h);
}

/// Height of the osculating torus over its tangent plane at the foot. `u`
/// runs along the tube direction (curvature `k_hi`) and `v` along the major
/// circle (curvature `k_lo`). The tube circle is revolved about an axis along
/// the tube direction, a signed distance `1/k_lo` along the normal. Both sags
/// are exact and subtraction-free, so flat and cylindrical limits are smooth:
///
///   h = sag(k_hi, u) + sag(k_lo / (1 − k_lo sag(k_hi, u)), v)
///
/// The 2-jet is ½(k_hi u² + k_lo v²). False past either silhouette.
inline bool osculating_torus_height(real k_hi, real k_lo, real u, real v, real &h) {
  real ht = 0;
  if (!circle_sag(k_hi, u, ht))
    return false;
  const real denom = real(1) - k_lo * ht;
  if (!(std::abs(denom) > real(1e-12)))
    return false;
  real hv = 0;
  if (!circle_sag(k_lo / denom, v, hv))
    return false;
  h = ht + hv;
  return std::isfinite(h);
}

/// Store shape at the foot: the two curvatures and their tangent directions
/// as coefficients on the chart's `(T, B)`. The tube is the larger `|κ|`.
struct torus_ring_shape {
  real k_hi = 0;
  real k_lo = 0;
  real hx = 1;
  real hy = 0;
  real lx = 0;
  real ly = 1;
};

/// Rigid pose of the store's torus against the ring: an offset of the foot
/// along the normal and a tilt of the normal toward `T` and `B`. The spin about
/// the normal is the store's principal angle and is not fitted.
struct torus_ring_pose {
  real offset = 0;
  real tilt_t = 0;
  real tilt_b = 0;
};

/// Tilted chart. `N` leans by `(tilt_t, tilt_b)`. `T` is re-projected and `B`
/// completes the frame, so the store's principal coefficients carry over.
struct torus_ring_frame {
  vec3 T = vec3::UnitX();
  vec3 B = vec3::UnitY();
  vec3 N = vec3::UnitZ();
  bool ok = false;
};

inline torus_ring_frame tilt_frame(const vec3 &T, const vec3 &B, const vec3 &N, const torus_ring_pose &pose) {
  torus_ring_frame f;
  vec3 n = N + pose.tilt_t * T + pose.tilt_b * B;
  if (!(n.squaredNorm() > real(1e-16)) || !n.allFinite())
    return f;
  n.normalize();
  vec3 t = T - n * T.dot(n);
  if (!(t.squaredNorm() > real(1e-16)))
    return f;
  t.normalize();
  f.N = n;
  f.T = t;
  f.B = n.cross(t);
  f.ok = f.B.allFinite();
  return f;
}

/// One ring point, measured from the foot, with its row scale.
struct torus_ring_point {
  vec3 d = vec3::Zero();
  real s = 0;
};

/// Height and its gradient. `∂sag/∂t = k t / s` and `∂sag/∂k = t² / (s (1 + s))`
/// with `s = sqrt(1 − k² t²)`, chained through the revolved curvature
/// `m = k_lo / (1 − k_lo sag(k_hi, u))`.
inline bool osculating_torus_height_grad(real k_hi, real k_lo, real u, real v, real &h, real &hu, real &hv) {
  real ht = 0;
  if (!circle_sag(k_hi, u, ht))
    return false;
  const real su = std::sqrt(std::max(real(0), real(1) - k_hi * k_hi * u * u));
  const real denom = real(1) - k_lo * ht;
  if (!(su > real(1e-9)) || !(std::abs(denom) > real(1e-12)))
    return false;
  const real m = k_lo / denom;
  real hv_sag = 0;
  if (!circle_sag(m, v, hv_sag))
    return false;
  const real sv = std::sqrt(std::max(real(0), real(1) - m * m * v * v));
  if (!(sv > real(1e-9)))
    return false;
  const real ht_u = k_hi * u / su;
  const real m_u = k_lo * k_lo * ht_u / (denom * denom);
  const real sag_m = v * v / (sv * (real(1) + sv));
  h = ht + hv_sag;
  hu = ht_u + sag_m * m_u;
  hv = m * v / sv;
  return std::isfinite(h) && std::isfinite(hu) && std::isfinite(hv);
}

/// Landing of `d` on the posed torus, measured from the unposed foot. The
/// default drops it along the posed normal and keeps its tangential
/// coordinates. `closest` takes the closest point instead, by foot-point
/// iteration on the height graph: drop onto the graph, project onto its
/// tangent plane there, repeat. About ten times the cost. `gap` is the signed distance along the
/// surface normal at the landing, so `d − land = gap · normal` and the pose fit
/// minimizes exactly the row energy of the landing.
inline bool land_on_posed_torus(const torus_ring_shape &shape, const torus_ring_frame &f, real offset,
                                const vec3 &d, vec3 &land, real &gap, bool closest = false) {
  const vec3 eh = shape.hx * f.T + shape.hy * f.B;
  if (!closest) {
    // Vertical drop along the posed normal: the tangential coordinates stay.
    const vec3 el0 = shape.lx * f.T + shape.ly * f.B;
    const vec3 p0 = d - offset * f.N;
    const real u0 = p0.dot(eh);
    const real v0 = p0.dot(el0);
    real h0 = 0;
    if (!osculating_torus_height(shape.k_hi, shape.k_lo, u0, v0, h0))
      return false;
    gap = p0.dot(f.N) - h0;
    land = offset * f.N + u0 * eh + v0 * el0 + h0 * f.N;
    return land.allFinite() && std::isfinite(gap);
  }
  const vec3 el = shape.lx * f.T + shape.ly * f.B;
  const vec3 p = d - offset * f.N;
  const real pu = p.dot(eh);
  const real pv = p.dot(el);
  const real pn = p.dot(f.N);
  real u = pu, v = pv, h = 0, hu = 0, hv = 0;
  for (int it = 0; it < 8; ++it) {
    if (!osculating_torus_height_grad(shape.k_hi, shape.k_lo, u, v, h, hu, hv))
      return false;
    const real inv = real(1) / std::sqrt(real(1) + hu * hu + hv * hv);
    const real dist = ((pu - u) * -hu + (pv - v) * -hv + (pn - h)) * inv;
    u = pu + dist * hu * inv;
    v = pv + dist * hv * inv;
  }
  if (!osculating_torus_height_grad(shape.k_hi, shape.k_lo, u, v, h, hu, hv))
    return false;
  const real inv = real(1) / std::sqrt(real(1) + hu * hu + hv * hv);
  gap = ((pu - u) * -hu + (pv - v) * -hv + (pn - h)) * inv;
  land = offset * f.N + u * eh + v * el + h * f.N;
  return land.allFinite() && std::isfinite(gap);
}

/// How a ring row treats the volume of its one-ring cone.
enum class ring_volume_mode {
  /// No volume term.
  off,
  /// Restore a fraction `weight / (1 + weight)` of the volume error at each
  /// projection.
  soft,
  /// Restore the volume exactly at each projection.
  hard,
  /// Variational: the volume is a term of the ring's own fit. Each projection
  /// minimizes the ring points' gaps, a pull of both curvatures toward the
  /// store (`store_weight`), and the one-ring medial volume's error (`weight`)
  /// together over the pose and both curvatures, then lands the ring on the
  /// fitted shape. The volume is measured on the model, between the landed
  /// patch and each point's curvature-sphere centre, so it is a cone at an
  /// umbilic and a wedge on a tube, and only a change of shape can meet it.
  /// See `ring_variational_fit.hpp`.
  fit,
  /// Corotational (ARAP / Müller): the row lands its ring as it does with
  /// the volume off, then the landings are adjusted so each fan prism, from a
  /// fan triangle to its points' tube-axis points, keeps its rest volume. The
  /// rest pose and skeleton are taken at the step's first projection and
  /// corotated onto the current ones by a Procrustes fit over the ring and
  /// skeleton together, so rigid motion is free. Weights in `corotational`;
  /// with its volume weight 0 the row is the bare row. See
  /// `corotational_ring_step`.
  corotational,
};

/// Volume of the cone from the foot's medial (spine) point to its one-ring fan.
///
/// The spine point sits `1 / k_hi` along the normal from the foot: the tube
/// centre of the osculating torus, the centre of the larger curvature sphere,
/// the local medial-axis point. The cone is one tetrahedron per fan triangle
/// `(foot, x_j, x_j+1)` with that apex, so its volume is the slice of solid
/// between this patch and the medial axis. On a tube of radius `r` it scales
/// like `r^2` per unit length, the tube's own volume.
///
/// The spine point is found at the step's first projection and then held
/// fixed in space for the step, the fixed spine reference: a store curvature
/// edit does not move the apex, so it cannot inflate or deflate the cone on its
/// own. The term holds the cone volume at its first-projection value by the
/// smallest change of the ring's targets and the foot's position in the row's
/// metric. As the surface flattens the apex runs off to infinity. It is capped at `apex_cap`
/// times the ring's reach, and in that limit the term becomes "no
/// area-weighted mean normal motion". The apex side is fixed at the first
/// projection, so a curvature crossing zero mid-step does not flip the cone.
struct ring_volume_config {
  ring_volume_mode mode = ring_volume_mode::off;
  /// Soft mode only. The fraction of the volume error restored at each
  /// projection is `weight / (1 + weight)`.
  real weight = 1;
  /// Land the ring on the fitted shape before the volume correction. False,
  /// with a volume `mode` on, keeps the ring points where they are and applies
  /// only the volume correction, with the spine read from the store at the
  /// vertex chart, and skips the fit. The row is then a pure local volume
  /// constraint, for use beside a curvature solver that already realizes the
  /// shape. With `mode` off the row always lands: "volume only" with no volume
  /// would be a row that does nothing.
  bool land = true;

  /// Whether the row lands on its shape (see `land`).
  bool lands() const { return land || mode == ring_volume_mode::off; }
  /// Farthest apex, in multiples of the ring's reach.
  real apex_cap = 10;
  /// Weight of the foot's absolute row, relative to a one-ring row (soft and
  /// hard modes; `fit` mode has no foot row). The ring
  /// rows are differences, so their forces cancel in sum and cannot move the
  /// patch as a whole toward or away from the spine; this row can, and the
  /// volume correction moves the foot along with the ring. Between
  /// corrections it holds the foot where it is, like a little extra mass. 0
  /// drops it, leaving only shape corrections relative to the foot.
  real foot_weight = 1;
  /// Hold a persistent reference instead of re-recording the volume at each
  /// step. The reference is the cone's volume per unit fan area, a local
  /// medial depth (about r / 3 on a tube) that survives remeshing: it is set
  /// the first time a vertex is seen and the demo carries it through the
  /// remesher as a vertex datum. Each step's target is that depth times the
  /// ring's area at the step's start. Without it the reference resets every
  /// step, so whatever a step loses (to remesh splits, or to the compromise
  /// with other terms in the solve) is kept. Needs the row to be given a depth
  /// store; without one it falls back to the per-step reference.
  bool persistent = true;
  /// `fit` mode only: pull of both curvatures toward the store, per one-ring
  /// point, in the height a mismatch bends a point at the reach. The volume
  /// term (`weight`, on the same per-point scale) has to win against it to
  /// change the shape.
  real store_weight = 10;
  /// `fit` mode only: the low curvature's store pull as a multiple of
  /// `store_weight`, so the volume acts through the tube radius rather than
  /// by bending the patch around the major circle.
  real lo_stiffness = 10;
  /// `fit` mode only, experimental: hold each one-ring point's tube-axis
  /// point in space for the step, so the fit's volume is the region between
  /// the moving patch and a fixed skeleton, a real volume, and free the foot
  /// offset to move the patch relative to it through the foot row. In the
  /// monge growth run this is unstable: every vertex sits in about seven
  /// rings' cells, each restores its whole volume error at every projection,
  /// and the tube inflates or tears (a weaker foot row makes it worse). Off,
  /// the axis points follow the model: stable, but the volume then measures
  /// the local tube radius, which an undulating surface can hold while its
  /// enclosed volume shrinks.
  bool held_medial = false;
  /// `corotational` mode: the spring weights of its local solve.
  corotational_config corotational{};
};

/// Signed volume of the cone from `c` over the fan `(f, d_j, d_j+1)`, for the
/// first `n` entries of `d` in cyclic order: one tetrahedron
/// `(c, f, d_j, d_j+1)` per fan triangle. `f` is the foot, the origin unless
/// given.
inline real fan_cone_volume(const std::vector<vec3> &d, size_t n, const vec3 &c,
                            const vec3 &f = vec3::Zero()) {
  real v = 0;
  const vec3 w = f - c;
  for (size_t j = 0; j < n; ++j)
    v += (d[j] - c).dot((d[(j + 1) % n] - c).cross(w));
  return v / real(6);
}

/// Area of the fan `(0, d_j, d_j+1)` over the first `n` entries of `d`.
inline real fan_area(const std::vector<vec3> &d, size_t n) {
  real a = 0;
  for (size_t j = 0; j < n; ++j)
    a += real(0.5) * d[j].cross(d[(j + 1) % n]).norm();
  return a;
}

/// Gradient of `fan_cone_volume` in each `d_j`, and in the foot `f` when
/// `g_foot` is given.
inline void fan_cone_gradient(const std::vector<vec3> &d, size_t n, const vec3 &c, std::vector<vec3> &g,
                              const vec3 &f = vec3::Zero(), vec3 *g_foot = nullptr) {
  g.assign(n, vec3::Zero());
  const vec3 w = f - c;
  vec3 gf = vec3::Zero();
  for (size_t j = 0; j < n; ++j) {
    const vec3 here = d[j] - c;
    const vec3 next = d[(j + 1) % n] - c;
    const vec3 prev = d[(j + n - 1) % n] - c;
    g[j] = (next.cross(w) + w.cross(prev)) / real(6);
    gf += here.cross(next) / real(6);
  }
  if (g_foot)
    *g_foot = gf;
}

/// Moves the first `n` targets (the cyclic one-ring, relative to the foot) and
/// the foot itself by the smallest change in the metric
/// `sum_j s_j^2 |dt_j|^2 + foot^2 |delta|^2` that brings the cone volume about
/// the fixed apex `c` to `v0` (`hard`), or a fraction `weight / (1 + weight)`
/// of the way (`soft`). The ring rows are relative to the foot, so a foot
/// shift `delta` carries the whole ring with it: the cone is evaluated over
/// the ring at `delta + t_j` with the foot at `delta`. Returns `delta`; with
/// `foot` 0 it stays zero. Linearized steps, at most eight, absorb the
/// volume's cubic dependence; a ring a tenth off converges in four or five.
inline vec3 hold_fan_cone_volume(std::vector<vec3> &t, size_t n, const std::vector<real> &scales,
                                 const vec3 &c, real v0, const ring_volume_config &cfg, real foot = 0) {
  vec3 delta = vec3::Zero();
  if (cfg.mode == ring_volume_mode::off || n < 3 || !std::isfinite(v0))
    return delta;
  const real fraction = cfg.mode == ring_volume_mode::hard
                            ? real(1)
                            : std::max(cfg.weight, real(0)) / (real(1) + std::max(cfg.weight, real(0)));
  if (!(fraction > real(0)))
    return delta;
  const real foot2 = foot > real(0) ? foot * foot : real(0);
  std::vector<vec3> P(n), g;
  auto volume = [&]() {
    for (size_t j = 0; j < n; ++j)
      P[j] = delta + t[j];
    return fan_cone_volume(P, n, c, delta);
  };
  const real now = volume();
  const real goal = now + fraction * (v0 - now);
  for (int it = 0; it < 8; ++it) {
    const real err = goal - volume();
    if (!(std::abs(err) > real(1e-14) * std::abs(goal)))
      break;
    vec3 g_foot;
    fan_cone_gradient(P, n, c, g, delta, &g_foot);
    // d/d delta moves the foot and every ring point together.
    vec3 g_delta = g_foot;
    for (size_t j = 0; j < n; ++j)
      g_delta += g[j];
    real denom = foot2 > real(0) ? g_delta.squaredNorm() / foot2 : real(0);
    for (size_t j = 0; j < n; ++j)
      denom += g[j].squaredNorm() / std::max(scales[j] * scales[j], real(1e-12));
    if (!(denom > real(1e-30)))
      break;
    const real lambda = err / denom;
    for (size_t j = 0; j < n; ++j)
      t[j] += lambda * g[j] / std::max(scales[j] * scales[j], real(1e-12));
    if (foot2 > real(0))
      delta += lambda * g_delta / foot2;
  }
  return delta;
}

/// Spine point relative to the unposed foot: the posed foot plus `1 / |k_hi|`
/// along the posed normal on `side`, capped at `cap`.
inline vec3 ring_spine_point(real k_hi, real offset, const vec3 &normal, real side, real cap) {
  real a = std::abs(k_hi) > real(0) ? real(1) / std::abs(k_hi) : cap;
  if (!(a < cap) || !std::isfinite(a))
    a = cap;
  return (offset + side * a) * normal;
}

/// Per-row state of the one-ring volume term: the one-ring's size at the front
/// of the ring points, the spine point held in space for the step, and the
/// cone volume it holds.
struct ring_volume_term {
  size_t n_ring1 = 0;
  bool apex_set = false;
  vec3 apex = vec3::Zero();
  real v0 = std::numeric_limits<real>::quiet_NaN();
  /// This vertex's entry in the persistent depth store (NaN until set), or
  /// null for a per-step reference.
  real *depth = nullptr;

  /// The spine point relative to the foot `xi`. The first call places it from
  /// the tube curvature (`1 / |k_hi|` along the posed normal, on the side of
  /// `k_hi`, capped) and holds it in space; later calls only re-express it
  /// relative to the moving foot. `d(k)` reads ring point `k` relative to the
  /// unposed foot.
  template <class D>
  vec3 spine(real k_hi, real offset, const vec3 &normal, D &&d, const ring_volume_config &cfg,
             const vec3 &xi) {
    if (!apex_set) {
      real reach2 = 0;
      for (size_t k = 0; k < n_ring1; ++k)
        reach2 = std::max(reach2, d(k).squaredNorm());
      const real side = k_hi < real(0) ? real(-1) : real(1);
      apex = xi + ring_spine_point(k_hi, offset, normal, side, cfg.apex_cap * std::sqrt(reach2));
      apex_set = true;
    }
    return apex - xi;
  }

  /// Records the cone volume of the points `d` on the first call, then moves
  /// the one-ring `target`s, and returns the foot shift, to hold it.
  template <class D>
  vec3 hold(std::vector<vec3> &target, D &&d, const std::vector<real> &scales, const vec3 &spine,
            const ring_volume_config &cfg) {
    if (cfg.mode == ring_volume_mode::off || cfg.mode == ring_volume_mode::fit ||
        cfg.mode == ring_volume_mode::corotational || n_ring1 < 3)
      return vec3::Zero();
    if (!std::isfinite(v0)) {
      std::vector<vec3> d0(n_ring1);
      for (size_t k = 0; k < n_ring1; ++k)
        d0[k] = d(k);
      v0 = fan_cone_volume(d0, n_ring1, spine);
      const real area = fan_area(d0, n_ring1);
      if (cfg.persistent && depth && area > real(0)) {
        if (std::isfinite(*depth))
          v0 = *depth * area;
        else
          *depth = v0 / area;
      }
    }
    return hold_fan_cone_volume(target, n_ring1, scales, spine, v0, cfg, cfg.foot_weight);
  }

  /// Whether the row carries the foot's absolute row.
  /// In `fit` mode only with `held_medial`. Without it the volume is
  /// intrinsic to the fitted shape, the ring's difference rows carry it, and
  /// the fit's normal offset is only the model's least-squares placement
  /// against the ring: moving the foot to it every projection ratchets the
  /// vertex along its normal.
  bool foot_row(const ring_volume_config &cfg) const {
    return cfg.mode != ring_volume_mode::off && (cfg.mode != ring_volume_mode::fit || cfg.held_medial) &&
           cfg.foot_weight > real(0) && n_ring1 >= 3;
  }

  /// `fit` mode with `held_medial`: the step's skeleton in world coordinates.
  bool medial_set = false;
  std::vector<vec3> medial_world;
  vec3 medial_foot_world = vec3::Zero();
  /// `corotational` mode: the step's rest pose and skeleton.
  corotational_rest corot;
};

struct osculating_torus_ring_config {
  /// Row scale of the 2-ring points relative to the 1-ring.
  real two_ring = 0.5;
  /// Gauss-Newton iterations on the pose per projection.
  int iterations = 4;
  /// Levenberg pull of the pose toward the chart pose (offset over the ring
  /// reach, and each tilt), and of fitted curvatures toward the store.
  real prior = 1e-3;
  /// The store drives changes, not absolute values. At the first projection of
  /// a step each row fits both curvatures and the principal angle to its ring
  /// and records the store's reading. Later projections target that fit plus
  /// the store's change since, so a bias in the store's discrete reading does
  /// not flow the mesh.
  bool relative = false;
  /// Scale on the store's change since the reference under `relative`. 0
  /// holds the frame-start fit.
  real store_gain = 1;
  /// Land on the closest point of the torus instead of dropping along the
  /// normal. About ten times the cost of a projection.
  bool closest_point = false;
  /// One-ring cone volume about the spine. See `ring_volume_config`.
  ring_volume_config volume{};
};

/// The principal directions turned by `a` about the normal.
inline torus_ring_shape spin_shape(torus_ring_shape s, real a) {
  const real c = std::cos(a), sn = std::sin(a);
  const real hx = c * s.hx - sn * s.hy, hy = sn * s.hx + c * s.hy;
  const real lx = c * s.lx - sn * s.ly, ly = sn * s.lx + c * s.ly;
  s.hx = hx;
  s.hy = hy;
  s.lx = lx;
  s.ly = ly;
  return s;
}

struct torus_ring_fit {
  torus_ring_pose pose;
  /// The fitted shape. Equal to the input shape unless curvatures were fitted.
  torus_ring_shape shape;
  /// Turn of the principal directions about the normal, when fitted.
  real spin = 0;
  real cost0 = 0;
  real cost = 0;
  bool ok = false;
};

/// Least-squares pose of a torus on the ring. With `fit_curvature` set, both
/// curvatures are solved too, and the principal directions stay the input's.
/// The point set that lands at the starting pose is kept fixed through the
/// solve, and a trial that drops one of them is rejected.
inline torus_ring_fit fit_torus_ring_pose(const torus_ring_shape &shape0, const vec3 &T, const vec3 &B,
                                          const vec3 &N, const std::vector<torus_ring_point> &pts,
                                          std::vector<char> &in, const osculating_torus_ring_config &cfg,
                                          bool fit_curvature = false) {
  torus_ring_fit fit;
  fit.shape = shape0;
  const size_t n = pts.size();
  in.assign(n, 0);
  real reach2 = 0;
  for (const torus_ring_point &pt : pts)
    reach2 = std::max(reach2, pt.d.squaredNorm());
  const real L = std::sqrt(reach2);
  if (!(L > real(1e-12)))
    return fit;
  const real prior = std::max(cfg.prior, real(0));
  const int np = fit_curvature ? 6 : 3;

  // Parameters: offset, tilt_t, tilt_b, then k_hi, k_lo and the spin of the
  // principal directions about the normal when the shape is fitted.
  using param = Eigen::Matrix<real, 6, 1>;
  auto unpack = [&](const param &x, torus_ring_pose &pose, torus_ring_shape &shape) {
    pose.offset = x[0];
    pose.tilt_t = x[1];
    pose.tilt_b = x[2];
    shape = shape0;
    if (fit_curvature) {
      shape.k_hi = x[3];
      shape.k_lo = x[4];
      shape = spin_shape(shape, x[5]);
    }
  };
  // Residual vector: one height gap per kept point, then the prior rows.
  auto residuals = [&](const param &x, std::vector<real> &r, bool select) -> bool {
    torus_ring_pose pose;
    torus_ring_shape shape;
    unpack(x, pose, shape);
    const torus_ring_frame f = tilt_frame(T, B, N, pose);
    if (!f.ok)
      return false;
    r.clear();
    for (size_t i = 0; i < n; ++i) {
      if (!select && !in[i])
        continue;
      vec3 land;
      real gap = 0;
      const bool ok = land_on_posed_torus(shape, f, pose.offset, pts[i].d, land, gap, cfg.closest_point);
      if (select) {
        in[i] = ok ? 1 : 0;
        if (!ok)
          continue;
      } else if (!ok) {
        return false;
      }
      r.push_back(pts[i].s * gap);
    }
    // Prior rows: the offset over the reach and each tilt, then each fitted
    // curvature in the height it moves a point at the reach.
    const real sp = std::sqrt(prior);
    r.push_back(sp * pose.offset / L);
    r.push_back(sp * pose.tilt_t);
    r.push_back(sp * pose.tilt_b);
    if (fit_curvature) {
      r.push_back(sp * real(0.5) * L * L * (shape.k_hi - shape0.k_hi));
      r.push_back(sp * real(0.5) * L * L * (shape.k_lo - shape0.k_lo));
      r.push_back(sp * real(0.5) * L * L * std::abs(shape0.k_hi - shape0.k_lo) * x[5]);
    }
    return true;
  };
  auto cost_of = [](const std::vector<real> &r) {
    real c = 0;
    for (const real v : r)
      c += v * v;
    return c;
  };

  param x = param::Zero();
  x[3] = shape0.k_hi;
  x[4] = shape0.k_lo;
  std::vector<real> r;
  if (!residuals(x, r, true))
    return fit;
  int kept = 0;
  for (const char c : in)
    kept += c ? 1 : 0;
  if (kept < np)
    return fit;
  fit.cost0 = cost_of(r);
  real cost = fit.cost0;

  param h;
  h << real(1e-6) * L, real(1e-6), real(1e-6), real(1e-6) / L, real(1e-6) / L, real(1e-6);
  std::vector<real> rp, rm, rt;
  for (int it = 0; it < std::max(cfg.iterations, 0); ++it) {
    const size_t m = r.size();
    Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> J(static_cast<Eigen::Index>(m), np);
    bool ok = true;
    for (int k = 0; k < np && ok; ++k) {
      param xp = x, xm = x;
      xp[k] += h[k];
      xm[k] -= h[k];
      ok = residuals(xp, rp, false) && residuals(xm, rm, false) && rp.size() == m && rm.size() == m;
      for (size_t i = 0; ok && i < m; ++i)
        J(static_cast<Eigen::Index>(i), k) = (rp[i] - rm[i]) / (real(2) * h[k]);
    }
    if (!ok)
      break;
    Eigen::Map<const vecX> rv(r.data(), static_cast<Eigen::Index>(m));
    const matX H = J.transpose() * J;
    const vecX g = J.transpose() * rv;
    const vecX step = -H.ldlt().solve(g);
    if (!step.allFinite())
      break;
    bool accepted = false;
    real a = 1;
    for (int ls = 0; ls < 5; ++ls, a *= real(0.5)) {
      param trial = x;
      for (int k = 0; k < np; ++k)
        trial[k] += a * step[k];
      if (!residuals(trial, rt, false))
        continue;
      const real c = cost_of(rt);
      if (c < cost) {
        x = trial;
        r = rt;
        cost = c;
        accepted = true;
        break;
      }
    }
    if (!accepted)
      break;
  }
  unpack(x, fit.pose, fit.shape);
  fit.spin = fit_curvature ? x[5] : real(0);
  fit.cost = cost;
  fit.ok = true;
  return fit;
}

/// Process-wide tallies over every `osculating_torus_ring::project` call, for
/// probes. `move` sums |landing − point| and `reach` sums |point|, so
/// `move / reach` is the mean relative pull. Reset with `reset()`.
struct osculating_torus_ring_stats {
  /// Off by default so the solve pays nothing. Probes turn it on.
  static inline std::atomic<bool> enabled{false};
  static inline std::atomic<long> calls{0};
  static inline std::atomic<long> no_frame{0};
  static inline std::atomic<long> no_fit{0};
  static inline std::atomic<long> points{0};
  static inline std::atomic<long> landed{0};
  static inline std::atomic<double> move{0};
  static inline std::atomic<double> reach{0};
  /// Store against mesh at the chart pose, over the 1-ring: Σ k_hi, Σ k_lo,
  /// and Σ h_mesh h_store, Σ h_store² for the regression slope of the mesh's
  /// normal heights on the store's.
  static inline std::atomic<long> shapes{0};
  static inline std::atomic<double> k_hi{0};
  static inline std::atomic<double> k_lo{0};
  static inline std::atomic<double> hm_hs{0};
  static inline std::atomic<double> hs_hs{0};
  static void reset() {
    shapes = 0;
    k_hi = 0;
    k_lo = 0;
    hm_hs = 0;
    hs_hs = 0;
    calls = 0;
    no_frame = 0;
    no_fit = 0;
    points = 0;
    landed = 0;
    move = 0;
    reach = 0;
  }
};

/// Osculating-torus projector on one vertex's 1-ring and 2-ring. The store
/// owns the shape: both curvatures and the principal angle come from the live
/// edge store at the foot. The ring owns the pose: the foot's normal offset and
/// the normal's tilt are fitted to the ring points. Each ring point is then
/// landed on that posed torus. Rows are differences `x_j − x_i`, so the
/// projection system carries the absolute position.
class osculating_torus_ring : public block_constraint {
public:
  DEFINE_CREATE_FUNC(osculating_torus_ring)

  /// `ids` is the foot, then the 1-ring, then the 2-ring. `scales` holds one
  /// row scale per ring point.
  /// `n_ring1` is how many of the ring points, from the front, are the cyclic
  /// one-ring the volume term's cone is built on.
  osculating_torus_ring(std::vector<index_t> ids, std::vector<real> scales, const real &w,
                        std::vector<sim_block::ptr> blocks, const std::vector<vertex_chart> *charts,
                        const kusama::weingarten::vertex_stars *stars, osculating_torus_ring_config cfg,
                        size_t n_ring1 = 0, real *depth = nullptr)
      : block_constraint(std::move(ids), w, std::move(blocks)), _scales(std::move(scales)), _charts(charts),
        _stars(stars), _cfg(cfg) {
    _volume.n_ring1 = std::min(n_ring1, _scales.size());
    _volume.depth = depth;
    _rest_initialized = true;
  }

  virtual std::string name() { return "osculating_torus_ring"; }

  virtual void project(const vecX &q, vecX &p) {
    if (_blocks.empty() || _ids.empty())
      return;
    const index_t vi = _ids[0];
    const vec3 xi = _blocks[0]->get_vec3(vi, q);
    const size_t n = _scales.size();
    std::vector<torus_ring_point> pts(n);
    for (size_t k = 0; k < n; ++k) {
      pts[k].d = _blocks[0]->get_vec3(_ids[k + 1], q) - xi;
      pts[k].s = _scales[k];
    }
    std::vector<vec3> target(n);
    for (size_t k = 0; k < n; ++k)
      target[k] = pts[k].d;

    shape_w W;
    vec3 T, B, N;
    _last_ok = false;
    vec3 spine = vec3::Zero();
    bool have_spine = false;
    bool variational = false;
    vec3 variational_foot = vec3::Zero();
    bool corot_done = false;
    vec3 corot_foot = vec3::Zero();
    auto ring_d = [&](size_t k) -> const vec3 & { return pts[k].d; };
    if (vertex_row_frame(_blocks, _charts, vi, q, W, T, B, N, _stars)) {
      W = shape_with_axes(W);
      if (W.ok) {
        const monge_axes ax = monge_axes_of(W);
        torus_ring_shape shape{
            .k_hi = ax.k_hi, .k_lo = ax.k_lo, .hx = ax.hx, .hy = ax.hy, .lx = ax.lx, .ly = ax.ly};
        if (_cfg.relative) {
          if (!_reference) {
            std::vector<char> in0;
            const torus_ring_fit ref = fit_torus_ring_pose(shape, T, B, N, pts, in0, _cfg, true);
            _reference = true;
            _ref_ok = ref.ok;
            _store0_hi = shape.k_hi;
            _store0_lo = shape.k_lo;
            _fit0_hi = ref.shape.k_hi;
            _fit0_lo = ref.shape.k_lo;
            _spin0 = ref.spin;
          }
          if (_ref_ok) {
            shape = spin_shape(shape, _spin0);
            shape.k_hi = _fit0_hi + _cfg.store_gain * (shape.k_hi - _store0_hi);
            shape.k_lo = _fit0_lo + _cfg.store_gain * (shape.k_lo - _store0_lo);
          }
        }
        if (osculating_torus_ring_stats::enabled.load(std::memory_order_relaxed)) {
          using stats = osculating_torus_ring_stats;
          double a = 0, b = 0;
          for (size_t k = 0; k < n; ++k) {
            if (_scales[k] != real(1))
              continue;
            real h = 0;
            if (!osculating_torus_height(shape.k_hi, shape.k_lo, pts[k].d.dot(ax.hx * T + ax.hy * B),
                                         pts[k].d.dot(ax.lx * T + ax.ly * B), h))
              continue;
            a += pts[k].d.dot(N) * h;
            b += h * h;
          }
          stats::shapes.fetch_add(1, std::memory_order_relaxed);
          stats::k_hi.fetch_add(shape.k_hi, std::memory_order_relaxed);
          stats::k_lo.fetch_add(shape.k_lo, std::memory_order_relaxed);
          stats::hm_hs.fetch_add(a, std::memory_order_relaxed);
          stats::hs_hs.fetch_add(b, std::memory_order_relaxed);
        }
        variational = _cfg.volume.mode == ring_volume_mode::fit && _cfg.volume.lands();
        if (variational && fit_variational(shape, T, B, N, xi, pts, target, variational_foot))
          _last_ok = true;
        if (!_cfg.volume.lands()) {
          spine = _volume.spine(shape.k_hi, 0, N, ring_d, _cfg.volume, xi);
          have_spine = true;
        }
        std::vector<char> in;
        const torus_ring_fit fit = _cfg.volume.lands() && !variational
                                       ? fit_torus_ring_pose(shape, T, B, N, pts, in, _cfg)
                                       : torus_ring_fit{};
        const torus_ring_frame f = tilt_frame(T, B, N, fit.pose);
        static const int debug_vert =
            std::getenv("GAUDI_TORUS_DEBUG") ? std::atoi(std::getenv("GAUDI_TORUS_DEBUG")) : -1;
        if (debug_vert == -2 && !_debug_printed) {
          _debug_printed = true;
          real worst = 0;
          for (size_t k = 0; k < n; ++k) {
            vec3 land;
            real gap = 0;
            if (land_on_posed_torus(shape, f, fit.pose.offset, pts[k].d, land, gap, _cfg.closest_point))
              worst = std::max(worst, (land - pts[k].d).norm());
          }
          std::cerr << "[torus_all] " << vi << " " << xi.transpose() << " " << shape.k_hi << " " << shape.k_lo
                    << " " << fit.cost << " " << worst << " " << n << "\n";
        }
        if (vi == debug_vert && !_debug_printed) {
          _debug_printed = true;
          std::cerr << "[torus_ring] v " << vi << " x " << xi.transpose() << " N " << N.transpose() << " T "
                    << T.transpose() << " k_hi " << shape.k_hi << " k_lo " << shape.k_lo << " ref fit "
                    << _fit0_hi << " " << _fit0_lo << " store0 " << _store0_hi << " " << _store0_lo
                    << " e_hi " << (ax.hx * T + ax.hy * B).transpose() << " pose " << fit.pose.offset << " "
                    << fit.pose.tilt_t << " " << fit.pose.tilt_b << " cost " << fit.cost0 << " -> "
                    << fit.cost << "\n";
          for (size_t k = 0; k < n; ++k) {
            vec3 land;
            real gap = 0;
            const bool ok =
                land_on_posed_torus(shape, f, fit.pose.offset, pts[k].d, land, gap, _cfg.closest_point);
            std::cerr << "  pt " << _ids[k + 1] << " s " << _scales[k] << " d " << pts[k].d.transpose()
                      << " ok " << ok << " gap " << gap << " move " << (land - pts[k].d).norm() << "\n";
          }
        }
        if (fit.ok && f.ok) {
          for (size_t k = 0; k < n; ++k) {
            if (!in[k])
              continue;
            vec3 land;
            real gap = 0;
            if (land_on_posed_torus(shape, f, fit.pose.offset, pts[k].d, land, gap, _cfg.closest_point))
              target[k] = land;
          }
          _last_ok = true;
          _last_fit = fit;
          spine = _volume.spine(shape.k_hi, fit.pose.offset, f.N, ring_d, _cfg.volume, xi);
          have_spine = true;
          if (_cfg.volume.mode == ring_volume_mode::corotational)
            corot_done = corotational_targets(shape, f, fit.pose.offset, pts, target, corot_foot);
        }
      }
    }
    vec3 foot_shift = variational ? variational_foot : vec3::Zero();
    if (have_spine)
      foot_shift = _volume.hold(target, ring_d, _scales, spine, _cfg.volume);
    if (corot_done)
      foot_shift = corot_foot;
    using stats = osculating_torus_ring_stats;
    if (stats::enabled.load(std::memory_order_relaxed)) {
      stats::calls.fetch_add(1, std::memory_order_relaxed);
      if (!_last_ok) {
        stats::no_fit.fetch_add(1, std::memory_order_relaxed);
      } else {
        long landed = 0;
        double mv = 0, rc = 0;
        for (size_t k = 0; k < n; ++k) {
          const real m = (target[k] - pts[k].d).norm();
          landed += m > real(0) ? 1 : 0;
          mv += m;
          rc += pts[k].d.norm();
        }
        stats::points.fetch_add(static_cast<long>(n), std::memory_order_relaxed);
        stats::landed.fetch_add(landed, std::memory_order_relaxed);
        stats::move.fetch_add(mv, std::memory_order_relaxed);
        stats::reach.fetch_add(rc, std::memory_order_relaxed);
      }
    }

    index_t row = _id0;
    for (size_t k = 0; k < n; ++k) {
      vec3 rowp = _w * _scales[k] * target[k];
      if (!rowp.allFinite())
        rowp = vec3::Zero();
      p.block(row, 0, 3, 1) = rowp;
      row += 3;
    }
    if (_volume.foot_row(_cfg.volume)) {
      vec3 rowp = _w * _cfg.volume.foot_weight * (xi + foot_shift);
      if (!rowp.allFinite())
        rowp = vec3::Zero();
      p.block(row, 0, 3, 1) = rowp;
      row += 3;
    }
  }

  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    if (_blocks.empty() || !_blocks[0] || _ids.empty())
      return;
    const index_t ii = _blocks[0]->get_offset_idx(_ids[0]);
    for (size_t k = 0; k < _scales.size(); ++k) {
      const index_t jj = _blocks[0]->get_offset_idx(_ids[k + 1]);
      const real s = _w * _scales[k];
      for (int ax = 0; ax < 3; ++ax) {
        triplets.push_back(trip(id0 + ax, ii + ax, -s));
        triplets.push_back(trip(id0 + ax, jj + ax, s));
      }
      id0 += 3;
    }
    if (_volume.foot_row(_cfg.volume)) {
      for (int ax = 0; ax < 3; ++ax)
        triplets.push_back(trip(id0 + ax, ii + ax, _w * _cfg.volume.foot_weight));
      id0 += 3;
    }
  }

  bool last_ok() const { return _last_ok; }
  const torus_ring_fit &last_fit() const { return _last_fit; }
  /// Cone volume the volume term holds (NaN before the first projection).
  real reference_volume() const { return _volume.v0; }

  /// `fit` volume mode: the variational fit over pose and both curvatures on
  /// the torus model (vertical landing, as the row's default), with the
  /// store's curvatures as the pull target. Writes the landed `target`s and
  /// the posed foot's offset from the unposed one.
  bool fit_variational(const torus_ring_shape &shape, const vec3 &T, const vec3 &B, const vec3 &N,
                       const vec3 &xi, const std::vector<torus_ring_point> &pts, std::vector<vec3> &target,
                       vec3 &foot_offset) {
    const size_t n = pts.size();
    thread_local std::vector<vec3> d;
    thread_local std::vector<real> sc;
    d.resize(n);
    sc.resize(n);
    for (size_t k = 0; k < n; ++k) {
      d[k] = pts[k].d;
      sc[k] = pts[k].s;
    }
    auto land = [&](const ring_fit_params &x, const vec3 &dk) {
      ring_model_point e;
      const torus_ring_frame fr = tilt_frame(T, B, N, {.offset = x[0], .tilt_t = x[1], .tilt_b = x[2]});
      if (!fr.ok)
        return e;
      const vec3 eh = shape.hx * fr.T + shape.hy * fr.B;
      const vec3 el = shape.lx * fr.T + shape.ly * fr.B;
      const vec3 p = dk - x[0] * fr.N;
      const real u = p.dot(eh), v = p.dot(el);
      real h = 0, hu = 0, hv = 0;
      if (!osculating_torus_height_grad(x[3], x[4], u, v, h, hu, hv))
        return e;
      e.gap = p.dot(fr.N) - h;
      e.land = x[0] * fr.N + u * eh + v * el + h * fr.N;
      e.normal = (fr.N - hu * eh - hv * el).normalized();
      e.ok = e.land.allFinite() && e.normal.allFinite();
      return e;
    };
    auto foot = [&](const ring_fit_params &x, vec3 &F, vec3 &NF) {
      const torus_ring_frame fr = tilt_frame(T, B, N, {.offset = x[0], .tilt_t = x[1], .tilt_b = x[2]});
      NF = fr.ok ? fr.N : N;
      F = x[0] * NF;
    };
    ring_fit_params start;
    start << 0, 0, 0, shape.k_hi, shape.k_lo;
    ring_variational_config vcfg{.w_store = _cfg.volume.store_weight,
                                 .lo_stiffness = _cfg.volume.lo_stiffness,
                                 .w_volume = _cfg.volume.weight,
                                 .prior = _cfg.prior,
                                 .iterations = _cfg.iterations,
                                 .apex_cap = _cfg.volume.apex_cap};
    // A held skeleton gives the foot offset a meaning: it moves the patch
    // relative to the skeleton. Without one it is only a nuisance and is pinned.
    if (_cfg.volume.held_medial)
      vcfg.offset_weight = 0;
    const size_t n1 = _volume.n_ring1;
    const real cap = vcfg.apex_cap * ring_reach(d, n1);
    // The step's skeleton, held in space, re-expressed about the moving foot.
    if (_cfg.volume.held_medial && !_volume.medial_set) {
      std::vector<vec3> M;
      vec3 MF;
      if (!ring_medial_images(start, d, n1, land, foot, cap, M, MF))
        return false;
      _volume.medial_world.resize(n1);
      for (size_t j = 0; j < n1; ++j)
        _volume.medial_world[j] = xi + M[j];
      _volume.medial_foot_world = xi + MF;
      _volume.medial_set = true;
    }
    ring_medial_hold hold;
    if (_volume.medial_set) {
      hold.M.resize(n1);
      for (size_t j = 0; j < n1; ++j)
        hold.M[j] = _volume.medial_world[j] - xi;
      hold.MF = _volume.medial_foot_world - xi;
    }
    const ring_medial_hold *held = _volume.medial_set ? &hold : nullptr;
    if (!std::isfinite(_volume.v0)) {
      real v = 0, a = 0;
      if (!ring_medial_volume(start, d, n1, land, foot, cap, v, a, held))
        return false;
      _volume.v0 = v;
      if (_cfg.volume.persistent && _volume.depth) {
        if (std::isfinite(*_volume.depth))
          _volume.v0 = *_volume.depth * a;
        else
          *_volume.depth = v / a;
      }
    }
    if (_have_theta)
      start = _theta;
    const ring_variational_fit fit =
        fit_ring_variational(start, shape.k_hi, shape.k_lo, d, sc, n1, _volume.v0, vcfg, land, foot, held);
    if (!fit.ok)
      return false;
    _theta = fit.theta;
    _have_theta = true;
    static const int debug_vert =
        std::getenv("GAUDI_RING_DEBUG") ? std::atoi(std::getenv("GAUDI_RING_DEBUG")) : -1;
    if (debug_vert == static_cast<int>(_ids[0]) && !_debug_printed) {
      _debug_printed = true;
      real v_start = 0, a_start = 0;
      ring_fit_params s0;
      s0 << 0, 0, 0, shape.k_hi, shape.k_lo;
      ring_medial_volume(s0, d, n1, land, foot, vcfg.apex_cap * ring_reach(d, n1), v_start, a_start);
      std::cerr << "[ring_fit] v " << _ids[0] << " store k " << shape.k_hi << " " << shape.k_lo << " depth "
                << (_volume.depth ? *_volume.depth : real(-1)) << " v0 " << _volume.v0 << " V@store "
                << v_start << " A " << a_start << " | fit k " << fit.theta[3] << " " << fit.theta[4]
                << " offset " << fit.theta[0] << " V " << fit.volume << " A " << fit.area << "\n";
    }
    for (size_t k = 0; k < n; ++k) {
      const ring_model_point e = land(fit.theta, d[k]);
      if (e.ok)
        target[k] = e.land;
    }
    vec3 NF;
    foot(fit.theta, foot_offset, NF);
    return true;
  }

  /// `corotational` mode: the skeleton from the posed torus's normals at the
  /// landings, then the corotational step on top of the row's landings.
  /// Rewrites `target` (relative to the foot) and gives the foot's shift.
  bool corotational_targets(const torus_ring_shape &shape, const torus_ring_frame &f, real offset,
                            const std::vector<torus_ring_point> &pts, std::vector<vec3> &target, vec3 &foot) {
    const size_t n = pts.size(), n1 = _volume.n_ring1;
    if (n1 < 3)
      return false;
    thread_local std::vector<vec3> d, S, y;
    d.resize(n);
    for (size_t k = 0; k < n; ++k)
      d[k] = pts[k].d;
    const real cap = _cfg.volume.apex_cap * ring_reach(d, n1);
    const vec3 eh = shape.hx * f.T + shape.hy * f.B;
    const vec3 el = shape.lx * f.T + shape.ly * f.B;
    S.resize(n1);
    for (size_t k = 0; k < n1; ++k) {
      const vec3 p = target[k] - offset * f.N;
      real h = 0, hu = 0, hv = 0;
      vec3 normal = f.N;
      if (osculating_torus_height_grad(shape.k_hi, shape.k_lo, p.dot(eh), p.dot(el), h, hu, hv))
        normal = (f.N - hu * eh - hv * el).normalized();
      S[k] = medial_image(target[k], normal, shape.k_hi, cap);
    }
    const vec3 SF = medial_image(offset * f.N, f.N, shape.k_hi, cap);
    real *depth = _cfg.volume.corotational.persistent ? _volume.depth : nullptr;
    const bool first = !_volume.corot.set;
    if (!corotational_ring_step(d, n1, target, S, SF, _volume.corot, _cfg.volume.corotational, depth, y))
      return false;
    static const int debug_vert =
        std::getenv("GAUDI_RING_DEBUG") ? std::atoi(std::getenv("GAUDI_RING_DEBUG")) : -1;
    if (debug_vert == static_cast<int>(_ids[0]) && (first || _corot_calls++ == 40)) {
      real v0 = 0, vl = 0, vy = 0;
      for (real v : _volume.corot.V0)
        v0 += v;
      auto prism_total = [&](const std::vector<vec3> &pts3, const vec3 &F) {
        real v = 0;
        for (size_t j = 0; j < n1; ++j) {
          const size_t k = (j + 1) % n1;
          v += tet_volume(F, pts3[j], pts3[k], SF) + tet_volume(pts3[j], pts3[k], SF, S[j]) +
               tet_volume(pts3[k], SF, S[j], S[k]);
        }
        return v;
      };
      std::vector<vec3> yr(y.begin() + 1, y.begin() + 1 + static_cast<long>(n1));
      vl = prism_total(target, vec3::Zero());
      vy = prism_total(yr, y[0]);
      real mean_dn = 0;
      for (size_t k = 0; k < n1; ++k)
        mean_dn += (y[k + 1] - target[k]).dot(f.N) / real(n1);
      std::cerr << "[corot] v " << _ids[0] << (first ? " first" : " later") << " k_hi " << shape.k_hi
                << " V0 " << v0 << " V(landing, current skeleton) " << vl << " V(targets) " << vy
                << " foot.N " << y[0].dot(f.N) << " ring dn " << mean_dn << " depth "
                << (depth ? *depth : real(-1)) << "\n";
    }
    for (size_t k = 0; k < n; ++k)
      target[k] = y[k + 1] - y[0];
    foot = y[0];
    return true;
  }

  static real ring_reach(const std::vector<vec3> &d, size_t n1) {
    real r2 = 0;
    for (size_t k = 0; k < n1 && k < d.size(); ++k)
      r2 = std::max(r2, d[k].squaredNorm());
    return std::sqrt(r2);
  }

  std::vector<real> _scales;
  const std::vector<vertex_chart> *_charts = nullptr;
  const kusama::weingarten::vertex_stars *_stars = nullptr;
  osculating_torus_ring_config _cfg;
  bool _last_ok = false;
  bool _debug_printed = false;
  /// Frame-start reference for `relative`: the ring's fitted curvatures and
  /// the store's reading at that first projection.
  bool _reference = false;
  bool _ref_ok = false;
  real _store0_hi = 0;
  real _store0_lo = 0;
  real _fit0_hi = 0;
  real _fit0_lo = 0;
  real _spin0 = 0;
  torus_ring_fit _last_fit;
  ring_volume_term _volume;
  int _corot_calls = 0;
  /// `fit` volume mode: the previous projection's fitted parameters, the warm
  /// start of the next.
  bool _have_theta = false;
  ring_fit_params _theta = ring_fit_params::Zero();
};

/// The 1-ring of `vi`, then the 2-ring, without the foot or repeats.
inline void torus_ring_ids(const asawa::shell::shell &shell, index_t vi, std::vector<index_t> &ring1,
                           std::vector<index_t> &ring2) {
  ring1.clear();
  ring2.clear();
  const int nV = shell.vert_count();
  if (vi < 0 || vi >= nV)
    return;
  std::vector<char> seen(static_cast<size_t>(nV), 0);
  seen[static_cast<size_t>(vi)] = 1;
  const auto one = shell.get_one_ring(asawa::shell::vert_id(vi));
  for (const asawa::shell::VertId nj : one) {
    const int j = static_cast<int>(nj);
    if (j < 0 || j >= nV || seen[static_cast<size_t>(j)])
      continue;
    seen[static_cast<size_t>(j)] = 1;
    ring1.push_back(j);
  }
  for (const index_t j : ring1) {
    for (const asawa::shell::VertId nk : shell.get_one_ring(asawa::shell::vert_id(j))) {
      const int k = static_cast<int>(nk);
      if (k < 0 || k >= nV || seen[static_cast<size_t>(k)])
        continue;
      seen[static_cast<size_t>(k)] = 1;
      ring2.push_back(k);
    }
  }
}

/// One osculating-torus row per vertex. `blocks` is positions, then the edge
/// store that `stars` reads.
inline void init_osculating_torus_ring(const asawa::shell::shell &shell,
                                       std::vector<projection_constraint::ptr> &constraints, const real &w,
                                       std::vector<sim_block::ptr> blocks,
                                       const std::vector<vertex_chart> *charts,
                                       const kusama::weingarten::vertex_stars *stars,
                                       osculating_torus_ring_config cfg = {},
                                       std::vector<real> *depth = nullptr) {
  if (!(w > real(0)) || !std::isfinite(w) || blocks.empty() || !blocks[0] || !charts)
    return;
  if (depth && depth->size() != static_cast<size_t>(shell.vert_count()))
    depth = nullptr;
  if (blocks.size() < 2 || !blocks[1])
    stars = nullptr;
  std::vector<index_t> ring1, ring2;
  for (int vi = 0; vi < shell.vert_count(); ++vi) {
    torus_ring_ids(shell, vi, ring1, ring2);
    if (ring1.size() < 3)
      continue;
    std::vector<index_t> ids{vi};
    std::vector<real> scales;
    ids.reserve(1 + ring1.size() + ring2.size());
    scales.reserve(ring1.size() + ring2.size());
    for (const index_t j : ring1) {
      ids.push_back(j);
      scales.push_back(real(1));
    }
    if (cfg.two_ring > real(0)) {
      for (const index_t j : ring2) {
        ids.push_back(j);
        scales.push_back(cfg.two_ring);
      }
    }
    real *entry = depth ? &(*depth)[static_cast<size_t>(vi)] : nullptr;
    constraints.push_back(osculating_torus_ring::create(std::move(ids), std::move(scales), w, blocks, charts,
                                                        stars, cfg, ring1.size(), entry));
  }
}

} // namespace block
} // namespace hepworth
} // namespace gaudi
