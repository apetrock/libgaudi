#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <variant>
#include <vector>

#include "gaudi/common.h"

namespace gaudi {
namespace hepworth {
namespace block {

/// Graph over the tangent plane, except the two cylinders, which are the
/// osculating circular tubes of the principal curvatures. `cylinder_minor` is
/// the larger-magnitude curvature (the tight tube). `cylinder_major` is the
/// smaller one (the loose tube). The other principal direction is the ruling.
enum class monge_shape {
  quadric,
  hyperboloid,
  cylinder_minor,
  cylinder_major,
  sphere,
  osculating_torus
};

/// Fairing projector: the spoke lands on the graph. Residual is a length,
/// x - Q(x). Not a curvature descent.
struct monge_fairing {
  monge_shape shape = monge_shape::quadric;
};

/// Closest-point landing on the one ring. The samples are the ring positions.
/// Normals are the fan triangles of those positions. Quadric, sphere, and
/// cyclide are the Calder fits (`calder::quadric`, `calder::constrained_sphere`,
/// `calder::darboux_cyclide_normal_constrained`).
enum class albers_shape { quadric, sphere, darboux_cyclide };

struct albers_fairing {
  albers_shape shape = albers_shape::quadric;
};

/// One position ARAP row. `angle` is signed radians from the minor tangent
/// toward the major, the larger-magnitude principal direction and the other.
/// 0 leaves `x` on the minor tangent and `y` on the major. `x` and `y` are
/// extra fractions on those axes. 0 keeps that axis. 0.1 makes it 1.1 times
/// as long. A negative value shrinks it. (0, 0) leaves the spoke unchanged.
/// `power` multiplies each fraction by `|κ|^power` of that axis. 0 leaves the
/// fractions. Negative pushes a flat direction more. The tangents are rebuilt
/// every projection.
/// One closed stretch line on a torus: turns around the ring and around the
/// tube. Both 0 keeps `stretch_target::angle` as given; otherwise the demo
/// fits the torus at startup and sets the angle so the line closes,
/// `atan2(ring * sqrt(R^2 - r^2), tube * r)` from the tube circle. The line
/// crosses the curvature lines at that same angle everywhere on a torus of
/// revolution. {.ring = 1, .tube = 1} is the Villarceau circle (the Hopf
/// fiber).
struct stretch_winding {
  int ring = 0;
  int tube = 0;
  /// +1 or -1: which of the two mirror families.
  int handedness = 1;
};

/// Where a stretch row's axis comes from each step.
///  angle    — `angle` (or the winding) off the live principal frame.
///  material — set from `angle` at the first step, then carried with the
///             surface (`kusama::material_direction`: turned with each one
///             ring, optionally leaned back toward the frame and smoothed).
///  stripes  — `angle` off the frame, replaced by the level lines of a
///             warm-started stripe pattern across it (`kusama::stripes`), so
///             the axis follows closed curves.
enum class direction_source { angle, material, stripes };

struct stretch_target {
  real angle = 0;
  real x = 0;
  real y = 0;
  /// Exponent on `|κ|` along the push axis. 0 leaves `x` and `y`. Negative
  /// pushes a flat direction more.
  real power = 0;
  /// Sets `angle` from winding numbers at startup (`stretch_winding`).
  stretch_winding winding{};
  direction_source direction = direction_source::angle;
};

/// Directional growth of the edge store's curvature
/// (`shape_distortion::growth`), the shape-space match of `stretch_target`:
/// each frame the normal curvature along the axis at `angle` (radians from
/// the frame-start tensor's tight direction k_hi toward k_lo) is scaled by
/// `1 + x` and across it by `1 + y` (`D h0 D`, `D` the square roots on the
/// axis frame). A negative fraction bends less that way. `power`
/// multiplies each fraction by `(|k| / k_bar)^power`, `k` the frame-start
/// normal curvature along that axis and `k_bar` the store's mean
/// `(|k1| + |k2|) / 2`: positive grows the more curved edges faster (the
/// growth feeds on its own curvature), negative the flatter ones; 0 leaves
/// the fractions. `.winding` sets the angle at startup like stretch.
struct curvature_growth_target {
  real angle = 0;
  real x = 0;
  real y = 0;
  real power = 0;
  stretch_winding winding{};
};

/// An area band relative to the step-start area: `{1, 1}` keeps it,
/// `area_band::free()` leaves it free.
struct area_band {
  real lo = 1;
  real hi = 1;
  static constexpr area_band free() { return {real(0), std::numeric_limits<real>::infinity()}; }
  constexpr bool contains(real a) const { return lo <= a && a <= hi; }
  constexpr real clamp(real a) const { return std::clamp(a, lo, hi); }
};

/// The per-triangle `area` rows: weight and band.
struct area_config {
  real w = 0.5;
  area_band band{0.5, 1.5};
};

/// Locks the one ring's covariant anisotropy: the spoke rms along the tight
/// principal tangent (k_hi) over that along the other is held at `ratio`
/// (1 is round). Skew and rotation are free; the ring's area stays within
/// `area` of its frame-start area (`axis_ratio`).
struct axis_ratio_target {
  real ratio = 1;
  area_band area{};
};

/// Area-preserving shear of the one ring each frame: stretch by
/// `1 + shear` along the direction at `angle` and `1 - shear` across it,
/// against the frame-start ring. That ratio is held and the ring's area kept
/// within `area` of its frame-start area; the skew of that
/// axis and rotation are free (`axis_ratio`), so at 45 degrees from the
/// principal tangents it is independent of `axis_ratio_target`.
struct shear_target {
  real shear = 0;
  /// Radians from the tight principal tangent (k_hi) toward the other, as
  /// `stretch_target::angle`. NaN is the 45-degree bisector.
  real angle = std::numeric_limits<real>::quiet_NaN();
  /// Sets `angle` from winding numbers at startup (`stretch_winding`).
  stretch_winding winding{};
  area_band area{};
};

/// The anisotropy of one point, from its two principal curvatures `k0`, `k1`
/// (red and green in the demo; any order, any signs).
///
/// The ratio is `r = |k_hi| / |k_lo| >= 1`: how many times harder the
/// surface bends one way than the other. 1 is a round point (sphere-like,
/// `k0 = k1`) or a minimal saddle (`k0 = -k1`); infinity is a parabolic
/// point (one direction flat, like a cylinder).
///
/// `q = (r - 1) / (r + 1) = (|k_hi| - |k_lo|) / (|k_hi| + |k_lo|)` is the
/// same thing squeezed into [0, 1]: 0 round or minimal, 1 parabolic. With
/// `H = (k0 + k1) / 2` and `D = (k0 - k1) / 2`, `q = min(|H|, |D|) /
/// max(|H|, |D|)`, and `log r = 2 atanh(q)`.
inline real curvature_q(real k0, real k1) {
  const real a = std::abs(k0), b = std::abs(k1);
  return a + b > real(0) ? std::abs(a - b) / (a + b) : real(0);
}

/// `r` from `q`; infinite at `q = 1`.
inline real ratio_of_q(real q) {
  return q < real(1) ? (real(1) + q) / (real(1) - q) : std::numeric_limits<real>::infinity();
}

/// `q` from `r`; `r` below 1 is read as `1 / r`, infinity gives 1.
inline real q_of_ratio(real r) {
  if (!(r > real(0)))
    return real(0);
  if (r < real(1))
    r = real(1) / r;
  return std::isfinite(r) ? (r - real(1)) / (r + real(1)) : real(1);
}

/// A per-point number a gain can be multiplied by. All are plain functions
/// of the ratio `r = |k_hi| / |k_lo|`:
///  aniso     — `q = (r - 1) / (r + 1)` in [0, 1]: 0 round or minimal saddle,
///              1 parabolic.
///  log_ratio — `log r` in [0, log 1000]: 0 round, capped at ratio 1000.
///  ratio     — `1 / r` in [0, 1]: 1 round, 0 parabolic.
enum class curvature_metric { aniso, log_ratio, ratio };

inline real curvature_measure(curvature_metric m, real k0, real k1) {
  if (!std::isfinite(k0) || !std::isfinite(k1))
    return real(0);
  const real q = curvature_q(k0, k1);
  switch (m) {
  case curvature_metric::aniso:
    return q;
  case curvature_metric::log_ratio:
    return std::min(real(2) * std::atanh(std::min(q, real(1) - real(1e-15))),
                    std::log(real(1000)));
  case curvature_metric::ratio:
    return (real(1) - q) / (real(1) + q);
  }
  return real(0);
}

/// One curvature edit on the edge store (`shape_distortion::anisotropy`).
/// Designated initializers must follow this order; leave out what you do not
/// use. Ratio means `|k_hi| / |k_lo|` (larger over smaller curvature).
struct anisotropy_target {
  /// Ratio to reach on every edge. 1 is round (saddles: minimal), 3 is one
  /// direction bending three times the other. NaN is off.
  real ratio = std::numeric_limits<real>::quiet_NaN();
  /// `ratio` is given as `log(ratio)`: 0 is round, 1.1 is ratio 3.
  bool log = false;
  /// `ratio` multiplies each edge's frame-start ratio, every frame (with
  /// `log`, it is added): 0.9 or log -0.1 is 10% less anisotropic per frame.
  bool relative = false;
  /// Mean curvature times this, per frame: 1.05 bends 5% tighter. NaN is
  /// off.
  real mean = std::numeric_limits<real>::quiet_NaN();
};

/// Shear of the edge store's curvature pattern (`shape_distortion::shear`),
/// the shape-space match of `shear_target`: each frame every vertex ring's
/// frame-start tensors are pushed through the area-preserving shear that
/// stretches by `sqrt((1 + shear) / (1 - shear))` along `angle` and its
/// inverse across (`II -> S^-T II S^-1`, the ring skewed the way the
/// position row skews its spokes), then turned by the ring's best common
/// turn (ARAP in shape space). Only the direction of the lines is held: the
/// mean curvature and the anisotropy are free, so it composes with
/// `anisotropy_target` instead of fighting its ratio.
struct curvature_shear_target {
  real shear = 0;
  /// Radians from the tight principal tangent (k_hi) toward the other, as
  /// `shear_target::angle`. NaN is the 45-degree bisector.
  real angle = std::numeric_limits<real>::quiet_NaN();
  /// Sets `angle` from winding numbers at startup (`stretch_winding`).
  stretch_winding winding{};
};

/// One feature-size band of a `willmore_target`: sizes from `radius` up to
/// the next band's radius (the last band runs to the whole surface), in
/// mean edge lengths.
struct willmore_band {
  real radius = 0;
  /// Fraction of the Willmore flow taken per frame in this band. 0 leaves
  /// the band free: the row ignores it, so position rows can bend it.
  real gain = 0;
};

/// Willmore flow on the edge store with a feature-size equalizer
/// (`shape_distortion::willmore`). Designated initializers must follow this
/// order.
struct willmore_target {
  /// Bands in ascending radius. A first radius above 0 leaves the finer
  /// sizes free. Empty is off.
  std::vector<willmore_band> bands;
  /// 0 flows only the anisotropic part (rounds); 1 also the mean curvature.
  real mean = 0;
};

/// Which principal tangent a position row turns about: `k0` (larger signed
/// curvature, red) or `k1` (green).
enum class principal_axis { k0, k1 };

/// How `directed_tilt` turns each frame-start spoke `d` of a vertex's ring.
/// `a` is the chosen principal tangent, `N` the outward normal, `b = N x a`,
/// `l` the ring's mean spoke length, `theta = gain * metric * |κ|^power`.
///  rigid — tilts the whole ring about `a` (see `tilt_sign` for which way).
///  spin  — turns the whole ring about `N` by `theta`, counter-clockwise
///          seen from outside. The same sense everywhere.
///  twist — turns each spoke about `a` by `theta (d . a) / l`: the ring
///          warps into a saddle along `a` (normal offset ~ `(d . a)(d . b)`).
///  bend  — turns each spoke about `a` by `theta (d . b) / l`: the ring
///          curls across `a` (normal offset ~ `(d . b)^2`).
/// `twist` and `bend` do not depend on which way `a` points.
enum class tilt_mode { rigid, spin, twist, bend };

/// Which way `rigid` tilts. A principal tangent is a line, so `a` and `-a`
/// are equally valid, and tilting about them goes opposite ways.
///  gradient  — tilt toward rising `metric`: the turn is `gain * l *
///              (grad metric . b)` about `a`, the same whichever way `a`
///              points, and it fades to 0 where the metric does not change
///              across the axis. Local, no seams.
///  propagate — one sign of `a` for the whole mesh, chosen so neighboring
///              axes agree; turn `theta`. Seams remain only around round
///              (umbilic) points, which no sign choice avoids.
enum class tilt_sign { gradient, propagate };

/// Position row that turns each vertex's ring (as it was at the start of the
/// frame) about a principal axis or the normal. It does not settle: the ring
/// turns by its angle every frame, and the row also holds the ring at its
/// frame-start pose otherwise.
struct tilt_target {
  principal_axis axis = principal_axis::k1;
  tilt_mode mode = tilt_mode::rigid;
  tilt_sign sign = tilt_sign::gradient;
  /// Radians per frame, times `metric` and `|κ|^power`. 0 is off.
  real gain = 0;
  /// What multiplies `gain` (default `aniso`, in [0, 1]).
  curvature_metric metric = curvature_metric::aniso;
  /// Exponent on `|κ|`, the normal curvature along the axis. 0 leaves it
  /// out; 1 turns curved directions more, negative turns flat ones more.
  real power = 0;
};

enum class spine_kind { none, cylinder, sphere, torus, ellipsoid };

/// Canal spines. Each listed shape has its own spine. The spoke is the ray
/// from the nearest spine point, and the ring keeps the smaller residual.
struct spine_fairing {
  spine_kind shape = spine_kind::torus;
};

using objective = std::variant<monge_fairing, albers_fairing, spine_fairing>;

} // namespace block
} // namespace hepworth
} // namespace gaudi
