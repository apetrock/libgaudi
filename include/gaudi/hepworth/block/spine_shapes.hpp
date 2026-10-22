#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include "gaudi/common.h"
#include "shape_objectives.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

/// Which spine won the ring. `hit` is the ray landing where `ok` is set.
struct spine_landing {
  spine_kind fit = spine_kind::none;
  std::vector<char> ok;
  std::vector<vec3> hit;
};

/// Principal frame at the foot. `tube` is the tight direction, `ruling` the
/// other. Spokes are already in the foot frame.
struct spine_curvature {
  vec3 tube = vec3::UnitX();
  vec3 ruling = vec3::UnitY();
  vec3 normal = vec3::UnitZ();
  real k_hi = 0;
  real k_lo = 0;
};

/// Osculating canal footing from neighborhood samples (typically a 2-ring).
/// Center lies on the normal line through the foot: `center = t * normal`,
/// `radius = |t|`.
struct osculating_spine {
  bool ready = false;
  real t = 0;
  real radius = 0;
  vec3 center = vec3::Zero();
};

namespace spine_detail {

struct spine_fit {
  bool ready = false;
  real l2 = 0;
  std::vector<char> ok;
  std::vector<vec3> hit;
};

/// Move `p` along the canal radial from spine point `ps` toward radius
/// `radius`. `step` ∈ [0, 1] is the fraction of the residual to take
/// (`1` = full land, `0` = hold `p`).
inline bool land_from_point(const vec3 &ps, real radius, const vec3 &p,
                            vec3 &out, real step = real(1)) {
  if (!ps.allFinite() || !p.allFinite() || !(radius > real(1e-12)))
    return false;
  const vec3 ray = p - ps;
  const real q = ray.norm();
  if (!(q > real(1e-12)) || !std::isfinite(q))
    return false;
  real a = std::isfinite(step) ? step : real(1);
  if (a < real(0))
    a = real(0);
  if (a > real(1))
    a = real(1);
  const real r = (real(1) - a) * q + a * radius;
  out = ps + r * (ray / q);
  return out.allFinite();
}

inline bool nearest_on_line(const vec3 &origin, const vec3 &dir, const vec3 &p,
                            vec3 &ps) {
  if (!(dir.squaredNorm() > real(1e-16)))
    return false;
  ps = origin + (p - origin).dot(dir) * dir;
  return ps.allFinite();
}

/// `radial` points from the circle's center toward `spine`.
inline bool nearest_on_circle(const vec3 &spine, const vec3 &radial,
                              const vec3 &ruling, real major, const vec3 &p,
                              vec3 &ps) {
  const vec3 v = p - spine;
  const real tang = v.dot(ruling);
  const real radial_off = v.dot(radial);
  const real along = major + radial_off;
  const real rho = std::hypot(along, tang);
  if (!std::isfinite(rho))
    return false;
  if (!(rho > real(1e-18))) {
    ps = spine;
    return ps.allFinite();
  }
  const real scale = major / rho;
  const real denom = along + rho;
  const real sag =
      (std::abs(denom) > real(1e-18)) ? (scale * tang * tang) / denom : real(0);
  ps = spine - sag * radial + (scale * tang) * ruling;
  return ps.allFinite();
}

inline bool land_ellipsoid(const vec3 &center, const vec3 &eh, const vec3 &el,
                           const vec3 &un, real a, real b, real c,
                           const vec3 &p, vec3 &out, real step = real(1)) {
  if (!(a > real(1e-12)) || !(b > real(1e-12)) || !(c > real(1e-12)))
    return false;
  const vec3 d = p - center;
  const real A = d.dot(eh) / a;
  const real B = d.dot(el) / b;
  const real C = d.dot(un) / c;
  const real Q = A * A + B * B + C * C;
  if (!(Q > real(1e-16)) || !std::isfinite(Q))
    return false;
  const vec3 hard = center + d / std::sqrt(Q);
  if (!hard.allFinite())
    return false;
  real s = std::isfinite(step) ? step : real(1);
  if (s < real(0))
    s = real(0);
  if (s > real(1))
    s = real(1);
  out = (real(1) - s) * p + s * hard;
  return out.allFinite();
}

template <typename Land>
spine_fit fit_spokes(const std::vector<vec3> &spokes, Land &&land) {
  spine_fit fit;
  const size_t n = spokes.size();
  fit.ok.assign(n, 0);
  fit.hit.assign(n, vec3::Zero());
  int finite = 0;
  for (size_t i = 0; i < n; ++i) {
    const vec3 &p = spokes[i];
    if (!p.allFinite())
      continue;
    ++finite;
    vec3 hit = vec3::Zero();
    if (!land(p, hit))
      return fit;
    fit.ok[i] = 1;
    fit.hit[i] = hit;
    fit.l2 += (hit - p).squaredNorm();
  }
  fit.ready = finite > 0;
  return fit;
}

inline void take_if_lower(const spine_fit &fit, spine_kind kind,
                          spine_landing &best, real &best_l2, bool &have,
                          real tol = 0) {
  if (!fit.ready)
    return;
  if (have && !(fit.l2 + tol < best_l2))
    return;
  best.fit = kind;
  best.ok = fit.ok;
  best.hit = fit.hit;
  best_l2 = fit.l2;
  have = true;
}

/// Circle through the foot with center on `N`: `C = t N`, `r = |t|`.
/// Each spoke gives `t = ‖q‖² / (2 N·q)`. Cylinder fits in the `(tube, N)`
/// plane and requires a tube-direction signal (a pure ruling offset is the
/// major circle, not the canal tube). Sphere uses the full spoke.
inline osculating_spine fit_osculating_normal_line(
    const vec3 &N, const std::vector<vec3> &spokes, bool cylinder_plane,
    const vec3 &tube) {
  osculating_spine out;
  if (!N.allFinite() || !(N.squaredNorm() > real(1e-16)))
    return out;
  const vec3 Nn = N.normalized();
  vec3 Tt = tube - Nn * tube.dot(Nn);
  if (cylinder_plane && !(Tt.squaredNorm() > real(1e-16)))
    return out;
  if (cylinder_plane)
    Tt.normalize();
  real sum_t = 0;
  real sum_w = 0;
  real tube_signal = 0;
  for (const vec3 &p : spokes) {
    if (!p.allFinite())
      continue;
    const real along_tube = cylinder_plane ? std::abs(p.dot(Tt)) : real(1);
    if (cylinder_plane)
      tube_signal = std::max(tube_signal, along_tube);
    // Major-circle samples (no tube component) are not canal observations.
    if (cylinder_plane && !(along_tube > real(1e-6)))
      continue;
    vec3 q = p;
    if (cylinder_plane)
      q = (p.dot(Tt)) * Tt + (p.dot(Nn)) * Nn;
    const real nd = Nn.dot(q);
    if (!(std::abs(nd) > real(1e-12)))
      continue;
    const real t = q.squaredNorm() / (real(2) * nd);
    if (!std::isfinite(t))
      continue;
    const real w = std::abs(nd);
    sum_t += w * t;
    sum_w += w;
  }
  // Ruling-only rings give no canal tube. Leave `ready` false.
  if (cylinder_plane && !(tube_signal > real(1e-6)))
    return out;
  if (!(sum_w > real(1e-12)))
    return out;
  out.t = sum_t / sum_w;
  out.radius = std::abs(out.t);
  if (!(out.radius > real(1e-12)))
    return out;
  out.center = out.t * Nn;
  out.ready = out.center.allFinite();
  return out;
}

/// Osculating tube-spine centerpoint at the foot. Kept under a ratio target;
/// the axis (ruling / major circle) is replaced from the new κ. Sliding this
/// point with `r` would pin the foot and leave spoke residuals O(θ²).
inline vec3 osculating_centerpoint(const osculating_spine &geom, const vec3 &N,
                                   real r_geom_fallback, real sgn) {
  const vec3 Nn = N.normalized();
  const vec3 fallback = sgn * r_geom_fallback * Nn;
  if (!geom.ready || !(geom.radius > real(1e-12)))
    return fallback;
  return geom.center;
}

} // namespace spine_detail

struct spine_set {
  bool cylinder = false;
  bool sphere = false;
  bool torus = false;
  bool ellipsoid = false;
  bool any() const { return cylinder || sphere || torus || ellipsoid; }
};

inline spine_set all_spine_kinds() {
  spine_set s;
  s.cylinder = true;
  s.sphere = true;
  s.torus = true;
  s.ellipsoid = true;
  return s;
}

inline spine_set spine_set_of(const std::vector<objective> &objectives) {
  spine_set s;
  for (const objective &o : objectives) {
    const auto *sp = std::get_if<spine_fairing>(&o);
    if (!sp)
      continue;
    switch (sp->shape) {
    case spine_kind::cylinder:
      s.cylinder = true;
      break;
    case spine_kind::sphere:
      s.sphere = true;
      break;
    case spine_kind::torus:
      s.torus = true;
      break;
    case spine_kind::ellipsoid:
      s.ellipsoid = true;
      break;
    case spine_kind::none:
      break;
    }
  }
  return s;
}

/// Rest canal centers from an osculating fit on neighborhood samples.
/// `C_cyl` is the cylinder-plane fit (tube spine); `C_sph` is the full-spoke
/// sphere fit. Each projector holds this rest center and lands with the goal
/// κ from the solve.
struct spine_rest {
  bool ready = false;
  vec3 C_cyl = vec3::Zero();
  vec3 C_sph = vec3::Zero();
  real r_geom = 0; ///< tube radius from the canal fit
};

/// Osculating normal-line canal fit on `fit_spokes` (typically a 2-ring).
inline spine_rest rest_centers_from_spokes(const spine_curvature &frame,
                                           const std::vector<vec3> &fit_spokes) {
  spine_rest out;
  vec3 normal = frame.normal;
  if (!normal.allFinite() || !(normal.squaredNorm() > real(1e-16)))
    return out;
  normal.normalize();
  vec3 tube = frame.tube - normal * frame.tube.dot(normal);
  const real sgn = (frame.k_hi >= real(0)) ? real(1) : real(-1);
  const real r_tgt =
      std::abs(frame.k_hi) > real(1e-8) ? real(1) / std::abs(frame.k_hi)
                                        : real(0);
  const osculating_spine cyl = spine_detail::fit_osculating_normal_line(
      normal, fit_spokes, /*cylinder_plane=*/true,
      tube.squaredNorm() > real(1e-16) ? tube.normalized() : vec3::UnitX());
  const osculating_spine sph = spine_detail::fit_osculating_normal_line(
      normal, fit_spokes, /*cylinder_plane=*/false, tube);
  const real r_fit =
      cyl.ready ? cyl.radius : (sph.ready ? sph.radius : r_tgt);
  out.C_cyl =
      spine_detail::osculating_centerpoint(cyl, normal, r_fit, sgn);
  out.C_sph =
      spine_detail::osculating_centerpoint(sph, normal, r_fit, sgn);
  out.r_geom = cyl.ready ? cyl.radius : out.C_cyl.norm();
  out.ready = out.C_cyl.allFinite() &&
              (cyl.ready || sph.ready || out.r_geom > real(1e-12));
  return out;
}

/// Rest centers from an osculating canal fit on `fit_spokes` (typically
/// 2-ring). Land `project_spokes` (1-ring) with the goal κ / axis from the
/// solve row `cr`. `step` ∈ [0, 1] is the canal-radial residual bite
/// (`1` = full land). Empty `fit_spokes` reuses `project_spokes` for the fit.
/// Vanishing major κ floors to a large-\(R\) torus (*Points as Tori*).
inline spine_landing choose_spine_shapes(
    const spine_curvature &cr, const std::vector<vec3> &project_spokes,
    spine_set kinds = all_spine_kinds(),
    const std::vector<vec3> &fit_spokes = {}, real step = real(1)) {
  spine_landing best;
  vec3 tube = cr.tube;
  vec3 ruling = cr.ruling;
  vec3 normal = cr.normal;
  if (!(std::abs(cr.k_hi) > real(1e-8)) || !normal.allFinite())
    return best;
  normal.normalize();
  tube -= normal * tube.dot(normal);
  ruling -= normal * ruling.dot(normal);
  if (!(tube.squaredNorm() > real(1e-16)) ||
      !(ruling.squaredNorm() > real(1e-16)))
    return best;
  tube.normalize();
  ruling = ruling - tube * ruling.dot(tube);
  if (!(ruling.squaredNorm() > real(1e-16)))
    return best;
  ruling.normalize();

  const real r_tgt = real(1) / std::abs(cr.k_hi);
  const real sgn = (cr.k_hi >= real(0)) ? real(1) : real(-1);
  constexpr real k_lo_floor = real(1e-8);
  const real abs_lo = std::max(std::abs(cr.k_lo), k_lo_floor);
  const bool curved = std::abs(cr.k_lo) > real(1e-12);
  const std::vector<vec3> &fit =
      fit_spokes.empty() ? project_spokes : fit_spokes;

  // Rest C from canal fit; goal κ (this `cr`) for radius / axis / major.
  spine_curvature fit_frame = cr;
  fit_frame.normal = normal;
  fit_frame.tube = tube;
  const spine_rest rest = rest_centers_from_spokes(fit_frame, fit);
  if (!rest.ready)
    return best;
  const vec3 C_cyl = rest.C_cyl;
  const vec3 C_sph = rest.C_sph.allFinite() ? rest.C_sph : rest.C_cyl;

  spine_detail::spine_fit cylinder;
  if (kinds.cylinder) {
    cylinder =
        spine_detail::fit_spokes(project_spokes, [&](const vec3 &p, vec3 &hit) {
          vec3 ps = vec3::Zero();
          return spine_detail::nearest_on_line(C_cyl, ruling, p, ps) &&
                 spine_detail::land_from_point(ps, r_tgt, p, hit, step);
        });
  }
  spine_detail::spine_fit sphere;
  if (kinds.sphere) {
    sphere =
        spine_detail::fit_spokes(project_spokes, [&](const vec3 &p, vec3 &hit) {
          return spine_detail::land_from_point(C_sph, r_tgt, p, hit, step);
        });
  }

  spine_detail::spine_fit torus;
  if (kinds.torus) {
    const real s_prod = (cr.k_hi * cr.k_lo >= real(0)) ? real(1) : real(-1);
    const real inv_lo = real(1) / abs_lo;
    const real major = std::max(real(0), inv_lo - s_prod * r_tgt);
    // New major axis through the held centerpoint: O = C − R radial, with
    // radial / ruling from the live (goal) principal frame.
    vec3 radial = normal - ruling * normal.dot(ruling);
    if (radial.squaredNorm() > real(1e-16) && std::isfinite(major) &&
        major > real(0)) {
      radial.normalize();
      const vec3 O_hint = sgn * s_prod * inv_lo * normal;
      if ((C_cyl - O_hint).dot(radial) < real(0))
        radial = -radial;
      torus = spine_detail::fit_spokes(
          project_spokes, [&](const vec3 &p, vec3 &hit) {
            vec3 ps = vec3::Zero();
            return spine_detail::nearest_on_circle(C_cyl, radial, ruling, major,
                                                  p, ps) &&
                   spine_detail::land_from_point(ps, r_tgt, p, hit, step);
          });
    }
  }

  spine_detail::spine_fit ellipsoid;
  if (kinds.ellipsoid && curved) {
    const real a = std::sqrt(r_tgt / std::abs(cr.k_hi));
    const real b = std::sqrt(r_tgt / std::abs(cr.k_lo));
    const vec3 un = sgn * normal;
    ellipsoid =
        spine_detail::fit_spokes(project_spokes, [&](const vec3 &p, vec3 &hit) {
          return spine_detail::land_ellipsoid(C_cyl, tube, ruling, un, a, b,
                                              r_tgt, p, hit, step);
        });
  }

  real best_l2 = 0;
  bool have = false;
  real reach2 = 0;
  for (const vec3 &p : project_spokes) {
    if (p.allFinite())
      reach2 = std::max(reach2, p.squaredNorm());
  }
  const real tol = reach2 * real(1e-8);
  // Torus before cylinder so large-R ≈ line keeps the torus.
  spine_detail::take_if_lower(torus, spine_kind::torus, best, best_l2, have);
  spine_detail::take_if_lower(sphere, spine_kind::sphere, best, best_l2, have,
                              tol);
  spine_detail::take_if_lower(ellipsoid, spine_kind::ellipsoid, best, best_l2,
                              have, tol);
  spine_detail::take_if_lower(cylinder, spine_kind::cylinder, best, best_l2,
                              have, tol);
  return best;
}

/// Per-ring spine drive: goal canal radius from κ vs rest radius from the
/// osculating canal fit, and mean spoke displacement after landing.
/// `fit_spokes` is the 2-ring; `project_spokes` are the landed 1-ring samples.
struct spine_ring_metric {
  bool ok = false;
  real r_tgt = 0;      ///< 1/|κ_hi| from goal S
  real r_geom = 0;     ///< tube radius from rest centers (fit S)
  real r_gap = 0;      ///< |r_tgt − r_geom|
  real R_tgt = 0;      ///< torus major ring: max(0, 1/|κ_lo| − s·r_tgt)
  real mean_hit_gap = 0; ///< mean ‖hit − spoke‖ over landed project spokes
  real max_hit_gap = 0;
  int n_land = 0;
  spine_kind fit = spine_kind::none;
  real k_hi = 0;
  real k_lo = 0;
};

inline spine_ring_metric measure_spine_ring(
    const spine_curvature &cr, const std::vector<vec3> &project_spokes,
    spine_set kinds = all_spine_kinds(),
    const std::vector<vec3> &fit_spokes = {}, real step = real(1)) {
  spine_ring_metric m;
  m.k_hi = cr.k_hi;
  m.k_lo = cr.k_lo;
  if (!(std::abs(cr.k_hi) > real(1e-8)))
    return m;
  m.r_tgt = real(1) / std::abs(cr.k_hi);
  {
    constexpr real k_lo_floor = real(1e-8);
    const real abs_lo = std::max(std::abs(cr.k_lo), k_lo_floor);
    const real s_prod = (cr.k_hi * cr.k_lo >= real(0)) ? real(1) : real(-1);
    const real inv_lo = real(1) / abs_lo;
    m.R_tgt = std::max(real(0), inv_lo - s_prod * m.r_tgt);
  }
  const std::vector<vec3> &fit =
      fit_spokes.empty() ? project_spokes : fit_spokes;
  const spine_rest rest = rest_centers_from_spokes(cr, fit);
  if (rest.ready)
    m.r_geom = rest.r_geom;
  m.r_gap = std::abs(m.r_tgt - m.r_geom);
  const spine_landing land =
      choose_spine_shapes(cr, project_spokes, kinds, fit, step);
  m.fit = land.fit;
  real sum = 0;
  int n = 0;
  for (size_t i = 0; i < project_spokes.size() && i < land.ok.size(); ++i) {
    if (!land.ok[i] || i >= land.hit.size())
      continue;
    const real g = (land.hit[i] - project_spokes[i]).norm();
    if (!std::isfinite(g))
      continue;
    sum += g;
    m.max_hit_gap = std::max(m.max_hit_gap, g);
    ++n;
  }
  m.n_land = n;
  if (n > 0)
    m.mean_hit_gap = sum / real(n);
  m.ok = true;
  return m;
}

} // namespace block
} // namespace hepworth
} // namespace gaudi
