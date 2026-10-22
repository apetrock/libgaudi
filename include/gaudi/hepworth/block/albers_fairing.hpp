#pragma once

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

#include "gaudi/albers/darboux_cyclide.hpp"
#include "gaudi/albers/osculating_torus.hpp"
#include "gaudi/albers/quadric.hpp"
#include "gaudi/albers/sphere.hpp"
#include "gaudi/calder/weight_functions.hpp"
#include "shape_objectives.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

/// One sample of a Calder fit, already in the foot frame. `w` is the positive
/// area weight (Calder's primitive area). The kernel is applied here.
struct albers_fit_sample {
  vec3 dp = vec3::Zero();
  vec3 N = vec3::UnitZ();
  real w = 0;
};

/// Value of the 10-coefficient quadric. Matches `mk_quad_A`'s value row and
/// `quadric_grad`: cross terms carry the factor two.
inline real albers_quadric_value(const albers::vec10 &Q, const vec3 &x) {
  return Q[0] * x[0] * x[0] + Q[1] * x[1] * x[1] + Q[2] * x[2] * x[2] +
         real(2) * Q[3] * x[0] * x[1] + real(2) * Q[4] * x[0] * x[2] +
         real(2) * Q[5] * x[1] * x[2] + real(2) * Q[6] * x[0] +
         real(2) * Q[7] * x[1] + real(2) * Q[8] * x[2] + Q[9];
}

/// Closest point on the implicit. Walk along the gradient from the spoke,
/// and give up if that walk leaves the one-ring.
template <typename Eval>
inline bool land_along_gradient(const vec3 &dp, Eval &&eval, vec3 &out) {
  if (!dp.allFinite() || !(dp.norm() > real(1e-8)))
    return false;
  const real spoke = dp.norm();
  vec3 x = dp;
  for (int it = 0; it < 8; ++it) {
    const auto fg = eval(x);
    const real f = fg.first;
    const vec3 &g = fg.second;
    const real g2 = g.squaredNorm();
    if (!std::isfinite(f) || !g.allFinite() || !(g2 > real(1e-18)))
      return false;
    if (std::abs(f) <= real(1e-10) * (real(1) + std::sqrt(g2)))
      break;
    const vec3 step = (f / g2) * g;
    if (!step.allFinite() || (x - step - dp).norm() > spoke)
      return false;
    x -= step;
  }
  if (!x.allFinite())
    return false;
  const auto fg = eval(x);
  const real gnorm = fg.second.norm();
  if (!(std::abs(fg.first) <= real(1e-4) * (spoke * gnorm + real(1))))
    return false;
  out = x;
  return (out - dp).norm() <= spoke + real(1e-8);
}

inline bool land_albers_quadric(const albers::vec10 &Q, const vec3 &d,
                                vec3 &out) {
  if (!Q.allFinite())
    return false;
  return land_along_gradient(
      d,
      [&](const vec3 &x) {
        return std::make_pair(albers_quadric_value(Q, x), albers::quadric_grad(Q, x));
      },
      out);
}

inline bool land_albers_sphere(const vec3 &center, real radius, const vec3 &dp,
                               vec3 &out) {
  if (!center.allFinite() || !dp.allFinite() || !std::isfinite(radius) ||
      !(radius > real(1e-8)))
    return false;
  const vec3 v = dp - center;
  const real vn = v.norm();
  if (!(vn > real(1e-8)))
    return false;
  out = center + radius * (v / vn);
  if (!out.allFinite())
    return false;
  const real spoke = dp.norm();
  return (out - dp).norm() <= spoke + real(1e-8);
}

/// Closest point on the tube, \(p_s + r\,\mathrm{normalize}(p - p_s)\).
/// A constructor torus measures \(p_s\) from the nearby spine. A hand-built
/// torus without that frame uses one SDF step, which is the same point when
/// \(\rho - R\) is accurate. A correction longer than the spoke is rejected.
inline bool land_albers_torus(const albers::osculating_torus &tor,
                              const vec3 &dp, vec3 &out) {
  if (!tor.valid || !dp.allFinite())
    return false;
  if (tor.tube) {
    vec3 ps;
    if (!albers::nearest_circle_spine(tor, dp, ps) ||
        !albers::land_from_spine(ps, tor.r, dp, out))
      return false;
  } else {
    const real s = albers::torus_sdf(dp, tor);
    const vec3 g = albers::torus_sdf_grad(dp, tor);
    if (!std::isfinite(s) || !g.allFinite() || !(g.norm() > real(1e-8)))
      return false;
    out = dp - s * g;
  }
  if (!out.allFinite())
    return false;
  const real spoke = dp.norm();
  if (!(spoke > real(1e-8)))
    return false;
  return (out - dp).norm() <= spoke + real(1e-8);
}

enum class minor_tube_fit { none, torus, cylinder };

/// One landing per spoke. `hit` is meaningful where `ok` is set.
struct minor_tube_choice {
  minor_tube_fit fit = minor_tube_fit::none;
  std::vector<char> ok;
  std::vector<vec3> hit;
};

/// Project every spoke onto the spine circle and onto the minor-radius line.
/// The residual is the sum of squared \(p_i - p_{\mathrm{proj}}\) over spokes
/// both landings accept. The smaller residual wins, and every spoke of that
/// fit is kept. Equal residuals keep the circle — a large-\(R\) torus is the
/// cylinder / plane limit (*Points as Tori*), not a drop to the line.
inline minor_tube_choice
choose_minor_tube_landing(const albers::osculating_torus &tor,
                          const std::vector<vec3> &spokes) {
  minor_tube_choice out;
  const size_t n = spokes.size();
  out.ok.assign(n, 0);
  out.hit.assign(n, vec3::Zero());
  if (!tor.valid && !tor.tube)
    return out;

  std::vector<char> tor_ok(n, 0);
  std::vector<char> cyl_ok(n, 0);
  std::vector<vec3> tor_hit(n, vec3::Zero());
  std::vector<vec3> cyl_hit(n, vec3::Zero());
  real tor_l2 = 0;
  real cyl_l2 = 0;
  real reach = 0;
  int both = 0;
  for (size_t i = 0; i < n; ++i) {
    const vec3 &p = spokes[i];
    if (!p.allFinite())
      continue;
    reach = std::max(reach, p.norm());
    vec3 ps = vec3::Zero();
    vec3 hit = vec3::Zero();
    if (tor.valid && albers::nearest_circle_spine(tor, p, ps) &&
        albers::land_from_spine(ps, tor.r, p, hit)) {
      tor_ok[i] = 1;
      tor_hit[i] = hit;
    }
    if (tor.tube && albers::nearest_line_spine(tor, p, ps) &&
        albers::land_from_spine(ps, tor.r, p, hit)) {
      cyl_ok[i] = 1;
      cyl_hit[i] = hit;
    }
    if (tor_ok[i] && cyl_ok[i]) {
      tor_l2 += (tor_hit[i] - p).squaredNorm();
      cyl_l2 += (cyl_hit[i] - p).squaredNorm();
      ++both;
    }
  }

  // Cylinder only when it beats the circle by more than float noise on a
  // large-R torus (circle ≈ line). Ties keep the torus.
  const real tol = (reach > real(0)) ? reach * reach * real(1e-8) : real(0);
  const bool want_cyl =
      tor.tube && (!tor.valid || (both > 0 && cyl_l2 + tol < tor_l2));
  if (want_cyl) {
    out.fit = minor_tube_fit::cylinder;
    out.ok = std::move(cyl_ok);
    out.hit = std::move(cyl_hit);
    return out;
  }
  if (tor.valid) {
    out.fit = minor_tube_fit::torus;
    out.ok = std::move(tor_ok);
    out.hit = std::move(tor_hit);
  }
  return out;
}

/// Closest point on the cyclide, along its gradient. A sample on the surface
/// stays there.
inline bool land_albers_cyclide(const albers::vec14 &Q, const vec3 &dp,
                                vec3 &out) {
  if (!Q.allFinite())
    return false;
  return land_along_gradient(
      dp,
      [&](const vec3 &x) {
        return std::make_pair(albers::eval_darboux(Q, x), albers::darboux_grad(Q, x));
      },
      out);
}

inline real albers_kernel(const vec3 &dp, real l0) {
  if (!(l0 > real(1e-8)))
    return 0;
  return calder::calc_inv_dist(dp, l0, real(3));
}

/// `calder::quadric`, plus a foot pin equal to the sample mass so the vertex
/// stays on the surface. `generic_fit` uses the same relative weight when
/// `foot_normal_w` is 1.
inline bool fit_albers_quadric(const std::vector<albers_fit_sample> &samples,
                               real l0, const vec3 &foot_N, albers::vec10 &Q) {
  if (!foot_N.allFinite() || !(foot_N.norm() > real(1e-8)))
    return false;
  albers::quadric fit;
  real mass = 0;
  int n = 0;
  for (const albers_fit_sample &s : samples) {
    if (!(s.w > real(0)) || !s.dp.allFinite() || !s.N.allFinite())
      continue;
    const real wt = s.w * albers_kernel(s.dp, l0);
    if (!(wt > real(1e-8)))
      continue;
    fit.accumulate(wt, s.dp, s.N);
    mass += wt;
    ++n;
  }
  if (n < 4 || !(mass > real(0)))
    return false;
  fit.accumulate(mass, vec3::Zero(), foot_N.normalized());
  Q = fit.solve();
  return Q.allFinite();
}

/// `calder::constrained_sphere` on the whole ring, plus the same foot pin.
/// The outward-side gate is a Barnes–Hut neighborhood test; on a one-ring it
/// drops spokes and the center flips.
inline bool fit_albers_sphere(const std::vector<albers_fit_sample> &samples,
                              real l0, const vec3 &foot_N, vec3 &center,
                              real &radius) {
  if (!foot_N.allFinite() || !(foot_N.norm() > real(1e-8)))
    return false;
  albers::constrained_sphere fit;
  real mass = 0;
  int n = 0;
  for (const albers_fit_sample &s : samples) {
    if (!(s.w > real(0)) || !s.dp.allFinite() || !s.N.allFinite())
      continue;
    const real wt = s.w * albers_kernel(s.dp, l0);
    if (!(wt > real(1e-8)))
      continue;
    fit.accumulate(wt, s.dp, s.N);
    mass += wt;
    ++n;
  }
  if (n < 4 || !(mass > real(0)))
    return false;
  fit.accumulate(mass, vec3::Zero(), foot_N.normalized());
  // `constrained_sphere::solve` completes the square as if the constant term
  // were −1. The rows do not enforce that, so read the center and radius off
  // the same polynomial: |x − c|² = |c|² − E/A.
  const albers::vec5 coef = fit.A.colPivHouseholderQr().solve(fit.b);
  const real A = coef[0];
  if (!coef.allFinite() || !(std::abs(A) > real(1e-12)))
    return false;
  const vec3 BCD = coef.segment(1, 3);
  center = -BCD / (real(2) * A);
  const real r2 = center.squaredNorm() - coef[4] / A;
  if (!center.allFinite() || !std::isfinite(r2) || !(r2 > real(1e-16)))
    return false;
  radius = std::sqrt(r2);
  const real cap = real(1e4) * std::max(l0, real(1e-8));
  return center.norm() < cap && radius < cap;
}

/// `calder::darboux_cyclide_normal_constrained`, with the foot pin at the
/// full sample mass. The 1e-2 pin used by the global fit does not hold a ring.
inline bool fit_albers_cyclide(const std::vector<albers_fit_sample> &samples,
                               real l0, const vec3 &foot_N, albers::vec14 &Q) {
  if (!foot_N.allFinite() || !(foot_N.norm() > real(1e-8)))
    return false;
  albers::normal_constrained_darboux_cyclide fit;
  real mass = 0;
  int n = 0;
  for (const albers_fit_sample &s : samples) {
    if (!(s.w > real(0)) || !s.dp.allFinite() || !s.N.allFinite())
      continue;
    const real wt = s.w * albers_kernel(s.dp, l0);
    if (!(wt > real(1e-8)))
      continue;
    fit.accumulate(wt, s.dp, s.N);
    mass += wt;
    ++n;
  }
  if (n < 4 || !(mass > real(0)))
    return false;
  fit.accumulate(mass, vec3::Zero(), foot_N.normalized());
  Q = fit.solve();
  return Q.allFinite();
}

struct albers_mask {
  bool quadric = false;
  bool sphere = false;
  bool cyclide = false;
  bool any() const { return quadric || sphere || cyclide; }
};

inline albers_mask albers_mask_of(const std::vector<objective> &objectives) {
  albers_mask m;
  for (const objective &o : objectives) {
    const auto *a = std::get_if<albers_fairing>(&o);
    if (a == nullptr)
      continue;
    switch (a->shape) {
    case albers_shape::quadric:
      m.quadric = true;
      break;
    case albers_shape::sphere:
      m.sphere = true;
      break;
    case albers_shape::darboux_cyclide:
      m.cyclide = true;
      break;
    }
  }
  return m;
}

inline bool objectives_have_monge(const std::vector<objective> &objectives) {
  for (const objective &o : objectives)
    if (std::get_if<monge_fairing>(&o) != nullptr)
      return true;
  return false;
}

} // namespace block
} // namespace hepworth
} // namespace gaudi
