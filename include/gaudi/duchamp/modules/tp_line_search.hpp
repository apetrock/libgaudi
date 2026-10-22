#ifndef GAUDI_DUCHAMP_TP_LINE_SEARCH_HPP
#define GAUDI_DUCHAMP_TP_LINE_SEARCH_HPP

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <vector>

#include "gaudi/common.h"
#include "gaudi/duchamp/demo_log.hpp"

namespace gaudi {
namespace duchamp {

struct tp_line_search_config {
  /// Cap on the search start: |δ|_∞ ≤ trust_frac_lavg · lavg.
  /// Collision time can only make the start shorter.
  real trust_frac_lavg = 0.25;
  /// Geometric shrink of Boyd & Vandenberghe Algorithm 9.2 (β in (0, 1)).
  real shrink = 0.5;
  /// How many times to multiply the step by shrink, starting at α0.
  int max_backtracks = 32;
  /// Bisections of the bracket around the best shrink sample.
  int max_refine = 4;

  bool verbose = false;
};

struct tp_line_search_result {
  real alpha = 0.0;
  real alpha_quad = 0.0;
  real g_d = 0.0;
  real H_dd = 0.0;
  real E0 = 0.0;
  real E_accepted = 0.0;
  real E_tp0 = 0.0;
  real E_smooth0 = 0.0;
  real peak_delta = 0.0;
  int n_backtracks = 0;
  bool used_fallback = false;
};

inline real tp_field_dot(const std::vector<vec3> &a, const std::vector<vec3> &b) {
  real s = 0.0;
  const size_t n = std::min(a.size(), b.size());
  for (size_t i = 0; i < n; ++i)
    s += a[i].dot(b[i]);
  return s;
}

inline real tp_field_peak(const std::vector<vec3> &v) {
  real m = 0.0;
  for (const vec3 &a : v) {
    const real n = a.norm();
    if (std::isfinite(n))
      m = std::max(m, n);
  }
  return m;
}

inline real tp_field_rms(const std::vector<vec3> &v) {
  if (v.empty())
    return real(0.0);
  return std::sqrt(tp_field_dot(v, v) / real(v.size()));
}

inline void tp_scale_field(std::vector<vec3> &v, real alpha) {
  for (vec3 &a : v)
    a *= alpha;
}

/// Peak-normalize so |d|_∞ = 1 (α then has length units). Returns prior peak.
inline real tp_peak_normalize(std::vector<vec3> &d, real eps = 1.0e-16) {
  const real peak = tp_field_peak(d);
  if (peak > eps) {
    for (vec3 &a : d)
      a /= peak;
  }
  return peak;
}

inline real tp_probe_length(real lavg, const tp_line_search_config &cfg) {
  const real L = std::max(lavg, real(1e-6));
  return std::max(cfg.trust_frac_lavg * L, real(0));
}

/// Shrink from α0, then bisect the bracket around the best sample.
/// Accept only when that energy is below E0 + tol. A step that never
/// beats E0 is held (α = 0).
/// Repulsive Curves §5.4 cites Boyd & Vandenberghe Algorithm 9.2 for the
/// shrink. max_refine fills the gap between those geometric samples.
inline tp_line_search_result tp_best_energy_line_search(
    real E0, real alpha0, const std::function<real(real)> &eval_E,
    const tp_line_search_config &cfg = {}, real E_tol = 0.0) {
  tp_line_search_result out;
  out.E0 = E0;
  out.alpha_quad = alpha0;
  const real amax = std::max(alpha0, real(0));
  const real E_accept = E0 + std::max(E_tol, real(0.0));
  const real beta = (cfg.shrink > real(0) && cfg.shrink < real(1))
                        ? cfg.shrink
                        : real(0.5);
  const int n_shrink = std::max(0, cfg.max_backtracks);
  const int n_refine = std::max(0, cfg.max_refine);

  if (cfg.verbose) {
    demo_log() << "[tp best-E] begin E0=" << E0 << " α0=" << amax
              << " shrink=" << beta << " backtracks=" << n_shrink
              << " refine=" << n_refine << std::endl;
  }

  int n_eval = 0;
  auto probe = [&](real alpha) -> real {
    const real E1 = eval_E(alpha);
    ++n_eval;
    if (cfg.verbose) {
      demo_log() << "[tp best-E] k=" << (n_eval - 1) << " alpha=" << alpha
                << " E=" << E1 << " dE=" << (E1 - E0) << std::endl;
    }
    return E1;
  };

  real a = amax;
  real best_a = real(0);
  real best_E = E0;
  bool have_best = false;
  real lo = real(0);
  real hi = amax;
  bool have_lo = false;
  bool have_hi = false;
  int geom_evals = 0;

  for (int k = 0; k <= n_shrink; ++k) {
    if (!(a > real(0)) || !std::isfinite(a))
      break;
    const real E = probe(a);
    ++geom_evals;
    const bool better =
        std::isfinite(E) && (!have_best || E < best_E);
    if (better) {
      if (have_best && best_a > a) {
        hi = best_a;
        have_hi = true;
      }
      best_a = a;
      best_E = E;
      have_best = true;
      have_lo = false;
    } else if (have_best && best_E < E_accept && a < best_a) {
      lo = a;
      have_lo = true;
      break;
    }
    const real next = a * beta;
    if (!(next < a))
      break;
    a = next;
  }

  out.n_backtracks = std::max(0, geom_evals - 1);
  if (!have_best || !(best_E < E_accept)) {
    out.used_fallback = true;
    out.alpha = 0.0;
    out.E_accepted = E0;
    out.peak_delta = 0.0;
    if (cfg.verbose) {
      demo_log() << "[tp best-E] hold α=0 α0=" << amax << " E0=" << E0
                << " samples=" << n_eval << " (no sample beat E0+tol)"
                << std::endl;
    }
    return out;
  }

  if (!have_hi)
    hi = std::max(best_a, amax);
  if (!have_lo)
    lo = best_a * beta;
  if (!(lo < best_a))
    lo = real(0);
  if (!(hi > best_a))
    hi = best_a;

  real aL = lo;
  real aR = hi;
  real aM = best_a;
  real EM = best_E;
  for (int r = 0; r < n_refine && aR > aL; ++r) {
    const real gapL = aM - aL;
    const real gapR = aR - aM;
    if (!(gapL > real(0)) && !(gapR > real(0)))
      break;
    const real aC = (gapL >= gapR) ? real(0.5) * (aL + aM)
                                   : real(0.5) * (aM + aR);
    if (!(aC > aL) || !(aC < aR))
      break;
    const real EC = probe(aC);
    const bool c_left = aC < aM;
    if (std::isfinite(EC) && EC < EM) {
      if (c_left)
        aR = aM;
      else
        aL = aM;
      aM = aC;
      EM = EC;
    } else if (c_left) {
      aL = aC;
    } else {
      aR = aC;
    }
  }

  out.alpha = aM;
  out.E_accepted = EM;
  out.peak_delta = out.alpha;
  if (cfg.verbose) {
    demo_log() << "[tp best-E] α*=" << out.alpha << " α0=" << amax
              << " E0=" << E0 << " E*=" << EM << " dE*=" << (EM - E0)
              << " samples=" << n_eval << std::endl;
  }
  return out;
}

inline std::vector<vec3> tp_scaled_direction(const std::vector<vec3> &d,
                                             real alpha) {
  std::vector<vec3> out = d;
  tp_scale_field(out, alpha);
  return out;
}

} // namespace duchamp
} // namespace gaudi

#endif // GAUDI_DUCHAMP_TP_LINE_SEARCH_HPP
