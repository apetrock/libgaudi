#ifndef GAUDI_DUCHAMP_TP_LINE_SEARCH_HPP
#define GAUDI_DUCHAMP_TP_LINE_SEARCH_HPP

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <vector>

#include "gaudi/common.h"

namespace gaudi {
namespace duchamp {

struct tp_line_search_config {
  /// After normalizing d so peak|d|=1, α has length units.
  real alpha_min = 1.0e-8;
  real alpha_max = 0.0; ///< if > trust probe, used as larger window probe
  real fallback_alpha = 0.0; ///< unused by window search (kept for compat)

  /// Probe / trust: |δ|_∞ ≤ trust_frac_lavg · lavg (big step for window).
  real trust_frac_lavg = 0.25;

  real shrink = 0.5;
  int max_backtracks = 32;

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
  real probe = cfg.trust_frac_lavg * L;
  if (cfg.alpha_max > probe)
    probe = cfg.alpha_max;
  return std::max(probe, cfg.alpha_min);
}

/// Large α₀, sample α,ρ α,ρ²α,… and pick lowest finite E.
inline tp_line_search_result tp_best_energy_line_search(
    real E0, real alpha0, const std::function<real(real)> &eval_E,
    const tp_line_search_config &cfg = {}, real E_tol = 0.0) {
  tp_line_search_result out;
  out.E0 = E0;
  out.alpha_quad = alpha0;
  const real amax = std::max(alpha0, cfg.alpha_min);
  out.alpha = std::clamp(alpha0, cfg.alpha_min, amax);

  real best_alpha = cfg.alpha_min;
  real best_E = E0;
  bool found = false;
  int n_eval = 0;
  int best_k = -1;
  const real E_accept = E0 + std::max(E_tol, real(0.0));

  if (cfg.verbose) {
    std::cout << "[tp best-E] begin E0=" << E0 << " α0=" << amax
              << " E_tol=" << E_tol << " shrink=" << cfg.shrink << std::endl;
  }

  real alpha = out.alpha;
  for (int k = 0; k <= cfg.max_backtracks; ++k) {
    const real E1 = eval_E(alpha);
    ++n_eval;
    out.n_backtracks = k;
    const real dE = E1 - E0;
    bool is_best = false;
    if (std::isfinite(E1) && (!found || E1 < best_E)) {
      best_E = E1;
      best_alpha = alpha;
      found = true;
      best_k = k;
      is_best = true;
    }
    if (cfg.verbose) {
      std::cout << "[tp best-E] k=" << k << " alpha=" << alpha << " E=" << E1
                << " dE=" << dE << " best=" << (is_best ? 1 : 0) << std::endl;
    }
    if (k == cfg.max_backtracks)
      break;
    alpha *= cfg.shrink;
    if (!(alpha > cfg.alpha_min))
      break;
    alpha = std::max(alpha, cfg.alpha_min);
  }

  if (!found || !(best_E < E_accept)) {
    out.used_fallback = true;
    out.alpha = 0.0;
    out.E_accepted = E0;
    out.peak_delta = 0.0;
    if (cfg.verbose) {
      std::cout << "[tp best-E] hold α=0 α0=" << amax << " E0=" << E0
                << " samples=" << n_eval << " (no sample beat E0+tol)"
                << std::endl;
    }
    return out;
  }

  out.alpha = best_alpha;
  out.E_accepted = best_E;
  out.peak_delta = out.alpha;
  if (cfg.verbose) {
    std::cout << "[tp best-E] α*=" << out.alpha << " α0=" << amax
              << " E0=" << E0 << " E*=" << best_E << " dE*=" << (best_E - E0)
              << " samples=" << n_eval << " best_k=" << best_k << std::endl;
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
