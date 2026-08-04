#ifndef GAUDI_DUCHAMP_MODULES_CCD_STEP_HPP
#define GAUDI_DUCHAMP_MODULES_CCD_STEP_HPP

#include <algorithm>
#include <cmath>
#include <vector>

#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/common.h"
#include "gaudi/duchamp/modules/cfl_neighborhood.hpp"

namespace gaudi {
namespace duchamp {

/// Geometry-only substep sizing. Scale drive to dx_max before estimate (see
/// ccd_scale_field_to_dx_max); dt comes from the scaled probe field.
struct ccd_config {
  /// Max |Δx| per substep. ≤0 ⇒ no peak cap.
  real dx_max = 0.0;
  /// Per-substep ceiling. ≤0 ⇒ uncapped.
  real dt_max = 0.0;
  real dt_min = 1.0e-8;
  int max_substeps = 128;
  real eps = 1.0e-16;
  cfl_neighborhood_config neighborhood{};
};

struct ccd_step_stats {
  int substeps = 0;
  real dt_accum = 0.0;
  real last_dti = 0.0;
  cfl_neighborhood_result last{};
};

struct ccd_field_stats {
  size_t n = 0;
  int n_nonfinite = 0;
  int i_peak = -1;
  real peak = 0.0;
  real rms = 0.0;
};

struct ccd_dt_breakdown {
  real dt_geom = 0.0;
  real dt_dx = 0.0; ///< peak cap; 0 ⇒ inactive / no binding cap stored
  real dt = 0.0;
  real drive_scale = 1.0; ///< dx_max / max|drive| when peak exceeded
  ccd_field_stats drive{};
  ccd_field_stats rod_x{};
};

inline ccd_field_stats ccd_field_stats_vec(const std::vector<vec3> &v) {
  ccd_field_stats s;
  s.n = v.size();
  if (v.empty())
    return s;
  real sum2 = 0.0;
  for (size_t i = 0; i < v.size(); ++i) {
    if (!v[i].array().isFinite().all()) {
      ++s.n_nonfinite;
      continue;
    }
    const real n = v[i].norm();
    sum2 += n * n;
    if (n > s.peak) {
      s.peak = n;
      s.i_peak = static_cast<int>(i);
    }
  }
  s.rms = std::sqrt(sum2 / std::max<size_t>(1, v.size()));
  return s;
}

inline void ccd_diagnose_rod(const asawa::rod::rod &rod, ccd_field_stats *x_out) {
  if (x_out)
    *x_out = ccd_field_stats_vec(rod.x());
}

inline real ccd_peak_norm(const std::vector<vec3> &v) {
  real peak = 0.0;
  for (const vec3 &a : v) {
    const real n = a.norm();
    if (std::isfinite(n))
      peak = std::max(peak, n);
  }
  return peak;
}

/// Scale v in place so max|v| ≤ dx_max. Returns factor applied (≤ 1).
inline real ccd_scale_field_to_dx_max(std::vector<vec3> &v, real dx_max,
                                      real eps = 1.0e-16) {
  if (!(dx_max > 0.0))
    return 1.0;
  const real peak = ccd_peak_norm(v);
  if (!(peak > eps) || !std::isfinite(peak) || peak <= dx_max)
    return 1.0;
  const real alpha = dx_max / peak;
  for (vec3 &a : v)
    a *= alpha;
  return alpha;
}

/// Scale f so max|h²f| ≤ dx_max at reference step h.
inline real ccd_scale_force_to_dx_max(std::vector<vec3> &f, real h_ref,
                                      real dx_max, real eps = 1.0e-16) {
  if (!(dx_max > 0.0) || !(h_ref > 0.0))
    return 1.0;
  const real h2 = h_ref * h_ref;
  real peak_u = 0.0;
  for (const vec3 &a : f) {
    const real n = (h2 * a).norm();
    if (std::isfinite(n))
      peak_u = std::max(peak_u, n);
  }
  if (!(peak_u > eps) || peak_u <= dx_max)
    return 1.0;
  const real alpha = dx_max / peak_u;
  for (vec3 &a : f)
    a *= alpha;
  return alpha;
}

inline real ccd_sanitize_dt(real dt, const ccd_config &cfg) {
  if (std::isfinite(dt))
    return dt;
  if (cfg.dt_max > 0.0)
    return cfg.dt_max;
  if (cfg.dt_min > 0.0)
    return cfg.dt_min;
  return real(0.0);
}

inline real ccd_residual_dt(real dti, real dt_step, real dt_accum) {
  if (!std::isfinite(dti))
    return std::max(real(0.0), dt_step - dt_accum);
  return std::max(real(0.0), std::min(dti, dt_step - dt_accum));
}

/// Displacement u at h_ref=1 (unit probe): safe dt from edge + envelope projector.
inline real ccd_displacement(const asawa::rod::rod &rod,
                             const std::vector<vec3> &u,
                             const ccd_config &cfg,
                             cfl_neighborhood_result *detail = nullptr) {
  if (u.size() != rod.x().size()) {
    cfl_neighborhood_result r;
    r.dt = ccd_sanitize_dt(0.0, cfg);
    if (detail)
      *detail = r;
    return r.dt;
  }
  cfl_neighborhood_config neigh = cfg.neighborhood;
  neigh.h_ref = 1.0;
  if (neigh.dt_max <= 0.0 && cfg.dt_max > 0.0)
    neigh.dt_max = cfg.dt_max;
  const cfl_neighborhood_result r = cfl_estimate_rod_from_u(rod, u, neigh);
  if (detail)
    *detail = r;
  return ccd_sanitize_dt(r.dt, cfg);
}

/// v = Δx per unit time (x += dt·v). u_probe = v.
inline real ccd_velocity(const asawa::rod::rod &rod, const std::vector<vec3> &v,
                         const ccd_config &cfg,
                         cfl_neighborhood_result *detail = nullptr,
                         ccd_dt_breakdown *breakdown = nullptr) {
  if (breakdown) {
    breakdown->drive = ccd_field_stats_vec(v);
    ccd_diagnose_rod(rod, &breakdown->rod_x);
  }

  real dt = ccd_displacement(rod, v, cfg, detail);
  const real dt_geom = dt;

  if (cfg.dt_max > 0.0)
    dt = std::min(dt, cfg.dt_max);
  if (cfg.dt_min > 0.0) {
    if (!(dt > 0.0))
      dt = cfg.dt_min;
    else
      dt = std::max(dt, cfg.dt_min);
  }

  real dt_dx = 0.0;
  if (cfg.dx_max > 0.0) {
    const real peak = ccd_peak_norm(v);
    if (peak > cfg.eps && std::isfinite(peak)) {
      dt_dx = cfg.dx_max / peak;
      dt = std::min(dt, dt_dx); // peak cap wins over dt_min
    }
  }
  dt = ccd_sanitize_dt(dt, cfg);

  if (breakdown) {
    breakdown->dt_geom = dt_geom;
    breakdown->dt_dx = dt_dx;
    breakdown->dt = dt;
  }
  return dt;
}

/// Force integrator: Δx = dt²·f. Probe u = h_ref²·f at h_ref.
inline real ccd_force(const asawa::rod::rod &rod, const std::vector<vec3> &f,
                      real h_ref, const ccd_config &cfg,
                      cfl_neighborhood_result *detail = nullptr) {
  if (f.size() != rod.x().size() || !(h_ref > 0.0) || !std::isfinite(h_ref)) {
    cfl_neighborhood_result r;
    r.dt = ccd_sanitize_dt(0.0, cfg);
    if (detail)
      *detail = r;
    return r.dt;
  }
  const size_t n = f.size();
  std::vector<vec3> u(n);
  const real h2 = h_ref * h_ref;
  for (size_t i = 0; i < n; ++i)
    u[i] = h2 * f[i];

  cfl_neighborhood_config neigh = cfg.neighborhood;
  neigh.h_ref = h_ref;
  if (neigh.dt_max <= 0.0 && cfg.dt_max > 0.0)
    neigh.dt_max = cfg.dt_max;
  const cfl_neighborhood_result r = cfl_estimate_rod_from_u(rod, u, neigh);
  if (detail)
    *detail = r;

  real dt = r.dt;
  if (cfg.dt_max > 0.0)
    dt = std::min(dt, cfg.dt_max);
  if (cfg.dt_min > 0.0) {
    if (!(dt > 0.0))
      dt = cfg.dt_min;
    else
      dt = std::max(dt, cfg.dt_min);
  }
  if (cfg.dx_max > 0.0) {
    const real peak = ccd_peak_norm(f);
    if (peak > cfg.eps && std::isfinite(peak))
      dt = std::min(dt, std::sqrt(cfg.dx_max / peak));
  }
  return ccd_sanitize_dt(dt, cfg);
}

/// Accumulate substeps until target_dt or max_substeps.
/// substep(dt_accum) runs one substep (simulate, estimate, retopo, …) and
/// returns the dti actually used.
template <typename SubstepFn>
ccd_step_stats ccd_step(real target_dt, SubstepFn &&substep,
                        const ccd_config &cfg = {}) {
  ccd_step_stats stats;
  if (!(target_dt > 0.0) || cfg.max_substeps <= 0)
    return stats;

  const real eps = cfg.eps;
  real accum = 0.0;
  while (accum + eps < target_dt && stats.substeps < cfg.max_substeps) {
    const real dti = substep(accum);
    if (!(dti > eps))
      break;
    accum += dti;
    stats.substeps += 1;
    stats.last_dti = dti;
  }
  stats.dt_accum = accum;
  return stats;
}

using ccd_neighborhood_config = cfl_neighborhood_config;
using ccd_neighborhood_result = cfl_neighborhood_result;
inline const char *ccd_limit_term_name(cfl_limit_term t) {
  return cfl_limit_term_name(t);
}

} // namespace duchamp
} // namespace gaudi

#endif
