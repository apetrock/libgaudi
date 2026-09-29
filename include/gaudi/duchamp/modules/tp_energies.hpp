#ifndef GAUDI_DUCHAMP_TP_ENERGIES_HPP
#define GAUDI_DUCHAMP_TP_ENERGIES_HPP

#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>
#include <vector>

#include "gaudi/asawa/rod/dynamic.hpp"
#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/calder/tangent_point_integrators.hpp"
#include "gaudi/common.h"
#include "gaudi/duchamp/modules/tp_config.hpp"
#include "gaudi/duchamp/modules/tp_line_search.hpp"
#include "gaudi/duchamp/modules/tp_sobolev_jacobi.hpp"

namespace gaudi {
namespace duchamp {

inline real tp_sum_energy(const std::vector<real> &e) {
  return std::accumulate(e.begin(), e.end(), real(0.0));
}

inline real tp_energy_regularized(asawa::rod::rod &rod,
                                  const asawa::rod::dynamic &dynamic,
                                  const tangent_point_force_config &cfg) {
  const real eps = dynamic._Cc;
  const std::vector<vec3> &x = rod.x();
  const std::vector<real> l = rod.l0();
  const std::vector<vec3> T = rod.N2c();
  const real e = tp_sum_energy(
      calder::tangent_point_energy(rod, x, l, T, cfg.l0 * eps, cfg.p));
  return cfg.w * e;
}

inline real tp_energy_soft(asawa::rod::rod &rod,
                           const soft_tangent_point_force_config &cfg_in) {
  soft_tangent_point_force_config cfg = cfg_in;
  resolve_soft_tangent_scales(cfg, rod._r);
  const std::vector<vec3> &x = rod.x();
  const std::vector<real> l = rod.l0();
  const std::vector<vec3> T = rod.N2c();
  const real e = tp_sum_energy(calder::tangent_point_energy_soft(
      rod, x, l, T, cfg.R_min, cfg.tau, cfg.p));
  return cfg.w * e;
}

inline real tp_energy_cauchy(asawa::rod::rod &rod,
                             const asawa::rod::dynamic &dynamic,
                             const cauchy_tangent_point_force_config &cfg,
                             real bh_eps = 0.25, real cone_eps = -1) {
  const real eps = dynamic._Cc;
  const std::vector<vec3> &x = rod.x();
  const std::vector<real> l = rod.l0();
  const std::vector<vec3> T = rod.N2c();
  const real e = tp_sum_energy(calder::tangent_point_energy_cauchy(
      rod, x, l, T, cfg.l0 * eps, cfg.p, bh_eps, cone_eps));
  return cfg.w * e;
}

/// ∇E from soft TP (before force negate), weighted.
inline std::vector<vec3> soft_tp_grad_energy(asawa::rod::rod &rod, real R_min,
                                             real tau, real p, real w) {
  const std::vector<vec3> &x = rod.x();
  const std::vector<real> l = rod.l0();
  const std::vector<vec3> T = rod.N2c();
  auto g = calder::tangent_point_gradient_soft(rod, x, l, T, R_min, tau, p,
                                               /*hess_alpha=*/0.0);
  if (w != 1.0) {
    for (auto &a : g)
      a *= w;
  }
  return g;
}

/// ∇E of one kernel (energy gradient, before force negate), weighted.
inline std::vector<vec3>
tp_grad_energy_one(asawa::rod::rod &rod, const asawa::rod::dynamic &dynamic,
                   const tangent_point_solver_config &cfg,
                   tangent_point_type type) {
  const std::vector<vec3> &x = rod.x();
  const std::vector<real> l = rod.l0();
  const std::vector<vec3> T = rod.N2c();
  const real eps = dynamic._Cc;

  switch (type) {
  case tangent_point_type::soft: {
    soft_tangent_point_force_config soft = cfg.soft;
    resolve_soft_tangent_scales(soft, rod._r);
    return soft_tp_grad_energy(rod, soft.R_min, soft.tau, soft.p, soft.w);
  }
  case tangent_point_type::cauchy: {
    auto g = calder::tangent_point_gradient_cauchy(
        rod, x, l, T, cfg.cauchy.l0 * eps, cfg.cauchy.p);
    if (cfg.cauchy.w != 1.0) {
      for (auto &a : g)
        a *= cfg.cauchy.w;
    }
    return g;
  }
  case tangent_point_type::sobolev:
    return sobolev_tp_gradient(rod, dynamic, cfg.sobolev);
  case tangent_point_type::regularized:
  default: {
    auto g = calder::tangent_point_gradient(rod, x, l, T,
                                            cfg.regularized.l0 * eps,
                                            cfg.regularized.p);
    if (cfg.regularized.w != 1.0) {
      for (auto &a : g)
        a *= cfg.regularized.w;
    }
    return g;
  }
  }
}

inline real tp_energy_one(asawa::rod::rod &rod,
                          const asawa::rod::dynamic &dynamic,
                          const tangent_point_solver_config &cfg,
                          tangent_point_type type) {
  switch (type) {
  case tangent_point_type::soft:
    return tp_energy_soft(rod, cfg.soft);
  case tangent_point_type::cauchy:
    return tp_energy_cauchy(rod, dynamic, cfg.cauchy);
  case tangent_point_type::sobolev: {
    const real l0 = cfg.sobolev.l0 * dynamic._Cc;
    const auto ev = sobolev_detail::eval_paper_tp(
        rod, l0, cfg.sobolev.p, cfg.sobolev.bh_eps, cfg.sobolev.cone_eps);
    return cfg.sobolev.w * ev.energy;
  }
  case tangent_point_type::regularized:
  default:
    return tp_energy_regularized(rod, dynamic, cfg.regularized);
  }
}

/// ∇E summed over cfg.type. A one-element list is that kernel alone.
inline std::vector<vec3>
tp_grad_energy(asawa::rod::rod &rod, const asawa::rod::dynamic &dynamic,
               const tangent_point_solver_config &cfg) {
  if (cfg.type.empty())
    return tp_grad_energy_one(rod, dynamic, cfg, tangent_point_type::soft);
  std::vector<vec3> sum =
      tp_grad_energy_one(rod, dynamic, cfg, cfg.type.front());
  for (size_t k = 1; k < cfg.type.size(); ++k) {
    const std::vector<vec3> g =
        tp_grad_energy_one(rod, dynamic, cfg, cfg.type[k]);
    const size_t n = std::min(sum.size(), g.size());
    for (size_t i = 0; i < n; ++i)
      sum[i] += g[i];
  }
  return sum;
}

/// Energy summed over cfg.type, matching the summed gradient.
inline real tp_energy(asawa::rod::rod &rod, const asawa::rod::dynamic &dynamic,
                      const tangent_point_solver_config &cfg) {
  if (cfg.type.empty())
    return tp_energy_one(rod, dynamic, cfg, tangent_point_type::soft);
  real e = 0.0;
  for (tangent_point_type type : cfg.type)
    e += tp_energy_one(rod, dynamic, cfg, type);
  return e;
}

struct tp_curvature_propose {
  tp_line_search_result ls{};
  std::vector<vec3> displacement; ///< α_accepted · d
  std::vector<vec3> direction;    ///< search direction d (peak-norm)
  real propose_ms = 0.0;
};

/// Naive gradient descent: d = −g, peak-norm, alpha window on E_tp.
inline tp_curvature_propose
tp_propose_curvature_step(asawa::rod::rod &rod,
                          const asawa::rod::dynamic &dynamic,
                          const tangent_point_solver_config &cfg) {
  using clock = std::chrono::steady_clock;
  const auto t0 = clock::now();

  tp_curvature_propose out;
  const tp_line_search_config &ls_cfg = cfg.line_search;

  const std::vector<vec3> g = tp_grad_energy(rod, dynamic, cfg);
  const size_t n = g.size();
  out.direction.resize(n);
  for (size_t i = 0; i < n; ++i)
    out.direction[i] = -g[i];
  tp_peak_normalize(out.direction);

  const std::vector<vec3> x0 = rod.x();
  const real E0 = tp_energy(rod, dynamic, cfg);
  out.ls.E0 = E0;
  out.ls.E_tp0 = E0;
  out.ls.g_d = tp_field_dot(g, out.direction);

  auto eval_E = [&](real alpha) -> real {
    for (size_t i = 0; i < n; ++i)
      rod.x()[i] = x0[i] + alpha * out.direction[i];
    rod._update_frames();
    const real E = tp_energy(rod, dynamic, cfg);
    for (size_t i = 0; i < n; ++i)
      rod.x()[i] = x0[i];
    rod._update_frames();
    return E;
  };

  const real probe = tp_probe_length(rod.lavg(), ls_cfg);
  out.ls = tp_best_energy_line_search(E0, probe, eval_E, ls_cfg);
  out.ls.E_tp0 = E0;
  out.ls.g_d = tp_field_dot(g, out.direction);
  out.displacement = tp_scaled_direction(out.direction, out.ls.alpha);

  out.propose_ms =
      real(std::chrono::duration<double, std::milli>(clock::now() - t0).count());
  return out;
}

} // namespace duchamp
} // namespace gaudi

#endif // GAUDI_DUCHAMP_TP_ENERGIES_HPP
