#ifndef GAUDI_DUCHAMP_TP_ENERGIES_HPP
#define GAUDI_DUCHAMP_TP_ENERGIES_HPP

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

inline real tp_energy_harmonic(asawa::rod::rod &rod,
                               const asawa::rod::dynamic &dynamic,
                               const harmonic_tangent_point_force_config &cfg_in) {
  harmonic_tangent_point_force_config cfg = cfg_in;
  resolve_harmonic_tangent_scales(cfg, rod._r, dynamic._Cc);
  const real eps = dynamic._Cc;
  const std::vector<vec3> &x = rod.x();
  const std::vector<real> l = rod.l0();
  const std::vector<vec3> T = rod.N2c();
  const real e = tp_sum_energy(calder::tangent_point_energy_harmonic(
      rod, x, l, T, cfg.l0 * eps, cfg.p0, cfg.smooth_l0, cfg.p1));
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
                             const cauchy_tangent_point_force_config &cfg) {
  const real eps = dynamic._Cc;
  const std::vector<vec3> &x = rod.x();
  const std::vector<real> l = rod.l0();
  const std::vector<vec3> T = rod.N2c();
  const real e = tp_sum_energy(calder::tangent_point_energy_cauchy(
      rod, x, l, T, cfg.l0 * eps, cfg.p));
  return cfg.w * e;
}

inline real tp_energy_cauchy_agree(
    asawa::rod::rod &rod, const asawa::rod::dynamic &dynamic,
    const cauchy_agree_tangent_point_force_config &cfg) {
  const real eps = dynamic._Cc;
  const std::vector<vec3> &x = rod.x();
  const std::vector<real> l = rod.l0();
  const std::vector<vec3> T = rod.N2c();
  const real e = tp_sum_energy(calder::tangent_point_energy_cauchy_agree(
      rod, x, l, T, cfg.l0 * eps, cfg.p, cfg.lam));
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

/// ∇E (energy gradient, before force negate), weighted.
inline std::vector<vec3>
tp_grad_energy(asawa::rod::rod &rod, const asawa::rod::dynamic &dynamic,
               const tangent_point_solver_config &cfg) {
  const std::vector<vec3> &x = rod.x();
  const std::vector<real> l = rod.l0();
  const std::vector<vec3> T = rod.N2c();
  const real eps = dynamic._Cc;

  switch (cfg.type) {
  case tangent_point_type::soft: {
    soft_tangent_point_force_config soft = cfg.soft;
    resolve_soft_tangent_scales(soft, rod._r);
    return soft_tp_grad_energy(rod, soft.R_min, soft.tau, soft.p, soft.w);
  }
  case tangent_point_type::harmonic: {
    harmonic_tangent_point_force_config harm = cfg.harmonic;
    resolve_harmonic_tangent_scales(harm, rod._r, eps);
    auto g = calder::tangent_point_harmonic_gradient(
        rod, x, l, T, harm.l0 * eps, harm.p0, harm.smooth_l0, harm.p1);
    if (harm.w != 1.0) {
      for (auto &a : g)
        a *= harm.w;
    }
    return g;
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
  case tangent_point_type::cauchy_agree: {
    auto g = calder::tangent_point_gradient_cauchy_agree(
        rod, x, l, T, cfg.cauchy_agree.l0 * eps, cfg.cauchy_agree.p,
        cfg.cauchy_agree.lam);
    if (cfg.cauchy_agree.w != 1.0) {
      for (auto &a : g)
        a *= cfg.cauchy_agree.w;
    }
    return g;
  }
  case tangent_point_type::surface_slide: {
    // Leaf = w0·slide + w1·cauchy (same composition as compute_* before −).
    const auto &ss = cfg.surface_slide;
    std::vector<vec3> g;
    if (ss.w1 != real(0.0)) {
      g = calder::tangent_point_gradient_cauchy(rod, x, l, T, ss.l0 * eps,
                                               ss.p);
      if (ss.w1 != real(1.0)) {
        for (auto &a : g)
          a *= ss.w1;
      }
    } else {
      g.assign(x.size(), vec3::Zero());
    }
    if (ss.w0 != real(0.0)) {
      auto slide = calder::tangent_point_surface_slide(rod, x, l, T);
      const real s = ss.w0;
      for (size_t i = 0; i < g.size(); ++i)
        g[i] += s * slide[i];
    }
    return g;
  }
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

inline real tp_energy(asawa::rod::rod &rod, const asawa::rod::dynamic &dynamic,
                      const tangent_point_solver_config &cfg) {
  switch (cfg.type) {
  case tangent_point_type::soft:
    return tp_energy_soft(rod, cfg.soft);
  case tangent_point_type::harmonic:
    return tp_energy_harmonic(rod, dynamic, cfg.harmonic);
  case tangent_point_type::cauchy:
    return tp_energy_cauchy(rod, dynamic, cfg.cauchy);
  case tangent_point_type::cauchy_agree:
    return tp_energy_cauchy_agree(rod, dynamic, cfg.cauchy_agree);
  case tangent_point_type::surface_slide:
    // Monitor with Cauchy core (w1); slide has no scalar energy.
    return tp_energy_cauchy(rod, dynamic,
                            {.w = cfg.surface_slide.w1,
                             .l0 = cfg.surface_slide.l0,
                             .p = cfg.surface_slide.p});
  case tangent_point_type::regularized:
  default:
    return tp_energy_regularized(rod, dynamic, cfg.regularized);
  }
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
