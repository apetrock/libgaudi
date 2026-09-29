#ifndef GAUDI_DUCHAMP_TP_CONFIG_HPP
#define GAUDI_DUCHAMP_TP_CONFIG_HPP

#include <vector>

#include "gaudi/common.h"
#include "gaudi/duchamp/modules/ccd_step.hpp"
#include "gaudi/duchamp/modules/tp_line_search.hpp"

namespace gaudi {
namespace duchamp {

struct tangent_point_force_config {
  real w = 0.0;
  real l0 = 1.0;
  real p = 6.0;
};

/// Softmax-floor TP force (no Cauchy ℓ₀). Classic density k=(1/R̃)^p.
struct soft_tangent_point_force_config {
  real w = 1.0e-6;
  real p = 6.0;
  real R_min_frac = 1.0;
  real tau_frac = 0.1;
  real R_min = 0.0;
  real tau = 0.0;
};

inline void resolve_soft_tangent_scales(soft_tangent_point_force_config &cfg,
                                        real rod_radius) {
  const real r = std::max(rod_radius, real(1e-6));
  if (!(cfg.R_min > 0.0))
    cfg.R_min = std::max(cfg.R_min_frac, real(1e-6)) * r;
  if (!(cfg.tau > 0.0))
    cfg.tau = std::max(cfg.tau_frac, real(1e-6)) * cfg.R_min;
}

/// Cauchy leaf: K = |N·dp|^p / (|dp|²+ε²)^p with ε = l0 · Cc (already regularized).
struct cauchy_tangent_point_force_config {
  real w = 1.0;
  real l0 = 1.0; ///< ε / Cc
  real p = 6.0;
};

/// Sobolev: f = Cauchy ∇K, then MG/Jacobi Ag=f with A=M+α_B B+α_0 B_0.
struct sobolev_tangent_point_force_config {
  real w = 1.0;   ///< Cauchy leaf weight
  real l0 = 0.01; ///< Cauchy ε / Cc (also metric length regularizer)
  real p = 6.0;   ///< Cauchy / metric power
  real alpha_B = 1.0;
  real alpha_0 = 1.0;
  real mass = 1.0; ///< diagonal M scale on sample masses
  int jacobi_iters = 8; ///< pre/post smooth + coarsest Jacobi budget
  int num_v_cycles = 2;
  int max_levels = 4; ///< incl. L0; ≤0 ⇒ full forest depth
  real bh_eps = 0.25;  ///< box ratio r_x < dist * bh_eps; ≤0 ⇒ leaves-only
  real cone_eps = 1.0; ///< tangent chord r_T < cone_eps; ≤0 ⇒ leaves-only
};

enum class tangent_point_type {
  /// Differentiate E0=|Pdp|^p/|dp|^{2p}, then Reg denoms (l0²).
  regularized,
  /// Softmin-floor composition; exact chain-rule ∇/H.
  soft,
  /// Cauchy K with ε in E; exact ∇K/HK.
  cauchy,
  /// Sobolev: f = Cauchy ∇K, then MG Ag=f (A=M+α_B B+α_0 B_0).
  sobolev
};

struct tangent_point_solver_config {
  /// One or more kernels. The propose field is the sum of their gradients.
  std::vector<tangent_point_type> type{tangent_point_type::soft};
  tangent_point_force_config regularized{};
  soft_tangent_point_force_config soft{};
  cauchy_tangent_point_force_config cauchy{};
  sobolev_tangent_point_force_config sobolev{};
  /// Shared alpha-window settings (all TP types).
  tp_line_search_config line_search{};

  /// Drive conditioning: Ĝ = G / (peak_scale · max|G|). ≤0 disables.
  real peak_scale = 1.0;
  /// Peak |Δx| / lavg per substep when ccd.dx_max ≤ 0.
  real ccd_dx_max_frac = 0.25;
  ccd_config ccd{
      .dx_max = 0.0, ///< 0 ⇒ ccd_dx_max_frac · lavg
      .dt_max = 0.05,
      .dt_min = 1.0e-8,
      .max_substeps = 16,
  };
  ccd_neighborhood_config neighborhood{
      .cfl_frac = 0.4,
      .agg = cfl_agg_mode::hard_min,
      .gate_frac = 0.1,
      .use_envelope = true,
      .chain_skip_hops = 8,
  };
};

inline tangent_point_type
tp_primary_type(const tangent_point_solver_config &cfg) {
  if (cfg.type.empty())
    return tangent_point_type::soft;
  return cfg.type.front();
}

} // namespace duchamp
} // namespace gaudi

#endif // GAUDI_DUCHAMP_TP_CONFIG_HPP
