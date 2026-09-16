#ifndef GAUDI_DUCHAMP_TP_CONFIG_HPP
#define GAUDI_DUCHAMP_TP_CONFIG_HPP

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

struct harmonic_tangent_point_force_config {
  real w = 1.0;
  /// Leaf = regularized TP density k; l0 scaled by Cc unless absolute set later.
  real l0 = 0.1;
  real l1 = 4.0; ///< Gaussian mollifier length / Cc
  real p0 = 6.0;
  real p1 = 2.0; ///< unused (Gaussian)
  real smooth_l0 = 0.0; ///< if >0, absolute mollifier length
};

inline void resolve_harmonic_tangent_scales(harmonic_tangent_point_force_config &cfg,
                                            real /*rod_radius*/, real Cc) {
  const real eps = std::max(Cc, real(1e-6));
  if (!(cfg.smooth_l0 > 0.0))
    cfg.smooth_l0 = std::max(cfg.l1, real(1e-6)) * eps;
}

/// Cauchy leaf: K = |N·dp|^p / (|dp|²+ε²)^p with ε = l0 · Cc (already regularized).
struct cauchy_tangent_point_force_config {
  real w = 1.0;
  real l0 = 1.0; ///< ε / Cc
  real p = 6.0;
};

/// Cauchy + λ circle-agreement: E = ½(Ki+Kj) + λ EΔ(|ci−cj|).
struct cauchy_agree_tangent_point_force_config {
  real w = 1.0;
  real l0 = 1.0; ///< ε / Cc
  real p = 6.0;
  real lam = 0.5; ///< weight on scale-free circle center mismatch EΔ
};

/// Compose: drive = −(w0·slide + w1·cauchy) = w0·t − w1·∇K (positive w0,w1).
struct surface_slide_tangent_point_force_config {
  real w0 = 1.0; ///< circle-slide weight
  real w1 = 1.0; ///< Cauchy core weight
  real l0 = 1.0; ///< Cauchy ε / Cc
  real p = 6.0;
};

enum class tangent_point_type {
  /// Differentiate E0=|Pdp|^p/|dp|^{2p}, then Reg denoms (l0²).
  regularized,
  /// Softmin-floor composition; exact chain-rule ∇/H.
  soft,
  /// E=M(T)=κ*T with T=k; exact leaf ∇T/HT then mollify.
  harmonic,
  /// Cauchy K with ε in E; exact ∇K/HK.
  cauchy,
  /// Cauchy + λ same-circle agreement (EΔ).
  cauchy_agree,
  /// w0·slide + w1·cauchy (composed from the two integrators).
  surface_slide
};

struct tangent_point_solver_config {
  tangent_point_type type = tangent_point_type::soft;
  tangent_point_force_config regularized{};
  soft_tangent_point_force_config soft{};
  harmonic_tangent_point_force_config harmonic{};
  cauchy_tangent_point_force_config cauchy{};
  cauchy_agree_tangent_point_force_config cauchy_agree{};
  surface_slide_tangent_point_force_config surface_slide{};
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

} // namespace duchamp
} // namespace gaudi

#endif // GAUDI_DUCHAMP_TP_CONFIG_HPP
