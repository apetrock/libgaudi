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

enum class sobolev_solver { jacobi, cg, chebyshev };

/// Smoother on rods above the coarse saddle. The leaf factors the KKT
/// system and does not use this choice.
struct sobolev_solver_config {
  sobolev_solver type = sobolev_solver::jacobi;
  int iter = 8;
  /// Chebyshev damps [chebyshev_lower · λ_max, 1.1 λ_max]. 0.3 is the
  /// usual multigrid window: the coarse rod keeps the modes below it.
  real chebyshev_lower = 0.3;
};

/// Sobolev: f = Cauchy ∇K, then multigrid Ag=f with A=α_B B+α_0 B_0.
struct sobolev_tangent_point_force_config {
  real w = 1.0;    ///< energy weight
  real l0 = 0.01;  ///< ε / Cc (energy and metric length regularizer)
  real p0 = 6.0;   ///< unused by the vertex energy; metric order is p1
  real p1 = 6.0;   ///< metric: B is k_0^{p1}, B0 is k_2^{p1+4}
  real tpe_alpha = 3; ///< energy |T×r|^α / |r|^β, reference tpeAlpha
  real tpe_beta = 6;  ///< reference tpeBeta
  real alpha_B = 1.0;
  real alpha_0 = 1.0;
  sobolev_solver_config solver{};
  /// Coarse rod with fewer corners than this is a dense all-pairs matvec
  /// plus Eigen CG to coarse_cg_tol. Larger rods use `solver`.
  /// ≤0 keeps the tree matvec on every rod.
  int dense_cg_threshold = 100;
  real coarse_cg_tol = 1e-8; ///< relative residual for that dense CG solve
  int num_v_cycles = 2;
  int max_levels = 4; ///< incl. L0; ≤0 ⇒ full forest depth
  real bh_eps = 0.25;  ///< box ratio r_x < dist * bh_eps; ≤0 ⇒ leaves-only
  real cone_eps = 1.0; ///< tangent chord r_T < cone_eps; ≤0 ⇒ leaves-only
  /// Coulomb load added to ∇Ê before Ag = f. 0 leaves the load alone.
  /// The vector is projected onto the binormal T × (x_prev − 2x + x_next)
  /// so the in-plane part does not fight the tangent-point gradient.
  /// Positive pulls toward other samples (the sign of `coulomb_force`).
  real coulomb_w = 0;
  real coulomb_p = 3;
  /// Vortex load, same binormal projection as Coulomb. 0 leaves it off.
  /// The kernel is −κ (dp × T), then only the curve binormal survives.
  real vortex_w = 0;
  real vortex_p = 4;
  real vortex_q = 1;
};

/// Harmonic: Cauchy leaf (l0, p0), then ∇(κ ⟨T_i, T_j⟩ Tp)
/// with Gaussian mollifier κ (l1, p1). Lengths are multiples of Cc.
/// p1 = 2 is exp(−|dp|²/l1²).
struct harmonic_tangent_point_force_config {
  real w = 1.0;
  real l0 = 0.1; ///< Cauchy ε / Cc
  real l1 = 4.0; ///< Gaussian length / Cc
  real p0 = 6.0; ///< Cauchy power
  real p1 = 2.0; ///< Gaussian power
};

enum class sobolev_constraint_type {
  edge_length,
  barycenter,
  c1,
  c2,
  c4,
  elastic
};

/// `w` scales a projective row (the Gram uses w²). `slack` is only read on
/// hard rows: 0 stays in the nullspace projector, nonzero stays out of it.
struct sobolev_constraint_spec {
  sobolev_constraint_type type = sobolev_constraint_type::edge_length;
  real w = 1.0;
  real slack = 0.0;
};

enum class tangent_point_type {
  /// Differentiate E0=|Pdp|^p/|dp|^{2p}, then Reg denoms (l0²).
  regularized,
  /// Softmin-floor composition; exact chain-rule ∇/H.
  soft,
  /// Cauchy K with ε in E; exact ∇K/HK.
  cauchy,
  /// Sobolev: f = Cauchy ∇K, then MG Ag=f (A=M+α_B B+α_0 B_0).
  sobolev,
  /// Mollified ∇(w Tp) = w′ Tp + w ∇Tp, w = κ ⟨T_i, T_j⟩.
  harmonic
};

struct tangent_point_solver_config {
  /// One or more kernels. The propose field is the sum of their gradients.
  std::vector<tangent_point_type> type{tangent_point_type::soft};
  tangent_point_force_config regularized{};
  soft_tangent_point_force_config soft{};
  cauchy_tangent_point_force_config cauchy{};
  sobolev_tangent_point_force_config sobolev{};
  harmonic_tangent_point_force_config harmonic{};
  /// Shared alpha-window settings (all TP types).
  tp_line_search_config line_search{};

  /// Drive conditioning: Ĝ = G / (peak_scale · max|G|). ≤0 disables.
  /// Unused when the step is the collision-capped line search.
  real peak_scale = 1.0;
  /// Peak |Δx| / lavg per substep when ccd.dx_max ≤ 0 and the
  /// position block is integrating force. Unused on the velocity path.
  real ccd_dx_max_frac = 0.25;
  /// Step clamps plus the geometric contact query (ccd.neighborhood).
  ccd_config ccd{
      .dx_max = 0.0, ///< 0 ⇒ ccd_dx_max_frac · lavg (force path only)
      .dt_max = 0.05,
      .dt_min = 1.0e-8,
      .max_substeps = 16,
  };
  /// Exact rows. Edge length and barycenter with slack 0 enter the
  /// nullspace projector and the coarse multiplier block.
  std::vector<sobolev_constraint_spec> hard_constraints{};
  /// Least-squares rows. Added as Cᵀ W C on A, with W = w². C1, C2, C4, and
  /// elastic live here. They never enter BBᵀ.
  std::vector<sobolev_constraint_spec> proj_constraints{};
  int constraint_cg = 8;
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
