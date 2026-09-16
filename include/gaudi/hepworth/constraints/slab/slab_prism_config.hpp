#ifndef GAUDI_HEPWORTH_SLAB_PRISM_CONFIG_HPP
#define GAUDI_HEPWORTH_SLAB_PRISM_CONFIG_HPP

#include "local_slab_field_fit.hpp"

#include "gaudi/common.h"

namespace gaudi {
namespace hepworth {
namespace block {

enum class field_live_update {
  /// Refit Q_k from current stencil geometry every project().
  refit_every_project = 0,
  /// Refit Q_k every field_refit_period project() calls per foot.
  refit_every_n = 1,
  /// Q_k stays equal to Q_r — ΔD ≡ 0, projection is a no-op (identity tests only).
  live_matches_rest = 2,
};

struct slab_prism_config {
  /// Shape-field bulk divergence row weight (legacy refit path).
  real w_field = 0.0;
  /// Incompressible phi material row weight (f_r ∘ φ_q).
  real w_phi = 0.0;
  /// Local RBF reduced elastic element (fixed-center RBF modes; gauge + ARAP strain).
  real w_rbf = 0.0;
  /// Relative dual-ring twist coupling.
  real w_twist = 1.0;
  /// Top-ring bulk Procrustes + normal-aligned volume pump.
  real w_bulk = 0.0;
  real L_min = 1e-4;
  /// Cap on effective row weight w/L_rest (0 = no cap).
  real w_eff_max = 400.0;
  bool twist_penalize_bend = false;

  /// Bulk Z-step: rest σ-product spring weight.
  real bulk_w_volume = 1.0;
  /// Bulk Z-step: per-axis rest σ spring weight.
  real bulk_w_length = 1.0;
  /// Max singular-value stretch in bulk Z-step (1 = rotation only).
  real bulk_max_stretch = 1.08;
  /// Max pump / spoke target change per project, as a fraction of rest pierce length.
  real bulk_target_step_frac = 0.2;
  /// Pump axis alignment uses live ring pole normal; false uses rest pole normal.
  bool bulk_pump_use_live_normal = true;

  slab_field_family field_family = slab_field_family::darboux;
  real lam_gap = 1.0;
  real field_foot_normal_w = 1.0;
  int field_gd_iters = 8;
  real field_gd_step_scale = 0.25;
  field_live_update field_live = field_live_update::refit_every_project;
  int field_refit_period = 8;

  /// Incompressible φ local NR iterations (0 = skip inner solve).
  int phi_nr_iters = 8;
  real phi_mu = 1.0;
  real phi_lam_gap = 1.0;

  /// RBF ROM: gauge weight μ on (f_r(x+u))².
  real rbf_mu = 1.0;
  /// RBF ROM: spoke ARAP strain weight η.
  real rbf_eta = 1.0;
  real rbf_lam_gap = 1.0;
  int rbf_nr_iters = 8;
  int rbf_n_centers = 4;
  real rbf_sigma_scale = 2.0;
};

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // GAUDI_HEPWORTH_SLAB_PRISM_CONFIG_HPP
