#ifndef __CALDER_MLS_JET_BOOTSTRAP__
#define __CALDER_MLS_JET_BOOTSTRAP__

#include "gaudi/albers/osculating_torus.hpp"
#include "gaudi/albers/quadric.hpp"
#include "gaudi/albers/shape_operator_constrained_fit.hpp"
#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/asawa/shell/shape_operator.hpp"
#include "gaudi/calder/least_squares_fit.hpp"
#include "gaudi/calder/shape_operator_weights.hpp"
#include "gaudi/common.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace gaudi {
namespace calder {

enum class mls_stage2_weight {
  torus_soft_gaussian,
  aniso_W1,
  soft_gaussian_only,
};

/// Ablation ladder for diagnosing bootstrap stages.
enum class mls_bootstrap_ablation {
  /// Single-stage normal-constrained Darboux with plain inv_dist (no convexity).
  plain_harmonic,
  /// two_ring jet → Green double-layer stage-1 → plain Gaussian(σ=R) stage-2.
  jet_green_plain_gaussian,
  /// Legacy: disc-aniso R into soft-convex Gaussian stage-2.
  radii_into_legacy_gaussian,
  // Soft-convex Darboux → foot-jet osculating torus → (Nj·∇T̂)×Gaussian(σ).
  soft_darboux_torus_normal_align,
  /// Soft-convex Darboux → foot-jet torus → Green max(0, dp·∇T̂)/r³.
  soft_darboux_torus_green_dp,
  /// Plain harmonic Darboux stage-1 → foot-jet torus → Green stage-2.
  harmonic_darboux_torus_green_dp,
  /// Plain harmonic stage-1 → foot-jet torus → soft-convex × (Nj·∇T̂) × G(σ=R).
  /// Degenerate torus prior → soft-convex × G (same as adaptive gaussian).
  harmonic_darboux_torus_normal_align,
  /// Plain harmonic stage-1 → radii from |κ|_max → soft-convex × Gaussian(σ).
  /// (Plain G without the half-space gate poisons the jet from the backside.)
  harmonic_plain_gaussian,
  /// Full torus-dist stage-2 (experimental; deferred).
  torus_stage2,
};

struct mls_jet_bootstrap_params {
  mls_bootstrap_ablation ablation =
      mls_bootstrap_ablation::soft_darboux_torus_green_dp;
  asawa::shell::face_curvature_stencil disc_stencil =
      asawa::shell::face_curvature_stencil::two_ring;
  bool use_disc_aniso = false;
  bool use_torus_prior = false;
  albers::torus_sign_mode torus_sign = albers::torus_sign_mode::pat_eq22;
  mls_stage2_weight stage2_weight = mls_stage2_weight::soft_gaussian_only;
  real fit_p = 3.0;
  real fit_w0 = 1.0;
  real fit_w_shape = 0.0;
  real radius_scale = 1.0;
  real aniso_eps = 1e-3;
  real torus_dist_p = 2.0;
};

inline const char *mls_bootstrap_ablation_name(mls_bootstrap_ablation a) {
  switch (a) {
  case mls_bootstrap_ablation::plain_harmonic:
    return "plainHarmonic(inv_dist^p)";
  case mls_bootstrap_ablation::jet_green_plain_gaussian:
    return "twoRingJetGreen→plainGaussian";
  case mls_bootstrap_ablation::radii_into_legacy_gaussian:
    return "anisoR→legacySoftGaussian";
  case mls_bootstrap_ablation::soft_darboux_torus_normal_align:
    return "softDarbouxTorus→(Nj·∇T̂)×G(σ)";
  case mls_bootstrap_ablation::soft_darboux_torus_green_dp:
    return "softDarbouxTorus→Green(dp·∇T)";
  case mls_bootstrap_ablation::harmonic_darboux_torus_green_dp:
    return "harmonicDarbouxTorus→Green(dp·∇T)";
  case mls_bootstrap_ablation::harmonic_darboux_torus_normal_align:
    return "harmonic→soft×(Nj·∇T̂)×G(σ=R)";
  case mls_bootstrap_ablation::harmonic_plain_gaussian:
    return "harmonic→softGaussian(σ=R)";
  case mls_bootstrap_ablation::torus_stage2:
    return "torus_stage2";
  }
  return "unknown";
}

inline mat3 shape_operator_from_quadric(const albers::vec10 &Q) {
  const vec3 g = albers::quadric_grad(Q, vec3::Zero());
  const mat3 H = albers::quadric_hessian(Q);
  mat3 W = mat3::Zero();
  albers::medial_generated::shape_operator_from_GH(g, H, W);
  return W;
}

inline std::vector<asawa::shell::face_height_jet>
local_jets_for_queries(asawa::shell::shell &M, const std::vector<vec3> &p_pov,
                       const std::vector<vec3> &N_pov,
                       const mls_jet_bootstrap_params &params) {
  const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);
  std::vector<asawa::shell::face_height_jet> jets(p_pov.size());

  const bool vertex_aligned =
      (p_pov.size() == x.size());
  if (vertex_aligned) {
    bool same = true;
    for (size_t i = 0; i < p_pov.size(); ++i) {
      if ((p_pov[i] - x[i]).squaredNorm() > 1e-20) {
        same = false;
        break;
      }
    }
    if (same) {
      for (size_t i = 0; i < p_pov.size(); ++i) {
        const vec3 n = (i < N_pov.size()) ? N_pov[i] : vec3::UnitZ();
        jets[i] = asawa::shell::vertex_height_jet_fit(
            M, x, asawa::shell::vert_id(static_cast<int>(i)), n,
            params.disc_stencil);
      }
      return jets;
    }
  }

  for (size_t i = 0; i < p_pov.size(); ++i) {
    int best = 0;
    real best_d2 = (p_pov[i] - x[0]).squaredNorm();
    for (size_t v = 1; v < x.size(); ++v) {
      const real d2 = (p_pov[i] - x[v]).squaredNorm();
      if (d2 < best_d2) {
        best_d2 = d2;
        best = static_cast<int>(v);
      }
    }
    const vec3 n = (i < N_pov.size()) ? N_pov[i] : vec3::UnitZ();
    jets[i] = asawa::shell::vertex_height_jet_fit(
        M, x, asawa::shell::vert_id(best), n, params.disc_stencil);
  }
  return jets;
}

/// Stage 1: quadric MLS with jet-Green double-layer weights. Returns W1.
inline std::vector<mat3>
quadric_jet_W_green(asawa::shell::shell &M, const std::vector<vec3> &p_pov,
                    const std::vector<vec3> &N_pov, real l0,
                    const std::vector<asawa::shell::face_height_jet> &jets,
                    const mls_jet_bootstrap_params &params) {
  const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);
  std::vector<real> areas = asawa::shell::face_areas(M, x);
  std::vector<vec3> Ns = asawa::shell::face_normals(M, x);
  for (int i = 0; i < static_cast<int>(Ns.size()); ++i) {
    Ns[i] = areas[i] * Ns[i];
  }

  auto weight = [&jets, &N_pov](
                    int i, int /*j*/,
                    const std::vector<calder::datum::ptr> & /*data*/,
                    shell_bundle::Sum_Type::Node_Type /*node_type*/,
                    const vec3 &dp, const vec3 & /*Ni*/, const vec3 & /*Nj*/,
                    real l0_, real /*p*/) -> real {
    if (i < 0 || static_cast<size_t>(i) >= jets.size()) {
      return 0.0;
    }
    const vec3 Ni =
        (static_cast<size_t>(i) < N_pov.size()) ? N_pov[static_cast<size_t>(i)]
                                                 : vec3::UnitZ();
    return jet_green_double_layer_weight(dp, l0_, jets[static_cast<size_t>(i)],
                                         Ni);
  };

  const std::vector<albers::vec10> Q =
      generic_fit<albers::quadric, shell_bundle>(
          M, Ns, p_pov, N_pov, l0, params.fit_p, weight, params.fit_w0);

  std::vector<mat3> W1(Q.size(), mat3::Zero());
  for (size_t i = 0; i < Q.size(); ++i) {
    W1[i] = shape_operator_from_quadric(Q[i]);
  }
  return W1;
}

/// Legacy stage-1: aniso from disc W0 (kept for A/B).
inline std::vector<mat3>
quadric_jet_W(asawa::shell::shell &M, const std::vector<vec3> &p_pov,
              const std::vector<vec3> &N_pov, real l0,
              const std::vector<mat3> &W0,
              const mls_jet_bootstrap_params &params) {
  const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);
  std::vector<real> areas = asawa::shell::face_areas(M, x);
  std::vector<vec3> Ns = asawa::shell::face_normals(M, x);
  for (int i = 0; i < static_cast<int>(Ns.size()); ++i) {
    Ns[i] = areas[i] * Ns[i];
  }

  std::vector<albers::principal_curvature_frame> frames(W0.size());
  for (size_t i = 0; i < W0.size(); ++i) {
    const vec3 n = (i < N_pov.size()) ? N_pov[i] : vec3::UnitZ();
    frames[i] = albers::principal_frame_from_shape_operator(W0[i], n);
  }

  const real aniso_eps = params.aniso_eps;
  const bool use_aniso = params.use_disc_aniso;
  auto weight = [frames, use_aniso, aniso_eps](
                    int i, int /*j*/,
                    const std::vector<calder::datum::ptr> & /*data*/,
                    shell_bundle::Sum_Type::Node_Type /*node_type*/,
                    const vec3 &dp, const vec3 & /*Ni*/, const vec3 & /*Nj*/,
                    real l0_, real p_) -> real {
    if (!use_aniso || i < 0 || static_cast<size_t>(i) >= frames.size()) {
      return calc_inv_dist(dp, l0_, p_);
    }
    return aniso_inv_dist_weight(dp, l0_, p_, frames[static_cast<size_t>(i)],
                                 aniso_eps);
  };

  const std::vector<albers::vec10> Q =
      generic_fit<albers::quadric, shell_bundle>(
          M, Ns, p_pov, N_pov, l0, params.fit_p, weight, params.fit_w0);

  std::vector<mat3> W1(Q.size(), mat3::Zero());
  for (size_t i = 0; i < Q.size(); ++i) {
    W1[i] = shape_operator_from_quadric(Q[i]);
  }
  return W1;
}

inline std::vector<real>
radii_from_shape_operators(const std::vector<mat3> &W1, real l0) {
  const real r_min = 0.25 * std::max(l0, real(1e-12));
  const real r_max = 64.0 * std::max(l0, real(1e-12));
  std::vector<real> radii(W1.size(), l0);
  for (size_t i = 0; i < W1.size(); ++i) {
    real R = radius_from_shape_operator(W1[i], l0);
    if (!std::isfinite(R) || R < 1e-12) {
      R = l0;
    }
    radii[i] = std::clamp(R, r_min, r_max);
  }
  return radii;
}

inline std::vector<mat3>
disc_shape_operators_for_queries(asawa::shell::shell &M,
                                 const std::vector<vec3> &p_pov,
                                 const mls_jet_bootstrap_params &params) {
  const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);
  const std::vector<mat3> W_vert =
      asawa::shell::vertex_shape_operators(M, x, params.disc_stencil);

  std::vector<mat3> W0(p_pov.size(), mat3::Zero());
  if (p_pov.size() == W_vert.size()) {
    bool same = true;
    for (size_t i = 0; i < p_pov.size(); ++i) {
      if ((p_pov[i] - x[i]).squaredNorm() > 1e-20) {
        same = false;
        break;
      }
    }
    if (same) {
      return W_vert;
    }
  }
  for (size_t i = 0; i < p_pov.size(); ++i) {
    int best = 0;
    real best_d2 = (p_pov[i] - x[0]).squaredNorm();
    for (size_t v = 1; v < x.size(); ++v) {
      const real d2 = (p_pov[i] - x[v]).squaredNorm();
      if (d2 < best_d2) {
        best_d2 = d2;
        best = static_cast<int>(v);
      }
    }
    W0[i] = W_vert[static_cast<size_t>(best)];
  }
  return W0;
}

/// Default: two_ring jet → Green stage-1 → plain Gaussian stage-2.
inline std::vector<albers::vec14>
darboux_fit_bootstrapped(asawa::shell::shell &M, const std::vector<vec3> &p_pov,
                         const std::vector<vec3> &N_pov, real l0,
                         const mls_jet_bootstrap_params &params_in =
                             mls_jet_bootstrap_params{}) {
  mls_jet_bootstrap_params params = params_in;

  if (params.ablation == mls_bootstrap_ablation::plain_harmonic) {
    // Single-stage: w = 1 / (|dp|^p + l0^p), p=fit_p (demo default 3).
    // No soft-convex / Green half-space gate.
    return darboux_cyclide_normal_constrained(M, p_pov, N_pov, l0, params.fit_p,
                                              params.fit_w0);
  }

  if (params.ablation == mls_bootstrap_ablation::harmonic_plain_gaussian) {
    // Stage-1 harmonic jet → R from |κ|_max → stage-2 soft-convex × G(σ=R).
    // Clamp to the mesh feature range, NOT edge-length multiples: with
    // l0≈avg_edge, 8*l0≈0.04 while true |κ|^{-1} is O(0.1–1) on a unit mesh,
    // so an edge clamp silently kills the cyclide prior.
    const std::vector<albers::vec14> Q0 = darboux_cyclide_normal_constrained(
        M, p_pov, N_pov, l0, params.fit_p, params.fit_w0);
    const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);
    vec3 lo = x.empty() ? vec3::Zero() : x.front();
    vec3 hi = lo;
    for (const vec3 &p : x) {
      lo = lo.cwiseMin(p);
      hi = hi.cwiseMax(p);
    }
    const real scene = std::max(real(0.5) * (hi - lo).norm(), l0);
    const real r_min = 2.0 * std::max(l0, real(1e-12));
    const real r_max = 0.5 * scene;
    std::vector<real> radii(Q0.size(), l0);
    for (size_t i = 0; i < Q0.size(); ++i) {
      real r = radius_from_cyclide_max_curvature(Q0[i], l0);
      if (!std::isfinite(r) || r < 1e-12) {
        r = l0;
      }
      radii[i] = std::clamp(r, r_min, r_max);
    }
    return darboux_cyclide_normal_constrained_adaptive_gaussian(
        M, p_pov, N_pov, l0, radii, params.fit_p, params.fit_w0,
        params.radius_scale);
  }

  if (params.ablation == mls_bootstrap_ablation::jet_green_plain_gaussian) {
    params.disc_stencil = asawa::shell::face_curvature_stencil::two_ring;
    params.use_disc_aniso = false;
    params.use_torus_prior = false;

    const auto jets = local_jets_for_queries(M, p_pov, N_pov, params);
    const std::vector<mat3> W1 =
        quadric_jet_W_green(M, p_pov, N_pov, l0, jets, params);
    const std::vector<real> radii = radii_from_shape_operators(W1, l0);
    return darboux_cyclide_normal_constrained_adaptive_gaussian_plain(
        M, p_pov, N_pov, l0, radii, params.fit_p, params.fit_w0,
        params.radius_scale);
  }

  if (params.ablation == mls_bootstrap_ablation::radii_into_legacy_gaussian) {
    params.use_disc_aniso = true;
    params.use_torus_prior = false;
    const std::vector<mat3> W0 =
        disc_shape_operators_for_queries(M, p_pov, params);
    const std::vector<mat3> W1 =
        quadric_jet_W(M, p_pov, N_pov, l0, W0, params);
    const std::vector<real> radii = radii_from_shape_operators(W1, l0);
    return darboux_cyclide_normal_constrained_adaptive_gaussian(
        M, p_pov, N_pov, l0, radii, params.fit_p, params.fit_w0,
        params.radius_scale);
  }

  if (params.ablation ==
          mls_bootstrap_ablation::soft_darboux_torus_normal_align ||
      params.ablation == mls_bootstrap_ablation::soft_darboux_torus_green_dp ||
      params.ablation ==
          mls_bootstrap_ablation::harmonic_darboux_torus_green_dp ||
      params.ablation ==
          mls_bootstrap_ablation::harmonic_darboux_torus_normal_align) {
    // Stage-1 jet → shape op → foot-local osculating torus.
    // soft_*: soft-convex inv_dist gate; harmonic_*: plain inv_dist^p.
    const bool harmonic_stage1 =
        (params.ablation ==
             mls_bootstrap_ablation::harmonic_darboux_torus_green_dp ||
         params.ablation ==
             mls_bootstrap_ablation::harmonic_darboux_torus_normal_align);
    const std::vector<albers::vec14> Q0 =
        harmonic_stage1
            ? darboux_cyclide_normal_constrained(M, p_pov, N_pov, l0,
                                                 params.fit_p, params.fit_w0)
            : darboux_cyclide_normal_constrained_convexity(
                  M, p_pov, N_pov, l0, params.fit_p, params.fit_w0);

    std::vector<albers::osculating_torus> tori(Q0.size());
    std::vector<real> sigmas(Q0.size(), l0);
    const std::vector<vec3> &x_bbox = asawa::const_get_vec_data(M, 0);
    vec3 lo = x_bbox.empty() ? vec3::Zero() : x_bbox.front();
    vec3 hi = lo;
    for (const vec3 &p : x_bbox) {
      lo = lo.cwiseMin(p);
      hi = hi.cwiseMax(p);
    }
    const real scene = std::max(real(0.5) * (hi - lo).norm(), l0);
    const real r_min = 2.0 * std::max(l0, real(1e-12));
    const real r_max = 0.5 * scene;
    for (size_t i = 0; i < Q0.size(); ++i) {
      const vec3 n = (i < N_pov.size()) ? N_pov[i] : vec3::UnitZ();
      const mat3 W = albers::shape_operator_at(Q0[i], vec3::Zero());
      tori[i] = albers::osculating_torus_from_shape_operator(
          vec3::Zero(), n, W, params.torus_sign, l0);
      // Feature scale from cyclide |κ|_max (same as adaptive-G path).
      // Torus supplies the align prior only — do not swap σ for raw tube r.
      real r_feat = radius_from_cyclide_max_curvature(Q0[i], l0);
      if (!std::isfinite(r_feat) || r_feat < 1e-12) {
        r_feat = l0;
      }
      r_feat = std::clamp(r_feat, r_min, r_max);
      sigmas[i] = std::max(params.radius_scale * r_feat, real(1e-12));
    }

    const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);
    std::vector<real> areas = asawa::shell::face_areas(M, x);
    std::vector<vec3> Ns = asawa::shell::face_normals(M, x);
    for (int i = 0; i < static_cast<int>(Ns.size()); ++i) {
      Ns[i] = areas[i] * Ns[i];
    }

    const real fit_p = params.fit_p;
    const bool use_green_dp =
        (params.ablation ==
             mls_bootstrap_ablation::soft_darboux_torus_green_dp ||
         params.ablation ==
             mls_bootstrap_ablation::harmonic_darboux_torus_green_dp);
    const bool use_align_gaussian =
        (params.ablation ==
             mls_bootstrap_ablation::soft_darboux_torus_normal_align ||
         params.ablation ==
             mls_bootstrap_ablation::harmonic_darboux_torus_normal_align);
    auto weight = [&tori, &sigmas, fit_p, use_green_dp, use_align_gaussian](
                      int i, int /*j*/,
                      const std::vector<calder::datum::ptr> & /*data*/,
                      shell_bundle::Sum_Type::Node_Type /*node_type*/,
                      const vec3 &dp, const vec3 & /*Ni*/, const vec3 &Nj,
                      real l0_, real /*p*/) -> real {
      if (i < 0 || static_cast<size_t>(i) >= tori.size()) {
        return 0.0;
      }
      const auto &T = tori[static_cast<size_t>(i)];
      if (use_green_dp) {
        return torus_green_dp_weight(dp, l0_, fit_p, Nj, T);
      }
      if (use_align_gaussian) {
        const real sigma =
            (static_cast<size_t>(i) < sigmas.size())
                ? sigmas[static_cast<size_t>(i)]
                : l0_;
        // w = soft-convex × (Nj·∇T̂) × G(σ=R); degenerate → soft×G.
        return soft_torus_align_gaussian_weight(dp, l0_, Nj, T, sigma);
      }
      return torus_inv_normal_align_weight(dp, l0_, fit_p, Nj, T);
    };

    return generic_fit<albers::normal_constrained_darboux_cyclide,
                       shell_bundle>(M, Ns, p_pov, N_pov, l0, params.fit_p,
                                     weight, params.fit_w0);
  }

  // Experimental torus-dist path.
  const std::vector<mat3> W0 =
      disc_shape_operators_for_queries(M, p_pov, params);
  const std::vector<mat3> W1 =
      quadric_jet_W(M, p_pov, N_pov, l0, W0, params);
  const std::vector<real> radii = radii_from_shape_operators(W1, l0);

  std::vector<real> sigmas(W1.size(), l0);
  std::vector<albers::osculating_torus> tori(W1.size());
  std::vector<albers::principal_curvature_frame> frames1(W1.size());
  for (size_t i = 0; i < W1.size(); ++i) {
    const vec3 n = (i < N_pov.size()) ? N_pov[i] : vec3::UnitZ();
    frames1[i] = albers::principal_frame_from_shape_operator(W1[i], n);
    sigmas[i] = std::max(params.radius_scale * radii[i], real(1e-12));
    if (params.use_torus_prior) {
      tori[i] = albers::osculating_torus_from_shape_operator(
          p_pov[i], n, W1[i], params.torus_sign, l0);
    }
  }

  const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);
  std::vector<real> areas = asawa::shell::face_areas(M, x);
  std::vector<vec3> Ns = asawa::shell::face_normals(M, x);
  for (int i = 0; i < static_cast<int>(Ns.size()); ++i) {
    Ns[i] = areas[i] * Ns[i];
  }

  const mls_stage2_weight wmode = params.stage2_weight;
  const bool use_torus = params.use_torus_prior;
  const real torus_p = params.torus_dist_p;
  const real aniso_eps = params.aniso_eps;

  auto weight = [=](int i, int /*j*/,
                    const std::vector<calder::datum::ptr> & /*data*/,
                    shell_bundle::Sum_Type::Node_Type /*node_type*/,
                    const vec3 &dp, const vec3 & /*Ni*/, const vec3 &Nj,
                    real /*l0*/, real /*p*/) -> real {
    if (i < 0 || static_cast<size_t>(i) >= sigmas.size()) {
      return 0.0;
    }
    const real sigma = sigmas[static_cast<size_t>(i)];
    if (wmode == mls_stage2_weight::aniso_W1) {
      return aniso_inv_dist_weight(dp, sigma, 2.0,
                                   frames1[static_cast<size_t>(i)], aniso_eps);
    }
    if (wmode == mls_stage2_weight::soft_gaussian_only || !use_torus ||
        !tori[static_cast<size_t>(i)].valid) {
      return calc_gaussian(dp, sigma);
    }
    return torus_soft_gaussian_weight(dp, Nj, sigma,
                                      p_pov[static_cast<size_t>(i)],
                                      tori[static_cast<size_t>(i)], torus_p);
  };

  return generic_fit<albers::normal_constrained_darboux_cyclide, shell_bundle>(
      M, Ns, p_pov, N_pov, l0, params.fit_p, weight, params.fit_w0);
}

} // namespace calder
} // namespace gaudi

#endif // __CALDER_MLS_JET_BOOTSTRAP__
