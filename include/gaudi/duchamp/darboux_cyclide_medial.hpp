#ifndef __GAUDI_DUCHAMP_DARBOUX_CYCLIDE_MEDIAL_HPP__
#define __GAUDI_DUCHAMP_DARBOUX_CYCLIDE_MEDIAL_HPP__

#include "gaudi/albers/darboux_cyclide.hpp"
#include "gaudi/asawa/datums.hpp"
#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/calder/least_squares_fit.hpp"
#include "gaudi/calder/mls_jet_bootstrap.hpp"
#include "gaudi/common.h"
#include "gaudi/geometry_logger.hpp"
#include "gaudi/kusama/cyclide_jet_smooth.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

namespace gaudi {
namespace duchamp {

enum class medial_axis_display {
  MedialAxis,
  SmoothedFitHessian,
};

struct cyclide_medial_params {
  real l0_scale = 0.1;
  real fit_p = 2.0;
  // Relative foot tip-in for MLS (quadric stage-1 and Darboux stage-2):
  //   w_foot = fit_w0 * Σ w_MLS   (1.0 ≈ match neighborhood mass)
  real fit_w0 = 1.0;
  real radius_scale = 1.0; // Gaussian σ = radius_scale * R
  // Stage-1 R = 1/|κ|_max from soft-convex MLS + foot tip-in.
  // Flip to quadric to A/B bandwidth without touching stage-2.
  calder::stage1_radius_model stage1_radius =
      calder::stage1_radius_model::quadric;
  // When true, use disc→aniso-quadric→torus-weighted Darboux bootstrap.
  bool use_mls_bootstrap = false;
  calder::mls_jet_bootstrap_params bootstrap;
  real normal_l0 = 0.35;
  real min_normal_alignment = -0.25;
  real medial_smooth_scale = 0.0;
  real medial_smooth_blend = 0.0;
  int max_iters = 120;
  real tol = 1e-8;
  real max_travel_scale = 12.0;
  real newton_step_scale = 1000.0;
  int probe_vertex = -1;
  kusama::cyclide_jet_smooth_params smooth;
  // Temporarily off while investigating divergent medial sites.
  bool enable_cyclide_smooth = false;
  medial_axis_display display = medial_axis_display::MedialAxis;
  //medial_axis_display display = medial_axis_display::SmoothedFitHessian;
  real hessian_frame_scale = 3.0;
  real hessian_line_radius = 0.004;
};

struct cyclide_medial_candidate {
  int vertex = -1;
  vec3 center_world = vec3::Zero();
  vec3 center_local = vec3::Zero();
  vec3 surface_world = vec3::Zero();
  vec3 surface_local = vec3::Zero();
  vec3 cylinder_axis_world = vec3::Zero();
  std::vector<vec3> trace_world;
  real travel = 0.0;
  real residual = 0.0;
  real surface_residual = 0.0;
  int iterations = 0;
  int clamped_steps = 0;
  albers::darboux_ridge_failure failure =
      albers::darboux_ridge_failure::not_converged;
  albers::darboux_ridge_failure projection_failure =
      albers::darboux_ridge_failure::none;
  bool projection_attempted = false;
  bool projection_converged = false;
  bool accepted = false;
};

struct cyclide_medial_stats {
  int total = 0;
  int accepted = 0;
  int rejected = 0;
  real avg_travel = 0.0;
  real avg_residual = 0.0;
};

inline bool cyclide_medial_candidate_valid(
    const albers::darboux_ridge_estimate &est) {
  return est.accepted && est.center.allFinite() && est.travel > 1e-12 &&
         std::isfinite(est.residual);
}

inline void smooth_pipe_fit_medial_axis(
    asawa::shell::shell &M, std::vector<cyclide_medial_candidate> &candidates,
    real smooth_scale, real blend) {
  if (smooth_scale <= 0.0 || candidates.empty()) {
    return;
  }

  const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);
  const real avg_len = asawa::shell::avg_length(M, x);
  const real l0 = std::max(smooth_scale * avg_len, real(1e-12));
  const real alpha = std::clamp(blend, real(0.0), real(1.0));

  std::vector<vec3> centers(x.size(), vec3::Zero());
  std::vector<vec3> axes(x.size(), vec3::Zero());
  std::vector<real> valid(x.size(), 0.0);
  for (const cyclide_medial_candidate &cand : candidates) {
    if (!cand.accepted || cand.vertex < 0 ||
        cand.vertex >= static_cast<int>(x.size())) {
      continue;
    }
    const size_t vi = static_cast<size_t>(cand.vertex);
    centers[vi] = cand.center_world;
    if (cand.cylinder_axis_world.squaredNorm() > 1e-12) {
      axes[vi] = cand.cylinder_axis_world.normalized();
    }
    valid[vi] = 1.0;
  }

  std::vector<vec3> axis_faces(M.face_count(), vec3::Zero());
  for (auto fi : M.get_face_range()) {
    vec3 ref = vec3::Zero();
    vec3 axis_sum = vec3::Zero();
    int axis_count = 0;
    M.const_for_each_face(asawa::shell::face_id(fi),
                          [&](asawa::shell::CornerId c0,
                              const asawa::shell::shell &shell) {
                            const vec3 axis = axes[shell.vert(c0)];
                            if (axis.squaredNorm() <= 1e-12) {
                              return;
                            }
                            if (ref.squaredNorm() <= 1e-12) {
                              ref = axis.normalized();
                            }
                            axis_sum += axis.dot(ref) < 0.0 ? -axis : axis;
                            ++axis_count;
                          });
    if (axis_count > 0 && axis_sum.squaredNorm() > 1e-12) {
      axis_faces[fi] = axis_sum.normalized();
    }
  }

  const std::vector<vec3> center_faces =
      asawa::shell::vert_to_face<vec3>(M, x, centers);
  const std::vector<real> valid_faces =
      asawa::shell::vert_to_face<real>(M, x, valid);
  const std::vector<vec3> smooth_centers =
      calder::mls_avg<vec3>(M, center_faces, x, l0, 2.0);
  const std::vector<mat3> smooth_axis_frames =
      calder::gaussian_covariant_vector_frame(M, axis_faces, x, l0);
  const std::vector<real> smooth_valid =
      calder::mls_avg<real>(M, valid_faces, x, l0, 2.0);

  for (cyclide_medial_candidate &cand : candidates) {
    if (!cand.accepted || cand.vertex < 0 ||
        cand.vertex >= static_cast<int>(x.size())) {
      continue;
    }
    const size_t vi = static_cast<size_t>(cand.vertex);
    if (smooth_valid[vi] <= 1e-8) {
      continue;
    }

    const vec3 center = smooth_centers[vi] / smooth_valid[vi];
    if (center.allFinite()) {
      cand.center_world = (1.0 - alpha) * cand.center_world + alpha * center;
      cand.center_local = cand.center_world - cand.surface_world;
      cand.travel = cand.center_local.norm();
    }

    if (smooth_axis_frames[vi].allFinite()) {
      vec3 axis = smooth_axis_frames[vi].col(2);
      if (!axis.allFinite() || axis.norm() <= 1e-12) {
        continue;
      }
      if (axis.dot(cand.cylinder_axis_world) < 0.0) {
        axis *= -1.0;
      }
      const vec3 blended_axis =
          (1.0 - alpha) * cand.cylinder_axis_world + alpha * axis.normalized();
      if (blended_axis.norm() > 1e-12) {
        cand.cylinder_axis_world = blended_axis.normalized();
      }
    }
  }
}

inline void dump_cyclide_medial_probe(const albers::vec14 &Q, int vertex,
                                      const vec3 &local_start,
                                      const vec3 &zero_dir_in, int max_iters,
                                      real tol, real max_travel,
                                      real max_newton_step) {
  std::cerr << "cyclide medial probe vertex=" << vertex
            << " local_start=" << local_start.transpose()
            << " D_start=" << albers::eval_darboux(Q, local_start)
            << std::endl;

  if (zero_dir_in.norm() < 1e-12) {
    std::cerr << "  zero: invalid direction" << std::endl;
    return;
  }

  const vec3 zero_dir = zero_dir_in.normalized();
  real t_zero = 0.0;
  vec3 x_zero = local_start;
  bool zero_ok = false;
  for (int iter = 0; iter < std::max(1, max_iters); ++iter) {
    x_zero = local_start + t_zero * zero_dir;
    const real D = albers::eval_darboux(Q, x_zero);
    const vec3 g = albers::darboux_grad(Q, x_zero);
    const real denom = g.dot(zero_dir);
    std::cerr << "  zero iter=" << iter << " t=" << t_zero << " D=" << D
              << " |g|=" << g.norm() << " dDdt=" << denom
              << " x=" << x_zero.transpose() << std::endl;
    if (!x_zero.allFinite() || !std::isfinite(D) || !g.allFinite()) {
      std::cerr << "  zero fail=nonfinite" << std::endl;
      return;
    }
    if (std::abs(D) <= std::max(tol, std::sqrt(tol))) {
      zero_ok = true;
      break;
    }
    if (!std::isfinite(denom) || std::abs(denom) < 1e-12) {
      std::cerr << "  zero fail=small_denominator" << std::endl;
      return;
    }
    const real dt = -D / denom;
    std::cerr << "    zero dt=" << dt << std::endl;
    t_zero += dt;
    if (!std::isfinite(dt) || std::abs(t_zero) > max_travel) {
      std::cerr << "  zero fail=max_travel/nonfinite t=" << t_zero
                << std::endl;
      return;
    }
  }

  if (!zero_ok) {
    std::cerr << "  zero fail=not_converged" << std::endl;
    return;
  }

  const vec3 g_zero = albers::darboux_grad(Q, x_zero);
  if (!g_zero.allFinite() || g_zero.norm() < 1e-12) {
    std::cerr << "  medial fail=invalid zero gradient |g|=" << g_zero.norm()
              << std::endl;
    return;
  }

  const vec3 medial_dir = -g_zero.normalized();
  real t_medial = 0.0;
  std::cerr << "  medial start x_zero=" << x_zero.transpose()
            << " D_zero=" << albers::eval_darboux(Q, x_zero)
            << " dir=" << medial_dir.transpose() << std::endl;
  for (int iter = 0; iter < std::max(1, max_iters); ++iter) {
    const vec3 x_medial = x_zero + t_medial * medial_dir;
    const real D = albers::eval_darboux(Q, x_medial);
    const real ridge =
        albers::darboux_directional_deriv(Q, x_medial, medial_dir);
    const real denom =
        albers::darboux_directional_second_deriv(Q, x_medial, medial_dir);
    std::cerr << "  medial iter=" << iter << " t=" << t_medial
              << " D=" << D << " ridge=" << ridge << " d2=" << denom
              << " x=" << x_medial.transpose() << std::endl;
    if (!x_medial.allFinite() || !std::isfinite(ridge) ||
        !std::isfinite(denom)) {
      std::cerr << "  medial fail=nonfinite" << std::endl;
      return;
    }
    if (std::abs(ridge) <= tol) {
      std::cerr << "  medial converged travel=" << t_medial << std::endl;
      return;
    }
    if (std::abs(denom) < 1e-12) {
      std::cerr << "  medial fail=small_denominator" << std::endl;
      return;
    }
    const real raw_dt = -ridge / denom;
    real dt = raw_dt;
    if (std::isfinite(max_newton_step) && max_newton_step > 0.0 &&
        std::abs(dt) > max_newton_step) {
      dt = std::copysign(max_newton_step, dt);
    }
    const real t_next = std::max(t_medial + dt, real(0.0));
    std::cerr << "    medial raw_dt=" << raw_dt << " dt=" << dt
              << " max_step=" << max_newton_step << " t_next=" << t_next
              << (t_medial + dt < 0.0 ? " clamped" : "") << std::endl;
    t_medial = t_next;
    if (!std::isfinite(raw_dt) || t_medial > max_travel) {
      std::cerr << "  medial fail=max_travel/nonfinite t=" << t_medial
                << std::endl;
      return;
    }
  }
  std::cerr << "  medial fail=not_converged final_t=" << t_medial
            << std::endl;
}

inline std::vector<cyclide_medial_candidate>
compute_single_fit_cyclide_medial_candidates(
    asawa::shell::shell &M, int fit_vertex,
    const cyclide_medial_params &params,
    cyclide_medial_stats *stats_out = nullptr,
    albers::vec14 *fit_out = nullptr) {
  const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);
  const std::vector<vec3> vertex_normals =
      asawa::shell::vertex_normals(M, x);
  const real avg_len = asawa::shell::avg_length(M, x);
  const real l0 = std::max(params.l0_scale * avg_len, real(1e-12));
  auto verts = M.get_vert_range();
  if (verts.empty()) {
    if (stats_out != nullptr) {
      *stats_out = cyclide_medial_stats{};
    }
    return {};
  }

  fit_vertex = std::clamp(fit_vertex, 0, static_cast<int>(x.size()) - 1);
  const std::vector<vec3> p_fit = {x[fit_vertex]};
  const std::vector<vec3> n_fit = {vertex_normals[fit_vertex]};
  const std::vector<albers::vec14> fits =
      calder::darboux_cyclide_tangent_plane(M, p_fit, n_fit, l0, params.fit_p);
  const albers::vec14 Q = fits.empty() ? albers::vec14::Zero() : fits[0];
  if (fit_out != nullptr) {
    *fit_out = Q;
  }

  const vec3 fit_origin = x[fit_vertex];
  const real max_travel = params.max_travel_scale * avg_len;
  const real max_newton_step =
      std::max(params.newton_step_scale * l0, real(1e-12));
  std::vector<cyclide_medial_candidate> out;
  cyclide_medial_stats stats;

  real travel_sum = 0.0;
  real residual_sum = 0.0;
  real grad_alignment_sum = 0.0;
  real projection_turn_sum = 0.0;
  real min_grad_alignment = 1.0;
  real max_projection_turn = 0.0;
  int grad_alignment_samples = 0;
  int projection_turn_samples = 0;
  int projection_turn_gt_60 = 0;
  int projection_attempted = 0;
  int projection_failed = 0;
  int projection_fail_max_travel = 0;
  int projection_fail_not_converged = 0;
  int projection_fail_other = 0;
  int reject_small_denominator = 0;
  int reject_nonfinite = 0;
  int reject_max_travel = 0;
  int reject_not_converged = 0;
  int reject_not_converged_clamped = 0;
  int reject_other = 0;
  real local_not_converged_residual_sum = 0.0;
  real local_not_converged_residual_max = 0.0;
  real local_not_converged_travel_sum = 0.0;
  real local_not_converged_travel_max = 0.0;
  int local_not_converged_samples = 0;
  real not_converged_residual_sum = 0.0;
  real not_converged_residual_max = 0.0;
  real not_converged_travel_sum = 0.0;
  real not_converged_travel_max = 0.0;
  int not_converged_samples = 0;

  for (int k = 0; k < static_cast<int>(verts.size()); ++k) {
    const int vi = static_cast<int>(verts[static_cast<size_t>(k)]);
    cyclide_medial_candidate cand;
    cand.vertex = vi;

    const vec3 local_start = x[vi] - fit_origin;
    const vec3 inward = -vertex_normals[vi].normalized();
    const real D_start = albers::eval_darboux(Q, local_start);
    const vec3 g_start = albers::darboux_grad(Q, local_start);
    if (g_start.allFinite() && g_start.norm() > 1e-12) {
      vec3 fit_normal = g_start.normalized();
      const vec3 mesh_normal = vertex_normals[vi].normalized();
      if (fit_normal.dot(mesh_normal) < 0.0) {
        fit_normal *= -1.0;
      }
      const real align = fit_normal.dot(mesh_normal);
      grad_alignment_sum += align;
      min_grad_alignment = std::min(min_grad_alignment, align);
      ++grad_alignment_samples;

      if (D_start > params.tol) {
        const vec3 first_step =
            -(D_start / g_start.squaredNorm()) * g_start;
        if (first_step.allFinite() && first_step.norm() > 1e-12) {
          const real c = std::clamp(first_step.normalized().dot(inward),
                                    real(-1.0), real(1.0));
          const real angle = std::acos(c) * real(180.0 / M_PI);
          projection_turn_sum += angle;
          max_projection_turn = std::max(max_projection_turn, angle);
          ++projection_turn_samples;
          if (angle > 60.0) {
            ++projection_turn_gt_60;
          }
        }
      }
    }
    std::vector<vec3> trace_local;
    const albers::darboux_ridge_estimate est = albers::estimate_center_ridge(
        Q, local_start, inward, params.max_iters, params.tol, &trace_local,
        max_travel);

    cand.center_world = fit_origin + est.center;
    cand.center_local = est.center - local_start;
    cand.trace_world.reserve(trace_local.size());
    for (const vec3 &p_local : trace_local) {
      cand.trace_world.push_back(fit_origin + p_local);
    }
    cand.travel = est.travel;
    cand.residual = est.residual;
    cand.surface_residual = est.surface_residual;
    cand.iterations = est.iterations;
    cand.clamped_steps = est.clamped_steps;
    cand.failure = est.failure;
    cand.projection_failure = est.projection_failure;
    cand.projection_attempted = est.projection_attempted;
    cand.projection_converged = est.projection_converged;
    cand.accepted = cyclide_medial_candidate_valid(est);

    if (cand.projection_attempted) {
      ++projection_attempted;
      if (!cand.projection_converged) {
        ++projection_failed;
        switch (cand.projection_failure) {
        case albers::darboux_ridge_failure::max_travel:
          ++projection_fail_max_travel;
          break;
        case albers::darboux_ridge_failure::not_converged:
          ++projection_fail_not_converged;
          break;
        default:
          ++projection_fail_other;
          break;
        }
      }
    }
    if (cand.accepted) {
      ++stats.accepted;
      travel_sum += cand.travel;
      residual_sum += cand.residual;
    } else {
      ++stats.rejected;
      switch (cand.failure) {
      case albers::darboux_ridge_failure::small_denominator:
        ++reject_small_denominator;
        break;
      case albers::darboux_ridge_failure::nonfinite:
        ++reject_nonfinite;
        break;
      case albers::darboux_ridge_failure::max_travel:
        ++reject_max_travel;
        break;
      case albers::darboux_ridge_failure::not_converged:
        ++reject_not_converged;
        ++not_converged_samples;
        not_converged_residual_sum += cand.residual;
        not_converged_residual_max =
            std::max(not_converged_residual_max, cand.residual);
        not_converged_travel_sum += cand.travel;
        not_converged_travel_max =
            std::max(not_converged_travel_max, cand.travel);
        if (cand.clamped_steps > 0) {
          ++reject_not_converged_clamped;
        }
        break;
      default:
        ++reject_other;
        break;
      }
    }
    out.push_back(cand);
  }

  if (stats.accepted > 0) {
    stats.avg_travel = travel_sum / real(stats.accepted);
    stats.avg_residual = residual_sum / real(stats.accepted);
  }
  stats.total = static_cast<int>(out.size());
  if (stats_out != nullptr) {
    *stats_out = stats;
  }
  std::cout << "cyclide medial newton diagnostics:"
            << " grad_align_avg="
            << (grad_alignment_samples > 0
                    ? grad_alignment_sum / real(grad_alignment_samples)
                    : 0.0)
            << " grad_align_min=" << min_grad_alignment
            << " projection_turn_avg="
            << (projection_turn_samples > 0
                    ? projection_turn_sum / real(projection_turn_samples)
                    : 0.0)
            << " projection_turn_max=" << max_projection_turn
            << " projection_turn_gt_60=" << projection_turn_gt_60
            << " projection_attempted=" << projection_attempted
            << " projection_failed=" << projection_failed
            << " projection_fail_max_travel=" << projection_fail_max_travel
            << " projection_fail_not_converged="
            << projection_fail_not_converged
            << " projection_fail_other=" << projection_fail_other
            << " reject_small_denominator=" << reject_small_denominator
            << " reject_nonfinite=" << reject_nonfinite
            << " reject_max_travel=" << reject_max_travel
            << " reject_not_converged=" << reject_not_converged
            << " reject_not_converged_clamped="
            << reject_not_converged_clamped
            << " not_converged_residual_avg="
            << (not_converged_samples > 0
                    ? not_converged_residual_sum /
                          real(not_converged_samples)
                    : 0.0)
            << " not_converged_residual_max=" << not_converged_residual_max
            << " not_converged_travel_avg="
            << (not_converged_samples > 0
                    ? not_converged_travel_sum / real(not_converged_samples)
                    : 0.0)
            << " not_converged_travel_max=" << not_converged_travel_max
            << " l0=" << l0 << " max_step=" << max_newton_step
            << " max_travel=" << max_travel
            << " reject_other=" << reject_other << std::endl;
  return out;
}

inline std::vector<cyclide_medial_candidate>
compute_single_fit_cyclide_zero_surface_candidates(
    asawa::shell::shell &M, int fit_vertex,
    const cyclide_medial_params &params,
    cyclide_medial_stats *stats_out = nullptr,
    albers::vec14 *fit_out = nullptr) {
  const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);
  const std::vector<vec3> vertex_normals =
      asawa::shell::vertex_normals(M, x);
  const real avg_len = asawa::shell::avg_length(M, x);
  const real l0 = std::max(params.l0_scale * avg_len, real(1e-12));
  auto verts = M.get_vert_range();
  if (verts.empty()) {
    if (stats_out != nullptr) {
      *stats_out = cyclide_medial_stats{};
    }
    return {};
  }

  fit_vertex = std::clamp(fit_vertex, 0, static_cast<int>(x.size()) - 1);
  const std::vector<vec3> p_fit = {x[fit_vertex]};
  const std::vector<vec3> n_fit = {vertex_normals[fit_vertex]};
  const std::vector<albers::vec14> fits =
      calder::darboux_cyclide_tangent_plane(M, p_fit, n_fit, l0, params.fit_p);
  const albers::vec14 Q = fits.empty() ? albers::vec14::Zero() : fits[0];
  if (fit_out != nullptr) {
    *fit_out = Q;
  }

  const vec3 fit_origin = x[fit_vertex];
  const real max_travel = params.max_travel_scale * avg_len;
  const real max_newton_step =
      std::max(params.newton_step_scale * l0, real(1e-12));
  std::vector<cyclide_medial_candidate> out;
  cyclide_medial_stats stats;

  real travel_sum = 0.0;
  real residual_sum = 0.0;
  real first_step_inward_angle_sum = 0.0;
  real first_step_inward_angle_max = 0.0;
  real first_step_root_angle_sum = 0.0;
  real first_step_root_angle_max = 0.0;
  int first_step_samples = 0;
  int first_step_inward_gt_60 = 0;
  int first_step_root_gt_60 = 0;
  int reject_nonfinite = 0;
  int reject_max_travel = 0;
  int reject_not_converged = 0;
  int reject_other = 0;
  int longest_medial_vertex = -1;
  real longest_medial_length = 0.0;
  int medial_negative_step_candidates = 0;
  int medial_negative_steps = 0;
  int medial_backward_clamp_candidates = 0;
  int medial_backward_clamps = 0;

  for (int k = 0; k < static_cast<int>(verts.size()); ++k) {
    const int vi = static_cast<int>(verts[static_cast<size_t>(k)]);
    cyclide_medial_candidate cand;
    cand.vertex = vi;

    const vec3 local_start = x[vi] - fit_origin;
    const vec3 inward = -vertex_normals[vi].normalized();
    std::vector<vec3> trace_local;
    trace_local.push_back(local_start);

    vec3 x_surface = local_start;
    int root_iters = 0;
    albers::darboux_ridge_failure failure =
        albers::darboux_ridge_failure::not_converged;
    const bool ok = albers::darboux_surface_point_line_newton(
        Q, local_start, inward, params.max_iters, params.tol, &x_surface,
        &root_iters, &trace_local, max_travel, &failure);

    if (trace_local.size() > 1) {
      const vec3 first_step = trace_local[1] - trace_local[0];
      const real first_step_norm = first_step.norm();
      const vec3 g_start = albers::darboux_grad(Q, local_start);
      const real d_start = albers::eval_darboux(Q, local_start);
      const real g2 = g_start.squaredNorm();
      if (first_step.allFinite() && first_step_norm > 1e-12 &&
          g_start.allFinite() && g2 > 1e-24) {
        const vec3 step_dir = first_step / first_step_norm;
        const real inward_cos =
            std::clamp(step_dir.dot(inward), real(-1.0), real(1.0));
        const real inward_angle = std::acos(inward_cos) * real(180.0 / M_PI);
        first_step_inward_angle_sum += inward_angle;
        first_step_inward_angle_max =
            std::max(first_step_inward_angle_max, inward_angle);
        if (inward_angle > 60.0) {
          ++first_step_inward_gt_60;
        }

        const vec3 root_step = -(d_start / g2) * g_start;
        const real root_step_norm = root_step.norm();
        if (root_step.allFinite() && root_step_norm > 1e-12) {
          const real root_cos = std::clamp(
              step_dir.dot(root_step / root_step_norm), real(-1.0), real(1.0));
          const real root_angle = std::acos(root_cos) * real(180.0 / M_PI);
          first_step_root_angle_sum += root_angle;
          first_step_root_angle_max =
              std::max(first_step_root_angle_max, root_angle);
          if (root_angle > 60.0) {
            ++first_step_root_gt_60;
          }
        }
        ++first_step_samples;
      }
    }

    cand.surface_world = fit_origin + x_surface;
    cand.surface_local = x_surface - local_start;
    cand.surface_residual = std::abs(albers::eval_darboux(Q, x_surface));
    cand.projection_attempted = true;
    cand.projection_converged = ok;
    cand.projection_failure = failure;

    std::vector<vec3> medial_trace_local;
    albers::darboux_ridge_estimate medial;
    if (ok) {
      const vec3 surface_grad = albers::darboux_grad(Q, x_surface);
      if (surface_grad.allFinite() && surface_grad.norm() > 1e-12) {
        const vec3 medial_dir = -surface_grad.normalized();
        medial = albers::estimate_center_ridge(
            Q, x_surface, medial_dir, params.max_iters, params.tol,
            &medial_trace_local, max_travel, max_newton_step);
      } else {
        medial.center = x_surface;
        medial.failure = albers::darboux_ridge_failure::invalid_direction;
      }
      cand.center_world = fit_origin + medial.center;
      cand.center_local = medial.center - local_start;
      cand.travel = cand.center_local.norm();
      cand.residual = medial.residual;
      cand.iterations = root_iters + medial.iterations;
      cand.clamped_steps = medial.clamped_steps;
      cand.failure = medial.failure;
      cand.accepted = cyclide_medial_candidate_valid(medial);
      if (medial.negative_steps > 0) {
        ++medial_negative_step_candidates;
        medial_negative_steps += medial.negative_steps;
      }
      if (medial.backward_clamps > 0) {
        ++medial_backward_clamp_candidates;
        medial_backward_clamps += medial.backward_clamps;
      }
    } else {
      cand.center_world = cand.surface_world;
      cand.center_local = cand.surface_local;
      cand.travel = cand.center_local.norm();
      cand.residual = cand.surface_residual;
      cand.iterations = root_iters;
      cand.failure = failure;
      cand.accepted = false;
    }

    cand.trace_world.reserve(trace_local.size() + medial_trace_local.size());
    for (const vec3 &p_local : trace_local) {
      cand.trace_world.push_back(fit_origin + p_local);
    }
    for (const vec3 &p_local : medial_trace_local) {
      cand.trace_world.push_back(fit_origin + p_local);
    }

    const real medial_length = (cand.center_world - cand.surface_world).norm();
    if (cand.projection_converged && cand.center_world.allFinite() &&
        std::isfinite(medial_length) && medial_length > longest_medial_length) {
      longest_medial_length = medial_length;
      longest_medial_vertex = vi;
    }

    if (cand.accepted) {
      ++stats.accepted;
      travel_sum += cand.travel;
      residual_sum += cand.residual;
    } else {
      ++stats.rejected;
      switch (cand.failure) {
      case albers::darboux_ridge_failure::nonfinite:
        ++reject_nonfinite;
        break;
      case albers::darboux_ridge_failure::max_travel:
        ++reject_max_travel;
        break;
      case albers::darboux_ridge_failure::not_converged:
        ++reject_not_converged;
        break;
      default:
        ++reject_other;
        break;
      }
    }
    out.push_back(cand);
  }

  if (stats.accepted > 0) {
    stats.avg_travel = travel_sum / real(stats.accepted);
    stats.avg_residual = residual_sum / real(stats.accepted);
  }
  stats.total = static_cast<int>(out.size());
  if (stats_out != nullptr) {
    *stats_out = stats;
  }
  std::cerr << "cyclide medial longest vertex=" << longest_medial_vertex
            << " zero_to_medial=" << longest_medial_length << std::endl;
  int probe_vertex = longest_medial_vertex;
  if (params.probe_vertex >= 0 &&
      params.probe_vertex < static_cast<int>(x.size())) {
    probe_vertex = params.probe_vertex;
    std::cerr << "cyclide medial cached probe vertex=" << probe_vertex
              << std::endl;
  }
  if (probe_vertex >= 0) {
    const vec3 local_start = x[static_cast<size_t>(probe_vertex)] - fit_origin;
    const vec3 zero_dir =
        -vertex_normals[static_cast<size_t>(probe_vertex)].normalized();
    dump_cyclide_medial_probe(Q, probe_vertex, local_start, zero_dir,
                              params.max_iters, params.tol, max_travel,
                              max_newton_step);
  }
  std::cout << "cyclide zero-surface diagnostics:"
            << " first_step_samples=" << first_step_samples
            << " first_step_inward_avg="
            << (first_step_samples > 0
                    ? first_step_inward_angle_sum / real(first_step_samples)
                    : 0.0)
            << " first_step_inward_max=" << first_step_inward_angle_max
            << " first_step_inward_gt_60=" << first_step_inward_gt_60
            << " first_step_root_avg="
            << (first_step_samples > 0
                    ? first_step_root_angle_sum / real(first_step_samples)
                    : 0.0)
            << " first_step_root_max=" << first_step_root_angle_max
            << " first_step_root_gt_60=" << first_step_root_gt_60
            << " reject_nonfinite=" << reject_nonfinite
            << " reject_max_travel=" << reject_max_travel
            << " reject_not_converged=" << reject_not_converged
            << " reject_other=" << reject_other
            << " medial_negative_step_candidates="
            << medial_negative_step_candidates
            << " medial_negative_steps=" << medial_negative_steps
            << " medial_backward_clamp_candidates="
            << medial_backward_clamp_candidates
            << " medial_backward_clamps=" << medial_backward_clamps
            << std::endl;
  return out;
}

inline std::vector<cyclide_medial_candidate>
compute_local_fit_cyclide_medial_candidates(
    asawa::shell::shell &M, const cyclide_medial_params &params,
    cyclide_medial_stats *stats_out = nullptr) {
  const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);
  const std::vector<vec3> vertex_normals =
      asawa::shell::vertex_normals(M, x);
  const real avg_len = asawa::shell::avg_length(M, x);
  const real l0 = std::max(params.l0_scale * avg_len, real(1e-12));
  auto verts = M.get_vert_range();
  if (verts.empty()) {
    if (stats_out != nullptr) {
      *stats_out = cyclide_medial_stats{};
    }
    return {};
  }

  std::vector<albers::vec14> fits;
  if (params.use_mls_bootstrap) {
    calder::mls_jet_bootstrap_params bp = params.bootstrap;
    bp.fit_p = params.fit_p;
    bp.fit_w0 = params.fit_w0;
    bp.radius_scale = params.radius_scale;
    fits = calder::darboux_fit_bootstrapped(M, x, vertex_normals, l0, bp);
  } else {
    fits = calder::darboux_cyclide_shell_fit(M, x, vertex_normals, l0,
                                             params.fit_p, params.fit_w0,
                                             params.radius_scale,
                                             params.stage1_radius);
  }

  const real max_travel = params.max_travel_scale * avg_len;
  const real max_newton_step =
      std::max(params.newton_step_scale * l0, real(1e-12));
  std::vector<cyclide_medial_candidate> out;
  cyclide_medial_stats stats;

  real travel_sum = 0.0;
  real residual_sum = 0.0;
  real grad_alignment_sum = 0.0;
  real min_grad_alignment = 1.0;
  int grad_alignment_samples = 0;
  int projection_attempted = 0;
  int projection_failed = 0;
  int projection_fail_max_travel = 0;
  int projection_fail_not_converged = 0;
  int projection_fail_other = 0;
  int reject_small_denominator = 0;
  int reject_nonfinite = 0;
  int reject_max_travel = 0;
  int reject_not_converged = 0;
  int reject_not_converged_clamped = 0;
  int reject_other = 0;
  real local_not_converged_residual_sum = 0.0;
  real local_not_converged_residual_max = 0.0;
  real local_not_converged_travel_sum = 0.0;
  real local_not_converged_travel_max = 0.0;
  int local_not_converged_samples = 0;
  std::vector<vec3> cylinder_pov(x.size(), vec3::Zero());
  std::vector<vec3> cylinder_normals(x.size(), vec3::UnitZ());
  std::vector<int> cylinder_valid(x.size(), 0);

  for (auto v_id : verts) {
    const int vi = static_cast<int>(v_id);
    if (vi < 0 || vi >= static_cast<int>(fits.size())) {
      continue;
    }

    cyclide_medial_candidate cand;
    cand.vertex = vi;

    const vec3 fit_origin = x[vi];
    const vec3 local_start = vec3::Zero();
    const albers::vec14 &Q = fits[vi];
    const vec3 inward = -vertex_normals[vi].normalized();

    const vec3 g_start = albers::darboux_grad(Q, local_start);
    real normal_alignment = -1.0;
    if (g_start.allFinite() && g_start.norm() > 1e-12) {
      vec3 fit_normal = g_start.normalized();
      const vec3 mesh_normal = vertex_normals[vi].normalized();
      if (fit_normal.dot(mesh_normal) < 0.0) {
        fit_normal *= -1.0;
      }
      const real align = fit_normal.dot(mesh_normal);
      normal_alignment = align;
      grad_alignment_sum += align;
      min_grad_alignment = std::min(min_grad_alignment, align);
      ++grad_alignment_samples;
    }

    std::vector<vec3> trace_local;
    albers::darboux_ridge_estimate est;
    vec3 medial_dir = vec3::Zero();
    if (g_start.allFinite() && g_start.norm() > 1e-12) {
      medial_dir = -g_start.normalized();
      est = albers::estimate_center_ridge(
          Q, local_start, medial_dir, params.max_iters, params.tol,
          &trace_local, max_travel, max_newton_step);
    } else {
      est.center = local_start;
      est.failure = albers::darboux_ridge_failure::invalid_direction;
    }

    cand.center_world = fit_origin + est.center;
    cand.center_local = est.center;
    cand.surface_world = fit_origin;
    cand.surface_local = vec3::Zero();
    cand.trace_world.reserve(trace_local.size());
    for (const vec3 &p_local : trace_local) {
      cand.trace_world.push_back(fit_origin + p_local);
    }
    cand.travel = est.travel;
    cand.residual = est.residual;
    cand.surface_residual = est.surface_residual;
    cand.iterations = est.iterations;
    cand.clamped_steps = est.clamped_steps;
    cand.failure = est.failure;
    cand.projection_failure = est.projection_failure;
    cand.projection_attempted = false;
    cand.projection_converged = true;
    cand.accepted = cyclide_medial_candidate_valid(est);
    if (normal_alignment < params.min_normal_alignment) {
      cand.accepted = false;
    }
    if (cand.accepted) {
      cylinder_pov[static_cast<size_t>(vi)] = cand.center_world;
      cylinder_normals[static_cast<size_t>(vi)] = inward;
      cylinder_valid[static_cast<size_t>(vi)] = 1;
    }
    if (params.probe_vertex == vi) {
      std::cerr << "local cyclide medial probe vertex=" << vi
                << " accepted=" << cand.accepted
                << " failure=" << static_cast<int>(cand.failure)
                << " travel=" << cand.travel
                << " residual=" << cand.residual
                << " normal_alignment=" << normal_alignment
                << " D_start=" << albers::eval_darboux(Q, local_start)
                << " |g_start|=" << g_start.norm()
                << " medial_dir=" << medial_dir.transpose() << std::endl;
      if (medial_dir.squaredNorm() > 1e-12) {
        const real ridge0 =
            albers::darboux_directional_deriv(Q, local_start, medial_dir);
        const real d2_0 =
            albers::darboux_directional_second_deriv(Q, local_start,
                                                     medial_dir);
        std::cerr << "  initial ridge=" << ridge0 << " d2=" << d2_0
                  << " raw_dt=" << (-ridge0 / d2_0) << std::endl;
        for (int ti = 0; ti < static_cast<int>(trace_local.size()); ++ti) {
          const vec3 &pt = trace_local[static_cast<size_t>(ti)];
          const real ridge =
              albers::darboux_directional_deriv(Q, pt, medial_dir);
          const real d2 =
              albers::darboux_directional_second_deriv(Q, pt, medial_dir);
          const real step =
              ti > 0 ? (pt - trace_local[static_cast<size_t>(ti - 1)]).norm()
                     : 0.0;
          std::cerr << "  trace[" << ti << "] x=" << pt.transpose()
                    << " D=" << albers::eval_darboux(Q, pt)
                    << " ridge=" << ridge << " d2=" << d2
                    << " step=" << step << std::endl;
        }
      }
    }

    if (cand.projection_attempted) {
      ++projection_attempted;
      if (!cand.projection_converged) {
        ++projection_failed;
        switch (cand.projection_failure) {
        case albers::darboux_ridge_failure::max_travel:
          ++projection_fail_max_travel;
          break;
        case albers::darboux_ridge_failure::not_converged:
          ++projection_fail_not_converged;
          break;
        default:
          ++projection_fail_other;
          break;
        }
      }
    }
    if (cand.accepted) {
      ++stats.accepted;
      travel_sum += cand.travel;
      residual_sum += cand.residual;
    } else {
      ++stats.rejected;
      switch (cand.failure) {
      case albers::darboux_ridge_failure::small_denominator:
        ++reject_small_denominator;
        break;
      case albers::darboux_ridge_failure::nonfinite:
        ++reject_nonfinite;
        break;
      case albers::darboux_ridge_failure::max_travel:
        ++reject_max_travel;
        break;
      case albers::darboux_ridge_failure::not_converged:
        ++reject_not_converged;
        ++local_not_converged_samples;
        local_not_converged_residual_sum += cand.residual;
        local_not_converged_residual_max =
            std::max(local_not_converged_residual_max, cand.residual);
        local_not_converged_travel_sum += cand.travel;
        local_not_converged_travel_max =
            std::max(local_not_converged_travel_max, cand.travel);
        if (cand.clamped_steps > 0) {
          ++reject_not_converged_clamped;
        }
        break;
      default:
        ++reject_other;
        break;
      }
    }
    out.push_back(cand);
  }

  const std::vector<vec6> cylinder_fits =
      calder::normal_aligned_line_convexity(M, cylinder_pov, cylinder_normals,
                                            l0, params.fit_p);
  for (cyclide_medial_candidate &cand : out) {
    if (!cand.accepted || cand.vertex < 0 ||
        cand.vertex >= static_cast<int>(cylinder_fits.size()) ||
        cylinder_valid[static_cast<size_t>(cand.vertex)] == 0) {
      continue;
    }
    vec3 axis = albers::plucker_line_direction(
        cylinder_fits[static_cast<size_t>(cand.vertex)]);
    if (!axis.allFinite() || axis.norm() < 1e-12) {
      continue;
    }
    if (axis.dot(cand.center_local) < 0.0) {
      axis *= -1.0;
    }
    cand.cylinder_axis_world = axis.normalized();
  }

  smooth_pipe_fit_medial_axis(M, out, params.medial_smooth_scale,
                              params.medial_smooth_blend);

  travel_sum = 0.0;
  residual_sum = 0.0;
  for (const cyclide_medial_candidate &cand : out) {
    if (!cand.accepted) {
      continue;
    }
    travel_sum += cand.travel;
    residual_sum += cand.residual;
  }

  if (stats.accepted > 0) {
    stats.avg_travel = travel_sum / real(stats.accepted);
    stats.avg_residual = residual_sum / real(stats.accepted);
  }
  stats.total = static_cast<int>(out.size());
  if (stats_out != nullptr) {
    *stats_out = stats;
  }
  std::cout << "cyclide local-fit medial diagnostics:"
            << " grad_align_avg="
            << (grad_alignment_samples > 0
                    ? grad_alignment_sum / real(grad_alignment_samples)
                    : 0.0)
            << " grad_align_min=" << min_grad_alignment
            << " projection_attempted=" << projection_attempted
            << " projection_failed=" << projection_failed
            << " projection_fail_max_travel=" << projection_fail_max_travel
            << " projection_fail_not_converged="
            << projection_fail_not_converged
            << " projection_fail_other=" << projection_fail_other
            << " reject_small_denominator=" << reject_small_denominator
            << " reject_nonfinite=" << reject_nonfinite
            << " reject_max_travel=" << reject_max_travel
            << " reject_not_converged=" << reject_not_converged
            << " reject_not_converged_clamped="
            << reject_not_converged_clamped
            << " not_converged_residual_avg="
            << (local_not_converged_samples > 0
                    ? local_not_converged_residual_sum /
                          real(local_not_converged_samples)
                    : 0.0)
            << " not_converged_residual_max="
            << local_not_converged_residual_max
            << " not_converged_travel_avg="
            << (local_not_converged_samples > 0
                    ? local_not_converged_travel_sum /
                          real(local_not_converged_samples)
                    : 0.0)
            << " not_converged_travel_max=" << local_not_converged_travel_max
            << " l0=" << l0 << " max_step=" << max_newton_step
            << " max_travel=" << max_travel
            << " reject_other=" << reject_other << std::endl;
  return out;
}

inline void print_cyclide_medial_stats(std::ostream &os,
                                       const cyclide_medial_stats &stats) {
  os << "cyclide medial: total=" << stats.total
     << " accepted=" << stats.accepted << " rejected=" << stats.rejected
     << " avg_travel=" << stats.avg_travel
     << " avg_residual=" << stats.avg_residual << std::endl;
}

inline cyclide_medial_params default_cyclide_medial_demo_params() {
  cyclide_medial_params params;
  params.l0_scale = 2.0;
  params.fit_p = 3.0;
  params.normal_l0 = 2.0;
  params.min_normal_alignment = -0.5;
  // Off while investigating raw ridge outliers (was MLS pipe smooth).
  params.medial_smooth_scale = 0.0;
  params.medial_smooth_blend = 0.0;
  params.max_iters = 100;
  params.tol = 1e-8;
  params.max_travel_scale = 1000.0;
  params.newton_step_scale = 1000.0;
  // Stage-1 harmonic → R from |κ|_max → soft-convex × (Nj·∇T̂) × G(σ=R).
  params.use_mls_bootstrap = true;
  params.fit_p = 3.0;
  params.fit_w0 = 1.0;
  params.radius_scale = 1.0;
  params.bootstrap.ablation =
      calder::mls_bootstrap_ablation::harmonic_darboux_torus_normal_align;
  params.bootstrap.torus_sign = albers::torus_sign_mode::pat_eq22;
  return params;
}

inline void normalize_cyclide_medial_demo_mesh(asawa::shell::shell &M,
                                               real target_radius = 1.5) {
  std::vector<vec3> &x = asawa::get_vec_data(M, 0);
  if (x.empty()) {
    return;
  }

  vec3 lo = x.front();
  vec3 hi = x.front();
  for (const vec3 &p : x) {
    lo = lo.cwiseMin(p);
    hi = hi.cwiseMax(p);
  }

  const vec3 center = 0.5 * (lo + hi);
  real radius = 0.0;
  for (const vec3 &p : x) {
    radius = std::max(radius, (p - center).norm());
  }
  if (radius < 1e-12) {
    return;
  }

  const real scale = target_radius / radius;
  for (vec3 &p : x) {
    p = scale * (p - center);
  }
  std::cerr << "normalized mesh radius=" << target_radius
            << " scale=" << scale << std::endl;
}

class darboux_cyclide_medial_demo {
public:
  using ptr = std::shared_ptr<darboux_cyclide_medial_demo>;

  static ptr create(cyclide_medial_params params =
                        default_cyclide_medial_demo_params()) {
    return std::make_shared<darboux_cyclide_medial_demo>(params);
  }

  explicit darboux_cyclide_medial_demo(cyclide_medial_params params)
      : __params(params) {
    __M = asawa::shell::load_bunny();
    asawa::shell::triangulate(*__M);
    normalize_cyclide_medial_demo_mesh(*__M);
    configure_scene_frame();
    __candidates =
        compute_local_fit_cyclide_medial_candidates(*__M, __params, &__stats);
    std::cerr << "bunny local-fit zero-to-medial" << std::endl;
    print_cyclide_medial_stats(std::cerr, __stats);
  }

  void step(int frame) {
    _frame = frame;
    draw_medial_rays();
  }

  int frame() const { return _frame; }
  const cyclide_medial_stats &stats() const { return __stats; }
  const std::vector<cyclide_medial_candidate> &candidates() const {
    return __candidates;
  }

  asawa::shell::shell::ptr __M;

private:
  void configure_scene_frame() {
    const std::vector<vec3> &x = asawa::const_get_vec_data(*__M, 0);
    if (x.empty()) {
      return;
    }

    vec3 lo = x.front();
    vec3 hi = x.front();
    for (const vec3 &p : x) {
      lo = lo.cwiseMin(p);
      hi = hi.cwiseMax(p);
    }

    __center = 0.5 * (lo + hi);
    __major_radius = 0.5 * (hi - lo).norm();
    __minor_radius = 4.0 * asawa::shell::avg_length(*__M, x);
    std::cerr << "loaded bunny verts=" << __M->vert_count()
              << " faces=" << __M->face_count()
              << " scene_radius=" << __major_radius
              << " slice_extent=" << __minor_radius << std::endl;
  }

  void draw_medial_rays() {
    const vec4 zero_color(0.0, 0.85, 1.0, 1.0);
    const vec4 medial_color(1.0, 0.65, 0.05, 1.0);
    const vec4 cylinder_axis_color(0.9, 0.15, 1.0, 1.0);
    const vec4 axis_color(0.35, 0.35, 0.35, 1.0);

    geometry_logger::line(__center - 1.7 * __major_radius * vec3::UnitZ(),
                          __center + 1.7 * __major_radius * vec3::UnitZ(),
                          axis_color);

    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    const real cylinder_axis_len = 8.0 * asawa::shell::avg_length(*__M, x);

    for (const auto &cand : __candidates) {
      if (cand.vertex < 0 || cand.vertex >= static_cast<int>(x.size())) {
        continue;
      }

      const vec3 p = x[static_cast<size_t>(cand.vertex)];
      if (cand.projection_converged) {
        if ((cand.surface_world - p).norm() > 1e-12) {
          geometry_logger::line(p, cand.surface_world, zero_color);
        } else {
          geometry_logger::point(cand.surface_world, zero_color);
        }
      }
      if (cand.projection_converged && cand.center_world.allFinite() &&
          (cand.center_world - cand.surface_world).norm() > 1e-12) {
        geometry_logger::line(cand.surface_world, cand.center_world,
                              medial_color);
      }
      if (cand.accepted && cand.cylinder_axis_world.squaredNorm() > 1e-12) {
        const vec3 dp =
            cylinder_axis_len * cand.cylinder_axis_world.normalized();
        geometry_logger::line(p - dp, p + dp, cylinder_axis_color);
      }
    }
  }

  cyclide_medial_params __params;
  cyclide_medial_stats __stats;
  std::vector<cyclide_medial_candidate> __candidates;
  vec3 __center = vec3::Zero();
  real __major_radius = 1.25;
  real __minor_radius = 0.35;
  int _frame = 0;
};

} // namespace duchamp
} // namespace gaudi

#endif // __GAUDI_DARBOUX_CYCLIDE_MEDIAL_HPP__
