#ifndef __GAUDI_DUCHAMP_DARBOUX_CYCLIDE_MEDIAL_HPP__
#define __GAUDI_DUCHAMP_DARBOUX_CYCLIDE_MEDIAL_HPP__

#include "gaudi/albers/darboux_cyclide.hpp"
#include "gaudi/asawa/datums.hpp"
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/calder/least_squares_fit.hpp"
#include "gaudi/common.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

namespace gaudi {
namespace duchamp {

struct cyclide_medial_params {
  real l0_scale = 0.1;
  real fit_p = 3.0;
  int max_iters = 12;
  real tol = 1e-8;
  real max_travel_scale = 12.0;
};

struct cyclide_medial_candidate {
  int vertex = -1;
  vec3 center_world = vec3::Zero();
  vec3 center_local = vec3::Zero();
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
      calder::darboux_cyclide(M, p_fit, n_fit, l0, params.fit_p);
  const albers::vec14 Q = fits.empty() ? albers::vec14::Zero() : fits[0];
  if (fit_out != nullptr) {
    *fit_out = Q;
  }

  const vec3 fit_origin = x[fit_vertex];
  const real max_travel = params.max_travel_scale * avg_len;
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

} // namespace duchamp
} // namespace gaudi

#endif // __GAUDI_DARBOUX_CYCLIDE_MEDIAL_HPP__
