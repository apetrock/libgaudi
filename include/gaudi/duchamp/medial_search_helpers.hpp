#ifndef __GAUDI_DUCHAMP_MEDIAL_SEARCH_HELPERS__
#define __GAUDI_DUCHAMP_MEDIAL_SEARCH_HELPERS__

#include "gaudi/albers/darboux_cyclide.hpp"
#include "gaudi/albers/darboux_medial_geometry.hpp"
#include "gaudi/albers/medial_shape_search.hpp"
#include "gaudi/duchamp/medial_result_types.hpp"

namespace gaudi {
namespace duchamp {

inline bool legacy_medial_accepted(const albers::darboux_ridge_estimate &est) {
  return est.accepted && est.center.allFinite() && est.travel > 1e-12 &&
         std::isfinite(est.residual);
}

inline medial_point search_medial_shape_energy(
    const albers::vec14 &Q, const vec3 &foot, const vec3 &mesh_normal,
    real max_travel, const albers::medial_shape_search_params &search) {
  medial_point result;
  const vec3 dir = albers::aligned_inward_ray_dir(Q, mesh_normal);
  if (dir.squaredNorm() < 1e-24) {
    return result;
  }
  const auto found =
      albers::search_medial_along_ray(Q, foot, dir, max_travel, search);
  result.point = found.center_world;
  result.accepted = found.converged && found.center_world.allFinite() &&
                    std::isfinite(found.energy) && found.travel > 1e-12;
  return result;
}

inline medial_point search_medial_legacy_ridge(
    const albers::vec14 &Q, const vec3 &foot, const vec3 &mesh_normal,
    real max_travel, int max_iters, real tol, real max_newton_step,
    real min_normal_alignment) {
  medial_point result;
  const vec3 local_start = vec3::Zero();
  const vec3 g_start = albers::darboux_grad(Q, local_start);
  real normal_alignment = -1.0;
  albers::darboux_ridge_estimate est;
  if (g_start.allFinite() && g_start.norm() > 1e-12) {
    vec3 fit_normal = g_start.normalized();
    const vec3 n = mesh_normal.normalized();
    if (fit_normal.dot(n) < 0.0) {
      fit_normal *= -1.0;
    }
    normal_alignment = fit_normal.dot(n);
    const vec3 medial_dir = -g_start.normalized();
    est = albers::estimate_center_ridge(Q, local_start, medial_dir, max_iters,
                                        tol, nullptr, max_travel,
                                        max_newton_step);
  } else {
    est.center = local_start;
    est.failure = albers::darboux_ridge_failure::invalid_direction;
  }

  result.point = foot + est.center;
  result.accepted = legacy_medial_accepted(est);
  if (normal_alignment < min_normal_alignment) {
    result.accepted = false;
  }
  return result;
}

} // namespace duchamp
} // namespace gaudi

#endif // __GAUDI_DUCHAMP_MEDIAL_SEARCH_HELPERS__
