#ifndef __GAUDI_DARBOUX_CYCLIDE_TESTS_HPP__
#define __GAUDI_DARBOUX_CYCLIDE_TESTS_HPP__

#include "gaudi/albers/darboux_cyclide.hpp"
#include "gaudi/asawa/datums.hpp"
#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/test/darboux_cyclide_torus_fixture.hpp"
#include "gaudi/calder/least_squares_fit.hpp"
#include "gaudi/duchamp/darboux_cyclide_medial.hpp"
#include "gaudi/test/test.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace gaudi {
namespace test {

GAUDI_TEST(calder_darboux_cyclide_fits_offset_torus) {
  const real major_radius = 1.25;
  const real minor_radius = 0.35;
  const TorusFrame frame =
      make_torus_frame(vec3(0.7, -0.45, 0.3), vec3(0.35, 0.65, 0.9));

  TorusMesh torus =
      make_offset_torus_shell(18, 12, major_radius, minor_radius, frame);
  GAUDI_ASSERT(torus.shell != nullptr);
  asawa::shell::shell &M = *torus.shell;
  const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);

  std::vector<vec3> face_centers = asawa::shell::face_centers(M, x);
  std::vector<real> face_areas = asawa::shell::face_areas(M, x);
  std::vector<vec3> face_normals = asawa::shell::face_normals(M, x);
  const auto face_range = M.get_face_range();
  GAUDI_ASSERT(!face_range.empty());

  std::vector<vec3> weighted_normals(M.face_count(), vec3::Zero());
  for (auto f : face_range) {
    weighted_normals[f] = face_areas[f] * face_normals[f];
  }

  const real l0 = 2.0 * asawa::shell::avg_length(M, x);
  using Cyclide = albers::darboux_cyclide;
  std::vector<Cyclide::coefficients> fits =
      calder::generic_fit<albers::darboux_cyclide, calder::shell_bundle>(
          M, weighted_normals, face_centers, face_normals, l0, 3.0,
          calder::shell_inv_dist_weight);
  GAUDI_ASSERT(fits.size() == face_centers.size());

  real max_surface_residual = 0.0;
  real max_center_residual = 0.0;
  real sum_surface_residual2 = 0.0;
  real sum_normal_alignment = 0.0;
  int checked_samples = 0;
  int checked_normals = 0;

  for (auto fi : face_range) {
    const Cyclide::coefficients &Q = fits[fi];
    GAUDI_EXPECT(Q.allFinite());

    const vec3 pi = face_centers[fi];
    max_center_residual =
        std::max(max_center_residual, std::abs(albers::eval_darboux(Q, vec3::Zero())));
    for (auto fj : face_range) {
      const vec3 dp = face_centers[fj] - pi;
      if (dp.norm() > l0) {
        continue;
      }

      const real residual = std::abs(albers::eval_darboux(Q, dp));
      max_surface_residual = std::max(max_surface_residual, residual);
      sum_surface_residual2 += residual * residual;

      const vec3 grad = albers::darboux_grad(Q, dp);
      if (grad.norm() > 1e-10) {
        const real alignment =
            std::abs(grad.normalized().dot(face_normals[fj]));
        sum_normal_alignment += alignment;
        ++checked_normals;
      }
      ++checked_samples;
    }
  }

  const real rms_surface_residual =
      std::sqrt(sum_surface_residual2 / real(std::max(checked_samples, 1)));
  const real mean_normal_alignment =
      sum_normal_alignment / real(std::max(checked_normals, 1));

  std::cerr << "\n[calder_darboux_cyclide_fits_offset_torus]"
            << " faces=" << face_range.size() << " l0=" << l0
            << " checked_samples=" << checked_samples
            << " checked_normals=" << checked_normals
            << " max_center_residual=" << max_center_residual
            << " max_surface_residual=" << max_surface_residual
            << " rms_surface_residual=" << rms_surface_residual
            << " mean_normal_alignment=" << mean_normal_alignment << "\n";

  GAUDI_ASSERT(checked_samples > 0);
  GAUDI_ASSERT(checked_normals > 0);
  GAUDI_EXPECT(max_center_residual < 7.5e-2);
  GAUDI_EXPECT(max_surface_residual < 1e-1);
  GAUDI_EXPECT(rms_surface_residual < 5e-2);
  GAUDI_EXPECT(mean_normal_alignment > 0.95);
}

GAUDI_TEST(calder_darboux_cyclide_fits_canonical_torus_vertices) {
  const real major_radius = 1.25;
  const real minor_radius = 0.35;
  const TorusFrame frame = make_torus_frame(vec3::Zero(), vec3::UnitZ());

  TorusMesh torus =
      make_offset_torus_shell(36, 20, major_radius, minor_radius, frame);
  GAUDI_ASSERT(torus.shell != nullptr);
  asawa::shell::shell &M = *torus.shell;
  const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);
  const std::vector<vec3> vertex_normals = asawa::shell::vertex_normals(M, x);
  const real l0 = 2.0 * asawa::shell::avg_length(M, x);

  using Cyclide = albers::darboux_cyclide;
  const std::vector<Cyclide::coefficients> fits =
      calder::darboux_cyclide(M, x, vertex_normals, l0, 3.0);
  GAUDI_ASSERT(fits.size() == x.size());

  real sum_alignment = 0.0;
  int checked = 0;
  auto vert_range = M.get_vert_range();
  for (auto vi : vert_range) {
    const vec3 mesh_n = vertex_normals[vi].normalized();
    vec3 cyclide_n = -fits[vi][9] * albers::darboux_grad(fits[vi], vec3::Zero());
    if (cyclide_n.norm() < 1e-10) {
      continue;
    }
    cyclide_n.normalize();
    if (cyclide_n.dot(mesh_n) < 0.0) {
      cyclide_n *= -1.0;
    }
    sum_alignment += cyclide_n.dot(mesh_n);
    ++checked;
  }

  const real mean_alignment = sum_alignment / real(std::max(checked, 1));
  std::cerr << "\n[calder_darboux_cyclide_fits_canonical_torus_vertices]"
            << " verts=" << vert_range.size() << " checked=" << checked
            << " mean_normal_alignment=" << mean_alignment << "\n";

  GAUDI_ASSERT(checked > 0);
  GAUDI_EXPECT(mean_alignment > 0.90);
}

GAUDI_TEST(duchamp_cyclide_medial_newton_offset_torus) {
  const TorusFrame frame =
      make_torus_frame(vec3(0.7, -0.45, 0.3), vec3(0.25, 0.65, 0.9));
  TorusMesh torus = make_offset_torus_shell(36, 20, 1.25, 0.35, frame);
  asawa::shell::shell &M = *torus.shell;
  const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);
  const std::vector<vec3> vertex_normals = asawa::shell::vertex_normals(M, x);
  const real avg_len = asawa::shell::avg_length(M, x);
  duchamp::cyclide_medial_params params;
  duchamp::cyclide_medial_stats stats;
  albers::vec14 Q;
  auto cands = duchamp::compute_single_fit_cyclide_medial_candidates(
      M, 0, params, &stats, &Q);

  int outside = 0;
  int outside_accepted = 0;
  int outside_rejected = 0;
  int crazy_travel = 0;
  const vec3 fit_origin = x[0];
  const real max_travel = params.max_travel_scale * avg_len;
  for (const auto &cand : cands) {
    const vec3 local_start = x[cand.vertex] - fit_origin;
    if (albers::eval_darboux(Q, local_start) > params.tol) {
      ++outside;
      if (cand.accepted) {
        ++outside_accepted;
      } else {
        ++outside_rejected;
      }
    }
    if (cand.accepted && cand.travel > max_travel) {
      ++crazy_travel;
    }
  }

  std::cerr << "\n[duchamp_cyclide_medial_newton_offset_torus]"
            << " accepted=" << stats.accepted << " rejected=" << stats.rejected
            << " avg_travel=" << stats.avg_travel
            << " avg_residual=" << stats.avg_residual
            << " outside=" << outside << " outside_accepted=" << outside_accepted
            << " outside_rejected=" << outside_rejected
            << " crazy_travel=" << crazy_travel << "\n";

  GAUDI_EXPECT(stats.accepted > 0);
  GAUDI_ASSERT(!cands.empty());
  GAUDI_EXPECT(crazy_travel == 0);
}

} // namespace test
} // namespace gaudi

#endif // __GAUDI_DARBOUX_CYCLIDE_TESTS_HPP__
