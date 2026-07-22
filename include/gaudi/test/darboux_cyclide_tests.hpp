#ifndef __GAUDI_DARBOUX_CYCLIDE_TESTS_HPP__
#define __GAUDI_DARBOUX_CYCLIDE_TESTS_HPP__

#include "gaudi/albers/darboux_cyclide.hpp"
#include "gaudi/asawa/datums.hpp"
#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/arp/tree_view.hpp"
#include "gaudi/test/darboux_cyclide_torus_fixture.hpp"
#include "gaudi/calder/integrators.hpp"
#include "gaudi/calder/least_squares_fit.hpp"
#include "gaudi/duchamp/darboux_cyclide_medial.hpp"
#include "gaudi/test/test.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace gaudi {
namespace test {

inline vec3 canonical_torus_normal(const vec3 &p, real major_radius) {
  vec3 radial(p[0], p[1], 0.0);
  if (radial.norm() < 1e-12) {
    return vec3::UnitZ();
  }
  radial.normalize();
  const vec3 tube_center = major_radius * radial;
  return (p - tube_center).normalized();
}

GAUDI_TEST(albers_darboux_cyclide_direct_fit_canonical_torus) {
  const real major_radius = 1.25;
  const real minor_radius = 0.35;
  const vec3 pov = torus_point(0.37, 1.13, major_radius, minor_radius,
                               make_torus_frame(vec3::Zero(), vec3::UnitZ()));

  albers::darboux_cyclide fit;
  int samples = 0;
  for (int i = 0; i < 48; ++i) {
    const real u = 2.0 * M_PI * real(i) / real(48);
    for (int j = 0; j < 24; ++j) {
      const real v = 2.0 * M_PI * real(j) / real(24);
      const vec3 p =
          torus_point(u, v, major_radius, minor_radius,
                      make_torus_frame(vec3::Zero(), vec3::UnitZ()));
      const vec3 n = canonical_torus_normal(p, major_radius);
      fit.accumulate(1.0, p - pov, n);
      ++samples;
    }
  }

  const albers::vec14 Q = fit.solve();
  GAUDI_EXPECT(Q.allFinite());

  real max_surface_residual = 0.0;
  real rms_surface_residual = 0.0;
  real mean_normal_alignment = 0.0;
  int checked = 0;
  for (int i = 0; i < 48; ++i) {
    const real u = 2.0 * M_PI * real(i) / real(48);
    for (int j = 0; j < 24; ++j) {
      const real v = 2.0 * M_PI * real(j) / real(24);
      const vec3 p =
          torus_point(u, v, major_radius, minor_radius,
                      make_torus_frame(vec3::Zero(), vec3::UnitZ()));
      const vec3 local = p - pov;
      const real residual = std::abs(albers::eval_darboux(Q, local));
      max_surface_residual = std::max(max_surface_residual, residual);
      rms_surface_residual += residual * residual;

      const vec3 g = albers::darboux_grad(Q, local);
      if (g.norm() > 1e-12) {
        mean_normal_alignment +=
            std::abs(g.normalized().dot(canonical_torus_normal(p, major_radius)));
        ++checked;
      }
    }
  }

  rms_surface_residual = std::sqrt(rms_surface_residual / real(samples));
  mean_normal_alignment /= real(std::max(checked, 1));
  std::cerr << "\n[albers_darboux_cyclide_direct_fit_canonical_torus]"
            << " samples=" << samples
            << " max_surface_residual=" << max_surface_residual
            << " rms_surface_residual=" << rms_surface_residual
            << " mean_normal_alignment=" << mean_normal_alignment << "\n";

  GAUDI_EXPECT(max_surface_residual < 1e-6);
  GAUDI_EXPECT(rms_surface_residual < 1e-7);
  GAUDI_EXPECT(mean_normal_alignment > 0.999);
}

GAUDI_TEST(calder_darboux_cyclide_identity_fits_canonical_torus) {
  const real major_radius = 1.25;
  const real minor_radius = 0.35;
  const TorusFrame frame = make_torus_frame(vec3::Zero(), vec3::UnitZ());

  TorusMesh torus =
      make_offset_torus_shell(48, 24, major_radius, minor_radius, frame);
  GAUDI_ASSERT(torus.shell != nullptr);
  asawa::shell::shell &M = *torus.shell;
  const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);

  std::vector<vec3> face_centers = asawa::shell::face_centers(M, x);
  std::vector<real> face_areas = asawa::shell::face_areas(M, x);
  const auto face_range = M.get_face_range();
  GAUDI_ASSERT(!face_range.empty());

  std::vector<vec3> analytic_normals(M.face_count(), vec3::Zero());
  std::vector<vec3> weighted_normals(M.face_count(), vec3::Zero());
  for (auto f : face_range) {
    analytic_normals[f] = canonical_torus_normal(face_centers[f], major_radius);
    weighted_normals[f] = face_areas[f] * analytic_normals[f];
  }

  using Cyclide = albers::darboux_cyclide;
  const std::vector<Cyclide::coefficients> fits =
      calder::generic_fit<albers::darboux_cyclide, calder::shell_bundle>(
          M, weighted_normals, face_centers, analytic_normals,
          2.0 * asawa::shell::avg_length(M, x), 3.0,
          calder::shell_identity_weight);
  GAUDI_ASSERT(fits.size() == face_centers.size());

  real max_surface_residual = 0.0;
  real rms_surface_residual = 0.0;
  real mean_normal_alignment = 0.0;
  int checked = 0;
  int checked_normals = 0;
  for (auto fi : face_range) {
    const Cyclide::coefficients &Q = fits[fi];
    GAUDI_EXPECT(Q.allFinite());
    const vec3 pi = face_centers[fi];
    for (auto fj : face_range) {
      const vec3 local = face_centers[fj] - pi;
      const real residual = std::abs(albers::eval_darboux(Q, local));
      max_surface_residual = std::max(max_surface_residual, residual);
      rms_surface_residual += residual * residual;

      const vec3 g = albers::darboux_grad(Q, local);
      if (g.norm() > 1e-12) {
        mean_normal_alignment +=
            std::abs(g.normalized().dot(analytic_normals[fj]));
        ++checked_normals;
      }
      ++checked;
    }
  }

  rms_surface_residual = std::sqrt(rms_surface_residual / real(checked));
  mean_normal_alignment /= real(std::max(checked_normals, 1));
  std::cerr << "\n[calder_darboux_cyclide_identity_fits_canonical_torus]"
            << " fits=" << fits.size() << " checked=" << checked
            << " max_surface_residual=" << max_surface_residual
            << " rms_surface_residual=" << rms_surface_residual
            << " mean_normal_alignment=" << mean_normal_alignment << "\n";

  GAUDI_EXPECT(max_surface_residual < 1e-2);
  GAUDI_EXPECT(rms_surface_residual < 1e-3);
  GAUDI_EXPECT(mean_normal_alignment > 0.999);
}

GAUDI_TEST(calder_shell_tree_node_data_matches_leaf_sums) {
  const real major_radius = 1.25;
  const real minor_radius = 0.35;
  const TorusFrame frame = make_torus_frame(vec3::Zero(), vec3::UnitZ());

  TorusMesh torus =
      make_offset_torus_shell(18, 12, major_radius, minor_radius, frame);
  GAUDI_ASSERT(torus.shell != nullptr);
  asawa::shell::shell &M = *torus.shell;
  const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);

  std::vector<index_t> face_vert_ids = M.get_face_vert_ids();
  auto face_ids_typed = M.get_face_range();
  std::vector<index_t> face_ids(face_ids_typed.begin(), face_ids_typed.end());
  std::vector<real> areas = asawa::shell::face_areas(M, x);
  std::vector<vec3> centers = asawa::shell::face_centers(M, x);
  std::vector<vec3> normals = asawa::shell::face_normals(M, x);
  std::vector<vec3> weighted_normals(M.face_count(), vec3::Zero());
  for (auto f : face_ids_typed) {
    weighted_normals[f] = areas[f] * normals[f];
  }

  calder::Shell_Tree_Type::ptr tree =
      calder::Shell_Tree_Type::create(face_vert_ids, x, 16);
  const arp::tree_view view = arp::make_tree_view(*tree);

  auto area_datum = calder::scalar_datum::create(face_ids, areas);
  auto normal_datum = calder::vec3_datum::create(face_ids, weighted_normals);
  auto com = calder::com_datum::create(face_ids, areas, centers);
  area_datum->pyramid(*tree);
  normal_datum->pyramid(*tree);
  com->pyramid(*tree);

  real max_area_error = 0.0;
  real max_normal_error = 0.0;
  real max_com_error = 0.0;
  int checked_nodes = 0;
  for (index_t node_id = 0;
       node_id < static_cast<index_t>(view.internal_count()); ++node_id) {
    real exact_area = 0.0;
    vec3 exact_normal = vec3::Zero();
    vec3 exact_com_numer = vec3::Zero();

    for (index_t leaf_id = 0;
         leaf_id < static_cast<index_t>(view.leaf_count()); ++leaf_id) {
      index_t parent = view.leaf_parent(leaf_id);
      bool in_subtree = false;
      while (parent != arp::UNULL) {
        if (parent == node_id) {
          in_subtree = true;
          break;
        }
        parent = view.internal_parent(parent);
      }
      if (!in_subtree) {
        continue;
      }

      const index_t orig = view.index(leaf_id);
      const index_t face = face_ids[orig];
      exact_area += areas[face];
      exact_normal += weighted_normals[face];
      exact_com_numer += areas[face] * centers[face];
    }

    if (exact_area <= 0.0) {
      continue;
    }

    const real node_area = area_datum->node_data()[node_id];
    const vec3 node_normal = normal_datum->node_data()[node_id];
    const vec3 node_com = com->get_node_com(node_id);
    const vec3 exact_com = exact_com_numer / exact_area;

    max_area_error = std::max(max_area_error, std::abs(node_area - exact_area));
    max_normal_error =
        std::max(max_normal_error, (node_normal - exact_normal).norm());
    max_com_error = std::max(max_com_error, (node_com - exact_com).norm());
    ++checked_nodes;
  }

  std::cerr << "\n[calder_shell_tree_node_data_matches_leaf_sums]"
            << " checked_nodes=" << checked_nodes
            << " max_area_error=" << max_area_error
            << " max_normal_error=" << max_normal_error
            << " max_com_error=" << max_com_error << "\n";

  GAUDI_ASSERT(checked_nodes > 0);
  GAUDI_EXPECT(max_area_error < 1e-12);
  GAUDI_EXPECT(max_normal_error < 1e-12);
  GAUDI_EXPECT(max_com_error < 1e-12);
}

GAUDI_TEST(calder_shell_leaf_normal_data_matches_face_indices) {
  const real major_radius = 1.25;
  const real minor_radius = 0.35;
  const TorusFrame frame = make_torus_frame(vec3::Zero(), vec3::UnitZ());

  TorusMesh torus =
      make_offset_torus_shell(18, 12, major_radius, minor_radius, frame);
  GAUDI_ASSERT(torus.shell != nullptr);
  asawa::shell::shell &M = *torus.shell;
  const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);

  std::vector<index_t> face_vert_ids = M.get_face_vert_ids();
  auto face_ids_typed = M.get_face_range();
  std::vector<index_t> face_ids(face_ids_typed.begin(), face_ids_typed.end());

  std::vector<real> synthetic_area(M.face_count(), 0.0);
  std::vector<vec3> synthetic_normal(M.face_count(), vec3::Zero());
  for (auto f : face_ids_typed) {
    synthetic_area[f] = real(f) + 0.25;
    synthetic_normal[f] =
        vec3(real(f + 1), real(2 * f + 3), real(-3 * f - 5));
  }

  calder::Shell_Tree_Type::ptr tree =
      calder::Shell_Tree_Type::create(face_vert_ids, x, 16);
  const arp::tree_view view = arp::make_tree_view(*tree);

  calder::Shell_Sum_Type sum(*tree);
  sum.bind(calder::vec3_datum::create(face_ids, synthetic_normal));
  sum.bind(calder::scalar_datum::create(face_ids, synthetic_area));
  for (auto &datum : sum.__data) {
    datum->pyramid(*tree);
  }

  real max_normal_error = 0.0;
  real max_area_error = 0.0;
  int checked = 0;
  for (index_t leaf_id = 0;
       leaf_id < static_cast<index_t>(view.leaf_count()); ++leaf_id) {
    const index_t orig = view.get_index(leaf_id);
    const index_t face = face_ids[orig];

    const vec3 got_normal =
        calder::get_data<vec3>(calder::LEAF, leaf_id, 0, sum.__data);
    const real got_area =
        calder::get_data<real>(calder::LEAF, leaf_id, 1, sum.__data);

    max_normal_error =
        std::max(max_normal_error, (got_normal - synthetic_normal[face]).norm());
    max_area_error =
        std::max(max_area_error, std::abs(got_area - synthetic_area[face]));
    ++checked;
  }

  std::cerr << "\n[calder_shell_leaf_normal_data_matches_face_indices]"
            << " checked=" << checked
            << " max_normal_error=" << max_normal_error
            << " max_area_error=" << max_area_error << "\n";

  GAUDI_ASSERT(checked == static_cast<int>(face_ids.size()));
  GAUDI_EXPECT(max_normal_error < 1e-12);
  GAUDI_EXPECT(max_area_error < 1e-12);
}

GAUDI_TEST(calder_shell_bh_traversal_covers_each_leaf_once) {
  const real major_radius = 1.25;
  const real minor_radius = 0.35;
  const TorusFrame frame = make_torus_frame(vec3::Zero(), vec3::UnitZ());

  TorusMesh torus =
      make_offset_torus_shell(18, 12, major_radius, minor_radius, frame);
  GAUDI_ASSERT(torus.shell != nullptr);
  asawa::shell::shell &M = *torus.shell;
  const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);

  std::vector<index_t> face_vert_ids = M.get_face_vert_ids();
  std::vector<vec3> centers = asawa::shell::face_centers(M, x);
  calder::Shell_Tree_Type::ptr tree =
      calder::Shell_Tree_Type::create(face_vert_ids, x, 16);
  const arp::tree_view view = arp::make_tree_view(*tree);

  auto mark_subtree = [&](index_t node_id, std::vector<int> &counts) {
    for (index_t leaf_id = 0;
         leaf_id < static_cast<index_t>(view.leaf_count()); ++leaf_id) {
      index_t parent = view.leaf_parent(leaf_id);
      while (parent != arp::UNULL) {
        if (parent == node_id) {
          counts[leaf_id] += 1;
          break;
        }
        parent = view.internal_parent(parent);
      }
    }
  };

  real eps = 0.25;
  int max_misses = 0;
  int max_duplicates = 0;
  for (size_t qi = 0; qi < std::min<size_t>(centers.size(), 16); ++qi) {
    std::vector<int> counts(view.leaf_count(), 0);
    calder::traverse_bh_opening(
        *tree, static_cast<index_t>(qi), centers[qi], eps,
        [&](index_t node_id, index_t, const vec3 &) {
          mark_subtree(node_id, counts);
        },
        [&](index_t leaf_id, index_t, index_t, const vec3 &) {
          counts[leaf_id] += 1;
        });

    int misses = 0;
    int duplicates = 0;
    for (int c : counts) {
      if (c == 0) {
        ++misses;
      } else if (c > 1) {
        duplicates += c - 1;
      }
    }
    max_misses = std::max(max_misses, misses);
    max_duplicates = std::max(max_duplicates, duplicates);
  }

  std::cerr << "\n[calder_shell_bh_traversal_covers_each_leaf_once]"
            << " max_misses=" << max_misses
            << " max_duplicates=" << max_duplicates << "\n";

  GAUDI_EXPECT(max_misses == 0);
  GAUDI_EXPECT(max_duplicates == 0);
}

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
    vec3 cyclide_n = albers::darboux_grad(fits[vi], vec3::Zero());
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

GAUDI_TEST(calder_darboux_cyclide_convexity_fit_w0_foot_normal) {
  const real major_radius = 1.25;
  const real minor_radius = 0.35;
  const TorusFrame frame = make_torus_frame(vec3::Zero(), vec3::UnitZ());
  TorusMesh torus =
      make_offset_torus_shell(48, 24, major_radius, minor_radius, frame);
  GAUDI_ASSERT(torus.shell != nullptr);
  asawa::shell::shell &M = *torus.shell;
  const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);
  const real l0 = 2.0 * asawa::shell::avg_length(M, x);

  const int probe = 17;
  const vec3 pov = x[probe];
  const vec3 true_normal = canonical_torus_normal(pov, major_radius);
  const vec3 perturbed_normal =
      (true_normal + 0.35 * frame.x_axis).normalized();
  const std::vector<vec3> p_fit = {pov};
  const std::vector<vec3> n_fit = {perturbed_normal};

  auto foot_alignment = [&](const albers::vec14 &Q,
                            const vec3 &target_normal) -> real {
    const vec3 g = albers::darboux_grad(Q, vec3::Zero());
    if (!g.allFinite() || g.norm() < 1e-12) {
      return -1.0;
    }
    vec3 gn = g.normalized();
    if (gn.dot(target_normal) < 0.0) {
      gn *= -1.0;
    }
    return gn.dot(target_normal);
  };

  const std::vector<albers::vec14> Q_none =
      calder::darboux_cyclide_normal_constrained_convexity(M, p_fit, n_fit, l0,
                                                           3.0, 0.0);
  const std::vector<albers::vec14> Q_tiny =
      calder::darboux_cyclide_normal_constrained_convexity(M, p_fit, n_fit, l0,
                                                           3.0, 1e-12);
  const std::vector<albers::vec14> Q_strong =
      calder::darboux_cyclide_normal_constrained_convexity(M, p_fit, n_fit, l0,
                                                           3.0, 1.0);

  GAUDI_ASSERT(Q_none.size() == 1);
  GAUDI_ASSERT(Q_tiny.size() == 1);
  GAUDI_ASSERT(Q_strong.size() == 1);
  GAUDI_EXPECT(Q_none[0].allFinite());
  GAUDI_EXPECT(Q_tiny[0].allFinite());
  GAUDI_EXPECT(Q_strong[0].allFinite());

  const real align_none = foot_alignment(Q_none[0], perturbed_normal);
  const real align_tiny = foot_alignment(Q_tiny[0], perturbed_normal);
  const real align_strong = foot_alignment(Q_strong[0], perturbed_normal);
  const real align_strong_true = foot_alignment(Q_strong[0], true_normal);

  std::cerr << "\n[calder_darboux_cyclide_convexity_fit_w0_foot_normal]"
            << " align_w0=0: " << align_none
            << " align_w0=1e-12: " << align_tiny
            << " align_w0=1: " << align_strong
            << " strong_vs_true: " << align_strong_true << "\n";

  GAUDI_EXPECT((Q_none[0] - Q_tiny[0]).norm() < 1e-10);
  GAUDI_EXPECT(align_strong > align_none + 1e-3);
  GAUDI_EXPECT(align_strong > 0.9);
}

} // namespace test
} // namespace gaudi

#endif // __GAUDI_DARBOUX_CYCLIDE_TESTS_HPP__
