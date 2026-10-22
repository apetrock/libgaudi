#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/arp/pairwise_tests.hpp"
#include "gaudi/arp/ray_simplex.hpp"
#include "gaudi/arp/ray_traversal.hpp"
#include "gaudi/duchamp/pierce_cylinder_medial.hpp"
#include "gaudi/hepworth/constraints/slab/slab_pierce_stencil.hpp"
#include "gaudi/hepworth/constraints/slab/slab_prism_init.hpp"

#include <algorithm>
#include <cstring>
#include <iostream>
#include <limits>
#include <vector>

namespace {

std::vector<gaudi::index_t>
demo_load_verts(const std::vector<gaudi::vec3> &x) {
  std::vector<gaudi::index_t> loads;
  if (x.empty())
    return loads;
  const int n = static_cast<int>(x.size());
  gaudi::index_t top = 0;
  gaudi::real y_min = x[0].y();
  gaudi::real y_max = x[0].y();
  for (int i = 0; i < n; ++i) {
    y_min = std::min(y_min, x[static_cast<size_t>(i)].y());
    y_max = std::max(y_max, x[static_cast<size_t>(i)].y());
    if (x[static_cast<size_t>(i)].y() > x[static_cast<size_t>(top)].y())
      top = static_cast<gaudi::index_t>(i);
  }
  loads.push_back(top);
  const gaudi::real y_cut = y_min + 0.72 * (y_max - y_min);
  gaudi::index_t side = top;
  gaudi::real best_r = -1.0;
  for (int i = 0; i < n; ++i) {
    const gaudi::vec3 &xi = x[static_cast<size_t>(i)];
    if (xi.y() < y_cut)
      continue;
    gaudi::vec3 r = xi;
    r.y() = 0.0;
    const gaudi::real rn = r.norm();
    if (rn > best_r && static_cast<gaudi::index_t>(i) != top) {
      best_r = rn;
      side = static_cast<gaudi::index_t>(i);
    }
  }
  if (side != top)
    loads.push_back(side);
  return loads;
}

void print_ray_diagnosis(
    const gaudi::hepworth::block::slab_ray_pierce_diagnosis &diag) {
  std::cout << "    fat_candidates=" << diag.fat_ray_candidates
            << " self_incident=" << diag.fat_ray_self_incident
            << " fat_valid=" << diag.fat_ray_valid
            << " brute_valid(-n/+n)=" << diag.brute_valid_minus_n << "/"
            << diag.brute_valid_plus_n << " miss_detail="
            << gaudi::hepworth::block::ray_miss_detail_name(diag.miss_detail)
            << " antipodal_vert=" << diag.antipodal_vert
            << " dot=" << diag.antipodal_dot
            << " spine=" << diag.antipodal_spine << "\n";
}

struct pierce_outlier_row {
  int foot = -1;
  gaudi::real travel = 0.0;
  gaudi::real ray_gap = 0.0;
  gaudi::real surface_gap = 0.0;
  gaudi::real claim_gap = 0.0;
  gaudi::real t_tie_ratio = 0.0;
  gaudi::index_t hit_face = -1;
};

void audit_pierce_outliers(gaudi::asawa::shell::shell &shell,
                           std::vector<gaudi::vec3> &x, const char *label,
                           bool demo_mesh) {
  (void)demo_mesh;
  gaudi::arp::bvh_tree<3>::ptr tree =
      gaudi::hepworth::block::make_shell_face_bvh(shell, x);
  if (!tree) {
    return;
  }
  const std::vector<gaudi::index_t> face_map = shell.get_face_map(true);
  const gaudi::hepworth::block::slab_pierce_config cfg;
  const gaudi::real avg = gaudi::asawa::shell::avg_length(shell, x);
  const gaudi::real t_min =
      std::max(cfg.t_min_scale * avg, gaudi::real(1e-8));
  const gaudi::real min_spine =
      std::max(cfg.self_skip_tol, cfg.t_min_scale * avg);
  const gaudi::real slop =
      std::max(cfg.fat_ray_radius_scale * avg, gaudi::real(1e-10));
  const gaudi::real lift =
      std::max(t_min * gaudi::real(2.0), cfg.origin_lift_frac * avg);

  std::vector<pierce_outlier_row> rows;

  for (int iv = 0; iv < static_cast<int>(shell.vert_count()); ++iv) {
    const auto foot = gaudi::asawa::shell::vert_id(iv);
    gaudi::vec3 n_out = gaudi::asawa::shell::vert_normal(shell, foot, x);
    if (n_out.squaredNorm() < 1e-24) {
      continue;
    }
    n_out.normalize();
    const gaudi::vec3 origin = x[static_cast<size_t>(iv)];
    const gaudi::arp::Ray ray(origin + n_out * lift, -n_out, t_min, cfg.t_max);

    const auto pierce = gaudi::hepworth::block::ray_pierce_slab_shell(
        shell, foot, origin, n_out, x, *tree, face_map, cfg, t_min, avg);
    if (!pierce) {
      continue;
    }

    pierce_outlier_row row;
    row.foot = iv;
    row.hit_face = static_cast<gaudi::index_t>(pierce->fi);
    row.travel = (pierce->hit.point - origin).norm();
    row.ray_gap = gaudi::arp::ray_distance_to_point(ray, pierce->hit.point);

    const std::array<gaudi::vec3, 1> query = {pierce->hit.point};
    const auto nearest = tree->get_nearest(query, avg * gaudi::real(10.0));
    if (!nearest.empty()) {
      const gaudi::index_t leaf = nearest.back();
      const auto tri = (*tree->permuted_data_view_)[leaf];
      row.surface_gap =
          gaudi::arp::test_point_tri(query, std::array<gaudi::vec3, 3>{
                                                tri[0], tri[1], tri[2]});
    }

    const auto claim_tri = shell.get_tri(pierce->fi);
    row.claim_gap = gaudi::arp::test_point_tri(
        query, std::array<gaudi::vec3, 3>{
                   x[static_cast<size_t>(claim_tri[0])],
                   x[static_cast<size_t>(claim_tri[1])],
                   x[static_cast<size_t>(claim_tri[2])]});

    const std::vector<gaudi::arp::ray_hit<3>> all_hits =
        gaudi::arp::ray_slop_pierce_tree(*tree, ray, slop);
    gaudi::real min_t = std::numeric_limits<gaudi::real>::infinity();
    for (const auto &h : all_hits) {
      min_t = std::min(min_t, h.t_ray);
    }
    gaudi::real min_travel = std::numeric_limits<gaudi::real>::infinity();
    for (const auto &h : all_hits) {
      if (h.t_ray > min_t + 1e-7) {
        continue;
      }
      const auto mapped =
          gaudi::hepworth::block::slab_pierce_hit_from_bvh(shell, face_map, h);
      if (!mapped) {
        continue;
      }
      if (gaudi::hepworth::block::slab_pierce_hit_is_self(
              shell, foot, *mapped, origin, x, min_spine, slop, ray)) {
        continue;
      }
      min_travel = std::min(
          min_travel, (mapped->hit.point - origin).norm());
    }
    if (std::isfinite(min_travel) && min_travel > 1e-12) {
      row.t_tie_ratio = row.travel / min_travel;
    }

    rows.push_back(row);
  }

  std::sort(rows.begin(), rows.end(),
            [](const pierce_outlier_row &a, const pierce_outlier_row &b) {
              if (a.t_tie_ratio != b.t_tie_ratio) {
                return a.t_tie_ratio > b.t_tie_ratio;
              }
              return a.ray_gap > b.ray_gap;
            });

  std::cout << label << " pierce_outliers n=" << rows.size()
            << " (sorted by t_tie_ratio, ray_gap)\n";
  std::cout << "  foot  travel   ray_gap  surf_gap claim_gap t_tie  hit_face\n";
  const int show = std::min(static_cast<int>(rows.size()), 12);
  for (int i = 0; i < show; ++i) {
    const pierce_outlier_row &r = rows[static_cast<size_t>(i)];
    std::cout << "  " << r.foot << "  " << r.travel << "  " << r.ray_gap << "  "
              << r.surface_gap << "  " << r.claim_gap << "  " << r.t_tie_ratio
              << "  " << r.hit_face << "\n";
  }
}

} // namespace

static void audit_mesh(const char *label, gaudi::asawa::shell::shell &shell,
                       std::vector<gaudi::vec3> &x, bool diagnose_failures) {
  gaudi::arp::bvh_tree<3>::ptr tree =
      gaudi::hepworth::block::make_shell_face_bvh(shell, x);
  if (!tree)
    return;
  const std::vector<gaudi::index_t> face_map = shell.get_face_map(true);
  const gaudi::hepworth::block::slab_pierce_config cfg;
  const gaudi::real avg = gaudi::asawa::shell::avg_length(shell, x);

  const gaudi::hepworth::block::slab_prism_init_stats stats =
      gaudi::hepworth::block::audit_slab_pierce_stencils(shell, x, *tree, face_map,
                                                         cfg, avg);

  size_t antipodal_ok = 0;
  for (int iv = 0; iv < static_cast<int>(shell.vert_count()); ++iv) {
    const gaudi::asawa::shell::VertId vid = gaudi::asawa::shell::vert_id(iv);
    const gaudi::real A = gaudi::asawa::shell::vert_area(shell, vid, x);
    if (A < 1e-8)
      continue;
    if (shell.get_one_ring(vid).empty())
      continue;
    const auto result = gaudi::hepworth::block::build_slab_pierce_stencil(
        shell, vid, x, *tree, face_map, cfg, avg);
    if (result.used_antipodal_fallback)
      ++antipodal_ok;
  }

  std::cout << label << " verts=" << shell.vert_count() << " avg_edge=" << avg
            << " active_feet=" << stats.active_feet
            << " slab_ok=" << stats.slab_count
            << " pierce_fail=" << stats.pierce_fail_count
            << " antipodal_ok=" << antipodal_ok << "\n";
  for (const auto &[reason, count] : stats.fail_by_reason)
    std::cout << "  "
              << gaudi::hepworth::block::pierce_failure_reason_name(reason)
              << "=" << count << "\n";

  if (!diagnose_failures)
    return;

  std::cout << "ray_miss diagnostics:\n";
  for (int iv = 0; iv < static_cast<int>(shell.vert_count()); ++iv) {
    const gaudi::asawa::shell::VertId vid = gaudi::asawa::shell::vert_id(iv);
    const gaudi::real A = gaudi::asawa::shell::vert_area(shell, vid, x);
    if (A < 1e-8)
      continue;
    if (shell.get_one_ring(vid).empty())
      continue;
    const auto result = gaudi::hepworth::block::build_slab_pierce_stencil(
        shell, vid, x, *tree, face_map, cfg, avg);
    if (result.stencil.has_slab)
      continue;
    if (result.reason != gaudi::hepworth::block::pierce_failure_reason::ray_miss)
      continue;
    const gaudi::vec3 n =
        gaudi::asawa::shell::vert_normal(shell, vid, x).normalized();
    std::cout << "  foot=" << iv << " y=" << x[static_cast<size_t>(iv)].y()
              << " valence=" << shell.get_one_ring(vid).size()
              << " miss_detail="
              << gaudi::hepworth::block::ray_miss_detail_name(result.miss_detail)
              << "\n";
    print_ray_diagnosis(gaudi::hepworth::block::diagnose_slab_ray_pierce(
        shell, vid, x, *tree, face_map, cfg, avg));
  }
}

int main(int argc, char **argv) {
  using namespace gaudi;
  const bool diagnose = argc > 1 && std::strcmp(argv[1], "diagnose") == 0;
  const bool outliers = argc > 1 && std::strcmp(argv[1], "outliers") == 0;
  const bool demo_mesh =
      outliers && argc > 2 && std::strcmp(argv[2], "demo") == 0;
  const bool sphere_only =
      (argc > 1 && std::strcmp(argv[1], "sphere") == 0) ||
      (argc > 2 && std::strcmp(argv[2], "sphere") == 0);
  if (sphere_only) {
    auto M = asawa::shell::load_sphere(1.0, 8, 6);
    if (!M)
      return 1;
    asawa::shell::triangulate(*M);
    std::vector<vec3> x = asawa::const_get_vec_data(*M, 0);
    audit_mesh("sphere", *M, x, diagnose);
    return 0;
  }

  if (outliers) {
    auto M = asawa::shell::load_bunny();
    if (!M) {
      return 1;
    }
    asawa::shell::triangulate(*M);
    std::vector<vec3> x = asawa::const_get_vec_data(*M, 0);
    if (demo_mesh) {
      duchamp::normalize_medial_demo_mesh(*M);
      x = asawa::const_get_vec_data(*M, 0);
      audit_pierce_outliers(*M, x, "demo bunny", true);
    } else {
      asawa::center(x);
      audit_pierce_outliers(*M, x, "centered bunny", false);
    }
    return 0;
  }

  auto M = asawa::shell::load_bunny();
  if (!M)
    return 1;
  asawa::shell::triangulate(*M);
  std::vector<vec3> x = asawa::const_get_vec_data(*M, 0);
  asawa::center(x);

  arp::bvh_tree<3>::ptr tree = hepworth::block::make_shell_face_bvh(*M, x);
  if (!tree)
    return 1;
  const std::vector<index_t> face_map = M->get_face_map(true);
  const hepworth::block::slab_pierce_config cfg;
  const gaudi::real avg = asawa::shell::avg_length(*M, x);

  const hepworth::block::slab_prism_init_stats stats =
      hepworth::block::audit_slab_pierce_stencils(*M, x, *tree, face_map, cfg,
                                                  avg);

  std::cout << "bunny verts=" << M->vert_count() << " avg_edge=" << avg
            << " active_feet=" << stats.active_feet
            << " slab_ok=" << stats.slab_count
            << " pierce_fail=" << stats.pierce_fail_count << "\n";
  for (const auto &[reason, count] : stats.fail_by_reason)
    std::cout << "  " << hepworth::block::pierce_failure_reason_name(reason)
              << "=" << count << "\n";

  const std::vector<index_t> loads = demo_load_verts(x);
  std::cout << "demo load verts:\n";
  for (index_t vi : loads) {
    const auto result = hepworth::block::build_slab_pierce_stencil(
        *M, asawa::shell::vert_id(static_cast<int>(vi)), x, *tree, face_map, cfg,
        avg);
    std::cout << "  load_vert=" << vi
              << " has_slab=" << (result.stencil.has_slab ? "yes" : "NO")
              << " reason=" << hepworth::block::pierce_failure_reason_name(result.reason)
              << " y=" << x[static_cast<size_t>(vi)].y()
              << " hit_face=" << result.hit_face << " ray_t=" << result.ray_t
              << " antipodal=" << (result.used_antipodal_fallback ? "yes" : "no")
              << "\n";
  }

  {
    hepworth::block::slab_pierce_config no_ant = cfg;
    no_ant.antipodal_fallback = false;
    const gaudi::real t_min =
        std::max(no_ant.t_min_scale * avg, gaudi::real(1e-8));
    int pierced = 0;
    gaudi::real travel_sum = 0, max_travel = 0;
    gaudi::real spine_sum = 0, max_spine = 0;
    for (int iv = 0; iv < static_cast<int>(M->vert_count()); ++iv) {
      const auto foot = gaudi::asawa::shell::vert_id(iv);
      gaudi::vec3 n_out = gaudi::asawa::shell::vert_normal(*M, foot, x);
      if (n_out.squaredNorm() < 1e-24)
        continue;
      n_out.normalize();
      const auto hit = hepworth::block::ray_pierce_slab_shell(
          *M, foot, x[static_cast<size_t>(iv)], n_out, x, *tree, face_map, no_ant,
          t_min, avg);
      if (!hit)
        continue;
      ++pierced;
      const gaudi::real travel =
          (hit->hit.point - x[static_cast<size_t>(iv)]).norm();
      travel_sum += travel;
      max_travel = std::max(max_travel, travel);
      const auto tri = M->get_tri(hit->fi);
      const gaudi::index_t pierce_v =
          hepworth::block::closest_triangle_vertex_excluding(
              tri, x, hit->hit.point, foot);
      if (pierce_v >= 0) {
        const gaudi::real spine =
            (x[static_cast<size_t>(pierce_v)] - x[static_cast<size_t>(iv)])
                .norm();
        spine_sum += spine;
        max_spine = std::max(max_spine, spine);
      }
    }
    std::cout << "no_antipodal ray pierce: pierced=" << pierced
              << " avg_travel="
              << (pierced > 0 ? travel_sum / gaudi::real(pierced) : 0)
              << " max_travel=" << max_travel
              << " avg_spine="
              << (pierced > 0 ? spine_sum / gaudi::real(pierced) : 0)
              << " max_spine=" << max_spine << "\n";
  }

  std::cout << "all pierce failures:\n";
  for (int iv = 0; iv < static_cast<int>(M->vert_count()); ++iv) {
    const asawa::shell::VertId vid = asawa::shell::vert_id(iv);
    const gaudi::real A = asawa::shell::vert_area(*M, vid, x);
    if (A < 1e-8)
      continue;
    if (M->get_one_ring(vid).empty())
      continue;
    const auto result = hepworth::block::build_slab_pierce_stencil(
        *M, vid, x, *tree, face_map, cfg, avg);
    if (result.stencil.has_slab)
      continue;
    const gaudi::vec3 n =
        asawa::shell::vert_normal(*M, vid, x).normalized();
    std::cout << "  foot=" << iv
              << " reason=" << hepworth::block::pierce_failure_reason_name(result.reason)
              << " y=" << x[static_cast<size_t>(iv)].y()
              << " n_in=" << (-n).transpose() << " valence="
              << M->get_one_ring(vid).size()
              << " miss_detail="
              << hepworth::block::ray_miss_detail_name(result.miss_detail)
              << "\n";
    if (diagnose &&
        result.reason == hepworth::block::pierce_failure_reason::ray_miss)
      print_ray_diagnosis(hepworth::block::diagnose_slab_ray_pierce(
          *M, vid, x, *tree, face_map, cfg, avg));
  }
  return 0;
}
