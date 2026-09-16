#ifndef __GAUDI_DUCHAMP_PIERCE_CYLINDER_MEDIAL_HPP__
#define __GAUDI_DUCHAMP_PIERCE_CYLINDER_MEDIAL_HPP__

#include "gaudi/albers/line_cylinder.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/calder/least_squares_fit.hpp"
#include "gaudi/hepworth/constraints/slab/slab_pierce_stencil.hpp"

#include <cmath>
#include <iostream>
#include <vector>

namespace gaudi {
namespace duchamp {

struct pierce_cylinder_medial_params {
  hepworth::block::slab_pierce_config pierce;
  real axis_display_scale = 8.0;
  bool draw_rejected = true;
  real miss_ray_scale = 0.35;
  /// Gaussian MLS σ = mls_sigma_pierce_scale * (travel / 2).
  real mls_sigma_pierce_scale = 1.0;
  /// Floor on σ relative to average edge length (tiny pierces).
  real mls_sigma_min_edge_scale = 1.0;
  real mls_fit_p = 3.0;
};

struct pierce_cylinder_medial_stats {
  int total = 0;
  int accepted = 0;
  int rejected = 0;
  int pierce_miss_count = 0;
  int fit_fail_count = 0;
  int antipodal_count = 0;
  int mls_fit_count = 0;
  int plucker_fallback_count = 0;
  real avg_travel = 0.0;
  real avg_center_mid_drift = 0.0;
  real avg_center_mid_drift_frac = 0.0;
};

enum class pierce_cylinder_status { pierce_miss, fit_fail, accepted };

struct pierce_cylinder_sample {
  int vertex = -1;
  int pierce_vert = -1;
  vec3 foot = vec3::Zero();
  vec3 foot_normal = vec3::Zero();
  vec3 pierce = vec3::Zero();
  vec3 midpoint = vec3::Zero();
  vec3 center = vec3::Zero();
  vec3 axis = vec3::Zero();
  real travel = 0.0;
  real radius = 0.0;
  pierce_cylinder_status status = pierce_cylinder_status::pierce_miss;
  bool accepted = false;
  bool used_antipodal = false;
  bool used_plucker_fallback = false;
};

inline pierce_cylinder_medial_params default_pierce_cylinder_medial_params() {
  pierce_cylinder_medial_params params;
  params.pierce.antipodal_fallback = false;
  return params;
}

inline real pierce_mls_sigma_floor(real avg_edge_length,
                                   const pierce_cylinder_medial_params &params) {
  return params.mls_sigma_min_edge_scale *
         std::max(avg_edge_length, real(1e-12));
}

/// Gaussian MLS length scale from half pierce travel (local thickness prior).
inline real pierce_mls_gaussian_sigma(real travel, real avg_edge_length,
                                      const pierce_cylinder_medial_params &params) {
  const real half_pierce = 0.5 * std::max(travel, real(0.0));
  const real sigma = params.mls_sigma_pierce_scale * half_pierce;
  return std::max(sigma, pierce_mls_sigma_floor(avg_edge_length, params));
}

inline void print_pierce_cylinder_medial_stats(
    std::ostream &os, const pierce_cylinder_medial_stats &stats) {
  os << "medial axis estimation: total=" << stats.total
     << " accepted=" << stats.accepted << " rejected=" << stats.rejected
     << " pierce_miss=" << stats.pierce_miss_count
     << " fit_fail=" << stats.fit_fail_count
     << " antipodal=" << stats.antipodal_count
     << " mls_fit=" << stats.mls_fit_count
     << " plucker_fallback=" << stats.plucker_fallback_count
     << " avg_travel=" << stats.avg_travel
     << " avg_center_mid_drift=" << stats.avg_center_mid_drift
     << " avg_center_mid_drift_frac=" << stats.avg_center_mid_drift_frac
     << std::endl;
}

inline void normalize_medial_demo_mesh(asawa::shell::shell &M,
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
    p = (p - center) * scale;
  }
  std::cerr << "normalized mesh radius=" << target_radius << " scale=" << scale
            << std::endl;
}

struct inward_pierce_vertex_result {
  index_t pierce_vert = -1;
  vec3 foot_pos = vec3::Zero();
  vec3 pierce_pos = vec3::Zero();
  real travel = 0.0;
  bool ok = false;
  bool used_antipodal = false;
};

inline vec3 pierce_midpoint_proxy(const vec3 &foot_pos, const vec3 &pierce_pos) {
  return 0.5 * (foot_pos + pierce_pos);
}

inline inward_pierce_vertex_result pierce_inward_vertex(
    const asawa::shell::shell &shell, asawa::shell::VertId foot,
    const std::vector<vec3> &x, const arp::bvh_tree<3> &face_tree,
    const std::vector<index_t> &face_map,
    const hepworth::block::slab_pierce_config &cfg, real avg_edge_length) {
  inward_pierce_vertex_result out;
  out.foot_pos = x[static_cast<size_t>(foot)];

  vec3 n_out = asawa::shell::vert_normal(shell, foot, x);
  if (n_out.squaredNorm() < 1e-24) {
    return out;
  }
  n_out.normalize();

  const real t_min =
      std::max(cfg.t_min_scale * avg_edge_length, real(1e-8));
  const vec3 origin = out.foot_pos;
  const real min_spine =
      std::max(cfg.self_skip_tol, cfg.t_min_scale * avg_edge_length);

  index_t pierce = -1;
  const auto pierce_hit = hepworth::block::ray_pierce_slab_shell(
      shell, foot, origin, n_out, x, face_tree, face_map, cfg, t_min,
      avg_edge_length);
  if (pierce_hit) {
    const auto tri = shell.get_tri(pierce_hit->fi);
    pierce = hepworth::block::closest_triangle_vertex_excluding(
        tri, x, pierce_hit->hit.point, foot);
    if (pierce < 0) {
      return out;
    }
    out.pierce_vert = pierce;
    out.pierce_pos = x[static_cast<size_t>(pierce)];
    out.travel = (out.pierce_pos - origin).norm();
  } else if (cfg.antipodal_fallback) {
    const auto ant = hepworth::block::antipodal_pierce_vertex(
        shell, foot, x, n_out, min_spine, cfg.antipodal_min_dot);
    if (!ant) {
      return out;
    }
    pierce = *ant;
    out.used_antipodal = true;
    out.pierce_vert = pierce;
    out.pierce_pos = x[static_cast<size_t>(pierce)];
    out.travel = (out.pierce_pos - origin).norm();
  } else {
    return out;
  }

  if (pierce == static_cast<index_t>(foot)) {
    return out;
  }

  if (out.travel < min_spine) {
    return out;
  }

  out.pierce_vert = pierce;
  out.ok = true;
  return out;
}

/// Optional drift gate (disabled in the live path — we keep raw MLS centers).
inline bool accept_medial_center(const vec3 &center, const vec3 &pov,
                                 const vec3 &pierce_pos, real spine_len) {
  if (!center.allFinite()) {
    return false;
  }
  const real mid_drift = (center - pov).norm();
  const real max_drift =
      spine_len > 1e-12 ? real(0.35) * spine_len : real(1e-3);
  if (mid_drift > max_drift) {
    return false;
  }
  if (spine_len > 1e-12) {
    const real pierce_drift = (center - pierce_pos).norm();
    if (mid_drift >= pierce_drift) {
      return false;
    }
  }
  return true;
}

inline bool extract_mls_cylinder_fit(const calder::vec7 &fit,
                                     const vec3 &spine_dir, vec3 *center_out,
                                     vec3 *axis_out, real *radius_out) {
  if (!center_out || !axis_out || !radius_out) {
    return false;
  }
  // calder::cylinder vec7 layout: [center, cylinder_axis, radius]
  const vec3 center = fit.segment(0, 3);
  vec3 axis = fit.segment(3, 3);
  const real radius = fit[6];
  if (!center.allFinite() || !axis.allFinite() || !std::isfinite(radius)) {
    return false;
  }
  if (axis.squaredNorm() < 1e-24 || radius <= 1e-8) {
    return false;
  }
  axis.normalize();
  if (spine_dir.squaredNorm() > 1e-24 &&
      axis.dot(spine_dir.normalized()) < 0.0) {
    axis = -axis;
  }
  *center_out = center;
  *axis_out = axis;
  *radius_out = radius;
  return true;
}

inline bool fit_pierce_cylinder_plucker(const asawa::shell::shell &shell,
                                        const std::vector<vec3> &x,
                                        asawa::shell::VertId foot,
                                        asawa::shell::VertId pierce,
                                        const vec3 &spine_dir, vec3 *center_out,
                                        vec3 *axis_out, real *radius_out) {
  if (!center_out || !axis_out || !radius_out) {
    return false;
  }

  const vec3 foot_pos = x[static_cast<size_t>(foot)];
  const vec3 pierce_pos = x[static_cast<size_t>(pierce)];

  vec3 n_foot = asawa::shell::vert_normal(shell, foot, x);
  vec3 n_pierce = asawa::shell::vert_normal(shell, pierce, x);
  if (n_foot.squaredNorm() < 1e-24 || n_pierce.squaredNorm() < 1e-24) {
    return false;
  }
  n_foot.normalize();
  n_pierce.normalize();

  const vec6 Q = albers::fit_normal_aligned_line(
      {foot_pos, pierce_pos}, {n_foot, n_pierce});
  vec3 axis = albers::plucker_line_direction(Q);
  if (!axis.allFinite() || axis.squaredNorm() < 1e-24) {
    return false;
  }
  axis.normalize();
  if (spine_dir.squaredNorm() > 1e-24 &&
      axis.dot(spine_dir.normalized()) < 0.0) {
    axis = -axis;
  }

  const vec3 center = albers::plucker_line_closest_point(Q, foot_pos);
  if (!center.allFinite()) {
    return false;
  }

  const real radius = (center - foot_pos).norm();
  if (!(radius > 1e-8) || !std::isfinite(radius)) {
    return false;
  }

  *center_out = center;
  *axis_out = axis;
  *radius_out = radius;
  return true;
}

struct inward_pierce_probe_result {
  vec3 foot_pos = vec3::Zero();
  vec3 foot_normal = vec3::Zero();
  vec3 pierce_pos = vec3::Zero();
  real travel = 0.0;
  bool ok = false;
};

inline inward_pierce_probe_result pierce_inward_probe(
    const asawa::shell::shell &shell, const vec3 &origin, const vec3 &n_out,
    asawa::shell::FaceId source_face, const std::vector<vec3> &x,
    const arp::bvh_tree<3> &face_tree, const std::vector<index_t> &face_map,
    const hepworth::block::slab_pierce_config &cfg, real avg_edge_length) {
  inward_pierce_probe_result out;
  out.foot_pos = origin;
  if (!n_out.allFinite() || n_out.squaredNorm() < 1e-24) {
    return out;
  }
  out.foot_normal = n_out.normalized();
  const real t_min = std::max(cfg.t_min_scale * avg_edge_length, real(1e-8));
  const real slop =
      std::max(cfg.fat_ray_radius_scale * avg_edge_length, real(1e-10));
  const real min_spine =
      std::max(cfg.self_skip_tol, cfg.t_min_scale * avg_edge_length);
  const real lift =
      std::max(t_min * real(2.0), cfg.origin_lift_frac * avg_edge_length);
  const arp::Ray ray(origin + out.foot_normal * lift, -out.foot_normal, t_min,
                     cfg.t_max);

  std::vector<arp::ray_hit<3>> hits =
      arp::ray_slop_pierce_tree(face_tree, ray, slop);
  arp::ray_hits_make_min_heap(hits, &ray);
  while (!hits.empty()) {
    const arp::ray_hit<3> h = arp::ray_hits_pop_min(hits, &ray);
    const auto mapped =
        hepworth::block::slab_pierce_hit_from_bvh(shell, face_map, h);
    if (!mapped) {
      continue;
    }
    if (source_face >= 0 && mapped->fi == source_face) {
      continue;
    }
    if ((mapped->hit.point - origin).norm() < min_spine) {
      continue;
    }
    if (arp::ray_distance_to_point(ray, mapped->hit.point) > slop + real(1e-4)) {
      continue;
    }
    out.pierce_pos = mapped->hit.point;
    out.travel = (out.pierce_pos - origin).norm();
    out.ok = out.travel >= min_spine;
    if (out.ok) {
      return out;
    }
  }
  return out;
}

inline std::vector<pierce_cylinder_sample> compute_pierce_cylinder_medial_at(
    asawa::shell::shell &shell, const std::vector<vec3> &p_pov,
    const std::vector<vec3> &n_pov, const pierce_cylinder_medial_params &params,
    pierce_cylinder_medial_stats *stats_out = nullptr) {
  std::vector<pierce_cylinder_sample> out;
  const std::vector<vec3> &x = asawa::const_get_vec_data(shell, 0);
  if (x.empty() || p_pov.empty() || p_pov.size() != n_pov.size()) {
    return out;
  }

  const real avg_len = asawa::shell::avg_length(shell, x);
  arp::bvh_tree<3>::ptr face_tree = hepworth::block::make_shell_face_bvh(shell, x);
  if (!face_tree) {
    return out;
  }
  const std::vector<index_t> face_map = shell.get_face_map(true);

  pierce_cylinder_medial_stats stats;
  stats.total = static_cast<int>(p_pov.size());
  out.resize(static_cast<size_t>(stats.total));

  struct pierce_fit_request {
    int vi = -1;
    vec3 mid = vec3::Zero();
    vec3 pierce_pos = vec3::Zero();
    vec3 spine_dir = vec3::Zero();
    vec3 n_probe = vec3::Zero();
    real spine_len = 0.0;
  };
  std::vector<pierce_fit_request> fit_requests;

  for (int vi = 0; vi < stats.total; ++vi) {
    pierce_cylinder_sample sample;
    sample.vertex = vi;
    sample.foot = p_pov[static_cast<size_t>(vi)];
    if (n_pov[static_cast<size_t>(vi)].squaredNorm() > 1e-24) {
      sample.foot_normal = n_pov[static_cast<size_t>(vi)].normalized();
    }
    const asawa::shell::FaceId source_face = asawa::shell::face_id(vi);
    const auto pierce = pierce_inward_probe(
        shell, sample.foot, sample.foot_normal, source_face, x, *face_tree,
        face_map, params.pierce, avg_len);
    if (!pierce.ok) {
      sample.status = pierce_cylinder_status::pierce_miss;
      out[static_cast<size_t>(vi)] = sample;
      ++stats.pierce_miss_count;
      continue;
    }

    sample.pierce = pierce.pierce_pos;
    sample.travel = pierce.travel;
    sample.midpoint = pierce_midpoint_proxy(sample.foot, sample.pierce);

    const vec3 spine = sample.pierce - sample.foot;
    pierce_fit_request req;
    req.vi = vi;
    req.mid = sample.midpoint;
    req.pierce_pos = sample.pierce;
    req.spine_dir = spine;
    req.n_probe = -sample.foot_normal;
    req.spine_len = sample.travel;
    fit_requests.push_back(req);
    out[static_cast<size_t>(vi)] = sample;
  }

  if (!fit_requests.empty()) {
    std::vector<vec3> mids;
    std::vector<vec3> fit_normals;
    std::vector<real> sigmas;
    mids.reserve(fit_requests.size());
    fit_normals.reserve(fit_requests.size());
    sigmas.reserve(fit_requests.size());
    for (const pierce_fit_request &req : fit_requests) {
      mids.push_back(req.mid);
      fit_normals.push_back(req.n_probe);
      sigmas.push_back(
          pierce_mls_gaussian_sigma(req.spine_len, avg_len, params));
    }

    const std::vector<calder::vec7> mls_fits = calder::cylinder_gaussian(
        shell, mids, fit_normals, sigmas, params.mls_fit_p);

    for (size_t fi = 0; fi < fit_requests.size(); ++fi) {
      const pierce_fit_request &req = fit_requests[fi];
      pierce_cylinder_sample &sample = out[static_cast<size_t>(req.vi)];

      vec3 center = vec3::Zero();
      vec3 axis = vec3::Zero();
      real radius = 0.0;
      const bool extracted = extract_mls_cylinder_fit(
          mls_fits[fi], req.spine_dir, &center, &axis, &radius);
      if (extracted) {
        // Keep raw MLS center/axis/radius — no drift miss filter.
        sample.center = center;
        sample.axis = axis;
        sample.radius = radius;
        sample.status = pierce_cylinder_status::accepted;
        sample.accepted = true;
        ++stats.mls_fit_count;
        ++stats.accepted;
        continue;
      }

      sample.status = pierce_cylinder_status::fit_fail;
      ++stats.fit_fail_count;
    }
  }

  real travel_sum = 0.0;
  real drift_sum = 0.0;
  real drift_frac_sum = 0.0;
  for (const pierce_cylinder_sample &sample : out) {
    if (!sample.accepted) {
      continue;
    }
    travel_sum += sample.travel;
    const real drift = (sample.center - sample.midpoint).norm();
    drift_sum += drift;
    if (sample.travel > 1e-12) {
      drift_frac_sum += drift / sample.travel;
    }
  }

  stats.rejected = stats.total - stats.accepted;
  if (stats.accepted > 0) {
    stats.avg_travel = travel_sum / real(stats.accepted);
    stats.avg_center_mid_drift = drift_sum / real(stats.accepted);
    stats.avg_center_mid_drift_frac = drift_frac_sum / real(stats.accepted);
  }
  if (stats_out) {
    *stats_out = stats;
  }
  return out;
}

inline std::vector<pierce_cylinder_sample> compute_pierce_cylinder_medial(
    asawa::shell::shell &shell, const pierce_cylinder_medial_params &params,
    pierce_cylinder_medial_stats *stats_out = nullptr) {
  std::vector<pierce_cylinder_sample> out;
  const std::vector<vec3> &x = asawa::const_get_vec_data(shell, 0);
  if (x.empty()) {
    return out;
  }

  const real avg_len = asawa::shell::avg_length(shell, x);
  arp::bvh_tree<3>::ptr face_tree =
      hepworth::block::make_shell_face_bvh(shell, x);
  if (!face_tree) {
    return out;
  }
  const std::vector<index_t> face_map = shell.get_face_map(true);

  pierce_cylinder_medial_stats stats;
  stats.total = static_cast<int>(shell.vert_count());
  out.resize(static_cast<size_t>(stats.total));

  struct pierce_fit_request {
    int vi = -1;
    vec3 mid = vec3::Zero();
    vec3 pierce_pos = vec3::Zero();
    vec3 spine_dir = vec3::Zero();
    vec3 n_pov = vec3::Zero();
    real spine_len = 0.0;
    asawa::shell::VertId pierce_vid = asawa::shell::vert_id(-1);
  };
  std::vector<pierce_fit_request> fit_requests;

  for (int vi = 0; vi < stats.total; ++vi) {
    pierce_cylinder_sample sample;
    sample.vertex = vi;
    sample.foot = x[static_cast<size_t>(vi)];

    vec3 n_foot = asawa::shell::vert_normal(shell, asawa::shell::vert_id(vi), x);
    if (n_foot.squaredNorm() > 1e-24) {
      sample.foot_normal = n_foot.normalized();
    }

    const auto pierce = pierce_inward_vertex(
        shell, asawa::shell::vert_id(vi), x, *face_tree, face_map,
        params.pierce, avg_len);
    if (!pierce.ok) {
      sample.status = pierce_cylinder_status::pierce_miss;
      out[static_cast<size_t>(vi)] = sample;
      ++stats.pierce_miss_count;
      continue;
    }

    sample.pierce = pierce.pierce_pos;
    sample.pierce_vert = static_cast<int>(pierce.pierce_vert);
    sample.travel = pierce.travel;
    sample.used_antipodal = pierce.used_antipodal;
    sample.midpoint = pierce_midpoint_proxy(pierce.foot_pos, pierce.pierce_pos);

    if (sample.foot_normal.squaredNorm() < 1e-24) {
      sample.status = pierce_cylinder_status::fit_fail;
      out[static_cast<size_t>(vi)] = sample;
      ++stats.fit_fail_count;
      continue;
    }

    const vec3 spine = pierce.pierce_pos - pierce.foot_pos;
    const vec3 n_pov = -sample.foot_normal;

    pierce_fit_request req;
    req.vi = vi;
    req.mid = sample.midpoint;
    req.pierce_pos = pierce.pierce_pos;
    req.spine_dir = spine;
    req.n_pov = n_pov;
    req.spine_len = pierce.travel;
    req.pierce_vid = asawa::shell::vert_id(pierce.pierce_vert);
    fit_requests.push_back(req);
    out[static_cast<size_t>(vi)] = sample;
  }

  if (!fit_requests.empty()) {
    std::vector<vec3> mids;
    std::vector<vec3> n_pov;
    std::vector<real> sigmas;
    mids.reserve(fit_requests.size());
    n_pov.reserve(fit_requests.size());
    sigmas.reserve(fit_requests.size());
    for (const pierce_fit_request &req : fit_requests) {
      mids.push_back(req.mid);
      n_pov.push_back(req.n_pov);
      sigmas.push_back(
          pierce_mls_gaussian_sigma(req.spine_len, avg_len, params));
    }

    const std::vector<calder::vec7> mls_fits =
        calder::cylinder_gaussian(shell, mids, n_pov, sigmas, params.mls_fit_p);

    for (size_t fi = 0; fi < fit_requests.size(); ++fi) {
      const pierce_fit_request &req = fit_requests[fi];
      pierce_cylinder_sample &sample = out[static_cast<size_t>(req.vi)];

      vec3 center = vec3::Zero();
      vec3 axis = vec3::Zero();
      real radius = 0.0;
      const bool extracted = extract_mls_cylinder_fit(
          mls_fits[fi], req.spine_dir, &center, &axis, &radius);
      if (extracted) {
        // Keep raw MLS center/axis/radius — no drift miss filter.
        sample.center = center;
        sample.axis = axis;
        sample.radius = radius;
        sample.status = pierce_cylinder_status::accepted;
        sample.accepted = true;
        ++stats.mls_fit_count;
        ++stats.accepted;
        if (sample.used_antipodal) {
          ++stats.antipodal_count;
        }
        continue;
      }

      vec3 p_center = vec3::Zero();
      vec3 p_axis = vec3::Zero();
      real p_radius = 0.0;
      if (fit_pierce_cylinder_plucker(shell, x, asawa::shell::vert_id(req.vi),
                                      req.pierce_vid, req.spine_dir, &p_center,
                                      &p_axis, &p_radius)) {
        sample.center = p_center;
        sample.axis = p_axis;
        sample.radius = p_radius;
        sample.used_plucker_fallback = true;
        sample.status = pierce_cylinder_status::accepted;
        sample.accepted = true;
        ++stats.plucker_fallback_count;
        ++stats.accepted;
        if (sample.used_antipodal) {
          ++stats.antipodal_count;
        }
        continue;
      }

      sample.status = pierce_cylinder_status::fit_fail;
      ++stats.fit_fail_count;
    }
  }

  real travel_sum = 0.0;
  real drift_sum = 0.0;
  real drift_frac_sum = 0.0;
  for (const pierce_cylinder_sample &sample : out) {
    if (!sample.accepted) {
      continue;
    }
    travel_sum += sample.travel;
    const real drift = (sample.center - sample.midpoint).norm();
    drift_sum += drift;
    if (sample.travel > 1e-12) {
      drift_frac_sum += drift / sample.travel;
    }
  }

  stats.rejected = stats.total - stats.accepted;
  if (stats.accepted > 0) {
    stats.avg_travel = travel_sum / real(stats.accepted);
    stats.avg_center_mid_drift = drift_sum / real(stats.accepted);
    stats.avg_center_mid_drift_frac = drift_frac_sum / real(stats.accepted);
  }
  if (stats_out) {
    *stats_out = stats;
  }
  return out;
}

} // namespace duchamp
} // namespace gaudi

#endif // __GAUDI_DUCHAMP_PIERCE_CYLINDER_MEDIAL_HPP__
