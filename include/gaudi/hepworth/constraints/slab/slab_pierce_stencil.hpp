#ifndef GAUDI_HEPWORTH_SLAB_PIERCE_STENCIL_HPP
#define GAUDI_HEPWORTH_SLAB_PIERCE_STENCIL_HPP

#include "gaudi/arp/arp.h"
#include "gaudi/arp/ray_traversal.hpp"
#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/smallvec.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace gaudi {
namespace hepworth {
namespace block {

struct slab_pierce_stencil {
  smallvec<index_t, 16> top_ids;
  smallvec<index_t, 16> bot_ids;
  smallvec<real, 16> top_weights;
  smallvec<real, 16> bot_weights;
  index_t foot = -1;
  index_t pierce_foot = -1;
  index_t hit_face = -1;
  bool has_slab = false;
};

struct slab_pierce_config {
  real t_min_scale = 1e-3;
  real t_max = 1e6;
  real self_skip_tol = 1e-6;
  /// Lift ray origin off the foot vertex along ±normal (fraction of avg edge).
  real origin_lift_frac = 0.02;
  /// When BVH pierce misses, pick the most opposite mesh vertex (closed shells).
  bool antipodal_fallback = false;
  real antipodal_min_dot = 0.65;
  /// Cylinder radius = fat_ray_radius_scale * avg_edge_length (ray slop).
  real fat_ray_radius_scale = 0.05;
};

enum class ray_miss_detail {
  none,
  bvh_no_hit,
  all_hits_self_incident,
  face_map_oob,
};

inline const char *ray_miss_detail_name(ray_miss_detail d) {
  switch (d) {
  case ray_miss_detail::none:
    return "none";
  case ray_miss_detail::bvh_no_hit:
    return "bvh_no_hit";
  case ray_miss_detail::all_hits_self_incident:
    return "all_hits_self_incident";
  case ray_miss_detail::face_map_oob:
    return "face_map_oob";
  }
  return "unknown";
}

enum class pierce_failure_reason {
  ok,
  empty_one_ring,
  top_cotan_mismatch,
  degenerate_normal,
  ray_miss,
  face_map_oob,
  self_face_hit,
  pierce_equals_foot,
  pierce_too_close,
  empty_bot_ring,
  bot_cotan_mismatch,
};

inline const char *pierce_failure_reason_name(pierce_failure_reason r) {
  switch (r) {
  case pierce_failure_reason::ok:
    return "ok";
  case pierce_failure_reason::empty_one_ring:
    return "empty_one_ring";
  case pierce_failure_reason::top_cotan_mismatch:
    return "top_cotan_mismatch";
  case pierce_failure_reason::degenerate_normal:
    return "degenerate_normal";
  case pierce_failure_reason::ray_miss:
    return "ray_miss";
  case pierce_failure_reason::face_map_oob:
    return "face_map_oob";
  case pierce_failure_reason::self_face_hit:
    return "self_face_hit";
  case pierce_failure_reason::pierce_equals_foot:
    return "pierce_equals_foot";
  case pierce_failure_reason::pierce_too_close:
    return "pierce_too_close";
  case pierce_failure_reason::empty_bot_ring:
    return "empty_bot_ring";
  case pierce_failure_reason::bot_cotan_mismatch:
    return "bot_cotan_mismatch";
  }
  return "unknown";
}

struct slab_pierce_result {
  slab_pierce_stencil stencil;
  pierce_failure_reason reason = pierce_failure_reason::ok;
  index_t hit_face = -1;
  real ray_t = -1.0;
  bool used_antipodal_fallback = false;
  ray_miss_detail miss_detail = ray_miss_detail::none;
};

struct slab_ray_pierce_diagnosis {
  int fat_ray_candidates = 0;
  int fat_ray_self_incident = 0;
  int fat_ray_valid = 0;
  int brute_valid_minus_n = 0;
  int brute_valid_plus_n = 0;
  ray_miss_detail miss_detail = ray_miss_detail::none;
  index_t antipodal_vert = -1;
  real antipodal_dot = -2.0;
  real antipodal_spine = 0.0;
};

inline bool shell_face_incident_to_vert(const asawa::shell::shell &shell,
                                        asawa::shell::FaceId fi,
                                        asawa::shell::VertId vi) {
  bool found = false;
  shell.const_for_each_face(fi, [&](asawa::shell::CornerId ci,
                                    const asawa::shell::shell &M) {
    if (M.vert(ci) == vi)
      found = true;
  });
  return found;
}

inline arp::bvh_tree<3>::ptr
make_shell_face_bvh(asawa::shell::shell &shell,
                    const std::vector<vec3> &x) {
  const std::vector<index_t> face_verts = shell.get_face_vert_ids(true);
  if (face_verts.empty())
    return {};
  arp::simplex_set<3> face_set(x, face_verts);
  return arp::bvh_tree<3>::create(face_set);
}

inline index_t closest_triangle_vertex(const std::array<asawa::shell::VertId, 3> &tri,
                                       const std::vector<vec3> &x,
                                       const vec3 &p) {
  index_t best = static_cast<index_t>(tri[0]);
  real best_d = (x[best] - p).squaredNorm();
  for (int k = 1; k < 3; ++k) {
    const index_t vi = static_cast<index_t>(tri[k]);
    const real d = (x[vi] - p).squaredNorm();
    if (d < best_d) {
      best_d = d;
      best = vi;
    }
  }
  return best;
}

inline index_t closest_triangle_vertex_excluding(
    const std::array<asawa::shell::VertId, 3> &tri, const std::vector<vec3> &x,
    const vec3 &p, asawa::shell::VertId exclude) {
  index_t best = -1;
  real best_d = std::numeric_limits<real>::infinity();
  for (int k = 0; k < 3; ++k) {
    const index_t vi = static_cast<index_t>(tri[k]);
    if (vi == static_cast<index_t>(exclude)) {
      continue;
    }
    const real d = (x[vi] - p).squaredNorm();
    if (d < best_d) {
      best_d = d;
      best = vi;
    }
  }
  return best;
}

struct slab_pierce_hit {
  asawa::shell::FaceId fi = asawa::shell::face_id(-1);
  arp::ray_hit<3> hit{};
};

inline std::optional<slab_pierce_hit>
slab_pierce_hit_from_bvh(const asawa::shell::shell &shell,
                         const std::vector<index_t> &face_map,
                         const arp::ray_hit<3> &hit) {
  if (!hit.valid || hit.id < 0 ||
      static_cast<size_t>(hit.id) >= face_map.size()) {
    return std::nullopt;
  }
  const asawa::shell::FaceId fi =
      asawa::shell::face_id(face_map[static_cast<size_t>(hit.id)]);
  return slab_pierce_hit{fi, hit};
}

inline bool slab_vert_in_one_ring(const asawa::shell::shell &shell,
                                  asawa::shell::VertId foot, index_t vi) {
  for (asawa::shell::VertId v : shell.get_one_ring(foot)) {
    if (static_cast<index_t>(v) == vi) {
      return true;
    }
  }
  return false;
}

/// Reject heap pop: foot-incident face, local cap, bad contact, micro spine.
inline bool slab_pierce_hit_is_self(const asawa::shell::shell &shell,
                                    asawa::shell::VertId foot,
                                    const slab_pierce_hit &mapped,
                                    const vec3 &origin,
                                    const std::vector<vec3> &x, real min_spine,
                                    real slop, const arp::Ray &ray) {
  if (shell_face_incident_to_vert(shell, mapped.fi, foot)) {
    return true;
  }
  if ((mapped.hit.point - origin).norm() < min_spine) {
    return true;
  }
  if (arp::ray_distance_to_point(ray, mapped.hit.point) > slop + real(1e-4)) {
    return true;
  }
  const auto tri = shell.get_tri(mapped.fi);
  const index_t pierce_v =
      closest_triangle_vertex_excluding(tri, x, mapped.hit.point, foot);
  if (pierce_v < 0 || pierce_v == static_cast<index_t>(foot)) {
    return true;
  }
  if (slab_vert_in_one_ring(shell, foot, pierce_v)) {
    return true;
  }
  if ((x[static_cast<size_t>(pierce_v)] - origin).norm() < min_spine) {
    return true;
  }
  return false;
}

inline std::optional<slab_pierce_hit>
ray_pierce_slab_shell(const asawa::shell::shell &shell,
                      asawa::shell::VertId foot, const vec3 &origin,
                      const vec3 &n_out, const std::vector<vec3> &x,
                      const arp::bvh_tree<3> &face_tree,
                      const std::vector<index_t> &face_map,
                      const slab_pierce_config &cfg, real t_min,
                      real avg_edge_length) {
  const real slop =
      std::max(cfg.fat_ray_radius_scale * avg_edge_length, real(1e-10));
  const real min_spine =
      std::max(cfg.self_skip_tol, cfg.t_min_scale * avg_edge_length);
  const real lift =
      std::max(t_min * real(2.0), cfg.origin_lift_frac * avg_edge_length);
  const arp::Ray ray(origin + n_out * lift, -n_out, t_min, cfg.t_max);

  std::vector<arp::ray_hit<3>> hits =
      arp::ray_slop_pierce_tree(face_tree, ray, slop);
  arp::ray_hits_make_min_heap(hits, &ray);
  while (!hits.empty()) {
    const arp::ray_hit<3> h = arp::ray_hits_pop_min(hits, &ray);
    const auto mapped = slab_pierce_hit_from_bvh(shell, face_map, h);
    if (!mapped) {
      continue;
    }
    if (slab_pierce_hit_is_self(shell, foot, *mapped, origin, x, min_spine,
                                slop, ray)) {
      continue;
    }
    return mapped;
  }
  return std::nullopt;
}

inline std::optional<index_t>
antipodal_pierce_vertex(const asawa::shell::shell &shell,
                        asawa::shell::VertId foot, const std::vector<vec3> &x,
                        const vec3 &n_out, real min_spine, real min_dot) {
  const vec3 o = x[static_cast<size_t>(foot)];
  index_t best = -1;
  real best_dot = min_dot;
  for (int iv = 0; iv < static_cast<int>(shell.vert_count()); ++iv) {
    if (iv == static_cast<int>(foot))
      continue;
    const vec3 d = x[static_cast<size_t>(iv)] - o;
    const real len = d.norm();
    if (len < min_spine)
      continue;
    const real dot = (-n_out).dot(d / len);
    if (dot > best_dot) {
      best_dot = dot;
      best = static_cast<index_t>(iv);
    }
  }
  if (best < 0)
    return std::nullopt;
  return best;
}

inline slab_ray_pierce_diagnosis
diagnose_slab_ray_pierce(const asawa::shell::shell &shell,
                         asawa::shell::VertId foot, const std::vector<vec3> &x,
                         const arp::bvh_tree<3> &face_tree,
                         const std::vector<index_t> &face_map,
                         const slab_pierce_config &cfg, real avg_edge_length) {
  slab_ray_pierce_diagnosis diag;
  vec3 n_out = asawa::shell::vert_normal(shell, foot, x);
  if (n_out.squaredNorm() < 1e-24)
    return diag;
  n_out.normalize();

  const real t_min =
      std::max(cfg.t_min_scale * avg_edge_length, real(1e-8));
  const real min_spine =
      std::max(cfg.self_skip_tol, cfg.t_min_scale * avg_edge_length);
  const vec3 origin = x[static_cast<size_t>(foot)];
  const real slop =
      std::max(cfg.fat_ray_radius_scale * avg_edge_length, real(1e-10));
  const real lift =
      std::max(t_min * real(2.0), cfg.origin_lift_frac * avg_edge_length);
  const arp::Ray ray(origin + n_out * lift, -n_out, t_min, cfg.t_max);

  const std::vector<arp::ray_hit<3>> raw =
      arp::ray_slop_pierce_tree(face_tree, ray, slop);
  diag.fat_ray_candidates = static_cast<int>(raw.size());
  for (const arp::ray_hit<3> &h : raw) {
    const auto mapped = slab_pierce_hit_from_bvh(shell, face_map, h);
    if (!mapped) {
      continue;
    }
    if (slab_pierce_hit_is_self(shell, foot, *mapped, origin, x, min_spine,
                                slop, ray)) {
      ++diag.fat_ray_self_incident;
      continue;
    }
    ++diag.fat_ray_valid;
  }

  for (size_t leaf = 0; leaf < face_map.size(); ++leaf) {
    const asawa::shell::FaceId fi =
        asawa::shell::face_id(face_map[leaf]);
    const auto tri = shell.get_tri(fi);
    for (int di = 0; di < 2; ++di) {
      const vec3 dir = (di == 0) ? -n_out : n_out;
      const arp::Ray ray(origin, dir, t_min, cfg.t_max);
      const auto h = arp::ray_pierce_triangle(
          ray, x[static_cast<size_t>(tri[0])], x[static_cast<size_t>(tri[1])],
          x[static_cast<size_t>(tri[2])], static_cast<index_t>(leaf));
      if (!h)
        continue;
      const slab_pierce_hit mapped{fi, *h};
      if (slab_pierce_hit_is_self(shell, foot, mapped, origin, x, min_spine,
                                  slop, ray))
        continue;
      if (di == 0)
        ++diag.brute_valid_minus_n;
      else
        ++diag.brute_valid_plus_n;
    }
  }

  if (const auto ant = antipodal_pierce_vertex(shell, foot, x, n_out, min_spine,
                                               cfg.antipodal_min_dot)) {
    diag.antipodal_vert = *ant;
    const vec3 d = x[static_cast<size_t>(*ant)] - origin;
    diag.antipodal_spine = d.norm();
    diag.antipodal_dot = (-n_out).dot(d / diag.antipodal_spine);
  }

  if (diag.fat_ray_valid > 0)
    diag.miss_detail = ray_miss_detail::none;
  else if (diag.fat_ray_candidates > 0 ||
           diag.fat_ray_self_incident > 0)
    diag.miss_detail = ray_miss_detail::all_hits_self_incident;
  else
    diag.miss_detail = ray_miss_detail::bvh_no_hit;
  return diag;
}

inline slab_pierce_result
build_slab_pierce_stencil(const asawa::shell::shell &shell,
                          asawa::shell::VertId foot,
                          const std::vector<vec3> &x,
                          const arp::bvh_tree<3> &face_tree,
                          const std::vector<index_t> &face_map,
                          const slab_pierce_config &cfg,
                          real avg_edge_length) {
  slab_pierce_result out;
  out.stencil.foot = static_cast<index_t>(foot);
  out.stencil.has_slab = false;

  const auto ring = shell.get_one_ring(foot);
  if (ring.empty()) {
    out.reason = pierce_failure_reason::empty_one_ring;
    return out;
  }

  out.stencil.top_ids.push_back(static_cast<index_t>(foot));
  for (asawa::shell::VertId v : ring)
    out.stencil.top_ids.push_back(static_cast<index_t>(v));

  const std::vector<real> top_w =
      asawa::shell::vert_cotan_weights(shell, foot, x);
  if (top_w.size() + 1 != out.stencil.top_ids.size()) {
    out.reason = pierce_failure_reason::top_cotan_mismatch;
    return out;
  }
  for (real w : top_w)
    out.stencil.top_weights.push_back(w);

  vec3 n_out = asawa::shell::vert_normal(shell, foot, x);
  if (n_out.squaredNorm() < 1e-24) {
    out.reason = pierce_failure_reason::degenerate_normal;
    return out;
  }
  n_out.normalize();

  const real t_min =
      std::max(cfg.t_min_scale * avg_edge_length, real(1e-8));
  const vec3 origin = x[static_cast<size_t>(foot)];

  index_t pierce = -1;
  const auto pierce_hit = ray_pierce_slab_shell(
      shell, foot, origin, n_out, x, face_tree, face_map, cfg, t_min,
      avg_edge_length);
  if (pierce_hit) {
    const asawa::shell::FaceId fi = pierce_hit->fi;
    const arp::ray_hit<3> &hit = pierce_hit->hit;
    out.hit_face = static_cast<index_t>(fi);
    out.ray_t = hit.t_ray;
    const auto tri = shell.get_tri(fi);
    pierce = closest_triangle_vertex_excluding(tri, x, hit.point, foot);
    if (pierce < 0) {
      out.reason = pierce_failure_reason::pierce_equals_foot;
      return out;
    }
  } else if (cfg.antipodal_fallback) {
    const real min_spine =
        std::max(cfg.self_skip_tol, cfg.t_min_scale * avg_edge_length);
    const auto ant = antipodal_pierce_vertex(shell, foot, x, n_out, min_spine,
                                             cfg.antipodal_min_dot);
    if (!ant) {
      out.reason = pierce_failure_reason::ray_miss;
      out.miss_detail = diagnose_slab_ray_pierce(shell, foot, x, face_tree,
                                                   face_map, cfg, avg_edge_length)
                            .miss_detail;
      return out;
    }
    pierce = *ant;
    out.used_antipodal_fallback = true;
    out.hit_face = -1;
    out.ray_t = (x[static_cast<size_t>(pierce)] - origin).norm();
  } else {
    out.reason = pierce_failure_reason::ray_miss;
    out.miss_detail = diagnose_slab_ray_pierce(shell, foot, x, face_tree,
                                                 face_map, cfg, avg_edge_length)
                          .miss_detail;
    return out;
  }

  if (pierce == out.stencil.foot) {
    out.reason = pierce_failure_reason::pierce_equals_foot;
    return out;
  }
  const real spine_len = (x[static_cast<size_t>(pierce)] - origin).norm();
  const real min_spine =
      std::max(cfg.self_skip_tol, cfg.t_min_scale * avg_edge_length);
  if (spine_len < min_spine) {
    out.reason = pierce_failure_reason::pierce_too_close;
    return out;
  }

  const asawa::shell::VertId pierce_vid =
      asawa::shell::vert_id(pierce);
  const auto bot_ring = shell.get_one_ring(pierce_vid);
  if (bot_ring.empty()) {
    out.reason = pierce_failure_reason::empty_bot_ring;
    return out;
  }

  out.stencil.bot_ids.push_back(pierce);
  for (asawa::shell::VertId v : bot_ring)
    out.stencil.bot_ids.push_back(static_cast<index_t>(v));

  const std::vector<real> bot_w =
      asawa::shell::vert_cotan_weights(shell, pierce_vid, x);
  if (bot_w.size() + 1 != out.stencil.bot_ids.size()) {
    out.reason = pierce_failure_reason::bot_cotan_mismatch;
    return out;
  }
  for (real w : bot_w)
    out.stencil.bot_weights.push_back(w);

  out.stencil.pierce_foot = pierce;
  out.stencil.hit_face = out.hit_face;
  out.stencil.has_slab = true;
  out.reason = pierce_failure_reason::ok;
  return out;
}

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // GAUDI_HEPWORTH_SLAB_PIERCE_STENCIL_HPP
