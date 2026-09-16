#ifndef GAUDI_ARP_RAY_SIMPLEX_HPP
#define GAUDI_ARP_RAY_SIMPLEX_HPP

#include "gaudi/arp/disk_triangle_hit.hpp"
#include "gaudi/arp/pairwise_tests.hpp"
#include "gaudi/common.h"
#include "gaudi/vec_addendum.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>

namespace gaudi {
namespace arp {

struct Ray {
  vec3 origin = vec3::Zero();
  vec3 dir = vec3::UnitZ();
  real t_min = 0.0;
  real t_max = std::numeric_limits<real>::infinity();

  Ray() = default;
  Ray(vec3 origin_in, vec3 dir_in, real t_min_in = 0.0,
      real t_max_in = std::numeric_limits<real>::infinity())
      : origin(std::move(origin_in)), dir(std::move(dir_in)), t_min(t_min_in),
        t_max(t_max_in) {
    normalize_dir();
  }

  vec3 at(real t) const { return origin + t * dir; }

  void normalize_dir() {
    const real n = dir.norm();
    if (n > 1e-12)
      dir /= n;
  }
};

struct point_local {};
struct line_local {
  real t_seg = 0.0;
};
struct tri_local {
  real u = 0.0;
  real v = 0.0;
  real w() const { return 1.0 - u - v; }
};

template <int N> struct simplex_local;

template <> struct simplex_local<1> {
  using type = point_local;
};
template <> struct simplex_local<2> {
  using type = line_local;
};
template <> struct simplex_local<3> {
  using type = tri_local;
};

template <int Ndata> struct ray_hit {
  index_t id = -1;
  real t_ray = 0.0;
  typename simplex_local<Ndata>::type local{};
  vec3 point = vec3::Zero();
  bool valid = false;
};

inline real ray_distance_to_point(const Ray &ray, const vec3 &p) {
  const vec3 w = p - ray.origin;
  const vec3 perp = w - ray.dir * w.dot(ray.dir);
  return perp.norm();
}

inline real ray_param_on_infinite_line(const Ray &ray, const vec3 &p) {
  return (p - ray.origin).dot(ray.dir);
}

inline std::optional<ray_hit<1>>
ray_point_within_tol(const Ray &ray, const vec3 &p, real tol,
                     index_t id = -1) {
  const real d = ray_distance_to_point(ray, p);
  if (d > tol)
    return std::nullopt;
  ray_hit<1> out;
  out.id = id;
  out.t_ray = ray_param_on_infinite_line(ray, p);
  out.point = p;
  out.valid = true;
  return out;
}

inline std::optional<ray_hit<2>>
ray_pierce_segment(const Ray &ray, const vec3 &p0, const vec3 &p1,
                   index_t id = -1, real eps = 1e-8) {
  const vec3 edge = p1 - p0;
  const vec3 h = ray.dir.cross(edge);
  const real det = h.squaredNorm();
  if (det < eps * eps)
    return std::nullopt;

  const vec3 s = ray.origin - p0;
  const real u = s.cross(edge).dot(h) / det;
  if (u < ray.t_min - eps || u > ray.t_max + eps)
    return std::nullopt;

  const real v = s.cross(ray.dir).dot(h) / det;
  if (v < -eps || v > 1.0 + eps)
    return std::nullopt;

  ray_hit<2> out;
  out.id = id;
  out.t_ray = u;
  out.local.t_seg = v;
  out.point = ray.at(u);
  out.valid = true;
  return out;
}

inline std::optional<ray_hit<3>>
ray_pierce_triangle(const Ray &ray, const vec3 &v0, const vec3 &v1,
                    const vec3 &v2, index_t id = -1, real eps = 1e-8) {
  const vec3 e1 = v1 - v0;
  const vec3 e2 = v2 - v0;
  const vec3 pvec = ray.dir.cross(e2);
  const real det = e1.dot(pvec);
  if (std::abs(det) < eps)
    return std::nullopt;

  const vec3 tvec = ray.origin - v0;
  const real u = tvec.dot(pvec) / det;
  if (u < -eps || u > 1.0 + eps)
    return std::nullopt;

  const vec3 qvec = tvec.cross(e1);
  const real v = ray.dir.dot(qvec) / det;
  if (v < -eps || u + v > 1.0 + eps)
    return std::nullopt;

  const real t = e2.dot(qvec) / det;
  if (t < ray.t_min - eps || t > ray.t_max + eps)
    return std::nullopt;

  ray_hit<3> out;
  out.id = id;
  out.t_ray = t;
  out.local.u = u;
  out.local.v = v;
  out.point = ray.at(t);
  out.valid = true;
  return out;
}

inline void ray_null_axes(const vec3 &dir, vec3 &e1, vec3 &e2) {
  const vec3 d = dir.normalized();
  const real ax = std::abs(d.x());
  const real ay = std::abs(d.y());
  const real az = std::abs(d.z());
  vec3 a = vec3(1, 0, 0);
  if (ay <= ax && ay <= az) {
    a = vec3(0, 1, 0);
  } else if (az <= ax && az <= ay) {
    a = vec3(0, 0, 1);
  }
  e1 = d.cross(a).normalized();
  e2 = d.cross(e1).normalized();
}

struct ray_slop_coords {
  real s = 0.0;
  real u = 0.0;
  real w = 0.0;
  real perp2() const { return u * u + w * w; }
};

inline ray_slop_coords ray_slop_coords_of(const Ray &ray, const vec3 &e1,
                                          const vec3 &e2, const vec3 &p) {
  const vec3 rel = p - ray.origin;
  return {rel.dot(ray.dir), rel.dot(e1), rel.dot(e2)};
}

/// Ray + slop (cylinder along ray axis): earliest t where the tube contacts triangle.
inline std::optional<ray_hit<3>>
ray_slop_intersects_triangle(const Ray &ray, real slop, const vec3 &v0,
                             const vec3 &v1, const vec3 &v2, index_t id = -1) {
  if (slop <= 0.0) {
    return ray_pierce_triangle(ray, v0, v1, v2, id);
  }

  vec3 e1, e2;
  ray_null_axes(ray.dir, e1, e2);
  const ray_slop_coords c0 = ray_slop_coords_of(ray, e1, e2, v0);
  const ray_slop_coords c1 = ray_slop_coords_of(ray, e1, e2, v1);
  const ray_slop_coords c2 = ray_slop_coords_of(ray, e1, e2, v2);

  const disk_triangle_hit_result dt = disk_triangle_hit_2d(
      slop, vec2(c0.u, c0.w), vec2(c1.u, c1.w), vec2(c2.u, c2.w), c0.s, c1.s,
      c2.s, ray.t_min, ray.t_max);
  if (!dt.hit) {
    return std::nullopt;
  }

  const real perp2 = dt.u * dt.u + dt.w * dt.w;
  const real h =
      std::sqrt(std::max(real(0), slop * slop - perp2));
  const vec3 contact = ray.origin + (dt.s + h) * ray.dir + dt.u * e1 + dt.w * e2;

  const vec3 e1_tri = v1 - v0;
  const vec3 e2_tri = v2 - v0;
  const vec3 tvec = contact - v0;
  const vec3 pvec = ray.dir.cross(e2_tri);
  const real det = e1_tri.dot(pvec);

  ray_hit<3> out;
  out.id = id;
  out.t_ray = dt.s;
  out.point = contact;
  if (std::abs(det) > 1e-12) {
    out.local.u = tvec.dot(pvec) / det;
    const vec3 qvec = tvec.cross(e1_tri);
    out.local.v = ray.dir.dot(qvec) / det;
  }
  out.valid = true;
  return out;
}

inline real ray_distance_to_segment(const Ray &ray, const vec3 &p0,
                                    const vec3 &p1) {
  const vec3 far = ray.origin + ray.dir * std::max(ray.t_max, real(1.0));
  const std::array<real, 3> d =
      va::distance_Segment_Segment(ray.origin, far, p0, p1);
  return d[0];
}

inline real ray_distance_to_triangle(const Ray &ray, const vec3 &v0,
                                     const vec3 &v1, const vec3 &v2) {
  vec3 cp;
  const std::array<vec3, 3> tri = {v0, v1, v2};
  const real d0 = va::distance_from_triangle(tri, ray.origin, cp);
  const real t_proj = ray_param_on_infinite_line(ray, cp);
  const vec3 on_ray = ray.at(t_proj);
  vec3 cp2;
  const real d1 = va::distance_from_triangle(tri, on_ray, cp2);
  return std::min(d0, d1);
}

template <int Ndata>
inline std::optional<ray_hit<Ndata>>
ray_pierce_simplex(const Ray &ray, const std::array<vec3, Ndata> &simplex,
                   index_t id = -1) {
  if constexpr (Ndata == 1) {
    if (ray_distance_to_point(ray, simplex[0]) > 1e-8)
      return std::nullopt;
    ray_hit<1> out;
    out.id = id;
    out.t_ray = ray_param_on_infinite_line(ray, simplex[0]);
    if (out.t_ray < ray.t_min || out.t_ray > ray.t_max)
      return std::nullopt;
    out.point = simplex[0];
    out.valid = true;
    return out;
  } else if constexpr (Ndata == 2) {
    return ray_pierce_segment(ray, simplex[0], simplex[1], id);
  } else {
    return ray_pierce_triangle(ray, simplex[0], simplex[1], simplex[2], id);
  }
}

template <int Ndata>
inline std::optional<ray_hit<Ndata>>
ray_slop_intersects_simplex(const Ray &ray, real slop,
                            const std::array<vec3, Ndata> &simplex,
                            index_t id = -1) {
  if constexpr (Ndata == 3) {
    return ray_slop_intersects_triangle(ray, slop, simplex[0], simplex[1],
                                        simplex[2], id);
  } else if constexpr (Ndata == 2) {
    if (slop <= 0.0) {
      return ray_pierce_segment(ray, simplex[0], simplex[1], id);
    }
    vec3 e1, e2;
    ray_null_axes(ray.dir, e1, e2);
    const auto c0 = ray_slop_coords_of(ray, e1, e2, simplex[0]);
    const auto c1 = ray_slop_coords_of(ray, e1, e2, simplex[1]);
    const disk_triangle_hit_result dt = disk_triangle_hit_2d(
        slop, vec2(c0.u, c0.w), vec2(c1.u, c1.w), vec2(c1.u, c1.w), c0.s, c1.s,
        c1.s, ray.t_min, ray.t_max);
    if (!dt.hit) {
      return std::nullopt;
    }
    const real perp2 = dt.u * dt.u + dt.w * dt.w;
    const real h =
        std::sqrt(std::max(real(0), slop * slop - perp2));
    ray_hit<2> out;
    out.id = id;
    out.t_ray = dt.s;
    out.point = ray.origin + (dt.s + h) * ray.dir + dt.u * e1 + dt.w * e2;
    out.valid = true;
    return out;
  } else {
    if (slop <= 0.0) {
      return ray_pierce_simplex<1>(ray, simplex, id);
    }
    return ray_point_within_tol(ray, simplex[0], slop, id);
  }
}

template <int Ndata>
inline real ray_distance_to_simplex(const Ray &ray,
                                    const std::array<vec3, Ndata> &simplex) {
  if constexpr (Ndata == 1) {
    return ray_distance_to_point(ray, simplex[0]);
  } else if constexpr (Ndata == 2) {
    return ray_distance_to_segment(ray, simplex[0], simplex[1]);
  } else {
    return ray_distance_to_triangle(ray, simplex[0], simplex[1], simplex[2]);
  }
}

template <int Ndata>
inline bool ray_collect_simplex_within_tol(
    const Ray &ray, const std::array<vec3, Ndata> &simplex, real tol,
    index_t id, std::vector<ray_hit<Ndata>> &out) {
  const real d = ray_distance_to_simplex<Ndata>(ray, simplex);
  if (d > tol)
    return false;
  ray_hit<Ndata> hit;
  hit.id = id;
  hit.t_ray = ray_param_on_infinite_line(ray, simplex[0]);
  hit.point = simplex[0];
  hit.valid = true;
  out.push_back(hit);
  return true;
}

template <int Ndata>
inline void ray_hit_keep_nearest(std::optional<ray_hit<Ndata>> &best,
                                 const std::optional<ray_hit<Ndata>> &candidate) {
  if (!candidate)
    return;
  if (!best || candidate->t_ray < best->t_ray)
    best = candidate;
}

} // namespace arp
} // namespace gaudi

#endif // GAUDI_ARP_RAY_SIMPLEX_HPP
