#ifndef GAUDI_ARP_DISK_TRIANGLE_HIT_HPP
#define GAUDI_ARP_DISK_TRIANGLE_HIT_HPP

#include "gaudi/common.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace gaudi {
namespace arp {

struct disk_triangle_hit_result {
  bool hit = false;
  real s = 0.0;
  real u = 0.0;
  real w = 0.0;
};

struct disk_triangle_best {
  real s = std::numeric_limits<real>::infinity();
  real u = 0.0;
  real w = 0.0;
};

inline bool disk_triangle_point_inside_2d(const vec2 &p0, const vec2 &p1,
                                          const vec2 &p2, const vec2 &p) {
  const vec2 v0 = p1 - p0;
  const vec2 v1 = p2 - p0;
  const vec2 v2 = p - p0;
  const real dot00 = v0.dot(v0);
  const real dot01 = v0.dot(v1);
  const real dot02 = v0.dot(v2);
  const real dot11 = v1.dot(v1);
  const real dot12 = v1.dot(v2);
  const real inv_denom = dot00 * dot11 - dot01 * dot01;
  if (std::abs(inv_denom) < 1e-24) {
    return false;
  }
  const real inv = real(1) / inv_denom;
  const real u = (dot11 * dot02 - dot01 * dot12) * inv;
  const real v = (dot00 * dot12 - dot01 * dot02) * inv;
  return u >= -1e-10 && v >= -1e-10 && u + v <= 1.0 + 1e-10;
}

inline real disk_triangle_bary_s_at_point(const vec2 &p0, const vec2 &p1,
                                          const vec2 &p2, real s0, real s1,
                                          real s2, const vec2 &p) {
  const vec2 v0 = p1 - p0;
  const vec2 v1 = p2 - p0;
  const vec2 v2 = p - p0;
  const real dot00 = v0.dot(v0);
  const real dot01 = v0.dot(v1);
  const real dot02 = v0.dot(v2);
  const real dot11 = v1.dot(v1);
  const real dot12 = v1.dot(v2);
  const real inv_denom = dot00 * dot11 - dot01 * dot01;
  if (std::abs(inv_denom) < 1e-24) {
    return s0;
  }
  const real inv = real(1) / inv_denom;
  const real u = (dot11 * dot02 - dot01 * dot12) * inv;
  const real v = (dot00 * dot12 - dot01 * dot02) * inv;
  const real w = real(1) - u - v;
  return w * s0 + u * s1 + v * s2;
}

inline void disk_triangle_consider_contact(real s_axial, real u, real w, real r2,
                                           real s_min, real s_max,
                                           disk_triangle_best &best) {
  const real perp2 = u * u + w * w;
  if (perp2 > r2 + 1e-12) {
    return;
  }
  const real h = std::sqrt(std::max(real(0), r2 - perp2));
  const real s_entry = s_axial - h;
  if (s_entry < s_min - 1e-10 || s_entry > s_max + 1e-10) {
    return;
  }
  if (s_entry < best.s) {
    best.s = s_entry;
    best.u = u;
    best.w = w;
  }
}

inline real disk_triangle_closest_alpha_segment(const vec2 &p0, const vec2 &p1) {
  const vec2 e = p1 - p0;
  const real len2 = e.squaredNorm();
  if (len2 < 1e-24) {
    return 0.0;
  }
  return std::clamp(-p0.dot(e) / len2, real(0), real(1));
}

inline void disk_triangle_consider_edge(real r2, real s0, real u0, real w0,
                                        real s1, real u1, real w1, real s_min,
                                        real s_max, disk_triangle_best &best) {
  const real ds = s1 - s0;
  const real du = u1 - u0;
  const real dw = w1 - w0;
  const vec2 p0(u0, w0);
  const vec2 p1(u1, w1);

  auto consider_alpha = [&](real alpha) {
    if (alpha < -1e-10 || alpha > 1.0 + 1e-10) {
      return;
    }
    const real a = std::clamp(alpha, real(0), real(1));
    const real u = u0 + a * du;
    const real w = w0 + a * dw;
    const real s = s0 + a * ds;
    disk_triangle_consider_contact(s, u, w, r2, s_min, s_max, best);
  };

  consider_alpha(0.0);
  consider_alpha(1.0);
  consider_alpha(disk_triangle_closest_alpha_segment(p0, p1));

  const real A = du * du + dw * dw;
  const real B = 2.0 * (u0 * du + w0 * dw);
  const real C = u0 * u0 + w0 * w0 - r2;
  if (std::abs(A) < 1e-14) {
    if (std::abs(B) > 1e-14) {
      consider_alpha(-C / B);
    }
    return;
  }
  const real disc = B * B - 4.0 * A * C;
  if (disc < 0.0) {
    return;
  }
  const real sq = std::sqrt(std::max(real(0), disc));
  consider_alpha((-B - sq) / (2.0 * A));
  consider_alpha((-B + sq) / (2.0 * A));
}

/// Disk at origin radius r in the uw plane; triangle verts (u,w) with along-ray s.
inline disk_triangle_hit_result
disk_triangle_hit_2d(real r, const vec2 &p0, const vec2 &p1, const vec2 &p2,
                     real s0, real s1, real s2, real s_min, real s_max) {
  disk_triangle_hit_result out;
  if (r <= 0.0) {
    return out;
  }

  const real r2 = r * r;
  disk_triangle_best best;

  if (disk_triangle_point_inside_2d(p0, p1, p2, vec2::Zero())) {
    const real s_origin =
        disk_triangle_bary_s_at_point(p0, p1, p2, s0, s1, s2, vec2::Zero());
    disk_triangle_consider_contact(s_origin, 0.0, 0.0, r2, s_min, s_max, best);
    if (std::isfinite(best.s)) {
      out.hit = true;
      out.s = best.s;
      out.u = best.u;
      out.w = best.w;
      return out;
    }
    best = disk_triangle_best{};
  }

  const real area2 =
      std::abs((p1.x() - p0.x()) * (p2.y() - p0.y()) -
               (p1.y() - p0.y()) * (p2.x() - p0.x()));
  if (area2 < 1e-20) {
    const real d01 = (p1 - p0).squaredNorm();
    const real d02 = (p2 - p0).squaredNorm();
    const real d12 = (p2 - p1).squaredNorm();
    if (d01 < 1e-20 && d02 < 1e-20) {
      disk_triangle_consider_contact(s0, p0.x(), p0.y(), r2, s_min, s_max, best);
    } else if (d01 >= d02 && d01 >= d12) {
      disk_triangle_consider_edge(r2, s1, p1.x(), p1.y(), s2, p2.x(), p2.y(),
                                  s_min, s_max, best);
    } else if (d02 >= d12) {
      disk_triangle_consider_edge(r2, s0, p0.x(), p0.y(), s2, p2.x(), p2.y(),
                                  s_min, s_max, best);
    } else {
      disk_triangle_consider_edge(r2, s0, p0.x(), p0.y(), s1, p1.x(), p1.y(),
                                  s_min, s_max, best);
    }
  } else {
    disk_triangle_consider_edge(r2, s0, p0.x(), p0.y(), s1, p1.x(), p1.y(),
                                s_min, s_max, best);
    disk_triangle_consider_edge(r2, s1, p1.x(), p1.y(), s2, p2.x(), p2.y(),
                                s_min, s_max, best);
    disk_triangle_consider_edge(r2, s2, p2.x(), p2.y(), s0, p0.x(), p0.y(),
                                s_min, s_max, best);
  }

  if (!std::isfinite(best.s)) {
    return out;
  }

  out.hit = true;
  out.s = best.s;
  out.u = best.u;
  out.w = best.w;
  return out;
}

} // namespace arp
} // namespace gaudi

#endif // GAUDI_ARP_DISK_TRIANGLE_HIT_HPP
