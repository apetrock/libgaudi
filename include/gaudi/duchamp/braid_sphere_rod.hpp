#ifndef GAUDI_DUCHAMP_BRAID_SPHERE_ROD_HPP
#define GAUDI_DUCHAMP_BRAID_SPHERE_ROD_HPP

#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/common.h"
#include "gaudi/duchamp/braid_planar_rod.hpp"
#include "gaudi/windychien/braid.hpp"
#include "gaudi/windychien/catalog.hpp"

#include <cmath>
#include <stdexcept>
#include <vector>

namespace gaudi {
namespace duchamp {

struct braid_sphere_params {
  real radius = 1.0;
  real lat_min = -0.45; // elevation from equator (radians)
  real lat_max = 0.45;
  /// r = radius + radial_gain * planar_z
  real radial_gain = 1.0;
  /// Circumferential subdivisions: circle-space dx = circumference / n_lat.
  /// Braid-space target edge before project: (N_crossing * dx) / n_lat.
  int n_lat = 36;
  /// Extra densification vs n_lat (braid max_len /= subdiv_safety).
  real subdiv_safety = 1.0;
  int max_subdiv_iters = 32;
};

namespace sphere_detail {

inline real wrap_x_toward(real x, real toward, real period) {
  while (x - toward > 0.5 * period) {
    x -= period;
  }
  while (toward - x > 0.5 * period) {
    x += period;
  }
  return x;
}

inline real mod_period(real x, real period) {
  real y = std::fmod(x, period);
  if (y < 0.0) {
    y += period;
  }
  return y;
}

inline vec3 midpoint_periodic(const vec3 &a, const vec3 &b, real period) {
  const real bx = wrap_x_toward(b[0], a[0], period);
  vec3 m = 0.5 * (a + vec3(bx, b[1], b[2]));
  m[0] = mod_period(m[0], period);
  return m;
}

inline real edge_length(const vec3 &a, const vec3 &b, real period) {
  if (period > 0.0) {
    const real bx = wrap_x_toward(b[0], a[0], period);
    return (vec3(bx, b[1], b[2]) - a).norm();
  }
  return (b - a).norm();
}

inline asawa::rod::CornerId insert_vert(asawa::rod::rod &R, const vec3 &p) {
  const asawa::rod::CornerId id = R.insert_edge();
  R.corner_verts().push_back(p);
  return id;
}

} // namespace sphere_detail

/// Split edges longer than max_len. period > 0 unwraps x on wrap edges.
inline void subdivide_rod_to_length(asawa::rod::rod &R, real max_len,
                                    real period = 0.0, int max_iters = 24) {
  if (max_len <= 0.0) {
    return;
  }
  for (int iter = 0; iter < max_iters; ++iter) {
    std::vector<asawa::rod::CornerId> to_split;
    for (size_t i = 0; i < R.corner_count(); ++i) {
      const auto ci = asawa::rod::corner_id(static_cast<int>(i));
      const auto cj = R.next(ci);
      if (cj < asawa::rod::corner_id(0)) {
        continue;
      }
      const real len = sphere_detail::edge_length(
          R.x()[i], R.x()[static_cast<size_t>(cj)], period);
      if (len > max_len) {
        to_split.push_back(ci);
      }
    }
    if (to_split.empty()) {
      break;
    }
    for (asawa::rod::CornerId c0 : to_split) {
      const asawa::rod::CornerId c1 = R.next(c0);
      if (c1 < asawa::rod::corner_id(0)) {
        continue;
      }
      const vec3 &a = R.x()[static_cast<size_t>(c0)];
      const vec3 &b = R.x()[static_cast<size_t>(c1)];
      const vec3 mid = (period > 0.0)
                           ? sphere_detail::midpoint_periodic(a, b, period)
                           : vec3(0.5 * (a + b));
      const asawa::rod::CornerId cnew = sphere_detail::insert_vert(R, mid);
      R.link(c0, cnew);
      R.link(cnew, c1);
    }
  }
  R._init_params();
}

/// Chart: x → lon [0,2π), y → lat [lat_min, lat_max], z → radial offset.
inline void project_rod_to_sphere(asawa::rod::rod &R,
                                  const braid_sphere_params &sp, real x_period,
                                  real y_min, real y_max) {
  const real dy = std::max(y_max - y_min, real(1e-12));
  const real dlat = sp.lat_max - sp.lat_min;
  for (vec3 &q : R.x()) {
    // Keep longitude continuous for wrap beads at x ∈ (0, x_period].
    real x = q[0];
    if (x_period > 0.0) {
      x = std::fmod(x, x_period);
      if (x < 0.0) {
        x += x_period;
      }
    }
    const real lon =
        (x_period > 0.0) ? (real(2) * real(M_PI) * (x / x_period)) : 0.0;
    const real lat = sp.lat_min + dlat * ((q[1] - y_min) / dy);
    const real r = sp.radius + sp.radial_gain * q[2];
    const real cl = std::cos(lat);
    q = vec3(r * cl * std::cos(lon), r * cl * std::sin(lon), r * std::sin(lat));
  }
  R._init_params();
}

/// Closed braid on a sphere: build in braid space → subdivide → project.
/// Target chart edge: (N_crossing * dx) / n_lat / subdiv_safety.
inline asawa::rod::rod::ptr
braid_to_sphere_rod(const windychien::braid &b,
                    braid_planar_params planar = {},
                    const braid_sphere_params &sphere = {}) {
  windychien::validate_braid(b);
  const int n = b.strands;
  const int L = windychien::braid_n_columns(b); // axial frames / columns
  if (L < 1) {
    throw std::runtime_error("braid_to_sphere_rod: empty word");
  }
  if (sphere.n_lat < 1) {
    throw std::runtime_error("braid_to_sphere_rod: n_lat must be >= 1");
  }

  // Stay in braid/lanyard chart units (do not replace dx with sphere arcs).
  planar.close = true;
  planar.center = false;
  if (planar.dx <= 0.0) {
    planar.dx = 1.0;
  }
  if (planar.dy <= 0.0) {
    planar.dy = planar.dx;
  }

  asawa::rod::rod::ptr R = braid_to_planar_rod(b, planar);

  const real x_period = planar.dx * real(L); // N_crossing * dx
  const real y_min = 0.0;
  const real y_max = planar.dy * real(std::max(n - 1, 0));
  const real safety = std::max(sphere.subdiv_safety, real(1e-6));
  // Circle-space dx = circumference / n_lat; same ratio in braid space:
  // max_len = period / n_lat = (N_crossing * dx) / n_lat.
  const real max_len =
      std::max(x_period / real(sphere.n_lat) / safety, real(1e-6));

  subdivide_rod_to_length(*R, max_len, x_period, sphere.max_subdiv_iters);
  project_rod_to_sphere(*R, sphere, x_period, y_min, y_max);
  return R;
}

inline asawa::rod::rod::ptr
load_braid_sphere_rod(const std::string &knot_name,
                      const braid_planar_params &planar = {},
                      const braid_sphere_params &sphere = {}) {
  return braid_to_sphere_rod(windychien::load_braid(knot_name), planar, sphere);
}

} // namespace duchamp
} // namespace gaudi

#endif
