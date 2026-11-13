#ifndef GAUDI_HEPWORTH_BLOCK_SURFACE_WALK_HPP
#define GAUDI_HEPWORTH_BLOCK_SURFACE_WALK_HPP

/// Straightest walks on a triangle mesh (scripts/sympy/cr_poisson/
/// curvature_advection.py `walk`): from a point in a face, along a tangent
/// direction, crossing edges into the neighbouring face with the direction
/// unfolded about the shared edge (its component along the edge kept, the
/// rest turned into the new face's plane). Samples are taken every `step`
/// of arc length. The walk stops at a boundary.

#include "gaudi/kusama/stencil_table.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace gaudi {
namespace hepworth {
namespace block {

struct walk_sample {
  int face = -1;
  vec3 bary = vec3::Zero(); ///< barycentric in the face's vertex order
  vec3 p = vec3::Zero();
  vec3 dir = vec3::Zero(); ///< unit, in the face's plane
  real s = 0;              ///< arc length from the start
};

namespace walk_detail {

inline vec3 face_normal(const kusama::stencil_table &t, const std::vector<vec3> &x, int f) {
  const auto &id = t.ids[static_cast<size_t>(f)];
  const vec3 n = (x[static_cast<size_t>(id[1])] - x[static_cast<size_t>(id[0])])
                     .cross(x[static_cast<size_t>(id[2])] - x[static_cast<size_t>(id[0])]);
  const real l = n.norm();
  return l > real(0) ? vec3(n / l) : vec3::Zero();
}

/// Barycentric coordinates of `p` (projected into the face's plane).
inline vec3 bary(const kusama::stencil_table &t, const std::vector<vec3> &x, int f, const vec3 &p) {
  const auto &id = t.ids[static_cast<size_t>(f)];
  const vec3 a = x[static_cast<size_t>(id[0])], e1 = x[static_cast<size_t>(id[1])] - a,
             e2 = x[static_cast<size_t>(id[2])] - a, d = p - a;
  const real g11 = e1.dot(e1), g12 = e1.dot(e2), g22 = e2.dot(e2), det = g11 * g22 - g12 * g12;
  if (!(std::abs(det) > real(0)))
    return vec3(1, 0, 0);
  const real b1 = (g22 * e1.dot(d) - g12 * e2.dot(d)) / det, b2 = (g11 * e2.dot(d) - g12 * e1.dot(d)) / det;
  return vec3(1 - b1 - b2, b1, b2);
}

/// The rate of change of the barycentric coordinates along a tangent `d`.
inline vec3 bary_dir(const kusama::stencil_table &t, const std::vector<vec3> &x, int f, const vec3 &d) {
  const auto &id = t.ids[static_cast<size_t>(f)];
  const vec3 a = x[static_cast<size_t>(id[0])];
  const vec3 b = bary(t, x, f, a + d);
  return b - vec3(1, 0, 0);
}

} // namespace walk_detail

/// Walk `length` from `p0` in face `f0` along `d0`, sampling every `step`
/// (the start included). Returns the samples in order.
inline std::vector<walk_sample> trace_on_surface(const kusama::stencil_table &t, const std::vector<vec3> &x, int f0,
                                                 const vec3 &p0, const vec3 &d0, real length, real step) {
  using namespace walk_detail;
  std::vector<walk_sample> out;
  if (f0 < 0 || f0 >= t.n_faces || !(step > real(0)))
    return out;
  int f = f0;
  vec3 n = face_normal(t, x, f);
  vec3 d = d0 - n * n.dot(d0);
  if (!(d.norm() > real(0)))
    return out;
  d.normalize();
  vec3 p = p0;
  real s = 0, next = 0;
  int guard = 0;
  while (guard++ < 100000) {
    const vec3 lam = bary(t, x, f, p);
    const vec3 dl = bary_dir(t, x, f, d);
    // Arc length to leave the face (a coordinate reaching zero while decreasing).
    real tau = std::numeric_limits<real>::infinity();
    int k_exit = -1;
    const real scale = dl.cwiseAbs().maxCoeff();
    for (int k = 0; k < 3; ++k) {
      if (!(dl[k] < -real(1e-9) * scale))
        continue;
      const real tk = -lam[k] / dl[k];
      if (tk < tau) {
        tau = std::max(tk, real(0));
        k_exit = k;
      }
    }
    // Samples inside this face.
    while (next <= length + real(1e-12) && next <= s + tau + real(1e-12)) {
      walk_sample w;
      w.face = f;
      w.p = p + (next - s) * d;
      w.bary = bary(t, x, f, w.p);
      w.dir = d;
      w.s = next;
      out.push_back(w);
      next += step;
    }
    if (s + tau >= length || k_exit < 0)
      break;
    // Cross the edge opposite vertex `k_exit`: local edge (k + 1) % 3.
    p += tau * d;
    s += tau;
    const int le = (k_exit + 1) % 3;
    const int e = t.face_edges[static_cast<size_t>(f)][static_cast<size_t>(le)];
    if (e < 0)
      break;
    const auto &ef = t.edge_faces[static_cast<size_t>(e)];
    const int g = ef[0] == f ? ef[1] : ef[0];
    if (g < 0)
      break;
    const auto &id = t.ids[static_cast<size_t>(f)];
    vec3 ea = x[static_cast<size_t>(id[static_cast<size_t>(le)])],
         eb = x[static_cast<size_t>(id[static_cast<size_t>((le + 1) % 3)])];
    vec3 ax = (eb - ea).normalized();
    const vec3 n2 = face_normal(t, x, g);
    // Unfold: keep the component along the edge, turn the rest into the new plane.
    const real along = d.dot(ax);
    vec3 perp = n.cross(ax);
    const real across = d.dot(perp);
    vec3 perp2 = n2.cross(ax);
    d = along * ax + across * perp2;
    d -= n2 * n2.dot(d);
    if (!(d.norm() > real(0)))
      break;
    d.normalize();
    // Nudge into the new face so the next exit search does not stall on the edge.
    p += real(1e-12) * d;
    f = g;
    n = n2;
  }
  return out;
}

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif
