#ifndef GAUDI_HEPWORTH_BLOCK_MOBIUS_GROWTH_HPP
#define GAUDI_HEPWORTH_BLOCK_MOBIUS_GROWTH_HPP

/// Moebius growth (scripts/sympy/cr_poisson/MATH.md 19.20): a Moebius map
/// is locally a conformal scale `e^u`; on a mesh each edge scales by its
/// endpoints' factors, `l_ij -> e^{(u_i + u_j) / 2} l_ij`, so a triangle
/// whose vertices carry different `u` changes shape (the kite at the triangle
/// level). The per-frame vertex log-scale `u` comes from the frame-start
/// charts' curvature tensors (`vertex_chart::h`, `S = -dN`, outward normals,
/// so convex `H < 0`), by `source`:
/// - `mean_curvature`: `u = alpha (-H - level)`: convex grows, concave
///   shrinks, saddles stay;
/// - `anisotropy`: `u = alpha (|z| - level)`, `z` the traceless part;
/// - `flat`: `u = alpha (1 - |k_lo| / |k_hi| - level)`: grows where one
///   principal curvature goes to zero (around parabolic lines), not at round
///   points.
/// With `balanced` the area-weighted mean is removed (redistributes scale, no
/// net growth). Per edge an `edge_strain` row with goal length
/// `e^{(u_i + u_j) / 2}` times its frame-start length, rebuilt each frame, so
/// the growth accumulates plastically. Positions only.

#include "gaudi/hepworth/block/cr_shape.hpp"
#include "gaudi/hepworth/block/shell_constraints.hpp"
#include "gaudi/kusama/stencil_table.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace gaudi {
namespace hepworth {
namespace block {

enum class mobius_source { mean_curvature, anisotropy, flat };

/// `alpha`: log-scale per unit of the source per frame; `max_u`: cap on the
/// per-frame `|u|`; `level`: the source value at which nothing changes;
/// `balanced`: remove the mean.
struct mobius_target {
  real alpha = 0.02;
  real max_u = 0.1;
  mobius_source source = mobius_source::mean_curvature;
  real level = 0;
  bool balanced = false;
};

namespace growth_detail {

/// Mean curvature, traceless magnitude, and the principal curvatures ordered
/// by magnitude (`k_hi` the tight one) of a chart tensor.
struct principals {
  real H = 0, z = 0, k_hi = 0, k_lo = 0;
  /// Unit principal directions in the chart's `(T, B)`, tight then loose,
  /// right-handed (`e_lo` is `e_hi` turned a quarter toward `B`).
  vec2 e_hi = vec2(1, 0), e_lo = vec2(0, 1);
};

inline principals principals_of(const vec3 &h) {
  principals p;
  p.H = real(0.5) * (h[0] + h[2]);
  const real dx = real(0.5) * (h[0] - h[2]);
  p.z = std::hypot(dx, h[1]);
  const real a = p.H + p.z, b = p.H - p.z;
  const real th = real(0.5) * std::atan2(h[1], dx);
  const vec2 ea(std::cos(th), std::sin(th)), eb(-std::sin(th), std::cos(th));
  if (std::abs(a) >= std::abs(b)) {
    p.k_hi = a;
    p.k_lo = b;
    p.e_hi = ea;
    p.e_lo = eb;
  } else {
    p.k_hi = b;
    p.k_lo = a;
    p.e_hi = eb;
    p.e_lo = -ea;
  }
  return p;
}

inline std::array<int, 2> edge_ends(const kusama::stencil_table &t, int e) {
  for (int s = 0; s < 2; ++s) {
    const int f = t.edge_faces[static_cast<size_t>(e)][static_cast<size_t>(s)];
    const int i = t.edge_local[static_cast<size_t>(e)][static_cast<size_t>(s)];
    if (f < 0 || i < 0)
      continue;
    const auto &id = t.ids[static_cast<size_t>(f)];
    return {id[static_cast<size_t>(i)], id[static_cast<size_t>((i + 1) % 3)]};
  }
  return {-1, -1};
}

inline std::vector<real> vertex_areas(const kusama::stencil_table &t, const std::vector<vec3> &x) {
  std::vector<real> a(static_cast<size_t>(t.n_verts), real(0));
  for (int f = 0; f < t.n_faces; ++f) {
    const auto &id = t.ids[static_cast<size_t>(f)];
    const real A = real(0.5) * (x[static_cast<size_t>(id[1])] - x[static_cast<size_t>(id[0])])
                                   .cross(x[static_cast<size_t>(id[2])] - x[static_cast<size_t>(id[0])])
                                   .norm();
    for (int k = 0; k < 3; ++k)
      a[static_cast<size_t>(id[static_cast<size_t>(k)])] += A / real(3);
  }
  return a;
}

} // namespace growth_detail

/// The per-frame log-scale from the charts, capped at `max_u` (0 where a
/// chart is missing).
inline std::vector<real> mobius_log_scale(const kusama::stencil_table &t, const std::vector<vec3> &x,
                                          const std::vector<vertex_chart> &charts, const mobius_target &spec) {
  const size_t n = static_cast<size_t>(std::max(t.n_verts, 0));
  std::vector<real> u(n, real(0));
  if (charts.size() != n || x.size() != n)
    return u;
  for (size_t i = 0; i < n; ++i) {
    if (!charts[i].ok)
      continue;
    const growth_detail::principals p = growth_detail::principals_of(charts[i].h);
    real s = 0;
    if (spec.source == mobius_source::mean_curvature)
      s = -p.H;
    else if (spec.source == mobius_source::anisotropy)
      s = p.z;
    else
      s = std::abs(p.k_hi) > real(0) ? real(1) - std::abs(p.k_lo) / std::abs(p.k_hi) : real(0);
    u[i] = spec.alpha * (s - spec.level);
  }
  if (spec.balanced) {
    const std::vector<real> a = growth_detail::vertex_areas(t, x);
    real su = 0, sa = 0;
    for (size_t i = 0; i < n; ++i) {
      su += a[i] * u[i];
      sa += a[i];
    }
    const real mean = su / std::max(sa, real(1e-300));
    for (real &v : u)
      v -= mean;
  }
  if (spec.max_u > real(0))
    for (real &v : u)
      v = std::clamp(v, -spec.max_u, spec.max_u);
  return u;
}

/// One `edge_strain` row per edge with goal length `e^{(u_i + u_j) / 2}`
/// times its frame-start length. Returns the log-scale used (empty when
/// nothing was added).
inline std::vector<real> init_mobius_growth(const kusama::stencil_table &t, const std::vector<vec3> &x,
                                            const std::vector<vertex_chart> &charts, const mobius_target &spec,
                                            real w, std::vector<projection_constraint::ptr> &constraints,
                                            std::vector<sim_block::ptr> blocks) {
  if (!(w > real(0)) || static_cast<int>(x.size()) != t.n_verts)
    return {};
  const std::vector<real> u = mobius_log_scale(t, x, charts, spec);
  for (int e = 0; e < t.n_edges; ++e) {
    const std::array<int, 2> v = growth_detail::edge_ends(t, e);
    if (v[0] < 0)
      continue;
    const real l0 = (x[static_cast<size_t>(v[1])] - x[static_cast<size_t>(v[0])]).norm();
    if (!(l0 > real(0)))
      continue;
    auto row = edge_strain::create(std::vector<index_t>{v[0], v[1]}, x, w, blocks);
    row->set_goal_length(l0 * std::exp(real(0.5) * (u[static_cast<size_t>(v[0])] + u[static_cast<size_t>(v[1])])));
    constraints.push_back(row);
  }
  return u;
}

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif
