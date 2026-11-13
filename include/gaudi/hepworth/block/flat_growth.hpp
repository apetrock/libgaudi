#ifndef GAUDI_HEPWORTH_BLOCK_FLAT_GROWTH_HPP
#define GAUDI_HEPWORTH_BLOCK_FLAT_GROWTH_HPP

/// Directional growth with a curvature policy. Per edge, goal length
/// `l0 (1 + g w_e ((1 - iso) (dhat0 . a_e)^2 + iso))`:
/// - `a_e` the growth axis of its ends, `cos(angle) e_lo + sin(angle) e_hi`
///   from the frame-start charts (sign-free; `angle` 0 the loose principal
///   direction, `pi / 2` the tight one);
/// - `w_e` the mean of its ends' weights `m_v^power` for `flatness` (already
///   in [0, 1]) and `(m_v / mean m)^power` for the dimensional measures, `m`
///   the chosen `measure` of the chart tensor (`power` 0: uniform);
/// - `iso` mixes in growth independent of direction (`iso = 1`,
///   `power = 0`: a plain uniform scale growth).
/// With `flatness` and a high power the growth sits where one principal
/// curvature vanishes: on a torus the top and bottom of the tube, where the
/// loose direction is the hoop, so those circles lengthen and ruffle.
/// Rebuilt each frame against the frame-start mesh and charts; positions
/// only.

#include "gaudi/hepworth/block/mobius_growth.hpp"
#include "gaudi/hepworth/block/shape_objectives.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace gaudi {
namespace hepworth {
namespace block {

/// What sets where the growth is:
/// - `flatness`: `1 - |k_lo| / |k_hi|` (1 where one curvature vanishes, 0 at
///   round points);
/// - `mean_curvature`: `|H|` (negative power favours flatter regions);
/// - `anisotropy`: `|z| = |k_hi - k_lo| / 2`;
/// - `k_lo`: `|k_lo|` (negative power favours a flat loose direction).
enum class growth_measure { flatness, mean_curvature, anisotropy, k_lo };

/// `g`: fractional length added per frame at full weight; `measure` and
/// `power`: the weight (`flatness^power`, or `(m / mean m)^power`); `angle`: the axis, radians from
/// the loose principal direction toward the tight one; `winding`: sets
/// `angle` at startup from a closed line on the fitted torus (resolved by the
/// caller, as the shear's); `iso`: the share of direction-free growth.
struct flat_growth_target {
  real g = 0.05;
  growth_measure measure = growth_measure::flatness;
  real power = 8;
  real angle = 0;
  stretch_winding winding{};
  real iso = 0;
};

/// The measure of one chart tensor.
inline real growth_measure_of(const growth_detail::principals &p, growth_measure m) {
  switch (m) {
  case growth_measure::mean_curvature:
    return std::abs(p.H);
  case growth_measure::anisotropy:
    return p.z;
  case growth_measure::k_lo:
    return std::abs(p.k_lo);
  case growth_measure::flatness:
  default:
    return std::abs(p.k_hi) > real(0) ? real(1) - std::abs(p.k_lo) / std::abs(p.k_hi) : real(0);
  }
}

/// Per vertex: the weight and the growth axis in world (zero where the chart
/// is missing).
inline void flat_growth_field(const std::vector<vertex_chart> &charts, const flat_growth_target &spec,
                              std::vector<real> &weight, std::vector<vec3> &dir) {
  const size_t n = charts.size();
  weight.assign(n, real(0));
  dir.assign(n, vec3::Zero());
  std::vector<real> m(n, real(0));
  real mbar = 0;
  int nm = 0;
  for (size_t i = 0; i < n; ++i) {
    if (!charts[i].ok)
      continue;
    const growth_detail::principals p = growth_detail::principals_of(charts[i].h);
    m[i] = growth_measure_of(p, spec.measure);
    mbar += m[i];
    ++nm;
    const vec2 e = std::cos(spec.angle) * p.e_lo + std::sin(spec.angle) * p.e_hi;
    dir[i] = (e[0] * charts[i].T + e[1] * charts[i].B).normalized();
  }
  mbar = nm > 0 ? mbar / nm : real(0);
  for (size_t i = 0; i < n; ++i) {
    if (!charts[i].ok)
      continue;
    real w = 1;
    if (spec.power != real(0)) {
      const real scale = spec.measure == growth_measure::flatness ? real(1) : mbar;
      if (!(scale > real(0)))
        continue;
      w = std::pow(std::max(m[i] / scale, real(1e-3)), spec.power);
    }
    weight[i] = std::isfinite(w) ? w : real(0);
  }
}

/// One `edge_strain` row per edge at its grown goal length.
inline void init_flat_growth(const kusama::stencil_table &t, const std::vector<vec3> &x,
                             const std::vector<vertex_chart> &charts, const flat_growth_target &spec, real w,
                             std::vector<projection_constraint::ptr> &constraints,
                             std::vector<sim_block::ptr> blocks) {
  if (!(w > real(0)) || static_cast<int>(x.size()) != t.n_verts || charts.size() != x.size())
    return;
  std::vector<real> wt;
  std::vector<vec3> dir;
  flat_growth_field(charts, spec, wt, dir);
  const real iso = std::clamp(spec.iso, real(0), real(1));
  for (int e = 0; e < t.n_edges; ++e) {
    const std::array<int, 2> v = growth_detail::edge_ends(t, e);
    if (v[0] < 0)
      continue;
    const size_t i = static_cast<size_t>(v[0]), j = static_cast<size_t>(v[1]);
    const vec3 d = x[j] - x[i];
    const real l0 = d.norm();
    if (!(l0 > real(0)))
      continue;
    const vec3 a = dir[i] + (dir[i].dot(dir[j]) < real(0) ? vec3(-dir[j]) : dir[j]);
    const real an = a.norm();
    const real c = an > real(0) ? d.dot(a) / (l0 * an) : real(0);
    const real we = real(0.5) * (wt[i] + wt[j]);
    auto row = edge_strain::create(std::vector<index_t>{v[0], v[1]}, x, w, blocks);
    row->set_goal_length(l0 * (real(1) + spec.g * we * ((real(1) - iso) * c * c + iso)));
    constraints.push_back(row);
  }
}

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif
