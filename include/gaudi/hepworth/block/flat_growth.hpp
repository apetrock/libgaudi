#ifndef GAUDI_HEPWORTH_BLOCK_FLAT_GROWTH_HPP
#define GAUDI_HEPWORTH_BLOCK_FLAT_GROWTH_HPP

/// Flat-direction growth: excess length along an axis set from the loose
/// (near zero) principal direction. Per edge, goal length
/// `l0 (1 + g w_e (dhat0 . a_e)^2)`: `a_e` the growth axis of its ends,
/// `cos(angle) e_lo + sin(angle) e_hi` (sign-free; `angle` 0 is the loose
/// direction, `pi / 2` the tight one), `w_e` the mean of its ends' policy
/// weights
/// `smoothstep(1 - |k_lo| / |k_hi|; lo, hi) (|k_hi| / mean |k_hi|)^power`.
/// Round points (`k_lo ~ k_hi`) get none; `power > 0` keeps nearly flat
/// regions, where the direction is noise, out. On a torus the loose
/// direction at the top and bottom of the tube is the hoop, so those circles
/// lengthen and the top ruffles (the vessel-rim buckle). Rebuilt each frame
/// against the frame-start mesh and charts; positions only.

#include "gaudi/hepworth/block/mobius_growth.hpp"
#include "gaudi/hepworth/block/shape_objectives.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace gaudi {
namespace hepworth {
namespace block {

/// `g`: fractional length added per frame along the growth axis at full
/// weight; `lo`, `hi`: the smoothstep on `1 - |k_lo| / |k_hi|`; `power`:
/// exponent on `|k_hi|` against its mean; `angle`: the axis, radians from
/// the loose principal direction toward the tight one; `winding`: sets
/// `angle` at startup from a closed line on the fitted torus (resolved by
/// the caller, as the shear's).
struct flat_growth_target {
  real g = 0.05;
  real lo = 0.5;
  real hi = 0.9;
  real power = 1;
  real angle = 0;
  stretch_winding winding{};
};

/// Per vertex: the policy weight and the loose direction in world (zero
/// where the chart is missing).
inline void flat_growth_field(const std::vector<vertex_chart> &charts, const flat_growth_target &spec,
                              std::vector<real> &weight, std::vector<vec3> &dir) {
  const size_t n = charts.size();
  weight.assign(n, real(0));
  dir.assign(n, vec3::Zero());
  real kbar = 0;
  int nk = 0;
  for (const vertex_chart &c : charts)
    if (c.ok) {
      kbar += std::abs(growth_detail::principals_of(c.h).k_hi);
      ++nk;
    }
  kbar = nk > 0 ? kbar / nk : real(0);
  auto smooth = [](real t, real a, real b) {
    if (!(b > a))
      return t >= a ? real(1) : real(0);
    const real s = std::clamp((t - a) / (b - a), real(0), real(1));
    return s * s * (real(3) - real(2) * s);
  };
  for (size_t i = 0; i < n; ++i) {
    const vertex_chart &c = charts[i];
    if (!c.ok)
      continue;
    const growth_detail::principals p = growth_detail::principals_of(c.h);
    if (!(std::abs(p.k_hi) > real(0)))
      continue;
    real w = smooth(real(1) - std::abs(p.k_lo) / std::abs(p.k_hi), spec.lo, spec.hi);
    if (spec.power != real(0) && kbar > real(0))
      w *= std::pow(std::abs(p.k_hi) / kbar, spec.power);
    weight[i] = std::isfinite(w) ? w : real(0);
    const vec2 e = std::cos(spec.angle) * p.e_lo + std::sin(spec.angle) * p.e_hi;
    dir[i] = (e[0] * c.T + e[1] * c.B).normalized();
  }
}

/// One `edge_strain` row per edge at its flat-direction goal length.
inline void init_flat_growth(const kusama::stencil_table &t, const std::vector<vec3> &x,
                             const std::vector<vertex_chart> &charts, const flat_growth_target &spec, real w,
                             std::vector<projection_constraint::ptr> &constraints,
                             std::vector<sim_block::ptr> blocks) {
  if (!(w > real(0)) || static_cast<int>(x.size()) != t.n_verts || charts.size() != x.size())
    return;
  std::vector<real> wt;
  std::vector<vec3> dir;
  flat_growth_field(charts, spec, wt, dir);
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
    row->set_goal_length(l0 * (real(1) + spec.g * we * c * c));
    constraints.push_back(row);
  }
}

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif
