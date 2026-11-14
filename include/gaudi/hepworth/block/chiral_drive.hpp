#ifndef GAUDI_HEPWORTH_BLOCK_CHIRAL_DRIVE_HPP
#define GAUDI_HEPWORTH_BLOCK_CHIRAL_DRIVE_HPP

/// Chiral curvature drive (`chiral_target`, MATH.md 19.22; Armon, Efrati,
/// Kupferman, Sharon 2011, `reference/armon2011-chiral-seed-pods.pdf`): per
/// directed spoke the normal curvature and twist of each vertex's target
/// tensor `b = h_ref + k0 w (e+ e+^T - e- e-^T)`, for `frame_rows`' spoke
/// target rows (`set_spoke_targets`, chord form). Edge-native, so no new row
/// type. Paired with a held metric (frame-0 edge lengths) it is Armon's
/// incompatible strip: a helicoid when narrow, a helix when wide.

#include "gaudi/hepworth/block/flat_growth.hpp"
#include "gaudi/hepworth/block/frame_rows.hpp"

#include <cmath>
#include <vector>

namespace gaudi {
namespace hepworth {
namespace block {

/// Per vertex the target tensor `(s11, s12, s22)` in its chart's `(T, B)`:
/// `h_ref` (the chart's own `h` when null or non-finite) plus the saddle,
/// the axis read off the reference. `h_first` (the first frame's tensor) is
/// what `k_max` caps the traceless change against in `rate` mode.
inline std::vector<vec3> chiral_tensors(const std::vector<vertex_chart> &charts, const std::vector<vec3> *h_ref,
                                        const chiral_target &spec, const std::vector<vec3> *h_first = nullptr) {
  const size_t n = charts.size();
  std::vector<vec3> out(n, vec3::Zero());
  std::vector<real> w;
  measure_weight(charts, spec.measure, spec.power, w);
  auto traceless = [](const vec3 &s) { return vec2(real(0.5) * (s[0] - s[2]), s[1]); };
  for (size_t v = 0; v < n; ++v) {
    const vertex_chart &c = charts[v];
    if (!c.ok || !c.h.allFinite()) {
      out[v] = c.h.allFinite() ? c.h : vec3::Zero();
      continue;
    }
    const vec3 ref = h_ref && v < h_ref->size() && (*h_ref)[v].allFinite() ? (*h_ref)[v] : c.h;
    const growth_detail::principals p = growth_detail::principals_of(ref);
    const real q = std::abs(p.H) + p.z > real(0) ? p.z / (std::abs(p.H) + p.z) : real(0);
    const bool principal =
        spec.axis == shear_axis::principal || (spec.axis == shear_axis::automatic && q >= spec.axis_q);
    const real axis = (principal ? std::atan2(p.e_hi[1], p.e_hi[0]) : real(0)) + spec.angle;
    vec3 S = vec3::Zero();
    for (int sgn : {1, -1}) {
      const real a = axis + real(sgn) * spec.opening;
      const real cs = std::cos(a), sn = std::sin(a);
      S += real(sgn) * vec3(cs * cs, cs * sn, sn * sn);
    }
    const real H = real(0.5) * (ref[0] + ref[2]);
    const vec2 z0 = traceless(ref);
    vec2 z = z0 + spec.k0 * w[v] * traceless(S);
    if (spec.keep_eigenvalues && z0.norm() > real(1e-3) * (std::abs(H) + z0.norm()) && z.norm() > real(0))
      z *= z0.norm() / z.norm();
    if (spec.mode == chiral_mode::rate && std::isfinite(spec.k_max) && h_first && v < h_first->size() &&
        (*h_first)[v].allFinite()) {
      const vec2 zf = traceless((*h_first)[v]);
      const vec2 d = z - zf;
      if (d.norm() > spec.k_max)
        z = zf + d * (spec.k_max / d.norm());
    }
    out[v] = vec3(H + z[0], z[1], H - z[0]);
  }
  return out;
}

/// Per directed spoke (`tp` order) the target normal curvature and twist of
/// its foot's chiral tensor along the spoke; spokes whose foot has no chart
/// keep their live reading (no drive).
inline void chiral_spoke_targets(const kusama::contact_net::topology &tp, const std::vector<vec3> &x,
                                 const std::vector<vertex_chart> &charts, const std::vector<vec3> *h_ref,
                                 const chiral_target &spec, std::vector<real> &k, std::vector<real> &tau,
                                 const std::vector<vec3> *h_first = nullptr) {
  std::vector<vec3> N(charts.size());
  for (size_t v = 0; v < charts.size(); ++v)
    N[v] = charts[v].N;
  frame_rows::spoke_readings(tp, x, N, k, tau);
  const std::vector<vec3> b = chiral_tensors(charts, h_ref, spec, h_first);
  for (size_t s = 0; s < tp.spokes.size(); ++s) {
    const size_t i = static_cast<size_t>(tp.spokes[s][0]), j = static_cast<size_t>(tp.spokes[s][1]);
    const vertex_chart &c = charts[i];
    if (!c.ok)
      continue;
    const vec3 dx = x[j] - x[i], nb = (N[i] + N[j]).normalized();
    const vec3 dxt = dx - nb * dx.dot(nb);
    real a = dxt.dot(c.T), bb = dxt.dot(c.B);
    const real r = std::hypot(a, bb);
    if (!(r > real(0)))
      continue;
    a /= r;
    bb /= r;
    const vec3 &S = b[i];
    k[s] = a * a * S[0] + 2 * a * bb * S[1] + bb * bb * S[2];
    tau[s] = -a * bb * S[0] + (a * a - bb * bb) * S[1] + a * bb * S[2];
  }
}

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif
