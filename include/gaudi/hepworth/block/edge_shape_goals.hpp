#ifndef GAUDI_HEPWORTH_BLOCK_EDGE_SHAPE_GOALS_HPP
#define GAUDI_HEPWORTH_BLOCK_EDGE_SHAPE_GOALS_HPP

/// Goal fields for `edge_shape_fair`, the edge ports of the vertex `H`
/// operators. Each goal is a local step re-projected from the live store
/// every PD iterate (as in ARAP) and returns a tensor per edge in that edge's
/// frame; each maps a store already at its target to itself.
/// - `goal_zero`: smoothing (with the Dirichlet and Dupin rows);
/// - `goal_fixed`: test targets;
/// - `goal_curvature_growth`: directional growth of the normal curvature;
/// - `goal_shape_arap`: the frame-start tensors with their lines turned by
///   the best turn of each vertex one-ring (ARAP in shape space);
/// - `goal_principal_mean`: each principal curvature to the neighbors'
///   normal curvature along its own direction (the old frame minimizer);
/// - `goal_mean_tension`: surface tension in curvature space, the mean
///   curvature to a constant or to the store's mean (`tension_target`).

#include "gaudi/hepworth/block/edge_shape_fair.hpp"
#include "gaudi/hepworth/block/shape_objectives.hpp"
#include "gaudi/hepworth/block/shell_constraints.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <memory>
#include <vector>

namespace gaudi {
namespace hepworth {
namespace block {

/// Normal curvature of edge `j`'s tensor along world direction `d` (projected
/// into `j`'s plane).
inline real edge_normal_curvature(const edge_frames &fr, const vec3 &hj, int j, const vec3 &d) {
  const vec3 &t = fr.t[static_cast<size_t>(j)], &b = fr.b[static_cast<size_t>(j)];
  const real u = d.dot(t), v = d.dot(b);
  const real l2 = u * u + v * v;
  if (!(l2 > real(1e-16)))
    return std::numeric_limits<real>::quiet_NaN();
  return (hj[0] * u * u + 2 * hj[1] * u * v + hj[2] * v * v) / l2;
}

/// Principal pair of one tensor ordered by magnitude: `k_hi` along
/// `(hx, hy)`, `|k_hi| >= |k_lo|`, so ratios are `|k_hi| / |k_lo| >= 1` and
/// do not flip with the sign of the curvature.
struct edge_principals {
  real k_hi = 0, k_lo = 0, hx = 1, hy = 0, lx = 0, ly = 1;
  bool ok = false;
};

inline edge_principals principals_of(const vec3 &h) {
  edge_principals p;
  const shape_w W = shape_w_from_coeffs(h);
  if (!W.ok)
    return p;
  const monge_axes a = monge_axes_of(W);
  p = {a.k_hi, a.k_lo, a.hx, a.hy, a.lx, a.ly, true};
  return p;
}

inline vec3 coeffs_of_principals(const edge_principals &p, const vec3 &fallback) {
  const shape_w out = shape_from_principals(p.k_hi, p.k_lo, p.hx, p.hy, p.lx, p.ly);
  return out.ok ? coeffs_of_shape(out) : fallback;
}

/// Signed principals `k0 >= k1` of `(k11, k12, k22)`.
inline void signed_principals(const vec3 &h, real &k0, real &k1) {
  const real H = real(0.5) * (h[0] + h[2]);
  const real r = std::hypot(real(0.5) * (h[0] - h[2]), h[1]);
  k0 = H + r;
  k1 = H - r;
}

/// `q = (r - 1) / (r + 1)` of one edge tensor (`curvature_q`).
inline real edge_q(const vec3 &h) {
  real k0, k1;
  signed_principals(h, k0, k1);
  return curvature_q(k0, k1);
}

/// `|k_hi| / |k_lo|` of one edge tensor; infinite on parabolic edges.
inline real edge_ratio(const vec3 &h) { return ratio_of_q(edge_q(h)); }

/// `(H, d, tau)`, the mean and the two deviator parts, set to anisotropy
/// `q`: the larger of `|H|` and `|D| = |(d, tau)|` stays and the smaller
/// becomes `q` times it, signs and axes kept. Elliptic (`|D| <= |H|`) keeps
/// `H`; saddle keeps `D`. A round point with `q > 0` has no axes and stays.
inline void set_anisotropy(real &H, real &d, real &tau, real q) {
  const real D = std::hypot(d, tau), aH = std::abs(H);
  if (D <= aH * (real(1) + real(1e-12))) {
    if (!(D > real(1e-12) * aH)) {
      if (q == real(0))
        d = tau = real(0);
      return;
    }
    const real s = q * aH / D;
    d *= s;
    tau *= s;
  } else {
    H = std::copysign(q * D, H);
  }
}

/// The tensor with principals `k0`, `k1` and the `k0` axis at `phi` radians
/// in the edge frame.
inline vec3 coeffs_at_angle(real k0, real k1, real phi) {
  const real H = real(0.5) * (k0 + k1), d = real(0.5) * (k0 - k1);
  const real c = std::cos(2 * phi), s = std::sin(2 * phi);
  return vec3(H + d * c, d * s, H - d * c);
}

/// The `q` an `anisotropy_target`'s ratio asks for on an edge whose
/// frame-start tensor is `h0`; NaN when the ratio is off.
inline real target_q(const anisotropy_target &spec, const vec3 &h0) {
  const real v = spec.ratio;
  const real nan = std::numeric_limits<real>::quiet_NaN();
  if (!std::isfinite(v))
    return nan;
  if (!spec.relative)
    return spec.log ? std::tanh(real(0.5) * std::abs(v)) : (v > real(0) ? q_of_ratio(v) : nan);
  // The frame-start ratio is capped at 1000 so parabolic edges respond too.
  const real r0 = std::min(ratio_of_q(edge_q(h0)), real(1000));
  if (spec.log)
    return std::tanh(real(0.5) * std::max(std::log(r0) + v, real(0)));
  if (!(v > real(0)))
    return nan;
  return q_of_ratio(std::max(r0 * v, real(1)));
}

/// Anisotropy knobs (`anisotropy_target`) on one tensor `h`, measured against
/// its frame-start value `h0` (the rest, as ARAP keeps rest spokes). The
/// deviator is `h`'s; the mean is `h`'s unless `mean` is set; then `ratio`
/// keeps the larger of the two. A tensor at its target maps to itself, so
/// re-projecting every iterate does not compound the edit, and the defaults
/// return `h`.
inline vec3 anisotropy_edge(const vec3 &h, const vec3 &h0, const anisotropy_target &spec) {
  if (!h.allFinite())
    return h;
  const vec3 r = h0.allFinite() ? h0 : h;
  real d = real(0.5) * (h[0] - h[2]), tau = h[1];
  real H = std::isfinite(spec.mean) ? spec.mean * real(0.5) * (r[0] + r[2])
                                    : real(0.5) * (h[0] + h[2]);
  const real q = target_q(spec, r);
  if (std::isfinite(q))
    set_anisotropy(H, d, tau, q);
  return vec3(H + d, tau, H - d);
}

inline edge_goal_fn goal_zero() {
  return [](const vecX &, std::vector<vec3> &g) { std::fill(g.begin(), g.end(), vec3::Zero()); };
}

/// The live store itself: a pin toward it is a proximal term, damping the
/// store's step toward its previous iterate, and zero once the iterates
/// settle.
inline edge_goal_fn goal_live(sim_block::ptr store) {
  return [store](const vecX &q, std::vector<vec3> &g) {
    for (size_t e = 0; e < g.size(); ++e)
      g[e] = store->get_vec3(static_cast<index_t>(e), q);
  };
}

inline edge_goal_fn goal_fixed(std::vector<vec3> targets) {
  return [targets = std::move(targets)](const vecX &, std::vector<vec3> &g) {
    for (size_t e = 0; e < g.size() && e < targets.size(); ++e)
      g[e] = targets[e];
  };
}

/// `h0` grown per `curvature_growth_target` (`angle` resolved): `D h0 D` in
/// the edge frame, `D = sqrt(1 + x) u u^T + sqrt(1 + y) v v^T`, `u` at `angle`
/// from `h0`'s k_hi axis toward its k_lo one and `v` across. With `power`,
/// each fraction is times `(|k| / k_bar)^power` of the normal curvature `k`
/// along its axis (`|k| / k_bar` floored at 1e-3 so a negative power stays
/// finite, as in `stretch_axis_scale`). An umbilic (no axes) scales by
/// `1 + (x + y) / 2` with `|H|` for `k`. Fractions at or below -1 are
/// clamped to bend nearly flat.
inline vec3 curvature_growth_tensor(const vec3 &h0, const curvature_growth_target &spec,
                                    real k_bar = 1) {
  if (!h0.allFinite())
    return h0;
  const bool fed = std::isfinite(spec.power) && spec.power != real(0) && k_bar > real(0) &&
                   std::isfinite(k_bar);
  auto gain = [&](real fraction, real k) {
    fraction = std::isfinite(fraction) ? fraction : real(0);
    if (fed)
      fraction *= std::pow(std::max(std::abs(k) / k_bar, real(1e-3)), spec.power);
    return std::max(real(1) + fraction, real(1e-6));
  };
  Eigen::Matrix<real, 2, 2> H, D;
  H << h0[0], h0[1], h0[1], h0[2];
  const shape_w axes = shape_with_axes(shape_w_from_coeffs(h0));
  if (!axes.ok || !(std::abs(axes.k1 - axes.k2) > real(1e-9) * (std::abs(axes.k1) + std::abs(axes.k2)))) {
    const real k = real(0.5) * (h0[0] + h0[2]);
    return real(0.5) * (gain(spec.x, k) + gain(spec.y, k)) * h0;
  }
  const monge_axes ax = monge_axes_of(axes);
  const real a = std::isfinite(spec.angle) ? spec.angle : real(0.25) * real(M_PI);
  const vec2 u = std::cos(a) * vec2(ax.hx, ax.hy) + std::sin(a) * vec2(ax.lx, ax.ly);
  const vec2 v(-u[1], u[0]);
  const real gx = gain(spec.x, u.dot(H * u));
  const real gy = gain(spec.y, v.dot(H * v));
  D = std::sqrt(gx) * u * u.transpose() + std::sqrt(gy) * v * v.transpose();
  const Eigen::Matrix<real, 2, 2> G = D * H * D;
  return vec3(G(0, 0), real(0.5) * (G(0, 1) + G(1, 0)), G(1, 1));
}

/// The store's mean `(|k1| + |k2|) / 2` over finite edges (0 if none).
inline real mean_abs_curvature(const std::vector<vec3> &h) {
  real s = 0;
  int n = 0;
  for (const vec3 &k : h) {
    if (!k.allFinite())
      continue;
    const real H = real(0.5) * (k[0] + k[2]);
    const real d = std::sqrt(real(0.25) * (k[0] - k[2]) * (k[0] - k[2]) + k[1] * k[1]);
    s += real(0.5) * (std::abs(H + d) + std::abs(H - d));
    ++n;
  }
  return n > 0 ? s / real(n) : real(0);
}

/// The grown frame-start store, held for the frame (`goal_fixed`), `power`
/// against the frame-start mean `|k|`. Use with tensor-metric pin rows.
inline edge_goal_fn goal_curvature_growth(const curvature_growth_target &spec,
                                          const std::vector<vec3> &h0) {
  const real k_bar = spec.power != real(0) ? mean_abs_curvature(h0) : real(1);
  std::vector<vec3> T(h0.size());
  for (size_t e = 0; e < h0.size(); ++e)
    T[e] = curvature_growth_tensor(h0[e], spec, k_bar);
  return goal_fixed(std::move(T));
}

/// Shape-space ARAP on vertex one-rings: each edge's frame-start tensor
/// `h0`, its curvature lines turned by the best turn of the rings it is in.
/// The anisotropic part as one complex number, `z = (k11 - k22) / 2 + i
/// k12`, has magnitude `(k0 - k1) / 2` and phase twice the `k0` axis angle;
/// turning the lines by `psi` multiplies it by `e^{2 i psi}`. The SVD (polar)
/// of `sum S S0^T` over a ring is then closed form: the turn is the phase of
/// `sum z conj(z0)`, which is the same in every edge frame, so the ring needs
/// no transport. Each vertex fits its ring; each edge takes the mean turn of
/// its two ends. The row resists any change of how anisotropic the surface
/// is and of the pattern of its lines, but lets a ring turn its lines
/// together. The mean curvature is `h0`'s when `hold_mean`, else the live
/// store's.
inline edge_goal_fn goal_shape_arap(sim_block::ptr store, const kusama::stencil_table *tb,
                                    std::vector<vec3> h0, bool hold_mean) {
  // The two vertices of each edge.
  std::vector<std::array<int, 2>> ends(static_cast<size_t>(tb->n_edges), {-1, -1});
  for (int f = 0; f < tb->n_faces; ++f)
    for (int i = 0; i < 3; ++i) {
      const int e = tb->face_edges[static_cast<size_t>(f)][static_cast<size_t>(i)];
      if (e >= 0 && e < tb->n_edges)
        ends[static_cast<size_t>(e)] = {tb->ids[static_cast<size_t>(f)][static_cast<size_t>(i)],
                                        tb->ids[static_cast<size_t>(f)][static_cast<size_t>((i + 1) % 3)]};
    }
  const int nV = tb->n_verts;
  return [store, hold_mean, nV, ends = std::move(ends),
          h0 = std::move(h0)](const vecX &q, std::vector<vec3> &g) {
    using cplx = std::complex<real>;
    const size_t nE = g.size();
    std::vector<vec3> h(nE);
    std::vector<cplx> z(nE), z0(nE), turn(nE, cplx(0, 0));
    std::vector<cplx> ring(static_cast<size_t>(std::max(nV, 0)), cplx(0, 0));
    for (size_t e = 0; e < nE; ++e) {
      h[e] = store->get_vec3(static_cast<index_t>(e), q);
      const vec3 &r = e < h0.size() ? h0[e] : h[e];
      z[e] = cplx(real(0.5) * (h[e][0] - h[e][2]), h[e][1]);
      z0[e] = cplx(real(0.5) * (r[0] - r[2]), r[1]);
      if (e < ends.size())
        for (int v : ends[e])
          if (v >= 0 && v < nV)
            ring[static_cast<size_t>(v)] += z[e] * std::conj(z0[e]);
    }
    for (cplx &s : ring)
      s = std::abs(s) > real(1e-30) ? s / std::abs(s) : cplx(0, 0);
    for (size_t e = 0; e < nE && e < ends.size(); ++e)
      for (int v : ends[e])
        if (v >= 0 && v < nV)
          turn[e] += ring[static_cast<size_t>(v)];
    for (size_t e = 0; e < nE; ++e) {
      const cplx u = std::abs(turn[e]) > real(1e-30) ? turn[e] / std::abs(turn[e]) : cplx(1, 0);
      const cplx zg = u * z0[e];
      const vec3 &r = e < h0.size() ? h0[e] : h[e];
      const real H = real(0.5) * (hold_mean ? r[0] + r[2] : h[e][0] + h[e][2]);
      g[e] = vec3(H + zg.real(), zg.imag(), H - zg.real());
    }
  };
}

/// Each principal value to the mean of the face neighbors' normal curvature
/// along the edge's own principal directions (the frame minimizer). At its
/// fixed point the goal is the store itself.
inline edge_goal_fn goal_principal_mean(sim_block::ptr store,
                                        std::shared_ptr<const edge_frames> fr) {
  return [store, fr](const vecX &q, std::vector<vec3> &g) {
    const size_t nE = g.size();
    std::vector<vec3> h(nE);
    for (size_t e = 0; e < nE; ++e)
      h[e] = store->get_vec3(static_cast<index_t>(e), q);
    for (size_t e = 0; e < nE; ++e) {
      g[e] = h[e];
      const shape_w W = shape_w_from_coeffs(h[e]);
      if (!W.ok)
        continue;
      const vec3 d1 = W.e1x * fr->t[e] + W.e1y * fr->b[e];
      const vec3 d2 = W.e2x * fr->t[e] + W.e2y * fr->b[e];
      real s1 = 0, s2 = 0;
      int n = 0;
      for (int j : fr->nbr[e]) {
        const real c1 = edge_normal_curvature(*fr, h[static_cast<size_t>(j)], j, d1);
        const real c2 = edge_normal_curvature(*fr, h[static_cast<size_t>(j)], j, d2);
        if (!std::isfinite(c1) || !std::isfinite(c2))
          continue;
        s1 += c1;
        s2 += c2;
        ++n;
      }
      if (n == 0)
        continue;
      const shape_w out = shape_from_principals(s1 / n, s2 / n, W.e1x, W.e1y, W.e2x, W.e2y);
      if (out.ok)
        g[e] = coeffs_of_shape(out);
    }
  };
}

/// Surface tension in curvature space. Area is not a curvature quantity, but
/// its first variation is the mean curvature; balanced against pressure it is
/// Young-Laplace, `H = dp / (2 gamma)`.
///   `constant`: every edge's mean curvature to `h0` (constant mean curvature,
///               bubbles; the store's sign: `k > 0` bends toward the normal,
///               so an outward-normal sphere of radius `r` has `H = -1 / r`).
///   `mean`:     to the store's mean over the edges (constant-mean-curvature
///               fairing: the curvature-space analogue of volume-preserving
///               mean-curvature flow).
struct tension_target {
  /// `volume` (frames mode only; the edge rows read it as `mean`): the
  /// frame mean flattened by `volume_gain` times the relative volume deficit
  /// against the first frame's.
  enum class mode { constant, mean, volume };
  mode kind = mode::mean;
  real h0 = 0;
  real volume_gain = 1;
};

/// Each edge tensor keeps its trace-free part and moves its mean curvature
/// `(k11 + k22) / 2` to the target: `g = h + (H* - H_e)(1, 0, 1)`. Linear in
/// the store.
inline edge_goal_fn goal_mean_tension(sim_block::ptr store, tension_target spec) {
  return [store, spec](const vecX &q, std::vector<vec3> &g) {
    const size_t nE = g.size();
    std::vector<vec3> h(nE);
    real sum = 0;
    size_t n = 0;
    for (size_t e = 0; e < nE; ++e) {
      h[e] = store->get_vec3(static_cast<index_t>(e), q);
      if (h[e].allFinite()) {
        sum += real(0.5) * (h[e][0] + h[e][2]);
        ++n;
      }
    }
    const real target = spec.kind == tension_target::mode::constant ? spec.h0 : (n ? sum / real(n) : real(0));
    for (size_t e = 0; e < nE; ++e) {
      g[e] = h[e];
      if (!h[e].allFinite())
        continue;
      const real d = target - real(0.5) * (h[e][0] + h[e][2]);
      g[e][0] += d;
      g[e][2] += d;
    }
  };
}

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif
