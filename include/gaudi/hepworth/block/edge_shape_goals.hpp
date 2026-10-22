#ifndef GAUDI_HEPWORTH_BLOCK_EDGE_SHAPE_GOALS_HPP
#define GAUDI_HEPWORTH_BLOCK_EDGE_SHAPE_GOALS_HPP

/// Goal fields for `edge_shape_fair`, the edge ports of the vertex `H`
/// operators. Each goal is a local step re-projected from the live store
/// every PD iterate (as in ARAP) and returns a tensor per edge in that edge's
/// frame; each maps a store already at its target to itself.
/// - `goal_zero`: smoothing (with `dirichlet` rows: the old `grad` and
///   `curvature_mean` rows);
/// - `goal_fixed`: test targets;
/// - `goal_anisotropy`: ratio and mean; `goal_curvature_shear`: torsion
///   (`anisotropy_target`), measured against the frame-start store;
/// - `goal_willmore_eq`: Willmore flow with a feature-size equalizer
///   (`willmore_target`), free bands left to the position rows;
/// - `goal_shape_arap`: the frame-start tensors with their lines turned by
///   the best turn of each vertex one-ring (ARAP in shape space);
/// - `goal_ratio_smooth`: each edge's anisotropy to the mean over its face
///   neighbors (the old `curvature_ratio`);
/// - `goal_principal_mean`: each principal curvature to the neighbors'
///   normal curvature along its own direction (the old frame minimizer).

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

/// The closest point to `h` of the set an `anisotropy_target` fixes, in the
/// tensor (Frobenius) metric `k11^2 + 2 k12^2 + k22^2 = 2 H^2 + 2 |delta|^2`
/// (`delta = ((k11 - k22) / 2, k12)`), which is round in `(H, delta)`. The set
/// is fixed for the frame (everything it reads comes from `h0`), so with
/// tensor-metric pin rows the goal is an exact projection and the row is
/// variational:
/// - `ratio`: the cones `|delta| = q |H|` and `|H| = q |delta|` with `q`
///   from `h0`; `delta`'s own axes are kept (so a torsion row turning them
///   is not resisted) and `(H, |delta|)` goes to the nearest of the rays
///   `(+-1, q)`, `(+-q, 1)`;
/// - `mean`: `H` is `mean H0`.
inline vec3 anisotropy_project(const vec3 &h, const vec3 &h0, const anisotropy_target &spec) {
  if (!h.allFinite())
    return h;
  const vec3 r = h0.allFinite() ? h0 : h;
  const bool fixed_mean = std::isfinite(spec.mean);
  const real q = target_q(spec, r);
  const real H0 = real(0.5) * (r[0] + r[2]);
  real H = fixed_mean ? spec.mean * H0 : real(0.5) * (h[0] + h[2]);
  vec2 delta(real(0.5) * (h[0] - h[2]), h[1]);
  const vec2 u0(real(0.5) * (r[0] - r[2]), r[1]);
  if (std::isfinite(q)) {
    // The deviator's own direction (h0's at an umbilic).
    vec2 dir = delta;
    if (!(dir.norm() > real(1e-12)))
      dir = u0.norm() > real(1e-12) ? u0 : vec2(1, 0);
    dir.normalize();
    const real s = delta.norm();
    // The set is both cones, |delta| = q |H| (|H| the larger) and
    // |H| = q |delta| (|delta| the larger): four rays in the (H, |delta|)
    // half-plane.
    real Ht = 0, Dt = 0;
    if (fixed_mean) {
      Ht = H;
      const real d1 = q * std::abs(H);
      const real d2 = q > real(0) ? std::abs(H) / q : std::numeric_limits<real>::infinity();
      Dt = std::abs(s - d1) <= std::abs(s - d2) ? d1 : d2;
    } else {
      const real len = std::sqrt(real(1) + q * q);
      const std::array<vec2, 4> rays{vec2(1, q) / len, vec2(-1, q) / len, vec2(q, 1) / len,
                                     vec2(-q, 1) / len};
      real best = 0;
      for (const vec2 &u : rays) {
        const real l = H * u[0] + s * u[1];
        if (l > best) {
          best = l;
          Ht = l * u[0];
          Dt = l * u[1];
        }
      }
    }
    H = Ht;
    delta = Dt * dir;
  }
  return vec3(H + delta[0], delta[1], H - delta[0]);
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

/// `h0` is the frame-start store (the live reading at prepare). Use with
/// tensor-metric pin rows (`edge_shape_fair(..., tensor_metric = true)`).
inline edge_goal_fn goal_anisotropy(sim_block::ptr store, anisotropy_target spec,
                                    std::vector<vec3> h0) {
  return [store, spec, h0 = std::move(h0)](const vecX &q, std::vector<vec3> &g) {
    for (size_t e = 0; e < g.size(); ++e) {
      const vec3 h = store->get_vec3(static_cast<index_t>(e), q);
      g[e] = anisotropy_project(h, e < h0.size() ? h0[e] : h, spec);
    }
  };
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

/// Shape-space ARAP under a shear (`curvature_shear_target`): each vertex
/// `v`'s ring of frame-start tensors `h0`, pushed through the shear `S_v`
/// in `v`'s step-start chart (`W -> S_v^-T W S_v^-1` on the world tensor,
/// read back in the edge frame), is turned by the ring's best common turn
/// (the phase of `sum z conj(z_s)`, as `goal_shape_arap`); each edge takes
/// the direction of the mean of its two ends' turned, pushed deviators, and
/// keeps its live deviator's magnitude and mean curvature, so the row holds
/// only the pattern of the lines and leaves the anisotropy to
/// `anisotropy_target`. A vertex without a chart pushes nothing (`S = I`).
/// `spec.angle` is the resolved axis (no winding lookup here).
inline edge_goal_fn goal_curvature_shear(sim_block::ptr store, const kusama::stencil_table *tb,
                                         const edge_frames &fr,
                                         const std::vector<vertex_chart> &charts,
                                         curvature_shear_target spec, std::vector<vec3> h0) {
  using cplx = std::complex<real>;
  const size_t nE = static_cast<size_t>(tb->n_edges);
  const int nV = tb->n_verts;
  std::vector<std::array<int, 2>> ends(nE, {-1, -1});
  for (int f = 0; f < tb->n_faces; ++f)
    for (int i = 0; i < 3; ++i) {
      const int e = tb->face_edges[static_cast<size_t>(f)][static_cast<size_t>(i)];
      if (e >= 0 && e < tb->n_edges)
        ends[static_cast<size_t>(e)] = {tb->ids[static_cast<size_t>(f)][static_cast<size_t>(i)],
                                        tb->ids[static_cast<size_t>(f)][static_cast<size_t>((i + 1) % 3)]};
    }
  // S_v^-1 per vertex (identity without a chart or a shear).
  const real g = std::isfinite(spec.shear) ? std::clamp(spec.shear, real(-0.99), real(0.99)) : real(0);
  const real sigma = std::sqrt((real(1) + g) / (real(1) - g));
  const real beta = std::isfinite(spec.angle) ? spec.angle : real(0.25) * real(M_PI);
  std::vector<mat3> S_inv(static_cast<size_t>(std::max(nV, 0)), mat3::Identity());
  for (int v = 0; v < nV && g != real(0); ++v) {
    shape_w W;
    vec3 T, B, N;
    if (!vertex_chart_at(&charts, v, W, T, B, N))
      continue;
    const shape_w axes = shape_with_axes(W);
    if (!axes.ok)
      continue;
    const monge_axes ax = monge_axes_of(axes);
    N.normalize();
    vec3 a = std::cos(beta) * (ax.hx * T + ax.hy * B) + std::sin(beta) * (ax.lx * T + ax.ly * B);
    a -= N * a.dot(N);
    if (!(a.squaredNorm() > real(1e-16)))
      continue;
    a.normalize();
    const vec3 ap = N.cross(a);
    S_inv[static_cast<size_t>(v)] =
        a * a.transpose() / sigma + sigma * ap * ap.transpose() + N * N.transpose();
  }
  // Each edge's frame-start deviator pushed through each end's shear.
  std::vector<std::array<cplx, 2>> zs(nE, {cplx(0, 0), cplx(0, 0)});
  for (size_t e = 0; e < nE && e < h0.size() && e < fr.t.size(); ++e) {
    const mat3 W0 = edge_world_tensor(h0[e], fr.t[e], fr.b[e]);
    for (int k = 0; k < 2; ++k) {
      const int v = ends[e][static_cast<size_t>(k)];
      const mat3 Si = v >= 0 && v < nV ? S_inv[static_cast<size_t>(v)] : mat3::Identity();
      const mat3 Wp = Si.transpose() * W0 * Si;
      const vec3 &t = fr.t[e], &b = fr.b[e];
      zs[e][static_cast<size_t>(k)] =
          cplx(real(0.5) * (t.dot(Wp * t) - b.dot(Wp * b)), t.dot(Wp * b));
    }
  }
  return [store, nV, ends = std::move(ends), zs = std::move(zs)](const vecX &q,
                                                                  std::vector<vec3> &goal) {
    const size_t n = std::min(goal.size(), zs.size());
    std::vector<vec3> h(n);
    std::vector<cplx> ring(static_cast<size_t>(std::max(nV, 0)), cplx(0, 0));
    for (size_t e = 0; e < n; ++e) {
      h[e] = store->get_vec3(static_cast<index_t>(e), q);
      const cplx z(real(0.5) * (h[e][0] - h[e][2]), h[e][1]);
      for (int k = 0; k < 2; ++k) {
        const int v = ends[e][static_cast<size_t>(k)];
        if (v >= 0 && v < nV)
          ring[static_cast<size_t>(v)] += z * std::conj(zs[e][static_cast<size_t>(k)]);
      }
    }
    for (cplx &r : ring)
      r = std::abs(r) > real(1e-30) ? r / std::abs(r) : cplx(1, 0);
    for (size_t e = 0; e < n; ++e) {
      cplx zg(0, 0);
      int c = 0;
      for (int k = 0; k < 2; ++k) {
        const int v = ends[e][static_cast<size_t>(k)];
        if (v >= 0 && v < nV) {
          zg += ring[static_cast<size_t>(v)] * zs[e][static_cast<size_t>(k)];
          ++c;
        }
      }
      const cplx z(real(0.5) * (h[e][0] - h[e][2]), h[e][1]);
      zg = c > 0 && std::abs(zg) > real(1e-30) ? zg / std::abs(zg) * std::abs(z) : z;
      const real H = real(0.5) * (h[e][0] + h[e][2]);
      goal[e] = vec3(H + zg.real(), zg.imag(), H - zg.real());
    }
  };
}

/// `n` passes of neighbor averaging of an edge tensor field: each edge's
/// world tensor averaged with its face neighbors', read back in its own
/// frame. `n` passes blur over about `0.45 sqrt(n)` edge lengths.
inline std::vector<vec3> blur_edge_tensors(const edge_frames &fr, std::vector<vec3> h, int n) {
  const size_t nE = h.size();
  std::vector<mat3> W(nE), Wn(nE);
  for (int it = 0; it < n; ++it) {
    for (size_t e = 0; e < nE; ++e)
      W[e] = edge_world_tensor(h[e], fr.t[e], fr.b[e]);
    for (size_t e = 0; e < nE; ++e) {
      mat3 s = W[e];
      int c = 1;
      for (int j : fr.nbr[e])
        if (j >= 0 && static_cast<size_t>(j) < nE) {
          s += W[static_cast<size_t>(j)];
          ++c;
        }
      Wn[e] = s / real(c);
    }
    for (size_t e = 0; e < nE; ++e) {
      const vec3 &t = fr.t[e], &b = fr.b[e];
      h[e] = vec3(t.dot(Wn[e] * t), t.dot(Wn[e] * b), b.dot(Wn[e] * b));
    }
  }
  return h;
}

/// Passes for a blur radius of `r` mean edge lengths (`blur_edge_tensors`).
inline int blur_passes(real r) {
  if (!(r > real(0)) || !std::isfinite(r))
    return 0;
  return std::min(2000, static_cast<int>(std::ceil(r * r / (real(0.45) * real(0.45)))));
}

/// Screened blurs of an edge tensor field in Crane's filter form,
/// `blur(x, r) = (M + sigma L)^{-1} M x`: `L` the CR spin-2 operator
/// (`build_spin2`), `M` each edge's dual area times the tensor inner product,
/// `sigma = (r * mean edge)^2`. One factorization per radius; the cost does
/// not grow with the radius. Fields with no spin-2 energy (constant in a
/// flat region) pass unchanged.
class edge_tensor_blur {
public:
  edge_tensor_blur(const kusama::stencil_table &t, const std::vector<vec3> &x,
                   const std::vector<real> &radii) {
    const int nE = t.n_edges;
    const std::vector<real> c = kusama::weingarten::coupling_factor(t, x);
    const Eigen::Matrix<real, 3, 3> Q = kusama::weingarten::tensor_inner();
    std::vector<trip> mt;
    real len = 0;
    int nl = 0;
    for (int f = 0; f < t.n_faces; ++f)
      for (int i = 0; i < 3; ++i) {
        len += (x[static_cast<size_t>(t.ids[f][(i + 1) % 3])] - x[static_cast<size_t>(t.ids[f][i])])
                   .norm();
        ++nl;
      }
    const real l = nl > 0 ? len / nl : real(1);
    for (int e = 0; e < nE; ++e) {
      const real A = c[static_cast<size_t>(e)] / kusama::weingarten::k_coupling;
      for (int r = 0; r < 3; ++r)
        for (int s = 0; s < 3; ++s)
          if (Q(r, s) != real(0))
            mt.emplace_back(3 * e + r, 3 * e + s, A * Q(r, s));
    }
    _M = matS(3 * nE, 3 * nE);
    _M.setFromTriplets(mt.begin(), mt.end());
    const matS L = kusama::weingarten::build_spin2(t, x);
    for (real r : radii) {
      auto ldlt = std::make_shared<Eigen::SimplicialLDLT<matS>>();
      const matS A = _M + (r * l) * (r * l) * L;
      ldlt->compute(A);
      _solve.push_back(ldlt);
    }
  }

  /// `blur(h, radii[k])`.
  std::vector<vec3> apply(size_t k, const std::vector<vec3> &h) const {
    const size_t nE = h.size();
    vecX b(3 * static_cast<index_t>(nE));
    for (size_t e = 0; e < nE; ++e)
      b.segment<3>(3 * static_cast<index_t>(e)) = h[e];
    const vecX y = _solve[k]->solve(_M * b);
    std::vector<vec3> out(nE);
    for (size_t e = 0; e < nE; ++e)
      out[e] = y.segment<3>(3 * static_cast<index_t>(e));
    return out;
  }

  size_t size() const { return _solve.size(); }

private:
  matS _M;
  std::vector<std::shared_ptr<Eigen::SimplicialLDLT<matS>>> _solve;
};

/// The Willmore flow of a store: `-Dev(h) - mean * H(h) I` per edge.
inline std::vector<vec3> willmore_flow(const std::vector<vec3> &h0, real mean) {
  std::vector<vec3> v(h0.size());
  for (size_t e = 0; e < h0.size(); ++e) {
    const vec3 &h = h0[e];
    const real H = real(0.5) * (h[0] + h[2]), d = real(0.5) * (h[0] - h[2]);
    v[e] = vec3(-d - mean * H, -h[1], d - mean * H);
  }
  return v;
}

/// A `willmore_target` as a store goal: the feature-size equalizer.
/// Bands `B_i(x) = blur(x, r_i) - blur(x, r_{i+1})` (the last `blur(x,
/// r_n)`, the first starting at radius 0, the identity) sum to the identity.
/// The frame target is `T = h0 + sum g_i B_i(v0)`; the goal is
/// `h0 - E(h0 - T)` with `E = sum c_i B_i`, `c_i = 1` where `g_i != 0`: the
/// active bands moved to `T`, the free ones at their frame-start values,
/// held for the frame. All gains equal (one band included) is plain Willmore,
/// a fixed pin toward `T` with no blur. Returns an empty function when off.
inline edge_goal_fn goal_willmore_eq(sim_block::ptr store, const kusama::stencil_table &t,
                                     const std::vector<vec3> &x, const std::vector<vec3> &h0,
                                     const willmore_target &spec) {
  std::vector<willmore_band> bands = spec.bands;
  std::sort(bands.begin(), bands.end(),
            [](const willmore_band &a, const willmore_band &b) { return a.radius < b.radius; });
  if (!bands.empty() && bands.front().radius > real(0))
    bands.insert(bands.begin(), willmore_band{0, 0});
  for (willmore_band &b : bands)
    if (!std::isfinite(b.gain))
      b.gain = 0;
  bool any = false, equal = true;
  for (const willmore_band &b : bands) {
    any = any || b.gain != real(0);
    equal = equal && b.gain == bands.front().gain;
  }
  if (!any)
    return {};
  const std::vector<vec3> v0 = willmore_flow(h0, std::isfinite(spec.mean) ? spec.mean : real(0));
  const size_t nE = h0.size();
  if (equal) {
    std::vector<vec3> T(nE);
    for (size_t e = 0; e < nE; ++e)
      T[e] = h0[e] + bands.front().gain * v0[e];
    return goal_fixed(std::move(T));
  }
  // Blurs at every band start but the first (radius 0, the identity).
  std::vector<real> radii;
  for (size_t i = 1; i < bands.size(); ++i)
    radii.push_back(bands[i].radius);
  auto blur = std::make_shared<edge_tensor_blur>(t, x, radii);
  const size_t n = bands.size();
  // b_k(x) = blur(x, r_k), b_0 = x; B_i = b_i - b_{i+1}, B_{n-1} = b_{n-1}.
  auto blurs = [blur, n](const std::vector<vec3> &f) {
    std::vector<std::vector<vec3>> b(n);
    b[0] = f;
    for (size_t k = 1; k < n; ++k)
      b[k] = blur->apply(k - 1, f);
    return b;
  };
  const std::vector<std::vector<vec3>> bv = blurs(v0);
  std::vector<vec3> T = h0;
  for (size_t i = 0; i < n; ++i)
    for (size_t e = 0; e < nE; ++e)
      T[e] += bands[i].gain * (bv[i][e] - (i + 1 < n ? bv[i + 1][e] : vec3::Zero()));
  // E(f) = sum_k d_k b_k(f), d_k = c_k - c_{k-1}: only where presence changes.
  std::vector<real> d(n);
  for (size_t k = 0; k < n; ++k) {
    const real ck = bands[k].gain != real(0) ? real(1) : real(0);
    const real cp = k > 0 && bands[k - 1].gain != real(0) ? real(1) : real(0);
    d[k] = ck - cp;
  }
  auto filter = [blur, d, n](const std::vector<vec3> &f) {
    std::vector<vec3> out(f.size(), vec3::Zero());
    for (size_t k = 0; k < n; ++k) {
      if (d[k] == real(0))
        continue;
      const std::vector<vec3> bk = k == 0 ? f : blur->apply(k - 1, f);
      for (size_t e = 0; e < f.size(); ++e)
        out[e] += d[k] * bk[e];
    }
    return out;
  };
  // The goal S - E(S - T) taken at the frame-start store and held: a fixed
  // pin, so the row is variational. Within a frame the free bands are held
  // at their frame-start values; each frame re-reads them.
  const std::vector<vec3> ET = filter(T), EH = filter(h0);
  std::vector<vec3> g0(nE);
  for (size_t e = 0; e < nE; ++e)
    g0[e] = h0[e] - EH[e] + ET[e];
  return goal_fixed(std::move(g0));
}

/// Each edge's anisotropy `q` (`curvature_q`, the ratio squeezed into
/// [0, 1]) to the mean over itself and its face neighbors, at its own mean
/// or deviator and axes (`set_anisotropy`).
inline edge_goal_fn goal_ratio_smooth(sim_block::ptr store,
                                      std::shared_ptr<const edge_frames> fr) {
  return [store, fr](const vecX &q, std::vector<vec3> &g) {
    const size_t nE = g.size();
    std::vector<vec3> h(nE);
    std::vector<real> Q(nE);
    for (size_t e = 0; e < nE; ++e) {
      h[e] = store->get_vec3(static_cast<index_t>(e), q);
      Q[e] = edge_q(h[e]);
    }
    for (size_t e = 0; e < nE; ++e) {
      g[e] = h[e];
      if (!h[e].allFinite())
        continue;
      real sum = Q[e];
      int n = 1;
      for (int j : fr->nbr[e])
        if (std::isfinite(Q[static_cast<size_t>(j)])) {
          sum += Q[static_cast<size_t>(j)];
          ++n;
        }
      real H = real(0.5) * (h[e][0] + h[e][2]), d = real(0.5) * (h[e][0] - h[e][2]),
           tau = h[e][1];
      set_anisotropy(H, d, tau, sum / n);
      g[e] = vec3(H + d, tau, H - d);
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

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif
