#pragma once

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <vector>

#include "gaudi/common.h"
#include "osculating_dupin_gap_generated.hpp"
#include "osculating_torus_ring.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

/// Special conformal transform about the origin,
///
///   S_b(x) = (x + b |x|²) / (1 + 2 b·x + |b|² |x|²),
///
/// an inversion, a translation by `b`, and the inversion again. It fixes the
/// origin with the identity differential, maps spheres to spheres, and the
/// transforms compose additively (`S_b ∘ S_c = S_{b+c}`), so `S_{-b}` inverts
/// `S_b`.
inline vec3 special_conformal(const vec3 &b, const vec3 &x) {
  const real x2 = x.squaredNorm();
  const real den = real(1) + real(2) * b.dot(x) + b.squaredNorm() * x2;
  return (x + b * x2) / den;
}

/// Scale factor `λ` of `S_b` at `x` (`dS_b = λ R`), its Jacobian, and
/// `∇σ` for `σ = log λ`.
inline real special_conformal_frame(const vec3 &b, const vec3 &x, mat3 &J, vec3 &grad_sigma) {
  const real x2 = x.squaredNorm();
  const real den = real(1) + real(2) * b.dot(x) + b.squaredNorm() * x2;
  const vec3 dden = real(2) * b + real(2) * b.squaredNorm() * x;
  const vec3 S = (x + b * x2) / den;
  J = (mat3::Identity() + real(2) * b * x.transpose()) / den - S * dden.transpose() / den;
  grad_sigma = -dden / den;
  return real(1) / den;
}

/// Local shape of a Dupin cyclide at the foot: the store's osculating torus
/// carried by `S_b` with `b` tangent at the foot. In the torus's principal
/// coordinates `(u, v, n)`, `b = (b_hi, b_lo, 0)`. To first order in `b` the
/// height gains exactly the two cross cubic terms of a Dupin cyclide,
///
///   δh = (k_lo − k_hi) b_hi u v² + (k_hi − k_lo) b_lo u² v,
///
/// and no ridge terms, so `b` is the cyclide's third-order freedom past the
/// torus. A normal `b` would only shift both curvatures, so it is not used.
struct dupin_ring_shape {
  torus_ring_shape torus;
  real b_hi = 0;
  real b_lo = 0;
  vec3 b() const { return vec3(b_hi, b_lo, real(0)); }
};

/// Landing of a point on the cyclide, in the foot's principal coordinates
/// `(u, v, n)`: pull it back by `S_{-b}`, drop it onto the torus graph, and
/// push the landing forward by `S_b`. `gap` is the signed offset along the
/// cyclide normal there.
inline bool dupin_land_local(const dupin_ring_shape &shape, const vec3 &p, vec3 &land, real &gap,
                             vec3 *normal = nullptr) {
  const vec3 b = shape.b();
  const vec3 q = special_conformal(-b, p);
  real h = 0, hu = 0, hv = 0;
  if (!osculating_torus_height_grad(shape.torus.k_hi, shape.torus.k_lo, q[0], q[1], h, hu, hv))
    return false;
  const vec3 qs(q[0], q[1], h);
  mat3 J;
  vec3 gs;
  special_conformal_frame(b, qs, J, gs);
  const vec3 n = (J * vec3(-hu, -hv, real(1))).normalized();
  land = special_conformal(b, qs);
  gap = (p - land).dot(n);
  if (normal)
    *normal = n;
  return land.allFinite() && std::isfinite(gap);
}

/// The cyclide's curvature tensor near the point `p` (principal coordinates),
/// as a 3×3 tangent tensor, with its unit normal there. The torus tensor at the
/// pulled-back point comes from the graph's Hessian, and the conformal law
/// carries it forward: principal directions map by `dS_b`, and each principal
/// curvature becomes `(κ − n·∇σ) / λ` in this convention, where a height
/// `½ κ t²` bends toward `+n`.
inline bool dupin_tensor_local(const dupin_ring_shape &shape, const vec3 &p, real step, mat3 &W,
                               vec3 &normal) {
  const vec3 b = shape.b();
  const vec3 q = special_conformal(-b, p);
  const real kh = shape.torus.k_hi, kl = shape.torus.k_lo;
  real h = 0, hu = 0, hv = 0;
  if (!osculating_torus_height_grad(kh, kl, q[0], q[1], h, hu, hv))
    return false;
  // Exact Hessian of the torus graph: the generated torus stage returns
  // (h, h_u, h_v) and its Jacobian, whose (h_u, h_v) rows over (u, v) are it.
  (void)step;
  Eigen::Matrix<real, 3, 1> th;
  Eigen::Matrix<real, 3, 4> Jt;
  dupin_generated::dupin_torus(Eigen::Matrix<real, 4, 1>(q[0], q[1], kh, kl), th, Jt);
  const real huu = Jt(1, 0);
  const real huv = real(0.5) * (Jt(1, 1) + Jt(2, 0));
  const real hvv = Jt(2, 1);
  if (!std::isfinite(huu) || !std::isfinite(huv) || !std::isfinite(hvv))
    return false;
  const real g = std::sqrt(real(1) + hu * hu + hv * hv);
  const vec3 n = vec3(-hu, -hv, real(1)) / g;
  // Orthonormal tangent frame at the torus point and the second form in it.
  // A tangent vector's first two coordinates are its graph parameters.
  const vec3 t1 = vec3(real(1), real(0), hu).normalized();
  const vec3 t2 = n.cross(t1);
  auto II = [&](const vec3 &a, const vec3 &c) {
    return (a[0] * (huu * c[0] + huv * c[1]) + a[1] * (huv * c[0] + hvv * c[1])) / g;
  };
  Eigen::Matrix<real, 2, 2> M;
  M << II(t1, t1), II(t1, t2), II(t2, t1), II(t2, t2);
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix<real, 2, 2>> es(M);
  if (es.info() != Eigen::Success)
    return false;
  const vec3 qs(q[0], q[1], h);
  mat3 J;
  vec3 gs;
  const real lambda = special_conformal_frame(b, qs, J, gs);
  normal = (J * n).normalized();
  const real shift = n.dot(gs);
  W.setZero();
  for (int i = 0; i < 2; ++i) {
    const vec3 d = es.eigenvectors()(0, i) * t1 + es.eigenvectors()(1, i) * t2;
    const vec3 dm = (J * d).normalized();
    const real k = (es.eigenvalues()[i] - shift) / lambda;
    W += k * dm * dm.transpose();
  }
  return W.allFinite() && normal.allFinite();
}

/// One edge of the foot's star: its end positions relative to the foot, its
/// store reading in the store's frame (edge direction low to high id, `b = n ×
/// t`), and its dual mass.
struct dupin_ring_edge {
  vec3 lo = vec3::Zero();
  vec3 hi = vec3::Zero();
  vec3 store = vec3::Zero();
  real mass = 0;
};

/// The store tensor a cyclide predicts on `edge`: its tensor at the edge
/// midpoint, read in the edge direction projected into the cyclide's tangent
/// plane there.
inline bool dupin_edge_reading(const dupin_ring_shape &shape, const torus_ring_frame &f, real offset,
                               const dupin_ring_edge &edge, real step, vec3 &reading) {
  const vec3 eh = shape.torus.hx * f.T + shape.torus.hy * f.B;
  const vec3 el = shape.torus.lx * f.T + shape.torus.ly * f.B;
  auto local = [&](const vec3 &d) { return vec3(d.dot(eh), d.dot(el), d.dot(f.N)); };
  const vec3 m = local(real(0.5) * (edge.lo + edge.hi) - offset * f.N);
  mat3 W;
  vec3 n;
  if (!dupin_tensor_local(shape, m, step, W, n))
    return false;
  vec3 t = local(edge.hi - edge.lo);
  t -= n * t.dot(n);
  if (!(t.squaredNorm() > real(1e-24)))
    return false;
  t.normalize();
  const vec3 bt = n.cross(t);
  reading = vec3(t.dot(W * t), t.dot(W * bt), bt.dot(W * bt));
  return reading.allFinite();
}

/// Landing point only, for the row targets: no gap, so no conformal Jacobian
/// and no normal.
inline bool land_point_on_posed_dupin(const dupin_ring_shape &shape, const torus_ring_frame &f, real offset,
                                      const vec3 &d, vec3 &land) {
  const vec3 eh = shape.torus.hx * f.T + shape.torus.hy * f.B;
  const vec3 el = shape.torus.lx * f.T + shape.torus.ly * f.B;
  const vec3 p = d - offset * f.N;
  const vec3 b = shape.b();
  const vec3 q = special_conformal(-b, vec3(p.dot(eh), p.dot(el), p.dot(f.N)));
  real h = 0;
  if (!osculating_torus_height(shape.torus.k_hi, shape.torus.k_lo, q[0], q[1], h))
    return false;
  const vec3 l = special_conformal(b, vec3(q[0], q[1], h));
  land = offset * f.N + l[0] * eh + l[1] * el + l[2] * f.N;
  return land.allFinite();
}

/// Landing of `d` (from the unposed foot) on the posed cyclide.
inline bool land_on_posed_dupin(const dupin_ring_shape &shape, const torus_ring_frame &f, real offset,
                                const vec3 &d, vec3 &land, real &gap) {
  const vec3 eh = shape.torus.hx * f.T + shape.torus.hy * f.B;
  const vec3 el = shape.torus.lx * f.T + shape.torus.ly * f.B;
  const vec3 p = d - offset * f.N;
  vec3 l;
  if (!dupin_land_local(shape, vec3(p.dot(eh), p.dot(el), p.dot(f.N)), l, gap))
    return false;
  land = offset * f.N + l[0] * eh + l[1] * el + l[2] * f.N;
  return land.allFinite();
}

/// Knobs of the osculating Dupin cyclide row. The row weights themselves
/// (`w` for the position rows, `w_shape` for the store rows) are passed to
/// `init_osculating_dupin_ring` separately.
struct osculating_dupin_ring_config {
  /// Weight of each 2-ring point relative to a 1-ring point, in both the fit
  /// and the position rows. The 2-ring is what sees third order in the
  /// positions, so it constrains `b` from the mesh side. 0 drops the 2-ring:
  /// `b` then comes from the store's star alone, the rows touch only the
  /// 1-ring, and the row is cheaper.
  real two_ring = 0.5;
  /// Gauss-Newton iterations per fit, both the frame-start shape fit and the
  /// pose refit at every later projection. The pose refit converges in two or
  /// three, the shape fit in about five.
  int iterations = 6;
  /// Levenberg pull toward the starting guess: the pose toward the vertex
  /// chart, the curvatures toward the store at the foot, `b` toward zero, the
  /// plain torus. Each row is scaled to the length it moves a point at the
  /// ring's reach, so this is relative to the squared gaps. Raise it when a
  /// noisy or flat ring makes `b` wander.
  real prior = 1e-6;
  /// Weight of the store's star of edge tensors against the ring points in the
  /// frame-start shape fit. At 1 a curvature mismatch counts like the height it
  /// bends a point at the reach. Higher trusts the store's shape over the
  /// mesh's; 0 fits the shape to the points alone.
  real fit_shape = 1;
  /// How much of the store's curvature change since the frame-start fit the
  /// cyclide follows during the step. 1 follows it fully, so shape edits on
  /// the store reach the mesh. 0 holds the fitted curvatures for the step, a
  /// pure smoother. `b` always holds for the step.
  real store_gain = 1;
  /// Ring-point rows of the fit's Jacobian from the sympy-generated gradient
  /// (`scripts/sympy/osculating_dupin_jacobian.py`). False takes central
  /// differences instead, to cross-check or compare cost.
  bool exact_gradient = true;
  /// Frame-start shape in closed form instead of the joint fit. The curvatures
  /// come from the store at the foot. The third-order pair comes from how the
  /// store's star of edge tensors changes across the ring: a least-squares
  /// slope of the low curvature along the tube direction, and of the high
  /// curvature across it, each divided by twice the curvature gap. That gap is
  /// regularized, so an umbilic gets a small `b` instead of an infinite one.
  /// About fifty times cheaper than the fit; `fit_shape` is unused when on.
  bool direct_shape = true;
  /// Use the cyclide's third-order pair `b`. False holds `b` at zero, so the
  /// row lands on the plain osculating torus, for comparison.
  bool third_order = true;
  /// Ridge on the third-order pair. A cyclide's cross cubic terms are the
  /// curvature gap times `b`, so near an umbilic a measured slope asks for a
  /// huge `b`, and `S_b` with a huge `b` is wildly nonlinear at the ring: those
  /// landings pulled a 2-ring row inward until the tube collapsed. The ridge
  /// solves `min (2 Δ b − slope)² + (2 ε b)²` with `ε = b_ridge (|k_hi| +
  /// |k_lo|)` plus a tiny flat-patch floor, so `b = Δ slope / (2 (Δ² + ε²))`: exact where the
  /// gap is large next to `ε`, and smoothly to zero at an umbilic, where no
  /// Dupin cyclide carries cross terms. The joint fit gets the same ridge in
  /// height units. 0 removes it.
  real b_ridge = 0.15;
  /// Backstop on `|b|` times the star's reach, past the ridge. The map `S_b`
  /// moves a point at distance `r` tangentially by about `|b| r^2`.
  real b_cap = 2.0;
  /// Gauss-Newton iterations for the pose refit at each projection, which
  /// starts from the previous projection's pose and stops early once its
  /// step is negligible. `iterations` still caps the joint shape fit.
  int pose_iterations = 3;
  /// Refit the pose at every this-many projections and reuse the last pose in
  /// between. With the warm start the refit usually stops at its first
  /// convergence check, so this rarely buys much; it is there for large meshes
  /// or many solver iterations. 1 refits at every projection.
  int pose_refit_every = 1;
  /// One-ring cone volume about the spine, the cyclide's tube centre at the
  /// foot. See `ring_volume_config`. With `land` false and no store rows the
  /// row skips both fits and is a pure volume constraint.
  ring_volume_config volume{};
};

/// Gap of `d` (from the unposed foot, in world coordinates) against the posed
/// cyclide, and its exact gradient in `(offset, tilt_t, tilt_b, k_hi, k_lo,
/// b_hi, b_lo)`, by chaining the sympy-generated stages. The point must land
/// (`land_on_posed_dupin`); the gap matches that function's.
inline void dupin_gap_and_gradient(const dupin_ring_shape &shape, const vec3 &T, const vec3 &B, const vec3 &N,
                                   const torus_ring_pose &pose, const vec3 &d, real &gap,
                                   Eigen::Matrix<real, 7, 1> &grad) {
  namespace g = dupin_generated;
  Eigen::Matrix<real, 7, 1> c;
  c << d.dot(T), d.dot(B), d.dot(N), shape.torus.hx, shape.torus.hy, shape.torus.lx, shape.torus.ly;
  // Pose: P and dP / d(offset, tilt_t, tilt_b).
  Eigen::Matrix<real, 3, 1> P;
  Eigen::Matrix<real, 3, 3> JP;
  g::dupin_pose(c, Eigen::Matrix<real, 3, 1>(pose.offset, pose.tilt_t, pose.tilt_b), P, JP);
  Eigen::Matrix<real, 3, 7> dP = Eigen::Matrix<real, 3, 7>::Zero();
  dP.leftCols<3>() = JP;
  // Pullback: q and its derivative through P, plus b directly.
  Eigen::Matrix<real, 5, 1> xq;
  xq << P, shape.b_hi, shape.b_lo;
  Eigen::Matrix<real, 2, 1> q;
  Eigen::Matrix<real, 2, 5> Jq;
  g::dupin_pullback(xq, q, Jq);
  Eigen::Matrix<real, 2, 7> dq = Jq.leftCols<3>() * dP;
  dq.rightCols<2>() += Jq.rightCols<2>();
  // Torus graph at q: (h, h_u, h_v), through q, plus the curvatures directly.
  Eigen::Matrix<real, 4, 1> xt(q[0], q[1], shape.torus.k_hi, shape.torus.k_lo);
  Eigen::Matrix<real, 3, 1> th;
  Eigen::Matrix<real, 3, 4> Jt;
  g::dupin_torus(xt, th, Jt);
  Eigen::Matrix<real, 3, 7> dt = Jt.leftCols<2>() * dq;
  dt.middleCols<2>(3) += Jt.rightCols<2>();
  // Gap of P against the pushed-forward torus point.
  Eigen::Matrix<real, 10, 1> xg;
  xg << P, shape.b_hi, shape.b_lo, q, th;
  Eigen::Matrix<real, 1, 1> y;
  Eigen::Matrix<real, 1, 10> Jg;
  g::dupin_gap(xg, y, Jg);
  gap = y[0];
  const Eigen::Matrix<real, 1, 7> row =
      Jg.leftCols<3>() * dP + Jg.middleCols<2>(5) * dq + Jg.rightCols<3>() * dt;
  grad = row.transpose();
  grad.tail<2>() += Jg.middleCols<2>(3).transpose();
}

/// The frame-start cyclide in closed form, from the store alone. `shape0`
/// carries the store's curvatures and principal directions at the foot. Each
/// star edge's stored tensor is read in the foot's principal frame at the
/// edge midpoint, and two weighted least-squares planes give the cross
/// derivatives: `c = d k_lo / d u` (tube direction) and `d = d k_hi / d v`.
/// A cyclide's cubic term over its torus, `(k_lo − k_hi) b_hi u v² +
/// (k_hi − k_lo) b_lo u² v`, has exactly those slopes, so
///
///   b_hi =  c / (2 Δ),   b_lo = −d / (2 Δ),   Δ = k_lo − k_hi,
///
/// with `1 / Δ` regularized as `Δ / (Δ² + ε²)`. The ridge slopes, which no
/// Dupin cyclide can hold, are left out. First order in `b`, which is small
/// at ring scale.
inline dupin_ring_shape direct_dupin_shape(const dupin_ring_shape &shape0, const vec3 &T, const vec3 &B,
                                           const vec3 &N, const std::vector<dupin_ring_edge> &edges,
                                           real b_ridge = real(0.15), real b_cap = real(2)) {
  dupin_ring_shape out = shape0;
  out.b_hi = 0;
  out.b_lo = 0;
  if (edges.size() < 3)
    return out;
  const vec3 eh = shape0.torus.hx * T + shape0.torus.hy * B;
  const vec3 el = shape0.torus.lx * T + shape0.torus.ly * B;
  Eigen::Matrix<real, 3, 3> A = Eigen::Matrix<real, 3, 3>::Zero();
  Eigen::Matrix<real, 3, 1> rl = Eigen::Matrix<real, 3, 1>::Zero(), rh = Eigen::Matrix<real, 3, 1>::Zero();
  real reach2 = 0;
  for (const dupin_ring_edge &e : edges) {
    vec3 t = e.hi - e.lo;
    t -= N * t.dot(N);
    if (!(t.squaredNorm() > real(1e-24)) || !(e.mass > real(0)) || !e.store.allFinite())
      continue;
    t.normalize();
    const vec3 bt = N.cross(t);
    const mat3 S = e.store[0] * t * t.transpose() + e.store[1] * (t * bt.transpose() + bt * t.transpose()) +
                   e.store[2] * bt * bt.transpose();
    const vec3 m = real(0.5) * (e.lo + e.hi);
    const Eigen::Matrix<real, 3, 1> row(real(1), m.dot(eh), m.dot(el));
    reach2 = std::max(reach2, m.squaredNorm());
    A += e.mass * row * row.transpose();
    rl += e.mass * row * el.dot(S * el);
    rh += e.mass * row * eh.dot(S * eh);
  }
  Eigen::LDLT<Eigen::Matrix<real, 3, 3>> ldlt(A);
  if (ldlt.info() != Eigen::Success || !(reach2 > real(0)))
    return out;
  const Eigen::Matrix<real, 3, 1> pl = ldlt.solve(rl), ph = ldlt.solve(rh);
  if (!pl.allFinite() || !ph.allFinite())
    return out;
  const real c = pl[1], d = ph[2];
  const real delta = shape0.torus.k_lo - shape0.torus.k_hi;
  const real k_scale = std::abs(shape0.torus.k_hi) + std::abs(shape0.torus.k_lo);
  const real eps = std::max(b_ridge, real(0)) * k_scale + real(1e-3) / std::sqrt(reach2);
  const real inv = delta / (delta * delta + eps * eps);
  out.b_hi = real(0.5) * c * inv;
  out.b_lo = real(-0.5) * d * inv;
  // `b` moves a point at the reach by about |b| r^3 |Δ|; keep that a small
  // correction even when the store's slopes are noisy.
  const real cap = std::max(b_cap, real(0)) / std::sqrt(reach2);
  out.b_hi = std::clamp(out.b_hi, -cap, cap);
  out.b_lo = std::clamp(out.b_lo, -cap, cap);
  return out;
}

struct dupin_ring_fit {
  torus_ring_pose pose;
  dupin_ring_shape shape;
  real cost0 = 0;
  real cost = 0;
  bool ok = false;
};

/// Least-squares fit on one ring. With `fit_shape` the unknowns are the pose,
/// both curvatures and `b`, and the store's edge tensors enter beside the
/// ring points. Otherwise only the pose is solved, against the points.
inline dupin_ring_fit fit_dupin_ring(const dupin_ring_shape &shape0, const vec3 &T, const vec3 &B,
                                     const vec3 &N, const std::vector<torus_ring_point> &pts,
                                     const std::vector<dupin_ring_edge> &edges, std::vector<char> &in,
                                     const osculating_dupin_ring_config &cfg, bool fit_shape,
                                     const torus_ring_pose *start = nullptr) {
  dupin_ring_fit fit;
  fit.shape = shape0;
  const size_t n = pts.size();
  in.assign(n, 0);
  real reach2 = 0;
  for (const torus_ring_point &pt : pts)
    reach2 = std::max(reach2, pt.d.squaredNorm());
  const real L = std::sqrt(reach2);
  if (!(L > real(1e-12)))
    return fit;
  const real sp = std::sqrt(std::max(cfg.prior, real(0)));
  const real sw = std::sqrt(std::max(cfg.fit_shape, real(0)));
  // Ridge on b in the height its cubic moves a point at the reach: the data
  // see b through the gap, `|Δ| L^3`; the ridge adds `2 ε L^3` that does not
  // vanish at an umbilic.
  const real k_scale0 = std::abs(shape0.torus.k_hi) + std::abs(shape0.torus.k_lo);
  const real eps_b = std::max(cfg.b_ridge, real(0)) * k_scale0 + real(1e-3) / L;
  const real b_prior =
      std::hypot(sp * std::abs(shape0.torus.k_hi - shape0.torus.k_lo), real(2) * eps_b) * L * L * L;
  const real step = real(1e-4) * L;
  const int np = fit_shape ? 7 : 3;
  real mass = 0;
  for (const dupin_ring_edge &e : edges)
    mass += e.mass;
  mass = edges.empty() ? real(1) : mass / real(edges.size());

  using param = Eigen::Matrix<real, 7, 1>;
  auto unpack = [&](const param &x, torus_ring_pose &pose, dupin_ring_shape &shape) {
    pose = {.offset = x[0], .tilt_t = x[1], .tilt_b = x[2]};
    shape = shape0;
    if (fit_shape) {
      shape.torus.k_hi = x[3];
      shape.torus.k_lo = x[4];
      shape.b_hi = x[5];
      shape.b_lo = x[6];
    }
  };
  auto residuals = [&](const param &x, std::vector<real> &r, bool select) -> bool {
    torus_ring_pose pose;
    dupin_ring_shape shape;
    unpack(x, pose, shape);
    const torus_ring_frame f = tilt_frame(T, B, N, pose);
    if (!f.ok)
      return false;
    r.clear();
    for (size_t i = 0; i < n; ++i) {
      if (!select && !in[i])
        continue;
      vec3 land;
      real gap = 0;
      const bool ok = land_on_posed_dupin(shape, f, pose.offset, pts[i].d, land, gap);
      if (select) {
        in[i] = ok ? 1 : 0;
        if (!ok)
          continue;
      } else if (!ok) {
        return false;
      }
      r.push_back(pts[i].s * gap);
    }
    if (fit_shape && sw > real(0)) {
      // Store tensors in the height they bend a point at the reach.
      const real hs = sw * real(0.5) * L * L;
      for (const dupin_ring_edge &e : edges) {
        vec3 reading;
        if (!dupin_edge_reading(shape, f, pose.offset, e, step, reading))
          return false;
        const real s = hs * std::sqrt(e.mass / mass);
        r.push_back(s * (reading[0] - e.store[0]));
        r.push_back(s * std::sqrt(real(2)) * (reading[1] - e.store[1]));
        r.push_back(s * (reading[2] - e.store[2]));
      }
    }
    r.push_back(sp * pose.offset);
    r.push_back(sp * L * pose.tilt_t);
    r.push_back(sp * L * pose.tilt_b);
    if (fit_shape) {
      r.push_back(sp * real(0.5) * L * L * (shape.torus.k_hi - shape0.torus.k_hi));
      r.push_back(sp * real(0.5) * L * L * (shape.torus.k_lo - shape0.torus.k_lo));
      r.push_back(b_prior * shape.b_hi);
      r.push_back(b_prior * shape.b_lo);
    }
    return true;
  };
  auto cost_of = [](const std::vector<real> &r) {
    real c = 0;
    for (const real v : r)
      c += v * v;
    return c;
  };

  param x = param::Zero();
  if (start) {
    x[0] = start->offset;
    x[1] = start->tilt_t;
    x[2] = start->tilt_b;
  }
  x[3] = shape0.torus.k_hi;
  x[4] = shape0.torus.k_lo;
  x[5] = shape0.b_hi;
  x[6] = shape0.b_lo;
  // Scratch reused across calls on this thread: the fit runs at every
  // projection of every vertex.
  thread_local std::vector<real> r;
  if (!residuals(x, r, true))
    return fit;
  int kept = 0;
  for (const char c : in)
    kept += c ? 1 : 0;
  if (kept < 3)
    return fit;
  fit.cost0 = cost_of(r);
  real cost = fit.cost0;
  param h;
  h << real(1e-6) * L, real(1e-6), real(1e-6), real(1e-6) / L, real(1e-6) / L, real(1e-6) / L, real(1e-6) / L;
  thread_local std::vector<real> rp, rm, rt;
  // Jacobian rows in residual order: ring points exactly from the generated
  // gradient, store edges by central differences (frame-start shape fit only),
  // and the prior rows, which are diagonal.
  auto jacobian = [&](const param &x, size_t m, matX &Jm) -> bool {
    torus_ring_pose pose;
    dupin_ring_shape shape;
    unpack(x, pose, shape);
    Jm.setZero(static_cast<Eigen::Index>(m), np);
    Eigen::Index row = 0;
    for (size_t i = 0; i < n; ++i) {
      if (!in[i])
        continue;
      real gap = 0;
      Eigen::Matrix<real, 7, 1> grad;
      if (cfg.exact_gradient) {
        dupin_gap_and_gradient(shape, T, B, N, pose, pts[i].d, gap, grad);
      } else {
        for (int k = 0; k < np; ++k) {
          param xp = x, xm = x;
          xp[k] += h[k];
          xm[k] -= h[k];
          torus_ring_pose pp, pm;
          dupin_ring_shape sp_, sm;
          unpack(xp, pp, sp_);
          unpack(xm, pm, sm);
          vec3 land;
          real gp = 0, gm = 0;
          if (!land_on_posed_dupin(sp_, tilt_frame(T, B, N, pp), pp.offset, pts[i].d, land, gp) ||
              !land_on_posed_dupin(sm, tilt_frame(T, B, N, pm), pm.offset, pts[i].d, land, gm))
            return false;
          grad[k] = (gp - gm) / (real(2) * h[k]);
        }
      }
      if (!grad.head(np).allFinite())
        return false;
      Jm.row(row++) = pts[i].s * grad.head(np).transpose();
    }
    const Eigen::Index n_edge_rows = fit_shape && sw > real(0) ? Eigen::Index(3 * edges.size()) : 0;
    if (n_edge_rows > 0) {
      auto edge_rows = [&](const param &xe, std::vector<real> &out) -> bool {
        torus_ring_pose pe;
        dupin_ring_shape se;
        unpack(xe, pe, se);
        const torus_ring_frame fe = tilt_frame(T, B, N, pe);
        if (!fe.ok)
          return false;
        out.clear();
        const real hs = sw * real(0.5) * L * L;
        for (const dupin_ring_edge &e : edges) {
          vec3 reading;
          if (!dupin_edge_reading(se, fe, pe.offset, e, step, reading))
            return false;
          const real s = hs * std::sqrt(e.mass / mass);
          out.push_back(s * reading[0]);
          out.push_back(s * std::sqrt(real(2)) * reading[1]);
          out.push_back(s * reading[2]);
        }
        return true;
      };
      for (int k = 0; k < np; ++k) {
        param xp = x, xm = x;
        xp[k] += h[k];
        xm[k] -= h[k];
        if (!edge_rows(xp, rp) || !edge_rows(xm, rm))
          return false;
        for (Eigen::Index e = 0; e < n_edge_rows; ++e)
          Jm(row + e, k) = (rp[static_cast<size_t>(e)] - rm[static_cast<size_t>(e)]) / (real(2) * h[k]);
      }
      row += n_edge_rows;
    }
    Jm(row++, 0) = sp;
    Jm(row++, 1) = sp * L;
    Jm(row++, 2) = sp * L;
    if (fit_shape) {
      Jm(row++, 3) = sp * real(0.5) * L * L;
      Jm(row++, 4) = sp * real(0.5) * L * L;
      Jm(row++, 5) = b_prior;
      Jm(row++, 6) = b_prior;
    }
    return row == static_cast<Eigen::Index>(m);
  };
  const int max_it = std::max(fit_shape ? cfg.iterations : cfg.pose_iterations, 0);
  for (int it = 0; it < max_it; ++it) {
    const size_t m = r.size();
    thread_local matX Jm;
    const bool ok = jacobian(x, m, Jm);
    if (!ok)
      break;
    Eigen::Map<const vecX> rv(r.data(), static_cast<Eigen::Index>(m));
    const vecX stepv = -(Jm.transpose() * Jm).ldlt().solve(Jm.transpose() * rv);
    if (!stepv.allFinite())
      break;
    // Converged: the pose would move a point at the reach by a hair.
    if (std::abs(stepv[0]) + L * (std::abs(stepv[1]) + std::abs(stepv[2])) < real(1e-10) * L &&
        (!fit_shape || stepv.tail(np - 3).norm() * L * L < real(1e-10) * L))
      break;
    bool accepted = false;
    real a = 1;
    for (int ls = 0; ls < 5; ++ls, a *= real(0.5)) {
      param trial = x;
      for (int k = 0; k < np; ++k)
        trial[k] += a * stepv[k];
      if (!residuals(trial, rt, false))
        continue;
      const real c = cost_of(rt);
      if (c < cost) {
        x = trial;
        r = rt;
        cost = c;
        accepted = true;
        break;
      }
    }
    if (!accepted)
      break;
  }
  unpack(x, fit.pose, fit.shape);
  fit.cost = cost;
  fit.ok = true;
  return fit;
}

/// Principal axes of a later store reading, matched to the frame-start ones.
/// The tube is the larger `|κ|` at each reading, so on a saddle with nearly
/// equal and opposite curvatures growth or noise can swap the roles mid-step.
/// The roles follow the reference's tube direction instead. Each direction's
/// sign is then aligned with the reference's, because the eigenvector sign is
/// arbitrary and `b` is stored along those directions.
inline monge_axes align_axes_to_reference(monge_axes ax, const torus_ring_shape &ref) {
  const real stay = ax.hx * ref.hx + ax.hy * ref.hy;
  const real swap = ax.lx * ref.hx + ax.ly * ref.hy;
  if (std::abs(swap) > std::abs(stay)) {
    std::swap(ax.k_hi, ax.k_lo);
    std::swap(ax.hx, ax.lx);
    std::swap(ax.hy, ax.ly);
  }
  if (ax.hx * ref.hx + ax.hy * ref.hy < real(0)) {
    ax.hx = -ax.hx;
    ax.hy = -ax.hy;
  }
  if (ax.lx * ref.lx + ax.ly * ref.ly < real(0)) {
    ax.lx = -ax.lx;
    ax.ly = -ax.ly;
  }
  return ax;
}

/// Process-wide tallies over `osculating_dupin_ring::project`, for probes. Off
/// by default so the solve pays nothing.
struct osculating_dupin_ring_stats {
  static inline std::atomic<bool> enabled{false};
  static inline std::atomic<long> calls{0};
  static inline std::atomic<long> no_frame{0};
  static inline std::atomic<long> no_fit{0};
  static inline std::atomic<long> points{0};
  static inline std::atomic<long> landed{0};
  static inline std::atomic<double> move{0};
  static inline std::atomic<double> reach{0};
  static void reset() {
    calls = 0;
    no_frame = 0;
    no_fit = 0;
    points = 0;
    landed = 0;
    move = 0;
    reach = 0;
  }
};

/// Osculating Dupin cyclide projector on one vertex's 1-ring and 2-ring, with
/// rows on both the positions and the edge store.
///
/// At the first projection of a step it fits the cyclide's shape (both
/// curvatures and `b`) and pose to the ring points and the store's star of
/// edge tensors together. Later projections fit only the pose; the curvatures
/// follow the store's change since that fit and `b` holds, so the cyclide is
/// carried with the surface through the step.
///
/// Position rows are differences `x_j − x_i` landed on the posed cyclide.
/// Shape rows pull each star edge's stored tensor toward the cyclide's tensor
/// at that edge, at the edge's dual mass, like the coupling rows.
class osculating_dupin_ring : public block_constraint {
public:
  DEFINE_CREATE_FUNC(osculating_dupin_ring)

  /// `ids` is the foot, then the 1-ring, then the 2-ring, with one row scale
  /// per ring point in `scales`. `blocks` is positions, then the edge store.
  /// `w` weighs the position rows and `w_shape` the store rows.
  osculating_dupin_ring(std::vector<index_t> ids, std::vector<real> scales, const real &w, real w_shape,
                        std::vector<sim_block::ptr> blocks, const std::vector<vertex_chart> *charts,
                        const kusama::weingarten::vertex_stars *stars, osculating_dupin_ring_config cfg,
                        size_t n_ring1 = 0, real *depth = nullptr)
      : block_constraint(std::move(ids), w, std::move(blocks)), _scales(std::move(scales)), _w_shape(w_shape),
        _charts(charts), _stars(stars), _cfg(cfg) {
    _volume.n_ring1 = std::min(n_ring1, _scales.size());
    _volume.depth = depth;
    if (_stars && !_ids.empty() && _ids[0] < static_cast<index_t>(_stars->at.size()))
      _star = _stars->at[static_cast<size_t>(_ids[0])];
    if (_blocks.size() < 2 || !_blocks[1])
      _star.clear();
    _rest_initialized = true;
  }

  virtual std::string name() { return "osculating_dupin_ring"; }

  virtual void project(const vecX &q, vecX &p) {
    if (_blocks.empty() || _ids.empty())
      return;
    const index_t vi = _ids[0];
    const vec3 xi = _blocks[0]->get_vec3(vi, q);
    const size_t n = _scales.size();
    auto &pts = _pts;
    auto &target = _target;
    pts.resize(n);
    target.resize(n);
    for (size_t k = 0; k < n; ++k) {
      pts[k] = {.d = _blocks[0]->get_vec3(_ids[k + 1], q) - xi, .s = _scales[k]};
      target[k] = pts[k].d;
    }
    auto &edges = _edges;
    auto &edge_target = _edge_target;
    edges.resize(_star.size());
    edge_target.resize(_star.size());
    for (size_t k = 0; k < _star.size(); ++k) {
      const auto &sp = _star[k];
      edges[k] = {.lo = _blocks[0]->get_vec3(sp.lo, q) - xi,
                  .hi = _blocks[0]->get_vec3(sp.hi, q) - xi,
                  .store = _blocks[1]->get_vec3(sp.edge, q),
                  .mass = sp.w};
      edge_target[k] = edges[k].store;
    }

    shape_w W;
    vec3 T, B, N;
    vec3 spine = vec3::Zero();
    bool have_spine = false;
    auto ring_d = [&](size_t k) -> const vec3 & { return pts[k].d; };
    // Volume-only rows need neither fit: the spine comes from the store.
    const bool fits = _cfg.volume.lands() || _w_shape > real(0);
    bool frame_ok = false, fit_ok = false;
    bool variational = false;
    vec3 variational_foot = vec3::Zero();
    bool corot_done = false;
    vec3 corot_foot = vec3::Zero();
    if (vertex_row_frame(_blocks, _charts, vi, q, W, T, B, N, _stars)) {
      W = shape_with_axes(W);
      frame_ok = W.ok;
      if (W.ok && !fits) {
        spine = _volume.spine(monge_axes_of(W).k_hi, 0, N, ring_d, _cfg.volume, xi);
        have_spine = true;
      }
      if (W.ok && fits) {
        const monge_axes ax =
            _reference ? align_axes_to_reference(monge_axes_of(W), _fit0.torus) : monge_axes_of(W);
        dupin_ring_shape shape{
            .torus = {.k_hi = ax.k_hi, .k_lo = ax.k_lo, .hx = ax.hx, .hy = ax.hy, .lx = ax.lx, .ly = ax.ly}};
        if (!_reference) {
          _reference = true;
          _store0_hi = shape.torus.k_hi;
          _store0_lo = shape.torus.k_lo;
          if (_cfg.direct_shape) {
            _fit0 = direct_dupin_shape(shape, T, B, N, edges, _cfg.b_ridge, _cfg.b_cap);
            _ref_ok = true;
          } else {
            std::vector<char> in0;
            const dupin_ring_fit ref = fit_dupin_ring(shape, T, B, N, pts, edges, in0, _cfg, true);
            _ref_ok = ref.ok;
            _fit0 = ref.shape;
          }
          if (!_cfg.third_order) {
            _fit0.b_hi = 0;
            _fit0.b_lo = 0;
          }
        }
        if (_ref_ok) {
          shape.torus.k_hi = _fit0.torus.k_hi + _cfg.store_gain * (shape.torus.k_hi - _store0_hi);
          shape.torus.k_lo = _fit0.torus.k_lo + _cfg.store_gain * (shape.torus.k_lo - _store0_lo);
          shape.b_hi = _fit0.b_hi;
          shape.b_lo = _fit0.b_lo;
        }
        std::vector<char> in;
        dupin_ring_fit fit;
        const bool refit = !_have_pose || ++_since_refit >= std::max(_cfg.pose_refit_every, 1);
        variational = _cfg.volume.mode == ring_volume_mode::fit && _cfg.volume.lands();
        if (variational) {
          // Pose and both curvatures from the variational fit; b holds.
          ring_fit_params theta;
          if (fit_variational(shape, T, B, N, xi, pts, theta)) {
            shape.torus.k_hi = theta[3];
            shape.torus.k_lo = theta[4];
            fit.pose = {.offset = theta[0], .tilt_t = theta[1], .tilt_b = theta[2]};
            fit.shape = shape;
            fit.ok = true;
            in.assign(n, 1);
            variational_foot = theta[0] * tilt_frame(T, B, N, fit.pose).N;
          }
        } else if (refit) {
          fit = fit_dupin_ring(shape, T, B, N, pts, edges, in, _cfg, false,
                               _have_pose ? &_last_fit.pose : nullptr);
          _since_refit = 0;
        } else {
          // Reuse the last pose; every point is a landing candidate.
          fit = _last_fit;
          fit.shape = shape;
          in.assign(n, 1);
        }
        const torus_ring_frame f = tilt_frame(T, B, N, fit.pose);
        if (fit.ok && f.ok) {
          fit_ok = true;
          _last_fit = fit;
          _have_pose = true;
          spine = _volume.spine(shape.torus.k_hi, fit.pose.offset, f.N, ring_d, _cfg.volume, xi);
          have_spine = true;
          for (size_t k = 0; _cfg.volume.lands() && k < n; ++k) {
            vec3 land;
            if (in[k] && land_point_on_posed_dupin(shape, f, fit.pose.offset, pts[k].d, land))
              target[k] = land;
          }
          if (_cfg.volume.mode == ring_volume_mode::corotational)
            corot_done = corotational_targets(shape, f, fit.pose.offset, pts, target, corot_foot);
          if (_w_shape > real(0)) {
            real reach2 = 0;
            for (const torus_ring_point &pt : pts)
              reach2 = std::max(reach2, pt.d.squaredNorm());
            const real step = real(1e-4) * std::sqrt(reach2);
            for (size_t k = 0; k < edges.size(); ++k) {
              vec3 reading;
              if (dupin_edge_reading(shape, f, fit.pose.offset, edges[k], step, reading))
                edge_target[k] = reading;
            }
          }
        }
      }
    }

    vec3 foot_shift = variational ? variational_foot : vec3::Zero();
    if (have_spine)
      foot_shift = _volume.hold(target, ring_d, _scales, spine, _cfg.volume);
    if (corot_done)
      foot_shift = corot_foot;
    if (osculating_dupin_ring_stats::enabled.load(std::memory_order_relaxed)) {
      using stats = osculating_dupin_ring_stats;
      stats::calls.fetch_add(1, std::memory_order_relaxed);
      if (!frame_ok)
        stats::no_frame.fetch_add(1, std::memory_order_relaxed);
      else if (!fit_ok)
        stats::no_fit.fetch_add(1, std::memory_order_relaxed);
      long landed = 0;
      double mv = 0, rc = 0;
      for (size_t k = 0; k < n; ++k) {
        const real m = (target[k] - pts[k].d).norm();
        landed += m > real(0) ? 1 : 0;
        mv += m;
        rc += pts[k].d.norm();
      }
      stats::points.fetch_add(static_cast<long>(n), std::memory_order_relaxed);
      stats::landed.fetch_add(landed, std::memory_order_relaxed);
      stats::move.fetch_add(mv, std::memory_order_relaxed);
      stats::reach.fetch_add(rc, std::memory_order_relaxed);
    }

    index_t row = _id0;
    for (size_t k = 0; k < n; ++k) {
      vec3 rowp = _w * _scales[k] * target[k];
      if (!rowp.allFinite())
        rowp = vec3::Zero();
      p.block(row, 0, 3, 1) = rowp;
      row += 3;
    }
    if (_volume.foot_row(_cfg.volume)) {
      vec3 rowp = _w * _cfg.volume.foot_weight * (xi + foot_shift);
      if (!rowp.allFinite())
        rowp = vec3::Zero();
      p.block(row, 0, 3, 1) = rowp;
      row += 3;
    }
    if (_w_shape > real(0)) {
      for (size_t k = 0; k < _star.size(); ++k) {
        vec3 rowp = _w_shape * _star[k].w * edge_target[k];
        if (!rowp.allFinite())
          rowp = vec3::Zero();
        p.block(row, 0, 3, 1) = rowp;
        row += 3;
      }
    }
  }

  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    if (_blocks.empty() || !_blocks[0] || _ids.empty())
      return;
    const index_t ii = _blocks[0]->get_offset_idx(_ids[0]);
    for (size_t k = 0; k < _scales.size(); ++k) {
      const index_t jj = _blocks[0]->get_offset_idx(_ids[k + 1]);
      const real s = _w * _scales[k];
      for (int ax = 0; ax < 3; ++ax) {
        triplets.push_back(trip(id0 + ax, ii + ax, -s));
        triplets.push_back(trip(id0 + ax, jj + ax, s));
      }
      id0 += 3;
    }
    if (_volume.foot_row(_cfg.volume)) {
      for (int ax = 0; ax < 3; ++ax)
        triplets.push_back(trip(id0 + ax, ii + ax, _w * _cfg.volume.foot_weight));
      id0 += 3;
    }
    if (_w_shape > real(0)) {
      for (const auto &sp : _star) {
        const index_t ee = _blocks[1]->get_offset_idx(sp.edge);
        for (int c = 0; c < 3; ++c)
          triplets.push_back(trip(id0 + c, ee + c, _w_shape * sp.w));
        id0 += 3;
      }
    }
  }

  const dupin_ring_fit &last_fit() const { return _last_fit; }
  const dupin_ring_shape &reference_shape() const { return _fit0; }
  /// Cone volume the volume term holds (NaN before the first projection).
  real reference_volume() const { return _volume.v0; }

  /// `corotational` mode: the skeleton from the posed cyclide's normals at the
  /// landings, then the corotational step on top of the row's landings.
  bool corotational_targets(const dupin_ring_shape &shape, const torus_ring_frame &f, real offset,
                            const std::vector<torus_ring_point> &pts, std::vector<vec3> &target, vec3 &foot) {
    const size_t n = pts.size(), n1 = _volume.n_ring1;
    if (n1 < 3)
      return false;
    thread_local std::vector<vec3> d, S, y;
    d.resize(n);
    real r2 = 0;
    for (size_t k = 0; k < n; ++k) {
      d[k] = pts[k].d;
      if (k < n1)
        r2 = std::max(r2, d[k].squaredNorm());
    }
    const real cap = _cfg.volume.apex_cap * std::sqrt(r2);
    const vec3 eh = shape.torus.hx * f.T + shape.torus.hy * f.B;
    const vec3 el = shape.torus.lx * f.T + shape.torus.ly * f.B;
    S.resize(n1);
    for (size_t k = 0; k < n1; ++k) {
      const vec3 p = target[k] - offset * f.N;
      vec3 l, nl;
      real gap = 0;
      vec3 normal = f.N;
      if (dupin_land_local(shape, vec3(p.dot(eh), p.dot(el), p.dot(f.N)), l, gap, &nl))
        normal = (nl[0] * eh + nl[1] * el + nl[2] * f.N).normalized();
      S[k] = medial_image(target[k], normal, shape.torus.k_hi, cap);
    }
    const vec3 SF = medial_image(offset * f.N, f.N, shape.torus.k_hi, cap);
    real *depth = _cfg.volume.corotational.persistent ? _volume.depth : nullptr;
    if (!corotational_ring_step(d, n1, target, S, SF, _volume.corot, _cfg.volume.corotational, depth, y))
      return false;
    for (size_t k = 0; k < n; ++k)
      target[k] = y[k + 1] - y[0];
    foot = y[0];
    return true;
  }

  /// `fit` volume mode: the variational fit over pose and both curvatures on
  /// the cyclide model, `b` held from the frame-start shape, with the store's
  /// curvatures in `shape` as the pull target.
  bool fit_variational(const dupin_ring_shape &shape, const vec3 &T, const vec3 &B, const vec3 &N,
                       const vec3 &xi, const std::vector<torus_ring_point> &pts, ring_fit_params &theta) {
    const size_t n = pts.size();
    thread_local std::vector<vec3> d;
    thread_local std::vector<real> sc;
    d.resize(n);
    sc.resize(n);
    for (size_t k = 0; k < n; ++k) {
      d[k] = pts[k].d;
      sc[k] = pts[k].s;
    }
    auto land = [&](const ring_fit_params &x, const vec3 &dk) {
      ring_model_point e;
      const torus_ring_frame fr = tilt_frame(T, B, N, {.offset = x[0], .tilt_t = x[1], .tilt_b = x[2]});
      if (!fr.ok)
        return e;
      dupin_ring_shape sh = shape;
      sh.torus.k_hi = x[3];
      sh.torus.k_lo = x[4];
      const vec3 eh = sh.torus.hx * fr.T + sh.torus.hy * fr.B;
      const vec3 el = sh.torus.lx * fr.T + sh.torus.ly * fr.B;
      const vec3 p = dk - x[0] * fr.N;
      vec3 l, nl;
      if (!dupin_land_local(sh, vec3(p.dot(eh), p.dot(el), p.dot(fr.N)), l, e.gap, &nl))
        return e;
      e.land = x[0] * fr.N + l[0] * eh + l[1] * el + l[2] * fr.N;
      e.normal = (nl[0] * eh + nl[1] * el + nl[2] * fr.N).normalized();
      e.ok = e.land.allFinite() && e.normal.allFinite();
      return e;
    };
    auto foot = [&](const ring_fit_params &x, vec3 &F, vec3 &NF) {
      const torus_ring_frame fr = tilt_frame(T, B, N, {.offset = x[0], .tilt_t = x[1], .tilt_b = x[2]});
      NF = fr.ok ? fr.N : N;
      F = x[0] * NF;
    };
    ring_fit_params start;
    start << 0, 0, 0, shape.torus.k_hi, shape.torus.k_lo;
    ring_variational_config vcfg{.w_store = _cfg.volume.store_weight,
                                 .lo_stiffness = _cfg.volume.lo_stiffness,
                                 .w_volume = _cfg.volume.weight,
                                 .prior = _cfg.prior,
                                 .iterations = _cfg.pose_iterations,
                                 .apex_cap = _cfg.volume.apex_cap};
    // A held skeleton gives the foot offset a meaning: it moves the patch
    // relative to the skeleton. Without one it is only a nuisance and is pinned.
    if (_cfg.volume.held_medial)
      vcfg.offset_weight = 0;
    const size_t n1 = _volume.n_ring1;
    real r2 = 0;
    for (size_t k = 0; k < n1 && k < n; ++k)
      r2 = std::max(r2, d[k].squaredNorm());
    const real cap = vcfg.apex_cap * std::sqrt(r2);
    // The step's tube-axis skeleton, held in space, re-expressed about the
    // moving foot.
    if (_cfg.volume.held_medial && !_volume.medial_set) {
      std::vector<vec3> M;
      vec3 MF;
      if (!ring_medial_images(start, d, n1, land, foot, cap, M, MF))
        return false;
      _volume.medial_world.resize(n1);
      for (size_t j = 0; j < n1; ++j)
        _volume.medial_world[j] = xi + M[j];
      _volume.medial_foot_world = xi + MF;
      _volume.medial_set = true;
    }
    ring_medial_hold hold;
    if (_volume.medial_set) {
      hold.M.resize(n1);
      for (size_t j = 0; j < n1; ++j)
        hold.M[j] = _volume.medial_world[j] - xi;
      hold.MF = _volume.medial_foot_world - xi;
    }
    const ring_medial_hold *held = _volume.medial_set ? &hold : nullptr;
    if (!std::isfinite(_volume.v0)) {
      real v = 0, a = 0;
      if (!ring_medial_volume(start, d, n1, land, foot, cap, v, a, held))
        return false;
      _volume.v0 = v;
      if (_cfg.volume.persistent && _volume.depth) {
        if (std::isfinite(*_volume.depth))
          _volume.v0 = *_volume.depth * a;
        else
          *_volume.depth = v / a;
      }
    }
    if (_have_theta)
      start = _theta;
    const ring_variational_fit fit = fit_ring_variational(start, shape.torus.k_hi, shape.torus.k_lo, d, sc,
                                                          n1, _volume.v0, vcfg, land, foot, held);
    if (!fit.ok)
      return false;
    _theta = fit.theta;
    _have_theta = true;
    theta = fit.theta;
    return true;
  }

private:
  std::vector<real> _scales;
  real _w_shape = 0;
  const std::vector<vertex_chart> *_charts = nullptr;
  const kusama::weingarten::vertex_stars *_stars = nullptr;
  std::vector<kusama::weingarten::vertex_stars::spoke> _star;
  osculating_dupin_ring_config _cfg;
  bool _reference = false;
  bool _ref_ok = false;
  real _store0_hi = 0;
  real _store0_lo = 0;
  dupin_ring_shape _fit0;
  dupin_ring_fit _last_fit;
  /// The pose refit starts from the previous projection's pose once there is one.
  bool _have_pose = false;
  /// `fit` volume mode: the previous projection's fitted parameters.
  bool _have_theta = false;
  ring_fit_params _theta = ring_fit_params::Zero();
  /// Per-projection scratch, reused so a projection does not allocate.
  std::vector<torus_ring_point> _pts;
  std::vector<vec3> _target;
  std::vector<dupin_ring_edge> _edges;
  std::vector<vec3> _edge_target;
  int _since_refit = 0;
  ring_volume_term _volume;
};

/// One osculating-cyclide row per vertex. `blocks` is positions, then the edge
/// store that `stars` reads. `w` weighs the position rows and `w_shape` the
/// store rows.
inline void init_osculating_dupin_ring(const asawa::shell::shell &shell,
                                       std::vector<projection_constraint::ptr> &constraints, const real &w,
                                       real w_shape, std::vector<sim_block::ptr> blocks,
                                       const std::vector<vertex_chart> *charts,
                                       const kusama::weingarten::vertex_stars *stars,
                                       osculating_dupin_ring_config cfg = {},
                                       std::vector<real> *depth = nullptr) {
  if (!(w > real(0) || w_shape > real(0)) || blocks.empty() || !blocks[0] || !charts)
    return;
  if (depth && depth->size() != static_cast<size_t>(shell.vert_count()))
    depth = nullptr;
  if (blocks.size() < 2 || !blocks[1])
    stars = nullptr;
  std::vector<index_t> ring1, ring2;
  for (int vi = 0; vi < shell.vert_count(); ++vi) {
    torus_ring_ids(shell, vi, ring1, ring2);
    if (ring1.size() < 3)
      continue;
    std::vector<index_t> ids{vi};
    std::vector<real> scales;
    ids.reserve(1 + ring1.size() + ring2.size());
    for (const index_t j : ring1) {
      ids.push_back(j);
      scales.push_back(real(1));
    }
    if (cfg.two_ring > real(0)) {
      for (const index_t j : ring2) {
        ids.push_back(j);
        scales.push_back(cfg.two_ring);
      }
    }
    real *entry = depth ? &(*depth)[static_cast<size_t>(vi)] : nullptr;
    constraints.push_back(osculating_dupin_ring::create(std::move(ids), std::move(scales), w, w_shape, blocks,
                                                        charts, stars, cfg, ring1.size(), entry));
  }
}

} // namespace block
} // namespace hepworth
} // namespace gaudi
