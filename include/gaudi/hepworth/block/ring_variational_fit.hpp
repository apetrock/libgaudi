#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

#include "gaudi/common.h"

namespace gaudi {
namespace hepworth {
namespace block {

/// One ring point evaluated on a posed local model: its landing, the model's
/// unit normal there, and its signed gap along the model.
struct ring_model_point {
  bool ok = false;
  vec3 land = vec3::Zero();
  vec3 normal = vec3::UnitZ();
  real gap = 0;
};

/// Fit parameters shared by the torus and cyclide rows: the pose (normal
/// offset of the foot and tilt of its normal) and both curvatures.
using ring_fit_params = Eigen::Matrix<real, 5, 1>;

/// Signed volume of the tetrahedron `(a, b, c, d)`.
inline real tet_volume(const vec3 &a, const vec3 &b, const vec3 &c, const vec3 &d) {
  return (b - a).dot((c - a).cross(d - a)) / real(6);
}

/// Volume between the fan `(F, L_j, L_j+1)` and its medial images `(M_F, M_j,
/// M_j+1)`: one prism per fan triangle, three tetrahedra each, for the first
/// `n` points in cyclic order. Each medial image is the centre of that point's
/// larger curvature sphere, so the region is a cone where the patch is
/// umbilic (the images coincide), a wedge on a tube (they lie on the axis),
/// and between the two in general: the slice of solid between the patch and
/// its local medial set.
inline real medial_prism_volume(const vec3 &F, const vec3 &MF, const std::vector<vec3> &L,
                                const std::vector<vec3> &M, size_t n) {
  real v = 0;
  for (size_t j = 0; j < n; ++j) {
    const vec3 &a = L[j], &b = L[(j + 1) % n];
    const vec3 &ma = M[j], &mb = M[(j + 1) % n];
    v += tet_volume(F, a, b, MF) + tet_volume(a, b, MF, ma) + tet_volume(b, MF, ma, mb);
  }
  return v;
}

/// Area of the fan `(F, L_j, L_j+1)`.
inline real fan_area_about(const vec3 &F, const std::vector<vec3> &L, size_t n) {
  real a = 0;
  for (size_t j = 0; j < n; ++j)
    a += real(0.5) * (L[j] - F).cross(L[(j + 1) % n] - F).norm();
  return a;
}

/// Medial image of a model point: `p + n / k_hi`, the centre of the sphere of
/// its larger (tube) curvature, in the convention where a height `½ k t²`
/// bends toward `+n`. On a tube or torus these points lie on the tube's axis,
/// so a ring's images line up along it and the region they bound is a
/// cylindrical wedge; only at an umbilic do they coincide in one point. The distance is capped at `cap` so a flat patch stays finite.
inline vec3 medial_image(const vec3 &p, const vec3 &normal, real k_hi, real cap) {
  real a = std::abs(k_hi) > real(0) ? real(1) / std::abs(k_hi) : cap;
  if (!(a < cap) || !std::isfinite(a))
    a = cap;
  return p + (k_hi < real(0) ? -a : a) * normal;
}

struct ring_variational_config {
  /// Pull of both curvatures toward the store's, in the height a mismatch
  /// bends a point at the ring's reach, per one-ring point. High keeps the
  /// store's shape; the volume term then has to win against it to change it.
  real w_store = 10;
  /// Extra stiffness of the low curvature's store pull, as a multiple of
  /// `w_store`. A ring's medial volume depends on both curvatures, and the
  /// low one is the cheaper lever: bending the patch around the major circle
  /// changes how its curvature-sphere centres spread. Left equal, the fit
  /// meets the volume that way while the tube radius follows the store.
  /// Stiffer here makes the volume act through the tube curvature.
  real lo_stiffness = 10;
  /// Weight of the volume residual, `(V − V₀) / A`, a mean normal height,
  /// per one-ring point. Around 1 it trades evenly with the store pull per
  /// point; 1e4 makes it effectively hard.
  real w_volume = 1;
  /// Pose pull toward the chart, relative to the squared gaps.
  real prior = 1e-6;
  /// Pull of the model's foot onto the vertex, in units of a one-ring gap. On
  /// a one-ring every point is about the same distance from the foot, so a
  /// free normal offset absorbs any mean-curvature mismatch: the mean ring
  /// height is `offset + ½ k L²`. A ring's mean curvature lives in how high
  /// the centre sits on its ring, so the model has to pass through its own
  /// vertex for the curvatures, and so the volume, to reach the mesh. Large
  /// pins it there; 0 frees it (the pose-only fit's behaviour).
  real offset_weight = 1e4;
  /// Gauss-Newton iterations; the fit starts from the previous projection's
  /// parameters and stops early once its step is negligible.
  int iterations = 4;
  /// Farthest medial image, in multiples of the ring's reach.
  real apex_cap = 10;
};

/// A medial set held in space for the step: each one-ring point's and the
/// foot's curvature-sphere centre at the step's first projection, given here
/// relative to the current unposed foot. With it the fit's volume is the
/// region between the moving patch and a fixed skeleton, a real volume: a
/// patch moving toward its skeleton loses it whether or not it curves. Without
/// it the images follow the model, and the volume only measures the local
/// curvature radius.
struct ring_medial_hold {
  std::vector<vec3> M;
  vec3 MF = vec3::Zero();
};

struct ring_variational_fit {
  ring_fit_params theta = ring_fit_params::Zero();
  real volume = 0;
  real area = 0;
  bool ok = false;
};

/// The ring's medial volume and fan area for parameters `theta`. `land(theta,
/// d)` evaluates a ring point on the posed model; `foot(theta, F, NF)` gives
/// the posed foot and its normal. Only the first `n1` points (the cyclic
/// one-ring) enter the volume.
template <class Land, class Foot>
bool ring_medial_volume(const ring_fit_params &theta, const std::vector<vec3> &d, size_t n1, Land &&land,
                        Foot &&foot, real cap, real &volume, real &area,
                        const ring_medial_hold *hold = nullptr) {
  if (n1 < 3)
    return false;
  thread_local std::vector<vec3> L, M;
  L.resize(n1);
  M.resize(n1);
  for (size_t j = 0; j < n1; ++j) {
    const ring_model_point e = land(theta, d[j]);
    if (!e.ok)
      return false;
    L[j] = e.land;
    M[j] = hold ? hold->M[j] : medial_image(e.land, e.normal, theta[3], cap);
  }
  vec3 F, NF;
  foot(theta, F, NF);
  const vec3 MF = hold ? hold->MF : medial_image(F, NF, theta[3], cap);
  volume = medial_prism_volume(F, MF, L, M, n1);
  area = fan_area_about(F, L, n1);
  return std::isfinite(volume) && area > real(0);
}

/// Curvature-sphere centres of the first `n1` ring points and the foot on the
/// model at `theta`, relative to the unposed foot.
template <class Land, class Foot>
bool ring_medial_images(const ring_fit_params &theta, const std::vector<vec3> &d, size_t n1, Land &&land,
                        Foot &&foot, real cap, std::vector<vec3> &M, vec3 &MF) {
  M.resize(n1);
  for (size_t j = 0; j < n1; ++j) {
    const ring_model_point e = land(theta, d[j]);
    if (!e.ok)
      return false;
    M[j] = medial_image(e.land, e.normal, theta[3], cap);
  }
  vec3 F, NF;
  foot(theta, F, NF);
  MF = medial_image(F, NF, theta[3], cap);
  return MF.allFinite();
}

/// Variational fit of one ring: the minimizer over `theta` of
///
///   Σ_j s_j² gap_j²                                      (ring points on the model)
///   + w_store n1 (½ L²)² |k(theta) − k_store|²           (store pull, both curvatures)
///   + w_volume n1 ((V(theta) − v0) / A)²                 (one-ring medial volume)
///   + prior (offset² + L² tilt²)
///
/// by Gauss-Newton with central-difference Jacobians from `start`. With `v0`
/// not finite the volume term is off. `d` and `s` are the ring points from
/// the unposed foot and their row scales.
template <class Land, class Foot>
ring_variational_fit fit_ring_variational(const ring_fit_params &start, real k_store_hi, real k_store_lo,
                                          const std::vector<vec3> &d, const std::vector<real> &s, size_t n1,
                                          real v0, const ring_variational_config &cfg, Land &&land,
                                          Foot &&foot, const ring_medial_hold *hold = nullptr) {
  ring_variational_fit out;
  const size_t n = d.size();
  real reach2 = 0;
  for (size_t j = 0; j < n1 && j < n; ++j)
    reach2 = std::max(reach2, d[j].squaredNorm());
  const real L = std::sqrt(reach2);
  if (!(L > real(1e-12)) || n1 < 3)
    return out;
  const real cap = cfg.apex_cap * L;
  const bool with_volume = std::isfinite(v0) && cfg.w_volume > real(0);
  const real sk = std::sqrt(std::max(cfg.w_store, real(0)) * real(n1)) * real(0.5) * L * L;
  const real sk_lo = sk * std::sqrt(std::max(cfg.lo_stiffness, real(0)));
  const real sv = std::sqrt(std::max(cfg.w_volume, real(0)) * real(n1));
  const real sp = std::sqrt(std::max(cfg.prior, real(0)));

  // Points that land at the start stay the residual set for the whole fit.
  std::vector<char> in(n, 0);
  for (size_t j = 0; j < n; ++j)
    in[j] = land(start, d[j]).ok ? 1 : 0;

  auto residuals = [&](const ring_fit_params &x, std::vector<real> &r) -> bool {
    r.clear();
    for (size_t j = 0; j < n; ++j) {
      if (!in[j])
        continue;
      const ring_model_point e = land(x, d[j]);
      if (!e.ok)
        return false;
      r.push_back(s[j] * e.gap);
    }
    r.push_back(sk * (x[3] - k_store_hi));
    r.push_back(sk_lo * (x[4] - k_store_lo));
    if (with_volume) {
      real v = 0, a = 0;
      if (!ring_medial_volume(x, d, n1, land, foot, cap, v, a, hold))
        return false;
      r.push_back(sv * (v - v0) / a);
    }
    r.push_back(std::sqrt(std::max(cfg.offset_weight, real(0)) * real(n1)) * x[0] + sp * x[0]);
    r.push_back(sp * L * x[1]);
    r.push_back(sp * L * x[2]);
    return true;
  };
  auto cost_of = [](const std::vector<real> &r) {
    real c = 0;
    for (const real v : r)
      c += v * v;
    return c;
  };

  ring_fit_params x = start;
  thread_local std::vector<real> r, rp, rm, rt;
  if (!residuals(x, r))
    return out;
  real cost = cost_of(r);
  ring_fit_params h;
  h << real(1e-7) * L, real(1e-7), real(1e-7), real(1e-7) / L, real(1e-7) / L;
  for (int it = 0; it < std::max(cfg.iterations, 0); ++it) {
    const Eigen::Index m = static_cast<Eigen::Index>(r.size());
    Eigen::Matrix<real, Eigen::Dynamic, 5> J(m, 5);
    bool ok = true;
    for (int k = 0; k < 5 && ok; ++k) {
      ring_fit_params xp = x, xm = x;
      xp[k] += h[k];
      xm[k] -= h[k];
      ok = residuals(xp, rp) && residuals(xm, rm) && static_cast<Eigen::Index>(rp.size()) == m &&
           static_cast<Eigen::Index>(rm.size()) == m;
      for (Eigen::Index i = 0; ok && i < m; ++i)
        J(i, k) = (rp[static_cast<size_t>(i)] - rm[static_cast<size_t>(i)]) / (real(2) * h[k]);
    }
    if (!ok)
      break;
    Eigen::Map<const vecX> rv(r.data(), m);
    const ring_fit_params step = -(J.transpose() * J).ldlt().solve(J.transpose() * rv);
    if (!step.allFinite())
      break;
    if (std::abs(step[0]) + L * (std::abs(step[1]) + std::abs(step[2])) +
            real(0.5) * L * L * (std::abs(step[3]) + std::abs(step[4])) <
        real(1e-10) * L)
      break;
    bool accepted = false;
    real a = 1;
    for (int ls = 0; ls < 5; ++ls, a *= real(0.5)) {
      const ring_fit_params trial = x + a * step;
      if (!residuals(trial, rt))
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
  out.theta = x;
  out.ok = ring_medial_volume(x, d, n1, land, foot, cap, out.volume, out.area, hold);
  return out;
}

/// Gradients of `tet_volume(a, b, c, d)` in each vertex.
inline void tet_volume_grad(const vec3 &a, const vec3 &b, const vec3 &c, const vec3 &d, vec3 &ga, vec3 &gb,
                            vec3 &gc, vec3 &gd) {
  gb = (c - a).cross(d - a) / real(6);
  gc = (d - a).cross(b - a) / real(6);
  gd = (b - a).cross(c - a) / real(6);
  ga = -(gb + gc + gd);
}

/// Rest state of one ring for the corotational volume mode, recorded at the
/// step's first projection, relative to the foot then.
struct corotational_rest {
  bool set = false;
  /// Ring points (all of them) and the foot, the foot at the origin.
  std::vector<vec3> P0;
  /// Tube-axis points of the one-ring and of the foot.
  std::vector<vec3> S0;
  vec3 SF0 = vec3::Zero();
  /// Rest volume of each fan prism.
  std::vector<real> V0;
  real A0 = 0;
};

/// Weighted Procrustes: the rotation `R` and translation `t` taking `src` to
/// `dst` best in `sum_i w_i |R src_i + t − dst_i|²`, as a proper rotation.
inline void weighted_procrustes(const std::vector<vec3> &src, const std::vector<vec3> &dst,
                                const std::vector<real> &w, mat3 &R, vec3 &t) {
  real W = 0;
  vec3 cs = vec3::Zero(), cd = vec3::Zero();
  for (size_t i = 0; i < src.size(); ++i) {
    W += w[i];
    cs += w[i] * src[i];
    cd += w[i] * dst[i];
  }
  R = mat3::Identity();
  t = vec3::Zero();
  if (!(W > real(0)))
    return;
  cs /= W;
  cd /= W;
  mat3 H = mat3::Zero();
  for (size_t i = 0; i < src.size(); ++i)
    H += w[i] * (src[i] - cs) * (dst[i] - cd).transpose();
  Eigen::JacobiSVD<mat3> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
  mat3 D = mat3::Identity();
  if ((svd.matrixV() * svd.matrixU().transpose()).determinant() < real(0))
    D(2, 2) = real(-1);
  R = svd.matrixV() * D * svd.matrixU().transpose();
  t = cd - R * cs;
}

/// Weights of the corotational volume mode's local solve. With `w_volume` 0
/// and `w_rest` and `w_stay` 0 the targets are the row's own landings, so the
/// mode is the bare row exactly; the volume and rest springs are added to it.
struct corotational_config {
  /// Spring of each point onto the row's own landing on its torus or cyclide.
  real w_proj = 1;
  /// Spring of each point to its corotated rest position: the ARAP / Müller
  /// shape-matching goal. 0 leaves the shape to the row's landing.
  real w_rest = 0;
  /// Spring of each point to where it is now.
  real w_stay = 0;
  /// Spring of each fan prism's volume to its rest volume, in the mean
  /// normal height its error stands for.
  real w_volume = 1;
  /// Weight of the tube-axis points in the Procrustes fit, relative to a ring
  /// point. They anchor the rotation out of the ring's nearly flat plane.
  real w_spine = 1;
  /// Gauss-Newton steps of the volume solve.
  int iterations = 3;
  /// Stiffness of the foot in the volume solve, relative to a ring point.
  /// The foot belongs to every prism of its ring, so with equal stiffness the
  /// smallest correction moves it most; its neighbours, whose rings hold it
  /// as an ordinary point, then push it back, and the mesh checkerboards. As
  /// stiff as the ring is wide spreads the correction over the ring instead.
  /// Negative uses the ring's size.
  real foot_stiffness = -1;
  /// Scale the rest volumes by the vertex's persistent medial depth (volume
  /// per fan area) instead of taking them as measured at the step's start.
  /// Off by default: the depth is a local tube radius, so holding it across
  /// steps holds each vertex's radius at its first value, which pushes a
  /// growing tube inward. The per-step rest holds the volume through each
  /// step and lets the radius follow the store between steps.
  bool persistent = false;
};

/// Corotational volume step of one ring, on top of the row's own landings.
///
/// The rest pose is the foot, the ring points and the model's tube-axis
/// points at the step's first projection, with one prism per fan triangle
/// between the patch and that skeleton. Each projection:
///
/// 1. a weighted Procrustes fit takes the rest points and rest skeleton onto
///    the current points and the current model's skeleton (`R`, `t`): the
///    corotated rest goals `g` and the corotated skeleton;
/// 2. the targets solve, by Gauss-Newton,
///      min Σ_k (w_proj |y_k − ℓ_k|² + w_rest |y_k − g_k|² + w_stay |y_k − x_k|²)
///          + w_volume Σ_j ((V_j(y) − V0_j) / a)²
///    with `ℓ_k` the row's landing (the foot's is where it is) and `a` the
///    mean rest prism area, against the corotated skeleton held fixed.
///
/// Positions are relative to the current foot; `d` and `l` are the ring
/// points and their landings, `S` and `SF` the current one-ring and foot
/// tube-axis points. Returns targets for the foot (`y[0]`) and the ring
/// points. `depth`, when given, makes the rest volumes persistent: their sum
/// is `*depth` times the rest fan area.
inline bool corotational_ring_step(const std::vector<vec3> &d, size_t n1, const std::vector<vec3> &l,
                                   const std::vector<vec3> &S, const vec3 &SF, corotational_rest &rest,
                                   const corotational_config &cc, real *depth, std::vector<vec3> &y) {
  const size_t n = d.size();
  if (n1 < 3 || n1 > n || l.size() != n || S.size() < n1)
    return false;
  if (!rest.set) {
    rest.P0.assign(n + 1, vec3::Zero());
    for (size_t k = 0; k < n; ++k)
      rest.P0[k + 1] = d[k];
    rest.S0.assign(S.begin(), S.begin() + static_cast<long>(n1));
    rest.SF0 = SF;
    rest.V0.resize(n1);
    real total = 0;
    for (size_t j = 0; j < n1; ++j) {
      const size_t k = (j + 1) % n1;
      const vec3 &A = rest.P0[j + 1], &B = rest.P0[k + 1];
      rest.V0[j] = tet_volume(vec3::Zero(), A, B, rest.SF0) + tet_volume(A, B, rest.SF0, rest.S0[j]) +
                   tet_volume(B, rest.SF0, rest.S0[j], rest.S0[k]);
      total += rest.V0[j];
    }
    std::vector<vec3> ring(rest.P0.begin() + 1, rest.P0.begin() + 1 + static_cast<long>(n1));
    rest.A0 = fan_area_about(vec3::Zero(), ring, n1);
    if (depth && rest.A0 > real(0) && std::abs(total) > real(0)) {
      if (std::isfinite(*depth)) {
        const real scale = *depth * rest.A0 / total;
        for (real &v : rest.V0)
          v *= scale;
      } else {
        *depth = total / rest.A0;
      }
    }
    rest.set = true;
  }
  // 1. Corotate the rest pose and skeleton onto the current ones.
  std::vector<vec3> src, dst;
  std::vector<real> w;
  src.push_back(vec3::Zero());
  dst.push_back(vec3::Zero());
  w.push_back(1);
  for (size_t k = 0; k < n1; ++k) {
    src.push_back(rest.P0[k + 1]);
    dst.push_back(d[k]);
    w.push_back(1);
  }
  src.push_back(rest.SF0);
  dst.push_back(SF);
  w.push_back(cc.w_spine);
  for (size_t k = 0; k < n1; ++k) {
    src.push_back(rest.S0[k]);
    dst.push_back(S[k]);
    w.push_back(cc.w_spine);
  }
  mat3 R;
  vec3 t;
  weighted_procrustes(src, dst, w, R, t);
  std::vector<vec3> Sr(n1);
  for (size_t k = 0; k < n1; ++k)
    Sr[k] = R * rest.S0[k] + t;
  const vec3 SFr = R * rest.SF0 + t;

  // 2. The springs' blend, then the volume springs by Gauss-Newton.
  const real wsum = cc.w_proj + cc.w_rest + cc.w_stay;
  if (!(wsum > real(0)))
    return false;
  std::vector<vec3> ystar(n + 1);
  ystar[0] = cc.w_rest * t / wsum; // the foot's landing and current position are the origin
  for (size_t k = 0; k < n; ++k)
    ystar[k + 1] = (cc.w_proj * l[k] + cc.w_rest * (R * rest.P0[k + 1] + t) + cc.w_stay * d[k]) / wsum;
  y = ystar;
  const real a = rest.A0 / real(n1);
  const real sv = std::sqrt(std::max(cc.w_volume, real(0))) / std::max(a, real(1e-30));
  if (!(sv > real(0)))
    return true;
  real reach2 = 0;
  for (size_t k = 0; k < n1; ++k)
    reach2 = std::max(reach2, d[k].squaredNorm());
  const Eigen::Index m = static_cast<Eigen::Index>(3 * (n1 + 1));
  const real foot_k = cc.foot_stiffness < real(0) ? real(n1) : std::max(cc.foot_stiffness, real(1e-6));
  for (int it = 0; it < std::max(cc.iterations, 0); ++it) {
    // Normal equations of  Σ_i m_i |Y_i − Y*_i|² + Σ_j (sv (V_j − V0_j))²,
    // m_i = wsum for a ring point and the foot's stiffness times it.
    matX H = matX::Identity(m, m) * wsum;
    H.block<3, 3>(0, 0) *= foot_k;
    vecX rhs = vecX::Zero(m);
    for (size_t i = 0; i <= n1; ++i)
      rhs.segment<3>(static_cast<Eigen::Index>(3 * i)) =
          (i == 0 ? foot_k : real(1)) * wsum * (ystar[i] - y[i]);
    for (size_t j = 0; j < n1; ++j) {
      const size_t k = (j + 1) % n1;
      const vec3 &F = y[0], &A = y[j + 1], &B = y[k + 1];
      const real v =
          tet_volume(F, A, B, SFr) + tet_volume(A, B, SFr, Sr[j]) + tet_volume(B, SFr, Sr[j], Sr[k]);
      vec3 gF = vec3::Zero(), gA = vec3::Zero(), gB = vec3::Zero(), g1, g2, g3, g4;
      tet_volume_grad(F, A, B, SFr, g1, g2, g3, g4);
      gF += g1;
      gA += g2;
      gB += g3;
      tet_volume_grad(A, B, SFr, Sr[j], g1, g2, g3, g4);
      gA += g1;
      gB += g2;
      tet_volume_grad(B, SFr, Sr[j], Sr[k], g1, g2, g3, g4);
      gB += g1;
      const std::array<std::pair<size_t, vec3>, 3> rows{{{0, sv * gF}, {j + 1, sv * gA}, {k + 1, sv * gB}}};
      const real res = sv * (v - rest.V0[j]);
      for (const auto &[ia, ja] : rows) {
        rhs.segment<3>(static_cast<Eigen::Index>(3 * ia)) -= ja * res;
        for (const auto &[ib, jb] : rows)
          H.block<3, 3>(static_cast<Eigen::Index>(3 * ia), static_cast<Eigen::Index>(3 * ib)) +=
              ja * jb.transpose();
      }
    }
    const vecX step = H.ldlt().solve(rhs);
    if (!step.allFinite())
      break;
    for (size_t i = 0; i <= n1; ++i)
      y[i] += step.segment<3>(static_cast<Eigen::Index>(3 * i));
    if (step.norm() < real(1e-12) * std::sqrt(reach2))
      break;
  }
  return true;
}

} // namespace block
} // namespace hepworth
} // namespace gaudi
