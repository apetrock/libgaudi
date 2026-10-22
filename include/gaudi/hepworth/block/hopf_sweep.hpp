#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "shape_objectives.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

/// Jacobi cn, sn, dn. m = 0 is cos, m = 1 is sech.
/// The Hopf profile does not use this at m = 0: cn(s, 0) = cos(s), and the
/// catalogue circle is constant curvature.
inline void jacobi_sn_cn_dn(real uu, real m, real &sn, real &cn, real &dn) {
  if (!(m > real(0))) {
    sn = std::sin(uu);
    cn = std::cos(uu);
    dn = 1;
    return;
  }
  if (!(m < real(1))) {
    sn = std::tanh(uu);
    cn = real(1) / std::cosh(uu);
    dn = cn;
    return;
  }
  real emc = real(1) - m;
  real a = 1;
  real c = 0;
  real em[14];
  real en[14];
  int l = 0;
  dn = 1;
  for (int i = 1; i <= 13; ++i) {
    l = i;
    em[i] = a;
    en[i] = (emc = std::sqrt(std::max(emc, real(0))));
    c = real(0.5) * (a + emc);
    if (std::abs(a - emc) <= real(1e-15) * a)
      break;
    emc *= a;
    a = c;
  }
  const real u = uu * c;
  sn = std::sin(u);
  cn = std::cos(u);
  if (sn != 0) {
    a = cn / sn;
    c *= a;
    for (int ii = l; ii >= 1; --ii) {
      const real b = em[ii];
      a *= c;
      c *= dn;
      dn = (en[ii] + a) / (b + a);
      a = c / b;
    }
    a = real(1) / std::sqrt(c * c + real(1));
    sn = sn >= 0 ? a : -a;
    cn = c * sn;
  }
}

/// Twist shifts m off the circular section. Nonpositive twist stays at 0.
inline real hopf_modulus(real twist) {
  if (!std::isfinite(twist) || !(twist > real(0)))
    return 0;
  const real cap = real(1) - real(1e-6);
  return twist < cap ? twist : cap;
}

/// Circle of curvature kappa0 in the (tube, normal) plane. Arc length 0 is
/// the foot. Positive kappa0 bends toward +N.
inline void hopf_circle_at(real kappa0, real s, real &x, real &y) {
  const real th = kappa0 * s;
  x = std::sin(th) / kappa0;
  y = (real(1) - std::cos(th)) / kappa0;
}

inline bool hopf_circle_closest(real kappa0, real px, real py, real &qx,
                                real &qy) {
  if (!(std::abs(kappa0) > real(1e-8)))
    return false;
  const real r = real(1) / kappa0;
  const real vx = px;
  const real vy = py - r;
  const real vn = std::sqrt(vx * vx + vy * vy);
  if (!(vn > real(1e-14)))
    return false;
  const real radius = std::abs(r);
  qx = radius * vx / vn;
  qy = r + radius * vy / vn;
  return std::isfinite(qx) && std::isfinite(qy);
}

struct hopf_table {
  static constexpr int N = 64;
  static constexpr int count = 2 * N + 1;
  real s[count]{};
  real x[count]{};
  real y[count]{};
  real th[count]{};
  real kappa0 = 0;
  real m = 0;
  real phase = 0;
  real frequency = 1;
  real S = 0;
  bool ok = false;
};

inline void hopf_rk4_step(real kappa0, real m, real phase, real frequency,
                          real h, real &s, real &x, real &y, real &th) {
  const real scale = std::abs(kappa0) * frequency;
  auto kappa_at = [&](real ss) {
    real sn = 0, cn = 1, dn = 1;
    jacobi_sn_cn_dn(scale * ss + phase, m, sn, cn, dn);
    return kappa0 * cn;
  };
  const real th1 = th;
  const real k1 = kappa_at(s);
  const real dx1 = std::cos(th1);
  const real dy1 = std::sin(th1);

  const real th2 = th + real(0.5) * h * k1;
  const real k2 = kappa_at(s + real(0.5) * h);
  const real dx2 = std::cos(th2);
  const real dy2 = std::sin(th2);

  const real th3 = th + real(0.5) * h * k2;
  const real k3 = kappa_at(s + real(0.5) * h);
  const real dx3 = std::cos(th3);
  const real dy3 = std::sin(th3);

  const real th4 = th + h * k3;
  const real k4 = kappa_at(s + h);
  const real dx4 = std::cos(th4);
  const real dy4 = std::sin(th4);

  s += h;
  x += (h / real(6)) * (dx1 + real(2) * dx2 + real(2) * dx3 + dx4);
  y += (h / real(6)) * (dy1 + real(2) * dy2 + real(2) * dy3 + dy4);
  th += (h / real(6)) * (k1 + real(2) * k2 + real(2) * k3 + k4);
}

/// Profile κ(s) = κ0 cn(frequency |κ0| s + phase, m). frequency 1 is the
/// elastic scaling. A larger frequency packs the wave into a shorter arc.
inline hopf_table hopf_make_table(real kappa0, real m, real phase = 0,
                                  real frequency = 1) {
  hopf_table t;
  t.kappa0 = kappa0;
  t.m = m;
  t.phase = phase;
  t.frequency = std::isfinite(frequency) && std::abs(frequency) > real(1e-8)
                    ? frequency
                    : real(1);
  const real scale = std::abs(kappa0);
  if (!(scale > real(1e-8)) || !(m > real(0)) || !(m < real(1)))
    return t;
  t.S = real(6) / scale;
  const real h = t.S / real(hopf_table::N);
  const int mid = hopf_table::N;
  for (int i = 0; i < hopf_table::N; ++i) {
    const int fwd = mid + i;
    t.s[fwd + 1] = t.s[fwd];
    t.x[fwd + 1] = t.x[fwd];
    t.y[fwd + 1] = t.y[fwd];
    t.th[fwd + 1] = t.th[fwd];
    hopf_rk4_step(kappa0, m, phase, t.frequency, h, t.s[fwd + 1], t.x[fwd + 1],
                  t.y[fwd + 1], t.th[fwd + 1]);
    const int back = mid - i;
    t.s[back - 1] = t.s[back];
    t.x[back - 1] = t.x[back];
    t.y[back - 1] = t.y[back];
    t.th[back - 1] = t.th[back];
    hopf_rk4_step(kappa0, m, phase, t.frequency, -h, t.s[back - 1],
                  t.x[back - 1], t.y[back - 1], t.th[back - 1]);
  }
  t.ok = true;
  return t;
}

struct hopf_profile {
  vec3 axis = vec3::Zero();
  vec3 tube = vec3::Zero();
  vec3 N = vec3::Zero();
  real kappa0 = 0;
  real m = 0;
  real phase = 0;
  hopf_table table;
  bool spin = false;
  bool ok = false;
};

/// One spoke of the wave fit. `sigma` is arc length along the profile,
/// `alpha` along the fiber. `kappa` is the spoke's normal curvature.
struct hopf_wave_sample {
  real sigma = 0;
  real alpha = 0;
  real kappa = 0;
};

/// Slope of profile curvature along the wave. A constant principal jet gives
/// 0, so the foot stays at phase 0.
inline real hopf_wave_slope(real k_hi, real k_lo, const hopf_wave_sample *samples,
                            int n) {
  if (samples == nullptr || n <= 0 || !std::isfinite(k_hi) || !std::isfinite(k_lo))
    return 0;
  real num = 0;
  real den = 0;
  for (int i = 0; i < n; ++i) {
    const real sigma = samples[i].sigma;
    const real alpha = samples[i].alpha;
    const real L2 = sigma * sigma + alpha * alpha;
    if (!(L2 > real(1e-16)) || !std::isfinite(samples[i].kappa))
      continue;
    const real pred = (k_hi * sigma * sigma + k_lo * alpha * alpha) / L2;
    const real c = (sigma * sigma * sigma) / (real(2) * L2);
    const real r = samples[i].kappa - pred;
    if (!std::isfinite(c) || !std::isfinite(r))
      continue;
    num += c * r;
    den += c * c;
  }
  if (!(den > real(1e-20)))
    return 0;
  const real slope = num / den;
  return std::isfinite(slope) ? slope : 0;
}

/// Jacobi phase at the foot. `kappa` is the foot curvature, `slope` is
/// dκ/ds along the wave, `m` is the modulus. Phase 0 is the crest. The
/// amplitude is kappa / cn(phase), so the foot curvature stays `kappa`.
inline real hopf_wave_phase(real kappa, real slope, real m,
                            real frequency = 1) {
  if (!std::isfinite(kappa) || !std::isfinite(slope) || !std::isfinite(m))
    return 0;
  if (!(std::abs(kappa) > real(1e-8)) || !(std::abs(slope) > real(1e-8)))
    return 0;
  if (!(m > real(0)) || !(m < real(1)))
    return 0;
  const real freq = std::isfinite(frequency) && std::abs(frequency) > real(1e-8)
                        ? frequency
                        : real(1);
  auto predicted = [&](real phi) {
    real sn = 0, cn = 1, dn = 1;
    jacobi_sn_cn_dn(phi, m, sn, cn, dn);
    if (!(std::abs(cn) > real(0.25)))
      return std::numeric_limits<real>::quiet_NaN();
    const real amp = kappa / cn;
    return -amp * std::abs(amp) * freq * sn * dn;
  };
  real best = 0;
  real best_err = std::abs(slope);
  constexpr int n = 48;
  for (int i = 0; i <= n; ++i) {
    const real phi = real(-1) + real(2) * real(i) / real(n);
    const real p = predicted(phi);
    if (!std::isfinite(p))
      continue;
    const real err = std::abs(p - slope);
    if (err < best_err) {
      best_err = err;
      best = phi;
    }
  }
  return best;
}

inline hopf_profile hopf_profile_make(const vec3 &axis, const vec3 &tube,
                                      const vec3 &N, real kappa0, real twist,
                                      real phase = 0, real frequency = 1) {
  hopf_profile hp;
  if (!axis.allFinite() || !tube.allFinite() || !N.allFinite())
    return hp;
  if (!(axis.norm() > real(1e-12)) || !(tube.norm() > real(1e-12)) ||
      !(N.norm() > real(1e-12)))
    return hp;
  hp.axis = axis.normalized();
  hp.tube = tube.normalized();
  hp.N = N.normalized();
  hp.kappa0 = kappa0;
  hp.m = hopf_modulus(twist);
  hp.ok = std::isfinite(hp.kappa0) && std::abs(hp.kappa0) > real(1e-8);
  if (!(hp.ok && hp.m > real(0)))
    return hp;
  real amp = hp.kappa0;
  if (std::isfinite(phase) && std::abs(phase) > real(1e-10)) {
    real sn = 0, cn = 1, dn = 1;
    jacobi_sn_cn_dn(phase, hp.m, sn, cn, dn);
    if (std::abs(cn) > real(0.25)) {
      amp = hp.kappa0 / cn;
      hp.phase = phase;
    }
  }
  if (!(std::abs(amp) > real(1e-8)) || !std::isfinite(amp)) {
    hp.ok = false;
    return hp;
  }
  hp.table = hopf_make_table(amp, hp.m, hp.phase, frequency);
  return hp;
}

/// Curvature of the unit sphere in Pinkall's shape operator, the off-diagonal
/// of W. On a mesh this is a world-unit curvature.
inline constexpr real k_hopf_sphere_curvature = real(1);

/// κ(s) = A cn(frequency |A| s + phase, m), with A chosen so κ(0) is H.
/// m outside (0, 1) stays at the measured mean curvature.
inline real hopf_elastic_kappa(real H, real m, real frequency, real phase,
                               real s) {
  if (!std::isfinite(H))
    return std::numeric_limits<real>::quiet_NaN();
  if (!(m > real(0)) || !(m < real(1)))
    return H;
  real phi = 0;
  real amp = H;
  if (std::isfinite(phase) && std::abs(phase) > real(1e-10)) {
    real sn = 0, cn = 1, dn = 1;
    jacobi_sn_cn_dn(phase, m, sn, cn, dn);
    if (std::isfinite(cn) && std::abs(cn) > real(0.25)) {
      amp = H / cn;
      phi = phase;
    }
  }
  if (!std::isfinite(amp))
    return H;
  const real freq = std::isfinite(frequency) ? frequency : real(1);
  const real arg = freq * std::abs(amp) * s + phi;
  real sn = 0, cn = 1, dn = 1;
  jacobi_sn_cn_dn(arg, m, sn, cn, dn);
  const real kappa = amp * cn;
  return std::isfinite(kappa) ? kappa : H;
}

/// Height of Pinkall's quadric. In the frame (across the fiber, along it),
/// W = [[2κ, σ], [σ, 0]] and the height is ½ tᵀ W t = κ u² + σ u v.
inline real hopf_pinkall_height(real kappa, real u, real v) {
  if (!std::isfinite(kappa) || !std::isfinite(u) || !std::isfinite(v))
    return std::numeric_limits<real>::quiet_NaN();
  return kappa * u * u + k_hopf_sphere_curvature * u * v;
}

/// Larger-magnitude eigenvalue of Pinkall's W. This is the tube curvature.
inline real hopf_tube_curvature(real kappa) {
  if (!std::isfinite(kappa))
    return std::numeric_limits<real>::quiet_NaN();
  const real sig = k_hopf_sphere_curvature;
  const real s = std::sqrt(kappa * kappa + sig * sig);
  const real lp = kappa + s;
  const real lm = kappa - s;
  return std::abs(lp) >= std::abs(lm) ? lp : lm;
}

/// Below this |scale|, the circle center would run away. The floor limits
/// |curvature| and keeps the sign, which is the side of the tangent plane.
inline constexpr real k_hopf_coeff_floor = real(0.2);

/// Project (u, h) onto the circle of curvature `k_hi * coeff`. The center is
/// 1/(k_hi * coeff): the sign of that curvature is the side of the tangent
/// plane, and the radius moves with the center. The circle stays tangent at
/// the foot, so a looser radius does not throw the foot through the plane.
/// coeff = 1 is the measured circle.
inline bool hopf_cylinder_land(real k_hi, real coeff, real u, real h, real &uo,
                               real &ho) {
  if (!std::isfinite(k_hi) || !std::isfinite(coeff) || !std::isfinite(u) ||
      !std::isfinite(h) || !(std::abs(k_hi) > real(1e-8)))
    return false;
  const real mag = std::max(std::abs(coeff), k_hopf_coeff_floor);
  const real signed_coeff = std::copysign(mag, coeff);
  const real k = k_hi * signed_coeff;
  if (!(std::abs(k) > real(1e-8)))
    return false;
  const real rho = real(1) / std::abs(k);
  const real ch = real(1) / k;
  const real du = u;
  const real dh = h - ch;
  const real n2 = du * du + dh * dh;
  if (!(n2 > real(1e-16)))
    return false;
  const real inv = rho / std::sqrt(n2);
  uo = du * inv;
  ho = ch + dh * inv;
  return std::isfinite(uo) && std::isfinite(ho);
}

/// Project (u, h) onto the circle of radius 1/|λ+| about the measured tube
/// center (0, 1/k_hi). The center stays. The radius changes.
inline bool hopf_radial_land(real k_hi, real kappa, real u, real h, real &uo,
                             real &ho) {
  if (!std::isfinite(k_hi) || !std::isfinite(u) || !std::isfinite(h) ||
      !(std::abs(k_hi) > real(1e-8)))
    return false;
  const real lam = hopf_tube_curvature(kappa);
  if (!std::isfinite(lam) || !(std::abs(lam) > real(1e-8)))
    return false;
  const real rho = real(1) / std::abs(lam);
  const real ch = real(1) / k_hi;
  const real du = u;
  const real dh = h - ch;
  const real n2 = du * du + dh * dh;
  if (!(n2 > real(1e-16)))
    return false;
  const real inv = rho / std::sqrt(n2);
  uo = du * inv;
  ho = ch + dh * inv;
  return std::isfinite(uo) && std::isfinite(ho);
}

/// One Hopf jet. `axis` is the fiber, `tube` is across it. `H` is the mean
/// curvature at the foot. `k_hi` places the measured tube center.
struct hopf_jet {
  vec3 axis = vec3::Zero();
  vec3 tube = vec3::Zero();
  real H = 0;
  real k_hi = 0;
  real m = 0;
  real phase = 0;
  real frequency = 1;
  bool ok = false;
};

/// Rotate the fiber about the surface normal. angle 0 leaves the ruling on
/// `axis0` and the wave on `tube0`.
inline void hopf_fiber_axes(const vec3 &axis0, const vec3 &tube0, real angle,
                            vec3 &axis, vec3 &tube) {
  const real c = std::cos(std::isfinite(angle) ? angle : real(0));
  const real s = std::sin(std::isfinite(angle) ? angle : real(0));
  axis = c * axis0 + s * tube0;
  tube = c * tube0 - s * axis0;
}

/// One neighbor of the connection eigenproblem. `t` is the tangent
/// displacement from the foot. `spoke_z` is the spin about the normal from
/// the foot frame. `link_z` is the spin from this neighbor to the next one
/// in circulation order; the last entry closes the ring.
struct hopf_ring_node {
  vec3 t = vec3::Zero();
  real spoke_z = 0;
  real link_z = 0;
};

/// A tangent and its opposite are the same ruling, so a spin of π is no spin.
/// The result lies in [-π/2, π/2].
inline real hopf_fold_line_spin(real omega) {
  if (!std::isfinite(omega))
    return 0;
  return real(0.5) * std::atan2(std::sin(real(2) * omega), std::cos(real(2) * omega));
}

/// Angle of `transport` in the frame built from the one-ring connection
/// Laplacian. The connection is the quaternion spin about the normal, folded
/// by π. The smallest eigenvector's argument steps define a tangent V. Phase
/// 0 means the transport line lies along V. A flat or negligible spin
/// returns 0.
inline real hopf_ring_phase(const vec3 &N, const vec3 &transport,
                            const hopf_ring_node *nodes, int n) {
  if (nodes == nullptr || n < 3 || !N.allFinite() || !transport.allFinite())
    return 0;
  if (!(N.squaredNorm() > real(1e-16)))
    return 0;
  const int m = n + 1;
  Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> A =
      Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic>::Zero(m, m);
  Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> B =
      Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic>::Zero(m, m);
  auto add_edge = [&](int i, int j, real omega) {
    if (!std::isfinite(omega))
      return;
    const real c = std::cos(omega);
    const real s = std::sin(omega);
    A(i, i) += 1;
    A(j, j) += 1;
    A(i, j) -= c;
    A(j, i) -= c;
    B(i, j) += s;
    B(j, i) -= s;
  };
  real spoke_len = 0;
  int spoke_n = 0;
  for (int k = 0; k < n; ++k) {
    add_edge(0, k + 1, hopf_fold_line_spin(nodes[k].spoke_z));
    add_edge(k + 1, ((k + 1) % n) + 1, hopf_fold_line_spin(nodes[k].link_z));
    const real len = nodes[k].t.norm();
    if (std::isfinite(len) && len > real(0)) {
      spoke_len += len;
      ++spoke_n;
    }
  }
  Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> H(2 * m, 2 * m);
  H.topLeftCorner(m, m) = A;
  H.topRightCorner(m, m) = -B;
  H.bottomLeftCorner(m, m) = B;
  H.bottomRightCorner(m, m) = A;
  Eigen::SelfAdjointEigenSolver<
      Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic>>
      es(H);
  if (es.info() != Eigen::Success)
    return 0;
  const Eigen::Matrix<real, Eigen::Dynamic, 1> ev = es.eigenvectors().col(0);
  const real u0 = ev(0);
  const real v0 = ev(m);
  const real n0 = std::hypot(u0, v0);
  if (!(n0 > real(1e-14)))
    return 0;
  const real c0 = u0 / n0;
  const real s0 = v0 / n0;
  const vec3 Nn = N.normalized();
  vec3 E1 = Nn.unitOrthogonal();
  vec3 E2 = Nn.cross(E1);
  real ata00 = 0, ata01 = 0, ata11 = 0, atd0 = 0, atd1 = 0;
  for (int k = 0; k < n; ++k) {
    const real uk = ev(k + 1);
    const real vk = ev(m + k + 1);
    const real u = uk * c0 + vk * s0;
    const real v = vk * c0 - uk * s0;
    const real delta = std::atan2(v, u);
    const vec3 tk = nodes[k].t - nodes[k].t.dot(Nn) * Nn;
    const real x = tk.dot(E1);
    const real y = tk.dot(E2);
    ata00 += x * x;
    ata01 += x * y;
    ata11 += y * y;
    atd0 += x * delta;
    atd1 += y * delta;
  }
  const real det = ata00 * ata11 - ata01 * ata01;
  if (!(std::abs(det) > real(1e-16)))
    return 0;
  const real vx = (ata11 * atd0 - ata01 * atd1) / det;
  const real vy = (ata00 * atd1 - ata01 * atd0) / det;
  vec3 V = vx * E1 + vy * E2;
  const real mean_len =
      spoke_n > 0 ? spoke_len / real(spoke_n) : real(1);
  // |V| is radians per length. Below this, the argument steps are noise.
  if (!(V.norm() * mean_len > real(1e-3)))
    return 0;
  V.normalize();
  vec3 L = transport - transport.dot(Nn) * Nn;
  if (!(L.squaredNorm() > real(1e-16)))
    return 0;
  L.normalize();
  const real phase = std::atan2(L.dot(Nn.cross(V)), L.dot(V));
  return std::isfinite(phase) ? phase : 0;
}

/// One one-ring frame. `t` is the tangent displacement from the foot. `q` is
/// that neighbor's quaternion. Samples are in circulation order.
struct hopf_section_sample {
  vec3 t = vec3::Zero();
  quat q = quat::Identity();
};

/// Angle of a normal about `fiber`, in the plane perpendicular to it.
/// The zero is the projection of world Z into that plane, or world X when the
/// fiber is parallel to Z.
inline real hopf_normal_angle(const vec3 &Ns_in, const vec3 &fiber) {
  if (!Ns_in.allFinite() || !fiber.allFinite() ||
      !(Ns_in.squaredNorm() > real(1e-16)) ||
      !(fiber.squaredNorm() > real(1e-16)))
    return 0;
  const vec3 R = fiber.normalized();
  vec3 ref = vec3::UnitZ() - R * R.dot(vec3::UnitZ());
  if (!(ref.squaredNorm() > real(1e-8)))
    ref = vec3::UnitX() - R * R.dot(vec3::UnitX());
  if (!(ref.squaredNorm() > real(1e-16)))
    return 0;
  ref.normalize();
  const vec3 Ns = Ns_in.normalized();
  const real ang = std::atan2(Ns.dot(R.cross(ref)), Ns.dot(ref));
  return std::isfinite(ang) ? ang : 0;
}

inline real hopf_normal_angle(const quat &q_in, const vec3 &fiber) {
  quat q = q_in;
  const real n2 = q.squaredNorm();
  if (!(n2 > real(1e-16)) || !std::isfinite(n2))
    return 0;
  q.normalize();
  return hopf_normal_angle(q.toRotationMatrix().col(2), fiber);
}

/// Major tangent Z × (x, y, 0) and the tube perpendicular to it. On a torus of
/// revolution the tube is the minor-circle direction, one smooth axis.
inline bool hopf_revolution_axes(const vec3 &point, const vec3 &Nin, vec3 &major,
                                 vec3 &tube) {
  if (!point.allFinite() || !Nin.allFinite() || !(Nin.squaredNorm() > real(1e-16)))
    return false;
  const vec3 N = Nin.normalized();
  const vec3 radial(point.x(), point.y(), real(0));
  vec3 M = vec3::UnitZ().cross(radial);
  if (!(M.squaredNorm() > real(1e-8)))
    M = vec3::UnitZ().cross(N);
  if (!(M.squaredNorm() > real(1e-16)))
    return false;
  M.normalize();
  M -= N * N.dot(M);
  if (!(M.squaredNorm() > real(1e-16)))
    return false;
  M.normalize();
  const vec3 T = N.cross(M);
  if (!(T.squaredNorm() > real(1e-16)))
    return false;
  major = M;
  tube = T.normalized();
  return true;
}

/// One neighbor used to read normal curvature. `t` is the displacement from
/// the foot. `N` is that neighbor's normal.
struct hopf_turn_sample {
  vec3 t = vec3::Zero();
  vec3 N = vec3::Zero();
};

/// Normal curvature in direction `dir`: Σ (t·dir) ((N_j − N)·dir) / Σ (t·dir)².
/// The sign is the side of the tangent plane the center sits on when this
/// value is used as `k` in the circle landing.
inline real hopf_direction_curvature(const vec3 &N, const vec3 &dir,
                                     const hopf_turn_sample *samples, int n) {
  if (samples == nullptr || n < 1 || !N.allFinite() || !dir.allFinite() ||
      !(dir.squaredNorm() > real(1e-16)))
    return 0;
  const vec3 D = dir.normalized();
  real num = 0;
  real den = 0;
  for (int k = 0; k < n; ++k) {
    const vec3 &t = samples[k].t;
    const vec3 &Nj = samples[k].N;
    if (!t.allFinite() || !Nj.allFinite())
      continue;
    const real ds = t.dot(D);
    num += ds * (Nj - N).dot(D);
    den += ds * ds;
  }
  if (!(den > real(1e-16)))
    return 0;
  const real kappa = num / den;
  return std::isfinite(kappa) ? kappa : 0;
}

/// A principal axis has no sign. Flip `fiber`, and `tube` with it, so the fiber
/// agrees with the major tangent Z × (x, y, 0). The angle about the opposite
/// fiber is the negative phase. When the foot sits on the Z axis, fall back to
/// Z × N.
inline void hopf_orient_fiber(const vec3 &point, const vec3 &N, vec3 &fiber,
                              vec3 &tube) {
  if (!point.allFinite() || !fiber.allFinite() || !tube.allFinite())
    return;
  const vec3 radial(point.x(), point.y(), real(0));
  vec3 ref = vec3::UnitZ().cross(radial);
  if (!(ref.squaredNorm() > real(1e-8)) && N.allFinite())
    ref = vec3::UnitZ().cross(N);
  if (!(ref.squaredNorm() > real(1e-8)) || !(fiber.squaredNorm() > real(1e-16)))
    return;
  if (fiber.dot(ref) < real(0)) {
    fiber = -fiber;
    tube = -tube;
  }
}

/// Shift the osculating angle so `twist` is the new zero.
inline real hopf_phase_shift(real phase, real twist) {
  if (!std::isfinite(phase))
    return 0;
  const real t = std::isfinite(twist) ? twist : real(0);
  return phase - t;
}

inline real hopf_wrap_angle(real a) {
  if (!std::isfinite(a))
    return 0;
  return std::atan2(std::sin(a), std::cos(a));
}

/// Angle of the position about the world Z axis.
inline real hopf_major_angle(const vec3 &point) {
  if (!point.allFinite())
    return 0;
  return std::atan2(point.y(), point.x());
}

/// Tube radius of the measured section. Zero when the curvature vanishes.
inline real hopf_tube_radius(real k_hi) {
  if (!std::isfinite(k_hi) || !(std::abs(k_hi) > real(1e-8)))
    return 0;
  return real(1) / std::abs(k_hi);
}

/// Phase constant along the ruling rotated by `fiber` from the principal axis.
/// `minor` is the angle about that axis and `major` is `hopf_major_angle`.
/// `radius` is the measured tube radius and `rho` is the horizontal distance
/// from the Z axis, so the gradient follows the rotated tube. `twist` shifts
/// the zero. The result lies in [-π, π].
inline real hopf_aligned_phase(real minor, real major, real fiber, real rho,
                               real radius, real twist) {
  if (!std::isfinite(minor))
    return 0;
  const real th = std::isfinite(fiber) ? fiber : real(0);
  const real maj = std::isfinite(major) ? major : real(0);
  real scale = 1;
  if (std::isfinite(rho) && std::isfinite(radius) && radius > real(1e-8))
    scale = rho / radius;
  const real mixed = std::cos(th) * minor - std::sin(th) * scale * maj;
  return hopf_wrap_angle(hopf_phase_shift(mixed, twist));
}

/// Normal curvature of the section perpendicular to the Hopf ruling.
/// `fiber` is the rotation of that ruling off the low-curvature axis, so the
/// tube direction is the same angle off the high-curvature axis.
inline real hopf_section_curvature(real k_hi, real k_lo, real fiber) {
  const real th = std::isfinite(fiber) ? fiber : real(0);
  const real c = std::cos(th);
  const real s = std::sin(th);
  const real kh = std::isfinite(k_hi) ? k_hi : real(0);
  const real kl = std::isfinite(k_lo) ? k_lo : real(0);
  return kh * c * c + kl * s * s;
}

/// Move from the measured section circle to the phase-scaled one. A scale of
/// 1 leaves the spoke where it is. The displacement lies in the plane
/// perpendicular to the Hopf ruling.
inline bool hopf_wave_move(real k, real scale, real u, real h, real &du,
                           real &dh) {
  real u1 = 0;
  real h1 = 0;
  real us = 0;
  real hs = 0;
  if (!hopf_cylinder_land(k, 1, u, h, u1, h1))
    return false;
  if (!hopf_cylinder_land(k, scale, u, h, us, hs))
    return false;
  du = us - u1;
  dh = hs - h1;
  return std::isfinite(du) && std::isfinite(dh);
}

/// Signed fraction of the measured section curvature. `sin = ±1` targets
/// radius ±R, where R = 1/|k| is the radius already on the surface. The sign
/// is the side of the tangent plane. Near a zero the landing floor keeps the
/// radius from running away.
inline constexpr real k_hopf_wave_amplitude = real(1);

inline real hopf_wave_scale(real frequency, real phase) {
  const real freq = std::isfinite(frequency) ? frequency : real(1);
  const real ph = std::isfinite(phase) ? phase : real(0);
  return k_hopf_wave_amplitude * std::sin(freq * ph);
}

/// Project (u, h) onto the circle of curvature `k`. The center is 1/k, so the
/// circle stays tangent at the foot. The sign of `k` is the side of the
/// tangent plane.
inline bool hopf_curvature_land(real k, real u, real h, real &uo, real &ho) {
  if (!std::isfinite(k) || !std::isfinite(u) || !std::isfinite(h) ||
      !(std::abs(k) > real(1e-8)))
    return false;
  const real rho = real(1) / std::abs(k);
  const real ch = real(1) / k;
  const real du = u;
  const real dh = h - ch;
  const real n2 = du * du + dh * dh;
  if (!(n2 > real(1e-16)))
    return false;
  const real inv = rho / std::sqrt(n2);
  uo = du * inv;
  ho = ch + dh * inv;
  return std::isfinite(uo) && std::isfinite(ho);
}

/// Displacement from the circle of curvature `k0` to the circle of curvature
/// `k1`. Equal curvatures leave the spoke where it is.
inline bool hopf_curvature_move(real k0, real k1, real u, real h, real &du,
                                real &dh) {
  real u0 = 0;
  real h0 = 0;
  real u1 = 0;
  real h1 = 0;
  if (!hopf_curvature_land(k0, u, h, u0, h0))
    return false;
  if (!hopf_curvature_land(k1, u, h, u1, h1))
    return false;
  du = u1 - u0;
  dh = h1 - h0;
  return std::isfinite(du) && std::isfinite(dh);
}

/// Target section curvature. `x` is the vertex world coordinate.
inline real hopf_curvature_plus_sin(real k0, real x) {
  if (!std::isfinite(k0) || !std::isfinite(x))
    return std::numeric_limits<real>::quiet_NaN();
  return k0 + std::sin(x);
}

/// Spoke offset along the normal so the cotangent sum changes by
/// `4 * area * sin(x) * N`. The row is Σ K (x_j − x_i) = 4 H A N, with
/// Δ = Σ K (·) / (2A), so this is H += sin(x). Positive x adds the outward
/// normal. The same offset goes on every spoke; the sum is the only quantity
/// the row keeps.
inline real hopf_curvature_normal_offset(real area, real sum_k, real x) {
  if (!std::isfinite(area) || !std::isfinite(sum_k) || !std::isfinite(x))
    return 0;
  if (!(area > real(0)) || !(std::abs(sum_k) > real(1e-12)))
    return 0;
  const real beta = real(4) * area * std::sin(x) / sum_k;
  return std::isfinite(beta) ? beta : 0;
}

/// Phase in the osculating cross-section. The fiber is oriented about world Z
/// first. It then falls between a consecutive pair of one-ring frames; the
/// opposite direction falls between another pair. Each pair is slerped to the
/// fiber, and the phase is the average of those normals' angles about the fiber.
inline real hopf_section_phase(const vec3 &point, const vec3 &N, const vec3 &fiber_in,
                               const vec3 &tube_in, const hopf_section_sample *samples,
                               int n) {
  if (samples == nullptr || n < 2 || !N.allFinite() || !fiber_in.allFinite() ||
      !tube_in.allFinite())
    return 0;
  if (!(N.squaredNorm() > real(1e-16)) ||
      !(fiber_in.squaredNorm() > real(1e-16)) ||
      !(tube_in.squaredNorm() > real(1e-16)))
    return 0;
  vec3 fiber = fiber_in;
  vec3 tube = tube_in;
  hopf_orient_fiber(point, N, fiber, tube);
  const vec3 Nn = N.normalized();
  const vec3 R = fiber.normalized();
  const vec3 T = tube.normalized();
  auto wrap = [](real a) { return std::atan2(std::sin(a), std::cos(a)); };
  std::vector<real> ang(static_cast<size_t>(n));
  for (int k = 0; k < n; ++k) {
    const vec3 tk = samples[k].t - samples[k].t.dot(Nn) * Nn;
    ang[static_cast<size_t>(k)] = std::atan2(tk.dot(T), tk.dot(R));
  }
  auto wedge = [&](int i, int j, real target, real &fraction) {
    const real dab = wrap(ang[static_cast<size_t>(j)] - ang[static_cast<size_t>(i)]);
    const real dat = wrap(target - ang[static_cast<size_t>(i)]);
    if (!(std::abs(dab) > real(1e-8)))
      return false;
    if (!(dat * dab >= real(0)) || std::abs(dat) > std::abs(dab) + real(1e-6))
      return false;
    fraction = dat / dab;
    return std::isfinite(fraction);
  };
  auto side = [&](real target, real &phase) {
    real best = std::numeric_limits<real>::infinity();
    int ia = -1;
    int ib = -1;
    real frac = 0;
    for (int k = 0; k < n; ++k) {
      const int next = (k + 1) % n;
      real f = 0;
      if (!wedge(k, next, target, f))
        continue;
      const real width =
          std::abs(wrap(ang[static_cast<size_t>(next)] - ang[static_cast<size_t>(k)]));
      if (width < best) {
        best = width;
        ia = k;
        ib = next;
        frac = f;
      }
    }
    if (ia < 0)
      return false;
    quat qa = samples[ia].q;
    quat qb = samples[ib].q;
    const real na = qa.squaredNorm();
    const real nb = qb.squaredNorm();
    if (!(na > real(1e-16)) || !(nb > real(1e-16)))
      return false;
    qa.normalize();
    qb.normalize();
    const quat qs = qa.slerp(std::clamp(frac, real(0), real(1)), qb);
    phase = hopf_normal_angle(qs, R);
    return std::isfinite(phase);
  };
  real fore = 0;
  real back = 0;
  const bool have_fore = side(0, fore);
  const bool have_back = side(real(M_PI), back);
  if (have_fore && have_back)
    return real(0.5) * (fore + back);
  if (have_fore)
    return fore;
  if (have_back)
    return back;
  return 0;
}

inline void hopf_advance_to(const hopf_table &table, int i0, real s_target,
                            real &x, real &y) {
  real s = table.s[i0];
  x = table.x[i0];
  y = table.y[i0];
  real th = table.th[i0];
  const real span = s_target - s;
  if (!(std::abs(span) > real(1e-14)))
    return;
  const int n = 4;
  const real h = span / real(n);
  for (int i = 0; i < n; ++i)
    hopf_rk4_step(table.kappa0, table.m, table.phase, table.frequency, h, s, x,
                  y, th);
}

inline bool hopf_table_closest(const hopf_table &table, real px, real py,
                               real &qx, real &qy) {
  if (!table.ok)
    return false;
  int best = 0;
  real best_d2 = std::numeric_limits<real>::infinity();
  for (int i = 0; i < hopf_table::count; ++i) {
    const real dx = table.x[i] - px;
    const real dy = table.y[i] - py;
    const real d2 = dx * dx + dy * dy;
    if (d2 < best_d2) {
      best_d2 = d2;
      best = i;
    }
  }
  const int left = std::max(0, best - 1);
  const int right = std::min(hopf_table::count - 1, best + 1);
  real a = table.s[left];
  real b = table.s[right];
  if (!(b > a)) {
    qx = table.x[best];
    qy = table.y[best];
    return true;
  }
  const real inv_phi = real(0.5) * (std::sqrt(real(5)) - real(1));
  real c = b - inv_phi * (b - a);
  real d = a + inv_phi * (b - a);
  auto dist2 = [&](real s) {
    real x = 0, y = 0;
    hopf_advance_to(table, left, s, x, y);
    const real dx = x - px;
    const real dy = y - py;
    return dx * dx + dy * dy;
  };
  real fc = dist2(c);
  real fd = dist2(d);
  for (int iter = 0; iter < 18; ++iter) {
    if (fc < fd) {
      b = d;
      d = c;
      fd = fc;
      c = b - inv_phi * (b - a);
      fc = dist2(c);
    } else {
      a = c;
      c = d;
      fc = fd;
      d = a + inv_phi * (b - a);
      fd = dist2(d);
    }
  }
  const real s = real(0.5) * (a + b);
  real rx = 0, ry = 0;
  hopf_advance_to(table, left, s, rx, ry);
  const real dx = rx - px;
  const real dy = ry - py;
  if (dx * dx + dy * dy < best_d2) {
    qx = rx;
    qy = ry;
  } else {
    qx = table.x[best];
    qy = table.y[best];
  }
  return std::isfinite(qx) && std::isfinite(qy);
}

/// Body spin about N larger than this is a frame flip, not a Hopf phase.
/// The tube probe's quiet spokes sit near 0.02–0.11; the flips start near 1.
inline constexpr real k_hopf_spin_max = real(0.5);

/// Roll a Hopf landing about the fiber through the foot by this spoke's ωz.
/// That is the relative fiber phase. The wave coordinate is unchanged in
/// direction: a spoke along the ruling stays. A flip is returned unchanged.
inline vec3 hopf_project_fiber(const vec3 &hit, const vec3 &axis, real wz) {
  if (!hit.allFinite() || !axis.allFinite() || !std::isfinite(wz))
    return hit;
  if (!(std::abs(wz) > real(1e-8)) || std::abs(wz) > k_hopf_spin_max)
    return hit;
  const real a2 = axis.squaredNorm();
  if (!(a2 > real(1e-16)))
    return hit;
  const vec3 a = axis / std::sqrt(a2);
  const real c = std::cos(wz);
  const real s = std::sin(wz);
  const real along = a.dot(hit);
  return c * hit + s * a.cross(hit) + (real(1) - c) * along * a;
}

/// Turn a Hopf landing about the normal by this spoke's body spin ωz.
/// The height stays. The tangential part swings. A flip is returned unchanged.
/// `hopf_sweep::spin` gates the call. Set that false to undo the trial.
inline vec3 hopf_project_spin(const vec3 &hit, const vec3 &N, real wz) {
  if (!hit.allFinite() || !N.allFinite() || !std::isfinite(wz))
    return hit;
  if (!(std::abs(wz) > real(1e-8)) || std::abs(wz) > k_hopf_spin_max)
    return hit;
  const real n2 = N.squaredNorm();
  if (!(n2 > real(1e-16)))
    return hit;
  const vec3 n = N / std::sqrt(n2);
  const real c = std::cos(wz);
  const real s = std::sin(wz);
  const real hn = hit.dot(n);
  const vec3 t = hit - hn * n;
  return c * t + s * n.cross(t) + hn * n;
}

/// Closest point on the swept profile. The component along `axis` stays.
/// m = 0 is the circle of curvature kappa0. A correction longer than the
/// spoke is rejected.
inline bool land_hopf(const hopf_profile &hp, const vec3 &dp, vec3 &out) {
  if (!hp.ok || !dp.allFinite())
    return false;
  const real axial = dp.dot(hp.axis);
  const real px = dp.dot(hp.tube);
  const real py = dp.dot(hp.N);
  real qx = 0, qy = 0;
  const bool hit = hp.m > real(0)
                       ? hopf_table_closest(hp.table, px, py, qx, qy)
                       : hopf_circle_closest(hp.kappa0, px, py, qx, qy);
  if (!hit)
    return false;
  out = axial * hp.axis + qx * hp.tube + qy * hp.N;
  if (!out.allFinite())
    return false;
  if ((out - dp).norm() > dp.norm() + real(1e-8))
    return false;
  return true;
}

} // namespace block
} // namespace hepworth
} // namespace gaudi
