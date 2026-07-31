#pragma once

#include <cmath>

#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>

namespace gaudi {
namespace vermeer {
namespace dq {

// Dual quaternion for SE(3): real = rotation, dual encodes translation.
struct dual_quat {
  glm::quat real = glm::quat(1, 0, 0, 0);
  glm::quat dual = glm::quat(0, 0, 0, 0);
};

inline dual_quat from_rigid(const glm::quat &r, const glm::vec3 &t) {
  dual_quat dq;
  dq.real = glm::normalize(r);
  // dual = 0.5 * t_quat * real
  const glm::quat tq(0.0f, t.x, t.y, t.z);
  dq.dual = 0.5f * (tq * dq.real);
  return dq;
}

inline glm::vec3 translation(const dual_quat &dq) {
  const glm::quat tq = 2.0f * (dq.dual * glm::conjugate(dq.real));
  return glm::vec3(tq.x, tq.y, tq.z);
}

inline dual_quat conjugate(const dual_quat &dq) {
  return {glm::conjugate(dq.real), glm::conjugate(dq.dual)};
}

inline dual_quat mul(const dual_quat &a, const dual_quat &b) {
  return {a.real * b.real, a.real * b.dual + a.dual * b.real};
}

inline float real_dot(const dual_quat &a, const dual_quat &b) {
  return glm::dot(a.real, b.real);
}

inline dual_quat normalize(const dual_quat &dq) {
  const float n = glm::length(dq.real);
  if (n < 1.0e-12f)
    return dual_quat{};
  dual_quat out;
  out.real = dq.real / n;
  out.dual = dq.dual / n;
  // Enforce <real, dual> = 0
  const float d = glm::dot(out.real, out.dual);
  out.dual -= out.real * d;
  return out;
}

// Screw linear interpolation (ScLERP) between unit dual quaternions.
inline dual_quat sclerp(dual_quat a, dual_quat b, float t) {
  a = normalize(a);
  b = normalize(b);
  if (real_dot(a, b) < 0.0f) {
    b.real = -b.real;
    b.dual = -b.dual;
  }

  // diff = a^* * b
  dual_quat diff = normalize(mul(conjugate(a), b));

  // Pow(diff, t) via screw parameters from real part.
  glm::quat qr = diff.real;
  float cos_half = glm::clamp(qr.w, -1.0f, 1.0f);
  float sin_half = std::sqrt(std::max(0.0f, 1.0f - cos_half * cos_half));
  float angle = 2.0f * std::atan2(sin_half, cos_half);

  glm::vec3 axis(qr.x, qr.y, qr.z);
  const float axis_n = glm::length(axis);
  if (axis_n > 1.0e-8f)
    axis /= axis_n;
  else
    axis = glm::vec3(0, 0, 1);

  // Pitch from dual part along screw.
  const glm::vec3 tr = translation(diff);
  const float pitch = glm::dot(tr, axis);

  const float at = angle * t;
  const float pt = pitch * t;
  const float sh = std::sin(0.5f * at);
  const float ch = std::cos(0.5f * at);

  dual_quat powt;
  powt.real = glm::quat(ch, axis.x * sh, axis.y * sh, axis.z * sh);
  // dual for pure screw: 0.5 * (pt * axis) * real  (approx for screw motion)
  const glm::vec3 t_axis = pt * axis;
  // Also include part of translation orthogonal to axis via linear blend fallback
  // when angle is tiny.
  if (std::abs(angle) < 1.0e-5f) {
    const glm::vec3 t_lerp = tr * t;
    return normalize(from_rigid(
        glm::slerp(a.real, b.real, t),
        glm::mix(translation(a), translation(b), t)));
  }

  // Full dual for screw: dual = 0.5 * (t_quat) * real with t = pitch*axis +
  // moment; use standard formula:
  // For unit dual quat, log/exp ScLERP:
  powt = from_rigid(powt.real, t_axis);
  // Blend in orthogonal translation component linearly (moment).
  const glm::vec3 t_orth = tr - pitch * axis;
  const glm::vec3 t_full = t_axis + t_orth * t;
  powt = from_rigid(powt.real, t_full);

  return normalize(mul(a, powt));
}

} // namespace dq
} // namespace vermeer
} // namespace gaudi
