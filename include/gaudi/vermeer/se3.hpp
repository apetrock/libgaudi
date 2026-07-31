#pragma once

#include <cmath>

#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>

namespace gaudi {
namespace vermeer {
namespace se3 {

struct twist {
  glm::vec3 w{0.0f}; // rotational (so(3))
  glm::vec3 v{0.0f}; // translational
};

struct rigid {
  glm::quat R{1, 0, 0, 0};
  glm::vec3 t{0.0f};
};

inline glm::mat3 hat(const glm::vec3 &w) {
  return glm::mat3(0.0f, -w.z, w.y, w.z, 0.0f, -w.x, -w.y, w.x, 0.0f);
}

inline rigid inverse(const rigid &T) {
  const glm::quat Ri = glm::conjugate(glm::normalize(T.R));
  return {Ri, -(Ri * T.t)};
}

inline rigid mul(const rigid &A, const rigid &B) {
  const glm::quat Ra = glm::normalize(A.R);
  const glm::quat Rb = glm::normalize(B.R);
  return {glm::normalize(Ra * Rb), Ra * B.t + A.t};
}

inline glm::vec3 so3_log(glm::quat q) {
  q = glm::normalize(q);
  if (q.w < 0.0f)
    q = -q;
  const float w = glm::clamp(q.w, -1.0f, 1.0f);
  const float xyz_n = std::sqrt(q.x * q.x + q.y * q.y + q.z * q.z);
  if (xyz_n < 1.0e-8f)
    return glm::vec3(0.0f);
  const float theta = 2.0f * std::atan2(xyz_n, w);
  return (theta / xyz_n) * glm::vec3(q.x, q.y, q.z);
}

inline glm::quat so3_exp(const glm::vec3 &w) {
  const float theta = glm::length(w);
  if (theta < 1.0e-8f)
    return glm::quat(1, 0.5f * w.x, 0.5f * w.y, 0.5f * w.z);
  const glm::vec3 a = w / theta;
  const float s = std::sin(0.5f * theta);
  return glm::quat(std::cos(0.5f * theta), a.x * s, a.y * s, a.z * s);
}

inline twist log(const rigid &T) {
  twist xi;
  xi.w = so3_log(T.R);
  const float theta = glm::length(xi.w);
  if (theta < 1.0e-6f) {
    xi.v = T.t;
    return xi;
  }
  const glm::mat3 W = hat(xi.w);
  // V^{-1} = I - 0.5 W + (1/theta^2)(1 - theta/2 * cot(theta/2)) W^2
  const float half = 0.5f * theta;
  const float cot = std::cos(half) / std::max(std::sin(half), 1.0e-8f);
  const float a = (1.0f - half * cot) / (theta * theta);
  const glm::mat3 Vinv =
      glm::mat3(1.0f) - 0.5f * W + a * (W * W);
  xi.v = Vinv * T.t;
  return xi;
}

inline rigid exp(const twist &xi) {
  rigid T;
  T.R = so3_exp(xi.w);
  const float theta = glm::length(xi.w);
  if (theta < 1.0e-6f) {
    T.t = xi.v;
    return T;
  }
  const glm::mat3 W = hat(xi.w);
  const float A = std::sin(theta) / theta;
  const float B = (1.0f - std::cos(theta)) / (theta * theta);
  const float C = (1.0f - A) / (theta * theta);
  const glm::mat3 V =
      glm::mat3(1.0f) + B * W + C * (W * W);
  T.t = V * xi.v;
  return T;
}

inline twist log_relative(const rigid &from, const rigid &to) {
  return log(mul(inverse(from), to));
}

inline rigid apply_twist(const rigid &from, const twist &xi) {
  return mul(from, exp(xi));
}

inline twist scale(const twist &xi, float s) {
  return {xi.w * s, xi.v * s};
}

// Geodesic step: T * exp(α log(T^{-1} T_to))
inline rigid chase(const rigid &from, const rigid &to, float alpha) {
  alpha = std::max(0.0f, std::min(1.0f, alpha));
  if (alpha <= 0.0f)
    return from;
  if (alpha >= 1.0f)
    return to;
  return apply_twist(from, scale(log_relative(from, to), alpha));
}

inline void twist_to_array(const twist &xi, float out[6]) {
  out[0] = xi.w.x;
  out[1] = xi.w.y;
  out[2] = xi.w.z;
  out[3] = xi.v.x;
  out[4] = xi.v.y;
  out[5] = xi.v.z;
}

inline twist twist_from_array(const float in[6]) {
  return {glm::vec3(in[0], in[1], in[2]), glm::vec3(in[3], in[4], in[5])};
}

} // namespace se3
} // namespace vermeer
} // namespace gaudi
