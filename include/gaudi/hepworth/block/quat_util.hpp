#ifndef GAUDI_HEPWORTH_BLOCK_QUAT_UTIL_HPP
#define GAUDI_HEPWORTH_BLOCK_QUAT_UTIL_HPP

/// Quaternions as `vec4` `(w, x, y, z)`.

#include "gaudi/common.h"

#include <cmath>

namespace gaudi {
namespace hepworth {
namespace block {
namespace quat_ops {

inline vec4 qmul(const vec4 &a, const vec4 &b) {
  return vec4(a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
              a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
              a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1],
              a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0]);
}

inline vec4 qconj(const vec4 &a) { return vec4(a[0], -a[1], -a[2], -a[3]); }

/// `v` rotated by `q` normalized.
inline vec3 qrot(const vec4 &q, const vec3 &v) {
  const real l = q.norm();
  const vec4 u = l > real(0) ? vec4(q / l) : vec4(1, 0, 0, 0);
  const vec4 r = qmul(qmul(u, vec4(0, v[0], v[1], v[2])), qconj(u));
  return vec3(r[1], r[2], r[3]);
}

/// The unit quaternion of a rotation matrix.
inline vec4 quat_from_matrix(const mat3 &R) {
  vec4 q;
  const real tr = R.trace();
  if (tr > 0) {
    const real s = real(2) * std::sqrt(tr + 1);
    q = vec4(real(0.25) * s, (R(2, 1) - R(1, 2)) / s, (R(0, 2) - R(2, 0)) / s, (R(1, 0) - R(0, 1)) / s);
  } else {
    int i = 0;
    if (R(1, 1) > R(i, i))
      i = 1;
    if (R(2, 2) > R(i, i))
      i = 2;
    const int j = (i + 1) % 3, k = (i + 2) % 3;
    const real s = real(2) * std::sqrt(1 + R(i, i) - R(j, j) - R(k, k));
    vec3 v;
    v[i] = real(0.25) * s;
    v[j] = (R(j, i) + R(i, j)) / s;
    v[k] = (R(k, i) + R(i, k)) / s;
    q = vec4((R(k, j) - R(j, k)) / s, v[0], v[1], v[2]);
  }
  return q.normalized();
}

} // namespace quat_ops
} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif
