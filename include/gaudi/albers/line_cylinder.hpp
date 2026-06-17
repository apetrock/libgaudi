#ifndef ALBERS_LINE_CYLINDER_H
#define ALBERS_LINE_CYLINDER_H

#include <Eigen/Dense>

#include "gaudi/common.h"

namespace gaudi {
namespace albers {

TYPEDEF_MAT_NM(2, 6)

inline mat3 cross_matrix(const vec3 &x) {
  mat3 X;
  X << 0.0, -x[2], x[1], x[2], 0.0, -x[0], -x[1], x[0], 0.0;
  return X;
}

inline mat26 mk_normal_aligned_line_A(const vec3 &x, const vec3 &N) {
  // Plucker line q = [d, m], m = c x d.
  // s = x x d - m equals radial x d. If radial is parallel to N, then:
  //   d dot N = 0
  //   s dot N = 0
  // Both constraints are linear in [d, m].
  mat26 A = mat26::Zero();
  A.block<1, 3>(0, 0) = N.transpose();
  A.block<1, 3>(1, 0) = (N.transpose() * cross_matrix(x));
  A.block<1, 3>(1, 3) = -N.transpose();
  return A;
}

inline vec3 plucker_line_direction(const vec6 &Q) {
  const vec3 d = Q.segment(0, 3);
  return d.norm() > 1e-12 ? d.normalized() : vec3::Zero();
}

inline vec3 plucker_line_point(const vec6 &Q) {
  const vec3 d = Q.segment(0, 3);
  vec3 m = Q.segment(3, 3);
  const real d2 = d.squaredNorm();
  if (d2 < 1e-24) {
    return vec3::Zero();
  }
  m -= (d.dot(m) / d2) * d;
  return d.cross(m) / d2;
}

inline vec3 plucker_line_closest_point(const vec6 &Q, const vec3 &x) {
  const vec3 d = plucker_line_direction(Q);
  if (d.squaredNorm() < 1e-24) {
    return vec3::Zero();
  }
  const vec3 c = plucker_line_point(Q);
  return c + d * d.dot(x - c);
}

class normal_aligned_line {
public:
  using coefficients = vec6;

  normal_aligned_line() {
    A = mat6::Zero();
    b = vec6::Zero();
  }

  void accumulate(real w, const vec3 &x, const vec3 &N) {
    mat26 Ai = mk_normal_aligned_line_A(x, N);
    A += w * Ai.transpose() * Ai;
  }

  vec6 solve() {
    Eigen::SelfAdjointEigenSolver<mat6> es(A);
    if (es.info() != Eigen::Success) {
      return vec6::Zero();
    }
    vec6 q = es.eigenvectors().col(0);
    if (q.segment(0, 3).norm() < 1e-12) {
      return vec6::Zero();
    }
    return q;
  }

  mat6 A;
  vec6 b;
};

} // namespace albers
} // namespace gaudi

#endif // ALBERS_LINE_CYLINDER_H
