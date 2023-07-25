/*
 *  geometry_types.hpp
 *  Manifold
 *
 *  Created by John Delaney on 8/10/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __GEOMETRY_TYPES__
#define __GEOMETRY_TYPES__

#include "gaudi/common.h"
#include "gaudi/vec_addendum.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <iostream>

template <typename T, typename CTYPE> struct swept_triangle {
  static const int size = 3;
  CTYPE p[3];
  CTYPE v[3];
  T dt;
  swept_triangle(){};

  swept_triangle(CTYPE p0, CTYPE p1, CTYPE p2, CTYPE v0, CTYPE v1, CTYPE v2,
                 T Dt) {
    p[0] = p0;
    p[1] = p1;
    p[2] = p2;
    v[0] = v0;
    v[1] = v1;
    v[2] = v2;
  };

  CTYPE &operator[](int i) { return p[i]; }
  CTYPE operator[](int i) const { return p[i]; }

  CTYPE normal() { return va::calculate_normal(p[0], p[1], p[2]); }

  CTYPE center() { return 0.33333 * (p[0] + p[1] + p[2]); }

  void getExtents(CTYPE &min, CTYPE &max) {
    min = this->center();
    max = min;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        min[i] = p[j][i] < min[i] ? p[j][i] : min[i];
        max[i] = p[j][i] > max[i] ? p[j][i] : max[i];
      }
      for (int j = 0; j < 3; j++) {
        min[i] =
            p[j][i] + dt * v[j][i] < min[i] ? p[j][i] + dt * v[j][i] : min[i];
        max[i] =
            p[j][i] + dt * v[j][i] > max[i] ? p[j][i] + dt * v[j][i] : max[i];
      }
    }
  }

  /*
  void draw(){
    glBegin(GL_POLYGON);
    CTYPE norm = this->normal();
    glNormal3f(norm[0], norm[1], norm[2]);
    glVertex3d(p[0][0],p[0][1],p[0][2]);
    glVertex3d(p[1][0],p[1][1],p[1][2]);
    glVertex3d(p[2][0],p[2][1],p[2][2]);
    glEnd();
  }
  */
};

template <typename T, typename CTYPE> struct bounding_box {
public:
  CTYPE center;
  CTYPE half;
  bounding_box(){};
  bounding_box(const CTYPE &cen, const CTYPE &h) : center(cen), half(h){};

  void setBounds(const CTYPE &min, const CTYPE &max) {
    // should probably assert that min < max
    center = 0.5 * (min + max);
    half = 0.5 * (max - min);
  }

  void inflate(const CTYPE &tol) { half += tol; }

  void expandBy(const T &r) { half = va::max(half, r); }

  void expandBy(const bounding_box<T, CTYPE> &other) {
    CTYPE omin = other.center - other.half;
    CTYPE omax = other.center + other.half;
    CTYPE tmin = this->center - this->half;
    CTYPE tmax = this->center + this->half;

    tmin = va::min(tmin, omin);
    tmax = va::max(tmax, omax);

    this->setBounds(tmin, tmax);
  }

  bool overlap(const bounding_box<T, CTYPE> &other) const {
    CTYPE omin = other.center - other.half;
    CTYPE omax = other.center + other.half;

    CTYPE tmin = this->center - this->half;
    CTYPE tmax = this->center + this->half;

    if (va::greater_than(omin, tmax))
      return false;
    if (va::less_than(omax, tmin))
      return false;
    /*
    if (omin[0] > tmax[0])
      return false;
    else if (omax[0] < tmin[0])
      return false;

    else if (omin[1] > tmax[1])
      return false;
    else if (omax[1] < tmin[1])
      return false;

    else if (omin[2] > tmax[2])
      return false;
    else if (omax[2] < tmin[2])
      return false;
*/
    return true;
  }
};

template <typename T, typename CTYPE>
bounding_box<T, CTYPE> makeBoxMinMax(const CTYPE &min, const CTYPE &max) {
  bounding_box<T, CTYPE> bb;
  bb.setBounds(min, max);
  return bb;
}

template <typename T, typename CTYPE> struct line {

public:
  static const int size = 3;
  CTYPE p[2];
  bounding_box<T, CTYPE> _bbox;
  bool _bbox_init = false;

  line(){};

  line(const CTYPE &p0, const CTYPE &p1) { setP(p0, p1); };

  void setP(const CTYPE &p0, const CTYPE &p1) {
    p[0] = p0;
    p[1] = p1;
    update_bbox();
  }

  CTYPE &operator[](int i) { return p[i]; }
  CTYPE operator[](int i) const { return p[i]; }

  CTYPE center() const { return 0.5 * (p[0] + p[1]); }

  T distanceFrom(CTYPE point) { return distance_from_line(p, point); }

  T avgSqdDist(const line &B) const {

    const CTYPE &ca0 = this->p[0];
    const CTYPE &ca1 = this->p[1];
    const CTYPE &cb0 = B.p[0];
    const CTYPE &cb1 = B.p[1];
    T d0 = 1.0 / 2.0 * ((cb0 - ca0).squaredNorm() + (cb1 - ca1).squaredNorm());
    T d1 = 1.0 / 2.0 * ((cb0 - ca1).squaredNorm() + (cb1 - ca0).squaredNorm());

    return min(d0, d1);
  };

  void update_bbox() {
    CTYPE min, max;
    this->getExtents(min, max);
    _bbox = makeBoxMinMax<T, CTYPE>(min, max);
  }

  const bounding_box<T, CTYPE> &bbox() const { return _bbox; }

  void getExtents(CTYPE &min, CTYPE &max) {
    min = va::min(p[0], p[1]);
    max = va::max(p[0], p[1]);
  }
};

template <typename T, typename CTYPE> struct triangle {

public:
  static const int size = 3;
  CTYPE p[3];

  triangle(){};

  triangle(CTYPE p0, CTYPE p1, CTYPE p2) {
    p[0] = p0;
    p[1] = p1;
    p[2] = p2;
  };

  CTYPE &operator[](int i) { return p[i]; }
  CTYPE operator[](int i) const { return p[i]; }

  CTYPE normal() const {

    CTYPE e1 = p[1] - p[0];
    CTYPE e2 = p[2] - p[0];
    CTYPE N = va::cross(e1, e2);
    N.normalize();
    return N;
  }

  CTYPE center() const { return 0.33333 * (p[0] + p[1] + p[2]); }

  T area() const {
    T out = 0;
    CTYPE c0 = p[0];
    CTYPE c1 = p[1];
    CTYPE c2 = p[2];
    CTYPE c10 = c1 - c0;
    CTYPE c20 = c2 - c0;
    CTYPE n = va::cross(c10, c20);
    out += n.norm() * 0.5;
    return out;
  }

  bounding_box<T, CTYPE> bbox() {
    CTYPE min, max;
    this->getExtents(min, max);
    return makeBoxMinMax<T, CTYPE>(min, max);
  }

  void getExtents(CTYPE &min, CTYPE &max) {
    min = p[0], max = p[0];
    min = va::min(min, p[1]);
    min = va::min(min, p[2]);
    max = va::max(max, p[1]);
    max = va::max(max, p[2]);
  }

  T solidAngle(CTYPE pi) {
    CTYPE A = p[0] - pi;
    CTYPE B = p[1] - pi;
    CTYPE C = p[2] - pi;
    T a = va::norm(A);
    T b = va::norm(B);
    T c = va::norm(C);

    A /= a;
    B /= b;
    C /= c;

    T divisor = va::dot(A, B) + va::dot(B, C) + va::dot(C, A) + T(1);

    T det = va::determinant(A, B, C);

    if (det == 0)
      return T(0);
    return T(2) * atan2(det, divisor);
  }

  T determinant(CTYPE pi) {
    CTYPE A = p[0] - pi;
    CTYPE B = p[1] - pi;
    CTYPE C = p[2] - pi;
    /*
        T a = va::norm(A);
        T b = va::norm(B);
        T c = va::norm(C);

        A /= a;
        B /= b;
        C /= c;
    */
    return va::determinant(A, B, C);
  }

  T distanceFrom(CTYPE point) { return va::distance_from_triangle(p, point); }

  bool rayIntersect(const CTYPE &r0, const CTYPE &r1, T &d) {
    CTYPE pi;
    return va::ray_triangle_intersect<T>(pi, r0, r1, p[0], p[1], p[2], d);
  }

  T angle(const triangle &B) const {
    CTYPE NA = this->normal();
    CTYPE NB = B.normal();
    return va::dot(NA, NB);
  }

  T avgSqdDist(const triangle &B) const {

    const CTYPE &ca0 = this->p[0];
    const CTYPE &ca1 = this->p[1];
    const CTYPE &ca2 = this->p[2];
    // assume the other coordinate rotates in opposite direction since they
    // are facing opposite directions
    const CTYPE &cb0 = B.p[2];
    const CTYPE &cb1 = B.p[1];
    const CTYPE &cb2 = B.p[0];
    T d0 = 1.0 / 3.0 *
           ((cb0 - ca0).squaredNorm() + (cb1 - ca1).squaredNorm() +
            (cb2 - ca2).squaredNorm());
    T d1 = 1.0 / 3.0 *
           ((cb0 - ca1).squaredNorm() + (cb1 - ca2).squaredNorm() +
            (cb2 - ca0).squaredNorm());
    T d2 = 1.0 / 3.0 *
           ((cb0 - ca2).squaredNorm() + (cb1 - ca0).squaredNorm() +
            (cb2 - ca1).squaredNorm());

    T dmin = min(min(d0, d1), d2);
    return dmin;
  };
};

template <typename T, typename CTYPE> struct swept_point {
public:
  CTYPE p;
  CTYPE v;
  T dt;
  swept_point(CTYPE p0, CTYPE vi) {
    p = p0;
    v = vi;
  };

  CTYPE center() { return 0.5 * (p + v); }

  void getExtents(CTYPE &min, CTYPE &max) {
    min = this->center();
    max = min;
    for (int i = 0; i < 3; i++) {
      min[i] = p[i] < min[i] ? p[i] : min[i];
      max[i] = p[i] > max[i] ? p[i] : max[i];
      min[i] = p[i] + dt * v[i] < min[i] ? p[i] + dt * v[i] : min[i];
      max[i] = p[i] + dt * v[i] > max[i] ? p[i] + dt * v[i] : max[i];
    }
  }
  /*
  void draw(){
    glBegin(GL_LINE);
    glVertex3d(p[0],p[1],p[2]);
    glVertex3d(p[0] + dt*v[0],
               p[1] + dt*v[1],
               p[2] + dt*v[2]);
    glEnd();
    }
*/
};

template <typename T, typename CTYPE> struct box {
  CTYPE center;
  CTYPE half;

  box(){};

  box(CTYPE cen, CTYPE h) {
    center = cen;
    half = h;
  };
};

template <typename T, typename CTYPE> struct ls_sphere {
  T //
      Sx = 0.0,
      Sy = 0.0, Sz = 0.0,                 //
      Sxx = 0.0, Syy = 0.0, Szz = 0.0,    //
      Sxy = 0.0, Sxz = 0.0, Syz = 0.0,    //
      Sxxx = 0.0, Syyy = 0.0, Szzz = 0.0, //
      Sxyy = 0.0, Sxzz = 0.0, Sxxy = 0.0, //
      Sxxz = 0.0, Syyz = 0.0, Syzz = 0.0;
  T W = 0;

  void add_vec(CTYPE xi, T w) {
    T x = xi[0];
    T y = xi[1];
    T z = xi[2];
    Sxx += x * x, Syy += y * y, Szz += z * z;
    Sxy += x * y, Sxz += x * z, Syz += y * z;

    Sxxx += x * x * x, Syyy += y * y * y, Szzz += z * z * z;
    Sxyy += x * y * y, Sxzz += x * z * z, Sxxy += x * x * y;
    Sxxz += x * x * z, Syyz += y * y * z, Syzz += y * z * z;
    W += w;
  }

  ls_sphere<T, CTYPE> &operator+=(const ls_sphere<T, CTYPE> &rhs) {
    Sxx += rhs.Sxx, Syy += rhs.Syy, Szz += rhs.Szz;
    Sxy += rhs.Sxy, Sxz += rhs.Sxz, Syz += rhs.Syz;

    Sxxx += rhs.Sxxx, Syyy += rhs.Syyy, Szzz += rhs.Szzz;
    Sxyy += rhs.Sxyy, Sxzz += rhs.Sxzz, Sxxy += rhs.Sxxy;
    Sxxz += rhs.Sxxz, Syyz += rhs.Syyz, Syzz += rhs.Syzz;
    W += rhs.W;

    return *this;
  }

  void calc(CTYPE &X, T &R) {
    T A1 = Sxx + Syy + Szz;
    T a = 2.0 * Sx * Sx - 2.0 * W * Sxx;
    T b = 2.0 * Sx * Sy - 2.0 * W * Sxy;
    T c = 2.0 * Sx * Sz - 2.0 * W * Sxz;
    T d = -W * (Sxxx + Sxyy + Sxzz) + A1 * Sx;

    T e = 2.0 * Sx * Sy - 2.0 * W * Sxy;
    T f = 2.0 * Sy * Sy - 2.0 * W * Syy;
    T g = 2.0 * Sy * Sz - 2.0 * W * Syz;
    T h = -W * (Sxxy + Syyy + Syzz) + A1 * Sy;
    T j = 2.0 * Sx * Sz - 2.0 * W * Sxz;
    T k = 2.0 * Sy * Sz - 2.0 * W * Syz;
    T l = 2.0 * Sz * Sz - 2.0 * W * Szz;
    T m = -W * (Sxxz + Syyz + Szzz) + A1 * Sz;

    T delta = a * (f * l - g * k) - e * (b * l - c * k) + j * (b * g - c * f);
    X[0] = (d * (f * l - g * k) - h * (b * l - c * k) + m * (b * g - c * f)) /
           delta;
    X[1] = (a * (h * l - m * g) - e * (d * l - m * c) + j * (d * g - h * c)) /
           delta;
    X[2] = (a * (f * m - h * k) - e * (b * m - d * k) + j * (b * h - d * f)) /
           delta;

    R = sqrt(X[0] * X[0] + X[1] * X[1] + X[2] * X[2] +
             (A1 - 2.0 * (X[0] * Sx + X[1] * Sy + X[2] * Sz)) / W);
  }
};

namespace gaudi {
namespace z {

template <typename T> T zero() { return T(0.0); }
template <typename T> T one() { return T(1.0); }

template <> inline double zero<double>() { return 0.0; }
template <> inline double one<double>() { return 1.0; }

template <>
inline Eigen::Matrix<double, 3, 1> zero<Eigen::Matrix<double, 3, 1>>() {
  return Eigen::Matrix<double, 3, 1>(0, 0, 0);
}
template <>
inline Eigen::Matrix<double, 3, 1> one<Eigen::Matrix<double, 3, 1>>() {
  return Eigen::Matrix<double, 3, 1>(1, 1, 1);
}

template <>
inline Eigen::Matrix<double, 3, 3> zero<Eigen::Matrix<double, 3, 3>>() {
  Eigen::Matrix<double, 3, 3> m = Eigen::Matrix<double, 3, 3>::Zero();
  // m << 0, 0, 0, 0, 0, 0, 0, 0, 0;
  return m;
}

template <>
inline Eigen::Matrix<double, 3, 3> one<Eigen::Matrix<double, 3, 3>>() {
  Eigen::Matrix<double, 3, 3> m;
  m << 1, 0, 0, 0, 1, 0, 0, 0, 1;
  return m;
}

template <>
inline Eigen::Matrix<double, 4, 3> zero<Eigen::Matrix<double, 4, 3>>() {
  // Eigen::Matrix<double, 4, 3> m;
  Eigen::Matrix<double, 4, 3> m = Eigen::Matrix<double, 4, 3>::Zero();
  // m << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
  return m;
}
} // namespace z

namespace ext {
typedef std::array<vec3, 2> extents_t;

bool overlap(const extents_t &A, const extents_t &B) {
  if (va::greater_than(A[0], B[1]))
    return false;
  if (va::less_than(A[1], B[0]))
    return false;

  return true;
}

extents_t inflate(extents_t e, real eps) {
  vec3 deps(eps, eps, eps);
  e[0] -= deps;
  e[1] += deps;
  return e;
}

extents_t expand(const extents_t &e, const vec3 &x) {
  extents_t eout;
  eout[0] = va::min(e[0], x);
  eout[1] = va::max(e[1], x);
  return eout;
}

extents_t expand(const extents_t &eA, const extents_t &eB) {
  extents_t eout;
  eout = expand(eA, eB[0]);
  eout = expand(eout, eB[1]);
  return eout;
}

real distance(const extents_t &e, const vec3 &x) {
  extents_t eout;
  vec3 c = 0.5 * (e[0] + e[1]);
  return (x - c).norm();
}
} // namespace ext
} // namespace gaudi
#endif
