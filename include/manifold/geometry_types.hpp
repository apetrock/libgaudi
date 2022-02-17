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

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "vec_addendum.h"
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

  CTYPE normal() { return cross(p[1] - p[0], p[2] - p[0]).normalize(); }

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

  void expandBy(const T &r) { half = m2::va::max(half, r); }

  void expandBy(const bounding_box<T, CTYPE> &other) {
    CTYPE omin = other.center - other.half;
    CTYPE omax = other.center + other.half;
    CTYPE tmin = this->center - this->half;
    CTYPE tmax = this->center + this->half;

    tmin = m2::va::min(tmin, omin);
    tmax = m2::va::max(tmax, omax);

    this->setBounds(tmin, tmax);
  }

  bool overlap(const bounding_box<T, CTYPE> &other) {
    CTYPE omin = other.center - other.half;
    CTYPE omax = other.center + other.half;

    CTYPE tmin = this->center - this->half;
    CTYPE tmax = this->center + this->half;

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

    return true;
  }
};

template <typename T, typename CTYPE>
bounding_box<T, CTYPE> makeBoxMinMax(const CTYPE &min, const CTYPE &max) {
  bounding_box<T, CTYPE> bb;
  bb.setBounds(min, max);
  return bb;
}

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
    CTYPE N = m2::va::cross(e1, e2);
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
    CTYPE n = m2::va::cross(c10, c20);
    out += n.norm() * 0.5;
    return out;
  }

  T solidAngle(CTYPE pi) {
    CTYPE A = p[0] - pi;
    CTYPE B = p[1] - pi;
    CTYPE C = p[2] - pi;
    T a = m2::va::norm(A);
    T b = m2::va::norm(B);
    T c = m2::va::norm(C);

    A /= a;
    B /= b;
    C /= c;

    T divisor =
        m2::va::dot(A, B) + m2::va::dot(B, C) + m2::va::dot(C, A) + T(1);

    T det = m2::va::determinant(A, B, C);

    if (det == 0)
      return T(0);
    return T(2) * atan2(det, divisor);
  }

  T distanceFrom(CTYPE point) { return distance_from_triangle(p, point); }

  T angle(const triangle &B) const {
    CTYPE NA = this->normal();
    CTYPE NB = B.normal();
    return m2::va::dot(NA, NB);
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
            (cb2 - ca0).norm());
    T d2 = 1.0 / 3.0 *
           ((cb0 - ca2).squaredNorm() + (cb1 - ca0).squaredNorm() +
            (cb2 - ca1).squaredNorm());

    T dmin = min(min(d0, d1), d2);
    return dmin;
  };

  bounding_box<T, CTYPE> bbox() {
    CTYPE min, max;
    this->getExtents(min, max);
    return makeBoxMinMax<T, CTYPE>(min, max);
  }

  void getExtents(CTYPE &min, CTYPE &max) {
    min = p[0], max = p[0];
    min = m2::va::min(min, p[1]);
    min = m2::va::min(min, p[2]);
    max = m2::va::max(max, p[1]);
    max = m2::va::max(max, p[2]);
  }
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

template <typename T, typename CTYPE> struct line {
  CTYPE p[2];
  line(CTYPE p0, CTYPE p1) {
    p[0] = p0;
    p[1] = p1;
  };

  CTYPE &operator[](int i) { return p[i]; }
  CTYPE operator[](int i) const { return p[i]; }

  CTYPE center() { return 0.5 * (p[0], p[1]); }
  /*
  void draw(){
    glBegin(GL_LINE);
    glVertex3d(p[0][0],p[0][1],p[0][2]);
    glVertex3d(p[1][0],p[1][1],p[1][2]);
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

// template <typename T, typename CTYPE, typename PRIMITIVE>
// bool boxIntersect(PRIMITIVE& prim, box<T,CTYPE> & b);

// template<typename T, typename vectype>
// class ray{};

template <typename T> class euclidean_space {
public:
  typedef T double_type;
  typedef T real;
  typedef std::complex<T> complex;

  // always use homogeneous coordinates, provides decent error checking
  typedef Eigen::Matrix<T, 3, 1> coordinate_type;

  typedef line<T, coordinate_type> line_type;
  typedef triangle<T, coordinate_type> triangle_type;
  typedef swept_point<T, coordinate_type> swept_point_type;
  typedef swept_triangle<T, coordinate_type> swept_triangle_type;
  typedef bounding_box<T, coordinate_type> box_type;
  typedef Eigen::Matrix<T, 3, 3> mat3;
  typedef Eigen::Matrix<T, 4, 4> mat4;
  typedef unsigned short ushort;
  typedef unsigned int uint;
  typedef unsigned long ulong;

  typedef double double_t;
  typedef float float_t;

  typedef Eigen::Quaternion<T> quat;
  typedef Eigen::Matrix<T, 2, 1> vec2;
  typedef Eigen::Matrix<T, 3, 1> vec3;
  typedef Eigen::Matrix<T, 4, 1> vec4;

  typedef Eigen::Matrix<uint, 2, 1> uint2;
  typedef Eigen::Matrix<uint, 2, 1> uint4;

  typedef Eigen::Matrix<T, 2, 1> int2;
  typedef Eigen::Matrix<T, 4, 1> int3;
  typedef Eigen::Matrix<T, 4, 1> int4;

  enum class edge_index { MAXINDEX = 0 };

  enum class vertex_index { COORDINATE = 0, MAXINDEX = 1 };

  enum class face_index { NORMAL = 0, CENTER = 1, AREA = 2, MAXINDEX = 3 };

  enum class face_vertex_index { MAXINDEX = 0 };
};

typedef euclidean_space<double> space3;

namespace z {

template <typename T> T zero() { return T(0.0); }
template <typename T> T one() { return T(1.0); }

template <> double zero<double>() { return 0.0; }
template <> double one<double>() { return 1.0; }

template <> Eigen::Matrix<double, 3, 1> zero<Eigen::Matrix<double, 3, 1>>() {
  return Eigen::Matrix<double, 3, 1>(0, 0, 0);
}
template <> Eigen::Matrix<double, 3, 1> one<Eigen::Matrix<double, 3, 1>>() {
  return Eigen::Matrix<double, 3, 1>(1, 1, 1);
}

} // namespace z
#endif
