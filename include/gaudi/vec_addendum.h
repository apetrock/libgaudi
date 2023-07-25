/*
 *  vec_addendum.h
 *  Manifold
 *
 *  Created by John Delaney on 3/23/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>
#include <array>
#include <cmath>
#include <complex>
#include <iostream>
#include <math.h>
#include <vector>
#ifndef __VECADD__
#define __VECADD__

using namespace std;

namespace va {
template <typename T> using VEC3 = Eigen::Matrix<T, 3, 1>;
template <typename T> using VEC4 = Eigen::Matrix<T, 4, 1>;
template <typename T> using QUAT = Eigen::Quaternion<T>;

template <int N, typename T> using VEC = Eigen::Matrix<T, N, 1>;

template <typename T>
using MAT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template <typename T> using MAT3 = Eigen::Matrix<T, 3, 3>;

template <typename T> T sgn(T val) { return T(T(0) < val) - (val < T(0)); }
template <typename T> T sgn(const VEC3<T> &v0, const VEC3<T> &v1) {
  return sgn(v0.dot(v1));
}

template <class Tf, class Tv> inline Tv linear(Tf f, const Tv &x, const Tv &y) {
  return (y - x) * f + x;
}

template <class Tf> inline bool mix(Tf f, const bool &x0, const bool &x1) {
  if (x0 && x1)
    return true;
  else if (!x0 && !x1)
    return false;
  if (f < 0.5)
    return x0;
  else
    return x1;
}

template <class Tf, class Tv> inline Tv mix(Tf f, const Tv &x0, const Tv &x1) {
  return (1.0 - f) * x0 + f * x1;
}

template <class T> inline T clamp(const T &v, const T &lb, const T &ub) {
  return std::max(lb, std::min(v, ub));
}

template <typename T>
VEC3<T> catmull_rom(const VEC3<T> &p0, const VEC3<T> &p1, const VEC3<T> &p2,
                    const VEC3<T> &p3, T t /* between 0 and 1 */,
                    T alpha = .5f /* between 0 and 1 */) {
  T tension = 0.0;
  T t0 = 0.0f;
  T t01 = pow((p1 - p0).norm(), alpha);
  T t12 = pow((p2 - p1).norm(), alpha);
  T t23 = pow((p3 - p2).norm(), alpha);

  VEC3<T> m1 = (1.0f - tension) *
               (p2 - p1 + t12 * ((p1 - p0) / t01 - (p2 - p0) / (t01 + t12)));
  VEC3<T> m2 = (1.0f - tension) *
               (p2 - p1 + t12 * ((p3 - p2) / t23 - (p3 - p1) / (t12 + t23)));
  VEC3<T> a = 2.0f * (p1 - p2) + m1 + m2;
  VEC3<T> b = -3.0f * (p1 - p2) - m1 - m1 - m2;
  VEC3<T> c = m1;
  VEC3<T> d = p1;
  return a * t * t * t + b * t * t + c * t + d;
}

template <int N, typename T> struct vecComp {
  bool operator()(const T &x, const T &y) { return ((*x)[N] < (*y)[N]); }
};

// template <typename T>
// inline Eigen::Matrix<T,3,1> cross(const Eigen::Matrix<T,3,1>& vecA, const
// Eigen::Matrix<T,3,1> vecB){ 	Eigen::Matrix<T,3,1> out; 	out[0] =
// vecA[1]*vecB[2]
//- vecA[2]*vecB[1]; 	out[1] = vecA[2]*vecB[0] - vecA[0]*vecB[2]; 	out[2] =
// vecA[0]*vecB[1] - vecA[1]*vecB[0]; 	return out;
//}

template <typename T>
inline VEC3<T> ray_point_intersect(const VEC3<T> &l0, const VEC3<T> &l1,
                                   const VEC3<T> &p0, T &dist) {
  // calculate the intersection of a ray through to a perpendicular plane
  // passing through a point This is for detecting if a point lies in a box or
  // not.
  VEC3<T> l = (l1 - l0);
  l.normalize();
  T den = dot(l, l);
  VEC3<T> vect = p0 - l0;
  T d = dot(vect, l) / den;
  VEC3<T> p_out = l0 + l * d;
  dist = d;
  return p_out;
}

// functions-------------------------------
template <typename T>
inline VEC3<T> line_plane_intersect(const VEC3<T> &r0, const VEC3<T> &r1,
                                    const VEC3<T> &v0, const VEC3<T> &v1,
                                    const VEC3<T> &v2) {
  // calculate the intersection between a line and a plane
  VEC3<T> N = cross((v1 - v0), (v2 - v0));
  VEC3<T> n = N.normalize();

  VEC3<T> r = (r1 - r0).normalize();
  T d = dot((v0 - r0), n) / dot(r, n);
  VEC3<T> p_out = r0 + r * d;
  return p_out;
}

// template <typename T>
// inline Eigen::Matrix<T,3,1> cross(const Eigen::Matrix<T,3,1>& a, const
// Eigen::Matrix<T,3,1>& b){ 	Eigen::Matrix<T,3,1> r; 	r[0] = a[1]*b[2]
// -
// a[2]*b[1]; 	r[1] = a[2]*b[0] - a[0]*b[2]; 	r[2] = a[0]*b[1] -
// a[1]*b[0]; r[3] = 1.0; 	return r;
//}

/*
template <typename T>
const Eigen::Matrix<T,4,1> cross(const Eigen::Matrix<T,4,1>& a,
                                 const Eigen::Matrix<T,4,1>& b){
  return Eigen::Matrix<T,4,1>(a[1]*b[2] - a[2]*b[1],
                              a[2]*b[0] - a[0]*b[2],
                              a[0]*b[1] - a[1]*b[0], 0.0);
}
*/

template <typename T> inline T dist(const VEC3<T> &A, const VEC3<T> &B) {
  return (A - B).norm();
}

template <int N, typename T>
inline T dist(const VEC<N, T> &A, const VEC<N, T> &B) {
  return (A - B).norm();
}

template <int N, typename T>
inline VEC<N, T> max(const VEC<N, T> &A, const VEC<N, T> &B) {
  return A.cwiseMax(B);
}

template <int N, typename T>
inline VEC<N, T> min(const VEC<N, T> &A, const VEC<N, T> &B) {
  return A.cwiseMin(B);
}

template <int N, typename T>
inline VEC<N, T> max(const VEC<N, T> &A, const T &s) {
  return A.cwiseMax(s);
}

template <int N, typename T>
inline VEC<N, T> min(const VEC<N, T> &A, const T &s) {
  return A.cwiseMin(s);
}

template <int N, typename T>
inline T dot(const VEC<N, T> &A, const VEC<N, T> &B) {
  return A.dot(B);
}
template <typename T> inline T dot(const VEC3<T> &A, const VEC3<T> &B) {
  return A.dot(B);
}

inline double norm(double a) { return a; };

template <typename T> T norm(VEC3<T> a) { return a.norm(); };

template <typename T> double norm(std::complex<T> a) {
  a = a * a;
  return a.real() - a.imag();
};

template <typename T> T norm2(VEC3<T> a) {
  return a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
};

template <typename T> VEC3<T> normalize(VEC3<T> a) { return a / norm(a); };
template <typename T> VEC4<T> normalize(VEC4<T> a) { return a / norm(a); };

template <typename T>
inline VEC3<T> cross(const VEC3<T> &vecA, const VEC3<T> &vecB) {
  VEC3<T> out;
  out[0] = vecA[1] * vecB[2] - vecA[2] * vecB[1];
  out[1] = vecA[2] * vecB[0] - vecA[0] * vecB[2];
  out[2] = vecA[0] * vecB[1] - vecA[1] * vecB[0];
  return out;
}

template <typename T>
inline MAT3<T> outer(const VEC3<T> &vecA, const VEC3<T> &vecB) {
  // return vecA.transpose() * vecB;
  return vecA * vecB.transpose();
}

template <typename T>
inline T cotan(const VEC3<T> &c0, const VEC3<T> &c1, const VEC3<T> &c2) {

  VEC3<T> dc10 = c1 - c0;
  VEC3<T> dc20 = c2 - c0;
  // T denom = abs(dc10)*abs(dc20);
  T cosP = dot(dc10, dc20);
  T sinP = norm(cross(dc10, dc20));
  T cotP = cosP / sinP;
  if (sinP > 1e-8)
    return cotP;
  else
    return 0.0;
}

template <typename T>
inline T abs_cos(const VEC3<T> &c0, const VEC3<T> &c1, const VEC3<T> &c2) {

  T e1 = norm(VEC3<T>(c1 - c0));
  T e2 = norm(VEC3<T>(c2 - c0));
  T e3 = norm(VEC3<T>(c2 - c1));

  if (e1 <= 1.0e-8 || e2 <= 1.0e-8)
    return 0.0;

  T cos_alpha = fabs((e1 * e1 + e2 * e2 - e3 * e3) / (2.0f * e1 * e2));

  return cos_alpha;
}

template <typename T>
inline T abs_cotan(const VEC3<T> &c0, const VEC3<T> &c1, const VEC3<T> &c2) {

  T cos_alpha = abs_cos<T>(c0, c1, c2);

  if (cos_alpha > 0.98)
    return 4.0;

  T cotan1 = cos_alpha / sqrt(1.0f - cos_alpha * cos_alpha);
  cotan1 = std::min<T>(cotan1, T(8.0));
  // std::cout << cotan1 << std::endl;
  if (!std::isfinite(cotan1)) {
    return 4.0;
    std::cout << cotan1 << " " << cos_alpha << std::endl;
    std::cout << c0.transpose() << " - " << c1.transpose() << " - "
              << c2.transpose() << std::endl;
    abort();
  }
  return cotan1;
}

template <typename T>
inline bool ray_triangle_intersect(VEC3<T> &pi, const VEC3<T> &r0,
                                   const VEC3<T> &r1, const VEC3<T> &v0,
                                   const VEC3<T> &v1, const VEC3<T> &v2,
                                   T &dist) {
  // Adapted from:
  // Copyright 2001, softSurfer (www.softsurfer.com)
  // using paramatric coordinates, V(s,t) = v0 + su +tv
  // where u = v1-v0
  // and   v = v2-v0
  // if s+t > 1 then we return false
  // we update p with the intersection in the plane

  VEC3<T> u = v1 - v0;
  VEC3<T> v = v2 - v0;

  VEC3<T> n = cross(u, v);
  n = n.normalized();
  T eps = 0.00001;
  VEC3<T> dir = (r1 - r0);
  T b = dot(dir, n);
  if (fabs(b) < eps) {
    // dist = -1;
    return 0;
  }
  dist = (v0 - r0).dot(n) / b;
  pi = r0 + dir * dist;

  VEC3<T> w = pi - v0;

  T duu = dot(u, u);
  T dvv = dot(v, v);
  T duv = dot(u, v);
  T dwu = dot(w, u);
  T dwv = dot(w, v);
  T denom = 1. / (duv * duv - duu * dvv);
  T s0 = (duv * dwv - dvv * dwu) * denom;
  if (s0 < -eps || s0 > 1.0 + eps) {
    // dist = -1;
    return 0; // I is outside T
  }
  T t0 = (duv * dwu - duu * dwv) * denom;
  if (t0 < -eps || (s0 + t0) > 1.0 + eps) {
    // dist = -1;
    return 0; // I is outside T
  }

  return 1; // I is in T
}

template <typename T>
inline bool ray_triangle_intersectII(VEC3<T> &pi, const VEC3<T> &r0,
                                     const VEC3<T> &r1, const VEC3<T> &v0,
                                     const VEC3<T> &v1, const VEC3<T> &v2,
                                     T &dist) {
  // Fast Minimum Storage Ray-Triangle Intersection
  //	Tomas Moller
  //	Ben Trumbore
  // http://www.cs.virginia.edu/~gfx/Courses/2003/ImageSynthesis/papers/Acceleration/Fast%20MinimumStorage%20RayTriangle%20Intersection.pdf

  VEC3<T> e1 = v1 - v0;
  VEC3<T> e2 = v2 - v0;
  VEC3<T> dir = (r1 - r0);

  VEC3<T> p = cross(dir, e2);
  T det = dot(p, e1);
  T eps = 0.0001;
  bool out = 0;
  if (fabs(det) > eps) {

    VEC3<T> tt = r0 - v0;
    VEC3<T> Q = cross(tt, e1);
    T u = dot(p, tt) / det;
    if (u > -eps && u < 1.0 + eps) {
      T v = dot(Q, dir) / det;

      if (v > -eps && u + v < 1.0 + eps) {
        T t = dot(Q, e2) / det;
        dist = t;
        pi = r0 + dir * t;
        out = 1;
      }
    }
  }
  return out;
}

template <typename T>
inline bool point_in_bounds(const VEC3<T> &pt, const T &xmin, const T &xmax,
                            const T &ymin, const T &ymax, const T &zmin,
                            const T &zmax) {
  bool out = true;
  if (pt[0] >= xmax)
    out = false;
  if (pt[0] <= xmin)
    out = false;
  if (pt[1] >= ymax)
    out = false;
  if (pt[1] <= ymin)
    out = false;
  if (pt[2] >= zmax)
    out = false;
  if (pt[2] <= zmin)
    out = false;
  return out;
}

template <typename T, typename POINT_TYPE>
inline bool box_intersection(VEC3<T> cen, VEC3<T> r0, VEC3<T> r1) {
  // this little guy is pretty cool.  So what I'm doing is I'm calculating the
  // plane normal to the ray that intersects with the center of the box.  Then
  // if the point intersecting with that ray and the plane is is outside the
  // box we return false.  This way I don't have to calculate every
  // intersection withe every side.  Bleh, that sounds terrible 	T eps =
  // 0.1; 	T r = o.radius() + eps;
  T r = 0.0;
  T xmax = cen[0] + r, xmin = cen[0] - r;
  T ymax = cen[1] + r, ymin = cen[1] - r;
  T zmax = cen[2] + r, zmin = cen[2] - r;
  bool out = true;

  VEC3<T> d1 = ray_point_intersect(r0, r1, cen);
  out = point_in_bounds(d1, xmin, xmax, ymin, ymax, zmin, zmax);
  return out;
};

template <typename T>
inline T determinant(const VEC3<T> &A, const VEC3<T> &B, const VEC3<T> &C) {
  /*if the volume of the determinant is negative then
    the vectors are clockwise, else anti-clockwise*/
  /*  a b c		A
  // d e f ->	B
  // g h i		C

  detA	= aei + bfg + cdh - afh - bdi - ceg
  ////	= Ax*By*Cz + Ay*Bz*Cx + Az*Bx*Cy
  ////	- Ax*Bz*Cy - Ay*Bx*Cz - Az*By*Cx
  */

  T Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz;

  Ax = A[0];
  Ay = A[1];
  Az = A[2];
  Bx = B[0];
  By = B[1];
  Bz = B[2];
  Cx = C[0];
  Cy = C[1];
  Cz = C[2];

  T out = Ax * By * Cz + Ay * Bz * Cx + Az * Bx * Cy - Ax * Bz * Cy -
          Ay * Bx * Cz - Az * By * Cx;

  return out;
};

template <typename T>
T solidAngle(const VEC3<T> &pi, const VEC3<T> &p0, const VEC3<T> &p1,
             const VEC3<T> &p2) {
  VEC3<T> A = p0 - pi;
  VEC3<T> B = p1 - pi;
  VEC3<T> C = p2 - pi;
  T a = norm(A);
  T b = norm(B);
  T c = norm(C);

  A /= a;
  B /= b;
  C /= c;

  T divisor = dot(A, B) + dot(B, C) + dot(C, A) + T(1);

  T det = determinant(A, B, C);

  if (det == 0)
    return T(0);
  return T(2) * atan2(det, divisor);
}

template <typename T>
T sgn(const VEC3<T> &N0, //
      const VEC3<T> &N1, //
      const VEC3<T> &e) {
  // 2.0 tan(thet/2)
  return sgn(va::determinant(N0, N1, e));
};

template <typename T>
T sin(const VEC3<T> &N0, //
      const VEC3<T> &N1, //
      const VEC3<T> &e) {

  return sgn(N0, N1, e) * norm(VEC3<T>(N1 - N0));
};

template <typename T>
T cos(const VEC3<T> &N0, //
      const VEC3<T> &N1, //
      const VEC3<T> &e) {
  return va::norm(VEC3<T>(N1 + N0));
};

template <typename T>
inline T angle_from_vectors(const VEC3<T> &angA, const VEC3<T> angB) {

  T AB = dot(angA, angB);

  T ABab, magA, magB;
  magA = angA.mag();
  magB = angB.mag();
  ABab = AB / magA / magB;
  T out = acos(ABab);
  return out;
};

template <typename T>
inline VEC3<T> project_on_line(const VEC3<T> &v0, const VEC3<T> &v1,
                               const VEC3<T> &pt, bool clamp = true) {
  VEC3<T> dx = v1 - v0;
  T s = (pt - v0).dot(dx) / dx.dot(dx);
  if (clamp) {
    if (s > 1) {
      s = 1;
    }
    if (s < 0) {
      s = 0;
    }
  }
  VEC3<T> ptl = v0 + s * dx;
  return ptl;
}

template <typename T>
inline T distance_from_line(const VEC3<T> &v0, const VEC3<T> &v1,
                            const VEC3<T> &pt) {
  VEC3<T> ptl = project_on_line(v0, v1, pt);
  T d = (pt - ptl).norm();
  return d;
}

template <typename T>
inline bool contained_in_line(const VEC3<T> &v0, const VEC3<T> &v1,
                              const VEC3<T> &pt) {
  VEC3<T> dx = v1 - v0;
  T s = (pt - v0).dot(dx) / dx.dot(dx);
  if (s > 1 || s < 0) {
    // std::cout << "s out of bounds: " << s << std::endl;
    return false;
  }
  return true;
}

template <typename T>
inline T distance_line_line(const VEC3<T> &x0, const VEC3<T> &x1,
                            const VEC3<T> &x2, const VEC3<T> &x3) {
  VEC3<T> a = x1 - x0;
  VEC3<T> b = x3 - x2;
  VEC3<T> c = x2 - x0;
  VEC3<T> axb = cross(a, b);
  return dot(c, axb) / norm(axb);
}

template <typename T>
inline T distance_from_line(Eigen::Matrix<T, 2, 1> &x1,
                            Eigen::Matrix<T, 2, 1> &x2,
                            Eigen::Matrix<T, 2, 1> &pt) {

  Eigen::Matrix<T, 2, 1> dl = x2 - x1;

  T d =
      ((x2[0] - x1[0]) * (x1[1] - pt[1]) - (x1[0] - pt[0]) * (x2[1] - x1[1])) /
      dl.mag();
  return d;
}

template <typename T>
inline T distance_from_plane(const VEC3<T> &v0, const VEC3<T> &v1,
                             const VEC3<T> &v2, const VEC3<T> &r0) {
  VEC3<T> u = v1 - v0;
  VEC3<T> v = v2 - v0;
  VEC3<T> n = cross(u, v);
  T num = dot(n, r0 - v0);
  T d = num / sqrt(dot(n, n));

  return d;
}

template <typename T>
inline std::array<T, 3>
distance_Segment_Segment(const VEC3<T> &s00, const VEC3<T> &s01,
                         const VEC3<T> &s10, const VEC3<T> &s11,
                         bool clamp_s = true, bool clamp_t = true) {
  // http://geomalgorithms.com/a07-_distance.html#dist3D_Segment_to_Segment()
  //    Input:  two 3D line segments S1 and S2
  //    Return: the shortest distance between S1 and S2
  VEC3<T> u = s01 - s00;
  VEC3<T> v = s11 - s10;
  VEC3<T> w = s00 - s10;
  T a = dot(u, u); // always >= 0
  T b = dot(u, v);
  T c = dot(v, v); // always >= 0
  T d = dot(u, w);
  T e = dot(v, w);
  T D = a * c - b * b; // always >= 0
  T sc, sN, sD = D;    // sc = sN / sD, default sD = D >= 0
  T tc, tN, tD = D;    // tc = tN / tD, default tD = D >= 0

  // compute the line parameters of the two closest points
  if (D < 1e-12) { // the lines are almost parallel
    sN = 0.0;      // force using point P0 on segment S1
    sD = 1.0;      // to prevent possible division by 0.0 later
    tN = e;
    tD = c;
  } else { // get the closest points on the infinite lines
    sN = (b * e - c * d);
    tN = (a * e - b * d);
    if (sN < 0.0 && clamp_s) { // sc < 0 => the s=0 edge is visible
      sN = 0.0;
      tN = e;
      tD = c;
    } else if (sN > sD && clamp_s) { // sc > 1  => the s=1 edge is visible
      sN = sD;
      tN = e + b;
      tD = c;
    }
  }

  if (tN < 0.0 && clamp_t) { // tc < 0 => the t=0 edge is visible
    tN = 0.0;
    // recompute sc for this edge
    if (-d < 0.0 && clamp_s)
      sN = 0.0;
    else if (-d > a && clamp_s)
      sN = sD;
    else {
      sN = -d;
      sD = a;
    }
  } else if (tN > tD && clamp_t) { // tc > 1  => the t=1 edge is visible
    tN = tD;
    // recompute sc for this edge
    if ((-d + b) < 0.0 && clamp_s)
      sN = 0;
    else if ((-d + b) > a && clamp_s)
      sN = sD;
    else {
      sN = (-d + b);
      sD = a;
    }
  }
  // finally do the division to get sc and tc
  sc = (abs(sN) < 1e-12 ? 0.0 : sN / sD);
  tc = (abs(tN) < 1e-12 ? 0.0 : tN / tD);

  // get the difference of the two closest points
  VEC3<T> dP = w + (sc * u) - (tc * v); // =  S1(sc) - S2(tc)
  return {norm2(dP), sc, tc};           // return the closest distance
}

template <typename T>
inline T distance_from_triangle(const std::array<VEC3<T>, 3> tri, VEC3<T> r0,
                                VEC3<T> &pt) {
  // then makes sure its in the direction of the plane.
  const VEC3<T> &v0 = tri[0];
  const VEC3<T> &v1 = tri[1];
  const VEC3<T> &v2 = tri[2];

  VEC3<T> u = v1 - v0;
  VEC3<T> v = v2 - v0;
  VEC3<T> N = cross(u, v);
  T iN2 = 1.0 / (dot(N, N));
  VEC3<T> w = r0 - v0;
  T b10 = dot(cross(u, w), N) * iN2;
  T b20 = dot(cross(w, v), N) * iN2;
  T b12 = 1.0 - b10 - b20;
  b10 = clamp<T>(b10, 0, 1.0);
  b20 = clamp<T>(b20, 0, 1.0);
  b12 = clamp<T>(b12, 0, 1.0);

  pt = b10 * v2 + b20 * v1 + b12 * v0;
  return (r0 - pt).norm();
}

template <typename T>
inline std::array<T, 4> closest_point(const std::array<VEC3<T>, 3> tri,
                                      const VEC3<T> &p) {
  // https://github.com/juj/MathGeoLib/blob/master/src/Geometry/Triangle.cpp
  /** The code for Triangle-float3 test is from Christer Ericson's Real-Time
   * Collision Detection, pp. 141-142. */

  // Check if P is in vertex region outside A.
  VEC3<T> a = tri[0];
  VEC3<T> b = tri[1];
  VEC3<T> c = tri[2];

  VEC3<T> ab = b - a;
  VEC3<T> ac = c - a;
  VEC3<T> ap = p - a;
  T d1 = dot(ab, ap);
  T d2 = dot(ac, ap);
  if (d1 <= 0.f && d2 <= 0.f)
    // Barycentric coordinates are (1,0,0).
    return {(p - a).norm(), 1.0, 0.0, 0.0};

  // Check if P is in vertex region outside B.
  VEC3<T> bp = p - b;
  T d3 = dot(ab, bp);
  T d4 = dot(ac, bp);
  if (d3 >= 0.f && d4 <= d3)
    // Barycentric coordinates are (0,1,0).
    return {(p - b).norm(), 0.0, 1.0, 0.0};

  // Check if P is in edge region of AB, and if so, return the projection of P
  // onto AB.
  T vc = d1 * d4 - d3 * d2;
  if (vc <= 0.f && d1 >= 0.f && d3 <= 0.f) {
    T v = d1 / (d1 - d3);
    // The barycentric coordinates are (1-v, v, 0).
    return {(p - (a + v * ab)).norm(), 1.0 - v, v, 0.0};
  }

  // Check if P is in vertex region outside C.
  VEC3<T> cp = p - c;
  T d5 = dot(ab, cp);
  T d6 = dot(ac, cp);
  if (d6 >= 0.f && d5 <= d6)
    return {(p - c).norm(), 0.0, 0.0, 1.0};

  // Check if P is in edge region of AC, and if so, return the projection of P
  // onto AC.
  T vb = d5 * d2 - d1 * d6;
  if (vb <= 0.f && d2 >= 0.f && d6 <= 0.f) {
    T w = d2 / (d2 - d6);
    // The barycentric coordinates are (1-w, 0, w).
    return {(p - (a + w * ac)).norm(), 1.0 - w, 0.0, w};
  }

  // Check if P is in edge region of BC, and if so, return the projection of P
  // onto BC.
  T va = d3 * d6 - d5 * d4;
  if (va <= 0.f && d4 - d3 >= 0.f && d5 - d6 >= 0.f) {
    T w = (d4 - d3) / (d4 - d3 + d5 - d6);
    // The barycentric coordinates are (0, 1-w, w).
    return {(p - (b + w * (c - b))).norm(), 0.0, 1.0 - w, w};
  }

  // P must be inside the face region. Compute the closest point through its
  // barycentric coordinates (u,v,w).
  T denom = 1.f / (va + vb + vc);
  T v = vb * denom;
  T w = vc * denom;
  return {(p - (a + ab * v + ac * w)).norm(), va * denom, v, w};
}

template <typename T>
VEC3<T> project_to_nullspace(const VEC3<T> &v, const VEC3<T> &n) {
  VEC3<T> u = n.normalized();
  VEC3<T> w = v - v.dot(u) * u;

  return v - w;
}

template <typename T> MAT3<T> rejection_matrix(const VEC3<T> &N) {
  return MAT3<T>::Identity() - N * N.transpose();
}

template <typename T> MAT3<T> projection_matrix(const VEC3<T> &N) {
  return N * N.transpose();
}

template <typename T> MAT3<T> skew_symmetric_matrix(const VEC3<T> &x) {
  MAT3<T> X = MAT3<T>::Zero();
  X.col(0) = VEC3<T>(0, x[2], -x[1]);
  X.col(1) = VEC3<T>(-x[2], 0, x[0]);
  X.col(2) = VEC3<T>(x[1], -x[0], 0);
  return X;
}

template <typename T>
inline VEC3<T> reject(const VEC3<T> &N, const VEC3<T> &A) {
  // N has to be normalized
  T dist = dot(N, A);
  VEC3<T> out = A - dist * N;
  return out;
};

template <typename T>
inline vector<VEC3<T>> reject(const VEC3<T> &norm,
                              const vector<VEC3<T>> &verts) {

  T Nxx, Nxy, Nxz, Nyy, Nyz, Nzz;

  Nxx = norm[0] * norm[0] - 1;
  Nxy = norm[0] * norm[1];
  Nxz = norm[0] * norm[2];
  Nyy = norm[1] * norm[1] - 1;
  Nyz = norm[1] * norm[2];
  Nzz = norm[2] * norm[2] - 1;

  vector<VEC3<T>> out;
  typename vector<VEC3<T>>::iterator itb = verts.begin();
  for (itb; itb != verts.end(); itb++) {
    VEC3<T> cur = *itb;
    VEC3<T> nv;
    nv[0] = Nxx * cur[0] + Nxy * cur[1] + Nxz * cur[2];
    nv[1] = Nxy * cur[0] + Nyy * cur[1] + Nyz * cur[2];
    nv[2] = Nxz * cur[0] + Nyz * cur[1] + Nzz * cur[2];
    out.push_back(nv);
  }

  return out;
}

template <typename T>
inline VEC3<T> reflect(const VEC3<T> &norm, const VEC3<T> &x) {

  T Nxx, Nxy, Nxz, Nyy, Nyz, Nzz;
  VEC3<T> rx;

  Nxx = 1 - 2 * norm[0] * norm[0];
  Nxy = 2 * norm[0] * norm[1];
  Nxz = 2 * norm[0] * norm[2];
  Nyy = 1 - 2 * norm[1] * norm[1];
  Nyz = 2 * norm[1] * norm[2];
  Nzz = 1 - 2 * norm[2] * norm[2];

  rx[0] = Nxx * x[0] + Nxy * x[1] + Nxz * x[2];
  rx[1] = Nxy * x[0] + Nyy * x[1] + Nyz * x[2];
  rx[2] = Nxz * x[0] + Nyz * x[1] + Nzz * x[2];

  return rx;
};

template <typename T> inline VEC3<T> project(const VEC3<T> &N, const T &A) {
  return A;
};

template <typename T>
inline VEC3<T> project(const VEC3<T> &N, const VEC3<T> &A) {
  // N has to be normalized
  T dist = dot(N, A);
  T N2 = dot(N, N);
  VEC3<T> out = (dist / N2) * N;
  return out;
};

// type distance_from_line(Eigen::Matrix<T,3,1> la, Eigen::Matrix<T,3,1> lb,
// Eigen::Matrix<T,3,1> pt){ 	Eigen::Matrix<T,3,1> out; 	T m = lb - la;
// T p = pt -
// la
//
//	return d;
//}

template <typename T>
inline T angle_from_plane(const VEC3<T> &norm, const VEC3<T> &vertA,
                          const VEC3<T> &vertB) {
  T Nxx, Nxy, Nxz, Nyy, Nyz, Nzz;
  VEC3<T> A, B, AC, BC;

  Nxx = norm[0] * norm[0] - 1;
  Nxy = norm[0] * norm[1];
  Nxz = norm[0] * norm[2];
  Nyy = norm[1] * norm[1] - 1;
  Nyz = norm[1] * norm[2];
  Nzz = norm[2] * norm[2] - 1;

  AC = vertA;
  BC = vertB;

  A[0] = Nxx * AC[0] + Nxy * AC[1] + Nxz * AC[2];
  A[1] = Nxy * AC[0] + Nyy * AC[1] + Nyz * AC[2];
  A[2] = Nxz * AC[0] + Nyz * AC[1] + Nzz * AC[2];

  B[0] = Nxx * BC[0] + Nxy * BC[1] + Nxz * BC[2];
  B[1] = Nxy * BC[0] + Nyy * BC[1] + Nyz * BC[2];
  B[2] = Nxz * BC[0] + Nyz * BC[1] + Nzz * BC[2];

  return angle_from_vectors(A, B);
};

template <typename T>
inline T calculate_area(const VEC3<T> &p0, //
                        const VEC3<T> &p1, //
                        const VEC3<T> &p2) {
  return 0.5 * (p1 - p0).cross(p2 - p0).norm();
}

template <typename T>
inline VEC3<T> calculate_normal(const VEC3<T> &p0, //
                                const VEC3<T> &p1, //
                                const VEC3<T> &p2) {
  return (p1 - p0).cross(p2 - p0).normalized();
}

template <typename T> inline VEC3<T> calculate_normal(vector<VEC3<T>> vt_list) {
  // using newell's method
  VEC3<T> tNormal;
  size_t size = vt_list.size();

  for (int i = 0; i < vt_list.size(); i++) {
    VEC3<T> &curr = vt_list[i];
    VEC3<T> &next = vt_list[(i + 1) % size];
    tNormal[0] += (curr[1] - next[1]) * (curr[2] + next[2]);
    tNormal[1] += (curr[2] - next[2]) * (curr[0] + next[0]);
    tNormal[2] += (curr[0] - next[0]) * (curr[1] + next[1]);
  }
  tNormal.normalize();
  return tNormal;
}

template <typename T>
inline VEC3<T> calculate_average(vector<VEC3<T>> vt_list) {
  VEC3<T> tCenter;
  tCenter.zero();
  T size = (T)vt_list.size();
  typename vector<VEC3<T>>::iterator itb = vt_list.begin();
  typename vector<VEC3<T>>::iterator ite = vt_list.end();
  tCenter.zero();
  VEC3<T> coord;
  T inv = 1 / size;
  while (itb != ite) {
    coord = (*itb);
    tCenter += coord * inv;
    ++itb;
  }
  return tCenter;
}

template <typename T>
inline VEC3<T> fit_sphere(const vector<VEC3<T>> &X, vector<VEC3<T>> &Xc, T &r) {

  T N = T(X.size());
  T Sx = 0.0, Sy = 0.0, Sz = 0.0,         //
      Sxx = 0.0, Syy = 0.0, Szz = 0.0,    //
      Sxy = 0.0, Sxz = 0.0, Syz = 0.0,    //
      Sxxx = 0.0, Syyy = 0.0, Szzz = 0.0, //
      Sxyy = 0.0, Sxzz = 0.0, Sxxy = 0.0, //
      Sxxz = 0.0, Syyz = 0.0, Syzz = 0.0;

  for (auto xi : X) {
    T x = xi[0];
    T y = xi[1];
    T z = xi[2];
    Sxx += x * x, Syy += y * y, Szz += z * z;
    Sxy += x * y, Sxz += x * z, Syz += y * z;

    Sxxx += x * x * x, Syyy += y * y * y, Szzz += z * z * z;
    Sxyy += x * y * y, Sxzz += x * z * z, Sxxy += x * x * y;
    Sxxz += x * x * z, Syyz += y * y * z, Syzz += y * z * z;
  }

  T A1 = Sxx + Syy + Szz;
  T a = 2 * Sx * Sx - 2 * N * Sxx;
  T b = 2 * Sx * Sy - 2 * N * Sxy;
  T c = 2 * Sx * Sz - 2 * N * Sxz;
  T d = -N * (Sxxx + Sxyy + Sxzz) + A1 * Sx;

  T e = 2 * Sx * Sy - 2 * N * Sxy;
  T f = 2 * Sy * Sy - 2 * N * Syy;
  T g = 2 * Sy * Sz - 2 * N * Syz;
  T h = -N * (Sxxy + Syyy + Syzz) + A1 * Sy;
  T j = 2 * Sx * Sz - 2 * N * Sxz;
  T k = 2 * Sy * Sz - 2 * N * Syz;
  T l = 2 * Sz * Sz - 2 * N * Szz;
  T m = -N * (Sxxz + Syyz + Szzz) + A1 * Sz;
  T delta = a * (f * l - g * k) - e * (b * l - c * k) + j * (b * g - c * f);

  T xc =
      (d * (f * l - g * k) - h * (b * l - c * k) + m * (b * g - c * f)) / delta;
  T yc =
      (a * (h * l - m * g) - e * (d * l - m * c) + j * (d * g - h * c)) / delta;
  T zc =
      (a * (f * m - h * k) - e * (b * m - d * k) + j * (b * h - d * f)) / delta;
  T R = sqrt(xc ^ 2 + yc ^ 2 + zc ^
             2 + (A1 - 2 * (xc * Sx + yc * Sy + zc * Sz)) / N);
}

template <typename T>
void estimate_3D_circle(VEC3<T> p1, VEC3<T> p2, VEC3<T> p3, VEC3<T> &c,
                        T &radius) {
  // https: // github.com/sergarrido/random/blob/master/circle3d/circle3d.cpp
  VEC3<T> v1 = p2 - p1;
  VEC3<T> v2 = p3 - p1;
  T v1v1, v2v2, v1v2;
  v1v1 = v1.dot(v1);
  v2v2 = v2.dot(v2);
  v1v2 = v1.dot(v2);

  T base = 0.5 / (v1v1 * v2v2 - v1v2 * v1v2);
  T k1 = base * v2v2 * (v1v1 - v1v2);
  T k2 = base * v1v1 * (v2v2 - v1v2);
  c = p1 + v1 * k1 + v2 * k2; // center

  radius = (c - p1).norm();
}

inline float Q_rsqrt(float number) {
  long i;
  float x2, y;
  const float threehalfs = 1.5F;

  x2 = number * 0.5F;
  y = number;
  i = *(long *)&y;           // evil floating point bit level hacking [sic]
  i = 0x5f3759df - (i >> 1); // what the fuck? [sic]
  y = *(float *)&i;
  y = y * (threehalfs - (x2 * y * y)); // 1st iteration
  y = y * (threehalfs - (x2 * y * y)); // 2nd iteration, this can be removed
  return y;
}

template <typename T> inline VEC3<T> fnormalize(VEC3<T> vecin) {
  T mag = vecin[0] * vecin[0] + vecin[1] * vecin[1] + vecin[2] * vecin[2];
  T imag = Q_rsqrt(mag);
  return vecin * imag;
}

template <typename T>
inline VEC3<T> calc_bary(VEC3<T> c, std::vector<VEC3<T>> vertices) {
  MAT3<T> V;
  int i = 0;

  for (auto v : vertices) {
    if (i > 2)
      break;
    V.block(0, i, 3, 1) = v;
    i++;
  }

  // VEC3<T> l = V.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(c );
  VEC3<T> l = V.colPivHouseholderQr().solve(c);
  // std::cout << l.transpose() << std::endl;
  return l;
}

template <typename T> bool greater_than(const VEC3<T> &A, const VEC3<T> &B) {
  return (A.array() > B.array()).sum() > 0;
};

template <typename T> bool less_than(const VEC3<T> &A, const VEC3<T> &B) {
  return (A.array() < B.array()).sum() > 0;
};

template <typename T> QUAT<T> slerp0(QUAT<T> a, QUAT<T> b, float t) {
  T dotAB = a.dot(b);
  dotAB = abs(dotAB);
  T s0, s1;
  if (dotAB < 0.9995f) {
    T s = sqrt(1.0 - dotAB * dotAB); //   Sine of relative angle
    T a = atan2(s, dotAB);
    T c = std::cos(t * a);

    s1 = sqrtf(1 - c * c) / s;
    s0 = c - dotAB * s1;
  } else {
    s0 = 1.0f - t;
    s1 = t;
  }

  return QUAT<T>(s0 * a.coeffs() + s1 * b.coeffs());
}

template <typename T> QUAT<T> slerp(QUAT<T> a, QUAT<T> b, float t) {
  T dotAB = a.dot(b);
  dotAB = abs(dotAB);

  if (dotAB < 0.0) {
    b.coeffs() *= -1.0;
    dotAB *= 1.0;
  }

  // b.coeffs() *= sgn(dotAB);
  // dotAB *= sgn(dotAB);
  if (dotAB < 0.9995f) {
    // T s_ab = sqrt(1.0 - dotAB * dotAB); //   Sine of relative angle
    // T theta = atan2(s_ab, dotAB);

    T theta = acosf(dotAB);

    T sinTheta = sinf(theta);
    T af = sinf((1.0f - t) * theta) / sinTheta;
    T bf = sinf(t * theta) / sinTheta;
    return QUAT<T>(af * a.coeffs() + bf * b.coeffs());
  } else {
    T af = 1.0f - t;
    T bf = t;
    return QUAT<T>(af * a.coeffs() + bf * b.coeffs());
  }
}

template <typename T> QUAT<T> Reciprocal(QUAT<T> q) {
  return QUAT<T>(q.conjugate().coeffs() * (1.0f / q.dot(q)));
}

template <typename T> QUAT<T> Exp(QUAT<T> q) {
  double b = q.vec().norm();

  if (fabs(b) <= 1.0e-14 * abs(q.w()))
    return QUAT<T>(exp(q.w()), 0.0f, 0.0f, 0.0f);
  else {
    T e = exp(q.w());
    T f = std::sin(b) / b;

    return QUAT<T>(e * std::cos(b), //
                   e * f * q.x(),   //
                   e * f * q.y(),   //
                   e * f * q.z());
  }
}

template <typename T> QUAT<T> slerp_far(QUAT<T> q0, QUAT<T> q1, float t) {
  QUAT<T> q01 = Reciprocal(q0) * q1;
  QUAT<T> Vdash(0.0, q01.x(), q01.y(), q01.z());
  QUAT<T> V = QUAT<T>(Vdash.coeffs() * (1.0f / Vdash.norm()));
  T theta = 2.0f * atan2f(Vdash.norm(), q01.w());
  T CompTheta = theta - 2.0f * M_PI;

  return q1 * Exp(QUAT<T>(t * CompTheta * V.coeffs() * 0.5f));
}

} // namespace va
#endif
