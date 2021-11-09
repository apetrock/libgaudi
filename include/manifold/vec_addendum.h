/*
 *  vec_addendum.h
 *  Manifold
 *
 *  Created by John Delaney on 3/23/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
#include <iostream>
#include <math.h>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#ifndef __VECADD__
#define __VECADD__

using namespace std;

namespace m2 {
namespace va {
template <typename T> using VEC3 = Eigen::Matrix<T, 3, 1>;
template <typename T> using VEC4 = Eigen::Matrix<T, 4, 1>;
template <int N, typename T> using VEC = Eigen::Matrix<T, N, 1>;

template <class Tf, class Tv> inline Tv linear(Tf f, const Tv &x, const Tv &y) {
  return (y - x) * f + x;
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
inline VEC3<T> ray_point_intersect(const VEC3<T> &l0, const VEC3<T> &l1, const VEC3<T> &p0,
                                T &dist) {
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
inline VEC3<T> line_plane_intersect(const VEC3<T> &r0, const VEC3<T> &r1, const VEC3<T> &v0,
                                 const VEC3<T> &v1, const VEC3<T> &v2) {
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
inline T dot(const VEC<N, T> &A, const VEC<N, T> &B) {
  return A.dot(B);
}
template <typename T> inline T dot(const VEC3<T> &A, const VEC3<T> &B) {
  return A.dot(B);
}

template <typename T> T norm(VEC3<T> a) { return a.norm(); };
template <typename T> T norm2(VEC3<T> a) {
  return a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
};

template <typename T> VEC3<T> normalize(VEC3<T> a) { 
  return a / norm(a);
};
template <typename T> VEC4<T> normalize(VEC4<T> a) { return a / norm(a); };

template <typename T> inline VEC3<T> cross(const VEC3<T> &vecA, const VEC3<T> &vecB) {
  VEC3<T> out;
  out[0] = vecA[1] * vecB[2] - vecA[2] * vecB[1];
  out[1] = vecA[2] * vecB[0] - vecA[0] * vecB[2];
  out[2] = vecA[0] * vecB[1] - vecA[1] * vecB[0];
  return out;
}


template <typename T>
inline T cotan(const VEC3<T> &c0, const VEC3<T> &c1, const VEC3<T> &c2) {

  VEC3<T> dc10 = c1 - c0;
  VEC3<T> dc20 = c2 - c0;
  // T denom = abs(dc10)*abs(dc20);
  T cosP = dot(dc10, dc20);
  T sinP = norm(cross(dc10, dc20));
  T cotP = cosP / sinP;
  if (sinP > 1e-12)
    return cotP;
  else
    return 1.0;
}

template <typename T>
inline bool ray_triangle_intersect(VEC3<T> &pi, const VEC3<T> &r0, const VEC3<T> &r1,
                                   const VEC3<T> &v0, const VEC3<T> &v1,
                                   const VEC3<T> &v2, T &dist) {
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
  n = n.normalize();
  T eps = 0.00001;
  VEC3<T> dir = (r1 - r0);
  T b = dot(dir, n);
  if (fabs(b) < eps) {
    return 0;
  }
  dist = dot((v0 - r0), n) / b;
  pi = r0 + dir * dist;

  VEC3<T> w = pi - v0;

  T duu = dot(u, u);
  T dvv = dot(v, v);
  T duv = dot(u, v);
  T dwu = dot(w, u);
  T dwv = dot(w, v);
  T denom = 1. / (duv * duv - duu * dvv);
  T s0 = (duv * dwv - dvv * dwu) * denom;
  if (s0 < -eps || s0 > 1.0 + eps)
    return 0; // I is outside T
  T t0 = (duv * dwu - duu * dwv) * denom;
  if (t0 < -eps || (s0 + t0) > 1.0 + eps)
    return 0; // I is outside T

  return 1; // I is in T
}

template <typename T>
inline bool ray_triangle_intersectII(VEC3<T> &pi, const VEC3<T> &r0, const VEC3<T> &r1,
                                     const VEC3<T> &v0, const VEC3<T> &v1,
                                     const VEC3<T> &v2, T &dist) {
  // Fast Minimum Storage RayTriangle Intersection
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
inline T angle_from_vectors(const VEC3<T> &angA, const VEC3<T> angB) {

  T AB = dot(angA, angB);

  T ABab, magA, magB;
  magA = angA.mag();
  magB = angB.mag();
  ABab = AB / magA / magB;
  T out = acos(ABab);
  return out;
};

template <int N, typename T>
inline bool project_on_line(const Eigen::Matrix<T, N, 1> &x0,
                            const Eigen::Matrix<T, N, 1> &x1,
                            const Eigen::Matrix<T, N, 1> &pt,
                            Eigen::Matrix<T, N, 1> &pr, T eps) {
  Eigen::Matrix<T, N, 1> Q = pt - x0;
  Eigen::Matrix<T, N, 1> dx = x1 - x0;
  T mdx = dx.mag();
  Eigen::Matrix<T, N, 1> ndx = dx / mdx;
  T t = dot(Q, ndx) / mdx;
  pr = x0 + t * dx;
  if (t > -eps && t <= 1 + eps) {
    return true;
  } else {
    return false;
  }
}

template <typename T>
inline T distance_from_line(const VEC3<T> &v0, const VEC3<T> &v1, const VEC3<T> &pt) {
  VEC3<T> dx = v1 - v0;
  T s = dot((pt - v0), dx) / (dot(dx, dx));
  if (s >= 1)
    s = 1;
  if (s <= 0)
    s = 0;
  VEC3<T> ptl = v0 + s * dx;

  T d = norm2(pt - ptl);
  return d;
}

template <typename T>
inline T distance_line_line(const VEC3<T> &x0, const VEC3<T> &x1, const VEC3<T> &x2,
                            const VEC3<T> &x3) {
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
inline T distance_from_plane(const VEC3<T> &v0, const VEC3<T> &v1, const VEC3<T> &v2,
                             const VEC3<T> &r0) {
  VEC3<T> u = v1 - v0;
  VEC3<T> v = v2 - v0;
  VEC3<T> n = cross(u, v);
  T num = dot(n, r0 - v0);
  T d = num / sqrt(dot(n, n));

  return d;
}

template <typename T>
inline T distance_from_triangle(const VEC3<T> *tri, VEC3<T> r0) {
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

  if (b10 >= 0.0 && b20 >= 0.0 && b12 >= 0.0) {
    VEC3<T> c = b10 * v2 + b20 * v1 + b12 * v0;
    return norm2(r0 - c);
  } else {
#if 1
    if (b10 <= 0) {
      return distance_from_line(v0, v1, r0);
    } else if (b20 <= 0) {
      return distance_from_line(v0, v2, r0);
    } else {
      return distance_from_line(v1, v2, r0);
    }
#endif
  }
}

template <typename T>
inline vector<VEC3<T>> orthogonal_project(const VEC3<T> &norm,
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

template <typename T> inline VEC3<T> reflect(const VEC3<T> &norm, const VEC3<T> &x) {

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

template <typename T>
inline VEC3<T> orthogonal_project(const VEC3<T> &N, const VEC3<T> &A) {
  // N has to be normalized
  T dist = dot(N, A);
  VEC3<T> out = A - dist * N;
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

template <typename T> inline VEC3<T> calculate_average(vector<VEC3<T>> vt_list) {
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

template <typename T> inline void radial_test() {
  VEC3<T> A;
  VEC3<T> B;

  T imax = 16;
  for (int i = 0; i < int(imax); i++) {
    A.set(0, 1, 0);
    B.set(sin(2 * 3.14 * T(i) / imax), cos(2 * 3.14 * T(i) / imax), 0);

    T magA = A.mag();
    T magB = B.mag();
    VEC3<T> cAB = cross(A, B);
    T magcAB = cAB.mag();

    T dotAB = dot(A, B);

    T ndotAB = dotAB / magA / magB;
    T nmagcAB = magcAB / magA / magB;

    T theta1 = acos(ndotAB);
    T theta2 = asin(nmagcAB);
    cout << i << ": " << 2 * 3.14 * T(i) / imax << endl;
    cout << "   " << ndotAB << ", " << nmagcAB << endl;
    cout << "   " << theta1 << ", " << theta2 << endl;
  }
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
} // namespace va
} // namespace m2
#endif
