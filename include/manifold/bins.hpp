#ifndef __M2BINS__
#define __M2BINS__

#include "conj_grad.hpp"
#include "m2Includes.h"

#include "quartic/cubic.hpp"
#include "tribox_test.hpp"
#include <stack>

#include "TIMER.h"

#include <cmath>
#include <limits>

bool boxOverlap(space3::triangle_type &tri, space3::box_type &b) {
  m2::va::tri_box<space3::double_type, space3::coordinate_type> tribox;
  bool inBox = tribox.triBoxOverlap(b.center, b.half, tri.p);
  return inBox;
};

namespace m2 {
template <typename SPACE> struct line_tests {
  M2_TYPEDEFS;

  inline bool ray_triangle_intersect(const coordinate_type &r0,
                                     const coordinate_type &r1,
                                     const coordinate_type &v0,
                                     const coordinate_type &v1,
                                     const coordinate_type &v2,
                                     coordinate_type &pi, T &dist) {
    // Adapted from:
    // Copyright 2001, softSurfer (www.softsurfer.com)
    // using paramatric coordinates, V(s,t) = v0 + su +tv
    // where u = v1-v0
    // and   v = v2-v0
    // if s+t > 1 then we return false
    // we update p with the intersection in the plane

    coordinate_type u = v1 - v0;
    coordinate_type v = v2 - v0;

    coordinate_type n = cross(u, v);
    n = n.normalize();
    T eps = 1e-12;
    coordinate_type dir = (r1 - r0);
    T b = dot(dir, n);
    if (fabs(b) < eps)
      return 0;
    dist = dot((v0 - r0), n) / b;
    pi = r0 + dir * dist;

    coordinate_type w = pi - v0;

    T duu = dot(u, u);
    T dvv = dot(v, v);
    T duv = dot(u, v);
    T dwu = dot(w, u);
    T dwv = dot(w, v);
    T denom = (duv * duv - duu * dvv);
    if (fabs(denom) < eps)
      return 0;
    T idenom = 1. / denom;
    T s0 = (duv * dwv - dvv * dwu) * idenom;
    T t0 = (duv * dwu - duu * dwv) * idenom;
    if (s0 < -eps || s0 > 1.0 + eps)
      return 0; // I is outside T
    if (t0 < -eps || (s0 + t0) > 1.0 + eps)
      return 0; // I is outside T
    return 1;   // I is in T
  }

  inline bool lineTriangle(coordinate_type e0, coordinate_type e1,
                           coordinate_type p0, coordinate_type p1,
                           coordinate_type p2) {
    T t;
    coordinate_type p;
    bool hit = ray_triangle_intersect(e0, e1, p0, p1, p2, p, t);
    if (!hit)
      return false;

    T eps = 1e-10;
    if (t < eps)
      return false;
    if (t > 1.0 - eps)
      return false;

    std::ostringstream ss;
    ss << t;
    std::string s(ss.str());
    return true;
  }

  bool rayIntersectBox(coordinate_type ro, coordinate_type rd,
                       coordinate_type boxmin, coordinate_type boxmax, T *tnear,
                       T *tfar) {
    // compute intersection of ray with all six bbox planes
    coordinate_type invR =
        coordinate_type(1.0 / rd[0], 1.0 / rd[1], 1.0 / rd[2]);
    coordinate_type tbot = invR * (boxmin - ro);
    coordinate_type ttop = invR * (boxmax - ro);

    // re-order intersections to find smallest and largest on each axis
    coordinate_type tmin = min(ttop, tbot);
    coordinate_type tmax = max(ttop, tbot);

    // find the largest tmin and the smallest tmax
    T largest_tmin = max(max(tmin[0], tmin[1]), max(tmin[0], tmin[2]));
    T smallest_tmax = min(min(tmax[0], tmax[1]), min(tmax[0], tmax[2]));

    *tnear = largest_tmin;
    *tfar = smallest_tmax;
    return smallest_tmax > largest_tmin;
  }

  inline bool lineBox(coordinate_type L1, coordinate_type L2,
                      coordinate_type B1, coordinate_type B2) {
    T tnear, tfar;
    if (norm(L2 - L1) < 1e-12)
      return false;
    coordinate_type Ld = L2 - L1;
    bool hit = rayIntersectBox(L1, Ld, B1, B2, &tnear, &tfar);
    if (!hit)
      return false;
    T dist = norm(L2 - L1);
    if (tnear > 0 && tnear > dist)
      return false;
    if (tnear < 0 && tfar < 0)
      return false;
    if (tnear > 1 && tfar > 1)
      return false;
    // std::cout << tnear << " " << tfar << std::endl;
    coordinate_type p0 = L1 + Ld * tnear;
    coordinate_type p1 = L1 + Ld * tfar;

    return true;
  }
};

template <typename SPACE> class geometry_calculator {
public:
  M2_TYPEDEFS;
  coordinate_type center(coordinate_type p) { return p; }

  coordinate_type center(line_type p) { return 0.5 * (p[0] + p[1]); }

  coordinate_type center(triangle_type p) {
    return 0.33333 * (p[0] + p[1] + p[2]);
  }

  void getExtents(coordinate_type p, coordinate_type &min,
                  coordinate_type &max) {
    min = p;
    max = p;
  }

  void getExtents(line_type p, coordinate_type &min, coordinate_type &max) {
    min = this->center(p);
    max = min;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 2; j++) {
        min[i] = p[j][i] < min[i] ? p[j][i] : min[i];
        max[i] = p[j][i] > max[i] ? p[j][i] : max[i];
      }
    }
  }

  void getExtents(swept_point_type p, coordinate_type &min,
                  coordinate_type &max) {
    min = p.p;
    max = min;
    T dt = p.dt;
    for (int i = 0; i < 3; i++) {
      min[i] = p.p[i] < min[i] ? p.p[i] : min[i];
      max[i] = p.p[i] > max[i] ? p.p[i] : max[i];
      min[i] = (p.p[i] + dt * p.v[i]) < min[i] ? p.p[i] + dt * p.v[i] : min[i];
      max[i] = (p.p[i] + dt * p.v[i]) > max[i] ? p.p[i] + dt * p.v[i] : max[i];
    }
  }

  void getExtents(triangle_type p, coordinate_type &min, coordinate_type &max) {
    min = this->center(p);
    max = min;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        min[i] = p[j][i] < min[i] ? p[j][i] : min[i];
        max[i] = p[j][i] > max[i] ? p[j][i] : max[i];
      }
    }
  }
};

template <typename SPACE> class distance_calculator {
public:
  M2_TYPEDEFS;

  T calcEdgeEdgeDistance(edge_ptr e0, edge_ptr e1) {
    coordinate_type x00 = get_coordinate(e0->v1());
    coordinate_type x01 = get_coordinate(e0->v2());
    coordinate_type x10 = get_coordinate(e1->v1());
    coordinate_type x11 = get_coordinate(e1->v2());
    return va::distance_Segment_Segment(x00, x01, x10, x11);
  }

  T distance(line_type e0, line_type e1) {
    return va::distance_Segment_Segment(e0[0], e0[1], e1[0], e1[1]);
  }

  T distance(triangle_type t, coordinate_type c) {
    return va::distance_from_triangle(t.p, c);
  }

  inline void dimPolynomial(coordinate_type a, coordinate_type va,
                            coordinate_type b, coordinate_type vb,
                            coordinate_type c, coordinate_type vc, int i, int j,
                            int k, T &o, T &p, T &q, T &r) {
    o = va[i] * vb[j] * vc[k];
    p = a[i] * vb[j] * vc[k] + va[i] * b[j] * vc[k] + va[i] * vb[j] * c[k];
    q = va[i] * b[j] * c[k] + a[i] * vb[j] * c[k] + a[i] * b[j] * vc[k];
    r = a[i] * b[j] * c[k];
    // std::cout << o << " " << p << " " << q << " " << r << std::endl;
  }

  void buildPolynomial(coordinate_type a, coordinate_type va, coordinate_type b,
                       coordinate_type vb, coordinate_type c,
                       coordinate_type vc, T &p, T &q, T &r) {
    T o0, o1, o2, o3, o4, o5;
    T p0, p1, p2, p3, p4, p5;
    T q0, q1, q2, q3, q4, q5;
    T r0, r1, r2, r3, r4, r5;
    dimPolynomial(a, va, b, vb, c, vc, 1, 2, 0, o0, p0, q0, r0);
    dimPolynomial(a, va, b, vb, c, vc, 2, 1, 0, o1, p1, q1, r1);
    dimPolynomial(a, va, b, vb, c, vc, 2, 0, 1, o2, p2, q2, r2);
    dimPolynomial(a, va, b, vb, c, vc, 0, 2, 1, o3, p3, q3, r3);
    dimPolynomial(a, va, b, vb, c, vc, 0, 1, 2, o4, p4, q4, r4);
    dimPolynomial(a, va, b, vb, c, vc, 1, 0, 2, o5, p5, q5, r5);
    T o = o0 - o1 + o2 - o3 + o4 - o5;
    p = p0 - p1 + p2 - p3 + p4 - p5;
    q = q0 - q1 + q2 - q3 + q4 - q5;
    r = r0 - r1 + r2 - r3 + r4 - r5;
    // std::cout << o << " " << p << " " << q << " " << r << std::endl;
    // o = o < 1e-12 ? 1e-12 : o;
    o += 1e-12;
    p /= o;
    q /= o;
    r /= o;
    // std::cout << o << " " << p << " " << q << " " << r << std::endl;
  }
#if 1
  T distance(swept_triangle_type t, swept_point_type c) {
    // TIMER function//TIMER(__FUNCTION__);
    // int pj = tree.permutation[j];
    T dt = t.dt;
    coordinate_type ci = c.p;
    coordinate_type vi = c.v;
    coordinate_type c0 = t[0];
    coordinate_type c1 = t[1];
    coordinate_type c2 = t[2];
    coordinate_type v0 = t.v[0];
    coordinate_type v1 = t.v[1];
    coordinate_type v2 = t.v[2];
    T p, q, r, roots[3];
    roots[0] = 0.0;
    roots[1] = 0.0;
    roots[2] = 0.0;
    buildPolynomial(c1 - c0, v1 - v0, c2 - c0, v2 - v0, ci - c0, vi - v0, p, q,
                    r);
    int rootcount =
        magnet::math::cubicSolve(p, q, r, roots[0], roots[1], roots[2]);
    T dist = 999;
    for (int k = 0; k < rootcount; k++) {

      if (roots[k] > 0 && roots[k] < dt) {
        T ti = roots[k];
        if (ti > dt)
          continue;
        if (ti < 0.0)
          continue;
        m2::distance_calculator<SPACE> calc;
        triangle_type tri(c0 + ti * v0, c1 + ti * v1, c2 + ti * v2);
        T tdist = this->distance(tri, ci + ti * vi);
        dist = dist < tdist ? dist : tdist;
        // std::cout << k << " : " << roots[k] << " " << tdist << " " << dist <<
        // " "
        //	    << p << " " << q << " " << r <<  std::endl;
      }
    }
    return dist;
  }
#endif
#if 0
    T distance(swept_triangle_type t, swept_point_type c){
      //int pj = tree.permutation[j];
      T dt = t.dt;
      coordinate_type e0 = c.p;
      coordinate_type e1 = e0 + dt*c.v;
      coordinate_type c00 = t[0]; 
      coordinate_type c01 = t[1];
      coordinate_type c02 = t[2];
      coordinate_type c10 = c00 + dt*t.v[0];
      coordinate_type c11 = c11 + dt*t.v[1];
      coordinate_type c12 = c12 + dt*t.v[2];
      T dist = 999.0, tdist = 999.0;
      coordinate_type pi; 
      bool hit = false;
      m2::line_tests<SPACE> tester;
      
      hit = tester.ray_triangle_intersect(e0, e1, c00, c01, c02, pi, tdist);
      dist = (hit && tdist < dist) ? tdist : dist;
      hit = tester.ray_triangle_intersect(e0, e1, c10, c11, c12, pi, tdist);
      dist = (hit && tdist < dist) ? tdist : dist;

      hit = tester.ray_triangle_intersect(e0, e1, c00, c01, c10, pi, tdist);
      dist = (hit && tdist < dist) ? tdist : dist;
      hit = tester.ray_triangle_intersect(e0, e1, c01, c02, c11, pi, tdist);
      dist = (hit && tdist < dist) ? tdist : dist;
      hit = tester.ray_triangle_intersect(e0, e1, c02, c00, c12, pi, tdist);
      dist = (hit && tdist < dist) ? tdist : dist;

      hit = tester.ray_triangle_intersect(e0, e1, c10, c11, c01, pi, tdist);
      dist = (hit && tdist < dist) ? tdist : dist;
      hit = tester.ray_triangle_intersect(e0, e1, c11, c12, c02, pi, tdist);
      dist = (hit && tdist < dist) ? tdist : dist;
      hit = tester.ray_triangle_intersect(e0, e1, c12, c10, c00, pi, tdist);
      dist = (hit && tdist < dist) ? tdist : dist;
      return dist;
    }
#endif
};

template <typename SPACE> struct face_bin {
  M2_TYPEDEFS;

  face_bin() {}

  face_bin(face_bin &other) { *this = other; }

  face_bin(surf_ref mesh, T dx) { this->binCenters(mesh, dx); }

  ~face_bin() {}

  face_bin &operator=(face_bin &rhs) {
    if (this != &rhs) {
      _binnedFaces = rhs.binnedFaces();
      _binStart = rhs.binStart();
      _binCounter = rhs.binnedCounter();
      _faceLoc = rhs.faceLoc();
      _faces = rhs.faces();
      _xRes = rhs.xRes();
      _yRes = rhs.yRes();
      _zRes = rhs.zRes();
      _dx = rhs.dx();
      _center = rhs.center();
      _lengths = rhs.lengths();
    }
    return *this;
  }
  // void binTriangles(surf_ref mesh, T dx){}
#if 1
  void binTriangles(surf_ref mesh, T dx) {
    face_array &faces = mesh.get_faces();
    coordinate_type gmin = mesh.calc_min();
    coordinate_type gmax = mesh.calc_max();
    _dx = dx;
    T idx = 1.0 / dx;
    coordinate_type gcen = 0.5 * (gmin + gmax);
    coordinate_type glengths = gmax - gmin;

    _center = gcen;
    //_lengths = glengths;
    // T maxl = glengths[0];
    // maxl = maxl > lengths[1] ? maxl : glengths[1];
    // maxl = maxl > lengths[2] ? maxl : glengths[2];

    _xRes = ceil(glengths[0] / _dx);
    _yRes = ceil(glengths[1] / _dx);
    _zRes = ceil(glengths[2] / _dx);
    _lengths = coordinate_type(_xRes * _dx, _yRes * _dx, _zRes * _dx);

    _binStart.resize(_xRes * _yRes * _zRes + 1, 0);
    _binCounter.resize(_xRes * _yRes * _zRes, 0);

    for (int c = 0; c < faces.size(); c++) {
      face_ptr f = faces[c];
      if (f) {
        face_vertex_ptr fv1 = f->fbegin();
        face_vertex_ptr fv2 = fv1->next();
        face_vertex_ptr fv3 = fv2->next();
        coordinate_type tri[3] = {get_coordinate(fv1), get_coordinate(fv2),
                                  get_coordinate(fv3)};

        coordinate_type cen(0, 0, 0), min, max;
        for (int j = 0; j < 3; j++) {
          cen += tri[j];
        }
        cen = tri[0] + tri[1] + tri[2];
        cen /= 3.0;

        min = cen;
        max = cen;
        for (int j = 0; j < 3; j++) {
          min[0] = tri[j][0] < min[0] ? tri[j][0] : min[0];
          min[1] = tri[j][1] < min[1] ? tri[j][1] : min[1];
          min[2] = tri[j][2] < min[2] ? tri[j][2] : min[2];
          max[0] = tri[j][0] > max[0] ? tri[j][0] : max[0];
          max[1] = tri[j][1] > max[1] ? tri[j][1] : max[1];
          max[2] = tri[j][2] > max[2] ? tri[j][2] : max[2];
        }
        int mini[3];
        nearestBin(min, mini);
        coordinate_type boxSize = idx * (max - min);
        int bRes[3] = {ceil(boxSize[0]), ceil(boxSize[1]), ceil(boxSize[2])};
        bRes[0] = bRes[0] == 0 ? 1 : bRes[0];
        bRes[1] = bRes[1] == 0 ? 1 : bRes[1];
        bRes[2] = bRes[2] == 0 ? 1 : bRes[2];

        for (int i = 0; i < bRes[0]; i++) {
          for (int j = 0; j < bRes[1]; j++) {
            for (int k = 0; k < bRes[2]; k++) {
              int ix = mini[0] + i;
              int iy = mini[1] + j;
              int iz = mini[2] + k;
              if (ix < _xRes && ix > -1 && iy < _yRes && iy > -1 &&
                  iz < _zRes && iz > -1) {
                int index = ix + iy * _xRes + iz * _xRes * _yRes;
                _binCounter[index]++;
              }
            }
          }
        }
      }
    }

    _binStart[0] = 0;
    for (int i = 1; i < _binStart.size(); i++) {
      _binStart[i] = _binStart[i - 1] + _binCounter[i - 1];
      _binCounter[i - 1] = 0;
    }
    _faceLoc.resize(faces.size());
    _binnedFaces.resize(_binStart.back());
    for (int c = 0; c < faces.size(); c++) {
      face_ptr f = faces[c];
      if (f) {
        face_vertex_ptr fv1 = f->fbegin();
        face_vertex_ptr fv2 = fv1->next();
        face_vertex_ptr fv3 = fv2->next();
        coordinate_type tri[3] = {get_coordinate(fv1), get_coordinate(fv2),
                                  get_coordinate(fv3)};

        coordinate_type cen(0, 0, 0), min, max;
        for (int j = 0; j < 3; j++) {
          cen += tri[j];
        }
        cen = tri[0] + tri[1] + tri[2];
        cen /= 3.0;

        min = cen;
        max = cen;
        for (int j = 0; j < 3; j++) {
          min[0] = tri[j][0] < min[0] ? tri[j][0] : min[0];
          min[1] = tri[j][1] < min[1] ? tri[j][1] : min[1];
          min[2] = tri[j][2] < min[2] ? tri[j][2] : min[2];
          max[0] = tri[j][0] > max[0] ? tri[j][0] : max[0];
          max[1] = tri[j][1] > max[1] ? tri[j][1] : max[1];
          max[2] = tri[j][2] > max[2] ? tri[j][2] : max[2];
        }

        int mini[3];
        nearestBin(min, mini);
        coordinate_type boxSize = idx * (max - min);
        int bRes[3] = {ceil(boxSize[0]), ceil(boxSize[1]), ceil(boxSize[2])};
        bRes[0] = bRes[0] == 0 ? 1 : bRes[0];
        bRes[1] = bRes[1] == 0 ? 1 : bRes[1];
        bRes[2] = bRes[2] == 0 ? 1 : bRes[2];
        coordinate_type hbox = 0.5 * coordinate_type(dx, dx, dx);
        for (int i = 0; i < bRes[0]; i++) {
          for (int j = 0; j < bRes[1]; j++) {
            for (int k = 0; k < bRes[2]; k++) {
              int ix = mini[0] + i;
              int iy = mini[1] + j;
              int iz = mini[2] + k;
              coordinate_type ncen =
                  _center - 0.5 * _lengths +
                  coordinate_type(ix * dx, iy * dx, iz * dx) +
                  0.5 * coordinate_type(dx, dx, dx);
              if (ix < _xRes && ix > -1 && iy < _yRes && iy > -1 &&
                  iz < _zRes && iz > -1) {
                // tri_box<SPACE> tribox;
                // bool inBox = tribox.triBoxOverlap(cen,hbox,tri);
                // box_type tbox(cen,hbox);
                bool inBox = boxOverlap(box_type(cen, hbox), tri);
                if (inBox) {
                  int index = ix + iy * _xRes + iz * _xRes * _yRes;
                  int bStart = _binStart[index];
                  int bCount = _binCounter[index];

                  _binnedFaces[bStart + bCount] = c;
                  _binCounter[index]++;
                  _faceLoc[c] =
                      index; // makes life a little easier to store the id;
                }
              }
            }
          }
        }
        _faces.push_back(f);
      }
    }
  }
#endif

#if 1
  void binCenters(surf_ref mesh, T dx) {
    // future home of barycenters???
    face_array &faces = mesh.get_faces();
    _dx = dx;
    T idx = 1.0 / dx;
    box_type bb = ci::bound<SPACE>(&mesh);
    coordinate_type gcen = bb.center;
    coordinate_type glengths = 2.0 * bb.half;

    _center = gcen;
    //_lengths = glengths;
    // T maxl = glengths[0];
    // maxl = maxl > lengths[1] ? maxl : glengths[1];
    // maxl = maxl > lengths[2] ? maxl : glengths[2];

    _xRes = ceil(glengths[0] / _dx);
    _yRes = ceil(glengths[1] / _dx);
    _zRes = ceil(glengths[2] / _dx);
    _lengths = coordinate_type(_xRes * _dx, _yRes * _dx, _zRes * _dx);

    _binStart.resize(_xRes * _yRes * _zRes + 1, 0);
    _binCounter.resize(_xRes * _yRes * _zRes, 0);

    for (int c = 0; c < faces.size(); c++) {
      face_ptr f = faces[c];
      if (f) {
        coordinate_type p = ci::center<SPACE>(f);
        int b[3];
        nearestBin(p, b);
        int index = b[0] + b[1] * _xRes + b[2] * _xRes * _yRes;
        _binCounter[index]++;
      }
    }

    _binStart[0] = 0;
    for (int i = 1; i < _binStart.size(); i++) {
      _binStart[i] = _binStart[i - 1] + _binCounter[i - 1];
      _binCounter[i - 1] = 0;
    }

    std::cout << _binStart.back() << std::endl;
    _binnedFaces.resize(_binStart.back(), 0);
    for (int c = 0; c < faces.size(); c++) {
      face_ptr f = faces[c];
      if (f) {
        coordinate_type p = ci::center<SPACE>(f);
        int b[3];
        nearestBin(p, b);
        int index = b[0] + b[1] * _xRes + b[2] * _xRes * _yRes;
        int bStart = _binStart[index];
        int bCount = _binCounter[index];

        _binnedFaces[bStart + bCount] = c;
        _binCounter[index]++;
      }
    }
  }
#endif

  void clear() {
    _binnedFaces.clear();
    _binStart.clear();
    _binCounter.clear();
    _faceLoc.clear();
    _faces.clear();
  }

  int xRes() { return _xRes; }
  int yRes() { return _yRes; }
  int zRes() { return _zRes; }
  T dx() { return _dx; }
  coordinate_type center() { return _center; }
  coordinate_type lengths() { return _lengths; }

  coordinate_type center(int x, int y, int z) {
    coordinate_type cen = _center - 0.5 * _lengths +
                          coordinate_type(x * _dx, y * _dx, z * _dx) +
                          0.5 * coordinate_type(_dx, _dx, _dx);
    return cen;
  }

  void nearestBin(coordinate_type pos, int *bini) {
    coordinate_type binf =
        pos - (_center - 0.5 * _lengths + 0.5 * coordinate_type(_dx, _dx, _dx));
    bini[0] = binf[0] / _dx + 0.5;
    bini[1] = binf[1] / _dx + 0.5;
    bini[2] = binf[2] / _dx + 0.5;
  }

  vector<face_ptr> &faces() { return _faces; }

  vector<int> &binStart() { return _binStart; }
  vector<int> &faceLoc() { return _faceLoc; }
  vector<int> &binnedFaces() { return _binnedFaces; }
  vector<int> &binnedCounter() { return _binCounter; }
  vector<face_ptr> _faces;

  // binning variables
  vector<int> _binnedFaces;
  vector<int> _binStart;
  vector<int> _binCounter;
  vector<int> _faceLoc;

  // grouping variables
  int _numGroups;

  int _xRes;
  int _yRes;
  int _zRes;
  T _dx;
  coordinate_type _center;
  coordinate_type _lengths;
};

template <typename SPACE> struct edge_bin {
  M2_TYPEDEFS;

  edge_bin() {}

  edge_bin(edge_bin &other) { *this = other; }

  edge_bin(surf_ref mesh, T dx) { this->binEdges(mesh, dx); }

  ~edge_bin() {}

  edge_bin &operator=(edge_bin &rhs) {
    if (this != &rhs) {
      _binnedEdges = rhs.binnedEdges();
      _binStart = rhs.binStart();
      _binCounter = rhs.binnedCounter();
      _edgeLoc = rhs.edgeLoc();
      _edges = rhs.edges();
      _xRes = rhs.xRes();
      _yRes = rhs.yRes();
      _zRes = rhs.zRes();
      _dx = rhs.dx();
      _center = rhs.center();
      _lengths = rhs.lengths();
    }
    return *this;
  }

  inline void incrementDDA(T &tMaxX, T &tMaxY, T &tMaxZ, int &x, int &y, int &z,
                           coordinate_type &tDelta, int &stepX, int &stepY,
                           int &stepZ, int &xRes, int &yRes, int &zRes) {

    if (tMaxX < tMaxY) {
      if (tMaxX < tMaxZ) {
        x += stepX;
        tMaxX += tDelta[0];
      } else {
        z += stepZ;
        tMaxZ += tDelta[2];
      }
    } else if (tMaxY < tMaxZ) {
      y += stepY;
      tMaxY += tDelta[1];
    } else {
      z += stepZ;
      tMaxZ += tDelta[2];
    }
  }

  void DDA_3D(int ei, coordinate_type r0, coordinate_type r1,
              coordinate_type bmin, coordinate_type bmax, int phase) {

    // puts this in voxel coordinates 0:N
    coordinate_type Nf(_xRes, _yRes, _zRes, 1);
    coordinate_type dbox = bmax - bmin;
    coordinate_type r_op = ((r0 - bmin).array() * Nf.array()) / dbox.array();
    coordinate_type r_dp = ((r1 - r0).array() * Nf.array()) / dbox.array();

    coordinate_type r0p = ((r0 - bmin).array() * Nf.array()) / dbox.array();
    coordinate_type r1p = ((r1 - bmin).array() * Nf.array()) / dbox.array();

    T eps = 1e-12;
    r_dp[0] = fabs(r_dp[0]) < eps ? eps : r_dp[0];
    r_dp[1] = fabs(r_dp[1]) < eps ? eps : r_dp[1];
    r_dp[2] = fabs(r_dp[2]) < eps ? eps : r_dp[2];

    int stepX = r_dp[0] > 0 ? 1 : -1, stepY = r_dp[1] > 0 ? 1 : -1,
        stepZ = r_dp[2] > 0 ? 1 : -1;

    coordinate_type tDelta(stepX / r_dp[0], stepY / r_dp[1], stepZ / r_dp[2],
                           0);

    int x = floor(r0p[0]);
    int y = floor(r0p[1]);
    int z = floor(r0p[2]);
    int x0 = x;
    int y0 = y;
    int z0 = z;

    int xMax = abs(floor(r_dp[0] + 0.5));
    int yMax = abs(floor(r_dp[1] + 0.5));
    int zMax = abs(floor(r_dp[2] + 0.5));

    coordinate_type voxmax(x + (stepX > 0 ? 1 : 0), y + (stepY > 0 ? 1 : 0),
                           z + (stepZ > 0 ? 1 : 0), 0);

    coordinate_type vtMax = (voxmax - r_op).array() / r_dp.array();
    T tMaxX = vtMax[0], tMaxY = vtMax[1], tMaxZ = vtMax[2];

    while (abs(x0 - x) < xMax || abs(y0 - y) < yMax || abs(z0 - z) < zMax) {

      incrementDDA(tMaxX, tMaxY, tMaxZ, x, y, z, tDelta, stepX, stepY, stepZ,
                   _xRes, _yRes, _zRes);

      if (x < 0 || y < 0 || z < 0)
        continue;
      if (x > _xRes - 1 || y > _yRes - 1 || z > _zRes - 1)
        continue;
      int index = x + y * _xRes + z * _xRes * _yRes;
      if (phase == 0) {
        _binCounter[index]++;
      } else if (phase == 1) {
        int bStart = _binStart[index];
        int bCount = _binCounter[index];
        /// std::cout << index << " " << bStart << " " << bCount << std::endl;
        _binnedEdges[bStart + bCount] = ei;
        _binCounter[index]++;
        _edgeLoc[ei] = index;
      }
    }
  }
  // void binTriangles(surf_ref mesh, T dx){}
#if 1
  void binEdges(surf_ref mesh, T dx) {
    edge_array &edges = mesh.get_edges();
    coordinate_type gmin = mesh.calc_min();
    coordinate_type gmax = mesh.calc_max();
    _dx = dx;
    T idx = 1.0 / dx;
    coordinate_type gcen = 0.5 * (gmin + gmax);
    coordinate_type glengths = gmax - gmin;

    _center = gcen;
    //_lengths = glengths;
    // T maxl = glengths[0];
    // maxl = maxl > lengths[1] ? maxl : glengths[1];
    // maxl = maxl > lengths[2] ? maxl : glengths[2];

    _xRes = ceil(0.25 * glengths[0] / _dx);
    _yRes = ceil(0.25 * glengths[1] / _dx);
    _zRes = ceil(0.25 * glengths[2] / _dx);
    _lengths = coordinate_type(_xRes * _dx, _yRes * _dx, _zRes * _dx, 0);

    _binStart.resize(_xRes * _yRes * _zRes + 1, 0);
    _binCounter.resize(_xRes * _yRes * _zRes, 0);

    for (int c = 0; c < edges.size(); c++) {
      edge_ptr e = edges[c];
      if (!e)
        continue;

      face_vertex_ptr fv1 = e->v1();
      face_vertex_ptr fv2 = e->v2();
      coordinate_type p1 = get_coordinate(fv1);
      coordinate_type p2 = get_coordinate(fv2);
      DDA_3D(c, p1, p2, gmin, gmax, 0);
    }

    _binStart[0] = 0;
    for (int i = 1; i < _binStart.size(); i++) {
      _binStart[i] = _binStart[i - 1] + _binCounter[i - 1];
      _binCounter[i - 1] = 0;
    }
    _edgeLoc.resize(edges.size());
    _binnedEdges.resize(_binStart.back());

    for (int c = 0; c < edges.size(); c++) {
      edge_ptr e = edges[c];
      if (!e)
        continue;

      face_vertex_ptr fv1 = e->v1();
      face_vertex_ptr fv2 = e->v2();
      coordinate_type p1 = get_coordinate(fv1);
      coordinate_type p2 = get_coordinate(fv2);
      DDA_3D(c, p1, p2, gmin, gmax, 1);
    }
  }

#endif
  void clear() {
    _binnedEdges.clear();
    _binStart.clear();
    _binCounter.clear();
    _edgeLoc.clear();
    _edges.clear();
  }

  int xRes() { return _xRes; }
  int yRes() { return _yRes; }
  int zRes() { return _zRes; }
  T dx() { return _dx; }
  coordinate_type center() { return _center; }
  coordinate_type lengths() { return _lengths; }

  coordinate_type center(int x, int y, int z) {
    coordinate_type cen = _center - 0.5 * _lengths +
                          coordinate_type(x * _dx, y * _dx, z * _dx) +
                          0.5 * coordinate_type(_dx, _dx, _dx);
    return cen;
  }

  void nearestBin(coordinate_type pos, int *bini) {
    coordinate_type binf =
        pos - (_center - 0.5 * _lengths + 0.5 * coordinate_type(_dx, _dx, _dx));
    bini[0] = binf[0] / _dx + 0.5;
    bini[1] = binf[1] / _dx + 0.5;
    bini[2] = binf[2] / _dx + 0.5;
  }

  vector<edge_ptr> &edges() { return _edges; }

  vector<int> &binStart() { return _binStart; }
  vector<int> &edgeLoc() { return _edgeLoc; }
  vector<int> &binnedEdges() { return _binnedEdges; }
  vector<int> &binnedCounter() { return _binCounter; }
  vector<edge_ptr> _edges;

  // binning variables
  vector<int> _binnedEdges;
  vector<int> _binStart;
  vector<int> _binCounter;
  vector<int> _edgeLoc;

  // grouping variables
  int _numGroups;

  int _xRes;
  int _yRes;
  int _zRes;
  T _dx;
  coordinate_type _center;
  coordinate_type _lengths;
};

template <typename SPACE> struct vertex_bin {
  M2_TYPEDEFS;

  vertex_bin() {}

  vertex_bin(vertex_bin &other) { *this = other; }

  vertex_bin(surf_ref mesh, T dx) { this->binPoints(mesh, dx); }

  ~vertex_bin() {}

  vertex_bin &operator=(vertex_bin &rhs) {
    if (this != &rhs) {
      _binnedVerts = rhs.binnedVerts();
      _binStart = rhs.binStart();
      _binCounter = rhs.binnedCounter();

      _verts = rhs.vertices();
      _xRes = rhs.xRes();
      _yRes = rhs.yRes();
      _zRes = rhs.zRes();
      _dx = rhs.dx();
      _center = rhs.center();
      _lengths = rhs.lengths();
    }
    return *this;
  }

#if 1
  void binPoints(surf_ref mesh, T dx) {
    vertex_array &verts = mesh.get_vertices();
    coordinate_type gmin = mesh.calc_min();
    coordinate_type gmax = mesh.calc_max();
    _dx = dx;
    T idx = 1.0 / dx;
    coordinate_type gcen = 0.5 * (gmin + gmax);
    coordinate_type glengths = gmax - gmin;

    _center = gcen;
    //_lengths = glengths;
    // T maxl = glengths[0];
    // maxl = maxl > lengths[1] ? maxl : glengths[1];
    // maxl = maxl > lengths[2] ? maxl : glengths[2];

    _xRes = ceil(glengths[0] / _dx);
    _yRes = ceil(glengths[1] / _dx);
    _zRes = ceil(glengths[2] / _dx);
    _lengths = coordinate_type(_xRes * _dx, _yRes * _dx, _zRes * _dx, 0.0);

    _binStart.resize(_xRes * _yRes * _zRes + 1, 0);
    _binCounter.resize(_xRes * _yRes * _zRes, 0);

    for (int c = 0; c < verts.size(); c++) {
      vertex_ptr v = verts[c];
      if (v) {
        coordinate_type p = get_coordinate(v);
        int b[3];
        nearestBin(p, b);
        int index = b[0] + b[1] * _xRes + b[2] * _xRes * _yRes;
        _binCounter[index]++;
      }
    }

    _binStart[0] = 0;
    for (int i = 1; i < _binStart.size(); i++) {
      _binStart[i] = _binStart[i - 1] + _binCounter[i - 1];
      _binCounter[i - 1] = 0;
    }

    std::cout << _binStart.back() << std::endl;
    _binnedVerts.resize(_binStart.back(), 0);
    for (int c = 0; c < verts.size(); c++) {
      vertex_ptr v = verts[c];
      if (v) {
        coordinate_type p = get_coordinate(v);
        int b[3];
        nearestBin(p, b);
        int index = b[0] + b[1] * _xRes + b[2] * _xRes * _yRes;
        int bStart = _binStart[index];
        int bCount = _binCounter[index];

        _binnedVerts[bStart + bCount] = c;
        _binCounter[index]++;
      }
    }
  }
#endif
  void clear() {
    _binnedVerts.clear();
    _binStart.clear();
    _binCounter.clear();
    _verts.clear();
  }

  int xRes() { return _xRes; }
  int yRes() { return _yRes; }
  int zRes() { return _zRes; }
  T dx() { return _dx; }
  coordinate_type center() { return _center; }
  coordinate_type lengths() { return _lengths; }

  coordinate_type center(int x, int y, int z) {
    coordinate_type cen = _center - 0.5 * _lengths +
                          coordinate_type(x * _dx, y * _dx, z * _dx, 0.0) +
                          0.5 * coordinate_type(_dx, _dx, _dx, 0.0);
    return cen;
  }

  void nearestBin(coordinate_type pos, int *bini) {
    coordinate_type binf = pos - (_center - 0.5 * _lengths +
                                  0.5 * coordinate_type(_dx, _dx, _dx, 0.0));
    bini[0] = floor(binf[0] / _dx + 0.5);
    bini[1] = floor(binf[1] / _dx + 0.5);
    bini[2] = floor(binf[2] / _dx + 0.5);
  }

  vector<vertex_ptr> &vertices() { return _verts; }

  vector<int> &binStart() { return _binStart; }
  vector<int> &binnedVerts() { return _binnedVerts; }
  vector<int> &binnedCounter() { return _binCounter; }
  vector<vertex_ptr> _verts;

  // binning variables
  vector<int> _binnedVerts;
  vector<int> _binStart;
  vector<int> _binCounter;

  // grouping variables
  int _numGroups;

  int _xRes;
  int _yRes;
  int _zRes;
  T _dx;
  coordinate_type _center;
  coordinate_type _lengths;
};

template <typename SPACE> struct pole_node {
public:
  M2_TYPEDEFS;
  int id;
  int begin;
  int size;
  int level;
  int parent;
  coordinate_type center;
  coordinate_type centerOfMass;
  coordinate_type half;
  int children[8];
  // int neighbors[6]; to be implemented later

  pole_node() {
    id = -1;
    begin = -1;
    size = -1;
    level = -1;
    parent = -1;
    for (int i = 0; i < 8; i++)
      children[i] = -1;
  }

  ~pole_node() {}

  pole_node(const pole_node &rhs) {
    center = rhs.center;
    centerOfMass = rhs.centerOfMass;
    half = rhs.half;
    id = rhs.id;
    size = rhs.size;
    begin = rhs.begin;
    parent = rhs.parent;
    level = rhs.level;
    for (int i = 0; i < 8; i++)
      children[i] = rhs.children[i];
  }

  pole_node &operator=(const pole_node &rhs) {
    // this = new pole_node();
    if (this != &rhs) {
      // children = new int[8];
      center = rhs.center;
      centerOfMass = rhs.centerOfMass;
      half = rhs.half;
      id = rhs.id;
      size = rhs.size;
      begin = rhs.begin;
      parent = rhs.parent;
      level = rhs.level;
      for (int i = 0; i < 8; i++)
        children[i] = rhs.children[i];
    }
    return *this;
  }
};

template <typename SPACE> struct pole_tree {
  // specialized for points, this makes things really easy: basically quicksort
  M2_TYPEDEFS;

public:
  typedef pole_node<SPACE> node_type;

  pole_tree() {}
  ~pole_tree() {}

  pole_tree(const pole_tree &other) { *this = other; }

  pole_tree(vector<coordinate_type> &points) { this->build(points, 20); }

  pole_tree &operator=(const pole_tree &rhs) {
    if (this != &rhs) {
      nodes = rhs.nodes;
      permutation = rhs.permutation;
      leafNodes = rhs.leafNodes;
    }
    return *this;
  }

  void calcHalfCenter(coordinate_type &half, coordinate_type &cen,
                      coordinate_type &com, vector<coordinate_type> &points,
                      vector<int> &permutation, int beg, int N) {

    com = coordinate_type(0, 0, 0);
    for (int i = beg; i < beg + N; i++)
      com += points[permutation[i]];
    com /= (T)N;
    coordinate_type min = com;
    coordinate_type max = com;
    for (int i = beg; i < beg + N; i++) {
      coordinate_type p = points[permutation[i]];
      for (int j = 0; j < 3; j++)
        min[j] = min[j] < p[j] ? min[j] : p[j];
      for (int j = 0; j < 3; j++)
        max[j] = max[j] > p[j] ? max[j] : p[j];
    }
    half = 0.5 * (max - min);
    cen = 0.5 * (max + min);
  }

  int getBin(const coordinate_type &p, const coordinate_type &c) {
    int flag = 0;
    if (p[0] > c[0])
      flag |= 1;
    if (p[1] > c[1])
      flag |= 2;
    if (p[2] > c[2])
      flag |= 4;
    return flag;
  }

  void build(vector<coordinate_type> &points, int maxLevel) {
    // TIMER function//TIMER(__FUNCTION__);
    // inititalize permutation
    if (points.empty())
      return;

    permutation.resize(points.size());
    leafIds.resize(points.size());
    // nodes.reserve(points.size()*points.size());
    // permutation.reserve(points.size());

    for (int i = 0; i < permutation.size(); i++)
      permutation[i] = i;
    node_type root;
    root.begin = 0;
    root.level = 0;
    root.size = points.size();
    root.id = nodes.size();
    root.parent = -1;
    calcHalfCenter(root.half, root.center, root.centerOfMass, points,
                   permutation, root.begin, root.size);
    stack<int> stack;

    stack.push(nodes.size());

    int size = points.size();
    size *= log(8 * size);
    nodes.reserve(size);

    nodes.push_back(root);

    while (stack.size() > 0) {
      int pNodeId = stack.top();
      stack.pop();
      node_type pNode = nodes[pNodeId];

      int beg = pNode.begin;
      int N = pNode.size;

      int cN[8], cCounter[8], cAccum[8];
      int *lPerm = new int[N];

      for (int j = 0; j < 8; j++) {
        cN[j] = 0;
        cCounter[j] = 0;
        cAccum[0] = 0;
      }

      for (int i = beg; i < beg + N; i++) {
        coordinate_type p = points[permutation[i]];
        int bin = getBin(p, pNode.centerOfMass);
        cN[bin]++;
      }

      for (int j = 1; j < 8; j++)
        cAccum[j] = cAccum[j - 1] + cN[j - 1];

      for (int i = beg; i < beg + N; i++) {
        coordinate_type p = points[permutation[i]];
        int bin = getBin(p, pNode.centerOfMass);
        lPerm[cAccum[bin] + cCounter[bin]] = permutation[i];
        cCounter[bin]++;
      }

      int ii = 0;

      for (int i = beg; i < beg + N; i++) {
        // update the global permutation with the local permutation
        permutation[i] = lPerm[ii];
        ii++;
      }

      for (int j = 0; j < 8; j++) {
        if (cN[j] == 0)
          continue;
        int cNodeId = nodes.size();
        node_type cNode;

        cNode.level = nodes[pNodeId].level + 1;
        cNode.begin = beg + cAccum[j];
        cNode.size = cN[j];
        cNode.id = cNodeId;
        cNode.parent = pNodeId;

        nodes[pNodeId].children[j] = cNodeId;

        calcHalfCenter(cNode.half, cNode.center, cNode.centerOfMass, points,
                       permutation, cNode.begin, cNode.size);
        nodes.push_back(cNode);
        if (cNode.size > 1 && cNode.level < maxLevel) {
          stack.push(cNodeId);
        } else if (cNode.size == 1 || cNode.level == maxLevel) {
          leafNodes.push_back(cNodeId);
        } // leafIds.push_back(cNodeId);
      }
      // delete lPerm;
    }
  }

  vector<node_type> nodes;
  vector<int> leafNodes;
  vector<int> leafIds;
  vector<int> permutation;
};

template <typename SPACE> struct aabb_node {
public:
  M2_TYPEDEFS;
  int dim;
  int id;
  int begin;
  int size;
  int level;
  int parent;
  int children[2];
  box_type bbox;
  // coordinate_type centerOfMass;
  // int neighbors[6]; to be implemented later

  aabb_node() {
    dim = 0;
    id = -1;
    begin = -1;
    size = -1;
    level = -1;
    parent = -1;
    children[0] = -1;
    children[1] = -1;
  }

  ~aabb_node() {}

  aabb_node(const aabb_node &rhs) {
    bbox = rhs.bbox;
    // centerOfMass     = rhs.centerOfMass;
    dim = rhs.dim;
    id = rhs.id;
    begin = rhs.begin;
    size = rhs.size;
    size = rhs.size;
    parent = rhs.parent;
    level = rhs.level;

    children[0] = rhs.children[0];
    children[1] = rhs.children[1];
  }

  aabb_node &operator=(const aabb_node &rhs) {
    // this = new aabb_node();
    if (this != &rhs) {

      bbox = rhs.bbox;
      // centerOfMass     = rhs.centerOfMass;
      dim = rhs.dim;
      id = rhs.id;
      begin = rhs.begin;
      size = rhs.size;
      parent = rhs.parent;
      level = rhs.level;

      children[0] = rhs.children[0];
      children[1] = rhs.children[1];
    }
    return *this;
  }

  int getNumChildren() { return 2; }

  bool isLeaf() const { return children[0] < 0 && children[1] < 0; }
};

template <typename SPACE, typename PRIMITIVE> struct aabb_tree {

  M2_TYPEDEFS;

public:
  typedef aabb_node<SPACE> node_type;

  aabb_tree() {}
  ~aabb_tree() {}

  aabb_tree(const aabb_tree &other) { *this = other; }

  aabb_tree(vector<PRIMITIVE> &points) { this->build(points, 24); }

  aabb_tree &operator=(const aabb_tree &rhs) {
    if (this != &rhs) {
      nodes = rhs.nodes;
      leafNodes = rhs.leafNodes;
      permutation = rhs.permutation;
    }
    return *this;
  }

  void calcHalfCenter(coordinate_type &half, coordinate_type &cen,
                      vector<PRIMITIVE> &primitives,
                      const vector<int> &permutation, int beg, int N) {

    if (permutation.empty())
      return;

    box_type bb = primitives[permutation[beg]].bbox();
    for (int i = beg; i < beg + N; i++) {
      PRIMITIVE p = primitives[permutation[i]];
      bb.expandBy(p.bbox());
    }
    half = bb.half;
    cen = bb.center;
  }

  void build(vector<PRIMITIVE> &primitives, int maxLevel) {
    // TIMER function//TIMER(__FUNCTION__);

    // inititalize permutation
    permutation.resize(primitives.size());
    leafIds.resize(primitives.size());
    // nodes.reserve(points.size()*points.size());
    // permutation.reserve(points.size());

    for (int i = 0; i < permutation.size(); i++)
      permutation[i] = i;

    node_type root;
    root.begin = 0;
    root.level = 0;
    root.size = primitives.size();
    root.id = nodes.size();
    root.parent = -1;
    if (primitives.empty())
      return;
    calcHalfCenter(root.bbox.half, root.bbox.center, primitives, permutation,
                   root.begin, root.size);

    stack<int> stack;
    nodes.reserve(log(primitives.size()) * primitives.size());
    stack.push(nodes.size());
    nodes.push_back(root);
    while (stack.size() > 0) {
      int pNodeId = stack.top();
      stack.pop();
      node_type pNode = nodes[pNodeId];

      int beg = pNode.begin;
      int N = pNode.size;
      int dim = pNode.dim;

      coordinate_type center = pNode.bbox.center;

      int cN[2] = {0, 0}, cCounter[2] = {0, 0}, cAccum[2] = {0, 0};
      std::vector<int> lPerm(N);

      for (int i = beg; i < beg + N; i++) {
        PRIMITIVE p = primitives[permutation[i]];
        int bin = (p.center()[dim] < center[dim]) ? 0 : 1;
        cN[bin]++;
      }

      for (int j = 1; j < 2; j++)
        cAccum[j] = cAccum[j - 1] + cN[j - 1];

      for (int i = beg; i < beg + N; i++) {
        PRIMITIVE p = primitives[permutation[i]];
        int bin = (p.center()[dim] < center[dim]) ? 0 : 1;
        lPerm[cAccum[bin] + cCounter[bin]] = permutation[i];
        cCounter[bin]++;
      }

      int ii = 0;

      for (int i = beg; i < beg + N; i++) {
        // update the global permutation with the local permutation
        permutation[i] = lPerm[ii];
        ii++;
      }

      for (int j = 0; j < 2; j++) {
        if (cN[j] == 0)
          continue;

        int cNodeId = nodes.size();
        node_type cNode;

        cNode.level = nodes[pNodeId].level + 1;
        cNode.begin = beg + cAccum[j];
        cNode.size = cN[j];
        cNode.id = cNodeId;
        cNode.parent = pNodeId;
        cNode.dim = (dim + 1) % 3;

        nodes[pNodeId].children[j] = cNodeId;

        calcHalfCenter(cNode.bbox.half, cNode.bbox.center, primitives,
                       permutation, cNode.begin, cNode.size);

        nodes.push_back(cNode);

        if (cNode.size < pNode.size && cNode.level < maxLevel)
          stack.push(cNodeId);
        else if (cNode.size == pNode.size || cNode.size == 1 ||
                 cNode.level == maxLevel)
          leafNodes.push_back(cNodeId);

        // leafIds.push_back(cNodeId);
      }
      // delete lPerm;
    }
  }

  vector<node_type> nodes;
  vector<int> leafIds;
  vector<int> leafNodes;
  vector<int> permutation;
};

template <typename SPACE, typename PRIMITIVE_A, typename PRIMITIVE_B>
PRIMITIVE_A
getNearest(PRIMITIVE_B &primB, const aabb_tree<SPACE, PRIMITIVE_A> &faceTree,
           const vector<PRIMITIVE_A> &primitives,
           std::function<typename SPACE::double_type(const PRIMITIVE_A &a,
                                                     const PRIMITIVE_B &b)>
               testAB,
           typename SPACE::double_type tol) {
  M2_TYPEDEFS;
  // TIMER function//TIMER(__FUNCTION__);
  typedef aabb_tree<SPACE, PRIMITIVE_A> tree_type;
  typedef typename tree_type::node_type Node;

  PRIMITIVE_A primMin;
  T dmin = std::numeric_limits<T>::infinity();

  const Node &root = faceTree.nodes[0];
  std::stack<int> cstack;
  cstack.push(0);
  bool hit = false;
  // T tol = 0.05;
  box_type boxB = primB.bbox();
  boxB.inflate(coordinate_type(tol, tol, tol));

  while (cstack.size() > 0) {
    line_tests<SPACE> test;
    int cId = cstack.top();
    cstack.pop();
    const Node &cnode = faceTree.nodes[cId];

    if (cnode.children[0] == -1 && cnode.children[1] == -1) {

      for (int k = cnode.begin; k < cnode.begin + cnode.size; k++) {

        const PRIMITIVE_A &primA = primitives[faceTree.permutation[k]];
        const box_type &boxA = primA.bbox();

        if (!boxA.overlap(boxB)) {
          continue;
        }

        T dist = testAB(primA, primB);

        if (dist < dmin && dist < std::numeric_limits<T>::infinity()) {
          dmin = dist;
          primMin = primA;
        }
      }
    }

    for (int i = 0; i < 2; i++) {
      if (cnode.children[i] > -1) {
        box_type boxA = faceTree.nodes[cnode.children[i]].bbox;
        if (boxA.overlap(boxB)) {
          cstack.push(cnode.children[i]);
        } else {
          // std::cout << cId << ": nover " << cnode.children[i] << std::endl;
          continue;
        }
      }
    }
  }
  return primMin;
};

template <typename SPACE, typename PRIMITIVE> struct dynamic_octnode {
public:
  M2_TYPEDEFS;
  int id;
  int level;
  int parent;
  bool isLeaf;
  box_type bbox;
  int *children;
  list<int> data; // index into the global primitive array
  // int neighbors[6]; to be implemented later

  dynamic_octnode() {
    isLeaf = false;
    children = new int[8];
    id = -1;
    level = -1;
    parent = -1;
    for (int i = 0; i < 8; i++)
      children[i] = -1;
  }

  ~dynamic_octnode() { delete children; }

  dynamic_octnode(const dynamic_octnode &rhs) {
    isLeaf = rhs.isLeaf;
    bbox = rhs.bbox;
    id = rhs.id;
    parent = rhs.parent;
    level = rhs.level;
    children = new int[8];
    data = rhs.data;
    for (int i = 0; i < 8; i++)
      children[i] = rhs.children[i];
  }

  dynamic_octnode &operator=(const dynamic_octnode &rhs) {
    // this = new dynamic_octnode();
    if (this != &rhs) {
      isLeaf = rhs.isLeaf;
      bbox = rhs.bbox;
      children = new int[8];
      id = rhs.id;
      parent = rhs.parent;
      level = rhs.level;
      data = rhs.data;
      for (int i = 0; i < 8; i++)
        children[i] = rhs.children[i];
    }
    return *this;
  }
};

template <typename SPACE, typename PRIMITIVE> struct dynamic_octree {
  // specialized for faces, but...
  M2_TYPEDEFS;

public:
  ///////////////////////
  // typedefs
  ///////////////////////
  typedef dynamic_octnode<SPACE, PRIMITIVE> node_type;
  typedef list<int> data_list;
  typedef typename data_list::iterator data_iterator;

  ///////////////////////
  // constructor/destructors
  ///////////////////////
  dynamic_octree() { this->maxLeafSize = 8; }
  ~dynamic_octree() {}

  dynamic_octree(const dynamic_octree &other) { *this = other; }

  dynamic_octree(vector<coordinate_type> &points) { this->build(points, 20); }

  ///////////////////////
  // operators
  ///////////////////////
  dynamic_octree &operator=(const dynamic_octree &rhs) {
    if (this != &rhs) {
      nodes = rhs.nodes;
      leafNodes = rhs.leafNodes;
    }
    return *this;
  }

  void remove(vector<PRIMITIVE> &primitives, int p, int maxLevel) {}

  void insert(vector<PRIMITIVE> &primitives, int p, int maxLevel) {
    typedef std::pair<int, int> frame_type; // node, primitive
    stack<frame_type> stack;
    // store the node and point associated with it
    if (nodes.size() == 0) {
      node_type root;
      root.bbox = this->bbox;
      root.level = 0;
      root.parent = -1;
      root.id = nodes.size();
      root.isLeaf = true;
      nodes.push_back(root);
    }
    frame_type frame;
    frame.first = 0;
    frame.second = p;
    stack.push(frame);

    while (stack.size() > 0) {
      frame_type cFrame = stack.top();
      stack.pop();
      int pNodeId = cFrame.first;
      int pi = cFrame.second;

      PRIMITIVE &pPrim = primitives[pi];
      node_type pNode = nodes[pNodeId];
      // if the maximum number of primitives in the node are exceeded,
      // then the node must be resubdivided, flushed, and all the data
      // needs to be binned further downward.

      // if not a leaf push stack until a leaf is found
      if (nodes[pNodeId].isLeaf == false) {
        for (int i = 0; i < 8; i++) {
          int childId = nodes[pNodeId].children[i];
          bool overlap = boxOverlap(primitives[pi], nodes[childId].bbox);
          if (overlap) {
            frame_type newFrame(childId, pi);
            stack.push(newFrame);
          }
        }
        continue;
      }
      // if the current node is a leaf then push the current frame into the
      // current node.
      nodes[pNodeId].data.push_back(pi);
      // if the current node is at the maximum level, continue, never subdivide
      if (nodes[pNodeId].level == maxLevel)
        continue;

      // if the current leaf is too full and less than the maximum level.
      // Subdivide it and push all of its stored data onto the stack

      if (nodes[pNodeId].data.size() > maxLeafSize) {
        coordinate_type cen = pNode.bbox.center, half = pNode.bbox.half;

        node_type cNodes[8];
        coordinate_type offset;
        for (int i = 0; i < 8; i++) {
          T xh = 0.5 * half[0];
          T yh = 0.5 * half[1];
          T zh = 0.5 * half[2];
          if (i == 0)
            offset = coordinate_type(-xh, -yh, -zh);
          if (i == 1)
            offset = coordinate_type(xh, -yh, -zh);
          if (i == 2)
            offset = coordinate_type(-xh, yh, -zh);
          if (i == 3)
            offset = coordinate_type(xh, yh, -zh);
          if (i == 4)
            offset = coordinate_type(-xh, -yh, zh);
          if (i == 5)
            offset = coordinate_type(xh, -yh, zh);
          if (i == 6)
            offset = coordinate_type(-xh, yh, zh);
          if (i == 7)
            offset = coordinate_type(xh, yh, zh);

          box_type cbox(cen + offset, 0.5 * half);
          cNodes[i].bbox = cbox;
          cNodes[i].level = pNode.level + 1;
          cNodes[i].parent = pNode.id;
          cNodes[i].id = nodes.size();
          cNodes[i].isLeaf = true;
          nodes.push_back(cNodes[i]);
          nodes[pNodeId].children[i] = cNodes[i].id;
        }
        data_iterator itb = nodes[pNodeId].data.begin();
        data_iterator ite = nodes[pNodeId].data.end();
        int d = 0;
        while (itb != ite) {
          for (int i = 0; i < 8; i++) {
            bool overlap = boxOverlap(primitives[*itb], cNodes[i].bbox);
            if (overlap) {
              frame_type newFrame(cNodes[i].id, *itb);
              stack.push(newFrame);
            }
          }
          ++itb;
          ++d;
        }
        nodes[pNodeId].isLeaf = false;
        nodes[pNodeId].data.clear();
      }
    }
  };

  void build(vector<PRIMITIVE> &prims, int maxLevel) {
    // TIMER function//TIMER(__FUNCTION__);
    for (int i = 0; i < prims.size(); i++) {
      this->insert(prims, i, maxLevel);
    }
  }

  void draw() {
    // std::cout << " number of bins: " << nodes.size() << std::endl;
    glDisable(GL_LIGHTING);
    for (int i = 0; i < nodes.size(); i++) {
      node_type &node = nodes[i];

      if (node.data.size() > 0) {
        glColor3f(0.5, 0.5, 0.5);
        node.bbox.draw();
        node.bbox.drawCenter(3.0);
      }
    }
  }

  int maxLeafSize;
  box_type bbox;
  vector<node_type> nodes;
  vector<int> leafNodes;
  vector<int> leafIds;
};

template <typename SPACE, typename PRIMITIVE>
inline bool intersectLineTest(typename SPACE::coordinate_type e0,
                              typename SPACE::coordinate_type e1,
                              dynamic_octree<SPACE, PRIMITIVE> &faceTree,
                              vector<PRIMITIVE> &corners) {
  M2_TYPEDEFS;

  typedef dynamic_octree<SPACE, PRIMITIVE> tree_type;
  typedef typename tree_type::node_type node;
  typedef typename tree_type::data_iterator data_iterator;
  node &root = faceTree.nodes[0];
  std::stack<int> cstack;
  cstack.push(0);
  bool hit = false;

  if (corners.size() == 0)
    return false;

  while (cstack.size() > 0) {
    line_tests<SPACE> test;
    int cId = cstack.top();
    cstack.pop();
    node &cnode = faceTree.nodes[cId];

    if (cnode.isLeaf) {
      // in order to make this generic, I need to turn this piece into a
      // function or create iterators for the leaves
      data_iterator itb = faceTree.nodes[cId].data.begin();
      data_iterator ite = faceTree.nodes[cId].data.end();
      int d = 0;

      while (itb != ite) {
        int ci = *itb;
        hit = test.lineTriangle(e0, e1, corners[ci].p[0], corners[ci].p[1],
                                corners[ci].p[2]);
        itb++;
      }
    }

    for (int i = 0; i < 8; i++) {
      if (cnode.children[i] > -1) {
        coordinate_type c0 = faceTree.nodes[cnode.children[i]].bbox.center;
        coordinate_type h0 = faceTree.nodes[cnode.children[i]].bbox.half;

        if (test.lineBox(e0, e1, c0 - h0, c0 + h0)) {
          cstack.push(cnode.children[i]);
        }
      }
    }
  }
  return hit;
}

template <typename SPACE, typename PRIMITIVE_A, typename PRIMITIVE_B>
inline void getAllNearest(PRIMITIVE_B &primB,
                          aabb_tree<SPACE, PRIMITIVE_A> &faceTree,
                          vector<PRIMITIVE_A> &corners,
                          std::vector<int> &collectedPrimitives,
                          typename SPACE::double_type tol) {
  M2_TYPEDEFS;
  // TIMER function//TIMER(__FUNCTION__);
  typedef aabb_tree<SPACE, PRIMITIVE_A> tree_type;
  typedef typename tree_type::node_type node;
  node &root = faceTree.nodes[0];
  std::stack<int> cstack;
  cstack.push(0);
  bool hit = false;
  // T tol = 0.05;
  while (cstack.size() > 0) {
    line_tests<SPACE> test;
    int cId = cstack.top();
    cstack.pop();
    node &cnode = faceTree.nodes[cId];

    if (cnode.children[0] == -1 && cnode.children[1] == -1) {
      // in order to make this generic, I need to turn this piece into a
      // function or create iterators for the leaves
      int beg = faceTree.nodes[cId].begin;
      int end = beg + faceTree.nodes[cId].size;
      int d = 0;
      // m2::Debugger& debug = m2::Debugger::get_instance();
      // debug.DebugLines0.push_back(cnode.bbox.center);
      // debug.DebugLines0.push_back(primB);
      while (beg != end) {
        PRIMITIVE_A primA = corners[faceTree.permutation[beg]];
        m2::distance_calculator<SPACE> calc;
        T dist = calc.distance(primA, primB);
        // if(dist < 10) std::cout << dist << std::endl;
        if (dist < tol * tol)
          collectedPrimitives.push_back(faceTree.permutation[beg]);
        // m2::Debugger& debug = m2::Debugger::get_instance();
        // debug.DebugLines0.push_back(primB);
        // debug.DebugLines0.push_back(primA.center());
        beg++;
      }
    }

    for (int i = 0; i < 2; i++) {
      if (cnode.children[i] > -1) {
        coordinate_type c = faceTree.nodes[cnode.children[i]].bbox.center;
        coordinate_type h = faceTree.nodes[cnode.children[i]].bbox.half;
        coordinate_type minA = c - h;
        coordinate_type maxA = c + h;
        coordinate_type minB;
        coordinate_type maxB;
        geometry_calculator<SPACE> calc;
        calc.getExtents(primB, minB, maxB);
        // m2::Debugger& debug = m2::Debugger::get_instance();
        // debug.DebugBoxes.push_back(c);
        // debug.DebugBoxes.push_back(h);
        if (minA[0] > maxB[0] + tol || minB[0] - tol > maxA[0] ||
            minA[1] > maxB[1] + tol || minB[1] - tol > maxA[1] ||
            minA[2] > maxB[2] + tol || minB[2] - tol > maxA[2]) {
          continue;
        } else
          cstack.push(cnode.children[i]);
      }
    }
  }
};

template <typename SPACE>
inline void
getAllNearestTriPoint(typename SPACE::coordinate_type &p,
                      aabb_tree<SPACE, typename SPACE::triangle_type> &faceTree,
                      vector<typename SPACE::triangle_type> &corners,
                      std::vector<int> &collectedPrimitives,
                      typename SPACE::double_type tol) {
  M2_TYPEDEFS;
  // TIMER function//TIMER(__FUNCTION__);
  typedef aabb_tree<SPACE, triangle_type> tree_type;
  typedef typename tree_type::node_type node;
  node &root = faceTree.nodes[0];
  std::stack<int> cstack;
  cstack.push(0);
  bool hit = false;
  // T tol = 0.05;
  while (cstack.size() > 0) {
    line_tests<SPACE> test;
    int cId = cstack.top();
    cstack.pop();
    node &cnode = faceTree.nodes[cId];

    if (cnode.children[0] == -1 && cnode.children[1] == -1) {
      // in order to make this generic, I need to turn this piece into a
      // function or create iterators for the leaves
      int beg = faceTree.nodes[cId].begin;
      int end = beg + faceTree.nodes[cId].size;
      int d = 0;
      // m2::Debugger& debug = m2::Debugger::get_instance();
      // debug.DebugLines0.push_back(cnode.bbox.center);
      // debug.DebugLines0.push_back(p);
      while (beg != end) {
        triangle_type primA = corners[faceTree.permutation[beg]];
        m2::distance_calculator<SPACE> calc;
        T dist = calc.distance(primA, p);
        // if(dist < 10) std::cout << dist << std::endl;
        if (dist < tol * tol)
          collectedPrimitives.push_back(faceTree.permutation[beg]);
        // m2::Debugger& debug = m2::Debugger::get_instance();
        // debug.DebugLines0.push_back(p);
        // debug.DebugLines0.push_back(primA.center());
        beg++;
      }
    }

    for (int i = 0; i < 2; i++) {
      if (cnode.children[i] > -1) {
        coordinate_type c = faceTree.nodes[cnode.children[i]].bbox.center;
        coordinate_type h = faceTree.nodes[cnode.children[i]].bbox.half;
        coordinate_type minA = c - h;
        coordinate_type maxA = c + h;
        coordinate_type minB;
        coordinate_type maxB;
        geometry_calculator<SPACE> calc;
        calc.getExtents(p, minB, maxB);
        // m2::Debugger& debug = m2::Debugger::get_instance();
        // debug.DebugBoxes.push_back(c);
        // debug.DebugBoxes.push_back(h);
        /*if(minA[0] > maxB[0] + tol || minB[0] - tol > maxA[0] ||
             minA[1] > maxB[1] + tol || minB[1] - tol > maxA[1] ||
             minA[2] > maxB[2] + tol || minB[2] - tol > maxA[2]){
            continue;
          }*/
        T d = 0; // sphere box intersect the dumb unoptimized way.
        for (int k = 0; k < 3; k++) {
          if (p[k] < minA[k]) {
            T e = p[k] - minA[k];
            d += e * e;
          } else if (p[k] > maxA[k]) {
            T e = p[k] - maxA[k];
            d += e * e;
          }
        }
        if (d > tol * tol) {
          continue;
        } else
          cstack.push(cnode.children[i]);
      }
    }
  }
};

template <typename SPACE, typename PRIMITIVE_A>
inline int getNearestInRange(typename SPACE::coordinate_type &p,
                             aabb_tree<SPACE, PRIMITIVE_A> &faceTree,
                             vector<PRIMITIVE_A> &corners,
                             typename SPACE::double_type tol) {
  M2_TYPEDEFS;
  // TIMER function//TIMER(__FUNCTION__);
  typedef aabb_tree<SPACE, PRIMITIVE_A> tree_type;
  typedef typename tree_type::node_type node;
  node &root = faceTree.nodes[0];
  std::stack<int> cstack;
  cstack.push(0);
  bool hit = false;
  int closestPrim = 0;
  T min = 9999;

  while (cstack.size() > 0) {
    line_tests<SPACE> test;
    int cId = cstack.top();
    cstack.pop();
    node &cnode = faceTree.nodes[cId];
    // debug.DebugBoxes.push_back(cnode.bbox.center);
    // debug.DebugBoxes.push_back(cnode.bbox.half);

    if (cnode.children[0] == -1 && cnode.children[1] == -1) {
      // in order to make this generic, I need to turn this piece into a
      // function or create iterators for the leaves
      int beg = faceTree.nodes[cId].begin;
      int end = beg + faceTree.nodes[cId].size;
      int d = 0;
      // m2::Debugger& debug = m2::Debugger::get_instance();
      // debug.DebugLines0.push_back(cnode.bbox.center);
      // debug.DebugLines0.push_back(p);
      // std::cout << beg << " " << end << std::endl;
      while (beg != end) {
        PRIMITIVE_A primA = corners[faceTree.permutation[beg]];
        m2::distance_calculator<SPACE> calc;
        T dist = calc.distance(primA, p);
        if (dist < tol + 1e-16 && dist < min) {
          min = dist;
          closestPrim = faceTree.permutation[beg];
        };
        beg++;
      }
      // debug.DebugLines0.push_back(cnode.bbox.center);
      // debug.DebugLines0.push_back(p);
    }

    for (int i = 0; i < 2; i++) {
      if (cnode.children[i] > -1) {
        coordinate_type c = faceTree.nodes[cnode.children[i]].bbox.center;
        coordinate_type h = faceTree.nodes[cnode.children[i]].bbox.half;
        coordinate_type minA = c - h;
        coordinate_type maxA = c + h;
        T d = 0; // sphere box intersect the dumb unoptimized way.
        for (int k = 0; k < 3; k++) {
          if (p[k] < minA[k]) {
            T e = p[k] - minA[k];
            d += e * e;
          } else if (p[k] > maxA[k]) {
            T e = p[k] - maxA[k];
            d += e * e;
          }
        }
        if (d > tol * tol) {
          continue;
        } else
          cstack.push(cnode.children[i]);
      }
    }
  }
  return closestPrim;
};

template <typename SPACE, typename PRIMITIVE_A>
typename SPACE::double_type
getNearestRange(typename SPACE::coordinate_type &p,
                aabb_tree<SPACE, PRIMITIVE_A> &faceTree,
                vector<PRIMITIVE_A> &corners) {
  M2_TYPEDEFS;
  ////TIMER functionTimer(__FUNCTION__);
  typedef aabb_tree<SPACE, PRIMITIVE_A> tree_type;
  typedef typename tree_type::node_type node;
  node &root = faceTree.nodes[0];
  std::stack<int> cstack;
  cstack.push(0);
  bool hit = false;
  // T tol = 0.05;
  // first we find the closest leaf as the best first guess
  T min = 9999.0; // minimum calculated distance
  int minId = 0;  // the current closest node
  int minTri = 0;

  while (cstack.size() > 0) {
    line_tests<SPACE> test;
    int cId = cstack.top();
    cstack.pop();
    node &cnode = faceTree.nodes[cId];
    // debug.DebugBoxes.push_back(cnode.bbox.center);
    // debug.DebugBoxes.push_back(cnode.bbox.half);

    if (cnode.children[0] == -1 && cnode.children[1] == -1) {
      // in order to make this generic, I need to turn this piece into a
      // function or create iterators for the leaves
      int beg = faceTree.nodes[cId].begin;
      int end = beg + faceTree.nodes[cId].size;
      int d = 0;
      // m2::Debugger& debug = m2::Debugger::get_instance();
      // debug.DebugLines0.push_back(cnode.bbox.center);
      // debug.DebugLines0.push_back(primB);
      while (beg != end) {
        PRIMITIVE_A primA = corners[faceTree.permutation[beg]];
        m2::distance_calculator<SPACE> calc;
        T dist = calc.distance(primA, p);
        if (dist < min) {
          minId = cId;
          minTri = faceTree.permutation[beg];
          min = dist;
        }
        beg++;
      }
    } else {

      if (cnode.children[0] == -1)
        cstack.push(cnode.children[1]);
      else if (cnode.children[1] == -1)
        cstack.push(cnode.children[0]);
      else {
        int pDim = cnode.dim;
        int cDim = (pDim + 1) % 3;
        coordinate_type c0 = faceTree.nodes[cnode.children[0]].bbox.center;
        coordinate_type c0h = faceTree.nodes[cnode.children[0]].bbox.half;

        coordinate_type c1 = faceTree.nodes[cnode.children[1]].bbox.center;
        coordinate_type c1h = faceTree.nodes[cnode.children[1]].bbox.half;
        if (p[pDim] < c0[pDim] + c0h[pDim])
          cstack.push(cnode.children[0]);
        if (p[pDim] > c1[pDim] - c1h[pDim])
          cstack.push(cnode.children[1]);
      }
    }
  }
  // coordinate_type cen = faceTree.nodes[minId].bbox.center;

  /*
    PRIMITIVE_A primA = corners[minTri];
    m2::distance_calculator<SPACE> calc;
    T dist = calc.distance(primA,p);
    std::cout << dist << " " << norm(p-cen) << std::endl;
    */

  // debug.DebugLines1.push_back(p);
  // debug.DebugLines1.push_back(cen);

  return sqrt(min);
};

template <typename SPACE, typename PRIMITIVE_A>
inline int getNearest(typename SPACE::coordinate_type &p,
                      aabb_tree<SPACE, PRIMITIVE_A> &faceTree,
                      vector<PRIMITIVE_A> &corners) {
  M2_TYPEDEFS;
  T min = getNearestRange<SPACE, PRIMITIVE_A>(p, faceTree, corners);
  // std::cout << min << std::endl;
  return getNearestInRange<SPACE, PRIMITIVE_A>(p, faceTree, corners, min);
};

template <typename SPACE> struct ccd {
  M2_TYPEDEFS;

  T edgeEdge(line_type e12, line_type e34) {
    // Ec = b
    Eigen::MatrixXd E(2, 2);
    E = Eigen::MatrixXd::Zero(2, 2);
    coordinate_type x21 = e12[1] - e12[0];
    coordinate_type x43 = e34[1] - e34[0];
    coordinate_type x31 = e34[0] - e12[0];

    E(0, 0) = dot(x21, x21);
    E(0, 1) = -dot(x21, x43);
    E(1, 0) = -dot(x21, x43);
    E(1, 1) = dot(x43, x43);

    Eigen::VectorXd b(2);
    b(0) = dot(x21, x31);
    b(1) = -dot(x43, x31);
    Eigen::VectorXd c = E.ldlt().solve(b);
    return coordinate_type(b(0), b(1), 0.0);
  }

  T pointTriangle(coordinate_type p, triangle_type tri, coordinate_type w) {
    // Ec = b
    Eigen::MatrixXd E(2, 2);
    E = Eigen::MatrixXd::Zero(2, 2);
    coordinate_type x13 = tri[0] - tri[2];
    coordinate_type x23 = tri[1] - tri[2];
    coordinate_type x43 = p - tri[2];

    E(0, 0) = dot(x13, x13);
    E(0, 1) = dot(x13, x23);
    E(1, 0) = dot(x13, x23);
    E(1, 1) = dot(x23, x23);

    Eigen::VectorXd b(2);
    b(0) = dot(x13, x43);
    b(1) = -dot(x23, x43);
    Eigen::VectorXd c = E.ldlt().solve(b);
    return coordinate_type(b(0), b(1), 0.0);
  }

  T timeToCollision(coordinate_type x1, coordinate_type v1, coordinate_type x2,
                    coordinate_type v2, coordinate_type x3, coordinate_type v3,
                    coordinate_type x4, coordinate_type v4) {}
};

} // namespace m2
#endif
