//
//  modify.hpp
//  Manifold
//
//  Created by John Delaney on 5/21/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#ifndef __M2MODIFY__
#define __M2MODIFY__
#include "m2Includes.h"
#include <cmath>

namespace m2 {
template <typename SPACE> class affine : public default_interface<SPACE> {
  M2_TYPEDEFS;

public:
  affine() {}
  ~affine() {}

  void translate(surf_ptr obj, T x, T y, T z) {
    vertex_array &verts = obj->get_vertices();
    long sz = verts.size();
    for (long i = 0; i < sz; i++) {
      coordinate_type cd = this->coordinate(verts[i]);
      this->coordinate(cd + coordinate_type(x, y, z), verts[i]);
    }
    obj->update_all();
  }

  void move_to(surf_ptr obj, coordinate_type newCenter) {
    vertex_array &verts = obj->get_vertices();

    coordinate_type oldCenter = center(obj);
    coordinate_type dC = newCenter - oldCenter;
    long sz = verts.size();
    for (long i = 0; i < verts.size(); i++) {
      coordinate_type ci = this->coordinate(verts[i]);
      this->coordinate(ci + dC, verts[i]);
    }
    obj->update_all();
  }

  void scale(surf_ptr obj, T x, T y, T z) {
    vertex_array &verts = obj->get_vertices();
    long sz = verts.size();
    for (long i = 0; i < sz; i++) {
      coordinate_type ci = this->coordinate(verts[i]);
      this->coordinate(ci.array() * coordinate_type(x, y, z).array(), verts[i]);
    }
    obj->update_all();
  }

  void rotate(surf_ptr obj, mat4 R) {
    vertex_array &verts = obj->get_vertices();
    long sz = verts.size();
    for (long i = 0; i < sz; i++) {
      coordinate_type ci = this->coordinate(verts[i]);
      this->coordinate(R.transform(ci), verts[i]);
    }
    obj->update_all();
  }

  void center(surf_ref in) {

    vertex_array &tVerts = in.get_vertices();
    int fc = 0;
    coordinate_type gcen(0.0, 0.0, 0.0);
    //	for (long i = fc; i < fc+1; i++) {
    for (long i = 0; i < tVerts.size(); i++) {
      vertex_ptr v = tVerts[i];
      coordinate_type c = this->coordinate(v);
      gcen += c;
    }

    gcen *= 1.0 / (T)tVerts.size();
    coordinate_type gmin(gcen);
    coordinate_type gmax(gcen);

    for (long i = 0; i < tVerts.size(); i++) {
      
      vertex_ptr v = tVerts[i];
      coordinate_type c = this->coordinate(v);
      gmin[0] = c[0] < gmin[0] ? c[0] : gmin[0];
      gmin[1] = c[1] < gmin[1] ? c[1] : gmin[1];
      gmin[2] = c[2] < gmin[2] ? c[2] : gmin[2];
      gmax[0] = c[0] > gmax[0] ? c[0] : gmax[0];
      gmax[1] = c[1] > gmax[1] ? c[1] : gmax[1];
      gmax[2] = c[2] > gmax[2] ? c[2] : gmax[2];
    }

    coordinate_type dl = gmax - gmin;
    T maxl = dl[0];
    maxl = maxl > dl[1] ? maxl : dl[1];
    maxl = maxl > dl[2] ? maxl : dl[2];
    T s = 1.0 / maxl;
    m2::affine<space3> mod;
    std::cout << "offset: " << gcen << std::endl;
    std::cout << "scale: " << s << std::endl;
    mod.translate(&in, -gcen[0], -gcen[1], -gcen[2]);
    mod.scale(&in, s, s, s);
  }

  void centerGeometry(surf_ref in) {
    vector<vertex_ptr> &tVerts = in.get_vertices();
    int fc = 0;

    box_type bb = this->bound(&in);
    coordinate_type cen = bb.center;

    coordinate_type dl = 2.0 * bb.half;
    T maxl = dl[0];
    maxl = maxl > dl[1] ? maxl : dl[1];
    maxl = maxl > dl[2] ? maxl : dl[2];
    T s = 2.0 / maxl;
    m2::affine<SPACE> mod;
    std::cout << "offset: " << cen.transpose() << std::endl;
    std::cout << "scale: " << s << std::endl;
    this->translate(&in, -cen[0], -cen[1], -cen[2]);
    this->scale(&in, s, s, s);
  };

  void scale_face(face_ptr f1, T scale) {
    f1->update_center();
    face_vertex_ptr fvb = f1->fbegin();
    face_vertex_ptr fve = f1->fend();
    bool it = true;
    while (it) {
      it = fvb != fve;
      Eigen::Matrix<T, 3, 1> ins = this->coordinate(fvb) - this->center(f1);
      ins.normalize();
      this->coordinate(fvb) += ins * scale;
      fvb = fvb->next();
    }
    f1->update_normal();
  }

  void translate_face_along_vector(face_ptr f1, coordinate_type c1, T offset) {
    face_vertex_ptr fvb = f1->fbegin();
    face_vertex_ptr fve = f1->fend();
    bool it = true;
    while (it) {
      it = fvb != fve;
      this->coordinate(fvb) += c1 * offset;
      fvb = fvb->next();
    }
  };
}; // Class Modify
}; // namespace m2
#endif
