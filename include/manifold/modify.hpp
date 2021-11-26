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
template <typename SPACE> class modify {
  M2_TYPEDEFS;

public:
  modify() {}
  ~modify() {}

  void translate(surf_ptr obj, T x, T y, T z) {
    vertex_array &verts = obj->get_vertices();
    long sz = verts.size();
    for (long i = 0; i < sz; i++) {
      coordinate_type &cd = verts[i]->coordinate();
      cd[0] += x;
      cd[1] += y;
      cd[2] += z;
    }
    obj->update_all();
  }

  void move_to(surf_ptr obj, coordinate_type newCenter) {
    vertex_array &verts = obj->get_vertices();

    coordinate_type oldCenter = obj->calc_center();
    coordinate_type dC = newCenter - oldCenter;
    long sz = verts.size();
    for (long i = 0; i < verts.size(); i++) {
      verts[i]->coordinate()[0] += dC[0];
      verts[i]->coordinate()[1] += dC[1];
      verts[i]->coordinate()[2] += dC[2];
    }
    obj->update_all();
  }

  void scale(surf_ptr obj, T x, T y, T z) {
    vertex_array &verts = obj->get_vertices();
    long sz = verts.size();
    for (long i = 0; i < sz; i++) {
      coordinate_type &cd = verts[i]->coordinate();
      cd[0] *= x;
      cd[1] *= y;
      cd[2] *= z;
    }
    obj->update_all();
  }

  void rotate(surf_ptr obj, mat4 R) {
    vertex_array &verts = obj->get_vertices();
    long sz = verts.size();
    for (long i = 0; i < sz; i++) {
      coordinate_type &cd = verts[i]->coordinate();
      coordinate_type cd4 = R.transform(cd);
      cd = coordinate_type(cd4[0], cd4[1], cd4[2]);
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
      gcen += coordinate_type(v->coordinate()[0], v->coordinate()[1],
                              v->coordinate()[2]);
    }
    gcen *= 1.0 / (T)tVerts.size();
    coordinate_type gmin(gcen);
    coordinate_type gmax(gcen);

    for (long i = 0; i < tVerts.size(); i++) {
      vertex_ptr v = tVerts[i];
      gmin[0] = v->coordinate()[0] < gmin[0] ? v->coordinate()[0] : gmin[0];
      gmin[1] = v->coordinate()[1] < gmin[1] ? v->coordinate()[1] : gmin[1];
      gmin[2] = v->coordinate()[2] < gmin[2] ? v->coordinate()[2] : gmin[2];
      gmax[0] = v->coordinate()[0] > gmax[0] ? v->coordinate()[0] : gmax[0];
      gmax[1] = v->coordinate()[1] > gmax[1] ? v->coordinate()[1] : gmax[1];
      gmax[2] = v->coordinate()[2] > gmax[2] ? v->coordinate()[2] : gmax[2];
    }

    coordinate_type dl = gmax - gmin;
    T maxl = dl[0];
    maxl = maxl > dl[1] ? maxl : dl[1];
    maxl = maxl > dl[2] ? maxl : dl[2];
    T s = 1.0 / maxl;
    m2::modify<space3> mod;
    std::cout << "offset: " << gcen << std::endl;
    std::cout << "scale: " << s << std::endl;
    mod.translate(&in, -gcen[0], -gcen[1], -gcen[2]);
    mod.scale(&in, s, s, s);
  }

  void centerGeometry(surf_ref in) {
    M2_TYPEDEFS
    vector<vertex_ptr> &tVerts = in.get_vertices();
    int fc = 0;
    coordinate_type cen = in.calc_center();
    box_type bb = in.calc_bbox();

    coordinate_type dl = 2.0 * bb.half;
    T maxl = dl[0];
    maxl = maxl > dl[1] ? maxl : dl[1];
    maxl = maxl > dl[2] ? maxl : dl[2];
    T s = 2.0 / maxl;
    m2::modify<SPACE> mod;
    std::cout << "offset: " << cen.transpose() << std::endl;
    std::cout << "scale: " << s << std::endl;
    mod.translate(&in, -cen[0], -cen[1], -cen[2]);
    mod.scale(&in, s, s, s);
  };

  void scale_face(face_ptr f1, T scale) {
    f1->update_center();
    face_vertex_ptr fvb = f1->fbegin();
    face_vertex_ptr fve = f1->fend();
    bool it = true;
    while (it) {
      it = fvb != fve;
      Eigen::Matrix<T, 3, 1> ins = fvb->coordinate() - f1->center();
      ins.normalize();
      fvb->coordinate() += ins * scale;
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
      fvb->coordinate() += c1 * offset;
      fvb = fvb->next();
    }
    f1->update_normal();
  };

  edge_ptr flip_edge(edge_ptr e1) {
    // specialized for triangle meshes
    face_vertex_ptr fv10 = e1->v1();
    face_vertex_ptr fv11 = fv10->prev();
    face_vertex_ptr fv12 = fv11->prev();

    vertex_ptr v1 = fv10->vertex();
    vertex_ptr v2 = fv11->vertex();

    face_vertex_ptr fv20 = e1->v2();
    face_vertex_ptr fv21 = fv20->prev();
    face_vertex_ptr fv22 = fv21->prev();

    vertex_ptr v3 = fv20->vertex();
    vertex_ptr v4 = fv21->vertex();

    v1->remove_face_vertex(fv10);
    v2->remove_face_vertex(fv11);
    v3->remove_face_vertex(fv12);
    v3->remove_face_vertex(fv20);
    v4->remove_face_vertex(fv21);
    v1->remove_face_vertex(fv22);

    v2->add_face_vertex(fv10);
    v3->add_face_vertex(fv11);
    v4->add_face_vertex(fv12);

    v4->add_face_vertex(fv20);
    v1->add_face_vertex(fv21);
    v2->add_face_vertex(fv22);

    edge_ptr e11 = fv11->edge();
    edge_ptr e12 = fv12->edge();
    edge_ptr e21 = fv21->edge();
    edge_ptr e22 = fv22->edge();

    face_vertex_ptr fv11t = e11->other(fv11);
    face_vertex_ptr fv21t = e21->other(fv21);
    e11->other(fv11) = e12->other(fv12);
    e11->other(fv11)->edge() = e11;
    e21->other(fv21) = e22->other(fv22);
    e21->other(fv21)->edge() = e21;
    e12->other(fv12) = fv21t;
    fv21t->edge() = e12;
    e22->other(fv22) = fv11t;
    fv11t->edge() = e22;

    //			face_vertex_ptr itb1 = fv1;
    //			face_vertex_ptr ite1 = fv1->next();
    //			bool iterating = true;
    //			while (iterating) {
    //				iterating = itb1 != ite1;
    //
    //				itb1 = itb1->prev();
    //			}

    face_ptr f1 = fv10->face();
    face_ptr f2 = fv20->face();

    f1->update_all();
    f2->update_all();
    return e1;
  }
}; // Class Modify
}; // namespace m2
#endif
