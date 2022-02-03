/*
 *  untitled.h
 *  Phase Vocoder
 *
 *  Created by John Delaney on 1/23/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __TWOMANIFOLDFUNCTIONS__
#define __TWOMANIFOLDFUNCTIONS__
#include "m2Includes.h"

namespace m2 {

template <typename SPACE>
class geometry_helper : public default_interface<SPACE> {
public:
  M2_TYPEDEFS;

  inline coordinate_type periodicTangent(vector<coordinate_type> &path, int i) {

    int N = path.size();
    int im2 = i - 2;
    im2 = im2 < 0 ? im2 + N : im2;
    int im1 = i - 1;
    im1 = im1 < 0 ? im1 + N : im1;
    int ip1 = i + 1;
    ip1 = ip1 > N - 1 ? ip1 - N : ip1;
    int ip2 = i + 2;
    ip2 = ip2 > N - 1 ? ip2 - N : ip2;

    coordinate_type t =
        -path[ip2] + 8.0 * path[ip1] - 8.0 * path[im1] + path[im2];
    return t;
  }

  inline coordinate_type tangent(vector<coordinate_type> &path, int i) {
    coordinate_type t;
    int N = path.size();
#if 1
    if (i == 0) {
      // t = -path[i+2] + 8.0*path[i+1] - 8.0*path[i-1] + path[i-2];
      t = -path[i + 3] + 8.0 * path[i + 2] - 8.0 * path[i + 1] + path[i];
    } else if (i == N - 1) {
      t = path[N - 1] - path[N - 2];
    } else if (i < N - 3 && i > 1) {
      int N = path.size();
      t = -path[i + 2] + 8.0 * path[i + 1] - 8.0 * path[i - 1] + path[i - 2];
    } else {
      t = path[i + 1] - path[i - 1];
    }
#endif
    return t;
  }

  inline quat getRotation(coordinate_type t0, coordinate_type t1) {

    coordinate_type N = cross(t0, t1);
    if (N.mag() < 1e-10)
      return quat(0, 1, 0, 0);
    double cost = t0[0] * t1[0] + t0[1] * t1[1] + t0[2] * t1[2];
    double sint = (cross(t0, t1)).mag();
    double thet = 0.5 * acos(cost);

    N.normalize();
    N = sin(thet) * N;
    quat q;
    q.set(cos(thet), N[0], N[1], N[2]);
    return q;
  }

  inline quat getRotation(vector<coordinate_type> &path, int i) {
    coordinate_type t0 = tangent(path, i), t1 = tangent(path, i + 1);
    t0.normalize();
    t1.normalize();

    return getRotation(t0, t1);
  }

  inline quat getPeriodicRotation(vector<coordinate_type> &path, int i) {
    int N = path.size();
    coordinate_type t0 = periodicTangent(path, i),
                    t1 = periodicTangent(path, (i + 1) % N);
    t0.normalize();
    t1.normalize();

    return getRotation(t0, t1);
  }

  inline quat parallelTransport(vector<coordinate_type> &path, quat &q0,
                                int i) {

    quat q = getRotation(path, i);
    // std::cout << q[0] << " " << q[1] << " " << q[2] << " " << q[3] <<
    // std::endl;
    quat q1 = q * q0;
    q1.normalize();
    return q1;
  }

  inline quat rotorFromFrame(const coordinate_type &t, const coordinate_type &n,
                             const coordinate_type &b) {

    // T trace = m[0]+m[5]+m[10];
    // w = sqrt(1. + trace)*0.5;

    // if(trace > 0.) {
    // 	x = (m[9] - m[6])/(4.*w);
    // 	y = (m[2] - m[8])/(4.*w);
    // 	z = (m[4] - m[1])/(4.*w);
    // }
    // else {
    // 	if(m[0] > m[5] && m[0] > m[10]) {
    // 		// m[0] is greatest
    // 		x = sqrt(1. + m[0]-m[5]-m[10])*0.5;
    // 		w = (m[9] - m[6])/(4.*x);
    // 		y = (m[4] + m[1])/(4.*x);
    // 		z = (m[8] + m[2])/(4.*x);
    // 	}
    // 	else if(m[5] > m[0] && m[5] > m[10]) {
    // 		// m[1] is greatest
    // 		y = sqrt(1. + m[5]-m[0]-m[10])*0.5;
    // 		w = (m[2] - m[8])/(4.*y);
    // 		x = (m[4] + m[1])/(4.*y);
    // 		z = (m[9] + m[6])/(4.*y);
    // 	}
    // 	else { //if(m[10] > m[0] && m[10] > m[5]) {
    // 		// m[2] is greatest
    // 		z = sqrt(1. + m[10]-m[0]-m[5])*0.5;
    // 		w = (m[4] - m[1])/(4.*z);
    // 		x = (m[8] + m[2])/(4.*z);
    // 		y = (m[9] + m[6])/(4.*z);
    // 	}
    // }

    // M[0] = b[0];
    // M[1] = b[1];
    // M[2] = b[2];
    // M[3] = 0.0;

    // M[4] = n[0];
    // M[5] = n[1];
    // M[6] = n[2];
    // M[7] = 0.0;

    // M[8]  = t1[0];
    // M[9]  = t1[1];
    // M[10] = t1[2];
    // M[11] = 0.0;

    // M[12] = 0;
    // M[13] = 0;
    // M[14] = 0;
    // M[15] = 1.0;
    T x, y, z, w;
    T trace = b[0] + n[1] + t[2];
    w = sqrt(1. + trace) * 0.5;

    if (trace > 0.) {
      x = (t[1] - n[2]) / (4. * w);
      y = (b[2] - t[0]) / (4. * w);
      z = (n[0] - b[1]) / (4. * w);
    } else {
      if (b[0] > n[1] && b[0] > t[2]) {
        // m[0] is greatest
        x = sqrt(1. + b[0] - n[1] - b[2]) * 0.5;
        w = (t[1] - n[2]) / (4. * x);
        y = (n[0] + b[1]) / (4. * x);
        z = (t[0] + b[2]) / (4. * x);
      } else if (n[1] > b[0] && n[1] > t[2]) {
        // m[1] is greatest
        y = sqrt(1. + n[1] - b[0] - t[2]) * 0.5;
        w = (b[2] - t[0]) / (4. * y);
        x = (n[0] + b[1]) / (4. * y);
        z = (t[1] + n[2]) / (4. * y);
      } else { // if(m[10] > m[0] && m[10] > m[5]) {
        // m[2] is greatest
        z = sqrt(1. + t[2] - b[0] - n[1]) * 0.5;
        w = (n[0] - b[1]) / (4. * z);
        x = (t[0] + b[2]) / (4. * z);
        y = (t[1] + n[2]) / (4. * z);
      }
    }
    return quat(w, x, y, z);
  }

}; // geometry helper

template <typename SPACE> class construct : public default_interface<SPACE> {
  M2_TYPEDEFS;

public:
  construct() {}
  ~construct() {}

  edge_ptr insert_edge(surf_ptr obj_in, face_vertex_ptr v1,
                       face_vertex_ptr v2) {
    face_ptr f1 = v1->face();
    face_ptr f2 = v2->face();
    bool same_face = f1 == f2;
    if (f1 == f2) {
      edge_ptr e1 = split_face(obj_in, v1, v2);
      e1->v1()->face()->update_all();
      e1->v2()->face()->update_all();
      obj_in->push_edge(e1);
      return e1;
    } else {
      face_ptr face1 = v1->face();
      face_ptr face2 = v2->face();
      edge_ptr e1 = connect_face(obj_in, v1, v2);
      return e1;
    }
  }

  edge_ptr split_face(surf_ptr obj_in, face_vertex_ptr v1, face_vertex_ptr v2) {

    face_ptr f1 = v1->face();
    face_ptr f2 = new face_type();

    // f1->print();

    edge_ptr e1 = new edge_type();

    face_vertex_ptr v1c = v1->add_prev();
    face_vertex_ptr v2c = v2->add_prev();

    v1->prev() = v2c;
    v2c->next() = v1;
    v2->prev() = v1c;
    v1c->next() = v2;

    e1->set(v1c, v2c);
    f1->fbegin() = v1c;
    f2->fbegin() = v2c;
    f1->update_all();
    f2->update_all();
    int s1 = 0, s2 = 0;

    face_vertex_ptr fvb = v1c->vertex()->fbegin();
    face_vertex_ptr fve = v1c->vertex()->fend();

    // bool iterating = true;
    // for(fvb; iterating;  fvb = fvb->vnext()){
    // 	iterating = fvb != fve;
    // 	s1++;
    // }
    // fvb = v2c->vertex()->fbegin();
    // fve = v2c->vertex()->fend();
    // iterating = true;
    // for(fvb; iterating;  fvb = fvb->vnext()){
    // 	iterating = fvb != fve;
    // 	s2++;
    // }

    obj_in->push_face(f2);
    return e1;
  }

  edge_ptr connect_face(surf_ptr obj_in, face_vertex_ptr v1,
                        face_vertex_ptr v2) {

    face_ptr face1 = v1->face();
    face_ptr face2 = v2->face();
    edge_ptr e1 = new edge_type();
    size_t flag = 0;

    if (face1->size() == 1) {
      flag += 1;
    }

    if (face2->size() == 1) {
      flag += 2;
    }

    if (face1->size() > 1 && face2->size() > 1) {
      // neither are face_vertices
      face_vertex_ptr tmp1 = v1->add_prev();
      face_vertex_ptr tmp2 = v2->add_prev();

      tmp1->next() = v2;
      v2->prev() = tmp1;

      tmp2->next() = v1;
      v1->prev() = tmp2;

      e1->set(tmp1, tmp2);
      face1->fbegin() = v1;

      face1->update_all();
      e1->verify();
      face1->verify();
    }

    else if (face1->size() == 1 && face2->size() > 1) {
      face_vertex_ptr tmp = v2->add_prev();
      tmp->next() = v1;
      v1->prev() = tmp;
      v1->next() = v2;
      v2->prev() = v1;
      face1->fbegin() = v1;
      e1->set(v1, tmp);
    }

    else if (face1->size() > 1 && face2->size() == 1) {
      face_vertex_ptr tmp = v1->add_prev();
      tmp->next() = v2;
      v2->prev() = tmp;
      v2->next() = v1;
      v1->prev() = v2;
      face1->fbegin() = v2;
      e1->set(v2, tmp);
    }

    else {
      v1->next() = v2;
      v2->prev() = v1;
      v1->prev() = v2;
      v2->next() = v1;
      e1->set(v1, v2);
      v1->face() = face1;
      v2->face() = face1;
    }

    face1->update_all();
    obj_in->remove_face(face2->position_in_set());
    obj_in->push_edge(e1);
    face1->update_all();

    bool check = true;
    return e1;
  }

  face_ptr delete_edge(surf_ptr obj_in, edge_ptr e) {

    face_vertex_ptr v1 = e->v1();
    face_vertex_ptr v2 = e->v2();

    v1->vertex()->front() = v1->vnext();
    v2->vertex()->front() = v2->vnext();

    face_ptr f1 = v1->face();
    face_ptr f2 = v2->face();
    bool same_face = f1 == f2;
    int sz1 = f1->size();
    int sz2 = f2->size();
    face_ptr fout;
    int sw = 0;

    //std::cout << e << " " << f1 << " " << f2 << std::endl;

    if (f1 == f2) {
      face_ptr f1 = disconnect_face(obj_in, e);
      fout = f1;
    }

    else {
      face_ptr f1 = join_face(obj_in, e);
      fout = f1;
    }
    if (fout) {
      fout->update_all();
    }

    return fout;
  }

  face_ptr join_face(surf_ptr obj_in, edge_ptr e) {

    // this removes an edge
    // we gather up all the face vertices

    face_vertex_ptr v1a = e->v1();
    face_vertex_ptr v1n = v1a->next();
    face_vertex_ptr v1p = v1a->prev();
    vertex_ptr v1 = v1a->vertex();

    face_vertex_ptr v2a = e->v2();
    face_vertex_ptr v2n = v2a->next();
    face_vertex_ptr v2p = v2a->prev();
    vertex_ptr v2 = v2a->vertex();
    face_ptr f1 = v1a->face();
    face_ptr f2 = v2a->face();

#if 0
      std::cout << __FUNCTION__ << std::endl;
       f1->print();
       f2->print();
       std::cout <<" 1: "<< v1 << " " << v1a << " " << v1n << " " << v1p << std::endl;
       std::cout <<" 2: "<< v2 << " " << v2a << " " << v2n << " " << v2p << std::endl;

#endif
    if (v1p == v1a && v2p == v2a) {
      obj_in->remove_face(f1->position_in_set());
      obj_in->remove_face(f2->position_in_set());
      obj_in->remove_edge(e->position_in_set());

      v1->remove_face_vertex(v1a);
      v2->remove_face_vertex(v2a);
      delete v1a;
      delete v2a;
      return NULL;
    }

    else if (v1p == v1a) {
      v2p->next() = v2n;
      v2n->prev() = v2p;

      f1->fbegin() = v2n;
      f1->update_all();
      v2->front() = v2n;

      obj_in->remove_face(f2->position_in_set());
      obj_in->remove_edge(e->position_in_set());
      v1->remove_face_vertex(v1a);
      v2->remove_face_vertex(v2a);
      delete v1a;
      delete v2a;
    }

    else if (v2p == v2a) {
      v1p->next() = v1n;
      v1n->prev() = v1p;
      f1->fbegin() = v1n;
      f1->update_all();
      v1->front() = v1n;

      obj_in->remove_face(f2->position_in_set());
      obj_in->remove_edge(e->position_in_set());
      v1->remove_face_vertex(v1a);
      v2->remove_face_vertex(v2a);
      delete v1a;
      delete v2a;
    }

    else {
      //	 std::cout << e->position_in_set() << " " <<
      // v1->position_in_set() << " " << v2->position_in_set() << " "
      //	 << v1a << " " << v2a << " " << " " << v1p << " " << v2n <<
      // std::endl;
      v1p->next() = v2n;
      v2n->prev() = v1p;
      v2p->next() = v1n;
      v1n->prev() = v2p;
      f1->fbegin() = v1p;

      f1->update_all();
      v1->front() = v2n;
      v2->front() = v1n;

      obj_in->remove_face(f2->position_in_set());
      obj_in->remove_edge(e->position_in_set());

      v1->remove_face_vertex(v1a);
      v2->remove_face_vertex(v2a);

      delete v1a;
      delete v2a;
    }

    // else{

    // }
    return f1;
  }

  face_ptr disconnect_face(surf_ptr obj_in, edge_ptr e) {
    // std::cout << "face 1 == face 2" << std::endl;
    // we gather up all the face vertices

    face_vertex_ptr v1a = e->v1();
    face_vertex_ptr v1n = v1a->next();
    face_vertex_ptr v1p = v1a->prev();
    vertex_ptr v1 = v1a->vertex();

    face_vertex_ptr v2a = e->v2();
    face_vertex_ptr v2n = v2a->next();
    face_vertex_ptr v2p = v2a->prev();
    vertex_ptr v2 = v2a->vertex();
#if 0
      std::cout << __FUNCTION__ << std::endl;
      face_ptr f1 = v1a->face();
      f1->print();
      std::cout <<" 1: "<< v1 << " " << v1a << " " << v1n << " " << v1p << std::endl;
      std::cout <<" 2: "<< v2 << " " << v2a << " " << v2n << " " << v2p << std::endl;
#endif

    if (v1n == v2a && v2n == v1a) { // two face vertices
      face_ptr f1 = v1a->face();

      if (v1a != v2a) {
        v1->remove_face_vertex(v1a);
        v2->remove_face_vertex(v2a);

        obj_in->remove_face(f1->position_in_set());
        obj_in->remove_edge(e->position_in_set());
        delete v1a;
        delete v2a;
      } else {
        v1->remove_face_vertex(v1a);

        obj_in->remove_face(f1->position_in_set());
        obj_in->remove_edge(e->position_in_set());
        delete v1a;
      }
      return NULL;
    }

    else if (v1n == v2a) { // dangling point
      face_ptr f1 = v1a->face();

      // connect up stuff.
      v1p->next() = v2n;
      v2n->prev() = v1p;

      f1->fbegin() = v2n;
      f1->update_vertex_faces();
      f1->update_all();

      if (v1a != v2a) {

        obj_in->remove_edge(e->position_in_set());
        v1->remove_face_vertex(v1a);
        v2->remove_face_vertex(v2a);
        delete v1a;
        delete v2a;
      } else {

        obj_in->remove_edge(e->position_in_set());
        v1->remove_face_vertex(v1a);
        delete v1a;
      }

      v1->front() = v2n;

      return f1;
    }

    else if (v2n == v1a) { // dangling point
      face_ptr f1 = v2a->face();
      face_vertex_ptr fvb = f1->fbegin();
      face_vertex_ptr fve = f1->fend();
      bool iterating = true;

      v2p->next() = v1n;
      v1n->prev() = v2p;

      f1->fbegin() = v1n;
      f1->update_vertex_faces();
      f1->update_all();

      if (v1a != v2a) {

        obj_in->remove_edge(e->position_in_set());
        v1->remove_face_vertex(v1a);
        v2->remove_face_vertex(v2a);

        delete v1a;
        delete v2a;
      } else {
        v1->remove_face_vertex(v1a);
        obj_in->remove_edge(e->position_in_set());
        delete v1a;
      }

      v2->front() = v1n;
      return f1;
    } else { // full pipe

      face_ptr f1 = v1a->face();

      if (v1a != v2a) {
        face_ptr nf = new face_type();
        f1->fbegin() = v1n;
        nf->fbegin() = v2n;
        v2p->next() = v1n;
        v1n->prev() = v2p;
        v1p->next() = v2n;
        v2n->prev() = v1p;

        f1->update_vertex_faces();
        nf->update_vertex_faces();
        f1->update_all();
        nf->update_all();

        obj_in->remove_edge(e->position_in_set());
        v1->remove_face_vertex(v1a);
        v2->remove_face_vertex(v2a);
        delete v1a;
        delete v2a;
        v1->front() = v2n;
        v2->front() = v1n;
        obj_in->push_face(nf);
        return nf;
      } else {
        f1->fbegin() = v1n;
        v2p->next() = v1n;
        v1n->prev() = v2p;
        v1p->next() = v2n;
        v2n->prev() = v1p;

        f1->update_vertex_faces();
        f1->update_all();
        v1->remove_face_vertex(v1a);
        v1->front() = v1n;
        obj_in->remove_edge(e->position_in_set());
        delete v1a;
        return f1;
      }
    }
  }

  vertex_ptr subdivide_edge(surf_ptr obj_in, edge_ptr edge_in) {

    face_vertex_ptr fv1 = edge_in->v1();
    face_vertex_ptr fv2 = edge_in->v2();

    coordinate_type c1 = this->coordinate(fv1);
    coordinate_type c2 = this->coordinate(fv2);
    coordinate_type cn = 0.5 * (c1 + c2);

    vertex_ptr v1 = fv1->vertex();
    vertex_ptr v2 = fv2->vertex();

    vertex_ptr vn = new vertex_type();
    this->coordinate(cn, vn);
    edge_ptr e1 = edge_in;
    edge_ptr e2 = new edge_type();
    face_vertex_ptr fv1n = fv1->add_next();
    face_vertex_ptr fv2n = fv2->add_next();

    fv1->vertex()->front() = fv1;
    fv2->vertex()->front() = fv2;

    v1->remove_face_vertex(fv1n);
    v2->remove_face_vertex(fv2n);

    fv1n->face() = fv1->face();
    fv2n->face() = fv2->face();

    vn->add_face_vertex(fv1n);
    vn->add_face_vertex(fv2n);
    vn->front() = fv1n;

    e1->set(fv1, fv2n);
    e2->set(fv1n, fv2);

    e1->flag = 1;
    e2->flag = 1;
    vn->flag = 1;

    obj_in->push_vertex(vn);
    obj_in->push_edge(e2);

    fv1n->flag = 1;
    fv2n->flag = 1;
    if (v1->pinned == true && v2->pinned == true) {
      vn->pinned = true;
    } else
      vn->pinned = false;

    return vn;
  }

  list<edge_ptr> pipe_face(surf_ptr obj_in, face_ptr f1, face_ptr f2) {
    list<edge_ptr> new_edges;
    face_vertex_ptr fv1 = f1->fbegin(), fv1e = f1->fend(), fv2 = f2->fbegin(),
                    fv2e = f2->fend(), fv1t = fv1, fv2t = fv2;

    bool it1 = true, it2 = true;
    coordinate_type c0 = this->coordinate(fv1);
    T mz1 = c0[2];
    T mz2 = c0[2];

    it1 = true;
    it2 = true;
    coordinate_type c1, c2;

    while (it1) {
      it1 = fv1 != fv1e;
      c1 = this->coordinate(fv1);

      if (c1[2] > mz1) {
        mz1 = c1[2];
        fv1t = fv1;
      }
      fv1 = fv1->next();
    }

    while (it2) {
      it2 = fv2 != fv2e;
      c2 = this->coordinate(fv2);

      if (c2[2] > mz2) {
        mz2 = c2[2];
        fv2t = fv2;
      }
      fv2 = fv2->next();
    }

    vector<face_vertex_ptr> vfv1;
    vector<face_vertex_ptr> vfv2;

    fv1 = fv1t;
    fv2 = fv2t;
    fv1e = fv1->prev();
    int i = 0;
    it1 = true;
    while (it1) {
      it1 = fv1 != fv1e;
      i++;
      vfv1.push_back(fv1);
      vfv2.push_back(fv2);

      fv1 = fv1->next();
      fv2 = fv2->prev();
    }

    typename vector<face_vertex_ptr>::iterator itb1 = vfv1.begin();
    typename vector<face_vertex_ptr>::iterator itb2 = vfv2.begin();
    while (itb1 != vfv1.end()) {
      edge_ptr ne = insert_edge(obj_in, *itb1, *itb2);
      itb2++;
      itb1++;
      new_edges.push_back(ne);
    }
    return new_edges;
  }

  list<edge_ptr> stitch_unlike_faces(surf_ptr obj_in, face_ptr f1,
                                     face_ptr f2) {
    list<edge_ptr> new_edges;
    face_vertex_ptr fv1 = f1->fbegin(), fv1e = f1->fend(),

                    fv2 = f2->fbegin(), fv2e = f2->fend(),

                    fv1t = fv1, fv2t = fv2;

    bool it1 = true, it2 = true;

    coordinate_type c0 = this->coordinate(fv1);
    T mz1 = c0[2];
    T mz2 = c0[2];

    it1 = true;
    it2 = true;
    coordinate_type c1, c2;

    while (it1) {
      it1 = fv1 != fv1e;
      c1 = this->coordinate(fv1);

      if (c1[2] > mz1) {
        mz1 = c1[2];
        fv1t = fv1;
      }
      fv1 = fv1->next();
    }

    while (it2) {
      it2 = fv2 != fv2e;
      c2 = this->coordinate(fv2);

      if (c2[2] > mz2) {
        mz2 = c2[2];
        fv2t = fv2;
      }
      fv2 = fv2->next();
    }

    vector<face_vertex_ptr> vfv1;
    vector<face_vertex_ptr> vfv2;

    fv1 = fv1t;
    fv2 = fv2t;
    fv1e = fv1->prev();
    int i = 0;
    it1 = true;
    while (it1) {
      it1 = fv1 != fv1e;
      i++;
      vfv1.push_back(fv1);
      vfv2.push_back(fv2);

      fv1 = fv1->next();
      fv2 = fv2->prev();
    }

    typename vector<face_vertex_ptr>::iterator itb1 = vfv1.begin();
    typename vector<face_vertex_ptr>::iterator itb2 = vfv2.begin();
    while (itb1 != vfv1.end()) {
      edge_ptr ne = insert_edge(obj_in, *itb1, *itb2);
      itb2++;
      itb1++;
      new_edges.push_back(ne);
    }
    return new_edges;
  }

  void cleanup_zero_edges(surf_ptr &obj_in, vertex_ptr &v) {
    bool cleanup = true;
    while (cleanup) {
      face_vertex_ptr fvb = v->fbegin(), fve = v->fend();

      bool iterating = true;
      bool hasZeroEdge = false;
      edge_ptr zeroEdge;
      while (iterating) {
        iterating = fvb != fve;
        if (fvb->edge()->v1()->vertex() == fvb->edge()->v2()->vertex()) {
          hasZeroEdge = true;
          zeroEdge = fvb->edge();
        }
        fvb = fvb->vnext();
      }
      if (hasZeroEdge) {
        this->delete_edge(obj_in, zeroEdge);
      }
      cleanup = hasZeroEdge;
    }
  }

  face_ptr delete_vertex_primitive(surf_ptr obj_in, vertex_ptr v) {
    // we gather up all the face vertices
    // vertex_ptr v = obj_in->get_vertices()[i];
    if (v->size() < 1) {
      obj_in->remove_vertex(v->position_in_set());
      return NULL;
    }

    else {

      // cleanup_zero_edges(obj_in,v);
      face_vertex_ptr fvb = v->fbegin(), fve = v->fend();
      bool iterating = true;
      // clear any existing edge flags should they be there
      while (iterating) {
        iterating = fvb != fve;
        fvb->edge()->flag = 0;
        fvb = fvb->vnext();
      }
      iterating = true;
      fvb = v->fbegin();
      fve = v->fend();
      // iterate and grab edge pointers flagging them so as to not delete them
      // twice
      vector<edge_ptr> edgesToDelete;
      while (iterating) {
        iterating = fvb != fve;
        if (fvb->edge()->delete_flag == 0)
          edgesToDelete.push_back(fvb->edge());
        fvb->edge()->delete_flag += 1;
        fvb = fvb->vnext();
      }

      face_ptr nf;

      for (int i = 0; i < edgesToDelete.size(); i++) {
        edgesToDelete[i]->verify();
        nf = this->delete_edge(obj_in, edgesToDelete[i]);
      }

      obj_in->remove_vertex(v->position_in_set());
      return nf;
    }
  }

  face_ptr delete_vertex(surf_ptr obj_in, vertex_ptr v) {

    face_vertex_ptr fvb = v->fbegin(), fve = v->fend();
    bool iterating = true;

    face_vertex_ptr this_wing,
        inco_edge; // we'll need this after connecting everything up
    // 1. we need to loop around and gather the current vertex's next and the
    // next vertex's previous;
    face_ptr nf = new face_type();
    int j = 0;
    iterating = true;
    vector<face_vertex_ptr> chk;
    while (iterating) {
      // the corners we're grabbing should make som illogical sense...
      iterating = fvb != fve;

      this_wing = fvb->next();
      inco_edge = fvb->vnext()->prev()->prev();

      this_wing->prev() = inco_edge;
      inco_edge->next() = this_wing;

      this_wing->face() = nf;
      inco_edge->face() = nf;
      inco_edge->position_in_face() = j;
      j++;
      this_wing->position_in_face() = j;
      j++;
      this_wing->vertex()->front() = this_wing;
      chk.push_back(this_wing);
      chk.push_back(inco_edge);
      fvb = fvb->vnext();
    }
    iterating = true;
    nf->fbegin() = v->fbegin()->next();

    nf->size() = j;
    nf->update_all();
    obj_in->push_face(nf);

    // 3. now the nasty part about deletion;
    fvb = nf->fbegin();
    fve = nf->fend();

    iterating = true;
    while (iterating) {
      iterating = fvb != fve;
      fvb = fvb->next();
    }
    fvb = v->fbegin();
    fve = v->fend();
    iterating = true;
    while (iterating) {
      iterating = fvb != fve;

      face_ptr cf = fvb->face();
      edge_ptr ce = fvb->edge();
      face_vertex_ptr fvp1 = ce->v1();
      face_vertex_ptr fvp2 = ce->v2();
      // uhh... make sure we advance before we delete...
      fvb = fvb->vnext();
      obj_in->remove_face(cf->position_in_set());
      obj_in->remove_edge(ce->position_in_set());
      fvp1->vertex()->remove_face_vertex(fvp1);
      fvp2->vertex()->remove_face_vertex(fvp2);
    }
    obj_in->remove_vertex(v->position_in_set());
    bool here = true;
    return nf;
  }

  bool delete_degenerates(surf_ptr obj, vertex_ptr v, real tol = 1e-8) {
    vector<vertex_ptr> vertDegenerates;
    if(v->calc_size() < 3){
      //std::cout << " deleting: " << v->position_in_set() << std::endl;
      this->delete_vertex_primitive(obj, v);
      return true;
    }
    if (m2::ci::area<SPACE>(v) < tol) {
      this->delete_vertex_primitive(obj, v);
      return true;
    }
    return false;
  };

  bool delete_degenerates(surf_ptr obj, edge_ptr e) {
    if (m2::ci::length<SPACE>(e) == 0) {
      // std::cout << " deleting 1: " << e->position_in_set() << std::endl;
      this->delete_edge(obj, e);
      return true;
    }

    if(e->v1()->face()->size() < 3){
      //std::cout << " deleting 1: " << e->position_in_set() << std::endl;
      this->delete_edge(obj, e);
      return true;
    }
    if (e->v2()->face()->size() < 3) {
      //std::cout << " deleting 2: " << e->position_in_set() << std::endl;
      this->delete_edge(obj, e);
      return true;
    }
    return false;

  };

  vertex_ptr collapse_edge_triangle(surf_ptr obj, edge_ptr e) {
    // TODO:: move to triangle_ops, compare with joinFaceOneSharedEdge
    if (e->v1()->vertex() == e->v2()->vertex())
      return e->v1()->vertex();

    face_vertex_ptr va0 = e->v1();
    face_vertex_ptr vb0 = e->v2();
    vertex_ptr v1 = va0->vertex();
    vertex_ptr v2 = vb0->vertex();

    coordinate_type c1 = this->coordinate(v1);
    coordinate_type c2 = this->coordinate(v2);
    coordinate_type cp = 0.5 * (c1 + c2);
    this->coordinate(cp, v1);

    if (va0->face()->size() == 1 || va0->face()->size() == 2) {
      std::cout << "e pos (v2): " << e->position_in_set() << std::endl;
      v2->print();
      delete_edge(obj, e);
      return v1;
    }

    if (vb0->face()->size() == 1 || vb0->face()->size() == 2) {
      std::cout << "e pos (v2): " << e->position_in_set() << std::endl;
      v2->print();
      delete_edge(obj, e);
      return v2;
    }

    if (v1->calc_size() == 1 || v1->calc_size() == 2) {
      std::cout << "e pos (v1): " << e->position_in_set() << std::endl;
      v1->print();
      // std::cout << " deleting v1" << std::endl;
      face_ptr f = delete_vertex_primitive(obj, v1);
      return v2;
    }

    if (v2->calc_size() == 1 || v2->calc_size() == 2) {
      std::cout << "e pos (v2): " << e->position_in_set() << std::endl;
      v2->print();
      face_ptr f = delete_vertex_primitive(obj, v2);
      //	std::cout << "output 3: " << std::endl;
      // v1->print();
      return v1;
    }

    if (v1->calc_size() == 1 || v1->calc_size() == 2) {
      std::cout << "e pos (v1): " << e->position_in_set() << std::endl;
      v1->print();
      // std::cout << " deleting v1" << std::endl;
      face_ptr f = delete_vertex_primitive(obj, v1);
      return v2;
    }

    face_vertex_ptr va1 = va0->next();
    face_vertex_ptr va2 = va0->prev();
    edge_ptr ea1 = va1->edge();
    edge_ptr ea2 = va2->edge();
    face_vertex_ptr va2p = ea1->other(va1);
    face_vertex_ptr va0p = ea2->other(va2);

    face_vertex_ptr vb1 = vb0->next();
    face_vertex_ptr vb2 = vb0->prev();
    edge_ptr eb1 = vb1->edge();
    edge_ptr eb2 = vb2->edge();
    face_vertex_ptr vb2p = eb1->other(vb1);
    face_vertex_ptr vb0p = eb2->other(vb2);

    if (e->position_in_set() == 57456) {

      std::cout << " e: " << e->position_in_set() << std::endl;
      std::cout << " ea2: " << ea2->position_in_set() << std::endl;
      std::cout << " eb2: " << eb2->position_in_set() << std::endl;

      std::cout << " e: " << e->v1() << " " << e->v2() << std::endl;
      std::cout << va0 << " " << va1 << " " << va2 << " " << vb0 << " " << vb1
                << " " << vb2 << std::endl;

      std::cout << va0p << " " << va2p << " " << vb0p << " " << vb2p
                << std::endl;

      bool iterating = true;
      face_vertex_ptr fvb = v1->fbegin();
      face_vertex_ptr fve = v1->fend();
      while (iterating) {
        iterating = fvb != fve;
        std::cout << fvb->edge()->position_in_set() << " ";
        fvb = fvb->vnext();
      }
      std::cout << std::endl;

      iterating = true;
      fvb = v1->fbegin();
      fve = v1->fend();
      while (iterating) {
        iterating = fvb != fve;
        std::cout << fvb->edge()->v1() << " ";
        fvb = fvb->vnext();
      }
      std::cout << std::endl;

      iterating = true;
      fvb = v1->fbegin();
      fve = v1->fend();
      while (iterating) {
        iterating = fvb != fve;
        std::cout << fvb->edge()->v2() << " ";
        //<< fvb->next()->vertex()->calc_size() <<  " ";
        fvb = fvb->vnext();
      }
      std::cout << std::endl;
    }

    ea1->set(va0p, va2p);
    eb1->set(vb0p, vb2p);

    va0p->vertex() = v1;

    va0p->vertex()->front() = va0p;
    va2p->vertex()->front() = va2p;
    vb0p->vertex()->front() = vb0p;
    vb2p->vertex()->front() = vb2p;

    face_ptr fa = va0->face();
    face_ptr fb = vb0->face();

    v1->color.r = 1.0;

    obj->remove_face(fa->position_in_set());
    obj->remove_face(fb->position_in_set());

    obj->remove_edge(e->position_in_set());
    obj->remove_edge(ea2->position_in_set());
    obj->remove_edge(eb2->position_in_set());
    obj->remove_vertex(v2->position_in_set());
    // these gotta go last
    if (va0)
      delete va0;
    if (va1)
      delete va1;
    if (va2)
      delete va2;
    if (vb0)
      delete vb0;
    if (vb1)
      delete vb1;
    if (vb2)
      delete vb2;

    for_each(v1, [&v1](face_vertex_ptr fv) -> void {
      fv->vertex() = v1;
      fv->next()->vertex()->front() = fv->next();
    });

    // while(this->delete_degenerates(obj,v1));
    return v1;
  }

  vertex_ptr collapse_edge(surf_ptr obj, edge_ptr e) {
    if (e->v1()->size() == 0)
      return e->v2()->vertex();
    if (e->v2()->size() == 0)
      return e->v1()->vertex();

    if (e->v1()->vertex() == e->v2()->vertex())
      return e->v1()->vertex();

    face_vertex_ptr fv1 = e->v1();
    face_vertex_ptr fv2 = e->v2();
    vertex_ptr v1 = fv1->vertex();
    vertex_ptr v2 = fv2->vertex();

    // std::cout << "input: " << std::endl;
    // v1->print();
    // v2->print();

    coordinate_type c1 = this->coordinate(v1);
    coordinate_type c2 = this->coordinate(v2);
    coordinate_type cp = 0.5 * (c1 + c2);
    this->coordinate(cp, v1);
    
    /*
    if (v2->size() == 3) {

      face_ptr f = delete_vertex_primitive(obj, v2);
      //	std::cout << "output 3: " << std::endl;
      // v1->print();
      return v1;
    }
    */
    std::cout << v1->position_in_set() << ":" << v1->size() << " "
              << v2->position_in_set() << ":" << v2->size() << " " << std::endl;

    /*
    if (v2->size() == 2) {
      edge_ptr e1n = fv1->next()->edge();
      edge_ptr e2n = fv2->next()->edge();
      this->delete_edge(obj, e);
      this->delete_edge(obj, e1n);
      this->delete_edge(obj, e2n);
      // std::cout << "output 2: " << std::endl;
      // v1->print();
      // v2->print();
      obj->remove_vertex(v2->position_in_set());
      return v1;
    }
    */
    if (fv1->face()->size() == 3)
      this->delete_edge(obj, fv1->next()->edge());
    if (fv2->face()->size() == 3)
      this->delete_edge(obj, fv2->prev()->edge());

    face_vertex_ptr v1p = fv1->prev();
    face_vertex_ptr v1n = fv1->next();
    face_vertex_ptr v2p = fv2->prev();
    face_vertex_ptr v2n = fv2->next();

    fv1->face()->fbegin() = v1p;
    fv2->face()->fbegin() = v2p;
    v1p->next() = v1n;
    v1n->prev() = v1p;
    v2p->next() = v2n;
    v2n->prev() = v2p;
    fv1->face()->update_all();
    fv2->face()->update_all();

    bool iterating = true;
    face_vertex_ptr fvb = fv2->vnext();
    face_vertex_ptr fve = fv2->vprev();


    while (iterating) {
      iterating = fvb != fve;
      v2->remove_face_vertex(fvb);
      v1->add_face_vertex(fvb);
      fvb = fvb->vnext();
    }

    v1->front() = fvb;

    v1->remove_face_vertex(fv1);
    delete fv1;
    v2->remove_face_vertex(fv2);
    delete fv2;
    // std::cout << "output: " << std::endl;
    // v1->print();
    // v2->print();

    obj->remove_edge(e->position_in_set());
    obj->remove_vertex(v2->position_in_set());

    return v1;
  }

  vertex_ptr collapse_edge_primitive(surf_ptr obj, edge_ptr e) {

    face_vertex_ptr fv1 = e->v1();
    face_vertex_ptr fv2 = e->v2();

    vertex_ptr v1 = fv1->vertex();
    vertex_ptr v2 = fv2->vertex();
    coordinate_type c1 = coordinate(v1);
    coordinate_type c2 = coordinate(v2);
    coordinate_type cp = 0.5 * (c1 + c2);
    coordinate(v1, cp);

    face_vertex_ptr fv1p = fv1->vnext();
    // if(fv1->face()->size() < 3 || fv2->face()->size() < 3) {
    // 	//handle the troubling edge and abort the collapse
    // 	this->delete_edge(obj,e); return;
    // };

    face_ptr f = delete_vertex_primitive(obj, v2);
    if (!f)
      return v1;
    if (f->size() < 4)
      return v1;

    // face_vertex_ptr fvb = fv1p;
    face_vertex_ptr fvb = f->get_corner_on_vertex(v1);

    fvb->face()->print();

    bool iterating = true;
    int i = 0;
    while (iterating) {
      face_vertex_ptr fvnn = fvb->next()->next();

      if (fvb->vertex() != fvnn->vertex()) {
        edge_ptr ei = insert_edge(obj, fvb, fvnn);
        fvb = ei->v1();
      } else {
        return v1;
      }

      face_vertex_ptr fvnnn = fvb->next()->next()->next();

      iterating = fvb != fvnnn;
      i++;
    }
    return v1;
  }

  void flip_edge(surf_ptr &obj, edge_ptr &e) {

#if 1
    // valid if both faces are of size three.

    face_vertex_ptr fv1 = e->v1();
    face_vertex_ptr fv2 = e->v2();
    face_vertex_ptr fv1n = fv1->next();
    face_vertex_ptr fv2n = fv2->next();
    face_vertex_ptr fv1p = fv1->prev();
    face_vertex_ptr fv2p = fv2->prev();
    face_vertex_ptr fv1pp = fv1p->prev();
    face_vertex_ptr fv2pp = fv2p->prev();

    e->v1() = fv1p;
    e->v2() = fv2p;

    edge_ptr e1p = fv1p->edge();
    edge_ptr e2p = fv2p->edge();

    fv1p->edge() = e;
    fv2p->edge() = e;

    if (e1p->v1() == fv1p)
      e1p->v1() = fv2;
    else
      e1p->v2() = fv2;

    if (e2p->v1() == fv2p)
      e2p->v1() = fv1;
    else
      e2p->v2() = fv1;

    fv1->edge() = e2p;
    fv2->edge() = e1p;

    face_ptr f1 = fv1->face();
    face_ptr f2 = fv2->face();

    vertex_ptr v1 = fv1->vertex();
    vertex_ptr v2 = fv2->vertex();
    vertex_ptr v1p = fv1p->vertex();
    vertex_ptr v2p = fv2p->vertex();

    v2->remove_face_vertex(fv2);
    v1->remove_face_vertex(fv1);
    v1p->add_face_vertex(fv2);
    v2p->add_face_vertex(fv1);

    v1p->front() = fv2;
    v2p->front() = fv1;

    v1->front() = fv2n;
    v2->front() = fv1n;

    f1->update_all();
    f2->update_all();

#endif
  }

  bool stitch_faces(surf_ptr obj_in, face_ptr f1, face_ptr f2, T tol) {

    face_vertex_ptr fvb2 = f2->fbegin();
    face_vertex_ptr fve2 = f2->fend();
    coordinate_type fc = coordinate(f1->fbegin());

    T normv = 1000., normt = 1000.;
    bool it1 = true;

    while (it1) {
      it1 = fvb2 != fve2;
      normv = dot(fc, coordinate(fvb2));
      if (normv < normt) {
        normt = normv;
        f2->fbegin() = fvb2;
      }
      fvb2 = fvb2->next();
    }

    it1 = true;
    face_vertex_ptr fvb1 = f1->fbegin();
    fvb2 = f2->fbegin();
    face_vertex_ptr fve1 = f1->fend();

    while (it1) {
      it1 = fvb1 != fve1;
      face_vertex_ptr fv1p = fvb1->vprev();
      face_vertex_ptr fv2n = fvb2->next()->vprev();
      edge_ptr e1 = fv1p->edge();
      edge_ptr e2 = fv2n->edge();

      e1->set_other(fv1p, fv2n);
      fv2n->edge() = e1;

      obj_in->remove_edge(e2->position_in_set());
      fvb2->vertex()->add_face_vertex(fv1p);

      // fv1p->vertex() = fvb2->vertex();
      fv2n->vertex()->add_face_vertex(fv1p->next());

      fvb1 = fvb1->next();
      fvb2 = fvb2->prev();
    }

    fvb1 = f1->fbegin();
    fvb2 = f2->fbegin();

    it1 = true;
    while (it1) {
      it1 = fvb1 != fve1;
      obj_in->remove_vertex(fvb1->vertex()->position_in_set());
      fvb2->vertex()->remove_face_vertex(fvb2->position_in_vertex());

      fvb1 = fvb1->next();
      fvb2 = fvb2->next();
    }

    obj_in->remove_face(f1->position_in_set());
    obj_in->remove_face(f2->position_in_set());
    return true;
  }

  bool stitch_faces_diff_size(surf_ptr obj_in, face_ptr f1, face_ptr f2,
                              T tol) {

    face_vertex_ptr fvb2 = f2->fbegin();
    face_vertex_ptr fve2 = f2->fend();
    coordinate_type fc = coordinate(f1->fbegin());

    T normv = 1000., normt = 1000.;
    bool it1 = true;

    while (it1) {
      it1 = fvb2 != fve2;
      normv = dot(fc, coordinate(fvb2));
      if (normv < normt) {
        normt = normv;
        f2->fbegin() = fvb2;
      }
      fvb2 = fvb2->next();
    }

    it1 = true;
    face_vertex_ptr fvb1 = f1->fbegin();
    fvb2 = f2->fbegin();
    face_vertex_ptr fve1 = f1->fend();

    while (it1) {
      it1 = fvb1 != fve1;
      face_vertex_ptr fv1p = fvb1->vprev();
      face_vertex_ptr fv2n = fvb2->next()->vprev();
      edge_ptr e1 = fv1p->edge();
      edge_ptr e2 = fv2n->edge();

      e1->set_other(fv1p, fv2n);
      fv2n->edge() = e1;

      obj_in->remove_edge(e2->position_in_set());
      fvb2->vertex()->add_face_vertex(fv1p);

      // fv1p->vertex() = fvb2->vertex();
      fv2n->vertex()->add_face_vertex(fv1p->next());

      fvb1 = fvb1->next();
      fvb2 = fvb2->prev();
    }

    fvb1 = f1->fbegin();
    fvb2 = f2->fbegin();

    it1 = true;
    while (it1) {
      it1 = fvb1 != fve1;
      obj_in->remove_vertex(fvb1->vertex()->position_in_set());
      fvb2->vertex()->remove_face_vertex(fvb2->position_in_vertex());

      fvb1 = fvb1->next();
      fvb2 = fvb2->next();
    }

    obj_in->remove_face(f1->position_in_set());
    obj_in->remove_face(f2->position_in_set());
    return true;
  }

  void bevel(surf_ptr obj_in, const int &i1, const T &offset, const T &inset) {
    face_ptr f1 = obj_in->face(i1);
    this->bevel_face(obj_in, f1, offset, inset);
  }

  bool bevel_face(surf_ptr obj_in, face_ptr f1, const T &offset,
                  const T &inset) {

    face_vertex_ptr fvb = f1->fbegin();
    face_vertex_ptr fve = fvb->prev();
    bool it1 = true;
    coordinate_type N = this->normal(fvb->face());
    coordinate_type cen = this->center(fvb->face());

    vector<vertex_ptr> nverts;
    vector<vertex_ptr> overts;

    vector<face_vertex_ptr> topfv;
    vector<face_vertex_ptr> botfv;

    vector<coordinate_type> in_cont;

    int i = 0;
    bool set_head = true;
    // iterate through the face and replace the old face vertices with the
    // new offset ones

    while (it1) {
      it1 = fvb != fve;
      coordinate_type c = this->coordinate(fvb);

      coordinate_type C = c - cen;
      C.normalize();
      coordinate_type cn = c + N * offset + inset * C;

      vertex_ptr nv = new vertex_type();
      this->coordinate(cn, nv);
      vertex_ptr ov = fvb->vertex();
      // push back the old vertices and add the new ofset vertices
      nverts.push_back(nv);
      overts.push_back(ov);

      // now we have to remove the old face vertex from the vertex set
      ov->remove_face_vertex(fvb);
      // replace the vertex with the new offset one
      fvb->vertex() = nv;
      // and add this fv to the vertex, and in turn add the new vertex to the
      // set
      nv->add_face_vertex(fvb);
      obj_in->push_vertex(nv);
      // move through to the next face vertex
      fvb = fvb->next();
    }

    // now we loop around the face again and we have all these dangly edges that
    // still point to the offset face, these edges need to be replaced with
    // these dangly edges need to be replaced with faces, or rather we add
    // new faces and attach those edges to the new faces

    fvb = f1->fbegin();
    i = 0;
    it1 = true;
    while (it1) {
      it1 = fvb != fve;
      vertex_ptr nv0, nv1, nv2, nv3;
      face_vertex_ptr fv0, fv1, fv2, fv3;
      // now we grab the top and bottom vertices for this face,
      // we're circulating around the lists
      nv3 = overts[i];
      nv2 = nverts[i];
      nv1 = nverts[(i + 1) % nverts.size()];
      nv0 = overts[(i + 1) % nverts.size()];

      // now allocate a face and four new face vertices
      face_ptr nf = new face_type();
      fv0 = new face_vertex_type();
      fv1 = new face_vertex_type();
      fv2 = new face_vertex_type();
      fv3 = new face_vertex_type();
      nf->fbegin() = fv0; // set this new face to begin at fv0;

      // we connect up our new face so it circulates properly
      fv0->next() = fv1;
      fv1->prev() = fv0;
      fv1->next() = fv2;
      fv2->prev() = fv1;
      fv2->next() = fv3;
      fv3->prev() = fv2;
      fv3->next() = fv0;
      fv0->prev() = fv3;

      /*now we grab the dangly wangly edgy thingy
        and we set it to point to a new face vertex.
        this way the edge now points to the bottom of the
        new face rather than the old face
        _    _
        |\|->|_|
        \
      */
      edge_ptr ed = fvb->edge();
      if (ed->v1() == fvb) {
        ed->v1() = fv3;
      } else {
        ed->v2() = fv3;
      }

      fv3->edge() = ed;

      // now we add these new face vertices to their respective vertices.
      fv0->vertex() = nv0;
      nv0->add_face_vertex(fv0);
      fv0->face() = nf;
      fv1->vertex() = nv1;
      nv1->add_face_vertex(fv1);
      fv1->face() = nf;
      fv2->vertex() = nv2;
      nv2->add_face_vertex(fv2);
      fv2->face() = nf;
      fv3->vertex() = nv3;
      nv3->add_face_vertex(fv3);
      fv3->face() = nf;
      nf->size() = 4;
      // in order to connect up the sides, we'll have to loop around the top and
      // bottom vertices
      topfv.push_back(fv2);
      topfv.push_back(fv1);
      botfv.push_back(fv3);
      botfv.push_back(fv0);

      obj_in->push_face(nf);
      fvb = fvb->next();
      i++;
    }

    int sz = topfv.size();
    for (int i = 0; i < sz; i += 2) {
      edge_ptr ne = new edge_type();
      face_vertex_ptr fvtp = topfv[(i)];
      face_vertex_ptr fvbm = botfv[(sz + i - 1) % sz];
      ne->v1() = fvtp;
      fvtp->edge() = ne;
      ne->v2() = fvbm;
      fvbm->edge() = ne;
      obj_in->push_edge(ne);
    }

    fvb = f1->fbegin();
    i = 0;
    it1 = true;
    while (it1) {
      it1 = fvb != fve;
      face_vertex_ptr fvo = topfv[(i + 1) % sz];
      edge_ptr ne = new edge_type();
      ne->v2() = fvb;
      fvb->edge() = ne;
      ne->v1() = fvo;
      fvo->edge() = ne;
      fvb = fvb->next();
      i += 2;
      obj_in->push_edge(ne);
    }
    return true;
  }

  bool incremental_patch_hole(face_vertex_ptr &fv, surf_ptr obj) {
    face_vertex_ptr fvb = fv;
    face_vertex_ptr fve = fv->vprev();
    if (!fve) {
      fv->face()->color.r = 0;
      fv->face()->color.g = 0;
      fv->face()->color.b = 1;
      // obj->remove_face(fvb->face()->position_in_set());
      return 0;
    }
    if (!fvb->next()->vnext()) {
      fvb->face()->color.r = 0;
      fvb->face()->color.g = 1;
      fvb->face()->color.b = 0;
      return 0;
    }
#if 1
    bool atHead = false;
    bool isLoop = false;
    int maxIterations = 1000;
    int loopSize = 0;
    int iterations = 0;
    face_vertex_ptr fv0 = fvb;
    ;
    face_vertex_ptr fv1;
    face_vertex_ptr fv2;
    if (!fv0->next()->edge())
      fv1 = fv0->next();
    else
      fv1 = fv0->next()->vnext();
    if (!fv1->next()->edge())
      fv2 = fv1->next();
    else
      fv2 = fv1->next()->vnext();

    face_ptr nf = new face_type();
    face_vertex_ptr nfv0 = new face_vertex_type();
    nfv0->face() = nf;
    nfv0->vertex() = fv0->vertex();
    fv0->vertex()->add_face_vertex(nfv0);

    face_vertex_ptr nfv1 = new face_vertex_type();
    nfv1->face() = nf;
    nfv1->vertex() = fv1->vertex();
    fv1->vertex()->add_face_vertex(nfv1);

    face_vertex_ptr nfv2 = new face_vertex_type();
    nfv2->face() = nf;
    nfv2->vertex() = fv2->vertex();
    fv2->vertex()->add_face_vertex(nfv2);

    edge_ptr ne0 = new edge_type(*fv0->next(), *nfv0);
    obj->push_edge(ne0);
    edge_ptr ne1 = new edge_type(*fv1->next(), *nfv1);
    obj->push_edge(ne1);

    nfv0->next() = nfv1;
    nfv1->prev() = nfv0;
    nfv1->next() = nfv2;
    nfv2->prev() = nfv1;
    nfv2->next() = nfv0;
    nfv0->prev() = nfv2;

    nfv0->vertex()->color.r = 1.0;
    nfv0->vertex()->color.g = 0.0;
    nfv0->vertex()->color.b = 0.0;

    nf->color.r = 1.0;
    nf->color.g = 0.0;
    nf->color.b = 0.0;

    nf->fbegin() = nfv0;
    nf->update_all();
    obj->push_face(nf);

#endif
  }

  bool patch_hole(face_vertex_ptr fv, surf_ptr obj) {
    face_vertex_ptr fvb = fv;
    face_vertex_ptr fve = fv->vprev();
    if (!fve) {
      fv->face()->color.r = 0;
      fv->face()->color.g = 0;
      fv->face()->color.b = 1;
      // obj->remove_face(fvb->face()->position_in_set());
      return 0;
    }
    if (!fvb->next()->vnext()) {
      fvb->face()->color.r = 0;
      fvb->face()->color.g = 1;
      fvb->face()->color.b = 0;
      return 0;
    }
#if 1
    bool atHead = false;
    bool isLoop = false;
    int maxIterations = 1000;
    int loopSize = 0;
    int iterations = 0;
    while (!atHead && iterations < maxIterations) {
      if (!fvb->edge()) {
        // std::cout << fvb << " " << fvb->next()->vnext() << " " <<
        // fvb->vertex() <<" " << fve  << std::endl;
        fvb->face()->color.r = 0;
        fvb->face()->color.g = 0.5;
        fvb->face()->color.b = 0.5;
        if (!fvb->next()->edge())
          fvb = fvb->next();
        else
          fvb = fvb->next()
                    ->vnext(); // i wish I could draw pictures of all of this.
        loopSize++;
      } else
        fvb = fvb->vnext();
      atHead = fvb == fve;
      iterations++;
    }

    if (iterations < maxIterations && loopSize > 1) {
      std::cout << "found unclosed loop! loop size: " << loopSize << std::endl;
      face_ptr nf = new face_type();
      vector<face_vertex_ptr> nfvs;
      vector<face_vertex_ptr> ofvs;
      fvb = fv;
      fve = fv->vprev();

      atHead = false;
      while (!atHead) {
        fvb->face()->color.r = 0.5;
        fvb->face()->color.g = 0.5;
        fvb->face()->color.b = 0.0;
        std::cout << " fvb " << fvb << " has edge: " << fvb->edge()
                  << std::endl;
        face_vertex_ptr nfv = new face_vertex_type();
        nfv->face() = nf;
        nfv->vertex() = fvb->vertex();
        fvb->vertex()->add_face_vertex(nfv);
        nfvs.push_back(nfv);
        ofvs.push_back(fvb);
        if (!fvb->next()->edge())
          fvb = fvb->next();
        else {
          fvb = fvb->next(); // i wish I could draw pictures of all of this.
          bool hasEdge = true;
          while (hasEdge) {
            fvb = fvb->vnext();
            hasEdge = fvb->edge() != NULL;
          }
        }
        atHead = fvb == fv;
      }
      if (nfvs.size() > 1) {
        // std::cout << "done collecting face vertices" << std::endl;
        int N = nfvs.size();
        for (int i = 0; i < nfvs.size(); i++) {
          int in = (i + 1 < N) ? i + 1 : N - i - 1;
          int ip = (i - 1 >= 0) ? i - 1 : N - 1 + i;
          std::cout << ip << " " << i << " " << in << " " << std::endl;
          nfvs[i]->next() = nfvs[ip];
          nfvs[ip]->prev() = nfvs[i];
          // nfvs[i]->prev()  = nfvs[ip];
          // nfvs[ip]->next() = nfvs[i];

          edge_ptr ne = new edge_type(nfvs[i], ofvs[ip]);

          obj->push_edge(ne);

          std::cout << ne->position_in_set() << " "
                    << ne->fv1->vertex()->position_in_set() << " "
                    << ne->fv2->vertex()->position_in_set() << " "
                    << ne->fv1->next()->vertex()->position_in_set() << " "
                    << ne->fv2->next()->vertex()->position_in_set() << " "
                    << ne->fv1 << " " << ne->fv2 << " " << std::endl;
        }
        nf->fbegin() = nfvs[0];
        nf->update_all();
        obj->push_face(nf);
        // std::cout << "done patching face" << std::endl;
        return 1;
      } else
        return 0;
    } else
      return 0;

#endif
  }

  void fill_holes(surf_ptr obj, vector<list<face_vertex_ptr>> corners) {
    vector<vertex_ptr> &verts = obj->get_vertices();
    int kk = 0;
    int max = 40;
    for (int i = 0; i < verts.size(); i++) {
      // verts[i]->pack();
      list<face_vertex_ptr> &rFaceVertices = corners[i];
      fvl_iterator jtb = rFaceVertices.begin();
      fvl_iterator jte = rFaceVertices.end();
      for (jtb; jtb != jte; jtb++) {
        face_vertex_ptr fvj = *jtb;
        if (fvj && fvj->edge() == NULL) {
          fvj->face()->color.r = 1.0;
          fvj->face()->color.g = 0.0;
          fvj->face()->color.b = 0.0;
          bool patch = false;
          patch = patch_hole(fvj, obj);
          if (patch == 1) {
            kk++;
          }
          // if(kk > max) {std::cout << " returning" << std::endl; return;};
        }
      }
    }
  }

  void incremental_fill_holes(surf_ptr obj) {
    vector<vertex_ptr> &verts = obj->get_vertices();
    int kk = 0;
    int max = 40;
    int edgeLessfvs = 1;
    while (edgeLessfvs > 0 && kk < 1) {
      for (int i = 0; i < verts.size(); i++) {
        verts[i]->pack();
        vector<face_vertex_ptr> &rFaceVertices = verts[i]->get_face_vertices();
        for (int j = 0; j < rFaceVertices.size(); j++) {
          face_vertex_ptr fvj = rFaceVertices[j];
          if (fvj && fvj->edge() == NULL) {
            fvj->face()->color.r = 1.0;
            fvj->face()->color.g = 0.0;
            fvj->face()->color.b = 0.0;
            bool patch = false;
            patch = incremental_patch_hole(fvj, obj);
            if (patch == 1) {
            }
          }
        }
      }

      edgeLessfvs = 0;
      for (int i = 0; i < verts.size(); i++) {
        verts[i]->pack();
        vector<face_vertex_ptr> &rFaceVertices = verts[i]->get_face_vertices();
        for (int j = 0; j < rFaceVertices.size(); j++) {
          face_vertex_ptr fvj = rFaceVertices[j];
          if (fvj && fvj->edge() == NULL) {
            edgeLessfvs++;
          }
        }
      }
      kk++;
      std::cout << "remaining edges: " << edgeLessfvs << std::endl;
    }
  }
};

} // namespace m2
#endif
