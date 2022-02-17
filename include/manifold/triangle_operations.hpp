//
//  modify.hpp
//  Manifold
//
//  Created by John Delaney on 5/21/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#ifndef __M2TRIOPS__
#define __M2TRIOPS__
#include "m2Includes.h"
#include "manifold/m2.hpp"
#include <algorithm>
#include <cmath>

namespace m2 {

template <typename SPACE>

// forward declaration
class edge_collapser;

template <typename SPACE> class edge_split {
  M2_TYPEDEFS;

public:
  vertex_ptr operator()(surf_ptr obj_in, edge_ptr edge_in) {

    face_vertex_ptr fv1 = edge_in->v1();
    face_vertex_ptr fv2 = edge_in->v2();

    vertex_ptr v1 = fv1->vertex();
    vertex_ptr v2 = fv2->vertex();

    vertex_ptr vn = new vertex_type();

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
};

template <typename SPACE> class edge_collapse {
  M2_TYPEDEFS;

public:
#if 1
  vertex_ptr operator()(surf_ptr obj, edge_ptr e) {
    // TODO:: move to triangle_ops, compare with joinFaceOneSharedEdge
    if (e->v1()->vertex() == e->v2()->vertex())
      return e->v1()->vertex();

    face_vertex_ptr va0 = e->v1();
    face_vertex_ptr vb0 = e->v2();
    vertex_ptr v1 = va0->vertex();
    vertex_ptr v2 = vb0->vertex();
    construct<SPACE> cons;
    if (va0->face()->size() == 1 || va0->face()->size() == 2) {
      std::cout << "e pos (v2): " << e->position_in_set() << std::endl;
      v2->print();
      cons.delete_edge(obj, e);
      return v1;
    }

    if (vb0->face()->size() == 1 || vb0->face()->size() == 2) {
      std::cout << "e pos (v2): " << e->position_in_set() << std::endl;
      v2->print();
      cons.delete_edge(obj, e);
      return v2;
    }

    if (v1->calc_size() == 1 || v1->calc_size() == 2) {
      std::cout << "e pos (v1): " << e->position_in_set() << std::endl;
      v1->print();
      // std::cout << " deleting v1" << std::endl;
      face_ptr f = cons.delete_vertex_primitive(obj, v1);
      return v2;
    }

    if (v2->calc_size() == 1 || v2->calc_size() == 2) {
      std::cout << "e pos (v2): " << e->position_in_set() << std::endl;
      v2->print();
      face_ptr f = cons.delete_vertex_primitive(obj, v2);
      //	std::cout << "output 3: " << std::endl;
      // v1->print();
      return v1;
    }

    if (v1->calc_size() == 1 || v1->calc_size() == 2) {
      std::cout << "e pos (v1): " << e->position_in_set() << std::endl;
      v1->print();
      // std::cout << " deleting v1" << std::endl;
      face_ptr f = cons.delete_vertex_primitive(obj, v1);
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

    m2::for_each_vertex<SPACE>(v1, [&v1](face_vertex_ptr fv) -> void {
      fv->vertex() = v1;
      fv->next()->vertex()->front() = fv->next();
    });

    // while(this->delete_degenerates(obj,v1));
    return v1;
  }
#else
  vertex_ptr operator()(surf_ptr obj, edge_ptr e) {
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

    construct<SPACE> cons;
    if (fv1->face()->size() == 3) {
      std::cout << "dm1 " << fv2->next()->edge()->flag[0] << std::endl;
      cons.delete_edge(obj, fv2->next()->edge());
    }
    if (fv2->face()->size() == 3) {
      std::cout << "dm2" << fv2->next()->edge()->flag[0] << std::endl;
      cons.delete_edge(obj, fv2->next()->edge());
    }

    v1->front() = fv1->vprev();
    v2->front() = fv2->vprev();

    std::cout << "v1: " << v1->front() << " " << fv1->vnext() << std::endl;
    std::cout << "v2: " << v2->front() << " " << fv2->vnext() << std::endl;

    std::cout << v1->position_in_set() << ":" << v1->size() << " "
              << v2->position_in_set() << ":" << v2->size() << " " << std::endl;

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

    // bool iterating = true;
    // face_vertex_ptr fvb = fv2->vnext();
    // face_vertex_ptr fve = fv2->vprev();

    v1->update_all();
    v2->update_all();
    /*
    while (iterating) {
      iterating = fvb != fve;
      v2->remove_face_vertex(fvb);
      v1->add_face_vertex(fvb);
      fvb = fvb->vnext();
    }

    v1->front() = fvb;
    */
    std::cout << e << " " << fv1 << " " << fv2 << std::endl;

    std::cout << "d0" << std::endl;
    v1->remove_face_vertex(fv1);
    delete fv1;

    std::cout << "d1" << std::endl;
    v2->remove_face_vertex(fv2);
    delete fv2;

    std::cout << "update_v1"
              << " " << v1->front() << std::endl;
    v1->update_all();

    std::cout << "update_v2"
              << " " << v2->front() << std::endl;
    v2->update_all();
    // std::cout << "output: " << std::endl;
    // v1->print();
    // v2->print();

    std::cout << "d2" << std::endl;
    obj->remove_edge(e->position_in_set());

    std::cout << "d3" << std::endl;
    obj->remove_vertex(v2->position_in_set());

    return v1;
  }
#endif
};

template <typename SPACE> class edge_flip {
  M2_TYPEDEFS;

public:
  edge_flip() {}
  ~edge_flip() {}

  edge_ptr operator()(edge_ptr e1) {
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

template <typename SPACE>
class triangle_operations_base : public default_interface<SPACE> {
public:
  M2_TYPEDEFS;

  triangle_operations_base(const surf_ptr &surf, real max)
      : default_interface<SPACE>(), _thresh(max), _surf(surf) {}

  void reset_flags() {
    //_surf->reset_flags();

    edge_array edges = _surf->get_edges();
    for (int i = 0; i < edges.size(); i++) {
      if (!edges[i])
        continue;
      edges[i]->flags[0] = 0;
      edges[i]->flags[1] = 0;
    }

    vertex_array vertices = _surf->get_vertices();
    for (int i = 0; i < vertices.size(); i++) {
      if (!vertices[i])
        continue;
      vertices[i]->topologyChange = -1;
    }

    std::cout << "done:" << std::endl;
  }

  struct {
    bool operator()(edge_ptr ei, edge_ptr ej) {
      return (ci::length<SPACE>(ei) > ci::length<SPACE>(ej));
    }
  } mEdgeSorterGreater;

  struct {
    bool operator()(edge_ptr ei, edge_ptr ej) {
      return (ci::length<SPACE>(ei) > ci::length<SPACE>(ej));
    }
  } mEdgeSorterLesser;

  surf_ptr _surf;
  real _thresh = 0.1;
};

template <typename SPACE>
class edge_flipper : public triangle_operations_base<SPACE> {
  M2_TYPEDEFS;

public:
  edge_flipper(const surf_ptr &surf)
      : triangle_operations_base<SPACE>(surf, 0.0) {}

  bool flip() {
    // TIMER function//TIMER(__FUNCTION__);
    edge_array &edges = this->_surf->get_edges();

    bool flipped = false;
    edge_array permEdges;
    for (int i = 0; i < edges.size(); i++) {
      if (!this->_surf->has_edge(i))
        continue;
      edge_ptr e = edges[i];
      if (e->v1()->face()->size() == 3 && e->v2()->face()->size() == 3)
        permEdges.push_back(edges[i]);
    }

    for (int i = 0; i < permEdges.size(); i++) {
      int card = rand() % permEdges.size();
      edge_ptr et = permEdges[i];
      permEdges[i] = permEdges[card];
      permEdges[card] = et;
    }

    for (int i = 0; i < permEdges.size(); i++) {
      // std::cout << "A" << std::endl;
      edge_ptr e = permEdges[i];
      bool pinned = e->v1()->vertex()->pinned == true &&
                    e->v2()->vertex()->pinned == true;
      bool notPinned = e->v1()->vertex()->pinned != true &&
                       e->v2()->vertex()->pinned != true;
      // std::cout << "B" << std::endl;
      if (pinned)
        continue;
      if (notPinned)
        continue;

      // if(e->v1()->vertex()->size() < 4) continue;
      // if(e->v2()->vertex()->size() < 4) continue;

      face_vertex_ptr v0 = e->v1();
      face_vertex_ptr v1 = v0->prev();
      face_vertex_ptr v2 = e->v2();
      face_vertex_ptr v3 = v2->prev();

      coordinate_type c0 = this->coordinate(v0);
      coordinate_type c1 = this->coordinate(v1);
      coordinate_type c2 = this->coordinate(v2);
      coordinate_type c3 = this->coordinate(v3);

      face_vertex_ptr fvEdge = NULL;
      T m01 = 1.0 / (c0 - c1).norm();
      T m12 = 1.0 / (c1 - c2).norm();
      T m23 = 1.0 / (c2 - c3).norm();
      T m30 = 1.0 / (c3 - c0).norm();

      T cos0 = (c1 - c0).dot(c3 - c0) * m01 * m30;
      T cos1 = (c0 - c1).dot(c2 - c1) * m01 * m12;
      T cos2 = (c1 - c2).dot(c3 - c2) * m12 * m23;
      T cos3 = (c0 - c3).dot(c2 - c3) * m30 * m23;
      // half angle cos^2(2a) = 0.5*(1+cos(a))
      T cSame = acos(cos1) + acos(cos3); // corresponds to flipped edge
      T cFlip = acos(cos0) + acos(cos2); // corresponds to flipped edge
      // T cSame =  (M_PI - (acos(cos1) + acos(cos3))); //corresponds to flipped
      // edge

      T eFlip = cFlip * cFlip;
      T eSame = cSame * cSame;

      // if(cSame > M_PI){
      if (cFlip < cSame) {
        edge_flip<SPACE> flip;
        flip(e);
        // e->v1()->data *= -1;
        // e->v2()->data *= -1;
        flipped = true;
      }
    }
    return flipped;
  }

}; // set_operations

#if 1
template <typename SPACE>
class face_merger : public triangle_operations_base<SPACE> {
public:
  M2_TYPEDEFS;

  face_merger(const surf_ptr &surf, real max)
      : triangle_operations_base<SPACE>(surf, max) {}

  void joinFaces(face_ptr fA, face_ptr fB, surf_ptr mesh) {
    // std::cout << "no shared_edge: " << fA->position_in_set() << " "
    //          << fB->position_in_set() << std::endl;
    face_vertex_ptr fvb0 = fB->fbegin()->prev();
    std::vector<face_vertex_ptr> fv_to_delete;
    std::vector<edge_ptr> e_to_delete;

    face_vertex_ptr fvA0 = fA->fbegin();
    face_vertex_ptr fvA1 = fvA0->next();
    face_vertex_ptr fvA2 = fvA1->next();

    face_vertex_ptr fvB0 = fB->fbegin()->prev();
    face_vertex_ptr fvB1 = fvB0->prev();
    face_vertex_ptr fvB2 = fvB1->prev();

    vertex_ptr vA0 = fvA0->vertex();
    vertex_ptr vA1 = fvA1->vertex();
    vertex_ptr vA2 = fvA2->vertex();
    vertex_ptr vB0 = fvB0->vertex();
    vertex_ptr vB1 = fvB1->vertex();
    vertex_ptr vB2 = fvB2->vertex();

    std::cout << " vA0: " << vA0->position_in_set() << " "
              << vB0->position_in_set() << std::endl;
    std::cout << " vA1: " << vA1->position_in_set() << " "
              << vB1->position_in_set() << std::endl;
    std::cout << " vA2: " << vA2->position_in_set() << " "
              << vB2->position_in_set() << std::endl;

    m2::for_each_face<SPACE>(fA, [&](face_vertex_ptr fva0) {
      face_vertex_ptr fvb1 = fvb0->next();

      vertex_ptr vA0 = fva0->vertex();
      vertex_ptr vB0 = fvb0->vertex();
      vertex_ptr vB1 = fvb1->vertex();

      edge_ptr ea = fva0->edge();
      edge_ptr eb = fvb0->edge();

      if (ea != eb) {
        e_to_delete.push_back(eb);
        vA0->front() = fvB0->coedge();
        vB0->remove_face_vertex(fva0);
        vA0->remove_face_vertex(fvb0);

        fv_to_delete.push_back(fva0);
        fv_to_delete.push_back(fvb0);
        std::cout << fva0 << " " << fva0->coedge() << std::endl;
        std::cout << fvb0 << " " << fvb0->coedge() << std::endl;
        ea->set(fva0->coedge(), fvb0->coedge());
      }

      fvb0 = fvb0->prev();
    });

    edge_ptr eA0 = fvA0->edge();
    edge_ptr eA1 = fvA1->edge();
    edge_ptr eA2 = fvA2->edge();

    vA0->update_all();
    vA1->update_all();
    vA2->update_all();

    vA0->print();
    vA1->print();
    vA2->print();

    std::for_each(e_to_delete.begin(), e_to_delete.end(), [&mesh](edge_ptr e) {
      mesh->remove_edge(e->position_in_set());
    });
    std::for_each(fv_to_delete.begin(), fv_to_delete.end(),
                  [&mesh](face_vertex_ptr fv) { delete fv; });

    mesh->remove_face(fA->position_in_set());
    mesh->remove_face(fB->position_in_set());

    std::cout << "vA0 ver" << std::endl;
    vA0->verify();
    std::cout << "vA1 ver" << std::endl;
    vA1->verify();
    std::cout << "vA2 ver" << std::endl;
    vA2->verify();

    eA0->v1()->face()->verify();
    eA0->v2()->face()->verify();
    eA1->v1()->face()->verify();
    eA1->v2()->face()->verify();
    eA2->v1()->face()->verify();
    eA2->v2()->face()->verify();
  }

  std::vector<edge_ptr> get_bridged_edges(face_ptr fA, face_ptr fB) {
    std::vector<edge_ptr> collected;

    m2::set_vertex_flag<SPACE>(fA, 0, 0);
    m2::set_vertex_flag<SPACE>(fB, 0, 1);

    std::cout << std::endl;
    m2::for_each_face<SPACE>(fA, [&](face_vertex_ptr fv) {
      vertex_ptr v = fv->vertex();

      m2::for_each_vertex<SPACE>(v, [&](face_vertex_ptr fvv) {
        if (fvv->coedge()->vertex()->flags[0]) {
          collected.push_back(fvv->edge());
        }
      });
      std::cout << std::endl;
    });

    m2::set_vertex_flag<SPACE>(fA, 0, 0);
    m2::set_vertex_flag<SPACE>(fB, 0, 0);
    return collected;
  }

  void alignFaces(face_ptr fA, face_ptr fB) {
    fA->update_all();
    fB->update_all();

    face_vertex_ptr fvA0 = fA->fbegin();
    face_vertex_ptr fvA1 = fvA0->next();
    face_vertex_ptr fvA2 = fvA1->next();

    face_vertex_ptr fvB0 = fB->fbegin();
    face_vertex_ptr fvB1 = fvB0->prev();
    face_vertex_ptr fvB2 = fvB1->prev();

    /*
    vertex_ptr vA0 = fvA0->vertex();
    vertex_ptr vA1 = fvA1->vertex();
    vertex_ptr vA2 = fvA2->vertex();

    vertex_ptr vB0 = fvB0->vertex();
    vertex_ptr vB1 = fvB1->vertex();
    vertex_ptr vB2 = fvB2->vertex();
    */
    coordinate_type ca0 = this->coordinate(fvA0);
    coordinate_type ca1 = this->coordinate(fvA1);
    coordinate_type ca2 = this->coordinate(fvA2);
    // assume the other coordinate rotates in opposite direction since they
    // are facing opposite directions
    coordinate_type cb0 = this->coordinate(fvB0);
    coordinate_type cb1 = this->coordinate(fvB1);
    coordinate_type cb2 = this->coordinate(fvB2);

    T d0 = 1.0 / 3.0 *
           ((cb0 - ca0).squaredNorm() + (cb1 - ca1).squaredNorm() +
            (cb2 - ca2).squaredNorm());
    T d1 = 1.0 / 3.0 *
           ((cb0 - ca1).squaredNorm() + (cb1 - ca2).squaredNorm() +
            (cb2 - ca0).norm());
    T d2 = 1.0 / 3.0 *
           ((cb0 - ca2).squaredNorm() + (cb1 - ca0).squaredNorm() +
            (cb2 - ca1).squaredNorm());

    // std::cout << " distance: " << d0 << " " << d1 << " " << d2 << std::endl;

    face_vertex_ptr fvAb = fvA0;
    T db = d0;

    fvAb = (db < d1) ? fvAb : fvA1;
    db = (db < d1) ? db : d1;
    fvAb = (db < d2) ? fvAb : fvA2;
    db = (db < d2) ? db : d2;

    fA->front() = fvAb;
    fB->front() = fvB0;
  }

  void averageVerts(face_ptr fA, face_ptr fB) {
    face_vertex_ptr fvA = fA->fbegin();
    face_vertex_ptr fvAe = fA->fend();
    face_vertex_ptr fvB = fB->fbegin();
    bool itf = true;
    while (itf) {
      itf = fvA != fvAe;
      vertex_ptr vA = fvA->vertex();
      vertex_ptr vB = fvB->vertex();
      std::cout << vA->position_in_set() << " " << vB->position_in_set()
                << std::endl;
      if (vA != vB) {
        coordinate_type cA = this->coordinate(vA);
        coordinate_type cB = this->coordinate(vB);
        coordinate_type c = 0.5 * (cA + cB);
        this->coordinate(c, vA);
      }
      fvA = fvA->next();
      fvB = fvB->prev();
    }
  }

  int min_valence(face_ptr f) {
    int mv = 999;
    m2::for_each_face<SPACE>(f, [&mv](face_vertex_ptr fv) {
      typename SPACE::real li = 0.0;
      int v = fv->vertex()->calc_size();
      mv = v < mv ? v : mv;
    });
    return mv;
  }

  std::vector<edge_ptr> get_shared_edges(face_ptr fA, face_ptr fB) {
    std::vector<edge_ptr> shared;
    face_vertex_ptr fvA = fA->fbegin();
    face_vertex_ptr fvAe = fA->fend();
    bool itA = true;
    while (itA) {
      itA = fvA != fvAe;
      edge_ptr eA = fvA->edge();
      face_vertex_ptr fvB = fB->fbegin();
      face_vertex_ptr fvBe = fB->fend();
      bool itB = true;
      while (itB) {
        itB = fvB != fvBe;
        edge_ptr eB = fvB->edge();
        if (eA == eB) {
          shared.push_back(eB);
        };
        fvB = fvB->next();
      }
      fvA = fvA->next();
    }
    return shared;
  }

  vector<triangle_pair> get_face_pairs_to_merge(std::vector<face_ptr> &faces) {
    std::vector<triangle_type> tris = m2::ci::get_tris<SPACE>(
        faces); // get tris is a good thing to pass to namespace ci
    m2::aabb_tree<SPACE, triangle_type> tree(tris);

    std::cout << "tree.nodes.size(): " << tree.nodes.size() << std::endl;

    auto triDist = [](const triangle_type &A, const triangle_type &B) -> T {
      return A.dist(B);
    };

    T tol = this->_thresh;

    vector<triangle_pair> collected;
    // for (int i = 0; i < 1; i++) {
    // std::cout << " foo !" << std::endl;
    for (int i = 0; i < tris.size(); i++) {
      triangle_type triNear =
          m2::getNearest<SPACE, triangle_type, triangle_type>(
              tris[i], tree, tris, triDist, tol);
      if (triNear.faceId > -1) {
        collected.push_back(triangle_pair(tris[i], triNear));
      }
    }

    std::sort(collected.begin(), collected.end(), std::less<triangle_pair>());

    std::cout << " full list: " << collected.size() << std::endl;
    auto it = std::unique(collected.begin(), collected.end());

    collected.erase(it, collected.end());
    std::cout << " unique list: " << collected.size() << std::endl;

    collected.erase(std::remove_if(collected.begin(), collected.end(),
                                   [faces, tol](auto p) {
                                     face_ptr fA = faces[p.A.faceId];
                                     face_ptr fB = faces[p.B.faceId];

                                     T d = p.dist();
                                     if (d > tol)
                                       return true;

                                     if (fA->flags[0] > 0)
                                       return true;
                                     if (fB->flags[0] > 0)
                                       return true;

                                     fA->flags[0] = 1;
                                     fB->flags[0] = 1;
                                     return false;
                                   }),
                    collected.end());

    std::cout << " merging " << collected.size() << " tris!" << std::endl;
    return collected;
  }

  std::vector<triangle_pair>
  collapse_bridges(std::vector<triangle_pair> &collected,
                   std::vector<face_ptr> &faces) {
    std::vector<edge_ptr> edgeToCollapse;
    // collect bridges
    for (auto p : collected) {
      face_ptr fA = faces[p.A.faceId];
      face_ptr fB = faces[p.B.faceId];

      std::vector<edge_ptr> i_bridges = get_bridged_edges(fA, fB);

      // copy edges to collection if they haven't been already been flagged
      std::copy_if(i_bridges.begin(), i_bridges.end(),
                   std::back_inserter(edgeToCollapse), [](edge_ptr e) {
                     // this is nice symmetry that lets you skip every other edge
                     bool insert = !(e->flags[0] | e->v1()->face()->flags[1] |
                                     e->v2()->face()->flags[1]);
                     if (insert) {
                       e->v1()->face()->flags[1] = 1;
                       e->v2()->face()->flags[1] = 1;
                       e->flags[0] = 1;
                     }

                     return insert;
                   });
    }

    std::vector<triangle_pair> newCollect;
    std::copy_if(collected.begin(), collected.end(),
                 std::back_inserter(newCollect), [faces](triangle_pair p) {
                   std::cout << "A/B id: " << p.A.faceId << " " << p.B.faceId
                             << std::endl;
                   face_ptr fA = faces[p.A.faceId];
                   face_ptr fB = faces[p.B.faceId];
                   return (!(fA->flags[1] | fB->flags[1]));
                 });

    // clear flags
    std::for_each(edgeToCollapse.begin(), edgeToCollapse.end(),
                  [](edge_ptr &e) {
                    e->flags[0] = 0;
                    e->v1()->face()->flags[1] = 0;
                    e->v2()->face()->flags[1] = 0;
                  });

    std::vector<coordinate_type> coords(edgeToCollapse.size());
    int i = 0;
    for (auto e : edgeToCollapse) {
      vertex_ptr v0 = e->v1()->vertex();
      vertex_ptr v1 = e->v2()->vertex();
      coordinate_type c0 = m2::ci::get_coordinate<SPACE>(v0);
      coordinate_type c1 = m2::ci::get_coordinate<SPACE>(v1);
      coords[i] = 0.5 * c0 + 0.5 * c1;
      i++;
    }

    m2::edge_collapser<SPACE> collapser(this->_surf, 0.0);

    std::cout << "shared edges: " << edgeToCollapse.size() << std::endl;
    collapser.collapse_edges(edgeToCollapse);
    vertex_array vertices = this->_surf->get_vertices();

    std::cout << "assigning calc'd vertex values" << std::endl;

    for (auto v : vertices) {
      if (!v)
        continue;
      int topologyChange = v->topologyChange;
      if (topologyChange < 0)
        continue;
      m2::ci::set_coordinate<SPACE>(coords[topologyChange], v);
    }

    std::cout << "new collect size: " << newCollect.size() << std::endl;
    return newCollect;
  }

  bool merge_collected(std::vector<triangle_pair> &collected,
                       std::vector<face_ptr> &faces) {

    bool merging = true;
    int phase = 0;
    while (merging && collected.size()) {

      collected = collapse_bridges(collected, faces);

      std::sort(collected.begin(), collected.end(), std::less<triangle_pair>());

      std::cout << "merging " << phase
                << ", triangle pair count: " << collected.size() << std::endl;

      std::vector<triangle_pair> new_collect;

      for (auto p : collected) {

        face_ptr fA = faces[p.A.faceId];
        face_ptr fB = faces[p.B.faceId];

        int shared_vertices = fA->count_shared_vertices(fB);
        int shared_edges = fA->count_shared_edges(fB);

        // std::cout << " shared: " << shared << std::end
        coordinate_type NA = ci::normal<SPACE>(fA);
        coordinate_type NB = ci::normal<SPACE>(fB);
        real angle = va::dot(NA, NB);
        T d = p.dist();
        std::cout << "dist/angle: " << d << " / " << angle << std::endl;

        alignFaces(fA, fB);
        std::cout << "popping hole " << shared_vertices << std::endl;
        averageVerts(fA, fB);
        fA->flag_edges(2);
        fB->flag_edges(2);
        joinFaces(fA, fB, this->_surf);
        merging = true;
      }
      collected = new_collect;
    }

    this->_surf->verify();
    this->_surf->update_all();
    this->_surf->pack();

    this->reset_flags();
    return merging;
  }

  bool merge() {
    std::vector<face_ptr> faces = this->_surf->get_faces();
    vector<triangle_pair> collected = get_face_pairs_to_merge(faces);
    return this->merge_collected(collected, faces);
  }

  void reset_flags() {

    for (auto f : this->_surf->get_faces()) {
      f->flags[0] = 0;
    }

    for (auto v : this->_surf->get_vertices()) {
      v->flags[0] = 0;
      v->flags[1] = 0;
    }

    for (auto e : this->_surf->get_edges()) {
      e->flags[2] = 0;
    }
  }
};
#endif

#if 1
template <typename SPACE>
class edge_splitter : public triangle_operations_base<SPACE> {
public:
  M2_TYPEDEFS;

  edge_splitter(const surf_ptr &surf, real max)
      : triangle_operations_base<SPACE>(surf, max) {}

  edge_ptr subdivideFace(face_vertex_ptr fv) {

    // coordinate_type data = fv->face()->data;
    // coordinate_type data2 = fv->face()->data2;

    // TIMER function//TIMER(__FUNCTION__);
    // assuming that the face vertex is the newly inserted one.
    face_vertex_ptr fv1 = fv; // this is the
    face_vertex_ptr fv2 = fv->next();
    face_vertex_ptr fv3 = fv2->next();
    face_vertex_ptr fv4 = fv3->next();

    m2::construct<SPACE> cons;
    edge_ptr enew = cons.insert_edge(this->_surf, fv1, fv3);

    return enew;
  }

  void split_edge(int i, edge_ptr e) {
    // TIMER function//TIMER(__FUNCTION__);
    face_vertex_ptr fv1 = e->v1();
    face_vertex_ptr fv2 = e->v2();
    T circ1 = fv1->data;
    T circ2 = fv2->data;
    edge_split<SPACE> split;
    // coordinate_type c = getButterflyWeight(e);
    vertex_ptr nv = split(this->_surf, e);
    nv->topologyChange = i;
    // std::cout << "nsplit:" << i << " " << nv->position_in_set()
    //          << std::endl;

    // nv->coordinate() = c;
    fv1->edge()->flags[0] = 1;
    fv2->edge()->flags[0] = 1;

    edge_ptr e11 = subdivideFace(fv1->next());
    edge_ptr e21 = subdivideFace(fv2->next());

    e11->flags[1] = 0;
    e21->flags[1] = 0;
  }

  edge_array get_edges_to_split() {
    edge_array &edges = this->_surf->get_edges();

    edge_array edgesToSplit;
#if 1
    for (int i = 0; i < edges.size(); i++) {
      if (!this->_surf->has_edge(i))
        continue;
      edge_ptr ei = edges[i];
      bool pinned = ei->v1()->vertex()->pinned == true &&
                    ei->v2()->vertex()->pinned == true;
      if (pinned)
        continue;
      T l = this->length(ei);
      if (l > this->_thresh) {
        edgesToSplit.push_back(ei);
        continue;
      }
    }
#endif

    std::cout << " - sorting " << edgesToSplit.size() << " edges" << std::endl;
    std::sort(edgesToSplit.begin(), edgesToSplit.end(),
              this->mEdgeSorterGreater);
    return edgesToSplit;
  }

  bool split_edges(edge_array edgesToSplit) {
    // TIMER function//TIMER(__FUNCTION__);

    std::cout << "splitting " << edgesToSplit.size() << " edges" << std::endl;

    for (int i = 0; i < edgesToSplit.size(); i++) {
      this->split_edge(i, edgesToSplit[i]);
    }

    return edgesToSplit.size() > 0;
  }
};
#endif

#if 1
template <typename SPACE>
class edge_collapser : public triangle_operations_base<SPACE> {
public:
  M2_TYPEDEFS;

  edge_collapser(const surf_ptr &surf, real max)
      : triangle_operations_base<SPACE>(surf, max) {}

  edge_array get_edges_to_collapse() {
    edge_array &edges = this->_surf->get_edges();

    edge_array edgeToCollapse;
#if 1
    for (int i = 0; i < edges.size(); i++) {
      if (!this->_surf->has_edge(i))
        continue;
      edge_ptr ei = edges[i];

      face_vertex_ptr fv1 = ei->v1();
      face_vertex_ptr fv2 = ei->v2();
      fv1->next()->edge()->flags[0] = true;
      fv2->next()->edge()->flags[0] = true;

      bool pinned = ei->v1()->vertex()->pinned == true &&
                    ei->v2()->vertex()->pinned == true;
      if (pinned)
        continue;

      T l = this->length(ei);
      if (l < this->_thresh) {
        edgeToCollapse.push_back(ei);
        continue;
      }
    }
#endif

    std::cout << " - sorting " << edgeToCollapse.size() << " edges"
              << std::endl;
    std::sort(edgeToCollapse.begin(), edgeToCollapse.end(),
              this->mEdgeSorterLesser);

    edgeToCollapse.erase(std::remove_if(edgeToCollapse.begin(),
                                        edgeToCollapse.end(),
                                        [](edge_ptr e) { return e->flags[0]; }),
                         edgeToCollapse.end());

    return edgeToCollapse;
  }

  bool collapse_edges(edge_array edgeToCollapse) {
    // TIMER function//TIMER(__FUNCTION__);
    std::cout << " - deleting: " << edgeToCollapse.size() << " Tiny edges"
              << std::endl;

    for (int i = 0; i < edgeToCollapse.size(); i++) {
      edge_collapse<SPACE> collapser;
      std::cout << " i: " << i << std::endl;
      if (!edgeToCollapse[i])
        continue;
      edge_ptr e = edgeToCollapse[i];
      vertex_ptr nv = collapser(this->_surf, e);
      if (nv)
        nv->topologyChange = i;
    }
    this->_surf->pack();
    return edgeToCollapse.size() > 0;
  }
};
#endif
}; // namespace m2
#endif
