//
//  modify.hpp
//  Manifold
//
//  Created by John Delaney on 5/21/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#ifndef __M2TRIOPS__
#define __M2TRIOPS__
#include "construct.hpp"
#include "coordinate_interface.hpp"
#include "m2.hpp"

#include "manifold/asawa/remesh.hpp"

#include "manifold/bontecou/laplacian.hpp"

#include "manifold/bins.hpp"
#include "manifold/diffuse.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <set>
#include <vector>

namespace asawa {

template <typename SPACE>

// forward declaration
class edge_collapser;

template <typename SPACE> class edge_split {
  M2_TYPEDEFS;

public:
  vertex_ptr operator()(surf_ptr obj_in, edge_ptr edge_in, edge_ptr &e1,
                        edge_ptr &e2) {

    face_vertex_ptr fv1 = edge_in->v1();
    face_vertex_ptr fv2 = edge_in->v2();

    vertex_ptr v1 = fv1->vertex();
    vertex_ptr v2 = fv2->vertex();

    vertex_ptr vn = new vertex_type();

    e1 = edge_in;
    e2 = new edge_type();
    face_vertex_ptr fv1n = fv1->add_next();
    face_vertex_ptr fv2n = fv2->add_next();

    fv1->vertex()->set_front(fv1);
    fv2->vertex()->set_front(fv2);

    v1->remove_face_vertex(fv1n);
    v2->remove_face_vertex(fv2n);

    fv1n->face() = fv1->face();
    fv2n->face() = fv2->face();

    vn->add_face_vertex(fv1n);
    vn->add_face_vertex(fv2n);
    vn->set_front(fv1n);

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

    vn->update();
    v1->update();
    v2->update();
    return vn;
  }
};

template <typename SPACE> class edge_merge {
  M2_TYPEDEFS;

public:
  void operator()(surf_ptr obj_in, edge_ptr eA, edge_ptr eB) {
    alignEdges(eA, eB);

    if (!compatible(eA, eB))
      return;

    //    std::cout << " joining: " << eA->position_in_set() << " "
    //              << eB->position_in_set() << std::endl;

    if (eA->v1()->face() == eB->v1()->face())
      return;
    if (eA->v2()->face() == eB->v2()->face())
      return;

    if (eA->v2()->face() == eB->v1()->face())
      return;
    if (eA->v1()->face() == eB->v2()->face())
      return;

    face_vertex_ptr cA0 = eA->v1();
    face_vertex_ptr cA1 = eA->v2();
    face_vertex_ptr cB0 = eB->v1();
    face_vertex_ptr cB1 = eB->v2();

    vertex_ptr vA0 = cA0->vertex();
    vertex_ptr vA1 = cA1->vertex();
    vertex_ptr vB0 = cB0->vertex();
    vertex_ptr vB1 = cB1->vertex();
    face_vertex_ptr cA0n = cA1->next();
    face_vertex_ptr cA1n = cA0->next();
    /*
    vA0->print();
    vA1->print();
    vB0->print();
    vB1->print();

    std::cout << " Av: " << vA0->position_in_set() << " "
              << vA1->position_in_set() << std::endl;
    std::cout << " Af: " << cA0->face()->position_in_set() << " "
              << cA1->face()->position_in_set() << std::endl;
    std::cout << " Bv: " << vB0->position_in_set() << " "
              << vB1->position_in_set() << std::endl;
    std::cout << " Bf: " << cB0->face()->position_in_set() << " "
              << cB1->face()->position_in_set() << std::endl;
    std::cout << " cc: " << cA0n << " " << cA1n << std::endl;
    */
    eA->set(cA0, cB1);
    eB->set(cB0, cA1);
    vA0->set_front(cA0);
    vA1->set_front(cA1);

    if (vA0 == vB0) {
      // std::cout << " A.0 " << std::endl;
      vertex_ptr vN0 = new vertex_type();
      obj_in->push_vertex(vN0);

      coordinate_type cN = ci::get_coordinate<SPACE>(vA0);
      ci::set_coordinate<SPACE>(cN, vN0);

      vN0->set_front(cA0n); // this sets vertex
      vA0->update();
      vN0->update();

      // vA0->print();
      // vN0->print();
    } else {
      // std::cout << " A.1 " << std::endl;
      vA0->update();
      // vA0->print();
      obj_in->remove_vertex(vB0->position_in_set());
    }

    if (vA1 == vB1) {
      // std::cout << " B.0 " << std::endl;
      vertex_ptr vN1 = new vertex_type();
      obj_in->push_vertex(vN1);

      coordinate_type cN = ci::get_coordinate<SPACE>(vA1);
      ci::set_coordinate<SPACE>(cN, vN1);

      vN1->set_front(cA1n); // this sets vertex
      vA1->update();
      vN1->update();
      // vA1->print();
      // vN1->print();

    } else {
      // std::cout << " B.1 " << std::endl;
      vA1->update();
      // vA1->print();
      obj_in->remove_vertex(vB1->position_in_set());
    }

    vA0->template set<real>(SPACE::vertex_index::SMOOTH, 1.0);
    vA1->template set<real>(SPACE::vertex_index::SMOOTH, 1.0);

    // this->_surf->verify();
  }

  bool compatible(edge_ptr ea, edge_ptr eb) {
    bool ahasb1 = ea->v1()->vertex()->has_vertex(eb->v1()->vertex());
    bool ahasb2 = ea->v2()->vertex()->has_vertex(eb->v2()->vertex());
    return (!ahasb1 && !ahasb2);
  }

  void alignEdges(edge_ptr eA, edge_ptr eB) {

    if (eA->v1()->vertex() == eB->v1()->vertex()) {
      return;
    }
    if (eA->v2()->vertex() == eB->v2()->vertex()) {
      return;
    }

    if (eA->v1()->vertex() == eB->v2()->vertex()) {
      eB->swap_corners();
      return;
    }

    if (eA->v2()->vertex() == eB->v1()->vertex()) {
      eB->swap_corners();
      return;
    }
    coordinate_type ca0 = ci::get_coordinate<SPACE>(eA->v1());
    coordinate_type ca1 = ci::get_coordinate<SPACE>(eA->v2());
    coordinate_type cb0 = ci::get_coordinate<SPACE>(eB->v1());
    coordinate_type cb1 = ci::get_coordinate<SPACE>(eB->v2());

    T d0 = 1.0 / 3.0 * ((cb0 - ca0).squaredNorm() + (cb1 - ca1).squaredNorm());
    T d1 = 1.0 / 3.0 * ((cb0 - ca1).squaredNorm() + (cb1 - ca0).squaredNorm());

    if (d1 < d0)
      eB->swap_corners();
    /*
    std::cout << eA->v1()->vertex()->position_in_set() << " "
              << eA->v2()->vertex()->position_in_set() << std::endl;
    std::cout << eB->v1()->vertex()->position_in_set() << " "
              << eB->v2()->vertex()->position_in_set() << std::endl;
    */
  }

  void averageVerts(edge_ptr eA, edge_ptr eB) {
    auto avg = [this](face_vertex_ptr fvA, face_vertex_ptr fvB) {
      vertex_ptr vA = fvA->vertex();
      vertex_ptr vB = fvB->vertex();

      coordinate_type cA = ci::get_coordinate<SPACE>(vA);
      coordinate_type cB = ci::get_coordinate<SPACE>(vB);

      coordinate_type c = 0.5 * (cA + cB);

      ci::set_coordinate<SPACE>(c, vA);
      ci::set_coordinate<SPACE>(c, vB);
    };

    face_vertex_ptr fvA0 = eA->v1();
    face_vertex_ptr fvA1 = eA->v2();
    face_vertex_ptr fvB0 = eB->v1();
    face_vertex_ptr fvB1 = eB->v2();
    avg(eA->v1(), eB->v1());
    avg(eA->v2(), eB->v2());
  }
};

template <typename SPACE> class edge_collapse {
  M2_TYPEDEFS;

public:
  vertex_ptr base_collapse(surf_ptr obj, edge_ptr e, bool del2 = true) {
    // std::cout << "collapsing: " << e->position_in_set() << std::endl;
    face_vertex_ptr va0 = e->v1();
    face_vertex_ptr vb0 = e->v2();
    vertex_ptr v1 = va0->vertex();
    vertex_ptr v2 = vb0->vertex();

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
    // v1->print();
    // v2->print();
    /*
    std::cout << " sz: " << v1->size() << " " << v2->size() << std::endl;

    std::cout << " a: " << va0->vertex()->position_in_set() << " "
              << va1->vertex()->position_in_set() << " "
              << va2->vertex()->position_in_set() << std::endl;
    std::cout << " b: " << vb0->vertex()->position_in_set() << " "
              << vb1->vertex()->position_in_set() << " "
              << vb2->vertex()->position_in_set() << std::endl;
    */
    // std::cout << 1 << std::endl;
    ea1->set(va0p, va2p);
    eb1->set(vb0p, vb2p);

    v1->set_front(va0p);
    va0p->vertex()->set_front(va0p);
    va2p->vertex()->set_front(va2p);
    vb0p->vertex()->set_front(vb0p);
    vb2p->vertex()->set_front(vb2p);

    va0p->vertex()->update();
    va2p->vertex()->update();
    vb0p->vertex()->update();
    vb2p->vertex()->update();

    face_ptr fa = va0->face();
    face_ptr fb = vb0->face();

    v1->color.r = 1.0;

    obj->remove_face(fa->position_in_set());
    obj->remove_face(fb->position_in_set());

    obj->remove_edge(e->position_in_set());
    obj->remove_edge(ea2->position_in_set());
    obj->remove_edge(eb2->position_in_set());

    if (v1 != v2 && del2 == true) {
      obj->remove_vertex(v2->position_in_set());
    }

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
    v1->update_all();
    // obj->verify();
    return v1;
  }

  vertex_ptr degenerate_collapse(surf_ptr obj, edge_ptr e) {
    face_vertex_ptr fva = e->v1();
    face_vertex_ptr fvb = e->v2();
    vertex_ptr v1 = fva->vertex();
    vertex_ptr v2 = fvb->vertex();
    vertex_ptr v1p = fva->prev()->vertex();
    vertex_ptr v2p = fvb->prev()->vertex();

    construct<SPACE> cons;
    // std::cout << " collapse; " << v1->size() << " " << v2->size() <<
    // std::endl; std::cout << " collapse; " << v1->position_in_set() << " "
    //           << v2->position_in_set() << std::endl;

    if (v1->is_degenerate()) {
      // std::cout << " deg collapse 1: " << v1->position_in_set() << " "
      //           << v1->size() << std::endl;
      return NULL;
      // cons.delete_vertex_primitive(obj, v1);
    }
    if (v2->is_degenerate()) {
      // std::cout << " deg collapse 2: " << v2->position_in_set() << " "
      //           << v2->size() << std::endl;
      return NULL;
      cons.delete_vertex_primitive(obj, v2);
    }

    std::vector<edge_ptr> edges = v1->get_shared_edges(v2);
    if (edges.size() == 1) { // not a pinch point
      // std::cout << " edge collapse 1: " << v2->position_in_set() << " "
      //           << v2->size() << std::endl;

      return base_collapse(obj, e, true);
    } else if (edges.size() == 2) {
      // std::cout << " edge collapse 2: " << edges[0]->position_in_set() << " "
      //           << edges[1]->position_in_set() << std::endl;
      edge_merge<SPACE> merge;
      merge(obj, edges[0], edges[1]);
      return NULL;
    } else {
      return NULL;
    }
    /*
    else if (edges.size() > 1) { // pinch point collapse both eges
      std::cout << " pinch collapse! e: " << e->position_in_set()
                << " , " << edges.size()
                << " v1: " << v1->position_in_set()
                << " v2: " << v2->position_in_set() << std::endl;
      std::cout << 3 << std::endl;
      // make sure edges point to opposite vertices
      int i = 0;
      for (auto ei : edges) {
        if (ei->v1()->vertex() == v2 && i == 0){
          ei->swap_corners();
          std::cout << 4 << std::endl;
        }
        else if (ei->v1()->vertex() == v1 && i != 0){
          ei->swap_corners();
          std::cout << 5 << std::endl;
        }
        // do the collapse
        std::cout << i << " " << 6 << std::endl;
        base_collapse(obj, ei, false);
        i++;
      }
    }
    */
  }

  vertex_ptr operator()(surf_ptr obj, edge_ptr e) {
    return degenerate_collapse(obj, e);
  }
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
    /*
    std::cout << "===============" << std::endl;
    v1->print();
    std::cout << "before flip: v1" << std::endl;
    v1->print();
    std::cout << "before flip: v2" << std::endl;
    v2->print();
    std::cout << "before flip: v3" << std::endl;
    v3->print();
    std::cout << "before flip: v4" << std::endl;
    v4->print();
    */
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

    v1->set_front(fv21);
    v2->set_front(fv10);
    v3->set_front(fv11);
    v4->set_front(fv20);

    edge_ptr e11 = fv11->edge();
    edge_ptr e12 = fv12->edge();
    edge_ptr e21 = fv21->edge();
    edge_ptr e22 = fv22->edge();

    e11->set_this(fv11, fv22);
    e21->set_this(fv21, fv12);
    e12->set_this(fv12, fv11);
    e22->set_this(fv22, fv21);

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
    v1->update_all();
    v2->update_all();
    v3->update_all();
    v4->update_all();
    /*
    std::cout << "after flip: v1" << std::endl;
    v1->print();
    std::cout << "after flip: v2" << std::endl;
    v2->print();
    std::cout << "after flip: v3" << std::endl;
    v3->print();
    std::cout << "after flip: v4" << std::endl;
    v4->print();
    */
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
      vertices[i]->topologyChangeId = -1;
    }
  }

  virtual edge_array get_edges() { return edge_array(); }
  virtual bool op(edge_array edges) { return true; }
  virtual bool op() { return true; }

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

  bool skip_edge(edge_ptr e) {
    if (e->is_degenerate())
      return true;

    if (e->v1()->vertex()->is_degenerate())
      return true;

    if (e->v2()->vertex()->is_degenerate())
      return true;

    if (e->v1()->prev()->vertex() == e->v2()->prev()->vertex())
      return true;

    if (e->v1()->face()->size() != 3)
      return true;
    if (e->v2()->face()->size() != 3)
      return true;

    return false;
  }

  virtual bool op() {
    // TIMER function//TIMER(__FUNCTION__);
    edge_array &edges = this->_surf->get_edges();

    bool flipped = false;
    edge_array permEdges;
    for (int i = 0; i < edges.size(); i++) {
      if (!this->_surf->has_edge(i))
        continue;
      edge_ptr e = edges[i];
      if (skip_edge(e))
        continue;
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

      if (skip_edge(e))
        continue;

      bool pinned = e->v1()->vertex()->pinned == true &&
                    e->v2()->vertex()->pinned == true;
      bool notPinned = e->v1()->vertex()->pinned != true &&
                       e->v2()->vertex()->pinned != true;
      // std::cout << "B" << std::endl;
      // if (pinned)
      //  continue;
      // if (notPinned)
      //  continue;

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
      T A0 = va::calculate_area<T>(c0, c1, c2);
      T A1 = va::calculate_area<T>(c2, c3, c0);
      edge_flip<SPACE> flip;
      if (A0 / A1 > 4.0) {
        flip(e);
        flipped = true;
        continue;
      } else if (A1 / A0 > 4) {
        flip(e);
        flipped = true;
        continue;
      }
      // std::cout << A0 << "  " << A1 << " " << A0 / A1 << " " << A1 / A0
      //           << std::endl;
      //  interior angles
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
                                         // surface angles

      // current normals
      coordinate_type N00 = va::calculate_normal(c1, c0, c2);
      coordinate_type N01 = va::calculate_normal(c3, c2, c0);
      // new normals
      coordinate_type N10 = va::calculate_normal(c0, c3, c1);
      coordinate_type N11 = va::calculate_normal(c2, c1, c3);
      T cosN0 = va::dot(N00, N01);
      T tSame = M_PI - acos(cosN0);
      T cosN1 = va::dot(N10, N11);
      T tFlip = M_PI - acos(cosN1);
      T nFlip = tFlip;

      T eFlip = cFlip * cFlip + 0.5 * tFlip * tFlip;
      T eSame = cSame * cSame + 0.0 * tSame * tSame;

      // if (false) {
      if (eFlip < 0.5 * eSame) {

        flip(e);
        flipped = true;
      }
    }
    return flipped;
  }

}; // set_operations
#if 1
template <typename SPACE>
class edge_merger : public triangle_operations_base<SPACE> {
public:
  M2_TYPEDEFS;

  edge_merger(const surf_ptr &surf, real max)
      : triangle_operations_base<SPACE>(surf, max) {}

  void clipEar(edge_ptr e, surf_ptr mesh) {

    face_ptr fA = e->v1()->face(), fB = e->v2()->face();
    // std::cout << "shared_edge: " << fA->position_in_set() << " "
    //          << fB->position_in_set() << std::endl;
    face_vertex_ptr fvA0 = e->v1();
    face_vertex_ptr fvA1 = fvA0->next();
    face_vertex_ptr fvA2 = fvA1->next();

    face_vertex_ptr fvB1 = e->v2();
    face_vertex_ptr fvB2 = fvB1->prev();
    face_vertex_ptr fvB0 = fvB2->prev();

    vertex_ptr vA0 = fvA0->vertex();
    vertex_ptr vA1 = fvA1->vertex();
    vertex_ptr vA2 = fvA2->vertex();

    vertex_ptr vB0 = fvB0->vertex();
    vertex_ptr vB1 = fvB1->vertex();
    vertex_ptr vB2 = fvB2->vertex();

    if (vA0 == vA1)
      return;
    if (vA0 == vA2)
      return;
    if (vA1 == vA2)
      return;

    if (vB0 == vB1)
      return;
    if (vB0 == vB2)
      return;
    if (vB1 == vB2)
      return;

    fvA0->vertex()->remove_face_vertex(fvA0);
    fvA1->vertex()->remove_face_vertex(fvA1);
    fvA2->vertex()->remove_face_vertex(fvA2);

    fvB0->vertex()->remove_face_vertex(fvB0);
    fvB1->vertex()->remove_face_vertex(fvB1);
    fvB2->vertex()->remove_face_vertex(fvB2);
    /*
    std::cout << " va0: " << vA0->position_in_set() << " "
              << fvA2->coedge()->vertex()->position_in_set() << std::endl;
    std::cout << " va1: " << vA1->position_in_set() << " "
              << fvA0->coedge()->vertex()->position_in_set() << std::endl;
    std::cout << " va2: " << vA2->position_in_set() << " "
              << fvA1->coedge()->vertex()->position_in_set() << std::endl;
    */
    edge_ptr eA1 = fvA1->edge();
    edge_ptr eA2 = fvA2->edge();

    edge_ptr eB0 = fvB0->edge();
    edge_ptr eB1 = fvB1->edge();
    edge_ptr eB2 = fvB2->edge();

    if (vA2 != vB2) {
      std::cout << " not adding " << std::endl;
      mesh->remove_vertex(vB2->position_in_set());
    } else {
      std::cout << " adding new vertex" << std::endl;
      vertex_ptr vB2p = new vertex_type(*vA2);
      vB2p->front() = fvB0->coedge();
      vB2p->update_all();
      mesh->push_vertex(vB2p);
    }
  }

  vector<line_pair> get_edge_pairs_to_merge(std::vector<edge_ptr> &edges) {

    std::vector<edge_ptr> filtered_edges;
    for (int i = 0; i < edges.size(); i++) {
      if (!this->_surf->has_edge(i))
        continue;
      if (edges[i]->is_degenerate())
        continue;
      filtered_edges.push_back(edges[i]);
    }

    std::vector<line_type> lines = asawa::ci::get_lines<SPACE>(filtered_edges);
    bins::aabb_tree<SPACE, line_type> tree(lines);

    // std::cout << __FUNCTION__ << ": edges.size(): " << edges.size()
    //           << std::endl;
    // std::cout << __FUNCTION__ << ": tree.nodes.size(): " << tree.nodes.size()
    //           << std::endl;

    auto lineDist = [](const line_type &A, const line_type &B) -> T {
      return A.dist(B);
    };

    T tol = this->_thresh;

    vector<line_pair> collected(lines.size());
// for (int i = 0; i < 1; i++) {
// std::cout << " foo !" << std::endl;
#pragma omp parallel for
    for (int i = 0; i < lines.size(); i++) {
      line_type lineNear = bins::getNearest<SPACE, line_type, line_type>(
          lines[i], tree, lines, lineDist, tol);
      collected[i] = line_pair(lines[i], lineNear);
      // if (lineNear.edgeId > -1) {
      //   collected.push_back(line_pair(lines[i], lineNear));
      // }
    }

    std::sort(collected.begin(), collected.end(), std::less<line_pair>());

    // std::cout << " full list: " << collected.size() << std::endl;
    auto it = std::unique(collected.begin(), collected.end());

    collected.erase(it, collected.end());
    // std::cout << " unique list: " << collected.size() << std::endl;

    collected.erase(std::remove_if(collected.begin(), collected.end(),
                                   [edges, tol](auto p) {
                                     if (p.A.edgeId < 0)
                                       return true;
                                     if (p.B.edgeId < 0)
                                       return true;
                                     edge_ptr eA = edges[p.A.edgeId];
                                     edge_ptr eB = edges[p.B.edgeId];
                                     T d = p.dist();

                                     if (d > tol)
                                       return true;

                                     if (eA->flags[0] > 0)
                                       return true;
                                     if (eB->flags[0] > 0)
                                       return true;

                                     eA->flags[0] = 1;
                                     eB->flags[0] = 1;
                                     return false;
                                   }),
                    collected.end());

    // std::cout << " merging " << collected.size() << " lines!" << std::endl;
    return collected;
  }

  bool merge_collected(std::vector<line_pair> &collected,
                       std::vector<edge_ptr> &edges) {

    bool merging = true;
    int phase = 0;
    while (merging && collected.size()) {

      // collected = collapse_bridges(collected, faces);

      std::sort(collected.begin(), collected.end(), std::less<line_pair>());

      // std::cout << "merging phase" << phase
      //           << ", line pair count: " << collected.size() << std::endl;

      std::vector<line_pair> new_collect;

      for (auto p : collected) {
        if (!this->_surf->has_edge(p.A.edgeId) ||
            !this->_surf->has_edge(p.B.edgeId)) {
          // std::cout << " moving on... " << std::endl;
          continue;
        };

        edge_ptr eA = edges[p.A.edgeId];
        edge_ptr eB = edges[p.B.edgeId];

        if (eA->is_degenerate() || eB->is_degenerate()) {
          continue;
        };
        // std::cout << " shared: " << shared << std::end
        coordinate_type NA = ci::normal<SPACE>(eA);
        coordinate_type NB = ci::normal<SPACE>(eB);
        real angle = va::dot(NA, NB);
        T d = p.dist();

        if (angle < -0.65) {

          edge_merge<SPACE> merge;
          merge.averageVerts(eA, eB);
          merge(this->_surf, eA, eB);
        }
        merging = true;
      }
      collected = new_collect;
    }

    this->_surf->update_all();
    this->_surf->verify();
    this->reset_flags();
    return merging;
  }

  virtual bool op() {
    std::vector<edge_ptr> edges = this->_surf->get_edges();
    vector<line_pair> collected = get_edge_pairs_to_merge(edges);
    if (collected.empty())
      return false;
    return this->merge_collected(collected, edges);
    return true;
  }

  void reset_flags() {

    for (auto e : this->_surf->get_edges()) {
      if (!e)
        continue;
      e->flags[0] = 0;
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

    asawa::construct<SPACE> cons;
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
    edge_ptr e0, e1;
    vertex_ptr nv = split(this->_surf, e, e0, e1);
    nv->topologyChangeId = i;

    // nv->coordinate() = c;
    fv1->edge()->flags[0] = 1;
    fv2->edge()->flags[0] = 1;

    edge_ptr e11 = subdivideFace(fv1->next());
    edge_ptr e21 = subdivideFace(fv2->next());

    e0->topologyChangeId = 4 * i + 0;
    e1->topologyChangeId = 4 * i + 1;
    e11->topologyChangeId = 4 * i + 2;
    e21->topologyChangeId = 4 * i + 3;

    e11->flags[1] = 0;
    e21->flags[1] = 0;
  }

  virtual edge_array get_edges() {
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

    // std::cout << " - sorting " << edgesToSplit.size() << " edges" <<
    // std::endl;
    std::sort(edgesToSplit.begin(), edgesToSplit.end(),
              this->mEdgeSorterGreater);
    return edgesToSplit;
  }

  virtual bool op(edge_array edgesToSplit) {
    // TIMER function//TIMER(__FUNCTION__);

    // std::cout << "splitting " << edgesToSplit.size() << " edges" <<
    // std::endl;

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

  virtual edge_array get_edges() {
    edge_array &edges = this->_surf->get_edges();

    edge_array edgeToCollapse;
#if 1
    for (int i = 0; i < edges.size(); i++) {
      if (!this->_surf->has_edge(i))
        continue;
      edge_ptr ei = edges[i];

      face_vertex_ptr fv1 = ei->v1();
      face_vertex_ptr fv2 = ei->v2();

      if (fv1->vertex()->size() < 3 || fv2->vertex()->size() < 3) {
        continue;
      }
      if (ei->flags[0])
        continue;

      if (fv1->prev()->vertex()->size() < 3 ||
          fv2->prev()->vertex()->size() < 3) {
        continue;
      }

      bool pinned = ei->v1()->vertex()->pinned == true &&
                    ei->v2()->vertex()->pinned == true;
      if (pinned)
        continue;

      T l = this->length(ei);

      if (l < this->_thresh) {
        edgeToCollapse.push_back(ei);
        fv1->next()->edge()->flags[0] = true;
        fv2->next()->edge()->flags[0] = true;
        continue;
      }
    }
#endif

    // std::cout << " - sorting " << edgeToCollapse.size() << " edges"
    //           << std::endl;
    std::sort(edgeToCollapse.begin(), edgeToCollapse.end(),
              this->mEdgeSorterLesser);
    /*
    edgeToCollapse.erase(std::remove_if(edgeToCollapse.begin(),
                                        edgeToCollapse.end(),
                                        [](edge_ptr e) { return e->flags[0]; }),
                         edgeToCollapse.end());
    */
    // std::cout << " collapsing: " << edgeToCollapse.size() << std::endl;
    return edgeToCollapse;
  }

  virtual bool op(edge_array edgeToCollapse) {
    // TIMER function//TIMER(__FUNCTION__);
    std::cout << " - deleting: " << edgeToCollapse.size() << " edges"
              << std::endl;
    std::vector<size_t> indicesToCollapse;
    for (int i = 0; i < edgeToCollapse.size(); i++) {
      indicesToCollapse.push_back(edgeToCollapse[i]->position_in_set());
    }

    edge_array edges = this->_surf->get_edges();

    for (int i = 0; i < indicesToCollapse.size(); i++) {
      if (!this->_surf->has_edge(indicesToCollapse[i]))
        continue;
      edge_ptr e = edges[indicesToCollapse[i]];
      edge_collapse<SPACE> collapser;
      vertex_ptr nv = collapser(this->_surf, e);
      if (nv)
        nv->topologyChangeId = i;
    }
    return edgeToCollapse.size() > 0;
  }
};

#endif
#if 1
template <typename SPACE>
class degenerate_deleter : public triangle_operations_base<SPACE> {
public:
  M2_TYPEDEFS;

  degenerate_deleter(const surf_ptr &surf)
      : triangle_operations_base<SPACE>(surf, 0.0) {}

  bool op(real eps) {
    asawa::construct<SPACE> cons;

    bool testing = true;
    int i;

    while (testing) {
      i = 0;
      bool deleting = false;
      for (auto e : this->_surf->get_edges()) {
        if (!this->_surf->has_edge(i++))
          continue;

        deleting |= cons.delete_degenerates(this->_surf, e);
      }
      testing = deleting;

      i = 0;
      deleting = false;
      for (auto v : this->_surf->get_vertices()) {
        if (!this->_surf->has_vertex(i++))
          continue;
        deleting = cons.delete_degenerates(this->_surf, v);
      }

      i = 0;
      deleting = false;
      asawa::edge_collapse<SPACE> collapser;
      for (auto f : this->_surf->get_faces()) {
        if (!this->_surf->has_face(i++))
          continue;
        deleting = cons.delete_degenerates(this->_surf, f, eps);
        if (f->is_degenerate()) {
          if (f->get_front()->edge()) {
            collapser.degenerate_collapse(this->_surf, f->get_front()->edge());
            testing = true;
          }
        }

        real area = asawa::ci::area<SPACE>(f);
        if (area < eps) {
          if (f->get_front()->edge()) {
            collapser.degenerate_collapse(this->_surf, f->get_front()->edge());
            testing = true;
          }
        }

        return false;
      }

      testing |= deleting;
    }

    i = 0;
    for (auto f : this->_surf->get_faces()) {
      if (!this->_surf->has_face(i++))
        continue;
      face_vertex_ptr fv = f->get_front();
      if (f->is_null() && fv->vertex()->size() > 1) {
        this->_surf->remove_face(f->position_in_set());
        delete (fv);
      }
    }

    asawa::remesh<SPACE> rem;
    rem.triangulate(this->_surf);
    return true;
  }
};
#endif

#if 1
enum op_type { split, merge };

template <typename SPACE> class vertex_policy {
public:
  M2_TYPEDEFS;

  vertex_policy(typename SPACE::vertex_index id) : _id(id) {}
  virtual void calc(int i, edge_ptr &e, op_type op) = 0;
  virtual void apply(int i, vertex_ptr &v) = 0;
  virtual void reserve(long N) = 0;
  virtual void clear() = 0;
  typename SPACE::vertex_index _id;
};

template <typename SPACE, typename TYPE>
class vertex_policy_t : public vertex_policy<SPACE> {
public:
  M2_TYPEDEFS;

  vertex_policy_t(typename SPACE::vertex_index id) : vertex_policy<SPACE>(id) {
    std::cout << " id: " << int(this->_id) << std::endl;
  }

  virtual void calc(int i, edge_ptr &e, op_type op) {

    face_vertex_ptr v0 = e->v1();
    face_vertex_ptr v1 = v0->prev();
    face_vertex_ptr v2 = e->v2();
    face_vertex_ptr v3 = v2->prev();

    TYPE t0 = v0->vertex()->template get<TYPE>(this->_id);
    TYPE t1 = v1->vertex()->template get<TYPE>(this->_id);
    TYPE t2 = v2->vertex()->template get<TYPE>(this->_id);
    TYPE t3 = v3->vertex()->template get<TYPE>(this->_id);
    _vals[i] = 0.5 * (t0 + t2);
    return;

    coordinate_type c0 = ci::get_coordinate<SPACE>(v0);
    coordinate_type c1 = ci::get_coordinate<SPACE>(v1);
    coordinate_type c2 = ci::get_coordinate<SPACE>(v2);
    coordinate_type c3 = ci::get_coordinate<SPACE>(v3);

    T cos1 = va::abs_cos<T>(c1, c0, c2);
    T cos3 = va::abs_cos<T>(c3, c2, c0);
    T sin1 = sqrt(1.0 - cos1 * cos1);
    T sin3 = sqrt(1.0 - cos3 * cos3);
    T hcot1 = sin1 / (1.0 - cos1);
    T hcot3 = sin3 / (1.0 - cos3);

    real w0 = hcot1 + hcot3;
    real w2 = hcot1 + hcot3;

    real w1 = va::abs_cotan(c0, c1, c2) + va::abs_cotan(c2, c1, c0);
    real w3 = va::abs_cotan(c0, c3, c2) + va::abs_cotan(c2, c3, c0);

    _vals[i] = (w0 * t0 + w1 * t1 + w2 * t2 + w3 * t3) / (w0 + w1 + w2 + w3);
  }

  virtual void apply(int i, vertex_ptr &v) {
    v->template set<TYPE>(this->_id, _vals[i]);
  }

  virtual void reserve(long N) { _vals = std::vector<TYPE>(N); }
  virtual void clear() { _vals.clear(); }
  std::vector<TYPE> _vals;
};

template <typename SPACE> class edge_policy {
public:
  M2_TYPEDEFS;

  edge_policy(typename SPACE::edge_index id) : _id(id) {}
  virtual void calc(int i, edge_ptr &e, op_type op) = 0;
  virtual void apply(int i, edge_ptr &e) = 0;
  virtual void reserve(long N) = 0;
  virtual void clear() = 0;
  typename SPACE::edge_index _id;
};

template <typename SPACE, typename TYPE>
class edge_policy_t : public edge_policy<SPACE> {
public:
  M2_TYPEDEFS;

  edge_policy_t(typename SPACE::edge_index id) : edge_policy<SPACE>(id) {}

  virtual void calc(int i, edge_ptr &e, op_type op) {
    //_vals[4 * i + 0] = blah;
    //_vals[4 * i + 1] = blah;
    //_vals[4 * i + 2] = blah;
    //_vals[4 * i + 3] = blah;
    std::cout << _vals.size() << " " << 4 * i << std::endl;
  }

  virtual void apply(int i, edge_ptr &e) {

    e->template set<TYPE>(this->_id, _vals[i]);
  }
  virtual void reserve(long N) { _vals = std::vector<TYPE>(4 * N); }
  virtual void clear() { _vals.clear(); }

  std::vector<TYPE> _vals;
};

template <typename SPACE> class surf_integrator {
public:
  M2_TYPEDEFS;

  surf_integrator(const surf_ptr &surf, real min = 0.5, real max = 3.0,
                  real merge = 1.0)
      : _surf(surf) {

    this->add_default_vertex_policy<coordinate_type>(
        SPACE::vertex_index::COORDINATE);

    double l0 = asawa::ci::geometric_mean_length<SPACE>(_surf);
    std::cout << " integrator avg length: " << l0 << std::endl;
    _min = min * l0;
    _max = max * _min;
    _merge = merge * _min;
  }

  void set_mesh(const surf_ptr &surf) { this->_surf = surf; }

  void add_vertex_policy(vertex_policy<SPACE> *policy) {
    _vertex_policies.push_back(policy);
  }

  template <typename TYPE>
  void add_default_vertex_policy(typename SPACE::vertex_index id) {
    _vertex_policies.push_back(new vertex_policy_t<SPACE, TYPE>(id));
  }

  void add_edge_policy(edge_policy<SPACE> *policy) {
    _edge_policies.push_back(policy);
  }

  void cacheBary(asawa::surf<SPACE> *surf) {
    M2_TYPEDEFS;
    int i = 0;
    for (auto f : surf->get_faces()) {

      if (!surf->has_face(i++))
        continue;

      typename SPACE::coordinate_type c = asawa::ci::center<SPACE>(f);
      typename SPACE::coordinate_type l = asawa::ci::point_to_bary<SPACE>(f, c);

      assert(!isnan(l[0]));
      assert(!isnan(l[1]));
      assert(!isnan(l[2]));

      int i = 0;
      asawa::for_each_face<SPACE>(f, [l, i](face_vertex_ptr fv) mutable {
        typename SPACE::real li = 0.0;
        if (i < 3)
          li = l[i];
        fv->template set<typename SPACE::real>(SPACE::face_vertex_index::BARY,
                                               li);
        i++;
      });
    }
  }

  void smoothMesh(typename SPACE::real C, int N) {

    // return;
    vertex_array &vertices = this->_surf->get_vertices();
    int i = 0;
    // asawa::area_laplacian_0<SPACE, coordinate_type> M(this->_surf);
    asawa::laplacian3<SPACE> M(this->_surf);

    std::cout << "ugly smoothing " << std::flush;
    coordinate_array coords = asawa::ci::get_coordinates<SPACE>(this->_surf);
    coordinate_array normals =
        asawa::ci::get_vertex_normals<SPACE>(this->_surf);

    std::vector<real> sm =
        ci::get<SPACE, real>(_surf, SPACE::vertex_index::SMOOTH);
    for (int k = 0; k < N; k++) {
      std::cout << "." << std::flush;
      M.build();
      coords = M.smooth(coords, C, C + 3e-5);
    }

    asawa::ci::set_coordinates<SPACE>(coords, this->_surf);
    std::cout << "done!" << std::endl;
  }

  bool operate_edges(triangle_operations_base<SPACE> &op, op_type opt) {

    op.reset_flags();
    edge_array edges_to_op = op.get_edges();

    if (edges_to_op.empty())
      return false;

    auto calc_policy = [&opt](const edge_array &edges_to_op, auto &&policies) {
      for (auto p : policies) {
        p->reserve(edges_to_op.size());
      }
      int i = 0;
      for (auto e : edges_to_op) {
        for (auto p : policies) {
          p->calc(i, e, opt);
        }
        i++;
      }
    };

    calc_policy(edges_to_op, _vertex_policies);
    calc_policy(edges_to_op, _edge_policies);

    op.op(edges_to_op);

    // std::cout << "assigning calc'd vertex values" << std::endl;
    auto apply_policy = [](auto &&primitives, auto &&policies) {
      for (auto prim : primitives) {
        if (!prim)
          continue;
        int topologyChangeId = prim->topologyChangeId;
        if (topologyChangeId < 0)
          continue;

        for (auto p : policies) {
          p->apply(topologyChangeId, prim);
        }
        prim->topologyChangeId = -1;
      }

      for (auto p : policies) {
        p->clear();
      }
    };

    vertex_array vertices = _surf->get_vertices();
    apply_policy(vertices, _vertex_policies);
    edge_array edges = _surf->get_edges();
    apply_policy(edges, _edge_policies);

    return true;
  }

  bool splitEdges() {
    asawa::edge_splitter<SPACE> splitter(_surf, _max);
    return this->operate_edges(splitter, split);
  }

  bool collapsEdges() {
    using namespace asawa;
    M2_TYPEDEFS;
    // std::cout << "-- collapsing --" << std::endl;
    asawa::edge_collapser<SPACE> collapser(_surf, _min);
    return this->operate_edges(collapser, merge);
  }

  bool mergeEdges() {
    // std::cout << "-- merging --" << std::endl;
    edge_merger<SPACE> merger(this->_surf, _merge);
    return merger.op();
  }

  void flipEdges() {
    edge_flipper<SPACE> flipper(this->_surf);
    flipper.op();
  }

  void delete_degenerates() {

    // std::cout << "-- delete_degenerates --" << std::endl;
    degenerate_deleter<SPACE> zapper(this->_surf);
    zapper.op(0.25 * _min * _min);
    // std::cout << " done! " << std::endl;
  }

  bool integrate() { // mergeTris(_meshGraph, rxA, rxB);
    bool relaxing = true;
    int i = 0;
    while (relaxing) {
      if (i++ > 5)
        break;
      relaxing = splitEdges();
      delete_degenerates();
      relaxing |= collapsEdges();
      delete_degenerates();
      relaxing |= mergeEdges();
      delete_degenerates();
      flipEdges();
      delete_degenerates();

      this->_surf->update_all();
      this->_surf->pack();
      // relaxing = false;
    }

    // diffuse edge flags
    this->cacheBary(this->_surf);

    smoothMesh(0.001, 10);

    std::vector<real> x =
        ci::get<SPACE, real>(_surf, SPACE::vertex_index::SMOOTH);

    asawa::diffuse<SPACE> diff(_surf);
    std::vector<real> fs(x);

    for (auto &xi : x) {
      xi = std::min(std::max(xi, 0.0), 1.0);
      xi *= 0.5;
    }
    /*
    for (auto &f : fs)
      f = 0.0;

    x = diff.second_order(x, fs, 0.01, _max);
    double Sm = 0.0;

    for (auto &xi : x) {
      Sm += xi;
      xi = clamp(xi, 0.0, 1.0);
    }

    std::cout << "   smooth norm: " << Sm / double(x.size()) << std::endl;
    */
    ci::set<SPACE, real>(_surf, x, SPACE::vertex_index::SMOOTH);

    return true;
  }
  asawa::surf<SPACE> *_surf;
  real _min = 0.5, _max = 3.0, _merge = 1.0;

  std::vector<vertex_policy<SPACE> *> _vertex_policies;
  std::vector<edge_policy<SPACE> *> _edge_policies;
};

#endif
}; // namespace asawa
#endif
