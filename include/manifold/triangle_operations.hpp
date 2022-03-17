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
#include "manifold/construct.hpp"
#include "manifold/laplacian.hpp"
#include "manifold/m2.hpp"
#include "manifold/remesh.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <set>
#include <vector>

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

template <typename SPACE> class edge_collapse {
  M2_TYPEDEFS;

public:
  vertex_ptr base_collapse(surf_ptr obj, edge_ptr e) {

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
    /*
    std::cout << " sz: " << v1->size() << " " << v2->size() << std::endl;

    std::cout << " a: " << va0->vertex()->position_in_set() << " "
              << va1->vertex()->position_in_set() << " "
              << va2->vertex()->position_in_set() << std::endl;
    std::cout << " b: " << vb0->vertex()->position_in_set() << " "
              << vb1->vertex()->position_in_set() << " "
              << vb2->vertex()->position_in_set() << std::endl;
    */
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

    if (v1 != v2) {
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
    vertex_ptr v1 = fva->prev()->vertex();
    vertex_ptr v2 = fvb->prev()->vertex();

    construct<SPACE> cons;
    // std::cout << " collapse; " << v1->size() << " " << v2->size() <<
    // std::endl;
    if (v1->is_degenerate()) {
      // std::cout << " deg collapse 1" << v1->position_in_set() << std::endl;
      cons.delete_vertex_primitive(obj, v1);
    }
    if (v2->is_degenerate()) {
      // std::cout << " deg collapse 2" << v2->position_in_set() << std::endl;
      cons.delete_vertex_primitive(obj, v2);
    }

    return base_collapse(obj, e);
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
    std::cout << "after flip: v1" << std::endl;
    v1->print();
    std::cout << "after flip: v2" << std::endl;
    v2->print();
    std::cout << "after flip: v3" << std::endl;
    v3->print();
    std::cout << "after flip: v4" << std::endl;
    v4->print();

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

    std::cout << "done:" << std::endl;
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
    ;
    if (e->v2()->vertex()->is_degenerate())
      return true;

    if (e->v1()->prev()->vertex() == e->v2()->prev()->vertex())
      return true;

    if (e->v1()->face()->size() != 3)
      return true;
    if (e->v2()->face()->size() != 3)
      return true;
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
      if (eFlip < eSame) {
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
class edge_merger : public triangle_operations_base<SPACE> {
public:
  M2_TYPEDEFS;

  edge_merger(const surf_ptr &surf, real max)
      : triangle_operations_base<SPACE>(surf, max) {}

  void joinEdges(edge_ptr eA, edge_ptr eB, surf_ptr mesh) {
    std::cout << " joining: " << eA->position_in_set() << " "
              << eB->position_in_set() << std::endl;

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

    eA->set(cA0, cB1);
    eB->set(cB0, cA1);
    vA0->set_front(cA0);
    vA1->set_front(cA1);

    if (vA0 == vB0) {
      std::cout << " A.0 " << std::endl;
      vertex_ptr vN0 = new vertex_type();
      this->_surf->push_vertex(vN0);
      coordinate_type cN = this->coordinate(vA0);
      this->coordinate(cN, vN0);
      vN0->set_front(cA0n); // this sets vertex
      vA0->update();
      vN0->update();

      vA0->print();
      vN0->print();
    } else {
      std::cout << " A.1 " << std::endl;
      vA0->update();
      vA0->print();
      mesh->remove_vertex(vB0->position_in_set());
    }

    if (vA1 == vB1) {
      std::cout << " B.0 " << std::endl;
      vertex_ptr vN1 = new vertex_type();
      this->_surf->push_vertex(vN1);
      coordinate_type cN = this->coordinate(vA1);
      this->coordinate(cN, vN1);
      vN1->set_front(cA1n); // this sets vertex
      vA1->update();
      vN1->update();
      vA1->print();
      vN1->print();

    } else {
      std::cout << " B.1 " << std::endl;
      vA1->update();
      vA1->print();
      mesh->remove_vertex(vB1->position_in_set());
    }

    // this->_surf->verify();
  }

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

    coordinate_type ca0 = this->coordinate(eA->v1());
    coordinate_type ca1 = this->coordinate(eA->v2());

    coordinate_type cb0 = this->coordinate(eB->v1());
    coordinate_type cb1 = this->coordinate(eB->v2());

    T d0 = 1.0 / 3.0 * ((cb0 - ca0).squaredNorm() + (cb1 - ca1).squaredNorm());
    T d1 = 1.0 / 3.0 * ((cb0 - ca1).squaredNorm() + (cb1 - ca0).squaredNorm());

    if (d1 < d0)
      eB->swap_corners();
    std::cout << eA->v1()->vertex()->position_in_set() << " "
              << eA->v2()->vertex()->position_in_set() << std::endl;
    std::cout << eB->v1()->vertex()->position_in_set() << " "
              << eB->v2()->vertex()->position_in_set() << std::endl;
  }

  void averageVerts(edge_ptr eA, edge_ptr eB) {
    auto avg = [this](face_vertex_ptr fvA, face_vertex_ptr fvB) {
      vertex_ptr vA = fvA->vertex();
      vertex_ptr vB = fvB->vertex();
      coordinate_type cA = this->coordinate(vA);
      coordinate_type cB = this->coordinate(vB);
      coordinate_type c = 0.5 * (cA + cB);
      this->coordinate(c, vA);
      this->coordinate(c, vB);
    };

    face_vertex_ptr fvA0 = eA->v1();
    face_vertex_ptr fvA1 = eA->v2();
    face_vertex_ptr fvB0 = eB->v1();
    face_vertex_ptr fvB1 = eB->v2();
    avg(eA->v1(), eB->v1());
    avg(eA->v2(), eB->v2());
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

    std::vector<line_type> lines = m2::ci::get_lines<SPACE>(filtered_edges);
    m2::aabb_tree<SPACE, line_type> tree(lines);

    std::cout << "tree.nodes.size(): " << tree.nodes.size() << std::endl;

    auto lineDist = [](const line_type &A, const line_type &B) -> T {
      return A.dist(B);
    };

    T tol = this->_thresh;

    vector<line_pair> collected;
    // for (int i = 0; i < 1; i++) {
    // std::cout << " foo !" << std::endl;
    for (int i = 0; i < lines.size(); i++) {
      line_type lineNear = m2::getNearest<SPACE, line_type, line_type>(
          lines[i], tree, lines, lineDist, tol);
      if (lineNear.edgeId > -1) {
        collected.push_back(line_pair(lines[i], lineNear));
      }
    }

    std::sort(collected.begin(), collected.end(), std::less<line_pair>());

    std::cout << " full list: " << collected.size() << std::endl;
    auto it = std::unique(collected.begin(), collected.end());

    collected.erase(it, collected.end());
    std::cout << " unique list: " << collected.size() << std::endl;

    collected.erase(std::remove_if(collected.begin(), collected.end(),
                                   [edges, tol](auto p) {
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

    std::cout << " merging " << collected.size() << " lines!" << std::endl;
    return collected;
  }

  bool compatible(edge_ptr ea, edge_ptr eb) {
    bool ahasb1 = ea->v1()->vertex()->has_vertex(eb->v1()->vertex());
    bool ahasb2 = ea->v2()->vertex()->has_vertex(eb->v2()->vertex());
    return (!ahasb1 && !ahasb2);
  }

  bool merge_collected(std::vector<line_pair> &collected,
                       std::vector<edge_ptr> &edges) {

    bool merging = true;
    int phase = 0;
    while (merging && collected.size()) {

      // collected = collapse_bridges(collected, faces);

      std::sort(collected.begin(), collected.end(), std::less<line_pair>());

      std::cout << "merging phase" << phase
                << ", line pair count: " << collected.size() << std::endl;

      std::vector<line_pair> new_collect;

      for (auto p : collected) {
        if (!this->_surf->has_edge(p.A.edgeId) ||
            !this->_surf->has_edge(p.B.edgeId)) {
          std::cout << " moving on... " << std::endl;
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

        alignEdges(eA, eB);
        if (!compatible(eA, eB))
          continue;

        if (angle < -0.75) {

          averageVerts(eA, eB);
          joinEdges(eA, eB, this->_surf);
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
    return this->merge_collected(collected, edges);
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
class face_merger : public triangle_operations_base<SPACE> {
public:
  M2_TYPEDEFS;

  face_merger(const surf_ptr &surf, real max)
      : triangle_operations_base<SPACE>(surf, max) {}

  void joinFaces(face_ptr fA, face_ptr fB, surf_ptr mesh) {

    face_vertex_ptr fvb0 = fB->fbegin();
    std::vector<edge_ptr> e_border;

    std::set<face_vertex_ptr> fv_all;
    std::set<edge_ptr> e_all;
    std::set<vertex_ptr> v_del;
    std::set<vertex_ptr> v_new;
    std::set<edge_ptr> f_del;

    face_vertex_ptr fvA0 = fA->fbegin();
    face_vertex_ptr fvB0 = fB->fbegin()->prev();

    // std::cout << "===========" << std::endl;
    // std::cout << "- fA -" << std::endl;

    // fA->print();
    m2::for_each_face<SPACE>(fA, [](face_vertex_ptr fv) {
      // fv->vertex()->print();
      fv->vertex()->print_adj_sz();
    });

    // std::cout << "===========" << std::endl;
    // std::cout << "- fB -" << std::endl;

    // fB->print();
    // m2::for_each_face<SPACE>(fB, [](face_vertex_ptr fv) {
    //   fv->vertex()->print();
    //   fv->vertex()->print_adj_sz();
    // });

    std::set<vertex_ptr> vd;
    m2::for_each_face<SPACE>(fA, [&vd](face_vertex_ptr fv) {
      if (fv->next()->vertex()->size() < 4) {
        vd.insert(fv->next()->vertex());
      }
    });

    m2::for_each_face<SPACE>(fB, [&vd](face_vertex_ptr fv) {
      if (fv->next()->vertex()->size() < 4) {
        vd.insert(fv->next()->vertex());
      }
    });

    bool ret = false;
    for (auto v : vd) {
      construct<SPACE> cons;
      // std::cout << " pre del: " << v->position_in_set() << v->size()
      //            << std::endl;
      for_each_vertex<SPACE>(v, [&fA, &fB, &ret](face_vertex_ptr fv) {
        if (fv->edge()) {
          ret |= fv->edge()->v1()->face() == fA;
          ret |= fv->edge()->v2()->face() == fA;
          ret |= fv->edge()->v1()->face() == fB;
          ret |= fv->edge()->v2()->face() == fB;
        }
      });
      cons.delete_vertex_primitive(mesh, v);
    }

    if (ret) {
      return;
    }

    int i = 0;
    m2::for_each_face<SPACE>(fA, [&](face_vertex_ptr fva0) {
      std::cout << i++ << std::endl;
      face_vertex_ptr fva1 = fva0->next();
      face_vertex_ptr fvb1 = fvb0->prev();
      face_vertex_ptr fva2 = fva0->prev();
      face_vertex_ptr fvb2 = fvb0->next();

      vertex_ptr vA0 = fva0->vertex();
      vertex_ptr vB0 = fvb0->vertex();
      vertex_ptr vA1 = fva1->vertex();

      vertex_ptr vB1 = fvb1->vertex();
      vertex_ptr vA2 = fva2->vertex();
      vertex_ptr vB2 = fvb2->vertex();

      edge_ptr ea0 = fva0->edge();
      edge_ptr eb0 = fvb0->edge();

      edge_ptr ea1 = fva1->edge();
      edge_ptr eb1 = fvb1->edge();

      edge_ptr ea2 = fva2->edge();
      edge_ptr eb2 = fvb2->edge();

      std::cout << " vA: " << vA0->position_in_set() << " "
                << vA1->position_in_set() << " " << vA2->position_in_set()
                << std::endl;
      std::cout << " vB: " << vB0->position_in_set() << " "
                << vB1->position_in_set() << " " << vB2->position_in_set()
                << std::endl;
      std::cout << " eA: " << ea0->position_in_set() << " "
                << ea1->position_in_set() << " " << ea2->position_in_set()
                << std::endl;
      std::cout << " eB: " << eb0->position_in_set() << " "
                << eb1->position_in_set() << " " << eb2->position_in_set()
                << std::endl;

      e_all.insert(ea0);
      e_all.insert(eb1);
      fv_all.insert(fva0);
      fv_all.insert(fvb1);
      v_del.insert(vB0);

      if (ea0 != eb1) {
        // std::cout << " -don't share " << vA0->position_in_set() <<
        // std::endl;

        vA0->set_front(fva0->vnext());
        e_border.push_back(ea0);
        face_vertex_ptr cva = fvb1->coedge();
        face_vertex_ptr cvb = fva0->coedge();
        ea0->set(cvb, cva);
        cva->vertex() = vA0;
        cvb->vertex() = vA1;

        v_new.insert(vA0);
        // std::cout << " foo: " << vA0->position_in_set() << std::endl;

        if (vA0 == vB0 && eb0 != ea2) {
          vertex_ptr vN = new vertex_type();
          this->_surf->push_vertex(vN);
          coordinate_type cN = this->coordinate(vA0);
          this->coordinate(cN, vN);
          face_vertex_ptr fvb0n = fvb0->vnext();
          vN->set_front(fvb0n); // this sets vertex
          e_border.push_back(fvb0->vnext()->edge());

          std::cout << " -ins v " << vN->position_in_set() << " "
                    << vN->get_front()->next()->vertex()->position_in_set()
                    << " "
                    << vA0->get_front()->next()->vertex()->position_in_set()
                    << " " << std::endl;
        }
      } else {

        if (fva0->vnext()->face() == fB)
          vA0->set_front(fva0->vprev());
        if (fva0->vprev()->face() == fB)
          vA0->set_front(fva0->vnext());

        e_border.push_back(ea0);
      }

      fvb0 = fvb0->prev();
    });

    std::cout << "es to keep: " << std::endl;
    std::for_each(e_border.begin(), e_border.end(),
                  [&e_all, &v_del, &fv_all, &v_new](edge_ptr e) {
                    e_all.erase(e);
                    face_vertex_ptr fv0 = e->v1();
                    face_vertex_ptr fv1 = e->v2();
                    vertex_ptr v0 = fv0->vertex();
                    vertex_ptr v1 = fv1->vertex();

                    v_del.erase(fv0->vertex());
                    v_del.erase(fv1->vertex());
                    v_new.insert(fv0->vertex());
                    v_new.insert(fv1->vertex());
                    fv_all.erase(fv0);
                    fv_all.erase(fv1);
                  });

    std::for_each(v_new.begin(), v_new.end(), [&mesh, &v_del](vertex_ptr v) {
      v_del.erase(v); // gaurd from deleting
      if (v->get_front()->vertex() != v) {
        std::cout << " bummer " << v->position_in_set() << " "
                  << v->get_front()->vertex()->position_in_set() << std::endl;
        v_del.insert(v);
      } else {
        v->update();
        std::cout << " cool.. " << v->position_in_set() << " "
                  << v->get_front()->vertex()->position_in_set() << std::endl;
      }
    });

    std::cout << "vr1: " << std::endl;
    std::for_each(v_del.begin(), v_del.end(), [&mesh, &v_new](vertex_ptr v) {
      std::cout << " : " << v->position_in_set() << std::endl;
    });

    std::cout << "fvs: " << std::endl;
    std::for_each(fv_all.begin(), fv_all.end(), [&mesh](face_vertex_ptr fv) {
      std::cout << " : " << fv << " " << fv->vertex()->position_in_set() << " "
                << fv->edge()->position_in_set() << " " << std::endl;
      delete fv;
    });

    std::cout << "er: " << std::endl;
    std::for_each(e_all.begin(), e_all.end(), [&mesh](edge_ptr e) {
      std::cout << " : " << e->position_in_set() << " " << e << std::endl;
      mesh->remove_edge(e->position_in_set());
    });

    std::cout << "vr: " << std::endl;
    std::for_each(v_del.begin(), v_del.end(), [&mesh, &v_new](vertex_ptr v) {
      std::cout << " : " << v->position_in_set() << std::endl;
      v_new.erase(v);
      mesh->remove_vertex(v->position_in_set());
    });

    // std::for_each(v_new.begin(), v_new.end(), [&mesh](vertex_ptr v) {
    //   // std::cout << "after" << std::endl;
    //   v->print();
    // });
    mesh->remove_face(fA->position_in_set());
    mesh->remove_face(fB->position_in_set());
  }

  std::vector<edge_ptr> get_bridged_edges(face_ptr fA, face_ptr fB) {
    std::vector<edge_ptr> collected;

    m2::set_vertex_flag<SPACE>(fA, 0, 0);
    m2::set_vertex_flag<SPACE>(fB, 0, 1);

    int idA = fA->position_in_set();
    int idB = fB->position_in_set();

    m2::for_each_face<SPACE>(fA, [&](face_vertex_ptr fv) {
      vertex_ptr v = fv->vertex();

      m2::for_each_vertex<SPACE>(v, [&](face_vertex_ptr fvv) {
        edge_ptr e = fvv->edge();
        face_ptr f0 = e->v1()->face();
        face_ptr f1 = e->v2()->face();

        if (fvv->edge()->is_degenerate()) {
          collected.push_back(fvv->edge());
          return;
        }

        if (fvv->next()->vertex()->flags[0] && //
            f0 != fA && f0 != fB &&            //
            f1 != fA && f1 != fB) {
          collected.push_back(fvv->edge());
          return;
        }
      });
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

    if (fvA0->edge() == fvB0->edge()) {
      fA->set_front(fvA0);
      fB->set_front(fvB2);
      return;
    }

    if (fvA0->edge() == fvB1->edge()) {
      fA->set_front(fvA0);
      fB->set_front(fvB0);
      return;
    }

    if (fvA0->edge() == fvB2->edge()) {
      fA->set_front(fvA0);
      fB->set_front(fvB1);
      return;
    }

    if (fvA1->edge() == fvB0->edge()) {
      fA->set_front(fvA1);
      fB->set_front(fvB2);
      return;
    }

    if (fvA1->edge() == fvB1->edge()) {
      fA->set_front(fvA1);
      fB->set_front(fvB0);
      return;
    }

    if (fvA1->edge() == fvB2->edge()) {
      fA->set_front(fvA1);
      fB->set_front(fvB1);
      return;
    }

    if (fvA2->edge() == fvB0->edge()) {
      fA->set_front(fvA2);
      fB->set_front(fvB2);
      return;
    }

    if (fvA2->edge() == fvB1->edge()) {
      fA->set_front(fvA2);
      fB->set_front(fvB0);
      return;
    }

    if (fvA2->edge() == fvB2->edge()) {
      fA->set_front(fvA2);
      fB->set_front(fvB1);
      return;
    }
    //-----------

    if (fvA0->vertex() == fvB0->vertex() && fvA1->vertex() == fvB1->vertex()) {
      fA->set_front(fvA0);
      fB->set_front(fvB0);
      return;
    }

    if (fvA0->vertex() == fvB1->vertex() && fvA1->vertex() == fvB2->vertex()) {
      fA->set_front(fvA0);
      fB->set_front(fvB1);
      return;
    }

    if (fvA0->vertex() == fvB2->vertex() && fvA1->vertex() == fvB0->vertex()) {
      fA->set_front(fvA0);
      fB->set_front(fvB2);
      return;
    }

    //-----------
    if (fvA1->vertex() == fvB1->vertex() && fvA2->vertex() == fvB2->vertex()) {
      fA->set_front(fvA1);
      fB->set_front(fvB1);
      return;
    }

    if (fvA1->vertex() == fvB2->vertex() && fvA2->vertex() == fvB0->vertex()) {
      fA->set_front(fvA1);
      fB->set_front(fvB2);
      return;
    }

    if (fvA1->vertex() == fvB0->vertex() && fvA2->vertex() == fvB1->vertex()) {
      fA->set_front(fvA1);
      fB->set_front(fvB0);
      return;
    }

    //-----------
    if (fvA2->vertex() == fvB2->vertex() && fvA0->vertex() == fvB0->vertex()) {
      fA->set_front(fvA2);
      fB->set_front(fvB2);
      return;
    }

    if (fvA2->vertex() == fvB0->vertex() && fvA0->vertex() == fvB1->vertex()) {
      fA->set_front(fvA2);
      fB->set_front(fvB0);
      return;
    }

    if (fvA2->vertex() == fvB1->vertex() && fvA0->vertex() == fvB2->vertex()) {
      fA->set_front(fvA2);
      fB->set_front(fvB1);
      return;
    }

    //-----------
    if (fvA0->vertex() == fvB0->vertex()) {
      fA->set_front(fvA0);
      fB->set_front(fvB0);
      return;
    }

    if (fvA0->vertex() == fvB1->vertex()) {
      fA->set_front(fvA0);
      fB->set_front(fvB1);
      return;
    }

    if (fvA0->vertex() == fvB2->vertex()) {
      fA->set_front(fvA0);
      fB->set_front(fvB2);
      return;
    }

    //-----------
    if (fvA1->vertex() == fvB0->vertex()) {
      fA->set_front(fvA1);
      fB->set_front(fvB0);
      return;
    }

    if (fvA1->vertex() == fvB1->vertex()) {
      fA->set_front(fvA1);
      fB->set_front(fvB1);
      return;
    }

    if (fvA1->vertex() == fvB2->vertex()) {
      fA->set_front(fvA1);
      fB->set_front(fvB2);
      return;
    }

    //-----------
    if (fvA2->vertex() == fvB0->vertex()) {
      fA->set_front(fvA2);
      fB->set_front(fvB0);
      return;
    }

    if (fvA2->vertex() == fvB1->vertex()) {
      fA->set_front(fvA2);
      fB->set_front(fvB1);
      return;
    }

    if (fvA2->vertex() == fvB2->vertex()) {
      fA->set_front(fvA2);
      fB->set_front(fvB2);
      return;
    }
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

    face_vertex_ptr fvAb = fvA0;
    T db = d0;

    fvAb = (db < d1) ? fvAb : fvA1;
    db = (db < d1) ? db : d1;
    fvAb = (db < d2) ? fvAb : fvA2;
    db = (db < d2) ? db : d2;

    fA->set_front(fvAb);
    fB->set_front(fvB0);
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
      int v = fv->vertex()->size();
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

    std::vector<face_ptr> filtered_faces;
    for (int i = 0; i < faces.size(); i++) {
      if (!this->_surf->has_face(i))
        continue;
      if (faces[i]->size() < 3)
        continue;
      filtered_faces.push_back(faces[i]);
    }

    std::vector<triangle_type> tris = m2::ci::get_tris<SPACE>(filtered_faces);
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

  void flag_edge(edge_ptr e) {

    face_vertex_ptr fva = e->v1();
    face_vertex_ptr fvb = e->v2();
    face_vertex_ptr va2 = fva->prev();
    edge_ptr ea2 = va2->edge();
    face_vertex_ptr vb2 = fvb->prev();
    edge_ptr eb2 = vb2->edge();
    e->flags[0] = 1;
    ea2->flags[0] = 1;
    eb2->flags[0] = 1;
    fva->face()->flags[1] = 1;
    fvb->face()->flags[1] = 1;

    vertex_ptr v1p = fva->prev()->vertex();
    vertex_ptr v2p = fvb->prev()->vertex();

    if (v1p->is_degenerate()) {
      for_each_vertex<SPACE>(v1p, [](face_vertex_ptr fv) {
        std::cout << " flagging deg" << fv->face()->position_in_set()
                  << std::endl;
        fv->face()->flags[1] = 1;
        fv->edge()->flags[1] = 1;
      });
    }

    if (v2p->is_degenerate()) {
      for_each_vertex<SPACE>(v2p, [](face_vertex_ptr fv) {
        std::cout << " flagging deg" << fv->face()->position_in_set()
                  << std::endl;
        fv->face()->flags[1] = 1;
        fv->edge()->flags[1] = 1;
      });
    }
  }

  std::vector<triangle_pair>
  collapse_bridges(std::vector<triangle_pair> &collected,
                   std::vector<face_ptr> &faces) {

    // collect bridges

    std::vector<triangle_pair> oldCollect = collected;

    bool topology_change = true;
    while (topology_change) {
      std::vector<edge_ptr> edgeToCollapse;

      for (auto p : oldCollect) {
        std::cout << " getting these dudes: " << p.A.faceId << " " << p.B.faceId
                  << std::endl;
        face_ptr fA = faces[p.A.faceId];
        face_ptr fB = faces[p.B.faceId];

        if (!fA || !fB)
          continue;

        std::vector<edge_ptr> i_bridges = get_bridged_edges(fA, fB);

        // copy edges to collection if they haven't been already been flagged
        std::copy_if(i_bridges.begin(), i_bridges.end(),
                     std::back_inserter(edgeToCollapse),
                     [fA, fB, this](edge_ptr e) {
                       bool flags = !(e->flags[0] | e->v1()->face()->flags[1] |
                                      e->v2()->face()->flags[1]);
                       this->flag_edge(e);
                       return flags;
                     });
      }

      topology_change = edgeToCollapse.size() > 0;
      std::cout << " edges to collapse " << edgeToCollapse.size() << std::endl;
      if (edgeToCollapse.size() > 0) {

        std::vector<triangle_pair> newCollect;
        std::copy_if(oldCollect.begin(), oldCollect.end(),
                     std::back_inserter(newCollect), [faces](triangle_pair p) {
                       face_ptr fA = faces[p.A.faceId];
                       face_ptr fB = faces[p.B.faceId];
                       return (!(fA->flags[1] | fB->flags[1]));
                     });

        std::cout << " new collect size: " << newCollect.size() << std::endl;
        oldCollect = newCollect;

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
          int changeId = v->topologyChangeId;
          if (changeId < 0)
            continue;
          m2::ci::set_coordinate<SPACE>(coords[changeId], v);
          v->topologyChangeId = -1;
        }

        std::cout << "done" << std::endl;
      }
    }

    std::cout << "new collect size: " << oldCollect.size() << std::endl;
    return oldCollect;
  }

  bool merge_collected(std::vector<triangle_pair> &collected,
                       std::vector<face_ptr> &faces) {

    bool merging = true;
    int phase = 0;
    while (merging && collected.size()) {

      // collected = collapse_bridges(collected, faces);

      std::sort(collected.begin(), collected.end(), std::less<triangle_pair>());

      std::cout << "merging " << phase
                << ", triangle pair count: " << collected.size() << std::endl;

      std::vector<triangle_pair> new_collect;

      for (auto p : collected) {
        if (!this->_surf->has_face(p.A.faceId) ||
            !this->_surf->has_face(p.B.faceId)) {
          std::cout << " moving on... " << std::endl;
          continue;
        };

        face_ptr fA = faces[p.A.faceId];
        face_ptr fB = faces[p.B.faceId];

        if (fA->is_degenerate() || fB->is_degenerate()) {
          continue;
        };
        // std::cout << " shared: " << shared << std::end
        coordinate_type NA = ci::normal<SPACE>(fA);
        coordinate_type NB = ci::normal<SPACE>(fB);
        real angle = va::dot(NA, NB);
        T d = p.dist();
        alignFaces(fA, fB);
        averageVerts(fA, fB);
        fA->flag_edges(2);
        fB->flag_edges(2);

        if (angle < -0.95) {

          joinFaces(fA, fB, this->_surf);
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

  bool merge() {
    std::vector<face_ptr> faces = this->_surf->get_faces();

    vector<triangle_pair> collected = get_face_pairs_to_merge(faces);
    return this->merge_collected(collected, faces);
  }

  void reset_flags() {

    for (auto f : this->_surf->get_faces()) {
      if (!f)
        continue;
      f->flags[0] = 0;
      f->flags[1] = 0;
    }

    for (auto v : this->_surf->get_vertices()) {
      if (!v)
        continue;
      v->flags[0] = 0;
      v->flags[1] = 0;
    }

    for (auto e : this->_surf->get_edges()) {
      if (!e)
        continue;
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

    std::cout << __FUNCTION__ << " " << e->v1()->face()->size() << " "
              << e->v2()->face()->size() << std::endl;
#if 1
    std::cout << "  " << __FUNCTION__ << std::endl;
    e->print();
#endif

    face_vertex_ptr fv1 = e->v1();
    face_vertex_ptr fv2 = e->v2();
    T circ1 = fv1->data;
    T circ2 = fv2->data;
    edge_split<SPACE> split;
    // coordinate_type c = getButterflyWeight(e);
    vertex_ptr nv = split(this->_surf, e);
    nv->topologyChangeId = i;

    // nv->coordinate() = c;
    fv1->edge()->flags[0] = 1;
    fv2->edge()->flags[0] = 1;

    edge_ptr e11 = subdivideFace(fv1->next());
    edge_ptr e21 = subdivideFace(fv2->next());

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

    std::cout << " - sorting " << edgesToSplit.size() << " edges" << std::endl;
    std::sort(edgesToSplit.begin(), edgesToSplit.end(),
              this->mEdgeSorterGreater);
    return edgesToSplit;
  }

  virtual bool op(edge_array edgesToSplit) {
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
      fv1->next()->edge()->flags[0] = true;
      fv2->next()->edge()->flags[0] = true;

      if (fv1->vertex()->size() < 3 || fv2->vertex()->size() < 3) {
        continue;
      }

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

  virtual bool op(edge_array edgeToCollapse) {
    // TIMER function//TIMER(__FUNCTION__);
    std::cout << " - deleting: " << edgeToCollapse.size() << " edges"
              << std::endl;

    for (int i = 0; i < edgeToCollapse.size(); i++) {
      edge_collapse<SPACE> collapser;
      if (!edgeToCollapse[i])
        continue;
      edge_ptr e = edgeToCollapse[i];
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

  bool op() {
    m2::construct<SPACE> cons;

    bool testing = true;
    int i;

    while (testing) {
      i = 0;
      bool deleting = false;
      for (auto e : this->_surf->get_edges()) {
        if (!this->_surf->has_edge(i++))
          continue;
        // std::cout << i << std::endl;
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

    m2::remesh<SPACE> rem;
    rem.triangulate(this->_surf);
  }
};
#endif

#if 1

template <typename SPACE> class vertex_policy {
public:
  M2_TYPEDEFS;

  vertex_policy(typename SPACE::vertex_index id) : _id(id) {}
  virtual void calc(int i, vertex_ptr &v0, vertex_ptr &v1) = 0;
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

  virtual void calc(int i, vertex_ptr &v0, vertex_ptr &v1) {
    TYPE t0 = v0->template get<TYPE>(this->_id);
    TYPE t1 = v1->template get<TYPE>(this->_id);
    _vals[i] = 0.5 * t0 + 0.5 * t1;
  }

  virtual void apply(int i, vertex_ptr &v) {
    v->template set<TYPE>(this->_id, _vals[i]);
  }

  virtual void reserve(long N) { _vals = std::vector<TYPE>(N); }

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

    double l0 = m2::ci::geometric_mean_length<SPACE>(_surf);
    l0 *= 0.5;
    std::cout << " integrator avg length: " << l0 << std::endl;
    _min = min * l0;
    _max = max * _min;
    _merge = merge * min;
  }

  void add_vertex_policy(vertex_policy<SPACE> *policy) {
    _policies.push_back(policy);
  }

  template <typename TYPE>
  void add_default_vertex_policy(typename SPACE::vertex_index id) {
    _policies.push_back(new vertex_policy_t<SPACE, TYPE>(id));
  }

  void cacheBary(m2::surf<SPACE> *surf) {
    M2_TYPEDEFS;
    int i = 0;
    for (auto f : surf->get_faces()) {

      if (!surf->has_face(i++))
        continue;

      if (f->size() < 3) {
        std::cout << " bary_size: " << std::endl;
        std::cout << "   face: " << f->size() << std::endl;
        std::cout << "   vert: " << f->fbegin()->vertex()->size() << std::endl;
        std::cout << "   edge: " << f->fbegin()->edge() << std::endl;
        f->print();
        f->fbegin()->vertex()->print();
        continue;
      }

      typename SPACE::coordinate_type c = m2::ci::center<SPACE>(f);
      typename SPACE::coordinate_type l = m2::ci::point_to_bary<SPACE>(f, c);

      assert(!isnan(l[0]));
      assert(!isnan(l[1]));
      assert(!isnan(l[2]));

      int i = 0;
      m2::for_each_face<SPACE>(f, [l, i](face_vertex_ptr fv) mutable {
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

    this->cacheBary(this->_surf);

    m2::area_laplacian<SPACE, coordinate_type> M(this->_surf);
    int i = 0;

    for (int k = 0; k < N; k++) {
      coordinate_array coords = m2::ci::get_coordinates<SPACE>(this->_surf);
      coordinate_array normals = m2::ci::get_vertex_normals<SPACE>(this->_surf);
      coordinate_array ncoords = M.mult(coords);

      for (int i = 0; i < coords.size(); i++) {
        coordinate_type cp = va::orthogonal_project(normals[i], ncoords[i]);

        if (isnan(cp[0])) {
          std::cout << "cp: " << cp.transpose() << std::endl;
          std::cout << "no: " << normals[i].transpose() << std::endl;
          std::cout << "nc: " << ncoords[i].transpose() << std::endl;
          std::cout << " c: " << coords[i].transpose() << std::endl;
        }

        coords[i] = coords[i] + C * cp;
        assert(!isnan(coords[i][0]));
        assert(!isnan(coords[i][1]));
        assert(!isnan(coords[i][2]));
      }
      m2::ci::set_coordinates<SPACE>(coords, this->_surf);
    }
  }

  void operate_edges(triangle_operations_base<SPACE> &op) {

    op.reset_flags();
    edge_array edges_to_op = op.get_edges();

    std::cout << "calculating new vals" << std::endl;
    for (auto p : _policies) {
      p->reserve(edges_to_op.size());
    }
    int i = 0;
    for (auto e : edges_to_op) {
      vertex_ptr v0 = e->v1()->vertex();
      vertex_ptr v1 = e->v2()->vertex();
      for (auto p : _policies) {
        p->calc(i, v0, v1);
      }
      i++;
    }
    std::cout << "op: " << edges_to_op.size() << " edges" << std::endl;

    op.op(edges_to_op);

    vertex_array vertices = _surf->get_vertices();

    std::cout << "assigning calc'd vertex values" << std::endl;

    for (auto v : vertices) {
      if (!v)
        continue;
      int topologyChangeId = v->topologyChangeId;
      if (topologyChangeId < 0)
        continue;

      for (auto p : _policies) {
        p->apply(topologyChangeId, v);
      }
      v->topologyChangeId = -1;
    }

    for (auto p : _policies) {
      p->clear();
    }
  }

  void splitEdges() {

    m2::edge_splitter<SPACE> splitter(_surf, _max);
    this->operate_edges(splitter);
  }

  void collapsEdges() {
    using namespace m2;
    M2_TYPEDEFS;
    std::cout << "-- collapsing --" << std::endl;
    m2::edge_collapser<SPACE> collapser(_surf, _min);
    this->operate_edges(collapser);
  }

  void mergeEdges() {
    std::cout << "-- merging --" << std::endl;
    edge_merger<SPACE> merger(this->_surf, _merge);
    bool merging = true;
    merging = merger.op();
  }

  void flipEdges() {
    edge_flipper<SPACE> flipper(this->_surf);
    flipper.op();
  }

  void delete_degenerates() {
    m2::construct<SPACE> cons;
    degenerate_deleter<SPACE> zapper(this->_surf);
    zapper.op();
  }

  bool integrate() { // mergeTris(_meshGraph, rxA, rxB);

    // mergeEdges();
    delete_degenerates();
    splitEdges();
    delete_degenerates();
    collapsEdges();
    delete_degenerates();
    flipEdges();
    delete_degenerates();
    std::cout << "ugly smoothing " << std::endl;
    this->_surf->pack();
    smoothMesh(0.01, 10);
    std::cout << "done " << std::endl;
  }
  m2::surf<SPACE> *_surf;
  real _min = 0.5, _max = 3.0, _merge = 1.0;

  std::vector<vertex_policy<SPACE> *> _policies;
};

#endif
}; // namespace m2
#endif
