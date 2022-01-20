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
#include <cmath>

namespace m2 {

template <typename SPACE> class edge_splitter {
  M2_TYPEDEFS;

public:
  vertex_ptr split(surf_ptr obj_in, edge_ptr edge_in) {

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

template <typename SPACE> class edge_flipper {
  M2_TYPEDEFS;

public:
  edge_flipper() {}
  ~edge_flipper() {}

  edge_ptr flip(edge_ptr e1) {
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

#if 0
template <typename SPACE> class merge_faces : public default_interface<SPACE> {
public:
  M2_TYPEDEFS;

  std::vector<edge_ptr> getSharedEdges(face_ptr fA, face_ptr fB) {

    std::vector<edge_ptr> edges;

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
        fvB = fvB->next();
        edge_ptr eB = fvB->edge();
        if (eA == eB)
          edges.push_back(eA);
      }
      fvA = fvA->next();
    }
    return edges;
  }

  void joinFaceOneSharedEdge(face_ptr fA, face_ptr fB, surf_ptr mesh) {
    std::vector<edge_ptr> edges = getSharedEdges(fA, fB);
    if (edges.empty())
      return;
    edge_ptr e = edges[0];

    std::cout << "shared_edge: " << fA->position_in_set() << " "
              << fB->position_in_set() << std::endl;
    face_vertex_ptr fvA0 = e->v1();
    face_vertex_ptr fvA1 = fvA0->next();
    face_vertex_ptr fvA2 = fvA1->next();

    face_vertex_ptr fvB0 = e->v2();
    face_vertex_ptr fvB1 = fvB0->prev();
    face_vertex_ptr fvB2 = fvB1->prev();

    fvA0->vertex()->remove_face_vertex(fvA0);
    fvA1->vertex()->remove_face_vertex(fvA1);
    fvA2->vertex()->remove_face_vertex(fvA2);

    fvB0->vertex()->remove_face_vertex(fvB0);
    fvB1->vertex()->remove_face_vertex(fvB1);
    fvB2->vertex()->remove_face_vertex(fvB2);

    vertex_ptr vA0 = fvA0->vertex();
    vertex_ptr vA2 = fvA2->vertex();
    vertex_ptr vB0 = fvB0->vertex();
    vertex_ptr vB1 = fvB1->vertex();

    vA0->front() = fvA0->vprev();
    vB0->front() = fvB0->vprev();

    vA2->add_face_vertex(fvB2->coedge());
    vA2->update_all();

    edge_ptr eA1 = fvA1->edge();
    edge_ptr eA2 = fvA2->edge();
    edge_ptr eB1 = fvB1->edge();
    edge_ptr eB2 = fvB2->edge();

    eA1->set(fvA1->coedge(), fvB1->coedge());
    eA2->set(fvA2->coedge(), fvB2->coedge());

    delete fvA0;
    delete fvA1;
    delete fvA2;
    delete fvB0;
    delete fvB1;
    delete fvB2;

    mesh->remove_edge(e->position_in_set());
    mesh->remove_edge(eB1->position_in_set());
    mesh->remove_edge(eB2->position_in_set());
    mesh->remove_face(fA->position_in_set());
    mesh->remove_face(fB->position_in_set());
    mesh->remove_vertex(vB1->position_in_set());

    vA0->verify();
    vA2->verify();
    vB0->verify();
    eA1->v1()->face()->verify();
    eA1->v2()->face()->verify();
    eA2->v1()->face()->verify();
    eA2->v2()->face()->verify();
  }

  void joinFaceNoSharedEdges(face_ptr fA, face_ptr fB, surf_ptr mesh) {
    std::cout << "no shared_edge: " << fA->position_in_set() << " "
              << fB->position_in_set() << std::endl;

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

    vA0->remove_face_vertex(fvA0);
    vA1->remove_face_vertex(fvA1);
    vA2->remove_face_vertex(fvA2);

    vB0->remove_face_vertex(fvB0);
    vB1->remove_face_vertex(fvB1);
    vB2->remove_face_vertex(fvB2);

    vA0->add_face_vertex(fvB0->coedge());
    vA1->add_face_vertex(fvB1->coedge());
    vA2->add_face_vertex(fvB2->coedge());

    vA0->update_all();
    vA1->update_all();
    vA2->update_all();

    edge_ptr eA0 = fvA0->edge();
    edge_ptr eA1 = fvA1->edge();
    edge_ptr eA2 = fvA2->edge();

    edge_ptr eB0 = fvB0->edge();
    edge_ptr eB1 = fvB1->edge();
    edge_ptr eB2 = fvB2->edge();

    eA0->set(fvA0->coedge(), fvB0->coedge());
    eA1->set(fvA1->coedge(), fvB1->coedge());
    eA2->set(fvA2->coedge(), fvB2->coedge());

    delete fvA0;
    delete fvA1;
    delete fvA2;
    delete fvB0;
    delete fvB1;
    delete fvB2;

    mesh->remove_edge(eB0->position_in_set());
    mesh->remove_edge(eB1->position_in_set());
    mesh->remove_edge(eB2->position_in_set());

    mesh->remove_vertex(vB0->position_in_set());
    mesh->remove_vertex(vB1->position_in_set());
    mesh->remove_vertex(vB2->position_in_set());

    mesh->remove_face(fA->position_in_set());
    mesh->remove_face(fB->position_in_set());

    eA0->v1()->face()->verify();
    eA0->v2()->face()->verify();
    eA1->v1()->face()->verify();
    eA1->v2()->face()->verify();
    eA2->v1()->face()->verify();
    eA2->v2()->face()->verify();
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

    fA->print_vert_ids();
    fB->print_vert_ids(true);

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

    std::cout << " distance: " << d0 << " " << d1 << " " << d2 << std::endl;

    face_vertex_ptr fvAb = fvA0;
    fvAb = (d1 < d0) && (d1 < d2) ? fvA1 : fvAb;
    fvAb = (d2 < d0) && (d2 < d1) ? fvA2 : fvAb;

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
      if (vA != vB) {
        coordinate_type cA = this->coordinate(vA);
        coordinate_type cB = this->coordinate(vB);
        this->coordinate(0.5 * (cA + cB), vA);
      }
      fvA = fvA->next();
      fvB = fvB->prev();
    }
  }

  void getNearest(surf_ptr mesh) {

    // using namespace nanogui;

    std::vector<face_ptr> faces = mesh->get_faces();
    face_ptr f = faces[0];

    std::vector<triangle_type> tris;

    std::for_each(faces.begin(), faces.end(), [&tris](const face_ptr &f) {
      std::vector<triangle_type> ftris = m2::ci::get_tris<SPACE>(f);
      tris.insert(tris.end(), ftris.begin(), ftris.end());
    });

    m2::aabb_tree<SPACE, triangle_type> tree(tris);

    std::cout << "tree.nodes.size(): " << tree.nodes.size() << std::endl;

    auto triDist = [](const triangle_type &A, const triangle_type &B) -> T {
      return A.dist(B);
    };

    T tol = 0.025;

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

    bool trimming = true;
    while (trimming) {
      T d = collected.back().dist();
      if (d > tol)
        collected.pop_back();
      else
        trimming = false;
    }
    std::cout << " trimmed list: " << collected.size() << std::endl;

    collected.erase(std::remove_if(collected.begin(), collected.end(),
                                   [faces](auto p) {
                                     face_ptr fA = faces[p.A.faceId];
                                     face_ptr fB = faces[p.B.faceId];

                                     if (fA->flag > 0)
                                       return true;
                                     if (fB->flag > 0)
                                       return true;
                                     fA->flag = 1;
                                     fB->flag = 1;
                                     return false;
                                   }),
                    collected.end());

    std::cout << " flagged duplicates list: " << collected.size() << std::endl;

    for (int i = 0; i < collected.size(); i++) {
      // /for (int i = 0; i < 3; i++) {
      auto p = collected[i];
      face_ptr fA = faces[p.A.faceId];
      face_ptr fB = faces[p.B.faceId];

      std::cout << i << " " << fA->size() << " " << fB->size() << std::endl;

      int shared = fA->count_shared_vertices(fB);

      if (shared == 1) {
        std::cout << " shared == 1, punting!" << std::endl;
        continue;
      }
      std::cout << " shared: " << shared << std::endl;

      alignFaces(fA, fB);
      averageVerts(fA, fB);

      // continue;

      fA->print_vert_ids();
      fB->print_vert_ids(true);
      fA->print_shared_edge();
      fB->print_shared_edge();

      fA->flag_vertices();
      if (shared == 0) {
        int adjacent_shared = fA->count_adjacent_shared_vertices(fB);
        std::cout << " adjacent shared: " << shared << std::endl;
        if (adjacent_shared == 0)
          joinFaceNoSharedEdges(fA, fB, mesh);

      } else if (shared == 2)
        joinFaceOneSharedEdge(fA, fB, mesh);
    };

    mesh->verify();
    mesh->update_all();
    mesh->pack();
  }
};

template <typename SPACE> class split_edges : public default_interface<SPACE> {
public:
  M2_TYPEDEFS;

  edge_ptr subdivideFace(face_vertex_ptr fv) {

    // coordinate_type data = fv->face()->data;
    // coordinate_type data2 = fv->face()->data2;

    // TIMER function//TIMER(__FUNCTION__);
    // assuming that the face vertex is the newly inserted one.
    face_vertex_ptr fv1 = fv; // this is the
    face_vertex_ptr fv2 = fv->next();
    face_vertex_ptr fv3 = fv2->next();
    face_vertex_ptr fv4 = fv3->next();
#if 1
    coordinate_type e21 = this->coordinate(fv2) - this->coordinate(fv1);
    coordinate_type e32 = this->coordinate(fv3) - this->coordinate(fv2);
    coordinate_type e43 = this->coordinate(fv4) - this->coordinate(fv3);
    coordinate_type e14 = this->coordinate(fv1) - this->coordinate(fv4);
    coordinate_type e13 = this->coordinate(fv1) - this->coordinate(fv3);

    T e1 = fv1->data;
    T e2 = fv2->data;
    T e3 = fv3->data;
    T e4 = fv4->data;
    coordinate_type gamma0 = e2 * e21 + e3 * e32 + e4 * e43 + e1 * e14;
    Eigen::MatrixXd E(5, 6);
    E = Eigen::MatrixXd::Zero(5, 6);
    for (int i = 0; i < 3; i++) {
      E(i, 0) = e21[i];
      E(i, 1) = e13[i];
      E(i, 2) = e32[i];
      E(i, 3) = e14[i];
      E(i, 4) = e43[i];
      E(i, 5) = -e13[i];
    }
    E(3, 0) = 1;
    E(3, 1) = 1;
    E(3, 2) = 1;
    E(4, 3) = 1;
    E(4, 4) = 1;
    E(4, 5) = 1;
    Eigen::VectorXd b(5);
    b(0) = gamma0[0];
    b(1) = gamma0[1];
    b(2) = gamma0[2];
    b(3) = 0;
    b(4) = 0;
    // b = E.transpose()*b;
    // Eigen::MatrixXd EtE(3,3);
    // EtE = E.transpose()*E;
    // Eigen::VectorXd g = EtE.ldlt().solve(b);
    Eigen::VectorXd g =
        E.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);

    m2::construct<SPACE> cons;
    edge_ptr enew = cons.insert_edge(mMesh, fv1, fv3);
    enew->v1()->data = g(1);
    enew->v2()->data = g(5);
    fv3->data = g(0);
    fv2->data = g(2);
    fv1->data = g(3);
    fv4->data = g(4);
#else
    m2::construct<SPACE> cons;
    edge_ptr enew = cons.insert_edge(mMesh, fv1, fv3);

    // enew->v1()->face()->data = e1;
    // enew->v2()->face()->data = e2;

    // enew->v1()->face()->data = e3;
    // enew->v2()->face()->data = e4;
    enew->v1()->data = 0.0;
    enew->v2()->data = 0.0;
#endif
    return enew;
  }

  void split_edge(edge_ptr e) {
    // TIMER function//TIMER(__FUNCTION__);
    face_vertex_ptr fv1 = e->v1();
    face_vertex_ptr fv2 = e->v2();
    T circ1 = fv1->data;
    T circ2 = fv2->data;
    edge_splitter<SPACE> split;
    // coordinate_type c = getButterflyWeight(e);
    vertex_ptr nv = split.split(mMesh, e);

    // nv->coordinate() = c;
    fv1->data = circ2;
    fv2->data = circ1;
    fv1->next()->data = circ2;
    fv2->next()->data = circ1;

    nv->winding =
        0.5 * (e->v1()->vertex()->winding + e->v2()->vertex()->winding);
    subdivideFace(fv1->next());
    subdivideFace(fv2->next());
  }

  struct edge_sort {
    bool operator()(edge_ptr ei, edge_ptr ej) {
      return (ci::length<SPACE>(ei) < ci::length<SPACE>(ej));
    }
  } mEdgeSorter;

  bool split_edges() {
    // TIMER function//TIMER(__FUNCTION__);
    edge_array &edges = mMesh->get_edges();
    vector<vertex_ptr> &verts = mMesh->get_vertices();

    face_array &faces = mMesh->get_faces();
    int N = faces.size();

    bool topology_change = false;
    vector<T> edgeWeights;
    vector<T> vertexWeights;
    mMesh->update_all();
    mMesh->reset_flags();
    mMesh->pack();

    vector<edge_ptr> edgesToSplit;
#if 1
    for (int i = 0; i < edges.size(); i++) {
      if (!mMesh->has_edge(i))
        continue;
      edge_ptr ei = edges[i];
      bool pinned = ei->v1()->vertex()->pinned == true &&
                    ei->v2()->vertex()->pinned == true;
      if (pinned)
        continue;
      T l = this->length(ei);
      if (l > 1.75 * minLength) {
        edgesToSplit.push_back(ei);
        continue;
      }
    }
#endif

    std::cout << " - sorting " << edgesToSplit.size() << " edges, ";
    std::sort(edgesToSplit.begin(), edgesToSplit.end(), mEdgeSorter);
    std::cout << "splitting " << edgesToSplit.size() << " edges" << std::endl;
    for (int i = edgesToSplit.size(); i > 0; i--) {
      this->split_edge(edgesToSplit[i - 1]);
    }

    return topology_change;
  }
};

template <typename SPACE>
class collapse_edge : public default_interface<SPACE> {
public:
  M2_TYPEDEFS;

  bool collapse_edges() {
    // TIMER function//TIMER(__FUNCTION__);
    bool topology_change = false;
    edge_array collectedEdges;
    vector<edge_ptr> &edges = mMesh->get_edges();
    mMesh->reset_flags();
    for (int i = 0; i < edges.size(); i++) {

      if (!mMesh->has_edge(i))
        continue;
      edge_ptr e = edges[i];
      if (e->flag == 1)
        continue;
      if (e->v1()->vertex()->flag == 1)
        continue;
      if (e->v2()->vertex()->flag == 1)
        continue;
      if (e->v2()->vertex()->pinned)
        continue;
      T dist = this->length(e);
      if (dist < minCollapseLength) {
        e->v1()->vertex()->flag = 1;
        e->v2()->vertex()->flag = 1;
        collectedEdges.push_back(e);
      }
    }

    std::cout << " - deleting: " << collectedEdges.size() << " Tiny edges"
              << std::endl;

    for (int i = 0; i < collectedEdges.size(); i++) {
      construct<SPACE> cons;
      if (!collectedEdges[i])
        continue;
      edge_ptr e = collectedEdges[i];
      coordinate_type avg = 0.5 * (this->coordinate(e->v1()->vertex()) +
                                   this->coordinate(e->v2()->vertex()));
      this->coordinate(avg, e->v1()->vertex());
      // face_ptr nf = cons.delete_vertex_primitive(mMesh,e->v2()->vertex());
      // cons.collapse_edge_primitive(mMesh,e);
      cons.collapse_edge(mMesh, e);
      topology_change = true;
    }
    if (topology_change) {
      m2::remesh<SPACE> rem;
      rem.triangulate(mMesh);
      mMesh->pack();
      mMesh->update_all();
      mMesh->reset_flags();
    }
    return topology_change;
  }
};
#endif

}; // namespace m2
#endif
