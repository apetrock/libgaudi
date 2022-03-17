//
//  remesh.hpp
//  Manifold
//
//  Created by John Delaney on 5/3/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

//
//  subdivide.h
//  Manifold
//
//  Created by John Delaney on 5/3/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#ifndef __REMESH__
#define __REMESH__
#include "m2Includes.h"
#include <cmath>

namespace m2 {
template <typename SPACE> class remesh : public default_interface<SPACE> {
  M2_TYPEDEFS;

public:
  remesh() {}
  ~remesh() {}
  void stellate(surf_ref control_in) {

    face_array &rFaces = control_in.get_faces();
    fa_iterator fitb = rFaces.begin();
    fa_iterator fite = rFaces.end();

    long sz = rFaces.size();
    for (long i = 0; i < sz; i++) {
      face_ptr fi = rFaces[i];
      if (fi)
        fi->flag = 12;
    }

    for (long i = 0; i < sz; i++) {
      face_ptr fi = rFaces[i];
      if (fi) {
        if (fi->flag == 12) {
          this->stellate_face_center(control_in, *fi);
          // control_in.remove_face(fi->position_in_set());
        }
      }
    }

    // control_in->reset_flags();
    control_in.update_all();
  }

  vertex_ptr stellate_face_center(surf_ptr obj_in, face_ptr face_in) {

    // 1. nitialize lists
    vertex_list vlist;
    face_vertex_list flist;

    // 2. calculate center vertex
    coordinate_array vl = face_in->flagged_coordinate_trace(0);
    coordinate_type cc = calculate_average(vl);
    return this->stellate_face_generic(obj_in, face_in, cc);
  }

  vertex_ptr stellate_face_generic(surf_ptr obj_in, face_ptr face_in,
                                   coordinate_type cen) {

    // 1. nitialize lists
    vertex_list vlist;
    face_vertex_list flist;

    // 2. calculate center vertex
    coordinate_type cc = cen;
    vertex_ptr v0 = new vertex_type();
    this->coordinate(cc,v0);
    obj_in->push_vertex(v0);

    face_vertex_ptr itb = face_in->fbegin();
    face_vertex_ptr ite = itb;

    int i = 0;

    std::vector<face_ptr> new_faces;
    bool iterating = true;
    while (iterating) {
      std::cout << " itb: " << itb << " " << ite << std::endl;
      face_vertex_ptr fv1 = itb;
      face_vertex_ptr fv2 = fv1->next();
      itb = itb->next();
      iterating = itb != ite;
      face_ptr nf;

      vertex_ptr v1 = fv1->vertex();
      vertex_ptr v2 = fv2->vertex();
      face_vertex_ptr p1 = fv1;
      face_vertex_ptr p0 = v0->make_face_vertex();
      face_vertex_ptr p2 = v2->make_face_vertex();

      nf = new face_type(p0, p1, p2);
      new_faces.push_back(nf);
      obj_in->push_face(nf);
    }

    for (int i = 0; i < new_faces.size(); i++) {
      int ip = (i + 1) % new_faces.size();
      face_ptr f0 = new_faces[i];
      face_ptr f1 = new_faces[ip];
      face_vertex_ptr p0 = f0->fbegin()->prev();
      face_vertex_ptr p1 = f1->fbegin();
      edge_ptr e = new edge_type(p0, p1);
      obj_in->push_edge(e);
    }

    for (int i = 0; i < new_faces.size(); i++) {
      face_ptr f0 = new_faces[i];
      f0->update_all();
    }

    obj_in->remove_face(face_in->position_in_set());

    return v0;
  }

  void root2(surf_ref control_in) {

    edge_array &rEdges = control_in.get_edges();
    for (long i = 0; i < rEdges.size(); i++) {
      if (rEdges[i])
        rEdges[i]->flag = 7;
    }

    this->stellate(control_in);

    for (long i = 0; i < rEdges.size(); i++) {
      edge_ptr ei = rEdges[i];
      if (rEdges[i]) {
        if (ei->flag == 7) {
          construct<SPACE> con;
          con.delete_non_cofacial_edge(control_in, *ei);
        }
      }
    }

    control_in.reset_flags();
    control_in.update_all();
  }

  bool merge_adjacent_planar_faces(surf_ptr obj_in, face_ptr face_in, T tol) {
    face_vertex_ptr fvb = face_in->fbegin();
    face_vertex_ptr fve = face_in->fend();
    coordinate_type n0 = face_in->normal();
    bool iterating = true;
    while (iterating) {
      iterating = fvb != fve;
      face_ptr coface = fvb->coface();
      coordinate_type n1 = coface->normal();
      T d = dist(n0, n1);
      if (d < tol) {
        edge_ptr e = fvb->edge();
        construct<SPACE> cons;
        face_ptr nf = cons.delete_non_cofacial_edge(*obj_in, *e);
        nf->color.g = 1;
        nf->color.b = 0;
        nf->color.r = 0;
        if (nf != face_in) {
          merge_adjacent_planar_faces(obj_in, nf, tol);
        }
        iterating = false;
      }
      fvb = fvb->next();
    }

    return true;
  }

  bool merge_all_adjacent_planar_faces(surf_ptr obj_in, T tol) {
    obj_in->reset_flags();
    face_array faces = obj_in->get_faces();
    for (int i = 0; i < faces.size(); i++) {
      //				face_array tFaces;
      if (obj_in->has_face(i)) {
        merge_adjacent_planar_faces(obj_in, faces[i], tol);
      }
      //				remesh<T> rem;
      //				face_ptr nf = merge_faces(obj_in,
      // tFaces); 				if (faces[i]->flag != 5) {
      // faces[i]->flag = 5; 					nf->flag = 5;
      //				}
    }
    obj_in->reset_flags();
    return true;
  }

  bool find_group_boundary(face_vertex_ptr cfv, face_ptr &tface,
                           vector<face_vertex_ptr> &outer_ring,
                           vector<face_vertex_ptr> &to_remove) {
    face_vertex_ptr fvb = cfv;
    face_vertex_ptr fve = cfv->prev();
    bool iterating = true;
    while (iterating) {
      iterating = fvb != fve;
      if (fvb->flag != 3) {
        fvb->flag = 3;
        face_ptr coface = fvb->coface();
        if (coface->flag == 3) {
          face_vertex_ptr fvn = fvb->vnext();
          to_remove.push_back(fvb);
          find_group_boundary(fvn, coface, outer_ring, to_remove);
        } else {
          outer_ring.push_back(fvb);
        }
      }
      fvb = fvb->next();
    }
    return true;
  }

  bool mark_adjacent_edges(vector<face_ptr> &faces_to_remove,
                           vector<edge_ptr> &edges_to_mark, int face_mark) {
    for (int i = 0; i < faces_to_remove.size(); i++) {
      face_vertex_ptr fvb = faces_to_remove[i]->fbegin();
      face_vertex_ptr fve = fvb->prev();
      bool iterating = true;
      while (iterating) {
        iterating = fvb != fve;
        if (fvb->flag != face_mark) {
          face_ptr coface = fvb->coface();
          if (coface->flag == face_mark) {
            std::cout << "marking adjacent edge: " << fvb->edge()->ID()
                      << std::endl;
            face_vertex_ptr fvn = fvb->vnext();
            fvb->edge()->flag = 0;

            edges_to_mark.push_back(fvb->edge());
          }
        }
        fvb = fvb->next();
      }
    }

    return true;
  }

  face_ptr merge_faces(surf_ptr obj_in, vector<face_ptr> &faces_to_remove,
                       int face_flag) {
    // why don't I loop through everything, then tag them locally.  Then find
    // the outside edges and do a recursive search to connect them up right now
    // this will assume all edges are connected, however detecting and mergeing
    // seperate loops will be in order eventually

    if (faces_to_remove.size() == 1) {
      return faces_to_remove.front();
    }

    for (long i = 0; i < faces_to_remove.size(); i++) {
      faces_to_remove[i]->flag = 3;
    }

    vector<edge_ptr> edges_to_remove;

    face_ptr nf = faces_to_remove[0];
    int edge_flag = 5;
    mark_adjacent_edges(faces_to_remove, edges_to_remove, face_flag);
    obj_in->toggle_clean_up(); // we need to make sure that pointers persist
                               // during this operation, so we'll clean up
                               // afterwards.
    for (long i = 0; i < edges_to_remove.size(); i++) {
      // obj_in->remove_vertex(edges_to_remove[i]->vertex()->position_in_set());
      edge_ptr e = edges_to_remove[i];
      if (e->flag != edge_flag) {
        std::cout << i << " ";
        e->flag = edge_flag;
        std::cout << "removing edge: " << e->ID() << std::endl;

        construct<SPACE> cons;
        nf = cons.delete_non_cofacial_edge(*obj_in, *e);
      }
    }
    nf->update_all();
    obj_in->clean_up();
    return nf;
    //			if (edge_ring.size() > 0) {
    //			}
  }

  surf_ref dual(surf_ref rhs) {
    rhs.pack();
    surf_ptr out = new surf_type();
    vertex_array &cverts = rhs.get_vertices();
    face_array &cfaces = rhs.get_faces();
    edge_array &cedges = rhs.get_edges();
    vertex_array nverts;
    face_array nfaces;
    edge_array nedges;

    nfaces.resize(cverts.size());
    nverts.resize(cfaces.size());
    nedges.resize(cedges.size());

    for (long i = 0; i < cverts.size(); i++) {
      if (cverts[i]->position_in_set() != i) {
        bool here = true;
      }
    }

    for (long i = 0; i < cfaces.size(); i++) {
      if (cfaces[i]->position_in_set() != i) {
        bool here = true;
      }
    }

    for (long i = 0; i < cedges.size(); i++) {
      if (cedges[i]->position_in_set() != i) {
        bool here = true;
      }
    }

    for (long i = 0; i < cverts.size(); i++) {
      nfaces[i] = new face_type();
      nfaces[i]->position_in_set() = i;
    }

    for (long i = 0; i < cfaces.size(); i++) {
      nverts[i] = new vertex_type(cfaces[i]->calculate_center());
      nverts[i]->position_in_set() = i;
    }

    for (long i = 0; i < nedges.size(); i++) {
      nedges[i] = new edge_type();
      nedges[i]->position_in_set() = i;
    }
    std::cout << "finished allocating new arrays" << std::endl;
    for (long i = 0; i < cverts.size(); i++) {

      face_ptr nf = nfaces[i];

      vertex_ptr ov = cverts[i];
      nf->size() = ov->size();
      face_vertex_ptr itb = ov->fbegin();
      face_vertex_ptr ite = ov->fbegin()->vnext();
      face_vertex_ptr fv0, fv1;

      fv0 = new face_vertex_type();
      vertex_ptr nv0 = nverts[ite->face()->position_in_set()];
      fv0->vertex() = nv0;
      nv0->add_face_vertex(fv0);

      edge_ptr ne0 = nedges[ite->edge()->position_in_set()];
      fv0->edge() = ne0;
      fv0->face() = nf;
      nf->set_front(fv0);
      if (!ne0->v1())
        ne0->v1() = fv0;
      else
        ne0->v2() = fv0;
      std::cout << "iterating through each face" << std::endl;
      bool iterating = true;
      while (iterating) {
        iterating = itb != ite->vnext();
        fv1 = new face_vertex_type();
        vertex_ptr nv1 = nverts[itb->face()->position_in_set()];
        fv1->vertex() = nv1;
        nv1->add_face_vertex(fv1);
        edge_ptr ne1 = nedges[itb->edge()->position_in_set()];
        fv1->edge() = ne1;
        fv1->face() = nf;

        if (!ne1->v1())
          ne1->v1() = fv1;
        else
          ne1->v2() = fv1;

        fv1->next() = fv0;
        fv0->prev() = fv1;
        fv0 = fv1;
        itb = itb->vprev();
      }
      fv1 = nf->fbegin();
      fv1->next() = fv0;
      fv0->prev() = fv1;
    }

    out->get_faces() = nfaces;
    out->get_edges() = nedges;
    out->get_vertices() = nverts;
    return *out;
  }

  void split4(surf_ptr &mMesh, face_ptr f) {
    face_vertex_ptr itb = f->fbegin();
    face_vertex_ptr ite = f->fend();
    while (itb->flag != 1) {
      itb = itb->next();
    }
    m2::construct<SPACE> cons;
    edge_ptr ne = cons.insert_edge(mMesh, itb, itb->next()->next());
    itb->flag = 0;
  }

  void split5(surf_ptr &mMesh, face_ptr f) {
    face_vertex_ptr itb0 = f->fbegin();
    while (itb0->flag != 1) {
      itb0 = itb0->next();
    }

    m2::construct<SPACE> cons;
    if (itb0->next()->next()->flag == 1) {
      face_vertex_ptr itb1 = itb0->next()->next();
      itb0->flag = 0;
      itb1->flag = 0;
      edge_ptr e1 = cons.insert_edge(mMesh, itb0, itb1);
      // itb1 = itb1->vprev();
      face_vertex_ptr itb2 = itb1->next()->next();
      edge_ptr e2 = cons.insert_edge(mMesh, itb1, itb2);
    } else {
      face_vertex_ptr itbt = itb0->prev()->prev();
      face_vertex_ptr itb1 = itb0;
      itb0 = itbt;
      itb0->flag = 0;
      itb1->flag = 0;
      edge_ptr e1 = cons.insert_edge(mMesh, itb0, itb1);
      // itb1 = itb1->vprev();
      face_vertex_ptr itb2 = itb1->next()->next();
      edge_ptr e2 = cons.insert_edge(mMesh, itb1, itb2);
    }
  }

  void split6(surf_ptr &mMesh, face_ptr f) {
    face_vertex_ptr itb0 = f->fbegin();
    while (itb0->flag != 1) {
      itb0 = itb0->next();
    }
    face_vertex_ptr itb1 = itb0->next()->next();
    face_vertex_ptr itb2 = itb1->next()->next();
    m2::construct<SPACE> cons;
    itb0->flag = 0;
    itb1->flag = 0;
    itb2->flag = 0;

    edge_ptr e1 = cons.insert_edge(mMesh, itb0, itb1);
    edge_ptr e2 = cons.insert_edge(mMesh, itb1, itb2);
    face_vertex_ptr itb0p = itb2->next()->next();
    edge_ptr e3 = cons.insert_edge(mMesh, itb0p, itb2);
  }

  void postSplitTriangulateFace(surf_ptr &mMesh, face_ptr f) {

    T maxArea = 0;
    T minDist = 99999;
    T maxDist = 0;
    bool iterating = true;
    face_vertex_ptr itb = f->fbegin();
    face_vertex_ptr ite = f->fend();
    face_vertex_ptr fvEdge = NULL;
#if 1
    int numFlags = 0;
    while (iterating) {
      iterating = itb != ite;
      if (itb->flag == 1)
        numFlags++;
      itb = itb->next();
    }

    if (numFlags == 1) {
      this->split4(mMesh, f);
    }
    if (numFlags == 2) {
      this->split5(mMesh, f);
    }
    if (numFlags == 3) {
      this->split6(mMesh, f);
    }
#endif
  }

  void triangulate_quad(surf_ptr &mMesh, face_ptr f) {
    f->update_all();

    T maxArea = 0;
    T minDist = 99999;
    T maxDist = 0;
    bool iterating = true;
    face_vertex_ptr v0 = f->fbegin();
    face_vertex_ptr v1 = v0->next();
    face_vertex_ptr v2 = v1->next();
    face_vertex_ptr v3 = v2->next();

    face_vertex_ptr fvEdge = NULL;
    T cot13 = this->cotan(v1) + this->cotan(v3); // corresponds to e13
    T cot02 = this->cotan(v0) + this->cotan(v2);
    if (cot13 * cot13 < cot02 * cot02) {
      m2::construct<SPACE> cons;
      edge_ptr ne = cons.insert_edge(mMesh, v1, v3);
    } else {
      m2::construct<SPACE> cons;
      edge_ptr ne = cons.insert_edge(mMesh, v0, v2);
    }
  }

  void slice_and_dice_face(surf_ptr &mMesh, face_ptr f) {
    f->update_all();

    T maxArea = 0;
    T minDist = 99999;
    T maxDist = 0;
    bool iterating = true;
    face_vertex_ptr p0 = f->fbegin();
    face_vertex_ptr p1 = p0->next()->next();
    face_vertex_array fvPairs;

    std::cout << f->size() << std::endl;

    while (iterating) {

      bool t1 = bool(p0->prev() == p1->next());
      bool t2 = bool(p0->prev() == p1);
      bool t3 = bool(p1->next() == p0);
      if (t1 || t2 || t3)
        iterating = false;
      if (p0->vertex() != p1->vertex()) {
        fvPairs.push_back(p0);
        fvPairs.push_back(p1);
      }
      p0 = p0->prev();
      p1 = p1->next();
    }
    m2::construct<SPACE> cons;
    for (int i = 0; i < fvPairs.size(); i += 2) {
      edge_ptr ne = cons.insert_edge(mMesh, fvPairs[i], fvPairs[i + 1]);
      std::cout << ne->v1()->face()->size() << " " << ne->v2()->face()->size()
                << std::endl;
    }
    // std::cout << minDist << std::endl;
  }

  void triangulate_face(surf_ptr mMesh, face_ptr f) {
    f->update_all();

    T maxArea = 0;
    T minDist = 99999;
    T maxDist = 0;
    bool iterating = true;
    face_vertex_ptr itb = f->fbegin();
    face_vertex_ptr ite = f->fend();
    face_vertex_ptr fvEdge = NULL;
    face_array fStack;
    fStack.push_back(f);
    int i = 0;
    while (fStack.size() > 0) {
      face_ptr fi = fStack.back();
      fStack.pop_back();
      if (fi->size() > 4) {
        m2::construct<SPACE> cons;
        face_vertex_ptr p0 = fi->fbegin();
        face_vertex_ptr p1 = p0->next()->next();
        if (i % 2 == 0)
          std::swap(p0, p1);
        edge_ptr ne = cons.insert_edge(mMesh, p0, p1);
        if (ne->v1()->face()->size() > 4) {
          fStack.push_back(ne->v1()->face());
        }
        if (ne->v2()->face()->size() > 4) {
          fStack.push_back(ne->v2()->face());
        }
        i++;
      }
    }
  }

  void triangulate(surf_ptr mMesh) {
    bool triangulating = true;
    int k = 0;
    while (triangulating && k < 100) {
      bool stillTriangulating = false;
      face_array &faces = mMesh->get_faces();
      int N = faces.size();
      for (int i = 0; i < N; i++) {
        if (faces[i]) {
          face_ptr f = faces[i];
          if (f->size() > 3) {
            if (f->size() == 4)
              this->triangulate_quad(mMesh, f);
            else
              this->triangulate_face(mMesh, f);
            // this->slice_and_dice_face(mMesh,f);
            stillTriangulating = true;
          }
        }
      }
      k++;
      triangulating = stillTriangulating;
    }
  }

  void flip_edges(surf_ptr obj) {
    bool triangulating = true;
    int k = 0;

    edge_array &edges = obj->get_edges();
    int N = edges.size();
    for (int i = 0; i < N; i++) {
      edge_ptr e = edges[i];
      if (e) {
        face_vertex_ptr v0 = e->v1();
        face_vertex_ptr v1 = v0->prev();
        face_vertex_ptr v2 = e->v2();
        face_vertex_ptr v3 = v2->prev();

        coordinate_type c0 = this->coordinate(v0);
        coordinate_type c1 = this->coordinate(v1);
        coordinate_type c2 = this->coordinate(v2);
        coordinate_type c3 = this->coordinate(v3);

        face_vertex_ptr fvEdge = NULL;
        T cos0 = dot((c1 - c0), (c3 - c0));
        T cos1 = dot((c0 - c1), (c2 - c1));
        T cos2 = dot((c1 - c2), (c3 - c2));
        T cos3 = dot((c0 - c3), (c2 - c3));

        T eFlip = cos0 + cos2; // corresponds to flipped edge
        T eSame = cos1 + cos3; // corresponds to flipped edge

        if (eFlip * eFlip > eSame * eSame) {
          m2::construct<SPACE> cons;

          e->v1()->vertex()->color.b = 0.75;
          e->v2()->vertex()->color.b = 0.75;
          cons.flip_edge(obj, e);
        }
      }
    }
  }

  void reverse_edge(edge_ptr &in) {
    face_vertex_ptr fv1 = in->v1();
    face_vertex_ptr fv2 = in->v2();
    face_vertex_ptr fv1n = fv1->next();
    face_vertex_ptr fv2n = fv2->next();
    in->v1() = fv2n;
    fv2n->edge() = in;
    in->v2() = fv1n;
    fv1n->edge() = in;
  }

  void reverse_face(face_ptr &in) {
    face_vertex_ptr fvb = in->fbegin();
    face_vertex_ptr fve = in->fend();
    bool iterating = true;
    while (iterating) {
      iterating = fvb != fve;
      face_vertex_ptr next = fvb->next();
      face_vertex_ptr prev = fvb->prev();

      fvb->next() = prev;
      fvb->prev() = next;
      fvb = next;
    }
  }

  void reverse(surf_ref rhs) {
    edge_array edges = rhs.get_edges();
    for (long i = 0; i < edges.size(); i++) {
      reverse_edge(edges[i]);
    }
    face_array faces = rhs.get_faces();
    for (long i = 0; i < faces.size(); i++) {
      reverse_face(faces[i]);
    }
  }
};
} // namespace m2
#endif
