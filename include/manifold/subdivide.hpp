//
//  subdivide.h
//  Manifold
//
//  Created by John Delaney on 5/3/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#ifndef __SUBDIVIDE__
#define __SUBDIVIDE__

#include "construct.hpp"
#include "m2Includes.h"
#include "vec_addendum.h"

#include <cmath>
namespace m2 {
template <typename SPACE> class subdivide : public default_interface<SPACE> {
  M2_TYPEDEFS

public:
  subdivide() {}
  ~subdivide() {}

  void subdivide_catmull_clark(surf_ptr &obj_in) {
    coordinate_array control;
    coordinate_array edges;
    coordinate_array centers;
    vertex_array control_vertices;
    vertex_array edge_vertices;
    vertex_array center_vertices;

    catmull_crawl(obj_in, control, edges, centers);
    subdivide_control(obj_in, control_vertices, edge_vertices, center_vertices);
    //            update_manifold(control, edges, centers,
    //                            control_vertices, edge_vertices,
    //                            center_vertices);

    // now update normals
    face_array &rFaces = obj_in->get_faces();
    fa_iterator fitb = rFaces.begin();
    fa_iterator fite = rFaces.end();

    while (fitb != fite) {
      face_ptr fi = (*fitb);
      fi->update_normal();
      fitb++;
    }

    obj_in->print();
  }

  coordinate_type calc_vertex(vertex_ptr &vertex_in) {
    coordinate_type out(0, 0, 0);

    T B = 3.0 / 2.0;
    T G = 1.0 / 4.0;
    face_vertex_ptr fvb = vertex_in->fbegin();
    face_vertex_ptr fve = vertex_in->fend();

    coordinate_type cc = this->coordinate(vertex_in);

    // out = cc * (1. - B - G);
    bool iterating = true;
    if (fvb == NULL)
      iterating == false;
    int i = 0;
    int maxIt = 100;
    T k = 0.0;
    while (iterating && i < maxIt) {
      iterating = fvb != fve;
      iterating &= fvb != NULL;
      if (fvb != NULL) {
        face_vertex_ptr fvn = fvb->prev();
        coordinate_type adj = this->coordinate(fvn);
        coordinate_type opp = this->coordinate(fvn->prev());
        out += (B * adj + G * opp);
        fvb = fvb->vnext();
        k += 1.0;
        i++;
      }
    }
    out /= k * k;
    out += (1.0 - B / k - G / k) * cc;
    /*
    if (i == maxIt) {
      return 0;
      fvb->face()->color.r = 1.0;
      fvb->face()->color.g = 0.0;
      fvb->face()->color.b = 0.0;
    } else
      return 1;
    */
    return out;
  }

  coordinate_type calc_edge(edge_ptr &edge_in) {
    face_vertex_ptr fv1 = edge_in->v1();
    face_vertex_ptr fv2 = edge_in->v2();
    // c11     c21
    // c1  o   c2
    // c12     c22
    coordinate_type c1 = this->coordinate(fv1);
    coordinate_type c2 = this->coordinate(fv2);

    coordinate_type c11 = this->coordinate(fv1->vnext()->next());
    coordinate_type c12 = this->coordinate(fv1->prev());

    coordinate_type c21 = this->coordinate(fv2->vnext()->next());
    coordinate_type c22 = this->coordinate(fv2->prev());

    coordinate_type out =
        (c11 + c12 + c21 + c22) * 1. / 16. + (c1 + c2) * 3. / 8.;
    edge_in->flag = 1;

    return out;
  }

  coordinate_type calc_center(face_ptr &face_in) {
    return this->center(face_in);
  }

  surf_ref subdivide_control(surf_ref control_in) {
    std::cout << " ===subdividing=== " << std::endl;
    // now we make a new face for each edge, and look to the parent to see if
    // the edge has been cracked or not we've preallocated our new arrays, so we
    // need to put all the new pieces back where they should be vN = vi if
    // original vertex location vN = vs + ei      if edge vertex vN = vs + es +
    // fi if face_center

    // eN = 2*ei + j where j is the 0 or 1 depending on which edge it
    // corresponds to eN = 2*ei + fvi where fi is a counter of the face vertices
    // as loop through

    // fN = fvi where fvi is a running face vertex index

    // we need to repack the control in, to get rid of null pointers after
    // topological changes. Shouldn't be a problem with the shader version, as
    // its patch based and only renders non-null patchs

    control_in.pack();
    control_in.update_all();
    vertex_array &cverts = control_in.get_vertices();
    edge_array &cedges = control_in.get_edges();
    face_array &cfaces = control_in.get_faces();

    // because some faces have more than four vertices, we can't
    // just multiply by four/two for everything. So for full
    // two manifold objects we can use the EP characteristic:
    // Euler-Poincare characteristic for 2-manifold objects:
    // V-E+F = 2;

    // during a subdivide we'll get a new vertex for every element, an edge
    // yields a new vertex, so does a face and we need to copy over the old
    // vertics
    long vsize = cverts.size() + cedges.size() + cfaces.size();
    // we know that an edge doubles itself, but then has edges emmanating from
    // itself to the center soo:
    long esize = 4 * cedges.size();
    // now for number of faces
    //            long fsize = 2 - (vsize - esize);
    long fsize = 0;
    for (long i = 0; i < cfaces.size(); i++) {
      fsize += cfaces[i]->size();
    }
    vertex_array nverts;
    nverts.resize(vsize);
    edge_array nedges;
    nedges.resize(esize);
    face_array nfaces;
    nfaces.resize(fsize);

    // loop through vertics and create all new vertices

    for (long i = 0; i < cverts.size(); i++) {
      vertex_ptr nv = new vertex_type();
      coordinate_type c = calc_vertex(cverts[i]);
      this->coordinate(c, nv);
      // if (err == 0)
      //  return control_in;

      nv->position_in_set() = i;
      
      nverts[i] = nv;
    }

    // loop through edges and create all new vertices
    for (long i = 0; i < cedges.size(); i++) {
      vertex_ptr nv = new vertex_type();

      coordinate_type c = calc_edge(cedges[i]);
      this->coordinate(c, nv);

      long setpos = cverts.size() + i;
      nv->position_in_set() = setpos;
      nverts[setpos] = nv;

      // now we push back some subdivided edges
      // we'll keep the pointers null for now.
      edge_ptr e1 = new edge_type();
      edge_ptr e2 = new edge_type();
      nedges[2 * i] = e1;
      e1->position_in_set() = 2 * i;
      nedges[2 * i + 1] = e2;
      e2->position_in_set() = 2 * i + 1;
    }

    long fcntr = 0;
    // loop through faces and create all new center face vertices
    for (long i = 0; i < cfaces.size(); i++) {
      face_ptr cf = cfaces[i];

      if (cf) {
        vertex_ptr nv = new vertex_type();

        coordinate_type c = calc_center(cf);
        this->coordinate(c, nv);
        long setpos = cverts.size() + cedges.size() + i;
        nv->position_in_set() = setpos;
        nverts[setpos] = nv;
      }
      // now we push back some MORE subdivided edges
      // we'll keep the pointers null for now, these correspond to the four face
      // edges. damn need iterator, don't like iterators, but it would be nice;
      if (cf) {
        for (int j = 0; j < cf->size(); j++) {

          edge_ptr e1 = new edge_type();
          long setpos = 2 * cedges.size() + fcntr;
          e1->position_in_set() = setpos;
          nedges[setpos] = e1;
          fcntr++;
        }
      }
    }

    fcntr = 0;
    // loop through faces and create new sub faces
    for (long i = 0; i < cfaces.size(); i++) {
      bool iterating = true;
      face_ptr cf = cfaces[i];
      if (cf) {
        face_vertex_ptr ftb = cf->fbegin();
        face_vertex_ptr fte = cf->fend();
        vertex_ptr vcent = nverts[cverts.size() + cedges.size() + i];
        long j = 0;

        int fcntbg = fcntr;
        while (iterating) {
          iterating = ftb != fte;

          face_ptr nf = new face_type();

          vertex_ptr vnext =
              nverts[cverts.size() + ftb->edge()->position_in_set()];
          vertex_ptr vcurr = nverts[ftb->vertex()->position_in_set()];
          vertex_ptr vprev =
              nverts[cverts.size() + ftb->prev()->edge()->position_in_set()];

          face_vertex_ptr fv0 = new face_vertex_type();
          fv0->vertex() = vcurr;
          fv0->face() = nf;
          vcurr->add_face_vertex(fv0);
          face_vertex_ptr fv1 = new face_vertex_type();
          fv1->vertex() = vnext;
          fv1->face() = nf;
          vnext->add_face_vertex(fv1);
          face_vertex_ptr fv2 = new face_vertex_type();
          fv2->vertex() = vcent;
          fv2->face() = nf;
          vcent->add_face_vertex(fv2);
          face_vertex_ptr fv3 = new face_vertex_type();
          fv3->vertex() = vprev;
          fv3->face() = nf;
          vprev->add_face_vertex(fv3);
          fv0->next() = fv1;
          fv1->prev() = fv0;
          fv1->next() = fv2;
          fv2->prev() = fv1;
          fv2->next() = fv3;
          fv3->prev() = fv2;
          fv3->next() = fv0;
          fv0->prev() = fv3;
          nf->fbegin() = fv0;
          // now lets find our edges: e0 and e3 come from the edge pile
          // we need to allocate our interior vertices e1 and e2
          // edges out:
          long es = cedges.size();
          long ei0, ei3;

          int e0num = cedges[ftb->edge()->position_in_set()]->vnum(ftb);
          int e3num =
              cedges[ftb->prev()->edge()->position_in_set()]->vnum(ftb->prev());

          ei0 = 2 * (ftb->edge()->position_in_set());
          ei3 = 2 * (ftb->prev()->edge()->position_in_set());

          edge_ptr e0;
          if (e0num == 1) {
            e0 = nedges[ei0 + 1];
            e0->v1() = fv0;
          } else {
            e0 = nedges[ei0];
            e0->v2() = fv0;
          }

          fv0->edge() = e0;

          edge_ptr e1 = nedges[2 * es + fcntr];
          if (e1->v1())
            e1->v2() = fv1;
          else
            e1->v1() = fv1;
          fv1->edge() = e1;

          int fcntnxt = fcntr - 1;
          if (j - 1 < 0)
            fcntnxt = fcntbg + cf->size() - 1;
          edge_ptr e2 = nedges[2 * es + fcntnxt];
          if (e2->v1())
            e2->v2() = fv2;
          else
            e2->v1() = fv2;
          fv2->edge() = e2;

          edge_ptr e3;
          if (e3num == 1) {
            e3 = nedges[ei3];
            e3->v1() = fv3;
          } else {
            e3 = nedges[ei3 + 1];
            e3->v2() = fv3;
          }

          fv3->edge() = e3;
          nf->size() = 4;

          nf->position_in_set() = fcntr;
          nfaces[fcntr] = nf;
          j++;
          fcntr++;
          ftb = ftb->next();
        }
      }
    }

    surf_ptr outptr = new surf_type();
    surf_ref out = *outptr;
    //			cfaces.clear();
    //			cedges.clear();
    //			cverts.clear();
    //			cfaces = nfaces;
    //			cedges = nedges;
    //			cverts = nverts;
    //            out->get_faces().clear();
    //            out->get_edges().clear();
    //            out->get_vertices().clear();
    //            swap(out->get_faces(),nfaces);
    //            swap(out->get_edges(),nedges);
    //            swap(out->get_vertices(),nverts);
    out.get_faces() = nfaces;
    out.get_edges() = nedges;
    out.get_vertices() = nverts;
    out.update_all();
    return out;
  }

  void subdivide_edges(surf_ptr obj_in) {
    obj_in->pack();
    edge_array &E = obj_in->get_edges();
    long sz = E.size();
    for (long i = 0; i < sz; i++) {
      m2::construct<SPACE> cons;
      cons.subdivide_edge(obj_in, E[i]);
    }
  }

}; // class subdivide
} // namespace m2
#endif
