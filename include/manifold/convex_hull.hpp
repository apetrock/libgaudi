/*
 *  convex_hull.hpp
 *  Manifold
 *
 *  Created by John Delaney on 7/14/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __TWOMANICONVEXHULL__
#define __TWOMANICONVEXHULL__

#include "construct.hpp"
#include "m2Includes.h"
#include "triangle_operations.hpp"
#include "remesh.hpp"
#include <time.h>

// functions-------------------------------

namespace m2 {
template <typename SPACE> class convex_hull {
  // class designed to produce primitives, such as convex hulls or loading .obj,
  // etc;
  /*TODO: doesn't work perfectly yet, there is still a bug somewhere, it might
   * be a precision thing */
  M2_TYPEDEFS
public:
  struct face_bundle {
    face_ptr cface;
    coordinate_list clist;
  };

  convex_hull(){};
  ~convex_hull(){};

  surf_ptr tetrahedron(coordinate_type cxu, coordinate_type cxl,
                       coordinate_type cyu, coordinate_type cyl,
                       coordinate_type czu, coordinate_type czl) {
    /*  czu
         /\
      cxl--cxu
         \/
         czl */
    surf_ptr ttet = new surf_type();
    construct<SPACE> cons;

    vector<coordinate_type> max_pts;
    max_pts.push_back(cxu);
    max_pts.push_back(cxl);
    max_pts.push_back(cyu);
    max_pts.push_back(cyl);
    max_pts.push_back(czu);
    max_pts.push_back(czl);

    ttet->insert_vertex(cxl);
    ttet->insert_vertex(cyu);
    ttet->insert_vertex(cxu);
    ttet->insert_vertex(cyl);

    cons.insert_edge(ttet, ttet->vertex(0), ttet->vertex(1));
    cons.insert_edge(ttet, ttet->vertex(1), ttet->vertex(2));
    cons.insert_edge(ttet, ttet->vertex(2), ttet->vertex(3));
    cons.insert_edge(ttet, ttet->vertex(3), ttet->vertex(0));
    remesh<SPACE> rem;
    ttet->pack();
    ttet->update_all();

    face_ptr f0 = ttet->face(0);
    face_ptr f1 = ttet->face(1);

    // std::cout << "stellating f0" << std::endl;
    // rem.stellate_face_generic(*ttet, *f0, czu);
    // ttet->pack(); ttet->update_all();
    // std::cout << "stellating f1" << std::endl;
    // rem.stellate_face_generic(*ttet, *f1, czl);
    // ttet->pack(); ttet->update_all();
    rem.triangulate(ttet);
    // ttet->pack();
    // ttet->update_all();

    std::cout << cxu << "\n"
              << cxl << "\n"
              << cyu << "\n"
              << cyl << "\n"
              << czu << "\n"
              << czl << std::endl;

    std::cout << ttet->get_vertices().size() << std::endl;

    std::cout << ttet->get_vertices().size() << std::endl;
    this->concave_flip(ttet);
    ttet->pack();
    ttet->update_all();
    return ttet;
  }

  void collect_points(face_ptr f, coordinate_list &clist,
                      coordinate_list &to_fill) {
    size_t s = clist.size();
    if (!clist.empty()) {
      typename coordinate_list::iterator itb = clist.begin();
      typename coordinate_list::iterator ite = clist.end();
      face_vertex_ptr fv1 = f->fbegin();
      face_vertex_ptr fv2 = fv1->next();
      face_vertex_ptr fv3 = fv2->next();
      std::cout << fv1 << " " << fv2 << " " << fv3 << std::endl;
      for (; itb != ite; itb++) {
        coordinate_type c4 = *itb;
        coordinate_type v41 = get_coordinate(fv1) - c4;
        coordinate_type v42 = get_coordinate(fv2) - c4;
        coordinate_type v43 = get_coordinate(fv3) - c4;
        T det = determinant(v41, v42, v43);
        if (det < 0.0) {
          to_fill.push_back(c4);
          itb = clist.erase(itb);
        }
      }
    }
  }

  void concave_flip(surf_ptr obj) {
    vector<edge_ptr> &edges = obj->get_edges();
    for (int i = 0; i < edges.size(); i++) {
      T det = edge_concave(edges[i]);
      std::cout << det << std::endl;
      if (det < 0) {
        std::cout << "flipping edges" << std::endl;
        edge_ptr e = edges[i];
        m2::triangle_operations<SPACE> triops;
        e = triops.flip(e);
      }
    }
  }

  T edge_concave(edge_ptr e) {
    coordinate_type v1 = get_coordinate(e->v1()->vertex());
    coordinate_type v2 = get_coordinate(e->v1()->prev()->vertex());
    coordinate_type v3 = get_coordinate(e->v2()->vertex());
    coordinate_type v4 = get_coordinate(e->v2()->prev()->vertex());

    T det = determinant(v2 - v1, v2 - v3, v2 - v4);
    return det;
  }

  void find_horizon(face_ptr &face_in, coordinate_type cen,
                    vector<face_ptr> &faces_to_remove, int face_flag) {
    vector<face_ptr> f_stack;
    face_in->flag = face_flag;
    f_stack.push_back(face_in);

    bool searching = true;
    while (f_stack.size() > 0) {
      face_ptr cface = f_stack.back();
      f_stack.pop_back();
      face_vertex_ptr fvb = cface->fbegin();
      face_vertex_ptr fve = cface->fend();

      faces_to_remove.push_back(cface);

      bool iterating = true;
      while (iterating) {
        iterating = fvb != fve;
        face_ptr oface = fvb->edge()->coface(fvb);
        face_vertex_ptr fv1 = oface->fbegin();
        face_vertex_ptr fv2 = fv1->next();
        ;
        face_vertex_ptr fv3 = fv2->next();
        ;
        coordinate_type c4 = cen;
        coordinate_type v41 = get_coordinate(fv1) - c4;
        coordinate_type v42 = get_coordinate(fv2) - c4;
        coordinate_type v43 = get_coordinate(fv3) - c4;

        T det1 = determinant(v41, v42, v43);
        if (det1 < 0) {
          if (oface->flag != face_flag) {
            oface->flag = face_flag;
            f_stack.push_back(oface);
          }
        }
        fvb = fvb->next();
      }
    }
  }

  surf_ptr quick_hull(std::vector<coordinate_type> clist) {
    coordinate_type czu, czl, cxu, cxl, cyu, cyl, avg;

    for (long i = 0; i < clist.size(); i++) {
      avg += clist[i];
    }

    avg /= (T)clist.size();

    cxu = avg;
    cxl = avg;
    cyu = avg;
    cyl = avg;
    czu = avg;
    czl = avg;
    for (long i = 0; i < clist.size(); i++) {
      if (clist[i][0] > cxu[0]) {
        cxu = clist[i];
      }
      if (clist[i][1] > cyu[1]) {
        cyu = clist[i];
      }
      if (clist[i][2] > czu[2]) {
        czu = clist[i];
      }

      if (clist[i][0] < cxl[0]) {
        cxl = clist[i];
      }
      if (clist[i][1] < cyl[1]) {
        cyl = clist[i];
      }
      if (clist[i][2] < czl[2]) {
        czl = clist[i];
      }
    }
    surf_ptr out = tetrahedron(cxu, cxl, cyu, cyl, czu, czl);

    // return out;
    //			man.merge_vertices(out);
    //			out->pack();

    // std::list<face_ptr>        face_stack;
    // std::list<coordinate_list> point_stack;
    std::list<face_bundle> face_stack;

    coordinate_list tlist; // temporary list of points
    for (long i = 0; i < clist.size(); i++) {
      tlist.push_back(clist[i]);
    }
    std::cout << "collecting points per face" << std::endl;
    for (int i = 0; i < out->get_faces().size(); i++) {
      coordinate_list tarr;
      face_bundle fb;
      fb.cface = out->face(i);
      collect_points(out->face(i), tlist, fb.clist);
      if (!fb.clist.empty()) {
        std::cout << "pushing back face bundle" << std::endl;
        face_stack.push_back(fb);
      }
    }

    long n = 0;
    std::cout << "entering traversal" << std::endl;
    while (face_stack.size() > 0) {
      //  		while (n < 0) {
      face_bundle cnode = face_stack.front();
      face_ptr cface = cnode.cface;
      coordinate_list &carr = cnode.clist;
      std::cout << cface << std::endl;
      // face_stack.pop_back();
      // point_stack.pop_back();

      face_vertex_ptr fv1 = cface->fbegin();
      face_vertex_ptr fv2 = fv1->next();
      face_vertex_ptr fv3 = fv2->next();
      coordinate_type vmax = carr.front();
      T dmax = 0;
      typename coordinate_list::iterator itb = carr.begin();
      typename coordinate_list::iterator ite = carr.end();
      std::cout << "looping through next coordinate list" << std::endl;
      while (itb != ite) {
        coordinate_type c1 = *itb;
        coordinate_type p1 = get_coordinate(fv1);
        coordinate_type p2 = get_coordinate(fv2);
        coordinate_type p3 = get_coordinate(fv3);
        T d = distance_from_plane(p1, p2, p3, c1);

        if (d > dmax) {
          dmax = d;
          vmax = *itb;
        }
        itb++;
      }

      remesh<SPACE> rem;
      vector<face_ptr> faces_in_horizon;
      int face_flag = 3;
      this->find_horizon(cface, vmax, faces_in_horizon, face_flag);
      std::cout << "found some faces to remove" << std::endl;
      if (faces_in_horizon.size() > 0) {

        // remove the old old faces from the set and merge their points with the
        // current face's
        std::cout << "removing faces and merging points size: "
                  << faces_in_horizon.size() << std::endl;
        for (int i = 0; i < faces_in_horizon.size(); i++) {
          face_ptr to_remove = faces_in_horizon[i];

          typename std::list<face_bundle>::iterator itb = face_stack.begin();
          typename std::list<face_bundle>::iterator ite = face_stack.end();

          bool iterating = true;
          while (iterating) {
            iterating = itb != ite;
            long t = (*itb).cface->ID();
            long t0 = to_remove->ID();
            std::cout << t << " " << t0 << std::endl;

            if (to_remove->ID() == (*itb).cface->ID()) {
              std::cout << "removing face" << std::endl;
              // coordinate_list tarr = (*itb);
              typename coordinate_list::iterator tatb = (*itb).clist.begin();
              typename coordinate_list::iterator tate = (*itb).clist.end();

              while (tatb != tate) {
                carr.push_back(*tatb);
                tatb++;
              }

              face_stack.erase(itb);
              iterating = false;
            }
            itb++;
          }
          // for(;itb!=ite;itb++){

          // }
        }

        face_ptr nf = rem.merge_faces(out, faces_in_horizon, face_flag);
        out->pack();
        out->update_all();
        vertex_ptr vc = NULL;
        vc = rem.stellate_face_generic(*out, *nf, vmax);
        out->pack();
        out->update_all();

        face_vertex_ptr fvb1 = vc->fbegin();
        face_vertex_ptr fve1 = vc->fend();
        bool iterating1 = true;

        while (
            iterating1) { // we need to capture all the faces we just stellated
          iterating1 = fvb1 != fve1;
          iterating1 &= fvb1 != NULL;
          if (fvb1 != NULL) {
            face_bundle fb;
            std::cout << "face: " << fvb1->face()->ID() << std::endl;
            fb.cface = fvb1->face();
            collect_points(fvb1->face(), carr, fb.clist);

            if (!fb.clist.empty()) {
              // fvb1->face()->color.b = 1;
              // fvb1->face()->color.r = 0;
              // fvb1->face()->color.g = 0;

              face_stack.push_back(fb);
            }
            fvb1 = fvb1->vnext();
          }
        }
      };
      out->reset_flags();
      std::cout << "next iteration" << std::endl;
      n++;
    }
    out->pack();
    out->update_all();
    return out;
  }
};
}; // namespace m2
#endif
