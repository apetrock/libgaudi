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
template <typename SPACE> class triangle_operations {
  M2_TYPEDEFS;

public:
  triangle_operations() {}
  ~triangle_operations() {}

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
