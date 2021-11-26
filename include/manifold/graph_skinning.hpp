/*
 *  graph_skinning.hpp
 *  Manifold
 *
 *  Created by John Delaney on 7/14/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __TWOMANIFOLDGRAPHSKIN__
#define __TWOMANIFOLDGRAPHSKIN__

#include <time.h>

#include "m2Includes.h"
#include "m2Operators.h"
#include "m1Includes.h"

//functions-------------------------------

namespace m2 {
  template <typename SPACE>
  class graph_skinning  {
    //class designed to produce primitives, such as convex hulls or loading .obj, etc;		
    M2_TYPEDEFS;
    public:
    T radius;
		
    graph_skinning(){};
    ~graph_skinning(){};
    surf_ptr build(m1::control<T>* graph_in){
      surf_ptr out = new surf_type();
      vector< face_list > faces_on_vertex;
      faces_on_vertex.resize(graph_in->get_vertices().size()); //some micromanagement here is necessary.
			
      int N2 = graph_in->size_edge();
      for(int i = 0; i < N2; i++){		
				
	coordinate_type de;
	coordinate_type vb = graph_in->E(i)->v1()->coordinate();		
	de = graph_in->E(i)->v2()->coordinate() - graph_in->E(i)->v1()->coordinate();
	typename m1::control<T>::edge_ptr ce = graph_in->E(i);
	ce->update_orientation();
				
	coordinate_type c1,c2,c3,c4,cp1,cp2,cp3,cp4;
	quat q1 = ce->orientation();
	quat q2;
				
	q2.w = 0;	q2.x = 0;	q2.y = 0;	q2.z = 1;		
	quat dq = q1-q2;
	dq.normalize();
	T r = radius;
	//T r = radius*al::rnd::uniform(1.0) + radius*0.5;				
	c1[0] = -r; c1[1] =  r; c1[2] = 0.0;
	c2[0] =  r; c2[1] =  r; c2[2] = 0.0;
	c3[0] =  r; c3[1] = -r; c3[2] = 0.0;
	c4[0] = -r; c4[1] = -r; c4[2] = 0.0;		
				
	cp1 = dq.rotate(c1);
	cp2 = dq.rotate(c2);
	cp3 = dq.rotate(c3);
	cp4 = dq.rotate(c4);
				
	T	jr = radius*8;
	//T	jr = 0.25;				
	coordinate_type dn = de;
	dn.normalize(); 
	T el = de.mag();
	//T el = 1.0;				
	cp1 += vb + dn*jr;	cp2 += vb + dn*jr;
	cp3 += vb + dn*jr;	cp4 += vb + dn*jr;
				
	m2::construct<SPACE>	cons;
	surf_ptr seg = new m2::surf<SPACE>();
	seg->insert_vertex(cp1);
	seg->insert_vertex(cp2);
	seg->insert_vertex(cp3);
	seg->insert_vertex(cp4);
				
	cons.insert_edge(seg, seg->vertex(0), seg->vertex(1));
	cons.insert_edge(seg, seg->vertex(1), seg->vertex(2));
	cons.insert_edge(seg, seg->vertex(2), seg->vertex(3));
	cons.insert_edge(seg, seg->vertex(3), seg->vertex(0));
				
	seg->pack();
				
	face_ptr f1 = seg->face(0);
	face_ptr f2 = seg->face(1);
				
	cons.bevel_face(seg, f1, 0, 0);
	m2::modify<SPACE> mod;
	/*
	  to finish out this portion, we'll have to put the beginning face and end face
	  into an array
	*/
	mod.translate_face_along_vector(f1, dn, el - jr);
	seg->update_all();
	out->merge(*seg);				
				
	//now we have to push this back into a list of edges
	long i1 = ce->v1()->position_in_set();
	long i2 = ce->v2()->position_in_set();
				
	face_list& fov1 = faces_on_vertex[i1]; //list because we will be pushing back data
	face_list& fov2 = faces_on_vertex[i2]; //list because we will be pushing back data
	//now all the faces are pushed into their respective vertices, then we can do some psuedo convex hulling
	//		vector<m1::edge<T>*> ev = graph_in->get_edges();
	//		vector<m1::vertex<T>*> vv = graph_in->get_vertices();		
	fov1.push_back(f2);
	fov2.push_back(f1);		
      }
			
      //for(long i = 1; i < 0; i++){
			for(long i = 0; i < faces_on_vertex.size(); i++){
	face_list& fov = faces_on_vertex[i];
	if (fov.size() > 1) {
	  fl_iterator itb = fov.begin();
	  fl_iterator ie  = fov.end();
					
	  std::vector<coordinate_type > jverts;		
	  while (itb != ie) {
	    face_ptr fi = *itb;
	    face_vertex_ptr fvb = fi->fbegin();
	    face_vertex_ptr fve = fi->fend();
	    bool iterating = true;
	    while (iterating) {
	      iterating = fvb != fve;
	      jverts.push_back(fvb->coordinate());
	      fvb = fvb->next();
	    }
	    itb++;
	  }
					
	  surf_ptr joint;
	  m2::convex_hull<SPACE> ch;
	  joint = ch.quick_hull(jverts);
					
	  remesh<SPACE> rem;
	  rem.merge_all_adjacent_planar_faces(joint, 0.005);
	  face_array jfaces = joint->get_faces();					
	  out->merge(*joint);
	  //					
	  std::list<std::pair<face_ptr,face_ptr> > jpairs;
	  itb = fov.begin();
	  int sfov = fov.size();
	  while (itb != ie) {
	    face_ptr fi = *itb;
	    for (long j = 0; j<jfaces.size(); j++) {
	      face_ptr fj = jfaces[j];
	      if (fi && fj) {
		bool is_match = opposing_face_match(fi, fj, 0.005);
		//								bool is_match = false;								
		if (is_match) {
		  pair<face_ptr,face_ptr> cpair(fi,fj);
		  jpairs.push_back(cpair);
		}
	      }
	    }
	    itb++;
	  }
					
	  typename std::list<std::pair<face_ptr,face_ptr> >::iterator pb = jpairs.begin();
	  typename std::list<std::pair<face_ptr,face_ptr> >::iterator pe = jpairs.end();
	  int N3 = 0;
	  while (pb != pe) {
	    face_ptr fi = (*pb).first;
	    face_ptr fj = (*pb).second;	
						
	    stitch_faces(out, fi, fj, 0.001);
	    pb++; N3++;// if(N3 > 0) pb = pe;
						
	    //pb++;
	  }
	  bool here = true;
	}
				
      }
      //remesh<T> rem;
      //rem.merge_all_adjacent_planar_faces(out, 0.25);
      vertex_array& verts = out->get_vertices();
      for(long i = 0; i < verts.size(); i++){
	if (verts[i]) {
	  if (verts[i]->get_face_vertices().size() == 0) {
	    out->remove_vertex(verts[i]->position_in_set());
	  }
	}
      }
			
      return out;
    }
		
    bool opposing_face_match(face_ptr f1, face_ptr f2, T tol){
      f1->update_all();
      f2->update_all();			
      coordinate_type n1 = f1->normal();
      coordinate_type n2 = f2->normal();
      n1.normalize();
      n2.normalize();
			
      bool
	norm_match = false,
	vert_match = false,
	size_match = false;
      T d = dist(n1,n2);
      if (d >= 2. - tol) {
	norm_match = true;
	if (f1->size() == f2->size()) {
	  size_match = true;
					
	  face_vertex_ptr fvb1 = f1->fbegin();
	  face_vertex_ptr fve1 = f1->fend();
	  face_vertex_ptr fvb2 = f2->fbegin();			
	  face_vertex_ptr fve2 = f2->fend();								
	  T normv = 0,normt = 0;
	  bool it1 = true;
					
	  while (it1) {
	    it1 = fvb2 != fve2;
	    normv = dist(fvb1->coordinate(),fvb2->coordinate());
	    if (normv < tol) it1 = false;
	    else fvb2 = fvb2->next();						
	  }
					
	  it1 = true;
	  while (it1) {
	    it1 = fvb1 != fve1;
	    normt += dist(fvb1->coordinate(),fvb2->coordinate());
	    fvb1 = fvb1->next();				
	    fvb2 = fvb2->prev();
	  }
					
	  if (normt/(T)f1->size() <= tol) {
	    vert_match = true;
	  }
	}				
				
      }
			
      bool out = false;
      if (norm_match && size_match && vert_match) {
	out = true;
      } else out = false;	
      return out;
    }
		
    bool stitch_faces(surf_ptr obj_in, face_ptr f1, face_ptr f2, T tol){
			
      face_vertex_ptr fvb2 = f2->fbegin();			
      face_vertex_ptr fve2 = f2->fend();								
      coordinate_type fc = f1->fbegin()->coordinate();
			
      T normv = 1000., normt = 1000.;
      bool it1 = true;
			
      while (it1) {
	it1 = fvb2 != fve2;
	normv = dist(fc,fvb2->coordinate());
	if (normv < normt)	{
	  normt = normv;
	  f2->fbegin() = fvb2;
	}
	fvb2	= fvb2->next();
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
				
	e1->set_other(fv1p,fv2n);	fv2n->edge() = e1;
				
	obj_in->remove_edge(e2->position_in_set());
	fvb2->vertex()->add_face_vertex(fv1p);
				
	//fv1p->vertex() = fvb2->vertex();
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
  };
}; //end m2
#endif
