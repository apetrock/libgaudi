//
//  modify.hpp
//  Manifold
//
//  Created by John Delaney on 5/21/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//


#ifndef __M2EDGEPOINTSUBDIVIDE__
#define __M2EDGEPOINTSUBDIVIDE__
#include "m2Includes.h"
#include <cmath>

//oh man, 350 lines of code and nary a comment, I'll have to comment my way through this, blech!
namespace m2 {
  template <typename SPACE>
  class epSubdivide{
    M2_TYPEDEFS;
  public:
    void inset_vertex(control_ptr in, vertex_ptr v){
      face_ptr nf = new face_type();
      face_vertex_ptr fvb = v->fbegin();
      face_vertex_ptr fve = v->fend();
      std::vector<vertex_ptr> va;
      std::vector<face_vertex_ptr> fva;	
      std::vector<face_vertex_ptr> nfva;				
      bool iterating = true;
			
      while (iterating) {
	iterating = fvb != fve;
	coordinate_type de = fvb->next()->coordinate() - v->coordinate();
	coordinate_type nvc = v->coordinate() + 0.45*de;  //right now do 1/3 in next make it variable
	vertex_ptr nv = new vertex_type();
	nv->coordinate() = nvc;
	face_vertex_ptr nfv = new face_vertex_type();				
	va.push_back(nv);
	fva.push_back(fvb); //I pick this up in an array because topological changes ruin connectivity
	nfva.push_back(nfv); //I pick this up in an array because topological changes ruin connectivity				
	fvb = fvb->vnext();
      }
			
      int N = fva.size();
      for (int i = 0; i < N; i++) {
				
	int ip = i-1 > -1 ? i-1:N-1;				
	vertex_ptr Vertex0    = va[ip];
	vertex_ptr Vertex1    = va[i%N];
	vertex_ptr Vertex2    = va[(i+1)%N];				
				
	face_vertex_ptr newFace0 = nfva[ip];
	face_vertex_ptr newFace1 = nfva[i];	
	face_vertex_ptr newFace2 = nfva[(i+1)%N];					
				
	face_vertex_ptr oldFace0 = fva[i];
	face_vertex_ptr oldFace1 = new face_vertex_type();
				
	edge_ptr e1 = oldFace0->edge();
	edge_ptr e2 = new edge_type();				
				
	e2->v1() = oldFace0;
	e2->v2() = newFace2;
				
	oldFace0->edge() = e2;								
	newFace2->edge() = e2;
				
	in->push_edge(e2);
				
	e1->this_fv(oldFace0) = oldFace1;
	oldFace1->edge() = e1;
				
	newFace0->face() = nf;
	oldFace1->face() = oldFace0->face();
	face_vertex_ptr oldFace2 = oldFace0->next();
	oldFace1->next()    = oldFace2;	oldFace2->prev() = oldFace1;
	oldFace0->next()    = oldFace1;	oldFace1->prev() = oldFace0;
	newFace1->next()    = newFace0;	newFace0->prev() = newFace1;
				
	Vertex0->add_face_vertex(oldFace0);
	Vertex1->add_face_vertex(oldFace1);
	Vertex0->add_face_vertex(newFace1);
	oldFace0->face()->update_all();
				
	in->push_vertex(Vertex1);				
      }
      nf->fbegin() = nfva[0];
      nf->update_all();
      in->push_face(nf);
      in->remove_vertex(v->position_in_set());
      //in->update_all();			
    }
		
    void inset_vertex_triangle_collapse(control_ptr in, vertex_ptr v){
      face_ptr nf = new face_type();
      face_vertex_ptr fvb = v->fbegin();
      face_vertex_ptr fve = v->fend();
      std::vector<vertex_ptr> va;
      std::vector<face_vertex_ptr> fva;
      std::vector<face_vertex_ptr> nfva;				
      bool iterating = true;
			
      while (iterating) {
	iterating = fvb != fve;
	if (fvb->coface()->size() > 3) {
	  coordinate_type de = (fvb->next()->coordinate() + fvb->vnext()->next()->coordinate())*0.5 - v->coordinate();
	  coordinate_type nvc = v->coordinate() + 0.1*de;  //right now do 1/3 in next make it variable
	  vertex_ptr nv = new vertex_type();
	  nv->coordinate() = nvc;
	  face_vertex_ptr nfv = new face_vertex_type();
	  va.push_back(nv);
	  nfva.push_back(nfv); //I pick this up in an array because topological changes ruin connectivity				
	}
	fva.push_back(fvb);
	fvb = fvb->vnext();
      }
			
      int N  = fva.size();
      int Nj = va.size();
      int j = 0;
      for (int i = 0; i < N; i++) {
				
	int ip = i-1 > -1 ? i-1:N-1;
	int jp = j-1 > -1 ? j-1:Nj-1;					
	face_ptr nfi = new face_type();
				
	face_vertex_ptr oldFace0 = fva[i];
	face_vertex_ptr oldFace1 = fva[(i+1)%N];
	vertex_ptr Vertex0    = va[jp];
	vertex_ptr Vertex1    = va[j];
				
	if (oldFace0->coface()->size() > 3) {
				
	  vertex_ptr VertexN    = oldFace0->next()->vertex();
	  face_vertex_ptr newFace0 = nfva[jp];
	  face_vertex_ptr newFace1 = nfva[j];
	  face_vertex_ptr newFace2 = nfva[(j+1)%Nj];
					
	  face_vertex_ptr newFacei0 = new face_vertex_type();
	  face_vertex_ptr newFacei1 = new face_vertex_type();
	  face_vertex_ptr newFacei2 = new face_vertex_type();
					
	  edge_ptr e1 = oldFace0->edge();
	  edge_ptr e2 = new edge_type();
	  edge_ptr e3 = new edge_type();				
					
	  e2->v1() = newFace2;
	  e2->v2() = newFacei0;
	  newFace2->edge() = e2;								
	  newFacei0->edge() = e2;
	  in->push_edge(e2);
					
	  e3->v1() = oldFace1->prev();
	  e3->v2() = newFacei1;
	  oldFace1->prev()->edge() = e3;
	  newFacei1->edge() = e3;
	  in->push_edge(e3);
					
	  e1->other(oldFace0) = newFacei2;
	  newFacei2->edge() = e1;
					
	  oldFace0->face()->update_all();					
	  newFace1->next() = newFace0;	newFace0->prev() = newFace1;
					
	  newFacei0->next() = newFacei1;	newFacei1->prev() = newFacei0;
	  newFacei1->next() = newFacei2;	newFacei2->prev() = newFacei1;
	  newFacei2->next() = newFacei0;	newFacei0->prev() = newFacei2;
					
	  Vertex0->add_face_vertex(oldFace0);
	  Vertex0->add_face_vertex(newFace1);					
	  Vertex0->add_face_vertex(newFacei0);
	  Vertex1->add_face_vertex(newFacei1);
	  VertexN->add_face_vertex(newFacei2);
					
	  nfi->fbegin() = newFacei0;
	  nfi->update_all();
	  in->push_face(nfi);
	  in->push_vertex(Vertex0);
					
	  if (newFacei2->coface()->size() == 3) {
	    construct<T> cons;
	    cons.delete_non_cofacial_edge(in, newFacei2->edge());					
	  }
					
	  j++;
	}
	else {
	  Vertex0->add_face_vertex(oldFace0);
	}
				
      }
      nf->fbegin() = nfva[0];
      nf->update_all();
      in->push_face(nf);
      in->remove_vertex(v->position_in_set());
      //in->update_all();			
    }
		
    void inset_vertex_flag_collapse(control_ptr in, vertex_ptr v, vector<coordinate_type>& lc, long lcb){
      face_ptr nf = new face_type();
      face_vertex_ptr fvb = v->fbegin();
      face_vertex_ptr fve = v->fend();
      std::vector<vertex_ptr> va;
      std::vector<face_vertex_ptr> fva;
      std::vector<face_vertex_ptr> nfva;				
      bool iterating = true;			
			
			
      while (iterating) {
	iterating = fvb != fve;
	if (fvb->coface()->flag != 3) {
	  coordinate_type nvc = lc[lcb];
	  vertex_ptr nv = new vertex_type();
	  nv->coordinate() = nvc;
	  face_vertex_ptr nfv = new face_vertex_type();
	  va.push_back(nv);
	  nfva.push_back(nfv); //I pick this up in an array because topological changes ruin connectivity				
	  lcb++;
	}
	fva.push_back(fvb);
	fvb = fvb->vnext();
      }
			
      int N  = fva.size();
      int Nj = va.size();
      int j = 0;
      for (int i = 0; i < N; i++) {
				
	int ip = i-1 > -1 ? i-1:N-1;
	int jp = j-1 > -1 ? j-1:Nj-1;					
	face_ptr nfi = new face_type();
				
	face_vertex_ptr oldFace0 = fva[i];
	face_vertex_ptr oldFace1 = fva[(i+1)%N];
	vertex_ptr Vertex0    = va[jp];
	vertex_ptr Vertex1    = va[j];
				
	if (oldFace0->coface()->flag != 3) {
					
	  vertex_ptr VertexN    = oldFace0->next()->vertex();
	  face_vertex_ptr newFace0 = nfva[jp];
	  face_vertex_ptr newFace1 = nfva[j];
	  face_vertex_ptr newFace2 = nfva[(j+1)%Nj];
					
	  face_vertex_ptr newFacei0 = new face_vertex_type();
	  face_vertex_ptr newFacei1 = new face_vertex_type();
	  face_vertex_ptr newFacei2 = new face_vertex_type();
					
	  edge_ptr e1 = oldFace0->edge();
	  edge_ptr e2 = new edge_type();
	  edge_ptr e3 = new edge_type();				
					
	  e2->v1() = newFace2;
	  e2->v2() = newFacei0;
	  newFace2->edge() = e2;								
	  newFacei0->edge() = e2;
	  in->push_edge(e2);
					
	  e3->v1() = oldFace1->prev();
	  e3->v2() = newFacei1;
	  oldFace1->prev()->edge() = e3;
	  newFacei1->edge() = e3;
	  in->push_edge(e3);
					
	  e1->other(oldFace0) = newFacei2;
	  newFacei2->edge() = e1;					
	  nfi->flag = 3;
					
	  oldFace0->face()->update_all();					
	  newFace1->next() = newFace0;	newFace0->prev() = newFace1;
					
	  newFacei0->next() = newFacei1;	newFacei1->prev() = newFacei0;
	  newFacei1->next() = newFacei2;	newFacei2->prev() = newFacei1;
	  newFacei2->next() = newFacei0;	newFacei0->prev() = newFacei2;
					
	  Vertex0->add_face_vertex(oldFace0);
	  Vertex0->add_face_vertex(newFace1);					
	  Vertex0->add_face_vertex(newFacei0);
	  Vertex1->add_face_vertex(newFacei1);
	  VertexN->add_face_vertex(newFacei2);
					
	  nfi->fbegin() = newFacei0;
	  nfi->update_all();
	  in->push_face(nfi);
	  in->push_vertex(Vertex0);
					
	  if (newFacei2->coface()->flag == 3) {
	    construct<T> cons;
	    cons.delete_non_cofacial_edge(in, newFacei2->edge());					
	  }
					
	  j++;
	}
	else {
	  Vertex0->add_face_vertex(oldFace0);
	}
				
      }
      nf->fbegin() = nfva[0];
      nf->update_all();
      in->push_face(nf);
      in->remove_vertex(v->position_in_set());
      //in->update_all();			
    }
		
    bool subdivide_set(control_ptr in, list<vertex_ptr>& lv){
      vector<coordinate_type> cvlist;
      vector<long> cvbegin;			
      typename list<vertex_ptr>::iterator vi = lv.begin();
      typename list<vertex_ptr>::iterator ve = lv.end();
      bool iterating = true;
      while (vi != ve) {
	vertex_ptr v = *vi;
	face_vertex_ptr fvb = v->fbegin();
	face_vertex_ptr fve = v->fend();
	bool iterating2 = true;
	int N = 0;
	list<coordinate_type> cl;
	cvbegin.push_back(cvlist.size());
	while (iterating2) {
	  iterating2 = fvb != fve;
	  if (fvb->coface()->flag != 3) {
	    ////coordinate_type norm = fvb->face()->normal(); 
	    //coordinate_type de = (fvb->next()->coordinate()*0.10 + 
	    //				   fvb->vnext()->next()->coordinate()*0.90) - v->coordinate() - norm*0.05;
	    coordinate_type c1 = fvb->next()->coordinate();
	    coordinate_type c2 = fvb->vnext()->next()->coordinate();
	    coordinate_type de = (c1*0.5 + c2*0.5) - v->coordinate();
						
	    coordinate_type nvc = v->coordinate() + 0.15*de;  //right now do 1/3 in next make it variable
	    cvlist.push_back(nvc);
	  }
	  fvb = fvb->vnext();
	  N++;
	}
	vi++;
      }
      iterating = true;
      vi = lv.begin();
      int i = 0;
      while (vi != ve) {
	//			while (i<1) {	
	vertex_ptr v = *vi;
	inset_vertex_flag_collapse(in, v, cvlist, cvbegin[i]);
	vi++; i++;
      }
      in->reset_flags();
      return true;
    }
  };
}
#endif
