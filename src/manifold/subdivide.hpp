//
//  subdivide.h
//  Manifold
//
//  Created by John Delaney on 5/3/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#ifndef __SUBDIVIDE__
#define __SUBDIVIDE__

#include "m2Includes.h"
#include "octree.hpp"
#include <cmath>
namespace m2 {
  template <typename SPACE>
  class subdivide{
    M2_TYPEDEFS
        
    public:
    subdivide(){}
    ~subdivide(){}
        
    void subdivide_catmull_clark(control_ptr& obj_in){
      coordinate_array	control;
      coordinate_array	edges;
      coordinate_array	centers;
      vertex_array		control_vertices;
      vertex_array		edge_vertices;
      vertex_array		center_vertices;
            
      catmull_crawl(obj_in, control, edges, centers);
      subdivide_control(obj_in, control_vertices, edge_vertices, center_vertices);
      //            update_manifold(control, edges, centers,
      //                            control_vertices, edge_vertices, center_vertices);
            
      //now update normals
      face_array&	rFaces	= obj_in->get_faces();
      fa_iterator fitb	= rFaces.begin();
      fa_iterator fite	= rFaces.end();
            
      while (fitb != fite) {
	face_ptr fi = (*fitb);
	fi->update_normal();
	fitb++;
      }
            
      obj_in->print();
            
    }
        
    int calc_vertex(vertex_ptr& vertex_in,
		    coordinate_type & out){
      size_t k = vertex_in->size(); //calculate the valence
      T B = 3./2./(T)k;
      T G = 1./4./(T)k;
      face_vertex_ptr fvb = vertex_in->fbegin();
      face_vertex_ptr fve = vertex_in->fend();
            
      coordinate_type cc = vertex_in->coordinate();

      out = cc*(1. - B - G);
      bool iterating = true;
      if (fvb == NULL) iterating == false; 
      int i = 0; int maxIt = 100;
      while (iterating && i < maxIt) {
	iterating = fvb != fve;
	iterating &= fvb != NULL;                
	if (fvb != NULL) {
	  face_vertex_ptr fvn = fvb->prev();
	  coordinate_type adj = fvn->coordinate();				
	  coordinate_type opp = fvn->prev()->coordinate();
	  out += adj*B/(T)k + opp*G/(T)k;
	  fvb = fvb->vnext();
	  i++;
	}
        
      }
      out[3] = 1;
      if(i == maxIt) {
	return 0;
	fvb->face()->color.r = 1.0;
	fvb->face()->color.g = 0.0;
	fvb->face()->color.b = 0.0;
      }
      else return 1;         
    }
        
    coordinate_type calc_edge(edge_ptr& edge_in){
      face_vertex_ptr  fv1 = edge_in->v1();
      face_vertex_ptr  fv2 = edge_in->v2();
      //c11     c21
      //c1  o   c2
      //c12     c22
      coordinate_type c1 = fv1->coordinate();
      coordinate_type c2 = fv2->coordinate();

            
      coordinate_type c11 = fv1->vnext()->next()->coordinate();
      coordinate_type c12 = fv1->prev()->coordinate();
            
      coordinate_type c21 = fv2->vnext()->next()->coordinate();
      coordinate_type c22 = fv2->prev()->coordinate();
            
      coordinate_type out  = (c11+c12+c21+c22)*1./16. + (c1 + c2)*3./8.;
      out[3] = 1;
      edge_in->flag = 1;
      return out;
            
    }
        
    coordinate_type calc_center(face_ptr& face_in){
      return face_in->calc_center();
      /*
      face_vertex_ptr  fvb = face_in->fbegin();
      face_vertex_ptr  fve = face_in->fend();
      T size = (T)face_in->size();
      bool iterating = true;
      coordinate_type out;
      while (iterating) {
	iterating = fvb != fve;
                
	face_vertex_ptr fvn = fvb->next();
	coordinate_type adj = fvn->coordinate();
	out += adj/size;
                
	fvb = fvb->next();				
      }		
      out[3] = 1.0;
      return out;
      */
    }
        
    coordinate_type calc_normal(coordinate_type& c0,coordinate_type& c1,
				coordinate_type& c2,coordinate_type& c3){
      coordinate_type out;
      out[0] += (c0[1] - c1[1])*(c0[2] + c1[2]);
      out[1] += (c0[2] - c1[2])*(c0[0] + c1[0]);
      out[2] += (c0[0] - c1[0])*(c0[1] + c1[1]);
            
      out[0] += (c1[1] - c2[1])*(c1[2] + c2[2]);
      out[1] += (c1[2] - c2[2])*(c1[0] + c2[0]);
      out[2] += (c1[0] - c2[0])*(c1[1] + c2[1]);
            
      out[0] += (c2[1] - c3[1])*(c2[2] + c3[2]);
      out[1] += (c2[2] - c3[2])*(c2[0] + c3[0]);
      out[2] += (c2[0] - c3[0])*(c2[1] + c3[1]);
            
      out[0] += (c3[1] - c0[1])*(c3[2] + c0[2]);
      out[1] += (c3[2] - c0[2])*(c3[0] + c0[0]);
      out[2] += (c3[0] - c0[0])*(c3[1] + c0[1]);
      T inv = Q_rsqrt(out[0]*out[0] + out[1]*out[1] + out[2]*out[2]);
      //            out = fnormalize(out);
      out = out*inv;

      return out;
    }
        
    control_ref subdivide_control(control_ref	control_in){
      std::cout << " ===subdividing=== " << std::endl;
      //now we make a new face for each edge, and look to the parent to see if the edge has been cracked or not
      //we've preallocated our new arrays, so we need to put all the new pieces back where they should be
      //vN = vi if original vertex location
      //vN = vs + ei      if edge vertex
      //vN = vs + es + fi if face_center
            
      //eN = 2*ei + j where j is the 0 or 1 depending on which edge it corresponds to
      //eN = 2*ei + fvi where fi is a counter of the face vertices as loop through
            
      //fN = fvi where fvi is a running face vertex index
            
      //we need to repack the control in, to get rid of null pointers after topological changes.
      //Shouldn't be a problem with the shader version, as its patch based and only renders non-null
      //patchs
            
      control_in.pack();
      control_in.update_all();			
      vertex_array & cverts = control_in.get_vertices();
      edge_array &   cedges = control_in.get_edges();
      face_array &   cfaces = control_in.get_faces();
            
      //because some faces have more than four vertices, we can't 
      //just multiply by four/two for everything. So for full
      //two manifold objects we can use the EP characteristic:
      //Euler-Poincare characteristic for 2-manifold objects:
      //V-E+F = 2;
            
      //during a subdivide we'll get a new vertex for every element, an edge yields a new vertex,
      //so does a face and we need to copy over the old vertics
      long vsize = cverts.size() + cedges.size() + cfaces.size();
      //we know that an edge doubles itself, but then has edges emmanating from itself to the center soo:
      long esize = 4*cedges.size();
      //now for number of faces
      //            long fsize = 2 - (vsize - esize);
      long fsize = 0;
      for (long i = 0; i < cfaces.size(); i++) {
	fsize+=cfaces[i]->size();
      }
      vertex_array nverts; nverts.resize(vsize);
      edge_array   nedges; nedges.resize(esize);
      face_array   nfaces; nfaces.resize(fsize);
            
      //loop through vertics and create all new vertices
            
      for(long i = 0; i < cverts.size(); i++){
	vertex_ptr nv = new vertex_type();
	//                nv->coordinate() = cverts[i]->coordinate();
	int err = calc_vertex(cverts[i], nv->coordinate());
        if(err == 0) return control_in;
				
	nv->position_in_set() = i;
	nverts[i] = nv;
      }
            
      //loop through edges and create all new vertices
      for(long i = 0; i < cedges.size(); i++){
	vertex_ptr nv = new vertex_type();
	//                nv->coordinate() = (cedges[i]->v1()->coordinate() + cedges[i]->v2()->coordinate())*0.5;
	nv->coordinate() = calc_edge(cedges[i]);
        nv->coordinate()[3] = 1;
        //std::cout << nv->coordinate().transpose() << std::endl;

	long setpos = cverts.size() + i;
	nv->position_in_set() = setpos;
	nverts[setpos] = nv;
                
	//now we push back some subdivided edges 
	//we'll keep the pointers null for now.
	edge_ptr e1 = new edge_type();
	edge_ptr e2 = new edge_type();
	nedges[2*i]     = e1;
	e1->position_in_set() = 2*i;
	nedges[2*i + 1] = e2;
	e2->position_in_set() = 2*i + 1;
      }
            
            
      long fcntr = 0;
      //loop through faces and create all new center face vertices
      for(long i = 0; i < cfaces.size(); i++){
	face_ptr cf = cfaces[i];
                
	if(cf){
	  vertex_ptr nv = new vertex_type();
	  // nv->coordinate() = cf->center();
	  nv->coordinate() = calc_center(cf);
          //std::cout << nv->coordinate().transpose() << std::endl;
          long setpos = cverts.size() + cedges.size() + i;
	  nv->position_in_set() = setpos;
	  nverts[setpos] = nv;                    
	}
	//now we push back some MORE subdivided edges 
	//we'll keep the pointers null for now, these correspond to the four face edges.
	//damn need iterator, don't like iterators, but it would be nice;                
	if (cf) { 
	  for (int j = 0; j < cf->size(); j++) {   

	    edge_ptr e1 = new edge_type();
	    long setpos = 2*cedges.size() + fcntr;
	    e1->position_in_set() = setpos;
	    nedges[setpos] = e1; fcntr++;
	  }
	}
      }
            
      fcntr = 0;
      //loop through faces and create new sub faces
      for(long i = 0; i < cfaces.size(); i++){
	bool iterating = true;
	face_ptr cf = cfaces[i];
	if (cf) {
	  face_vertex_ptr ftb = cf->fbegin();
	  face_vertex_ptr fte = cf->fend();
	  vertex_ptr vcent = nverts[cverts.size() + cedges.size() + i];
	  long j = 0;
                    
	  int fcntbg = fcntr;
	  while(iterating){ iterating = ftb != fte;
                        
	    face_ptr nf = new face_type();
                        
	    vertex_ptr vnext = nverts[cverts.size() + ftb->edge()->position_in_set()];
	    vertex_ptr vcurr = nverts[ftb->vertex()->position_in_set()];
	    vertex_ptr vprev = nverts[cverts.size() + ftb->prev()->edge()->position_in_set()];
                        
	    face_vertex_ptr fv0 = new face_vertex_type(); 
	    fv0->vertex() = vcurr; fv0->face() = nf; vcurr->add_face_vertex(fv0);
	    face_vertex_ptr fv1 = new face_vertex_type(); 
	    fv1->vertex() = vnext; fv1->face() = nf; vnext->add_face_vertex(fv1);
	    face_vertex_ptr fv2 = new face_vertex_type(); 
	    fv2->vertex() = vcent; fv2->face() = nf; vcent->add_face_vertex(fv2);
	    face_vertex_ptr fv3 = new face_vertex_type(); 
	    fv3->vertex() = vprev; fv3->face() = nf; vprev->add_face_vertex(fv3);
                        
	    fv0->next() = fv1; fv1->prev() = fv0;
	    fv1->next() = fv2; fv2->prev() = fv1;
	    fv2->next() = fv3; fv3->prev() = fv2;
	    fv3->next() = fv0; fv0->prev() = fv3;
	    nf->fbegin() = fv0;
	    //now lets find our edges: e0 and e3 come from the edge pile
	    //we need to allocate our interior vertices e1 and e2
	    //edges out:
	    long es = cedges.size();
	    long ei0,ei3;
                        
	    int e0num = cedges[ftb->edge()->position_in_set()]->vnum(ftb);
	    int e3num = cedges[ftb->prev()->edge()->position_in_set()]->vnum(ftb->prev());
                        
	    ei0 = 2*(ftb->edge()->position_in_set());
	    ei3 = 2*(ftb->prev()->edge()->position_in_set());
                        
	    edge_ptr e0;
	    if (e0num == 1) {
	      e0 = nedges[ei0+1];
	      e0->v1() = fv0;
	    } 
	    else {
	      e0 = nedges[ei0];
	      e0->v2() = fv0; 
	    }
                        
	    fv0->edge() = e0;
                        
	    edge_ptr e1 = nedges[2*es + fcntr];                
	    if (e1->v1()) e1->v2() = fv1; else e1->v1() = fv1;
	    fv1->edge() = e1;                        
                        
	    int fcntnxt = fcntr - 1;
	    if (j-1 < 0) fcntnxt = fcntbg + cf->size()-1;
	    edge_ptr e2 = nedges[2*es + fcntnxt];
	    if (e2->v1()) e2->v2() = fv2; else e2->v1() = fv2;
	    fv2->edge() = e2;
                        
	    edge_ptr e3; 
	    if (e3num == 1) {
	      e3 = nedges[ei3];
	      e3->v1() = fv3;
	    } 
	    else {
	      e3 = nedges[ei3+1];
	      e3->v2() = fv3; 
	    }
                        
	    fv3->edge() = e3;
	    nf->size() = 4;
                        
	    //                        nf->update_normal();
	    nf->normal() = calc_normal(fv0->coordinate(),
				       fv1->coordinate(),
				       fv2->coordinate(),
				       fv3->coordinate());

	    nf->position_in_set() = fcntr;
	    nfaces[fcntr] = nf; j++;
	    fcntr++;                       
	    ftb = ftb->next();                        
	  }
	}
      }
            
      control_ptr outptr = new control_type();
      control_ref out = *outptr;
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
    void subdivide_edges(control_ptr obj_in){
      obj_in->pack();
      edge_array& E = obj_in->get_edges();
      long sz = E.size();
      for (long i = 0; i < sz; i++) {
	subdivide_edge(obj_in, E[i]);
      }
    }
		
    vertex_ptr subdivide_edge(control_ptr	obj_in,
			      edge_ptr		edge_in){

      face_vertex_ptr fv1 =  edge_in->v1();
      face_vertex_ptr fv2 =  edge_in->v2();
      
      coordinate_type c1 = fv1->coordinate();
      coordinate_type c2 = fv2->coordinate();			
      coordinate_type cn = 0.5*(c1 + c2);
            
      vertex_ptr v1 = fv1->vertex();
      vertex_ptr v2 = fv2->vertex();
            
      vertex_ptr vn = new vertex_type(cn);
            
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

      e1->set(fv1,fv2n);
      e2->set(fv1n,fv2);

      e1->flag = 1; e2->flag = 1; vn->flag = 1;
            
      obj_in->push_vertex(vn);
      obj_in->push_edge(e2);

      fv1n->flag = 1;
      fv2n->flag = 1;
      if(v1->pinned ==true && v2->pinned == true){ 
      	vn->pinned = true;
      }
      else vn->pinned = false;

      return vn;
    }
  }; //class subdivide
}  //namespace M2
#endif
