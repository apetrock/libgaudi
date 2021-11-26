/*
 *  m1Construct.cpp
 *  Manifold
 *
 *  Created by John Delaney on 5/29/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
#include "m1.hpp"
#include "vec_addendum.h"
#include "m1Spring.hpp"
#ifndef __M1COSTRUCT__
#define __M1COSTRUCT__

namespace m1 {
  template <typename T>
  class construct {
    M1_TYPEDEFS
    public:
    void add_moment_support(surf_ptr& in){
      vertex_array& V = in->get_vertices();
      edge_array& E = in->get_edges();			
      vector<edge_ptr> collector;
      for (long i = 0; i < E.size(); i++){
	E[i]->k = 60.;
      }
      for (long i = 0; i < V.size(); i++) {
	vertex_ptr v = V[i];
				
	//				vector<coordinate_type> cv;
	//				for (vertex_iterator<T> itb = v->begin(); itb < v->end(); itb++) {
	//					coordinate_type tmp = V[i]->coordinate() - (*itb)->coordinate();
	//					cv.push_back(tmp.normalize());
	//				}
				
	//				coordinate_type N = calculate_average(cv);
	//				N.normalize();
				
	for (vertex_iterator<T> itb = v->begin(); itb < v->end(); itb++) {
	  for (vertex_iterator<T> jtb = v->begin(); jtb < v->end(); jtb++) {
	    //bool chk = jtb < jte;
	    vertex_ptr vi = *itb;
	    vertex_ptr vj = *jtb;
	    coordinate_type ci = vi->coordinate();
	    coordinate_type cj = vj->coordinate();
	    coordinate_type ct = V[i]->coordinate();
	    //						coordinate_type cti = orthogonal_project(N,ct-ci);
	    //						coordinate_type ctj = orthogonal_project(N,ct-cj);
	    coordinate_type cti = ct-ci;
	    coordinate_type ctj = ct-cj;						
	    T dotij = dot((cti).normalize(),(ctj).normalize());
	    if (vi != vj) {						
	      if (dotij > - 1.0 && dotij < -0.8){
		edge_ptr en = new edge_type();
		en->k = 60;
		en->v1() = vi; en->v2() = vj;
								
		collector.push_back(en);
	      }
	      //								if (dotij >   0.0 && dotij < 1.0){
	      //									edge_ptr en = new edge_type();
	      //									en->v1() = vi; en->v2() = vj;
	      //									collector.push_back(en);
	      //								}
							
	    }
	  }
	}
				
      }
      for (long i = 0; i < collector.size(); i++) {
	edge_ptr e = collector[i];
	e->v1()->push_singular_edge(e);
	e->v2()->push_singular_edge(e);
	in->push_edge(e);
      }	
			
    }
		
    void insert_vertex(surf_ptr ci, long ei){
      edge_ptr eo = ci->edge(ei);
      vertex_ptr v1 = eo->v1();
      vertex_ptr v2 = eo->v2();
      vertex_ptr vn = new vertex_type();
      coordinate_type cn = (v1->coordinate() + v2->coordinate())*0.5;
      vn->coordinate() = cn;
      edge_ptr en = new edge_type();
			
      eo->v2()  = vn;
			
      en->v1() = vn;
      en->v2() = v2;
      vn->push_edge(eo);
      vn->push_edge(en);
      ci->push_edge(en);
      ci->push_vertex(vn);
    }
		
    surf_ptr rand_graph(int N){
      surf_ptr out = new surf_type();
  
      for (int i = 0; i < N; i++){
	coordinate_type c;
	c[0]=(T)rand()/(T)RAND_MAX;
	c[1]=(T)rand()/(T)RAND_MAX;
	c[2]=(T)rand()/(T)RAND_MAX;
	m1::vertex<T>*  v = new m1::vertex<T>();
	v->coordinate() = c;
	out->push_vertex(v);
      }
			
			
      int N1 = out->size_vertex();
      for (int i = 0;i < N1; i++){
	int c1 = rand()%N1; int c2 = rand()%N1;
	if (rand()%100 > 70) {
	  out->connect(i,(i+1)%N1);
	}
      }
			
      vector<m1::edge<T>*> edc = out->get_edges();
      for(long i = 0; i < edc.size(); i++){
	m1::construct<T> cons;
	cons.insert_vertex(out, i);
      }
			
      N1 = out->size_vertex();
      for (int i = 0;i < N1; i++){
	int c1 = rand()%N1; int c2 = rand()%N1;
	if (rand()%100 > 11) {
	  if (c1!=c2)	out->connect(c1,c2);
	}
      }
			
			
      //	out->remove_duplicate_edges();
      N1 = out->size_vertex();
      for (long i = 0; i < N1; i++) {
	if (out->N(i)->size_edge() == 0) {
	  out->remove_vertex(i);
	}
      }
			
      edc = out->get_edges();
      for(long i = 0; i < edc.size(); i++){
	m1::construct<T> cons;
	cons.insert_vertex(out, i);
      }			
      out->pack();
			
      m1::spring_electric<T>* spr = new m1::spring_electric<T>(out);
      for(int i = 0; i < 400; i++){
	if (i%100 == 0) {
	  spr->reset_t();
	}
	spr->update_positions();
      }
      return out;
			
    }	
  };
};
#endif
