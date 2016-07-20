/*
 *  m2control.hpp
 *  Phase Vocoder
 *
 *  Created by John Delaney on 1/8/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef __TWOMANIFOLDCONTROL__
#define __TWOMANIFOLDCONTROL__
#include "m2Common.hpp"
#include "m2Face.hpp"
#include "m2FaceVertex.hpp"
#include "m2Vertex.hpp"
#include "m2Edge.hpp"

namespace m2 {	
  template <typename SPACE>
  class control {				
    M2_TYPEDEFS
		
    public:
    control(){
      manual_clean_up = false;
    }
		
    control(control* rhs){
      manual_clean_up = false;
      //mFaces.clear();
      //mVertices.clear();
      //mEdges.clear();			
      mFaces.resize(rhs->mFaces.size());
      mVertices.resize(rhs->mVertices.size());          
      mEdges.resize(rhs->mEdges.size());
            
      for(long i = 0; i < rhs->mVertices.size(); i++){
	mVertices[i] = new vertex_type(rhs->mVertices[i]->coordinate());
	mVertices[i]->position_in_set() = i;
      }
            
      for(long i = 0; i < rhs->mEdges.size(); i++){
	mEdges[i] = new edge_type();
	mEdges[i]->position_in_set() = i;
      }
            
      for(long i = 0; i < rhs->mFaces.size(); i++){
	if (rhs->mFaces[i]) {
	  face_ptr nf =  new face_type();
	  nf->position_in_set() = i;
                    
	  face_ptr of = rhs->mFaces[i];
	  nf->size() = of->size();
	  face_vertex_ptr itb = of->fbegin();
	  face_vertex_ptr ite = of->fend();
	  face_vertex_ptr fv0,fv1;
                    
	  fv0 = new face_vertex_type(); 
	  vertex_ptr nv0 = mVertices[ite->vertex()->position_in_set()];
	  fv0->vertex() = nv0;
	  nv0->add_face_vertex(fv0);
                    
	  edge_ptr ne0 = mEdges[ite->edge()->position_in_set()];                
	  fv0->edge() = ne0;
	  fv0->face() = nf;
	  nf->fbegin() = fv0;
	  if(!ne0->v1()) ne0->v1() = fv0; else ne0->v2() = fv0;
                    
	  bool iterating = true;
	  while (iterating) {                                        
	    iterating = itb != ite->prev();
	    fv1 = new face_vertex_type(); 
	    vertex_ptr nv1 = mVertices[itb->vertex()->position_in_set()];
	    fv1->vertex() = nv1;
	    nv1->add_face_vertex(fv1);
	    edge_ptr ne1 = mEdges[itb->edge()->position_in_set()];                
	    fv1->edge() = ne1;
	    fv1->face() = nf;
                        
	    if(!ne1->v1()) ne1->v1() = fv1; else ne1->v2() = fv1;
                        
	    fv0->next() = fv1; fv1->prev() = fv0;
	    fv0 = fv1;
	    itb = itb->next();
	  }
	  fv1 = nf->fbegin();
	  fv0->next() = fv1; fv1->prev() = fv0;                    
	  mFaces[i] = nf;
	}
      }
      mFaceRecycle = rhs->mFaceRecycle;
      mVertexRecycle = rhs->mVertexRecycle;
      mEdgeRecycle = rhs->mEdgeRecycle;
      bool here = true;
    }
        
    ~control(){}		
		
    //		face<T>& get_face(size_t ind){
    //			return *mFaces[ind];
    //		}
		
    face_ptr&	face(size_t ind)	{return mFaces[ind];}		
    edge_ptr&	edge(size_t ind)	{return mEdges[ind];}		
    vertex_ptr&	vertex(size_t ind)	{return mVertices[ind];}
    bool has_face(size_t ind){if(mFaces[ind])return true; else return false;}
    bool has_edge(size_t ind){if(mEdges[ind])return true; else return false;}
    bool has_vertex(size_t ind){if(mVertices[ind])return true; else return false;}
		
    face_array&   get_faces(){return mFaces;}
    edge_array&   get_edges(){return mEdges;}
    vertex_array& get_vertices(){return mVertices;}
    void operator=(control& rhs){
      mFaces.clear();
      mVertices.clear();            
      mEdges.clear();
            
      for(long i = 0; i < rhs.mFaces.size(); i++){
	if (rhs.mFaces[i]) {
	  this->push_face(rhs.mFaces[i]);
	}
      }
      for(long i = 0; i < rhs.mEdges.size(); i++){
	if (rhs.mEdges[i]) {
	  this->push_edge(rhs.mEdges[i]);
	}
      }
      for(long i = 0; i < rhs.mVertices.size(); i++){
	if (rhs.mVertices[i]) {
	  this->push_vertex(rhs.mVertices[i]);
	}
      }
    }
	
    void merge(control_ref other){
      this->pack();
      other.pack();
      for(long i = 0; i < other.mFaces.size(); i++){
	if (other.mFaces[i]) {
	  this->push_face(other.mFaces[i]);
	}
      }
      for(long i = 0; i < other.mEdges.size(); i++){
	if (other.mEdges[i]) {
	  this->push_edge(other.mEdges[i]);
	}
      }
      for(long i = 0; i < other.mVertices.size(); i++){
	if (other.mVertices[i]) {
	  this->push_vertex(other.mVertices[i]);
	}
      }
      other.mFaces.clear();
      other.mEdges.clear();
      other.mVertices.clear();			
    }
		
    void push_vertex(vertex_ptr& in){
      if (mVertexRecycle.size() > 0 && !manual_clean_up){
	long i = mVertexRecycle.back();
	in->position_in_set() = i;
	mVertexRecycle.pop_back();
	mVertices[i] = in;
      }
      else{
	mVertices.push_back(in);
	in->position_in_set() = mVertices.size()-1;
      }
    }
		
    void push_edge(edge_ptr& in){			
      if (mEdgeRecycle.size() > 0 && !manual_clean_up){
	long i = mEdgeRecycle.back();
	in->position_in_set() = i;
	mEdgeRecycle.pop_back();
	mEdges[i] = in;
      }
      else{
	mEdges.push_back(in);
	in->position_in_set() = mEdges.size()-1;
      }
    }
		
    void push_face(face_ptr& in){
      if (mFaceRecycle.size() > 0 && manual_clean_up){
	long i = mFaceRecycle.back();
	in->position_in_set() = i;
	mFaceRecycle.pop_back();
	mFaces[i] = in;
      }
      else{
	mFaces.push_back(in);
	in->position_in_set() = mFaces.size()-1;
      }
    }
		
    void remove_vertex(long i){			
      vertex_ptr v = mVertices[i];
      mVertexRecycle.push_back(i);
      if(!manual_clean_up) {
	mVertices[i] = NULL; 
	delete v;
      }
    }
		
    void remove_edge(long i){
      edge_ptr e = mEdges[i];
      mEdgeRecycle.push_back(i);
      if(!manual_clean_up) {
	mEdges[i] = NULL; 
	delete e;
      }
    }
		
    void remove_face(long i){
      face_ptr f = mFaces[i];
      mFaceRecycle.push_back(i);
      if(!manual_clean_up) {
	mFaces[i] = NULL; 
	delete f; 
      };
    }
		
    void toggle_clean_up(){manual_clean_up ^= true;}

    vertex_ptr insert_vertex(coordinate_type in){
      return this->insert_vertex(in[0], in[1], in[2]);
    }
    


    vertex_ptr insert_vertex(T x,T y,T z)
    {
      vertex_ptr	new_vert = new vertex_type(x,y,z);
      new_vert->init();			
      face_ptr	new_face = new_vert->front()->face();
      this->push_face(new_face);
      this->push_vertex(new_vert);			
      return new_vert;
    }
		
    void clean_up(){
      //cleanup globally deletes pointers after an operation that needs them to exist.
      for(int i = 0; i < mFaceRecycle.size(); i++){
	int ii = mFaceRecycle[i];
	delete mFaces[ii]; mFaces[ii] = NULL;
      }
      for(int i = 0; i < mEdgeRecycle.size(); i++){
	int ii = mEdgeRecycle[i];
	delete mEdges[ii]; mEdges[ii] = NULL;
      }
      for(int i = 0; i < mVertexRecycle.size(); i++){
	int ii = mVertexRecycle[i];
	delete mVertices[ii]; mVertices[ii] = NULL;
      }
      this->toggle_clean_up();
    }

    void pack(){
      //TODO: safer pack, iterating from mRecycle[i] to mRecycle[i+1]
      if (mFaceRecycle.size() > 0) {
	face_array tFaces;
	long j = 0;
	for (long i = 0; i < mFaces.size(); i++) {
	  if(mFaces[i]){ 
	    tFaces.push_back(mFaces[i]);
	    tFaces.back()->position_in_set() = j;
	    j++;
	  }
	}
	swap(mFaces,tFaces);
      }
      mFaceRecycle.clear();
            
      if (mVertexRecycle.size() > 0) {
	vertex_array tVertices;
	long j = 0;
	for (long i = 0; i < mVertices.size(); i++) {
	  if(mVertices[i]){ 
	    //mVertices[i]->pack();
	    tVertices.push_back(mVertices[i]);
	    tVertices.back()->position_in_set() = j; j++;
	  }
	}
	swap(mVertices,tVertices);              
      }
      mVertexRecycle.clear();
            
      if (mEdgeRecycle.size() > 0) {
	edge_array tEdges;
	long j = 0;
	for (long i = 0; i < mEdges.size(); i++) {
	  if(mEdges[i]){
	    tEdges.push_back(mEdges[i]);
	    tEdges.back()->position_in_set() = j; j++;
	  }
	}
	swap(mEdges,tEdges);    
      }
      mEdgeRecycle.clear();
      this->pack_vertices();
    }
		
    void pack_vertices(){
      for(long i = 0; i < this->mVertices.size(); i++){
	//mVertices[i]->pack();
      }
    }
		
    void draw(){            
      //glDisable(GL_BLEND);
      //glEnable(GL_DEPTH | GL_DOUBLE | GLUT_RGB);	// Enables Depth Testing			
      this->draw_faces();
      //glDisable(GL_DEPTH_TEST);					// Enables Depth Testing 
      //glEnable(GL_BLEND);
    }
		
    bool draw(T off){
      glEnable(GL_COLOR_MATERIAL);
      this->draw_edges(off);			
      this->draw_faces(off);
      return true;
    }
		
    void draw_faces(){
      fa_iterator it_b = mFaces.begin();
      fa_iterator it_e = mFaces.end();
      while (it_b != it_e) {
	if(*it_b){
	  (*it_b)->draw();                
	}
	it_b++;
      }
    }
		
    void draw_faces(T off){
      fa_iterator it_b = mFaces.begin();
      fa_iterator it_e = mFaces.end();
      while (it_b != it_e) {
	if(*it_b){
	  (*it_b)->draw(off);
	  (*it_b)->draw_normal(off);
	}				
	it_b++;
      }
    }
		
		
    void draw_edges(T off){
      ea_iterator it_b = mEdges.begin();
      ea_iterator it_e = mEdges.end();
      size_t size = mEdges.size();
      it_b = mEdges.begin();
			
      while (it_b != it_e) {
	edge_ptr cur_edge = *it_b;				
	if(*it_b){ cur_edge->draw(off);
	}
	it_b++;
      }
			
    }

    void draw_vertex_colors(){
      fa_iterator it_b = mFaces.begin();
      fa_iterator it_e = mFaces.end();
      while (it_b != it_e) {
	if(*it_b){
	  (*it_b)->draw_vertex_colors();
	}				
	it_b++;
      }
    }

    void draw_vertices(T off){
      va_iterator it_b = mVertices.begin();
      va_iterator it_e = mVertices.end();
      while (it_b != it_e) {
	//(*it_b)->update_normal();
	(*it_b)->draw();
	//(*it_b)->draw_label();
	it_b++;
      }
    }
        
    void reset_flags(){
      for(long i = 0; i < mFaces.size(); i++){
	if (mFaces[i]) {
	  mFaces[i]->flag = 0;
	}
      }
            
      for(long i = 0; i < mEdges.size(); i++){
	if(mEdges[i]){ 
	  mEdges[i]->flag = 0;		
	}
      }
            
      for(long i = 0; i < mVertices.size(); i++){
	if(mVertices[i]){
	  mVertices[i]->flag = 0;
	}
      }			
    }

    void color_dead_pointers(){
      for(long i = 0; i < mFaces.size(); i++){
	if (mFaces[i]) {
	  if(!mFaces[i]->fend() || !mFaces[i]->fbegin()){
	    mFaces[i]->color.r = 1.0;
	    mFaces[i]->color.g = 0.0;
	    mFaces[i]->color.b = 0.0;
	    std::cout << "bad bad face" << std::endl;		  
	  }
	}
      }
            
      for(long i = 0; i < mEdges.size(); i++){
	if(mEdges[i]){ 
	  if(!mEdges[i]->v2()){
	    mEdges[i]->v1()->face()->color.r = 0.0;
	    mEdges[i]->v1()->face()->color.g = 0.5;
	    mEdges[i]->v1()->face()->color.b = 0.5;
	    std::cout << "bad edge v2 pntr" << std::endl;		  
	  }
	  if(!mEdges[i]->v1()){
	    mEdges[i]->v2()->face()->color.r = 0.0;
	    mEdges[i]->v2()->face()->color.g = 0.5;
	    mEdges[i]->v2()->face()->color.b = 0.0;
	    std::cout << "bad edge v1 pntr" << std::endl;		  
	  }
	}
      }
      for(long i = 0; i < mVertices.size(); i++){
	if(mVertices[i]){ 
	  vector<face_vertex_ptr>& fva = mVertices[i]->get_face_vertices();
	  for(int j = 0; j < fva.size(); j++){
	    face_vertex_ptr & fv = fva[j];
#if 1
	      if(!fv){
		// fv->face()->color.r = 0.20;
		// fv->face()->color.g = 0.10;
		// fv->face()->color.b = 0.5;
		break;
		}
	      face_vertex_ptr fve = fv->vprev();
	      if(!fve){
		// fv->face()->color.r = 0.30;
		// fv->face()->color.g = 0.40;
		// fv->face()->color.b = 0.1;
		break;
	      }
	      int i = 0; int maxIt = 100;
	      while(fv != fve && i < maxIt){
		if(!fv){
		  break;
		}
		if(!fv->prev()){
		  fv->face()->color.r = 1.0;
		  fv->face()->color.g = 0.5;
		  fv->face()->color.b = 0.0;
		  std::cout << "bad prev pntr" << std::endl;		  
		  break;
		}
		if(!fv->next()){
		  fv->face()->color.r = 1.0;
		  fv->face()->color.g = 0.0;
		  fv->face()->color.b = 0.5;
		  std::cout << "bad next pntr" << std::endl;		  
		  break;
		}
		if(!fv->vnext()){
		  fv->face()->color.r = 0.30;
		  fv->face()->color.g = 0.10;
		  fv->face()->color.b = 0.50;
		  std::cout << "bad vnext pntr" << std::endl;		  
		  break;
		}
		if(!fv->vprev()){
		  fv->face()->color.r = 0.6;
		  fv->face()->color.g = 0.5;
		  fv->face()->color.b = 0.2;
		  std::cout << "bad vprev pntr" << std::endl;		  
		  break;
		}
		if(i>maxIt/2){
		  fv->face()->color.r = 0.1;
		  fv->face()->color.g = 0.4;
		  fv->face()->color.b = 0.5;
		  fv->vertex()->color.r = 0.6;
		  fv->vertex()->color.g = 0.4;
		  fv->vertex()->color.b = 0.1;
		  std::cout << "vertex exceeded maximum iterations" << std::endl;
		  break;
		}
		fv = fv->vnext();
		i++;

	      }

#endif
	    }
	  }
	}
			
    }
    
        
    void update_all(){
      fa_iterator fit_b = mFaces.begin();
      fa_iterator fit_e = mFaces.end();
      for(long i = 0; i < mFaces.size(); i++){
	if (mFaces[i]) {
	  mFaces[i]->update_all();
	}
      }
            
      //			ea_iterator eit_b = mEdges.begin();
      //			ea_iterator eit_e = mEdges.end();
      //			size_t size = mEdges.size();
      //			while (eit_b != eit_e) {
      //				edge_ptr cur_edge = *eit_b;
      //				cur_edge->flag = 0;				
      //				eit_b++;
      //			}
            
      va_iterator vit_b = mVertices.begin();
      va_iterator vit_e = mVertices.end();
      while (vit_b !=vit_e) {
	if (*vit_b) {
	  (*vit_b)->update_normal();
	}
	vit_b++;
      }
            
    }
        
    void print(){
      cout << "----begin dump----" << endl;
      cout << "number of faces: " << mFaces.size() << endl;
      cout << "number of edges: " << mEdges.size() << endl;
      cout << "number of vertices: " << mVertices.size() << endl;
            
            
    }
        
    void print_stack(){
      cout << "----begin dump----" << endl;
      cout << "number of faces: " << mFaces.size() << endl;
      cout << "number of edges: " << mEdges.size() << endl;
      cout << "number of vertices: " << mVertices.size() << endl;
      size_t i = 0;
      // fa_iterator itb = mFaces.begin();
      // fa_iterator ite = mFaces.end();
      // while (itb != ite) {
      // 	if(*itb){
      // 	  cout << "Stack Number: "<< i << endl;
      // 	  i++;							
      // 	  face_ptr cur_face = *itb;
      // 	  cout << cur_face << endl;
      // 	  cur_face->print();
      // 	}
      // 	++itb;
      // }

      va_iterator itb = mVertices.begin();
      va_iterator ite = mVertices.end();
      while (itb != ite) {
	if(*itb){
	  cout << "Vertex Number: "<< i << endl;
	  i++;							
	  cout << (*itb)->coordinate() << endl;
	}
	++itb;
      }
      cout << "----end dump----" << endl;		
    }
        
    void print_edge(){
      cout << "----begin dump----" << endl;
      cout << "number of faces: " << mFaces.size() << endl;
      cout << "number of edges: " << mEdges.size() << endl;
      cout << "number of vertices: " << mVertices.size() << endl;
            
      ea_iterator itb = mEdges.begin();
      ea_iterator ite = mEdges.end();
      while (itb != ite) {
	//cout <<"pointer addres: "<< *itb << endl;
	face_ptr cur_face = (*itb)->v1()->face();
	cout << cur_face << endl;
	cur_face->print();
	++itb;
      }
      cout << "----end dump----" << endl;	
    }
        
  protected:
    bool manual_clean_up;
    face_array		mFaces;     vector<int> mFaceRecycle;
    edge_array		mEdges;     vector<int> mEdgeRecycle;
    vertex_array	mVertices;  vector<int>  mVertexRecycle;
        
  };	
}
#endif
