/*
 *  m2Vertex.h
 *  Phase Vocoder
 *
 *  Created by John Delaney on 1/8/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef __TWOMANIFOLDVERTEX__
#define __TWOMANIFOLDVERTEX__
#include "m2Common.hpp"
#include "m2Control.hpp"
#include "m2Face.hpp"
#include "m2FaceVertex.hpp"
#include "m2Edge.hpp"
#include <list>

namespace m2 {
  template <typename SPACE>
  class vertex{
    M2_TYPEDEFS
    public:
    colorRGB color;

    vertex(){
      mNormal.zero();			
      m2::ID& manager = m2::ID::get_instance();
      mID = manager.new_vertex_id();
      mFaceVertices.clear();
      mRecycle.clear();
      flag = 0;
      color.r = 0.4 + randd(0.1);
      color.g = 0.4 + randd(0.1);
      color.b = 0.4 + randd(0.1);
      color.a = 1.0;
      pinned = false;
    }
		
    vertex(T x, T y, T z){
      mCoordinate[0] = x;
      mCoordinate[1] = y;
      mCoordinate[2] = z;
      mNormal.zero();
      mRecycle.clear();
      m2::ID& manager = m2::ID::get_instance();
      mID = manager.new_vertex_id();
      mFaceVertices.clear();
      flag = 0;
      color.r = 0.4 + randd(0.1);
      color.g = 0.4 + randd(0.1);
      color.b = 0.4 + randd(0.1);
      color.a = 1.0;
      pinned = false;
    }
		
    vertex(coordinate_type co){
      mCoordinate[0] = co[0];
      mCoordinate[1] = co[1];
      mCoordinate[2] = co[2];
      mNormal.zero();
      mRecycle.clear();
      m2::ID& manager = m2::ID::get_instance();
      mID = manager.new_vertex_id();
			
      mFaceVertices.clear();
      flag = 0;
      color.r = 0.4 + randd(0.1);
      color.g = 0.4 + randd(0.1);
      color.b = 0.4 + randd(0.1);
      color.a = 1.0;
      pinned = false;
    }		      
		
    void init(){
      mFaceVertices.clear();
			
      face_ptr		new_face = new face_type(*this);
      face_vertex_ptr new_fv = new_face->fbegin();
      size_t sz = mFaceVertices.size();
      new_fv->face()	= new_face;
      this->add_face_vertex(new_fv);
      pinned = false;
    }
		
    long& position_in_set(){return mSetPosition;}
		
    T & x()		  {return mCoordinate[0];}
    T   x()	const {return mCoordinate[0];}
    T & y()		  {return mCoordinate[1];}
    T   y()	const {return mCoordinate[1];}
    T & z()		  {return mCoordinate[2];}
    T   z()	const {return mCoordinate[2];}
    T & operator[](int i)		{return mCoordinate[i];}
    T   operator[](int i) const {return mCoordinate[i];}
		
    long	&	ID()			{	return mSetPosition;}
    long		ID()	const	{	return mSetPosition;}
		
    size_t size()	{	return mFaceVertices.size();	}
		
    coordinate_type	normal() const	{	return mNormal;}
		
    coordinate_type&	coordinate()		{	return mCoordinate;}
    coordinate_type		coordinate() const	{	return mCoordinate;}
        
    void add_face_vertex(face_vertex_ptr& new_fv){
      mFront = new_fv;
      mFaceVertices.push_back(new_fv);
      fvl_iterator end = mFaceVertices.end();
      end--;
      new_fv->position_in_vertex() = end;
      // if (mRecycle.size()>0) {
      // 	long ins = mRecycle.back();
      // 	mRecycle.pop_back();
      // 	mFaceVertices[ins] = new_fv;
      // 	new_fv->position_in_vertex() = ins;
      // } else {
      // 	mFaceVertices.push_back(new_fv);
      // 	new_fv->position_in_vertex() = mFaceVertices.size()-1;
      // }
      new_fv->vertex() = this;
      //mFaceVertices.push_back(new_fv);
    }        

    void remove_face_vertex(fvl_iterator & it){
      if(*it && mFaceVertices.size() > 0){
	fvl_iterator itb = mFaceVertices.begin();
	fvl_iterator ite = mFaceVertices.end();
	bool inset = false;
	for(;itb != ite; ++itb){
	  if(*it == *itb) inset = true;
	}
	if(inset){
	  mFaceVertices.erase(it);
	  //mFront = mFaceVertices.back();
	}
      }
      // if (mFaceVertices[in]) {
      // 	face_vertex_ptr fv = mFaceVertices[in];
      // 	mFaceVertices[in] = NULL;
      // 	mRecycle.push_back(in);
      // }
    }

    // void pack(){
    //   if (mRecycle.size() > 0) {
    // 	vector<face_vertex_ptr> tFaceVertices;
    // 	long j = 0;
    // 	for (long i = 0; i < mFaceVertices.size(); i++) {
    // 	  if(mFaceVertices[i]){ 
    // 	    tFaceVertices.push_back(mFaceVertices[i]);
    // 	    tFaceVertices.back()->position_in_vertex() = j;
    // 	    j++;
    // 	  }
    // 	}
    // 	swap(mFaceVertices,tFaceVertices);
    //   mRecycle.clear();
    //   }
    // }        
        
    face_vertex_ptr & front(){
      return mFront;
    }
		
    face_vertex_ptr fbegin(){return this->front();}
    face_vertex_ptr fend()  {return this->front()->vprev();}
        
    list<face_vertex_ptr>& get_face_vertices(){ return mFaceVertices;}
        
    face_vertex_ptr get_insertion_face_vertex(vertex_ptr& that){
      if (mFaceVertices.size() > 1) {
	return find_insertion_face_vertex(that);
      }
      else {
	return mFaceVertices.front();
      }
            
    }
        
    face_vertex_ptr find_insertion_face_vertex(vertex_ptr& that){
            
      face_vertex_ptr itb = mFaceVertices.front();			
      face_vertex_ptr ite = itb->vprev();
            
      coordinate_type this_point = this->coordinate();
      coordinate_type that_point = that->coordinate();
      coordinate_type next_pt, xc,xp;
      face_vertex_ptr out = itb; T d;
      size_t sz = mFaceVertices.size();
      xp = cross(that_point - this_point, ite->next()->coordinate() - this_point);
      next_pt = itb->next()->coordinate();				
      T dt = dist(that_point, next_pt); d = dt;
      while (itb != ite) {
	next_pt = itb->next()->coordinate();				
	xc = cross(that_point - this_point, next_pt - this_point);
	T  s = xc[0]*xp[0] + xc[1]*xp[1] + xc[2]*xp[2];
                
	T dt = dist(that_point, next_pt);
	if (s < 0) {
	  if (dt < d){
	    d = dt;
	    out = itb;
	  }
	}
	xp = xc;
	itb = itb->vnext();
      }			
      return out;
    }
        
    coordinate_array normal_trace(){
      coordinate_array tOut;
      fvl_iterator itb = mFaceVertices.begin();
      fvl_iterator ite = mFaceVertices.end();
      size_t fv_size = mFaceVertices.size();
            
      while (itb != ite) {
	face_vertex_ptr fv = (*itb);
	if (fv) {
	  face_ptr f = fv->face();
	  coordinate_type vert = f->normal();
	  tOut.push_back(vert);
	}
	++itb;
      }
      return tOut;
    }
		
        
    coordinate_array vertex_trace(){
      coordinate_array tOut;
      fvl_iterator itb = mFaceVertices.begin();
      fvl_iterator ite = mFaceVertices.end();
      size_t fv_size = mFaceVertices.size();
            
      while (itb != ite) {
	coordinate_type other_vert = (*itb)->next()->coordinate();
	tOut.push_back(other_vert-mCoordinate);
	//tOut.push_back((*itb)->face()->normal());
	++itb;
      }
      return tOut;
    }
        
    coordinate_array rel_vertex_trace(){
      coordinate_array tOut;
      fvl_iterator itb = mFaceVertices.begin();
      fvl_iterator ite = mFaceVertices.end();
      size_t fv_size = mFaceVertices.size();
            
      while (itb != ite) {
	coordinate_type other_vert = (*itb)->next()->coordinate();
	tOut.push_back(other_vert-mCoordinate);
	//tOut.push_back((*itb)->face()->normal());
	++itb;
      }
      return tOut;
    }
        
    coordinate_array ordered_normal_trace(){
      coordinate_array tOut;
      face_vertex_ptr itb = mFaceVertices.front();
      face_vertex_ptr ite = itb->vprev();
      size_t fv_size = mFaceVertices.size();
            
      while (itb != ite) {
	coordinate_type other_vert = (*itb)->vnext()->face()->normal();
	tOut.push_back(other_vert);
	//tOut.push_back((*itb)->face()->normal());
	itb = itb->vnext();
      }
      return tOut;
    }
        
    void update_normal(){
      coordinate_array vt_list;
      //coordinate_type nNormal;						
      vt_list = this->normal_trace();
      mNormal = calculate_average(vt_list);							
      mNormal.normalize();
      calc_normal = false;
    }
        
    void update_center_of_mass(){
      coordinate_array vt_list = this->vertex_trace();
      vt_list.push_back(mCoordinate);
      mCentMass = calculate_average(vt_list);
    }
        
    void draw(){
            
      this->draw_vertex_point();
      //this->draw_normal();
    }
        
    void draw_vertex_point(){
      glPushMatrix();
      glColor3f(color.r,color.g,color.b);
      glPointSize(5.0f);
      glBegin(GL_POINTS);
      glVertex3f(mCoordinate[0],mCoordinate[1], mCoordinate[2]);
      glEnd();
      glPopMatrix();
    }

    void draw_vertex(){
      glPushMatrix();
      glColor3f(0,0.0,0.0);
      glPointSize(5.0f);
      glBegin(GL_POINTS);
      glVertex3f(mCoordinate[0],mCoordinate[1], mCoordinate[2]);
      glEnd();
      glPopMatrix();
			
      glLineWidth(1.0);
      glBegin(GL_LINE_STRIP);
      T off = 0.05;
      for (int i = 0; i < 3; i++) {
	bool iterating = true;
	face_vertex_ptr fvb = this->fbegin();
	face_vertex_ptr fve = this->fend();				
	while (iterating) {
	  iterating = fvb != fve;
	  coordinate_type cn = fvb->edge()->other(fvb)->vertex()->coordinate();
	  coordinate_type dcn = mCoordinate + off*(cn-mCoordinate);
	  glColor3f(0.0, 0.0, 0.0);
	  glVertex3f(dcn[0],dcn[1],dcn[2]);
	  fvb = fvb->vnext();
	}
	off+=0.05;
      }

      glEnd();
    }
        
    void draw_label(){
      char str[10];
      if (mRecycle.size()>0) {
	this->pack();
      }
      sprintf(str,"%d",position_in_set());
      glColor3f(0., 0., 0.);
      renderModStringX(mCoordinate[0],
		       mCoordinate[1],
		       mCoordinate[2], str);
    }
        
    void draw_normal(){
      if (calc_normal) {
	this->update_normal();
      }
            
      T nx = mCoordinate[0] + mNormal[0]*0.25;
      T ny = mCoordinate[1] + mNormal[1]*0.25;
      T nz = mCoordinate[2] + mNormal[2]*0.25;
            
      glColor4f((1 + mNormal[0])/2, (1 + mNormal[1])/2, (1+mNormal[2])/2, 0.5);
      glLineWidth(0.25);
      glBegin(GL_LINES);
      glVertex3f(mCoordinate[0], mCoordinate[1], mCoordinate[2]);
      glVertex3f(nx, ny, nz);
      glEnd();
    }
        
        
  protected:
    coordinate_type mCoordinate;
    coordinate_type mCentMass;
    coordinate_type mNormal;
    face_vertex_ptr mFront;
    std::list<face_vertex_ptr> mFaceVertices;
    vector<long> mRecycle;
    size_t mID;
    bool calc_normal;
    long mSetPosition;
    
  public:
    int pinned;
    unsigned int flag;
        
    void blank_iteration_block(){
      fvl_iterator itb = mFaceVertices.begin();
      fvl_iterator ite = mFaceVertices.end();
            
      while (itb != ite) {
      }			
    }
        
    bool find_rotation_point_by_angle(fvl_iterator& itb, coordinate_type new_v, T& angle_diff){
            
      coordinate_type vertex_in, vertex_prev, vertex_next, vertex_head;
      T angle_whole, angle_half1, angle_half2;
            
      vertex_in = new_v - mCoordinate;
      bool found_match = false;
      mNormal.normalize();
            
      vertex_head = this->get_head_face_vertex()->next()->get_coordinate() - mCoordinate;
      face_vertex_ptr fvn = (*itb)->next();
      face_vertex_ptr fvp = (*itb)->prev();
      vertex_next = fvn->get_coordinate() - mCoordinate;
      vertex_prev = fvp->get_coordinate() - mCoordinate;
      angle_whole = angle_from_plane( mNormal, vertex_next, vertex_prev);
      angle_half1 = angle_from_plane( mNormal, vertex_in, vertex_next);
      angle_half2 = angle_from_plane( mNormal, vertex_in, vertex_prev);
            
      //			bool next_head	= vertex_next == vertex_head;
            
      bool is_between = (angle_half1+angle_half2) <= (angle_whole + .01);
            
      //bool less_180	= (angle_half1+angle_half2) < PI;
            
      if (is_between) {
	return found_match = true;
      }
            
      else {
	itb++;
      }
            
      return found_match;			
    }
        
        
    void find_rotation_point_by_turn(fvl_iterator& it,coordinate_type new_v){
            
      coordinate_type vec_next, vec_prev, vec_in, vec_cen;
      T det_prev, det_next;
      vec_in = new_v	- mCoordinate;
            
      fvl_iterator& itb = mFaceVertices.begin();
      for (itb = mFaceVertices.begin(); itb!=mFaceVertices.end(); itb++) {
	face_vertex_ptr cur = (*itb);
	vec_next = cur->next()->coordinate() - mCoordinate;
	vec_prev = cur->prev()->coordinate() - mCoordinate;
	//vec_prev2 = (*itb)->prev()->coordinate() - mCoordinate;					
	det_next = determinant(vec_prev, vec_in, vec_next);										
	if (det_next*det_prev< 0) {
	  //if the determinant changes direction then we want to 
	  //insert between the two vectors
	  it = itb;
	}
      }			
    }
        
    bool find_rotation_point_by_pointer(fvl_iterator& it,face_vertex_ptr& new_fv){
      bool inserted = false;			
      vertex_ptr nxt = new_fv->next()->vertex();
      vertex_ptr cur;
      fvl_iterator& itb = mFaceVertices.begin();
      for (itb; itb !=mFaceVertices.end(); itb++) {
	cur = (*itb)->next()->vertex();
	if (cur->ID() == nxt->ID()) {
	  it = itb;
	  inserted = true;
	}
      }
      return inserted;
    }
        
    bool pointer_insert(face_vertex_ptr& new_fv){
      fvl_iterator itb;			
      bool inserted = find_rotation_point_by_turn(itb,new_fv);						
      mFaceVertices.insert(itb,new_fv);
      return inserted;
    }
        
    bool radial_insert(face_vertex_ptr& new_fv){
      fvl_iterator itb;			
      bool inserted = find_rotation_point_by_angle(itb,new_fv->coordinate());						
      mFaceVertices.insert(itb,new_fv);
      return inserted;
    }
  };
}
#endif
