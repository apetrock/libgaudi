/*
 *  m2Face.hpp
 *  Phase Vocoder
 *
 *  Created by John Delaney on 1/8/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef __TWOMANIFOLDFACE__
#define __TWOMANIFOLDFACE__
#include "m2Common.hpp"
#include "m2Control.hpp"
#include "m2FaceVertex.hpp"
#include "m2Vertex.hpp"
#include "m2Edge.hpp"
#include "vec_addendum.h"

//#include "al_Mesh.hpp"

namespace m2 {
  template <typename SPACE>
  class face{
		
    //typedef macro
    M2_TYPEDEFS
		
    public:
    face(){
      m2::ID& manager = m2::ID::get_instance();      
      mID = manager.new_face_id();
      flag = 0;
      color.r = 0.4 + randd(0.1);
      color.g = 0.4 + randd(0.1);
      color.b = 0.4 + randd(0.1);
      color.a = 1.0;
    }
		
    face(vertex_ref pnt){
      //this makes a the face in a hypersphere from one point
      m2::ID& manager = m2::ID::get_instance();
      mID = manager.new_face_id();
			
      face_vertex_ptr nfv = new face_vertex_type();			
      mSize = 1;
			
      nfv->vertex() = &pnt;
      nfv->face() = this;
			
      nfv->prev() = nfv;
      nfv->next() = nfv;
			
      fHead = nfv;
      this->renumber_vertex_IDs();
			
      flag = 0;
			
      color.r = 0.8 + randd(0.2);
      color.g = 0.0 + randd(0.1);
      color.b = 0.0 + randd(0.1);
      color.a = 1.0;
    }
		
		
    ~face(){
      //std::cout << "erasing face: " << ID() << std::endl;
    }
		
    size_t& size()		 { return mSize;}
    size_t  size() const { return mSize;}
		
    long& position_in_set(){
      return mSetPosition;
    }
		
    face_vertex_ptr& fbegin(){
      return fHead;
    }
		
    face_vertex_ptr& fend(){
      return fHead->prev();
    }
		
    face_ref operator=(const face_ref rhs){
      face_ptr out = new face_type();
      m2::ID& manager = m2::ID::get_instance();
      out->mID = manager.new_face_id();
      out->mNormal = rhs.mNormal;
      out->mCenter = rhs.mCenter;
      out->mVertices = rhs.mVertices;
			
      return *out;
    }
		
    coordinate_array flagged_coordinate_trace(unsigned long flag_num){
      face_vertex_ptr itb = fHead;
      face_vertex_ptr ite = fHead->prev();
      coordinate_array array_out;
			
      bool at_head = false;
      while (!at_head) {
	if (itb == ite) {
	  at_head = true;
	}
	long cflag = itb->vertex()->flag;
	if (cflag == flag_num) {
	  coordinate_type cd= itb->vertex()->coordinate();
	  array_out.push_back(cd);
	  cflag = itb->vertex()->flag;
	}
	itb = itb->next();
      }
      return array_out;
    }				

    face_vertex_list vertex_trace(){
      face_vertex_ptr itb = fHead;
      face_vertex_ptr ite = fHead->prev();
      face_vertex_list array_out;
			
      bool at_head = false;
      while (!at_head) {
	if (itb == ite) {
	  at_head = true;
	}
	array_out.push_back(itb);
	itb = itb->next();
      }
      return array_out;
    }
		
    face_vertex_list flagged_vertex_trace(unsigned long flag_num){
      face_vertex_ptr itb = fHead;
      face_vertex_ptr ite = fHead->prev();
      face_vertex_list array_out;
			
      bool at_head = false;
      while (!at_head) {
	if (itb == ite) {
	  at_head = true;
	}
	if (itb->vertex()->flag == flag_num) {
	  array_out.push_back(itb);
	  //itb->vertex()->flag = 0;
	}
	itb = itb->next();
      }
      return array_out;
    }
		
    edge_list flagged_edge_trace(unsigned long flag_num){
      face_vertex_ptr itb = fbegin();
      face_vertex_ptr ite = fend();
      edge_list array_out;
			
      bool at_head = false;
      while (!at_head) {
	if (itb == ite) {
	  at_head = true;
	}
	if (itb->edge()->flag == flag_num) {
	  array_out.push_back(itb->edge());
	  //itb->vertex()->flag = 0;
	}
	itb = itb->next();
      }
      return array_out;
    }
		
    coordinate_array coordinate_trace(){
      face_vertex_ptr itb = fHead;
      face_vertex_ptr ite = fHead->prev();
      coordinate_array array_out;
			
      bool at_head = false;
      while (!at_head) {
	if (itb == ite) {
	  at_head = true;
	}
	coordinate_type cd = itb->coordinate();
	array_out.push_back(cd);
	itb = itb->next();
      }
      return array_out;
    }
		
    edge_list edge_trace(){
      face_vertex_ptr itb = fHead;
      face_vertex_ptr ite = fHead->prev();
      edge_list list_out;
			
      bool at_head = false;
      while (!at_head) {
	if (itb == ite) {
	  at_head = true;
	}
	edge_ptr ep = itb->edge();
	list_out.push_back(ep);
	itb = itb->next();
      }
      return list_out;
    }
		
    void update_all(){
      this->update_vertex_faces();
      this->renumber_vertex_IDs();
      this->update_normal();
      this->update_center();
      //this->update_vertex_normals();			
      mArea = this->calc_area();
    }
		
    void renumber_vertex_IDs(){
      face_vertex_ptr itb = fbegin();
      face_vertex_ptr ite = fend();
      size_t new_id = 0;
      bool at_head = false;
      bool iterating = true;
      while (iterating) {
	iterating = itb != ite;
	itb->position_in_face() = new_id;
	itb = itb->next(); new_id++;
      }
      this->mSize = new_id;
    }
		
    void update_vertex_faces(){
      face_vertex_ptr itb = fbegin();
      face_vertex_ptr ite = fend();
      int i = 0;
      bool at_head = false;
      while (!at_head) {
	at_head = itb == ite;
	itb->face() = this;
	itb = itb->next(); i++;
      }
    }        

    void update_vertex_normals(){
      face_vertex_ptr itb = fHead;
      face_vertex_ptr ite = fHead->prev();
			
      bool at_head = false;
      while (!at_head) {
	if (itb == ite) {
	  at_head = true;
	}
	itb->vertex()->update_normal();
	itb = itb->next();
      }
    }
        
    coordinate_type calculate_center(){
      face_vertex_ptr itb = fHead;
      face_vertex_ptr ite = fHead->prev();
      coordinate_type out;
      bool at_head = false;
      while (!at_head) {
	at_head = itb == ite;
	out += itb->coordinate();
	itb = itb->next();
      }
      out /= (T)this->size();
      return out;
    }

		
    void delete_covertices(){
      face_vertex_ptr itb = fHead->next();
      face_vertex_ptr ite = fHead;
      face_vertex_ptr nxt;
			
      bool at_head = false;
      while (!at_head) {
	at_head = itb==ite;
	nxt = itb->next();
	if (itb->vertex_ID() == nxt->vertex_ID()) {
	  if (itb->next() == fHead) {
	    fHead = itb;
	    ite = fHead;
	  }
	  itb->delete_next();
	}
	itb = itb->next();
      }
      this->renumber_vertex_IDs();
    }
		
    T calc_area(){
      face_vertex_ptr it0 = fHead;
      face_vertex_ptr it1 = it0->next();
      face_vertex_ptr it2 = it1->next();
      face_vertex_ptr ite = fHead->prev();

      coordinate_type c0 = it0->coordinate();
      T out;
      bool iterating = true;
      while (iterating) {
	iterating = it2 != ite;
	coordinate_type c1 = it1->coordinate();
	coordinate_type c2 = it2->coordinate();                
	coordinate_type n = cross(c1-c0,c2-c0);
	out += n.mag()*0.5;                
	it1 = it1->next();
	it2 = it2->next();
      }
      return out;
    }
        
    void update_normal(){
            
      face_vertex_ptr itb = fHead;
      face_vertex_ptr ite = fHead->prev();
      mNormal = 0;
      bool iterating = true;
      while (iterating) {
                
	iterating = itb != ite;
	coordinate_type& curr = itb->coordinate();
	coordinate_type& next = itb->next()->coordinate();
	mNormal[0] += (curr[1] - next[1])*(curr[2] + next[2]);
	mNormal[1] += (curr[2] - next[2])*(curr[0] + next[0]);
	mNormal[2] += (curr[0] - next[0])*(curr[1] + next[1]);
	itb = itb->next();
      }
      mNormal.normalize();
      //            mNormal = fnormalize(mNormal);
      ncolor.r = 0.5 + mNormal[0]*0.5;
      ncolor.g = 0.5 + mNormal[1]*0.5;
      ncolor.b = 0.5 + mNormal[2]*0.5;
      ncolor.a = 1.0;
      calc_normal = false;
			
			
    }
        
		
    void update_center(){
      mCenter = this->calculate_center();
      calc_center = false;
    }
		
    void draw(){
      this->draw_face();
      //this->draw_normal(0.0);
    }
		
    void draw(T off){
      //			this->draw_normal(off);
      this->draw_face(off);
      //			this->draw_face_vertices(off);		
    }
		
    void draw_vertex_colors(){
            if (calc_normal) {
	this->update_normal();
      }
			
      face_vertex_ptr itb = fHead;
      face_vertex_ptr c0 = fHead;
      face_vertex_ptr c1 = c0->next();
      face_vertex_ptr c2 = c1->next();
      face_vertex_ptr ite = fHead->prev();
			
      coordinate_type& vect0 = c0->coordinate();

      glDisable(GL_BLEND);
      glEnable(GL_DEPTH_TEST);					// Enables Depth Testing 
      glEnable(GL_LIGHTING);
      glEnable(GL_COLOR_MATERIAL);
      glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
      glColorMaterial(GL_FRONT, GL_SHININESS);
      glShadeModel(GL_SMOOTH);	
			
      glBegin(GL_TRIANGLE_FAN);
      glNormal3f(mNormal[0], mNormal[1], mNormal[2]);
      gl_set_color(c0->vertex()->color);
      glVertex3f(vect0[0],
		 vect0[1],
		 vect0[2]);
      bool at_head;
      while (!at_head) {
	at_head = itb==ite;
				
	coordinate_type& vect1 = c1->coordinate();
	coordinate_type& vect2 = c2->coordinate();				
	bool vNULL = &vect0 != NULL;
	if (vNULL) {
	  gl_set_color(c1->vertex()->color);
	  glVertex3f(vect1[0],
		     vect1[1],
		     vect1[2]);
	  gl_set_color(c2->vertex()->color);
	  glVertex3f(vect2[0],
		     vect2[1],
		     vect2[2]);
	}
	c1 = c2;
	c2 = c2->next();
	itb = itb->next();
      }			
      glEnd();					
      glEnable(GL_BLEND);
      glDisable(GL_LIGHTING);
    }
    
    void draw_face(){			
      if (calc_normal) {
	this->update_normal();
      }
			
      face_vertex_ptr itb = fHead;
      face_vertex_ptr c0 = fHead;
      face_vertex_ptr c1 = c0->next();
      face_vertex_ptr c2 = c1->next();
      face_vertex_ptr ite = fHead->prev();
			
      coordinate_type& vect0 = c0->coordinate();

      glDisable(GL_BLEND);
      glEnable(GL_DEPTH_TEST);					// Enables Depth Testing 
      glEnable(GL_LIGHTING);
      glEnable(GL_COLOR_MATERIAL);
      glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
      glColorMaterial(GL_FRONT, GL_SHININESS);
      glShadeModel(GL_SMOOTH);	
			
      glBegin(GL_TRIANGLE_FAN);
      glNormal3f(mNormal[0], mNormal[1], mNormal[2]);
      gl_set_color(color);

      if(c0->vertex()->pinned)
	glColor3f(1,0,0);

      glVertex3f(vect0[0],
		 vect0[1],
		 vect0[2]);
      bool at_head;
      while (!at_head) {
	at_head = itb==ite;
				
	coordinate_type& vect1 = c1->coordinate();
	coordinate_type& vect2 = c2->coordinate();				
	bool vNULL = &vect0 != NULL;
	if (vNULL) {
	  //glColor3f(randd(1.0), randd(1.0), randd(1.0));
	  if(c1->vertex()->pinned)
	    glColor3f(1,0,0);
	  else
	    gl_set_color(color);
	  glVertex3f(vect1[0],
		     vect1[1],
		     vect1[2]);
	  //glColor3f(randd(1.0), randd(1.0), randd(1.0));
	  if(c2->vertex()->pinned)
	    glColor3f(1,0,0);
	  else
	    gl_set_color(color);

	  glVertex3f(vect2[0],
		     vect2[1],
		     vect2[2]);
	}
	c1 = c2;
	c2 = c2->next();
	itb = itb->next();
      }			
      glEnd();					
      glEnable(GL_BLEND);
      glDisable(GL_LIGHTING);
    }
		
    void draw_face(T off){
			
      if (calc_normal) {
	this->update_normal();
      }
			
      face_vertex_ptr itb = fHead;
      face_vertex_ptr ite = fHead->prev();
			
			
      glBegin(GL_POLYGON);
      gl_set_color(color);			
      glNormal3f(mNormal[0], mNormal[1], mNormal[2]);
      bool at_head = false;
      while (!at_head) {
	at_head = itb==ite;
	coordinate_type& vect = itb->coordinate();
	if (&vect != NULL) {
	  glVertex3f(vect[0] + off*mNormal[0],
		     vect[1] + off*mNormal[1],
		     vect[2] + off*mNormal[2]);
	}
	itb = itb->next();
      }
			
      glEnd();					
    }	
		
    void draw_face_vertices(T off){
      face_vertex_ptr itb = fHead;
      face_vertex_ptr ite = fHead->prev();
			
      bool at_head = false;
      while (!at_head) {
	if (itb == ite) {
	  at_head = true;
	}
	itb->draw(off);
	itb = itb->next();
      }
			
      glEnd();					
    }
		
    void draw_normal(T off){
      if (calc_normal) {
	this->update_normal();
      }
			
      if (calc_center) {
	this->update_center();
      }
			
      T norm_dist = 0.001;
      T n0x = mCenter[0] + mNormal[0]*off;
      T n0y = mCenter[1] + mNormal[1]*off;
      T n0z = mCenter[2] + mNormal[2]*off;
      T n1x = n0x + mNormal[0]*norm_dist;
      T n1y = n0y + mNormal[1]*norm_dist;
      T n1z = n0z + mNormal[2]*norm_dist;
			
      gl_set_color(ncolor);
      glLineWidth(0.25);
      glBegin(GL_LINES);
      glVertex3f(n0x, n0y, n0z);
      glVertex3f(n1x, n1y, n1z);
      glEnd();
    }
		
    void print(){
      cout << "face " << mID << ": number of vertices: " << this->size() << endl;
      cout << "vertices: " << endl;
			
      face_vertex_ptr itb = fHead;
      face_vertex_ptr ite = fHead->prev();
      size_t new_id = 0;
      bool at_head = false;
      while (!at_head) {
	if (itb == ite) {
	  at_head = true;
	}
	cout << itb->vertex_ID() << ", ";
	itb = itb->next();
      }
      cout << endl;
			
    }
    T area(){
      this->calc_area();
      return mArea;
    }					
    coordinate_type & normal(){return mNormal;}			
    coordinate_type & center(){return mCenter;}
		
    size_t ID() const	{return this->mID;}
    //		size_t	size()		const	{return this->fHead->prev()->face_ID() + 1;}

    T& x() {return mCenter[0];}
    T& y() {return mCenter[1];}
    T& z() {return mCenter[2];}
		
    T  x() const {return mCenter[0];}
    T  y() const {return mCenter[1];}
    T  z() const {return mCenter[2];}
		
  protected:
    //face_vertex_array	mVertices;
		
    face_vertex_ptr		fHead;
    coordinate_type		mCenter;
    coordinate_type		mNormal;
    T					mArea;
    size_t				mID;
    size_t				mSize;
		
    long mSetPosition;
		
    bool calc_normal;
    bool calc_center;
		
  public:
    colorRGB		color;
    colorRGB		ncolor;
    unsigned int flag;
  };	
}
#endif
