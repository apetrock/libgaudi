/*
 *  m2Edge.hpp
 *  Phase Vocoder
 *
 *  Created by John Delaney on 1/8/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef __TWOMANIFOLDEDGE__
#define __TWOMANIFOLDEDGE__
#include "m2Common.hpp"
#include "m2Control.hpp"
#include "m2Face.hpp"
#include "m2FaceVertex.hpp"
#include "m2Vertex.hpp"


namespace m2 {
  template <typename SPACE>
  class edge{
		
    M2_TYPEDEFS		
		
    public:
    edge(){
      fv1 = NULL;
      fv2 = NULL;
      m2::ID& manager = m2::ID::get_instance();
      mID = manager.new_edge_id();
      flag = 0;
    };		
		
    edge(face_vertex_ref in1, face_vertex_ref in2){
      fv1 = &in1;
      fv2 = &in2;
      in1.set_edge(*this);
      in2.set_edge(*this);
      flag = 0;
    }	
		
    edge_ref operator=(edge_ref rhs){
      edge_ptr out = new edge_type(*rhs.fv1,*rhs.fv2);
      m2::ID& manager = m2::ID::get_instance();
      out->mID = manager.new_edge_id();	
      return *out;
    }
		
		
    ~edge(){
    };
		
    long&  position_in_set(){
      return mSetPosition;
    }
    
    size_t	ID()		const	{return this->mID;}
    face_vertex_ptr	&	vertex_1(){	return fv1;}
    face_vertex_ptr	&	vertex_2(){	return fv2;}	
    face_vertex_ptr& v1(){		
      return fv1;
    }	
		
    face_vertex_ptr& v2(){
      return fv2;
    }
		
    face_vertex_ptr	&	other(face_vertex_ptr cv){
      if (cv == fv1) {
	return fv2;
      }
      else return fv1;
    }
		
    face_vertex_ptr	&	other(face_ptr cv){
      if (cv->face() == fv1->face()) {
	return fv2;
      }
      else return fv1;
    }
		
    face_vertex_ptr	&	other(vertex_ptr cv){
      if (cv == fv1->vertex()) {
	return fv2;
      }
      else return fv1;
    }
		
    face_vertex_ptr	&	return_this(vertex_ptr cv){
      if (cv != fv1->vertex()) {
	return fv2;
      }
      else return fv1;
    }
		
    face_vertex_ptr	&	this_fv(face_vertex_ptr cv){
      if (cv != fv1) {
	return fv2;
      }
      else return fv1;
    }
		
    void	set_this(face_vertex_ptr this_vertex, face_vertex_ptr new_vertex){
      if (this_vertex == fv1) fv1 = new_vertex;
      else					fv2 = new_vertex;
      new_vertex->edge() = this;
    }
		
    void	set_other(face_vertex_ptr this_vertex, face_vertex_ptr new_vertex){
      if (this_vertex == fv1) fv2 = new_vertex;
      else					fv1 = new_vertex;
      new_vertex->edge() = this;
    }
		
    void	set_other(face_ptr cv, face_vertex_ptr ov){
      if (cv->face() == fv1->face()) {
	return fv2 = ov;
      }
      else fv1 = ov;
    }
						
    void update_vertex(face_vertex_ptr old_, face_vertex_ptr new_){
      if (fv1 == old_) {
	fv1 = new_;
      }
      else {
	fv2 = new_;
      }
    }
		
    face_ptr coface(face_vertex_ptr this_vert){
      if (this_vert == fv1) {
	return fv2->face();
      }
      else {
	return fv1->face();
      }
    }
		
    void update_face_vertices(){
      fv1->edge() = this;
      fv2->edge() = this;
    }
        
    int vnum(face_vertex_ptr v){            
      if(v == fv1) return 1;
      else if(v == fv2) return 2;
      else return 0;
    }
        
    T dist(){
      coordinate_type c0 = this->v1()->coordinate();
      coordinate_type c1 = this->v2()->coordinate();	
      return (c0 - c1).mag();
    }

    T length(){
      coordinate_type c0 = this->v1()->coordinate();
      coordinate_type c1 = this->v2()->coordinate();	
      return (c0 - c1).mag();
    }
		
    void draw(){
      this->draw(0.0);
    }
		
    void draw(T off){
      if (fv1 && fv2) {
	coordinate_type& vect1 = fv1->coordinate();
	coordinate_type& vect2 = fv2->coordinate();
                
	face_ptr f1 = fv1->face();
	face_ptr f2 = fv2->face();
                
	//if(f1)f1->update_normal();
	//if(f2)f2->update_normal();			
                
	//               glPushMatrix();
	glLineWidth(0.50f);
	glBegin(GL_LINES);
	coordinate_type& n1 = f1->normal();
	coordinate_type& n2 = f2->normal();
	gl_set_color_opaque(f1->color);
	glVertex3f(vect1[0] + off*n1[0],
		   vect1[1] + off*n1[1], 
		   vect1[2] + off*n1[2]);			
                
	gl_set_color_opaque(f2->color);
	glVertex3f(vect2[0] + off*n2[0],
		   vect2[1] + off*n2[1],
		   vect2[2] + off*n2[2]);
	glEnd();

	//                glPopMatrix();

      }
    }
    face_vertex_ptr fv1;
    face_vertex_ptr fv2;		
  protected:

    size_t mID;
    long mSetPosition;

		
  public:
    unsigned int flag;
		
	
  };

}
#endif
