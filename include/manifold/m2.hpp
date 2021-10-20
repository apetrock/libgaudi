/*
 *  m2Common.h
 *  Phase Vocoder
 *
 *  Created by John Delaney on 1/8/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __TWOMANIFOLDCOMMON__
#define __TWOMANIFOLDCOMMON__

#include <iostream>
#include <list>
#include <vector>
#include <algorithm>
#include <math.h>
#include <assert.h>
#include <stack>

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
//#include <GLUT/glut.h>
#else
#ifdef _WIN32
  #include <windows.h>
#endif
//#include <GL/gl.h>
//#include <GL/glu.h>
//#include <GL/glew.h>
//#include <GL/glut.h>
#endif

#include <memory>

#include "manifold_singleton.h"
#include "geometry_types.hpp"
#include "vec_addendum.h"


// this typedef list is ugly but useful!!!
#define M2_TYPEDEFS\
  typedef typename SPACE::coordinate_type  coordinate_type;\
  typedef typename SPACE::line_type       line_type;\
  typedef typename SPACE::triangle_type    triangle_type;\  
  typedef typename SPACE::swept_point_type       swept_point_type;\  
  typedef typename SPACE::swept_triangle_type    swept_triangle_type;\
  typedef typename SPACE::box_type         box_type;\
  typedef typename SPACE::quat          quat;\
  typedef typename SPACE::mat3          mat3;\
  typedef typename SPACE::mat4          mat4;\
  typedef typename SPACE::double_type	 	   T;\
  typedef	m2::face<SPACE>		    face_type;\
  typedef	m2::edge<SPACE>             edge_type;\
  typedef	m2::vertex<SPACE>	    vertex_type;\
  typedef	m2::face_vertex<SPACE>      face_vertex_type;\
  typedef	m2::control<SPACE>	    control_type;\
  typedef       m2::face<SPACE>*	    face_ptr;\
  typedef       m2::edge<SPACE>*	    edge_ptr;\
  typedef       m2::vertex<SPACE>*      vertex_ptr;    \
  typedef       m2::face_vertex<SPACE>* face_vertex_ptr; \
  typedef       m2::control<SPACE>*     control_ptr;   \
							\
  typedef	m2::face<SPACE>&	    face_ref;	                \
  typedef	m2::edge<SPACE>&	    edge_ref;	                \
  typedef	m2::vertex<SPACE>&          vertex_ref;			\
  typedef	m2::face_vertex<SPACE>&     face_vertex_ref;		\
  typedef	m2::control<SPACE>&         control_ref;		\
									\
  typedef vector<coordinate_type >		coordinate_array;	\
  typedef	vector<face_ptr>		face_array;             \
  typedef	vector<edge_ptr>		edge_array;             \
  typedef	vector<vertex_ptr>		vertex_array;           \
  typedef	vector<face_vertex_ptr>		face_vertex_array;      \
									\
  typedef list<coordinate_type >		coordinate_list;        \
  typedef list<face_ptr>			face_list;              \
  typedef list<edge_ptr>			edge_list;              \
  typedef list<vertex_ptr>			vertex_list;	        \
  typedef list<face_vertex_ptr>			face_vertex_list;	\
									\
  typedef typename coordinate_array::iterator		ca_iterator;	\
  typedef	typename vector<face_ptr>::iterator	fa_iterator;    \
  typedef	typename vector<edge_ptr>::iterator	ea_iterator;    \
  typedef	typename vector<vertex_ptr>::iterator	va_iterator;	\
  typedef typename vector<face_vertex_ptr>::iterator	fva_iterator;	\
									\
  typedef typename list<face_ptr>::iterator			fl_iterator; \
  typedef typename list<edge_ptr>::iterator			el_iterator; \
  typedef typename list<vertex_ptr>::iterator			vl_iterator; \
  typedef typename list<face_vertex_ptr>::iterator	fvl_iterator;		

using namespace std;

namespace m2 {


  template <typename SPACE>	
  class face_vertex;
	
  template <typename SPACE>
  class face;	
	
  template <typename SPACE>
  class edge;
	
  template <typename SPACE>
  class vertex;
	
  template <typename SPACE>
  class control;
	
  template <typename SPACE>
  class construct;
	
  template <typename SPACE>
  class subdivide;
  
  template <typename SPACE>  
  struct bounding_box{
    M2_TYPEDEFS;
  public:
    coordinate_type min;
    coordinate_type max;
    bounding_box(const coordinate_type & mn,
		 const coordinate_type & mx)
      :min(mn), max(mx){};
  };

  struct colorRGB {
    double r,g,b,a;
  };

  /*
  inline void gl_set_color(colorRGB in){
    glColor4f(in.r, in.g, in.b, in.a);
  }
  inline void gl_set_color_opaque(colorRGB in){
    glColor4f(in.r, in.g, in.b, 1.0);
  }	
  */
  /*
  template <typename T>
  ostream &operator<<( ostream &out, const al::Vec<3,T>& in ) {
    out << " " << in[0] << " " <<  in[1] << " " <<  in[2] << " ";
    return out;
  }
  */
	
  template <typename T>
  inline	T randd(T range_) {	
    T output =  (T)(rand()/(T(RAND_MAX)+1.0))*range_;
    return output;
  }
  /*
  template <typename T>	
  inline Eigen::Matrix<T,4,1> cross(Eigen::Matrix<T,4,1> a,Eigen::Matrix<T,4,1> b){
    return Eigen::Matrix<T,4,1>( a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x, 1.0);
  }
  */
  template <typename SPACE>
  class edge{		
    M2_TYPEDEFS;	
  public:
    edge(){
      fv1 = NULL;
      fv2 = NULL;
      mSetPosition = -1;
      m2::ID& manager = m2::ID::get_instance();
      mID = manager.new_edge_id();
      flag = 0;
      delete_flag = 0;
    };		
		
    edge(face_vertex_ref in1, face_vertex_ref in2){
      mSetPosition = -1;
      fv1 = &in1;
      fv2 = &in2;
      in1.set_edge(*this);
      in2.set_edge(*this);
      flag = 0;
      delete_flag = 0;
    }	
		
    edge_ref operator=(edge_ref rhs){
      edge_ptr out = new edge_type(*rhs.fv1,*rhs.fv2);
      m2::ID& manager = m2::ID::get_instance();
      out->mID = manager.new_edge_id();	
      return *out;
    }
		
		
    ~edge(){};
		
    int&  position_in_set()       {return mSetPosition;}
    int   position_in_set() const { return mSetPosition;}

    int    group() const {return mGroup;}
    int &  group()       {return mGroup;}

    size_t ID() const {return this->mID;}

    face_vertex_ptr& v1()       { return fv1;}	
    face_vertex_ptr  v1() const { return fv1;}	
    face_vertex_ptr& v2()       {return fv2;}
    face_vertex_ptr  v2() const {return fv2;}
		
    face_vertex_ptr other(const face_vertex_ptr & cv) const {
      if (cv == fv1) {
	return fv2;
      }
      else return fv1;
    }

#if 0		
    face_vertex_ptr  other(face_ptr cv){
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
#endif
		
    face_vertex_ptr return_this(vertex_ptr cv) const {
      if (cv != fv1->vertex()) {
	return fv2;
      }
      else return fv1;
    }
		
    face_vertex_ptr this_fv(face_vertex_ptr cv) const {
      if (cv != fv1) {
	return fv2;
      }
      else return fv1;
    }

    void set(face_vertex_ptr nfv1, face_vertex_ptr nfv2){
      fv1 = nfv1;
      fv2 = nfv2;
      nfv1->edge() = this;
      nfv2->edge() = this;
    }
    
    void set_this(face_vertex_ptr this_vertex, face_vertex_ptr new_vertex){
      if (this_vertex == fv1) fv1 = new_vertex;
      else					fv2 = new_vertex;
      new_vertex->edge() = this;
    }
		
    void set_other(face_vertex_ptr this_vertex, face_vertex_ptr new_vertex){
      if (this_vertex == fv1) fv2 = new_vertex;
      else					fv1 = new_vertex;
      new_vertex->edge() = this;
    }
		
    void set_other(face_ptr cv, face_vertex_ptr ov){
      if (cv->face() == fv1->face()) {
	return fv2 = ov;
      }
      else fv1 = ov;
    }
    
    coordinate_type normal(){
      coordinate_type n1 = fv1->face()->normal();
      coordinate_type n2 = fv2->face()->normal();
      coordinate_type na = n1+n2;
      na.normalize();
      return na;
    }
						
    void update_vertex(face_vertex_ptr old_, face_vertex_ptr new_){
      if (fv1 == old_) {
	fv1 = new_;
      }
      else {
	fv2 = new_;
      }
    }
		
    face_ptr coface(face_vertex_ptr this_vert) const {
      if (this_vert == fv1) {
	return fv2->face();
      }
      else {
	return fv1->face();
      }
    }
		
    void verify(){
      assert(fv1->edge() == this);
      assert(fv2->edge() == this);
    };

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
      return (c0 - c1).norm();
    }

    T length(){
      coordinate_type c0 = this->v1()->coordinate();
      coordinate_type c1 = this->v2()->coordinate();	
      return (c0 - c1).norm();
    }
		
    void draw(){
      this->draw(0.0);
    }
    /*	
    void draw(T off){
      if (fv1 && fv2) {
	coordinate_type& vect1 = fv1->coordinate();
	coordinate_type& vect2 = fv2->coordinate();
                
	face_ptr f1 = fv1->face();
	face_ptr f2 = fv2->face();
                
	//if(f1)f1->update_normal();
	//if(f2)f2->update_normal();			
                
	//               glPushMatrix();
	glLineWidth(4.0f);
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
*/
    face_vertex_ptr fv1;
    face_vertex_ptr fv2;		
  protected:

    size_t mID;
    int mSetPosition;
    int mGroup;
    		
  public:
    unsigned int idata;
    unsigned int flag;
    unsigned int delete_flag;
    
  };

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
      color.a = 0.75;
      fHead = NULL;
      mArea = 0.0;
      mSetPosition = -1;
      mSize = 0;
      data3 = 1.0;
      mNormal = coordinate_type(0,0,0,1);
    }
		
    face(vertex_ref pnt){
      //this makes a the face in a hypersphere from one point
      m2::ID& manager = m2::ID::get_instance();
      mID = manager.new_face_id();
			
      face_vertex_ptr nfv = new face_vertex_type();			
      mSize = 1;
      mArea = 0.0;
			
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
      mSetPosition = -1;
      
      mNormal = coordinate_type(0,0,0,1);
      data3 = 1.0;
    }
		
		
    ~face(){
      //std::cout << "erasing face: " << ID() << std::endl;
    }
		
    int & size()       { return mSize;}
    int   size() const { return mSize;}
		
    int& position_in_set()       {return mSetPosition;}
    int  position_in_set() const {return mSetPosition;}

    void setHead(face_vertex_ptr head){
      fHead = head;
    }

    face_vertex_ptr& fbegin(){
      return fHead;
    }

    face_vertex_ptr fbegin() const {
      return fHead;
    }
		
    face_vertex_ptr fend(){
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

    
    coordinate_type calc_center(){
      face_vertex_ptr itb = fHead;
      face_vertex_ptr ite = fHead->prev();
      coordinate_type cen = coordinate_type(0.0,0.0,0.0,0.0); 
      float n = 0.0;
      bool at_head = false;
      
      while (!at_head) {
	at_head = itb == ite;
	cen += itb->coordinate();
        
        n += 1.0;
	itb = itb->next();
      }
      cen *= 1.0/n;
      return cen;
    }
    
    bounding_box<SPACE> calc_bbox(){
      face_vertex_ptr itb = fHead;
      face_vertex_ptr ite = fHead->prev();
      
      coordinate_type 
	min = this->calc_center(),
	max = min;
      
      bool at_head = false;
      while (!at_head) {
	at_head = itb == ite;
	
	coordinate_type p = itb->coordinate();
	for(int j = 0; j < 3; j++){
	  max[j] = max[j] > p[j] ? max[j] : p[j];
	  min[j] = min[j] < p[j] ? min[j] : p[j];
	}
	itb = itb->next();
      }
      return bounding_box<SPACE>(min,max);
    }
		
    coordinate_array flagged_coordinate_trace(unsigned int flag_num){
      face_vertex_ptr itb = fHead;
      face_vertex_ptr ite = fHead->prev();
      coordinate_array array_out;
			
      bool at_head = false;
      while (!at_head) {
	if (itb == ite) {
	  at_head = true;
	}
	int cflag = itb->vertex()->flag;
	if (cflag == flag_num) {
	  coordinate_type cd= itb->vertex()->coordinate();
	  array_out.push_back(cd);
	  cflag = itb->vertex()->flag;
	}
	itb = itb->next();
      }
      return array_out;
    }				

    face_vertex_ptr get_corner_on_vertex(vertex_ptr v){
      face_vertex_ptr itb = fHead;
      face_vertex_ptr ite = fHead->prev();

      face_vertex_ptr ito = fHead;
      bool at_head = false;
      while (!at_head) {
	if (itb == ite) {
	  at_head = true;
	}
	if (itb->vertex() == v) ito = itb;
	itb = itb->next();
      }
      return ito;
    }

    face_vertex_array face_vertex_trace(){
      face_vertex_ptr itb = fHead;
      face_vertex_ptr ite = fHead->prev();
      face_vertex_array array_out;
			
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
		
    face_vertex_list flagged_vertex_trace(unsigned int flag_num){
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
		
    edge_list flagged_edge_trace(unsigned int flag_num){
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

    T thinness(){
    face_vertex_ptr fv0 = this->fbegin();
    face_vertex_ptr fv1 = fv0->next();
      face_vertex_ptr fv2 = fv1->next();
      T d0 = fv0->edge()->length();
      T d1 = fv1->edge()->length();
      T d2 = fv2->edge()->length();
      T tEps = 0.95;
      T thin0 = d0/(d1 + d2);
      T thin1 = d1/(d0 + d2);
      T thin2 = d2/(d0 + d1);
      T thin = thin0 > thin1 ? thin0 : thin1;
      thin = thin > thin2 ? thin : thin2;
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
      int i = 0;
      while (iterating && i < 200) {
	iterating = itb != ite;
	itb->position_in_face() = new_id;
	itb = itb->next(); new_id++; i++;      
      }
      this->mSize = new_id;
    }
		
    void update_vertex_faces(){
      face_vertex_ptr itb = fbegin();
      face_vertex_ptr ite = fend();
      int i = 0;
      bool at_head = false;
      while (!at_head) {
	if(itb->vertex() == NULL) std::cout << "NULL VERTEX: " << std::endl;
	if(itb->next() == NULL) std::cout << "NULL NEXT: " << std::endl;
	if(itb->prev() == NULL) std::cout << "NULL PREV: " << std::endl;
	at_head = itb == ite;
	itb->face() = this;
	itb = itb->next(); i++;
      }
    }        

    void update_vertex_normals(){
      face_vertex_ptr itb = fHead;
      face_vertex_ptr ite = fHead->prev();
      int i = 0;
      bool at_head = false;
      while (!at_head && i < 200) {
	if (itb == ite) {
	  at_head = true;
	}
	itb->vertex()->update_normal();
	itb = itb->next(); i++; 
      }
    }
    /*   
    coordinate_type calculate_center(){
      face_vertex_ptr itb = fHead;
      face_vertex_ptr ite = fHead->prev();
      coordinate_type out(0,0,0,0);
      int i = 0;
      bool at_head = false;
      while (!at_head && i < 200) {
	at_head = itb == ite;
	out += itb->coordinate();
	itb = itb->next(); i++;
      }
      out /= (T)this->size();
      return out;
    }

    */		
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
      T out = 0;
      bool iterating = true;
      while (iterating) {
	iterating = it2 != ite;
	coordinate_type c1 = it1->coordinate();
	coordinate_type c2 = it2->coordinate();  
	coordinate_type c10 = c1-c0;
	coordinate_type c20 = c2-c0;
	coordinate_type n = cross(c10,c20);
	out += n.norm()*0.5;                
	it1 = it1->next();
	it2 = it2->next();
      }
      return out;
    }
        
    void update_normal(){
            
      face_vertex_ptr itb = fHead;
      face_vertex_ptr ite = fHead->prev();
      mNormal = coordinate_type(0,0,0,0);
      bool iterating = true;

      while (iterating) {
	iterating = itb != ite;
	coordinate_type curr = itb->coordinate();
	coordinate_type next = itb->next()->coordinate();
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
      mCenter = this->calc_center();
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
    /*
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
    */
    /*
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
    */
    /*
    void draw_face(T off){
			
      if (calc_normal) {
	this->update_normal();
      }
			
      face_vertex_ptr itb = fHead;
      face_vertex_ptr ite = fHead->prev();
			
			
      glBegin(GL_POLYGON);
      gl_set_color(color);
      glColor4f(color.r,color.g,color.b,data3);
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
    */
    /*
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
    */

    /*
    void draw_normal(T off){
      
      this->update_normal();
      this->update_center();
			
      T norm_dist = 0.01;
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
    */
    
    bool has_vertex(vertex_ptr v){
      face_vertex_ptr itb = fbegin();
      face_vertex_ptr ite = fend();      
      bool iterating = true;
      while (iterating) {
	iterating = itb != ite;
	if(itb->next()->vertex() == v) return true;
	itb = itb->next();
      }
      return false;
    }
    
    void verify(){
      face_vertex_ptr itb = fbegin();
      face_vertex_ptr ite = fend();
      bool iterating = true;
      while (iterating) {
	iterating = itb != ite;
	assert(itb->prev()->next() == itb);
	assert(itb->next()->prev() == itb);
	assert(itb->face() == this);
	itb = itb->next();
      }
    }
		
    void print(){
      cout << "face " << mSetPosition << ": number of vertices: " << this->size() << endl;
      cout << "vertices: " << endl;
			
      face_vertex_ptr itb = fbegin();
      face_vertex_ptr ite = fend();
      size_t new_id = 0;
      bool iterating = true;
      while (iterating) {
	iterating = itb != ite;
	cout << itb->vertex()->position_in_set() << ", ";
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
    int    group() const {return mGroup;}
    int &  group()       {return mGroup;}
  protected:
    //face_vertex_array	mVertices;
    int mGroup;

    face_vertex_ptr fHead;
    coordinate_type mCenter;
    coordinate_type mNormal;
    T		     mArea;
    int mID;
    int mSize;
		
    int mSetPosition;
		
    bool calc_normal;
    //bool calc_center;
		
  public:
    colorRGB	 color;
    colorRGB	 ncolor;
    unsigned int flag;
    coordinate_type data;
    coordinate_type data2;
    T data3;
  };

  template <typename SPACE>
  class face_vertex{
		
    M2_TYPEDEFS
		
    public:
		
    face_vertex(){
      mVertex	= NULL;
      mEdge	= NULL;
      mFace	= NULL;

      m2::ID& manager = m2::ID::get_instance();
      mID = manager.new_face_vertex_id();
      fID = 0;
      flag = 0;
      data = 0;
    }    
		
    face_vertex(const face_vertex_ref rhs){
      this->mEdge	= rhs.mEdge;
      //			mEdge->set_this(&rhs,this);
      this->mFace	= rhs.mFace;
      this->mVertex	= rhs.mVertex;
      this->nxt_face	= rhs.nxt_face;
      this->prv_face	= rhs.prv_face;	
      this->data = rhs.data;
      fID = 0;
			
      m2::ID& manager = m2::ID::get_instance();
      mID = manager.new_face_vertex_id();
			
      flag = 0;
    }
		
    bool operator==(const face_vertex_ref rhs){
      if (mID == rhs.mID) {
	return true;
      }
      else return false;
    }
		
    ~face_vertex(){
      //mVertex->remove_face_vertex(mVertexPosition);
    };
		
    face_vertex_ref operator=(const face_vertex_ref rhs){
      face_vertex_ptr out = new face_vertex_type(rhs);			
      out->mEdge->update_vertex(&rhs,this);
      return *out;
    }
    
    T angle(){
      coordinate_type ci = this->coordinate();
      coordinate_type ca = this->next()->coordinate();
      coordinate_type cb = this->prev()->coordinate();
      coordinate_type cai = ca - ci;
      coordinate_type cbi = cb - ci;
      T maga = sqrt(cai[0]*cai[0] + cai[1]*cai[1] + cai[2]*cai[2]);
      T magb = sqrt(cbi[0]*cbi[0] + cbi[1]*cbi[1] + cbi[2]*cbi[2]);
      T dotab = cai[0]*cbi[0] + cai[1]*cbi[1] + cai[2]*cbi[2];
      return acos(dotab/(maga*magb));
    }

    face_ref coface(){
      return mEdge->return_coface();
    }
    
    face_vertex_ptr coedge(){
      return mEdge->other(this);
    }

    face_ref get_face(){
      return *mFace;
    }
		
    face_vertex_ptr add_next(){
      face_vertex_ptr	out = new face_vertex(*this);
      face_vertex_ptr nxt = this->nxt_face;	
      out->next()  = nxt;
      nxt->prev()  = out;
      out->prev()  = this;
      this->next() = out;			
			
      mFace->size() += 1;

      out->face() = this->face();
      out->face()->fbegin() = out;
      this->vertex()->add_face_vertex(out);	
      this->vertex()->front() = out;
			
      return out;
    }
		
		
    face_vertex_ptr add_prev(){
      face_vertex_ptr out = new face_vertex(*this);
      face_vertex_ptr prv = this->prv_face;
      out->prev()  = prv;
      prv->next()  = out;
      out->next()  = this;
      this->prev() = out;

      //prv->edge()->set_this(this,out);
      mFace->size() += 1;

      out->face() = this->face();	
      out->face()->fbegin() = out;	
      this->vertex()->add_face_vertex(out);	
      this->vertex()->front() = out;
      return out;
    }

    face_vertex_ptr & next()       {return nxt_face;}
    face_vertex_ptr   next() const {return nxt_face;}		
    face_vertex_ptr & prev()       {return prv_face;}
    face_vertex_ptr   prev() const {return prv_face;}
		
    face_vertex_ptr  vnext() {
      if (mEdge == NULL) {
	return NULL;
      }
      else {
	face_vertex_ptr out = mEdge->other(this);
	return out->next();
      }			
    }
		
    face_vertex_ptr vprev() {
      if (this->prev() == NULL || this->prev()->mEdge == NULL) {
	return NULL;
      }
      else {
	face_vertex_ptr out = this->prev();
	if (out->mEdge)	return out->mEdge->other(out);
	else return NULL;
      }
    }
		
    void draw(T off){
      this->draw_vertex(off);
      this->draw_tail(off);
    }
    /*
    void draw_vertex(T off){

      coordinate_type n1 = mFace->normal();
			
      T
	t0x = n1[0]*off + this->x(),
	t0y = n1[1]*off + this->y(),
	t0z = n1[2]*off + this->z();
			
      //			glPushMatrix();
      glPointSize(2.0f);
      glBegin(GL_POINTS);
      glColor4f(0.0,0.0,0.0,0.0);			
      glVertex3f(t0x,t0y,t0z);
      glEnd();
      //			glPopMatrix();
    }
    */		
    /*
    void draw_tail(T off){

      coordinate_type n1 = mFace->normal();			
      T 
	t0x = n1[0]*off + this->x(),
	t0y = n1[1]*off + this->y(),
	t0z = n1[2]*off + this->z(),
			
	t1x = t0x - (this->x() - nxt_face->x())*0.3,
	t1y = t0y - (this->y() - nxt_face->y())*0.3,
	t1z = t0z - (this->z() - nxt_face->z())*0.3;

      //			glPushMatrix();
      glLineWidth(0.5f);			
      glBegin(GL_LINES);			
      glColor4f(0.1,0.1,0.1,0.0);			
      glVertex3f(t0x, t0y, t0z);					
      glVertex3f(t1x, t1y, t1z);						
      glEnd();					
      //			glPopMatrix();
    }
    */		
		
    edge_ptr & edge()		{return	mEdge;}
    edge_ptr   edge() const     {return	mEdge;}
    face_ptr & face()		{return	mFace;}
    face_ptr   face()   const   {return	mFace;}
    face_ptr   coface()	const   {return mEdge->other(this)->face();}
    vertex_ptr & vertex()	{return	mVertex;}
    vertex_ptr   vertex() const {return	mVertex;}

    int	vertex_ID()	const	{return mVertex->ID();}
    int	ID()		const	{return this->mID;}
    int	face_ID()	const	{return this->mFacePosition;}
    int&	face_ID()			{return this->mFacePosition;}

    int& position_in_face()       {return mFacePosition;}
    int  position_in_face() const {return mFacePosition;}
    void	set_edge(edge_ref input)	{mEdge	= &input;};
    void	set_face(face_ref input)	{mFace	= &input;};
    void	set_vertex(vertex_ref input){mVertex	= &input;};
		
    coordinate_type& coordinate(){return mVertex->coordinate();}
		
    T & x()	      {return mVertex->coordinate()[0];}
    T   x()	const {return mVertex->coordinate()[0];}
    T & y()	      {return mVertex->coordinate()[1];}
    T   y()	const {return mVertex->coordinate()[1];}
    T & z()	      {return mVertex->coordinate()[2];}
    T   z()	const {return mVertex->coordinate()[2];}
    int &  group()       {return mGroup;}
    int    group() const {return mGroup;}
    T& operator[](int i)        {return mVertex->coordinate()[i];}
    T  operator[](int i) const  {return mVertex->coordinate()[i];}

  protected:
    int                 mGroup;
    int			mID;
    int			fID;
    
    face_vertex_ptr nxt_face;
    face_vertex_ptr prv_face;

    edge_ptr		mEdge;
    face_ptr		mFace;
    vertex_ptr		mVertex;

    int mFacePosition;
  public:
    unsigned int flag;
    T data;
  };


  template <typename SPACE, typename F> 
  void for_each(m2::vertex<SPACE> * v, F& func){
    M2_TYPEDEFS;
    face_vertex_ptr fvb = v->fbegin();
    face_vertex_ptr fve = v->fend();
    bool iterating = true;
    while(iterating){
      iterating = fvb != fve;
      func(fvb);
      fvb = fvb->vnext();
    }
  }


  template <typename SPACE>
  class vertex{
    M2_TYPEDEFS;
    public:
    int graphColor;
    bool smooth;
    int isDegenerate;
    colorRGB color;

    vertex(){
      int graphColor = -1;
      //mNormal.zero();			
      mNormal = coordinate_type(0.0,0.0,0.0,1.0);
      mCoordinate = coordinate_type(0.0,0.0,0.0,1.0);
      m2::ID& manager = m2::ID::get_instance();
      mID = manager.new_vertex_id();

      flag = 0;
      mSize = 0;
      color.r = 0.4 + randd(0.1);
      color.g = 0.4 + randd(0.1);
      color.b = 0.4 + randd(0.1);
      color.a = 1.0;
      pinned = false;
      mFront = NULL;
      mSetPosition = -1;
      winding = 0;
      data2 = 0;
      isDegenerate = -1;
    }
		
    vertex(T x, T y, T z){
      int graphColor = -1;

      mCoordinate = coordinate_type(0.0,0.0,0.0,1.0);
      mCoordinate[0] = x;
      mCoordinate[1] = y;
      mCoordinate[2] = z;
      mNormal = coordinate_type(0.0,0.0,0.0,1.0);
      
      m2::ID& manager = m2::ID::get_instance();
      mID = manager.new_vertex_id();

      flag = 0;
      mSize = 0;
      color.r = 0.4 + randd(0.1);
      color.g = 0.4 + randd(0.1);
      color.b = 0.4 + randd(0.1);
      color.a = 1.0;
      pinned = false;
      mFront = NULL;
      mSetPosition = -1;
      winding = 0;
      data2 = 0;
      isDegenerate = -1;
    }
		
    vertex(coordinate_type co){
      int graphColor = -1;
      mCoordinate[0] = co[0];
      mCoordinate[1] = co[1];
      mCoordinate[2] = co[2];
      mCoordinate[3] = co[3];
      mNormal.setZero();

      m2::ID& manager = m2::ID::get_instance();
      mID = manager.new_vertex_id();
			
      flag = 0;
      mSize = 0;
      color.r = 0.4 + randd(0.1);
      color.g = 0.4 + randd(0.1);
      color.b = 0.4 + randd(0.1);
      color.a = 1.0;
      pinned = false;
      mFront = NULL;
      mSetPosition = -1;
      winding = 0;
      data2 = 0;
      isDegenerate = -1;
    }		      
		
    void init(){
      int graphColor = -1;
      face_ptr		new_face = new face_type(*this);
      face_vertex_ptr new_fv = new_face->fbegin();
      new_fv->face()	= new_face;
      this->add_face_vertex(new_fv);
      pinned = false;
      data2 = 0;
    }
		
    int& position_in_set(){return mSetPosition;}
		
    T & x()	      {return mCoordinate[0];}
    T   x()	const {return mCoordinate[0];}
    T & y()	      {return mCoordinate[1];}
    T   y()	const {return mCoordinate[1];}
    T & z()	      {return mCoordinate[2];}
    T   z()	const {return mCoordinate[2];}
    T & operator[](int i)       {return mCoordinate[i];}
    T   operator[](int i) const {return mCoordinate[i];}
		
    int & group()       {return mGroup;}
    int   group() const {return mGroup;}

    int & ID()          {return mSetPosition;}
    int   ID()	const	{return mSetPosition;}

    int calc_size(){
      face_vertex_ptr fvb = this->fbegin();
      face_vertex_ptr fve = this->fend();
      bool iterating = true;
      int sz = 0;
      while(iterating){
	iterating = fvb != fve;
	fvb = fvb->vnext();
	sz++;
      }
      return sz;
    }

    int size() {return mSize;}
    coordinate_type normal() const{return mNormal;}
    coordinate_type& coordinate()       {return mCoordinate;}
    coordinate_type  coordinate() const	{return mCoordinate;}
        
    void add_face_vertex(face_vertex_ptr new_fv){
      mFront = new_fv;
      mSize ++;
      // mFaceVertices.push_back(new_fv);
      // fvl_iterator end = mFaceVertices.end();
      // end--;
      // new_fv->position_in_vertex() = end;
      new_fv->vertex() = this;
    }        

    void remove_face_vertex(face_vertex_ptr fv){
      mSize--;
      //mFaceVertices.erase(it);

      // if (mFaceVertices[in]) {
      // 	face_vertex_ptr fv = mFaceVertices[in];
      // 	mFaceVertices[in] = NULL;
      // 	mRecycle.push_back(in);
      // }
    }

    // void pack(){
    //   if (mRecycle.size() > 0) {
    // 	vector<face_vertex_ptr> tFaceVertices;
    // 	int j = 0;
    // 	for (int i = 0; i < mFaceVertices.size(); i++) {
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
        
    face_vertex_ptr & front(){return mFront;}
    face_vertex_ptr fbegin() {return this->front();}
    face_vertex_ptr fend()   {return this->front()->vprev();}
    
    T thinness(){
      if(mSize == 0) return 0.0;
	vertex_ptr v = this;	
	face_vertex_ptr fvb = v->fbegin();
	face_vertex_ptr fve = v->fend();
	bool iterating = true;
	T vthin = 0; int s = 0;
	while(iterating){
	  iterating = fvb != fve;
	  vthin += fvb->face()->thinness();
	  fvb = fvb->vnext(); s++;
	}
	vthin /= (T)s;

	return vthin;
    }

    bool shares_edge_with(vertex_ptr vi){
      if(mSize == 0) return false;
	vertex_ptr v = this;
	face_vertex_ptr fvb = v->fbegin();
	face_vertex_ptr fve = v->fend();
	bool iterating = true;
	bool share = false;
	while(iterating){
	  iterating = fvb != fve;
	  if(fvb->next()->vertex() == vi) share = true;
	  fvb = fvb->vnext();
	}
	return share;
    }

    edge_ptr get_shared_edge(vertex_ptr vi){
      if(mSize == 0) return false;
	vertex_ptr v = this;
	face_vertex_ptr fvb = v->fbegin();
	face_vertex_ptr fve = v->fend();
	bool iterating = true;
	bool share = false;
	while(iterating){
	  iterating = fvb != fve;
	  if(fvb->coedge()->vertex() == vi) return fvb->edge();
	  fvb = fvb->vnext();
	}
	return NULL;
    }

    face_vertex_ptr & get_insertion_face_vertex(vertex_ptr that){
      if (size() > 1) {
	return find_insertion_face_vertex(that);
      }
      else {
	return mFront;
      }       
    }
        
    face_vertex_ptr & find_insertion_face_vertex(vertex_ptr that){
            
      face_vertex_ptr itb = mFront;			
      face_vertex_ptr ite = itb->vprev();
            
      coordinate_type this_point = this->coordinate();
      coordinate_type that_point = that->coordinate();
      coordinate_type next_pt, xc,xp;
      face_vertex_ptr out = itb; T d;
      size_t sz = size();
      coordinate_type dp  = that_point - this_point;
      coordinate_type dnp = ite->next()->coordinate() - this_point;
      xp = cross(dp, dnp);
      next_pt = itb->next()->coordinate();				
      T dt = dist(that_point, next_pt); d = dt;
      while (itb != ite) {
	next_pt = itb->next()->coordinate();				
	//coordinate_type dp  = that_point - this_point;
        coordinate_type dnp = next_pt - this_point;
        xc = cross(dp, dnp);
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
        
        
    // coordinate_array vertex_trace(){
    //   coordinate_array tOut;
    //   fvl_iterator itb = mFaceVertices.begin();
    //   fvl_iterator ite = mFaceVertices.end();
    //   size_t fv_size = mFaceVertices.size();
            
    //   while (itb != ite) {
    // 	coordinate_type other_vert = (*itb)->next()->coordinate();
    // 	tOut.push_back(other_vert-mCoordinate);
    // 	//tOut.push_back((*itb)->face()->normal());
    // 	++itb;
    //   }
    //   return tOut;
    // }
        
    // coordinate_array one_ring_trace(){
    //   coordinate_array tOut;
    //   fvl_iterator itb = mFaceVertices.begin();
    //   fvl_iterator ite = mFaceVertices.end();
    //   size_t fv_size = mFaceVertices.size();
            
    //   while (itb != ite) {
    // 	coordinate_type other_vert = (*itb)->next()->coordinate();
    // 	tOut.push_back(other_vert-mCoordinate);
    // 	//tOut.push_back((*itb)->face()->normal());
    // 	++itb;
    //   }
    //   return tOut;
    // }
        
    // coordinate_array ordered_normal_trace(){
    //   coordinate_array tOut;
    //   face_vertex_ptr itb = mFaceVertices.front();
    //   face_vertex_ptr ite = itb->vprev();
    //   size_t fv_size = mFaceVertices.size();
            
    //   while (itb != ite) {
    // 	coordinate_type other_vert = (*itb)->vnext()->face()->normal();
    // 	tOut.push_back(other_vert);
    // 	//tOut.push_back((*itb)->face()->normal());
    // 	itb = itb->vnext();
    //   }
    //   return tOut;
    // }
        
    void update_normal(){
      if(mSize == 0) return;
      coordinate_type N(0,0,0,0);
      int n = 0;
      face_vertex_ptr itb = mFront;
      face_vertex_ptr ite = mFront->vprev();
      bool iterating = true;

      while(iterating) {
	iterating = itb != ite;
	face_ptr f = itb->face();
	if(f->area() > 1e-12){
	  coordinate_type Ni = f->normal();
	  N += Ni;  n++;
	}
	itb = itb->vnext();
      }

      mNormal = 1.0/(T)n*N;
      mNormal.normalize();
      calc_normal = false;
    }

    void verify(){
      if(mSize == 0) return;
      face_vertex_ptr itb = this->fbegin();
      face_vertex_ptr ite = this->fend();
      bool iterating = true;
      while(iterating) {
	iterating = itb != ite;
	assert(itb->vertex() == this);
	itb = itb->vnext();
      }
    }

    void print(){
      if(mSize == 0) return;
      face_vertex_ptr itb = this->fbegin();
      face_vertex_ptr ite = this->fend();
      bool iterating = true;
      std::cout << " - vertex: " << mSetPosition << ", size: "<< this->mSize << std::endl;
      std::cout << " edges: ";
      while(iterating) {
	iterating = itb != ite;
	//std::cout << itb->next()->vertex()->position_in_set() << ", ";
	std::cout << itb->edge()->position_in_set() << ", ";
	itb = itb->vnext();
      }
      std::cout << std::endl;
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
    /*  
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
    */  
        
  protected:
    coordinate_type mCoordinate;
    coordinate_type mCentMass;
    coordinate_type mNormal;
    face_vertex_ptr mFront;

    size_t mID;
    bool calc_normal;
    int mSetPosition;
    int mSize;
    int mGroup;
  public:
    int pinned;
    unsigned int flag;
    coordinate_type data;
    T data2;
    T winding;
  };


  template <typename SPACE>
  class control {				
    public:
    M2_TYPEDEFS

    
    struct faceDistFromPointSorter {
      bool operator() (face_ptr fi, face_ptr fj) { 
	coordinate_type ci = fi->center();
	coordinate_type cj = fj->center();
	T di = norm2(ci - mPoint);
	T dj = norm2(cj - mPoint);
	return (di > dj);
      }
      coordinate_type mPoint;
    } mFaceSorter;
		
    control(){
      maxGraphColor = 0;
      manual_clean_up = false;
    }

    control(const control_ref rhs){
      maxGraphColor = 0;
      manual_clean_up = false;
      //mFaces.clear();
      //mVertices.clear();
      //mEdges.clear();			
      mFaces.resize(rhs.mFaces.size());
      mVertices.resize(rhs.mVertices.size());          
      mEdges.resize(rhs.mEdges.size());
            
      for(int i = 0; i < rhs.mVertices.size(); i++){
	coordinate_type nc(rhs.mVertices[i]->coordinate());
	mVertices[i] = new vertex_type(nc);
	mVertices[i]->position_in_set() = i;
      }
            
      for(int i = 0; i < rhs.mEdges.size(); i++){
	mEdges[i] = new edge_type();
	mEdges[i]->position_in_set() = i;
      }
            
      for(int i = 0; i < rhs.mFaces.size(); i++){
	if (rhs.mFaces[i]) {
	  face_ptr nf =  new face_type();
	  nf->position_in_set() = i;
                    
	  face_ptr of = rhs.mFaces[i];
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
      mFaceRecycle = rhs.mFaceRecycle;
      mVertexRecycle = rhs.mVertexRecycle;
      mEdgeRecycle = rhs.mEdgeRecycle;
      bool here = true;
    }

    control_ref operator=(const control_ref rhs){
      std::cout << "deep copy equals" << std::endl;
      //deep copy
      maxGraphColor = 0;
      manual_clean_up = false;
      mFaces.clear();
      mVertices.clear();
      mEdges.clear();			
      mFaces.resize(rhs.mFaces.size());
      mVertices.resize(rhs.mVertices.size());          
      mEdges.resize(rhs.mEdges.size());
            
      for(int i = 0; i < rhs.mVertices.size(); i++){
	coordinate_type nc(rhs.mVertices[i]->coordinate());
	mVertices[i] = new vertex_type(nc);
	mVertices[i]->position_in_set() = i;
      }
            
      for(int i = 0; i < rhs.mEdges.size(); i++){
	// std::cout << rhs.mEdges[i]->v1()->face()->position_in_set() << " "
	// 	  << rhs.mEdges[i]->v2()->face()->position_in_set() << std::endl;
	mEdges[i] = new edge_type();
	mEdges[i]->position_in_set() = i;
      }
            
      for(int i = 0; i < rhs.mFaces.size(); i++){
	if (rhs.mFaces[i]) {
	  face_ptr nf =  new face_type();
	  nf->position_in_set() = i;

	  face_ptr of = rhs.mFaces[i];
	  nf->size() = of->size();
	  face_vertex_ptr itb = of->fbegin();
	  face_vertex_ptr ite = of->fend();                    
	  vector<face_vertex_ptr> tmpFaceArray; tmpFaceArray.reserve(3);
	  bool iterating = true;
	  int fs = 0;
	  while (iterating) {
	    iterating = itb != ite;
	    face_vertex_ptr fv1 = new face_vertex_type(); 
	    edge_ptr   oe = itb->edge();
	    edge_ptr   ne = mEdges[oe->position_in_set()];  
	    vertex_ptr nv = mVertices[itb->vertex()->position_in_set()];       

	    //std::cout << i << " " << oe->position_in_set() << std::endl;

	    if(oe->v1() == itb) ne->v1() = fv1;
	    else                ne->v2() = fv1;

	    //if(!ne1->v1()) ne1->v1() = fv1; else ne1->v2() = fv1;
	    fv1->vertex() = nv;
	    nv->add_face_vertex(fv1);

	    fv1->edge() = ne;
	    fv1->face() = nf;
	    
	    tmpFaceArray.push_back(fv1);
	    fs++;
	    itb = itb->next();
	  }
	  for(int j = 0; j < tmpFaceArray.size(); j++){
	    face_vertex_ptr fvi = tmpFaceArray[j];
	    face_vertex_ptr fvn = tmpFaceArray[(j+1)%fs];
	    fvi->next() = fvn; fvn->prev() = fvi;
	  }
	  nf->fbegin() = tmpFaceArray[0];
	  mFaces[i] = nf;
	}

      }
      mFaceRecycle = rhs.mFaceRecycle;
      mVertexRecycle = rhs.mVertexRecycle;
      mEdgeRecycle = rhs.mEdgeRecycle;
      bool here = true;
      this->update_all();
      return *this;
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

    vector<coordinate_type> get_coordinates(){
      vector<coordinate_type> out; out.resize(mVertices.size());
      for(int i = 0; i < mVertices.size(); i++){
	if(!mVertices[i]) continue;
	out[i] = mVertices[i]->coordinate();
      }
      return out;
    }

    vector<coordinate_type> get_normals(){
      vector<coordinate_type> out; out.resize(mVertices.size());
      for(int i = 0; i < mVertices.size(); i++){
	if(!mVertices[i]) continue;
	mVertices[i]->update_normal();
	out[i] = mVertices[i]->normal();
      }
      return out;
    }

    void assign_coordinates(const vector<coordinate_type> & in){
      for(int i = 0; i < mVertices.size(); i++){
	if(!mVertices[i]) continue;
	if(mVertices[i]->pinned) continue;
	mVertices[i]->coordinate() = in[i];
      }
    }

    void assign_offset_coordinates(const vector<coordinate_type> & in, T amt){
      for(int i = 0; i < mVertices.size(); i++){
	if(!mVertices[i]) continue;
	if(mVertices[i]->pinned) continue;
	mVertices[i]->coordinate() += amt*in[i];
      }
    }

    coordinate_type calc_center(){
      coordinate_type avg(0,0,0,0);
      for(int i = 0; i < mVertices.size(); i++){
	if(!mVertices[i]) continue;
	avg += mVertices[i]->coordinate();
      }
      avg /= (T)mVertices.size();
      return avg;
    }

    coordinate_type calc_min(){
      coordinate_type min = this->calc_center();
      for(int i = 0; i < mVertices.size(); i++){
	coordinate_type p = mVertices[i]->coordinate();
	for(int j = 0; j < 3; j++)
	  min[j] = min[j] < p[j] ? min[j] : p[j];
      }
      return min;
    }

    coordinate_type calc_max(){
      coordinate_type max = this->calc_center();
      for(int i = 0; i < mVertices.size(); i++){
	coordinate_type p = mVertices[i]->coordinate();
	for(int j = 0; j < 3; j++)
	  max[j] = max[j] > p[j] ? max[j] : p[j];
      }
      return max;
    }

    bounding_box<SPACE> calc_bbox(){
      
      coordinate_type max = this->calc_center();
      coordinate_type min = max;
     
      for(int i = 0; i < mVertices.size(); i++){
	coordinate_type p = mVertices[i]->coordinate();
	for(int j = 0; j < 3; j++){
	  max[j] = max[j] > p[j] ? max[j] : p[j];
	  min[j] = min[j] < p[j] ? min[j] : p[j];
	}
      }
      
      return bounding_box<SPACE>(min,max);
    }

    void sortFacesByPoint(coordinate_type c){
      this->pack();
      mFaceSorter.mPoint = c;
      std::sort (mFaces.begin(), mFaces.end(), mFaceSorter);
      for(int i = 0; i < mFaces.size(); i++){
	//std::cout << norm2(mFaces[i]->center() - c) << std::endl;
	mFaces[i]->position_in_set() = i;
      }
    }
	
    void merge(control_ref other){
      this->pack();
      other.pack();
      for(int i = 0; i < other.mFaces.size(); i++){
	if (other.mFaces[i]) {
	  this->push_face(other.mFaces[i]);
	}
      }
      for(int i = 0; i < other.mEdges.size(); i++){
	if (other.mEdges[i]) {
	  this->push_edge(other.mEdges[i]);
	}
      }
      for(int i = 0; i < other.mVertices.size(); i++){
	if (other.mVertices[i]) {
	  this->push_vertex(other.mVertices[i]);
	}
      }
      other.mFaces.clear();
      other.mEdges.clear();
      other.mVertices.clear();			
    }
		
    void push_vertex(vertex_ptr in){
      if (mVertexRecycle.size() > 0 && !manual_clean_up){
	int i = mVertexRecycle.back();
	in->position_in_set() = i;
	mVertexRecycle.pop_back();
	mVertices[i] = in;
      }
      else{
	mVertices.push_back(in);
	in->position_in_set() = mVertices.size()-1;
      }
    }
		
    void push_edge(edge_ptr  in){			
      if (mEdgeRecycle.size() > 0 && !manual_clean_up){
	int i = mEdgeRecycle.back();
	in->position_in_set() = i;
	mEdgeRecycle.pop_back();
	mEdges[i] = in;
      }
      else{
	in->position_in_set() = mEdges.size();
	mEdges.push_back(in);
      }
    }
		
    void push_face(face_ptr in){
      if (mFaceRecycle.size() > 0 && manual_clean_up){
	int i = mFaceRecycle.back();
	in->position_in_set() = i;
	mFaceRecycle.pop_back();
	mFaces[i] = in;
      }
      else{
	mFaces.push_back(in);
	in->position_in_set() = mFaces.size()-1;
      }
    }
		
    void remove_vertex(int i){			
      vertex_ptr v = mVertices[i];
      //std::cout << " removing vert: " << i << " " << v  << " " << v->size() << std::endl;
      mVertexRecycle.push_back(i);
      if(!manual_clean_up) {
	//if(v->size() > 1) throw("fuck you fucktard");
	mVertices[i] = NULL; 	
	delete v;
      }
    }
		
    void remove_edge(int i){
      edge_ptr e = mEdges[i];
      //if(i == 175548) throw('gdb here we come!');
      if(this->mEdgeDeleteFunc)
	this->mEdgeDeleteFunc(e);
      mEdgeRecycle.push_back(i);
      
      if(!manual_clean_up) {
	delete mEdges[i];
	mEdges[i] = NULL; 
      }
    }
		
    void remove_face(int i){
      face_ptr f = mFaces[i];
      mFaceRecycle.push_back(i);
      if(!manual_clean_up) {
	mFaces[i] = NULL; 
	delete f; 
      };
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
		
    void toggle_clean_up(){manual_clean_up ^= true;}

    vertex_ptr insert_vertex(coordinate_type in){
      return this->insert_vertex(in[0], in[1], in[2]);
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
	int j = 0;
	for (int i = 0; i < mFaces.size(); i++) {
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
	int j = 0;
	for (int i = 0; i < mVertices.size(); i++) {
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
	int j = 0;
	for (int i = 0; i < mEdges.size(); i++) {
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
      for(int i = 0; i < this->mVertices.size(); i++){
	//mVertices[i]->pack();
      }
    }
    
    void validateFaceVertices(){
      for(int i = 0; i < mVertices.size(); i++){
	vertex_ptr v = mVertices[i];
	if(v){
	  face_vertex_ptr fvb = v->fbegin();
	  face_vertex_ptr fve = v->fend();
	  bool iterating = true;
	  while(iterating){
	    if(fvb->vertex() != v){
	      std::cout << "face vertex: " << fvb 
			<< " has improperly assigned vertex!" << std::endl;
	    }
	    iterating = fvb != fve;
	    fvb = fvb->vnext();
	  }
	}
      }

      for(int i = 0; i < mFaces.size(); i++){
	face_ptr f = mFaces[i];
	if(f){
	  face_vertex_ptr fvb = f->fbegin();
	  face_vertex_ptr fve = f->fend();
	  bool iterating = true;
	  while(iterating){
	    if(fvb->face() != f){
	      std::cout << "face vertex: " << fvb 
			<< " has improperly assigned face!" << std::endl;
	    }
	    iterating = fvb != fve;
	    fvb = fvb->next();
	  }
	}
      }
      for(int i = 0; i < mEdges.size(); i++){
	edge_ptr e = mEdges[i];
	if(e){
	  if(e->v1()->edge() != e){
	    std::cout << "face vertex: " << e->v1()->edge()
		      << " has improperly assigned edge!" << std::endl;
	  }
	  if(e->v2()->edge() != e){
	    std::cout << "face vertex: " << e->v2()->edge()
		      << " has improperly assigned edge!" << std::endl;
	  }
	}
      }
    }
#if 1
    void groupElements(){

      bool coloring = true;
      maxGraphColor = 0;
      for(int i = 0; i < mVertices.size(); i++){
	if(mVertices[i])
	  mVertices[i]->group() = -1;
      }
      for(int i = 0; i < mEdges.size(); i++){
	if(mEdges[i]){
	  mEdges[i]->group() = -1;
	  mEdges[i]->v1()->group() = -1;
	  mEdges[i]->v2()->group() = -1;
	}
      }
      for(int i = 0; i < mFaces.size(); i++){
	if(mFaces[i])
	  mFaces[i]->group() = -1;
      }
      
      int currentGroup = 0;
      T r = (T)rand()/(T)RAND_MAX;
      T g = (T)rand()/(T)RAND_MAX;
      T b = (T)rand()/(T)RAND_MAX;
      mGroupColors.push_back(coordinate_type(r,g,b));

      for(int i = 0; i < mVertices.size(); i++){
	vertex_ptr vi = mVertices[i];
	if(vi->group() == -1){
	  currentGroup++;
	  T r = (T)rand()/(T)RAND_MAX;
	  T g = (T)rand()/(T)RAND_MAX;
	  T b = (T)rand()/(T)RAND_MAX;
	  if(currentGroup > mGroupColors.size()-1)
	    mGroupColors.push_back(coordinate_type(r,g,b));
	}
	else continue;
	std::stack<int> stack;
	stack.push(i);
	while(stack.size() > 0){
	  int pId = stack.top(); stack.pop();
	  vertex_ptr v = mVertices[pId];
	  v->group() = currentGroup;
	  if(v && v->size() > 0){
	    face_vertex_ptr fvb = v->fbegin();
	    face_vertex_ptr fve = v->fend();
	    bool iterating = true;
	    while(iterating){	    
	      iterating = fvb != fve;
	      fvb->group()  = currentGroup;
	      fvb->edge()->group() = currentGroup;
	      fvb->face()->group() = currentGroup;
	      if(fvb->next()->vertex()->group() == -1){ 
		stack.push(fvb->next()->vertex()->position_in_set());
	      }
	      //fvb->face()->color.r = mGroupColors[currentGroup][0];
	      //fvb->face()->color.g = mGroupColors[currentGroup][1];
	      //fvb->face()->color.b = mGroupColors[currentGroup][2];
	      fvb = fvb->vnext();
	    }	  
	  }
	}
      }
    }
#endif

    void colorVertices(){
      //greedy coloring
      bool coloring = true;
      maxGraphColor = 0;
      for(int i = 0; i < mVertices.size(); i++){
	if(mVertices[i])
	  mVertices[i]->graphColor = -1;
      }

      vertex_array permVerts = mVertices;
      for(int i = 0; i < permVerts.size(); i++){
	int card = rand()%permVerts.size();
	vertex_ptr vt = permVerts[i];
	permVerts[i] = permVerts[card];
	permVerts[card] = vt;
      }

      for(int i = 0; i < permVerts.size(); i++){
	vertex_ptr v = permVerts[i];
	if(v && v->size() > 0){
	  int maxNeighborColor = 0;
	  int minNeighborColor = 0;
	  face_vertex_ptr fvb = v->fbegin();
	  face_vertex_ptr fve = v->fend();
	  bool iterating = true;
	  while(iterating){
	    iterating = fvb != fve;
	    int neighborColor = fvb->next()->vertex()->graphColor;
	    maxNeighborColor = neighborColor > maxNeighborColor ?
	      neighborColor : maxNeighborColor;
	    if(neighborColor > 0)
	      minNeighborColor = neighborColor < minNeighborColor ?
		neighborColor : minNeighborColor;
	    fvb = fvb->vnext();
	  }
	  //std::cout << minNeighborColor << std::endl;
	  //std::cout << maxNeighborColor << std::endl;
	  if(minNeighborColor - 1 > -1){
	    v->graphColor = minNeighborColor - 1;
	  }
	  else{
	    v->graphColor = maxNeighborColor + 1;
	  }
	  maxGraphColor = 
	    maxGraphColor > maxNeighborColor + 1 ? 
	    maxGraphColor : maxNeighborColor + 1;
	  
	}
      }
    }

    void setFaceFlags(int k){
    for(int i = 0; i < mFaces.size(); i++){
	if(mFaces[i])
	  mFaces[i]->flag = k;
      }
    }
    /*	
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
		
    void draw_normals(){
      fa_iterator it_b = mFaces.begin();
      fa_iterator it_e = mFaces.end();
      while (it_b != it_e) {
	if(*it_b){
	  (*it_b)->draw_normal(0.0);
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
    */
        
    void reset_flags(){
      for(int i = 0; i < mFaces.size(); i++){
	if (mFaces[i]) {
	  mFaces[i]->flag = 0;
	}
      }
            
      for(int i = 0; i < mEdges.size(); i++){
	if(mEdges[i]){ 
	  mEdges[i]->flag = 0;
	  mEdges[i]->v1()->flag = 0;
	  mEdges[i]->v2()->flag = 0;
	}
      }
            
      for(int i = 0; i < mVertices.size(); i++){
	if(mVertices[i]){
	  mVertices[i]->flag = 0;
	}
      }			
    }

    void color_dead_pointers(){
      for(int i = 0; i < mFaces.size(); i++){
	if (mFaces[i]) {
	  if(!mFaces[i]->fend() || !mFaces[i]->fbegin()){
	    mFaces[i]->color.r = 1.0;
	    mFaces[i]->color.g = 0.0;
	    mFaces[i]->color.b = 0.0;
	    std::cout << "bad bad face" << std::endl;		  
	  }
	}
      }
            
      for(int i = 0; i < mEdges.size(); i++){
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
      for(int i = 0; i < mVertices.size(); i++){
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
      for(int i = 0; i < mFaces.size(); i++){
	if (mFaces[i]) {
	  mFaces[i]->update_all();
	}
      }
            
      for(int i = 0; i < mVertices.size(); i++){
	if (mVertices[i]) {
	  mVertices[i]->update_normal();
	}
      }     
    }
    
    void set_edge_delete_func(std::function<void(edge_ptr)> func){
      this->mEdgeDeleteFunc = func;
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
    bool                manual_clean_up;
    int                 numGroups;
    vector<int>         groupSizes;
    face_array		mFaces;     vector<int> mFaceRecycle;
    edge_array		mEdges;     vector<int> mEdgeRecycle;
    vertex_array	mVertices;  vector<int>  mVertexRecycle;
    coordinate_array mGroupColors;
    std::function<void(edge_ptr)>  mEdgeDeleteFunc;
  public:
    int maxGraphColor;
  };	

}
//#undef TYPEDEF_LIST
#endif
