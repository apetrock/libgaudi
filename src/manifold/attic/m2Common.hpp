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
#include <GLUT/glut.h>
#include <OpenGL/OpenGL.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>	
#include "al_Vec.hpp"
#include "al_Quat.hpp"

// this typedef list is ugly but useful!!!

#define M2_TYPEDEFS							\
  typedef typename SPACE::coordinate_type  coordinate_type;		\
  typedef typename SPACE::rotor_type       rotor_type;	         	\
  typedef typename SPACE::type			 T;			\
  typedef typename SPACE::Real		    Real;			\
  typedef	m2::face<SPACE>		    face_type;                  \
  typedef	m2::edge<SPACE>             edge_type;			\
  typedef	m2::vertex<SPACE>	    vertex_type;	        \
  typedef	m2::face_vertex<SPACE>      face_vertex_type;		\
  typedef	m2::surf<SPACE>	    surf_type;	        \
									\
  typedef	m2::face<SPACE>*	    face_ptr;	                \
  typedef	m2::edge<SPACE>*	    edge_ptr;	                \
  typedef	m2::vertex<SPACE>*          vertex_ptr;			\
  typedef	m2::face_vertex<SPACE>*     face_vertex_ptr;		\
  typedef	m2::surf<SPACE>*         surf_ptr;		\
									\
  typedef	m2::face<SPACE>&	    face_ref;	                \
  typedef	m2::edge<SPACE>&	    edge_ref;	                \
  typedef	m2::vertex<SPACE>&          vertex_ref;			\
  typedef	m2::face_vertex<SPACE>&     face_vertex_ref;		\
  typedef	m2::surf<SPACE>&         surf_ref;		\
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

  /*TODO: 
    1) Need to add a place for storing arbitrary data in SPACE so that there 
       is an optional data rider associated with each part of the graph
    2) There needs to be a face iterator class and a vertex iterator class
    3) An iterator that will iterate between two directed pieces of geometry
       ie move across directed edges or faces with four sides, so one can implement
       better collection strategies
  */

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
	
  struct colorRGB {
    double r,g,b,a;
  };
	
  inline void gl_set_color(colorRGB in){
    glColor4f(in.r, in.g, in.b, in.a);
  }
  inline void gl_set_color_opaque(colorRGB in){
    glColor4f(in.r, in.g, in.b, 1.0);
  }	

  template <typename T>
  inline ostream &operator<<( ostream &out, const al::Vec<3,T>& in ) {
    out << " " << in[0] << " " <<  in[1] << " " <<  in[2] << " ";
    return out;
  }
	
  template <typename T>
  inline	T randd(T range_) {	
    T output =  (T)(rand()/(T(RAND_MAX)+1.0))*range_;
    return output;
  }
  template <typename T>	
  inline al::Vec<4,T> cross(al::Vec<4,T> a,al::Vec<4,T> b){
    return al::Vec<4,T>( a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x, 1.0);
  }
	

}
//#undef TYPEDEF_LIST
#endif
