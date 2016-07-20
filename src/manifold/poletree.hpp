/*
 *  GLOtree.cpp
 *  OTree
 *
 *  Created by John Delaney on 4/5/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */



#ifndef __POLETREE__
#define __POLETREE__

#include <iostream>
#include <math.h>
#include "al_Vec.hpp"


// COTD Entry submitted by Paul Nettle [midnight@FluidStudios.com]
// Corresponds with an Ask MidNight response (http://www.flipcode.com/askmid/)

// -----------------------------------------------------------------------------
// This defines a traverse_callback for traversal
// -----------------------------------------------------------------------------

//forward declaration
template <typename T, typename POINT_TYPE>
class   Octree;


template <typename T, typename POINT_TYPE>
class octree_traversal_base {
public:
  octree_traversal_base(){}
  ~octree_traversal_base(){}	
  virtual bool operator()(const Octree<T,POINT_TYPE> &o, void *data){ return true; };
};

template <typename T, typename POINT_TYPE>
class octree_build_base {
public:	
  octree_build_base(){}
  ~octree_build_base(){}	
  virtual bool operator()(const Octree<T,POINT_TYPE> &o, void *data){return true;}
};

template <typename T, typename POINT_TYPE>
class draw_leaf_node : public octree_traversal_base<T,POINT_TYPE> {
public:	
  draw_leaf_node(){};
  ~draw_leaf_node(){};
  bool operator()(const Octree<T,POINT_TYPE> &o, void *data){
    bool out = true;
    T r = o.radius();
    al::Vec<3,T> cen = o.center();
		
    if (o.is_leaf()) {
      glLineWidth(0.25);
      glColor4f(0.50 , 0.50, 0.70, 0.6);
      draw_box(cen[0]-r, cen[0]+r, 
	       cen[1]-r, cen[1]+r,
	       cen[2]-r, cen[2]+r);
      vector<POINT_TYPE*> cl = o.get_points();
			
      //			pl_iterator itb = cl.begin();
      //			pl_iterator ite = cl.end();
			
      //		while (itb != ite) {
      //			POINT_TYPE& pt = **itb;
      //			pt.draw();
      //			++itb;
      //		}
    }
    return out;
  }
};

template <typename T, typename POINT_TYPE>
bool get_child_points(const Octree<T,POINT_TYPE> &o, void *data){
  vector<POINT_TYPE*>& collect = *(vector<POINT_TYPE*>*)data;
	
  if (o.is_leaf()) {
    vector<POINT_TYPE*> opts = o.points();
    collect.merge(opts);
    return false;
  }
  else {
    return true;
  }
}

template <typename T, typename POINT_TYPE>
bool draw_node(const Octree<T,POINT_TYPE> &o, void *data){
  typedef typename vector<POINT_TYPE*>::iterator pl_iterator;
  bool out = true;
  T r = o.radius();
  al::Vec<3,T> cen = o.center();
	
  glLineWidth(0.25);
  glColor4f(0.50 , 0.50, 0.70, 0.6);
  draw_box(cen[0]-r, cen[0]+r, 
	   cen[1]-r, cen[1]+r,
	   cen[2]-r, cen[2]+r);
	
  //	glColor4f(0.25, 0.25, 0.25, 0.0);
  //	
  //	if (o.pointCount() > 0) {
  //		list<POINT_TYPE*> cl = o.get_points();
  //		pl_iterator itb = cl.begin();
  //		pl_iterator ite = cl.end();
  //		glPointSize(5.0);
  //		glBegin(GL_POINTS);		
  //		
  //		while (itb != ite) {
  //			POINT_TYPE& pt = **itb;
  //			glVertex3d(pt[0], pt[1], pt[2]);
  //			itb++;
  //		}
  //		glEnd();
  //		
  //		itb = cl.begin();
  //		ite = cl.end();
  //		
  //		while (itb != ite) {
  //			glBegin(GL_LINES);
  //			POINT_TYPE& pt = **itb;
  //			al::Vec<3,T> c = o.center();
  //			glVertex3d(c[0], c[1], c[2]);
  //			glVertex3d(pt[0], pt[1], pt[2]);
  //			glEnd();
  //			itb++;
  //		}
  //	}
  //	if (o.pointCount() == 1) {
  //		out = false;
  //	}
  return out;
}

template <typename OCT, typename PROC, typename DATA>
const bool traverse(OCT& O, PROC& proc, DATA *data)
{
  // Call the traverse_callback for this node (if the traverse_callback returns false, then
  // stop traversing.
  bool out = true;
  if (!proc(O, data)) out = false;
		
  // If I'm a node, recursively traverse my children
  size_t cur_count = O.get_points().size();
  if (!cur_count)
    {
      for (unsigned int i = 0; i < 8; i++)
	{
	  // We store incomplete trees (i.e. we're not guaranteed
	  // that a node has all 8 children)
	  bool here = true;
	  if (!O.has_child(i)) continue;
				
	  if (!traverse(O.child(i),proc, data)) out = false;
	}
    }
		
  return out;
}

template <typename T, typename POINT_TYPE>
class   Octree
{
public:
  // Construction/Destruction
	
  // Accessors
  typedef int            (*traverse_callback)(const Octree &o, void *data);
  typedef bool			(*build_callback)(const Octree &o, 
						  POINT_TYPE*	p, 
						  al::Vec<3,T>	c, T r,
						  void			*data);
  inline  vector<POINT_TYPE*>		points()	 const {return _points;}
  inline  const   unsigned int    pointCount() const {return _points.size();}
	
  typedef Octree<T,POINT_TYPE> octree_node;
  typedef typename list<POINT_TYPE*>::iterator pl_iterator;
  typedef typename vector<POINT_TYPE*>::iterator pv_iterator;
  // Implementation
	
	
  // -----------------------------------------------------------------------------
  // Construction -- Just "nullify" the class
  // -----------------------------------------------------------------------------
	
  Octree()
    :  _center(0,0,0), _radius(0.0)
  {
    _points.clear();
    memset(_child, 0, sizeof(_child));
  }
	
  // -----------------------------------------------------------------------------
  // Destruction -- free up memory
  // -----------------------------------------------------------------------------
	
  ~Octree()
  {
  }
	
  // -----------------------------------------------------------------------------
  // Build the octree
  // -----------------------------------------------------------------------------
	
  const   bool    build(vector<POINT_TYPE*>  points,
			const unsigned int threshold,
			const unsigned int maximumDepth,
			const al::Vec<3,T> icenter,
			const T			 iradius,	
			const unsigned int currentDepth)
  {
		
    // You know you're a leaf when...
    //
    // 1. The number of points is <= the threshold
    // 2. We've recursed too deep into the tree
    //    (currentDepth >= maximumDepth)
    //
    //    NOTE: We specifically use ">=" for the depth comparison so that we
    //          can set the maximumDepth depth to 0 if we want a tree with
    //          no depth.
    this->_center = icenter;
    this->_radius = iradius;
		
		
    unsigned int count = points.size();
    //		if (count <= threshold || currentDepth >= maximumDepth)		
    size_t check = points.size();
    if (points.size() <= threshold || currentDepth >= maximumDepth)
      {
	// Just store the points in the node, making it a leaf						
	_points = points;		
			
	return true;
      }
		
    // We'll need this (see comment further down in this source)
    // set everything to zero, this isn't guaranteed.
    unsigned int    childPointCounts[8];
    for (int i = 0; i < 8; i++) {
      childPointCounts[i]=0;
    }
		
    // Classify each point to a child node
    pv_iterator itb = points.begin();
    pv_iterator ite = points.end();
    vector<int> flags;
		
    while (itb != ite) {
      // Current point		
      POINT_TYPE   *p = *itb;
			
      // Center of this node
			
      const al::Vec<3,T> &c = icenter;
			
      // Here, we need to know which child each point belongs to. To
      // do this, we build an index into the _child[] array using the
      // relative position of the point to the center of the current
      // node
      //			bool inbounds = this->check_bounds(p);
      //			if (!inbounds){
      //				bool stop_here= true;
      //			}
			
      int flag = 0;
      if ((*p)[0] > c[0]) flag |= 1;
      if ((*p)[1] > c[1]) flag |= 2;
      if ((*p)[2] > c[2]) flag |= 4;
			
      // We'll need to keep track of how many points get stuck in each
      // child so we'll just keep track of it here, since we have the
      // information handy.
      flags.push_back(flag);
      childPointCounts[flag]++;
      itb++;
    }
		
    // Recursively call build() for each of the 8 children		
    // Generate a new bounding volume -- We do this with a touch of
    // trickery...
    //
    // We use a table of offsets. These offsets determine where a
    // node is, relative to it's parent. So, for example, if want to
    // generate the bottom-left-rear (-x, -y, -z) child for a node,
    // we use (-1, -1, -1).
    // 
    // However, since the radius of a child is always half of its
    // parent's, we use a table of 0.5, rather than 1.0.
    // 
    // These values are stored the following table. Note that this
    // won't compile because it assumes Points are structs, but you
    // get the idea.
		
    T bounds_offset[] = 		
      {
	-0.5, -0.5, -0.5,			
	+0.5, -0.5, -0.5,
	-0.5, +0.5, -0.5,
	+0.5, +0.5, -0.5,
	-0.5, -0.5, +0.5,
	+0.5, -0.5, +0.5,
	-0.5, +0.5, +0.5,
	+0.5, +0.5, +0.5
      };
		
    for (unsigned int i = 0; i < 8; i++)
      {
	// Don't bother going any further if there aren't any points for
	// this child
			
	if (!childPointCounts[i]) continue;
			
	// Allocate the child
			
	_child[i] = new Octree;
			
	// Allocate a   of points that were coded JUST for this child
	// only
			
			
	// Go through the input list of points and copy over the points
	// that were coded for this child
			
	vector<POINT_TYPE*>   newList;
	newList.clear();
	itb = points.begin();
	ite = points.end();
	vector<int>::iterator ftb = flags.begin();
			
	while (itb!=ite) {
	  POINT_TYPE* pt = *itb;
	  int flag = *ftb;
	  if ((flag) == i)						
	    {
	      //itb = points.erase(itb);
	      newList.push_back(pt);
	      (*ftb) = 10;
	    }
	  itb++; ftb++;
	}
			
	// Calculate our offset from the center of the parent's node to
	// the center of the child's node
	al::Vec<3,T> offset;
	offset[0] = bounds_offset[3*i+0]*iradius;
	offset[1] = bounds_offset[3*i+1]*iradius;
	offset[2] = bounds_offset[3*i+2]*iradius;
			
	// Create a new Bounds, with the center offset and half the
	// radius
	T				nradius = iradius * 0.5;
	al::Vec<3,T>	ncenter = icenter + offset;
			
	// Recurse
			
	_child[i]->build(newList, threshold, maximumDepth,
			 ncenter, nradius, currentDepth+1);
	_child[i]->parent() = this;
	_child[i]->quadrant() = i;
			
      }
		
    return true;
  }
	
  bool is_leaf() const {
    return _points.size();
  }
	
  const   bool    kbuild(vector<POINT_TYPE*>  points,
			 build_callback	  proc,
			 const unsigned int threshold,
			 const unsigned int maximumDepth,
			 const al::Vec<3,T> icenter,
			 const T			 iradius,	
			 const unsigned int currentDepth)
  {
		
		
    this->_center = icenter;
    this->_radius = iradius;
		
    unsigned int count = points.size();
    //		if (count <= threshold || currentDepth >= maximumDepth)		
    size_t check = points.size();
    if (points.size() <= threshold || currentDepth >= maximumDepth)
      {
	// Just store the points in the node, making it a leaf						
	_points = points;		
			
	return true;
      }
		
    unsigned int    childPointCounts[8];
    for (int i = 0; i < 8; i++) {
      childPointCounts[i]=0;
    }
		
    pv_iterator itb = points.begin();
    pv_iterator ite = points.end();
		
    vector<int>			flags;
    vector<POINT_TYPE*> temp_list;
    //		inline unsigned int IX(unsigned int in,unsigned int  jn,unsigned int  kn){
    //			return ((in) + (im+2)*(jn) + (im+2)*(jm+2)*kn);
    //		}
    T bounds_offset[] = 		
      {						//zyx
	-0.5, -0.5, -0.5,	//000
	+0.5, -0.5, -0.5,	//001
	-0.5, +0.5, -0.5,	//010
	+0.5, +0.5, -0.5,	//011
	-0.5, -0.5, +0.5,	//100
	+0.5, -0.5, +0.5,	//101
	-0.5, +0.5, +0.5,	//110
	+0.5, +0.5, +0.5	//111
      };
		
		
    while (itb != ite) {
      POINT_TYPE   *p = *itb;
      for (int j = 0; j<8; j++) {
				
	al::Vec<3,T> offset;
	offset[0] = bounds_offset[3*j+0]*iradius;
	offset[1] = bounds_offset[3*j+1]*iradius;
	offset[2] = bounds_offset[3*j+2]*iradius;
	T				nradius = iradius * 0.5;
	al::Vec<3,T>	ncenter = icenter + offset;
				
	bool intx = false;
	intx = proc(*this,p,ncenter,nradius, NULL);
	if (intx) {
	  flags.push_back(j);
	  temp_list.push_back(p);
	  childPointCounts[j]++;					
	}
				
      }
			
      itb++;
    }
		
    for (int i = 0; i < 8; i++)
      {
	if (!childPointCounts[i]) continue;
			
	_child[i] = new Octree;
			
	vector<POINT_TYPE*>   newList;
	newList.clear();
	pv_iterator fitb = temp_list.begin();
	pv_iterator fite = temp_list.end();
	vector<int>::iterator ftb = flags.begin();
			
	while (fitb!=fite) {
	  POINT_TYPE* pt = *fitb;
	  int flag = *ftb;
	  //cout << flag << endl;
	  if ((flag) == i)					
	    {
	      //itb = points.erase(itb);
	      newList.push_back(pt);
	    }
	  fitb++; ftb++;				
	}
			
	// Calculate our offset from the center of the parent's node to
	// the center of the child's node
	al::Vec<3,T> offset;
	offset[0] = bounds_offset[3*i+0]*iradius;
	offset[1] = bounds_offset[3*i+1]*iradius;
	offset[2] = bounds_offset[3*i+2]*iradius;
			
	// Create a new Bounds, with the center offset and half the
	// radius
	T				nradius = iradius * 0.5;
	al::Vec<3,T>	ncenter = icenter + offset;
			
	// Recurse
			
	_child[i]->kbuild(newList, proc, threshold, maximumDepth,
			  ncenter, nradius, currentDepth+1);
	_child[i]->parent() = this;
	_child[i]->quadrant() = i;
			
      }
		
    return true;
  }
	
  // -----------------------------------------------------------------------------
  // Determine the [cubic] bounding volume for a set of points
  // -----------------------------------------------------------------------------
	
  T* calculate_cubic_bounds(vector<POINT_TYPE*>& vl )
  {
		
    size_t n = vl.size();
		
    size_t i = 0;
    T* out = new T[4];
		
    T
      xh=0,yh=0,zh=0,
      xl=0,yl=0,zl=0;
		
		
    typename vector<POINT_TYPE*>::iterator itb = vl.begin();
    typename vector<POINT_TYPE*>::iterator ite = vl.end();			
		
    while (itb != ite) {
      T xt, yt, zt;
      xt = (**itb)[0];
      yt = (**itb)[1];
      zt = (**itb)[2];
      //collect the extreme points these will give us the box radius to use
      if (xt > xh) xh = xt;
      if (yt > yh) yh = yt;
      if (zt > zh) zh = zt;
      if (xt < xl) xl = xt;
      if (yt < yl) yl = yt;
      if (zt < zl) zl = zt;
			
      i++; itb++;
    }
		
    out[0] = xl + 0.5*(xh-xl);
    out[1] = yl + 0.5*(yh-yl);
    out[2] = zl + 0.5*(zh-zl);
		
    out[3] = fabs(xh - xl);
    if (fabs(yh - yl) > out[3]) {
      out[3] = fabs(yh - yl);
    }
    if (fabs(zh - zl) > out[3]) {
      out[3] = fabs(zh - zl);
    }
    out[3] /= 2.;
		
    return out;			
		
  }
  //	
  // -----------------------------------------------------------------------------
  // Generic tree traversal
  // -----------------------------------------------------------------------------
	
  const bool traverse(octree_traversal_base<T,POINT_TYPE>& proc, void *data) const
  {
    // Call the traverse_callback for this node (if the traverse_callback returns false, then
    // stop traversing.
    bool out = true;
    if (!proc(*this, data)) out = false;
		
    // If I'm a node, recursively traverse my children
    size_t cur_count = _points.size();
    if (!cur_count)
      {
	for (unsigned int i = 0; i < 8; i++)
	  {
	    // We store incomplete trees (i.e. we're not guaranteed
	    // that a node has all 8 children)
	    bool here = true;
	    if (!_child[i]) continue;
				
	    if (!_child[i]->traverse(proc, data)) out = false;
	  }
      }
		
    return out;
  }
	
	
  al::Vec<3,T>&	center()		{		return _center;	} 
  al::Vec<3,T>	center() const	{		return _center;	} 
	
  T&				radius()		{	return _radius;	} 
  T				radius() const	{	return _radius;	}	
  octree_node*&			parent()			{	return _parent;	} 
  octree_node*			parent()	const	{	return _parent;	} 	
  int						quadrant()	const	{	return _quadrant;	}
  int&					quadrant()			{	return _quadrant;	} 

  vector<POINT_TYPE*>& get_points(){
    return _points;
  }
	
  vector<POINT_TYPE*> get_points() const {
    return _points;
  }
	
  Octree& child(int i){
    return *_child[i];
  }
	
  bool has_child(int i){
    if (_child[i]) return true;
    else		   return false;
  }	
	
  void draw() const{
    this->traverse(draw_node,NULL);
  }
	
  void draw_leaves() const{
    draw_leaf_node<T,POINT_TYPE> trav;
    this->traverse(trav,NULL);
  }
	
protected:
  Octree *_parent; 
  int	  _quadrant;	
  Octree  *_child[8];
  vector<POINT_TYPE*>	_points;
  al::Vec<3,T>		_center;
  al::Vec<3,T>		_centerOfMass;
  al::Vec<3,T>		_avgPotential;
  T			_radius;
};





#endif
