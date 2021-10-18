
/*
 *  m2FaceVertex.h
 *  Phase Vocoder
 *
 *  Created by John Delaney on 1/8/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef __TWOMANIFOLDFACEVERTEX__
#define __TWOMANIFOLDFACEVERTEX__
#include "m2Common.hpp"
#include "m2Control.hpp"
#include "m2Face.hpp"
#include "m2Vertex.hpp"
#include "m2Edge.hpp"

namespace m2 {
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
    }    
		
    face_vertex(const face_vertex_ref rhs){
      this->mEdge		= rhs.mEdge;
      //			mEdge->set_this(&rhs,this);
      this->mFace		= rhs.mFace;
      this->mVertex	= rhs.mVertex;
      this->nxt_face	= rhs.nxt_face;
      this->prv_face	= rhs.prv_face;			
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

    face_ref get_coface(){
      return mEdge->return_coface();
    }
    face_ref get_face(){
      return *mFace;
    }
		
    face_vertex_ptr add_next(){
      face_vertex_ptr	out = new face_vertex(*this);
      face_vertex_ptr nxt = this->nxt_face;
			
      out->next() = nxt;
      nxt->prev() = out;
			
      out->prev() = this;
      this->next() = out;
			
      this->edge()->set_this(this,out);
			
      mFace->size() += 1;
			
      return out;
    }
		
		
    face_vertex_ptr add_prev(){
      face_vertex_ptr	out = new face_vertex(*this);
      face_vertex_ptr prv = this->prv_face;
			
      out->next() = this;
      this->prev() = out;
			
      out->prev() = prv;
      prv->next() = out;
			
      prv->edge()->set_this(this,out);
			
      mFace->size() += 1;
			
      return out;
    }
		
    face_vertex_ptr add_next(vertex_ptr& pt){
      face_vertex_ptr	out = new face_vertex(*this);
      out->mVertex = pt;
      face_vertex_ptr nxt = this->nxt_face;
			
      out->next() = nxt;
      nxt->prev() = out;
			
      out->prev() = this;
      this->next() = out;
			
      mFace->size() += 1;
			
      //this->edge()->set_this(this,out);
			
      return out;
    }
		
		
    face_vertex_ptr add_prev(vertex_ptr& pt){
      face_vertex_ptr	out = new face_vertex(*this);
      out->mVertex = pt;
      face_vertex_ptr prv = this->prv_face;
			
      out->next() = this;
      this->prev() = out;
			
      out->prev() = prv;
      prv->next() = out;
			
      mFace->size() += 1;
			
      //prv->edge()->set_this(this,out);
			
      return out;
    }
		
    void insert_next(face_vertex_ptr& aft){
      face_vertex_ptr nxt = this->nxt_face;
			
      aft->next() = nxt;
      nxt->prev() = aft;
			
      aft->prev() = this;
      this->next() = aft;
			
    }
		
    void insert_previous(face_vertex_ptr& bef){
      face_vertex_ptr prv = this->prv_face;
			
      prv->next() = bef;
      bef->prev() = prv;
			
      bef->next() = this;
      this->prev() = bef;
		
    }
		
    void delete_next(){
      face_vertex_ptr nxt = this->nxt_face;
      face_vertex_ptr nxtnxt = nxt->next();
      this->next() = nxtnxt;
      nxtnxt->prev() = this;			  
			
    }
		
    face_vertex_ptr delete_this(){
      face_vertex_ptr prv = this->prv_face;
      face_vertex_ptr nxt = this->nxt_face;
      nxt->prv_face = prv;
      prv->nxt_face = nxt;
      return prv;
			
    }
		
    void remove_next(){}
    face_vertex_ptr & next() {return nxt_face;}		
    face_vertex_ptr & prev() {return prv_face;}
		
    face_vertex_ptr  vnext() {
      face_vertex_ptr out = this;
      if (mEdge == NULL) {
	return NULL;
      }
      else {
	out = mEdge->other(this);                    			
	return out->next();

      }			
    }
		
    face_vertex_ptr vprev()	{
      face_vertex_ptr out = this;
      if (out->prev() == NULL || out->prev()->mEdge == NULL) {
	return NULL;
      }
      else {
	out = out->prev();
	if (out->mEdge)	return out->mEdge->other(out);
	else return NULL;

      }
    }
		
    void draw(T off){
      this->draw_vertex(off);
      this->draw_tail(off);
    }

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
		
		
    edge_ptr & edge()		{return		mEdge;}
    face_ptr & face()		{return		mFace;}
    face_ptr & coface()	{return		mEdge->other(this)->face();}
    vertex_ptr & vertex()	{return		mVertex;}
		
    //		edge_type	edge()		const	{return	*	mEdge;}
    //		face_type	face()		const	{return	*	mFace;}
    //		vertex_type	vertex()	const	{return	*	mVertex;}
    long	vertex_ID()	const	{return mVertex->ID();}
    long	ID()		const	{return this->mID;}
    long	face_ID()	const	{return this->mFacePosition;}
    long&	face_ID()			{return this->mFacePosition;}
    fvl_iterator & position_in_vertex()  {return mVertexPosition;}

    long& position_in_face()    {return mFacePosition;}
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
    T& operator[](int i)		{return mVertex->coordinate()[i];}
    T  operator[](int i) const  {return mVertex->coordinate()[i];}

  protected:
		
    long			mID;
    long			fID;
		
    face_vertex_ptr nxt_face;
    face_vertex_ptr prv_face;

    fvl_iterator mVertexPosition;

    edge_ptr		mEdge;
    face_ptr		mFace;
    vertex_ptr		mVertex;
    //long mVertexPosition;
    long mFacePosition;
  public:
    unsigned int flag;
  };
	
}
#endif
