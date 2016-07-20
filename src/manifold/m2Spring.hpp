//
//  spring.h
//  Manifold
//
//  Created by John Delaney on 5/19/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#ifndef __TWOMANIFOLDSPRING__
#define __TWOMANIFOLDSPRING__

#include "EIGEN.h"
#include <stack> 
#include "m2Includes.h"
#include "conj_grad.hpp"
#include "modify.hpp"
#include "Octree.hpp"
#include <cmath>

namespace m2 {
    
  template <typename SPACE>
  class spring{ M2_TYPEDEFS;
  typedef coordinate_type							point_type ;
    typedef typename list<coordinate_type >::iterator	pl_iterator;
  public:
    spring(control_ptr obj_in){
      obj_in->pack();
      mks	= 15.0;//for now uniform k spring
      mdt	= 0.1;
      mMesh = obj_in;
      vertex_array vl = mMesh->get_vertices();
            
      this->initialize();
    }
        
    void initialize(){
      mM.resize(mMesh->get_vertices().size());
      pos0.resize(mMesh->get_vertices().size());
      vel0.resize(mMesh->get_vertices().size());
      vel1.resize(mMesh->get_vertices().size());
      mF.resize(mMesh->get_vertices().size());
      mFext.resize(mMesh->get_vertices().size());
      mFint.resize(mMesh->get_vertices().size());
			
      mIsPinned.resize(mMesh->get_vertices().size());
      ml0.resize(mMesh->get_edges().size());
      mP0.resize(mMesh->get_vertices().size());
			
      mKr.resize(6*mMesh->get_vertices().size());			
      mK.resize(6*mMesh->get_edges().size());
      vector<edge_ptr> ev = mMesh->get_edges();
			
      for (long i = 0; i < ev.size(); i++) {
	edge_ptr e = ev[i];
	coordinate_type c1 = e->v1()->coordinate();
	coordinate_type c2 = e->v2()->coordinate();
	T l = dist(c1,c2);
	ml0[i] = l;
      }
			
      for (long i = 0; i < mM.size(); i++) {
	mM[i] = 1.0;
      }
			
      for (long i = 0; i < mIsPinned.size(); i++) {
	mIsPinned[i] = false;
      }
			
      vector<vertex_ptr> V = mMesh->get_vertices();
      for (long i = 0; i < pos0.size(); i++) {				
	mP0[i] = V[i]->coordinate();
	pos0[i] = V[i]->coordinate();
      }
    }
		
		
    al::Vec<3,T>& get_force(long i){
      return mFint[i];
    }
		
    void calculate_forces_and_stiffness(){
      edge_array E = mMesh->get_edges();
      for (long i = 0; i < mFint.size(); i++) {
	mFint[i] = 0;
      }
			
      for (long i = 0; i < E.size(); i++) {
	edge_ptr e = E[i];
	coordinate_type c1 = e->v1()->coordinate();
	coordinate_type c2 = e->v2()->coordinate();
	long v1 = e->v1()->vertex()->position_in_set();
	long v2 = e->v2()->vertex()->position_in_set();
				
	T l = dist(c1,c2);
	T l0 = ml0[i];
				
	//T kd = mks/l0;
	coordinate_type fs,fd;
	fs = mks*(c2-c1)*(l-l0)/l;				
	//fd = kd*(vel0[v1]-vel0[v2])*(c1-c2)/l;				
	fd = 0;
				
	mFint[v1] += mIsPinned[v1] ? 0 : (fs+fd);
	mFint[v2] -= mIsPinned[v2] ? 0 : (fs+fd);
				
				
	T dx = c2[0] - c1[0];
	T dy = c2[1] - c1[1];
	T dz = c2[2] - c1[2];				
				
	T linv  = 1./l;
	T
	  k00 = mks*(-1. + l0*linv*(1. -(dx*dx*linv*linv))),
	  k11 = mks*(-1. + l0*linv*(1. -(dy*dy*linv*linv))),
	  k22 = mks*(-1. + l0*linv*(1. -(dz*dz*linv*linv))),				
	  k01 = -mks*l0*linv*(dx*dy*linv*linv),
	  k02 = -mks*l0*linv*(dx*dz*linv*linv),
	  k12 = -mks*l0*linv*(dy*dz*linv*linv);
				
	mK[6*i + 0] = k00;
	mK[6*i + 1] = k11;
	mK[6*i + 2] = k22;
	mK[6*i + 3] = k01;
	mK[6*i + 4] = k02;
	mK[6*i + 5] = k12;
      }
    }
		
    void set_pin(int i, bool isPinned){
      mIsPinned[i] = isPinned;
      pinID.push_back(i);
    }
		
    void update_pinned_positions(){
      vertex_array& V = mMesh->get_vertices();
      for (long i = 0; i < pinID.size(); i++) {
	V[pinID[i]]->coordinate() = mP0[V[pinID[i]]->position_in_set()];
      }
    }
		
    virtual void update_positions(){
      this->solve();
      this->update_pinned_positions();
      //fill mesh data structure, update lengths
    }
		
    virtual void solve(){			
      //this->explicite_euler();
      this->implicite_euler();			
      //this->verlet();
      mMesh->update_all();
    }
		
    virtual void explicite_euler(){
      this->calculate_forces_and_stiffness();
      vector<vertex_ptr> V = mMesh->get_vertices();		
      for (long i = 0; i < vel1.size(); i++) {
	vel1[i] = mdt/mM[i]*(mFint[i] + mFext[i]) + vel0[i];
	V[i]->coordinate() += (vel1[i])*mdt;				
      }
      swap(vel0, vel1);
    }
		
    virtual void verlet(){
			
      vector<vertex_ptr> V = mMesh->get_vertices();		
      for (long i = 0; i < pos0.size(); i++) {				
	pos0[i] = V[i]->coordinate();
      }
      this->calculate_forces_and_stiffness();
      for (long i = 0; i < vel1.size(); i++) {
	V[i]->coordinate() = pos0[i] + vel0[i]*mdt + (mFint[i]+mFext[i])*mdt*mdt/mM[i];
	vel1[i] = (V[i]->coordinate() - pos0[i])/mdt;
      }
      swap(vel0, vel1);
			
    }
		
    virtual void implicite_euler(){			
			
      vector<vertex_ptr> V = mMesh->get_vertices();
			
      for (long i = 0; i < V.size(); i++) {				
	pos0[i] = V[i]->coordinate();
      }
			
      int N = 0;
      while (N < 3) {
				
	this->calculate_forces_and_stiffness();				
	for(long i = 0; i < vel0.size(); i++){
	  mF[i] = mM[i]*vel0[i] + mdt*mFint[i] + mdt*mFext[i];
	}				
	cgsolve(*this,vel1,mF);
				
	for (long i = 0; i < V.size(); i++) {				
	  V[i]->coordinate() = pos0[i] + mdt*vel1[i];//*0.5 + mdt*vel0[i]*0.5;
	}				
				
				
	swap(vel0, vel1);
	N++;
      }			
			
    }
		
    virtual void mult(vector<coordinate_type>& x, vector<coordinate_type>& b){
			
      vector<edge_ptr> E = mMesh->get_edges();		
      vector<vertex_ptr> V = mMesh->get_vertices();			
			
      for (long i = 0; i < x.size(); i++) {
	vertex_ptr v = V[i];
	coordinate_type& cii = x[i];
				
	if (!mIsPinned[i]){b[i] = cii*mM[i];}
	else {b[i] = cii;}							
				
	face_vertex_ptr fvb = mMesh->vertex(i)->fbegin();
	face_vertex_ptr fve = mMesh->vertex(i)->fend();
	coordinate_type fadj = 0;
				
	while (fvb != fve) {
	  edge_ptr e = fvb->edge();
	  int ei = e->position_in_set();					
	  long vi1 = e->v1()->vertex()->position_in_set();
	  long vi2 = e->v2()->vertex()->position_in_set();
					
	  long eij = (e->v1() == fvb) ? vi2 : vi1;
					
	  coordinate_type& cij = x[eij];
	  coordinate_type  bij = 0;
	  T k00,k11,k22,k01,k02,k12;
					
					
	  if (!mIsPinned[i]) {
	    k00 = mK[6*ei + 0];
	    k11 = mK[6*ei + 1];
	    k22 = mK[6*ei + 2];
	    k01 = mK[6*ei + 3];
	    k02 = mK[6*ei + 4];
	    k12 = mK[6*ei + 5];
	  }
					
	  else {
	    k00 = 0; k11 = 0; k22 = 0;
	    k01 = 0; k02 = 0; k12 = 0;
	  }
					
	  //					b[i][0] -= k00*(cii[0]) + k01*(cii[1]) + k02*(cii[2]);
	  //					b[i][1] -= k01*(cii[0]) + k11*(cii[1]) + k12*(cii[2]);
	  //					b[i][2] -= k02*(cii[0]) + k12*(cii[1]) + k22*(cii[2]);
					
	  bij[0]  += k00*(cij[0]) + k01*(cij[1]) + k02*(cij[2]);
	  bij[1]  += k01*(cij[0]) + k11*(cij[1]) + k12*(cij[2]);
	  bij[2]  += k02*(cij[0]) + k12*(cij[1]) + k22*(cij[2]);;
					
	  b[i]  -= mdt*mdt*bij; //because kij,kji is negative and Kij = M[i,j] - dt*dt*K[i,j] and M[i,j] = 0;										
	  fvb = fvb->vnext();
	}
      }
    }
		
		
		
    void randomize(){
      vertex_array V = mMesh->get_vertices();
      for (long i = 0; i < V.size(); i++) {
	V[i]->coordinate()[0] += randd(0.1);
	V[i]->coordinate()[1] += randd(0.1);
	V[i]->coordinate()[2] += randd(0.1);
	pos0[i] = V[i]->coordinate();
      }
    }
		
    void draw(){
    }
		
    void reset_t(){
      mdt = 0.1;
    }
		
  private:
    control_ptr		 mMesh;
    T mdt,meps,mks,mkrs;
    vector<coordinate_type>	pos0;	//new velocity		
    vector<coordinate_type>	vel1;	//new velocity
    vector<coordinate_type>	vel0;	//old velocity			
    vector<coordinate_type>	mFext;		//external forces		
    vector<coordinate_type>	mFint;		//external forces				
    vector<coordinate_type>	mF;			//external forces
		
    //pins
    vector<bool>			mIsPinned;
    vector<long>			pinID;
    vector<coordinate_type> mP0;		
		
    vector<T>	mM;					//mass matrix		
    vector<T>	mK;					//tangential stiffness matrix one for each edge
    vector<T>	mKr;				//tangential stiffness matrix from rest position
    vector<T>   ml0;				//initial lengths
  };
	
  //--------------------------------------------------------------------------------------
	
  template <typename SPACE>
  class spring_electric{ 
    M2_TYPEDEFS;
    typedef coordinate_type							point_type ;
    typedef typename list<coordinate_type >::iterator	pl_iterator;

  public:
    spring_electric(control_ptr obj_in){
      K_	= 1.0;
      C_	= 1.0;
      eps_= 0.1;
      dt_	= 0.1;
      mMesh = obj_in;
      vertex_array vl = mMesh->get_vertices();
			
      va_iterator itb = vl.begin();
      va_iterator ite = vl.end();			
      while (itb != ite) {
	pos1.push_back((*itb)->coordinate());
	++itb;
      }			
    }
    virtual void update_positions(){
      vertex_array vl = mMesh->get_vertices();
      size_t gs = vl.size();
			
      pl_iterator p0tb = pos1.begin();
      pl_iterator p0te = pos1.end();
			
      va_iterator itb = vl.begin();
      va_iterator ite = vl.end();
			
      while (itb != ite) {
	double ax, ay, az;
	ax = 0.0; ay = 0.0; az = 0.0;
	double a_norm = 0.0; 
	coordinate_type &posi = (*itb)->coordinate();
				
	bool nool = true;
				
	va_iterator jtb = vl.begin();
	va_iterator jte = vl.end();
				
	coordinate_type ni = (*itb)->normal();
	ni.normalize();
	while (jtb != jte) {
	  vertex_ptr vi = *itb;
	  vertex_ptr vj = *jtb;
	  coordinate_type &posj = (*jtb)->coordinate();
	  if (vi->ID() != vj->ID()) {
	    coordinate_type nj = (*jtb)->normal();
	    double dx, dy, dz;
	    dx = 0.0; dy = 0.0; dz = 0.0;
						
	    nj.normalize();
						
	    coordinate_type dij = cross(posi, posj);
						
	    coordinate_type nij;
	    nij[0] = (nj[0] - ni[0]);
	    nij[1] = (nj[1] - ni[1]);
	    nij[2] = (nj[2] - ni[2]);
						
	    dij.normalize();
	    nij.normalize();
	    T d_ = 0.0, n_ = -0.0;
						
	    dx = (posj[0] - posi[0]);
	    dy = (posj[1] - posi[1]);
	    dz = (posj[2] - posi[2]);						
						
	    T inv1 = 1/(dx*dx + dy*dy + dz*dz + eps_);
	    T f = -C_*inv1*inv1;//  - 10.*inv1*inv1*inv1*inv1;
	    ax += f*dx + d_*inv1*dij[0] + n_*inv1*nij[0];
	    ay += f*dy + d_*inv1*dij[1] + n_*inv1*nij[1];
	    az += f*dz + d_*inv1*dij[2] + n_*inv1*nij[2];						
	  }
	  ++jtb;					
	}
				
				
	face_vertex_ptr ktb = (*itb)->fbegin();
	face_vertex_ptr kte = (*itb)->fend();
	bool iterating = true;
	while (iterating) {
	  iterating = ktb != kte;
	  coordinate_type &posj = ktb->edge()->other(*itb)->vertex()->coordinate();
	  double dx, dy, dz;
	  dx = 0.0; dy = 0.0; dz = 0.0;
	  dx = posj[0] - posi[0];
	  dy = posj[1] - posi[1];
	  dz = posj[2] - posi[2];
	  double dist = sqrt(dx*dx + dy*dy + dz*dz + eps_);
	  double f = dist*dist/K_;
	  ax += f*dx;
	  ay += f*dy;
	  az += f*dz;
	  ktb = ktb->vnext();
	}
				
				
	a_norm = sqrt(ax*ax + ay*ay + az*az);
				
	coordinate_type &pos1i = (*p0tb);
	pos1i[0] = posi[0] + dt_*ax/a_norm;
	pos1i[1] = posi[1] + dt_*ay/a_norm;
	pos1i[2] = posi[2] + dt_*az/a_norm;
				
				
	++itb; ++p0tb;
      }
			
      itb = vl.begin();
      p0tb = pos1.begin();
			
      while (itb != ite) {
	coordinate_type &posi = (*itb)->coordinate();
	coordinate_type &pos1i = *p0tb;
				
	swap(posi, pos1i);
				
	itb++;
	p0tb++;
      }
			
      mMesh->update_all();
      //dt_ = 0.98*dt_;
    }
    void randomize(){
      vertex_array vl = mMesh->get_vertices();
      va_iterator itb = vl.begin();
      va_iterator ite = vl.end();
      while (itb != ite) {
	(*itb)->coordinate()[0] = randd(0.1);
	(*itb)->coordinate()[1] = randd(0.1);
	(*itb)->coordinate()[2] = randd(0.1);
				
	itb++;
      }
    }
    void reset_t(){
      dt_ = 0.1;
    }
		
  private:
    control_ptr		 mMesh;
    T dt_,eps_,C_,K_;
    list<coordinate_type>	pos1;        
  };
}//end namespace
#endif
