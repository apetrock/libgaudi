//
//  spring.h
//  Manifold
//
//  Created by John Delaney on 5/19/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#ifndef __M1SPRING__
#define __M1SPRING__

#include "conj_grad.hpp"
#include "m1.hpp"
#include "m2Common.hpp"

#include "vec_addendum.h"
#include <cmath>
//this will become a generic file, however, for now for simplicity, copy and paste.
namespace m1 {
	
  template <typename T>
  class spring;
  template <typename T>
  class vel_mult;
  template <typename T>
  class pos_mult;
	
  template <typename T>
  class vel_mult { M1_TYPEDEFS		
  public:
    m1::spring<T>* As;
    virtual void mult(vector<coordinate_type>& x, vector<coordinate_type>& b){
      vector<edge_ptr>   E = As->sp_->get_edges();		
      vector<vertex_ptr> V = As->sp_->get_vertices();			
      T dt = As->mdt*0.5;
      T m  = As->mM;
      int nthreads, tid, i, chunk;
      chunk = 10;
      //#pragma omp parallel private(i,tid)	
      {
	tid = omp_get_thread_num();
	if (tid == 0)
	  {
	    nthreads = omp_get_num_threads();
	  }
	//#pragma omp for schedule(dynamic,chunk)
				
	for (long i = 0; i < b.size(); i++) {
	  b[i] = x[i];
	  //b[i] = 0;
	}
				
	for (long i = 0; i < E.size(); i++) {
	  edge_ptr e = E[i];					
					
	  long v1 = e->v1()->position_in_set();
	  long v2 = e->v2()->position_in_set();
	  coordinate_type c1 = As->pos0[v1];
	  coordinate_type c2 = As->pos0[v2];

	  T l0 = As->ml0[i];
	  T l = norm(c1-c2);
					
	  coordinate_type dv = x[v1] - x[v2];
	  coordinate_type dc = c1-c2;
	  coordinate_type dn = dc/l;
	  T k = e->k;
	  T dotdv = dv*dn;
	  T A = dt*dt*dotdv*k/l0/m;
	  coordinate_type bij = dn*A;
	  bool pinned = As->mIsPinned[v1];
	  b[v1] +=  As->mIsPinned[v1] ? 0 : bij;
	  b[v2] -=  As->mIsPinned[v2] ? 0 : bij;
	}
      }
    }
		
  };

	
  template <typename T>
  class spring{ M1_TYPEDEFS
  typedef coordinate_type							point_type ;
    typedef typename list<coordinate_type >::iterator	pl_iterator;
    friend class vel_mult<T>;
    friend class pos_mult<T>;
  public:
    spring(control_ptr obj_in){

      mbegin = true;
      obj_in->pack();
      mks	= 0.0;//for now uniform k spring
      mdt	= 0.1;
      sp_ = obj_in;
      vertex_array vl = sp_->get_vertices();
            
      this->initialize();
    }
        
		
    void initialize(){
      mM = 1.0;
      pos0.resize(sp_->get_vertices().size());
      pos1.resize(sp_->get_vertices().size());			
      vel0.resize(sp_->get_vertices().size());
      vel1.resize(sp_->get_vertices().size());
      mF.resize(sp_->get_vertices().size());
      mF0.resize(sp_->get_vertices().size());			
      mFext.resize(sp_->get_vertices().size());
      mFint.resize(sp_->get_vertices().size());
			
      mIsPinned.resize(sp_->get_vertices().size());
      ml0.resize(sp_->get_edges().size());
      ml.resize(sp_->get_edges().size());			
      mP0.resize(sp_->get_vertices().size());
			
      mKr.resize(6*sp_->get_vertices().size());			
      mK.resize(6*sp_->get_edges().size());
      vector<edge_ptr> ev = sp_->get_edges();
			
      for (long i = 0; i < ev.size(); i++) {
	edge_ptr e = ev[i];
	coordinate_type c1 = e->v1()->coordinate();
	coordinate_type c2 = e->v2()->coordinate();
	T l = norm(c1-c2);
	ml0[i] = l;
      }
			
      //			for (long i = 0; i < mM.size(); i++) {
      //                mM = 1.0;
      //            }
			
      for (long i = 0; i < mIsPinned.size(); i++) {
	mIsPinned[i] = false;
      }
			
      vector<vertex_ptr> V = sp_->get_vertices();
      for (long i = 0; i < pos0.size(); i++) {
	coordinate_type nc0(V[i]->coordinate());
	coordinate_type nc1(V[i]->coordinate());
				
	mP0[i] = V[i]->coordinate();
	pos0[i] = nc0;
	pos1[i] = nc1;				
      }
    }
		
		
    al::Vec<3,T>& get_force(long i){
      return mFint[i];
    }
		
    al::Vec<3,T>& external_force(long i){
      return mFext[i];
    }
		
    size_t size(){
      return mFint.size();
    }
		
    m1::control<T>*&  graph()       {return sp_;}
    m1::control<T>*   graph() const {return sp_;}		
		
    T &M()       {return mM;}
    T  M() const {return mM;}
    T &dh()       {return mdt;}
    T  dh() const {return mdt;}
		
    coordinate_type calc_Fint(coordinate_type c1, coordinate_type c2, T l, T l0, T k){
      coordinate_type dc = c1 - c2;
      coordinate_type dn = dc/l;
      coordinate_type fs = 0;
			
      T A = (dot(dc,dn) - l0)*k/l0;					
      fs = dn*A;
      return fs;
    }
		
    void rhsVelocity(){
      edge_array E = sp_->get_edges();
			
      for (long i = 0; i < mFint.size() ; i++) {
	mFint[i] = 0;
      }
      int nthreads, tid, i, chunk;
      chunk = 10;
			
      //#pragma omp parallel private(i,tid)	
      {
	tid = omp_get_thread_num();
	if (tid == 0)
	  {
	    nthreads = omp_get_num_threads();
	  }
	//#pragma omp for schedule(dynamic,chunk)
				
				
	for (long i = 0; i < E.size(); i++) {
	  edge_ptr e = E[i];
	  long v1 = e->v1()->position_in_set();
	  long v2 = e->v2()->position_in_set();			
	  //					coordinate_type c1 = e->v1()->coordinate();
	  //					coordinate_type c2 = e->v2()->coordinate();
	  coordinate_type c1 = pos0[v1];
	  coordinate_type c2 = pos0[v2];
					
	  T k = e->k;
	  T l = norm(c1-c2);
	  T l0 = ml0[i];
	  coordinate_type fs = calc_Fint(c1,c2,l,l0,k);
	  mFint[v1] -= mIsPinned[v1] ? 0 : fs;
	  mFint[v2] += mIsPinned[v2] ? 0 : fs;					
	}
      }
    }
		
    void rhsPosition(){
      edge_array E = sp_->get_edges();
			
      for (long i = 0; i < mFint.size() ; i++) {
	mFint[i] = 0;
      }
			
      int nthreads, tid, i, chunk;
      chunk = 10;
			
      //#pragma omp parallel private(i,tid)	
      {
	tid = omp_get_thread_num();
	if (tid == 0)
	  {
	    nthreads = omp_get_num_threads();
	  }
	//#pragma omp for schedule(dynamic,chunk)

	T dt = mdt*0.5;
	for (long i = 0; i < E.size(); i++) {
	  edge_ptr e = E[i];
	  long v1 = e->v1()->position_in_set();
	  long v2 = e->v2()->position_in_set();					
	  coordinate_type c1 = pos0[v1];
	  coordinate_type c2 = pos0[v2];
	  //					coordinate_type c1 = e->v1()->coordinate();
	  //					coordinate_type c2 = e->v2()->coordinate();					
	  T l = dist(c1,c2);
	  T l0 = ml0[i];
	  coordinate_type dc = c1 - c2;				
	  coordinate_type dn = dc/l;
	  T dotDcn = dot(dc,dn);
	  coordinate_type xs,xp;
	  T k = e->k;

	  xs = calc_Fint(c1, c2, l, l0, k);
	  xp = dotDcn*k/l0;
	  coordinate_type xe = xs-xp;
	  mFint[v1] -= mIsPinned[v1] ? 0 : xe;
	  mFint[v2] += mIsPinned[v2] ? 0 : xe;
	}
      }
    }
		
    void set_pin(int i, bool isPinned){
      mIsPinned[i] = isPinned;
      pinID.push_back(i);			
    }
		
    bool is_pinned(int i){
      return mIsPinned[i];
    }
		
    void update_pinned_positions(){
      vertex_array& V = sp_->get_vertices();
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
      //sp_->update_all();
    }
		
    virtual void explicite_euler(){
      this->rhsVelocity();
      vector<vertex_ptr> V = sp_->get_vertices();		
      for (long i = 0; i < vel1.size(); i++) {
	vel1[i] = mdt/mM*(mFint[i] + mFext[i]) + vel0[i];
	V[i]->coordinate() += vel1[i]*mdt;
	pos0[i] = V[i]->coordinate();
      }
      swap(vel0, vel1);
    }
		
    virtual void verlet(){
			
      vector<vertex_ptr> V = sp_->get_vertices();		
      for (long i = 0; i < pos0.size(); i++) {				
	pos0[i] = V[i]->coordinate();
      }
      this->rhsVelocity();
      for (long i = 0; i < vel1.size(); i++) {
	pos1[i] = pos0[i] + vel0[i]*mdt + (mFint[i]+mFext[i])*mdt*mdt/mM;
	vel1[i] = (V[i]->coordinate() - pos0[i])/mdt;
      }
      swap(pos0, pos1);
      swap(vel0, vel1);
			
    }
		
    virtual void implicite_euler(){			
      vel_mult<T> Avel;			
      Avel.As = this;			
      T dt = mdt/2.;
			
      vector<vertex_ptr> V = sp_->get_vertices();
      //			if (mbegin) {
      //				this->explicite_euler();
      //				mbegin = false;
      //			}
			
      for (long i = 0; i < pos0.size(); i++) {
	pos0[i] = V[i]->coordinate();
      }
			
      //this->verlet();
      this->rhsVelocity();				
      for(long i = 0; i < vel0.size(); i++){
	vel0[i] = vel0[i] + (mFint[i] + mFext[i])*dt/mM;				
      }					
      cgsolve(Avel,vel1,vel0);
      for (long i = 0; i < pos0.size(); i++) {				
	pos0[i] += vel1[i]*dt;
      }
			
      vel0.swap(vel1);			
			
      this->rhsVelocity();				
      for(long i = 0; i < vel0.size(); i++){
	vel0[i] = vel0[i] + (mFint[i] + mFext[i])*dt/mM;		
      }	
      cgsolve(Avel,vel1,vel0);
				
      for (long i = 0; i < pos0.size(); i++) {
				
	V[i]->coordinate() += vel1[i]*mdt;
      }
			
      pos0.swap(pos1);
      vel0.swap(vel1);
			
      //			
      //			this->rhsPosition();			
      //			for (long i = 0; i < pos0.size(); i++) {					
      //				pos0[i] = pos0[i] + vel0[i]*dt + (mFint[i] + mFext[i])*dt*dt/mM;	//update for dt/2
      //			}
      //			
      //			cgsolve(Avel,pos1,pos0);
      //			
      //			for (long i = 0; i < V.size(); i++) {
      //				vel1[i] = (pos1[i] - pos0[i])/2./dt;
      //				//pos1[i] += vel1[i]*dt;
      //				V[i]->coordinate() = pos1[i];
      //			}
      //			
      //			swap(vel0, vel1);
      //			swap(pos0, pos1);				
    }		
		
		
		
    void randomize(){
      vertex_array V = sp_->get_vertices();
      for (long i = 0; i < V.size(); i++) {
	V[i]->coordinate()[0] += randd((T)0.1);
	V[i]->coordinate()[1] += randd((T)0.1);
	V[i]->coordinate()[2] += randd((T)0.1);
	pos0[i] = V[i]->coordinate();
      }
    }
		
    void draw(){
      glDisable(GL_DEPTH_TEST);					// Enables Depth Testing 
      glEnable(GL_BLEND);					
      glLineWidth(1.0);
      glBegin(GL_LINES);			
      for (long i = 0; i < mFint.size(); i++) {
	glColor3f(0.25, 0.25, 0.25);
	glVertex3f(sp_->vertex(i)->coordinate()[0],
		   sp_->vertex(i)->coordinate()[1],
		   sp_->vertex(i)->coordinate()[2]);
	glVertex3f(sp_->vertex(i)->coordinate()[0] + mFint[i][0],
		   sp_->vertex(i)->coordinate()[1] + mFint[i][1],
		   sp_->vertex(i)->coordinate()[2] + mFint[i][2]);
      }
      glEnd();
			
      glPointSize(2.0);
      glBegin(GL_POINTS);			
      for (long i = 0; i < mIsPinned.size(); i++) {
	glColor3f(1.0, 0.0, 0.0);
	if (mIsPinned[i]) {
	  glVertex3f(sp_->vertex(i)->coordinate()[0],
		     sp_->vertex(i)->coordinate()[1],
		     sp_->vertex(i)->coordinate()[2]);
	}
      }
      glEnd();
    }
		
    void reset_t(){
      mdt = 0.1;
    }
		
  private:
    bool mbegin;
    control_ptr		 sp_;
    T mdt,meps,mks,mkrs;
    vector<coordinate_type>	pos0;	//new velocity
    vector<coordinate_type>	pos1;	//new velocity		
    vector<coordinate_type>	vel1;	//new velocity
    vector<coordinate_type>	vel0;	//old velocity			
    vector<coordinate_type>	mFext;		//external forces		
    vector<coordinate_type>	mFint;		//external forces				
    vector<coordinate_type>	mF;			//external forces
    vector<coordinate_type>	mF0;			//external forces
		
    //pins
    vector<bool>			mIsPinned;
    vector<long>			pinID;
    vector<coordinate_type> mP0;		
		
    T			mM;					//mass matrix		
    vector<T>	mK;					//tangential stiffness matrix one for each edge
    vector<T>	mKr;				//tangential stiffness matrix from rest position
    vector<T>   ml0;				//initial lengths
    vector<T>   ml;				//initial lengths
  };
	
  //--------------------------------------------------------------------------------------
    
  template <typename T>
  class spring_non_linear{ M1_TYPEDEFS
  typedef coordinate_type							point_type ;
    typedef typename list<coordinate_type >::iterator	pl_iterator;
  public:
    spring_non_linear(control_ptr obj_in){
      obj_in->pack();
      mks	= 500.0;//for now uniform k spring_non_linear
      mdt	= 0.1;
      sp_ = obj_in;
      vertex_array vl = sp_->get_vertices();
            
      this->initialize();
    }
        
    void initialize(){
      mM = 50.0;
      pos0.resize(sp_->get_vertices().size());
      vel0.resize(sp_->get_vertices().size());
      vel1.resize(sp_->get_vertices().size());
      mF.resize(sp_->get_vertices().size());
      mFext.resize(sp_->get_vertices().size());
      mFint.resize(sp_->get_vertices().size());
			
      mIsPinned.resize(sp_->get_vertices().size());
      ml0.resize(sp_->get_edges().size());
      mP0.resize(sp_->get_vertices().size());
			
      mKr.resize(6*sp_->get_vertices().size());			
      mK.resize(6*sp_->get_edges().size());
      vector<edge_ptr> ev = sp_->get_edges();
			
      for (long i = 0; i < ev.size(); i++) {
	edge_ptr e = ev[i];
	coordinate_type c1 = e->v1()->coordinate();
	coordinate_type c2 = e->v2()->coordinate();
	T l = dist(c1,c2);
	ml0[i] = l;
      }
			
      //			for (long i = 0; i < mM.size(); i++) {
      //                mM = 1.0;
      //            }
			
      for (long i = 0; i < mIsPinned.size(); i++) {
	mIsPinned[i] = false;
      }
			
      vector<vertex_ptr> V = sp_->get_vertices();
      for (long i = 0; i < pos0.size(); i++) {				
	mP0[i] = V[i]->coordinate();
	pos0[i] = V[i]->coordinate();
      }
    }
		
		
    al::Vec<3,T>& get_force(long i){
      return mFint[i];
    }
		
    T &M()       {return mM;}
    T  M() const {return mM;}		
		
    void calculate_forces_and_stiffness(){
      edge_array E = sp_->get_edges();
      for (long i = 0; i < mFint.size(); i++) {
	mFint[i] = 0;
      }
			
      int nthreads, tid, i, chunk;
      chunk = 10;
#pragma omp parallel private(i,tid)	
      //#pragma omp parallel shared(b,R,xt,x,p,nthreads) private(i,tid)
      {
	tid = omp_get_thread_num();
	if (tid == 0)
	  {
	    nthreads = omp_get_num_threads();
	  }
#pragma omp for schedule(dynamic,chunk)
				
				
	for (long i = 0; i < E.size(); i++) {
	  edge_ptr e = E[i];
	  coordinate_type c1 = e->v1()->coordinate();
	  coordinate_type c2 = e->v2()->coordinate();
	  long v1 = e->v1()->position_in_set();
	  long v2 = e->v2()->position_in_set();
					
	  T l = dist(c1,c2);
	  T l0 = ml0[i];
					
	  T kd = mks/l0;
	  coordinate_type fs,fd;
	  fs = e->k*(c2-c1)*(l-l0)/l;				
	  //fd = kd*(vel0[v1]-vel0[v2])*(c1-c2)/l;				
	  //fd = 0;
					
	  mFint[v1] += mIsPinned[v1] ? 0 : (fs+fd);
	  mFint[v2] -= mIsPinned[v2] ? 0 : (fs+fd);
					
					
	  T dx = c2[0] - c1[0];
	  T dy = c2[1] - c1[1];
	  T dz = c2[2] - c1[2];				
					
	  T linv  = 1./l;
	  T
	    k00 = e->k*(-1. + l0*linv*(1. -(dx*dx*linv*linv))),
	    k11 = e->k*(-1. + l0*linv*(1. -(dy*dy*linv*linv))),
	    k22 = e->k*(-1. + l0*linv*(1. -(dz*dz*linv*linv))),				
	    k01 = -e->k*l0*linv*(dx*dy*linv*linv),
	    k02 = -e->k*l0*linv*(dx*dz*linv*linv),
	    k12 = -e->k*l0*linv*(dy*dz*linv*linv);                    
					
	  mK[6*i + 0] = k00;
	  mK[6*i + 1] = k11;
	  mK[6*i + 2] = k22;
	  mK[6*i + 3] = k01;
	  mK[6*i + 4] = k02;
	  mK[6*i + 5] = k12;
	}
      }
    }
		
    void set_pin(int i, bool isPinned){
      mIsPinned[i] = isPinned;
      pinID.push_back(i);
    }
		
    void update_pinned_positions(){
      vertex_array& V = sp_->get_vertices();
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
      //sp_->update_all();
    }
		
    virtual void explicite_euler(){
      this->calculate_forces_and_stiffness();
      vector<vertex_ptr> V = sp_->get_vertices();		
      for (long i = 0; i < vel1.size(); i++) {
	vel1[i] = mdt/mM*(mFint[i] + mFext[i]) + vel0[i];
	V[i]->coordinate() += (vel1[i])*mdt;				
      }
      swap(vel0, vel1);
    }
		
    virtual void verlet(){
			
      vector<vertex_ptr> V = sp_->get_vertices();		
      for (long i = 0; i < pos0.size(); i++) {				
	pos0[i] = V[i]->coordinate();
      }
      this->calculate_forces_and_stiffness();
      for (long i = 0; i < vel1.size(); i++) {
	V[i]->coordinate() = pos0[i] + vel0[i]*mdt + (mFint[i]+mFext[i])*mdt*mdt/mM;
	vel1[i] = (V[i]->coordinate() - pos0[i])/mdt;
      }
      swap(vel0, vel1);
			
    }
		
    virtual void implicite_euler(){			
			
      vector<vertex_ptr> V = sp_->get_vertices();
			
      for (long i = 0; i < V.size(); i++) {				
	pos0[i] = V[i]->coordinate();
      }
			
      int N = 0;
      while (N < 4) {
				
	this->calculate_forces_and_stiffness();				
	for(long i = 0; i < vel0.size(); i++){
	  mF[i] = mM*vel0[i] + mdt*mFint[i] + mdt*mFext[i];
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
			
      vector<edge_ptr> E = sp_->get_edges();		
      vector<vertex_ptr> V = sp_->get_vertices();
			
      int nthreads, tid, i, chunk;
      chunk = 10;
#pragma omp parallel private(i,tid)	
      //#pragma omp parallel shared(b,R,xt,x,p,nthreads) private(i,tid)
      {
	tid = omp_get_thread_num();
	if (tid == 0)
	  {
	    nthreads = omp_get_num_threads();
	  }
#pragma omp for schedule(dynamic,chunk)
				
                
	for (long i = 0; i < x.size(); i++) {
	  vertex_ptr v = V[i];
	  coordinate_type& cii = x[i];
					
	  if (!mIsPinned[i]){b[i] = cii*mM;}
	  else {b[i] = cii;}							
					
	  vertex_iterator<T> fvb = sp_->vertex(i)->begin();
	  vertex_iterator<T> fve = sp_->vertex(i)->end();
	  coordinate_type fadj = 0;
					
	  while (fvb != fve) {
	    edge_ptr e = fvb.edge();
	    int ei = e->position_in_set();					
	    long vi1 = e->v1()->position_in_set();
	    long vi2 = e->v2()->position_in_set();
						
	    long eij = (e->v1() == *fvb) ? vi2 : vi1;
						
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
	    fvb++;
	  }
	}
      }
    }
		
		
		
    void randomize(){
      vertex_array V = sp_->get_vertices();
      for (long i = 0; i < V.size(); i++) {
	V[i]->coordinate()[0] += randd((T)0.1);
	V[i]->coordinate()[1] += randd((T)0.1);
	V[i]->coordinate()[2] += randd((T)0.1);
	pos0[i] = V[i]->coordinate();
      }
    }
		
    void draw(){
      glDisable(GL_DEPTH_TEST);					// Enables Depth Testing 
      glEnable(GL_BLEND);					
      glLineWidth(2.0);
      glBegin(GL_LINES);			
      for (long i = 0; i < mFint.size(); i++) {
	glColor3f(0.25, 0.25, 0.25);
	glVertex3f(sp_->vertex(i)->coordinate()[0],
		   sp_->vertex(i)->coordinate()[1],
		   sp_->vertex(i)->coordinate()[2]);
	glVertex3f(sp_->vertex(i)->coordinate()[0] + mFint[i][0],
		   sp_->vertex(i)->coordinate()[1] + mFint[i][1],
		   sp_->vertex(i)->coordinate()[2] + mFint[i][2]);
      }
      glEnd();
			
      glPointSize(2.0);
      glBegin(GL_POINTS);			
      for (long i = 0; i < mIsPinned.size(); i++) {
	glColor3f(1.0, 0.0, 0.0);
	if (mIsPinned[i]) {
	  glVertex3f(sp_->vertex(i)->coordinate()[0],
		     sp_->vertex(i)->coordinate()[1],
		     sp_->vertex(i)->coordinate()[2]);
	}
      }
      glEnd();
    }
		
    void reset_t(){
      mdt = 0.1;
    }
		
  private:
    control_ptr		 sp_;
    T mdt,meps,mks,mkrs;
    vector<coordinate_type>	pos0;	//new velocity		
    vector<coordinate_type>	vel1;	//new velocity
    vector<coordinate_type>	vel0;	//old velocity			
    vector<coordinate_type>	mFext;	//external forces		
    vector<coordinate_type>	mFint;	//external forces				
    vector<coordinate_type>	mF;	//external forces
		
    //pins
    vector<bool>			mIsPinned;
    vector<long>			pinID;
    vector<coordinate_type> mP0;		
		
    T			mM;			//mass matrix		
    vector<T>	mK;				//tangential stiffness matrix one for each edge
    vector<T>	mKr;				//tangential stiffness matrix from rest position
    vector<T>   ml0;				//initial lengths
  };
    
	
  template <typename T>
  class spring_electric{ M1_TYPEDEFS
  typedef coordinate_type							point_type ;
    typedef typename list<coordinate_type >::iterator	pl_iterator;
  public:
    spring_electric(control_ptr obj_in){
      K_	= 0.1;
      C_	= 0.5;
      eps_= 0.01;
      dt_	= 0.1;
      sp_ = obj_in;
      vertex_array vl = sp_->get_vertices();
			
      energy0 = 10000000;
      energy  = 0;
      progress = 0;
			
      va_iterator itb = vl.begin();
      va_iterator ite = vl.end();
      while (itb != ite) {
	pos1.push_back((*itb)->coordinate());
	++itb;
      }
    }
    virtual void update_positions(){
      vertex_array vl = sp_->get_vertices();
      vector<coordinate_type> flist;
      va_iterator itb = vl.begin();
      va_iterator ite = vl.end();
      energy0 = energy; energy = 0;
      for (long i = 0; i < vl.size(); i++) {
				
	coordinate_type posf;
				
	vertex_ptr vi = vl[i];
	for (long j = 0; j < vl.size(); j++) {
	  vertex_ptr vj = vl[j];
	  if (i != j) {
	    coordinate_type dpos = (vj->coordinate() - vi->coordinate());
	    T dist = dpos.magSqr();
	    T inv1 = 1./(dist + eps_);
	    T f = -C_*K_*K_*inv1;
	    posf += dpos*f;
	  }								
	}
	flist.push_back(posf);
				
	//				edge_array va = vi->get_edges();
	//				for (long j = 0; j < va.size(); j++) {
	//					vertex_ptr vj = va[j]->other(vi);
	//					coordinate_type dpos = (vj->coordinate() - vi->coordinate());					
	//					
	//					T dist = dpos.mag();
	//					T f   = dist/K_;
	//					posf += dpos*f;
	//				}
				
	//				T magf = posf.mag();
	//				if (magf > 0.0001) {
	//					pos1[i]  = vi->coordinate() + posf/magf*dt_;
	//				}
	//				
	//				T magf2 = posf.magSqr();
	//				energy += magf2;
				
      }
			
      edge_array el = sp_->get_edges();
      for (long i = 0; i < el.size(); i++) {
	vertex_ptr vi = el[i]->v1();
	vertex_ptr vj = el[i]->v2();
	coordinate_type dpos = (vj->coordinate() - vi->coordinate());					
				
	T dist = dpos.mag();
	T f   = dist/K_;
	coordinate_type posf = dpos*f;
	flist[vi->position_in_set()] += posf;
	flist[vj->position_in_set()] -= posf;
      }
			
      for (long i = 0; i < vl.size(); i++) {
	coordinate_type posf = flist[i];
				
	T magf = posf.mag();
	T magf2 = posf.magSqr();
	energy += magf2;
	vl[i]->coordinate()+= posf/magf*dt_;
	if(magf > 0.0001) pos1[i] += posf/magf*dt_;
				
	coordinate_type &posi = vl[i]->coordinate();
	coordinate_type &pos1i = pos1[i];				
	//swap(posi, pos1i);
      }
			
      this->update_t();
    }
		
    void update_t(){
      T t = 0.90;
      if (energy<energy0) {
	progress += 1;
	if(progress>=5){
	  dt_ /= t;
	  progress = 0;
	}
      }
      else {
	progress = 0;
	dt_ *= t;
      }
			
			
    }
    void randomize(){
      vertex_array vl = sp_->get_vertices();
      va_iterator itb = vl.begin();
      va_iterator ite = vl.end();
      while (itb != ite) {
	(*itb)->coordinate()[0] = randd((T)0.1);
	(*itb)->coordinate()[1] = randd((T)0.1);
	(*itb)->coordinate()[2] = randd((T)0.1);
				
	itb++;
      }
    }
    void reset_t(){
      dt_ = 0.1;
    }
		
  private:
    control_ptr		 sp_;
    T dt_,eps_,C_,K_, energy0, energy;
    int progress;
    vector<coordinate_type>	pos1;  
  };
	
}
#endif
