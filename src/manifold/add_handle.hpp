/*
 *  add_handle.hpp
 *  libGaudi
 *
 *  Created by John Delaney on 10/8/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
#include <iostream>
//#include "al_Quat.hpp"
//#include "al_Vec.hpp"
#include "bezier.hpp"
#include "graph_skinning.hpp"

template <typename T>
Eigen::Matrix<T,4,1> rotate_to_plane(Eigen::Matrix<T,4,1> normP, 
				     Eigen::Matrix<T,4,1> norm, 
				     Eigen::Matrix<T,4,1> cen, 
				     Eigen::Matrix<T,4,1> vec){
  norm.normalize();
  normP.normalize();	
  Eigen::Matrix<T,4,1> xe = cross(norm,normP).normalize();
  T r1 = norm.mag();
  T rP = normP.mag();
  T o = acos(dot(normP,norm)/r1/rP) - M_PI;
  T o0 = acos(dot(normP,norm)/r1/rP);
  Eigen::Quaternion<T> qP(normP);
  Eigen::Quaternion<T> qi(norm);
  //	Eigen::Quaternion<T>  dq = (qP-qi).normalize();
  Eigen::Quaternion<T> dq(cos(o/2.), sin(o/2.)*xe[0], sin(o/2.)*xe[1], sin(o/2.)*xe[2]);		
  //	Eigen::Quaternion<T> dq = qP*qi.conj()*qP;
  Eigen::Matrix<T,4,1> c = dq.normalize().rotate(vec - cen);
  return c;
}

template <typename T>
Eigen::Matrix<T,4,1> rotate_to_xy(Eigen::Matrix<T,4,1> norm, 
				  Eigen::Matrix<T,4,1> cen, 
				  Eigen::Matrix<T,4,1> vec){
  Eigen::Matrix<T,4,1> vz(0,0,1);
  norm.normalize();
  Eigen::Matrix<T,4,1> xe = cross(norm,vz).normalize();
  T r1 = norm.mag();
  T o = acos(dot(vz,norm)/r1);
  Eigen::Quaternion<T> dq(cos(o/2.), sin(o/2.)*xe[0], sin(o/2.)*xe[1], sin(o/2.)*xe[2]);	
  Eigen::Matrix<T,4,1> c = dq.rotate(vec - cen);
  return c;
}

template <typename T>
Eigen::Matrix<T,4,1> rotate_from_xy(Eigen::Matrix<T,4,1> norm, 
				    Eigen::Matrix<T,4,1> cen, 
				    Eigen::Matrix<T,4,1> vec){
  Eigen::Matrix<T,4,1> vz(0,0,1);
  norm.normalize();
  Eigen::Matrix<T,4,1> xe = -cross(norm,vz).normalize();
  T r1 = norm.mag();
  T o = acos(dot(vz,norm)/r1);
  Eigen::Quaternion<T> dq(cos(o/2.), sin(o/2.)*xe[0], sin(o/2.)*xe[1], sin(o/2.)*xe[2]);	
  Eigen::Matrix<T,4,1> c = dq.rotate(vec);
  c += cen;
  return c;
}

template<typename T>
Eigen::Quaternion<T> naive_slerp(const Eigen::Quaternion<T>& input, const Eigen::Quaternion<T>& target, T amt){
  Eigen::Quaternion<T> result;
	
  if (amt==T(0)) {
    return input;
  } else if (amt==T(1)) {
    return target;
  }
	
  int bflip = 0;
  T dot_prod = input.dot(target);
  T a, b;
	
  //clamp
  dot_prod = (dot_prod < -1) ? -1 : ((dot_prod > 1) ? 1 : dot_prod);
	
  // if B is on opposite hemisphere from A, use -B instead
  if (dot_prod < 0.0) {
    dot_prod = -dot_prod;
    bflip = 1;
  }
	
  T cos_angle = acos(dot_prod);
  if(ABS(cos_angle) > QUAT_EPSILON) {
    T sine = sin(cos_angle);
    T inv_sine = 1./sine;
		
    a = sin(cos_angle*(1.-amt)) * inv_sine;
    b = sin(cos_angle*amt) * inv_sine;
		
    if (bflip) { b = -b; }
  } else {
    // nearly the same;
    // approximate without trigonometry
    a = amt;
    b = 1.-amt;
  }
	
  result.w = a*input.w + b*target.w;
  result.x = a*input.x + b*target.x;
  result.y = a*input.y + b*target.y;
  result.z = a*input.z + b*target.z;
	
  result.normalize();
  return result;
}


template <typename T>
Eigen::Matrix<T,4,1> blend(T t, Eigen::Matrix<T,4,1> p00, Eigen::Matrix<T,4,1> pNorm,
		   Eigen::Matrix<T,4,1> p10){
  T r0 = p00.mag();
  T r1 = p10.mag();
  T rs = linear(t,r0,r1);
  T dotp = dot(p00,p10);
  dotp = abs(dotp);
  T theta = acos(dotp/r0/r1);
  Eigen::Quaternion<T> dq0;
  p00.normalize();
  p10.normalize();
  Eigen::Quaternion<T> q00(p00);
  Eigen::Quaternion<T> q10(p10);	
  if (theta < 0.01) {
    //just give it a little nudge in the right direction
    T o = 0.10;
    Eigen::Quaternion<T> dq(cos(o/2.), sin(o/2.)*pNorm[0], sin(o/2.)*pNorm[1], sin(o/2.)*pNorm[2]);	
    Eigen::Matrix<T,4,1> p01 = dq.rotate(p00);
    Eigen::Quaternion<T> q01(p01);		
    dq0 = Eigen::Quaternion<T>::slerp(q01,q10,t);	
  } 
  else {
    dq0 = Eigen::Quaternion<T>::slerp(q00,q10,t);	
  }	
  T dott = p00.dot(p10);
  T dotq0 = q00.dot(q10);	
  T dotq1 = dq0.dot(q10);
  return Eigen::Matrix<T,4,1>(rs*dq0[1],rs*dq0[2],rs*dq0[3]);	
  //return Eigen::Matrix<T,4,1>(p01);
}

//template <typename T>
//Eigen::Matrix<T,4,1> blend(T t, Eigen::Matrix<T,4,1> p0, Eigen::Matrix<T,4,1> p1){
//	T tol = 0.00001;
//
//	T r0 = p0.mag();
//	T r1 = p1.mag();

//	if ((p0-p1).magSqr() < tol) {
//		return (1.-t)*p0 + t*p1;
//	}
//	return sin((1.-t)*theta)/sin(theta)*p0 + sin(t*theta)/sin(theta)*p1;
//}



namespace m2 {
  template <typename SPACE>
  class add_handle{
    M2_TYPEDEFS;
  public:
    bool add(control_ptr ob, long i, long j, T teni, T tenj){
      //	m2::subdivide<T> sub;
      //	m2Ch = sub.subdivide_control(m2Ch);
      //bezier_curve<point_space<3,T> > mCurve(3);
      face_ptr fi = ob->face(i);
      face_ptr fj = ob->face(j);
			
      bezier_curve<point_space<3,typename SPACE::type> > mCurve(4);	
      coordinate_type c0 = fi->calculate_center();	
      coordinate_type c0t = c0 + teni*ob->face(i)->normal();
      coordinate_type c1 = fj->calculate_center();	
      coordinate_type c1t = c1 + tenj*ob->face(j)->normal();			
			
      if (fi->size() != fj->size() ) {
	return false;
      }
			
      mCurve.push_point(c0);
      mCurve.push_point(c0t);	
      mCurve.push_point(c1t);		
      mCurve.push_point(c1);
			
      mCurve.update_knot_vector();
			
			
      Eigen::Quaternion<T>  q1(ob->face(j)->normal());q1.normalize();			
			
      vector<coordinate_type > f0 = fi->coordinate_trace();
      vector<coordinate_type > f1 = fj->coordinate_trace();
			
      long Ni = mCurve.mKnots.size();
      T dt = T(1./T(Ni-1));
      T t = dt*0.5;
      //						for(long ii = 0; ii < 3; ++ii){			
      for(long ii = 0; ii < Ni-1; ++ii){
	coordinate_type Ni = fi->normal();
	Eigen::Quaternion<T>  q0(Ni); q0.normalize();
	coordinate_type de = mCurve.mKnots[ii]-mCurve.mKnots[ii+1];
	de.normalize();
	Eigen::Quaternion<T> qi(de); qi.normalize();
	//
	coordinate_type xe = cross(Ni,-de).normalize();
	T r0 = de.mag();
	T r1 = Ni.mag();
	T o = acos(dot(-de,Ni)/r0/r1);
				
	Eigen::Quaternion<T> dq(cos(o/2.), sin(o/2.)*xe[0], sin(o/2.)*xe[1], sin(o/2.)*xe[2]);
	//				Eigen::Quaternion<T> dq = (qi*q0); //dq[0]*=M_PI;
	//dq.normalize();			
	coordinate_type cen = fi->calculate_center();
	face_vertex_ptr itb = fi->fbegin();
	face_vertex_ptr ite = fi->fend();
				
	construct<SPACE> bevel;
	bevel.bevel_face(ob, fi, 0.0, 0.0);			
				
	bool at_head = false;
	while (!at_head) {
	  at_head = itb == ite;
	  coordinate_type c = itb->coordinate();
	  coordinate_type cp0 = dq.rotate(c-cen);
	  itb->coordinate() = cp0 + (mCurve.mKnots[ii+1]);
	  itb = itb->next();					
	}
				
	fi->update_normal();
	fi->update_center();				
	this->slerp_face(t, fi, fj);				
	//				cen = fi->calculate_center();
	//				itb = fi->fbegin();
	//				at_head = false;
	//				while (!at_head) {
	//					at_head = itb == ite;
	//					coordinate_type c = itb->coordinate();
	//					coordinate_type N = fi->normal();
	//					T o = 0.5*t*M_PI;
	//					Eigen::Quaternion<T> qb(cos(o*0.5), sin(o*0.5)*N[0], sin(o*0.5)*N[1], sin(o*0.5)*N[2]);
	//					qb.normalize();
	//					Eigen::Matrix<T,4,1> cp1 = qb.rotate(c-cen);
	//					//					Eigen::Matrix<T,4,1> cp1 = cp0;
	//					//					itb->coordinate() = cp0 + (mCurve.mKnots[ii] + mCurve.mKnots[ii + 1])*0.5;
	//					itb->coordinate() = cp1 + (mCurve.mKnots[ii+1]);
	//					itb = itb->next();					
	//				}
				
	t += dt;
      }
      graph_skinning<SPACE> gs; //TODO: stitch_faces needs to get moved to a surgery object
      gs.stitch_faces(ob, fi, fj, 0.01);
      return true;
    }
    bool slerp_face(T t, face_ptr fi, face_ptr fj){
      coordinate_type Ni = fi->normal();
      coordinate_type Nj = fj->normal();
      coordinate_type ceni = fi->calculate_center();
      coordinate_type cenj = fj->calculate_center();			
      face_vertex_ptr itb = fi->fbegin(), ite = fi->fend();
      face_vertex_ptr jtb = fj->fbegin(), jte = fj->fend();
      bool at_head = false;
      while (!at_head) {
	at_head = itb == ite;
				
	coordinate_type ci = itb->coordinate()-ceni;
	coordinate_type cip = itb->next()->coordinate()-ceni;				
	coordinate_type cj = rotate_to_plane(Ni,Nj,cenj,jtb->coordinate());
	//mOrient.push_back(ci + ceni);
	//mOrient.push_back(cj + ceni);
	coordinate_type cij = blend(t,ci,Ni,cj);				
	itb->coordinate() = cij + ceni;
				
	//				coordinate_type ci = rotate_to_xy(Ni,ceni,itb->coordinate()); 
	//				coordinate_type cip = rotate_to_xy(Ni,ceni,itb->next()->coordinate()); 				
	//				coordinate_type cj = rotate_to_xy(Nj,cenj,jtb->coordinate()); 
	//				coordinate_type cjp = rotate_to_xy(Nj,cenj,jtb->next()->coordinate()); 								
	//
	//				mOrient.push_back(ci + ceni);
	//				mOrient.push_back(cj + ceni);
	//				
	//				coordinate_type cij = blend(t,ci,Eigen::Matrix<T,4,1>(0,0,1),cj);				
	//				itb->coordinate() = rotate_from_xy(Ni,ceni,cij);
				
				
	itb = itb->next();
	jtb = jtb->prev();				
      }
      fi->update_normal();
      fi->update_center();
      cout << endl;
      cout << "---------" << endl;
      return true;
			
    }
  };	
}
