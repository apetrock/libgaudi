/*
 *  bezier.hpp
 *  Manifold
 *
 *  Created by John Delaney on 8/17/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <math.h>
#include <vector>
#include <iostream>
#include "nearest_point.hpp"
#include "geometry_types.hpp"

//#include "al_Functions.hpp" //pull some of that stuff here, and make it work with eigen
//#include "al_Quat.hpp"
//#include "al_Interpolation.hpp"
#ifndef __BEZIER__
#define __BEZIER__


template <class Tf, class Tv>
class slerp_interpolator{
public:
  virtual inline Tv operator()(Tf t, const Tv& p0, const Tv& p1){
    return slerp(p0,p1,t);
  }
};

template <class Tf, class Tv>
class lerp_interpolator{
public:
  virtual inline Tv operator()(Tf t, const Tv& p0, const Tv& p1){
    //return al::ipl::linear(t,p0,p1);
  }
};

template <int N, typename T>
class point_space {
public:
  typedef T										double_t;
  typedef Eigen::Matrix<T,N,1>						point_t;
  typedef lerp_interpolator<double_t,point_t>	interpolator_t;
  typedef unsigned int							uint;
  typedef unsigned long							ulong;
};

template <int N, typename T>
class rotor_space {
public:	
  typedef T										double_t;
  typedef Eigen::Quaternion<T>								point_t;
  typedef slerp_interpolator<double_t,point_t>	interpolator_t;
  typedef unsigned int							uint;
  typedef unsigned long							ulong;
};

template<typename SPACE> 
class bezier_curve{
public:
	
  typedef typename SPACE::double_t		double_t;
  typedef typename SPACE::point_t		point_t;
  typedef typename SPACE::interpolator_t	interpolator_t;
  typedef typename SPACE::uint			uint;
  typedef typename SPACE::ulong			ulong;
	
  bezier_curve(int chordsize){
    mChordSize = 1 << chordsize;
    mBaseKnotOrder = 4;
  }
  ~bezier_curve(){
  }
	
  int push_point(point_t& P){
    mControl.push_back(P);
    return mControl.size()-1;
  }
	
  void delete_point(ulong id){
    vector<point_t> nControl;
    for (ulong i = 0; i < mControl.size(); i++) {
      if (i != id) {
	nControl.push_back(mControl[i]);
      }
    }
    mControl = nControl;
  }
	
  void update_knot_vector(){
    mKnots.clear();
    this->calc_cubic_knot(mControl, mKnots);
  }

	
  void calc_knot(vector<point_t>& P, vector<point_t>& U){
    uint ks = mBaseKnotOrder;
    for(ulong i = 0; i < P.size()-1; i += (ks-1)){
      vector<point_t> pi;
      if (P.size() - i < ks) {
	ks = P.size() - i;
      }
      for(ulong j = i; j < i+ks; j++){
	pi.push_back(P[j]);
      }
      calc_chord(pi, U);
    }
  }
	
  void calc_chord(vector<point_t>& P, vector<point_t>& U){
    long s = mChordSize;
    vector<point_t> Q = P;		
    for(ulong c = 0; c < s; c++){
      ulong n = P.size();
      double_t t = double_t(c)/double_t(s);
      for(ulong k = 1; k < n; k++){					
	for(ulong i = 0; i < n-k; i++){
	  Q[i] = (1-t)*Q[i] + t*Q[i+1];
	  point_t q = Q[i];						
	}
      }
      U.push_back(Q[0]);
    }
  }

  void calc_cubic_knot(vector<point_t>& P, vector<point_t>& U){
    uint ks = mBaseKnotOrder;
    for(ulong i = 0; i < P.size()-1; i += 3){
      calc_cubic_chord(i,P, U);
    }
  }
	
  void calc_cubic_chord(uint beg, vector<point_t>& P, vector<point_t>& U){
    ulong s = mChordSize;
    vector<point_t> Q = P;
    double_t t = 0; double_t dt = 1./double_t(s);
    point_t P0 = P[beg];
    point_t P1 = P[beg+1];
    point_t P2 = P[beg+2];
    point_t P3 = P[beg+3];
    interpolator_t lerp;
    for(ulong i = 0; i < s+1; i++){
      point_t b01  = lerp(t,P0,P1);
      point_t b12  = lerp(t,P1,P2);
      point_t b23  = lerp(t,P2,P3);			
      point_t b012 = lerp(t,b01,b12);
      point_t b123 = lerp(t,b12,b23);
      point_t b0123 = lerp(t,b012,b123);			
      U.push_back(b0123);
      t+=dt;
    }
  }
	
  point_t find_nearest(point_t pi, double_t& t, int& knot){
    //cubic knots only;
    double_t tn;
    int ks = mBaseKnotOrder; double_t dist, ndist; point_t c;
    dist = 1e8;
    vector<point_t>& P = mControl;
    for(ulong i = 0; i < P.size()-1; i += (ks-1)){
      point_t* Pi = new point_t[ks];
      if (P.size() - i < ks) {
	ks = P.size() - i;
      }
      int is = 0;
      for(ulong j = i; j < i+ks; j++){
	Pi[is] = P[j]; is++;
      }
      point_t ci = find_nearest_chord(pi,Pi,ndist, tn);
      if (ndist < dist) {
	dist = ndist; c = ci; t = tn; knot = i;
      }
    }
    return c;
  }
	
  point_t find_nearest_chord(point_t Pi, point_t *P, double_t&dist, double_t& t){
    bezier_nearest_point<double_t> np;
    return np.NearestPointOnCurve(Pi, P, dist,t);
  }
	
  void insert(int knot, double_t t){
    vector<point_t> ct,cn;
    for (long i = knot; i < knot+mBaseKnotOrder; i++) {
      ct.push_back(mControl[i]);
    }		
    this->subdivide_knot(ct, t);
    cout << "knot" << knot << endl;
    cout << "original" << endl;
    for(long i = 0; i < mControl.size(); i++){
      cout << i << ": " << mControl[i][0] << " "<< mControl[i][1] << " ";
      cout << endl;
    }
    cout << endl;		
    for (long i = 0; i < knot; i++) {
      cn.push_back(mControl[i]);
    }
    for (long i = 0; i < ct.size(); i++) {
      cn.push_back(ct[i]);
    }
    for (long i = knot+mBaseKnotOrder; i < mControl.size(); i++) {
      cn.push_back(mControl[i]);
    }
    cout << "insert" << endl;
    for(long i = 0; i < ct.size(); i++){
      cout << i << ": " << ct[i][0] << " "<< ct[i][1] << " ";
      cout << endl;
    }
    cout << endl;
    cout << "inserted" << endl;
    for(long i = 0; i < cn.size(); i++){
      cout << i << ": " << cn[i][0] << " "<< cn[i][1] << " ";
      cout << endl;
    }
    mControl = cn;
  }
	
  void subdivide_knot(vector<point_t >& Po, double_t t){
    vector<point_t > Q = Po;
    vector<point_t > Pn;
    long n = Po.size();
    long ip = 1, kp = n;
    Pn.push_back(Po.back());			
    for(long k = 1; k < n; k++){
      for(long i = 0; i < n-k; i++){
	Q[i] = (1-t)*Q[i] + t*Q[i+1];
	point_t q = Q[i];						
      }
      point_t t = Q[0];
      Po[ip] = Q[0];
      Pn.push_back(Q[n-k]);
      ip++;
    }
		
    for(long i = Pn.size()-1; i > 0; i--){
      Po.push_back(Pn[i]);
    }
		
    //nQ[nN - 1] = P[N];
  }
  int mChordSize;
  int mBaseKnotOrder;
  vector<point_t>		mKnots;	
  vector<point_t>		mControl;
  vector<double_t>	mDistTree;		
};

#endif
