/*
 *  control.hpp
 *  Phase Vocoder
 *
 *  Created by John Delaney on 12/16/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef __TWOMANIFOLDMAKE__
#define __TWOMANIFOLDMAKE__

#include "m2Includes.h"
#include "m2Operators.h"
#include <time.h>

//functions-------------------------------

namespace m2 {
  template <typename SPACE>
  class make  {
    //class designed to produce primitives, such as convex hulls or loading .obj, etc;		
    M2_TYPEDEFS
    public:
    make(){};
    ~make(){};
        
    control_ptr cube(T x, T y, T z){
      control_ptr tcube = new control_type();
      construct<SPACE>	cons;
      subdivide<SPACE>	subd;
      tcube->insert_vertex(-x*0.5,-y*0.5,-z*0.5);
      tcube->insert_vertex( x*0.5,-y*0.5,-z*0.5);
      tcube->insert_vertex( x*0.5, y*0.5,-z*0.5);
      tcube->insert_vertex(-x*0.5, y*0.5,-z*0.5);
            
      cons.insert_edge(tcube, tcube->vertex(0), tcube->vertex(1));
      cons.insert_edge(tcube, tcube->vertex(1), tcube->vertex(2));
      cons.insert_edge(tcube, tcube->vertex(2), tcube->vertex(3));
      cons.insert_edge(tcube, tcube->vertex(3), tcube->vertex(0));            
      cons.bevel(tcube, 0.0, z, 0.0);            
      return tcube;
    }

    control_ptr extruded_knot(vector<coordinate_type> & knot, T r){
      control_ptr cyl = new control_type();
      construct<SPACE>	cons;
      geometry_helper<SPACE> rotator;
      int N = 64;
      int kN = knot.size();
      coordinate_type t0, t1, tN, n, b;
      tN = rotator.periodicTangent(knot,kN-1);
      t0 = rotator.periodicTangent(knot,0);
      t1 = rotator.periodicTangent(knot,1);
      n = cross(tN,t1); n.normalize();	
      b = cross(n, t0); b.normalize();
      std::cout << n << " " << b << std::endl;
      T dThet = 2.0/(T)N*M_PI;
      T thet = 2.0*M_PI;
      for(int i = 0; i < N; i++){
	T x = knot[0][0] + r*(b[0]*cos(thet) + n[0]*sin(thet));
	T y = knot[0][1] + r*(b[1]*cos(thet) + n[1]*sin(thet));
	T z = knot[0][2] + r*(b[2]*cos(thet) + n[2]*sin(thet));
	cyl->insert_vertex(x,y, z);
	thet -= dThet;
      }

      for(int i = 0; i < N; i++){
	int ip = (i+1)%N;
	  cons.insert_edge(cyl, cyl->vertex(i), cyl->vertex(ip));
      }

      face_ptr fi = cyl->get_faces()[N];
      face_ptr f0 = cyl->get_faces()[0];
      for(int i = 0; i < knot.size() - 1; i++){
	cons.bevel_face(cyl,fi, 0.0, 0.0);
	face_vertex_ptr fvb= fi->fbegin();
	face_vertex_ptr fve= fi->fend();
	bool iterating = true;

	//q0 = rotator.parallelTransport(knot,q0,i);
	quat q0 = rotator.getPeriodicRotation(knot,i);

	coordinate_type cen = fi->calculate_center();
	while(iterating){
	  iterating = fvb != fve;
	  coordinate_type ci = fvb->vertex()->coordinate();
	  coordinate_type dc = ci - cen;
	  dc = q0.rotate(dc);
	  coordinate_type cp = dc + cen;

	  int im1 = i-1; im1 = im1 < 0 ? im1 + kN : im1;
	  int ip1 = i+1; ip1 = ip1 > kN-1 ? ip1 - kN : ip1;
	  fvb->vertex()->coordinate() = cp + 0.5*(knot[ip1] - knot[im1]);
	  fvb = fvb->next();
	}
      }
      cons.pipe_face(cyl, fi, f0);
      cyl->pack();
      return cyl;
    }

    control_ptr trefoil(T r){

      // The (p,q)-torus knot can be given by the parametrization
      // x = r*cos(p*phi)
      // y = r*sin(p*phi)
      // z = -sin(q*phi)
      // where r = cos(q*phi)+2 and 0 < phi< 2*pi. This lies on the surface of the
      // torus given by (r-2)^2 + z^2 = 1 (in cylindrical coordinates).

      vector<coordinate_type> knot;
      int knotN = 128;
 
      T dt = 2.0*M_PI/(T)knotN;
      T t = 0;
      for(int i = 0; i < knotN; i++){
	coordinate_type ki(sin(t) + 2.0*sin(2.0*t),
			   cos(t) - 2.0*cos(2.0*t),
			   -sin(3.0*t));
	knot.push_back(2.0*r*ki);
	t+=dt;
      }
      return extruded_knot(knot, r);
    }

    control_ptr torus(T r0, T r1){
      vector<coordinate_type> knot;
      int knotN = 64;
 
      T dt = 2.0*M_PI/(T)knotN;
      T t = 0;
      for(int i = 0; i < knotN; i++){
	coordinate_type ki(sin(t),cos(t),0);
	knot.push_back(r1*ki);
	t+=dt;
      }
      return extruded_knot(knot, r0);
    }

    control_ptr cylinder(T h, T r){
      control_ptr cyl = new control_type();
      construct<SPACE>	cons;
      subdivide<SPACE>	subd;
      
      int N = 32; 
      T dThet = 2.0/(T)N*M_PI;
      T thet = 0;
      for(int i = 0; i < N; i++){
	T x = r*cos(thet);
	T z = r*sin(thet);

	cyl->insert_vertex(x,-h*0.5, z);
	thet += dThet;
      }
      for(int i = 0; i < N; i++){
	int ip = (i+1)%N;
	  cons.insert_edge(cyl, cyl->vertex(i), cyl->vertex(ip));
      }

      cons.bevel(cyl,N, 0.25*h, 0.);
      cons.bevel(cyl,N, 0.25*h, 0.);
      cons.bevel(cyl,N, 0.25*h, 0.);
      cons.bevel(cyl,N, 0.25*h, 0.);

      cons.bevel(cyl,0, 0.0, -0.25*h);
      cons.bevel(cyl,0, 0.0, -0.25*h);
      cons.bevel(cyl,0, 0.0, -0.25*h);
      cons.bevel(cyl,0, 0.0, -0.25*h);

      cons.bevel(cyl,N, 0.0, -0.25*h);
      cons.bevel(cyl,N, 0.0, -0.25*h);
      cons.bevel(cyl,N, 0.0, -0.25*h);
      cons.bevel(cyl,N, 0.0, -0.25*h);

      cyl->pack();
      return cyl;
    }
                		
    control_ptr max_tet(std::vector<coordinate_type> clist){
      control_ptr ttet = new control_type();
      construct<SPACE>	cons;
      coordinate_type	avg = calc_avg(clist);
      T dmax = 0;
      coordinate_type vmax;
      for (int i = 0; i < clist.size(); i++) {
	T d = dist(clist[i], avg);
	if (d>dmax) {
	  vmax = clist[i];
	}
      }
      return ttet;
    }
				
  };
}; //end m2
#endif
