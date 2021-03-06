/********************************************************/
/* AABB-triangle overlap test code                      */
/* by Tomas Akenine-Möller                              */
/* Function: int triBoxOverlap(float boxcenter[3],      */
/*          float boxhalfsize[3],float triverts[3][3]); */
/* History:                                             */
/*   2001-03-05: released the code in its first version */
/*   2001-06-18: changed the order of the tests, faster */
/*                                                      */
/* Acknowledgement: Many thanks to Pierre Terdiman for  */
/* suggestions and discussions on how to optimize code. */
/* Thanks to David Hunt for finding a ">="-bug!         */
/********************************************************/
#ifndef __TRIBOXX_TEST__
#define __TRIBOXX_TEST__

#include <math.h>
#include <stdio.h>

//#include "vec_addendum.h"
//#include "allocore/al_Allocore.hpp"

template <typename T, typename CTYPE>
class tri_box{

#define X 0
#define Y 1
#define Z 2

#define FINDMINMAX(x0,x1,x2,min,max)		\
  min = max = x0;				\
  if(x1<min) min=x1;				\
  if(x1>max) max=x1;				\
  if(x2<min) min=x2;				\
  if(x2>max) max=x2;
public:
  int planeBoxOverlap(CTYPE normal,CTYPE vert, CTYPE maxbox)
  {
    int q;
    CTYPE vmin(0,0,0,0),vmax(0,0,0,0);
    float v;
    for(q=X;q<=Z;q++)
      {
        v=vert[q];
	if(normal[q]>0.0)
	  {
	    vmin[q]=-maxbox[q]-v;
	    vmax[q]=maxbox[q]-v;
	  }
	else
	  {
	    vmin[q]=maxbox[q]-v;
	    vmax[q]=-maxbox[q]-v;
	  }
      }
    //std::cout << " min/max: " << vmin.transpose() << " " << vmax.transpose() << std::endl;
    //std::cout << " normal: " << normal.transpose() << std::endl;
    
    if(dot(normal,vmin) > 0.0) return false;
    if(dot(normal,vmax) >=0.0) return true;
	
    return false;
  }


  /*======================== X-tests ========================*/
#define AXISTEST_X01(a, b, fa, fb)			\
  p0 = a*v0[Y] - b*v0[Z];				\
  p2 = a*v2[Y] - b*v2[Z];				\
  if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;}	\
  rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];	\
  if(min>rad || max<-rad) return false;

#define AXISTEST_X2(a, b, fa, fb)			\
  p0 = a*v0[Y] - b*v0[Z];				\
  p1 = a*v1[Y] - b*v1[Z];				\
  if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;}	\
  rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];	\
  if(min>rad || max<-rad) return false;

  /*======================== Y-tests ========================*/
#define AXISTEST_Y02(a, b, fa, fb)			\
  p0 = -a*v0[X] + b*v0[Z];				\
  p2 = -a*v2[X] + b*v2[Z];				\
  if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;}	\
  rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];	\
  if(min>rad || max<-rad) return false;

#define AXISTEST_Y1(a, b, fa, fb)			\
  p0 = -a*v0[X] + b*v0[Z];				\
  p1 = -a*v1[X] + b*v1[Z];				\
  if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;}	\
  rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];	\
  if(min>rad || max<-rad) return false;

  /*======================== Z-tests ========================*/

#define AXISTEST_Z12(a, b, fa, fb)			\
  p1 = a*v1[X] - b*v1[Y];				\
  p2 = a*v2[X] - b*v2[Y];				\
  if(p2<p1) {min=p2; max=p1;} else {min=p1; max=p2;}	\
  rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];	\
  if(min>rad || max<-rad) return false;

#define AXISTEST_Z0(a, b, fa, fb)			\
  p0 = a*v0[X] - b*v0[Y];				\
  p1 = a*v1[X] - b*v1[Y];				\
  if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;}	\
  rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];	\
  if(min>rad || max<-rad) return false;

  bool triBoxOverlap(const CTYPE boxcenter,
		     const CTYPE boxhalfsize,
		     const CTYPE triverts[3])
  {
    /*    use separating axis theorem to test overlap between triangle and box */
    /*    need to test for overlap in these directions: */
    /*    1) the {x,y,z}-directions (actually, since we use the AABB of the triangle */
    /*       we do not even need to test these) */
    /*    2) normal of the triangle */
    /*    3) crossproduct(edge from tri, {x,y,z}-directin) */
    /*       this gives 3x3=9 more tests */
    CTYPE v0,v1,v2;
    CTYPE normal,e0,e1,e2;
    CTYPE axis;
    T min,max,d,p0,p1,p2,rad,fex,fey,fez;  
	
    /* This is the fastest branch on Sun */
    /* move everything so that the boxcenter is in (0,0,0) */
    v0 = triverts[0] - boxcenter;
    v1 = triverts[1] - boxcenter;
    v2 = triverts[2] - boxcenter;
   
    /* compute triangle edges */
    e0 = v1 - v0;      /* tri edge 0 */
    e1 = v2 - v1;      /* tri edge 1 */
    e2 = v0 - v2;      /* tri edge 2 */
	
    /* Bullet 3:  */
    /*  test the 9 tests first (this was faster) */
    fex = fabs(e0[X]);
    fey = fabs(e0[Y]);
    fez = fabs(e0[Z]);
    AXISTEST_X01(e0[Z], e0[Y], fez, fey);
    AXISTEST_Y02(e0[Z], e0[X], fez, fex);
    AXISTEST_Z12(e0[Y], e0[X], fey, fex);
	
    fex = fabs(e1[X]);
    fey = fabs(e1[Y]);
    fez = fabs(e1[Z]);
    AXISTEST_X01(e1[Z], e1[Y], fez, fey);
    AXISTEST_Y02(e1[Z], e1[X], fez, fex);
    AXISTEST_Z0(e1[Y], e1[X], fey, fex);
	
    fex = fabs(e2[X]);
    fey = fabs(e2[Y]);
    fez = fabs(e2[Z]);
    AXISTEST_X2(e2[Z], e2[Y], fez, fey);
    AXISTEST_Y1(e2[Z], e2[X], fez, fex);
    AXISTEST_Z12(e2[Y], e2[X], fey, fex);
	
    /* Bullet 1: */
    /*  first test overlap in the {x,y,z}-directions */
    /*  find min, max of the triangle each direction, and test for overlap in */
    /*  that direction -- this is equivalent to testing a minimal AABB around */
    /*  the triangle against the AABB */
	
    /* test in X-direction */
    FINDMINMAX(v0[X],v1[X],v2[X],min,max);
    if(min>boxhalfsize[X] || max<-boxhalfsize[X]) return false;
	
    /* test in Y-direction */
    FINDMINMAX(v0[Y],v1[Y],v2[Y],min,max);
    if(min>boxhalfsize[Y] || max<-boxhalfsize[Y]) return false;
	
    /* test in Z-direction */
    FINDMINMAX(v0[Z],v1[Z],v2[Z],min,max);
    if(min>boxhalfsize[Z] || max<-boxhalfsize[Z]) return false;
	
    /* Bullet 2: */
    /*  test if the box intersects the plane of the triangle */
    /*  compute plane equation of triangle: normal*x+d=0 */
    //	normal = ( e0.y*e1.z - e0.z*e1.y, e0.z*e1.x - e0.x*e1.z, e0.x*e1.y - e0.y*e1.x );
    normal = cross(e0,e1);
    d = dot(normal,v0);  /* plane eq: normal.x+d=0 */
    //std::cout << " norm, v0: "  << normal.transpose() << " " << v0.transpose() << std::endl;
    if(!planeBoxOverlap(normal,v0,boxhalfsize)) return false;
	
    return true;   /* box and triangle overlaps */
  }
#undef X
#undef Y
#undef Z
#undef FINDMINMAX
#undef AXISTEST_X01
#undef AXISTEST_X2
#undef AXISTEST_Y02
#undef AXISTEST_Y1
#undef AXISTEST_Z12
#undef AXISTEST_Z0
};

#endif
