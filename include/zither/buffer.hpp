/*
 *  buffer.h
 *  Phase Vocoder
 *
 *  Created by John Delaney on 4/14/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


#ifndef __BUFFER__
#define __BUFFER__
//#include <OpenGL/gl.h>				// Header File For The OpenGL32 Library
//#include <OpenGL/glu.h>			// Header File For The GLu32 Library
//#include <GLUT/glut.h>			// Header File For The GLUT Library
#include <iostream>
#include <math.h>

using namespace std;
namespace zither{

  template <class T, int S>
  class buffer {	
  protected:
    std::array<T,S> mArray;
	
  public:
    buffer(){}

    buffer(const buffer & in) : 
      mArray(in.mArray)
    {
    }
	
    ~buffer(){
    }
    
    size_t size(){return S;}
    
    T& operator[](const int i){  return mArray[i];}	
    T operator[](const int i) const{ return mArray[i];}
    T operator()(T x){
      int xi = (int)(std::ceil(x));
      int x0 = xi%S;
      int x1 = (x0 + 1)%S;		
      T y0 = mArray[x0];
      T y1 = mArray[x1];
      T C = (T)xi - x;
      T out = C*y0 + (1.0-C)*y1;
      return out;
    }

    void print() const{
		
      cout << "buffer: " << endl;
      for(int i=0; i < this->capacity(); i++){
        cout << bufferHead[i] << " ";
      }
      cout << endl;
    }
	

    //	void draw(){
    //		int bufferLength=this->capacity()-1;
    //		T lbuffer = 0; T cbuffer = 0;
    //		T ci, li = 0;
    //		
    //		glPushMatrix();
    //		glLoadIdentity();
    //		glMatrixMode(GL_PROJECTION);
    //		//gluOrtho2D ( 0.0, 1.0, 0.0, 1.0 );
    //		glPushMatrix();
    //		glLoadIdentity ();
    //
    //		glTranslatef(-0.85, -0.65, 0.0);
    //		glScalef(1.0, 1.0, 1.0);
    //
    //		glColor4f(1.0f, 1.0f, 1.0f, 0.5f);
    //		glBegin(GL_QUADS);
    //		glVertex3f( 0 , 0.25, 0.0f);
    //		glVertex3f( 1 , 0.25, 0.0f);
    //		glVertex3f( 1 , -0.25, 0.0f);
    //		glVertex3f( 0 , -0.25, 0.0f);
    //		glEnd();
    //		
    //		glColor4f(0.3f, 0.3f, 0.3f, 0.15f);
    //		for (int i=0; i<bufferLength; i+=(bufferLength/50) ){
    //			
    //			cbuffer = fabs(bufferHead[i])/5;
    //			ci = (T)i/(T)bufferLength;
    //			if (cbuffer > 1.0) cbuffer = 1.0;
    //			
    //			glBegin(GL_QUADS);
    //			glVertex3f( li , lbuffer,0.0f);
    //			glVertex3f( ci , cbuffer,0.0f);
    //			glVertex3f( ci ,-cbuffer,0.0f);
    //			glVertex3f( li ,-lbuffer,0.0f);
    //			glEnd();
    //			lbuffer = cbuffer; li = ci;
    //		}
    //		
    //		glPopMatrix();
    //		glMatrixMode(GL_MODELVIEW);
    //		glPopMatrix();
    //	}
  };
}//zither
#endif
