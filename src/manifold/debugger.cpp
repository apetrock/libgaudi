/*
 *  manifold_singleton.cpp
 *  Phase Vocoder
 *
 *  Created by John Delaney on 12/29/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

//#include <gl/glut.h>
#include "manifold/debugger.h"

namespace m2 {

  //////////////////////////
  //debugger
  //////////////////////////

  bool      Debugger::instance_flag = false;
  Debugger* Debugger::global_instance = NULL;

  Debugger& Debugger::get_instance(){
    static Debugger* debug;
    if (!debug->initialized()) {
      debug = new Debugger();
      debug->labelId = 0;
      debug->initialized() = true;
    }
    return *debug;
  }

  void Debugger::add_label(float x, float y, float z, string label){
    DebugLabels.push_back(label);
    labelId++;
    DebugLabelPoints.push_back(coordinate_type(x,y,z,1));
  }

  void Debugger::draw_labels(){
    glPointSize(5.0);

    for(int i = 0; i < DebugLabelPoints.size(); i++){
      coordinate_type c1 = DebugLabelPoints[i];

      glBegin(GL_POINTS);
      glColor3f(1.0,0.0,0.5);
      glVertex3d(c1[0],c1[1],c1[2]);
      glEnd();

      char* str = new char[DebugLabels[i].size()];
      for(int j = 0; j < DebugLabels[i].size(); j++){
	str[j] = DebugLabels[i][j];
      }
      renderString(c1[0],c1[1],c1[2],0.00005,str);
      delete str;
    }
    glPointSize(1.0);

  }

  void Debugger::add_id(float x, float y, float z){
    DebugIds.push_back(labelId);
    DebugIdPoints.push_back(coordinate_type(x,y,z,1));
    labelId++;
  }

  void Debugger::draw_ids(){

    glColor3f(1.0,0.0,0.5);
    for(int i = 0; i < DebugIdPoints.size(); i++){
      coordinate_type c1 = DebugIdPoints[i];
      char str[10];
      sprintf(str,"%d",DebugIds[i]);
      glPointSize(5.0);
      glBegin(GL_POINTS);
      glVertex3d(c1[0],c1[1],c1[2]);
      glEnd();
      renderString(c1[0],c1[1],c1[2],0.0001,str);
    }
    glPointSize(1.0);
  }

  void Debugger::draw_points(){
    glPointSize(10.0);
    glBegin(GL_POINTS);
    glColor3f(1.0,0.0,0.5);
    for(int i = 0; i < DebugPoints.size(); i++){
      coordinate_type c1 = DebugPoints[i];
      glVertex3d(c1[0],c1[1],c1[2]);
    }
    glEnd();
    glPointSize(1.0);
  }

  void Debugger::draw_points_0(){
    glPointSize(10.0);
    glBegin(GL_POINTS);
    glColor3f(0.1,1.0,0.0);
    for(int i = 0; i < DebugPoints0.size(); i++){
      coordinate_type c1 = DebugPoints0[i];
      glVertex3d(c1[0],c1[1],c1[2]);
    }
    glEnd();
    glPointSize(1.0);
  }

  void Debugger::draw_mag_points(){
    glPointSize(2.0);

    for(int i = 0; i < DebugMagPoints.size(); i++){
      double mag = DebugMags[i];
      coordinate_type c1 = DebugMagPoints[i];
      double amag = abs(mag)*10.;
      if(mag < 0.0)
	glColor3f(amag,0,0);
      else
	glColor3f(0,amag,0);

      glBegin(GL_POINTS);
      glVertex3d(c1[0],c1[1],c1[2]);
      glEnd();
    }
  }
    
  void Debugger::draw_lines(){
    glLineWidth(2.0);
    glBegin(GL_LINES);
    glColor3f(1.0,0.0,0.0);
    for(int i = 0; i < DebugLines.size(); i+=2){
      coordinate_type c0 = DebugLines[i];
      coordinate_type c1 = DebugLines[i+1];
      glVertex3d(c0[0],c0[1],c0[2]);
      glVertex3d(c1[0],c1[1],c1[2]);
    }
    glEnd();
  }

  void Debugger::draw_lines_0(){
    glLineWidth(2.0);
    glBegin(GL_LINES);
    glColor3f(0.75,0.75,0.0);
    for(int i = 0; i < DebugLines0.size(); i+=2){
      coordinate_type c0 = DebugLines0[i];
      coordinate_type c1 = DebugLines0[i+1];
      glVertex3d(c0[0],c0[1],c0[2]);
      glVertex3d(c1[0],c1[1],c1[2]);
    }
    glEnd();
  }

  void Debugger::draw_lines_1(){
    glLineWidth(2.0);
    glBegin(GL_LINES);
    glColor3f(0.75,0.0,0.75);
    for(int i = 0; i < DebugLines1.size(); i+=2){
      coordinate_type c0 = DebugLines1[i];
      coordinate_type c1 = DebugLines1[i+1];
      glVertex3d(c0[0],c0[1],c0[2]);
      glVertex3d(c1[0],c1[1],c1[2]);
    }
    glEnd();
  }

  void Debugger::draw_triangles(){      
    for(int i = 0; i < DebugTriangles.size(); i+=3){
      glLineWidth(2.0);
      glBegin(GL_POLYGON);
      coordinate_type c0 = DebugTriangles[i+0];
      coordinate_type c1 = DebugTriangles[i+1];
      coordinate_type c2 = DebugTriangles[i+2];
      glColor3f(0.0,1.0,0.0);
      glVertex3d(c0[0],c0[1],c0[2]);
      glVertex3d(c1[0],c1[1],c1[2]);
      glVertex3d(c2[0],c2[1],c2[2]);
      glEnd();
    }
    glEnd();
  }

  void Debugger::draw_triangles_0(){      
    for(int i = 0; i < DebugTriangles0.size(); i+=3){
      glLineWidth(5.0);
      glBegin(GL_POLYGON);
      coordinate_type c0 = DebugTriangles0[i+0];
      coordinate_type c1 = DebugTriangles0[i+1];
      coordinate_type c2 = DebugTriangles0[i+2];
      glColor3f(0.75,0.75,0.75);
      glVertex3d(c0[0],c0[1],c0[2]);
      glVertex3d(c1[0],c1[1],c1[2]);
      glVertex3d(c2[0],c2[1],c2[2]);
      glEnd();
    }
    glEnd();
  }

  void Debugger::draw_cached_triangles(){      
    for(int i = 0; i < CachedTriangles.size(); i+=3){
      glLineWidth(1.0);
      glBegin(GL_POLYGON);
      coordinate_type c0 = CachedTriangles[i+0];
      coordinate_type c1 = CachedTriangles[i+1];
      coordinate_type c2 = CachedTriangles[i+2];
      glColor3f(0.5,0.5,0.5);
      glVertex3d(c0[0],c0[1],c0[2]);
      glVertex3d(c1[0],c1[1],c1[2]);
      glVertex3d(c2[0],c2[1],c2[2]);
      glEnd();
    }
    glEnd();
  }

  void Debugger::draw_cached_triangle_lines(){      
    for(int i = 0; i < CachedTriangles.size(); i+=3){
      glLineWidth(1.0);
      glBegin(GL_LINES);
      coordinate_type c0 = CachedTriangles[i+0];
      coordinate_type c1 = CachedTriangles[i+1];
      coordinate_type c2 = CachedTriangles[i+2];
      glColor3f(0.5,0.5,0.5);
      glVertex3d(c0[0],c0[1],c0[2]);
      glVertex3d(c1[0],c1[1],c1[2]);
      glVertex3d(c1[0],c1[1],c1[2]);
      glVertex3d(c2[0],c2[1],c2[2]);
      glVertex3d(c2[0],c2[1],c2[2]);
      glVertex3d(c0[0],c0[1],c0[2]);
      glEnd();
    }
    glEnd();
  }
	

  void Debugger::draw_boxes(){
    glLineWidth(1.0);
    glColor3f(0.5,0.5,0.8);
    for(int i = 0; i < DebugBoxes.size(); i+=2){
      coordinate_type p = DebugBoxes[i];
      coordinate_type h = DebugBoxes[i+1];
      glBegin(GL_LINES);
      glVertex3d(p[0] - h[0],p[1] - h[1],p[2]-h[2]);
      glVertex3d(p[0] + h[0],p[1] - h[1],p[2]-h[2]);

      glVertex3d(p[0] + h[0],p[1] - h[1],p[2]-h[2]);
      glVertex3d(p[0] + h[0],p[1] + h[1],p[2]-h[2]);

      glVertex3d(p[0] + h[0],p[1] + h[1],p[2]-h[2]);
      glVertex3d(p[0] - h[0],p[1] + h[1],p[2]-h[2]);

      glVertex3d(p[0] - h[0],p[1] + h[1],p[2]-h[2]);
      glVertex3d(p[0] - h[0],p[1] - h[1],p[2]-h[2]);


      glVertex3d(p[0] - h[0],p[1] - h[1],p[2]+h[2]);
      glVertex3d(p[0] + h[0],p[1] - h[1],p[2]+h[2]);

      glVertex3d(p[0] + h[0],p[1] - h[1],p[2]+h[2]);
      glVertex3d(p[0] + h[0],p[1] + h[1],p[2]+h[2]);

      glVertex3d(p[0] + h[0],p[1] + h[1],p[2]+h[2]);
      glVertex3d(p[0] - h[0],p[1] + h[1],p[2]+h[2]);

      glVertex3d(p[0] - h[0],p[1] + h[1],p[2]+h[2]);
      glVertex3d(p[0] - h[0],p[1] - h[1],p[2]+h[2]);


      glVertex3d(p[0] - h[0],p[1] - h[1],p[2] - h[2]);
      glVertex3d(p[0] - h[0],p[1] - h[1],p[2] + h[2]);

      glVertex3d(p[0] + h[0],p[1] - h[1],p[2] - h[2]);
      glVertex3d(p[0] + h[0],p[1] - h[1],p[2] + h[2]);

      glVertex3d(p[0] + h[0],p[1] + h[1],p[2] - h[2]);
      glVertex3d(p[0] + h[0],p[1] + h[1],p[2] + h[2]);

      glVertex3d(p[0] - h[0],p[1] + h[1],p[2] - h[2]);
      glVertex3d(p[0] - h[0],p[1] + h[1],p[2] + h[2]);

      glVertex3d(p[0] - h[0],p[1] - h[1],p[2] - h[2]);
      glVertex3d(p[0] - h[0],p[1] - h[1],p[2] + h[2]);
      glEnd();
    }

  }
}//end m2
