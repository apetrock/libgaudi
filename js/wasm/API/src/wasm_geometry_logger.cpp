/*
 *  wasm_geometry_logger.cpp
 *  WASM Geometry Logger Implementation
 *
 *  Created by John Delaney on 12/29/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

// #include <gl/glut.h>
#include "../include/wasm_geometry_logger.h"
#include <numeric>

namespace gaudi {

//////////////////////////
// debugger
//////////////////////////

bool wasm_geometry_logger::instance_flag = false;
wasm_geometry_logger *wasm_geometry_logger::global_instance = NULL;

// Convenience functions for appending vectors to float arrays
inline void appendVec3(const vec3 &v, std::vector<float> &floatArray) {
  floatArray.push_back(v[0]);
  floatArray.push_back(v[1]);
  floatArray.push_back(v[2]);
}

inline void appendVec4(const vec4 &v, std::vector<float> &floatArray) {
  floatArray.push_back(v[0]);
  floatArray.push_back(v[1]);
  floatArray.push_back(v[2]);
  floatArray.push_back(v[3]);
}

wasm_geometry_logger &wasm_geometry_logger::get_instance() {
  static wasm_geometry_logger logger;
  if (!logger.initialized()) {
    logger.initialized() = true;
  }

  return logger;
}

void wasm_geometry_logger::render() {
  wasm_geometry_logger &logger = wasm_geometry_logger::get_instance();
  // Implementation would depend on how you want to render the float data
  // This is a placeholder for the rendering logic
}

void wasm_geometry_logger::clear() {
  wasm_geometry_logger &logger = wasm_geometry_logger::get_instance();
  logger._lines.clear();
  logger._line_colors.clear();
  logger._points.clear();
  logger._point_colors.clear();
}

void wasm_geometry_logger::point(const vec3 &p0, const vec4 &color) {
  wasm_geometry_logger &logger = wasm_geometry_logger::get_instance();
  
  // Append point coordinates (3 floats)
  appendVec3(p0, logger._points);
  
  // Append point color (4 floats)
  appendVec4(color, logger._point_colors);
}

void wasm_geometry_logger::line(const vec3 &p0, const vec3 &p1, const vec4 &color) {
  wasm_geometry_logger &logger = wasm_geometry_logger::get_instance();
  
  // Append line start coordinates (3 floats)
  appendVec3(p0, logger._lines);
  
  // Append line end coordinates (3 floats)
  appendVec3(p1, logger._lines);
  
  // Append line color (4 floats) - same color for both endpoints
  appendVec4(color, logger._line_colors);
  appendVec4(color, logger._line_colors);
}

void wasm_geometry_logger::box(const vec3 &cen, const vec3 &h, const vec4 &color) {
  wasm_geometry_logger &logger = wasm_geometry_logger::get_instance();
  
  // Calculate box corners
  vec3 min = cen - h;
  vec3 max = cen + h;

  // Draw the 12 edges of the box
  // Bottom face
  logger.line(vec3(min[0], min[1], min[2]), 
              vec3(max[0], min[1], min[2]), color);
  logger.line(vec3(max[0], min[1], min[2]), 
              vec3(max[0], min[1], max[2]), color);
  logger.line(vec3(max[0], min[1], max[2]), 
              vec3(min[0], min[1], max[2]), color);
  logger.line(vec3(min[0], min[1], max[2]), 
              vec3(min[0], min[1], min[2]), color);
  
  // Top face
  logger.line(vec3(min[0], max[1], min[2]), 
              vec3(max[0], max[1], min[2]), color);
  logger.line(vec3(max[0], max[1], min[2]), 
              vec3(max[0], max[1], max[2]), color);
  logger.line(vec3(max[0], max[1], max[2]), 
              vec3(min[0], max[1], max[2]), color);
  logger.line(vec3(min[0], max[1], max[2]), 
              vec3(min[0], max[1], min[2]), color);
  
  // Vertical edges
  logger.line(vec3(min[0], min[1], min[2]), 
              vec3(min[0], max[1], min[2]), color);
  logger.line(vec3(max[0], min[1], min[2]), 
              vec3(max[0], max[1], min[2]), color);
  logger.line(vec3(max[0], min[1], max[2]), 
              vec3(max[0], max[1], max[2]), color);
  logger.line(vec3(min[0], min[1], max[2]), 
              vec3(min[0], max[1], max[2]), color);
}

void wasm_geometry_logger::ext(const vec3 &mn, const vec3 &mx, const vec4 &color) {
  wasm_geometry_logger &logger = wasm_geometry_logger::get_instance();

  vec3 cen = 0.5 * (mx + mn);
  vec3 h = 0.5 * (mx - mn);

  logger.box(cen, h, color);
}

vec3 _rainbow(double d) {
  double r = 0.5 + 0.5 * cos(2.0 * M_PI * (d + 0.000));
  double g = 0.5 + 0.5 * cos(2.0 * M_PI * (d + 0.333));
  double b = 0.5 + 0.5 * cos(2.0 * M_PI * (d + 0.666));
  return vec3(r, g, b);
}

void wasm_geometry_logger::frame(mat3 M, vec3 c, double C) {
  vec3 t0 = M.col(0);
  vec3 t1 = M.col(1);
  vec3 t2 = M.col(2);

  wasm_geometry_logger::line(c - 0.5 * C * t0, c + C * t0, vec4(1.0, 0.0, 0.0, 1.0));
  wasm_geometry_logger::line(c - 0.5 * C * t1, c + C * t1, vec4(0.0, 1.0, 0.0, 1.0));
  wasm_geometry_logger::line(c - 0.5 * C * t2, c + C * t2, vec4(0.0, 0.0, 1.0, 1.0));
}

} // namespace gg
