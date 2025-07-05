/*
 *  manifold_singleton.h
 *  Phase Vocoder
 *
 *  Created by John Delaney on 12/29/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include "gaudi/common.h"
#include <emscripten/val.h>

#ifndef __WASM_DEBUG_INTERFACE__
#define __WASM_DEBUG_INTERFACE__

namespace gaudi {

enum PresetColor {
  grey,
  red,
  green,
  blue,
  rainbow,
};

inline vec4 rainbow4(double d) {
  double r = 0.5 + 0.5 * cos(2.0 * M_PI * (d + 0.000));
  double g = 0.5 + 0.5 * cos(2.0 * M_PI * (d + 0.333));
  double b = 0.5 + 0.5 * cos(2.0 * M_PI * (d + 0.666));
  return vec4(r, g, b, 1.0);
}

inline vec4 sdf4(double d) {
  vec4 inside(0.0, 1.0, 0.0, 1.0);
  vec4 outside(1.0, 0.0, 0.0, 1.0);
  if (d < 0)
    return abs(d) * inside;
  else
    return abs(d) * outside;
}

class wasm_logger {

public:

  static wasm_logger &get_instance();
  static void point(const vec3 &p0, const vec4 &color);
  static void line(const vec3 &p0, const vec3 &p1, const vec4 &color);

  static void box(const vec3 &cen, const vec3 &h, const vec4 &col);
  static void ext(const vec3 &mn, const vec3 &mx, const vec4 &col);

  static void frame(mat3 M, vec3 c, double C);
  
  // Add missing declarations
  static void render();
  static void clear();

  bool &initialized() { return instance_flag; }
  bool initialized() const { return instance_flag; }
  std::vector<float> _lines;
  std::vector<float> _line_colors;
  std::vector<float> _points;
  std::vector<float> _point_colors;

private:
  wasm_logger() {
  }

  wasm_logger(const wasm_logger &);
  wasm_logger &operator=(const wasm_logger &);

  static wasm_logger *global_instance;
  static bool instance_flag;
};

// Global logger API functions for WASM modules
emscripten::val get_lines_data();
emscripten::val get_line_colors_data();
emscripten::val get_points_data();
emscripten::val get_point_colors_data();
int get_line_count();
int get_point_count();
int get_line_color_count();
int get_point_color_count();

} // namespace gaudi

#endif
