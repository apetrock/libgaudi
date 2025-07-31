#ifdef __EMSCRIPTEN__
#include <emscripten/emscripten.h>
#else
#define EMSCRIPTEN_KEEPALIVE
#endif

#ifdef __EMSCRIPTEN__
#include <emscripten/emscripten.h>
#else
#define EMSCRIPTEN_KEEPALIVE
#endif

#include <emscripten/val.h>
#include <emscripten/bind.h>

#include "../include/wasm_geometry_logger.h"
#include "gaudi/logger.hpp"

using namespace gaudi;

emscripten::val vector_string_to_js_array(const std::vector<std::string> &vec)
{
  emscripten::val jsArray = emscripten::val::array();
  for (const auto &str : vec)
  {
    jsArray.call<void>("push", emscripten::val(str));
  }
  return jsArray;
}

emscripten::val vector_float_to_js_array(const std::vector<float> &data)
{
  emscripten::val view{emscripten::typed_memory_view(data.size(), data.data())};
  auto result = emscripten::val::global("Float32Array").new_(data.size());
  result.call<void>("set", view);
  return result;
}

emscripten::val vector_int_to_js_array(const std::vector<int> &data)
{
  emscripten::val view{emscripten::typed_memory_view(data.size(), data.data())};
  auto result = emscripten::val::global("Int32Array").new_(data.size());
  result.call<void>("set", view);
  return result;
}

// WASM Logger API functions
namespace gaudi {

// Get the lines data as Float32Array
emscripten::val get_lines_data() {
  wasm_geometry_logger &logger = wasm_geometry_logger::get_instance();
  return vector_float_to_js_array(logger._lines);
}

// Get the line colors data as Float32Array
emscripten::val get_line_colors_data() {
  wasm_geometry_logger &logger = wasm_geometry_logger::get_instance();
  return vector_float_to_js_array(logger._line_colors);
}

// Get the points data as Float32Array
emscripten::val get_points_data() {
  wasm_geometry_logger &logger = wasm_geometry_logger::get_instance();
  return vector_float_to_js_array(logger._points);
}

// Get the point colors data as Float32Array
emscripten::val get_point_colors_data() {
  wasm_geometry_logger &logger = wasm_geometry_logger::get_instance();
  return vector_float_to_js_array(logger._point_colors);
}

// Get line count (number of line segments)
int get_line_count() {
  wasm_geometry_logger &logger = wasm_geometry_logger::get_instance();
  return logger._lines.size() / 6; // 6 floats per line (2 points * 3 coordinates)
}

// Get point count
int get_point_count() {
  wasm_geometry_logger &logger = wasm_geometry_logger::get_instance();
  return logger._points.size() / 3; // 3 floats per point
}

// Get line color count (number of color entries for lines)
int get_line_color_count() {
  wasm_geometry_logger &logger = wasm_geometry_logger::get_instance();
  return logger._line_colors.size() / 4; // 4 floats per color
}

// Get point color count
int get_point_color_count() {
  wasm_geometry_logger &logger = wasm_geometry_logger::get_instance();
  return logger._point_colors.size() / 4; // 4 floats per color
}

// Wrapper functions for the static methods to work with emscripten
void wasm_geometry_logger_point(float x, float y, float z, float r, float g, float b, float a) {
  wasm_geometry_logger::point(vec3(x, y, z), vec4(r, g, b, a));
}

void wasm_geometry_logger_line(float x0, float y0, float z0, float x1, float y1, float z1, 
                      float r, float g, float b, float a) {
  wasm_geometry_logger::line(vec3(x0, y0, z0), vec3(x1, y1, z1), vec4(r, g, b, a));
}

void wasm_geometry_logger_box(float cx, float cy, float cz, float hx, float hy, float hz,
                     float r, float g, float b, float a) {
  wasm_geometry_logger::box(vec3(cx, cy, cz), vec3(hx, hy, hz), vec4(r, g, b, a));
}

void wasm_geometry_logger_ext(float min_x, float min_y, float min_z, 
                     float max_x, float max_y, float max_z,
                     float r, float g, float b, float a) {
  wasm_geometry_logger::ext(vec3(min_x, min_y, min_z), vec3(max_x, max_y, max_z), vec4(r, g, b, a));
}

void wasm_geometry_logger_frame_wrapper(const emscripten::val& M, float cx, float cy, float cz, float C) {
  // Convert JS array to mat3 - assuming M is a 3x3 matrix as array
  mat3 matrix;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      matrix(i, j) = M[i * 3 + j].as<double>();
    }
  }
  wasm_geometry_logger::frame(matrix, vec3(cx, cy, cz), C);
}

void wasm_geometry_logger_clear() {
  wasm_geometry_logger::clear();
}

void wasm_geometry_logger_render() {
  wasm_geometry_logger::render();
}

} // namespace gaudi

// Emscripten bindings
EMSCRIPTEN_BINDINGS(wasm_geometry_logger_module) {
  using namespace emscripten;
  
  // Register vector types
  register_vector<float>("vector_float");
  register_vector<int>("vector_int");
  register_vector<std::string>("vector_string");
  
  // Bind the convenience functions
  function("vector_string_to_js_array", &vector_string_to_js_array);
  function("vector_float_to_js_array", &vector_float_to_js_array);
  function("vector_int_to_js_array", &vector_int_to_js_array);
  
  // Bind the wasm_geometry_logger data access functions
  function("get_lines_data", &gaudi::get_lines_data);
  function("get_line_colors_data", &gaudi::get_line_colors_data);
  function("get_points_data", &gaudi::get_points_data);
  function("get_point_colors_data", &gaudi::get_point_colors_data);
  function("get_line_count", &gaudi::get_line_count);
  function("get_point_count", &gaudi::get_point_count);
  function("get_line_color_count", &gaudi::get_line_color_count);
  function("get_point_color_count", &gaudi::get_point_color_count);
  
  // Bind the wasm_geometry_logger drawing functions
  function("wasm_geometry_logger_point", &gaudi::wasm_geometry_logger_point);
  function("wasm_geometry_logger_line", &gaudi::wasm_geometry_logger_line);
  function("wasm_geometry_logger_box", &gaudi::wasm_geometry_logger_box);
  function("wasm_geometry_logger_ext", &gaudi::wasm_geometry_logger_ext);
  function("wasm_geometry_logger_frame", &gaudi::wasm_geometry_logger_frame_wrapper);
  function("wasm_geometry_logger_clear", &gaudi::wasm_geometry_logger_clear);
  function("wasm_geometry_logger_render", &gaudi::wasm_geometry_logger_render);
}