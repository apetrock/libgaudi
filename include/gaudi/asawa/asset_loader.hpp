
#ifndef __ASAWA_PRIM_LOAD__
#define __ASAWA_PRIM_LOAD__

#include "faceloader.hpp"
#include "objloader_refactor.hpp"
#include "primitive_objects.hpp"
#include "shell.hpp"

namespace gaudi {
namespace asawa {

shell::ptr load_obj(const std::string &name) {
  std::string file(name);
  // std::string file("assets/skeleton.obj");
  std::vector<vec3> vertices;
  std::vector<std::vector<int>> faces;
  loadObjfile(file, vertices, faces);
  // make_cube(vertices, faces);
  // faces.pop_back();
  std::vector<index_t> corners_next, corners_vert, corners_face;
  assemble_table(vertices, faces, corners_next, corners_vert, corners_face);
  shell::ptr M = shell::create(corners_next, corners_vert, corners_face);
  datum_t<vec3>::ptr vdata = datum_t<vec3>::create(prim_type::VERTEX, vertices);
  M->insert_datum(vdata);

  return M;
}
shell::ptr load_bunny() { return load_obj("assets/bunny.obj"); }
shell::ptr load_skeleton() { return load_obj("assets/skeleton.obj"); }

shell::ptr load_heart() { return load_obj("assets/heart.obj"); }

shell::ptr load_cube() {
  std::vector<vec3> vertices;
  std::vector<std::vector<int>> faces;
  make_cube(vertices, faces);
  std::vector<index_t> corners_next, corners_vert, corners_face;
  assemble_table(vertices, faces, corners_next, corners_vert, corners_face);
  shell::ptr M = shell::create(corners_next, corners_vert, corners_face);
  datum_t<vec3>::ptr vdata = datum_t<vec3>::create(prim_type::VERTEX, vertices);
  M->insert_datum(vdata);

  return M;
}
} // namespace asawa
} // namespace gaudi
#endif