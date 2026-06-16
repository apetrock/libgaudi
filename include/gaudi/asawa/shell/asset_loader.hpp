
#ifndef __ASAWA_PRIM_LOAD__
#define __ASAWA_PRIM_LOAD__

#include "../datums.hpp"
#include "../faceloader.hpp"
#include "../objloader_refactor.hpp"
#include "../primitive_objects.hpp"
#include "gaudi/paths.hpp"
#include "shell.hpp"

namespace gaudi {
namespace asawa {
namespace shell {

shell::ptr load_obj(const std::string &name) {
  const std::filesystem::path resolved = gaudi::resolve_path(name);
  std::string file = resolved.string();
  // std::string file("assets/skeleton.obj");
  std::vector<vec3> vertices;
  std::vector<std::vector<int>> faces;
  loadObjfile(file, vertices, faces);
  if (vertices.empty() || faces.empty()) {
    throw std::runtime_error("OBJ asset '" + file +
                             "' loaded no geometry");
  }
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
shell::ptr load_messer() { return load_obj("assets/messer.obj"); }
shell::ptr load_crab() { return load_obj("assets/blue-crab-low.obj"); }

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

shell::ptr load_tet() {
  std::vector<vec3> vertices;
  std::vector<std::vector<int>> faces;
  make_tet(vertices, faces);
  std::vector<index_t> corners_next, corners_vert, corners_face;
  assemble_table(vertices, faces, corners_next, corners_vert, corners_face);
  shell::ptr M = shell::create(corners_next, corners_vert, corners_face);
  datum_t<vec3>::ptr vdata = datum_t<vec3>::create(prim_type::VERTEX, vertices);
  M->insert_datum(vdata);

  return M;
}

shell::ptr load_sphere(real radius = 1.0, int u_segments = 24,
                       int v_segments = 16) {
  std::vector<vec3> vertices;
  std::vector<std::vector<int>> faces;
  make_sphere(vertices, faces, radius, u_segments, v_segments);
  std::vector<index_t> corners_next, corners_vert, corners_face;
  assemble_table(vertices, faces, corners_next, corners_vert, corners_face);
  shell::ptr M = shell::create(corners_next, corners_vert, corners_face);
  datum_t<vec3>::ptr vdata = datum_t<vec3>::create(prim_type::VERTEX, vertices);
  M->insert_datum(vdata);

  return M;
}

shell::ptr load_torus(real major_radius = 1.0, real minor_radius = 0.33,
                      int major_segments = 48, int minor_segments = 24) {
  std::vector<vec3> vertices;
  std::vector<std::vector<int>> faces;
  make_torus(vertices, faces, major_radius, minor_radius, major_segments,
             minor_segments);
  std::vector<index_t> corners_next, corners_vert, corners_face;
  assemble_table(vertices, faces, corners_next, corners_vert, corners_face);
  shell::ptr M = shell::create(corners_next, corners_vert, corners_face);
  datum_t<vec3>::ptr vdata = datum_t<vec3>::create(prim_type::VERTEX, vertices);
  M->insert_datum(vdata);

  return M;
}
} // namespace shell
} // namespace asawa
} // namespace gaudi
#endif