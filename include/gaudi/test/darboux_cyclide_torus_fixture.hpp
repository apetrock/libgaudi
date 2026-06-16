#ifndef __GAUDI_DARBOUX_CYCLIDE_TORUS_FIXTURE_HPP__
#define __GAUDI_DARBOUX_CYCLIDE_TORUS_FIXTURE_HPP__

#include "gaudi/asawa/datums.hpp"
#include "gaudi/asawa/faceloader.hpp"
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/common.h"

#include <cmath>
#include <vector>

namespace gaudi {
namespace test {

struct TorusFrame {
  vec3 center = vec3::Zero();
  vec3 x_axis = vec3::UnitX();
  vec3 y_axis = vec3::UnitY();
  vec3 z_axis = vec3::UnitZ();
};

inline TorusFrame make_torus_frame(const vec3 &center, const vec3 &axis) {
  TorusFrame frame;
  frame.center = center;
  frame.z_axis = axis.normalized();
  vec3 seed = std::abs(frame.z_axis.dot(vec3::UnitZ())) > 0.9 ? vec3::UnitX()
                                                              : vec3::UnitZ();
  frame.x_axis = seed.cross(frame.z_axis).normalized();
  frame.y_axis = frame.z_axis.cross(frame.x_axis).normalized();
  return frame;
}

struct TorusMesh {
  asawa::shell::shell::ptr shell;
  std::vector<vec3> vertices;
};

inline vec3 torus_point(real u, real v, real major_radius, real minor_radius,
                        const TorusFrame &frame) {
  const real cu = std::cos(u);
  const real su = std::sin(u);
  const real cv = std::cos(v);
  const real sv = std::sin(v);
  const real rho = major_radius + minor_radius * cv;
  return frame.center + rho * cu * frame.x_axis + rho * su * frame.y_axis +
         minor_radius * sv * frame.z_axis;
}

inline TorusMesh make_offset_torus_shell(int major_segments, int minor_segments,
                                         real major_radius, real minor_radius,
                                         const TorusFrame &frame) {
  TorusMesh out;
  out.vertices.reserve(static_cast<size_t>(major_segments * minor_segments));
  for (int i = 0; i < major_segments; ++i) {
    const real u = 2.0 * M_PI * real(i) / real(major_segments);
    for (int j = 0; j < minor_segments; ++j) {
      const real v = 2.0 * M_PI * real(j) / real(minor_segments);
      out.vertices.push_back(
          torus_point(u, v, major_radius, minor_radius, frame));
    }
  }

  std::vector<std::vector<int>> faces;
  faces.reserve(static_cast<size_t>(2 * major_segments * minor_segments));
  auto vid = [major_segments, minor_segments](int i, int j) {
    return ((i % major_segments + major_segments) % major_segments) *
               minor_segments +
           ((j % minor_segments + minor_segments) % minor_segments);
  };

  for (int i = 0; i < major_segments; ++i) {
    for (int j = 0; j < minor_segments; ++j) {
      const int i1 = (i + 1) % major_segments;
      const int j1 = (j + 1) % minor_segments;
      faces.push_back({vid(i, j), vid(i1, j), vid(i1, j1)});
      faces.push_back({vid(i, j), vid(i1, j1), vid(i, j1)});
    }
  }

  std::vector<index_t> corners_next;
  std::vector<index_t> corners_vert;
  std::vector<index_t> corners_face;
  asawa::assemble_table(out.vertices, faces, corners_next, corners_vert,
                        corners_face);
  out.shell =
      asawa::shell::shell::create(corners_next, corners_vert, corners_face);
  out.shell->insert_datum(
      asawa::datum_t<vec3>::create(asawa::prim_type::VERTEX, out.vertices));
  return out;
}

} // namespace test
} // namespace gaudi

#endif // __GAUDI_DARBOUX_CYCLIDE_TORUS_FIXTURE_HPP__
