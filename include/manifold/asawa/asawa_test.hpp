#include "m2.hpp"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "manifold/vec_addendum.h"

#include "GaudiGraphics/geometry_logger.h"

#include "dynamic_surface.hpp"
#include "m2_refactor.hpp"
#include "primitive_objects.hpp"
#include "primitive_operations.hpp"
#include <vector>
#include <zlib.h>

#ifndef __M2REFACTOR_TEST__
#define __M2REFACTOR_TEST__
namespace asawa {
typedef double real;
typedef int index_t;
typedef Eigen::Matrix<real, 3, 1> vec3;
typedef Eigen::Matrix<real, 4, 1> vec4;

void debug_manifold(manifold &M, const std::vector<vec3> verts) {
  for (int i = 0; i < M.__corners_next.size(); i += 2) {
    if (M.__corners_next[i] < 0)
      continue;
    int i0 = i;
    int i1 = M.other(i0);
    int v0 = M.vert(i0);
    int v1 = M.vert(i1);
    gg::geometry_logger::line(verts[v0], verts[v1], vec4(0.5, 0.5, 0.5, 1.0));
  }
}

real avg_length(manifold &M, const std::vector<vec3> &coords) {
  real accum = 0.0;
  for (int i = 0; i < M.__corners_next.size(); i += 2) {
    if (M.__corners_next[i] < 0)
      continue;
    int i0 = i;
    int i1 = M.other(i0);
    int v0 = M.vert(i0);
    int v1 = M.vert(i1);
    accum += (coords[M.vert(i0)] - coords[M.vert(i1)]).norm();
  }
  return 0.5 * accum / real(M.corner_count());
}

manifold::ptr build_cube() {
  std::vector<vec3> vertices;
  std::vector<std::vector<int>> faces;
  make_cube(vertices, faces);
  std::vector<index_t> corners_next, corners_prev, corners_vert, corners_face;
  assemble_table(vertices, faces, corners_next, corners_prev, corners_vert,
                 corners_face);
  manifold::ptr M = manifold::create(corners_next, corners_vert, corners_face);
  datum_t<vec3>::ptr vdata = datum_t<vec3>::create(prim_type::VERTEX, vertices);
  M->insert_datum(vdata);

  return M;
}

void test() {
  manifold::ptr M = build_cube();
  triangulate(*M);

  vec3_datum::ptr v_datum = static_pointer_cast<vec3_datum>(M->get_datum(0));
  const std::vector<vec3> &coords = v_datum->data();
  real l0 = avg_length(*M, v_datum->data());

  dynamic_surface::ptr surf =
      dynamic_surface::create(M, 0.5 * l0, 1.5 * l0, 0.5 * l0);

  remove_vertex(*M, 0);
  remove_vertex(*M, 1);
  remove_vertex(*M, 2);
  remove_vertex(*M, 3);
  remove_vertex(*M, 4);

  //  remove_vertex(*M, 2);
  //   subdivide_edges(*M);
  //   collapse_edges(*M);

  //  gather_edges_parallel(*M, vdata->data(), 1.0);
  //   merge_face(*M, 3, M->other(3));
  debug_manifold(*M, v_datum->data());
}

class asawa_dynamic_test {
public:
  typedef std::shared_ptr<asawa_dynamic_test> ptr;

  static ptr create() { return std::make_shared<asawa_dynamic_test>(); }

  asawa_dynamic_test() {
    __M = build_cube();
    triangulate(*__M);

    vec3_datum::ptr c_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(0));

    const std::vector<vec3> &coords = c_datum->data();
    real l0 = 0.25 * avg_length(*__M, c_datum->data());

    __surf = dynamic_surface::create(__M, 0.5 * l0, 3.0 * l0, 0.5 * l0);
  };

  void step(int frame) {

    vec3_datum::ptr v_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(1));
    std::vector<vec3> &dx = v_datum->data();
    std::cout << "frame: " << frame << std::endl;

    if (frame == 1) {

      dx[0] = vec3(0.0, 0.01, -0.005);
      dx[2] = vec3(0.01, 0.00, -0.005);
      dx[4] = vec3(-0.01, 0.0, -0.005);
      dx[6] = vec3(0.00, -0.01, -0.005);

      dx[1] = vec3(0.01, 0.0, 0.005);
      dx[3] = vec3(0.00, -0.01, 0.005);
      dx[5] = vec3(0.00, 0.01, 0.005);
      dx[7] = vec3(-0.01, 0.00, 0.005);
    }
    __surf->step(dx);

    vec3_datum::ptr c_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(0));
    std::vector<vec3> &x = c_datum->data();

    debug_manifold(*__M, c_datum->data());
  }

  manifold::ptr __M;
  dynamic_surface::ptr __surf;
};

} // namespace asawa
#endif