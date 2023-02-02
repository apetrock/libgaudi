#ifndef __MESH_HELPER__
#define __MESH_HELPER__

#include <GaudiGraphics/buffers.hpp>
#include <GaudiGraphics/viewer.hpp>
#include <GaudiMath/typedefs.hpp>
#include <cassert>

#include <gaudi/asawa/asawa.h>

#include <gaudi/asawa/manifold.hpp>

#include <gaudi/vec_addendum.h>

#include <nanogui/glutil.h>

#include <iostream>
#include <string>
#include <vector>

namespace gg {
using namespace gaudi;
//////////////////////////////////////////////////
// refactor
//////////////////////////////////////////////////
struct colorRGB {
  double r, g, b, a;
  colorRGB() {
    r = 0.5;
    g = 0.5;
    b = 0.5;
    a = 1.0;
  }
  colorRGB(double rr, double gg, double bb, double aa) {
    r = rr;
    g = gg;
    b = bb;
    a = aa;
  }
};
void fillBuffer_ref(asawa::manifold &M, gg::BufferObjectPtr obj,
                    colorRGB col = colorRGB(1.0, 0.5, 0.5, 1.0)) {
  using namespace nanogui;

  asawa::vec3_datum::ptr x_datum =
      static_pointer_cast<asawa::vec3_datum>(M.get_datum(0));
  const std::vector<gaudi::vec3> &x = x_datum->data();

  int numVerts = 0;
  int numIndices = 0;

  std::vector<std::vector<int>> faces;
  for (int i = 0; i < M.face_count(); i++) {
    std::vector<int> face;
    if (M.fbegin(i) < 0)
      continue;
    if (M.fsize(i) != 3)
      continue;
    M.for_each_face(
        i, [&face](int ci, asawa::manifold &M) { face.push_back(M.vert(ci)); });

    numIndices += face.size();
    faces.push_back(face);
  }
  numVerts = x.size();

  obj->fillBuffer([&](gg::BufferObject &o) -> void {
    o.allocateVerts(faces.size(), numVerts);
    auto &indices = o.indices();
    auto &positions = o.positions();
    auto &colors = o.colors();
    for (int i = 0; i < faces.size(); i++) {
      for (int j = 0; j < faces[i].size(); j++) {
        indices.col(i)[j] = faces[i][j];
      }
    }

    for (int i = 0; i < x.size(); i++) {
      for (int j = 0; j < 3; j++) {
        positions.col(i)[j] = x[i][j];
      }
      colors.col(i)[0] = col.r;
      colors.col(i)[1] = col.g;
      colors.col(i)[2] = col.b;
    }
    // std::cout << indices << std::endl;
  });
};
} // namespace gg

#endif