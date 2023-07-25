#ifndef __MESH_HELPER__
#define __MESH_HELPER__

#include "gaudi/common.h"
#include <GaudiGraphics/buffers.hpp>
#include <GaudiGraphics/viewer.hpp>
#include <GaudiMath/typedefs.hpp>
#include <cassert>

#include <cmath>
#include <gaudi/asawa/asawa.h>

#include <gaudi/asawa/shell/shell.hpp>

#include <gaudi/asawa/rod/rod.hpp>

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
void fillBuffer_ref(asawa::shell::shell &M, gg::BufferObjectPtr obj,
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
    // std::cout << i << " " << M.fsize(i) << " ";
    M.for_each_face(i, [&face](int ci, asawa::shell::shell &M) {
      // std::cout << ci << " " << M.vert(ci) << " ";
      face.push_back(M.vert(ci));
    });
    // std::cout << std::endl;

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

void fillBuffer_ref(asawa::rod::rod &R, gg::BufferObjectPtr obj,
                    colorRGB col = colorRGB(1.0, 0.5, 0.5, 0.5)) {
  using namespace nanogui;

  const std::vector<gaudi::vec3> &x = R.__x;

  std::vector<gaudi::vec3> circ;
  int Nc = 32;
  matX x_sect(3, Nc);
  double r = 0.5 * R._r;
  for (int i = 0; i < Nc; i++) {
    double thet = double(i) / double(Nc) * 2.0 * M_PI;
    double x = r * cos(thet);
    double y = r * sin(thet);
    x_sect.col(i) << gaudi::vec3(x, y, 0.0);
    circ.push_back(gaudi::vec3(x, y, 0.0));
  }

  std::vector<gaudi::vec3> verts;
  // gaudi::quat q0 = R.calc_frenet_frame(0);
  //  gaudi::quat q0 = gaudi::quat::Identity();

  for (int i = 0; i < R.corner_count(); i++) {
    // q0.normalize();
    // gaudi::quat qi = R.get_frenet(i);
    // if (R.next(i) < 0)
    //  continue;

    gaudi::vec3 x0 = x[i];
    auto idx = R.consec(i);

    gaudi::quat qi = R.__u[i];
    // gaudi::quat dq = qi * q0;

    Eigen::Matrix3d Q = qi.toRotationMatrix();
    /*
    gg::geometry_logger::line(x0, x0 + 0.025 * Q.col(0),
                              vec4(1.0, 0.0, 0.0, 1.0));
    gg::geometry_logger::line(x0, x0 + 0.025 * Q.col(1),
                              vec4(0.0, 1.0, 0.0, 1.0));
    gg::geometry_logger::line(x0, x0 + 0.025 * Q.col(2),
                              vec4(0.0, 0.0, 1.0, 1.0));
    */
    matX x_sect_r = Q * x_sect;

    for (int j0 = 0; j0 < x_sect_r.cols(); j0++) {
      int j1 = (j0 + 1) % x_sect_r.cols();
      gaudi::vec3 c0 = x0 + x_sect_r.col(j0);
      gaudi::vec3 c1 = x0 + x_sect_r.col(j1);
      // gg::geometry_logger::line(c0, c1, vec4(0.0, 0.75, 0.95, 1.0));
      verts.push_back(c0);
    }

    // q0 = qi * dq;
  }

  std::vector<std::vector<int>> faces;
  for (int i0 = 0; i0 < R.corner_count(); i0++) {
    if (R.next(i0) < 0)
      continue;
    int i1 = R.next(i0);
    for (int j0 = 0; j0 < Nc; j0++) {
      int j1 = (j0 + 1) % Nc;
      faces.push_back({Nc * i0 + j0, Nc * i0 + j1, Nc * i1 + j1});
      faces.push_back({Nc * i0 + j0, Nc * i1 + j1, Nc * i1 + j0});
      // faces.push_back({Nc * i0 + j1, Nc * i0 + j0, Nc * i1 + j1});
    }
  }

  obj->fillBuffer([&](gg::BufferObject &o) -> void {
    o.allocateVerts(faces.size(), verts.size());
    auto &indices = o.indices();
    auto &positions = o.positions();
    auto &colors = o.colors();
    for (int i = 0; i < faces.size(); i++) {
      for (int j = 0; j < faces[i].size(); j++) {
        indices.col(i)[j] = faces[i][j];
      }
    }

    for (int i = 0; i < verts.size(); i++) {
      for (int j = 0; j < 3; j++) {

        positions.col(i)[j] = verts[i][j];
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