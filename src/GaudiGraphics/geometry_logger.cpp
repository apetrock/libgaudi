/*
 *  manifold_singleton.cpp
 *  Phase Vocoder
 *
 *  Created by John Delaney on 12/29/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

//#include <gl/glut.h>
#include "GaudiGraphics/geometry_logger.h"
#include <numeric>

namespace gg {

//////////////////////////
// debugger
//////////////////////////

bool geometry_logger::instance_flag = false;
geometry_logger *geometry_logger::global_instance = NULL;

geometry_logger &geometry_logger::get_instance() {
  static geometry_logger logger;
  if (!logger.initialized()) {
    logger.initialized() = true;
  }

  return logger;
}

void geometry_logger::render() {
  geometry_logger &logger = geometry_logger::get_instance();
  logger.debugLines->renderLines();
}

void geometry_logger::clear() {
  geometry_logger &logger = geometry_logger::get_instance();
  logger.debugLines->clear();
}

void geometry_logger::line4(const Vec4 &p0, const Vec4 &p1, const Vec4 &color) {
  geometry_logger &logger = geometry_logger::get_instance();
  logger.debugLines->pushLine(p0, p1, color);
}
void geometry_logger::line(const Vec3d &p0, const Vec3d &p1,
                           const Vec4d &color) {
  geometry_logger &logger = geometry_logger::get_instance();
  logger.debugLines->pushLine(xyzw(p0), xyzw(p1),
                              Vec4(color[0], color[1], color[2], color[3]));
}

/*
template <typename SPACE>
void lines(const std::vector<Vec3d> &p0, const std::vector<Vec3d> &p1,
           const Vec4d &color, double D = 0.05) {

  double mx = std::accumulate(p1.begin(), p1.end(), 0.0, [](double a, auto &c) {
    return max(a, c.norm());
  });

  for (int i = 0; i < p0.size(); i++) {

    const auto &p = p0[i];
    const auto &a = p1[i];

    auto pa = p + D * a / mx;
    gg::geometry_logger::line(Vec4(p0[0], p0[1], p0[2], 1.0),
                              Vec4(p1[0], p1[1], p1[2], 1.0),
                              Vec4(c[0], c[1], c[2], 1.0));
  }
}
*/

void geometry_logger::frame(Mat3d M, Vec3d c, double C = 1.0) {

  Vec3d t0 = M.block(0, 0, 3, 1);
  Vec3d t1 = M.block(0, 1, 3, 1);
  Vec3d t2 = M.block(0, 2, 3, 1);

  gg::geometry_logger::line(c - C * t0, c + C * t0, Vec4d(1.0, 0.0, 0.0, 1.0));
  gg::geometry_logger::line(c - C * t1, c + C * t1, Vec4d(0.0, 1.0, 0.0, 1.0));
  gg::geometry_logger::line(c - C * t2, c + C * t2, Vec4d(0.0, 0.0, 1.0, 1.0));
}

} // namespace gg
