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

} // namespace gg
