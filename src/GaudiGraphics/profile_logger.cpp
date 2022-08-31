/*
 *  manifold_singleton.cpp
 *  Phase Vocoder
 *
 *  Created by John Delaney on 12/29/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

//#include <gl/glut.h>
#include "GaudiGraphics/profile_logger.h"
#include <numeric>

namespace gg {

//////////////////////////
// debugger
//////////////////////////

bool profile_logger::instance_flag = false;
profile_logger *geometry_logger::global_instance = NULL;

profile_logger &profile_logger::get_instance() {
  static profile_logger logger;
  if (!logger.initialized()) {
    logger.initialized() = true;
  }

  return logger;
}

void geometry_logger::clear() {}

void geometry_logger::frame() {

  // gg::geometry_logger::line(c - C * t0, c + C * t0, Vec4d(1.0, 0.0,
  // 0.0, 1.0)); gg::geometry_logger::line(c - C * t1, c + C * t1,
  // Vec4d(0.0, 1.0, 0.0, 1.0)); gg::geometry_logger::line(c - C * t2, c + C *
  // t2, Vec4d(0.0, 0.0, 1.0, 1.0));
}

} // namespace gg
