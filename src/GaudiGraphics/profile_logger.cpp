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

} // namespace gg
