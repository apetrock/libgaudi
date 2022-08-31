/*
 *  manifold_singleton.h
 *  Phase Vocoder
 *
 *  Created by John Delaney on 12/29/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "buffers.hpp"
#include <iostream>

#ifndef __GG_DEBUGGER__
#define __GG_DEBUGGER__

namespace gg {
class geometry_logger;

class geometry_logger {

public:
  gg::DebugBufferPtr debugLines = NULL;

  static geometry_logger &get_instance();
  static void render();
  static void clear();
  static void line4(const Vec4 &p0, const Vec4 &p1, const Vec4 &color);
  static void line(const Vec3d &p0, const Vec3d &p1, const Vec4d &color);
  void frame(Mat3d M, Vec3d c, double C);

  bool &initialized() { return instance_flag; }
  bool initialized() const { return instance_flag; }

private:
  geometry_logger() {
    global_instance = this;

    debugLines = gg::DebugBuffer::create();
    debugLines->init();
  }

  geometry_logger(const geometry_logger &);
  geometry_logger &operator=(const geometry_logger &);

  static geometry_logger *global_instance;
  static bool instance_flag;
};
} // namespace gg

#endif
