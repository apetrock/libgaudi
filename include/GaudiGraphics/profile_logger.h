/*
 *  manifold_singleton.h
 *  Phase Vocoder
 *
 *  Created by John Delaney on 12/29/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>

#ifndef __GG_PROFILER__
#define __GG_PROFILER__

namespace gg {

class timer {

  timer(geometry_logger *logger) { _logger = logger; }
  timer ~() { _logger = NULL; }

  geometry_logger *_logger;
};

class profile_logger {

public:
  static profile_logger &get_instance();

  static void clear();
  static void frame();

  bool &initialized() { return instance_flag; }
  bool initialized() const { return instance_flag; }

private:
  profile_logger() { global_instance = this; }

  profile_logger(const profile_logger &);
  profile_logger &operator=(const profile_logger &);

  static profile_logger *global_instance;
  static bool instance_flag;
};
} // namespace gg

#endif
