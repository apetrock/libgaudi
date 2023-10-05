
#ifndef __DUCHAMP_MODULE_BASE__
#define __DUCHAMP_MODULE_BASE__
#include "gaudi/define_create_func.h"
#include <vector>
// bare bones module base class, not sure how these will interact together
namespace gaudi {
namespace duchamp {

class module_base {
public:
  DEFINE_CREATE_FUNC(module_base)
  module_base(){};
  virtual ~module_base(){};
  virtual void step(real h){};

  // this should be virtual, one of the many things a module might return
  // virtual const std::vector<vec3> &forces() { return std::vector<vec3>(); }
};
} // namespace duchamp
} // namespace gaudi
#endif