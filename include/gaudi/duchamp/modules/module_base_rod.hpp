
#ifndef __DUCHAMP_MODULE_ROD_BASE__
#define __DUCHAMP_MODULE_ROD_BASE__
#include "gaudi/asawa/rod/rod.hpp"
#include "module_base.hpp"

namespace gaudi {
namespace duchamp {

class module_base_rod : public module_base {
public:
  DEFINE_CREATE_FUNC(module_base_rod)
  module_base_rod(asawa::rod::rod::ptr R) : _R(R){};
  virtual ~module_base_rod(){};
  asawa::rod::rod::ptr _R;
};

} // namespace duchamp
} // namespace gaudi
#endif
