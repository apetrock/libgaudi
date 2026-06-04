
#ifndef __DUCHAMP_MODULE_ROD_BASE__
#define __DUCHAMP_MODULE_ROD_BASE__

#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/duchamp/modules/rx/detail/rx_shell_smooth.hpp"
#include "module_base.hpp"

namespace gaudi {
namespace duchamp {

class module_base_shell : public module_base {
public:
  DEFINE_CREATE_FUNC(module_base_shell)
  module_base_shell(asawa::shell::shell::ptr M) : _M(M){};
  virtual ~module_base_shell(){};
  asawa::shell::shell::ptr _M;
};

} // namespace duchamp
} // namespace gaudi
#endif
