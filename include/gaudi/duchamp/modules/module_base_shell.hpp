
#ifndef __DUCHAMP_MODULE_ROD_BASE__
#define __DUCHAMP_MODULE_ROD_BASE__

#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/bontecou/laplacian.hpp"

#include "module_base.hpp"

namespace gaudi {
namespace duchamp {

void diffuse(asawa::shell::shell::ptr M, const std::vector<vec3> &x,
             std::vector<real> &f, real dt) {

  bontecou::laplacian L(M, x);
  std::vector<real> f_comp = asawa::shell::compress_to_vert_range<real>(*M, f);

  std::vector<real> d = L.diffuse2(f_comp, dt);
  std::vector<real> d_exp = asawa::shell::expand_from_vert_range<real>(*M, d);

  for (int k = 0; k < f.size(); k++) {
    f[k] = d_exp[k];
  }
}

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
