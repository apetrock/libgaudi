
#ifndef __DUCHAMP_MODULE_ROD_BASE__
#define __DUCHAMP_MODULE_ROD_BASE__

#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/bontecou/laplacian.hpp"

#include "module_base.hpp"
#include <functional>
#include <optional>

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

/// One implicit cotan Crank–Nicolson–style diffusion step on vertex scalars.
/// Callers can substitute an anisotropic Laplacian by capturing their own
/// `laplacian` + `diffuse2` inside the callback.
using rx_smooth_fn = std::function<void(std::vector<real> &field, real dt)>;

/// Apply optional custom smooth; otherwise default cotan `diffuse`.
/// \p input_scale multiplies \p dt when set (Laplacian “input scale”).
inline void rx_apply_smooth(asawa::shell::shell::ptr M,
                            const std::vector<vec3> &x,
                            const std::optional<rx_smooth_fn> &custom_smooth,
                            std::vector<real> &f, real dt,
                            const std::optional<real> &input_scale) {
  real dtm = dt;
  if (input_scale.has_value())
    dtm *= *input_scale;
  if (custom_smooth.has_value() && *custom_smooth)
    (*custom_smooth)(f, dtm);
  else
    diffuse(M, x, f, dtm);
}

/// Per-vertex multiplicative scale; \p n is vertex slot count (e.g. f.size()).
inline real rx_effect_coeff_at(const std::vector<real> *effect_coeff, int i,
                               int n) {
  (void)n;
  if (!effect_coeff || i < 0 || i >= static_cast<int>(effect_coeff->size()))
    return 1.0;
  return (*effect_coeff)[i];
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
