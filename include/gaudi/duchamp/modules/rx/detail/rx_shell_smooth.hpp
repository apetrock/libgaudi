#ifndef GAUDI_DUCHAMP_RX_DETAIL_RX_SHELL_SMOOTH_HPP
#define GAUDI_DUCHAMP_RX_DETAIL_RX_SHELL_SMOOTH_HPP

#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/kusama/laplacian.hpp"
#include <functional>
#include <optional>
#include <vector>

namespace gaudi {
namespace duchamp {

inline void diffuse(asawa::shell::shell::ptr M, const std::vector<vec3> &x,
                    std::vector<real> &f, real dt) {

  kusama::laplacian L(M, x);
  std::vector<real> f_comp = asawa::shell::compress_to_vert_range<real>(*M, f);

  std::vector<real> d = L.diffuse2(f_comp, dt);
  std::vector<real> d_exp = asawa::shell::expand_from_vert_range<real>(*M, d);

  for (int k = 0; k < f.size(); k++) {
    f[k] = d_exp[k];
  }
}

/// One implicit cotan Crank–Nicolson–style diffusion step on vertex scalars.
using rx_smooth_fn = std::function<void(std::vector<real> &field, real dt)>;

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

inline real rx_effect_coeff_at(const std::vector<real> *effect_coeff, int i,
                               int n) {
  (void)n;
  if (!effect_coeff || i < 0 || i >= static_cast<int>(effect_coeff->size()))
    return 1.0;
  return (*effect_coeff)[i];
}

} // namespace duchamp
} // namespace gaudi

#endif
