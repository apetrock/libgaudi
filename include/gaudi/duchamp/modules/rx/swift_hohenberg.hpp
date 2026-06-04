#ifndef __GAUDI_DUCHAMP_RX_SWIFT_HOHENBERG_HPP__
#define __GAUDI_DUCHAMP_RX_SWIFT_HOHENBERG_HPP__

#include "gaudi/asawa/datums.hpp"
#include "gaudi/duchamp/modules/module_base_shell.hpp"
#include "gaudi/duchamp/modules/rx/detail/rx_step_modules.hpp"
#include "gaudi/duchamp/modules/rx/rx_pipeline.hpp"
#include "gaudi/kusama/rx/swift_hohenberg.hpp"
#include <optional>
#include <vector>

namespace gaudi {
namespace duchamp {
namespace rx {

using sh_reaction_step_mode = kusama::rx::swift_hohenberg::reaction_step_mode;

/// Scalar Swift–Hohenberg–style split step (pattern driver).
class swift_hohenberg : public module_base_shell {
public:
  DEFINE_CREATE_FUNC(swift_hohenberg)

  swift_hohenberg(asawa::shell::shell::ptr M, real eps0 = 0.04, real g0 = 1.0)
      : module_base_shell(M), _eps0(eps0), _g0(g0) {
    _iu = asawa::init_vert_datum<real>(*_M, 0.0);
    std::vector<real> &u = get_u();
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    for (int i = 0; i < u.size(); ++i) {
      u[i] = 0.1 * dis(gen);
      if (dis(gen) > 0.97)
        u[i] += 0.5;
    }
  }

  void set_diffuse_smooth(std::optional<rx_smooth_fn> fn) { _smooth = std::move(fn); }
  void set_diffuse_input_scale(std::optional<real> s) { _input_scale = s; }
  void set_effect_coeff_field(const std::vector<real> *p) { _effect_coeff = p; }
  void set_reaction_mode(sh_reaction_step_mode m) { _reaction_mode = m; }

  std::vector<real> &get_u() { return asawa::get_real_data(*_M, _iu); }

  void step(real h, const std::vector<real> &epsilon,
            const std::vector<real> &g_cubic,
            const std::vector<real> &lap_smooth_dt,
            const std::vector<real> *effect_coeff = nullptr) {

    detail::rx_cotan_diffuser diff{_M, _smooth, _input_scale};
    detail::swift_hohenberg_stepper<detail::rx_cotan_diffuser> stepper(
        _M, _iu, std::move(diff), _reaction_mode, _eps0, _g0);
    stepper.step(h, epsilon, g_cubic, lap_smooth_dt, effect_coeff, _effect_coeff);
  }

  virtual void step(real h) override {
    const int n = static_cast<int>(_M->vert_count());
    std::vector<real> eps(n, _eps0), g(n, _g0), lap(n, 1e-3);
    step(h, eps, g, lap, nullptr);
  }

  rx_pipeline make_rx_pipeline(
      real h, const std::vector<real> &epsilon, const std::vector<real> &g_cubic,
      const std::vector<real> &lap_smooth_dt, const std::vector<real> *effect_coeff) {
    rx_pipeline p;
    p.push_back([this, h, &epsilon, &g_cubic, &lap_smooth_dt, effect_coeff]() {
      step(h, epsilon, g_cubic, lap_smooth_dt, effect_coeff);
    });
    return p;
  }

private:
  index_t _iu = -1;
  real _eps0 = 0.04;
  real _g0 = 1.0;
  sh_reaction_step_mode _reaction_mode = sh_reaction_step_mode::newton;
  std::optional<rx_smooth_fn> _smooth;
  std::optional<real> _input_scale;
  const std::vector<real> *_effect_coeff = nullptr;
};

} // namespace rx
} // namespace duchamp
} // namespace gaudi

#endif
