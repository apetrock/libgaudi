#ifndef __GAUDI_DUCHAMP_RX_DETAIL_RX_STEP_MODULES_HPP__
#define __GAUDI_DUCHAMP_RX_DETAIL_RX_STEP_MODULES_HPP__

#include "gaudi/duchamp/modules/module_base_shell.hpp"
#include "gaudi/kusama/rx/grey_scott.hpp"
#include "gaudi/kusama/rx/swift_hohenberg.hpp"
#include <cmath>
#include <optional>
#include <vector>

namespace gaudi {
namespace duchamp {
namespace rx {
namespace detail {

/// Shared scalar diffusion contract for Grey–Scott and Swift–Hohenberg steppers.
struct rx_cotan_diffuser {
  asawa::shell::shell::ptr M;
  std::optional<rx_smooth_fn> smooth;
  std::optional<real> input_scale;

  void operator()(std::vector<real> &field, real dt) const {
    std::vector<vec3> &x = asawa::get_vec_data(*M, 0);
    rx_apply_smooth(M, x, smooth, field, dt, input_scale);
  }
};

/// Grey–Scott: local kinetics, then the same @p Diffuser on `u` and `v`.
template <class Diffuser>
class grey_scott_stepper {
public:
  grey_scott_stepper(asawa::shell::shell::ptr M, index_t u_idx, index_t v_idx,
                     Diffuser diffuse)
      : _M(std::move(M)), _u_idx(u_idx), _v_idx(v_idx), _diffuse(std::move(diffuse)) {}

  void step_anisotropic(real h, const std::vector<real> &f, const std::vector<real> &k,
                        real d_a, real d_b, const std::vector<real> *effect_coeff,
                        const std::vector<real> *effect_field) {

    std::vector<real> &rxa = asawa::get_real_data(*_M, _u_idx);
    std::vector<real> &rxb = asawa::get_real_data(*_M, _v_idx);
    const std::vector<real> *ec = effect_coeff ? effect_coeff : effect_field;
    const int n = static_cast<int>(rxa.size());

    for (int i = 0; i < n; i++) {
      real c = rx_effect_coeff_at(ec, i, n);
      real fi = c * f[i];
      real ki = c * k[i];
      auto p = kusama::rx::grey_scott::trapezoid({rxa[i], rxb[i]}, fi, ki, h);
      rxa[i] = p[0];
      rxb[i] = p[1];
    }
    _diffuse(rxa, h * d_a);
    _diffuse(rxb, h * d_b);
  }

  void step_isotropic(real h, real fi, real ki, real d_a, real d_b,
                      const std::vector<real> *effect_field) {
    const int n = static_cast<int>(asawa::get_real_data(*_M, _u_idx).size());
    std::vector<real> f(n, fi), k(n, ki);
    step_anisotropic(h, f, k, d_a, d_b, nullptr, effect_field);
  }

private:
  asawa::shell::shell::ptr _M;
  index_t _u_idx;
  index_t _v_idx;
  Diffuser _diffuse;
};

/// Swift–Hohenberg: scalar reaction then two @p Diffuser passes.
template <class Diffuser>
class swift_hohenberg_stepper {
public:
  swift_hohenberg_stepper(
      asawa::shell::shell::ptr M, index_t u_idx, Diffuser diffuse,
      kusama::rx::swift_hohenberg::reaction_step_mode mode =
          kusama::rx::swift_hohenberg::reaction_step_mode::newton,
      real eps0 = 0.04, real g0 = 1.0)
      : _M(std::move(M)), _u_idx(u_idx), _diffuse(std::move(diffuse)), _mode(mode),
        _eps0(eps0), _g0(g0) {}

  void set_reaction_mode(kusama::rx::swift_hohenberg::reaction_step_mode m) {
    _mode = m;
  }

  void step(real h, const std::vector<real> &epsilon, const std::vector<real> &g_cubic,
            const std::vector<real> &lap_smooth_dt, const std::vector<real> *effect_coeff,
            const std::vector<real> *effect_field) {

    std::vector<real> &u = asawa::get_real_data(*_M, _u_idx);
    const std::vector<real> *ec = effect_coeff ? effect_coeff : effect_field;
    const int n = static_cast<int>(u.size());

    for (int i = 0; i < n; ++i) {
      real c = rx_effect_coeff_at(ec, i, n);
      real eps = i < static_cast<int>(epsilon.size()) ? epsilon[i] : _eps0;
      real g = i < static_cast<int>(g_cubic.size()) ? g_cubic[i] : _g0;
      real ui = u[i];
      if (!std::isfinite(ui))
        ui = 0.0;
      ui = std::max(-50.0, std::min(50.0, ui));
      if (_mode == kusama::rx::swift_hohenberg::reaction_step_mode::newton) {
        const real u0 = ui;
        u[i] = kusama::rx::swift_hohenberg::reaction_newton_solve(u0, c, eps, g, h);
      } else {
        u[i] = ui + h * c * (eps * ui - g * ui * ui * ui);
      }
    }

    for (int pass = 0; pass < 2; ++pass) {
      real mean_dt = 0.0;
      int cnt = 0;
      for (int i = 0; i < n && i < static_cast<int>(lap_smooth_dt.size()); ++i) {
        mean_dt += lap_smooth_dt[i];
        ++cnt;
      }
      if (cnt > 0)
        mean_dt /= static_cast<real>(cnt);
      else
        mean_dt = 1e-3;
      _diffuse(u, mean_dt);
    }

    for (int i = 0; i < n; ++i) {
      if (!std::isfinite(u[i]))
        u[i] = 0.0;
      u[i] = std::max(-50.0, std::min(50.0, u[i]));
    }
  }

private:
  asawa::shell::shell::ptr _M;
  index_t _u_idx;
  Diffuser _diffuse;
  kusama::rx::swift_hohenberg::reaction_step_mode _mode;
  real _eps0;
  real _g0;
};

} // namespace detail
} // namespace rx
} // namespace duchamp
} // namespace gaudi

#endif
