#ifndef __GAUDI_DUCHAMP_SWIFT_HOHENBERG__
#define __GAUDI_DUCHAMP_SWIFT_HOHENBERG__

#include "gaudi/asawa/datums.hpp"
#include "gaudi/duchamp/modules/module_base_shell.hpp"
#include <cmath>
#include <optional>
#include <vector>

namespace gaudi {
namespace duchamp {

/// Scalar Swift–Hohenberg–style split step (pattern driver):
/// reaction `∂u ≈ h·c·(ε u − g u³)` then two implicit cotan (or custom) smooths
/// with timestep `lap_smooth_dt` each to approximate `(1+λΔ)²` in operator-split form.
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

  std::vector<real> &get_u() { return asawa::get_real_data(*_M, _iu); }

  /// \p epsilon, \p g_cubic, \p lap_smooth_dt per vertex (same layout as shell verts).
  /// Each implicit smooth uses timestep `lap_smooth_dt[i]` averaged over the edge
  /// endpoints of the vertex star — here we use vertex value directly.
  void step(real h, const std::vector<real> &epsilon,
            const std::vector<real> &g_cubic,
            const std::vector<real> &lap_smooth_dt,
            const std::vector<real> *effect_coeff = nullptr) {

    std::vector<vec3> &x = asawa::get_vec_data(*_M, 0);
    std::vector<real> &u = get_u();
    const std::vector<real> *ec = effect_coeff ? effect_coeff : _effect_coeff;
    const int n = static_cast<int>(u.size());

    for (int i = 0; i < n; ++i) {
      real c = rx_effect_coeff_at(ec, i, n);
      real eps = i < static_cast<int>(epsilon.size()) ? epsilon[i] : _eps0;
      real g = i < static_cast<int>(g_cubic.size()) ? g_cubic[i] : _g0;
      real ui = u[i];
      if (!std::isfinite(ui))
        ui = 0.0;
      ui = std::max(-50.0, std::min(50.0, ui));
      u[i] = ui + h * c * (eps * ui - g * ui * ui * ui);
    }

    for (int pass = 0; pass < 2; ++pass) {
      // One global smooth per pass: use mean(lap_smooth_dt) so implicit solve stays well-posed.
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
      rx_apply_smooth(_M, x, _smooth, u, mean_dt, _input_scale);
    }

    for (int i = 0; i < n; ++i) {
      if (!std::isfinite(u[i]))
        u[i] = 0.0;
      u[i] = std::max(-50.0, std::min(50.0, u[i]));
    }
  }

  virtual void step(real h) override {
    const int n = _M->vert_count();
    std::vector<real> eps(n, _eps0), g(n, _g0), lap(n, 1e-3);
    step(h, eps, g, lap, nullptr);
  }

private:
  index_t _iu = -1;
  real _eps0 = 0.04;
  real _g0 = 1.0;
  std::optional<rx_smooth_fn> _smooth;
  std::optional<real> _input_scale;
  const std::vector<real> *_effect_coeff = nullptr;
};

} // namespace duchamp
} // namespace gaudi

#endif
