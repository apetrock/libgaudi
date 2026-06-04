#ifndef __GAUDI_DUCHAMP_RX_GINZBURG_LANDAU_HPP__
#define __GAUDI_DUCHAMP_RX_GINZBURG_LANDAU_HPP__

#include "gaudi/asawa/datums.hpp"
#include "gaudi/duchamp/modules/module_base_shell.hpp"
#include "gaudi/duchamp/modules/rx/rx_pipeline.hpp"
#include "gaudi/kusama/complex_laplacian.hpp"
#include "gaudi/kusama/laplacian.hpp"
#include "gaudi/kusama/rx/ginzburg_landau.hpp"
#include <cmath>
#include <optional>
#include <vector>

namespace gaudi {
namespace duchamp {
namespace rx {

enum class gl_reaction_step_mode { forward_euler, newton };

/// Complex Ginzburg--Landau: `u = Re A`, `v = Im A` with operator split
/// (reaction, then implicit (1+iα) Δ via @ref kusama::rx::cgle::linear_operator).
class ginzburg_landau : public module_base_shell {
public:
  using sparse_mat = Eigen::SparseMatrix<real>;

  DEFINE_CREATE_FUNC(ginzburg_landau)

  ginzburg_landau(asawa::shell::shell::ptr M, real alpha0 = 1.5, real beta0 = 1.0)
      : module_base_shell(M), _alpha0(alpha0), _beta0(beta0) {
    _iu = asawa::init_vert_datum<real>(*_M, 0.0);
    _iv = asawa::init_vert_datum<real>(*_M, 0.0);
    std::vector<real> &u = get_u();
    std::vector<real> &v = get_v();
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    for (int i = 0; i < u.size(); ++i) {
      u[i] = 0.1 * dis(gen);
      v[i] = 0.1 * dis(gen);
    }
  }

  void set_dispersive_stiffness(const sparse_mat *C) { _dispersive_C = C; }
  void set_cgle_linear_config(kusama::rx::cgle::linear_config c) {
    _cgle_cfg = std::move(c);
  }
  void set_reaction_mode(gl_reaction_step_mode m) { _reaction_mode = m; }
  void set_effect_coeff_field(const std::vector<real> *p) { _effect_coeff = p; }

  std::vector<real> &get_u() { return asawa::get_real_data(*_M, _iu); }
  std::vector<real> &get_v() { return asawa::get_real_data(*_M, _iv); }

  std::vector<real> amplitude_abs() {
    const std::vector<real> &u = get_u();
    const std::vector<real> &v = get_v();
    std::vector<real> out(u.size());
    for (std::size_t i = 0; i < u.size(); ++i)
      out[i] = std::sqrt(u[i] * u[i] + v[i] * v[i]);
    return out;
  }

  std::vector<real> amplitude_phase() {
    const std::vector<real> &u = get_u();
    const std::vector<real> &v = get_v();
    std::vector<real> out(u.size());
    for (std::size_t i = 0; i < u.size(); ++i)
      out[i] = std::atan2(v[i], u[i]);
    return out;
  }

  void step_reaction(const std::vector<real> &alpha, const std::vector<real> &beta,
                     real h, const std::vector<real> *effect_coeff) {
    (void)alpha;
    std::vector<real> &u = get_u();
    std::vector<real> &v = get_v();
    const std::vector<real> *ec = effect_coeff ? effect_coeff : _effect_coeff;
    const int n = static_cast<int>(u.size());
    for (int i = 0; i < n; ++i) {
      real c = rx_effect_coeff_at(ec, i, n);
      real b = i < static_cast<int>(beta.size()) ? beta[i] : _beta0;
      real ui = u[i], vi = v[i];
      if (!std::isfinite(ui))
        ui = 0.0;
      if (!std::isfinite(vi))
        vi = 0.0;
      ui = std::max(-50.0, std::min(50.0, ui));
      vi = std::max(-50.0, std::min(50.0, vi));
      if (_reaction_mode == gl_reaction_step_mode::newton) {
        const real u0 = ui, v0 = vi;
        kusama::rx::ginzburg_landau::reaction_newton_2d(ui, vi, u0, v0, c, b, h);
        u[i] = ui;
        v[i] = vi;
      } else {
        real r2 = ui * ui + vi * vi;
        real gu = ui - r2 * (ui - b * vi);
        real gv = vi - r2 * (vi + b * ui);
        u[i] = ui + h * c * gu;
        v[i] = vi + h * c * gv;
      }
    }
    for (int i = 0; i < n; ++i) {
      if (!std::isfinite(u[i]))
        u[i] = 0.0;
      if (!std::isfinite(v[i]))
        v[i] = 0.0;
      u[i] = std::max(-50.0, std::min(50.0, u[i]));
      v[i] = std::max(-50.0, std::min(50.0, v[i]));
    }
  }

  void step_linear(const std::vector<real> &alpha, real h,
                   const std::vector<real> *effect_coeff) {
    (void)effect_coeff;
    std::vector<vec3> &x = asawa::get_vec_data(*_M, 0);
    std::vector<real> &u = get_u();
    std::vector<real> &v = get_v();
    ::gaudi::kusama::laplacian L(_M, x);
    if (_dispersive_C != nullptr &&
        static_cast<index_t>(_dispersive_C->rows()) == _M->vert_count() &&
        static_cast<index_t>(_dispersive_C->cols()) == _M->vert_count()) {
      L.set_stiffness(*_dispersive_C);
    }
    kusama::rx::cgle::linear_operator::apply(L, alpha, h, u, v, _cgle_cfg);
  }

  void step(real h, const std::vector<real> &alpha, const std::vector<real> &beta,
            const std::vector<real> *effect_coeff = nullptr) {
    step_reaction(alpha, beta, h, effect_coeff);
    step_linear(alpha, h, effect_coeff);
  }

  virtual void step(real h) override {
    const int n = static_cast<int>(_M->vert_count());
    std::vector<real> alpha(n, _alpha0), beta(n, _beta0);
    step(h, alpha, beta, nullptr);
  }

  rx_pipeline make_rx_pipeline(real h, const std::vector<real> &alpha,
                               const std::vector<real> &beta,
                               const std::vector<real> *effect_coeff) {
    rx_pipeline p;
    p.push_back([this, h, &alpha, &beta, effect_coeff]() {
      step_reaction(alpha, beta, h, effect_coeff);
    });
    p.push_back([this, h, &alpha, effect_coeff]() {
      step_linear(alpha, h, effect_coeff);
    });
    return p;
  }

private:
  index_t _iu = -1, _iv = -1;
  real _alpha0 = 1.5;
  real _beta0 = 1.0;
  const sparse_mat *_dispersive_C = nullptr;
  const std::vector<real> *_effect_coeff = nullptr;
  kusama::rx::cgle::linear_config _cgle_cfg;
  gl_reaction_step_mode _reaction_mode = gl_reaction_step_mode::forward_euler;
};

} // namespace rx
} // namespace duchamp
} // namespace gaudi

#endif
