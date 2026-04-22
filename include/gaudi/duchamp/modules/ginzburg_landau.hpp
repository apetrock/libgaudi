#ifndef __GAUDI_DUCHAMP_GINZBURG_LANDAU__
#define __GAUDI_DUCHAMP_GINZBURG_LANDAU__

#include "gaudi/asawa/datums.hpp"
#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/bontecou/laplacian.hpp"
#include "gaudi/duchamp/modules/module_base_shell.hpp"
#include <Eigen/Sparse>
#include <cmath>
#include <optional>
#include <vector>

namespace gaudi {
namespace duchamp {

/// Complex Ginzburg–Landau on two real components `u = Re A`, `v = Im A`:
/// `∂t A = A + (1+iα) Δ A − (1+iβ) |A|² A` with explicit split:
/// (1) forward Euler on `G = A − (1+iβ)|A|²A` (real/imag as in plan),
/// (2) explicit `(1+iα)Δ` via cotan stiffness `multC` (optionally custom sparse C).
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

  /// Optional anisotropic (or otherwise custom) cotan replacement for `multC`.
  void set_dispersive_stiffness(const sparse_mat *C) { _dispersive_C = C; }

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

  void step(real h, const std::vector<real> &alpha,
            const std::vector<real> &beta,
            const std::vector<real> *effect_coeff = nullptr) {

    std::vector<vec3> &x = asawa::get_vec_data(*_M, 0);
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
      real r2 = ui * ui + vi * vi;
      real gu = ui - r2 * (ui - b * vi);
      real gv = vi - r2 * (vi + b * ui);
      u[i] = ui + h * c * gu;
      v[i] = vi + h * c * gv;
    }

    std::vector<real> Du = mult_c_on_vertices(_M, x, u, _dispersive_C);
    std::vector<real> Dv = mult_c_on_vertices(_M, x, v, _dispersive_C);

    for (int i = 0; i < n; ++i) {
      real c = rx_effect_coeff_at(ec, i, n);
      real a = i < static_cast<int>(alpha.size()) ? alpha[i] : _alpha0;
      u[i] += h * c * (Du[i] - a * Dv[i]);
      v[i] += h * c * (Dv[i] + a * Du[i]);
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

  virtual void step(real h) override {
    const int n = _M->vert_count();
    std::vector<real> alpha(n, _alpha0), beta(n, _beta0);
    step(h, alpha, beta, nullptr);
  }

private:
  static std::vector<real>
  mult_c_on_vertices(asawa::shell::shell::ptr M, const std::vector<vec3> &x,
                     const std::vector<real> &f, const sparse_mat *override_C) {
    bontecou::laplacian L(M, x);
    if (override_C != nullptr &&
        static_cast<index_t>(override_C->rows()) == M->vert_count() &&
        static_cast<index_t>(override_C->cols()) == M->vert_count())
      L.set_stiffness(*override_C);
    std::vector<real> fc = asawa::shell::compress_to_vert_range<real>(*M, f);
    std::vector<real> m = L.multC(fc);
    return asawa::shell::expand_from_vert_range<real>(*M, m);
  }

  index_t _iu = -1, _iv = -1;
  real _alpha0 = 1.5;
  real _beta0 = 1.0;
  const sparse_mat *_dispersive_C = nullptr;
  const std::vector<real> *_effect_coeff = nullptr;
};

} // namespace duchamp
} // namespace gaudi

#endif
