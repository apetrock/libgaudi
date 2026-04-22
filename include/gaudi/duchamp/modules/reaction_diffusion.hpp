#ifndef __DUCHAMP_RXDIFF_MODULE__
#define __DUCHAMP_RXDIFF_MODULE__

#include "gaudi/asawa/datums.hpp"
#include "gaudi/bontecou/laplacian.hpp"
#include "module_base_shell.hpp"
#include <optional>
#include <vector>

namespace gaudi {
namespace duchamp {
class reaction_diffusion : public module_base_shell {
public:
  DEFINE_CREATE_FUNC(reaction_diffusion)
  reaction_diffusion(asawa::shell::shell::ptr M, real f, real k, real da,
                     real db)
      : module_base_shell(M), _f(f), _k(k), _da(da), _db(db) {
    init_rx();
  };

  virtual ~reaction_diffusion() {}

  void set_diffuse_smooth(std::optional<rx_smooth_fn> fn) { _smooth = std::move(fn); }
  void set_diffuse_input_scale(std::optional<real> s) { _input_scale = s; }
  void set_effect_coeff_field(const std::vector<real> *p) { _effect_coeff = p; }

  // pretty good library function, actually,
  index_t _init_datum() {
    return gaudi::asawa::init_vert_datum<real>(*_M, 0.0);
  }

  void init_rx() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0f, 1.0f);

    _irxa = _init_datum();
    _irxb = _init_datum();
    std::vector<real> &rxa = asawa::get_real_data(*_M, _irxa);
    std::vector<real> &rxb = asawa::get_real_data(*_M, _irxb);
    for (int i = 0; i < rxa.size(); i++) {
      real ta = dis(gen);
      real tb = dis(gen);
      // rxb[i] = 1.0;
      rxa[i] = 1.0;
      rxb[i] = 0.0;

      if (tb > 0.975) {
        rxb[i] = 1.0;
      }
    }
  }

  /// \p effect_coeff: optional per-vertex scale on `f` and `k` in the reaction
  /// (nullptr uses \ref set_effect_coeff_field pointer, else all ones).
  virtual void step_anisotropic(real h, //
                                const std::vector<real> &f,
                                const std::vector<real> &k,
                                const std::vector<real> *effect_coeff = nullptr) {

    std::vector<vec3> &x = asawa::get_vec_data(*_M, 0);
    std::vector<real> &rxa = get_rxa();
    std::vector<real> &rxb = get_rxb();
    const std::vector<real> *ec =
        effect_coeff ? effect_coeff : _effect_coeff;
    const int n = static_cast<int>(rxa.size());

    for (int i = 0; i < n; i++) {
      real c = rx_effect_coeff_at(ec, i, n);
      real fi = c * f[i];
      real ki = c * k[i];
      auto [rxai, rxbi] =
          bontecou::grey_scott_2({rxa[i], rxb[i]}, fi, ki, h);
      rxa[i] = rxai, rxb[i] = rxbi;
    }

    rx_apply_smooth(_M, x, _smooth, rxa, h * _da, _input_scale);
    rx_apply_smooth(_M, x, _smooth, rxb, h * _db, _input_scale);
  }

  virtual void step_isotropic(real h, real fi, real ki) {
    std::vector<real> f(get_rxa().size(), fi);
    std::vector<real> k(get_rxa().size(), ki);
    step_anisotropic(h, f, k, nullptr);
  }

  virtual void step(real h) { step_isotropic(h, _f, _k); }

  std::vector<real> &get_rxa() { return asawa::get_real_data(*_M, _irxa); }
  std::vector<real> &get_rxb() { return asawa::get_real_data(*_M, _irxb); }

  index_t _irxa = -1, _irxb = -1;

  real _f = 0.025, _k = 0.535;
  real _da = 5.00e-4, _db = 0.4 * _da;

  std::optional<rx_smooth_fn> _smooth;
  std::optional<real> _input_scale;
  const std::vector<real> *_effect_coeff = nullptr;
};

} // namespace duchamp
} // namespace gaudi

#endif
