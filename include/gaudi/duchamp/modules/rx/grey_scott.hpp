#ifndef __GAUDI_DUCHAMP_RX_GREY_SCOTT_HPP__
#define __GAUDI_DUCHAMP_RX_GREY_SCOTT_HPP__

#include "gaudi/asawa/datums.hpp"
#include "gaudi/duchamp/modules/module_base_shell.hpp"
#include "gaudi/duchamp/modules/rx/detail/rx_step_modules.hpp"
#include "gaudi/duchamp/modules/rx/rx_pipeline.hpp"
#include "gaudi/kusama/laplacian.hpp"
#include <optional>
#include <vector>

namespace gaudi {
namespace duchamp {
namespace rx {

/// Grey–Scott on a shell. Default \p _f, \p _k sit in a well-studied spot of the
/// *f*–*k* parameter plane (spots / labyrinthine-type behavior; see e.g. Pearson,
/// *Pattern formation* / graphics references for full phase diagrams).
class grey_scott : public module_base_shell {
public:
  DEFINE_CREATE_FUNC(grey_scott)
  grey_scott(asawa::shell::shell::ptr M, real f, real k, real da, real db)
      : module_base_shell(M), _f(f), _k(k), _da(da), _db(db) {
    init_rx();
  };

  virtual ~grey_scott() {}

  void set_diffuse_smooth(std::optional<rx_smooth_fn> fn) { _smooth = std::move(fn); }
  void set_diffuse_input_scale(std::optional<real> s) { _input_scale = s; }
  void set_effect_coeff_field(const std::vector<real> *p) { _effect_coeff = p; }

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
      real tb = dis(gen);
      rxa[i] = 1.0;
      rxb[i] = 0.0;

      if (tb > 0.975) {
        rxb[i] = 1.0;
      }
    }
  }

  virtual void step_anisotropic(real h, //
                                const std::vector<real> &f,
                                const std::vector<real> &k,
                                const std::vector<real> *effect_coeff = nullptr) {

    detail::rx_cotan_diffuser diff{_M, _smooth, _input_scale};
    detail::grey_scott_stepper<detail::rx_cotan_diffuser> stepper(
        _M, _irxa, _irxb, std::move(diff));
    stepper.step_anisotropic(h, f, k, _da, _db, effect_coeff, _effect_coeff);
  }

  virtual void step_isotropic(real h, real fi, real ki) {
    std::vector<real> f(get_rxa().size(), fi);
    std::vector<real> k(get_rxa().size(), ki);
    step_anisotropic(h, f, k, nullptr);
  }

  virtual void step(real h) { step_isotropic(h, _f, _k); }

  rx_pipeline make_step_pipeline(real h, const std::vector<real> &f,
                                 const std::vector<real> &k,
                                 const std::vector<real> *effect_coeff) {
    rx_pipeline p;
    p.push_back([=, this]() { step_anisotropic(h, f, k, effect_coeff); });
    return p;
  }

  std::vector<real> &get_rxa() { return asawa::get_real_data(*_M, _irxa); }
  std::vector<real> &get_rxb() { return asawa::get_real_data(*_M, _irxb); }

  index_t _irxa = -1, _irxb = -1;

  real _f = 0.025, _k = 0.535;
  real _da = 5.00e-4, _db = 0.4 * _da;

  std::optional<rx_smooth_fn> _smooth;
  std::optional<real> _input_scale;
  const std::vector<real> *_effect_coeff = nullptr;
};

} // namespace rx
} // namespace duchamp
} // namespace gaudi

#endif
