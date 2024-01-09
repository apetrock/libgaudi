#ifndef __DUCHAMP_RXDIFF_MODULE__
#define __DUCHAMP_RXDIFF_MODULE__

#include "gaudi/asawa/datums.hpp"
#include "gaudi/bontecou/laplacian.hpp"
#include "module_base_shell.hpp"
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

  virtual ~reaction_diffusion(){

  };

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

  virtual void step_anisotropic(real h, //
                                const std::vector<real> &f,
                                const std::vector<real> &k) {

    std::vector<vec3> &x = asawa::get_vec_data(*_M, 0);
    std::vector<real> &rxa = get_rxa();
    std::vector<real> &rxb = get_rxb();

    for (int i = 0; i < rxa.size(); i++) {
      real u0 = rxa[i];
      real v0 = rxb[i];
      vec2 un0 = vec2(u0, v0);
      auto [rxai, rxbi] =
          bontecou::grey_scott_2({rxa[i], rxb[i]}, f[i], k[i], h);
      rxa[i] = rxai, rxb[i] = rxbi;
    }

    diffuse(_M, x, rxa, h * _da);
    diffuse(_M, x, rxb, h * _db);
  }

  virtual void step_isotropic(real h, real fi, real ki) {
    std::vector<vec3> &x = asawa::get_vec_data(*_M, 0);
    std::vector<real> f(x.size(), fi);
    std::vector<real> k(x.size(), ki);
    step_anisotropic(h, f, k);
  }

  virtual void step(real h) { step_isotropic(h, _f, _k); }

  std::vector<real> &get_rxa() { return asawa::get_real_data(*_M, _irxa); }
  std::vector<real> &get_rxb() { return asawa::get_real_data(*_M, _irxb); }

  index_t _irxa = -1, _irxb = -1;

  real _f = 0.025, _k = 0.535;
  real _da = 5.00e-4, _db = 0.4 * _da;
};

} // namespace duchamp
} // namespace gaudi

#endif