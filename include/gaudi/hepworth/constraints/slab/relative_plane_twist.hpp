#ifndef GAUDI_HEPWORTH_RELATIVE_PLANE_TWIST_HPP
#define GAUDI_HEPWORTH_RELATIVE_PLANE_TWIST_HPP

#include "slab_dual_ring_arap_constraint.hpp"
#include "slab_prism_config.hpp"
#include "slab_prism_operators.hpp"

#include "gaudi/define_create_func.h"

namespace gaudi {
namespace hepworth {
namespace block {

class relative_plane_twist : public slab_dual_ring_arap_constraint {
public:
  DEFINE_CREATE_FUNC(relative_plane_twist)

  relative_plane_twist(const slab_pierce_stencil &stencil,
                       const slab_prism_config &cfg, const real &w,
                       std::vector<sim_block::ptr> blocks)
      : slab_dual_ring_arap_constraint(
            top_ids_from_stencil(stencil), weights_from_stencil(stencil),
            bot_ids_from_stencil(stencil), bot_weights_from_stencil(stencil), w,
            blocks),
        _stencil(stencil), _cfg(cfg), _w_user(w) {}

  virtual std::string name() override { return "relative_plane_twist"; }

  void init_rest(const vecX &q, solver_context &ctx) override {
    (void)ctx;
    if (_rest_initialized)
      return;
    const auto pos_at = [&](index_t vi) { return _blocks[0]->get_vec3(vi, q); };
    init_spoke_rest(pos_at);
    _rest = build_slab_prism_rest(_stencil, pos_at);
    _w = effective_slab_weight(_w_user, _rest.L_rest, _cfg.L_min, _cfg.w_eff_max);
    _rest_initialized = true;
  }

  virtual void project(const vecX &q, vecX &p) override {
    const vec3 q0 = _blocks[0]->get_vec3(_stencil.foot, q);
    mat3 R0 = fit_ring_procrustes(q0, _top_ids, _top_weights, _d_top_rest, q);

    const vec3 c_bot = _blocks[0]->get_vec3(_stencil.pierce_foot, q);
    mat3 R1 = fit_ring_procrustes(c_bot, _bot_ids, _bot_weights, _d_bot_rest, q);

    vec3 d_hat = c_bot - q0;
    const real L = std::max(d_hat.norm(), _cfg.L_min);
    d_hat /= L;

    apply_twist_to_top_frame(R0, R1, _rest.R0_rest, _rest.R1_rest, d_hat,
                             _cfg.twist_penalize_bend);

    for (size_t k = 0; k < _top_weights.size(); ++k) {
      const real wk = _top_weights[k];
      p.block(_id0_top + static_cast<index_t>(3 * k), 0, 3, 1) =
          _w * wk * (R0 * _d_top_rest[k]);
    }
    for (size_t k = 0; k < _bot_weights.size(); ++k) {
      const real wk = _bot_weights[k];
      p.block(_id0_bot + static_cast<index_t>(3 * k), 0, 3, 1) =
          _w * wk * (R1 * _d_bot_rest[k]);
    }
  }

private:
  static std::vector<index_t>
  top_ids_from_stencil(const slab_pierce_stencil &stencil) {
    return std::vector<index_t>(stencil.top_ids.begin(), stencil.top_ids.end());
  }

  static std::vector<real>
  weights_from_stencil(const slab_pierce_stencil &stencil) {
    return std::vector<real>(stencil.top_weights.begin(),
                             stencil.top_weights.end());
  }

  static std::vector<index_t>
  bot_ids_from_stencil(const slab_pierce_stencil &stencil) {
    return std::vector<index_t>(stencil.bot_ids.begin(), stencil.bot_ids.end());
  }

  static std::vector<real>
  bot_weights_from_stencil(const slab_pierce_stencil &stencil) {
    return std::vector<real>(stencil.bot_weights.begin(),
                             stencil.bot_weights.end());
  }

  slab_pierce_stencil _stencil;
  slab_prism_rest _rest;
  slab_prism_config _cfg;
  real _w_user = 0.0;
};

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // GAUDI_HEPWORTH_RELATIVE_PLANE_TWIST_HPP
