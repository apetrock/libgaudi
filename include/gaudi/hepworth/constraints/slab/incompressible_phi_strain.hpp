#ifndef GAUDI_HEPWORTH_INCOMPRESSIBLE_PHI_STRAIN_HPP
#define GAUDI_HEPWORTH_INCOMPRESSIBLE_PHI_STRAIN_HPP

#include "incompressible_phi_energy.hpp"
#include "level_set_field_bundle.hpp"
#include "local_slab_field_fit.hpp"
#include "phi_div_free_basis.hpp"
#include "slab_dual_ring_arap_constraint.hpp"
#include "slab_prism_config.hpp"
#include "slab_prism_operators.hpp"

#include "gaudi/define_create_func.h"

namespace gaudi {
namespace hepworth {
namespace block {

/// Incompressible phi material: fixed quadric Q_r, live f = f_r ∘ phi_q (GD on q).
/// Local ROM grounds preferred ring x*; emit classic ARAP p = w R d_rest with
/// R = Procrustes(x*, d_rest) so neighboring feet transport via hub-spoke votes.
class incompressible_phi_strain : public slab_dual_ring_arap_constraint {
public:
  DEFINE_CREATE_FUNC(incompressible_phi_strain)

  incompressible_phi_strain(const slab_pierce_stencil &stencil,
                            const slab_prism_config &cfg, const real &w,
                            std::vector<sim_block::ptr> blocks,
                            real avg_edge_length)
      : slab_dual_ring_arap_constraint(
            top_ids_from_stencil(stencil), weights_from_stencil(stencil),
            bot_ids_from_stencil(stencil), bot_weights_from_stencil(stencil), w,
            blocks),
        _stencil(stencil), _cfg(cfg), _avg_edge_length(avg_edge_length) {}

  incompressible_phi_strain(const slab_pierce_stencil &stencil,
                            const slab_prism_rest &rest,
                            const slab_prism_config &cfg, const real &w,
                            std::vector<sim_block::ptr> blocks,
                            slab_field_coeffs Q_rest, real avg_edge_length)
      : slab_dual_ring_arap_constraint(
            top_ids_from_stencil(stencil), weights_from_stencil(stencil),
            rest.d_top_rest, bot_ids_from_stencil(stencil),
            bot_weights_from_stencil(stencil), rest.d_bot_rest, w, blocks),
        _stencil(stencil), _rest(rest), _cfg(cfg), _Q_rest(std::move(Q_rest)),
        _avg_edge_length(avg_edge_length) {
    ensure_basis_from_rest();
  }

  virtual std::string name() override { return "incompressible_phi_strain"; }

  void init_rest(const vecX &q, solver_context &ctx) override {
    (void)ctx;
    if (_rest_initialized)
      return;
    const auto pos_at = [&](index_t vi) { return _blocks[0]->get_vec3(vi, q); };
    init_spoke_rest(pos_at);
    _rest = build_slab_prism_rest(_stencil, pos_at);
    local_field_fit_config fit_cfg;
    fit_cfg.foot_normal_weight = _cfg.field_foot_normal_w;
    _Q_rest = fit_phi_quadric_on_stencil(_stencil, pos_at, fit_cfg);
    _basis.build(_stencil, pos_at);
    _q_phi = vecX::Zero(_basis.dim_q);
    _rest_initialized = true;
  }

  virtual void project(const vecX &q, vecX &p) override {
    ensure_basis_from_rest();
    const auto pos_at = [&](index_t vi) -> vec3 {
      return _blocks[0]->get_vec3(vi, q);
    };
    const slab_ring_positions ring =
        slab_ring_positions::from_stencil(_stencil, pos_at);

    if (_basis.dim_q > 0 && _cfg.phi_nr_iters > 0) {
      _q_phi = resolve_incompressible_phi_gd(
          _stencil, _Q_rest.quadric, ring, _basis, _basis.foot_rest, _q_phi,
          _cfg.phi_lam_gap, _cfg.phi_nr_iters, 0.25, _cfg.phi_mu);
    }

    if (_basis.dim_q <= 0 || !_q_phi.allFinite())
      return;

    const slab_ring_positions x_star =
        apply_phi_to_ring(_stencil, ring, _basis, _q_phi);

    const auto x_at = [&](index_t vi) -> vec3 {
      const vec3 t = x_star.at(vi);
      return t.allFinite() ? t : ring.at(vi);
    };
    write_top_arap_from_resolved(x_at, p);
    write_bot_arap_from_resolved(x_at, p);
  }

private:
  static std::vector<index_t>
  top_ids_from_stencil(const slab_pierce_stencil &stencil) {
    return std::vector<index_t>(stencil.top_ids.begin(), stencil.top_ids.end());
  }

  static std::vector<real> weights_from_stencil(const slab_pierce_stencil &stencil) {
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

  void ensure_basis_from_rest() {
    if (_basis.dim_q > 0)
      return;
    if (!_rest_initialized)
      return;
    const auto pos_at = [&](index_t vi) -> vec3 {
      if (vi == _stencil.foot)
        return _rest.c_top_rest;
      if (vi == _stencil.pierce_foot)
        return _rest.c_bot_rest;
      for (size_t k = 0; k < _stencil.top_ids.size() - 1; ++k) {
        if (_stencil.top_ids[k + 1] == vi)
          return _rest.c_top_rest + _rest.d_top_rest[k];
      }
      for (size_t k = 0; k < _stencil.bot_ids.size() - 1; ++k) {
        if (_stencil.bot_ids[k + 1] == vi)
          return _rest.c_bot_rest + _rest.d_bot_rest[k];
      }
      return vec3::Zero();
    };
    _basis.build(_stencil, pos_at);
    if (_q_phi.size() != _basis.dim_q)
      _q_phi = vecX::Zero(_basis.dim_q);
  }

  slab_pierce_stencil _stencil;
  slab_prism_rest _rest;
  slab_prism_config _cfg;
  slab_field_coeffs _Q_rest;
  phi_div_free_basis _basis;
  real _avg_edge_length = 1.0;
  vecX _q_phi;
};

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // GAUDI_HEPWORTH_INCOMPRESSIBLE_PHI_STRAIN_HPP
