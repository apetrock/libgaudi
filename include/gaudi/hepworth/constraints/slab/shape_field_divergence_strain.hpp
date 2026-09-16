#ifndef GAUDI_HEPWORTH_SHAPE_FIELD_DIVERGENCE_STRAIN_HPP
#define GAUDI_HEPWORTH_SHAPE_FIELD_DIVERGENCE_STRAIN_HPP

#include "level_set_field_bundle.hpp"
#include "slab_field_divergence_energy.hpp"
#include "slab_prism_config.hpp"
#include "slab_prism_operators.hpp"
#include "slab_top_ring_arap_constraint.hpp"

#include "gaudi/define_create_func.h"

#include <cassert>

namespace gaudi {
namespace hepworth {
namespace block {

/// Pierce stencil builds Q_r / Q_k; projection DOFs are the foot one-ring only.
///
/// Local step freezes Q_r and Q_k and places top-ring samples inside the frozen
/// ΔD gauge (bulk |∇ΔD|² + λ(ΔD)²). Live refit of Q_k is outer-loop only.
class shape_field_divergence_strain : public slab_top_ring_arap_constraint {
public:
  DEFINE_CREATE_FUNC(shape_field_divergence_strain)

  shape_field_divergence_strain(const slab_pierce_stencil &stencil,
                                const slab_prism_config &cfg, const real &w,
                                std::vector<sim_block::ptr> blocks,
                                real avg_edge_length)
      : slab_top_ring_arap_constraint(top_ids_from_stencil(stencil),
                                      weights_from_stencil(stencil), w, blocks),
        _stencil(stencil), _cfg(cfg), _avg_edge_length(avg_edge_length) {}

  shape_field_divergence_strain(
      const slab_pierce_stencil &stencil, const slab_prism_rest &rest,
      const slab_prism_config &cfg, const real &w,
      std::vector<sim_block::ptr> blocks, slab_field_coeffs Q_rest,
      real avg_edge_length)
      : slab_top_ring_arap_constraint(top_ids_from_stencil(stencil),
                                      weights_from_stencil(stencil),
                                      rest.d_top_rest, w, blocks),
        _stencil(stencil), _rest(rest), _cfg(cfg), _Q_rest(std::move(Q_rest)),
        _avg_edge_length(avg_edge_length) {
    _Q_live = _Q_rest;
    if (_cfg.field_live != field_live_update::live_matches_rest)
      _Q_live_valid = false;
  }

  virtual std::string name() override { return "shape_field_divergence_strain"; }

  void init_rest(const vecX &q, solver_context &ctx) override {
    (void)ctx;
    if (_rest_initialized)
      return;
    assert(!_ids.empty() && _ids[0] == _stencil.foot);
    const auto pos_at = [&](index_t vi) { return _blocks[0]->get_vec3(vi, q); };
    init_spoke_rest(pos_at);
    _rest = build_slab_prism_rest(_stencil, pos_at);
    local_field_fit_config fit_cfg;
    fit_cfg.foot_normal_weight = _cfg.field_foot_normal_w;
    fit_cfg.max_rbf_centers = slab_field_max_rbf_centers;
    // Fit uses the full pierce stencil (top + bottom rings).
    _Q_rest = fit_field_coeffs_on_stencil(_cfg.field_family, _stencil, pos_at,
                                          fit_cfg, _avg_edge_length);
    _Q_live = _Q_rest;
    _Q_live_valid = _cfg.field_live == field_live_update::live_matches_rest;
    _project_count = 0;
    _rest_initialized = true;
  }

  /// Refit Q_k from current pierce stencil. Q_r stays fixed.
  void refit_live_field(const vecX &q) {
    if (_cfg.field_live == field_live_update::live_matches_rest)
      return;
    const auto pos_at = [&](index_t vi) -> vec3 {
      return _blocks[0]->get_vec3(vi, q);
    };
    local_field_fit_config fit_cfg;
    fit_cfg.foot_normal_weight = _cfg.field_foot_normal_w;
    fit_cfg.max_rbf_centers = slab_field_max_rbf_centers;
    _Q_live = fit_field_coeffs_on_stencil(_Q_rest.family, _stencil, pos_at,
                                          fit_cfg, _avg_edge_length, &_Q_rest);
    _Q_live_valid = true;
  }

  virtual void project(const vecX &q, vecX &p) override {
    // Stale pierce ids after remesh pack → OOB q.block and heap corruption.
    const index_t nq = static_cast<index_t>(q.size());
    for (index_t vi : _ids) {
      const index_t off = _blocks[0]->get_offset_idx(vi);
      if (off < 0 || off + 2 >= nq)
        return;
    }
    for (index_t vi : _stencil.bot_ids) {
      const index_t off = _blocks[0]->get_offset_idx(vi);
      if (off < 0 || off + 2 >= nq)
        return;
    }

    const bool refit_live = [&]() {
      if (_cfg.field_live == field_live_update::live_matches_rest)
        return false;
      switch (_cfg.field_live) {
      case field_live_update::refit_every_project:
        return true;
      case field_live_update::refit_every_n:
        return !_Q_live_valid ||
               (_project_count % std::max(_cfg.field_refit_period, 1)) == 0;
      default:
        return false;
      }
    }();

    if (refit_live)
      refit_live_field(q);
    ++_project_count;

    const auto pos_at = [&](index_t vi) -> vec3 {
      return _blocks[0]->get_vec3(vi, q);
    };

    top_ring_positions top =
        top_ring_positions::from_stencil(_stencil, pos_at);

    const bool skip_gd =
        _cfg.field_live == field_live_update::live_matches_rest;

    top_ring_positions x_star = top;
    if (!skip_gd && _cfg.field_gd_iters > 0) {
      // Q frozen for this local placement step (Green / bulk ΔD gauge).
      const real e0 = eval_shape_field_divergence_energy(
          _stencil, _Q_rest, _Q_live, top, _cfg.lam_gap);
      if (std::isfinite(e0) && e0 > 0.0) {
        x_star = resolve_shape_field_divergence_gd(
            _stencil, _Q_rest, _Q_live, top, _cfg.lam_gap, _cfg.field_gd_iters,
            _cfg.field_gd_step_scale);

        for (const vec3 &xi : x_star.x) {
          if (!xi.allFinite()) {
            x_star = top;
            break;
          }
        }
      }
    }

    write_spoke_targets(x_star.at(_stencil.foot),
                        [&](index_t vi) { return x_star.at(vi); }, p);
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

  slab_pierce_stencil _stencil;
  slab_prism_rest _rest;
  slab_prism_config _cfg;
  slab_field_coeffs _Q_rest;
  slab_field_coeffs _Q_live;
  real _avg_edge_length = 1.0;
  int _project_count = 0;
  bool _Q_live_valid = false;
};

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // GAUDI_HEPWORTH_SHAPE_FIELD_DIVERGENCE_STRAIN_HPP
