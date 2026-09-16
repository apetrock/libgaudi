#ifndef GAUDI_HEPWORTH_SLAB_DUAL_RING_VOLUME_PROCRUSTES_HPP
#define GAUDI_HEPWORTH_SLAB_DUAL_RING_VOLUME_PROCRUSTES_HPP

#include "gaudi/hepworth/block/block_constraint.hpp"
#include "slab_foot_projection_probe.hpp"
#include "slab_prism_config.hpp"
#include "gaudi/hepworth/block/slab_prism_log.hpp"
#include "slab_prism_operators.hpp"

#include "gaudi/define_create_func.h"

#include <vector>

namespace gaudi {
namespace hepworth {
namespace block {

/// Top-ring bulk: σ/volume springs in project(); pump along the live bulk SVD axis
/// that best aligns with the one-ring outward normal.
class slab_dual_ring_volume_procrustes : public block_constraint {
public:
  DEFINE_CREATE_FUNC(slab_dual_ring_volume_procrustes)

  slab_dual_ring_volume_procrustes(const slab_pierce_stencil &stencil,
                                   const slab_prism_config &cfg, const real &w,
                                   std::vector<sim_block::ptr> blocks)
      : block_constraint(top_ids_from_stencil(stencil), w, blocks),
        _stencil(stencil), _cfg(cfg), _w_user(w) {}

  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) override {
    _id0 = id0;
    const index_t i_foot = _blocks[0]->get_offset_idx(_stencil.foot);
    for (size_t k = 0; k < _top_ring_samples.size(); ++k) {
      const real wk = _w * _top_ring_samples[k].weight;
      const index_t jj =
          _blocks[0]->get_offset_idx(_top_ring_samples[k].vi);
      const index_t row = id0 + static_cast<index_t>(3 * k);
      for (int ax = 0; ax < 3; ++ax) {
        triplets.push_back(trip(row + ax, i_foot + ax, -wk));
        triplets.push_back(trip(row + ax, jj + ax, wk));
      }
    }
    id0 += static_cast<index_t>(3 * _top_ring_samples.size());
  }

  virtual std::string name() override {
    return "slab_dual_ring_volume_procrustes";
  }

  void init_rest(const vecX &q, solver_context &ctx) override {
    (void)ctx;
    if (_rest_initialized)
      return;
    const auto pos_at = [&](index_t vi) { return _blocks[0]->get_vec3(vi, q); };
    vec3 foot;
    _top_ring_samples = build_top_ring_fit_samples(_stencil, pos_at, foot);
    _d_rest = pos_at(_stencil.pierce_foot) - foot;
    _bulk_rest = build_foot_frame_bulk_rest(_top_ring_samples, _cfg.L_min);
    _L_spine = std::max(_d_rest.norm(), _cfg.L_min);
    _w = effective_slab_weight(_w_user, _L_spine, _cfg.L_min, _cfg.w_eff_max);
    const std::vector<index_t> top_ids(_stencil.top_ids.begin(),
                                         _stencil.top_ids.end());
    _n_rest = outward_ring_normal(
        top_ring_pole_normal(foot, top_ids, pos_at, _cfg.L_min), _d_rest,
        _cfg.L_min);
    _rest_initialized = true;
  }

  virtual void project(const vecX &q, vecX &p) override {
    const auto pos_at = [&](index_t vi) { return _blocks[0]->get_vec3(vi, q); };
    const vec3 hub = pos_at(_stencil.foot);
    const std::vector<index_t> top_ids(_stencil.top_ids.begin(),
                                         _stencil.top_ids.end());

    const mat3 M_cross =
        foot_frame_cross_covariance(_stencil.foot, _top_ring_samples, pos_at);
    const mat3 Cov_live =
        foot_frame_live_covariance(_stencil.foot, _top_ring_samples, pos_at);
    const bulk_ring_map m = bulk_procrustes_map(
        M_cross, Cov_live, _bulk_rest, _cfg.bulk_max_stretch, _cfg.bulk_w_volume,
        _cfg.bulk_w_length);

    const vec3 d_spine = pos_at(_stencil.pierce_foot) - hub;
    vec3 n_ring = _n_rest;
    if (_cfg.bulk_pump_use_live_normal) {
      n_ring = outward_ring_normal(
          top_ring_pole_normal(hub, top_ids, pos_at, _cfg.L_min), d_spine,
          _cfg.L_min);
    }
    // Live principal axes from Cov_live (not cross-covariance).
    const vec3 pump_axis =
        bulk_live_axis_nearest_normal(Cov_live, n_ring, _cfg.L_min);

    const real V_live = m.l0[0] * m.l0[1] * m.l0[2];
    const real V_target = m.l[0] * m.l[1] * m.l[2];
    const real max_pump =
        _cfg.bulk_target_step_frac * std::max(_L_spine, _cfg.L_min);
    const real delta = bulk_volume_pump_delta(
        V_live, V_target, _bulk_rest.sigma_rest, _bulk_rest.V_rest, _L_spine,
        _cfg.L_min, max_pump);

    if (slab_foot_probe_active(_stencil.foot)) {
      slab_prism_log(
          "bulk foot=" + std::to_string(_stencil.foot) +
          " |l0-lt|=" +
          std::to_string((m.l0 - _bulk_rest.sigma_rest).norm()) +
          " Z~I=" + std::to_string((m.Z - mat3::Identity()).norm()) +
          " dV=" + std::to_string(V_target - V_live) +
          " pump=" + std::to_string(delta) +
          " axis·n=" + std::to_string(pump_axis.dot(n_ring)));
    }

    const real max_spoke_step = max_pump;
    for (size_t k = 0; k < _top_ring_samples.size(); ++k) {
      const slab_foot_local_sample &s = _top_ring_samples[k];
      const real wk = _w * s.weight;
      // ARAP-style transport: Z maps rest spokes (R from live, stretch from
      // Cov_live vs σ_rest); pump is the volume correction along the normal.
      vec3 tgt = m.Z * s.local_rest + delta * pump_axis;
      const vec3 cur = pos_at(s.vi) - hub;
      tgt = clamp_slab_spoke_target_delta(tgt, cur, max_spoke_step);
      p.block(_id0 + static_cast<index_t>(3 * k), 0, 3, 1) = wk * tgt;
    }
  }

  size_t fit_sample_count() const { return _top_ring_samples.size(); }

private:
  static std::vector<index_t>
  top_ids_from_stencil(const slab_pierce_stencil &stencil) {
    return std::vector<index_t>(stencil.top_ids.begin(), stencil.top_ids.end());
  }

  slab_pierce_stencil _stencil;
  slab_prism_config _cfg;
  real _w_user = 0.0;
  std::vector<slab_foot_local_sample> _top_ring_samples;
  foot_frame_bulk_rest _bulk_rest;
  vec3 _d_rest = vec3::Zero();
  vec3 _n_rest = vec3::UnitZ();
  real _L_spine = 1.0;
  index_t _id0 = 0;
};

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // GAUDI_HEPWORTH_SLAB_DUAL_RING_VOLUME_PROCRUSTES_HPP
