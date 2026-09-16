#ifndef GAUDI_HEPWORTH_SLAB_DUAL_RING_ARAP_CONSTRAINT_HPP
#define GAUDI_HEPWORTH_SLAB_DUAL_RING_ARAP_CONSTRAINT_HPP

#include "gaudi/hepworth/block/block_constraint.hpp"
#include "slab_pierce_stencil.hpp"

#include <Eigen/SVD>

#include <functional>
#include <vector>

namespace gaudi {
namespace hepworth {
namespace block {

/// Top + bottom hub-and-spoke rows (two ARAP rings per slab stencil).
class slab_dual_ring_arap_constraint : public block_constraint {
public:
  slab_dual_ring_arap_constraint(const std::vector<index_t> &top_ids,
                                 const std::vector<real> &top_weights,
                                 const std::vector<index_t> &bot_ids,
                                 const std::vector<real> &bot_weights,
                                 const real &w, std::vector<sim_block::ptr> blocks)
      : block_constraint(top_ids, w, blocks), _top_ids(top_ids),
        _top_weights(top_weights), _bot_ids(bot_ids), _bot_weights(bot_weights) {
    assert(top_ids.size() >= 2);
    assert(top_weights.size() + 1 == top_ids.size());
    assert(bot_ids.size() >= 2);
    assert(bot_weights.size() + 1 == bot_ids.size());
    _d_top_rest.resize(top_weights.size());
    _d_bot_rest.resize(bot_weights.size());
  }

  slab_dual_ring_arap_constraint(const std::vector<index_t> &top_ids,
                                 const std::vector<real> &top_weights,
                                 const std::vector<vec3> &d_top_rest,
                                 const std::vector<index_t> &bot_ids,
                                 const std::vector<real> &bot_weights,
                                 const std::vector<vec3> &d_bot_rest,
                                 const real &w, std::vector<sim_block::ptr> blocks)
      : block_constraint(top_ids, w, blocks), _top_ids(top_ids),
        _top_weights(top_weights), _d_top_rest(d_top_rest), _bot_ids(bot_ids),
        _bot_weights(bot_weights), _d_bot_rest(d_bot_rest) {
    assert(top_ids.size() >= 2);
    assert(top_weights.size() + 1 == top_ids.size());
    assert(d_top_rest.size() == top_weights.size());
    assert(bot_ids.size() >= 2);
    assert(bot_weights.size() + 1 == bot_ids.size());
    assert(d_bot_rest.size() == bot_weights.size());
    _rest_initialized = true;
  }

  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) override {
    _id0_top = id0;
    fill_ring_A(_top_ids, _top_weights, _id0_top, triplets);
    id0 += static_cast<index_t>(3 * _top_weights.size());

    _id0_bot = id0;
    fill_ring_A(_bot_ids, _bot_weights, _id0_bot, triplets);
    id0 += static_cast<index_t>(3 * _bot_weights.size());

    _id0 = _id0_top;
  }

  size_t top_spoke_count() const { return _top_weights.size(); }
  size_t bot_spoke_count() const { return _bot_weights.size(); }
  size_t total_spoke_count() const {
    return top_spoke_count() + bot_spoke_count();
  }

protected:
  void init_spoke_rest(const std::function<vec3(index_t)> &pos_at) {
    const vec3 c_top = pos_at(_top_ids[0]);
    for (size_t k = 0; k < _top_weights.size(); ++k)
      _d_top_rest[k] = pos_at(_top_ids[k + 1]) - c_top;
    const vec3 c_bot = pos_at(_bot_ids[0]);
    for (size_t k = 0; k < _bot_weights.size(); ++k)
      _d_bot_rest[k] = pos_at(_bot_ids[k + 1]) - c_bot;
  }

  void write_top_spoke_targets(const vec3 &x_hub,
                               const std::function<vec3(index_t)> &target_at,
                               vecX &p) const {
    write_ring_spoke_targets(_id0_top, _top_ids, _top_weights, x_hub, target_at,
                             p);
  }

  void write_bot_spoke_targets(const vec3 &x_hub,
                               const std::function<vec3(index_t)> &target_at,
                               vecX &p) const {
    write_ring_spoke_targets(_id0_bot, _bot_ids, _bot_weights, x_hub, target_at,
                             p);
  }

  /// Fit R from a resolved ring configuration x*, then emit classic ARAP
  /// targets p = w R d_rest (not raw x*-hub). Grounds slab/ROM energy in the
  /// same hub-spoke transport as shell ARAP.
  void write_top_arap_from_resolved(const std::function<vec3(index_t)> &x_at,
                                    vecX &p) const {
    write_ring_arap_from_resolved(_id0_top, _top_ids, _top_weights, _d_top_rest,
                                  x_at, p);
  }

  void write_bot_arap_from_resolved(const std::function<vec3(index_t)> &x_at,
                                    vecX &p) const {
    write_ring_arap_from_resolved(_id0_bot, _bot_ids, _bot_weights, _d_bot_rest,
                                  x_at, p);
  }

  mat3 fit_ring_procrustes(const vec3 &hub, const std::vector<index_t> &ring_ids,
                           const std::vector<real> &weights,
                           const std::vector<vec3> &d_rest,
                           const vecX &q) const {
    mat3 M = mat3::Zero();
    for (size_t k = 0; k < weights.size(); ++k) {
      const real wk = weights[k];
      const vec3 dq = _blocks[0]->get_vec3(ring_ids[k + 1], q) - hub;
      M += wk * dq * d_rest[k].transpose();
    }
    return procrustes_R_from_M(M);
  }

  index_t _id0_top = 0;
  index_t _id0_bot = 0;
  std::vector<index_t> _top_ids;
  std::vector<real> _top_weights;
  std::vector<vec3> _d_top_rest;
  std::vector<index_t> _bot_ids;
  std::vector<real> _bot_weights;
  std::vector<vec3> _d_bot_rest;

private:
  static mat3 procrustes_R_from_M(const mat3 &M) {
    if (M.squaredNorm() < 1e-24)
      return mat3::Identity();
    Eigen::JacobiSVD<mat3> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
    mat3 U = svd.matrixU();
    mat3 V = svd.matrixV();
    mat3 R = U * V.transpose();
    if (R.determinant() < real(0.0)) {
      U.col(2) *= real(-1.0);
      R = U * V.transpose();
    }
    return R;
  }

  void fill_ring_A(const std::vector<index_t> &ring_ids,
                   const std::vector<real> &weights, index_t id0,
                   std::vector<trip> &triplets) const {
    const index_t ii = _blocks[0]->get_offset_idx(ring_ids[0]);
    for (size_t k = 0; k < weights.size(); ++k) {
      const real wk = weights[k];
      const index_t jj = _blocks[0]->get_offset_idx(ring_ids[k + 1]);
      const index_t row = id0 + static_cast<index_t>(3 * k);
      for (int ax = 0; ax < 3; ax++)
        triplets.push_back(trip(row + ax, ii + ax, -_w * wk));
      for (int ax = 0; ax < 3; ax++)
        triplets.push_back(trip(row + ax, jj + ax, _w * wk));
    }
  }

  void write_ring_spoke_targets(index_t row0, const std::vector<index_t> &ring_ids,
                                const std::vector<real> &weights,
                                const vec3 &x_hub,
                                const std::function<vec3(index_t)> &target_at,
                                vecX &p) const {
    for (size_t k = 0; k < weights.size(); ++k) {
      const real wk = weights[k];
      const vec3 target =
          _w * wk * (target_at(ring_ids[k + 1]) - x_hub);
      p.block(row0 + static_cast<index_t>(3 * k), 0, 3, 1) = target;
    }
  }

  void write_ring_arap_targets(index_t row0, const std::vector<real> &weights,
                               const std::vector<vec3> &d_rest, const mat3 &R,
                               vecX &p) const {
    for (size_t k = 0; k < weights.size(); ++k) {
      const real wk = weights[k];
      const vec3 target = _w * wk * (R * d_rest[k]);
      p.block(row0 + static_cast<index_t>(3 * k), 0, 3, 1) = target;
    }
  }

  void write_ring_arap_from_resolved(
      index_t row0, const std::vector<index_t> &ring_ids,
      const std::vector<real> &weights, const std::vector<vec3> &d_rest,
      const std::function<vec3(index_t)> &x_at, vecX &p) const {
    if (weights.empty() || d_rest.size() != weights.size() ||
        ring_ids.size() < weights.size() + 1) {
      return;
    }
    const vec3 hub = x_at(ring_ids[0]);
    mat3 M = mat3::Zero();
    for (size_t k = 0; k < weights.size(); ++k) {
      const vec3 dq = x_at(ring_ids[k + 1]) - hub;
      if (!dq.allFinite())
        continue;
      M += weights[k] * dq * d_rest[k].transpose();
    }
    write_ring_arap_targets(row0, weights, d_rest, procrustes_R_from_M(M), p);
  }
};

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // GAUDI_HEPWORTH_SLAB_DUAL_RING_ARAP_CONSTRAINT_HPP
