#ifndef GAUDI_HEPWORTH_SLAB_TOP_RING_ARAP_CONSTRAINT_HPP
#define GAUDI_HEPWORTH_SLAB_TOP_RING_ARAP_CONSTRAINT_HPP

#include "gaudi/hepworth/block/block_constraint.hpp"
#include "slab_pierce_stencil.hpp"

#include <Eigen/SVD>

#include <cassert>
#include <functional>
#include <vector>

namespace gaudi {
namespace hepworth {
namespace block {

/// Shared top-ring hub-and-spokes row layout (matches shell arap::fill_A).
class slab_top_ring_arap_constraint : public block_constraint {
public:
  slab_top_ring_arap_constraint(const std::vector<index_t> &top_ids,
                                const std::vector<real> &edge_weights,
                                const real &w, std::vector<sim_block::ptr> blocks)
      : block_constraint(top_ids, w, blocks), _edge_weights(edge_weights) {
    assert(top_ids.size() >= 2);
    assert(edge_weights.size() + 1 == top_ids.size());
    _d_rest.assign(edge_weights.size(), vec3::Zero());
  }

  slab_top_ring_arap_constraint(const std::vector<index_t> &top_ids,
                                const std::vector<real> &edge_weights,
                                const std::vector<vec3> &d_rest, const real &w,
                                std::vector<sim_block::ptr> blocks)
      : block_constraint(top_ids, w, blocks), _edge_weights(edge_weights),
        _d_rest(d_rest) {
    assert(top_ids.size() >= 2);
    assert(edge_weights.size() + 1 == top_ids.size());
    assert(d_rest.size() == edge_weights.size());
    _rest_initialized = true;
  }

  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) override {
    _id0 = id0;
    assert(!_ids.empty());
    assert(_edge_weights.size() + 1 == _ids.size());
    const index_t ii = _blocks[0]->get_offset_idx(_ids[0]);
    assert(ii >= 0);
    for (size_t k = 0; k < _edge_weights.size(); ++k) {
      const real wk = _edge_weights[k];
      const index_t jj = _blocks[0]->get_offset_idx(_ids[k + 1]);
      assert(jj >= 0);
      const index_t row = id0;
      for (int ax = 0; ax < 3; ax++)
        triplets.push_back(trip(row + ax, ii + ax, -_w * wk));
      for (int ax = 0; ax < 3; ax++)
        triplets.push_back(trip(row + ax, jj + ax, _w * wk));
      id0 += 3;
    }
  }

  size_t top_spoke_count() const { return _edge_weights.size(); }

protected:
  void init_spoke_rest(const std::function<vec3(index_t)> &pos_at) {
    const vec3 c = pos_at(_ids[0]);
    for (size_t k = 0; k < _edge_weights.size(); ++k)
      _d_rest[k] = pos_at(_ids[k + 1]) - c;
  }

  void write_spoke_targets(const vec3 &x_hub,
                           const std::function<vec3(index_t)> &target_at,
                           vecX &p) {
    assert(_id0 + static_cast<index_t>(3 * _edge_weights.size()) <=
           static_cast<index_t>(p.size()));
    for (size_t k = 0; k < _edge_weights.size(); ++k) {
      const real wk = _edge_weights[k];
      const vec3 target = _w * wk * (target_at(_ids[k + 1]) - x_hub);
      p.block(_id0 + static_cast<index_t>(3 * k), 0, 3, 1) = target;
    }
  }

  mat3 fit_top_procrustes(const vec3 &q0, const vecX &q) const {
    mat3 M = mat3::Zero();
    for (size_t k = 0; k < _edge_weights.size(); ++k) {
      const real wk = _edge_weights[k];
      const vec3 dq = _blocks[0]->get_vec3(_ids[k + 1], q) - q0;
      M += wk * dq * _d_rest[k].transpose();
    }
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

  std::vector<real> _edge_weights;
  std::vector<vec3> _d_rest;
};

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // GAUDI_HEPWORTH_SLAB_TOP_RING_ARAP_CONSTRAINT_HPP
