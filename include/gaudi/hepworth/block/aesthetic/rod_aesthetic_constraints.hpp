#ifndef __HEP_ROD_AESTHETIC_CONSTRAINTS__
#define __HEP_ROD_AESTHETIC_CONSTRAINTS__

// Geometric / regularization priors that are not Cosserat stretch–shear /
// bend–twist physics. Kept out of rod_constraints.hpp so the physics surface
// stays clear.

#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include "gaudi/common.h"
#include "gaudi/hepworth/block/block_constraint.hpp"
#include "gaudi/hepworth/block/sim_block.hpp"
#include "gaudi/vec_addendum.h"

namespace gaudi {
namespace hepworth {
namespace block {

class uniform_step : public block_constraint {
public:
  typedef std::shared_ptr<uniform_step> ptr;
  static constexpr int kStencil = 5;

  static ptr create(const std::vector<index_t> &ids, const real &w,
                    std::vector<sim_block::ptr> blocks,
                    real induced_twist = 0.0) {
    return std::make_shared<uniform_step>(ids, w, blocks, induced_twist);
  }

  uniform_step(const std::vector<index_t> &ids, const real &w,
               std::vector<sim_block::ptr> blocks, real induced_twist = 0.0)
      : block_constraint(ids, w, blocks), _induced_twist(induced_twist) {
    assert(ids.size() >= kStencil);
  }

  virtual std::string name() { return typeid(*this).name(); }

  // Markley average: largest eigenvector of Σ q qᵀ (hemisphere-aligned).
  static quat markley_average(const quat us[kStencil]) {
    mat4 A = mat4::Zero();
    vec4 ref = us[0].coeffs();
    if (ref.squaredNorm() < 1.0e-32)
      return quat::Identity();
    for (int i = 0; i < kStencil; ++i) {
      vec4 c = us[i].coeffs();
      if (c.dot(ref) < 0.0)
        c = -c;
      A += c * c.transpose();
    }
    Eigen::SelfAdjointEigenSolver<mat4> es(A);
    if (es.info() != Eigen::Success)
      return us[0];
    quat mean(es.eigenvectors().col(3).data());
    mean.normalize();
    if (mean.coeffs().dot(ref) < 0.0)
      mean.coeffs() = -mean.coeffs();
    return mean;
  }

  virtual void project(const vecX &q, vecX &p) {
    quat u[kStencil];
    for (int i = 0; i < kStencil; ++i)
      u[i] = _blocks[0]->get_quat(this->_ids[i], q).normalized();

    quat mean = markley_average(u);
    if (std::abs(_induced_twist) > 1.0e-16) {
      mean = (mean * quat(Eigen::AngleAxis<real>(_induced_twist, vec3::UnitZ())))
                 .normalized();
    }

    const vec4 t(mean.coeffs().data());
    for (int i = 0; i < kStencil; ++i) {
      if (mean.coeffs().hasNaN()) {
        std::cout << __PRETTY_FUNCTION__ << " mean is nan" << std::endl;
        exit(0);
      }
      p.block(_id0 + 4 * i, 0, 4, 1) = _w * t;
    }
  }

  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    for (int i = 0; i < kStencil; ++i) {
      index_t ii = _blocks[0]->get_offset_idx(this->_ids[i]);
      for (int ax = 0; ax < 4; ++ax)
        triplets.push_back(trip(_id0 + 4 * i + ax, ii + ax, _w));
    }
    id0 += 4 * kStencil;
  }

  real _induced_twist = 0.0;
};
// Local tube prior from ONE local SVD of mean-centered stencil points.
//
// Citation:
//   J. Álvarez-Vizoso, R. Arn, M. Kirby, C. Peterson, B. Draper,
//   "Geometry of Curves in R^n from the Local Singular Value Decomposition",
//   Linear Algebra and its Applications 571 (2019) 180–210.
//   doi:10.1016/j.laa.2019.02.006
//   preprint: arXiv:1511.05008
//
// Method (Theorems 3.1–3.2): local covariance of curve samples → SVD;
//   singular vectors ≈ Frenet frame: u1≈T, u2≈N, u3≈B (descending σ);
//   κ₁ = (√20 / 3) · σ₂ / σ₁² ,   R = 1/κ₁.
//
// PD residual is differential (not a world-space point snap, not a smoother):
//   A encodes δ²p = p₋ + p₊ − 2 p₀  (same sparsity as smooth)
//   project target: δ² → h² κ N = (h²/R) N   (Frenet: γ'' = κ N)
// So this enforces osculating-cylinder curvature, not δ² → 0.
//
// ids layout (init_helicity): {ci, ip3, ip2, ip1, ip0, ci, in0, in1, in2, in3}
//   full stencil → SVD;  (ids[4], ids[0], ids[6]) = (p₋, p₀, p₊) → δ².
class helicitiy : public block_constraint {
public:
  typedef std::shared_ptr<helicitiy> ptr;

  static ptr create(const std::vector<index_t> &ids, const real &w,
                    std::vector<sim_block::ptr> blocks) {
    return std::make_shared<helicitiy>(ids, w, blocks);
  }

  helicitiy(const std::vector<index_t> &ids, const real &w,
            std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks) {}
  virtual std::string name() { return typeid(*this).name(); }

  struct tube_frame_t {
    vec3 tangent = vec3::UnitX();
    vec3 normal = vec3::UnitY(); // principal normal, toward center of curvature
    vec3 binormal = vec3::UnitZ();
    vec3 center = vec3::Zero();
    real kappa = 0.0;
    real radius = -1.0;
    bool ok = false;
  };

  // Ordered polyline: ip3..ip0, ci, in0..in3 (init_helicity layout).
  void gather_polyline(const vecX &q, vec3 &q0, std::vector<vec3> &xs) const {
    q0 = _blocks[0]->get_vec3(this->_ids[0], q);
    xs.clear();
    if (this->_ids.size() >= 10) {
      const index_t order[] = {1, 2, 3, 4, 0, 6, 7, 8, 9};
      xs.reserve(9);
      for (index_t k : order)
        xs.push_back(_blocks[0]->get_vec3(this->_ids[k], q));
      return;
    }
    xs.push_back(q0);
    for (size_t i = 1; i < this->_ids.size(); ++i) {
      if (this->_ids[i] == this->_ids[0])
        continue;
      xs.push_back(_blocks[0]->get_vec3(this->_ids[i], q));
    }
  }

  // Immediate neighbors for δ² (requires init_helicity 10-id layout).
  bool gather_triple(const vecX &q, vec3 &qm, vec3 &q0, vec3 &qp) const {
    if (this->_ids.size() < 7)
      return false;
    q0 = _blocks[0]->get_vec3(this->_ids[0], q);
    qm = _blocks[0]->get_vec3(this->_ids[4], q); // ip0
    qp = _blocks[0]->get_vec3(this->_ids[6], q); // in0
    return true;
  }

  // Single local SVD / eigh of mean-centered point covariance.
  static tube_frame_t tube_from_local_svd(const std::vector<vec3> &xs,
                                          const vec3 &q0) {
    tube_frame_t out;
    if (xs.size() < 3)
      return out;

    vec3 mean = vec3::Zero();
    for (const vec3 &x : xs)
      mean += x;
    mean /= real(xs.size());

    mat3 C = mat3::Zero();
    for (const vec3 &x : xs) {
      const vec3 d = x - mean;
      C += d * d.transpose();
    }
    C /= real(xs.size());

    Eigen::SelfAdjointEigenSolver<mat3> es(C);
    if (es.info() != Eigen::Success)
      return out;

    // Ascending eigenvalues λ0≤λ1≤λ2 → descending singular values σ_i=√λ.
    // Frenet: u1↔T (largest), u2↔N, u3↔B (smallest).
    const real l0 = std::max(real(0.0), es.eigenvalues()[0]);
    const real l1 = std::max(real(0.0), es.eigenvalues()[1]);
    const real l2 = std::max(real(0.0), es.eigenvalues()[2]);
    const real s1 = std::sqrt(l2); // σ1
    const real s2 = std::sqrt(l1); // σ2
    if (s1 < 1.0e-14 || s2 < 1.0e-14)
      return out;

    out.tangent = es.eigenvectors().col(2).normalized();
    out.normal = es.eigenvectors().col(1).normalized();
    out.binormal = es.eigenvectors().col(0).normalized();

    // Orient T along the polyline chord through the center vertex.
    if (xs.size() >= 3) {
      const size_t mid = xs.size() / 2;
      vec3 chord = xs[std::min(mid + 1, xs.size() - 1)] -
                   xs[mid > 0 ? mid - 1 : 0];
      if (chord.squaredNorm() > 1.0e-16 && chord.dot(out.tangent) < 0.0)
        out.tangent = -out.tangent;
    }

    // Right-handed Frenet: B = T × N, N = B × T after fixing.
    out.binormal = out.tangent.cross(out.normal);
    if (out.binormal.squaredNorm() < 1.0e-16)
      out.binormal = es.eigenvectors().col(0).normalized();
    else
      out.binormal.normalize();
    out.normal = out.binormal.cross(out.tangent).normalized();

    // κ₁ = (√20 / 3) σ₂ / σ₁²  (Álvarez et al. 2019, Thm. 3.2); R = 1/κ₁.
    constexpr real k_pre = 1.4907119849998597; // √20 / 3
    const real kappa = k_pre * s2 / (s1 * s1);
    if (kappa < 1.0e-14)
      return out;
    out.kappa = kappa;
    out.radius = 1.0 / kappa;

    // Principal normal toward interior (stencil mean ≈ toward axis).
    vec3 N = out.normal;
    const vec3 to_mean = mean - q0;
    const vec3 rad = to_mean - to_mean.dot(out.binormal) * out.binormal;
    if (rad.squaredNorm() > 1.0e-16 && rad.dot(N) < 0.0)
      N = -N;
    out.normal = N;
    out.center = q0 + out.radius * N; // C = γ + R N
    out.ok = true;
    return out;
  }

  virtual void project(const vecX &q, vecX &p) {
    vec3 qm, q0, qp;
    if (!gather_triple(q, qm, q0, qp)) {
      p.block(_id0, 0, 3, 1).setZero();
      return;
    }

    const vec3 d2 = qm + qp - real(2.0) * q0;

    std::vector<vec3> xs;
    vec3 q0b;
    gather_polyline(q, q0b, xs);
    const tube_frame_t tube = tube_from_local_svd(xs, q0);
    if (!tube.ok) {
      // Identity in δ²-space (no force) rather than a position snap.
      p.block(_id0, 0, 3, 1) = _w * d2;
      return;
    }

    // Align N with discrete acceleration (both should point to center).
    vec3 N = tube.normal;
    if (d2.squaredNorm() > 1.0e-16 && d2.dot(N) < 0.0)
      N = -N;

    const real h =
        real(0.5) * ((q0 - qm).norm() + (qp - q0).norm());
    if (h < 1.0e-14) {
      p.block(_id0, 0, 3, 1) = _w * d2;
      return;
    }

    // Frenet: γ'' = κ N  ⇒  δ²p → h² κ N.
    const vec3 target = (h * h * tube.kappa) * N;
    p.block(_id0, 0, 3, 1) = _w * target;
  }

  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    // δ² = p₋ + p₊ − 2 p₀  (ids: center=0, prev=ip0=4, next=in0=6)
    index_t i0 = _blocks[0]->get_offset_idx(this->_ids[0]);
    index_t im = _blocks[0]->get_offset_idx(this->_ids[4]);
    index_t ip = _blocks[0]->get_offset_idx(this->_ids[6]);

    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + ax, i0 + ax, -real(2.0) * _w));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + ax, im + ax, _w));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + ax, ip + ax, _w));
    id0 += 3;
  }
};

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // __HEP_ROD_AESTHETIC_CONSTRAINTS__
