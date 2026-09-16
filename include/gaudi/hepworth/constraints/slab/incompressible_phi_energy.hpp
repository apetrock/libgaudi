#ifndef GAUDI_HEPWORTH_INCOMPRESSIBLE_PHI_ENERGY_HPP
#define GAUDI_HEPWORTH_INCOMPRESSIBLE_PHI_ENERGY_HPP

#include "implicit_field.hpp"
#include "phi_div_free_basis.hpp"
#include "slab_pierce_stencil.hpp"
#include "slab_prism_operators.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

/// Gauge sample ΔD = f_r(φ_q(x)) targeting zero (not f_r(x) − f_r(φ_q(x))).
/// @p foot_rest  Fixed foot origin used when Q_r was fitted (not the live hub).
inline field_eval eval_phi_gauge_sample(const albers::vec10 &Q, const vec3 &local,
                                        const phi_div_free_basis &basis,
                                        const vecX &q, index_t vi) {
  const vec3 u = basis.eval_u(local, vi, q);
  return eval_quadric_field(Q, local + u);
}

/// Ring spoke samples used for gauge energy (exclude hub verts — not in Q_r fit).
inline void for_each_phi_gauge_sample(
    const slab_pierce_stencil &stencil, const slab_ring_positions &ring,
    const std::function<void(index_t vi)> &fn) {
  for (index_t vi : ring.ids) {
    if (vi == stencil.foot || vi == stencil.pierce_foot)
      continue;
    fn(vi);
  }
}

/// Bulk density for φ gauge: (f_r(φ_q(x)))² — value only (grad is nonzero on the rest surface).
inline real phi_gauge_bulk_density(const field_eval &F, real lam_gap) {
  return lam_gap * F.f * F.f;
}

inline real eval_incompressible_phi_gauge_energy(
    const slab_pierce_stencil &stencil, const albers::vec10 &Q,
    const slab_ring_positions &ring, const phi_div_free_basis &basis,
    const vec3 &foot_rest, const vecX &q, real lam_gap, real mu = 1.0) {
  if (basis.dim_q <= 0)
    return 0.0;
  real E = 0.0;
  for_each_phi_gauge_sample(stencil, ring, [&](index_t vi) {
    const vec3 local = ring.at(vi) - foot_rest;
    const field_eval delta = eval_phi_gauge_sample(Q, local, basis, q, vi);
    if (!std::isfinite(delta.f))
      return;
    E += mu * phi_gauge_bulk_density(delta, lam_gap);
  });
  return E;
}

inline vecX eval_incompressible_phi_gauge_grad_q_analytic(
    const slab_pierce_stencil &stencil, const albers::vec10 &Q,
    const slab_ring_positions &ring, const phi_div_free_basis &basis,
    const vec3 &foot_rest, const vecX &q, real lam_gap, real mu = 1.0) {
  const int dim = basis.dim_q;
  vecX g = vecX::Zero(dim);
  if (dim <= 0)
    return g;

  for_each_phi_gauge_sample(stencil, ring, [&](index_t vi) {
    const vec3 local = ring.at(vi) - foot_rest;
    const vec3 u = basis.eval_u(local, vi, q);
    const field_eval F = eval_quadric_field(Q, local + u);
    if (!std::isfinite(F.f))
      return;
    for (int j = 0; j < dim; ++j) {
      const vec3 bj = basis.eval_basis_field(j, local, vi);
      const real dF = F.grad.dot(bj);
      g[j] += mu * lam_gap * 2.0 * F.f * dF;
    }
  });
  return g;
}

inline vecX eval_incompressible_phi_gauge_grad_q_fd(
    const slab_pierce_stencil &stencil, const albers::vec10 &Q,
    const slab_ring_positions &ring, const phi_div_free_basis &basis,
    const vec3 &foot_rest, const vecX &q, real lam_gap, real mu = 1.0,
    real eps = 1e-6) {
  const int dim = basis.dim_q;
  vecX g = vecX::Zero(dim);
  for (int i = 0; i < dim; ++i) {
    vecX qp = q;
    vecX qm = q;
    qp[i] += eps;
    qm[i] -= eps;
    const real ep = eval_incompressible_phi_gauge_energy(stencil, Q, ring, basis,
                                                       foot_rest, qp, lam_gap, mu);
    const real em = eval_incompressible_phi_gauge_energy(stencil, Q, ring, basis,
                                                       foot_rest, qm, lam_gap, mu);
    g[i] = (ep - em) / (2.0 * eps);
  }
  return g;
}

inline slab_ring_positions apply_phi_to_ring(
    const slab_pierce_stencil &stencil, const slab_ring_positions &ring,
    const phi_div_free_basis &basis, const vecX &q) {
  slab_ring_positions out = ring;
  for (index_t vi : ring.ids) {
    const vec3 local = ring.at(vi) - basis.foot_rest;
    const vec3 u = basis.eval_u(local, vi, q);
    if (!u.allFinite())
      continue;
    out.set(vi, ring.at(vi) + u);
  }
  return out;
}

inline vecX resolve_incompressible_phi_gd(
    const slab_pierce_stencil &stencil, const albers::vec10 &Q,
    const slab_ring_positions &ring, const phi_div_free_basis &basis,
    const vec3 &foot_rest, vecX q0, real lam_gap, int iters,
    real step_scale = 0.25, real mu = 1.0) {
  const int dim = basis.dim_q;
  vecX q = q0;
  if (q.size() < dim)
    q = vecX::Zero(dim);
  if (!q.allFinite())
    q = vecX::Zero(dim);
  if (iters <= 0 || dim <= 0)
    return q;

  const vecX q_zero = vecX::Zero(dim);
  const real e_rest = eval_incompressible_phi_gauge_energy(
      stencil, Q, ring, basis, foot_rest, q_zero, lam_gap, mu);
  if (e_rest < 1e-20)
    return q_zero;

  for (int iter = 0; iter < iters; ++iter) {
    const vecX grad = eval_incompressible_phi_gauge_grad_q_analytic(
        stencil, Q, ring, basis, foot_rest, q, lam_gap, mu);
    if (!grad.allFinite())
      break;

    real g2 = 0.0;
    for (int i = 0; i < grad.size(); ++i)
      g2 += grad[i] * grad[i];
    if (g2 < 1e-32)
      break;

    const real e0 = eval_incompressible_phi_gauge_energy(stencil, Q, ring, basis,
                                                         foot_rest, q, lam_gap, mu);
    if (!std::isfinite(e0))
      break;

    real alpha = step_scale / std::sqrt(g2);
    vecX trial = q;
    bool improved = false;
    for (int ls = 0; ls < 8; ++ls) {
      trial = q;
      for (int i = 0; i < trial.size(); ++i)
        trial[i] -= alpha * grad[i];
      if (!trial.allFinite())
        break;
      const real e1 = eval_incompressible_phi_gauge_energy(stencil, Q, ring, basis,
                                                         foot_rest, trial, lam_gap, mu);
      if (std::isfinite(e1) && e1 < e0) {
        q = trial;
        improved = true;
        break;
      }
      alpha *= 0.5;
    }
    if (!improved)
      break;
    if (e0 < 1e-24)
      break;
  }
  if (!q.allFinite())
    q = vecX::Zero(dim);
  return q;
}

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // GAUDI_HEPWORTH_INCOMPRESSIBLE_PHI_ENERGY_HPP
