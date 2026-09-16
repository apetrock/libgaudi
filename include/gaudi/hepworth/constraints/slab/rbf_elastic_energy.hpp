#ifndef GAUDI_HEPWORTH_RBF_ELASTIC_ENERGY_HPP
#define GAUDI_HEPWORTH_RBF_ELASTIC_ENERGY_HPP

#include "implicit_field.hpp"
#include "rbf_displacement_basis.hpp"
#include "slab_pierce_stencil.hpp"
#include "slab_prism_operators.hpp"

#include <cmath>
#include <functional>

namespace gaudi {
namespace hepworth {
namespace block {

inline void for_each_rbf_gauge_sample(const slab_pierce_stencil &stencil,
                                      const slab_ring_positions &ring,
                                      const std::function<void(index_t)> &fn) {
  for (index_t vi : ring.ids) {
    if (vi == stencil.foot || vi == stencil.pierce_foot)
      continue;
    fn(vi);
  }
}

inline slab_ring_positions apply_rbf_to_ring(const slab_ring_positions &ring,
                                             const rbf_displacement_basis &basis,
                                             const vecX &q) {
  slab_ring_positions out = ring;
  for (index_t vi : ring.ids) {
    const vec3 u = basis.eval_u_at_vertex(vi, q);
    if (!u.allFinite())
      continue;
    out.set(vi, ring.at(vi) + u);
  }
  return out;
}

/// Rest-surface gauge: Σ (f_r(x_i + u_i))² on ring spokes (fixed Q_r, foot_rest).
inline real eval_rbf_gauge_energy(const slab_pierce_stencil &stencil,
                                  const albers::vec10 &Q,
                                  const slab_ring_positions &ring,
                                  const rbf_displacement_basis &basis,
                                  const vec3 &foot_rest, const vecX &q,
                                  real lam_gap, real mu) {
  if (basis.dim_q <= 0 || !(mu > 0.0))
    return 0.0;
  real E = 0.0;
  for_each_rbf_gauge_sample(stencil, ring, [&](index_t vi) {
    const vec3 local = ring.at(vi) - foot_rest;
    const vec3 u = basis.eval_u_at_vertex(vi, q);
    const field_eval F = eval_quadric_field(Q, local + u);
    if (!std::isfinite(F.f))
      return;
    E += mu * lam_gap * F.f * F.f;
  });
  return E;
}

inline vecX eval_rbf_gauge_grad_q(const slab_pierce_stencil &stencil,
                                  const albers::vec10 &Q,
                                  const slab_ring_positions &ring,
                                  const rbf_displacement_basis &basis,
                                  const vec3 &foot_rest, const vecX &q,
                                  real lam_gap, real mu) {
  const int dim = basis.dim_q;
  vecX g = vecX::Zero(dim);
  if (dim <= 0 || !(mu > 0.0))
    return g;

  for_each_rbf_gauge_sample(stencil, ring, [&](index_t vi) {
    const vec3 local = ring.at(vi) - foot_rest;
    const vec3 u = basis.eval_u_at_vertex(vi, q);
    const field_eval F = eval_quadric_field(Q, local + u);
    if (!std::isfinite(F.f))
      return;
    for (int j = 0; j < dim; ++j) {
      const vec3 bj = basis.eval_basis_field_at_vertex(j, vi);
      g[j] += mu * lam_gap * 2.0 * F.f * F.grad.dot(bj);
    }
  });
  return g;
}

/// Spoke ARAP strain of reconstructed ring vs rest spokes (Procrustes per ring).
inline real eval_rbf_strain_energy(const slab_pierce_stencil &stencil,
                                   const slab_prism_rest &rest,
                                   const slab_ring_positions &x_star, real eta) {
  if (!(eta > 0.0))
    return 0.0;

  const auto ring_energy = [&](index_t hub, const auto &ids,
                               const std::vector<vec3> &d_rest,
                               const std::vector<real> &wts) -> real {
    if (d_rest.empty() || wts.size() != d_rest.size())
      return 0.0;
    if (ids.size() < d_rest.size() + 1)
      return 0.0;
    const vec3 c = x_star.at(hub);
    mat3 M = mat3::Zero();
    for (size_t k = 0; k < d_rest.size(); ++k) {
      const index_t vi = ids[k + 1];
      const vec3 d = x_star.at(vi) - c;
      M += wts[k] * d * d_rest[k].transpose();
    }
    const mat3 R = (M.squaredNorm() < 1e-24) ? mat3::Identity()
                                             : fit_procrustes_rotation(M);
    real E = 0.0;
    for (size_t k = 0; k < d_rest.size(); ++k) {
      const index_t vi = ids[k + 1];
      const vec3 d = x_star.at(vi) - c;
      const vec3 r = d - R * d_rest[k];
      E += wts[k] * r.squaredNorm();
    }
    return E;
  };

  real E = 0.0;
  E += ring_energy(stencil.foot, stencil.top_ids, rest.d_top_rest, rest.top_w);
  E += ring_energy(stencil.pierce_foot, stencil.bot_ids, rest.d_bot_rest,
                   rest.bot_w);
  return eta * E;
}

inline vecX eval_rbf_strain_grad_q_fd(const slab_pierce_stencil &stencil,
                                      const slab_prism_rest &rest,
                                      const slab_ring_positions &ring,
                                      const rbf_displacement_basis &basis,
                                      const vecX &q, real eta,
                                      real eps = 1e-6) {
  const int dim = basis.dim_q;
  vecX g = vecX::Zero(dim);
  if (dim <= 0 || !(eta > 0.0))
    return g;
  for (int i = 0; i < dim; ++i) {
    vecX qp = q;
    vecX qm = q;
    qp[i] += eps;
    qm[i] -= eps;
    const real ep = eval_rbf_strain_energy(
        stencil, rest, apply_rbf_to_ring(ring, basis, qp), eta);
    const real em = eval_rbf_strain_energy(
        stencil, rest, apply_rbf_to_ring(ring, basis, qm), eta);
    g[i] = (ep - em) / (2.0 * eps);
  }
  return g;
}

inline real eval_rbf_elastic_energy(const slab_pierce_stencil &stencil,
                                    const albers::vec10 &Q,
                                    const slab_prism_rest &rest,
                                    const slab_ring_positions &ring,
                                    const rbf_displacement_basis &basis,
                                    const vec3 &foot_rest, const vecX &q,
                                    real lam_gap, real mu, real eta) {
  const real Eg =
      eval_rbf_gauge_energy(stencil, Q, ring, basis, foot_rest, q, lam_gap, mu);
  const slab_ring_positions x_star = apply_rbf_to_ring(ring, basis, q);
  const real Es = eval_rbf_strain_energy(stencil, rest, x_star, eta);
  return Eg + Es;
}

inline vecX eval_rbf_elastic_grad_q(const slab_pierce_stencil &stencil,
                                    const albers::vec10 &Q,
                                    const slab_prism_rest &rest,
                                    const slab_ring_positions &ring,
                                    const rbf_displacement_basis &basis,
                                    const vec3 &foot_rest, const vecX &q,
                                    real lam_gap, real mu, real eta) {
  vecX g = eval_rbf_gauge_grad_q(stencil, Q, ring, basis, foot_rest, q, lam_gap,
                                 mu);
  g += eval_rbf_strain_grad_q_fd(stencil, rest, ring, basis, q, eta);
  return g;
}

inline vecX resolve_rbf_elastic_gd(const slab_pierce_stencil &stencil,
                                   const albers::vec10 &Q,
                                   const slab_prism_rest &rest,
                                   const slab_ring_positions &ring,
                                   const rbf_displacement_basis &basis,
                                   const vec3 &foot_rest, vecX q0, real lam_gap,
                                   real mu, real eta, int iters,
                                   real step_scale = 0.25) {
  const int dim = basis.dim_q;
  vecX q = q0;
  if (q.size() < dim)
    q = vecX::Zero(dim);
  if (!q.allFinite())
    q = vecX::Zero(dim);
  if (iters <= 0 || dim <= 0)
    return q;

  for (int iter = 0; iter < iters; ++iter) {
    const vecX grad = eval_rbf_elastic_grad_q(stencil, Q, rest, ring, basis,
                                              foot_rest, q, lam_gap, mu, eta);
    if (!grad.allFinite())
      break;
    real g2 = 0.0;
    for (int i = 0; i < grad.size(); ++i)
      g2 += grad[i] * grad[i];
    if (g2 < 1e-32)
      break;

    const real e0 = eval_rbf_elastic_energy(stencil, Q, rest, ring, basis,
                                            foot_rest, q, lam_gap, mu, eta);
    if (!std::isfinite(e0))
      break;

    real alpha = step_scale / std::sqrt(g2);
    bool improved = false;
    for (int ls = 0; ls < 8; ++ls) {
      vecX trial = q;
      for (int i = 0; i < trial.size(); ++i)
        trial[i] -= alpha * grad[i];
      if (!trial.allFinite())
        break;
      const real e1 = eval_rbf_elastic_energy(stencil, Q, rest, ring, basis,
                                              foot_rest, trial, lam_gap, mu, eta);
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

#endif // GAUDI_HEPWORTH_RBF_ELASTIC_ENERGY_HPP
