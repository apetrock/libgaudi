#ifndef GAUDI_HEPWORTH_SLAB_FIELD_DIVERGENCE_ENERGY_HPP
#define GAUDI_HEPWORTH_SLAB_FIELD_DIVERGENCE_ENERGY_HPP

#include "implicit_field.hpp"
#include "level_set_field_bundle.hpp"
#include "slab_pierce_stencil.hpp"
#include "slab_prism_operators.hpp"

#include <cassert>
#include <functional>
#include <type_traits>
#include <unordered_map>
#include <vector>

namespace gaudi {
namespace hepworth {
namespace block {

namespace detail {

inline field_eval eval_field_delta_samples(const field_sample &rest,
                                           const field_sample &live) {
  field_eval out;
  out.f = rest.f - live.f;
  out.grad = rest.grad - live.grad;
  if (rest.has_hess && live.has_hess) {
    out.hess = rest.hess - live.hess;
    out.has_hess = true;
  }
  return out;
}

template <FieldBundle B>
inline real bulk_sample_energy_generic(const typename B::q_type &Qr,
                                       const typename B::q_type &Qk,
                                       const vec3 &x_local, real lam_gap) {
  const B br{Qr};
  const B bk{Qk};
  const field_sample er = br.eval(x_local);
  const field_sample ek = bk.eval(x_local);
  const field_eval delta = eval_field_delta_samples(er, ek);
  return slab_field_bulk_density(delta, lam_gap);
}

template <FieldBundle B>
inline vec3 bulk_sample_grad_generic(const typename B::q_type &Qr,
                                     const typename B::q_type &Qk,
                                     const vec3 &x_local, real lam_gap) {
  const B br{Qr};
  const B bk{Qk};
  const field_sample er = br.eval(x_local);
  const field_sample ek = bk.eval(x_local);
  const field_eval delta = eval_field_delta_samples(er, ek);
  return slab_field_bulk_density_grad_wrt_x(delta, lam_gap);
}

template <FieldBundle B>
inline real bulk_sample_energy(const typename B::q_type &Qr,
                                 const typename B::q_type &Qk,
                                 const vec3 &x_local, real lam_gap) {
  if constexpr (std::is_same_v<B, darboux_field_bundle>) {
    return darboux_field_bundle::bulk_sample_energy(Qr, Qk, x_local, lam_gap);
  } else {
    return bulk_sample_energy_generic<B>(Qr, Qk, x_local, lam_gap);
  }
}

template <FieldBundle B>
inline vec3 bulk_sample_grad(const typename B::q_type &Qr,
                             const typename B::q_type &Qk,
                             const vec3 &x_local, real lam_gap) {
  if constexpr (std::is_same_v<B, darboux_field_bundle>) {
    return darboux_field_bundle::bulk_sample_grad(Qr, Qk, x_local, lam_gap);
  } else {
    return bulk_sample_grad_generic<B>(Qr, Qk, x_local, lam_gap);
  }
}

} // namespace detail

/// Sample bulk |∇ΔD|² + λ(ΔD)² on a vertex set with foot-local coords.
/// Q_rest / Q_live are frozen; only sample positions move.
template <FieldBundle B>
inline real eval_shape_field_divergence_energy_on(
    const index_t foot, const typename B::q_type &Q_rest,
    const typename B::q_type &Q_live, const std::vector<index_t> &sample_ids,
    const std::function<vec3(index_t)> &pos_at, real lam_gap) {
  const vec3 x_foot = pos_at(foot);
  real E = 0.0;
  for (index_t vi : sample_ids) {
    const vec3 local = pos_at(vi) - x_foot;
    E += detail::bulk_sample_energy<B>(Q_rest, Q_live, local, lam_gap);
  }
  return E;
}

template <FieldBundle B>
inline real eval_shape_field_divergence_energy(
    const slab_pierce_stencil &stencil, const typename B::q_type &Q_rest,
    const typename B::q_type &Q_live, const slab_ring_positions &ring,
    real lam_gap) {
  return eval_shape_field_divergence_energy_on<B>(
      stencil.foot, Q_rest, Q_live, ring.ids,
      [&](index_t vi) { return ring.at(vi); }, lam_gap);
}

/// Foot one-ring only: quadrature / DOFs for local placement inside frozen ΔD.
template <FieldBundle B>
inline real eval_shape_field_divergence_energy(
    const slab_pierce_stencil &stencil, const typename B::q_type &Q_rest,
    const typename B::q_type &Q_live, const top_ring_positions &top,
    real lam_gap) {
  return eval_shape_field_divergence_energy_on<B>(
      stencil.foot, Q_rest, Q_live, top.ids,
      [&](index_t vi) { return top.at(vi); }, lam_gap);
}

template <FieldBundle B>
inline void shape_field_divergence_grad_on(
    const index_t foot, const typename B::q_type &Q_rest,
    const typename B::q_type &Q_live, const std::vector<index_t> &sample_ids,
    const std::unordered_map<index_t, size_t> &idx,
    const std::function<vec3(index_t)> &pos_at, real lam_gap,
    std::vector<vec3> &grad) {
  grad.assign(sample_ids.size(), vec3::Zero());
  const vec3 x_foot = pos_at(foot);

  const auto add = [&](index_t vi, const vec3 &g) {
    const auto it = idx.find(vi);
    if (it == idx.end())
      return;
    assert(it->second < grad.size());
    grad[it->second] += g;
  };

  vec3 g_foot = vec3::Zero();
  for (index_t vi : sample_ids) {
    const vec3 local = pos_at(vi) - x_foot;
    const vec3 g_local =
        detail::bulk_sample_grad<B>(Q_rest, Q_live, local, lam_gap);
    add(vi, g_local);
    g_foot -= g_local;
  }
  add(foot, g_foot);
}

template <FieldBundle B>
inline void shape_field_divergence_grad(
    const slab_pierce_stencil &stencil, const typename B::q_type &Q_rest,
    const typename B::q_type &Q_live, const slab_ring_positions &ring,
    real lam_gap, std::vector<vec3> &grad) {
  shape_field_divergence_grad_on<B>(
      stencil.foot, Q_rest, Q_live, ring.ids, ring.idx,
      [&](index_t vi) { return ring.at(vi); }, lam_gap, grad);
}

template <FieldBundle B>
inline void shape_field_divergence_grad(
    const slab_pierce_stencil &stencil, const typename B::q_type &Q_rest,
    const typename B::q_type &Q_live, const top_ring_positions &top,
    real lam_gap, std::vector<vec3> &grad) {
  shape_field_divergence_grad_on<B>(
      stencil.foot, Q_rest, Q_live, top.ids, top.idx,
      [&](index_t vi) { return top.at(vi); }, lam_gap, grad);
}

template <typename RingPos, FieldBundle B>
inline RingPos resolve_shape_field_divergence_gd_ring(
    const slab_pierce_stencil &stencil, const typename B::q_type &Q_rest,
    const typename B::q_type &Q_live, const RingPos &x0, real lam_gap,
    int iters, real step_scale = 0.25) {
  RingPos x = x0;
  if (iters <= 0)
    return x;

  for (int iter = 0; iter < iters; ++iter) {
    std::vector<vec3> grad;
    shape_field_divergence_grad<B>(stencil, Q_rest, Q_live, x, lam_gap, grad);
    real g2 = 0.0;
    bool grad_ok = true;
    for (const vec3 &g : grad) {
      if (!g.allFinite()) {
        grad_ok = false;
        break;
      }
      g2 += g.squaredNorm();
    }
    if (!grad_ok)
      return x0;
    if (!(g2 > 1e-32) || !std::isfinite(g2))
      break;

    RingPos trial = x;
    real alpha = step_scale / std::sqrt(g2);
    const real e0 =
        eval_shape_field_divergence_energy<B>(stencil, Q_rest, Q_live, x, lam_gap);
    if (!std::isfinite(e0))
      return x0;
    bool accepted = false;
    for (int ls = 0; ls < 8; ++ls) {
      bool trial_ok = true;
      for (size_t i = 0; i < x.x.size(); ++i) {
        trial.x[i] = x.x[i] - alpha * grad[i];
        if (!trial.x[i].allFinite()) {
          trial_ok = false;
          break;
        }
      }
      if (!trial_ok) {
        alpha *= 0.5;
        continue;
      }
      const real e1 = eval_shape_field_divergence_energy<B>(
          stencil, Q_rest, Q_live, trial, lam_gap);
      if (std::isfinite(e1) && e1 < e0) {
        x = trial;
        accepted = true;
        break;
      }
      alpha *= 0.5;
    }
    if (!accepted || e0 < 1e-24)
      break;
  }
  return x;
}

template <FieldBundle B>
inline slab_ring_positions resolve_shape_field_divergence_gd(
    const slab_pierce_stencil &stencil, const typename B::q_type &Q_rest,
    const typename B::q_type &Q_live, const slab_ring_positions &x0,
    real lam_gap, int iters, real step_scale = 0.25) {
  return resolve_shape_field_divergence_gd_ring<slab_ring_positions, B>(
      stencil, Q_rest, Q_live, x0, lam_gap, iters, step_scale);
}

template <FieldBundle B>
inline top_ring_positions resolve_shape_field_divergence_gd(
    const slab_pierce_stencil &stencil, const typename B::q_type &Q_rest,
    const typename B::q_type &Q_live, const top_ring_positions &x0,
    real lam_gap, int iters, real step_scale = 0.25) {
  return resolve_shape_field_divergence_gd_ring<top_ring_positions, B>(
      stencil, Q_rest, Q_live, x0, lam_gap, iters, step_scale);
}

// Non-template entry points for mixed-family slab_field_coeffs storage.

namespace detail {

inline real eval_shape_field_divergence_energy_generic(
    const index_t foot, const slab_field_coeffs &Q_rest,
    const slab_field_coeffs &Q_live, const std::vector<index_t> &sample_ids,
    const std::function<vec3(index_t)> &pos_at, real lam_gap) {
  const vec3 x_foot = pos_at(foot);
  real E = 0.0;
  for (index_t vi : sample_ids) {
    const vec3 local = pos_at(vi) - x_foot;
    const field_eval er = eval_slab_field(Q_rest, local);
    const field_eval ek = eval_slab_field(Q_live, local);
    E += slab_field_bulk_density(eval_field_delta(er, ek), lam_gap);
  }
  return E;
}

inline void shape_field_divergence_grad_generic(
    const index_t foot, const slab_field_coeffs &Q_rest,
    const slab_field_coeffs &Q_live, const std::vector<index_t> &sample_ids,
    const std::unordered_map<index_t, size_t> &idx,
    const std::function<vec3(index_t)> &pos_at, real lam_gap,
    std::vector<vec3> &grad) {
  grad.assign(sample_ids.size(), vec3::Zero());
  const vec3 x_foot = pos_at(foot);
  const auto add = [&](index_t vi, const vec3 &g) {
    const auto it = idx.find(vi);
    if (it == idx.end())
      return;
    assert(it->second < grad.size());
    grad[it->second] += g;
  };
  vec3 g_foot = vec3::Zero();
  for (index_t vi : sample_ids) {
    const vec3 local = pos_at(vi) - x_foot;
    const field_eval er = eval_slab_field(Q_rest, local);
    const field_eval ek = eval_slab_field(Q_live, local);
    const vec3 g_local = slab_field_bulk_density_grad_wrt_x(
        eval_field_delta(er, ek), lam_gap);
    add(vi, g_local);
    g_foot -= g_local;
  }
  add(foot, g_foot);
}

template <typename RingPos>
inline real eval_shape_field_energy_coeffs(const slab_pierce_stencil &stencil,
                                           const slab_field_coeffs &Q_rest,
                                           const slab_field_coeffs &Q_live,
                                           const RingPos &ring, real lam_gap) {
  switch (Q_rest.family) {
  case slab_field_family::darboux:
    return eval_shape_field_divergence_energy<darboux_field_bundle>(
        stencil, Q_rest.darboux, Q_live.darboux, ring, lam_gap);
  case slab_field_family::rbf_mls:
    return eval_shape_field_divergence_energy<rbf_field_bundle>(
        stencil, Q_rest.rbf, Q_live.rbf, ring, lam_gap);
  case slab_field_family::quadric:
  default:
    return eval_shape_field_divergence_energy_generic(
        stencil.foot, Q_rest, Q_live, ring.ids,
        [&](index_t vi) { return ring.at(vi); }, lam_gap);
  }
}

template <typename RingPos>
inline void shape_field_grad_coeffs(const slab_pierce_stencil &stencil,
                                    const slab_field_coeffs &Q_rest,
                                    const slab_field_coeffs &Q_live,
                                    const RingPos &ring, real lam_gap,
                                    std::vector<vec3> &grad) {
  switch (Q_rest.family) {
  case slab_field_family::darboux:
    shape_field_divergence_grad<darboux_field_bundle>(
        stencil, Q_rest.darboux, Q_live.darboux, ring, lam_gap, grad);
    return;
  case slab_field_family::rbf_mls:
    shape_field_divergence_grad<rbf_field_bundle>(
        stencil, Q_rest.rbf, Q_live.rbf, ring, lam_gap, grad);
    return;
  case slab_field_family::quadric:
  default:
    shape_field_divergence_grad_generic(
        stencil.foot, Q_rest, Q_live, ring.ids, ring.idx,
        [&](index_t vi) { return ring.at(vi); }, lam_gap, grad);
    return;
  }
}

template <typename RingPos>
inline RingPos resolve_shape_field_gd_coeffs(const slab_pierce_stencil &stencil,
                                             const slab_field_coeffs &Q_rest,
                                             const slab_field_coeffs &Q_live,
                                             const RingPos &x0, real lam_gap,
                                             int iters, real step_scale) {
  switch (Q_rest.family) {
  case slab_field_family::darboux:
    return resolve_shape_field_divergence_gd<darboux_field_bundle>(
        stencil, Q_rest.darboux, Q_live.darboux, x0, lam_gap, iters, step_scale);
  case slab_field_family::rbf_mls:
    return resolve_shape_field_divergence_gd<rbf_field_bundle>(
        stencil, Q_rest.rbf, Q_live.rbf, x0, lam_gap, iters, step_scale);
  case slab_field_family::quadric:
  default: {
    RingPos x = x0;
    if (iters <= 0)
      return x;
    for (int iter = 0; iter < iters; ++iter) {
      std::vector<vec3> grad;
      shape_field_grad_coeffs(stencil, Q_rest, Q_live, x, lam_gap, grad);
      real g2 = 0.0;
      bool grad_ok = true;
      for (const vec3 &g : grad) {
        if (!g.allFinite()) {
          grad_ok = false;
          break;
        }
        g2 += g.squaredNorm();
      }
      if (!grad_ok)
        return x0;
      if (!(g2 > 1e-32) || !std::isfinite(g2))
        break;

      RingPos trial = x;
      real alpha = step_scale / std::sqrt(g2);
      const real e0 =
          eval_shape_field_energy_coeffs(stencil, Q_rest, Q_live, x, lam_gap);
      if (!std::isfinite(e0))
        return x0;
      bool accepted = false;
      for (int ls = 0; ls < 8; ++ls) {
        bool trial_ok = true;
        for (size_t i = 0; i < x.x.size(); ++i) {
          trial.x[i] = x.x[i] - alpha * grad[i];
          if (!trial.x[i].allFinite()) {
            trial_ok = false;
            break;
          }
        }
        if (!trial_ok) {
          alpha *= 0.5;
          continue;
        }
        const real e1 = eval_shape_field_energy_coeffs(stencil, Q_rest, Q_live,
                                                       trial, lam_gap);
        if (std::isfinite(e1) && e1 < e0) {
          x = trial;
          accepted = true;
          break;
        }
        alpha *= 0.5;
      }
      if (!accepted || e0 < 1e-24)
        break;
    }
    return x;
  }
  }
}

} // namespace detail

inline real eval_shape_field_divergence_energy(
    const slab_pierce_stencil &stencil, const slab_field_coeffs &Q_rest,
    const slab_field_coeffs &Q_live, const slab_ring_positions &ring,
    real lam_gap) {
  return detail::eval_shape_field_energy_coeffs(stencil, Q_rest, Q_live, ring,
                                                lam_gap);
}

inline real eval_shape_field_divergence_energy(
    const slab_pierce_stencil &stencil, const slab_field_coeffs &Q_rest,
    const slab_field_coeffs &Q_live, const top_ring_positions &top,
    real lam_gap) {
  return detail::eval_shape_field_energy_coeffs(stencil, Q_rest, Q_live, top,
                                                lam_gap);
}

inline void shape_field_divergence_grad(const slab_pierce_stencil &stencil,
                                        const slab_field_coeffs &Q_rest,
                                        const slab_field_coeffs &Q_live,
                                        const slab_ring_positions &ring,
                                        real lam_gap,
                                        std::vector<vec3> &grad) {
  detail::shape_field_grad_coeffs(stencil, Q_rest, Q_live, ring, lam_gap, grad);
}

inline void shape_field_divergence_grad(const slab_pierce_stencil &stencil,
                                        const slab_field_coeffs &Q_rest,
                                        const slab_field_coeffs &Q_live,
                                        const top_ring_positions &top,
                                        real lam_gap,
                                        std::vector<vec3> &grad) {
  detail::shape_field_grad_coeffs(stencil, Q_rest, Q_live, top, lam_gap, grad);
}

inline slab_ring_positions resolve_shape_field_divergence_gd(
    const slab_pierce_stencil &stencil, const slab_field_coeffs &Q_rest,
    const slab_field_coeffs &Q_live, const slab_ring_positions &x0,
    real lam_gap, int iters, real step_scale = 0.25) {
  return detail::resolve_shape_field_gd_coeffs(stencil, Q_rest, Q_live, x0,
                                               lam_gap, iters, step_scale);
}

inline top_ring_positions resolve_shape_field_divergence_gd(
    const slab_pierce_stencil &stencil, const slab_field_coeffs &Q_rest,
    const slab_field_coeffs &Q_live, const top_ring_positions &x0,
    real lam_gap, int iters, real step_scale = 0.25) {
  return detail::resolve_shape_field_gd_coeffs(stencil, Q_rest, Q_live, x0,
                                               lam_gap, iters, step_scale);
}

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // GAUDI_HEPWORTH_SLAB_FIELD_DIVERGENCE_ENERGY_HPP
