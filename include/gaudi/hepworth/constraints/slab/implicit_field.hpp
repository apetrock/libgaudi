#ifndef GAUDI_HEPWORTH_IMPLICIT_FIELD_HPP
#define GAUDI_HEPWORTH_IMPLICIT_FIELD_HPP

#include "gaudi/albers/darboux_cyclide.hpp"
#include "gaudi/albers/quadric.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

/// Value + first (+ optional second) derivatives w.r.t. evaluation point (foot-local coords).
struct field_eval {
  real f = 0.0;
  vec3 grad = vec3::Zero();
  mat3 hess = mat3::Zero();
  bool has_hess = false;
};

inline real eval_quadric(const albers::vec10 &Q, const vec3 &x) {
  return Q[0] * x[0] * x[0] + Q[1] * x[1] * x[1] + Q[2] * x[2] * x[2] +
         Q[3] * x[0] * x[1] + Q[4] * x[0] * x[2] + Q[5] * x[1] * x[2] +
         2.0 * Q[6] * x[0] + 2.0 * Q[7] * x[1] + 2.0 * Q[8] * x[2] + Q[9];
}

inline field_eval eval_darboux_field(const albers::vec14 &Q, const vec3 &x) {
  field_eval out;
  out.f = albers::eval_darboux(Q, x);
  out.grad = albers::darboux_grad(Q, x);
  out.hess = albers::darboux_hessian(Q, x);
  out.has_hess = true;
  return out;
}

inline field_eval eval_quadric_field(const albers::vec10 &Q, const vec3 &x) {
  field_eval out;
  out.f = eval_quadric(Q, x);
  out.grad = albers::quadric_grad(Q, x);
  out.hess = 2.0 * albers::quadric_hessian(Q);
  out.has_hess = true;
  return out;
}

inline field_eval eval_field_delta(const field_eval &rest, const field_eval &live) {
  field_eval out;
  out.f = rest.f - live.f;
  out.grad = rest.grad - live.grad;
  if (rest.has_hess && live.has_hess) {
    out.hess = rest.hess - live.hess;
    out.has_hess = true;
  }
  return out;
}

/// Bulk density |∇(ΔD)|² + λ_gap (ΔD)² at one sample (ΔD = f_rest − f_live).
inline real slab_field_bulk_density(const field_eval &delta, real lam_gap) {
  return delta.grad.squaredNorm() + lam_gap * delta.f * delta.f;
}

/// ∂density/∂x at the sample point (x is foot-local).
inline vec3 slab_field_bulk_density_grad_wrt_x(const field_eval &delta,
                                               real lam_gap) {
  if (!delta.has_hess)
    return 2.0 * lam_gap * delta.f * delta.grad;
  return 2.0 * delta.hess * delta.grad + 2.0 * lam_gap * delta.f * delta.grad;
}

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // GAUDI_HEPWORTH_IMPLICIT_FIELD_HPP
