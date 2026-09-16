#ifndef GAUDI_HEPWORTH_LEVEL_SET_FIELD_BUNDLE_HPP
#define GAUDI_HEPWORTH_LEVEL_SET_FIELD_BUNDLE_HPP

#include "implicit_field.hpp"
#include "local_slab_field_fit.hpp"
#include "slab_pierce_stencil.hpp"

#include "slab_field_energy_generated.hpp"

#include <concepts>
#include <functional>

namespace gaudi {
namespace hepworth {
namespace block {

/// Value + derivatives of D(Q, x) w.r.t. sample point x (foot-local).
struct field_sample {
  real f = 0.0;
  vec3 grad = vec3::Zero();
  mat3 hess = mat3::Zero();
  bool has_hess = false;
};

inline field_sample from_field_eval(const field_eval &e) {
  field_sample s;
  s.f = e.f;
  s.grad = e.grad;
  s.hess = e.hess;
  s.has_hess = e.has_hess;
  return s;
}

/// Rest/live Darboux coefficients on a slab stencil.
struct darboux_field_bundle {
  using q_type = albers::vec14;

  albers::vec14 Q = albers::vec14::Zero();

  static constexpr slab_field_family family() {
    return slab_field_family::darboux;
  }

  static darboux_field_bundle from_coeffs(const slab_field_coeffs &c) {
    darboux_field_bundle b;
    b.Q = c.darboux;
    return b;
  }

  slab_field_coeffs to_coeffs() const {
    slab_field_coeffs c;
    c.family = slab_field_family::darboux;
    c.darboux = Q;
    return c;
  }

  static darboux_field_bundle
  fit_Q(const slab_pierce_stencil &stencil,
        const std::function<vec3(index_t)> &pos_at,
        const local_field_fit_config &cfg, real /*avg_edge_length*/) {
    darboux_field_bundle b;
    b.Q = fit_darboux_on_stencil(stencil, pos_at, cfg);
    return b;
  }

  field_sample eval(const vec3 &x_local) const {
    return from_field_eval(eval_darboux_field(Q, x_local));
  }

  static real bulk_sample_energy(const q_type &Qr, const q_type &Qk,
                                 const vec3 &x_local, real lam_gap) {
    return albers::medial_generated::slab_bulk_energy(Qr, Qk, x_local, lam_gap);
  }

  static vec3 bulk_sample_grad(const q_type &Qr, const q_type &Qk,
                               const vec3 &x_local, real lam_gap) {
    return albers::medial_generated::slab_bulk_energy_grad(Qr, Qk, x_local,
                                                          lam_gap);
  }
};

/// Local Gaussian RBF + linear polynomial MLS on stencil samples.
struct rbf_field_bundle {
  using q_type = local_rbf_mls_field;

  local_rbf_mls_field Q;

  static constexpr slab_field_family family() {
    return slab_field_family::rbf_mls;
  }

  static rbf_field_bundle from_coeffs(const slab_field_coeffs &c) {
    rbf_field_bundle b;
    b.Q = c.rbf;
    return b;
  }

  slab_field_coeffs to_coeffs() const {
    slab_field_coeffs c;
    c.family = slab_field_family::rbf_mls;
    c.rbf = Q;
    return c;
  }

  static rbf_field_bundle
  fit_Q(const slab_pierce_stencil &stencil,
        const std::function<vec3(index_t)> &pos_at,
        const local_field_fit_config &cfg, real avg_edge_length,
        const local_rbf_mls_field *rest_basis = nullptr) {
    rbf_field_bundle b;
    b.Q = fit_rbf_mls_on_stencil(
        stencil, pos_at, cfg, avg_edge_length,
        rest_basis ? &rest_basis->centers : nullptr,
        rest_basis ? rest_basis->sigma : 0.0);
    return b;
  }

  field_sample eval(const vec3 &x_local) const {
    return from_field_eval(Q.eval(x_local));
  }
};

template <typename B>
concept FieldBundle = requires(
    const B &b, const typename B::q_type &qr, const typename B::q_type &qk,
    const vec3 &x, const slab_pierce_stencil &stencil,
    const std::function<vec3(index_t)> &pos_at, const local_field_fit_config &cfg,
    real avg_edge) {
  typename B::q_type;
  { B::family() } -> std::same_as<slab_field_family>;
  { B::fit_Q(stencil, pos_at, cfg, avg_edge) } -> std::same_as<B>;
  { b.eval(x) } -> std::same_as<field_sample>;
  { b.to_coeffs() } -> std::same_as<slab_field_coeffs>;
};

static_assert(FieldBundle<darboux_field_bundle>);
static_assert(FieldBundle<rbf_field_bundle>);

inline slab_field_coeffs fit_field_coeffs_on_stencil(
    slab_field_family family, const slab_pierce_stencil &stencil,
    const std::function<vec3(index_t)> &pos_at,
    const local_field_fit_config &cfg, real avg_edge_length,
    const slab_field_coeffs *rest_basis = nullptr) {
  switch (family) {
  case slab_field_family::darboux:
    return darboux_field_bundle::fit_Q(stencil, pos_at, cfg, avg_edge_length)
        .to_coeffs();
  case slab_field_family::quadric: {
    slab_field_coeffs c;
    c.family = slab_field_family::quadric;
    c.quadric = fit_quadric_on_stencil(stencil, pos_at, cfg);
    return c;
  }
  case slab_field_family::rbf_mls:
    return rbf_field_bundle::fit_Q(
               stencil, pos_at, cfg, avg_edge_length,
               rest_basis && rest_basis->family == slab_field_family::rbf_mls
                   ? &rest_basis->rbf
                   : nullptr)
        .to_coeffs();
  }
  return {};
}

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // GAUDI_HEPWORTH_LEVEL_SET_FIELD_BUNDLE_HPP
