#ifndef GAUDI_ALBERS_SHAPE_OPERATOR_CONSTRAINED_FIT_HPP
#define GAUDI_ALBERS_SHAPE_OPERATOR_CONSTRAINED_FIT_HPP

#include "gaudi/albers/darboux_cyclide.hpp"
#include "gaudi/albers/ncls.hpp"
#include "gaudi/albers/quadric.hpp"
#include "gaudi/common.h"
#include "shape_operator_fit_generated.hpp"

namespace gaudi {
namespace albers {

/// Normal-constrained quadric + optional linear shape-operator tip-in:
/// \(\mathrm{vech}(P_n H P_n) \approx \mathrm{vech}(W_\star)\) via generated \(M(n)\).
class shape_operator_constrained_quadric {
public:
  using coefficients = vec10;

  shape_operator_constrained_quadric() {
    A = mat10::Zero();
    b = vec10::Zero();
  }

  void set_shape_weight(real w) { shape_weight_ = w; }
  real shape_weight() const { return shape_weight_; }

  void accumulate(real w, const vec3 &x, const vec3 &N,
                  const mat3 &W_star = mat3::Zero()) {
    const mat410 Ab = mk_quad_A(x);
    const vec4 Nb = mk_N(N);
    A += w * Ab.transpose() * Ab;
    b += w * Ab.transpose() * Nb;

    if (shape_weight_ <= 0.0 || W_star.isZero(0.0)) {
      return;
    }
    const vec3 n = N.normalized();
    if (n.norm() < 1e-12) {
      return;
    }
    const Eigen::Matrix<real, 6, 10> M =
        medial_generated::mk_quad_W_rows(n);
    accumulate_hessian_block3(shape_weight_, w, M, W_star, A, b);
  }

  vec10 solve() { return A.colPivHouseholderQr().solve(b); }

  mat10 A;
  vec10 b;

private:
  real shape_weight_ = 0.0;
};

/// Normal-constrained Darboux + optional linear \(W\) tip-in at sample \(x\).
class shape_operator_constrained_darboux_cyclide {
public:
  using coefficients = vec14;

  shape_operator_constrained_darboux_cyclide() {
    A = mat14::Zero();
    b = vec14::Zero();
  }

  void set_shape_weight(real w) { shape_weight_ = w; }
  real shape_weight() const { return shape_weight_; }

  void accumulate(real w, const vec3 &x, const vec3 &N,
                  const mat3 &W_star = mat3::Zero()) {
    const mat414 Ab = mk_darboux_A(x);
    const vec4 Nb = mk_N(N);
    A += w * Ab.transpose() * Ab;
    b += w * Ab.transpose() * Nb;

    if (shape_weight_ <= 0.0 || W_star.isZero(0.0)) {
      return;
    }
    const vec3 n = N.normalized();
    if (n.norm() < 1e-12) {
      return;
    }
    const Eigen::Matrix<real, 6, 14> M =
        medial_generated::mk_darboux_W_rows(n, x);
    accumulate_hessian_block3(shape_weight_, w, M, W_star, A, b);
  }

  vec14 solve() { return A.colPivHouseholderQr().solve(b); }

  mat14 A;
  vec14 b;

private:
  real shape_weight_ = 0.0;
};

} // namespace albers
} // namespace gaudi

#endif // GAUDI_ALBERS_SHAPE_OPERATOR_CONSTRAINED_FIT_HPP
