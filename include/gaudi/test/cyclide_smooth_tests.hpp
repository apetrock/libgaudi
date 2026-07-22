#ifndef __GAUDI_TEST_CYCLIDE_SMOOTH_TESTS_HPP__
#define __GAUDI_TEST_CYCLIDE_SMOOTH_TESTS_HPP__

#include <cmath>

#include "cyclide_smooth_generated.hpp"
#include "gaudi/albers/darboux_cyclide.hpp"
#include "gaudi/test/test.hpp"

namespace gaudi {
namespace test {

namespace {

inline void fill_cyclide_smooth_probe(albers::vec14 &Qi, albers::vec14 &Q0,
                                      albers::vec14 &Qj) {
  for (int i = 0; i < 14; ++i) {
    Qi[i] = real(0.3) + real(0.02) * real(i);
    Q0[i] = real(0.07) * real(i);
    Qj[i] = real(0.25) + real(0.015) * real(i);
  }
}

inline real central_diff_anchor_energy(const albers::vec14 &Qi,
                                       const albers::vec14 &Q0, real wi, int k,
                                       real eps = real(1e-6)) {
  albers::vec14 up = Qi;
  albers::vec14 dn = Qi;
  up[k] += eps;
  dn[k] -= eps;
  return (albers::medial_generated::cyclide_smooth_anchor_energy(up, Q0, wi) -
          albers::medial_generated::cyclide_smooth_anchor_energy(dn, Q0, wi)) /
         (real(2) * eps);
}

inline real central_diff_neighbor_energy(const albers::vec14 &Qi,
                                         const albers::vec14 &Qj,
                                         const vec3 &x_j_in_i, real wi,
                                         real w_ij, real alpha_G, real alpha_H,
                                         int k, real eps = real(1e-6)) {
  albers::vec14 up = Qi;
  albers::vec14 dn = Qi;
  up[k] += eps;
  dn[k] -= eps;
  return (albers::medial_generated::cyclide_smooth_neighbor_energy(
              up, Qj, x_j_in_i, wi, w_ij, alpha_G, alpha_H) -
          albers::medial_generated::cyclide_smooth_neighbor_energy(
              dn, Qj, x_j_in_i, wi, w_ij, alpha_G, alpha_H)) /
         (real(2) * eps);
}

} // namespace

GAUDI_TEST(cyclide_smooth_generated_anchor_finite) {
  albers::vec14 Qi;
  albers::vec14 Q0;
  albers::vec14 Qj;
  fill_cyclide_smooth_probe(Qi, Q0, Qj);

  const real wi = real(0.85);
  const real E =
      albers::medial_generated::cyclide_smooth_anchor_energy(Qi, Q0, wi);
  GAUDI_ASSERT(std::isfinite(E));

  albers::vec14 grad = albers::vec14::Zero();
  albers::medial_generated::cyclide_smooth_anchor_grad(Qi, Q0, wi, grad);
  GAUDI_ASSERT(grad.allFinite());
}

GAUDI_TEST(cyclide_smooth_generated_neighbor_finite) {
  albers::vec14 Qi;
  albers::vec14 Q0;
  albers::vec14 Qj;
  fill_cyclide_smooth_probe(Qi, Q0, Qj);

  const vec3 x_j_in_i(real(0.12), real(-0.05), real(0.08));
  const real wi = real(0.85);
  const real w_ij = real(0.5);
  const real alpha_G = real(1.0);
  const real alpha_H = real(1.0);

  const real E = albers::medial_generated::cyclide_smooth_neighbor_energy(
      Qi, Qj, x_j_in_i, wi, w_ij, alpha_G, alpha_H);
  GAUDI_ASSERT(std::isfinite(E));

  albers::vec14 grad = albers::vec14::Zero();
  albers::medial_generated::cyclide_smooth_neighbor_grad(
      Qi, Qj, x_j_in_i, wi, w_ij, alpha_G, alpha_H, grad);
  GAUDI_ASSERT(grad.allFinite());
}

GAUDI_TEST(cyclide_smooth_generated_anchor_grad_numeric) {
  albers::vec14 Qi;
  albers::vec14 Q0;
  albers::vec14 Qj;
  fill_cyclide_smooth_probe(Qi, Q0, Qj);

  const real wi = real(0.85);
  albers::vec14 grad = albers::vec14::Zero();
  albers::medial_generated::cyclide_smooth_anchor_grad(Qi, Q0, wi, grad);

  for (int k = 0; k < 14; ++k) {
    const real num = central_diff_anchor_energy(Qi, Q0, wi, k);
    GAUDI_ASSERT(std::abs(num - grad[k]) < real(1e-4));
  }
}

GAUDI_TEST(cyclide_smooth_generated_neighbor_grad_numeric) {
  albers::vec14 Qi;
  albers::vec14 Q0;
  albers::vec14 Qj;
  fill_cyclide_smooth_probe(Qi, Q0, Qj);

  const vec3 x_j_in_i(real(0.12), real(-0.05), real(0.08));
  const real wi = real(0.85);
  const real w_ij = real(0.5);
  const real alpha_G = real(1.0);
  const real alpha_H = real(1.0);

  albers::vec14 grad = albers::vec14::Zero();
  albers::medial_generated::cyclide_smooth_neighbor_grad(
      Qi, Qj, x_j_in_i, wi, w_ij, alpha_G, alpha_H, grad);

  for (int k = 0; k < 14; ++k) {
    const real num = central_diff_neighbor_energy(Qi, Qj, x_j_in_i, wi, w_ij,
                                                  alpha_G, alpha_H, k);
    GAUDI_ASSERT(std::abs(num - grad[k]) < real(2e-3));
  }
}

} // namespace test
} // namespace gaudi

#endif // __GAUDI_TEST_CYCLIDE_SMOOTH_TESTS_HPP__
