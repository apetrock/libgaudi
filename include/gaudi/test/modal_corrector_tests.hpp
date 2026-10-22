#ifndef GAUDI_TEST_MODAL_CORRECTOR_TESTS_HPP
#define GAUDI_TEST_MODAL_CORRECTOR_TESTS_HPP

#include "gaudi/hepworth/block/modal_corrector.hpp"
#include "gaudi/test/test.hpp"

#include <Eigen/Dense>

#include <random>
#include <vector>

namespace gaudi {
namespace test {

namespace modal_detail {

/// A sparse SPD preconditioner: a 1D Laplacian plus a mass diagonal.
inline matS laplace_mass(int n, real mass) {
  std::vector<Eigen::Triplet<real>> t;
  for (int i = 0; i < n; ++i) {
    t.emplace_back(i, i, 2 + mass);
    if (i + 1 < n) {
      t.emplace_back(i, i + 1, -1);
      t.emplace_back(i + 1, i, -1);
    }
  }
  matS P(n, n);
  P.setFromTriplets(t.begin(), t.end());
  return P;
}

/// `k` vectors orthonormal in the `P` inner product.
inline Eigen::MatrixXd p_orthonormal(const Eigen::MatrixXd &P, int k, unsigned seed) {
  std::mt19937 rng(seed);
  std::normal_distribution<real> g(0, 1);
  Eigen::MatrixXd U(P.rows(), k);
  for (int j = 0; j < k; ++j) {
    vecX u(P.rows());
    for (Eigen::Index i = 0; i < u.size(); ++i)
      u[i] = g(rng);
    for (int m = 0; m < j; ++m)
      u -= U.col(m).dot(P * u) * U.col(m);
    U.col(j) = u / std::sqrt(u.dot(P * u));
  }
  return U;
}

/// `H = P - sum (1 - lambda_j) (P u_j)(P u_j)^T`: generalized eigenvalue
/// `lambda_j` on each `u_j`, 1 on their P-complement.
inline Eigen::MatrixXd with_modes(const Eigen::MatrixXd &P, const Eigen::MatrixXd &U,
                                  const std::vector<real> &lambda) {
  Eigen::MatrixXd H = P;
  for (size_t j = 0; j < lambda.size(); ++j) {
    const vecX Pu = P * U.col(static_cast<Eigen::Index>(j));
    H -= (real(1) - lambda[j]) * Pu * Pu.transpose();
  }
  return H;
}

} // namespace modal_detail

/// The modes from the LDLT-transformed Lanczos match a dense generalized
/// eigensolve, and come out P-orthonormal.
GAUDI_TEST(modal_corrector_matches_dense_generalized_eigs) {
  using namespace modal_detail;
  const int n = 80;
  const matS P = laplace_mass(n, 0.1);
  const Eigen::MatrixXd Pd(P);
  const Eigen::MatrixXd U = p_orthonormal(Pd, 3, 3);
  const Eigen::MatrixXd H = with_modes(Pd, U, {-0.4, 0.05, 0.3});
  matS Pc = P;
  m_solver S(Pc);
  const hepworth::block::hessian_vector_fn hv = [&](const vecX &v) { return vecX(H * v); };
  hepworth::block::modal_config cfg;
  cfg.modes = 4;
  hepworth::block::modal_corrector mc;
  GAUDI_ASSERT(mc.refresh(S, hv, cfg));
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges(H, Pd);
  std::vector<real> got(mc.lambdas().data(), mc.lambdas().data() + mc.lambdas().size());
  std::sort(got.begin(), got.end());
  real worst = 0;
  for (size_t i = 0; i < got.size(); ++i)
    worst = std::max(worst, std::abs(got[i] - ges.eigenvalues()[static_cast<Eigen::Index>(i)]));
  const Eigen::MatrixXd G = mc.modes().transpose() * Pd * mc.modes();
  const real ortho = (G - Eigen::MatrixXd::Identity(G.rows(), G.cols())).norm();
  console_logger::info << "modal: lambda " << got[0] << " " << got[1] << " " << got[2] << " "
                       << got[3] << ", worst gap to dense " << worst << ", P-orthonormality "
                       << ortho << std::endl;
  GAUDI_EXPECT(worst < real(1e-6));
  GAUDI_EXPECT(ortho < real(1e-6));
  GAUDI_EXPECT(mc.last().negative == 1);
}

/// On a quadratic whose preconditioner overestimates three directions
/// (generalized eigenvalues 0.01, 0.05, 0.2), the corrected iteration
/// converges in far fewer iterates than the plain one.
GAUDI_TEST(modal_corrector_accelerates_slow_modes) {
  using namespace modal_detail;
  const int n = 80;
  const matS P = laplace_mass(n, 0.1);
  const Eigen::MatrixXd Pd(P);
  const Eigen::MatrixXd U = p_orthonormal(Pd, 3, 5);
  const Eigen::MatrixXd H = with_modes(Pd, U, {0.01, 0.05, 0.2});
  const vecX f = vecX::Ones(n);
  auto energy = [&](const vecX &q) { return real(0.5) * q.dot(H * q) - f.dot(q); };
  matS Pc = P;
  m_solver S(Pc);
  const hepworth::block::hessian_vector_fn hv = [&](const vecX &v) { return vecX(H * v); };
  auto run = [&](bool modal) {
    hepworth::block::modal_config cfg;
    cfg.modes = 3;
    hepworth::block::modal_corrector mc;
    if (modal)
      mc.refresh(S, hv, cfg);
    vecX q = vecX::Zero(n);
    for (int k = 0; k < 5000; ++k) {
      const vecX r = H * q - f;
      if (r.norm() < real(1e-10) * f.norm())
        return k;
      vecX b = Pd * q - r;
      vecX qn = S.solve(b);
      if (modal)
        qn = mc.correct(qn, r, energy, energy(qn), real(1), cfg);
      q = qn;
    }
    return 5000;
  };
  const int plain = run(false), modal = run(true);
  console_logger::info << "modal: iterates to 1e-10, plain " << plain << ", modal " << modal
                       << std::endl;
  GAUDI_EXPECT(modal * 10 < plain);
}

/// A symmetric double well along one mode (generalized eigenvalue -0.5 at
/// the saddle; two stable modes at 0.6 and 0.8 keep the spectrum from
/// collapsing): the plain iteration started at the saddle stays there; the
/// corrector steps out and settles in a well.
GAUDI_TEST(modal_corrector_escapes_a_saddle) {
  using namespace modal_detail;
  const int n = 60;
  const matS P = laplace_mass(n, 0.1);
  const Eigen::MatrixXd Pd(P);
  const Eigen::MatrixXd U = p_orthonormal(Pd, 3, 7);
  const Eigen::MatrixXd H = with_modes(Pd, U, {-0.5, 0.6, 0.8});
  const vecX Pu = Pd * U.col(0);
  const real beta = 2;
  auto energy = [&](const vecX &q) {
    const real t = Pu.dot(q);
    return real(0.5) * q.dot(H * q) + real(0.25) * beta * t * t * t * t;
  };
  auto grad = [&](const vecX &q) {
    const real t = Pu.dot(q);
    return vecX(H * q + beta * t * t * t * Pu);
  };
  matS Pc = P;
  m_solver S(Pc);
  const hepworth::block::hessian_vector_fn hv = [&](const vecX &v) { return vecX(H * v); };
  auto run = [&](bool modal) {
    hepworth::block::modal_config cfg;
    cfg.modes = 2;
    hepworth::block::modal_corrector mc;
    if (modal)
      mc.refresh(S, hv, cfg);
    vecX q = vecX::Zero(n);
    for (int k = 0; k < 200; ++k) {
      const vecX r = grad(q);
      vecX b = Pd * q - r;
      vecX qn = S.solve(b);
      if (modal && k == 0)
        qn = mc.correct(qn, r, energy, energy(qn), real(1), cfg);
      q = qn;
    }
    return q;
  };
  const vecX q_plain = run(false), q_modal = run(true);
  const real t_well = std::sqrt(real(0.5) / beta);
  console_logger::info << "modal: saddle amplitude plain " << Pu.dot(q_plain) << ", modal "
                       << Pu.dot(q_modal) << " (well " << t_well << "), energy "
                       << energy(q_plain) << " -> " << energy(q_modal) << std::endl;
  GAUDI_EXPECT(std::abs(Pu.dot(q_plain)) < real(1e-12));
  GAUDI_EXPECT(std::abs(std::abs(Pu.dot(q_modal)) - t_well) < real(1e-6));
  GAUDI_EXPECT(energy(q_modal) < energy(q_plain));
}

} // namespace test
} // namespace gaudi

#endif
