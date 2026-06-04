/*
 * Sparse Laplacian spectrum via Spectra.
 *
 * Symmetrize the weak cotan stiffness \f$L\f$ first, then **negate** to obtain a
 * **positive** semi-definite operator \f$-L\f$ when \f$L\f$ is negative semi-definite.
 * Regularize \f$A = -L_{\mathrm{sym}} + \varepsilon I\f$ and use **SmallestAlge** for
 * low-frequency modes (near-null of \f$L\f$, i.e. smallest eigenvalues of \f$A\f$).
 *
 * Do **not** apply **LargestAlge** directly to \f$L+\varepsilon I\f$ when \f$L\f$ may
 * be indefinite (signed cotan): the algebraic maximum is then unrelated to smooth modes.
 *
 * **Largest magnitude** on the same \f$A\f$ targets the high-frequency tail of \f$-L\f$.
 *
 * **Shift–invert** with \f$\sigma\f$ and `LargestMagn` on \f$(A-\sigma I)^{-1}\f$ finds
 * eigenvalues of \f$A\f$ **closest** to \f$\sigma\f$ (mid-band / interior spectrum).
 *
 * Implementation: thin wrappers around `spectrum.hpp` (`sparse_sym_eigs_*`).
 */

#ifndef __GAUDI_KUSAMA_LAPLACE_SPECTRUM__
#define __GAUDI_KUSAMA_LAPLACE_SPECTRUM__

#include "gaudi/kusama/spectrum.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace gaudi {
namespace kusama {

// `real` / `index_t` come from `spectrum.hpp`.

/// Build \f$L + \varepsilon I\f$ (sparse).
inline Eigen::SparseMatrix<real>
regularize_stiffness(const Eigen::SparseMatrix<real> &L, real epsilon) {
  Eigen::SparseMatrix<real> A = L;
  const index_t n = static_cast<index_t>(A.rows());
  for (index_t i = 0; i < n; ++i)
    A.coeffRef(i, i) += epsilon;
  A.makeCompressed();
  return A;
}

/// Spectra `SparseSymMatProd` multiplies using `selfadjointView<Eigen::Lower>()` only.
/// Cotan assembly may store asymmetric off-diagonals; symmetrize so matvecs match the
/// intended weak Laplacian.
inline Eigen::SparseMatrix<real>
symmetrize_sparse(const Eigen::SparseMatrix<real> &A) {
  Eigen::SparseMatrix<real> At = A.transpose();
  Eigen::SparseMatrix<real> S = A + At;
  S *= 0.5;
  S.makeCompressed();
  return S;
}

/// Log eigenvalue residual and coefficient variation (stderr). Use to verify solves.
/// Loose bound: \f$\|A\|_2 \le \|A\|_\infty\f$ (max absolute row sum) for symmetric \f$A\f$.
inline real sparse_sym_max_row_sum_abs(const Eigen::SparseMatrix<real> &A) {
  const index_t n = static_cast<index_t>(A.rows());
  if (n <= 0 || A.cols() != n)
    return real(0);
  Eigen::VectorXd row_sum = Eigen::VectorXd::Zero(n);
  for (index_t j = 0; j < A.outerSize(); ++j) {
    for (Eigen::SparseMatrix<real>::InnerIterator it(A, j); it; ++it)
      row_sum(it.row()) += std::abs(it.value());
  }
  return row_sum.maxCoeff();
}

inline void log_laplace_eigen_stats(const Eigen::SparseMatrix<real> &L, index_t mode_j,
                                    real lambda, const Eigen::VectorXd &x) {
  const index_t n = static_cast<index_t>(x.size());
  if (n == 0 || L.rows() != n) {
    std::cerr << "[laplace_spectrum] mode " << mode_j << " bad sizes L=" << L.rows()
              << " x=" << n << "\n";
    return;
  }
  Eigen::VectorXd r = L * x - lambda * x;
  const real res_inf = r.cwiseAbs().maxCoeff();
  const real x_inf = std::max(x.cwiseAbs().maxCoeff(), real(1e-30));
  const real mn = x.minCoeff(), mx = x.maxCoeff();
  const real mean = x.mean();
  const real var =
      (x.array() - mean).square().sum() / static_cast<real>(std::max(index_t(1), n - 1));
  const real sigma = std::sqrt(std::max(real(0), var));
  const real l2 = x.norm();
  const real span = mx - mn;
  const real rms_uniform = real(1) / std::sqrt(static_cast<real>(n));
  std::cerr << "[laplace_spectrum] mode " << mode_j << "  lambda=" << lambda
            << "  |phi|_2=" << l2 << "  phi_span=" << span << "  phi_sigma=" << sigma
            << "  ~1/sqrt(n)=" << rms_uniform << "  phi[" << mn << "," << mx << "]"
            << "  rel_res_inf=" << (res_inf / x_inf) << "\n";
}

/// Compute \f$k\f$ eigenpairs of \f$L_{\mathrm{reg}}\f$ with **largest magnitude**
/// eigenvalues. For the repo cotan stiffness (negative semi-definite before shift),
/// those correspond to the most oscillatory modes in the batch.
///
/// On success, `evals` has size `k` (descending magnitude order as returned by
/// Spectra), `evecs` is `n × k` with columns as eigenvectors.
///
/// @param ncv  Krylov size; if `<= 0`, uses `min(n, max(2*k+1, k+8))`.
/// @param warm_start  Optional Lanczos start (length `n`); see `spectrum.hpp`.
/// @return     `true` if `CompInfo::Successful`.
inline bool laplace_eigs_largest_magnitude(const Eigen::SparseMatrix<real> &L_reg,
                                           index_t k, index_t ncv,
                                           Eigen::VectorXd &evals,
                                           Eigen::MatrixXd &evecs,
                                           const Eigen::VectorXd *warm_start = nullptr) {
  return sparse_sym_eigs_largest_magnitude(L_reg, k, ncv, evals, evecs, warm_start);
}

/// Low-frequency batch: \f$k\f$ eigenpairs with **smallest algebraic** eigenvalues of
/// \f$A = -L_{\mathrm{sym}} + \varepsilon I\f$ (PSD when \f$L\f$ is NS-definite cotan).
/// Eigenvectors match those of \f$L_{\mathrm{sym}}\f$; eigenvalues are \f$\varepsilon-\mu\f$
/// for Laplacian eigenvalues \f$\mu\f$.
///
/// Same `evals` / `evecs` layout as `laplace_eigs_largest_magnitude`.
inline bool laplace_eigs_low_frequency(const Eigen::SparseMatrix<real> &A_reg, index_t k,
                                       index_t ncv, Eigen::VectorXd &evals,
                                       Eigen::MatrixXd &evecs,
                                       const Eigen::VectorXd *warm_start = nullptr) {
  return sparse_sym_eigs_smallest_algebraic(A_reg, k, ncv, evals, evecs, warm_start);
}

/// \f$k\f$ eigenpairs whose eigenvalues of \f$A\f$ lie **nearest** to \f$\sigma\f$
/// (shift–invert + `LargestMagn` on \f$(A-\sigma I)^{-1}\f$). \f$A\f$ must be symmetric.
/// Fails if \f$A-\sigma I\f$ is singular or factorization breaks.
///
/// @param ncv  if `<= 0`, uses `min(n, max(3*k+2, 2*k+12))` (shift mode wants larger Krylov).
/// @param warm_start  Optional Lanczos start (length `n`).
inline bool laplace_eigs_shift_invert_nearest(const Eigen::SparseMatrix<real> &A_reg,
                                              real sigma, index_t k, index_t ncv,
                                              Eigen::VectorXd &evals,
                                              Eigen::MatrixXd &evecs,
                                              const Eigen::VectorXd *warm_start = nullptr) {
  return sparse_sym_eigs_shift_invert_nearest(A_reg, sigma, k, ncv, evals, evecs,
                                              warm_start);
}

/// Back-compat alias for `laplace_eigs_low_frequency`.
inline bool laplace_eigs_smallest_magnitude(const Eigen::SparseMatrix<real> &L_reg,
                                            index_t k, index_t ncv,
                                            Eigen::VectorXd &evals,
                                            Eigen::MatrixXd &evecs,
                                            const Eigen::VectorXd *warm_start = nullptr) {
  return laplace_eigs_low_frequency(L_reg, k, ncv, evals, evecs, warm_start);
}

} // namespace kusama
} // namespace gaudi

#endif
