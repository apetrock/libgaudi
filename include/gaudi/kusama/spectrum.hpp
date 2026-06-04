/*
 * Generic sparse **symmetric** eigenproblems via Spectra.
 *
 * Callers supply a symmetric `Eigen::SparseMatrix<real>` (PSD / shifted as needed).
 * No mesh or Laplacian assumptions — see `laplace_spectrum.hpp` for cotan-specific docs.
 */

#ifndef GAUDI_KUSAMA_SPECTRUM_HPP
#define GAUDI_KUSAMA_SPECTRUM_HPP

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/MatOp/SparseSymShiftSolve.h>
#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/SymEigsSolver.h>

namespace gaudi {
namespace kusama {

using real = double;
using index_t = int;

/// Smallest **algebraic** eigenvalues of symmetric \f$A\f$.
/// @param warm_start  If non-null and `warm_start->size() == A.rows()`, seeds Lanczos;
///                    otherwise Spectra uses a random start vector.
inline bool sparse_sym_eigs_smallest_algebraic(const Eigen::SparseMatrix<real> &A,
                                               index_t k, index_t ncv,
                                               Eigen::VectorXd &evals,
                                               Eigen::MatrixXd &evecs,
                                               const Eigen::VectorXd *warm_start = nullptr) {
  const index_t n = static_cast<index_t>(A.rows());
  if (n < 2 || k <= 0 || A.cols() != n)
    return false;
  k = std::min(k, n - 1);
  if (ncv <= 0)
    ncv = std::min(n, std::max(2 * k + 1, k + 8));
  ncv = std::min(std::max(ncv, k + 1), n);

  Spectra::SparseSymMatProd<real> op(A);
  Spectra::SymEigsSolver<Spectra::SparseSymMatProd<real>> eigs(op, k, ncv);
  if (warm_start && warm_start->size() == n)
    eigs.init(warm_start->data());
  else
    eigs.init();
  eigs.compute(Spectra::SortRule::SmallestAlge);
  if (eigs.info() != Spectra::CompInfo::Successful)
    return false;
  evals = eigs.eigenvalues();
  evecs = eigs.eigenvectors();
  return true;
}

/// Eigenvalues with **largest magnitude** (Spectra `LargestMagn` on \f$A\f$).
inline bool sparse_sym_eigs_largest_magnitude(const Eigen::SparseMatrix<real> &A,
                                              index_t k, index_t ncv,
                                              Eigen::VectorXd &evals,
                                              Eigen::MatrixXd &evecs,
                                              const Eigen::VectorXd *warm_start = nullptr) {
  const index_t n = static_cast<index_t>(A.rows());
  if (n < 2 || k <= 0 || A.cols() != n)
    return false;
  k = std::min(k, n - 1);
  if (ncv <= 0)
    ncv = std::min(n, std::max(2 * k + 1, k + 8));
  ncv = std::min(std::max(ncv, k + 1), n);

  Spectra::SparseSymMatProd<real> op(A);
  Spectra::SymEigsSolver<Spectra::SparseSymMatProd<real>> eigs(op, k, ncv);
  if (warm_start && warm_start->size() == n)
    eigs.init(warm_start->data());
  else
    eigs.init();
  eigs.compute(Spectra::SortRule::LargestMagn);
  if (eigs.info() != Spectra::CompInfo::Successful)
    return false;
  evals = eigs.eigenvalues();
  evecs = eigs.eigenvectors();
  return true;
}

/// Shift–invert: eigenvalues of \f$A\f$ **nearest** to \f$\sigma\f$.
inline bool sparse_sym_eigs_shift_invert_nearest(const Eigen::SparseMatrix<real> &A,
                                                 real sigma, index_t k, index_t ncv,
                                                 Eigen::VectorXd &evals,
                                                 Eigen::MatrixXd &evecs,
                                                 const Eigen::VectorXd *warm_start = nullptr) {
  const index_t n = static_cast<index_t>(A.rows());
  if (n < 2 || k <= 0 || A.cols() != n)
    return false;
  k = std::min(k, n - 1);
  if (ncv <= 0)
    ncv = std::min(n, std::max(3 * k + 2, 2 * k + 12));
  ncv = std::min(std::max(ncv, k + 1), n);

  try {
    Spectra::SparseSymShiftSolve<real> op(A);
    Spectra::SymEigsShiftSolver<Spectra::SparseSymShiftSolve<real>> eigs(op, k, ncv,
                                                                         sigma);
    if (warm_start && warm_start->size() == n)
      eigs.init(warm_start->data());
    else
      eigs.init();
    eigs.compute(Spectra::SortRule::LargestMagn);
    if (eigs.info() != Spectra::CompInfo::Successful)
      return false;
    evals = eigs.eigenvalues();
    evecs = eigs.eigenvectors();
    return true;
  } catch (const std::invalid_argument &e) {
    std::cerr << "[spectrum] shift-invert: " << e.what() << "\n";
    return false;
  } catch (const std::exception &e) {
    std::cerr << "[spectrum] shift-invert: " << e.what() << "\n";
    return false;
  }
}

} // namespace kusama
} // namespace gaudi

#endif
