#ifndef GAUDI_BONTECOU_VECTOR_DIRICHLET_SPECTRUM_HPP
#define GAUDI_BONTECOU_VECTOR_DIRICHLET_SPECTRUM_HPP

/// Partial spectrum of the Stein et al. **vector Dirichlet** stiffness
/// (`build_vector_dirichlet_energy`, \f$2|E|\times 2|E|\f$).
///
/// The assembled matrix is symmetric positive semi-definite; we add \f$\varepsilon I\f$
/// so smallest-algebraic / largest-magnitude Spectra modes are well-defined.

#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/bontecou/laplace_spectrum.hpp"
#include "gaudi/bontecou/spectrum.hpp"
#include "gaudi/bontecou/vector_dirichlet.hpp"
#include "gaudi/common.h"

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace gaudi {
namespace bontecou {

/// Which batch of eigenpairs to request (same regimes as vertex Laplace helpers).
enum class vector_dirichlet_spectrum_band {
  smallest_algebraic, ///< `sparse_sym_eigs_smallest_algebraic` on \f$L+\varepsilon I\f$.
  largest_magnitude,  ///< `sparse_sym_eigs_largest_magnitude` on \f$L+\varepsilon I\f$.
  shift_invert_nearest ///< `sparse_sym_eigs_shift_invert_nearest` on \f$L+\varepsilon I\f$.
};

/// Build \f$L+\varepsilon I\f$ from [`build_vector_dirichlet_energy`](vector_dirichlet.hpp),
/// then compute `k` eigenpairs. Sets `num_edge_dofs` to \f$n_E\f$ (half the matrix size).
///
/// @return `true` if the Spectra solve succeeded.
inline bool vector_dirichlet_partial_spectrum(
    asawa::shell::shell &M, const std::vector<vec3> &x, double epsilon_shift,
    vector_dirichlet_spectrum_band band, double shift_sigma, index_t k, index_t ncv,
    Eigen::VectorXd &evals, Eigen::MatrixXd &evecs, int *num_edge_dofs = nullptr,
    const Eigen::VectorXd *warm_start = nullptr) {
  int nE = 0;
  Eigen::SparseMatrix<real> L = build_vector_dirichlet_energy(M, x, &nE);
  if (num_edge_dofs)
    *num_edge_dofs = nE;
  const int n = 2 * nE;
  if (n < 2 || k <= 0)
    return false;
  Eigen::SparseMatrix<real> A = regularize_stiffness(L, static_cast<real>(epsilon_shift));
  switch (band) {
  case vector_dirichlet_spectrum_band::smallest_algebraic:
    return sparse_sym_eigs_smallest_algebraic(A, k, ncv, evals, evecs, warm_start);
  case vector_dirichlet_spectrum_band::largest_magnitude:
    return sparse_sym_eigs_largest_magnitude(A, k, ncv, evals, evecs, warm_start);
  case vector_dirichlet_spectrum_band::shift_invert_nearest:
    return sparse_sym_eigs_shift_invert_nearest(
        A, static_cast<real>(shift_sigma), k, ncv, evals, evecs, warm_start);
  }
  return false;
}

} // namespace bontecou
} // namespace gaudi

#endif
