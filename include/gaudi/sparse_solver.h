#ifndef __LIBGAUDI_SPARSE_SOLVER__
#define __LIBGAUDI_SPARSE_SOLVER__

#include <algorithm>

#ifdef CHOLMOD_H
#define USE_CHOLMOD 1
#else
#define USE_CHOLMOD 0 // I guess I need to force it..
#endif

#if USE_CHOLMOD
#include <Eigen/CholmodSupport>
#endif

#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <iostream>
#include "gaudi/duchamp/demo_log.hpp"
#include "gaudi/logger.hpp"

typedef double real;
typedef Eigen::SparseMatrix<double> matS;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vecX;
namespace gaudi {

namespace {

// Direct sparse LDLT fill-in / memory grows superlinearly; above this size use
// CG. The coupled curvature systems (fourth-order in positions, a massless
// store) defeat diagonal-preconditioned CG, so stay direct well past 4M
// nonzeros.
inline bool prefer_iterative_solve(const matS &A) {
  const Eigen::Index n = A.rows();
  if (n <= 0)
    return false;
  const Eigen::Index nnz = A.nonZeros();
  return n > 3200000 || nnz > 40000000;
}

using CGSolver =
    Eigen::ConjugateGradient<matS, Eigen::Lower | Eigen::Upper,
                             Eigen::DiagonalPreconditioner<double>>;

inline void setup_cg(CGSolver &cg, const matS &A) {
  cg.compute(A);
  const int n = static_cast<int>(A.rows());
  // Large normal equations need many iterations; cap to avoid pathological runtimes.
  cg.setMaxIterations(std::max(500, std::min(50000, std::max(n / 4, 2000))));
  cg.setTolerance(1e-5);
}

} // namespace

class m_solver {
public:
  m_solver(matS &A) { this->compute(A); }

  void compute(matS &A) {
    _decomposed = false;
    _use_iterative = false;

    gaudi::duchamp::demo_log() << "[m_solver] compute: " << A.rows() << "x" << A.cols()
              << " nnz=" << A.nonZeros() << std::endl;
    gaudi::duchamp::demo_log().flush();

    gaudi::duchamp::demo_log() << "[m_solver] copying matrix for storage..." << std::endl;
    gaudi::duchamp::demo_log().flush();
    _matrix = A;
    gaudi::duchamp::demo_log() << "[m_solver] copy done" << std::endl;
    gaudi::duchamp::demo_log().flush();

    if (prefer_iterative_solve(A)) {
      gaudi::duchamp::demo_log() << "[m_solver] large system: using CG+Diagonal (skip direct LDLT)"
                << std::endl;
      gaudi::duchamp::demo_log().flush();
      gaudi::duchamp::demo_log() << "[m_solver] setup_cg: compute()..." << std::endl;
      gaudi::duchamp::demo_log().flush();
      setup_cg(_cg_iter, _matrix);
      gaudi::duchamp::demo_log() << "[m_solver] setup_cg: done info=" << _cg_iter.info()
                << std::endl;
      gaudi::duchamp::demo_log().flush();
      _use_iterative = true;
      _decomposed = (_cg_iter.info() == Eigen::Success);
      if (!_decomposed)
        gaudi::logger::error << "[m_solver] CG preconditioner setup failed"
                             << std::endl;
      gaudi::duchamp::demo_log() << "[m_solver] compute() returning (iterative path)" << std::endl;
      gaudi::duchamp::demo_log().flush();
      return;
    }

    gaudi::duchamp::demo_log() << "[m_solver] matrix copied, factorizing (LDLT)..." << std::endl;
    gaudi::duchamp::demo_log().flush();

#if USE_CHOLMOD
    gaudi::logger::info << "solving with cholmod" << std::endl;
#endif

    gaudi::duchamp::demo_log() << "[m_solver] calling SimplicialLDLT::compute (this can take RAM)..."
              << std::endl;
    gaudi::duchamp::demo_log().flush();
    __solver.compute(A);
    gaudi::duchamp::demo_log() << "[m_solver] factorization complete info=" << __solver.info();
#if !USE_CHOLMOD
    if (__solver.info() == Eigen::Success) {
      const double l_nnz = static_cast<double>(__solver.matrixL().nestedExpression().nonZeros());
      gaudi::duchamp::demo_log() << " L_nnz=" << l_nnz << " (" << l_nnz / std::max<double>(A.nonZeros(), 1)
                                 << "x A, ~" << l_nnz * (sizeof(real) + sizeof(int)) / 1048576.0
                                 << " MB)";
    }
#endif
    gaudi::duchamp::demo_log() << std::endl;
    gaudi::duchamp::demo_log().flush();

    if (__solver.info() == Eigen::Success) {
      _decomposed = true;
    } else {
      gaudi::logger::error << ".....decomposition error! trying CG fallback"
                           << std::endl;
      setup_cg(_cg_iter, _matrix);
      _use_iterative = true;
      _decomposed = (_cg_iter.info() == Eigen::Success);
    }
  }

  vecX solve(vecX &b) {
    ++_solve_call;
    const int sc = _solve_call;
    const bool trace = (sc <= 3) || (sc % 10 == 0);
    const real b_norm = b.norm();
    if (trace) {
      gaudi::duchamp::demo_log() << "[m_solver] solve #" << sc << " begin b.rows=" << b.size()
                << " b_norm=" << b_norm << std::endl;
      gaudi::duchamp::demo_log().flush();
    }
    if (_use_iterative) {
      vecX x = _cg_iter.solve(b);
      if (trace) {
        const real x_norm = x.norm();
        const real x_over_b = (b_norm > 0.0) ? x_norm / b_norm : 0.0;
        gaudi::duchamp::demo_log() << "[m_solver] solve #" << sc << " end (CG) info="
                  << _cg_iter.info() << " x_norm=" << x_norm
                  << " x/b=" << x_over_b
                  << " iters=" << _cg_iter.iterations() << std::endl;
        gaudi::duchamp::demo_log().flush();
      }
      if (_cg_iter.info() == Eigen::Success)
        return x;
      // An unconverged CG step is garbage; factor directly instead and stay
      // direct for this matrix.
      gaudi::logger::error << ".....iterative solve did not converge; factoring directly"
                           << std::endl;
      __solver.compute(_matrix);
      if (__solver.info() != Eigen::Success) {
        gaudi::logger::error << ".....direct fallback failed; returning the CG result"
                             << std::endl;
        return x;
      }
      _use_iterative = false;
    }
    vecX x = __solver.solve(b);
    if (trace) {
      const real x_norm = x.norm();
      const real x_over_b = (b_norm > 0.0) ? x_norm / b_norm : 0.0;
      gaudi::duchamp::demo_log() << "[m_solver] solve #" << sc << " end (LDLT) info="
                << __solver.info() << " x_norm=" << x_norm
                << " x/b=" << x_over_b << std::endl;
      gaudi::duchamp::demo_log().flush();
    }
    if (__solver.info() != Eigen::Success) {
      gaudi::logger::error << ".....solve error! " << std::endl;
    }
    return x;
  }

  bool success() { return _decomposed; }
  /// The direct LDLT factor `__solver` is valid (not the CG path).
  bool direct() const { return _decomposed && !_use_iterative; }
  // macro to see if cholamod header is available

#if USE_CHOLMOD
  Eigen::CholmodSupernodalLLT<matS> __solver;
#else
  // AMDOrdering cuts fill-in vs NaturalOrdering on mesh Laplacians / AᵀA.
  Eigen::SimplicialLDLT<matS, Eigen::Lower, Eigen::AMDOrdering<int>> __solver;
#endif

  CGSolver _cg_iter;

  bool _decomposed = false;
  bool _use_iterative = false;
  int _solve_call = 0;
  matS _matrix; // copy for iterative path & fallback
};

vecX solve(matS &A, vecX &b) {

  // Eigen::ConjugateGradient<matS, Eigen::Upper> solver;
  Eigen::SimplicialLDLT<matS> solver;

  solver.compute(A);

  if (solver.info() != Eigen::Success) {
    // decomposition failed
    gaudi::logger::error << ".....decomposition error! " << std::endl;
  }
  vecX x = solver.solve(b);
  if (solver.info() != Eigen::Success) {
    // solving failed
    gaudi::logger::error << ".....solve error! " << std::endl;
  }

  return x;
}

} // namespace gaudi
#endif
