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
#include "gaudi/logger.hpp"

typedef double real;
typedef Eigen::SparseMatrix<double> matS;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vecX;
namespace gaudi {

namespace {

// Direct sparse LDLT fill-in / memory grows superlinearly; above this size use CG.
inline bool prefer_iterative_solve(const matS &A) {
  const Eigen::Index n = A.rows();
  if (n <= 0)
    return false;
  const Eigen::Index nnz = A.nonZeros();
  return n > 3200000 || nnz > 4000000;
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

    std::cerr << "[m_solver] compute: " << A.rows() << "x" << A.cols()
              << " nnz=" << A.nonZeros() << std::endl;
    std::cerr.flush();

    std::cerr << "[m_solver] copying matrix for storage..." << std::endl;
    std::cerr.flush();
    _matrix = A;
    std::cerr << "[m_solver] copy done" << std::endl;
    std::cerr.flush();

    if (prefer_iterative_solve(A)) {
      std::cerr << "[m_solver] large system: using CG+Diagonal (skip direct LDLT)"
                << std::endl;
      std::cerr.flush();
      std::cerr << "[m_solver] setup_cg: compute()..." << std::endl;
      std::cerr.flush();
      setup_cg(_cg_iter, _matrix);
      std::cerr << "[m_solver] setup_cg: done info=" << _cg_iter.info()
                << std::endl;
      std::cerr.flush();
      _use_iterative = true;
      _decomposed = (_cg_iter.info() == Eigen::Success);
      if (!_decomposed)
        gaudi::logger::error << "[m_solver] CG preconditioner setup failed"
                             << std::endl;
      std::cerr << "[m_solver] compute() returning (iterative path)" << std::endl;
      std::cerr.flush();
      return;
    }

    std::cerr << "[m_solver] matrix copied, factorizing (LDLT)..." << std::endl;
    std::cerr.flush();

#if USE_CHOLMOD
    gaudi::logger::info << "solving with cholmod" << std::endl;
#endif

    std::cerr << "[m_solver] calling SimplicialLDLT::compute (this can take RAM)..."
              << std::endl;
    std::cerr.flush();
    __solver.compute(A);
    std::cerr << "[m_solver] factorization complete info=" << __solver.info()
              << std::endl;
    std::cerr.flush();

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
    if (trace) {
      std::cerr << "[m_solver] solve #" << sc << " begin b.rows=" << b.size()
                << " b_norm=" << b.norm() << std::endl;
      std::cerr.flush();
    }
    if (_use_iterative) {
      vecX x = _cg_iter.solve(b);
      if (trace) {
        std::cerr << "[m_solver] solve #" << sc << " end (CG) info="
                  << _cg_iter.info() << " x_norm=" << x.norm()
                  << " iters=" << _cg_iter.iterations() << std::endl;
        std::cerr.flush();
      }
      if (_cg_iter.info() != Eigen::Success) {
        gaudi::logger::error << ".....iterative solve error! " << std::endl;
      }
      return x;
    }
    vecX x = __solver.solve(b);
    if (trace) {
      std::cerr << "[m_solver] solve #" << sc << " end (LDLT) info="
                << __solver.info() << " x_norm=" << x.norm() << std::endl;
      std::cerr.flush();
    }
    if (__solver.info() != Eigen::Success) {
      gaudi::logger::error << ".....solve error! " << std::endl;
    }
    return x;
  }

  bool success() { return _decomposed; }
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
