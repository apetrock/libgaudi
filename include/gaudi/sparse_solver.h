#ifndef __LIBGAUDI_SPARSE_SOLVER__
#define __LIBGAUDI_SPARSE_SOLVER__

#ifdef CHOLMOD_H
#define USE_CHOLMOD 1
#else
#define USE_CHOLMOD 0 // I guess I need to force it..
#endif

#if USE_CHOLMOD
#include <Eigen/CholmodSupport>
#endif

#include <Eigen/Sparse>
#include "gaudi/logger.hpp"

typedef double real;
typedef Eigen::SparseMatrix<double> matS;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vecX;
namespace gaudi {

class m_solver {
public:
  m_solver(matS &A) { this->compute(A); }

  void compute(matS &A) {
    _decomposed = false;
    _use_iterative = false;
    
    // Store matrix for potential iterative solver use
    _matrix = A;
    
    // Try direct factorization first
    __solver.compute(A);

#if USE_CHOLMOD
    gaudi::logger::info << "solving with cholmod" << std::endl;
#endif

    if (__solver.info() == Eigen::Success) {
      _decomposed = true;
    } else {
      gaudi::logger::error << ".....decomposition error! " << std::endl;
      // Fallback to iterative solver for WebAssembly
      _use_iterative = true;
      _decomposed = true; // Mark as successful so we can use iterative solver
    }
  }
  
  vecX solve(vecX &b) {
    if (_use_iterative) {
      // Use iterative solver as fallback
      Eigen::ConjugateGradient<matS, Eigen::Lower | Eigen::Upper> cg;
      cg.compute(_matrix);
      cg.setMaxIterations(1000);
      cg.setTolerance(1e-6);
      vecX x = cg.solve(b);
      if (cg.info() != Eigen::Success) {
        gaudi::logger::error << ".....iterative solve error! " << std::endl;
      }
      return x;
    } else {
      vecX x = __solver.solve(b);
      if (__solver.info() != Eigen::Success) {
        // solving failed
        gaudi::logger::error << ".....solve error! " << std::endl;
      }
      return x;
    }
  }

  bool success() { return _decomposed; }
  // macro to see if cholamod header is available

#if USE_CHOLMOD
  Eigen::CholmodSupernodalLLT<matS> __solver;
#else
  Eigen::SimplicialLDLT<matS, Eigen::Lower, Eigen::NaturalOrdering<int>> __solver;
#endif

  bool _decomposed = false;
  bool _use_iterative = false;
  matS _matrix; // Store matrix for iterative solver
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