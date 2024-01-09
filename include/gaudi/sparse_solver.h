

#ifndef __LIBGAUDI_SPARSE_SOLVER__
#define __LIBGAUDI_SPARSE_SOLVER__

#ifdef CHOLMOD_H
#define USE_CHOLMOD 1
#else
#define USE_CHOLMOD 1 // I guess I need to force it..
#endif

#if USE_CHOLMOD
#include <Eigen/CholmodSupport>
#endif

#include <Eigen/Sparse>

typedef double real;
typedef Eigen::SparseMatrix<double> matS;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vecX;
namespace gaudi {

class m_solver {
public:
  m_solver(matS &A) { this->compute(A); }

  void compute(matS &A) {
    _decomposed = false;
    __solver.compute(A);

#if USE_CHOLMOD
    std::cout << "solving with cholmod" << std::endl;
#endif

    if (__solver.info() == Eigen::Success) {
      _decomposed = true;
    } else {
      std::cout << ".....decomposition error! " << std::endl;
    }
  }
  vecX solve(vecX &b) {
    vecX x = __solver.solve(b);
    if (__solver.info() != Eigen::Success) {
      // solving failed
      std::cout << ".....solve error! " << std::endl;
    }
    return x;
  }

  bool success() { return _decomposed; }
  // macro to see if cholamod header is available

#if USE_CHOLMOD
  Eigen::CholmodSupernodalLLT<matS> __solver;
#else
  Eigen::SimplicialLDLT<matS> __solver;
#endif

  bool _decomposed = false;
};

vecX solve(matS &A, vecX &b) {

  // Eigen::ConjugateGradient<matS, Eigen::Upper> solver;
  Eigen::SimplicialLDLT<matS> solver;
  solver;

  solver.compute(A);

  if (solver.info() != Eigen::Success) {
    // decomposition failed
    std::cout << ".....decomposition error! " << std::endl;
  }
  vecX x = solver.solve(b);
  if (solver.info() != Eigen::Success) {
    // solving failed
    std::cout << ".....solve error! " << std::endl;
  }

  return x;
}

} // namespace gaudi
#endif