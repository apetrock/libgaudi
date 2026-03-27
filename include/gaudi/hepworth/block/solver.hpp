
#ifndef __HEP_ROD_SOLVER__
#define __HEP_ROD_SOLVER__

#include "gaudi/sparse_solver.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cxxabi.h>

#include <iostream>
#include <memory.h>
#include <numeric>
#include <ostream>
#include <stdio.h>
#include <type_traits>
#include <vector>

#include "block_constraint.hpp"
#include "gaudi/common.h"
#include "sim_block.hpp"

///////////////////////////////////////////////////////
// solver
///////////////////////////////////////////////////////

namespace gaudi {
namespace hepworth {
namespace block {
class projection_solver {
public:
  projection_solver() {}

  void set_mass(std::vector<sim_block::ptr> &blocks) {

    vecX M;
    for (auto &block : blocks) {
      block->map_mass(M);
    }
    std::vector<trip> triplets;
    for (int i = 0; i < M.rows(); i++) {
      real Mi = M[i] < 1e-8 ? 1e-8 : M[i];
      triplets.push_back(trip(i, i, M[i]));
    }
    __M = matS(M.size(), M.size());
    __M.setFromTriplets(triplets.begin(), triplets.end());
  }

  void step(std::vector<sim_block::ptr> &blocks, const real &h = 0.01,
            const real &damping = 0.5, const int &ITS = 50) {

    vecX q;
    vecX s;
    set_mass(blocks);

    for (auto &block : blocks) {
      block->map_to_x(q);
      block->integrate_inertia(h);
      block->map_to_x(s);
    }
    int Nm = q.size();

    std::vector<trip> triplets;
    index_t id0 = 0;

    for (int i = 0; i < _constraints.size(); i++) {
      auto &constraint = _constraints[i];
      constraint->fill_A(id0, triplets);
    }

    matS A(id0, Nm);
    matS &M = __M;
    M *= 1.0 / h / h;
    std::cout << "Nm: " << Nm << std::endl;
    std::cerr << "[projection_solver] setFromTriplets..." << std::endl;
    std::cerr.flush();
    A.setFromTriplets(triplets.begin(), triplets.end());
    std::cerr << "[projection_solver] A nnz=" << A.nonZeros() << std::endl;
    std::cerr.flush();
    std::cerr << "[projection_solver] computing AtA = A^T * A (heavy)..."
              << std::endl;
    std::cerr.flush();
    matS AtA = A.transpose() * A;
    std::cerr << "[projection_solver] AtA nnz=" << AtA.nonZeros() << std::endl;
    std::cerr.flush();
    std::cerr << "[projection_solver] MAtA = M + AtA..." << std::endl;
    std::cerr.flush();
    matS MAtA = M + AtA;
    std::cerr << "[projection_solver] MAtA nnz=" << MAtA.nonZeros() << std::endl;
    std::cerr.flush();

#if 0
    if (!MAtA.isApprox(MAtA.transpose())) {
      throw std::runtime_error("Possibly non semi-positive definitie matrix!");
    }
    vecX z = vecX::Zero(MAtA.rows());
    if ((MAtA * z).hasNaN()) {
      std::cout << "MAtA has NaN" << std::endl;
      exit(0);
    }
    for (int i = 0; i < MAtA.rows(); i++) {
      if (MAtA.coeff(i, i) < 1e-9) {
        std::cout << "MAtA has zero diagonal" << std::endl;
        exit(0);
      }
    }
#endif

    std::cout << "A   sum: " << A.sum() << std::endl;
    std::cout << "AtA sum: " << AtA.sum() << std::endl;
    std::cout << "M sum: " << M.sum() << std::endl;
    std::cout << "MAtA sum: " << MAtA.sum() << std::endl;
    std::cout << "A: " << A.rows() << "x" << A.cols()
              << " nnz=" << A.nonZeros() << std::endl;
    std::cout << "MAtA: " << MAtA.rows() << "x" << MAtA.cols()
              << " nnz=" << MAtA.nonZeros() << std::endl;
    std::cout.flush();
    std::cerr << "[projection_solver] constructing m_solver(MAtA)..."
              << std::endl;
    std::cerr.flush();

    m_solver S(MAtA);

    std::cerr << "[projection_solver] m_solver ready, success=" << S.success()
              << std::endl;
    std::cerr.flush();
    if (!S.success()) {
      std::cout << "Solve failed, attempting normalization" << std::endl;
      Eigen::SparseMatrix<double> I(MAtA.rows(), MAtA.cols());
      MAtA += 1e-6 * I;
      S.compute(MAtA);
    }
    if (!S.success()) {
      std::cout << "Solve failed" << std::endl;
      return;
    }

    vecX p = vecX::Zero(id0);
    vecX q0 = q;
    for (int k = 0; k < ITS; k++) {
      if (k == 0 || k % 10 == 0) {
        std::cerr << "[projection_solver] outer iter k=" << k << "/"
                  << (ITS - 1) << std::endl;
        std::cerr.flush();
      }
      p.setZero();
      int ii = 0;
//#pragma omp parallel for
      if (k == 0) {
        std::cerr << "[projection_solver] project: " << _constraints.size()
                  << " constraints..." << std::endl;
        std::cerr.flush();
      }
      for (int i = 0; i < _constraints.size(); i++) {
        auto &constraint = _constraints[i];
        constraint->project(q, p);
#if 0
      if (i < _constraints.size() - 1) {
          index_t id0 = _constraints[i]->_id0;
          index_t id1 = _constraints[i + 1]->_id0;
          vecX pi = p.segment(id0, id1 - id0);
          if (pi.hasNaN()) {
            std::cout << i << " " << _constraints.size() << std::endl;
            std::cout << "q has NaN: "
                      << "ii - " << id0 << " " << id1 - id0 << " " << p.size()
                      << " " << constraint->name() << std::endl;
            exit(0);
          }
        }
#endif
        ii++;
      }
      if (k == 0) {
        std::cerr << "[projection_solver] project done" << std::endl;
        std::cerr.flush();
      }
      if (q.hasNaN()) {
        std::cout << "q has NaN" << std::endl;
      }

      if (p.hasNaN()) {
        std::cout << "p has NaN" << std::endl;
      }

      if (k % 10 == 0)
        std::cout << "k: " << k << " -pnorm: " << p.norm() << std::endl;

      if (k == 0) {
        std::cerr << "[projection_solver] forming rhs b = M*s + A^T*p..."
                  << std::endl;
        std::cerr.flush();
      }
      vecX b = M * s + A.transpose() * p;
      if (k == 0) {
        std::cerr << "[projection_solver] calling S.solve (k=" << k << ")..."
                  << std::endl;
        std::cerr.flush();
      }

      q = S.solve(b);

      if (k == 0) {
        std::cerr << "[projection_solver] S.solve returned k=" << k << std::endl;
        std::cerr.flush();
      }

      real dq = (q - q0).norm();
      if (dq < 1e-8)
        break;
      q0 = q;

      if (k % 10 == 0)
        std::cout << "k: " << k << " -norms: "
                  << " q-q0: " << dq << std::endl;

      // q = qi + dq.min(bnd).max(-bnd);
    }

    std::cerr << "[projection_solver] outer loop finished, map_from_x..."
              << std::endl;
    std::cerr.flush();
    for (auto &block : blocks) {
      block->map_from_x(q, h, damping);
    }
    std::cerr << "[projection_solver] step() complete" << std::endl;
    std::cerr.flush();
  }

  void
  set_constraints(const std::vector<projection_constraint::ptr> &constraints) {
    _constraints = constraints;
  }

  matS __M;
  std::vector<projection_constraint::ptr> _constraints;
}; // class projection_solver
} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif