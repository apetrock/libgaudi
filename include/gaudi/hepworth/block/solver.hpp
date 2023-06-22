
#ifndef __HEP_ROD_SOLVER__
#define __HEP_ROD_SOLVER__

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cxxabi.h>
#include <execinfo.h>
#include <iostream>
#include <memory.h>
#include <numeric>
#include <ostream>
#include <stdio.h>
#include <type_traits>
#include <vector>
#include <zlib.h>

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

  void set_mass(std::vector<sim_block::ptr> &blocks) {

    vecX M;
    for (auto &block : blocks) {
      block->map_mass(M);
    }
    std::vector<trip> triplets;
    for (int i = 0; i < M.rows(); i++) {
      triplets.push_back(trip(i, i, M[i]));
    }
    __M = matS(M.size(), M.size());
    __M.setFromTriplets(triplets.begin(), triplets.end());
  }

  void step(std::vector<sim_block::ptr> &blocks, const real &h = 0.01,
            const int &ITS = 50) {

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
    for (auto &constraint : _constraints) {
      constraint->fill_A(id0, triplets);
    }

    matS A(id0, Nm);
    matS &M = __M;
    M *= 1.0 / h / h;
    std::cout << "Nm: " << Nm << std::endl;
    A.setFromTriplets(triplets.begin(), triplets.end());
    matS AtA = A.transpose() * A;
    matS MAtA = M + AtA;
    m_solver S(MAtA);
    std::cout << "A   sum: " << A.sum() << std::endl;
    std::cout << "AtA sum: " << AtA.sum() << std::endl;
    std::cout << "M sum: " << M.sum() << std::endl;

    vecX p = vecX::Zero(id0);

    for (int k = 0; k < ITS; k++) {
      p.setZero();
      for (auto &constraint : _constraints) {
        constraint->project(q, p);
      }

      if (k % 10 == 0)
        std::cout << "k: " << k << " -pnorm: " << p.norm() << std::endl;

      vecX b = M * s + A.transpose() * p;
      q = S.solve(b);

      // q = qi + dq.min(bnd).max(-bnd);
    }

    for (auto &block : blocks) {
      block->map_from_x(q, h, 0.5);
    }
  }

  void
  set_constraints(const std::vector<projection_constraint::ptr> &constraints) {
    _constraints = constraints;
  }

  matS __M;
  std::vector<projection_constraint::ptr> _constraints;
};
} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif