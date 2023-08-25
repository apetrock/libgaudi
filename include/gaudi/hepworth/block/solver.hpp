
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
    A.setFromTriplets(triplets.begin(), triplets.end());
    matS AtA = A.transpose() * A;
    matS MAtA = M + AtA;

#if 1
    vecX z = vecX::Zero(MAtA.rows());
    if ((MAtA * z).hasNaN()) {
      std::cout << "MAtA has NaN" << std::endl;
      exit(0);
    }
#endif

    m_solver S(MAtA);
    if (!S.success()) {
      Eigen::SparseMatrix<double> I(MAtA.rows(), MAtA.cols());
      MAtA += 1e-6 * I;
      S.compute(MAtA);
    }
    if (!S.success()) {
      std::cout << "Solve failed" << std::endl;
      return;
    }

    std::cout << "A   sum: " << A.sum() << std::endl;
    std::cout << "AtA sum: " << AtA.sum() << std::endl;
    std::cout << "M sum: " << M.sum() << std::endl;

    vecX p = vecX::Zero(id0);

    for (int k = 0; k < ITS; k++) {
      p.setZero();
      int ii = 0;
      // #pragma omp parallel for
      for (int i = 0; i < _constraints.size(); i++) {
        auto &constraint = _constraints[i];
        constraint->project(q, p);
#if 1
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
      if (q.hasNaN()) {
        std::cout << "q has NaN" << std::endl;
      }

      if (p.hasNaN()) {
        std::cout << "p has NaN" << std::endl;
      }

      if (k % 10 == 0)
        std::cout << "k: " << k << " -pnorm: " << p.norm() << std::endl;

      vecX b = M * s + A.transpose() * p;

      q = S.solve(b);

      if (k % 10 == 0)
        std::cout << "k: " << k << " -norms: "
                  << " p-" << p.norm() //
                  << " q-" << q.norm() //
                  << " s-" << s.norm() << std::endl;

      if (q.hasNaN()) {
        std::cout << "q has NaN" << std::endl;
        exit(0);
      }
      // q = qi + dq.min(bnd).max(-bnd);
    }

    for (auto &block : blocks) {
      block->map_from_x(q, h, damping);
    }
  }

  void
  set_constraints(const std::vector<projection_constraint::ptr> &constraints) {
    _constraints = constraints;
  }

  matS __M;
  std::vector<projection_constraint::ptr> _constraints;
}; // namespace block
} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif