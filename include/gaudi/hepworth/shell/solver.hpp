
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

#include "constraints.hpp"
#include "gaudi/common.h"

///////////////////////////////////////////////////////
// solver
///////////////////////////////////////////////////////

namespace gaudi {
namespace hepworth {
namespace shell {
class projection_solver {
public:
  projection_solver() {}

  void set_mass(const std::vector<vec3> &m) {

    vecX M = to(m);

    std::vector<trip> triplets;
    for (int i = 0; i < M.rows(); i++) {
      triplets.push_back(trip(i, i, M[i]));
      // triplets.push_back(trip(i, i, 0.15));
    }
    __M = matS(M.size(), M.size());
    __M.setFromTriplets(triplets.begin(), triplets.end());
  }

  void update_velocity(const real &h, std::vector<vec3> &x_,
                       const std::vector<vec3> &v_,
                       const std::vector<vec3> &f_) {
    for (size_t i = 0; i < v_.size(); i++) {
      x_[i] += h * v_[i] + h * h * f_[i];
    }
  }

  void step(std::vector<vec3> &x_,
            std::vector<vec3> &v_,       // velocity;
            const std::vector<vec3> &f_, // velocity;
            const real &h = 0.01) {

    int Ni = 3 * x_.size();

    vecX q0 = to(x_);
    vecX v = to(v_);
    vecX f = to(f_);

    vecX q = q0 + h * v + h * h * f;

    std::vector<trip> triplets;
    index_t id0 = 0;
    for (auto &constraint : _constraints) {
      constraint->fill_A(id0, triplets);
    }

    // std::cout << A.rows() << " " << A.cols() << " " << triplets.size()
    //          << std::endl;
    matS A(id0, Ni);
    matS &M = __M;
    M *= 1.0 / h / h;

    A.setFromTriplets(triplets.begin(), triplets.end());
    matS AtA = A.transpose() * A;
    matS MAtA = M + AtA;
    m_solver S(MAtA);
    std::cout << "A   sum: " << A.sum() << std::endl;
    std::cout << "AtA sum: " << AtA.sum() << std::endl;
    std::cout << "M sum: " << M.sum() << std::endl;

    // vecX x0 = x;
    vecX s = q;
    vecX p = vecX::Zero(id0);

    for (int k = 0; k < 10; k++) {
      p.setZero();
      for (auto &constraint : _constraints) {
        constraint->project(q, p);
      }

      if (k % 1 == 0) {
        std::cout << "pnorm: " << p.norm() << std::endl;
        std::cout << "ATpnorm: " << (A.transpose() * p).norm() << std::endl;
      }

      vecX b = M * s + A.transpose() * p;

      q = S.solve(b);

      // q = qi + dq.min(bnd).max(-bnd);
    }

    std::cout << " x norm: " << (q - q0).norm() << std::endl;
    real damp = 0.05;
    v = (1.0 - damp) / h * (q - q0);
    from(v_, v);
    from(x_, q);
  }

  void
  set_constraints(const std::vector<projection_constraint::ptr> &constraints) {
    _constraints = constraints;
  }

  matS __M;
  std::vector<projection_constraint::ptr> _constraints;
};
} // namespace shell
} // namespace hepworth
} // namespace gaudi

#endif