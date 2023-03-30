
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
namespace rod {
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

  void set_mass(const std::vector<vec3> &M_, const std::vector<vec4> &J_) {

    vecX M = to(M_);
    vecX J = to(J_);

    vecX MJ = concat(M, J);
    std::vector<trip> triplets;
    for (int i = 0; i < MJ.rows(); i++) {
      triplets.push_back(trip(i, i, MJ[i]));
    }
    __M = matS(MJ.size(), MJ.size());
    __M.setFromTriplets(triplets.begin(), triplets.end());
  }

  void update_rotational_inertia(const real &h, std::vector<quat> &u_,
                                 const std::vector<quat> &o_ // omega
  ) {
    for (size_t i = 0; i < o_.size(); i++) {
      quat O = o_[i];
      quat u = u_[i];

      // real J = J;
      quat sO = O; // torque terms + h / J
      quat su = u;
      su.coeffs() += 0.5 * h * (u * sO).coeffs();
      u_[i] = su;
      u_[i].normalize();
    }
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
            std::vector<quat> &u_,       // quats,
            std::vector<quat> &o_,       // omega
            const real &h = 0.01) {

    int Ni = x_.size();
    int Nm = 3 * Ni + 4 * Ni;

    vecX x0 = to(x_);
    vecX u0 = to(u_);

    update_rotational_inertia(h, u_, o_);
    update_velocity(h, x_, v_, f_);

    vecX x = to(x_);
    vecX v = to(v_);
    vecX u = to(u_);
    vecX o = to(o_);

    std::vector<trip> triplets;
    index_t id0 = 0;
    for (auto &constraint : _constraints) {
      constraint->fill_A(id0, triplets);
    }

    matS A(id0, Nm);
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

    vecX q = concat(x0, u0);
    vecX s = concat(x, u);

    vecX p = vecX::Zero(id0);

    for (int k = 0; k < 50; k++) {
      p.setZero();
      for (auto &constraint : _constraints) {
        constraint->project(q, p);
      }

      if (k % 10 == 0)
        std::cout << "pnorm: " << p.norm() << std::endl;

      vecX b = M * s + A.transpose() * p;
      // std::cout << p << std::endl;
      q = S.solve(b);
      // q = qi + dq.min(bnd).max(-bnd);
    }

    split(q, x, u);
    std::cout << " x norm: " << (x - x0).norm() << std::endl;
    std::cout << " u norm: " << (u - u0).norm() << std::endl;

    real damp = 0.5;
    v = (1.0 - damp) / h * (x - x0);
    from(v_, v);
    from(x_, x);
    from(u_, u);
    std::vector<quat> u0_(u_);
    from(u0_, u0);
    for (int i = 0; i < u_.size(); i++) {
      o_[0] = u0_[i].conjugate() * u_[i];
      o_[0].coeffs() *= 2.0 * (1.0 - damp) / h;
    }
  }

  void
  set_constraints(const std::vector<projection_constraint::ptr> &constraints) {
    _constraints = constraints;
  }

  matS __M;
  std::vector<projection_constraint::ptr> _constraints;
};
} // namespace rod
} // namespace hepworth
} // namespace gaudi

#endif