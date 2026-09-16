
#ifndef __HEP_ROD_SOLVER__
#define __HEP_ROD_SOLVER__

#include "gaudi/sparse_solver.h"

#include "projection_diagnostics.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cxxabi.h>
#include <functional>
#include <memory>
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

struct solver_context;

namespace detail {

inline size_t hash_triplets(const std::vector<trip> &triplets) {
  size_t h = triplets.size();
  for (const trip &t : triplets) {
    h ^= std::hash<int>()(static_cast<int>(t.row())) + 0x9e3779b9 + (h << 6) +
         (h >> 2);
    h ^= std::hash<int>()(static_cast<int>(t.col())) + 0x9e3779b9 + (h << 6) +
         (h >> 2);
    h ^= std::hash<double>()(t.value()) + 0x9e3779b9 + (h << 6) + (h >> 2);
  }
  return h;
}

inline real mass_diagonal_sum(const matS &M) {
  real sum = 0.0;
  for (int k = 0; k < M.outerSize(); ++k) {
    for (matS::InnerIterator it(M, k); it; ++it) {
      if (it.row() == it.col())
        sum += it.value();
    }
  }
  return sum;
}

struct linear_system_cache_key {
  int nm = 0;
  index_t a_rows = 0;
  size_t constraint_count = 0;
  size_t triplet_hash = 0;
  real h = 0.0;
  real mass_sum = 0.0;

  bool operator==(const linear_system_cache_key &o) const {
    return nm == o.nm && a_rows == o.a_rows &&
           constraint_count == o.constraint_count &&
           triplet_hash == o.triplet_hash && h == o.h && mass_sum == o.mass_sum;
  }
};

} // namespace detail

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
      triplets.push_back(trip(i, i, Mi));
    }
    __M = matS(M.size(), M.size());
    __M.setFromTriplets(triplets.begin(), triplets.end());
  }

  void step(std::vector<sim_block::ptr> &blocks, const real &h = 0.01,
            const real &damping = 0.5, const int &ITS = 50,
            solver_context *rest_ctx = nullptr) {

    vecX q;
    vecX s;
    set_mass(blocks);

    for (auto &block : blocks) {
      block->map_to_x(q);
    }

    if (rest_ctx) {
      for (int i = 0; i < static_cast<int>(_constraints.size()); i++) {
        if (auto bc =
                std::dynamic_pointer_cast<block_constraint>(_constraints[i])) {
          if (!bc->rest_initialized())
            bc->init_rest(q, *rest_ctx);
        }
      }
    }

    for (auto &block : blocks) {
      block->integrate_inertia(h);
      block->map_to_x(s);
    }
    // Soler et al. Alg.1 line 6: start local/global from the momentum prediction.
    // Without this, project() runs on q^t and never sees torque-tipped frames, so
    // bend/twist springs cannot propagate twist to neighbors.
    q = s;
    int Nm = q.size();

    std::vector<trip> triplets;
    index_t id0 = 0;

    for (int i = 0; i < _constraints.size(); i++) {
      auto &constraint = _constraints[i];
      constraint->fill_A(id0, triplets);
    }

    const detail::linear_system_cache_key key{
        Nm,
        id0,
        _constraints.size(),
        detail::hash_triplets(triplets),
        h,
        detail::mass_diagonal_sum(__M),
    };

    const bool diag = projection_diag_enabled();
    matS A;
    m_solver *S = nullptr;

    if (_cache_valid && key == _cache_key) {
      if (diag)
        projection_diag_log("linear_system cache hit");
      A = _cache_A;
      S = _cache_solver.get();
    } else {
      if (diag)
        projection_diag_log("linear_system cache miss");
      A = matS(id0, Nm);
      matS M = __M;
      M *= 1.0 / h / h;
      if (diag) {
        projection_diag_log("setFromTriplets...");
        std::cout << "Nm: " << Nm << std::endl;
      }
      A.setFromTriplets(triplets.begin(), triplets.end());
      if (diag)
        log_matrix_report("A", diagnose_sparse(A));
      matS AtA = A.transpose() * A;
      if (diag)
        log_matrix_report("AtA", diagnose_sparse(AtA));
      matS MAtA = M + AtA;
      if (diag)
        log_matrix_report("MAtA", diagnose_sparse(MAtA));

      _cache_solver = std::make_unique<m_solver>(MAtA);
      _cache_A = A;
      _cache_AtA = AtA;
      _cache_MAtA = MAtA;
      _cache_key = key;
      _cache_valid = true;
      S = _cache_solver.get();
    }

    if (!S->success()) {
      std::cout << "Solve failed, attempting normalization" << std::endl;
      matS MAtA = _cache_MAtA;
      Eigen::SparseMatrix<double> I(MAtA.rows(), MAtA.cols());
      MAtA += 1e-6 * I;
      _cache_solver = std::make_unique<m_solver>(MAtA);
      _cache_MAtA = MAtA;
      S = _cache_solver.get();
    }
    if (!S->success()) {
      std::cout << "Solve failed" << std::endl;
      return;
    }

    vecX p = vecX::Zero(id0);
    vecX q0 = q;
    for (int k = 0; k < ITS; k++) {
      p.setZero();
      for (int i = 0; i < _constraints.size(); i++) {
        auto &constraint = _constraints[i];
        constraint->project(q, p);
      }
      if (q.hasNaN()) {
        std::cout << "q has NaN" << std::endl;
      }

      if (p.hasNaN()) {
        std::cout << "p has NaN" << std::endl;
      }

      if (diag && (k == 0 || k % 10 == 0))
        std::cout << "k: " << k << " -pnorm: " << p.norm() << std::endl;

      matS M = __M;
      M *= 1.0 / h / h;
      vecX b = M * s + A.transpose() * p;

      q = S->solve(b);

      real dq = (q - q0).norm();
      if (dq < 1e-8)
        break;
      q0 = q;

      if (diag && k % 10 == 0)
        std::cout << "k: " << k << " -norms: "
                  << " q-q0: " << dq << std::endl;
    }

    for (auto &block : blocks) {
      block->map_from_x(q, h, damping);
    }
  }

  void
  set_constraints(const std::vector<projection_constraint::ptr> &constraints) {
    _cache_valid = false;
    _constraints = constraints;
  }

  matS __M;
  std::vector<projection_constraint::ptr> _constraints;

private:
  bool _cache_valid = false;
  detail::linear_system_cache_key _cache_key;
  matS _cache_A;
  matS _cache_AtA;
  matS _cache_MAtA;
  std::unique_ptr<m_solver> _cache_solver;
}; // class projection_solver
} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif
