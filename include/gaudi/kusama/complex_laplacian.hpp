#ifndef __GAUDI_KUSAMA_COMPLEX_LAPLACIAN_HPP__
#define __GAUDI_KUSAMA_COMPLEX_LAPLACIAN_HPP__

/// Implicit Crank--Nicolson for the linear (1+iα)Δ part of the complex Ginzburg--Landau
/// operator, in real 2×2 block form on (u, v). See block structure in implementation.

#include "gaudi/kusama/laplacian.hpp"
#include "gaudi/sparse_solver.h" // vecX, real
#include <cmath>
#include <vector>

#include <Eigen/SparseLU>

namespace gaudi {
namespace kusama {
namespace rx {
namespace cgle {

enum class linear_solve_mode { block_gauss_seidel, sparse_lu };

struct linear_config {
  linear_solve_mode mode = linear_solve_mode::block_gauss_seidel;
  int max_gs_iter = 80;
  real gs_tol = 1e-9;
  /// If 0, use the mean of `alpha_verts` for a scalar block `h α C` on off-diagonals.
  real alpha_scalar = 0.0;
};

namespace detail {

using sparmat = ::gaudi::kusama::laplacian::sparmat;
using triplet = Eigen::Triplet<real>;

inline void assemble_block_lu(int n, const sparmat &A, const sparmat &C, real h,
                              real alpha, sparmat *out) {
  std::vector<triplet> t;
  t.reserve(static_cast<size_t>(2) * (A.nonZeros() + 2 * C.nonZeros()));
  for (int k = 0; k < A.outerSize(); ++k) {
    for (Eigen::SparseMatrix<real>::InnerIterator it(A, k); it; ++it) {
      const int r = static_cast<int>(it.row());
      const int c = static_cast<int>(it.col());
      t.emplace_back(r, c, it.value());
      t.emplace_back(r + n, c + n, it.value());
    }
  }
  for (int k = 0; k < C.outerSize(); ++k) {
    for (Eigen::SparseMatrix<real>::InnerIterator it(C, k); it; ++it) {
      const int r = static_cast<int>(it.row());
      const int c = static_cast<int>(it.col());
      const real v = it.value();
      t.emplace_back(r, c + n, h * alpha * v);
      t.emplace_back(r + n, c, -h * alpha * v);
    }
  }
  out->resize(2 * n, 2 * n);
  out->setFromTriplets(t.begin(), t.end());
  out->makeCompressed();
}

} // namespace detail

/// One Crank--Nicolson step for the \((1+i\alpha)\Delta\) linear term (CGLE linear part).
inline void linear_crank_nicolson(::gaudi::kusama::laplacian &L,
                                  const std::vector<real> &alpha_verts, real h,
                                  std::vector<real> &u, std::vector<real> &v,
                                  const linear_config &cfg) {

  using sparmatL = ::gaudi::kusama::laplacian::sparmat;
  const sparmatL &M = *L.mass_ptr();
  const sparmatL &C = *L.stiffness_ptr();
  const sparmatL A = L.cn_diffusion_matrix(h);

  real alpha = cfg.alpha_scalar;
  if (alpha == 0.0) {
    if (alpha_verts.empty())
      alpha = 1.0;
    else {
      alpha = 0.0;
      for (real a : alpha_verts)
        alpha += a;
      alpha /= static_cast<real>(alpha_verts.size());
    }
  }

  vecX u0e = L.to_eigen(u);
  vecX v0e = L.to_eigen(v);
  const int n = static_cast<int>(u0e.size());

  vecX b_u = 2.0 * (M * u0e) + h * (C * u0e) - h * alpha * (C * v0e);
  vecX b_v = 2.0 * (M * v0e) + h * (C * v0e) + h * alpha * (C * u0e);

  if (cfg.mode == linear_solve_mode::sparse_lu) {
    sparmatL B(2 * n, 2 * n);
    detail::assemble_block_lu(n, A, C, h, alpha, &B);
    vecX bfull(2 * n);
    bfull.head(n) = b_u;
    bfull.tail(n) = b_v;
    Eigen::SparseLU<sparmatL, Eigen::COLAMDOrdering<int>> slu;
    slu.compute(B);
    if (slu.info() == Eigen::Success) {
      vecX sol = slu.solve(bfull);
      if (slu.info() == Eigen::Success) {
        u = L.from_eigen(sol.head(n));
        v = L.from_eigen(sol.tail(n));
        return;
      }
    }
  }

  // block Gauss–Seidel: A u = b_u - hαC v,  A v = b_v + hαC u
  auto solve_line = [&](const vecX &r) {
    return L.to_eigen(
        L.solve_system_copy(::gaudi::kusama::laplacian::sparmat(A), L.from_eigen(r)));
  };

  vecX u_k = u0e;
  vecX v_k = v0e;
  vecX u_last = u_k, v_last = v_k;
  for (int it = 0; it < cfg.max_gs_iter; ++it) {
    u_k = solve_line(b_u - h * alpha * (C * v_k));
    v_k = solve_line(b_v + h * alpha * (C * u_k));
    real ch = (u_k - u_last).norm() + (v_k - v_last).norm();
    u_last = u_k;
    v_last = v_k;
    if (it > 0 && ch < cfg.gs_tol * (1.0 + u_k.norm() + v_k.norm()))
      break;
  }
  u = L.from_eigen(u_k);
  v = L.from_eigen(v_k);
}

/// Ginzburg–Landau **linear** substep only (no pluggable real scalar diffuser).
struct linear_operator {
  static void apply(::gaudi::kusama::laplacian &L, const std::vector<real> &alpha, real h,
                    std::vector<real> &u, std::vector<real> &v,
                    const linear_config &cfg) {
    linear_crank_nicolson(L, alpha, h, u, v, cfg);
  }
};

} // namespace cgle
} // namespace rx
} // namespace kusama
} // namespace gaudi

#endif
