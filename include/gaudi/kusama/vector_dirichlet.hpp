#ifndef GAUDI_KUSAMA_VECTOR_DIRICHLET_HPP
#define GAUDI_KUSAMA_VECTOR_DIRICHLET_HPP

/// Stein et al., "A Simple Discretization of the Vector Dirichlet Energy" (CGF 2020).
/// Assembly matches libigl \c cr_vector_laplacian_intrinsic from the authors' reference
/// repository (Crouzeix–Raviart / edge-based tangential DOFs, \f$2|E|\times 2|E|\f$).

#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/common.h"

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace gaudi {
namespace kusama {

/// Contiguous dof index for each physical edge slot \c c/2 (\c c from \c get_edge_range).
inline int build_compact_edge_dof_map(asawa::shell::shell &M,
                                      std::vector<int> &edge_slot_to_dof) {
  const int ns = static_cast<int>(M.corner_count() / 2);
  edge_slot_to_dof.assign(ns, -1);
  int nE = 0;
  for (asawa::shell::CornerId c : M.get_edge_range()) {
    const int slot = static_cast<int>(c) / 2;
    if (edge_slot_to_dof[static_cast<size_t>(slot)] < 0)
      edge_slot_to_dof[static_cast<size_t>(slot)] = nE++;
  }
  return nE;
}

inline int edge_orientation_sign(const asawa::shell::shell &M,
                                 asawa::shell::CornerId c) {
  const int ia = static_cast<int>(M.vert(c));
  const int ib = static_cast<int>(M.vert(M.other(c)));
  return ia <= ib ? 1 : -1;
}

namespace laplace {
namespace vector_dirichlet {

/// One triangle's contribution (6 x 6), local index `comp * 3 + edge`
/// (`comp` 0 parallel, 1 perpendicular; local edge `e = (X_e, X_{e+1})`).
/// `sign[e]` is +1 when local edge `e` runs from the lower vertex id to the
/// higher (the edge dof's canonical orientation), -1 otherwise. A degenerate
/// triangle returns zero.
inline Eigen::Matrix<real, 6, 6> face_element(const vec3 &X0, const vec3 &X1,
                                              const vec3 &X2,
                                              const std::array<int, 3> &sign) {
  Eigen::Matrix<real, 6, 6> K = Eigen::Matrix<real, 6, 6>::Zero();
  const vec3 X[3] = {X0, X1, X2};
  double lsq[3];
  for (int e = 0; e < 3; ++e)
    lsq[e] = static_cast<double>((X[(e + 1) % 3] - X[e]).squaredNorm());
  const double dA = (X1 - X0).cross(X2 - X0).norm();
  if (dA < 1e-30)
    return K;
  for (int e = 0; e < 3; ++e) {
    const double eij = lsq[e];
    const double ejk = lsq[(e + 1) % 3];
    const double eki = lsq[(e + 2) % 3];
    if (!(eij > 1e-24) || !(eki > 1e-24))
      continue;
    const double lens = std::sqrt(eij * eki);
    const int j = (e + 2) % 3;
    const int o = sign[static_cast<size_t>(e)] * sign[static_cast<size_t>(j)];
    const double diag = 2.0 / dA * eij;
    const double num = eij - ejk + eki;
    const double Dijki = static_cast<double>(o) * num * num / (2.0 * eij * eki * dA) * lens;
    const double Dperp = -static_cast<double>(o) * num / (eij * eki) * lens;
    // A sliver divides by a vanishing area and writes Inf. The normal
    // equations then turn that into NaN.
    if (!std::isfinite(diag) || !std::isfinite(Dijki) || !std::isfinite(Dperp))
      continue;
    K(e, e) += diag;
    K(3 + e, 3 + e) += diag;
    K(e, j) += Dijki;
    K(j, e) += Dijki;
    K(3 + e, 3 + j) += Dijki;
    K(3 + j, 3 + e) += Dijki;
    K(e, 3 + j) += Dperp;
    K(3 + j, e) += Dperp;
    K(3 + e, j) -= Dperp;
    K(j, 3 + e) -= Dperp;
  }
  return K;
}

/// Sparse symmetric stiffness for \f$\tfrac12 u^\top L u\f$ (two scalar DOFs per edge).
inline Eigen::SparseMatrix<real> energy(asawa::shell::shell &M,
                                        const std::vector<vec3> &x,
                                        int *out_num_edge_dofs = nullptr) {
  if (!M.verts_are_dense_packed()) {
    throw std::runtime_error(
        "kusama::laplace::vector_dirichlet::energy: shell vertices must be "
        "dense-packed (see kusama::build_lap)");
  }
  if (x.size() != static_cast<size_t>(M.vert_count())) {
    throw std::runtime_error(
        "kusama::laplace::vector_dirichlet::energy: position count mismatch");
  }
  std::vector<int> slot_map;
  const int nE = build_compact_edge_dof_map(M, slot_map);
  if (nE <= 0) {
    if (out_num_edge_dofs)
      *out_num_edge_dofs = 0;
    return Eigen::SparseMatrix<real>(0, 0);
  }

  std::vector<Eigen::Triplet<real>> triplets;
  auto frange = M.get_face_range(true);
  triplets.reserve(frange.size() * 48);

  for (asawa::shell::FaceId fi : frange) {
    if (M.fsize(fi) != 3)
      continue;
    bool first = true; // a triangle's fan repeats once, degenerate
    M.const_for_each_face_tri(
        fi, [&](asawa::shell::CornerId c0, asawa::shell::CornerId c1,
                asawa::shell::CornerId c2, const asawa::shell::shell &Ms) {
          if (!first)
            return;
          first = false;
          const asawa::shell::CornerId ce[3] = {c0, c1, c2};
          std::array<int, 3> sign{};
          int dof[3];
          for (int e = 0; e < 3; ++e) {
            sign[static_cast<size_t>(e)] = edge_orientation_sign(Ms, ce[e]);
            dof[e] = slot_map[static_cast<size_t>(static_cast<int>(ce[e]) / 2)];
            if (dof[e] < 0)
              return;
          }
          const Eigen::Matrix<real, 6, 6> K = face_element(
              x[Ms.vert(c0)], x[Ms.vert(c1)], x[Ms.vert(c2)], sign);
          for (int a = 0; a < 6; ++a)
            for (int b = 0; b < 6; ++b)
              if (K(a, b) != real(0))
                triplets.emplace_back(dof[a % 3] + (a / 3) * nE,
                                      dof[b % 3] + (b / 3) * nE, K(a, b));
        });
  }

  Eigen::SparseMatrix<real> L(2 * nE, 2 * nE);
  L.setFromTriplets(triplets.begin(), triplets.end());
  if (out_num_edge_dofs)
    *out_num_edge_dofs = nE;
  return L;
}

/// Pack a length \f$2 n_E\f$ DOF vector into \f$n_E\f$ per-edge `vec2` values.
inline std::vector<vec2> pack_dof_to_vec2(const Eigen::VectorXd &u,
                                          int num_edge_dofs) {
  if (u.size() != 2 * num_edge_dofs) {
    throw std::runtime_error(
        "kusama::laplace::vector_dirichlet::pack_dof_to_vec2: expected u.size() == 2 * "
        "num_edge_dofs (got " +
        std::to_string(static_cast<int>(u.size())) + " vs 2*" +
        std::to_string(num_edge_dofs) + ")");
  }
  std::vector<vec2> out(static_cast<size_t>(num_edge_dofs));
  for (int i = 0; i < num_edge_dofs; ++i)
    out[static_cast<size_t>(i)] = vec2(u(i), u(i + num_edge_dofs));
  return out;
}

} // namespace vector_dirichlet
} // namespace laplace

/// @deprecated Prefer `kusama::laplace::vector_dirichlet::energy`.
inline Eigen::SparseMatrix<real>
build_vector_dirichlet_energy(asawa::shell::shell &M,
                              const std::vector<vec3> &x,
                              int *out_num_edge_dofs = nullptr) {
  return laplace::vector_dirichlet::energy(M, x, out_num_edge_dofs);
}

/// @deprecated Prefer `kusama::laplace::vector_dirichlet::pack_dof_to_vec2`.
inline std::vector<vec2> pack_vector_dirichlet_dof_to_vec2(const Eigen::VectorXd &u,
                                                           int num_edge_dofs) {
  return laplace::vector_dirichlet::pack_dof_to_vec2(u, num_edge_dofs);
}

} // namespace kusama
} // namespace gaudi

#endif
