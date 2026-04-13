#ifndef GAUDI_BONTECOU_VECTOR_DIRICHLET_HPP
#define GAUDI_BONTECOU_VECTOR_DIRICHLET_HPP

/// Stein et al., "A Simple Discretization of the Vector Dirichlet Energy" (CGF 2020).
/// Assembly matches libigl \c cr_vector_laplacian_intrinsic from the authors' reference
/// repository (Crouzeix–Raviart / edge-based tangential DOFs, \f$2|E|\times 2|E|\f$).

#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/common.h"

#include <Eigen/Sparse>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace gaudi {
namespace bontecou {

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
  const int ib = static_cast<int>(M.vert(M.next(c)));
  return ia < ib ? 1 : -1;
}

/// Sparse symmetric stiffness for \f$\tfrac12 u^\top L u\f$ (two scalar DOFs per edge).
inline Eigen::SparseMatrix<real>
build_vector_dirichlet_energy(asawa::shell::shell &M,
                              const std::vector<vec3> &x,
                              int *out_num_edge_dofs = nullptr) {
  if (!M.verts_are_dense_packed()) {
    throw std::runtime_error(
        "bontecou::build_vector_dirichlet_energy: shell vertices must be "
        "dense-packed (see bontecou::build_lap)");
  }
  if (x.size() != static_cast<size_t>(M.vert_count())) {
    throw std::runtime_error(
        "bontecou::build_vector_dirichlet_energy: position count mismatch");
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
    M.const_for_each_face_tri(
        fi, [&](asawa::shell::CornerId c0, asawa::shell::CornerId c1,
                asawa::shell::CornerId c2, const asawa::shell::shell &Ms) {
          const asawa::shell::CornerId ce[3] = {c0, c1, c2};
          double lsq[3];
          for (int e = 0; e < 3; ++e) {
            const real len =
                asawa::shell::edge_length(Ms, ce[static_cast<size_t>(e)], x);
            lsq[static_cast<size_t>(e)] = static_cast<double>(len * len);
          }
          const vec3 X0 = x[Ms.vert(c0)];
          const vec3 X1 = x[Ms.vert(c1)];
          const vec3 X2 = x[Ms.vert(c2)];
          const double dA = (X1 - X0).cross(X2 - X0).norm();
          if (dA < 1e-30)
            return;

          for (int e = 0; e < 3; ++e) {
            const double eij = lsq[static_cast<size_t>(e)];
            const double ejk = lsq[static_cast<size_t>((e + 1) % 3)];
            const double eki = lsq[static_cast<size_t>((e + 2) % 3)];
            const double lens = std::sqrt(eij * eki);
            const int o = edge_orientation_sign(Ms, ce[static_cast<size_t>(e)]) *
                          edge_orientation_sign(
                              Ms, ce[static_cast<size_t>((e + 2) % 3)]);
            const int Ei =
                slot_map[static_cast<size_t>(static_cast<int>(ce[static_cast<size_t>(e)]) /
                                             2)];
            const int Ej = slot_map[static_cast<size_t>(
                static_cast<int>(ce[static_cast<size_t>((e + 2) % 3)]) / 2)];
            if (Ei < 0 || Ej < 0)
              return;

            const double diag = 2.0 / dA * eij;
            triplets.emplace_back(Ei, Ei, diag);
            triplets.emplace_back(Ei + nE, Ei + nE, diag);

            const double num = eij - ejk + eki;
            const double Dijki =
                static_cast<double>(o) * num * num / (2.0 * eij * eki * dA) * lens;
            triplets.emplace_back(Ei, Ej, Dijki);
            triplets.emplace_back(Ej, Ei, Dijki);
            triplets.emplace_back(Ei + nE, Ej + nE, Dijki);
            triplets.emplace_back(Ej + nE, Ei + nE, Dijki);

            const double Dperp =
                -static_cast<double>(o) * num / (eij * eki) * lens;
            triplets.emplace_back(Ei, Ej + nE, Dperp);
            triplets.emplace_back(Ej + nE, Ei, Dperp);
            triplets.emplace_back(Ei + nE, Ej, -Dperp);
            triplets.emplace_back(Ej, Ei + nE, -Dperp);
          }
        });
  }

  Eigen::SparseMatrix<real> L(2 * nE, 2 * nE);
  L.setFromTriplets(triplets.begin(), triplets.end());
  if (out_num_edge_dofs)
    *out_num_edge_dofs = nE;
  return L;
}

} // namespace bontecou
} // namespace gaudi

#endif
