#include "gaudi/asawa/datums.hpp"
#include "gaudi/asawa/faceloader.hpp"
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/kusama/vector_dirichlet.hpp"
#include "gaudi/kusama/vector_dirichlet_guided.hpp"
#include "gaudi/common.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

int main() {
  using namespace gaudi;
  std::vector<vec3> V = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}};
  std::vector<std::vector<index_t>> F = {{0, 1, 2}, {1, 3, 2}};
  std::vector<index_t> cn, cv, cf;
  asawa::assemble_table(V, F, cn, cv, cf);
  asawa::shell::shell::ptr M = asawa::shell::shell::create(cn, cv, cf);
  asawa::init_vert_datum(*M, vec3(0, 0, 0));
  std::vector<vec3> &x = asawa::get_vec_data(*M, 0);
  x = V;

  int nE = 0;
  Eigen::SparseMatrix<double> L =
      kusama::build_vector_dirichlet_energy(*M, x, &nE);
  assert(nE > 0);
  assert(L.rows() == 2 * nE);
  assert(L.cols() == 2 * nE);
  assert(L.nonZeros() > 0);
  std::cout << "[vector_dirichlet_smoke] nE=" << nE << " nnz=" << L.nonZeros()
            << "\n";

  // edge_barycentric_dual_mass = 0.5 * (vert_area(v0) + vert_area(v1))
  {
    for (asawa::shell::CornerId c : M->get_edge_range()) {
      double m_e = asawa::shell::edge_barycentric_dual_mass(*M, c, x);
      double a0 = asawa::shell::vert_area(*M, M->vert(c), x);
      double a1 = asawa::shell::vert_area(*M, M->vert(M->next(c)), x);
      assert(std::abs(m_e - 0.5 * (a0 + a1)) < 1e-10);
      break;
    }
  }

  // Guided solve residual ||(L+λM)u - λMg||
  {
    std::vector<int> slot_map;
    int nE2 = kusama::build_compact_edge_dof_map(*M, slot_map);
    assert(nE2 == nE);
    std::vector<asawa::shell::CornerId> dof_to_corner;
    kusama::build_dof_to_corner(*M, slot_map, nE, dof_to_corner);
    std::vector<std::vector<asawa::shell::FaceId>> edge_inc;
    kusama::build_edge_incident_faces(*M, slot_map, nE, edge_inc);
    Eigen::VectorXd g = kusama::build_edge_curvature_guidance(
        *M, x, slot_map, nE, dof_to_corner, edge_inc,
        asawa::shell::face_curvature_stencil::one_ring);
    Eigen::VectorXd mass =
        kusama::build_edge_barycentric_mass_diagonal(*M, x, nE, dof_to_corner);
    const double lambda = 1e6;
    Eigen::VectorXd u =
        kusama::solve_guided_vector_dirichlet(L, mass, g, lambda);
    Eigen::SparseMatrix<double> A = L;
    for (int i = 0; i < 2 * nE; ++i)
      A.coeffRef(i, i) += lambda * mass(i);
    Eigen::VectorXd b = lambda * mass.cwiseProduct(g);
    Eigen::VectorXd r = A * u - b;
    const double bn = std::max(1.0, b.norm());
    assert(r.norm() < 1e-5 * bn);
    std::cout << "[vector_dirichlet_smoke] guided residual ratio "
              << (r.norm() / bn) << "\n";
  }

  // Sign coherence idempotent on oriented guidance
  {
    std::vector<int> slot_map;
    int nE2 = kusama::build_compact_edge_dof_map(*M, slot_map);
    (void)nE2;
    std::vector<asawa::shell::CornerId> dof_to_corner;
    kusama::build_dof_to_corner(*M, slot_map, nE, dof_to_corner);
    std::vector<std::vector<asawa::shell::FaceId>> edge_inc;
    kusama::build_edge_incident_faces(*M, slot_map, nE, edge_inc);
    Eigen::VectorXd g = kusama::build_edge_curvature_guidance(
        *M, x, slot_map, nE, dof_to_corner, edge_inc,
        asawa::shell::face_curvature_stencil::one_ring);
    std::vector<vec3> edge_n(static_cast<size_t>(nE), vec3::UnitZ());
    for (int e = 0; e < nE; ++e)
      edge_n[static_cast<size_t>(e)] =
          kusama::edge_average_normal(*M, x, edge_inc[static_cast<size_t>(e)]);
    std::vector<std::vector<kusama::edge_triangle_link>> adj;
    kusama::build_edge_triangle_adjacency(*M, slot_map, nE, adj);
    kusama::orient_edge_guidance_sign_coherence(*M, x, nE, dof_to_corner,
                                                  edge_n, adj, g);
    Eigen::VectorXd once = g;
    kusama::orient_edge_guidance_sign_coherence(*M, x, nE, dof_to_corner,
                                                  edge_n, adj, g);
    assert((g - once).norm() < 1e-12);
    std::cout << "[vector_dirichlet_smoke] sign coherence idempotent ok\n";
  }

  return 0;
}
