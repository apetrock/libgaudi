/// Headless checks: iso fallback (σu=σv=1) matches isotropic cotan for random \f$g\f$;
/// sign flip invariance; conductance vs fem_d Frobenius distance at anisotropic weights.

#include "gaudi/asawa/datums.hpp"
#include "gaudi/asawa/faceloader.hpp"
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/kusama/laplacian.hpp"
#include "gaudi/kusama/laplacian_anisotropic.hpp"
#include "gaudi/kusama/vector_dirichlet.hpp"

#include <Eigen/Sparse>
#include <cassert>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

namespace {

gaudi::kusama::real sparse_frob_diff(const Eigen::SparseMatrix<gaudi::kusama::real> &A,
                             const Eigen::SparseMatrix<gaudi::kusama::real> &B) {
  assert(A.rows() == B.rows() && A.cols() == B.cols());
  const int n = static_cast<int>(A.rows());
  Eigen::MatrixXd Da = Eigen::MatrixXd::Zero(n, n);
  Eigen::MatrixXd Db = Eigen::MatrixXd::Zero(n, n);
  for (int k = 0; k < A.outerSize(); ++k) {
    for (Eigen::SparseMatrix<gaudi::kusama::real>::InnerIterator it(A, k); it;
         ++it)
      Da(it.row(), it.col()) = it.value();
  }
  for (int k = 0; k < B.outerSize(); ++k) {
    for (Eigen::SparseMatrix<gaudi::kusama::real>::InnerIterator it(B, k); it;
         ++it)
      Db(it.row(), it.col()) = it.value();
  }
  return (Da - Db).norm();
}

} // namespace

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

  const Eigen::SparseMatrix<gaudi::kusama::real> L_iso = kusama::build_lap<1>(
      *M, x,
      [](asawa::shell::shell &Ms, asawa::shell::CornerId c,
         const std::vector<vec3> &xs) -> gaudi::kusama::real {
        return asawa::shell::weak_cotan(Ms, Ms.prev(c), xs) +
               asawa::shell::weak_cotan(Ms, Ms.prev(Ms.other(c)), xs);
      });

  std::vector<int> slot_map;
  const int nE = kusama::build_compact_edge_dof_map(*M, slot_map);
  assert(nE > 0);

  std::mt19937 rng(42);
  std::normal_distribution<double> nd(0.0, 1.0);
  Eigen::VectorXd g = Eigen::VectorXd::Zero(2 * nE);
  for (int i = 0; i < 2 * nE; ++i)
    g[i] = static_cast<gaudi::kusama::real>(nd(rng));

  const Eigen::SparseMatrix<gaudi::kusama::real> L_c =
      kusama::build_anisotropic_cotan_laplacian(
          *M, x, g, slot_map, nE, gaudi::kusama::real(1), gaudi::kusama::real(1),
          kusama::anisotropic_laplacian_kind::conductance);
  const Eigen::SparseMatrix<gaudi::kusama::real> L_f =
      kusama::build_anisotropic_cotan_laplacian(
          *M, x, g, slot_map, nE, gaudi::kusama::real(1), gaudi::kusama::real(1),
          kusama::anisotropic_laplacian_kind::fem_d);

  const gaudi::kusama::real d_iso_c = sparse_frob_diff(L_iso, L_c);
  const gaudi::kusama::real d_iso_f = sparse_frob_diff(L_iso, L_f);
  assert(d_iso_c < static_cast<gaudi::kusama::real>(1e-10));
  assert(d_iso_f < static_cast<gaudi::kusama::real>(1e-10));
  std::cout << "[anisotropic_cotan_smoke] iso match conductance Frobenius diff " << d_iso_c
            << ", fem_d " << d_iso_f << "\n";

  Eigen::VectorXd g_neg = -g;
  const Eigen::SparseMatrix<gaudi::kusama::real> L_c_neg =
      kusama::build_anisotropic_cotan_laplacian(
          *M, x, g_neg, slot_map, nE, gaudi::kusama::real(1),
          gaudi::kusama::real(0.01),
          kusama::anisotropic_laplacian_kind::conductance);
  const Eigen::SparseMatrix<gaudi::kusama::real> L_c_pos =
      kusama::build_anisotropic_cotan_laplacian(
          *M, x, g, slot_map, nE, gaudi::kusama::real(1),
          gaudi::kusama::real(0.01),
          kusama::anisotropic_laplacian_kind::conductance);
  assert(sparse_frob_diff(L_c_neg, L_c_pos) <
         static_cast<gaudi::kusama::real>(1e-10));
  std::cout << "[anisotropic_cotan_smoke] sign flip invariance (conductance, aniso) ok\n";

  const Eigen::SparseMatrix<gaudi::kusama::real> L_f_neg =
      kusama::build_anisotropic_cotan_laplacian(
          *M, x, g_neg, slot_map, nE, gaudi::kusama::real(1),
          gaudi::kusama::real(0.01), kusama::anisotropic_laplacian_kind::fem_d);
  const Eigen::SparseMatrix<gaudi::kusama::real> L_f_pos =
      kusama::build_anisotropic_cotan_laplacian(
          *M, x, g, slot_map, nE, gaudi::kusama::real(1),
          gaudi::kusama::real(0.01), kusama::anisotropic_laplacian_kind::fem_d);
  assert(sparse_frob_diff(L_f_neg, L_f_pos) <
         static_cast<gaudi::kusama::real>(1e-10));
  std::cout << "[anisotropic_cotan_smoke] sign flip invariance (fem_d, aniso) ok\n";

  const Eigen::SparseMatrix<gaudi::kusama::real> L_f_an =
      kusama::build_anisotropic_cotan_laplacian(
          *M, x, g, slot_map, nE, gaudi::kusama::real(1),
          gaudi::kusama::real(0.01),
          kusama::anisotropic_laplacian_kind::fem_d);
  const gaudi::kusama::real d_kinds = sparse_frob_diff(L_c_pos, L_f_an);
  const gaudi::kusama::real n_iso = L_iso.norm();
  std::cout << "[anisotropic_cotan_smoke] ‖L_cond - L_fem‖_F / ‖L_iso‖_F = "
            << (n_iso > static_cast<gaudi::kusama::real>(1e-30) ? d_kinds / n_iso
                                                                  : d_kinds)
            << "\n";

  std::cout << "[anisotropic_cotan_smoke] all passed\n";
  return 0;
}
