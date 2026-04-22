#include "gaudi/asawa/datums.hpp"
#include "gaudi/asawa/faceloader.hpp"
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/bontecou/laplace_spectrum.hpp"
#include "gaudi/bontecou/spectrum.hpp"
#include "gaudi/bontecou/vector_dirichlet.hpp"
#include "gaudi/bontecou/vector_dirichlet_spectrum.hpp"
#include "gaudi/common.h"
#include "gaudi/duchamp/spectral_projection_integrator.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

int main() {
  using namespace gaudi;

  // --- spectrum.hpp on tiny SPD diagonal ---
  {
    const int n = 8;
    std::vector<Eigen::Triplet<double>> tr;
    for (int i = 0; i < n; ++i)
      tr.emplace_back(i, i, static_cast<double>(i + 1));
    Eigen::SparseMatrix<double> A(n, n);
    A.setFromTriplets(tr.begin(), tr.end());
    A.makeCompressed();
    Eigen::VectorXd evals;
    Eigen::MatrixXd evecs;
    const bool ok =
        bontecou::sparse_sym_eigs_smallest_algebraic(A, 3, 0, evals, evecs, nullptr);
    assert(ok);
    assert(evals.size() == 3);
    std::vector<double> got(3);
    for (int i = 0; i < 3; ++i)
      got[static_cast<size_t>(i)] = evals[i];
    std::sort(got.begin(), got.end());
    assert(std::abs(got[0] - 1.0) < 1e-4);
    assert(std::abs(got[1] - 2.0) < 1e-4);
    assert(std::abs(got[2] - 3.0) < 1e-4);
    std::cout << "[spectral_smoke] sparse_sym_eigs_smallest_algebraic ok\n";
  }

  // --- vector Dirichlet pack ---
  {
    int nE = 0;
    std::vector<vec3> V = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}};
    std::vector<std::vector<index_t>> F = {{0, 1, 2}, {1, 3, 2}};
    std::vector<index_t> cn, cv, cf;
    asawa::assemble_table(V, F, cn, cv, cf);
    asawa::shell::shell::ptr M = asawa::shell::shell::create(cn, cv, cf);
    asawa::init_vert_datum(*M, vec3(0, 0, 0));
    std::vector<vec3> &x = asawa::get_vec_data(*M, 0);
    x = V;
    Eigen::SparseMatrix<double> L = bontecou::build_vector_dirichlet_energy(*M, x, &nE);
    (void)L;
    Eigen::VectorXd u = Eigen::VectorXd::Ones(2 * nE);
    std::vector<vec2> pv = bontecou::pack_vector_dirichlet_dof_to_vec2(u, nE);
    assert(static_cast<int>(pv.size()) == nE);
    for (int i = 0; i < nE; ++i) {
      assert(std::abs(pv[static_cast<size_t>(i)][0] - 1.0) < 1e-12);
      assert(std::abs(pv[static_cast<size_t>(i)][1] - 1.0) < 1e-12);
    }
    std::cout << "[spectral_smoke] pack_vector_dirichlet_dof_to_vec2 ok\n";
  }

  // --- spectral_projection_integrator on 2-tri square ---
  {
    std::vector<vec3> V = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}};
    std::vector<std::vector<index_t>> F = {{0, 1, 2}, {1, 3, 2}};
    std::vector<index_t> cn, cv, cf;
    asawa::assemble_table(V, F, cn, cv, cf);
    asawa::shell::shell::ptr M = asawa::shell::shell::create(cn, cv, cf);
    asawa::init_vert_datum(*M, vec3(0, 0, 0));
    std::vector<vec3> &x = asawa::get_vec_data(*M, 0);
    x = V;

    duchamp::spectral_projection_config cfg;
    cfg.requested_modes = 4;
    cfg.epsilon_shift = 1e-6;
    cfg.band = duchamp::laplace_modes_band::low_frequency;
    duchamp::spectral_projection_integrator I(cfg);
    I.set_mesh(M);
    assert(I.rebuild_operator());
    assert(I.compute_modes(nullptr));
    assert(I.mode_count() >= 1);
    const Eigen::VectorXd f = I.vertex_field();
    I.project_field_onto_current_basis(f);
    assert((I.coeffs() - I.evecs().transpose() * f).norm() < 1e-8);

    I.remember_embedding_before_deformation();
    x[0][0] += 1e-4;
    assert(I.rebuild_operator());
    assert(I.compute_modes(nullptr, false));
    I.project_embedding_after_new_spectrum();
    const Eigen::VectorXd f2 = I.vertex_field();
    assert(f2.size() == f.size());
    assert((f2 - f).norm() < 5e-3);
    std::cout << "[spectral_smoke] spectral_projection_integrator projection ok\n";
  }

  // --- injected vertex operator builder (diagonal perturbation of default) ---
  {
    std::vector<vec3> V = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}};
    std::vector<std::vector<index_t>> F = {{0, 1, 2}, {1, 3, 2}};
    std::vector<index_t> cn, cv, cf;
    asawa::assemble_table(V, F, cn, cv, cf);
    asawa::shell::shell::ptr M = asawa::shell::shell::create(cn, cv, cf);
    asawa::init_vert_datum(*M, vec3(0, 0, 0));
    std::vector<vec3> &x = asawa::get_vec_data(*M, 0);
    x = V;

    duchamp::spectral_projection_config cfg;
    cfg.requested_modes = 4;
    cfg.epsilon_shift = 1e-6;
    cfg.band = duchamp::laplace_modes_band::low_frequency;

    auto builder = [](asawa::shell::shell::ptr Mp, std::vector<vec3> &xp,
                      const duchamp::spectral_projection_config &c) -> Eigen::SparseMatrix<double> {
      Eigen::SparseMatrix<double> A = duchamp::make_default_cotan_vertex_operator(Mp, xp, c);
      const int n = static_cast<int>(A.rows());
      for (int i = 0; i < n; ++i)
        A.coeffRef(i, i) += 1e-4;
      A.makeCompressed();
      return A;
    };

    duchamp::spectral_projection_integrator I(cfg, builder, {});
    I.set_mesh(M);
    assert(I.rebuild_operator());
    assert(I.compute_modes(nullptr));
    assert(I.mode_count() >= 1);
    std::cout << "[spectral_smoke] injected operator builder ok\n";
  }

  // --- vector Dirichlet (Stein) partial spectrum on 2-tri square ---
  {
    std::vector<vec3> V = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}};
    std::vector<std::vector<index_t>> F = {{0, 1, 2}, {1, 3, 2}};
    std::vector<index_t> cn, cv, cf;
    asawa::assemble_table(V, F, cn, cv, cf);
    asawa::shell::shell::ptr M = asawa::shell::shell::create(cn, cv, cf);
    asawa::init_vert_datum(*M, vec3(0, 0, 0));
    std::vector<vec3> &x = asawa::get_vec_data(*M, 0);
    x = V;
    int nE = 0;
    Eigen::VectorXd evals;
    Eigen::MatrixXd evecs;
    const bool ok = bontecou::vector_dirichlet_partial_spectrum(
        *M, x, 1e-8, bontecou::vector_dirichlet_spectrum_band::smallest_algebraic, 0.0, 4, 0,
        evals, evecs, &nE, nullptr);
    assert(ok);
    assert(nE > 0);
    assert(static_cast<int>(evals.size()) >= 3);
    assert(evecs.rows() == 2 * nE);
    std::vector<vec2> line = bontecou::pack_vector_dirichlet_dof_to_vec2(evecs.col(0), nE);
    assert(static_cast<int>(line.size()) == nE);
    std::cout << "[spectral_smoke] vector_dirichlet_partial_spectrum ok (nE=" << nE << ")\n";
  }

  // --- advance_projection_step ---
  {
    std::vector<vec3> V = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}};
    std::vector<std::vector<index_t>> F = {{0, 1, 2}, {1, 3, 2}};
    std::vector<index_t> cn, cv, cf;
    asawa::assemble_table(V, F, cn, cv, cf);
    asawa::shell::shell::ptr M = asawa::shell::shell::create(cn, cv, cf);
    asawa::init_vert_datum(*M, vec3(0, 0, 0));
    std::vector<vec3> &x = asawa::get_vec_data(*M, 0);
    x = V;

    duchamp::spectral_projection_config cfg;
    cfg.requested_modes = 4;
    cfg.epsilon_shift = 1e-6;
    cfg.band = duchamp::laplace_modes_band::low_frequency;
    duchamp::spectral_projection_integrator I(cfg);
    I.set_mesh(M);
    assert(I.rebuild_operator());
    assert(I.compute_modes(nullptr));
    const Eigen::VectorXd f0 = I.vertex_field();
    assert(I.advance_projection_step(1e-5, 1.0, true));
    const Eigen::VectorXd f1 = I.vertex_field();
    assert(f0.size() == f1.size());
    assert(f1.norm() < 1e3);
    std::cout << "[spectral_smoke] advance_projection_step ok\n";
  }

  // --- Hungarian matches brute-force optimum on small overlap costs ---
  {
    const int k = 4;
    Eigen::MatrixXd prev = Eigen::MatrixXd::Random(10, k);
    Eigen::MatrixXd raw = Eigen::MatrixXd::Random(10, k);
    const duchamp::mode_alignment h = duchamp::compute_mode_alignment_hungarian(prev, raw);
    double hung = 0;
    for (int j = 0; j < k; ++j) {
      const int i = static_cast<int>(h.raw_col_for_tracked[static_cast<size_t>(j)]);
      hung += std::abs(prev.col(j).dot(raw.col(i)));
    }
    std::vector<int> p(static_cast<size_t>(k));
    for (int i = 0; i < k; ++i)
      p[static_cast<size_t>(i)] = i;
    double best = -1e300;
    do {
      double s = 0;
      for (int j = 0; j < k; ++j)
        s += std::abs(prev.col(j).dot(raw.col(p[static_cast<size_t>(j)])));
      best = std::max(best, s);
    } while (std::next_permutation(p.begin(), p.end()));
    assert(std::abs(hung - best) < 1e-8);
    std::cout << "[spectral_smoke] hungarian alignment (vs brute perm) ok\n";
  }

  // --- truncate_coeffs_by_magnitude ---
  {
    Eigen::VectorXd c(5);
    c << 1.0, 5.0, 2.0, 4.0, 3.0;
    duchamp::truncate_coeffs_by_magnitude(c, 2);
    int nz = 0;
    for (int i = 0; i < 5; ++i)
      if (std::abs(c[i]) > 1e-12)
        ++nz;
    assert(nz == 2);
    assert(std::abs(c[1] - 5.0) < 1e-12 && std::abs(c[3] - 4.0) < 1e-12);
    std::cout << "[spectral_smoke] truncate_coeffs_by_magnitude ok\n";
  }

  // --- greedy alignment (permutation recovery) ---
  {
    Eigen::MatrixXd prev(4, 2);
    prev << 1, 0, 0, 1, 0, 0, 0, 0;
    Eigen::MatrixXd raw(4, 2);
    raw << 0, 1, 1, 0, 0, 0, 0, 0;
    const duchamp::mode_alignment al =
        duchamp::compute_mode_alignment_greedy(prev, raw);
    Eigen::MatrixXd E = raw;
    Eigen::VectorXd evals(2);
    evals << 0.1, 0.2;
    Eigen::VectorXd c(2);
    c = raw.transpose() * Eigen::VectorXd::Ones(4);
    duchamp::apply_mode_alignment(al, E, evals, c);
    assert((E.col(0) - prev.col(0)).norm() < 1e-10 ||
           (E.col(0) + prev.col(0)).norm() < 1e-10);
    std::cout << "[spectral_smoke] mode_alignment ok\n";
  }

  std::cout << "[spectral_smoke] all passed\n";
  return 0;
}
