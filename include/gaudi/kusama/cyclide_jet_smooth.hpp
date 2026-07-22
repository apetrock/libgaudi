#ifndef GAUDI_KUSAMA_CYCLIDE_JET_SMOOTH_HPP
#define GAUDI_KUSAMA_CYCLIDE_JET_SMOOTH_HPP

#include <algorithm>
#include <stdexcept>
#include <vector>

#include "cyclide_smooth_generated.hpp"
#include "gaudi/albers/darboux_cyclide.hpp"
#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/asawa/shell/shell.hpp"

namespace gaudi {
namespace kusama {

struct cyclide_jet_smooth_params {
  real wi = 0.9;
  real alpha_G = 0.5;
  real alpha_H = 0.25;
  bool use_cotan_weights = true;
  // Jacobi sweeps: each vert solves locally vs neighbors from the previous
  // iterate (sweep 0 uses the raw fit). Anchor always stays on the original Q0.
  int sweeps = 1;
};

inline std::vector<albers::vec14> cyclide_jet_smooth(
    asawa::shell::shell &M, const std::vector<vec3> &pov,
    const std::vector<albers::vec14> &Q_global,
    const cyclide_jet_smooth_params &params) {

  const size_t n = pov.size();
  if (n != Q_global.size()) {
    throw std::runtime_error("cyclide_jet_smooth: pov/Q_global size mismatch");
  }
  if (static_cast<size_t>(M.vert_count()) != n) {
    throw std::runtime_error("cyclide_jet_smooth: vert_count mismatch");
  }

  std::vector<real> edge_cotan;
  if (params.use_cotan_weights) {
    edge_cotan = asawa::shell::edge_cotan_weights(M, pov);
  }

  const int n_sweeps = std::max(1, params.sweeps);
  std::vector<albers::vec14> Q_prev = Q_global;
  std::vector<albers::vec14> Q_smooth = Q_global;
  const albers::vec14 Qi_zero = albers::vec14::Zero();
  albers::mat14 H_acc;
  albers::mat14 H_part;
  albers::vec14 b_acc;
  albers::vec14 g_part;

  for (int sweep = 0; sweep < n_sweeps; ++sweep) {
    for (size_t i = 0; i < n; ++i) {
      const asawa::shell::VertId vi =
          asawa::shell::vert_id(static_cast<asawa::shell::index_t>(i));
      const albers::vec14 &Q0_i = Q_global[i];

      H_acc.setZero();
      b_acc.setZero();

      albers::medial_generated::cyclide_smooth_anchor_hess(params.wi, H_acc);
      albers::medial_generated::cyclide_smooth_anchor_grad(Qi_zero, Q0_i,
                                                           params.wi, g_part);
      b_acc += g_part;

      M.for_each_vertex(vi, [&](asawa::shell::CornerId c,
                                asawa::shell::shell &Ms) {
        const int j = Ms.vert(Ms.next(c));
        if (j < 0 || static_cast<size_t>(j) >= n) {
          return;
        }
        const vec3 x_j_in_i = pov[static_cast<size_t>(j)] - pov[i];
        real w_ij = 1.0;
        if (params.use_cotan_weights) {
          const size_t eid = static_cast<size_t>(c) / 2;
          if (eid < edge_cotan.size()) {
            w_ij = edge_cotan[eid];
          }
        }
        if (w_ij <= real(0)) {
          w_ij = real(1e-12);
        }

        // Neighbors from previous Jacobi iterate (raw fit on sweep 0).
        const albers::vec14 &Q_j = Q_prev[static_cast<size_t>(j)];

        albers::medial_generated::cyclide_smooth_neighbor_hess(
            Q_j, x_j_in_i, params.wi, w_ij, params.alpha_G, params.alpha_H,
            H_part);
        H_acc += H_part;

        albers::medial_generated::cyclide_smooth_neighbor_grad(
            Qi_zero, Q_j, x_j_in_i, params.wi, w_ij, params.alpha_G,
            params.alpha_H, g_part);
        b_acc += g_part;
      });

      const Eigen::LDLT<albers::mat14> solver(H_acc);
      if (solver.info() != Eigen::Success) {
        continue;
      }
      const albers::vec14 Qi_sol = solver.solve(-b_acc);
      if (Qi_sol.allFinite()) {
        Q_smooth[i] = Qi_sol;
      }
    }
    Q_prev = Q_smooth;
  }

  return Q_smooth;
}

} // namespace kusama
} // namespace gaudi

#endif // GAUDI_KUSAMA_CYCLIDE_JET_SMOOTH_HPP
