#ifndef GAUDI_ENGMAN_STORE_METRIC_HPP
#define GAUDI_ENGMAN_STORE_METRIC_HPP

/// The store's step metric (`scripts/sympy/cr_poisson/MATH.md` section 17.6):
/// the split spin-2 operator `lambda_H L_tau + lambda_z L_z` added to the
/// store block of the Engman Hessian. Built from `edge_shape_fair` rows
/// `dirichlet_trace` and `dirichlet_traceless` toward the live store
/// (`goal_live`): their residual and gradient are zero at every iterate, so
/// they shape the steps (smooth mean curvature and principal directions) and
/// leave the solution where it is.

#include "gaudi/hepworth/block/edge_shape_fair.hpp"
#include "gaudi/hepworth/block/edge_shape_goals.hpp"

#include <cmath>
#include <vector>

namespace gaudi {
namespace engman {

/// Energies for `solver::set_energies` (append to the others). Row weights
/// are `sqrt(lambda)`, so `lambda_H`, `lambda_z` weigh the energies.
inline std::vector<hepworth::projection_constraint::ptr> store_metric(const kusama::stencil_table *table,
                                                                      const std::vector<vec3> &x,
                                                                      hepworth::sim_block::ptr store,
                                                                      real lambda_H, real lambda_z) {
  namespace hb = hepworth::block;
  std::vector<hepworth::projection_constraint::ptr> out;
  if (lambda_H > real(0))
    out.push_back(hb::edge_shape_fair::create(table, x, hb::edge_fair_rows::dirichlet_trace, std::sqrt(lambda_H),
                                              std::vector<hepworth::sim_block::ptr>{store}, hb::goal_live(store)));
  if (lambda_z > real(0))
    out.push_back(hb::edge_shape_fair::create(table, x, hb::edge_fair_rows::dirichlet_traceless,
                                              std::sqrt(lambda_z), std::vector<hepworth::sim_block::ptr>{store},
                                              hb::goal_live(store)));
  return out;
}

} // namespace engman
} // namespace gaudi

#endif
