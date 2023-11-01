//
//  m2Includes.h
//  Manifold
//
//  Created by John Delaney on 5/22/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//
// includes files that are required for the data structure only, not derived
// files that USE the data structure
#ifndef __MLS_TORUS_ESTIMATION__
#define __MLS_TORUS_ESTIMATION__

#include "rod_integrators.hpp"
#include "shell_integrators.hpp"
#include <cmath>
#include <cstdlib>
#include <ostream>
#include <vector>

namespace gaudi {

namespace calder {

std::vector<vec3> calc_center(asawa::rod::rod &R, const std::vector<vec3> &Nr,
                              const std::vector<vec3> &p_pov, real l0,
                              real p = 3.0) {
  std::vector<real> weights = R.l0();
  std::vector<real> sums(p_pov.size(), 0.0);
  std::vector<vec3> us = integrate_over_rod<vec3>(
      R, p_pov,
      [&Nr, &weights](const std::vector<index_t> &edge_ids, Rod_Sum_Type &sum) {
        sum.bind(calder::vec3_datum::create(edge_ids, Nr));
        sum.bind(calder::scalar_datum::create(edge_ids, weights));
      },
      [&](const index_t i, const index_t j, //
          const vec3 &pi, const vec3 &pj,
          const std::vector<calder::datum::ptr> &data,
          Rod_Sum_Type::Node_Type node_type, //
          const Rod_Sum_Type::Node &node,    //
          const Rod_Sum_Type::Tree &tree) -> vec3 {
        const calder::vec3_datum::ptr F_datum =
            static_pointer_cast<calder::vec3_datum>(data[0]);
        vec3 Nr = get_data<vec3>(node_type, j, 0, data);
        real w = get_data<real>(node_type, j, 1, data);
        Nr.normalize();
        vec3 dp = pj - pi;
        real ndp = dp.squaredNorm();

        real nPdp = (Nr * Nr.transpose() * dp).norm();
        real r = 0.5 * ndp / (nPdp + 1.0 * l0);

        real dist = (pj - pi).norm();
        real kappa = computeK(dist, l0, p);
        sums[i] += w * kappa;
        return w * kappa * dp;
      });

  for (int i = 0; i < us.size(); i++) {
    if (sums[i] == 0.0)
      continue;
    us[i] /= sums[i];
  }

  return us;
}

} // namespace calder
} // namespace gaudi
#endif