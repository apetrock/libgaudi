//
//  m2Includes.h
//  Manifold
//
//  Created by John Delaney on 5/22/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//
// includes files that are required for the data structure only, not derived
// files that USE the data structure
#ifndef __QUADRIC_POINT_INTEGRATOR__
#define __QUADRIC_POINT_INTEGRATOR__

#include "gaudi/common.h"
#include "rod_integrators.hpp"
#include "shell_integrators.hpp"
#include <cmath>
#include <cstdlib>
#include <ostream>
#include <vector>

namespace gaudi {

namespace calder {

std::vector<vec10> quadric(asawa::rod::rod &R, const std::vector<vec3> &Nr,
                           const std::vector<vec3> &p_pov,
                           const std::vector<vec3> &N_pov, real l0,
                           real p = 3.0, real w0 = 1e-3) {
  std::vector<real> weights = R.l0();

  auto mkA = [](vec3 dx) {
    mat410 A;
    real x = dx[0];
    real y = dx[1];
    real z = dx[2];
    real xx = x * x;
    real yy = y * y;
    real zz = z * z;
    real xy = 2.0 * x * y;
    real xz = 2.0 * x * z;
    real yz = 2.0 * y * z;

    A.col(0) = vec4(xx, 2.0 * x, 0.0 * y, 0.0 * z);
    A.col(1) = vec4(yy, 0.0 * x, 2.0 * y, 0.0 * z);
    A.col(2) = vec4(zz, 0.0 * x, 0.0 * y, 2.0 * z);

    A.col(3) = vec4(xy, 2.0 * y, 2.0 * x, 0.0 * z);
    A.col(4) = vec4(xz, 2.0 * z, 0.0 * y, 2.0 * x);
    A.col(5) = vec4(yz, 0.0 * x, 2.0 * z, 2.0 * y);

    A.col(6) = vec4(2.0 * x, 2.0, 0.0, 0.0);
    A.col(7) = vec4(2.0 * y, 0.0, 2.0, 0.0);
    A.col(8) = vec4(2.0 * z, 0.0, 0.0, 2.0);

    A.col(9) = vec4(1.0, 0.0, 0.0, 0.0);
    return A;
  };

  auto mkN = [](vec3 N) { return vec4(0.0, N[0], N[1], N[2]); };
  std::vector<mat10> As(p_pov.size(), mat10::Zero());
  std::vector<vec10> ns(p_pov.size(), vec10::Zero());

  for (int i = 0; i < p_pov.size(); i++) {
    mat410 A = mkA(vec3::Zero());
    vec4 n = mkN(N_pov[i]);
    ns[i] += w0 * A.transpose() * n;
    As[i] += w0 * A.transpose() * A;
  }

  std::vector<real> us = integrate_over_rod<real>(
      R, p_pov,
      [&Nr](const std::vector<index_t> &edge_ids, Rod_Sum_Type &sum) {
        sum.bind(calder::vec3_datum::create(edge_ids, Nr));
      },
      [&](const index_t i, const index_t j, //
          const vec3 &pi, const vec3 &pj,
          const std::vector<calder::datum::ptr> &data,
          Rod_Sum_Type::Node_Type node_type, //
          const Rod_Sum_Type::Node &node,    //
          const Rod_Sum_Type::Tree &tree) -> real {
        const calder::vec3_datum::ptr F_datum =
            static_pointer_cast<calder::vec3_datum>(data[0]);
        vec3 Nr = get_data<vec3>(node_type, j, 0, data);
        real w = Nr.norm();
        Nr /= w;

        vec3 dp = pj - pi;
        real dist = dp.norm();
        real kappa = computeK(dist, l0, p);
        // real kappa = computeKg(dist, l0, p);

        vec3 dpN = dp / dist;

        real NtN1 = 1.0 + Nr.dot(N_pov[i]);
        NtN1 = pow(0.5 * NtN1, 0.5);
        real Ndp1 = 1.0 + Nr.dot(dpN);
        Ndp1 = pow(0.5 * Ndp1, 0.5);
        real NNdp = 1.0 - dp.dot(Nr.cross(N_pov[i]));
#if 1
        real Nrdp = Nr.dot(dp);
        real NtN = Nr.dot(N_pov[i]);
        if (Nrdp < 0 && dist < 5.0 * l0) {
          // quadric doesn't compute torus, so we reflect the normal
          // to make it look more like a torus... we convexify the point
          vec3 pj_t = pj - 0.5 * Nrdp * Nr;
          dp = pj_t - pi;
          Nr *= -1.0;
          dp = vec3::Zero();
        }
#endif

        // real w_total = w * kappa;
        real w_total = w * kappa;

        mat410 A = mkA(dp);
        vec4 n = mkN(Nr);
        ns[i] += w_total * A.transpose() * n;
        As[i] += w_total * A.transpose() * A;
        return 0.0;
      });

  std::vector<vec10> out(p_pov.size(), vec10::Zero());
  for (int i = 0; i < As.size(); i++) {
    mat10 &A = As[i];
    vec10 &n = ns[i];
    out[i] = A.colPivHouseholderQr().solve(n);
  }

  return out;
}

std::vector<vec3> quadric_grad(asawa::rod::rod &R, const std::vector<vec3> &Nr,
                               const std::vector<vec3> &p_pov,
                               const std::vector<vec3> &N_pov, real l0,
                               real p = 3.0, bool normalize = true) {
  std::vector<vec10> Q = quadric(R, Nr, p_pov, N_pov, l0, p, normalize);
  std::vector<vec3> out(p_pov.size(), vec3::Zero());
  for (int i = 0; i < Q.size(); i++) {
    out[i] = -Q[i][9] * Q[i].segment(6, 3);
    if (out[i].hasNaN()) {
      out[i] = vec3::Zero();
    }
  }
  return out;
}

std::vector<real> quadric_sdf(asawa::rod::rod &R, const std::vector<vec3> &Nr,
                              const std::vector<vec3> &p_pov,
                              const std::vector<vec3> &N_pov, real l0,
                              real p = 3.0, bool normalize = true) {
  std::vector<vec10> Q = quadric(R, Nr, p_pov, N_pov, l0, p, normalize);
  std::vector<real> out(p_pov.size(), 0.0);
  for (int i = 0; i < Q.size(); i++) {
    out[i] = Q[i][9]; //* q.segment(6, 3);
    if (std::isnan(out[i])) {
      out[i] = 0.0;
    }
  }
  return out;
}

} // namespace calder
} // namespace gaudi
#endif