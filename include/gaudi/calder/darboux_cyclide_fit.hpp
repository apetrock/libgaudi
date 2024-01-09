//
//  m2Includes.h
//  Manifold
//
//  Created by John Delaney on 5/22/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//
// includes files that are required for the data structure only, not derived
// files that USE the data structure
#ifndef __DARBOUX_CYCLIDE_FIT__
#define __DARBOUX_CYCLIDE_FIT__

#include "gaudi/common.h"
#include "rod_integrators.hpp"
#include "shell_integrators.hpp"
#include <cmath>
#include <cstdlib>
#include <ostream>
#include <vector>

namespace gaudi {

namespace calder {

TYPEDEF_VEC(14)
TYPEDEF_MAT(14)
TYPEDEF_MAT_NM(4, 14)

std::vector<vec14> darboux_cyclide(asawa::rod::rod &R,
                                   const std::vector<vec3> &Nr,
                                   const std::vector<vec3> &p_pov,
                                   const std::vector<vec3> &N_pov, real l0,
                                   real p = 3.0, real w0 = 1e-2) {
  std::vector<real> weights = R.l0();

  auto mkA = [](vec3 dx) {
    mat414 A;
    real x = dx[0];
    real y = dx[1];
    real z = dx[2];

    real xx = x * x;
    real yy = y * y;
    real zz = z * z;
    real xy = 2.0 * x * y;
    real xz = 2.0 * x * z;
    real yz = 2.0 * y * z;

    real xyz = xx + yy + zz;
    real xyz2 = xyz * xyz;
    real xxyz = x * xyz;
    real yxyz = y * xyz;
    real zxyz = z * xyz;
    A.col(0) = vec4(xx, 2.0 * x, 0.0 * y, 0.0 * z); // A
    A.col(1) = vec4(yy, 0.0 * x, 2.0 * y, 0.0 * z); // B
    A.col(2) = vec4(zz, 0.0 * x, 0.0 * y, 2.0 * z); // C

    A.col(3) = vec4(xy, 2.0 * y, 2.0 * x, 0.0 * z); // D
    A.col(4) = vec4(xz, 2.0 * z, 0.0 * y, 2.0 * x); // E
    A.col(5) = vec4(yz, 0.0 * x, 2.0 * z, 2.0 * y); // F

    A.col(6) = vec4(2.0 * x, 2.0, 0.0, 0.0); // G
    A.col(7) = vec4(2.0 * y, 0.0, 2.0, 0.0); // H
    A.col(8) = vec4(2.0 * z, 0.0, 0.0, 2.0); // I

    A.col(9) = vec4(1.0, 0.0, 0.0, 0.0); // J

    A.col(10) = vec4(xyz2, 4.0 * xxyz, 4.0 * yxyz, 4.0 * zxyz);    // lambda
    A.col(11) = vec4(x * xyz, 2.0 * xx + xyz, 2.0 * xy, 2.0 * xz); // mu
    A.col(12) = vec4(y * xyz, 2.0 * xy, 2.0 * yy + xyz, 2.0 * yz); // nu
    A.col(13) = vec4(z * xyz, 2.0 * xz, 2.0 * yz, 2.0 * zz + xyz); // kappa

    return A;
  };

  auto mkN = [](vec3 N) { return vec4(0.0, N[0], N[1], N[2]); };
  std::vector<mat14> As(p_pov.size(), mat14::Zero());
  std::vector<vec14> ns(p_pov.size(), vec14::Zero());

  for (int i = 0; i < p_pov.size(); i++) {
    mat414 A = mkA(vec3::Zero());
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
        vec3 Ns = N_pov[i];
        real w = Nr.norm();
        Nr /= w;

        vec3 dp = pj - pi;
        real dist = dp.norm();
        real kappa = computeK(dist, l0, p);
        // real kappa = computeKg(dist, l0, p);
        real w_total = w * kappa;

        // Nr = dp.normalized();
        mat414 A = mkA(dp);
        vec4 n = mkN(Nr);

        ns[i] += w_total * A.transpose() * n;
        As[i] += w_total * A.transpose() * A;
        return 0.0;
      });

  std::vector<vec14> out(p_pov.size(), vec14::Zero());
  for (int i = 0; i < As.size(); i++) {
    mat14 &A = As[i];
    vec14 &n = ns[i];
    out[i] = A.colPivHouseholderQr().solve(n);
  }

  return out;
}

std::vector<vec3> darboux_cyclide_grad(asawa::rod::rod &R,
                                       const std::vector<vec3> &Nr,
                                       const std::vector<vec3> &p_pov,
                                       const std::vector<vec3> &N_pov, real l0,
                                       real p = 3.0) {
  std::vector<vec14> Q = darboux_cyclide(R, Nr, p_pov, N_pov, l0, p, 1e-3);
  std::vector<vec3> out(p_pov.size(), vec3::Zero());
  for (int i = 0; i < Q.size(); i++) {
    out[i] = -Q[i][9] * Q[i].segment(6, 3);
    if (out[i].hasNaN()) {
      out[i] = vec3::Zero();
    }
  }
  return out;
}

std::vector<real> darboux_cyclide_sdf(asawa::rod::rod &R,
                                      const std::vector<vec3> &Nr,
                                      const std::vector<vec3> &p_pov,
                                      const std::vector<vec3> &N_pov, real l0,
                                      real p = 3.0) {
  std::vector<vec14> Q = darboux_cyclide(R, Nr, p_pov, N_pov, l0, p, 1e-4);
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