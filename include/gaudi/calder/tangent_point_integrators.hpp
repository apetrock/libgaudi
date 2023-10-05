//
//  m2Includes.h
//  Manifold
//
//  Created by John Delaney on 5/22/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//
// includes files that are required for the data structure only, not derived
// files that USE the data structure
#ifndef __TANGENT_POINT_INTEGRATOR__
#define __TANGENT_POINT_INTEGRATOR__

#include "gaudi/logger.hpp"
#include "integrators.hpp"
#include "rod_integrators.hpp"
#include "shell_integrators.hpp"
#include <cmath>
#include <cstdlib>
#include <ostream>
#include <vector>

namespace gaudi {

namespace calder {

real compute_tangent_point_radius(const vec3 &dp, const vec3 &N) {
  real ndp = dp.squaredNorm();
  real nPdp = (N * N.transpose() * dp).norm();
  return 0.5 * ndp / nPdp;
};

real compute_tangent_point_inverse_radius(const vec3 &dp, const vec3 &N,
                                          const real &l0, const double &p) {
  // computes the inverse radius of the tangent point to power p
  real ndp = dp.norm();
  real nPdp = (N * N.transpose() * dp).norm();
  real lp = pow(l0, p);
  real k = pow(nPdp, p) / (pow(ndp, 2.0 * p) + lp);
  return k;
};

vec3 compute_tangent_point_radius_grad_0(const vec3 &dp0, const vec3 &N,
                                         const real &l0, const real &p) {
  real Ndp = dp0.dot(N);

  real ndp = dp0.norm();
  vec3 dp = ndp * dp0.normalized();

  mat3 P = N * N.transpose();
  vec3 Pdp = P * dp;

  real nPdp = (Pdp).norm();
  real lp = pow(l0, p);
  real l2 = pow(l0, 2.0);

  real k = pow(nPdp, p) / (pow(ndp, 2.0 * p) + lp);
  vec3 dk0 = P * Pdp / (nPdp * nPdp + l2);
  vec3 dk1 = 2.0 * dp / (ndp * ndp + l2);
  vec3 dk = p * k * (dk0 - dk1);
  // vec3 dk = p * k * (dk0);

  return dk;
};

vec3 compute_tangent_point_radius_grad_1(const vec3 &dp0, const vec3 &N,
                                         const real &l0, const real &p) {
  real Ndp = dp0.dot(N);

  real fx = dp0.norm();
  vec3 dfx = fx * dp0.normalized();

  mat3 P = N * N.transpose();
  vec3 Pdp = P * dp0;

  real fPx = (Pdp).norm();
  vec3 dfPx = fPx * Pdp.normalized();

  real lp = pow(l0, p);

  real denom = (pow(fx, 2.0 * p) + lp);
  real k = pow(fPx, p) / denom;
  vec3 dk0 = P * pow(fPx, p - 1.0) * dfPx / denom;
  vec3 dk1 = 2.0 * k * pow(fx, 2.0 * p - 1.0) * dfx / denom;
  vec3 dk = p * (dk0 - dk1);
  // vec3 dk = p * k * (dk0);

  return -dk;
};

vec3 compute_tangent_point_radius_grad(const vec3 &dp, const vec3 &N,
                                       const real &l0, const real &p) {

  real fx = dp.norm();
  mat3 P = N * N.transpose();
  vec3 Px = P * dp;

  real fPx = Px.norm();
  // vec3 dfPx = fPx * Pdp.normalized();

  real lp = pow(l0, p);
  real l2 = pow(l0, 2.0);

  real k = pow(fPx, p) / (pow(fx, 2.0 * p) + lp);

  //  P = N*Nt => P * P = N*Nt*N*Nt = N*Nt = P
  vec3 dk = p * k * (Px / (Px.dot(Px) + l2) - 2.0 * dp / (dp.dot(dp) + l2));

  // vec3 dk = p * k * (dk0);

  return dk;
};

std::vector<vec3> tangent_point_gradient(asawa::rod::rod &R,
                                         const std::vector<vec3> &p_pov,
                                         const std::vector<real> &w_pov,
                                         const std::vector<vec3> &T_pov,
                                         real l0, real p = 3.0) {
  std::vector<real> weights = R.l0();
  std::vector<vec3> Tc = R.N2c();

  std::vector<vec3> us = integrate_over_rod<vec3>(
      R, p_pov,
      [&weights, &Tc](const std::vector<index_t> &edge_ids, Rod_Sum_Type &sum) {
        sum.bind(calder::scalar_datum::create(edge_ids, weights));
        sum.bind(calder::vec3_datum::create(edge_ids, Tc));
      },
      [l0, &w_pov, &T_pov, p](const index_t i, const index_t j, //
                              const vec3 &pi, const vec3 &pj,
                              const std::vector<calder::datum::ptr> &data,
                              Rod_Sum_Type::Node_Type node_type, //
                              const Rod_Sum_Type::Node &node,    //
                              const Rod_Sum_Type::Tree &tree) -> vec3 {
        real wi = w_pov[i];
        vec3 Ti = T_pov[i];
        real wj = get_data<real>(node_type, j, 0, data);
        vec3 Tj = get_data<vec3>(node_type, j, 1, data);
        Ti.normalize();
        Tj.normalize();

        vec3 dp = pj - pi;
        // vec3 Ni = va::rejection_matrix(Ti) * dp;
        // vec3 Nj = va::rejection_matrix(Tj) * dp;
        vec3 Bi = Ti.cross(dp).normalized();
        vec3 Bj = Tj.cross(dp).normalized();
        vec3 Ni = Bi.cross(Ti).normalized();
        vec3 Nj = Bj.cross(Tj).normalized();

        vec3 gj = compute_tangent_point_radius_grad(dp, Nj, l0, p);
        vec3 gi = compute_tangent_point_radius_grad(-dp, Ni, l0, p);

        return 0.5 * (wi * gi - wj * gj);
      });
  return us;
}

std::vector<real> tangent_point_energy(asawa::shell::shell &M,
                                       const std::vector<vec3> &p_pov,
                                       const std::vector<real> &w_pov,
                                       const std::vector<vec3> &N_pov, //
                                       real l0, real p = 3.0) {
  std::vector<vec3> x = asawa::get_vec_data(M, 0);
  std::vector<real> weights = asawa::shell::face_areas(M, x);
  std::vector<vec3> Nc = asawa::shell::face_normals(M, x);

  std::vector<real> us = integrate_over_shell<real>(
      M, p_pov,
      [&weights, &Nc](const std::vector<index_t> &face_ids,
                      Shell_Sum_Type &sum) {
        sum.bind(calder::scalar_datum::create(face_ids, weights));
        sum.bind(calder::vec3_datum::create(face_ids, Nc));
      },
      [l0, &w_pov, &N_pov, p](const index_t i, const index_t j, //
                              const vec3 &pi, const vec3 &pj,
                              const std::vector<calder::datum::ptr> &data,
                              Shell_Sum_Type::Node_Type node_type, //
                              const Shell_Sum_Type::Node &node,    //
                              const Shell_Sum_Type::Tree &tree) -> real {
        real wi = w_pov[i];
        vec3 Ni = N_pov[i];

        real wj = get_data<real>(node_type, j, 0, data);
        vec3 Nj = get_data<vec3>(node_type, j, 1, data);

        wi = std::max(wi, 1e-6);
        wj = std::max(wj, 1e-6);

        Ni.normalize();
        Nj.normalize();

        vec3 dp = pj - pi;

        real jr = compute_tangent_point_inverse_radius(dp, Nj, l0, p);
        real ir = compute_tangent_point_inverse_radius(-dp, Ni, l0, p);
        // not sure if symmetric?
        return wj * jr;
      });
  return us;
}

std::vector<vec3> tangent_point_gradient(asawa::shell::shell &M,
                                         const std::vector<vec3> &p_pov,
                                         const std::vector<real> &w_pov,
                                         const std::vector<vec3> &N_pov,
                                         real l0, real p = 3.0) {
  std::vector<vec3> x = asawa::get_vec_data(M, 0);
  std::vector<real> weights = asawa::shell::face_areas(M, x);
  std::vector<vec3> Nc = asawa::shell::face_normals(M, x);

  std::vector<vec3> us = integrate_over_shell<vec3>(
      M, p_pov,
      [&weights, &Nc](const std::vector<index_t> &face_ids,
                      Shell_Sum_Type &sum) {
        sum.bind(calder::scalar_datum::create(face_ids, weights));
        sum.bind(calder::vec3_datum::create(face_ids, Nc));
      },
      [l0, &w_pov, &N_pov, p](const index_t i, const index_t j, //
                              const vec3 &pi, const vec3 &pj,
                              const std::vector<calder::datum::ptr> &data,
                              Shell_Sum_Type::Node_Type node_type, //
                              const Shell_Sum_Type::Node &node,    //
                              const Shell_Sum_Type::Tree &tree) -> vec3 {
        real wi = w_pov[i];
        vec3 Ni = N_pov[i];
        real wj = get_data<real>(node_type, j, 0, data);
        vec3 Nj = get_data<vec3>(node_type, j, 1, data);
        wi = std::max(wi, 1e-6);
        wj = std::max(wj, 1e-6);

        vec3 dp = pj - pi;
        // if (i == 0) {
        //   logger::line(pi, pj, vec4(0.0, 1.0, 0.5, 1.0));
        //   logger::line(pj, pj + 0.1 * Nj, vec4(0.0, 1.0, 0.5, 1.0));
        // }
        Ni.normalize();
        Nj.normalize();
        real pk = p;

        vec3 gj = compute_tangent_point_radius_grad(dp, Nj, l0, pk);
        vec3 gi = compute_tangent_point_radius_grad(-dp, Ni, l0, pk);

        return wi * gi - wj * gj;
        // return wi * gi;

        // return -wj * gj;
      });
  return us;
}

std::vector<mat3> tangent_point_gradient_frame(asawa::shell::shell &M,
                                               const std::vector<vec3> &p_pov,
                                               const std::vector<real> &w_pov,
                                               const std::vector<vec3> &N_pov,
                                               real l0, real p = 3.0) {
  std::vector<vec3> x = asawa::get_vec_data(M, 0);
  std::vector<real> weights = asawa::shell::face_areas(M, x);
  std::vector<vec3> Nc = asawa::shell::face_normals(M, x);

  std::vector<real> sums(x.size(), 0.0);
  std::vector<mat3> us = integrate_over_shell<mat3>(
      M, p_pov,
      [&weights, &Nc](const std::vector<index_t> &face_ids,
                      Shell_Sum_Type &sum) {
        sum.bind(calder::scalar_datum::create(face_ids, weights));
        sum.bind(calder::vec3_datum::create(face_ids, Nc));
      },
      [l0, &w_pov, &N_pov, p](const index_t i, const index_t j, //
                              const vec3 &pi, const vec3 &pj,
                              const std::vector<calder::datum::ptr> &data,
                              Shell_Sum_Type::Node_Type node_type, //
                              const Shell_Sum_Type::Node &node,    //
                              const Shell_Sum_Type::Tree &tree) -> mat3 {
        real wi = w_pov[i];
        vec3 Ni = N_pov[i];
        real wj = get_data<real>(node_type, j, 0, data);
        vec3 Nj = get_data<vec3>(node_type, j, 1, data);
        wi = std::max(wi, 1e-6);
        wj = std::max(wj, 1e-6);

        vec3 dp = pj - pi;
        // if (i == 0) {
        //   logger::line(pi, pj, vec4(0.0, 1.0, 0.5, 1.0));
        //   logger::line(pj, pj + 0.1 * Nj, vec4(0.0, 1.0, 0.5, 1.0));
        // }
        Ni.normalize();
        Nj.normalize();
        real pk = p;

        vec3 gj = compute_tangent_point_radius_grad(dp, Nj, l0, pk);
        vec3 gi = compute_tangent_point_radius_grad(-dp, Ni, l0, pk);
        vec3 wg = wi * gi - wj * gj;
        // sums[i] += wi + wj;
        return wg * wg.transpose();
        // return wi * gi;

        // return -wj * gj;
      });
  /*
  for (int i = 0; i < us.size(); i++) {
    us[i] /= sums[i];
  }
  */
  return us;
}

} // namespace calder
} // namespace gaudi
#endif