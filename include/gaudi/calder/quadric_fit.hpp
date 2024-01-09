//
//  m2Includes.h
//  Manifold
//
//  Created by John Delaney on 5/22/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//
// includes files that are required for the data structure only, not derived
// files that USE the data structure
#ifndef __QUADRIC_FIT__
#define __QUADRIC_FIT__

#include "gaudi/common.h"
#include "integrators.hpp"
#include "rod_integrators.hpp"
#include "shell_integrators.hpp"
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <ostream>
#include <unistd.h>
#include <utility>
#include <vector>

/*
TODO: Merge this file with the darboux file, rename regressions.

*/
namespace gaudi {

namespace calder {

TYPEDEF_VEC(7)
TYPEDEF_VEC(10)
TYPEDEF_MAT(10)
TYPEDEF_MAT_NM(4, 10)

vec3 quadric_grad(const vec10 &Q, const vec3 x) {

  vec3 g = {
      2.0 * Q[0] * x[0] + 2.0 * Q[3] * x[1] + 2.0 * Q[4] * x[2] + 2.0 * Q[6], //
      2.0 * Q[3] * x[0] + 2.0 * Q[1] * x[1] + 2.0 * Q[5] * x[2] + 2.0 * Q[7], //
      2.0 * Q[4] * x[0] + 2.0 * Q[5] * x[1] + 2.0 * Q[2] * x[2] + 2.0 * Q[8]  //
  };
  return g;
}

mat410 mk_quad_A(vec3 dx) {
  mat410 A;
  real x = dx[0];
  real y = dx[1];
  real z = dx[2];
  real x2 = 2.0 * x;
  real y2 = 2.0 * y;
  real z2 = 2.0 * z;

  real xx = x * x;
  real yy = y * y;
  real zz = z * z;
  real xy = 2.0 * x * y;
  real xz = 2.0 * x * z;
  real yz = 2.0 * y * z;

  A.row(0) = vec10(xx, yy, zz, xy, xz, yz, x2, y2, z2, 1.0);
  // A = X^2 * I + sort of Skew(X) + 1
  A.row(1) = vec10(x2, 0.0, 0.0, y2, z2, 0.0, 2.0, 0.0, 0.0, 0.0);
  A.row(2) = vec10(0.0, y2, 0.0, x2, 0.0, z2, 0.0, 2.0, 0.0, 0.0);
  A.row(3) = vec10(0.0, 0.0, z2, 0.0, x2, y2, 0.0, 0.0, 2.0, 0.0);

  return A;
}

vec4 mk_quad_N(const vec3 &N) { return vec4(0.0, N[0], N[1], N[2]); }

mat3 quadric_hessian(const vec10 &Q) {
  mat3 A = mat3::Zero();
  A.row(0) = vec3(Q[0], Q[3], Q[4]);
  A.row(1) = vec3(Q[3], Q[1], Q[5]);
  A.row(2) = vec3(Q[4], Q[5], Q[2]);
  return A;
}

vec3 quadric_center(const vec10 &Q) {
  mat3 A = quadric_hessian(Q);
  vec3 b = -vec3(Q[6], Q[7], Q[8]);
  vec3 c = A.colPivHouseholderQr().solve(b);
  Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(A);
  int rank = lu_decomp.rank();
  if (rank < 3) {
    std::cout << "rank: " << rank << std::endl;
  }
  // vec3 c = -A.inverse() * b;
#if 0
  real det = A.determinant();
  Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(A);
  int rank = lu_decomp.rank();
  std::cout << "Q: " << Q.transpose() << std::endl;
  std::cout << "rank: " << rank << std::endl;
  std::cout << "det: " << det << std::endl;

#endif
  std::cout << "calc'd grad: " << quadric_grad(Q, c).norm() << std::endl;
  return c;
}

std::vector<vec10> quadric(asawa::rod::rod &R, const std::vector<vec3> &Nr,
                           const std::vector<vec3> &p_pov,
                           const std::vector<vec3> &N_pov, real l0,
                           real p = 3.0, real w0 = 1e-2, real w_r = 6.0,
                           real k_r = 0.5) {
  std::vector<real> weights = R.l0();

  std::vector<mat10> As(p_pov.size(), mat10::Zero());
  std::vector<vec10> ns(p_pov.size(), vec10::Zero());

  for (int i = 0; i < p_pov.size(); i++) {
    mat410 A = mk_quad_A(vec3::Zero());
    vec4 n = mk_quad_N(N_pov[i]);
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

        vec3 dpN = dp / dist;

        real NtN1 = 1.0 + Nr.dot(Ns);
        NtN1 = pow(0.5 * NtN1, 0.5);
        real Ndp1 = 1.0 + Nr.dot(dpN);
        Ndp1 = pow(0.5 * Ndp1, 0.5);
        real NNdp = 1.0 - dp.dot(Nr.cross(Ns));
#if 1
        real Nrdp = Nr.dot(dpN);
        real Nsdp = Ns.dot(dpN);
        real NtN = Nr.dot(Ns);
        if (Nrdp < 0 && Nsdp < 0 && dist < w_r * l0) {

          // quadric doesn't compute torus, so we reflect the normal
          // to make it look more like a torus... we convexify the point
          // we vary the convexification based on the distance from the
          // edge of the sphere.
          real t_r = (w_r * l0 - dist) / dist;
          vec3 pj_t = pj + t_r * k_r * Nrdp * Nsdp * dist * Nr;
          dp = pj_t - pi;
          Nr *= -1.0;
          dp = vec3::Zero();
        }
#endif

        // real w_total = w * kappa;
        real w_total = w * kappa;

        mat410 A = mk_quad_A(dp);
        vec4 n = mk_quad_N(Nr);
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
                               real p = 3.0, real w_r = 6.0, real k_r = 0.5) {
  std::vector<vec10> Q = quadric(R, Nr, p_pov, N_pov, l0, p, 1e-5, w_r, k_r);
  std::vector<vec3> out(p_pov.size(), vec3::Zero());
  for (int i = 0; i < Q.size(); i++) {
    out[i] = quadric_grad(Q[i], p_pov[i]);
    if (out[i].hasNaN()) {
      out[i] = vec3::Zero();
    }
  }
  return out;
}

std::vector<real> quadric_sdf(asawa::rod::rod &R, const std::vector<vec3> &Nr,
                              const std::vector<vec3> &p_pov,
                              const std::vector<vec3> &N_pov, real l0,
                              real p = 3.0, real w_r = 6.0, real k_r = 0.5) {
  std::vector<vec10> Q = quadric(R, Nr, p_pov, N_pov, l0, p, 1e-2, w_r, k_r);
  std::vector<real> out(p_pov.size(), 0.0);
  for (int i = 0; i < Q.size(); i++) {
    out[i] = Q[i][9]; //* q.segment(6, 3);
    if (std::isnan(out[i])) {
      out[i] = 0.0;
    }
  }
  return out;
}

std::vector<vec10> quadric(asawa::shell::shell &M,
                           const std::vector<vec3> &p_pov,
                           const std::vector<vec3> &n_pov, real l0,
                           real p = 3.0) {

  const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);

  std::vector<real> weights = asawa::shell::face_areas(M, x);
  std::vector<vec3> Ns = asawa::shell::face_normals(M, x);

  for (int i = 0; i < Ns.size(); i++) {
    Ns[i] = weights[i] * Ns[i];
  }

  std::vector<mat10> As(p_pov.size(), mat10::Zero());
  std::vector<vec10> ns(p_pov.size(), vec10::Zero());
#if 0
  real w0 = 0.1;
  for (int i = 0; i < p_pov.size(); i++) {
    mat410 A = mk_quad_A(vec3::Zero());
    vec4 n = mkN(n_pov[i]);
    ns[i] += w0 * A.transpose() * n;
    As[i] += w0 * A.transpose() * A;
  }
#endif

  std::vector<real> us = integrate_over_shell<real>(
      M, p_pov,
      [&Ns](const std::vector<index_t> &edge_ids, Shell_Sum_Type &sum) {
        sum.bind(calder::vec3_datum::create(edge_ids, Ns));
      },
      [&](const index_t i, const index_t j, //
          const vec3 &pi, const vec3 &pj,
          const std::vector<calder::datum::ptr> &data,
          Shell_Sum_Type::Node_Type node_type, //
          const Shell_Sum_Type::Node &node,    //
          const Shell_Sum_Type::Tree &tree) -> real {
        const calder::vec3_datum::ptr F_datum =
            static_pointer_cast<calder::vec3_datum>(data[0]);
        /*
        if (node_type == Shell_Sum_Type::Node_Type::BRANCH) {
          return 0.0;
        }
        */
        vec3 Nj = get_data<vec3>(node_type, j, 0, data);
        vec3 Ni = n_pov[i];

        real wj = Nj.norm();
        if (wj < 1e-8) {
          return 0.0;
        }
        Nj /= wj;

        vec3 dp = pj - pi;
        vec3 dpN = dp.normalized();
        real dpdN = dpN.dot(Ni);
        real NdN = Nj.dot(Ni);

        if (NdN < 0) {
          return 0.0;
        }

        real dist = dp.norm();
        /*
        if (i == 1000) {
          logger::line(pi, pj, vec4(0.0, 1.0, 0.0, 1.0));
          logger::line(pj, pj + 0.1 * Nj, vec4(0.0, 0.0, 1.0, 1.0));
        }
*/
        real kappa = computeK(dist, l0, p);
        // real kappa = computeKg(dist, l0, 2.0);

        // real w_total = w * kappa;
        real w_total = wj * kappa * NdN;

        mat410 A = mk_quad_A(dp);
        vec4 n = mk_quad_N(Nj);
        ns[i] += w_total * A.transpose() * n;
        As[i] += w_total * A.transpose() * A;
        return 0.0;
      });

  std::vector<vec10> out(p_pov.size(), vec10::Zero());
  for (int i = 0; i < As.size(); i++) {
    const mat10 &A = As[i];
    const vec10 &n = ns[i];
    out[i] = A.colPivHouseholderQr().solve(n);
  }

  return out;
}

std::vector<std::pair<vec3, mat3>>
quadric_hessian(asawa::shell::shell &M, const std::vector<vec3> &p_pov,
                const std::vector<vec3> &n_pov, real l0, real p = 3.0) {
  std::vector<vec10> Q = quadric(M, p_pov, n_pov, l0, p);
  auto z = std::make_pair(vec3::Zero(), mat3::Zero());
  std::vector<std::pair<vec3, mat3>> out(p_pov.size(), z);
  for (int i = 0; i < Q.size(); i++) {
    vec3 x = p_pov[i];
    vec10 Qi = Q[i];
    vec3 g = quadric_grad(Q[i], p_pov[i]);

    // vec3 Ni = g.normalized();
    vec3 Ni = n_pov[i];
    // logger::line(x, x + 0.05 * Ni, vec4(1.0, 0.0, 1.0, 1.0));
    mat3 P = Ni * Ni.transpose();
    mat3 R = va::rejection_matrix(Ni);
    out[i] = std::make_pair(g, quadric_hessian(Q[i]));
  }
  return out;
}

std::vector<mat3> quadric_curvature(asawa::shell::shell &M,
                                    const std::vector<vec3> &p_pov,
                                    const std::vector<vec3> &n_pov, real l0,
                                    real p = 3.0) {
  std::vector<std::pair<vec3, mat3>> H =
      quadric_hessian(M, p_pov, n_pov, l0, p);
  std::vector<mat3> out(p_pov.size(), mat3::Zero());
  for (int i = 0; i < H.size(); i++) {
    vec3 x = p_pov[i];
    mat3 Hi = H[i].second;
    vec3 Ni = n_pov[i];
    // logger::line(x, x + 0.05 * Ni, vec4(1.0, 0.0, 1.0, 1.0));
    mat3 P = Ni * Ni.transpose();
    mat3 R = va::rejection_matrix(Ni);

    mat3 RH = R * Hi;
    vec3 c0 = RH.col(0);
    vec3 c1 = RH.col(1);
    vec3 c2 = RH.col(2);

    Eigen::JacobiSVD<mat3> svd(RH, Eigen::ComputeFullU | Eigen::ComputeFullV);
    mat3 U = svd.matrixU();
    mat3 V = svd.matrixV();
    vec3 s = svd.singularValues();
    mat3 S = mat3::Zero();
    S.col(0) = s[0] * U.col(0);
    S.col(1) = s[1] * U.col(1);
    S.col(2) = s[2] * U.col(2);
    if (s.hasNaN()) {
      continue;
    }
    /*
    logger::line(x, x + 0.01 * S.col(0), vec4(1.0, 0.0, 0.0, 1.0));
    logger::line(x, x + 0.01 * S.col(1), vec4(0.0, 1.0, 0.0, 1.0));
    logger::line(x, x + 0.01 * S.col(2), vec4(0.0, 0.0, 1.0, 1.0));
*/

    out[i] = S;
  }
  return out;
}

std::vector<vec3> quadric_center(asawa::shell::shell &M,
                                 const std::vector<vec3> &p_pov,
                                 const std::vector<vec3> &n_pov, real l0,
                                 real p = 3.0) {
  std::vector<vec10> Q = quadric(M, p_pov, n_pov, l0, p);
  std::vector<vec3> out(p_pov.size(), vec3::Zero());
  for (int i = 0; i < Q.size(); i++) {
    vec3 x = p_pov[i];
    vec10 Qi = Q[i];
    vec3 cen = quadric_center(Q[i]);
    vec3 g = quadric_grad(Q[i], p_pov[i]);

    out[i] = cen;
  }
  return out;
}

TYPEDEF_VEC(18)
vec18 mk_sphere_v(const vec3 &xi) {
  vec18 S = vec18::Zero();
  real x = xi[0];
  real y = xi[1];
  real z = xi[2];
  S[0] = x;
  S[1] = y;
  S[2] = z;
  S[3] = x * x, S[4] = y * y, S[5] = z * z;
  S[6] = x * y, S[7] = x * z, S[8] = y * z;
  S[9] = x * x * x, S[10] = x * x * y, S[11] = x * x * z;
  S[12] = y * y * x, S[13] = y * y * y, S[14] = y * y * z;
  S[15] = z * z * x, S[16] = z * z * y, S[17] = z * z * z;
  return S;
}

vec4 calc_sphere(real W, const vec18 &S) {
  real Sx = S[0], Sy = S[1], Sz = S[2];
  real Sxx = S[3], Syy = S[4], Szz = S[5];
  real Sxy = S[6], Sxz = S[7], Syz = S[8];
  real Sxxx = S[9], Sxxy = S[10], Sxxz = S[11];
  real Syyx = S[12], Syyy = S[13], Syyz = S[14];
  real Szzx = S[15], Szzy = S[16], Szzz = S[17];

  real A1 = Sxx + Syy + Szz;
  real a = 2.0 * Sx * Sx - 2.0 * W * Sxx;
  real b = 2.0 * Sx * Sy - 2.0 * W * Sxy;
  real c = 2.0 * Sx * Sz - 2.0 * W * Sxz;
  real d = -W * (Sxxx + Syyx + Szzx) + A1 * Sx;

  real e = 2.0 * Sx * Sy - 2.0 * W * Sxy;
  real f = 2.0 * Sy * Sy - 2.0 * W * Syy;
  real g = 2.0 * Sy * Sz - 2.0 * W * Syz;
  real h = -W * (Sxxy + Syyy + Szzy) + A1 * Sy;
  real j = 2.0 * Sx * Sz - 2.0 * W * Sxz;
  real k = 2.0 * Sy * Sz - 2.0 * W * Syz;
  real l = 2.0 * Sz * Sz - 2.0 * W * Szz;
  real m = -W * (Sxxz + Syyz + Szzz) + A1 * Sz;

  real delta = a * (f * l - g * k) - e * (b * l - c * k) + j * (b * g - c * f);
  vec3 X = vec3::Zero();
  X[0] =
      (d * (f * l - g * k) - h * (b * l - c * k) + m * (b * g - c * f)) / delta;
  X[1] =
      (a * (h * l - m * g) - e * (d * l - m * c) + j * (d * g - h * c)) / delta;
  X[2] =
      (a * (f * m - h * k) - e * (b * m - d * k) + j * (b * h - d * f)) / delta;

  real R = sqrt(X[0] * X[0] + X[1] * X[1] + X[2] * X[2] +
                (A1 - 2.0 * (X[0] * Sx + X[1] * Sy + X[2] * Sz)) / W);
  return vec4(X[0], X[1], X[2], R);
}

std::vector<vec4> sphere(asawa::shell::shell &M, const std::vector<vec3> &p_pov,
                         const std::vector<vec3> &n_pov, real l0,
                         real p = 3.0) {

  const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);

  std::vector<real> weights = asawa::shell::face_areas(M, x);

  std::vector<vec18> S(p_pov.size(), vec18::Zero());
  std::vector<real> W = std::vector<real>(p_pov.size(), 0.0);
  std::vector<vec3> Ns = asawa::shell::face_normals(M, x);
  for (int i = 0; i < Ns.size(); i++) {
    Ns[i] = weights[i] * Ns[i];
  }

  std::vector<real> us = integrate_over_shell<real>(
      M, p_pov,
      [&weights, &Ns](const std::vector<index_t> &edge_ids,
                      Shell_Sum_Type &sum) {
        sum.bind(calder::scalar_datum::create(edge_ids, weights));
        sum.bind(calder::vec3_datum::create(edge_ids, Ns));
      },
      [&](const index_t i, const index_t j, //
          const vec3 &pi, const vec3 &pj,
          const std::vector<calder::datum::ptr> &data,
          Shell_Sum_Type::Node_Type node_type, //
          const Shell_Sum_Type::Node &node,    //
          const Shell_Sum_Type::Tree &tree) -> real {
        real wj = get_data<real>(node_type, j, 0, data);
        vec3 Nj = get_data<vec3>(node_type, j, 1, data);
        vec3 Ni = n_pov[i];

        vec3 dp = pj - pi;
        if (dp.dot(Ni) > 0) {
          return 0.0;
        }

        real dist = dp.norm();
        // real kappa = computeK(dist, l0, p);
        real kappa = computeKg(dist, l0, p);
        //  real kappa = computeKg(dist, l0, 2.0);

        real w_total = wj * kappa;
        // real w_total = 1.0;

        // S[i] += mk_sphere_v(w_total * dp);
        S[i] += w_total * mk_sphere_v(dp);
        W[i] += w_total;
        return 0.0;
      });

  std::vector<vec4> out(p_pov.size(), vec4::Zero());
  for (int i = 0; i < S.size(); i++) {
    const vec18 &Si = S[i];
    out[i] = calc_sphere(W[i], Si);
    vec3 x = p_pov[i];
    vec3 dc = out[i].segment(0, 3);
    out[i].segment(0, 3) = x + dc;
    //  logger::line(x, x + dc, vec4(0.0, 1.0, 0.0, 1.0));
  }

  return out;
}

mat3 mk_cyl_A(vec2 dx, const vec2 &N) {
  mat3 A;
  real x = dx[0];
  real y = dx[1];
  real x2 = 2.0 * x;
  real y2 = 2.0 * y;
  real Nx = N[0];
  real Ny = N[1];

  real Nx1 = sqrt(1.0 - Nx * Nx);
  real Ny1 = sqrt(1.0 - Ny * Ny);

  A.row(0) = vec3(x2, y2, 1.0);
  // A = X^2 * I + sort of Skew(X) + 1
  A.row(1) = vec3(Nx1, -Nx, 0.0);
  A.row(2) = vec3(-Ny, Ny1, 0.0);

  return A;
}

vec3 mk_cyl_B(const vec2 dx, const vec2 &N) {
  real x = dx[0];
  real y = dx[1];
  real xx = x * x;
  real yy = y * y;
  real Nx = N[0];
  real Ny = N[1];

  real Nx1 = sqrt(1.0 - Nx * Nx);
  real Ny1 = sqrt(1.0 - Ny * Ny);
  real b0 = xx + yy;
  real b1 = Nx1 * x - Nx * y;
  real b2 = -Ny * x + Ny1 * y;
  return vec3(b0, b1, b2);
}

vec3 mk_cyl_A(vec2 dx) {
  mat3 A;
  real x = dx[0];
  real y = dx[1];
  real x2 = 2.0 * x;
  real y2 = 2.0 * y;
  return vec3(x2, y2, 1.0);
}

real mk_cyl_B(const vec2 dx) {
  real x = dx[0];
  real y = dx[1];
  real xx = x * x;
  real yy = y * y;
  return xx + yy;
}

std::vector<vec7> cylinder(asawa::shell::shell &M,
                           const std::vector<vec3> &p_pov,
                           const std::vector<vec3> &n_pov, real l0,
                           real p = 3.0) {

  const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);
  std::vector<vec3> Nv = asawa::shell::vertex_normals(M, x);

  std::vector<real> weights = asawa::shell::face_areas(M, x);
  // what would be cool would be to use the normal frame to compute
  // a an angle measure between the observer null space and the face null space
  // this would require smearing the normal frame over the face or computing

  std::vector<mat3> Fn = normal_covariant_frame(M, p_pov, n_pov, l0, p);
  std::vector<mat3> Fnv = normal_covariant_frame(M, x, Nv, l0, p);
  std::vector<vec3> Yv(x.size(), vec3::Zero());
  for (int i = 0; i < x.size(); i++) {
    Yv[i] = Fnv[i].col(2);
  }
  std::vector<vec3> Yf = asawa::shell::vert_to_face<vec3>(M, x, Yv);

  std::vector<mat3> As = std::vector<mat3>(p_pov.size(), mat3::Zero());
  std::vector<vec3> bs = std::vector<vec3>(p_pov.size(), vec3::Zero());

  std::vector<vec3> Ns = asawa::shell::face_normals(M, x);
  for (int i = 0; i < Ns.size(); i++) {
    Ns[i] = weights[i] * Ns[i];
    Yf[i] = weights[i] * Yf[i].normalized();
  }

  std::vector<real> us = integrate_over_shell<real>(
      M, p_pov,
      [&](const std::vector<index_t> &face_ids, Shell_Sum_Type &sum) {
        sum.bind(calder::scalar_datum::create(face_ids, weights));
        sum.bind(calder::vec3_datum::create(face_ids, Ns));
        sum.bind(calder::vec3_datum::create(face_ids, Yf));
      },
      [&](const index_t i, const index_t j, //
          const vec3 &pi, const vec3 &pj,
          const std::vector<calder::datum::ptr> &data,
          Shell_Sum_Type::Node_Type node_type, //
          const Shell_Sum_Type::Node &node,    //
          const Shell_Sum_Type::Tree &tree) -> real {
        real wj = get_data<real>(node_type, j, 0, data);
        vec3 Nj = get_data<vec3>(node_type, j, 1, data);
        vec3 Yj = get_data<vec3>(node_type, j, 2, data);
        Yj = Yj.normalized();

        Nj = Nj.normalized();

        vec3 Ni = n_pov[i];

        mat3 Fi = Fn[i];
        mat23 P;
        P.row(0) = Fi.col(0);
        P.row(1) = Fi.col(1);
        vec3 Yi = Fi.col(2);
        vec3 dp = pj - pi;
        real Ndp = dp.dot(Nj);

        vec2 dp2 = P * dp;
        vec2 Nj2 = P * Nj;
        vec2 Ni2 = P * Ni;

        Nj2 = Nj2.normalized();
        Ni2 = Ni2.normalized();
        if (dp2.dot(Ni2) > 0) {
          return 0.0;
        }
        mat3 Ai = mk_cyl_A(dp2, Nj2);
        vec3 bi = mk_cyl_B(dp2, Nj2);

        // vec3 Ai = mk_cyl_A(dp2);
        // real bi = mk_cyl_B(dp2);

        real dist = dp.norm();
        real kappa = computeK(dist, l0, p);
        real Yij = pow(Yi.dot(Yj), 2.0);
        real w = Yij * wj * kappa;
    // real kappa = computeKg(dist, l0, p);
    //   real kappa = computeKg(dist, l0, 2.0);
#if 0
        if (i == 0) {
          vec3 dpUp = P.transpose() * dp2;
          logger::line(pi, pi + dpUp, vec4(0.0, 1.0, 0.0, 1.0));
          logger::line(pi, pi + 0.01 * Fi.col(2), vec4(1.0, 0.0, 0.0, 1.0));
          logger::line(vec3::Zero(), vec3(dp2[0], dp2[1], 0.0),
                       vec4(0.0, 0.0, 1.0, 1.0));
        }
#endif
        // As[i] += wj * kappa * Ai * Ai.transpose();
        // bs[i] += wj * kappa * bi * Ai.transpose();
        As[i] += w * Ai.transpose() * Ai;
        bs[i] += w * Ai.transpose() * bi;

        return 0.0;
      });

  std::vector<vec7> out(p_pov.size(), vec7::Zero());
  for (int i = 0; i < As.size(); i++) {
    vec3 xi = p_pov[i];
    mat3 &A = As[i];
    vec3 &b = bs[i];
    mat3 Fi = Fn[i];
    vec3 Ni = Fi.col(2);
    mat23 P;
    P.row(0) = Fi.col(0);
    P.row(1) = Fi.col(1);

    vec3 z = A.colPivHouseholderQr().solve(b);
    vec2 z2 = z.segment(0, 2);
    vec3 cen = xi + P.transpose() * z2;

    real r = sqrt(z[0] * z[0] + z[1] * z[1] + z[2]);

    vec3 xp = va::project_on_line(cen, vec3(cen + Ni), xi);
    vec3 dpn = (xp - xi).normalized();
#if 0
    logger::line(xi - 0.01 * Ni, xi + 0.01 * Ni, vec4(0.0, 1.0, 1.0, 1.0));
    // logger::line(xi, xp, vec4(1.0, 0.0, 1.0, 1.0));
    logger::line(cen, cen - r * dpn, vec4(1.0, 0.0, 1.0, 1.0));
    // logger::line(vec3::Zero(), vec3(z2[0], z2[1], 0.0),
    //              vec4(1.0, 0.0, 1.0, 1.0));
#endif
    out[i].segment(0, 3) = cen;
    out[i].segment(3, 3) = Ni;
    out[i][6] = r;
  }

  return out;
}

} // namespace calder
} // namespace gaudi
#endif