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
#include "gaudi/asawa/datums.hpp"
#include "gaudi/asawa/shell/datum_x.hpp"
#include "rod_integrators.hpp"
#include "shell_integrators.hpp"
#include "weight_functions.hpp"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <ostream>
#include <vector>
#include "gaudi/geometry_logger.hpp"

namespace gaudi
{

  namespace calder
  {

    std::vector<vec3> tangent_point_gradient(asawa::rod::rod &R,
                                             const std::vector<vec3> &p_pov,
                                             const std::vector<real> &w_pov,
                                             const std::vector<vec3> &T_pov,
                                             real l0, real p = 3.0)
    {
      std::vector<real> weights = R.l0();
      std::vector<vec3> Tc = R.N2c();

      auto bind = [&weights, &Tc](const std::vector<index_t> &edge_ids,
                                  cov_rod_sum &sum) {
        sum.bind(calder::scalar_datum::create(edge_ids, weights));
        sum.bind(calder::vec3_datum::create(edge_ids, Tc));
      };
      auto kernel = [l0, &w_pov, &T_pov, p](
                        const index_t i, const index_t j, const vec3 &pi,
                        const vec3 &pj, const std::vector<calder::datum::ptr> &data,
                        Node_Type node_type,
                        const arp::edge_cov_tree<2> &) -> vec3 {
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
            vec3 gj = calc_tangent_point_radius_grad(dp, Nj, l0, p);
            vec3 gi = calc_tangent_point_radius_grad(-dp, Ni, l0, p);
            //vec3 g = 0.5 * (wi * gi - wj * gj);
            //std::cout << "g: " << g.transpose() << std::endl;
            return 0.5 * (wi * gi - wj * gj);
          };
      return integrate_over_rod_cov<vec3>(R, p_pov, bind, kernel, kernel);
    }

    /// Cauchy leaf gradient: exact ∇ of K = |N·dp|^p / (|dp|²+ε²)^p, ε=l0.
    std::vector<vec3> tangent_point_gradient_cauchy(
        asawa::rod::rod &R, const std::vector<vec3> &p_pov,
        const std::vector<real> &w_pov, const std::vector<vec3> &T_pov,
        real l0, real p = 3.0, real bh_eps = 0.25,
        const arp::edge_cov_tree<2> *cov = nullptr, real cone_eps = -1) {
      std::vector<real> weights = R.l0();
      std::vector<vec3> Tc = R.N2c();

      auto bind = [&weights, &Tc](const std::vector<index_t> &edge_ids,
                                  cov_rod_sum &sum) {
        sum.bind(calder::scalar_datum::create(edge_ids, weights));
        sum.bind(calder::vec3_datum::create(edge_ids, Tc));
      };
      auto kernel = [l0, &w_pov, &T_pov, p](
                        const index_t i, const index_t j, const vec3 &pi,
                        const vec3 &pj, const std::vector<calder::datum::ptr> &data,
                        Node_Type node_type,
                        const arp::edge_cov_tree<2> &) -> vec3 {
            real wi = w_pov[i];
            vec3 Ti = T_pov[i];
            real wj = get_data<real>(node_type, j, 0, data);
            vec3 Tj = get_data<vec3>(node_type, j, 1, data);
            Ti.normalize();
            Tj.normalize();
            vec3 dp = pj - pi;
            vec3 ci = Ti.cross(dp);
            vec3 cj = Tj.cross(dp);
            if (ci.squaredNorm() < real(1.0e-32) ||
                cj.squaredNorm() < real(1.0e-32))
              return vec3::Zero();
            vec3 Ni = ci.normalized().cross(Ti).normalized();
            vec3 Nj = cj.normalized().cross(Tj).normalized();
            const vec3 gj =
                calc_tangent_point_radius_grad_cauchy(dp, Nj, l0, p);
            const vec3 gi =
                calc_tangent_point_radius_grad_cauchy(-dp, Ni, l0, p);
            return real(0.5) * (wi * gi - wj * gj);
          };
      if (cov)
        return integrate_over_cov_tree<vec3>(*cov, rod_cov_edge_ids(R), p_pov,
                                             bind, kernel, kernel, bh_eps, 1,
                                             cone_eps);
      return integrate_over_rod_cov<vec3>(R, p_pov, bind, kernel, kernel,
                                          bh_eps, 1, cone_eps);
    }

    /// Matrix-free Hess·v for Cauchy TP (exact leaf Hess, same pair assembly).
    std::vector<vec3> tangent_point_hessian_vec_cauchy(
        asawa::rod::rod &R, const std::vector<vec3> &p_pov,
        const std::vector<real> &w_pov, const std::vector<vec3> &T_pov,
        const std::vector<vec3> &v_pov, real l0, real p = 3.0) {
      std::vector<real> weights = R.l0();
      std::vector<vec3> Tc = R.N2c();

      return integrate_over_rod<vec3>(
          R, p_pov,
          [&weights, &Tc, &v_pov](const std::vector<index_t> &edge_ids,
                                  Rod_Sum_Type &sum) {
            sum.bind(calder::scalar_datum::create(edge_ids, weights));
            sum.bind(calder::vec3_datum::create(edge_ids, Tc));
            sum.bind(calder::vec3_datum::create(edge_ids, v_pov));
          },
          [l0, &w_pov, &T_pov, &v_pov, p](
              const index_t i, const index_t j, const vec3 &pi, const vec3 &pj,
              const std::vector<calder::datum::ptr> &data,
              Rod_Sum_Type::Node_Type node_type,
              const Rod_Sum_Type::Tree & /*tree*/) -> vec3 {
            real wi = w_pov[i];
            vec3 Ti = T_pov[i];
            vec3 vi = v_pov[i];
            real wj = get_data<real>(node_type, j, 0, data);
            vec3 Tj = get_data<vec3>(node_type, j, 1, data);
            vec3 vj = get_data<vec3>(node_type, j, 2, data);
            Ti.normalize();
            Tj.normalize();
            vec3 dp = pj - pi;
            vec3 ci = Ti.cross(dp);
            vec3 cj = Tj.cross(dp);
            if (ci.squaredNorm() < real(1.0e-32) ||
                cj.squaredNorm() < real(1.0e-32))
              return vec3::Zero();
            vec3 Ni = ci.normalized().cross(Ti).normalized();
            vec3 Nj = cj.normalized().cross(Tj).normalized();
            const vec3 dv = vj - vi;
            const mat3 Hj =
                calc_tangent_point_radius_hessian_cauchy(dp, Nj, l0, p);
            const mat3 Hi =
                calc_tangent_point_radius_hessian_cauchy(-dp, Ni, l0, p);
            const vec3 gj = Hj * dv;
            const vec3 gi = Hi * (-dv);
            vec3 out = real(0.5) * (wi * gi - wj * gj);
            if (!std::isfinite(out[0]) || !std::isfinite(out[1]) ||
                !std::isfinite(out[2]))
              return vec3::Zero();
            return out;
          });
    }

    /// Matrix-free Hess·v for regularized TP: Jac(Reg(∇E0)) · v.
    std::vector<vec3> tangent_point_hessian_vec(
        asawa::rod::rod &R, const std::vector<vec3> &p_pov,
        const std::vector<real> &w_pov, const std::vector<vec3> &T_pov,
        const std::vector<vec3> &v_pov, real l0, real p = 3.0) {
      std::vector<real> weights = R.l0();
      std::vector<vec3> Tc = R.N2c();

      std::vector<vec3> us = integrate_over_rod<vec3>(
          R, p_pov,
          [&weights, &Tc, &v_pov](const std::vector<index_t> &edge_ids,
                                  Rod_Sum_Type &sum) {
            sum.bind(calder::scalar_datum::create(edge_ids, weights));
            sum.bind(calder::vec3_datum::create(edge_ids, Tc));
            sum.bind(calder::vec3_datum::create(edge_ids, v_pov));
          },
          [l0, &w_pov, &T_pov, &v_pov, p](
              const index_t i, const index_t j, const vec3 &pi, const vec3 &pj,
              const std::vector<calder::datum::ptr> &data,
              Rod_Sum_Type::Node_Type node_type,
              const Rod_Sum_Type::Tree & /*tree*/) -> vec3 {
            real wi = w_pov[i];
            vec3 Ti = T_pov[i];
            vec3 vi = v_pov[i];
            real wj = get_data<real>(node_type, j, 0, data);
            vec3 Tj = get_data<vec3>(node_type, j, 1, data);
            vec3 vj = get_data<vec3>(node_type, j, 2, data);
            Ti.normalize();
            Tj.normalize();

            vec3 dp = pj - pi;
            vec3 ci = Ti.cross(dp);
            vec3 cj = Tj.cross(dp);
            if (ci.squaredNorm() < real(1.0e-32) ||
                cj.squaredNorm() < real(1.0e-32))
              return vec3::Zero();
            vec3 Bi = ci.normalized();
            vec3 Bj = cj.normalized();
            vec3 Ni = Bi.cross(Ti).normalized();
            vec3 Nj = Bj.cross(Tj).normalized();
            const vec3 dv = vj - vi;
            const mat3 Hj =
                calc_tangent_point_radius_hessian_reg(dp, Nj, l0, p);
            const mat3 Hi =
                calc_tangent_point_radius_hessian_reg(-dp, Ni, l0, p);
            const vec3 gj = Hj * dv;
            const vec3 gi = Hi * (-dv);

            vec3 out = real(0.5) * (wi * gi - wj * gj);
            if (!std::isfinite(out[0]) || !std::isfinite(out[1]) ||
                !std::isfinite(out[2]))
              return vec3::Zero();
            return out;
          });
      return us;
    }

    /// Matrix-free Hess·v for soft TP density (same pair assembly as gradient).
    /// For each pair: dv = vj-vi, contrib mirrors 0.5*(wi Hi(-dv) - wj Hj dv).
    std::vector<vec3> tangent_point_hessian_vec_soft(
        asawa::rod::rod &R, const std::vector<vec3> &p_pov,
        const std::vector<real> &w_pov, const std::vector<vec3> &T_pov,
        const std::vector<vec3> &v_pov, real R_min, real tau, real p = 3.0) {
      std::vector<real> weights = R.l0();
      std::vector<vec3> Tc = R.N2c();

      std::vector<vec3> us = integrate_over_rod<vec3>(
          R, p_pov,
          [&weights, &Tc, &v_pov](const std::vector<index_t> &edge_ids,
                                  Rod_Sum_Type &sum) {
            sum.bind(calder::scalar_datum::create(edge_ids, weights));
            sum.bind(calder::vec3_datum::create(edge_ids, Tc));
            sum.bind(calder::vec3_datum::create(edge_ids, v_pov));
          },
          [R_min, tau, &w_pov, &T_pov, &v_pov, p](
              const index_t i, const index_t j, const vec3 &pi, const vec3 &pj,
              const std::vector<calder::datum::ptr> &data,
              Rod_Sum_Type::Node_Type node_type,
              const Rod_Sum_Type::Tree & /*tree*/) -> vec3 {
            real wi = w_pov[i];
            vec3 Ti = T_pov[i];
            vec3 vi = v_pov[i];
            real wj = get_data<real>(node_type, j, 0, data);
            vec3 Tj = get_data<vec3>(node_type, j, 1, data);
            vec3 vj = get_data<vec3>(node_type, j, 2, data);
            Ti.normalize();
            Tj.normalize();

            vec3 dp = pj - pi;
            vec3 ci = Ti.cross(dp);
            vec3 cj = Tj.cross(dp);
            if (ci.squaredNorm() < real(1.0e-32) ||
                cj.squaredNorm() < real(1.0e-32))
              return vec3::Zero();
            vec3 Bi = ci.normalized();
            vec3 Bj = cj.normalized();
            vec3 Ni = Bi.cross(Ti).normalized();
            vec3 Nj = Bj.cross(Tj).normalized();
            const vec3 dv = vj - vi;
            const mat3 Hj =
                calc_tangent_point_radius_hessian_soft(dp, Nj, R_min, tau, p);
            const mat3 Hi =
                calc_tangent_point_radius_hessian_soft(-dp, Ni, R_min, tau, p);
            const vec3 gj = Hj * dv;
            const vec3 gi = Hi * (-dv);

            vec3 out = real(0.5) * (wi * gi - wj * gj);
            if (!std::isfinite(out[0]) || !std::isfinite(out[1]) ||
                !std::isfinite(out[2]))
              return vec3::Zero();
            return out;
          });
      return us;
    }

    /// Softmax-floor TP force (no Cauchy ℓ₀). See weight_functions.hpp.
    /// Optional probe: g_eff = g + α H dp, dp = pj − pi (leaf or branch center).
    std::vector<vec3> tangent_point_gradient_soft(
        asawa::rod::rod &R, const std::vector<vec3> &p_pov,
        const std::vector<real> &w_pov, const std::vector<vec3> &T_pov,
        real R_min, real tau, real p = 3.0, real hess_alpha = 0.0) {
      std::vector<real> weights = R.l0();
      std::vector<vec3> Tc = R.N2c();
      const bool use_hess = hess_alpha != real(0.0);

      auto bind = [&weights, &Tc](const std::vector<index_t> &edge_ids,
                                  cov_rod_sum &sum) {
        sum.bind(calder::scalar_datum::create(edge_ids, weights));
        sum.bind(calder::vec3_datum::create(edge_ids, Tc));
      };
      auto kernel = [R_min, tau, &w_pov, &T_pov, p, hess_alpha, use_hess](
                        const index_t i, const index_t j, const vec3 &pi,
                        const vec3 &pj, const std::vector<calder::datum::ptr> &data,
                        Node_Type node_type,
                        const arp::edge_cov_tree<2> &) -> vec3 {
            real wi = w_pov[i];
            vec3 Ti = T_pov[i];
            real wj = get_data<real>(node_type, j, 0, data);
            vec3 Tj = get_data<vec3>(node_type, j, 1, data);
            Ti.normalize();
            Tj.normalize();

            vec3 dp = pj - pi;
            vec3 ci = Ti.cross(dp);
            vec3 cj = Tj.cross(dp);
            if (ci.squaredNorm() < real(1.0e-32) ||
                cj.squaredNorm() < real(1.0e-32))
              return vec3::Zero();
            vec3 Bi = ci.normalized();
            vec3 Bj = cj.normalized();
            vec3 Ni = Bi.cross(Ti).normalized();
            vec3 Nj = Bj.cross(Tj).normalized();
            vec3 gj =
                calc_tangent_point_radius_gradient_soft(dp, Nj, R_min, tau, p);
            vec3 gi =
                calc_tangent_point_radius_gradient_soft(-dp, Ni, R_min, tau, p);
            if (use_hess) {
              const mat3 Hj = calc_tangent_point_radius_hessian_soft(
                  dp, Nj, R_min, tau, p);
              const mat3 Hi = calc_tangent_point_radius_hessian_soft(
                  -dp, Ni, R_min, tau, p);
              gj -= hess_alpha * (Hj * (dp));
              gi += hess_alpha * (Hi * (dp));
            }
            vec3 out = 0.5 * (wi * gi - wj * gj);
            if (!std::isfinite(out[0]) || !std::isfinite(out[1]) ||
                !std::isfinite(out[2]))
              return vec3::Zero();
            return out;
          };
      return integrate_over_rod_cov<vec3>(R, p_pov, bind, kernel, kernel);
    }

    /// Regularized TP energy matching regularized gradient:
    ///   E_ij = ½ (wi f(-dp,Ni) + wj f(dp,Nj)),  ∇_i ⇒ ½(wi gi − wj gj).
    std::vector<real> tangent_point_energy(asawa::rod::rod &R,
                                           const std::vector<vec3> &p_pov,
                                           const std::vector<real> &w_pov,
                                           const std::vector<vec3> &T_pov,
                                           real l0, real p = 3.0)
    {
      std::vector<real> weights = R.l0();
      std::vector<vec3> Tc = R.N2c();

      auto bind = [&weights, &Tc](const std::vector<index_t> &edge_ids,
                                  cov_rod_sum &sum) {
        sum.bind(calder::scalar_datum::create(edge_ids, weights));
        sum.bind(calder::vec3_datum::create(edge_ids, Tc));
      };
      auto kernel = [l0, &w_pov, &T_pov, p](
                        const index_t i, const index_t j, const vec3 &pi,
                        const vec3 &pj, const std::vector<calder::datum::ptr> &data,
                        Node_Type node_type,
                        const arp::edge_cov_tree<2> &) -> real {
            real wi = w_pov[i];
            vec3 Ti = T_pov[i];
            real wj = get_data<real>(node_type, j, 0, data);
            vec3 Tj = get_data<vec3>(node_type, j, 1, data);
            wi = std::max(wi, real(1e-6));
            wj = std::max(wj, real(1e-6));
            Ti.normalize();
            Tj.normalize();

            vec3 dp = pj - pi;
            vec3 ci = Ti.cross(dp);
            vec3 cj = Tj.cross(dp);
            if (ci.squaredNorm() < real(1e-30) || cj.squaredNorm() < real(1e-30))
              return real(0.0);
            vec3 Ni = ci.normalized().cross(Ti).normalized();
            vec3 Nj = cj.normalized().cross(Tj).normalized();
            const real fi =
                calc_tangent_point_inverse_radius(-dp, Ni, l0, p);
            const real fj =
                calc_tangent_point_inverse_radius(dp, Nj, l0, p);
            return real(0.5) * (wi * fi + wj * fj);
          };
      return integrate_over_rod_cov<real>(R, p_pov, bind, kernel, kernel);
    }

    /// Cauchy TP energy: K = |N·dp|^p / (|dp|²+ε²)^p, ε=l0.
    /// Pair form matches gradient: E_ij = ½(wi fi + wj fj).
    std::vector<real> tangent_point_energy_cauchy(
        asawa::rod::rod &R, const std::vector<vec3> &p_pov,
        const std::vector<real> &w_pov, const std::vector<vec3> &T_pov,
        real l0, real p = 3.0, real bh_eps = 0.25, real cone_eps = -1) {
      std::vector<real> weights = R.l0();
      std::vector<vec3> Tc = R.N2c();

      auto bind = [&weights, &Tc](const std::vector<index_t> &edge_ids,
                                  cov_rod_sum &sum) {
        sum.bind(calder::scalar_datum::create(edge_ids, weights));
        sum.bind(calder::vec3_datum::create(edge_ids, Tc));
      };
      auto kernel = [l0, &w_pov, &T_pov, p](
                        const index_t i, const index_t j, const vec3 &pi,
                        const vec3 &pj, const std::vector<calder::datum::ptr> &data,
                        Node_Type node_type,
                        const arp::edge_cov_tree<2> &) -> real {
            real wi = w_pov[i];
            vec3 Ti = T_pov[i];
            real wj = get_data<real>(node_type, j, 0, data);
            vec3 Tj = get_data<vec3>(node_type, j, 1, data);
            wi = std::max(wi, real(1e-6));
            wj = std::max(wj, real(1e-6));
            Ti.normalize();
            Tj.normalize();
            vec3 dp = pj - pi;
            vec3 ci = Ti.cross(dp);
            vec3 cj = Tj.cross(dp);
            if (ci.squaredNorm() < real(1e-30) || cj.squaredNorm() < real(1e-30))
              return real(0.0);
            vec3 Ni = ci.normalized().cross(Ti).normalized();
            vec3 Nj = cj.normalized().cross(Tj).normalized();
            const real fi =
                calc_tangent_point_inverse_radius_cauchy(-dp, Ni, l0, p);
            const real fj =
                calc_tangent_point_inverse_radius_cauchy(dp, Nj, l0, p);
            return real(0.5) * (wi * fi + wj * fj);
          };
      return integrate_over_rod_cov<real>(R, p_pov, bind, kernel, kernel, bh_eps,
                                          1, cone_eps);
    }

    /// Softmax-floor TP energy matching soft gradient:
    ///   E_ij = ½ (wi f(-dp,Ni) + wj f(dp,Nj)),  ∇_i ⇒ ½(wi gi − wj gj).
    std::vector<real> tangent_point_energy_soft(
        asawa::rod::rod &R, const std::vector<vec3> &p_pov,
        const std::vector<real> &w_pov, const std::vector<vec3> &T_pov,
        real R_min, real tau, real p = 3.0) {
      std::vector<real> weights = R.l0();
      std::vector<vec3> Tc = R.N2c();

      auto bind = [&weights, &Tc](const std::vector<index_t> &edge_ids,
                                  cov_rod_sum &sum) {
        sum.bind(calder::scalar_datum::create(edge_ids, weights));
        sum.bind(calder::vec3_datum::create(edge_ids, Tc));
      };
      auto kernel = [R_min, tau, &w_pov, &T_pov, p](
                        const index_t i, const index_t j, const vec3 &pi,
                        const vec3 &pj, const std::vector<calder::datum::ptr> &data,
                        Node_Type node_type,
                        const arp::edge_cov_tree<2> &) -> real {
            real wi = w_pov[i];
            vec3 Ti = T_pov[i];
            real wj = get_data<real>(node_type, j, 0, data);
            vec3 Tj = get_data<vec3>(node_type, j, 1, data);
            wi = std::max(wi, real(1e-6));
            wj = std::max(wj, real(1e-6));
            Ti.normalize();
            Tj.normalize();
            vec3 dp = pj - pi;
            vec3 ci = Ti.cross(dp);
            vec3 cj = Tj.cross(dp);
            if (ci.squaredNorm() < real(1e-30) || cj.squaredNorm() < real(1e-30))
              return real(0.0);
            vec3 Ni = ci.normalized().cross(Ti).normalized();
            vec3 Nj = cj.normalized().cross(Tj).normalized();
            const real fi =
                calc_tangent_point_radius_soft(-dp, Ni, R_min, tau, p);
            const real fj =
                calc_tangent_point_radius_soft(dp, Nj, R_min, tau, p);
            return real(0.5) * (wi * fi + wj * fj);
          };
      return integrate_over_rod_cov<real>(R, p_pov, bind, kernel, kernel);
    }

    std::vector<real> tangent_point_energy(asawa::shell::shell &M,
                                           const std::vector<vec3> &p_pov,
                                           const std::vector<real> &w_pov,
                                           const std::vector<vec3> &N_pov, //
                                           real l0, real p = 3.0)
    {
      std::vector<vec3> x = asawa::get_vec_data(M, 0);
      std::vector<real> weights = asawa::shell::face_areas(M, x);
      std::vector<vec3> Nc = asawa::shell::face_normals(M, x);

      std::vector<real> us = integrate_over_shell<real>(
          M, p_pov,
          [&weights, &Nc](const std::vector<index_t> &face_ids,
                          Shell_Sum_Type &sum)
          {
            sum.bind(calder::scalar_datum::create(face_ids, weights));
            sum.bind(calder::vec3_datum::create(face_ids, Nc));
          },
          [l0, &w_pov, &N_pov, p](const index_t i, const index_t j, //
                                  const vec3 &pi, const vec3 &pj,
                                  const std::vector<calder::datum::ptr> &data,
                                  Shell_Sum_Type::Node_Type node_type, //
                                  const Shell_Sum_Type::Tree &tree) -> real
          {
            real wi = w_pov[i];
            vec3 Ni = N_pov[i];

            real wj = get_data<real>(node_type, j, 0, data);
            vec3 Nj = get_data<vec3>(node_type, j, 1, data);

            wi = std::max(wi, 1e-6);
            wj = std::max(wj, 1e-6);

            Ni.normalize();
            Nj.normalize();

            vec3 dp = pj - pi;

            const real fj = calc_tangent_point_inverse_radius(dp, Nj, l0, p);
            const real fi = calc_tangent_point_inverse_radius(-dp, Ni, l0, p);
            // Match shell gradient assembly: ½(wi fi + wj fj).
            return real(0.5) * (wi * fi + wj * fj);
          });
      return us;
    }

    std::vector<vec3> tangent_point_gradient(asawa::shell::shell &M,
                                             const std::vector<vec3> &p_pov,
                                             const std::vector<real> &w_pov,
                                             const std::vector<vec3> &N_pov,
                                             real l0, real p = 3.0)
    {
      std::vector<vec3> x = asawa::get_vec_data(M, 0);
      std::vector<real> weights = asawa::shell::face_areas(M, x);
      std::vector<vec3> Nc = asawa::shell::face_normals(M, x);

      std::vector<vec3> us = integrate_over_shell<vec3>(
          M, p_pov,
          [&weights, &Nc](const std::vector<index_t> &face_ids,
                          Shell_Sum_Type &sum)
          {
            sum.bind(calder::scalar_datum::create(face_ids, weights));
            sum.bind(calder::vec3_datum::create(face_ids, Nc));
          },
          [l0, &w_pov, &N_pov, p](const index_t i, const index_t j, //
                                  const vec3 &pi, const vec3 &pj,
                                  const std::vector<calder::datum::ptr> &data,
                                  Shell_Sum_Type::Node_Type node_type, //
                                  const Shell_Sum_Type::Tree &tree) -> vec3
          {
            real wi = w_pov[i];
            vec3 Ni = N_pov[i];
            real wj = get_data<real>(node_type, j, 0, data);
            vec3 Nj = get_data<vec3>(node_type, j, 1, data);
            wi = std::max(wi, 1e-6);
            wj = std::max(wj, 1e-6);

            vec3 dp = pj - pi;
            Ni.normalize();
            Nj.normalize();
            real pk = p;
            vec3 gj = calc_tangent_point_radius_grad(dp, Nj, l0, pk);
            vec3 gi = calc_tangent_point_radius_grad(-dp, Ni, l0, pk);

            return wi * gi - wj * gj;
          });
      return us;
    }

    /// Axis-restricted (soft-slice) shell TP gradient — Cauchy core, same
    /// energy family as braid_planar_demo (`tangent_point_type::cauchy`):
    ///   K = |N·dp|^p / (|dp|² + ε²)^p,  ε = l0.
    /// Evaluated in the null plane Πᵢ through pᵢ (project dp and normals),
    /// with Gaussian window κ = calc_gaussian(axial·e₁, σ).
    /// Optionally project the force back into Πᵢ.
    std::vector<vec3> tangent_point_gradient_axis_slice(
        asawa::shell::shell &M, const std::vector<vec3> &p_pov,
        const std::vector<real> &w_pov, const std::vector<vec3> &N_pov,
        const std::vector<vec3> &A_pov, const std::vector<real> &sigma_par,
        real l0, real p = 6.0, bool project_force_to_plane = true) {
      assert(w_pov.size() == p_pov.size());
      assert(N_pov.size() == p_pov.size());
      assert(A_pov.size() == p_pov.size());
      assert(sigma_par.size() == p_pov.size());

      std::vector<vec3> x = asawa::get_vec_data(M, 0);
      std::vector<real> weights = asawa::shell::face_areas(M, x);
      std::vector<vec3> Nc = asawa::shell::face_normals(M, x);

      return integrate_over_shell<vec3>(
          M, p_pov,
          [&weights, &Nc](const std::vector<index_t> &face_ids,
                          Shell_Sum_Type &sum) {
            sum.bind(calder::scalar_datum::create(face_ids, weights));
            sum.bind(calder::vec3_datum::create(face_ids, Nc));
          },
          [l0, p, project_force_to_plane, &w_pov, &N_pov, &A_pov, &sigma_par](
              const index_t i, const index_t j, const vec3 &pi, const vec3 &pj,
              const std::vector<calder::datum::ptr> &data,
              Shell_Sum_Type::Node_Type node_type,
              const Shell_Sum_Type::Tree & /*tree*/) -> vec3 {
            const size_t ii = static_cast<size_t>(i);
            vec3 ai = A_pov[ii];
            if (!ai.allFinite() || ai.squaredNorm() < real(1e-24))
              return vec3::Zero();
            ai.normalize();

            real wi = std::max(w_pov[ii], real(1e-6));
            vec3 Ni = N_pov[ii];
            real wj = get_data<real>(node_type, j, 0, data);
            vec3 Nj = get_data<vec3>(node_type, j, 1, data);
            wj = std::max(wj, real(1e-6));
            if (!Ni.allFinite() || !Nj.allFinite() || !std::isfinite(wj))
              return vec3::Zero();
            Ni.normalize();
            Nj.normalize();

            const vec3 dp = pj - pi;
            if (!dp.allFinite())
              return vec3::Zero();
            const real axial = dp.dot(ai);
            const real sigma = std::max(sigma_par[ii], real(1e-12));
            const real kappa =
                calc_gaussian(vec3(axial, real(0.0), real(0.0)), sigma);
            if (!(kappa > real(1e-8)))
              return vec3::Zero();

            const vec3 dp_plane = dp - axial * ai;
            if (dp_plane.squaredNorm() < real(1e-32))
              return vec3::Zero();

            vec3 Ni_p = Ni - Ni.dot(ai) * ai;
            vec3 Nj_p = Nj - Nj.dot(ai) * ai;
            if (Ni_p.squaredNorm() < real(1e-24) ||
                Nj_p.squaredNorm() < real(1e-24))
              return vec3::Zero();
            Ni_p.normalize();
            Nj_p.normalize();

            // Cauchy core (braid / rod TP): ∇ K with ε = l0.
            const vec3 gj =
                calc_tangent_point_radius_grad_cauchy(dp_plane, Nj_p, l0, p);
            const vec3 gi =
                calc_tangent_point_radius_grad_cauchy(-dp_plane, Ni_p, l0, p);
            if (!gi.allFinite() || !gj.allFinite())
              return vec3::Zero();
            // Match rod cauchy assembly: ½(wi gi − wj gj).
            vec3 g = kappa * real(0.5) * (wi * gi - wj * gj);
            if (!g.allFinite())
              return vec3::Zero();
            if (project_force_to_plane)
              g = g - g.dot(ai) * ai;
            return g;
          });
    }

    std::vector<mat3> tangent_point_gradient_frame(asawa::shell::shell &M,
                                                   const std::vector<vec3> &p_pov,
                                                   const std::vector<real> &w_pov,
                                                   const std::vector<vec3> &N_pov,
                                                   real l0, real p = 3.0)
    {
      std::vector<vec3> x = asawa::get_vec_data(M, 0);
      std::vector<real> weights = asawa::shell::face_areas(M, x);
      std::vector<vec3> Nc = asawa::shell::face_normals(M, x);

      std::vector<real> sums(x.size(), 0.0);
      std::vector<mat3> us = integrate_over_shell<mat3>(
          M, p_pov,
          [&weights, &Nc](const std::vector<index_t> &face_ids,
                          Shell_Sum_Type &sum)
          {
            sum.bind(calder::scalar_datum::create(face_ids, weights));
            sum.bind(calder::vec3_datum::create(face_ids, Nc));
          },
          [l0, &w_pov, &N_pov, p](const index_t i, const index_t j, //
                                  const vec3 &pi, const vec3 &pj,
                                  const std::vector<calder::datum::ptr> &data,
                                  Shell_Sum_Type::Node_Type node_type, //
                                  const Shell_Sum_Type::Tree &tree) -> mat3
          {
            real wi = w_pov[i];
            vec3 Ni = N_pov[i];
            real wj = get_data<real>(node_type, j, 0, data);
            vec3 Nj = get_data<vec3>(node_type, j, 1, data);
            wi = std::max(wi, 1e-6);
            wj = std::max(wj, 1e-6);

            vec3 dp = pj - pi;
            // if (i == 0) {
            //   geometry_logger::line(pi, pj, vec4(0.0, 1.0, 0.5, 1.0));
            //   geometry_logger::line(pj, pj + 0.1 * Nj, vec4(0.0, 1.0, 0.5, 1.0));
            // }
            Ni.normalize();
            Nj.normalize();
            real pk = p;
            vec3 gj = calc_tangent_point_radius_grad(dp, Nj, l0, pk);
            vec3 gi = calc_tangent_point_radius_grad(-dp, Ni, l0, pk);
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

    std::vector<vec3> tangent_point_center(asawa::shell::shell &M,
                                           const std::vector<vec3> &p_pov,
                                           const std::vector<real> &w_pov,
                                           const std::vector<vec3> &N_pov, real l0,
                                           real p = 3.0)
    {
      std::vector<vec3> x = asawa::get_vec_data(M, 0);
      std::vector<real> weights = asawa::shell::face_areas(M, x);
      std::vector<vec3> Nc = asawa::shell::face_normals(M, x);
      std::vector<vec3> acc(p_pov.size(), vec3::Zero());
      std::vector<real> acc_w(p_pov.size(), 0.0);

      std::vector<vec3> us = integrate_over_shell<vec3>(
          M, p_pov,
          [&weights, &Nc](const std::vector<index_t> &face_ids,
                          Shell_Sum_Type &sum)
          {
            sum.bind(calder::scalar_datum::create(face_ids, weights));
            sum.bind(calder::vec3_datum::create(face_ids, Nc));
          },
          [l0, &w_pov, &N_pov, &acc, &acc_w,
           p](const index_t i, const index_t j, //
              const vec3 &pi, const vec3 &pj,
              const std::vector<calder::datum::ptr> &data,
              Shell_Sum_Type::Node_Type node_type, //
              const Shell_Sum_Type::Tree &tree) -> vec3
          {
            real wi = w_pov[i];
            vec3 Ni = N_pov[i];
            real wj = get_data<real>(node_type, j, 0, data);
            vec3 Nj = get_data<vec3>(node_type, j, 1, data);
            wi = std::max(wi, 1e-6);
            wj = std::max(wj, 1e-6);

            vec3 dp = pj - pi;
            real dpNj = dp.dot(Nj);
            if (dp.dot(Ni) > 0)
            {
              return vec3::Zero();
            }

            //  if (i == 0) {
            //    geometry_logger::line(pi, pj, vec4(0.0, 1.0, 0.5, 1.0));
            //    geometry_logger::line(pj, pj + 0.1 * Nj, vec4(0.0, 1.0, 0.5, 1.0));
            //  }
            Ni.normalize();
            Nj.normalize();
            real pk = p;
            real ri = calc_tangent_point_radius(dp, Ni);
            real rj = calc_tangent_point_radius(-dp, Nj);
            vec3 ci = pi - 0.5 * ri * Ni;
            vec3 cj = pj - 0.5 * rj * Nj;
            real kij = calc_inv_dist(dp, l0, p);
            acc[i] += kij * pj;
            acc_w[i] += kij;
            return vec3::Zero();
            // return wi * gi;

            // return -wj * gj;
          });
      for (int i = 0; i < acc.size(); i++)
      {
        us[i] = acc[i] / acc_w[i];
      }
      return us;
    }

    /// Unit-tangent dot. Zero when either tangent is degenerate.
    inline real tp_tangent_align(const vec3 &Ti, const vec3 &Tj) {
      const real ni = Ti.norm();
      const real nj = Tj.norm();
      if (!(ni > real(1e-12)) || !(nj > real(1e-12)))
        return real(0);
      return Ti.dot(Tj) / (ni * nj);
    }

    /// Sample of w Tp with w = κ ⟨T_i, T_j⟩.
    /// κ is the Gaussian mollifier of length l1 and power p1.
    /// ⟨T_i, T_j⟩ is the Repulsive Curves cross-term coefficient, held fixed
    /// (Yu, Schumacher, Crane, reference/crane2021-repulsive-curves.pdf, §5.2.2).
    inline real harmonic_tp_density(const vec3 &dp, const vec3 &Ti, const vec3 &Tj,
                                    real Tp, real l1, real p1 = 2.0) {
      if (!(l1 > real(1e-16)))
        return real(0);
      return calc_gaussian(dp, l1, p1) * tp_tangent_align(Ti, Tj) * Tp;
    }

    /// ∇(w Tp) = w′ Tp + w ∇Tp, w = κ ⟨T_i, T_j⟩.
    /// κ(l1, p1) is differentiated; ⟨T_i, T_j⟩ is the same frozen coefficient.
    /// ∇_i κ = −∂κ/∂dp. Matches the shell product rule
    /// (scripts/sympy/tp_harmonic_compose.py) with the extra tangent weight.
    inline vec3 harmonic_tp_pair(const vec3 &dp, const vec3 &Ti, const vec3 &Tj,
                                 real Tp, const vec3 &grad_Tp, real l1,
                                 real p1 = 2.0) {
      if (!(l1 > real(1e-16)))
        return vec3::Zero();
      const real align = tp_tangent_align(Ti, Tj);
      const real kappa = calc_gaussian(dp, l1, p1);
      const vec3 dkappa = calc_d_gaussian(dp, l1, p1);
      return (-dkappa) * (align * Tp) + (kappa * align) * grad_Tp;
    }

    /// Rod harmonic potential Φ = (κ ⟨T_i, T_j⟩) ∗ Tp.
    /// Tp is the Cauchy density (l0, p0). κ is the Gaussian (l1, p1).
    inline std::vector<real> tangent_point_harmonic_potential(
        asawa::rod::rod &R, const std::vector<vec3> &p_pov,
        const std::vector<real> &w_pov, const std::vector<vec3> &T_pov, real l0,
        real p0, real l1, real p1) {
      const std::vector<real> Ev =
          tangent_point_energy_cauchy(R, p_pov, w_pov, T_pov, l0, p0);
      std::vector<vec3> Tc = R.N2c();
      return integrate_over_rod<real>(
          R, p_pov,
          [&Ev, &Tc](const std::vector<index_t> &edge_ids, Rod_Sum_Type &sum) {
            sum.bind(calder::scalar_datum::create(edge_ids, Ev));
            sum.bind(calder::vec3_datum::create(edge_ids, Tc));
          },
          [&T_pov, l1, p1](const index_t i, const index_t j, const vec3 &pi,
                           const vec3 &pj,
                           const std::vector<calder::datum::ptr> &data,
                           Rod_Sum_Type::Node_Type node_type,
                           const Rod_Sum_Type::Tree &) -> real {
            const real Tp = get_data<real>(node_type, j, 0, data);
            vec3 Tj = get_data<vec3>(node_type, j, 1, data);
            const vec3 Ti = T_pov[static_cast<size_t>(i)];
            return harmonic_tp_density(pj - pi, Ti, Tj, Tp, l1, p1);
          });
    }

    /// Rod harmonic gradient ∇(w Tp) = w′ Tp + w ∇Tp, w = κ ⟨T_i, T_j⟩.
    /// Leaf Tp / ∇Tp are the Cauchy density and its gradient (l0, p0).
    /// κ is the Gaussian mollifier (l1, p1).
    inline std::vector<vec3> tangent_point_harmonic_gradient(
        asawa::rod::rod &R, const std::vector<vec3> &p_pov,
        const std::vector<real> &w_pov, const std::vector<vec3> &T_pov, real l0,
        real p0, real l1, real p1) {
      const std::vector<vec3> Gv =
          tangent_point_gradient_cauchy(R, p_pov, w_pov, T_pov, l0, p0);
      const std::vector<real> Ev =
          tangent_point_energy_cauchy(R, p_pov, w_pov, T_pov, l0, p0);
      std::vector<vec3> Tc = R.N2c();
      return integrate_over_rod<vec3>(
          R, p_pov,
          [&Gv, &Ev, &Tc](const std::vector<index_t> &edge_ids,
                          Rod_Sum_Type &sum) {
            sum.bind(calder::vec3_datum::create(edge_ids, Gv));
            sum.bind(calder::scalar_datum::create(edge_ids, Ev));
            sum.bind(calder::vec3_datum::create(edge_ids, Tc));
          },
          [&T_pov, l1, p1](const index_t i, const index_t j, const vec3 &pi,
                           const vec3 &pj,
                           const std::vector<calder::datum::ptr> &data,
                           Rod_Sum_Type::Node_Type node_type,
                           const Rod_Sum_Type::Tree &) -> vec3 {
            const vec3 grad_Tp = get_data<vec3>(node_type, j, 0, data);
            const real Tp = get_data<real>(node_type, j, 1, data);
            vec3 Tj = get_data<vec3>(node_type, j, 2, data);
            const vec3 Ti = T_pov[static_cast<size_t>(i)];
            return harmonic_tp_pair(pj - pi, Ti, Tj, Tp, grad_Tp, l1, p1);
          });
    }

    /// Mollified / harmonic TP gradient via convolution Φ = κ ∗ T.
    /// Leaf T uses regularized density energy; shell force path uses
    /// post-diff Reg leaf grad (same as tangent_point_gradient).
    /// Shell: vertex query, face-averaged G/E then mollify.
    inline std::vector<vec3>
    tangent_point_harmonic_gradient(asawa::shell::shell &M,
                                    const std::vector<vec3> &p_pov,
                                    const std::vector<real> &w_pov,
                                    const std::vector<vec3> &N_pov, real l0,
                                    real p0, real smooth_l0, real p1 = 2.0)
    {
      std::vector<vec3> &x = asawa::get_vec_data(M, 0);
      // Shell demo leaf: regularized + mollify (force flow, not Armijo).
      std::vector<vec3> Gv =
          tangent_point_gradient(M, p_pov, w_pov, N_pov, l0, p0);
      std::vector<real> Ev =
          tangent_point_energy(M, p_pov, w_pov, N_pov, l0, p0);
      std::vector<vec3> Gf = asawa::shell::vert_to_face<vec3>(M, x, Gv);
      std::vector<real> Ef = asawa::shell::vert_to_face<real>(M, x, Ev);
      std::vector<vec3> Gs = smoothed_gradient(M, p_pov, Gf, smooth_l0, p1);
      std::vector<vec3> Ks = gradient_scalar(M, p_pov, Ef, smooth_l0, p1);
      std::vector<vec3> G(Gs.size(), vec3::Zero());
      for (size_t i = 0; i < G.size(); ++i)
        G[i] = Gs[i] + Ks[i];
      return G;
    }

  } // namespace calder
} // namespace gaudi
#endif
