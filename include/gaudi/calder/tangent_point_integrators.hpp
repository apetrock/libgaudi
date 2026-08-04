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

      std::vector<vec3> us = integrate_over_rod<vec3>(
          R, p_pov,
          [&weights, &Tc](const std::vector<index_t> &edge_ids, Rod_Sum_Type &sum)
          {
            sum.bind(calder::scalar_datum::create(edge_ids, weights));
            sum.bind(calder::vec3_datum::create(edge_ids, Tc));
          },
          [l0, &w_pov, &T_pov, p](const index_t i, const index_t j, //
                                  const vec3 &pi, const vec3 &pj,
                                  const std::vector<calder::datum::ptr> &data,
                                  Rod_Sum_Type::Node_Type node_type, //
                                  const Rod_Sum_Type::Tree &tree) -> vec3
          {
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
#if 0
            // Rosenhead/Cauchy core: K = |N·dp|^p / (|dp|² + ε²)^p, ε = l0
            vec3 gj = calc_tangent_point_radius_grad_cauchy(dp, Nj, l0, p);
            vec3 gi = calc_tangent_point_radius_grad_cauchy(-dp, Ni, l0, p);
#else
            vec3 gj = calc_tangent_point_radius_grad(dp, Nj, l0, p);
            vec3 gi = calc_tangent_point_radius_grad(-dp, Ni, l0, p);
#endif
            //vec3 g = 0.5 * (wi * gi - wj * gj);
            //std::cout << "g: " << g.transpose() << std::endl;
            return 0.5 * (wi * gi - wj * gj);
          });
      return us;
    }

    /// Matrix-free Hess·v for soft TP density (same pair assembly as gradient).
    /// For each pair: dv = vj−vi, contrib mirrors 0.5*(wi Hi(-dv) − wj Hj dv).
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

            vec3 out = (wi * gi - wj * gj);
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
      hess_alpha = 1.0;
      const bool use_hess = hess_alpha != real(0.0);

      std::vector<vec3> us = integrate_over_rod<vec3>(
          R, p_pov,
          [&weights, &Tc](const std::vector<index_t> &edge_ids,
                          Rod_Sum_Type &sum) {
            sum.bind(calder::scalar_datum::create(edge_ids, weights));
            sum.bind(calder::vec3_datum::create(edge_ids, Tc));
          },
          [R_min, tau, &w_pov, &T_pov, p, hess_alpha, use_hess](
              const index_t i, const index_t j, const vec3 &pi, const vec3 &pj,
              const std::vector<calder::datum::ptr> &data,
              Rod_Sum_Type::Node_Type node_type,
              const Rod_Sum_Type::Tree & /*tree*/) -> vec3 {
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
              //std::cout << "gj" << gj.transpose() << std::endl;;
              //std::cout << "hj" << (Hj * (dp)).transpose() << " " << hess_alpha << std::endl;
              gj -= hess_alpha * (Hj * (dp));
              gi += hess_alpha * (Hi * (dp));
            }
            //float wpar = 1.0 - dp.normalized().dot(Tj);
            vec3 out = 0.5 * (wi * gi - wj * gj);
            if (!std::isfinite(out[0]) || !std::isfinite(out[1]) ||
                !std::isfinite(out[2]))
              return vec3::Zero();
            return out;
          });
      return us;
    }

    std::vector<real> tangent_point_energy(asawa::rod::rod &R,
                                           const std::vector<vec3> &p_pov,
                                           const std::vector<real> &w_pov,
                                           const std::vector<vec3> &T_pov,
                                           real l0, real p = 3.0)
    {
      (void)w_pov;
      (void)T_pov;
      std::vector<real> weights = R.l0();
      std::vector<vec3> Tc = R.N2c();

      return integrate_over_rod<real>(
          R, p_pov,
          [&weights, &Tc](const std::vector<index_t> &edge_ids, Rod_Sum_Type &sum)
          {
            sum.bind(calder::scalar_datum::create(edge_ids, weights));
            sum.bind(calder::vec3_datum::create(edge_ids, Tc));
          },
          [l0, p](const index_t /*i*/, const index_t j, //
                  const vec3 &pi, const vec3 &pj,
                  const std::vector<calder::datum::ptr> &data,
                  Rod_Sum_Type::Node_Type node_type, //
                  const Rod_Sum_Type::Tree & /*tree*/) -> real
          {
            real wj = get_data<real>(node_type, j, 0, data);
            vec3 Tj = get_data<vec3>(node_type, j, 1, data);
            wj = std::max(wj, 1e-6);
            Tj.normalize();

            vec3 dp = pj - pi;
            vec3 Bj = Tj.cross(dp);
            if (Bj.squaredNorm() < 1e-30)
              return 0.0;
            vec3 Nj = Bj.normalized().cross(Tj).normalized();
            #if 0
              return wj * calc_tangent_point_inverse_radius_cauchy(dp, Nj, l0, p);
            #else
              return wj * calc_tangent_point_inverse_radius(dp, Nj, l0, p);
#endif
          });
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

            real jr = calc_tangent_point_inverse_radius(dp, Nj, l0, p);
            real ir = calc_tangent_point_inverse_radius(-dp, Ni, l0, p);
            // not sure if symmetric?
            return wj * jr;
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

    /// Product-rule / harmonic TP force: Gs - Ks for f = w·g.
    /// Shell: vertex query, face-averaged G/E then MLS smooth / ∇scalar.
    inline std::vector<vec3>
    tangent_point_harmonic_gradient(asawa::shell::shell &M,
                                    const std::vector<vec3> &p_pov,
                                    const std::vector<real> &w_pov,
                                    const std::vector<vec3> &N_pov, real l0,
                                    real p0, real smooth_l0, real p1 = 2.0)
    {
      std::vector<vec3> &x = asawa::get_vec_data(M, 0);
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
        G[i] = Gs[i] - Ks[i];
      return G;
    }

    /// Rod product-rule / harmonic TP force (corner-valued G/E, no face lift).
    inline std::vector<vec3>
    tangent_point_harmonic_gradient(asawa::rod::rod &R,
                                    const std::vector<vec3> &p_pov,
                                    const std::vector<real> &w_pov,
                                    const std::vector<vec3> &T_pov, real l0,
                                    real p0, real smooth_l0, real p1 = 2.0)
    {
      std::vector<vec3> Gv =
          tangent_point_gradient(R, p_pov, w_pov, T_pov, l0, p0);
      std::vector<real> Ev =
          tangent_point_energy(R, p_pov, w_pov, T_pov, l0, p0);
      std::vector<vec3> Gs = smoothed_gradient(R, p_pov, Gv, smooth_l0, p1);
      std::vector<vec3> Ks = gradient_scalar(R, p_pov, Ev, smooth_l0, p1);
      std::vector<vec3> sKs = smoothed_gradient(R, p_pov, Ks, smooth_l0, p1);
      std::vector<vec3> G(Gs.size(), vec3::Zero());

      for (size_t i = 0; i < G.size(); ++i)
        G[i] = Gs[i] - sKs[i];

      return G;
    }

  } // namespace calder
} // namespace gaudi
#endif
