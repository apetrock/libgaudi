//
//  m2Includes.h
//  Manifold
//
//  Created by John Delaney on 5/22/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//
// includes files that are required for the data structure only, not derived
// files that USE the data structure
#ifndef __CALDER_LEAST_SQUARES_FIT__
#define __CALDER_LEAST_SQUARES_FIT__

#include "gaudi/common.h"
#include "integrators.hpp"
#include "rod_integrators.hpp"
#include "shell_integrators.hpp"
#include "gaudi/albers/quadric.hpp"
#include "gaudi/albers/sphere.hpp"
#include "gaudi/albers/darboux_cyclide.hpp"
#include "gaudi/albers/circle.hpp"

#include <chrono>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <ostream>
#include <unistd.h>
#include <utility>
#include <vector>

#define GENERATE_WEIGHT_FUNCS(func)                      \
  weight_func<rod_bundle> rod_##func = func<rod_bundle>; \
  weight_func<shell_bundle> shell_##func = func<shell_bundle>;

namespace gaudi
{

  namespace calder
  {

    TYPEDEF_VEC(7)
    TYPEDEF_VEC(10)
    TYPEDEF_MAT(10)
    TYPEDEF_MAT_NM(4, 10)

    using rod_bundle = rod_integration_bundle<real>;
    using shell_bundle = shell_integration_bundle<real>;
    using rod_node_type = rod_bundle::Sum_Type::Node_Type;
    using shell_node_type = shell_bundle::Sum_Type::Node_Type;
    using data_vector = std::vector<calder::datum::ptr>;

    template <typename M_TYPE>
    using weight_func = real (*)(int, int, const std::vector<calder::datum::ptr> &, typename M_TYPE::Sum_Type::Node_Type, const vec3 &, const vec3 &, const vec3 &, real, real);

    template <typename M_TYPE>
    real identity_weight(int i, int j, //
                         const std::vector<calder::datum::ptr> &data,
                         typename M_TYPE::Sum_Type::Node_Type node_type, //
                         const vec3 &dp, const vec3 &Ni, const vec3 &Nj, real l0, real p = 3.0)
    {
      return 1.0;
    }

    template <typename M_TYPE>
    real inv_dist_weight(int i, int j, //
                         const std::vector<calder::datum::ptr> &data,
                         typename M_TYPE::Sum_Type::Node_Type node_type, //
                         const vec3 &dp, const vec3 &Ni, const vec3 &Nj, real l0, real p = 3.0)
    {
      return calc_inv_dist(dp, l0, p);
    }

    template <typename M_TYPE>
    real inv_rad_weight(int i, int j, //
                        const std::vector<calder::datum::ptr> &data,
                        typename M_TYPE::Sum_Type::Node_Type node_type, //
                        const vec3 &dp, const vec3 &Ni, const vec3 &Nj, real l0, real p = 3.0)
    {
      return /*pow(Ni.dot(Nj), p) */ calc_tangent_point_inverse_radius(dp, Nj, l0, p);
    }

    template <typename M_TYPE>
    real convex_weight(int i, int j, //
                       const std::vector<calder::datum::ptr> &data,
                       typename M_TYPE::Sum_Type::Node_Type node_type, //
                       const vec3 &dp, const vec3 &Ni, const vec3 &Nj, real l0, real p = 3.0)
    {
      return dp.dot(Nj) > 0 ? 1.0 : 0.0;
    }

    GENERATE_WEIGHT_FUNCS(identity_weight)
    GENERATE_WEIGHT_FUNCS(inv_dist_weight)
    GENERATE_WEIGHT_FUNCS(inv_rad_weight)
    GENERATE_WEIGHT_FUNCS(convex_weight)

    template <typename F_TYPE, typename M_TYPE>
    std::vector<typename F_TYPE::coefficients> generic_fit(typename M_TYPE::Manifold_Type &M,
                                                           const std::vector<vec3> &Ns,
                                                           const std::vector<vec3> &p_pov,
                                                           const std::vector<vec3> &n_pov, real l0,
                                                           real p = 3.0, //
                                                           weight_func<M_TYPE>
                                                               weight_func = identity_weight<M_TYPE>)
    {

      using type = typename M_TYPE::type;

      std::vector<F_TYPE> accumulators(p_pov.size());

      std::vector<type> us = M_TYPE::integrate(
          M, p_pov,
          [&Ns](const std::vector<index_t> &edge_ids, typename M_TYPE::Sum_Type &sum)
          {
            sum.bind(calder::vec3_datum::create(edge_ids, Ns));
          },
          [&](const index_t i, const index_t j, //
              const vec3 &pi, const vec3 &pj,
              const std::vector<calder::datum::ptr> &data,
              typename M_TYPE::Sum_Type::Node_Type node_type, //
              const typename M_TYPE::Sum_Type::Node &node,    //
              const typename M_TYPE::Sum_Type::Tree &tree) -> type
          {
            vec3 Nj = get_data<vec3>(node_type, j, 0, data);
            vec3 Ni = n_pov[i];

            real wj = Nj.norm();

            Nj /= wj;

            vec3 dp = pj - pi;
            real wf = weight_func(i, j, data, node_type, dp, Ni, Nj, l0, p);
            real kappa = calc_inv_dist(dp, l0, p);
            real w_total = wj * wf;

            if (w_total < 1e-8)
            { // prevent negative weights
              return 0.0;
            }

            accumulators[i].accumulate(w_total, dp, Nj);

            return 0.0;
          });

      std::vector<typename F_TYPE::coefficients> out(p_pov.size(), F_TYPE::coefficients::Zero());
      for (int i = 0; i < accumulators.size(); i++)
      {
        out[i] = accumulators[i].solve();
      }

      return out;
    }

    std::vector<vec10> quadric(asawa::rod::rod &R, const std::vector<vec3> &Nr,
                               const std::vector<vec3> &p_pov,
                               const std::vector<vec3> &N_pov, real l0,
                               real p = 3.0)
    {

      std::vector<real> weights = R.l0();

      std::vector<vec3> Ns = Nr;
      for (int i = 0; i < Ns.size(); i++)
      {
        Ns[i] = weights[i] * Nr[i];
      }

      return generic_fit<albers::quadric, rod_bundle>(
          R, Ns, p_pov, N_pov, l0, p, rod_inv_dist_weight);
    }

    std::vector<vec3> quadric_grad(asawa::rod::rod &R, const std::vector<vec3> &Nr,
                                   const std::vector<vec3> &p_pov,
                                   const std::vector<vec3> &N_pov, real l0,
                                   real p = 3.0)
    {
      std::vector<vec10> Q = quadric(R, Nr, p_pov, N_pov, l0, p);
      std::vector<vec3> out(p_pov.size(), vec3::Zero());
      for (int i = 0; i < Q.size(); i++)
      {
        out[i] = albers::quadric_grad(Q[i], p_pov[i]);
        if (out[i].hasNaN())
        {
          out[i] = vec3::Zero();
        }
      }
      return out;
    }

    std::vector<real> quadric_sdf(asawa::rod::rod &R, const std::vector<vec3> &Nr,
                                  const std::vector<vec3> &p_pov,
                                  const std::vector<vec3> &N_pov, real l0,
                                  real p = 3.0)
    {
      std::vector<vec10> Q = quadric(R, Nr, p_pov, N_pov, l0, p);
      std::vector<real> out(p_pov.size(), 0.0);
      for (int i = 0; i < Q.size(); i++)
      {
        out[i] = Q[i][9]; //* q.segment(6, 3);
        if (std::isnan(out[i]))
        {
          out[i] = 0.0;
        }
      }
      return out;
    }

    std::vector<vec3> quadric_center(asawa::rod::rod &R,
                                     const std::vector<vec3> &Nr,
                                     const std::vector<vec3> &p_pov,
                                     const std::vector<vec3> &n_pov, real l0,
                                     real p = 3.0)
    {
      std::vector<vec10> Q = quadric(R, Nr, p_pov, n_pov, l0, p);
      std::vector<vec3> out(p_pov.size(), vec3::Zero());
      for (int i = 0; i < Q.size(); i++)
      {
        vec3 x = p_pov[i];
        vec10 Qi = Q[i];
        vec3 cen = albers::quadric_center(Q[i]);
        out[i] = cen;
      }
      return out;
    }

    std::vector<vec10> quadric(asawa::shell::shell &M,
                               const std::vector<vec3> &p_pov,
                               const std::vector<vec3> &n_pov, real l0,
                               real p = 3.0)
    {
      const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);
      std::vector<real> weights = asawa::shell::face_areas(M, x);
      std::vector<vec3> Ns = asawa::shell::face_normals(M, x);
      for (int i = 0; i < Ns.size(); i++)
      {
        Ns[i] = weights[i] * Ns[i];
      }
      return generic_fit<albers::quadric, shell_bundle>(
          M, Ns, p_pov, n_pov, l0, p, shell_inv_dist_weight);
    }

    std::vector<std::pair<vec3, mat3>>
    quadric_hessian(asawa::shell::shell &M, const std::vector<vec3> &p_pov,
                    const std::vector<vec3> &n_pov, real l0, real p = 3.0)
    {
      std::vector<vec10> Q = quadric(M, p_pov, n_pov, l0, p);
      auto z = std::make_pair(vec3::Zero(), mat3::Zero());
      std::vector<std::pair<vec3, mat3>> out(p_pov.size(), z);
      for (int i = 0; i < Q.size(); i++)
      {
        vec3 x = p_pov[i];
        vec10 Qi = Q[i];
        vec3 g = albers::quadric_grad(Q[i], p_pov[i]);

        // vec3 Ni = g.normalized();
        vec3 Ni = n_pov[i];
        // logger::line(x, x + 0.05 * Ni, vec4(1.0, 0.0, 1.0, 1.0));
        mat3 P = Ni * Ni.transpose();
        mat3 R = va::rejection_matrix(Ni);
        out[i] = std::make_pair(g, albers::quadric_hessian(Q[i]));
      }
      return out;
    }

    std::vector<mat3> quadric_curvature(asawa::shell::shell &M,
                                        const std::vector<vec3> &p_pov,
                                        const std::vector<vec3> &n_pov, real l0,
                                        real p = 3.0)
    {
      std::vector<std::pair<vec3, mat3>> H =
          quadric_hessian(M, p_pov, n_pov, l0, p);
      std::vector<mat3> out(p_pov.size(), mat3::Zero());
      for (int i = 0; i < H.size(); i++)
      {
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
        if (s.hasNaN())
        {
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
                                     real p = 3.0)
    {
      std::vector<vec10> Q = quadric(M, p_pov, n_pov, l0, p);
      std::vector<vec3> out(p_pov.size(), vec3::Zero());
      for (int i = 0; i < Q.size(); i++)
      {
        vec3 x = p_pov[i];
        vec10 Qi = Q[i];
        vec3 cen = albers::quadric_center(Q[i]);
        vec3 g = albers::quadric_grad(Q[i], p_pov[i]);

        out[i] = cen;
      }
      return out;
    }

    std::vector<vec4> sphere(asawa::shell::shell &M, const std::vector<vec3> &p_pov,
                             const std::vector<vec3> &n_pov, real l0,
                             real p = 3.0)
    {

      const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);
      std::vector<vec3> Ns = asawa::shell::face_normals(M, x);
      std::vector<real> weights = asawa::shell::face_areas(M, x);

      for (int i = 0; i < Ns.size(); i++)
        Ns[i] = weights[i] * Ns[i];

      weight_func<shell_bundle> convex_weight =
          [](int i, int j, const std::vector<calder::datum::ptr> &data, shell_node_type node_type, const vec3 &dp, const vec3 &Ni, const vec3 &Nj, real l0, real p)
      {
        return dp.dot(Nj) > 0 ? 1.0 : 0.0;
      };
      std::vector<albers::sphere> accum(p_pov.size());
      std::vector<vec4> S = generic_fit<albers::sphere, shell_bundle>(
          M, asawa::shell::face_normals(M, x), p_pov, n_pov, l0, p,
          convex_weight);

      for (int i = 0; i < S.size(); i++)
      {
        vec3 x = p_pov[i];
        vec3 dc = S[i].segment(0, 3);
        S[i].segment(0, 3) = x + dc;
        //  logger::line(x, x + dc, vec4(0.0, 1.0, 0.0, 1.0));
      }

      return S;
    }

    std::vector<vec7> cylinder(asawa::shell::shell &M,
                               const std::vector<vec3> &p_pov,
                               const std::vector<vec3> &n_pov, real l0,
                               real p = 3.0)
    {
      // 2D circle regression from:
      // https://www.scribd.com/document/14819165/Regressions-coniques-quadriques-circulaire-spherique
      // with added normal constraint

      const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);
      std::vector<vec3> Nv = asawa::shell::vertex_normals(M, x);

      std::vector<real> weights = asawa::shell::face_areas(M, x);
      // what would be cool would be to use the normal frame to compute
      // a an angle measure between the observer null space and the face null space
      // this would require smearing the normal frame over the face or computing

      std::vector<mat3> Fn = normal_covariant_frame(M, p_pov, n_pov, l0, p);
      std::vector<mat3> Fnv = normal_covariant_frame(M, x, Nv, l0, p);
      std::vector<vec3> Yv(x.size(), vec3::Zero());
      for (int i = 0; i < x.size(); i++)
      {
        Yv[i] = Fnv[i].col(2);
      }
      std::vector<vec3> Yf = asawa::shell::vert_to_face<vec3>(M, x, Yv);

      std::vector<albers::constrained_circle> accum(p_pov.size());
      std::vector<vec3> Ns = asawa::shell::face_normals(M, x);
      for (int i = 0; i < Ns.size(); i++)
      {
        Ns[i] = weights[i] * Ns[i];
        Yf[i] = weights[i] * Yf[i].normalized();
      }

      std::vector<real> us = integrate_over_shell<real>(
          M, p_pov,
          [&](const std::vector<index_t> &face_ids, Shell_Sum_Type &sum)
          {
            sum.bind(calder::scalar_datum::create(face_ids, weights));
            sum.bind(calder::vec3_datum::create(face_ids, Ns));
            sum.bind(calder::vec3_datum::create(face_ids, Yf));
          },
          [&](const index_t i, const index_t j, //
              const vec3 &pi, const vec3 &pj,
              const std::vector<calder::datum::ptr> &data,
              Shell_Sum_Type::Node_Type node_type, //
              const Shell_Sum_Type::Node &node,    //
              const Shell_Sum_Type::Tree &tree) -> real
          {
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
            if (dp2.dot(Ni2) > 0)
            {
              return 0.0;
            }

            real dist = dp.norm();
            real kappa = calc_inv_dist(dp, l0, p);
            real Yij = pow(Yi.dot(Yj), 2.0);
            real w = Yij * wj * kappa;

            accum[i].accumulate(w, dp2, Nj2);

            return 0.0;
          });

      std::vector<vec7> out(p_pov.size(), vec7::Zero());
      for (int i = 0; i < accum.size(); i++)
      {
        vec3 xi = p_pov[i];

        mat3 Fi = Fn[i];
        vec3 Ni = Fi.col(2);
        mat23 P;
        P.row(0) = Fi.col(0);
        P.row(1) = Fi.col(1);

        vec3 z = accum[i].solve();
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

    TYPEDEF_VEC(14)
    TYPEDEF_MAT(14)
    TYPEDEF_MAT_NM(4, 14)

    std::vector<vec14> darboux_cyclide(asawa::rod::rod &R,
                                       const std::vector<vec3> &Nr,
                                       const std::vector<vec3> &p_pov,
                                       const std::vector<vec3> &N_pov, real l0,
                                       real p = 3.0)
    {
      std::vector<real> weights = R.l0();
      std::vector<vec3> Ns = Nr;
      for (int i = 0; i < Ns.size(); i++)
      {
        Ns[i] = weights[i] * Nr[i];
      }
      return generic_fit<albers::darboux_cyclide, rod_bundle>(
          R, Ns, p_pov, N_pov, l0, p, rod_inv_dist_weight);
    }

    std::vector<vec3> darboux_cyclide_grad(asawa::rod::rod &R,
                                           const std::vector<vec3> &Nr,
                                           const std::vector<vec3> &p_pov,
                                           const std::vector<vec3> &N_pov, real l0,
                                           real p = 3.0)
    {
      std::vector<vec14> Q = darboux_cyclide(R, Nr, p_pov, N_pov, l0, p);
      std::vector<vec3> out(p_pov.size(), vec3::Zero());
      for (int i = 0; i < Q.size(); i++)
      {
        out[i] = -Q[i][9] * albers::darboux_grad(Q[i], vec3::Zero()); // double check, this is right...
        if (out[i].hasNaN())
        {
          out[i] = vec3::Zero();
        }
      }
      return out;
    }

    std::vector<real> darboux_cyclide_sdf(asawa::rod::rod &R,
                                          const std::vector<vec3> &Nr,
                                          const std::vector<vec3> &p_pov,
                                          const std::vector<vec3> &N_pov, real l0,
                                          real p = 3.0)
    {
      std::vector<vec14> Q = darboux_cyclide(R, Nr, p_pov, N_pov, l0, p);
      std::vector<real> out(p_pov.size(), 0.0);
      for (int i = 0; i < Q.size(); i++)
      {
        out[i] = Q[i][9]; //* q.segment(6, 3);
        if (std::isnan(out[i]))
        {
          out[i] = 0.0;
        }
      }
      return out;
    }

    std::vector<real> darboux_cyclide_sdf_quadric_centers(asawa::rod::rod &R,
                                                          const std::vector<vec3> &Nr,
                                                          const std::vector<vec3> &p_pov,
                                                          const std::vector<vec3> &N_pov, real l0,
                                                          real p = 3.0)
    {
      std::vector<vec3> cens = quadric_center(R, Nr, p_pov, N_pov, l0, p);
      std::vector<vec14> Q = darboux_cyclide(R, Nr, cens, N_pov, l0, p);
      std::vector<real> out(p_pov.size(), 0.0);
      for (int i = 0; i < Q.size(); i++)
      {
        out[i] = Q[i][9]; //* q.segment(6, 3);
        vec3 dp = p_pov[i] - cens[i];
        out[i] = albers::eval_darboux(Q[i], dp);
        if (std::isnan(out[i]))
        {
          out[i] = 0.0;
        }
      }
      return out;
    }

    std::vector<vec14> darboux_cyclide(asawa::shell::shell &M,
                                       const std::vector<vec3> &p_pov,
                                       const std::vector<vec3> &N_pov, real l0,
                                       real p = 3.0, real w0 = 1e-2)
    {
      const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);
      std::vector<real> weights = asawa::shell::face_areas(M, x);
      std::vector<vec3> Ns = asawa::shell::face_normals(M, x);
      for (int i = 0; i < Ns.size(); i++)
      {
        Ns[i] = weights[i] * Ns[i];
      }
      return generic_fit<albers::darboux_cyclide, shell_bundle>(
          M, Ns, p_pov, N_pov, l0, p, shell_inv_rad_weight);
      return generic_fit<albers::darboux_cyclide, shell_bundle>(
          M, Ns, p_pov, N_pov, l0, p, shell_inv_dist_weight);
    }

    std::vector<vec3> darboux_cyclide_grad(asawa::shell::shell &M,
                                           const std::vector<vec3> &p_pov,
                                           const std::vector<vec3> &N_pov, real l0,
                                           real p = 3.0)
    {
      std::vector<vec14> Q = darboux_cyclide(M, p_pov, N_pov, l0, p, 1e-3);
      std::vector<vec3> out(p_pov.size(), vec3::Zero());
      for (int i = 0; i < Q.size(); i++)
      {
        out[i] = -Q[i][9] * Q[i].segment(6, 3); // double check, this is right...
        if (out[i].hasNaN())
        {
          out[i] = vec3::Zero();
        }
      }
      return out;
    }

    std::vector<real> darboux_cyclide_sdf(asawa::shell::shell &M,
                                          const std::vector<vec3> &p_pov,
                                          const std::vector<vec3> &N_pov, real l0,
                                          real p = 3.0)
    {
      std::vector<vec14> Q = darboux_cyclide(M, p_pov, N_pov, l0, p, 0.0);
      std::vector<real> out(p_pov.size(), 0.0);
      for (int i = 0; i < Q.size(); i++)
      {
        out[i] = Q[i][9]; //* q.segment(6, 3);
        if (std::isnan(out[i]))
        {
          out[i] = 0.0;
        }
      }
      return out;
    }

    std::vector<std::pair<vec3, mat3>>
    darboux_hessian(asawa::shell::shell &M, const std::vector<vec3> &p_pov,
                    const std::vector<vec3> &n_pov, real l0, real p = 3.0)
    {
      std::vector<vec14> Q = darboux_cyclide(M, p_pov, n_pov, l0, p);
      auto z = std::make_pair(vec3::Zero(), mat3::Zero());
      std::vector<std::pair<vec3, mat3>> out(p_pov.size(), z);
      for (int i = 0; i < Q.size(); i++)
      {
        vec3 x = p_pov[i];
        vec14 Qi = Q[i];
        vec3 g = albers::darboux_grad(Q[i], p_pov[i]);

        // vec3 Ni = g.normalized();
        vec3 Ni = n_pov[i];
        // logger::line(x, x + 0.05 * Ni, vec4(1.0, 0.0, 1.0, 1.0));
        mat3 P = Ni * Ni.transpose();
        mat3 R = va::rejection_matrix(Ni);
        out[i] = std::make_pair(g, albers::darboux_hessian(Q[i], vec3::Zero()));
      }
      return out;
    }

    std::vector<mat3> darboux_curvature(asawa::shell::shell &M,
                                        const std::vector<vec3> &p_pov,
                                        const std::vector<vec3> &n_pov, real l0,
                                        real p = 3.0)
    {
      std::vector<std::pair<vec3, mat3>> H =
          darboux_hessian(M, p_pov, n_pov, l0, p);
      std::vector<mat3> out(p_pov.size(), mat3::Zero());
      for (int i = 0; i < H.size(); i++)
      {
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
        if (s.hasNaN())
        {
          continue;
        }
#if 0        
        logger::line(x, x + 0.01 * S.col(0), vec4(1.0, 0.0, 0.0, 1.0));
        logger::line(x, x + 0.01 * S.col(1), vec4(0.0, 1.0, 0.0, 1.0));
        logger::line(x, x + 0.01 * S.col(2), vec4(0.0, 0.0, 1.0, 1.0));
#endif
        out[i] = S;
      }
      return out;
    }

  } // namespace calder
} // namespace gaudi
#endif