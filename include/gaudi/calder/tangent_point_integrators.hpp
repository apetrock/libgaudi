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
#include "rod_integrators.hpp"
#include "shell_integrators.hpp"
#include <cmath>
#include <cstdlib>
#include <ostream>
#include <vector>

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
                                  const Rod_Sum_Type::Node &node,    //
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

            vec3 gj = calc_tangent_point_radius_grad(dp, Nj, l0, p);
            vec3 gi = calc_tangent_point_radius_grad(-dp, Ni, l0, p);

            return 0.5 * (wi * gi - wj * gj);
          });
      return us;
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
                                  const Shell_Sum_Type::Node &node,    //
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
                                  const Shell_Sum_Type::Node &node,    //
                                  const Shell_Sum_Type::Tree &tree) -> vec3
          {
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

            vec3 gj = calc_tangent_point_radius_grad(dp, Nj, l0, pk);
            vec3 gi = calc_tangent_point_radius_grad(-dp, Ni, l0, pk);

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
                                  const Shell_Sum_Type::Node &node,    //
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
            //   logger::line(pi, pj, vec4(0.0, 1.0, 0.5, 1.0));
            //   logger::line(pj, pj + 0.1 * Nj, vec4(0.0, 1.0, 0.5, 1.0));
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
              const Shell_Sum_Type::Node &node,    //
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
            //    logger::line(pi, pj, vec4(0.0, 1.0, 0.5, 1.0));
            //    logger::line(pj, pj + 0.1 * Nj, vec4(0.0, 1.0, 0.5, 1.0));
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

  } // namespace calder
} // namespace gaudi
#endif