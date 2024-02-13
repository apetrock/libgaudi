//
//  m2Includes.h
//  Manifold
//
//  Created by John Delaney on 5/22/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//
// includes files that are required for the data structure only, not derived
// files that USE the data structure
#ifndef __SHELL_INTEGRATOR__
#define __SHELL_INTEGRATOR__

#include "integrators.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <vector>

namespace gaudi
{

  namespace calder
  {

    using Shell_Tree_Type = arp::T3;
    using Shell_Sum_Type = calder::fast_summation<Shell_Tree_Type>;

    using Shell_Bind_Fcn =
        std::function<void(const std::vector<index_t> &, Shell_Sum_Type &)>;

    template <typename Q>
    using Shell_Compute_Fcn =
        std::function<Q(const index_t &i, const index_t &j, const vec3 &,
                        const vec3 &, const std::vector<datum::ptr> &,
                        Shell_Sum_Type::Node_Type &, const Shell_Tree_Type::node &,
                        const Shell_Tree_Type &)>;

    template <typename T>
    T get_data(Shell_Sum_Type::Node_Type node_type, index_t j, index_t data_id,
               const std::vector<calder::datum::ptr> &data)
    {
      const typename calder::datum_t<T>::ptr F_datum =
          static_pointer_cast<typename calder::datum_t<T>>(data[data_id]);
      if (node_type == Shell_Sum_Type::Node_Type::LEAF)
      {
        return F_datum->leaf_data()[j];
      }
      else
      {
        return F_datum->node_data()[j];
      }
    }

    template <typename T>
    std::vector<T>
    integrate_over_shell(asawa::shell::shell &M, const std::vector<vec3> &p_pov,
                         Shell_Bind_Fcn bind_fcn = nullptr,
                         Shell_Compute_Fcn<T> compute_fcn = nullptr)
    {

      std::vector<vec3> &x = asawa::get_vec_data(M, 0);
      std::vector<index_t> face_vert_ids = M.get_face_vert_ids();
      std::vector<index_t> face_map = M.get_face_map();
      std::vector<index_t> face_ids = M.get_face_range();
      // std::cout << __PRETTY_FUNCTION__ << std::endl;
      std::cout << "summing" << std::endl;
      std::cout << " -n_faces: " << face_ids.size() << std::endl;
      std::cout << " -create: " << std::endl;

      Shell_Tree_Type::ptr face_tree = arp::T3::create(face_vert_ids, x, 16);
      std::cout << " -sum: " << std::endl;
      Shell_Sum_Type sum(*face_tree);
      bind_fcn(face_ids, sum);

      std::cout << " -compute: " << std::endl;
      std::vector<T> us = sum.calc<T>(
          p_pov,
          [&compute_fcn](const index_t &i, const index_t &j, const vec3 &pi,
                         const std::vector<calder::datum::ptr> &data,
                         Shell_Sum_Type::Node_Type node_type,
                         const Shell_Sum_Type::Node &node,
                         const Shell_Sum_Type::Tree &tree) -> T
          {
            vec3 x0 = tree.vert(3 * j + 0);
            vec3 x1 = tree.vert(3 * j + 1);
            vec3 x2 = tree.vert(3 * j + 2);
            vec3 pj;
            real dist = va::distance_from_triangle({x0, x1, x2}, pi, pj);
            pj = 0.333 * (x0 + x1 + x2);
            return compute_fcn(i, j, pi, pj, data, node_type, node, tree);
          },
          [&compute_fcn](const index_t &i, const index_t &j, const vec3 &pi,
                         const std::vector<calder::datum::ptr> &data,
                         Shell_Sum_Type::Node_Type node_type,
                         const Shell_Sum_Type::Node &node,
                         const Shell_Sum_Type::Tree &tree) -> T
          {
            vec3 pj = node.center();
            return compute_fcn(i, j, pi, pj, data, node_type, node, tree);
          },
          0.25, false);
      return us;
    }

    real calc_w(vec3 dx, real eps, real p)
    {
      // std laplace kernel
      real dist = dx.norm();
      real distp = pow(dist, p);
      real kappa = 1.0 / (distp + eps);
      return kappa;
    };

    template <typename T>
    class shell_integration_bundle
    {
    public:
      using type = T;
      using Manifold_Type = asawa::shell::shell;
      using Tree_Type = Shell_Tree_Type;
      using Sum_Type = Shell_Sum_Type;
      using Bind_Fcn = Shell_Bind_Fcn;
      using Compute_Fcn = Shell_Compute_Fcn<T>;

      static std::vector<T> integrate(Manifold_Type &M, const std::vector<vec3> &p_pov,
                       Bind_Fcn bind_fcn = nullptr,
                       Compute_Fcn compute_fcn = nullptr)
      {
        return integrate_over_shell<T>(M, p_pov, bind_fcn, compute_fcn);
      }
      
    };

    vec3 calc_dw(vec3 dx, real eps, real p)
    {
      real dist = dx.norm();
      real distpm1 = pow(dist, p - 1);
      real distp = pow(dist, p);
      real distp_eps_2 = pow(distp + eps, 2.0);
      return -p * distpm1 / distp_eps_2 * dx / dist;
    };

    // calc w/dw but using gaussian kernel instead of std laplace
    real calc_kg(vec3 dx, real l)
    {
      real dist = dx.norm();
      real distp = pow(dist, 2.0);
      real lp = pow(l, 2.0);
      real C = 1.0 / std::sqrt(2.0 * M_PI);
      real expx2 = exp(-distp / lp);
      real kappa = C / l * expx2;
      // real kappa = exp(-distp / eps);
      return kappa;
    };

    vec3 calc_dkg(vec3 dx, real l)
    {
      real dist = dx.norm();
      real distp = pow(dist, 2.0);
      real lp = pow(l, 2.0);
      real C = 1.0 / std::sqrt(2.0 * M_PI);
      real expx2 = exp(-distp / lp);

      return -2.0 * dx * C / l / lp * expx2;
    };

    std::vector<vec3> smoothed_gradient(asawa::shell::shell &M,
                                        const std::vector<vec3> &p_pov,
                                        const std::vector<vec3> &omega, real l0,
                                        real p = 3.0)
    {
      std::vector<vec3> x = asawa::get_vec_data(M, 0);
      std::vector<real> weights = asawa::shell::face_areas(M, x);
      std::vector<vec3> us = integrate_over_shell<vec3>(
          M, p_pov,
          [&omega, &weights](const std::vector<index_t> &face_ids,
                             Shell_Sum_Type &sum)
          {
            sum.bind(calder::scalar_datum::create(face_ids, weights));
            sum.bind(calder::vec3_datum::create(face_ids, omega));
          },
          [l0, p](const index_t i, const index_t j, //
                  const vec3 &pi, const vec3 &pj,
                  const std::vector<calder::datum::ptr> &data,
                  Shell_Sum_Type::Node_Type node_type, //
                  const Shell_Sum_Type::Node &node,    //
                  const Shell_Sum_Type::Tree &tree) -> vec3
          {
            real wj = get_data<real>(node_type, j, 0, data);
            vec3 w = get_data<vec3>(node_type, j, 1, data);
            vec3 dp = pj - pi;
            // real kappa = calc_w(dp, l0, p);
            real kappa = calc_kg(dp, l0);

            return kappa * w;
          });
      return us;
    }

    std::vector<vec3> gradient_scalar(asawa::shell::shell &M,
                                      const std::vector<vec3> &p_pov,
                                      const std::vector<real> &omega, real l0,
                                      real p = 3.0)
    {

      std::vector<vec3> us = integrate_over_shell<vec3>(
          M, p_pov,
          [&omega](const std::vector<index_t> &face_ids, Shell_Sum_Type &sum)
          {
            sum.bind(calder::scalar_datum::create(face_ids, omega));
          },
          [l0, p](const index_t i, const index_t j, //
                  const vec3 &pi, const vec3 &pj,
                  const std::vector<calder::datum::ptr> &data,
                  Shell_Sum_Type::Node_Type node_type, //
                  const Shell_Sum_Type::Node &node,    //
                  const Shell_Sum_Type::Tree &tree) -> vec3
          {
            real w = get_data<real>(node_type, j, 0, data);
            vec3 dp = pj - pi;
            // vec3 dkappa = calc_dw(dp, l0, p);
            vec3 dkappa = calc_dkg(dp, l0);
            return -dkappa * w;
          });
      return us;
    }

    void log_v(vec3 pi, vec3 e)
    {
      logger::line(pi, pi + 10.0 * e, vec4(0.0, 0.3, 1.0, 1.0));
    }

    void log_v(vec3 pi, real e) {}

    template <typename T>
    std::vector<T> mls_avg(asawa::shell::shell &M, const std::vector<T> &v,
                           const std::vector<vec3> &p_pov, real l0, real p = 3.0)
    {
      // takes a surface and a vector field, v defined on the faces of the surface
      // and returns the mls average of v at the points p_pov
      const std::vector<vec3> &x = asawa::get_vec_data(M, 0);

      std::vector<real> sums(p_pov.size(), 0.0);
      std::vector<real> weights = asawa::shell::face_areas(M, x);
      std::vector<T> wV(v);
      for (int i = 0; i < v.size(); i++)
      {
        wV[i] *= weights[i];
      }

      std::vector<T> us = integrate_over_shell<T>(
          M, p_pov,
          [&wV, &weights](const std::vector<index_t> &face_ids,
                          Shell_Sum_Type &sum)
          {
            sum.bind(calder::datum_t<T>::create(face_ids, wV));
            sum.bind(calder::datum_t<real>::create(face_ids, weights));
          },
          [l0, &sums, p](const index_t i, const index_t j, //
                         const vec3 &pi, const vec3 &pj,
                         const std::vector<calder::datum::ptr> &data,
                         Shell_Sum_Type::Node_Type node_type, //
                         const Shell_Sum_Type::Node &node,    //
                         const Shell_Sum_Type::Tree &tree) -> T
          {
            T e = get_data<T>(node_type, j, 0, data);
            real w = get_data<real>(node_type, j, 1, data);

            real dist = (pj - pi).norm();

            real kappa = computeKg(dist, l0, p);

            // real kappa = computeK(dist, l0, p);
            // if (i == 1250) {
            // logger::line(pi, pj, vec4(1.0, 0.3, 0.3, 1.0));
            // log_v(pj, kappa * e);
            // logger::line(pj, pj + 0.1 * vec3(e), vec4(0.0, 0.3, 1.0, 1.0));
            //}
            sums[i] += w * kappa;
            return kappa * e;
          });
#if 1
      for (int i = 0; i < p_pov.size(); i++)
      {
        if (sums[i] < 1e-16)
          continue;
        us[i] /= sums[i];
      }
#endif
      return us;
    }

#if 0
std::vector<vec3> collision_filter(asawa::shell::shell &M,
                                   const std::vector<vec3> &v, ,
                                   const std::vector<vec3> &v_pov,
                                   const std::vector<vec3> &p_pov, real l0,
                                   real p = 3.0) {
  const std::vector<vec3> &x = asawa::get_vec_data(M, 0);

  std::vector<real> sums(p_pov.size(), 0.0);
  std::vector<real> weights = asawa::shell::face_areas(M, x);
  std::vector<T> wV(v);
  for (int i = 0; i < v.size(); i++) {
    wV[i] *= weights[i];
  }

  std::vector<vec3> us = integrate_over_shell<T>(
      M, p_pov,
      [&wV, &weights, &v_pov](const std::vector<index_t> &face_ids,
                              Shell_Sum_Type &sum) {
        sum.bind(calder::datum_t<T>::create(face_ids, wV));
        sum.bind(calder::datum_t<real>::create(face_ids, weights));
      },
      [l0, &sums, p](const index_t i, const index_t j, //
                     const vec3 &pi, const vec3 &pj,
                     const std::vector<calder::datum::ptr> &data,
                     Shell_Sum_Type::Node_Type node_type, //
                     const Shell_Sum_Type::Node &node,    //
                     const Shell_Sum_Type::Tree &tree) -> vec3 {
        vec3 vj = get_data<T>(node_type, j, 0, data);
        vec3 vi = v_pov[i];

        real dotvjvi = vj.normalized().dot(vi.normalized());
        real w = get_data<real>(node_type, j, 1, data);

        real dist = (pj - pi).norm();

        real kappa = computeKg(dist, l0, p);
        // real kappa = computeK(dist, l0, p);
        // if (i == 1250) {
        // logger::line(pi, pj, vec4(1.0, 0.3, 0.3, 1.0));
        // log_v(pj, kappa * e);
        // logger::line(pj, pj + 0.1 * vec3(e), vec4(0.0, 0.3, 1.0, 1.0));
        //}
        sums[i] += w * kappa;
        return kappa * e;
      });
#if 1
  for (int i = 0; i < p_pov.size(); i++) {
    if (sums[i] < 1e-16)
      continue;
    us[i] /= sums[i];
  }
#endif
  return us;
}
#endif

    std::vector<vec3> vortex_force(asawa::shell::shell &M,
                                   const std::vector<vec3> &p_pov,
                                   const std::vector<vec3> &omega, real l0,
                                   real p = 3.0)
    {

      std::vector<vec3> us = integrate_over_shell<vec3>(
          M, p_pov,
          [&omega](const std::vector<index_t> &edge_ids, Shell_Sum_Type &sum)
          {
            sum.bind(calder::vec3_datum::create(edge_ids, omega));
          },
          [l0, p](const index_t i, const index_t j, //
                  const vec3 &pi, const vec3 &pj,
                  const std::vector<calder::datum::ptr> &data,
                  Shell_Sum_Type::Node_Type node_type, //
                  const Shell_Sum_Type::Node &node,    //
                  const Shell_Sum_Type::Tree &tree) -> vec3
          {
            vec3 w = get_data<vec3>(node_type, j, 0, data);
            vec3 dp = pj - pi;
            real dist = dp.norm();
            // dist = std::max(dist, l0);
            real kappa = computeKm(dist, l0, p);

            return -kappa * dp.cross(w);
          });
      return us;
    }

#if 1
    // use Taubin curvature
    std::vector<mat3> covariant_frame(asawa::shell::shell &M,
                                      const std::vector<vec3> &p_pov, //
                                      real l0, real p = 3.0)
    {
      std::vector<vec3> x = asawa::get_vec_data(M, 0);
      std::vector<real> weights = asawa::shell::face_areas(M, x);
      std::vector<vec3> N = asawa::shell::face_normals(M, x);
      for (int i = 0; i < N.size(); i++)
      {
        if (weights[i] < 1e-16)
          continue;
        N[i] = weights[i] * N[i];
      }
      std::vector<real> sums(p_pov.size(), 0.0);
      std::vector<mat3> us = integrate_over_shell<mat3>(
          M, p_pov,
          [&weights, &N](const std::vector<index_t> &face_ids,
                         Shell_Sum_Type &sum)
          {
            sum.bind(calder::vec3_datum::create(face_ids, N));
          },
          [l0, p, &N, &sums](const index_t i, const index_t j, //
                             const vec3 &pi, const vec3 &pj,
                             const std::vector<calder::datum::ptr> &data,
                             Shell_Sum_Type::Node_Type node_type, //
                             const Shell_Sum_Type::Node &node,    //
                             const Shell_Sum_Type::Tree &tree) -> mat3
          {
            // vec3 Ni = N[i].normalized();
            vec3 Nj = get_data<vec3>(node_type, j, 0, data);
            real wN = Nj.norm();
            Nj /= wN;
            vec3 dp = pj - pi;
            real dist = dp.norm();
            // real w = calc_w(dp, l0, p);
            real w = computeK(dist, l0, p);
            sums[i] += w * wN;
            return w * wN * dp * dp.transpose();
          });

      for (int i = 0; i < us.size(); i++)
      {
        if (sums[i] < 1e-16)
        {
          us[i] = mat3::Zero();
          continue;
        }
        mat3 H = us[i] / sums[i];

        Eigen::JacobiSVD<mat3> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
        mat3 U = svd.matrixU();
        mat3 V = svd.matrixV();
        vec3 s = svd.singularValues();

        us[i] = U * s.asDiagonal();
      }

      return us;
    }
#endif

#if 1
    std::vector<mat3> normal_covariant_frame(asawa::shell::shell &M,
                                             const std::vector<vec3> &p_pov, //
                                             const std::vector<vec3> &N_pov, //
                                             real l0, real p = 3.0)
    {
      std::vector<vec3> x = asawa::get_vec_data(M, 0);
      std::vector<real> weights = asawa::shell::face_areas(M, x);
      std::vector<vec3> N = asawa::shell::face_normals(M, x);
      for (int i = 0; i < N.size(); i++)
      {
        if (weights[i] < 1e-16)
          continue;
        N[i] = weights[i] * N[i];
      }
      std::vector<real> sums(p_pov.size(), 0.0);
      std::vector<mat3> us = integrate_over_shell<mat3>(
          M, p_pov,
          [&weights, &N](const std::vector<index_t> &face_ids,
                         Shell_Sum_Type &sum)
          {
            sum.bind(calder::vec3_datum::create(face_ids, N));
          },
          [l0, p, &N_pov, &sums](const index_t i, const index_t j, //
                                 const vec3 &pi, const vec3 &pj,
                                 const std::vector<calder::datum::ptr> &data,
                                 Shell_Sum_Type::Node_Type node_type, //
                                 const Shell_Sum_Type::Node &node,    //
                                 const Shell_Sum_Type::Tree &tree) -> mat3
          {
            vec3 Ni = N_pov[i];
            vec3 Nj = get_data<vec3>(node_type, j, 0, data);
            real wN = Nj.norm();
            Nj /= wN;

            vec3 dp = pj - pi;
            real Nidp = Ni.dot(dp);
            if (Nidp > 0)
              return mat3::Zero();

            real dist = dp.norm();

            real w = computeK(dist, l0, p);
            sums[i] += w * wN;
            return w * wN * Nj * Nj.transpose();
            // return w * wN * dp * dp.transpose();
          });

      for (int i = 0; i < us.size(); i++)
      {
        if (sums[i] < 1e-16)
        {
          us[i] = mat3::Zero();
          continue;
        }
        mat3 H = us[i] / sums[i];

        Eigen::JacobiSVD<mat3> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
        mat3 U = svd.matrixU();
        mat3 V = svd.matrixV();
        vec3 s = svd.singularValues();

        us[i] = U; // always assume U.col[2] is represents the null vector of the
                   // subspace
        // us[i] = U * s.asDiagonal();
      }

      return us;
    }
#endif

#if 1
    // use Taubin curvature
    std::vector<mat3> taubin_curvature(asawa::shell::shell &M,
                                       const std::vector<vec3> &p_pov, //
                                       const std::vector<vec3> &N_pov, //
                                       real l0, real p = 3.0)
    {
      std::vector<vec3> x = asawa::get_vec_data(M, 0);
      std::vector<real> weights = asawa::shell::face_areas(M, x);
      std::vector<vec3> N = asawa::shell::face_normals(M, x);
      for (int i = 0; i < N.size(); i++)
      {
        if (weights[i] < 1e-16)
          continue;
        N[i] = weights[i] * N[i];
      }
      std::vector<real> sums(p_pov.size(), 0.0);
      std::vector<mat3> us = integrate_over_shell<mat3>(
          M, p_pov,
          [&weights, &N](const std::vector<index_t> &face_ids,
                         Shell_Sum_Type &sum)
          {
            sum.bind(calder::vec3_datum::create(face_ids, N));
          },
          [l0, p, &N_pov, &sums](const index_t i, const index_t j, //
                                 const vec3 &pi, const vec3 &pj,
                                 const std::vector<calder::datum::ptr> &data,
                                 Shell_Sum_Type::Node_Type node_type, //
                                 const Shell_Sum_Type::Node &node,    //
                                 const Shell_Sum_Type::Tree &tree) -> mat3
          {
            vec3 Ni = N_pov[i];
            vec3 Nj = get_data<vec3>(node_type, j, 0, data);
            real wN = Nj.norm();
            Nj /= wN;
            mat3 R = va::rejection_matrix(Ni);
            // mat3 R = mat3::Identity() - Ni * Nj.transpose();
            real Nij = Ni.dot(Nj);
            if (Nij < 0)
              return mat3::Zero();

            vec3 dp = pj - pi;
            vec3 Rdp = R * dp;
            Rdp.normalize();

            real kij = Ni.dot(dp) / dp.dot(dp);

            real dist = dp.norm();
            sums[i] += wN;

            return wN * kij * Rdp * Rdp.transpose();
            // return w * wN * dp * dp.transpose();
          });

      for (int i = 0; i < us.size(); i++)
      {
        if (sums[i] < 1e-16)
        {
          us[i] = mat3::Zero();
          continue;
        }
        mat3 H = us[i] / sums[i];

        Eigen::JacobiSVD<mat3> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
        mat3 U = svd.matrixU();
        mat3 V = svd.matrixV();
        vec3 s = svd.singularValues();

        us[i] = U * s.asDiagonal();
      }

      return us;
    }

#endif

  } // namespace calder
} // namespace gaudi
#endif