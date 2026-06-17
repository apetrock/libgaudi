#ifndef __SHELL_INTEGRATOR__
#define __SHELL_INTEGRATOR__

#include "gaudi/arp/arp.h"
#include "gaudi/arp/hash_tree.hpp"
#include "gaudi/arp/datums.hpp"
#include "integrators.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <vector>
#include "gaudi/geometry_logger.hpp"

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
                        Shell_Sum_Type::Node_Type,
                        const Shell_Tree_Type &)>;

    template <typename T>
    std::vector<T>
    integrate_over_shell(asawa::shell::shell &M, const std::vector<vec3> &p_pov,
                         Shell_Bind_Fcn bind_fcn = nullptr,
                         Shell_Compute_Fcn<T> compute_fcn = nullptr)
    {

      std::vector<vec3> &x = asawa::get_vec_data(M, 0);
      std::vector<index_t> face_vert_ids = M.get_face_vert_ids();
      std::vector<index_t> face_map = M.get_face_map();
      auto face_ids_typed = M.get_face_range();
      std::vector<index_t> face_ids(face_ids_typed.begin(), face_ids_typed.end());
      std::cout << "summing" << std::endl;
      std::cout << " -n_faces: " << face_ids.size() << std::endl;
      std::cout << " -create: " << std::endl;

      Shell_Tree_Type::ptr face_tree = arp::T3::create(face_vert_ids, x, 16);
      std::cout << " -sum: " << std::endl;
      Shell_Sum_Type sum(*face_tree);
      bind_fcn(face_ids, sum);

      std::vector<real> areas = asawa::shell::face_areas(M, x);
      std::vector<vec3> centroids = asawa::shell::face_centers(M, x);
      auto com = calder::com_datum::create(face_ids, areas, centroids);
      com->pyramid(*face_tree);

      std::cout << " -compute: " << std::endl;
      std::vector<T> us = sum.template calc<T>(
          p_pov,
          [&compute_fcn, &com](const index_t &i, const index_t &j, const vec3 &pi,
                         const std::vector<calder::datum::ptr> &data,
                         Shell_Sum_Type::Node_Type node_type,
                         const Shell_Sum_Type::Tree &tree) -> T
          {
            vec3 pj = com->get_leaf_com(j);
            return compute_fcn(i, j, pi, pj, data, node_type, tree);
          },
          [&compute_fcn, &com](const index_t &i, const index_t &j, const vec3 &pi,
                         const std::vector<calder::datum::ptr> &data,
                         Shell_Sum_Type::Node_Type node_type,
                         const Shell_Sum_Type::Tree &tree) -> T
          {
            vec3 pj = com->get_node_com(j);
            return compute_fcn(i, j, pi, pj, data, node_type, tree);
          },
          0.25, false);
      return us;
    }

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
          [l0, p](const index_t i, const index_t j,
                  const vec3 &pi, const vec3 &pj,
                  const std::vector<calder::datum::ptr> &data,
                  Shell_Sum_Type::Node_Type node_type,
                  const Shell_Sum_Type::Tree &tree) -> vec3
          {
            real wj = get_data<real>(node_type, j, 0, data);
            vec3 w = get_data<vec3>(node_type, j, 1, data);
            vec3 dp = pj - pi;
            real kappa = calc_gaussian(dp, l0);

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
          [l0, p](const index_t i, const index_t j,
                  const vec3 &pi, const vec3 &pj,
                  const std::vector<calder::datum::ptr> &data,
                  Shell_Sum_Type::Node_Type node_type,
                  const Shell_Sum_Type::Tree &tree) -> vec3
          {
            real w = get_data<real>(node_type, j, 0, data);
            vec3 dp = pj - pi;
            vec3 dkappa = calc_d_gaussian(dp, l0);
            return -dkappa * w;
          });
      return us;
    }

    void log_v(vec3 pi, vec3 e)
    {
      geometry_logger::line(pi, pi + 10.0 * e, vec4(0.0, 0.3, 1.0, 1.0));
    }

    void log_v(vec3 pi, real e) {}

    template <typename T>
    std::vector<T> mls_avg(asawa::shell::shell &M, const std::vector<T> &v,
                           const std::vector<vec3> &p_pov, real l0, real p = 3.0)
    {
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
          [l0, &sums, p](const index_t i, const index_t j,
                         const vec3 &pi, const vec3 &pj,
                         const std::vector<calder::datum::ptr> &data,
                         Shell_Sum_Type::Node_Type node_type,
                         const Shell_Sum_Type::Tree &tree) -> T
          {
            T e = get_data<T>(node_type, j, 0, data);
            real w = get_data<real>(node_type, j, 1, data);

            vec3 dp = pj - pi;
            real kappa = calc_gaussian(dp, l0);

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

    template <typename KERNEL_FUNC>
    std::vector<mat3> covariant_vector_frame(
        asawa::shell::shell &M, const std::vector<vec3> &v,
        const std::vector<vec3> &p_pov, real l0, KERNEL_FUNC kernel_func)
    {
      const std::vector<vec3> &x = asawa::get_vec_data(M, 0);
      std::vector<real> weights = asawa::shell::face_areas(M, x);
      std::vector<mat3> frames(v.size(), mat3::Zero());
      std::vector<real> active(v.size(), 0.0);
      for (int i = 0; i < static_cast<int>(v.size()); i++)
      {
        if (weights[i] < 1e-16 || v[i].norm() < 1e-12)
        {
          continue;
        }
        const vec3 axis = v[i].normalized();
        frames[i] = weights[i] * axis * axis.transpose();
        active[i] = weights[i];
      }

      std::vector<real> sums(p_pov.size(), 0.0);
      std::vector<mat3> us = integrate_over_shell<mat3>(
          M, p_pov,
          [&frames, &active](const std::vector<index_t> &face_ids,
                             Shell_Sum_Type &sum)
          {
            sum.bind(calder::datum_t<mat3>::create(face_ids, frames));
            sum.bind(calder::datum_t<real>::create(face_ids, active));
          },
          [l0, &sums, kernel_func](const index_t i, const index_t j,
                                   const vec3 &pi, const vec3 &pj,
                                   const std::vector<calder::datum::ptr> &data,
                                   Shell_Sum_Type::Node_Type node_type,
                                   const Shell_Sum_Type::Tree &tree) -> mat3
          {
            (void)tree;
            const mat3 frame = get_data<mat3>(node_type, j, 0, data);
            const real w = get_data<real>(node_type, j, 1, data);
            if (w < 1e-16)
            {
              return mat3::Zero();
            }
            const vec3 dp = pj - pi;
            const real kappa = kernel_func(dp, l0);
            sums[i] += w * kappa;
            return kappa * frame;
          });

      for (int i = 0; i < static_cast<int>(us.size()); i++)
      {
        if (sums[i] < 1e-16)
        {
          us[i] = mat3::Zero();
          continue;
        }
        const mat3 C = us[i] / sums[i];
        Eigen::SelfAdjointEigenSolver<mat3> es(C);
        us[i] = es.info() == Eigen::Success ? es.eigenvectors() : mat3::Zero();
      }
      return us;
    }

    std::vector<mat3> gaussian_covariant_vector_frame(
        asawa::shell::shell &M, const std::vector<vec3> &v,
        const std::vector<vec3> &p_pov, real l0)
    {
      return covariant_vector_frame(
          M, v, p_pov, l0,
          [](const vec3 &dp, real l0) -> real
          {
            return calc_gaussian(dp, l0);
          });
    }

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
          [l0, p](const index_t i, const index_t j,
                  const vec3 &pi, const vec3 &pj,
                  const std::vector<calder::datum::ptr> &data,
                  Shell_Sum_Type::Node_Type node_type,
                  const Shell_Sum_Type::Tree &tree) -> vec3
          {
            vec3 w = get_data<vec3>(node_type, j, 0, data);
            vec3 dp = pj - pi;
            real kappa = calc_mollified(dp, l0, p);

            return -kappa * dp.cross(w);
          });
      return us;
    }

    std::vector<mat3> covariant_frame(asawa::shell::shell &M,
                                      const std::vector<vec3> &p_pov,
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
            sum.bind(calder::scalar_datum::create(face_ids, weights));
          },
          [l0, p, &N, &sums](const index_t i, const index_t j,
                             const vec3 &pi, const vec3 &pj,
                             const std::vector<calder::datum::ptr> &data,
                             Shell_Sum_Type::Node_Type node_type,
                             const Shell_Sum_Type::Tree &tree) -> mat3
          {
            vec3 Nj = get_data<vec3>(node_type, j, 0, data); // sum(area * n)
            real wN = get_data<real>(node_type, j, 1, data); // sum(area)
            if (Nj.norm() < 1e-12)
              return mat3::Zero();
            Nj.normalize();
            vec3 dp = pj - pi;
            real w = calc_inv_dist(dp, l0, p);
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
        vec3 s = svd.singularValues();

        us[i] = U * s.asDiagonal();
      }

      return us;
    }

    std::vector<mat3> normal_covariant_frame(asawa::shell::shell &M,
                                             const std::vector<vec3> &p_pov,
                                             const std::vector<vec3> &N_pov,
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
            sum.bind(calder::scalar_datum::create(face_ids, weights));
          },
          [l0, p, &N_pov, &sums](const index_t i, const index_t j,
                                 const vec3 &pi, const vec3 &pj,
                                 const std::vector<calder::datum::ptr> &data,
                                 Shell_Sum_Type::Node_Type node_type,
                                 const Shell_Sum_Type::Tree &tree) -> mat3
          {
            vec3 Ni = N_pov[i];
            vec3 Nj = get_data<vec3>(node_type, j, 0, data); // sum(area * n)
            real wN = get_data<real>(node_type, j, 1, data); // sum(area)
            if (Nj.norm() < 1e-12)
              return mat3::Zero();
            Nj.normalize();

            vec3 dp = pj - pi;
            real Nidp = Ni.dot(dp);
            if (Nidp > 0)
              return mat3::Zero();

            real w = calc_inv_dist(dp, l0, p);
            sums[i] += w * wN;
            return w * wN * Nj * Nj.transpose();
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
        vec3 s = svd.singularValues();

        us[i] = U;
      }

      return us;
    }

    std::vector<mat3> taubin_curvature(asawa::shell::shell &M,
                                       const std::vector<vec3> &p_pov,
                                       const std::vector<vec3> &N_pov,
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
            sum.bind(calder::scalar_datum::create(face_ids, weights));
          },
          [l0, p, &N_pov, &sums](const index_t i, const index_t j,
                                 const vec3 &pi, const vec3 &pj,
                                 const std::vector<calder::datum::ptr> &data,
                                 Shell_Sum_Type::Node_Type node_type,
                                 const Shell_Sum_Type::Tree &tree) -> mat3
          {
            vec3 Ni = N_pov[i];
            vec3 Nj = get_data<vec3>(node_type, j, 0, data); // sum(area * n)
            real wN = get_data<real>(node_type, j, 1, data); // sum(area)
            if (Nj.norm() < 1e-12)
              return mat3::Zero();
            Nj.normalize();
            mat3 R = va::rejection_matrix(Ni);
            real Nij = Ni.dot(Nj);
            if (Nij < 0)
              return mat3::Zero();

            vec3 dp = pj - pi;
            vec3 Rdp = R * dp;
            Rdp.normalize();

            real kij = Ni.dot(dp) / dp.dot(dp);

            sums[i] += wN;

            return wN * kij * Rdp * Rdp.transpose();
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
        vec3 s = svd.singularValues();

        us[i] = U * s.asDiagonal();
      }

      return us;
    }

    void visualize_shell_bvh(asawa::shell::shell &M,
                             const std::vector<vec3> &queries,
                             real eps = 0.25,
                             const vec4 &far_color = vec4(0.5, 0.5, 0.1, 1.0),
                             const vec4 &near_color = vec4(0.1, 0.8, 0.2, 1.0),
                             const vec4 &morton_color = vec4(1.0, 0.5, 0.0, 1.0))
    {
      std::vector<vec3> &x = asawa::get_vec_data(M, 0);
      std::vector<index_t> face_vert_ids = M.get_face_vert_ids();
      auto face_ids_typed = M.get_face_range();
      std::vector<index_t> face_ids(face_ids_typed.begin(), face_ids_typed.end());

      Shell_Tree_Type::ptr face_tree = arp::T3::create(face_vert_ids, x, 16);

      // Step 1: Morton-sorted COM polyline (iterates all leaves in sorted order)
      const auto &coms = face_tree->coms_;
      for (size_t i = 0; i + 1 < coms.size(); i++) {
        const vec3 &a = std::get<1>(coms[i]);
        const vec3 &b = std::get<1>(coms[i + 1]);
        geometry_logger::line(a, b, morton_color);
      }

      // Step 2: BH bbox visualization via production sum.calc pathway
      Shell_Sum_Type sum(*face_tree);
      std::vector<real> dummy_weights(face_ids.size(), 1.0);
      sum.bind(calder::scalar_datum::create(face_ids, dummy_weights));

      sum.template calc<real>(
          queries,
          [&near_color](const index_t &i, const index_t &j, const vec3 &pi,
                        const std::vector<calder::datum::ptr> &data,
                        Shell_Sum_Type::Node_Type node_type,
                        const Shell_Sum_Type::Tree &tree) -> real
          {
            const ext::extents_t &e = tree.bvh_.leaf[j];
            geometry_logger::ext(e[0], e[1], near_color);
            return 0.0;
          },
          [&far_color](const index_t &i, const index_t &j, const vec3 &pi,
                       const std::vector<calder::datum::ptr> &data,
                       Shell_Sum_Type::Node_Type node_type,
                       const Shell_Sum_Type::Tree &tree) -> real
          {
            const ext::extents_t &e = tree.bvh_.internal[j];
            geometry_logger::ext(e[0], e[1], far_color);
            return 0.0;
          },
          eps, false);
    }

  } // namespace calder
} // namespace gaudi
#endif
