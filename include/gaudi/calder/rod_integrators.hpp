#ifndef __ROD_INTEGRATOR__
#define __ROD_INTEGRATOR__

#include "gaudi/common.h"
#include "gaudi/arp/hash_tree.hpp"
#include "integrators.hpp"
#include "gaudi/geometry_logger.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ostream>
#include <vector>

namespace gaudi
{

  namespace calder
  {

    using Rod_Tree_Type = arp::T2;
    using Rod_Sum_Type = calder::fast_summation<Rod_Tree_Type>;

    using Rod_Bind_Fcn =
        std::function<void(const std::vector<index_t> &, Rod_Sum_Type &)>;
    template <typename Q>
    using Rod_Compute_Fcn =
        std::function<Q(const index_t &i, const index_t &j, const vec3 &,
                        const vec3 &, const std::vector<datum::ptr> &,
                        Rod_Sum_Type::Node_Type,
                        const Rod_Tree_Type &)>;

    template <typename T>
    std::vector<T> integrate_over_rod(asawa::rod::rod &R,
                                      const std::vector<vec3> &p_pov,
                                      Rod_Bind_Fcn bind_fcn = nullptr,
                                      Rod_Compute_Fcn<T> compute_fcn = nullptr)
    {
      std::vector<vec3> x = R.xc();

      std::vector<index_t> edge_verts = R.get_edge_vert_ids();
      std::vector<index_t> edge_map = R.get_edge_map();
      auto rverts = R.get_vert_range();
      std::vector<index_t> edge_ids(rverts.begin(), rverts.end());

      std::cout << "summing" << std::endl;
      std::cout << " -n_faces: " << edge_ids.size() << std::endl;
      std::cout << " -create: " << std::endl;
      Rod_Tree_Type::ptr edge_tree = arp::T2::create(edge_verts, x, 12);

      Rod_Sum_Type sum(*edge_tree);
      bind_fcn(edge_ids, sum);
      std::cout << " -compute: " << std::endl;
      std::vector<T> us = sum.template calc<T>(
          p_pov,
          [&compute_fcn](const index_t &i, const index_t &j, const vec3 &pi,
                         const std::vector<calder::datum::ptr> &data,
                         Rod_Sum_Type::Node_Type node_type,
                         const Rod_Tree_Type &tree) -> T
          {
            auto simplex = tree.leaf_simplex(j);
            vec3 x0 = simplex[0], x1 = simplex[1];
            vec3 pj = va::project_on_line(x0, x1, pi);
            return compute_fcn(i, j, pi, pj, data, node_type, tree);
          },
          [&compute_fcn](const index_t &i, const index_t &j, const vec3 &pi,
                         const std::vector<calder::datum::ptr> &data,
                         Rod_Sum_Type::Node_Type node_type,
                         const Rod_Tree_Type &tree) -> T
          {
            const ext::extents_t &ext = tree.bvh_.internal[j];
            vec3 pj = 0.5 * (ext[0] + ext[1]);
            return compute_fcn(i, j, pi, pj, data, node_type, tree);
          },
          0.25, false);
      return us;
    }

    template <typename T>
    class rod_integration_bundle
    {
    public:
      using type = T;
      using Manifold_Type = asawa::rod::rod;
      using Tree_Type = Rod_Tree_Type;
      using Sum_Type = Rod_Sum_Type;
      using Bind_Fcn = Rod_Bind_Fcn;
      using Compute_Fcn = Rod_Compute_Fcn<T>;

      static std::vector<T> integrate(Manifold_Type &M, const std::vector<vec3> &p_pov,
                                      Bind_Fcn bind_fcn = nullptr,
                                      Compute_Fcn compute_fcn = nullptr)
      {
        return integrate_over_rod<T>(M, p_pov, bind_fcn, compute_fcn);
      }
    };

    template <typename T>
    std::vector<T> mls_avg(asawa::rod::rod &R, const std::vector<T> &v,
                           const std::vector<vec3> &p_pov, real l0, real p = 3.0)
    {

      std::vector<real> sums(p_pov.size(), 0.0);
      std::vector<real> weights = R.l0();
      std::vector<T> wV(v);
      for (int i = 0; i < v.size(); i++)
      {
        wV[i] *= weights[i];
      }

      std::vector<T> us = integrate_over_rod<T>(
          R, p_pov,
          [&wV, &v, &weights](const std::vector<index_t> &edge_ids,
                              Rod_Sum_Type &sum)
          {
            sum.bind(calder::vec3_datum::create(edge_ids, wV));
            sum.bind(calder::datum_t<real>::create(edge_ids, weights));
          },
          [l0, &sums, p](const index_t i, const index_t j,
                         const vec3 &pi, const vec3 &pj,
                         const std::vector<calder::datum::ptr> &data,
                         Rod_Sum_Type::Node_Type node_type,
                         const Rod_Sum_Type::Tree &tree) -> T
          {
            T e = get_data<T>(node_type, j, 0, data);
            real w = get_data<real>(node_type, j, 1, data);

            vec3 dp = pj - pi;
            real kappa = calc_gaussian(dp, l0);
            sums[i] += w * kappa;
            return kappa * e;
          });
      int max_count = 0;
      int max_count_i = 0;
      for (int i = 0; i < p_pov.size(); i++)
      {
        if (sums[i] < 1e-16)
          continue;
        us[i] /= sums[i];
      }
      return us;
    }

    std::vector<vec3> coulomb_force(asawa::rod::rod &R,
                                    const std::vector<vec3> &p_pov, real l0,
                                    real p = 3.0)
    {
      std::vector<real> weights = R.l0();
      std::vector<real> sums(p_pov.size(), 0.0);
      std::vector<vec3> us = integrate_over_rod<vec3>(
          R, p_pov,
          [&weights](const std::vector<index_t> &edge_ids, Rod_Sum_Type &sum)
          {
            sum.bind(calder::scalar_datum::create(edge_ids, weights));
          },
          [l0, &sums, p](const index_t i, const index_t j,
                         const vec3 &pi, const vec3 &pj,
                         const std::vector<calder::datum::ptr> &data,
                         Rod_Sum_Type::Node_Type node_type,
                         const Rod_Sum_Type::Tree &tree) -> vec3
          {
            real w = get_data<real>(node_type, j, 0, data);
            vec3 dp = pj - pi;
            real kappa = calc_inv_dist(dp, l0, p);
            return w * kappa * dp;
          });
      return us;
    }

    std::vector<vec3> null_coulomb_force(asawa::rod::rod &R,
                                         const std::vector<vec3> &p_pov, real l0,
                                         real p = 3.0)
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
          [l0, p](const index_t i, const index_t j,
                  const vec3 &pi, const vec3 &pj,
                  const std::vector<calder::datum::ptr> &data,
                  Rod_Sum_Type::Node_Type node_type,
                  const Rod_Sum_Type::Tree &tree) -> vec3
          {
            real w = get_data<real>(node_type, j, 0, data);
            vec3 T = get_data<vec3>(node_type, j, 1, data);
            vec3 dp = pj - pi;
            T.normalize();
            vec3 N = va::rejection_matrix(T) * dp;
            real kappa = calc_inv_dist(dp, l0, p);
            return w * kappa * N;
          });
      return us;
    }

    std::vector<vec3> vortex_force(asawa::rod::rod &R,
                                   const std::vector<vec3> &p_pov,
                                   const std::vector<real> &phi, real l0 = 1e-2,
                                   real p = 4.0)
    {

      std::vector<vec3> &x = R.x();
      std::vector<vec3> T = R.dirs();
      std::vector<real> weights = R.l0();

      index_t i = 0;
      for (auto &t : T)
      {
        t *= phi[i++];
      }

      std::vector<vec3> us = integrate_over_rod<vec3>(
          R, p_pov,
          [&weights, &T](const std::vector<index_t> &edge_ids, Rod_Sum_Type &sum)
          {
            sum.bind(calder::scalar_datum::create(edge_ids, weights));
            sum.bind(calder::vec3_datum::create(edge_ids, T));
          },
          [l0, p](const index_t i, const index_t j,
                  const vec3 &pi, const vec3 &pj,
                  const std::vector<calder::datum::ptr> &data,
                  Rod_Sum_Type::Node_Type node_type,
                  const Rod_Sum_Type::Tree &tree) -> vec3
          {
            real w = get_data<real>(node_type, j, 0, data);
            vec3 T = get_data<vec3>(node_type, j, 1, data);
            vec3 dp = pj - pi;
            real kappa = calc_inv_dist(dp, l0, p);

            return -w * kappa * dp.cross(T);
          });
      return us;
    }

    std::vector<mat3> covariance(asawa::rod::rod &R, const std::vector<vec3> &Nr,
                                 const std::vector<vec3> &p_pov, real r, real l0,
                                 real p = 4.0, bool normalize = true)
    {

      std::vector<vec3> &x = R.__x;
      std::vector<quat> &q = R.__u;

      vec3 u0 = vec3(0.0, 0.0, 1.0);
      std::vector<vec3> ue;
      std::vector<real> ls;
      ue.reserve(q.size());
      for (int i = 0; i < q.size(); i++)
      {
        auto ci = asawa::rod::corner_id(i);
        ue[i] = vec3::Zero();
        if (R.next(ci) == asawa::rod::corner_id(-1))
          continue;
        asawa::rod::consec_t ids = R.consec(ci);
        real l = (x[ids[2]] - x[ids[1]]).norm();
        ls.push_back(l);
        vec3 ui = q[i].normalized() * (l * u0);
        ue.push_back(ui);
      }

      std::vector<index_t> edge_verts = R.get_edge_vert_ids();
      std::vector<index_t> edge_map = R.get_edge_map();
      auto rverts_cov = R.get_vert_range();
      std::vector<index_t> edge_ids(rverts_cov.begin(), rverts_cov.end());

      Rod_Tree_Type::ptr edge_tree = arp::T2::create(edge_verts, x, 12);

      calder::fast_summation<Rod_Tree_Type> sum(*edge_tree);
      sum.bind(calder::edge_frame_datum::create(edge_ids, ue));
      sum.bind(calder::scalar_datum::create(edge_ids, ls));

      std::vector<real> sums(p_pov.size(), 0.0);
      std::vector<mat3> u = sum.calc<mat3>(
          p_pov,
          [&](const index_t &i, const index_t &j, const vec3 &pi,
              const std::vector<calder::datum::ptr> &data,
              Rod_Sum_Type::Node_Type node_type,
              const Rod_Sum_Type::Tree &tree) -> mat3
          {
            const calder::edge_frame_datum::ptr F_datum =
                static_pointer_cast<calder::edge_frame_datum>(data[0]);

            const vec3 &e = F_datum->sorted_leaf_data()[j];

            auto simplex = tree.leaf_simplex(j);
            vec3 x0 = simplex[0], x1 = simplex[1];
            vec3 pj = va::project_on_line(x0, x1, pi);
            pj -= r * Nr[i];
            vec3 dp = pj - pi;
            real dist = va::norm(dp);

            real kappa = calc_inv_dist(dp, l0, p);
            real w = (x0 - x1).norm();
            sums[i] += w * kappa;
            return kappa * w * e * e.transpose();
          },
          [&](const index_t &i, const index_t &j, const vec3 &pi,
              const std::vector<calder::datum::ptr> &data,
              Rod_Sum_Type::Node_Type node_type,
              const Rod_Sum_Type::Tree &tree) -> mat3
          {
            const calder::edge_frame_datum::ptr F_datum =
                static_pointer_cast<calder::edge_frame_datum>(data[0]);
            const mat3 &E = F_datum->node_data()[j];
            const calder::scalar_datum::ptr R_datum =
                static_pointer_cast<calder::scalar_datum>(data[1]);
            const real &w = R_datum->node_data()[j];

            const ext::extents_t &ext = tree.bvh_.internal[j];
            vec3 pj = 0.5 * (ext[0] + ext[1]);
            pj -= r * Nr[i];
            vec3 dp = pj - pi;
            real dist = va::norm(dp);
            real kappa = calc_inv_dist(dp, l0, p);

            sums[i] += w * kappa;
            return kappa * E;
          });

      std::vector<mat3> Us(p_pov.size());
      for (int i = 0; i < p_pov.size(); i++)
      {
        Eigen::JacobiSVD<mat3> svd(u[i], Eigen::ComputeFullU | Eigen::ComputeFullV);
        mat3 U = svd.matrixU();
        vec3 s = svd.singularValues();
        mat3 S = mat3::Zero();

        S.col(0) = s[0] * U.col(0);
        S.col(1) = s[1] * U.col(1);
        S.col(2) = s[2] * U.col(2);
        Us[i] = S;
      }

      return Us;
    }

    void visualize_rod_bvh(asawa::rod::rod &R,
                            const std::vector<vec3> &queries,
                            real eps = 0.25,
                            const vec4 &far_color = vec4(0.5, 0.5, 0.1, 1.0),
                            const vec4 &near_color = vec4(0.1, 0.8, 0.2, 1.0),
                            const vec4 &morton_color = vec4(1.0, 0.5, 0.0, 1.0))
    {
      std::vector<vec3> x = R.xc();
      std::vector<index_t> edge_verts = R.get_edge_vert_ids();
      auto rverts = R.get_vert_range();
      std::vector<index_t> edge_ids(rverts.begin(), rverts.end());

      Rod_Tree_Type::ptr edge_tree = arp::T2::create(edge_verts, x, 12);

      const auto &coms = edge_tree->coms_;
      for (size_t i = 0; i + 1 < coms.size(); i++) {
        const vec3 &a = std::get<1>(coms[i]);
        const vec3 &b = std::get<1>(coms[i + 1]);
        geometry_logger::line(a, b, morton_color);
      }

      Rod_Sum_Type sum(*edge_tree);
      std::vector<real> dummy_weights(edge_ids.size(), 1.0);
      sum.bind(calder::scalar_datum::create(edge_ids, dummy_weights));

      sum.template calc<real>(
          queries,
          [&near_color](const index_t &i, const index_t &j, const vec3 &pi,
                        const std::vector<calder::datum::ptr> &data,
                        Rod_Sum_Type::Node_Type node_type,
                        const Rod_Sum_Type::Tree &tree) -> real
          {
            const ext::extents_t &e = tree.bvh_.leaf[j];
            geometry_logger::ext(e[0], e[1], near_color);
            return 0.0;
          },
          [&far_color](const index_t &i, const index_t &j, const vec3 &pi,
                       const std::vector<calder::datum::ptr> &data,
                       Rod_Sum_Type::Node_Type node_type,
                       const Rod_Sum_Type::Tree &tree) -> real
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
