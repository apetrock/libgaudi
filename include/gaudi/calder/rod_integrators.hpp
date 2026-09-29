#ifndef __ROD_INTEGRATOR__
#define __ROD_INTEGRATOR__

#include "gaudi/common.h"
#include "gaudi/arp/arp.h"
#include "gaudi/arp/edge_cov_tree.hpp"
#include "gaudi/arp/hash_tree.hpp"
#include "integrators.hpp"
#include <memory>
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
                                      Rod_Compute_Fcn<T> compute_fcn = nullptr,
                                      real bh_eps = 0.25)
    {
      // Use vertex positions (__x), not edge midpoints (xc). Adjacency is
      // (i, next(i)), so leaves are the true rod edges that attributes bind to.
      // bh_eps <= 0 forces leaves-only (opening never accepts BRANCH).
      const std::vector<vec3> &x = R.x();

      std::vector<index_t> edge_verts = R.get_edge_vert_ids();
      std::vector<index_t> edge_map = R.get_edge_map();
      auto rverts = R.get_vert_range();
      std::vector<index_t> edge_ids(rverts.begin(), rverts.end());

      Rod_Tree_Type::ptr edge_tree = arp::T2::create(edge_verts, x, 12);

      Rod_Sum_Type sum(*edge_tree);
      bind_fcn(edge_ids, sum);
      std::vector<T> us = sum.template calc<T>(
          p_pov,
          [&compute_fcn](const index_t &i, const index_t &j, const vec3 &pi,
                         const std::vector<calder::datum::ptr> &data,
                         Rod_Sum_Type::Node_Type node_type,
                         const Rod_Tree_Type &tree) -> T
          {
            auto simplex = tree.leaf_simplex(j);
            vec3 x0 = simplex[0], x1 = simplex[1];

            if((x0 - x1).norm() < 1e-8) return z::zero<T>();

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
          bh_eps, false);
      return us;
    }

    // Datum binder for the covariance walk. bind() matches fast_summation.
    struct cov_rod_sum {
      std::vector<datum::ptr> data;
      void bind(const datum::ptr &x) { data.push_back(x); }
    };

    inline arp::edge_cov_tree<2>::ptr build_rod_cov_tree(asawa::rod::rod &R) {
      const std::vector<vec3> &x = R.x();
      const std::vector<index_t> edges = R.get_edge_vert_ids();
      return arp::edge_cov_tree<2>::create(edges, x);
    }

    inline std::vector<index_t> rod_cov_edge_ids(const asawa::rod::rod &R) {
      auto rverts = R.get_vert_range();
      return std::vector<index_t>(rverts.begin(), rverts.end());
    }

    // One count per integrate_over_cov_tree call, not per query.
    inline int &cov_walk_count_slot() {
      static int n = 0;
      return n;
    }
    inline int cov_walk_count() { return cov_walk_count_slot(); }

    // cone_eps < 0 uses the box ratio for the tangent chord as well.
    inline real cov_cone_eps(real box_eps, real cone_eps) {
      return cone_eps < real(0) ? box_eps : cone_eps;
    }

    // Point queries against an existing rod edge_cov_tree. Near and far are
    // the pair kernels. The helper fills the sample: a leaf is the query
    // projected onto that edge; a far node is μ_c, with the tangent datum's
    // node value set to μ_T so the kernel reads the axis. A non-positive box
    // or cone threshold never opens. The tree is not rebuilt.
    template <typename T, class Bind, class Near, class Far>
    std::vector<T>
    integrate_over_cov_tree(const arp::edge_cov_tree<2> &tree,
                            const std::vector<index_t> &edge_ids,
                            const std::vector<vec3> &p_pov, Bind &&bind_fcn,
                            Near &&near, Far &&far, real bh_eps = 0.25,
                            int tangent_slot = 1, real cone_eps = -1) {
      std::vector<T> us(p_pov.size(), z::zero<T>());
      if (tree.internal_nodes_.empty())
        return us;

      cov_rod_sum sum;
      bind_fcn(edge_ids, sum);
      for (datum::ptr &d : sum.data)
        d->pyramid(tree);
      if (tangent_slot >= 0 &&
          tangent_slot < static_cast<int>(sum.data.size())) {
        auto tangent = std::dynamic_pointer_cast<vec3_datum>(
            sum.data[static_cast<size_t>(tangent_slot)]);
        if (tangent && tangent->node_data().size() == tree.mu_T_.size())
          tangent->node_data() = tree.mu_T_;
      }

      const std::vector<datum::ptr> &data = sum.data;
      const real cone = cov_cone_eps(bh_eps, cone_eps);
      ++cov_walk_count_slot();
#pragma omp parallel for
      for (int i = 0; i < static_cast<int>(p_pov.size()); ++i) {
        T acc = z::zero<T>();
        arp::traverse_cov_opening(
            tree, static_cast<index_t>(i), p_pov[static_cast<size_t>(i)],
            bh_eps, cone,
            [&](index_t node_id, index_t qi, const vec3 &pi) {
              acc += far(qi, node_id, pi, tree.mu_c(node_id), data, BRANCH,
                         tree);
            },
            [&](index_t leaf_id, index_t, index_t qi, const vec3 &pi) {
              const index_t prim = tree.get_index(static_cast<size_t>(leaf_id));
              const index_t ia = tree.adjacency_[static_cast<size_t>(2 * prim)];
              const index_t ib =
                  tree.adjacency_[static_cast<size_t>(2 * prim + 1)];
              const vec3 &x0 = tree.data_[static_cast<size_t>(ia)];
              const vec3 &x1 = tree.data_[static_cast<size_t>(ib)];
              if ((x0 - x1).norm() < real(1e-8))
                return;
              const vec3 pj = va::project_on_line(x0, x1, pi);
              acc += near(qi, leaf_id, pi, pj, data, LEAF, tree);
            });
        us[static_cast<size_t>(i)] = acc;
      }
      return us;
    }

    template <typename T, class Bind, class Near, class Far>
    std::vector<T> integrate_over_rod_cov(asawa::rod::rod &R,
                                          const std::vector<vec3> &p_pov,
                                          Bind &&bind_fcn, Near &&near,
                                          Far &&far, real bh_eps = 0.25,
                                          int tangent_slot = 1,
                                          real cone_eps = -1) {
      arp::edge_cov_tree<2>::ptr tree = build_rod_cov_tree(R);
      if (!tree)
        return std::vector<T>(p_pov.size(), z::zero<T>());
      return integrate_over_cov_tree<T>(*tree, rod_cov_edge_ids(R), p_pov,
                                        std::forward<Bind>(bind_fcn),
                                        std::forward<Near>(near),
                                        std::forward<Far>(far), bh_eps,
                                        tangent_slot, cone_eps);
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
                                      Compute_Fcn compute_fcn = nullptr,
                                      real bh_eps = 0.25)
      {
        return integrate_over_rod<T>(M, p_pov, bind_fcn, compute_fcn, bh_eps);
      }
    };

    /// MLS-smoothed vector field (shell `smoothed_gradient` rod port).
    inline std::vector<vec3> smoothed_gradient(asawa::rod::rod &R,
                                               const std::vector<vec3> &p_pov,
                                               const std::vector<vec3> &omega,
                                               real l0, real p = 3.0) {
      (void)p;
      std::vector<real> weights = R.l0();
      return integrate_over_rod<vec3>(
          R, p_pov,
          [&omega, &weights](const std::vector<index_t> &edge_ids,
                             Rod_Sum_Type &sum) {
            sum.bind(calder::scalar_datum::create(edge_ids, weights));
            sum.bind(calder::vec3_datum::create(edge_ids, omega));
          },
          [l0](const index_t /*i*/, const index_t j, const vec3 &pi,
               const vec3 &pj, const std::vector<calder::datum::ptr> &data,
               Rod_Sum_Type::Node_Type node_type,
               const Rod_Sum_Type::Tree & /*tree*/) -> vec3 {
            (void)get_data<real>(node_type, j, 0, data);
            vec3 w = get_data<vec3>(node_type, j, 1, data);
            vec3 dp = pj - pi;
            real kappa = calc_gaussian(dp, l0);
            // Match shell smoothed_gradient: kappa * omega (area bound for pyramid).
            return kappa * w;
          });
    }

    /// MLS-smoothed scalar field (κ ∗ ω).
    inline std::vector<real> smoothed_scalar(asawa::rod::rod &R,
                                             const std::vector<vec3> &p_pov,
                                             const std::vector<real> &omega,
                                             real l0, real p = 3.0) {
      (void)p;
      return integrate_over_rod<real>(
          R, p_pov,
          [&omega](const std::vector<index_t> &edge_ids, Rod_Sum_Type &sum) {
            sum.bind(calder::scalar_datum::create(edge_ids, omega));
          },
          [l0](const index_t /*i*/, const index_t j, const vec3 &pi,
               const vec3 &pj, const std::vector<calder::datum::ptr> &data,
               Rod_Sum_Type::Node_Type node_type,
               const Rod_Sum_Type::Tree & /*tree*/) -> real {
            real w = get_data<real>(node_type, j, 0, data);
            vec3 dp = pj - pi;
            return calc_gaussian(dp, l0) * w;
          });
    }

    /// One BH walk: f_k = W(ℓ_k) ∗ ω for ℓ_k = l0 · 2^k, k = 0..n_octaves-1.
    /// Unnormalized (raw κ∗ω). Prefer mls_scalar_octaves for partition-of-unity.
    /// Returns per-POV Eigen vec8 (only first n_octaves entries used; n ≤ 8).
    inline std::vector<Eigen::Matrix<real, 8, 1>>
    smoothed_scalar_octaves(asawa::rod::rod &R, const std::vector<vec3> &p_pov,
                            const std::vector<real> &omega, real l0,
                            int n_octaves, real bh_eps = 0.25) {
      using oct8 = Eigen::Matrix<real, 8, 1>;
      const int n_oct = std::clamp(n_octaves, 1, 8);
      std::vector<real> lengths(static_cast<size_t>(n_oct));
      real lk = std::max(l0, real(1e-16));
      for (int k = 0; k < n_oct; ++k) {
        lengths[static_cast<size_t>(k)] = lk;
        lk *= real(2.0);
      }

      return integrate_over_rod<oct8>(
          R, p_pov,
          [&omega](const std::vector<index_t> &edge_ids, Rod_Sum_Type &sum) {
            sum.bind(calder::scalar_datum::create(edge_ids, omega));
          },
          [n_oct, &lengths](const index_t /*i*/, const index_t j, const vec3 &pi,
                            const vec3 &pj,
                            const std::vector<calder::datum::ptr> &data,
                            Rod_Sum_Type::Node_Type node_type,
                            const Rod_Sum_Type::Tree & /*tree*/) -> oct8 {
            const real w = get_data<real>(node_type, j, 0, data);
            const vec3 dp = pj - pi;
            oct8 out = oct8::Zero();
            for (int k = 0; k < n_oct; ++k)
              out[k] = calc_gaussian(dp, lengths[static_cast<size_t>(k)]) * w;
            return out;
          },
          bh_eps);
    }

    /// MLS / partition-of-unity octaves:
    ///   f_k[i] = (κ_k ∗ (ℓ₀ ω))_i / (κ_k ∗ ℓ₀)_i
    /// so each scale is a true local average (no spatially varying gain).
    inline std::vector<Eigen::Matrix<real, 8, 1>>
    mls_scalar_octaves(asawa::rod::rod &R, const std::vector<vec3> &p_pov,
                       const std::vector<real> &omega, real l0, int n_octaves,
                       real bh_eps = 0.25) {
      using oct8 = Eigen::Matrix<real, 8, 1>;
      const int n_oct = std::clamp(n_octaves, 1, 8);
      std::vector<real> lengths(static_cast<size_t>(n_oct));
      real lk = std::max(l0, real(1e-16));
      for (int k = 0; k < n_oct; ++k) {
        lengths[static_cast<size_t>(k)] = lk;
        lk *= real(2.0);
      }

      std::vector<real> weights = R.l0();
      std::vector<real> w_omega(omega);
      const size_t n_w = std::min(w_omega.size(), weights.size());
      for (size_t i = 0; i < n_w; ++i)
        w_omega[i] *= weights[i];

      std::vector<oct8> mass(p_pov.size(), oct8::Zero());
      std::vector<oct8> us = integrate_over_rod<oct8>(
          R, p_pov,
          [&w_omega, &weights](const std::vector<index_t> &edge_ids,
                               Rod_Sum_Type &sum) {
            sum.bind(calder::scalar_datum::create(edge_ids, w_omega));
            sum.bind(calder::scalar_datum::create(edge_ids, weights));
          },
          [n_oct, &lengths, &mass](const index_t i, const index_t j,
                                   const vec3 &pi, const vec3 &pj,
                                   const std::vector<calder::datum::ptr> &data,
                                   Rod_Sum_Type::Node_Type node_type,
                                   const Rod_Sum_Type::Tree & /*tree*/) -> oct8 {
            const real wo = get_data<real>(node_type, j, 0, data);
            const real w = get_data<real>(node_type, j, 1, data);
            const vec3 dp = pj - pi;
            oct8 out = oct8::Zero();
            for (int k = 0; k < n_oct; ++k) {
              const real kappa =
                  calc_gaussian(dp, lengths[static_cast<size_t>(k)]);
              mass[static_cast<size_t>(i)][k] += kappa * w;
              out[k] = kappa * wo;
            }
            return out;
          },
          bh_eps);

      for (size_t i = 0; i < us.size(); ++i) {
        for (int k = 0; k < n_oct; ++k) {
          const real m = mass[i][k];
          us[i][k] = (m > real(1e-16)) ? us[i][k] / m : real(0.0);
        }
      }
      return us;
    }

    /// Band residuals R_k = rms(f_k − f_{k+1}) over POVs; size n_octaves-1.
    inline std::vector<real>
    octave_band_rms(const std::vector<Eigen::Matrix<real, 8, 1>> &f,
                    int n_octaves) {
      const int n_oct = std::clamp(n_octaves, 1, 8);
      std::vector<real> R(static_cast<size_t>(std::max(n_oct - 1, 0)),
                          real(0.0));
      if (f.empty() || n_oct < 2)
        return R;
      const real inv_n = real(1.0) / real(f.size());
      for (int k = 0; k < n_oct - 1; ++k) {
        real acc = 0.0;
        for (const auto &fi : f) {
          const real d = fi[k] - fi[k + 1];
          acc += d * d;
        }
        R[static_cast<size_t>(k)] = std::sqrt(acc * inv_n);
      }
      return R;
    }

    /// ∇ of an MLS-smoothed scalar (shell `gradient_scalar` rod port).
    inline std::vector<vec3> gradient_scalar(asawa::rod::rod &R,
                                             const std::vector<vec3> &p_pov,
                                             const std::vector<real> &omega,
                                             real l0, real p = 3.0) {
      (void)p;
      return integrate_over_rod<vec3>(
          R, p_pov,
          [&omega](const std::vector<index_t> &edge_ids, Rod_Sum_Type &sum) {
            sum.bind(calder::scalar_datum::create(edge_ids, omega));
          },
          [l0](const index_t /*i*/, const index_t j, const vec3 &pi,
               const vec3 &pj, const std::vector<calder::datum::ptr> &data,
               Rod_Sum_Type::Node_Type node_type,
               const Rod_Sum_Type::Tree & /*tree*/) -> vec3 {
            real w = get_data<real>(node_type, j, 0, data);
            vec3 dp = pj - pi;
            vec3 dkappa = calc_d_gaussian(dp, l0);
            return -dkappa * w;
          });
    }

    /// Convolution Hess·v extras from κ(π_i−π_j) T_j / κ G_j (leaf κ only).
    ///   (Hess_κ·δdp)∗T + (δκ)∗G + (∇_i κ)(G·v)
    /// with δdp=v_j−v_i, δκ=(∂κ/∂dp)·δdp, ∇_i κ=−∂κ/∂dp.
    /// Caller adds κ∗(H_T v) via smoothed_gradient.
    inline std::vector<vec3>
    mollifier_hessian_vec(asawa::rod::rod &R, const std::vector<vec3> &p_pov,
                          const std::vector<real> &T,
                          const std::vector<vec3> &G,
                          const std::vector<vec3> &v_pov, real l0) {
      return integrate_over_rod<vec3>(
          R, p_pov,
          [&T, &G, &v_pov](const std::vector<index_t> &edge_ids,
                           Rod_Sum_Type &sum) {
            sum.bind(calder::scalar_datum::create(edge_ids, T));
            sum.bind(calder::vec3_datum::create(edge_ids, G));
            sum.bind(calder::vec3_datum::create(edge_ids, v_pov));
          },
          [l0, &v_pov](const index_t i, const index_t j, const vec3 &pi,
                       const vec3 &pj,
                       const std::vector<calder::datum::ptr> &data,
                       Rod_Sum_Type::Node_Type node_type,
                       const Rod_Sum_Type::Tree & /*tree*/) -> vec3 {
            const real Tj = get_data<real>(node_type, j, 0, data);
            const vec3 Gj = get_data<vec3>(node_type, j, 1, data);
            const vec3 vj = get_data<vec3>(node_type, j, 2, data);
            const vec3 vi = v_pov[i];
            const vec3 dp = pj - pi;
            const vec3 ddp = vj - vi;
            const vec3 dkappa = calc_d_gaussian(dp, l0);
            const mat3 Hkappa = calc_d2_gaussian(dp, l0);
            // ∇_i κ = −∂κ/∂dp;  δ(∇_i κ) = −H_κ δdp
            const vec3 d_grad_kappa = -(Hkappa * ddp);
            const real d_kappa = dkappa.dot(ddp);
            const real dT = Gj.dot(vj); // δT_j ≈ G_j · v_j
            return d_grad_kappa * Tj + (-dkappa) * dT + d_kappa * Gj;
          });
    }

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
                                   real p = 4.0, real q = 1.0)
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
          [l0, p, q](const index_t i, const index_t j,
                  const vec3 &pi, const vec3 &pj,
                  const std::vector<calder::datum::ptr> &data,
                  Rod_Sum_Type::Node_Type node_type,
                  const Rod_Sum_Type::Tree &tree) -> vec3
          {
            real w = get_data<real>(node_type, j, 0, data);
            vec3 T = get_data<vec3>(node_type, j, 1, data);
            vec3 dp = pj - pi;
            vec3 Xt = dp.cross(T);
            const real geom = dp.squaredNorm() * T.squaredNorm();
            if (geom < 1e-24) {
              return vec3::Zero();
            }
            const real sin2 = Xt.squaredNorm() / geom;
            const real angular = q == 0.0 ? 1.0 : std::pow(sin2, q/2.0);
            real kappa = calc_inv_dist(dp, l0, p);
            return -angular * w * kappa * Xt;
          });
      return us;
    }

    std::vector<mat3> covariance(asawa::rod::rod &R, const std::vector<vec3> &Nr,
                                 const std::vector<vec3> &p_pov, real r, real l0,
                                 real p = 4.0, bool normalize = true)
    {

      std::vector<vec3> &x = R.x();
      std::vector<quat> &q = R.u();

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
      const std::vector<vec3> &x = R.x();
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
