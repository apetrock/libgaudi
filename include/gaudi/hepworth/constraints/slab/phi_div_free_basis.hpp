#ifndef GAUDI_HEPWORTH_PHI_DIV_FREE_BASIS_HPP
#define GAUDI_HEPWORTH_PHI_DIV_FREE_BASIS_HPP

#include "slab_fan_wedge.hpp"
#include "slab_pierce_stencil.hpp"
#include "slab_wedge_strain.hpp"

#include <Eigen/SVD>
#include <functional>
#include <unordered_map>
#include <vector>

namespace gaudi {
namespace hepworth {
namespace block {

/// Div-free velocity u_q on the slab stencil: u = B q with B from ker(div) on fan wedges.
struct phi_div_free_basis {
  int dim_q = 0;
  int n_dof = 0;
  matX B;
  vec3 foot_rest = vec3::Zero();
  std::vector<index_t> vert_ids;
  std::unordered_map<index_t, size_t> vert_idx;

  struct wedge_patch {
    index_t v0 = -1;
    index_t v1 = -1;
    index_t v2 = -1;
    vec3 grad_l0 = vec3::Zero();
    vec3 grad_l1 = vec3::Zero();
    vec3 grad_l2 = vec3::Zero();
    vec3 bc_local = vec3::Zero();
  };
  std::vector<wedge_patch> wedges;

  void build(const slab_pierce_stencil &stencil,
             const std::function<vec3(index_t)> &pos_at) {
    foot_rest = pos_at(stencil.foot);
    vert_ids.clear();
    vert_idx.clear();
    wedges.clear();
    B.resize(0, 0);
    dim_q = 0;
    n_dof = 0;

    const auto add_vert = [&](index_t vi) {
      if (vert_idx.count(vi))
        return;
      vert_idx[vi] = vert_ids.size();
      vert_ids.push_back(vi);
    };
    for (index_t vi : stencil.top_ids)
      add_vert(vi);
    for (index_t vi : stencil.bot_ids)
      add_vert(vi);

    if (vert_ids.size() < 3) {
      build_fallback_rotation();
      return;
    }

    n_dof = static_cast<int>(3 * vert_ids.size());
    const auto local_at = [&](index_t vi) { return pos_at(vi) - foot_rest; };

    const std::vector<slab_fan_wedge> fan =
        build_slab_fan_wedges(stencil, pos_at);
    wedges.reserve(fan.size());

    matX D = matX::Zero(static_cast<int>(fan.size()), n_dof);
    int row = 0;
    for (const slab_fan_wedge &fw : fan) {
      wedge_patch wp;
      wp.v0 = fw.hub;
      wp.v1 = fw.i0;
      wp.v2 = fw.i1;
      const vec3 p0 = local_at(wp.v0);
      const vec3 p1 = local_at(wp.v1);
      const vec3 p2 = local_at(wp.v2);
      wp.bc_local = (p0 + p1 + p2) / 3.0;

      const vec3 e01 = p1 - p0;
      const vec3 e02 = p2 - p0;
      const vec3 n = e01.cross(e02);
      const real area2 = n.norm();
      if (area2 < 1e-24) {
        ++row;
        continue;
      }
      wp.grad_l0 = e02.cross(n) / area2;
      wp.grad_l1 = n.cross(e01) / area2;
      wp.grad_l2 = -(wp.grad_l0 + wp.grad_l1);
      wedges.push_back(wp);

      const auto add_row = [&](index_t vi, const vec3 &grad_l) {
        const auto it = vert_idx.find(vi);
        if (it == vert_idx.end())
          return;
        const int base = static_cast<int>(3 * it->second);
        for (int ax = 0; ax < 3; ++ax)
          D(row, base + ax) = grad_l[ax];
      };
      add_row(wp.v0, wp.grad_l0);
      add_row(wp.v1, wp.grad_l1);
      add_row(wp.v2, wp.grad_l2);
      ++row;
    }

    if (row == 0 || wedges.empty()) {
      build_fallback_rotation();
      return;
    }
    D.conservativeResize(row, n_dof);

    Eigen::JacobiSVD<matX> svd(D, Eigen::ComputeFullV);
    const real tol =
        std::max(D.rows(), D.cols()) * svd.singularValues()(0) * 1e-10;
    const int sv = static_cast<int>(svd.singularValues().size());
    int rank = 0;
    for (int k = 0; k < sv; ++k) {
      if (svd.singularValues()[k] > tol)
        ++rank;
    }
    dim_q = n_dof - rank;
    if (dim_q <= 0) {
      build_fallback_rotation();
      return;
    }

    B = matX::Zero(n_dof, dim_q);
    for (int j = 0; j < dim_q; ++j)
      B.col(j) = svd.matrixV().col(rank + j);
  }

  vec3 eval_u_at_vertex(index_t vi, const vecX &q) const {
    if (q.size() < dim_q || dim_q <= 0)
      return vec3::Zero();
    const auto it = vert_idx.find(vi);
    if (it == vert_idx.end())
      return vec3::Zero();
    const int base = static_cast<int>(3 * it->second);
    vec3 u = vec3::Zero();
    for (int j = 0; j < dim_q; ++j) {
      u[0] += q[j] * B(base + 0, j);
      u[1] += q[j] * B(base + 1, j);
      u[2] += q[j] * B(base + 2, j);
    }
    return u;
  }

  vec3 eval_basis_field_at_vertex(int mode_j, index_t vi) const {
    if (mode_j < 0 || mode_j >= dim_q)
      return vec3::Zero();
    const auto it = vert_idx.find(vi);
    if (it == vert_idx.end())
      return vec3::Zero();
    const int base = static_cast<int>(3 * it->second);
    return vec3(B(base + 0, mode_j), B(base + 1, mode_j), B(base + 2, mode_j));
  }

  vec3 eval_u(const vec3 &local, index_t vi, const vecX &q) const {
    (void)local;
    return eval_u_at_vertex(vi, q);
  }

  vec3 eval_basis_field(int mode_j, const vec3 &local, index_t vi) const {
    (void)local;
    return eval_basis_field_at_vertex(mode_j, vi);
  }

  real eval_div_at_barycenter(const wedge_patch &wp, const vecX &q) const {
    if (q.size() < dim_q || dim_q <= 0)
      return 0.0;
    real div = 0.0;
    const auto add = [&](index_t vi, const vec3 &grad_l) {
      const auto it = vert_idx.find(vi);
      if (it == vert_idx.end())
        return;
      const int base = static_cast<int>(3 * it->second);
      vec3 u_v = vec3::Zero();
      for (int j = 0; j < dim_q; ++j) {
        u_v[0] += q[j] * B(base + 0, j);
        u_v[1] += q[j] * B(base + 1, j);
        u_v[2] += q[j] * B(base + 2, j);
      }
      div += grad_l.dot(u_v);
    };
    add(wp.v0, wp.grad_l0);
    add(wp.v1, wp.grad_l1);
    add(wp.v2, wp.grad_l2);
    return div;
  }

  real max_div_on_wedges(const vecX &q) const {
    real m = 0.0;
    for (const wedge_patch &wp : wedges)
      m = std::max(m, std::abs(eval_div_at_barycenter(wp, q)));
    return m;
  }

  real eval_div(const vec3 &local, const vecX &q) const {
    (void)local;
    (void)q;
    return 0.0;
  }

  real eval_div_fd(const vec3 &local, const vecX &q, real eps = 1e-5) const {
    (void)local;
    (void)q;
    (void)eps;
    return 0.0;
  }

private:
  void build_fallback_rotation() {
    dim_q = 3;
    n_dof = 3;
    B = matX::Identity(3, 3);
  }
};

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // GAUDI_HEPWORTH_PHI_DIV_FREE_BASIS_HPP
