#ifndef GAUDI_HEPWORTH_RBF_DISPLACEMENT_BASIS_HPP
#define GAUDI_HEPWORTH_RBF_DISPLACEMENT_BASIS_HPP

#include "local_slab_field_fit.hpp"
#include "slab_pierce_stencil.hpp"

#include <cmath>
#include <functional>
#include <unordered_map>
#include <vector>

namespace gaudi {
namespace hepworth {
namespace block {

/// Fixed-center RBF displacement ROM on a slab stencil.
///
/// Modes (chosen once at rest): for each center c and axis a,
///   ψ_{c,a}(x) = φ(|x − c| / σ) ê_a
/// with Gaussian φ. Live DOFs are coeffs q; u = Ψ q.
struct rbf_displacement_basis {
  int dim_q = 0;
  int n_centers = 0;
  real sigma = 1.0;
  vec3 foot_rest = vec3::Zero();
  std::vector<vec3> centers; // foot-local
  std::vector<index_t> vert_ids;
  std::unordered_map<index_t, size_t> vert_idx;
  /// Dense mode matrix: rows = 3 * vert_ids, cols = dim_q (= 3 * n_centers).
  matX B;

  void clear() {
    dim_q = 0;
    n_centers = 0;
    sigma = 1.0;
    foot_rest = vec3::Zero();
    centers.clear();
    vert_ids.clear();
    vert_idx.clear();
    B.resize(0, 0);
  }

  void build(const slab_pierce_stencil &stencil,
             const std::function<vec3(index_t)> &pos_at, real avg_edge_length,
             int max_centers = 4, real sigma_scale = 2.0) {
    clear();
    foot_rest = pos_at(stencil.foot);
    sigma = std::max(sigma_scale * avg_edge_length, real(1e-6));

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
      build_fallback_translation();
      return;
    }

    const auto local_at = [&](index_t vi) { return pos_at(vi) - foot_rest; };
    const vec3 pierce_local = local_at(stencil.pierce_foot);

    centers.clear();
    centers.push_back(vec3::Zero());
    if (pierce_local.squaredNorm() > 1e-24)
      centers.push_back(pierce_local);

    std::vector<vec3> spoke_locals;
    spoke_locals.reserve(stencil.top_ids.size() + stencil.bot_ids.size());
    for (size_t k = 1; k < stencil.top_ids.size(); ++k)
      spoke_locals.push_back(local_at(stencil.top_ids[k]));
    for (size_t k = 1; k < stencil.bot_ids.size(); ++k)
      spoke_locals.push_back(local_at(stencil.bot_ids[k]));

    const int budget = std::clamp(max_centers, 1, slab_field_max_rbf_centers);
    if (!spoke_locals.empty() && static_cast<int>(centers.size()) < budget) {
      const int need = budget - static_cast<int>(centers.size());
      const int n = static_cast<int>(spoke_locals.size());
      for (int i = 0; i < need; ++i) {
        const int idx = (i * std::max(n, 1)) / std::max(need, 1);
        centers.push_back(spoke_locals[static_cast<size_t>(idx % n)]);
      }
    }

    n_centers = static_cast<int>(centers.size());
    if (n_centers <= 0) {
      build_fallback_translation();
      return;
    }

    dim_q = 3 * n_centers;
    const int n_dof = static_cast<int>(3 * vert_ids.size());
    B = matX::Zero(n_dof, dim_q);

    for (size_t vi = 0; vi < vert_ids.size(); ++vi) {
      const vec3 x = local_at(vert_ids[vi]);
      const int row = static_cast<int>(3 * vi);
      for (int c = 0; c < n_centers; ++c) {
        const real phi =
            detail::rbf_phi(x, centers[static_cast<size_t>(c)], sigma);
        B(row + 0, 3 * c + 0) = phi;
        B(row + 1, 3 * c + 1) = phi;
        B(row + 2, 3 * c + 2) = phi;
      }
    }
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

  /// API parity with phi_div_free_basis (local unused — modes tabulated on verts).
  vec3 eval_u(const vec3 &local, index_t vi, const vecX &q) const {
    (void)local;
    return eval_u_at_vertex(vi, q);
  }

  vec3 eval_basis_field(int mode_j, const vec3 &local, index_t vi) const {
    (void)local;
    return eval_basis_field_at_vertex(mode_j, vi);
  }

private:
  void build_fallback_translation() {
    n_centers = 0;
    centers.clear();
    dim_q = 3;
    const int n_dof = std::max(3, static_cast<int>(3 * vert_ids.size()));
    B = matX::Zero(n_dof, 3);
    const int n_v = std::max(1, static_cast<int>(vert_ids.size()));
    for (int vi = 0; vi < n_v; ++vi) {
      B(3 * vi + 0, 0) = 1.0;
      B(3 * vi + 1, 1) = 1.0;
      B(3 * vi + 2, 2) = 1.0;
    }
    if (vert_ids.empty()) {
      vert_ids.push_back(-1);
      vert_idx[-1] = 0;
    }
  }
};

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // GAUDI_HEPWORTH_RBF_DISPLACEMENT_BASIS_HPP
