#ifndef GAUDI_HEPWORTH_LOCAL_SLAB_FIELD_FIT_HPP
#define GAUDI_HEPWORTH_LOCAL_SLAB_FIELD_FIT_HPP

#include "implicit_field.hpp"
#include "slab_pierce_stencil.hpp"

#include "gaudi/albers/darboux_cyclide.hpp"
#include "gaudi/albers/quadric.hpp"
#include "gaudi/smallvec.hpp"

#include <algorithm>
#include <array>
#include <functional>

namespace gaudi {
namespace hepworth {
namespace block {

/// Max RBF centers (subsampled ring spokes); stack matrices sized from these.
inline constexpr int slab_field_max_rbf_centers = 8;
inline constexpr int slab_field_max_ring_samples = 32;
inline constexpr int slab_field_rbf_poly_dim = 4; // linear a0 + a1 x + a2 y + a3 z

enum class slab_field_family {
  darboux = 0,
  quadric = 1,
  /// Local Gaussian RBF + linear polynomial (two-ring hub normals only).
  rbf_mls = 2,
};

struct local_rbf_mls_field {
  vec3 foot_normal = vec3::UnitZ();
  real sigma = 1.0;
  vec4 poly = vec4::Zero(); // a0 + a1 x + a2 y + a3 z
  smallvec<vec3, slab_field_max_rbf_centers> centers;
  smallvec<real, slab_field_max_rbf_centers> weights;

  field_eval eval(const vec3 &x_local) const;
};

struct slab_field_coeffs {
  slab_field_family family = slab_field_family::darboux;
  albers::vec14 darboux = albers::vec14::Zero();
  albers::vec10 quadric = albers::vec10::Zero();
  local_rbf_mls_field rbf;
};

struct local_field_fit_config {
  real sample_weight = 1.0;
  /// Row weight for ∇f(0)·n_top = 1 (fixes homogeneous scale).
  real foot_normal_weight = 1.0;
  real rbf_sigma_scale = 2.0;
  int max_rbf_centers = slab_field_max_rbf_centers;
  /// Tikhonov on coefficients when the MLS system is rank-deficient.
  real rbf_reg = 1e-10;
};

/// Cotan-weighted one-ring hub normal from spoke positions only (no full mesh).
inline vec3 ring_hub_normal(
    const vec3 &hub, const smallvec<index_t, 16> &ring_ids,
    const smallvec<real, 16> &ring_weights,
    const std::function<vec3(index_t)> &pos_at) {
  vec3 n = vec3::Zero();
  const size_t nsp = ring_weights.size();
  if (nsp == 0)
    return vec3::UnitZ();
  for (size_t k = 0; k < nsp; ++k) {
    const vec3 a = pos_at(ring_ids[k + 1]) - hub;
    const vec3 b = pos_at(ring_ids[(k + 1) % nsp + 1]) - hub;
    n += ring_weights[k] * a.cross(b);
  }
  if (n.squaredNorm() < 1e-24)
    return vec3::UnitZ();
  return n.normalized();
}

struct slab_stencil_normals {
  vec3 n_top = vec3::UnitZ();
  vec3 n_bot = vec3::UnitZ();
};

inline slab_stencil_normals
compute_stencil_ring_normals(const slab_pierce_stencil &stencil,
                             const std::function<vec3(index_t)> &pos_at) {
  slab_stencil_normals out;
  const vec3 top = pos_at(stencil.foot);
  const vec3 bot = pos_at(stencil.pierce_foot);
  out.n_top =
      ring_hub_normal(top, stencil.top_ids, stencil.top_weights, pos_at);
  out.n_bot =
      ring_hub_normal(bot, stencil.bot_ids, stencil.bot_weights, pos_at);
  return out;
}

/// Normal for MLS samples: hub normals on the two rings only (no per-spoke synthesis).
inline vec3 stencil_sample_normal(const slab_pierce_stencil &stencil,
                                  const slab_stencil_normals &ns, index_t vi) {
  if (vi == stencil.foot)
    return ns.n_top;
  if (vi == stencil.pierce_foot)
    return ns.n_bot;
  for (index_t id : stencil.top_ids)
    if (id == vi)
      return ns.n_top;
  return ns.n_bot;
}

inline field_eval eval_slab_field(const slab_field_coeffs &c, const vec3 &x) {
  switch (c.family) {
  case slab_field_family::darboux:
    return eval_darboux_field(c.darboux, x);
  case slab_field_family::quadric:
    return eval_quadric_field(c.quadric, x);
  case slab_field_family::rbf_mls:
    return c.rbf.eval(x);
  }
  return field_eval{};
}

/// Normal-aligned Darboux MLS on stencil verts (foot-local coords, origin at foot).
inline albers::vec14
fit_darboux_on_stencil(const slab_pierce_stencil &stencil,
                       const std::function<vec3(index_t)> &pos_at,
                       const local_field_fit_config &cfg) {
  const vec3 foot = pos_at(stencil.foot);
  const slab_stencil_normals ns =
      compute_stencil_ring_normals(stencil, pos_at);
  albers::normal_constrained_darboux_cyclide acc;
  real w_sum = 0.0;

  const auto add = [&](index_t vi, real w) {
    if (w <= 0.0)
      return;
    const vec3 dp = pos_at(vi) - foot;
    acc.accumulate(w * cfg.sample_weight, dp,
                   stencil_sample_normal(stencil, ns, vi));
    w_sum += w;
  };

  for (size_t k = 0; k < stencil.top_weights.size(); ++k)
    add(stencil.top_ids[k + 1], stencil.top_weights[k]);
  for (size_t k = 0; k < stencil.bot_weights.size(); ++k)
    add(stencil.bot_ids[k + 1], stencil.bot_weights[k]);

  if (cfg.foot_normal_weight > 0.0 && w_sum > 0.0) {
    const real w_foot = cfg.foot_normal_weight * w_sum;
    acc.accumulate(w_foot, vec3::Zero(), ns.n_top);
  }
  albers::vec14 Q = acc.solve();
  if (!Q.allFinite())
    return albers::vec14::Zero();
  return Q;
}

inline albers::vec10
fit_quadric_on_stencil(const slab_pierce_stencil &stencil,
                       const std::function<vec3(index_t)> &pos_at,
                       const local_field_fit_config &cfg) {
  const vec3 foot = pos_at(stencil.foot);
  const slab_stencil_normals ns =
      compute_stencil_ring_normals(stencil, pos_at);
  albers::quadric acc;
  real w_sum = 0.0;

  const auto add = [&](index_t vi, real w) {
    if (w <= 0.0)
      return;
    const vec3 dp = pos_at(vi) - foot;
    acc.accumulate(w * cfg.sample_weight, dp,
                   stencil_sample_normal(stencil, ns, vi));
    w_sum += w;
  };

  for (size_t k = 0; k < stencil.top_weights.size(); ++k)
    add(stencil.top_ids[k + 1], stencil.top_weights[k]);
  for (size_t k = 0; k < stencil.bot_weights.size(); ++k)
    add(stencil.bot_ids[k + 1], stencil.bot_weights[k]);

  if (cfg.foot_normal_weight > 0.0 && w_sum > 0.0) {
    const real w_foot = cfg.foot_normal_weight * w_sum;
    acc.accumulate(w_foot, vec3::Zero(), ns.n_top);
  }
  return acc.solve();
}

/// Foot-local LS quadric with f=0 on spoke ring samples + ∇f(0)·n_top = 1.
/// Used by incompressible φ so gauge energy vanishes at rest (Q_r matches stencil verts).
inline albers::vec10
fit_quadric_level_zero_on_stencil(
    const slab_pierce_stencil &stencil,
    const std::function<vec3(index_t)> &pos_at,
    const local_field_fit_config &cfg) {
  const vec3 foot = pos_at(stencil.foot);
  const slab_stencil_normals ns =
      compute_stencil_ring_normals(stencil, pos_at);

  Eigen::Matrix<real, Eigen::Dynamic, 10> A;
  Eigen::Matrix<real, Eigen::Dynamic, 1> b;
  int n_rows = 0;
  real w_sum = 0.0;

  const auto count_spoke = [&](index_t vi, real w) {
    if (w <= 0.0 || vi == stencil.foot || vi == stencil.pierce_foot)
      return;
    ++n_rows;
    w_sum += w;
  };
  for (size_t k = 0; k < stencil.top_weights.size(); ++k)
    count_spoke(stencil.top_ids[k + 1], stencil.top_weights[k]);
  for (size_t k = 0; k < stencil.bot_weights.size(); ++k)
    count_spoke(stencil.bot_ids[k + 1], stencil.bot_weights[k]);

  if (cfg.foot_normal_weight > 0.0 && w_sum > 0.0)
    ++n_rows;

  if (n_rows <= 0)
    return albers::vec10::Zero();

  A.resize(n_rows, 10);
  b.setZero(n_rows);
  int row = 0;

  const auto add_f_zero = [&](index_t vi, real w) {
    if (w <= 0.0 || vi == stencil.foot || vi == stencil.pierce_foot)
      return;
    const vec3 dp = pos_at(vi) - foot;
    const real sw = std::sqrt(w * cfg.sample_weight);
    A.row(row) = sw * albers::mk_quad_A(dp).row(0);
    ++row;
  };

  for (size_t k = 0; k < stencil.top_weights.size(); ++k)
    add_f_zero(stencil.top_ids[k + 1], stencil.top_weights[k]);
  for (size_t k = 0; k < stencil.bot_weights.size(); ++k)
    add_f_zero(stencil.bot_ids[k + 1], stencil.bot_weights[k]);

  if (cfg.foot_normal_weight > 0.0 && w_sum > 0.0) {
    const real sw = std::sqrt(cfg.foot_normal_weight * w_sum);
    albers::vec10 grad_row = albers::vec10::Zero();
    grad_row[6] = 2.0 * ns.n_top.x();
    grad_row[7] = 2.0 * ns.n_top.y();
    grad_row[8] = 2.0 * ns.n_top.z();
    A.row(row) = sw * grad_row;
    b(row) = sw;
  }

  return A.colPivHouseholderQr().solve(b);
}

inline slab_field_coeffs fit_phi_quadric_on_stencil(
    const slab_pierce_stencil &stencil,
    const std::function<vec3(index_t)> &pos_at,
    const local_field_fit_config &cfg) {
  slab_field_coeffs out;
  out.family = slab_field_family::quadric;
  out.quadric = fit_quadric_level_zero_on_stencil(stencil, pos_at, cfg);
  return out;
}

namespace detail {

inline constexpr int rbf_max_rows =
    slab_field_max_ring_samples + 2; // samples + two hub normal rows
inline constexpr int rbf_max_cols =
    slab_field_max_rbf_centers + slab_field_rbf_poly_dim;

using rbf_fit_mat = Eigen::Matrix<real, rbf_max_rows, rbf_max_cols>;
using rbf_fit_vec = Eigen::Matrix<real, rbf_max_rows, 1>;
using rbf_sol_vec = Eigen::Matrix<real, rbf_max_cols, 1>;

struct rbf_ring_sample {
  vec3 local = vec3::Zero();
  real weight = 0.0;
};

inline real rbf_phi(const vec3 &a, const vec3 &b, real sigma) {
  const real r2 = (a - b).squaredNorm();
  return std::exp(-r2 / (sigma * sigma));
}

inline vec3 rbf_grad_phi(const vec3 &x, const vec3 &c, real sigma) {
  const real p = rbf_phi(x, c, sigma);
  return (-2.0 / (sigma * sigma)) * (x - c) * p;
}

inline void add_rbf_normal_row(rbf_fit_mat &A, rbf_fit_vec &b, int row, const vec3 &x_eval,
                               const vec3 &n, real row_weight, real target,
                               const smallvec<vec3, slab_field_max_rbf_centers> &centers,
                               real sigma, int n_rbf) {
  for (int k = 0; k < n_rbf; ++k) {
    const vec3 gphi = rbf_grad_phi(x_eval, centers[static_cast<size_t>(k)], sigma);
    A(row, k) = row_weight * n.dot(gphi);
  }
  A(row, n_rbf + 0) = 0.0;
  A(row, n_rbf + 1) = row_weight * n.x();
  A(row, n_rbf + 2) = row_weight * n.y();
  A(row, n_rbf + 3) = row_weight * n.z();
  b(row) = row_weight * target;
}

inline rbf_sol_vec solve_rbf_least_squares(const rbf_fit_mat &A, const rbf_fit_vec &b,
                                           int n_rows, int n_cols, real reg) {
  const Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> block =
      A.topLeftCorner(n_rows, n_cols);
  const Eigen::Matrix<real, Eigen::Dynamic, 1> rhs = b.head(n_rows);
  Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> ata =
      block.transpose() * block;
  if (reg > 0.0)
    ata.diagonal().array() += reg;
  const Eigen::Matrix<real, Eigen::Dynamic, 1> atb = block.transpose() * rhs;
  Eigen::Matrix<real, Eigen::Dynamic, 1> sol_dyn =
      ata.ldlt().solve(atb);
  rbf_sol_vec sol = rbf_sol_vec::Zero();
  sol.head(n_cols) = sol_dyn;
  return sol;
}

} // namespace detail

inline bool rbf_fields_near_equal(const local_rbf_mls_field &a,
                                  const local_rbf_mls_field &b, real tol) {
  if (a.centers.size() != b.centers.size())
    return false;
  if ((a.foot_normal - b.foot_normal).norm() > tol)
    return false;
  if (std::abs(a.sigma - b.sigma) > tol * std::max(a.sigma, real(1.0)))
    return false;
  for (size_t k = 0; k < a.centers.size(); ++k) {
    if ((a.centers[k] - b.centers[k]).norm() > tol)
      return false;
    if (std::abs(a.weights[k] - b.weights[k]) > tol)
      return false;
  }
  return (a.poly - b.poly).norm() <= tol;
}

inline local_rbf_mls_field
fit_rbf_mls_on_stencil(const slab_pierce_stencil &stencil,
                       const std::function<vec3(index_t)> &pos_at,
                       const local_field_fit_config &cfg, real avg_edge_length,
                       const smallvec<vec3, slab_field_max_rbf_centers> *fixed_centers =
                           nullptr,
                       real fixed_sigma = 0.0) {
  local_rbf_mls_field out;
  const slab_stencil_normals ns =
      compute_stencil_ring_normals(stencil, pos_at);
  out.foot_normal = ns.n_top;

  const vec3 foot = pos_at(stencil.foot);
  const vec3 pierce_local = pos_at(stencil.pierce_foot) - foot;

  std::array<detail::rbf_ring_sample, slab_field_max_ring_samples> samples{};
  int n_samples = 0;

  const auto add_spoke = [&](index_t vi, real w) {
    if (w <= 0.0 || vi == stencil.foot || vi == stencil.pierce_foot)
      return;
    if (n_samples >= slab_field_max_ring_samples)
      return;
    samples[static_cast<size_t>(n_samples)] = {pos_at(vi) - foot,
                                               w * cfg.sample_weight};
    ++n_samples;
  };

  for (size_t k = 0; k < stencil.top_weights.size(); ++k)
    add_spoke(stencil.top_ids[k + 1], stencil.top_weights[k]);
  for (size_t k = 0; k < stencil.bot_weights.size(); ++k)
    add_spoke(stencil.bot_ids[k + 1], stencil.bot_weights[k]);

  const int max_centers =
      std::clamp(cfg.max_rbf_centers, 0, slab_field_max_rbf_centers);
  int n_rbf = 0;
  if (fixed_centers && !fixed_centers->empty()) {
    out.centers = *fixed_centers;
    n_rbf = static_cast<int>(out.centers.size());
    out.sigma = fixed_sigma > 0.0
                    ? fixed_sigma
                    : std::max(cfg.rbf_sigma_scale * avg_edge_length, real(1e-6));
  } else {
    n_rbf = std::max(
        0, std::min({n_samples, max_centers, n_samples + 2 - slab_field_rbf_poly_dim}));
    out.sigma = std::max(cfg.rbf_sigma_scale * avg_edge_length, real(1e-6));
    out.centers.clear();
    for (int k = 0; k < n_rbf; ++k)
      out.centers.push_back(samples[static_cast<size_t>(k)].local);
  }

  const int n_poly = slab_field_rbf_poly_dim;
  const int n_cols = n_rbf + n_poly;
  const int n_rows = n_samples + 2;

  detail::rbf_fit_mat A = detail::rbf_fit_mat::Zero();
  detail::rbf_fit_vec b = detail::rbf_fit_vec::Zero();

  for (int i = 0; i < n_samples; ++i) {
    const vec3 xi = samples[static_cast<size_t>(i)].local;
    const real wi = samples[static_cast<size_t>(i)].weight;
    for (int k = 0; k < n_rbf; ++k)
      A(i, k) = wi * detail::rbf_phi(xi, out.centers[static_cast<size_t>(k)],
                                     out.sigma);
    A(i, n_rbf + 0) = wi * 1.0;
    A(i, n_rbf + 1) = wi * xi.x();
    A(i, n_rbf + 2) = wi * xi.y();
    A(i, n_rbf + 3) = wi * xi.z();
  }

  const real n_w =
      cfg.foot_normal_weight > 0.0 ? cfg.foot_normal_weight : real(1.0);
  // Inhomogeneous gauge: ∇f(0)·n_top = 1 (TDD rest orientation).
  detail::add_rbf_normal_row(A, b, n_samples, vec3::Zero(), ns.n_top, n_w, 1.0,
                             out.centers, out.sigma, n_rbf);
  detail::add_rbf_normal_row(A, b, n_samples + 1, pierce_local, ns.n_bot, n_w, 0.0,
                             out.centers, out.sigma, n_rbf);

  detail::rbf_sol_vec sol = detail::solve_rbf_least_squares(
      A, b, n_rows, n_cols, cfg.rbf_reg);

  out.weights.clear();
  for (int k = 0; k < n_rbf; ++k)
    out.weights.push_back(sol(k));
  out.poly = vec4(sol(n_rbf + 0), sol(n_rbf + 1), sol(n_rbf + 2), sol(n_rbf + 3));

  // Single sign fix if the normal constraint was satisfied with wrong hemisphere.
  field_eval e0 = out.eval(vec3::Zero());
  if (e0.grad.dot(out.foot_normal) < 0.0) {
    for (real &wk : out.weights)
      wk = -wk;
    out.poly = -out.poly;
  }
  return out;
}

inline field_eval local_rbf_mls_field::eval(const vec3 &x_local) const {
  field_eval out;
  out.f = poly[0] + poly[1] * x_local.x() + poly[2] * x_local.y() +
          poly[3] * x_local.z();
  out.grad = vec3(poly[1], poly[2], poly[3]);
  out.hess = mat3::Zero();
  out.has_hess = true;

  for (size_t k = 0; k < centers.size(); ++k) {
    const real w = weights[k];
    const vec3 c = centers[k];
    const real p = detail::rbf_phi(x_local, c, sigma);
    const vec3 g = detail::rbf_grad_phi(x_local, c, sigma);
    out.f += w * p;
    out.grad += w * g;
    const vec3 d = x_local - c;
    const real s2 = sigma * sigma;
    out.hess += w * ((-2.0 / s2) * p * mat3::Identity() +
                     (4.0 / (s2 * s2)) * p * d * d.transpose());
  }
  return out;
}

inline slab_field_coeffs
fit_slab_field_on_stencil(slab_field_family family,
                          const slab_pierce_stencil &stencil,
                          const std::function<vec3(index_t)> &pos_at,
                          const local_field_fit_config &cfg,
                          real avg_edge_length = 1.0) {
  slab_field_coeffs out;
  out.family = family;
  switch (family) {
  case slab_field_family::darboux:
    out.darboux = fit_darboux_on_stencil(stencil, pos_at, cfg);
    break;
  case slab_field_family::quadric:
    out.quadric = fit_quadric_on_stencil(stencil, pos_at, cfg);
    break;
  case slab_field_family::rbf_mls:
    out.rbf = fit_rbf_mls_on_stencil(stencil, pos_at, cfg, avg_edge_length);
    break;
  }
  return out;
}

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // GAUDI_HEPWORTH_LOCAL_SLAB_FIELD_FIT_HPP