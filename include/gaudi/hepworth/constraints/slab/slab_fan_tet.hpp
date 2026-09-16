#ifndef GAUDI_HEPWORTH_SLAB_FAN_TET_HPP
#define GAUDI_HEPWORTH_SLAB_FAN_TET_HPP

#include "phi_div_free_basis.hpp"
#include "slab_pierce_stencil.hpp"
#include "slab_wedge_strain.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <functional>
#include <unordered_map>
#include <vector>

namespace gaudi {
namespace hepworth {
namespace block {

/// One fan tet: anchor + three hub-relative edges, volume strain at rest baked in.
struct slab_fan_tet {
  index_t anchor = -1;
  index_t v1 = -1;
  index_t v2 = -1;
  index_t v3 = -1;
  real weight = 1.0;
  mat3 E_rest = mat3::Identity();
  mat3 E_rest_inv = mat3::Identity();
  real V_rest = 0.0;
  index_t row0 = 0;
  bool bottom_cap = false;
};

struct slab_fan_gauge_cache {
  phi_div_free_basis basis;
  vecX q_rest;
  real align_rms = 0.0;
};

inline real stencil_vert_weight(const slab_pierce_stencil &stencil, index_t vi) {
  if (vi == stencil.foot)
    return 0.0;
  for (size_t k = 0; k < stencil.top_weights.size(); ++k) {
    if (stencil.top_ids[k + 1] == vi)
      return stencil.top_weights[k];
  }
  for (size_t k = 0; k < stencil.bot_weights.size(); ++k) {
    if (stencil.bot_ids[k + 1] == vi)
      return stencil.bot_weights[k];
  }
  if (vi == stencil.pierce_foot) {
    real w = 0.0;
    for (real wb : stencil.bot_weights)
      w += wb;
    return w;
  }
  return 1.0;
}

/// Weighted LS: q_rest = argmin ||W^{1/2}(x_local - B q)||.
inline vecX fit_slab_fan_nullspace_q(const phi_div_free_basis &basis,
                                     const slab_pierce_stencil &stencil,
                                     const std::function<vec3(index_t)> &pos_at,
                                     real &align_rms_out) {
  align_rms_out = 0.0;
  if (basis.dim_q <= 0 || basis.B.size() == 0)
    return vecX::Zero(0);

  vecX x_local = vecX::Zero(basis.n_dof);
  vecX w_dof = vecX::Zero(basis.n_dof);
  for (size_t k = 0; k < basis.vert_ids.size(); ++k) {
    const index_t vi = basis.vert_ids[k];
    const vec3 xl = pos_at(vi) - basis.foot_rest;
    const real wv = stencil_vert_weight(stencil, vi);
    const int base = static_cast<int>(3 * k);
    for (int ax = 0; ax < 3; ++ax) {
      x_local[base + ax] = xl[ax];
      w_dof[base + ax] = std::max(wv, real(0.0));
    }
  }

  const vecX w_sqrt = w_dof.cwiseSqrt();
  const matX Bw = w_sqrt.asDiagonal() * basis.B;
  const vecX xw = w_sqrt.cwiseProduct(x_local);
  const matX normal = Bw.transpose() * Bw;
  const vecX rhs = Bw.transpose() * xw;
  vecX q_rest = vecX::Zero(basis.dim_q);
  if (normal.rows() > 0) {
    Eigen::LDLT<matX> ldlt(normal);
    if (ldlt.info() == Eigen::Success)
      q_rest = ldlt.solve(rhs);
  }

  const vecX resid = x_local - basis.B * q_rest;
  const real num = (w_sqrt.cwiseProduct(resid)).norm();
  const real den = std::max(w_sqrt.cwiseProduct(x_local).norm(), real(1e-12));
  align_rms_out = num / den;
  return q_rest;
}

inline void bake_fan_tet_rest(slab_fan_tet &t,
                              const std::function<vec3(index_t)> &pos_at,
                              gaudi::real w_scale) {
  const vec3 anchor = pos_at(t.anchor);
  const vec3 e0 = pos_at(t.v1) - anchor;
  const vec3 e1 = pos_at(t.v2) - anchor;
  const vec3 e2 = pos_at(t.v3) - anchor;

  mat3 E = mat3::Zero();
  E.col(0) = e0;
  E.col(1) = e1;
  E.col(2) = e2;

  gaudi::real volume_det = e0.dot(e1.cross(e2));
  if (!std::isfinite(volume_det)) {
    t.E_rest = mat3::Identity();
    t.E_rest_inv = mat3::Identity();
    t.V_rest = 0.0;
    t.weight = w_scale * t.weight;
    return;
  }

  if (volume_det < 0.0) {
    std::swap(t.v1, t.v2);
    E.col(0) = pos_at(t.v1) - anchor;
    E.col(1) = pos_at(t.v2) - anchor;
    volume_det = -volume_det;
  }

  if (volume_det < 1e-24) {
    t.E_rest = mat3::Identity();
    t.E_rest_inv = mat3::Identity();
    t.V_rest = 0.0;
    t.weight = w_scale * t.weight;
    return;
  }

  t.E_rest = E;
  t.E_rest_inv = E.inverse();
  t.V_rest = volume_det / 6.0;
  t.weight = w_scale * t.weight * std::cbrt(volume_det);
}

inline std::vector<slab_fan_tet>
build_slab_fan_tets(const slab_pierce_stencil &stencil,
                    const std::function<vec3(index_t)> &pos_at, gaudi::real w_scale = 1.0) {
  std::vector<slab_fan_tet> out;
  const size_t n_top = stencil.top_ids.size() - 1;
  const size_t n_bot = stencil.bot_ids.size() - 1;
  out.reserve(n_top + n_bot);

  for (size_t k = 0; k < n_top; ++k) {
    slab_fan_tet t;
    t.anchor = stencil.foot;
    t.v1 = stencil.top_ids[k + 1];
    t.v2 = stencil.top_ids[(k + 1) % n_top + 1];
    t.v3 = stencil.pierce_foot;
    t.weight = stencil.top_weights[k];
    t.bottom_cap = false;
    bake_fan_tet_rest(t, pos_at, w_scale);
    out.push_back(t);
  }

  for (size_t k = 0; k < n_bot; ++k) {
    slab_fan_tet t;
    t.anchor = stencil.pierce_foot;
    t.v1 = stencil.bot_ids[k + 1];
    t.v2 = stencil.bot_ids[(k + 1) % n_bot + 1];
    t.v3 = stencil.foot;
    t.weight = stencil.bot_weights[k];
    t.bottom_cap = true;
    bake_fan_tet_rest(t, pos_at, w_scale);
    out.push_back(t);
  }
  return out;
}

inline mat3 clamp_tet_F_volume_preserving(const mat3 &F) {
  if (F.squaredNorm() < 1e-24)
    return mat3::Identity();
  Eigen::JacobiSVD<mat3> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
  vec3 S = svd.singularValues();
  const real det_sigma = S(0) * S(1) * S(2);
  if (std::abs(det_sigma) < 1e-24)
    return mat3::Identity();
  const real lam = std::cbrt(real(1.0) / det_sigma);
  S *= lam;
  return svd.matrixU() * S.asDiagonal() * svd.matrixV().transpose();
}

/// Local 3×3 volume projection for one fan tet (live edges in, weighted target cols out).
inline bool project_fan_tet_volume(const vec3 &e0_cur, const vec3 &e1_cur,
                                   const vec3 &e2_cur, const mat3 &E_rest,
                                   const mat3 &E_rest_inv, real w, vec3 &out0,
                                   vec3 &out1, vec3 &out2) {
  mat3 E_cur;
  E_cur.col(0) = e0_cur;
  E_cur.col(1) = e1_cur;
  E_cur.col(2) = e2_cur;

  if (std::abs(E_rest.col(0).dot(E_rest.col(1).cross(E_rest.col(2)))) < 1e-24 ||
      std::abs(E_cur.col(0).dot(E_cur.col(1).cross(E_cur.col(2)))) < 1e-24) {
    out0 = w * E_rest.col(0);
    out1 = w * E_rest.col(1);
    out2 = w * E_rest.col(2);
    return false;
  }

  const mat3 F_vol = clamp_tet_F_volume_preserving(E_cur * E_rest_inv);
  const mat3 E_star = F_vol * E_rest;
  out0 = w * E_star.col(0);
  out1 = w * E_star.col(1);
  out2 = w * E_star.col(2);
  return true;
}

/// Build fan tets after optional div-nullspace gauge fit (tests / φ path only — not used in w_tet hot init).
inline std::vector<slab_fan_tet>
build_slab_fan_tets_aligned(const slab_pierce_stencil &stencil,
                            const std::function<vec3(index_t)> &pos_at,
                            slab_fan_gauge_cache &gauge_out, gaudi::real w_scale = 1.0) {
  gauge_out.basis.build(stencil, pos_at);
  gauge_out.q_rest =
      fit_slab_fan_nullspace_q(gauge_out.basis, stencil, pos_at, gauge_out.align_rms);
  // Gauge records ring alignment into ker(div); rest tet geometry uses raw hub-relative edges.
  return build_slab_fan_tets(stencil, pos_at, w_scale);
}

inline real eval_fan_tet_volume(const vec3 &e0, const vec3 &e1, const vec3 &e2) {
  return e0.dot(e1.cross(e2)) / 6.0;
}

inline real eval_fan_tet_volume_at(const slab_fan_tet &t,
                                   const std::function<vec3(index_t)> &pos_at) {
  const vec3 a = pos_at(t.anchor);
  return eval_fan_tet_volume(pos_at(t.v1) - a, pos_at(t.v2) - a, pos_at(t.v3) - a);
}

inline real eval_fan_tet_volume_energy(
    const std::vector<slab_fan_tet> &tets, const slab_ring_positions &ring,
    const std::function<vec3(index_t)> &fallback_at) {
  real E = 0.0;
  const auto pos_at = [&](index_t vi) -> vec3 {
    const auto it = ring.idx.find(vi);
    if (it != ring.idx.end())
      return ring.x[it->second];
    return fallback_at(vi);
  };
  for (const slab_fan_tet &t : tets) {
    const vec3 a = pos_at(t.anchor);
    const real V = eval_fan_tet_volume(pos_at(t.v1) - a, pos_at(t.v2) - a,
                                       pos_at(t.v3) - a);
    const real r = V - t.V_rest;
    E += t.weight * r * r;
  }
  return E;
}

inline void fan_tet_grad_slab(const slab_fan_tet &t, const vec3 &e0, const vec3 &e1,
                              const vec3 &e2,
                              const std::function<void(index_t, const vec3 &)> &add) {
  const vec3 g0 = e1.cross(e2) / 6.0;
  const vec3 g1 = e2.cross(e0) / 6.0;
  const vec3 g2 = e0.cross(e1) / 6.0;
  add(t.v1, g0);
  add(t.v2, g1);
  add(t.v3, g2);
  add(t.anchor, -(g0 + g1 + g2));
}

/// Gauss–Newton on Σ w_k (V_k − V_k^rest)² over stencil ring DOFs.
inline slab_ring_positions resolve_fan_tet_volume_gauss_newton(
    const slab_pierce_stencil &stencil, const std::vector<slab_fan_tet> &tets,
    const slab_ring_positions &x0,
    const std::function<vec3(index_t)> &fallback_at, int iters,
    real lm_damp = 1e-8) {
  slab_ring_positions x = x0;
  const int n = static_cast<int>(x.x.size());
  const int nd = 3 * n;
  if (nd <= 0 || tets.empty() || iters <= 0)
    return x;
  (void)stencil;

  const real e_start = eval_fan_tet_volume_energy(tets, x, fallback_at);
  if (e_start < 1e-24)
    return x;

  for (int iter = 0; iter < iters; ++iter) {
    const auto pos_at = [&](index_t vi) -> vec3 {
      const auto it = x.idx.find(vi);
      if (it != x.idx.end())
        return x.x[it->second];
      return fallback_at(vi);
    };

    Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> JTJ =
        Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic>::Zero(nd, nd);
    Eigen::Matrix<real, Eigen::Dynamic, 1> JTr =
        Eigen::Matrix<real, Eigen::Dynamic, 1>::Zero(nd);

    const auto add_grad = [&](index_t vi, const vec3 &g,
                              std::vector<vec3> &slot) {
      const auto it = x.idx.find(vi);
      if (it != x.idx.end())
        slot[it->second] += g;
    };

    for (const slab_fan_tet &t : tets) {
      const vec3 a = pos_at(t.anchor);
      const vec3 e0 = pos_at(t.v1) - a;
      const vec3 e1 = pos_at(t.v2) - a;
      const vec3 e2 = pos_at(t.v3) - a;
      const real r = eval_fan_tet_volume(e0, e1, e2) - t.V_rest;
      std::vector<vec3> gk(x.x.size(), vec3::Zero());
      fan_tet_grad_slab(t, e0, e1, e2, [&](index_t vi, const vec3 &g) {
        add_grad(vi, g, gk);
      });
      for (int a_i = 0; a_i < n; ++a_i) {
        for (int ax = 0; ax < 3; ++ax) {
          const int ia = 3 * a_i + ax;
          const real ga = gk[static_cast<size_t>(a_i)][ax];
          JTr(ia) += t.weight * r * ga;
          for (int b = 0; b < n; ++b) {
            for (int bx = 0; bx < 3; ++bx) {
              const int ib = 3 * b + bx;
              JTJ(ia, ib) += t.weight * ga * gk[static_cast<size_t>(b)][bx];
            }
          }
        }
      }
    }

    JTJ.diagonal().array() += lm_damp;
    const Eigen::Matrix<real, Eigen::Dynamic, 1> delta = -JTJ.ldlt().solve(JTr);
    if (!delta.allFinite())
      break;

    slab_ring_positions trial = x;
    real alpha = 1.0;
    const real e0 = eval_fan_tet_volume_energy(tets, x, fallback_at);
    for (int ls = 0; ls < 8; ++ls) {
      for (int a_i = 0; a_i < n; ++a_i) {
        for (int ax = 0; ax < 3; ++ax)
          trial.x[static_cast<size_t>(a_i)][ax] =
              x.x[static_cast<size_t>(a_i)][ax] + alpha * delta(3 * a_i + ax);
      }
      const real e1 = eval_fan_tet_volume_energy(tets, trial, fallback_at);
      if (e1 < e0) {
        x = trial;
        break;
      }
      alpha *= 0.5;
    }
    if (e0 < 1e-24)
      break;
  }
  return x;
}

inline real eval_fan_tet_volume_sum(
    const std::vector<slab_fan_tet> &tets,
    const std::function<vec3(index_t)> &pos_at) {
  real sum = 0.0;
  for (const slab_fan_tet &t : tets)
    sum += eval_fan_tet_volume_at(t, pos_at);
  return sum;
}

/// Unpack SL3 spoke targets p into absolute ring positions (foot at c_foot).
inline std::unordered_map<index_t, vec3>
unpack_sl3_spoke_targets(const slab_pierce_stencil &stencil,
                         const std::vector<slab_foot_local_sample> &spokes,
                         const vec3 &c_foot, const vecX &p, index_t row0, real w) {
  std::unordered_map<index_t, vec3> out;
  out[stencil.foot] = c_foot;
  for (size_t k = 0; k < spokes.size(); ++k) {
    const real wk = w * spokes[k].weight;
    if (std::abs(wk) < 1e-24)
      continue;
    const vec3 edge = p.block(row0 + static_cast<index_t>(3 * k), 0, 3, 1) / wk;
    out[spokes[k].vi] = c_foot + edge;
  }
  return out;
}

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // GAUDI_HEPWORTH_SLAB_FAN_TET_HPP
