#ifndef GAUDI_HEPWORTH_SLAB_WEDGE_STRAIN_HPP
#define GAUDI_HEPWORTH_SLAB_WEDGE_STRAIN_HPP

#include "bivector3.hpp"
#include "slab_cap_solid_angle_generated.hpp"
#include "slab_fan_wedge.hpp"
#include "slab_prism_operators.hpp"
#include "trivector3.hpp"

#include <Eigen/Dense>
#include <cmath>

namespace gaudi {
namespace hepworth {
namespace block {

inline vec3 cap_hub_fan_axis(const vec3 &c_hub,
                             const std::vector<slab_fan_wedge> &wedges,
                             const slab_ring_positions &ring, bool bottom_cap) {
  vec3 sum = vec3::Zero();
  for (const slab_fan_wedge &w : wedges) {
    if (w.bottom_cap != bottom_cap)
      continue;
    const vec3 e0 = ring.at(w.i0) - c_hub;
    const vec3 e1 = ring.at(w.i1) - c_hub;
    sum += e0.cross(e1);
  }
  const real n2 = sum.squaredNorm();
  if (n2 < 1e-24)
    return vec3::Zero();
  return sum / std::sqrt(n2);
}

inline real cap_wedge_solid_angle_signed(const vec3 &e0, const vec3 &e1,
                                         const vec3 &fan_axis) {
  const vec3 cr = e0.cross(e1);
  const real den = e0.dot(e1) + e0.norm() * e1.norm();
  if (den <= 1e-24)
    return 0.0;
  if (fan_axis.squaredNorm() < 1e-24)
    return cap_generated::cap_wedge_solid_angle(e0, e1);
  return 2.0 * std::atan2(cr.dot(fan_axis), den);
}

inline real cap_wedge_solid_angle_oriented(const vec3 &e0, const vec3 &e1,
                                           const vec3 &fan_axis) {
  const real s = cap_wedge_solid_angle_signed(e0, e1, fan_axis);
  return s >= 0.0 ? s : -s;
}

inline real cap_wedge_solid_angle_at_hub(
    const vec3 &c_hub, const std::vector<slab_fan_wedge> &wedges,
    const slab_ring_positions &ring, bool bottom_cap, const vec3 &e0,
    const vec3 &e1) {
  const vec3 axis = cap_hub_fan_axis(c_hub, wedges, ring, bottom_cap);
  return cap_wedge_solid_angle_oriented(e0, e1, axis);
}

inline trivector3 eval_fan_wedge_trivector(index_t foot, const vec3 &c_hub,
                                           index_t i0, index_t i1,
                                           const vec3 &e0, const vec3 &e1,
                                           const vec3 &d, bool bottom_cap) {
  const bivector3 B = bivector3::wedge(e0, e1);
  const trivector3 tau = trivector3::wedge(B, bottom_cap ? -d : d);
  (void)foot;
  (void)c_hub;
  (void)i0;
  (void)i1;
  return tau;
}

inline std::vector<slab_fan_wedge>
build_slab_fan_wedges(const slab_pierce_stencil &stencil,
                      const std::function<vec3(index_t)> &pos_at) {
  std::vector<slab_fan_wedge> out;
  const vec3 c_top = pos_at(stencil.foot);
  const vec3 c_bot = pos_at(stencil.pierce_foot);
  const vec3 d = c_bot - c_top;

  const size_t n_top = stencil.top_ids.size() - 1;
  out.reserve(n_top + stencil.bot_ids.size());
  for (size_t k = 0; k < n_top; ++k) {
    slab_fan_wedge w;
    w.hub = stencil.foot;
    w.i0 = stencil.top_ids[k + 1];
    w.i1 = stencil.top_ids[(k + 1) % n_top + 1];
    w.weight = stencil.top_weights[k];
    w.bottom_cap = false;
    const vec3 e0 = pos_at(w.i0) - c_top;
    const vec3 e1 = pos_at(w.i1) - c_top;
    w.tau_rest =
        eval_fan_wedge_trivector(stencil.foot, c_top, w.i0, w.i1, e0, e1, d, false);
    out.push_back(w);
  }

  const size_t n_bot = stencil.bot_ids.size() - 1;
  for (size_t k = 0; k < n_bot; ++k) {
    slab_fan_wedge w;
    w.hub = stencil.pierce_foot;
    w.i0 = stencil.bot_ids[k + 1];
    w.i1 = stencil.bot_ids[(k + 1) % n_bot + 1];
    w.weight = stencil.bot_weights[k];
    w.bottom_cap = true;
    const vec3 e0 = pos_at(w.i0) - c_bot;
    const vec3 e1 = pos_at(w.i1) - c_bot;
    w.tau_rest =
        eval_fan_wedge_trivector(stencil.foot, c_bot, w.i0, w.i1, e0, e1, d, true);
    out.push_back(w);
  }

  const auto pos_ring = [&](index_t vi) { return pos_at(vi); };
  const slab_ring_positions ring =
      slab_ring_positions::from_stencil(stencil, pos_ring);
  for (slab_fan_wedge &w : out) {
    const vec3 c = w.bottom_cap ? c_bot : c_top;
    const vec3 e0 = pos_at(w.i0) - c;
    const vec3 e1 = pos_at(w.i1) - c;
    w.omega_rest =
        cap_wedge_solid_angle_at_hub(c, out, ring, w.bottom_cap, e0, e1);
  }
  return out;
}

inline std::vector<slab_fan_wedge>
build_slab_fan_wedges(const slab_pierce_stencil &stencil,
                      const std::vector<vec3> &x_rest) {
  const auto pos_at = [&](index_t vi) { return x_rest[static_cast<size_t>(vi)]; };
  return build_slab_fan_wedges(stencil, pos_at);
}

inline void fan_wedge_grad_slab(
    index_t foot, index_t pierce_foot, const slab_fan_wedge &w, const vec3 &e0,
    const vec3 &e1, const vec3 &d,
    const std::function<void(index_t, const vec3 &)> &add) {
  const vec3 ed = w.bottom_cap ? -d : d;
  const vec3 g0 = e1.cross(ed);
  const vec3 g1 = ed.cross(e0);
  const vec3 cross_e = e0.cross(e1);
  add(w.i0, g0);
  add(w.i1, g1);
  if (!w.bottom_cap) {
    add(foot, -(g0 + g1) - cross_e);
    add(pierce_foot, cross_e);
  } else {
    add(pierce_foot, -(g0 + g1) - cross_e);
    add(foot, cross_e);
  }
}

inline real eval_wedge_strain_energy(
    const slab_pierce_stencil &stencil,
    const std::vector<slab_fan_wedge> &wedges, const slab_ring_positions &ring,
    const std::function<vec3(index_t)> &fallback_at) {
  const vec3 c_top = ring.at(stencil.foot);
  const vec3 c_bot = ring.at(stencil.pierce_foot);
  const vec3 d = c_bot - c_top;

  const auto pos_at = [&](index_t vi) -> vec3 {
    const auto it = ring.idx.find(vi);
    if (it != ring.idx.end())
      return ring.x[it->second];
    return fallback_at(vi);
  };

  real E = 0.0;
  for (const slab_fan_wedge &w : wedges) {
    const vec3 c = w.bottom_cap ? c_bot : c_top;
    const vec3 e0 = pos_at(w.i0) - c;
    const vec3 e1 = pos_at(w.i1) - c;
    const trivector3 tau =
        eval_fan_wedge_trivector(stencil.foot, c, w.i0, w.i1, e0, e1, d, w.bottom_cap);
    const real r = tau.coeff - w.tau_rest.coeff;
    E += w.weight * r * r;
  }
  return E;
}

inline void wedge_strain_grad(
    const slab_pierce_stencil &stencil,
    const std::vector<slab_fan_wedge> &wedges, const slab_ring_positions &ring,
    const std::function<vec3(index_t)> &fallback_at, std::vector<vec3> &grad) {
  grad.assign(ring.x.size(), vec3::Zero());
  const vec3 c_top = ring.at(stencil.foot);
  const vec3 c_bot = ring.at(stencil.pierce_foot);
  const vec3 d = c_bot - c_top;

  const auto pos_at = [&](index_t vi) -> vec3 {
    const auto it = ring.idx.find(vi);
    if (it != ring.idx.end())
      return ring.x[it->second];
    return fallback_at(vi);
  };

  const auto add = [&](index_t vi, const vec3 &g) {
    const auto it = ring.idx.find(vi);
    if (it != ring.idx.end())
      grad[it->second] += g;
  };

  for (const slab_fan_wedge &w : wedges) {
    const vec3 c = w.bottom_cap ? c_bot : c_top;
    const vec3 e0 = pos_at(w.i0) - c;
    const vec3 e1 = pos_at(w.i1) - c;
    const trivector3 tau =
        eval_fan_wedge_trivector(stencil.foot, c, w.i0, w.i1, e0, e1, d, w.bottom_cap);
    const real r = tau.coeff - w.tau_rest.coeff;
    std::vector<vec3> wg(ring.x.size(), vec3::Zero());
    const auto add_w = [&](index_t vi, const vec3 &g) {
      const auto it = ring.idx.find(vi);
      if (it != ring.idx.end())
        wg[it->second] += g;
    };
    fan_wedge_grad_slab(stencil.foot, stencil.pierce_foot, w, e0, e1, d, add_w);
    for (size_t i = 0; i < grad.size(); ++i)
      grad[i] += real(2.0) * w.weight * r * wg[i];
  }
}

inline real wedge_strain_energy_scale(
    const std::vector<slab_fan_wedge> &wedges) {
  real scale = 0.0;
  for (const slab_fan_wedge &w : wedges) {
    const real t = w.tau_rest.coeff;
    scale += w.weight * t * t;
  }
  return std::max(scale, real(1e-24));
}

/// Gauss–Newton on Σ w_k (τ_k − τ_k^rest)² (top + bottom ring DOFs).
inline slab_ring_positions resolve_wedge_strain_gauss_newton(
    const slab_pierce_stencil &stencil, const std::vector<slab_fan_wedge> &wedges,
    const slab_ring_positions &x0,
    const std::function<vec3(index_t)> &fallback_at, int iters,
    real lm_damp = 1e-8, real strain_tol_frac = 1e-10) {
  slab_ring_positions x = x0;
  const int n = static_cast<int>(x.x.size());
  const int nd = 3 * n;
  if (nd <= 0 || wedges.empty() || iters <= 0)
    return x;

  const real strain_tol =
      strain_tol_frac * wedge_strain_energy_scale(wedges);
  const real e_init = eval_wedge_strain_energy(stencil, wedges, x, fallback_at);
  if (e_init <= strain_tol)
    return x;

  Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> JTJ(nd, nd);
  Eigen::Matrix<real, Eigen::Dynamic, 1> JTr(nd);
  Eigen::Matrix<real, Eigen::Dynamic, 1> delta(nd);
  std::vector<vec3> gk(static_cast<size_t>(n));

  for (int iter = 0; iter < iters; ++iter) {
    const vec3 c_top = x.at(stencil.foot);
    const vec3 c_bot = x.at(stencil.pierce_foot);
    const vec3 d = c_bot - c_top;

    const auto pos_at = [&](index_t vi) -> vec3 {
      const auto it = x.idx.find(vi);
      if (it != x.idx.end())
        return x.x[it->second];
      return fallback_at(vi);
    };

    JTJ.setZero();
    JTr.setZero();

    const auto add_grad = [&](index_t vi, const vec3 &g) {
      const auto it = x.idx.find(vi);
      if (it != x.idx.end())
        gk[static_cast<size_t>(it->second)] += g;
    };

    for (const slab_fan_wedge &w : wedges) {
      for (vec3 &g : gk)
        g.setZero();
      const vec3 c = w.bottom_cap ? c_bot : c_top;
      const vec3 e0 = pos_at(w.i0) - c;
      const vec3 e1 = pos_at(w.i1) - c;
      const trivector3 tau =
          eval_fan_wedge_trivector(stencil.foot, c, w.i0, w.i1, e0, e1, d,
                                   w.bottom_cap);
      const real r = tau.coeff - w.tau_rest.coeff;
      fan_wedge_grad_slab(stencil.foot, stencil.pierce_foot, w, e0, e1, d,
                          add_grad);
      for (int a = 0; a < n; ++a) {
        for (int ax = 0; ax < 3; ++ax) {
          const int ia = 3 * a + ax;
          const real ga = gk[static_cast<size_t>(a)][ax];
          JTr(ia) += w.weight * r * ga;
          for (int b = 0; b < n; ++b) {
            for (int bx = 0; bx < 3; ++bx) {
              const int ib = 3 * b + bx;
              JTJ(ia, ib) += w.weight * ga * gk[static_cast<size_t>(b)][bx];
            }
          }
        }
      }
    }

    JTJ.diagonal().array() += lm_damp;
    delta = -JTJ.ldlt().solve(JTr);
    if (!delta.allFinite())
      break;

    slab_ring_positions trial = x;
    real alpha = 1.0;
    const real e0 = eval_wedge_strain_energy(stencil, wedges, x, fallback_at);
    for (int ls = 0; ls < 8; ++ls) {
      for (int a = 0; a < n; ++a) {
        for (int ax = 0; ax < 3; ++ax)
          trial.x[static_cast<size_t>(a)][ax] =
              x.x[static_cast<size_t>(a)][ax] +
              alpha * delta(3 * a + ax);
      }
      const real e1 = eval_wedge_strain_energy(stencil, wedges, trial, fallback_at);
      if (e1 < e0) {
        x = trial;
        break;
      }
      alpha *= 0.5;
    }
    if (e0 <= strain_tol)
      break;
  }
  return x;
}

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // GAUDI_HEPWORTH_SLAB_WEDGE_STRAIN_HPP
