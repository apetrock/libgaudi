#ifndef GAUDI_PRINCIPAL_CURVATURE_TESTS_HPP
#define GAUDI_PRINCIPAL_CURVATURE_TESTS_HPP

#include "gaudi/asawa/primitive_objects.hpp"
#include "gaudi/kusama/principal_curvature_solver.hpp"
#include "gaudi/test/test.hpp"

#include <cmath>
#include <limits>
#include <random>
#include <vector>

namespace gaudi {
namespace test {
namespace {

using kusama::principal_curvature_case;
using kusama::principal_curvatures;
using kusama::principal_curvature_at;
using kusama::principal_curvature_jacobian;
using kusama::principal_curvature_solver;
using kusama::principal_curvature_stencil;
using kusama::tri3;

/// `asawa::make_torus` (outward winding) as a triangle list. `phi[i]` is the
/// tube angle of vertex `i`.
inline void eig_torus(real R, real r, int major, int minor,
                      std::vector<vec3> &x, std::vector<tri3> &tris,
                      std::vector<real> &phi) {
  std::vector<std::vector<int>> faces;
  asawa::make_torus(x, faces, R, r, major, minor);
  tris.clear();
  for (const auto &f : faces)
    tris.push_back({f[0], f[1], f[2]});
  phi.resize(x.size());
  for (size_t i = 0; i < x.size(); ++i)
    phi[i] = 2.0 * M_PI * real(static_cast<int>(i) % minor) / real(minor);
}

/// Torus principal curvatures with outward normals: tube `1/r`, and
/// `cos(phi) / (R + r cos(phi))` along the ring.
inline std::vector<vec2> eig_torus_kappa(real R, real r,
                                         const std::vector<real> &phi) {
  std::vector<vec2> k(phi.size());
  for (size_t i = 0; i < phi.size(); ++i)
    k[i] = vec2(1.0 / r, std::cos(phi[i]) / (R + r * std::cos(phi[i])));
  return k;
}

/// Ring and tube radius of a torus about +Z.
inline vec2 eig_fit_torus(const std::vector<vec3> &x) {
  real rho_min = 1e30, rho_max = 0;
  for (const vec3 &p : x) {
    const real rho = std::hypot(p[0], p[1]);
    rho_min = std::min(rho_min, rho);
    rho_max = std::max(rho_max, rho);
  }
  const real Rc = 0.5 * (rho_min + rho_max);
  real tube = 0;
  for (const vec3 &p : x)
    tube += std::hypot(std::hypot(p[0], p[1]) - Rc, p[2]);
  return vec2(Rc, tube / real(x.size()));
}

/// Small random offsets so no vertex sits on an exact symmetry.
inline void eig_jitter(std::vector<vec3> &x, real amp, unsigned seed) {
  std::mt19937 rng(seed);
  std::uniform_real_distribution<real> u(-amp, amp);
  for (vec3 &p : x)
    p += vec3(u(rng), u(rng), u(rng));
}

/// Runs LM until it stops making progress. Returns accepted iterations.
inline int eig_solve(principal_curvature_solver &s, std::vector<vec3> &x,
                     const std::vector<vec2> &kh, int max_it,
                     bool *monotone = nullptr) {
  int accepted = 0;
  int stalls = 0;
  for (int it = 0; it < max_it && stalls < 4; ++it) {
    const kusama::principal_curvature_report rep = s.iterate(x, kh);
    if (monotone && rep.accepted && !(rep.energy1 < rep.energy0))
      *monotone = false;
    if (!rep.accepted) {
      ++stalls;
      continue;
    }
    stalls = 0;
    ++accepted;
    if (rep.energy0 - rep.energy1 <= 1e-12 * rep.energy0)
      break;
  }
  return accepted;
}

GAUDI_TEST(principal_curvature_torus_exact) {
  const real R = 1.0, r = 0.33;
  std::vector<vec3> x;
  std::vector<tri3> tris;
  std::vector<real> phi;
  eig_torus(R, r, 48, 24, x, tris, phi);
  principal_curvature_solver s;
  s.reset(static_cast<int>(x.size()), tris, x);
  const std::vector<vec2> k = s.principal(x);
  const std::vector<vec2> ka = eig_torus_kappa(R, r, phi);
  real e1 = 0, e2 = 0;
  for (size_t i = 0; i < k.size(); ++i) {
    GAUDI_ASSERT(k[i].allFinite());
    e1 += std::abs(k[i][0] - ka[i][0]) / ka[i][0];
    e2 += std::abs(k[i][1] - ka[i][1]);
  }
  e1 /= real(k.size());
  e2 /= real(k.size());
  console_logger::info << "torus kappa: mean rel err k1 " << e1
                       << ", mean abs err k2 " << e2 << std::endl;
  GAUDI_EXPECT(e1 < 0.03);
  GAUDI_EXPECT(e2 < 0.05);
}

GAUDI_TEST(principal_curvature_dihedral_gradient_matches_fd) {
  std::vector<vec3> x;
  std::vector<tri3> tris;
  std::vector<real> phi;
  eig_torus(1.0, 0.33, 24, 12, x, tris, phi);
  eig_jitter(x, 0.01, 3);
  const kusama::principal_curvature_mesh m =
      kusama::principal_curvature_mesh::build(static_cast<int>(x.size()), tris);
  const real h = 1e-6;
  real worst = 0;
  for (int v = 0; v < m.n; v += 7)
    for (const kusama::principal_curvature_spoke &s : m.spokes[v]) {
      const int c = kusama::principal_curvature_third(m.tris[s.f0], v, s.j);
      const int d = kusama::principal_curvature_third(m.tris[s.f1], v, s.j);
      const auto g = kusama::principal_curvature_dihedral_gradient(x[v], x[s.j], x[c],
                                                         x[d]);
      const int ids[4] = {v, s.j, c, d};
      const vec3 ga[4] = {g.ga, g.gb, g.gc, g.gd};
      for (int q = 0; q < 4; ++q)
        for (int ax = 0; ax < 3; ++ax) {
          const real x0 = x[ids[q]][ax];
          x[ids[q]][ax] = x0 + h;
          const real bp = kusama::principal_curvature_dihedral(x, m, v, s.j, s.f0, s.f1);
          x[ids[q]][ax] = x0 - h;
          const real bm = kusama::principal_curvature_dihedral(x, m, v, s.j, s.f0, s.f1);
          x[ids[q]][ax] = x0;
          const real fd = (bp - bm) / (2 * h);
          worst = std::max(worst, std::abs(fd - ga[q][ax]) / (1 + std::abs(fd)));
        }
    }
  console_logger::info << "dihedral gradient worst rel err " << worst
                       << std::endl;
  GAUDI_EXPECT(worst < 1e-6);
}

GAUDI_TEST(principal_curvature_tensor_jacobian_matches_fd) {
  std::vector<vec3> x;
  std::vector<tri3> tris;
  std::vector<real> phi;
  eig_torus(1.0, 0.33, 24, 12, x, tris, phi);
  eig_jitter(x, 0.01, 5);
  const kusama::principal_curvature_mesh m =
      kusama::principal_curvature_mesh::build(static_cast<int>(x.size()), tris);
  real worst = 0, scale = 0;
  for (int v = 0; v < m.n; v += 5) {
    const std::vector<int> st = principal_curvature_stencil(m, v);
    std::vector<mat3> da, df;
    GAUDI_ASSERT(kusama::principal_curvature_tensor_jacobian(x, m, v, st, da));
    GAUDI_ASSERT(kusama::principal_curvature_tensor_jacobian_fd(x, m, v, st, 1e-6, df));
    for (size_t i = 0; i < da.size(); ++i) {
      worst = std::max(worst, (da[i] - df[i]).cwiseAbs().maxCoeff());
      scale = std::max(scale, df[i].cwiseAbs().maxCoeff());
    }
  }
  console_logger::info << "dK worst abs err " << worst << " (scale " << scale
                       << ")" << std::endl;
  GAUDI_EXPECT(worst < 1e-6 * (1 + scale));
}

/// The whole stacked Jacobian (curvature, conformal, areal, closeness rows):
/// closed forms against central differences.
GAUDI_TEST(principal_curvature_solver_jacobian_matches_fd) {
  std::vector<vec3> x;
  std::vector<tri3> tris;
  std::vector<real> phi;
  eig_torus(1.0, 0.33, 24, 12, x, tris, phi);
  principal_curvature_solver s;
  s.weights.k_areal = 1.0;
  s.reset(static_cast<int>(x.size()), tris, x);
  eig_jitter(x, 0.01, 9);
  const std::vector<vec2> kh = eig_torus_kappa(1.2, 0.4, phi);
  vecX ra, rf;
  matS Ja, Jf;
  s.analytic = true;
  s.linearize(x, kh, ra, Ja);
  s.analytic = false;
  s.linearize(x, kh, rf, Jf);
  GAUDI_ASSERT(ra.size() == rf.size());
  GAUDI_EXPECT((ra - rf).cwiseAbs().maxCoeff() < 1e-12);
  const matS D = Ja - Jf;
  real worst = 0, scale = 0;
  for (int k = 0; k < D.outerSize(); ++k)
    for (matS::InnerIterator it(D, k); it; ++it)
      worst = std::max(worst, std::abs(it.value()));
  for (int k = 0; k < Jf.outerSize(); ++k)
    for (matS::InnerIterator it(Jf, k); it; ++it)
      scale = std::max(scale, std::abs(it.value()));
  console_logger::info << "solver J worst abs err " << worst << " (scale "
                       << scale << ")" << std::endl;
  GAUDI_EXPECT(worst < 1e-5 * (1 + scale));
}

/// A.4 against differences of the read-out itself, at distinct-eigenvalue
/// vertices (the cylindric / umbilic cases are approximations by design).
GAUDI_TEST(principal_curvature_derivative_matches_fd) {
  std::vector<vec3> x;
  std::vector<tri3> tris;
  std::vector<real> phi;
  eig_torus(1.0, 0.33, 48, 24, x, tris, phi);
  principal_curvature_solver s;
  s.reset(static_cast<int>(x.size()), tris, x);
  const auto &m = s.mesh();
  const real h = 1e-6;
  real worst = 0;
  int checked = 0;
  for (int v = 5; v < m.n; v += 37) {
    const principal_curvatures pr = principal_curvature_at(x, m, v);
    GAUDI_ASSERT(pr.ok);
    if (pr.kind != principal_curvature_case::distinct ||
        std::abs(pr.k2) < 0.05 * std::abs(pr.k1))
      continue;
    ++checked;
    const std::vector<int> st = principal_curvature_stencil(m, v);
    std::vector<mat3> dK;
    GAUDI_ASSERT(kusama::principal_curvature_tensor_jacobian(x, m, v, st, dK));
    Eigen::MatrixXd J;
    GAUDI_ASSERT(principal_curvature_jacobian(pr, dK, J));
    for (int si = 0; si < static_cast<int>(st.size()); ++si)
      for (int c = 0; c < 3; ++c) {
        const real x0 = x[st[si]][c];
        x[st[si]][c] = x0 + h;
        const principal_curvatures pp = principal_curvature_at(x, m, v);
        x[st[si]][c] = x0 - h;
        const principal_curvatures pm = principal_curvature_at(x, m, v);
        x[st[si]][c] = x0;
        const real d1 = (pp.k1 - pm.k1) / (2 * h);
        const real d2 = (pp.k2 - pm.k2) / (2 * h);
        worst = std::max(worst, std::abs(d1 - J(0, 3 * si + c)) /
                                    (1.0 + std::abs(d1)));
        worst = std::max(worst, std::abs(d2 - J(1, 3 * si + c)) /
                                    (1.0 + std::abs(d2)));
      }
  }
  console_logger::info << "kappa derivative worst rel err " << worst
                       << " over " << checked << " vertices" << std::endl;
  GAUDI_EXPECT(checked > 10);
  GAUDI_EXPECT(worst < 1e-4);
}

/// The Eigensatz energy (closeness off) is rigid invariant, so an exact
/// gradient carries no net force or torque. Reports how far ours is off
/// (cylindric and umbilic vertices use the paper's approximate derivative),
/// and checks `principal_curvature_remove_rigid` removes it.
GAUDI_TEST(principal_curvature_gradient_rigid_moments) {
  std::vector<vec3> x;
  std::vector<tri3> tris;
  std::vector<real> phi;
  eig_torus(1.0, 0.33, 48, 24, x, tris, phi);
  principal_curvature_solver s;
  s.weights.k_close = 0;
  s.reset(static_cast<int>(x.size()), tris, x);
  eig_jitter(x, 0.005, 21);
  const std::vector<vec2> kh = eig_torus_kappa(1.2, 0.396, phi);
  vecX r;
  matS J;
  s.linearize(x, kh, r, J);
  vecX g = J.transpose() * r;
  const std::vector<real> m(x.size(), 1.0);
  vec3 F, T;
  kusama::principal_curvature_rigid_moments(x, m, g, F, T);
  const real gn = g.norm();
  console_logger::info << "principal_curvature gradient |g| " << gn << "  |sum g| "
                       << F.norm() << "  |torque| " << T.norm() << std::endl;
  kusama::principal_curvature_remove_rigid(x, m, g);
  kusama::principal_curvature_rigid_moments(x, m, g, F, T);
  console_logger::info << "after remove_rigid: |sum g| " << F.norm()
                       << "  |torque| " << T.norm() << "  |g| " << g.norm()
                       << std::endl;
  GAUDI_EXPECT(F.norm() < 1e-9 * gn);
  GAUDI_EXPECT(T.norm() < 1e-9 * gn);
}

GAUDI_TEST(principal_curvature_idle_when_targets_match) {
  std::vector<vec3> x;
  std::vector<tri3> tris;
  std::vector<real> phi;
  eig_torus(1.0, 0.33, 48, 24, x, tris, phi);
  principal_curvature_solver s;
  s.reset(static_cast<int>(x.size()), tris, x);
  const std::vector<vec2> kh = s.principal(x);
  const std::vector<vec3> x0 = x;
  const kusama::principal_curvature_report rep = s.iterate(x, kh);
  GAUDI_EXPECT(rep.energy0 < 1e-20);
  real moved = 0;
  for (size_t i = 0; i < x.size(); ++i)
    moved = std::max(moved, (x[i] - x0[i]).norm());
  GAUDI_EXPECT(moved < 1e-12);
}

/// Targets of the torus scaled by `scale`. A similarity keeps every angle,
/// so the conformal term agrees with the targets and the answer is exact.
inline vec2 eig_scale_solve(real scale, bool prefactor, int *accepted,
                            bool *monotone) {
  const real R = 1.0, r = 0.33;
  std::vector<vec3> x;
  std::vector<tri3> tris;
  std::vector<real> phi;
  eig_torus(R, r, 48, 24, x, tris, phi);
  principal_curvature_solver s;
  s.prefactor = prefactor;
  s.reset(static_cast<int>(x.size()), tris, x);
  const std::vector<vec2> kh = eig_torus_kappa(R * scale, r * scale, phi);
  *accepted = eig_solve(s, x, kh, 400, monotone);
  return eig_fit_torus(x);
}

GAUDI_TEST(principal_curvature_torus_scale_change) {
  int accepted = 0;
  bool monotone = true;
  const vec2 fit = eig_scale_solve(1.2, true, &accepted, &monotone);
  console_logger::info << "torus x1.2: R " << fit[0] << " tube " << fit[1]
                       << " accepted " << accepted << std::endl;
  GAUDI_EXPECT(monotone);
  GAUDI_EXPECT(std::abs(fit[0] - 1.2) / 1.2 < 0.02);
  GAUDI_EXPECT(std::abs(fit[1] - 0.396) / 0.396 < 0.02);
}

GAUDI_TEST(principal_curvature_torus_shrink) {
  int accepted = 0;
  bool monotone = true;
  const vec2 fit = eig_scale_solve(0.9, true, &accepted, &monotone);
  console_logger::info << "torus x0.9: R " << fit[0] << " tube " << fit[1]
                       << " accepted " << accepted << std::endl;
  GAUDI_EXPECT(monotone);
  GAUDI_EXPECT(std::abs(fit[0] - 0.9) / 0.9 < 0.02);
  GAUDI_EXPECT(std::abs(fit[1] - 0.297) / 0.297 < 0.02);
}

GAUDI_TEST(principal_curvature_prefactor_matches_full_lm) {
  int a_pre = 0, a_full = 0;
  bool mono = true;
  const vec2 pre = eig_scale_solve(1.2, true, &a_pre, &mono);
  const vec2 full = eig_scale_solve(1.2, false, &a_full, &mono);
  console_logger::info << "prefactor R " << pre[0] << " tube " << pre[1]
                       << " (" << a_pre << " its); full LM R " << full[0]
                       << " tube " << full[1] << " (" << a_full << " its)"
                       << std::endl;
  GAUDI_EXPECT(mono);
  GAUDI_EXPECT(std::abs(pre[0] - full[0]) / full[0] < 0.005);
  GAUDI_EXPECT(std::abs(pre[1] - full[1]) / full[1] < 0.005);
}

} // namespace
} // namespace test
} // namespace gaudi

#endif
