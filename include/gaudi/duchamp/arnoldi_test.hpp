#ifndef __ARNOLDI_TEST__
#define __ARNOLDI_TEST__

#include "gaudi/duchamp/arnoldi.hpp"
#include "gaudi/vec_addendum.h"

#include "GaudiGraphics/geometry_logger.h"

#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/dynamic.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/asawa/shell/shell.hpp"

#include "gaudi/asawa/primitive_objects.hpp"
#include "gaudi/common.h"

#include "gaudi/kusama/laplacian.hpp"
#include "gaudi/calder/least_squares_fit.hpp"
#include "gaudi/calder/shell_integrators.hpp"
#include "gaudi/calder/tangent_point_integrators.hpp"

#include "gaudi/geometry_logger.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

namespace gaudi {
namespace duchamp {
using namespace asawa;

// Krylov primitives live in gaudi/duchamp/arnoldi.hpp.

enum energy_mode {
  TANGENT_POINT,
  QUADRIC_NORMAL,
  QUADRIC_GRAD,
  CYCLIDE_NORMAL
};

enum solver_mode { POWER_ITERATION, FULL_ARNOLDI };

// ── arnoldi_test class ──────────────────────────────────────────────────

class arnoldi_test {
public:
  typedef std::shared_ptr<arnoldi_test> ptr;

  static ptr create() { return std::make_shared<arnoldi_test>(); }

  arnoldi_test() { load_shell(); };

  void load_shell() {
    __M = shell::load_bunny();
    shell::triangulate(*__M);

    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    asawa::center(x, 2.0);

    real l0 = asawa::shell::avg_length(*__M, x);

    real C = 1.0;
    __surf = shell::dynamic::create(__M, C * l0, 2.5 * C * l0, 0.5 * C * l0);
    //if we do it all at once, the mesh will collapse anisotropically
    while (C <= 5.0) {
      __surf->set_collapse_threshold(C * l0);
      __surf->set_stretch_threshold(2.5 * C * l0);
      __surf->set_bridge_threshold(0.25 * C * l0);
      __surf->step();
      C += 0.25;
    }

    _eps = asawa::shell::avg_length(*__M, x);
    _fit_l0 = 2.0 * _eps;

    __dir_datum_id = asawa::init_vert_datum<vec3>(*__M, vec3::Zero());
  }

  // ── energy gradient methods ─────────────────────────────────────────

  // Mollified TP gradient: product rule / IBP ⇒ G = Gs + Ks (not Gs − Ks).
  //   Gs = κ ∗ ∇E_reg,  Ks = (∇κ) ∗ E_reg;  under IBP, Ks ≈ Gs.
  std::vector<vec3>
  compute_gradient_tangent_point(asawa::shell::shell &M,
                                 const std::vector<vec3> &x) {
    asawa::get_vec_data(M, 0) = x;

    std::vector<vec3> Nv = asawa::shell::vertex_normals(M, x);
    std::vector<real> w = asawa::shell::vertex_areas(M, x);
    real p0 = 6.0;
    real p1 = 2.0;

    std::vector<vec3> Gv =
        calder::tangent_point_gradient(M, x, w, Nv, 1.0 * _eps, p0);
    std::vector<vec3> Gf = asawa::shell::vert_to_face<vec3>(M, x, Gv);
    std::vector<vec3> Gs =
        calder::smoothed_gradient(M, x, Gf, 4.0 * _eps, p1);

    std::vector<real> Kv =
        calder::tangent_point_energy(M, x, w, Nv, 1.0 * _eps, p0);
    std::vector<real> Kf = asawa::shell::vert_to_face<real>(M, x, Kv);
    std::vector<vec3> Ks = calder::gradient_scalar(M, x, Kf, 4.0 * _eps, p1);

    int Nv_count = static_cast<int>(Ks.size());
    std::vector<vec3> G(Nv_count, vec3::Zero());
    for (int i = 0; i < Nv_count; i++) {
      G[i] = Gs[i] + Ks[i];
    }
    return G;
  }

  std::vector<vec3>
  compute_gradient_quadric(asawa::shell::shell &M,
                           const std::vector<vec3> &x,
                           bool use_implicit_grad) {
    asawa::get_vec_data(M, 0) = x;

    std::vector<vec3> Nv = asawa::shell::vertex_normals(M, x);
    int N = static_cast<int>(x.size());

    auto Q = calder::quadric(M, x, Nv, _fit_l0, _fit_p);

    std::vector<vec3> G(N, vec3::Zero());
    for (int i = 0; i < N; i++) {
      if (use_implicit_grad) {
        G[i] = albers::quadric_grad(Q[i], vec3::Zero());
      } else {
        G[i] = Q[i][9] * Nv[i];
      }
      if (G[i].hasNaN())
        G[i] = vec3::Zero();
    }
    return G;
  }

  std::vector<vec3>
  compute_gradient_cyclide(asawa::shell::shell &M,
                           const std::vector<vec3> &x) {
    asawa::get_vec_data(M, 0) = x;

    std::vector<vec3> Nv = asawa::shell::vertex_normals(M, x);
    int N = static_cast<int>(x.size());

    auto Qq = calder::quadric(M, x, Nv, _fit_l0, _fit_p);

    std::vector<vec3> cens(N);
    for (int i = 0; i < N; i++) {
      vec3 dc = albers::quadric_center(Qq[i]);
      if (dc.hasNaN()) {
        cens[i] = x[i];
        continue;
      }
      real lc = dc.norm();
      lc = std::min(lc, 8.0 * _fit_l0);
      if (lc > 1e-10)
        dc = dc.normalized() * lc;
      cens[i] = x[i] + dc;
    }

    auto Qd = calder::darboux_cyclide(M, cens, Nv, _fit_l0, _fit_p);

    std::vector<vec3> G(N, vec3::Zero());
    for (int i = 0; i < N; i++) {
      vec3 dp = x[i] - cens[i];
      real sdf = albers::eval_darboux(Qd[i], dp);
      if (std::isnan(sdf))
        sdf = 0.0;
      G[i] = sdf * Nv[i];
      if (G[i].hasNaN())
        G[i] = vec3::Zero();
    }
    return G;
  }

  std::vector<vec3> compute_gradient(asawa::shell::shell &M,
                                     const std::vector<vec3> &x) {
    switch (_energy_mode) {
    case TANGENT_POINT:
      return compute_gradient_tangent_point(M, x);
    case QUADRIC_NORMAL:
      return compute_gradient_quadric(M, x, false);
    case QUADRIC_GRAD:
      return compute_gradient_quadric(M, x, true);
    case CYCLIDE_NORMAL:
      return compute_gradient_cyclide(M, x);
    default:
      return compute_gradient_tangent_point(M, x);
    }
  }

  // ── logging & visualization ─────────────────────────────────────────

  void log_spectrum(const std::vector<ritz_pair> &pairs,
                    int active_idx) const {
    int np = static_cast<int>(pairs.size());
    real max_abs = 1e-10;
    for (int i = 0; i < np; i++)
      max_abs = std::max(max_abs, std::abs(pairs[i].eigenvalue));

    std::cout << "  arnoldi [k=" << _k << ", m=" << np << "]:" << std::endl;
    for (int i = 0; i < np; i++) {
      real lam = pairs[i].eigenvalue;

      const char *label;
      if (lam < -0.1 * max_abs)
        label = "unstable";
      else if (std::abs(lam) < 0.05 * max_abs)
        label = "soft";
      else if (lam > 0.5 * max_abs)
        label = "stiff";
      else
        label = "moderate";

      std::cout << "    " << (i == active_idx ? ">>>" : "   ") << " lambda_"
                << i << " = " << std::setw(12) << std::fixed
                << std::setprecision(4) << lam << "  (" << label << ")"
                << std::endl;
    }

    const char *mode_name[] = {"DOMINANT_RITZ", "UNSTABLE_RITZ",
                               "REWEIGHTED_GRADIENT"};
    std::cout << "  mode: " << mode_name[_mode]
              << "  target: " << _target_mode << "  active: " << active_idx
              << std::endl;
  }

  // ── main step ───────────────────────────────────────────────────────

  void step(int frame) {
    std::cout << "frame: " << frame << std::endl;
    std::cout << "  -verts: " << __M->vert_count()
              << "  faces: " << __M->face_count() << std::endl;

    __surf->step(false, false);

    asawa::shell::shell &M = *__M;
    const std::vector<vec3> &x0 = asawa::get_vec_data(M, 0);
    int N = static_cast<int>(x0.size());

    std::vector<vec3> G0 = compute_gradient(M, x0);

    real g_norm = krylov_norm(G0);
    if (g_norm < 1e-8) {
      std::cout << "  gradient near zero, skipping" << std::endl;
      return;
    }

    // ── read previous direction from datum, renormalize after remeshing ──
    std::vector<vec3> &d = asawa::get_vec_data(M, __dir_datum_id);
    real d_norm = krylov_norm(d);
    bool have_d = (d_norm > 1e-10);
    if (have_d)
      krylov_scale(1.0 / d_norm, d);

    if (!have_d) {
      for (int i = 0; i < N; i++)
        d[i] = G0[i] / g_norm;
      have_d = true;
    }

    // ── finite-difference Hessian matvec (shared by both solvers) ────────
    real fd_eps = _eps * _fd_eps_scale;
    auto matvec = [&](const std::vector<vec3> &v) -> std::vector<vec3> {
      std::vector<vec3> x_pert(N);
      for (int i = 0; i < N; i++)
        x_pert[i] = x0[i] + fd_eps * v[i];
      auto G_pert = compute_gradient(M, x_pert);
      std::vector<vec3> Hv(N);
      for (int i = 0; i < N; i++)
        Hv[i] = (G_pert[i] - G0[i]) / fd_eps;
      return Hv;
    };

    // ── dispatch to solver ───────────────────────────────────────────────
    mode_estimate est;

    if (_solver_mode == POWER_ITERATION) {
      est = power_iteration_step(matvec, d, _sigma);
      for (int pi = 1; pi < _power_iters; pi++)
        est = power_iteration_step(matvec, est.direction, _sigma);
      std::cout << "  power iter: lambda=" << std::fixed
                << std::setprecision(4) << est.eigenvalue
                << "  sigma=" << _sigma
                << "  iters=" << _power_iters << std::endl;
    } else {
      auto log_fn = [this](const std::vector<ritz_pair> &pairs, int idx) {
        log_spectrum(pairs, idx);
      };
      est = arnoldi_step(matvec, d, have_d, G0, _k, _mode, _target_mode,
                         _alpha, log_fn);
    }

    asawa::get_vec_data(M, 0) = x0;

    // ── sign-resolve and write direction back to datum ────────────────────
    real sign = (krylov_dot(est.direction, d) >= 0) ? 1.0 : -1.0;
    for (int i = 0; i < N; i++)
      d[i] = sign * est.direction[i];

    _lambda_est = est.eigenvalue;

    // ── smooth and apply displacement ────────────────────────────────────
    std::vector<vec3> dir_f =
        asawa::shell::vert_to_face<vec3>(M, x0, d);
    std::vector<vec3> dir_s =
        calder::smoothed_gradient(M, x0, dir_f, 4.0 * _eps, 2.0);

    real dir_norm = krylov_norm(dir_s);
    if (dir_norm < 1e-14)
      return;
    krylov_scale(1.0 / dir_norm, dir_s);

    std::vector<vec3> &x = asawa::get_vec_data(M, 0);
    for (int i = 0; i < N; i++) {
      geometry_logger::line(x[i], x[i] + 1.0 * dir_s[i],
                            vec4(1.0, 0.0, 0.0, 1.0));
      x[i] += 1.0 * _h * dir_s[i];
      //x[i] -= 0.1 * _h * G0[i] / g_norm;
    }
  }

  // ── tunable parameters ────────────────────────────────────────────────

  //solver_mode _solver_mode = FULL_ARNOLDI;
  solver_mode _solver_mode = POWER_ITERATION;
  

  // shared
  index_t __dir_datum_id = -1;
  real _lambda_est = 0.0;
  real _fd_eps_scale = 1e-4;
  real _h = 5e-2;
  real _eps;

  // power iteration
  real _sigma = 5000.0;
  int _power_iters = 8;

  // arnoldi
  step_mode _mode = UNSTABLE_RITZ;
  int _k = 8;
  int _target_mode = 7;
  real _alpha = 1.0;

  // energy
  //energy_mode _energy_mode = QUADRIC_NORMAL;
  energy_mode _energy_mode = TANGENT_POINT;

  // quadric/cyclide fitting parameters
  real _fit_l0;
  real _fit_p = 3.0;

  shell::shell::ptr __M;
  shell::dynamic::ptr __surf;
};

} // namespace duchamp
} // namespace gaudi

#endif
