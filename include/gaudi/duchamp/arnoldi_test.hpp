#ifndef __ARNOLDI_TEST__
#define __ARNOLDI_TEST__

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

// ── vec3-array linear algebra for Krylov methods ────────────────────────

inline real krylov_dot(const std::vector<vec3> &a,
                       const std::vector<vec3> &b) {
  real s = 0.0;
  for (size_t i = 0; i < a.size(); i++)
    s += a[i].dot(b[i]);
  return s;
}

inline real krylov_norm(const std::vector<vec3> &a) {
  return std::sqrt(krylov_dot(a, a));
}

inline void krylov_axpy(real alpha, const std::vector<vec3> &x,
                        std::vector<vec3> &y) {
  for (size_t i = 0; i < x.size(); i++)
    y[i] += alpha * x[i];
}

inline void krylov_scale(real alpha, std::vector<vec3> &x) {
  for (size_t i = 0; i < x.size(); i++)
    x[i] *= alpha;
}

// ── Ritz pair: eigenvalue + eigenvector ─────────────────────────────────

struct ritz_pair {
  real eigenvalue;
  std::vector<vec3> vector;
};

// ── Generic Arnoldi iteration ───────────────────────────────────────────
//
// MatVecFn: callable with signature
//   std::vector<vec3>(const std::vector<vec3>&)
//
// q0: starting vector (typically the gradient)
// k:  subspace dimension (each step costs one matvec call)
//
// Returns Ritz pairs sorted by eigenvalue, largest first.

template <typename MatVecFn>
std::vector<ritz_pair> arnoldi(MatVecFn &&matvec,
                               const std::vector<vec3> &q0, int k) {
  int N = static_cast<int>(q0.size());
  real q0_norm = krylov_norm(q0);
  if (q0_norm < 1e-14)
    return {};

  std::vector<std::vector<vec3>> Q(k + 1);
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(k + 1, k);

  Q[0].resize(N);
  for (int i = 0; i < N; i++)
    Q[0][i] = q0[i] / q0_norm;

  int m = k;
  for (int j = 0; j < k; j++) {
    std::vector<vec3> w = matvec(Q[j]);

    for (int i = 0; i <= j; i++) {
      real h = krylov_dot(w, Q[i]);
      H(i, j) = h;
      krylov_axpy(-h, Q[i], w);
    }

    real h_next = krylov_norm(w);
    H(j + 1, j) = h_next;

    if (h_next < 1e-12) {
      m = j + 1;
      std::cout << "  arnoldi: lucky breakdown at j=" << j << std::endl;
      break;
    }

    Q[j + 1].resize(N);
    for (int i = 0; i < N; i++)
      Q[j + 1][i] = w[i] / h_next;
  }

  Eigen::MatrixXd Hm = H.topLeftCorner(m, m);
  Eigen::EigenSolver<Eigen::MatrixXd> es(Hm);
  auto evals = es.eigenvalues();
  auto evecs = es.eigenvectors();

  std::vector<ritz_pair> pairs(m);
  for (int i = 0; i < m; i++) {
    pairs[i].eigenvalue = evals(i).real();
    pairs[i].vector.resize(N, vec3::Zero());
    for (int j = 0; j < m; j++) {
      real coeff = evecs(j, i).real();
      krylov_axpy(coeff, Q[j], pairs[i].vector);
    }
    real vn = krylov_norm(pairs[i].vector);
    if (vn > 1e-14)
      krylov_scale(1.0 / vn, pairs[i].vector);
  }

  std::sort(pairs.begin(), pairs.end(),
            [](const ritz_pair &a, const ritz_pair &b) {
              return a.eigenvalue > b.eigenvalue;
            });

  return pairs;
}

// ── Stepping modes ──────────────────────────────────────────────────────
//
// DOMINANT_RITZ:      largest |eigenvalue| — stiffest direction
// UNSTABLE_RITZ:      most negative eigenvalue — buckling candidate
// REWEIGHTED_GRADIENT: project G into subspace, weight by |lambda|^alpha

enum step_mode { DOMINANT_RITZ, UNSTABLE_RITZ, REWEIGHTED_GRADIENT };

enum energy_mode {
  TANGENT_POINT,
  QUADRIC_NORMAL,
  QUADRIC_GRAD,
  CYCLIDE_NORMAL
};

enum solver_mode { POWER_ITERATION, FULL_ARNOLDI };

// ── mode_estimate: shared return type for both solvers ──────────────────

struct mode_estimate {
  std::vector<vec3> direction;
  real eigenvalue;
};

// ── Shared primitive: one Hessian-vector product ────────────────────────

struct hv_result {
  std::vector<vec3> w;
  real rayleigh;
};

template <typename MatVecFn>
hv_result compute_Hv(MatVecFn &&matvec, const std::vector<vec3> &v) {
  auto w = matvec(v);
  real rayleigh = krylov_dot(v, w);
  return {std::move(w), rayleigh};
}

// ── Shifted power iteration step ────────────────────────────────────────

template <typename MatVecFn>
mode_estimate power_iteration_step(MatVecFn &&matvec,
                                   const std::vector<vec3> &d, real sigma) {
  auto [w, rayleigh] = compute_Hv(matvec, d);
  int N = static_cast<int>(d.size());
  for (int i = 0; i < N; i++)
    w[i] -= sigma * d[i];

  real sign = (krylov_dot(w, d) >= 0) ? 1.0 : -1.0;
  real w_norm = krylov_norm(w);
  if (w_norm > 1e-14)
    krylov_scale(sign / w_norm, w);

  return {std::move(w), rayleigh};
}

// ── Full Arnoldi step: k matvecs, Ritz extraction, overlap tracking ─────

template <typename MatVecFn>
mode_estimate arnoldi_step(MatVecFn &&matvec, const std::vector<vec3> &d,
                           bool have_d, const std::vector<vec3> &G0, int k,
                           step_mode mode, int target_mode, real alpha,
                           std::function<void(const std::vector<ritz_pair> &,
                                              int)> log_fn) {
  int N = static_cast<int>(G0.size());
  const auto &q0 = have_d ? d : G0;
  auto pairs = arnoldi(matvec, q0, k);

  int active_idx = -1;
  int np = static_cast<int>(pairs.size());

  if (have_d && !pairs.empty()) {
    real best = 0.0;
    for (int i = 0; i < np; i++) {
      real ov = std::abs(krylov_dot(d, pairs[i].vector));
      if (ov > best) {
        best = ov;
        active_idx = i;
      }
    }
  } else {
    if (np == 0) {
      active_idx = -1;
    } else if (target_mode >= 0 && target_mode < np) {
      active_idx = target_mode;
    } else {
      switch (mode) {
      case DOMINANT_RITZ: {
        active_idx = 0;
        real max_abs = 0.0;
        for (int i = 0; i < np; i++) {
          if (std::abs(pairs[i].eigenvalue) > max_abs) {
            max_abs = std::abs(pairs[i].eigenvalue);
            active_idx = i;
          }
        }
        break;
      }
      case UNSTABLE_RITZ: {
        active_idx = np - 1;
        for (int i = 0; i < np; i++) {
          if (pairs[i].eigenvalue < 0) {
            active_idx = i;
            break;
          }
        }
        break;
      }
      case REWEIGHTED_GRADIENT:
        active_idx = -1;
        break;
      }
    }
  }

  if (log_fn)
    log_fn(pairs, active_idx);

  std::vector<vec3> dir;
  if (mode == REWEIGHTED_GRADIENT) {
    dir.assign(N, vec3::Zero());
    if (!pairs.empty()) {
      for (int i = 0; i < np; i++) {
        real coeff = krylov_dot(G0, pairs[i].vector);
        real weight = std::pow(std::abs(pairs[i].eigenvalue) + 1e-10, alpha);
        krylov_axpy(coeff * weight, pairs[i].vector, dir);
      }
    } else {
      dir = G0;
    }
  } else if (active_idx >= 0 && active_idx < np) {
    dir = pairs[active_idx].vector;
  } else {
    dir = G0;
  }

  real dn = krylov_norm(dir);
  if (dn > 1e-14)
    krylov_scale(1.0 / dn, dir);

  real lambda = (active_idx >= 0 && active_idx < np)
                    ? pairs[active_idx].eigenvalue
                    : 0.0;

  return {std::move(dir), lambda};
}

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

  // E = w(x) * F(x)
  // dE/dx = dw/dx * F(x) + w(x) * dF/dx
  //       = Ks            + Gs
  // G = Gs - Ks
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
      G[i] = Gs[i] - Ks[i];
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
