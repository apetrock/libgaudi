#ifndef GAUDI_DUCHAMP_ARNOLDI_HPP
#define GAUDI_DUCHAMP_ARNOLDI_HPP

/// Arnoldi on a vec3 field. The matvec returns a vector of the same length.
/// Ritz pairs come back sorted by eigenvalue, largest first.

#include "gaudi/common.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>

namespace gaudi {
namespace duchamp {

inline real krylov_dot(const std::vector<vec3> &a, const std::vector<vec3> &b) {
  real s = 0.0;
  const size_t n = std::min(a.size(), b.size());
  for (size_t i = 0; i < n; ++i)
    s += a[i].dot(b[i]);
  return s;
}

inline real krylov_norm(const std::vector<vec3> &a) {
  return std::sqrt(krylov_dot(a, a));
}

inline void krylov_axpy(real alpha, const std::vector<vec3> &x,
                        std::vector<vec3> &y) {
  const size_t n = std::min(x.size(), y.size());
  for (size_t i = 0; i < n; ++i)
    y[i] += alpha * x[i];
}

inline void krylov_scale(real alpha, std::vector<vec3> &x) {
  for (vec3 &v : x)
    v *= alpha;
}

struct ritz_pair {
  real eigenvalue = 0;
  std::vector<vec3> vector;
};

/// k matvecs. An invariant subspace returns fewer than k pairs.
template <typename MatVecFn>
std::vector<ritz_pair> arnoldi(MatVecFn &&matvec, const std::vector<vec3> &q0,
                               int k) {
  const int N = static_cast<int>(q0.size());
  const real q0_norm = krylov_norm(q0);
  if (N <= 0 || k <= 0 || !(q0_norm > real(1e-14)))
    return {};

  std::vector<std::vector<vec3>> Q(static_cast<size_t>(k) + 1);
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(k + 1, k);

  Q[0].resize(static_cast<size_t>(N));
  for (int i = 0; i < N; ++i)
    Q[0][static_cast<size_t>(i)] = q0[static_cast<size_t>(i)] / q0_norm;

  int m = k;
  for (int j = 0; j < k; ++j) {
    std::vector<vec3> w = matvec(Q[static_cast<size_t>(j)]);
    if (static_cast<int>(w.size()) != N)
      w.assign(static_cast<size_t>(N), vec3::Zero());

    for (int i = 0; i <= j; ++i) {
      const real h = krylov_dot(w, Q[static_cast<size_t>(i)]);
      H(i, j) = h;
      krylov_axpy(-h, Q[static_cast<size_t>(i)], w);
    }

    const real h_next = krylov_norm(w);
    H(j + 1, j) = h_next;
    if (!(h_next > real(1e-12))) {
      m = j + 1;
      break;
    }

    Q[static_cast<size_t>(j) + 1].resize(static_cast<size_t>(N));
    for (int i = 0; i < N; ++i)
      Q[static_cast<size_t>(j) + 1][static_cast<size_t>(i)] =
          w[static_cast<size_t>(i)] / h_next;
  }

  const Eigen::MatrixXd Hm = H.topLeftCorner(m, m);
  Eigen::EigenSolver<Eigen::MatrixXd> es(Hm);
  const auto evals = es.eigenvalues();
  const auto evecs = es.eigenvectors();

  std::vector<ritz_pair> pairs(static_cast<size_t>(m));
  for (int i = 0; i < m; ++i) {
    pairs[static_cast<size_t>(i)].eigenvalue = evals(i).real();
    pairs[static_cast<size_t>(i)].vector.assign(static_cast<size_t>(N),
                                                vec3::Zero());
    for (int j = 0; j < m; ++j) {
      const real coeff = evecs(j, i).real();
      krylov_axpy(coeff, Q[static_cast<size_t>(j)],
                  pairs[static_cast<size_t>(i)].vector);
    }
    const real vn = krylov_norm(pairs[static_cast<size_t>(i)].vector);
    if (vn > real(1e-14))
      krylov_scale(real(1) / vn, pairs[static_cast<size_t>(i)].vector);
  }

  std::sort(pairs.begin(), pairs.end(),
            [](const ritz_pair &a, const ritz_pair &b) {
              return a.eigenvalue > b.eigenvalue;
            });
  return pairs;
}

enum step_mode { DOMINANT_RITZ, UNSTABLE_RITZ, REWEIGHTED_GRADIENT };

struct mode_estimate {
  std::vector<vec3> direction;
  real eigenvalue = 0;
};

struct hv_result {
  std::vector<vec3> w;
  real rayleigh = 0;
};

template <typename MatVecFn>
hv_result compute_Hv(MatVecFn &&matvec, const std::vector<vec3> &v) {
  auto w = matvec(v);
  return {std::move(w), krylov_dot(v, w)};
}

template <typename MatVecFn>
mode_estimate power_iteration_step(MatVecFn &&matvec, const std::vector<vec3> &d,
                                   real sigma) {
  auto [w, rayleigh] = compute_Hv(matvec, d);
  const int N = static_cast<int>(d.size());
  if (static_cast<int>(w.size()) != N)
    w.assign(static_cast<size_t>(std::max(N, 0)), vec3::Zero());
  for (int i = 0; i < N; ++i)
    w[static_cast<size_t>(i)] -= sigma * d[static_cast<size_t>(i)];

  const real sign = (krylov_dot(w, d) >= 0) ? real(1) : real(-1);
  const real w_norm = krylov_norm(w);
  if (w_norm > real(1e-14))
    krylov_scale(sign / w_norm, w);
  return {std::move(w), rayleigh};
}

template <typename MatVecFn>
mode_estimate arnoldi_step(MatVecFn &&matvec, const std::vector<vec3> &d,
                           bool have_d, const std::vector<vec3> &G0, int k,
                           step_mode mode, int target_mode, real alpha,
                           std::function<void(const std::vector<ritz_pair> &,
                                              int)> log_fn) {
  const int N = static_cast<int>(G0.size());
  const auto &q0 = have_d ? d : G0;
  auto pairs = arnoldi(matvec, q0, k);

  int active_idx = -1;
  const int np = static_cast<int>(pairs.size());

  if (have_d && !pairs.empty()) {
    real best = 0.0;
    for (int i = 0; i < np; ++i) {
      const real ov = std::abs(krylov_dot(d, pairs[static_cast<size_t>(i)].vector));
      if (ov > best) {
        best = ov;
        active_idx = i;
      }
    }
  } else if (np == 0) {
    active_idx = -1;
  } else if (target_mode >= 0 && target_mode < np) {
    active_idx = target_mode;
  } else {
    switch (mode) {
    case DOMINANT_RITZ: {
      active_idx = 0;
      real max_abs = 0.0;
      for (int i = 0; i < np; ++i) {
        if (std::abs(pairs[static_cast<size_t>(i)].eigenvalue) > max_abs) {
          max_abs = std::abs(pairs[static_cast<size_t>(i)].eigenvalue);
          active_idx = i;
        }
      }
      break;
    }
    case UNSTABLE_RITZ: {
      active_idx = np - 1;
      for (int i = 0; i < np; ++i) {
        if (pairs[static_cast<size_t>(i)].eigenvalue < 0) {
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

  if (log_fn)
    log_fn(pairs, active_idx);

  std::vector<vec3> dir;
  if (mode == REWEIGHTED_GRADIENT) {
    dir.assign(static_cast<size_t>(std::max(N, 0)), vec3::Zero());
    if (!pairs.empty()) {
      for (int i = 0; i < np; ++i) {
        const real coeff = krylov_dot(G0, pairs[static_cast<size_t>(i)].vector);
        const real weight =
            std::pow(std::abs(pairs[static_cast<size_t>(i)].eigenvalue) +
                         real(1e-10),
                     alpha);
        krylov_axpy(coeff * weight, pairs[static_cast<size_t>(i)].vector, dir);
      }
    } else {
      dir = G0;
    }
  } else if (active_idx >= 0 && active_idx < np) {
    dir = pairs[static_cast<size_t>(active_idx)].vector;
  } else {
    dir = G0;
  }

  const real dn = krylov_norm(dir);
  if (dn > real(1e-14))
    krylov_scale(real(1) / dn, dir);

  const real lambda = (active_idx >= 0 && active_idx < np)
                          ? pairs[static_cast<size_t>(active_idx)].eigenvalue
                          : real(0);
  return {std::move(dir), lambda};
}

} // namespace duchamp
} // namespace gaudi

#endif
