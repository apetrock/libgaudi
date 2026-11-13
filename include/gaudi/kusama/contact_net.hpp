#ifndef GAUDI_KUSAMA_CONTACT_NET_HPP
#define GAUDI_KUSAMA_CONTACT_NET_HPP

/// The principal contact element net as a solver (scripts/sympy/cr_poisson/
/// MATH.md 19.10 to 19.12; prototype `contact_net.py`).
///
/// State: one contact element per vertex, positions `x` and normals `n` both
/// unknown. Each edge holds one sphere, the common sphere of its endpoint
/// contact elements, through the discrete Rodrigues relation with twist,
///
///     dn = -k dx_t - tau (nbar x dx),   nbar = normalize(n_i + n_j),
///
/// `dx_t` the part of `dx` normal to `nbar`. `k` is the curvature along the
/// edge, `tau` the twist (0 on curvature-line edges); with unit normals the
/// two are the whole of `dn`. Per directed spoke the normal offset
/// `c_ij = n_i . d / |d|`. Targets `(k, tau, c)` are pose free.
///
/// Two solvers on the same rows:
/// - `gn_solver`: Levenberg-Marquardt on the residuals (Rodrigues, scalar
///   tangency both ways, unit length, a length memory that fades with the
///   merit), generated Jacobians (`contact_net_generated.hpp`), the factor
///   kept `gn_refresh` iterates. 5 to 15 iterates to 1e-9 on the shared suite.
/// - `pd_solver`: projective dynamics on one constant factor per coordinate
///   (cross terms lagged into the right-hand side, normals fitted to their
///   rings), with Anderson (Peng et al. 2018) or Chebyshev (Wang 2015)
///   acceleration. A reference `k_ref` in the factor lets targets change
///   without refactoring.

#include "gaudi/common.h"
#include "gaudi/kusama/stencil_table.hpp"

#include "contact_net_generated.hpp"

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include <algorithm>
#include <array>
#include <cmath>
#include <deque>
#include <limits>
#include <vector>

namespace gaudi {
namespace kusama {
namespace contact_net {

using trip = Eigen::Triplet<real>;
using matS = Eigen::SparseMatrix<real>;

/// Edges `(i < j)` by compact edge id, and directed ring spokes `(i, j)`
/// grouped by `i`.
struct topology {
  int n = 0;
  std::vector<std::array<int, 2>> edges;
  std::vector<std::array<int, 2>> spokes;
  std::vector<int> spoke_begin; // n + 1
  std::vector<std::array<int, 3>> tris;
};

inline topology build_topology(const stencil_table &t) {
  topology tp;
  tp.n = t.n_verts;
  tp.edges.assign(static_cast<size_t>(t.n_edges), {-1, -1});
  for (int f = 0; f < t.n_faces; ++f) {
    const auto &id = t.ids[static_cast<size_t>(f)];
    tp.tris.push_back({id[0], id[1], id[2]});
    for (int i = 0; i < 3; ++i) {
      const int e = t.face_edges[static_cast<size_t>(f)][static_cast<size_t>(i)];
      if (e < 0)
        continue;
      const int a = id[static_cast<size_t>(i)], b = id[static_cast<size_t>((i + 1) % 3)];
      tp.edges[static_cast<size_t>(e)] = {std::min(a, b), std::max(a, b)};
    }
  }
  std::vector<std::vector<int>> nbr(static_cast<size_t>(tp.n));
  for (const auto &e : tp.edges) {
    if (e[0] < 0)
      continue;
    nbr[static_cast<size_t>(e[0])].push_back(e[1]);
    nbr[static_cast<size_t>(e[1])].push_back(e[0]);
  }
  tp.spoke_begin.assign(static_cast<size_t>(tp.n) + 1, 0);
  for (int i = 0; i < tp.n; ++i) {
    std::sort(nbr[static_cast<size_t>(i)].begin(), nbr[static_cast<size_t>(i)].end());
    tp.spoke_begin[static_cast<size_t>(i)] = static_cast<int>(tp.spokes.size());
    for (int j : nbr[static_cast<size_t>(i)])
      tp.spokes.push_back({i, j});
  }
  tp.spoke_begin[static_cast<size_t>(tp.n)] = static_cast<int>(tp.spokes.size());
  return tp;
}

/// Area-weighted vertex normals (unit).
inline std::vector<vec3> vertex_normals(const topology &tp, const std::vector<vec3> &x) {
  std::vector<vec3> N(static_cast<size_t>(tp.n), vec3::Zero());
  for (const auto &t : tp.tris) {
    const vec3 c = (x[static_cast<size_t>(t[1])] - x[static_cast<size_t>(t[0])])
                       .cross(x[static_cast<size_t>(t[2])] - x[static_cast<size_t>(t[0])]);
    for (int k = 0; k < 3; ++k)
      N[static_cast<size_t>(t[static_cast<size_t>(k)])] += c;
  }
  for (vec3 &v : N)
    if (v.norm() > real(0))
      v.normalize();
  return N;
}

struct targets {
  std::vector<real> k, tau; // per edge
  std::vector<real> c;      // per spoke
};

/// `(k, tau)` per edge and `c` per spoke of the contact elements `(x, n)`.
inline targets read(const topology &tp, const std::vector<vec3> &x, const std::vector<vec3> &n) {
  targets T;
  const size_t nE = tp.edges.size();
  T.k.assign(nE, 0);
  T.tau.assign(nE, 0);
  for (size_t e = 0; e < nE; ++e) {
    const int i = tp.edges[e][0], j = tp.edges[e][1];
    if (i < 0)
      continue;
    const vec3 ni = n[static_cast<size_t>(i)].normalized(), nj = n[static_cast<size_t>(j)].normalized();
    const vec3 dx = x[static_cast<size_t>(j)] - x[static_cast<size_t>(i)], dn = nj - ni;
    vec3 nb = ni + nj;
    nb.normalize();
    const vec3 dxt = dx - dx.dot(nb) * nb, b = nb.cross(dx);
    const real q = std::max(dxt.squaredNorm(), real(1e-300));
    T.k[e] = -dn.dot(dxt) / q;
    T.tau[e] = -dn.dot(b) / q;
  }
  T.c.assign(tp.spokes.size(), 0);
  for (size_t s = 0; s < tp.spokes.size(); ++s) {
    const int i = tp.spokes[s][0], j = tp.spokes[s][1];
    const vec3 d = x[static_cast<size_t>(j)] - x[static_cast<size_t>(i)];
    T.c[s] = n[static_cast<size_t>(i)].normalized().dot(d) / std::max(d.norm(), real(1e-300));
  }
  return T;
}

enum class kind { gn, pd };
enum class accel { none, anderson, chebyshev, anderson_chebyshev };

struct config {
  kind solver = kind::gn;
  int iterations = 30;
  real w_rodrigues = 1, w_tangency = 1, w_unit = 1;
  /// Length memory: reseeded every iterate (a proximal damping). Gauss-Newton
  /// fades it with `merit / merit_start` (a proximal term on every step
  /// dominates the soft bending modes and makes convergence linear).
  real w_memory = 0.01;
  bool memory_fade = true;
  int gn_refresh = 5;
  real lam0 = 1e-4;
  /// Projective dynamics: normal fit pull, anchor, accelerator.
  real normal_mu = 1e-2;
  real eps = 1e-6;
  accel acceleration = accel::anderson;
  int anderson_m = 5;
  int cheb_delay = 30, cheb_window = 20;
  real tol = 1e-12;
  /// `contact_net_coupling` only: inner solves per frame (the store settles
  /// over the outer iterates; each move past `goal_tol` relative re-solves,
  /// warm started from the last net), the last result held after the budget.
  int solves_per_frame = 3;
  real goal_tol = 1e-3;
};

struct result {
  std::vector<vec3> x, n;
  std::vector<real> merit; // per iterate, [0] the start
  int iterations = 0, factorizations = 0, backsolves = 0;
};

/// RMS of the target mismatch: `k` and `tau` relative to the mean `|k*|`, `c`
/// absolute.
inline real merit(const topology &tp, const targets &T, const std::vector<vec3> &x, const std::vector<vec3> &n) {
  const targets C = read(tp, x, n);
  real kbar = 0;
  for (real v : T.k)
    kbar += std::abs(v);
  kbar = std::max(kbar / std::max<real>(real(T.k.size()), 1), real(1e-300));
  real a = 0, b = 0, c = 0;
  for (size_t e = 0; e < T.k.size(); ++e) {
    a += std::pow((C.k[e] - T.k[e]) / kbar, 2);
    b += std::pow((C.tau[e] - T.tau[e]) / kbar, 2);
  }
  for (size_t s = 0; s < T.c.size(); ++s)
    c += std::pow(C.c[s] - T.c[s], 2);
  const real nE = std::max<real>(real(T.k.size()), 1), nS = std::max<real>(real(T.c.size()), 1);
  return std::sqrt(a / nE + b / nE + c / nS);
}

inline real mean_edge(const topology &tp, const std::vector<vec3> &x) {
  real s = 0;
  int m = 0;
  for (const auto &e : tp.edges)
    if (e[0] >= 0) {
      s += (x[static_cast<size_t>(e[1])] - x[static_cast<size_t>(e[0])]).norm();
      ++m;
    }
  return m ? s / m : real(1);
}

// ------------------------------------------------------------ Gauss-Newton

/// Levenberg-Marquardt on the un-lagged rows. Unknowns `(x, n)`, `6V`.
class gn_solver {
public:
  gn_solver(const topology &tp, const config &cfg) : _tp(tp), _cfg(cfg) {}

  result solve(const std::vector<vec3> &x0, const std::vector<vec3> &n0, const targets &T) {
    const int n = _tp.n;
    result out;
    out.x = x0;
    out.n = n0;
    _h = mean_edge(_tp, x0);
    out.merit.push_back(merit(_tp, T, out.x, out.n));
    const real m0 = std::max(out.merit[0], real(1e-300));
    real lam = _cfg.lam0, scale = 1;
    int age = 0;
    bool have = false;
    Eigen::SimplicialLDLT<matS> ldlt;
    matS H;
    for (int it = 0; it < _cfg.iterations; ++it) {
      if (_cfg.memory_fade && (!have || age >= _cfg.gn_refresh))
        scale = std::min(real(1), out.merit.back() / m0);
      std::vector<real> L(_tp.edges.size(), 0);
      for (size_t e = 0; e < L.size(); ++e)
        if (_tp.edges[e][0] >= 0)
          L[e] = (out.x[static_cast<size_t>(_tp.edges[e][1])] - out.x[static_cast<size_t>(_tp.edges[e][0])]).norm();
      vecX r;
      matS J;
      assemble(out.x, out.n, T, L, scale, r, &J);
      const vecX g = J.transpose() * r;
      const real cost = real(0.5) * r.squaredNorm();
      auto refactor = [&]() {
        H = matS(J.transpose() * J);
        real dmean = 0;
        for (Eigen::Index i = 0; i < H.rows(); ++i)
          dmean += H.coeff(i, i);
        dmean /= std::max<real>(real(H.rows()), 1);
        matS D(H.rows(), H.cols());
        std::vector<trip> dt;
        for (Eigen::Index i = 0; i < H.rows(); ++i)
          dt.emplace_back(i, i, lam * H.coeff(i, i) + real(1e-12) * dmean);
        D.setFromTriplets(dt.begin(), dt.end());
        ldlt.compute(matS(H + D));
        ++out.factorizations;
        age = 0;
        have = true;
      };
      if (!have || age >= _cfg.gn_refresh)
        refactor();
      bool ok = false;
      std::vector<vec3> xn, nn;
      for (int tries = 0; tries < 10; ++tries) {
        const vecX dz = -ldlt.solve(g);
        ++out.backsolves;
        xn = out.x;
        nn = out.n;
        for (int v = 0; v < n; ++v) {
          xn[static_cast<size_t>(v)] += dz.segment<3>(3 * v);
          nn[static_cast<size_t>(v)] += dz.segment<3>(3 * n + 3 * v);
        }
        vecX rn;
        assemble(xn, nn, T, L, scale, rn, nullptr);
        if (rn.allFinite() && real(0.5) * rn.squaredNorm() < cost) {
          ok = true;
          break;
        }
        lam *= real(4);
        refactor();
      }
      if (!ok)
        break;
      lam = std::max(lam / real(3), real(1e-12));
      ++age;
      real move = 0;
      for (int v = 0; v < n; ++v)
        move = std::max(move, (xn[static_cast<size_t>(v)] - out.x[static_cast<size_t>(v)]).cwiseAbs().maxCoeff());
      out.x = std::move(xn);
      out.n = std::move(nn);
      out.merit.push_back(merit(_tp, T, out.x, out.n));
      ++out.iterations;
      if (move < _cfg.tol)
        break;
    }
    for (vec3 &v : out.n)
      if (v.norm() > real(0))
        v.normalize();
    return out;
  }

private:
  /// Residual `r` and (optionally) Jacobian over `(x, n)`.
  void assemble(const std::vector<vec3> &x, const std::vector<vec3> &n, const targets &T,
                const std::vector<real> &L, real scale, vecX &r, matS *J) const {
    const int nv = _tp.n;
    const size_t nE = _tp.edges.size(), nS = _tp.spokes.size();
    r.resize(static_cast<Eigen::Index>(3 * nE + nS + nv + nE));
    std::vector<trip> jt;
    if (J)
      jt.reserve(3 * nE * 12 + nS * 9 + 3 * nv + 6 * nE);
    const real sr = std::sqrt(_cfg.w_rodrigues), st = std::sqrt(_cfg.w_tangency), su = std::sqrt(_cfg.w_unit);
    const real sm = std::sqrt(_cfg.w_memory * scale) / _h;
    auto xc = [](int v) { return 3 * v; };
    auto nc = [nv](int v) { return 3 * nv + 3 * v; };
    Eigen::Index row = 0;
    for (size_t e = 0; e < nE; ++e, row += 3) {
      const int i = _tp.edges[e][0], j = _tp.edges[e][1];
      if (i < 0) {
        r.segment<3>(row).setZero();
        continue;
      }
      Eigen::Matrix<real, 12, 1> in;
      in << x[static_cast<size_t>(i)], x[static_cast<size_t>(j)], n[static_cast<size_t>(i)], n[static_cast<size_t>(j)];
      Eigen::Matrix<real, 2, 1> c(T.k[e], T.tau[e]);
      Eigen::Matrix<real, 3, 1> y;
      Eigen::Matrix<real, 3, 12> Jl;
      contact_generated::contact_rodrigues(c, in, y, Jl);
      r.segment<3>(row) = sr * y;
      if (J) {
        const int cols[4] = {xc(i), xc(j), nc(i), nc(j)};
        for (int a = 0; a < 3; ++a)
          for (int s = 0; s < 4; ++s)
            for (int q = 0; q < 3; ++q)
              jt.emplace_back(row + a, cols[s] + q, sr * Jl(a, 3 * s + q));
      }
    }
    for (size_t s = 0; s < nS; ++s, ++row) {
      const int i = _tp.spokes[s][0], j = _tp.spokes[s][1];
      Eigen::Matrix<real, 9, 1> in;
      in << x[static_cast<size_t>(i)], x[static_cast<size_t>(j)], n[static_cast<size_t>(i)];
      Eigen::Matrix<real, 1, 1> c(T.c[s]);
      Eigen::Matrix<real, 1, 1> y;
      Eigen::Matrix<real, 1, 9> Jl;
      contact_generated::contact_tangency(c, in, y, Jl);
      r[row] = st * y[0];
      if (J) {
        const int cols[3] = {xc(i), xc(j), nc(i)};
        for (int b = 0; b < 3; ++b)
          for (int q = 0; q < 3; ++q)
            jt.emplace_back(row, cols[b] + q, st * Jl(0, 3 * b + q));
      }
    }
    for (int v = 0; v < nv; ++v, ++row) {
      const vec3 &nv3 = n[static_cast<size_t>(v)];
      const real len = nv3.norm();
      r[row] = su * (len - real(1));
      if (J && len > real(0))
        for (int q = 0; q < 3; ++q)
          jt.emplace_back(row, nc(v) + q, su * nv3[q] / len);
    }
    for (size_t e = 0; e < nE; ++e, ++row) {
      const int i = _tp.edges[e][0], j = _tp.edges[e][1];
      if (i < 0) {
        r[row] = 0;
        continue;
      }
      const vec3 d = x[static_cast<size_t>(j)] - x[static_cast<size_t>(i)];
      const real len = d.norm();
      r[row] = sm * (len - L[e]);
      if (J && len > real(0))
        for (int q = 0; q < 3; ++q) {
          jt.emplace_back(row, xc(j) + q, sm * d[q] / len);
          jt.emplace_back(row, xc(i) + q, -sm * d[q] / len);
        }
    }
    if (J) {
      J->resize(row, 6 * nv);
      J->setFromTriplets(jt.begin(), jt.end());
    }
  }

  const topology &_tp;
  config _cfg;
  real _h = 1;
};

// ------------------------------------------------- projective dynamics

/// One constant factor per coordinate on `(x_c, n_c)` (`2V`): Rodrigues rows
/// with `k_ref` in the matrix and `(k - k_ref) dx`, `k (dx . nbar) nbar` and the
/// twist lagged; tangency spokes; unit-normal rows toward the normal fitted to
/// the ring's offsets; length memory reseeded every iterate; anchor.
class pd_solver {
public:
  pd_solver(const topology &tp, const std::vector<vec3> &x0, const std::vector<real> &k_ref, const config &cfg)
      : _tp(tp), _cfg(cfg), _kmat(k_ref) {
    const int n = _tp.n;
    _h = mean_edge(_tp, x0);
    std::vector<trip> at;
    Eigen::Index r = 0;
    const real sr = std::sqrt(_cfg.w_rodrigues), st = std::sqrt(_cfg.w_tangency) / _h,
               su = std::sqrt(_cfg.w_unit), sm = std::sqrt(_cfg.w_memory) / _h;
    for (size_t e = 0; e < _tp.edges.size(); ++e, ++r) {
      const int i = _tp.edges[e][0], j = _tp.edges[e][1];
      if (i < 0)
        continue;
      at.emplace_back(r, j, sr * _kmat[e]);
      at.emplace_back(r, i, -sr * _kmat[e]);
      at.emplace_back(r, n + j, sr);
      at.emplace_back(r, n + i, -sr);
    }
    _r_tan = r;
    for (const auto &s : _tp.spokes) {
      at.emplace_back(r, s[1], st);
      at.emplace_back(r, s[0], -st);
      ++r;
    }
    _r_unit = r;
    for (int v = 0; v < n; ++v)
      at.emplace_back(r++, n + v, su);
    _r_mem = r;
    for (const auto &e : _tp.edges) {
      if (e[0] >= 0) {
        at.emplace_back(r, e[1], sm);
        at.emplace_back(r, e[0], -sm);
      }
      ++r;
    }
    _A.resize(r, 2 * n);
    _A.setFromTriplets(at.begin(), at.end());
    matS K = matS(_A.transpose() * _A);
    real dmean = 0;
    for (int v = 0; v < n; ++v)
      dmean += K.coeff(v, v);
    dmean /= std::max(n, 1);
    _dmean = dmean;
    std::vector<trip> dt;
    for (int v = 0; v < 2 * n; ++v)
      dt.emplace_back(v, v, _cfg.eps * dmean);
    matS D(2 * n, 2 * n);
    D.setFromTriplets(dt.begin(), dt.end());
    _ldlt.compute(matS(K + D));
  }

  bool ok() const { return _ldlt.info() == Eigen::Success; }

  result solve(const std::vector<vec3> &x0, const std::vector<vec3> &n0, const targets &T) {
    const int n = _tp.n;
    result out;
    out.x = x0;
    out.n = n0;
    out.factorizations = 0;
    out.merit.push_back(merit(_tp, T, out.x, out.n));
    const Eigen::Index dim = 6 * n;
    auto pack = [&](const std::vector<vec3> &X, const std::vector<vec3> &N) {
      vecX z(dim);
      for (int v = 0; v < n; ++v) {
        z.segment<3>(3 * v) = X[static_cast<size_t>(v)];
        z.segment<3>(3 * n + 3 * v) = N[static_cast<size_t>(v)];
      }
      return z;
    };
    auto unpack = [&](const vecX &z, std::vector<vec3> &X, std::vector<vec3> &N) {
      X.resize(static_cast<size_t>(n));
      N.resize(static_cast<size_t>(n));
      for (int v = 0; v < n; ++v) {
        X[static_cast<size_t>(v)] = z.segment<3>(3 * v);
        N[static_cast<size_t>(v)] = z.segment<3>(3 * n + 3 * v);
      }
    };
    std::deque<vecX> hf, hg;
    vecX zprev;
    bool cheb_on = false;
    real rho = 0, omega = 1, mprev = out.merit[0];
    for (int it = 0; it < _cfg.iterations; ++it) {
      std::vector<vec3> Xg, Ng;
      G(out.x, out.n, T, Xg, Ng);
      out.backsolves += 3;
      const vecX z = pack(out.x, out.n), g = pack(Xg, Ng);
      std::vector<vec3> Xn = Xg, Nn = Ng;
      const accel a = _cfg.acceleration;
      const bool use_anderson = a == accel::anderson || a == accel::anderson_chebyshev;
      const bool use_cheb = a == accel::chebyshev || (a == accel::anderson_chebyshev && hf.size() < 2);
      bool took_anderson = false;
      if (use_anderson) {
        hf.push_back(g - z);
        hg.push_back(g);
        while (static_cast<int>(hf.size()) > _cfg.anderson_m + 1) {
          hf.pop_front();
          hg.pop_front();
        }
        if (hf.size() > 1) {
          const int m = static_cast<int>(hf.size()) - 1;
          Eigen::MatrixXd dF(dim, m), dG(dim, m);
          for (int i = 0; i < m; ++i) {
            dF.col(i) = hf[static_cast<size_t>(i + 1)] - hf[static_cast<size_t>(i)];
            dG.col(i) = hg[static_cast<size_t>(i + 1)] - hg[static_cast<size_t>(i)];
          }
          const Eigen::VectorXd gam = dF.colPivHouseholderQr().solve(hf.back());
          std::vector<vec3> Xa, Na;
          unpack(g - dG * gam, Xa, Na);
          bool finite = true;
          for (const vec3 &v : Xa)
            finite = finite && v.allFinite();
          if (finite && merit(_tp, T, Xa, Na) <= merit(_tp, T, Xg, Ng)) {
            Xn = std::move(Xa);
            Nn = std::move(Na);
            took_anderson = true;
          } else {
            hf.erase(hf.begin(), hf.end() - 1);
            hg.erase(hg.begin(), hg.end() - 1);
          }
        }
      }
      if (use_cheb && !took_anderson) {
        if (!cheb_on && it >= _cfg.cheb_delay && static_cast<int>(out.merit.size()) > _cfg.cheb_window) {
          const real h1 = out.merit.back(), h0 = out.merit[out.merit.size() - static_cast<size_t>(_cfg.cheb_window)];
          rho = std::clamp(std::pow(h1 / std::max(h0, real(1e-300)), real(1) / real(_cfg.cheb_window - 1)), real(0),
                           real(0.9999));
          cheb_on = true;
          omega = 1;
        }
        if (cheb_on && zprev.size() == dim) {
          omega = omega == real(1) ? real(2) / (real(2) - rho * rho) : real(4) / (real(4) - rho * rho * omega);
          unpack(omega * (g - zprev) + zprev, Xn, Nn);
        }
      }
      real mn = merit(_tp, T, Xn, Nn);
      if (cheb_on && use_cheb && !took_anderson && mn > mprev) {
        Xn = Xg;
        Nn = Ng;
        omega = 1;
        mn = merit(_tp, T, Xn, Nn);
      }
      zprev = z;
      real move = 0;
      for (int v = 0; v < n; ++v)
        move = std::max(move, (Xn[static_cast<size_t>(v)] - out.x[static_cast<size_t>(v)]).cwiseAbs().maxCoeff());
      out.x = std::move(Xn);
      out.n = std::move(Nn);
      out.merit.push_back(mn);
      mprev = mn;
      ++out.iterations;
      if (move < _cfg.tol || !std::isfinite(mn))
        break;
    }
    for (vec3 &v : out.n)
      if (v.norm() > real(0))
        v.normalize();
    return out;
  }

private:
  /// The fixed-point map: right-hand side at `(x, n)`, three backsolves.
  void G(const std::vector<vec3> &x, const std::vector<vec3> &N, const targets &T, std::vector<vec3> &Xo,
         std::vector<vec3> &No) const {
    const int n = _tp.n;
    const real sr = std::sqrt(_cfg.w_rodrigues), st = std::sqrt(_cfg.w_tangency) / _h,
               su = std::sqrt(_cfg.w_unit), sm = std::sqrt(_cfg.w_memory) / _h;
    Eigen::Matrix<real, Eigen::Dynamic, 3> rows(_A.rows(), 3);
    rows.setZero();
    for (size_t e = 0; e < _tp.edges.size(); ++e) {
      const int i = _tp.edges[e][0], j = _tp.edges[e][1];
      if (i < 0)
        continue;
      const vec3 dx = x[static_cast<size_t>(j)] - x[static_cast<size_t>(i)];
      vec3 nb = N[static_cast<size_t>(i)] + N[static_cast<size_t>(j)];
      nb.normalize();
      const vec3 rod = T.k[e] * dx.dot(nb) * nb - (T.k[e] - _kmat[e]) * dx - T.tau[e] * nb.cross(dx);
      rows.row(static_cast<Eigen::Index>(e)) = sr * rod.transpose();
    }
    for (size_t s = 0; s < _tp.spokes.size(); ++s) {
      const int i = _tp.spokes[s][0], j = _tp.spokes[s][1];
      const vec3 d = x[static_cast<size_t>(j)] - x[static_cast<size_t>(i)];
      const vec3 ni = N[static_cast<size_t>(i)].normalized();
      const vec3 ds = d - ni.dot(d) * ni + T.c[s] * d.norm() * ni;
      rows.row(_r_tan + static_cast<Eigen::Index>(s)) = st * ds.transpose();
    }
    const std::vector<vec3> nf = fitted_normals(x, N, T);
    for (int v = 0; v < n; ++v)
      rows.row(_r_unit + v) = su * nf[static_cast<size_t>(v)].transpose();
    for (size_t e = 0; e < _tp.edges.size(); ++e) {
      const int i = _tp.edges[e][0], j = _tp.edges[e][1];
      if (i < 0)
        continue;
      rows.row(_r_mem + static_cast<Eigen::Index>(e)) =
          sm * (x[static_cast<size_t>(j)] - x[static_cast<size_t>(i)]).transpose();
    }
    Eigen::Matrix<real, Eigen::Dynamic, 3> b = _A.transpose() * rows;
    for (int v = 0; v < n; ++v) {
      b.row(v) += _cfg.eps * _dmean * x[static_cast<size_t>(v)].transpose();
      b.row(n + v) += _cfg.eps * _dmean * N[static_cast<size_t>(v)].transpose();
    }
    Xo.resize(static_cast<size_t>(n));
    No.resize(static_cast<size_t>(n));
    for (int c = 0; c < 3; ++c) {
      const vecX z = _ldlt.solve(vecX(b.col(c)));
      for (int v = 0; v < n; ++v) {
        Xo[static_cast<size_t>(v)][c] = z[v];
        No[static_cast<size_t>(v)][c] = z[n + v];
      }
    }
  }

  /// Normals follow positions: per vertex `min |D n - c|^2 + mu |n - n_i|^2`,
  /// normalized.
  std::vector<vec3> fitted_normals(const std::vector<vec3> &x, const std::vector<vec3> &N, const targets &T) const {
    std::vector<vec3> out(N.size());
    for (int v = 0; v < _tp.n; ++v) {
      const vec3 cur = N[static_cast<size_t>(v)].normalized();
      mat3 A = _cfg.normal_mu * mat3::Identity();
      vec3 rhs = _cfg.normal_mu * cur;
      for (int s = _tp.spoke_begin[static_cast<size_t>(v)]; s < _tp.spoke_begin[static_cast<size_t>(v) + 1]; ++s) {
        vec3 d = x[static_cast<size_t>(_tp.spokes[static_cast<size_t>(s)][1])] - x[static_cast<size_t>(v)];
        d /= std::max(d.norm(), real(1e-300));
        A += d * d.transpose();
        rhs += T.c[static_cast<size_t>(s)] * d;
      }
      const vec3 nn = A.ldlt().solve(rhs);
      out[static_cast<size_t>(v)] = nn.norm() > real(0) ? vec3(nn.normalized()) : cur;
    }
    return out;
  }

  const topology &_tp;
  config _cfg;
  std::vector<real> _kmat;
  real _h = 1, _dmean = 1;
  matS _A;
  Eigen::Index _r_tan = 0, _r_unit = 0, _r_mem = 0;
  Eigen::SimplicialLDLT<matS> _ldlt;
};

/// Either solver per `cfg.solver` (PD factored at `x0` with `k_ref = T.k`).
inline result solve(const topology &tp, const std::vector<vec3> &x0, const std::vector<vec3> &n0, const targets &T,
                    const config &cfg) {
  if (cfg.solver == kind::pd) {
    pd_solver s(tp, x0, T.k, cfg);
    result r = s.solve(x0, n0, T);
    r.factorizations = 1;
    return r;
  }
  gn_solver s(tp, cfg);
  return s.solve(x0, n0, T);
}

} // namespace contact_net
} // namespace kusama
} // namespace gaudi

#endif
