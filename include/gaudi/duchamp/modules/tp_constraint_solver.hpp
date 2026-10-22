#ifndef GAUDI_DUCHAMP_TP_CONSTRAINT_SOLVER_HPP
#define GAUDI_DUCHAMP_TP_CONSTRAINT_SOLVER_HPP

#include <Eigen/QR>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>

#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/common.h"
#include "gaudi/duchamp/demo_log.hpp"
#include "gaudi/duchamp/modules/tp_config.hpp"
#include "gaudi/duchamp/modules/tp_sobolev_jacobi.hpp"

namespace gaudi {
namespace duchamp {

/// One entry of a constraint row: vertex index and the 3-vector coefficient.
struct sobolev_constraint_term {
  int i = -1;
  vec3 c = vec3::Zero();
};

/// One scalar row of C. `slack` is 1/w (0 = hard) for rows that stay
/// out of the projector. `null_ridge` is the diagonal added to BBᵀ for a
/// row that stays in the projector. 0 there is an exact nullspace row.
struct sobolev_constraint_row {
  std::vector<sobolev_constraint_term> terms;
  real slack = 0.0;
  real null_ridge = 0.0;
};

/// Object the saddle calls, in the same spirit as a Hepworth constraint.
/// `residual` and `jacobian` append the same number of rows, in the same order.
class sobolev_constraint {
public:
  using ptr = std::shared_ptr<sobolev_constraint>;
  virtual ~sobolev_constraint() = default;
  virtual void init_rest(const asawa::rod::rod &) {}
  /// Projective scale. The Gram stores w · row and w · residual.
  virtual real weight() const { return 0; }
  virtual void residual(const asawa::rod::rod &rod,
                        std::vector<real> &phi) const = 0;
  virtual void jacobian(const asawa::rod::rod &rod,
                        std::vector<sobolev_constraint_row> &rows) const = 0;
};

/// Hard edge length. Rest lives on the rod (rest_l) and follows split/collapse.
/// init_rest writes it once. Each call rebuilds rows from the current edges.
class sobolev_edge_length_constraint : public sobolev_constraint {
public:
  explicit sobolev_edge_length_constraint(real slack = 0)
      : slack_(std::max(slack, real(0))) {}

  static ptr create(real slack = 0) {
    return std::make_shared<sobolev_edge_length_constraint>(slack);
  }

  void init_rest(const asawa::rod::rod &rod) override {
    asawa::rod::rod &R = const_cast<asawa::rod::rod &>(rod);
    auto &L = R.rest_l();
    if (L.size() < R.corner_count())
      L.resize(R.corner_count(), real(0));
    const auto E = sobolev_detail::build_edges(R);
    for (size_t e = 0; e < E.n; ++e) {
      const int i0 = E.ep0[e];
      if (i0 < 0 || static_cast<size_t>(i0) >= L.size())
        continue;
      L[static_cast<size_t>(i0)] = E.ell[e];
    }
    have_rest_ = true;
  }

  void residual(const asawa::rod::rod &rod,
                std::vector<real> &phi) const override {
    if (!have_rest_)
      return;
    asawa::rod::rod &R = const_cast<asawa::rod::rod &>(rod);
    const auto &L = R.rest_l();
    const auto E = sobolev_detail::build_edges(R);
    for (size_t e = 0; e < E.n; ++e) {
      const int i0 = E.ep0[e];
      if (i0 < 0 || static_cast<size_t>(i0) >= L.size())
        continue;
      phi.push_back(E.ell[e] - L[static_cast<size_t>(i0)]);
    }
  }

  void jacobian(const asawa::rod::rod &rod,
                std::vector<sobolev_constraint_row> &rows) const override {
    if (!have_rest_)
      return;
    asawa::rod::rod &R = const_cast<asawa::rod::rod &>(rod);
    const auto &L = R.rest_l();
    const auto E = sobolev_detail::build_edges(R);
    for (size_t e = 0; e < E.n; ++e) {
      const int i0 = E.ep0[e];
      if (i0 < 0 || static_cast<size_t>(i0) >= L.size())
        continue;
      sobolev_constraint_row row;
      row.null_ridge = slack_;
      row.terms.push_back({i0, -E.T[e]});
      row.terms.push_back({E.ep1[e], E.T[e]});
      rows.push_back(std::move(row));
    }
  }

private:
  bool have_rest_ = false;
  real slack_ = 0.0;
};

/// Barycenter. slack 0 is hard. The mean position is captured once in init_rest.
class sobolev_barycenter_constraint : public sobolev_constraint {
public:
  explicit sobolev_barycenter_constraint(real slack = 0)
      : slack_(std::max(slack, real(0))) {}

  static ptr create(real slack = 0) {
    return std::make_shared<sobolev_barycenter_constraint>(slack);
  }

  void init_rest(const asawa::rod::rod &rod) override {
    rest_ = mean_of(rod);
    have_rest_ = true;
  }

  void residual(const asawa::rod::rod &rod,
                std::vector<real> &phi) const override {
    if (!have_rest_ || rod.x().empty())
      return;
    const vec3 d = mean_of(rod) - rest_;
    phi.push_back(d[0]);
    phi.push_back(d[1]);
    phi.push_back(d[2]);
  }

  void jacobian(const asawa::rod::rod &rod,
                std::vector<sobolev_constraint_row> &rows) const override {
    const auto &x = rod.x();
    if (!have_rest_ || x.empty())
      return;
    const real w = real(1) / real(x.size());
    for (int ax = 0; ax < 3; ++ax) {
      sobolev_constraint_row row;
      row.slack = slack_;
      vec3 e = vec3::Zero();
      e[ax] = w;
      row.terms.reserve(x.size());
      for (size_t i = 0; i < x.size(); ++i)
        row.terms.push_back({static_cast<int>(i), e});
      rows.push_back(std::move(row));
    }
  }

private:
  static vec3 mean_of(const asawa::rod::rod &rod) {
    const auto &x = rod.x();
    vec3 m = vec3::Zero();
    if (x.empty())
      return m;
    for (const vec3 &p : x)
      m += p;
    return m / real(x.size());
  }

  vec3 rest_ = vec3::Zero();
  bool have_rest_ = false;
  real slack_ = 0.0;
};

/// Elastic edge length. Row is ±T, residual is ℓ − ℓ_rest. The scale is w.
/// This is a projective row: the saddle adds w² CᵀC, and the load adds
/// w² (ℓ − ℓ_rest) Cᵀ.
class sobolev_elastic_constraint : public sobolev_constraint {
public:
  explicit sobolev_elastic_constraint(real w) : w_(std::max(w, real(0))) {}

  static ptr create(real w) {
    return std::make_shared<sobolev_elastic_constraint>(w);
  }

  real weight() const override { return w_; }

  void init_rest(const asawa::rod::rod &rod) override {
    asawa::rod::rod &R = const_cast<asawa::rod::rod &>(rod);
    auto &L = R.rest_l();
    if (L.size() < R.corner_count())
      L.resize(R.corner_count(), real(0));
    const auto E = sobolev_detail::build_edges(R);
    for (size_t e = 0; e < E.n; ++e) {
      const int i0 = E.ep0[e];
      if (i0 < 0 || static_cast<size_t>(i0) >= L.size())
        continue;
      L[static_cast<size_t>(i0)] = E.ell[e];
    }
    have_rest_ = true;
  }

  void residual(const asawa::rod::rod &rod,
                std::vector<real> &phi) const override {
    if (!(w_ > real(0)) || !have_rest_)
      return;
    asawa::rod::rod &R = const_cast<asawa::rod::rod &>(rod);
    const auto &L = R.rest_l();
    const auto E = sobolev_detail::build_edges(R);
    for (size_t e = 0; e < E.n; ++e) {
      const int i0 = E.ep0[e];
      if (i0 < 0 || static_cast<size_t>(i0) >= L.size())
        continue;
      phi.push_back(E.ell[e] - L[static_cast<size_t>(i0)]);
    }
  }

  void jacobian(const asawa::rod::rod &rod,
                std::vector<sobolev_constraint_row> &rows) const override {
    if (!(w_ > real(0)))
      return;
    asawa::rod::rod &R = const_cast<asawa::rod::rod &>(rod);
    const auto E = sobolev_detail::build_edges(R);
    for (size_t e = 0; e < E.n; ++e) {
      sobolev_constraint_row row;
      row.terms.push_back({E.ep0[e], -E.T[e]});
      row.terms.push_back({E.ep1[e], E.T[e]});
      rows.push_back(std::move(row));
    }
  }

private:
  bool have_rest_ = false;
  real w_ = 0;
};

/// Weighted C1. Residual is the second difference x_prev - 2 x + x_next.
/// Slack on each scalar row is 1/w. w <= 0 adds no rows.
class sobolev_c1_constraint : public sobolev_constraint {
public:
  explicit sobolev_c1_constraint(real w) : w_(w) {}

  static ptr create(real w) {
    return std::make_shared<sobolev_c1_constraint>(w);
  }

  real weight() const override { return w_; }

  void residual(const asawa::rod::rod &rod,
                std::vector<real> &phi) const override {
    if (!(w_ > real(0)))
      return;
    const auto &x = rod.x();
    const size_t n = x.size();
    for (size_t i = 0; i < n; ++i) {
      const auto ci = asawa::rod::corner_id(static_cast<int>(i));
      const auto cp = rod.prev(ci);
      const auto cn = rod.next(ci);
      if (cp < asawa::rod::corner_id(0) || cn < asawa::rod::corner_id(0))
        continue;
      const size_t ip = static_cast<size_t>(static_cast<int>(cp));
      const size_t in = static_cast<size_t>(static_cast<int>(cn));
      if (ip >= n || in >= n)
        continue;
      const vec3 d = x[ip] - real(2) * x[i] + x[in];
      phi.push_back(d[0]);
      phi.push_back(d[1]);
      phi.push_back(d[2]);
    }
  }

  void jacobian(const asawa::rod::rod &rod,
                std::vector<sobolev_constraint_row> &rows) const override {
    if (!(w_ > real(0)))
      return;
    const auto &x = rod.x();
    const size_t n = x.size();
    const real slack = real(1) / w_;
    for (size_t i = 0; i < n; ++i) {
      const auto ci = asawa::rod::corner_id(static_cast<int>(i));
      const auto cp = rod.prev(ci);
      const auto cn = rod.next(ci);
      if (cp < asawa::rod::corner_id(0) || cn < asawa::rod::corner_id(0))
        continue;
      const int ip = static_cast<int>(cp);
      const int in = static_cast<int>(cn);
      if (ip < 0 || in < 0 || static_cast<size_t>(ip) >= n ||
          static_cast<size_t>(in) >= n)
        continue;
      for (int ax = 0; ax < 3; ++ax) {
        vec3 ep = vec3::Zero();
        vec3 ec = vec3::Zero();
        vec3 en = vec3::Zero();
        ep[ax] = real(1);
        ec[ax] = real(-2);
        en[ax] = real(1);
        sobolev_constraint_row row;
        row.slack = slack;
        row.terms.push_back({ip, ep});
        row.terms.push_back({static_cast<int>(i), ec});
        row.terms.push_back({in, en});
        rows.push_back(std::move(row));
      }
    }
  }

private:
  real w_ = 1.0;
};

/// Weighted C2. Residual is the third difference
/// -x_{i-2} + 2 x_{i-1} - 2 x_{i+1} + x_{i+2} (center coefficient 0).
/// Slack on each scalar row is 1/w. w <= 0 adds no rows.
class sobolev_c2_constraint : public sobolev_constraint {
public:
  explicit sobolev_c2_constraint(real w) : w_(w) {}

  static ptr create(real w) {
    return std::make_shared<sobolev_c2_constraint>(w);
  }

  real weight() const override { return w_; }

  void residual(const asawa::rod::rod &rod,
                std::vector<real> &phi) const override {
    if (!(w_ > real(0)))
      return;
    const auto &x = rod.x();
    int im2 = 0, im1 = 0, ip1 = 0, ip2 = 0;
    for (size_t i = 0; i < x.size(); ++i) {
      if (!stencil(rod, static_cast<int>(i), im2, im1, ip1, ip2))
        continue;
      const vec3 d = -x[static_cast<size_t>(im2)] +
                     real(2) * x[static_cast<size_t>(im1)] -
                     real(2) * x[static_cast<size_t>(ip1)] +
                     x[static_cast<size_t>(ip2)];
      phi.push_back(d[0]);
      phi.push_back(d[1]);
      phi.push_back(d[2]);
    }
  }

  void jacobian(const asawa::rod::rod &rod,
                std::vector<sobolev_constraint_row> &rows) const override {
    if (!(w_ > real(0)))
      return;
    const real slack = real(1) / w_;
    int im2 = 0, im1 = 0, ip1 = 0, ip2 = 0;
    for (size_t i = 0; i < rod.x().size(); ++i) {
      if (!stencil(rod, static_cast<int>(i), im2, im1, ip1, ip2))
        continue;
      const int id[4] = {im2, im1, ip1, ip2};
      const real coeff[4] = {real(-1), real(2), real(-2), real(1)};
      for (int ax = 0; ax < 3; ++ax) {
        sobolev_constraint_row row;
        row.slack = slack;
        for (int k = 0; k < 4; ++k) {
          vec3 e = vec3::Zero();
          e[ax] = coeff[k];
          row.terms.push_back({id[k], e});
        }
        rows.push_back(std::move(row));
      }
    }
  }

private:
  static bool stencil(const asawa::rod::rod &rod, int i, int &im2, int &im1,
                      int &ip1, int &ip2) {
    const size_t n = rod.x().size();
    const auto ci = asawa::rod::corner_id(i);
    const auto cim1 = rod.prev(ci);
    const auto cip1 = rod.next(ci);
    if (cim1 < asawa::rod::corner_id(0) || cip1 < asawa::rod::corner_id(0))
      return false;
    const auto cim2 = rod.prev(cim1);
    const auto cip2 = rod.next(cip1);
    if (cim2 < asawa::rod::corner_id(0) || cip2 < asawa::rod::corner_id(0))
      return false;
    im2 = static_cast<int>(cim2);
    im1 = static_cast<int>(cim1);
    ip1 = static_cast<int>(cip1);
    ip2 = static_cast<int>(cip2);
    if (im2 < 0 || im1 < 0 || ip1 < 0 || ip2 < 0)
      return false;
    if (static_cast<size_t>(im2) >= n || static_cast<size_t>(im1) >= n ||
        static_cast<size_t>(ip1) >= n || static_cast<size_t>(ip2) >= n)
      return false;
    return true;
  }

  real w_ = 1.0;
};

/// Wider second difference, four edges across: x_{i-2} - 2 x_i + x_{i+2}.
/// The residual points at the chord of that window. w <= 0 adds no rows.
class sobolev_c4_constraint : public sobolev_constraint {
public:
  explicit sobolev_c4_constraint(real w) : w_(w) {}

  static ptr create(real w) {
    return std::make_shared<sobolev_c4_constraint>(w);
  }

  real weight() const override { return w_; }

  void residual(const asawa::rod::rod &rod,
                std::vector<real> &phi) const override {
    if (!(w_ > real(0)))
      return;
    const auto &x = rod.x();
    int im2 = 0, ip2 = 0;
    for (size_t i = 0; i < x.size(); ++i) {
      if (!ends(rod, static_cast<int>(i), im2, ip2))
        continue;
      const vec3 d = x[static_cast<size_t>(im2)] - real(2) * x[i] +
                     x[static_cast<size_t>(ip2)];
      phi.push_back(d[0]);
      phi.push_back(d[1]);
      phi.push_back(d[2]);
    }
  }

  void jacobian(const asawa::rod::rod &rod,
                std::vector<sobolev_constraint_row> &rows) const override {
    if (!(w_ > real(0)))
      return;
    const real slack = real(1) / w_;
    int im2 = 0, ip2 = 0;
    for (size_t i = 0; i < rod.x().size(); ++i) {
      if (!ends(rod, static_cast<int>(i), im2, ip2))
        continue;
      for (int ax = 0; ax < 3; ++ax) {
        vec3 em = vec3::Zero();
        vec3 ec = vec3::Zero();
        vec3 ep = vec3::Zero();
        em[ax] = real(1);
        ec[ax] = real(-2);
        ep[ax] = real(1);
        sobolev_constraint_row row;
        row.slack = slack;
        row.terms.push_back({im2, em});
        row.terms.push_back({static_cast<int>(i), ec});
        row.terms.push_back({ip2, ep});
        rows.push_back(std::move(row));
      }
    }
  }

private:
  static bool ends(const asawa::rod::rod &rod, int i, int &im2, int &ip2) {
    const size_t n = rod.x().size();
    const auto ci = asawa::rod::corner_id(i);
    const auto cim1 = rod.prev(ci);
    const auto cip1 = rod.next(ci);
    if (cim1 < asawa::rod::corner_id(0) || cip1 < asawa::rod::corner_id(0))
      return false;
    const auto cim2 = rod.prev(cim1);
    const auto cip2 = rod.next(cip1);
    if (cim2 < asawa::rod::corner_id(0) || cip2 < asawa::rod::corner_id(0))
      return false;
    im2 = static_cast<int>(cim2);
    ip2 = static_cast<int>(cip2);
    if (im2 < 0 || ip2 < 0)
      return false;
    if (static_cast<size_t>(im2) >= n || static_cast<size_t>(ip2) >= n)
      return false;
    return true;
  }

  real w_ = 1.0;
};

/// Stacks constraint objects and solves the saddle with the fractional metric.
/// Repulsive Curves §5.3. The V-cycle is A^{-1}. The solver does not name
/// the constraints it holds.
class sobolev_constraint_solver {
public:
  int max_cg = 8;

  void add(sobolev_constraint::ptr c) {
    if (c)
      constraints_.push_back(std::move(c));
  }

  void add_proj(sobolev_constraint::ptr c) {
    if (c)
      proj_.push_back(std::move(c));
  }

  bool empty() const { return constraints_.empty() && proj_.empty(); }

  void init_rest(const asawa::rod::rod &rod) {
    for (const auto &c : constraints_)
      c->init_rest(rod);
    for (const auto &c : proj_)
      c->init_rest(rod);
  }

  /// g solves (A + C_pᵀ W C_p) g + C_hᵀ λ = f + C_pᵀ W Φ.
  std::vector<vec3>
  project_gradient(asawa::rod::rod &rod, const std::vector<vec3> &f,
                   const sobolev_tangent_point_force_config &cfg,
                   real l0_abs) const {
    if (f.empty())
      return {};
    if (empty())
      return sobolev_mg_solve(rod, f, cfg, l0_abs);

    std::vector<vec3> load = f;
    sobolev_gram Gg;
    fill_gram(rod, Gg);
    Gg.add_force(load);

    std::vector<sobolev_constraint_row> rows;
    stack_jacobian(rod, rows);
    if (rows.empty())
      return apply_A_proj(rod, load, cfg, l0_abs);

    std::vector<sobolev_constraint_row> soft;
    for (const auto &row : rows) {
      if (row.slack != real(0))
        soft.push_back(row);
    }
    const std::vector<vec3> g0 = apply_A_proj(rod, load, cfg, l0_abs);
    if (soft.empty())
      return g0;

    std::vector<real> rhs;
    apply_C(soft, g0, rhs);
    std::vector<real> lambda(soft.size(), real(0));
    schur_proj(rod, soft, rhs, cfg, l0_abs, lambda);

    std::vector<vec3> ct;
    apply_CT(soft, lambda, f.size(), ct);
    const std::vector<vec3> corr = apply_A_proj(rod, ct, cfg, l0_abs);
    std::vector<vec3> g = g0;
    for (size_t i = 0; i < g.size() && i < corr.size(); ++i)
      g[i] -= corr[i];
    return g;
  }

  /// Paper §5.3.2: repeat the Newton correction until max |Φ| < 10^{-4}.
  /// Returns that max. Eight steps is the cap; the paper rarely needs three.
  static constexpr real backproj_tol = real(1e-4);
  static constexpr int backproj_newton = 8;

  real project_positions(asawa::rod::rod &rod,
                         const sobolev_tangent_point_force_config &cfg,
                         real l0_abs, bool log = false) const {
    if (empty() || rod.x().empty())
      return real(0);
    real viol = real(0);
    int used = 0;
    for (int it = 0; it < backproj_newton; ++it) {
      std::vector<sobolev_constraint_row> rows;
      std::vector<real> phi;
      stack(rod, rows, phi);
      viol = phi_abs_max(phi);
      if (rows.empty() || phi.size() != rows.size() || !(viol >= backproj_tol))
        break;
      const real before = viol;
      const std::vector<vec3> x_pre = rod.x();
      if (!newton_backproject_once(rod, cfg, l0_abs, rows, phi, log && used == 0))
        break;
      std::vector<real> after_phi;
      std::vector<sobolev_constraint_row> after_rows;
      stack(rod, after_rows, after_phi);
      const real after = phi_abs_max(after_phi);
      if (!(after < before)) {
        rod.x() = x_pre;
        rod._update_frames();
        viol = before;
        break;
      }
      ++used;
      viol = after;
    }
    if (log) {
      demo_log() << "[sobolev proj] newton=" << used
                << " |Phi|_max=" << viol << std::endl;
    }
    return viol;
  }

private:
  static real phi_abs_max(const std::vector<real> &phi) {
    real m = 0.0;
    for (real a : phi) {
      const real b = std::abs(a);
      if (std::isfinite(b))
        m = std::max(m, b);
    }
    return m;
  }

  /// One step of Eq. 23. z = C†(−Φ) is a particular solution. y solves
  /// P A P y = P A z, so x += z − y is the A-smallest correction.
  bool newton_backproject_once(asawa::rod::rod &rod,
                               const sobolev_tangent_point_force_config &cfg,
                               real l0_abs,
                               const std::vector<sobolev_constraint_row> &rows,
                               const std::vector<real> &phi, bool log) const {
    const int m = static_cast<int>(rows.size());
    if (m == 0 || static_cast<int>(phi.size()) != m)
      return false;
    matX G = matX::Zero(m, m);
    const int nvert = static_cast<int>(rod.x().size());
    std::vector<std::vector<std::pair<int, vec3>>> at(
        static_cast<size_t>(std::max(nvert, 0)));
    for (int a = 0; a < m; ++a) {
      for (const auto &t : rows[static_cast<size_t>(a)].terms) {
        if (t.i < 0 || t.i >= nvert)
          continue;
        at[static_cast<size_t>(t.i)].push_back({a, t.c});
      }
      G(a, a) += rows[static_cast<size_t>(a)].slack;
    }
    for (const auto &touch : at)
      for (const auto &ta : touch)
        for (const auto &tb : touch)
          G(ta.first, tb.first) += ta.second.dot(tb.second);
    vecX rhs(m);
    for (int a = 0; a < m; ++a)
      rhs(a) = -phi[static_cast<size_t>(a)];
    Eigen::CompleteOrthogonalDecomposition<matX> cod(G);
    const vecX lambda = cod.solve(rhs);
    std::vector<real> lam(static_cast<size_t>(m));
    for (int a = 0; a < m; ++a)
      lam[static_cast<size_t>(a)] = lambda(a);
    std::vector<vec3> z;
    apply_CT(rows, lam, rod.x().size(), z);
    std::vector<vec3> Az;
    sobolev_detail::matvec_A(rod, cfg.alpha_B, cfg.alpha_0, cfg.p1, l0_abs,
                             cfg.bh_eps, z, Az);
    sobolev_gram Gram;
    fill_gram(rod, Gram);
    Gram.apply(z, Az);
    const std::vector<vec3> y = apply_A_proj(rod, Az, cfg, l0_abs);
    if (log)
      log_position_projection(rod, cfg, l0_abs, phi, z, Az, y);
    auto &x = rod.x();
    for (size_t i = 0; i < x.size() && i < z.size(); ++i) {
      const vec3 yi = (i < y.size()) ? y[i] : vec3::Zero();
      x[i] += z[i] - yi;
    }
    rod._update_frames();
    return true;
  }

  static real vec_l2(const std::vector<vec3> &v) {
    real s = 0.0;
    for (const vec3 &a : v)
      s += a.squaredNorm();
    return std::sqrt(s);
  }

  static real vec_peak(const std::vector<vec3> &v) {
    real p = 0.0;
    for (const vec3 &a : v) {
      const real n = a.norm();
      if (std::isfinite(n))
        p = std::max(p, n);
    }
    return p;
  }

  /// |P A P z - P A P y| against the load actually solved, before x += z - y.
  void log_position_projection(asawa::rod::rod &rod,
                               const sobolev_tangent_point_force_config &cfg,
                               real l0_abs, const std::vector<real> &phi,
                               const std::vector<vec3> &z,
                               const std::vector<vec3> &Az,
                               const std::vector<vec3> &y) const {
    real phi2 = 0.0;
    for (real a : phi)
      phi2 += a * a;
    std::vector<vec3> dxy(z.size(), vec3::Zero());
    for (size_t i = 0; i < z.size(); ++i) {
      const vec3 yi = (i < y.size()) ? y[i] : vec3::Zero();
      dxy[i] = z[i] - yi;
    }

    sobolev_projector P;
    sobolev_gram Gunused;
    build_level(rod, P, Gunused,
                const_cast<sobolev_constraint_solver *>(this));
    std::vector<vec3> Paz = Az;
    P.apply(Paz);
    std::vector<vec3> Py = y;
    P.apply(Py);
    std::vector<vec3> APy;
    sobolev_detail::matvec_A(rod, cfg.alpha_B, cfg.alpha_0, cfg.p1, l0_abs,
                             cfg.bh_eps, Py, APy);
    P.apply(APy);
    real r2 = 0.0;
    const size_t n = std::min(Paz.size(), APy.size());
    for (size_t i = 0; i < n; ++i)
      r2 += (Paz[i] - APy[i]).squaredNorm();

    std::vector<vec3> Bz, B0z, By, B0y;
    sobolev_detail::matvec_B(rod, z, cfg.p1, l0_abs, cfg.bh_eps, Bz);
    sobolev_detail::matvec_B0(rod, z, cfg.p1, l0_abs, cfg.bh_eps, B0z);
    sobolev_detail::matvec_B(rod, y, cfg.p1, l0_abs, cfg.bh_eps, By);
    sobolev_detail::matvec_B0(rod, y, cfg.p1, l0_abs, cfg.bh_eps, B0y);
    const real ymax = vec_peak(y);
    demo_log() << "[sobolev proj] |Phi|=" << std::sqrt(phi2)
              << " |z|_max=" << vec_peak(z) << " |y|_max=" << ymax
              << " |z-y|_max=" << vec_peak(dxy) << " |Az|=" << vec_l2(Az)
              << " |Bz|=" << vec_l2(Bz) << " |B0z|=" << vec_l2(B0z)
              << " |PAz|=" << vec_l2(Paz) << " |By|=" << vec_l2(By)
              << " |B0y|=" << vec_l2(B0y) << " |r|=" << std::sqrt(r2)
              << std::endl;
    if (ymax > real(10))
      dump_projection_fields(rod.x(), y, z);
  }

  static void write_field(std::ostream &out, const char *name,
                          const std::vector<vec3> &v, size_t n) {
    out << "# " << name << " " << n << "\n";
    for (size_t i = 0; i < n; ++i) {
      const vec3 a = (i < v.size()) ? v[i] : vec3::Zero();
      out << a[0] << " " << a[1] << " " << a[2] << "\n";
    }
  }

  /// Pose before x += z - y, plus the fields that produced the jump.
  static void dump_projection_fields(const std::vector<vec3> &x,
                                     const std::vector<vec3> &y,
                                     const std::vector<vec3> &z) {
    std::error_code ec;
    std::filesystem::create_directories("dump", ec);
    const std::filesystem::path path = "dump/sobolev_blowup.txt";
    std::ofstream out(path);
    if (!out) {
      demo_log() << "[sobolev proj] dump failed " << path << std::endl;
      return;
    }
    const size_t n = x.size();
    write_field(out, "x", x, n);
    write_field(out, "y", y, n);
    write_field(out, "z", z, n);
    demo_log() << "[sobolev proj] dump " << std::filesystem::absolute(path)
              << " n=" << n << std::endl;
  }

  void stack_jacobian(const asawa::rod::rod &rod,
                      std::vector<sobolev_constraint_row> &rows) const {
    for (const auto &c : constraints_)
      c->jacobian(rod, rows);
  }

  void stack(const asawa::rod::rod &rod,
             std::vector<sobolev_constraint_row> &rows,
             std::vector<real> &phi) const {
    for (const auto &c : constraints_) {
      const size_t n0 = rows.size();
      const size_t p0 = phi.size();
      c->residual(rod, phi);
      c->jacobian(rod, rows);
      if (rows.size() - n0 != phi.size() - p0) {
        rows.resize(n0);
        phi.resize(p0);
      }
    }
  }

  static void apply_CT(const std::vector<sobolev_constraint_row> &rows,
                       const std::vector<real> &lambda, size_t n,
                       std::vector<vec3> &out) {
    out.assign(n, vec3::Zero());
    for (size_t k = 0; k < rows.size() && k < lambda.size(); ++k) {
      for (const auto &t : rows[k].terms) {
        if (t.i < 0 || static_cast<size_t>(t.i) >= n)
          continue;
        out[static_cast<size_t>(t.i)] += lambda[k] * t.c;
      }
    }
  }

  static void apply_C(const std::vector<sobolev_constraint_row> &rows,
                      const std::vector<vec3> &y, std::vector<real> &out) {
    out.assign(rows.size(), real(0));
    for (size_t k = 0; k < rows.size(); ++k) {
      real s = 0.0;
      for (const auto &t : rows[k].terms) {
        if (t.i < 0 || static_cast<size_t>(t.i) >= y.size())
          continue;
        s += t.c.dot(y[static_cast<size_t>(t.i)]);
      }
      out[k] = s;
    }
  }

  static void build_level(asawa::rod::rod &rod, sobolev_projector &P,
                          sobolev_gram &G, void *ctx) {
    auto *self = static_cast<sobolev_constraint_solver *>(ctx);
    self->prepare_level_rod(rod);
    std::vector<sobolev_constraint_row> rows;
    self->stack_jacobian(rod, rows);
    P.clear();
    P.ridge.reserve(rows.size());
    for (const auto &row : rows) {
      if (row.slack != real(0))
        continue;
      std::vector<sobolev_projector::term> terms;
      terms.reserve(row.terms.size());
      for (const auto &t : row.terms)
        terms.push_back({t.i, t.c});
      if (terms.empty())
        continue;
      P.rows.push_back(std::move(terms));
      P.ridge.push_back(std::max(row.null_ridge, real(0)));
    }
    P.factor();
    self->fill_gram(rod, G);
  }

  void prepare_level_rod(asawa::rod::rod &rod) const {
    bool edges = false;
    for (const auto &c : constraints_) {
      if (dynamic_cast<const sobolev_edge_length_constraint *>(c.get()))
        edges = true;
    }
    if (!edges)
      return;
    auto &L = rod.rest_l();
    if (L.size() >= rod.corner_count())
      return;
    L.assign(rod.corner_count(), real(0));
    const auto E = sobolev_detail::build_edges(rod);
    for (size_t e = 0; e < E.n; ++e) {
      const int i0 = E.ep0[e];
      if (i0 < 0 || static_cast<size_t>(i0) >= L.size())
        continue;
      L[static_cast<size_t>(i0)] = E.ell[e];
    }
  }

  void fill_gram(const asawa::rod::rod &rod, sobolev_gram &G) const {
    G.clear();
    for (const auto &c : proj_) {
      std::vector<sobolev_constraint_row> rows;
      std::vector<real> phi;
      c->residual(rod, phi);
      c->jacobian(rod, rows);
      const real w = c->weight();
      if (!(w > real(0)) || rows.empty())
        continue;
      const bool have_phi = phi.size() == rows.size();
      for (size_t k = 0; k < rows.size(); ++k) {
        sobolev_gram::row gr;
        gr.phi = have_phi ? w * phi[k] : real(0);
        gr.terms.reserve(rows[k].terms.size());
        for (const auto &t : rows[k].terms)
          gr.terms.push_back({t.i, w * t.c});
        if (!gr.terms.empty())
          G.rows.push_back(std::move(gr));
      }
    }
  }

  sobolev_mg_hook proj_hook() const {
    sobolev_mg_hook hook;
    hook.build = &sobolev_constraint_solver::build_level;
    hook.ctx = const_cast<sobolev_constraint_solver *>(this);
    return hook;
  }

  std::vector<vec3>
  apply_A_proj(asawa::rod::rod &rod, const std::vector<vec3> &b,
               const sobolev_tangent_point_force_config &cfg,
               real l0_abs) const {
    return sobolev_mg_solve(rod, b, cfg, l0_abs, {}, proj_hook());
  }

  static std::vector<vec3>
  apply_A_inv(asawa::rod::rod &rod, const std::vector<vec3> &b,
              const sobolev_tangent_point_force_config &cfg, real l0_abs) {
    return sobolev_mg_solve(rod, b, cfg, l0_abs);
  }

  /// (C A^{-1} Cᵀ + diag(slack)) λ = rhs.
  void schur_proj(asawa::rod::rod &rod,
                const std::vector<sobolev_constraint_row> &rows,
                const std::vector<real> &rhs,
                const sobolev_tangent_point_force_config &cfg, real l0_abs,
                std::vector<real> &lambda) const {
    const size_t m = rows.size();
    if (m == 0 || rhs.size() != m)
      return;
    lambda.assign(m, real(0));

    auto apply_S = [&](const std::vector<real> &v, std::vector<real> &out) {
      std::vector<vec3> ct;
      apply_CT(rows, v, rod.x().size(), ct);
      const std::vector<vec3> y = apply_A_proj(rod, ct, cfg, l0_abs);
      apply_C(rows, y, out);
      for (size_t k = 0; k < m; ++k)
        out[k] += rows[k].slack * v[k];
    };

    std::vector<real> r = rhs;
    std::vector<real> p = r;
    real rr = 0.0;
    for (real a : r)
      rr += a * a;
    const real rr0 = rr;
    const real tol = real(1e-10) * std::max(rr0, real(1));
    const int niter = std::max(1, max_cg);
    for (int it = 0; it < niter; ++it) {
      if (!(rr > tol))
        break;
      std::vector<real> Sp;
      apply_S(p, Sp);
      real pSp = 0.0;
      for (size_t k = 0; k < m; ++k)
        pSp += p[k] * Sp[k];
      if (!(std::abs(pSp) > real(1e-30)))
        break;
      const real alpha = rr / pSp;
      real rr_new = 0.0;
      for (size_t k = 0; k < m; ++k) {
        lambda[k] += alpha * p[k];
        r[k] -= alpha * Sp[k];
        rr_new += r[k] * r[k];
      }
      if (!(rr_new > tol))
        break;
      const real beta = rr_new / rr;
      for (size_t k = 0; k < m; ++k)
        p[k] = r[k] + beta * p[k];
      rr = rr_new;
    }
  }

  std::vector<sobolev_constraint::ptr> constraints_;
  std::vector<sobolev_constraint::ptr> proj_;
};

inline sobolev_constraint_solver
make_sobolev_constraint_solver(const tangent_point_solver_config &cfg) {
  sobolev_constraint_solver solver;
  solver.max_cg = std::max(1, cfg.constraint_cg);
  auto add_hard = [&](const sobolev_constraint_spec &spec) {
    if (spec.w == real(0))
      return;
    switch (spec.type) {
    case sobolev_constraint_type::edge_length:
      solver.add(sobolev_edge_length_constraint::create(spec.slack));
      break;
    case sobolev_constraint_type::barycenter:
      solver.add(sobolev_barycenter_constraint::create(spec.slack));
      break;
    case sobolev_constraint_type::c1:
      solver.add(sobolev_c1_constraint::create(spec.w));
      break;
    case sobolev_constraint_type::c2:
      solver.add(sobolev_c2_constraint::create(spec.w));
      break;
    case sobolev_constraint_type::c4:
      solver.add(sobolev_c4_constraint::create(spec.w));
      break;
    case sobolev_constraint_type::elastic:
      break;
    }
  };
  auto add_proj = [&](const sobolev_constraint_spec &spec) {
    if (spec.w == real(0))
      return;
    switch (spec.type) {
    case sobolev_constraint_type::c1:
      solver.add_proj(sobolev_c1_constraint::create(spec.w));
      break;
    case sobolev_constraint_type::c2:
      solver.add_proj(sobolev_c2_constraint::create(spec.w));
      break;
    case sobolev_constraint_type::c4:
      solver.add_proj(sobolev_c4_constraint::create(spec.w));
      break;
    case sobolev_constraint_type::elastic:
    case sobolev_constraint_type::edge_length:
      solver.add_proj(sobolev_elastic_constraint::create(spec.w));
      break;
    case sobolev_constraint_type::barycenter:
      break;
    }
  };
  for (const sobolev_constraint_spec &spec : cfg.hard_constraints)
    add_hard(spec);
  for (const sobolev_constraint_spec &spec : cfg.proj_constraints)
    add_proj(spec);
  return solver;
}

} // namespace duchamp
} // namespace gaudi

#endif
