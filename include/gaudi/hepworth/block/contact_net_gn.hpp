#ifndef GAUDI_HEPWORTH_BLOCK_CONTACT_NET_GN_HPP
#define GAUDI_HEPWORTH_BLOCK_CONTACT_NET_GN_HPP

/// The contact element net's rows inside `projection_solver`, as Gauss-Newton
/// rows like `weingarten_gn` (scripts/sympy/cr_poisson/MATH.md 19.10 to
/// 19.13; `kusama/contact_net.hpp` holds the standalone solvers). One solve
/// per frame, no inner solve: the outer factor is refreshed every
/// `gn_refresh_every` iterates.
///
/// Blocks: `blocks[0]` positions, `blocks[1]` the per-vertex normals
/// (`shell_cr_block::normal_dof`), `blocks[2]` the edge store
/// (`shell_cr_block::shape_dof`). Rows:
/// - Rodrigues, per edge (3): `k (dx - (dx . nbar) nbar) + dn + tau (nbar x dx)`,
///   `k` and `tau` the store's change of `(k11, k12)` in its edge frame added
///   to the frame-start read-out (increment form: no edit is an exact fixed
///   point), so the store is coupled natively (its columns `dr/dk = dx_t`,
///   `dr/dtau = nbar x dx`);
/// - tangency, per directed spoke (1): `n_i . d / |d| - c`, `c` the frame-start
///   offset plus `dk |d0| / 2`;
/// - unit length, per vertex (1): `|n_i| - 1`.
/// Rows scaled `w k_contact_gn A_e / |e|`: the Rodrigues and tangency
/// residuals are curvature times length, so the scaled rows are integrated
/// curvatures, the unit of `weingarten_gn` and `bend`; `k_contact_gn` is
/// calibrated against `bend` (`contact_net_gn_matches_bend_at_equal_weight`).
/// (Lie sphere rows `<s_ij, P_j>`, `<s_ij, T_j>` were tried in place of
/// Rodrigues / tangency and fail off curvature-line edges: the twist leaves
/// the neighbour free to slide on the sphere; MATH.md 19.13.)

#include "gaudi/hepworth/block/block_constraint.hpp"
#include "gaudi/kusama/contact_net.hpp"
#include "gaudi/kusama/stencil_table.hpp"
#include "gaudi/kusama/weingarten.hpp"

#include "contact_net_generated.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <vector>

namespace gaudi {
namespace hepworth {
namespace block {

/// Calibration of the rows against the cotan `bend` row at equal weights.
/// `contact_net_gn_matches_bend_at_equal_weight`: at 1 the rows' mean
/// position stiffness is 0.0047 / 0.0052 of `bend` on the two test tori, so
/// `sqrt(1 / 0.005)`.
constexpr real k_contact_gn = real(14.2);

class contact_net_gn : public block_constraint {
public:
  DEFINE_CREATE_FUNC(contact_net_gn)

  /// `c_off`: optional per-spoke tangency offsets `c - k |d| / 2` to hold
  /// (default: the frame start's).
  contact_net_gn(const kusama::stencil_table *table, const std::vector<vec3> &x0, const std::vector<vec3> &n0,
                 const std::vector<vec3> &h0, real w, std::vector<sim_block::ptr> blocks, real k = k_contact_gn,
                 const std::vector<real> *c_off = nullptr)
      : block_constraint({}, w, std::move(blocks)), _t(table), _h0(h0) {
    _tp = kusama::contact_net::build_topology(*_t);
    _T0 = kusama::contact_net::read(_tp, x0, n0);
    // Per edge `A_e / |e|`: a third of the incident faces' area over length.
    const size_t nE = _tp.edges.size();
    std::vector<real> area(nE, 0);
    for (int f = 0; f < _t->n_faces; ++f) {
      const auto &id = _t->ids[static_cast<size_t>(f)];
      const real a = real(0.5) * (x0[static_cast<size_t>(id[1])] - x0[static_cast<size_t>(id[0])])
                                     .cross(x0[static_cast<size_t>(id[2])] - x0[static_cast<size_t>(id[0])])
                                     .norm();
      for (int i = 0; i < 3; ++i) {
        const int e = _t->face_edges[static_cast<size_t>(f)][static_cast<size_t>(i)];
        if (e >= 0)
          area[static_cast<size_t>(e)] += a / real(3);
      }
    }
    _s.assign(nE, 0);
    real mean = 0;
    for (size_t e = 0; e < nE; ++e) {
      const int i = _tp.edges[e][0], j = _tp.edges[e][1];
      if (i < 0)
        continue;
      const real l = (x0[static_cast<size_t>(j)] - x0[static_cast<size_t>(i)]).norm();
      _s[e] = l > real(0) ? k * area[e] / l : real(0);
      mean += _s[e];
    }
    _s_unit = mean / std::max<real>(real(nE), 1);
    std::map<std::pair<int, int>, int> edge_of;
    for (size_t e = 0; e < nE; ++e)
      edge_of[{_tp.edges[e][0], _tp.edges[e][1]}] = static_cast<int>(e);
    const size_t nS = _tp.spokes.size();
    _spoke_edge.resize(nS);
    _spoke_len.resize(nS);
    _c_off.resize(nS);
    _kappa0.resize(nS);
    _kappa_from.assign(nS, {});
    // The one-sided normal curvature of spoke `(i, j)`: vertex `i`'s tensor
    // (`kusama::weingarten::vertex_tensor`, the weighted mean of its edge
    // tensors placed in its plane) along the spoke, `d^T S_i d`. Linear in the
    // store; its coefficients are frozen at the frame start.
    const kusama::weingarten::vertex_stars stars = kusama::weingarten::build_vertex_stars(*_t, x0);
    for (size_t s = 0; s < nS; ++s) {
      const int i = _tp.spokes[s][0], j = _tp.spokes[s][1];
      _spoke_edge[s] = edge_of[{std::min(i, j), std::max(i, j)}];
      const vec3 d0 = x0[static_cast<size_t>(j)] - x0[static_cast<size_t>(i)];
      _spoke_len[s] = d0.norm();
      _c_off[s] = _T0.c[s] - real(0.5) * _T0.k[static_cast<size_t>(_spoke_edge[s])] * _spoke_len[s];
      _kappa0[s] = _spoke_len[s] > real(0) ? real(2) * _T0.c[s] / _spoke_len[s] : real(0);
      const vec3 Ni = n0[static_cast<size_t>(i)].normalized();
      vec3 dh = d0 - Ni * d0.dot(Ni);
      if (!(dh.norm() > real(0)))
        continue;
      dh.normalize();
      real W = 0;
      for (const auto &sp : stars.at[static_cast<size_t>(i)]) {
        vec3 te = x0[static_cast<size_t>(sp.hi)] - x0[static_cast<size_t>(sp.lo)];
        te -= Ni * te.dot(Ni);
        const real l = te.norm();
        if (!(l > real(1e-20)) || !(sp.w > real(0)))
          continue;
        te /= l;
        const vec3 be = Ni.cross(te);
        const real a = te.dot(dh), b = be.dot(dh);
        _kappa_from[s].push_back({sp.edge, sp.w * vec3(a * a, real(2) * a * b, b * b)});
        W += sp.w;
      }
      for (auto &[e, coef] : _kappa_from[s])
        coef /= std::max(W, real(1e-300));
    }
    if (c_off && c_off->size() == _c_off.size())
      _c_off = *c_off;
    _rest_initialized = true;
  }

  /// The tangency target: `vertex` (default) the one-sided normal curvature
  /// from the store's vertex tensor (follows every edit, all three store
  /// components); `edge` the edge's `k` plus the frame start's one-sided
  /// offset (holds the old shape's offsets under large edits).
  enum class tangency { vertex, edge };
  tangency tangency_target = tangency::vertex;

  std::string name() override { return "contact_net_gn"; }
  bool gauss_newton() const override { return true; }

  void fill_J(const vecX &q, index_t &id0, std::vector<trip> &triplets) override {
    vecX r;
    linearize(q, r, &triplets, id0);
    id0 += static_cast<index_t>(r.size());
  }

  void add_gradient(const vecX &q, vecX &g) override {
    vecX r;
    std::vector<trip> jt;
    linearize(q, r, &jt, 0);
    for (const trip &t : jt)
      g[t.col()] += t.value() * r[t.row()];
  }

  real residual_sq(const vecX &q) override {
    vecX r;
    linearize(q, r, nullptr, 0);
    return r.squaredNorm();
  }

  index_t rows() const {
    return static_cast<index_t>(3 * _tp.edges.size() + _tp.spokes.size()) + _tp.n;
  }

private:
  /// The weighted residual `r` (rows from 0) and, optionally, its Jacobian
  /// triplets with rows offset by `id0` and the solver's global columns.
  void linearize(const vecX &q, vecX &r, std::vector<trip> *jt, index_t id0) const {
    const int n = _tp.n;
    const size_t nE = _tp.edges.size(), nS = _tp.spokes.size();
    std::vector<vec3> x(static_cast<size_t>(n)), N(static_cast<size_t>(n));
    for (int v = 0; v < n; ++v) {
      x[static_cast<size_t>(v)] = _blocks[0]->get_vec3(v, q);
      N[static_cast<size_t>(v)] = _blocks[1]->get_vec3(v, q);
    }
    std::vector<real> dk(nE, 0), dt(nE, 0);
    for (size_t e = 0; e < nE; ++e) {
      const vec3 h = _blocks[2]->get_vec3(static_cast<index_t>(e), q);
      dk[e] = h[0] - _h0[e][0];
      dt[e] = h[1] - _h0[e][1];
    }
    r.resize(rows());
    auto xc = [&](int v) { return _blocks[0]->get_offset_idx(v); };
    auto nc = [&](int v) { return _blocks[1]->get_offset_idx(v); };
    auto sc = [&](size_t e) { return _blocks[2]->get_offset_idx(static_cast<index_t>(e)); };
    Eigen::Index row = 0;
    for (size_t e = 0; e < nE; ++e, row += 3) {
      const int i = _tp.edges[e][0], j = _tp.edges[e][1];
      if (i < 0 || !(_s[e] > real(0))) {
        r.segment<3>(row).setZero();
        continue;
      }
      const real s = _w * _s[e];
      Eigen::Matrix<real, 14, 1> in;
      in << x[static_cast<size_t>(i)], x[static_cast<size_t>(j)], N[static_cast<size_t>(i)], N[static_cast<size_t>(j)],
          _T0.k[e] + dk[e], _T0.tau[e] + dt[e];
      Eigen::Matrix<real, 1, 1> unused(0);
      Eigen::Matrix<real, 3, 1> y;
      Eigen::Matrix<real, 3, 14> J;
      kusama::contact_generated::contact_rodrigues_kt(unused, in, y, J);
      r.segment<3>(row) = s * y;
      if (jt) {
        const index_t cols[4] = {xc(i), xc(j), nc(i), nc(j)};
        for (int a = 0; a < 3; ++a) {
          for (int b = 0; b < 4; ++b)
            for (int c = 0; c < 3; ++c)
              jt->emplace_back(id0 + row + a, cols[b] + c, s * J(a, 3 * b + c));
          jt->emplace_back(id0 + row + a, sc(e) + 0, s * J(a, 12));
          jt->emplace_back(id0 + row + a, sc(e) + 1, s * J(a, 13));
        }
      }
    }
    // Tangency as the Lie "point on sphere" row: `n_i . d / |d| = k |d| / 2`
    // for the sphere of curvature `k` through `x_i` toward `x_j` (exact on a
    // circle), at the live `|d|`, with the frame-start mismatch of the
    // one-sided offset against the edge's `k` held as a constant
    // (`_c_off = c0 - k0 |d0| / 2`), so the target follows edge stretch.
    const bool vertex = tangency_target == tangency::vertex;
    for (size_t sp = 0; sp < nS; ++sp, ++row) {
      const int i = _tp.spokes[sp][0], j = _tp.spokes[sp][1];
      const size_t e = static_cast<size_t>(_spoke_edge[sp]);
      const real s = _w * _s[e];
      const vec3 d = x[static_cast<size_t>(j)] - x[static_cast<size_t>(i)];
      const real len = d.norm();
      real k = _T0.k[e] + dk[e], off = _c_off[sp];
      if (vertex) {
        k = _kappa0[sp];
        for (const auto &[ee, coef] : _kappa_from[sp])
          k += coef.dot(_blocks[2]->get_vec3(ee, q) - _h0[static_cast<size_t>(ee)]);
        off = 0;
      }
      Eigen::Matrix<real, 9, 1> in;
      in << x[static_cast<size_t>(i)], x[static_cast<size_t>(j)], N[static_cast<size_t>(i)];
      Eigen::Matrix<real, 1, 1> c(off + real(0.5) * k * len);
      Eigen::Matrix<real, 1, 1> y;
      Eigen::Matrix<real, 1, 9> J;
      kusama::contact_generated::contact_tangency(c, in, y, J);
      r[row] = s * y[0];
      if (jt && s != real(0)) {
        const vec3 dh = len > real(0) ? vec3(d / len) : vec3::Zero();
        const index_t cols[3] = {xc(i), xc(j), nc(i)};
        for (int b = 0; b < 3; ++b)
          for (int cc = 0; cc < 3; ++cc) {
            real v = J(0, 3 * b + cc);
            if (b == 0)
              v += real(0.5) * k * dh[cc]; // -d(k |d| / 2)/dx_i
            else if (b == 1)
              v -= real(0.5) * k * dh[cc];
            jt->emplace_back(id0 + row, cols[b] + cc, s * v);
          }
        if (vertex) {
          for (const auto &[ee, coef] : _kappa_from[sp])
            for (int cc = 0; cc < 3; ++cc)
              if (coef[cc] != real(0))
                jt->emplace_back(id0 + row, sc(static_cast<size_t>(ee)) + cc, -s * real(0.5) * len * coef[cc]);
        } else {
          jt->emplace_back(id0 + row, sc(e) + 0, -s * real(0.5) * len);
        }
      }
    }
    for (int v = 0; v < n; ++v, ++row) {
      const vec3 &nv = N[static_cast<size_t>(v)];
      const real len = nv.norm();
      const real s = _w * _s_unit;
      r[row] = s * (len - real(1));
      if (jt && len > real(0))
        for (int c = 0; c < 3; ++c)
          jt->emplace_back(id0 + row, nc(v) + c, s * nv[c] / len);
    }
  }

  const kusama::stencil_table *_t = nullptr;
  kusama::contact_net::topology _tp;
  kusama::contact_net::targets _T0;
  std::vector<vec3> _h0;
  std::vector<real> _s, _spoke_len, _c_off, _kappa0;
  std::vector<std::vector<std::pair<int, vec3>>> _kappa_from;
  std::vector<int> _spoke_edge;
  real _s_unit = 1;
};

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif
