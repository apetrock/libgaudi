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
/// Frames mode (MATH.md 19.21, E3): `blocks[1]` a `vec4_block` of per-vertex
/// quaternions `Y` (`shell_cr_block::frame_dof`), the normal
/// `n = Y e3 conj(Y) / |Y|^2` with analytic `dn/dY`; the unit rows hold
/// `|Y| = 1` and a gauge row per vertex (weight `w_frame_gauge`) holds the
/// twist of `conj(Y0) Y` about the body normal (its `z` component, linear in
/// `Y`) at zero. A gauge carried by the shortest arc from the frame-start
/// normal stalled cube -> sphere: under large turns of the normal the arc's
/// holonomy fights the frames.
/// Store-free (frames mode with no third block): the curvature is the
/// frames' own, so there are no Rodrigues rows; each directed spoke's
/// tangency reads the frames' normal curvature,
/// `n_i . d / |d| - k_ij |d| / 2 - c0`, `k_ij = -(n_j - n_i) . t / |d_t|`
/// (`t` the spoke tangent against the mean normal), `c0` the frame start's
/// offset; the native rows (`frame_rows`) act on the frames.
/// Rows scaled `w k_contact_gn A_e / |e|`: the Rodrigues and tangency
/// residuals are curvature times length, so the scaled rows are integrated
/// curvatures, the unit of `weingarten_gn` and `bend`; `k_contact_gn` is
/// calibrated against `bend` (`contact_net_gn_matches_bend_at_equal_weight`).
/// (Lie sphere rows `<s_ij, P_j>`, `<s_ij, T_j>` were tried in place of
/// Rodrigues / tangency and fail off curvature-line edges: the twist leaves
/// the neighbour free to slide on the sphere; MATH.md 19.13.)

#include "gaudi/hepworth/block/block_constraint.hpp"
#include "gaudi/hepworth/block/quat_util.hpp"
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

/// `n = M(Y) / |Y|^2` with `M` the quadratic `Y e3 conj(Y)`, and `dn/dY`.
inline vec3 frame_normal(const vec4 &Y, Eigen::Matrix<real, 3, 4> *dn = nullptr) {
  const real w = Y[0], x = Y[1], y = Y[2], z = Y[3];
  const real l2 = Y.squaredNorm();
  const vec3 M(2 * (x * z + w * y), 2 * (y * z - w * x), w * w - x * x - y * y + z * z);
  if (dn && l2 > real(0)) {
    Eigen::Matrix<real, 3, 4> dM;
    dM << 2 * y, 2 * z, 2 * w, 2 * x, -2 * x, -2 * w, 2 * z, 2 * y, 2 * w, -2 * x, -2 * y, 2 * z;
    *dn = dM / l2 - (real(2) / (l2 * l2)) * M * Y.transpose();
  }
  return l2 > real(0) ? vec3(M / l2) : vec3::UnitZ();
}

class contact_net_gn : public block_constraint {
public:
  DEFINE_CREATE_FUNC(contact_net_gn)

  /// `c_off`: optional per-spoke tangency offsets `c - k |d| / 2` to hold
  /// (default: the frame start's).
  contact_net_gn(const kusama::stencil_table *table, const std::vector<vec3> &x0, const std::vector<vec3> &n0,
                 const std::vector<vec3> &h0, real w, std::vector<sim_block::ptr> blocks, real k = k_contact_gn,
                 const std::vector<real> *c_off = nullptr)
      : block_constraint({}, w, std::move(blocks)), _t(table), _h0(h0) {
    if (auto fb = std::dynamic_pointer_cast<vec4_block>(_blocks[1])) {
      _frames = true;
      _Y0 = fb->_x;
      _store_free = _blocks.size() < 3 || !_blocks[2];
    }
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
    if (_store_free) {
      std::vector<vec3> N0(_Y0.size());
      for (size_t v = 0; v < _Y0.size(); ++v)
        N0[v] = frame_normal(_Y0[v]);
      _c_free.resize(nS);
      for (size_t s = 0; s < nS; ++s) {
        const int i = _tp.spokes[s][0], j = _tp.spokes[s][1];
        _c_free[s] = free_tangency(x0[static_cast<size_t>(i)], x0[static_cast<size_t>(j)], N0[static_cast<size_t>(i)],
                                   N0[static_cast<size_t>(j)]);
      }
    }
    _rest_initialized = true;
  }

  /// The tangency target: `vertex` (default) the one-sided normal curvature
  /// from the store's vertex tensor (follows every edit, all three store
  /// components); `edge` the edge's `k` plus the frame start's one-sided
  /// offset (holds the old shape's offsets under large edits).
  enum class tangency { vertex, edge };
  tangency tangency_target = tangency::vertex;
  /// Frames mode: weight of the gauge rows (relative to the unit rows).
  real w_frame_gauge = real(1e-3);

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
    return static_cast<index_t>((_store_free ? 0 : 3) * _tp.edges.size() + _tp.spokes.size()) +
           (_frames ? 2 : 1) * _tp.n;
  }

private:
  /// Store-free tangency, `n_i . d / |d| - k_ij |d| / 2`, before the offset.
  static real free_tangency(const vec3 &xi, const vec3 &xj, const vec3 &ni, const vec3 &nj) {
    const vec3 d = xj - xi;
    const real len = d.norm();
    const vec3 nb = (ni + nj).normalized();
    const vec3 dt = d - nb * d.dot(nb);
    const real L = dt.norm();
    if (!(len > real(0)) || !(L > real(0)))
      return real(0);
    const real k = -(nj - ni).dot(dt / L) / L;
    return ni.dot(d) / len - real(0.5) * k * len;
  }

  /// The weighted residual `r` (rows from 0) and, optionally, its Jacobian
  /// triplets with rows offset by `id0` and the solver's global columns.
  void linearize(const vecX &q, vecX &r, std::vector<trip> *jt, index_t id0) const {
    const int n = _tp.n;
    const size_t nE = _tp.edges.size(), nS = _tp.spokes.size();
    std::vector<vec3> x(static_cast<size_t>(n)), N(static_cast<size_t>(n));
    std::vector<Eigen::Matrix<real, 3, 4>> dN(_frames ? static_cast<size_t>(n) : 0);
    std::vector<vec4> Y(_frames ? static_cast<size_t>(n) : 0);
    for (int v = 0; v < n; ++v) {
      x[static_cast<size_t>(v)] = _blocks[0]->get_vec3(v, q);
      if (_frames) {
        Y[static_cast<size_t>(v)] = get_elem<vec4>(v, q, _blocks[1]->_offset);
        N[static_cast<size_t>(v)] = frame_normal(Y[static_cast<size_t>(v)], &dN[static_cast<size_t>(v)]);
      } else {
        N[static_cast<size_t>(v)] = _blocks[1]->get_vec3(v, q);
      }
    }
    std::vector<real> dk(nE, 0), dt(nE, 0);
    for (size_t e = 0; e < (_store_free ? 0 : nE); ++e) {
      const vec3 h = _blocks[2]->get_vec3(static_cast<index_t>(e), q);
      dk[e] = h[0] - _h0[e][0];
      dt[e] = h[1] - _h0[e][1];
    }
    r.resize(rows());
    auto xc = [&](int v) { return _blocks[0]->get_offset_idx(v); };
    auto nc = [&](int v) { return _blocks[1]->get_offset_idx(v); };
    // A row's derivative `g` in vertex `v`'s normal, onto the normal's or
    // the frame's columns.
    auto emit_n = [&](index_t r, int v, const vec3 &g) {
      if (_frames) {
        const Eigen::Matrix<real, 1, 4> gy = g.transpose() * dN[static_cast<size_t>(v)];
        for (int c = 0; c < 4; ++c)
          if (gy[c] != real(0))
            jt->emplace_back(r, nc(v) + c, gy[c]);
      } else {
        for (int c = 0; c < 3; ++c)
          jt->emplace_back(r, nc(v) + c, g[c]);
      }
    };
    auto sc = [&](size_t e) { return _blocks[2]->get_offset_idx(static_cast<index_t>(e)); };
    Eigen::Index row = 0;
    for (size_t e = 0; e < (_store_free ? 0 : nE); ++e, row += 3) {
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
        const index_t cols[2] = {xc(i), xc(j)};
        for (int a = 0; a < 3; ++a) {
          for (int b = 0; b < 2; ++b)
            for (int c = 0; c < 3; ++c)
              jt->emplace_back(id0 + row + a, cols[b] + c, s * J(a, 3 * b + c));
          emit_n(id0 + row + a, i, s * vec3(J(a, 6), J(a, 7), J(a, 8)));
          emit_n(id0 + row + a, j, s * vec3(J(a, 9), J(a, 10), J(a, 11)));
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
    for (size_t sp = 0; sp < (_store_free ? nS : 0); ++sp, ++row) {
      const int i = _tp.spokes[sp][0], j = _tp.spokes[sp][1];
      const size_t ii = static_cast<size_t>(i), jj = static_cast<size_t>(j);
      const real s = _w * _s[static_cast<size_t>(_spoke_edge[sp])];
      r[row] = s * (free_tangency(x[ii], x[jj], N[ii], N[jj]) - _c_free[sp]);
      if (!jt || s == real(0))
        continue;
      // Positions: central differences; frames: through `dn/dY`.
      for (int side = 0; side < 2; ++side) {
        const size_t v = side == 0 ? ii : jj;
        vec3 gn;
        for (int c = 0; c < 3; ++c) {
          const real h = real(1e-7) * std::max(real(1), std::abs(x[v][c]));
          std::vector<vec3> xp{x[ii], x[jj]}, xm{x[ii], x[jj]};
          xp[side][c] += h;
          xm[side][c] -= h;
          jt->emplace_back(id0 + row, xc(static_cast<int>(v)) + c,
                           s * (free_tangency(xp[0], xp[1], N[ii], N[jj]) - free_tangency(xm[0], xm[1], N[ii], N[jj])) /
                               (real(2) * h));
          const real hn = real(1e-7);
          vec3 np = N[v], nm = N[v];
          np[c] += hn;
          nm[c] -= hn;
          gn[c] = side == 0 ? (free_tangency(x[ii], x[jj], np, N[jj]) - free_tangency(x[ii], x[jj], nm, N[jj]))
                            : (free_tangency(x[ii], x[jj], N[ii], np) - free_tangency(x[ii], x[jj], N[ii], nm));
          gn[c] /= real(2) * hn;
        }
        emit_n(id0 + row, static_cast<int>(v), s * gn);
      }
    }
    for (size_t sp = 0; sp < (_store_free ? 0 : nS); ++sp, ++row) {
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
        const index_t cols[2] = {xc(i), xc(j)};
        for (int b = 0; b < 2; ++b)
          for (int cc = 0; cc < 3; ++cc) {
            real v = J(0, 3 * b + cc);
            if (b == 0)
              v += real(0.5) * k * dh[cc]; // -d(k |d| / 2)/dx_i
            else
              v -= real(0.5) * k * dh[cc];
            jt->emplace_back(id0 + row, cols[b] + cc, s * v);
          }
        emit_n(id0 + row, i, s * vec3(J(0, 6), J(0, 7), J(0, 8)));
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
    const real su = _w * _s_unit;
    for (int v = 0; v < n; ++v, ++row) {
      if (_frames) {
        const vec4 &yv = Y[static_cast<size_t>(v)];
        const real len = yv.norm();
        r[row] = su * (len - real(1));
        if (jt && len > real(0))
          for (int c = 0; c < 4; ++c)
            jt->emplace_back(id0 + row, nc(v) + c, su * yv[c] / len);
        continue;
      }
      const vec3 &nv = N[static_cast<size_t>(v)];
      const real len = nv.norm();
      r[row] = su * (len - real(1));
      if (jt && len > real(0))
        for (int c = 0; c < 3; ++c)
          jt->emplace_back(id0 + row, nc(v) + c, su * nv[c] / len);
    }
    if (!_frames)
      return;
    const real sg = su * std::sqrt(std::max(w_frame_gauge, real(0)));
    for (int v = 0; v < n; ++v, ++row) {
      // `(conj(Y0) Y)_z = w0 z - z0 w - x0 y + y0 x`.
      const vec4 &a0 = _Y0[static_cast<size_t>(v)];
      const vec4 g(-a0[3], a0[2], -a0[1], a0[0]);
      r[row] = sg * g.dot(Y[static_cast<size_t>(v)]);
      if (jt && sg != real(0))
        for (int c = 0; c < 4; ++c)
          jt->emplace_back(id0 + row, nc(v) + c, sg * g[c]);
    }
  }

  const kusama::stencil_table *_t = nullptr;
  kusama::contact_net::topology _tp;
  kusama::contact_net::targets _T0;
  std::vector<vec3> _h0;
  bool _frames = false;
  std::vector<vec4> _Y0;
  bool _store_free = false;
  std::vector<real> _c_free;
  std::vector<real> _s, _spoke_len, _c_off, _kappa0;
  std::vector<std::vector<std::pair<int, vec3>>> _kappa_from;
  std::vector<int> _spoke_edge;
  real _s_unit = 1;
};

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif
