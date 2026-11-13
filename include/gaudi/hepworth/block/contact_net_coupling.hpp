#ifndef GAUDI_HEPWORTH_BLOCK_CONTACT_NET_COUPLING_HPP
#define GAUDI_HEPWORTH_BLOCK_CONTACT_NET_COUPLING_HPP

/// The contact element net (`kusama/contact_net.hpp`; scripts/sympy/
/// cr_poisson/MATH.md 19.10 to 19.12) as a projection row on positions
/// (`blocks[0]`), its targets read from the edge store (`blocks[1]`, read
/// only), in the pattern of `exponent_coupling`.
///
/// Rows: the CR derivative `w sqrt(k_gather M_e) D` at the frame start, the
/// same rows and calibration against `bend` as the `exponent` gather, so one
/// `w_coupling` means the same under every curvature solver. `project` runs
/// the net's inner solve when the store's goal moves past `goal_tol` (the
/// store settles over the solver's first iterates), at most
/// `solves_per_frame` times, each warm started from the last net, and
/// writes `D X*`, `X*` the net's positions; otherwise it holds `X*`.
///
/// Targets in increment form, so the row is idle at the frame start: the
/// net's own read-out of the frame-start contact elements plus the store's
/// change since then, in the store's edge frame (`t` lower to higher vertex,
/// `b = n x t`, the net's convention): `k += dh_0`, `tau += dh_1`, and each
/// spoke's normal offset `c += dk |d| / 2` (second order). With no edit the
/// frame is an exact fixed point.

#include "gaudi/hepworth/block/block_constraint.hpp"
#include "gaudi/hepworth/block/exponent_coupling.hpp"
#include "gaudi/kusama/contact_net.hpp"
#include "gaudi/kusama/stencil_table.hpp"
#include "gaudi/kusama/weingarten.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <vector>

namespace gaudi {
namespace hepworth {
namespace block {

class contact_net_coupling : public block_constraint {
public:
  DEFINE_CREATE_FUNC(contact_net_coupling)

  contact_net_coupling(const kusama::stencil_table *table, const std::vector<vec3> &x0,
                       const std::vector<vec3> &h0, real w, std::vector<sim_block::ptr> blocks,
                       kusama::contact_net::config cfg)
      : block_constraint({}, w, std::move(blocks)), _t(table), _cfg(cfg), _h0(h0), _x0(x0) {
    _tp = kusama::contact_net::build_topology(*_t);
    _T0 = kusama::contact_net::read(_tp, x0, kusama::contact_net::vertex_normals(_tp, x0));
    _D = kusama::weingarten::build_cr_derivative(*_t, x0);
    std::map<std::pair<int, int>, int> edge_of;
    for (size_t e = 0; e < _tp.edges.size(); ++e)
      edge_of[{_tp.edges[e][0], _tp.edges[e][1]}] = static_cast<int>(e);
    _spoke_edge.resize(_tp.spokes.size());
    _spoke_len.resize(_tp.spokes.size());
    for (size_t s = 0; s < _tp.spokes.size(); ++s) {
      const int i = _tp.spokes[s][0], j = _tp.spokes[s][1];
      _spoke_edge[s] = edge_of[{std::min(i, j), std::max(i, j)}];
      _spoke_len[s] = (x0[static_cast<size_t>(j)] - x0[static_cast<size_t>(i)]).norm();
    }
    _held = derivative_of(x0);
    _rest_initialized = true;
  }

  std::string name() override { return "contact_net_coupling"; }

  void fill_A(index_t &id0, std::vector<trip> &triplets) override {
    _id0 = id0;
    for (size_t e = 0; e < _D.size(); ++e) {
      const real s = scale(e);
      for (int j = 0; j < 2; ++j)
        for (int c = 0; c < 3; ++c) {
          for (int k = 0; k < 4; ++k) {
            const int v = _D[e].v[static_cast<size_t>(k)];
            if (v < 0)
              continue;
            const real coef = j == 0 ? _D[e].par[static_cast<size_t>(k)] : _D[e].perp[static_cast<size_t>(k)];
            if (coef != real(0))
              triplets.push_back(trip(id0, _blocks[0]->get_offset_idx(v) + c, s * coef));
          }
          ++id0;
        }
    }
  }

  void project(const vecX &q, vecX &p) override {
    const int n = _t->n_verts;
    std::vector<vec3> x(static_cast<size_t>(n));
    for (int v = 0; v < n; ++v)
      x[static_cast<size_t>(v)] = _blocks[0]->get_vec3(v, q);
    kusama::contact_net::targets T = _T0;
    if (_blocks.size() > 1 && _blocks[1]) {
      std::vector<real> dk(T.k.size(), 0);
      for (size_t e = 0; e < T.k.size(); ++e) {
        const vec3 h = _blocks[1]->get_vec3(static_cast<index_t>(e), q);
        dk[e] = h[0] - _h0[e][0];
        T.k[e] += dk[e];
        T.tau[e] += h[1] - _h0[e][1];
      }
      for (size_t s = 0; s < T.c.size(); ++s)
        T.c[s] += real(0.5) * dk[static_cast<size_t>(_spoke_edge[s])] * _spoke_len[s];
    }
    real moved = 0, sc = 0;
    for (size_t e = 0; e < T.k.size(); ++e) {
      const real lk = _goal_k.empty() ? _T0.k[e] : _goal_k[e];
      const real lt = _goal_tau.empty() ? _T0.tau[e] : _goal_tau[e];
      moved = std::max({moved, std::abs(T.k[e] - lk), std::abs(T.tau[e] - lt)});
      sc = std::max(sc, std::abs(T.k[e]));
    }
    if (moved > _cfg.goal_tol * std::max(sc, real(1e-12)) && _solves < _cfg.solves_per_frame) {
      _goal_k = T.k;
      _goal_tau = T.tau;
      // Warm start from the last net (already near the moved goal); the
      // first solve starts at the live positions.
      const std::vector<vec3> &xs = _xs.empty() ? x : _xs;
      const std::vector<vec3> ns = _ns.empty() ? kusama::contact_net::vertex_normals(_tp, x) : _ns;
      const kusama::contact_net::result r = kusama::contact_net::solve(_tp, xs, ns, T, _cfg);
      _xs = r.x;
      _ns = r.n;
      _held = derivative_of(r.x);
      _last_merit = r.merit.empty() ? real(0) : r.merit.back();
      _solves += 1;
      _iterations += r.iterations;
      _factorizations += r.factorizations;
    }
    write(_held, p);
  }

  real last_merit() const { return _last_merit; }
  int solves() const { return _solves; }
  int iterations() const { return _iterations; }
  int factorizations() const { return _factorizations; }

private:
  real scale(size_t e) const { return _w * std::sqrt(exponent::k_gather * std::max(_D[e].mass, real(0))); }

  /// `D x` with the frame-start derivative (the rows `fill_A` writes).
  std::vector<vec3> derivative_of(const std::vector<vec3> &x) const {
    const size_t nE = _D.size();
    std::vector<vec3> out(2 * nE, vec3::Zero());
    for (size_t e = 0; e < nE; ++e)
      for (int k = 0; k < 4; ++k) {
        const int v = _D[e].v[static_cast<size_t>(k)];
        if (v < 0)
          continue;
        out[e] += _D[e].par[static_cast<size_t>(k)] * x[static_cast<size_t>(v)];
        out[nE + e] += _D[e].perp[static_cast<size_t>(k)] * x[static_cast<size_t>(v)];
      }
    return out;
  }

  void write(const std::vector<vec3> &G, vecX &p) const {
    const size_t nE = _D.size();
    index_t row = _id0;
    for (size_t e = 0; e < nE; ++e) {
      const real s = scale(e);
      for (int j = 0; j < 2; ++j)
        for (int c = 0; c < 3; ++c)
          p[row++] = s * G[static_cast<size_t>(j) * nE + e][c];
    }
  }

  const kusama::stencil_table *_t = nullptr;
  kusama::contact_net::config _cfg;
  kusama::contact_net::topology _tp;
  kusama::contact_net::targets _T0;
  std::vector<kusama::weingarten::cr_edge_rows> _D;
  std::vector<vec3> _h0, _x0, _held, _xs, _ns;
  std::vector<int> _spoke_edge;
  std::vector<real> _spoke_len, _goal_k, _goal_tau;
  real _last_merit = 0;
  int _solves = 0, _iterations = 0, _factorizations = 0;
};

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif
