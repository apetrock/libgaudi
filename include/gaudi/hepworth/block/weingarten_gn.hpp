#ifndef GAUDI_HEPWORTH_BLOCK_WEINGARTEN_GN_HPP
#define GAUDI_HEPWORTH_BLOCK_WEINGARTEN_GN_HPP

/// The Weingarten solve (`kusama/weingarten.hpp`) as Gauss–Newton rows of the
/// projection solve, weighted like every other row: `A = w * factor`, with
/// the factor the row's own geometry.
/// - `weingarten_gn`: the read-out row, the store `S` against the mean of
///   each edge's two faces' live CR read-outs (positions and store columns);
/// - `dirichlet_memory_gn`: the intrinsic vector Dirichlet memory on the
///   positions, rest at construction (the frame start: frames are restarts).

#include "gaudi/hepworth/block/block_constraint.hpp"
#include "gaudi/kusama/stencil_table.hpp"
#include "gaudi/kusama/weingarten.hpp"

#include <cmath>
#include <vector>

namespace gaudi {
namespace hepworth {
namespace block {

namespace detail {

inline void gather_positions(const sim_block &pos, const vecX &q, int n,
                             std::vector<vec3> &x) {
  x.resize(static_cast<size_t>(n));
  for (int i = 0; i < n; ++i)
    x[static_cast<size_t>(i)] = pos.get_vec3(i, q);
}

inline real mean_face_area(const kusama::stencil_table &t, const std::vector<vec3> &x) {
  real s = 0;
  for (int f = 0; f < t.n_faces; ++f) {
    const auto &id = t.ids[static_cast<size_t>(f)];
    s += real(0.5) * (x[static_cast<size_t>(id[1])] - x[static_cast<size_t>(id[0])])
                         .cross(x[static_cast<size_t>(id[2])] - x[static_cast<size_t>(id[0])])
                         .norm();
  }
  return t.n_faces > 0 ? s / t.n_faces : real(0);
}

/// Zeroes the residual rows in `drop` and removes their Jacobian entries,
/// along with any non-finite entry, so one degenerate stencil cannot poison
/// the system.
inline void drop_rows(const std::vector<char> &drop, Eigen::VectorXd &r, std::vector<trip> *jt) {
  for (index_t i = 0; i < r.size(); ++i)
    if ((static_cast<size_t>(i) < drop.size() && drop[static_cast<size_t>(i)]) ||
        !std::isfinite(r[i]))
      r[i] = 0;
  if (!jt)
    return;
  jt->erase(std::remove_if(jt->begin(), jt->end(),
                           [&](const trip &t) {
                             const size_t row = static_cast<size_t>(t.row());
                             return (row < drop.size() && drop[row]) || !std::isfinite(t.value());
                           }),
            jt->end());
}

} // namespace detail

/// The coupling row between position space and shape space, like the rod
/// rows that keep positions and frames in sync. Rows `3 e + c`:
/// `w c_e (mean_e(R(x)) - S_e)`, `c_e = kusama::weingarten::coupling_factor`
/// (`k_coupling` times the edge's dual area) at `x0`.
/// `blocks[0]` positions, `blocks[1]` the edge store
/// (`shell_cr_block::shape_dof`). Shape rows on the store with the same
/// `c_e` (`edge_shape_fair`) combine with it in the least-squares solve:
/// pins give `S = (w^2 R + sum w_k^2 g_k) / (w^2 + sum w_k^2)`.
class weingarten_gn : public block_constraint {
public:
  DEFINE_CREATE_FUNC(weingarten_gn)

  weingarten_gn(const kusama::stencil_table *table, const std::vector<vec3> &x0, real w,
                std::vector<sim_block::ptr> blocks)
      : block_constraint({}, w, std::move(blocks)), _t(table) {
    // `edge_coupling` scales edge rows by `sqrt(mass)`: pass `c_e^2`.
    _c2 = kusama::weingarten::coupling_factor(*_t, x0);
    for (real &c : _c2)
      c *= c;
    // Edges touching a degenerate face at the frame start get no rows.
    const kusama::weingarten::degenerate_set bad = kusama::weingarten::find_degenerate(*_t, x0);
    _drop.assign(3 * static_cast<size_t>(_t->n_edges), 0);
    for (int e = 0; e < _t->n_edges; ++e)
      if (bad.edge[static_cast<size_t>(e)]) {
        _c2[static_cast<size_t>(e)] = 0;
        for (int c = 0; c < 3; ++c)
          _drop[3 * static_cast<size_t>(e) + c] = 1;
      }
    if (bad.n_faces > 0)
      std::cout << "[weingarten_gn] " << bad.n_faces << " degenerate face(s) (worst " << bad.worst_face
                << ", quality " << bad.worst << "): " << bad.n_edges
                << " edge(s) left out of the coupling" << std::endl;
    _rest_initialized = true;
  }

  std::string name() override { return "weingarten_gn"; }
  bool gauss_newton() const override { return true; }

  void fill_J(const vecX &q, index_t &id0, std::vector<trip> &triplets) override {
    vecX r;
    std::vector<trip> jt;
    linearize(q, r, &jt);
    for (const trip &t : jt)
      triplets.emplace_back(id0 + t.row(), column(t.col()), _w * t.value());
    id0 += static_cast<index_t>(r.size());
  }

  void add_gradient(const vecX &q, vecX &g) override {
    vecX r;
    std::vector<trip> jt;
    linearize(q, r, &jt);
    for (const trip &t : jt)
      g[column(t.col())] += _w * _w * t.value() * r[t.row()];
  }

  real residual_sq(const vecX &q) override {
    vecX r;
    linearize(q, r, nullptr);
    return _w * _w * r.squaredNorm();
  }

private:
  /// Solver-local columns: positions `3 v + c`, then the store `3 (n + e) + c`.
  void linearize(const vecX &q, vecX &r, std::vector<trip> *jt) const {
    std::vector<vec3> x, hs(static_cast<size_t>(_t->n_edges));
    detail::gather_positions(*_blocks[0], q, _t->n_verts, x);
    for (int e = 0; e < _t->n_edges; ++e)
      hs[static_cast<size_t>(e)] = _blocks[1]->get_vec3(e, q);
    Eigen::VectorXd rc;
    kusama::weingarten::edge_coupling(*_t, x, _c2, &hs, rc, jt, 0);
    detail::drop_rows(_drop, rc, jt);
    r = rc;
    if (!jt)
      return;
    const index_t hcol = 3 * static_cast<index_t>(_t->n_verts);
    for (int e = 0; e < _t->n_edges; ++e) {
      const real s = -std::sqrt(std::max(_c2[static_cast<size_t>(e)], real(0)));
      if (s != real(0))
        for (int c = 0; c < 3; ++c)
          jt->emplace_back(3 * e + c, hcol + 3 * e + c, s);
    }
  }

  index_t column(index_t col) const {
    const index_t n3 = 3 * static_cast<index_t>(_t->n_verts);
    if (col < n3)
      return _blocks[0]->get_offset_idx(col / 3) + col % 3;
    return _blocks[1]->get_offset_idx((col - n3) / 3) + (col - n3) % 3;
  }

  const kusama::stencil_table *_t = nullptr;
  std::vector<real> _c2;
  std::vector<char> _drop;
};

/// Per-face rows `K_f^{1/2} c_f` (12) plus per-edge mass rows of the
/// intrinsic vector Dirichlet memory, rest at `x0`, times `w sqrt(A)` (`A` the
/// mean face area at `x0`: the dimensionless frame change times a length).
/// `blocks[0]` positions.
class dirichlet_memory_gn : public block_constraint {
public:
  DEFINE_CREATE_FUNC(dirichlet_memory_gn)

  dirichlet_memory_gn(const kusama::stencil_table *table, const std::vector<vec3> &x0,
                      real w, std::vector<sim_block::ptr> blocks, real mu_rel = real(1e-3))
      : block_constraint({}, w, std::move(blocks)), _t(table),
        _rest(kusama::weingarten::build_memory_rest(*table, x0, mu_rel)),
        _l(std::sqrt(detail::mean_face_area(*table, x0))) {
    // Rows: 12 per face, then 3 mass rows per edge (`memory_residual`). An
    // edge's frame uses its diamond and wing vertices (`_rest.D[e].v`): drop
    // the rows of every edge whose stencil reaches a degenerate face, and the
    // rows of every face holding such an edge.
    const kusama::weingarten::degenerate_set bad = kusama::weingarten::find_degenerate(*_t, x0);
    const size_t nF = static_cast<size_t>(_t->n_faces), nE = static_cast<size_t>(_t->n_edges);
    _drop.assign(12 * nF + 3 * nE, 0);
    if (bad.n_faces > 0) {
      std::vector<char> edge_bad(nE, 0);
      for (size_t e = 0; e < nE && e < _rest.D.size(); ++e)
        for (int v : _rest.D[e].v)
          if (v >= 0 && static_cast<size_t>(v) < bad.vert.size() && bad.vert[static_cast<size_t>(v)])
            edge_bad[e] = 1;
      for (size_t e = 0; e < nE; ++e)
        if (edge_bad[e] || bad.edge[e])
          for (int c = 0; c < 3; ++c)
            _drop[12 * nF + 3 * e + c] = 1;
      for (size_t f = 0; f < nF; ++f) {
        bool any = bad.face[f];
        for (int i = 0; i < 3 && !any; ++i) {
          const int e = _t->face_edges[f][static_cast<size_t>(i)];
          any = e >= 0 && static_cast<size_t>(e) < nE && edge_bad[static_cast<size_t>(e)];
        }
        if (any)
          for (int k = 0; k < 12; ++k)
            _drop[12 * f + k] = 1;
      }
    }
    _rest_initialized = true;
  }

  std::string name() override { return "dirichlet_memory_gn"; }
  bool gauss_newton() const override { return true; }

  void fill_J(const vecX &q, index_t &id0, std::vector<trip> &triplets) override {
    vecX r;
    std::vector<trip> jt;
    linearize(q, r, &jt);
    const real s = _w * _l;
    for (const trip &t : jt)
      triplets.emplace_back(id0 + t.row(), column(t.col()), s * t.value());
    id0 += static_cast<index_t>(r.size());
  }

  void add_gradient(const vecX &q, vecX &g) override {
    vecX r;
    std::vector<trip> jt;
    linearize(q, r, &jt);
    const real s2 = _w * _w * _l * _l;
    for (const trip &t : jt)
      g[column(t.col())] += s2 * t.value() * r[t.row()];
  }

  real residual_sq(const vecX &q) override {
    vecX r;
    linearize(q, r, nullptr);
    return _w * _w * _l * _l * r.squaredNorm();
  }

private:
  void linearize(const vecX &q, vecX &r, std::vector<trip> *jt) const {
    std::vector<vec3> x;
    detail::gather_positions(*_blocks[0], q, _t->n_verts, x);
    Eigen::VectorXd rm;
    kusama::weingarten::memory_residual(*_t, _rest, x, real(1), rm, jt, 0);
    detail::drop_rows(_drop, rm, jt);
    r = rm;
  }

  index_t column(index_t col) const {
    return _blocks[0]->get_offset_idx(col / 3) + col % 3;
  }

  const kusama::stencil_table *_t = nullptr;
  kusama::weingarten::memory_rest _rest;
  real _l = 1;
  std::vector<char> _drop;
};

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif
