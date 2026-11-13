#ifndef GAUDI_ENGMAN_CONSTRAINT_HPP
#define GAUDI_ENGMAN_CONSTRAINT_HPP

/// Hard constraints of the Engman solver (`scripts/sympy/cr_poisson/MATH.md`
/// section 17): residuals `C(q) = 0` met by multipliers, not weighed against
/// the energies. Columns are the solver's `q` (the hepworth blocks' layout).

#include "gaudi/hepworth/block/sim_block.hpp"
#include "gaudi/kusama/stencil_table.hpp"
#include "gaudi/kusama/weingarten.hpp"

#include <Eigen/Sparse>

#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace gaudi {
namespace engman {

using trip = Eigen::Triplet<real>;

/// A block of constraint rows. `lambda` is its multiplier, kept here so a
/// caller can carry it between steps (sized by the solver when it does not
/// match `rows()`).
class constraint {
public:
  using ptr = std::shared_ptr<constraint>;
  virtual ~constraint() = default;
  virtual std::string name() const = 0;
  virtual index_t rows() const = 0;
  /// `C(q)` in `C` (`rows()`), and with `J` its Jacobian, rows from 0,
  /// columns in `q`.
  virtual void evaluate(const vecX &q, vecX &C, std::vector<trip> *J) const = 0;
  /// Rows (local, the first of a 3-row block) each paired with the first
  /// `q` column of a 3-column block only that row block reaches, for the
  /// Schur elimination (`kkt_system`). None by default.
  virtual void pairs(std::vector<std::pair<index_t, index_t>> &out) const { (void)out; }

  vecX lambda;
};

/// The Weingarten coupling as a constraint: per edge
/// `C_e = c_e (mean_e R(x) - S_e)`, `c_e = kusama::weingarten::coupling_factor`
/// at `x0`, the `weingarten_gn` row without its weight. `blocks[0]` positions,
/// `blocks[1]` the edge store. Edges touching a degenerate face at `x0` keep
/// zero rows.
class coupling : public constraint {
public:
  coupling(const kusama::stencil_table *table, const std::vector<vec3> &x0,
           std::vector<hepworth::sim_block::ptr> blocks)
      : _t(table), _blocks(std::move(blocks)) {
    _c2 = kusama::weingarten::coupling_factor(*_t, x0);
    for (real &c : _c2)
      c *= c;
    const kusama::weingarten::degenerate_set bad = kusama::weingarten::find_degenerate(*_t, x0);
    _drop.assign(static_cast<size_t>(_t->n_edges), 0);
    for (int e = 0; e < _t->n_edges; ++e)
      if (bad.edge[static_cast<size_t>(e)]) {
        _drop[static_cast<size_t>(e)] = 1;
        _c2[static_cast<size_t>(e)] = 0;
      }
  }
  static ptr create(const kusama::stencil_table *table, const std::vector<vec3> &x0,
                    std::vector<hepworth::sim_block::ptr> blocks) {
    return std::make_shared<coupling>(table, x0, std::move(blocks));
  }

  std::string name() const override { return "engman_coupling"; }
  /// Each edge's rows with its store block (`J_y = -c_e I`).
  void pairs(std::vector<std::pair<index_t, index_t>> &out) const override {
    for (int e = 0; e < _t->n_edges; ++e)
      if (!_drop[static_cast<size_t>(e)] && _c2[static_cast<size_t>(e)] > real(0))
        out.emplace_back(3 * static_cast<index_t>(e), _blocks[1]->get_offset_idx(e));
  }
  index_t rows() const override { return 3 * static_cast<index_t>(_t->n_edges); }

  void evaluate(const vecX &q, vecX &C, std::vector<trip> *J) const override {
    const int n = _t->n_verts, nE = _t->n_edges;
    std::vector<vec3> x(static_cast<size_t>(n)), hs(static_cast<size_t>(nE));
    for (int i = 0; i < n; ++i)
      x[static_cast<size_t>(i)] = _blocks[0]->get_vec3(i, q);
    for (int e = 0; e < nE; ++e)
      hs[static_cast<size_t>(e)] = _blocks[1]->get_vec3(e, q);
    Eigen::VectorXd r;
    std::vector<trip> jt;
    kusama::weingarten::edge_coupling(*_t, x, _c2, &hs, r, J ? &jt : nullptr, 0);
    C = r;
    for (int e = 0; e < nE; ++e)
      if (_drop[static_cast<size_t>(e)] || !C.segment<3>(3 * e).allFinite())
        C.segment<3>(3 * e).setZero();
    if (!J)
      return;
    for (const trip &t : jt) {
      const int e = static_cast<int>(t.row() / 3);
      if (_drop[static_cast<size_t>(e)] || !std::isfinite(t.value()))
        continue;
      J->emplace_back(t.row(), _blocks[0]->get_offset_idx(t.col() / 3) + t.col() % 3, t.value());
    }
    for (int e = 0; e < nE; ++e) {
      const real s = -std::sqrt(std::max(_c2[static_cast<size_t>(e)], real(0)));
      if (s != real(0))
        for (int c = 0; c < 3; ++c)
          J->emplace_back(3 * e + c, _blocks[1]->get_offset_idx(e) + c, s);
    }
  }

private:
  const kusama::stencil_table *_t = nullptr;
  std::vector<hepworth::sim_block::ptr> _blocks;
  std::vector<real> _c2;
  std::vector<char> _drop;
};

/// Positions held exactly: `C = x_i - x_i^goal` for the listed vertices of
/// `block`.
class position_pin : public constraint {
public:
  position_pin(std::vector<index_t> ids, std::vector<vec3> goals, hepworth::sim_block::ptr block)
      : _ids(std::move(ids)), _goals(std::move(goals)), _block(std::move(block)) {}
  static ptr create(std::vector<index_t> ids, std::vector<vec3> goals, hepworth::sim_block::ptr block) {
    return std::make_shared<position_pin>(std::move(ids), std::move(goals), std::move(block));
  }

  std::string name() const override { return "engman_position_pin"; }
  index_t rows() const override { return 3 * static_cast<index_t>(_ids.size()); }

  void evaluate(const vecX &q, vecX &C, std::vector<trip> *J) const override {
    C.resize(rows());
    for (size_t k = 0; k < _ids.size(); ++k) {
      C.segment<3>(3 * static_cast<index_t>(k)) = vec3(_block->get_vec3(_ids[k], q)) - _goals[k];
      if (J)
        for (int c = 0; c < 3; ++c)
          J->emplace_back(3 * static_cast<index_t>(k) + c, _block->get_offset_idx(_ids[k]) + c, real(1));
    }
  }

private:
  std::vector<index_t> _ids;
  std::vector<vec3> _goals;
  hepworth::sim_block::ptr _block;
};

/// Per-face area held exactly: `C_f = (A_f(x) - a*_f) / a_ref`, `a_ref` the
/// mean target area (so `rho_rel` means the same at any resolution). About
/// `2V` rows on `3V` positions, so it does not lock. Its multiplier is the
/// local surface tension.
class face_area : public constraint {
public:
  face_area(std::vector<std::array<int, 3>> tris, std::vector<real> targets, hepworth::sim_block::ptr block)
      : _tris(std::move(tris)), _a(std::move(targets)), _block(std::move(block)) {
    real s = 0;
    for (real a : _a)
      s += a;
    _ref = _a.empty() ? real(1) : std::max(s / real(_a.size()), real(1e-300));
  }
  static ptr create(std::vector<std::array<int, 3>> tris, std::vector<real> targets,
                    hepworth::sim_block::ptr block) {
    return std::make_shared<face_area>(std::move(tris), std::move(targets), std::move(block));
  }

  std::string name() const override { return "engman_face_area"; }
  index_t rows() const override { return static_cast<index_t>(_tris.size()); }

  void evaluate(const vecX &q, vecX &C, std::vector<trip> *J) const override {
    C.resize(rows());
    for (size_t f = 0; f < _tris.size(); ++f) {
      const auto &t = _tris[f];
      const vec3 X[3] = {_block->get_vec3(t[0], q), _block->get_vec3(t[1], q), _block->get_vec3(t[2], q)};
      const vec3 n2 = (X[1] - X[0]).cross(X[2] - X[0]);
      const real A2 = n2.norm();
      C[static_cast<index_t>(f)] = (real(0.5) * A2 - _a[f]) / _ref;
      if (!J || !(A2 > real(0)))
        continue;
      const vec3 n = n2 / A2;
      for (int i = 0; i < 3; ++i) {
        // dA / dX_i = 1/2 n x (X_{i+2} - X_{i+1})
        const vec3 g = real(0.5) * n.cross(X[(i + 2) % 3] - X[(i + 1) % 3]) / _ref;
        for (int c = 0; c < 3; ++c)
          J->emplace_back(static_cast<index_t>(f), _block->get_offset_idx(t[static_cast<size_t>(i)]) + c, g[c]);
      }
    }
  }

  /// The face areas of `x` over `tris`.
  static std::vector<real> areas(const std::vector<std::array<int, 3>> &tris, const std::vector<vec3> &x) {
    std::vector<real> a(tris.size());
    for (size_t f = 0; f < tris.size(); ++f) {
      const auto &t = tris[f];
      a[f] = real(0.5) * (x[static_cast<size_t>(t[1])] - x[static_cast<size_t>(t[0])])
                             .cross(x[static_cast<size_t>(t[2])] - x[static_cast<size_t>(t[0])])
                             .norm();
    }
    return a;
  }

private:
  std::vector<std::array<int, 3>> _tris;
  std::vector<real> _a;
  hepworth::sim_block::ptr _block;
  real _ref = 1;
};

} // namespace engman
} // namespace gaudi

#endif
