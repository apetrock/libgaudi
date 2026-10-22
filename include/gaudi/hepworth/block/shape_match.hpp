#ifndef GAUDI_HEPWORTH_BLOCK_SHAPE_MATCH_HPP
#define GAUDI_HEPWORTH_BLOCK_SHAPE_MATCH_HPP

/// Shape matching (Müller, Heidelberger, Teschner, Gross 2005) as one
/// projective-dynamics constraint (Bouaziz et al. 2014). A goal function
/// emits goal positions from the live state; the goal is aligned onto the
/// current vertices (none / rigid / linear) and written as pin rows or
/// edge-difference rows.

#include "gaudi/hepworth/block/block_constraint.hpp"

#include <Eigen/SVD>

#include <array>
#include <functional>
#include <utility>
#include <vector>

namespace gaudi {
namespace hepworth {
namespace block {

/// Emits goal positions `g[k]` for `ids[k]` from the live state `q`. Blocks
/// it reads are captured, the way constraints hold their blocks. Returning
/// false means no pull this iterate.
using shape_goal_fn = std::function<bool(
    const vecX &q, const std::vector<index_t> &ids, std::vector<vec3> &g)>;

/// Classic shape matching: the goal is a stored rest shape.
inline shape_goal_fn rest_shape_goal(std::vector<vec3> rest) {
  return [rest = std::move(rest)](const vecX &, const std::vector<index_t> &ids,
                                  std::vector<vec3> &g) {
    g.resize(ids.size());
    for (size_t k = 0; k < ids.size(); ++k) {
      if (ids[k] < 0 || static_cast<size_t>(ids[k]) >= rest.size())
        return false;
      g[k] = rest[static_cast<size_t>(ids[k])];
    }
    return true;
  };
}

/// Müller §4: `none` keeps the goal as emitted; `rigid` removes its rigid
/// motion relative to the current vertices (polar factor); `linear` blends
/// the volume-normalized best-fit linear map with that rotation by `beta`.
enum class shape_align { none, rigid, linear };

/// `pin`: one row per vertex. `edge`: one row per edge, translation free.
enum class shape_rows { pin, edge };

/// Shape-matching knobs, independent of the goal and the mesh.
struct shape_match_config {
  /// Stiffness relative to inertia: rows are `sqrt(w m_i) / h`, so 1 pulls
  /// as hard as the mass term `M / h^2`, whatever the step.
  real w = 1.0;
  shape_rows rows = shape_rows::pin;
  /// `rigid` is Procrustes (polar factor); `linear` blends in the
  /// best-fit linear map by `beta`.
  shape_align align = shape_align::rigid;
  real beta = 0;
};

class shape_match : public block_constraint {
public:
  DEFINE_CREATE_FUNC(shape_match)

  /// `blocks[0]` is the positions block the rows act on. `h` is the solver
  /// step, which makes `cfg.w` inertia-relative. `edges` are pairs of local
  /// indices into `ids`, used by `shape_rows::edge`.
  shape_match(const std::vector<index_t> &ids, const shape_match_config &cfg,
              real h, std::vector<sim_block::ptr> blocks, shape_goal_fn goal,
              std::vector<std::array<int, 2>> edges = {})
      : block_constraint(ids, cfg.w / (h * h), std::move(blocks)),
        _goal(std::move(goal)), _rows(cfg.rows), _align(cfg.align),
        _edges(std::move(edges)), _beta(cfg.beta) {
    _rest_initialized = true;
  }

  std::string name() override { return "shape_match"; }

  void fill_A(index_t &id0, std::vector<trip> &triplets) override {
    _id0 = id0;
    const size_t n = _ids.size();
    _m.assign(n, real(1));
    if (auto pb = std::dynamic_pointer_cast<vec3_block>(_blocks[0]))
      for (size_t k = 0; k < n; ++k)
        if (_ids[k] >= 0 && static_cast<size_t>(_ids[k]) < pb->_M.size())
          _m[k] = pb->_M[static_cast<size_t>(_ids[k])][0];
    if (_rows == shape_rows::pin) {
      _s.resize(n);
      for (size_t k = 0; k < n; ++k) {
        _s[k] = std::sqrt(_w * _m[k]);
        const index_t c = _blocks[0]->get_offset_idx(_ids[k]);
        for (int ax = 0; ax < 3; ++ax)
          triplets.push_back(trip(id0 + 3 * k + ax, c + ax, _s[k]));
      }
      id0 += 3 * static_cast<index_t>(n);
      return;
    }
    _s.resize(_edges.size());
    for (size_t e = 0; e < _edges.size(); ++e) {
      const int a = _edges[e][0], b = _edges[e][1];
      _s[e] = std::sqrt(_w * real(0.5) * (_m[a] + _m[b]));
      const index_t ca = _blocks[0]->get_offset_idx(_ids[a]);
      const index_t cb = _blocks[0]->get_offset_idx(_ids[b]);
      for (int ax = 0; ax < 3; ++ax) {
        triplets.push_back(trip(id0 + 3 * e + ax, ca + ax, -_s[e]));
        triplets.push_back(trip(id0 + 3 * e + ax, cb + ax, _s[e]));
      }
    }
    id0 += 3 * static_cast<index_t>(_edges.size());
  }

  void project(const vecX &q, vecX &p) override {
    const size_t n = _ids.size();
    std::vector<vec3> x(n);
    for (size_t k = 0; k < n; ++k)
      x[k] = _blocks[0]->get_vec3(_ids[k], q);
    std::vector<vec3> g;
    const bool ok = _goal && _goal(q, _ids, g) && g.size() == n;
    if (ok)
      align(x, g);
    const std::vector<vec3> &t = ok ? g : x;
    if (_rows == shape_rows::pin) {
      for (size_t k = 0; k < n; ++k)
        p.block(_id0 + 3 * k, 0, 3, 1) = _s[k] * t[k];
      return;
    }
    for (size_t e = 0; e < _edges.size(); ++e)
      p.block(_id0 + 3 * e, 0, 3, 1) =
          _s[e] * (t[_edges[e][1]] - t[_edges[e][0]]);
  }

private:
  void align(const std::vector<vec3> &x, std::vector<vec3> &g) const {
    if (_align == shape_align::none)
      return;
    real M = 0;
    vec3 cx = vec3::Zero(), cg = vec3::Zero();
    for (size_t k = 0; k < x.size(); ++k) {
      M += _m[k];
      cx += _m[k] * x[k];
      cg += _m[k] * g[k];
    }
    if (!(M > real(0)))
      return;
    cx /= M;
    cg /= M;
    mat3 Apq = mat3::Zero(), Aqq = mat3::Zero();
    for (size_t k = 0; k < x.size(); ++k) {
      const vec3 pk = x[k] - cx, qk = g[k] - cg;
      Apq += _m[k] * pk * qk.transpose();
      Aqq += _m[k] * qk * qk.transpose();
    }
    Eigen::JacobiSVD<mat3> svd(Apq, Eigen::ComputeFullU | Eigen::ComputeFullV);
    mat3 D = mat3::Identity();
    D(2, 2) = (svd.matrixU() * svd.matrixV().transpose()).determinant() < 0
                  ? real(-1)
                  : real(1);
    const mat3 R = svd.matrixU() * D * svd.matrixV().transpose();
    mat3 T = R;
    if (_align == shape_align::linear && _beta > real(0)) {
      mat3 A = Apq * Aqq.inverse();
      const real det = A.determinant();
      if (det > real(0) && A.allFinite()) {
        A /= std::cbrt(det);
        T = _beta * A + (real(1) - _beta) * R;
      }
    }
    for (size_t k = 0; k < g.size(); ++k)
      g[k] = T * (g[k] - cg) + cx;
  }

  shape_goal_fn _goal;
  shape_rows _rows;
  shape_align _align;
  std::vector<std::array<int, 2>> _edges;
  real _beta;
  std::vector<real> _m;
  std::vector<real> _s;
};

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif
