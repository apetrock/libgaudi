#ifndef GAUDI_HEPWORTH_BLOCK_AMAP_GN_HPP
#define GAUDI_HEPWORTH_BLOCK_AMAP_GN_HPP

/// As-Moebius-as-possible rows (Vaxman, Mueller, Weber 2015, *Conformal Mesh
/// Deformations with Moebius Transformations*, section 6.2, the 3D "edge
/// deviation" form; `reference/vaxman2015-conformal-mobius-deformations.pdf`).
///
/// Per vertex a quaternion `Y_i` (not unit: its modulus is the local scale),
/// per edge `(i, k)` the Gauss-Newton residual
///
///     w_ik - conj(Y_i) q_ik Y_k
///
/// with `w_ik = x_k - x_i` the live edge and `q_ik` the rest edge (frame start;
/// both pure quaternions; 4 rows). Equal `Y_i = Y_k` is a similarity of the
/// edge (rotation and scale); unequal ones are the Moebius (inversion) part,
/// damped by `E_INV = alpha sum |Y_i - Y_k|^2` (4 rows per edge). At the
/// frame start (`Y = 1`, `w = q`) every row is zero. Unlike
/// `triangle_conformal` (each face's deformation gradient projected onto its
/// closest similarity), this is vertex-based, couples neighbours through the
/// shared `Y`, and lets Moebius maps through exactly.
///
/// Rows scaled `w sqrt(A_e) / |q_ik|` (the residual is a length; per unit
/// length it is dimensionless, times the root of the edge's dual area like
/// `triangle_conformal`'s `sqrt(A)`). Blocks: `blocks[0]` positions,
/// `blocks[1]` the quaternions (`vec4_block`, `(w, x, y, z)`).

#include "gaudi/hepworth/block/block_constraint.hpp"
#include "gaudi/kusama/stencil_table.hpp"

#include <array>
#include <cmath>
#include <vector>

namespace gaudi {
namespace hepworth {
namespace block {

namespace amap_detail {

/// Hamilton product, `(w, x, y, z)`.
inline vec4 qmul(const vec4 &a, const vec4 &b) {
  return vec4(a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
              a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
              a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1],
              a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0]);
}
inline vec4 qconj(const vec4 &a) { return vec4(a[0], -a[1], -a[2], -a[3]); }
inline vec4 pure(const vec3 &v) { return vec4(0, v[0], v[1], v[2]); }
inline vec4 unit(int c) {
  vec4 e = vec4::Zero();
  e[c] = 1;
  return e;
}

} // namespace amap_detail

class amap_gn : public block_constraint {
public:
  DEFINE_CREATE_FUNC(amap_gn)

  /// `x0`: the rest (frame-start) positions; `alpha_inv`: the inversion
  /// damping `alpha_INV`, relative to the AMAP rows.
  amap_gn(const kusama::stencil_table *table, const std::vector<vec3> &x0, real w, real alpha_inv,
          std::vector<sim_block::ptr> blocks)
      : block_constraint({}, w, std::move(blocks)), _alpha(alpha_inv) {
    const size_t nE = static_cast<size_t>(table->n_edges);
    _ends.assign(nE, {-1, -1});
    std::vector<real> area(nE, real(0));
    for (int f = 0; f < table->n_faces; ++f) {
      const auto &id = table->ids[static_cast<size_t>(f)];
      const real a = real(0.5) * (x0[static_cast<size_t>(id[1])] - x0[static_cast<size_t>(id[0])])
                                     .cross(x0[static_cast<size_t>(id[2])] - x0[static_cast<size_t>(id[0])])
                                     .norm();
      for (int i = 0; i < 3; ++i) {
        const int e = table->face_edges[static_cast<size_t>(f)][static_cast<size_t>(i)];
        if (e < 0 || e >= static_cast<int>(nE))
          continue;
        _ends[static_cast<size_t>(e)] = {id[static_cast<size_t>(i)], id[static_cast<size_t>((i + 1) % 3)]};
        area[static_cast<size_t>(e)] += a / real(3);
      }
    }
    _q.assign(nE, vec3::Zero());
    _s.assign(nE, real(0));
    for (size_t e = 0; e < nE; ++e) {
      const auto [i, k] = _ends[e];
      if (i < 0)
        continue;
      _q[e] = x0[static_cast<size_t>(k)] - x0[static_cast<size_t>(i)];
      const real l = _q[e].norm();
      _s[e] = l > real(0) ? std::sqrt(area[e]) / l : real(0);
    }
    _rest_initialized = true;
  }

  std::string name() override { return "amap_gn"; }
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

private:
  /// Residual rows (already weighted) and, if `jt`, their Jacobian at global
  /// columns, rows offset by `id0`.
  void linearize(const vecX &q, vecX &r, std::vector<trip> *jt, index_t id0) const {
    using namespace amap_detail;
    const size_t nE = _ends.size();
    const bool inv = _alpha > real(0);
    const size_t per = inv ? 8 : 4;
    r = vecX::Zero(static_cast<index_t>(per * nE));
    const sim_block &P = *_blocks[0];
    const sim_block &Yb = *_blocks[1];
    const real sa = std::sqrt(std::max(_alpha, real(0)));
    for (size_t e = 0; e < nE; ++e) {
      const auto [i, k] = _ends[e];
      if (i < 0 || !(_s[e] > real(0)))
        continue;
      const real s = _w * _s[e];
      const vec4 Yi = get_elem<vec4>(i, q, Yb._offset), Yk = get_elem<vec4>(k, q, Yb._offset);
      const vec4 Q = pure(_q[e]);
      const vec4 W = pure(P.get_vec3(k, q) - P.get_vec3(i, q));
      const vec4 QYk = qmul(Q, Yk), cYiQ = qmul(qconj(Yi), Q);
      const vec4 res = W - qmul(qconj(Yi), QYk);
      const index_t row = static_cast<index_t>(per * e);
      for (int c = 0; c < 4; ++c)
        r[row + c] = s * res[c];
      if (jt) {
        const index_t xi = P.get_offset_idx(i), xk = P.get_offset_idx(k);
        for (int a = 0; a < 3; ++a) {
          jt->emplace_back(id0 + row + 1 + a, xk + a, s);
          jt->emplace_back(id0 + row + 1 + a, xi + a, -s);
        }
        const index_t yi = Yb.get_offset_idx(i), yk = Yb.get_offset_idx(k);
        for (int c = 0; c < 4; ++c) {
          // d(conj(Y_i) Q Y_k) / dY_i,c = conj(e_c) Q Y_k; / dY_k,c = conj(Y_i) Q e_c.
          const vec4 di = qmul(qconj(unit(c)), QYk), dk = qmul(cYiQ, unit(c));
          for (int m = 0; m < 4; ++m) {
            if (di[m] != real(0))
              jt->emplace_back(id0 + row + m, yi + c, -s * di[m]);
            if (dk[m] != real(0))
              jt->emplace_back(id0 + row + m, yk + c, -s * dk[m]);
          }
        }
      }
      if (inv) {
        // `Y` is dimensionless: weight by the root of the dual area alone.
        const real si = sa * _w * _s[e] * _q[e].norm();
        const vec4 d = Yi - Yk;
        for (int c = 0; c < 4; ++c)
          r[row + 4 + c] = si * d[c];
        if (jt) {
          const index_t yi = Yb.get_offset_idx(i), yk = Yb.get_offset_idx(k);
          for (int c = 0; c < 4; ++c) {
            jt->emplace_back(id0 + row + 4 + c, yi + c, si);
            jt->emplace_back(id0 + row + 4 + c, yk + c, -si);
          }
        }
      }
    }
  }

  real _alpha = 0;
  std::vector<std::array<int, 2>> _ends;
  std::vector<vec3> _q;
  std::vector<real> _s;
};

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif
