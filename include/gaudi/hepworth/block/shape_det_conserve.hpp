#ifndef GAUDI_HEPWORTH_BLOCK_SHAPE_DET_CONSERVE_HPP
#define GAUDI_HEPWORTH_BLOCK_SHAPE_DET_CONSERVE_HPP

/// Conservation of the edge store's Gauss-image area over each vertex ring:
/// `sum_e c_e det(h_e) = C0`, `det h = k11 k22 - k12^2` (the Gaussian
/// curvature `K`), `c_e` the edge's dual area (`coupling_factor`), `C0` from
/// the frame-start store. It is the shape-space match of the `area` row: a
/// ring keeps its total `K dA` and redistributes it between its edges and
/// between their mean curvature and anisotropy. `K` is intrinsic, so the
/// held quantity is the one bending without stretching keeps.
///
/// `vertex_gauss_conserve` holds pointwise `K` instead: each vertex's
/// Gaussian curvature `det h_v` (the vertex tensor of its spokes) at a rest
/// `K0` the caller keeps (first frame, relaxed or not). Pointwise `K` scales
/// as `1 / s^2`, so a held rest blocks a uniform shrink, which `K dA` does
/// not.

#include "gaudi/hepworth/block/block_constraint.hpp"
#include "gaudi/kusama/stencil_table.hpp"
#include "gaudi/kusama/weingarten.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace gaudi {
namespace hepworth {
namespace block {

/// The multiplier `l` in `(-1, 1)` with `sH / (1 - l)^2 - sD / (1 + l)^2 = C0`
/// (`sH`, `sD` the weighted sums of `H^2` and `|delta|^2`): the closest
/// tensors with `sum det = C0` are `H / (1 - l)`, `delta / (1 + l)`. The left
/// side increases with `l`, so Newton from 0 safeguarded by the bracket
/// finds it; a root past the bracket (no mean part, or no deviator) stops
/// at its end.
inline real det_multiplier(real sH, real sD, real C0) {
  auto f = [&](real l) {
    return sH / ((real(1) - l) * (real(1) - l)) - sD / ((real(1) + l) * (real(1) + l)) - C0;
  };
  auto df = [&](real l) {
    return real(2) * sH / ((real(1) - l) * (real(1) - l) * (real(1) - l)) +
           real(2) * sD / ((real(1) + l) * (real(1) + l) * (real(1) + l));
  };
  real lo = real(-1) + real(1e-6), hi = real(1) - real(1e-6);
  if (!(f(lo) < real(0)))
    return lo;
  if (!(f(hi) > real(0)))
    return hi;
  const real tol = real(1e-14) * (std::abs(sH) + std::abs(sD) + std::abs(C0) + real(1e-300));
  real l = real(0);
  for (int it = 0; it < 100; ++it) {
    const real fl = f(l);
    if (std::abs(fl) <= tol)
      break;
    (fl < real(0) ? lo : hi) = l;
    const real d = df(l);
    real next = d > real(0) ? l - fl / d : real(0.5) * (lo + hi);
    if (!(next > lo && next < hi))
      next = real(0.5) * (lo + hi);
    if (std::abs(next - l) <= real(1e-16))
      break;
    l = next;
  }
  return l;
}

/// The closest tensor to `h` (tensor metric, `k11^2 + 2 k12^2 + k22^2`) with
/// `det = K0`: the eigenvectors kept and the eigenvalues `H +- |delta|`
/// moved to the nearest point of `l1 l2 = K0`, signs kept (saddles stay on
/// their branch).
inline vec3 closest_det_tensor(const vec3 &h, real K0) {
  const real H = real(0.5) * (h[0] + h[2]);
  const real d = real(0.5) * (h[0] - h[2]);
  const real l = det_multiplier(H * H, d * d + h[1] * h[1], K0);
  const real Hn = H / (real(1) - l), dn = d / (real(1) + l);
  return vec3(Hn + dn, h[1] / (real(1) + l), Hn - dn);
}

/// One ring's conservation row on the edge store (`blocks[0]`, one vec3 per
/// edge). Rows are the tensor metric `k11^2 + 2 k12^2 + k22^2 = 2 H^2 + 2
/// |delta|^2` weighted by `sqrt(c_e)`. The projection is the closest ring
/// (in that metric) with `sum c_e (H^2 - |delta|^2) = C0`: by Lagrange every
/// edge takes `H / (1 - lambda)` and `delta / (1 + lambda)` with one
/// `lambda` in `(-1, 1)`, the root of a function increasing in `lambda`
/// (bisection, then the endpoints are never reached). `C0` is fixed for the
/// step, so the row is variational.
class ring_det_conserve : public block_constraint {
public:
  DEFINE_CREATE_FUNC(ring_det_conserve)

  ring_det_conserve(std::vector<index_t> edges, std::vector<real> c, real C0, real w,
                    std::vector<sim_block::ptr> blocks)
      : block_constraint(std::move(edges), w, std::move(blocks)), _c(std::move(c)), _C0(C0) {
    _rest_initialized = true;
  }

  std::string name() override { return "ring_det_conserve"; }

  /// The ring's `sum c_e det(h_e)` for tensors `h`.
  static real ring_det(const std::vector<vec3> &h, const std::vector<real> &c) {
    real s = 0;
    for (size_t k = 0; k < h.size() && k < c.size(); ++k)
      s += c[k] * (h[k][0] * h[k][2] - h[k][1] * h[k][1]);
    return s;
  }

  /// The projected ring of `h` (the closest with `ring_det = C0`).
  static std::vector<vec3> project_ring(const std::vector<vec3> &h, const std::vector<real> &c,
                                        real C0) {
    real sH = 0, sD = 0;
    std::vector<real> H(h.size());
    for (size_t k = 0; k < h.size(); ++k) {
      H[k] = real(0.5) * (h[k][0] + h[k][2]);
      const real d = real(0.5) * (h[k][0] - h[k][2]);
      sH += c[k] * H[k] * H[k];
      sD += c[k] * (d * d + h[k][1] * h[k][1]);
    }
    const real l = det_multiplier(sH, sD, C0);
    std::vector<vec3> out(h.size());
    for (size_t k = 0; k < h.size(); ++k) {
      const real Hn = H[k] / (real(1) - l);
      const real dn = real(0.5) * (h[k][0] - h[k][2]) / (real(1) + l);
      out[k] = vec3(Hn + dn, h[k][1] / (real(1) + l), Hn - dn);
    }
    return out;
  }

  void project(const vecX &q, vecX &p) override {
    std::vector<vec3> h(_ids.size());
    for (size_t k = 0; k < _ids.size(); ++k)
      h[k] = _blocks[0]->get_vec3(_ids[k], q);
    const std::vector<vec3> g = project_ring(h, _c, _C0);
    for (size_t k = 0; k < _ids.size(); ++k) {
      const real s = _w * std::sqrt(std::max(_c[k], real(0)));
      p[_id0 + static_cast<index_t>(3 * k) + 0] = s * g[k][0];
      p[_id0 + static_cast<index_t>(3 * k) + 1] = s * std::sqrt(real(2)) * g[k][1];
      p[_id0 + static_cast<index_t>(3 * k) + 2] = s * g[k][2];
    }
  }

  void fill_A(index_t &id0, std::vector<trip> &triplets) override {
    _id0 = id0;
    for (size_t k = 0; k < _ids.size(); ++k) {
      const index_t i = _blocks[0]->get_offset_idx(_ids[k]);
      const real s = _w * std::sqrt(std::max(_c[k], real(0)));
      triplets.push_back(trip(id0 + 0, i + 0, s));
      triplets.push_back(trip(id0 + 1, i + 1, s * std::sqrt(real(2))));
      triplets.push_back(trip(id0 + 2, i + 2, s));
      id0 += 3;
    }
  }

  real held() const { return _C0; }

private:
  std::vector<real> _c;
  real _C0 = 0;
};

/// The edges incident on each vertex of the table.
inline std::vector<std::vector<index_t>> vertex_edge_rings(const kusama::stencil_table &tb) {
  std::vector<std::vector<index_t>> rings(static_cast<size_t>(std::max(tb.n_verts, 0)));
  for (int f = 0; f < tb.n_faces; ++f)
    for (int i = 0; i < 3; ++i) {
      const int e = tb.face_edges[static_cast<size_t>(f)][static_cast<size_t>(i)];
      for (int v : {tb.ids[static_cast<size_t>(f)][static_cast<size_t>(i)],
                    tb.ids[static_cast<size_t>(f)][static_cast<size_t>((i + 1) % 3)]})
        if (v >= 0 && v < tb.n_verts && e >= 0)
          rings[static_cast<size_t>(v)].push_back(e);
    }
  for (auto &r : rings) {
    std::sort(r.begin(), r.end());
    r.erase(std::unique(r.begin(), r.end()), r.end());
  }
  return rings;
}

/// Each vertex's chart `(T, B, N)` (`N` the area-weighted face normal, `T`
/// the first spoke projected) and its spokes: edge, weight `a = w_e / W`
/// (`kusama::weingarten::vertex_stars`), and `cos`, `sin` of the projected
/// edge's angle from `T`, so `h_v = sum a R(phi) h_e` is
/// `kusama::weingarten::vertex_tensor` at the step start, linear in the store.
struct vertex_star_charts {
  struct spoke {
    index_t edge = -1;
    real a = 0, c = 1, s = 0;
  };
  std::vector<std::vector<spoke>> at;
  std::vector<vec3> T, B, N;
};

/// `h` from a spoke's frame into one turned by `(c, s)`: the spoke tangent
/// is `c T + s B`.
inline vec3 turn_tensor(const vec3 &h, real c, real s) {
  return vec3(h[0] * c * c - real(2) * h[1] * c * s + h[2] * s * s,
              (h[0] - h[2]) * c * s + h[1] * (c * c - s * s),
              h[0] * s * s + real(2) * h[1] * c * s + h[2] * c * c);
}

inline vertex_star_charts build_vertex_star_charts(const kusama::stencil_table &tb,
                                                   const std::vector<vec3> &x) {
  vertex_star_charts out;
  const size_t nv = static_cast<size_t>(std::max(tb.n_verts, 0));
  out.at.resize(nv);
  out.T.assign(nv, vec3::Zero());
  out.B.assign(nv, vec3::Zero());
  out.N.assign(nv, vec3::Zero());
  if (x.size() < nv)
    return out;
  for (int f = 0; f < tb.n_faces; ++f) {
    const auto &id = tb.ids[static_cast<size_t>(f)];
    const vec3 n = (x[static_cast<size_t>(id[1])] - x[static_cast<size_t>(id[0])])
                       .cross(x[static_cast<size_t>(id[2])] - x[static_cast<size_t>(id[0])]);
    for (int i = 0; i < 3; ++i)
      out.N[static_cast<size_t>(id[static_cast<size_t>(i)])] += n;
  }
  const kusama::weingarten::vertex_stars stars = kusama::weingarten::build_vertex_stars(tb, x);
  for (size_t v = 0; v < nv; ++v) {
    const real nn = out.N[v].norm();
    if (!(nn > real(0)))
      continue;
    const vec3 N = out.N[v] / nn;
    out.N[v] = N;
    real W = 0;
    for (const auto &sp : stars.at[v])
      W += sp.w;
    if (!(W > real(0)))
      continue;
    for (const auto &sp : stars.at[v]) {
      vec3 t = x[static_cast<size_t>(sp.hi)] - x[static_cast<size_t>(sp.lo)];
      t -= N * t.dot(N);
      const real l = t.norm();
      if (!(l > real(1e-20)) || !(sp.w > real(0)))
        continue;
      t /= l;
      if (out.T[v].squaredNorm() == real(0)) {
        out.T[v] = t;
        out.B[v] = N.cross(t);
      }
      out.at[v].push_back({static_cast<index_t>(sp.edge), sp.w / W, t.dot(out.T[v]),
                           t.dot(out.B[v])});
    }
  }
  return out;
}

/// The vertex tensors of store `h` in the charts (NaN where a vertex has no
/// spokes).
inline std::vector<vec3> vertex_tensors(const vertex_star_charts &ch, const std::vector<vec3> &h) {
  std::vector<vec3> out(ch.at.size(), vec3::Constant(std::numeric_limits<real>::quiet_NaN()));
  for (size_t v = 0; v < ch.at.size(); ++v) {
    if (ch.at[v].empty())
      continue;
    vec3 hv = vec3::Zero();
    for (const auto &sp : ch.at[v])
      hv += sp.a * turn_tensor(h[static_cast<size_t>(sp.edge)], sp.c, sp.s);
    out[v] = hv;
  }
  return out;
}

/// Each vertex's Gaussian curvature `det h_v` of store `h`.
inline std::vector<real> vertex_gauss(const vertex_star_charts &ch, const std::vector<vec3> &h) {
  std::vector<real> K;
  K.reserve(ch.at.size());
  for (const vec3 &hv : vertex_tensors(ch, h))
    K.push_back(hv[0] * hv[2] - hv[1] * hv[1]);
  return K;
}

/// One vertex's Gaussian-curvature row on the edge store (`blocks[0]`): its
/// spoke tensors `h_e` with `det h_v = K0`, `h_v = sum a R(phi) h_e`. Rows
/// `w sqrt(c_e) (1, sqrt2, 1)`, `c_e` the coupling factor. The projection is
/// the closest spokes in that metric: turns keep the tensor metric, so the
/// vertex moves to `closest_det_tensor(h_v, K0)` and each spoke takes the
/// least-norm share `a / (c_e sigma) R(-phi) (g - h_v)`,
/// `sigma = sum a^2 / c_e`. Mean curvature and anisotropy are free; `K0` is
/// fixed for the step, so the row is variational.
class vertex_gauss_conserve : public block_constraint {
public:
  DEFINE_CREATE_FUNC(vertex_gauss_conserve)

  vertex_gauss_conserve(std::vector<vertex_star_charts::spoke> spokes, std::vector<real> c, real K0,
                        real w, std::vector<sim_block::ptr> blocks)
      : block_constraint(edges_of(spokes), w, std::move(blocks)), _spokes(std::move(spokes)),
        _c(std::move(c)), _K0(K0) {
    _rest_initialized = true;
  }

  std::string name() override { return "vertex_gauss_conserve"; }

  /// The vertex tensor of the spokes `h`.
  vec3 vertex_tensor(const std::vector<vec3> &h) const {
    vec3 hv = vec3::Zero();
    for (size_t k = 0; k < _spokes.size(); ++k)
      hv += _spokes[k].a * turn_tensor(h[k], _spokes[k].c, _spokes[k].s);
    return hv;
  }

  /// The projected spokes of `h`.
  std::vector<vec3> project_spokes(const std::vector<vec3> &h) const {
    real sigma = 0;
    for (size_t k = 0; k < _spokes.size(); ++k)
      sigma += _spokes[k].a * _spokes[k].a / _c[k];
    std::vector<vec3> out = h;
    if (!(sigma > real(0)))
      return out;
    const vec3 hv = vertex_tensor(h);
    const vec3 d = closest_det_tensor(hv, _K0) - hv;
    for (size_t k = 0; k < _spokes.size(); ++k)
      out[k] += (_spokes[k].a / (_c[k] * sigma)) * turn_tensor(d, _spokes[k].c, -_spokes[k].s);
    return out;
  }

  void project(const vecX &q, vecX &p) override {
    std::vector<vec3> h(_ids.size());
    for (size_t k = 0; k < _ids.size(); ++k)
      h[k] = _blocks[0]->get_vec3(_ids[k], q);
    const std::vector<vec3> g = project_spokes(h);
    for (size_t k = 0; k < _ids.size(); ++k) {
      const real s = _w * std::sqrt(_c[k]);
      p[_id0 + static_cast<index_t>(3 * k) + 0] = s * g[k][0];
      p[_id0 + static_cast<index_t>(3 * k) + 1] = s * std::sqrt(real(2)) * g[k][1];
      p[_id0 + static_cast<index_t>(3 * k) + 2] = s * g[k][2];
    }
  }

  void fill_A(index_t &id0, std::vector<trip> &triplets) override {
    _id0 = id0;
    for (size_t k = 0; k < _ids.size(); ++k) {
      const index_t i = _blocks[0]->get_offset_idx(_ids[k]);
      const real s = _w * std::sqrt(_c[k]);
      triplets.push_back(trip(id0 + 0, i + 0, s));
      triplets.push_back(trip(id0 + 1, i + 1, s * std::sqrt(real(2))));
      triplets.push_back(trip(id0 + 2, i + 2, s));
      id0 += 3;
    }
  }

  real held() const { return _K0; }

private:
  static std::vector<index_t> edges_of(const std::vector<vertex_star_charts::spoke> &sp) {
    std::vector<index_t> e;
    e.reserve(sp.size());
    for (const auto &s : sp)
      e.push_back(s.edge);
    return e;
  }

  std::vector<vertex_star_charts::spoke> _spokes;
  std::vector<real> _c;
  real _K0 = 0;
};

/// One `vertex_gauss_conserve` per vertex with a finite `K0[v]`, over its
/// spokes with positive coupling. `w` 0 leaves the rows off.
inline void init_vertex_gauss(const kusama::stencil_table &tb, const std::vector<vec3> &x,
                              const vertex_star_charts &ch, const std::vector<real> &K0, real w,
                              sim_block::ptr store,
                              std::vector<projection_constraint::ptr> &constraints) {
  if (!(w > real(0)) || !std::isfinite(w) || !store || K0.size() < ch.at.size())
    return;
  const std::vector<real> c = kusama::weingarten::coupling_factor(tb, x);
  for (size_t v = 0; v < ch.at.size(); ++v) {
    if (ch.at[v].empty() || !std::isfinite(K0[v]))
      continue;
    std::vector<vertex_star_charts::spoke> sp;
    std::vector<real> cs;
    for (const auto &s : ch.at[v]) {
      const real ce = c[static_cast<size_t>(s.edge)];
      if (!(ce > real(0)))
        continue;
      sp.push_back(s);
      cs.push_back(ce);
    }
    if (!sp.empty())
      constraints.push_back(vertex_gauss_conserve::create(std::move(sp), std::move(cs), K0[v], w,
                                                          std::vector<sim_block::ptr>{store}));
  }
}

/// One `ring_det_conserve` per vertex, holding its ring's `sum c_e det(h0_e)`
/// from the frame-start store `h0`. `w` 0 leaves the rows off.
inline void init_ring_det(const kusama::stencil_table &tb, const std::vector<vec3> &x,
                          const std::vector<vec3> &h0, real w, sim_block::ptr store,
                          std::vector<projection_constraint::ptr> &constraints) {
  if (!(w > real(0)) || !std::isfinite(w) || !store ||
      static_cast<int>(h0.size()) != tb.n_edges)
    return;
  const std::vector<real> c = kusama::weingarten::coupling_factor(tb, x);
  for (const std::vector<index_t> &ring : vertex_edge_rings(tb)) {
    if (ring.empty())
      continue;
    std::vector<real> cr(ring.size());
    std::vector<vec3> hr(ring.size());
    for (size_t k = 0; k < ring.size(); ++k) {
      cr[k] = c[static_cast<size_t>(ring[k])];
      hr[k] = h0[static_cast<size_t>(ring[k])];
    }
    constraints.push_back(ring_det_conserve::create(
        ring, cr, ring_det_conserve::ring_det(hr, cr), w, std::vector<sim_block::ptr>{store}));
  }
}

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif
