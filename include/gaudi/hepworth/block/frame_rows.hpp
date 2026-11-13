#ifndef GAUDI_HEPWORTH_BLOCK_FRAME_ROWS_HPP
#define GAUDI_HEPWORTH_BLOCK_FRAME_ROWS_HPP

/// Native curvature rows on positions and per-vertex quaternion frames
/// (scripts/sympy/cr_poisson/MATH.md 19.21, E3; prototype `e3_native.py`).
///
/// Blocks: `blocks[0]` positions, `blocks[1]` the frames `Y`
/// (`shell_cr_block::frame_dof`), as the contact net's frames mode. The frame
/// gives `(T, B, N) = Y (e1, e2, e3) conj(Y) / |Y|^2`. Each vertex tensor
/// `S_v = (s11, s12, s22)` in `v`'s own `(T, B)` is the least-squares fit of
/// the spokes' Rodrigues normal curvature `k = t^T S t` and twist
/// `tau = t^T S (N x t)`, `t` the spoke tangent, `k = -dn . t / |dx_t|`,
/// `tau = -dn . (nbar x t) / |dx_t|`. Rows that compare vertices read `S`
/// along frozen world directions projected into each tangent plane, or carry
/// `z = (s11 - s22) / 2 + i s12` by the frames' relative turn, so the turn of
/// `Y` about `N` drops out.
///
/// Rows (weight 0 leaves a row out), each scaled `w k_frame_rows A`, `A` the
/// vertex (or edge) dual area, so the residuals are integrated curvatures:
/// - memory: `z_v` held at its frame-start value turned by the ring's best
///   common turn (AMAP-style shape memory; the turn frozen per iterate), and
///   `H_v` at its frame start's with `memory_mean`;
/// - tension: `H_v` to a per-frame constant (`frame_tension`);
/// - minimizer: per directed spoke, `k1` and `k2` (along the frame start's
///   principal directions of the spoke's foot) differenced across their own
///   lines (across-line Dirichlet);
/// - dupin: the same differences along their lines; with `dupin_reach > 0`
///   (and the vertex-tensor block) instead the traced-line form: from each
///   face's centroid both principal lines are walked `dupin_reach` mean edge
///   lengths each way (one mean edge length per step, snapped to the reached
///   face's principal direction), each sample's normal curvature is the
///   barycentric blend of its face's vertex tensors along the line, and the
///   rows hold the residual of a degree-`dupin_degree` fit in arc length.
///   Traced at construction and frozen (samples, weights, projections), so
///   the rows are linear in the tensors;
/// - growth: the normal curvature along `growth_dir` times `1 + x`, across it
///   times `1 + y`, the torsion along it held;
/// - det: each vertex's sum of its faces' normal solid angles (the Gauss
///   image area) at its frame-start value (scaled by the mean edge length);
/// - gradtl: per edge, `z_i - z_j e^{2 i phi}`, `phi` the frames' relative
///   turn;
/// - twist: the geodesic torsion along `twist_dir` at its frame-start value
///   plus `tau0`.
/// Jacobians by central differences over each row's local unknowns.
///
/// Every row reads the vertex tensors as functions of the ring's positions
/// and frames (the fit), with the Jacobian chained through it (`dS` per
/// vertex, cached once per linearization): the rows act on the frames
/// directly and the contact net resolves frames against positions. (A tensor
/// block tied to the frames by per-spoke fit rows was tried and removed: the
/// fit rows, about 12 per vertex against 3 unknowns, could not be satisfied
/// and damped every ring ripple, so any frame row removed all buckling;
/// MATH.md 19.21.)
///
/// - target: each directed spoke's Rodrigues normal curvature and twist from
///   the frames pinned to given values (`set_spoke_targets`, the target
///   mesh's readings), in the chord form (times the live tangential length,
///   over the frame-start one); edge-native, so a row touches the spoke's
///   two ends.

#include "gaudi/hepworth/block/block_constraint.hpp"
#include "gaudi/hepworth/block/quat_util.hpp"
#include "gaudi/hepworth/block/surface_walk.hpp"
#include "gaudi/kusama/contact_net.hpp"
#include "gaudi/kusama/stencil_table.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <functional>
#include <limits>
#include <vector>

namespace gaudi {
namespace hepworth {
namespace block {

/// Row scale, as `k_contact_gn` for the contact net's rows.
constexpr real k_frame_rows = real(14.2);

enum class frame_tension { constant, frame_mean, volume };

struct frame_rows_config {
  real w_memory = 0;
  bool memory_mean = true;
  real w_tension = 0;
  frame_tension tension = frame_tension::frame_mean;
  /// `constant` target (outward normals: convex `H < 0`).
  real h0 = 0;
  /// `volume`: the frame-mean target flattened by `volume_gain` times the
  /// relative volume deficit against `volume_ref` (NaN: no correction).
  real volume_gain = 1;
  real volume_ref = std::numeric_limits<real>::quiet_NaN();
  real w_minimizer = 0;
  real w_dupin = 0;
  /// Traced-line Dupin rows: reach each way in mean edge lengths (0: the
  /// per-spoke rows) and the fit degree along the line.
  real dupin_reach = 0;
  int dupin_degree = 0;
  real w_growth = 0;
  real growth_x = 0, growth_y = 0;
  real w_det = 0;
  real w_gradtl = 0;
  real w_twist = 0;
  real tau0 = 0;
  /// Spoke targets (`set_spoke_targets`).
  real w_target = 0;
};

namespace frame_detail {

struct frame3 {
  vec3 T, B, N;
};

inline frame3 frame_of(const vec4 &Y) {
  const real l = Y.norm();
  const vec4 q = l > real(0) ? vec4(Y / l) : vec4(1, 0, 0, 0);
  const real w = q[0], x = q[1], y = q[2], z = q[3];
  return {vec3(1 - 2 * (y * y + z * z), 2 * (x * y + w * z), 2 * (x * z - w * y)),
          vec3(2 * (x * y - w * z), 1 - 2 * (x * x + z * z), 2 * (y * z + w * x)),
          vec3(2 * (x * z + w * y), 2 * (y * z - w * x), 1 - 2 * (x * x + y * y))};
}

/// Normal curvature and geodesic torsion of `s` along `w` projected into `f`.
inline void along(const vec3 &s, const frame3 &f, const vec3 &w, real &k, real &tau) {
  real a = w.dot(f.T), b = w.dot(f.B);
  const real r = std::sqrt(a * a + b * b);
  if (!(r > real(0))) {
    k = tau = 0;
    return;
  }
  a /= r;
  b /= r;
  k = a * a * s[0] + 2 * a * b * s[1] + b * b * s[2];
  tau = -a * b * s[0] + (a * a - b * b) * s[1] + a * b * s[2];
}

inline std::complex<real> z_of(const vec3 &s) { return {real(0.5) * (s[0] - s[2]), s[1]}; }
inline real H_of(const vec3 &s) { return real(0.5) * (s[0] + s[2]); }

/// Angle of `fj.T`, carried onto `fi`'s plane by the shortest arc, in `fi`.
inline real relative_turn(const frame3 &fi, const frame3 &fj) {
  const vec3 c = fj.N.cross(fi.N);
  const real d = fj.N.dot(fi.N);
  const vec3 t = fj.T * d + c.cross(fj.T) + c * c.dot(fj.T) / std::max(real(1) + d, real(1e-12));
  return std::atan2(t.dot(fi.B), t.dot(fi.T));
}

inline real solid_angle(const vec3 &a, const vec3 &b, const vec3 &c) {
  return real(2) * std::atan2(a.dot(b.cross(c)), real(1) + a.dot(b) + b.dot(c) + c.dot(a));
}

/// `S_v = (s11, s12, s22)` in `F[v]`'s `(T, B)`: the least-squares fit of the
/// spokes' Rodrigues normal curvature and twist from positions and frames.
inline vec3 fit_vertex_tensor(const kusama::contact_net::topology &tp, int v, const std::vector<vec3> &x,
                              const std::vector<frame3> &F) {
  Eigen::Matrix<real, 3, 3> M = Eigen::Matrix<real, 3, 3>::Identity() * real(1e-12);
  vec3 rhs = vec3::Zero();
  const size_t vi = static_cast<size_t>(v);
  for (int s = tp.spoke_begin[vi]; s < tp.spoke_begin[vi + 1]; ++s) {
    const size_t j = static_cast<size_t>(tp.spokes[static_cast<size_t>(s)][1]);
    const vec3 dx = x[j] - x[vi], dn = F[j].N - F[vi].N;
    const vec3 nb = (F[vi].N + F[j].N).normalized();
    const vec3 dxt = dx - nb * dx.dot(nb);
    const real L = dxt.norm();
    if (!(L > real(0)))
      continue;
    const vec3 t = dxt / L, b = nb.cross(t);
    const real k = -dn.dot(t) / L, tau = -dn.dot(b) / L;
    real a = t.dot(F[vi].T), c = t.dot(F[vi].B);
    const real r = std::sqrt(a * a + c * c);
    if (!(r > real(0)))
      continue;
    a /= r;
    c /= r;
    const vec3 r1(a * a, 2 * a * c, c * c), r2(-a * c, a * a - c * c, a * c);
    M += r1 * r1.transpose() + r2 * r2.transpose();
    rhs += r1 * k + r2 * tau;
  }
  return M.ldlt().solve(rhs);
}

} // namespace frame_detail

class frame_rows : public block_constraint {
public:
  DEFINE_CREATE_FUNC(frame_rows)

  /// `x0`, `Y0`: the frame start. `growth_dir`, `twist_dir`: per vertex world
  /// directions (empty: the frame start's tight principal direction).
  frame_rows(const kusama::stencil_table *table, const std::vector<vec3> &x0, const std::vector<vec4> &Y0,
             frame_rows_config cfg, std::vector<sim_block::ptr> blocks, std::vector<real> growth_angle = {},
             std::vector<real> twist_angle = {})
      : block_constraint({}, real(1), std::move(blocks)), _t(table), _cfg(cfg) {
    namespace fd = frame_detail;
    _tp = kusama::contact_net::build_topology(*_t);
    const int n = _tp.n;
    const size_t nv = static_cast<size_t>(n);
    _faces_of.assign(nv, {});
    for (size_t f = 0; f < _tp.tris.size(); ++f)
      for (int k = 0; k < 3; ++k)
        _faces_of[static_cast<size_t>(_tp.tris[f][static_cast<size_t>(k)])].push_back(static_cast<int>(f));
    _A.assign(nv, 0);
    real len = 0;
    int nl = 0;
    for (const auto &tr : _tp.tris) {
      const vec3 &a = x0[static_cast<size_t>(tr[0])], &b = x0[static_cast<size_t>(tr[1])],
                 &c = x0[static_cast<size_t>(tr[2])];
      const real ar = real(0.5) * (b - a).cross(c - a).norm();
      for (int k = 0; k < 3; ++k)
        _A[static_cast<size_t>(tr[static_cast<size_t>(k)])] += ar / real(3);
    }
    for (const auto &e : _tp.edges)
      if (e[0] >= 0) {
        len += (x0[static_cast<size_t>(e[1])] - x0[static_cast<size_t>(e[0])]).norm();
        ++nl;
      }
    _hl = nl > 0 ? len / nl : real(1);
    _F0.resize(nv);
    for (size_t v = 0; v < nv; ++v)
      _F0[v] = fd::frame_of(Y0[v]);
    std::vector<vec3> N0(nv);
    for (size_t v = 0; v < nv; ++v)
      N0[v] = _F0[v].N;
    _S0.resize(nv);
    for (int v = 0; v < n; ++v)
      _S0[static_cast<size_t>(v)] = fit(v, x0, _F0);
    // Principal directions (world), tight (larger |k|) first.
    _e1.resize(nv);
    _e2.resize(nv);
    for (size_t v = 0; v < nv; ++v) {
      const vec3 &s = _S0[v];
      const real th = real(0.5) * std::atan2(2 * s[1], s[0] - s[2]);
      vec3 e1 = std::cos(th) * _F0[v].T + std::sin(th) * _F0[v].B;
      vec3 e2 = _F0[v].N.cross(e1);
      real k1, k2, t;
      fd::along(s, _F0[v], e1, k1, t);
      fd::along(s, _F0[v], e2, k2, t);
      if (std::abs(k2) > std::abs(k1)) {
        const vec3 tmp = e1;
        e1 = e2;
        e2 = -tmp;
      }
      _e1[v] = e1;
      _e2[v] = e2;
    }
    auto dir = [&](const std::vector<real> &ang, size_t v) {
      const real a = v < ang.size() ? ang[v] : real(0);
      return vec3(std::cos(a) * _e1[v] + std::sin(a) * _e2[v]);
    };
    _ga.resize(nv);
    _gk.resize(nv);
    _tw.resize(nv);
    _tw0.resize(nv);
    _z0.resize(nv);
    _H0.resize(nv);
    for (size_t v = 0; v < nv; ++v) {
      _z0[v] = fd::z_of(_S0[v]);
      _H0[v] = fd::H_of(_S0[v]);
      _ga[v] = dir(growth_angle, v);
      real ka, ta, kb, tb;
      fd::along(_S0[v], _F0[v], _ga[v], ka, ta);
      fd::along(_S0[v], _F0[v], _F0[v].N.cross(_ga[v]), kb, tb);
      _gk[v] = vec3((1 + _cfg.growth_x) * ka, (1 + _cfg.growth_y) * kb, ta);
      _tw[v] = dir(twist_angle, v);
      real kt;
      fd::along(_S0[v], _F0[v], _tw[v], kt, _tw0[v]);
    }
    // Memory: the neighbours' frame-start z carried into each foot's frame.
    _zring0.resize(_tp.spokes.size());
    for (size_t s = 0; s < _tp.spokes.size(); ++s) {
      const size_t i = static_cast<size_t>(_tp.spokes[s][0]), j = static_cast<size_t>(_tp.spokes[s][1]);
      _zring0[s] = _z0[j] * std::polar(real(1), real(2) * fd::relative_turn(_F0[i], _F0[j]));
    }
    // Tension target.
    real Hm = 0, At = 0;
    for (size_t v = 0; v < nv; ++v) {
      Hm += _A[v] * _H0[v];
      At += _A[v];
    }
    Hm /= std::max(At, real(1e-30));
    _Ht = Hm;
    if (_cfg.tension == frame_tension::constant) {
      _Ht = _cfg.h0;
    } else if (_cfg.tension == frame_tension::volume && std::isfinite(_cfg.volume_ref)) {
      real V = 0;
      for (const auto &tr : _tp.tris)
        V += x0[static_cast<size_t>(tr[0])].dot(
                 x0[static_cast<size_t>(tr[1])].cross(x0[static_cast<size_t>(tr[2])])) /
             real(6);
      _Ht = Hm + _cfg.volume_gain * (_cfg.volume_ref - V) / std::abs(_cfg.volume_ref) * std::abs(Hm);
    }
    _Om0.assign(nv, 0);
    for (int v = 0; v < n; ++v)
      _Om0[static_cast<size_t>(v)] = gauss_area(v, N0);
    if (_cfg.w_dupin > 0 && _cfg.dupin_reach > 0)
      build_dupin_lines(x0);
    _spoke_l0.assign(_tp.spokes.size(), real(0));
    for (size_t sp = 0; sp < _tp.spokes.size(); ++sp) {
      const size_t i = static_cast<size_t>(_tp.spokes[sp][0]), j = static_cast<size_t>(_tp.spokes[sp][1]);
      const vec3 dx = x0[j] - x0[i], nb = (_F0[i].N + _F0[j].N).normalized();
      _spoke_l0[sp] = (dx - nb * dx.dot(nb)).norm();
    }
    _rest_initialized = true;
  }

  std::string name() override { return "frame_rows"; }
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

  /// The frame-start vertex tensor in the frame-start frame.
  const std::vector<vec3> &start_tensors() const { return _S0; }
  /// Per directed spoke (`kusama::contact_net::build_topology` order) the
  /// target normal curvature and twist for the `w_target` rows.
  void set_spoke_targets(std::vector<real> k, std::vector<real> tau) {
    _kt = std::move(k);
    _taut = std::move(tau);
  }
  /// Directed spoke `(k, tau)` readings of positions and frames, in the
  /// topology's spoke order (`n` the frames' normals).
  static void spoke_readings(const kusama::contact_net::topology &tp, const std::vector<vec3> &x,
                             const std::vector<vec3> &n, std::vector<real> &k, std::vector<real> &tau) {
    k.assign(tp.spokes.size(), real(0));
    tau.assign(tp.spokes.size(), real(0));
    for (size_t s = 0; s < tp.spokes.size(); ++s) {
      const size_t i = static_cast<size_t>(tp.spokes[s][0]), j = static_cast<size_t>(tp.spokes[s][1]);
      read_spoke(x[i], x[j], n[i], n[j], k[s], tau[s]);
    }
  }
  static void read_spoke(const vec3 &xi, const vec3 &xj, const vec3 &ni, const vec3 &nj, real &k, real &tau) {
    const vec3 dx = xj - xi, dn = nj - ni;
    const vec3 nb = (ni + nj).normalized();
    const vec3 dxt = dx - nb * dx.dot(nb);
    const real L = dxt.norm();
    k = tau = 0;
    if (!(L > real(0)))
      return;
    const vec3 t = dxt / L, b = nb.cross(t);
    k = -dn.dot(t) / L;
    tau = -dn.dot(b) / L;
  }
  /// The tension target in use this frame.
  real tension_target() const { return _Ht; }

private:
  using frame3 = frame_detail::frame3;

  /// `S_v` from positions and frames.
  vec3 fit(int v, const std::vector<vec3> &x, const std::vector<frame3> &F) const {
    return frame_detail::fit_vertex_tensor(_tp, v, x, F);
  }

  real gauss_area(int v, const std::vector<vec3> &N) const {
    real om = 0;
    for (int f : _faces_of[static_cast<size_t>(v)]) {
      const auto &tr = _tp.tris[static_cast<size_t>(f)];
      om += frame_detail::solid_angle(N[static_cast<size_t>(tr[0])], N[static_cast<size_t>(tr[1])],
                                      N[static_cast<size_t>(tr[2])]);
    }
    return om;
  }

  /// A traced Dupin line: the vertices whose tensors it reads and the rows
  /// `C S` (stacked `(s11, s12, s22)` per vertex, `w_dupin` applied later).
  struct dupin_line {
    std::vector<int> verts;
    matX C;
  };

  /// The face's principal directions from its corners' frame-start tensors
  /// (world tensors averaged, read in the face plane).
  bool face_dirs(const std::vector<vec3> &x, int f, vec3 &d1, vec3 &d2) const {
    const auto &id = _t->ids[static_cast<size_t>(f)];
    const vec3 a = x[static_cast<size_t>(id[0])], b = x[static_cast<size_t>(id[1])],
               c = x[static_cast<size_t>(id[2])];
    vec3 nf = (b - a).cross(c - a);
    if (!(nf.norm() > real(0)))
      return false;
    nf.normalize();
    mat3 W = mat3::Zero();
    for (int k = 0; k < 3; ++k) {
      const size_t u = static_cast<size_t>(id[static_cast<size_t>(k)]);
      const vec3 &S = _S0[u];
      const vec3 &T = _F0[u].T, &B = _F0[u].B;
      W += S[0] * T * T.transpose() + S[1] * (T * B.transpose() + B * T.transpose()) + S[2] * B * B.transpose();
    }
    const vec3 t1 = (b - a).normalized(), t2 = nf.cross(t1);
    const real p = t1.dot(W * t1), q = t1.dot(W * t2), r = t2.dot(W * t2);
    const real th = real(0.5) * std::atan2(2 * q, p - r);
    d1 = std::cos(th) * t1 + std::sin(th) * t2;
    d2 = nf.cross(d1);
    return d1.allFinite() && d2.allFinite();
  }

  void build_dupin_lines(const std::vector<vec3> &x0) {
    const kusama::stencil_table &t = *_t;
    const int steps = std::max(1, static_cast<int>(std::lround(_cfg.dupin_reach)));
    const int deg = std::max(0, _cfg.dupin_degree);
    const real step = _hl;
    auto walk = [&](int f, const vec3 &p0, const vec3 &d0) {
      std::vector<walk_sample> out;
      int g = f;
      vec3 p = p0, d = d0;
      for (int k = 0; k < steps; ++k) {
        const std::vector<walk_sample> w = trace_on_surface(t, x0, g, p, d, step, step);
        if (w.empty())
          break;
        walk_sample ws = w.back();
        vec3 e1, e2;
        if (face_dirs(x0, ws.face, e1, e2)) {
          const vec3 &best = std::abs(e1.dot(ws.dir)) >= std::abs(e2.dot(ws.dir)) ? e1 : e2;
          ws.dir = best.dot(ws.dir) >= 0 ? best : vec3(-best);
        }
        ws.s = real(k + 1) * step;
        out.push_back(ws);
        g = ws.face;
        p = ws.p;
        d = ws.dir;
      }
      return out;
    };
    for (int f = 0; f < t.n_faces; ++f) {
      vec3 d1, d2;
      if (!face_dirs(x0, f, d1, d2))
        continue;
      const auto &id = t.ids[static_cast<size_t>(f)];
      const vec3 c =
          (x0[static_cast<size_t>(id[0])] + x0[static_cast<size_t>(id[1])] + x0[static_cast<size_t>(id[2])]) / real(3);
      const real Af = real(0.5) * (x0[static_cast<size_t>(id[1])] - x0[static_cast<size_t>(id[0])])
                                      .cross(x0[static_cast<size_t>(id[2])] - x0[static_cast<size_t>(id[0])])
                                      .norm();
      for (const vec3 &d : {d1, d2}) {
        const std::vector<walk_sample> fwd = walk(f, c, d), bwd = walk(f, c, -d);
        std::vector<walk_sample> line;
        for (auto it = bwd.rbegin(); it != bwd.rend(); ++it) {
          walk_sample w = *it;
          w.s = -w.s;
          line.push_back(w);
        }
        walk_sample w0;
        w0.face = f;
        w0.bary = vec3::Constant(real(1) / real(3));
        w0.p = c;
        w0.dir = d;
        w0.s = 0;
        line.push_back(w0);
        line.insert(line.end(), fwd.begin(), fwd.end());
        const int m = static_cast<int>(line.size());
        if (m < deg + 2)
          continue;
        dupin_line dl;
        std::vector<int> slot(static_cast<size_t>(_tp.n), -1);
        for (const walk_sample &w : line)
          for (int k = 0; k < 3; ++k) {
            const int u = t.ids[static_cast<size_t>(w.face)][static_cast<size_t>(k)];
            if (slot[static_cast<size_t>(u)] < 0) {
              slot[static_cast<size_t>(u)] = static_cast<int>(dl.verts.size());
              dl.verts.push_back(u);
            }
          }
        const Eigen::Index nv3 = 3 * static_cast<Eigen::Index>(dl.verts.size());
        matX K = matX::Zero(m, nv3);
        matX V(m, deg + 1);
        const real sc = std::max(_cfg.dupin_reach * step, real(1e-12));
        for (int i = 0; i < m; ++i) {
          const walk_sample &w = line[static_cast<size_t>(i)];
          for (int k = 0; k < 3; ++k) {
            const size_t u = static_cast<size_t>(t.ids[static_cast<size_t>(w.face)][static_cast<size_t>(k)]);
            real a = w.dir.dot(_F0[u].T), b = w.dir.dot(_F0[u].B);
            const real r = std::sqrt(a * a + b * b);
            if (!(r > real(0)))
              continue;
            a /= r;
            b /= r;
            const Eigen::Index col = 3 * slot[u];
            K(i, col + 0) += w.bary[k] * a * a;
            K(i, col + 1) += w.bary[k] * 2 * a * b;
            K(i, col + 2) += w.bary[k] * b * b;
          }
          for (int pw = 0; pw <= deg; ++pw)
            V(i, pw) = std::pow(w.s / sc, real(pw));
        }
        const Eigen::HouseholderQR<matX> qr(V);
        const matX Q = qr.householderQ() * matX::Identity(m, deg + 1);
        const matX P = matX::Identity(m, m) - Q * Q.transpose();
        dl.C = (k_frame_rows * Af / std::sqrt(real(m))) * P * K;
        _lines.push_back(std::move(dl));
      }
    }
  }

  /// One row unit: `m` residuals from the cached vertex tensors `Sc` of
  /// `sv` and the live positions and frames of `dv`.
  struct unit {
    int m = 0;
    std::vector<int> sv, dv;
    std::function<void(real *)> f;
  };

  void ring_of(int v, std::vector<int> &out) const {
    out.push_back(v);
    const size_t vi = static_cast<size_t>(v);
    for (int s = _tp.spoke_begin[vi]; s < _tp.spoke_begin[vi + 1]; ++s)
      out.push_back(_tp.spokes[static_cast<size_t>(s)][1]);
  }

  index_t col_of(int v, int comp) const {
    return comp < 3 ? _blocks[0]->get_offset_idx(v) + comp : _blocks[1]->get_offset_idx(v) + comp - 3;
  }

  void linearize(const vecX &q, vecX &r, std::vector<trip> *jt, index_t id0) const {
    namespace fd = frame_detail;
    const int n = _tp.n;
    const size_t nv = static_cast<size_t>(n);
    std::vector<vec3> X(nv);
    std::vector<vec4> Y(nv);
    std::vector<frame3> F(nv);
    for (int v = 0; v < n; ++v) {
      X[static_cast<size_t>(v)] = _blocks[0]->get_vec3(v, q);
      Y[static_cast<size_t>(v)] = get_elem<vec4>(v, q, _blocks[1]->_offset);
      F[static_cast<size_t>(v)] = fd::frame_of(Y[static_cast<size_t>(v)]);
    }
    const auto &c = _cfg;
    const bool need_s = c.w_memory > 0 || c.w_tension > 0 || c.w_growth > 0 || c.w_twist > 0 || c.w_minimizer > 0 ||
                        c.w_dupin > 0 || c.w_gradtl > 0;
    // Vertex tensors and, with a Jacobian, their derivatives in the ring's
    // positions and frames (`dS[v]` is 3 x 7 |ring|, columns in `ring_of` order).
    std::vector<vec3> Sc(nv, vec3::Zero());
    std::vector<Eigen::Matrix<real, 3, Eigen::Dynamic>> dS(jt && need_s ? nv : 0);
    if (need_s)
      for (int v = 0; v < n; ++v) {
        const size_t vi = static_cast<size_t>(v);
        Sc[vi] = fit(v, X, F);
        if (!jt)
          continue;
        std::vector<int> ring;
        ring_of(v, ring);
        dS[vi].resize(3, 7 * static_cast<Eigen::Index>(ring.size()));
        for (size_t a = 0; a < ring.size(); ++a) {
          const size_t u = static_cast<size_t>(ring[a]);
          for (int comp = 0; comp < 7; ++comp) {
            const bool pos = comp < 3;
            real &slot = pos ? X[u][comp] : Y[u][comp - 3];
            const real keep = slot, h = real(1e-7) * std::max(real(1), std::abs(keep));
            slot = keep + h;
            if (!pos)
              F[u] = fd::frame_of(Y[u]);
            const vec3 sp = fit(v, X, F);
            slot = keep - h;
            if (!pos)
              F[u] = fd::frame_of(Y[u]);
            const vec3 sm = fit(v, X, F);
            slot = keep;
            if (!pos)
              F[u] = fd::frame_of(Y[u]);
            dS[vi].col(7 * static_cast<Eigen::Index>(a) + comp) = (sp - sm) / (real(2) * h);
          }
        }
      }
    const real k = k_frame_rows;
    std::vector<unit> units;
    std::vector<std::complex<real>> rot;
    if (c.w_memory > 0) {
      // The ring's common turn from the live state, frozen for the linearization.
      rot.assign(nv, {1, 0});
      for (int v = 0; v < n; ++v) {
        const size_t vi = static_cast<size_t>(v);
        std::complex<real> acc = fd::z_of(Sc[vi]) * std::conj(_z0[vi]);
        for (int s = _tp.spoke_begin[vi]; s < _tp.spoke_begin[vi + 1]; ++s) {
          const size_t j = static_cast<size_t>(_tp.spokes[static_cast<size_t>(s)][1]);
          acc += fd::z_of(Sc[j]) * std::polar(real(1), real(2) * fd::relative_turn(F[vi], F[j])) *
                 std::conj(_zring0[static_cast<size_t>(s)]);
        }
        rot[vi] = std::abs(acc) > real(1e-30) ? acc / std::abs(acc) : std::complex<real>(1, 0);
      }
    }
    for (int v = 0; v < n; ++v) {
      const size_t vi = static_cast<size_t>(v);
      const real sv = k * _A[vi];
      if (c.w_memory > 0) {
        const std::complex<real> g = rot[vi] * _z0[vi];
        const bool mean = c.memory_mean;
        const real s = c.w_memory * sv;
        units.push_back({mean ? 3 : 2, {v}, {}, [&, g, s, mean, vi](real *o) {
                           const std::complex<real> d = fd::z_of(Sc[vi]) - g;
                           o[0] = s * d.real();
                           o[1] = s * d.imag();
                           if (mean)
                             o[2] = s * (fd::H_of(Sc[vi]) - _H0[vi]);
                         }});
      }
      if (c.w_tension > 0) {
        const real s = c.w_tension * sv;
        units.push_back({1, {v}, {}, [&, s, vi](real *o) { o[0] = s * (fd::H_of(Sc[vi]) - _Ht); }});
      }
      if (c.w_growth > 0) {
        const real s = c.w_growth * sv;
        units.push_back({3, {v}, {v}, [&, s, vi](real *o) {
                           real ka, ta, kb, tb;
                           fd::along(Sc[vi], F[vi], _ga[vi], ka, ta);
                           fd::along(Sc[vi], F[vi], F[vi].N.cross(_ga[vi]), kb, tb);
                           o[0] = s * (ka - _gk[vi][0]);
                           o[1] = s * (kb - _gk[vi][1]);
                           o[2] = s * (ta - _gk[vi][2]);
                         }});
      }
      if (c.w_twist > 0) {
        const real s = c.w_twist * sv;
        units.push_back({1, {v}, {v}, [&, s, vi](real *o) {
                           real kt, tt;
                           fd::along(Sc[vi], F[vi], _tw[vi], kt, tt);
                           o[0] = s * (tt - _tw0[vi] - c.tau0);
                         }});
      }
      if (c.w_det > 0) {
        const real s = c.w_det * k * _hl;
        std::vector<int> verts;
        for (int f : _faces_of[vi])
          for (int kk = 0; kk < 3; ++kk)
            verts.push_back(_tp.tris[static_cast<size_t>(f)][static_cast<size_t>(kk)]);
        std::sort(verts.begin(), verts.end());
        verts.erase(std::unique(verts.begin(), verts.end()), verts.end());
        units.push_back({1, {}, verts, [&, s, vi](real *o) {
                           real om = 0;
                           for (int f : _faces_of[vi]) {
                             const auto &tr = _tp.tris[static_cast<size_t>(f)];
                             om += fd::solid_angle(F[static_cast<size_t>(tr[0])].N, F[static_cast<size_t>(tr[1])].N,
                                                   F[static_cast<size_t>(tr[2])].N);
                           }
                           o[0] = s * (om - _Om0[vi]);
                         }});
      }
    }
    const bool lines = c.w_dupin > 0 && c.dupin_reach > 0;
    if (lines)
      for (const dupin_line &dl : _lines) {
        const real wd = c.w_dupin;
        units.push_back({static_cast<int>(dl.C.rows()), dl.verts, {}, [&, wd](real *o) {
                           vecX Sv(3 * static_cast<Eigen::Index>(dl.verts.size()));
                           for (size_t a = 0; a < dl.verts.size(); ++a)
                             Sv.segment<3>(3 * static_cast<Eigen::Index>(a)) = Sc[static_cast<size_t>(dl.verts[a])];
                           const vecX r = wd * (dl.C * Sv);
                           for (Eigen::Index a = 0; a < r.size(); ++a)
                             o[a] = r[a];
                         }});
      }
    if (c.w_minimizer > 0 || (c.w_dupin > 0 && !lines)) {
      const real wm = c.w_minimizer, wd = lines ? real(0) : c.w_dupin;
      const int m = (wm > 0 ? 2 : 0) + (wd > 0 ? 2 : 0);
      for (size_t sp = 0; sp < _tp.spokes.size(); ++sp) {
        const int v = _tp.spokes[sp][0], w = _tp.spokes[sp][1];
        const size_t vi = static_cast<size_t>(v), wi = static_cast<size_t>(w);
        const real s = k * _A[vi];
        units.push_back({m, {v, w}, {v, w}, [&, vi, wi, s, wm, wd](real *o) {
                           real k1v, k1w, k2v, k2w, t;
                           fd::along(Sc[vi], F[vi], _e1[vi], k1v, t);
                           fd::along(Sc[wi], F[wi], _e1[vi], k1w, t);
                           fd::along(Sc[vi], F[vi], _e2[vi], k2v, t);
                           fd::along(Sc[wi], F[wi], _e2[vi], k2w, t);
                           const vec3 dh = (X[wi] - X[vi]).normalized();
                           const real g1 = dh.dot(_e1[vi]), g2 = dh.dot(_e2[vi]);
                           int at = 0;
                           if (wm > 0) {
                             o[at++] = wm * s * g2 * (k1w - k1v);
                             o[at++] = wm * s * g1 * (k2w - k2v);
                           }
                           if (wd > 0) {
                             o[at++] = wd * s * g1 * (k1w - k1v);
                             o[at++] = wd * s * g2 * (k2w - k2v);
                           }
                         }});
      }
    }
    if (c.w_target > 0 && _kt.size() == _tp.spokes.size() && _taut.size() == _tp.spokes.size()) {
      for (int v = 0; v < n; ++v) {
        const size_t vi = static_cast<size_t>(v);
        const int deg = _tp.spoke_begin[vi + 1] - _tp.spoke_begin[vi];
        const real s = c.w_target * k * _A[vi] / std::sqrt(std::max(real(deg), real(1)));
        for (int sp = _tp.spoke_begin[vi]; sp < _tp.spoke_begin[vi + 1]; ++sp) {
          const size_t si = static_cast<size_t>(sp);
          const int j = _tp.spokes[si][1];
          const size_t jj = static_cast<size_t>(j);
          // Chord form, as the net's Rodrigues rows: `(k - k_t) |dx_t|` over the
          // frame-start `|dx_t|`, no division by the live length (a spoke
          // turning toward the normal would blow up the reading).
          const real L0 = _spoke_l0[si] > real(0) ? _spoke_l0[si] : real(1);
          units.push_back({2, {}, {v, j}, [&, vi, jj, si, s, L0](real *o) {
                             const vec3 dx = X[jj] - X[vi], dn = F[jj].N - F[vi].N;
                             const vec3 nb = (F[vi].N + F[jj].N).normalized();
                             const vec3 dxt = dx - nb * dx.dot(nb);
                             const real L = dxt.norm();
                             o[0] = o[1] = 0;
                             if (!(L > real(0)))
                               return;
                             const vec3 t = dxt / L, b = nb.cross(t);
                             o[0] = s * (-dn.dot(t) - _kt[si] * L) / L0;
                             o[1] = s * (-dn.dot(b) - _taut[si] * L) / L0;
                           }});
        }
      }
    }
    if (c.w_gradtl > 0) {
      for (const auto &e : _tp.edges) {
        if (e[0] < 0)
          continue;
        const int i = e[0], j = e[1];
        const size_t ii = static_cast<size_t>(i), jj = static_cast<size_t>(j);
        const real s = c.w_gradtl * k * real(0.5) * (_A[ii] + _A[jj]);
        units.push_back({2, {i, j}, {i, j}, [&, ii, jj, s](real *o) {
                           const std::complex<real> d =
                               fd::z_of(Sc[ii]) -
                               fd::z_of(Sc[jj]) * std::polar(real(1), real(2) * fd::relative_turn(F[ii], F[jj]));
                           o[0] = s * d.real();
                           o[1] = s * d.imag();
                         }});
      }
    }
    index_t total = 0;
    for (const unit &u : units)
      total += u.m;
    r.resize(total);
    index_t row = 0;
    int mmax = 1;
    for (const unit &u : units)
      mmax = std::max(mmax, u.m);
    std::vector<real> bufv(static_cast<size_t>(mmax)), bpv(static_cast<size_t>(mmax)), bmv(static_cast<size_t>(mmax));
    real *buf = bufv.data(), *bp = bpv.data(), *bm = bmv.data();
    std::vector<std::vector<std::pair<index_t, real>>> acc(static_cast<size_t>(mmax));
    auto add = [&](int a, index_t col, real val) {
      if (val == real(0))
        return;
      for (auto &pr : acc[a])
        if (pr.first == col) {
          pr.second += val;
          return;
        }
      acc[a].push_back({col, val});
    };
    for (const unit &u : units) {
      u.f(buf);
      for (int a = 0; a < u.m; ++a)
        r[row + a] = buf[a];
      if (jt) {
        for (int a = 0; a < u.m; ++a)
          acc[a].clear();
        // Through the cached tensors.
        for (int v : u.sv) {
          const size_t vi = static_cast<size_t>(v);
          std::vector<int> ring;
          ring_of(v, ring);
          for (int cc = 0; cc < 3; ++cc) {
            const real keep = Sc[vi][cc], h = real(1e-7) * std::max(real(1), std::abs(keep));
            Sc[vi][cc] = keep + h;
            u.f(bp);
            Sc[vi][cc] = keep - h;
            u.f(bm);
            Sc[vi][cc] = keep;
            for (int a = 0; a < u.m; ++a) {
              const real dr = (bp[a] - bm[a]) / (real(2) * h);
              if (dr == real(0))
                continue;
              for (size_t b = 0; b < ring.size(); ++b)
                for (int comp = 0; comp < 7; ++comp)
                  add(a, col_of(ring[b], comp), dr * dS[vi](cc, 7 * static_cast<Eigen::Index>(b) + comp));
            }
          }
        }
        // Direct reads of positions and frames (tensors held).
        for (int v : u.dv) {
          const size_t vi = static_cast<size_t>(v);
          for (int comp = 0; comp < 7; ++comp) {
            const bool pos = comp < 3;
            real &slot = pos ? X[vi][comp] : Y[vi][comp - 3];
            const real keep = slot, h = real(1e-7) * std::max(real(1), std::abs(keep));
            slot = keep + h;
            if (!pos)
              F[vi] = fd::frame_of(Y[vi]);
            u.f(bp);
            slot = keep - h;
            if (!pos)
              F[vi] = fd::frame_of(Y[vi]);
            u.f(bm);
            slot = keep;
            if (!pos)
              F[vi] = fd::frame_of(Y[vi]);
            for (int a = 0; a < u.m; ++a)
              add(a, col_of(v, comp), (bp[a] - bm[a]) / (real(2) * h));
          }
        }
        for (int a = 0; a < u.m; ++a)
          for (const auto &pr : acc[a])
            jt->emplace_back(id0 + row + a, pr.first, pr.second);
      }
      row += u.m;
    }
  }

  const kusama::stencil_table *_t = nullptr;
  frame_rows_config _cfg;
  kusama::contact_net::topology _tp;
  std::vector<std::vector<int>> _faces_of;
  std::vector<real> _A;
  real _hl = 1;
  std::vector<frame3> _F0;
  std::vector<vec3> _S0, _e1, _e2, _ga, _gk, _tw;
  std::vector<real> _tw0, _H0, _Om0;
  std::vector<std::complex<real>> _z0, _zring0;
  real _Ht = 0;
  std::vector<real> _kt, _taut, _spoke_l0;
  std::vector<dupin_line> _lines;
};

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif
