#ifndef GAUDI_HEPWORTH_BLOCK_EDGE_SHAPE_FAIR_HPP
#define GAUDI_HEPWORTH_BLOCK_EDGE_SHAPE_FAIR_HPP

/// Projection-based fairing of the edge curvature store
/// (`shell_cr_block::shape_dof`): one constraint type for every curvature
/// operator. The rows `A` are fixed for the frame; `project` writes
/// `p = A g(q)` for a goal field `g` (`edge_shape_goals.hpp`), so the solve
/// pulls the store toward the goal's shape under that metric:
/// - `pin`: `c_e I` per edge;
/// - `dirichlet`: per face, `E_f^{1/2}` of the CR spin-2 element over its three
///   edge tensors (only the goal's variation; blind to a uniform change);
/// - `dirichlet_trace`: the same on the trace `k11 + k22` only;
/// - `dirichlet_traceless`: the same on the traceless part `z` only
///   (`spin2_traceless_face_element`): smooths the principal directions
///   weighted by anisotropy, mean curvature free. `dirichlet_trace` plus
///   `dirichlet_traceless` is `dirichlet` exactly;
/// - `h1`: `pin` stacked with `sqrt(lambda)` times `dirichlet`: a uniform
///   change costs its pin, a rough one more;
/// - `dupin`: per face, each principal curvature's derivative along its own
///   direction, `e_1 . grad k_1` and `e_2 . grad k_2` (the Dupin condition:
///   curvature lines are circles, as on a torus or any Dupin cyclide). The
///   face's principal frame is frozen from the frame-start store `h0`, so
///   `k_d = d^T S d` is linear in the store and the energy is a fixed
///   quadratic (variational); use it with `goal_zero`. With a `dupin_scale`
///   `s` the rows hold `s` times `h0`'s own gradients instead of zero: 1
///   keeps each line's curvature gradient (its skew) at its frame-start
///   value, above 1 grows it each frame, 0 is the Dupin condition.
/// - `dupin_line`: the Dupin condition on a long stencil. Per face and
///   principal direction (frozen from `h0`), the curvature line through the
///   centroid is traced `reach` mean edge lengths both ways
///   (`trace_on_surface`), the normal curvature along the transported
///   direction sampled every edge length (linear in the store: each
///   sample's face edges through the Crouzeix-Raviart basis), and the rows
///   are the residual of a least-squares polynomial fit of degree `degree`
///   along the line (0: constant curvature, the line a circle; 1: also lets
///   the curvature change linearly). A fixed quadratic: use it with
///   `goal_zero`.
/// Rows are `w` times the row's geometry, as every PD row: the coupling
/// factor `kusama::weingarten::coupling_factor` per edge (the face's mean on
/// the face rows), the same factor as the coupling row `weingarten_gn`, so a
/// pin at `w` and the coupling at `w` are an even least-squares blend of the
/// store.

#include "gaudi/hepworth/block/block_constraint.hpp"
#include "gaudi/hepworth/block/surface_walk.hpp"
#include "gaudi/kusama/stencil_table.hpp"
#include "gaudi/kusama/weingarten.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <algorithm>
#include <array>
#include <functional>
#include <memory>
#include <vector>

namespace gaudi {
namespace hepworth {
namespace block {

/// World frame of each edge: `t` from the lower vertex id to the higher, `n`
/// the mean of the incident face normals, `b = n x t`; and the other edges of
/// the edge's faces.
struct edge_frames {
  std::vector<vec3> t, b, n;
  std::vector<std::vector<int>> nbr;
};

inline std::shared_ptr<const edge_frames> build_edge_frames(const kusama::stencil_table &tb,
                                                            const std::vector<vec3> &x) {
  auto fr = std::make_shared<edge_frames>();
  const size_t nE = static_cast<size_t>(tb.n_edges);
  fr->t.assign(nE, vec3::UnitX());
  fr->b.assign(nE, vec3::UnitY());
  fr->n.assign(nE, vec3::Zero());
  fr->nbr.assign(nE, {});
  for (int f = 0; f < tb.n_faces; ++f) {
    const auto &id = tb.ids[static_cast<size_t>(f)];
    const vec3 N = (x[static_cast<size_t>(id[1])] - x[static_cast<size_t>(id[0])])
                       .cross(x[static_cast<size_t>(id[2])] - x[static_cast<size_t>(id[0])]);
    const vec3 Nn = N.norm() > real(0) ? vec3(N.normalized()) : vec3::Zero();
    for (int i = 0; i < 3; ++i) {
      const size_t e = static_cast<size_t>(tb.face_edges[static_cast<size_t>(f)][static_cast<size_t>(i)]);
      const int a = id[static_cast<size_t>(i)], b = id[static_cast<size_t>((i + 1) % 3)];
      fr->t[e] = x[static_cast<size_t>(std::max(a, b))] - x[static_cast<size_t>(std::min(a, b))];
      fr->n[e] += Nn;
      for (int j = 1; j < 3; ++j)
        fr->nbr[e].push_back(tb.face_edges[static_cast<size_t>(f)][static_cast<size_t>((i + j) % 3)]);
    }
  }
  for (size_t e = 0; e < nE; ++e) {
    vec3 &n = fr->n[e];
    vec3 &t = fr->t[e];
    if (n.norm() > real(0))
      n.normalize();
    t -= n * t.dot(n);
    if (t.norm() > real(0))
      t.normalize();
    fr->b[e] = n.cross(t);
  }
  return fr;
}

/// `(k11, k12, k22)` in `(t, b)` as a world tensor.
inline mat3 edge_world_tensor(const vec3 &h, const vec3 &t, const vec3 &b) {
  return h[0] * t * t.transpose() + h[1] * (t * b.transpose() + b * t.transpose()) +
         h[2] * b * b.transpose();
}

enum class edge_fair_rows { pin, dirichlet, dirichlet_trace, h1, dupin, dirichlet_traceless, dupin_line };

/// Rows `dupin_line`: how far each curvature line is traced (mean edge
/// lengths each way) and the degree of the fit whose residual the rows are.
struct dupin_line_spec {
  real reach = 3;
  int degree = 0;
};

/// Goal tensor per edge from the live state.
using edge_goal_fn = std::function<void(const vecX &q, std::vector<vec3> &g)>;

class edge_shape_fair : public block_constraint {
public:
  DEFINE_CREATE_FUNC(edge_shape_fair)

  /// `blocks[0]` is the edge store; the goal is identity when `goal` is
  /// empty. `lambda` weighs the Dirichlet part of rows `h1`. Rows `dupin`
  /// read their principal frames from `h0` (the frame-start store).
  /// `tensor_metric` weighs the pin's `k12` by `sqrt(2)`, the tensor norm
  /// `k11^2 + 2 k12^2 + k22^2`, in which `goal_anisotropy` is an exact
  /// projection. `dupin_scale` sets rows `dupin` to hold that multiple of
  /// `h0`'s gradients (0: the Dupin zero).
  edge_shape_fair(const kusama::stencil_table *table, const std::vector<vec3> &x,
                  edge_fair_rows rows, real w, std::vector<sim_block::ptr> blocks,
                  edge_goal_fn goal = {}, real lambda = real(1),
                  const std::vector<vec3> *h0 = nullptr, bool tensor_metric = false,
                  real dupin_scale = real(0), dupin_line_spec line = {})
      : block_constraint({}, w, std::move(blocks)), _nE(table->n_edges),
        _goal(std::move(goal)), _tensor_metric(tensor_metric), _line(line) {
    build(*table, x, rows, lambda, h0);
    if (rows == edge_fair_rows::dupin && dupin_scale != real(0) && std::isfinite(dupin_scale) &&
        h0 && static_cast<int>(h0->size()) == _nE) {
      _target.assign(static_cast<size_t>(_n_rows), real(0));
      for (const trip &t : _rows)
        _target[static_cast<size_t>(t.row())] +=
            dupin_scale * t.value() * (*h0)[static_cast<size_t>(t.col() / 3)][t.col() % 3];
    }
    _rest_initialized = true;
  }

  std::string name() override { return "edge_shape_fair"; }

  void fill_A(index_t &id0, std::vector<trip> &triplets) override {
    _id0 = id0;
    for (const trip &t : _rows)
      triplets.emplace_back(id0 + t.row(),
                            _blocks[0]->get_offset_idx(t.col() / 3) + t.col() % 3,
                            _w * t.value());
    id0 += _n_rows;
  }

  void project(const vecX &q, vecX &p) override {
    std::vector<vec3> g(static_cast<size_t>(_nE));
    if (_goal)
      _goal(q, g);
    else
      for (int e = 0; e < _nE; ++e)
        g[static_cast<size_t>(e)] = _blocks[0]->get_vec3(e, q);
    for (vec3 &v : g)
      if (!v.allFinite())
        v = vec3::Zero();
    for (const trip &t : _rows)
      p[_id0 + t.row()] += _w * t.value() * g[static_cast<size_t>(t.col() / 3)][t.col() % 3];
    for (size_t r = 0; r < _target.size(); ++r)
      p[_id0 + static_cast<index_t>(r)] += _w * _target[r];
  }

  index_t row_count() const { return _n_rows; }

private:
  /// Rows in local columns `3 e + c`.
  void build(const kusama::stencil_table &t, const std::vector<vec3> &x, edge_fair_rows rows,
             real lambda, const std::vector<vec3> *h0) {
    namespace wg = kusama::weingarten;
    _c = wg::coupling_factor(t, x);
    _n_rows = 0;
    if (rows == edge_fair_rows::dupin) {
      if (h0 && static_cast<int>(h0->size()) == t.n_edges)
        dupin_rows(t, x, *h0);
      return;
    }
    if (rows == edge_fair_rows::dupin_line) {
      if (h0 && static_cast<int>(h0->size()) == t.n_edges)
        dupin_line_rows(t, x, *h0);
      return;
    }
    if (rows == edge_fair_rows::pin || rows == edge_fair_rows::h1)
      pin_rows();
    if (rows == edge_fair_rows::pin)
      return;
    const real sl = rows == edge_fair_rows::h1 ? std::sqrt(std::max(lambda, real(0))) : real(1);
    for (int f = 0; f < t.n_faces; ++f) {
      const auto &id = t.ids[static_cast<size_t>(f)];
      const vec3 &X0 = x[static_cast<size_t>(id[0])];
      const vec3 &X1 = x[static_cast<size_t>(id[1])];
      const vec3 &X2 = x[static_cast<size_t>(id[2])];
      const auto &fe = t.face_edges[static_cast<size_t>(f)];
      const real s = sl * (_c[fe[0]] + _c[fe[1]] + _c[fe[2]]) / real(3);
      if (rows != edge_fair_rows::dirichlet_trace) {
        const Eigen::Matrix<real, 9, 9> E = rows == edge_fair_rows::dirichlet_traceless
                                                ? wg::spin2_traceless_face_element(X0, X1, X2)
                                                : wg::spin2_face_element(X0, X1, X2);
        if (!E.allFinite())
          continue;
        const Eigen::Matrix<real, 9, 9> A = psd_sqrt(E);
        for (int r = 0; r < 9; ++r)
          for (int a = 0; a < 3; ++a)
            for (int c = 0; c < 3; ++c) {
              const real v = s * A(r, 3 * a + c);
              if (v != real(0))
                _rows.emplace_back(_n_rows + r,
                                   3 * t.face_edges[static_cast<size_t>(f)][static_cast<size_t>(a)] + c, v);
            }
        _n_rows += 9;
      } else {
        // `1/2 tau^T K tau` matches the trace part of the spin-2 inner product.
        const Eigen::Matrix<real, 3, 3> K = wg::cr_face_stiffness(X0, X1, X2);
        if (!K.allFinite())
          continue;
        const Eigen::Matrix<real, 3, 3> A = psd_sqrt(Eigen::Matrix<real, 3, 3>(real(0.5) * K));
        for (int r = 0; r < 3; ++r)
          for (int a = 0; a < 3; ++a) {
            const real v = s * A(r, a);
            if (v == real(0))
              continue;
            const int e = t.face_edges[static_cast<size_t>(f)][static_cast<size_t>(a)];
            _rows.emplace_back(_n_rows + r, 3 * e + 0, v);
            _rows.emplace_back(_n_rows + r, 3 * e + 2, v);
          }
        _n_rows += 3;
      }
    }
  }

  void dupin_rows(const kusama::stencil_table &t, const std::vector<vec3> &x,
                  const std::vector<vec3> &h0) {
    const auto fr = build_edge_frames(t, x);
    for (int f = 0; f < t.n_faces; ++f) {
      const auto &id = t.ids[static_cast<size_t>(f)];
      const auto &fe = t.face_edges[static_cast<size_t>(f)];
      const std::array<vec3, 3> X{x[static_cast<size_t>(id[0])], x[static_cast<size_t>(id[1])],
                                  x[static_cast<size_t>(id[2])]};
      const vec3 Nf = (X[1] - X[0]).cross(X[2] - X[0]);
      const real A = real(0.5) * Nf.norm();
      if (!(A > real(0)) || !std::isfinite(A))
        continue;
      const vec3 n = Nf.normalized();
      // The face's principal frame: the mean of its edges' tensors in the face plane.
      mat3 T = mat3::Zero();
      for (int i = 0; i < 3; ++i) {
        const size_t e = static_cast<size_t>(fe[static_cast<size_t>(i)]);
        T += edge_world_tensor(h0[e], fr->t[e], fr->b[e]) / real(3);
      }
      const vec3 u = (X[1] - X[0]).normalized(), v = n.cross(u);
      Eigen::Matrix<real, 2, 2> S2;
      S2 << u.dot(T * u), u.dot(T * v), v.dot(T * u), v.dot(T * v);
      Eigen::SelfAdjointEigenSolver<Eigen::Matrix<real, 2, 2>> es(real(0.5) * (S2 + S2.transpose()));
      const real k1 = es.eigenvalues()[0], k2 = es.eigenvalues()[1];
      // Principal directions are undefined at an umbilic.
      if (!(std::abs(k1 - k2) > real(1e-3) * (std::abs(k1) + std::abs(k2))))
        continue;
      const std::array<vec3, 2> dir{es.eigenvectors()(0, 0) * u + es.eigenvectors()(1, 0) * v,
                                    es.eigenvectors()(0, 1) * u + es.eigenvectors()(1, 1) * v};
      const real s = std::sqrt(A) * (_c[fe[0]] + _c[fe[1]] + _c[fe[2]]) / real(3);
      for (const vec3 &d : dir) {
        for (int i = 0; i < 3; ++i) {
          // Crouzeix-Raviart basis of local edge i (vertex i to i + 1).
          const vec3 grad = -n.cross(X[static_cast<size_t>((i + 1) % 3)] - X[static_cast<size_t>(i)]) / A;
          const int e = fe[static_cast<size_t>(i)];
          const real a = d.dot(fr->t[static_cast<size_t>(e)]), b = d.dot(fr->b[static_cast<size_t>(e)]);
          const real l2 = a * a + b * b;
          if (!(l2 > real(1e-16)))
            continue;
          const real g = s * d.dot(grad) / l2;
          _rows.emplace_back(_n_rows, 3 * e + 0, g * a * a);
          _rows.emplace_back(_n_rows, 3 * e + 1, g * real(2) * a * b);
          _rows.emplace_back(_n_rows, 3 * e + 2, g * b * b);
        }
        ++_n_rows;
      }
    }
  }

  /// The face's principal directions from `h0` (none at an umbilic).
  static bool face_principals(const kusama::stencil_table &t, const std::vector<vec3> &x,
                              const std::vector<vec3> &h0, const edge_frames &fr, int f,
                              std::array<vec3, 2> &dir, vec3 &n, real &A) {
    const auto &id = t.ids[static_cast<size_t>(f)];
    const auto &fe = t.face_edges[static_cast<size_t>(f)];
    const std::array<vec3, 3> X{x[static_cast<size_t>(id[0])], x[static_cast<size_t>(id[1])],
                                x[static_cast<size_t>(id[2])]};
    const vec3 Nf = (X[1] - X[0]).cross(X[2] - X[0]);
    A = real(0.5) * Nf.norm();
    if (!(A > real(0)) || !std::isfinite(A))
      return false;
    n = Nf.normalized();
    mat3 T = mat3::Zero();
    for (int i = 0; i < 3; ++i) {
      const size_t e = static_cast<size_t>(fe[static_cast<size_t>(i)]);
      T += edge_world_tensor(h0[e], fr.t[e], fr.b[e]) / real(3);
    }
    const vec3 u = (X[1] - X[0]).normalized(), v = n.cross(u);
    Eigen::Matrix<real, 2, 2> S2;
    S2 << u.dot(T * u), u.dot(T * v), v.dot(T * u), v.dot(T * v);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<real, 2, 2>> es(real(0.5) * (S2 + S2.transpose()));
    const real k1 = es.eigenvalues()[0], k2 = es.eigenvalues()[1];
    if (!(std::abs(k1 - k2) > real(1e-3) * (std::abs(k1) + std::abs(k2))))
      return false;
    dir = {es.eigenvectors()(0, 0) * u + es.eigenvectors()(1, 0) * v,
           es.eigenvectors()(0, 1) * u + es.eigenvectors()(1, 1) * v};
    return true;
  }

  void dupin_line_rows(const kusama::stencil_table &t, const std::vector<vec3> &x,
                       const std::vector<vec3> &h0) {
    const auto fr = build_edge_frames(t, x);
    real hm = 0;
    int ne = 0;
    for (int f = 0; f < t.n_faces; ++f)
      for (int i = 0; i < 3; ++i) {
        const auto &id = t.ids[static_cast<size_t>(f)];
        hm += (x[static_cast<size_t>(id[static_cast<size_t>((i + 1) % 3)])] - x[static_cast<size_t>(id[static_cast<size_t>(i)])]).norm();
        ++ne;
      }
    hm /= std::max(ne, 1);
    const real reach = std::max(_line.reach, real(0)) * hm;
    const int deg = std::max(_line.degree, 0);
    for (int f = 0; f < t.n_faces; ++f) {
      std::array<vec3, 2> dir;
      vec3 n;
      real A = 0;
      if (!face_principals(t, x, h0, *fr, f, dir, n, A))
        continue;
      const auto &id = t.ids[static_cast<size_t>(f)];
      const auto &fe = t.face_edges[static_cast<size_t>(f)];
      const vec3 c = (x[static_cast<size_t>(id[0])] + x[static_cast<size_t>(id[1])] + x[static_cast<size_t>(id[2])]) / real(3);
      const real s = std::sqrt(A) * (_c[fe[0]] + _c[fe[1]] + _c[fe[2]]) / real(3);
      // A curvature line: one edge length at a time, the direction snapped
      // after each step to the closer principal direction of the face reached
      // (a straightest walk alone drifts off lines that are not geodesics).
      auto trace_line = [&](vec3 d0) {
        std::vector<walk_sample> out;
        int g = f;
        vec3 p = c, d = d0;
        real s0 = 0;
        while (s0 < reach - real(1e-9)) {
          std::vector<walk_sample> seg = trace_on_surface(t, x, g, p, d, std::min(hm, reach - s0), hm);
          if (seg.empty())
            break;
          if (out.empty())
            out.push_back(seg.front());
          if (seg.size() < 2)
            break;
          walk_sample w = seg.back();
          w.s += s0;
          std::array<vec3, 2> pd;
          vec3 ng;
          real Ag = 0;
          if (face_principals(t, x, h0, *fr, w.face, pd, ng, Ag)) {
            const vec3 &best = std::abs(pd[0].dot(w.dir)) >= std::abs(pd[1].dot(w.dir)) ? pd[0] : pd[1];
            w.dir = best.dot(w.dir) >= 0 ? best : vec3(-best);
          }
          out.push_back(w);
          g = w.face;
          p = w.p;
          d = w.dir;
          s0 = w.s;
        }
        return out;
      };
      for (const vec3 &d : dir) {
        // The line both ways, as samples ordered by signed arc length.
        std::vector<walk_sample> fwd = trace_line(d);
        std::vector<walk_sample> bwd = trace_line(-d);
        std::vector<walk_sample> line;
        for (size_t k = bwd.size(); k-- > 1;) {
          walk_sample w = bwd[k];
          w.s = -w.s;
          line.push_back(w);
        }
        for (const walk_sample &w : fwd)
          line.push_back(w);
        const int m = static_cast<int>(line.size());
        if (m < deg + 2)
          continue;
        // Each sample's normal curvature along its direction, as a row over store columns.
        std::vector<std::vector<std::pair<int, real>>> C(static_cast<size_t>(m));
        for (int j = 0; j < m; ++j) {
          const walk_sample &w = line[static_cast<size_t>(j)];
          const auto &ge = t.face_edges[static_cast<size_t>(w.face)];
          for (int i = 0; i < 3; ++i) {
            // Crouzeix-Raviart basis of local edge i (vertices i, i + 1): 1 - 2 lambda_{i + 2}.
            const real phi = real(1) - real(2) * w.bary[(i + 2) % 3];
            const int e = ge[static_cast<size_t>(i)];
            const real a = w.dir.dot(fr->t[static_cast<size_t>(e)]), b = w.dir.dot(fr->b[static_cast<size_t>(e)]);
            const real l2 = a * a + b * b;
            if (!(l2 > real(1e-16)))
              continue;
            C[static_cast<size_t>(j)].push_back({3 * e + 0, phi * a * a / l2});
            C[static_cast<size_t>(j)].push_back({3 * e + 1, phi * real(2) * a * b / l2});
            C[static_cast<size_t>(j)].push_back({3 * e + 2, phi * b * b / l2});
          }
        }
        // Residual of the degree-`deg` fit: the orthonormal complement of the
        // polynomial columns (arc length in mean edge lengths).
        Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> V(m, deg + 1);
        for (int j = 0; j < m; ++j)
          for (int p = 0; p <= deg; ++p)
            V(j, p) = std::pow(line[static_cast<size_t>(j)].s / hm, p);
        Eigen::HouseholderQR<Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic>> qr(V);
        const Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> Q =
            qr.householderQ() * Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic>::Identity(m, m);
        // Curvature deviations scaled like the face rows' gradients: per
        // edge length, normalized by the number of samples.
        const real sc = s / (hm * std::sqrt(real(m)));
        for (int r = deg + 1; r < m; ++r) {
          for (int j = 0; j < m; ++j) {
            const real qj = Q(j, r);
            if (qj == real(0))
              continue;
            for (const auto &[col, v] : C[static_cast<size_t>(j)])
              _rows.emplace_back(_n_rows, col, sc * qj * v);
          }
          ++_n_rows;
        }
      }
    }
  }

  void pin_rows() {
    for (int e = 0; e < _nE; ++e) {
      const real s = _c[static_cast<size_t>(e)];
      for (int c = 0; c < 3; ++c)
        _rows.emplace_back(_n_rows + c, 3 * e + c,
                           c == 1 && _tensor_metric ? s * std::sqrt(real(2)) : s);
      _n_rows += 3;
    }
  }

  template <int N> static Eigen::Matrix<real, N, N> psd_sqrt(const Eigen::Matrix<real, N, N> &E) {
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<real, N, N>> es(real(0.5) * (E + E.transpose()));
    // Round-off in the null space would leave sqrt(1e-16) rows.
    const real floor = real(1e-12) * es.eigenvalues().cwiseAbs().maxCoeff();
    Eigen::Matrix<real, N, 1> ev = es.eigenvalues();
    for (int i = 0; i < N; ++i)
      ev[i] = ev[i] > floor ? std::sqrt(ev[i]) : real(0);
    return es.eigenvectors() * ev.asDiagonal() * es.eigenvectors().transpose();
  }

  int _nE = 0;
  index_t _n_rows = 0;
  std::vector<real> _c;
  std::vector<trip> _rows;
  /// Constant row values held for the step (rows `dupin` with a scale).
  std::vector<real> _target;
  edge_goal_fn _goal;
  bool _tensor_metric = false;
  dupin_line_spec _line{};
};

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif
