#ifndef GAUDI_ENGMAN_ENERGIES_HPP
#define GAUDI_ENGMAN_ENERGIES_HPP

/// Soft energies for the Engman solver (`engman/solver.hpp`), in the two row
/// forms it takes: projection rows (`fill_A`, `project`; constant `A^T A`)
/// and Gauss-Newton rows (`fill_J`, `add_gradient`, `residual_sq`).
///
/// - `as_energy`: any `engman::constraint` as a weighted Gauss-Newton energy,
///   `1/2 w^2 |C(q)|^2`, the soft side of a hard/soft switch;
/// - `store_bending`: bending on the edge store, `plate` `w int |W|^2` or
///   `helfrich` `w int (2H - c0)^2`, rows `sqrt(w A_e)`;
/// - `torsion_target`: a target on the geodesic torsion along an axis only
///   (the Helfrich-Prost chiral term), the curvatures free; `read_torsion`
///   measures it;
/// - `deviator_hold`: holds each edge's anisotropy `|z_e|` at a fixed value
///   (its phase free), an exact projection, so a twist cannot be met by
///   erasing the anisotropy;
/// - `twisted_spin2`: the traceless spin-2 CR Dirichlet energy with its
///   connection shifted by `omega` along an axis field `a`,
///   `w int |grad_a z - i omega z|^2` (plus the across-axis gradient), zero on
///   a principal field that turns at rate `omega` along `a` against parallel
///   transport;
/// - `twist_rate`: the same target as a scale-free Gauss-Newton energy on
///   the measured rate, `w sum A_f |z0_f|^2 (omega_f(S) - omega)^2`, which
///   anisotropy loss cannot satisfy (the default form);
/// - `twist_readout`: the measured twist rate `Im(conj z d_a z) / |z|^2` per
///   face and its `A |z|`-weighted mean.
///
/// Store tensors are `(k11, k12, k22)` in each edge's frame `(t_e, N x t_e)`;
/// `z = ((k11 - k22) / 2, k12)` is the traceless part as a complex number.

#include "gaudi/engman/constraint.hpp"
#include "gaudi/hepworth/block/edge_shape_fair.hpp"
#include "gaudi/hepworth/block/sim_block.hpp"
#include "gaudi/hepworth/projection_constraint.hpp"
#include "gaudi/kusama/stencil_table.hpp"
#include "gaudi/kusama/weingarten.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace gaudi {
namespace engman {

// ------------------------------------------------------------------ helpers

namespace detail {

inline real face_area(const vec3 &a, const vec3 &b, const vec3 &c) {
  return real(0.5) * (b - a).cross(c - a).norm();
}

/// Dual area per edge: the incident face areas over 3.
inline std::vector<real> edge_dual_area(const kusama::stencil_table &t, const std::vector<vec3> &x) {
  std::vector<real> A(static_cast<size_t>(t.n_edges), real(0));
  for (int f = 0; f < t.n_faces; ++f) {
    const auto &id = t.ids[static_cast<size_t>(f)];
    const real a = face_area(x[static_cast<size_t>(id[0])], x[static_cast<size_t>(id[1])],
                             x[static_cast<size_t>(id[2])]);
    for (int i = 0; i < 3; ++i)
      A[static_cast<size_t>(t.face_edges[static_cast<size_t>(f)][static_cast<size_t>(i)])] += a / real(3);
  }
  return A;
}

/// Turns `z` of `(k11, k12, k22)` by `beta` (`z -> e^{i beta} z`), the trace kept.
inline Eigen::Matrix<real, 3, 3> rotate_z(real beta) {
  const real c = std::cos(beta), s = std::sin(beta);
  Eigen::Matrix<real, 3, 3> R;
  R << (1 + c) / 2, -s, (1 - c) / 2, //
      s / 2, c, -s / 2,              //
      (1 - c) / 2, s, (1 + c) / 2;
  return R;
}

template <int N>
Eigen::Matrix<real, N, N> psd_sqrt(const Eigen::Matrix<real, N, N> &E) {
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix<real, N, N>> es(real(0.5) * (E + E.transpose()));
  const real floor = real(1e-12) * std::max(es.eigenvalues().cwiseAbs().maxCoeff(), real(1e-300));
  Eigen::Matrix<real, N, 1> ev = es.eigenvalues();
  for (int i = 0; i < N; ++i)
    ev[i] = ev[i] > floor ? std::sqrt(ev[i]) : real(0);
  return es.eigenvectors() * ev.asDiagonal() * es.eigenvectors().transpose();
}

/// Face frame `(T, B, N)` as `kusama::weingarten::spin2_face_element_q` builds it.
inline void face_frame(const vec3 &X0, const vec3 &X1, const vec3 &X2, vec3 &T, vec3 &B, vec3 &N) {
  N = (X1 - X0).cross(X2 - X0).normalized();
  T = (X1 - X0).normalized();
  B = N.cross(T);
}

} // namespace detail

/// Axis field: a tangent direction per face from the face's centroid and
/// normal (projected into the face plane and normalized by the energy).
using axis_field = std::function<vec3(const vec3 &centroid, const vec3 &normal)>;

// ---------------------------------------------------------------- as_energy

/// An `engman::constraint` as the Gauss-Newton energy `1/2 w^2 |C(q)|^2`:
/// the soft side of a hard/soft switch. At large `w` it approaches the hard
/// constraint as `1/w^2`.
class as_energy : public hepworth::projection_constraint {
public:
  using ptr = std::shared_ptr<as_energy>;
  static ptr create(constraint::ptr c, real w) { return std::make_shared<as_energy>(std::move(c), w); }
  as_energy(constraint::ptr c, real w) : hepworth::projection_constraint({}, w), _c(std::move(c)) {}

  std::string name() override { return "soft_" + _c->name(); }
  bool gauss_newton() const override { return true; }

  void fill_J(const vecX &q, index_t &id0, std::vector<trip> &triplets) override {
    vecX C;
    std::vector<trip> J;
    _c->evaluate(q, C, &J);
    for (const trip &t : J)
      triplets.emplace_back(id0 + t.row(), t.col(), _w * t.value());
    id0 += _c->rows();
  }
  void add_gradient(const vecX &q, vecX &g) override {
    vecX C;
    std::vector<trip> J;
    _c->evaluate(q, C, &J);
    for (const trip &t : J)
      g[t.col()] += _w * _w * t.value() * C[t.row()];
  }
  real residual_sq(const vecX &q) override {
    vecX C;
    _c->evaluate(q, C, nullptr);
    return _w * _w * C.squaredNorm();
  }

private:
  constraint::ptr _c;
};

// ------------------------------------------------------- store edge energies

/// Per-edge projection rows on the store: `rows_e` (k x 3) and target
/// `p_e(h)`, both fixed by the subclass at construction.
class store_edge_energy : public hepworth::projection_constraint {
public:
  store_edge_energy(int n_edges, int k, real w, hepworth::sim_block::ptr store)
      : hepworth::projection_constraint({}, w), _nE(n_edges), _k(k), _store(std::move(store)),
        _rows(static_cast<size_t>(n_edges), Eigen::MatrixXd::Zero(k, 3)) {}

  void fill_A(index_t &id0, std::vector<trip> &triplets) override {
    _id0 = id0;
    for (int e = 0; e < _nE; ++e) {
      const Eigen::MatrixXd &R = _rows[static_cast<size_t>(e)];
      for (int r = 0; r < _k; ++r)
        for (int c = 0; c < 3; ++c)
          if (R(r, c) != real(0))
            triplets.emplace_back(id0 + _k * e + r, _store->get_offset_idx(e) + c, R(r, c));
    }
    id0 += static_cast<index_t>(_k) * _nE;
  }

  void project(const vecX &q, vecX &p) override {
    for (int e = 0; e < _nE; ++e) {
      const vec3 h = _store->get_vec3(e, q);
      const Eigen::VectorXd pe = target(e, h);
      if (pe.allFinite())
        p.segment(_id0 + static_cast<index_t>(_k) * e, _k) += pe;
    }
  }

protected:
  /// `p_e` for the live tensor `h` (rows times the goal).
  virtual Eigen::VectorXd target(int e, const vec3 &h) const = 0;

  int _nE = 0, _k = 0;
  hepworth::sim_block::ptr _store;
  std::vector<Eigen::MatrixXd> _rows;
};

enum class bending_mode { off, plate, helfrich };

/// Bending as a store energy. `plate`: `w sum A_e |S_e|^2_F`
/// (`k11^2 + 2 k12^2 + k22^2`, the thin-shell `int |W|^2`). `helfrich`:
/// `w sum A_e (tau_e - c0)^2`, `tau = 2H`. Quadratic in the store: under the
/// hard coupling it reaches positions only through the coupling.
class store_bending : public store_edge_energy {
public:
  using ptr = std::shared_ptr<store_bending>;
  static ptr create(const kusama::stencil_table &t, const std::vector<vec3> &x, bending_mode mode, real w,
                    real c0, hepworth::sim_block::ptr store) {
    return std::make_shared<store_bending>(t, x, mode, w, c0, std::move(store));
  }
  store_bending(const kusama::stencil_table &t, const std::vector<vec3> &x, bending_mode mode, real w,
                real c0, hepworth::sim_block::ptr store)
      : store_edge_energy(t.n_edges, mode == bending_mode::helfrich ? 1 : 3, w, std::move(store)),
        _mode(mode), _c0(c0) {
    const std::vector<real> A = detail::edge_dual_area(t, x);
    _s.resize(A.size());
    for (size_t e = 0; e < A.size(); ++e) {
      _s[e] = std::sqrt(std::max(w * A[e], real(0)));
      Eigen::MatrixXd &R = _rows[e];
      if (_mode == bending_mode::helfrich) {
        R << _s[e], 0, _s[e];
      } else {
        R.setZero();
        R(0, 0) = _s[e];
        R(1, 1) = std::sqrt(real(2)) * _s[e];
        R(2, 2) = _s[e];
      }
    }
  }
  std::string name() override { return _mode == bending_mode::helfrich ? "helfrich" : "plate_bending"; }

protected:
  Eigen::VectorXd target(int e, const vec3 &) const override {
    if (_mode == bending_mode::helfrich)
      return Eigen::VectorXd::Constant(1, _s[static_cast<size_t>(e)] * _c0);
    return Eigen::VectorXd::Zero(3);
  }

private:
  bending_mode _mode;
  real _c0 = 0;
  std::vector<real> _s;
};

/// Holds `|z_e|` at `mag_e` with `z`'s phase and the trace free:
/// `w sum A_e | z_e - mag_e z_e / |z_e| |^2`. The set is fixed for the step,
/// so the goal is the exact closest point and the row is variational.
class deviator_hold : public store_edge_energy {
public:
  using ptr = std::shared_ptr<deviator_hold>;
  static ptr create(const kusama::stencil_table &t, const std::vector<vec3> &x, std::vector<real> mag, real w,
                    hepworth::sim_block::ptr store) {
    return std::make_shared<deviator_hold>(t, x, std::move(mag), w, std::move(store));
  }
  deviator_hold(const kusama::stencil_table &t, const std::vector<vec3> &x, std::vector<real> mag, real w,
                hepworth::sim_block::ptr store)
      : store_edge_energy(t.n_edges, 2, w, std::move(store)), _mag(std::move(mag)) {
    const std::vector<real> A = detail::edge_dual_area(t, x);
    _s.resize(A.size());
    for (size_t e = 0; e < A.size(); ++e) {
      _s[e] = std::sqrt(std::max(w * A[e], real(0)));
      _rows[e] << _s[e] / 2, 0, -_s[e] / 2, //
          0, _s[e], 0;
    }
  }
  std::string name() override { return "deviator_hold"; }

  /// `|z|` of every edge of `h`.
  static std::vector<real> magnitudes(const std::vector<vec3> &h) {
    std::vector<real> m(h.size());
    for (size_t e = 0; e < h.size(); ++e)
      m[e] = std::hypot((h[e][0] - h[e][2]) / 2, h[e][1]);
    return m;
  }

protected:
  Eigen::VectorXd target(int e, const vec3 &h) const override {
    const real za = (h[0] - h[2]) / 2, zb = h[1];
    const real n = std::hypot(za, zb);
    const real m = _mag[static_cast<size_t>(e)];
    Eigen::VectorXd p(2);
    if (n > real(1e-12) * (m + real(1e-300)))
      p << _s[static_cast<size_t>(e)] * m * za / n, _s[static_cast<size_t>(e)] * m * zb / n;
    else
      p << _s[static_cast<size_t>(e)] * za, _s[static_cast<size_t>(e)] * zb; // no direction: leave it
    return p;
  }

private:
  std::vector<real> _mag, _s;
};

namespace detail {

/// Midpoint of edge `e` from either of its faces.
inline vec3 edge_midpoint(const kusama::stencil_table &t, const std::vector<vec3> &x, int e) {
  const size_t es = static_cast<size_t>(e);
  const bool first = t.edge_faces[es][0] >= 0;
  const int f = first ? t.edge_faces[es][0] : t.edge_faces[es][1];
  const int i = first ? t.edge_local[es][0] : t.edge_local[es][1];
  const auto &id = t.ids[static_cast<size_t>(f)];
  return real(0.5) * (x[static_cast<size_t>(id[static_cast<size_t>(i)])] +
                      x[static_cast<size_t>(id[static_cast<size_t>((i + 1) % 3)])]);
}

/// The row of `tau_g = W(a, N x a)` on `(k11, k12, k22)` in edge `e`'s frame:
/// `(-sin 2phi / 2, cos 2phi, sin 2phi / 2)`, `phi` the angle of `a` from
/// `t_e`. Zero when the axis has no tangent part.
inline Eigen::RowVector3d torsion_row(const hepworth::block::edge_frames &fr, size_t e, const vec3 &axis) {
  vec3 a = axis - fr.n[e] * axis.dot(fr.n[e]);
  if (!(a.norm() > real(0)))
    return Eigen::RowVector3d::Zero();
  a.normalize();
  const real c = a.dot(fr.t[e]), sn = a.dot(fr.b[e]);
  const real s2 = real(2) * sn * c, c2 = c * c - sn * sn;
  return Eigen::RowVector3d(-s2 / 2, c2, s2 / 2);
}

} // namespace detail

/// A torsion target on the store, only on the geodesic torsion along an axis
/// field: `w sum A_e (tau_g(S_e) - tau)^2`, `tau_g = W(a, N x a) =
/// ((k0 - k1) / 2) sin 2 delta`. One row per edge, linear in the store: the
/// trace, the anisotropy magnitude and the principal angle are all free, so
/// the solve can meet the twist by growing the anisotropy, turning the
/// principal frame, or buckling, whichever the coupling and the other terms
/// make cheapest. Its linear part is the Helfrich-Prost chiral term
/// `-2 w tau int tau_g`; the sign of `tau` is the handedness.
class torsion_target : public store_edge_energy {
public:
  using ptr = std::shared_ptr<torsion_target>;
  static ptr create(const kusama::stencil_table &t, const std::vector<vec3> &x, const axis_field &axis,
                    real tau, real w, hepworth::sim_block::ptr store) {
    return std::make_shared<torsion_target>(t, x, axis, tau, w, std::move(store));
  }
  torsion_target(const kusama::stencil_table &t, const std::vector<vec3> &x, const axis_field &axis, real tau,
                 real w, hepworth::sim_block::ptr store)
      : store_edge_energy(t.n_edges, 1, w, std::move(store)) {
    const std::vector<real> A = detail::edge_dual_area(t, x);
    const auto fr = hepworth::block::build_edge_frames(t, x);
    _p.assign(A.size(), real(0));
    for (size_t e = 0; e < A.size(); ++e) {
      const real sc = std::sqrt(std::max(w * A[e], real(0)));
      const Eigen::RowVector3d r =
          detail::torsion_row(*fr, e, axis(detail::edge_midpoint(t, x, static_cast<int>(e)), fr->n[e]));
      _rows[e] = sc * r;
      _p[e] = r.squaredNorm() > real(0) ? sc * tau : real(0);
    }
  }
  std::string name() override { return "torsion_target"; }

protected:
  Eigen::VectorXd target(int e, const vec3 &) const override {
    return Eigen::VectorXd::Constant(1, _p[static_cast<size_t>(e)]);
  }

private:
  std::vector<real> _p;
};

/// Geodesic torsion of the store along an axis field, per edge
/// `tau_g = W(a, N x a) = k12 cos 2phi - (k11 - k22) / 2 sin 2phi`, as the
/// area mean and relative to the anisotropy (`sum A tau_g / sum A |z|`, the
/// mean `sin 2 delta`: 0 aligned, +-1 at 45 degrees).
struct torsion_reading {
  real mean = 0;
  real relative = 0;
};

inline torsion_reading read_torsion(const kusama::stencil_table &t, const std::vector<vec3> &x,
                                    const std::vector<vec3> &h, const axis_field &axis) {
  const std::vector<real> A = detail::edge_dual_area(t, x);
  const auto fr = hepworth::block::build_edge_frames(t, x);
  real st = 0, sz = 0, sa = 0;
  for (int e = 0; e < t.n_edges; ++e) {
    const size_t es = static_cast<size_t>(e);
    const vec3 &he = h[es];
    if (!he.allFinite())
      continue;
    const Eigen::RowVector3d r =
        detail::torsion_row(*fr, es, axis(detail::edge_midpoint(t, x, e), fr->n[es]));
    if (!(r.squaredNorm() > real(0)))
      continue;
    st += A[es] * r.dot(he);
    sz += A[es] * std::hypot((he[0] - he[2]) / 2, he[1]);
    sa += A[es];
  }
  return {.mean = sa > real(0) ? st / sa : real(0), .relative = sz > real(0) ? st / sz : real(0)};
}

// ------------------------------------------------------- twisted spin-2

/// One triangle's twisted traceless spin-2 element (9 x 9 over its three edge
/// tensors): the CR element `K_ab (R_a M_a)^T Q_z (R_b M_b)` with each edge's
/// `z` demodulated by `R_a = rotate_z(-omega a . (m_a - c))`, `m_a` the edge
/// midpoint, `c` the centroid, `a` the unit axis in the face plane. At
/// `omega = 0` it is `kusama::weingarten::spin2_traceless_face_element`.
inline Eigen::Matrix<real, 9, 9> twisted_spin2_face_element(const vec3 &X0, const vec3 &X1, const vec3 &X2,
                                                            const vec3 &axis, real omega) {
  namespace wg = kusama::weingarten;
  const vec3 X[3] = {X0, X1, X2};
  const Eigen::Matrix<real, 3, 3> K = wg::cr_face_stiffness(X0, X1, X2);
  const Eigen::Matrix<real, 3, 3> Q = wg::tensor_inner_traceless();
  vec3 T, B, N;
  detail::face_frame(X0, X1, X2, T, B, N);
  vec3 a = axis - N * axis.dot(N);
  a = a.norm() > real(0) ? vec3(a.normalized()) : vec3::Zero();
  const vec3 c = (X0 + X1 + X2) / real(3);
  std::array<Eigen::Matrix<real, 3, 3>, 3> M;
  for (int i = 0; i < 3; ++i) {
    const vec3 v = (X[(i + 1) % 3] - X[i]).normalized();
    const vec3 m = real(0.5) * (X[i] + X[(i + 1) % 3]);
    M[static_cast<size_t>(i)] = detail::rotate_z(-omega * a.dot(m - c)) * wg::tensor_map(v.dot(T), v.dot(B));
  }
  Eigen::Matrix<real, 9, 9> E;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      E.block<3, 3>(3 * i, 3 * j) =
          K(i, j) * M[static_cast<size_t>(i)].transpose() * Q * M[static_cast<size_t>(j)];
  return E;
}

/// `w int |grad z - i omega a z|^2` over the store, rows `sqrt(w E_f)`, target
/// zero. Fixed for the step (geometry and axis at construction).
class twisted_spin2 : public hepworth::projection_constraint {
public:
  using ptr = std::shared_ptr<twisted_spin2>;
  static ptr create(const kusama::stencil_table &t, const std::vector<vec3> &x, const axis_field &axis,
                    real omega, real w, hepworth::sim_block::ptr store) {
    return std::make_shared<twisted_spin2>(t, x, axis, omega, w, std::move(store));
  }
  twisted_spin2(const kusama::stencil_table &t, const std::vector<vec3> &x, const axis_field &axis,
                real omega, real w, hepworth::sim_block::ptr store)
      : hepworth::projection_constraint({}, w), _store(std::move(store)) {
    const real sw = std::sqrt(std::max(w, real(0)));
    for (int f = 0; f < t.n_faces; ++f) {
      const auto &id = t.ids[static_cast<size_t>(f)];
      const vec3 &X0 = x[static_cast<size_t>(id[0])], &X1 = x[static_cast<size_t>(id[1])],
                 &X2 = x[static_cast<size_t>(id[2])];
      const vec3 Nf = (X1 - X0).cross(X2 - X0);
      if (!(Nf.norm() > real(0)))
        continue;
      const vec3 a = axis((X0 + X1 + X2) / real(3), Nf.normalized());
      const Eigen::Matrix<real, 9, 9> E = twisted_spin2_face_element(X0, X1, X2, a, omega);
      if (!E.allFinite())
        continue;
      const Eigen::Matrix<real, 9, 9> S = detail::psd_sqrt(E);
      for (int r = 0; r < 9; ++r)
        for (int i = 0; i < 3; ++i)
          for (int c = 0; c < 3; ++c) {
            const real v = sw * S(r, 3 * i + c);
            if (v != real(0))
              _rows.emplace_back(_n + r, 3 * t.face_edges[static_cast<size_t>(f)][static_cast<size_t>(i)] + c,
                                 v);
          }
      _n += 9;
    }
  }
  std::string name() override { return "twisted_spin2"; }

  void fill_A(index_t &id0, std::vector<trip> &triplets) override {
    _id0 = id0;
    for (const trip &t : _rows)
      triplets.emplace_back(id0 + t.row(), _store->get_offset_idx(t.col() / 3) + t.col() % 3, t.value());
    id0 += _n;
  }
  void project(const vecX &, vecX &) override {} // target zero

private:
  hepworth::sim_block::ptr _store;
  std::vector<trip> _rows;
  index_t _n = 0;
};

// ------------------------------------------------------- twist rate energy

/// The twist rate as a scale-free Gauss-Newton energy on the store: per face
/// `r_f = sqrt(w A_f) |z0_f| (omega_f(S) - omega)`, `omega_f =
/// Im(conj z d_a z) / |z|^2` from the live store (the readout of
/// `read_twist`), `|z0_f|` frozen at construction. `omega_f` does not change
/// when `z` is scaled, so the energy cannot be met by erasing the
/// anisotropy; the weight `|z0_f|` keeps near-umbilic faces (where the
/// principal directions are undefined) out, and `|z|^2` is regularized by
/// `(0.1 |z0_f|)^2`. Geometry (face frames, CR gradients, axis) is fixed for
/// the step.
class twist_rate : public hepworth::projection_constraint {
public:
  using ptr = std::shared_ptr<twist_rate>;
  static ptr create(const kusama::stencil_table &t, const std::vector<vec3> &x, const std::vector<vec3> &h0,
                    const axis_field &axis, real omega, real w, hepworth::sim_block::ptr store,
                    real z_floor = real(0.05)) {
    return std::make_shared<twist_rate>(t, x, h0, axis, omega, w, std::move(store), z_floor);
  }
  twist_rate(const kusama::stencil_table &t, const std::vector<vec3> &x, const std::vector<vec3> &h0,
             const axis_field &axis, real omega, real w, hepworth::sim_block::ptr store, real z_floor)
      : hepworth::projection_constraint({}, w), _omega(omega), _store(std::move(store)) {
    namespace wg = kusama::weingarten;
    Eigen::Matrix<real, 2, 3> P;
    P << real(0.5), 0, real(-0.5), 0, 1, 0;
    real zsum = 0, asum = 0;
    std::vector<face> all;
    for (int f = 0; f < t.n_faces; ++f) {
      const auto &id = t.ids[static_cast<size_t>(f)];
      const vec3 X[3] = {x[static_cast<size_t>(id[0])], x[static_cast<size_t>(id[1])],
                         x[static_cast<size_t>(id[2])]};
      vec3 T, B, N;
      detail::face_frame(X[0], X[1], X[2], T, B, N);
      const real A = detail::face_area(X[0], X[1], X[2]);
      if (!(A > real(0)) || !N.allFinite())
        continue;
      vec3 a = axis((X[0] + X[1] + X[2]) / real(3), N);
      a -= N * a.dot(N);
      if (!(a.norm() > real(0)))
        continue;
      a.normalize();
      face fc;
      fc.U.setZero();
      fc.D.setZero();
      for (int i = 0; i < 3; ++i) {
        fc.e[static_cast<size_t>(i)] = t.face_edges[static_cast<size_t>(f)][static_cast<size_t>(i)];
        const vec3 v = (X[(i + 1) % 3] - X[i]).normalized();
        const Eigen::Matrix<real, 2, 3> Z = P * wg::tensor_map(v.dot(T), v.dot(B));
        const vec3 grad = -N.cross(X[(i + 1) % 3] - X[i]) / A;
        fc.U.block<2, 3>(0, 3 * i) = Z / real(3);
        fc.D.block<2, 3>(0, 3 * i) = a.dot(grad) * Z;
      }
      Eigen::Matrix<real, 9, 1> hf;
      for (int i = 0; i < 3; ++i)
        hf.segment<3>(3 * i) = h0[static_cast<size_t>(fc.e[static_cast<size_t>(i)])];
      fc.z0 = (fc.U * hf).norm();
      fc.s = std::sqrt(std::max(w * A, real(0))) * fc.z0;
      zsum += A * fc.z0;
      asum += A;
      all.push_back(fc);
    }
    // Near-umbilic faces (principal directions undefined) get no row.
    const real zmin = z_floor * (asum > real(0) ? zsum / asum : real(0));
    for (const face &fc : all)
      if (fc.z0 > zmin && fc.s > real(0))
        _faces.push_back(fc);
  }

  /// Faces with a row.
  size_t face_count() const { return _faces.size(); }

  std::string name() override { return "twist_rate"; }
  bool gauss_newton() const override { return true; }

  void fill_J(const vecX &q, index_t &id0, std::vector<trip> &triplets) override {
    for (size_t k = 0; k < _faces.size(); ++k) {
      Eigen::Matrix<real, 1, 9> J;
      residual(_faces[k], q, &J);
      for (int i = 0; i < 3; ++i)
        for (int c = 0; c < 3; ++c)
          if (J(3 * i + c) != real(0))
            triplets.emplace_back(id0 + static_cast<index_t>(k),
                                  _store->get_offset_idx(_faces[k].e[static_cast<size_t>(i)]) + c,
                                  J(3 * i + c));
    }
    id0 += static_cast<index_t>(_faces.size());
  }
  void add_gradient(const vecX &q, vecX &g) override {
    for (const face &fc : _faces) {
      Eigen::Matrix<real, 1, 9> J;
      const real r = residual(fc, q, &J);
      for (int i = 0; i < 3; ++i)
        g.segment<3>(_store->get_offset_idx(fc.e[static_cast<size_t>(i)])) +=
            r * J.segment<3>(3 * i).transpose();
    }
  }
  real residual_sq(const vecX &q) override {
    real s = 0;
    for (const face &fc : _faces) {
      const real r = residual(fc, q, nullptr);
      s += r * r;
    }
    return s;
  }

private:
  struct face {
    std::array<int, 3> e{};
    Eigen::Matrix<real, 2, 9> U, D; ///< z and d_a z at the centroid, per store component
    real z0 = 0;                    ///< |z| at construction
    real s = 0;                     ///< sqrt(w A) |z0|
  };

  real residual(const face &fc, const vecX &q, Eigen::Matrix<real, 1, 9> *J) const {
    Eigen::Matrix<real, 9, 1> hf;
    for (int i = 0; i < 3; ++i)
      hf.segment<3>(3 * i) = _store->get_vec3(fc.e[static_cast<size_t>(i)], q);
    const real u = fc.U.row(0).dot(hf), v = fc.U.row(1).dot(hf);
    const real p = fc.D.row(0).dot(hf), qd = fc.D.row(1).dot(hf);
    // Regularized at 10% of the face's starting |z| (a 1% bias at the start):
    // a face passing near zero anisotropy mid-solve must not blow up.
    const real n2 = u * u + v * v + real(0.01) * fc.z0 * fc.z0;
    const real num = u * qd - v * p;
    const real om = num / n2;
    if (J) {
      const Eigen::Matrix<real, 1, 9> dnum =
          qd * fc.U.row(0) + u * fc.D.row(1) - p * fc.U.row(1) - v * fc.D.row(0);
      const Eigen::Matrix<real, 1, 9> dn2 = real(2) * (u * fc.U.row(0) + v * fc.U.row(1));
      *J = fc.s * (dnum / n2 - num * dn2 / (n2 * n2));
    }
    return fc.s * (om - _omega);
  }

  real _omega = 0;
  hepworth::sim_block::ptr _store;
  std::vector<face> _faces;
};

// ------------------------------------------------------- twist readout

struct twist_reading {
  std::vector<real> omega;  ///< per face, rad / length along the axis (NaN where |z| vanishes)
  std::vector<real> weight; ///< per face, A_f |z_f|
  real mean = 0;            ///< weighted mean of omega
  real mean_z = 0;          ///< area-weighted mean |z|
};

/// Per face: the three edge tensors in the face frame, `z` at the centroid
/// (their mean) and `d_a z` from the CR basis gradients, so
/// `omega_f = Im(conj z d_a z) / |z|^2`. Rows of `|z|` below `z_floor` times
/// the mean are skipped.
inline twist_reading read_twist(const kusama::stencil_table &t, const std::vector<vec3> &x,
                                const std::vector<vec3> &h, const axis_field &axis,
                                real z_floor = real(1e-3)) {
  namespace wg = kusama::weingarten;
  using cplx = std::complex<real>;
  twist_reading out;
  out.omega.assign(static_cast<size_t>(t.n_faces), std::numeric_limits<real>::quiet_NaN());
  out.weight.assign(static_cast<size_t>(t.n_faces), real(0));
  std::vector<cplx> zf(static_cast<size_t>(t.n_faces)), dz(static_cast<size_t>(t.n_faces));
  std::vector<real> Af(static_cast<size_t>(t.n_faces), real(0));
  real zsum = 0, asum = 0;
  for (int f = 0; f < t.n_faces; ++f) {
    const auto &id = t.ids[static_cast<size_t>(f)];
    const vec3 X[3] = {x[static_cast<size_t>(id[0])], x[static_cast<size_t>(id[1])],
                       x[static_cast<size_t>(id[2])]};
    vec3 T, B, N;
    detail::face_frame(X[0], X[1], X[2], T, B, N);
    const real A = detail::face_area(X[0], X[1], X[2]);
    if (!(A > real(0)) || !N.allFinite())
      continue;
    vec3 a = axis((X[0] + X[1] + X[2]) / real(3), N);
    a -= N * a.dot(N);
    if (!(a.norm() > real(0)))
      continue;
    a.normalize();
    cplx z0(0, 0), d(0, 0);
    for (int i = 0; i < 3; ++i) {
      const int e = t.face_edges[static_cast<size_t>(f)][static_cast<size_t>(i)];
      const vec3 v = (X[(i + 1) % 3] - X[i]).normalized();
      const vec3 hf = wg::tensor_map(v.dot(T), v.dot(B)) * h[static_cast<size_t>(e)];
      const cplx zi((hf[0] - hf[2]) / 2, hf[1]);
      // CR basis of local edge i (1 at its midpoint, -1 at the opposite vertex).
      const vec3 grad = -N.cross(X[(i + 1) % 3] - X[i]) / A;
      z0 += zi / real(3);
      d += zi * a.dot(grad);
    }
    zf[static_cast<size_t>(f)] = z0;
    dz[static_cast<size_t>(f)] = d;
    Af[static_cast<size_t>(f)] = A;
    zsum += A * std::abs(z0);
    asum += A;
  }
  out.mean_z = asum > real(0) ? zsum / asum : real(0);
  real wsum = 0, osum = 0;
  for (int f = 0; f < t.n_faces; ++f) {
    const cplx z0 = zf[static_cast<size_t>(f)];
    const real n2 = std::norm(z0);
    if (!(std::sqrt(n2) > z_floor * out.mean_z))
      continue;
    const real om = std::imag(std::conj(z0) * dz[static_cast<size_t>(f)]) / n2;
    out.omega[static_cast<size_t>(f)] = om;
    const real w = Af[static_cast<size_t>(f)] * std::sqrt(n2);
    out.weight[static_cast<size_t>(f)] = w;
    wsum += w;
    osum += w * om;
  }
  out.mean = wsum > real(0) ? osum / wsum : real(0);
  return out;
}

} // namespace engman
} // namespace gaudi

#endif
