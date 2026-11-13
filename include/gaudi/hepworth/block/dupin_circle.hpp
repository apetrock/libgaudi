#ifndef GAUDI_HEPWORTH_BLOCK_DUPIN_CIRCLE_HPP
#define GAUDI_HEPWORTH_BLOCK_DUPIN_CIRCLE_HPP

/// Position-space Dupin rows: four consecutive points along a curvature line
/// fitted to one planar circle and projected onto it (ShapeUp-style local
/// step, one projective row per group). Dupin cyclides are exactly the
/// surfaces whose curvature lines in both families are circles
/// (`reference/2312.04341v1.pdf`, Bobenko and Fairley, section 4); any three
/// points lie on a circle, so the constraint needs four, and four points
/// being concyclic is Moebius-invariant (their cross ratio is real, Vaxman et
/// al. 2015 section 4.1).
///
/// Groups: for an edge `(i, j)` within `max_angle` of one family's principal
/// direction at both ends, `a` the neighbour of `i` best aligned with `-t` and
/// `b` the neighbour of `j` best aligned with `+t` (each within `max_angle`);
/// weighted by the alignment and the endpoints' anisotropy (principal
/// directions are undefined at umbilics, so those vertices drop out).
///
/// Every iterate fits the live circle (least-squares plane, then the
/// algebraic circle fit in it; nearly collinear groups fit a line) and the
/// mode chooses each point's target on it:
///   closest     each live point to its closest point on the live circle
///               (the plain Dupin row);
///   slide       a rotation of the circle (elliptic Moebius): every point
///               advances `delta` mean edge lengths of arc per frame;
///   even        ends kept, the inner points at thirds of the arc;
///   graded      ends kept, the three arcs in ratio 1 : g : g^2;
///   hyperbolic  a hyperbolic Moebius map of the circle fixing the group's
///               mean point and its antipode: `x* = lambda^sigma x`,
///               `x = tan((theta - m0) / 2)`; keeps the cross ratio;
///   parabolic   a parabolic Moebius map fixing the antipode: `x* = x + sigma
///               beta`; keeps the cross ratio;
///   loxodromic  the circle moved on its family's curvature sphere at `i` by
///               the loxodromic map fixing `x_i` and its antipode
///               (stereographic `z* = lambda e^{i psi} z`): turns the
///               curvature line about `x_i` by `psi`, scales it toward `x_i`.
/// Drive modes (all but `closest`) compute seed points once, at
/// construction (the frame start), and each iterate targets the seeds'
/// closest points on the live circle: a drive acts once per frame and does
/// not drift, while the live refit keeps the line circular. `sigma` (the
/// direction along the line) comes from a reference tangent field; modes
/// that need it skip groups without one.

#include "gaudi/hepworth/block/block_constraint.hpp"
#include "gaudi/hepworth/block/shell_constraints.hpp"
#include "gaudi/kusama/stencil_table.hpp"

#include <Eigen/Eigenvalues>
#include <array>
#include <cmath>
#include <complex>
#include <vector>

namespace gaudi {
namespace hepworth {
namespace block {

struct dupin_circle_spec {
  enum class mode { closest, slide, even, graded, hyperbolic, parabolic, loxodromic };
  mode kind = mode::closest;
  /// 0 both families, 1 the tight one (larger |k|), 2 the loose one.
  int family = 0;
  real w = 1;
  real delta = 0.1;  ///< slide: arc per frame, in mean edge lengths
  real g = 1.2;      ///< graded: arc ratio
  real lambda = 1.1; ///< hyperbolic / loxodromic multiplier
  real beta = 0.05;  ///< parabolic shift (stereographic coordinate)
  real psi = 0.05;   ///< loxodromic turn per frame (radians)
  real max_angle = real(25.0 * M_PI / 180.0);
  real min_aniso = 0.05;
};

namespace dupin_circle_detail {

struct fit4 {
  vec3 m = vec3::Zero();   ///< centroid
  vec3 c = vec3::Zero();   ///< circle centre
  vec3 nu = vec3::UnitZ(); ///< plane normal
  vec3 dir = vec3::UnitX(); ///< line direction (largest principal axis)
  real r = 0;
  bool line = false;
  bool ok = false;
};

inline fit4 fit(const std::array<vec3, 4> &p, real r_max = real(1e3)) {
  fit4 f;
  f.m = real(0.25) * (p[0] + p[1] + p[2] + p[3]);
  mat3 C = mat3::Zero();
  real size = 0;
  for (const vec3 &q : p) {
    C += (q - f.m) * (q - f.m).transpose();
    size = std::max(size, (q - f.m).norm());
  }
  if (!(size > real(0)) || !C.allFinite())
    return f;
  Eigen::SelfAdjointEigenSolver<mat3> es(C);
  f.nu = es.eigenvectors().col(0);
  f.dir = es.eigenvectors().col(2);
  const vec3 e1 = f.dir, e2 = f.nu.cross(e1);
  Eigen::Matrix<real, 4, 3> A;
  Eigen::Matrix<real, 4, 1> b;
  for (int k = 0; k < 4; ++k) {
    const real x = (p[k] - f.m).dot(e1), y = (p[k] - f.m).dot(e2);
    A.row(k) << x, y, real(1);
    b[k] = -(x * x + y * y);
  }
  const vec3 s = A.colPivHouseholderQr().solve(b);
  const real cx = -s[0] / 2, cy = -s[1] / 2, r2 = cx * cx + cy * cy - s[2];
  f.ok = true;
  if (!(r2 > real(0)) || !std::isfinite(r2) || std::sqrt(r2) > r_max * size) {
    f.line = true;
    return f;
  }
  f.r = std::sqrt(r2);
  f.c = f.m + cx * e1 + cy * e2;
  return f;
}

inline vec3 closest(const fit4 &f, const vec3 &q) {
  if (f.line)
    return f.m + f.dir * (q - f.m).dot(f.dir);
  vec3 d = q - f.nu * (q - f.c).dot(f.nu) - f.c;
  const real l = d.norm();
  if (!(l > real(1e-300)))
    return q;
  return f.c + f.r * d / l;
}

/// Mean distance of the four points to their fitted circle.
inline real residual(const std::array<vec3, 4> &p) {
  const fit4 f = fit(p);
  if (!f.ok)
    return 0;
  real s = 0;
  for (const vec3 &q : p)
    s += (closest(f, q) - q).norm();
  return s / 4;
}

/// The loxodromic map of the sphere (centre `C`, radius `R`) fixing `x` (on
/// it) and its antipode, in stereographic coordinates from the antipode on
/// the tangent plane at `x` with basis `(e, n x e)`: `z* = lambda e^{i psi} z`.
inline vec3 loxodromic(const vec3 &p, const vec3 &C, real R, const vec3 &x, const vec3 &n, real lambda, real psi) {
  const vec3 m = (x - C).normalized();
  const vec3 A = C - R * m;
  vec3 e = (std::abs(m[0]) < 0.9 ? vec3::UnitX() : vec3::UnitY());
  e = (e - m * e.dot(m)).normalized();
  const vec3 f = n.cross(e).normalized();
  const real den = (p - A).dot(m);
  if (!(std::abs(den) > real(1e-300)))
    return p;
  const vec3 z = A + (2 * R / den) * (p - A);
  std::complex<real> w((z - x).dot(e), (z - x).dot(f));
  w *= lambda * std::exp(std::complex<real>(0, psi));
  const vec3 zp = x + w.real() * e + w.imag() * f;
  const vec3 d = zp - A;
  const real dd = d.squaredNorm();
  if (!(dd > real(1e-300)))
    return p;
  return A + (2 * R * m.dot(d) / dd) * d;
}

} // namespace dupin_circle_detail

/// One group `(a, i, j, b)`. Rows: the four positions minus their centroid,
/// times `w s` (12 rows); targets the mode's points on the live circle.
class circle4 : public block_constraint {
public:
  DEFINE_CREATE_FUNC(circle4)

  /// `x0`: frame-start positions; `s`: the group's weight (alignment times
  /// anisotropy); `sigma`: the direction along the line (+1, -1, or 0 when
  /// unknown); `sphere_c`, `sphere_r`, `n_i`: the family's curvature sphere
  /// at `i` and the normal there (loxodromic only).
  circle4(const std::array<index_t, 4> &ids, const std::vector<vec3> &x0, const dupin_circle_spec &spec, real s,
          int sigma, const vec3 &sphere_c, real sphere_r, const vec3 &n_i, std::vector<sim_block::ptr> blocks)
      : block_constraint(std::vector<index_t>(ids.begin(), ids.end()), spec.w, std::move(blocks)),
        _s(s), _kind(spec.kind) {
    using namespace dupin_circle_detail;
    std::array<vec3, 4> p;
    for (int k = 0; k < 4; ++k)
      p[k] = x0[static_cast<size_t>(ids[k])];
    _rest_initialized = true;
    if (_kind == dupin_circle_spec::mode::closest)
      return;
    if (!seed(p, spec, sigma, sphere_c, sphere_r, n_i))
      _kind = dupin_circle_spec::mode::closest;
  }

  std::string name() override { return "dupin_circle"; }

  void project(const vecX &q, vecX &p) override {
    using namespace dupin_circle_detail;
    std::array<vec3, 4> x;
    for (int k = 0; k < 4; ++k)
      x[k] = _blocks[0]->get_vec3(_ids[k], q);
    const fit4 f = fit(x);
    std::array<vec3, 4> t = x;
    if (f.ok)
      for (int k = 0; k < 4; ++k)
        t[k] = closest(f, _kind == dupin_circle_spec::mode::closest ? x[k] : _seed[k]);
    const vec3 tm = real(0.25) * (t[0] + t[1] + t[2] + t[3]);
    // The circle condition is scale-free, so the row keeps the group's size:
    // the targets rescaled about their centroid to the live points' RMS
    // spread. Without it the algebraic fit's small-radius bias on short arcs
    // shrinks the surface frame after frame.
    const vec3 xm = real(0.25) * (x[0] + x[1] + x[2] + x[3]);
    real sx = 0, st = 0;
    for (int k = 0; k < 4; ++k) {
      sx += (x[k] - xm).squaredNorm();
      st += (t[k] - tm).squaredNorm();
    }
    const real scale = st > real(0) ? std::sqrt(sx / st) : real(1);
    for (int k = 0; k < 4; ++k)
      p.block(_id0 + 3 * k, 0, 3, 1) = _w * _s * scale * (t[k] - tm);
  }

  void fill_A(index_t &id0, std::vector<trip> &triplets) override {
    _id0 = id0;
    for (int k = 0; k < 4; ++k)
      for (int l = 0; l < 4; ++l) {
        const real c = _w * _s * ((k == l ? real(1) : real(0)) - real(0.25));
        const index_t col = _blocks[0]->get_offset_idx(_ids[l]);
        for (int ax = 0; ax < 3; ++ax)
          triplets.push_back(trip(id0 + 3 * k + ax, col + ax, c));
      }
    id0 += 12;
  }

  const std::array<vec3, 4> &seeds() const { return _seed; }
  dupin_circle_spec::mode kind() const { return _kind; }

private:
  /// Seed points from the frame-start fit; false falls back to `closest`.
  bool seed(const std::array<vec3, 4> &p, const dupin_circle_spec &spec, int sigma, const vec3 &sc, real sr,
            const vec3 &n_i) {
    using namespace dupin_circle_detail;
    using mode = dupin_circle_spec::mode;
    const fit4 f = fit(p);
    if (!f.ok || f.line)
      return false;
    std::array<vec3, 4> on;
    for (int k = 0; k < 4; ++k)
      on[k] = closest(f, p[k]);
    if (spec.kind == mode::loxodromic) {
      if (!(sr > real(0)) || !std::isfinite(sr))
        return false;
      const vec3 xi = sc + sr * (p[1] - sc).normalized();
      for (int k = 0; k < 4; ++k) {
        const vec3 onS = sc + sr * (on[k] - sc).normalized();
        _seed[k] = loxodromic(onS, sc, sr, xi, n_i, spec.lambda, spec.psi);
      }
      return _seed[0].allFinite() && _seed[1].allFinite() && _seed[2].allFinite() && _seed[3].allFinite();
    }
    // Angles on the frame-start circle from `i`, oriented so `j` follows.
    vec3 u = (on[1] - f.c).normalized(), nu = f.nu;
    vec3 v = nu.cross(u);
    auto ang = [&](const vec3 &q) { return std::atan2((q - f.c).dot(v), (q - f.c).dot(u)); };
    if (ang(on[2]) < 0) {
      nu = -nu;
      v = -v;
    }
    std::array<real, 4> th{ang(on[0]), real(0), ang(on[2]), ang(on[3])};
    if (th[3] < th[2])
      th[3] += 2 * real(M_PI);
    if (!(th[0] < 0 && th[2] > 0 && th[3] > th[2]))
      return false;
    std::array<real, 4> ts = th;
    const real span = th[3] - th[0];
    const real l = ((p[1] - p[0]).norm() + (p[2] - p[1]).norm() + (p[3] - p[2]).norm()) / 3;
    switch (spec.kind) {
    case mode::slide:
      if (sigma == 0)
        return false;
      for (real &t : ts)
        t += sigma * spec.delta * l / f.r;
      break;
    case mode::even:
      ts = {th[0], th[0] + span / 3, th[0] + 2 * span / 3, th[3]};
      break;
    case mode::graded: {
      if (sigma == 0)
        return false;
      const real g = std::max(spec.g, real(1e-6));
      const std::array<real, 3> a = sigma > 0 ? std::array<real, 3>{1, g, g * g} : std::array<real, 3>{g * g, g, 1};
      const real d = span / (a[0] + a[1] + a[2]);
      ts = {th[0], th[0] + a[0] * d, th[0] + (a[0] + a[1]) * d, th[3]};
      break;
    }
    case mode::hyperbolic:
    case mode::parabolic: {
      if (sigma == 0)
        return false;
      const real m0 = real(0.25) * (th[0] + th[1] + th[2] + th[3]);
      for (real &t : ts) {
        real x = std::tan((t - m0) / 2);
        x = spec.kind == mode::hyperbolic ? std::pow(spec.lambda, real(sigma)) * x : x + sigma * spec.beta;
        t = m0 + 2 * std::atan(x);
      }
      break;
    }
    default:
      return false;
    }
    for (int k = 0; k < 4; ++k)
      _seed[k] = f.c + f.r * (std::cos(ts[k]) * u + std::sin(ts[k]) * v);
    return true;
  }

  real _s = 1;
  dupin_circle_spec::mode _kind = dupin_circle_spec::mode::closest;
  std::array<vec3, 4> _seed{};
};

/// A curvature-line group and its data.
struct dupin_group {
  std::array<index_t, 4> ids{};
  int family = 1;
  real s = 0;
  int sigma = 0;
};

/// The groups of `x` for `family` (1 tight, 2 loose), from the vertex charts.
/// `ref` (optional, per vertex): the reference tangent field for `sigma`.
inline std::vector<dupin_group> select_dupin_groups(const kusama::stencil_table &t, const std::vector<vec3> &x,
                                                    const std::vector<vertex_chart> &charts, int family,
                                                    real max_angle, real min_aniso,
                                                    const std::vector<vec3> *ref = nullptr) {
  const size_t n = x.size();
  std::vector<std::vector<int>> nbr(n);
  std::vector<std::array<int, 2>> edges(static_cast<size_t>(t.n_edges), {-1, -1});
  for (int f = 0; f < t.n_faces; ++f) {
    const auto &id = t.ids[static_cast<size_t>(f)];
    for (int k = 0; k < 3; ++k) {
      const int a = id[static_cast<size_t>(k)], b = id[static_cast<size_t>((k + 1) % 3)];
      if (a < 0 || b < 0 || static_cast<size_t>(a) >= n || static_cast<size_t>(b) >= n)
        continue;
      const int e = t.face_edges[static_cast<size_t>(f)][static_cast<size_t>(k)];
      if (e >= 0 && e < t.n_edges)
        edges[static_cast<size_t>(e)] = {a, b};
      if (std::find(nbr[a].begin(), nbr[a].end(), b) == nbr[a].end())
        nbr[a].push_back(b);
      if (std::find(nbr[b].begin(), nbr[b].end(), a) == nbr[b].end())
        nbr[b].push_back(a);
    }
  }
  std::vector<vec3> dir(n, vec3::Zero()), N(n, vec3::Zero());
  std::vector<real> q(n, 0);
  for (size_t v = 0; v < n && v < charts.size(); ++v) {
    shape_w W;
    vec3 T, B, Nv;
    if (!vertex_chart_at(&charts, static_cast<index_t>(v), W, T, B, Nv))
      continue;
    const shape_w ax = shape_with_axes(W);
    if (!ax.ok)
      continue;
    const monge_axes m = monge_axes_of(ax);
    dir[v] = family == 1 ? vec3(m.hx * T + m.hy * B) : vec3(m.lx * T + m.ly * B);
    N[v] = Nv.normalized();
    q[v] = std::abs(m.k_hi - m.k_lo) / std::max(std::abs(m.k_hi) + std::abs(m.k_lo), real(1e-300));
  }
  const real ca = std::cos(max_angle);
  auto unit_t = [&](int a, int b, int at) {
    vec3 d = x[static_cast<size_t>(b)] - x[static_cast<size_t>(a)];
    d -= N[static_cast<size_t>(at)] * d.dot(N[static_cast<size_t>(at)]);
    const real l = d.norm();
    return l > real(0) ? vec3(d / l) : vec3::Zero();
  };
  std::vector<dupin_group> out;
  for (const auto &e : edges) {
    const int i = e[0], j = e[1];
    if (i < 0 || q[static_cast<size_t>(i)] < min_aniso || q[static_cast<size_t>(j)] < min_aniso)
      continue;
    const vec3 ti = unit_t(i, j, i), tj = unit_t(i, j, j);
    const real ali = std::abs(ti.dot(dir[static_cast<size_t>(i)])), alj = std::abs(tj.dot(dir[static_cast<size_t>(j)]));
    if (ali < ca || alj < ca)
      continue;
    int a = -1, b = -1;
    real best_a = ca, best_b = ca;
    for (int c : nbr[static_cast<size_t>(i)]) {
      if (c == j)
        continue;
      const real s = unit_t(i, c, i).dot(-ti);
      if (s > best_a) {
        best_a = s;
        a = c;
      }
    }
    for (int c : nbr[static_cast<size_t>(j)]) {
      if (c == i)
        continue;
      const real s = unit_t(j, c, j).dot(tj);
      if (s > best_b) {
        best_b = s;
        b = c;
      }
    }
    if (a < 0 || b < 0 || a == b)
      continue;
    dupin_group g;
    g.ids = {a, i, j, b};
    g.family = family;
    g.s = real(0.5) * (ali + alj) * std::min(q[static_cast<size_t>(i)], q[static_cast<size_t>(j)]);
    if (ref && static_cast<size_t>(i) < ref->size()) {
      const real d = ti.dot((*ref)[static_cast<size_t>(i)]);
      g.sigma = std::abs(d) > real(0.1) ? (d > 0 ? 1 : -1) : 0;
    }
    out.push_back(g);
  }
  return out;
}

/// The rows of every spec. Returns the number of groups added.
inline int init_dupin_circles(const kusama::stencil_table &t, const std::vector<vec3> &x,
                              const std::vector<vertex_chart> &charts, const std::vector<dupin_circle_spec> &specs,
                              real scale, const std::vector<vec3> *ref,
                              std::vector<projection_constraint::ptr> &constraints,
                              std::vector<sim_block::ptr> blocks) {
  int count = 0;
  if (!(scale > real(0)))
    return 0;
  for (dupin_circle_spec spec : specs) {
    spec.w *= scale;
    if (!(spec.w > real(0)))
      continue;
    for (int fam = 1; fam <= 2; ++fam) {
      if (spec.family != 0 && spec.family != fam)
        continue;
      for (const dupin_group &g : select_dupin_groups(t, x, charts, fam, spec.max_angle, spec.min_aniso, ref)) {
        vec3 sc = vec3::Zero(), ni = vec3::Zero();
        real sr = 0;
        if (spec.kind == dupin_circle_spec::mode::loxodromic) {
          shape_w W;
          vec3 T, B, Nv;
          if (!vertex_chart_at(&charts, g.ids[1], W, T, B, Nv))
            continue;
          const monge_axes m = monge_axes_of(shape_with_axes(W));
          const real k = fam == 1 ? m.k_hi : m.k_lo;
          if (!(std::abs(k) > real(1e-6)))
            continue;
          ni = Nv.normalized();
          sc = x[static_cast<size_t>(g.ids[1])] + ni / k;
          sr = real(1) / std::abs(k);
        }
        constraints.push_back(circle4::create(g.ids, x, spec, g.s, g.sigma, sc, sr, ni, blocks));
        ++count;
      }
    }
  }
  return count;
}

/// Mean circle residual of the curvature-line groups (both families), over
/// the mean edge length: how far the curvature lines are from circles.
inline real dupin_circle_residual(const kusama::stencil_table &t, const std::vector<vec3> &x,
                                  const std::vector<vertex_chart> &charts, real max_angle = real(25.0 * M_PI / 180.0),
                                  real min_aniso = real(0.05)) {
  real sum = 0, len = 0;
  size_t n = 0, ne = 0;
  for (int fam = 1; fam <= 2; ++fam)
    for (const dupin_group &g : select_dupin_groups(t, x, charts, fam, max_angle, min_aniso)) {
      std::array<vec3, 4> p;
      for (int k = 0; k < 4; ++k)
        p[k] = x[static_cast<size_t>(g.ids[k])];
      sum += dupin_circle_detail::residual(p);
      ++n;
      len += (p[2] - p[1]).norm();
      ++ne;
    }
  if (n == 0 || !(len > real(0)))
    return std::numeric_limits<real>::quiet_NaN();
  return (sum / real(n)) / (len / real(ne));
}

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif
