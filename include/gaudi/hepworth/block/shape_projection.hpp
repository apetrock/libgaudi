#ifndef GAUDI_HEPWORTH_SHAPE_PROJECTION_HPP
#define GAUDI_HEPWORTH_SHAPE_PROJECTION_HPP

#include <vector>

#include "gaudi/hepworth/block/cr_dirichlet.hpp"
#include "gaudi/hepworth/block/edge_shape_fair.hpp"
#include "gaudi/hepworth/block/shape_objectives.hpp"
#include "gaudi/hepworth/block/shape_operator_datum.hpp"
#include "gaudi/hepworth/block/shell_constraints.hpp"
#include "gaudi/hepworth/block/shell_constraints_init.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

inline int vec3_block_index(const std::vector<sim_block::ptr> &blocks, int skip) {
  for (int i = 0; i < static_cast<int>(blocks.size()); ++i) {
    if (i == skip)
      continue;
    if (std::dynamic_pointer_cast<vec3_block>(blocks[static_cast<size_t>(i)]))
      return i;
  }
  return -1;
}

struct shape_fan {
  index_t i = -1;
  std::vector<fan_face> faces;
  std::vector<fan_face> support;
  std::vector<fan_edge> edges;
  std::vector<int> tangent;
};

/// Calls `row(ids, weights, i)` for each live vertex ring: `ids` the vertex
/// then its spokes' far ends, `weights` the cotan spoke weights (positive
/// ones only). Nothing without dense-packed vertices.
template <typename RowFn>
inline void for_each_spoke_ring(const asawa::shell::shell &shell, const std::vector<vec3> &x,
                                const std::vector<char> &valid, RowFn row) {
  if (!shell.verts_are_dense_packed())
    return;
  for (const fair_ring &ring : build_frame_rings(shell, x, valid)) {
    if (ring.i < 0 || ring.i >= static_cast<int>(x.size()))
      continue;
    std::vector<index_t> ids{ring.i};
    std::vector<real> weights;
    for (const fair_spoke &sp : ring.spokes) {
      if (!(sp.K > real(0)) || !sp.ok || sp.j < 0 || sp.j >= static_cast<int>(x.size()))
        continue;
      ids.push_back(sp.j);
      weights.push_back(sp.K);
    }
    if (!weights.empty())
      row(ids, weights, ring.i);
  }
}

/// Entry `i` of `field` (zero when `field` is null or short).
inline vec3 field_at(const std::vector<vec3> *field, int i) {
  return field && i >= 0 && static_cast<size_t>(i) < field->size() ? (*field)[static_cast<size_t>(i)]
                                                                   : vec3::Zero();
}

/// Position ARAP. Each entry is one row. `angle` 0 puts `x` on the minor
/// tangent and `y` on the major. Those are the live principal directions, the
/// same world vectors as the fairing tube and ruling, rebuilt every
/// projection. Positive `angle` turns from minor toward major. `x` and `y`
/// are extra fractions on those axes, each multiplied by `|κ|^power` of that
/// axis. `power` 0 leaves the fractions. `dirs[t]`, when set, is target
/// `t`'s stretch axis per vertex (world, held for the step) in place of its
/// `angle`. (0, 0) adds no row. `w` 0 or an empty list leaves the row off.
/// Columns are the position block only.
inline void init_directed_stretch(const asawa::shell::shell &shell,
                                 std::vector<projection_constraint::ptr> &constraints,
                                 const std::vector<vec3> &x,
                                 const std::vector<char> &valid,
                                 std::vector<sim_block::ptr> blocks, real w,
                                 const std::vector<stretch_target> &targets,
                                 const std::vector<vertex_chart> *charts = nullptr,
                                 const std::vector<const std::vector<vec3> *> &dirs = {}) {
  if (!(w > real(0)) || !std::isfinite(w) || targets.empty() || !charts || blocks.empty() ||
      !blocks[0])
    return;
  const std::vector<sim_block::ptr> held{blocks[0]};
  for_each_spoke_ring(shell, x, valid, [&](const std::vector<index_t> &ids,
                                           const std::vector<real> &weights, int i) {
    for (size_t t = 0; t < targets.size(); ++t) {
      const stretch_target &spec = targets[t];
      if (!std::isfinite(spec.x) || !std::isfinite(spec.y))
        continue;
      if (!(std::abs(spec.x) > real(1e-16) || std::abs(spec.y) > real(1e-16)))
        continue;
      constraints.push_back(directed_stretch::create(
          ids, weights, x, fan_local{}, 0, spec.angle, spec.x, spec.y, spec.power, w, held,
          charts, field_at(t < dirs.size() ? dirs[t] : nullptr, i)));
    }
  });
}

/// `directed_tilt` per vertex and target (`tilt_target`), cotan spoke
/// weights. A `rigid` / `propagate` target shares one `tilt_orientation`
/// over the one-ring graph across its rows. A zero gain adds no row. `w` 0 or
/// an empty list leaves the row off. Columns are the position block only.
inline void init_directed_tilt(const asawa::shell::shell &shell,
                               std::vector<projection_constraint::ptr> &constraints,
                               const std::vector<vec3> &x, const std::vector<char> &valid,
                               std::vector<sim_block::ptr> blocks, real w,
                               const std::vector<tilt_target> &targets,
                               const std::vector<vertex_chart> *charts) {
  if (!(w > real(0)) || !std::isfinite(w) || targets.empty() || !charts)
    return;
  if (!shell.verts_are_dense_packed() || blocks.empty() || !blocks[0])
    return;
  std::vector<fair_ring> rings = build_frame_rings(shell, x, valid);
  std::vector<std::shared_ptr<tilt_orientation>> fields(targets.size());
  for (size_t t = 0; t < targets.size(); ++t) {
    if (targets[t].mode != tilt_mode::rigid || targets[t].sign != tilt_sign::propagate)
      continue;
    auto f = std::make_shared<tilt_orientation>();
    f->axis = targets[t].axis;
    f->nbr.assign(x.size(), {});
    for (const fair_ring &ring : rings)
      if (ring.i >= 0 && ring.i < static_cast<int>(x.size()))
        for (const fair_spoke &sp : ring.spokes)
          if (sp.j >= 0 && sp.j < static_cast<int>(x.size()))
            f->nbr[static_cast<size_t>(ring.i)].push_back(sp.j);
    fields[t] = f;
  }
  std::vector<sim_block::ptr> held{blocks[0]};
  for (const fair_ring &ring : rings) {
    if (ring.i < 0 || ring.i >= static_cast<int>(x.size()))
      continue;
    std::vector<index_t> ids{ring.i};
    std::vector<real> weights;
    for (const fair_spoke &sp : ring.spokes) {
      if (!(sp.K > real(0)) || !sp.ok || sp.j < 0 || sp.j >= static_cast<int>(x.size()))
        continue;
      ids.push_back(sp.j);
      weights.push_back(sp.K);
    }
    if (weights.empty())
      continue;
    for (size_t t = 0; t < targets.size(); ++t) {
      const tilt_target &spec = targets[t];
      if (!std::isfinite(spec.gain) || spec.gain == real(0))
        continue;
      constraints.push_back(
          directed_tilt::create(ids, weights, x, spec, w, held, charts, fields[t]));
    }
  }
}

/// World direction per vertex that a row reads from `angle`: from the
/// chart's k_hi tangent toward its k_lo tangent (`monge_axes_of`; NaN is 45
/// degrees). The sense of that turn about the normal follows the
/// eigenvector order, so it is mirrored where the larger-magnitude
/// curvature is the smaller eigenvalue. Unit, sign per vertex; zero without
/// a chart.
inline std::vector<vec3> chart_directions(const std::vector<vertex_chart> &charts, real angle) {
  const real a = std::isfinite(angle) ? angle : real(0.25) * real(M_PI);
  std::vector<vec3> out(charts.size(), vec3::Zero());
  for (size_t i = 0; i < charts.size(); ++i) {
    shape_w W;
    vec3 T, B, N;
    if (!vertex_chart_at(&charts, static_cast<index_t>(i), W, T, B, N))
      continue;
    const shape_w axes = shape_with_axes(W);
    if (!axes.ok)
      continue;
    const monge_axes ax = monge_axes_of(axes);
    vec3 u = std::cos(a) * (ax.hx * T + ax.hy * B) + std::sin(a) * (ax.lx * T + ax.ly * B);
    const vec3 n = N.normalized();
    u -= n * u.dot(n);
    if (!(u.squaredNorm() > real(1e-16)) || !u.allFinite())
      continue;
    out[i] = u.normalized();
  }
  return out;
}

/// The angle from the tube circle (the tight tangent) of a loxodrome winding
/// `w` on a torus of radii `R > r`: `handedness * atan2(ring * sqrt(R^2 -
/// r^2), tube * r)`, the line straight in the torus's conformal chart that
/// closes after `ring` turns around the ring and `tube` around the tube.
/// `fallback` without windings or without `R > r > 0`.
inline real loxodrome_angle(const stretch_winding &w, real R, real r, real fallback) {
  if ((w.ring == 0 && w.tube == 0) || !(r > real(0)) || !(R > r) || !std::isfinite(R))
    return fallback;
  const real hand = w.handedness < 0 ? real(-1) : real(1);
  return hand * std::atan2(real(w.ring) * std::sqrt(R * R - r * r), real(w.tube) * r);
}

/// `axis_ratio` per vertex and target: the ring's spoke rms along k_hi over
/// that across held at `ratio`, area within the target's band, skew and
/// rotation free. A
/// non-positive ratio adds no row. `w` 0 or an empty list leaves the row
/// off. Columns are the position block only.
inline void init_axis_ratio(const asawa::shell::shell &shell,
                            std::vector<projection_constraint::ptr> &constraints,
                            const std::vector<vec3> &x, const std::vector<char> &valid,
                            std::vector<sim_block::ptr> blocks, real w,
                            const std::vector<axis_ratio_target> &targets,
                            const std::vector<vertex_chart> *charts) {
  if (!(w > real(0)) || !std::isfinite(w) || targets.empty() || !charts || blocks.empty() ||
      !blocks[0])
    return;
  const std::vector<sim_block::ptr> held{blocks[0]};
  for_each_spoke_ring(shell, x, valid, [&](const std::vector<index_t> &ids,
                                           const std::vector<real> &weights, int) {
    for (const axis_ratio_target &spec : targets)
      if (spec.ratio > real(0) && std::isfinite(spec.ratio))
        constraints.push_back(axis_ratio::create(ids, weights, x, spec.ratio, true, real(0), w,
                                                 held, charts, "axis_ratio", spec.area));
  });
}

/// `axis_ratio` per vertex and target as a shear: the frame-start ring
/// stretched by `(1 + shear) / (1 - shear)` along `angle` relative to across
/// it, area, that axis's skew and rotation free. Zero shear adds no row.
/// `w` 0 or an empty list leaves the row off. Columns are the position block
/// only.
inline void init_shear(const asawa::shell::shell &shell,
                       std::vector<projection_constraint::ptr> &constraints,
                       const std::vector<vec3> &x, const std::vector<char> &valid,
                       std::vector<sim_block::ptr> blocks, real w,
                       const std::vector<shear_target> &targets,
                       const std::vector<vertex_chart> *charts) {
  if (!(w > real(0)) || !std::isfinite(w) || targets.empty() || !charts || blocks.empty() ||
      !blocks[0])
    return;
  const std::vector<sim_block::ptr> held{blocks[0]};
  for_each_spoke_ring(shell, x, valid, [&](const std::vector<index_t> &ids,
                                           const std::vector<real> &weights, int) {
    for (const shear_target &spec : targets) {
      if (!std::isfinite(spec.shear) || spec.shear == real(0))
        continue;
      const real g = std::clamp(spec.shear, real(-0.99), real(0.99));
      const real angle = std::isfinite(spec.angle) ? spec.angle : real(0.25) * real(M_PI);
      constraints.push_back(axis_ratio::create(ids, weights, x, (real(1) + g) / (real(1) - g),
                                               false, angle, w, held, charts, "shear",
                                               spec.area));
    }
  });
}

/// The in-plane stretch of one triangle whose curvature tensor `S` (2x2,
/// face basis) is to become `Sh^-T S Sh^-1`, `Sh` the area-preserving shear
/// by `sqrt((1 + g) / (1 - g))` along the axis at `angle` from `S`'s k_hi
/// direction toward its k_lo one. With the edge normals held, an in-plane
/// map `F` of the triangle takes its CR tensor to `S F^-1`, so
/// `F = S*^-1 S`; the result is `F`'s stretch `U` (`F = R U`), `det U = 1`.
/// False near parabolic (`|det S|` under `min_det` of `|S|^2`) or when `F`
/// flips the triangle.
inline bool face_shear_stretch(const Eigen::Matrix<real, 2, 2> &S, real g, real angle,
                               Eigen::Matrix<real, 2, 2> &U, real min_det = real(1e-2),
                               Eigen::Matrix<real, 2, 2> *F_out = nullptr,
                               Eigen::Matrix<real, 2, 2> *S_star_out = nullptr) {
  using m2 = Eigen::Matrix<real, 2, 2>;
  using v2 = Eigen::Matrix<real, 2, 1>;
  const m2 Ss = real(0.5) * (S + S.transpose());
  if (!Ss.allFinite() || !(std::abs(Ss.determinant()) > min_det * Ss.squaredNorm()))
    return false;
  Eigen::SelfAdjointEigenSolver<m2> es(Ss);
  const bool hi0 = std::abs(es.eigenvalues()[0]) >= std::abs(es.eigenvalues()[1]);
  const v2 h = es.eigenvectors().col(hi0 ? 0 : 1), l = es.eigenvectors().col(hi0 ? 1 : 0);
  const v2 a = std::cos(angle) * h + std::sin(angle) * l, ap(-a[1], a[0]);
  g = std::clamp(g, real(-0.99), real(0.99));
  const real sigma = std::sqrt((real(1) + g) / (real(1) - g));
  const m2 Sh_inv = a * a.transpose() / sigma + sigma * ap * ap.transpose();
  const m2 S_star = Sh_inv.transpose() * Ss * Sh_inv;
  const m2 F = S_star.inverse() * Ss;
  if (!F.allFinite() || !(F.determinant() > real(0)))
    return false;
  if (F_out)
    *F_out = F;
  if (S_star_out)
    *S_star_out = S_star;
  Eigen::JacobiSVD<m2> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
  U = svd.matrixV() * svd.singularValues().asDiagonal() * svd.matrixV().transpose();
  U /= std::sqrt(U.determinant());
  return U.allFinite();
}

/// `face_stretch_target` per triangle and target: the face's curvature
/// tensor, the mean of its three edges' frame-start store tensors `h0` in
/// the face plane, sheared per `curvature_shear_target` (`angle` resolved),
/// and the triangle stretched so its CR tensor follows (`face_shear_stretch`).
/// Faces near parabolic add no row. `w` 0 or an empty list leaves the row
/// off. Columns are the position block only.
inline void init_face_shear(const kusama::stencil_table &tb, const edge_frames &fr,
                            const std::vector<vec3> &h0,
                            std::vector<projection_constraint::ptr> &constraints,
                            const std::vector<vec3> &x, std::vector<sim_block::ptr> blocks, real w,
                            const std::vector<curvature_shear_target> &targets) {
  if (!(w > real(0)) || !std::isfinite(w) || targets.empty() || blocks.empty() || !blocks[0] ||
      static_cast<int>(h0.size()) != tb.n_edges || static_cast<int>(fr.t.size()) != tb.n_edges)
    return;
  const std::vector<sim_block::ptr> held{blocks[0]};
  for (int f = 0; f < tb.n_faces; ++f) {
    const auto &id = tb.ids[static_cast<size_t>(f)];
    const auto &fe = tb.face_edges[static_cast<size_t>(f)];
    if (id[0] < 0 || id[1] < 0 || id[2] < 0)
      continue;
    const vec3 X0 = x[static_cast<size_t>(id[0])];
    const vec3 e1 = x[static_cast<size_t>(id[1])] - X0, e2 = x[static_cast<size_t>(id[2])] - X0;
    const vec3 n = e1.cross(e2);
    if (!(n.squaredNorm() > real(1e-30)))
      continue;
    mat3 T = mat3::Zero();
    for (int i = 0; i < 3; ++i) {
      const size_t e = static_cast<size_t>(fe[static_cast<size_t>(i)]);
      T += edge_world_tensor(h0[e], fr.t[e], fr.b[e]) / real(3);
    }
    const vec3 u = e1.normalized(), v = n.normalized().cross(u);
    Eigen::Matrix<real, 2, 2> S;
    S << u.dot(T * u), u.dot(T * v), v.dot(T * u), v.dot(T * v);
    for (const curvature_shear_target &spec : targets) {
      if (!std::isfinite(spec.shear) || spec.shear == real(0))
        continue;
      const real angle = std::isfinite(spec.angle) ? spec.angle : real(0.25) * real(M_PI);
      Eigen::Matrix<real, 2, 2> U;
      if (!face_shear_stretch(S, spec.shear, angle, U))
        continue;
      constraints.push_back(face_stretch_target::create(
          std::array<index_t, 3>{id[0], id[1], id[2]}, x, U, w, held, "face_shear"));
    }
  }
}

/// Same SVD polar as `polar_rotation` (deformation `Σ cur restᵀ`, not covariance).
inline bool diamond_rotation(const mat3 &M, mat3 &R) {
  return polar_rotation(M, R);
}

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif
