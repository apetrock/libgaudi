#ifndef GAUDI_HEPWORTH_SHELL_SHAPE_BLOCK_HPP
#define GAUDI_HEPWORTH_SHELL_SHAPE_BLOCK_HPP

#include <exception>
#include <memory>
#include <vector>

#include "gaudi/asawa/datums.hpp"
#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/hepworth/block/cr_dirichlet.hpp"
#include "gaudi/hepworth/block/sim_block.hpp"
#include "gaudi/hepworth/block/solver_composition.hpp"
#include "gaudi/hepworth/blocks/block_wrapper.hpp"
#include "gaudi/hepworth/block/shape_operator_datum.hpp"
#include "gaudi/kusama/shape_operator.hpp"
#include "gaudi/kusama/vector_dirichlet_guided.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

/// World tangent of each compact edge in `H`. `H` is the Dirichlet unknown:
/// parallel coefficients, then perpendicular. This is the vector the logger draws.
inline std::vector<vec3> extract_tan(asawa::shell::shell &shell,
                                     const std::vector<vec3> &x,
                                     const std::vector<real> &H) {
  std::vector<int> slot_map;
  const int nE = kusama::build_compact_edge_dof_map(shell, slot_map);
  if (nE <= 0 || static_cast<int>(H.size()) < 2 * nE)
    return {};
  std::vector<asawa::shell::CornerId> dof_to_corner;
  kusama::build_dof_to_corner(shell, slot_map, nE, dof_to_corner);
  std::vector<std::vector<asawa::shell::FaceId>> incident;
  kusama::build_edge_incident_faces(shell, slot_map, nE, incident);
  std::vector<vec3> tan(static_cast<size_t>(nE), vec3::Zero());
  for (int e = 0; e < nE; ++e) {
    const asawa::shell::CornerId c = dof_to_corner[static_cast<size_t>(e)];
    if (static_cast<int>(c) < 0)
      continue;
    const vec3 n = kusama::edge_average_normal(
        shell, x, incident[static_cast<size_t>(e)]);
    const vec3 t = kusama::edge_guidance_vector_3d(
        shell, c, x, n, H[static_cast<size_t>(e)],
        H[static_cast<size_t>(e + nE)]);
    if (t.allFinite())
      tan[static_cast<size_t>(e)] = t;
  }
  return tan;
}

/// Polar of the incident edge tangents. `frames[i].col(0)` is the vertex tangent.
inline std::vector<mat3> vertex_tangent_frames(asawa::shell::shell &shell,
                                        const std::vector<vec3> &x,
                                        const std::vector<vec3> &tangents) {
  const int n = shell.vert_count();
  std::vector<mat3> frames(static_cast<size_t>(std::max(n, 0)), mat3::Zero());
  if (n <= 0 || static_cast<int>(x.size()) != n || tangents.empty())
    return frames;
  std::vector<int> slot_map;
  const int nE = kusama::build_compact_edge_dof_map(shell, slot_map);
  if (nE <= 0 || static_cast<int>(tangents.size()) < nE)
    return frames;
  for (int i = 0; i < n; ++i) {
    const asawa::shell::VertId vi = asawa::shell::vert_id(i);
    if (static_cast<int>(shell.vbegin(vi)) < 0)
      continue;
    const std::vector<real> weights =
        asawa::shell::vert_cotan_weights(shell, vi, x);
    std::vector<vec3> ring;
    ring.reserve(weights.size());
    shell.const_for_each_vertex(
        vi, [&](asawa::shell::CornerId ci, const asawa::shell::shell &Ms) {
          (void)Ms;
          const int slot = static_cast<int>(ci) / 2;
          const int e =
              slot >= 0 && slot < static_cast<int>(slot_map.size())
                  ? slot_map[static_cast<size_t>(slot)]
                  : -1;
          ring.push_back(e >= 0 && e < nE ? tangents[static_cast<size_t>(e)]
                                          : vec3::Zero());
        });
    mat3 R;
    if (asawa::shell::edge_to_vert(ring, weights, R))
      frames[static_cast<size_t>(i)] = R;
  }
  return frames;
}

/// Vertex operator in the polar frame of the incident edge tangents.
/// `charts`, when set, receives that frame and `(k11, k12, k22)`.
/// `g_par` and `g_perp`, when both set, are the cached coordinate gradients.
inline bool vertex_weingarten(asawa::shell::shell &shell,
                             const std::vector<vec3> &x,
                             const std::vector<real> &coeff,
                             std::vector<vec3> &H, std::vector<vec3> &tang,
                             std::vector<vec3> &axes, std::vector<char> &wrote,
                             std::vector<vertex_chart> *charts = nullptr,
                             const std::vector<vec3> *g_par = nullptr,
                             const std::vector<vec3> *g_perp = nullptr) {
  const int n = shell.vert_count();
  if (static_cast<int>(x.size()) != n)
    return false;
  std::vector<int> slot_map;
  const int nE = kusama::build_compact_edge_dof_map(shell, slot_map);
  if (nE <= 0 || static_cast<int>(coeff.size()) < 2 * nE)
    return false;
  std::vector<asawa::shell::CornerId> dof_to_corner;
  kusama::build_dof_to_corner(shell, slot_map, nE, dof_to_corner);
  std::vector<cr_tri> tris;
  std::vector<int> face_local(static_cast<size_t>(shell.face_count()), -1);
  for (asawa::shell::FaceId fi : shell.get_face_range(true)) {
    if (shell.fsize(fi) != 3)
      continue;
    cr_tri tri;
    bool filled = false;
    shell.const_for_each_face_tri(
        fi, [&](asawa::shell::CornerId c0, asawa::shell::CornerId c1,
                asawa::shell::CornerId c2, const asawa::shell::shell &Ms) {
          if (filled)
            return;
          filled = true;
          const asawa::shell::CornerId c[3] = {c0, c1, c2};
          for (int i = 0; i < 3; ++i) {
            tri.v[static_cast<size_t>(i)] = static_cast<int>(Ms.vert(c[i]));
            const int slot = static_cast<int>(c[i]) / 2;
            tri.e[static_cast<size_t>(i)] =
                slot >= 0 && slot < static_cast<int>(slot_map.size())
                    ? slot_map[static_cast<size_t>(slot)]
                    : -1;
          }
        });
    face_local[idx(fi)] =
        static_cast<int>(tris.size());
    tris.push_back(tri);
  }
  if (tris.empty())
    return false;
  std::vector<cr_edge> edges(static_cast<size_t>(nE));
  for (int e = 0; e < nE; ++e) {
    const asawa::shell::CornerId c = dof_to_corner[static_cast<size_t>(e)];
    if (static_cast<int>(c) < 0)
      continue;
    edges[static_cast<size_t>(e)].v0 = static_cast<int>(shell.vert(c));
    edges[static_cast<size_t>(e)].v1 =
        static_cast<int>(shell.vert(shell.other(c)));
  }
  const std::vector<cr_face_operator> faceS =
      cr_face_operators(x, tris, edges, g_par, g_perp);
  const std::vector<vec3> edge_tan = extract_tan(shell, x, coeff);
  const std::vector<mat3> frames = vertex_tangent_frames(shell, x, edge_tan);
  wrote.assign(static_cast<size_t>(n), 0);
  if (charts)
    charts->assign(static_cast<size_t>(n), vertex_chart{});
  for (int i = 0; i < n; ++i) {
    const asawa::shell::VertId vi = asawa::shell::vert_id(i);
    if (static_cast<int>(shell.vbegin(vi)) < 0)
      continue;
    std::vector<int> faces;
    shell.const_for_each_vertex(
        vi, [&](asawa::shell::CornerId ci, const asawa::shell::shell &Ms) {
          const int lj = static_cast<int>(Ms.face(ci));
          if (lj >= 0 && lj < static_cast<int>(face_local.size())) {
            const int local = face_local[static_cast<size_t>(lj)];
            if (local >= 0)
              faces.push_back(local);
          }
        });
    vec3 N = asawa::shell::vert_normal(shell, vi, x);
    if (!N.allFinite() || !(N.squaredNorm() > real(1e-16)))
      continue;
    N.normalize();
    if (i >= static_cast<int>(frames.size()))
      continue;
    const mat3 &R = frames[static_cast<size_t>(i)];
    if (!R.allFinite() || !(R.col(0).squaredNorm() > real(1e-16)))
      continue;
    vec3 T = R.col(0);
    T -= N * T.dot(N);
    if (!(T.squaredNorm() > real(1e-16)) || !T.allFinite())
      continue;
    T.normalize();
    vec3 B = N.cross(T);
    if (!(B.squaredNorm() > real(1e-16)))
      continue;
    B.normalize();
    const mat3 S = cr_vertex_operator(faceS, faces);
    const vec3 h(T.dot(S * T), T.dot(S * B), B.dot(S * B));
    if (!h.allFinite())
      continue;
    real k1, k2, x1, y1;
    shape_sym2(h[0], h[1], h[2], k1, k2, x1, y1);
    vec3 axis = x1 * T + y1 * B;
    if (!(axis.squaredNorm() > real(1e-16)))
      axis = T;
    else
      axis.normalize();
    H[static_cast<size_t>(i)] = h;
    tang[static_cast<size_t>(i)] = T;
    axes[static_cast<size_t>(i)] = axis;
    wrote[static_cast<size_t>(i)] = 1;
    if (charts) {
      vertex_chart &chart = (*charts)[static_cast<size_t>(i)];
      chart.T = T;
      chart.B = B;
      chart.N = N;
      chart.h = h;
      chart.ok = true;
    }
  }
  return true;
}

/// Write Dirichlet world tangents onto the CR edge datum (remesh carry).
inline void commit_dirichlet_to_edge_datum(asawa::shell::shell &mesh,
                                           shell_cr_block &cr,
                                           const std::vector<real> &coeff) {
  if (!cr.tangents)
    return;
  std::vector<int> slot_map;
  const int nE = kusama::build_compact_edge_dof_map(mesh, slot_map);
  if (nE <= 0 || static_cast<int>(coeff.size()) < 2 * nE)
    return;
  std::vector<asawa::shell::CornerId> dof_to_corner;
  kusama::build_dof_to_corner(mesh, slot_map, nE, dof_to_corner);
  const std::vector<vec3> &x = asawa::const_get_vec_data(mesh, 0);
  const std::vector<vec3> tan = extract_tan(mesh, x, coeff);
  if (static_cast<int>(tan.size()) < nE)
    return;
  auto &T = cr.tangents->data();
  for (int e = 0; e < nE; ++e) {
    const asawa::shell::CornerId c = dof_to_corner[static_cast<size_t>(e)];
    if (static_cast<int>(c) < 0)
      continue;
    const vec3 &tangent = tan[static_cast<size_t>(e)];
    if (!tangent.allFinite() || !(tangent.squaredNorm() > real(1e-16)))
      continue;
    const int slot = static_cast<int>(c) / 2;
    if (slot < 0)
      continue;
    if (static_cast<size_t>(slot) >= T.size())
      T.resize(static_cast<size_t>(slot) + 1, vec3::Zero());
    T[static_cast<size_t>(slot)] = tangent;
  }
}

/// Square-one / fill: Dirichlet → edge datum + `u` → `vertex_frames`. No
/// parallel shape datum.
inline bool seed_cr_from_dirichlet(asawa::shell::shell &mesh, shell_cr_block &cr,
                                   bool reseed = true, real lambda = real(100),
                                   real frame_angle = 0) {
  if (!mesh.verts_are_dense_packed() || !(lambda > real(0)))
    return false;
  // Remesh changes edge count; `store_solved_tangents` needs a live topology.
  {
    solver_context ctx;
    cr.prepare(ctx);
  }
  if (cr.topology.nE <= 0)
    return false;
  const std::vector<vec3> &x = asawa::const_get_vec_data(mesh, 0);
  const int n = mesh.vert_count();
  if (static_cast<int>(x.size()) != n)
    return false;
  // Carry rebuilds the charts from the remesh-carried tangents on the live
  // topology; old charts are stale indices.
  if (!reseed && !cr.guidance_tangents &&
      static_cast<int>(cr.u.size()) >= 2 * cr.topology.nE) {
    std::vector<vec3> H(static_cast<size_t>(n), vec3::Zero());
    std::vector<vec3> tang(static_cast<size_t>(n), vec3::Zero());
    std::vector<vec3> axes(static_cast<size_t>(n), vec3::Zero());
    std::vector<char> wrote;
    std::vector<vertex_chart> charts;
    if (vertex_weingarten(mesh, x, cr.u, H, tang, axes, wrote, &charts)) {
      cr.vertex_frames = std::move(charts);
      bool open = false;
      for (int i = 0; i < n; ++i)
        if (!wrote[static_cast<size_t>(i)])
          open = true;
      if (!open)
        return true;
    }
  }
  Eigen::VectorXd field;
  try {
    field = kusama::solve_curvature_guided_vector_dirichlet(
        mesh, x, lambda, asawa::shell::face_curvature_stencil::one_ring);
  } catch (const std::exception &) {
    return false;
  }
  if (field.size() < 2)
    return false;
  std::vector<real> coeff(static_cast<size_t>(field.size()));
  for (int i = 0; i < field.size(); ++i)
    coeff[static_cast<size_t>(i)] = field(i);
  commit_dirichlet_to_edge_datum(mesh, cr, coeff);
  cr.store_solved_tangents(coeff, frame_angle, /*replace_all=*/reseed);
  std::vector<vec3> solved(static_cast<size_t>(n), vec3::Zero());
  std::vector<vec3> tang(static_cast<size_t>(n), vec3::Zero());
  std::vector<vec3> axes(static_cast<size_t>(n), vec3::Zero());
  std::vector<char> wrote(static_cast<size_t>(n), 0);
  std::vector<vertex_chart> charts;
  if (!vertex_weingarten(mesh, x, coeff, solved, tang, axes, wrote, &charts))
    return false;
  if (reseed || static_cast<int>(cr.vertex_frames.size()) != n)
    cr.vertex_frames.assign(static_cast<size_t>(n), vertex_chart{});
  for (int i = 0; i < n; ++i) {
    if (!wrote[static_cast<size_t>(i)])
      continue;
    if (!reseed && cr.vertex_frames[static_cast<size_t>(i)].ok)
      continue;
    cr.vertex_frames[static_cast<size_t>(i)] = charts[static_cast<size_t>(i)];
  }
  return true;
}

/// Vertex shape operator `(k11, k12, k22)` in that vertex's tangent frame.
/// Mass is vertex area, the same diagonal the position block contributes.
class shell_shape_block : public dof_block_wrapper {
public:
  using ptr = std::shared_ptr<shell_shape_block>;

  explicit shell_shape_block(asawa::shell::shell::ptr mesh)
      : mesh(std::move(mesh)) {
    const int n = this->mesh ? this->mesh->vert_count() : 0;
    _op = shape_operator_datum::create(
        asawa::VERTEX, std::vector<vec3>(static_cast<size_t>(n), vec3::Zero()));
    _shape_id = this->mesh->insert_datum(_op);
    _vel_id = asawa::init_vert_datum<vec3>(*this->mesh, vec3::Zero());
  }

  asawa::shell::shell::ptr mesh;
  std::vector<vec3> mass;
  std::vector<vec3> forces;

  std::vector<vec3> &shape() { return _op->data(); }

  /// Pin target. Reseed rewrites this with the new bake; otherwise it is
  /// written once when a vertex is first seeded.
  std::vector<vec3> &rest() { return _op->rest; }

  /// Chart tangent for `shape()` / `rest`. Same basis the polar cache uses.
  std::vector<vec3> &chart_T() { return _op->chart_T; }

  std::vector<vec3> &velocity() { return asawa::get_vec_data(*mesh, _vel_id); }

  std::vector<real> &seeded() { return _op->mark; }

  shape_operator_datum::ptr descriptor() const { return _op; }

  /// Copy polar-frame tangents onto the datum so remesh can transport H in
  /// the same charts the projection constraints use. When the chart turns,
  /// `shape()` is pulled into the new basis — never leave coeffs under a
  /// replaced `chart_T`.
  void sync_charts_from(const std::vector<vertex_chart> &charts) {
    if (!_op || !mesh)
      return;
    const int n = mesh->vert_count();
    _op->fit_sides(n);
    std::vector<vec3> &T = chart_T();
    std::vector<vec3> &H = shape();
    if (static_cast<int>(H.size()) < n)
      H.resize(static_cast<size_t>(n), vec3::Zero());
    for (int i = 0; i < n; ++i) {
      if (static_cast<size_t>(i) >= charts.size())
        break;
      const vertex_chart &c = charts[static_cast<size_t>(i)];
      if (!c.ok || !c.T.allFinite() || !(c.T.squaredNorm() > real(1e-16)) ||
          !c.N.allFinite() || !(c.N.squaredNorm() > real(1e-16)))
        continue;
      const vec3 &Told = T[static_cast<size_t>(i)];
      const vec3 &Hold = H[static_cast<size_t>(i)];
      if (Told.allFinite() && Told.squaredNorm() > real(1e-16) &&
          Hold.allFinite() &&
          std::abs(Told.normalized().dot(c.T.normalized())) < real(0.999)) {
        vec3 T0, B0, T1, B1;
        if (surface_chart(Told, c.N, T0, B0) &&
            surface_chart(c.T, c.N, T1, B1)) {
          if (T1.dot(T0) < real(0)) {
            T1 = -T1;
            B1 = -B1;
          }
          const mat3 S = shape_world_from_coeffs(Hold, T0, B0);
          if (S.allFinite()) {
            const vec3 h1 = shape_coeffs_from_world(S, T1, B1);
            if (h1.allFinite())
              H[static_cast<size_t>(i)] = h1;
          }
        }
      }
      T[static_cast<size_t>(i)] = c.T;
    }
  }

  /// Dirichlet bake onto edge tangents, then vertex frames / shape operators
  /// from that field (same `vertex_weingarten` the projection cache uses).
  /// `reseed`: square one — replace carried `u`, rewrite every framed vertex's
  /// `shape()`, `rest`, `chart_T`, and `cr.vertex_frames`. Off fills only
  /// unmarked zeros and leaves an existing rest.
  void seed_from_operator(shell_cr_block &cr, bool reseed = false,
                          real lambda = real(100), real frame_angle = 0) {
    if (!mesh || !_op)
      return;
    if (!mesh->verts_are_dense_packed())
      return;
    if (!(lambda > real(0)))
      return;
    const std::vector<vec3> &x = asawa::const_get_vec_data(*mesh, 0);
    std::vector<vec3> &H = shape();
    std::vector<vec3> &R = rest();
    std::vector<real> &mark = seeded();
    const int n = mesh->vert_count();
    _op->fit_sides(n);
    if (static_cast<int>(H.size()) < n)
      H.resize(static_cast<size_t>(n), vec3::Zero());
    if (x.size() != static_cast<size_t>(n))
      return;
    if (!reseed) {
      bool open = false;
      for (int i = 0; i < n; ++i) {
        if (mark[static_cast<size_t>(i)] > real(0.5))
          continue;
        if (H[static_cast<size_t>(i)].squaredNorm() > real(1e-24)) {
          mark[static_cast<size_t>(i)] = real(1);
          continue;
        }
        open = true;
      }
      if (!open)
        return;
    }
    Eigen::VectorXd field;
    try {
      field = kusama::solve_curvature_guided_vector_dirichlet(
          *mesh, x, lambda, asawa::shell::face_curvature_stencil::one_ring);
    } catch (const std::exception &) {
      return;
    }
    if (field.size() < 2)
      return;
    std::vector<real> coeff(static_cast<size_t>(field.size()));
    for (int i = 0; i < field.size(); ++i)
      coeff[static_cast<size_t>(i)] = field(i);
    store_dirichlet_edge_tangents(cr, coeff);
    // Reseed must replace carried remesh `u` so frames and H share one field.
    cr.store_solved_tangents(coeff, frame_angle, /*replace_all=*/reseed);
    std::vector<vec3> solved(static_cast<size_t>(n), vec3::Zero());
    std::vector<vec3> tang(static_cast<size_t>(n), vec3::Zero());
    std::vector<vec3> axes(static_cast<size_t>(n), vec3::Zero());
    std::vector<char> wrote(static_cast<size_t>(n), 0);
    std::vector<vertex_chart> charts;
    if (!vertex_weingarten(*mesh, x, coeff, solved, tang, axes, wrote, &charts))
      return;
    std::vector<vec3> &Tc = chart_T();
    std::vector<vec3> &vel = velocity();
    if (reseed) {
      cr.vertex_frames = std::move(charts);
      if (static_cast<int>(vel.size()) < n)
        vel.assign(static_cast<size_t>(n), vec3::Zero());
    }
    for (int i = 0; i < n; ++i) {
      if (!wrote[static_cast<size_t>(i)])
        continue;
      if (!reseed) {
        if (mark[static_cast<size_t>(i)] > real(0.5))
          continue;
        if (H[static_cast<size_t>(i)].squaredNorm() > real(1e-24))
          continue;
      }
      H[static_cast<size_t>(i)] = solved[static_cast<size_t>(i)];
      if (tang[static_cast<size_t>(i)].allFinite() &&
          tang[static_cast<size_t>(i)].squaredNorm() > real(1e-16))
        Tc[static_cast<size_t>(i)] = tang[static_cast<size_t>(i)];
      if (reseed) {
        // Pin target must match the new bake; leaving an old rest is what
        // made reseed look noisier than carry.
        R[static_cast<size_t>(i)] = solved[static_cast<size_t>(i)];
        mark[static_cast<size_t>(i)] = real(1);
        vel[static_cast<size_t>(i)] = vec3::Zero();
        continue;
      }
      if (mark[static_cast<size_t>(i)] > real(0.5))
        continue;
      R[static_cast<size_t>(i)] = solved[static_cast<size_t>(i)];
      mark[static_cast<size_t>(i)] = real(1);
    }
  }

  void prepare(solver_context &ctx) override {
    (void)ctx;
    if (!mesh)
      return;
    const int n = mesh->vert_count();
    auto fit = [n](std::vector<vec3> &v) {
      if (static_cast<int>(v.size()) < n)
        v.resize(static_cast<size_t>(n), vec3::Zero());
      else if (static_cast<int>(v.size()) > n)
        v.resize(static_cast<size_t>(n));
    };
    fit(shape());
    fit(rest());
    fit(velocity());
    std::vector<real> &mark = seeded();
    if (static_cast<int>(mark.size()) < n)
      mark.resize(static_cast<size_t>(n), real(0));
    else if (static_cast<int>(mark.size()) > n)
      mark.resize(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
      vec3 &h = shape()[static_cast<size_t>(i)];
      vec3 &v = velocity()[static_cast<size_t>(i)];
      if (h.allFinite() && v.allFinite())
        continue;
      h = vec3::Zero();
      v = vec3::Zero();
      rest()[static_cast<size_t>(i)] = vec3::Zero();
      mark[static_cast<size_t>(i)] = real(0);
    }
    const std::vector<vec3> &x = asawa::const_get_vec_data(*mesh, 0);
    mass = asawa::shell::vertex_areas_3(*mesh, x);
    if (static_cast<int>(mass.size()) < n)
      mass.resize(static_cast<size_t>(n), vec3::Zero());
    forces.assign(static_cast<size_t>(n), vec3::Zero());
  }

  void apply_external_input(const std::vector<vec3> &input, real scale = 1.0) override {
    (void)input;
    (void)scale;
  }

  void flush_external_inputs() override {
    forces.assign(forces.size(), vec3::Zero());
  }

  void emit_dof_blocks(std::vector<sim_block::ptr> &out) override {
    _vec3 = vec3_block::create(mass, shape(), velocity(), forces);
    out.push_back(_vec3);
  }

  vec3_block::ptr dof_vec3() const { return _vec3; }

private:
  void store_dirichlet_edge_tangents(shell_cr_block &cr,
                                     const std::vector<real> &coeff) {
    if (!mesh || !cr.tangents)
      return;
    std::vector<int> slot_map;
    const int nE = kusama::build_compact_edge_dof_map(*mesh, slot_map);
    if (nE <= 0 || static_cast<int>(coeff.size()) < 2 * nE)
      return;
    std::vector<asawa::shell::CornerId> dof_to_corner;
    kusama::build_dof_to_corner(*mesh, slot_map, nE, dof_to_corner);
    const std::vector<vec3> &x = asawa::const_get_vec_data(*mesh, 0);
    const std::vector<vec3> tan = extract_tan(*mesh, x, coeff);
    if (static_cast<int>(tan.size()) < nE)
      return;
    auto &T = cr.tangents->data();
    for (int e = 0; e < nE; ++e) {
      const asawa::shell::CornerId c = dof_to_corner[static_cast<size_t>(e)];
      if (static_cast<int>(c) < 0)
        continue;
      const vec3 &tangent = tan[static_cast<size_t>(e)];
      if (!tangent.allFinite() || !(tangent.squaredNorm() > real(1e-16)))
        continue;
      const int slot = static_cast<int>(c) / 2;
      if (slot < 0)
        continue;
      if (static_cast<size_t>(slot) >= T.size())
        T.resize(static_cast<size_t>(slot) + 1, vec3::Zero());
      T[static_cast<size_t>(slot)] = tangent;
    }
  }

  index_t _shape_id = -1;
  index_t _vel_id = -1;
  shape_operator_datum::ptr _op;
  vec3_block::ptr _vec3;
};

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif
