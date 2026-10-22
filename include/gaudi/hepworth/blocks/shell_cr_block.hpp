#ifndef GAUDI_HEPWORTH_SHELL_CR_BLOCK_HPP
#define GAUDI_HEPWORTH_SHELL_CR_BLOCK_HPP

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

#include "gaudi/asawa/datums.hpp"
#include "gaudi/asawa/shell/stencil_walker.hpp"
#include "gaudi/hepworth/block/cr_field.hpp"
#include "gaudi/hepworth/block/cr_shape.hpp"
#include "gaudi/hepworth/block/sim_block.hpp"
#include "gaudi/hepworth/block/solver_composition.hpp"
#include "gaudi/hepworth/blocks/block_wrapper.hpp"
#include "gaudi/kusama/stencil_table.hpp"
#include "gaudi/kusama/vector_dirichlet_guided.hpp"
#include "gaudi/kusama/weingarten.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

/// Turn one edge `(gp, gq)` in its plane. Positive `theta` turns parallel
/// toward perpendicular. The pair length stays.
inline void rotate_edge_pair(real &gp, real &gq, real theta) {
  if (!std::isfinite(theta) || !(std::abs(theta) > real(0)) ||
      !std::isfinite(gp) || !std::isfinite(gq))
    return;
  const real c = std::cos(theta);
  const real s = std::sin(theta);
  const real p = gp;
  const real q = gq;
  gp = p * c - q * s;
  gq = p * s + q * c;
}

/// Turn every edge pair. `u` is parallel, then perpendicular, `nE` of each.
inline void rotate_edge_coefficients(std::vector<real> &u, int nE, real theta) {
  if (nE <= 0 || static_cast<int>(u.size()) < 2 * nE)
    return;
  for (int e = 0; e < nE; ++e)
    rotate_edge_pair(u[static_cast<size_t>(e)],
                     u[static_cast<size_t>(e + nE)], theta);
}

/// Two Stein coefficients per edge. Index `e` is parallel, `e + nE` perpendicular.
class cr_coeff_block : public sim_block {
public:
  using ptr = std::shared_ptr<cr_coeff_block>;

  static ptr create(std::vector<real> &u, std::vector<real> &v,
                    std::vector<real> &mass, int nE) {
    return std::make_shared<cr_coeff_block>(u, v, mass, nE);
  }

  cr_coeff_block(std::vector<real> &u, std::vector<real> &v,
                 std::vector<real> &mass, int nE)
      : u(u), v(v), mass(mass), nE(nE) {}

  void map_to_x(vecX &q) override {
    _offset = q.size();
    q = concat(q, to(u));
  }

  void map_mass(vecX &q) override { q = concat(q, to(mass)); }

  void map_from_x(const vecX &q, const real &h, const real &damp) override {
    if (u.empty() || !(h > real(0)))
      return;
    std::vector<real> u1(u.size(), real(0));
    for (size_t i = 0; i < u.size(); ++i)
      u1[i] = q(_offset + static_cast<index_t>(i));
    const std::vector<real> &u0 = u_prev.size() == u.size() ? u_prev : u;
    for (size_t i = 0; i < u.size(); ++i)
      v[i] = (real(1) - damp) / h * (u1[i] - u0[i]);
    u = std::move(u1);
  }

  void integrate_inertia(const real &h) override {
    u_prev = u;
    const size_t n = std::min(u.size(), v.size());
    for (size_t i = 0; i < n; ++i)
      u[i] += h * v[i];
  }

  index_t get_offset_idx(index_t ii) const override { return _offset + ii; }

  std::vector<real> &u;
  std::vector<real> &v;
  std::vector<real> &mass;
  int nE = 0;
  std::vector<real> u_prev;
};

/// Edge tangents for the shell projection. Owned here so the field survives remeshing
/// until the edge count changes, then it is seeded from the curvature guidance.
class shell_cr_block : public dof_block_wrapper {
public:
  using ptr = std::shared_ptr<shell_cr_block>;

  explicit shell_cr_block(asawa::shell::shell::ptr mesh)
      : mesh(std::move(mesh)) {
    if (!this->mesh)
      return;
    std::vector<vec3> carried(this->mesh->edge_count(), vec3::Zero());
    tangents = edge_tangent_datum::create(asawa::prim_type::EDGE, carried);
    this->mesh->insert_datum(tangents);
  }

  asawa::shell::shell::ptr mesh;
  cr_topology topology;
  edge_tangent_datum::ptr tangents;
  std::vector<real> u;
  std::vector<real> v;
  std::vector<real> mass;
  /// Live stencil table of the triangles; edge ids are the CR edge dofs.
  kusama::stencil_table table;
  /// Incident edges per vertex: vertex charts read `h` from the edge store.
  kusama::weingarten::vertex_stars stars;
  /// The curvature store: `(k11, k12, k22)` per edge in its unfolded frame
  /// `(t_e, N x t_e)`, `t_e` from the lower vertex id to the higher. Set to
  /// the live edge-average reading each prepare, so the Weingarten coupling
  /// is idle at the frame start. Massless: rows determine it. Not carried
  /// across a remesh.
  std::vector<vec3> shape;
  /// Edge coordinate gradients. A cache: not a solver unknown.
  std::vector<vec3> g_par;
  std::vector<vec3> g_perp;
  /// Vertex polar frame (edge-tangent polar) and `h`, the local polar of the
  /// edge store at the step start. A cache: not a solver unknown.
  std::vector<vertex_chart> vertex_frames;
  /// `(k11, k12, k22)` per face, in the face's first edge, for the face-site
  /// operators (`one_ring_operators` off). A solver unknown.
  std::vector<vec3> face_shape;
  std::vector<vec3> shape_velocity;
  std::vector<vec3> shape_force;
  std::vector<vec3> shape_mass;
  std::vector<vec3> face_velocity;
  std::vector<vec3> face_force;
  std::vector<vec3> face_mass;
  /// Applied once, when a Dirichlet solution fills an edge.
  real frame_angle = 0;
  /// `prepare` loaded every edge from raw guidance. The next Dirichlet
  /// solution replaces that placeholder.
  bool guidance_tangents = false;

  /// Write Dirichlet coefficients into `u` and rebuild edge operators.
  /// Empty slots and guidance placeholders are always filled. A carried pair
  /// is left alone unless `replace_all` (reseed) or the prepare was guidance.
  /// `angle` turns each newly written pair once. Velocity is dropped.
  void store_solved_tangents(const std::vector<real> &coeff, real angle,
                             bool replace_all = false) {
    const int nE = topology.nE;
    if (nE <= 0 || static_cast<int>(coeff.size()) < 2 * nE)
      return;
    if (static_cast<int>(u.size()) < 2 * nE)
      u.assign(static_cast<size_t>(2 * nE), real(0));
    const bool replace = replace_all || guidance_tangents;
    guidance_tangents = false;
    for (int e = 0; e < nE; ++e) {
      real &gp = u[static_cast<size_t>(e)];
      real &gq = u[static_cast<size_t>(e + nE)];
      if (!replace && gp * gp + gq * gq > real(1e-16))
        continue;
      gp = coeff[static_cast<size_t>(e)];
      gq = coeff[static_cast<size_t>(e + nE)];
      rotate_edge_pair(gp, gq, angle);
    }
    v.assign(u.size(), real(0));
  }

  /// Write the solved coefficients back onto the edge slots before a remesh.
  void commit_edge_tangents() {
    if (!mesh || !tangents || topology.nE <= 0)
      return;
    if (static_cast<int>(u.size()) < 2 * topology.nE)
      return;
    const std::vector<vec3> &x = asawa::const_get_vec_data(*mesh, 0);
    auto &T = tangents->data();
    for (int e = 0; e < topology.nE; ++e) {
      const int c = static_cast<int>(topology.dof_to_corner[static_cast<size_t>(e)]);
      const int slot = c / 2;
      if (slot < 0)
        continue;
      vec3 e_par, e_perp;
      if (!cr_edge_frame(*mesh, x, topology.dof_to_corner[static_cast<size_t>(e)],
                         topology.edge_incident[static_cast<size_t>(e)], e_par,
                         e_perp))
        continue;
      const vec3 dir =
          u[static_cast<size_t>(e)] * e_par +
          u[static_cast<size_t>(e + topology.nE)] * e_perp;
      if (static_cast<size_t>(slot) >= T.size())
        T.resize(static_cast<size_t>(slot) + 1, vec3::Zero());
      T[static_cast<size_t>(slot)] = dir;
    }
  }

  int count() const { return topology.nE; }

  void prepare(solver_context &ctx) override {
    (void)ctx;
    cr_topology next;
    guidance_tangents = false;
    if (!mesh || !assemble_cr_topology(*mesh, next)) {
      topology = cr_topology{};
      table = kusama::stencil_table{};
      stars = kusama::weingarten::vertex_stars{};
      u.clear();
      v.clear();
      mass.clear();
      shape.clear();
      shape_velocity.clear();
      shape_force.clear();
      shape_mass.clear();
      g_par.clear();
      g_perp.clear();
      face_shape.clear();
      face_velocity.clear();
      face_force.clear();
      face_mass.clear();
      return;
    }
    const std::vector<vec3> &x = asawa::const_get_vec_data(*mesh, 0);
    topology = std::move(next);
    u.assign(static_cast<size_t>(2 * topology.nE), real(0));
    v.assign(u.size(), real(0));
    bool carried = false;
    if (tangents) {
      const auto &T = tangents->data();
      for (int e = 0; e < topology.nE; ++e) {
        const int c =
            static_cast<int>(topology.dof_to_corner[static_cast<size_t>(e)]);
        const int slot = c / 2;
        if (slot < 0 || static_cast<size_t>(slot) >= T.size())
          continue;
        const vec3 dir = T[static_cast<size_t>(slot)];
        if (!(dir.squaredNorm() > real(1e-16)))
          continue;
        vec3 e_par, e_perp;
        if (!cr_edge_frame(*mesh, x, topology.dof_to_corner[static_cast<size_t>(e)],
                           topology.edge_incident[static_cast<size_t>(e)], e_par,
                           e_perp))
          continue;
        u[static_cast<size_t>(e)] = dir.dot(e_par);
        u[static_cast<size_t>(e + topology.nE)] = dir.dot(e_perp);
        carried = true;
      }
    }
    // Placeholder until the Dirichlet presolve. A carried datum is the field,
    // and a split spoke stays empty so that presolve can fill it.
    if (!carried && x.size() == static_cast<size_t>(mesh->vert_count())) {
      const Eigen::VectorXd g = cr_curvature_guidance(*mesh, x, topology);
      for (int i = 0; i < g.size() && i < static_cast<int>(u.size()); ++i)
        u[static_cast<size_t>(i)] = g(i);
      guidance_tangents = true;
    }
    const Eigen::VectorXd md = kusama::build_edge_barycentric_mass_diagonal(
        *mesh, x, topology.nE, topology.dof_to_corner);
    mass.resize(static_cast<size_t>(2 * topology.nE));
    for (int i = 0; i < md.size() && i < static_cast<int>(mass.size()); ++i)
      mass[static_cast<size_t>(i)] = std::max(md(i), real(1e-8));
    const int nE = topology.nE;
    table = kusama::build_stencil_table(asawa::shell::shell_walker(*mesh));
    const bool live = table.n_edges == nE && x.size() == static_cast<size_t>(table.n_verts);
    shape = live ? kusama::weingarten::edge_readings(table, x)
                 : std::vector<vec3>(static_cast<size_t>(nE), vec3::Zero());
    stars = live ? kusama::weingarten::build_vertex_stars(table, x)
                 : kusama::weingarten::vertex_stars{};
    for (vec3 &s : shape)
      if (!s.allFinite())
        s = vec3::Zero();
    // A degenerate face's read-out divides by a vanishing area: its edges
    // start at zero curvature (the coupling leaves them out this frame).
    if (live) {
      const kusama::weingarten::degenerate_set bad = kusama::weingarten::find_degenerate(table, x);
      for (int e = 0; e < nE && e < static_cast<int>(bad.edge.size()); ++e)
        if (bad.edge[static_cast<size_t>(e)])
          shape[static_cast<size_t>(e)] = vec3::Zero();
    }
    shape_velocity.assign(static_cast<size_t>(nE), vec3::Zero());
    shape_force.assign(static_cast<size_t>(nE), vec3::Zero());
    // No inertia: the store is the weighted least-squares blend of its rows
    // (the read-out `weingarten_gn` and the target rows), which keep it
    // determined. The floor only keeps the diagonal positive.
    shape_mass.assign(static_cast<size_t>(nE), vec3::Constant(real(1e-8)));
    g_par.assign(static_cast<size_t>(nE), vec3::Zero());
    g_perp.assign(static_cast<size_t>(nE), vec3::Zero());
    const int nF = mesh->face_count();
    face_shape.assign(static_cast<size_t>(std::max(nF, 0)), vec3::Zero());
    face_velocity.assign(face_shape.size(), vec3::Zero());
    face_force.assign(face_shape.size(), vec3::Zero());
    face_mass.assign(face_shape.size(), vec3(real(1e-8), real(1e-8), real(1e-8)));
  }

  void apply_external_input(const std::vector<vec3> &input, real scale = 1.0) override {
    (void)input;
    (void)scale;
  }

  void flush_external_inputs() override {}

  void emit_dof_blocks(std::vector<sim_block::ptr> &out) override {
    _dof = cr_coeff_block::create(u, v, mass, topology.nE);
    _shape_dof = vec3_block::create(shape_mass, shape, shape_velocity, shape_force);
    out.push_back(_shape_dof);
    _face_dof = nullptr;
    if (face_site) {
      _face_dof =
          vec3_block::create(face_mass, face_shape, face_velocity, face_force);
      out.push_back(_face_dof);
    }
  }

  /// Emit the per-face store (`face_shape`) for the face-site rows
  /// (`face_stretch`, `face_fair`, the Dirichlet `w_pin`). Off drops a third
  /// of the unknowns when only the edge store is used.
  bool face_site = true;

  cr_coeff_block::ptr dof() const { return _dof; }
  vec3_block::ptr shape_dof() const { return _shape_dof; }
  vec3_block::ptr face_dof() const { return _face_dof; }

private:
  cr_coeff_block::ptr _dof;
  vec3_block::ptr _shape_dof;
  vec3_block::ptr _face_dof;
};

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif
