#ifndef GAUDI_HEPWORTH_BLOCKS_SHELL_E3_BLOCK_HPP
#define GAUDI_HEPWORTH_BLOCKS_SHELL_E3_BLOCK_HPP

/// The contact net's vertex quaternion frames (scripts/sympy/cr_poisson/
/// MATH.md 19.21, E3) as a solver block wrapper: one quaternion `Y` per
/// vertex (`n = Y e3 conj(Y) / |Y|^2`), and the frame-start
/// charts the position rows read (`T, B, N` from `Y`, `h` the vertex tensor
/// fitted from the frames). Reseeded every frame from the live mesh: normals
/// from the area normals, tangents from a carried vertex datum (the frames'
/// turn about the normal, kept across frames and remeshes).

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

#include "gaudi/asawa/datums.hpp"
#include "gaudi/asawa/shell/stencil_walker.hpp"
#include "gaudi/hepworth/block/cr_shape.hpp"
#include "gaudi/hepworth/block/frame_rows.hpp"
#include "gaudi/hepworth/block/quat_util.hpp"
#include "gaudi/hepworth/block/sim_block.hpp"
#include "gaudi/hepworth/block/solver_composition.hpp"
#include "gaudi/hepworth/blocks/block_wrapper.hpp"
#include "gaudi/kusama/contact_net.hpp"
#include "gaudi/kusama/stencil_table.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

class shell_e3_block : public dof_block_wrapper {
public:
  using ptr = std::shared_ptr<shell_e3_block>;

  explicit shell_e3_block(asawa::shell::shell::ptr mesh) : mesh(std::move(mesh)) {
    if (this->mesh)
      _gauge_slot = asawa::init_vert_datum<vec3>(*this->mesh, vec3::Zero());
  }

  asawa::shell::shell::ptr mesh;
  /// Live stencil table of the triangles.
  kusama::stencil_table table;
  /// Area-weighted unit vertex normals at the frame start.
  std::vector<vec3> normal;
  /// Per-vertex quaternion frames, a solver unknown.
  std::vector<vec4> frame;
  std::vector<vec4> frame_mass;
  /// Frame mass per unit vertex dual area (0: the floor only).
  real frame_mass_scale = 0;
  /// Frame-start charts: `T, B, N` from the frames, `h` the fitted tensor.
  std::vector<vertex_chart> vertex_frames;

  vec4_block::ptr frame_dof() const { return _frame_dof; }

  void prepare(solver_context &ctx) override {
    (void)ctx;
    table = kusama::stencil_table{};
    if (!mesh)
      return;
    const std::vector<vec3> &x = asawa::const_get_vec_data(*mesh, 0);
    table = kusama::build_stencil_table(asawa::shell::shell_walker(*mesh));
    const size_t n = x.size();
    const bool live = static_cast<size_t>(table.n_verts) == n;
    normal.assign(n, vec3::Zero());
    std::vector<real> area(n, real(0));
    if (live)
      for (int f = 0; f < table.n_faces; ++f) {
        const auto &id = table.ids[static_cast<size_t>(f)];
        const vec3 c = (x[static_cast<size_t>(id[1])] - x[static_cast<size_t>(id[0])])
                           .cross(x[static_cast<size_t>(id[2])] - x[static_cast<size_t>(id[0])]);
        for (int k = 0; k < 3; ++k) {
          normal[static_cast<size_t>(id[static_cast<size_t>(k)])] += c;
          area[static_cast<size_t>(id[static_cast<size_t>(k)])] += c.norm() / real(6);
        }
      }
    for (vec3 &nv : normal)
      if (nv.norm() > real(0))
        nv.normalize();
    const std::vector<vec3> *gauge =
        _gauge_slot >= 0 ? &asawa::const_get_vec_data(*mesh, _gauge_slot) : nullptr;
    frame.assign(n, vec4(1, 0, 0, 0));
    frame_mass.resize(n);
    for (size_t i = 0; i < n; ++i) {
      const vec3 nv = normal[i].norm() > real(0) ? normal[i] : vec3::UnitZ();
      vec3 t = gauge && i < gauge->size() ? (*gauge)[i] : vec3::Zero();
      t -= nv * nv.dot(t);
      if (!(t.norm() > real(1e-9))) {
        t = std::abs(nv[0]) < real(0.9) ? vec3::UnitX() : vec3::UnitY();
        t -= nv * nv.dot(t);
      }
      t.normalize();
      mat3 R;
      R.col(0) = t;
      R.col(1) = nv.cross(t);
      R.col(2) = nv;
      frame[i] = quat_ops::quat_from_matrix(R);
      frame_mass[i] = vec4::Constant(std::max(frame_mass_scale * area[i], real(1e-8)));
    }
    update_charts(x);
  }

  /// Charts from the current frames and positions.
  void update_charts(const std::vector<vec3> &x) {
    const size_t n = x.size();
    vertex_frames.assign(n, vertex_chart{});
    if (static_cast<size_t>(table.n_verts) != n || frame.size() != n)
      return;
    const kusama::contact_net::topology tp = kusama::contact_net::build_topology(table);
    std::vector<frame_detail::frame3> F(n);
    for (size_t i = 0; i < n; ++i)
      F[i] = frame_detail::frame_of(frame[i]);
    for (size_t i = 0; i < n; ++i) {
      vertex_chart &c = vertex_frames[i];
      c.T = F[i].T;
      c.B = F[i].B;
      c.N = F[i].N;
      c.h = frame_detail::fit_vertex_tensor(tp, static_cast<int>(i), x, F);
      c.ok = c.h.allFinite() && tp.spoke_begin[i + 1] - tp.spoke_begin[i] >= 3;
    }
  }

  /// After a solve: the solved frames' tangents become the carried gauge.
  void commit_gauge() {
    if (!mesh || _gauge_slot < 0)
      return;
    std::vector<vec3> &g = asawa::get_vec_data(*mesh, _gauge_slot);
    if (g.size() != frame.size())
      return;
    for (size_t i = 0; i < frame.size(); ++i)
      g[i] = quat_ops::qrot(frame[i], vec3::UnitX());
  }

  void apply_external_input(const std::vector<vec3> &input, real scale = 1.0) override {
    (void)input;
    (void)scale;
  }

  void flush_external_inputs() override {}

  void emit_dof_blocks(std::vector<sim_block::ptr> &out) override {
    _frame_dof = vec4_block::create(frame_mass, frame);
    out.push_back(_frame_dof);
  }

private:
  index_t _gauge_slot = -1;
  vec4_block::ptr _frame_dof;
};

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif
