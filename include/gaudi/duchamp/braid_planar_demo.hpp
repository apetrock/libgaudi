#ifndef GAUDI_DUCHAMP_BRAID_PLANAR_DEMO_HPP
#define GAUDI_DUCHAMP_BRAID_PLANAR_DEMO_HPP

#include <algorithm>
#include <memory>
#include <optional>
#include <string>

#include "gaudi/asawa/rod/dynamic.hpp"
#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/common.h"
#include "gaudi/duchamp/braid_planar_rod.hpp"
#include "gaudi/duchamp/braid_sphere_rod.hpp"
#include "gaudi/duchamp/demo_snapshots.hpp"
#include "gaudi/duchamp/demo_trait.hpp"
#include "gaudi/duchamp/modules/rod_forces.hpp"
#include "gaudi/hepworth/blocks/rod_position_block.hpp"
#include "gaudi/hepworth/constraints/bundles.hpp"
#include "gaudi/hepworth/nodes/block_solver_node.hpp"
#include "gaudi/hepworth/nodes/solver_builder.hpp"
#include "gaudi/windychien/braid.hpp"
#include "gaudi/windychien/catalog.hpp"
#include "liblombardi/graph_context.hpp"

namespace gaudi {
namespace duchamp {

struct braid_planar_demo_config {
  /// If true: plain weave (even frames 0-1,2-3,…; odd 1-2,3-4,…).
  /// If false: load `knot_name` from the catalog.
  bool use_plain_weave = true;
  int weave_strands = 50;
  int weave_frames = 50;
  // Knot Atlas Hoste–Thistlethwaite K11a359 (not Stoimenov 11_359).
  std::string knot_name = "K11a359";
  braid_planar_params planar{};
  // Dense circumferential subdiv (higher than weave_frames for smooth tubes).
  braid_sphere_params sphere{.n_lat = 256};
  vec3 rod_color = vec3(0.95, 0.55, 0.15);
  vec3 planar_color = vec3(0.35, 0.75, 1.0);
  real rod_radius = 0.01;
  /// Place planar chart beside the sphere for comparison.
  bool show_planar_graph = true;

  /// Position-only solver: rest-length + fairing + collisions (+ TP force).
  /// No Cosserat quat DOFs; tubes use parallel-transport frames.
  real edge_stretch_w = 5.0e-2;
  real smooth_w = 1.0e-3;
  real min_kink_w = 1.0e-4;
  real collision_w = 1.0;
  real dt = 0.05;
  real damping = 0.1;
  int iterations = 20;
  tangent_point_force_config tangent{.w = 1.0e-8, .l0 = 3.0, .p = 6.0};
};

inline const braid_planar_demo_config k_braid_planar_demo_default{};

/// Closed braid → sphere (+ optional planar graph debug overlay).
/// 1-block position solver: edge_stretch + smooth + min_kink + collisions,
/// driven by tangent-point forces. Tube frames via `_update_frames`.
class braid_planar_demo {
public:
  using ptr = std::shared_ptr<braid_planar_demo>;

  static ptr create(const braid_planar_demo_config &cfg =
                        k_braid_planar_demo_default) {
    return std::make_shared<braid_planar_demo>(cfg);
  }

  explicit braid_planar_demo(
      const braid_planar_demo_config &cfg = k_braid_planar_demo_default)
      : _cfg(cfg) {
    const windychien::braid b =
        _cfg.use_plain_weave
            ? windychien::make_plain_weave(_cfg.weave_strands,
                                           _cfg.weave_frames)
            : windychien::load_braid(_cfg.knot_name);

    braid_planar_params planar = _cfg.planar;
    planar.close = true;
    planar.center = false;
    // Mesh tube radius = rod._r, so diameter = 2*rod_radius.
    // eps_z = radius + 0.1*r ⇒ small air gap at crossings (flush was eps_z = r).
    planar.eps_z = _cfg.rod_radius + real(0.1) * _cfg.rod_radius;
    // Braid/lanyard chart units (subdivide here, then project).
    if (planar.dx <= 0.0) {
      planar.dx = 1.0;
    }
    if (planar.dy <= 0.0) {
      planar.dy = planar.dx;
    }

    // Open chart for graph viz (no wrap-around chords through the diagram).
    braid_planar_params planar_viz = planar;
    planar_viz.close = false;
    __R_planar = braid_to_planar_rod(b, planar_viz);

    __R = braid_to_sphere_rod(b, planar, _cfg.sphere);
    __R->_r = _cfg.rod_radius;
    __R_planar->_r = _cfg.rod_radius;

    if (_cfg.show_planar_graph) {
      // Center in x–y; keep over/under as z; park below the sphere.
      real x_min = __R_planar->x()[0][0], x_max = x_min;
      real y_min = __R_planar->x()[0][1], y_max = y_min;
      for (const vec3 &q : __R_planar->x()) {
        x_min = std::min(x_min, q[0]);
        x_max = std::max(x_max, q[0]);
        y_min = std::min(y_min, q[1]);
        y_max = std::max(y_max, q[1]);
      }
      const vec3 cen(0.5 * (x_min + x_max), 0.5 * (y_min + y_max), 0.0);
      const vec3 shift(0.0, 0.0, -_cfg.sphere.radius - 0.75);
      for (vec3 &q : __R_planar->x()) {
        q[0] -= cen[0];
        q[1] -= cen[1];
        q += shift;
      }
      __R_planar->_init_params();
    }

    const real lavg = std::max(__R->lavg(), real(1e-6));
    __Rd = asawa::rod::dynamic::create(__R, 0.25 * lavg, 2.5 * lavg,
                                       0.25 * lavg);

    _rod_pos = std::make_shared<hepworth::block::rod_position_block>(__R, __Rd);

    _solver_config =
        hepworth::block::block_solver_builder<
            hepworth::block::rod_position_block>::create()
            .with_blocks(_rod_pos)
            .with_recompute(
                hepworth::block::make_rod_edge_stretch_recompute<0>(
                    __R, _cfg.edge_stretch_w))
            .with_recompute(hepworth::block::make_rod_smooth_recompute<0>(
                __R, _cfg.smooth_w))
            .with_recompute(hepworth::block::make_rod_min_kink_recompute<0>(
                __R, _cfg.min_kink_w))
            .with_recompute(
                hepworth::block::make_rod_collisions_recompute<0>(
                    __R, __Rd, _cfg.collision_w))
            .dt(_cfg.dt)
            .damping(_cfg.damping)
            .iterations(_cfg.iterations)
            .build();

    _tangent = _graph.create_node<tangent_point_gradient_node>(
        __R, __Rd, _cfg.tangent);
    _solver_node =
        _graph.create_node<hepworth::block::block_solver_node<
            hepworth::block::rod_position_block>>(_solver_config);
    _graph.link(_tangent->output(), _solver_node->input_at<0>());
  }

  const braid_planar_demo_config &config() const { return _cfg; }
  const asawa::rod::rod &rod() const { return *__R; }
  const asawa::rod::rod &planar_rod() const { return *__R_planar; }
  asawa::rod::rod::ptr rod_ptr() const { return __R; }

  void step(int /*frame*/) {
    _graph.run();
    __Rd->step();
    // Cosserat quats are not DOFs; rebuild tube frames by parallel transport.
    __R->_update_frames();
  }

private:
  braid_planar_demo_config _cfg;
  asawa::rod::rod::ptr __R_planar;
  asawa::rod::rod::ptr __R;
  asawa::rod::dynamic::ptr __Rd;

  hepworth::block::rod_position_block::ptr _rod_pos;
  hepworth::block::block_solver_config<hepworth::block::rod_position_block>
      _solver_config;
  liblombardi::GraphContext _graph;
  tangent_point_gradient_node::ptr _tangent;
  hepworth::block::block_solver_node<hepworth::block::rod_position_block>::ptr
      _solver_node;
};

class braid_planar_demo_adapter : public demo_adapter<braid_planar_demo> {
public:
  static ptr create(const braid_planar_demo_config &cfg =
                        k_braid_planar_demo_default) {
    return std::make_shared<braid_planar_demo_adapter>(cfg);
  }

  explicit braid_planar_demo_adapter(
      const braid_planar_demo_config &cfg = k_braid_planar_demo_default)
      : demo_adapter<braid_planar_demo>(braid_planar_demo::create(cfg),
                                        "braid_sphere") {}

  std::optional<mesh_snapshot> rod_mesh() const override {
    return make_rod_mesh_snapshot(demo()->rod(), demo()->config().rod_color);
  }

  std::optional<rod_snapshot> rod_polyline() const override {
    if (!demo()->config().show_planar_graph) {
      return std::nullopt;
    }
    return make_rod_snapshot(demo()->planar_rod(),
                             demo()->config().planar_color);
  }
};

} // namespace duchamp
} // namespace gaudi

#endif
