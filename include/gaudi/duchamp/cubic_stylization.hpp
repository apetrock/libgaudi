#pragma once

#include <memory>
#include <string>
#include <vector>

#include "gaudi/asawa/datums.hpp"
#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/asawa/shell/dynamic.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/duchamp/fields.hpp"
#include "gaudi/hepworth/blocks/shell_position_block.hpp"
#include "gaudi/hepworth/constraints/bundles.hpp"
#include "gaudi/hepworth/nodes/solver_builder.hpp"

namespace gaudi {
namespace duchamp {

/// Halloween cubification: spoke ARAP + hard-quantized normal tip-in
/// (Liu–Jacobson style local step; Hepworth ADMM as the outer loop).
/// N_rest rides remesh as a vert datum; N_goal is re-quantized each step
/// from current vertex normals.
class cubic_stylization {
public:
  using ptr = std::shared_ptr<cubic_stylization>;

  static ptr create(const std::string &mesh_path = {}) {
    return std::make_shared<cubic_stylization>(mesh_path);
  }

  explicit cubic_stylization(const std::string &mesh_path = {})
      : _mesh_path(mesh_path) {
    reset_scene();
  }

  void reset_scene() {
    if (_mesh_path.empty())
      __M = asawa::shell::load_bunny();
    else
      __M = asawa::shell::load_obj(_mesh_path);

    asawa::shell::triangulate(*__M);
    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    asawa::center(x);

    asawa::init_vert_datum<vec3>(*__M, vec3::Zero()); // velocities @ 1
    auto xs = std::make_shared<shell_vert_positions>(__M, 0);
    auto vs = std::make_shared<shell_vert_velocities>(__M, 1);
    _shell = std::make_shared<hepworth::block::shell_position_block>(__M, xs, vs);

    const std::vector<vec3> N0 = asawa::shell::vertex_normals(*__M, xs->get());
    _iN_rest = asawa::init_vert_datum<vec3>(*__M, vec3::Zero());
    _iN_goal = asawa::init_vert_datum<vec3>(*__M, vec3::Zero());
    asawa::get_vec_data(*__M, _iN_rest) = N0;
    asawa::get_vec_data(*__M, _iN_goal) =
        hepworth::block::quantize_normals_axis(N0);

    const real l0 = asawa::shell::avg_length(*__M, xs->get());
    __surf = asawa::shell::dynamic::create(__M, l0, 1.25 * l0, 0.5 * l0);

    _stretch_w = 0.1;
    _cubify_w = 4.0;
    _lambda = 1.0e-1;
    _pin_w = 1.0e-1;
    _solver_dt = 0.01;
    _solver_iters = 16;

    _config =
        hepworth::block::block_solver_builder<hepworth::block::shell_position_block>::create()
            .with_blocks(_shell)
            .with_recompute(hepworth::block::make_shell_stretch_recompute<0>(
                _shell, _stretch_w))
            .with_bundle(hepworth::block::make_shell_cubification_bundle<0>(
                _shell, _iN_rest, _iN_goal, _cubify_w, _lambda))
            .with_bundle(
                hepworth::block::make_shell_pin_bundle<0>(_shell, _pin_w))
            .dt(_solver_dt)
            .damping(0.2)
            .iterations(_solver_iters)
            .build();
  }

  void step(int frame) {
    (void)frame;
    refresh_goal_normals();
    hepworth::block::run_solver_step(_config, _solver);
    __surf->step(true);
  }

  asawa::shell::shell::ptr __M;
  asawa::shell::dynamic::ptr __surf;

private:
  void refresh_goal_normals() {
    const std::vector<vec3> &x = _shell->xs->get();
    const std::vector<vec3> N =
        asawa::shell::vertex_normals(*__M, x);
    asawa::get_vec_data(*__M, _iN_goal) =
        hepworth::block::quantize_normals_axis(N);
  }

  std::string _mesh_path;
  hepworth::block::shell_position_block::ptr _shell;
  index_t _iN_rest = -1;
  index_t _iN_goal = -1;

  real _stretch_w = 0.5;
  real _cubify_w = 1.0;
  real _lambda = 0.1;
  real _pin_w = 1.0e-2;
  real _solver_dt = 0.01;
  int _solver_iters = 16;

  hepworth::block::block_solver_config<hepworth::block::shell_position_block>
      _config;
  hepworth::block::projection_solver _solver;
};

} // namespace duchamp
} // namespace gaudi
