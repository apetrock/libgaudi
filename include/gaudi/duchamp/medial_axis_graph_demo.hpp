#ifndef __GAUDI_DUCHAMP_MEDIAL_AXIS_GRAPH_DEMO__
#define __GAUDI_DUCHAMP_MEDIAL_AXIS_GRAPH_DEMO__

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>

#include "gaudi/albers/darboux_medial_geometry.hpp"
#include "gaudi/albers/line_cylinder.hpp"
#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/duchamp/body_datum.hpp"
#include "gaudi/duchamp/calder_graph_nodes.hpp"
#include "gaudi/duchamp/darboux_cyclide_medial.hpp"
#include "gaudi/duchamp/field_graph_nodes.hpp"
#include "gaudi/duchamp/medial_graph_nodes.hpp"
#include "gaudi/duchamp/medial_result_types.hpp"
#include "gaudi/geometry_logger.hpp"

#include "liblombardi/graph_context.hpp"

namespace gaudi {
namespace duchamp {

enum class medial_viz_mode {
  HessianFrame,
  CylinderFromMesh,
  CylinderFromMedialPoints
};

template <medial_viz_mode VizMode = medial_viz_mode::HessianFrame>
class medial_axis_graph_demo {
public:
  using ptr = std::shared_ptr<medial_axis_graph_demo>;

  static ptr create(cyclide_medial_params params =
                        default_cyclide_medial_demo_params()) {
    return std::make_shared<medial_axis_graph_demo>(params);
  }

  explicit medial_axis_graph_demo(cyclide_medial_params params)
      : _params(params) {
    _mesh = asawa::shell::load_bunny();
    asawa::shell::triangulate(*_mesh);
    normalize_cyclide_medial_demo_mesh(*_mesh);
    configure_scene_frame();
    wire_graph();
    _graph.run();
    summarize_results();
  }

  void step(int frame) {
    _frame = frame;
    draw_medial_viz();
  }

  int frame() const { return _frame; }
  const cyclide_medial_stats &stats() const { return _stats; }
  asawa::shell::shell::ptr shell() const { return _mesh; }

private:
  void wire_graph() {
    const std::vector<vec3> &x = asawa::const_get_vec_data(*_mesh, 0);
    const real avg_len = asawa::shell::avg_length(*_mesh, x);
    _l0 = std::max(_params.l0_scale * avg_len, real(1e-12));
    _max_travel = _params.max_travel_scale * avg_len;

    auto body = _body = _graph.create_node<body_constant_node>(
        make_shell_body(_mesh));
    auto positions = _positions =
        _graph.create_node<position_snapshot_node>();
    auto normals = _normals =
        _graph.create_node<vertex_normals_snapshot_node>();
    auto cyclide_fit = _cyclide_fit =
        _graph.create_node<darboux_cyclide_fit_node>(_l0, _params.fit_p,
                                                     _params.fit_w0);
    auto cyclide_smooth = _cyclide_smooth =
        _graph.create_node<darboux_cyclide_smooth_node>(_params.smooth);
    auto medial_search = _medial_search =
        _graph.create_node<medial_shape_energy_search_node>(_max_travel);
    auto hessian = _hessian = _graph.create_node<hessian_frame_node>();
    auto cyl_mesh = _cyl_mesh =
        _graph.create_node<cylinder_fit_mesh_node>(_params.normal_l0,
                                                   _params.fit_p);
    auto cyl_points = _cyl_points =
        _graph.create_node<cylinder_fit_points_node>();

    _graph.link(body->output(), positions->body());
    _graph.link(body->output(), normals->body());
    _graph.link(body->output(), cyclide_fit->body());
    _graph.link(positions->output(), cyclide_fit->pov());
    _graph.link(normals->output(), cyclide_fit->n_pov());
    _graph.link(body->output(), cyclide_smooth->body());
    _graph.link(positions->output(), cyclide_smooth->pov());
    _graph.link(cyclide_fit->output(), cyclide_smooth->cyclide_in());
    _graph.link(positions->output(), medial_search->pov());
    if (_params.enable_cyclide_smooth) {
      _graph.link(cyclide_smooth->cyclide_out(), medial_search->cyclide());
      _graph.link(cyclide_smooth->cyclide_out(), hessian->cyclide());
    } else {
      _graph.link(cyclide_fit->output(), medial_search->cyclide());
      _graph.link(cyclide_fit->output(), hessian->cyclide());
    }
    _graph.link(normals->output(), medial_search->n_pov());
    _graph.link(positions->output(), hessian->pov());
    _graph.link(medial_search->output(), hessian->position());
    _graph.link(body->output(), cyl_mesh->body());
    _graph.link(positions->output(), cyl_mesh->pov());
    _graph.link(normals->output(), cyl_mesh->n_pov());
    _graph.link(positions->output(), cyl_points->pov());
    _graph.link(normals->output(), cyl_points->n_pov());
    _graph.link(medial_search->output(), cyl_points->points());
  }

  void summarize_results() {
    const auto &mask =
        _medial_search
            ->get_datum<medial_shape_energy_search_node::MaskPortDef>()
            ->data();
    const auto &points =
        _medial_search
            ->get_datum<medial_shape_energy_search_node::OutputPortDef>()
            ->data();
    const auto &pov =
        _positions->get_datum<position_snapshot_node::OutputPortDef>()->data();
    _stats.total = static_cast<int>(mask.size());
    _stats.accepted = count_medial_accepted(mask);
    _stats.rejected = _stats.total - _stats.accepted;
    real travel_sum = 0.0;
    for (size_t i = 0; i < mask.size(); ++i) {
      if (mask[i] <= real(0.5) || i >= points.size() || i >= pov.size()) {
        continue;
      }
      travel_sum += (points[i] - pov[i]).norm();
    }
    if (_stats.accepted > 0) {
      _stats.avg_travel = travel_sum / real(_stats.accepted);
    }
    std::cerr << "medial graph demo (ShapeEnergy) accepted=" << _stats.accepted
              << " rejected=" << _stats.rejected
              << " avg_travel=" << _stats.avg_travel << std::endl;
  }

  void configure_scene_frame() {
    const std::vector<vec3> &x = asawa::const_get_vec_data(*_mesh, 0);
    if (x.empty()) {
      return;
    }
    vec3 lo = x.front();
    vec3 hi = x.front();
    for (const vec3 &p : x) {
      lo = lo.cwiseMin(p);
      hi = hi.cwiseMax(p);
    }
    _center = 0.5 * (lo + hi);
    _major_radius = 0.5 * (hi - lo).norm();
    _minor_radius = 4.0 * asawa::shell::avg_length(*_mesh, x);
  }

  void draw_medial_viz() {
    const vec4 axis_color(0.35, 0.35, 0.35, 1.0);
    geometry_logger::line(_center - 1.7 * _major_radius * vec3::UnitZ(),
                          _center + 1.7 * _major_radius * vec3::UnitZ(),
                          axis_color);

    switch (_params.display) {
    case medial_axis_display::MedialAxis:
      draw_medial_axis_overlay();
      break;
    case medial_axis_display::SmoothedFitHessian:
      draw_smoothed_fit_hessian();
      break;
    }
  }

  void draw_medial_axis_overlay() {
    const vec4 surface_color(0.0, 0.85, 1.0, 1.0);
    const vec4 medial_color(1.0, 0.65, 0.05, 1.0);
    const vec4 cylinder_color(0.9, 0.15, 1.0, 1.0);
    const vec4 rejected_color(0.35, 0.35, 0.35, 0.45);

    const auto &pov =
        _positions->get_datum<position_snapshot_node::OutputPortDef>()->data();
    const auto &medial =
        _medial_search
            ->get_datum<medial_shape_energy_search_node::OutputPortDef>()
            ->data();
    const auto &mask =
        _medial_search
            ->get_datum<medial_shape_energy_search_node::MaskPortDef>()
            ->data();
    const real axis_len = 8.0 * asawa::shell::avg_length(*_mesh, pov);

    for (size_t i = 0; i < pov.size(); ++i) {
      if (i >= mask.size()) {
        continue;
      }
      const vec3 p = pov[i];
      if (mask[i] <= real(0.5)) {
        geometry_logger::point(p, rejected_color);
        continue;
      }
      const vec3 m = medial[i];
      const real travel = (m - p).norm();
      if (travel > 1e-12) {
        geometry_logger::line(p, m, surface_color);
      } else {
        geometry_logger::point(m, surface_color);
      }
      geometry_logger::point(m, medial_color);

      if constexpr (VizMode == medial_viz_mode::HessianFrame) {
        const auto &frames =
            _hessian->get_datum<hessian_frame_node::OutputPortDef>()->data();
        if (i >= frames.size()) {
          continue;
        }
        const mat3 &F = frames[i];
        const real len = 0.35 * travel;
        for (int ax = 0; ax < 3; ++ax) {
          const vec3 dir = F.col(ax).normalized();
          geometry_logger::line(m - len * dir, m + len * dir, medial_color);
        }
      } else if constexpr (VizMode == medial_viz_mode::CylinderFromMesh) {
        const auto &lines =
            _cyl_mesh->get_datum<cylinder_fit_mesh_node::OutputPortDef>()
                ->data();
        if (i >= lines.size()) {
          continue;
        }
        const vec3 d = albers::plucker_line_direction(lines[i]);
        if (d.squaredNorm() > 1e-12) {
          const vec3 dp = axis_len * d.normalized();
          geometry_logger::line(p - dp, p + dp, cylinder_color);
        }
      } else if constexpr (VizMode ==
                           medial_viz_mode::CylinderFromMedialPoints) {
        const auto &lines =
            _cyl_points->get_datum<cylinder_fit_points_node::OutputPortDef>()
                ->data();
        if (i >= lines.size()) {
          continue;
        }
        const vec3 d = albers::plucker_line_direction(lines[i]);
        if (d.squaredNorm() > 1e-12) {
          const vec3 dp = axis_len * d.normalized();
          geometry_logger::line(m - dp, m + dp, cylinder_color);
        }
      }
    }
  }

  void draw_smoothed_fit_hessian() {
    const vec4 foot_color(0.85, 0.85, 0.85, 0.65);
    const vec4 axis_colors[3] = {vec4(1.0, 0.0, 0.0, 1.0),
                                  vec4(0.0, 1.0, 0.0, 1.0),
                                  vec4(0.0, 0.0, 1.0, 1.0)};

    const auto &pov =
        _positions->get_datum<position_snapshot_node::OutputPortDef>()->data();
    const auto &Q =
        _params.enable_cyclide_smooth
            ? _cyclide_smooth
                  ->get_datum<darboux_cyclide_smooth_node::CyclideOutPortDef>()
                  ->data()
            : _cyclide_fit->get_datum<darboux_cyclide_fit_node::OutputPortDef>()
                  ->data();

    const real len =
        _params.hessian_frame_scale * asawa::shell::avg_length(*_mesh, pov);
    const real line_radius = _params.hessian_line_radius;
    const size_t n = std::min(pov.size(), Q.size());

    for (size_t i = 0; i < n; ++i) {
      const vec3 p = pov[i];
      geometry_logger::point(p, foot_color);

      const mat3 W = albers::shape_operator_at(Q[i], vec3::Zero());
      if (!W.allFinite()) {
        continue;
      }

      Eigen::SelfAdjointEigenSolver<mat3> es(W);
      if (es.info() != Eigen::Success) {
        continue;
      }

      const mat3 &V = es.eigenvectors();
      for (int ax = 0; ax < 3; ++ax) {
        const vec3 dir = V.col(ax);
        if (dir.squaredNorm() < 1e-24) {
          continue;
        }
        const vec3 d = dir.normalized();
        geometry_logger::line(p - len * d, p + len * d, axis_colors[ax],
                              line_radius);
      }
    }
  }

  cyclide_medial_params _params;
  cyclide_medial_stats _stats;
  asawa::shell::shell::ptr _mesh;
  liblombardi::GraphContext _graph;
  body_constant_node::ptr _body;
  position_snapshot_node::ptr _positions;
  vertex_normals_snapshot_node::ptr _normals;
  darboux_cyclide_fit_node::ptr _cyclide_fit;
  darboux_cyclide_smooth_node::ptr _cyclide_smooth;
  medial_shape_energy_search_node::ptr _medial_search;
  hessian_frame_node::ptr _hessian;
  cylinder_fit_mesh_node::ptr _cyl_mesh;
  cylinder_fit_points_node::ptr _cyl_points;
  real _l0 = 1.0;
  real _max_travel = 1.0;
  vec3 _center = vec3::Zero();
  real _major_radius = 1.0;
  real _minor_radius = 0.35;
  int _frame = 0;
};

} // namespace duchamp
} // namespace gaudi

#endif // __GAUDI_DUCHAMP_MEDIAL_AXIS_GRAPH_DEMO__
