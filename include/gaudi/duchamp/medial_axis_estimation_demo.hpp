#ifndef __GAUDI_DUCHAMP_MEDIAL_AXIS_ESTIMATION_DEMO_HPP__
#define __GAUDI_DUCHAMP_MEDIAL_AXIS_ESTIMATION_DEMO_HPP__

#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/duchamp/pierce_cylinder_medial.hpp"
#include "gaudi/geometry_logger.hpp"

#include <iostream>
#include <memory>

namespace gaudi {
namespace duchamp {

class medial_axis_estimation_demo {
public:
  using ptr = std::shared_ptr<medial_axis_estimation_demo>;

  struct cached_debug_seg {
    vec3 p0 = vec3::Zero();
    vec3 p1 = vec3::Zero();
    vec4 color = vec4::Zero();
    real radius = geometry_logger::k_default_line_radius;
  };

  struct cached_debug_point {
    vec3 p = vec3::Zero();
    vec4 color = vec4::Zero();
  };

  static ptr create(pierce_cylinder_medial_params params =
                        default_pierce_cylinder_medial_params()) {
    return std::make_shared<medial_axis_estimation_demo>(params);
  }

  explicit medial_axis_estimation_demo(pierce_cylinder_medial_params params)
      : _params(params) {
    _mesh = asawa::shell::load_bunny();
    asawa::shell::triangulate(*_mesh);
    normalize_medial_demo_mesh(*_mesh);
    configure_scene_frame();
    _samples = compute_pierce_cylinder_medial(*_mesh, _params, &_stats);
    print_pierce_cylinder_medial_stats(std::cerr, _stats);
    build_debug_cache();
  }

  void step(int frame) { _frame = frame; }

  void draw_debug() {
    for (const cached_debug_point &pt : _point_cache) {
      geometry_logger::point(pt.p, pt.color);
    }
    for (const cached_debug_seg &seg : _debug_cache) {
      geometry_logger::line(seg.p0, seg.p1, seg.color, seg.radius);
    }
  }

  int frame() const { return _frame; }
  const pierce_cylinder_medial_stats &stats() const { return _stats; }
  const std::vector<pierce_cylinder_sample> &samples() const {
    return _samples;
  }
  asawa::shell::shell::ptr shell() const { return _mesh; }

private:
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
    std::cerr << "loaded bunny verts=" << _mesh->vert_count()
              << " faces=" << _mesh->face_count()
              << " scene_radius=" << _major_radius << std::endl;
  }

  void append_seg(const vec3 &p0, const vec3 &p1, const vec4 &color,
                  real radius) {
    if ((p1 - p0).norm() < 1e-12) {
      return;
    }
    _debug_cache.push_back({p0, p1, color, radius});
  }

  void build_debug_cache() {
    const vec4 midpoint_color(0.0, 0.85, 1.0, 1.0);
    const vec4 cylinder_axis_color(0.9, 0.15, 1.0, 1.0);
    const vec4 rejected_axis_color(1.0, 0.85, 0.05, 1.0);
    const vec4 pierce_miss_color(1.0, 0.15, 0.15, 1.0);
    constexpr real k_radius = geometry_logger::k_default_line_radius;

    _debug_cache.clear();
    _point_cache.clear();

    const std::vector<vec3> &x = asawa::const_get_vec_data(*_mesh, 0);
    const real avg_len = asawa::shell::avg_length(*_mesh, x);
    const real miss_len = _params.miss_ray_scale * avg_len;
    const real axis_len = _params.axis_display_scale * avg_len;

    for (const pierce_cylinder_sample &sample : _samples) {
      const vec3 &p = sample.foot;

      if (sample.status == pierce_cylinder_status::pierce_miss) {
        // Always draw pierce misses when a foot normal exists (no reject filter).
        if (sample.foot_normal.squaredNorm() < 1e-24) {
          continue;
        }
        append_seg(p, p - miss_len * sample.foot_normal, pierce_miss_color,
                   k_radius);
        continue;
      }

      _point_cache.push_back({sample.midpoint, midpoint_color});

      const bool draw_axis =
          sample.axis.squaredNorm() > 1e-24 && sample.center.allFinite();
      if (!draw_axis) {
        continue;
      }

      const vec3 dp = axis_len * sample.axis.normalized();
      const vec4 &axis_col =
          sample.accepted ? cylinder_axis_color : rejected_axis_color;
      append_seg(sample.center - dp, sample.center + dp, axis_col, k_radius);
    }
  }

  pierce_cylinder_medial_params _params;
  pierce_cylinder_medial_stats _stats;
  std::vector<pierce_cylinder_sample> _samples;
  std::vector<cached_debug_seg> _debug_cache;
  std::vector<cached_debug_point> _point_cache;
  asawa::shell::shell::ptr _mesh;
  vec3 _center = vec3::Zero();
  real _major_radius = 1.0;
  int _frame = 0;
};

} // namespace duchamp
} // namespace gaudi

#endif // __GAUDI_DUCHAMP_MEDIAL_AXIS_ESTIMATION_DEMO_HPP__
