#pragma once

/// An Engman scene with drawing. The scene and its terms live in
/// `gaudi/engman/scene.hpp` and `gaudi/engman/terms.hpp`; this file only draws.
///
/// Pick a preset with `ENGMAN_PRESET` (`torsion` default, `torsion_hopf`,
/// `torsion_elastic`, `twist_rate`, `relax`), then edit `engman_demo_config`.

#include <array>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

#include "gaudi/common.h"
#include "gaudi/duchamp/demo_trait.hpp"
#include "gaudi/engman/scene.hpp"
#include "gaudi/geometry_logger.hpp"

namespace gaudi {
namespace duchamp {

/// The configuration the demo starts from: the preset named by
/// `ENGMAN_PRESET`, then any edits here.
inline engman::scene_config engman_demo_config() {
  const char *name = std::getenv("ENGMAN_PRESET");
  engman::scene_config cfg = engman::preset(name ? std::string(name) : std::string());
  // Edits on top of the preset, for example:
  //   engman::find<engman::terms::torsion>(cfg.terms)->tau = -1; // other hand
  //   cfg.terms.push_back(engman::terms::face_area{.m = engman::mode::soft});
  // or a whole scene (fields in declaration order):
  //   cfg = {.rest = {.relax = 0},
  //          .terms = {engman::terms::coupling{},
  //                    engman::terms::memory{},
  //                    engman::terms::bending{.w = 0.1},
  //                    engman::terms::torsion{.tau = 1, .axis = {.kind = engman::axis_kind::winding}},
  //                    engman::terms::face_area{},
  //                    engman::terms::conformal{.mu = 1}}};
  return cfg;
}

class engman_demo {
public:
  using ptr = std::shared_ptr<engman_demo>;
  static ptr create() { return create(engman_demo_config()); }
  static ptr create(const engman::scene_config &cfg) { return std::make_shared<engman_demo>(cfg); }

  explicit engman_demo(const engman::scene_config &cfg) : scene(engman::scene::create(cfg)) {}

  void step(int frame) {
    (void)frame;
    scene->step();
  }

  void draw_debug() const {
    const engman::draw_config &d = scene->cfg.draw;
    if (d.frames) {
      std::vector<engman::scene::vertex_frame> fr;
      scene->vertex_frames(fr);
      real kbar = 0;
      for (const auto &f : fr)
        kbar += (std::abs(f.k0) + std::abs(f.k1)) / 2;
      kbar = fr.empty() ? real(1) : std::max(kbar / real(fr.size()), real(1e-12));
      const real L = d.frame_length * scene->mean_edge;
      const vec4 red(1, 0, 0, 1), green(0, 1, 0, 1), blue(0, 0, 1, 1);
      for (const auto &f : fr) {
        const real l0 = d.frames_curvature ? L * std::abs(f.k0) / kbar : L;
        const real l1 = d.frames_curvature ? L * std::abs(f.k1) / kbar : L;
        geometry_logger::line(f.x - l0 * f.d0, f.x + l0 * f.d0, red);
        geometry_logger::line(f.x - l1 * f.d1, f.x + l1 * f.d1, green);
        geometry_logger::line(f.x - L * f.n, f.x + L * f.n, blue);
      }
    }
    if (d.principal_lines) {
      std::vector<std::array<vec3, 2>> seg;
      scene->principal_segments(seg);
      const vec4 col(0.1, 0.1, 0.1, 1.0);
      const real radius = real(0.08) * scene->mean_edge;
      for (const auto &sg : seg)
        geometry_logger::line(sg[0], sg[1], col, radius);
    }
  }

  mesh_snapshot snapshot() const {
    mesh_snapshot m;
    m.positions = scene->x;
    m.indices.reserve(3 * scene->tris.size());
    for (const auto &t : scene->tris)
      for (int k = 0; k < 3; ++k)
        m.indices.push_back(static_cast<uint32_t>(t[static_cast<size_t>(k)]));
    m.colors = scene->colors();
    return m;
  }

  engman::scene::ptr scene;
};

} // namespace duchamp
} // namespace gaudi
