#ifndef GAUDI_DUCHAMP_BRAID_CIRCLE_BUNDLE_HPP
#define GAUDI_DUCHAMP_BRAID_CIRCLE_BUNDLE_HPP

#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/common.h"
#include "gaudi/duchamp/braid_planar_rod.hpp"
#include "gaudi/duchamp/braid_sphere_rod.hpp"
#include "gaudi/windychien/braid.hpp"
#include "gaudi/windychien/catalog.hpp"

#include <string>

namespace gaudi {
namespace duchamp {

/// Braid knot + matching sphere shell, built in the same ambient circle.
struct braid_circle_config {
  /// If true: plain weave (even frames 0-1,2-3,…; odd 1-2,3-4,…).
  /// If false: load `knot_name` from the catalog.
  bool use_plain_weave = true;
  int weave_strands = 50;
  int weave_frames = 50;
  std::string knot_name = "K11a359";
  braid_planar_params planar{};
  // Keep dense circumferential subdiv even with fewer weave frames.
  braid_sphere_params sphere{.n_lat = 256};
  real rod_radius = 0.01;
  /// Shell tessellation (u = longitude, v = latitude stacks).
  int shell_u = 256;
  int shell_v = 128;
};

struct braid_circle_bundle {
  asawa::shell::shell::ptr shell;
  asawa::rod::rod::ptr rod;
  braid_circle_config config;
};

/// Build sphere shell at `sphere.radius` and closed braid rod on that circle.
inline braid_circle_bundle
make_braid_circle_bundle(const braid_circle_config &cfg = {}) {
  braid_planar_params planar = cfg.planar;
  planar.close = true;
  planar.center = false;
  // Flush rod radius + small air gap at crossings.
  planar.eps_z = cfg.rod_radius + real(0.1) * cfg.rod_radius;
  if (planar.dx <= 0.0) {
    planar.dx = 1.0;
  }
  if (planar.dy <= 0.0) {
    planar.dy = planar.dx;
  }

  const windychien::braid b =
      cfg.use_plain_weave
          ? windychien::make_plain_weave(cfg.weave_strands, cfg.weave_frames)
          : windychien::load_braid(cfg.knot_name);

  braid_circle_bundle out;
  out.config = cfg;
  out.shell = asawa::shell::load_sphere(cfg.sphere.radius, cfg.shell_u,
                                       cfg.shell_v);
  out.rod = braid_to_sphere_rod(b, planar, cfg.sphere);
  out.rod->_r = cfg.rod_radius;
  return out;
}

} // namespace duchamp
} // namespace gaudi

#endif
