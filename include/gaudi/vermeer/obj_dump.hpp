#pragma once

#include <string>

#include "gaudi/duchamp/obj_export.hpp"
#include "gaudi/vermeer/duchamp_playback.hpp"
#include "gaudi/vermeer/scene_frame.hpp"

namespace gaudi {
namespace vermeer {

// If playback requested an export, write dump/<demo>_f<N>.{obj,mtl} from frame.
inline void maybe_dump_scene_obj(duchamp_playback &pb, const SceneFrame &frame,
                                 const std::string &demo_name) {
  if (!pb.export_obj_once.exchange(false))
    return;

  const std::string stem =
      "dump/" + demo_name + "_f" + std::to_string(frame.sim_frame);
  duchamp::write_scene_meshes_obj(stem, frame.shell, frame.rod);
}

} // namespace vermeer
} // namespace gaudi
