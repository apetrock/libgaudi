#pragma once

#include "gaudi/vermeer/duchamp_project.hpp"
#include "gaudi/vermeer/vermeer_config.hpp"

namespace gaudi {
namespace vermeer {

class medial_axis_project : public duchamp_project {
public:
  explicit medial_axis_project(duchamp::demo_trait::ptr demo)
      : duchamp_project(std::move(demo),
                        vermeer_config{.renderer_type = render_graph_type::ssao_deferred}) {}
};

} // namespace vermeer
} // namespace gaudi
