#ifndef __GAUDI_HEPWORTH_DOF_INPUTS__
#define __GAUDI_HEPWORTH_DOF_INPUTS__

#include "gaudi/common.h"
#include "gaudi/duchamp/field_nodes.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

struct forces_tag {};
struct torques_tag {};
struct velocity_tag {};

template <typename DatumT, typename Tag>
struct dof_input {
  using datum_type = DatumT;
  using value_type = typename DatumT::value_type;
  using tag = Tag;
};

using vec3_forces_port = dof_input<duchamp::field_datum<vec3>, forces_tag>;
using vec3_torques_port = dof_input<duchamp::field_datum<vec3>, torques_tag>;
using vec3_velocity_port = dof_input<duchamp::field_datum<vec3>, velocity_tag>;

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // __GAUDI_HEPWORTH_DOF_INPUTS__
