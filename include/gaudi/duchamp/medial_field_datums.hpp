#ifndef __GAUDI_DUCHAMP_MEDIAL_FIELD_DATUMS__
#define __GAUDI_DUCHAMP_MEDIAL_FIELD_DATUMS__

#include "gaudi/albers/darboux_cyclide.hpp"
#include "gaudi/duchamp/field_nodes.hpp"

namespace gaudi {
namespace duchamp {

using cyclide_field_datum = field_datum<albers::vec14>;
using line_field_datum = field_datum<vec6>;
using frame_field_datum = field_datum<mat3>;

} // namespace duchamp
} // namespace gaudi

#endif // __GAUDI_DUCHAMP_MEDIAL_FIELD_DATUMS__
