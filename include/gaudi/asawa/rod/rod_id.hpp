#ifndef GAUDI_ASAWA_ROD_ROD_ID_HPP
#define GAUDI_ASAWA_ROD_ROD_ID_HPP

#include "gaudi/strong_id.hpp"

namespace gaudi {
namespace asawa {
namespace rod {

enum class RodIdKind : unsigned char { corner };

using CornerId = Id<RodIdKind::corner, int>;

inline constexpr CornerId corner_id(int i) { return CornerId{i}; }

} // namespace rod
} // namespace asawa
} // namespace gaudi

#endif
