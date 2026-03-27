#ifndef GAUDI_ASAWA_SHELL_SHELL_ID_HPP
#define GAUDI_ASAWA_SHELL_SHELL_ID_HPP

#include "gaudi/strong_id.hpp"

namespace gaudi {
namespace asawa {
namespace shell {

enum class ShellIdKind : unsigned char { corner, face, vert };

using CornerId = Id<ShellIdKind::corner, int>;
using FaceId = Id<ShellIdKind::face, int>;
using VertId = Id<ShellIdKind::vert, int>;

/// Explicit construction at loaders and raw index boundaries only.
inline constexpr CornerId corner_id(int i) { return CornerId{i}; }
inline constexpr FaceId face_id(int i) { return FaceId{i}; }
inline constexpr VertId vert_id(int i) { return VertId{i}; }

/// Shell uses negative sentinel ids (e.g. -1) for invalid vertices.
inline constexpr bool vert_id_valid(VertId v) noexcept {
  return v >= vert_id(0);
}

} // namespace shell
} // namespace asawa
} // namespace gaudi

#endif
