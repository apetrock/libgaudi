#pragma once

#include <type_traits>
#include <cstddef>

namespace liblombardi {

using uint = unsigned int;

/// Port Definition Template
/// Defines a port with its element type and static ID
template <typename T, auto PortId>
struct PortDef {
    using buffer_type = T;          // Element type (int, vec3, etc.)
    static constexpr auto id = PortId;  // Local ID within node

    // Compile-time assertions to ensure valid port definition
    static_assert(std::is_same_v<T, std::remove_reference_t<T>>,
                  "PortDef buffer_type must not be a reference");
    static_assert(!std::is_void_v<T>, "PortDef buffer_type cannot be void");
};

/// Convert a PortDef::id (enum or integral) to a local port index.
template <auto Id>
inline constexpr uint port_index_v =
    static_cast<uint>(static_cast<std::underlying_type_t<decltype(Id)>>(Id));

/// Same as port_index_v but accepts a PortDef type (works with dependent types).
template <typename TPortDef>
inline constexpr uint port_index_from_def_v =
    static_cast<uint>(static_cast<std::underlying_type_t<decltype(TPortDef::id)>>(TPortDef::id));

/// Bundles a node reference with its port definition for type-safe linking.
template <typename TNode, typename TPortDef>
struct PortRef {
    TNode& node;
    using port_def = TPortDef;
};

} // namespace liblombardi
