#pragma once

#include "datum_pool.hpp"

namespace liblombardi {

/// Port direction
enum class PortDirection {
    Input,
    Output
};

/// Handle to a datum in the pool
struct PortHandle {
    datum_index_t index;
    bool is_valid() const { return index != static_cast<datum_index_t>(-1); }

    bool operator==(const PortHandle& other) const {
        return index == other.index;
    }

    bool operator!=(const PortHandle& other) const {
        return index != other.index;
    }
};

/// Default invalid handle
inline constexpr PortHandle invalid_handle{static_cast<datum_index_t>(-1)};

/// Port wrapper that holds a handle and direction
template <typename T>
struct Port {
    PortHandle handle;

    Port() : handle(invalid_handle) {}
    Port(datum_index_t idx) : handle({idx}) {}
    Port(const PortHandle& h) : handle(h) {}

    /// Check if handle is valid
    bool is_valid() const { return handle.is_valid(); }

    /// Get datum index
    datum_index_t get_index() const { return handle.index; }

    /// Set handle
    void set_handle(datum_index_t idx) { handle.index = idx; }

    /// Set from PortHandle
    void set_handle(const PortHandle& h) { handle = h; }
};

/// Alias for input port
template <typename T>
using InputPort = Port<T>;

/// Alias for output port
template <typename T>
using OutputPort = Port<T>;

} // namespace liblombardi
