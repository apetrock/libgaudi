#ifndef GAUDI_STRONG_ID_HPP
#define GAUDI_STRONG_ID_HPP

#include <compare>
#include <functional>
#include <ostream>

namespace gaudi {

/// Distinct wrapper per \p K (usually an \c enum class enumerator).
/// C++17 \c auto non-type template parameter; domain enums stay in subsystem
/// headers (e.g. shell_id.hpp).
///
/// Implicit conversion to \c T supports indexing and arithmetic at use sites;
/// construction from \c T remains explicit via \c Id{...} or subsystem
/// factories (e.g. \c vert_id(i)).
template <auto K, typename T = int> struct Id {
private:
  T value_{};

public:
  explicit constexpr Id(T v) noexcept : value_(v) {}

  constexpr operator T() const noexcept { return value_; }

  constexpr auto operator<=>(const Id &) const noexcept = default;

  friend constexpr std::ostream &operator<<(std::ostream &os, Id id) {
    return os << static_cast<T>(id);
  }
};

} // namespace gaudi

namespace std {

template <auto K, typename T> struct hash<gaudi::Id<K, T>> {
  std::size_t operator()(gaudi::Id<K, T> id) const noexcept {
    return std::hash<T>{}(static_cast<T>(id));
  }
};

} // namespace std

#endif
