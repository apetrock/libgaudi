#ifndef __GAUDI_DUCHAMP_BODY_DATUM__
#define __GAUDI_DUCHAMP_BODY_DATUM__

#include <memory>
#include <stdexcept>
#include <vector>

#include "gaudi/asawa/datums.hpp"
#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/common.h"

#include "liblombardi/datum_pool.hpp"

namespace gaudi {
namespace duchamp {

enum class body_kind { rod, shell };

struct body_handle {
  virtual ~body_handle() = default;
  virtual body_kind kind() const = 0;
};

struct rod_body_handle : body_handle {
  asawa::rod::rod::ptr rod;

  explicit rod_body_handle(asawa::rod::rod::ptr r) : rod(std::move(r)) {}

  body_kind kind() const override { return body_kind::rod; }
  asawa::rod::rod &ref() const { return *rod; }
};

struct shell_body_handle : body_handle {
  asawa::shell::shell::ptr shell;

  explicit shell_body_handle(asawa::shell::shell::ptr m) : shell(std::move(m)) {}

  body_kind kind() const override { return body_kind::shell; }
  asawa::shell::shell &ref() const { return *shell; }
};

inline std::shared_ptr<body_handle>
make_rod_body(asawa::rod::rod::ptr rod) {
  return std::make_shared<rod_body_handle>(std::move(rod));
}

inline std::shared_ptr<body_handle>
make_shell_body(asawa::shell::shell::ptr shell) {
  return std::make_shared<shell_body_handle>(std::move(shell));
}

// Unbound position snapshot at the current body state (rod corner / shell vertex).
inline std::vector<vec3> snapshot_positions(const body_handle &body) {
  switch (body.kind()) {
  case body_kind::rod:
    return static_cast<const rod_body_handle &>(body).ref().xc();
  case body_kind::shell: {
    const auto &M = static_cast<const shell_body_handle &>(body).ref();
    return asawa::const_get_vec_data(M, 0);
  }
  }
  throw std::runtime_error("snapshot_positions: unknown body kind");
}

inline std::vector<vec3> snapshot_vertex_normals(const body_handle &body) {
  switch (body.kind()) {
  case body_kind::shell: {
    const auto &M = static_cast<const shell_body_handle &>(body).ref();
    const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);
    return asawa::shell::vertex_normals(M, x);
  }
  case body_kind::rod:
    throw std::runtime_error(
        "snapshot_vertex_normals: rod body not supported");
  }
  throw std::runtime_error("snapshot_vertex_normals: unknown body kind");
}

// Graph port datum: opaque body reference resolved at compute time.
struct body_datum : public liblombardi::Datum {
  std::shared_ptr<body_handle> handle;

  body_datum() = default;
  explicit body_datum(std::shared_ptr<body_handle> h) : handle(std::move(h)) {}

  void resize(size_t) override {}
  size_t size() const override { return handle ? 1u : 0u; }
  void clear() override { handle.reset(); }
  void *get_data() override { return handle.get(); }
  const void *get_data() const override { return handle.get(); }
};

inline asawa::rod::rod &require_rod(const body_handle &body) {
  if (body.kind() != body_kind::rod) {
    throw std::runtime_error("require_rod: expected rod body");
  }
  return static_cast<const rod_body_handle &>(body).ref();
}

inline asawa::shell::shell &require_shell(const body_handle &body) {
  if (body.kind() != body_kind::shell) {
    throw std::runtime_error("require_shell: expected shell body");
  }
  return static_cast<const shell_body_handle &>(body).ref();
}

} // namespace duchamp
} // namespace gaudi

#endif // __GAUDI_DUCHAMP_BODY_DATUM__
