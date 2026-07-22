#ifndef __GAUDI_DUCHAMP_FIELD_GRAPH_NODES__
#define __GAUDI_DUCHAMP_FIELD_GRAPH_NODES__

#include <memory>
#include <stdexcept>
#include <vector>

#include "gaudi/calder/rod_integrators.hpp"
#include "gaudi/calder/shell_integrators.hpp"
#include "gaudi/duchamp/body_datum.hpp"
#include "gaudi/duchamp/field_nodes.hpp"

#include "liblombardi/node_base.hpp"
#include "liblombardi/port_def.hpp"

namespace gaudi {
namespace duchamp {

// Emits a fixed body handle each compute (for wiring body into downstream nodes).
class body_constant_node : public liblombardi::Node {
public:
  using ptr = std::shared_ptr<body_constant_node>;

  enum class PortId { Output = 0 };
  using OutputPortDef = liblombardi::PortDef<body_datum, PortId::Output>;

  explicit body_constant_node(std::shared_ptr<body_handle> body)
      : _body(std::move(body)) {}

  void compute() override { get_datum<OutputPortDef>()->handle = _body; }

  unsigned int port_count() const override { return 1; }

  liblombardi::PortRef<body_constant_node, OutputPortDef> output() {
    return {*this};
  }

private:
  std::shared_ptr<body_handle> _body;
};

// Copies body positions into an unbound field_datum<vec3> at this pipeline stage.
class position_snapshot_node : public liblombardi::Node {
public:
  using ptr = std::shared_ptr<position_snapshot_node>;

  enum class PortId { Body = 0, Output = 1 };
  using BodyPortDef = liblombardi::PortDef<body_datum, PortId::Body>;
  using OutputPortDef = liblombardi::PortDef<field_datum<vec3>, PortId::Output>;

  void compute() override {
    const auto body = get_datum<BodyPortDef>();
    auto out = get_datum<OutputPortDef>();
    if (!body->handle) {
      throw std::runtime_error("position_snapshot_node: missing body handle");
    }
    out->data() = snapshot_positions(*body->handle);
  }

  unsigned int port_count() const override { return 2; }

  liblombardi::PortRef<position_snapshot_node, BodyPortDef> body() {
    return {*this};
  }
  liblombardi::PortRef<position_snapshot_node, OutputPortDef> output() {
    return {*this};
  }
};

// Copies shell vertex normals into an unbound field_datum<vec3>.
class vertex_normals_snapshot_node : public liblombardi::Node {
public:
  using ptr = std::shared_ptr<vertex_normals_snapshot_node>;

  enum class PortId { Body = 0, Output = 1 };
  using BodyPortDef = liblombardi::PortDef<body_datum, PortId::Body>;
  using OutputPortDef = liblombardi::PortDef<field_datum<vec3>, PortId::Output>;

  void compute() override {
    const auto body = get_datum<BodyPortDef>();
    auto out = get_datum<OutputPortDef>();
    if (!body->handle) {
      throw std::runtime_error(
          "vertex_normals_snapshot_node: missing body handle");
    }
    out->data() = snapshot_vertex_normals(*body->handle);
  }

  unsigned int port_count() const override { return 2; }

  liblombardi::PortRef<vertex_normals_snapshot_node, BodyPortDef> body() {
    return {*this};
  }
  liblombardi::PortRef<vertex_normals_snapshot_node, OutputPortDef> output() {
    return {*this};
  }
};

// Emits a fixed unbound POV list (field_datum<vec3>).
class pov_constant_node : public liblombardi::Node {
public:
  using ptr = std::shared_ptr<pov_constant_node>;

  enum class PortId { Output = 0 };
  using OutputPortDef = liblombardi::PortDef<field_datum<vec3>, PortId::Output>;

  explicit pov_constant_node(std::vector<vec3> pov) : _pov(std::move(pov)) {}

  void compute() override { get_datum<OutputPortDef>()->data() = _pov; }

  unsigned int port_count() const override { return 1; }

  liblombardi::PortRef<pov_constant_node, OutputPortDef> output() {
    return {*this};
  }

private:
  std::vector<vec3> _pov;
};

struct rod_mls_avg_tag {
  using value_type = vec3;
};

struct shell_mls_avg_tag {
  using value_type = real;
};

namespace detail {

template <typename T>
inline std::vector<T>
dispatch_mls_avg(rod_mls_avg_tag, const body_handle &body,
                 const std::vector<T> &source_field,
                 const std::vector<vec3> &p_pov, real l0, real p) {
  return calder::mls_avg(require_rod(body), source_field, p_pov, l0, p);
}

template <typename T>
inline std::vector<T>
dispatch_mls_avg(shell_mls_avg_tag, const body_handle &body,
                 const std::vector<T> &source_field,
                 const std::vector<vec3> &p_pov, real l0, real p) {
  return calder::mls_avg(require_shell(body), source_field, p_pov, l0, p);
}

} // namespace detail

// Calder eval node: body + unbound POV in, sampled field_datum<T> out.
template <typename IntegratorTag, typename T = typename IntegratorTag::value_type>
class nbody_eval_node : public liblombardi::Node {
public:
  using ptr = std::shared_ptr<nbody_eval_node>;

  enum class PortId { Body = 0, Pov = 1, Output = 2 };
  using BodyPortDef = liblombardi::PortDef<body_datum, PortId::Body>;
  using PovPortDef = liblombardi::PortDef<field_datum<vec3>, PortId::Pov>;
  using OutputPortDef = liblombardi::PortDef<field_datum<T>, PortId::Output>;

  nbody_eval_node(std::vector<T> source_field, real l0, real p = 3.0)
      : _source_field(std::move(source_field)), _l0(l0), _p(p) {}

  void compute() override {
    const auto body = get_datum<BodyPortDef>();
    const auto pov = get_datum<PovPortDef>();
    auto out = get_datum<OutputPortDef>();
    if (!body->handle) {
      throw std::runtime_error("nbody_eval_node: missing body handle");
    }
    out->data() = detail::dispatch_mls_avg(IntegratorTag{}, *body->handle,
                                           _source_field, pov->data(), _l0, _p);
  }

  unsigned int port_count() const override { return 3; }

  liblombardi::PortRef<nbody_eval_node, BodyPortDef> body() { return {*this}; }
  liblombardi::PortRef<nbody_eval_node, PovPortDef> pov() { return {*this}; }
  liblombardi::PortRef<nbody_eval_node, OutputPortDef> output() {
    return {*this};
  }

private:
  std::vector<T> _source_field;
  real _l0 = 1.0;
  real _p = 3.0;
};

} // namespace duchamp
} // namespace gaudi

#endif // __GAUDI_DUCHAMP_FIELD_GRAPH_NODES__
