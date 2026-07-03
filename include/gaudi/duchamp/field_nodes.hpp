#ifndef __GAUDI_DUCHAMP_FIELD_NODES__
#define __GAUDI_DUCHAMP_FIELD_NODES__

#include <memory>
#include <vector>

#include "gaudi/common.h"
#include "gaudi/duchamp/fields.hpp"

#include "liblombardi/datum_pool.hpp"
#include "liblombardi/graph_context.hpp"
#include "liblombardi/node_base.hpp"
#include "liblombardi/port_def.hpp"

namespace gaudi {
namespace duchamp {

// Bridge datum: a liblombardi pool datum that holds a std::vector<T>, so field
// data can flow through the liblombardi graph. Mirrors DatumImpl but is its own
// type so PortDef<field_datum<T>, ...> stays distinct from the pool's defaults.
template <typename T>
struct field_datum : public liblombardi::Datum {
  using value_type = T;
  std::vector<T> _data;

  field_datum() = default;
  explicit field_datum(const std::vector<T> &data) : _data(data) {}

  void resize(size_t n) override { _data.resize(n); }
  size_t size() const override { return _data.size(); }
  void clear() override { _data.clear(); }
  void *get_data() override { return _data.data(); }
  const void *get_data() const override { return _data.data(); }

  std::vector<T> &data() { return _data; }
  const std::vector<T> &data() const { return _data; }
  std::vector<T> &values() { return _data; }
  const std::vector<T> &values() const { return _data; }
};

// Single explicit march step: out[i] = in[i] + step * normal[i].
//
// ONE input port (positions) and ONE output port (marched positions). The
// normal field is a constructor parameter passed polymorphically as
// field_base<...>::ptr, so the node both exercises the field system and stays a
// clean 1-in / 1-out node for the graph.
class march_node : public liblombardi::Node {
public:
  using ptr = std::shared_ptr<march_node>;
  using normal_field =
      field_base<asawa::shell::shell, vec3, asawa::prim_type::VERTEX>;

  enum class PortId { Input = 0, Output = 1 };

  using InputPortDef = liblombardi::PortDef<field_datum<vec3>, PortId::Input>;
  using OutputPortDef = liblombardi::PortDef<field_datum<vec3>, PortId::Output>;

  march_node() = default;
  march_node(normal_field::ptr normals, real step)
      : _normals(std::move(normals)), _step(step) {}

  void set_normals(normal_field::ptr normals) { _normals = std::move(normals); }
  void set_step(real step) { _step = step; }

  void compute() override {
    auto in = get_datum<InputPortDef>();
    auto out = get_datum<OutputPortDef>();
    const std::vector<vec3> &x = in->data();
    std::vector<vec3> &y = out->data();
    y.resize(x.size());

    if (!_normals) {
      y = x;
      return;
    }
    const std::vector<vec3> &n = _normals->get();
    for (size_t i = 0; i < x.size(); ++i) {
      const vec3 ni = i < n.size() ? n[i] : vec3::Zero();
      y[i] = x[i] + _step * ni;
    }
  }

  unsigned int port_count() const override { return 2; }

  liblombardi::PortRef<march_node, InputPortDef> input() { return {*this}; }
  liblombardi::PortRef<march_node, OutputPortDef> output() { return {*this}; }

private:
  normal_field::ptr _normals;
  real _step = 0.0;
};

} // namespace duchamp
} // namespace gaudi

#endif // __GAUDI_DUCHAMP_FIELD_NODES__
