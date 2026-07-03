#ifndef __GAUDI_DUCHAMP_ROD_FORCE_NODES__
#define __GAUDI_DUCHAMP_ROD_FORCE_NODES__

#include <memory>

#include "gaudi/duchamp/field_nodes.hpp"
#include "gaudi/duchamp/rod_constraints_dynamics.hpp"
#include "liblombardi/junction_node.hpp"

namespace gaudi {
namespace duchamp {

class boundary_gradient_node : public liblombardi::Node {
public:
  using ptr = std::shared_ptr<boundary_gradient_node>;

  enum class PortId { Output = 0 };
  using OutputPortDef = liblombardi::PortDef<field_datum<vec3>, PortId::Output>;

  boundary_gradient_node(asawa::rod::rod::ptr rod, sdf_base::ptr sdf0, sdf_base::ptr sdf1)
      : _rod(std::move(rod)), _sdf0(std::move(sdf0)), _sdf1(std::move(sdf1)) {}

  void set_frame(int frame) { _frame = frame; }

  void compute() override {
    auto out = get_datum<OutputPortDef>();
    const auto sdf = select_rod_sdf(_frame, _sdf0, _sdf1);
    out->data() = compute_boundary_gradients(*_rod, *sdf);
  }

  unsigned int port_count() const override { return 1; }

  liblombardi::PortRef<boundary_gradient_node, OutputPortDef> output() { return {*this}; }

private:
  asawa::rod::rod::ptr _rod;
  sdf_base::ptr _sdf0;
  sdf_base::ptr _sdf1;
  int _frame = 0;
};

class tangent_point_gradient_node : public liblombardi::Node {
public:
  using ptr = std::shared_ptr<tangent_point_gradient_node>;

  enum class PortId { Output = 0 };
  using OutputPortDef = liblombardi::PortDef<field_datum<vec3>, PortId::Output>;

  tangent_point_gradient_node(asawa::rod::rod::ptr rod, asawa::rod::dynamic::ptr dynamic,
                              real scale = 1.0)
      : _rod(std::move(rod)), _dynamic(std::move(dynamic)), _scale(scale) {}

  void compute() override {
    auto out = get_datum<OutputPortDef>();
    auto gradient = compute_tangent_point_gradient(*_rod, *_dynamic);
    if (_scale != 1.0) {
      for (auto &g : gradient) {
        g *= _scale;
      }
    }
    out->data() = std::move(gradient);
  }

  unsigned int port_count() const override { return 1; }

  liblombardi::PortRef<tangent_point_gradient_node, OutputPortDef> output() { return {*this}; }

private:
  asawa::rod::rod::ptr _rod;
  asawa::rod::dynamic::ptr _dynamic;
  real _scale = 1.0;
};

class vortex_force_node : public liblombardi::Node {
public:
  using ptr = std::shared_ptr<vortex_force_node>;

  enum class PortId { Output = 0 };
  using OutputPortDef = liblombardi::PortDef<field_datum<vec3>, PortId::Output>;

  vortex_force_node(asawa::rod::rod::ptr rod, asawa::rod::dynamic::ptr dynamic,
                    real scale = 1.0, real p = 4.0, real q = 1.0)
      : _rod(std::move(rod)), _dynamic(std::move(dynamic)), _scale(scale), _p(p), _q(q) {}

  void compute() override {
    auto out = get_datum<OutputPortDef>();
    auto forces = compute_vortex_force(*_rod, *_dynamic, _p, _q);
    if (_scale != 1.0) {
      for (auto &f : forces) {
        f *= _scale;
      }
    }
    out->data() = std::move(forces);
  }

  unsigned int port_count() const override { return 1; }

  liblombardi::PortRef<vortex_force_node, OutputPortDef> output() { return {*this}; }

private:
  asawa::rod::rod::ptr _rod;
  asawa::rod::dynamic::ptr _dynamic;
  real _scale = 1.0;
  real _p = 4.0;
  real _q = 1.0;
};

template <int N>
using vec3_junction_node =
    liblombardi::junction_node<N, field_datum<vec3>, liblombardi::add_op<vec3>>;

} // namespace duchamp
} // namespace gaudi

#endif
