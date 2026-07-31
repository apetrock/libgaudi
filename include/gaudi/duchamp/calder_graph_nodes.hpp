#ifndef __GAUDI_DUCHAMP_CALDER_GRAPH_NODES__
#define __GAUDI_DUCHAMP_CALDER_GRAPH_NODES__

#include <stdexcept>
#include <utility>

#include "gaudi/albers/line_cylinder.hpp"
#include "gaudi/calder/least_squares_fit.hpp"
#include "gaudi/calder/mls_jet_bootstrap.hpp"
#include "gaudi/duchamp/body_datum.hpp"
#include "gaudi/duchamp/field_graph_nodes.hpp"
#include "gaudi/duchamp/medial_field_datums.hpp"
#include "gaudi/kusama/cyclide_jet_smooth.hpp"

#include "liblombardi/node_base.hpp"
#include "liblombardi/port_def.hpp"

namespace gaudi {
namespace duchamp {

class darboux_cyclide_fit_node : public liblombardi::Node {
public:
  using ptr = std::shared_ptr<darboux_cyclide_fit_node>;

  enum class PortId { Body = 0, Pov = 1, NPov = 2, Output = 3 };
  using BodyPortDef = liblombardi::PortDef<body_datum, PortId::Body>;
  using PovPortDef = liblombardi::PortDef<field_datum<vec3>, PortId::Pov>;
  using NPovPortDef = liblombardi::PortDef<field_datum<vec3>, PortId::NPov>;
  using OutputPortDef = liblombardi::PortDef<cyclide_field_datum, PortId::Output>;

  darboux_cyclide_fit_node(real l0, real fit_p = 3.0, real fit_w0 = 1.0,
                           real radius_scale = 1.0,
                           calder::stage1_radius_model stage1 =
                               calder::stage1_radius_model::cyclide,
                           bool use_mls_bootstrap = false,
                           calder::mls_jet_bootstrap_params bootstrap = {})
      : _l0(l0), _fit_p(fit_p), _fit_w0(fit_w0), _radius_scale(radius_scale),
        _stage1(stage1), _use_mls_bootstrap(use_mls_bootstrap),
        _bootstrap(std::move(bootstrap)) {}

  void compute() override {
    const auto body = get_datum<BodyPortDef>();
    const auto pov = get_datum<PovPortDef>();
    const auto n_pov = get_datum<NPovPortDef>();
    auto out = get_datum<OutputPortDef>();
    if (!body->handle) {
      throw std::runtime_error("darboux_cyclide_fit_node: missing body");
    }
    if (body->handle->kind() != body_kind::shell) {
      throw std::runtime_error("darboux_cyclide_fit_node: shell body required");
    }
    auto &M = static_cast<shell_body_handle &>(*body->handle).ref();
    if (_use_mls_bootstrap) {
      calder::mls_jet_bootstrap_params bp = _bootstrap;
      bp.fit_p = _fit_p;
      bp.fit_w0 = _fit_w0;
      bp.radius_scale = _radius_scale;
      out->data() = calder::darboux_fit_bootstrapped(M, pov->data(),
                                                     n_pov->data(), _l0, bp);
    } else {
      out->data() = calder::darboux_cyclide_shell_fit(
          M, pov->data(), n_pov->data(), _l0, _fit_p, _fit_w0, _radius_scale,
          _stage1);
    }
  }

  unsigned int port_count() const override { return 4; }

  liblombardi::PortRef<darboux_cyclide_fit_node, BodyPortDef> body() {
    return {*this};
  }
  liblombardi::PortRef<darboux_cyclide_fit_node, PovPortDef> pov() {
    return {*this};
  }
  liblombardi::PortRef<darboux_cyclide_fit_node, NPovPortDef> n_pov() {
    return {*this};
  }
  liblombardi::PortRef<darboux_cyclide_fit_node, OutputPortDef> output() {
    return {*this};
  }

private:
  real _l0 = 1.0;
  real _fit_p = 3.0;
  real _fit_w0 = 1.0; // w_foot = _fit_w0 * Σ w_MLS
  real _radius_scale = 1.0;
  calder::stage1_radius_model _stage1 = calder::stage1_radius_model::cyclide;
  bool _use_mls_bootstrap = false;
  calder::mls_jet_bootstrap_params _bootstrap;
};

class darboux_cyclide_smooth_node : public liblombardi::Node {
public:
  using ptr = std::shared_ptr<darboux_cyclide_smooth_node>;

  enum class PortId { Body = 0, Pov = 1, CyclideIn = 2, CyclideOut = 3 };
  using BodyPortDef = liblombardi::PortDef<body_datum, PortId::Body>;
  using PovPortDef = liblombardi::PortDef<field_datum<vec3>, PortId::Pov>;
  using CyclideInPortDef =
      liblombardi::PortDef<cyclide_field_datum, PortId::CyclideIn>;
  using CyclideOutPortDef =
      liblombardi::PortDef<cyclide_field_datum, PortId::CyclideOut>;

  explicit darboux_cyclide_smooth_node(
      kusama::cyclide_jet_smooth_params params = {})
      : _params(params) {}

  void compute() override {
    const auto body = get_datum<BodyPortDef>();
    const auto pov = get_datum<PovPortDef>();
    const auto cyclide_in = get_datum<CyclideInPortDef>();
    auto out = get_datum<CyclideOutPortDef>();
    if (!body->handle) {
      throw std::runtime_error("darboux_cyclide_smooth_node: missing body");
    }
    if (body->handle->kind() != body_kind::shell) {
      throw std::runtime_error(
          "darboux_cyclide_smooth_node: shell body required");
    }
    auto &M = static_cast<shell_body_handle &>(*body->handle).ref();
    out->data() = kusama::cyclide_jet_smooth(M, pov->data(),
                                             cyclide_in->data(), _params);
  }

  unsigned int port_count() const override { return 4; }

  liblombardi::PortRef<darboux_cyclide_smooth_node, BodyPortDef> body() {
    return {*this};
  }
  liblombardi::PortRef<darboux_cyclide_smooth_node, PovPortDef> pov() {
    return {*this};
  }
  liblombardi::PortRef<darboux_cyclide_smooth_node, CyclideInPortDef>
  cyclide_in() {
    return {*this};
  }
  liblombardi::PortRef<darboux_cyclide_smooth_node, CyclideOutPortDef>
  cyclide_out() {
    return {*this};
  }

private:
  kusama::cyclide_jet_smooth_params _params;
};

class cylinder_fit_mesh_node : public liblombardi::Node {
public:
  using ptr = std::shared_ptr<cylinder_fit_mesh_node>;

  enum class PortId { Body = 0, Pov = 1, NPov = 2, Output = 3 };
  using BodyPortDef = liblombardi::PortDef<body_datum, PortId::Body>;
  using PovPortDef = liblombardi::PortDef<field_datum<vec3>, PortId::Pov>;
  using NPovPortDef = liblombardi::PortDef<field_datum<vec3>, PortId::NPov>;
  using OutputPortDef = liblombardi::PortDef<line_field_datum, PortId::Output>;

  cylinder_fit_mesh_node(real l0, real fit_p = 3.0)
      : _l0(l0), _fit_p(fit_p) {}

  void compute() override {
    const auto body = get_datum<BodyPortDef>();
    const auto pov = get_datum<PovPortDef>();
    const auto n_pov = get_datum<NPovPortDef>();
    auto out = get_datum<OutputPortDef>();
    if (!body->handle) {
      throw std::runtime_error("cylinder_fit_mesh_node: missing body");
    }
    if (body->handle->kind() != body_kind::shell) {
      throw std::runtime_error("cylinder_fit_mesh_node: shell body required");
    }
    auto &M = static_cast<shell_body_handle &>(*body->handle).ref();
    out->data() = calder::normal_aligned_line_convexity(
        M, pov->data(), n_pov->data(), _l0, _fit_p);
  }

  unsigned int port_count() const override { return 4; }

  liblombardi::PortRef<cylinder_fit_mesh_node, BodyPortDef> body() {
    return {*this};
  }
  liblombardi::PortRef<cylinder_fit_mesh_node, PovPortDef> pov() {
    return {*this};
  }
  liblombardi::PortRef<cylinder_fit_mesh_node, NPovPortDef> n_pov() {
    return {*this};
  }
  liblombardi::PortRef<cylinder_fit_mesh_node, OutputPortDef> output() {
    return {*this};
  }

private:
  real _l0 = 1.0;
  real _fit_p = 3.0;
};

class cylinder_fit_points_node : public liblombardi::Node {
public:
  using ptr = std::shared_ptr<cylinder_fit_points_node>;

  enum class PortId { Pov = 0, NPov = 1, Points = 2, Output = 3 };
  using PovPortDef = liblombardi::PortDef<field_datum<vec3>, PortId::Pov>;
  using NPovPortDef = liblombardi::PortDef<field_datum<vec3>, PortId::NPov>;
  using PointsPortDef = liblombardi::PortDef<field_datum<vec3>, PortId::Points>;
  using OutputPortDef = liblombardi::PortDef<line_field_datum, PortId::Output>;

  void compute() override {
    const auto pov = get_datum<PovPortDef>();
    const auto n_pov = get_datum<NPovPortDef>();
    const auto points = get_datum<PointsPortDef>();
    auto out = get_datum<OutputPortDef>();

    const size_t n =
        std::min({pov->data().size(), n_pov->data().size(), points->data().size()});
    out->data().resize(n);
    for (size_t i = 0; i < n; ++i) {
      out->data()[i] = albers::fit_normal_aligned_line(
          {points->data()[i]}, {n_pov->data()[i]});
    }
  }

  unsigned int port_count() const override { return 4; }

  liblombardi::PortRef<cylinder_fit_points_node, PovPortDef> pov() {
    return {*this};
  }
  liblombardi::PortRef<cylinder_fit_points_node, NPovPortDef> n_pov() {
    return {*this};
  }
  liblombardi::PortRef<cylinder_fit_points_node, PointsPortDef> points() {
    return {*this};
  }
  liblombardi::PortRef<cylinder_fit_points_node, OutputPortDef> output() {
    return {*this};
  }
};

} // namespace duchamp
} // namespace gaudi

#endif // __GAUDI_DUCHAMP_CALDER_GRAPH_NODES__
