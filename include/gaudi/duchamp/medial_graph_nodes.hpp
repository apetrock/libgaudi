#ifndef __GAUDI_DUCHAMP_MEDIAL_GRAPH_NODES__
#define __GAUDI_DUCHAMP_MEDIAL_GRAPH_NODES__

#include <cmath>
#include <limits>
#include <vector>

#include "gaudi/albers/darboux_cyclide.hpp"
#include "gaudi/duchamp/field_graph_nodes.hpp"
#include "gaudi/duchamp/medial_field_datums.hpp"
#include "gaudi/duchamp/medial_search_helpers.hpp"

#include "liblombardi/node_base.hpp"
#include "liblombardi/port_def.hpp"

namespace gaudi {
namespace duchamp {

// Shared port layout for both medial search backends:
//   inputs:  pov, cyclide, n_pov
//   outputs: points (vec3), mask (real; 1=accepted, 0=rejected)
class medial_shape_energy_search_node : public liblombardi::Node {
public:
  using ptr = std::shared_ptr<medial_shape_energy_search_node>;

  enum class PortId { Pov = 0, Cyclide = 1, NPov = 2, Output = 3, Mask = 4 };
  using PovPortDef = liblombardi::PortDef<field_datum<vec3>, PortId::Pov>;
  using CyclidePortDef =
      liblombardi::PortDef<cyclide_field_datum, PortId::Cyclide>;
  using NPovPortDef = liblombardi::PortDef<field_datum<vec3>, PortId::NPov>;
  using OutputPortDef = liblombardi::PortDef<field_datum<vec3>, PortId::Output>;
  using MaskPortDef = liblombardi::PortDef<field_datum<real>, PortId::Mask>;

  explicit medial_shape_energy_search_node(
      real max_travel, albers::medial_shape_search_params search = {})
      : _max_travel(max_travel), _search(std::move(search)) {}

  void compute() override {
    const auto pov = get_datum<PovPortDef>();
    const auto cyclide = get_datum<CyclidePortDef>();
    const auto n_pov = get_datum<NPovPortDef>();
    auto out = get_datum<OutputPortDef>();
    auto mask = get_datum<MaskPortDef>();

    const size_t n = std::min(
        {pov->data().size(), cyclide->data().size(), n_pov->data().size()});
    out->data().resize(n, vec3::Zero());
    mask->data().resize(n, real(0.0));

    for (size_t i = 0; i < n; ++i) {
      const medial_point r = search_medial_shape_energy(
          cyclide->data()[i], pov->data()[i], n_pov->data()[i], _max_travel,
          _search);
      out->data()[i] = r.point;
      mask->data()[i] = r.accepted ? real(1.0) : real(0.0);
    }
  }

  unsigned int port_count() const override { return 5; }

  liblombardi::PortRef<medial_shape_energy_search_node, PovPortDef> pov() {
    return {*this};
  }
  liblombardi::PortRef<medial_shape_energy_search_node, CyclidePortDef>
  cyclide() {
    return {*this};
  }
  liblombardi::PortRef<medial_shape_energy_search_node, NPovPortDef> n_pov() {
    return {*this};
  }
  liblombardi::PortRef<medial_shape_energy_search_node, OutputPortDef>
  output() {
    return {*this};
  }
  liblombardi::PortRef<medial_shape_energy_search_node, MaskPortDef> mask() {
    return {*this};
  }

private:
  real _max_travel = 1.0;
  albers::medial_shape_search_params _search;
};

class medial_legacy_ridge_search_node : public liblombardi::Node {
public:
  using ptr = std::shared_ptr<medial_legacy_ridge_search_node>;

  enum class PortId { Pov = 0, Cyclide = 1, NPov = 2, Output = 3, Mask = 4 };
  using PovPortDef = liblombardi::PortDef<field_datum<vec3>, PortId::Pov>;
  using CyclidePortDef =
      liblombardi::PortDef<cyclide_field_datum, PortId::Cyclide>;
  using NPovPortDef = liblombardi::PortDef<field_datum<vec3>, PortId::NPov>;
  using OutputPortDef = liblombardi::PortDef<field_datum<vec3>, PortId::Output>;
  using MaskPortDef = liblombardi::PortDef<field_datum<real>, PortId::Mask>;

  medial_legacy_ridge_search_node(
      real max_travel, int max_iters = 100, real tol = 1e-8,
      real max_newton_step = std::numeric_limits<real>::infinity(),
      real min_normal_alignment = -0.25)
      : _max_travel(max_travel), _max_iters(max_iters), _tol(tol),
        _max_newton_step(max_newton_step),
        _min_normal_alignment(min_normal_alignment) {}

  void compute() override {
    const auto pov = get_datum<PovPortDef>();
    const auto cyclide = get_datum<CyclidePortDef>();
    const auto n_pov = get_datum<NPovPortDef>();
    auto out = get_datum<OutputPortDef>();
    auto mask = get_datum<MaskPortDef>();

    const size_t n = std::min(
        {pov->data().size(), cyclide->data().size(), n_pov->data().size()});
    out->data().resize(n, vec3::Zero());
    mask->data().resize(n, real(0.0));

    for (size_t i = 0; i < n; ++i) {
      const medial_point r = search_medial_legacy_ridge(
          cyclide->data()[i], pov->data()[i], n_pov->data()[i], _max_travel,
          _max_iters, _tol, _max_newton_step, _min_normal_alignment);
      out->data()[i] = r.point;
      mask->data()[i] = r.accepted ? real(1.0) : real(0.0);
    }
  }

  unsigned int port_count() const override { return 5; }

  liblombardi::PortRef<medial_legacy_ridge_search_node, PovPortDef> pov() {
    return {*this};
  }
  liblombardi::PortRef<medial_legacy_ridge_search_node, CyclidePortDef>
  cyclide() {
    return {*this};
  }
  liblombardi::PortRef<medial_legacy_ridge_search_node, NPovPortDef> n_pov() {
    return {*this};
  }
  liblombardi::PortRef<medial_legacy_ridge_search_node, OutputPortDef>
  output() {
    return {*this};
  }
  liblombardi::PortRef<medial_legacy_ridge_search_node, MaskPortDef> mask() {
    return {*this};
  }

private:
  real _max_travel = 1.0;
  int _max_iters = 100;
  real _tol = 1e-8;
  real _max_newton_step = std::numeric_limits<real>::infinity();
  real _min_normal_alignment = -0.25;
};

class hessian_frame_node : public liblombardi::Node {
public:
  using ptr = std::shared_ptr<hessian_frame_node>;

  enum class PortId { Cyclide = 0, Pov = 1, Position = 2, Output = 3 };
  using CyclidePortDef =
      liblombardi::PortDef<cyclide_field_datum, PortId::Cyclide>;
  using PovPortDef = liblombardi::PortDef<field_datum<vec3>, PortId::Pov>;
  using PositionPortDef =
      liblombardi::PortDef<field_datum<vec3>, PortId::Position>;
  using OutputPortDef = liblombardi::PortDef<frame_field_datum, PortId::Output>;

  void compute() override {
    const auto cyclide = get_datum<CyclidePortDef>();
    const auto pov = get_datum<PovPortDef>();
    const auto position = get_datum<PositionPortDef>();
    auto out = get_datum<OutputPortDef>();

    const size_t n = std::min(
        {cyclide->data().size(), pov->data().size(), position->data().size()});
    out->data().resize(n, mat3::Identity());

    for (size_t i = 0; i < n; ++i) {
      const vec3 x_local = position->data()[i] - pov->data()[i];
      const mat3 W =
          albers::shape_operator_at(cyclide->data()[i], x_local);
      if (!W.allFinite()) {
        out->data()[i].setIdentity();
        continue;
      }
      Eigen::SelfAdjointEigenSolver<mat3> es(W);
      if (es.info() != Eigen::Success) {
        out->data()[i].setIdentity();
        continue;
      }
      out->data()[i] = es.eigenvectors();
    }
  }

  unsigned int port_count() const override { return 4; }

  liblombardi::PortRef<hessian_frame_node, CyclidePortDef> cyclide() {
    return {*this};
  }
  liblombardi::PortRef<hessian_frame_node, PovPortDef> pov() { return {*this}; }
  liblombardi::PortRef<hessian_frame_node, PositionPortDef> position() {
    return {*this};
  }
  liblombardi::PortRef<hessian_frame_node, OutputPortDef> output() {
    return {*this};
  }
};

} // namespace duchamp
} // namespace gaudi

#endif // __GAUDI_DUCHAMP_MEDIAL_GRAPH_NODES__
