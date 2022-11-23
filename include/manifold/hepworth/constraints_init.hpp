
#ifndef __PENKO_CONSTRAINTS_INIT__
#define __PENKO_CONSTRAINTS_INIT__

#include <manifold/asawa/coordinate_interface.hpp>
#include <manifold/asawa/m2.hpp>

#include "constraints.hpp"
#include "objective_function.hpp"

template <typename SPACE>
void init_stretch_constraints(
    asawa::surf<SPACE> *surf,
    typename hepworth::constraint_set<SPACE>::ptr constraints,
    double w = 1e-3) {
  M2_TYPEDEFS;

  std::vector<edge_ptr> edges = surf->get_edges();

  for (auto e : edges) {
    vertex_ptr v0 = e->v1()->vertex();
    vertex_ptr v1 = e->v2()->vertex();
    vertex_ptr v2 = e->v1()->prev()->vertex();
    vertex_ptr v3 = e->v2()->prev()->vertex();

    size_t i0 = v0->position_in_set();
    size_t i1 = v1->position_in_set();
    size_t i2 = v2->position_in_set();
    size_t i3 = v3->position_in_set();
#if 1
    typename hepworth::edge_stretch<SPACE>::ptr stretch01 =
        hepworth::edge_stretch<SPACE>::create(i0, i1);
    stretch01->set_weight(w);

    constraints->add_constraint(stretch01);
#endif
  }
}

template <typename SPACE>
void init_cross_constraints(
    asawa::surf<SPACE> *surf,
    typename hepworth::constraint_set<SPACE>::ptr constraints, const double &w,
    const double &k) {
  M2_TYPEDEFS;

  std::vector<edge_ptr> edges = surf->get_edges();

  for (auto e : edges) {
    vertex_ptr v0 = e->v1()->vertex();
    vertex_ptr v1 = e->v2()->vertex();
    vertex_ptr v2 = e->v1()->prev()->vertex();
    vertex_ptr v3 = e->v2()->prev()->vertex();

    size_t i0 = v0->position_in_set();
    size_t i1 = v1->position_in_set();
    size_t i2 = v2->position_in_set();
    size_t i3 = v3->position_in_set();

#if 0
      typename hepworth::edge_length<SPACE>::ptr length =
          hepworth::edge_length<SPACE>::create(i0, i1);
      // length->set_weight(1e-4 * ci::cotan<SPACE>(e));
      constraints->add_constraint(length);

#endif
#if 1
    typename hepworth::cross_ratio<SPACE>::ptr cross =
        hepworth::cross_ratio<SPACE>::create(i0, i1, i2, i3, k);
    // cross->set_weight(1e-6 * ci::cotan<SPACE>(e));
    cross->set_weight(w);

    constraints->add_constraint(cross);
#endif
#if 0
      typename hepworth::cross_ratio<SPACE>::ptr cross =
          hepworth::cross_ratio<SPACE>::create(i0, i1, i2, i3);
      cross->set_weight(5e-5);
      constraints->add_constraint(cross);
#endif
  }
}

template <typename SPACE>
void init_bend_constraints(
    asawa::surf<SPACE> *surf,
    typename hepworth::constraint_set<SPACE>::ptr constraints,
    double w = 1e-3) {
  M2_TYPEDEFS;

  std::vector<edge_ptr> edges = surf->get_edges();

  for (auto e : edges) {
    vertex_ptr v1 = e->v1()->vertex();
    vertex_ptr v2 = e->v2()->vertex();
    vertex_ptr v0 = e->v1()->prev()->vertex();
    vertex_ptr v3 = e->v2()->prev()->vertex();

    size_t i0 = v0->position_in_set();
    size_t i1 = v1->position_in_set();
    size_t i2 = v2->position_in_set();
    size_t i3 = v3->position_in_set();
#if 1
    typename hepworth::edge_bend<SPACE>::ptr bend =
        hepworth::edge_bend<SPACE>::create(i0, i1, i2, i3);
    bend->set_weight(w);
    constraints->add_constraint(bend);
#endif
  }
}

template <typename SPACE>
void init_mem_bend_constraints(
    asawa::surf<SPACE> *surf,
    typename hepworth::constraint_set<SPACE>::ptr constraints, double w0 = 1e-3,
    double w1 = 1e-3) {
  M2_TYPEDEFS;

  std::vector<edge_ptr> edges = surf->get_edges();

  for (auto e : edges) {
    vertex_ptr v1 = e->v1()->vertex();
    vertex_ptr v2 = e->v2()->vertex();
    vertex_ptr v0 = e->v1()->prev()->vertex();
    vertex_ptr v3 = e->v2()->prev()->vertex();

    size_t i0 = v0->position_in_set();
    size_t i1 = v1->position_in_set();
    size_t i2 = v2->position_in_set();
    size_t i3 = v3->position_in_set();

    coordinate_type c1 = ci::get_coordinate<SPACE>(e->v1());
    coordinate_type c2 = ci::get_coordinate<SPACE>(e->v2());

    coordinate_type de = (c2 - c1).normalized();
    coordinate_type N0 = ci::normal<SPACE>(e->v1()->face());
    coordinate_type N1 = ci::normal<SPACE>(e->v2()->face());

    real phi0 = e->template get<real>(SPACE::edge_index::PHI0);
    real phi1 = hepworth::tan_psi<SPACE>::static_calc_phi(N0, N1, de);
    if (phi0 < -9998) {
      phi0 = hepworth::tan_psi<SPACE>::static_calc_phi(N0, N1, de);
    }
    real k = 1.0;
    phi0 = (1.0 - k) * phi0 + k * phi1;
    e->template set<real>(SPACE::edge_index::PHI0, phi0);

#if 1
    typename hepworth::edge_mem_bend<SPACE>::ptr bend =
        hepworth::edge_mem_bend<SPACE>::create(i0, i1, i2, i3, phi0);
    bend->set_weight(w0, w1);
    constraints->add_constraint(bend);
#endif
  }
}

template <typename SPACE>
void init_willmore_constraints(
    asawa::surf<SPACE> *surf,
    typename hepworth::constraint_set<SPACE>::ptr constraints,
    double w = 1e-7) {
  M2_TYPEDEFS;

  std::vector<edge_ptr> edges = surf->get_edges();

  for (auto e : edges) {
    vertex_ptr v1 = e->v1()->vertex();
    vertex_ptr v2 = e->v2()->vertex();
    vertex_ptr v0 = e->v1()->prev()->vertex();
    vertex_ptr v3 = e->v2()->prev()->vertex();

    size_t i0 = v0->position_in_set();
    size_t i1 = v1->position_in_set();
    size_t i2 = v2->position_in_set();
    size_t i3 = v3->position_in_set();
#if 1
    typename hepworth::willmore_bend<SPACE>::ptr bend =
        hepworth::willmore_bend<SPACE>::create(i0, i1, i2, i3);
    bend->set_weight(w);
    constraints->add_constraint(bend);
#endif
  }
}
#endif