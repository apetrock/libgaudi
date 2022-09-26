
#ifndef __PENKO_CONSTRAINTS_INIT__
#define __PENKO_CONSTRAINTS_INIT__

#include <manifold/coordinate_interface.hpp>
#include <manifold/m2.hpp>

#include "constraints.hpp"
#include "objective_function.hpp"

template <typename SPACE>
void init_stretch_constraints(
    m2::surf<SPACE> *surf,
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
    m2::surf<SPACE> *surf,
    typename hepworth::constraint_set<SPACE>::ptr constraints) {
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
        hepworth::cross_ratio<SPACE>::create(i0, i1, i2, i3, 1.1);
    cross->set_weight(1e-6 * ci::cotan<SPACE>(e));

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
    m2::surf<SPACE> *surf,
    typename hepworth::constraint_set<SPACE>::ptr constraints,
    double w = 1e-3) {
  M2_TYPEDEFS;

  std::vector<edge_ptr> edges = surf->get_edges();
  std::vector<size_t> indices;
  std::vector<real> lengths;

  for (auto e : edges) {
    vertex_ptr v1 = e->v1()->vertex();
    vertex_ptr v2 = e->v2()->vertex();
    vertex_ptr v0 = e->v1()->prev()->vertex();
    vertex_ptr v3 = e->v2()->prev()->vertex();

    size_t i0 = v0->position_in_set();
    size_t i1 = v1->position_in_set();
    size_t i2 = v2->position_in_set();
    size_t i3 = v3->position_in_set();

    indices.push_back(i0);
    indices.push_back(i1);
    indices.push_back(i2);
    indices.push_back(i3);

    coordinate_type c0 = m2::ci::get_coordinate<SPACE>(v0);
    coordinate_type c1 = m2::ci::get_coordinate<SPACE>(v1);
    coordinate_type c2 = m2::ci::get_coordinate<SPACE>(v2);
    coordinate_type c3 = m2::ci::get_coordinate<SPACE>(v3);
    real l01 = m2::va::norm(coordinate_type(c0 - c1));
    real l23 = m2::va::norm(coordinate_type(c3 - c2));
    lengths.push_back(l01);
    lengths.push_back(l23);
  }

  int i = 0;
  int N = 0.5 * lengths.size();

  for (i = 0; i < N; i++) {
    size_t i0 = indices[4 * i + 0];
    size_t i1 = indices[4 * i + 1];
    size_t i2 = indices[4 * i + 2];
    size_t i3 = indices[4 * i + 3];

    real l01 = lengths[2 * i + 0];
    real l23 = lengths[2 * i + 1];

#if 1
    typename hepworth::edge_bend<SPACE>::ptr bend =
        hepworth::edge_bend<SPACE>::create(i0, i1, i2, i3);
    bend->set_weight(w);
    constraints->add_constraint(bend);
#endif
  }
}

#endif