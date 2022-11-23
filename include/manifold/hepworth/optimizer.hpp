
#ifndef __PENKO_OPTIMIZER__
#define __PENKO_OPTIMIZER__

#include <bits/types/wint_t.h>
#include <vector>

#include <manifold/asawa/m2.hpp>
#include <manifold/calder/harmonic_integrators.hpp>

#include "solver.hpp"

namespace hepworth {
template <typename SPACE> class optimizer {

  M2_TYPEDEFS;

public:
  typedef Eigen::Matrix<real, 2, 1> vec2;
  typedef Eigen::Matrix<real, 6, 1> vec6;
  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 3, 2> mat32;
  typedef Eigen::Matrix<real, 6, 6> mat66;
  typedef Eigen::Matrix<real, Eigen::Dynamic, 1> vecX;
  typedef Eigen::SparseMatrix<real> sparmat;

  optimizer() {}

  void test_rosenbrock() {
    typename rosenbrock_function<SPACE>::ptr fcn =
        rosenbrock_function<SPACE>::create();
    std::cout << "--- rosenbrock test ---" << std::endl;
    vecX x(2);
    x[0] = 12.0;
    x[1] = 0.2125;

    fcn->set_x(x);
    // gradient_descent_solver<SPACE> solver;
    newton_raphson_solver<SPACE> solver;
    solver.solve(fcn);
    // std::cout << fcn->get_x() << std::endl;
    // std::cout << "--- --------------- ---" << std::endl;
  }

  void update(typename objective_function<SPACE>::ptr objective_function) {
    test_rosenbrock();
#if 0
    newton_raphson_solver<SPACE> solver;
#else
    gradient_descent_solver<SPACE> solver;
#endif
    solver.solve(objective_function);
  }

}; // class optimizer

template <typename SPACE> class position_optimizer {

  M2_TYPEDEFS;

public:
  typedef Eigen::Matrix<real, 2, 1> vec2;
  typedef Eigen::Matrix<real, 6, 1> vec6;
  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 3, 2> mat32;
  typedef Eigen::Matrix<real, 6, 6> mat66;
  typedef Eigen::Matrix<real, Eigen::Dynamic, 1> vecX;
  typedef Eigen::SparseMatrix<real> sparmat;

  position_optimizer(const std::vector<vec3> &p0, const std::vector<vec3> &p1,
                     const std::vector<real> &weights)
      : _p0(p0), _p1(p1), _weights(weights) {}

#if 0

  virtual vecX project_gradient(const vecX &g) {
    coordinate_array positions = get_positions();
    coordinate_array gc(positions);
    vec_interface<SPACE>::from(gc, g);
    coordinate_array gf =
        asawa::ci::verts_to_faces<SPACE, coordinate_type>(gc, _surf);
    asawa::mesh_calculator<SPACE> calc;
    // this needs to be reworked to take a surface topology and a position set
    coordinate_array gn = calc.harmonicAvg(_surf, gf, positions, 1.0 * _reg);

    return vec_interface<SPACE>::to(gn);
  };
#else
  coordinate_array project(const coordinate_array &p0,
                           const coordinate_array &p1,
                           typename constraint_set<SPACE>::ptr constraint_set) {
    surf_ptr surf = constraint_set->_surf;
    real reg = 3.0 * asawa::ci::geometric_mean_length<SPACE>(surf);

    coordinate_array pp(p0);
    coordinate_array dp(p0);

    for (int i = 0; i < p0.size(); i++) {
      dp[i] = p1[i] - p0[i];
    }

    coordinate_array dpf =
        asawa::ci::verts_to_faces<SPACE, coordinate_type>(dp, surf);

    calder::mesh_calculator<SPACE> calc;
    // this needs to be reworked to take a surface topology and a position set
    coordinate_array dpn = calc.harmonicAvg(surf, dpf, p0, 1.0 * reg);

    for (int i = 0; i < p0.size(); i++) {
      pp[i] += dpn[i];
    }
    return pp;
  };
#endif

  void update(typename constraint_set<SPACE>::ptr constraint_set) {
    const std::vector<vec3> &positions = constraint_set->get_positions();
    // do stuff with the positions.
    std::vector<typename point<SPACE>::ptr> position_constraints;

    for (int i = 0; i < _p1.size(); i++) {
      vec3 target = _p1[i];
      real weight = _weights[i];
      typename point<SPACE>::ptr pos = point<SPACE>::create(i, target, weight);
      position_constraints.push_back(pos);
      constraint_set->add_constraint(pos);
    }

    constraint_set->set_positions(_p0);

#if 0
    newton_raphson_solver<SPACE> solver;
#else
    gradient_descent_solver<SPACE> solver;
#endif
    solver.solve(constraint_set);

    coordinate_array pp = project(_p0, _p1, constraint_set);
    constraint_set->set_positions(pp);
  }
  std::vector<real> _weights;
  std::vector<vec3> _p1;
  std::vector<vec3> _p0;
}; // class position_optimizer

template <typename SPACE> class velocity_optimizer {

  M2_TYPEDEFS;

public:
  typedef Eigen::Matrix<real, 2, 1> vec2;
  typedef Eigen::Matrix<real, 6, 1> vec6;
  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 3, 2> mat32;
  typedef Eigen::Matrix<real, 6, 6> mat66;
  typedef Eigen::Matrix<real, Eigen::Dynamic, 1> vecX;
  typedef Eigen::SparseMatrix<real> sparmat;

  velocity_optimizer(const std::vector<vec3> &p0, const std::vector<vec3> &p1)
      : _p0(p0), _p1(p1) {}

#if 0

  virtual vecX project_gradient(const vecX &g) {
    coordinate_array positions = get_positions();
    coordinate_array gc(positions);
    vec_interface<SPACE>::from(gc, g);
    coordinate_array gf =
        asawa::ci::verts_to_faces<SPACE, coordinate_type>(gc, _surf);
    asawa::mesh_calculator<SPACE> calc;
    // this needs to be reworked to take a surface topology and a position set
    coordinate_array gn = calc.harmonicAvg(_surf, gf, positions, 1.0 * _reg);

    return vec_interface<SPACE>::to(gn);
  };
#else
  coordinate_array project(const coordinate_array &p0,
                           const coordinate_array &p1,
                           typename constraint_set<SPACE>::ptr constraint_set) {
    surf_ptr surf = constraint_set->_surf;
    real reg = 3.0 * asawa::ci::geometric_mean_length<SPACE>(surf);

    coordinate_array pp(p0);
    coordinate_array dp(p0);

    for (int i = 0; i < p0.size(); i++) {
      dp[i] = p1[i] - p0[i];
    }

    coordinate_array dpf =
        asawa::ci::verts_to_faces<SPACE, coordinate_type>(dp, surf);

    calder::mesh_calculator<SPACE> calc;
    // this needs to be reworked to take a surface topology and a position set
    coordinate_array dpn = calc.harmonicAvg(surf, dpf, p0, 1.0 * reg);

    for (int i = 0; i < p0.size(); i++) {
      pp[i] += dpn[i];
    }
    return pp;
  };
#endif

  void update(typename constraint_set<SPACE>::ptr constraint_set,
              int iter = 1) {
    real h = 1.0 / real(iter);

    const std::vector<vec3> &positions = constraint_set->get_positions();
    std::vector<vec3> vel(_p1);
    for (int i = 0; i < _p1.size(); i++) {
      vel[i] = _p1[i] - _p0[i];
    }

    constraint_set->set_positions(_p0);
    constraint_set->init_rest();

    backward_euler<SPACE> solver;
    solver.solve(constraint_set, vec_interface<SPACE>::to(vel), h);
    _p1 = constraint_set->get_positions();
    //    coordinate_array pp = project(_p0, _p1, constraint_set);
    //    constraint_set->set_positions(pp);
  }

  std::vector<vec3> _p1;
  std::vector<vec3> _p0;
}; // class position_optimizer
} // namespace hepworth
#endif