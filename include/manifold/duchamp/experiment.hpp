#ifndef __EXPERIMENT_HARNESS__
#define __EXPERIMENT_HARNESS__

#include "manifold/asawa/asawa.h"

#include "manifold/hepworth/constraints_init.hpp"
#include "manifold/hepworth/objective_function.hpp"
#include "manifold/hepworth/optimizer.hpp"
namespace duchamp {

template <typename SPACE>
class edge_init_policy
    : public asawa::edge_policy_t<SPACE, typename SPACE::real> {
public:
  M2_TYPEDEFS;
  edge_init_policy(typename SPACE::edge_index id)
      : edge_policy_t<SPACE, real>(id) {}

  virtual void calc(int i, edge_ptr &e, asawa::op_type op) {
    this->_vals[4 * i + 0] = -9999;
    this->_vals[4 * i + 1] = -9999;
    this->_vals[4 * i + 2] = -9999;
    this->_vals[4 * i + 3] = -9999;
    return;
  }
};

using namespace asawa;
template <typename SPACE> class simple_experiment_base {
public:
  M2_TYPEDEFS;

  typedef std::shared_ptr<simple_experiment_base<SPACE>> ptr;
  static ptr create() {
    return std::make_shared<simple_experiment_base<SPACE>>();
  }

  virtual void init() {
    asawa::obj_loader<SPACE> load;
    asawa::affine<SPACE> aff;
    asawa::remesh<SPACE> rem;

    _surf = &load("assets/bunny.obj");
    rem.triangulate(_surf);
    _surf->update_all();
    _surf->pack();

    aff.centerGeometry(*_surf);

    init_dynamic_mesh();
    this->_init();
  }

  virtual void init_checkpoint(const std::string &start_frame = "") {

    if (!start_frame.empty()) {
      FILE *file;
      file = fopen(start_frame.c_str(), "rb");
      if (file != NULL) {
        this->load_gaudi(start_frame);
      }
    }
    asawa::affine<SPACE> aff;
    aff.centerGeometry(*_surf);
    init_dynamic_mesh();
  }

  virtual void init_dynamic_mesh() {
    _integrator = new asawa::surf_integrator<SPACE>(_surf, 0.5, 2.5, 0.75);
    _integrator->template add_default_vertex_policy<typename SPACE::real>(
        SPACE::vertex_index::SMOOTH);
  }

  virtual void step(const int &frame) {
    _surf->update_all();
    _surf->reset_flags();
    _surf->pack();
    _integrator->integrate();
    this->_step(frame);
  }

  virtual void _init() = 0;
  virtual void _step(const int &frame) = 0;

  void init_phi(asawa::surf<SPACE> *surf) {
    M2_TYPEDEFS;
    std::vector<edge_ptr> edges = surf->get_edges();
    for (auto e : edges) {
      e->template set<real>(SPACE::edge_index::PHI0, -9999);
    }
  }

  virtual asawa::surf<SPACE> *get_surf() { return _surf; }

  void calc_stdev(const std::vector<real> &K, real &mean, real &stdev, real &mn,
                  real &mx) {
    real sum = std::accumulate(std::begin(K), std::end(K), 0.0);
    mean = sum / K.size();
    real accum = 0.0;
    std::for_each(std::begin(K), std::end(K),
                  [&](const double d) { accum += (d - mean) * (d - mean); });
    stdev = sqrt(accum / (K.size() - 1));

    mn = K[0];
    mx = K[0];

    std::for_each(std::begin(K), std::end(K), [&](const double d) {
      mn = std::min(mn, d);
      mx = std::max(mx, d);
    });
    std::cout << " K mean: " << mean << " stdev: " << stdev << std::endl;
  }

  vector<asawa::colorRGB> getColor() {
    M2_TYPEDEFS;
    asawa::cotan_curvature<SPACE> curve(_surf);
    std::vector<real> K = curve();
    real mean, stdev, mn, mx;
    calc_stdev(K, mean, stdev, mn, mx);

    std::vector<asawa::colorRGB> vert_colors(_surf->get_vertices().size());

    int i = 0;
    for (auto v : _surf->get_vertices()) {
      typename SPACE::real k = (K[i] - mn) / (mx - mn);
      typename SPACE::real N = 0;

      typename SPACE::coordinate_type colorS(1.0, 0.0, 0.0);
      typename SPACE::coordinate_type colorC(0.56, 0.60, 0.40);
      typename SPACE::coordinate_type colorH(0.5, 0.0, 0.5);
      typename SPACE::coordinate_type colorK(0.75, 0.75, 0.0);
      typename SPACE::coordinate_type mx = colorC;
      // mx = va::mix(k, colorK, mx);
      vert_colors[i] = asawa::colorRGB(mx[0], mx[1], mx[2], 1.0);
      i++;
    }
    return vert_colors;
  }

  virtual void load_gaudi(std::string file_name) {
    std::cout << " loading" << std::endl;
    asawa::flattened_surf<SPACE> fsurf;
    fsurf.clear();
    fsurf.read(file_name);
    _surf = fsurf.to_surf();
    //_integrator->set_mesh(_meshGraph);
  }

  virtual void dump_gaudi(std::string fname, int frame = 0) {
    std::cout << " dumping" << std::endl;
    asawa::flattened_surf<SPACE> fsurf(_surf);
    fsurf.write(fname + "." + std::to_string(frame) + ".gaudi");
  }

  virtual void dump_obj(std::string fname, int frame = 0) {
    _surf->pack();
    std::stringstream ss;
    ss << fname << "." << frame << ".obj";
    asawa::write_obj<SPACE>(*_surf, ss.str());
  }

  asawa::surf_integrator<SPACE> *_integrator;
  asawa::surf<SPACE> *_surf;
};

template <typename SPACE>
class bend_experiment : public simple_experiment_base<SPACE> {
public:
  M2_TYPEDEFS;

  typedef std::shared_ptr<bend_experiment<SPACE>> ptr;
  static ptr create() { return std::make_shared<bend_experiment<SPACE>>(); }

  virtual void _init() {
    this->init_phi(this->_surf);
    this->_integrator->add_edge_policy(
        new edge_init_policy<SPACE>(SPACE::edge_index::PHI0));
  };

  void willmore(const std::vector<typename SPACE::vec3> &p0) {

    typename hepworth::constraint_set<SPACE>::ptr constraints =
        hepworth::constraint_set<SPACE>::create(this->_surf);
    std::cout << "adding willmore constraints " << std::endl;
    init_willmore_constraints<SPACE>(this->_surf, constraints, 1.0e-8);

    std::vector<typename SPACE::vec3> p1(p0);
    hepworth::velocity_optimizer<SPACE> opt(p0, p1);
    opt.update(constraints);

    p1 = constraints->get_positions();
    asawa::ci::set_coordinates<SPACE>(p1, this->_surf);
  }

  void bending(const std::vector<typename SPACE::vec3> &p0,
               const std::vector<typename SPACE::vec3> &p1) {
    std::cout << "building constraints " << std::endl;
    typename hepworth::constraint_set<SPACE>::ptr constraints =
        hepworth::constraint_set<SPACE>::create(this->_surf);

    std::cout << "adding constraints " << std::endl;
    init_stretch_constraints<SPACE>(this->_surf, constraints, 4.0e-4);
    init_bend_constraints<SPACE>(this->_surf, constraints, 4.0e-6);
    // init_mem_bend_constraints<growth>(this->_surf, constraints, 1.0e-6,
    //                                    1.0e-5);

    constraints->add_constraint(hepworth::internal_collisions<SPACE>::create(
        this->_surf, this->_integrator->_max));

    hepworth::velocity_optimizer<SPACE> opt(p0, p1);
    opt.update(constraints);
    std::vector<typename SPACE::vec3> pp = constraints->get_positions();

    asawa::ci::set_coordinates<SPACE>(pp, this->_surf);
  }

  void update(const std::vector<typename SPACE::vec3> &p0) {
    asawa::ci::set_coordinates<SPACE>(p0, this->_surf);
  }

  virtual void _step(const int &frame){};
};

} // namespace duchamp
#endif