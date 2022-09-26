//#include "nanoguiincludes.h"

#include "manifold/coordinate_interface.hpp"

#include "manifold/hepworth/constraints_init.hpp"
#include "manifold/hepworth/objective_function.hpp"

#include "manifold/m2.hpp"
#include <algorithm>
#include <exception>
#include <vector>
#if defined(WIN32)
#include <windows.h>
#endif

#include <stdio.h> /* defines FILENAME_MAX */
// #define WINDOWS  /* uncomment this line to use it for windows.*/
#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

#include <complex>
#include <iostream>
#include <random>
#include <string>

#include "GaudiGraphics/buffers.hpp"
#include "GaudiGraphics/geometry_logger.h"
#include "GaudiGraphics/mesh_helper.hpp"
#include "GaudiGraphics/viewer.hpp"

#include "manifold/conj_grad.hpp"
#include "manifold/diffuse.hpp"
#include "manifold/laplacian.hpp"

#include "manifold/bins.hpp"
#include "manifold/m2Includes.h"
#include "manifold/m2Operators.h"
#include "manifold/make.hpp"
#include "manifold/objloader.hpp"

#include "manifold/hepworth/optimizer.hpp"

#include "manifold/harmonic_integrators.hpp"
#include "manifold/vec_addendum.h"

#include "manifold/triangle_operations.hpp"

//#include "m2Operators.h"

#define TRACKBALLSIZE (0.8f)
#define RENORMCOUNT 97

using std::cerr;
using std::cout;
using std::endl;

using namespace GaudiMath;

class Scene;
using ScenePtr = std::shared_ptr<Scene>;

class Scene : public gg::Scene {

public:
  static ScenePtr create() { return std::make_shared<Scene>(); }

  Scene() : gg::Scene() {
    initScene();
    // initUI();
  }

  void initScene() {

    m2::obj_loader<stretch> load;
    m2::subdivide<stretch> sub;
    m2::make<stretch> mk;
    m2::convex_hull<stretch> ch;

    m2::construct<stretch> bevel;
    m2::affine<stretch> mod;
    std::string start_frame = "";
    // std::string start_frame = "stretch.46.gaudi";
    if (!start_frame.empty()) {
      FILE *file;
      file = fopen(start_frame.c_str(), "rb");
      if (file != NULL) {
        this->load_gaudi(start_frame);
      }
    } else {

      _meshGraph = &load("assets/bunny.obj");
      //_meshGraph = &load("assets/close.obj");

      //_meshGraph = &load("assets/icosahedron.obj");
      //_meshGraph = &load("assets/sphere.obj");

      //_meshGraph = &load("assets/messer.obj");
      //_meshGraph = &load("assets/tet.obj");
      // std::cout << "--make cube" << std::endl;
      //_meshGraph = mk.cube(1.0, 1.0, 1.0);
      //_meshGraph = mk.tet();

      //_meshGraph = &sub.subdivide_control(*_meshGraph);
      //_meshGraph = &sub.subdivide_control(*_meshGraph);
      //_meshGraph = &sub.subdivide_control(*_meshGraph);
      //_meshGraph = &sub.subdivide_control(*_meshGraph);
      //_meshGraph = &sub.subdivide_control(*_meshGraph);
      //_meshGraph = &sub.subdivide_control(*_meshGraph);
      //_meshGraph = &sub.subdivide_control(*_meshGraph);

      m2::remesh<stretch> rem;
      rem.triangulate(_meshGraph);
      //_meshGraph = &sub.subdivide_control(*_meshGraph);

      std::cout << "--center" << std::endl;
      std::cout << "--update_all" << std::endl;
      _meshGraph->update_all();
      std::cout << "--pack" << std::endl;
      _meshGraph->pack();

      mod.centerGeometry(*_meshGraph);
    }

    int N = 0;
    _integrator = new m2::surf_integrator<stretch>(_meshGraph, 0.5, 2.5, 0.75);
    //_integrator = new m2::surf_integrator<stretch>(_meshGraph, 0.1, 3.0, 0.5);
    _integrator->add_default_vertex_policy<typename stretch::real>(
        stretch::vertex_index::SMOOTH);
    _integrator->add_default_vertex_policy<typename stretch::real>(
        stretch::vertex_index::ACT);
    _integrator->add_vertex_policy(
        new action_policy<stretch>(stretch::vertex_index::ACT));
    _max = _integrator->_max;
    _min = _integrator->_min;

    std::cout << "--init rx" << std::endl;

    std::cout << "creating buffer" << std::endl;

    _obj = gg::BufferObject::create();
    _obj->init();
    mSceneObjects.push_back(_obj);

    mSceneObjects.push_back(gg::geometry_logger::get_instance().debugLines);
  }

  template <typename SPACE>
  vector<m2::colorRGB> getColor(m2::surf<SPACE> *surf) {
    M2_TYPEDEFS;
    auto smooth = m2::ci::get<SPACE, real>(surf, SPACE::vertex_index::SMOOTH);
    auto hot = m2::ci::get<SPACE, real>(surf, SPACE::vertex_index::ACT);

    m2::cotan_curvature<SPACE> curve(surf);
    std::vector<real> K = curve();
    real mn = K[0];
    real mx = K[0];
    std::for_each(std::begin(K), std::end(K), [&](const double d) {
      mn = std::min(mn, d);
      mx = std::max(mx, d);
    });

    std::vector<m2::colorRGB> vert_colors(smooth.size());
    int i = 0;
    for (auto v : surf->get_vertices()) {
      typename SPACE::real k = (K[i] - mn) / (mx - mn);
      typename SPACE::real N = 0;
      typename SPACE::real s = smooth[i];
      typename SPACE::real h = hot[i];

      typename SPACE::coordinate_type colorS(1.0, 0.0, 0.0);
      typename SPACE::coordinate_type colorC(0.56, 0.50, 0.60);
      typename SPACE::coordinate_type colorH(0.5, 0.0, 0.5);
      typename SPACE::coordinate_type colorK(0.5, 0.5, 0.0);
      typename SPACE::coordinate_type mx = colorC;
      mx = m2::va::mix(k, colorK, mx);
      mx = m2::va::mix(s, colorS, mx);
      mx = m2::va::mix(h, colorH, mx);
      mx = m2::va::mix(k, colorC, mx);

      vert_colors[i] = m2::colorRGB(mx[0], mx[1], mx[2], 1.0);
      i++;
    }
    return vert_colors;
  }

  virtual void onAnimate(int frame) {

    _meshGraph->update_all();
    _meshGraph->reset_flags();
    _meshGraph->pack();
    // if (frame == 1)
    _integrator->integrate();
    //    if (frame > 1)
    //      return;

    std::cout << "frame: " << frame << std::endl;
    std::cout << "====== " << std::endl;
    std::cout << " verts: " << _meshGraph->get_vertices().size() << std::endl;
    std::cout << " edges: " << _meshGraph->get_edges().size() << std::endl;
    std::cout << " faces: " << _meshGraph->get_faces().size() << std::endl;
    std::cout << " mean edge length: "
              << m2::ci::geometric_mean_length<stretch>(_meshGraph)
              << std::endl;

    double dt = _params.dt;
    auto colors = getColor(_meshGraph);
#if 1

    std::vector<stretch::vec3> p0 =
        m2::ci::get_coordinates<stretch>(_meshGraph);
    std::vector<stretch::vec3> p1(p0);

    // build constraints to capture current config
    std::cout << "====== " << std::endl;
    std::cout << "integrating " << std::endl;
    double Nn = 0, NRxn = 0;
    using edge_ptr = typename m2::surf<stretch>::edge *;
    using vertex_ptr = typename m2::surf<stretch>::vertex *;

    std::cout << "building constraints " << std::endl;
    hepworth::constraint_set<stretch>::ptr constraints =
        hepworth::constraint_set<stretch>::create(_meshGraph);

    std::cout << "adding constraints " << std::endl;
    // init_stretch_constraints<stretch>(_meshGraph, constraints, 1e-5);
    //       init_cross_constraints<stretch>(_meshGraph, constraints);
    init_bend_constraints<stretch>(_meshGraph, constraints, 5e-5);

    constraints->add_constraint(
        hepworth::internal_collisions<stretch>::create(_meshGraph, _max));

    std::vector<stretch::vec3> normals =
        m2::ci::get_vertex_normals<stretch>(_meshGraph);

    std::vector<stretch::vec3> forces(p0.size(), stretch::vec3::Zero());
    std::vector<stretch::vec3> wN(p0);
    int i = 0;
#if 0
    // perturb config
    int iii = 0;
    forces[iii] = 0.01 * normals[iii];
    p1[iii] += forces[iii];

#elif 1
    std::mt19937_64 rng;
    // initialize the random number generator with time-dependent seed
    // rng.seed(419770629407);
    // initialize a uniform distribution between 0 and 1
    std::uniform_real_distribution<double> unif(-1.0, 1.0);
    std::uniform_real_distribution<double> unif01(0.0, 1.0);
    auto &verts = _meshGraph->get_vertices();

    m2::cotan_curvature<stretch> curve(_meshGraph);
    std::vector<double> K = curve();
    double sum = std::accumulate(std::begin(K), std::end(K), 0.0);
    double mean = sum / K.size();
    double accum = 0.0;
    std::for_each(std::begin(K), std::end(K),
                  [&](const double d) { accum += (d - mean) * (d - mean); });
    double stdev = sqrt(accum / (K.size() - 1));

    std::cout << " K mean: " << mean << " stdev: " << stdev << std::endl;

    double N = double(verts.size());
    double d0 = 0.0;
    double dp = 0.0;
    double dn = 0.0;
    for (auto &v : verts) {
      double a = v->template get<double>(stretch::vertex_index::ACT);
      d0 += sqrt(a * a) / N;
    }

    std::cout << "current density: " << d0 << std::endl;
    double max_d = 0.01;
    double prob_d = 1.0 / N / max_d;
    double d1 = d0;
    std::cout << "max density: " << max_d << ", prob density: " << prob_d
              << std::endl;

    for (auto &v : verts) {
      double a = v->template get<double>(stretch::vertex_index::ACT);
      forces[i].setZero();

      if (unif01(rng) < prob_d && d1 < max_d && a < 1e-1) {
        if (dp < 16.0 * dn) {
          a = 1.0;
          dp += fabs(a) / N;
        } else {
          a = -1.0;
          dn += fabs(a) / N;
        }
        d1 += sqrt(a * a) / N;
      }

      if (a * a > 0) {
        // if (unif01(rng) > 0.99)
        //   a = 0.0;
        a *= (1.0 - 0.002 * unif01(rng));
      }
      // std::cout << forces[i].transpose() << " - " << a << std::endl;
      // double mag = a < 0 ? -a : a;
      // p1[i] += a * 0.03 * normals[i];
      wN[i] = stretch::vec3(a * normals[i]);
      v->template set<double>(stretch::vertex_index::ACT, a);

      i++;
    }

#elif 0
    for (int i = 0; i < p1.size(); i++) {
      forces[i] = 0.01 * normals[i];
      p1[i] += forces[i];
    }

#endif
#if 1

    i = 0;
    double C_N = 0.01;
    for (auto v : verts) {
      auto f = wN[i].norm();
      auto sn = va::sgn(wN[i].dot(normals[i]));
      for_each_vertex<stretch>(
          v, [i, &p0, &p1, &normals, f, sn,
              C_N](typename m2::surf<stretch>::face_vertex_ptr fv) {
            int j = fv->next()->vertex()->position_in_set();
            auto Nj = normals[j];
            p1[j] += C_N * f * sn * Nj;
          });
      p1[i] -= 0.5 * C_N * f * sn * normals[i];
      i++;
    }
#endif
#if 0
    i = 0;

    double C_C = C_N;
    for (auto v : verts) {
      auto f = wN[i].norm();
      auto sn = va::sgn(wN[i].dot(normals[i]));

      for_each_vertex<stretch>(
          v, [i, &p0, &p1, &normals, f, sn,
              C_C](typename m2::surf<stretch>::face_vertex_ptr fv) {
            int j = fv->next()->vertex()->position_in_set();
            auto Ni = normals[i];
            auto pi = p0[i];
            auto pj = p0[j];
            typename stretch::coordinate_type dp = pj - pi;
            dp.normalize();
            auto Ncp = dp.cross(Ni);
            p1[j] += 0.5 * C_C * f * sn * Ncp;
          });
      i++;
    }
#endif
#if 1
    std::vector<stretch::vec3> fN =
        calcPotential<stretch>(_meshGraph, wN, p0, 2.0 * _max);
    i = 0;
    for (auto f : fN) {
      p1[i] += 0.25 * f;
      i++;
    }
    gg::geometry_logger::field(p0, fN, 1.0);

#endif
    // hepworth::position_optimizer<stretch> opt(p0, p1, weights);

    // gg::geometry_logger::field(p0, wN, 0.1, gg::PresetColor::red);

    hepworth::velocity_optimizer<stretch> opt(p0, p1);
    opt.update(constraints);
    p1 = constraints->get_positions();
    m2::ci::set_coordinates<stretch>(p1, _meshGraph);

#if 1

    // std::cout << "print vecs" << std::endl;
    std::cout << "rendering debug" << std::endl;
    gg::geometry_logger::render();
#endif

    _meshGraph->print();

    gg::fillBuffer(_meshGraph, _obj, colors);
#endif
  }

  virtual void save(int frame) {
    _meshGraph->pack();
    std::stringstream ss;
    ss << "stretch." << frame << ".obj";
    m2::write_obj<stretch>(*_meshGraph, ss.str());
  }

  virtual void load_gaudi(std::string file_name) {
    std::cout << " loading" << std::endl;
    m2::flattened_surf<stretch> fsurf;
    fsurf.clear();
    fsurf.read(file_name);
    _meshGraph = fsurf.to_surf();
    //_integrator->set_mesh(_meshGraph);
  }

  virtual void dump_gaudi(int frame = 0) {
    std::cout << " dumping" << std::endl;
    m2::flattened_surf<stretch> fsurf(_meshGraph);
    fsurf.write("stretch." + std::to_string(frame) + ".gaudi");
  }

  virtual void onDraw(gg::Viewer &viewer) {

    std::for_each(mSceneObjects.begin(), mSceneObjects.end(),
                  [&](gg::DrawablePtr obj) mutable {
                    if (obj->isVisible)
                      obj->draw(viewer.getProjection(), viewer.getModelView());
                  });

    gg::geometry_logger::clear();
  }

  struct {
    double dt = 0.01;
  } _params;

private:
  double _max = 0.0;
  double _min = 0.0;

  std::vector<gg::DrawablePtr> mSceneObjects;

  gg::BufferObjectPtr _obj = NULL;
  m2::surf<stretch> *_meshGraph;
  m2::surf_integrator<stretch> *_integrator;
};

std::string GetCurrentWorkingDir(void) {
  char buff[FILENAME_MAX];
  GetCurrentDir(buff, FILENAME_MAX);
  std::string current_working_dir(buff);
  return current_working_dir;
}

class App;
using AppPtr = std::shared_ptr<App>;

class App : public gg::SimpleApp {
public:
  static AppPtr create(std::string file) { return std::make_shared<App>(file); }

  typedef double Real;

  App(std::string file) : gg::SimpleApp() {
    this->setScene(scene = Scene::create());
    this->initUI();
  }

  void initUI() {
    using namespace nanogui;
    int w = 256;
    performLayout();
    // window->center();
  }

  ~App() {}

  ScenePtr scene;
};

int main(int argc, char *argv[]) {
  try {
    cout << "You have entered " << argc << " arguments:"
         << "\n";

    for (int i = 0; i < argc; ++i)
      cout << argv[i] << "\n";

    nanogui::init();

    AppPtr app = App::create(std::string(argv[0]));

    // app->setScene(Scene::create());
    app->drawAll();
    app->setVisible(true);
    nanogui::mainloop();
    // delete app;
    nanogui::shutdown();

  } catch (const std::runtime_error &e) {
    std::string error_msg =
        std::string("Caught a fatal error: ") + std::string(e.what());

#if defined(WIN32)
    MessageBoxA(nullptr, error_msg.c_str(), NULL, MB_ICONERROR | MB_OK);
#else
    std::cerr << error_msg << endl;
#endif

    return -1;
  }

  return 0;
}
