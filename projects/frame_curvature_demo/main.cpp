#include <algorithm>
#include <exception>
#include <iostream>
#include <string>
#include <vector>
#if defined(WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#ifdef ERROR
#undef ERROR
#endif
#endif

#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/asawa/shell/operations.hpp" // triangulate
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/asawa/datums.hpp"
#include "gaudi/kusama/vector_dirichlet_guided.hpp"
#include "gaudi/geometry_logger.hpp"

#include <Eigen/Dense>

#include "GaudiGraphics/buffers.hpp"
#include "GaudiGraphics/geometry_logger.h"
#include "GaudiGraphics/mesh_helper.hpp"
#include "GaudiGraphics/viewer.hpp"

using std::cerr;
using std::cout;
using std::endl;

namespace {

void sync_gaudi_debug_to_gg() {
  const auto &lines = gaudi::geometry_logger::get_lines();
  const auto &line_cols = gaudi::geometry_logger::get_line_colors();
  for (size_t i = 0; i + 1 < lines.size(); i += 2) {
    gg::geometry_logger::line(lines[i], lines[i + 1], line_cols[i]);
  }
}

} // namespace

class Scene;
using ScenePtr = std::shared_ptr<Scene>;

class Scene : public gg::Scene {
public:
  static ScenePtr create(int stride, bool butterfly_stencil, bool guided_edges) {
    return std::make_shared<Scene>(stride, butterfly_stencil, guided_edges);
  }

  Scene(int stride, bool butterfly_stencil, bool guided_edges)
      : gg::Scene(), __stride(stride), __butterfly(butterfly_stencil),
        __guided_edges(guided_edges) {
    __M = gaudi::asawa::shell::load_bunny();
    gaudi::asawa::shell::triangulate(*__M);
    // Vector Dirichlet / guided solve require vertex ids 0..N-1 ≡ indices into x.
    if (!__M->verts_are_dense_packed())
      gaudi::asawa::shell::pack(*__M);
    std::vector<gaudi::vec3> &x = gaudi::asawa::get_vec_data(*__M, 0);
    gaudi::asawa::center(x, 2.0);
    _objs.resize(1);
    _objs[0] = gg::BufferObject::create();
    _objs[0]->init();
    mSceneObjects.push_back(_objs[0]);
    mSceneObjects.push_back(gg::geometry_logger::get_instance().debugLines);
    _mesh_dirty = true;
  }

  void refresh_mesh_buffer() {
    gg::fillBuffer_ref(*__M, _objs[0],
                       gg::colorRGB(0.72, 0.74, 0.78, 1.0));
  }

  void draw_curvature_frames() {
    gaudi::geometry_logger::clear();
    gaudi::asawa::shell::shell &M = *__M;
    std::vector<gaudi::vec3> &x = gaudi::asawa::get_vec_data(M, 0);
    const gaudi::real scale =
        gaudi::asawa::shell::avg_length(M, x) * 3.0;
    auto fr = M.get_face_range(true);
    const auto stencil =
        __butterfly
            ? gaudi::asawa::shell::face_curvature_stencil::butterfly
            : gaudi::asawa::shell::face_curvature_stencil::one_ring;
    if (!__guided_edges) {
      for (gaudi::asawa::shell::FaceId fi : fr) {
        gaudi::asawa::shell::face_curvature_frame fc =
            gaudi::asawa::shell::face_curvature_frame_fit(M, x, fi, stencil);
        gaudi::vec3 c = gaudi::asawa::shell::face_center(M, fi, x);
        gaudi::vec4 cmin(0.95, 0.35, 0.2, 1.0);
        gaudi::geometry_logger::line(c, c + scale * fc.t_min, cmin);
      }
    }

    if (__guided_edges) {
      ensure_guided_edge_field(M, x, stencil);
      if (__u_edge.size() == 2 * __nE && __nE > 0) {
        gaudi::vec4 c_edge(0.2, 0.85, 0.55, 1.0);
        for (gaudi::asawa::shell::CornerId c : M.get_edge_range()) {
          const int slot = static_cast<int>(c) / 2;
          const int e = __slot_map[static_cast<size_t>(slot)];
          if (e < 0 || e >= __nE)
            continue;
          gaudi::vec3 mid =
              0.5 * (x[M.vert(c)] + x[M.vert(M.next(c))]);
          gaudi::vec3 n_e = gaudi::kusama::edge_average_normal(
              M, x, __edge_inc[static_cast<size_t>(e)]);
          const double gp = __u_edge(e);
          const double gq = __u_edge(e + __nE);
          gaudi::vec3 dir = gaudi::kusama::edge_guidance_vector_3d(
              M, c, x, n_e, static_cast<gaudi::real>(gp),
              static_cast<gaudi::real>(gq));
          if (dir.norm() > 1e-12)
            dir.normalize();
          gaudi::geometry_logger::line(mid, mid + scale * dir, c_edge);
        }
      }
    }
    sync_gaudi_debug_to_gg();
  }

  void ensure_guided_edge_field(
      gaudi::asawa::shell::shell &M, std::vector<gaudi::vec3> &x,
      gaudi::asawa::shell::face_curvature_stencil stencil) {
    if (__guided_solve_done)
      return;
    __guided_solve_done = true;
    try {
      __nE = gaudi::kusama::build_compact_edge_dof_map(M, __slot_map);
      if (__nE <= 0)
        return;
      gaudi::kusama::build_edge_incident_faces(M, __slot_map, __nE,
                                                 __edge_inc);
      __u_edge = gaudi::kusama::solve_curvature_guided_vector_dirichlet(
          M, x, 3.0, stencil, true);
    } catch (const std::exception &ex) {
      cerr << "[frame_curvature_demo] guided field solve skipped: " << ex.what()
           << "\n";
      __u_edge.resize(0);
      __nE = 0;
    }
  }

  void mark_mesh_dirty() { _mesh_dirty = true; }

  virtual void onAnimate(int /*frame*/) {
    draw_curvature_frames();
    gg::geometry_logger::render();
  }

  virtual void onDraw(gg::Viewer &viewer) {
    if (_mesh_dirty) {
      refresh_mesh_buffer();
      _mesh_dirty = false;
    }
    std::for_each(mSceneObjects.begin(), mSceneObjects.end(),
                  [&](gg::DrawablePtr obj) mutable {
                    if (obj->isVisible)
                      obj->draw(viewer.getProjection(), viewer.getModelView());
                  });
    gg::geometry_logger::clear();
  }

  gaudi::asawa::shell::shell::ptr __M;
  int __stride;
  bool __butterfly;
  bool __guided_edges = false;
  bool __guided_solve_done = false;
  int __nE = 0;
  std::vector<int> __slot_map;
  std::vector<std::vector<gaudi::asawa::shell::FaceId>> __edge_inc;
  Eigen::VectorXd __u_edge;

private:
  std::vector<gg::DrawablePtr> mSceneObjects;
  std::vector<gg::BufferObjectPtr> _objs;
  bool _mesh_dirty = true;
};

class App;
using AppPtr = std::shared_ptr<App>;

class App : public gg::SimpleApp {
public:
  static AppPtr create(int width, int height, int stride, bool butterfly,
                       bool guided) {
    return std::make_shared<App>(width, height, stride, butterfly, guided);
  }

  App(int width, int height, int stride, bool butterfly, bool guided)
      : gg::SimpleApp(width, height, 4.0, false, "frame_curvature_") {
    scene = Scene::create(stride, butterfly, guided);
    this->setScene(scene);
    this->initUI();
  }

  void initUI() { performLayout(); }

  ~App() override = default;

  ScenePtr scene;
};

int main(int argc, char *argv[]) {
  try {
    int stride = 40;
    bool butterfly = false;
    bool guided = false;
    for (int i = 1; i < argc; ++i) {
      std::string a(argv[i]);
      if (a == "--stride" && i + 1 < argc) {
        stride = std::max(1, std::atoi(argv[++i]));
      } else if (a == "--butterfly") {
        butterfly = true;
      } else if (a == "--guided") {
        guided = true;
      }
    }

    cout << "[frame_curvature_demo] default: per-face t_min | --butterfly "
            "stencil | --guided: smoothed edge vectors only (solve once)\n";

    nanogui::init();
    AppPtr app = App::create(1280, 740, stride, butterfly, guided);
    app->drawAll();
    app->setVisible(true);
    nanogui::mainloop();
    nanogui::shutdown();
  } catch (const std::runtime_error &e) {
    std::string error_msg =
        std::string("Caught a fatal error: ") + std::string(e.what());
    std::cerr << error_msg << endl;
#if defined(WIN32)
    MessageBoxA(nullptr, error_msg.c_str(), nullptr, MB_ICONERROR | MB_OK);
#endif
    return -1;
  }
  return 0;
}
