#include <algorithm>
#include <cmath>
#include <cstring>
#include <exception>
#include <iostream>
#include <limits>
#include <numeric>
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

#include "GaudiGraphics/buffers.hpp"
#include "GaudiGraphics/geometry_logger.h"
#include "GaudiGraphics/mesh_helper.hpp"
#include "GaudiGraphics/viewer.hpp"

#include "gaudi/asawa/datums.hpp"
#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/calder/least_squares_fit.hpp"
#include "gaudi/duchamp/darboux_cyclide_medial.hpp"
#include "gaudi/geometry_logger.hpp"
#include "gaudi/test/darboux_cyclide_torus_fixture.hpp"

namespace {

class Scene;
using ScenePtr = std::shared_ptr<Scene>;

gaudi::duchamp::cyclide_medial_params demo_params() {
  gaudi::duchamp::cyclide_medial_params params;
  params.l0_scale = 2.0;
  params.fit_p = 3.0;
  params.normal_l0 = 2.0;
  params.min_normal_alignment = -0.5;
  params.medial_smooth_scale = 8.0;
  params.medial_smooth_blend = 1.0;
  params.max_iters = 100;
  params.tol = 1e-8;
  params.max_travel_scale = 1000.0;
  params.newton_step_scale = 1000.0;
  return params;
}

void normalize_mesh(gaudi::asawa::shell::shell &M,
                    gaudi::real target_radius = 1.5) {
  std::vector<gaudi::vec3> &x = gaudi::asawa::get_vec_data(M, 0);
  if (x.empty()) {
    return;
  }

  gaudi::vec3 lo = x.front();
  gaudi::vec3 hi = x.front();
  for (const gaudi::vec3 &p : x) {
    lo = lo.cwiseMin(p);
    hi = hi.cwiseMax(p);
  }

  const gaudi::vec3 center = 0.5 * (lo + hi);
  gaudi::real radius = 0.0;
  for (const gaudi::vec3 &p : x) {
    radius = std::max(radius, (p - center).norm());
  }
  if (radius < 1e-12) {
    return;
  }

  const gaudi::real scale = target_radius / radius;
  for (gaudi::vec3 &p : x) {
    p = scale * (p - center);
  }
  std::cerr << "normalized mesh radius=" << target_radius
            << " scale=" << scale << std::endl;
}

void inspect_exact_convex_fit(gaudi::asawa::shell::shell &M, int vertex,
                              const gaudi::duchamp::cyclide_medial_params &params) {
  const std::vector<gaudi::vec3> &x = gaudi::asawa::const_get_vec_data(M, 0);
  const std::vector<gaudi::vec3> vertex_normals =
      gaudi::asawa::shell::vertex_normals(M, x);
  const std::vector<gaudi::vec3> face_centers =
      gaudi::asawa::shell::face_centers(M, x);
  const std::vector<gaudi::vec3> face_normals =
      gaudi::asawa::shell::face_normals(M, x);
  const std::vector<gaudi::real> face_areas =
      gaudi::asawa::shell::face_areas(M, x);
  const gaudi::real avg_len = gaudi::asawa::shell::avg_length(M, x);
  const gaudi::real l0 = std::max(params.l0_scale * avg_len, gaudi::real(1e-12));

  vertex = std::clamp(vertex, 0, static_cast<int>(x.size()) - 1);
  const gaudi::vec3 pov = x[static_cast<size_t>(vertex)];
  const gaudi::vec3 Ni = vertex_normals[static_cast<size_t>(vertex)].normalized();

  gaudi::albers::normal_constrained_darboux_cyclide exact_acc;
  gaudi::real w_sum = 0.0;
  gaudi::real area_sum = 0.0;
  gaudi::real max_w = 0.0;
  gaudi::real min_dist = std::numeric_limits<gaudi::real>::infinity();
  gaudi::real max_dist = 0.0;
  gaudi::vec3 weighted_center = gaudi::vec3::Zero();
  gaudi::vec3 weighted_normal = gaudi::vec3::Zero();
  int accepted_faces = 0;

  for (int fi = 0; fi < static_cast<int>(face_centers.size()); ++fi) {
    const gaudi::vec3 dp = face_centers[static_cast<size_t>(fi)] - pov;
    const gaudi::real dist = dp.norm();
    if (dist < 1e-12) {
      continue;
    }
    const gaudi::vec3 Nj = face_normals[static_cast<size_t>(fi)].normalized();
    const gaudi::real convexity =
        std::max(gaudi::real(0.0), dp.normalized().dot(Nj));
    const gaudi::real wf =
        convexity * gaudi::calder::calc_inv_dist(dp, l0, params.fit_p);
    const gaudi::real w_total = face_areas[static_cast<size_t>(fi)] * wf;
    if (w_total < 1e-8) {
      continue;
    }
    exact_acc.accumulate(w_total, dp, Nj);
    w_sum += w_total;
    area_sum += face_areas[static_cast<size_t>(fi)];
    max_w = std::max(max_w, w_total);
    min_dist = std::min(min_dist, dist);
    max_dist = std::max(max_dist, dist);
    weighted_center += w_total * face_centers[static_cast<size_t>(fi)];
    weighted_normal += w_total * Nj;
    ++accepted_faces;
  }

  const gaudi::albers::vec14 exact_Q = exact_acc.solve();
  const std::vector<gaudi::vec3> p_fit = {pov};
  const std::vector<gaudi::vec3> n_fit = {Ni};
  const std::vector<gaudi::albers::vec14> bh_Qs =
      gaudi::calder::darboux_cyclide_shell_fit(M, p_fit, n_fit, l0, params.fit_p);
  const gaudi::albers::vec14 bh_Q =
      bh_Qs.empty() ? gaudi::albers::vec14::Zero() : bh_Qs.front();

  const auto print_fit = [&](const char *label, const gaudi::albers::vec14 &Q) {
    const gaudi::vec3 g = gaudi::albers::darboux_grad(Q, gaudi::vec3::Zero());
    gaudi::real align = -1.0;
    if (g.allFinite() && g.norm() > 1e-12) {
      gaudi::vec3 gn = g.normalized();
      if (gn.dot(Ni) < 0.0) {
        gn *= -1.0;
      }
      align = gn.dot(Ni);
    }
    std::cerr << "  " << label
              << " D0=" << gaudi::albers::eval_darboux(Q, gaudi::vec3::Zero())
              << " |g0|=" << g.norm()
              << " grad_align=" << align << std::endl;
  };

  if (w_sum > 0.0) {
    weighted_center /= w_sum;
    weighted_normal /= w_sum;
  }
  std::cerr << "exact convex support vertex=" << vertex
            << " faces=" << accepted_faces
            << " w_sum=" << w_sum
            << " max_w_frac=" << (w_sum > 0.0 ? max_w / w_sum : 0.0)
            << " area_sum=" << area_sum
            << " min_dist=" << min_dist
            << " max_dist=" << max_dist
            << " centroid_offset=" << (weighted_center - pov).transpose()
            << " mean_normal_align="
            << (weighted_normal.norm() > 1e-12
                    ? weighted_normal.normalized().dot(Ni)
                    : 0.0)
            << std::endl;
  print_fit("barnes_hut", bh_Q);
  print_fit("exact_faces", exact_Q);
}

int run_headless_probe() {
  const auto params = demo_params();
  gaudi::asawa::shell::shell::ptr mesh = gaudi::asawa::shell::load_bunny();
  gaudi::asawa::shell::shell &M = *mesh;
  normalize_mesh(M);

  gaudi::duchamp::cyclide_medial_stats stats;
  const std::vector<gaudi::duchamp::cyclide_medial_candidate> candidates =
      gaudi::duchamp::compute_local_fit_cyclide_medial_candidates(M, params,
                                                                  &stats);
  const std::vector<gaudi::vec3> &x = gaudi::asawa::const_get_vec_data(M, 0);
  const gaudi::real avg_len = gaudi::asawa::shell::avg_length(M, x);
  std::vector<int> ranked;
  ranked.reserve(candidates.size());
  for (int i = 0; i < static_cast<int>(candidates.size()); ++i) {
    if (candidates[static_cast<size_t>(i)].accepted) {
      ranked.push_back(i);
    }
  }
  std::sort(ranked.begin(), ranked.end(), [&](int a, int b) {
    return candidates[static_cast<size_t>(a)].travel >
           candidates[static_cast<size_t>(b)].travel;
  });
  const int n_report = std::min<int>(8, ranked.size());
  std::cerr << "longest accepted medial rays:" << std::endl;
  for (int rank = 0; rank < n_report; ++rank) {
    const auto &cand = candidates[static_cast<size_t>(ranked[rank])];
    const gaudi::vec3 p = x[static_cast<size_t>(cand.vertex)];
    std::cerr << "  rank=" << rank << " vertex=" << cand.vertex
              << " travel=" << cand.travel
              << " travel/avg_edge=" << cand.travel / avg_len
              << " residual=" << cand.residual
              << " iterations=" << cand.iterations
              << " clamped_steps=" << cand.clamped_steps
              << " p=" << p.transpose()
              << " center=" << cand.center_world.transpose() << std::endl;
  }
  if (!ranked.empty()) {
    const auto &cand = candidates[static_cast<size_t>(ranked.front())];
    std::cerr << "worst medial trace vertex=" << cand.vertex
              << " samples=" << cand.trace_world.size() << std::endl;
    for (int i = 0; i < static_cast<int>(cand.trace_world.size()); ++i) {
      const gaudi::vec3 prev =
          i > 0 ? cand.trace_world[static_cast<size_t>(i - 1)]
                : cand.trace_world[static_cast<size_t>(i)];
      const gaudi::vec3 cur = cand.trace_world[static_cast<size_t>(i)];
      std::cerr << "  trace[" << i << "] p=" << cur.transpose()
                << " step=" << (cur - prev).norm() << std::endl;
    }
    gaudi::duchamp::cyclide_medial_params probe_params = params;
    probe_params.probe_vertex = cand.vertex;
    gaudi::duchamp::cyclide_medial_stats probe_stats;
    (void)gaudi::duchamp::compute_local_fit_cyclide_medial_candidates(
        M, probe_params, &probe_stats);
    inspect_exact_convex_fit(M, cand.vertex, params);
  }
  std::cerr << "headless bunny local-fit zero-to-medial" << std::endl;
  gaudi::duchamp::print_cyclide_medial_stats(std::cerr, stats);
  return 0;
}

void sync_gaudi_debug_to_gg() {
  const auto &lines = gaudi::geometry_logger::get_lines();
  const auto &line_cols = gaudi::geometry_logger::get_line_colors();
  for (size_t i = 0; i + 1 < lines.size(); i += 2) {
    gg::geometry_logger::line(lines[i], lines[i + 1], line_cols[i]);
  }

  const auto &points = gaudi::geometry_logger::get_points();
  const auto &point_cols = gaudi::geometry_logger::get_point_colors();
  for (size_t i = 0; i < points.size(); ++i) {
    gg::geometry_logger::point(points[i], point_cols[i]);
  }
}

class Scene : public gg::Scene {
public:
  static ScenePtr create(const gaudi::duchamp::cyclide_medial_params &params) {
    return std::make_shared<Scene>(params);
  }

  Scene(const gaudi::duchamp::cyclide_medial_params &params)
      : gg::Scene(), __params(params) {

    __M = gaudi::asawa::shell::load_bunny();

    gaudi::asawa::shell::shell &M = *__M;
    normalize_mesh(M);
    configure_scene_frame(M);
    __candidates =
        gaudi::duchamp::compute_local_fit_cyclide_medial_candidates(
            M, __params, &__stats);
    std::cerr << "bunny local-fit zero-to-medial" << std::endl;
    gaudi::duchamp::print_cyclide_medial_stats(std::cerr, __stats);

    _objs.resize(1);
    _objs[0] = gg::BufferObject::create();
    _objs[0]->init();
    refresh_mesh_buffer();
    mSceneObjects.push_back(_objs[0]);
    mSceneObjects.push_back(gg::geometry_logger::get_instance().debugLines);
  }

  void refresh_mesh_buffer() {
    gg::fillBuffer_ref(*__M, _objs[0], gg::colorRGB(0.72, 0.74, 0.78, 1.0));
  }

  void configure_scene_frame(gaudi::asawa::shell::shell &M) {
    const std::vector<gaudi::vec3> &x = gaudi::asawa::const_get_vec_data(M, 0);
    if (x.empty()) {
      return;
    }

    gaudi::vec3 lo = x.front();
    gaudi::vec3 hi = x.front();
    for (const gaudi::vec3 &p : x) {
      lo = lo.cwiseMin(p);
      hi = hi.cwiseMax(p);
    }

    __frame = gaudi::test::make_torus_frame(0.5 * (lo + hi), gaudi::vec3::UnitZ());
    __major_radius = 0.5 * (hi - lo).norm();
    __minor_radius = 4.0 * gaudi::asawa::shell::avg_length(M, x);
    std::cerr << "loaded bunny verts=" << M.vert_count()
              << " faces=" << M.face_count()
              << " scene_radius=" << __major_radius
              << " slice_extent=" << __minor_radius << std::endl;
  }

  void draw_medial_rays() {
    const gaudi::vec4 zero_color(0.0, 0.85, 1.0, 1.0);
    const gaudi::vec4 medial_color(1.0, 0.65, 0.05, 1.0);
    const gaudi::vec4 cylinder_axis_color(0.9, 0.15, 1.0, 1.0);
    const gaudi::vec4 axis_color(0.35, 0.35, 0.35, 1.0);

    gg::geometry_logger::line(
        __frame.center - 1.7 * __major_radius * __frame.z_axis,
        __frame.center + 1.7 * __major_radius * __frame.z_axis, axis_color);

    gaudi::asawa::shell::shell &M = *__M;
    std::vector<gaudi::vec3> &x = gaudi::asawa::get_vec_data(M, 0);
    const gaudi::real cylinder_axis_len =
        8.0 * gaudi::asawa::shell::avg_length(M, x);

    for (const auto &cand : __candidates) {
      const gaudi::vec3 p = x[cand.vertex];
      if (cand.projection_converged) {
        if ((cand.surface_world - p).norm() > 1e-12) {
          gg::geometry_logger::line(p, cand.surface_world, zero_color);
        } else {
          gg::geometry_logger::point(cand.surface_world, zero_color);
        }
      }
      if (cand.projection_converged && cand.center_world.allFinite() &&
          (cand.center_world - cand.surface_world).norm() > 1e-12) {
        gg::geometry_logger::line(cand.surface_world, cand.center_world,
                                  medial_color);
      }
      if (cand.accepted && cand.cylinder_axis_world.squaredNorm() > 1e-12) {
        gaudi::vec3 dp = cylinder_axis_len * cand.cylinder_axis_world.normalized();
        gg::geometry_logger::line(
            p - dp, p + dp,
            cylinder_axis_color);
      }
    }
  }

  virtual void onAnimate(int /*frame*/) {}

  virtual void onDraw(gg::Viewer &viewer) {
    draw_medial_rays();
    sync_gaudi_debug_to_gg();
    gg::geometry_logger::render();
    std::for_each(mSceneObjects.begin(), mSceneObjects.end(),
                  [&](gg::DrawablePtr obj) mutable {
                    if (obj->isVisible)
                      obj->draw(viewer.getProjection(), viewer.getModelView());
                  });
    gg::geometry_logger::clear();
  }

private:
  gaudi::asawa::shell::shell::ptr __M;
  gaudi::test::TorusFrame __frame;
  gaudi::duchamp::cyclide_medial_params __params;
  gaudi::duchamp::cyclide_medial_stats __stats;
  std::vector<gaudi::duchamp::cyclide_medial_candidate> __candidates;
  gaudi::real __major_radius = 1.25;
  gaudi::real __minor_radius = 0.35;
  std::vector<gg::DrawablePtr> mSceneObjects;
  std::vector<gg::BufferObjectPtr> _objs;
};

class App;
using AppPtr = std::shared_ptr<App>;

class App : public gg::SimpleApp {
public:
  static AppPtr create(int width, int height,
                       const gaudi::duchamp::cyclide_medial_params &params) {
    return std::make_shared<App>(width, height, params);
  }

  App(int width, int height,
      const gaudi::duchamp::cyclide_medial_params &params)
      : gg::SimpleApp(width, height, 4.0, false, "darboux_cyclide_medial_") {
    scene = Scene::create(params);
    this->setScene(scene);
    this->initUI();
  }

  void initUI() { performLayout(); }

  ~App() override = default;

  ScenePtr scene;
};

} // namespace

int main(int argc, char **argv) {
  if (argc > 1 && std::strcmp(argv[1], "--probe") == 0) {
    return run_headless_probe();
  }

  try {
    const auto params = demo_params();
    nanogui::init();
    AppPtr app = App::create(1280, 740, params);
    app->drawAll();
    app->setVisible(true);
    nanogui::mainloop();
    nanogui::shutdown();
  } catch (const std::runtime_error &e) {
    std::cerr << "Caught a fatal error: " << e.what() << std::endl;
    return -1;
  }
  return 0;
}
