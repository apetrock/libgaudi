#pragma once

#include <atomic>
#include <chrono>
#include <memory>
#include <thread>
#include <utility>

#include "gaudi/duchamp/demo_trait.hpp"
#include "gaudi/geometry_logger.hpp"
#include "gaudi/vermeer/debug_line_utils.hpp"
#include "gaudi/vermeer/duchamp_playback.hpp"
#include "gaudi/vermeer/scene_frame.hpp"
#include "lewitt/geometry_logger.h"

namespace gaudi {
namespace vermeer {

// Owns the sim worker thread. demo_trait is touched only on that thread after
// start(); the render thread consumes SceneFrameChannel only.
class duchamp_sim_runtime {
public:
  duchamp_sim_runtime(duchamp::demo_trait::ptr demo,
                      std::shared_ptr<duchamp_playback> playback)
      : _demo(std::move(demo)), _playback(std::move(playback)) {}

  ~duchamp_sim_runtime() { stop(); }

  duchamp_sim_runtime(const duchamp_sim_runtime &) = delete;
  duchamp_sim_runtime &operator=(const duchamp_sim_runtime &) = delete;

  SceneFrameChannel &channel() { return _channel; }
  const SceneFrameChannel &channel() const { return _channel; }

  duchamp::demo_trait::ptr demo() const { return _demo; }

  // Call once on the main thread before start() (e.g. after GPU init).
  void reset_demo() {
    if (_demo)
      _demo->reset();
  }

  void start() {
    if (_running.exchange(true))
      return;
    _worker = std::thread([this] { worker_loop(); });
  }

  void stop() {
    if (!_running.exchange(false))
      return;
    _channel.shutdown();
    if (_worker.joinable())
      _worker.join();
  }

private:
  static scene_debug_line to_scene_line(
      const lewitt::doables::lineable::exported_line &line) {
    scene_debug_line out;
    out.p0 = glm::vec3(line.p0.x, line.p0.y, line.p0.z);
    out.p1 = glm::vec3(line.p1.x, line.p1.y, line.p1.z);
    out.color = glm::vec3(line.color.x, line.color.y, line.color.z);
    out.radius = line.radius;
    return out;
  }

  std::shared_ptr<const SceneFrame> capture_frame(int frame) {
    auto frame_ptr = std::make_shared<SceneFrame>();
    frame_ptr->sim_frame = frame;
    if (_demo) {
      frame_ptr->shell = _demo->shell_mesh();
      frame_ptr->rod = _demo->rod_mesh();
      frame_ptr->rod_polyline = _demo->rod_polyline();
    }

    auto exported = lewitt::logger::geometry::steal_lines();
    frame_ptr->debug_lines.reserve(exported.size());
    for (const auto &line : exported)
      frame_ptr->debug_lines.push_back(to_scene_line(line));

    auto spheres = lewitt::logger::geometry::steal_spheres();
    frame_ptr->debug_spheres.reserve(spheres.size());
    for (const auto &s : spheres) {
      frame_ptr->debug_spheres.push_back(
          {glm::vec3(s.center.x, s.center.y, s.center.z), s.radius,
           glm::vec3(s.color.x, s.color.y, s.color.z)});
    }

    auto tori = lewitt::logger::geometry::steal_tori();
    frame_ptr->debug_tori.reserve(tori.size());
    for (const auto &t : tori) {
      frame_ptr->debug_tori.push_back(
          {glm::vec3(t.center.x, t.center.y, t.center.z),
           glm::vec3(t.axis.x, t.axis.y, t.axis.z), t.major_radius,
           t.minor_radius, glm::vec3(t.color.x, t.color.y, t.color.z)});
    }

    // Rod polyline is display data; bake into debug lines so render stays dump.
    if (frame_ptr->rod_polyline) {
      auto poly = lines_from_rod_snapshot(*frame_ptr->rod_polyline);
      for (const auto &line : poly) {
        frame_ptr->debug_lines.push_back(
            {line.p0, line.p1, line.color, line.radius});
      }
    }
    return frame_ptr;
  }

  void worker_loop() {
    using namespace std::chrono_literals;

    // Publish post-reset geometry before any step so init / first present see
    // the starting scene without waiting on (or skipping past) step 0.
    if (_running.load(std::memory_order_acquire)) {
      geometry_logger::clear();
      const int frame0 =
          _playback ? _playback->sim_frame.load(std::memory_order_relaxed)
                    : _fallback_frame;
      _channel.publish(capture_frame(frame0));
    }

    while (_running.load(std::memory_order_acquire)) {
      // Don't start the next (expensive) step until render has taken the last
      // frame — otherwise the worker busy-loops and starves the main thread.
      if (!_channel.wait_to_produce())
        break;
      if (!_running.load(std::memory_order_acquire))
        break;

      const bool paused = _playback && _playback->paused;
      const bool step_once = _playback && _playback->step_once;
      if (paused && !step_once) {
        std::this_thread::sleep_for(2ms);
        continue;
      }

      const int frame =
          _playback ? _playback->sim_frame.load(std::memory_order_relaxed)
                    : _fallback_frame;

      geometry_logger::clear();
      if (_demo)
        _demo->step(frame);

      _channel.publish(capture_frame(frame));

      if (_playback) {
        _playback->sim_frame.fetch_add(1, std::memory_order_relaxed);
        if (step_once)
          _playback->step_once = false;
      } else {
        ++_fallback_frame;
      }
    }
  }

  duchamp::demo_trait::ptr _demo;
  std::shared_ptr<duchamp_playback> _playback;
  SceneFrameChannel _channel;
  std::atomic<bool> _running{false};
  std::thread _worker;
  int _fallback_frame = 0;
};

} // namespace vermeer
} // namespace gaudi
