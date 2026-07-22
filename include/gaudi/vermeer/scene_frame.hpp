#pragma once

#include <atomic>
#include <condition_variable>
#include <memory>
#include <mutex>
#include <optional>
#include <utility>
#include <vector>

#include <glm/glm.hpp>

#include "gaudi/duchamp/demo_trait.hpp"

namespace gaudi {
namespace vermeer {

struct scene_debug_line {
  glm::vec3 p0{0.0f};
  glm::vec3 p1{0.0f};
  glm::vec3 color{1.0f};
  float radius = 0.01f;
};

// CPU snapshot published by the sim worker; render thread uploads only this.
struct SceneFrame {
  int sim_frame = 0;
  std::optional<duchamp::mesh_snapshot> shell;
  std::optional<duchamp::mesh_snapshot> rod;
  std::optional<duchamp::rod_snapshot> rod_polyline;
  std::vector<scene_debug_line> debug_lines;
};

// Single-slot handoff: sim waits for a free slot, then publishes; render take
// clears the slot. Keeps the worker from busy-looping ahead of the GPU.
class SceneFrameChannel {
public:
  // Block until the previous frame was taken (or shutdown). Call before step().
  bool wait_to_produce() {
    std::unique_lock<std::mutex> lock(_mutex);
    _cv_consumed.wait(lock, [this] {
      return !_pending.load(std::memory_order_relaxed) ||
             _shutdown.load(std::memory_order_relaxed);
    });
    return !_shutdown.load(std::memory_order_relaxed);
  }

  void publish(std::shared_ptr<const SceneFrame> frame) {
    std::unique_lock<std::mutex> lock(_mutex);
    if (_shutdown.load(std::memory_order_relaxed))
      return;
    _latest = std::move(frame);
    _pending.store(true, std::memory_order_relaxed);
    _generation.fetch_add(1, std::memory_order_relaxed);
    lock.unlock();
    _cv_produced.notify_one();
  }

  // Returns the newest frame if newer than seen_generation; clears pending so
  // the producer can run again.
  std::shared_ptr<const SceneFrame> take_latest(uint64_t &seen_generation) {
    std::lock_guard<std::mutex> lock(_mutex);
    const uint64_t gen = _generation.load(std::memory_order_relaxed);
    if (!_latest || gen == seen_generation)
      return nullptr;
    seen_generation = gen;
    _pending.store(false, std::memory_order_relaxed);
    auto out = _latest;
    _cv_consumed.notify_one();
    return out;
  }

  void shutdown() {
    {
      std::lock_guard<std::mutex> lock(_mutex);
      _shutdown.store(true, std::memory_order_relaxed);
      _pending.store(false, std::memory_order_relaxed);
    }
    _cv_consumed.notify_all();
    _cv_produced.notify_all();
  }

private:
  std::mutex _mutex;
  std::condition_variable _cv_consumed;
  std::condition_variable _cv_produced;
  std::shared_ptr<const SceneFrame> _latest;
  std::atomic<uint64_t> _generation{0};
  std::atomic<bool> _pending{false};
  std::atomic<bool> _shutdown{false};
};

} // namespace vermeer
} // namespace gaudi
