#pragma once

#include <atomic>
#include <string>

#include <GLFW/glfw3.h>

namespace gaudi {
namespace vermeer {

// Shared run-harness playback state (Space toggles pause, '.' steps one frame,
// 'O' OBJ dump, 'C' toggles manual vs auto camera framing).
struct duchamp_playback {
  std::atomic<bool> paused{false};
  bool space_was_down = false;
  bool period_was_down = false;
  bool o_was_down = false;
  bool c_was_down = false;
  std::atomic<bool> step_once{false};
  std::atomic<bool> export_obj_once{false};
  std::atomic<bool> toggle_camera_manual{false};
  bool title_shows_paused = false;
  std::atomic<int> sim_frame{0};

  std::string window_title = "Vermeer";
  std::string paused_title;

  void set_window_title(const std::string &title) {
    window_title = title;
    paused_title = title + " [Paused]";
    title_shows_paused = false;
  }
};

inline void poll_playback_input(duchamp_playback &pb, GLFWwindow *window) {
  if (!window)
    return;

  const bool space = glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS;
  if (space && !pb.space_was_down)
    pb.paused = !pb.paused.load();
  pb.space_was_down = space;

  if (pb.paused.load()) {
    const bool period = glfwGetKey(window, GLFW_KEY_PERIOD) == GLFW_PRESS;
    if (period && !pb.period_was_down)
      pb.step_once = true;
    pb.period_was_down = period;
  } else {
    pb.period_was_down = false;
  }

  const bool o_key = glfwGetKey(window, GLFW_KEY_O) == GLFW_PRESS;
  if (o_key && !pb.o_was_down)
    pb.export_obj_once = true;
  pb.o_was_down = o_key;

  const bool c_key = glfwGetKey(window, GLFW_KEY_C) == GLFW_PRESS;
  if (c_key && !pb.c_was_down)
    pb.toggle_camera_manual = true;
  pb.c_was_down = c_key;

  const bool paused_now = pb.paused.load();
  if (paused_now != pb.title_shows_paused) {
    pb.title_shows_paused = paused_now;
    glfwSetWindowTitle(window,
                       paused_now ? pb.paused_title.c_str() : pb.window_title.c_str());
  }
}

} // namespace vermeer
} // namespace gaudi
