#pragma once

#include <algorithm>
#include <cmath>

#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>

#include "gaudi/common.h"
#include "gaudi/vermeer/scene_frame.hpp"
#include "gaudi/vermeer/scene_obb.hpp"
#include "lewitt/camera.hpp"

namespace gaudi {
namespace vermeer {

struct bbox_framing_config {
  bool auto_fit = true;
  /// Capture cadence and min-jerk tween length (sim seconds @ 30fps frames).
  float period_s = 8.0f;
  float fov_y = lewitt::camera::k_fov_y;
  /// Inertia below this is ignored; above treats coast as user control.
  float coast_velocity = 2.0e-4f;
  /// max/min OBB extent. Near 1 ⇒ ball-like; skip retarget.
  float anisotropy_min = 1.35f;
  /// If |dot(view, preferred axis)| ≥ this, current view is good enough.
  float keep_alignment = 0.82f;
  /// Require target score ≥ live score + this (score = align * (aniso − 1)).
  float min_score_improvement = 0.12f;
  /// Live vs expected-tween divergence ⇒ treat as user intervene (scroll, etc.).
  float pose_intervene_eps = 0.06f;
};

struct framing_pose {
  glm::vec3 look_at{0.0f};
  glm::quat rotation{1, 0, 0, 0}; // R * (0,0,1) = look_at → eye
  float distance = 1.0f;
};

inline glm::vec3 to_glm(const vec3 &v) {
  return glm::vec3(static_cast<float>(v.x()), static_cast<float>(v.y()),
                   static_cast<float>(v.z()));
}

inline glm::quat rotation_look_offset(const glm::vec3 &offset_dir,
                                      const glm::vec3 &world_up = glm::vec3(0, 0, 1)) {
  glm::vec3 f = glm::normalize(offset_dir);
  glm::vec3 r = glm::cross(world_up, f);
  if (glm::dot(r, r) < 1.0e-10f)
    r = glm::cross(glm::vec3(0, 1, 0), f);
  r = glm::normalize(r);
  glm::vec3 u = glm::normalize(glm::cross(f, r));
  return glm::normalize(glm::quat_cast(glm::mat3(r, u, f)));
}

inline framing_pose pose_from_camera(const lewitt::camera &cam) {
  framing_pose p;
  p.look_at = cam.look_at();
  p.rotation = rotation_look_offset(cam.offset_direction());
  p.distance = cam.distance();
  return p;
}

inline framing_pose target_pose_from_obb(const oriented_bbox &obb, float aspect,
                                         float fov_y,
                                         const framing_pose &current) {
  framing_pose t;
  t.look_at = to_glm(obb.com);

  glm::vec3 view = to_glm(obb.axes.col(0));
  const glm::vec3 cur_off =
      glm::normalize(current.rotation * glm::vec3(0, 0, 1));
  if (glm::dot(view, cur_off) < 0.0f)
    view = -view;
  t.rotation = rotation_look_offset(view);

  const float r = static_cast<float>(obb.sphere_radius);
  const float tan_y = std::tan(0.5f * fov_y);
  const float tan_x = tan_y * std::max(aspect, 1.0e-3f);
  t.distance = std::max(r / std::max(tan_y, 1.0e-6f),
                        r / std::max(tan_x, 1.0e-6f));
  t.distance = std::max(t.distance, 0.05f);
  return t;
}

inline void apply_pose_to_camera(lewitt::camera &cam, const framing_pose &p) {
  const glm::vec3 off = glm::normalize(p.rotation * glm::vec3(0, 0, 1));
  const lewitt::vec2 angles = lewitt::camera::angles_from_offset_dir(
      lewitt::vec3(off.x, off.y, off.z));
  cam.set_framing(lewitt::vec3(p.look_at.x, p.look_at.y, p.look_at.z), angles,
                  -std::log(std::max(p.distance, 1.0e-4f)));
}

// Flip quat / view so the new sample is nearest the previous pose.
inline void enforce_pose_continuity(framing_pose &pose, const framing_pose &prev) {
  if (glm::dot(pose.rotation, prev.rotation) < 0.0f)
    pose.rotation = -pose.rotation;
  const glm::vec3 a = glm::normalize(pose.rotation * glm::vec3(0, 0, 1));
  const glm::vec3 b = glm::normalize(prev.rotation * glm::vec3(0, 0, 1));
  if (glm::dot(a, b) < 0.0f) {
    pose.rotation = rotation_look_offset(-a);
    if (glm::dot(pose.rotation, prev.rotation) < 0.0f)
      pose.rotation = -pose.rotation;
  }
}

// Classic min-jerk ease on [0,1]: zero vel/accel at ends.
inline float min_jerk_ease(float s) {
  s = glm::clamp(s, 0.0f, 1.0f);
  const float s2 = s * s;
  const float s3 = s2 * s;
  const float s4 = s3 * s;
  const float s5 = s4 * s;
  return 10.0f * s3 - 15.0f * s4 + 6.0f * s5;
}

inline framing_pose lerp_pose(const framing_pose &a, const framing_pose &b,
                              float u) {
  u = glm::clamp(u, 0.0f, 1.0f);
  framing_pose p;
  p.look_at = glm::mix(a.look_at, b.look_at, u);
  glm::quat qb = b.rotation;
  if (glm::dot(a.rotation, qb) < 0.0f)
    qb = -qb;
  p.rotation = glm::normalize(glm::slerp(a.rotation, qb, u));
  const float la = std::log(std::max(a.distance, 1.0e-4f));
  const float lb = std::log(std::max(b.distance, 1.0e-4f));
  p.distance = std::exp(glm::mix(la, lb, u));
  return p;
}

inline float obb_anisotropy(const oriented_bbox &obb) {
  const real emin = std::max(obb.extents.minCoeff(), real(1.0e-8));
  return static_cast<float>(obb.extents.maxCoeff() / emin);
}

inline float view_axis_alignment(const framing_pose &pose,
                                 const oriented_bbox &obb) {
  const glm::vec3 view =
      glm::normalize(pose.rotation * glm::vec3(0.0f, 0.0f, 1.0f));
  const glm::vec3 pref = glm::normalize(to_glm(obb.axes.col(0)));
  return std::abs(glm::dot(view, pref));
}

/// Higher = better framing for an elongated OBB. Ball-like ⇒ ~0.
inline float framing_score(const framing_pose &pose, const oriented_bbox &obb) {
  const float aniso = obb_anisotropy(obb);
  return view_axis_alignment(pose, obb) * std::max(0.0f, aniso - 1.0f);
}

inline float pose_divergence(const framing_pose &a, const framing_pose &b) {
  const float scale = std::max(a.distance, 1.0e-3f);
  const float d_look = glm::length(a.look_at - b.look_at) / scale;
  const float d_rot = 1.0f - std::abs(glm::dot(a.rotation, b.rotation));
  const float d_dist = std::abs(std::log(std::max(a.distance, 1.0e-4f)) -
                                std::log(std::max(b.distance, 1.0e-4f)));
  return d_look + d_rot + d_dist;
}

inline bool should_retarget(const framing_pose &live, const framing_pose &target,
                            const oriented_bbox &obb,
                            const bbox_framing_config &cfg) {
  const float aniso = obb_anisotropy(obb);
  if (aniso < cfg.anisotropy_min)
    return false;
  if (view_axis_alignment(live, obb) >= cfg.keep_alignment)
    return false;
  const float live_s = framing_score(live, obb);
  const float target_s = framing_score(target, obb);
  return target_s >= live_s + cfg.min_score_improvement;
}

// Capture → capture PCA targets with min-jerk eased interpolation.
// User orbit/pan/coast/scroll marks interrupt; on release, ease from the live
// pose toward a new optimal only if it is tunably better (elliptical OBB +
// poor current alignment). Otherwise hold the live view — no snap-back.
class bbox_framing_controller {
public:
  explicit bbox_framing_controller(bbox_framing_config cfg = {}) : _cfg(cfg) {
    if (_cfg.period_s < 1.0e-3f)
      _cfg.period_s = 1.0e-3f;
  }

  bool manual() const { return _manual; }
  void toggle_manual() { _manual = !_manual; }

  // Call once per new SceneFrame. sim_time_s = captured_frame / 30 (nominal 30fps).
  void on_capture(lewitt::camera &cam, const SceneFrame &frame, float sim_time_s) {
    if (!_cfg.auto_fit || _manual)
      return;

    const float vel = glm::length(cam._drag_state.velocity);
    const bool user_control = cam._drag_state.active || cam._pan_state.active ||
                              vel > _cfg.coast_velocity;
    if (user_control) {
      _interrupted = true;
      return;
    }

    const framing_pose live = pose_from_camera(cam);

    // Scroll / external nudge mid-tween: stop fighting the user.
    if (_initialized && !_holding) {
      const float s_chk =
          glm::clamp((sim_time_s - _segment_t0) / _cfg.period_s, 0.0f, 1.0f);
      const framing_pose expected =
          lerp_pose(_from, _to, min_jerk_ease(s_chk));
      if (pose_divergence(live, expected) > _cfg.pose_intervene_eps)
        _interrupted = true;
    }

    const bool resume = _interrupted;
    _interrupted = false;

    const bool due =
        !_initialized || (sim_time_s - _segment_t0 >= _cfg.period_s);
    if (due || resume) {
      const oriented_bbox obb =
          oriented_bbox_from_snapshots(frame.shell, frame.rod);
      if (!obb.ok)
        return;

      framing_pose next =
          target_pose_from_obb(obb, cam.aspect(), _cfg.fov_y, live);
      enforce_pose_continuity(next, live);

      // Always anchor from the live camera (never an old _to) so mouse-up
      // glides instead of snapping back along a stale segment.
      const bool go =
          !_initialized || should_retarget(live, next, obb, _cfg);
      _from = live;
      if (go) {
        _to = next;
        _holding = false;
      } else {
        _to = live;
        _holding = true;
      }
      _segment_t0 = sim_time_s;
      _initialized = true;
    }

    if (!_initialized || _holding)
      return;

    const float s =
        glm::clamp((sim_time_s - _segment_t0) / _cfg.period_s, 0.0f, 1.0f);
    apply_pose_to_camera(cam, lerp_pose(_from, _to, min_jerk_ease(s)));
  }

private:
  bbox_framing_config _cfg;
  bool _manual = false;
  bool _initialized = false;
  bool _interrupted = false;
  bool _holding = false; // live view accepted; wait until retarget is worth it
  float _segment_t0 = 0.0f;
  framing_pose _from;
  framing_pose _to;
};

} // namespace vermeer
} // namespace gaudi
