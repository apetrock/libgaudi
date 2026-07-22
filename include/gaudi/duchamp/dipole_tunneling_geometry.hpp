#pragma once

#include <array>
#include <cmath>
#include <limits>
#include <set>
#include <tuple>
#include <vector>

#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/common.h"
#include "gaudi/geometry_types.hpp"
#include "gaudi/vec_addendum.h"

namespace gaudi {
namespace duchamp {
namespace dipole_tunneling {

enum class dipole_zone { outside, upper_clearance, lower_contact };

// (xr0, xr1, N, T, B)
using segment_frame = std::tuple<vec3, vec3, vec3, vec3, vec3>;

inline constexpr real k_dipole_eps2 = 1e-24;

inline bool is_finite_vec(const vec3 &v) {
  return std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]);
}

inline vec3 safe_normalized(const vec3 &v, const vec3 &fallback = vec3::UnitY()) {
  const real n2 = v.squaredNorm();
  if (!(n2 > k_dipole_eps2) || !is_finite_vec(v))
    return fallback;
  return v / std::sqrt(n2);
}

// Foot on segment; collapsed edges return v0 (no /0).
inline vec3 safe_project_on_line(const vec3 &v0, const vec3 &v1, const vec3 &pt) {
  const vec3 dx = v1 - v0;
  const real dx2 = dx.squaredNorm();
  if (!(dx2 > k_dipole_eps2))
    return v0;
  real s = (pt - v0).dot(dx) / dx2;
  s = va::clamp(s, 0.0, 1.0);
  return v0 + s * dx;
}

// Inf when dp ⟂ N or dp≈0 so gather/force radius gates reject cleanly.
inline real safe_tangent_point_radius(const vec3 &dp, const vec3 &N) {
  const real nPdp = (N * N.transpose() * dp).norm();
  if (!(nPdp > 1e-12) || !std::isfinite(nPdp))
    return std::numeric_limits<real>::infinity();
  const real R = 0.5 * dp.squaredNorm() / nPdp;
  return std::isfinite(R) ? R : std::numeric_limits<real>::infinity();
}

// Slerp rod normals along a segment: R(t) = I.slerp(t, FromTwoVectors(N0,N1)).
inline vec3 interpolate_rod_normal(const vec3 &N0, const vec3 &N1, real t) {
  const vec3 n0 = safe_normalized(N0, vec3::Zero());
  const vec3 n1 = safe_normalized(N1, vec3::Zero());
  if (n0.squaredNorm() < k_dipole_eps2 && n1.squaredNorm() < k_dipole_eps2)
    return vec3::UnitY();
  if (n0.squaredNorm() < k_dipole_eps2)
    return n1;
  if (n1.squaredNorm() < k_dipole_eps2 || (n0 - n1).squaredNorm() < k_dipole_eps2)
    return n0;
  if ((n0 + n1).squaredNorm() < k_dipole_eps2)
    return n0; // antipodal: slerp undefined; keep N0
  const quat q = quat::Identity().slerp(t, quat::FromTwoVectors(n0, n1));
  const vec3 out = q * n0;
  return safe_normalized(out, n0);
}

// Returns false if the rod edge is collapsed (caller should skip).
inline bool try_make_segment_frame(const vec3 &xr0, const vec3 &xr1,
                                   const vec3 &Ng, segment_frame &out) {
  const vec3 t_raw = xr1 - xr0;
  if (!(t_raw.squaredNorm() > k_dipole_eps2) || !is_finite_vec(t_raw))
    return false;
  const vec3 t = t_raw.normalized();
  const vec3 ng = safe_normalized(Ng, t.unitOrthogonal());
  vec3 b = ng.cross(t);
  if (!(b.squaredNorm() > k_dipole_eps2))
    b = t.unitOrthogonal(); // Nr ∥ T
  else
    b.normalize();
  const vec3 n = safe_normalized(t.cross(b), ng);
  if (!is_finite_vec(n) || !is_finite_vec(t) || !is_finite_vec(b))
    return false;
  out = {xr0, xr1, n, t, b};
  return true;
}

inline segment_frame make_segment_frame(const vec3 &xr0, const vec3 &xr1,
                                        const vec3 &Ng) {
  segment_frame f;
  if (try_make_segment_frame(xr0, xr1, Ng, f))
    return f;
  // Degenerate edge: stable orthonormal placeholder (N from Ng if possible).
  const vec3 t = vec3::UnitX();
  const vec3 n = safe_normalized(Ng, vec3::UnitY());
  const vec3 b = safe_normalized(n.cross(t), t.unitOrthogonal());
  const vec3 nn = safe_normalized(t.cross(b), n);
  return {xr0, xr1, nn, t, b};
}

inline segment_frame make_segment_frame(const vec3 &xr0, const vec3 &xr1,
                                        const vec3 &N0, const vec3 &N1,
                                        real t = 0.5) {
  return make_segment_frame(xr0, xr1, interpolate_rod_normal(N0, N1, t));
}

inline segment_frame make_segment_frame(index_t sid,
                                        const std::vector<vec3> &x,
                                        const vec3 &N) {
  return make_segment_frame(x[sid], x[sid + 1], N);
}

inline segment_frame make_segment_frame(index_t sid,
                                        const std::vector<vec3> &x,
                                        const std::vector<vec3> &Nr,
                                        real t = 0.5) {
  return make_segment_frame(x[sid], x[sid + 1], Nr[sid], Nr[sid + 1], t);
}

inline vec2 project_to_NB(const vec3 &dx, const vec3 &N, const vec3 &B) {
  return vec2(dx.dot(B), dx.dot(N));
}

inline real signed_disk_distance(const vec2 &pt, const vec2 &cen, real radius) {
  return (pt - cen).norm() - radius;
}

inline bool inside_disk(const vec2 &pt, const vec2 &cen, real radius) {
  return signed_disk_distance(pt, cen, radius) < 0.0;
}

inline dipole_zone classify_dipole_disk(const vec3 &x, const vec3 &xr,
                                        const segment_frame &f, real h,
                                        real r) {
  const auto &[xr0, xr1, N, T, B] = f;
  (void)xr0;
  (void)xr1;
  (void)T;
  const vec2 pt = project_to_NB(x - xr, N, B);
  if (inside_disk(pt, vec2(0.0, h), r))
    return dipole_zone::upper_clearance;
  if (inside_disk(pt, vec2(0.0, -h), r))
    return dipole_zone::lower_contact;
  return dipole_zone::outside;
}

inline int tunnel_sign(dipole_zone zone) {
  return zone == dipole_zone::upper_clearance ? 1 : -1;
}

inline std::array<vec3, 4> cap_disc_corners(const vec3 &p, const segment_frame &f,
                                            real r, int sign) {
  const auto &[xr0, xr1, N, T, B] = f;
  (void)xr0;
  (void)xr1;
  (void)T;
  const real s = real(sign);
  const vec3 C = p + s * r * N;
  return {C + r * B, C - r * B, p, p + 2.0 * s * r * N};
}

inline std::array<vec3, 8>
make_tangent_tunnel_corners(const vec3 &p0, const vec3 &p1,
                            const segment_frame &f, real r, dipole_zone zone) {
  const int sign = tunnel_sign(zone);
  const auto c0 = cap_disc_corners(p0, f, r, sign);
  const auto c1 = cap_disc_corners(p1, f, r, sign);
  return {c0[0], c0[1], c0[2], c0[3], c1[0], c1[1], c1[2], c1[3]};
}

inline ext::extents_t extents_from_corners(const std::array<vec3, 8> &corners) {
  return ext::calc_extents(corners);
}

inline ext::extents_t make_tangent_tunnel_aabb(const vec3 &p0, const vec3 &p1,
                                               const segment_frame &f, real r,
                                               dipole_zone zone) {
  return extents_from_corners(make_tangent_tunnel_corners(p0, p1, f, r, zone));
}

inline std::array<ext::extents_t, 2>
make_dual_tunnel_aabbs(const segment_frame &f, real r) {
  const auto &[xr0, xr1, N, T, B] = f;
  (void)N;
  (void)T;
  (void)B;
  return {make_tangent_tunnel_aabb(xr0, xr1, f, r, dipole_zone::upper_clearance),
          make_tangent_tunnel_aabb(xr0, xr1, f, r, dipole_zone::lower_contact)};
}

inline bool inside_extents(const ext::extents_t &box, const vec3 &p,
                           real eps = 1e-12) {
  return ext::dist(box, p) <= eps;
}

inline real dipole_margin(real r) { return 0.35 * r; }

// Match triangle_dipole_tunneling::calc_proj_N (#if 0 = tangent-point R, else fixed r).
inline vec3 dipole_sphere_center(const vec3 &xr, const vec3 &dp, const vec3 &N,
                                 real r) {
  const real sg = va::sgn(N.dot(dp));
  const real s = (sg == 0.0) ? 1.0 : sg;
#if 0
  const real R = safe_tangent_point_radius(dp, N);
  if (!std::isfinite(R))
    return xr + s * r * N;
  return xr + s * R * N;
#else
  return xr + s * r * N;
#endif
}

// Project p onto the dipole sphere: cen + R * (p - cen).normalized().
// Side from live geometry: sgn(N·(p - xr)), same as dipole_sphere_center.
inline vec3 proj_to_dipole(const vec3 &p_mesh, const vec3 &xr, const vec3 &Nr,
                           real R) {
  if (!(R > 0.0) || !std::isfinite(R) || !is_finite_vec(p_mesh) ||
      !is_finite_vec(xr))
    return p_mesh;
  const vec3 N = safe_normalized(Nr);
  const vec3 dp = p_mesh - xr;
  const real sg = va::sgn(N.dot(dp));
  const real s = (sg == 0.0) ? 1.0 : sg;
  const vec3 cen = xr + s * R * N;
  const vec3 dpcen = p_mesh - cen;
  const vec3 n =
      dpcen.squaredNorm() > k_dipole_eps2 ? safe_normalized(dpcen, s * N) : s * N;
  const vec3 out = cen + R * n;
  return is_finite_vec(out) ? out : p_mesh;
}

inline vec3 proj_to_dipole(const vec3 &p_mesh, const vec3 &xr,
                           const segment_frame &f, real R) {
  const auto &[xr0, xr1, N, T, B] = f;
  (void)xr0;
  (void)xr1;
  (void)T;
  (void)B;
  return proj_to_dipole(p_mesh, xr, N, R);
}

inline std::pair<std::vector<asawa::shell::FaceId>, std::vector<asawa::shell::FaceId>>
gather_faces_in_segment(const asawa::shell::shell &M,
                        const std::vector<vec3> &x,
                        const segment_frame &f, real r) {
  const real margin = dipole_margin(r);
  const real r_gather = r + margin;
  const auto &[upper_box, lower_box] = make_dual_tunnel_aabbs(f, r_gather);
  const auto &[xr0, xr1, N, T, B] = f;
  (void)T;
  (void)B;

  std::vector<asawa::shell::FaceId> upper_out;
  std::vector<asawa::shell::FaceId> lower_out;
  std::set<int> upper_seen;
  std::set<int> lower_seen;

  auto try_add = [&](asawa::shell::FaceId face, bool hit_upper, bool hit_lower) {
    const vec3 xf = asawa::shell::face_center(M, face, x);
    const vec3 xr = safe_project_on_line(xr0, xr1, xf);
    const vec3 dp = xf - xr;
    const real R = safe_tangent_point_radius(dp, N);
    if (!(R <= r_gather + 1e-9))
      return;
    const int fid = static_cast<int>(face);
    if (hit_upper && N.dot(dp) >= 0.0 && upper_seen.insert(fid).second)
      upper_out.emplace_back(face);
    if (hit_lower && N.dot(dp) < 0.0 && lower_seen.insert(fid).second)
      lower_out.emplace_back(face);
  };

  for (const auto fi : M.get_face_range()) {
    const asawa::shell::FaceId face = asawa::shell::face_id(fi);
    bool hit_upper = false;
    bool hit_lower = false;
    M.const_for_each_face_tri(
        face, [&](asawa::shell::CornerId c0, asawa::shell::CornerId c1,
                  asawa::shell::CornerId c2, const asawa::shell::shell &) {
          if (hit_upper && hit_lower)
            return;
          const vec3 p0 = x[M.vert(c0)];
          const vec3 p1 = x[M.vert(c1)];
          const vec3 p2 = x[M.vert(c2)];
          if (!hit_upper && ext::overlap(upper_box, p0, p1, p2))
            hit_upper = true;
          if (!hit_lower && ext::overlap(lower_box, p0, p1, p2))
            hit_lower = true;
        });
    if (!hit_upper && !hit_lower)
      continue;
    try_add(face, hit_upper, hit_lower);
  }
  return {upper_out, lower_out};
}

// Dipole classify for a prefiltered face (no full-mesh AABB sweep).
// Returns +1 upper, -1 lower, 0 outside.
inline int classify_face_on_segment(const asawa::shell::shell &M,
                                    asawa::shell::FaceId face,
                                    const std::vector<vec3> &x,
                                    const segment_frame &f, real r) {
  const real r_gather = r + dipole_margin(r);
  const auto &[xr0, xr1, N, T, B] = f;
  (void)T;
  (void)B;
  const vec3 xf = asawa::shell::face_center(M, face, x);
  const vec3 xr = safe_project_on_line(xr0, xr1, xf);
  const vec3 dp = xf - xr;
  const real R = safe_tangent_point_radius(dp, N);
  if (!(R <= r_gather + 1e-9))
    return 0;
  return N.dot(dp) >= 0.0 ? 1 : -1;
}

// Project p onto the dipole cylinder by an SDF step along the mesh normal:
//   side = sgn(N_rod · (p − xr))   (+1 upper, −1 lower)
//   p'   = p + (−side) · (R − dist) · N_mesh
// Upper half-plane → negative mesh-N; lower → positive mesh-N.
// Cylinder is rod-relative: axis through dipole_sphere_center along T.
inline vec3 proj_to_dipole_cylinder_along_normal(const vec3 &p_mesh,
                                                const vec3 &xr,
                                                const segment_frame &f, real R,
                                                const vec3 &N_mesh) {
  if (!(R > 0.0) || !std::isfinite(R) || !is_finite_vec(p_mesh) ||
      !is_finite_vec(xr))
    return p_mesh;
  const auto &[xr0, xr1, N, T, B] = f;
  (void)xr0;
  (void)xr1;
  (void)B;
  if (!is_finite_vec(N) || !is_finite_vec(T))
    return p_mesh;

  const vec3 dp = p_mesh - xr;
  const real sg = va::sgn(N.dot(dp));
  const real side = (sg == 0.0) ? 1.0 : sg;

  const vec3 cen = dipole_sphere_center(xr, dp, N, R);
  if (!is_finite_vec(cen))
    return p_mesh;

  const vec3 dpcen = p_mesh - cen;
  const vec3 radial = dpcen - dpcen.dot(T) * T;
  const real dist = radial.norm();
  if (!std::isfinite(dist) || !(dist > k_dipole_eps2))
    return p_mesh;

  const vec3 n = safe_normalized(N_mesh, vec3::Zero());
  if (!(n.squaredNorm() > k_dipole_eps2))
    return p_mesh;

  const vec3 out = p_mesh + (-side) * (R - dist) * n;
  return is_finite_vec(out) ? out : p_mesh;
}

// Displacement along the mesh normal onto the dipole tunnel cylinder.
// Same target as proj_to_dipole_cylinder_along_normal; returns dx = p_proj - x.
// Inside → push out; outside (within margin) → pull in. Caller scales by
// w/(dt*dt) so w=1 lands on the cylinder under x += dt*v + dt*dt*f.
inline vec3 dipole_radial_force_at(const vec3 &x, const segment_frame &f, real r,
                                  real /*disk_h*/, bool /*above*/, real /*C*/,
                                  const vec3 &N_tri) {
  const auto &[xr0, xr1, N, T, B] = f;
  (void)B;
  if (!(r > 0.0) || !std::isfinite(r) || !is_finite_vec(x) || !is_finite_vec(N) ||
      !is_finite_vec(T))
    return vec3::Zero();
  if ((xr1 - xr0).squaredNorm() <= k_dipole_eps2)
    return vec3::Zero();

  const real margin = dipole_margin(r);
  const vec3 xr = safe_project_on_line(xr0, xr1, x);
  const vec3 dp = x - xr;
  const real Rtp = safe_tangent_point_radius(dp, N);
  if (!(Rtp <= r + margin + 1e-9))
    return vec3::Zero();

  const vec3 cen = dipole_sphere_center(xr, dp, N, r);
  if (!is_finite_vec(cen))
    return vec3::Zero();
  const vec3 dpcen = x - cen;
  const real dist = (dpcen - dpcen.dot(T) * T).norm();
  if (!std::isfinite(dist) || !(dist > k_dipole_eps2))
    return vec3::Zero();
  if (dist > r + margin + 1e-9)
    return vec3::Zero();

  const vec3 p_proj =
      proj_to_dipole_cylinder_along_normal(x, xr, f, r, N_tri);
  const vec3 out = p_proj - x;
  return (out.squaredNorm() > 0.0 && is_finite_vec(out)) ? out : vec3::Zero();
}

inline vec3 dipole_face_radial_force(const vec3 &xf, const segment_frame &f, real r,
                                     real h, bool above, real C, const vec3 &N_tri) {
  return dipole_radial_force_at(xf, f, r, h, above, C, N_tri);
}

// edge_verts: flat pairs (v0,v1,...) from rod::get_edge_vert_ids().
// Do NOT iterate x_r[i]/x_r[i+1] — after remesh next(i) may not be i+1, and
// that creates phantom chords through the knot (force "trails").
// force_w is unused (kept for API); scale displacements by w/(dt*dt) at the
// call site so w=1 lands on the dipole cylinder under x += dt*v + dt*dt*f.
inline void accumulate_dipole_tunnel_forces(
    const asawa::shell::shell &M, const std::vector<vec3> &x,
    const std::vector<vec3> &x_r, const std::vector<vec3> &Nr,
    const std::vector<index_t> &edge_verts, real r, real h, real force_w,
    std::vector<vec3> &forces_out) {
  (void)force_w;
  forces_out.assign(M.vert_count(), vec3::Zero());
  assert(Nr.size() >= x_r.size());
  assert(edge_verts.size() % 2 == 0);

  for (size_t e = 0; e + 1 < edge_verts.size(); e += 2) {
    const index_t v0 = edge_verts[e];
    const index_t v1 = edge_verts[e + 1];
    if (v0 < 0 || v1 < 0 || static_cast<size_t>(v0) >= x_r.size() ||
        static_cast<size_t>(v1) >= x_r.size() ||
        static_cast<size_t>(v0) >= Nr.size() ||
        static_cast<size_t>(v1) >= Nr.size())
      continue;
    segment_frame f;
    if (!try_make_segment_frame(x_r[v0], x_r[v1],
                                interpolate_rod_normal(Nr[v0], Nr[v1], 0.5), f))
      continue;
    const auto [upper, lower] = gather_faces_in_segment(M, x, f, r);

    auto add_faces = [&](const std::vector<asawa::shell::FaceId> &faces,
                         bool above) {
      for (const asawa::shell::FaceId fid : faces) {
        if (M.fsize(fid) != 3)
          continue;
        const auto tri = M.get_tri(fid);
        for (index_t vi : tri) {
          const vec3 N_tri =
              asawa::shell::vert_normal(M, asawa::shell::vert_id(vi), x);
          const real s = [&]() {
            const vec3 &r0 = x_r[v0];
            const vec3 &r1 = x_r[v1];
            const vec3 dr = r1 - r0;
            const real dr2 = dr.squaredNorm();
            real t = dr2 > k_dipole_eps2 ? (x[vi] - r0).dot(dr) / dr2 : 0.5;
            return va::clamp(t, 0.0, 1.0);
          }();
          segment_frame fi;
          if (!try_make_segment_frame(
                  x_r[v0], x_r[v1],
                  interpolate_rod_normal(Nr[v0], Nr[v1], s), fi))
            continue;
          const vec3 fv =
              dipole_radial_force_at(x[vi], fi, r, h, above, /*C=*/1.0, N_tri);
          if (fv.squaredNorm() > 0.0 && is_finite_vec(fv))
            forces_out[vi] += fv;
        }
      }
    };

    add_faces(upper, true);
    add_faces(lower, false);
  }
}

inline void accumulate_dipole_tunnel_forces(
    const asawa::shell::shell &M, const std::vector<vec3> &x,
    const std::vector<vec3> &x_r, const vec3 &dipole_N,
    const std::vector<index_t> &edge_verts, real r, real h, real force_w,
    std::vector<vec3> &forces_out) {
  accumulate_dipole_tunnel_forces(
      M, x, x_r, std::vector<vec3>(x_r.size(), dipole_N), edge_verts, r, h,
      force_w, forces_out);
}

} // namespace dipole_tunneling
} // namespace duchamp
} // namespace gaudi
