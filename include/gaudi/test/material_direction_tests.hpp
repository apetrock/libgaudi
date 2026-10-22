#ifndef GAUDI_TEST_MATERIAL_DIRECTION_TESTS_HPP
#define GAUDI_TEST_MATERIAL_DIRECTION_TESTS_HPP

#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/hepworth/block/shape_projection.hpp"
#include "gaudi/kusama/material_direction.hpp"
#include "gaudi/test/test.hpp"

#include <Eigen/Geometry>

#include <cmath>
#include <random>
#include <vector>

namespace gaudi {
namespace test {

namespace material_detail {

/// Unit hoop direction about the z axis at each vertex (tangent to the
/// built-in torus).
inline std::vector<vec3> hoop(const std::vector<vec3> &x) {
  std::vector<vec3> u(x.size());
  for (size_t i = 0; i < x.size(); ++i)
    u[i] = vec3(-x[i].y(), x[i].x(), 0).normalized();
  return u;
}

inline real mean_angle(const std::vector<vec3> &a, const std::vector<vec3> &b) {
  real s = 0;
  for (size_t i = 0; i < a.size(); ++i)
    s += std::acos(std::clamp(std::abs(a[i].dot(b[i])), real(0), real(1)));
  return s / real(a.size());
}

} // namespace material_detail

/// `directed_stretch` given the world axis `chart_directions` reads from
/// `angle` matches the row reading `angle` itself (with curvature scaling),
/// either sign of the axis.
GAUDI_TEST(stretch_matches_its_angle_given_the_direction) {
  const std::vector<vec3> spokes{vec3(2, 0, real(0.4)), vec3(-2, 0, 0), vec3(0, 1, 0),
                                 vec3(0, -1, 0), vec3(1, 1, 0)};
  const std::vector<real> weights(spokes.size(), real(1));
  std::vector<vec3> rest_x{vec3::Zero()};
  rest_x.insert(rest_x.end(), spokes.begin(), spokes.end());
  // A turned live ring, so the polar is not the identity.
  std::vector<vec3> live_x = rest_x;
  const mat3 turn = Eigen::AngleAxis<real>(real(0.3), vec3(1, 2, 3).normalized()).matrix();
  for (vec3 &p : live_x)
    p = turn * p;
  std::vector<vec3> mass(rest_x.size(), vec3::Ones()), vel(rest_x.size(), vec3::Zero()),
      force(rest_x.size(), vec3::Zero());
  auto pos = hepworth::vec3_block::create(mass, live_x, vel, force);
  const std::vector<hepworth::sim_block::ptr> blocks{pos};
  std::vector<hepworth::block::vertex_chart> charts(1);
  charts[0].ok = true;
  charts[0].T = vec3::UnitX();
  charts[0].B = vec3::UnitY();
  charts[0].N = vec3::UnitZ();
  std::vector<index_t> ids(rest_x.size());
  for (size_t i = 0; i < ids.size(); ++i)
    ids[i] = static_cast<index_t>(i);
  vecX q;
  pos->map_to_x(q);
  auto projected = [&](auto row) {
    index_t id0 = 0;
    std::vector<trip> trips;
    row->fill_A(id0, trips);
    vecX p = vecX::Zero(id0);
    row->project(q, p);
    return row->rest();
  };
  auto stretch = [&](real angle, const vec3 &dir) {
    return projected(hepworth::block::directed_stretch::create(
        ids, weights, rest_x, hepworth::block::fan_local{}, 0, angle, real(0.2), real(-0.05),
        real(0.5), real(1), blocks, &charts, dir));
  };
  // k_hi the larger eigenvalue, then the smaller (negative curvatures).
  for (const vec3 &h : {vec3(real(3), real(0.4), real(1)), vec3(real(-3), real(0.4), real(-1))})
    for (const real angle : {real(0), real(0.4), real(-1.1)}) {
      charts[0].h = h;
      const vec3 dir = hepworth::block::chart_directions(charts, angle)[0];
      GAUDI_ASSERT(std::abs(dir.norm() - real(1)) < real(1e-12));
      const std::vector<vec3> a = stretch(angle, vec3::Zero());
      const std::vector<vec3> b = stretch(real(0), dir);
      const std::vector<vec3> c = stretch(real(0), -dir);
      GAUDI_ASSERT(a.size() == b.size() && a.size() == c.size());
      for (size_t k = 0; k < a.size(); ++k) {
        GAUDI_EXPECT((a[k] - b[k]).norm() < real(1e-10));
        GAUDI_EXPECT((a[k] - c[k]).norm() < real(1e-10));
      }
    }
}

GAUDI_TEST(material_direction_turns_with_a_rigid_motion) {
  asawa::shell::shell::ptr M = asawa::shell::load_torus(1.0, 0.33, 32, 16);
  const std::vector<vec3> x0 = asawa::get_vec_data(*M, 0);
  const mat3 R = Eigen::AngleAxis<real>(real(0.7), vec3(1, -2, 0.5).normalized()).matrix();
  std::vector<vec3> x1 = x0;
  for (vec3 &p : x1)
    p = R * p + vec3(0.3, -1, 2);
  // Tangent to the mesh's own vertex normals, which tilt slightly off the
  // smooth torus's.
  std::vector<vec3> u0 = material_detail::hoop(x0);
  for (size_t i = 0; i < u0.size(); ++i)
    GAUDI_ASSERT(kusama::material_direction::tangent_unit(
        u0[i], kusama::material_direction::unit_vert_normal(*M, static_cast<int>(i), x0)));
  std::vector<vec3> u = u0;
  kusama::material_direction::transport(*M, x0, x1, u);
  real err = 0;
  for (size_t i = 0; i < u.size(); ++i)
    err = std::max(err, (u[i] - R * u0[i]).norm());
  std::cout << "  rigid transport error " << err << std::endl;
  GAUDI_EXPECT(err < real(1e-9));
}

GAUDI_TEST(material_direction_orient_agrees_in_sign) {
  asawa::shell::shell::ptr M = asawa::shell::load_torus(1.0, 0.33, 32, 16);
  const std::vector<vec3> &x = asawa::get_vec_data(*M, 0);
  const std::vector<vec3> u0 = material_detail::hoop(x);
  std::vector<vec3> u = u0;
  std::mt19937 rng(7);
  for (vec3 &v : u)
    if (rng() & 1u)
      v = -v;
  kusama::material_direction::orient(*M, u);
  const real s = u[0].dot(u0[0]) > 0 ? real(1) : real(-1);
  int bad = 0;
  for (size_t i = 0; i < u.size(); ++i)
    bad += (s * u[i] - u0[i]).norm() > real(1e-12);
  GAUDI_EXPECT(bad == 0);
}

GAUDI_TEST(material_direction_lean_is_sign_free) {
  const vec3 a = vec3::UnitX(), b = vec3::UnitY();
  for (const real s : {real(1), real(-1)}) {
    std::vector<vec3> u{a}, g{s * b};
    kusama::material_direction::lean_toward(u, g, real(0));
    GAUDI_EXPECT((u[0] - a).norm() < real(1e-15));
    kusama::material_direction::lean_toward(u, g, real(0.5));
    GAUDI_EXPECT((u[0] - (a + s * b).normalized()).norm() < real(1e-12));
    u = {a};
    kusama::material_direction::lean_toward(u, g, real(1));
    GAUDI_EXPECT((u[0] - s * b).norm() < real(1e-12));
  }
}

GAUDI_TEST(material_direction_smooth_removes_noise_and_keeps_the_field) {
  asawa::shell::shell::ptr M = asawa::shell::load_torus(1.0, 0.33, 32, 16);
  const std::vector<vec3> &x = asawa::get_vec_data(*M, 0);
  const std::vector<vec3> u0 = material_detail::hoop(x);
  const real l = asawa::shell::avg_length(*M, x);
  const real t = std::pow(real(1.5) * l, 2);

  std::vector<vec3> clean = u0;
  kusama::material_direction::smooth(*M, x, clean, t);
  const real drift = material_detail::mean_angle(clean, u0);
  std::cout << "  clean drift " << drift << " rad" << std::endl;
  GAUDI_EXPECT(drift < real(0.05));

  std::vector<vec3> noisy = u0;
  std::mt19937 rng(3);
  std::uniform_real_distribution<real> pick(real(-0.5), real(0.5));
  for (size_t i = 0; i < noisy.size(); ++i) {
    const vec3 n = kusama::material_direction::unit_vert_normal(*M, static_cast<int>(i), x);
    noisy[i] = Eigen::AngleAxis<real>(pick(rng), n).matrix() * u0[i];
  }
  const real before = material_detail::mean_angle(noisy, u0);
  kusama::material_direction::smooth(*M, x, noisy, t);
  const real after = material_detail::mean_angle(noisy, u0);
  std::cout << "  noise " << before << " -> " << after << " rad" << std::endl;
  // Independent noise per vertex; a heat length of 1.5 edges averages about
  // seven vertices, so roughly a factor sqrt(7).
  GAUDI_EXPECT(after < real(0.5) * before);
  for (const vec3 &v : noisy)
    GAUDI_EXPECT(std::abs(v.norm() - real(1)) < real(1e-9));
}

} // namespace test
} // namespace gaudi

#endif
