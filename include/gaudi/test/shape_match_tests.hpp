#ifndef GAUDI_SHAPE_MATCH_TESTS_HPP
#define GAUDI_SHAPE_MATCH_TESTS_HPP

#include "gaudi/asawa/primitive_objects.hpp"
#include "gaudi/hepworth/block/shape_match.hpp"
#include "gaudi/hepworth/block/solver.hpp"
#include "gaudi/test/test.hpp"

#include <Eigen/Geometry>

#include <cmath>
#include <random>
#include <vector>

namespace gaudi {
namespace test {
namespace {

using hepworth::block::rest_shape_goal;
using hepworth::block::shape_align;
using hepworth::block::shape_match;
using hepworth::block::shape_rows;

inline std::vector<vec3> sm_torus(std::vector<std::array<int, 2>> *edges) {
  std::vector<vec3> x;
  std::vector<std::vector<int>> faces;
  asawa::make_torus(x, faces, 1.0, 0.33, 24, 12);
  if (edges) {
    edges->clear();
    for (const auto &f : faces)
      for (int k = 0; k < 3; ++k) {
        const int a = f[k], b = f[(k + 1) % 3];
        if (a < b)
          edges->push_back({a, b});
        else
          edges->push_back({b, a});
      }
    std::sort(edges->begin(), edges->end());
    edges->erase(std::unique(edges->begin(), edges->end()), edges->end());
  }
  return x;
}

inline std::vector<index_t> sm_all(size_t n) {
  std::vector<index_t> ids(n);
  for (size_t i = 0; i < n; ++i)
    ids[i] = static_cast<index_t>(i);
  return ids;
}

/// Max |A q - p| for one constraint over the block state `x`.
inline real sm_force(std::vector<vec3> &x, shape_rows rows, shape_align align,
                     const std::vector<vec3> &rest,
                     const std::vector<std::array<int, 2>> &edges) {
  std::vector<vec3> vel(x.size(), vec3::Zero()), force = vel;
  std::vector<vec3> mass(x.size(), vec3::Ones());
  auto pos = hepworth::vec3_block::create(mass, x, vel, force);
  vecX q;
  pos->map_to_x(q);
  hepworth::block::shape_match_config cfg;
  cfg.rows = rows;
  cfg.align = align;
  auto c = shape_match::create(sm_all(x.size()), cfg, real(1),
                               std::vector<hepworth::sim_block::ptr>{pos},
                               rest_shape_goal(rest), edges);
  index_t id0 = 0;
  std::vector<trip> trips;
  c->fill_A(id0, trips);
  matS A(id0, q.size());
  A.setFromTriplets(trips.begin(), trips.end());
  vecX p = vecX::Zero(id0);
  c->project(q, p);
  return (A * q - p).cwiseAbs().maxCoeff();
}

GAUDI_TEST(shape_match_rigid_idle) {
  const std::vector<vec3> rest = sm_torus(nullptr);
  const mat3 R =
      Eigen::AngleAxis<real>(0.7, vec3(1, 2, 3).normalized()).toRotationMatrix();
  std::vector<vec3> x(rest.size());
  for (size_t i = 0; i < rest.size(); ++i)
    x[i] = R * rest[i] + vec3(0.3, -1.2, 2.0);
  GAUDI_EXPECT(sm_force(x, shape_rows::pin, shape_align::rigid, rest, {}) <
               1e-10);
  GAUDI_EXPECT(sm_force(x, shape_rows::pin, shape_align::none, rest, {}) >
               1e-2);
}

GAUDI_TEST(shape_match_edge_translation_free) {
  std::vector<std::array<int, 2>> edges;
  const std::vector<vec3> rest = sm_torus(&edges);
  std::vector<vec3> x(rest.size());
  for (size_t i = 0; i < rest.size(); ++i)
    x[i] = rest[i] + vec3(5.0, -2.0, 1.0);
  GAUDI_EXPECT(sm_force(x, shape_rows::edge, shape_align::none, rest, edges) <
               1e-10);
}

GAUDI_TEST(shape_match_restores_rest) {
  const std::vector<vec3> rest = sm_torus(nullptr);
  std::vector<vec3> x = rest;
  std::mt19937 rng(7);
  std::uniform_real_distribution<real> u(-0.03, 0.03);
  const mat3 R =
      Eigen::AngleAxis<real>(0.4, vec3(0, 1, 1).normalized()).toRotationMatrix();
  for (vec3 &p : x)
    p = R * (p + vec3(u(rng), u(rng), u(rng))) + vec3(1, 0, 0);
  std::vector<vec3> vel(x.size(), vec3::Zero()), force = vel;
  std::vector<vec3> mass(x.size(), vec3::Ones());
  auto pos = hepworth::vec3_block::create(mass, x, vel, force);
  std::vector<hepworth::sim_block::ptr> blocks{pos};
  const real h = 0.1;
  hepworth::block::shape_match_config cfg;
  cfg.w = 10; // ten times inertia
  std::vector<hepworth::projection_constraint::ptr> constraints{
      shape_match::create(sm_all(x.size()), cfg, h, blocks,
                          rest_shape_goal(rest))};
  hepworth::block::projection_solver solver;
  solver.set_constraints(constraints);
  for (int s = 0; s < 10; ++s)
    solver.step(blocks, h, real(0.5), 10);
  // Residual after removing the best rigid motion from rest onto x.
  real err = 0;
  {
    std::vector<vec3> g = rest;
    vec3 cx = vec3::Zero(), cg = vec3::Zero();
    for (size_t i = 0; i < x.size(); ++i) {
      cx += x[i];
      cg += g[i];
    }
    cx /= real(x.size());
    cg /= real(x.size());
    mat3 Apq = mat3::Zero();
    for (size_t i = 0; i < x.size(); ++i)
      Apq += (x[i] - cx) * (g[i] - cg).transpose();
    Eigen::JacobiSVD<mat3> svd(Apq, Eigen::ComputeFullU | Eigen::ComputeFullV);
    const mat3 Rf = svd.matrixU() * svd.matrixV().transpose();
    for (size_t i = 0; i < x.size(); ++i)
      err = std::max(err, (Rf * (g[i] - cg) + cx - x[i]).norm());
  }
  console_logger::info << "shape_match restores rest: max err " << err
                       << std::endl;
  GAUDI_EXPECT(err < 1e-3);
}

} // namespace
} // namespace test
} // namespace gaudi

#endif
