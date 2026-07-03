#ifndef __GAUDI_TEST_FIELD_TESTS_HPP__
#define __GAUDI_TEST_FIELD_TESTS_HPP__

#include <memory>
#include <vector>

#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/duchamp/field_nodes.hpp"
#include "gaudi/duchamp/fields.hpp"
#include "gaudi/test/test.hpp"

#include "liblombardi/graph_context.hpp"

namespace gaudi {
namespace test {

GAUDI_TEST(field_bound_access) {
  auto M = asawa::shell::load_cube(); // 8 verts, positions at datum 0
  auto xs = std::make_shared<duchamp::shell_vert_positions>(M, 0);

  GAUDI_ASSERT(xs->size() == 8);

  // bound_field reads/writes mesh-owned storage in place.
  const vec3 original = (*xs)[0];
  (*xs)[0] = vec3(7.0, 8.0, 9.0);
  GAUDI_ASSERT(((*xs)[0] - vec3(7.0, 8.0, 9.0)).norm() < 1e-12);

  // It is a view onto the mesh datum, not a copy.
  const std::vector<vec3> &raw = asawa::const_get_vec_data(*M, 0);
  GAUDI_ASSERT((raw[0] - vec3(7.0, 8.0, 9.0)).norm() < 1e-12);
  (*xs)[0] = original;
}

GAUDI_TEST(field_surface_following_validation) {
  auto M = asawa::shell::load_cube(); // 8 verts
  auto params = std::make_shared<duchamp::shell_vert_normals>(M);

  params->resize(8);
  GAUDI_ASSERT(params->valid());

  params->resize(7);
  GAUDI_ASSERT(!params->valid());
}

GAUDI_TEST(field_polymorphic_algorithm) {
  auto M = asawa::shell::load_cube();
  auto xs = std::make_shared<duchamp::shell_vert_positions>(M, 0);

  auto ridges = std::make_shared<duchamp::shell_ridge_points>(M);
  ridges->resize(3);
  (*ridges)[0] = vec3(1, 0, 0);
  (*ridges)[1] = vec3(0, 1, 0);
  (*ridges)[2] = vec3(0, 0, 1);

  // Same algorithm consumes a bound field and a sampled field via field_base.
  vec3 avg_xs =
      duchamp::compute_average<duchamp::shell_mesh, vec3,
                               asawa::prim_type::VERTEX>(xs);
  vec3 avg_ridges =
      duchamp::compute_average<duchamp::shell_mesh, vec3,
                               asawa::prim_type::VERTEX>(ridges);

  GAUDI_ASSERT((avg_ridges - vec3(1, 1, 1) / 3.0).norm() < 1e-12);
  // cube centroid is (0.5, 0.5, 0.5) for the unit-cube fixture.
  GAUDI_ASSERT((avg_xs - vec3(0.5, 0.5, 0.5)).norm() < 1e-12);
}

GAUDI_TEST(field_sampled_independence) {
  auto ridges = std::make_shared<duchamp::shell_ridge_points>(nullptr);
  ridges->resize(10);

  GAUDI_ASSERT(ridges->size() == 10);
  GAUDI_ASSERT(ridges->mesh == nullptr);
}

GAUDI_TEST(field_march_node_one_step) {
  auto M = asawa::shell::load_sphere(1.0, 24, 16);
  const std::vector<vec3> &x = asawa::const_get_vec_data(*M, 0);

  auto normals = std::make_shared<duchamp::shell_vert_normals>(M);
  normals->get() = asawa::shell::vertex_normals(*M, x);
  GAUDI_ASSERT(normals->valid());

  const real step = 0.1;
  liblombardi::GraphContext ctx;
  auto node = ctx.create_node<duchamp::march_node>(normals, step);

  // Seed the single input port with the sphere positions.
  node->get_datum<duchamp::march_node::InputPortDef>()->data() = x;
  ctx.run();
  const std::vector<vec3> &out =
      node->get_datum<duchamp::march_node::OutputPortDef>()->data();

  GAUDI_ASSERT(out.size() == x.size());

  // Node mechanics: out = in + step * normal, exactly.
  const std::vector<vec3> &n = normals->get();
  real max_err = 0.0;
  real max_disp_err = 0.0;
  real mean_radius = 0.0;
  for (size_t i = 0; i < x.size(); ++i) {
    const vec3 expect = x[i] + step * n[i];
    max_err = std::max(max_err, (out[i] - expect).norm());
    // vert_normal is unit length, so each point moves exactly `step`.
    max_disp_err = std::max(max_disp_err, std::abs((out[i] - x[i]).norm() - step));
    mean_radius += out[i].norm();
  }
  mean_radius /= real(std::max<size_t>(out.size(), 1));

  GAUDI_ASSERT(max_err < 1e-10);
  GAUDI_ASSERT(max_disp_err < 1e-10);
  // make_sphere winds outward, so a +step march grows the radius to ~1 + step.
  GAUDI_ASSERT(std::abs(mean_radius - (1.0 + step)) < 5e-2);
}

} // namespace test
} // namespace gaudi

#endif // __GAUDI_TEST_FIELD_TESTS_HPP__
