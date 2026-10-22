#ifndef GAUDI_ROD_DATUM_STORAGE_TESTS_HPP
#define GAUDI_ROD_DATUM_STORAGE_TESTS_HPP

#include "gaudi/asawa/rod/dynamic.hpp"
#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/asawa/rod/rod_datums.hpp"
#include "gaudi/test/test.hpp"

#include <cmath>
#include <vector>

namespace gaudi {
namespace test {

GAUDI_TEST(rod_datum_pool_install_and_sizes) {
  std::vector<vec3> verts = {
      vec3(0, 0, 0),
      vec3(1, 0, 0),
      vec3(2, 0, 0),
      vec3(3, 0, 0),
  };
  auto R = asawa::rod::rod::create(verts, false);
  GAUDI_EXPECT(R->get_data().size() == 9);
  GAUDI_EXPECT(R->x().size() == R->corner_count());
  GAUDI_EXPECT(R->v().size() == R->corner_count());
  GAUDI_EXPECT(R->u().size() == R->corner_count());
  GAUDI_EXPECT(R->o().size() == R->corner_count());
  GAUDI_EXPECT(R->M().size() == R->corner_count());
  GAUDI_EXPECT(R->J().size() == R->corner_count());
  GAUDI_EXPECT(R->l0().size() == R->corner_count());
  GAUDI_EXPECT(R->t().size() == R->corner_count());
  GAUDI_EXPECT(R->rest_l().size() == R->corner_count());
  GAUDI_EXPECT((R->x()[1] - vec3(1, 0, 0)).norm() < 1e-12);
  GAUDI_EXPECT(R->get_datum(0)->slot_id == 0);
  GAUDI_EXPECT(R->channels().get("x") != nullptr);
  GAUDI_EXPECT(R->channels().velocity_slot() == asawa::rod::datum_id::v);
  asawa::velocity(R->channels())[0] = vec3(1, 2, 3);
  GAUDI_EXPECT((R->v()[0] - vec3(1, 2, 3)).norm() < 1e-12);
}

GAUDI_TEST(rod_datum_split_midpoint_and_interp_l) {
  std::vector<vec3> verts = {
      vec3(0, 0, 0),
      vec3(2, 0, 0),
      vec3(4, 0, 0),
  };
  auto R = asawa::rod::rod::create(verts, false);
  auto Rd = asawa::rod::dynamic::create(R, 0.1, 10.0, 0.0);

  const real l0_before = R->l0()[0];
  GAUDI_EXPECT(std::abs(l0_before - 2.0) < 1e-10);

  // Split first edge [0]->[1]
  Rd->split_edge(asawa::rod::corner_id(0));

  GAUDI_EXPECT(R->corner_count() == 4);
  // New corner is last: Catmull/midpoint of edge 0-1 → near (1,0,0)
  const vec3 &xn = R->x().back();
  GAUDI_EXPECT((xn - vec3(1, 0, 0)).norm() < 1e-8);

  // Rest lengths on the two half-edges sum to original
  const real sum_l = R->l0()[0] + R->l0().back();
  GAUDI_EXPECT(std::abs(sum_l - l0_before) < 1e-8);
}

GAUDI_TEST(rod_datum_split_quat_slerp) {
  std::vector<vec3> verts = {
      vec3(0, 0, 0),
      vec3(1, 0, 0),
      vec3(2, 0, 0),
  };
  auto R = asawa::rod::rod::create(verts, false);
  R->u()[0] = quat::Identity();
  R->u()[1] = quat(Eigen::AngleAxisd(M_PI, vec3::UnitZ()));
  R->u()[2] = quat::Identity();

  auto Rd = asawa::rod::dynamic::create(R, 0.1, 10.0, 0.0);
  Rd->split_edge(asawa::rod::corner_id(0));

  const quat &un = R->u().back();
  const quat expect =
      va::slerp(quat::Identity(), quat(Eigen::AngleAxisd(M_PI, vec3::UnitZ())),
                0.5f);
  GAUDI_EXPECT(un.angularDistance(expect) < 1e-6);
}

GAUDI_TEST(rod_datum_collapse_and_pack) {
  std::vector<vec3> verts = {
      vec3(0, 0, 0),
      vec3(1, 0, 0),
      vec3(2, 0, 0),
      vec3(3, 0, 0),
  };
  auto R = asawa::rod::rod::create(verts, false);
  auto Rd = asawa::rod::dynamic::create(R, 0.1, 10.0, 0.0);

  const size_t n0 = R->corner_count();
  Rd->collapse_edge(asawa::rod::corner_id(1));
  R->pack();

  GAUDI_EXPECT(R->corner_count() == n0 - 1);
  GAUDI_EXPECT(R->x().size() == R->corner_count());
  GAUDI_EXPECT(R->l0().size() == R->corner_count());
  GAUDI_EXPECT(R->rest_l().size() == R->corner_count());
  GAUDI_EXPECT(R->u().size() == R->corner_count());
  for (size_t i = 0; i < R->corner_count(); ++i) {
    if (R->next(asawa::rod::corner_id(static_cast<int>(i))) < 0 &&
        R->prev(asawa::rod::corner_id(static_cast<int>(i))) < 0)
      GAUDI_ASSERT(false); // pack should drop tombstones
  }
}

} // namespace test
} // namespace gaudi

#endif
