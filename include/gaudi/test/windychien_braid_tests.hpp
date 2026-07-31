#ifndef GAUDI_TEST_WINDYCHIEN_BRAID_TESTS_HPP
#define GAUDI_TEST_WINDYCHIEN_BRAID_TESTS_HPP

#include "gaudi/duchamp/braid_planar_rod.hpp"
#include "gaudi/duchamp/braid_sphere_rod.hpp"
#include "gaudi/test/test.hpp"
#include "gaudi/windychien/braid.hpp"
#include "gaudi/windychien/catalog.hpp"

#include <cmath>
#include <vector>

namespace gaudi {
namespace test {

GAUDI_TEST(windychien_catalog_loads_full_table) {
  const auto cat = windychien::load_braid_catalog();
  GAUDI_EXPECT(cat.size() >= 2977);
  const windychien::braid b = windychien::find_braid(cat, "3_1");
  GAUDI_EXPECT(b.strands == 2);
  GAUDI_EXPECT(b.word.size() == 3);
  GAUDI_EXPECT(std::abs(b.word[0]) == 1 && std::abs(b.word[1]) == 1 &&
               std::abs(b.word[2]) == 1);

  const windychien::braid t = windychien::find_braid(cat, "T8_7");
  GAUDI_EXPECT(t.strands == 7);
  GAUDI_EXPECT(t.word.size() == 48);
}

GAUDI_TEST(windychien_rejects_bad_generator) {
  windychien::braid b;
  b.name = "bad";
  b.strands = 2;
  b.word = {2};
  bool threw = false;
  try {
    windychien::validate_braid(b);
  } catch (const std::exception &) {
    threw = true;
  }
  GAUDI_EXPECT(threw);
}

GAUDI_TEST(duchamp_braid_planar_crossing_mids_separated) {
  // Coplanar columns + ±eps beads; beads share XY with opposite z.
  for (const char *name : {"3_1", "4_1", "T8_7", "12_477"}) {
    const windychien::braid b = windychien::load_braid(name);
    duchamp::braid_planar_params p;
    p.center = false;
    p.close = false;
    p.eps_z = 0.05;
    p.dx = 1.0;
    p.dy = 1.0;
    const asawa::rod::rod::ptr R = duchamp::braid_to_planar_rod(b, p);
    GAUDI_EXPECT(R != nullptr);

    const int L = static_cast<int>(b.word.size());
    GAUDI_EXPECT(static_cast<int>(R->x().size()) ==
                 b.strands * (L + 1) + 2 * L);

    // Pair opposite beads at the same XY midpoint.
    int pairs = 0;
    std::vector<char> used(R->x().size(), 0);
    for (size_t i = 0; i < R->x().size(); ++i) {
      if (used[i] || std::abs(std::abs(R->x()[i][2]) - p.eps_z) > 1e-12) {
        continue;
      }
      for (size_t j = i + 1; j < R->x().size(); ++j) {
        if (used[j] || std::abs(std::abs(R->x()[j][2]) - p.eps_z) > 1e-12) {
          continue;
        }
        const vec3 &a = R->x()[i];
        const vec3 &bpt = R->x()[j];
        if (std::hypot(a[0] - bpt[0], a[1] - bpt[1]) > 1e-9) {
          continue;
        }
        GAUDI_EXPECT(a[2] * bpt[2] < 0.0);
        used[i] = used[j] = 1;
        ++pairs;
        break;
      }
    }
    GAUDI_EXPECT(pairs == L);
  }
}

GAUDI_TEST(duchamp_braid_planar_rod_3_1_open) {
  const windychien::braid b = windychien::load_braid("3_1");
  duchamp::braid_planar_params p;
  p.center = false;
  p.close = false;
  p.eps_z = 0.05;
  const asawa::rod::rod::ptr R = duchamp::braid_to_planar_rod(b, p);
  GAUDI_EXPECT(R != nullptr);

  const int L = static_cast<int>(b.word.size());
  GAUDI_EXPECT(static_cast<int>(R->x().size()) == b.strands * (L + 1) + 2 * L);

  int open_starts = 0;
  for (size_t i = 0; i < R->corner_count(); ++i) {
    if (R->prev(asawa::rod::corner_id(static_cast<int>(i))) < 0) {
      ++open_starts;
    }
  }
  GAUDI_EXPECT(open_starts == b.strands);
}

GAUDI_TEST(duchamp_braid_planar_rod_3_1_closed) {
  const windychien::braid b = windychien::load_braid("3_1");
  duchamp::braid_planar_params p;
  p.close = true;
  p.center = false;
  const asawa::rod::rod::ptr R = duchamp::braid_to_planar_rod(b, p);
  GAUDI_EXPECT(R != nullptr);

  const int L = static_cast<int>(b.word.size());
  GAUDI_EXPECT(static_cast<int>(R->x().size()) == b.strands * L + 2 * L);

  int open_ends = 0;
  for (size_t i = 0; i < R->corner_count(); ++i) {
    const auto ci = asawa::rod::corner_id(static_cast<int>(i));
    if (R->prev(ci) < 0 || R->next(ci) < 0) {
      ++open_ends;
    }
  }
  GAUDI_EXPECT(open_ends == 0);
  // Closure of (σ1 σ2)^3 is one knot component (3-cycle braid perm).
  GAUDI_EXPECT(duchamp::detail::collect_strand_paths(*R, true).size() == 1);
}

GAUDI_TEST(duchamp_braid_planar_rod_T8_7_closed_one_component) {
  const windychien::braid b = windychien::load_braid("T8_7");
  duchamp::braid_planar_params p;
  p.close = true;
  p.center = false;
  const asawa::rod::rod::ptr R = duchamp::braid_to_planar_rod(b, p);
  GAUDI_EXPECT(R != nullptr);
  int open_ends = 0;
  for (size_t i = 0; i < R->corner_count(); ++i) {
    const auto ci = asawa::rod::corner_id(static_cast<int>(i));
    if (R->prev(ci) < 0 || R->next(ci) < 0) {
      ++open_ends;
    }
  }
  GAUDI_EXPECT(open_ends == 0);
  // T(8,7) = closure of (σ1…σ6)^8; gcd(8,7)=1 → single knot.
  GAUDI_EXPECT(duchamp::detail::collect_strand_paths(*R, true).size() == 1);
}

GAUDI_TEST(duchamp_braid_closed_wrap_bead_on_seam) {
  // Wrap crossing beads must sit at x = (L - 1/2)*dx, not mid-chart.
  const windychien::braid b = windychien::load_braid("3_1");
  duchamp::braid_planar_params p;
  p.close = true;
  p.center = false;
  p.dx = 1.0;
  p.dy = 1.0;
  p.eps_z = 0.05;
  const asawa::rod::rod::ptr R = duchamp::braid_to_planar_rod(b, p);
  const int L = static_cast<int>(b.word.size());
  const real seam_x = p.dx * (real(L) - 0.5);
  int seam_beads = 0;
  for (const vec3 &q : R->x()) {
    if (std::abs(std::abs(q[2]) - p.eps_z) > 1e-12) {
      continue;
    }
    if (std::abs(q[0] - seam_x) < 1e-12) {
      ++seam_beads;
    }
  }
  GAUDI_EXPECT(seam_beads == 2); // over + under on the wrap generator
}

GAUDI_TEST(duchamp_braid_sphere_rod_3_1) {
  const windychien::braid b = windychien::load_braid("3_1");
  duchamp::braid_sphere_params sp;
  sp.radius = 1.0;
  sp.n_lat = 32;
  const asawa::rod::rod::ptr R = duchamp::braid_to_sphere_rod(b, {}, sp);
  GAUDI_EXPECT(R != nullptr);
  // Closed + chart-subdivided: denser than the raw crossing skeleton.
  GAUDI_EXPECT(static_cast<int>(R->x().size()) >= b.strands * sp.n_lat / 2);

  int open_ends = 0;
  real r_min = 1e9;
  real r_max = -1e9;
  for (size_t i = 0; i < R->corner_count(); ++i) {
    const auto ci = asawa::rod::corner_id(static_cast<int>(i));
    if (R->prev(ci) < 0 || R->next(ci) < 0) {
      ++open_ends;
    }
    const real rn = R->x()[i].norm();
    r_min = std::min(r_min, rn);
    r_max = std::max(r_max, rn);
  }
  GAUDI_EXPECT(open_ends == 0);
  // Chart-subdivide then project: samples stay near r = R + gain·z.
  GAUDI_EXPECT(std::isfinite(r_min) && std::isfinite(r_max));
  GAUDI_EXPECT(r_min > 0.9);
  GAUDI_EXPECT(r_max < 1.1);
}

GAUDI_TEST(duchamp_braid_sphere_rod_T8_7) {
  const asawa::rod::rod::ptr R = duchamp::load_braid_sphere_rod("T8_7");
  GAUDI_EXPECT(R != nullptr);
  GAUDI_EXPECT(R->x().size() > 7 * 48);

  int open_ends = 0;
  for (size_t i = 0; i < R->corner_count(); ++i) {
    const auto ci = asawa::rod::corner_id(static_cast<int>(i));
    if (R->prev(ci) < 0 || R->next(ci) < 0) {
      ++open_ends;
    }
  }
  GAUDI_EXPECT(open_ends == 0);
}

} // namespace test
} // namespace gaudi

#endif
