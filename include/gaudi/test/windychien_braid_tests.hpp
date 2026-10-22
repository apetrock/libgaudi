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

GAUDI_TEST(duchamp_braid_planar_slide_crosses_several_rows) {
  windychien::braid b;
  b.name = "slide";
  b.strands = 4;
  b.word = {1, 2, 3};
  b.layer_counts = {3};
  duchamp::braid_planar_params p;
  p.center = false;
  p.close = false;
  p.dx = 1.0;
  p.dy = 1.0;
  p.eps_z = 0.05;
  const asawa::rod::rod::ptr R = duchamp::braid_to_planar_rod(b, p);
  GAUDI_EXPECT(static_cast<int>(R->x().size()) == 4 * 2 + 2 * 3);

  int over = 0;
  real y_prev = -1.0;
  for (const vec3 &q : R->x()) {
    if (std::abs(q[2] - p.eps_z) > 1e-12) {
      continue;
    }
    GAUDI_EXPECT(q[1] > y_prev);
    y_prev = q[1];
    ++over;
  }
  GAUDI_EXPECT(over == 3);
  GAUDI_EXPECT(y_prev > 2.0);

  // One weave in the middle of the upward pass: the middle crossing is under.
  b.word = {1, -2, 3};
  const asawa::rod::rod::ptr woven = duchamp::braid_to_planar_rod(b, p);
  GAUDI_EXPECT(static_cast<int>(woven->x().size()) == 4 * 2 + 2 * 3);
  int woven_over = 0;
  int woven_under = 0;
  for (const vec3 &q : woven->x()) {
    if (std::abs(std::abs(q[2]) - p.eps_z) > 1e-12) {
      continue;
    }
    if (q[2] > 0.0) {
      ++woven_over;
    } else {
      ++woven_under;
    }
  }
  GAUDI_EXPECT(woven_over == 3);
  GAUDI_EXPECT(woven_under == 3);
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

GAUDI_TEST(windychien_make_knot_weave_or_catalog) {
  const windychien::braid weave =
      windychien::make_knot({.name = "plain_weave", .strands = 4, .frames = 3});
  const windychien::braid plain = windychien::make_plain_weave(4, 3);
  GAUDI_EXPECT(weave.strands == plain.strands);
  GAUDI_EXPECT(weave.word == plain.word);
  GAUDI_EXPECT(weave.layer_counts == plain.layer_counts);

  const windychien::braid alias =
      windychien::make_knot({.name = "plane_weave", .strands = 4, .frames = 3});
  GAUDI_EXPECT(alias.word == plain.word);

  const windychien::braid cat = windychien::make_knot({.name = "3_1"});
  const windychien::braid loaded = windychien::load_braid("3_1");
  GAUDI_EXPECT(cat.strands == loaded.strands);
  GAUDI_EXPECT(cat.word == loaded.word);

  const windychien::braid kt = windychien::make_knot({.name = "K11n42"});
  GAUDI_EXPECT(kt.strands == 4);
  GAUDI_EXPECT(kt.word.size() == 13u);
  const windychien::braid conway = windychien::load_braid("K11n34");
  GAUDI_EXPECT(conway.strands == 4);
  GAUDI_EXPECT(conway.word.size() == 11u);

  // spherogram braid of the Knot Atlas DT code, not a Stoimenov serial index.
  const windychien::braid k11a17 = windychien::load_braid("K11a17");
  GAUDI_EXPECT(k11a17.strands == 6);
  GAUDI_EXPECT(windychien::braid_component_count(k11a17) == 1);

  const windychien::braid rnd =
      windychien::make_knot({.name = "random", .strands = 8, .frames = 24, .seed = 3});
  const windychien::braid again = windychien::make_random_braid(8, 24, 3);
  GAUDI_EXPECT(rnd.word == again.word);
  GAUDI_EXPECT(rnd.strands == 8);
  // 8 strands need an odd crossing count for one component, so 24 becomes 25.
  GAUDI_EXPECT(rnd.word.size() == 25u);
  GAUDI_EXPECT(!rnd.layer_counts.empty());
  int layer_sum = 0;
  int mid_spans = 0;
  std::vector<char> span_seen(8, 0);
  for (int c : rnd.layer_counts) {
    GAUDI_EXPECT(c >= 1 && c <= 7);
    span_seen[static_cast<size_t>(c)] = 1;
    layer_sum += c;
    if (c >= 2 && c <= 6) {
      ++mid_spans;
    }
  }
  int distinct = 0;
  for (char s : span_seen) {
    if (s) {
      ++distinct;
    }
  }
  GAUDI_EXPECT(layer_sum == static_cast<int>(rnd.word.size()));
  GAUDI_EXPECT(distinct >= 3);
  GAUDI_EXPECT(mid_spans >= 1);
  GAUDI_EXPECT(windychien::braid_component_count(rnd) == 1);
  std::vector<char> used(7, 0);
  for (int g : rnd.word) {
    const int i = std::abs(g);
    GAUDI_EXPECT(i >= 1 && i <= 7);
    used[static_cast<size_t>(i - 1)] = 1;
  }
  for (char u : used) {
    GAUDI_EXPECT(u == 1);
  }
  const windychien::braid other = windychien::make_random_braid(8, 24, 4);
  GAUDI_EXPECT(other.word != rnd.word);

  int off = 0;
  for (int c : rnd.layer_counts) {
    if (c >= 2) {
      bool pos = false;
      bool neg = false;
      for (int k = 0; k < c; ++k) {
        if (rnd.word[static_cast<size_t>(off + k)] > 0) {
          pos = true;
        } else {
          neg = true;
        }
      }
      GAUDI_EXPECT(pos && neg);
    }
    std::vector<int> inv(static_cast<size_t>(c));
    for (int k = 0; k < c; ++k) {
      inv[static_cast<size_t>(k)] =
          -rnd.word[static_cast<size_t>(off + c - 1 - k)];
    }
    int off2 = 0;
    for (int c2 : rnd.layer_counts) {
      if (c2 == c && off2 != off) {
        bool same = true;
        for (int k = 0; k < c; ++k) {
          if (rnd.word[static_cast<size_t>(off2 + k)] != inv[static_cast<size_t>(k)]) {
            same = false;
            break;
          }
        }
        GAUDI_EXPECT(!same);
      }
      off2 += c2;
    }
    off += c;
  }
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
