#ifndef GAUDI_TEST_SHELL_QC_TESTS_HPP
#define GAUDI_TEST_SHELL_QC_TESTS_HPP

#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/dynamic.hpp"
#include "gaudi/asawa/shell/element_qc.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/asawa/shell/shell_invariants.hpp"
#include "gaudi/test/test.hpp"

#include <iostream>
#include <random>
#include <vector>

namespace gaudi {
namespace test {

namespace shell_qc_detail {

using namespace asawa::shell;

inline bool expect_valid(const shell &M, const char *where, bool triangles = true) {
  const invariant_report r = check_invariants(M, triangles);
  if (!r.ok)
    std::cout << "[shell_qc] " << where << ": " << r.which << " id " << r.id
              << ": " << r.detail << std::endl;
  return r.ok;
}

inline int euler(const shell &M) {
  int V = 0, E = 0, F = 0;
  for (int v = 0; v < static_cast<int>(M.vert_count()); ++v)
    V += M.vbegin(vert_id(v)) >= 0;
  for (int c = 0; c < static_cast<int>(M.corner_count()); c += 2)
    E += M.next(corner_id(c)) >= 0;
  for (int f = 0; f < static_cast<int>(M.face_count()); ++f)
    F += M.fbegin(face_id(f)) >= 0;
  return V - E + F;
}

struct vef {
  int V = 0, E = 0, F = 0;
};

inline vef counts(const shell &M) {
  vef r;
  for (int v = 0; v < static_cast<int>(M.vert_count()); ++v)
    r.V += M.vbegin(vert_id(v)) >= 0;
  for (int c = 0; c < static_cast<int>(M.corner_count()); c += 2)
    r.E += M.next(corner_id(c)) >= 0;
  for (int f = 0; f < static_cast<int>(M.face_count()); ++f)
    r.F += M.fbegin(face_id(f)) >= 0;
  return r;
}

inline real ring_min_quality(const shell &M, VertId v, const std::vector<vec3> &x) {
  real lo = 1;
  M.const_for_each_vertex(v, [&](CornerId s, const shell &m) {
    lo = std::min(lo, face_quality(x[m.vert(s)], x[m.vert(m.next(s))], x[m.vert(m.prev(s))]));
  });
  return lo;
}

/// Collapse `c` the way the remesher does: datum hooks, then topology.
inline void collapse_with_hooks(shell &M, CornerId c) {
  for (auto d : M.get_data())
    d->collapse(M, c);
  collapse_edge(M, c);
}

/// Torus with tangential-ish noise so remeshing has uneven edges to work on.
inline shell::ptr noisy_torus(int maj, int min, real noise, unsigned seed = 7) {
  auto M = load_torus(1.0, 0.33, maj, min);
  std::vector<vec3> &x = asawa::get_vec_data(*M, 0);
  std::mt19937 rng(seed);
  std::uniform_real_distribution<real> u(-1, 1);
  for (vec3 &p : x)
    p += noise * vec3(u(rng), u(rng), u(rng));
  return M;
}

} // namespace shell_qc_detail

GAUDI_TEST(shell_invariants_hold_on_loaded_meshes) {
  using namespace shell_qc_detail;
  GAUDI_EXPECT(expect_valid(*load_tet(), "tet"));
  GAUDI_EXPECT(expect_valid(*load_torus(1.0, 0.33, 12, 8), "torus"));
  GAUDI_EXPECT(euler(*load_tet()) == 2);
  GAUDI_EXPECT(euler(*load_torus(1.0, 0.33, 12, 8)) == 0);
}

GAUDI_TEST(shell_invariants_catch_broken_links_and_labels) {
  using namespace shell_qc_detail;
  {
    auto M = load_torus(1.0, 0.33, 12, 8);
    M->corners_prev()[5] = M->corners_prev()[7];
    const invariant_report r = check_invariants(*M);
    GAUDI_EXPECT(!r.ok && r.which == "I1");
  }
  {
    auto M = load_torus(1.0, 0.33, 12, 8);
    M->corners_face()[5] = M->corners_face()[5] == 0 ? 1 : 0;
    const invariant_report r = check_invariants(*M);
    GAUDI_EXPECT(!r.ok && r.which == "I4");
  }
  {
    auto M = load_torus(1.0, 0.33, 12, 8);
    M->set_vbegin(vert_id(3), corner_id(-1));
    const invariant_report r = check_invariants(*M);
    GAUDI_EXPECT(!r.ok && r.which == "I6");
  }
}

/// The remesher's existing operations on a noisy torus: invariants and the
/// Euler characteristic after each sub-step and after a full `step`.
GAUDI_TEST(shell_invariants_hold_through_dynamic_step) {
  using namespace shell_qc_detail;
  auto M = noisy_torus(24, 12, 0.02);
  const real l0 = avg_length(*M, asawa::const_get_vec_data(*M, 0));
  auto dyn = dynamic::create(M, 0.5 * l0, 1.25 * l0, 0.25 * l0);
  dyn->update_trees();
  dyn->subdivide_edges();
  GAUDI_EXPECT(expect_valid(*M, "after subdivide"));
  GAUDI_EXPECT(euler(*M) == 0);
  dyn->collapse_edges();
  GAUDI_EXPECT(expect_valid(*M, "after collapse"));
  GAUDI_EXPECT(euler(*M) == 0);
  dyn->flip_edges();
  GAUDI_EXPECT(expect_valid(*M, "after flip"));
  GAUDI_EXPECT(euler(*M) == 0);
  for (int k = 0; k < 3; ++k) {
    dyn->update_trees();
    dyn->step(true);
    GAUDI_EXPECT(expect_valid(*M, "after step"));
    GAUDI_EXPECT(euler(*M) == 0);
  }
}

GAUDI_TEST(shell_collapse_lands_on_the_midpoint) {
  using namespace shell_qc_detail;
  auto M = load_torus(1.0, 0.33, 12, 8);
  std::vector<vec3> &x = asawa::get_vec_data(*M, 0);
  const CornerId c = corner_id(10);
  GAUDI_ASSERT(can_collapse(*M, c));
  const VertId v0 = M->vert(c), v1 = M->vert(M->other(c));
  const vec3 mid = real(0.5) * (x[v0] + x[v1]);
  for (auto d : M->get_data())
    d->collapse(*M, c);
  collapse_edge(*M, c);
  GAUDI_EXPECT(M->vbegin(v0) < 0);
  GAUDI_EXPECT(M->vbegin(v1) >= 0);
  GAUDI_EXPECT((x[v1] - mid).norm() < 1e-12);
  GAUDI_EXPECT(expect_valid(*M, "after collapse"));
  GAUDI_EXPECT(euler(*M) == 0);
}

/// On a 3 x 3 torus every vertex has valence 6 and most vertex pairs are
/// already joined: many collapses would pinch, and after one legal flip
/// (edge 8) six flips would duplicate an edge. A 12 x 8 torus has neither.
GAUDI_TEST(shell_flip_and_collapse_refuse_duplicates_and_pinches) {
  using namespace shell_qc_detail;
  auto small = load_torus(1.0, 0.33, 3, 3);
  GAUDI_ASSERT(expect_valid(*small, "3x3 torus"));
  auto dyn = dynamic::create(small, 0.1, 10.0, 0.1);
  GAUDI_ASSERT(!dyn->skip_flip(*small, corner_id(8)));
  flip_edge(*small, corner_id(8));
  GAUDI_ASSERT(expect_valid(*small, "3x3 torus after one flip"));
  int dup = 0, pinch = 0;
  for (CornerId c : small->get_edge_range()) {
    const VertId a = small->vert(small->prev(c));
    const VertId b = small->vert(small->prev(small->other(c)));
    if (small->find_edge_from_verts(a, b) >= 0) {
      ++dup;
      GAUDI_EXPECT(dyn->skip_flip(*small, c));
    }
    if (!collapse_link_ok(*small, c)) {
      ++pinch;
      GAUDI_EXPECT(collapse_edge(*small, c) == c);
    }
  }
  GAUDI_EXPECT(dup > 0);
  GAUDI_EXPECT(pinch > 0);
  GAUDI_EXPECT(expect_valid(*small, "3x3 torus after refusals"));
  GAUDI_EXPECT(euler(*small) == 0);

  auto big = load_torus(1.0, 0.33, 12, 8);
  auto dyn_big = dynamic::create(big, 0.1, 10.0, 0.1);
  for (CornerId c : big->get_edge_range()) {
    GAUDI_EXPECT(collapse_link_ok(*big, c));
    GAUDI_EXPECT(!dyn_big->skip_flip(*big, c));
  }
}

GAUDI_TEST(shell_best_triangulation_beats_the_fan) {
  using namespace shell_qc_detail;
  std::mt19937 rng(3);
  std::uniform_real_distribution<real> u(-0.15, 0.15);
  const vec3 n(0, 0, 1);
  for (int m = 4; m <= 12; ++m) {
    for (int jitter = 0; jitter < 2; ++jitter) {
      std::vector<vec3> p;
      for (int i = 0; i < m; ++i) {
        const real t = 2 * M_PI * i / m, r = 1 + (jitter ? u(rng) : 0);
        p.push_back(vec3(r * std::cos(t), r * std::sin(t), jitter ? u(rng) : 0));
      }
      real q_best = 0;
      const auto tri = best_triangulation(p, n, face_quality, [](int, int) { return true; }, &q_best);
      GAUDI_ASSERT(static_cast<int>(tri.size()) == m - 2);
      real q_fan = 1;
      for (int k = 1; k + 1 < m; ++k)
        q_fan = std::min(q_fan, face_quality(p[0], p[k], p[k + 1]));
      GAUDI_EXPECT(q_best >= q_fan - 1e-12);
      for (const auto &t : tri)
        GAUDI_EXPECT((p[t[1]] - p[t[0]]).cross(p[t[2]] - p[t[0]]).dot(n) > 0);
      GAUDI_EXPECT(best_triangulation(p, n, face_quality, [](int, int) { return false; }).empty());
    }
  }
}

/// A valence-6 torus vertex and a valence-3 one (a subdivided edge's vertex
/// after one spoke flip) both dissolve to a valid triangle mesh with the
/// Euler table's counts and no new slots.
GAUDI_TEST(shell_dissolve_vertex_counts_and_validity) {
  using namespace shell_qc_detail;
  for (int valence : {6, 3}) {
    auto M = load_torus(1.0, 0.33, 12, 8);
    std::vector<vec3> &x = asawa::get_vec_data(*M, 0);
    VertId v = vert_id(20);
    if (valence == 3) {
      const CornerId e = corner_id(40);
      const vec3 mid = real(0.5) * (x[M->vert(e)] + x[M->vert(M->other(e))]);
      // subdivide_edge returns a corner leaving the new vertex.
      v = M->vert(subdivide_edge(*M, e));
      x.resize(M->vert_count());
      x[v] = mid;
      auto dyn = dynamic::create(M, 0.1, 10.0, 0.1);
      GAUDI_ASSERT(M->vsize(v) == 4);
      const CornerId spoke = M->vbegin(v);
      GAUDI_ASSERT(!dyn->skip_flip(*M, spoke));
      flip_edge(*M, spoke);
      GAUDI_ASSERT(M->vsize(v) == 3);
      GAUDI_ASSERT(expect_valid(*M, "valence-3 setup"));
    }
    const vef before = counts(*M);
    const size_t nc = M->corner_count(), nf = M->face_count();
    GAUDI_ASSERT(dissolve_vertex(*M, v, x, face_quality));
    const vef after = counts(*M);
    GAUDI_EXPECT(expect_valid(*M, "after dissolve"));
    GAUDI_EXPECT(after.V == before.V - 1 && after.E == before.E - 3 && after.F == before.F - 2);
    GAUDI_EXPECT(M->corner_count() == nc && M->face_count() == nf);
    GAUDI_EXPECT(M->vbegin(v) < 0);
    GAUDI_EXPECT(euler(*M) == 0);
  }
}

/// Repeated collapses into one vertex raise its valence; dissolving it
/// leaves a valid mesh whose ring is no worse than the fan was.
GAUDI_TEST(shell_dissolve_high_valence_vertex) {
  using namespace shell_qc_detail;
  auto M = load_torus(1.0, 0.33, 24, 12);
  const std::vector<vec3> &x = asawa::const_get_vec_data(*M, 0);
  const VertId hub = vert_id(100);
  for (int it = 0; it < 6 && M->vsize(hub) < 10; ++it) {
    std::vector<CornerId> spokes;
    M->for_each_vertex(hub, [&](CornerId s, shell &) { spokes.push_back(s); });
    bool done = false;
    for (CornerId s : spokes) {
      const CornerId into = M->other(s); // removes the neighbor, keeps hub
      if (can_collapse(*M, into)) {
        collapse_with_hooks(*M, into);
        done = true;
        break;
      }
    }
    GAUDI_ASSERT(done);
  }
  GAUDI_ASSERT(M->vsize(hub) >= 10);
  GAUDI_ASSERT(expect_valid(*M, "after building the hub"));
  const real q_fan = ring_min_quality(*M, hub, x);
  const std::vector<VertId> ring = M->get_one_ring(hub);
  GAUDI_ASSERT(dissolve_vertex(*M, hub, x, face_quality, q_fan));
  GAUDI_EXPECT(expect_valid(*M, "after dissolving the hub"));
  GAUDI_EXPECT(euler(*M) == 0);
  // The faces filling the hole are the ones with all three vertices on the
  // old ring; the dissolve guaranteed them at least the fan's minimum.
  auto on_ring = [&](VertId w) { return std::find(ring.begin(), ring.end(), w) != ring.end(); };
  int hole = 0;
  for (FaceId f : M->get_face_range()) {
    const auto t = M->get_tri(f);
    if (!on_ring(t[0]) || !on_ring(t[1]) || !on_ring(t[2]))
      continue;
    ++hole;
    GAUDI_EXPECT(face_quality(x[t[0]], x[t[1]], x[t[2]]) >= q_fan - 1e-12);
  }
  GAUDI_EXPECT(hole == static_cast<int>(ring.size()) - 2);
}

GAUDI_TEST(shell_predictions_match_the_operations) {
  using namespace shell_qc_detail;
  auto M = noisy_torus(12, 8, 0.03);
  std::vector<vec3> &x = asawa::get_vec_data(*M, 0);
  {
    const CornerId c = corner_id(30);
    GAUDI_ASSERT(can_collapse(*M, c));
    const collapse_prediction pc = predict_collapse(*M, c, x, face_quality);
    const VertId keep = M->vert(M->other(c));
    collapse_with_hooks(*M, c);
    GAUDI_EXPECT(std::abs(ring_min_quality(*M, keep, x) - pc.q_after) < 1e-12);
    GAUDI_EXPECT(static_cast<int>(M->vsize(keep)) == pc.valence);
    GAUDI_EXPECT(!pc.inverts);
  }
  {
    const CornerId c = corner_id(90);
    auto dyn = dynamic::create(M, 0.1, 10.0, 0.1);
    GAUDI_ASSERT(!dyn->skip_flip(*M, c));
    const flip_prediction pf = predict_flip(*M, c, x, face_quality);
    flip_edge(*M, c);
    const FaceId f0 = M->face(c), f1 = M->face(M->other(c));
    auto fq = [&](FaceId f) {
      const CornerId b = M->fbegin(f);
      return face_quality(x[M->vert(b)], x[M->vert(M->next(b))], x[M->vert(M->prev(b))]);
    };
    GAUDI_EXPECT(std::abs(std::min(fq(f0), fq(f1)) - pf.q_after) < 1e-12);
    GAUDI_EXPECT(expect_valid(*M, "after flip"));
  }
}

GAUDI_TEST(shell_retriangulate_face_restores_triangles) {
  using namespace shell_qc_detail;
  auto M = load_torus(1.0, 0.33, 12, 8);
  const std::vector<vec3> &x = asawa::const_get_vec_data(*M, 0);
  const vef before = counts(*M);
  const CornerId c = corner_id(12);
  const FaceId f = merge_face(*M, c, M->other(c));
  GAUDI_ASSERT(face_len(*M, M->fbegin(f)) == 4);
  GAUDI_ASSERT(expect_valid(*M, "quad", false));
  GAUDI_ASSERT(retriangulate_face(*M, f, x, face_quality));
  GAUDI_EXPECT(expect_valid(*M, "after retriangulate"));
  const vef after = counts(*M);
  GAUDI_EXPECT(after.V == before.V && after.E == before.E && after.F == before.F);
}

GAUDI_TEST(shell_refused_collapses_move_nothing) {
  using namespace shell_qc_detail;
  auto M = noisy_torus(12, 8, 0.03);
  const std::vector<vec3> x0 = asawa::const_get_vec_data(*M, 0);
  auto dyn = dynamic::create(M, 10.0, 100.0, 0.1); // every edge is "short"
  dyn->set_collapse_pred([](shell &, CornerId) { return true; });
  dyn->collapse_edges();
  const std::vector<vec3> &x = asawa::const_get_vec_data(*M, 0);
  real moved = 0;
  for (size_t i = 0; i < x0.size(); ++i)
    moved = std::max(moved, (x[i] - x0[i]).norm());
  GAUDI_EXPECT(moved == 0);
  GAUDI_EXPECT(expect_valid(*M, "after refused collapses"));
}

namespace shell_qc_detail {

inline real min_face_quality(const shell &M, const std::vector<vec3> &x) {
  real lo = 1;
  for (FaceId f : M.get_face_range()) {
    const auto t = M.get_tri(f);
    lo = std::min(lo, face_quality(x[t[0]], x[t[1]], x[t[2]]));
  }
  return lo;
}

inline qc_report run_qc(dynamic &dyn) {
  qc_report rep;
  dyn.qc_pass(rep);
  dynamic::log_qc(rep);
  return rep;
}

} // namespace shell_qc_detail

/// A needle (an edge shrunk to 1% of its length) and a cap (a vertex pushed
/// onto its opposite edge) are both repaired, and nothing else breaks.
GAUDI_TEST(shell_qc_repairs_needle_and_cap) {
  using namespace shell_qc_detail;
  for (int kind = 0; kind < 2; ++kind) {
    auto M = load_torus(1.0, 0.33, 24, 12);
    std::vector<vec3> &x = asawa::get_vec_data(*M, 0);
    const CornerId c = corner_id(200);
    const VertId v0 = M->vert(c), v1 = M->vert(M->other(c));
    if (kind == 0) {
      x[v0] = x[v1] + real(0.01) * (x[v0] - x[v1]);
    } else {
      const VertId a = M->vert(M->next(M->next(c))), b = M->vert(M->next(c));
      x[v0] = real(0.5) * (x[a] + x[b]) + real(0.002) * (x[v0] - real(0.5) * (x[a] + x[b]));
    }
    GAUDI_ASSERT(min_face_quality(*M, x) < real(0.05));
    auto dyn = dynamic::create(M, 0.01, 10.0, 0.01);
    const qc_report rep = run_qc(*dyn);
    GAUDI_EXPECT(expect_valid(*M, kind == 0 ? "needle repaired" : "cap repaired"));
    GAUDI_EXPECT(euler(*M) == 0);
    GAUDI_EXPECT(rep.final.faces == 0);
    GAUDI_EXPECT(min_face_quality(*M, x) >= dyn->qc().q_warn);
    GAUDI_EXPECT(rep.collapsed + rep.flipped + rep.dissolved > 0);
  }
}

GAUDI_TEST(shell_qc_brings_high_valence_into_band) {
  using namespace shell_qc_detail;
  auto M = load_torus(1.0, 0.33, 24, 12);
  const VertId hub = vert_id(100);
  for (int it = 0; it < 6 && M->vsize(hub) < 11; ++it) {
    std::vector<CornerId> spokes;
    M->for_each_vertex(hub, [&](CornerId s, shell &) { spokes.push_back(s); });
    for (CornerId s : spokes)
      if (can_collapse(*M, M->other(s))) {
        collapse_with_hooks(*M, M->other(s));
        break;
      }
  }
  GAUDI_ASSERT(M->vsize(hub) >= 10);
  auto dyn = dynamic::create(M, 0.01, 10.0, 0.01);
  const qc_report rep = run_qc(*dyn);
  GAUDI_EXPECT(expect_valid(*M, "valence repaired"));
  GAUDI_EXPECT(euler(*M) == 0);
  GAUDI_EXPECT(rep.final.val_hi == 0);
  GAUDI_EXPECT(rep.final.hard == 0);
}

/// Full remesh steps with QC on: valid, triangular, no hard faces, and the
/// report agrees with an independent census.
GAUDI_TEST(shell_qc_step_invariants) {
  using namespace shell_qc_detail;
  auto M = noisy_torus(24, 12, 0.03);
  const real l0 = avg_length(*M, asawa::const_get_vec_data(*M, 0));
  auto dyn = dynamic::create(M, 0.5 * l0, 1.25 * l0, 0.25 * l0);
  for (int k = 0; k < 4; ++k) {
    std::vector<vec3> &x = asawa::get_vec_data(*M, 0);
    std::mt19937 rng(11 + k);
    std::uniform_real_distribution<real> u(-1, 1);
    for (vec3 &p : x)
      p += real(0.15) * l0 * vec3(u(rng), u(rng), u(rng));
    dyn->update_trees();
    dyn->step(true);
    const qc_report &rep = dyn->last_qc();
    GAUDI_EXPECT(expect_valid(*M, "after qc step"));
    GAUDI_EXPECT(euler(*M) == 0);
    GAUDI_EXPECT(rep.final.hard == 0);
    GAUDI_EXPECT(rep.final.ngons == 0);
    const qc_census now = dyn->qc_census_now();
    GAUDI_EXPECT(now.faces == rep.final.faces && now.val_hi == rep.final.val_hi);
  }
}

GAUDI_TEST(shell_corner_angle_covers_obtuse) {
  using namespace shell_qc_detail;
  auto M = load_tet();
  std::vector<vec3> &x = asawa::get_vec_data(*M, 0);
  const CornerId c = M->fbegin(face_id(0));
  for (real deg : {real(60), real(90), real(120)}) {
    const real t = deg * real(M_PI) / real(180);
    x[M->vert(c)] = vec3::Zero();
    x[M->vert(M->next(c))] = vec3(1, 0, 0);
    x[M->vert(M->prev(c))] = vec3(std::cos(t), std::sin(t), 0);
    GAUDI_EXPECT(std::abs(angle(*M, c, x) - t) < 1e-6);
  }
}

} // namespace test
} // namespace gaudi

#endif
