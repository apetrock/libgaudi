#ifndef GAUDI_CR_DIRICHLET_TESTS_HPP
#define GAUDI_CR_DIRICHLET_TESTS_HPP

#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/hepworth/block/cr_dirichlet.hpp"
#include "gaudi/hepworth/block/shell_constraints.hpp"
#include "gaudi/hepworth/block/shell_constraints_init.hpp"
#include "gaudi/hepworth/blocks/shell_cr_block.hpp"
#include "gaudi/hepworth/blocks/shell_shape_block.hpp"
#include "gaudi/test/test.hpp"

#include <cmath>
#include <memory>
#include <vector>

namespace gaudi {
namespace test {
namespace {

using hepworth::block::fair_ring;
using hepworth::block::fair_spoke;
using hepworth::block::fit_ring_frame;
using hepworth::block::frame_quat_from_axes;
using hepworth::block::tangent_basis;

quat normal_frame(const vec3 &Nin) {
  vec3 T, B;
  const vec3 N = Nin.normalized();
  tangent_basis(N, T, B);
  return frame_quat_from_axes(T, B, N, quat::Identity());
}

struct cr_soup {
  std::vector<vec3> x;
  std::vector<hepworth::block::cr_tri> tris;
  std::vector<hepworth::block::cr_edge> edges;

  int edge_of(int a, int b) {
    for (int i = 0; i < static_cast<int>(edges.size()); ++i) {
      const auto &e = edges[static_cast<size_t>(i)];
      if ((e.v0 == a && e.v1 == b) || (e.v0 == b && e.v1 == a))
        return i;
    }
    hepworth::block::cr_edge e;
    e.v0 = a;
    e.v1 = b;
    edges.push_back(e);
    return static_cast<int>(edges.size()) - 1;
  }

  void face(int a, int b, int c) {
    hepworth::block::cr_tri t;
    t.v = {a, b, c};
    t.e = {edge_of(a, b), edge_of(b, c), edge_of(c, a)};
    const int fi = static_cast<int>(tris.size());
    for (int k = 0; k < 3; ++k)
      edges[static_cast<size_t>(t.e[static_cast<size_t>(k)])].faces.push_back(fi);
    tris.push_back(std::move(t));
  }
};

} // namespace

/// Unit graph z = -½(x²+y²). Stored shape operator is −dN, so H = −1.
GAUDI_TEST(mean_curvature_is_normal_turning_in_the_frame) {
  const real a = real(0.05);
  const real expect = real(-1) / std::sqrt(real(1) + a * a);

  std::vector<vec3> mass(5, vec3::Ones());
  std::vector<vec3> x(5, vec3::Zero());
  std::vector<vec3> vel(5, vec3::Zero());
  std::vector<vec3> force(5, vec3::Zero());
  x[1] = vec3(a, 0, 0);
  x[2] = vec3(-a, 0, 0);
  x[3] = vec3(0, a, 0);
  x[4] = vec3(0, -a, 0);

  std::vector<quat> frames(5, quat::Identity());
  frames[1] = normal_frame(vec3(a, 0, 1));
  frames[2] = normal_frame(vec3(-a, 0, 1));
  frames[3] = normal_frame(vec3(0, a, 1));
  frames[4] = normal_frame(vec3(0, -a, 1));
  std::vector<quat> omega(5, quat::Identity());
  std::vector<vec4> J(5, vec4::Ones());

  auto pos = hepworth::vec3_block::create(mass, x, vel, force);
  auto frame = hepworth::quat_block::create(J, frames, omega);
  vecX q;
  pos->map_to_x(q);
  frame->map_to_x(q);

  fair_ring ring;
  ring.i = 0;
  for (int j = 1; j <= 4; ++j) {
    fair_spoke sp;
    sp.j = j;
    sp.K = 1;
    sp.ok = 1;
    ring.spokes.push_back(sp);
  }
  const auto rf = fit_ring_frame(*pos, *frame, q, ring);
  GAUDI_EXPECT(rf.axes);
  GAUDI_EXPECT(rf.W.ok);
  const real H = real(0.5) * (rf.W.k1 + rf.W.k2);
  GAUDI_EXPECT(std::abs(H - expect) < real(1e-6));
}

/// Edge rows only. The frame is not an unknown.
GAUDI_TEST(cr_dirichlet_rows_sit_with_the_vertex_blocks) {
  asawa::shell::shell::ptr M = asawa::shell::load_torus(1.0, 0.33, 8, 6);
  std::vector<vec3> &x = asawa::get_vec_data(*M, 0);
  const int nv = static_cast<int>(x.size());
  std::vector<vec3> mass(static_cast<size_t>(nv), vec3::Ones());
  std::vector<vec3> vel(static_cast<size_t>(nv), vec3::Zero());
  std::vector<vec3> force(static_cast<size_t>(nv), vec3::Zero());
  auto pos = hepworth::vec3_block::create(mass, x, vel, force);
  auto cr = std::make_shared<hepworth::block::shell_cr_block>(M);
  hepworth::block::solver_context ctx;
  cr->prepare(ctx);

  std::vector<hepworth::sim_block::ptr> emitted{pos};
  cr->emit_dof_blocks(emitted);
  GAUDI_ASSERT(cr->dof());
  std::vector<hepworth::sim_block::ptr> blocks{pos, cr->dof()};

  std::vector<hepworth::projection_constraint::ptr> constraints;
  hepworth::block::init_cr_dirichlet(*M, constraints, x, real(1), real(100), blocks);
  GAUDI_ASSERT(cr->count() > 0);
  GAUDI_ASSERT(static_cast<int>(constraints.size()) == cr->count());

  vecX q;
  for (const auto &b : blocks)
    b->map_to_x(q);
  index_t id0 = 0;
  std::vector<trip> trips;
  for (const auto &c : constraints)
    c->fill_A(id0, trips);
  GAUDI_EXPECT(id0 == 2 * cr->count());

  vecX p = vecX::Zero(id0);
  for (const auto &c : constraints)
    c->project(q, p);
  GAUDI_EXPECT(p.allFinite());
  GAUDI_EXPECT(p.norm() > real(0));

  hepworth::block::cr_complex cx;
  GAUDI_ASSERT(hepworth::block::build_cr_complex(*M, cx));
  auto at = [&](int i) { return x[static_cast<size_t>(i)]; };
  int matched = 0;
  for (const auto &c : constraints) {
    auto row = std::dynamic_pointer_cast<hepworth::block::cr_dirichlet>(c);
    GAUDI_ASSERT(row);
    const int e = row->dof();
    if (e < 0 || e >= static_cast<int>(cx.edges.size()))
      continue;
    const auto &ed = cx.edges[static_cast<size_t>(e)];
    real gp = 0;
    real gq = 0;
    hepworth::block::unoriented_edge_guidance(
        at, ed.v0, ed.v1, hepworth::block::dirichlet_faces_of(cx, e), gp, gq);
    const index_t r = row->_id0;
    const real mp = e < static_cast<int>(cr->mass.size()) ? cr->mass[static_cast<size_t>(e)]
                                                           : real(0);
    const real mq = e + cr->count() < static_cast<int>(cr->mass.size())
                        ? cr->mass[static_cast<size_t>(e + cr->count())]
                        : real(0);
    GAUDI_EXPECT(std::abs(p(r) - real(100) * mp * gp) < real(1e-8));
    GAUDI_EXPECT(std::abs(p(r + 1) - real(100) * mq * gq) < real(1e-8));
    ++matched;
  }
  GAUDI_EXPECT(matched == cr->count());
}

GAUDI_TEST(cr_shape_is_zero_on_a_flat_patch) {
  cr_soup mesh;
  mesh.x = {vec3(0, 0, 0), vec3(1, 0, 0), vec3(0, 1, 0), vec3(1, 1, 0)};
  mesh.face(0, 1, 2);
  mesh.face(1, 3, 2);
  const auto faces =
      hepworth::block::cr_face_operators(mesh.x, mesh.tris, mesh.edges);
  GAUDI_ASSERT(faces.size() == 2);
  for (const auto &f : faces) {
    GAUDI_EXPECT(f.S.allFinite());
    GAUDI_EXPECT(f.S.norm() < real(1e-8));
  }
}

GAUDI_TEST(cr_shape_trace_matches_the_graph_and_ignores_tangent_sign) {
  cr_soup mesh;
  mesh.x.push_back(vec3(0, 0, 0));
  const int n = 8;
  const real a = real(0.2);
  for (int i = 0; i < n; ++i) {
    const real th = real(2) * M_PI * real(i) / real(n);
    const real px = a * std::cos(th);
    const real py = a * std::sin(th);
    mesh.x.push_back(vec3(px, py, real(-0.5) * (px * px + py * py)));
  }
  for (int i = 0; i < n; ++i)
    mesh.face(0, 1 + i, 1 + (i + 1) % n);
  std::vector<int> around(static_cast<size_t>(n));
  for (int i = 0; i < n; ++i)
    around[static_cast<size_t>(i)] = i;
  const auto faces =
      hepworth::block::cr_face_operators(mesh.x, mesh.tris, mesh.edges);
  const mat3 S = hepworth::block::cr_vertex_operator(faces, around);
  vec3 N;
  GAUDI_ASSERT(hepworth::block::cr_fan_normal(mesh.x, mesh.tris, around, 0, N));
  vec3 T = vec3(1, 0, 0) - N * N.x();
  T.normalize();
  vec3 B = N.cross(T).normalized();
  const auto W = hepworth::block::shape_w_from_operator(S, T, B);
  const auto Wflip = hepworth::block::shape_w_from_operator(S, -T, -B);
  GAUDI_EXPECT(W.ok);
  const real H = real(0.5) * (W.k1 + W.k2);
  const real Hflip = real(0.5) * (Wflip.k1 + Wflip.k2);
  GAUDI_EXPECT(std::abs(H - Hflip) < real(1e-8));
  GAUDI_EXPECT(H < real(-0.2));
}

GAUDI_TEST(cr_shape_on_a_cylinder_separates_the_axial_curvature) {
  cr_soup mesh;
  const int n = 16;
  const int rings = 3;
  const real h = real(0.5);
  for (int k = 0; k < rings; ++k) {
    for (int i = 0; i < n; ++i) {
      const real th = real(2) * M_PI * real(i) / real(n);
      mesh.x.push_back(vec3(std::cos(th), std::sin(th), h * real(k)));
    }
  }
  auto id = [&](int k, int i) { return k * n + (i % n); };
  for (int k = 0; k < rings - 1; ++k) {
    for (int i = 0; i < n; ++i) {
      mesh.face(id(k, i), id(k, i + 1), id(k + 1, i));
      mesh.face(id(k, i + 1), id(k + 1, i + 1), id(k + 1, i));
    }
  }
  const int mid = n;
  std::vector<int> around;
  for (int f = 0; f < static_cast<int>(mesh.tris.size()); ++f) {
    const auto &t = mesh.tris[static_cast<size_t>(f)];
    if (t.v[0] == mid || t.v[1] == mid || t.v[2] == mid)
      around.push_back(f);
  }
  const auto faces =
      hepworth::block::cr_face_operators(mesh.x, mesh.tris, mesh.edges);
  const mat3 S = hepworth::block::cr_vertex_operator(faces, around);
  vec3 N;
  GAUDI_ASSERT(hepworth::block::cr_fan_normal(mesh.x, mesh.tris, around, mid, N));
  vec3 axial(0, 0, 1);
  axial -= N * axial.dot(N);
  axial.normalize();
  vec3 circ = N.cross(axial).normalized();
  const real k_axial = axial.dot(S * axial);
  const real k_circ = circ.dot(S * circ);
  GAUDI_EXPECT(std::abs(k_axial) < real(0.25));
  GAUDI_EXPECT(std::abs(std::abs(k_circ) - real(1)) < real(0.35));
  GAUDI_EXPECT(std::abs(k_axial) < std::abs(k_circ));
  const real k_flip = (-axial).dot(S * -axial);
  GAUDI_EXPECT(std::abs(k_axial - k_flip) < real(1e-8));
}

GAUDI_TEST(cr_shape_fairing_row_uses_the_operator) {
  asawa::shell::shell::ptr M = asawa::shell::load_torus(1.0, 0.33, 8, 6);
  std::vector<vec3> &x = asawa::get_vec_data(*M, 0);
  const int nv = static_cast<int>(x.size());
  std::vector<vec3> mass(static_cast<size_t>(nv), vec3::Ones());
  std::vector<vec3> vel(static_cast<size_t>(nv), vec3::Zero());
  std::vector<vec3> force(static_cast<size_t>(nv), vec3::Zero());
  auto pos = hepworth::vec3_block::create(mass, x, vel, force);
  auto cr = std::make_shared<hepworth::block::shell_cr_block>(M);
  hepworth::block::solver_context ctx;
  cr->prepare(ctx);
  std::vector<hepworth::sim_block::ptr> emitted{pos};
  cr->emit_dof_blocks(emitted);
  std::vector<hepworth::sim_block::ptr> blocks{pos, cr->dof()};
  std::vector<char> live(static_cast<size_t>(nv), 1);
  std::vector<hepworth::projection_constraint::ptr> constraints;
  hepworth::block::init_frame_fair(*M, constraints, x, live, real(1), blocks, {});
  GAUDI_ASSERT(constraints.size() > 1);
  vecX q;
  for (const auto &b : blocks)
    b->map_to_x(q);
  index_t id0 = 0;
  std::vector<trip> trips;
  for (const auto &c : constraints)
    c->fill_A(id0, trips);
  GAUDI_EXPECT(id0 > 0);
  GAUDI_EXPECT(id0 % 3 == 0);
  vecX p = vecX::Zero(id0);
  for (const auto &c : constraints)
    c->project(q, p);
  GAUDI_EXPECT(p.allFinite());
  GAUDI_EXPECT(p.norm() > real(0));

  hepworth::block::cr_complex cx;
  GAUDI_ASSERT(hepworth::block::build_cr_complex(*M, cx));
  const int nE = static_cast<int>(cx.edges.size());
  std::vector<real> coeff(static_cast<size_t>(2 * nE), real(0));
  const index_t coff = blocks[1]->get_offset_idx(0);
  for (int i = 0; i < 2 * nE; ++i)
    coeff[static_cast<size_t>(i)] = q(coff + i);
  const auto faceS = hepworth::block::cr_face_operators(x, cx.tris, cx.edges);
  int spoke_matched = 0;
  for (const auto &c : constraints) {
    auto row = std::dynamic_pointer_cast<hepworth::block::frame_fair>(c);
    if (!row || row->_ids.size() < 2)
      continue;
    const int v = static_cast<int>(row->_ids[0]);
    std::vector<int> faces;
    std::vector<int> edges;
    hepworth::block::stamp_cr_vertex(*M, cx, v, faces, edges);
    vec3 N;
    if (!hepworth::block::cr_fan_normal(x, cx.tris, faces, v, N))
      continue;
    vec3 T;
    if (!hepworth::block::cr_dirichlet_tangent(x, cx.edges, coeff, edges, N, T, nullptr))
      continue;
    vec3 B = N.cross(T);
    if (!(B.squaredNorm() > real(1e-24)))
      continue;
    B.normalize();
    T = B.cross(N);
    const auto W = hepworth::block::shape_w_from_operator(
        hepworth::block::cr_vertex_operator(faceS, faces), T, B);
    if (!W.ok)
      continue;
    for (int k = 0; k < static_cast<int>(row->_edge_weights.size()); ++k) {
      const vec3 d = x[static_cast<size_t>(row->_ids[static_cast<size_t>(k + 1)])] -
                     x[static_cast<size_t>(v)];
      vec3 hit = d;
      const std::vector<hepworth::block::objective> none;
      if (!hepworth::block::select_monge_spoke(none, d, W, nullptr, T, B, N, hit))
        hit = d;
      const index_t r = row->_id0 + 3 * k;
      const vec3 got(p(r), p(r + 1), p(r + 2));
      const vec3 expect = row->_edge_weights[static_cast<size_t>(k)] * hit;
      GAUDI_EXPECT((got - expect).norm() < real(1e-6) * (real(1) + expect.norm()));
      ++spoke_matched;
    }
    if (spoke_matched > 0)
      break;
  }
  GAUDI_EXPECT(spoke_matched > 0);
}

GAUDI_TEST(diamond_split_keeps_the_halves_and_turns_the_spokes) {
  const vec3 p0(real(-1), real(0), real(0));
  const vec3 p1(real(1), real(0), real(0));
  const vec3 a(real(0), real(1), real(0));
  const vec3 b(real(0), real(-1), real(0));
  const vec3 parent = vec3::UnitX();
  vec3 out[4];
  hepworth::block::split_diamond_tangents(p0, p1, a, b, real(0.5), parent, out);
  for (int i = 0; i < 2; ++i) {
    const vec3 u = out[i].normalized();
    GAUDI_EXPECT(std::abs(u.dot(vec3::UnitX()) - real(1)) < real(1e-5));
    GAUDI_EXPECT(std::abs(u.dot(vec3::UnitZ())) < real(1e-5));
  }
  for (int i = 2; i < 4; ++i) {
    GAUDI_EXPECT(out[i].allFinite());
    GAUDI_EXPECT(out[i].norm() > real(0.5));
    GAUDI_EXPECT(std::abs(std::abs(out[i].normalized().dot(vec3::UnitY())) - real(1)) <
                 real(1e-4));
    GAUDI_EXPECT(std::abs(out[i].normalized().dot(vec3::UnitZ())) < real(1e-5));
  }
}

GAUDI_TEST(edge_tangent_split_fills_the_diamond) {
  asawa::shell::shell::ptr M = asawa::shell::load_torus(1.0, 0.33, 8, 6);
  const auto edges = M->get_edge_range();
  GAUDI_ASSERT(!edges.empty());
  const index_t c = edges[0];
  const index_t slot = c / 2;
  std::vector<vec3> t(static_cast<size_t>(M->edge_count()), vec3::Zero());
  const vec3 parent = vec3(1, 0, 0);
  t[static_cast<size_t>(slot)] = parent;
  auto d = hepworth::block::edge_tangent_datum::create(asawa::prim_type::EDGE, t);
  const index_t nE = static_cast<index_t>(M->edge_count());
  d->subdivide(*M, slot, nE, nE + 1, nE + 2, real(0.5), c);
  for (index_t i : {slot, nE, nE + 1, nE + 2}) {
    const vec3 got = d->data()[static_cast<size_t>(i)];
    GAUDI_EXPECT(got.allFinite());
    GAUDI_EXPECT(got.norm() > real(1e-6));
  }
  GAUDI_EXPECT(d->data()[static_cast<size_t>(slot)].normalized().dot(parent) > real(0.5));
  GAUDI_EXPECT(d->data()[static_cast<size_t>(nE)].normalized().dot(parent) > real(0.5));
}

/// The face frame is the sign-aligned average of `vertex_weingarten` at the
/// three vertices. On a torus that tight axis is the tube.
GAUDI_TEST(face_shape_in_the_stored_tangent_follows_the_torus) {
  const real R = 1;
  const real r = 0.33;
  asawa::shell::shell::ptr M = asawa::shell::load_torus(R, r, 24, 12);
  const std::vector<vec3> &x = asawa::const_get_vec_data(*M, 0);
  hepworth::block::cr_topology topo;
  GAUDI_ASSERT(hepworth::block::assemble_cr_topology(*M, topo));
  const Eigen::VectorXd g = hepworth::block::cr_curvature_guidance(*M, x, topo);
  std::vector<real> coeff(static_cast<size_t>(g.size()), real(0));
  for (int i = 0; i < g.size(); ++i)
    coeff[static_cast<size_t>(i)] = g(i);
  const int nV = M->vert_count();
  std::vector<vec3> H(static_cast<size_t>(nV), vec3::Zero());
  std::vector<vec3> tang(static_cast<size_t>(nV), vec3::Zero());
  std::vector<vec3> axes(static_cast<size_t>(nV), vec3::Zero());
  std::vector<char> wrote(static_cast<size_t>(nV), 0);
  GAUDI_ASSERT(hepworth::block::vertex_weingarten(*M, x, coeff, H, tang, axes, wrote));
  int n = 0;
  real max_tight_deg = 0;
  const auto deg = [](const vec3 &a, const vec3 &b) {
    const real c =
        std::min(real(1), std::abs(a.normalized().dot(b.normalized())));
    return std::acos(c) * real(180) / std::acos(real(-1));
  };
  for (const asawa::shell::FaceId f : M->get_face_range(true)) {
    if (M->fsize(f) != 3)
      continue;
    const std::array<asawa::shell::VertId, 3> tri = M->get_tri(f);
    const int v[3] = {static_cast<int>(tri[0]), static_cast<int>(tri[1]),
                      static_cast<int>(tri[2])};
    const vec3 &x0 = x[static_cast<size_t>(v[0])];
    const vec3 &x1 = x[static_cast<size_t>(v[1])];
    const vec3 &x2 = x[static_cast<size_t>(v[2])];
    vec3 N = (x1 - x0).cross(x2 - x0);
    if (!(N.squaredNorm() > real(1e-16)))
      continue;
    N.normalize();
    vec3 axis = vec3::Zero();
    vec3 anchor = vec3::Zero();
    real k1 = 0;
    real k2 = 0;
    int got = 0;
    for (int i = 0; i < 3; ++i) {
      if (!wrote[static_cast<size_t>(v[i])])
        continue;
      vec3 a = axes[static_cast<size_t>(v[i])];
      a -= N * a.dot(N);
      if (!(a.squaredNorm() > real(1e-16)))
        continue;
      a.normalize();
      if (got == 0)
        anchor = a;
      else if (a.dot(anchor) < real(0))
        a = -a;
      axis += a;
      real ek1 = 0;
      real ek2 = 0;
      real ex = 0;
      real ey = 0;
      hepworth::block::shape_sym2(H[static_cast<size_t>(v[i])].x(),
                                  H[static_cast<size_t>(v[i])].y(),
                                  H[static_cast<size_t>(v[i])].z(), ek1, ek2, ex,
                                  ey);
      k1 += ek1;
      k2 += ek2;
      ++got;
    }
    if (got == 0 || !(axis.squaredNorm() > real(1e-16)))
      continue;
    axis.normalize();
    vec3 other = N.cross(axis).normalized();
    const vec3 tight = std::abs(k1) >= std::abs(k2) ? axis : other;
    const vec3 origin = (x0 + x1 + x2) / real(3);
    const real u = std::atan2(origin.y(), origin.x());
    const real rho = std::hypot(origin.x(), origin.y());
    const real vv = std::atan2(origin.z(), rho - R);
    const vec3 tube(-std::cos(u) * std::sin(vv), -std::sin(u) * std::sin(vv),
                    std::cos(vv));
    max_tight_deg = std::max(max_tight_deg, deg(tight, tube));
    ++n;
  }
  GAUDI_EXPECT(n > 24 * 12);
  GAUDI_EXPECT(max_tight_deg < real(8));
}

} // namespace test
} // namespace gaudi

#endif
