#ifndef GAUDI_WEINGARTEN_TESTS_HPP
#define GAUDI_WEINGARTEN_TESTS_HPP

#include "gaudi/asawa/datums.hpp"
#include "gaudi/asawa/primitive_objects.hpp"
#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/asawa/shell/stencil_walker.hpp"
#include "gaudi/hepworth/block/shape_det_conserve.hpp"
#include "gaudi/hepworth/block/variational_audit.hpp"
#include "gaudi/hepworth/block/edge_shape_fair.hpp"
#include "gaudi/hepworth/block/edge_shape_goals.hpp"
#include "gaudi/hepworth/block/generic_constraints.hpp"
#include "gaudi/hepworth/block/shell_constraints_init.hpp"
#include "gaudi/hepworth/block/solver.hpp"
#include "gaudi/hepworth/block/weingarten_gn.hpp"
#include "gaudi/hepworth/blocks/shell_cr_block.hpp"
#include "gaudi/hepworth/blocks/shell_shape_block.hpp"
#include "gaudi/kusama/stencil_table.hpp"
#include "gaudi/kusama/vector_dirichlet.hpp"
#include "gaudi/kusama/weingarten.hpp"
#include "gaudi/test/test.hpp"

#include "cr_curvature_edge_generated.hpp"
#include "cr_curvature_stencil_generated.hpp"

#include <cmath>
#include <map>
#include <random>
#include <set>
#include <vector>

namespace gaudi {
namespace test {
namespace {

inline std::vector<std::array<int, 3>> shell_tris(const asawa::shell::shell &M) {
  std::vector<std::array<int, 3>> tris;
  for (asawa::shell::FaceId f : M.get_face_range(true)) {
    const auto t = M.get_tri(f);
    tris.push_back({static_cast<int>(t[0]), static_cast<int>(t[1]),
                    static_cast<int>(t[2])});
  }
  return tris;
}

/// Open strip as a triangle list (boundary edges on all four sides).
inline void strip_tris(int nx, int ny, std::vector<vec3> &x,
                       std::vector<std::array<int, 3>> &tris) {
  x.clear();
  tris.clear();
  for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j)
      x.push_back(vec3(real(i) / (nx - 1), real(j) / (ny - 1),
                       real(0.1) * std::sin(real(i + 2 * j))));
  for (int i = 0; i + 1 < nx; ++i)
    for (int j = 0; j + 1 < ny; ++j) {
      const int a = i * ny + j, b = (i + 1) * ny + j, c = (i + 1) * ny + j + 1,
                d = i * ny + j + 1;
      tris.push_back({a, b, c});
      tris.push_back({a, c, d});
    }
}

GAUDI_TEST(weingarten_shell_walker_matches_triangle_walker) {
  asawa::shell::shell::ptr M = asawa::shell::load_torus(1.0, 0.33, 12, 8);
  const asawa::shell::shell_walker sw(*M);
  // Same triangles, rotated to start where the walker starts (`fend`), so
  // local indices line up. `get_tri` starts at `fbegin`, one corner later.
  std::vector<std::array<int, 3>> tris = shell_tris(*M);
  GAUDI_ASSERT(static_cast<int>(tris.size()) == sw.face_count());
  for (int f = 0; f < sw.face_count(); ++f) {
    const auto &t = tris[static_cast<size_t>(f)];
    GAUDI_EXPECT(t[0] == sw.vertex(f, 1) && t[1] == sw.vertex(f, 2) &&
                 t[2] == sw.vertex(f, 0));
    tris[static_cast<size_t>(f)] = {sw.vertex(f, 0), sw.vertex(f, 1), sw.vertex(f, 2)};
  }
  const asawa::shell::triangle_walker tw(static_cast<int>(M->vert_count()), tris);
  GAUDI_ASSERT(sw.face_count() == tw.face_count());
  GAUDI_EXPECT(sw.edge_count() == tw.edge_count());
  GAUDI_EXPECT(2 * sw.edge_count() == 3 * sw.face_count()); // closed
  int bad_v = 0, bad_o = 0, boundary = 0;
  std::set<std::pair<int, int>> pair_ids; // (shell edge id, list edge id)
  for (int f = 0; f < sw.face_count(); ++f)
    for (int i = 0; i < 3; ++i) {
      bad_v += sw.vertex(f, i) != tw.vertex(f, i);
      bad_o += sw.opposite(f, i) != tw.opposite(f, i);
      boundary += sw.opposite(f, i) < 0;
      pair_ids.insert({sw.edge(f, i), tw.edge(f, i)});
    }
  GAUDI_EXPECT(bad_v == 0);
  GAUDI_EXPECT(bad_o == 0);
  GAUDI_EXPECT(boundary == 0);
  // The two edge numberings are a bijection.
  GAUDI_EXPECT(static_cast<int>(pair_ids.size()) == sw.edge_count());
  // Shell edge ids follow the CR dof numbering.
  std::vector<int> slot_map;
  const int nE = kusama::build_compact_edge_dof_map(*M, slot_map);
  GAUDI_EXPECT(nE == sw.edge_count());
  int bad_e = 0;
  for (int f = 0; f < sw.face_count(); ++f) {
    bool first = true; // a triangle's fan visits it once, then repeats degenerate
    M->const_for_each_face_tri(
        sw.source_face(f),
        [&](asawa::shell::CornerId c0, asawa::shell::CornerId c1,
            asawa::shell::CornerId c2, const asawa::shell::shell &) {
          if (!first)
            return;
          first = false;
          const asawa::shell::CornerId c[3] = {c0, c1, c2};
          for (int k = 0; k < 3; ++k)
            bad_e += slot_map[static_cast<size_t>(static_cast<int>(c[k]) / 2)] !=
                     sw.edge(f, k);
        });
  }
  GAUDI_EXPECT(bad_e == 0);
}

GAUDI_TEST(weingarten_stencil_table_edges_and_boundary) {
  std::vector<vec3> x;
  std::vector<std::array<int, 3>> tris;
  strip_tris(6, 5, x, tris);
  const asawa::shell::triangle_walker tw(static_cast<int>(x.size()), tris);
  const kusama::stencil_table t = kusama::build_stencil_table(tw);
  int interior = 0, boundary = 0, bad = 0;
  for (int e = 0; e < t.n_edges; ++e) {
    const auto &fs = t.edge_faces[static_cast<size_t>(e)];
    const auto &ls = t.edge_local[static_cast<size_t>(e)];
    if (fs[1] < 0) {
      ++boundary;
      bad += !t.boundary(fs[0], ls[0]);
      continue;
    }
    ++interior;
    for (int s = 0; s < 2; ++s)
      bad += t.face_edges[static_cast<size_t>(fs[static_cast<size_t>(s)])]
                         [static_cast<size_t>(ls[static_cast<size_t>(s)])] != e;
  }
  GAUDI_EXPECT(bad == 0);
  GAUDI_EXPECT(boundary == 2 * (5 + 4)); // perimeter of a 5 x 4 cell grid
  GAUDI_EXPECT(interior + boundary == t.n_edges);
  // Reflected dummy on a boundary slot, real neighbor elsewhere.
  vec3 h;
  int checked = 0;
  for (int f = 0; f < t.n_faces; ++f) {
    const std::array<vec3, 6> P = kusama::stencil_positions(t, x, f, h);
    for (int i = 0; i < 3; ++i) {
      const vec3 want = t.boundary(f, i)
                            ? vec3(P[i] + P[(i + 1) % 3] - P[(i + 2) % 3])
                            : x[static_cast<size_t>(t.ids[static_cast<size_t>(f)][3 + i])];
      GAUDI_EXPECT((P[3 + i] - want).norm() < real(1e-14));
      GAUDI_EXPECT(h[i] == (t.boundary(f, i) ? real(0) : real(1)));
      ++checked;
    }
  }
  GAUDI_EXPECT(checked == 3 * t.n_faces);
}

/// Central differences of the generated read-outs over the 18 stencil
/// coordinates; the scaling identity `J p = -s` (degree -1 homogeneity).
GAUDI_TEST(weingarten_generated_readouts_match_fd) {
  asawa::shell::shell::ptr M = asawa::shell::load_torus(1.0, 0.33, 12, 8);
  std::vector<vec3> x = asawa::get_vec_data(*M, 0);
  std::mt19937 rng(3);
  std::normal_distribution<real> nd(0, 0.02);
  for (vec3 &p : x)
    p += vec3(nd(rng), nd(rng), nd(rng));
  const kusama::stencil_table t =
      kusama::build_stencil_table(asawa::shell::shell_walker(*M));
  real worst_s = 0, worst_r = 0, worst_euler = 0, worst_rot = 0;
  const real eps = 1e-6;
  for (int f = 0; f < t.n_faces; f += 7) {
    vec3 h;
    const std::array<vec3, 6> P = kusama::stencil_positions(t, x, f, h);
    vec3 s;
    Eigen::Matrix<real, 3, 18> Js;
    hepworth::block::cr_generated::cr_curvature_stencil(P, h, s, Js);
    Eigen::Matrix<real, 9, 1> R;
    Eigen::Matrix<real, 9, 18> Jr;
    hepworth::block::cr_generated::cr_curvature_edge_readings(P, h, R, Jr);
    for (int d = 0; d < 18; ++d) {
      auto Pp = P, Pm = P;
      Pp[d / 3][d % 3] += eps;
      Pm[d / 3][d % 3] -= eps;
      vec3 sp, sm;
      Eigen::Matrix<real, 3, 18> Jd;
      hepworth::block::cr_generated::cr_curvature_stencil(Pp, h, sp, Jd);
      hepworth::block::cr_generated::cr_curvature_stencil(Pm, h, sm, Jd);
      worst_s = std::max(worst_s, (Js.col(d) - (sp - sm) / (2 * eps)).cwiseAbs().maxCoeff());
      Eigen::Matrix<real, 9, 1> Rp, Rm;
      Eigen::Matrix<real, 9, 18> Jd9;
      hepworth::block::cr_generated::cr_curvature_edge_readings(Pp, h, Rp, Jd9);
      hepworth::block::cr_generated::cr_curvature_edge_readings(Pm, h, Rm, Jd9);
      worst_r = std::max(worst_r, (Jr.col(d) - (Rp - Rm) / (2 * eps)).cwiseAbs().maxCoeff());
    }
    Eigen::Matrix<real, 18, 1> p;
    for (int k = 0; k < 6; ++k)
      p.segment<3>(3 * k) = P[k];
    worst_euler = std::max(worst_euler, (Jr * p + R).cwiseAbs().maxCoeff());
    // The face's reading in edge 0's frame is the face read-out rotated by
    // the angle of edge 0 in the face frame, which is 0: they agree.
    worst_rot = std::max(worst_rot, (R.head<3>() - s).cwiseAbs().maxCoeff());
  }
  console_logger::info << "weingarten generated read-outs: J_s " << worst_s << " J_R "
                       << worst_r << " J p + R " << worst_euler << " R_0 - s "
                       << worst_rot << std::endl;
  GAUDI_EXPECT(worst_s < real(1e-6));
  GAUDI_EXPECT(worst_r < real(1e-6));
  GAUDI_EXPECT(worst_euler < real(1e-10));
  GAUDI_EXPECT(worst_rot < real(1e-12));
}

// ------------------------------------------------------------------- M2

namespace wg = kusama::weingarten;
namespace hb = hepworth::block;

/// Right-triangle grid on the unit square (CCW, +z normals).
inline void unit_grid(int n, std::vector<vec3> &x, std::vector<std::array<int, 3>> &tris) {
  x.clear();
  tris.clear();
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      x.push_back(vec3(real(i) / (n - 1), real(j) / (n - 1), 0));
  for (int i = 0; i + 1 < n; ++i)
    for (int j = 0; j + 1 < n; ++j) {
      const int a = i * n + j, b = (i + 1) * n + j, c = (i + 1) * n + j + 1, d = i * n + j + 1;
      tris.push_back({a, b, c});
      tris.push_back({a, c, d});
    }
}

/// `asawa::make_torus` as a triangle list.
inline void torus_tris(real R, real r, int major, int minor, std::vector<vec3> &x,
                       std::vector<std::array<int, 3>> &tris) {
  std::vector<std::vector<int>> faces;
  asawa::make_torus(x, faces, R, r, major, minor);
  tris.clear();
  for (const auto &f : faces)
    tris.push_back({f[0], f[1], f[2]});
}

/// Max distance after the best rigid alignment of `x` onto `y` (Kabsch).
inline real rigid_error(const std::vector<vec3> &x, const std::vector<vec3> &y) {
  vec3 cx = vec3::Zero(), cy = vec3::Zero();
  for (size_t i = 0; i < x.size(); ++i) {
    cx += x[i];
    cy += y[i];
  }
  cx /= real(x.size());
  cy /= real(y.size());
  mat3 C = mat3::Zero();
  for (size_t i = 0; i < x.size(); ++i)
    C += (y[i] - cy) * (x[i] - cx).transpose();
  Eigen::JacobiSVD<mat3> svd(C, Eigen::ComputeFullU | Eigen::ComputeFullV);
  mat3 U = svd.matrixU();
  if ((U * svd.matrixV().transpose()).determinant() < 0)
    U.col(2) *= -1;
  const mat3 Rm = U * svd.matrixV().transpose();
  real e = 0;
  for (size_t i = 0; i < x.size(); ++i)
    e = std::max(e, (Rm * (x[i] - cx) + cy - y[i]).norm());
  return e;
}

/// The spin-1 version of the spin-2 construction (CR stiffness times the
/// rotation between edge frames) reproduces Stein's vector Dirichlet matrix,
/// which also checks the `face_element` refactor of `energy()`.
GAUDI_TEST(weingarten_spin1_construction_matches_vector_dirichlet) {
  asawa::shell::shell::ptr M = asawa::shell::load_torus(1.0, 0.33, 10, 8);
  std::vector<vec3> &x = asawa::get_vec_data(*M, 0);
  const kusama::stencil_table t = kusama::build_stencil_table(asawa::shell::shell_walker(*M));
  int nE = 0;
  const wg::matS L = kusama::laplace::vector_dirichlet::energy(*M, x, &nE);
  GAUDI_ASSERT(nE == t.n_edges);
  Eigen::MatrixXd C = Eigen::MatrixXd::Zero(2 * nE, 2 * nE);
  for (int f = 0; f < t.n_faces; ++f) {
    const auto &id = t.ids[static_cast<size_t>(f)];
    const vec3 X[3] = {x[id[0]], x[id[1]], x[id[2]]};
    const Eigen::Matrix<real, 3, 3> K = wg::cr_face_stiffness(X[0], X[1], X[2]);
    const vec3 N = (X[1] - X[0]).cross(X[2] - X[0]).normalized();
    const vec3 T = (X[1] - X[0]).normalized();
    const vec3 B = N.cross(T);
    Eigen::Matrix<real, 2, 2> Rm[3];
    for (int a = 0; a < 3; ++a) {
      const int lo = std::min(id[a], id[(a + 1) % 3]), hi = std::max(id[a], id[(a + 1) % 3]);
      const vec3 v = (x[hi] - x[lo]).normalized(); // canonical edge direction
      const real c = v.dot(T), s = v.dot(B);
      Rm[a] << c, -s, s, c;
    }
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b) {
        const Eigen::Matrix<real, 2, 2> blk = K(a, b) * Rm[a].transpose() * Rm[b];
        const int ea = t.face_edges[f][a], eb = t.face_edges[f][b];
        for (int r = 0; r < 2; ++r)
          for (int c = 0; c < 2; ++c)
            C(ea + r * nE, eb + c * nE) += blk(r, c);
      }
  }
  const Eigen::MatrixXd Ld = Eigen::MatrixXd(L);
  const real diff = (C - Ld).cwiseAbs().maxCoeff() / Ld.cwiseAbs().maxCoeff();
  console_logger::info << "weingarten spin-1 construction vs vector_dirichlet " << diff
                       << std::endl;
  GAUDI_EXPECT(diff < real(1e-12));
}

/// Edge-frame components `(t, N x t)` of a planar tensor field at edge midpoints.
template <typename Field>
inline std::vector<vec3> sample_edges_planar(const kusama::stencil_table &t,
                                             const std::vector<vec3> &x, Field field) {
  std::vector<vec3> h(static_cast<size_t>(t.n_edges), vec3::Zero());
  for (int f = 0; f < t.n_faces; ++f)
    for (int i = 0; i < 3; ++i) {
      const int a = t.ids[f][i], b = t.ids[f][(i + 1) % 3];
      const vec3 tv = (x[std::max(a, b)] - x[std::min(a, b)]).normalized();
      const vec3 bv = vec3::UnitZ().cross(tv);
      const mat3 S = field(real(0.5) * (x[a] + x[b]));
      h[t.face_edges[f][i]] = vec3(tv.dot(S * tv), tv.dot(S * bv), bv.dot(S * bv));
    }
  return h;
}

GAUDI_TEST(weingarten_spin2_null_space_and_python_parity) {
  auto energy = [](int n, auto field) {
    std::vector<vec3> x;
    std::vector<std::array<int, 3>> tris;
    unit_grid(n, x, tris);
    const kusama::stencil_table t =
        kusama::build_stencil_table(asawa::shell::triangle_walker(static_cast<int>(x.size()), tris));
    const std::vector<vec3> h = sample_edges_planar(t, x, field);
    Eigen::VectorXd c(3 * t.n_edges);
    for (int e = 0; e < t.n_edges; ++e)
      c.segment<3>(3 * e) = h[e];
    return real(c.dot(wg::build_spin2(t, x) * c));
  };
  auto constant = [](const vec3 &) {
    mat3 S;
    S << 1.0, 0.3, 0, 0.3, -0.5, 0, 0, 0, 0;
    return S;
  };
  auto planar = [](const vec3 &p) {
    const real a = 1 + 0.5 * std::sin(2 * M_PI * p[1]);
    const real z0 = a * std::cos(2 * M_PI * p[0]), z1 = a * std::sin(2 * M_PI * p[0]);
    const real tau = 1 + 0.3 * p[0];
    mat3 S;
    S << tau / 2 + z0, z1, 0, z1, tau / 2 - z0, 0, 0, 0, 0;
    return S;
  };
  const real Ec = 5 * M_PI * M_PI + 0.045;
  const real e0 = energy(7, constant);
  const real r11 = (energy(11, planar) - Ec) / Ec;
  const real r21 = (energy(21, planar) - Ec) / Ec;
  console_logger::info << "weingarten spin-2: constant " << e0 << " rel err n11 " << r11
                       << " n21 " << r21 << std::endl;
  GAUDI_EXPECT(std::abs(e0) < real(1e-12));
  // MATH.md section 14, right grid: -8.19e-3 (n 11), -2.05e-3 (n 21).
  GAUDI_EXPECT(std::abs(r11 - real(-8.19e-3)) < real(5e-5));
  GAUDI_EXPECT(std::abs(r21 - real(-2.05e-3)) < real(5e-5));
}

GAUDI_TEST(weingarten_memory_and_coupling_derivatives) {
  std::vector<vec3> x0, xt;
  std::vector<std::array<int, 3>> tris;
  torus_tris(2.0, 0.7, 16, 8, x0, tris);
  const kusama::stencil_table t =
      kusama::build_stencil_table(asawa::shell::triangle_walker(static_cast<int>(x0.size()), tris));
  const wg::memory_rest m = wg::build_memory_rest(t, x0);
  Eigen::VectorXd r;
  wg::memory_residual(t, m, x0, real(1), r, nullptr);
  GAUDI_EXPECT(r.norm() < real(1e-12)); // zero at rest
  const mat3 Q = Eigen::AngleAxis<real>(0.7, vec3(1, 2, 3).normalized()).toRotationMatrix();
  std::vector<vec3> xr(x0.size());
  for (size_t i = 0; i < x0.size(); ++i)
    xr[i] = Q * x0[i] + vec3(0.3, -1, 2);
  wg::memory_residual(t, m, xr, real(1), r, nullptr);
  GAUDI_EXPECT(r.norm() < real(1e-10)); // rigid-invariant
  // Coupling idle at the readings of its own shape.
  std::vector<real> mass(static_cast<size_t>(t.n_edges));
  for (int e = 0; e < t.n_edges; ++e)
    mass[e] = m.D[e].mass;
  const std::vector<vec3> h0 = wg::edge_readings(t, x0);
  wg::edge_coupling(t, x0, mass, &h0, r, nullptr);
  GAUDI_EXPECT(r.norm() < real(1e-12));
  // Gradients against central differences at a perturbed shape.
  std::mt19937 rng(11);
  std::normal_distribution<real> nd(0, 0.03);
  std::vector<vec3> x = x0, d(x0.size());
  for (size_t i = 0; i < x.size(); ++i) {
    x[i] += vec3(nd(rng), nd(rng), nd(rng));
    d[i] = vec3(nd(rng), nd(rng), nd(rng));
  }
  torus_tris(2.0, 0.9, 16, 8, xt, tris);
  const std::vector<vec3> ht = wg::edge_readings(t, xt);
  auto check = [&](auto res, const char *name) {
    Eigen::VectorXd r0;
    std::vector<wg::trip> jt;
    res(x, r0, &jt);
    wg::matS J(r0.size(), 3 * static_cast<int>(x.size()));
    J.setFromTriplets(jt.begin(), jt.end());
    Eigen::VectorXd dv(3 * x.size());
    for (size_t i = 0; i < x.size(); ++i)
      dv.segment<3>(3 * static_cast<index_t>(i)) = d[i];
    const real an = (J.transpose() * r0).dot(dv);
    const real h = 1e-6;
    std::vector<vec3> xp = x, xm = x;
    for (size_t i = 0; i < x.size(); ++i) {
      xp[i] += h * d[i];
      xm[i] -= h * d[i];
    }
    Eigen::VectorXd rp, rm;
    res(xp, rp, nullptr);
    res(xm, rm, nullptr);
    const real fd = (real(0.5) * rp.squaredNorm() - real(0.5) * rm.squaredNorm()) / (2 * h);
    console_logger::info << "weingarten " << name << " grad " << an << " fd " << fd << std::endl;
    GAUDI_EXPECT(std::abs(an - fd) < real(1e-6) * (1 + std::abs(fd)));
  };
  check([&](const std::vector<vec3> &y, Eigen::VectorXd &rr, std::vector<wg::trip> *jt) {
          wg::memory_residual(t, m, y, real(1), rr, jt);
        },
        "memory");
  check([&](const std::vector<vec3> &y, Eigen::VectorXd &rr, std::vector<wg::trip> *jt) {
          wg::edge_coupling(t, y, mass, &ht, rr, jt);
        },
        "coupling");
}

/// Python parity (MATH.md sections 11, 14): torus tube 0.7 -> 0.9 with the
/// edge-average coupling, the intrinsic vector Dirichlet memory and restarts
/// reaches the target shape.
GAUDI_TEST(weingarten_solver_torus_tube_restarts) {
  std::vector<vec3> x0, xt;
  std::vector<std::array<int, 3>> tris;
  torus_tris(2.0, 0.7, 24, 12, x0, tris);
  torus_tris(2.0, 0.9, 24, 12, xt, tris);
  const kusama::stencil_table t =
      kusama::build_stencil_table(asawa::shell::triangle_walker(static_cast<int>(x0.size()), tris));
  wg::weingarten_solver s;
  s.reset(t, x0);
  const std::vector<vec3> h = s.readout(xt);
  std::vector<vec3> x = x0;
  const real e_start = rigid_error(x0, xt);
  const real E0 = s.energy(x, h);
  s.solve(x, h, 6, 30);
  const real e_end = rigid_error(x, xt);
  console_logger::info << "weingarten torus tube: rigid err " << e_start << " -> " << e_end
                       << "  E0 " << E0 << std::endl;
  GAUDI_EXPECT(e_end < real(1e-5));
}

/// Cube `[-1, 1]^3`, `k x k` quads per side, two outward triangles each
/// (`closed_targets.cube`).
inline void cube_tris(int k, std::vector<vec3> &x, std::vector<std::array<int, 3>> &tris) {
  x.clear();
  tris.clear();
  std::map<std::array<long, 3>, int> index;
  auto vid = [&](const vec3 &p) {
    const std::array<long, 3> key{std::lround(p[0] * 1e6), std::lround(p[1] * 1e6),
                                  std::lround(p[2] * 1e6)};
    auto it = index.find(key);
    if (it != index.end())
      return it->second;
    index[key] = static_cast<int>(x.size());
    x.push_back(p);
    return static_cast<int>(x.size()) - 1;
  };
  for (int axis = 0; axis < 3; ++axis)
    for (real side : {real(-1), real(1)}) {
      const int u = (axis + 1) % 3, v = (axis + 2) % 3;
      for (int i = 0; i < k; ++i)
        for (int j = 0; j < k; ++j) {
          std::array<int, 4> q{};
          const int di[4] = {0, 1, 1, 0}, dj[4] = {0, 0, 1, 1};
          for (int c = 0; c < 4; ++c) {
            vec3 p = vec3::Zero();
            p[axis] = side;
            p[u] = real(-1) + real(2) * (i + di[c]) / k;
            p[v] = real(-1) + real(2) * (j + dj[c]) / k;
            q[c] = vid(p);
          }
          tris.push_back({q[0], q[1], q[2]});
          tris.push_back({q[0], q[2], q[3]});
        }
    }
  for (auto &f : tris) {
    const vec3 n = (x[f[1]] - x[f[0]]).cross(x[f[2]] - x[f[0]]);
    if (n.dot(x[f[0]] + x[f[1]] + x[f[2]]) < 0)
      std::swap(f[1], f[2]);
  }
}

/// Python parity (MATH.md section 14, end to end, clean): cube to the
/// sphere of equal area, 8 frames of 30 iterates, memory weight 1e-3.
GAUDI_TEST(weingarten_solver_cube_to_sphere_restarts) {
  std::vector<vec3> x0;
  std::vector<std::array<int, 3>> tris;
  cube_tris(6, x0, tris);
  const real R = std::sqrt(real(24) / (4 * M_PI));
  std::vector<vec3> xt(x0.size());
  for (size_t i = 0; i < x0.size(); ++i)
    xt[i] = R * x0[i].normalized();
  const kusama::stencil_table t =
      kusama::build_stencil_table(asawa::shell::triangle_walker(static_cast<int>(x0.size()), tris));
  wg::weingarten_solver s;
  s.w_memory = real(1e-3);
  s.reset(t, x0);
  const std::vector<vec3> h = s.readout(xt);
  std::vector<vec3> x = x0;
  s.solve(x, h, 8, 30);
  real dist = 0;
  vec3 c = vec3::Zero();
  for (const vec3 &p : x)
    c += p;
  c /= real(x.size());
  for (const vec3 &p : x)
    dist = std::max(dist, std::abs((p - c).norm() - R));
  console_logger::info << "weingarten cube -> sphere: max |r - R| / size " << dist / (2 * R)
                       << "  rigid err " << rigid_error(x, xt) << std::endl;
  GAUDI_EXPECT(dist / (2 * R) < real(1e-4));
}

// ------------------------------------------------------------------- M4

/// Edge-derived vertex charts (local polar of the edge store) against the
/// face-operator charts on an untouched torus, and both against the analytic
/// principal curvatures (tube `1 / r`, hoop `cos(phi) / (R + r cos(phi))`).
GAUDI_TEST(weingarten_edge_charts_match_vertex_charts_on_torus) {
  const real R = 1.0, r = 0.33;
  asawa::shell::shell::ptr M = asawa::shell::load_torus(R, r, 48, 24);
  auto cr = std::make_shared<hepworth::block::shell_cr_block>(M);
  hepworth::block::solver_context ctx;
  cr->prepare(ctx);
  std::vector<vec3> &x = asawa::get_vec_data(*M, 0);
  const int n = static_cast<int>(x.size());
  std::vector<vec3> H(n, vec3::Zero()), tang(n, vec3::Zero()), axes(n, vec3::Zero());
  std::vector<char> wrote(n, 0);
  std::vector<hepworth::block::vertex_chart> charts;
  hepworth::block::cr_complex cx;
  GAUDI_ASSERT(hepworth::block::build_cr_complex(*M, cx));
  const hepworth::block::cr_edge_gradient G =
      hepworth::block::cr_coordinate_gradients(x, cx.tris, cx.edges);
  GAUDI_ASSERT(hepworth::block::vertex_weingarten(*M, x, cr->u, H, tang, axes, wrote, &charts,
                                                  &G.par, &G.perp));
  auto pos = [&](int i) { return x[static_cast<size_t>(i)]; };
  auto store = [&](int e) { return cr->shape[static_cast<size_t>(e)]; };
  real d_old = 0, d_new = 0, d_pair = 0, norm = 0;
  int count = 0;
  for (int i = 0; i < n; ++i) {
    const auto &c = charts[static_cast<size_t>(i)];
    if (!c.ok)
      continue;
    const vec3 he = wg::vertex_tensor(cr->stars, i, pos, store, c.T, c.B, c.N);
    GAUDI_ASSERT(he.allFinite());
    const vec3 &p = x[static_cast<size_t>(i)];
    const real rho = std::hypot(p[0], p[1]);
    const real cphi = (rho - R) / r;
    const real k_tube = 1 / r, k_hoop = cphi / (R + r * cphi);
    auto invariants = [](const vec3 &h) {
      return vec2(h[0] + h[2], h[0] * h[2] - h[1] * h[1]);
    };
    const vec2 exact(k_tube + k_hoop, k_tube * k_hoop);
    const vec2 io = invariants(c.h), ie = invariants(he);
    // The store's sign convention is the generated read-out's; compare |H|.
    d_old += std::pow(std::abs(io[0]) - exact[0], 2);
    d_new += std::pow(std::abs(ie[0]) - exact[0], 2);
    d_pair += (c.h - (io[0] * ie[0] < 0 ? -he : he)).squaredNorm();
    norm += exact[0] * exact[0];
    ++count;
  }
  const real e_old = std::sqrt(d_old / norm), e_new = std::sqrt(d_new / norm);
  const real e_pair = std::sqrt(d_pair / norm);
  console_logger::info << "weingarten charts on torus: rms |H| err face-operator " << e_old
                       << " edge-store " << e_new << "; edge vs face tensors " << e_pair
                       << " over " << count << " vertices" << std::endl;
  GAUDI_EXPECT(count == n);
  GAUDI_EXPECT(e_new < real(0.02));
  GAUDI_EXPECT(e_new < real(2) * e_old + real(1e-3));
}

/// Projection solves on the edge store alone (one block) with fairing rows.
struct store_fixture {
  kusama::stencil_table t;
  std::vector<vec3> x;
  std::vector<vec3> clean;
  std::vector<vec3> hs, hv, hf, hm;
  hepworth::vec3_block::ptr store;
  std::shared_ptr<const hb::edge_frames> frames;

  explicit store_fixture(real noise = 0) {
    std::vector<std::array<int, 3>> tris;
    torus_tris(2.0, 0.7, 24, 12, x, tris);
    t = kusama::build_stencil_table(asawa::shell::triangle_walker(static_cast<int>(x.size()), tris));
    clean = wg::edge_readings(t, x);
    hs = clean;
    std::mt19937 rng(5);
    std::normal_distribution<real> nd(0, 1);
    real rms = 0;
    for (const vec3 &h : clean)
      rms += h.squaredNorm();
    rms = std::sqrt(rms / (3 * clean.size()));
    for (vec3 &h : hs)
      h += noise * rms * vec3(nd(rng), nd(rng), nd(rng));
    const std::vector<wg::cr_edge_rows> D = wg::build_cr_derivative(t, x);
    hm.resize(hs.size());
    for (size_t e = 0; e < hs.size(); ++e)
      hm[e] = vec3::Constant(D[e].mass);
    hv.assign(hs.size(), vec3::Zero());
    hf = hv;
    store = hepworth::vec3_block::create(hm, hs, hv, hf);
    frames = hb::build_edge_frames(t, x);
  }

  void run(std::vector<hepworth::projection_constraint::ptr> cs, int frames_n = 10) {
    hb::projection_solver solver;
    solver.set_constraints(cs);
    std::vector<hepworth::sim_block::ptr> blocks{store};
    for (int k = 0; k < frames_n; ++k)
      solver.step(blocks, real(0.1), real(1), 20);
  }

  real err(const std::vector<vec3> &ref) const {
    real d = 0, n = 0;
    for (size_t e = 0; e < hs.size(); ++e) {
      d += (hs[e] - ref[e]).squaredNorm();
      n += ref[e].squaredNorm();
    }
    return std::sqrt(d / n);
  }
};

GAUDI_TEST(weingarten_edge_fair_rows_idle_at_rest) {
  store_fixture fx;
  vecX q;
  fx.store->map_to_x(q);
  for (hb::edge_fair_rows rows :
       {hb::edge_fair_rows::pin, hb::edge_fair_rows::dirichlet, hb::edge_fair_rows::dirichlet_trace}) {
    auto c = hb::edge_shape_fair::create(&fx.t, fx.x, rows, real(1),
                                         std::vector<hepworth::sim_block::ptr>{fx.store});
    index_t id0 = 0;
    std::vector<trip> trips;
    c->fill_A(id0, trips);
    wg::matS A(id0, q.size());
    A.setFromTriplets(trips.begin(), trips.end());
    vecX p = vecX::Zero(id0);
    c->project(q, p);
    const real r = (A * q - p).norm() / std::max(p.norm(), real(1e-30));
    console_logger::info << "edge_shape_fair rows " << int(rows) << " identity goal residual " << r
                         << " rows " << id0 << std::endl;
    GAUDI_EXPECT(r < real(1e-12));
  }
  // The spin-2 rows hold a constant tensor field with a zero goal.
  std::vector<vec3> xg;
  std::vector<std::array<int, 3>> tris;
  unit_grid(9, xg, tris);
  const kusama::stencil_table tg =
      kusama::build_stencil_table(asawa::shell::triangle_walker(static_cast<int>(xg.size()), tris));
  auto c = hb::edge_shape_fair::create(&tg, xg, hb::edge_fair_rows::dirichlet, real(1),
                                       std::vector<hepworth::sim_block::ptr>{fx.store});
  std::vector<vec3> hc = sample_edges_planar(tg, xg, [](const vec3 &) {
    mat3 S;
    S << 1.0, 0.3, 0, 0.3, -0.5, 0, 0, 0, 0;
    return S;
  });
  std::vector<vec3> hm(hc.size(), vec3::Ones()), hv(hc.size(), vec3::Zero()), hf = hv;
  auto blk = hepworth::vec3_block::create(hm, hc, hv, hf);
  auto c2 = hb::edge_shape_fair::create(&tg, xg, hb::edge_fair_rows::dirichlet, real(1),
                                        std::vector<hepworth::sim_block::ptr>{blk}, hb::goal_zero());
  vecX q2;
  blk->map_to_x(q2);
  index_t id0 = 0;
  std::vector<trip> trips;
  c2->fill_A(id0, trips);
  wg::matS A(id0, q2.size());
  A.setFromTriplets(trips.begin(), trips.end());
  vecX r2 = vecX::Random(q2.size()) * (q2.norm() / std::sqrt(real(q2.size())));
  const real rel = (A * q2).norm() / (A * r2).norm();
  console_logger::info << "edge_shape_fair spin-2 rows on a constant field, relative to random "
                       << rel << std::endl;
  GAUDI_EXPECT(rel < real(1e-8));
}

GAUDI_TEST(weingarten_edge_goals_move_toward_targets) {
  using rows = hb::edge_fair_rows;
  auto blocks = [](const store_fixture &fx) { return std::vector<hepworth::sim_block::ptr>{fx.store}; };
  {
    // Fixed target: pin rows reach it.
    store_fixture fx;
    std::vector<vec3> target = fx.clean;
    for (vec3 &h : target)
      h *= real(1.3);
    const real e0 = fx.err(target);
    fx.run({hb::edge_shape_fair::create(&fx.t, fx.x, rows::pin, real(100), blocks(fx),
                                        hb::goal_fixed(target))});
    const real e1 = fx.err(target);
    console_logger::info << "goal fixed: err " << e0 << " -> " << e1 << std::endl;
    GAUDI_EXPECT(e1 < real(1e-3) * e0);
  }
  {
    // Zero knobs are idle on either row kind.
    store_fixture fx;
    fx.run({hb::edge_shape_fair::create(&fx.t, fx.x, rows::pin, real(100), blocks(fx),
                                        hb::goal_anisotropy(fx.store, hb::anisotropy_target{}, fx.hs)),
            hb::edge_shape_fair::create(&fx.t, fx.x, rows::h1, real(100), blocks(fx),
                                        hb::goal_anisotropy(fx.store, hb::anisotropy_target{}, fx.hs))},
           3);
    console_logger::info << "goal anisotropy (zero knobs): drift " << fx.err(fx.clean)
                         << std::endl;
    GAUDI_EXPECT(fx.err(fx.clean) < real(1e-10));
  }
  {
    // Ratio 2 on every edge that has axes, saddle and parabolic ones included.
    store_fixture fx;
    hb::anisotropy_target spec;
    spec.ratio = 2;
    std::vector<char> has_axes(fx.hs.size());
    for (size_t e = 0; e < fx.hs.size(); ++e)
      has_axes[e] = hb::edge_q(fx.hs[e]) > real(1e-6);
    fx.run({hb::edge_shape_fair::create(&fx.t, fx.x, rows::pin, real(100), blocks(fx),
                                        hb::goal_anisotropy(fx.store, spec, fx.hs))});
    real worst = 0;
    for (size_t e = 0; e < fx.hs.size(); ++e)
      if (has_axes[e])
        worst = std::max(worst, std::abs(std::log(hb::edge_ratio(fx.hs[e]) / real(2))));
    console_logger::info << "goal anisotropy ratio 2: worst |log(r / 2)| " << worst << std::endl;
    GAUDI_EXPECT(worst < real(1e-2));
  }
  for (int which = 0; which < 4; ++which) {
    // Smoothers on a noisy store reduce the error to the clean field.
    store_fixture fx(real(0.1));
    const real e0 = fx.err(fx.clean);
    hepworth::projection_constraint::ptr c;
    const char *name = "";
    if (which == 0) {
      name = "dirichlet zero (grad)";
      c = hb::edge_shape_fair::create(&fx.t, fx.x, rows::dirichlet, real(1),
                                      blocks(fx), hb::goal_zero());
    } else if (which == 1) {
      name = "dirichlet_trace zero (curvature_mean)";
      c = hb::edge_shape_fair::create(&fx.t, fx.x, rows::dirichlet_trace, real(1),
                                      blocks(fx), hb::goal_zero());
    } else if (which == 2) {
      name = "ratio smooth (curvature_ratio)";
      c = hb::edge_shape_fair::create(&fx.t, fx.x, rows::pin, real(10), blocks(fx),
                                      hb::goal_ratio_smooth(fx.store, fx.frames));
    } else {
      name = "principal step (frame minimizer)";
      c = hb::edge_shape_fair::create(&fx.t, fx.x, rows::pin, real(10), blocks(fx),
                                      hb::goal_principal_mean(fx.store, fx.frames));
    }
    // Ratio smoothing is not a denoiser of the tensor: it is judged by the
    // roughness of the log ratio across neighboring edges.
    auto roughness = [&]() {
      real s = 0;
      for (size_t e = 0; e < fx.hs.size(); ++e) {
        const real Le = hb::edge_q(fx.hs[e]);
        for (int j : fx.frames->nbr[e]) {
          const real Lj = hb::edge_q(fx.hs[static_cast<size_t>(j)]);
          s += (Le - Lj) * (Le - Lj);
        }
      }
      return s;
    };
    const real r0 = roughness();
    fx.run({c}, 3);
    const real e1 = fx.err(fx.clean), r1 = roughness();
    console_logger::info << "goal " << name << ": err to clean " << e0 << " -> " << e1
                         << ", log-ratio roughness " << r0 << " -> " << r1 << std::endl;
    // Ratio smoothing and the minimizer smooth the field rather than denoise
    // it (the clean torus is not their fixed point): judge them by roughness.
    GAUDI_EXPECT(which >= 2 ? r1 < r0 : e1 < e0);
  }
}

// ------------------------------------------------ one scale, edits drive

/// Each knob's goal maps a tensor already at its target to itself, so
/// re-projecting every iterate does not compound the edit.
GAUDI_TEST(weingarten_anisotropy_goals_are_idempotent) {
  const vec3 h0(2.0, 0.3, -0.5);
  // A live tensor away from the rest: the goal must not depend on how far
  // the store has got.
  const vec3 h(1.7, -0.2, 0.1);
  real worst = 0;
  auto check = [&](const hb::anisotropy_target &spec) {
    // One application reaches the target; a second must not move it.
    const vec3 g1 = hb::anisotropy_edge(h, h0, spec);
    const vec3 g2 = hb::anisotropy_edge(g1, h0, spec);
    worst = std::max(worst, (g2 - g1).norm() / g1.norm());
    return g1;
  };
  hb::anisotropy_target mean, ratio, logr, rel, rel_log, mixed;
  mean.mean = 1.2;
  ratio.ratio = 1.5;
  logr.ratio = 0.3;
  logr.log = true;
  rel.ratio = 0.8;
  rel.relative = true;
  rel_log.ratio = -0.2;
  rel_log.log = true;
  rel_log.relative = true;
  mixed.mean = 1.1;
  mixed.ratio = 0.1;
  mixed.log = true;
  mixed.relative = true;
  for (const auto *s : {&mean, &ratio, &logr, &rel, &rel_log, &mixed})
    check(*s);
  console_logger::info << "anisotropy goals: worst second-application change " << worst << std::endl;
  GAUDI_EXPECT(worst < real(1e-12));
  GAUDI_EXPECT((hb::anisotropy_edge(h, h0, hb::anisotropy_target{}) - h).norm() < real(1e-14));
}

/// The shape-space ARAP shear: at zero shear the goal is the live store
/// (an identity ARAP); with a shear it turns the store's lines and leaves
/// the mean curvature and the anisotropy where they are.
GAUDI_TEST(weingarten_curvature_shear_goal) {
  store_fixture fx;
  // Analytic charts of the fixture torus (R 2, r 0.7, about z).
  std::vector<hb::vertex_chart> charts(fx.x.size());
  for (size_t i = 0; i < fx.x.size(); ++i) {
    const vec3 &p = fx.x[i];
    const real rho = std::hypot(p.x(), p.y());
    const vec3 c(p.x() / rho * real(2), p.y() / rho * real(2), 0);
    const vec3 n = (p - c).normalized();
    const vec3 hoop(-p.y() / rho, p.x() / rho, 0);
    charts[i].ok = true;
    charts[i].N = n;
    charts[i].T = n.cross(hoop).normalized();
    charts[i].B = n.cross(charts[i].T);
    charts[i].h = vec3(real(1) / real(0.7), 0, ((rho - real(2)) / real(0.7)) / rho);
  }
  vecX q;
  fx.store->map_to_x(q);
  auto goal_of = [&](real shear) {
    hb::curvature_shear_target spec;
    spec.shear = shear;
    std::vector<vec3> g(fx.hs.size());
    hb::goal_curvature_shear(fx.store, &fx.t, *fx.frames, charts, spec, fx.hs)(q, g);
    return g;
  };
  real idle = 0, moved = 0, dH = 0, dD = 0;
  const std::vector<vec3> g0 = goal_of(0), g1 = goal_of(0.3);
  auto D = [](const vec3 &h) { return std::hypot(real(0.5) * (h[0] - h[2]), h[1]); };
  for (size_t e = 0; e < fx.hs.size(); ++e) {
    idle = std::max(idle, (g0[e] - fx.hs[e]).norm());
    moved += (g1[e] - fx.hs[e]).norm();
    dH = std::max(dH, std::abs((g1[e][0] + g1[e][2]) - (fx.hs[e][0] + fx.hs[e][2])));
    dD = std::max(dD, std::abs(D(g1[e]) - D(fx.hs[e])));
  }
  moved /= real(fx.hs.size());
  console_logger::info << "curvature shear: zero-shear drift " << idle << ", mean move at 0.3 "
                       << moved << ", worst mean-curvature change " << dH
                       << ", worst anisotropy change " << dD << std::endl;
  GAUDI_EXPECT(idle < real(1e-10));
  GAUDI_EXPECT(moved > real(1e-2));
  GAUDI_EXPECT(dH < real(1e-12));
  GAUDI_EXPECT(dD < real(1e-12));
}

/// What each knob does to one tensor, on elliptic, saddle and parabolic
/// tensors (signed `k0 >= k1`, ratio `r = |k_hi| / |k_lo|`).
GAUDI_TEST(weingarten_anisotropy_knob_identities) {
  auto H_of = [](const vec3 &g) { return real(0.5) * (g[0] + g[2]); };
  auto D_of = [](const vec3 &g) { return std::hypot(real(0.5) * (g[0] - g[2]), g[1]); };
  auto axis_of = [](const vec3 &g) { return real(0.5) * std::atan2(g[1], real(0.5) * (g[0] - g[2])); };
  auto ratio_of = [](const vec3 &g) {
    real k0, k1;
    hb::signed_principals(g, k0, k1);
    return std::max(std::abs(k0), std::abs(k1)) / std::min(std::abs(k0), std::abs(k1));
  };
  struct tcase {
    vec3 h;
    bool saddle;
  };
  const std::vector<tcase> cases{{hb::coeffs_at_angle(2.0, 0.5, 0.4), false},   // elliptic
                                 {hb::coeffs_at_angle(3.0, -1.0, -0.7), true},  // saddle
                                 {hb::coeffs_at_angle(1.0, -1.0, 0.3), true},   // minimal
                                 {hb::coeffs_at_angle(1.5, 0.0, 1.1), false},   // parabolic
                                 {hb::coeffs_at_angle(-0.3, -2.0, 0.2), false}}; // concave
  real worst = 0;
  int check = 0, which = -1;
  auto gap = [&](real a, real b) {
    ++check;
    if (std::abs(a - b) > worst) {
      worst = std::abs(a - b);
      which = check;
    }
  };
  auto knob = [](real ratio, bool log = false, bool relative = false) {
    hb::anisotropy_target s;
    s.ratio = ratio;
    s.log = log;
    s.relative = relative;
    return s;
  };
  for (const tcase &c : cases) {
    const vec3 &h = c.h;
    // Ratio 2.5: measured |k_hi| / |k_lo| is 2.5; elliptic keeps H, saddle
    // keeps |D|; axes kept.
    const vec3 g = hb::anisotropy_edge(h, h, knob(2.5));
    gap(ratio_of(g), 2.5);
    gap(c.saddle ? D_of(g) : H_of(g), c.saddle ? D_of(h) : H_of(h));
    gap(std::remainder(axis_of(g) - axis_of(h), real(M_PI)), 0);
    // The same target as a log, and below 1 read as its inverse.
    gap((hb::anisotropy_edge(h, h, knob(std::log(2.5), true)) - g).norm(), 0);
    gap((hb::anisotropy_edge(h, h, knob(0.4)) - g).norm(), 0);
    // Ratio 1: round at the same H, or a minimal saddle at the same |D|.
    const vec3 g1 = hb::anisotropy_edge(h, h, knob(1));
    gap(c.saddle ? H_of(g1) : D_of(g1), 0);
    gap(c.saddle ? D_of(g1) : H_of(g1), c.saddle ? D_of(h) : H_of(h));
    gap((hb::anisotropy_edge(h, h, knob(0, true)) - g1).norm(), 0);
    // Relative: a huge negative log stops at ratio 1; times 0.5 halves r.
    gap((hb::anisotropy_edge(h, h, knob(-50, true, true)) - g1).norm(), 0);
    if (std::isfinite(ratio_of(h)) && ratio_of(h) > real(2.5)) {
      const vec3 gh = hb::anisotropy_edge(h, h, knob(0.5, false, true));
      gap(ratio_of(gh), real(0.5) * ratio_of(h));
    }
    // Mean: H scaled, deviator kept.
    hb::anisotropy_target mn;
    mn.mean = 1.5;
    const vec3 gm = hb::anisotropy_edge(h, h, mn);
    gap(H_of(gm), real(1.5) * H_of(h));
    gap(D_of(gm), D_of(h));
  }
  // Ratio leaves a round tensor (no axes) alone.
  const vec3 um(0.8, 0.0, 0.8);
  gap((hb::anisotropy_edge(um, um, knob(3)) - um).norm(), 0);
  // The metrics are functions of r: q, log r, 1 / r.
  gap(hb::curvature_measure(hb::curvature_metric::aniso, 3.0, 1.0), 0.5);
  gap(hb::curvature_measure(hb::curvature_metric::aniso, 3.0, -1.0), 0.5);
  gap(hb::curvature_measure(hb::curvature_metric::log_ratio, 3.0, 1.0), std::log(real(3)));
  gap(hb::curvature_measure(hb::curvature_metric::ratio, -1.0, 3.0), real(1) / 3);
  console_logger::info << "anisotropy knob identities: worst gap " << worst << " (check " << which
                       << " of " << check << ")" << std::endl;
  GAUDI_EXPECT(worst < real(1e-10));
}

/// One degenerate triangle (two vertices welded, as a remesh can leave)
/// must not poison the coupling and memory rows: they leave its edges out,
/// so every Jacobian entry stays finite and of ordinary size.
GAUDI_TEST(weingarten_rows_skip_degenerate_faces) {
  std::vector<vec3> x;
  std::vector<std::array<int, 3>> tris;
  torus_tris(2.0, 0.7, 24, 12, x, tris);
  const kusama::stencil_table t =
      kusama::build_stencil_table(asawa::shell::triangle_walker(static_cast<int>(x.size()), tris));
  auto max_entry = [&](const std::vector<vec3> &xs) {
    const int n = static_cast<int>(xs.size());
    std::vector<vec3> xp = xs, v(n, vec3::Zero()), f = v, m(n, vec3::Ones());
    std::vector<vec3> hs = wg::edge_readings(t, xs), hv(hs.size(), vec3::Zero()), hf = hv,
                      hm(hs.size(), vec3::Constant(real(1e-8)));
    for (vec3 &h : hs)
      if (!h.allFinite())
        h = vec3::Zero();
    auto pos = hepworth::vec3_block::create(m, xp, v, f);
    auto store = hepworth::vec3_block::create(hm, hs, hv, hf);
    std::vector<hepworth::sim_block::ptr> blocks{pos, store};
    vecX q;
    pos->map_to_x(q);
    store->map_to_x(q);
    real worst = 0;
    bool finite = true;
    for (hepworth::projection_constraint::ptr c :
         {hepworth::projection_constraint::ptr(hb::weingarten_gn::create(&t, xs, real(1), blocks)),
          hepworth::projection_constraint::ptr(hb::dirichlet_memory_gn::create(
              &t, xs, real(1), std::vector<hepworth::sim_block::ptr>{pos}))}) {
      index_t id0 = 0;
      std::vector<trip> trips;
      c->fill_J(q, id0, trips);
      for (const trip &tr : trips) {
        finite = finite && std::isfinite(tr.value());
        worst = std::max(worst, std::abs(tr.value()));
      }
      finite = finite && std::isfinite(c->residual_sq(q));
    }
    return std::pair{finite, worst};
  };
  const auto [ok0, clean] = max_entry(x);
  // Weld the second vertex of the first face onto the first.
  std::vector<vec3> xb = x;
  xb[static_cast<size_t>(tris[0][1])] = xb[static_cast<size_t>(tris[0][0])];
  const kusama::weingarten::degenerate_set bad = kusama::weingarten::find_degenerate(t, xb);
  const auto [ok1, welded] = max_entry(xb);
  console_logger::info << "degenerate guard: clean max |J| " << clean << ", welded max |J| "
                       << welded << " (" << bad.n_faces << " degenerate faces, " << bad.n_edges
                       << " edges dropped)" << std::endl;
  GAUDI_EXPECT(ok0 && ok1);
  GAUDI_EXPECT(bad.n_faces >= 2);
  GAUDI_EXPECT(welded < real(100) * clean);
}

/// Directional curvature growth: the normal curvature along the axis
/// scales by `1 + x` and across it by `1 + y`; zero growth is the identity.
GAUDI_TEST(weingarten_curvature_growth_goal) {
  store_fixture fx;
  hb::curvature_growth_target zero, grow;
  grow.x = 0.3;
  grow.y = -0.1;
  real idle = 0, along = 0, across = 0;
  int n = 0;
  for (const vec3 &h0 : fx.hs) {
    idle = std::max(idle, (hb::curvature_growth_tensor(h0, zero) - h0).norm());
    const hb::shape_w ax = hb::shape_with_axes(hb::shape_w_from_coeffs(h0));
    if (!ax.ok || !(std::abs(ax.k1 - ax.k2) > real(1e-6)))
      continue;
    const hb::monge_axes m = hb::monge_axes_of(ax);
    const vec2 u(m.hx, m.hy), v(-m.hy, m.hx);
    auto kn = [](const vec3 &h, const vec2 &d) {
      return h[0] * d[0] * d[0] + real(2) * h[1] * d[0] * d[1] + h[2] * d[1] * d[1];
    };
    const vec3 g = hb::curvature_growth_tensor(h0, grow);
    along = std::max(along, std::abs(kn(g, u) - real(1.3) * kn(h0, u)));
    across = std::max(across, std::abs(kn(g, v) - real(0.9) * kn(h0, v)));
    ++n;
  }
  // `power`: 0 is the plain growth whatever k_bar; 1 with k_bar 1 scales
  // each fraction by that axis's |k|.
  hb::curvature_growth_target fed = grow, flat = grow;
  fed.power = 1;
  flat.power = 0;
  real plain = 0, fed_along = 0, fed_across = 0;
  for (const vec3 &h0 : fx.hs) {
    plain = std::max(plain, (hb::curvature_growth_tensor(h0, flat, real(3)) -
                             hb::curvature_growth_tensor(h0, grow)).norm());
    const hb::shape_w ax = hb::shape_with_axes(hb::shape_w_from_coeffs(h0));
    if (!ax.ok || !(std::abs(ax.k1 - ax.k2) > real(1e-6)))
      continue;
    const hb::monge_axes m = hb::monge_axes_of(ax);
    const vec2 u(m.hx, m.hy), v(-m.hy, m.hx);
    auto kn = [](const vec3 &h, const vec2 &d) {
      return h[0] * d[0] * d[0] + real(2) * h[1] * d[0] * d[1] + h[2] * d[1] * d[1];
    };
    const real ku = kn(h0, u), kv = kn(h0, v);
    const vec3 g = hb::curvature_growth_tensor(h0, fed, real(1));
    const real fu = std::max(std::abs(ku), real(1e-3)), fv = std::max(std::abs(kv), real(1e-3));
    fed_along = std::max(fed_along, std::abs(kn(g, u) - (real(1) + real(0.3) * fu) * ku));
    fed_across = std::max(fed_across, std::abs(kn(g, v) - (real(1) - real(0.1) * fv) * kv));
  }
  console_logger::info << "curvature growth: idle " << idle << ", along gap " << along
                       << ", across gap " << across << " over " << n << " edges; power 0 gap "
                       << plain << ", power 1 along " << fed_along << " across " << fed_across
                       << std::endl;
  GAUDI_EXPECT(n > 0);
  GAUDI_EXPECT(idle < real(1e-14));
  GAUDI_EXPECT(along < real(1e-12));
  GAUDI_EXPECT(across < real(1e-12));
  GAUDI_EXPECT(plain < real(1e-14));
  GAUDI_EXPECT(fed_along < real(1e-12));
  GAUDI_EXPECT(fed_across < real(1e-12));
}

/// The ring determinant row: idle at the frame-start store, its projection
/// lands every ring (elliptic, saddle and mixed) exactly on its held
/// `sum c det`, it is variational, and with a growth pin it keeps the rings'
/// totals that the growth alone would change.
GAUDI_TEST(weingarten_ring_det_conserves_and_redistributes) {
  store_fixture fx;
  std::vector<hepworth::projection_constraint::ptr> rows;
  hb::init_ring_det(fx.t, fx.x, fx.hs, real(1), fx.store, rows);
  GAUDI_ASSERT(!rows.empty());
  vecX q;
  fx.store->map_to_x(q);
  real idle = 0;
  for (auto &r : rows) {
    index_t id0 = 0;
    std::vector<trip> trips;
    r->fill_A(id0, trips);
    vecX Aq = vecX::Zero(id0), p = vecX::Zero(id0);
    for (const trip &t : trips)
      Aq[t.row()] += t.value() * q[t.col()];
    r->project(q, p);
    idle = std::max(idle, (Aq - p).norm() / std::max(Aq.norm(), real(1e-30)));
  }
  // A perturbed ring lands on its held total; count the mixed-sign rings.
  const std::vector<real> c = kusama::weingarten::coupling_factor(fx.t, fx.x);
  const auto rings = hb::vertex_edge_rings(fx.t);
  std::mt19937 rng(9);
  std::normal_distribution<real> nd(0, 0.2);
  real worst = 0;
  int mixed = 0;
  for (const auto &ring : rings) {
    std::vector<vec3> h0(ring.size()), h(ring.size());
    std::vector<real> cr(ring.size());
    bool pos = false, neg = false;
    for (size_t k = 0; k < ring.size(); ++k) {
      h0[k] = fx.hs[static_cast<size_t>(ring[k])];
      cr[k] = c[static_cast<size_t>(ring[k])];
      const real det = h0[k][0] * h0[k][2] - h0[k][1] * h0[k][1];
      pos |= det > 0;
      neg |= det < 0;
      h[k] = h0[k] + h0[k].norm() * vec3(nd(rng), nd(rng), nd(rng));
    }
    mixed += pos && neg;
    const real C0 = hb::ring_det_conserve::ring_det(h0, cr);
    const real got = hb::ring_det_conserve::ring_det(hb::ring_det_conserve::project_ring(h, cr, C0), cr);
    worst = std::max(worst, std::abs(got - C0) / std::max(std::abs(C0), real(1e-12)));
  }
  // Audit a perturbed store on a sample of rows.
  std::vector<vec3> hp = fx.hs;
  for (vec3 &h : hp)
    h += h.norm() * vec3(nd(rng), nd(rng), nd(rng));
  std::vector<vec3> hv(hp.size(), vec3::Zero()), hf = hv;
  auto store2 = hepworth::vec3_block::create(fx.hm, hp, hv, hf);
  std::vector<hepworth::projection_constraint::ptr> rows2;
  hb::init_ring_det(fx.t, fx.x, fx.hs, real(1), store2, rows2);
  std::vector<hepworth::projection_constraint::ptr> sample;
  for (size_t i = 0; i < rows2.size(); i += 17)
    sample.push_back(rows2[i]);
  vecX q2;
  store2->map_to_x(q2);
  real audit = 0;
  for (const auto &a : hb::audit_variational(sample, q2))
    audit = std::max(audit, a.max);
  // Growth pin with and without the det rows, store only.
  auto run = [&](bool det) {
    store_fixture f2;
    // The store block holds `f2.hs` by reference: keep the start.
    const std::vector<vec3> start = f2.hs;
    hb::curvature_growth_target g;
    g.x = 0.3;
    std::vector<hepworth::projection_constraint::ptr> cs{hb::edge_shape_fair::create(
        &f2.t, f2.x, hb::edge_fair_rows::pin, real(1), std::vector<hepworth::sim_block::ptr>{f2.store},
        hb::goal_curvature_growth(g, f2.hs), real(1), nullptr, true)};
    if (det)
      hb::init_ring_det(f2.t, f2.x, f2.hs, real(10), f2.store, cs);
    f2.run(cs, 5);
    const std::vector<vec3> &hs = f2.hs;
    real num = 0, den = 0;
    for (const auto &ring : rings) {
      std::vector<vec3> a(ring.size()), b(ring.size());
      std::vector<real> cr(ring.size());
      for (size_t k = 0; k < ring.size(); ++k) {
        a[k] = start[static_cast<size_t>(ring[k])];
        b[k] = hs[static_cast<size_t>(ring[k])];
        cr[k] = c[static_cast<size_t>(ring[k])];
      }
      num += std::abs(hb::ring_det_conserve::ring_det(b, cr) - hb::ring_det_conserve::ring_det(a, cr));
      den += std::abs(hb::ring_det_conserve::ring_det(a, cr));
    }
    return num / den;
  };
  const real drift_free = run(false), drift_held = run(true);
  console_logger::info << "ring det: idle " << idle << ", projected worst " << worst << " ("
                       << mixed << " mixed of " << rings.size() << " rings), audit " << audit
                       << ", growth drift free " << drift_free << " held " << drift_held << std::endl;
  GAUDI_EXPECT(idle < real(1e-12));
  GAUDI_EXPECT(worst < real(1e-9));
  GAUDI_EXPECT(mixed > 0);
  GAUDI_EXPECT(audit < real(1e-4));
  GAUDI_EXPECT(drift_held < real(0.2) * drift_free);
}

/// The vertex Gauss row: its vertex tensor is `vertex_tensor`'s, it is idle
/// at the baked `K0`, its projection lands each vertex (elliptic, saddle,
/// near-parabolic) exactly on `K0`, it is variational, and with a growth
/// pin it keeps the vertices' `K` that the growth alone would change.
GAUDI_TEST(weingarten_vertex_gauss_conserves) {
  store_fixture fx;
  const hb::vertex_star_charts ch = hb::build_vertex_star_charts(fx.t, fx.x);
  const kusama::weingarten::vertex_stars stars = kusama::weingarten::build_vertex_stars(fx.t, fx.x);
  const std::vector<vec3> hv = hb::vertex_tensors(ch, fx.hs);
  real chart_gap = 0;
  for (int v = 0; v < fx.t.n_verts; ++v) {
    const vec3 ref = kusama::weingarten::vertex_tensor(
        stars, v, [&](int i) { return fx.x[static_cast<size_t>(i)]; },
        [&](int e) { return fx.hs[static_cast<size_t>(e)]; }, ch.T[static_cast<size_t>(v)],
        ch.B[static_cast<size_t>(v)], ch.N[static_cast<size_t>(v)]);
    chart_gap = std::max(chart_gap, (ref - hv[static_cast<size_t>(v)]).norm());
  }
  const std::vector<real> K0 = hb::vertex_gauss(ch, fx.hs);
  std::vector<hepworth::projection_constraint::ptr> rows;
  hb::init_vertex_gauss(fx.t, fx.x, ch, K0, real(1), fx.store, rows);
  GAUDI_ASSERT(!rows.empty());
  vecX q;
  fx.store->map_to_x(q);
  real idle = 0;
  for (auto &r : rows) {
    index_t id0 = 0;
    std::vector<trip> trips;
    r->fill_A(id0, trips);
    vecX Aq = vecX::Zero(id0), p = vecX::Zero(id0);
    for (const trip &t : trips)
      Aq[t.row()] += t.value() * q[t.col()];
    r->project(q, p);
    idle = std::max(idle, (Aq - p).norm() / std::max(Aq.norm(), real(1e-30)));
  }
  // Perturbed spokes land on K0.
  const std::vector<real> c = kusama::weingarten::coupling_factor(fx.t, fx.x);
  std::mt19937 rng(11);
  std::normal_distribution<real> nd(0, 0.2);
  real Kmean = 0;
  for (real k : K0)
    Kmean += std::abs(k) / real(K0.size());
  real worst = 0;
  int ell = 0, sad = 0, flat = 0;
  for (size_t v = 0; v < ch.at.size(); ++v) {
    if (ch.at[v].empty())
      continue;
    std::vector<real> cs;
    std::vector<vec3> h;
    for (const auto &sp : ch.at[v]) {
      cs.push_back(c[static_cast<size_t>(sp.edge)]);
      const vec3 &h0 = fx.hs[static_cast<size_t>(sp.edge)];
      h.push_back(h0 + h0.norm() * vec3(nd(rng), nd(rng), nd(rng)));
    }
    hb::vertex_gauss_conserve row(ch.at[v], cs, K0[v], real(1),
                                  std::vector<hepworth::sim_block::ptr>{fx.store});
    const vec3 g = row.vertex_tensor(row.project_spokes(h));
    worst = std::max(worst, std::abs(g[0] * g[2] - g[1] * g[1] - K0[v]) /
                                std::max(std::abs(K0[v]), real(1e-3) * Kmean));
    (std::abs(K0[v]) < real(0.05) * Kmean ? flat : K0[v] > 0 ? ell : sad)++;
  }
  // Audit a perturbed store on a sample of rows.
  std::vector<vec3> hp = fx.hs;
  for (vec3 &h : hp)
    h += h.norm() * vec3(nd(rng), nd(rng), nd(rng));
  std::vector<vec3> hvel(hp.size(), vec3::Zero()), hf = hvel;
  auto store2 = hepworth::vec3_block::create(fx.hm, hp, hvel, hf);
  std::vector<hepworth::projection_constraint::ptr> rows2, sample;
  hb::init_vertex_gauss(fx.t, fx.x, ch, K0, real(1), store2, rows2);
  for (size_t i = 0; i < rows2.size(); i += 17)
    sample.push_back(rows2[i]);
  vecX q2;
  store2->map_to_x(q2);
  real audit = 0;
  for (const auto &a : hb::audit_variational(sample, q2))
    audit = std::max(audit, a.max);
  // Isotropic growth pin with and without the rows, store only.
  auto run = [&](bool held) {
    store_fixture f2;
    hb::curvature_growth_target g;
    g.x = g.y = 0.2;
    std::vector<hepworth::projection_constraint::ptr> cs{hb::edge_shape_fair::create(
        &f2.t, f2.x, hb::edge_fair_rows::pin, real(1), std::vector<hepworth::sim_block::ptr>{f2.store},
        hb::goal_curvature_growth(g, f2.hs), real(1), nullptr, true)};
    if (held)
      hb::init_vertex_gauss(f2.t, f2.x, ch, K0, real(10), f2.store, cs);
    f2.run(cs, 5);
    const std::vector<real> K = hb::vertex_gauss(ch, f2.hs);
    real num = 0, den = 0;
    for (size_t v = 0; v < K.size(); ++v) {
      num += std::abs(K[v] - K0[v]);
      den += std::abs(K0[v]);
    }
    return num / den;
  };
  const real drift_free = run(false), drift_held = run(true);
  console_logger::info << "vertex gauss: chart gap " << chart_gap << ", idle " << idle
                       << ", projected worst " << worst << " (" << ell << " elliptic, " << sad
                       << " saddle, " << flat << " near flat), audit " << audit
                       << ", growth drift free " << drift_free << " held " << drift_held << std::endl;
  GAUDI_EXPECT(chart_gap < real(1e-12));
  GAUDI_EXPECT(idle < real(1e-12));
  GAUDI_EXPECT(worst < real(1e-9));
  GAUDI_EXPECT(ell > 0 && sad > 0);
  GAUDI_EXPECT(audit < real(1e-4));
  GAUDI_EXPECT(drift_held < real(0.2) * drift_free);
}

/// Dupin rows with a target scale hold that multiple of the frame-start
/// gradients: at the frame-start store the residual is zero at scale 1 and
/// plus or minus the plain Dupin residual at scales 0 and 2.
GAUDI_TEST(weingarten_skew_rows_hold_a_scaled_gradient) {
  store_fixture fx;
  vecX q;
  fx.store->map_to_x(q);
  auto residual = [&](real scale) {
    auto row = hb::edge_shape_fair::create(
        &fx.t, fx.x, hb::edge_fair_rows::dupin, real(1),
        std::vector<hepworth::sim_block::ptr>{fx.store}, hb::goal_zero(), real(1), &fx.hs, false,
        scale);
    index_t id0 = 0;
    std::vector<trip> trips;
    row->fill_A(id0, trips);
    vecX Aq = vecX::Zero(id0), p = vecX::Zero(id0);
    for (const trip &t : trips)
      Aq[t.row()] += t.value() * q[t.col()];
    row->project(q, p);
    return vecX(Aq - p);
  };
  const vecX r0 = residual(0), r1 = residual(1), r2 = residual(2);
  console_logger::info << "skew rows: |r| at scale 0 " << r0.norm() << ", 1 " << r1.norm()
                       << ", 2 " << r2.norm() << std::endl;
  GAUDI_EXPECT(r0.norm() > real(1e-6));
  GAUDI_EXPECT(r1.norm() < real(1e-12) * r0.norm() + real(1e-14));
  GAUDI_EXPECT((r2 + r0).norm() < real(1e-12) * r0.norm() + real(1e-14));
}

/// The Willmore equalizer: one band is the plain shrink of the anisotropic
/// part; the screened blur leaves a flat constant field alone and damps
/// edge-to-edge noise; a free fine band leaves noise nearly unpenalized,
/// while the active coarse band still holds the constant part.
/// Dupin rows vanish (to discretization) on a torus, where both principal
/// curvatures are constant along their own lines, and not on a torus whose
/// tube radius varies around and along it.
GAUDI_TEST(weingarten_dupin_rows_hold_a_torus) {
  auto dupin_ratio = [](const std::vector<vec3> &x, const std::vector<std::array<int, 3>> &tris) {
    const kusama::stencil_table t =
        kusama::build_stencil_table(asawa::shell::triangle_walker(static_cast<int>(x.size()), tris));
    std::vector<vec3> h = wg::edge_readings(t, x);
    std::vector<vec3> hm(h.size(), vec3::Ones()), hv(h.size(), vec3::Zero()), hf = hv;
    auto blk = hepworth::vec3_block::create(hm, h, hv, hf);
    const std::vector<hepworth::sim_block::ptr> blocks{blk};
    auto dupin = hb::edge_shape_fair::create(&t, x, hb::edge_fair_rows::dupin, real(1), blocks,
                                             hb::goal_zero(), real(1), &h);
    auto pin = hb::edge_shape_fair::create(&t, x, hb::edge_fair_rows::pin, real(1), blocks,
                                           hb::goal_zero());
    vecX q;
    blk->map_to_x(q);
    auto rows_norm = [&](hb::edge_shape_fair &c) {
      index_t id0 = 0;
      std::vector<trip> trips;
      c.fill_A(id0, trips);
      wg::matS A(id0, q.size());
      A.setFromTriplets(trips.begin(), trips.end());
      return std::make_pair((A * q).norm(), id0);
    };
    const auto [rd, nd] = rows_norm(*dupin);
    const auto [rp, np] = rows_norm(*pin);
    console_logger::info << "dupin rows " << nd << " |A h| " << rd << " pin |h| " << rp << std::endl;
    return rd / std::max(rp, real(1e-30));
  };
  const real R = 2.0, r = 0.7;
  std::vector<vec3> x;
  std::vector<std::array<int, 3>> tris;
  torus_tris(R, r, 48, 24, x, tris);
  const real clean = dupin_ratio(x, tris);
  std::vector<vec3> bumpy = x;
  for (vec3 &p : bumpy) {
    const real phi = std::atan2(p[1], p[0]);
    const real rho = std::hypot(p[0], p[1]) - R;
    const real theta = std::atan2(p[2], rho);
    const real rt = r * (real(1) + real(0.2) * std::sin(real(3) * phi) * std::cos(theta));
    p = vec3((R + rt * std::cos(theta)) * std::cos(phi), (R + rt * std::cos(theta)) * std::sin(phi),
             rt * std::sin(theta));
  }
  const real rough = dupin_ratio(bumpy, tris);
  console_logger::info << "dupin residual / |h|: torus " << clean << " bumpy " << rough << std::endl;
  GAUDI_EXPECT(clean < real(0.25) * rough);
}

GAUDI_TEST(weingarten_willmore_equalizer) {
  auto goal_at = [](const hb::edge_goal_fn &goal, hepworth::sim_block::ptr store,
                    const std::vector<vec3> &h) {
    vecX q;
    store->map_to_x(q);
    for (size_t e = 0; e < h.size(); ++e)
      q.segment<3>(store->get_offset_idx(static_cast<index_t>(e))) = h[e];
    std::vector<vec3> g(h.size());
    goal(q, g);
    return g;
  };
  // Plain: one band at gain 0.2 shrinks the anisotropic part by 20%.
  store_fixture fx;
  hb::willmore_target plain;
  plain.bands = {{0.0, 0.2}};
  const std::vector<vec3> g = goal_at(hb::goal_willmore_eq(fx.store, fx.t, fx.x, fx.hs, plain),
                                      fx.store, fx.hs);
  real worst = 0;
  for (size_t e = 0; e < g.size(); ++e) {
    const vec3 &h = fx.hs[e];
    const real H = real(0.5) * (h[0] + h[2]), d = real(0.5) * (h[0] - h[2]);
    const vec3 want(H + real(0.8) * d, real(0.8) * h[1], H - real(0.8) * d);
    worst = std::max(worst, (g[e] - want).norm());
  }
  // Flat grid: constant field and edge noise.
  std::vector<vec3> xg;
  std::vector<std::array<int, 3>> tris;
  unit_grid(9, xg, tris);
  const kusama::stencil_table tg =
      kusama::build_stencil_table(asawa::shell::triangle_walker(static_cast<int>(xg.size()), tris));
  const std::vector<vec3> hc = sample_edges_planar(tg, xg, [](const vec3 &) {
    mat3 S;
    S << 1.0, 0.3, 0, 0.3, -0.5, 0, 0, 0, 0;
    return S;
  });
  std::mt19937 rng(7);
  std::uniform_real_distribution<real> u(-1, 1);
  std::vector<vec3> noise(hc.size());
  for (vec3 &n : noise)
    n = vec3(u(rng), u(rng), u(rng));
  auto norm = [](const std::vector<vec3> &a) {
    real s = 0;
    for (const vec3 &v : a)
      s += v.squaredNorm();
    return std::sqrt(s);
  };
  const hb::edge_tensor_blur blur(tg, xg, {real(3)});
  std::vector<vec3> d_const = blur.apply(0, hc), d_noise = blur.apply(0, noise);
  for (size_t e = 0; e < hc.size(); ++e)
    d_const[e] -= hc[e];
  const real keep = norm(d_const) / norm(hc), damp = norm(d_noise) / norm(noise);
  // Balloon preset: fine band free, coarse band active. Noise added to the
  // store moves the residual (goal - S) far less than the noise itself.
  std::vector<vec3> hm(hc.size(), vec3::Constant(real(1e-8))), hv(hc.size(), vec3::Zero()),
      hf = hv, hs = hc;
  auto store = hepworth::vec3_block::create(hm, hs, hv, hf);
  hb::willmore_target balloon;
  balloon.bands = {{0.0, 0.0}, {3.0, 0.2}};
  const hb::edge_goal_fn gb = hb::goal_willmore_eq(store, tg, xg, hc, balloon);
  std::vector<vec3> noisy = hc;
  for (size_t e = 0; e < hc.size(); ++e)
    noisy[e] += noise[e];
  const std::vector<vec3> g0 = goal_at(gb, store, hc), g1 = goal_at(gb, store, noisy);
  std::vector<vec3> r0(hc.size()), r1(hc.size());
  for (size_t e = 0; e < hc.size(); ++e) {
    r0[e] = g0[e] - hc[e];
    r1[e] = g1[e] - noisy[e];
  }
  // The goal is taken at the frame-start store and held (a fixed pin, so the
  // row is variational): it does not depend on the live store, so within a
  // frame the free bands are held at their frame-start values too; each
  // frame re-reads them.
  std::vector<vec3> dg(hc.size());
  for (size_t e = 0; e < hc.size(); ++e)
    dg[e] = g1[e] - g0[e];
  const real moves = norm(dg) / norm(g0);
  console_logger::info << "willmore equalizer: plain vs shrink " << worst
                       << ", blur keeps a flat constant to " << keep << ", damps noise to "
                       << damp << ", goal change with the live store " << moves << std::endl;
  GAUDI_EXPECT(worst < real(1e-12));
  GAUDI_EXPECT(keep < real(1e-6));
  GAUDI_EXPECT(damp < real(0.3));
  GAUDI_EXPECT(moves < real(1e-12));
}

/// Shape ARAP on the torus store: idle at rest, free under a uniform turn of
/// the lines, and it restores the frame-start anisotropy and mean.
GAUDI_TEST(weingarten_shape_arap_goal) {
  store_fixture fx;
  const std::vector<vec3> h0 = fx.hs;
  auto goal = hb::goal_shape_arap(fx.store, &fx.t, h0, true);
  vecX q;
  fx.store->map_to_x(q);
  auto run = [&](const std::vector<vec3> &h) {
    for (size_t e = 0; e < h.size(); ++e)
      q.segment<3>(fx.store->get_offset_idx(static_cast<index_t>(e))) = h[e];
    std::vector<vec3> g(h.size());
    goal(q, g);
    return g;
  };
  auto gap = [](const std::vector<vec3> &a, const std::vector<vec3> &b) {
    real w = 0, n = 0;
    for (size_t e = 0; e < a.size(); ++e) {
      w = std::max(w, (a[e] - b[e]).norm());
      n = std::max(n, b[e].norm());
    }
    return w / std::max(n, real(1e-30));
  };
  // At rest the goal is the rest.
  const real idle = gap(run(h0), h0);
  // Every edge's lines turned by 0.4: the goal follows the turn.
  std::vector<vec3> turned = h0;
  for (vec3 &h : turned) {
    const real H = real(0.5) * (h[0] + h[2]);
    const std::complex<real> z =
        std::complex<real>(real(0.5) * (h[0] - h[2]), h[1]) * std::polar(real(1), real(0.8));
    h = vec3(H + z.real(), z.imag(), H - z.real());
  }
  const real turn = gap(run(turned), turned);
  // Anisotropic part and mean scaled: the goal is the rest again.
  std::vector<vec3> scaled = h0;
  for (vec3 &h : scaled)
    h *= real(1.4);
  const real restore = gap(run(scaled), h0);
  console_logger::info << "shape arap: idle " << idle << ", uniform turn " << turn
                       << ", scaled back to rest " << restore << std::endl;
  GAUDI_EXPECT(idle < real(1e-12));
  GAUDI_EXPECT(turn < real(1e-12));
  GAUDI_EXPECT(restore < real(1e-12));
}

/// `directed_tilt` turns the frame-start ring by `gain * metric` and holds
/// that target for the frame: about `N` for spin, about the chart's principal
/// axis for the others. Twist, bend and the gradient-signed rigid tilt do
/// not depend on which way the axis points; the propagated sign agrees
/// across neighbors.
GAUDI_TEST(weingarten_directed_tilt_turns_the_ring) {
  std::vector<vec3> x{vec3::Zero()};
  for (int k = 0; k < 6; ++k)
    x.push_back(vec3(std::cos(k * M_PI / 3), std::sin(k * M_PI / 3), 0));
  const int n = static_cast<int>(x.size());
  std::vector<vec3> v(n, vec3::Zero()), f = v, m(n, vec3::Ones());
  auto pos = hepworth::vec3_block::create(m, x, v, f);
  // k0 = 1 on T = x, k1 on B = y; k1 grows with x on the ring, so the
  // anisotropy falls toward +x.
  auto charts_with = [&](const vec3 &T, const vec3 &B) {
    std::vector<hb::vertex_chart> c(n);
    for (int i = 0; i < n; ++i) {
      const real k1 = real(0.2) + real(0.1) * x[static_cast<size_t>(i)][0];
      // The tensor in (T, B): k0 on world x, k1 on world y.
      const real cx = T.dot(vec3::UnitX()), cy = T.dot(vec3::UnitY());
      const real H = real(0.5) * (1 + k1), d = real(0.5) * (1 - k1);
      const real c2 = cx * cx - cy * cy, s2 = 2 * cx * cy;
      c[static_cast<size_t>(i)] = {T, B, vec3::UnitZ(), vec3(H + d * c2, d * s2, H - d * c2), true};
    }
    return c;
  };
  const std::vector<hb::vertex_chart> charts = charts_with(vec3::UnitX(), vec3::UnitY());
  const std::vector<hb::vertex_chart> flipped = charts_with(-vec3::UnitX(), -vec3::UnitY());
  std::vector<index_t> ids{0, 1, 2, 3, 4, 5, 6};
  std::vector<real> wk(6, real(1));
  const real gain = 0.3, w = 0.5;
  const real theta = gain * (real(1) - real(0.2)) / (real(1) + real(0.2));
  auto targets = [&](hb::tilt_target spec, const std::vector<hb::vertex_chart> &ch,
                     std::shared_ptr<hb::tilt_orientation> field = nullptr) {
    auto c = hb::directed_tilt::create(ids, wk, x, spec, w,
                                       std::vector<hepworth::sim_block::ptr>{pos}, &ch, field);
    index_t id0 = 0;
    std::vector<trip> trips;
    c->fill_A(id0, trips);
    vecX q;
    pos->map_to_x(q);
    vecX p = vecX::Zero(id0);
    c->project(q, p);
    // Moving the vertices must not move a held target.
    vecX q2 = q + vecX::Random(q.size()) * real(0.1);
    vecX p2 = vecX::Zero(id0);
    c->project(q2, p2);
    GAUDI_EXPECT((p2 - p).norm() < real(1e-14));
    return p;
  };
  real worst = 0;
  const vec3 a = vec3::UnitY(), b = vec3::UnitZ().cross(a);
  for (hb::tilt_mode mode : {hb::tilt_mode::spin, hb::tilt_mode::twist, hb::tilt_mode::bend}) {
    hb::tilt_target spec;
    spec.mode = mode;
    spec.gain = gain;
    const vecX p = targets(spec, charts);
    for (int k = 0; k < 6; ++k) {
      const vec3 d = x[static_cast<size_t>(k + 1)];
      real angle = theta;
      vec3 axis = a;
      if (mode == hb::tilt_mode::spin)
        axis = vec3::UnitZ();
      else if (mode == hb::tilt_mode::twist)
        angle = theta * d.dot(a);
      else
        angle = theta * d.dot(b);
      const vec3 want = w * (Eigen::AngleAxis<real>(angle, axis) * d);
      worst = std::max(worst, (p.segment<3>(3 * k) - want).norm());
    }
    // The same ring with the chart turned a half turn about N.
    worst = std::max(worst, (targets(spec, flipped) - p).norm());
  }
  // Gradient-signed rigid: tilts about the k1 axis, and the flipped chart
  // gives the same targets.
  hb::tilt_target rigid;
  rigid.gain = gain;
  const vecX pr = targets(rigid, charts);
  real moved = 0;
  for (int k = 0; k < 6; ++k)
    moved = std::max(moved, (pr.segment<3>(3 * k) - w * x[static_cast<size_t>(k + 1)]).norm());
  const real flip_gap = (targets(rigid, flipped) - pr).norm();
  // Propagated sign on a 3-vertex chain whose charts alternate a half turn:
  // every neighboring pair agrees after orientation.
  std::vector<hb::vertex_chart> chain = charts_with(vec3::UnitX(), vec3::UnitY());
  chain[1] = flipped[1];
  auto field = std::make_shared<hb::tilt_orientation>();
  field->nbr = {{1}, {0, 2}, {1}, {}, {}, {}, {}};
  field->orient(&chain);
  real agree = 1;
  for (int i : {0, 1})
    agree = std::min(agree, (field->sign[static_cast<size_t>(i)] *
                             hb::tilt_orientation::chart_axis(&chain, i, hb::principal_axis::k1))
                                .dot(field->sign[static_cast<size_t>(i + 1)] *
                                     hb::tilt_orientation::chart_axis(&chain, i + 1,
                                                                      hb::principal_axis::k1)));
  console_logger::info << "directed_tilt: worst target gap " << worst << " at theta " << theta
                       << ", rigid moved " << moved << " flip gap " << flip_gap
                       << ", propagated neighbor agreement " << agree << std::endl;
  GAUDI_EXPECT(worst < real(1e-12));
  GAUDI_EXPECT(moved > real(1e-3));
  GAUDI_EXPECT(flip_gap < real(1e-12));
  GAUDI_EXPECT(agree > real(0.99));
}

/// Coupled frames on the torus, as the demo runs them: the massless store
/// starts each frame at the live reading, `weingarten_gn` (the read-out row,
/// weight `w_c`) and the memory couple positions to it, and `rows(t, x,
/// store, h0)` adds the target rows. Returns the positions after each frame.
template <typename Rows>
inline std::vector<std::vector<vec3>> run_edit_frames(std::vector<vec3> x0,
                                                     const std::vector<std::array<int, 3>> &tris,
                                                     int frames, real w_c, Rows rows) {
  const kusama::stencil_table t =
      kusama::build_stencil_table(asawa::shell::triangle_walker(static_cast<int>(x0.size()), tris));
  const int n = static_cast<int>(x0.size());
  const real h = 0.1;
  std::vector<vec3> x = std::move(x0), v(n, vec3::Zero()), f = v, m(n, vec3::Ones());
  auto pos = hepworth::vec3_block::create(m, x, v, f);
  std::vector<vec3> hs, hv, hf, hm;
  auto store = hepworth::vec3_block::create(hm, hs, hv, hf);
  std::vector<hepworth::sim_block::ptr> blocks{pos, store};
  hb::projection_solver solver;
  solver.gn_refresh_every = 0;
  solver.gn_backtracks = 8;
  std::vector<std::vector<vec3>> out{x};
  for (int k = 0; k < frames; ++k) {
    const std::vector<wg::cr_edge_rows> D = wg::build_cr_derivative(t, x);
    hs = wg::edge_readings(t, x);
    hv.assign(hs.size(), vec3::Zero());
    hf.assign(hs.size(), vec3::Zero());
    hm.assign(hs.size(), vec3::Constant(real(1e-8)));
    std::vector<hepworth::projection_constraint::ptr> cs{
        hb::weingarten_gn::create(&t, x, w_c, blocks),
        hb::dirichlet_memory_gn::create(&t, x, real(1), std::vector<hepworth::sim_block::ptr>{pos})};
    for (auto &c : rows(t, x, store, hs))
      cs.push_back(c);
    solver.set_constraints(cs);
    solver.step(blocks, h, real(1), 30);
    out.push_back(x);
  }
  return out;
}

inline real mean_displacement(const std::vector<vec3> &a, const std::vector<vec3> &b) {
  real d = 0;
  for (size_t i = 0; i < a.size(); ++i)
    d += (a[i] - b[i]).norm();
  return d / real(a.size());
}

GAUDI_TEST(weingarten_shape_edits_drive_the_mesh) {
  std::vector<vec3> x0;
  std::vector<std::array<int, 3>> tris;
  torus_tris(2.0, 0.7, 24, 12, x0, tris);
  const kusama::stencil_table t =
      kusama::build_stencil_table(asawa::shell::triangle_walker(static_cast<int>(x0.size()), tris));
  using rows = hb::edge_fair_rows;
  using cvec = std::vector<hepworth::projection_constraint::ptr>;
  const real w_c = 1e3;
  // Target weights beside the read-out's: `w_c` is an even blend.
  auto knob = [w_c](rows r, hb::anisotropy_target spec, real w = real(-1)) {
    const real wt = w > 0 ? w : w_c;
    return [r, spec, wt](const kusama::stencil_table &tt, const std::vector<vec3> &x,
                         hepworth::sim_block::ptr store, const std::vector<vec3> &h0) {
      return cvec{hb::edge_shape_fair::create(&tt, x, r, wt,
                                              std::vector<hepworth::sim_block::ptr>{store},
                                              hb::goal_anisotropy(store, spec, h0))};
    };
  };
  {
    // Zero knobs: nothing moves.
    const auto xs = run_edit_frames(x0, tris, 3, w_c, knob(rows::h1, hb::anisotropy_target{}));
    const real d = mean_displacement(xs.back(), x0);
    console_logger::info << "edits: zero knobs, mean displacement over 3 frames " << d << std::endl;
    GAUDI_EXPECT(d < real(1e-8));
  }
  // Mean q on elliptic edges (on saddles q = |H| / |D| rises as D shrinks).
  auto mean_aniso = [&](const std::vector<vec3> &x) {
    real s = 0;
    int c = 0;
    for (const vec3 &h : wg::edge_readings(t, x)) {
      real k0, k1;
      hb::signed_principals(h, k0, k1);
      if (k0 * k1 > 0) {
        s += hb::edge_q(h);
        ++c;
      }
    }
    return s / real(std::max(c, 1));
  };
  // Mean |log r - target| (log r capped at log 1000).
  auto mean_log_gap = [&](const std::vector<vec3> &x, real target) {
    real s = 0;
    const std::vector<vec3> hs = wg::edge_readings(t, x);
    for (const vec3 &h : hs) {
      real k0, k1;
      hb::signed_principals(h, k0, k1);
      s += std::abs(hb::curvature_measure(hb::curvature_metric::log_ratio, k0, k1) - target);
    }
    return s / real(std::max<size_t>(hs.size(), 1));
  };
  {
    // Ratio times 0.8 per frame: the measured anisotropy drops. (Compounding
    // far past what a torus can carry runs away here: the harness holds the
    // positions with the memory alone.)
    hb::anisotropy_target spec;
    spec.ratio = 0.8;
    spec.relative = true;
    const auto xs = run_edit_frames(x0, tris, 3, w_c, knob(rows::h1, spec));
    for (size_t k = 1; k < xs.size(); ++k)
      console_logger::info << "  frame " << k << " aniso " << mean_aniso(xs[k]) << " |dx| "
                           << mean_displacement(xs[k], xs[k - 1]) << std::endl;
    const real a0 = mean_aniso(x0), a1 = mean_aniso(xs.back());
    console_logger::info << "edits: ratio x0.8 per frame, mean aniso " << a0 << " -> " << a1
                         << ", displacement " << mean_displacement(xs.back(), x0) << std::endl;
    GAUDI_EXPECT(a1 < real(0.95) * a0);
  }
  {
    // Ratio 1.5: the measured log ratio moves toward log 1.5. (The harness
    // keeps full velocity between frames, so past the target it overshoots.)
    hb::anisotropy_target spec;
    spec.ratio = 1.5;
    const real L = std::log(real(1.5));
    const auto xs = run_edit_frames(x0, tris, 2, w_c, knob(rows::h1, spec));
    for (size_t k = 1; k < xs.size(); ++k)
      console_logger::info << "  frame " << k << " |log r - log 1.5| " << mean_log_gap(xs[k], L)
                           << " |dx| " << mean_displacement(xs[k], xs[k - 1]) << std::endl;
    const real g0 = mean_log_gap(x0, L), g1 = mean_log_gap(xs.back(), L);
    console_logger::info << "edits: ratio 1.5, mean |log r - log 1.5| " << g0 << " -> " << g1
                         << ", displacement " << mean_displacement(xs.back(), x0) << std::endl;
    GAUDI_EXPECT(g1 < real(0.8) * g0);
  }
  {
    // h1 rows against Dirichlet rows on the same edit, across weights.
    hb::anisotropy_target spec;
    spec.ratio = 0.8;
    spec.relative = true;
    for (real w : {real(0.1) * w_c, w_c, real(10) * w_c}) {
      const real d_h1 =
          mean_displacement(run_edit_frames(x0, tris, 1, w_c, knob(rows::h1, spec, w))[1], x0);
      const real d_dir =
          mean_displacement(run_edit_frames(x0, tris, 1, w_c, knob(rows::dirichlet, spec, w))[1], x0);
      console_logger::info << "edits: ratio x0.8 one frame, w " << w << " displacement h1 " << d_h1
                           << " dirichlet " << d_dir << std::endl;
      // A relative ratio edit varies over the torus, so the
      // Dirichlet rows see it as well as the pin part of h1.
      GAUDI_EXPECT(d_h1 > real(1e-2) && d_dir > real(1e-2));
    }
  }
}

/// With the positions held, the least-squares solve puts the store at the
/// combination of the read-out row and a target row with the same factor,
/// `S = (w_s^2 R + w_t^2 g) / (w_s^2 + w_t^2)` per edge (rows `w * factor`);
/// with no target row it is the reading.
GAUDI_TEST(weingarten_store_is_the_blend_of_its_rows) {
  std::vector<vec3> x0;
  std::vector<std::array<int, 3>> tris;
  torus_tris(2.0, 0.7, 24, 12, x0, tris);
  const kusama::stencil_table t =
      kusama::build_stencil_table(asawa::shell::triangle_walker(static_cast<int>(x0.size()), tris));
  const std::vector<vec3> R = wg::edge_readings(t, x0);
  std::vector<vec3> g = R;
  for (vec3 &v : g)
    v *= real(1.2);
  const real h = 0.1;
  for (auto [w_s, w_t] : {std::pair{real(1e3), real(0)}, std::pair{real(1e3), real(1e3)},
                          std::pair{real(1e3), real(3e3)}}) {
    const int n = static_cast<int>(x0.size());
    std::vector<vec3> x = x0, v(n, vec3::Zero()), f = v, m(n, vec3::Ones());
    auto pos = hepworth::vec3_block::create(m, x, v, f);
    std::vector<vec3> hs = R, hv(R.size(), vec3::Zero()), hf = hv,
                      hm(R.size(), vec3::Constant(real(1e-8)));
    auto store = hepworth::vec3_block::create(hm, hs, hv, hf);
    std::vector<hepworth::sim_block::ptr> blocks{pos, store};
    std::vector<hepworth::projection_constraint::ptr> cs{hb::weingarten_gn::create(&t, x, w_s, blocks)};
    for (int i = 0; i < n; ++i)
      cs.push_back(hb::pinned::create(std::vector<index_t>{i}, x0[i], real(1e6),
                                      std::vector<hepworth::sim_block::ptr>{pos}));
    if (w_t > 0)
      cs.push_back(hb::edge_shape_fair::create(&t, x, hb::edge_fair_rows::pin, w_t,
                                               std::vector<hepworth::sim_block::ptr>{store},
                                               hb::goal_fixed(g)));
    hb::projection_solver solver;
    solver.gn_refresh_every = 0;
    solver.gn_backtracks = 8;
    solver.set_constraints(cs);
    solver.step(blocks, h, real(1), 30);
    real err = 0, nrm = 0;
    for (size_t e = 0; e < R.size(); ++e) {
      const vec3 want = (w_s * w_s * R[e] + w_t * w_t * g[e]) / (w_s * w_s + w_t * w_t);
      err += (hs[e] - want).squaredNorm();
      nrm += want.squaredNorm();
    }
    const real rel = std::sqrt(err / nrm);
    console_logger::info << "store blend: w_s " << w_s << " w_t " << w_t
                         << " |S - (w_s^2 R + w_t^2 g) / (w_s^2 + w_t^2)| / |.| " << rel
                         << ", positions moved " << mean_displacement(x, x0) << std::endl;
    GAUDI_EXPECT(rel < real(1e-3));
  }
}

/// Positions feel a target in series with the read-out,
/// `w_s^2 w_t^2 / (w_s^2 + w_t^2)`: with the read-out weak against inertia,
/// tripling `w_t` from `w_s` raises the motion by about `0.9 / 0.5 = 1.8`.
GAUDI_TEST(weingarten_target_reaches_the_mesh_at_the_series_weight) {
  std::vector<vec3> x0;
  std::vector<std::array<int, 3>> tris;
  torus_tris(2.0, 0.7, 24, 12, x0, tris);
  const real w_s = 1;
  auto run = [&](real w_t) {
    const auto xs = run_edit_frames(
        x0, tris, 1, w_s,
        [w_t](const kusama::stencil_table &tt, const std::vector<vec3> &x,
              hepworth::sim_block::ptr store, const std::vector<vec3> &h0) {
          std::vector<vec3> g = h0;
          for (vec3 &v : g)
            v *= real(1.05);
          return std::vector<hepworth::projection_constraint::ptr>{hb::edge_shape_fair::create(
              &tt, x, hb::edge_fair_rows::pin, w_t,
              std::vector<hepworth::sim_block::ptr>{store}, hb::goal_fixed(g))};
        });
    return mean_displacement(xs[1], x0);
  };
  const real d1 = run(w_s), d3 = run(real(3) * w_s);
  console_logger::info << "series weight: displacement w_t = w_s " << d1 << ", w_t = 3 w_s " << d3
                       << ", ratio " << d3 / d1 << " (predicted 1.8)" << std::endl;
  GAUDI_EXPECT(d3 / d1 > real(1.6) && d3 / d1 < real(2.0));
}

GAUDI_TEST(weingarten_minimizer_smooths_the_mesh) {
  std::vector<vec3> x0;
  std::vector<std::array<int, 3>> tris;
  torus_tris(2.0, 0.7, 24, 12, x0, tris);
  const kusama::stencil_table t =
      kusama::build_stencil_table(asawa::shell::triangle_walker(static_cast<int>(x0.size()), tris));
  std::mt19937 rng(3);
  std::normal_distribution<real> nd(0, 0.01);
  for (vec3 &p : x0)
    p += vec3(nd(rng), nd(rng), nd(rng));
  std::shared_ptr<const hb::edge_frames> fr = hb::build_edge_frames(t, x0);
  auto roughness = [&](const std::vector<vec3> &x) {
    const std::vector<vec3> h = wg::edge_readings(t, x);
    real s = 0;
    int c = 0;
    for (size_t e = 0; e < h.size(); ++e)
      for (int j : fr->nbr[e]) {
        s += std::abs((h[e][0] + h[e][2]) - (h[static_cast<size_t>(j)][0] + h[static_cast<size_t>(j)][2]));
        ++c;
      }
    return s / std::max(c, 1);
  };
  const auto xs = run_edit_frames(
      x0, tris, 5, real(1e3),
      [](const kusama::stencil_table &tt, const std::vector<vec3> &x, hepworth::sim_block::ptr store,
         const std::vector<vec3> &) {
        return std::vector<hepworth::projection_constraint::ptr>{hb::edge_shape_fair::create(
            &tt, x, hb::edge_fair_rows::pin, real(1e3), std::vector<hepworth::sim_block::ptr>{store},
            hb::goal_principal_mean(store, hb::build_edge_frames(tt, x)))};
      });
  const real r0 = roughness(x0), r1 = roughness(xs.back());
  console_logger::info << "minimizer: mean-curvature roughness " << r0 << " -> " << r1
                       << ", displacement " << mean_displacement(xs.back(), x0) << std::endl;
  GAUDI_EXPECT(r1 < real(0.8) * r0);
}

// ------------------------------------------------------------------- M3

/// The coupling row is `w` times its geometry: on the store its coefficients
/// are exactly `-w c_e` (`coupling_factor`).
GAUDI_TEST(weingarten_gn_rows_are_weight_times_geometry) {
  std::vector<vec3> x0;
  std::vector<std::array<int, 3>> tris;
  torus_tris(2.0, 0.7, 16, 8, x0, tris);
  const kusama::stencil_table t =
      kusama::build_stencil_table(asawa::shell::triangle_walker(static_cast<int>(x0.size()), tris));
  const int n = static_cast<int>(x0.size());
  std::vector<vec3> x = x0, v(n, vec3::Zero()), f = v, m(n, vec3::Ones());
  auto pos = hepworth::vec3_block::create(m, x, v, f);
  std::vector<vec3> hs = wg::edge_readings(t, x0), hv(hs.size(), vec3::Zero()), hf = hv,
                    hm(hs.size(), vec3::Constant(real(1e-8)));
  auto store = hepworth::vec3_block::create(hm, hs, hv, hf);
  vecX q;
  pos->map_to_x(q);
  store->map_to_x(q);
  const real w = 2.5;
  auto c = hb::weingarten_gn::create(&t, x0, w, std::vector<hepworth::sim_block::ptr>{pos, store});
  index_t id0 = 0;
  std::vector<trip> trips;
  c->fill_J(q, id0, trips);
  const std::vector<real> ce = wg::coupling_factor(t, x0);
  real worst = 0;
  int seen = 0;
  for (const trip &tr : trips) {
    const index_t col = tr.col() - store->get_offset_idx(0);
    if (col < 0)
      continue;
    const int e = static_cast<int>(col / 3);
    worst = std::max(worst, std::abs(tr.value() + w * ce[static_cast<size_t>(e)]));
    ++seen;
  }
  console_logger::info << "coupling store coefficients vs -w c_e: worst " << worst << " over "
                       << seen << std::endl;
  GAUDI_EXPECT(seen == 3 * t.n_edges);
  GAUDI_EXPECT(worst < real(1e-14));
}

/// At equal weights the coupling row and `bend` pull on positions alike: the
/// mean `diag(A^T A)` over position columns at `w = 1` (demo torus).
GAUDI_TEST(weingarten_coupling_matches_bend_at_equal_weight) {
  for (auto [R, r, nu, nv] : {std::tuple{1.0, 0.33, 48, 24}, std::tuple{2.0, 0.7, 24, 12}}) {
  asawa::shell::shell::ptr M = asawa::shell::load_torus(R, r, nu, nv);
  std::vector<vec3> &x = asawa::get_vec_data(*M, 0);
  const int n = static_cast<int>(x.size());
  const kusama::stencil_table t = kusama::build_stencil_table(asawa::shell::shell_walker(*M));
  std::vector<vec3> v(n, vec3::Zero()), f = v, m(n, vec3::Ones());
  auto pos = hepworth::vec3_block::create(m, x, v, f);
  std::vector<vec3> hs = wg::edge_readings(t, x), hv(hs.size(), vec3::Zero()), hf = hv,
                    hm(hs.size(), vec3::Constant(real(1e-8)));
  auto store = hepworth::vec3_block::create(hm, hs, hv, hf);
  vecX q;
  pos->map_to_x(q);
  store->map_to_x(q);
  auto mean_pos_diag = [&](const std::vector<trip> &trips) {
    std::map<index_t, real> d;
    for (const trip &tr : trips)
      if (tr.col() < 3 * n)
        d[tr.col()] += tr.value() * tr.value();
    real sum = 0;
    for (const auto &[c, val] : d)
      sum += val;
    return d.empty() ? real(0) : sum / real(d.size());
  };
  std::vector<hepworth::projection_constraint::ptr> bend;
  hb::init_bending(*M, bend, x, real(1), {pos});
  std::vector<trip> tb;
  index_t id0 = 0;
  for (auto &c : bend)
    c->fill_A(id0, tb);
  auto coupling = hb::weingarten_gn::create(&t, x, real(1),
                                            std::vector<hepworth::sim_block::ptr>{pos, store});
  std::vector<trip> tc;
  id0 = 0;
  coupling->fill_J(q, id0, tc);
  const real kb = mean_pos_diag(tb), kc = mean_pos_diag(tc);
  console_logger::info << "w = 1 position stiffness on torus (" << R << ", " << r << ", " << nu
                       << " x " << nv << "): bend " << kb << " coupling " << kc << " ratio "
                       << kc / kb << std::endl;
  GAUDI_EXPECT(kc / kb > real(1) / 3 && kc / kb < real(3));
  }
}

GAUDI_TEST(weingarten_gn_idle_at_rest_on_shell_store) {
  asawa::shell::shell::ptr M = asawa::shell::load_torus(1.0, 0.33, 16, 8);
  auto cr = std::make_shared<hb::shell_cr_block>(M);
  hb::solver_context ctx;
  cr->prepare(ctx);
  std::vector<hepworth::sim_block::ptr> dofs;
  cr->emit_dof_blocks(dofs);
  GAUDI_ASSERT(cr->shape_dof() && cr->table.n_edges == cr->topology.nE);
  for (int e = 0; e < cr->topology.nE; ++e) {
    const asawa::shell::CornerId c = cr->topology.dof_to_corner[static_cast<size_t>(e)];
    std::set<int> a{static_cast<int>(M->vert(c)), static_cast<int>(M->vert(M->other(c)))};
    bool found = false;
    for (int f = 0; f < cr->table.n_faces && !found; ++f)
      for (int i = 0; i < 3; ++i)
        if (cr->table.face_edges[f][i] == e)
          found = a == std::set<int>{cr->table.ids[f][i], cr->table.ids[f][(i + 1) % 3]};
    GAUDI_ASSERT(found); // table edge ids are the CR edge dofs
  }
  std::vector<vec3> &x = asawa::get_vec_data(*M, 0);
  const int n = static_cast<int>(x.size());
  std::vector<vec3> v(n, vec3::Zero()), f = v, m(n, vec3::Ones());
  auto pos = hepworth::vec3_block::create(m, x, v, f);
  vecX q;
  pos->map_to_x(q);
  cr->shape_dof()->map_to_x(q);
  auto c = hb::weingarten_gn::create(&cr->table, x, real(10),
                                     std::vector<hepworth::sim_block::ptr>{pos, cr->shape_dof()});
  auto mem = hb::dirichlet_memory_gn::create(&cr->table, x, real(1),
                                             std::vector<hepworth::sim_block::ptr>{pos});
  index_t id0 = 0;
  std::vector<trip> trips;
  c->fill_J(q, id0, trips);
  mem->fill_J(q, id0, trips);
  vecX g = vecX::Zero(q.size());
  c->add_gradient(q, g);
  mem->add_gradient(q, g);
  console_logger::info << "weingarten_gn at rest: residual " << c->residual_sq(q) << " memory "
                       << mem->residual_sq(q) << " |g| " << g.norm() << std::endl;
  GAUDI_EXPECT(c->residual_sq(q) < real(1e-16));
  GAUDI_EXPECT(mem->residual_sq(q) < real(1e-16));
  GAUDI_EXPECT(g.norm() < real(1e-8));
}

/// Frames through `projection_solver`: each frame resets the store to the
/// live reading and the memory rest to the frame start (restarts), and a
/// test pin pulls the store to the target's reading.
inline real run_hepworth_frames(const std::vector<vec3> &x0, const std::vector<vec3> &xt,
                                const std::vector<std::array<int, 3>> &tris, int frames,
                                real w_c, real w_m, real w_pin, int backtracks = 8) {
  const kusama::stencil_table t =
      kusama::build_stencil_table(asawa::shell::triangle_walker(static_cast<int>(x0.size()), tris));
  const std::vector<vec3> target = wg::edge_readings(t, xt);
  const int n = static_cast<int>(x0.size());
  const real h = 0.1;
  std::vector<vec3> x = x0, v(n, vec3::Zero()), f = v, m(n, vec3::Ones());
  auto pos = hepworth::vec3_block::create(m, x, v, f);
  std::vector<vec3> hs, hv, hf, hm;
  auto store = hepworth::vec3_block::create(hm, hs, hv, hf);
  std::vector<hepworth::sim_block::ptr> blocks{pos, store};
  hb::projection_solver solver;
  solver.gn_refresh_every = 0;
  solver.gn_backtracks = backtracks;
  for (int k = 0; k < frames; ++k) {
    hs = wg::edge_readings(t, x);
    hv.assign(hs.size(), vec3::Zero());
    hf.assign(hs.size(), vec3::Zero());
    hm.assign(hs.size(), vec3::Constant(real(1e-8)));
    std::vector<hepworth::projection_constraint::ptr> cs{
        hb::weingarten_gn::create(&t, x, w_c, blocks),
        hb::dirichlet_memory_gn::create(&t, x, w_m, std::vector<hepworth::sim_block::ptr>{pos}),
        hb::edge_shape_fair::create(&t, x, hb::edge_fair_rows::pin, w_pin,
                                    std::vector<hepworth::sim_block::ptr>{store},
                                    hb::goal_fixed(target))};
    solver.set_constraints(cs);
    solver.step(blocks, h, real(1), 30);
  }
  return rigid_error(x, xt);
}

GAUDI_TEST(weingarten_gn_torus_tube_over_frames) {
  std::vector<vec3> x0, xt;
  std::vector<std::array<int, 3>> tris;
  torus_tris(2.0, 0.7, 24, 12, x0, tris);
  torus_tris(2.0, 0.9, 24, 12, xt, tris);
  const real e0 = rigid_error(x0, xt);
  const real e1 = run_hepworth_frames(x0, xt, tris, 20, real(1e3), real(1), real(1e4));
  console_logger::info << "weingarten_gn torus tube over frames: rigid err " << e0 << " -> " << e1
                       << std::endl;
  GAUDI_EXPECT(e1 < real(1e-4) * e0);
}

GAUDI_TEST(weingarten_gn_cube_to_sphere_over_frames) {
  std::vector<vec3> x0;
  std::vector<std::array<int, 3>> tris;
  cube_tris(6, x0, tris);
  const real R = std::sqrt(real(24) / (4 * M_PI));
  std::vector<vec3> xt(x0.size());
  for (size_t i = 0; i < x0.size(); ++i)
    xt[i] = R * x0[i].normalized();
  const real e0 = rigid_error(x0, xt);
  // Inertia is the proximal anchor, so low-frequency rounding advances a
  // fraction per frame; a stiffer read-out pushes more through.
  const real e1 = run_hepworth_frames(x0, xt, tris, 20, real(1e4), real(1), real(1e5));
  console_logger::info << "weingarten_gn cube -> sphere over frames: rigid err " << e0 << " -> "
                       << e1 << std::endl;
  GAUDI_EXPECT(e1 < real(2e-3) * e0);
}

} // namespace
} // namespace test
} // namespace gaudi

#endif
