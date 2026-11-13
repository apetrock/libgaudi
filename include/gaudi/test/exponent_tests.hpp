#ifndef GAUDI_EXPONENT_TESTS_HPP
#define GAUDI_EXPONENT_TESTS_HPP

/// The exponent curvature solver (hepworth/block/exponent_coupling.hpp;
/// MATH.md section 16). Included after weingarten_tests.hpp, whose helpers
/// (torus_tris, rigid_error) it uses.

#include "gaudi/hepworth/block/exponent_coupling.hpp"
#include "gaudi/kusama/principal_curvature_solver.hpp"
#include "gaudi/test/test.hpp"
#include "gaudi/test/weingarten_tests.hpp"

#include <chrono>
#include <cmath>
#include <random>
#include <vector>

namespace gaudi {
namespace test {
namespace {

namespace ex = hepworth::block::exponent;

/// The harness strip (cr_ops.strip_param): `ns x nz` grid over
/// `width x height`, two triangles per cell.
inline void strip_grid(int ns, int nz, real width, real height, std::vector<std::array<real, 2>> &P,
                       std::vector<std::array<int, 3>> &tris) {
  P.clear();
  tris.clear();
  for (int i = 0; i < ns; ++i)
    for (int j = 0; j < nz; ++j)
      P.push_back({-width / 2 + width * real(i) / (ns - 1), -height / 2 + height * real(j) / (nz - 1)});
  for (int i = 0; i + 1 < ns; ++i)
    for (int j = 0; j + 1 < nz; ++j) {
      const int a = i * nz + j, b = (i + 1) * nz + j;
      tris.push_back({a, b, b + 1});
      tris.push_back({a, b + 1, a + 1});
    }
}

/// The closed band of radius `r`, `ns` around and `nz` along.
inline void closed_band(int ns, int nz, real r, std::vector<vec3> &x, std::vector<std::array<int, 3>> &tris) {
  x.clear();
  tris.clear();
  for (int i = 0; i < ns; ++i)
    for (int j = 0; j < nz; ++j) {
      const real a = real(2) * real(M_PI) * real(i) / ns;
      x.push_back(vec3(r * std::cos(a), r * std::sin(a), -real(0.5) + real(j) / (nz - 1)));
    }
  for (int i = 0; i < ns; ++i)
    for (int j = 0; j + 1 < nz; ++j) {
      const int a = i * nz + j, b = ((i + 1) % ns) * nz + j;
      tris.push_back({a, b, b + 1});
      tris.push_back({a, b + 1, a + 1});
    }
}

inline kusama::stencil_table table_of(int n, const std::vector<std::array<int, 3>> &tris) {
  return kusama::build_stencil_table(asawa::shell::triangle_walker(n, tris));
}

/// The rebuilt face read-out `s(l, th)` (the stencil from the forms).
inline vec3 forms_readout(const std::array<real, 9> &l, const vec3 &th, const vec3 &h) {
  std::array<vec3, 6> P;
  Eigen::Matrix<real, 18, 12> dP;
  hb::cr_generated::cr_ff_rebuild(l, th, P, dP);
  for (int i = 0; i < 3; ++i)
    if (!(h[i] > 0))
      P[static_cast<size_t>(3 + i)] = P[static_cast<size_t>(i)] + P[static_cast<size_t>((i + 1) % 3)] -
                                      P[static_cast<size_t>((i + 2) % 3)];
  vec3 s;
  Eigen::Matrix<real, 3, 18> J;
  hb::cr_generated::cr_curvature_stencil(P, h, s, J);
  return s;
}

/// Ratio q_min / q_max of the height quadric fitted above the best plane
/// (1 dish, 0 cylinder), and |q|_max.
inline real quadric_ratio(const std::vector<vec3> &x, real &qmax) {
  vec3 c = vec3::Zero();
  for (const vec3 &p : x)
    c += p;
  c /= real(x.size());
  Eigen::MatrixXd D(x.size(), 3);
  for (size_t i = 0; i < x.size(); ++i)
    D.row(static_cast<Eigen::Index>(i)) = (x[i] - c).transpose();
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(D, Eigen::ComputeThinV);
  const vec3 e1 = svd.matrixV().col(0), e2 = svd.matrixV().col(1), n = svd.matrixV().col(2);
  Eigen::MatrixXd A(x.size(), 6);
  Eigen::VectorXd b(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    const vec3 d = x[i] - c;
    const real u = d.dot(e1), v = d.dot(e2);
    A.row(static_cast<Eigen::Index>(i)) << real(0.5) * u * u, u * v, real(0.5) * v * v, u, v, 1;
    b[static_cast<Eigen::Index>(i)] = d.dot(n);
  }
  const Eigen::VectorXd q = A.colPivHouseholderQr().solve(b);
  Eigen::Matrix<real, 2, 2> Q;
  Q << q[0], q[1], q[1], q[2];
  const Eigen::Matrix<real, 2, 1> ev = Eigen::SelfAdjointEigenSolver<Eigen::Matrix<real, 2, 2>>(Q).eigenvalues();
  const real big = std::abs(ev[0]) > std::abs(ev[1]) ? ev[0] : ev[1];
  const real small = std::abs(ev[0]) > std::abs(ev[1]) ? ev[1] : ev[0];
  qmax = std::abs(big);
  return small / big;
}

/// `J_ff` (rebuild x generated read-out) against central differences of the
/// rebuilt read-out, and the rebuilt read-out against the live one, on a
/// perturbed torus and an open strip.
GAUDI_TEST(exponent_forms_jacobian) {
  std::mt19937 rng(7);
  std::normal_distribution<real> nd(0, 1);
  for (int mesh = 0; mesh < 2; ++mesh) {
    std::vector<vec3> x;
    std::vector<std::array<int, 3>> tris;
    if (mesh == 0) {
      torus_tris(2.0, 0.7, 12, 8, x, tris);
    } else {
      std::vector<std::array<real, 2>> P;
      strip_grid(6, 4, 2.0, 1.0, P, tris);
      for (const auto &p : P)
        x.push_back(vec3(p[0], real(0.1) * std::sin(real(3) * p[0]), p[1]));
    }
    for (vec3 &p : x)
      p += real(0.03) * vec3(nd(rng), nd(rng), nd(rng));
    const kusama::stencil_table t = table_of(static_cast<int>(x.size()), tris);
    real e_val = 0, e_J = 0;
    for (int f = 0; f < t.n_faces; ++f) {
      vec3 h, s;
      const std::array<vec3, 6> P = kusama::stencil_positions(t, x, f, h);
      const Eigen::Matrix<real, 3, 12> J = ex::face_forms_jacobian(P, h, s);
      vec3 s_live;
      Eigen::Matrix<real, 3, 18> Jl;
      hb::cr_generated::cr_curvature_stencil(P, h, s_live, Jl);
      e_val = std::max(e_val, (s - s_live).norm() / (1 + s_live.norm()));
      std::array<real, 9> l;
      vec3 th;
      ex::face_forms(P, h, l, th);
      const real eps = 1e-6;
      for (int j = 0; j < 12; ++j) {
        if (j >= 3 && j < 9 && !(h[(j - 3) / 2] > 0))
          continue;
        if (j >= 9 && !(h[j - 9] > 0))
          continue;
        std::array<real, 9> lp = l, lm = l;
        vec3 tp = th, tm = th;
        if (j < 9) {
          lp[static_cast<size_t>(j)] *= std::exp(eps);
          lm[static_cast<size_t>(j)] *= std::exp(-eps);
        } else {
          tp[j - 9] += eps;
          tm[j - 9] -= eps;
        }
        const vec3 fd = (forms_readout(lp, tp, h) - forms_readout(lm, tm, h)) / (2 * eps);
        e_J = std::max(e_J, (J.col(j) - fd).norm() / (1 + fd.norm()));
      }
    }
    console_logger::info << "exponent forms: mesh " << mesh << " rebuilt vs live " << e_val
                         << ", J_ff vs FD " << e_J << std::endl;
    GAUDI_EXPECT(e_val < real(1e-11));
    GAUDI_EXPECT(e_J < real(1e-6));
  }
}

/// `R(omega) U(L)` against a series matrix exponential of `[omega]x + P L P^T`:
/// exact for pure rotation and pure stretch, second order in the mix.
GAUDI_TEST(exponent_face_deformation) {
  auto expm = [](const Eigen::Matrix<real, 3, 3> &A) {
    int k = 0;
    real n = A.norm();
    while (n > 0.5) {
      n *= 0.5;
      ++k;
    }
    const Eigen::Matrix<real, 3, 3> S = A / std::pow(real(2), k);
    Eigen::Matrix<real, 3, 3> E = Eigen::Matrix<real, 3, 3>::Identity(), term = E;
    for (int i = 1; i < 20; ++i) {
      term = term * S / real(i);
      E += term;
    }
    for (int i = 0; i < k; ++i)
      E = E * E;
    return E;
  };
  ex::mat32 P;
  P.col(0) = vec3(1, 0, 0);
  P.col(1) = vec3(0, 1, 0);
  const vec3 N(0, 0, 1);
  real e_rot = 0, e_str = 0, e_mix = 0;
  for (real s : {real(1e-3), real(0.3), real(1.2)}) {
    const vec3 w = s * vec3(0.3, -0.5, 0.7), L = s * vec3(0.2, -0.1, 0.15);
    e_rot = std::max(e_rot, (ex::face_deformation(w, vec3::Zero(), P, N) - expm(ex::skew(w))).norm());
    e_str = std::max(e_str, (ex::face_deformation(vec3::Zero(), L, P, N) -
                             expm(P * (Eigen::Matrix<real, 2, 2>() << L[0], L[1], L[1], L[2]).finished() * P.transpose()))
                                .norm());
    if (s < real(0.01))
      e_mix = (ex::face_deformation(w, L, P, N) -
               expm(ex::skew(w) + P * (Eigen::Matrix<real, 2, 2>() << L[0], L[1], L[1], L[2]).finished() * P.transpose()))
                  .norm() /
              (w.norm() * L.norm());
  }
  console_logger::info << "exponent deformation: rotation " << e_rot << ", stretch " << e_str
                       << ", mix / (|w| |L|) at 1e-3 " << e_mix << std::endl;
  GAUDI_EXPECT(e_rot < real(1e-12));
  GAUDI_EXPECT(e_str < real(1e-12));
  GAUDI_EXPECT(e_mix < real(1));
}

/// Targets equal to the live read-out: no step moves anything.
GAUDI_TEST(exponent_idle) {
  std::vector<vec3> x;
  std::vector<std::array<int, 3>> tris;
  torus_tris(2.0, 0.7, 16, 8, x, tris);
  const kusama::stencil_table t = table_of(static_cast<int>(x.size()), tris);
  const std::vector<vec3> goal = ex::face_readouts(t, x);
  std::vector<vec3> y = x;
  const ex::solve_report rep = ex::solve(t, y, goal, ex::config{}, 10);
  real move = 0;
  for (size_t i = 0; i < x.size(); ++i)
    move = std::max(move, (y[i] - x[i]).norm());
  console_logger::info << "exponent idle: merit " << rep.merit.front() << ", iterations " << rep.iterations
                       << ", max move " << move << std::endl;
  GAUDI_EXPECT(rep.merit.front() < real(1e-10));
  GAUDI_EXPECT(move < real(1e-9));
}

/// Realizable targets (read-outs of meshes with the same connectivity) are
/// reached: strip flat to cylinder (bending), torus tube 0.7 to 0.9 and a
/// closed band radius 1 to 1.3 (both need stretch). Python parity
/// (exponent_solver.py): shape 1e-10 to 1e-9 in 4 to 6 iterates.
GAUDI_TEST(exponent_parity) {
  struct row {
    const char *name;
    std::vector<vec3> x0, xt;
    std::vector<std::array<int, 3>> tris;
  };
  std::vector<row> rows(3);
  {
    std::vector<std::array<real, 2>> P;
    strip_grid(17, 6, 2.0, 1.0, P, rows[0].tris);
    rows[0].name = "strip flat -> cylinder 1.5";
    for (const auto &p : P) {
      rows[0].x0.push_back(vec3(p[0], 0, p[1]));
      const real a = p[0] / real(1.5) - real(0.5) * real(M_PI);
      rows[0].xt.push_back(vec3(real(1.5) * std::cos(a), real(1.5) * std::sin(a) + real(1.5), p[1]));
    }
  }
  rows[1].name = "torus tube 0.7 -> 0.9";
  torus_tris(2.0, 0.7, 24, 12, rows[1].x0, rows[1].tris);
  {
    std::vector<std::array<int, 3>> tt;
    torus_tris(2.0, 0.9, 24, 12, rows[1].xt, tt);
  }
  rows[2].name = "closed band 1 -> 1.3";
  closed_band(24, 6, 1.0, rows[2].x0, rows[2].tris);
  {
    std::vector<std::array<int, 3>> tt;
    closed_band(24, 6, 1.3, rows[2].xt, tt);
  }
  for (auto &rw : rows) {
    const kusama::stencil_table t = table_of(static_cast<int>(rw.x0.size()), rw.tris);
    const std::vector<vec3> goal = ex::face_readouts(t, rw.xt);
    std::vector<vec3> x = rw.x0;
    const ex::solve_report rep = ex::solve(t, x, goal, ex::config{}, 60);
    real size = 0;
    for (const vec3 &p : rw.xt)
      for (const vec3 &q : rw.xt)
        size = std::max(size, (p - q).cwiseAbs().maxCoeff());
    const real err = rigid_error(x, rw.xt) / size;
    console_logger::info << "exponent parity " << rw.name << ": iterations " << rep.iterations << ", merit "
                         << rep.merit.front() << " -> " << rep.merit.back() << ", rigid error / size " << err
                         << std::endl;
    GAUDI_EXPECT(rep.merit.back() < real(1e-10) * rep.merit.front());
    // An open patch's read-out does not pin its shape (MATH.md section 11);
    // closed surfaces land on the target.
    GAUDI_EXPECT(err < (&rw == &rows[0] ? real(1e-2) : real(1e-10)));
  }
}

/// Frustrated isotropic curvature on a flat sheet (a sphere cap's read-out):
/// soft growth stretches into a dish, an elastic stiff metric snaps to a
/// cylinder (the bimetal plate).
GAUDI_TEST(exponent_plate_buckles) {
  std::vector<std::array<real, 2>> P;
  std::vector<std::array<int, 3>> tris;
  strip_grid(13, 13, 2.0, 2.0, P, tris);
  std::vector<vec3> flat, cap;
  for (const auto &p : P) {
    flat.push_back(vec3(p[0], 0, p[1]));
    const real rr = std::hypot(p[0], p[1]);
    const real a = rr;  // cap radius 1
    const real dx = rr > 1e-12 ? p[0] / rr : 0, dz = rr > 1e-12 ? p[1] / rr : 0;
    cap.push_back(vec3(std::sin(a) * dx, -(1 - std::cos(a)), std::sin(a) * dz));
  }
  const kusama::stencil_table t = table_of(static_cast<int>(flat.size()), tris);
  const std::vector<vec3> goal = ex::face_readouts(t, cap);
  std::vector<std::array<vec3, 3>> rest(static_cast<size_t>(t.n_faces));
  for (int f = 0; f < t.n_faces; ++f)
    for (int k = 0; k < 3; ++k)
      rest[static_cast<size_t>(f)][static_cast<size_t>(k)] = flat[static_cast<size_t>(t.ids[static_cast<size_t>(f)][static_cast<size_t>(k)])];
  std::mt19937 rng(3);
  std::normal_distribution<real> nd(0, 1);
  std::vector<vec3> seed = flat;
  for (vec3 &p : seed)
    p[1] += real(1e-3) * nd(rng);
  real ratios[2], qs[2];
  const real kms[2] = {1.0, 1e4};
  for (int k = 0; k < 2; ++k) {
    ex::config cfg;
    cfg.km = kms[k];
    std::vector<vec3> x = seed;
    ex::solve(t, x, goal, cfg, 90, 30, real(1e-13), &rest);
    ratios[k] = quadric_ratio(x, qs[k]);
    console_logger::info << "exponent plate km " << kms[k] << ": ratio " << ratios[k] << ", |q| " << qs[k]
                         << std::endl;
  }
  GAUDI_EXPECT(ratios[0] > real(0.9));
  GAUDI_EXPECT(std::abs(ratios[1]) < real(0.1));
}

/// A factor reused on a rigidly turned torus steps as a fresh one there (the
/// body-frame map), and the increment-form gather is idle at a zero step
/// on positions the system was not built at.
GAUDI_TEST(exponent_frozen_factor) {
  std::vector<vec3> x;
  std::vector<std::array<int, 3>> tris;
  torus_tris(1.0, 0.33, 24, 12, x, tris);
  const kusama::stencil_table t = table_of(static_cast<int>(x.size()), tris);
  const ex::topology tp = ex::build_topology(t);
  const ex::config cfg;
  const auto F = std::make_shared<const ex::factor>(t, tp, x, cfg);
  const Eigen::Matrix<real, 3, 3> R =
      Eigen::AngleAxis<real>(real(0.7), vec3(1, 2, 3).normalized()).toRotationMatrix();
  std::vector<vec3> y(x.size());
  for (size_t i = 0; i < x.size(); ++i)
    y[i] = R * x[i] + vec3(0.3, -0.2, 0.1);
  const ex::system reused(t, tp, y, cfg, F), fresh(t, tp, y, cfg);
  GAUDI_EXPECT(!reused.fresh() && fresh.fresh());
  std::vector<vec3> goal = ex::face_readouts(t, y);
  for (vec3 &g : goal)
    g *= real(1.05);
  vecX r;
  ex::merit(t, y, goal, reused.weights(), r);
  // Compared as gathered shapes: the rigid part of a step is held only by
  // the tiny `kb`, so it carries the round-off.
  const vecX a = reused.solve(r), b = fresh.solve(r);
  const std::vector<vec3> ya = reused.gather(reused.frames(y, a, 1), y), yb = fresh.gather(fresh.frames(y, b, 1), y);
  real moved = 0;
  for (size_t i = 0; i < y.size(); ++i)
    moved = std::max(moved, (yb[i] - y[i]).cwiseAbs().maxCoeff());
  const real step = rigid_error(ya, yb) / moved;
  // The floor: two fresh factors, at x and at the turned y.
  const ex::system at_x(t, tp, x, cfg);
  const std::vector<vec3> yx = at_x.gather(at_x.frames(x, at_x.solve(r), 1), x);
  const real floor = rigid_error(yx, yb) / moved;

  std::vector<vec3> z = y;
  for (size_t i = 0; i < z.size(); ++i)
    z[i] += real(0.01) * vec3(std::sin(real(i)), std::cos(real(3 * i)), std::sin(real(7 * i)));
  const std::vector<vec3> g = reused.gather(reused.frames(z), z);
  real idle = 0;
  for (size_t i = 0; i < z.size(); ++i)
    idle = std::max(idle, (g[i] - z[i]).cwiseAbs().maxCoeff());
  console_logger::info << "exponent frozen factor: reused vs fresh shape / move " << step
                       << " (fresh at x vs fresh at y " << floor << ", raw step difference "
                       << (a - b).cwiseAbs().maxCoeff() / b.cwiseAbs().maxCoeff() << "), idle gather on moved x "
                       << idle << std::endl;
  GAUDI_EXPECT(step < real(1e-5) && step < real(2) * floor + real(1e-12));
  GAUDI_EXPECT(idle < real(1e-12));
}

/// `d (k1, k2) / d s` against central differences of the eigenvalues away
/// from umbilics; at an umbilic both rows are `dH` (Eigensatz A.6).
GAUDI_TEST(exponent_principal_rows) {
  std::mt19937 rng(11);
  std::normal_distribution<real> nd(0, 1);
  real e = 0;
  for (int k = 0; k < 50; ++k) {
    const vec3 s(nd(rng), nd(rng), nd(rng));
    const Eigen::Matrix<real, 3, 3> K = ex::principal_rows(s, real(1e-3));
    const real h = 1e-6;
    for (int j = 0; j < 3; ++j) {
      vec3 sp = s, sm = s;
      sp[j] += h;
      sm[j] -= h;
      const vec3 fd = (ex::principal_values(sp) - ex::principal_values(sm)) / (2 * h);
      e = std::max(e, (K.col(j) - fd).norm());
    }
  }
  const Eigen::Matrix<real, 3, 3> U = ex::principal_rows(vec3(0.7, 0, 0.7), real(1e-3));
  const real eu = std::max((U.row(0) - Eigen::Matrix<real, 1, 3>(0.5, 0, 0.5)).norm(),
                           (U.row(1) - Eigen::Matrix<real, 1, 3>(0.5, 0, 0.5)).norm());
  console_logger::info << "exponent principal rows: vs FD " << e << ", umbilic vs dH " << eu << std::endl;
  GAUDI_EXPECT(e < real(1e-8));
  GAUDI_EXPECT(eu < real(1e-15));
}

/// Principal-curvature targets: idle at the live ones; the read-out
/// eigenvalues of the 0.9 tube are reached from the 0.7 one (directions
/// free, so the torus lands on the target up to its symmetry).
GAUDI_TEST(exponent_principal_parity) {
  std::vector<vec3> x0, xt;
  std::vector<std::array<int, 3>> tris, tt;
  torus_tris(2.0, 0.7, 24, 12, x0, tris);
  torus_tris(2.0, 0.9, 24, 12, xt, tt);
  const kusama::stencil_table t = table_of(static_cast<int>(x0.size()), tris);
  ex::config cfg;
  cfg.kind = ex::residual::principal;
  auto eigs = [&](const std::vector<vec3> &x) {
    std::vector<vec3> k = ex::face_readouts(t, x);
    for (vec3 &v : k)
      v = ex::principal_values(v);
    return k;
  };
  std::vector<vec3> y = x0;
  const ex::solve_report idle = ex::solve(t, y, eigs(x0), cfg, 10);
  real move = 0;
  for (size_t i = 0; i < y.size(); ++i)
    move = std::max(move, (y[i] - x0[i]).norm());
  std::vector<vec3> x = x0;
  const ex::solve_report rep = ex::solve(t, x, eigs(xt), cfg, 60);
  real size = 0;
  for (const vec3 &p : xt)
    size = std::max(size, p.cwiseAbs().maxCoeff());
  const real err = rigid_error(x, xt) / (2 * size);
  console_logger::info << "exponent principal: idle merit " << idle.merit.front() << " move " << move
                       << "; tube 0.7 -> 0.9 iterations " << rep.iterations << ", merit " << rep.merit.front()
                       << " -> " << rep.merit.back() << ", rigid error / size " << err << std::endl;
  GAUDI_EXPECT(idle.merit.front() < real(1e-10) && move < real(1e-9));
  GAUDI_EXPECT(rep.merit.back() < real(1e-8) * rep.merit.front());
  GAUDI_EXPECT(err < real(1e-6));
}

/// Eigensatz's metric terms hold against the rest: the 0.7 -> 0.9 tube's
/// targets need anisotropic and areal stretch; `w_conformal` cuts the
/// trace-free part of the reached stretch, `w_area` its trace.
GAUDI_TEST(exponent_principal_metric_rows) {
  std::vector<vec3> x0, xt;
  std::vector<std::array<int, 3>> tris, tt;
  torus_tris(2.0, 0.7, 24, 12, x0, tris);
  torus_tris(2.0, 0.9, 24, 12, xt, tt);
  const kusama::stencil_table t = table_of(static_cast<int>(x0.size()), tris);
  std::vector<vec3> goal = ex::face_readouts(t, xt);
  for (vec3 &v : goal)
    v = ex::principal_values(v);
  std::vector<std::array<vec3, 3>> rest(static_cast<size_t>(t.n_faces));
  for (int f = 0; f < t.n_faces; ++f)
    for (int k = 0; k < 3; ++k)
      rest[static_cast<size_t>(f)][static_cast<size_t>(k)] =
          x0[static_cast<size_t>(t.ids[static_cast<size_t>(f)][static_cast<size_t>(k)])];
  auto stretch = [&](real wc, real wa, real &dev, real &tr) {
    ex::config cfg;
    cfg.kind = ex::residual::principal;
    cfg.w_conformal = wc;
    cfg.w_area = wa;
    std::vector<vec3> x = x0;
    ex::solve(t, x, goal, cfg, 60, 30, real(1e-13), &rest);
    dev = tr = 0;
    for (const vec3 &L : ex::current_log_stretch(t, x, rest)) {
      dev += std::hypot(real(0.5) * (L[0] - L[2]), L[1]);
      tr += std::abs(L[0] + L[2]);
    }
    dev /= real(t.n_faces);
    tr /= real(t.n_faces);
  };
  real d0, t0, dc, tc, da, ta;
  stretch(0, 0, d0, t0);
  stretch(10, 0, dc, tc);
  stretch(0, 10, da, ta);
  console_logger::info << "exponent principal metric: free |dev L| " << d0 << " tr " << t0
                       << "; conformal 10 |dev L| " << dc << "; area 10 |tr L| " << ta << std::endl;
  GAUDI_EXPECT(dc < real(0.5) * d0);
  GAUDI_EXPECT(ta < real(0.5) * t0);
}

/// The kusama Eigensatz solver's torus x1.2 case (vertex targets, its own
/// discretization) against principal rows on the exponent system (face
/// targets at the centroids): both land on R 1.2, tube 0.396.
GAUDI_TEST(exponent_principal_vs_kusama) {
  const real R = 1.0, r = 0.33, sc = 1.2;
  std::vector<vec3> x0;
  std::vector<std::array<int, 3>> tris;
  torus_tris(R, r, 48, 24, x0, tris);
  auto fit = [](const std::vector<vec3> &x) {
    real lo = 1e30, hi = 0;
    for (const vec3 &p : x) {
      lo = std::min(lo, std::hypot(p[0], p[1]));
      hi = std::max(hi, std::hypot(p[0], p[1]));
    }
    const real Rc = real(0.5) * (lo + hi);
    real tube = 0;
    for (const vec3 &p : x)
      tube += std::hypot(std::hypot(p[0], p[1]) - Rc, p[2]);
    return vec2(Rc, tube / real(x.size()));
  };
  auto kappa = [&](const vec3 &p) {
    const real phi = std::atan2(p[2], std::hypot(p[0], p[1]) - R);
    return vec2(real(1) / (sc * r), std::cos(phi) / (sc * (R + r * std::cos(phi))));
  };
  // kusama: Levenberg-Marquardt on vertex positions.
  std::vector<vec3> xk = x0;
  {
    kusama::principal_curvature_solver s;
    std::vector<kusama::tri3> kt(tris.begin(), tris.end());
    s.reset(static_cast<int>(xk.size()), kt, xk);
    std::vector<vec2> kh(xk.size());
    for (size_t i = 0; i < xk.size(); ++i)
      kh[i] = kappa(x0[i]);
    for (int it = 0, stalls = 0; it < 400 && stalls < 4; ++it) {
      const kusama::principal_curvature_report rep = s.iterate(xk, kh);
      if (!rep.accepted) {
        ++stalls;
        continue;
      }
      stalls = 0;
      if (rep.energy0 - rep.energy1 <= 1e-12 * rep.energy0)
        break;
    }
  }
  // exponent: the same targets at face centroids, in the read-out's sign.
  const kusama::stencil_table t = table_of(static_cast<int>(x0.size()), tris);
  const std::vector<vec3> s0 = ex::face_readouts(t, x0);
  std::vector<vec3> goal(static_cast<size_t>(t.n_faces));
  real corr = 0;
  for (int f = 0; f < t.n_faces; ++f) {
    const auto &id = t.ids[static_cast<size_t>(f)];
    const vec3 c = (x0[static_cast<size_t>(id[0])] + x0[static_cast<size_t>(id[1])] + x0[static_cast<size_t>(id[2])]) / 3;
    const vec2 k = kappa(c);
    goal[static_cast<size_t>(f)] = vec3(k[0], k[1], 0);
    corr += (s0[static_cast<size_t>(f)][0] + s0[static_cast<size_t>(f)][2]) * (k[0] + k[1]);
  }
  if (corr < 0)
    for (vec3 &g : goal)
      g = vec3(-g[1], -g[0], 0);
  ex::config cfg;
  cfg.kind = ex::residual::principal;
  std::vector<vec3> xe = x0;
  const ex::solve_report rep = ex::solve(t, xe, goal, cfg, 60);
  const vec2 fk = fit(xk), fe = fit(xe);
  console_logger::info << "exponent principal vs kusama (torus x1.2): kusama R " << fk[0] << " tube " << fk[1]
                       << "; exponent R " << fe[0] << " tube " << fe[1] << " (iterations " << rep.iterations
                       << ", merit " << rep.merit.front() << " -> " << rep.merit.back() << ")" << std::endl;
  for (const vec2 &f : {fk, fe}) {
    GAUDI_EXPECT(std::abs(f[0] - sc * R) / (sc * R) < real(0.02));
    GAUDI_EXPECT(std::abs(f[1] - sc * r) / (sc * r) < real(0.02));
  }
}

/// At equal weights the exponent coupling pulls on positions like `bend` and
/// the `weingarten` coupling (`weingarten_coupling_matches_bend_at_equal_weight`):
/// the mean `diag(A^T A)` over position columns at `w = 1`.
GAUDI_TEST(exponent_coupling_matches_bend_at_equal_weight) {
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
    auto coupling = hb::exponent_coupling::create(&t, x, hs, real(1),
                                                  std::vector<hepworth::sim_block::ptr>{pos, store}, ex::config{});
    std::vector<trip> tc;
    id0 = 0;
    coupling->fill_A(id0, tc);
    const real kb = mean_pos_diag(tb), kc = mean_pos_diag(tc);
    console_logger::info << "exponent w = 1 position stiffness on torus (" << R << ", " << r << ", " << nu << " x "
                         << nv << "): bend " << kb << " coupling " << kc << " ratio " << kc / kb << std::endl;
    GAUDI_EXPECT(kc / kb > real(1) / 3 && kc / kb < real(3));
  }
}

/// Cost of one system build and one step on the demo torus (48 x 24, 2304
/// faces); reported, not asserted.
GAUDI_TEST(exponent_cost) {
  std::vector<vec3> x;
  std::vector<std::array<int, 3>> tris;
  torus_tris(1.0, 0.33, 48, 24, x, tris);
  const kusama::stencil_table t = table_of(static_cast<int>(x.size()), tris);
  const ex::topology tp = ex::build_topology(t);
  const ex::config cfg;
  auto t0 = std::chrono::steady_clock::now();
  const ex::system sys(t, tp, x, cfg);
  auto t1 = std::chrono::steady_clock::now();
  std::vector<vec3> goal = ex::face_readouts(t, x);
  for (vec3 &g : goal)
    g *= real(1.05);
  vecX r;
  const real m0 = ex::merit(t, x, goal, sys.weights(), r);
  auto t2 = std::chrono::steady_clock::now();
  const vecX xi = sys.solve(r);
  auto t3 = std::chrono::steady_clock::now();
  const std::vector<vec3> xn = sys.gather(sys.frames(x, xi, 1), x);
  auto t4 = std::chrono::steady_clock::now();
  auto ms = [](auto a, auto b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
  console_logger::info << "exponent cost (" << t.n_faces << " faces): system " << ms(t0, t1) << " ms, merit "
                       << ms(t1, t2) << " ms, solve " << ms(t2, t3) << " ms, frames + gather " << ms(t3, t4)
                       << " ms (merit " << m0 << ")" << std::endl;
  GAUDI_EXPECT(sys.ok());
}

} // namespace
} // namespace test
} // namespace gaudi

#endif
