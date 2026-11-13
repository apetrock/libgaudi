#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

#include "gaudi/duchamp/monge_fairing_demo.hpp"

// Steps the monge demo on its built-in torus with startup noise and prints how
// round the tube stays, with the osculating-torus ring row off and on.
//
//   osculating_torus_probe [frames] [perturb]

namespace {
using real = gaudi::real;

struct tube_stats {
  real mean = 0;
  real rel_std = 0;
  real max_abs = 0;
  int verts = 0;
};

tube_stats measure(gaudi::asawa::shell::shell &M) {
  const std::vector<gaudi::vec3> &x0 = gaudi::asawa::get_vec_data(M, 0);
  tube_stats s;
  // Rigid frame from the cloud: centroid, and the axis of least spread.
  std::vector<gaudi::vec3> x;
  gaudi::vec3 c = gaudi::vec3::Zero();
  for (const gaudi::vec3 &p : x0)
    if (p.allFinite()) {
      x.push_back(p);
      c += p;
    }
  if (x.empty())
    return s;
  c /= real(x.size());
  gaudi::mat3 C = gaudi::mat3::Zero();
  for (const gaudi::vec3 &p : x)
    C += (p - c) * (p - c).transpose();
  Eigen::SelfAdjointEigenSolver<gaudi::mat3> es(C);
  const gaudi::vec3 axis = es.eigenvectors().col(0);
  auto split = [&](const gaudi::vec3 &p, real &radial, real &height) {
    const gaudi::vec3 d = p - c;
    height = d.dot(axis);
    radial = (d - height * axis).norm();
  };
  real R = 0;
  for (const gaudi::vec3 &p : x) {
    real rad, hgt;
    split(p, rad, hgt);
    R += rad;
  }
  R /= real(x.size());
  std::vector<real> rho;
  for (const gaudi::vec3 &p : x) {
    real rad, hgt;
    split(p, rad, hgt);
    rho.push_back(std::hypot(rad - R, hgt));
    s.max_abs = std::max(s.max_abs, (p - c).norm());
  }
  for (const real r : rho)
    s.mean += r;
  s.mean /= rho.size();
  real var = 0;
  for (const real r : rho)
    var += (r - s.mean) * (r - s.mean);
  s.rel_std = std::sqrt(var / rho.size()) / s.mean;
  s.verts = static_cast<int>(x0.size());
  return s;
}

void run(const char *label, gaudi::duchamp::monge_fairing_demo_config cfg, int frames) {
  auto demo = gaudi::duchamp::monge_fairing_demo::create(cfg);
  const tube_stats s0 = measure(*demo->__M);
  std::cout << label << "  start: tube r " << s0.mean << " roughness " << s0.rel_std << " verts " << s0.verts
            << "\n";
  for (int f = 0; f < frames; ++f) {
    const auto t0 = std::chrono::steady_clock::now();
    gaudi::hepworth::block::osculating_torus_ring_stats::reset();
    gaudi::hepworth::block::osculating_dupin_ring_stats::reset();
    demo->step(f, !std::getenv("PROBE_NOREMESH"));
    const real dt = std::chrono::duration<real>(std::chrono::steady_clock::now() - t0).count();
    const tube_stats s = measure(*demo->__M);
    using st = gaudi::hepworth::block::osculating_torus_ring_stats;
    if (st::calls > 0)
      std::cout << label << "      torus rows: calls " << st::calls << " no_fit " << st::no_fit << " landed "
                << st::landed << "/" << st::points << " pull/reach "
                << st::move / std::max(st::reach.load(), 1e-30) << " mean k_hi "
                << st::k_hi / std::max<long>(st::shapes, 1) << " mean k_lo "
                << st::k_lo / std::max<long>(st::shapes, 1) << " mesh/store height slope "
                << st::hm_hs / std::max(st::hs_hs.load(), 1e-30) << "\n";
    st::reset();
    using dst = gaudi::hepworth::block::osculating_dupin_ring_stats;
    if (dst::calls > 0)
      std::cout << label << "      dupin rows: calls " << dst::calls << " no_frame " << dst::no_frame
                << " no_fit " << dst::no_fit << " landed " << dst::landed << "/" << dst::points
                << " pull/reach " << dst::move / std::max(dst::reach.load(), 1e-30) << "\n";
    dst::reset();
    const auto &r = demo->last_solve();
    std::cout << label << "  [" << f << "] tube r " << s.mean << " roughness " << s.rel_std << " max|x| "
              << s.max_abs << " verts " << s.verts << " E " << r.E_first << " -> " << r.E_last << " s " << dt
              << "\n";
  }
}

} // namespace

int main(int argc, char **argv) {
  gaudi::hepworth::block::osculating_torus_ring_stats::enabled = true;
  gaudi::hepworth::block::osculating_dupin_ring_stats::enabled = true;
  const bool relative = std::getenv("PROBE_RELATIVE") != nullptr;
  const bool closest = std::getenv("PROBE_CLOSEST") != nullptr;
  const int frames = argc > 1 ? std::atoi(argv[1]) : 8;
  const gaudi::real perturb = argc > 2 ? std::atof(argv[2]) : 0.3;
  gaudi::duchamp::monge_fairing_demo_config base{};
  base.perturb = perturb;

  if (std::getenv("PROBE_ALL"))
    run("frame_fair", base, frames);

  auto both = base;
  both.project.w_osculating_torus = 1.0;
  if (std::getenv("PROBE_ALL"))
    run("frame_fair+torus", both, frames);

  if (std::getenv("PROBE_ITS"))
    base.step.iterations = std::atoi(std::getenv("PROBE_ITS"));
  if (std::getenv("PROBE_GAIN"))
    base.project.osculating_torus.store_gain = std::atof(std::getenv("PROBE_GAIN"));
  if (std::getenv("PROBE_TWO"))
    base.project.osculating_torus.two_ring = std::atof(std::getenv("PROBE_TWO"));
  if (std::getenv("PROBE_POSE_ITS"))
    base.project.osculating_torus.iterations = std::atoi(std::getenv("PROBE_POSE_ITS"));
  if (std::getenv("PROBE_PRIOR"))
    base.project.osculating_torus.prior = std::atof(std::getenv("PROBE_PRIOR"));
  if (std::getenv("PROBE_NO_GROWTH"))
    base.shape.w_growth = 0.0;
  if (std::getenv("PROBE_ALONE")) {
    base.shape.w_growth = 0.0;
    base.shape.w_det = 0.0;
    base.shape.w_anisotropy = 0.0;
    base.position.w_bend = 0.0;
    base.position.w_conformal = 0.0;
  }
  base.project.osculating_torus.relative = relative;
  base.project.osculating_torus.closest_point = closest;
  if (std::getenv("PROBE_WEINGARTEN"))
    base.project.curvature.kind = gaudi::duchamp::curvature_solver::weingarten;
  if (std::getenv("PROBE_DUPIN")) {
    base.project.w_osculating_torus = 0.0;
    base.project.w_frame_fair = 0.0;
    base.project.w_osculating_dupin = std::atof(std::getenv("PROBE_DUPIN"));
    base.project.w_osculating_dupin_shape =
        std::getenv("PROBE_DUPIN_SHAPE") ? std::atof(std::getenv("PROBE_DUPIN_SHAPE")) : 0.0;
  }
  if (std::getenv("PROBE_DUPIN_POSE_ITS"))
    base.project.osculating_dupin.pose_iterations = std::atoi(std::getenv("PROBE_DUPIN_POSE_ITS"));
  if (std::getenv("PROBE_DUPIN_TWO"))
    base.project.osculating_dupin.two_ring = std::atof(std::getenv("PROBE_DUPIN_TWO"));
  if (std::getenv("PROBE_DUPIN_GAIN"))
    base.project.osculating_dupin.store_gain = std::atof(std::getenv("PROBE_DUPIN_GAIN"));
  if (std::getenv("PROBE_DUPIN_JOINT"))
    base.project.osculating_dupin.direct_shape = false;
  if (std::getenv("PROBE_DUPIN_VW"))
    base.project.osculating_dupin.volume.weight = std::atof(std::getenv("PROBE_DUPIN_VW"));
  if (std::getenv("PROBE_DUPIN_NOB"))
    base.project.osculating_dupin.third_order = false;
  if (std::getenv("PROBE_DUPIN_BCAP"))
    base.project.osculating_dupin.b_cap = std::atof(std::getenv("PROBE_DUPIN_BCAP"));
  if (std::getenv("PROBE_DUPIN_RIDGE"))
    base.project.osculating_dupin.b_ridge = std::atof(std::getenv("PROBE_DUPIN_RIDGE"));
  if (std::getenv("PROBE_DUPIN_FD"))
    base.project.osculating_dupin.exact_gradient = false;
  if (const char *v = std::getenv("PROBE_VOLUME")) {
    using gaudi::hepworth::block::ring_volume_mode;
    const std::string name(v);
    const ring_volume_mode mode = name == "hard"    ? ring_volume_mode::hard
                                  : name == "fit"   ? ring_volume_mode::fit
                                  : name == "corot" ? ring_volume_mode::corotational
                                                    : ring_volume_mode::soft;
    const bool land = std::getenv("PROBE_VOLUME_ONLY") == nullptr;
    const gaudi::real weight = std::getenv("PROBE_VOLUME_W") ? std::atof(std::getenv("PROBE_VOLUME_W")) : 1.0;
    const gaudi::real store = std::getenv("PROBE_STORE_W") ? std::atof(std::getenv("PROBE_STORE_W")) : 10.0;
    for (auto *vol : {&base.project.osculating_torus.volume, &base.project.osculating_dupin.volume}) {
      *vol = {};
      vol->mode = mode;
      vol->weight = weight;
      vol->land = land;
      vol->store_weight = store;
      vol->held_medial = std::getenv("PROBE_FOLLOW") == nullptr;
      vol->corotational.w_volume = weight;
      if (std::getenv("PROBE_PERSIST"))
        vol->persistent = std::atoi(std::getenv("PROBE_PERSIST")) != 0;
      if (std::getenv("PROBE_FOOT_K"))
        vol->corotational.foot_stiffness = std::atof(std::getenv("PROBE_FOOT_K"));
      if (std::getenv("PROBE_REST_W"))
        vol->corotational.w_rest = std::atof(std::getenv("PROBE_REST_W"));
      if (std::getenv("PROBE_FOOT_W"))
        vol->foot_weight = std::atof(std::getenv("PROBE_FOOT_W"));
      if (std::getenv("PROBE_LO_STIFF"))
        vol->lo_stiffness = std::atof(std::getenv("PROBE_LO_STIFF"));
    }
  }
  auto torus = base;
  torus.project.w_frame_fair =
      std::getenv("PROBE_FRAME_FAIR") ? std::atof(std::getenv("PROBE_FRAME_FAIR")) : 0.0;
  // A cyclide run keeps the torus row off unless PROBE_W asks for it too.
  if (std::getenv("PROBE_W"))
    torus.project.w_osculating_torus = std::atof(std::getenv("PROBE_W"));
  else if (!std::getenv("PROBE_DUPIN"))
    torus.project.w_osculating_torus = 1.0;
  run("torus", torus, frames);
  return 0;
}
