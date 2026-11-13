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
    const auto &r = demo->last_solve();
    std::cout << label << "  [" << f << "] tube r " << s.mean << " roughness " << s.rel_std << " max|x| "
              << s.max_abs << " verts " << s.verts << " E " << r.E_first << " -> " << r.E_last << " s " << dt
              << "\n";
  }
}

} // namespace

int main(int argc, char **argv) {
  gaudi::hepworth::block::osculating_torus_ring_stats::enabled = true;
  const bool relative = std::getenv("PROBE_RELATIVE") != nullptr;
  const bool closest = std::getenv("PROBE_CLOSEST") != nullptr;
  const int frames = argc > 1 ? std::atoi(argv[1]) : 8;
  const gaudi::real perturb = argc > 2 ? std::atof(argv[2]) : 0.3;
  gaudi::duchamp::monge_fairing_demo_config base{};
  base.perturb = perturb;

  if (std::getenv("PROBE_ALL"))
    run("base", base, frames);

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
    base.position.w_bend = 0.0;
  }
  base.project.osculating_torus.relative = relative;
  base.project.osculating_torus.closest_point = closest;
  if (const char *v = std::getenv("PROBE_VOLUME")) {
    using gaudi::hepworth::block::ring_volume_mode;
    auto &vol = base.project.osculating_torus.volume;
    vol = {};
    vol.mode = std::string(v) == "off" ? ring_volume_mode::off : ring_volume_mode::corotational;
    if (std::getenv("PROBE_VOLUME_W"))
      vol.corotational.w_volume = std::atof(std::getenv("PROBE_VOLUME_W"));
    if (std::getenv("PROBE_FOOT_K"))
      vol.corotational.foot_stiffness = std::atof(std::getenv("PROBE_FOOT_K"));
    if (std::getenv("PROBE_REST_W"))
      vol.corotational.w_rest = std::atof(std::getenv("PROBE_REST_W"));
    if (std::getenv("PROBE_FOOT_W"))
      vol.foot_weight = std::atof(std::getenv("PROBE_FOOT_W"));
  }
  auto torus = base;
  torus.project.w_osculating_torus = std::getenv("PROBE_W") ? std::atof(std::getenv("PROBE_W")) : 1.0;
  run("torus", torus, frames);
  return 0;
}
