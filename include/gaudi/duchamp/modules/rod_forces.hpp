#ifndef __GAUDI_DUCHAMP_MODULES_ROD_FORCES__
#define __GAUDI_DUCHAMP_MODULES_ROD_FORCES__

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>

#include "gaudi/asawa/rod/dynamic.hpp"
#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/calder/rod_integrators.hpp"
#include "gaudi/calder/tangent_point_integrators.hpp"
#include "gaudi/common.h"
#include "gaudi/duchamp/field_nodes.hpp"
#include "gaudi/duchamp/utils/sdf.hpp"
#include "liblombardi/junction_node.hpp"

namespace gaudi {
namespace duchamp {

// --- Shared sub-configs (composed by demo configs) ---------------------------

struct tangent_point_force_config {
  real w = 0.0;
  real l0 = 1.0;
  real p = 6.0;
};

struct vortex_force_config {
  real w = 0.0;
  real p = 4.0;
  real q = 0.0;
};

struct boundary_force_config {
  real w = 1.0;
};

struct pca_curve_force_config {
  real w = 0.0;
  int half_width = 4;
};

struct world_space_moment_config {
  real Mx = 0.0;
  real My = 0.0;
  real Mz = 0.0;
};

struct eigen_moment_config {
  real w = 0.0;
  real p = 1.0;       // power on relative-Markley mode eigenvalue
  int half_width = 4;
  int mode = 2;       // 0..2 Markley evec on dquats (ascending); 2 = first away from mean
};

struct local_moment_config {
  real Mbn = 0.0;
  real Mn = 0.0;
  real Mt = 0.0;
};

struct rod_physics_config {
  /// Cosserat stretch–shear (positions + quats). Set 0 to omit those rows.
  real stretch_w = 1e-1;
  /// Position-only 1D rest-length spring (legacy edge_stretch). Set 0 to omit.
  real edge_stretch_w = 0.0;
  real bend_w = 1e-1;
  real twist_w = 1e-1;
  real collision_w = 1.0;
  real smooth_w = 0.0;
  real min_kink_w = 0.0;
  real squad_smooth_w = 0.0;
  real dt = 0.05;
  real damping = 0.1;
  int iterations = 1;
};

// --- Shared rod force helpers (used by graph nodes and demos) ----------------

inline std::vector<vec3> compute_boundary_gradients(const asawa::rod::rod &rod,
                                                    const sdf_base &sdf) {
  const std::vector<real> dists = sdf.distance(rod.__x);
  const std::vector<vec3> gdists = sdf.grad_distance(rod.__x);
  std::vector<vec3> f(rod.__x.size(), vec3::Zero());
  for (int i = 0; i < static_cast<int>(rod.__x.size()); ++i) {
    if (dists[i] > 0.0)
      f[i] = -dists[i] * gdists[i];
  }
  return f;
}

inline std::vector<vec3>
compute_tangent_point_gradient(asawa::rod::rod &rod,
                               const asawa::rod::dynamic &dynamic,
                               real l0_mul = 1.0, real p = 6.0) {
  const real eps = dynamic._Cc;
  const std::vector<vec3> &x = rod.x();
  const std::vector<real> l = rod.l0();
  const std::vector<vec3> T = rod.N2c();
  auto forces =
      calder::tangent_point_gradient(rod, x, l, T, l0_mul * eps, p);
  for (auto &f : forces)
    f *= -1.0;
  return forces;
}

inline std::vector<vec3>
compute_tangent_point_gradient(asawa::rod::rod &rod,
                               const asawa::rod::dynamic &dynamic,
                               const tangent_point_force_config &cfg) {
  auto forces = compute_tangent_point_gradient(rod, dynamic, cfg.l0, cfg.p);
  if (cfg.w != 1.0) {
    for (auto &f : forces)
      f *= cfg.w;
  }
  return forces;
}

// Self-induced vortex: -w*kappa*pow(sin^2,q)*(dp x T). q=0 disables angular filter.
inline std::vector<vec3> compute_vortex_force(asawa::rod::rod &rod,
                                              const asawa::rod::dynamic &dynamic,
                                              real p = 4.0, real q = 1.0) {
  const real eps = dynamic._Cc;
  const std::vector<vec3> &x = rod.x();
  const std::vector<real> phi = rod.l0();
  return calder::vortex_force(rod, x, phi, eps, p, q);
}

inline std::vector<vec3> compute_vortex_force(asawa::rod::rod &rod,
                                              const asawa::rod::dynamic &dynamic,
                                              const vortex_force_config &cfg) {
  auto forces = compute_vortex_force(rod, dynamic, cfg.p, cfg.q);
  if (cfg.w != 1.0) {
    for (auto &f : forces)
      f *= cfg.w;
  }
  return forces;
}

// Local PCA / Álvarez curvature: f_i = w_curve * kappa_i * N_i (not Cosserat).
namespace detail {
struct pca_tube_t {
  vec3 normal = vec3::UnitY();
  real kappa = 0.0;
  bool ok = false;
};

inline pca_tube_t tube_from_local_svd(const std::vector<vec3> &xs,
                                      const vec3 &q0) {
  pca_tube_t out;
  if (xs.size() < 3)
    return out;

  vec3 mean = vec3::Zero();
  for (const vec3 &x : xs)
    mean += x;
  mean /= real(xs.size());

  mat3 C = mat3::Zero();
  for (const vec3 &x : xs) {
    const vec3 d = x - mean;
    C += d * d.transpose();
  }
  C /= real(xs.size());

  Eigen::SelfAdjointEigenSolver<mat3> es(C);
  if (es.info() != Eigen::Success)
    return out;

  const real l0 = std::max(real(0.0), es.eigenvalues()[0]);
  const real l1 = std::max(real(0.0), es.eigenvalues()[1]);
  const real l2 = std::max(real(0.0), es.eigenvalues()[2]);
  const real s1 = std::sqrt(l2);
  const real s2 = std::sqrt(l1);
  if (s1 < 1.0e-14 || s2 < 1.0e-14)
    return out;

  vec3 tangent = es.eigenvectors().col(2).normalized();
  vec3 normal = es.eigenvectors().col(1).normalized();
  vec3 binormal = es.eigenvectors().col(0).normalized();

  if (xs.size() >= 3) {
    const size_t mid = xs.size() / 2;
    const vec3 chord = xs[std::min(mid + 1, xs.size() - 1)] -
                       xs[mid > 0 ? mid - 1 : 0];
    if (chord.squaredNorm() > 1.0e-16 && chord.dot(tangent) < 0.0)
      tangent = -tangent;
  }

  binormal = tangent.cross(normal);
  if (binormal.squaredNorm() < 1.0e-16)
    binormal = es.eigenvectors().col(0).normalized();
  else
    binormal.normalize();
  normal = binormal.cross(tangent).normalized();

  constexpr real k_pre = 1.4907119849998597; // √20 / 3
  const real kappa = k_pre * s2 / (s1 * s1);
  if (kappa < 1.0e-14)
    return out;

  const vec3 to_mean = mean - q0;
  const vec3 rad = to_mean - to_mean.dot(binormal) * binormal;
  if (rad.squaredNorm() > 1.0e-16 && rad.dot(normal) < 0.0)
    normal = -normal;

  out.normal = normal;
  out.kappa = kappa;
  out.ok = true;
  return out;
}

inline bool gather_curve_stencil(const asawa::rod::rod &rod,
                                 asawa::rod::CornerId c, int half_width,
                                 std::vector<vec3> &xs, vec3 &q0) {
  using asawa::rod::CornerId;
  if (half_width < 1)
    return false;
  const std::vector<vec3> &x = rod.x();
  q0 = x[static_cast<size_t>(c)];

  std::vector<CornerId> left;
  left.reserve(static_cast<size_t>(half_width));
  CornerId p = c;
  for (int k = 0; k < half_width; ++k) {
    p = rod.prev(p);
    if (p < 0)
      return false;
    left.push_back(p);
  }

  xs.clear();
  xs.reserve(static_cast<size_t>(2 * half_width + 1));
  for (int k = static_cast<int>(left.size()) - 1; k >= 0; --k)
    xs.push_back(x[static_cast<size_t>(left[static_cast<size_t>(k)])]);
  xs.push_back(q0);

  CornerId n = c;
  for (int k = 0; k < half_width; ++k) {
    n = rod.next(n);
    if (n < 0)
      return false;
    xs.push_back(x[static_cast<size_t>(n)]);
  }
  return true;
}

inline bool gather_quat_stencil(const asawa::rod::rod &rod,
                                asawa::rod::CornerId c, int half_width,
                                std::vector<quat> &us, quat &u0) {
  using asawa::rod::CornerId;
  if (half_width < 1)
    return false;
  const std::vector<quat> &u = rod.u();
  if (static_cast<size_t>(c) >= u.size())
    return false;
  u0 = u[static_cast<size_t>(c)].normalized();

  std::vector<CornerId> left;
  left.reserve(static_cast<size_t>(half_width));
  CornerId p = c;
  for (int k = 0; k < half_width; ++k) {
    p = rod.prev(p);
    if (p < 0)
      return false;
    left.push_back(p);
  }

  us.clear();
  us.reserve(static_cast<size_t>(2 * half_width + 1));
  for (int k = static_cast<int>(left.size()) - 1; k >= 0; --k)
    us.push_back(u[static_cast<size_t>(left[static_cast<size_t>(k)])].normalized());
  us.push_back(u0);

  CornerId n = c;
  for (int k = 0; k < half_width; ++k) {
    n = rod.next(n);
    if (n < 0)
      return false;
    us.push_back(u[static_cast<size_t>(n)].normalized());
  }
  return true;
}

inline quat markley_average(const std::vector<quat> &us) {
  if (us.empty())
    return quat::Identity();
  mat4 A = mat4::Zero();
  vec4 ref = us[0].coeffs();
  if (ref.squaredNorm() < 1.0e-32)
    return quat::Identity();
  for (const quat &qi : us) {
    vec4 c = qi.coeffs();
    if (c.dot(ref) < 0.0)
      c = -c;
    A += c * c.transpose();
  }
  Eigen::SelfAdjointEigenSolver<mat4> es(A);
  if (es.info() != Eigen::Success)
    return us[0];
  quat mean(es.eigenvectors().col(3).data());
  mean.normalize();
  if (mean.coeffs().dot(ref) < 0.0)
    mean.coeffs() = -mean.coeffs();
  return mean;
}

// Full Markley eigensystem of Σ q qᵀ (ascending evals). col(3) = mean.
struct markley_modes_t {
  quat mean = quat::Identity();
  mat4 evecs = mat4::Identity(); // columns = eigenvectors
  vec4 evals = vec4::Zero();
  bool ok = false;
};

inline markley_modes_t markley_modes(const std::vector<quat> &us) {
  markley_modes_t out;
  if (us.empty())
    return out;
  mat4 A = mat4::Zero();
  vec4 ref = us[0].coeffs();
  if (ref.squaredNorm() < 1.0e-32)
    return out;
  for (const quat &qi : us) {
    vec4 c = qi.coeffs();
    if (c.dot(ref) < 0.0)
      c = -c;
    A += c * c.transpose();
  }
  Eigen::SelfAdjointEigenSolver<mat4> es(A);
  if (es.info() != Eigen::Success)
    return out;
  out.evecs = es.eigenvectors();
  out.evals = es.eigenvalues();
  out.mean = quat(out.evecs.col(3).data());
  out.mean.normalize();
  if (out.mean.coeffs().dot(ref) < 0.0)
    out.mean.coeffs() = -out.mean.coeffs();
  out.ok = true;
  return out;
}

// Map Markley tangent δq ⊥ mean (R⁴) → world so(3) axis.
inline vec3 markley_tangent_to_world_axis(const quat &mean, const vec4 &dq) {
  quat d;
  d.coeffs() = dq;
  const quat body = mean.conjugate() * d;
  vec3 a = mean * (real(2.0) * body.vec());
  const real n2 = a.squaredNorm();
  if (n2 < 1.0e-32)
    return vec3::Zero();
  return a / std::sqrt(n2);
}

inline vec3 quat_log_vec(quat q) {
  q.normalize();
  if (q.w() < 0.0)
    q.coeffs() = -q.coeffs();
  const vec3 v = q.vec();
  const real vn = v.norm();
  if (vn < 1.0e-14)
    return vec3::Zero();
  const real w = std::max(real(-1.0), std::min(real(1.0), q.w()));
  const real angle = real(2.0) * std::atan2(vn, w);
  return (angle / vn) * v;
}

inline quat hemisphere_align(const quat &ref, quat q) {
  if (q.coeffs().dot(ref.coeffs()) < 0.0)
    q.coeffs() = -q.coeffs();
  return q;
}
} // namespace detail

inline std::vector<vec3> compute_pca_curve_force(const asawa::rod::rod &rod,
                                                 real w_curve,
                                                 int half_width = 4) {
  using asawa::rod::corner_id;
  std::vector<vec3> f(rod.x().size(), vec3::Zero());
  if (w_curve == 0.0 || half_width < 1)
    return f;

  std::vector<vec3> xs;
  for (int i = 0; i < rod.corner_count(); ++i) {
    vec3 q0;
    if (!detail::gather_curve_stencil(rod, corner_id(i), half_width, xs, q0))
      continue;
    const detail::pca_tube_t tube = detail::tube_from_local_svd(xs, q0);
    if (!tube.ok)
      continue;
    f[static_cast<size_t>(i)] = w_curve * tube.kappa * tube.normal;
  }
  return f;
}

inline std::vector<vec3>
compute_pca_curve_force(const asawa::rod::rod &rod,
                        const pca_curve_force_config &cfg) {
  return compute_pca_curve_force(rod, cfg.w, cfg.half_width);
}

// Full 4×4 Markley on relative dquats (Frenet-force analogue):
//   1) absolute Markley mean ū
//   2) δq_j = ū^{-1} u_j
//   3) 4×4 Markley on {δq} → mode dQ (default v2) + μ
//   4) torque axis from dQ (Im map in ū-body → world)
//   5) sign from center residual
inline std::vector<vec3>
compute_pca_curve_torque(const asawa::rod::rod &rod,
                         const eigen_moment_config &cfg,
                         const std::vector<vec3> *tau_tie = nullptr) {
  (void)tau_tie;
  using asawa::rod::corner_id;
  std::vector<vec3> tau(rod.u().size(), vec3::Zero());
  if (cfg.w == 0.0 || cfg.half_width < 1)
    return tau;

  const int mode = std::max(0, std::min(2, cfg.mode));
  std::vector<quat> us;
  std::vector<quat> dqs;

  for (int i = 0; i < rod.corner_count(); ++i) {
    quat u0;
    if (!detail::gather_quat_stencil(rod, corner_id(i), cfg.half_width, us, u0))
      continue;

    const quat mean_abs = detail::markley_average(us);

    dqs.clear();
    dqs.reserve(us.size());
    for (const quat &uj : us) {
      const quat uj_a = detail::hemisphere_align(mean_abs, uj);
      dqs.push_back(
          detail::hemisphere_align(quat::Identity(),
                                   mean_abs.conjugate() * uj_a));
    }

    const detail::markley_modes_t mk = detail::markley_modes(dqs);
    if (!mk.ok)
      continue;

    const real mu = std::max(real(0.0), mk.evals[mode]);
    const real N = real(dqs.size());
    const real lam = (N > 0.0) ? (mu / N) : 0.0;
    if (lam < 1.0e-16)
      continue;

    // dQ mode → axis in ū-body (relatives left-trivialized at mean_abs), then world.
    const vec3 a_body =
        detail::markley_tangent_to_world_axis(mk.mean, mk.evecs.col(mode));
    if (a_body.squaredNorm() < 1.0e-32)
      continue;
    vec3 axis = mean_abs * a_body;
    const real an2 = axis.squaredNorm();
    if (an2 < 1.0e-32)
      continue;
    axis /= std::sqrt(an2);

    const quat u0_a = detail::hemisphere_align(mean_abs, u0);
    const vec3 d_center =
        mean_abs.toRotationMatrix() *
        detail::quat_log_vec(mean_abs.conjugate() * u0_a);
    if (axis.dot(d_center) < 0.0)
      axis = -axis;

    tau[static_cast<size_t>(i)] = cfg.w * std::pow(lam, cfg.p) * axis;
  }
  return tau;
}

inline std::vector<vec3>
apply_world_space_moment(std::vector<vec3> tau,
                         const world_space_moment_config &cfg) {
  const vec3 M(cfg.Mx, cfg.My, cfg.Mz);
  if (M.squaredNorm() == 0.0)
    return tau;
  for (auto &t : tau)
    t += M;
  return tau;
}

inline sdf_base::ptr select_rod_sdf(int frame, const sdf_base::ptr &sdf0,
                                    const sdf_base::ptr &sdf1) {
  (void)frame;
  (void)sdf1;
  return sdf0;
}

inline void grow_rod_rest_lengths(asawa::rod::rod &rod, real factor) {
  std::vector<real> &l0 = rod.l0();
  for (auto &li : l0)
    li *= factor;
}

// --- Lombardi graph nodes ----------------------------------------------------

class boundary_gradient_node : public liblombardi::Node {
public:
  using ptr = std::shared_ptr<boundary_gradient_node>;

  enum class PortId { Output = 0 };
  using OutputPortDef = liblombardi::PortDef<field_datum<vec3>, PortId::Output>;

  boundary_gradient_node(asawa::rod::rod::ptr rod, sdf_base::ptr sdf0,
                         sdf_base::ptr sdf1,
                         const boundary_force_config &cfg = {})
      : _rod(std::move(rod)), _sdf0(std::move(sdf0)), _sdf1(std::move(sdf1)),
        _cfg(cfg) {}

  boundary_gradient_node(asawa::rod::rod::ptr rod, sdf_base::ptr sdf0,
                         sdf_base::ptr sdf1, real scale)
      : boundary_gradient_node(std::move(rod), std::move(sdf0), std::move(sdf1),
                               boundary_force_config{scale}) {}

  void set_frame(int frame) { _frame = frame; }

  void compute() override {
    auto out = get_datum<OutputPortDef>();
    const auto sdf = select_rod_sdf(_frame, _sdf0, _sdf1);
    auto forces = compute_boundary_gradients(*_rod, *sdf);
    if (_cfg.w != 1.0) {
      for (auto &g : forces)
        g *= _cfg.w;
    }
    out->data() = std::move(forces);
  }

  unsigned int port_count() const override { return 1; }

  liblombardi::PortRef<boundary_gradient_node, OutputPortDef> output() {
    return {*this};
  }

private:
  asawa::rod::rod::ptr _rod;
  sdf_base::ptr _sdf0;
  sdf_base::ptr _sdf1;
  boundary_force_config _cfg;
  int _frame = 0;
};

class tangent_point_gradient_node : public liblombardi::Node {
public:
  using ptr = std::shared_ptr<tangent_point_gradient_node>;

  enum class PortId { Output = 0 };
  using OutputPortDef = liblombardi::PortDef<field_datum<vec3>, PortId::Output>;

  tangent_point_gradient_node(asawa::rod::rod::ptr rod,
                              asawa::rod::dynamic::ptr dynamic,
                              const tangent_point_force_config &cfg = {})
      : _rod(std::move(rod)), _dynamic(std::move(dynamic)), _cfg(cfg) {}

  // scale: force weight (negative ⇒ clumping). l0_mul: kernel length / Cc
  tangent_point_gradient_node(asawa::rod::rod::ptr rod,
                              asawa::rod::dynamic::ptr dynamic, real scale = 1.0,
                              real l0_mul = 1.0, real p = 6.0)
      : tangent_point_gradient_node(
            std::move(rod), std::move(dynamic),
            tangent_point_force_config{scale, l0_mul, p}) {}

  void compute() override {
    auto out = get_datum<OutputPortDef>();
    out->data() = compute_tangent_point_gradient(*_rod, *_dynamic, _cfg);
  }

  unsigned int port_count() const override { return 1; }

  liblombardi::PortRef<tangent_point_gradient_node, OutputPortDef> output() {
    return {*this};
  }

private:
  asawa::rod::rod::ptr _rod;
  asawa::rod::dynamic::ptr _dynamic;
  tangent_point_force_config _cfg;
};

class vortex_force_node : public liblombardi::Node {
public:
  using ptr = std::shared_ptr<vortex_force_node>;

  enum class PortId { Output = 0 };
  using OutputPortDef = liblombardi::PortDef<field_datum<vec3>, PortId::Output>;

  vortex_force_node(asawa::rod::rod::ptr rod, asawa::rod::dynamic::ptr dynamic,
                    const vortex_force_config &cfg = {})
      : _rod(std::move(rod)), _dynamic(std::move(dynamic)), _cfg(cfg) {}

  vortex_force_node(asawa::rod::rod::ptr rod, asawa::rod::dynamic::ptr dynamic,
                    real scale = 1.0, real p = 4.0, real q = 1.0)
      : vortex_force_node(std::move(rod), std::move(dynamic),
                          vortex_force_config{scale, p, q}) {}

  void compute() override {
    auto out = get_datum<OutputPortDef>();
    out->data() = compute_vortex_force(*_rod, *_dynamic, _cfg);
  }

  unsigned int port_count() const override { return 1; }

  liblombardi::PortRef<vortex_force_node, OutputPortDef> output() {
    return {*this};
  }

private:
  asawa::rod::rod::ptr _rod;
  asawa::rod::dynamic::ptr _dynamic;
  vortex_force_config _cfg;
};

class pca_curve_force_node : public liblombardi::Node {
public:
  using ptr = std::shared_ptr<pca_curve_force_node>;

  enum class PortId { Output = 0 };
  using OutputPortDef = liblombardi::PortDef<field_datum<vec3>, PortId::Output>;

  pca_curve_force_node(asawa::rod::rod::ptr rod,
                       const pca_curve_force_config &cfg = {})
      : _rod(std::move(rod)), _cfg(cfg) {}

  pca_curve_force_node(asawa::rod::rod::ptr rod, real w_curve = 0.0,
                       int half_width = 4)
      : pca_curve_force_node(std::move(rod),
                             pca_curve_force_config{w_curve, half_width}) {}

  void compute() override {
    auto out = get_datum<OutputPortDef>();
    out->data() = compute_pca_curve_force(*_rod, _cfg);
  }

  unsigned int port_count() const override { return 1; }

  liblombardi::PortRef<pca_curve_force_node, OutputPortDef> output() {
    return {*this};
  }

private:
  asawa::rod::rod::ptr _rod;
  pca_curve_force_config _cfg;
};

template <int N>
using vec3_junction_node =
    liblombardi::junction_node<N, field_datum<vec3>, liblombardi::add_op<vec3>>;

} // namespace duchamp
} // namespace gaudi

#endif
