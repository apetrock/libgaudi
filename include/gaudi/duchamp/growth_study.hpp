#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "gaudi/vec_addendum.h"

#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/dynamic.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/asawa/shell/triangulation_debug.hpp"
#include "gaudi/asawa/datums.hpp"

#include "gaudi/kusama/laplacian.hpp"
#include "gaudi/kusama/laplacian_anisotropic.hpp"
#include "gaudi/calder/shell_integrators.hpp"
#include "gaudi/duchamp/modules/rx_diffuse.hpp"

#include "gaudi/hepworth/block/shell_constraints.hpp"
#include "gaudi/hepworth/block/shell_constraints_init.hpp"
#include "gaudi/hepworth/block/sim_block.hpp"
#include "gaudi/hepworth/block/solver.hpp"
#include "gaudi/hepworth/block/block_geometry.hpp"
#include "gaudi/hepworth/constraints/slab/relative_plane_twist.hpp"
#include "gaudi/hepworth/constraints/slab/rbf_elastic_strain.hpp"
#include "gaudi/hepworth/block/solver_composition.hpp"
#include "gaudi/hepworth/constraints/slab/slab_pierce_stencil.hpp"
#include "gaudi/hepworth/constraints/slab/slab_prism_config.hpp"

#include "gaudi/common.h"
#include "gaudi/duchamp/pierce_cylinder_medial.hpp"
#include "gaudi/calder/tangent_point_integrators.hpp"
#include "gaudi/geometry_logger.hpp"

#include "modules/rx/ginzburg_landau.hpp"
#include "modules/rx/grey_scott.hpp"
#include "modules/rx_colormap.hpp"
#include "modules/rx/swift_hohenberg.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <vector>

#ifndef __M2REFACTOR_TEST__
#define __M2REFACTOR_TEST__

namespace gaudi {
namespace duchamp {

using namespace asawa;

/// Selects which reaction–diffusion module backs \ref growth_study.
enum class growth_rx_model {
  grey_scott,
  swift_hohenberg,
  ginzburg_landau,
};

/// Growth study: RD colors the mesh; growth comes from pierce-medial axial/hoop
/// via anisotropic triangle strain + relative_plane_twist and/or rbf_elastic,
/// plus a mild axial force.
///
/// Live path per frame:
///   1. step_rx        — anisotropic RD (colors only for visualization)
///   2. step_dynamics  — heat distance, face-probe medial axes,
///                       strain + twist/RBF slab + pnt–tri collisions
///   3. __surf->step   — remesh / dynamic shell
class growth_study {
public:
  typedef std::shared_ptr<growth_study> ptr;

  static ptr create(growth_rx_model model = growth_rx_model::grey_scott) {
    return std::make_shared<growth_study>(model);
  }

  static bool should_trace_frame(index_t frame) {
    return frame < 3 || frame % 60 == 0;
  }

  void set_nan_probe(bool on, int frame_stride = 1) {
    _nan_probe = on;
    _nan_probe_stride = std::max(1, frame_stride);
  }

  /// Additive axial/hoop growth rates. Stretch = 1 + heat * _growth.
  /// e.g. (0.1, -0.5) => +10% axial, -50% hoop at heat=1; (0,0) is rest.
  void set_growth(const vec2 &g) { _growth = g; }
  const vec2 &growth() const { return _growth; }

  explicit growth_study(growth_rx_model model = growth_rx_model::grey_scott)
      : _rx_model(model) {
    std::cerr << "[growth_study_core] ctor start" << std::endl;
    __M = shell::load_bunny();
    std::cerr << "[growth_study_core] mesh loaded" << std::endl;

    {
      const std::vector<vec3> &xpre = asawa::const_get_vec_data(*__M, 0);
      asawa::shell::dump_pre_triangulation_report(*__M, std::cerr, &xpre);
    }
    shell::triangulate(*__M);
    std::cerr << "[growth_study_core] triangulated" << std::endl;
    for (int i = 0; i < __M->face_count(); i++) {
      asawa::shell::FaceId fi = asawa::shell::face_id(i);
      if (__M->fbegin(fi) >= 0) {
        assert(__M->fsize(fi) == 3);
      }
    }

    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    asawa::center(x);

    real l0 = 4.0 * asawa::shell::avg_length(*__M, x);
    _eps = l0;
    __surf = shell::dynamic::create(__M, 1.0 * l0, 2.5 * l0, 0.5 * l0);
    std::cerr << "[growth_study_core] shell dynamic created" << std::endl;

    for (int i = 0; i < 4; ++i) {
      std::cout << "[growth_study_core] stepping the surface for warmup: " << i
                << std::endl;
      __surf->step(true);
    }
    std::cerr << "[growth_study_core] warmup steps done" << std::endl;

    real f0 = 0.031, k0 = 0.0585;
    real da0 = 3.0e-4, db0 = 0.5 * da0;
    switch (_rx_model) {
    case growth_rx_model::grey_scott: {
      rx::grey_scott::ptr rx0 =
          rx::grey_scott::create(__M, f0, k0, da0, db0);
      _rx0 = std::dynamic_pointer_cast<module_base>(rx0);
      std::cerr << "[growth_study_core] reaction diffusion (Grey–Scott) created"
                << std::endl;
      break;
    }
    case growth_rx_model::swift_hohenberg: {
      rx::swift_hohenberg::ptr sh = rx::swift_hohenberg::create(__M, 0.04, 1.0);
      _rx0 = std::dynamic_pointer_cast<module_base>(sh);
      std::cerr << "[growth_study_core] Swift–Hohenberg created" << std::endl;
      break;
    }
    case growth_rx_model::ginzburg_landau: {
      rx::ginzburg_landau::ptr gl = rx::ginzburg_landau::create(__M, 1.2, 1.0);
      _rx0 = std::dynamic_pointer_cast<module_base>(gl);
      std::cerr << "[growth_study_core] Ginzburg–Landau created" << std::endl;
      break;
    }
    }

    init_origin();
    std::cerr << "[growth_study_core] ctor complete" << std::endl;
  }

  void init_origin() {
    const std::vector<vec3> &x = asawa::const_get_vec_data(*__M, 0);
    std::array<vec3, 2> ext = asawa::extents(x);
    _origin = vec3(0.5 * (ext[0][0] + ext[1][0]), ext[0][1],
                   0.5 * (ext[0][2] + ext[1][2]));
  }

  void set_mesh_color_gradient(std::function<vec4(real)> g) {
    _mesh_color_gradient = std::move(g);
  }

  std::vector<vec4> get_mesh_colors() {
    auto default_grad = [](real t) {
      return vec4(t, 0.15 + 0.85 * (1.0 - t), 1.0 - 0.4 * t, 1.0);
    };
    std::function<vec4(real)> grad =
        _mesh_color_gradient ? *_mesh_color_gradient : default_grad;

    switch (_rx_model) {
    case growth_rx_model::grey_scott: {
      std::vector<vec4> colors(__M->vert_count(), vec4(1.0, 0.0, 0.0, 1.0));
      std::vector<real> &rx0a =
          std::dynamic_pointer_cast<rx::grey_scott>(_rx0)->get_rxa();
      std::vector<real> &rx0b =
          std::dynamic_pointer_cast<rx::grey_scott>(_rx0)->get_rxb();
      vec4 col_a(1.0, 0.0, 1.0, 1.0);
      vec4 col_b(0.0, 1.0, 1.0, 1.0);
      for (int k = 0; k < __M->vert_count(); k++) {
        colors[k] = 1.0 * rx0a[k] * col_a + 3.0 * rx0b[k] * col_b;
      }
      return colors;
    }
    case growth_rx_model::swift_hohenberg: {
      std::vector<real> u =
          std::dynamic_pointer_cast<rx::swift_hohenberg>(_rx0)->get_u();
      return field_colors_minmax(u, grad);
    }
    case growth_rx_model::ginzburg_landau: {
      auto gl = std::dynamic_pointer_cast<rx::ginzburg_landau>(_rx0);
      std::vector<real> mag = gl->amplitude_abs();
      return field_colors_minmax(mag, grad);
    }
    }
    return std::vector<vec4>(__M->vert_count(), vec4(1, 0, 0, 1));
  }

  vec3 get_origin() { return _origin; }

  index_t get_new_origin(asawa::shell::shell &shell) {
    const std::vector<vec3> &x = asawa::const_get_vec_data(shell, 0);
    index_t imin = 0;
    real min = std::numeric_limits<real>::max();
    for (int i = 0; i < static_cast<int>(x.size()); i++) {
      real d = (x[i] - _origin).norm();
      if (d < min) {
        min = d;
        imin = i;
      }
    }
    _origin = x[imin];
    return imin;
  }

  std::vector<real> vertex_geodesic_weight(asawa::shell::shell &shell) {
    const std::vector<vec3> &x = asawa::const_get_vec_data(shell, 0);
    index_t imin = get_new_origin(shell);
    std::vector<real> f(shell.vert_count(), 0.0);
    f[imin] = 1.0;
    kusama::laplacian L(__M, x);
    return L.heatDist(f, 0.2);
  }

  void calc_sh_params(const std::vector<real> &t) {
    real e0 = 0.02, e1 = 0.07;
    real g0 = 0.8, g1 = 1.6;
    real l0 = 2.0e-4, l1 = 7.0e-4;
    _eps_sh.resize(t.size());
    _g_sh.resize(t.size());
    _lam_sh.resize(t.size());
    for (int i = 0; i < static_cast<int>(t.size()); i++) {
      _eps_sh[i] = va::mix(t[i], e0, e1);
      _g_sh[i] = va::mix(t[i], g0, g1);
      _lam_sh[i] = va::mix(t[i], l0, l1);
    }
  }

  void calc_gl_params(const std::vector<real> &t) {
    real a0 = 0.8, a1 = 1.8;
    real b0 = 0.6, b1 = 1.4;
    _alpha_gl.resize(t.size());
    _beta_gl.resize(t.size());
    for (int i = 0; i < static_cast<int>(t.size()); i++) {
      _alpha_gl[i] = va::mix(t[i], a0, a1);
      _beta_gl[i] = va::mix(t[i], b0, b1);
    }
  }

  void calc_kf(const std::vector<real> &t) {
    real f0 = 0.034, k0 = 0.059; // fingerprints
    real f1 = 0.0257, k1 = 0.0555; // mazy
    _k = std::vector<real>(t.size(), 0.0);
    _f = std::vector<real>(t.size(), 0.0);
    for (int i = 0; i < static_cast<int>(t.size()); i++) {
      _k[i] = va::mix(t[i], k0, k1);
      _f[i] = va::mix(t[i], f0, f1);
    }
  }

  void step_rx(int frame) {
    if (should_trace_frame(frame)) {
      std::cerr << "[growth_study_core] step_rx begin frame=" << frame
                << std::endl;
    }

    const std::vector<vec3> &x_pos = asawa::const_get_vec_data(*__M, 0);
    _C_aniso = kusama::build_curvature_aligned_laplacian(
        *__M, x_pos, _aniso_vd_lambda, _aniso_sigma_u, _aniso_sigma_v,
        asawa::shell::face_curvature_stencil::one_ring, _aniso_kind);

    if (_rx_model == growth_rx_model::grey_scott) {
      std::dynamic_pointer_cast<rx::grey_scott>(_rx0)->set_diffuse_smooth(
          std::make_optional<rx_smooth_fn>([this](std::vector<real> &f, real dt) {
            const std::vector<vec3> &xp = asawa::const_get_vec_data(*__M, 0);
            rx_diffuse_scalar_implicit_cotan(__M, xp, f, dt, &_C_aniso);
          }));
    } else if (_rx_model == growth_rx_model::swift_hohenberg) {
      std::dynamic_pointer_cast<rx::swift_hohenberg>(_rx0)->set_diffuse_smooth(
          std::make_optional<rx_smooth_fn>([this](std::vector<real> &f, real dt) {
            const std::vector<vec3> &xp = asawa::const_get_vec_data(*__M, 0);
            rx_diffuse_scalar_implicit_cotan(__M, xp, f, dt, &_C_aniso);
          }));
    } else {
      std::dynamic_pointer_cast<rx::ginzburg_landau>(_rx0)->set_dispersive_stiffness(
          &_C_aniso);
    }

    const std::vector<real> d = vertex_geodesic_weight(*__M);
    calc_kf(d);
    calc_sh_params(d);
    calc_gl_params(d);

    const int N_gs = 10;
    const real h_gs = 16.0 * _h;
    const int N_sh = 1;
    const real h_sh = 0.01;
    const int N_gl = 1;
    const real h_gl = 0.006;

    switch (_rx_model) {
    case growth_rx_model::grey_scott:
      for (int i = 0; i < N_gs; i++) {
        std::dynamic_pointer_cast<rx::grey_scott>(_rx0)->step_anisotropic(
            h_gs, _f, _k, nullptr);
      }
      break;
    case growth_rx_model::swift_hohenberg:
      for (int i = 0; i < N_sh; i++) {
        std::dynamic_pointer_cast<rx::swift_hohenberg>(_rx0)->step(
            h_sh, _eps_sh, _g_sh, _lam_sh, nullptr);
      }
      break;
    case growth_rx_model::ginzburg_landau:
      for (int i = 0; i < N_gl; i++) {
        std::dynamic_pointer_cast<rx::ginzburg_landau>(_rx0)->step(
            h_gl, _alpha_gl, _beta_gl, nullptr);
      }
      break;
    }
    if (should_trace_frame(frame)) {
      std::cerr << "[growth_study_core] step_rx end" << std::endl;
    }
  }

  /// Per-face pierce-cylinder field used for growth + trunkification.
  struct medial_face_field {
    std::vector<vec3> axial;  ///< faired tangential axis (unit)
    std::vector<vec3> center; ///< MLS cylinder center
    std::vector<vec3> axis;   ///< faired world axis (unit; for radial foot)
    std::vector<real> radius; ///< MLS cylinder radius
    std::vector<char> valid;  ///< true when center/axis/radius usable
  };

  /// Tangential projection of pierce-medial axis onto each sample normal.
  std::vector<vec3> sample_medial_axial_dirs(
      const std::vector<pierce_cylinder_sample> &medial,
      const std::vector<vec3> &normals) const {
    std::vector<vec3> axial(normals.size(), vec3::Zero());
    const size_t n = std::min(medial.size(), normals.size());
    for (size_t i = 0; i < n; ++i) {
      if (medial[i].axis.squaredNorm() < 1e-24)
        continue;
      vec3 n = normals[i];
      if (!n.allFinite() || n.squaredNorm() < 1e-24)
        continue;
      n.normalize();
      vec3 a = medial[i].axis - medial[i].axis.dot(n) * n;
      if (!a.allFinite() || a.squaredNorm() < 1e-24)
        continue;
      axial[i] = a.normalized();
    }
    return axial;
  }

  /// Single cylinder MLS at face centers, then:
  ///   1) SVD-fair axes
  ///   2) second cylinder MLS that re-fits (C, R) from neighbor prior (c,a,r)
  ///      with axis held fixed (no mesh point/normal rows).
  medial_face_field compute_faired_medial_face_field(
      const std::vector<vec3> &x, pierce_cylinder_medial_stats *stats_out) {
    medial_face_field out;
    const std::vector<vec3> x_face = asawa::shell::face_centers(*__M, x);
    const std::vector<vec3> N_face = asawa::shell::face_normals(*__M, x);
    const std::vector<real> areas = asawa::shell::face_areas(*__M, x);

    pierce_cylinder_medial_params mparams =
        default_pierce_cylinder_medial_params();
    mparams.mls_sigma_pierce_scale = _cylinder_mls_sigma_scale;
    std::vector<pierce_cylinder_sample> medial_face =
        compute_pierce_cylinder_medial_at(*__M, x_face, N_face, mparams,
                                          stats_out);
    std::vector<vec3> axial_raw =
        sample_medial_axial_dirs(medial_face, N_face);

    const size_t n = medial_face.size();
    std::vector<vec3> c_raw(n, vec3::Zero());
    std::vector<vec3> a_raw(n, vec3::Zero());
    std::vector<real> r_raw(n, 0.0);
    std::vector<real> w_face(n, 0.0);
    for (size_t i = 0; i < n; ++i) {
      const pierce_cylinder_sample &s = medial_face[i];
      if (!s.accepted || !s.center.allFinite() || !s.axis.allFinite() ||
          s.axis.squaredNorm() < 1e-24 || !std::isfinite(s.radius) ||
          s.radius <= 1e-8) {
        continue;
      }
      c_raw[i] = s.center;
      a_raw[i] = s.axis.normalized();
      r_raw[i] = s.radius;
      w_face[i] = i < areas.size() ? std::max(areas[i], real(0.0)) : 1.0;
    }

    const real avg_len = asawa::shell::avg_length(*__M, x);
    // Axis SVD still lives on the surface → edge-length bandwidth.
    const real l0_axis =
        std::max(_axial_fair_l0_scale * avg_len, real(1e-12));

    // Fair axes first (covariance eigen / SVD).
    const std::vector<mat3> frames = calder::gaussian_covariant_vector_frame(
        *__M, axial_raw, x_face, l0_axis);
    std::vector<vec3> a_fair(n, vec3::Zero());
    for (size_t i = 0; i < n; ++i) {
      if (w_face[i] <= 0.0)
        continue;
      vec3 a = a_raw[i];
      if (i < frames.size()) {
        const mat3 &F = frames[i];
        if (F.allFinite() && F.squaredNorm() > 1e-24)
          a = F.col(2);
      }
      if (a.squaredNorm() < 1e-24)
        a = a_raw[i];
      if (a_raw[i].squaredNorm() > 1e-24 && a.dot(a_raw[i]) < 0.0)
        a = -a;
      if (a.squaredNorm() > 1e-24)
        a_fair[i] = a.normalized();
    }

    // Second cylinder MLS in center space: σᵢ = radius_scale · rᵢ.
    const std::vector<calder::vec7> faired = calder::cylinder_fair_from_priors(
        *__M, x_face, a_fair, c_raw, r_raw, w_face, _cylinder_fair_radius_scale);

    out.axial.assign(n, vec3::Zero());
    out.center.assign(n, vec3::Zero());
    out.axis.assign(n, vec3::Zero());
    out.radius.assign(n, 0.0);
    out.valid.assign(n, 0);

    for (size_t i = 0; i < n; ++i) {
      vec3 c = vec3::Zero();
      vec3 a = a_fair[i];
      real r = 0.0;
      if (i < faired.size() && faired[i].allFinite() && faired[i][6] > 1e-8) {
        c = faired[i].segment(0, 3);
        if (faired[i].segment(3, 3).squaredNorm() > 1e-24)
          a = faired[i].segment(3, 3).normalized();
        r = faired[i][6];
      } else if (w_face[i] > 0.0) {
        c = c_raw[i];
        r = r_raw[i];
      } else {
        continue;
      }
      if (!c.allFinite() || !a.allFinite() || a.squaredNorm() < 1e-24 ||
          !std::isfinite(r) || r <= 1e-8)
        continue;

      out.center[i] = c;
      out.radius[i] = r;
      out.axis[i] = a;
      out.valid[i] = 1;

      vec3 a_t = a;
      if (i < N_face.size() && N_face[i].squaredNorm() > 1e-24) {
        const vec3 nn = N_face[i].normalized();
        a_t = a - a.dot(nn) * nn;
      }
      if (a_t.squaredNorm() < 1e-24)
        a_t = axial_raw[i];
      if (a_t.squaredNorm() > 1e-24) {
        a_t.normalize();
        if (axial_raw[i].squaredNorm() > 1e-24 && a_t.dot(axial_raw[i]) < 0.0)
          a_t = -a_t;
        out.axial[i] = a_t;
      }
    }
    return out;
  }

  /// Mild heat-modulated force along axial, oriented into the heat-flow half-space.
  /// Heat flow is the Crane heat-method direction X = −∇u (same as heatFrame col 0),
  /// not ∇d of the normalized distance — that field is noisier and can anti-align locally.
  std::vector<vec3> axial_face_forces(const std::vector<vec3> &axial_face,
                                      const std::vector<vec3> &heat_flow_face,
                                      const std::vector<real> &d_face) const {
    const size_t n = std::min(axial_face.size(),
                              std::min(heat_flow_face.size(), d_face.size()));
    std::vector<vec3> f_face(n, vec3::Zero());
    for (size_t i = 0; i < n; ++i) {
      vec3 a = axial_face[i];
      const vec3 &g = heat_flow_face[i];
      if (!a.allFinite() || a.squaredNorm() < 1e-24 || !std::isfinite(d_face[i]))
        continue;
      if (!g.allFinite() || g.squaredNorm() < 1e-20)
        continue;
      a.normalize();
      // Flip into the heat-flow half-space so a · X >= 0.
      if (a.dot(g) < 0.0)
        a = -a;
      f_face[i] = (_axial_force_weight * d_face[i]) * a;
    }
    return f_face;
  }

  /// Face-tangent heat flow X = −∇u from a diffused source (heat-method direction).
  std::vector<vec3> heat_flow_face_dirs(kusama::laplacian &L,
                                        const std::vector<real> &fh,
                                        const std::vector<vec3> &x) const {
    const std::vector<real> u = L.diffuse(fh, 0.2);
    std::vector<vec3> g = asawa::shell::gradient(*__M, u, x);
    const std::vector<vec3> N = asawa::shell::face_normals(*__M, x);
    for (size_t i = 0; i < g.size(); ++i) {
      if (!g[i].allFinite() || g[i].squaredNorm() < 1e-24) {
        g[i] = vec3::Zero();
        continue;
      }
      // X = −∇u, then project to the face tangent plane.
      vec3 flow = -g[i];
      if (i < N.size() && N[i].squaredNorm() > 1e-24) {
        const vec3 nn = N[i].normalized();
        flow = flow - flow.dot(nn) * nn;
      }
      if (flow.squaredNorm() < 1e-24) {
        g[i] = vec3::Zero();
        continue;
      }
      g[i] = flow.normalized();
    }
    return g;
  }

  /// Flip face axials into the heat-flow half-space (shared by force + strain).
  void orient_axials_to_heat_flow(std::vector<vec3> &axial_face,
                                  const std::vector<vec3> &heat_flow_face) const {
    const size_t n = std::min(axial_face.size(), heat_flow_face.size());
    for (size_t i = 0; i < n; ++i) {
      vec3 &a = axial_face[i];
      const vec3 &g = heat_flow_face[i];
      if (!a.allFinite() || a.squaredNorm() < 1e-24)
        continue;
      if (!g.allFinite() || g.squaredNorm() < 1e-20)
        continue;
      if (a.dot(g) < 0.0)
        a = -a;
    }
  }

  /// Trunkification: radial residual step (r − ρ) ρ̂, then project onto the
  /// face normal and displace only along N: f = ((step · N) N).
  std::vector<vec3> trunkification_face_forces(
      const std::vector<vec3> &x_face, const std::vector<vec3> &N_face,
      const medial_face_field &med) const {
    const size_t n = x_face.size();
    std::vector<vec3> f_face(n, vec3::Zero());
    if (_trunk_weight <= 0.0)
      return f_face;
    for (size_t i = 0; i < n; ++i) {
      if (i >= med.valid.size() || !med.valid[i])
        continue;
      const vec3 &c = med.center[i];
      const vec3 &a = med.axis[i];
      const real r = med.radius[i];
      if (!c.allFinite() || !a.allFinite() || a.squaredNorm() < 1e-24 ||
          !std::isfinite(r) || r <= 1e-8)
        continue;
      if (i >= N_face.size() || !N_face[i].allFinite() ||
          N_face[i].squaredNorm() < 1e-24)
        continue;
      const vec3 a_hat = a.normalized();
      const vec3 n_hat = N_face[i].normalized();
      // Foot on the infinite axis: c_a = c + ((p−c)·a) a
      const vec3 p = x_face[i];
      const vec3 c_a = c + a_hat * (p - c).dot(a_hat);
      const vec3 radial = p - c_a;
      const real rho = radial.norm();
      if (rho < 1e-12)
        continue;
      const vec3 rho_hat = radial / rho;
      const vec3 step = (_trunk_weight * (rho - r)) * rho_hat;
      f_face[i] = step.dot(n_hat) * n_hat;
    }
    return f_face;
  }

  /// Soft-slice TP at vertices: Cauchy core (same as braid_planar_demo) in the
  /// null plane of the world medial axis; unit-peak Gaussian on axial sep.
  /// Descent force (−weight · ∇E).
  std::vector<vec3> axis_slice_tp_vert_forces(
      const std::vector<vec3> &x, const medial_face_field &med,
      int frame) const {
    std::vector<vec3> f(x.size(), vec3::Zero());
    if (_slice_tp_weight <= 0.0)
      return f;

    // Area-weight only valid medial faces so we do not invent axes from zeros.
    std::vector<vec3> A_acc(__M->vert_count(), vec3::Zero());
    std::vector<real> r_acc(__M->vert_count(), 0.0);
    std::vector<real> w_acc(__M->vert_count(), 0.0);
    const std::vector<real> areas = asawa::shell::face_areas(*__M, x);
    for (auto fi : __M->get_face_range()) {
      const size_t fii = static_cast<size_t>(static_cast<int>(fi));
      if (fii >= med.valid.size() || !med.valid[fii])
        continue;
      if (fii >= med.axis.size() || med.axis[fii].squaredNorm() < 1e-24)
        continue;
      const real aw =
          fii < areas.size() ? std::max(areas[fii], real(0.0)) : real(1.0);
      if (!(aw > 0.0))
        continue;
      const vec3 a = med.axis[fii].normalized();
      const real r =
          fii < med.radius.size() && std::isfinite(med.radius[fii])
              ? med.radius[fii]
              : real(0.0);
      __M->const_for_each_face(
          asawa::shell::face_id(fi),
          [&](asawa::shell::CornerId c0, const asawa::shell::shell &M) {
            const size_t vid =
                static_cast<size_t>(static_cast<int>(M.vert(c0)));
            if (vid >= A_acc.size())
              return;
            A_acc[vid] += aw * a;
            r_acc[vid] += aw * r;
            w_acc[vid] += aw;
          });
    }

    std::vector<vec3> A_vert(x.size(), vec3::Zero());
    std::vector<real> r_vert(x.size(), 0.0);
    int n_axis = 0;
    for (size_t i = 0; i < A_vert.size(); ++i) {
      if (i >= w_acc.size() || !(w_acc[i] > 0.0))
        continue;
      vec3 a = A_acc[i] / w_acc[i];
      if (a.squaredNorm() < 1e-24)
        continue;
      A_vert[i] = a.normalized();
      r_vert[i] = r_acc[i] / w_acc[i];
      ++n_axis;
    }
    if (n_axis == 0)
      return f;

    const real avg_len = asawa::shell::avg_length(*__M, x);
    const real l0 = std::max(_slice_tp_l0_scale * avg_len, real(1e-12));
    const real sigma_floor =
        std::max(_slice_tp_sigma_par_scale * avg_len, real(1e-12));

    std::vector<real> sigma_par(x.size(), sigma_floor);
    for (size_t i = 0; i < sigma_par.size(); ++i) {
      if (A_vert[i].squaredNorm() < 1e-24)
        continue;
      if (std::isfinite(r_vert[i]) && r_vert[i] > 1e-12)
        sigma_par[i] =
            std::max(_slice_tp_sigma_par_scale * r_vert[i], sigma_floor);
    }

    std::vector<vec3> Nv = asawa::shell::vertex_normals(*__M, x);
    std::vector<real> w = asawa::shell::vertex_areas(*__M, x);
    std::vector<vec3> g = calder::tangent_point_gradient_axis_slice(
        *__M, x, w, Nv, A_vert, sigma_par, l0, _slice_tp_p,
        /*project_force_to_plane=*/true);

    real g_max = 0.0;
    real f_max = 0.0;
    int n_bad = 0;
    const size_t n = std::min(f.size(), g.size());
    for (size_t i = 0; i < n; ++i) {
      if (A_vert[i].squaredNorm() < 1e-24)
        continue;
      if (!g[i].allFinite()) {
        ++n_bad;
        continue;
      }
      g_max = std::max(g_max, g[i].norm());
      vec3 fi = -_slice_tp_weight * g[i];
      if (_slice_tp_force_clamp > 0.0 && fi.norm() > _slice_tp_force_clamp)
        fi *= _slice_tp_force_clamp / fi.norm();
      f[i] = fi;
      f_max = std::max(f_max, fi.norm());
    }
    if (should_trace_frame(frame)) {
      std::cerr << "[growth_study_core] slice_tp: n_axis=" << n_axis
                << " g_max=" << g_max << " f_max=" << f_max
                << " n_bad=" << n_bad << " l0=" << l0 << std::endl;
    }
    return f;
  }

  void step_dynamics(int frame) {
    if (should_trace_frame(frame)) {
      std::cerr << "[growth_study_core] step_dynamics begin frame=" << frame
                << std::endl;
    }
    hepworth::block::projection_solver solver;
    std::vector<hepworth::projection_constraint::ptr> constraints;

    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    std::vector<vec3> &v = asawa::get_vec_data(*__M, __surf->__vdatum_id);

    if (_nan_probe && (frame % _nan_probe_stride == 0)) {
      nan_probe_vec3("x(raw)", x);
      nan_probe_vec3("v(raw)", v);
    }

    for (index_t ii = 0; ii < static_cast<index_t>(x.size()); ++ii) {
      if (!x[ii].allFinite())
        x[ii] = vec3::Zero();
      if (!v[ii].allFinite())
        v[ii] = vec3::Zero();
    }

    if (!_enable_topology && !__M->verts_are_dense_packed()) {
      std::cerr << "[growth_study_core] WARNING: verts not dense-packed while "
                   "topology is off — absolute id layout is holey"
                << std::endl;
    }

    std::vector<vec3> M_area = asawa::shell::vertex_areas_3(*__M, x);
    const real m_min = 1e-14;
    for (auto &mi : M_area) {
      if (!mi.allFinite() || mi[0] < m_min)
        mi = vec3(m_min, m_min, m_min);
    }

    index_t imin = get_new_origin(*__M);
    std::vector<real> fh(__M->vert_count(), 0.0);
    fh[imin] = 1.0;
    kusama::laplacian L(__M, x);
    const std::vector<real> d_vert = L.heatDist(fh, 0.2);
    const std::vector<vec3> x_face = asawa::shell::face_centers(*__M, x);
    const std::vector<vec3> N_face = asawa::shell::face_normals(*__M, x);
    const std::vector<real> d_face =
        asawa::shell::vert_to_face<real>(*__M, x, d_vert);
    const std::vector<vec3> heat_flow = heat_flow_face_dirs(L, fh, x);

    pierce_cylinder_medial_stats mstats;
    medial_face_field medial = compute_faired_medial_face_field(x, &mstats);
    orient_axials_to_heat_flow(medial.axial, heat_flow);
    const std::vector<vec3> &axial_face = medial.axial;
    if (should_trace_frame(frame)) {
      std::cerr << "[growth_study_core] medial(face): accepted="
                << mstats.accepted << "/" << mstats.total
                << " fair_l0_scale=" << _axial_fair_l0_scale
                << " cyl_fair_r_scale=" << _cylinder_fair_radius_scale
                << " mls_sigma_scale=" << _cylinder_mls_sigma_scale
                << " trunk_w=" << _trunk_weight
                << " slice_tp_w=" << _slice_tp_weight << std::endl;
    }

    std::vector<vec3> f_axial =
        axial_face_forces(axial_face, heat_flow, d_face);
    std::vector<vec3> f_trunk =
        trunkification_face_forces(x_face, N_face, medial);
    cache_axial_force_debug(x_face, f_axial, d_face);
    cache_trunk_force_debug(x_face, f_trunk, medial);
    std::vector<vec3> f_face = f_axial;
    {
      const size_t n = std::min(f_face.size(), f_trunk.size());
      for (size_t i = 0; i < n; ++i)
        f_face[i] += f_trunk[i];
    }
    // Absolute-id scatter: f[vert_id] averages incident live face samples.
    // Same layout as x/v/mass for the integrator — never compressed range size.
    std::vector<vec3> f(x.size(), vec3::Zero());
    {
      std::vector<real> w_accum(x.size(), 0.0);
      for (auto vi : __M->get_vert_range()) {
        const size_t vid = static_cast<size_t>(static_cast<int>(vi));
        __M->const_for_each_vertex(
            asawa::shell::vert_id(vi), [&](asawa::shell::CornerId c0,
                                          const asawa::shell::shell &M) {
              const index_t fi = static_cast<index_t>(M.face(c0));
              if (fi < 0 || static_cast<size_t>(fi) >= f_face.size())
                return;
              if (M.fbegin(asawa::shell::face_id(fi)) < 0)
                return;
              f[vid] += f_face[static_cast<size_t>(fi)];
              w_accum[vid] += 1.0;
            });
      }
      for (size_t i = 0; i < f.size(); ++i) {
        if (w_accum[i] > 0.0)
          f[i] /= w_accum[i];
      }
    }

    std::vector<vec3> f_slice_tp =
        axis_slice_tp_vert_forces(x, medial, frame);
    cache_slice_tp_force_debug(x, f_slice_tp);
    {
      const size_t n = std::min(f.size(), f_slice_tp.size());
      for (size_t i = 0; i < n; ++i)
        f[i] += f_slice_tp[i];
    }
    cache_vert_force_debug(x, f);
    hepworth::vec3_block::ptr X = hepworth::vec3_block::create(M_area, x, v, f);
    if (should_trace_frame(frame))
      std::cerr << "[growth_study_core] sim block created" << std::endl;

    // Anisotropic triangle strain: stretch = 1 + heat * _growth along axial/hoop.
    {
      auto face_range = __M->get_face_range();
      for (auto fi : face_range) {
        asawa::shell::FaceId fid = asawa::shell::face_id(fi);
        if (__M->fsize(fid) != 3)
          continue;
        if (asawa::shell::face_area(*__M, fid, x) < 1e-8)
          continue;
        if (fi >= static_cast<int>(axial_face.size()))
          continue;
        const vec3 axial = axial_face[fi];
        if (!axial.allFinite() || axial.squaredNorm() < 1e-24)
          continue;
        const real heat = fi < static_cast<int>(d_face.size()) &&
                                  std::isfinite(d_face[fi])
                              ? std::clamp(d_face[fi], real(0.0), real(1.0))
                              : 0.0;
        const vec2 stretch(1.0 + heat * _growth.x(),
                           1.0 + heat * _growth.y());
        auto tri = __M->get_tri(fid);
        constraints.push_back(
            hepworth::block::anisotropic_triangle_strain::create(
                std::vector<index_t>({tri[0], tri[1], tri[2]}), x, axial,
                stretch, _strain_weight,
                std::vector<hepworth::sim_block::ptr>{X}));
      }
    }

    const real avg_len = asawa::shell::avg_length(*__M, x);

    // Usual one-ring cotan bending (knotted_surface default scale).
    if (_bending_weight > 0.0) {
      hepworth::block::init_bending(*__M, constraints, x, _bending_weight, {X});
    }

    // Willmore fairing / smoothing (knotted_surface uses ~1).
    if (_willmore_weight > 0.0) {
      hepworth::block::init_willmore(*__M, constraints, x, _willmore_weight,
                                     {X});
    }

    // Slab pierce feet: relative plane twist and/or RBF reduced elastic element.
    if (_twist_weight > 0.0 || _rbf_weight > 0.0) {
      hepworth::block::slab_pierce_config pierce_cfg;
      hepworth::block::slab_prism_config prism_cfg;
      prism_cfg.w_twist = _twist_weight;
      prism_cfg.w_rbf = _rbf_weight;
      prism_cfg.w_field = 0.0;
      prism_cfg.w_phi = 0.0;
      prism_cfg.w_bulk = 0.0;
      prism_cfg.rbf_mu = _rbf_mu;
      prism_cfg.rbf_eta = _rbf_eta;
      prism_cfg.rbf_nr_iters = _rbf_nr_iters;
      prism_cfg.rbf_n_centers = _rbf_n_centers;

      arp::bvh_tree<3>::ptr face_tree =
          hepworth::block::make_shell_face_bvh(*__M, x);
      const std::vector<index_t> face_map = __M->get_face_map(true);
      size_t n_twist = 0;
      size_t n_rbf = 0;
      if (face_tree) {
        for (int iv = 0; iv < static_cast<int>(__M->vert_count()); ++iv) {
          const asawa::shell::VertId vid = asawa::shell::vert_id(iv);
          const real A = asawa::shell::vert_area(*__M, vid, x);
          if (A < 1e-8)
            continue;
          if (__M->get_one_ring(vid).empty())
            continue;
          auto pierce = hepworth::block::build_slab_pierce_stencil(
              *__M, vid, x, *face_tree, face_map, pierce_cfg, avg_len);
          if (!pierce.stencil.has_slab)
            continue;
          if (_twist_weight > 0.0) {
            constraints.push_back(hepworth::block::relative_plane_twist::create(
                pierce.stencil, prism_cfg, prism_cfg.w_twist,
                std::vector<hepworth::sim_block::ptr>{X}));
            ++n_twist;
          }
          if (_rbf_weight > 0.0) {
            constraints.push_back(hepworth::block::rbf_elastic_strain::create(
                pierce.stencil, prism_cfg, prism_cfg.w_rbf,
                std::vector<hepworth::sim_block::ptr>{X}, avg_len));
            ++n_rbf;
          }
        }
      }
      if (should_trace_frame(frame)) {
        std::cerr << "[growth_study_core] slab pierce twist=" << n_twist
                  << " rbf=" << n_rbf << std::endl;
      }
    }

    // Self-collisions: point–triangle only. Edge–edge is off until
    // one-ring / shared-vertex filtering is reliable.
    {
      const size_t n0 = constraints.size();
      const real tri_width = std::max(avg_len, real(1e-12));
      const real r = _collision_query_scale * tri_width;
      const real barrier_eps = _collision_barrier_scale * tri_width;
      hepworth::block::init_pnt_tri_collisions(
          *__M, *__surf, constraints, x, r, barrier_eps, _collision_weight,
          {X, X}, _collision_hardness);
      const size_t n_pt = constraints.size() - n0;
      if (should_trace_frame(frame)) {
        std::cerr << "[growth_study_core] collisions pnt_tri=" << n_pt
                  << " r=" << r << " tri_width=" << tri_width
                  << " barrier_eps=" << barrier_eps
                  << " hardness=" << _collision_hardness << std::endl;
      }
    }

    {
      hepworth::block::solver_context ctx;
      ctx.shell_mesh = __M.get();
      hepworth::block::init_block_constraints_rest(constraints, {X}, ctx);
    }

    if (should_trace_frame(frame)) {
      std::cerr << "[growth_study_core] constraints ready count="
                << constraints.size() << std::endl;
    }

    solver.set_constraints(constraints);
    std::vector<hepworth::sim_block::ptr> blocks = {X};
    solver.step(blocks, _h, 0.5, 30);
    if (should_trace_frame(frame))
      std::cerr << "[growth_study_core] step_dynamics end" << std::endl;
  }

  void step(int frame) {
    if (should_trace_frame(frame)) {
      std::cerr << "[growth_study_core] step begin frame=" << frame
                << " topology=" << (_enable_topology ? "on" : "off")
                << std::endl;
    }
    step_rx(frame);
    step_dynamics(frame);
    if (_enable_topology)
      __surf->step(false);
    if (should_trace_frame(frame)) {
      std::cerr << "[growth_study_core] shell dynamic step end" << std::endl;
    }
  }

  void draw_debug() const {
    draw_axial_forces();
    draw_trunk_forces();
    draw_slice_tp_forces();
    draw_vert_forces();
  }

  /// Additive axial/hoop rates; stretch = 1 + heat * _growth.
  vec2 _growth = vec2(0.0, 0.01);
  /// Cylinder MLS Gaussian σ = scale * (pierce_travel / 2). Raise to widen
  /// the neighborhood (e.g. 2–4); 1 matches the historical half-pierce kernel.
  real _cylinder_mls_sigma_scale = 4.0;
  ///  Find push the vertices in the direction of the medial
  /// axis calculated at that point, a volumetric "growth" force
  real _axial_force_weight = 0.0;
  /// Trunkification: radial step projected onto face normal. 0 disables
  /// (negative weight is ignored — flip is in the residual, not the scale).
  real _trunk_weight = 0.0;
  /// Soft-slice TP descent weight (0 disables). Force = −weight · ∇E_slice.
  real _slice_tp_weight = 1.0e-4;
  /// Cauchy ε = l0_scale × mean edge (braid cauchy.l0 ~ 0.5 in rod units).
  real _slice_tp_l0_scale = 0.5;
  /// Cauchy power p (braid cauchy.p = 6).
  real _slice_tp_p = 2.0;
  /// Axial Gaussian σ: σ = scale × max(r_medial, edge).
  /// Window κ = exp(−(axial/σ)²) (unit peak).
  real _slice_tp_sigma_par_scale = 0.5;
  /// Per-vertex clamp on |f_slice| (0 = no clamp).
  real _slice_tp_force_clamp = 1.0;
  /// Axis SVD fairing bandwidth (× mean edge) on the surface.
  real _axial_fair_l0_scale = 4.0;
  /// Cylinder (C,R) fairing bandwidth in center space: σᵢ = scale · rᵢ.
  real _cylinder_fair_radius_scale = 8.0;

  /// When false, skip remesh/pack so we can isolate fairing vs topology noise.
  bool _enable_topology = false;

  /// Debug: draw face-center arrows for axial forces (heat × axial).
  bool _draw_axial_forces = false;
  /// Arrow length scale: visual = _axial_force_draw_scale * f_face.
  real _axial_force_draw_scale = 0.15;
  /// Debug: draw face-center arrows for trunkification forces.
  bool _draw_trunk_forces = false;
  /// Arrow length scale for trunk force debug.
  real _trunk_force_draw_scale = 0.15;
  /// Debug: draw vertex arrows for soft-slice TP forces.
  bool _draw_slice_tp_forces = false;
  real _slice_tp_force_draw_scale = 0.15;
  /// Debug: draw the exact vertex forces passed to the integrator (scatter check).
  bool _draw_vert_forces = false;
  real _vert_force_draw_scale = 0.15;
  real _strain_weight = 1.0e-2;
  real _twist_weight = 0.1;
  /// Local RBF reduced elastic (gauge + strain → ARAP R). 0 disables.
  /// Compare against _twist_weight to isolate relative-plane twist bias.
  real _rbf_weight = 0.0;
  real _rbf_mu = 1.0;
  real _rbf_eta = 1.0;
  int _rbf_nr_iters = 8;
  int _rbf_n_centers = 4;
  /// Usual shell bending; set 0 to disable. knotted_surface default is 5e-2.
  real _bending_weight = 0.0;
  /// Willmore smoothing; set 0 to disable. knotted_surface default is 1.
  real _willmore_weight = 1.0;
  real _collision_weight = 1.0;
  /// Amplifies barrier / deep-penetration pull-out in pnt_tri_collision.
  /// 1 = legacy soft; larger (e.g. 4–8) resists material growth better.
  real _collision_hardness = 1.0;
  /// Collision query radius = _collision_query_scale * mean edge length.
  real _collision_query_scale = 4.0;
  /// Barrier thickness = _collision_barrier_scale * mean edge length.
  real _collision_barrier_scale = 1.0;

  module_base::ptr _rx0;
  growth_rx_model _rx_model;
  std::optional<std::function<vec4(real)>> _mesh_color_gradient;

  vec3 _origin;
  real _h = 0.1;
  real _eps = 0.5;
  std::vector<real> _f;
  std::vector<real> _k;
  std::vector<real> _eps_sh, _g_sh, _lam_sh;
  std::vector<real> _alpha_gl, _beta_gl;

  kusama::laplacian::sparmat _C_aniso;
  real _aniso_vd_lambda = 3.0;
  real _aniso_sigma_u = 1.0;
  real _aniso_sigma_v = 0.001;
  kusama::anisotropic_laplacian_kind _aniso_kind =
      kusama::anisotropic_laplacian_kind::fem_d;

  shell::shell::ptr __M;
  shell::dynamic::ptr __surf;

private:
  bool _nan_probe = false;
  int _nan_probe_stride = 1;

  /// Last-frame axial force debug (face centers, forces, heat weight scalar).
  std::vector<vec3> _dbg_axial_p;
  std::vector<vec3> _dbg_axial_f;
  std::vector<real> _dbg_axial_w;

  /// Last-frame trunk force debug (face centers, forces, signed residual ρ−r).
  std::vector<vec3> _dbg_trunk_p;
  std::vector<vec3> _dbg_trunk_f;
  std::vector<real> _dbg_trunk_res;

  /// Soft-slice TP vertex forces.
  std::vector<vec3> _dbg_slice_tp_p;
  std::vector<vec3> _dbg_slice_tp_f;

  /// Integrator vertex forces (post face→vert scatter) for index diagnostics.
  std::vector<vec3> _dbg_vert_p;
  std::vector<vec3> _dbg_vert_f;

  void cache_vert_force_debug(const std::vector<vec3> &x,
                              const std::vector<vec3> &f) {
    const size_t n = std::min(x.size(), f.size());
    _dbg_vert_p.assign(x.begin(), x.begin() + static_cast<std::ptrdiff_t>(n));
    _dbg_vert_f.assign(f.begin(), f.begin() + static_cast<std::ptrdiff_t>(n));
  }

  void cache_axial_force_debug(const std::vector<vec3> &x_face,
                               const std::vector<vec3> &f_face,
                               const std::vector<real> &d_face) {
    const size_t n = std::min(x_face.size(), f_face.size());
    _dbg_axial_p = x_face;
    _dbg_axial_f = f_face;
    _dbg_axial_p.resize(n);
    _dbg_axial_f.resize(n);
    _dbg_axial_w.resize(n, 0.0);
    for (size_t i = 0; i < n; ++i) {
      const real heat =
          i < d_face.size() && std::isfinite(d_face[i]) ? d_face[i] : 0.0;
      _dbg_axial_w[i] = _axial_force_weight * heat;
    }
  }

  void cache_trunk_force_debug(const std::vector<vec3> &x_face,
                               const std::vector<vec3> &f_trunk,
                               const medial_face_field &med) {
    const size_t n = std::min(x_face.size(), f_trunk.size());
    _dbg_trunk_p = x_face;
    _dbg_trunk_f = f_trunk;
    _dbg_trunk_p.resize(n);
    _dbg_trunk_f.resize(n);
    _dbg_trunk_res.resize(n, 0.0);
    for (size_t i = 0; i < n; ++i) {
      if (i >= med.valid.size() || !med.valid[i]) {
        _dbg_trunk_res[i] = 0.0;
        continue;
      }
      const vec3 &c = med.center[i];
      const vec3 &a = med.axis[i];
      const real r = med.radius[i];
      if (!c.allFinite() || !a.allFinite() || a.squaredNorm() < 1e-24 ||
          !std::isfinite(r) || r <= 1e-8) {
        _dbg_trunk_res[i] = 0.0;
        continue;
      }
      const vec3 a_hat = a.normalized();
      const vec3 c_a = c + a_hat * (x_face[i] - c).dot(a_hat);
      const real rho = (x_face[i] - c_a).norm();
      _dbg_trunk_res[i] = rho - r; // >0 inflate, <0 contract
    }
  }

  void cache_slice_tp_force_debug(const std::vector<vec3> &x,
                                  const std::vector<vec3> &f) {
    const size_t n = std::min(x.size(), f.size());
    _dbg_slice_tp_p.assign(x.begin(),
                           x.begin() + static_cast<std::ptrdiff_t>(n));
    _dbg_slice_tp_f.assign(f.begin(),
                           f.begin() + static_cast<std::ptrdiff_t>(n));
  }

  void draw_axial_forces() const {
    if (!_draw_axial_forces || _dbg_axial_p.empty())
      return;
    const size_t n =
        std::min(_dbg_axial_p.size(),
                 std::min(_dbg_axial_f.size(), _dbg_axial_w.size()));
    real w_max = 1e-12;
    for (size_t i = 0; i < n; ++i) {
      if (std::isfinite(_dbg_axial_w[i]))
        w_max = std::max(w_max, std::abs(_dbg_axial_w[i]));
    }
    const real arrow_r = 0.003;
    for (size_t i = 0; i < n; ++i) {
      const vec3 &f = _dbg_axial_f[i];
      if (!f.allFinite() || f.squaredNorm() < 1e-24)
        continue;
      const real t =
          std::clamp(std::abs(_dbg_axial_w[i]) / w_max, real(0.0), real(1.0));
      const vec4 col(0.15 + 0.85 * t, 0.35 * (1.0 - t), 1.0 - 0.7 * t, 1.0);
      geometry_logger::line(_dbg_axial_p[i],
                            _dbg_axial_p[i] + _axial_force_draw_scale * f, col,
                            arrow_r);
    }
  }

  void draw_trunk_forces() const {
    if (!_draw_trunk_forces || _dbg_trunk_p.empty())
      return;
    const size_t n =
        std::min(_dbg_trunk_p.size(),
                 std::min(_dbg_trunk_f.size(), _dbg_trunk_res.size()));
    real res_max = 1e-12;
    for (size_t i = 0; i < n; ++i) {
      if (std::isfinite(_dbg_trunk_res[i]))
        res_max = std::max(res_max, std::abs(_dbg_trunk_res[i]));
    }
    const real arrow_r = 0.003;
    for (size_t i = 0; i < n; ++i) {
      const vec3 &f = _dbg_trunk_f[i];
      if (!f.allFinite() || f.squaredNorm() < 1e-24)
        continue;
      const real t =
          std::clamp(std::abs(_dbg_trunk_res[i]) / res_max, real(0.0), real(1.0));
      const bool inflate = _dbg_trunk_res[i] >= 0.0;
      const vec4 col =
          inflate ? vec4(0.1, 0.55 + 0.45 * t, 0.2, 1.0)
                  : vec4(1.0, 0.35 * (1.0 - t), 0.05, 1.0);
      geometry_logger::line(_dbg_trunk_p[i],
                            _dbg_trunk_p[i] + _trunk_force_draw_scale * f, col,
                            arrow_r);
    }
  }

  void draw_slice_tp_forces() const {
    if (!_draw_slice_tp_forces || _dbg_slice_tp_p.empty())
      return;
    const size_t n = std::min(_dbg_slice_tp_p.size(), _dbg_slice_tp_f.size());
    real f_max = 1e-12;
    for (size_t i = 0; i < n; ++i) {
      if (_dbg_slice_tp_f[i].allFinite())
        f_max = std::max(f_max, _dbg_slice_tp_f[i].norm());
    }
    const real arrow_r = 0.0025;
    for (size_t i = 0; i < n; ++i) {
      const vec3 &fv = _dbg_slice_tp_f[i];
      if (!fv.allFinite() || fv.squaredNorm() < 1e-24)
        continue;
      const real t = std::clamp(fv.norm() / f_max, real(0.0), real(1.0));
      // Cyan: soft-slice TP descent in Π ⊥ a.
      const vec4 col(0.05, 0.75 + 0.25 * t, 0.95, 1.0);
      geometry_logger::line(
          _dbg_slice_tp_p[i],
          _dbg_slice_tp_p[i] + _slice_tp_force_draw_scale * fv, col, arrow_r);
    }
  }

  void draw_vert_forces() const {
    if (!_draw_vert_forces || _dbg_vert_p.empty())
      return;
    const size_t n = std::min(_dbg_vert_p.size(), _dbg_vert_f.size());
    real f_max = 1e-12;
    for (size_t i = 0; i < n; ++i) {
      if (_dbg_vert_f[i].allFinite())
        f_max = std::max(f_max, _dbg_vert_f[i].norm());
    }
    const real arrow_r = 0.0025;
    for (size_t i = 0; i < n; ++i) {
      const vec3 &fv = _dbg_vert_f[i];
      if (!fv.allFinite() || fv.squaredNorm() < 1e-24)
        continue;
      const real t = std::clamp(fv.norm() / f_max, real(0.0), real(1.0));
      // Magenta: integrator input after face→vert scatter.
      const vec4 col(0.75 + 0.25 * t, 0.1, 0.85, 1.0);
      geometry_logger::line(_dbg_vert_p[i],
                            _dbg_vert_p[i] + _vert_force_draw_scale * fv, col,
                            arrow_r);
    }
  }

  void nan_probe_vec3(const char *tag, const std::vector<vec3> &a) const {
    if (!_nan_probe)
      return;
    for (index_t i = 0; i < static_cast<index_t>(a.size()); ++i) {
      if (!a[i].allFinite()) {
        std::cerr << "[growth_study nan_probe] " << tag << " first bad vert "
                  << i << " -> " << a[i].transpose() << std::endl;
        return;
      }
    }
  }
};

} // namespace duchamp
} // namespace gaudi
#endif
