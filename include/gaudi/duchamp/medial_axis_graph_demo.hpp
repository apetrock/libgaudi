#ifndef __GAUDI_DUCHAMP_MEDIAL_AXIS_GRAPH_DEMO__
#define __GAUDI_DUCHAMP_MEDIAL_AXIS_GRAPH_DEMO__

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>

#include "gaudi/albers/darboux_medial_geometry.hpp"
#include "gaudi/albers/osculating_torus.hpp"
#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/duchamp/body_datum.hpp"
#include "gaudi/duchamp/calder_graph_nodes.hpp"
#include "gaudi/duchamp/darboux_cyclide_medial.hpp"
#include "gaudi/duchamp/field_graph_nodes.hpp"
#include "gaudi/duchamp/medial_graph_nodes.hpp"
#include "gaudi/duchamp/medial_result_types.hpp"
#include "gaudi/geometry_logger.hpp"

#include "liblombardi/graph_context.hpp"

namespace gaudi {
namespace duchamp {

enum class medial_viz_mode {
  HessianFrame,
  CylinderFromMesh,
  CylinderFromMedialPoints
};

template <medial_viz_mode VizMode = medial_viz_mode::HessianFrame>
class medial_axis_graph_demo {
public:
  using ptr = std::shared_ptr<medial_axis_graph_demo>;

  static ptr create(cyclide_medial_params params =
                        default_cyclide_medial_demo_params()) {
    return std::make_shared<medial_axis_graph_demo>(params);
  }

  explicit medial_axis_graph_demo(cyclide_medial_params params)
      : _params(params) {
    _mesh = asawa::shell::load_bunny();
    asawa::shell::triangulate(*_mesh);
    normalize_cyclide_medial_demo_mesh(*_mesh);
    configure_scene_frame();
    wire_graph();
    _graph.run();
    summarize_results();
  }

  void step(int frame) {
    _frame = frame;
    draw_medial_viz();
  }

  int frame() const { return _frame; }
  const cyclide_medial_stats &stats() const { return _stats; }
  asawa::shell::shell::ptr shell() const { return _mesh; }

private:
  void wire_graph() {
    const std::vector<vec3> &x = asawa::const_get_vec_data(*_mesh, 0);
    const real avg_len = asawa::shell::avg_length(*_mesh, x);
    _l0 = std::max(_params.l0_scale * avg_len, real(1e-12));
    _max_travel = _params.max_travel_scale * avg_len;

    auto body = _body = _graph.create_node<body_constant_node>(
        make_shell_body(_mesh));
    auto positions = _positions =
        _graph.create_node<position_snapshot_node>();
    auto normals = _normals =
        _graph.create_node<vertex_normals_snapshot_node>();
    auto cyclide_fit = _cyclide_fit =
        _graph.create_node<darboux_cyclide_fit_node>(
            _l0, _params.fit_p, _params.fit_w0, _params.radius_scale,
            _params.stage1_radius, _params.use_mls_bootstrap,
            _params.bootstrap);
    auto cyclide_smooth = _cyclide_smooth =
        _graph.create_node<darboux_cyclide_smooth_node>(_params.smooth);
    auto medial_search = _medial_search =
        _graph.create_node<medial_shape_energy_search_node>(_max_travel);
    auto hessian = _hessian = _graph.create_node<hessian_frame_node>();
    auto cyl_mesh = _cyl_mesh =
        _graph.create_node<cylinder_fit_mesh_node>(_params.normal_l0,
                                                   _params.fit_p);
    auto cyl_points = _cyl_points =
        _graph.create_node<cylinder_fit_points_node>();

    _graph.link(body->output(), positions->body());
    _graph.link(body->output(), normals->body());
    _graph.link(body->output(), cyclide_fit->body());
    _graph.link(positions->output(), cyclide_fit->pov());
    _graph.link(normals->output(), cyclide_fit->n_pov());
    _graph.link(body->output(), cyclide_smooth->body());
    _graph.link(positions->output(), cyclide_smooth->pov());
    _graph.link(cyclide_fit->output(), cyclide_smooth->cyclide_in());
    _graph.link(positions->output(), medial_search->pov());
    if (_params.enable_cyclide_smooth) {
      _graph.link(cyclide_smooth->cyclide_out(), medial_search->cyclide());
      _graph.link(cyclide_smooth->cyclide_out(), hessian->cyclide());
    } else {
      _graph.link(cyclide_fit->output(), medial_search->cyclide());
      _graph.link(cyclide_fit->output(), hessian->cyclide());
    }
    _graph.link(normals->output(), medial_search->n_pov());
    _graph.link(positions->output(), hessian->pov());
    _graph.link(medial_search->output(), hessian->position());
    _graph.link(body->output(), cyl_mesh->body());
    _graph.link(positions->output(), cyl_mesh->pov());
    _graph.link(normals->output(), cyl_mesh->n_pov());
    _graph.link(positions->output(), cyl_points->pov());
    _graph.link(normals->output(), cyl_points->n_pov());
    _graph.link(medial_search->output(), cyl_points->points());
  }

  void summarize_results() {
    const auto &mask =
        _medial_search
            ->get_datum<medial_shape_energy_search_node::MaskPortDef>()
            ->data();
    const auto &points =
        _medial_search
            ->get_datum<medial_shape_energy_search_node::OutputPortDef>()
            ->data();
    const auto &pov =
        _positions->get_datum<position_snapshot_node::OutputPortDef>()->data();
    _stats.total = static_cast<int>(mask.size());
    _stats.accepted = count_medial_accepted(mask);
    _stats.rejected = _stats.total - _stats.accepted;
    real travel_sum = 0.0;
    for (size_t i = 0; i < mask.size(); ++i) {
      if (mask[i] <= real(0.5) || i >= points.size() || i >= pov.size()) {
        continue;
      }
      travel_sum += (points[i] - pov[i]).norm();
    }
    if (_stats.accepted > 0) {
      _stats.avg_travel = travel_sum / real(_stats.accepted);
    }
    std::cerr << "medial graph demo (ShapeEnergy) accepted=" << _stats.accepted
              << " rejected=" << _stats.rejected
              << " avg_travel=" << _stats.avg_travel;
    if (_params.use_mls_bootstrap) {
      std::cerr << " fit=bootstrap/"
                << calder::mls_bootstrap_ablation_name(
                       _params.bootstrap.ablation);
    } else {
      std::cerr << " fit=legacyAdaptiveGaussian";
    }
    std::cerr << std::endl;
    probe_cyclide_jet_degeneracy();
  }

  // Dump Q + foot-jet health where an osculating torus would go insane.
  // A healthy local ball/cyclide has: D(0)≈0, ∇D∥n, both κ same sign / small
  // radii, and torus normal at the foot aligned with the mesh normal.
  void probe_cyclide_jet_degeneracy() {
    const auto &pov =
        _positions->get_datum<position_snapshot_node::OutputPortDef>()->data();
    const auto &n_pov =
        _normals->get_datum<vertex_normals_snapshot_node::OutputPortDef>()
            ->data();
    const auto &Q =
        _params.enable_cyclide_smooth
            ? _cyclide_smooth
                  ->get_datum<darboux_cyclide_smooth_node::CyclideOutPortDef>()
                  ->data()
            : _cyclide_fit->get_datum<darboux_cyclide_fit_node::OutputPortDef>()
                  ->data();
    if (pov.empty() || Q.empty() || n_pov.empty()) {
      return;
    }

    const size_t n = std::min({pov.size(), Q.size(), n_pov.size()});
    const real l0 = std::max(asawa::shell::avg_length(*_mesh, pov), real(1e-12));
    const real scene = std::max(_major_radius, l0);

    struct jet_probe {
      size_t i = 0;
      real score = 0.0;
      real D0 = 0.0;
      real g_norm = 0.0;
      real align = 0.0; // g·n / (|g||n|)
      real k_min = 0.0;
      real k_max = 0.0;
      real R = 0.0;
      real r = 0.0;
      real quad_norm = 0.0;
      real quartic_norm = 0.0;
      real linear_LX_norm = 0.0;
      real W_frob = 0.0;
      int torus_sign = 0;
      bool hyperbolic = false;
      bool torus_valid = false;
      albers::vec14 Q = albers::vec14::Zero();
      vec3 g = vec3::Zero();
      vec3 n = vec3::Zero();
    };

    auto fill = [&](size_t i) -> jet_probe {
      jet_probe p;
      p.i = i;
      p.Q = Q[i];
      p.n = n_pov[i];
      if (p.n.squaredNorm() > 1e-20) {
        p.n.normalize();
      }
      p.D0 = albers::eval_darboux(p.Q, vec3::Zero());
      p.g = albers::darboux_grad(p.Q, vec3::Zero());
      p.g_norm = p.g.norm();
      if (p.g_norm > 1e-14 && p.n.squaredNorm() > 1e-20) {
        p.align = p.g.dot(p.n) / p.g_norm;
      }
      const mat3 W = albers::shape_operator_at(p.Q, vec3::Zero());
      p.W_frob = W.norm();
      const albers::principal_curvature_frame pc =
          albers::principal_frame_from_shape_operator(W, p.n);
      if (pc.valid) {
        p.k_min = pc.k_min;
        p.k_max = pc.k_max;
        p.hyperbolic = (pc.k_min * pc.k_max < 0.0);
      }
      // Coefficient block norms: quad (A..J), L·X (μ,ν,κ), quartic λ.
      p.quad_norm = p.Q.template head<10>().norm();
      p.linear_LX_norm = p.Q.template segment<3>(11).norm();
      p.quartic_norm = std::abs(p.Q[10]);

      const albers::osculating_torus T =
          albers::osculating_torus_from_shape_operator(
              vec3::Zero(), p.n, W, albers::torus_sign_mode::pat_eq22, l0);
      p.torus_valid = T.valid;
      p.R = T.R;
      p.r = T.r;
      p.torus_sign = T.sign;

      // Degeneracy score: bad normal align, tiny |κ|, hyperbolic, huge R,
      // foot off zero-set, quartic dominating the local jet.
      p.score = 0.0;
      p.score += 5.0 * (1.0 - std::abs(p.align)); // want |align|≈1
      p.score += 2.0 * std::abs(p.D0) / std::max(p.g_norm, real(1e-6));
      if (p.hyperbolic) {
        p.score += 8.0;
      }
      const real k_hi = std::max(std::abs(p.k_min), std::abs(p.k_max));
      if (k_hi < 1.0 / (2.0 * scene)) {
        p.score += 6.0; // near-flat → giant radius
      }
      if (p.torus_valid && p.R > scene) {
        p.score += p.R / scene;
      }
      if (p.quad_norm > 1e-12) {
        p.score += 3.0 * p.quartic_norm / p.quad_norm;
      }
      if (p.g_norm < 1e-6) {
        p.score += 10.0;
      }
      return p;
    };

    auto dump = [&](const char *label, const jet_probe &p) {
      const mat3 W = albers::shape_operator_at(p.Q, vec3::Zero());
      const albers::principal_curvature_frame pc =
          albers::principal_frame_from_shape_operator(W, p.n);
      const albers::osculating_torus T_pat =
          albers::osculating_torus_from_shape_operator(
              vec3::Zero(), p.n, W, albers::torus_sign_mode::pat_eq22, l0);
      const albers::osculating_torus T_in =
          albers::osculating_torus_from_shape_operator(
              vec3::Zero(), p.n, W, albers::torus_sign_mode::interior_biased,
              l0);

      auto torus_n_at_foot = [&](const albers::osculating_torus &T) -> vec3 {
        // Outer-equator contact: tube center is foot - sign*r*n_pc;
        // outward tube normal is from tube-center → foot.
        if (!T.valid || !pc.valid) {
          return vec3::Zero();
        }
        const vec3 tube_center =
            /*foot=*/vec3::Zero() -
            static_cast<real>(T.sign) * T.r * pc.n;
        const vec3 n_geom = (vec3::Zero() - tube_center);
        const vec3 n_sdf = albers::torus_sdf_grad(vec3::Zero(), T);
        std::cerr << "    geom_n=" << n_geom.normalized().transpose()
                  << "  sdf_grad=" << n_sdf.normalized().transpose()
                  << "  sdf(foot)=" << albers::torus_sdf(vec3::Zero(), T)
                  << "\n";
        return n_geom.squaredNorm() > 1e-20 ? n_geom.normalized()
                                            : n_sdf.normalized();
      };

      std::cerr << "---- jet probe [" << label << "] i=" << p.i
                << " score=" << p.score << " ----\n";
      std::cerr << "  D(0)=" << p.D0 << " |g|=" << p.g_norm
                << " align(g,n)=" << p.align
                << " angle_deg="
                << (180.0 / M_PI) *
                       std::acos(std::max(real(-1), std::min(real(1), p.align)))
                << "\n";
      std::cerr << "  k_min=" << p.k_min << " k_max=" << p.k_max
                << " hyperbolic=" << (p.hyperbolic ? 1 : 0)
                << " |W|_F=" << p.W_frob << "\n";
      std::cerr << "  |Q_quad|=" << p.quad_norm
                << " |λ|=" << p.quartic_norm
                << " |(μ,ν,κ)|=" << p.linear_LX_norm
                << " |Q|=" << p.Q.norm() << "\n";
      std::cerr << "  Q = [";
      for (int k = 0; k < 14; ++k) {
        if (k)
          std::cerr << ", ";
        std::cerr << p.Q[k];
      }
      std::cerr << "]\n";
      std::cerr << "  labels: A B C D E F G H I J  λ  μ  ν  κ\n";
      std::cerr << "  g=" << p.g.transpose() << "\n";
      std::cerr << "  n_mesh=" << p.n.transpose() << "\n";
      if (pc.valid) {
        std::cerr << "  n_pc=" << pc.n.transpose()
                  << "  n_pc·n_mesh=" << pc.n.dot(p.n) << "\n";
      }
      if (p.g_norm > 1e-14) {
        const vec3 g_hat = p.g / p.g_norm;
        std::cerr << "  g_hat=" << g_hat.transpose()
                  << "  g_hat·n_mesh=" << g_hat.dot(p.n) << "\n";
      }

      std::cerr << "  --- normal triad (pat_eq22, sign=" << T_pat.sign
                << " R=" << T_pat.R << " r=" << T_pat.r << ") ---\n";
      const vec3 n_pat = torus_n_at_foot(T_pat);
      if (p.g_norm > 1e-14) {
        const vec3 g_hat = p.g / p.g_norm;
        std::cerr << "    n_torus·n_mesh=" << n_pat.dot(p.n)
                  << "  n_torus·g_hat=" << n_pat.dot(g_hat)
                  << "  n_torus·n_pc=" << (pc.valid ? n_pat.dot(pc.n) : 0.0)
                  << "\n";
      }

      std::cerr << "  --- normal triad (interior_biased, sign=" << T_in.sign
                << " R=" << T_in.R << " r=" << T_in.r << ") ---\n";
      const vec3 n_in = torus_n_at_foot(T_in);
      if (p.g_norm > 1e-14) {
        const vec3 g_hat = p.g / p.g_norm;
        std::cerr << "    n_torus·n_mesh=" << n_in.dot(p.n)
                  << "  n_torus·g_hat=" << n_in.dot(g_hat)
                  << "  n_torus·n_pc=" << (pc.valid ? n_in.dot(pc.n) : 0.0)
                  << "\n";
      }
    };

    // Population stats.
    int n_hyp = 0, n_flat = 0, n_misalign = 0, n_huge_R = 0, n_bad_D0 = 0;
    int n_torus_flip_pat = 0, n_torus_flip_in = 0, n_g_flip = 0;
    real align_sum = 0.0;
    jet_probe worst{};
    jet_probe best_ball{};
    bool have_ball = false;
    worst.score = -1.0;

    for (size_t i = 0; i < n; ++i) {
      const jet_probe p = fill(i);
      align_sum += std::abs(p.align);
      if (p.hyperbolic)
        ++n_hyp;
      if (std::max(std::abs(p.k_min), std::abs(p.k_max)) < 1.0 / (2.0 * scene))
        ++n_flat;
      if (std::abs(p.align) < 0.85)
        ++n_misalign;
      if (p.torus_valid && p.R > scene)
        ++n_huge_R;
      if (std::abs(p.D0) > 1e-3 * std::max(p.g_norm, real(1e-6)))
        ++n_bad_D0;
      if (p.align < 0.0)
        ++n_g_flip;

      const mat3 W = albers::shape_operator_at(p.Q, vec3::Zero());
      const albers::principal_curvature_frame pc =
          albers::principal_frame_from_shape_operator(W, p.n);
      if (pc.valid) {
        const albers::osculating_torus T_pat =
            albers::osculating_torus_from_shape_operator(
                vec3::Zero(), p.n, W, albers::torus_sign_mode::pat_eq22, l0);
        const albers::osculating_torus T_in =
            albers::osculating_torus_from_shape_operator(
                vec3::Zero(), p.n, W, albers::torus_sign_mode::interior_biased,
                l0);
        // Outward torus normal at foot = sign * n_pc (see torus_n_at_foot).
        if (T_pat.valid && static_cast<real>(T_pat.sign) * pc.n.dot(p.n) < 0.0)
          ++n_torus_flip_pat;
        if (T_in.valid && static_cast<real>(T_in.sign) * pc.n.dot(p.n) < 0.0)
          ++n_torus_flip_in;
      }

      if (p.score > worst.score) {
        worst = p;
      }
      // "Healthy ball-like": elliptic, |align| high, moderate κ, small R.
      if (!p.hyperbolic && std::abs(p.align) > 0.95 && p.torus_valid &&
          p.R + p.r < 0.75 * scene && p.r > 0.25 * l0) {
        const real ball_score =
            std::abs(p.k_max - p.k_min) /
            std::max(std::abs(p.k_max), real(1e-9));
        if (!have_ball || ball_score < best_ball.score) {
          best_ball = p;
          best_ball.score = ball_score;
          have_ball = true;
        }
      }
    }

    std::cerr << "cyclide jet degeneracy census (n=" << n << "):\n"
              << "  mean|align(g,n)|=" << (align_sum / real(std::max(n, size_t(1))))
              << "  misalign(|c|<0.85)=" << n_misalign
              << "  g_opposite_n=" << n_g_flip
              << "  hyperbolic=" << n_hyp << "  near_flat=" << n_flat
              << "  huge_R(>scene)=" << n_huge_R << "  bad_D0=" << n_bad_D0
              << "\n"
              << "  torus_n opposite mesh_n: pat_eq22=" << n_torus_flip_pat
              << "  interior_biased=" << n_torus_flip_in << "\n";
    dump("worst", worst);
    if (have_ball) {
      dump("healthy_ball_like", best_ball);
    }

    // Side check: medial vs mesh normal at the worst vert.
    const auto &medial =
        _medial_search
            ->get_datum<medial_shape_energy_search_node::OutputPortDef>()
            ->data();
    if (worst.i < pov.size() && worst.i < medial.size()) {
      const vec3 to_m = medial[worst.i] - pov[worst.i];
      std::cerr << "  side check worst: (medial-foot)·n_mesh="
                << to_m.dot(worst.n) << " |medial-foot|=" << to_m.norm()
                << "  (expect <0 if n is outward)\n";
      std::cerr << "  side check worst: T.center·n_mesh="
                << albers::osculating_torus_from_shape_operator(
                       vec3::Zero(), worst.n,
                       albers::shape_operator_at(worst.Q, vec3::Zero()),
                       albers::torus_sign_mode::pat_eq22, l0)
                       .center.dot(worst.n)
                << "  (sign=+1 → expect <0)\n";
    }

    // A/B: active soft×torus×G vs soft×G (no torus) vs plain harmonic.
    {
      calder::mls_jet_bootstrap_params bp = _params.bootstrap;
      bp.fit_p = _params.fit_p;
      bp.fit_w0 = _params.fit_w0;
      bp.radius_scale = _params.radius_scale;

      bp.ablation = calder::mls_bootstrap_ablation::plain_harmonic;
      const std::vector<albers::vec14> Q_harm =
          calder::darboux_fit_bootstrapped(*_mesh, pov, n_pov, l0, bp);

      bp.ablation = calder::mls_bootstrap_ablation::harmonic_plain_gaussian;
      const std::vector<albers::vec14> Q_softG =
          calder::darboux_fit_bootstrapped(*_mesh, pov, n_pov, l0, bp);

      // Radii stats from stage-1 harmonic (shared prior).
      const std::vector<albers::vec14> Q0 =
          calder::darboux_cyclide_normal_constrained(*_mesh, pov, n_pov, l0,
                                                     _params.fit_p,
                                                     _params.fit_w0);
      const real r_min = 2.0 * std::max(l0, real(1e-12));
      const real r_max = 0.5 * std::max(scene, l0);
      real sig_sum = 0, sig_max = 0, r_raw_sum = 0;
      size_t n_sig_hi = 0, n_clamp_lo = 0, n_clamp_hi = 0;
      const size_t nR = Q0.size();
      for (size_t i = 0; i < nR; ++i) {
        real r_raw = calder::radius_from_cyclide_max_curvature(Q0[i], l0);
        if (!std::isfinite(r_raw) || r_raw < 1e-12) {
          r_raw = l0;
        }
        r_raw_sum += r_raw;
        const real r = std::clamp(r_raw, r_min, r_max);
        if (r_raw < r_min) {
          ++n_clamp_lo;
        }
        if (r_raw > r_max) {
          ++n_clamp_hi;
        }
        const real sig = _params.radius_scale * r;
        sig_sum += sig;
        sig_max = std::max(sig_max, sig);
        if (sig > 0.25 * scene) {
          ++n_sig_hi;
        }
      }

      auto census = [&](const std::vector<albers::vec14> &Qfit, const char *tag) {
        size_t n_mis = 0, n_hyp = 0, n_flat = 0, n_huge = 0, n_bad = 0;
        real align_sum = 0;
        std::vector<real> Rs;
        Rs.reserve(n);
        for (size_t i = 0; i < n; ++i) {
          const albers::vec14 &Qi = Qfit[i];
          vec3 ni = n_pov[i];
          if (ni.squaredNorm() > 1e-20) {
            ni.normalize();
          }
          const real D0 = albers::eval_darboux(Qi, vec3::Zero());
          const vec3 g = albers::darboux_grad(Qi, vec3::Zero());
          const real gn = g.norm();
          const real c =
              (gn > 1e-14 && ni.squaredNorm() > 1e-20) ? g.dot(ni) / gn : 0;
          align_sum += std::abs(c);
          if (std::abs(c) < 0.85) {
            ++n_mis;
          }
          if (std::abs(D0) > 1e-3 * std::max(gn, real(1e-6))) {
            ++n_bad;
          }
          const mat3 W = albers::shape_operator_at(Qi, vec3::Zero());
          const albers::principal_curvature_frame pc =
              albers::principal_frame_from_shape_operator(W, ni);
          real k_min = 0, k_max = 0;
          if (pc.valid) {
            k_min = pc.k_min;
            k_max = pc.k_max;
            if (k_min * k_max < 0.0) {
              ++n_hyp;
            }
          }
          if (std::max(std::abs(k_min), std::abs(k_max)) < 1.0 / (2.0 * scene)) {
            ++n_flat;
          }
          const albers::osculating_torus T =
              albers::osculating_torus_from_shape_operator(
                  vec3::Zero(), ni, W, albers::torus_sign_mode::pat_eq22, l0);
          Rs.push_back(T.valid ? T.R : 1e12);
          if (T.valid && T.R > scene) {
            ++n_huge;
          }
        }
        std::sort(Rs.begin(), Rs.end());
        const real R_med = Rs.empty() ? 0 : Rs[Rs.size() / 2];
        std::cerr << "  [" << tag << "] mean|align|="
                  << (align_sum / real(std::max(n, size_t(1))))
                  << "  misalign=" << n_mis << "  hyp=" << n_hyp
                  << "  flat=" << n_flat << "  huge_R=" << n_huge
                  << "  bad_D0=" << n_bad << "  median_R=" << R_med << "\n";
      };

      std::cerr << "cyclide jet A/B (l0=" << l0 << " scene=" << scene
                << " meanR_raw="
                << (r_raw_sum / real(std::max(nR, size_t(1))))
                << " meanσ="
                << (sig_sum / real(std::max(nR, size_t(1))))
                << " maxσ=" << sig_max << " clamp_lo=" << n_clamp_lo
                << " clamp_hi=" << n_clamp_hi << " σ>0.25scene=" << n_sig_hi
                << "):\n";
      census(Q, "soft×torus×G");
      census(Q_softG, "soft×G(σ=R)");
      census(Q_harm, "harmonic   ");
    }
  }

  void configure_scene_frame() {
    const std::vector<vec3> &x = asawa::const_get_vec_data(*_mesh, 0);
    if (x.empty()) {
      return;
    }
    vec3 lo = x.front();
    vec3 hi = x.front();
    for (const vec3 &p : x) {
      lo = lo.cwiseMin(p);
      hi = hi.cwiseMax(p);
    }
    _center = 0.5 * (lo + hi);
    _major_radius = 0.5 * (hi - lo).norm();
    _minor_radius = 4.0 * asawa::shell::avg_length(*_mesh, x);
  }

  void draw_medial_viz() {
    const vec4 axis_color(0.35, 0.35, 0.35, 1.0);
    geometry_logger::line(_center - 1.7 * _major_radius * vec3::UnitZ(),
                          _center + 1.7 * _major_radius * vec3::UnitZ(),
                          axis_color);

    switch (_params.display) {
    case medial_axis_display::MedialAxis:
      draw_medial_axis_overlay();
      break;
    case medial_axis_display::SmoothedFitHessian:
      draw_smoothed_fit_hessian();
      break;
    }
    draw_debug_osculating_torus();
  }

  void draw_debug_osculating_torus() {
    const auto &pov =
        _positions->get_datum<position_snapshot_node::OutputPortDef>()->data();
    const auto &n_pov =
        _normals->get_datum<vertex_normals_snapshot_node::OutputPortDef>()
            ->data();
    const auto &Q =
        _params.enable_cyclide_smooth
            ? _cyclide_smooth
                  ->get_datum<darboux_cyclide_smooth_node::CyclideOutPortDef>()
                  ->data()
            : _cyclide_fit->get_datum<darboux_cyclide_fit_node::OutputPortDef>()
                  ->data();
    if (pov.empty() || Q.empty() || n_pov.empty()) {
      return;
    }

    const size_t n = std::min({pov.size(), Q.size(), n_pov.size()});
    const real l0 = std::max(asawa::shell::avg_length(*_mesh, pov), real(1e-12));
    const real scene = std::max(_major_radius, l0);

    // Draw the worst degenerate foot-jet (near-flat κ_lo → huge R). Do not
    // filter it out — that is the failure mode under investigation.
    size_t best_i = n;
    real best_score = -1.0;
    albers::osculating_torus best_T{};
    albers::principal_curvature_frame best_pc{};
    for (size_t i = 0; i < n; ++i) {
      const vec3 ni = n_pov[i];
      if (ni.squaredNorm() < 1e-20) {
        continue;
      }
      const mat3 W = albers::shape_operator_at(Q[i], vec3::Zero());
      if (!W.allFinite()) {
        continue;
      }
      const albers::principal_curvature_frame pc =
          albers::principal_frame_from_shape_operator(W, ni);
      if (!pc.valid) {
        continue;
      }
      // Use pat_eq22 so the drawable matches the probe's "worst" torus.
      const albers::osculating_torus T =
          albers::osculating_torus_from_shape_operator(
              vec3::Zero(), ni, W, albers::torus_sign_mode::pat_eq22, l0);
      if (!T.valid) {
        continue;
      }

      const real D0 = albers::eval_darboux(Q[i], vec3::Zero());
      const vec3 g = albers::darboux_grad(Q[i], vec3::Zero());
      const real g_norm = g.norm();
      real align = 0.0;
      if (g_norm > 1e-14) {
        align = g.dot(ni.normalized()) / g_norm;
      }
      const real quad_norm = Q[i].template head<10>().norm();
      const real quartic_norm = std::abs(Q[i][10]);

      real score = 0.0;
      score += 5.0 * (1.0 - std::abs(align));
      score += 2.0 * std::abs(D0) / std::max(g_norm, real(1e-6));
      if (pc.k_min * pc.k_max < 0.0) {
        score += 8.0;
      }
      const real k_hi = std::max(std::abs(pc.k_min), std::abs(pc.k_max));
      if (k_hi < 1.0 / (2.0 * scene)) {
        score += 6.0;
      }
      if (T.R > scene) {
        score += T.R / scene;
      }
      if (quad_norm > 1e-12) {
        score += 3.0 * quartic_norm / quad_norm;
      }
      if (g_norm < 1e-6) {
        score += 10.0;
      }
      // Near-parabolic κ_lo is the failure we care about — weight it hard.
      const real k_lo = std::min(std::abs(pc.k_min), std::abs(pc.k_max));
      if (k_lo < 1e-3 / scene) {
        score += 1.0 / std::max(k_lo, real(1e-12));
      }

      if (score > best_score) {
        best_score = score;
        best_i = i;
        best_T = T;
        best_pc = pc;
      }
    }
    if (best_i >= n || !best_T.valid) {
      return;
    }

    const vec3 foot = pov[best_i];
    const vec3 ni = n_pov[best_i].normalized();
    const vec3 g = albers::darboux_grad(Q[best_i], vec3::Zero());
    const vec3 g_hat = g.norm() > 1e-14 ? g.normalized() : ni;
    const vec3 n_torus = static_cast<real>(best_T.sign) * best_pc.n;
    const real axis_len = 0.25 * scene;
    const vec3 tube_center =
        foot - static_cast<real>(best_T.sign) * best_T.r * best_pc.n;

    geometry_logger::point(foot, vec4(1.0, 0.15, 0.9, 1.0));
    geometry_logger::sphere(tube_center, 0.03 * scene,
                            vec4(1.0, 0.15, 0.9, 0.7));
    geometry_logger::line(foot, foot + axis_len * ni,
                          vec4(0.2, 0.55, 1.0, 0.9)); // mesh n
    geometry_logger::line(foot, foot + axis_len * g_hat,
                          vec4(1.0, 0.85, 0.1, 0.9)); // ∇D
    geometry_logger::line(foot, foot + axis_len * n_torus,
                          vec4(0.25, 1.0, 0.35, 0.9)); // torus n
    geometry_logger::line(foot, tube_center, vec4(1.0, 0.15, 0.9, 0.85));

    // R from κ_lo under the 1e-12 floor is ~1e12 — not representable in f32
    // instance attrs / not on camera. Draw the true torus when finite; otherwise
    // the jet is a cylinder: circle of radius r about tube_center in the
    // (n, e_hi) plane, plus generators along e_lo.
    const bool R_drawable =
        std::isfinite(best_T.R) && best_T.R > 0.0 &&
        best_T.R + best_T.r < real(1e4) * scene;
    if (R_drawable) {
      const vec3 center_w = foot + best_T.center;
      geometry_logger::torus(center_w, best_T.axis, best_T.R, best_T.r,
                             vec4(0.25, 1.0, 0.35, 1.0));
    } else {
      const vec3 e_lo = best_T.axis.normalized();
      const vec3 e_rad = n_torus.normalized();
      const vec3 e_bin = e_lo.cross(e_rad).normalized();
      const int N = 48;
      const real gen = 0.35 * scene;
      for (int k = 0; k < N; ++k) {
        const real a0 = real(2.0 * M_PI) * real(k) / real(N);
        const real a1 = real(2.0 * M_PI) * real(k + 1) / real(N);
        const vec3 p0 =
            tube_center + best_T.r * (std::cos(a0) * e_rad + std::sin(a0) * e_bin);
        const vec3 p1 =
            tube_center + best_T.r * (std::cos(a1) * e_rad + std::sin(a1) * e_bin);
        geometry_logger::line(p0, p1, vec4(0.25, 1.0, 0.35, 1.0));
      }
      geometry_logger::line(tube_center - gen * e_lo, tube_center + gen * e_lo,
                            vec4(0.25, 1.0, 0.35, 0.85));
    }

    // Reference: a healthy near-umbilic jet under the same fit, so there is
    // always a true drawable osculating torus on screen.
    size_t ball_i = n;
    real ball_score = std::numeric_limits<real>::infinity();
    albers::osculating_torus ball_T{};
    for (size_t i = 0; i < n; ++i) {
      const vec3 n = n_pov[i];
      if (n.squaredNorm() < 1e-20) {
        continue;
      }
      const mat3 W = albers::shape_operator_at(Q[i], vec3::Zero());
      if (!W.allFinite()) {
        continue;
      }
      const albers::principal_curvature_frame pc =
          albers::principal_frame_from_shape_operator(W, n);
      if (!pc.valid || pc.k_min * pc.k_max < 0.0) {
        continue;
      }
      const vec3 g_i = albers::darboux_grad(Q[i], vec3::Zero());
      if (g_i.norm() < 1e-14) {
        continue;
      }
      const real align = g_i.normalized().dot(n.normalized());
      if (std::abs(align) < 0.95) {
        continue;
      }
      const albers::osculating_torus T =
          albers::osculating_torus_from_shape_operator(
              vec3::Zero(), n, W, albers::torus_sign_mode::pat_eq22, l0);
      if (!T.valid || T.R + T.r > 0.75 * scene || T.r < 0.25 * l0) {
        continue;
      }
      const real s = std::abs(pc.k_max - pc.k_min) /
                     std::max(std::abs(pc.k_max), real(1e-9));
      if (s < ball_score) {
        ball_score = s;
        ball_i = i;
        ball_T = T;
      }
    }
    if (ball_i < n && ball_T.valid) {
      const vec3 c = pov[ball_i] + ball_T.center;
      geometry_logger::point(pov[ball_i], vec4(0.2, 0.85, 1.0, 1.0));
      geometry_logger::torus(c, ball_T.axis, ball_T.R, ball_T.r,
                             vec4(0.2, 0.85, 1.0, 1.0));
    }

    if (_frame == 0) {
      std::cerr << "debug osculating torus (worst jet): i=" << best_i
                << " score=" << best_score << " R=" << best_T.R
                << " r=" << best_T.r << " sign=" << best_T.sign
                << " k_min=" << best_pc.k_min << " k_max=" << best_pc.k_max
                << " drawable=" << (R_drawable ? 1 : 0)
                << " n_torus·n_mesh=" << n_torus.dot(ni)
                << " g·n_mesh=" << g_hat.dot(ni) << std::endl;
      if (ball_i < n) {
        std::cerr << "  healthy ball torus: i=" << ball_i
                  << " R=" << ball_T.R << " r=" << ball_T.r << std::endl;
      }
    }
  }

  void draw_medial_axis_overlay() {
    const vec4 surface_color(0.0, 0.85, 1.0, 1.0);
    const vec4 medial_color(1.0, 0.65, 0.05, 1.0);
    const vec4 cylinder_color(0.9, 0.15, 1.0, 1.0);
    const vec4 rejected_color(0.35, 0.35, 0.35, 0.45);

    const auto &pov =
        _positions->get_datum<position_snapshot_node::OutputPortDef>()->data();
    const auto &medial =
        _medial_search
            ->get_datum<medial_shape_energy_search_node::OutputPortDef>()
            ->data();
    const auto &mask =
        _medial_search
            ->get_datum<medial_shape_energy_search_node::MaskPortDef>()
            ->data();
    const real axis_len = 8.0 * asawa::shell::avg_length(*_mesh, pov);

    for (size_t i = 0; i < pov.size(); ++i) {
      if (i >= mask.size()) {
        continue;
      }
      const vec3 p = pov[i];
      if (mask[i] <= real(0.5)) {
        geometry_logger::point(p, rejected_color);
        continue;
      }
      const vec3 m = medial[i];
      const real travel = (m - p).norm();
      if (travel > 1e-12) {
        geometry_logger::line(p, m, surface_color);
      } else {
        geometry_logger::point(m, surface_color);
      }
      geometry_logger::point(m, medial_color);

      if constexpr (VizMode == medial_viz_mode::HessianFrame) {
        const auto &frames =
            _hessian->get_datum<hessian_frame_node::OutputPortDef>()->data();
        if (i >= frames.size()) {
          continue;
        }
        const mat3 &F = frames[i];
        const real len = 0.35 * travel;
        for (int ax = 0; ax < 3; ++ax) {
          const vec3 dir = F.col(ax).normalized();
          geometry_logger::line(m - len * dir, m + len * dir, medial_color);
        }
      } else if constexpr (VizMode == medial_viz_mode::CylinderFromMesh) {
        const auto &lines =
            _cyl_mesh->get_datum<cylinder_fit_mesh_node::OutputPortDef>()
                ->data();
        if (i >= lines.size()) {
          continue;
        }
        const vec3 d = albers::plucker_line_direction(lines[i]);
        if (d.squaredNorm() > 1e-12) {
          const vec3 dp = axis_len * d.normalized();
          geometry_logger::line(p - dp, p + dp, cylinder_color);
        }
      } else if constexpr (VizMode ==
                           medial_viz_mode::CylinderFromMedialPoints) {
        const auto &lines =
            _cyl_points->get_datum<cylinder_fit_points_node::OutputPortDef>()
                ->data();
        if (i >= lines.size()) {
          continue;
        }
        const vec3 d = albers::plucker_line_direction(lines[i]);
        if (d.squaredNorm() > 1e-12) {
          const vec3 dp = axis_len * d.normalized();
          geometry_logger::line(m - dp, m + dp, cylinder_color);
        }
      }
    }
  }

  void draw_smoothed_fit_hessian() {
    const vec4 foot_color(0.85, 0.85, 0.85, 0.65);
    const vec4 axis_colors[3] = {vec4(1.0, 0.0, 0.0, 1.0),
                                  vec4(0.0, 1.0, 0.0, 1.0),
                                  vec4(0.0, 0.0, 1.0, 1.0)};

    const auto &pov =
        _positions->get_datum<position_snapshot_node::OutputPortDef>()->data();
    const auto &Q =
        _params.enable_cyclide_smooth
            ? _cyclide_smooth
                  ->get_datum<darboux_cyclide_smooth_node::CyclideOutPortDef>()
                  ->data()
            : _cyclide_fit->get_datum<darboux_cyclide_fit_node::OutputPortDef>()
                  ->data();

    const real len =
        _params.hessian_frame_scale * asawa::shell::avg_length(*_mesh, pov);
    const real line_radius = _params.hessian_line_radius;
    const size_t n = std::min(pov.size(), Q.size());

    for (size_t i = 0; i < n; ++i) {
      const vec3 p = pov[i];
      geometry_logger::point(p, foot_color);

      const mat3 W = albers::shape_operator_at(Q[i], vec3::Zero());
      if (!W.allFinite()) {
        continue;
      }

      Eigen::SelfAdjointEigenSolver<mat3> es(W);
      if (es.info() != Eigen::Success) {
        continue;
      }

      const mat3 &V = es.eigenvectors();
      for (int ax = 0; ax < 3; ++ax) {
        const vec3 dir = V.col(ax);
        if (dir.squaredNorm() < 1e-24) {
          continue;
        }
        const vec3 d = dir.normalized();
        geometry_logger::line(p - len * d, p + len * d, axis_colors[ax],
                              line_radius);
      }
    }
  }

  cyclide_medial_params _params;
  cyclide_medial_stats _stats;
  asawa::shell::shell::ptr _mesh;
  liblombardi::GraphContext _graph;
  body_constant_node::ptr _body;
  position_snapshot_node::ptr _positions;
  vertex_normals_snapshot_node::ptr _normals;
  darboux_cyclide_fit_node::ptr _cyclide_fit;
  darboux_cyclide_smooth_node::ptr _cyclide_smooth;
  medial_shape_energy_search_node::ptr _medial_search;
  hessian_frame_node::ptr _hessian;
  cylinder_fit_mesh_node::ptr _cyl_mesh;
  cylinder_fit_points_node::ptr _cyl_points;
  real _l0 = 1.0;
  real _max_travel = 1.0;
  vec3 _center = vec3::Zero();
  real _major_radius = 1.0;
  real _minor_radius = 0.35;
  int _frame = 0;
};

} // namespace duchamp
} // namespace gaudi

#endif // __GAUDI_DUCHAMP_MEDIAL_AXIS_GRAPH_DEMO__
