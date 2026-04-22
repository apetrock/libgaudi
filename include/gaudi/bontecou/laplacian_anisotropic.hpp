#ifndef GAUDI_BONTECOU_LAPLACIAN_ANISOTROPIC_HPP
#define GAUDI_BONTECOU_LAPLACIAN_ANISOTROPIC_HPP

/// Field-guided anisotropic cotan Laplacian (two discretizations) built from per-edge
/// VD-style guidance \f$g \in \mathbb{R}^{2|E|}\f$. See plan / docs in repo.

#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/bontecou/laplacian.hpp"
#include "gaudi/bontecou/vector_dirichlet.hpp"
#include "gaudi/bontecou/vector_dirichlet_guided.hpp"

#include <Eigen/Sparse>
#include <cmath>
#include <vector>

namespace gaudi {
namespace bontecou {

enum class anisotropic_laplacian_kind { conductance, fem_d };

/// Per-(edge, face) guidance lifted into face \f$f\f$'s tangent plane.
struct edge_face_frame_t {
  vec3 w_alpha;
  vec3 w_alpha_perp;
  bool degenerate;
};

/// \p c_edge: corner whose edge is \f$(\mathrm{vert}(c),\mathrm{vert}(\mathrm{next}(c)))\f$.
/// \p f: face for the frame (usually \c M.face(c_edge)).
/// \p g: VD vector length \c 2*nE, layout [parallel block | perp block].
inline edge_face_frame_t
edge_face_frame(const asawa::shell::shell &M, const std::vector<vec3> &x,
                asawa::shell::CornerId c_edge, asawa::shell::FaceId f,
                const Eigen::VectorXd &g, const std::vector<int> &slot_map,
                int nE) {
  edge_face_frame_t out{vec3::Zero(), vec3::Zero(), true};

  const int slot = static_cast<int>(c_edge) / 2;
  const int dof = slot_map[static_cast<size_t>(slot)];
  if (dof < 0 || dof >= nE)
    return out;

  const vec3 n_f = asawa::shell::face_normal(M, f, x);
  const real gp = g(dof);
  const real gq = g(dof + nE);
  const vec3 W = edge_guidance_vector_3d(M, c_edge, x, n_f, gp, gq);

  const vec3 Pi_W = W - n_f * W.dot(n_f);
  const real n_pw = Pi_W.norm();
  constexpr real kEps = 1e-14;
  if (n_pw < kEps)
    return out;

  out.w_alpha = Pi_W / n_pw;
  out.w_alpha_perp = n_f.cross(out.w_alpha);
  out.degenerate = false;
  return out;
}

/// Kind A: \f$k(e,f)\cdot\f$ weak_cotan. \p c_opp is the corner opposite the edge in \c face(c_opp).
inline real weight_conductance(const asawa::shell::shell &M,
                                const std::vector<vec3> &x,
                                asawa::shell::CornerId c_opp,
                                const Eigen::VectorXd &g,
                                const std::vector<int> &slot_map, int nE, real su,
                                real sv) {
  const asawa::shell::CornerId c_edge = M.next(c_opp);
  const asawa::shell::FaceId f = M.face(c_opp);
  const real cot_g = asawa::shell::weak_cotan(M, c_opp, x);

  const edge_face_frame_t fr = edge_face_frame(M, x, c_edge, f, g, slot_map, nE);
  if (fr.degenerate)
    return static_cast<real>(0.5) * (su + sv) * cot_g;

  const vec3 e_ij =
      x[M.vert(M.prev(c_opp))] - x[M.vert(M.next(c_opp))];
  const real L = e_ij.norm();
  if (L < static_cast<real>(1e-20))
    return static_cast<real>(0.5) * (su + sv) * cot_g;
  const vec3 eh = e_ij / L;

  const real a = eh.dot(fr.w_alpha);
  const real b = eh.dot(fr.w_alpha_perp);
  return (su * a * a + sv * b * b) * cot_g;
}

/// Kind B: P1-style rank-1 sum over the three edges of the face.
inline real weight_fem_d(const asawa::shell::shell &M, const std::vector<vec3> &x,
                          asawa::shell::CornerId c_opp, const Eigen::VectorXd &g,
                          const std::vector<int> &slot_map, int nE, real su,
                          real sv) {
  const asawa::shell::FaceId f = M.face(c_opp);
  const real A_f = asawa::shell::face_area(M, f, x);
  if (A_f < static_cast<real>(1e-20))
    return static_cast<real>(0);

  const asawa::shell::VertId vk = M.vert(c_opp);
  const asawa::shell::VertId vi = M.vert(M.next(c_opp));
  const asawa::shell::VertId vj = M.vert(M.prev(c_opp));
  const vec3 e_i = x[vk] - x[vj];
  const vec3 e_j = x[vi] - x[vk];
  const real e_i_dot_e_j = e_i.dot(e_j);

  std::vector<asawa::shell::CornerId> face_corners;
  face_corners.reserve(8);
  M.const_for_each_face(f, [&](asawa::shell::CornerId cid, const asawa::shell::shell &) {
    face_corners.push_back(cid);
  });
  const int nfc = static_cast<int>(face_corners.size());
  if (nfc < 2)
    return static_cast<real>(0);

  std::vector<real> ell(static_cast<size_t>(nfc));
  real ell_sum = static_cast<real>(0);
  for (int t = 0; t < nfc; ++t) {
    ell[static_cast<size_t>(t)] =
        asawa::shell::edge_length(M, face_corners[static_cast<size_t>(t)], x);
    ell_sum += ell[static_cast<size_t>(t)];
  }
  if (ell_sum < static_cast<real>(1e-20))
    return static_cast<real>(0);

  real bracket_sum = static_cast<real>(0);
  for (int t = 0; t < nfc; ++t) {
    const real W_e = ell[static_cast<size_t>(t)] / ell_sum;
    const edge_face_frame_t fr = edge_face_frame(
        M, x, face_corners[static_cast<size_t>(t)], f, g, slot_map, nE);
    if (fr.degenerate) {
      bracket_sum += W_e * static_cast<real>(0.5) * (su + sv) * e_i_dot_e_j;
      continue;
    }
    const real wi_p = fr.w_alpha_perp.dot(e_i);
    const real wj_p = fr.w_alpha_perp.dot(e_j);
    const real wi_a = fr.w_alpha.dot(e_i);
    const real wj_a = fr.w_alpha.dot(e_j);
    bracket_sum += W_e * (su * wi_p * wj_p + sv * wi_a * wj_a);
  }

  return -bracket_sum / (static_cast<real>(2) * A_f);
}

inline Eigen::SparseMatrix<real>
build_anisotropic_cotan_laplacian(asawa::shell::shell &M,
                                  const std::vector<vec3> &x,
                                  const Eigen::VectorXd &g,
                                  const std::vector<int> &slot_map, int nE,
                                  real sigma_u, real sigma_v,
                                  anisotropic_laplacian_kind kind) {
  if (kind == anisotropic_laplacian_kind::conductance) {
    return build_lap<1>(
        M, x,
        [&](asawa::shell::shell &Ms, asawa::shell::CornerId c,
            const std::vector<vec3> &xs) -> real {
          const asawa::shell::CornerId c0p = Ms.prev(c);
          const asawa::shell::CornerId c1p = Ms.prev(Ms.other(c));
          return weight_conductance(Ms, xs, c0p, g, slot_map, nE, sigma_u,
                                    sigma_v) +
                 weight_conductance(Ms, xs, c1p, g, slot_map, nE, sigma_u,
                                    sigma_v);
        });
  }
  return build_lap<1>(
      M, x,
      [&](asawa::shell::shell &Ms, asawa::shell::CornerId c,
          const std::vector<vec3> &xs) -> real {
        const asawa::shell::CornerId c0p = Ms.prev(c);
        const asawa::shell::CornerId c1p = Ms.prev(Ms.other(c));
        return weight_fem_d(Ms, xs, c0p, g, slot_map, nE, sigma_u, sigma_v) +
               weight_fem_d(Ms, xs, c1p, g, slot_map, nE, sigma_u, sigma_v);
      });
}

inline Eigen::SparseMatrix<real>
build_curvature_aligned_laplacian(asawa::shell::shell &M,
                                  const std::vector<vec3> &x, real vd_lambda,
                                  real sigma_u, real sigma_v,
                                  asawa::shell::face_curvature_stencil stencil,
                                  anisotropic_laplacian_kind kind) {
  Eigen::VectorXd g = solve_curvature_guided_vector_dirichlet(
      M, x, vd_lambda, stencil, /*apply_sign_coherence=*/true);

  std::vector<int> slot_map;
  const int nE = build_compact_edge_dof_map(M, slot_map);
  if (nE <= 0 || g.size() != 2 * nE) {
    return build_lap<1>(
        M, x,
        [](asawa::shell::shell &Ms, asawa::shell::CornerId c,
           const std::vector<vec3> &xs) -> real {
          return asawa::shell::weak_cotan(Ms, Ms.prev(c), xs) +
                 asawa::shell::weak_cotan(Ms, Ms.prev(Ms.other(c)), xs);
        });
  }
  return build_anisotropic_cotan_laplacian(M, x, g, slot_map, nE, sigma_u,
                                           sigma_v, kind);
}

} // namespace bontecou
} // namespace gaudi

#endif
