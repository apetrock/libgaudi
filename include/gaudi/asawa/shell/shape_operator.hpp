#ifndef __ASAWA_SHELL_SHAPE_OPERATOR__
#define __ASAWA_SHELL_SHAPE_OPERATOR__

#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/common.h"

#include <vector>

namespace gaudi {
namespace asawa {
namespace shell {

/// Weingarten / shape operator from a discrete principal frame:
/// \(W = \kappa_{\min} t_{\min}t_{\min}^\top + \kappa_{\max} t_{\max}t_{\max}^\top\).
inline mat3 face_shape_operator(const face_curvature_frame &fc) {
  return fc.k_min * fc.t_min * fc.t_min.transpose() +
         fc.k_max * fc.t_max * fc.t_max.transpose();
}

/// Per-face discrete shape operators (dense face index = FaceId).
inline std::vector<mat3> face_shape_operators(
    const shell &M, const std::vector<vec3> &x,
    face_curvature_stencil stencil = face_curvature_stencil::one_ring) {
  const int nf = static_cast<int>(M.face_count());
  std::vector<mat3> W(static_cast<size_t>(nf), mat3::Zero());
  for (int fi = 0; fi < nf; ++fi) {
    const FaceId f = face_id(fi);
    if (M.fbegin(f) < 0 || M.fsize(f) != 3) {
      continue;
    }
    W[static_cast<size_t>(fi)] =
        face_shape_operator(face_curvature_frame_fit(M, x, f, stencil));
  }
  return W;
}

/// Average incident-face \(W\) at each vertex (for vertex MLS queries).
inline std::vector<mat3> vertex_shape_operators_from_faces(
    const shell &M, const std::vector<mat3> &face_W) {
  const int nv = static_cast<int>(M.vert_count());
  std::vector<mat3> W(static_cast<size_t>(nv), mat3::Zero());
  std::vector<real> counts(static_cast<size_t>(nv), 0.0);
  const int nf = static_cast<int>(face_W.size());
  for (int fi = 0; fi < nf; ++fi) {
    const FaceId f = face_id(fi);
    if (M.fbegin(f) < 0 || M.fsize(f) != 3) {
      continue;
    }
    const mat3 &Wf = face_W[static_cast<size_t>(fi)];
    if (!Wf.allFinite()) {
      continue;
    }
    M.const_for_each_face(f, [&](CornerId c, const shell &Mm) {
      const int vi = static_cast<int>(Mm.vert(c));
      if (vi < 0 || vi >= nv) {
        return;
      }
      W[static_cast<size_t>(vi)] += Wf;
      counts[static_cast<size_t>(vi)] += 1.0;
    });
  }
  for (int vi = 0; vi < nv; ++vi) {
    if (counts[static_cast<size_t>(vi)] > 0.5) {
      W[static_cast<size_t>(vi)] /= counts[static_cast<size_t>(vi)];
    }
  }
  return W;
}

inline std::vector<mat3> vertex_shape_operators(
    const shell &M, const std::vector<vec3> &x,
    face_curvature_stencil stencil = face_curvature_stencil::one_ring) {
  return vertex_shape_operators_from_faces(M, face_shape_operators(M, x, stencil));
}

} // namespace shell
} // namespace asawa
} // namespace gaudi

#endif // __ASAWA_SHELL_SHAPE_OPERATOR__
