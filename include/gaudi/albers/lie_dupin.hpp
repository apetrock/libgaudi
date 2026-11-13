#ifndef GAUDI_ALBERS_LIE_DUPIN_HPP
#define GAUDI_ALBERS_LIE_DUPIN_HPP

#include "gaudi/albers/lie_sphere.hpp"

#include <Eigen/Eigenvalues>
#include <Eigen/SVD>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

/// A Dupin cyclide as one coordinate: the involution sigma of R^{4,2}.
///
/// A Dupin cyclide (torus, cylinder and cone of revolution, and their Möbius
/// and Laguerre images) is the surface whose two families of curvature spheres
/// each lie in a 3-dimensional subspace: R^{4,2} = U + U^perp, both of
/// signature (2,1), one family on the conic of null vectors in U, the other in
/// U^perp (Cecil, *Lie Sphere Geometry*, section 5.4). The cyclide is the
/// involution
///
///     sigma = +1 on U,  -1 on U^perp,
///
/// which is eta-self-adjoint (`eta sigma` symmetric), squares to the identity
/// and has trace 0. Such sigma form the symmetric space
/// O(4,2) / (O(2,1) x O(2,1)) of dimension 15 - 3 - 3 = 9, the classical
/// 9-parameter family of cyclides. -sigma is the same surface with the two
/// families swapped.
///
/// What the coordinate buys:
///
/// - Membership is linear. A contact element L lies on the cyclide exactly
///   when sigma L = L (`lie_dupin_contact_residual`), so fitting sigma to a
///   ring of points with normals is a homogeneous linear least-squares problem
///   over eta-self-adjoint matrices (`lie_dupin_fit`). Four contact elements
///   on a cyclide determine it (three leave a pencil); six or more is well
///   conditioned. Generic contact elements lie on no common cyclide, so on
///   real data the fit is least squares and its residual measures how far
///   the data is from Dupin.
/// - The principal curvatures at a point are read off as the members of L
///   inside U and inside U^perp (`lie_dupin_curvatures`), with no 1/k: a
///   plane is a sphere, and flat or umbilic points are not singular.
/// - A Lie sphere transformation A moves it by conjugation, `A sigma A^-1`
///   (`lie_dupin_transform`): advecting the shape with the surface is one
///   matrix product, and never leaves the manifold.
/// - Projection back onto the manifold is the matrix sign function
///   (`lie_dupin_normalize`), the analogue of normalizing a quaternion.
/// - Geodesics, log and exp (`lie_dupin_transvection`, `lie_dupin_log`,
///   `lie_dupin_geodesic`) blend cyclides the way slerp blends rotations.
/// - A small correction toward new data lives in the 9-dimensional tangent
///   space of algebra elements that anticommute with sigma
///   (`lie_dupin_tangent_basis`, `lie_dupin_correct`).
///
/// The canonical torus (`lie_dupin_torus`) puts the tube spheres in U, so
/// `lie_dupin_curvatures(...).k_u` is the tube curvature.

namespace gaudi {
namespace albers {

typedef Eigen::Matrix<real, 21, 1> vec21;
typedef Eigen::Matrix<real, 21, 21> mat21;
typedef Eigen::Matrix<real, 6, 3> mat63;

struct lie_dupin {
  mat6 sigma = mat6::Identity();
};

/// eta-orthogonal projector onto the span of the columns of `B` (non-degenerate
/// under eta).
inline mat6 lie_projector(const mat63 &B) {
  const mat6 &eta = lie_metric();
  const mat3 G = B.transpose() * eta * B;
  return B * G.inverse() * B.transpose() * eta;
}

/// The cyclide whose U is the span of `B` (signature (2,1)).
inline lie_dupin lie_dupin_from_subspace(const mat63 &B) {
  return {real(2) * lie_projector(B) - mat6::Identity()};
}

/// Torus of revolution about the z axis through the origin, spine radius `R`
/// (> 0) and tube radius `r`, with outward normals (away from the spine).
/// U holds the tube spheres, centres on the spine circle, radius -r.
inline lie_dupin lie_dupin_torus(real R, real r) {
  mat63 B = mat63::Zero();
  B.col(0) << (real(1) + R * R - r * r) / real(2),
      (real(1) - R * R + r * r) / real(2), 0, 0, 0, -r;
  B(2, 1) = real(1);
  B(3, 2) = real(1);
  return lie_dupin_from_subspace(B);
}

/// Conjugation `A sigma A^-1` by a Lie sphere transformation.
inline lie_dupin lie_dupin_transform(const mat6 &A, const lie_dupin &d) {
  return {A * d.sigma * lie_inverse(A)};
}

/// Torus with centre `c`, axis `axis` (unit), radii `R`, `r`, outward normals.
inline lie_dupin lie_dupin_torus(const vec3 &c, const vec3 &axis, real R,
                                 real r) {
  const mat3 Rz =
      Eigen::Quaternion<real>::FromTwoVectors(vec3::UnitZ(), axis)
          .toRotationMatrix();
  return lie_dupin_transform(lie_rigid(Rz, c), lie_dupin_torus(R, r));
}

/// Largest of `|sigma^2 - I|`, the eta-asymmetry of sigma and `|tr sigma|`:
/// zero for a cyclide coordinate.
inline real lie_dupin_defect(const lie_dupin &d) {
  const mat6 &eta = lie_metric();
  const mat6 es = eta * d.sigma;
  const real inv = (d.sigma * d.sigma - mat6::Identity()).cwiseAbs().maxCoeff();
  const real sym = (es - es.transpose()).cwiseAbs().maxCoeff();
  return std::max({inv, sym, std::abs(d.sigma.trace())});
}

/// Signature (positive, negative) of eta on U (`sign = +1`) or U^perp
/// (`sign = -1`). (2, 1) on both for a cyclide.
inline Eigen::Vector2i lie_dupin_signature(const lie_dupin &d, real sign = 1,
                                           real tol = 1e-9) {
  // eta Pi is symmetric and has the inertia of eta restricted to range(Pi).
  const mat6 K = lie_metric() * (real(0.5) * (mat6::Identity() + sign * d.sigma));
  Eigen::SelfAdjointEigenSolver<mat6> es(real(0.5) * (K + K.transpose()));
  Eigen::Vector2i out(0, 0);
  const real scale = std::max(es.eigenvalues().cwiseAbs().maxCoeff(), real(1));
  for (int i = 0; i < 6; ++i) {
    if (es.eigenvalues()[i] > tol * scale)
      ++out[0];
    else if (es.eigenvalues()[i] < -tol * scale)
      ++out[1];
  }
  return out;
}

/// The two Euclidean-normalized basis vectors of a contact element and a basis
/// of its eta-orthogonal complement (4-dimensional, containing it).
struct lie_contact_frame {
  vec6 a[2];
  Eigen::Matrix<real, 6, 4> perp;
};

inline lie_contact_frame lie_contact_frame_of(const vec3 &p, const vec3 &n) {
  const lie_contact_element L = lie_contact(p, n);
  lie_contact_frame f;
  f.a[0] = L.point.normalized();
  f.a[1] = L.plane.normalized();
  Eigen::Matrix<real, 2, 6> A;
  A.row(0) = (lie_metric() * f.a[0]).transpose();
  A.row(1) = (lie_metric() * f.a[1]).transpose();
  Eigen::JacobiSVD<Eigen::Matrix<real, 2, 6>> svd(A, Eigen::ComputeFullV);
  f.perp = svd.matrixV().rightCols<4>();
  return f;
}

/// The eight membership residuals `<sigma a_i, w_k>` of a contact element:
/// all zero exactly when sigma maps L into L (`L^perp^perp = L`).
inline Eigen::Matrix<real, 8, 1>
lie_dupin_contact_rows(const mat6 &sigma, const lie_contact_frame &f) {
  Eigen::Matrix<real, 8, 1> r;
  const mat6 &eta = lie_metric();
  for (int i = 0; i < 2; ++i) {
    const vec6 sa = eta * (sigma * f.a[i]);
    for (int k = 0; k < 4; ++k)
      r[4 * i + k] = sa.dot(f.perp.col(k));
  }
  return r;
}

/// How far the contact element `(p, n)` is from lying on the cyclide: the
/// norm of its membership rows, 0 on the surface. Scale-free in the sense
/// that the rows use unit-norm vectors, but not a Euclidean distance.
inline real lie_dupin_contact_residual(const lie_dupin &d, const vec3 &p,
                                       const vec3 &n) {
  return lie_dupin_contact_rows(d.sigma, lie_contact_frame_of(p, n)).norm();
}

struct lie_dupin_curvature_pair {
  /// Normal curvature of the U family's sphere through `(p, n)`, and of the
  /// U^perp family's; infinite when that sphere is the point itself.
  real k_u = std::numeric_limits<real>::quiet_NaN();
  real k_v = std::numeric_limits<real>::quiet_NaN();
  /// Smallest singular value of each pencil solve, relative: 0 when the
  /// contact element meets U (resp. U^perp).
  real res_u = std::numeric_limits<real>::infinity();
  real res_v = std::numeric_limits<real>::infinity();
};

/// The principal curvatures of the cyclide at the contact element `(p, n)`:
/// the spheres `k P + T` of the contact element that lie in U and in U^perp.
/// Solved projectively, so a point sphere (k = inf) or the tangent plane
/// (k = 0) is no special case.
inline lie_dupin_curvature_pair lie_dupin_curvatures(const lie_dupin &d,
                                                     const vec3 &p,
                                                     const vec3 &n) {
  const lie_contact_element L = lie_contact(p, n);
  lie_dupin_curvature_pair out;
  for (int s = 0; s < 2; ++s) {
    const mat6 Q = mat6::Identity() - (s == 0 ? real(1) : real(-1)) * d.sigma;
    Eigen::Matrix<real, 6, 2> M;
    M.col(0) = Q * L.point;
    M.col(1) = Q * L.plane;
    Eigen::JacobiSVD<Eigen::Matrix<real, 6, 2>> svd(M, Eigen::ComputeFullV);
    const Eigen::Vector2d ab = svd.matrixV().col(1);
    const real k = std::abs(ab[1]) > real(1e-300)
                       ? ab[0] / ab[1]
                       : std::numeric_limits<real>::infinity();
    const real rel = svd.singularValues()[1] /
                     std::max(L.point.norm() + L.plane.norm(), real(1e-300));
    (s == 0 ? out.k_u : out.k_v) = k;
    (s == 0 ? out.res_u : out.res_v) = rel;
  }
  return out;
}

/// The nearest cyclide coordinate to an approximate one: symmetrize under eta,
/// then take the matrix sign by Newton's iteration `X <- (X + X^-1) / 2`,
/// which keeps eta-self-adjointness and converges quadratically. The analogue
/// of normalizing a quaternion. Does not force trace 0: a matrix whose sign
/// splits other than (3,3) is not a cyclide, and `lie_dupin_defect` reports it.
inline lie_dupin lie_dupin_normalize(const mat6 &S, int iterations = 30,
                                     real tol = 1e-14) {
  const mat6 &eta = lie_metric();
  mat6 X = real(0.5) * (S + eta * S.transpose() * eta);
  for (int it = 0; it < iterations; ++it) {
    const mat6 next = real(0.5) * (X + X.inverse());
    const real step = (next - X).cwiseAbs().maxCoeff();
    X = next;
    if (step < tol)
      break;
  }
  return {X};
}

/// Linear least-squares fit of a cyclide to weighted contact elements.
///
/// Writes `sigma = eta M` with M symmetric (21 unknowns), so every membership
/// row `<sigma a, w> = a^T M w` is linear in M. The rows' null space always
/// holds the identity (which fixes every contact element); a trace-0 row
/// removes it, the smallest eigenvector of the normal equations is the fit,
/// and `lie_dupin_normalize` puts it on the manifold.
struct lie_dupin_fit {
  mat21 N = mat21::Zero();
  int count = 0;

  static vec21 coefficient_row(const vec6 &a, const vec6 &w) {
    vec21 row;
    int k = 0;
    for (int i = 0; i < 6; ++i)
      for (int j = i; j < 6; ++j)
        row[k++] = i == j ? a[i] * w[i] : a[i] * w[j] + a[j] * w[i];
    return row;
  }

  static mat6 unpack(const vec21 &m) {
    mat6 M;
    int k = 0;
    for (int i = 0; i < 6; ++i)
      for (int j = i; j < 6; ++j) {
        M(i, j) = m[k];
        M(j, i) = m[k];
        ++k;
      }
    return M;
  }

  void accumulate(real weight, const vec3 &p, const vec3 &n) {
    const lie_contact_frame f = lie_contact_frame_of(p, n);
    for (int i = 0; i < 2; ++i)
      for (int k = 0; k < 4; ++k) {
        const vec21 row = coefficient_row(f.a[i], f.perp.col(k));
        N.noalias() += weight * row * row.transpose();
      }
    ++count;
  }

  struct result {
    lie_dupin cyclide;
    /// Ratio of the two smallest eigenvalues of the constrained normal
    /// equations: near 0 when the data picks one cyclide, near 1 when it
    /// leaves a pencil (fewer than four generic contact elements).
    real gap = 1;
    bool ok = false;
  };

  result solve(real gap_tol = 1e-3) const {
    result out;
    if (count < 4)
      return out;
    vec21 tr = vec21::Zero();
    {
      int k = 0;
      for (int i = 0; i < 6; ++i)
        for (int j = i; j < 6; ++j, ++k)
          if (i == j)
            tr[k] = lie_metric()(i, i);
    }
    const mat21 A = N + N.trace() * tr * tr.transpose();
    Eigen::SelfAdjointEigenSolver<mat21> es(A);
    if (es.info() != Eigen::Success)
      return out;
    const real l0 = std::max(es.eigenvalues()[0], real(0));
    const real l1 = es.eigenvalues()[1];
    out.gap = l1 > real(0) ? l0 / l1 : real(1);
    const mat6 S = lie_metric() * unpack(es.eigenvectors().col(0));
    const real s2 = (S * S).trace() / real(6);
    if (!(s2 > real(0)))
      return out;
    out.cyclide = lie_dupin_normalize(S / std::sqrt(s2));
    out.ok = out.gap < gap_tol && lie_dupin_defect(out.cyclide) < real(1e-8);
    return out;
  }
};

/// The transvection `g = (sigma_b sigma_a)^(1/2)`: the group element along
/// the geodesic from `a` to `b`, with `g sigma_a g^-1 = sigma_b`. Defined when
/// `sigma_b sigma_a` has no eigenvalues on the negative real axis (`a` and
/// `b` not antipodal). This is the edge transport `g_ij` between two vertex
/// cyclides.
inline mat6 lie_dupin_transvection(const lie_dupin &a, const lie_dupin &b) {
  const mat6 T = b.sigma * a.sigma;
  return T.sqrt();
}

/// Tangent vector at `a` pointing to `b`: `log g`, an algebra element that
/// anticommutes with `sigma_a`.
inline mat6 lie_dupin_log(const lie_dupin &a, const lie_dupin &b) {
  return lie_log(lie_dupin_transvection(a, b));
}

/// `exp(t xi) sigma_a exp(-t xi)` for a tangent vector `xi` at `a`.
inline lie_dupin lie_dupin_exp(const lie_dupin &a, const mat6 &xi,
                               real t = 1) {
  const mat6 g = lie_exp(t * xi);
  return {g * a.sigma * lie_inverse(g)};
}

/// The point at parameter `t` on the geodesic from `a` (t = 0) to `b`
/// (t = 1): the cyclide analogue of slerp.
inline lie_dupin lie_dupin_geodesic(const lie_dupin &a, const lie_dupin &b,
                                    real t) {
  return lie_dupin_exp(a, lie_dupin_log(a, b), t);
}

/// Holonomy around the face `(i, j, k)` of the edge transports:
/// `g_ki g_jk g_ij`. It fixes `sigma_i` (it lies in the stabilizer
/// O(2,1) x O(2,1) of the base), and is the identity when the three cyclides
/// lie on one geodesic or coincide.
inline mat6 lie_holonomy(const mat6 &g_ij, const mat6 &g_jk, const mat6 &g_ki) {
  return g_ki * g_jk * g_ij;
}

/// A basis of the 9-dimensional tangent space at `d`: the algebra elements
/// that anticommute with sigma, `(X - sigma X sigma) / 2` over the 15
/// generators, orthonormalized in the Frobenius inner product.
inline std::vector<mat6> lie_dupin_tangent_basis(const lie_dupin &d) {
  Eigen::Matrix<real, 36, 15> C;
  for (int i = 0, k = 0; i < 6; ++i)
    for (int j = i + 1; j < 6; ++j, ++k) {
      const mat6 X = lie_algebra_generator(i, j);
      const mat6 Y = real(0.5) * (X - d.sigma * X * d.sigma);
      C.col(k) = Eigen::Map<const Eigen::Matrix<real, 36, 1>>(Y.data());
    }
  Eigen::JacobiSVD<Eigen::Matrix<real, 36, 15>> svd(C, Eigen::ComputeThinU);
  std::vector<mat6> out;
  const real s0 = svd.singularValues()[0];
  for (int k = 0; k < 15; ++k) {
    if (svd.singularValues()[k] < real(1e-9) * s0)
      break;
    out.push_back(Eigen::Map<const mat6>(svd.matrixU().col(k).data()));
  }
  return out;
}

/// The smallest move of a cyclide toward weighted contact elements: damped
/// Gauss-Newton in the tangent space, `sigma <- c sigma c^-1` with `c` the
/// Cayley transform of `xi`, minimizing
///
///     sum_k w_k |rows_k(sigma(xi))|^2 + prior |xi|_F^2.
///
/// To first order `sigma(xi) = sigma + 2 xi sigma`. This is the corrector of
/// the advect-then-correct step: `prior` holds the advected state, and a
/// large one leaves it in place.
inline lie_dupin lie_dupin_correct(const lie_dupin &start,
                                   const std::vector<vec3> &p,
                                   const std::vector<vec3> &n,
                                   const std::vector<real> &w, real prior,
                                   int iterations = 4) {
  lie_dupin d = start;
  std::vector<lie_contact_frame> frames(p.size());
  for (size_t k = 0; k < p.size(); ++k)
    frames[k] = lie_contact_frame_of(p[k], n[k]);
  for (int it = 0; it < iterations; ++it) {
    const std::vector<mat6> basis = lie_dupin_tangent_basis(d);
    const int m = int(basis.size());
    matX H = matX::Identity(m, m) * prior;
    vecX g = vecX::Zero(m);
    std::vector<mat6> dsig(m);
    for (int b = 0; b < m; ++b)
      dsig[b] = real(2) * basis[b] * d.sigma;
    for (size_t k = 0; k < frames.size(); ++k) {
      const Eigen::Matrix<real, 8, 1> r0 = lie_dupin_contact_rows(d.sigma, frames[k]);
      Eigen::Matrix<real, 8, Eigen::Dynamic> J(8, m);
      for (int b = 0; b < m; ++b)
        J.col(b) = lie_dupin_contact_rows(dsig[b], frames[k]);
      H.noalias() += w[k] * J.transpose() * J;
      g.noalias() += w[k] * J.transpose() * r0;
    }
    const vecX x = H.ldlt().solve(-g);
    mat6 xi = mat6::Zero();
    for (int b = 0; b < m; ++b)
      xi += x[b] * basis[b];
    const mat6 c = lie_cayley(xi);
    d = lie_dupin_normalize(c * d.sigma * lie_inverse(c));
  }
  return d;
}

} // namespace albers
} // namespace gaudi

#endif
