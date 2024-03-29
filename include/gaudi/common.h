#ifndef __LIBGAUDI_COMMON_TYPEDEFS__
#define __LIBGAUDI_COMMON_TYPEDEFS__

#include "Eigen/src/Geometry/Quaternion.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <vector>

#define TYPEDEF_VEC(N) typedef Eigen::Matrix<real, N, 1> vec##N;
#define TYPEDEF_MAT(N) typedef Eigen::Matrix<real, N, N> mat##N;
#define TYPEDEF_MAT_NM(N, M) typedef Eigen::Matrix<real, N, M> mat##N##M;

namespace gaudi {
typedef int index_t;
typedef double real;
TYPEDEF_VEC(2)
TYPEDEF_VEC(3)
TYPEDEF_VEC(4)
TYPEDEF_VEC(6)
TYPEDEF_VEC(8)
TYPEDEF_VEC(10)
TYPEDEF_VEC(12)

typedef Eigen::Matrix<real, Eigen::Dynamic, 1> vecX;

TYPEDEF_MAT(2)
TYPEDEF_MAT(3)
TYPEDEF_MAT(4)
TYPEDEF_MAT(6)
TYPEDEF_MAT(8)
TYPEDEF_MAT(10)
TYPEDEF_MAT(12)

TYPEDEF_MAT_NM(3, 2)
TYPEDEF_MAT_NM(2, 3)
TYPEDEF_MAT_NM(3, 6)
TYPEDEF_MAT_NM(3, 9)
TYPEDEF_MAT_NM(4, 10)

typedef Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> matX;

typedef Eigen::Quaternion<real> quat;
typedef Eigen::SparseMatrix<real> matS;
typedef Eigen::Triplet<real> trip;

template <int S, typename VEC> const VEC from(const vecX &vals, size_t i) {
  return VEC(vals.data() + S * i);
};

// these will be function pointers with a default
// these should be non-const
template <int S, typename VEC> vecX to(const std::vector<VEC> &x) {
  std::vector<VEC> tmp(x);
  vecX out = Eigen::Map<vecX, Eigen::Unaligned>(
      reinterpret_cast<real *>(tmp.data()), S * x.size());
  return out;
}

template <int S, typename VEC>
void from(std::vector<VEC> &positions, const vecX &x) {
  for (int i = 0; i < positions.size(); i++)
    positions[i] = from<S, VEC>(x, i);
}

vecX to(const std::vector<vec3> &positions) { return to<3, vec3>(positions); }
void from(std::vector<vec3> &positions, const vecX &x) {
  from<3, vec3>(positions, x);
}

vecX to(const std::vector<vec4> &positions) { return to<4, vec4>(positions); }
void from(std::vector<vec4> &positions, const vecX &x) {
  from<4, vec4>(positions, x);
}

vecX to(const std::vector<quat> &positions) { return to<4, quat>(positions); }
void from(std::vector<quat> &positions, const vecX &x) {
  from<4, quat>(positions, x);
}

vecX to(const std::vector<real> &U) {
  vecX Ue = Eigen::Map<const vecX, Eigen::Unaligned>(U.data(), U.size());
  return Ue;
}

vecX concat(const vecX &x, const vecX &u) {
  vecX q(x.size() + u.size());
  q << x, u;
  return q;
}

void split(const vecX &q, vecX &s, vecX &u) {
  int Ns = s.size();
  int Nu = u.size();
  s = q.block(0, 0, Ns, 1);
  u = q.block(Ns, 0, Nu, 1);
}

/*
vecX concat(const std::vector<vec3> &xv, const std::vector<quat> &uv) {
  vecX s = to(xv);
  vecX u = to(uv);
  return concat(s, u);
}

void split(const vecX &q, std::vector<vec3> &sv, std::vector<quat> &uv) {
  vecX s;
  vecX u;
  split(q, s, u);
  from(sv, s);
  from(uv, u);
}
*/

std::vector<real> from(vecX U) {
  return std::vector<real>(U.data(), U.data() + U.rows() * U.cols());
}

template <typename T>
std::vector<T> operator+(const std::vector<T> &a, const std::vector<T> &b) {
  assert(a.size() == b.size());
  std::vector<T> result(a.size());
  std::transform(a.begin(), a.end(), b.begin(), result.begin(), std::plus<T>());

  return std::move(result);
}

template <typename T>
std::vector<T> operator-(const std::vector<T> &a, const std::vector<T> &b) {
  assert(a.size() == b.size());
  std::vector<T> result(a.size());
  std::transform(a.begin(), a.end(), b.begin(), result.begin(),
                 std::minus<T>());

  return std::move(result);
}

template <typename T>
std::vector<T> operator*(const real &a, const std::vector<T> &b) {
  std::vector<T> result(b.size());
  std::transform(b.begin(), b.end(), result.begin(),
                 [&a](const T &elem) { return a * elem; });

  return std::move(result);
}

} // namespace gaudi
#endif