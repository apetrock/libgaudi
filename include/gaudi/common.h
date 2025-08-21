#ifndef __LIBGAUDI_COMMON_TYPEDEFS__
#define __LIBGAUDI_COMMON_TYPEDEFS__

#include "gaudi/console_logger.hpp"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <tuple>
#include <type_traits>
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

// Old slice class removed - replaced by the concept-based version below

template <typename T>
concept ArrayType = requires(T a, size_t i) {
  //  requires Indexable<T>;
  { a[i] } -> std::convertible_to<typename T::value_type>;
  { a.size() } -> std::convertible_to<size_t>;
};


/*
foo = [A B C]
adj = [2 0 1  1 2 0]
perm = [1 0]

[fcn, params: a, b, ... z]
bar = [permute, foo, [permute, adj [spread<3>, perm]]] => [B C A  C A B]

bar[0] = B
bar[1] = C
bar[2] = A
bar[3] = C
bar[4] = A
bar[5] = B
*/

// Concept for any array-like type (indexable + has size)
template <typename T>
concept TypeArray = requires(T t) {
  { t.size() } -> std::convertible_to<size_t>;
  typename T::value_type;
};

template <typename T>
concept Vec3Array = TypeArray<T> && requires(T t) {
  { t[0] } -> std::convertible_to<vec3>;
  { t.size() } -> std::convertible_to<size_t>;
};

// Concept for index arrays (indexable with size_t values)
template <typename T>
concept IndexArray = requires(T t, size_t i) {
  { t[i] } -> std::convertible_to<size_t>;
  { t.size() } -> std::convertible_to<size_t>;
};

// Permuted array - templates on both array and index types
template <typename ArrayType, typename IndexType>
  requires TypeArray<ArrayType> && IndexArray<IndexType>
class permuted {
  const ArrayType *data_;
  const IndexType *indices_;

public:
  using value_type = typename ArrayType::value_type;

  permuted(const ArrayType &data, const IndexType &indices)
      : data_(&data), indices_(&indices) {
    // console_logger::debug << "permuted: data: " << data.size() << std::endl;
    // console_logger::debug << "permuted: indices: " << indices.size() <<
    // std::endl;
  }

  const value_type &operator[](size_t i) const {
    return (*data_)[(*indices_)[i]];
  }
  bool empty() const { return indices_->empty(); }
  size_t size() const { return indices_->size(); }
};

// Deduction guide for permuted
template <typename ArrayType, typename IndexType>
permuted(const ArrayType &, const IndexType &)
    -> permuted<ArrayType, IndexType>;

// Spread array - templates on the array type
template <int STRIDE, typename IndexType>
  requires IndexArray<IndexType>
class spread {
  const IndexType *indices_;

public:
  using value_type = typename IndexType::value_type;

  spread(const IndexType &indices) : indices_(&indices) {
    // console_logger::debug << "spread: indices: " << indices.size() <<
    // std::endl;
  }

  index_t operator[](size_t i) const {
    return STRIDE * (*indices_)[i / STRIDE] + i % STRIDE;
  }
  bool empty() const { return indices_->empty(); }
  size_t size() const { return STRIDE * indices_->size(); }
};

// No deduction guide for spread - STRIDE must be explicit
// Usage: spread<STRIDE>(indices)

// Slice array - templates on the array type

template <typename T> class const_view {
  const T *data_; // Points to long-lived data

public:
  // Take reference to long-lived data, store as pointer
  using value_type = typename T::value_type;
  const_view(const T &data) : data_(&data) {}

  // Forward operations to pointed data
  auto operator[](size_t i) const -> decltype(auto) { return (*data_)[i]; }
  bool empty() const { return data_->empty(); }
  size_t size() const { return data_->size(); }
};

template <typename T> class ref_view {
  T data_; // Owns the data via move

public:
  using value_type = typename T::value_type;
  // Take temporary by rvalue reference and move it
  ref_view(T &&temp_data) : data_(std::move(temp_data)) {}

  // Copy constructor for when you need to copy
  ref_view(const T &data) : data_(data) {}

  // Forward operations to owned data
  auto operator[](size_t i) const -> decltype(auto) { return data_[i]; }
  bool empty() const { return data_.empty(); }
  size_t size() const { return data_.size(); }
};

template <typename T> auto make_view(T &&data) {
  if constexpr (std::is_lvalue_reference_v<T>) {
    // Lvalue reference -> const_view with pointer
    // console_logger::debug << "make_view: lvalue reference" << data.size() <<
    // std::endl;
    return const_view<std::remove_reference_t<T>>(std::forward<T>(data));
  } else {
    // console_logger::debug << "make_view: rvalue reference" << data.size() <<
    // std::endl;
    //  Rvalue reference -> ref_view with move
    return ref_view<std::remove_reference_t<T>>(std::forward<T>(data));
  }
}

template <ArrayType ArrayType, IndexArray IndexType> class permuted_view {
  const ArrayType *data_;
  const IndexType *permutation_;
  const permuted<ArrayType, IndexType> p_data_;

public:
  using value_type = typename ArrayType::value_type;
  permuted_view(const ArrayType &data, const IndexType &permutation)
      : data_(&data),permutation_(&permutation), p_data_(data, permutation) {
    // console_logger::debug << "permuted_view: data: " << data.size() <<
    // std::endl; console_logger::debug << "permuted_view: indices: " <<
    // indices.size() << std::endl;
  }
  const value_type &operator[](size_t i) const { return p_data_[i]; }
  bool empty() const { return p_data_.empty(); }
  size_t size() const { return p_data_.size(); }
};

template <ArrayType ArrayType, IndexArray IndexType>
using adjacency_view = permuted_view<ArrayType, IndexType>;

template <int STRIDE, ArrayType ArrayType, IndexArray IndexType>
class permuted_adjacency_view {
  const ArrayType *data_;
  const IndexType *adjacency_;
  const IndexType *permutation_;
  const spread<STRIDE, IndexType> spread_;
  const permuted<IndexType, decltype(spread_)> p_adjacency_;
  const permuted<ArrayType, decltype(p_adjacency_)> p_data_;

public:
  using value_type = typename ArrayType::value_type;
  permuted_adjacency_view(const ArrayType &data,      //
                          const IndexType &adjacency, //
                          const IndexType &permutation)
      : data_(&data),                     //
        adjacency_(&adjacency),           //
        permutation_(&permutation),          //
        spread_(permutation),            //
        p_adjacency_(adjacency, spread_), //
        p_data_(data, p_adjacency_) {}    //
  const value_type &operator[](size_t i) const { return p_data_[i]; }
  bool empty() const { return p_data_.empty(); }
  size_t size() const { return p_data_.size(); }
};

template <int STRIDE, typename ArrayType>
  requires TypeArray<ArrayType>
class slice {
  const ArrayType *data_;
  size_t offset_;

public:
  using value_type = typename ArrayType::value_type;

  slice(const ArrayType &data, size_t offset) : data_(&data), offset_(offset) {
    // console_logger::debug << "slice: offset: " << offset_ << std::endl;
  }

  const value_type &operator[](size_t i) const {
    return (*data_)[STRIDE * offset_ + i];
  }

  size_t size() const { return STRIDE; }
};

// No deduction guide for slice - STRIDE must be explicit
// Usage: slice<STRIDE>(data, offset)

// these will be function pointers with a default
// these should be non-const
template <int S, typename VEC> inline vecX to(const std::vector<VEC> &x) {
  std::vector<VEC> tmp(x);
  vecX out = Eigen::Map<vecX, Eigen::Unaligned>(
      reinterpret_cast<real *>(tmp.data()), S * x.size());
  return out;
}

template <int S, typename VEC>
inline void from(std::vector<VEC> &positions, const vecX &x) {
  for (int i = 0; i < positions.size(); i++)
    positions[i] = from<S, VEC>(x, i);
}

inline vecX to(const std::vector<vec3> &positions) {
  return to<3, vec3>(positions);
}
inline void from(std::vector<vec3> &positions, const vecX &x) {
  from<3, vec3>(positions, x);
}

inline vecX to(const std::vector<vec4> &positions) {
  return to<4, vec4>(positions);
}
inline void from(std::vector<vec4> &positions, const vecX &x) {
  from<4, vec4>(positions, x);
}

inline vecX to(const std::vector<quat> &positions) {
  return to<4, quat>(positions);
}
inline void from(std::vector<quat> &positions, const vecX &x) {
  from<4, quat>(positions, x);
}

inline vecX to(const std::vector<real> &U) {
  vecX Ue = Eigen::Map<const vecX, Eigen::Unaligned>(U.data(), U.size());
  return Ue;
}

inline vecX concat(const vecX &x, const vecX &u) {
  vecX q(x.size() + u.size());
  q << x, u;
  return q;
}

inline void split(const vecX &q, vecX &s, vecX &u) {
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

inline std::vector<real> from(vecX U) {
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

// Specialization for tuple_size to work with slice
namespace std {
template <typename ArrayType, typename IndexType>
struct tuple_size<gaudi::permuted<ArrayType, IndexType>>
    : integral_constant<size_t, std::tuple_size_v<IndexType>> {};

template <int STRIDE, typename ArrayType>
struct tuple_size<gaudi::slice<STRIDE, ArrayType>>
    : integral_constant<size_t, STRIDE> {};

template <int STRIDE, typename IndexType>
struct tuple_size<gaudi::spread<STRIDE, IndexType>>
    : integral_constant<size_t, STRIDE * std::tuple_size_v<IndexType>> {};
} // namespace std

#endif