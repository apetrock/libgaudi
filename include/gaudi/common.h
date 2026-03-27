#ifndef __LIBGAUDI_COMMON_TYPEDEFS__
#define __LIBGAUDI_COMMON_TYPEDEFS__

#include "gaudi/console_logger.hpp"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <cassert>
#include <cmath>
#include <tuple>
#include <type_traits>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif
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

template <typename T>
concept ViewType = requires(T a, size_t i) {
  //  requires Indexable<T>;
  { a[i] } -> std::convertible_to<typename T::value_type>;
  { a.size() } -> std::convertible_to<size_t>;
  //{a.stride()} -> std::convertible_to<size_t>;
};

/*
template <int STRIDE, typename T>
concept SliceableViewType = requires(T a, size_t i) {
  //  requires Indexable<T>;
  { a[i] } -> std::convertible_to<typename T::value_type>;
  { a.size() } -> std::convertible_to<size_t>;
  {a.stride()} -> std::convertible_to<size_t>;
  {a.slice(i)} -> std::convertible_to<Slice<STRIDE, T>>;
};
*/        

// Concept for any array-like type (indexable + has size)
template <typename T>
concept TypeArray = requires(T t) {
  { t.size() } -> std::convertible_to<size_t>;
  typename T::value_type;
};

template <typename T>
concept Vec3View = TypeArray<T> && requires(T t) {
  { t[0] } -> std::convertible_to<vec3>;
  { t.size() } -> std::convertible_to<size_t>;
  //{ t.stride() } -> std::convertible_to<size_t>;
};

// Concept for index arrays (indexable with size_t values)
template <typename T>
concept IndexArray = requires(T t, size_t i) {
  { t[i] } -> std::convertible_to<size_t>;
  { t.size() } -> std::convertible_to<size_t>;
  //{ t.stride() } -> std::convertible_to<size_t>;
};

// Stride traits - extracts compile-time stride from view types
// Primary template: use T::stride member
template <typename T, typename = void>
struct view_stride {
  static constexpr size_t value = T::stride;
};

// Specialization for std::array
template <typename T, size_t N>
struct view_stride<std::array<T, N>> {
  static constexpr size_t value = N;
};

template <typename T>
inline constexpr size_t view_stride_v = view_stride<T>::value;

// Parameterized Vec3View concept with compile-time size checking
template <typename T, int N>
concept Vec3ViewN = Vec3View<T> && (view_stride_v<T> == N);

// Convenience aliases for geometric primitives
template <typename T>
concept PointView = Vec3ViewN<T, 1>;

template <typename T>
concept LineView = Vec3ViewN<T, 2>;

template <typename T>
concept TriView = Vec3ViewN<T, 3>;

// SimplexView concept - for types that return tuples of vec3s (e.g., simplex_view)
// Any type satisfying this concept works with bvh_tree<SimplexType>
// This includes views (simplex_view, permuted_simplex_view) and stored data
// (std::vector<std::array<vec3, N>>) for compile-time verification
template <typename T>
concept SimplexView = requires(T t, size_t i) {
  // Must have a stride member or be extractable via tuple_size
  { T::stride } -> std::convertible_to<size_t>;
  // Must return array-like tuple of vec3
  { t[i] } -> std::convertible_to<typename T::value_type>;
  { t.size() } -> std::convertible_to<size_t>;
} && requires {
  // value_type must be array-like (have tuple_size)
  typename T::value_type;
  { std::tuple_size<typename T::value_type>::value } -> std::convertible_to<size_t>;
};

// Helper to extract stride from SimplexView types
template <SimplexView T>
inline constexpr size_t simplex_stride_v = T::stride;

// Singulus: A single simplex that satisfies the SimplexView concept
// From Latin "singulus" - single, one at a time
// Used when passing a single simplex as a query to getNearest
template <size_t N>
struct Singulus {
  static constexpr size_t stride = N;
  using value_type = std::array<vec3, N>;
  
  const value_type* data_;
  
  Singulus() : data_(nullptr) {}
  explicit Singulus(const value_type& simplex) : data_(&simplex) {}
  
  value_type operator[](size_t i) const { 
    assert(i == 0 && "Singulus holds only a single simplex");
    return *data_; 
  }
  size_t size() const { return 1; }
};

inline index_t get_index(const std::vector<index_t> &indices, size_t i) {
  return indices[i];
}

template <IndexArray IndexType>
index_t get_index(const IndexType &indices, size_t i) {
  return static_cast<index_t>(indices[i]);
}

// Permuted array - templates on both array and index types
template <ViewType ViewType, IndexArray IndexType>
class permuted {
  const ViewType *data_;
  const IndexType *indices_;

public:
  using value_type = typename ViewType::value_type;

  permuted(const ViewType &data, const IndexType &indices)
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
template <ViewType ViewType, IndexArray IndexType>
permuted(const ViewType &, const IndexType &)
    -> permuted<ViewType, IndexType>;

// Spread array - templates on the array type
//spread is a view on an index array that adapts a permutation to work on 
//an adjacency list.
//spread<3>([1 0 2]) = [3 4 5 0 1 2 6 7 8]
//size(indices) = 1/STRIDE * size(indices)
//i_index = i_spread / STRIDE
//indices[i_index] += i_spread % STRIDE
//so basically we're taking a 


template <int STRIDE, IndexArray IndexType>
  requires IndexArray<IndexType>
class spread {
  const IndexType *indices_;

public:
  static constexpr size_t stride = STRIDE;
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
  size_t slice(size_t i) const { return i / STRIDE; }
};

// No deduction guide for spread - STRIDE must be explicit
// Usage: spread<STRIDE>(indices)

// Slice array - templates on the array type


template <int STRIDE, ViewType ViewType>
  requires TypeArray<ViewType>
class slice {
  const ViewType *data_;
  size_t offset_;

public:
  static constexpr size_t stride = STRIDE;
  using value_type = typename ViewType::value_type;

  slice(const ViewType &data, size_t offset) : data_(&data), offset_(offset) {
    // console_logger::debug << "slice: offset: " << offset_ << std::endl;
  }

  const value_type &operator[](size_t i) const {
    return (*data_)[STRIDE * offset_ + i];
  }

  size_t size() const { return STRIDE; }
};

template <int STRIDE, ViewType ViewType>
inline slice<STRIDE, ViewType> make_simplex(const ViewType &data,
                                            size_t offset) {
  return slice<STRIDE, ViewType>(data, offset);
}

// simplex_view - groups flat vec3 data into tuples of N vec3s
// Returns std::array<vec3, N> on access (stack-allocated, computed on-demand)
template <int N, Vec3View ViewType>
class simplex_view {
  const ViewType *data_;

public:
  static constexpr size_t stride = N;
  using value_type = std::array<vec3, N>;

  simplex_view(const ViewType &data) : data_(&data) {}

  // Copy/move constructors - views can be copied
  simplex_view(const simplex_view &other) = default;
  simplex_view(simplex_view &&other) = default;

  // Assignment operators deleted - views are immutable after construction
  simplex_view &operator=(const simplex_view &other) = delete;
  simplex_view &operator=(simplex_view &&other) = delete;

  ~simplex_view() = default;

  // Returns a stack-allocated tuple of N vec3s
  value_type operator[](size_t j) const {
    value_type result;
    for (int i = 0; i < N; ++i) {
      result[i] = (*data_)[N * j + i];
    }
    return result;
  }

  size_t size() const { return data_->size() / N; }
  bool empty() const { return data_->empty(); }
  
  // Access to underlying flat data
  const ViewType &data() const { return *data_; }
};

// Deduction guide for simplex_view - N must be explicit
// Usage: simplex_view<N>(data)

// Forward declare for use in permuted_simplex
template <SimplexView SimplexViewType, IndexArray IndexType>
class permuted_simplex;

// permuted_simplex - permutes a SimplexView (reorders the tuples)
// Returns tuples in permuted order
template <SimplexView SimplexViewType, IndexArray IndexType>
class permuted_simplex {
  const SimplexViewType *data_;
  const IndexType *permutation_;

public:
  static constexpr size_t stride = SimplexViewType::stride;
  using value_type = typename SimplexViewType::value_type;

  permuted_simplex(const SimplexViewType &data, const IndexType &permutation)
      : data_(&data), permutation_(&permutation) {}

  // Copy/move constructors - views can be copied
  permuted_simplex(const permuted_simplex &other) = default;
  permuted_simplex(permuted_simplex &&other) = default;

  // Assignment operators deleted - views are immutable after construction
  permuted_simplex &operator=(const permuted_simplex &other) = delete;
  permuted_simplex &operator=(permuted_simplex &&other) = delete;

  ~permuted_simplex() = default;

  // Returns the permutation[i]-th simplex
  value_type operator[](size_t i) const {
    return (*data_)[(*permutation_)[i]];
  }

  // Get the original index for a permuted index
  index_t get_index(size_t i) const {
    return static_cast<index_t>((*permutation_)[i]);
  }

  size_t size() const { return permutation_->size(); }
  bool empty() const { return permutation_->empty(); }

  // Access to underlying data
  const SimplexViewType &data() const { return *data_; }
  const IndexType &permutation() const { return *permutation_; }
};

// Deduction guide for permuted_simplex
template <SimplexView SimplexViewType, IndexArray IndexType>
permuted_simplex(const SimplexViewType &, const IndexType &)
    -> permuted_simplex<SimplexViewType, IndexType>;

// Free function to permute a SimplexView
template <SimplexView SimplexViewType, IndexArray IndexType>
inline auto permute(const SimplexViewType &data, const IndexType &permutation) {
  return permuted_simplex<SimplexViewType, IndexType>(data, permutation);
}

// permuted_simplex_view - convenience wrapper that composes:
//   1. permute(data, adjacency) - permute flat data by adjacency indices
//   2. simplex_view<N> - group into tuples
//   3. permute(simplex, permutation) - permute the tuples
// Returns tuples (std::array<vec3, N>) on access
// Drop-in replacement for permuted_adjacency_view (but returns tuples, not flat)
template <int N, Vec3View ViewType, IndexArray IndexType>
class permuted_simplex_view {
  const ViewType *data_;
  const IndexType *adjacency_;
  const IndexType *permutation_;

public:
  static constexpr size_t stride = N;
  using value_type = std::array<vec3, N>;

  permuted_simplex_view(const ViewType &data,
                        const IndexType &adjacency,
                        const IndexType &permutation)
      : data_(&data), adjacency_(&adjacency), permutation_(&permutation) {}

  // Copy/move constructors - views can be copied
  permuted_simplex_view(const permuted_simplex_view &other) = default;
  permuted_simplex_view(permuted_simplex_view &&other) = default;

  // Assignment operators deleted - views are immutable after construction
  permuted_simplex_view &operator=(const permuted_simplex_view &other) = delete;
  permuted_simplex_view &operator=(permuted_simplex_view &&other) = delete;

  ~permuted_simplex_view() = default;

  // Returns a stack-allocated tuple of N vec3s at permuted index i
  // Computes: for j in 0..N: result[j] = data[adjacency[permutation[i] * N + j]]
  value_type operator[](size_t i) const {
    value_type result;
    const size_t p_index = (*permutation_)[i];
    for (int j = 0; j < N; ++j) {
      result[j] = (*data_)[(*adjacency_)[N * p_index + j]];
    }
    return result;
  }

  // Get the tuple of original indices for a permuted index
  std::array<index_t, N> get_tuple_ids(index_t i) const {
    std::array<index_t, N> tuple;
    const size_t p_index = (*permutation_)[i];
    for (int j = 0; j < N; ++j) {
      tuple[j] = static_cast<index_t>((*adjacency_)[N * p_index + j]);
    }
    return tuple;
  }

  // Get the original simplex index for a permuted index
  index_t get_index(size_t i) const {
    return static_cast<index_t>((*permutation_)[i]);
  }

  size_t size() const { return permutation_->size(); }
  bool empty() const { return permutation_->empty(); }

  // Access to underlying data
  const ViewType &data() const { return *data_; }
  const IndexType &adjacency() const { return *adjacency_; }
  const IndexType &permutation() const { return *permutation_; }
};

// Deduction guide for permuted_simplex_view - N must be explicit
// Usage: permuted_simplex_view<N>(data, adjacency, permutation)

template <ViewType ViewType, IndexArray IndexType> class permuted_view {
  const ViewType *data_;
  const IndexType *permutation_;
  const permuted<ViewType, IndexType> p_data_;

public:
  using value_type = typename ViewType::value_type;
  permuted_view(const ViewType &data, const IndexType &permutation)
      : data_(&data),permutation_(&permutation), p_data_(data, permutation) {
    // console_logger::debug << "permuted_view: data: " << data.size() <<
    // std::endl; console_logger::debug << "permuted_view: indices: " <<
    // indices.size() << std::endl;
  }

  const index_t &get_index(size_t i) const { return (*permutation_)[i]; }
  const value_type &operator[](size_t i) const { return p_data_[i]; }
  bool empty() const { return p_data_.empty(); }
  size_t size() const { return p_data_.size(); }
};

//this will be an array wrapper for a view into 
//std::vector, maybe we'll do views into eigen matrices 

template <ViewType ViewType, IndexArray IndexType> class view {
  const ViewType *data_;
public:
  using value_type = typename ViewType::value_type;
  view(const ViewType &data)
      : data_(&data) {
        //assuming you'd use in same way as permuted_view
        //however no permutation, so get_index is just i
      }
  //copy constructor
  view(const view &other)
      : data_(other.data_) {}
  //assignment operator
  view &operator=(const view &other) {
    data_ = other.data_;
    return *this;
  }
  //move constructor
  view(view &&other)
      : data_(other.data_) {}
  //move assignment operator
  view &operator=(view &&other) {
    data_ = other.data_;
    return *this;
  }
  //destructor
  ~view() {}
  //get index
  index_t get_index(const size_t & i) const { return i; }
  //get value
  const value_type &operator[](const size_t & i) const { return (*data_)[i]; }
  //empty
  bool empty() const { return data_->empty(); }
  //size
  size_t size() const { return data_->size(); }
};

template <ViewType ViewType, IndexArray IndexType> class adjacency_view {
  const ViewType *data_;
  const IndexType *permutation_;
  const permuted<ViewType, IndexType> p_data_;

public:
  using value_type = typename ViewType::value_type;
  adjacency_view(const ViewType &data, const IndexType &permutation)
      : data_(&data),permutation_(&permutation), p_data_(data, permutation) {}
  
  // Copy/move constructors - views can be copied
  adjacency_view(const adjacency_view &other) = default;
  adjacency_view(adjacency_view &&other) = default;
  
  // Assignment operators deleted - views are immutable after construction
  adjacency_view &operator=(const adjacency_view &other) = delete;
  adjacency_view &operator=(adjacency_view &&other) = delete;
  
  ~adjacency_view() = default;
  //get index tuple from adjacency
  template <int N>
  std::array<index_t, N> get_tuple_ids(const index_t & i) const {

    std::array<index_t, N> tuple;
    for(int j = 0; j < N; j++){
      tuple[j] = (*permutation_)[N * i + j];
    }
    return tuple;
  }
  
  // get_index for unpermuted view is identity
  index_t get_index(const size_t & i) const { return i; }
  
  const value_type &operator[](const size_t & i) const { return p_data_[i]; }
  //empty
  bool empty() const { return p_data_.empty(); }
  //size
  size_t size() const { return p_data_.size(); }
};


//so what this does is it takes a data which is a set of points
//the adjacency is a flattened set of tuples of indices into the data
//the permutation then is designed to permute the adjacency to a new order
// DEPRECATED: Use permuted_simplex_view instead
// This class returns flat data and has circular dependency issues with slice
// Kept for backward compatibility - will be removed in future version
template <int STRIDE, ViewType ViewType, IndexArray IndexType>
class permuted_adjacency_view {
  const ViewType *data_;
  const IndexType *adjacency_;
  const IndexType *permutation_;
  const spread<STRIDE, IndexType> spread_;
  const permuted<IndexType, decltype(spread_)> p_adjacency_;
  const permuted<ViewType, decltype(p_adjacency_)> p_data_;

public:
  static constexpr size_t stride = STRIDE;
  using value_type = typename ViewType::value_type;
  permuted_adjacency_view(const ViewType &data,      //
                          const IndexType &adjacency, //
                          const IndexType &permutation)
      : data_(&data),                     //
        adjacency_(&adjacency),           //
        permutation_(&permutation),          //
        spread_(permutation),            //
        p_adjacency_(adjacency, spread_), //
        p_data_(data, p_adjacency_) {}
  
  // Copy/move constructors - views can be copied
  permuted_adjacency_view(const permuted_adjacency_view &other) = default;
  permuted_adjacency_view(permuted_adjacency_view &&other) = default;
  
  // Assignment operators deleted - views are immutable after construction
  permuted_adjacency_view &operator=(const permuted_adjacency_view &other) = delete;
  permuted_adjacency_view &operator=(permuted_adjacency_view &&other) = delete;
  
  ~permuted_adjacency_view() = default;
  //get index
  const value_type &operator[](const size_t & i) const { return p_data_[i]; }

  std::array<index_t, STRIDE> get_tuple_ids(const index_t & i){

    std::array<index_t, STRIDE> tuple;
    for(int j = 0; j < STRIDE; j++){
      const size_t p_index = (*permutation_)[i];
      tuple[j] = (*adjacency_)[STRIDE*p_index+j];
    }
    return tuple;
  }

  const index_t get_index(const size_t & i) const {
     return (*permutation_)[i]; }

  // NOTE: take() method removed to avoid circular dependency with slice<>
  // Use permuted_simplex_view instead for tuple access
  
  bool empty() const { return p_data_.empty(); }
  size_t size() const { return p_data_.size(); }
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
template <typename ViewType, typename IndexType>
struct tuple_size<gaudi::permuted<ViewType, IndexType>>
    : integral_constant<size_t, std::tuple_size_v<IndexType>> {};

template <int STRIDE, typename ViewType>
struct tuple_size<gaudi::slice<STRIDE, ViewType>>
    : integral_constant<size_t, STRIDE> {};

template <int STRIDE, typename IndexType>
struct tuple_size<gaudi::spread<STRIDE, IndexType>>
    : integral_constant<size_t, STRIDE * std::tuple_size_v<IndexType>> {};
} // namespace std

#endif