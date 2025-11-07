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

index_t get_index(const std::vector<index_t> &indices, size_t i) {
  return indices[i];
}

template <IndexArray IndexType>
index_t get_index(const IndexType &indices, size_t i) {
  return get_index(indices, i);
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
template <int STRIDE, IndexArray IndexType>
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
  size_t stride() const { return STRIDE; }
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
  using value_type = typename ViewType::value_type;

  slice(const ViewType &data, size_t offset) : data_(&data), offset_(offset) {
    // console_logger::debug << "slice: offset: " << offset_ << std::endl;
  }

  const value_type &operator[](size_t i) const {
    return (*data_)[STRIDE * offset_ + i];
  }

  size_t size() const { return STRIDE; }
};

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
  const value_type &operator[](const size_t & i) const { return data_[i]; }
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
      : data_(&data),permutation_(&permutation), p_data_(data, permutation) {
        //assuming you'd use in same way as permuted_view
        //however no permutation, so get_index is just i
      }
  //copy constructor
  adjacency_view(const adjacency_view &other)
      : data_(other.data_),
        permutation_(other.permutation_),
        p_data_(other.p_data_) {}
  //assignment operator
  adjacency_view &operator=(const adjacency_view &other) {
    data_ = other.data_;
    permutation_ = other.permutation_;
    p_data_ = other.p_data_;
    return *this;
  }
  //move constructor
  adjacency_view(adjacency_view &&other)
      : data_(other.data_),
        permutation_(other.permutation_),
        p_data_(other.p_data_) {}
  //move assignment operator
  adjacency_view &operator=(adjacency_view &&other) {
    data_ = other.data_;
    permutation_ = other.permutation_;
    p_data_ = other.p_data_;
    return *this;
  }
  //destructor
  ~adjacency_view() {}
  //get inde
  template <int N>
  std::array<index_t, N> get_tuple_ids(const index_t & i){
    const std::vector<index_t> &edge_verts = data_;
    const std::vector<index_t> &permuted = indices_;
    std::array<index_t, N> tuple;
    for(int j = 0; j < N; j++){
      tuple[j] = edge_verts[2*permuted[i]+j];
    }
    return tuple;
  }
  
  const value_type &operator[](const size_t & i) const { return p_data_[i]; }
  //empty
  bool empty() const { return p_data_.empty(); }
  //size
  size_t size() const { return p_data_.size(); }
};



template <int STRIDE, ViewType ViewType, IndexArray IndexType>
class permuted_adjacency_view {
  const ViewType *data_;
  const IndexType *adjacency_;
  const IndexType *permutation_;
  const spread<STRIDE, IndexType> spread_;
  const permuted<IndexType, decltype(spread_)> p_adjacency_;
  const permuted<ViewType, decltype(p_adjacency_)> p_data_;

public:
  using value_type = typename ViewType::value_type;
  permuted_adjacency_view(const ViewType &data,      //
                          const IndexType &adjacency, //
                          const IndexType &permutation)
      : data_(&data),                     //
        adjacency_(&adjacency),           //
        permutation_(&permutation),          //
        spread_(permutation),            //
        p_adjacency_(adjacency, spread_), //
        p_data_(data, p_adjacency_) {}    //
  //copy constructor
  permuted_adjacency_view(const permuted_adjacency_view &other)
      : data_(other.data_),
        adjacency_(other.adjacency_),
        permutation_(other.permutation_),
        spread_(other.spread_),
        p_adjacency_(other.p_adjacency_),
        p_data_(other.p_data_) {}
  //assignment operator
  permuted_adjacency_view &operator=(const permuted_adjacency_view &other) {
    data_ = other.data_;
    adjacency_ = other.adjacency_;
    permutation_ = other.permutation_;
    spread_ = other.spread_;
    p_adjacency_ = other.p_adjacency_;
    p_data_ = other.p_data_;
    return *this;
  }
  //move constructor
  permuted_adjacency_view(permuted_adjacency_view &&other)
      : data_(other.data_),
        adjacency_(other.adjacency_),
        permutation_(other.permutation_),
        spread_(other.spread_),
        p_adjacency_(other.p_adjacency_),
        p_data_(other.p_data_) {}
  //move assignment operator
  permuted_adjacency_view &operator=(permuted_adjacency_view &&other) {
    data_ = other.data_;
    adjacency_ = other.adjacency_;
    permutation_ = other.permutation_;
    spread_ = other.spread_;
    p_adjacency_ = other.p_adjacency_;
    p_data_ = other.p_data_;
    return *this;
  }
  //destructor
  ~permuted_adjacency_view() {}
  //get index
  const value_type &operator[](const size_t & i) const { return p_data_[i]; }

  std::array<index_t, STRIDE> get_tuple_ids(const index_t & i){
    const std::vector<index_t> &edge_verts = adjacency_;
    const std::vector<index_t> &permuted = indices_;
    std::array<index_t, STRIDE> tuple;
    for(int j = 0; j < STRIDE; j++){
      tuple[j] = edge_verts[2*permuted[i]+j];
    }
    return tuple;
  }

  const index_t get_index(const size_t & i) const {
     return (*permutation_)[i]; }
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