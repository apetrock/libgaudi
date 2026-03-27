#ifndef __GAUDI_ARP_MORTON__
#define __GAUDI_ARP_MORTON__

#include "gaudi/common.h"
#include "gaudi/geometry_types.hpp"
#include "gaudi/vec_addendum.h"
#include "gaudi/console_logger.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <bitset>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#ifdef _MSC_VER
#include <intrin.h>
#endif

namespace gaudi {
namespace arp {

enum class MortonAxis { X, Y, Z };

template <size_t NumWords>
struct MortonCode {
  std::array<uint32_t, NumWords> words{};

  static constexpr size_t total_bits = NumWords * 32;
  static constexpr size_t bits_per_axis = total_bits / 3;
  static constexpr uint64_t scale_range_u64 = (bits_per_axis < 63) ? (1ULL << bits_per_axis) : 0ULL;

  MortonCode() = default;
  explicit MortonCode(const vec3 &p) : words(MortonCode::from(p).words) {}
  explicit operator vec3() const { return to_vec3(); }

  uint32_t &operator[](size_t i) { return words[i]; }
  const uint32_t &operator[](size_t i) const { return words[i]; }

  bool operator==(const MortonCode &rhs) const = default;
  bool operator<(const MortonCode &rhs) const { return words < rhs.words; }
  bool operator<=(const MortonCode &rhs) const { return !(rhs < *this); }
  bool operator>(const MortonCode &rhs) const { return rhs < *this; }
  bool operator>=(const MortonCode &rhs) const { return !(*this < rhs); }

  static MortonCode from_uint64(uint64_t v) {
    MortonCode out;
    if constexpr (NumWords >= 2) {
      out.words[NumWords - 1] = static_cast<uint32_t>(v & 0xFFFFFFFFULL);
      out.words[NumWords - 2] = static_cast<uint32_t>((v >> 32) & 0xFFFFFFFFULL);
    } else {
      out.words[0] = static_cast<uint32_t>(v & 0xFFFFFFFFULL);
    }
    return out;
  }

  uint64_t to_uint64() const {
    if constexpr (NumWords == 1) {
      return static_cast<uint64_t>(words[0]);
    } else {
      return (static_cast<uint64_t>(words[NumWords - 2]) << 32) |
             static_cast<uint64_t>(words[NumWords - 1]);
    }
  }

  static MortonCode from_quantized(uint64_t x, uint64_t y, uint64_t z) {
    MortonCode out;
    for (size_t bit = 0; bit < bits_per_axis; ++bit) {
      const uint64_t xb = (x >> bit) & 1ULL;
      const uint64_t yb = (y >> bit) & 1ULL;
      const uint64_t zb = (z >> bit) & 1ULL;
      const size_t z_pos = 3 * bit;
      const size_t y_pos = z_pos + 1;
      const size_t x_pos = z_pos + 2;
      if (z_pos < total_bits && zb)
        out.set_lsb_bit(z_pos);
      if (y_pos < total_bits && yb)
        out.set_lsb_bit(y_pos);
      if (x_pos < total_bits && xb)
        out.set_lsb_bit(x_pos);
    }
    return out;
  }

  static MortonCode from(const vec3 &p) {
#ifndef NDEBUG
    assert(std::isfinite(p[0]) && std::isfinite(p[1]) && std::isfinite(p[2]));
    assert(p[0] >= 0.0 && p[0] <= 1.0);
    assert(p[1] >= 0.0 && p[1] <= 1.0);
    assert(p[2] >= 0.0 && p[2] <= 1.0);
#endif

    const auto quantize = [](real in) -> uint64_t {
      real x = std::min(std::max(in, 0.0), 1.0);
      real scaled = x * static_cast<real>(scale_range_u64) + 0.5;
      scaled = std::min(std::max(scaled, 0.0), static_cast<real>(scale_range_u64 - 1ULL));
      return static_cast<uint64_t>(scaled);
    };

    return from_quantized(quantize(p[0]), quantize(p[1]), quantize(p[2]));
  }

  uint64_t deinterleave(MortonAxis axis) const {
    uint64_t out = 0;
    const size_t axis_shift = (axis == MortonAxis::X) ? 2 : ((axis == MortonAxis::Y) ? 1 : 0);
    for (size_t bit = 0; bit < bits_per_axis; ++bit) {
      const size_t pos = 3 * bit + axis_shift;
      if (pos < total_bits && get_lsb_bit(pos))
        out |= (1ULL << bit);
    }
    return out;
  }

  vec3 to_vec3() const {
    const auto decode = [](uint64_t q) -> real {
      return (static_cast<real>(q) + 0.5) / static_cast<real>(scale_range_u64);
    };
    return vec3(decode(deinterleave(MortonAxis::X)),
                decode(deinterleave(MortonAxis::Y)),
                decode(deinterleave(MortonAxis::Z)));
  }

  MortonCode bump(MortonAxis axis, int delta) const {
    const int shift = (axis == MortonAxis::X) ? 2 : ((axis == MortonAxis::Y) ? 1 : 0);
    const int64_t delta_scaled = static_cast<int64_t>(delta) << shift;
    MortonCode out = *this;

    if (delta_scaled >= 0) {
      uint64_t carry = static_cast<uint64_t>(delta_scaled);
      for (size_t i = NumWords; i-- > 0;) {
        const uint64_t wi = static_cast<uint64_t>(out.words[i]);
        const uint64_t add = carry & 0xFFFFFFFFULL;
        const uint64_t sum = wi + add;
        out.words[i] = static_cast<uint32_t>(sum & 0xFFFFFFFFULL);
        carry = (carry >> 32) + (sum >> 32);
      }
    } else {
      uint64_t sub = static_cast<uint64_t>(-delta_scaled);
      uint64_t borrow = 0;
      for (size_t i = NumWords; i-- > 0;) {
        const uint64_t wi = static_cast<uint64_t>(out.words[i]);
        const uint64_t to_sub = (sub & 0xFFFFFFFFULL) + borrow;
        const uint64_t diff = wi - to_sub;
        out.words[i] = static_cast<uint32_t>(diff & 0xFFFFFFFFULL);
        borrow = (wi < to_sub) ? 1ULL : 0ULL;
        sub >>= 32;
      }
    }
    return out;
  }

  MortonCode bump_next() const { return bump(MortonAxis::Z, 1); }

private:
  bool get_lsb_bit(size_t bit_index) const {
    const size_t lsw_word = NumWords - 1 - (bit_index / 32);
    const size_t in_word = bit_index % 32;
    return ((words[lsw_word] >> in_word) & 1u) != 0u;
  }

  void set_lsb_bit(size_t bit_index) {
    const size_t lsw_word = NumWords - 1 - (bit_index / 32);
    const size_t in_word = bit_index % 32;
    words[lsw_word] |= (1u << in_word);
  }
};

using Morton32 = MortonCode<1>;
using Morton64 = MortonCode<2>;
using Morton128 = MortonCode<4>;
using morton_t = Morton64;

template <size_t N>
inline MortonCode<N> xor_morton(const MortonCode<N> &a, const MortonCode<N> &b) {
  MortonCode<N> out;
  for (size_t i = 0; i < N; ++i) {
    out.words[i] = a.words[i] ^ b.words[i];
  }
  return out;
}

template <size_t N>
inline bool is_zero(const MortonCode<N> &a) {
  for (size_t i = 0; i < N; ++i) {
    if (a.words[i] != 0u)
      return false;
  }
  return true;
}

template <size_t N>
inline int clz_morton(const MortonCode<N> &a) {
  for (size_t i = 0; i < N; ++i) {
    const uint32_t w = a.words[i];
    if (w != 0u) {
#if defined(__GNUC__) || defined(__clang__)
      const int word_clz = __builtin_clz(w);
#elif defined(_MSC_VER)
      unsigned long idx;
      _BitScanReverse(&idx, w);
      const int word_clz = 31 - static_cast<int>(idx);
#else
      int word_clz = 0;
      uint32_t v = w;
      while ((v & 0x80000000u) == 0u) {
        v <<= 1;
        ++word_clz;
      }
#endif
      return static_cast<int>(i * 32 + static_cast<size_t>(word_clz));
    }
  }
  return static_cast<int>(N * 32);
}

template <size_t N>
inline std::string dump_binary(const MortonCode<N> &code) {
  std::ostringstream oss;
  for (size_t i = 0; i < N; ++i) {
    oss << std::bitset<32>(code.words[i]).to_string();
  }
  return oss.str();
}

// Expands a 10-bit integer into 30 bits by inserting 2 zeros after each bit
inline uint32_t expandBits(uint32_t v) {
  v = (v * 0x00010001u) & 0xFF0000FFu;
  v = (v * 0x00000101u) & 0x0F00F00Fu;
  v = (v * 0x00000011u) & 0xC30C30C3u;
  v = (v * 0x00000005u) & 0x49249249u;
  return v;
}

// Normalizes input to [0, c-1] range for Morton encoding
inline uint32_t scale(real x, real c) {
  return std::min(std::max(x * c, 0.0), c - 1.0);
}

// Calculates a 30-bit Morton code for the given 3D point located within the
// unit cube [0,1]
inline uint32_t morton3D(real x, real y, real z) {
  x = scale(x, 1024.0f);
  y = scale(y, 1024.0f);
  z = scale(z, 1024.0f);
  uint32_t xx = expandBits((uint32_t)x);
  uint32_t yy = expandBits((uint32_t)y);
  uint32_t zz = expandBits((uint32_t)z);
  return xx * 4 + yy * 2 + zz;
}

// Binary string representation for debugging Morton codes
inline std::string dump_binary(uint32_t i) {
  return std::bitset<32>(i).to_string();
}

inline std::string dump_binary(uint64_t i) {
  return std::bitset<64>(i).to_string();
}

// ---------------------------------------------------------------------------
// Named constexpr masks for 32-bit expand (10-bit -> 30-bit spread)
//
// Each stage doubles the gap between source bits by masking then multiplying.
// Input: 10 contiguous bits  (bits 0-9)
// Output: 10 bits spread across 30 positions (every 3rd bit set)
//
//   stage 0:  v = -------- ------98 76543210  (10 bits packed)
//   stage 1:  v = ------98 -------- 76543210  (groups of 8 separated by 8)
//   stage 2:  v = ----98-- --7654-- --3210--  (groups of 4 separated by 4)
//   stage 3:  v = --98--76 --54--32 --10----  (groups of 2 separated by 2)
//   stage 4:  v = -9-8-7-6 -5-4-3-2 -1-0----  (every 3rd bit)
// ---------------------------------------------------------------------------
namespace morton32 {
constexpr uint32_t mask0 = 0x000003FFu; // input: low 10 bits
constexpr uint32_t mul1  = 0x00010001u;
constexpr uint32_t mask1 = 0xFF0000FFu;
constexpr uint32_t mul2  = 0x00000101u;
constexpr uint32_t mask2 = 0x0F00F00Fu;
constexpr uint32_t mul3  = 0x00000011u;
constexpr uint32_t mask3 = 0xC30C30C3u;
constexpr uint32_t mul4  = 0x00000005u;
constexpr uint32_t mask4 = 0x49249249u; // final: every 3rd bit
constexpr real     scale_range = 1024.0; // 2^10 values per axis
} // namespace morton32

// ---------------------------------------------------------------------------
// 64-bit Morton encoding: 21 bits per axis -> 63-bit code
//
// Uses shift-OR approach (no multiply) since the gaps are wider.
// Input: 21 contiguous bits  (bits 0-20)
// Output: 21 bits spread across 63 positions (every 3rd bit set)
//
//   stage 0:  mask to 21 bits
//   stage 1:  shift 32, OR, mask -> groups separated by 32 zeros
//   stage 2:  shift 16, OR, mask -> groups separated by 16 zeros
//   stage 3:  shift 8,  OR, mask -> groups separated by 8 zeros
//   stage 4:  shift 4,  OR, mask -> groups separated by 4 zeros
//   stage 5:  shift 2,  OR, mask -> every 3rd bit
// ---------------------------------------------------------------------------
namespace morton64 {
constexpr uint64_t mask0 = 0x1fffffULL;               // input: low 21 bits
constexpr uint64_t mask1 = 0x1f00000000ffffULL;
constexpr uint64_t mask2 = 0x1f0000ff0000ffULL;
constexpr uint64_t mask3 = 0x100f00f00f00f00fULL;
constexpr uint64_t mask4 = 0x10c30c30c30c30c3ULL;
constexpr uint64_t mask5 = 0x1249249249249249ULL;      // final: every 3rd bit
constexpr real     scale_range = 2097152.0;             // 2^21 values per axis
} // namespace morton64

namespace morton128 {
constexpr uint64_t axis_bits_mask = 0x1fffffULL;           // 21 bits
constexpr std::array<uint32_t, 4> mask0_words = {
    0x00000000u, 0x000003FFu, 0xFFFFFFFFu, 0xFFFFFFFFu};   // low 42 bits set
constexpr real scale_range = 4398046511104.0;              // 2^42 values per axis
} // namespace morton128

struct u128_parts {
  uint64_t hi = 0;
  uint64_t lo = 0;
};

inline u128_parts u128_or(const u128_parts &a, const u128_parts &b) {
  return {a.hi | b.hi, a.lo | b.lo};
}

inline u128_parts u128_shift_left(const u128_parts &v, unsigned shift) {
  if (shift >= 128)
    return {};
  if (shift == 0)
    return v;
  if (shift >= 64) {
    return {v.lo << (shift - 64), 0};
  }
  return {(v.hi << shift) | (v.lo >> (64 - shift)), v.lo << shift};
}

inline Morton128 from_u128(const u128_parts &v) {
  Morton128 out;
  out.words[0] = static_cast<uint32_t>(v.hi >> 32);
  out.words[1] = static_cast<uint32_t>(v.hi & 0xFFFFFFFFULL);
  out.words[2] = static_cast<uint32_t>(v.lo >> 32);
  out.words[3] = static_cast<uint32_t>(v.lo & 0xFFFFFFFFULL);
  return out;
}

template <size_t Bits>
inline uint64_t quantize_axis(real in) {
  constexpr uint64_t range = (1ULL << Bits);
  real x = std::min(std::max(in, 0.0), 1.0);
  real scaled = x * static_cast<real>(range) + 0.5;
  scaled = std::min(std::max(scaled, 0.0), static_cast<real>(range - 1ULL));
  return static_cast<uint64_t>(scaled);
}

// Expands a 21-bit integer into 63 bits by inserting 2 zeros after each bit
inline uint64_t expandBits64(uint64_t v) {
  v &= morton64::mask0;
  v = (v | (v << 32)) & morton64::mask1;
  v = (v | (v << 16)) & morton64::mask2;
  v = (v | (v <<  8)) & morton64::mask3;
  v = (v | (v <<  4)) & morton64::mask4;
  v = (v | (v <<  2)) & morton64::mask5;
  return v;
}

// Morton128 fast spread path:
// split 42-bit input into two 21-bit halves, spread each half with 64-bit mask
// stages, then merge at +63-bit offset.
inline u128_parts expandBits128_masked(uint64_t v42) {
  uint64_t low21 = v42 & morton128::axis_bits_mask;
  uint64_t high21 = (v42 >> 21) & morton128::axis_bits_mask;
  u128_parts low = {0, expandBits64(low21)};
  u128_parts high = u128_shift_left({0, expandBits64(high21)}, 63);
  return u128_or(low, high);
}

inline Morton64 morton3D_64_ref(real x, real y, real z) {
  return Morton64::from(vec3(x, y, z));
}

inline Morton64 morton3D_64_fast(real x, real y, real z) {
  uint64_t qx = quantize_axis<21>(x);
  uint64_t qy = quantize_axis<21>(y);
  uint64_t qz = quantize_axis<21>(z);
  uint64_t xx = expandBits64(qx);
  uint64_t yy = expandBits64(qy);
  uint64_t zz = expandBits64(qz);
  return Morton64::from_uint64((xx << 2) | (yy << 1) | zz);
}

// Calculates a 63-bit Morton code for the given 3D point in the unit cube [0,1]
inline Morton64 morton3D_64(real x, real y, real z) {
  return morton3D_64_fast(x, y, z);
}

inline Morton128 morton3D_128_ref(real x, real y, real z) {
  return Morton128::from(vec3(x, y, z));
}

inline Morton128 morton3D_128_fast(real x, real y, real z) {
  uint64_t qx = quantize_axis<42>(x);
  uint64_t qy = quantize_axis<42>(y);
  uint64_t qz = quantize_axis<42>(z);
  u128_parts xx = u128_shift_left(expandBits128_masked(qx), 2);
  u128_parts yy = u128_shift_left(expandBits128_masked(qy), 1);
  u128_parts zz = expandBits128_masked(qz);
  return from_u128(u128_or(u128_or(xx, yy), zz));
}

// Calculates a 126-bit Morton code (42 bits/axis) for the unit cube [0,1]
inline Morton128 morton3D_128(real x, real y, real z) {
  return morton3D_128_fast(x, y, z);
}

// Normalize data to [0,1] range for Morton encoding
inline std::vector<vec3> normalize_data(const std::vector<vec3> &data) {
  if (data.empty())
    return {};

  vec3 min_val = data[0];
  vec3 max_val = data[0];

  // Find bounds
  for (const auto &point : data) {
    min_val = va::min(min_val, point);
    max_val = va::max(max_val, point);
  }

  // Normalize to [0,1]
  std::vector<vec3> normalized_data;
  normalized_data.reserve(data.size());

  for (const auto &point : data) {
    vec3 normalized = (point - min_val).cwiseQuotient(max_val - min_val);
    // Handle degenerate case where all points are the same
    for (int i = 0; i < 3; ++i) {
      if (std::isnan(normalized[i]) || std::isinf(normalized[i])) {
        normalized[i] = 0.0;
      }
    }
    normalized_data.push_back(normalized);
  }

  return normalized_data;
}

// Generate Morton codes for a set of 3D points
inline std::pair<std::vector<uint32_t>, std::vector<index_t>>
make_hash_3d(const std::vector<vec3> &data) {
  if (data.empty()) {
    return {{}, {}};
  }

  // Normalize data to [0,1] range
  std::vector<vec3> normalized_data = normalize_data(data);

  // Generate Morton codes
  std::vector<uint32_t> hashes;
  std::vector<index_t> indices;
  hashes.reserve(data.size());
  indices.reserve(data.size());

  for (index_t i = 0; i < static_cast<index_t>(normalized_data.size()); ++i) {
    const vec3 &point = normalized_data[i];
    uint32_t hash = morton3D(point[0], point[1], point[2]);
    hashes.push_back(hash);
    indices.push_back(i);
  }

  // Sort by hash values
  std::sort(indices.begin(), indices.end(),
            [&hashes](index_t a, index_t b) { return hashes[a] < hashes[b]; });

  // Reorder hashes to match sorted indices
  std::vector<uint32_t> sorted_hashes;
  sorted_hashes.reserve(hashes.size());
  for (index_t idx : indices) {
    sorted_hashes.push_back(hashes[idx]);
  }

  // Ensure unique codes: bump duplicates so the sorted sequence is strictly
  // monotonic. With 30-bit Morton codes and typical meshes the runs are short
  // (<=7 observed), so the spatial distortion is negligible.
  for (size_t i = 1; i < sorted_hashes.size(); i++) {
    if (sorted_hashes[i] <= sorted_hashes[i - 1])
      sorted_hashes[i] = sorted_hashes[i - 1] + 1;
  }

  return {sorted_hashes, indices};
}

template <typename MortonT>
inline std::pair<std::vector<MortonT>, std::vector<index_t>>
make_hash_3d_t(const std::vector<vec3> &data) {
  if (data.empty()) {
    return {{}, {}};
  }

  std::vector<vec3> normalized_data = normalize_data(data);
  std::vector<MortonT> hashes;
  std::vector<index_t> indices;
  hashes.reserve(data.size());
  indices.reserve(data.size());

  for (index_t i = 0; i < static_cast<index_t>(normalized_data.size()); ++i) {
    const vec3 &point = normalized_data[i];
    hashes.push_back(MortonT::from(point));
    indices.push_back(i);
  }

  std::sort(indices.begin(), indices.end(),
            [&hashes](index_t a, index_t b) {
              return (hashes[a] < hashes[b]) ||
                     (hashes[a] == hashes[b] && a < b);
            });

  std::vector<MortonT> sorted_hashes;
  sorted_hashes.reserve(hashes.size());
  for (index_t idx : indices) {
    sorted_hashes.push_back(hashes[idx]);
  }

  // Enforce strict monotonicity for radix-tree invariants.
  for (size_t i = 1; i < sorted_hashes.size(); i++) {
    if (sorted_hashes[i] <= sorted_hashes[i - 1]) {
      sorted_hashes[i] = sorted_hashes[i - 1].bump_next();
    }
  }

  return {sorted_hashes, indices};
}

// Generate 64-bit Morton codes for a set of 3D points.
inline std::pair<std::vector<Morton64>, std::vector<index_t>>
make_hash_3d_64(const std::vector<vec3> &data) {
  return make_hash_3d_t<Morton64>(data);
}

// Generate 128-bit Morton codes for a set of 3D points.
inline std::pair<std::vector<Morton128>, std::vector<index_t>>
make_hash_3d_128(const std::vector<vec3> &data) {
  return make_hash_3d_t<Morton128>(data);
}

// Map function for flat Vec3View (legacy)
template <int N, typename O, Vec3View TTYPE>
inline std::vector<O> map(const TTYPE &data,
                          auto &&func,
                          O default_val) {
  if (data.empty()) {
    return {};
  }
  if (data.size() == 1) {
    slice<N, TTYPE> datum(data, 0);
    return {func(datum, default_val)};
  }
  if (data.size() % N != 0) {
    throw std::runtime_error("Data size must be a multiple of " +
                             std::to_string(N));
  }
  std::vector<O> mapped;
  mapped.reserve(data.size() / N);
  for (size_t i = 0; i < data.size(); i += N) {
    O sum = default_val;
    slice<N, TTYPE> datum(data, i / N);
    sum = func(datum, sum);
    mapped.push_back(sum);
  }
  return mapped;
}

// Map function for SimplexView (type-based)
// Extracts stride from the SimplexView type automatically
template <typename O, SimplexView STYPE>
inline std::vector<O> map(const STYPE &data,
                          auto &&func,
                          O default_val) {
  if (data.empty()) {
    return {};
  }
  std::vector<O> mapped;
  mapped.reserve(data.size());
  for (size_t i = 0; i < data.size(); ++i) {
    auto simplex = data[i];
    O result = func(simplex, default_val);
    mapped.push_back(result);
  }
  return mapped;
}

using MassPoint = std::tuple<real, vec3>;

// Helper to compute center of mass for a simplex tuple
template <size_t N>
inline MassPoint compute_com(const std::array<vec3, N> &simplex) {
  real mass;
  vec3 com;
  if constexpr (N == 1) {
    mass = 1.0;
    com = simplex[0];
  } else if constexpr (N == 2) {
    vec3 edge = simplex[1] - simplex[0];
    mass = 1.0;
    com = (simplex[0] + simplex[1]) * 0.5;
  } else if constexpr (N == 3) {
    vec3 edge1 = simplex[1] - simplex[0];
    vec3 edge2 = simplex[2] - simplex[0];
    mass = 0.5 * edge1.cross(edge2).norm();
    com = (simplex[0] + simplex[1] + simplex[2]) / 3.0;
  } else {
    // Generic N: uniform mass, centroid
    mass = 1.0;
    com = vec3::Zero();
    for (size_t i = 0; i < N; ++i) {
      com += simplex[i];
    }
    com /= static_cast<real>(N);
  }
  return {mass, com};
}

// Specialized mass calculations for N=1,2,3 (legacy flat Vec3View)
template <int N, Vec3View TTYPE>
inline std::vector<MassPoint> calc_com(const TTYPE &data) {

  if (data.size() % N != 0) {
    throw std::runtime_error("data size must be a multiple of " +
                             std::to_string(N));
  }

  std::vector<MassPoint> results;
  results.reserve(data.size() / N);

  for (size_t i = 0; i < data.size(); i += N) {
    real mass;
    vec3 com;
    const int j = i / N;
    //console_logger::debug << "simplex i: " << i << " j: " << j << std::endl; 

    slice<N, TTYPE> simplex(data, j); 
    if constexpr (N == 1) {
      mass = 1.0;
      com = simplex[0];
    } else if constexpr (N == 2) {
      vec3 edge = simplex[1] - simplex[0];
      // mass = edge.norm(); test with uniform mass
      mass = 1.0;
      com = (simplex[0] + simplex[1]) * 0.5;
    } else if constexpr (N == 3) {
      vec3 edge1 = simplex[1] - simplex[0];
      vec3 edge2 = simplex[2] - simplex[0];
      mass = 0.5 * edge1.cross(edge2).norm();
      com = (simplex[0] + simplex[1] + simplex[2]) / 3.0;
    }

    results.push_back({mass, com});
  }

  return results;
}

// Center of mass calculation for SimplexView (type-based)
// Extracts stride from the SimplexView type automatically
template <SimplexView STYPE>
inline std::vector<MassPoint> calc_com(const STYPE &data) {
  constexpr size_t N = simplex_stride_v<STYPE>;
  
  std::vector<MassPoint> results;
  results.reserve(data.size());

  for (size_t i = 0; i < data.size(); ++i) {
    auto simplex = data[i];
    results.push_back(compute_com<N>(simplex));
  }

  return results;
}

// Extents calculation for flat Vec3View (legacy)
template <int N, Vec3View TTYPE>
inline std::vector<ext::extents_t> calc_extents(const TTYPE &data) {

  const auto map_fcn = [&](const slice<N, TTYPE> &a, const ext::extents_t &b) {
    auto out = b;
    for (int i = 0; i < N; ++i) {
      out = ext::expand(out, a[i]);
    }
    return out;
  };

  const auto default_val = ext::init();
  return map<N, ext::extents_t>(data, map_fcn, default_val);
}

// Helper to compute extents for a simplex tuple
template <size_t N>
inline ext::extents_t compute_extents(const std::array<vec3, N> &simplex) {
  auto out = ext::init();
  for (size_t i = 0; i < N; ++i) {
    out = ext::expand(out, simplex[i]);
  }
  return out;
}

// Extents calculation for SimplexView (type-based)
// Extracts stride from the SimplexView type automatically
template <SimplexView STYPE>
inline std::vector<ext::extents_t> calc_extents(const STYPE &data) {
  constexpr size_t N = simplex_stride_v<STYPE>;
  
  std::vector<ext::extents_t> results;
  results.reserve(data.size());

  for (size_t i = 0; i < data.size(); ++i) {
    auto simplex = data[i];
    results.push_back(compute_extents<N>(simplex));
  }

  return results;
}

template <int N, Vec3View TTYPE>
inline std::vector<mat3> calc_outers(const TTYPE &data) {
  const auto map_fcn = [&](const vec3 &a, const mat3 &b) {
    return b + a * a.transpose();
  };
  const auto default_val = mat3::Zero();
  return map<N, mat3>(data, map_fcn, default_val);
}

// Explicit instantiation declarations - controlled by CMake option
#if defined(GAUDI_USE_EXPLICIT_INSTANTIATIONS) &&                              \
    GAUDI_USE_EXPLICIT_INSTANTIATIONS
// Explicit instantiation declarations for calc_com
extern template std::vector<MassPoint> calc_com<1>(const std::vector<vec3> &);
extern template std::vector<MassPoint> calc_com<2>(const std::vector<vec3> &);
extern template std::vector<MassPoint> calc_com<3>(const std::vector<vec3> &);

// Explicit instantiation declarations for calc_extents
extern template std::vector<ext::extents_t>
calc_extents<1>(const std::vector<vec3> &);
extern template std::vector<ext::extents_t>
calc_extents<2>(const std::vector<vec3> &);
extern template std::vector<ext::extents_t>
calc_extents<3>(const std::vector<vec3> &);
#endif

} // namespace arp
} // namespace gaudi

#endif // __GAUDI_ARP_MORTON__