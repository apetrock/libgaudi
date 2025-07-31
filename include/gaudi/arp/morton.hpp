#ifndef __GAUDI_ARP_MORTON__
#define __GAUDI_ARP_MORTON__

#include "gaudi/common.h"
#include "gaudi/vec_addendum.h"
#include "gaudi/geometry_types.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace gaudi {
namespace arp {

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

// Calculates a 30-bit Morton code for the given 3D point located within the unit cube [0,1]
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
  std::string s;
  for (int j = 31; j >= 0; j--) {
    s += ((i & (1 << j)) ? '1' : '0');
  }
  return s;
}

// Normalize data to [0,1] range for Morton encoding
inline std::vector<vec3> normalize_data(const std::vector<vec3>& data) {
  if (data.empty()) return {};
  
  vec3 min_val = data[0];
  vec3 max_val = data[0];
  
  // Find bounds
  for (const auto& point : data) {
    min_val = va::min(min_val, point);
    max_val = va::max(max_val, point);
  }
  
  // Normalize to [0,1]
  std::vector<vec3> normalized_data;
  normalized_data.reserve(data.size());
  
  for (const auto& point : data) {
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
make_hash_3d(const std::vector<vec3>& data) {
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
    const vec3& point = normalized_data[i];
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
  
  return {sorted_hashes, indices};
}

template <int N, typename O>
inline std::vector<O> map(const std::vector<vec3>& data, std::function<O(const vec3&, const O &)> func, O default_val)  {
  if(data.empty()) {
    return {};
  }
  if(data.size() == 1) {
    return {func(data[0], default_val)};
  }
  if(data.size()% N != 0) {
    throw std::runtime_error("Data size must be a multiple of " + std::to_string(N));
  }
  std::vector<O> mapped;
  mapped.reserve(data.size() / N);
  for (size_t i = 0; i < data.size(); i += N) {
    O sum = default_val;
    for (size_t j = 0; j < N; ++j) {
      sum = func(data[i + j], sum);
    }
    mapped.push_back(sum);
  }
  return mapped;
}

using MassPoint = std::tuple<real, vec3>;

// Specialized mass calculations for N=1,2,3
template <int N>
inline std::vector<MassPoint> calc_com(const std::vector<vec3>& data) {
  if (data.size() % N != 0) {
    throw std::runtime_error("Data size must be a multiple of " + std::to_string(N));
  }
  
  std::vector<MassPoint> results;
  results.reserve(data.size() / N);
  
  for (size_t i = 0; i < data.size(); i += N) {
    real mass;
    vec3 com;
    
    if constexpr (N == 1) {
      mass = 1.0;
      com = data[i];
    } else if constexpr (N == 2) {
      vec3 edge = data[i + 1] - data[i];
      //mass = edge.norm(); test with uniform mass
      mass = 1.0;
      com = (data[i] + data[i + 1]) * 0.5;
    } else if constexpr (N == 3) {
      vec3 edge1 = data[i + 1] - data[i];
      vec3 edge2 = data[i + 2] - data[i];
      mass = 0.5 * edge1.cross(edge2).norm();
      com = (data[i] + data[i + 1] + data[i + 2]) / 3.0;
    }
    
    results.push_back({mass, com});
  }
  
  return results;
}

template <int N>
inline std::vector<ext::extents_t> calc_extents(const std::vector<vec3>& data)  {
  
  const auto map_fcn = [&](const vec3& a, const ext::extents_t& b) {
    return ext::expand(b, a);
  };
  
  const auto default_val = ext::init();
  return map<N, ext::extents_t>(data, map_fcn, default_val);
}

template <int N>
inline std::vector<mat3> calc_outers(const std::vector<vec3>& data, const std::vector<index_t>& indices)  {
  const auto map_fcn = [&](const vec3 & a, const mat3 &b){
    return b + a * a.transpose();
  };
  const auto default_val = mat3::Zero();
  return map<N, mat3>(data, map_fcn, default_val);
}
} // namespace arp
} // namespace gaudi

#endif // __GAUDI_ARP_MORTON__ 