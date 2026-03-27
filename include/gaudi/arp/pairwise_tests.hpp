#ifndef __GAUDI_ARP_PAIRWISE_TESTS__
#define __GAUDI_ARP_PAIRWISE_TESTS__

#include "gaudi/common.h"
#include "gaudi/vec_addendum.h"
#include <array>
#include <limits>

namespace gaudi {
namespace arp {

// Thin inline wrappers around vec_addendum.h primitives for pairwise distance
// tests. These use typed concepts (PointView, LineView, TriView) to enforce
// correct primitive types at compile time.

// Point-to-point distance
template <PointView P0, PointView P1>
inline real test_point_point(const P0 &p0, const P1 &p1) {
  const vec3 &p00 = p0[0];
  const vec3 &p10 = p1[0];
  return (p00 - p10).norm();
}

// Point-to-line distance (wraps va::distance_from_line)
template <PointView P0, LineView P1>
inline real test_point_line(const P0 &pA, const P1 &pB) {
  const vec3 &pA0 = pA[0];
  const vec3 &pB0 = pB[0];
  const vec3 &pB1 = pB[1];
  return va::distance_from_line(pB0, pB1, pA0);
}

// Point-to-triangle distance (wraps va::closest_point)
template <PointView P0, TriView P1>
inline real test_point_tri(const P0 &pA, const P1 &pB) {
  const vec3 &point = pA[0];
  std::array<vec3, 3> tri = {pB[0], pB[1], pB[2]};
  auto result = va::closest_point(tri, point);
  return result[0]; // distance is first element
}

// Line-to-line distance (wraps va::distance_Segment_Segment)
template <LineView P0, LineView P1>
inline real test_line_line(const P0 &pA, const P1 &pB) {
  const vec3 &pA0 = pA[0];
  const vec3 &pA1 = pA[1];
  const vec3 &pB0 = pB[0];
  const vec3 &pB1 = pB[1];

  std::array<real, 3> d = va::distance_Segment_Segment(pA0, pA1, pB0, pB1);
  return d[0];
}

// Line-to-line distance with colinearity filter
// Filters out adjacent edges by checking if the connecting vector
// is too aligned with either edge direction
template <LineView P0, LineView P1>
inline real test_line_line_filtered(const P0 &pA, const P1 &pB,
                                    real angle_threshold = 0.35) {
  const vec3 &pA0 = pA[0];
  const vec3 &pA1 = pA[1];
  const vec3 &pB0 = pB[0];
  const vec3 &pB1 = pB[1];

  std::array<real, 3> d = va::distance_Segment_Segment(pA0, pA1, pB0, pB1);
  real s = d[1];
  real t = d[2];

  vec3 xA = va::mix(s, pA0, pA1);
  vec3 xB = va::mix(t, pB0, pB1);
  vec3 dA = (pA1 - pA0).normalized();
  vec3 dB = (pB1 - pB0).normalized();

  vec3 xAB = (xB - xA).normalized();

  // Filter out colinear cases (typically adjacent edges)
  if (std::abs(dA.dot(xAB)) > angle_threshold)
    return std::numeric_limits<real>::max();
  if (std::abs(dB.dot(xAB)) > angle_threshold)
    return std::numeric_limits<real>::max();

  return d[0];
}

// Line-to-triangle distance (placeholder - not yet implemented)
template <LineView P0, TriView P1>
inline real test_line_tri(const P0 &p0, const P1 &p1) {
  // TODO: Implement using ray-triangle intersection or edge-edge tests
  return std::numeric_limits<real>::max();
}

// Triangle-to-triangle distance (placeholder - not yet implemented)
template <TriView P0, TriView P1>
inline real test_tri_tri(const P0 &p0, const P1 &p1) {
  // TODO: Implement using separating axis theorem or GJK
  return std::numeric_limits<real>::max();
}

// ============================================================================
// Tuple-based test functions for SimplexView types
// These operate on std::array<vec3, N> directly
// ============================================================================

// Point-to-point distance (tuple version)
inline real test_point_point_tuple(const std::array<vec3, 1> &p0,
                                   const std::array<vec3, 1> &p1) {
  return (p0[0] - p1[0]).norm();
}

// Point-to-line distance (tuple version)
inline real test_point_line_tuple(const std::array<vec3, 1> &pA,
                                  const std::array<vec3, 2> &pB) {
  return va::distance_from_line(pB[0], pB[1], pA[0]);
}

// Point-to-triangle distance (tuple version)
inline real test_point_tri_tuple(const std::array<vec3, 1> &pA,
                                 const std::array<vec3, 3> &pB) {
  auto result = va::closest_point(pB, pA[0]);
  return result[0]; // distance is first element
}

// Line-to-line distance (tuple version)
inline real test_line_line_tuple(const std::array<vec3, 2> &pA,
                                 const std::array<vec3, 2> &pB) {
  std::array<real, 3> d = va::distance_Segment_Segment(pA[0], pA[1], pB[0], pB[1]);
  return d[0];
}

// Line-to-triangle distance (tuple version, placeholder)
inline real test_line_tri_tuple(const std::array<vec3, 2> &p0,
                                const std::array<vec3, 3> &p1) {
  // TODO: Implement using ray-triangle intersection or edge-edge tests
  return std::numeric_limits<real>::max();
}

// Triangle-to-triangle distance (tuple version, placeholder)
inline real test_tri_tri_tuple(const std::array<vec3, 3> &p0,
                               const std::array<vec3, 3> &p1) {
  // TODO: Implement using separating axis theorem or GJK
  return std::numeric_limits<real>::max();
}

} // namespace arp
} // namespace gaudi

#endif // __GAUDI_ARP_PAIRWISE_TESTS__

