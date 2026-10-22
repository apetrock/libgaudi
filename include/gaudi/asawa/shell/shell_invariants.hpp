#ifndef GAUDI_ASAWA_SHELL_SHELL_INVARIANTS_HPP
#define GAUDI_ASAWA_SHELL_SHELL_INVARIANTS_HPP

/// The shell's map invariants I1-I8 (`docs/dynamic_mesh_correctness_plan.md`
/// section 0), checked without aborting: walks are bounded by the corner
/// count and only follow links already verified, so a broken shell reports
/// instead of looping or reading out of range.

#include "gaudi/asawa/shell/shell.hpp"

#include <string>
#include <vector>

namespace gaudi {
namespace asawa {
namespace shell {

struct invariant_report {
  bool ok = true;
  /// "I1" .. "I8", or "tri" when a face is not a triangle.
  std::string which;
  /// The corner, face or vertex the violation names.
  int id = -1;
  std::string detail;
};

/// First violated invariant, or `ok`. `triangles` also requires every live
/// face to have exactly three corners (the remesher's state between ops).
inline invariant_report check_invariants(const shell &M, bool triangles = true) {
  const int n = static_cast<int>(M.corner_count());
  const int nf = static_cast<int>(M.face_count());
  const int nv = static_cast<int>(M.vert_count());
  auto fail = [](const char *which, int id, std::string detail) {
    return invariant_report{false, which, id, std::move(detail)};
  };
  auto in_range = [n](int c) { return c >= 0 && c < n; };
  auto live = [&](int c) { return in_range(c) && M.next(corner_id(c)) >= 0; };

  for (int i = 0; i < n; ++i) {
    const CornerId c = corner_id(i);
    if (!live(i)) {
      if (M.prev(c) != -1 || M.vert(c) != -1 || M.face(c) != -1)
        return fail("I3", i, "dead corner keeps a link or label");
      if (live(M.other(c)))
        return fail("I2", i, "twin of a dead corner is live");
      continue;
    }
    const int nx = M.next(c), pv = M.prev(c);
    if (!live(nx) || !live(pv))
      return fail("I1", i, "next or prev is not a live corner");
    if (M.prev(corner_id(nx)) != c || M.next(corner_id(pv)) != c)
      return fail("I1", i, "next and prev are not inverse");
    if (M.face(corner_id(nx)) != M.face(c))
      return fail("I4", i, "face label changes along next");
    if (M.vert(M.vnext(c)) != M.vert(c))
      return fail("I5", i, "vert label changes along vnext");
    if (M.vert(c) == M.vert(M.other(c)))
      return fail("I8", i, "loop edge");
    if (M.face(c) < 0 || M.face(c) >= nf || M.vert(c) < 0 || M.vert(c) >= nv)
      return fail("I6", i, "label out of range");
  }

  // One cycle per label (I4 / I5), each head on its own cycle (I6), face
  // lengths (I7). Every step below follows links verified above.
  std::vector<char> seen(static_cast<size_t>(n), 0);
  std::vector<int> face_cycles(static_cast<size_t>(nf), 0);
  for (int i = 0; i < n; ++i) {
    if (!live(i) || seen[static_cast<size_t>(i)])
      continue;
    int len = 0;
    CornerId c = corner_id(i);
    do {
      seen[static_cast<size_t>(static_cast<int>(c))] = 1;
      c = M.next(c);
      ++len;
    } while (c != corner_id(i) && len <= n);
    const int f = M.face(corner_id(i));
    if (++face_cycles[static_cast<size_t>(f)] > 1)
      return fail("I4", f, "two face cycles share one id");
    if (len < 3)
      return fail("I7", f, "face with " + std::to_string(len) + " corners");
    if (triangles && len != 3)
      return fail("tri", f, "face with " + std::to_string(len) + " corners");
  }
  std::fill(seen.begin(), seen.end(), 0);
  std::vector<int> vert_cycles(static_cast<size_t>(nv), 0);
  for (int i = 0; i < n; ++i) {
    if (!live(i) || seen[static_cast<size_t>(i)])
      continue;
    int len = 0;
    CornerId c = corner_id(i);
    do {
      seen[static_cast<size_t>(static_cast<int>(c))] = 1;
      c = M.vnext(c);
      ++len;
    } while (c != corner_id(i) && len <= n);
    const int v = M.vert(corner_id(i));
    if (++vert_cycles[static_cast<size_t>(v)] > 1)
      return fail("I5", v, "two vertex rings share one id");
  }
  for (int f = 0; f < nf; ++f) {
    const int b = M.fbegin(face_id(f));
    if ((b >= 0) != (face_cycles[static_cast<size_t>(f)] > 0))
      return fail("I6", f, "face head and face cycle disagree");
    if (b >= 0 && (!live(b) || M.face(corner_id(b)) != f))
      return fail("I6", f, "face head is not on its cycle");
  }
  for (int v = 0; v < nv; ++v) {
    const int b = M.vbegin(vert_id(v));
    if ((b >= 0) != (vert_cycles[static_cast<size_t>(v)] > 0))
      return fail("I6", v, "vertex head and vertex ring disagree");
    if (b >= 0 && (!live(b) || M.vert(corner_id(b)) != v))
      return fail("I6", v, "vertex head is not on its ring");
  }
  return {};
}

} // namespace shell
} // namespace asawa
} // namespace gaudi

#endif
