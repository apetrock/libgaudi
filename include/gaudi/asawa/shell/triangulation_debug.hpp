#ifndef GAUDI_ASAWA_SHELL_TRIANGULATION_DEBUG_HPP
#define GAUDI_ASAWA_SHELL_TRIANGULATION_DEBUG_HPP

#include "gaudi/common.h"
#include "gaudi/asawa/shell/shell.hpp"

#include <iostream>
#include <map>
#include <ostream>
#include <vector>

namespace gaudi {
namespace asawa {
namespace shell {

enum class FaceRingResult {
  Ok,
  InvalidFbegin,
  NonPositiveNext,
  StepLimit,
};

/// Walk the face ring using next() only (no 100-step cap like for_each_face).
/// On success, \p corner_order and \p vert_order have one entry per side.
inline FaceRingResult walk_face_ring(const shell &M, index_t face_index,
                                     std::vector<index_t> &corner_order,
                                     std::vector<index_t> &vert_order,
                                     int max_steps = 1000000) {
  corner_order.clear();
  vert_order.clear();
  FaceId fi = face_id(face_index);
  CornerId c0 = M.fbegin(fi);
  if (c0 < 0)
    return FaceRingResult::InvalidFbegin;
  CornerId cur = c0;
  int steps = 0;
  do {
    corner_order.push_back(cur);
    vert_order.push_back(M.vert(cur));
    if (++steps > max_steps)
      return FaceRingResult::StepLimit;
    cur = M.next(cur);
    if (cur < 0)
      return FaceRingResult::NonPositiveNext;
  } while (cur != c0);
  return FaceRingResult::Ok;
}

/// Summarize topology before triangulation: boundary edges, polygon sizes,
/// invalid fbegin, broken rings, and corner/face id mismatches.
inline void dump_pre_triangulation_report(
    const shell &M, std::ostream &os,
    const std::vector<vec3> *positions = nullptr,
    std::size_t max_examples = 32) {
  const std::size_t nf = M.face_count();
  os << "[shell pre_triangulate] faces=" << nf
     << " corners=" << M.corner_count() << " verts=" << M.vert_count()
     << "\n";

  std::size_t boundary_pairs = 0;
  for (index_t e = 0; e < (index_t)M.corner_count(); e += 2) {
    CornerId ce = corner_id(e);
    CornerId o = M.other(ce);
    if (M.next(o) < 0)
      boundary_pairs++;
  }
  os << "[shell pre_triangulate] open boundary halfedge pairs: "
     << boundary_pairs
     << " (expect 0 for a closed watertight mesh after hole fill)\n";

  std::map<int, int> sides_hist;
  std::size_t invalid_fbegin = 0;
  std::size_t ring_fail = 0;
  std::size_t face_mismatch = 0;

  std::vector<index_t> corners, verts;

  for (index_t fi = 0; fi < (index_t)nf; ++fi) {
    FaceRingResult wr =
        walk_face_ring(M, fi, corners, verts, 1000000);
    if (wr == FaceRingResult::InvalidFbegin) {
      invalid_fbegin++;
      if (invalid_fbegin <= max_examples)
        os << "[shell pre_triangulate] face " << fi
           << ": invalid fbegin (no corner head / removed face), fbegin="
           << M.fbegin(face_id(fi)) << "\n";
      continue;
    }
    if (wr != FaceRingResult::Ok) {
      ring_fail++;
      if (ring_fail <= max_examples) {
        os << "[shell pre_triangulate] face " << fi
           << ": ring walk failed (result=" << static_cast<int>(wr)
           << ") fbegin=" << M.fbegin(face_id(fi)) << "\n";
      }
      continue;
    }

    int n = static_cast<int>(verts.size());
    int bucket = n >= 127 ? 127 : n;
    sides_hist[bucket]++;

    bool mismatch = false;
    for (index_t c : corners) {
      if (M.face(corner_id(c)) != fi) {
        mismatch = true;
        if (face_mismatch < max_examples)
          os << "[shell pre_triangulate] face " << fi << " corner " << c
             << " has corners_face=" << M.face(corner_id(c))             << " (expected " << fi << ")\n";
        break;
      }
    }
    if (mismatch)
      face_mismatch++;
  }

  os << "[shell pre_triangulate] summary: invalid_fbegin=" << invalid_fbegin
     << " ring_walk_failures=" << ring_fail
     << " corner_face_mismatches=" << face_mismatch << "\n";

  os << "[shell pre_triangulate] polygon side histogram (sides -> count):\n";
  for (const auto &kv : sides_hist) {
    if (kv.first == 127)
      os << "  " << kv.first << "+ : " << kv.second << "\n";
    else
      os << "  " << kv.first << " : " << kv.second << "\n";
  }

  if (positions && !positions->empty()) {
    std::size_t shown = 0;
    for (index_t fi = 0; fi < (index_t)nf && shown < max_examples; ++fi) {
      FaceRingResult wr =
          walk_face_ring(M, fi, corners, verts, 1000000);
      if (wr != FaceRingResult::Ok)
        continue;
      if (verts.size() == 3)
        continue;
      os << "[shell pre_triangulate] non-triangle face " << fi
         << " verts (indices):";
      for (index_t v : verts)
        os << " " << v;
      os << "\n";
      os << "    first positions:";
      for (std::size_t k = 0; k < verts.size() && k < 6; ++k) {
        index_t v = verts[k];
        if (v >= 0 && (std::size_t)v < positions->size())
          os << " [" << k << "]=" << (*positions)[v].transpose();
      }
      os << "\n";
      shown++;
    }
  }
  os << std::flush;
}

} // namespace shell
} // namespace asawa
} // namespace gaudi

#endif
