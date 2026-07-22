
#ifndef __ASAWA_PRIM_OPS__
#define __ASAWA_PRIM_OPS__

#include "shell.hpp"

#include "../datums.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <numeric>
#include <ostream>
#include <vector>

namespace gaudi {
namespace asawa {
namespace shell {

using index_t = int;

CornerId split_edge(shell &M, CornerId corner_index,
                    VertId new_vertex = vert_id(-1),
                    CornerId new_corner = corner_id(-1)) {
  CornerId c0 = corner_index;
  CornerId c1 = M.other(c0);

  CornerId c0n = M.next(c0);
  CornerId c1p = M.prev(c1);

  VertId vn = new_vertex < 0 ? M.insert_vertex() : new_vertex;
  VertId v1 = M.vert(c1);

  CornerId c0i = new_corner < 0 ? M.insert_edge_pair() : new_corner;
  CornerId c1i = M.other(c0i);

  M.link(c0, c0i);
  M.link(c0i, c0n);
  M.link(c1p, c1i);
  M.link(c1i, c1);
  M.set_face(c0i, M.face(c0));
  M.set_face(c1i, M.face(c1));

  M.set_vert(c1, vn);
  M.set_vert(c0i, vn);
  M.set_vert(c1i, v1);

  return c0i;
}

int count_cycle(shell &M, CornerId corner) {
  CornerId c0 = corner;
  CornerId c1 = M.other(c0);
  VertId v0 = M.vert(c0);
  VertId v1 = M.vert(c1);
  int t = 0;
  M.for_each_vertex(v0, [&t, v1](CornerId ci, shell &M) {
    VertId vi = M.vert(M.next(ci));
    t += int(vi == v1);
  });
  return t;
}

CornerId merge_edge(shell &M, CornerId corner_index) {
  CornerId c0 = corner_index;
  CornerId c1 = M.other(c0);
  CornerId c0n = M.next(c0);
  CornerId c0p = M.prev(c0);

  CornerId c1n = M.next(c1);
  CornerId c1p = M.prev(c1);

  VertId v0 = M.vert(c0);
  VertId v1 = M.vert(c1);
  FaceId f0 = M.face(c0);
  FaceId f1 = M.face(c1);

  if (v0 != v1)
    M.remove_vertex(v0);
  M.remove_edge_pair(c0);

  M.link(c0p, c0n);
  M.link(c1p, c1n);

  M.set_vert(c1n, v1);

  M.set_face(c0n, f0);
  M.set_face(c1n, f1);
  M.vupdate(v1);

  return c1n;
}

CornerId split_face(shell &M, CornerId c0, CornerId c1,
                    CornerId new_corner = corner_id(-1),
                    FaceId new_face = face_id(-1)) {

  FaceId nf = new_face < 0 ? M.insert_face() : new_face;
  FaceId f0 = M.face(c0);
  VertId v0 = M.vert(c0);
  VertId v1 = M.vert(c1);

  CornerId c0p = M.prev(c0);
  CornerId c1p = M.prev(c1);

  CornerId c0i = new_corner < 0 ? M.insert_edge_pair() : new_corner;
  CornerId c1i = M.other(c0i);

  M.link(c0p, c0i);
  M.link(c0i, c1);

  M.link(c1p, c1i);
  M.link(c1i, c0);

  M.set_face(c0, nf);
  M.set_face(c1i, nf);
  M.set_face(c0i, f0);

  M.set_vert(c0i, v0);
  M.set_vert(c1i, v1);

  M.fupdate(nf);
  M.fupdate(f0);
  M.vupdate(v0);

  return c0i;
}

FaceId remove_dangling(shell &M, CornerId c0, CornerId c1) {
  CornerId c0p = M.prev(c0);
  CornerId c0n = M.next(c0);
  CornerId c1p = M.prev(c1);
  FaceId f0 = M.face(c0);
  VertId v0 = M.vert(c0);
  VertId v1 = M.vert(c1);

  M.remove_vertex(v0);
  M.remove_edge_pair(c0);
  M.set_vert(c0n, v1);
  M.link(c1p, c0n);
  M.set_face(c0n, f0);

  M.vupdate(v1);

  return f0;
}

FaceId remove_cruft(shell &M, CornerId c0, CornerId c1) {
  FaceId f0 = M.face(c0);
  FaceId f1 = M.face(c1);
  VertId v0 = M.vert(c0);
  VertId v1 = M.vert(c1);

  M.remove_vertex(v0);
  M.remove_vertex(v1);
  M.remove_face(f0);
  M.remove_face(f1);

  M.remove_edge_pair(c0);
  return face_id(-1);
}

FaceId merge_face(shell &M, CornerId c0, CornerId c1) {
  assert(c1 == M.other(c0));
  CornerId c0p = M.prev(c0);
  CornerId c0n = M.next(c0);
  CornerId c1p = M.prev(c1);
  CornerId c1n = M.next(c1);

  FaceId f0 = M.face(c0);
  FaceId f1 = M.face(c1);

  if (f0 == f1) {
    if (c0 == c1n && c0n == c1) {
      return remove_cruft(M, c0, c1);
    }
    if (c0 == c1n) {
      return remove_dangling(M, c0, c1);
    }
    if (c1 == c0n) {
      return remove_dangling(M, c1, c0);
    }
  }

  VertId v0 = M.vert(c0);
  VertId v1 = M.vert(c1);
  M.remove_face(f1);
  M.remove_edge_pair(c0);

  M.link(c0p, c1n);
  M.link(c1p, c0n);
  M.set_vert(c1n, v0);
  M.set_vert(c0n, v1);
  M.set_face(c1n, f0);

  M.fupdate(f0);

  return f0;
}

CornerId subdivide_edge(shell &M, CornerId corner_index,
                        VertId new_vert = vert_id(-1),
                        CornerId new_corner0 = corner_id(-1),
                        CornerId new_corner1 = corner_id(-1),
                        CornerId new_corner2 = corner_id(-1),
                        FaceId new_face0 = face_id(-1),
                        FaceId new_face1 = face_id(-1)) {

  CornerId c0 = corner_index;
  CornerId c1 = M.next(M.other(c0));

  c0 = split_edge(M, c0, new_vert, new_corner0);
  c1 = M.next(M.other(c0));

  split_face(M, c0, M.next(M.next(c0)), new_corner1, new_face0);
  split_face(M, c1, M.next(M.next(c1)), new_corner2, new_face1);

  return c0;
}

CornerId collapse_edge(shell &M, CornerId corner_index,
                       bool degenerate = false) {

  CornerId c0 = corner_index;
  CornerId c1 = M.other(c0);

  CornerId c0p = M.prev(c0);
  CornerId c1p = M.prev(c1);

  if (count_cycle(M, c0) > 1 && !degenerate)
    return corner_index;

  if (M.vsize(M.vert(c0)) < 3)
    return c0;
  if (M.vsize(M.vert(c1)) < 3)
    return c0;
  if (M.vsize(M.vert(c0p)) < 4)
    return c0;
  if (M.vsize(M.vert(c1p)) < 4)
    return c0;

  merge_face(M, c0p, M.other(c0p));
  merge_face(M, c1p, M.other(c1p));

  c0 = merge_edge(M, c0);

  assert(M.vert(c0) != M.vert(M.next(c0)));

  return c0;
}

CornerId flip_edge(shell &M, CornerId corner_index) {
  CornerId c0i = corner_index;
  CornerId c0p = M.prev(c0i);
  CornerId c0pp = M.prev(c0p);
  CornerId c0n = M.next(c0i);

  CornerId c1i = M.other(c0i);
  CornerId c1p = M.prev(c1i);
  CornerId c1pp = M.prev(c1p);
  CornerId c1n = M.next(c1i);

  FaceId f0 = M.face(c0i);
  FaceId f1 = M.face(c1i);

  VertId v0 = M.vert(c0i);
  VertId v1 = M.vert(c0p);
  VertId v2 = M.vert(c1i);
  VertId v3 = M.vert(c1p);

  M.link(c0pp, c0i);
  M.link(c0i, c1p);
  M.link(c1p, c0n);

  M.link(c1pp, c1i);
  M.link(c1i, c0p);
  M.link(c0p, c1n);

  M.set_face(c0i, f0);
  M.set_face(c1i, f1);

  M.set_vert(c0i, v1);
  M.set_vert(c1i, v3);

  M.set_vert(c1n, v0);
  M.set_vert(c0n, v2);

  M.vupdate(v1);
  M.vupdate(v3);
  M.fupdate(f0);
  M.fupdate(f1);

  if (M.fsize(f0) != 3)
    M.fprintv(f0);
  assert(M.fsize(f0) == 3);
  if (M.fsize(f1) != 3)
    M.fprintv(f1);
  assert(M.fsize(f1) == 3);

  return c0i;
}

FaceId remove_vertex(shell &M, VertId v) {

  std::vector<CornerId> corners;
  M.for_each_vertex(v, [&corners](CornerId cid, shell &m) {
    corners.push_back(cid);
  });
  FaceId f = face_id(-1);
  for (CornerId c : corners)
    f = merge_face(M, c, M.other(c));

  return f;
}

void triangulate_face(shell &M, FaceId fi) {
  CornerId c0 = M.fbegin(fi);
  if (c0 < 0)
    return;

  bool splitting = true;
  while (splitting) {

    CornerId c0n = M.next(c0);
    CornerId c0nn = M.next(c0n);
    CornerId c0p = M.prev(c0);

    splitting = c0nn != c0p;
    if (splitting) {
      c0 = split_face(M, c0n, c0p);
      c0 = M.other(c0);
    }
  }
}

void triangulate_face(shell &M, FaceId fi);

void remove_orphan_vertices(shell &M);

void triangulate(shell &M) {
  size_t Nf = M.face_count();
  for (int i = 0; i < static_cast<int>(Nf); i++) {
    FaceId fi = face_id(i);
    if (M.fbegin(fi) < 0)
      continue;
    triangulate_face(M, fi);
  }
  remove_orphan_vertices(M);
}

// Compact the vertex array by dropping vertices that no face references.
// Orphan vertices keep stale positions in the per-vertex data arrays, which
// lets downstream POV snapshots emit points with no surface neighborhood
// (and thus no meaningful medial axis). Remaps corner vertex ids and
// permutes every VERTEX datum so indices stay dense and consistent.
void remove_orphan_vertices(shell &M) {
  const index_t nv = static_cast<index_t>(M.vert_count());
  std::vector<index_t> old_to_new(nv, -1);
  std::vector<index_t> permute;
  permute.reserve(nv);
  index_t n = 0;
  for (index_t i = 0; i < nv; ++i) {
    if (M.vbegin(vert_id(i)) > -1) {
      old_to_new[i] = n;
      permute.push_back(i); // permute[new] = old
      ++n;
    }
  }
  if (n == nv) {
    return; // no orphans
  }

  std::vector<index_t> &cv = M.corners_vert();
  for (index_t &v : cv) {
    if (v >= 0 && v < nv) {
      v = old_to_new[v];
    }
  }

  for (datum_ptr &d : M.get_data()) {
    if (d && d->type() == prim_type::VERTEX) {
      d->permute(permute);
      d->resize(static_cast<size_t>(n));
    }
  }

  M.update_head();
  M.update_prev();
}

bool has_vert(shell &M, VertId vA, VertId vB) {

  bool hasB = false;
  M.for_each_vertex(vA, [&hasB, vB](CornerId ci, shell &M) {
    VertId vi = M.vert(M.next(ci));
    hasB |= int(vi == vB);
  });
  return hasB;
}

bool corner_in_ring(shell &M, VertId vA, CornerId cB) {

  bool hasB = false;
  M.for_each_vertex(vA, [&hasB, cB](CornerId ci, shell &M) {
    CornerId cBi = M.next(ci);
    hasB |= M.edge_equal(cBi, cB);
  });
  std::cout << std::endl;
  return hasB;
}

bool adjacent0(shell &M, CornerId cA0, CornerId cB0) {

  CornerId cA1 = M.other(cA0);
  CornerId cB1 = M.other(cB0);

  CornerId cA0n = M.next(cA0);
  CornerId cA1n = M.next(cA1);
  CornerId cB0n = M.next(cB0);
  CornerId cB1n = M.next(cB1);
  if (M.edge_equal(cA0n, cB0n))
    return true;
  if (M.edge_equal(cA0n, cB1n))
    return true;
  if (M.edge_equal(cA1n, cB0n))
    return true;
  if (M.edge_equal(cA1n, cB1n))
    return true;

  return false;
}

bool adjacent(shell &M, CornerId cA0, CornerId cB0) {

  CornerId cA1 = M.other(cA0);
  CornerId cB1 = M.other(cB0);
  if (M.vert(cA0) == M.vert(cB0) && //
      M.vert(cA1) == M.vert(cB1)) {
    return false;
  }
  bool a0ha0 = has_vert(M, M.vert(cA0), M.vert(cB0));
  bool a1ha1 = has_vert(M, M.vert(cA1), M.vert(cB1));
  if (a0ha0 || a1ha1)
    return true;

  return false;
}

bool share_faces(shell &M, CornerId cA0, CornerId cB0) {
  CornerId cA1 = M.other(cA0);
  CornerId cB1 = M.other(cB0);

  FaceId fA0 = M.face(cA0);
  FaceId fA1 = M.face(cA1);

  FaceId fB0 = M.face(cB0);
  FaceId fB1 = M.face(cB1);

  if (fA0 == fB0)
    return true;
  if (fA1 == fB1)
    return true;

  if (fA1 == fB0)
    return true;
  if (fA0 == fB1)
    return true;

  return false;
}

void weld_adajacent_edges(shell &M, CornerId cA0, CornerId cB0) {

  CornerId cA1 = M.other(cA0);
  CornerId cB1 = M.other(cB0);

  std::vector<CornerId> corners;
  M.for_each_vertex(M.vert(cA0), [&corners, cB0, cB1](CornerId ci, shell &M) {
    CornerId cBi = M.next(ci);
    if (M.edge_equal(cBi, cB0))
      corners.push_back(cBi);
    if (M.edge_equal(cBi, cB1))
      corners.push_back(cBi);
  });

  M.for_each_vertex(M.vert(cA1), [&corners, cB0, cB1](CornerId ci, shell &M) {
    CornerId cBi = M.next(ci);
    if (M.edge_equal(cBi, cB0))
      corners.push_back(cBi);
    if (M.edge_equal(cBi, cB1))
      corners.push_back(cBi);
  });

  for (CornerId c : corners) {
    if (M.next(c) > -1) {
      collapse_edge(M, c);
    }
  }
}

std::array<int, 4> merge_edge(shell &M, CornerId cA0, CornerId cB0,
                              VertId new_vert0 = vert_id(-1),
                              VertId new_vert1 = vert_id(-1)) {

  std::array<int, 4> out = {-1, -1, -1, -1};

  CornerId cA1 = M.other(cA0);
  CornerId cB1 = M.other(cB0);
  VertId vA0 = M.vert(cA0);
  VertId vA1 = M.vert(cA1);
  VertId vB0 = M.vert(cB0);
  VertId vB1 = M.vert(cB1);

  if (vA0 < 0 || vA1 < 0 || vB0 < 0 || vB1 < 0)
    return out;

  if (vA0 == vA1 || vB0 == vB1)
    return out;

  if (M.vsize(vA0) < 4 || M.vsize(vA1) < 4 || //
      M.vsize(vB0) < 4 || M.vsize(vB1) < 4) {
    return out;
  }

  if (adjacent(M, cA0, cB0)) {
    weld_adajacent_edges(M, cA0, cB0);
    return out;
  }

  M.swap_rows(cA1, cB1);
  M.set_vbegin(vA0, cA0);
  M.set_vbegin(vA1, cA1);
  out[0] = vA0;
  out[1] = vA1;

  if (vA0 == vB0) {
    VertId vN0 = new_vert0 < 0 ? M.insert_vertex() : new_vert0;

    out[2] = vN0;
    M.set_vbegin(vN0, cB0);

    M.vupdate(vA0);
    M.vupdate(vN0);
  } else {
    M.vupdate(vA0);
    M.remove_vertex(vB0);
  }

  if (vA1 == vB1) {
    VertId vN1 = new_vert1 < 0 ? M.insert_vertex() : new_vert1;

    out[3] = vN1;

    M.set_vbegin(vN1, cB1);

    M.vupdate(vA1);
    M.vupdate(vN1);

  } else {
    M.vupdate(vA1);

    M.remove_vertex(vB1);
  }
  return out;
}

template <int S>
std::vector<index_t> get_pack_permutation(std::vector<index_t> &indices) {
  std::vector<index_t> perm(indices.size());
  std::iota(perm.begin(), perm.end(), 0);

  int w = 0;

  for (int r = S; r < static_cast<int>(indices.size()); r += S) {
    if (indices[perm[r - S]] < 0 && indices[perm[r]] > -1 &&
        indices[perm[w]] < 0) {
      for (int i = 0; i < S; i++) {
        std::swap(perm[r + i], perm[w + i]);
      }
      w += S;
    } else if (indices[perm[w]] > -1) {
      w += S;
    }
  }
  return perm;
}

std::vector<index_t> inverse_permutation(const std::vector<index_t> &perm) {
  std::vector<index_t> iperm(perm.size(), -1);
  for (int i = 0; i < static_cast<int>(perm.size()); i++) {
    if (perm[i] > -1)
      iperm[perm[i]] = i;
  }
  return iperm;
}

void apply_permutation(const std::vector<index_t> &perm,
                       std::vector<index_t> &indices) {
  std::vector<index_t> n_indices(indices);
  for (int i = 0; i < static_cast<int>(indices.size()); i++) {
    n_indices[i] = indices[perm[i]];
  }
  indices = n_indices;
}

size_t calc_new_size(std::vector<index_t> &indices) {
  auto position = std::find(indices.begin(), indices.end(), -1);
  int index = static_cast<int>(position - indices.begin());
  return static_cast<size_t>(index);
}

void apply_inverse_permutation(const std::vector<index_t> &iperm,
                               std::vector<index_t> &indices) {
  for (int i = 0; i < static_cast<int>(indices.size()); i++) {
    if (indices[i] > -1)
      indices[i] = iperm[indices[i]];
  }
}

void pack(shell &M) {

  std::cout << "*--- packing ---*" << std::endl;

  auto debug = [](const std::vector<index_t> &indices, index_t s, index_t e,
                  const std::string &txt) {
    const index_t n = static_cast<index_t>(indices.size());
    if (n == 0) {
      std::cout << txt << ": (empty)" << std::endl;
      return;
    }
    if (s < 0)
      s = 0;
    if (e > n)
      e = n;
    if (s >= e) {
      std::cout << txt << ": (skip debug window; size=" << n << ")" << std::endl;
      return;
    }
    std::cout << txt << ": ";
    for (index_t i = s; i < e; i++) {
      if (indices[static_cast<size_t>(i)] < 0)
        std::cout << -1 << " ";
      else
        std::cout << indices[static_cast<size_t>(i)] << " ";
    }
    std::cout << std::endl;
  };

  std::vector<index_t> vperm = get_pack_permutation<1>(M.vert_begin());
  std::vector<index_t> viperm = inverse_permutation(vperm);

  apply_permutation(vperm, M.vert_begin());
  apply_inverse_permutation(viperm, M.corners_vert());
  size_t Nv = calc_new_size(M.vert_begin());
  M.vert_begin().resize(Nv);
  std::cout << "permute verts: " << Nv << std::endl;
  for (auto d : M.get_data()) {
    if (d->type() == asawa::VERTEX) {
      d->permute(vperm);
      d->resize(Nv);
    }
  }

  std::vector<index_t> fperm = get_pack_permutation<1>(M.face_begin());
  std::vector<index_t> fiperm = inverse_permutation(fperm);

  apply_permutation(fperm, M.face_begin());
  apply_inverse_permutation(fiperm, M.corners_face());
  size_t Nf = calc_new_size(M.face_begin());
  M.face_begin().resize(Nf);
  std::cout << "permute faces: " << Nf << std::endl;
  for (auto d : M.get_data()) {
    if (d->type() == FACE) {
      d->permute(fperm);
      d->resize(Nf);
    }
  }

#if 1

  std::vector<index_t> cperm = get_pack_permutation<2>(M.corners_next());
  std::vector<index_t> ciperm = inverse_permutation(cperm);
  // Debug window: was hardcoded 732–764 for huge meshes; clamp to valid range.
  const index_t ncorn = static_cast<index_t>(cperm.size());
  const index_t dbg_lo = ncorn > 32 ? ncorn - 32 : 0;
  const index_t dbg_hi = ncorn;
  debug(cperm, dbg_lo, dbg_hi, "cperm: ");
  debug(ciperm, dbg_lo, dbg_hi, "ciperm: ");

  debug(M.corners_face(), dbg_lo, dbg_hi, "corners_f, before: ");
  apply_permutation(cperm, M.corners_next());
  apply_permutation(cperm, M.corners_prev());
  apply_permutation(cperm, M.corners_face());
  apply_permutation(cperm, M.corners_vert());
  debug(M.corners_face(), dbg_lo, dbg_hi, "corners_f, after: ");

  apply_inverse_permutation(ciperm, M.corners_next());
  apply_inverse_permutation(ciperm, M.corners_prev());

  apply_inverse_permutation(ciperm, M.face_begin());
  apply_inverse_permutation(ciperm, M.vert_begin());
  size_t Nc = calc_new_size(M.corners_next());
  M.corners_next().resize(Nc);
  M.corners_prev().resize(Nc);
  M.corners_face().resize(Nc);
  M.corners_vert().resize(Nc);

  std::cout << "permute corners: " << Nc << std::endl;
  std::cout << cperm.size() << " " << ciperm.size() << std::endl;
  std::vector<index_t> cperm_half(cperm.size() / 2);

  for (size_t i = 0; i < cperm.size(); i += 2) {
    cperm_half[i / 2] = cperm[i] / 2;
  }

  for (auto d : M.get_data()) {

    if (d->type() == EDGE) {
      std::cout << "edge" << std::endl;
      d->permute(cperm_half);
      d->resize(Nc / 2);
    }

    if (d->type() == CORNER) {
      std::cout << "corner" << std::endl;
      d->permute(cperm);
      d->resize(Nc);
    }
  }
#endif

  std::cout << "*--- done ---*" << std::endl;
}

} // namespace shell
} // namespace asawa
} // namespace gaudi
#endif
