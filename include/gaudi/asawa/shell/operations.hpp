
#ifndef __ASAWA_PRIM_OPS__
#define __ASAWA_PRIM_OPS__

//#include "datums.hpp"
#include "shell.hpp"

#include <array>
#include <cassert>
#include <ostream>
#include <vector>

namespace gaudi {
namespace asawa {
namespace shell {
typedef int index_t;
index_t split_edge(shell &M,                       //
                   const index_t &corner_index,    //
                   const index_t &new_vertex = -1, //
                   const index_t &new_corner = -1) {
  // std::cout << __FUNCTION__ << " " << corner_index << std::endl;

  // originial topology:
  //    c1_<--c1p
  //   /     /
  //  c0 -->c0n

  // target topology:
  //   c1_<--c1i<--c1p
  //  /     /
  // c0_-->c0i-->c0n

  index_t c0 = corner_index;
  index_t c1 = M.other(c0);

  index_t c0n = M.next(c0);
  index_t c1p = M.prev(c1);

  index_t vn = new_vertex < 0 ? M.insert_vertex() : new_vertex;
  index_t v1 = M.vert(c1);

  index_t c0i = new_corner < 0 ? M.insert_edge_pair() : new_corner;
  index_t c1i = M.other(c0i);

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

int count_cycle(shell &M, index_t corner) {
  index_t c0 = corner;
  index_t c1 = M.other(c0);
  index_t v0 = M.vert(c0);
  index_t v1 = M.vert(c1);
  int t = 0;
  M.for_each_vertex(v0, [&t, v1](index_t ci, shell &M) {
    index_t vi = M.vert(M.next(ci));
    // std::cout << vi << " - " << v1 << std::endl;
    t += int(vi == v1);
  });
  return t;
}

index_t merge_edge(shell &M, //
                   const index_t &corner_index) {
  // std::cout << __FUNCTION__ << " " << corner_index << std::endl;
  //  originial topology:
  //     c1_<--c1p
  //    /     /
  //   c0 -->c0n

  // target topology:
  //   c1_<--c1i<--c1p
  //  /     /
  // c0_-->c0i-->c0n

  index_t c0 = corner_index;
  index_t c1 = M.other(c0);
  index_t c0n = M.next(c0);
  index_t c0p = M.prev(c0);

  index_t c1n = M.next(c1);
  index_t c1p = M.prev(c1);

  index_t v0 = M.vert(c0);
  index_t v1 = M.vert(c1);
  index_t f0 = M.face(c0);
  index_t f1 = M.face(c1);

  if (v0 != v1)
    M.remove_vertex(v0);
  M.remove_edge_pair(c0);

  M.link(c0p, c0n);
  M.link(c1p, c1n);

  M.set_vert(c1n, v1);

  M.set_face(c0n, f0);
  M.set_face(c1n, f1);
  //  assert(M.fsize(f0) == 3);
  //  assert(M.fsize(f1) == 3);
  M.vupdate(v1);

  return c1n;
}

index_t split_face(shell &M, //
                   const index_t &c0, const index_t &c1,
                   const index_t new_corner = -1, //
                   const index_t new_face = -1) {

  index_t nf = new_face < 0 ? M.insert_face() : new_face;
  index_t f0 = M.face(c0);
  index_t v0 = M.vert(c0);
  index_t v1 = M.vert(c1);

  index_t c0p = M.prev(c0);
  index_t c1p = M.prev(c1);

  index_t c0i = new_corner < 0 ? M.insert_edge_pair() : new_corner;
  index_t c1i = M.other(c0i);

  M.link(c0p, c0i);
  M.link(c0i, c1);

  M.link(c1p, c1i);
  M.link(c1i, c0);
  // M.link(c1, c1i);

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

index_t remove_dangling(shell &M, //
                        const index_t &c0, const index_t &c1) {
  // std::cout << __func__ << std::endl;
  // std::cout << __PRETTY_FUNCTION__ << std::endl;
  index_t c0p = M.prev(c0);
  index_t c0n = M.next(c0);
  index_t c1p = M.prev(c1);
  index_t c1n = M.next(c1);
  index_t f0 = M.face(c0);
  index_t v0 = M.vert(c0);
  index_t v1 = M.vert(c1);

  M.remove_vertex(v0);
  M.remove_edge_pair(c0);
  M.set_vert(c0n, v1);
  M.link(c1p, c0n);
  M.set_face(c0n, f0);

  M.vupdate(v1);

  // M.uber_assert();
  return f0;
}

index_t remove_cruft(shell &M, //
                     const index_t &c0, const index_t &c1) {
  // std::cout << __PRETTY_FUNCTION__ << std::endl;

  index_t f0 = M.face(c0);
  index_t f1 = M.face(c1);
  index_t v0 = M.vert(c0);
  index_t v1 = M.vert(c1);

  M.remove_vertex(v0);
  M.remove_vertex(v1);
  M.remove_face(f0);
  M.remove_face(f1);

  M.remove_edge_pair(c0);
  return -1;
}

index_t merge_face(shell &M, //
                   const index_t &c0, const index_t &c1) {
  assert(c1 == M.other(c0));
  index_t c0p = M.prev(c0);
  index_t c0n = M.next(c0);
  index_t c1p = M.prev(c1);
  index_t c1n = M.next(c1);

  index_t f0 = M.face(c0);
  index_t f1 = M.face(c1);

  index_t v0 = M.vert(c0);
  index_t v1 = M.vert(c1);

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

  // index_t v0 = M.vert(c0);
  // index_t v1 = M.vert(c1);

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

index_t subdivide_edge(shell &M,                        //
                       const index_t &corner_index,     //
                       const index_t &new_vert = -1,    //
                       const index_t &new_corner0 = -1, //
                       const index_t &new_corner1 = -1, //
                       const index_t &new_corner2 = -1, //
                       const index_t &new_face0 = -1,   //
                       const index_t &new_face1 = -1) {

  // std::cout << __FUNCTION__ << " " << corner_index << std::endl;
  index_t c0 = corner_index;
  index_t c1 = M.next(M.other(c0));

  index_t c0p = M.prev(c0);
  index_t c1p = M.prev(c1);

  /*
  if (M.vsize(M.vert(c0p)) > 15) {
    return c0;
  }
  if (M.vsize(M.vert(c1p)) > 15)
    return c0;
  */

  c0 = split_edge(M, c0, new_vert, new_corner0);
  c1 = M.next(M.other(c0));

  split_face(M, c0, M.next(M.next(c0)), new_corner1, new_face0);
  split_face(M, c1, M.next(M.next(c1)), new_corner2, new_face1);

  return c0;
}

index_t collapse_edge(shell &M, //
                      const index_t &corner_index, bool degenerate = false) {
  // std::cout << __FUNCTION__ << " " << corner_index << std::endl;
  // M.cprint(corner_index);

  index_t c0 = corner_index;
  index_t c1 = M.other(c0);

  // M.cprint(c0);
  //  M.cprint(c1);

  index_t c0p = M.prev(c0);
  index_t c1p = M.prev(c1);

  index_t v0 = M.vert(c0);
  index_t v1 = M.vert(c1);

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

  // M.vprintv(v0);
  // M.vprintv(v1);

  merge_face(M, c0p, M.other(c0p));
  merge_face(M, c1p, M.other(c1p));

  c0 = merge_edge(M, c0);

  assert(M.vert(c0) != M.vert(M.next(c0)));

  // M.vprintv(v1);
  // if (!degenerate)
  //  M.uber_assert();
  return c0;
}

index_t flip_edge(shell &M, //
                  const index_t &corner_index) {
  // this is kind of like an edge delete/insert packaged up into a primitive
  // function
  index_t c0i = corner_index;
  index_t c0p = M.prev(c0i);
  index_t c0pp = M.prev(c0p);
  index_t c0n = M.next(c0i);

  index_t c1i = M.other(c0i);
  index_t c1p = M.prev(c1i);
  index_t c1pp = M.prev(c1p);
  index_t c1n = M.next(c1i);

  index_t f0 = M.face(c0i);
  index_t f1 = M.face(c1i);

  index_t v0 = M.vert(c0i);
  index_t v1 = M.vert(c0p);
  index_t v2 = M.vert(c1i);
  index_t v3 = M.vert(c1p);

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

  if (M.fsize(f0) != 3) {
    // std::cout << "fsize 0: " << c0i << " " << M.fsize(f0) << std::endl;
    // std::cout << "fsize 1: " << c1i << " " << M.fsize(f1) << std::endl;
    M.fprintv(f0);
  }
  assert(M.fsize(f0) == 3);
  if (M.fsize(f1) != 3) {
    // std::cout << "fsize 1: " << c1i << " " << M.fsize(f1) << std::endl;
    // std::cout << "fsize 0: " << c0i << " " << M.fsize(f0) << std::endl;
    M.fprintv(f1);
  }
  assert(M.fsize(f1) == 3);

  assert(M.vert(c0i) != M.vert(c1i));
  // M.uber_assert();

  return c0i;
}

index_t remove_vertex(shell &M, //
                      const index_t &v) {

  // std::cout << __FUNCTION__ << " " << v << std::endl;
  // std::cout << "  vsize: " << M.vsize(v) << std::endl;

  std::vector<index_t> corners;
  M.for_each_vertex(v, [&corners](index_t cid, shell &m) {
    //    std::cout << m.vert(cid) << " ";
    // m.cprint(cid);
    corners.push_back(cid);
  });
  //  std::cout << std::endl;
  index_t f = -1;
  for (index_t c : corners)
    f = merge_face(M, c, M.other(c));

  // M.uber_assert();
  return f;
}

void triangulate_face(shell &M, //
                      index_t i) {
  index_t c0 = M.fbegin(i);

  bool splitting = true;
  while (splitting) {

    index_t c0n = M.next(c0);
    index_t c0nn = M.next(c0n);
    index_t c0p = M.prev(c0);

    splitting = c0nn != c0p;
    if (splitting) {
      c0 = split_face(M, c0n, c0p);
      c0 = M.other(c0);
    }
  }
}

void triangulate(shell &M) {
  size_t Nf = M.face_count();
  for (int i = 0; i < Nf; i++) {
    if (M.face(i) < 0)
      continue;
    triangulate_face(M, i);
  }
}

bool has_vert(shell &M, index_t vA, index_t vB) {

  bool hasB = false;
  M.for_each_vertex(vA, [&hasB, vB](index_t ci, shell &M) {
    index_t vi = M.vert(M.next(ci));
    hasB |= int(vi == vB);
  });
  return hasB;
}

bool corner_in_ring(shell &M, index_t vA, index_t cB) {

  bool hasB = false;
  M.for_each_vertex(vA, [&hasB, cB](index_t ci, shell &M) {
    index_t cBi = M.next(ci);
    //std::cout << cBi << " " << cB << " - ";
    hasB |= M.edge_equal(cBi, cB);
  });
  std::cout << std::endl;
  return hasB;
}

bool adjacent0(shell &M, index_t cA0, index_t cB0) {

  index_t cA1 = M.other(cA0);
  index_t cB1 = M.other(cB0);

  index_t cA0n = M.next(cA0);
  index_t cA1n = M.next(cA1);
  index_t cB0n = M.next(cB0);
  index_t cB1n = M.next(cB1);
  //std::cout << cA0n << " " << cA1n << " - " << cB0n << " " << cB1n << std::endl;
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

bool adjacent(shell &M, index_t cA0, index_t cB0) {

  index_t cA1 = M.other(cA0);
  index_t cB1 = M.other(cB0);
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

bool share_faces(shell &M, index_t cA0, index_t cB0) {
  index_t cA1 = M.other(cA0);
  index_t cB1 = M.other(cB0);

  index_t fA0 = M.face(cA0);
  index_t fA1 = M.face(cA1);

  index_t fB0 = M.face(cB0);
  index_t fB1 = M.face(cB1);

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

void weld_adajacent_edges(shell &M,    //
                          index_t cA0, //
                          index_t cB0) {

  index_t cA1 = M.other(cA0);
  index_t cB1 = M.other(cB0);

  std::vector<index_t> corners;
  M.for_each_vertex(M.vert(cA0), [&corners, cB0, cB1](index_t ci, shell &M) {
    index_t cBi = M.next(ci);
    if (M.edge_equal(cBi, cB0))
      corners.push_back(cBi);
    if (M.edge_equal(cBi, cB1))
      corners.push_back(cBi);
  });

  M.for_each_vertex(M.vert(cA1), [&corners, cB0, cB1](index_t ci, shell &M) {
    index_t cBi = M.next(ci);
    if (M.edge_equal(cBi, cB0))
      corners.push_back(cBi);
    if (M.edge_equal(cBi, cB1))
      corners.push_back(cBi);
  });

  for (auto c : corners) {
    if (M.next(c) > -1) {
      collapse_edge(M, c);
    }
  }
}

std::array<index_t, 4> merge_edge(shell &M,               //
                                  index_t cA0,            //
                                  index_t cB0,            //
                                  index_t new_vert0 = -1, //
                                  index_t new_vert1 = -1) {
  // alignEdges(eA, eB); function of moving mesh

  std::array<index_t, 4> out = {-1, -1, -1, -1};

  // if (share_faces(M, cA0, cB0))
  //   return out;

  index_t cA1 = M.other(cA0);
  index_t cB1 = M.other(cB0);
  index_t vA0 = M.vert(cA0);
  index_t vA1 = M.vert(cA1);
  index_t vB0 = M.vert(cB0);
  index_t vB1 = M.vert(cB1);

  if (vA0 < 0 || vA1 < 0 || vB0 < 0 || vB1 < 0)
    return out;

  /*
  if (count_cycle(M, cA0) > 1)
    return out;
  if (count_cycle(M, cB0) > 1)
    return out;
*/
  if (vA0 == vA1 || vB0 == vB1)
    return out;

  if (M.vsize(vA0) < 4 || M.vsize(vA1) < 4 || //
      M.vsize(vB0) < 4 || M.vsize(vB1) < 4) {
    return out;
  }

  if (adjacent(M, cA0, cB0)) {
    weld_adajacent_edges(M, cA0, cB0);
    // std::cout << "adjacent" << std::endl;
    return out;
  }

  M.swap_rows(cA1, cB1);
  M.set_vbegin(vA0, cA0);
  M.set_vbegin(vA1, cA1);
  out[0] = vA0;
  out[1] = vA1;

  if (vA0 == vB0) {
    // std::cout << " A.0 " << std::endl;s
    index_t vN0 = new_vert0 < 0 ? M.insert_vertex() : new_vert0;

    out[2] = vN0;
    M.set_vbegin(vN0, cB0);
    // M.set_vbegin(vN0, M.next(cA1));

    M.vupdate(vA0);
    M.vupdate(vN0);

    // M.vprintv(vA0);
    // M.vprintv(vN0);
    // std::cout << "vs: " << M.vsize(vN0) << " " << M.vsize(vA0) << std::endl;
  } else {
    // std::cout << "here A1!" << std::endl;
    M.vupdate(vA0);
    // M.vprintv(vA0);
    M.remove_vertex(vB0);
  }

  if (vA1 == vB1) {
    // std::cout << "here B0!" << std::endl;
    index_t vN1 = new_vert1 < 0 ? M.insert_vertex() : new_vert1;

    out[3] = vN1;

    M.set_vbegin(vN1, cB1);
    // M.set_vbegin(vN1, M.next(cA0));

    M.vupdate(vA1);
    M.vupdate(vN1);

    // M.vprintv(vA1);
    // M.vprintv(vN1);
    // std::cout << "vs: " << M.vsize(vN1) << " " << M.vsize(vA1) << std::endl;

  } else {
    // std::cout << "here B1!" << std::endl;
    M.vupdate(vA1);

    // M.vprintv(vA1);
    M.remove_vertex(vB1);
  }
  /*
  std::vector<index_t> corners;
  M.for_each_vertex(vA0, [&corners](index_t c0, shell &M) {
    index_t c1 = M.other(c0);
    if (M.vert(c0) == M.vert(c1))
      corners.push_back(c0);
  });

  M.for_each_vertex(vA1, [&corners](index_t c0, shell &M) {
    index_t c1 = M.other(c0);
    if (M.vert(c0) == M.vert(c1))
      corners.push_back(c0);
  });

  for (index_t c : corners) {
    if (M.next(c) < 0)
      continue;
    collapse_edge(M, c, true);
  }

  M.uber_assert();
*/
  return out;
}

std::vector<index_t> get_pack_permutation(std::vector<index_t> &indices) {
  std::vector<index_t> perm(indices.size());
  std::iota(perm.begin(), perm.end(), 0);

  int w = 0;
  //-----------xxxxxxx----x------
  //           |      |

  for (int r = 1; r < indices.size(); r++) {
    if (indices[perm[r - 1]] < 0 && indices[perm[r]] > -1 &&
        indices[perm[w]] < 0) {
      std::swap(perm[r], perm[w]);
      w++;
    } else if (indices[perm[w]] > -1) {
      w++;
    }
  }
  return perm;
}

std::vector<index_t> inverse_permutation(const std::vector<index_t> &perm) {
  std::vector<index_t> iperm(perm.size(), -1);
  for (int i = 0; i < perm.size(); i++) {
    if (perm[i] > -1)
      iperm[perm[i]] = i;
  }
  return iperm;
}

void apply_permutation(const std::vector<index_t> &perm,
                       std::vector<index_t> &indices) {
  std::vector<index_t> n_indices(indices);
  for (int i = 0; i < indices.size(); i++) {
    n_indices[i] = indices[perm[i]];
  }
  indices = n_indices;
}

size_t calc_new_size(std::vector<index_t> &indices) {
  std::vector<index_t>::iterator position =
      std::find(indices.begin(), indices.end(), -1);
  int index = position - indices.begin();
  return index;
}

void apply_inverse_permutation(const std::vector<index_t> &iperm,
                               std::vector<index_t> &indices) {
  for (int i = 0; i < indices.size(); i++) {
    if (indices[i] > -1)
      indices[i] = iperm[indices[i]];
  }
}

void pack(shell &M) {

  index_t s = 0;
  index_t e = s + 16;

  std::cout << "*--- packing ---*" << std::endl;

  auto debug = [s, e](const std::vector<index_t> indices, std::string txt) {
    std::cout << txt << ": ";
    for (int i = s; i < e /*__corners_next.size()*/; i++) {
      if (indices[i] < 0)
        std::cout << -1 << " ";
      else
        std::cout << indices[i] << " ";
    }
    std::cout << std::endl;
  };

  std::vector<index_t> vperm = get_pack_permutation(M.vert_begin());
  std::vector<index_t> viperm = inverse_permutation(vperm);

  apply_permutation(vperm, M.vert_begin());
  apply_inverse_permutation(viperm, M.corners_vert());
  size_t Nv = calc_new_size(M.vert_begin());
  M.vert_begin().resize(Nv);

  for (auto d : M.get_data()) {
    if (d->type() == asawa::VERTEX) {
      d->permute(vperm);
      d->resize(Nv);
    }
  }

  // debug(M.vert_begin(), std::string(" after v"));
  // debug(M.corners_vert(), std::string(" after c"));

  std::vector<index_t> fperm = get_pack_permutation(M.face_begin());
  std::vector<index_t> fiperm = inverse_permutation(fperm);

  apply_permutation(fperm, M.face_begin());
  apply_inverse_permutation(fiperm, M.corners_face());
  size_t Nf = calc_new_size(M.face_begin());
  M.face_begin().resize(Nf);

  for (auto d : M.get_data()) {
    if (d->type() == FACE) {
      d->permute(fperm);
      d->resize(Nf);
    }
  }

  // std::cout << std::flush;
}

} // namespace shell
} // namespace asawa
} // namespace gaudi
#endif