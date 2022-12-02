
#ifndef __ASAWA_PRIM_OPS__
#define __ASAWA_PRIM_OPS__

#include "m2_refactor.hpp"
#include "manifold/asawa/subdivide.hpp"
#include <cassert>

namespace asawa {
index_t split_edge(manifold &M, //
                   const index_t &corner_index) {

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

  index_t vn = M.insert_vertex();
  index_t v1 = M.vert(c1);

  index_t c0i = M.insert_edge_pair();
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

index_t merge_edge(manifold &M, //
                   const index_t &corner_index) {

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
  index_t c0p = M.prev(c0);

  index_t c1n = M.next(c1);
  index_t c1p = M.prev(c1);

  index_t v0 = M.vert(c0);
  index_t v1 = M.vert(c1);
  index_t f0 = M.face(c0);
  index_t f1 = M.face(c1);

  M.remove_vertex(v0);
  M.remove_edge_pair(c0);

  M.link(c0p, c0n);
  M.link(c1p, c1n);
  M.set_vert(c1n, v1);

  M.set_face(c1n, f0);
  M.set_face(c0n, f1);

  M.vupdate(v1);

  return c1n;
}

index_t split_face(manifold &M, //
                   const index_t &c0, const index_t &c1) {

  index_t nf = M.insert_face();
  index_t f0 = M.face(c0);
  index_t v0 = M.vert(c0);
  index_t v1 = M.vert(c1);

  index_t c0p = M.prev(c0);
  index_t c1p = M.prev(c1);

  index_t c0i = M.insert_edge_pair();
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

  return c0i;
}

index_t merge_face(manifold &M, //
                   const index_t &c0, const index_t &c1) {
  assert(c1 == M.other(c0));
  std::cout << c0 << " " << c1 << std::endl;
  index_t c0p = M.prev(c0);
  index_t c0n = M.next(c0);
  index_t c1p = M.prev(c1);
  index_t c1n = M.next(c1);

  index_t f0 = M.face(c0);
  index_t f1 = M.face(c1);
  index_t v0 = M.vert(c0);
  index_t v1 = M.vert(c1);

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

index_t subdivide_edge(manifold &M, //
                       const index_t &corner_index) {

  index_t c0 = split_edge(M, corner_index);
  index_t c1 = M.next(M.other(c0));
  split_face(M, c0, M.next(M.next(c0)));
  split_face(M, c1, M.next(M.next(c1)));

  return c0;
}

index_t collapse_edge(manifold &M, //
                      const index_t &corner_index) {

  index_t c0 = corner_index;
  index_t c1 = M.other(c0);
  index_t c0p = M.prev(c0);
  index_t c1p = M.prev(c1);
  merge_face(M, c0p, M.other(c0p));
  merge_face(M, c1p, M.other(c1p));
  c0 = merge_edge(M, c0);
  return c0;
}

void triangulate_face(manifold &M, //
                      index_t i) {
  index_t c0 = M.fbegin(i);
  index_t c1 = M.next(M.next(M.fbegin(i)));
  index_t ce = M.fend(i);

  std::cout << "f" << i << ": ";
  M.for_each_face(i, [](index_t cid, manifold &m) { std::cout << cid << " "; });
  std::cout << " end" << std::endl;

  while (c1 != ce) {
    index_t c1t = c1;
    c1 = M.next(c1);
    index_t nc = split_face(M, c0, c1t);
  }
}

void triangulate(manifold &M) {
  size_t Nf = M.face_count();
  for (int i = 0; i < Nf; i++) {
    if (M.face(i) < 0)
      continue;
    triangulate_face(M, i);
  }
}

} // namespace asawa
#endif