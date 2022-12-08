
#ifndef __ASAWA_PRIM_OPS__
#define __ASAWA_PRIM_OPS__

#include "m2_refactor.hpp"
#include <cassert>

namespace asawa {
typedef int index_t;
index_t split_edge(manifold &M, //
                   const index_t &corner_index) {
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
  M.vprint(v0);
  M.vprint(v1);
  M.fprint(f0);
  M.fprint(f1);

  M.remove_vertex(v0);
  M.remove_edge_pair(c0);

  M.link(c0p, c0n);
  M.link(c1p, c1n);
  M.set_vert(c1n, v1);

  M.set_face(c0n, f0);
  M.set_face(c1n, f1);

  std::cout << "v1: " << v1 << std::endl;
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

  M.fprint(nf);
  M.fprint(f0);
  M.vprint(v0);

  M.fupdate(nf);
  M.fupdate(f0);
  M.vupdate(v0);

  return c0i;
}

index_t remove_dangling(manifold &M, //
                        const index_t &c0, const index_t &c1) {
  std::cout << __func__ << std::endl;
  index_t c0p = M.prev(c0);
  index_t c0n = M.next(c0);
  index_t c1p = M.prev(c1);
  index_t c1n = M.next(c1);
  index_t f0 = M.face(c0);
  M.fprintv(f0);
  M.remove_vertex(M.vert(c0));
  M.remove_edge_pair(c0);
  M.link(c1p, c0n);
  M.set_face(c0n, f0);
  M.fprintv(f0);
  return f0;
}

index_t merge_face(manifold &M, //
                   const index_t &c0, const index_t &c1) {
  assert(c1 == M.other(c0));
  index_t c0p = M.prev(c0);
  index_t c0n = M.next(c0);
  index_t c1p = M.prev(c1);
  index_t c1n = M.next(c1);

  index_t f0 = M.face(c0);
  index_t f1 = M.face(c1);

  std::cout << "--bfoo--" << std::endl;
  std::cout << c0 << " " << c1 << std::endl;
  std::cout << f0 << " " << f1 << std::endl;

  M.vprint(M.vert(c0));
  M.vprint(M.vert(c1));
  std::cout << "--efoo--" << std::endl;

  if (f0 == f1) {
    // TODO: This isn't finished, yet

    if (c0 == c1n) {
      return remove_dangling(M, c0, c1);
    }
    if (c1 == c0n) {
      return remove_dangling(M, c1, c0);
    }
  }

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

  std::cout << __FUNCTION__ << " " << corner_index << std::endl;
  index_t c0 = split_edge(M, corner_index);

  index_t c1 = M.next(M.other(c0));
  split_face(M, c0, M.next(M.next(c0)));
  split_face(M, c1, M.next(M.next(c1)));

  return c0;
}

index_t collapse_edge(manifold &M, //
                      const index_t &corner_index) {
  std::cout << __FUNCTION__ << " " << corner_index << std::endl;

  index_t c0 = corner_index;
  index_t c1 = M.other(c0);

  index_t c0p = M.prev(c0);
  index_t c1p = M.prev(c1);
  std::cout << M.vert(c0) << " " << M.vert(c1) << std::endl;
  M.fprintv(M.face(c0));
  M.fprintv(M.face(c1));

  merge_face(M, c0p, M.other(c0p));
  merge_face(M, c1p, M.other(c1p));
  c0 = merge_edge(M, c0);

  return c0;
}

index_t flip_edge(manifold &M, //
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

  return c0i;
}

index_t remove_vertex(manifold &M, //
                      const index_t &v) {
  std::cout << " " << __FUNCTION__ << " " << v << std::endl;
  std::vector<index_t> corners;
  M.for_each_vertex(
      v, [&corners](index_t cid, manifold &m) { corners.push_back(cid); });
  index_t f = -1;
  for (index_t c : corners)
    f = merge_face(M, c, M.other(c));
  return f;
}

void triangulate_face(manifold &M, //
                      index_t i) {
  index_t c0 = M.fbegin(i);
  index_t c1 = M.next(M.next(M.fbegin(i)));
  index_t ce = M.fend(i);

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