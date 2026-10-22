#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#if defined(__GNUC__) || defined(__clang__)
#include <cxxabi.h>
#endif
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory.h>
#include <string>
#include <numeric>
#include <ostream>
#include <set>
#include <stdio.h>
#include <type_traits>
#include <vector>

#include "gaudi/asawa/geometry_channels.hpp"
#include "gaudi/asawa/shell/shell_id.hpp"
#include "gaudi/exception_diagnostics.hpp"

#ifndef __ASAWA_SHELL__
#define __ASAWA_SHELL__

#ifdef _MSC_VER
#define __PRETTY_FUNCTION__ __FUNCSIG__
#endif

namespace gaudi {
namespace asawa {

inline std::ofstream &shell_link_log() {
  static std::ofstream out;
  if (!out.is_open()) {
    std::error_code ec;
    std::filesystem::create_directories("dump", ec);
    out.open("dump/shell_link.log", std::ios::out | std::ios::trunc);
    out << std::unitbuf;
  }
  return out;
}

class datum;

namespace shell {

typedef int index_t;

typedef std::shared_ptr<datum> datum_ptr;

class shell {
public:
  typedef std::shared_ptr<shell> ptr;

  static ptr create(const std::vector<index_t> &corners_next,
                    const std::vector<index_t> &corners_vert,
                    const std::vector<index_t> &corners_face) {

    return std::make_shared<shell>(corners_next, corners_vert, corners_face);
  }

  shell(const std::vector<index_t> &corners_next,
        const std::vector<index_t> &corners_vert,
        const std::vector<index_t> &corners_face)
      : __corners_next(corners_next), __corners_vert(corners_vert),
        __corners_face(corners_face) {
    this->update_head();
    this->update_prev();
  }

  void update_head() {
    auto update = [](std::vector<index_t> &p_beg,
                     const std::vector<index_t> &c_ptr,
                     const char *label) {
      int size = *std::max_element(c_ptr.begin(), c_ptr.end());

      p_beg = std::vector<index_t>(size + 1, -1);

      for (int i = 0; i < c_ptr.size(); i++) {
        if (c_ptr[i] < 0) {
          std::cerr << "[shell update_head] WARNING: " << label
                    << " has invalid index -1 at corner " << i << "\n";
          continue;
        }
        p_beg[c_ptr[i]] = i;
      }
    };

    update(__face_begin, __corners_face, "corners_face");
    update(__vert_begin, __corners_vert, "corners_vert");
  }

  void update_prev() {
    // can't use the iterator because the iterator needs prev.
    __corners_prev = std::vector<index_t>(__corners_next.size(), -1);
    for (int i = 0; i < __face_begin.size(); i++) {
      int j_end = __face_begin[i];
      int j0 = __corners_next[j_end];
      bool it = true;
      while (it) {
        it = j0 != j_end;
        require_corner_ref("update_prev", __corners_next[j0], j0);
        __corners_prev[__corners_next[j0]] = j0;
        j0 = __corners_next[j0];
      }
    }
  }

  // Defined in datums.hpp once datum is complete (stamps slot_id).
  index_t insert_datum(datum_ptr datum);

  datum_ptr &get_datum(index_t i) { return __channels.get(i); }
  const datum_ptr &const_get_datum(index_t i) const {
    return __channels.get(i);
  }

  std::vector<datum_ptr> &get_data() { return __channels.pool(); }
  const std::vector<datum_ptr> &get_data() const { return __channels.pool(); }
  GeometryChannels &channels() { return __channels; }
  const GeometryChannels &channels() const { return __channels; }

  size_t vert_count() const { return __vert_begin.size(); }

  /// Every slot \c 0 .. vert_count()-1 is an active vertex (\c vbegin valid).
  /// Per-vertex matrices (Laplacian, etc.) assume vertex id ≡ index into \c x[v].
  /// Linear scan, \f$O(N)\f$; a cached flag could follow if hot paths need it.
  bool verts_are_dense_packed() const {
    for (int i = 0; i < static_cast<int>(vert_count()); ++i) {
      if (static_cast<int>(vbegin(vert_id(i))) < 0)
        return false;
    }
    return true;
  }

  size_t face_count() const { return __face_begin.size(); }
  size_t edge_count() const { return __corners_next.size() / 2; }
  size_t corner_count() const { return __corners_next.size(); }

  /// A corner `c` is the half-edge from `vert(c)` to `vert(next(c))` in
  /// `face(c)`. `other(c)` is its twin in the same slot pair `{2k, 2k + 1}`,
  /// running the other way, so `vert(other(c)) == vert(next(c))`.
  CornerId other(CornerId id) const {
    return corner_id(2 * (id / 2) + (id + 1) % 2);
  }
  bool edge_equal(CornerId c0, CornerId c1) {
    return c0 / 2 == c1 / 2;
  }
  /// The outgoing corners at `vert(c)` in `for_each_vertex` order:
  /// `vnext(c) = other(prev(c))`. `face(c)` lies between spoke `c` and spoke
  /// `vnext(c)`. Needs `prev`.
  CornerId vnext(CornerId id) const { return other(prev(id)); }
  /// Inverse of `vnext`: `vprev(c) = next(other(c))`. Needs only `next` and
  /// `other`, so it still walks a ring while `prev` links are stale.
  CornerId vprev(CornerId id) const { return next(other(id)); }

  /// The corner after / before `c` in `face(c)`.
  CornerId next(CornerId id) const { return corner_id(__corners_next[id]); }
  CornerId prev(CornerId id) const { return corner_id(__corners_prev[id]); }

  VertId vert(CornerId id) const { return vert_id(__corners_vert[id]); }
  FaceId face(CornerId id) const { return face_id(__corners_face[id]); }

  /// Two triangles sharing the undirected edge of \p edge_c0 (0–2 faces).
  std::vector<FaceId> dihedral_face_ids(CornerId edge_c0) const {
    std::vector<FaceId> out;
    CornerId c1 = other(edge_c0);
    if (next(edge_c0) < 0)
      return out;
    FaceId f0 = face(edge_c0);
    if (static_cast<int>(f0) >= 0)
      out.push_back(f0);
    if (next(c1) >= 0) {
      FaceId f1 = face(c1);
      if (static_cast<int>(f1) >= 0 && f1 != f0)
        out.push_back(f1);
    }
    return out;
  }

  /// Butterfly stencil (see \c duchamp/modules/cross.hpp): the two
  /// faces on \p edge_c0 plus up to four faces across \c prev/next wings (≤6 unique faces).
  std::vector<FaceId> butterfly_face_ids(CornerId edge_c0) const {
    std::set<FaceId> seen;
    auto add_f = [&seen](FaceId f) {
      if (static_cast<int>(f) >= 0)
        seen.insert(f);
    };
    CornerId c0 = edge_c0;
    CornerId c1 = other(c0);
    add_f(face(c0));
    if (next(c1) >= 0)
      add_f(face(c1));
    if (next(c0) < 0)
      return std::vector<FaceId>(seen.begin(), seen.end());
    CornerId ci = prev(c0);
    CornerId cj = next(c0);
    CornerId ck = prev(c1);
    CornerId cl = next(c1);
    auto add_across = [this, &add_f](CornerId c) {
      if (next(c) < 0)
        return;
      CornerId t = other(c);
      if (next(t) >= 0)
        add_f(face(t));
    };
    add_across(ci);
    add_across(cj);
    if (next(c1) >= 0) {
      add_across(ck);
      add_across(cl);
    }
    return std::vector<FaceId>(seen.begin(), seen.end());
  }

  /// Union of face stars of the three vertices of triangle \p f (face 1-ring).
  void face_one_ring_face_ids(FaceId f, std::vector<FaceId> &out) const {
    std::set<FaceId> seen;
    const_for_each_face(f, [&seen, this](CornerId c0, const shell &M) {
      VertId v = M.vert(c0);
      M.const_for_each_vertex(v, [&seen, &M](CornerId ci, const shell &Mm) {
        FaceId fi = Mm.face(ci);
        if (static_cast<int>(fi) >= 0)
          seen.insert(fi);
      });
    });
    out.assign(seen.begin(), seen.end());
  }

  std::vector<FaceId> face_one_ring_face_ids(FaceId f) const {
    std::vector<FaceId> out;
    face_one_ring_face_ids(f, out);
    return out;
  }

  /// Faces within graph distance ≤ 2 under face-1-ring adjacency.
  void face_two_ring_face_ids(FaceId f, std::vector<FaceId> &out) const {
    std::set<FaceId> seen;
    std::vector<FaceId> ring1 = face_one_ring_face_ids(f);
    for (FaceId f1 : ring1) {
      if (static_cast<int>(f1) >= 0)
        seen.insert(f1);
      std::vector<FaceId> ring2;
      face_one_ring_face_ids(f1, ring2);
      for (FaceId f2 : ring2) {
        if (static_cast<int>(f2) >= 0)
          seen.insert(f2);
      }
    }
    out.assign(seen.begin(), seen.end());
  }

  std::vector<FaceId> face_two_ring_face_ids(FaceId f) const {
    std::vector<FaceId> out;
    face_two_ring_face_ids(f, out);
    return out;
  }

  /// First / last corner `for_each_face` visits; -1 marks a dead face.
  CornerId fbegin(FaceId id) const { return corner_id(__face_begin[id]); }
  CornerId fend(FaceId id) const { return prev(corner_id(__face_begin[id])); }

  /// First / last corner `for_each_vertex` visits (`vnext(vend) == vbegin`);
  /// -1 marks a dead vertex.
  CornerId vbegin(VertId id) const { return corner_id(__vert_begin[id]); }
  CornerId vend(VertId id) const { return vprev(vbegin(id)); }

  void set_vbegin(VertId id, CornerId c) { __vert_begin[id] = c; }
  void set_fbegin(FaceId id, CornerId c) { __face_begin[id] = c; }

  void set_next(CornerId id, CornerId c) { __corners_next[id] = c; }

  void report_link_failure(const char *where, const std::string &msg) const {
    std::ofstream &out = shell_link_log();
    out << "shell link [" << where << "] " << msg << "\n";
    gaudi::print_native_backtrace(out);
    out.flush();
    std::cerr << "shell link [" << where << "] " << msg
              << "\n  wrote dump/shell_link.log\n"
              << std::flush;
    std::abort();
  }

  void require_corner_ref(const char *what, index_t slot, index_t value) const {
    const index_t n = static_cast<index_t>(corner_count());
    const bool slot_ok = slot >= 0 && slot < n;
    const bool value_ok = value == -1 || (value >= 0 && value < n);
    if (slot_ok && value_ok)
      return;
    report_link_failure(
        what, "store slot=" + std::to_string(slot) +
                  " value=" + std::to_string(value) + " n=" + std::to_string(n));
  }

  void set_prev(CornerId id, CornerId c) {
    require_corner_ref("set_prev", static_cast<index_t>(id),
                       static_cast<index_t>(c));
    __corners_prev[id] = c;
  }

  // Between topology phases. A live corner's prev must name another live
  // corner. A prev of -1, past the array, or aimed at a deleted corner is the
  // sew that later walks blow up on. The phase name is the operation.
  void audit_links(const char *where) const {
    const int n = static_cast<int>(corner_count());
    auto live = [&](int c) {
      return c >= 0 && c < n && __corners_next[static_cast<size_t>(c)] >= 0;
    };
    for (int i = 0; i < n; ++i) {
      const int nx = __corners_next[static_cast<size_t>(i)];
      if (nx < 0)
        continue;
      const int pv = __corners_prev[static_cast<size_t>(i)];
      const int v = __corners_vert[static_cast<size_t>(i)];
      std::string why;
      if (nx >= n)
        why = "next out of range";
      else if (pv < 0 || pv >= n)
        why = "prev out of range";
      else if (__corners_next[static_cast<size_t>(pv)] < 0)
        why = "prev names a deleted corner";
      if (why.empty())
        continue;
      report_link_failure(
          where, why + " corner=" + std::to_string(i) + " next=" +
                     std::to_string(nx) + " prev=" + std::to_string(pv) +
                     " vert=" + std::to_string(v) + " n=" + std::to_string(n));
    }
    audit_vbegins(where);
    std::ofstream &out = shell_link_log();
    out << "shell link [" << where << "] ok corners=" << n
        << " verts=" << vert_count() << "\n";
  }

  // vbegin of a live vertex must head a live corner. Cheap enough to call
  // after each topology primitive.
  void audit_vbegins(const char *where) const {
    const int n = static_cast<int>(corner_count());
    auto live = [&](int c) {
      return c >= 0 && c < n && __corners_next[static_cast<size_t>(c)] >= 0;
    };
    for (int i = 0; i < n; ++i) {
      if (!live(i))
        continue;
      const int o = static_cast<int>(other(corner_id(i)));
      const int nx = __corners_next[static_cast<size_t>(i)];
      if (!live(o) || !live(nx))
        continue;
      const int nv = __corners_vert[static_cast<size_t>(nx)];
      const int ov = __corners_vert[static_cast<size_t>(o)];
      if (nv == ov)
        continue;
      report_link_failure(
          where, "next vert != twin vert c=" + std::to_string(i) +
                     " next=" + std::to_string(nx) + " next_vert=" +
                     std::to_string(nv) + " twin=" + std::to_string(o) +
                     " twin_vert=" + std::to_string(ov) + " n=" +
                     std::to_string(n));
    }
    for (int vi = 0; vi < static_cast<int>(vert_count()); ++vi) {
      const int vb = __vert_begin[static_cast<size_t>(vi)];
      if (vb < 0)
        continue;
      if (live(vb) && __corners_vert[static_cast<size_t>(vb)] != vi) {
        report_link_failure(
            where, "vbegin vert mismatch v=" + std::to_string(vi) +
                       " vbegin=" + std::to_string(vb) + " vert=" +
                       std::to_string(__corners_vert[static_cast<size_t>(vb)]) +
                       " n=" + std::to_string(n));
      }
      if (!live(vb)) {
        const int nx = (vb >= 0 && vb < n)
                           ? __corners_next[static_cast<size_t>(vb)]
                           : -2;
        const int pv = (vb >= 0 && vb < n)
                           ? __corners_prev[static_cast<size_t>(vb)]
                           : -2;
        const int cv = (vb >= 0 && vb < n)
                           ? __corners_vert[static_cast<size_t>(vb)]
                           : -2;
        int still = 0;
        int sample = -1;
        for (int c = 0; c < n; ++c) {
          if (!live(c) || __corners_vert[static_cast<size_t>(c)] != vi)
            continue;
          if (sample < 0)
            sample = c;
          ++still;
        }
        report_link_failure(
            where, "vbegin not a live corner v=" + std::to_string(vi) +
                       " vbegin=" + std::to_string(vb) + " next=" +
                       std::to_string(nx) + " prev=" + std::to_string(pv) +
                       " vert=" + std::to_string(cv) + " live_corners=" +
                       std::to_string(still) + " sample=" +
                       std::to_string(sample) + " n=" + std::to_string(n));
      }
    }
  }

  /// Label `cid` and make it the head: also sets `vbegin(v) = cid`.
  void set_vert(CornerId cid, VertId v) {
    __corners_vert[cid] = v;
    set_vbegin(v, cid);
  }

  /// Label `cid` and make it the head: also sets `fbegin(f) = cid`.
  void set_face(CornerId cid, FaceId f) {
    __corners_face[cid] = f;
    set_fbegin(f, cid);
  }

  /// A corner from `v0` to `v1`, or -1 when no edge joins them. O(valence).
  CornerId find_edge_from_verts(VertId v0, VertId v1) const {
    CornerId c = corner_id(-1);
    const_for_each_vertex(v0, [&](CornerId ci, const shell &M) {
      CornerId cn = M.other(ci);
      if (M.vert(cn) == v1)
        c = ci;
    });
    return c;
  }

  void link(CornerId c0, CornerId c1) {
    // std::cout << __PRETTY_FUNCTION__ << c0 << " " << c1 << std::endl;
    /*
    if (vert(c0) == vert(c1)) {
      print_stacktrace();
      __builtin_frame_address(1);
    }
    assert(vert(c0) != vert(c1));
    */
    set_next(c0, c1);
    set_prev(c1, c0);
  }

  void uber_assert() {
    for (int i = 0; i < corner_count(); i++) {
      CornerId c0 = corner_id(i);
      CornerId c1 = other(c0);
      if (next(c0) < 0)
        continue;

      if (vert(c0) == vert(c1)) {
        std::cout << "    " << __PRETTY_FUNCTION__ << " c: " << c0 << " "
                  << c1 << " v: " << vert(c0) << " "
                  << vert(c1) //
                  << " vs: " << vsize(vert(c0)) << " " << vsize(vert(c1))
                  << std::endl;
      }
      if (vsize(vert(c0)) > 1)
        assert(vert(c0) != vert(c1));
    }

    for (int i = 0; i < vert_count(); i++) {
      VertId vi = vert_id(i);
      if (vbegin(vi) < 0)
        continue;

      if (next(vbegin(vi)) < 0) {
        std::cout << "-" << __PRETTY_FUNCTION__ << " v: " << i << " "
                  << vbegin(vi) << " " << std::endl
                  << std::flush;
        vprintv(vi);
      }
      assert(next(vbegin(vi)) >= 0);

      if (vert(vbegin(vi)) != vi) {
        std::cout << "-" << __PRETTY_FUNCTION__ << " v: " //
                  << i << " "                             //
                  << vbegin(vi) << " "               //
                  << vert(vbegin(vi)) << " "        //
                  << std::endl
                  << std::flush;

        for_each_vertex(
            vi, [](CornerId cid, shell &m) {
              std::cout << m.vert(cid) << " ";
            });
        std::cout << std::endl;
      }
      assert(vert(vbegin(vi)) == vi);
    }

    for (int i = 0; i < face_count(); i++) {
      FaceId fi = face_id(i);
      if (__face_begin[i] < 0)
        continue;
      if (fsize(fi) != 3) {
        std::cout << "-" << __PRETTY_FUNCTION__ << " f: " << i << " "
                  << __face_begin[i] << " " << fsize(fi) << std::endl;
        fprintv(fi);
      }
      assert(fsize(fi) < 4);
    }
  }

  void inflate_verts(size_t N) {
    __vert_begin.resize(__vert_begin.size() + N, -1);
  }

  void inflate_faces(size_t N) {
    __face_begin.resize(__face_begin.size() + N, -1);
  }

  void inflate_edge_pairs(size_t N) {
    size_t Ns = __corners_next.size() + 2 * N;
    __corners_next.resize(Ns, -1);
    __corners_prev.resize(Ns, -1);
    __corners_vert.resize(Ns, -1);
    __corners_face.resize(Ns, -1);
  }

  VertId insert_vertex() {
    __vert_begin.push_back(-1);
    return vert_id(static_cast<int>(__vert_begin.size() - 1));
  }

  FaceId insert_face() {
    __face_begin.push_back(-1);
    return face_id(static_cast<int>(__face_begin.size() - 1));
  }

  CornerId insert_edge_pair() {
    __corners_next.push_back(-1);
    __corners_next.push_back(-1);
    __corners_prev.push_back(-1);
    __corners_prev.push_back(-1);
    __corners_vert.push_back(-1);
    __corners_face.push_back(-1);
    __corners_vert.push_back(-1);
    __corners_face.push_back(-1);

    return corner_id(static_cast<int>(__corners_next.size() - 2));
  }

  void swap_rows(CornerId cA, CornerId cB) {

    int cAn = __corners_next[cA];
    int cAp = __corners_prev[cA];
    int vA = __corners_vert[cA];
    int fA = __corners_face[cA];

    int cBn = __corners_next[cB];
    int cBp = __corners_prev[cB];
    int vB = __corners_vert[cB];
    int fB = __corners_face[cB];

    set_next(prev(cA), cB);
    set_prev(next(cA), cB);
    set_vbegin(vert_id(vA), cB);
    set_fbegin(face_id(fA), cB);

    set_next(prev(cB), cA);
    set_prev(next(cB), cA);
    set_vbegin(vert_id(vB), cA);
    set_fbegin(face_id(fB), cA);

    __corners_next[cA] = cBn;
    require_corner_ref("swap_rows", static_cast<index_t>(cA), cBp);
    __corners_prev[cA] = cBp;
    __corners_vert[cA] = vB;
    __corners_face[cA] = fB;

    __corners_next[cB] = cAn;
    require_corner_ref("swap_rows", static_cast<index_t>(cB), cAp);
    __corners_prev[cB] = cAp;
    __corners_vert[cB] = vA;
    __corners_face[cB] = fA;
  }

  void flip_edge(CornerId c0) {
    CornerId c1 = other(c0);
    swap_rows(c0, c1);
  }

  void remove_vertex(VertId i) { __vert_begin[i] = -1; }

  void remove_face(FaceId i) { __face_begin[i] = -1; }

  void remove_edge_pair(CornerId i0) {
    CornerId i1 = other(i0);
    __corners_next[i0] = -1;
    __corners_next[i1] = -1;

    __corners_prev[i0] = -1;
    __corners_prev[i1] = -1;

    __corners_vert[i0] = -1;
    __corners_vert[i1] = -1;

    __corners_face[i0] = -1;
    __corners_face[i1] = -1;
  }

  /// The corners of face `i` in `next` order from `fbegin(i)` to `fend(i)`.
  /// Stops after 100 corners; `face_len` / `walk_face_ring` walk long or
  /// broken cycles without that cap.
  void for_each_face(FaceId i,
                     std::function<void(CornerId cid, shell &m)> func) {
    int j0 = this->fbegin(i);
    int j_end = this->fend(i);
    bool it = true;
    int k = 0;
    while (it && k++ < 100) {
      it = j0 != j_end;
      func(corner_id(j0), *this);
      j0 = this->next(corner_id(j0));
    }
  }

  void const_for_each_face(
      FaceId i, std::function<void(CornerId cid, const shell &m)> func) const {
    int j0 = this->fbegin(i);
    int j_end = this->fend(i);
    bool it = true;
    int k = 0;
    while (it && k++ < 100) {
      it = j0 != j_end;
      func(corner_id(j0), *this);
      j0 = this->next(corner_id(j0));
    }
  }

  void for_each_face_tri(
      FaceId i,
      std::function<void(CornerId c0, CornerId c1, CornerId c2, shell &m)>
          func) {
    int j1 = this->fbegin(i);
    int j2 = this->next(corner_id(j1));
    int j0 = this->fend(i);
    bool it = true;
    int k = 0;
    while (it && k++ < 100) {
      it = j2 != j0;
      func(corner_id(j0), corner_id(j1), corner_id(j2), *this);
      int jn = this->next(corner_id(j2));
      j1 = j2;
      j2 = jn;
    }
  }

  void const_for_each_face_tri(
      FaceId i,
      std::function<void(CornerId c0, CornerId c1, CornerId c2,
                         const shell &m)> func) const {
    int j1 = this->fbegin(i);
    int j2 = this->next(corner_id(j1));
    int j0 = this->fend(i);
    bool it = true;
    int k = 0;
    while (it && k++ < 100) {
      it = j2 != j0;
      func(corner_id(j0), corner_id(j1), corner_id(j2), *this);
      int jn = this->next(corner_id(j2));
      j1 = j2;
      j2 = jn;
    }
  }

  /// The outgoing corners of vertex `i` in `vnext` order from `vbegin(i)` to
  /// `vend(i)`. On corner `c` the ring neighbor is `vert(other(c))` and the
  /// edge opposite `i` is `next(c)`. Aborts when the ring does not close.
  void for_each_vertex(VertId i,
                       std::function<void(CornerId cid, shell &m)> func) {
    // Walk the whole ring. The old bound of 100 left corners of a merged
    // vertex labeled with the vertex merge_edge had already removed.
    const int j_start = this->vbegin(i);
    int j0 = j_start;
    const int n = static_cast<int>(corner_count());
    for (int k = 0; k < n; ++k) {
      func(corner_id(j0), *this);
      j0 = this->vnext(corner_id(j0));
      if (j0 == j_start)
        return;
    }
    report_link_failure("for_each_vertex",
                        "ring did not close v=" + std::to_string(static_cast<int>(i)) +
                            " vbegin=" + std::to_string(j_start));
  }

  void const_for_each_vertex(
      VertId i, std::function<void(CornerId cid, const shell &m)> func) const {
    const int j_start = this->vbegin(i);
    int j0 = j_start;
    const int n = static_cast<int>(corner_count());
    for (int k = 0; k < n; ++k) {
      func(corner_id(j0), *this);
      j0 = this->vnext(corner_id(j0));
      if (j0 == j_start)
        return;
    }
    report_link_failure("for_each_vertex",
                        "ring did not close v=" + std::to_string(static_cast<int>(i)) +
                            " vbegin=" + std::to_string(j_start));
  }

  std::vector<CornerId> get_edge_range() const {
    std::vector<CornerId> range;
    range.reserve(corner_count() / 2);
    for (int i = 0; i < corner_count(); i += 2) {
      if (__corners_next[i] > -1)
        range.push_back(corner_id(i));
    }
    return range;
  }

  std::vector<index_t> get_edge_range_2() const {
    std::vector<index_t> range;
    range.reserve(corner_count() / 2);
    // replace this with some c++isms
    index_t ii = 0;
    for (int i = 0; i < corner_count(); i += 2) {
      if (__corners_next[i] > -1)
        range.push_back(ii++);
    }
    return range;
  }

  std::vector<VertId> get_vert_range() const {
    std::vector<VertId> range;
    range.reserve(vert_count());
    for (int i = 0; i < vert_count(); i++) {
      if (__vert_begin[i] > -1)
        range.push_back(vert_id(i));
    }
    return range;
  }

  std::vector<FaceId> get_face_range(bool filter_tris = true) const {
    std::vector<FaceId> range;
    range.reserve(face_count());
    for (int i = 0; i < face_count(); i++) {
      FaceId fi = face_id(i);
      if (fbegin(fi) < 0)
        continue;
      if (filter_tris && fsize(fi) != 3)
        continue;
      range.push_back(fi);
    }
    return range;
  }

  std::vector<index_t> get_edge_vert_ids() {
    std::vector<index_t> range;
    range.reserve(corner_count());
    for (int i = 0; i < corner_count(); i += 2) {
      if (__corners_next[i] < 0)
        continue;
      CornerId c = corner_id(i);
      range.push_back(vert(c));
      range.push_back(vert(other(c)));
    }
    return range;
  }

  std::vector<index_t> get_face_vert_ids(bool filter_tris = true) {
    std::vector<int> faces;
    faces.reserve(3 * face_count());

    for (int i = 0; i < face_count(); i++) {
      FaceId fi = face_id(i);
      if (fbegin(fi) < 0)
        continue;

      if (filter_tris && fsize(fi) != 3)
        continue;

      for_each_face(fi, [&faces](CornerId ci, shell &M) {
        faces.push_back(M.vert(ci));
      });
    }
    return faces;
  }

  /// Ring neighbors of `iv` in `for_each_vertex` order.
  std::vector<VertId> get_one_ring(VertId iv) const {
    std::vector<VertId> range;
    this->const_for_each_vertex(iv, [&range](CornerId ci, const shell &M) {
      range.push_back(M.vert(M.next(ci)));
    });
    return range;
  }

  /// The link of `iv`: corner `next(c)` for each spoke `c`, in the same
  /// order. On a triangle fan these are the ring edges `w_k -> w_{k+1}`
  /// opposite `iv`, and they chain head to tail.
  std::vector<CornerId> link_corners(VertId iv) const {
    std::vector<CornerId> range;
    this->const_for_each_vertex(iv, [&range](CornerId ci, const shell &M) {
      range.push_back(M.next(ci));
    });
    return range;
  }

  std::array<VertId, 3> get_tri(FaceId i) const {
    std::array<VertId, 3> tri = {vert_id(-1), vert_id(-1), vert_id(-1)};
    int j = 0;
    this->const_for_each_face(
        i, [&tri, &j](CornerId ci, const shell &M) { tri[j++] = M.vert(ci); });
    return tri;
  }

  std::vector<index_t> get_range_map(const std::vector<index_t> indices,
                                     int stride, bool filter_tris = true) {
    std::vector<index_t> map;
    map.reserve(indices.size() / stride);
    for (int i = 0; i < static_cast<int>(indices.size()); i += stride) {
      if (indices[i] < 0)
        continue;

      if (filter_tris && fsize(face_id(i)) != 3)
        continue;

      map.push_back(i);
    }
    return map;
  }

  std::vector<index_t> get_edge_map() {
    return get_range_map(__corners_next, 2, false);
  }
  std::vector<index_t> get_vert_map() {
    return get_range_map(__vert_begin, 1, false);
  }
  std::vector<index_t> get_face_map(bool filter_tris = true) {
    return get_range_map(__face_begin, 1, filter_tris);
  }

  void fupdate(FaceId f) {
    for_each_face(f, [f](CornerId cid, shell &m) { m.set_face(cid, f); });
  }

  void vupdate(VertId v) {
    for_each_vertex(v, [v](CornerId cid, shell &m) { m.set_vert(cid, v); });
  }

  size_t fsize(FaceId f) const {
    size_t s = 0;
    const_for_each_face(f, [&s](CornerId cid, const shell &m) { s++; });
    return s;
  }

  size_t vsize(VertId v) const {
    size_t s = 0;
    const_for_each_vertex(v, [&s](CornerId cid, const shell &m) { s++; });
    return s;
  }

  void cprint(CornerId c0) {
    CornerId c1 = other(c0);
    CornerId c0p = prev(c0);
    CornerId c1p = prev(c1);

    VertId v0 = vert(c0);
    VertId v1 = vert(c1);
    VertId v0p = vert(c0p);
    VertId v1p = vert(c1p);
    std::cout << "// =========== " << std::endl;
    std::cout << "//    " << c0p << std::endl;
    std::cout << "//" << c0 << "---" << c1 << std::endl;
    std::cout << "//    " << c1p << std::endl;
    std::cout << "//    " << v0p << std::endl;
    std::cout << "//" << v0 << "---" << v1 << std::endl;
    std::cout << "//    " << v1p << std::endl;
    std::cout << "// vs: " << vsize(v0) << " " << vsize(v1) << std::endl;
    std::cout << "// vps: " << vsize(v0) << " " << vsize(v1) << std::endl;
    vprint_graph_viz(v0);
    vprint_graph_viz(v1);
    vprint_graph_viz(v0p);
    vprint_graph_viz(v1p);
    vprintvs(v0);
    vprintvs(v1);
    vprintvs(v0p);
    vprintvs(v1p);

    std::cout << "// =========== " << std::endl;
  }

  void fprint(FaceId f) {
    std::cout << "//f-" << f << ": ";
    for_each_face(f, [](CornerId cid, shell &m) {
      std::cout << cid << " ";
    });
    std::cout << std::endl;
  }

  void fprintv(FaceId f) {
    std::cout << "//fvs-" << f << ": " << std::flush;
    for_each_face(f, [this](CornerId cid, shell &m) {
      std::cout << this->vert(cid) << " " << std::flush;
    });
    std::cout << std::endl << std::flush;
  }

  void vprint(VertId v) {
    std::cout << "//v-" << v << ": " << std::flush;
    for_each_vertex(v, [](CornerId cid, shell &m) {
      std::cout << cid << " " << std::flush;
    });
    std::cout << std::endl << std::flush;
  }

  void vprintv(VertId v) {
    std::cout << "//vvs-" << v << ": " << std::flush;
    for_each_vertex(v, [this](CornerId cid, shell &m) {
      std::cout << this->vert(this->next(cid)) << " " << std::flush;
    });
    std::cout << std::endl << std::flush;
  }

  void vprintvs(VertId v) {
    std::cout << "//vvs-" << v << ": " << std::flush;
    for_each_vertex(v, [this](CornerId cid, shell &m) {
      std::cout << this->vsize(this->vert(this->next(cid))) << " "
                << std::flush;
    });
    std::cout << std::endl << std::flush;
  }

  void vprintfs(VertId v) {
    std::cout << "//vfs-" << v << ": ";
    for_each_vertex(v, [this](CornerId cid, shell &m) {
      std::cout << this->fsize(this->face(cid)) << " ";
    });
    std::cout << std::endl;
  }

  void vprint_graph_viz(VertId v) {
    std::cout << "//    ========" << std::endl;
    std::cout << "//    vg-" << v << std::endl;
    for_each_vertex(v, [this](CornerId cid, shell &m) {
      std::cout << "    " << this->vert(cid) << " -> "
                << this->vert(this->next(cid)) << std::endl;
    });
    std::cout << "//    ========" << std::endl;
  }

  // accessors
  std::vector<index_t> &corners_next() { return __corners_next; }
  std::vector<index_t> &corners_prev() { return __corners_prev; }
  std::vector<index_t> &corners_vert() { return __corners_vert; }
  std::vector<index_t> &corners_face() { return __corners_face; }
  std::vector<index_t> &vert_begin() { return __vert_begin; }
  std::vector<index_t> &face_begin() { return __face_begin; }

  std::vector<index_t> __corners_next;
  std::vector<index_t> __corners_prev;
  std::vector<index_t> __corners_vert;
  std::vector<index_t> __corners_face;

  std::vector<index_t> __vert_begin;
  std::vector<index_t> __face_begin;
  GeometryChannels __channels;
};
} // namespace shell
} // namespace asawa
} // namespace gaudi
#endif