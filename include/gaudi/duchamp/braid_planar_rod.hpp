#ifndef GAUDI_DUCHAMP_BRAID_PLANAR_ROD_HPP
#define GAUDI_DUCHAMP_BRAID_PLANAR_ROD_HPP

#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/common.h"
#include "gaudi/windychien/braid.hpp"
#include "gaudi/windychien/catalog.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace gaudi {
namespace duchamp {

struct braid_planar_params {
  real dx = 1.0;
  real dy = 1.0;
  real eps_z = 0.05; // over/under bead height at crossing midpoints
  bool center = true;
  /// If true: L columns; last generator links to column 0 (closed braid).
  /// Centering is skipped (breaks periodic x).
  bool close = false;
};

namespace detail {

inline asawa::rod::CornerId insert_vert(asawa::rod::rod &R, const vec3 &p) {
  const asawa::rod::CornerId id = R.insert_edge();
  R.x()[id] = p;
  return id;
}

inline std::vector<std::vector<asawa::rod::CornerId>>
collect_strand_paths(const asawa::rod::rod &R, bool closed) {
  std::vector<std::vector<asawa::rod::CornerId>> paths;
  std::vector<char> visited(R.corner_count(), 0);

  auto walk_from = [&](asawa::rod::CornerId start) {
    std::vector<asawa::rod::CornerId> path;
    asawa::rod::CornerId i = start;
    do {
      const int ii = static_cast<int>(i);
      if (visited[static_cast<size_t>(ii)]) {
        break;
      }
      visited[static_cast<size_t>(ii)] = 1;
      path.push_back(i);
      i = R.next(i);
    } while (i >= asawa::rod::corner_id(0) && i != start);
    if (!path.empty()) {
      paths.push_back(std::move(path));
    }
  };

  if (!closed) {
    for (size_t i = 0; i < R.corner_count(); ++i) {
      const auto ci = asawa::rod::corner_id(static_cast<int>(i));
      if (R.prev(ci) < asawa::rod::corner_id(0) &&
          R.next(ci) >= asawa::rod::corner_id(0)) {
        walk_from(ci);
      }
    }
  } else {
    for (size_t i = 0; i < R.corner_count(); ++i) {
      if (visited[i]) {
        continue;
      }
      const auto ci = asawa::rod::corner_id(static_cast<int>(i));
      if (R.next(ci) < asawa::rod::corner_id(0)) {
        continue;
      }
      walk_from(ci);
    }
  }
  return paths;
}

/// If z[i]==0 and z[i-1]==z[i+1] (nonzero), set z[i] to that elevation.
/// Index arithmetic uses the path buffer with %N when closed.
inline void elevate_flat_verts_between_equal_neighbors(
    asawa::rod::rod &R,
    const std::vector<std::vector<asawa::rod::CornerId>> &strand_paths,
    bool closed, real z_tol = 1e-14) {
  for (const auto &path : strand_paths) {
    const int N = static_cast<int>(path.size());
    if (N < 3) {
      continue;
    }
    bool changed = true;
    while (changed) {
      changed = false;
      for (int i = 0; i < N; ++i) {
        if (!closed && (i == 0 || i + 1 == N)) {
          continue;
        }
        const int im = (i - 1 + N) % N;
        const int ip = (i + 1) % N;
        if (!closed && (im < 0 || ip >= N)) {
          continue;
        }
        vec3 &qi =
            R.corner_verts()[static_cast<size_t>(path[static_cast<size_t>(i)])];
        if (std::abs(qi[2]) > z_tol) {
          continue;
        }
        const real zm =
            R.corner_verts()[static_cast<size_t>(path[static_cast<size_t>(im)])]
                [2];
        const real zp =
            R.corner_verts()[static_cast<size_t>(path[static_cast<size_t>(ip)])]
                [2];
        if (std::abs(zm - zp) > z_tol || std::abs(zm) <= z_tol) {
          continue;
        }
        qi[2] = zm;
        changed = true;
      }
    }
  }
}

/// One strand crossing several rows inside a single column.
/// `lo`/`hi` are inclusive 0-based rows. `up` sends the strand at `lo` to `hi`.
/// `gens` is the Artin word in travel order. A sign flip along the run is a
/// weave: the moving strand switches over/under at that row.
struct column_slide {
  int lo = 0;
  int hi = 0;
  bool up = true;
  std::vector<int> gens;
};

/// Ascending run: lower strand travels up. Descending run: upper strand
/// travels down. Signs may change along the run. A lone generator is the
/// adjacent swap (positive: lower strand over).
inline std::vector<column_slide>
parse_column_slides(const std::vector<int> &gens, int n) {
  std::vector<column_slide> slides;
  size_t i = 0;
  while (i < gens.size()) {
    const int a0 = std::abs(gens[i]);
    int step = 0;
    if (i + 1 < gens.size()) {
      const int a1 = std::abs(gens[i + 1]);
      if (a1 == a0 + 1 || a1 == a0 - 1) {
        step = a1 - a0;
      }
    }
    size_t k = i + 1;
    if (step != 0) {
      while (k < gens.size() &&
             std::abs(gens[k]) == std::abs(gens[k - 1]) + step) {
        ++k;
      }
    }
    const int a_last = std::abs(gens[k - 1]);
    column_slide s;
    s.gens.assign(gens.begin() + static_cast<std::ptrdiff_t>(i),
                  gens.begin() + static_cast<std::ptrdiff_t>(k));
    if (step >= 0) {
      s.lo = a0 - 1;
      s.hi = a_last;
      s.up = true;
    } else {
      s.hi = a0;
      s.lo = a_last - 1;
      s.up = false;
    }
    if (s.lo < 0 || s.hi >= n || s.hi <= s.lo ||
        static_cast<int>(s.gens.size()) != s.hi - s.lo) {
      throw std::runtime_error("braid_to_planar_rod: slide leaves the braid");
    }
    slides.push_back(std::move(s));
    i = k;
  }

  std::vector<char> occ(static_cast<size_t>(n), 0);
  for (const column_slide &s : slides) {
    for (int r = s.lo; r <= s.hi; ++r) {
      if (occ[static_cast<size_t>(r)]) {
        throw std::runtime_error(
            "braid_to_planar_rod: overlapping generators in one layer");
      }
      occ[static_cast<size_t>(r)] = 1;
    }
  }
  return slides;
}

} // namespace detail

/// Planar braid: coplanar column verts (z=0) + ±eps crossing beads,
/// then elevate flat verts sandwiched between equal elevations.
///
/// One axial column per layer (braid.layer_counts); empty layer_counts ⇒
/// one column per word entry (legacy Artin word).
/// A column may be one adjacent swap, or a run of consecutive generators
/// that slides one strand across several rows. The moving strand is a
/// straight diagonal; each crossed strand gets an over/under bead where
/// the diagonal meets it. A sign change along the run weaves that row the
/// other way. Positive ascending run: the lower strand travels up and over.
inline asawa::rod::rod::ptr
braid_to_planar_rod(const windychien::braid &b,
                    const braid_planar_params &p = {}) {
  windychien::validate_braid(b);
  const int n = b.strands;
  const int L = windychien::braid_n_columns(b);
  if (p.close && L < 1) {
    throw std::runtime_error(
        "braid_to_planar_rod: close requires nonempty word");
  }

  std::vector<std::vector<int>> layers;
  layers.reserve(static_cast<size_t>(L));
  if (b.layer_counts.empty()) {
    for (int g : b.word) {
      layers.push_back({g});
    }
  } else {
    size_t off = 0;
    for (int cnt : b.layer_counts) {
      std::vector<int> layer(
          b.word.begin() + static_cast<std::ptrdiff_t>(off),
          b.word.begin() + static_cast<std::ptrdiff_t>(off + cnt));
      detail::parse_column_slides(layer, n);
      layers.push_back(std::move(layer));
      off += static_cast<size_t>(cnt);
    }
  }

  asawa::rod::rod::ptr R = asawa::rod::rod::create();

  std::vector<asawa::rod::CornerId> curr(static_cast<size_t>(n),
                                         asawa::rod::corner_id(-1));
  std::vector<asawa::rod::CornerId> col0(static_cast<size_t>(n),
                                         asawa::rod::corner_id(-1));

  for (int r = 0; r < n; ++r) {
    const asawa::rod::CornerId id =
        detail::insert_vert(*R, vec3(0.0, p.dy * real(r), 0.0));
    curr[static_cast<size_t>(r)] = id;
    col0[static_cast<size_t>(r)] = id;
  }

  for (int c = 0; c < L; ++c) {
    const bool wrap = p.close && (c + 1 == L);
    const real x0 = p.dx * real(c);
    // Wrap destination is x=0 ≡ L*dx; beads sit on that seam interval.
    const real x1 = wrap ? p.dx * real(L) : p.dx * real(c + 1);

    std::vector<asawa::rod::CornerId> next_col(static_cast<size_t>(n),
                                               asawa::rod::corner_id(-1));
    if (wrap) {
      for (int r = 0; r < n; ++r) {
        next_col[static_cast<size_t>(r)] = col0[static_cast<size_t>(r)];
      }
    } else {
      for (int r = 0; r < n; ++r) {
        next_col[static_cast<size_t>(r)] =
            detail::insert_vert(*R, vec3(x1, p.dy * real(r), 0.0));
      }
    }

    std::vector<char> crossed(static_cast<size_t>(n), 0);
    const std::vector<detail::column_slide> slides =
        detail::parse_column_slides(layers[static_cast<size_t>(c)], n);
    for (const detail::column_slide &s : slides) {
      const int jsrc = s.up ? s.lo : s.hi;
      const int jdst = s.up ? s.hi : s.lo;
      const int span = s.hi - s.lo;
      std::vector<int> dest(static_cast<size_t>(n), -1);
      if (s.up) {
        dest[static_cast<size_t>(s.lo)] = s.hi;
        for (int r = s.lo + 1; r <= s.hi; ++r) {
          dest[static_cast<size_t>(r)] = r - 1;
        }
      } else {
        dest[static_cast<size_t>(s.hi)] = s.lo;
        for (int r = s.lo; r < s.hi; ++r) {
          dest[static_cast<size_t>(r)] = r + 1;
        }
      }

      struct hit {
        real t;
        int row;
      };
      std::vector<hit> hits;
      hits.reserve(static_cast<size_t>(span));
      for (int r = s.lo; r <= s.hi; ++r) {
        if (r == jsrc) {
          continue;
        }
        const real t = s.up ? real(r - s.lo) / real(span + 1)
                            : real(s.hi - r) / real(span + 1);
        hits.push_back(hit{t, r});
      }
      std::sort(hits.begin(), hits.end(),
                [](const hit &a, const hit &b) { return a.t < b.t; });

      const real ysrc = p.dy * real(jsrc);
      const real ydst = p.dy * real(jdst);
      asawa::rod::CornerId prev = curr[static_cast<size_t>(jsrc)];
      int hit_i = 0;
      for (const hit &h : hits) {
        const int g = s.gens[static_cast<size_t>(hit_i)];
        ++hit_i;
        const bool jumper_over = s.up == (g > 0);
        const real x = (real(1) - h.t) * x0 + h.t * x1;
        const real y = (real(1) - h.t) * ysrc + h.t * ydst;
        const asawa::rod::CornerId over =
            detail::insert_vert(*R, vec3(x, y, p.eps_z));
        const asawa::rod::CornerId under =
            detail::insert_vert(*R, vec3(x, y, -p.eps_z));
        const asawa::rod::CornerId jbead = jumper_over ? over : under;
        const asawa::rod::CornerId obead = jumper_over ? under : over;
        R->link(prev, jbead);
        prev = jbead;
        R->link(curr[static_cast<size_t>(h.row)], obead);
        R->link(obead, next_col[static_cast<size_t>(dest[static_cast<size_t>(h.row)])]);
        crossed[static_cast<size_t>(h.row)] = 1;
      }
      R->link(prev, next_col[static_cast<size_t>(jdst)]);
      crossed[static_cast<size_t>(jsrc)] = 1;
    }

    for (int r = 0; r < n; ++r) {
      if (crossed[static_cast<size_t>(r)]) {
        continue;
      }
      R->link(curr[static_cast<size_t>(r)], next_col[static_cast<size_t>(r)]);
    }

    curr = next_col;
  }

  const auto paths = detail::collect_strand_paths(*R, p.close);
  detail::elevate_flat_verts_between_equal_neighbors(*R, paths, p.close);

  R->_init_params();

  if (p.center && !p.close && !R->x().empty()) {
    vec3 cen = vec3::Zero();
    for (const vec3 &q : R->x()) {
      cen += q;
    }
    cen /= real(R->x().size());
    for (vec3 &q : R->x()) {
      q -= cen;
    }
  }

  return R;
}

inline asawa::rod::rod::ptr
load_braid_planar_rod(const std::string &knot_name,
                      const braid_planar_params &p = {}) {
  return braid_to_planar_rod(windychien::load_braid(knot_name), p);
}

} // namespace duchamp
} // namespace gaudi

#endif
