#ifndef GAUDI_DUCHAMP_BRAID_PLANAR_ROD_HPP
#define GAUDI_DUCHAMP_BRAID_PLANAR_ROD_HPP

#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/common.h"
#include "gaudi/windychien/braid.hpp"
#include "gaudi/windychien/catalog.hpp"

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
  R.corner_verts().push_back(p);
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

} // namespace detail

/// Planar braid: coplanar column verts (z=0) + ±eps crossing midpoints,
/// then elevate flat verts sandwiched between equal elevations.
///
/// One axial column per layer (braid.layer_counts); empty layer_counts ⇒
/// one column per word entry (legacy Artin word).
/// Positive generator: lower track through +eps bead.
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
      std::vector<char> used(static_cast<size_t>(n), 0);
      for (int g : layer) {
        const int a = std::abs(g) - 1;
        const int bb = a + 1;
        if (used[static_cast<size_t>(a)] || used[static_cast<size_t>(bb)]) {
          throw std::runtime_error(
              "braid_to_planar_rod: overlapping generators in one layer");
        }
        used[static_cast<size_t>(a)] = 1;
        used[static_cast<size_t>(bb)] = 1;
      }
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
    // Wrap destination is x=0 ≡ L*dx; midpoint must sit on the periodic seam.
    const real x1 = wrap ? p.dx * real(L) : p.dx * real(c + 1);
    const real xm = 0.5 * (x0 + x1);

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
    for (int g : layers[static_cast<size_t>(c)]) {
      const int track_a = std::abs(g) - 1;
      const int track_b = track_a + 1;
      const real ya = p.dy * real(track_a);
      const real yb = p.dy * real(track_b);
      const real ym = 0.5 * (ya + yb);

      const asawa::rod::CornerId over =
          detail::insert_vert(*R, vec3(xm, ym, p.eps_z));
      const asawa::rod::CornerId under =
          detail::insert_vert(*R, vec3(xm, ym, -p.eps_z));
      const asawa::rod::CornerId via_a = (g > 0) ? over : under;
      const asawa::rod::CornerId via_b = (g > 0) ? under : over;

      R->link(curr[static_cast<size_t>(track_a)], via_a);
      R->link(via_a, next_col[static_cast<size_t>(track_b)]);
      R->link(curr[static_cast<size_t>(track_b)], via_b);
      R->link(via_b, next_col[static_cast<size_t>(track_a)]);
      crossed[static_cast<size_t>(track_a)] = 1;
      crossed[static_cast<size_t>(track_b)] = 1;
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
