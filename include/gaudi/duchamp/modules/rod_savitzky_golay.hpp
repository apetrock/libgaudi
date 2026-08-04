#ifndef GAUDI_DUCHAMP_MODULES_ROD_SAVITZKY_GOLAY_HPP
#define GAUDI_DUCHAMP_MODULES_ROD_SAVITZKY_GOLAY_HPP

#include <algorithm>
#include <cmath>
#include <vector>

#include <Eigen/Dense>

#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/common.h"

namespace gaudi {
namespace duchamp {

struct rod_sg_filter_config {
  bool enable = false;
  /// Odd window length ≥ 3 (samples along each strand).
  int window = 5;
  /// Polynomial order; must satisfy poly_order < window.
  int poly_order = 2;
};

inline int rod_sg_clamp_window(int window) {
  window = std::max(window, 3);
  if (window % 2 == 0)
    ++window;
  return window;
}

/// Savitzky–Golay smoothing coefficients (0th derivative at window center).
inline std::vector<real> rod_sg_coefficients(int window, int poly_order) {
  window = rod_sg_clamp_window(window);
  poly_order = std::max(0, std::min(poly_order, window - 1));

  const int half = window / 2;
  Eigen::MatrixXd A(window, poly_order + 1);
  for (int j = 0; j < window; ++j) {
    const real x = real(j - half);
    real xp = 1.0;
    for (int i = 0; i <= poly_order; ++i) {
      A(j, i) = xp;
      xp *= x;
    }
  }
  const Eigen::MatrixXd pinv =
      (A.transpose() * A).ldlt().solve(A.transpose());
  std::vector<real> coeffs(static_cast<size_t>(window));
  for (int j = 0; j < window; ++j)
    coeffs[static_cast<size_t>(j)] = pinv(0, j);
  return coeffs;
}

inline void rod_strand_arcs(const asawa::rod::rod &rod, std::vector<int> &strand,
                            std::vector<int> &arc,
                            std::vector<int> &strand_len) {
  const size_t n = rod.x().size();
  strand.assign(n, -1);
  arc.assign(n, -1);
  strand_len.clear();
  std::vector<char> visited(n, 0);
  int sid = 0;
  for (size_t seed = 0; seed < n; ++seed) {
    if (visited[seed])
      continue;
    const asawa::rod::CornerId cseed =
        asawa::rod::corner_id(static_cast<int>(seed));
    if (rod.next(cseed) < asawa::rod::corner_id(0) &&
        rod.prev(cseed) < asawa::rod::corner_id(0)) {
      visited[seed] = 1;
      strand[seed] = sid;
      arc[seed] = 0;
      strand_len.push_back(1);
      ++sid;
      continue;
    }

    asawa::rod::CornerId start = cseed;
    asawa::rod::CornerId s = start;
    for (;;) {
      const asawa::rod::CornerId p = rod.prev(s);
      if (p < asawa::rod::corner_id(0) || p == start)
        break;
      s = p;
    }

    int a = 0;
    asawa::rod::CornerId i = s;
    do {
      const size_t ii = static_cast<size_t>(static_cast<int>(i));
      visited[ii] = 1;
      strand[ii] = sid;
      arc[ii] = a++;
      const asawa::rod::CornerId j = rod.next(i);
      if (j < asawa::rod::corner_id(0))
        break;
      i = j;
    } while (i != s && !visited[static_cast<size_t>(static_cast<int>(i))]);
    strand_len.push_back(a);
    ++sid;
  }
}

inline bool rod_strand_closed(const asawa::rod::rod &rod, int corner_seed) {
  const asawa::rod::CornerId c = asawa::rod::corner_id(corner_seed);
  const asawa::rod::CornerId n = rod.next(c);
  if (n < asawa::rod::corner_id(0))
    return false;
  asawa::rod::CornerId i = n;
  const asawa::rod::CornerId start = c;
  while (i != start) {
    const asawa::rod::CornerId j = rod.next(i);
    if (j < asawa::rod::corner_id(0))
      return false;
    i = j;
  }
  return true;
}

/// Filter vec3 samples along one strand (open ends clamped, closed loops wrap).
inline void rod_sg_filter_strand(const std::vector<int> &indices,
                                 const std::vector<real> &coeffs,
                                 bool closed, std::vector<vec3> &field) {
  const int window = static_cast<int>(coeffs.size());
  const int half = window / 2;
  const int L = static_cast<int>(indices.size());
  if (L < window)
    return;

  std::vector<vec3> src(L);
  for (int a = 0; a < L; ++a)
    src[static_cast<size_t>(a)] = field[static_cast<size_t>(indices[static_cast<size_t>(a)])];

  auto sample = [&](int a) -> vec3 {
    if (closed) {
      a = ((a % L) + L) % L;
      return src[static_cast<size_t>(a)];
    }
    a = std::max(0, std::min(L - 1, a));
    return src[static_cast<size_t>(a)];
  };

  for (int a = 0; a < L; ++a) {
    vec3 out = vec3::Zero();
    for (int j = 0; j < window; ++j)
      out += coeffs[static_cast<size_t>(j)] * sample(a + j - half);
    field[static_cast<size_t>(indices[static_cast<size_t>(a)])] = out;
  }
}

/// Savitzky–Golay smooth a per-corner vec3 field along each rod strand in place.
inline void rod_sg_filter_field(const asawa::rod::rod &rod,
                                std::vector<vec3> &field,
                                const rod_sg_filter_config &cfg) {
  if (!cfg.enable || field.size() != rod.x().size())
    return;

  const int window = rod_sg_clamp_window(cfg.window);
  const int poly = std::max(0, std::min(cfg.poly_order, window - 1));
  const std::vector<real> coeffs = rod_sg_coefficients(window, poly);

  std::vector<int> strand, arc, strand_len;
  rod_strand_arcs(rod, strand, arc, strand_len);

  std::vector<std::vector<int>> verts(strand_len.size());
  for (size_t i = 0; i < strand.size(); ++i) {
    const int s = strand[i];
    const int a = arc[i];
    if (s < 0 || a < 0 || static_cast<size_t>(s) >= verts.size())
      continue;
    if (static_cast<size_t>(a) >= verts[static_cast<size_t>(s)].size())
      verts[static_cast<size_t>(s)].resize(static_cast<size_t>(a + 1), -1);
    verts[static_cast<size_t>(s)][static_cast<size_t>(a)] = static_cast<int>(i);
  }

  for (size_t s = 0; s < verts.size(); ++s) {
    std::vector<int> &idx = verts[s];
    idx.erase(std::remove(idx.begin(), idx.end(), -1), idx.end());
    if (idx.empty())
      continue;
    const bool closed =
        idx.size() >= 3 && rod_strand_closed(rod, idx.front());
    rod_sg_filter_strand(idx, coeffs, closed, field);
  }
}

} // namespace duchamp
} // namespace gaudi

#endif
