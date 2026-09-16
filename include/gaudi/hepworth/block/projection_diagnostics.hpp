#ifndef GAUDI_HEPWORTH_PROJECTION_DIAGNOSTICS_HPP
#define GAUDI_HEPWORTH_PROJECTION_DIAGNOSTICS_HPP

#include "gaudi/common.h"
#include "gaudi/hepworth/block/slab_prism_log.hpp"
#include "gaudi/hepworth/projection_constraint.hpp"

#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

namespace gaudi {
namespace hepworth {
namespace block {

inline bool projection_diag_enabled() {
  const char *env = std::getenv("GAUDI_PROJECTION_DIAG");
  return env != nullptr && env[0] != '\0' && env[0] != '0';
}

inline void projection_diag_log(const std::string &msg) {
  std::cerr << "[projection_diag] " << msg << std::endl;
  std::cerr.flush();
}

class phase_timer {
public:
  explicit phase_timer(std::string label) : _label(std::move(label)) {
    _t0 = clock::now();
  }

  ~phase_timer() {
    if (!projection_diag_enabled() || _label.empty())
      return;
    const double ms =
        std::chrono::duration<double, std::milli>(clock::now() - _t0).count();
    projection_diag_log(_label + " ms=" + std::to_string(ms));
  }

private:
  using clock = std::chrono::steady_clock;
  std::string _label;
  clock::time_point _t0;
};

struct triplet_scan_report {
  size_t count = 0;
  size_t non_finite = 0;
};

inline triplet_scan_report scan_triplets(const std::vector<trip> &triplets) {
  triplet_scan_report out;
  out.count = triplets.size();
  for (const trip &t : triplets) {
    const double v = t.value();
    if (!std::isfinite(v))
      ++out.non_finite;
  }
  return out;
}

struct sparse_matrix_report {
  size_t rows = 0;
  size_t cols = 0;
  size_t nnz = 0;
  size_t non_finite_nnz = 0;
  size_t tiny_diag_count = 0;
  double min_diag = std::numeric_limits<double>::infinity();
  double max_diag = 0.0;
};

inline sparse_matrix_report diagnose_sparse(const matS &M,
                                            double tiny_thresh = 1e-9) {
  sparse_matrix_report out;
  out.rows = static_cast<size_t>(M.rows());
  out.cols = static_cast<size_t>(M.cols());
  out.nnz = static_cast<size_t>(M.nonZeros());
  for (int k = 0; k < M.outerSize(); ++k) {
    for (matS::InnerIterator it(M, k); it; ++it) {
      const double v = it.value();
      if (!std::isfinite(v))
        ++out.non_finite_nnz;
      if (it.row() == it.col()) {
        out.min_diag = std::min(out.min_diag, v);
        out.max_diag = std::max(out.max_diag, v);
        if (std::abs(v) < tiny_thresh)
          ++out.tiny_diag_count;
      }
    }
  }
  if (out.min_diag == std::numeric_limits<double>::infinity())
    out.min_diag = 0.0;
  return out;
}

inline void log_matrix_report(const std::string &label,
                              const sparse_matrix_report &r) {
  const double diag_ratio =
      r.min_diag > 0.0 ? r.max_diag / r.min_diag : 0.0;
  projection_diag_log(label + " rows=" + std::to_string(r.rows) + " cols=" +
                      std::to_string(r.cols) + " nnz=" + std::to_string(r.nnz) +
                      " non_finite_nnz=" + std::to_string(r.non_finite_nnz) +
                      " min_diag=" + std::to_string(r.min_diag) +
                      " max_diag=" + std::to_string(r.max_diag) +
                      " diag_ratio=" + std::to_string(diag_ratio) +
                      " tiny_diag=" + std::to_string(r.tiny_diag_count));
}

inline bool build_constraint_matrices(
    const std::vector<projection_constraint::ptr> &constraints, int nm,
    matS &A, matS &AtA, index_t &a_rows) {
  std::vector<trip> triplets;
  index_t id0 = 0;
  for (const auto &constraint : constraints) {
    if (constraint)
      constraint->fill_A(id0, triplets);
  }
  a_rows = id0;
  A = matS(a_rows, nm);
  A.setFromTriplets(triplets.begin(), triplets.end());
  AtA = A.transpose() * A;
  return a_rows > 0 && nm > 0;
}

inline void log_constraint_counts(
    const std::vector<projection_constraint::ptr> &constraints) {
  std::unordered_map<std::string, size_t> counts;
  for (const auto &c : constraints) {
    if (c)
      ++counts[c->name()];
  }
  for (const auto &[name, n] : counts)
    projection_diag_log("constraints " + name + "=" + std::to_string(n));
}

/// Outer ADMM index for projection diagnostics (constraints read this).
inline int &projection_diag_outer_k() {
  static int k = -1;
  return k;
}

inline bool projection_diag_sample_pass() {
  return projection_diag_enabled() && projection_diag_outer_k() == 0;
}

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // GAUDI_HEPWORTH_PROJECTION_DIAGNOSTICS_HPP
