#ifndef GAUDI_HEPWORTH_SLAB_PRISM_LOG_HPP
#define GAUDI_HEPWORTH_SLAB_PRISM_LOG_HPP

#include <atomic>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <string>

namespace gaudi {
namespace hepworth {
namespace block {

/// Coarse slab-prism progress. Set GAUDI_SLAB_PRISM_LOG=1 to enable.
inline bool slab_prism_log_enabled() {
  const char *env = std::getenv("GAUDI_SLAB_PRISM_LOG");
  if (env == nullptr)
    return false;
  return env[0] != '\0' && env[0] != '0';
}

inline void slab_prism_log(const std::string &msg) {
  if (!slab_prism_log_enabled())
    return;
  std::cerr << "[slab_prism] " << msg << std::endl;
  std::cerr.flush();
}

/// Unconditional stderr trace (flushed). Use sparingly for first-solve milestones.
inline void slab_prism_trace(const std::string &msg) {
  std::cerr << "[slab_prism] " << msg << std::endl;
  std::cerr.flush();
}

inline bool slab_prism_trace_claim_first_solve() {
  static std::atomic<bool> claimed{false};
  return !claimed.exchange(true);
}

class slab_prism_phase_timer {
public:
  explicit slab_prism_phase_timer(std::string label) : _label(std::move(label)) {
    if (slab_prism_log_enabled())
      _t0 = clock::now();
  }

  ~slab_prism_phase_timer() {
    if (!slab_prism_log_enabled() || _label.empty())
      return;
    const double ms =
        std::chrono::duration<double, std::milli>(clock::now() - _t0).count();
    slab_prism_log(_label + " ms=" + std::to_string(ms));
  }

private:
  using clock = std::chrono::steady_clock;
  std::string _label;
  clock::time_point _t0{};
};

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // GAUDI_HEPWORTH_SLAB_PRISM_LOG_HPP
