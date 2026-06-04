#ifndef __GAUDI_DUCHAMP_RX_RX_PIPELINE_HPP__
#define __GAUDI_DUCHAMP_RX_RX_PIPELINE_HPP__

#include <functional>
#include <vector>

namespace gaudi {
namespace duchamp {
namespace rx {

/// Composable operator-splitting steps. Each `run()` is one PDE time step: execute
/// stages in order (reaction blocks, then diffusion, etc.); callers capture mesh and
/// fields in lambdas or functors to author different growth models.
class rx_pipeline {
public:
  void clear() { _stages.clear(); }
  void push_back(std::function<void()> f) { _stages.push_back(std::move(f)); }
  void run() {
    for (auto &s : _stages)
      s();
  }
  size_t size() const { return _stages.size(); }
  const std::vector<std::function<void()>> &stages() const { return _stages; }
  std::vector<std::function<void()>> &stages() { return _stages; }

private:
  std::vector<std::function<void()>> _stages;
};

} // namespace rx
} // namespace duchamp
} // namespace gaudi

#endif
