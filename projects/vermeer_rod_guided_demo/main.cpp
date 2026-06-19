#include "gaudi/duchamp/rod_guided_demo_adapter.hpp"
#include "gaudi/vermeer/duchamp_host.hpp"

int main() {
  auto demo = gaudi::duchamp::rod_guided_demo_adapter::create();
  gaudi::vermeer::duchamp_host host(std::move(demo));
  return host.run();
}
