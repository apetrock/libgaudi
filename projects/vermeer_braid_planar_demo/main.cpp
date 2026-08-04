#include "gaudi/duchamp/braid_planar_demo.hpp"
#include "gaudi/vermeer/vermeer.hpp"

int main() {
  gaudi::duchamp::braid_planar_demo_config cfg;
  cfg.show_tp_gradients = true;
  return gaudi::vermeer::vermeer(
      gaudi::duchamp::braid_planar_demo_adapter::create(cfg));
}
