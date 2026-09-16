#include "gaudi/duchamp/braid_planar_demo.hpp"
#include "gaudi/vermeer/vermeer.hpp"

int main() {
  return gaudi::vermeer::vermeer(
      gaudi::duchamp::braid_planar_demo_adapter::create());
}
