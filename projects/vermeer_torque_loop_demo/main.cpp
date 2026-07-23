#include "gaudi/duchamp/torque_loop.hpp"
#include "gaudi/vermeer/vermeer.hpp"

int main() {
  return gaudi::vermeer::vermeer(
      gaudi::duchamp::torque_loop_adapter::create());
}
