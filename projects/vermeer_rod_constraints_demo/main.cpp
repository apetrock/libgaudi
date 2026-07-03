#include "gaudi/duchamp/rod_constraints_solver_adapter.hpp"
#include "gaudi/vermeer/vermeer.hpp"

int main() {
  return gaudi::vermeer::vermeer(
      gaudi::duchamp::rod_constraints_solver_adapter::create());
}
