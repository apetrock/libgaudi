#include "gaudi/duchamp/monge_fairing_demo.hpp"

// Measures body spin ωz against the shape-operator tilt on one tube of the
// built-in torus, after one smoothing pass and before remesh.
//
// The Hopf landing is then turned about N by that spoke's ωz when
// hopf_sweep::spin is set (hopf_project_spin). |ωz| above 0.5 is left
// unturned. Set the second argument false to undo it.

int main() {
  auto demo = gaudi::duchamp::monge_fairing_demo::create();
  demo->relax();
  demo->report_tube_spin(std::cout);
  return 0;
}
