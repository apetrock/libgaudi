#ifndef __GAUDI_DUCHAMP_RX_RESEARCH_NOTES_HPP__
#define __GAUDI_DUCHAMP_RX_RESEARCH_NOTES_HPP__

// Intentionally no runtime API: design reminders for future R&D (see project plan).
//
// 1) Semi-Lagrangian on a tri surface: back-tracing along characteristics needs geodesic
//    (or approximate shortest) paths; this repo has heat / geodesic solvers (e.g.
//    heat_geodesics_test, related code). Wiring SL to that path is a dedicated follow-up.
//
// 2) Madelung / (ρ, S) or polar splits vs. the coupled 2n Crank–Nicolson in
//    gaudi/kusama/complex_laplacian.hpp: alternative discretizations, not drop-in
//    replacements, until stability is analyzed.

#endif
