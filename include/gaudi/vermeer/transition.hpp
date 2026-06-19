#pragma once

// Rendering transition strategy for libgaudi demos.
//
// Decision: migrate demo hosts to Vermeer across the board. GaudiGraphics
// remains only as legacy code while existing hosts are replaced.
//
// - Shared simulation/demo code lives in `include/gaudi/duchamp/`.
// - Shared debug API lives in `include/gaudi/geometry_logger.hpp`.
// - Demo hosts should consume snapshots and debug geometry through Vermeer.
//
// Remaining GaudiGraphics projects should be replaced or retired as their demos
// move to Vermeer-native rendering and future render-graph work.

namespace gaudi {
namespace vermeer {
enum class render_backend { vermeer_webgpu };
} // namespace vermeer
} // namespace gaudi
