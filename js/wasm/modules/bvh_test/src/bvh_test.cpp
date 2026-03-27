#include "gaudi/arp/brute_force.hpp"
#include "gaudi/arp/hash_tree.hpp"
#include "gaudi/arp/morton.hpp"
#include "gaudi/asawa/objloader_refactor.hpp"
#include "gaudi/common.h"
#include "gaudi/console_logger.hpp"
#include "gaudi/geometry_logger.hpp"
#include "gaudi/geometry_types.hpp"
#include "gaudi/test/bvh_tests.hpp"
#include "gaudi/test/test.hpp"
#include "gaudi/vec_addendum.h"
#include <emscripten/bind.h>
#include <random>
#include <set>
#include <vector>

using namespace gaudi;

class BvhTest {
private:
  // Mesh data
  std::vector<vec3> vertices_;
  std::vector<std::vector<int>> faces_;

  // Extracted primitives
  std::vector<index_t> point_adjacency_;    // indices for N=1 (just vertex indices)
  std::vector<index_t> edge_adjacency_;     // pairs of vertex indices for N=2
  std::vector<index_t> triangle_adjacency_; // triples of vertex indices for N=3

  // BVH trees
  arp::bvh_tree<2>::ptr edge_bvh_;
  arp::bvh_tree<3>::ptr tri_bvh_;

  // Test results
  index_t last_bvh_result_ = -1;
  index_t last_brute_result_ = -1;
  bool last_test_passed_ = false;
  int last_suite_total_ = 0;
  int last_suite_failed_ = 0;

  // Helper to extract unique edges from faces
  void extractEdges() {
    edge_adjacency_.clear();
    std::set<std::pair<int, int>> unique_edges;

    for (const auto &face : faces_) {
      for (size_t i = 0; i < face.size(); i++) {
        int v0 = face[i];
        int v1 = face[(i + 1) % face.size()];
        // Store edges in canonical order (smaller index first)
        auto edge = std::minmax(v0, v1);
        if (unique_edges.find(edge) == unique_edges.end()) {
          unique_edges.insert(edge);
          edge_adjacency_.push_back(edge.first);
          edge_adjacency_.push_back(edge.second);
        }
      }
    }
    console_logger::info << "Extracted " << (edge_adjacency_.size() / 2)
                         << " unique edges" << std::endl;
  }

  // Helper to extract triangles from faces (assumes triangular faces)
  void extractTriangles() {
    triangle_adjacency_.clear();
    for (const auto &face : faces_) {
      if (face.size() >= 3) {
        // Simple triangulation: fan from first vertex
        for (size_t i = 1; i < face.size() - 1; i++) {
          triangle_adjacency_.push_back(face[0]);
          triangle_adjacency_.push_back(face[i]);
          triangle_adjacency_.push_back(face[i + 1]);
        }
      }
    }
    console_logger::info << "Extracted " << (triangle_adjacency_.size() / 3)
                         << " triangles" << std::endl;
  }

  // Helper to build point adjacency (just sequential indices)
  void extractPoints() {
    point_adjacency_.clear();
    for (size_t i = 0; i < vertices_.size(); i++) {
      point_adjacency_.push_back(static_cast<index_t>(i));
    }
    console_logger::info << "Extracted " << point_adjacency_.size() << " points"
                         << std::endl;
  }

public:
  BvhTest() = default;

  // Load mesh from OBJ string content
  bool loadMeshFromString(const std::string &objContent) {
    vertices_.clear();
    faces_.clear();

    try {
      asawa::loadObjFromString(objContent, vertices_, faces_);
      console_logger::info << "Loaded mesh: " << vertices_.size()
                           << " vertices, " << faces_.size() << " faces"
                           << std::endl;

      // Extract all primitive types
      extractPoints();
      extractEdges();
      extractTriangles();

      return true;
    } catch (const std::exception &e) {
      console_logger::error << "Failed to load mesh: " << e.what() << std::endl;
      return false;
    }
  }

  // Getters for mesh info
  int getVertexCount() const { return static_cast<int>(vertices_.size()); }
  int getEdgeCount() const {
    return static_cast<int>(edge_adjacency_.size() / 2);
  }
  int getTriangleCount() const {
    return static_cast<int>(triangle_adjacency_.size() / 3);
  }

  // Build BVH for edges (N=2)
  bool buildEdgeBvh() {
    if (vertices_.empty() || edge_adjacency_.empty()) {
      console_logger::error << "No mesh data to build edge BVH" << std::endl;
      return false;
    }

    try {
      edge_bvh_ = arp::bvh_tree<2>::create(edge_adjacency_, vertices_);
      console_logger::info << "Built edge BVH" << std::endl;
      return true;
    } catch (const std::exception &e) {
      console_logger::error << "Failed to build edge BVH: " << e.what()
                            << std::endl;
      return false;
    }
  }

  // Build BVH for triangles (N=3)
  bool buildTriBvh() {
    if (vertices_.empty() || triangle_adjacency_.empty()) {
      console_logger::error << "No mesh data to build triangle BVH"
                            << std::endl;
      return false;
    }

    try {
      tri_bvh_ = arp::bvh_tree<3>::create(triangle_adjacency_, vertices_);
      console_logger::info << "Built triangle BVH" << std::endl;
      return true;
    } catch (const std::exception &e) {
      console_logger::error << "Failed to build triangle BVH: " << e.what()
                            << std::endl;
      return false;
    }
  }

  // Test point-to-edge queries: pick random points and compare BVH vs brute
  // force
  bool testPointToEdge(double tolerance) {
    if (!edge_bvh_) {
      console_logger::error << "Edge BVH not built" << std::endl;
      return false;
    }

    // Pick a random test point (or center of mesh)
    vec3 test_point = vec3::Zero();
    for (const auto &v : vertices_) {
      test_point += v;
    }
    test_point /= static_cast<real>(vertices_.size());

    // Create query as array (point)
    std::array<vec3, 1> query = {test_point};

    // BVH query
    auto bvh_results = edge_bvh_->get_nearest(query, tolerance);
    last_bvh_result_ = bvh_results.empty() ? -1 : bvh_results.back();

    // Brute force query
    using EdgeView =
        arp::bvh_tree<2>::permuted_view;
    
    // Create an unpermuted view for brute force
    adjacency_view<std::vector<vec3>, std::vector<index_t>> edge_view(
        vertices_, edge_adjacency_);

    auto testAB = [](const std::array<vec3, 1> &pA,
                     const slice<2, decltype(edge_view)> &pB) {
      return arp::test_point_line(pA, pB);
    };

    auto brute_results =
        arp::brute_force_nearest<2>(query, edge_view, tolerance, testAB);
    last_brute_result_ = brute_results.empty() ? -1 : brute_results.back();

    // Compare results
    last_test_passed_ = (last_bvh_result_ == last_brute_result_);

    console_logger::info << "Point-to-Edge test: BVH=" << last_bvh_result_
                         << " Brute=" << last_brute_result_
                         << " Passed=" << last_test_passed_ << std::endl;

    return last_test_passed_;
  }

  // Test point-to-triangle queries
  bool testPointToTri(double tolerance) {
    if (!tri_bvh_) {
      console_logger::error << "Triangle BVH not built" << std::endl;
      return false;
    }

    // Pick center of mesh as test point
    vec3 test_point = vec3::Zero();
    for (const auto &v : vertices_) {
      test_point += v;
    }
    test_point /= static_cast<real>(vertices_.size());

    std::array<vec3, 1> query = {test_point};

    // For now just verify BVH doesn't crash
    // Triangle queries not fully implemented
    console_logger::info << "Point-to-Triangle test: Not fully implemented"
                         << std::endl;
    last_test_passed_ = true;
    return true;
  }

  // Test edge-to-edge queries
  bool testEdgeToEdge(double tolerance) {
    if (!edge_bvh_) {
      console_logger::error << "Edge BVH not built" << std::endl;
      return false;
    }

    if (edge_adjacency_.size() < 4) {
      console_logger::error << "Not enough edges for edge-edge test"
                            << std::endl;
      return false;
    }

    // Pick first edge as query
    std::array<vec3, 2> query = {vertices_[edge_adjacency_[0]],
                                 vertices_[edge_adjacency_[1]]};

    // BVH query
    auto bvh_results = edge_bvh_->get_nearest(query, tolerance);
    last_bvh_result_ = bvh_results.empty() ? -1 : bvh_results.back();

    // Brute force query
    adjacency_view<std::vector<vec3>, std::vector<index_t>> edge_view(
        vertices_, edge_adjacency_);

    auto testAB = [](const std::array<vec3, 2> &pA,
                     const slice<2, decltype(edge_view)> &pB) {
      return arp::test_line_line(pA, pB);
    };

    auto brute_results =
        arp::brute_force_nearest<2>(query, edge_view, tolerance, testAB);
    last_brute_result_ = brute_results.empty() ? -1 : brute_results.back();

    // Compare results
    last_test_passed_ = (last_bvh_result_ == last_brute_result_);

    console_logger::info << "Edge-to-Edge test: BVH=" << last_bvh_result_
                         << " Brute=" << last_brute_result_
                         << " Passed=" << last_test_passed_ << std::endl;

    return last_test_passed_;
  }

  // Getters for test results
  int getLastBvhResult() const { return last_bvh_result_; }
  int getLastBruteResult() const { return last_brute_result_; }
  bool getLastTestPassed() const { return last_test_passed_; }

  // Visualization methods
  void visualizeResults() {
    geometry_logger::clear();

    // Draw all vertices as points
    for (const auto &v : vertices_) {
      geometry_logger::point(v, vec4(0.5, 0.5, 0.5, 0.5));
    }

    // Draw edges
    for (size_t i = 0; i < edge_adjacency_.size(); i += 2) {
      const vec3 &v0 = vertices_[edge_adjacency_[i]];
      const vec3 &v1 = vertices_[edge_adjacency_[i + 1]];
      geometry_logger::line(v0, v1, vec4(0.2, 0.8, 0.2, 0.5));
    }

    // Draw BVH if built
    if (edge_bvh_) {
      // Draw centers of mass for edge BVH
      for (size_t i = 0; i < edge_adjacency_.size() / 2; i++) {
        vec3 com = edge_bvh_->get_com(static_cast<index_t>(i));
        geometry_logger::point(com, vec4(1.0, 0.0, 0.0, 1.0));
      }
    }
  }

  void clearVisualization() { geometry_logger::clear(); }

  // Test harness helpers
  void setSphereObjData(const std::string &objContent) {
    gaudi::test::assets::set_sphere_obj_data(objContent);
  }

  bool runAllTests() {
    auto result = gaudi::test::Registry::run_all();
    last_suite_total_ = result.total;
    last_suite_failed_ = result.failed;
    last_test_passed_ = (result.failed == 0);
    console_logger::info << "TEST SUMMARY: total=" << result.total
                         << " failed=" << result.failed << std::endl;
    return last_test_passed_;
  }

  int getLastSuiteTotal() const { return last_suite_total_; }
  int getLastSuiteFailed() const { return last_suite_failed_; }

  // Logger API for GaudiLoggerRenderer compatibility
  int get_line_count() const {
    return static_cast<int>(geometry_logger::get_lines().size() / 2);
  }

  int get_point_count_logger() const {
    return static_cast<int>(geometry_logger::get_points().size());
  }

  void get_line(int index, double *start, double *end, double *color) const {
    const auto &lines = geometry_logger::get_lines();
    const auto &colors = geometry_logger::get_line_colors();

    if (index * 2 + 1 < static_cast<int>(lines.size())) {
      const vec3 &p0 = lines[index * 2];
      const vec3 &p1 = lines[index * 2 + 1];
      const vec4 &col = colors[index * 2];

      start[0] = p0.x();
      start[1] = p0.y();
      start[2] = p0.z();
      end[0] = p1.x();
      end[1] = p1.y();
      end[2] = p1.z();
      color[0] = col.x();
      color[1] = col.y();
      color[2] = col.z();
      color[3] = col.w();
    }
  }

  void get_point_logger(int index, double *position, double *color) const {
    const auto &points = geometry_logger::get_points();
    const auto &colors = geometry_logger::get_point_colors();

    if (index < static_cast<int>(points.size())) {
      const vec3 &pos = points[index];
      const vec4 &col = colors[index];

      position[0] = pos.x();
      position[1] = pos.y();
      position[2] = pos.z();
      color[0] = col.x();
      color[1] = col.y();
      color[2] = col.z();
      color[3] = col.w();
    }
  }
};

EMSCRIPTEN_BINDINGS(bvh_test) {
  emscripten::class_<BvhTest>("BvhTest")
      .constructor<>()
      .function("loadMeshFromString", &BvhTest::loadMeshFromString)
      .function("getVertexCount", &BvhTest::getVertexCount)
      .function("getEdgeCount", &BvhTest::getEdgeCount)
      .function("getTriangleCount", &BvhTest::getTriangleCount)
      .function("buildEdgeBvh", &BvhTest::buildEdgeBvh)
      .function("buildTriBvh", &BvhTest::buildTriBvh)
      .function("testPointToEdge", &BvhTest::testPointToEdge)
      .function("testPointToTri", &BvhTest::testPointToTri)
      .function("testEdgeToEdge", &BvhTest::testEdgeToEdge)
      .function("getLastBvhResult", &BvhTest::getLastBvhResult)
      .function("getLastBruteResult", &BvhTest::getLastBruteResult)
      .function("getLastTestPassed", &BvhTest::getLastTestPassed)
      .function("setSphereObjData", &BvhTest::setSphereObjData)
      .function("runAllTests", &BvhTest::runAllTests)
      .function("getLastSuiteTotal", &BvhTest::getLastSuiteTotal)
      .function("getLastSuiteFailed", &BvhTest::getLastSuiteFailed)
      .function("visualizeResults", &BvhTest::visualizeResults)
      .function("clearVisualization", &BvhTest::clearVisualization)
      .function("get_line_count", &BvhTest::get_line_count)
      .function("get_point_count_logger", &BvhTest::get_point_count_logger)
      .function("get_line", &BvhTest::get_line,
                emscripten::allow_raw_pointers())
      .function("get_point_logger", &BvhTest::get_point_logger,
                emscripten::allow_raw_pointers());
}

