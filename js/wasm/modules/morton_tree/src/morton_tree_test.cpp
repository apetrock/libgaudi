#include "gaudi/arp/hash_tree.hpp"
#include "gaudi/arp/morton.hpp"
#include "gaudi/common.h"
#include "gaudi/console_logger.hpp"
#include "gaudi/geometry_logger.hpp"
#include "gaudi/geometry_types.hpp"
#include "gaudi/vec_addendum.h"
#include <emscripten/bind.h>
#include <random>
#include <vector>


// Ensure we use explicit instantiations if enabled
#if defined(GAUDI_USE_EXPLICIT_INSTANTIATIONS) &&                              \
    GAUDI_USE_EXPLICIT_INSTANTIATIONS
// This will use the pre-compiled versions from gaudi_library
#else
// This will instantiate templates inline
#endif

using namespace gaudi;

class MortonTreeTest {
private:
  vec3 near_test_point;
  std::vector<vec3> test_points;
  std::array<vec3, 1> test_point; // Add missing test_point variable
  std::vector<uint32_t> hashes;
  std::vector<index_t> indices;
  std::vector<gaudi::arp::radix_tree_node> internal_nodes;
  std::vector<gaudi::arp::radix_tree_node> leaf_nodes;
  bool auto_visualize = true;

  void redraw_visualizations() {
    geometry_logger::clear();
    if (!test_points.empty() && auto_visualize) {
      draw_points();
      draw_sorted_order_lines();
    }
  }

  void draw_points() {
    // console_logger::debug << "Drawing " << test_points.size() << " points" <<
    // std::endl;
    console_logger::debug << "calculating com for " << test_points.size()
                          << " points" << std::endl;
    auto masses = arp::calc_com<2>(test_points);
    console_logger::debug << "Drawing " << masses.size() << " points"
                          << std::endl;
    for (size_t i = 0; i < masses.size(); i++) {
      const vec3 &p0 = std::get<1>(masses[i]);
      geometry_logger::point(p0, vec4(1.0, 0.1f, 0.1f, 1.0));
    }
    console_logger::debug << "Drawing " << test_points.size() << " lines"
                          << std::endl;
    for (size_t i = 0; i < test_points.size(); i += 2) {
      const vec3 &p0 = test_points[i + 0];
      const vec3 &p1 = test_points[i + 1];
      geometry_logger::line(p0, p1, vec4(0.1, 1.0f, 0.1f, 1.0));
    }
  }

  void draw_sorted_order_lines() {
    if (indices.size() < 2)
      return;

    // Draw lines connecting points in Morton-sorted order
    auto masses = arp::calc_com<2>(test_points);
    for (size_t i = 0; i < indices.size() - 1; i++) {
      const vec3 &current = std::get<1>(masses[indices[i]]);
      const vec3 &next = std::get<1>(masses[indices[i + 1]]);
      geometry_logger::line(current, next, vec4(0, 0, 255, 0.7f)); // Blue lines
    }
  }

public:
  MortonTreeTest() {
    // Initialize with some test points
    test_points = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
                   {0.0, 0.0, 1.0}, {1.0, 1.0, 1.0}, {0.0, 0.0, 1.0},
                   {1.0, 0.0, 1.0}, {1.0, 1.0, 1.0}};
    update_hashes();
  }

  void update_hashes() {
    auto result = arp::make_hash_3d(test_points);
    hashes = result.first;
    indices = result.second;
    redraw_visualizations();
  }

  // Get binary representation for debugging
  std::string get_binary_string(uint32_t hash) {
    return arp::dump_binary(hash);
  }

  // Test hash generation for current points
  int get_hash_count() const { return static_cast<int>(hashes.size()); }

  uint32_t get_hash(int index) const {
    if (index >= 0 && index < static_cast<int>(hashes.size())) {
      return hashes[index];
    }
    return 0;
  }

  int get_index(int index) const {
    if (index >= 0 && index < static_cast<int>(indices.size())) {
      return static_cast<int>(indices[index]);
    }
    return -1;
  }

  // Add a new test point
  void add_point(double x, double y, double z) {
    test_points.push_back({x, y, z});
    update_hashes(); // This triggers redraw for single point additions
  }

  // Clear all points
  void clear_points() {
    test_points.clear();
    hashes.clear();
    indices.clear();
  }

  // Get number of points
  int get_point_count() const { return static_cast<int>(test_points.size()); }

  // Generate random points
  void generate_random_points(int count) {

    // Clear existing points but don't trigger redraw yet
    test_points.clear();
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dist(-0.5, 0.5);
    // N = l * l * l;
    const double scale = 1.0 / pow(double(count), 1.0 / 3.0) * 0.5;
    // Generate all points in batch without triggering redraws
    test_points.reserve(count); // Pre-allocate for efficiency
    for (int i = 0; i < count; i++) {
      const vec3 p0(dist(gen), dist(gen), dist(gen));
      const vec3 p1(dist(gen), dist(gen), dist(gen));
      const vec3 c(dist(gen), dist(gen), dist(gen));
      const vec3 dp = scale * (p1 - p0).normalized();
      test_points.push_back(c - dp);
      test_points.push_back(c + dp);
    }
    // Only update tree and redraw once at the end
    console_logger::debug << "generate_random_points: Calling update_tree..."
                          << std::endl;
    update_tree();
    console_logger::debug << "generate_random_points: update_tree completed"
                          << std::endl;
    console_logger::info << "=== generate_random_points END ===" << std::endl;
  }

  // Generate grid points
  void generate_grid_points(int grid_size) {
    // Clear existing points but don't trigger redraw yet
    test_points.clear();

    const double spacing = 2.0 / (grid_size - 1); // Grid spans [-1, 1]
    const double segment_length =
        spacing * 0.3; // Line segment length relative to grid spacing

    // Generate all points in batch without triggering redraws
    test_points.reserve(grid_size * grid_size * grid_size *
                        2); // 2 points per line segment

    for (int x = 0; x < grid_size; x++) {
      for (int y = 0; y < grid_size; y++) {
        for (int z = 0; z < grid_size; z++) {
          // Calculate center position in [-1, 1] range
          const double cx = -1.0 + x * spacing;
          const double cy = -1.0 + y * spacing;
          const double cz = -1.0 + z * spacing;
          const vec3 center(cx, cy, cz);

          // Create a small line segment centered at this point
          // Use a consistent direction for all segments (e.g., along x-axis)
          const vec3 direction(1.0, 0.0, 0.0);
          const vec3 offset = (segment_length / 2.0) * direction;

          test_points.push_back(center - offset);
          test_points.push_back(center + offset);
        }
      }
    }

    // Only update tree and redraw once at the end
    console_logger::debug << "generate_grid_points: Calling update_tree..."
                          << std::endl;
    update_tree();
    console_logger::debug << "generate_grid_points: update_tree completed"
                          << std::endl;
    console_logger::info << "=== generate_grid_points END ===" << std::endl;
  }

  // Update tree structure
  void update_tree() {
    printf("=== update_tree START ===\n");
    printf("update_tree: Processing %zu points\n", test_points.size());

    printf("update_tree: Calling arp::make_hash_N...\n");
    auto [tree_hashes, tree_indices, internal, leaf] =
        arp::make_hash_N<2>(test_points);

    printf("update_tree: arp::make_hash_tree completed - hashes: %zu, indices: "
           "%zu, internal: %zu, leaf: %zu\n",
           tree_hashes.size(), tree_indices.size(), internal.size(),
           leaf.size());

    hashes = tree_hashes;
    indices = tree_indices;
    internal_nodes = internal;
    leaf_nodes = leaf;

    printf("update_tree: Calling redraw_visualizations...\n");
    redraw_visualizations();
    printf("update_tree: redraw_visualizations completed\n");
    printf("=== update_tree END ===\n");
  }

  // Public method to rebuild hash tree (for explicit API calls)
  bool mk_hash_tree() {
    try {
      printf("=== mk_hash_tree START ===\n");

      if (test_points.empty()) {
        printf("mk_hash_tree: No points available for tree construction\n");
        return false;
      }

      printf("mk_hash_tree: Building tree for %zu points\n",
             test_points.size());

      // Let's manually do what update_tree() does with more debugging
      printf("mk_hash_tree: Calling arp::make_hash_tree...\n");
      auto coms = arp::calc_com<2>(test_points);
      auto centers = std::vector<vec3>(coms.size());
      for (size_t i = 0; i < coms.size(); i++) {
        centers[i] = std::get<1>(coms[i]);
      }
      auto [tree_hashes, tree_indices, internal, leaf] =
          arp::make_hash_tree(centers);

      printf("mk_hash_tree: arp::make_hash_tree returned - hashes: %zu, "
             "indices: %zu, internal: %zu, leaf: %zu\n",
             tree_hashes.size(), tree_indices.size(), internal.size(),
             leaf.size());

      printf("mk_hash_tree: Assigning results...\n");
      hashes = tree_hashes;
      indices = tree_indices;
      internal_nodes = internal;
      leaf_nodes = leaf;

      printf("mk_hash_tree: Calling redraw_visualizations...\n");
      redraw_visualizations();

      printf("mk_hash_tree: Tree built successfully - %zu hashes, %zu internal "
             "nodes, %zu leaf nodes\n",
             hashes.size(), internal_nodes.size(), leaf_nodes.size());
      printf("=== mk_hash_tree END ===\n");
      return true;
    } catch (const std::exception &e) {
      printf("mk_hash_tree: Exception caught: %s\n", e.what());
      return false;
    } catch (...) {
      printf("mk_hash_tree: Unknown exception caught\n");
      return false;
    }
  }

  // Log hierarchy visualization
  void log_hierarchy() {
    try {
      geometry_logger::clear();
      arp::log_hierarchy<2>(test_points, indices, internal_nodes, leaf_nodes);
    } catch (const std::exception &e) {
      printf("log_hierarchy: Exception caught: %s\n", e.what());
      throw; // Re-throw to get stack trace in JS
    } catch (...) {
      printf("log_hierarchy: Unknown exception caught\n");
      throw; // Re-throw to get stack trace in JS
    }
  }

  // Log BVH visualization
  void log_bvh() {
    try {
      if (test_points.empty()) {
        printf("log_bvh: No points available\n");
        return;
      }

      if (hashes.empty() || indices.empty()) {
        printf("log_bvh: Tree data not initialized. Points: %zu, Hashes: %zu, "
               "Indices: %zu\n",
               test_points.size(), hashes.size(), indices.size());
        return;
      }

      printf("log_bvh: Logging BVH with %zu points, %zu hashes, %zu internal "
             "nodes, %zu leaf nodes\n",
             test_points.size(), hashes.size(), internal_nodes.size(),
             leaf_nodes.size());

      geometry_logger::clear();
      arp::log_bvh<2>(test_points, indices, internal_nodes, leaf_nodes);
      printf("log_bvh: Successfully logged BVH\n");
    } catch (const std::exception &e) {
      printf("log_bvh: Exception caught: %s\n", e.what());
      throw; // Re-throw to get stack trace in JS
    } catch (...) {
      printf("log_bvh: Unknown exception caught\n");
      throw; // Re-throw to get stack trace in JS
    }
  }

  void log_nearest(float t) {
    try {
      if (test_points.empty()) {
        printf("log_nearest: No points available\n");
        return;
      }
      
      // Clear geometry logger
      geometry_logger::clear();
      
      // Update test point with complex trajectory
      const float C = 0.1f; // Speed constant
      near_test_point = vec3(
        sin(C * t), 
        cos(C * t), 
        sin(C * t) * cos(C * t)
      );
      test_point = {near_test_point};
      
      // Log the test point
      geometry_logger::point(near_test_point, vec4(1.0, 1.0, 0.0, 1.0)); // Yellow test point
      
      const auto bvh_result =
          arp::make_bvh<2>(test_points, indices, internal_nodes, leaf_nodes);
      auto result = arp::getNearest<std::array<vec3, 1>, 2>(
          test_point, test_points, indices, internal_nodes, leaf_nodes,
          bvh_result, 10000.0,
          [](const std::array<vec3, 1> &t_verts, const arp::near_array<2> &s_verts) {
            const vec3 &xA = t_verts[0];
            const vec3 &xB0 = s_verts[0];
            const vec3 &xB1 = s_verts[1];
            gaudi::real d = va::distance_from_line(xB0, xB1, xA);
            return d;
          });
      
      printf("log_nearest: Result: %zu\n", result.size());
      
      // Draw lines to the nearest leaf nodes
      auto masses = arp::calc_com<2>(test_points);
      for (size_t i = 0; i < result.size(); i++) {
        int leaf_id = result[i];
        if (leaf_id >= 0 && leaf_id < static_cast<int>(masses.size())) {
          const vec3 &leaf_point = std::get<1>(masses[leaf_id]);
          //const vec3 & leaf_point = test_points[leaf_id];
          geometry_logger::line(near_test_point, leaf_point, vec4(1.0, 0.0, 1.0, 0.8f)); // Magenta lines
        }
      }
    } catch (const std::exception &e) {
      printf("log_nearest: Exception caught: %s\n", e.what());
      throw; // Re-throw to get stack trace in JS
    } catch (...) {
      printf("log_nearest: Unknown exception caught\n");
      throw; // Re-throw to get stack trace in JS
    }
  }

  // Visualization control methods
  void set_auto_visualize(bool enable) {
    auto_visualize = enable;
    if (enable) {
      redraw_visualizations();
    } else {
      geometry_logger::clear();
    }
  }

  bool get_auto_visualize() const { return auto_visualize; }

  void clear_visualizations() { geometry_logger::clear(); }

  void draw_points_only() {
    geometry_logger::clear();
    draw_points();
  }

  void draw_sorted_lines_only() {
    geometry_logger::clear();
    draw_sorted_order_lines();
  }

  void draw_all_visualizations() { redraw_visualizations(); }

  // Logger API methods for compatibility with GaudiLoggerRenderer
  int get_line_count() const {
    return geometry_logger::get_lines().size() / 2; // Each line has 2 points
  }

  int get_point_count_logger() const {
    return geometry_logger::get_points().size();
  }

  void get_line(int index, double *start, double *end, double *color) const {
    const auto &lines = geometry_logger::get_lines();
    const auto &colors = geometry_logger::get_line_colors();

    if (index * 2 + 1 < lines.size()) {
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

    if (index < points.size()) {
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

  // Get point coordinates for visualization
  double get_point_x(int index) const {
    if (index >= 0 && index < static_cast<int>(test_points.size())) {
      return test_points[index][0];
    }
    return 0.0;
  }

  double get_point_y(int index) const {
    if (index >= 0 && index < static_cast<int>(test_points.size())) {
      return test_points[index][1];
    }
    return 0.0;
  }

  double get_point_z(int index) const {
    if (index >= 0 && index < static_cast<int>(test_points.size())) {
      return test_points[index][2];
    }
    return 0.0;
  }

  // Get leaf node point coordinates
  double get_leaf_point_x(int index) const {
    if (test_points.empty()) return 0.0;
    auto masses = arp::calc_com<2>(test_points);
    if (index >= 0 && index < static_cast<int>(masses.size())) {
      return std::get<1>(masses[index]).x();
    }
    return 0.0;
  }

  double get_leaf_point_y(int index) const {
    if (test_points.empty()) return 0.0;
    auto masses = arp::calc_com<2>(test_points);
    if (index >= 0 && index < static_cast<int>(masses.size())) {
      return std::get<1>(masses[index]).y();
    }
    return 0.0;
  }

  double get_leaf_point_z(int index) const {
    if (test_points.empty()) return 0.0;
    auto masses = arp::calc_com<2>(test_points);
    if (index >= 0 && index < static_cast<int>(masses.size())) {
      return std::get<1>(masses[index]).z();
    }
    return 0.0;
  }

  // Get number of leaf nodes
  int get_leaf_count() const {
    if (test_points.empty()) return 0;
    auto masses = arp::calc_com<2>(test_points);
    return static_cast<int>(masses.size());
  }
};

EMSCRIPTEN_BINDINGS(morton_tree_test) {
  emscripten::class_<MortonTreeTest>("MortonTreeTest")
      .constructor<>()
      .function("get_binary_string", &MortonTreeTest::get_binary_string)
      .function("get_hash_count", &MortonTreeTest::get_hash_count)
      .function("get_hash", &MortonTreeTest::get_hash)
      .function("get_index", &MortonTreeTest::get_index)
      .function("add_point", &MortonTreeTest::add_point)
      .function("clear_points", &MortonTreeTest::clear_points)
      .function("get_point_count", &MortonTreeTest::get_point_count)
      .function("generate_random_points",
                &MortonTreeTest::generate_random_points)
      .function("generate_grid_points", &MortonTreeTest::generate_grid_points)
      .function("mk_hash_tree", &MortonTreeTest::mk_hash_tree)
      .function("log_hierarchy", &MortonTreeTest::log_hierarchy)
      .function("log_bvh", &MortonTreeTest::log_bvh)
      .function("log_nearest", &MortonTreeTest::log_nearest)
      .function("get_point_x", &MortonTreeTest::get_point_x)
      .function("get_point_y", &MortonTreeTest::get_point_y)
      .function("get_point_z", &MortonTreeTest::get_point_z)
      .function("get_leaf_point_x", &MortonTreeTest::get_leaf_point_x)
      .function("get_leaf_point_y", &MortonTreeTest::get_leaf_point_y)
      .function("get_leaf_point_z", &MortonTreeTest::get_leaf_point_z)
      .function("get_leaf_count", &MortonTreeTest::get_leaf_count)
      // Visualization controls
      .function("set_auto_visualize", &MortonTreeTest::set_auto_visualize)
      .function("get_auto_visualize", &MortonTreeTest::get_auto_visualize)
      .function("clear_visualizations", &MortonTreeTest::clear_visualizations)
      .function("draw_points_only", &MortonTreeTest::draw_points_only)
      .function("draw_sorted_lines_only",
                &MortonTreeTest::draw_sorted_lines_only)
      .function("draw_all_visualizations",
                &MortonTreeTest::draw_all_visualizations)
      // Logger API compatibility
      .function("get_line_count", &MortonTreeTest::get_line_count)
      .function("get_point_count_logger",
                &MortonTreeTest::get_point_count_logger)
      .function("get_line", &MortonTreeTest::get_line,
                emscripten::allow_raw_pointers())
      .function("get_point_logger", &MortonTreeTest::get_point_logger,
                emscripten::allow_raw_pointers());
}