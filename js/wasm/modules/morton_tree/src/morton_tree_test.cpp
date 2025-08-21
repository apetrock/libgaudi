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
  std::vector<index_t> adjacency;
  std::vector<gaudi::arp::radix_tree_node> internal_nodes;
  std::vector<gaudi::arp::radix_tree_node> leaf_nodes;
  bool auto_visualize = true;

  // Helper method to generate random Gaussian points
  std::vector<vec3> generate_random_gaussian_points(int count) {
    std::vector<vec3> points;
    points.reserve(count);
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dist(-0.5, 0.5);
    const double scale = 1.0 / pow(double(count), 1.0 / 3.0) * 0.5;
    
    for (int i = 0; i < count; i++) {
      const vec3 p0(dist(gen), dist(gen), dist(gen));
      const vec3 p1(dist(gen), dist(gen), dist(gen));
      const vec3 c(dist(gen), dist(gen), dist(gen));
      const vec3 dp = scale * (p1 - p0).normalized();
      points.push_back(c - dp);
      points.push_back(c + dp);
    }
    
    return points;
  }

  // Helper method to generate individual random Gaussian points
  std::vector<vec3> generate_random_control_points(int count) {
    std::vector<vec3> points;
    points.reserve(count);
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dist(-0.5, 0.5);
    const double scale = 0.8; // Scale factor to fit in [-1, 1] range
    
    for (int i = 0; i < count; i++) {
      const vec3 point(scale * dist(gen), scale * dist(gen), scale * dist(gen));
      points.push_back(point);
    }
    
    return points;
  }

  // Helper method to evaluate cubic Bezier curve
  vec3 evaluate_bezier(const vec3& p0, const vec3& p1, const vec3& p2, const vec3& p3, double t) {
    const double t2 = t * t;
    const double t3 = t2 * t;
    const double mt = 1.0 - t;
    const double mt2 = mt * mt;
    const double mt3 = mt2 * mt;
    
    return mt3 * p0 + 3.0 * mt2 * t * p1 + 3.0 * mt * t2 * p2 + t3 * p3;
  }

  // Helper method to generate farthest-first permutation using stack
  std::vector<int> generate_farthest_first_permutation(const std::vector<vec3>& points, int sample_size = 5) {
    int num_points = static_cast<int>(points.size());
    std::vector<int> permutation(num_points, -1);
    
    // Create a stack of randomly shuffled indices
    std::vector<int> index_stack;
    index_stack.reserve(num_points);
    for (int i = 0; i < num_points; i++) {
      index_stack.push_back(i);
    }
    
    // Shuffle the stack
    std::random_device rd;
    std::mt19937 gen(rd());
    std::shuffle(index_stack.begin(), index_stack.end(), gen);
    
    // Lambda to find farthest point from a given point using the stack
    auto find_farthest_point = [&](const vec3& point, const std::vector<vec3>& points, 
                                   std::vector<int>& stack) -> int {
      if (stack.empty()) return -1;
      
      double max_dist = -1.0;
      int max_idx = -1;
      int max_stack_pos = -1;
      
      // Sample up to sample_size points from the stack
      int samples_to_take = std::min(sample_size, static_cast<int>(stack.size()));
      
      for (int i = 0; i < samples_to_take; i++) {
        int stack_idx = stack.size() - 1 - i; // Take from end of stack
        int point_idx = stack[stack_idx];
        const vec3& candidate = points[point_idx];
        
        double dist = (point - candidate).squaredNorm();
        if (dist > max_dist) {
          max_dist = dist;
          max_idx = point_idx;
          max_stack_pos = stack_idx;
        }
      }
      
      // Remove the farthest point from stack
      if (max_stack_pos >= 0) {
        stack.erase(stack.begin() + max_stack_pos);
      }
      
      return max_idx;
    };
    
    // Generate permutation by finding farthest point for each position
    for (int i = 0; i < num_points; i++) {
      permutation[i] = find_farthest_point(points[i], points, index_stack);
    }
    
    return permutation;
  }

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
    console_logger::debug << "draw_points: test_points.size(): " << test_points.size() << std::endl;
                   
    console_logger::debug << "adjacency: " << adjacency.size() << std::endl;
    console_logger::debug << "indices: " << indices.size() << std::endl;
    console_logger::debug << "building permuted points" << std::endl;
    //auto p_test_points = permuted(test_points, permuted(adjacency, spread<2, decltype(indices)>(indices)));
    auto p_test_points = permuted_adjacency_view<2, decltype(test_points), decltype(adjacency)>(test_points, adjacency, indices);

    auto masses = arp::calc_com<2>(p_test_points);
    
    for (size_t i = 0; i < masses.size(); i++) {
      const vec3 &p0 = std::get<1>(masses[i]);
      geometry_logger::point(p0, vec4(1.0, 0.1f, 0.1f, 1.0));
    }

    for (size_t i = 0; i < p_test_points.size(); i += 2) {
      const vec3 &p0 = p_test_points[i + 0];
      const vec3 &p1 = p_test_points[i + 1];
      geometry_logger::line(p0, p1, vec4(0.1, 1.0f, 0.1f, 1.0));
    }
  }

  void draw_sorted_order_lines() {
    if (indices.size() < 2)
      return;

    // Draw lines connecting points in Morton-sorted order
    console_logger::debug << "draw_sorted_order_lines: test_points.size(): " << test_points.size() << std::endl;
    console_logger::debug << "draw_sorted_order_lines: indices.size(): " << indices.size() << std::endl;
    auto p_test_points = permuted_adjacency_view<2, decltype(test_points), decltype(adjacency)>(test_points, adjacency, indices);
    auto coms = arp::calc_com<2>(p_test_points);
    for(size_t i = 0; i < coms.size() - 1; i++) {
      const vec3 &p0 = std::get<1>(coms[i]);
      const vec3 &p1 = std::get<1>(coms[i+1]);
      geometry_logger::line(p0, p1, vec4(0, 0, 255, 0.7f));
    }
  }

public:
  MortonTreeTest() {
    // Initialize with some test points
    test_points = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
                   {0.0, 0.0, 1.0}, {1.0, 1.0, 0.0}, {1.0, 0.0, 1.0},
                   {0.0, 1.0, 1.0}, {1.0, 1.0, 1.0}};
    adjacency = {0, 1, 2, 3, 4, 5, 6, 7};
    update_hashes();
  }

  void update_hashes() {

    console_logger::debug << "update_hashes: " << test_points.size() << std::endl;
    auto p_test_points = adjacency_view<decltype(test_points), decltype(adjacency)>(test_points, adjacency);
    auto [tree_hashes, tree_indices, tree_internal, tree_leaf] = arp::make_hash<2>(p_test_points);
    console_logger::debug << "done building hash tree" << std::endl;
    this->hashes = tree_hashes;
    this->indices = tree_indices;
    this->internal_nodes = tree_internal;
    this->leaf_nodes = tree_leaf;
    if(indices.size() < 128) {
      //dump indices:
      console_logger::debug << "update_hashes: indices: ";
      for(size_t i = 0; i < indices.size(); i++) {  
        console_logger::debug << indices[i] << " ";
      }
      console_logger::debug << std::endl;
    }
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
    adjacency.clear();
    // Use the helper method to generate random points
    test_points = generate_random_gaussian_points(count);
    for(size_t i = 0; i < test_points.size(); i++) {
      adjacency.push_back(i);
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
    adjacency.clear();
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
          const index_t idx = x * grid_size * grid_size + y * grid_size + z;
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
          adjacency.push_back(idx * 2 + 0);
          adjacency.push_back(idx * 2 + 1);
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

  // Generate trefoil knot points
  void generate_trefoil_knot(int num_segments) {
    // Clear existing points but don't trigger redraw yet
    test_points.clear();
    adjacency.clear();
    // Trefoil knot parametric equations
    // x = (2 + cos(3t)) * cos(2t)
    // y = (2 + cos(3t)) * sin(2t) 
    // z = sin(3t)
    // where t goes from 0 to 2π
    
    const double scale = 0.3; // Scale factor to fit in [-1, 1] range
    const double t_step = 2.0 * M_PI / num_segments;
    
    // Generate all points in batch without triggering redraws
    test_points.reserve(num_segments * 2); // 2 points per line segment
    
    for (int i = 0; i < num_segments; i++) {
      const double t1 = i * t_step;
      const double x1 = scale * (2.0 + cos(3.0 * t1)) * cos(2.0 * t1);
      const double y1 = scale * (2.0 + cos(3.0 * t1)) * sin(2.0 * t1);
      const double z1 = scale * sin(3.0 * t1);      
      test_points.push_back(vec3(x1, y1, z1));
      if(i > 0) {
        int i0 = i;
        int i1 = (i + 1) % num_segments;
        adjacency.push_back(i0);
        adjacency.push_back(i1);
      }
    }

    // Only update tree and redraw once at the end
    console_logger::debug << "generate_trefoil_knot: Calling update_tree..."
                          << std::endl;
    update_tree();
    console_logger::debug << "generate_trefoil_knot: update_tree completed"
                          << std::endl;
    console_logger::info << "=== generate_trefoil_knot END ===" << std::endl;
  }

  // Generate general knot points
  void generate_knot(int k, int num_segments) {
    // Clear existing points but don't trigger redraw yet
    test_points.clear();
    adjacency.clear();
    // General knot parametric equations
    // x = cos(u) [ 2 - cos(2 u/(2 k + 1)) ]
    // y = sin(u) [ 2 - cos(2 u/(2 k + 1)) ]
    // z = -sin(2 u/(2 k + 1))
    // where 0 < u < (4 k + 2) pi
    
    const double scale = 0.3; // Scale factor to fit in [-1, 1] range
    const double u_max = (4.0 * k + 2.0) * M_PI;
    const double u_step = u_max / num_segments;
    
    // Generate all points in batch without triggering redraws
    test_points.reserve(num_segments * 2); // 2 points per line segment
    
    for (int i = 0; i < num_segments; i++) {
      const double u1 = i * u_step;

      // Calculate first point
      const double cos_term1 = cos(2.0 * u1 / (2.0 * k + 1.0));
      const double x1 = scale * cos(u1) * (2.0 - cos_term1);
      const double y1 = scale * sin(u1) * (2.0 - cos_term1);
      const double z1 = scale * (-sin(2.0 * u1 / (2.0 * k + 1.0)));
      
      // Calculate second point
      
      test_points.push_back(vec3(x1, y1, z1));
      if(i > 0) {
        int i0 = i;
        int i1 = (i + 1) % num_segments;
        adjacency.push_back(i0);
        adjacency.push_back(i1);
      }
    }

    // Only update tree and redraw once at the end
    console_logger::debug << "generate_knot: Calling update_tree..."
                          << std::endl;
    update_tree();
    console_logger::debug << "generate_knot: update_tree completed"
                          << std::endl;
    console_logger::info << "=== generate_knot END ===" << std::endl;
  }

  // Generate random knot using Bezier splines
  void generate_random_knot(int num_control_points, int segments_per_chord) {
    // Clear existing points but don't trigger redraw yet
    test_points.clear();
    adjacency.clear();
    // Generate random control points using the helper function
    std::vector<vec3> control_points = generate_random_control_points(num_control_points);
    
    // Generate farthest-first permutation
    std::vector<int> permutation = generate_farthest_first_permutation(control_points);
    
    // Generate knot points by connecting control points with Bezier splines using permutation
    test_points.reserve(num_control_points * segments_per_chord); // 2 points per line segment
    
    for (int i = 0; i < num_control_points; i+=2) {
      // Get indices with wrap-around
      int im1 = (i - 1 + num_control_points) % num_control_points;
      int im0 = i;
      int ip0 = (i + 1) % num_control_points;
      int ip1 = (i + 2) % num_control_points;
      
      // Get control points
      const vec3& pm1 = control_points[permutation[im1]];
      const vec3& pm0 = control_points[permutation[im0]];
      const vec3& pp0 = control_points[permutation[ip0]];
      const vec3& pp1 = control_points[permutation[ip1]];
      
      // Calculate Bezier control points for this chord
      const vec3 pc0 = 0.5 * (pm1 + pm0);
      const vec3 pc1 = pm0;
      const vec3 pc2 = pp0;
      const vec3 pc3 = 0.5 * (pp0 + pp1);
      
      // Generate segments for this chord
      for (int j = 0; j < segments_per_chord; j++) {
        const double t1 = static_cast<double>(j) / segments_per_chord;
        
        // Evaluate Bezier curve at t1 and t2
        const vec3 p1 = evaluate_bezier(pc0, pc1, pc2, pc3, t1);
        test_points.push_back(p1);
      }
    }

    for(size_t i = 0; i < test_points.size(); i++) {
      int i0 = i;
      int i1 = (i + 1) % test_points.size();
      adjacency.push_back(i0);
      adjacency.push_back(i1);
    }

    // Only update tree and redraw once at the end
    console_logger::debug << "generate_random_knot: Calling update_tree..."
                          << std::endl;
    update_tree();
    console_logger::debug << "generate_random_knot: update_tree completed"
                          << std::endl;
    console_logger::info << "=== generate_random_knot END ===" << std::endl;
  }

  // Update tree structure
  void update_tree() {
    printf("=== update_tree START ===\n");
    printf("update_tree: Processing %zu points\n", test_points.size());

    printf("update_tree: Calling arp::make_hash_N...\n");
    auto p_test_points = permuted(test_points, adjacency);
    auto [tree_hashes, tree_indices, internal, leaf] =
        arp::make_hash<2>(p_test_points);

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
      console_logger::debug << "=== mk_hash_tree START ===" << std::endl;

      if (test_points.empty()) {
        console_logger::debug << "mk_hash_tree: No points available for tree construction" << std::endl;
        return false;
      }

      console_logger::debug << "mk_hash_tree: Building tree for " << test_points.size() << " points" << std::endl;

      // Let's manually do what update_tree() does with more debugging
      console_logger::debug << "mk_hash_tree: Calling arp::calc_com<2>(test_points)" << std::endl;
      auto coms = arp::calc_com<2>(test_points);
      auto centers = std::vector<vec3>(coms.size());
      for (size_t i = 0; i < coms.size(); i++) {
        centers[i] = std::get<1>(coms[i]);
      }
      auto [tree_hashes, tree_indices, internal, leaf] =
          arp::make_hash_tree(centers);
      console_logger::debug << "mk_hash_tree: arp::make_hash_tree returned - hashes: " << tree_hashes.size() << ", indices: " << tree_indices.size() << ", internal: " << internal.size() << ", leaf: " << leaf.size() << std::endl;
      console_logger::debug << "mk_hash_tree: centers: " << centers.size() << std::endl;
      console_logger::debug << "mk_hash_tree: coms: " << coms.size() << std::endl;
      console_logger::debug << "mk_hash_tree: tree_hashes: " << tree_hashes.size() << std::endl;
      console_logger::debug << "mk_hash_tree: tree_indices: " << tree_indices.size() << std::endl;
      console_logger::debug << "mk_hash_tree: internal: " << internal.size() << std::endl;
      console_logger::debug << "mk_hash_tree: leaf: " << leaf.size() << std::endl;

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
      auto p_test_points = permuted_adjacency_view<2, decltype(test_points), decltype(adjacency)>(test_points, adjacency, indices);
      arp::log_hierarchy<2>(p_test_points, internal_nodes, leaf_nodes);
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
      auto p_test_points = permuted_adjacency_view<2, decltype(test_points), decltype(adjacency)>(test_points, adjacency, indices);
      arp::log_bvh<2>(p_test_points, internal_nodes, leaf_nodes);
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
        sin(C* 2.0 * t), 
        cos(C * 7.0 * t), 
        sin(C * 11.0 * t) * cos(C * 13.0 * t)
      );
      test_point = {near_test_point};
      
      // Log the test point
      geometry_logger::point(near_test_point, vec4(1.0, 1.0, 0.0, 1.0)); // Yellow test point
      const auto p_test_points = permuted_adjacency_view<2, decltype(test_points), decltype(adjacency)>(test_points, adjacency, indices);
      
      using TTYPE = decltype(p_test_points);
      using PTYPE = decltype(test_point);

      const auto bvh_result =
          arp::make_bvh<2, TTYPE>(p_test_points, internal_nodes, leaf_nodes);
      
      auto result = arp::getNearest<2, PTYPE, TTYPE>(
          test_point, p_test_points, internal_nodes, leaf_nodes,
          bvh_result, 10000.0,
          [](
            const PTYPE &t_verts, 
            const slice<2, TTYPE> &s_verts) {
            
            const vec3 &xA = t_verts[0];
            const vec3 &xB0 = s_verts[0];
            const vec3 &xB1 = s_verts[1];
            gaudi::real d = va::distance_from_line(xB0, xB1, xA);
            return d;
          });
      
      printf("log_nearest: Result: %zu\n", result.size());
      
      // Draw lines to the nearest leaf nodes
      auto masses = arp::calc_com<2>(p_test_points);
      for (size_t i = 0; i < result.size(); i++) {
        int leaf_id = result[i];
        if (leaf_id >= 0 && leaf_id < static_cast<int>(masses.size())) {
          const vec3 &leaf_point = std::get<1>(masses[leaf_id]);
          //const vec3 & leaf_point = test_points[leaf_id];
          geometry_logger::line(near_test_point, leaf_point, vec4(1.0, 0.0, 1.0, 0.8f)); // Magenta lines
        }
      }
      /*
      for(size_t i = 0; i < test_points.size(); i++){
        const vec3 pnt = test_points[i];
        geometry_logger::point(pnt, vec4(0.0, 0.5, 0.25, 1.0));
      }
        */
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

  // Get number of leaf nodes
  int get_leaf_count() const {
    if (test_points.empty()) return 0;
    console_logger::debug << "get_leaf_count: test_points.size(): " << test_points.size() << std::endl;
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
      .function("generate_trefoil_knot", &MortonTreeTest::generate_trefoil_knot)
      .function("generate_knot", &MortonTreeTest::generate_knot)
      .function("generate_random_knot", &MortonTreeTest::generate_random_knot)
      .function("mk_hash_tree", &MortonTreeTest::mk_hash_tree)
      .function("log_hierarchy", &MortonTreeTest::log_hierarchy)
      .function("log_bvh", &MortonTreeTest::log_bvh)
      .function("log_nearest", &MortonTreeTest::log_nearest)
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