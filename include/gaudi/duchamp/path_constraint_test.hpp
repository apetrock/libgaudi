#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include "gaudi/common.h"
#include "gaudi/geometry_types.hpp"
#include "gaudi/logger.hpp"
#include "gaudi/vec_addendum.h"
#include "gaudi/duchamp/forward_declarations.hpp"

#include "gaudi/asawa/rod/dynamic.hpp"
#include "gaudi/asawa/rod/rod.hpp"

#include "gaudi/calder/tangent_point_integrators.hpp"

// Only include what we actually need for function calls
#include "gaudi/hepworth/block/generic_constraints.hpp"
#include "gaudi/hepworth/block/rod_constraints.hpp"
#include "gaudi/hepworth/block/rod_constraints_init.hpp"
#include "gaudi/hepworth/block/sim_block.hpp"
#include "gaudi/hepworth/block/solver.hpp"

#include "utils/sdf.hpp"
#include "gaudi/geometry_logger.hpp"

#ifndef __M2REFACTOR_TEST__
#define __M2REFACTOR_TEST__

namespace gaudi {
namespace duchamp {

using namespace asawa;

class block_test {
public:
  typedef std::shared_ptr<block_test> ptr;

  static ptr create() { return std::make_shared<block_test>(); }

  block_test() {
    // load_fib_rod();
    load_strands();
    int N = 13;
    real r1 = 1.5;
    real r11 = 0.5;
    real pi43 = 4.0 / 3.0 * M_PI;
    real v0 = real(13.0) * pi43 * pow(r11, 3.0);
    real r0 = std::pow(v0 / pi43, 1.0 / 3.0);
    std::cout << "radius 0" << r0 << std::endl;
    __sdf0 = sdf_sphere::create(vec3(0.0, 0.0, 0.0), r0);
    __sdf1 = sdf_multi_sphere::create(get_fib(r1, 13), r11);
    // load_sdf();

    // Log initial setup
    geometry_logger::clear();
    geometry_logger::line(vec3(0, 0, 0), vec3(2, 0, 0), vec4(1, 0, 0, 1)); // X axis
    geometry_logger::line(vec3(0, 0, 0), vec3(0, 2, 0), vec4(0, 1, 0, 1)); // Y axis
    geometry_logger::line(vec3(0, 0, 0), vec3(0, 0, 2), vec4(0, 0, 1, 1)); // Z axis
  };

  std::vector<vec3> get_fib(real r0, int N = 13) {

    real golden = 0.5 * (1.0 + sqrt(5));
    std::vector<vec3> cens(N, vec3::Zero());
    for (int i = 0; i < N; i++) {
      real theta = 2.0 * M_PI * i / golden;
      real phi = acos(1.0 - 2.0 * (i + 0.5) / real(N));
      vec3 ci =
          r0 * vec3(cos(theta) * sin(phi), sin(theta) * sin(phi), cos(phi));
      cens[i] = ci;
    }
    return cens;
  }

  void load_fib_rod() {
    __R = rod::rod::create(get_fib(1.5, 13));

    real lavg = 0.01 * __R->lavg();
    __Rd = rod::dynamic::create(__R, 0.25 * lavg, 2.5 * lavg, 0.25 * lavg);
  }

  void load_strands() {
    // Create a blank rod
    __R = rod::rod::create();
    
    // Set rod radius for spacing calculations
    real rod_radius = 0.05;  // Default rod radius
    __R->set_radius(rod_radius);
    
    // Create 16 strands: 4x4 grid at z=0 to circle at z=10
    for (int i = 0; i < 16; i++) {
      load_strand(*__R, i, rod_radius);
    }
    
    real lavg = __R->lavg();
    __Rd = rod::dynamic::create(__R, 0.25 * lavg, 2.5 * lavg, 0.25 * lavg);
  }

  void load_strand(rod::rod &R, int strand_index, real rod_radius) {
    // Calculate grid position (4x4 grid)
    int grid_row = strand_index / 4;  // 0-3
    int grid_col = strand_index % 4;  // 0-3
    
    // Grid spacing based on rod radius (4 rod diameters between centers)
    real grid_spacing = 8.0 * rod_radius;
    
    // Grid center at [x0, y0] = [-2.0, 0.0]
    vec3 grid_center(-2.0, 0.0, 0.0);
    
    // Start point: 4x4 grid at z=0
    vec3 start_point(
      grid_center.x() + (grid_col - 1.5) * grid_spacing,  // x: grid center + offset
      grid_center.y() + (grid_row - 1.5) * grid_spacing,  // y: grid center + offset
      0.0  // z=0
    );

    // Circle center at [x1, y1] = [2.0, 0.0]
    vec3 circle_center(2.0, 0.0, 0.5);
    real circle_radius = 2.0 * grid_spacing;  // Circle radius based on grid spacing
    real angle = 2.0 * M_PI * strand_index / 16.0;  // Distribute evenly around circle
    
    // End point: circle at z=0.5
    vec3 end_point(
      circle_center.x() + circle_radius * cos(angle),  // x: circle center + radius * cos
      circle_center.y() + circle_radius * sin(angle),  // y: circle center + radius * sin
      circle_center.z()  // z=0.5
    );

    std::cout << "Strand " << strand_index << " start: " << start_point << std::endl;
    std::cout << "Strand " << strand_index << " end: " << end_point << std::endl;

    // Create a straight strand between the two points
    int num_segments = 32; // Number of segments in the strand
    std::vector<vec3> strand_points;
    strand_points.reserve(num_segments + 1);

    for (int i = 0; i <= num_segments; i++) {
      real t = real(i) / real(num_segments);
      vec3 point = va::mix(t, start_point, end_point);
      strand_points.push_back(point);
    }

    // Compute total initial length
    _lt0 = 0.0;
    for (int i = 0; i < num_segments; i++) {
      _lt0 += (strand_points[i + 1] - strand_points[i]).norm();
    }

    // Insert the strand into the rod
    R.insert_strand(strand_points, false); // false = not a loop
  }

#if 1
  std::vector<vec3> compute_tangent_point_gradient() {
    real eps = __Rd->_Cc;
    std::vector<vec3> &x = __R->x();
    std::vector<real> l = __R->l0();
    std::vector<vec3> T = __R->N2c();
    std::vector<vec3> xc = __R->xc();

    std::vector<vec3> g0 =
        calder::tangent_point_gradient(*__R, x, l, T, 1.0 * eps, 6.0);

    // Log tangent point gradients
#if 0
    for (size_t i = 0; i < g0.size(); i++) {
      if (g0[i].norm() > 1e-6) {
        geometry_logger::line(xc[i], xc[i] + 0.1 * g0[i], vec4(1, 0.5, 0, 0.8));
      }
    }
#endif

    return g0;
  }
#endif

  sdf_base::ptr get_sdf(int frame) {
    if ((frame / 400) % 2 == 0) {
      return __sdf0;
    } else {
      return __sdf1;
    }
  }

  std::vector<real> compute_growth_weights(index_t frame) {

    auto sdf = get_sdf(frame);

    std::vector<real> dists = sdf->distance(__R->__x);
    std::vector<real> w(__R->__x.size(), 0);
    std::vector<vec3> xc = __R->xc();
    std::vector<vec3> N = __R->N0c();
    auto verts = __R->get_vert_range();
    for (auto i : verts) {
      asawa::rod::consec_t c = __R->consec(i);

      vec3 xi = xc[i];
      vec3 Ni = N[i];
      real di = dists[i];

#if 0
          vec4 c0 = vec4(0.0, 1.0, 0.0, 1.0);
        vec4 c1 = vec4(1.0, 0.0, 0.0, 1.0);
        if(di > 0.0){
            geometry_logger::line(xi, xi + 0.1*di * Ni, c0);
        } else{
            geometry_logger::line(xi, xi + 0.1*di * Ni, c1);
        }
#endif
      di = di < 0.0 ? -1.0 : 0.5 * di;
      w[i] = (1.0 - 0.28 * di);
      w[i] = va::clamp(w[i], 0.0, 4.00);
      // l0[i] = std::max(l0[i], 1e-8);
    }

    return std::move(w);
  }

  std::vector<vec3> compute_boundary_gradients(index_t frame) {

    auto sdf = get_sdf(frame);

    std::vector<real> dists = sdf->distance(__R->__x);
    std::vector<vec3> gdists = sdf->grad_distance(__R->__x);
    std::vector<vec3> f(__R->__x.size(), vec3::Zero());
    std::vector<vec3> xc = __R->xc();
    for (int i = 0; i < __R->__x.size(); i++) {
      if (dists[i] > 0.0) {
        f[i] = -dists[i] * gdists[i];
      }
#if 0
        vec4 c0 = vec4(0.0, 1.0, 0.0, 1.0);
        vec4 c1 = vec4(1.0, 0.0, 0.0, 1.0);
        if(dists[i] > 0.0){
            geometry_logger::line(xc[i], xc[i] + f[i], c0);
        } else{
            geometry_logger::line(xc[i], xc[i] + f[i], c1);   
        }
#endif
    }

    // Log boundary gradients
    for (size_t i = 0; i < f.size(); i++) {
      if (f[i].norm() > 1e-6) {
        vec4 color = dists[i] > 0.0 ? vec4(0, 1, 0, 0.8) : vec4(1, 0, 0, 0.8);
        geometry_logger::line(xc[i], xc[i] + 0.1 * f[i], color);
      }
    }

    return std::move(f);
  }


  void step_dynamics(int frame) {
    std::cout << "frame: " << frame << ", size: " << __R->__x.size()
              << std::endl;
    hepworth::block::projection_solver solver;

    std::vector<hepworth::projection_constraint::ptr> constraints;

    std::vector<real> &l0 = __R->__l0;

    real h = 0.05;
    std::vector<vec3> f(__R->__v.size(), vec3::Zero()); 
    // std::vector<vec3> fr = compute_coulomb_gradient();
    // std::vector<vec3> fr = compute_null_coulomb_gradient();
    //std::vector<vec3> fr = compute_tangent_point_gradient();
    //f = 1e-6 * fr;  // Add rotation torque to the forces
    hepworth::vec3_block::ptr x =
        hepworth::vec3_block::create(__R->__M, __R->__x, __R->__v, f);
    hepworth::quat_block::ptr u =
        hepworth::quat_block::create(__R->__J, __R->__u, __R->__o);
    pin_endpoints(constraints, {x});  // Commented out to allow free rotation
    // hepworth::rod::init_smooth(*__R, constraints, 0.2);

    //  hepworth::rod::init_smooth_bend(*__R, constraints, 0.01);


    // Use dual-weight stretch_shear: w1 for stretch/shear, w2 for rotation
    hepworth::block::init_stretch_shear(*__R, constraints, l0, 1e-2, 1e-6, {x, u});
    hepworth::block::init_bend_twist(*__R, constraints, 1e-3, {u}, false);

    hepworth::block::init_collisions(*__R, *__Rd, constraints, 1.0, {x, x});
    solver.set_constraints(constraints);

    // f[0][0] = 1.0;
    std::vector<hepworth::sim_block::ptr> blocks = {x, u};
    solver.step(blocks, h, 0.5);

    // Log constraint forces
    // log_constraint_forces(x, u);
  }

  void step(int frame) {

    _frame = frame;

    // Clear previous frame's debug lines
    geometry_logger::clear();

    // Log coordinate axes
    // geometry_logger::line(vec3(0, 0, 0), vec3(2, 0, 0), vec4(1, 0, 0, 1)); // X axis
    // geometry_logger::line(vec3(0, 0, 0), vec3(0, 2, 0), vec4(0, 1, 0, 1)); // Y axis
    // geometry_logger::line(vec3(0, 0, 0), vec3(0, 0, 2), vec4(0, 0, 1, 1)); // Z axis

    step_dynamics(frame);
    __Rd->step();

    // Log final rod geometry
    log_rod_geometry(vec4(0.8, 0.8, 0.8, 1.0));

    // Log pinned points (endpoints)
    log_pinned_points();
    
    // Log rotation-affected endpoints
    log_rotation_affected_points();

    if (frame > 3000)
      exit(0);
    //__R->debug();
  }
  // Log pinned points for visualization

  void
  pin_endpoints(std::vector<hepworth::projection_constraint::ptr> &constraints,
                std::vector<hepworth::sim_block::ptr> blocks) {
    if (!__R)
      return;
    const std::vector<vec3> &x = __R->x();
    std::vector<index_t> pinned = get_pinned_indices();
    for (int iv = 0; iv < pinned.size(); iv++) {
      constraints.push_back(hepworth::block::pinned::create(
          std::vector<index_t>({pinned[iv]}), x[pinned[iv]], 1.0, blocks));
    }
  }

  // Helper function to log rod geometry
  void log_rod_geometry(const vec4 &color) {
    if (!__R)
      return;

    // Use the existing debug() method which properly handles the adjacency
    // table
    __R->debug();
  }

  // Helper function to log constraint forces
  void log_constraint_forces(hepworth::vec3_block::ptr x,
                             hepworth::quat_block::ptr u) {
    if (!x || !__R)
      return;

    const std::vector<vec3> &positions = __R->x();
    const std::vector<vec3> &forces = x->_f;

    // Log forces for valid vertices
    for (int i = 0; i < __R->corner_count(); i++) {
      if (__R->next(asawa::rod::corner_id(i)) < 0)
        continue;

      if (i < forces.size() && forces[i].norm() > 1e-6) {
        geometry_logger::line(positions[i], positions[i] + 0.05 * forces[i],
                     vec4(1, 1, 0, 0.8));
      }
    }
  }

  // Get the pinned indices (endpoints of the strand)
  std::vector<index_t> get_pinned_indices() const {
    std::vector<index_t> pinned;
    if (!__R)
      return pinned;

    real z_threshold = 0.1;
    const std::vector<vec3>& positions = __R->x();
    
    // Use the new get_endpoints() convenience function
    std::vector<index_t> endpoints = __R->get_endpoints();
    
    // Only pin endpoints below the z-threshold
    for (index_t i : endpoints) {
      if (i < positions.size() && positions[i].z() < z_threshold) {
        pinned.push_back(i);
      }
    }
    return pinned;
  }

  // Log pinned points for visualization (now commented out since we're not pinning)
  void log_pinned_points() {
    if (!__R)
      return;

    // std::vector<index_t> pinned = get_pinned_indices();
    // const std::vector<vec3> &positions = __R->x();

    // for (index_t idx : pinned) {
    //   if (idx < positions.size()) {
    //     // Log pinned points as red spheres
    //     geometry_logger::point(positions[idx], vec4(1.0, 0.0, 0.0, 1.0));
    //   }
    // }
  }

  
  // Log rotation-affected endpoints for visualization
  void log_rotation_affected_points() {
    if (!__R) return;
    
    const std::vector<vec3>& positions = __R->x();
    real z_threshold_high = 0.1;
    real z_threshold_low = -0.1;
    
    // Use the new get_endpoints() convenience function
    std::vector<index_t> endpoints = __R->get_endpoints();
    
    for (index_t i : endpoints) {
      if (i < positions.size()) {
        vec3 pos = positions[i];
        
        // Mark endpoints above high threshold (clockwise rotation) as green
        if (pos.z() >= z_threshold_high) {
          geometry_logger::point(pos, vec4(0.0, 1.0, 0.0, 1.0));
        }
        // Mark endpoints below low threshold (counter-clockwise rotation) as red
        else if (pos.z() <= z_threshold_low) {
          geometry_logger::point(pos, vec4(1.0, 0.0, 0.0, 1.0));
        }
      }
    }
  }

  int _frame;
  real _lt0 = 1.0;
  sdf_base::ptr __sdf0;
  sdf_base::ptr __sdf1;

  rod::rod::ptr __R;
  rod::dynamic::ptr __Rd;
};

} // namespace duchamp
} // namespace gaudi
#endif