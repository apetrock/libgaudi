#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "gaudi/vec_addendum.h"

#include "GaudiGraphics/geometry_logger.h"

#include "gaudi/arp/arp.h"

#include "gaudi/bontecou/laplacian.hpp"

#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/dynamic.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/asawa/shell/shell.hpp"

#include "gaudi/asawa/primitive_objects.hpp"
#include "gaudi/common.h"

#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>
#include <zlib.h>

#ifndef __M2REFACTOR_TEST__
#define __M2REFACTOR_TEST__

namespace gaudi {
namespace duchamp {

using namespace asawa;

void debug_shell(shell::shell &M, const std::vector<vec3> verts) {
  for (int i = 0; i < M.__corners_next.size(); i += 2) {
    if (M.__corners_next[i] < 0)
      continue;
    int i0 = i;
    int i1 = M.other(i0);
    int v0 = M.vert(i0);
    int v1 = M.vert(i1);
    gg::geometry_logger::line(verts[v0], verts[v1], vec4(0.5, 0.5, 0.5, 1.0));
  }
}

void center(std::vector<vec3> &coords) {
  real accum = 0.0;
  vec3 min = coords[0];
  vec3 max = coords[0];

  for (auto &c : coords) {
    min = va::min(c, min);
    max = va::max(c, max);
  }

  vec3 dl = (max - min);
  real maxl = dl[0];
  maxl = maxl > dl[1] ? maxl : dl[1];
  maxl = maxl > dl[2] ? maxl : dl[2];
  real s = 2.0 / maxl;
  vec3 cen = (0.5 * (max + min) - min) * s;
  std::cout << " scale: " << s << std::endl;
  cen -= min;
  for (auto &c : coords) {
    c -= min;
    c = s * c;
    c -= cen;
  }
}

std::array<vec3, 2> extents(std::vector<vec3> &coords) {
  real accum = 0.0;
  vec3 min = coords[0];
  vec3 max = coords[0];

  for (auto &c : coords) {
    min = va::min(c, min);
    max = va::max(c, max);
  }
  return {min, max};
}

/*
std::vector<vec3> get_normals(shell &M, const std::vector<vec3> &x) {

  std::vector<std::vector<int>> faces;
  for (int i = 0; i < M.face_count(); i++) {
    std::vector<int> face;
    if (M.fbegin(i) < 0)
      continue;
    M.for_each_face(
        i, [&face](int ci, asawa::shell &M) { face.push_back(M.vert(ci)); });

    faces.push_back(face);
  }
}
*/
class spin_twist {
public:
  spin_twist(real r, vec3 cen, vec3 axis) : __r(r), __cen(cen), __axis(axis) {}
  vec3 op(const vec3 &p) {
    vec3 dp = p - __cen;
    vec3 axis = __axis.normalized();
    vec3 pp = va::orthogonal_project(axis, dp);

    vec3 ppN = pp.normalized();
    vec3 pC = __r * ppN;
    vec3 dpC = dp - pC;
    vec3 dpCN = dpC.normalized();

    // gg::geometry_logger::line(__cen, __cen + pC, vec4(1.0, 0.6, 0.5, 1.0));
    // gg::geometry_logger::line(p, p - dpC, vec4(0.8, 1.0, 0.35, 1.0));

    vec3 ppT = ppN.cross(axis);
    vec3 ppB = ppT.cross(dpCN);

    real tr = pp.norm();
    real br = dpC.norm();

    return tr * ppT + 2.0 * br * ppB;
  }

  real get_weight(const vec3 &p) {
    vec3 dx = op(p);
    vec2 rt = get_rthet(p, __axis, __axisX, __axisY);
    vec2 rtcen = 0.5 * (__rtmax + __rtmin);
    real dt = rtcen[1] - rt[1];
    return dt;
  }

  vec2 get_rthet(const vec3 &x, vec3 ax0, vec3 axX, vec3 axY) {

    vec3 dp = x - __cen;
    vec3 pp = va::orthogonal_project(ax0, dp);
    real ri = pp.norm();

    real px = axX.dot(pp);
    real py = axY.dot(pp);
    real ti = atan2(py, px);
    return vec2(ri, ti);
  }

  void get_bounds(const std::vector<vec3> &x) {

    auto polar = [this](const real &r, const real &t) {
      return r * (cos(t) * __axisX + sin(t) * __axisY);
    };

    vec3 axis = __axis.normalized();
    vec3 dc = vec3(0, 0, 0);

    for (int i = 0; i < x.size(); i++) {
      vec3 p0 = x[i];
      vec3 dp = p0 - __cen;
      vec3 pp = va::orthogonal_project(axis, dp);
      dc += pp;
    }

    dc /= real(x.size());
    __axisX = dc.normalized();
    __axisY = axis.cross(__axisX);

    vec2 rtcen(0, 0);
    for (int i = 0; i < x.size(); i++) {
      vec3 p0 = x[i];
      vec2 rti = get_rthet(p0, __axis, __axisX, __axisY);
      rtcen += rti;
    }
    rtcen /= real(x.size());
    vec3 prt = polar(rtcen[0], rtcen[1]);
    __axisX = prt.normalized();
    __axisY = axis.cross(__axisX);

    gg::geometry_logger::line(__cen, __cen + __axisX, vec4(1.0, 0.0, 0.2, 1.0));
    gg::geometry_logger::line(__cen, __cen + __axisY, vec4(0.0, 1.0, 0.2, 1.0));

    vec2 rtmin = rtcen;
    vec2 rtmax = rtcen;
    for (int i = 0; i < x.size(); i++) {
      vec3 p0 = x[i];
      vec2 rti = get_rthet(p0, __axis, __axisX, __axisY);

      rtmin = va::min(rtmin, rti);
      rtmax = va::max(rtmax, rti);
    }

    __rtmin = rtmin, __rtmax = rtmax;

    // vec3 prt = polar(rtcen[0], rtcen[1]);
    vec3 prt00 = polar(rtmin[0], rtmin[1]);
    vec3 prt01 = polar(rtmin[0], rtmax[1]);
    vec3 prt10 = polar(rtmax[0], rtmin[1]);
    vec3 prt11 = polar(rtmax[0], rtmax[1]);

    gg::geometry_logger::line(__cen, __cen + prt, vec4(1.0, 0.6, 0.5, 1.0));

    vec4 col = vec4(0.0, 0.6, 1.0, 1.0);
    gg::geometry_logger::line(__cen + prt00, __cen + prt10, col);
    gg::geometry_logger::line(__cen + prt01, __cen + prt11, col);

    int N = 128;

    real dt = (rtmax[1] - rtmin[1]) / real(N);
    for (int i = 0; i < N; i++) {
      real t0 = __rtmin[1] + real(i) * dt;
      real t1 = t0 + dt;

      vec3 prt00 = polar(rtmin[0], t0);
      vec3 prt01 = polar(rtmin[0], t1);
      vec3 prt10 = polar(rtmax[0], t0);
      vec3 prt11 = polar(rtmax[0], t1);

      gg::geometry_logger::line(__cen + prt00, __cen + prt01, col);
      gg::geometry_logger::line(__cen + prt10, __cen + prt11, col);
    }
  }

  void debug() {
    gg::geometry_logger::line(__cen - __axis, __cen + __axis,
                              vec4(1.0, 0.6, 0.5, 1.0));
  }

  vec2 __rtmin, __rtmax;
  real __r;
  vec3 __cen;
  vec3 __axis, __axisX, __axisY;
};

void test() {
  shell::shell::ptr M = shell::load_cube();
  shell::triangulate(*M);

  vec3_datum::ptr v_datum = static_pointer_cast<vec3_datum>(M->get_datum(0));
  const std::vector<vec3> &coords = v_datum->data();
  real l0 = 2.0 * shell::avg_length(*M, v_datum->data());

  shell::dynamic::ptr surf = shell::dynamic::create(M, l0, 3.0 * l0, 0.25 * l0);

  remove_vertex(*M, 0);
  remove_vertex(*M, 1);
  remove_vertex(*M, 2);
  remove_vertex(*M, 3);
  remove_vertex(*M, 4);

  //  remove_vertex(*M, 2);
  //   subdivide_edges(*M);
  //   collapse_edges(*M);

  //  gather_edges_parallel(*M, vdata->data(), 1.0);
  //   merge_face(*M, 3, M->other(3));
  debug_shell(*M, v_datum->data());
}

class dynamic_shell_test {
public:
  typedef std::shared_ptr<dynamic_shell_test> ptr;

  static ptr create() { return std::make_shared<dynamic_shell_test>(); }

  dynamic_shell_test() {
    //__M = load_cube();
    __M = shell::load_bunny();

    shell::triangulate(*__M);
    for (int i = 0; i < __M->face_count(); i++) {
      if (__M->fbegin(i) > 0) {
        assert(__M->fsize(i) == 3);
      }
    }

    vec3_datum::ptr c_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(0));

    std::vector<vec3> &coords = c_datum->data();
    center(coords);

    /////////
    // twist
    /////////

    std::array<vec3, 2> ext = extents(coords);
    vec3 cen = 0.5 * (ext[1] + ext[0]);
    vec3 de = (ext[1] - ext[0]);
    vec3 twist_axis = vec3(0, 0, 1);
    vec3 twist_cen = cen + vec3(0.35 * de[0], 0.35 * de[1], 0.0);
    real r = 0.75 * de[0];
    _twist = std::make_shared<spin_twist>(r, twist_cen, twist_axis);

    /////////
    // dynamic surface
    /////////

    real l0 = 1.0 * asawa::shell::avg_length(*__M, c_datum->data());
    __surf = shell::dynamic::create(__M, 1.0 * l0, 3.0 * l0, 1.0 * l0);

    /////////
    // weights
    /////////

    vec3_datum::ptr x_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(0));
    std::vector<vec3> &x = x_datum->data();

    std::vector<real> weights(__M->vert_count(), 0.0);
    _twist->get_bounds(x);
    for (int i = 0; i < x.size(); i++) {
      weights[i] = _twist->get_weight(x[i]);
    }

    datum_t<real>::ptr wdata =
        datum_t<real>::create(prim_type::VERTEX, weights);
    _iw = __M->insert_datum(wdata);
  };

  void test_twist() {
    vec3_datum::ptr x_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(0));
    std::vector<vec3> &x = x_datum->data();

    vec3_datum::ptr v_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(1));
    std::vector<vec3> &dx = v_datum->data();

    real_datum::ptr w_datum =
        static_pointer_cast<real_datum>(__M->get_datum(_iw));
    std::vector<real> &w = w_datum->data();

    _twist->get_bounds(x);
#if 1
    real dt = 0.01;
    std::cout << "vert_count: " << x.size() << std::endl;
    for (int i = 0; i < x.size(); i++) {
      vec3 p0 = x[i];
      vec3 v0 = (w[i] - 1.0) * _twist->op(p0);
      vec3 pp = p0 + 0.5 * dt * v0;
      vec3 v1 = (w[i] - 1.0) * _twist->op(pp);

      dx[i] = dt * v1;
    }
  }

  void step(int frame) {

    std::cout << "frame: " << frame << std::endl;
    test_twist();

    vec3_datum::ptr v_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(1));
    std::vector<vec3> &dx = v_datum->data();

    __surf->step(1.0, dx);

    _twist->debug();
#endif

    // debug_shell(*__M, x_datum->data());
  }

  index_t _iw;
  shell::shell::ptr __M;
  shell::dynamic::ptr __surf;
  std::shared_ptr<spin_twist> _twist;
};

} // namespace duchamp
} // namespace gaudi
#endif