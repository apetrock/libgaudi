/*
 *  manifold_singleton.cpp
 *  Phase Vocoder
 *
 *  Created by John Delaney on 12/29/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

//#include <gl/glut.h>
#include "GaudiGraphics/geometry_logger.h"
#include "GaudiMath/typedefs.hpp"
#include <numeric>

namespace gg {

//////////////////////////
// debugger
//////////////////////////

bool geometry_logger::instance_flag = false;
geometry_logger *geometry_logger::global_instance = NULL;

geometry_logger &geometry_logger::get_instance() {
  static geometry_logger logger;
  if (!logger.initialized()) {
    logger.initialized() = true;
  }

  return logger;
}

void geometry_logger::render() {
  geometry_logger &logger = geometry_logger::get_instance();
  logger.debugLines->renderLines();
}

void geometry_logger::clear() {
  geometry_logger &logger = geometry_logger::get_instance();
  logger.debugLines->clear();
}

void geometry_logger::line4(const Vec4 &p0, const Vec4 &p1, const Vec4 &color) {
  geometry_logger &logger = geometry_logger::get_instance();
  logger.debugLines->pushLine(p0, p1, color);
}
void geometry_logger::line(const Vec3d &p0, const Vec3d &p1,
                           const Vec4d &color) {
  geometry_logger &logger = geometry_logger::get_instance();
  logger.debugLines->pushLine(xyzw(p0), xyzw(p1),
                              Vec4(color[0], color[1], color[2], color[3]));
}

void geometry_logger::box(const Vec3d &cen, const Vec3d &h,
                          const Vec4d &color) {
  geometry_logger &logger = geometry_logger::get_instance();
  logger.debugLines->pushBox(xyzw(cen), xyzw(h),
                             Vec4(color[0], color[1], color[2], color[3]));
}

void geometry_logger::ext(const Vec3d &mn, const Vec3d &mx,
                          const Vec4d &color) {
  geometry_logger &logger = geometry_logger::get_instance();

  Vec3d cen = 0.5 * (mx + mn);
  Vec3d h = 0.5 * (mx - mn);

  logger.debugLines->pushBox(xyzw(cen), xyzw(h),
                             Vec4(color[0], color[1], color[2], color[3]));
}

Vec3d _rainbow(double d) {
  double r = 0.5 + 0.5 * cos(2.0 * M_PI * (d + 0.000));
  double g = 0.5 + 0.5 * cos(2.0 * M_PI * (d + 0.333));
  double b = 0.5 + 0.5 * cos(2.0 * M_PI * (d + 0.666));
  return Vec3d(r, g, b);
}

Vec3d _grey(double d) { return Vec3d(0.5, 0.5, 0.5); }
Vec3d _red(double d) { return Vec3d(1.0, 0.0, 0.0); }
Vec3d _green(double d) { return Vec3d(0.0, 1.0, 0.0); }
Vec3d _blue(double d) { return Vec3d(0.0, 0.0, 1.0); }

void geometry_logger::lines(const std::vector<Vec3d> &p0s,
                            const std::vector<Vec3d> &p1s,
                            const std::vector<Vec4d> &colors) {

  for (int i = 0; i < p0s.size(); i++) {

    const auto &p0 = p0s[i];
    const auto &p1 = p1s[i];
    const auto &c = colors[i];
    geometry_logger &logger = geometry_logger::get_instance();
    logger.debugLines->pushLine(xyzw(p0), xyzw(p1),
                                Vec4(c[0], c[1], c[2], c[3]));
  }
}

void geometry_logger::field(const std::vector<Vec3d> &p,
                            const std::vector<Vec3d> &dirs, double D,
                            PresetColor col) {

  double mx = std::accumulate(p.begin(), p.end(), 0.0, [](double a, auto &c) {
    return std::max(a, c.norm());
  });

  std::vector<Vec3d> p0(p.size(), Vec3d::Zero());
  std::vector<Vec3d> p1(p.size(), Vec3d::Zero());

  std::vector<Vec4d> colors(p.size(), Vec4d::Zero());

  for (int i = 0; i < p0.size(); i++) {

    const auto &dir = dirs[i];

    // std::cout << a.transpose() << std::endl;
    auto mag = dir.norm();
    auto c = _grey(mag);
    if (col == red)
      c = _red(mag);
    if (col == green)
      c = _green(mag);
    if (col == blue)
      c = _blue(mag);
    if (col == rainbow)
      c = _rainbow(mag);

    p1[i] = p[i] + D * dir / mx;
    colors[i] = Vec4d(c[0], c[1], c[2], mag / mx);
  }

  gg::geometry_logger::lines(p, p1, colors);
}

void geometry_logger::frame(Mat3d M, Vec3d c, double C = 1.0) {

  Vec3d t0 = M.block(0, 0, 3, 1);
  Vec3d t1 = M.block(0, 1, 3, 1);
  Vec3d t2 = M.block(0, 2, 3, 1);

  gg::geometry_logger::line(c - C * t0, c + C * t0, Vec4d(1.0, 0.0, 0.0, 1.0));
  gg::geometry_logger::line(c - C * t1, c + C * t1, Vec4d(0.0, 1.0, 0.0, 1.0));
  gg::geometry_logger::line(c - C * t2, c + C * t2, Vec4d(0.0, 0.0, 1.0, 1.0));
}

} // namespace gg
