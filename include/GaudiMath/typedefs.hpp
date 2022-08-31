#ifndef __GAUDI_MATH_DEFS__
#define __GAUDI_MATH_DEFS__

#include <iostream>
#include <memory>
#include <nanogui/glutil.h>
#include <string>

namespace GaudiMath {
// nano gui deals mostly with vec3s, which is annoying
// since I'm adhering to using homogeneous coordinates
// so far as I can tell, Eigen can grab blocks, but
// that syntax is laborious, so here are a few helper functions
using Real = double;

using Mat3f = Eigen::Matrix3f;
using Mat3d = Eigen::Matrix3d;
using Mat3 = Mat3f;

using Mat4f = Eigen::Matrix4f;
using Mat4d = Eigen::Matrix4d;
using Mat4 = Mat4f;

using Vec2i = Eigen::Vector2i;
using Vec3i = Eigen::Vector3i;

using Vec3f = Eigen::Vector3f;
using Vec3d = Eigen::Vector3d;
using Vec3 = Vec3f;

using Vec4f = Eigen::Vector4f;
using Vec4d = Eigen::Vector4d;
using Vec4 = Vec4f;

using Quatf = Eigen::Quaternion<float>;
using Quatd = Eigen::Quaternion<double>;
using Quat = Quatf;

// using  MatXu = Eigen::MatrixXu;
using MatXd = Eigen::MatrixXd;
using MatXf = Eigen::MatrixXf;
using MatX = MatXd;

/*
Vec3 xyz(Vec4 in){
  return Vec3(in[0],in[1],in[2]);
};
*/

inline Vec3f xyz(Vec4f in) { return Vec3f(in[0], in[1], in[2]); };

inline Vec3d xyz(Vec4d in) { return Vec3d(in[0], in[1], in[2]); };

inline Vec4 xyzw(Vec3d in) { return Vec4(in[0], in[1], in[2], 1); };

inline Vec4 normalize(Vec4 in) {
  Vec4 out = in;
  out.normalize();
  return out;
};

inline Vec4 hnormalize(Vec4 in) {
  // this function never gets used.  As it turns out
  // most vectors that need to get normalized are already
  // nonhomogenous
  in[3] = 0;
  in.normalize();
  in[3] = 1;
  return in;
};

inline Vec4 rayPlaneIntersect(Vec4 l0, Vec4 r, Vec4 N, Vec4 p) {
  // projects a pt onto a plane with a  normal, N and point c
  // pC[3] = 0;
  //  pt[3] = 0;
  N[3] = 0;
  N.normalize();
  Vec4 dpc = p - l0;
  hnormalize(N);                          // makeSure its Normalized
  Vec4 itx = N.dot(dpc) / (N.dot(r)) * r; // project vector onto the normal
  return l0 + itx; // return the addition of pt + projected vector
}

inline Vec4 projectToPlane(Vec4 pt, Vec4 N, Vec4 c) {
  N[3] = 0;
  N.normalize();
  Vec4 dpc = c - pt;
  hnormalize(N);                // makeSure its Normalized
  Vec4 projpc = N.dot(dpc) * N; // project vector onto the normal
  return pt + projpc;           // return the addition of pt + projected vector
}

inline Vec4 projectToLine(Vec4 p, Vec4 l0, Vec4 l1) {
  // very ver similar to project to plane...

  Vec4 dl = l0 - l1;
  Vec4 dlp = l0 - p;

  Vec4 projdlpdl = dl.dot(dlp) * normalize(dl); // project vector onto the
                                                // normal
  return l0 + projdlpdl; // return the addition of pt + projected vector
}

inline Mat4 RandomRotation() {
  auto fRand = [](double fMin, double fMax) -> double {
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
  };
  Quat Q = Quat(Eigen::AngleAxisf(
      fRand(0, 2.0) * M_PI, Vec3(fRand(-1, 1), fRand(-1, 1), fRand(-1, 1))));
  Q.normalize();
  Mat3 R3 = Q.matrix();
  Mat4 R4;
  R4.setIdentity();
  R4.block<3, 3>(0, 0) = R3;
  return R4;
}

} // namespace GaudiMath
#endif
