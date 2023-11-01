#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "gaudi/common.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <iterator>
#include <memory>
#include <set>
#include <vector>
#include <zlib.h>

#ifndef __POINTY_WOINTY__
#define __POINTY_WOINTY__

namespace gaudi {
namespace duchamp {

using namespace asawa;

std::vector<vec3> createPoints(int N, real std = 0.5) {
  auto randNormalVec = [](real mean, real std) {
    auto randomFunc =
        [distribution_ = std::normal_distribution<double>(mean, std),
         random_engine_ = std::mt19937{std::random_device{}()}]() mutable {
          return vec3(distribution_(random_engine_),
                      distribution_(random_engine_),
                      distribution_(random_engine_));
          ;
        };
    return randomFunc;
  };

  std::vector<vec3> points;
  std::generate_n(std::back_inserter(points), N, randNormalVec(0, std));
  return points;
}

std::vector<vec3> createPoints(int N, const std::vector<vec3> &x, real std) {
  auto randNormalVec = [](vec3 xi, real mean, real std) {
    auto randomFunc =
        [xi, distribution_ = std::normal_distribution<double>(mean, std),
         random_engine_ = std::mt19937{std::random_device{}()}]() mutable {
          return xi + vec3(distribution_(random_engine_),
                           distribution_(random_engine_),
                           distribution_(random_engine_));
          ;
        };
    return randomFunc;
  };
  std::vector<vec3> rpts;
  for (int i = 0; i < x.size(); i++) {
    std::generate_n(std::back_inserter(rpts), N, randNormalVec(x[i], 0, std));
  }
  return rpts;
}

} // namespace duchamp
} // namespace gaudi
#endif