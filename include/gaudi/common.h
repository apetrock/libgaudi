

#ifndef __LIBGAUDI_COMMON_TYPEDEFS__
#define __LIBGAUDI_COMMON_TYPEDEFS__

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
namespace gaudi {
typedef int index_t;
typedef double real;
typedef Eigen::Matrix<real, 2, 1> vec2;
typedef Eigen::Matrix<real, 3, 1> vec3;
typedef Eigen::Matrix<real, 4, 1> vec4;

typedef Eigen::Matrix<real, 3, 3> mat3;
typedef Eigen::Matrix<real, 4, 4> mat4;

} // namespace gaudi
#endif