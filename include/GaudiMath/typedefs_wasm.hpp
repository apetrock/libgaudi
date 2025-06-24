#pragma once

#include <Eigen/Dense>

namespace GaudiMath {

typedef Eigen::Matrix<double, 3, 1, 0, 3, 1> Vector3d;
typedef Eigen::Matrix<double, 4, 1, 0, 4, 1> Vector4d;
typedef Eigen::Matrix<double, -1, 1, 0, -1, 1> VectorXd;
typedef Eigen::Matrix<double, 3, 3, 0, 3, 3> Matrix3d;
typedef Eigen::Matrix<double, 4, 4, 0, 4, 4> Matrix4d;
typedef Eigen::Matrix<double, -1, -1, 0, -1, -1> MatrixXd;

}
