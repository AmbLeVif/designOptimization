#ifndef TOOLS_H
#define TOOLS_H

#include <vector>
#include <cmath>
#include <iostream>

#include "Eigen/Dense"


double degToRad(double degrees);
Eigen::Matrix4d eulerToHomogeneousMatrix(const Eigen::Vector3d& eulerAnglesDeg, const Eigen::Vector3d& translation);
std::vector<double> eigenToStdVector(const Eigen::VectorXd& eigenVec);
Eigen::VectorXd stdVectorToEigen(const std::vector<double>& stdVec);
double angleBetweenDirections(const Eigen::Matrix3d& m1, const Eigen::Matrix3d& m2);
Eigen::Vector3d getRainbow(double value);

#endif  // TOOLS_H
