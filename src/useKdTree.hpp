#ifndef USEKDTREE_H
#define USEKDTREE_H

#include <string>
#include <iostream>
#include <vector>

#include "open3d/Open3D.h"
#include "tools.hpp"
#include "kdtree.hpp"

using namespace Eigen;
using namespace std;
using namespace open3d;
using namespace Kdtree;

KdTree KdTreeFromPointCloud(std::vector<Eigen::Vector3d> points_int);

KdTree KdTreeFromMesh(geometry::TriangleMesh mesh);


#endif  // USEKDTREE_H
