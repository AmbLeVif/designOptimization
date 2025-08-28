#include "utils.h"

void RemplirPointCloud(open3d::geometry::PointCloud &pcd) {
    pcd.points_.emplace_back(0.0, 0.0, 0.0);
    pcd.points_.emplace_back(1.0, 0.0, 0.0);
    pcd.points_.emplace_back(0.0, 1.0, 0.0);
    pcd.points_.emplace_back(0.0, 0.0, 1.0);
}