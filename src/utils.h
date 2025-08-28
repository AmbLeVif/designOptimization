#ifndef UTILS_H
#define UTILS_H

#include <open3d/Open3D.h>
#include "CtrModel.h"

#include "kdtree.hpp"
#include "tools.hpp"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Side_of_triangle_mesh.h>

void RemplirPointCloud(open3d::geometry::PointCloud &pcd);

#endif
