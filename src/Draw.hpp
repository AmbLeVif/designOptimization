#ifndef DRAW_H
#define DRAW_H

#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <thread>
#include <memory>

#include "Eigen/Dense"
#include "CtrModel.h"
#include "loadParameters.h"
#include "open3d/Open3D.h"

#include "Target.h"
#include "useKdTree.hpp"

using namespace Eigen;
using namespace CtrLib;
using namespace std;
using namespace open3d;

std::shared_ptr<open3d::geometry::TriangleMesh> meshFromLines(open3d::geometry::LineSet Lines, double D, Eigen::Vector3d color);

geometry::LineSet* LineFromCtr(VectorXd YtotLastId, MatrixXd PosCtr, int i);

std::shared_ptr<geometry::TriangleMesh> plotTubes(VectorXd YtotLastId, MatrixXd posCtr);

std::shared_ptr<geometry::TriangleMesh> sigleTarget(Target tar, Vector3d clr, double r, bool isOriented = false);

bool isPointOnSameSide(Target target, Eigen::Vector3d point);

std::shared_ptr<geometry::LineSet> AfficheDistDura(MatrixXd posCtr, Kdtree::KdTree* KdTreeBrain, VectorXd YtotLastId, Target entreeCerveau);

std::shared_ptr< geometry::PointCloud > plotWorkspace(std::vector<MatrixXd> wrkSpace, bool colorPlot = false, double radius = 0.001);

#endif //DRAW_H