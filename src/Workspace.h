#ifndef WORKSPACE_H
#define WORKSPACE_H

#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <thread>

#include "Eigen/Dense"
#include "CtrModel.h"
#include "loadParameters.h"

using namespace Eigen;
using namespace CtrLib;
using namespace std;

struct WorkspaceData {

    std::vector<VectorXd> q;
    std::vector<MatrixXd> data;
    std::vector<Vector3d> YtotLastIdList;
    std::vector<Vector3d> endPos;
    std::vector<MatrixXd> orientations;
    std::vector<int> isStable;

    WorkspaceData()
        : q(), data(), YtotLastIdList(), endPos(), orientations(), isStable(){}
};

std::vector<Eigen::MatrixXd> genererRevolution(
    const std::vector<Eigen::MatrixXd>& pointsInitiaux,  // Points initiaux
    double discretAngle                                  // Discrétisation de la révolution (en radians)
    );

WorkspaceData buildWorkspace(CtrModel ctr, int discretAngle, int discretLength, VectorXd offset, bool isOriented, bool smart, int freqRevolution);


#endif //WORKSPACE_H