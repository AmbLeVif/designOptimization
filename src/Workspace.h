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

#include "useKdTree.hpp"
#include "open3d/Open3D.h"

using namespace Eigen;
using namespace CtrLib;
using namespace std;
using namespace open3d;

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

//-------------------------------------------------------------------------
//-------------------------------  RRT  -----------------------------------
//-------------------------------------------------------------------------

// Structure représentant un noeud de l'arbre
struct Node {
    VectorXd actuators;  // Position des actionneurs (Eigen VectorXd pour les doubles)
    Node* parent;               // Pointeur vers le parent du noeud
    double score;               // Score du noeud
    VectorXi exploredChildren;  // Liste des configurations suivantes explorées
    bool isFullyExplored;  
    CtrModel ctr;

    Node(CtrModel ctr, Vector<double, 2*NB_TUBES> actuators = Vector<double, 2*NB_TUBES>(), Node* parent = nullptr, double score = 0.0)
        : actuators(actuators), parent(parent), score(score), exploredChildren(VectorXi::Zero(actuators.rows())), isFullyExplored(false), ctr(ctr){}

    Node(const VectorXd& actuators = VectorXd(), Node* parent = nullptr, double score = 0.0)
        : actuators(actuators), parent(parent), score(score), exploredChildren(VectorXi::Zero(actuators.rows())), isFullyExplored(false), ctr(ctr){
            ctr = parent->ctr;
            ctr.Compute(actuators, opt_LOAD);
        }
};

bool isConfigOk(CtrModel ctr, Eigen::VectorXd actuators);

Node* generateRandomConfiguration(Node* parentNode, const Eigen::VectorXd& increments, geometry::TriangleMesh* mesh, Kdtree::KdTree* KdTreeSkull);

Node* selectNodeBasedOnScore(const std::vector<Node*>& tree);

double score(Vector3d RobotTip, Kdtree::KdTree* KdTreeBrain);

Node* growRRT(std::vector<Node*>& tree, const Eigen::VectorXd& increments, geometry::TriangleMesh* mesh, Kdtree::KdTree* KdTreeSkull, Kdtree::KdTree* KdTreeBrain);

std::vector<MatrixXd> createRRT(int iterations, const Eigen::VectorXd increments, CtrModel ctr, double deltaTubes, geometry::TriangleMesh mesh, Kdtree::KdTree* KdTreeSkull, Kdtree::KdTree* KdTreeBrain);

#endif //WORKSPACE_H