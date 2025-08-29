#ifndef PSO_H
#define PSO_H

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

// Déclaration de la structure Particle (particule)
struct Particle {
    vector<double> position;
    vector<double> velocity;
    vector<double> best_position; 
    vector<double> best_objective_values;
    vector<double> objective_values;
};

struct PsoParam {
    int nbIteration;
    int nbParticules;
    int nbDim;
    std::vector<double> lowLimit;
    std::vector<double> highLimit;

    int nbObjFunction;
    int closeTargets;
    double deltaTarget;
    double deltaDura;

    int discretLength;
    int discretAngle;

    int nbObjective_values;
    double goodDistToTarget;
    double goodDistToDura;
    std::vector<double> goodDist;
    PsoParam() 
        : nbIteration(1), nbParticules(1), 
        nbDim(1), 
        lowLimit(std::vector<double>(0)), highLimit(std::vector<double>(0)), 
        nbObjFunction(1), closeTargets(1), deltaTarget(0.0), deltaDura(0.0), 
        discretLength(0), discretAngle(0),
        nbObjective_values(1), goodDistToTarget(0), goodDistToDura(0), goodDist(std::vector<double>(0)) {}

    PsoParam(int nbIteration, int nbParticules, int nbDim, std::vector<double> low, std::vector<double> high, int closeTargets, double deltaTarget, double deltaDura) 
        : nbIteration(nbIteration), nbParticules(nbParticules), nbDim(nbDim), closeTargets(closeTargets), deltaTarget(deltaTarget), deltaDura(deltaDura) {
            lowLimit = low;
            highLimit = high;
        }
};

VectorXd ObjectiveFunctionTargetPosition(std::vector<MatrixXd> wrkSpace, Vector3d target, VectorXi* NclosestId);

VectorXd ObjectiveFunctionObstacle(std::vector<MatrixXd> wrkSpace, VectorXi* NclosestId, Kdtree::KdTree* KdTreeBrain, Target entreeCerveau);

VectorXd ObjectiveFunction(PsoParam PsoParam, std::vector<MatrixXd> wrkSpace, Target target, Kdtree::KdTree* KdTreeBrain, Target entreeCerveau);

std::shared_ptr< geometry::TriangleMesh > PlotValidRobot(PsoParam PSOparam, std::vector<MatrixXd> wrkSpace, std::vector<Vector3d> YtotLastId, std::vector<Target> targets, Kdtree::KdTree* KdTreeBrain, Target entreeCerveau, std::shared_ptr<geometry::TriangleMesh> mesh, Kdtree::KdTree* KdTreeSkull);

vector<double> setObjectiveValue(PsoParam PsoParam, vector<double> ctrParam, std::vector<Target> targets, Kdtree::KdTree* KdTreeBrain, Target entreeCerveau,  std::shared_ptr<geometry::TriangleMesh> mesh, Kdtree::KdTree* KdTreeSkull);

void initializeParticles(vector<Particle>& particles, PsoParam param, 
                         const vector<double>& vitesse_max, std::vector<Target> targets, Kdtree::KdTree* KdTreeBrain,
                         Target entreeCerveau, std::shared_ptr<geometry::TriangleMesh> mesh, Kdtree::KdTree* KdTreeSkull);

bool dominates(const std::vector<double>& a, const std::vector<double>& b);

std::vector<Particle> getParetoFront(std::vector<Particle> oldPareto, std::vector<Particle> particles);

void updateParticle(PsoParam param, Particle& p, const vector<double>& global_best_position, 
                    const vector<double>& vitesse_max, std::vector<Target> targets, 
                    Kdtree::KdTree* KdTreeBrain, Target entreeCerveau, 
                    std::shared_ptr<geometry::TriangleMesh> mesh, Kdtree::KdTree* KdTreeSkull);

bool isGoodParticle(PsoParam param, int nbTargets, Particle p);

void updateParticlesWithThreads(std::vector<Particle>& particles,
                                 std::vector<Particle>& result,
                                 const PsoParam& param,
                                 const std::vector<double>& global_best_position,
                                 const std::vector<double>& vitesse_max,
                                 const std::vector<Target>& targets,
                                 Kdtree::KdTree* KdTreeBrain,
                                 Target entreeCerveau,
                                 std::shared_ptr<geometry::TriangleMesh> mesh, 
                                 Kdtree::KdTree* KdTreeSkull);

vector<Particle> PSO(PsoParam param, std::vector<Target> targets, Kdtree::KdTree* KdTreeBrain,
                     Target entreeCerveau, std::shared_ptr<geometry::TriangleMesh> mesh, Kdtree::KdTree* KdTreeSkull);

#endif //PSO_H