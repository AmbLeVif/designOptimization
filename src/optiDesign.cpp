// ----------------------------------------------------------------------------
// -                        Open3D: www.open3d.org                            -
// ----------------------------------------------------------------------------
// The MIT License (MIT)
//
// Copyright (c) 2021 www.open3d.org
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.
// ----------------------------------------------------------------------------


#include <string>
#include <cstdlib>
#include <typeinfo>
#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include "CtrModel.h"
#include "loadParameters.h"
#include "pinv.h"
#include "plotCtr.h"
#include <future>  // std::async, std::future
#include <mutex>
#include <thread>
#include <optional>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Side_of_triangle_mesh.h>

#include "open3d/Open3D.h"
//#include "../include/kdtree.hpp"
#include "useKdTree.hpp"
#include "tools.hpp"
#include "Workspace.h"
#include "Draw.hpp"
#include "Target.h"
#include "PSO.hpp"
//#include "wrkspace.hpp"

//#include "optimisateurs.hpp"
//#include "RRT.hpp"

//#include "plotScene.cpp"

using namespace open3d;
using namespace CtrLib;
using namespace Eigen;

void plotVec(visualization::Visualizer* visualizer, std::vector<std::shared_ptr< geometry::Geometry3D> > vecDisplay)
{
    for (auto geometry : vecDisplay)
    {
        visualizer->AddGeometry(geometry);
    }
}

const Vector3d zero = Vector3d::Zero();

std::shared_ptr<geometry::PointCloud> NuageDePointsDepuisObj(const std::string& filename, Vector3d clr, Vector3d translation = Vector3d::Zero(), Vector3d orientation = Vector3d::Zero(), Vector3d* center = nullptr){
    std::ifstream file(filename);  // Ouvre le fichier
    geometry::PointCloud pointcloud;
    std::shared_ptr<geometry::PointCloud> pointcloud_ptr(
            new geometry::PointCloud);

    if (!file.is_open()) {
        std::cerr << "Impossible d'ouvrir le fichier." << std::endl;
        return pointcloud_ptr;  // Retourne un vecteur vide en cas d'erreur d'ouverture
    }

    std::string line;
    int cpt = 0;
    while (std::getline(file, line)) {
        // Vérifie si la ligne commence par 'v' ou 'V'
        if (line[0] == 'v') {
            std::stringstream ss(line.substr(1));  // Ignore le 'v' et lit le reste de la ligne
            double num1, num2, num3;
            // Essaie de lire trois nombres
            if (ss >> num1 >> num2 >> num3) {
                // Ajoute le groupe de trois nombres dans le vecteur de vecteurs
                pointcloud.points_.push_back(Vector3d(num1, num2, num3));
                
            }
            cpt++;
        }
    }

    file.close();  // Ferme le fichier

    
    *pointcloud_ptr = pointcloud;

    //process pointcloud
    
    if (center != nullptr)
    {
        *center = pointcloud_ptr->GetCenter();
        pointcloud_ptr->Translate(-1*(*center));
        std::cerr << "centre cerveau :" << *center << std::endl;
    }
    pointcloud_ptr->Rotate(eulerToHomogeneousMatrix(orientation,{0,0,0}).block<3, 3>(0, 0), Vector3d(0,0,0));
    pointcloud_ptr->Translate(translation);
    pointcloud_ptr->Scale(1.0/1000, zero);
    pointcloud_ptr->NormalizeNormals();
    pointcloud_ptr->PaintUniformColor(clr);

    return pointcloud_ptr;  // Retourne le vecteur de vecteurs contenant les groupes de nombres

}

std::shared_ptr<geometry::TriangleMesh> MeshFromFile(const std::string& filename, Vector3d clr, Vector3d move = Vector3d::Zero(), Vector3d turn = Vector3d::Zero(), Vector3d turnSkull = Vector3d::Zero(), Vector3d center = {0, 0, 0}){

    std::shared_ptr<geometry::TriangleMesh> mesh = io::CreateMeshFromFile(filename);

    mesh->Rotate(eulerToHomogeneousMatrix(turnSkull,{0,0,0}).block<3, 3>(0, 0), Vector3d(0,0,0));
    mesh->Translate(-1 * center);
    mesh->Rotate(eulerToHomogeneousMatrix(turn,{0,0,0}).block<3, 3>(0, 0), Vector3d(0,0,0));
    mesh->Translate(move);
    mesh->Scale(1.0/1000, zero);
    mesh->ComputeVertexNormals();
    mesh->PaintUniformColor(clr);

    return mesh;
}


/*bool isInMesh(geometry::TriangleMesh* mesh, Kdtree::KdTree* KdTreeSkull, Vector3d point)
{
  Kdtree::KdNodeVector result;
  KdTreeSkull->k_nearest_neighbors(eigenToStdVector(point), 1, &result);

  Vector3i idVertTriangle = mesh->triangles_[result[0].index];
  Vector3d normal = mesh->triangle_normals_[result[0].index];

  Vector3d triangleVert1 = mesh->vertices_[idVertTriangle[0]];
  Vector3d triangleVert2 = mesh->vertices_[idVertTriangle[1]];
  Vector3d triangleVert3 = mesh->vertices_[idVertTriangle[2]];

  Eigen::Vector3d v1_to_point = point - triangleVert1;// Calculer le vecteur du triangle (de v1 au point
  return normal.dot(v1_to_point) < 0; // Si le produit scalaire est positif, le point est du côté de la normale 
}*/


MatrixXd transformCtr(CtrModel ctr, Vector3d theta, Vector3d X){
  Matrix4d T = eulerToHomogeneousMatrix(theta, X);
  MatrixXd YTot = ctr.GetYTot();
  YTot = YTot.topRows(3);
  
  YTot.conservativeResize(YTot.rows() + 1, YTot.cols());
  VectorXd ones = VectorXd::Ones(YTot.cols()).transpose();
  YTot.row(YTot.rows()-1) = Eigen::VectorXd::Ones(YTot.cols());

  return T*YTot;

}

double DistMeanCtrDura(MatrixXd posCtr, Kdtree::KdTree* KdTreeBrain, VectorXd YtotLastId){

    double DistTot = 0.0;
    Kdtree::KdNodeVector result;
    for (int i = 0; i < YtotLastId[0]; i++)
    {
        Vector3d ptCtr = posCtr.col(i).topRows(3);
        KdTreeBrain->k_nearest_neighbors(eigenToStdVector(ptCtr), 1, &result);
        DistTot += (ptCtr - stdVectorToEigen(result[0].point)).norm();
    }

    return DistTot/posCtr.cols();
}

double DistMaxCtrDura(MatrixXd posCtr, Kdtree::KdTree* KdTreeBrain, Target entreeCerveau){

    double distMax = 0.0;
    Kdtree::KdNodeVector result;
    for (int i = 0; i < posCtr.cols(); i++)
    {
        Vector3d ptCtr = posCtr.col(i).topRows(3);
        if (isPointOnSameSide(entreeCerveau, ptCtr))
        {
            KdTreeBrain->k_nearest_neighbors(eigenToStdVector(ptCtr), 1, &result);
            double dist = (ptCtr - stdVectorToEigen(result[0].point)).norm();
            if(dist > distMax)
                distMax = dist;
        }
    }

    return distMax;
}

VectorXi GetObjectivePointIds(std::shared_ptr< geometry::PointCloud > points, std::string path)
{
    std::shared_ptr< geometry::PointCloud > objectif = NuageDePointsDepuisObj(path, {0,0,0});

    int lengthBrainSurface = points->points_.size();
    VectorXi result(lengthBrainSurface);
    for (int i = 0; i < lengthBrainSurface; ++i) {
        bool trouve = false;
        for (int j = 0; j < objectif->points_.size(); ++j) {
            // Comparaison des vecteurs 3D avec une tolérance
            if (points->points_[i].isApprox(objectif->points_[j], 1e-6)) {
                trouve = true;
                break;
            }
        }
        result[i] = trouve ? 1 : 0;  // Si trouvé, met à 1, sinon laisse à 0
    }
    std::ofstream file("../anat/idObjectif.txt");
    if (file.is_open()) {
        file << path << std::endl;
        // Enregistrement de la taille du vecteur (optionnel, si nécessaire pour la lecture)
        file << result.size() << std::endl;
        // Enregistrement des éléments du vecteur
        for (int i = 0; i < result.size(); ++i) {
            file << result[i] << " ";
        }
        file.close();
        std::cout << "Vecteur entier enregistré dans le fichier." << std::endl;
    } else {
        std::cerr << "Erreur lors de l'ouverture du fichier." << std::endl;
    }
    return result;
}

Eigen::VectorXi getSurfaceCible(std::shared_ptr< geometry::PointCloud > points, const std::string& objectiveFile, const std::string& saveObjective) {
    std::cout << saveObjective << " et " << objectiveFile << std::endl;
    std::ifstream file(saveObjective);
    int size;
    Eigen::VectorXi v;
    std::string objectiveFileNameInFile;

    if (file.is_open()) {
        // Lecture du nom de fichier precedemment utilise pour lire les objectifs
        file >> objectiveFileNameInFile;

        if (!objectiveFileNameInFile.compare(objectiveFile)) // if objectiveFileNameInFile == objectiveFile
        {
            // Lire la taille du vecteur (optionnel, mais pratique)
            file >> size;
            v.resize(size);  // Redimensionner le vecteur pour contenir les valeurs lues
            // Lire les éléments du vecteur
            for (int i = 0; i < size; ++i) {
                file >> v[i];
            }
            file.close();
        }
        else
        {

            v = GetObjectivePointIds(points, objectiveFile);
        }
        
    } else {
        std::cerr << "Erreur lors de l'ouverture du fichier." << std::endl;
    }

    return v;
}

VectorXi UpdateReachedObjective(Kdtree::KdTree* KdTree, VectorXi surfaceCible, MatrixXd pts, double r)
{
    Kdtree::KdNodeVector result;
    for (int i=0; i < pts.cols(); i++)
    {
        KdTree->range_nearest_neighbors(eigenToStdVector(pts.col(i)), r, &result);
        for (int j = 0; j < result.size(); j++) {
            if (surfaceCible[result[j].index] == 1)
            {
                surfaceCible[result[j].index] = 2;
            }
        }
    }
    return surfaceCible;
}

double CalculateReachedObjective(VectorXi surfaceCible){
    int totObjective = 0;
    int reachedObjective = 0;
    for (int i=0; i < surfaceCible.size(); i++)
    {
        if (surfaceCible[i]!=0)
            totObjective++;
            if (surfaceCible[i]==2)
                reachedObjective++;
    }
    if (totObjective==0)
    {
        std::cout << "No objective to reached in CalculateReachedObjective()" << std::endl;
        return 0;
    }
    return reachedObjective/totObjective;
}

bool isInMesh( std::shared_ptr<geometry::TriangleMesh> mesh, Kdtree::KdTree* KdTreeSkull, Vector3d point)
{
  Kdtree::KdNodeVector result;
  KdTreeSkull->k_nearest_neighbors(eigenToStdVector(point), 1, &result);

  Vector3i idVertTriangle = mesh->triangles_[result[0].index];
  Vector3d normal = mesh->triangle_normals_[result[0].index];

  Vector3d triangleVert1 = mesh->vertices_[idVertTriangle[0]];
  Vector3d triangleVert2 = mesh->vertices_[idVertTriangle[1]];
  Vector3d triangleVert3 = mesh->vertices_[idVertTriangle[2]];

  Eigen::Vector3d v1_to_point = point - triangleVert1;// Calculer le vecteur du triangle (de v1 au point
  return normal.dot(v1_to_point) < 0; // Si le produit scalaire est positif, le point est du côté de la normale 
}

//-------------------------------------IS ROBOT IN MESH ---------------------------------------------
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef K::Point_3 Point;
namespace PMP = CGAL::Polygon_mesh_processing;

std::vector<bool> arePointsInsideMesh(MatrixXd points, Mesh* mesh) {
    

    CGAL::Side_of_triangle_mesh<Mesh, K> pointInsideTester(*mesh);

    std::vector<bool> results;
    results.reserve(points.size());

    for (int i = 0; i < points.cols() ; i++) {
        Vector3d pt = points.col(i);
        Point p(pt.x(), pt.y(), pt.z());
        auto side = pointInsideTester(p);
        results.push_back(side == CGAL::ON_BOUNDED_SIDE || side == CGAL::ON_BOUNDARY);
    }

    return results;
}
//---------------------------------END IS ROBOT IN MESH ---------------------------------------------

bool isRobotContact(MatrixXd robot, std::shared_ptr<geometry::TriangleMesh> mesh, Kdtree::KdTree* KdTreeSkull)
{
    /*for(int j = 0 ;  j < robot.cols() ; j++)
    {
        Vector3d pt = robot.col(j);
        std::cout << "isInMesh = " << isInMesh(mesh, KdTreeSkull, pt) << " et pt = " << pt << std::endl;
        if(isInMesh(mesh, KdTreeSkull, pt))
        {
            return true;
        }
    }
    return false;*/
    bool isContact = false;
    for(int j = 0 ;  j < robot.cols() ; j++)
    {
        Vector3d pt = robot.col(j);
        //std::cout << "isInMesh = " << isInMesh(mesh, KdTreeSkull, pt)/* << " et pt = " << pt */<< std::endl;
        if(isInMesh(mesh, KdTreeSkull, pt))
        {
            isContact = true;
            break;
        }
    }
    return isContact;
}

WorkspaceData cleanWorkspace(WorkspaceData wrksp,  std::shared_ptr<geometry::TriangleMesh> mesh, Kdtree::KdTree* KdTreeSkull)
{
    int cpt = 0;
    WorkspaceData newWrksp;
    for(int i = 0 ; i < wrksp.data.size() ; i++)
    //for(int i = 1000 ; i < 1001 ; i++)
    {
        //std::cout << "i  = " << i<< std::endl;
        //std::cout << "data[i] = " << wrksp.data[i].cols() << " et size = " << wrksp.data[i].size()<< std::endl;
        if (!isRobotContact(wrksp.data[i], mesh, KdTreeSkull))
        {
            cpt++;
            std::cout << "clean = " << cpt << "/" << wrksp.data.size() << std::endl;
            /*wrksp.data.erase(wrksp.data.begin() + i);
            wrksp.YtotLastIdList.erase(wrksp.YtotLastIdList.begin() + i);
            wrksp.orientations.erase(wrksp.orientations.begin() + i);
            wrksp.isStable.erase(wrksp.isStable.begin() + i);*/
            newWrksp.data.push_back(wrksp.data[i]);
            newWrksp.YtotLastIdList.push_back(wrksp.YtotLastIdList[i]);
            newWrksp.orientations.push_back(wrksp.orientations[i]);
            newWrksp.isStable.push_back(wrksp.isStable[i]);
        }
    }
    std::cout << "cpt clean = " << cpt << std::endl;
    return newWrksp;
}

int main(int argc, char *argv[]) {

    std::cout << __cplusplus << std::endl;
    srand(static_cast<unsigned int>(time(0)));  // Initialisation du générateur de nombres aléatoires

    Vector3d lTubesNecessary = {0.1, 0.06, 0.03};
    double deltaTubes = 0.005;

    //Test robot ------------------------------
    std::vector<parameters> vParameters;
    //if(loadParameters("../parameters/parametersMaison.csv",vParameters) != 0){
    if(loadParameters("../parameters/parametersMaison2.csv",vParameters) != 0){
    //if(loadParameters("../parameters/parameters.csv",vParameters) != 0){
      return -1;
    }
    parameters &pNominal =  vParameters[0];

    //pNominal.l = setLength(lTubesNecessary, pNominal.offset, deltaTubes);
    double l1 = pNominal.l[0];
    double l2 = pNominal.l[1];
    double l3 = pNominal.l[2];
    std::cout << "l1 : " << pNominal.l[0] << " l2 : " << pNominal.l[1] << " l3 : " << pNominal.l[2] << std::endl;

    CtrModel ctr(pNominal);


    Vector<double,NB_Q> q;
    //q << -0.400, -0.290, -0.180, -pi/2, pi/2, pi/2; // coord. articulaires (angle en rad) pi/4
    //q << -0.06,-0.04,-0.02, 0, pi, 0; // coord. articulaires (angle en rad) pi/4
    q << -l1+0.003,-l2+0.002,-l3+0.001, 0, pi, 0;
    // Compute model
    ctr.Compute(q, opt_LOAD);

    //MatrixXd posCtr = transformCtr(ctr, {0, 0, 0}, {0,0,0});
    MatrixXd posCtr = ctr.GetYTot().topRows(3);

    VectorXd YtotLastId(NB_TUBES);
    for (int i = 0; i<NB_TUBES; i++){YtotLastId[i]=getYtotIndexFromIend(ctr.segmented.iEnd(i)) + 1;}
    
    // Get end-effector position
    Vector3d boutDuRobot = posCtr.col(posCtr.cols() - 1);

    //Test affich -----------------------------

    std::vector<std::shared_ptr< geometry::Geometry3D> > vecDisplay;

    visualization::Visualizer visualizer;
    if (argc == 2) {
        std::string option(argv[1]);
        if (option == "--skip-for-unit-test") {
            open3d::utility::LogInfo("Skiped for unit test.");
            return 0;
        }
    }
    Kdtree::KdNodeVector KdVecBrain;
    //----------------------------- Decalage ---------------------------------------------
    //Vector3d brainRotation = {-90, 0, 0}; Vector3d brainTranslation = {0, -20, 70}; // paramètre boule
    //Vector3d brainRotation = {-90, 0, 0}; Vector3d brainTranslation = {-55, -60, 170}; // paramètre boule temporaire
    Vector3d brainRotation = {210, 0, 0}; Vector3d brainTranslation = {-6, 2, 50}; // {rouge, vert, bleu} // paramètre minipig
    //Vector3d brainRotation = {0, 0, 0}; Vector3d brainTranslation = {0, 0, 0}; // 0

    //----------------------------- Points cerveau ---------------------------------------

    Eigen::Vector3d centre_cerveau = {0,0,0};

    std::shared_ptr< geometry::PointCloud > pts_objectif = NuageDePointsDepuisObj("../anat/minipig/aires_cible.obj", {1, 0.6, 0.3}, brainTranslation, brainRotation);
    std::shared_ptr< geometry::PointCloud > points_int = NuageDePointsDepuisObj("../anat/minipig/brain_minipig.obj", {1, 0.8, 0.2}, brainTranslation, brainRotation, &centre_cerveau);
    pts_objectif->Translate((0.001*centre_cerveau));
    pts_objectif->Translate({0.011, -0.005, -0.012});
    pts_objectif->Rotate(eulerToHomogeneousMatrix({-3, 0, 0}, {0,0,0}).block<3, 3>(0, 0), Vector3d(0,0,0));
    //pts_objectif->Translate({0.01,-0.01,-0.01});
    //std::shared_ptr< geometry::PointCloud > points_int = NuageDePointsDepuisObj("../anat/boule10mm.obj");
    //VectorXi surfaceCible = getSurfaceCible(points_int, "../anat/minipig/minipig_pts_objective_realistic2.obj", "../anat/idObjectif.txt");
    VectorXi surfaceCible = getSurfaceCible(points_int, "../anat/minipig/objective.obj", "../anat/idObjectif.txt");
    
    Kdtree::KdTree KdTreeBrain = KdTreeFromPointCloud(points_int->points_);

    std::cout << "bout du robot : " << boutDuRobot << std::endl;
    Kdtree::KdNodeVector result;

    KdTreeBrain.k_nearest_neighbors(eigenToStdVector(boutDuRobot), 1, &result);
    //print_nodes(result);
    std::cout << "id resultat : " << result[0].index << std::endl;

    //------------------------------- Origine -----------------------------------------
    auto origin = geometry::TriangleMesh::CreateCoordinateFrame(); //repère
    origin->Scale(1.0/100, zero);
    vecDisplay.push_back(origin);

    //-------------------------------- Cible ------------------------------------------

    Target cible1({0.01, 0.016, 0.047}, {0, 0, 0});
    Target cible2({0.0095, -0.009, 0.032}, {90, 0, 0});

    Target cibleTest({-1.79804e-16, -0.014368, 0.0962135}, {0, 0, 0});

    std::vector<Target> targets = {cible1, cible2};
    std::cout << "dist à la cible = " << (cible1.position - boutDuRobot).norm() << std::endl;
    std::cout << "rot entre les 2 cibles = " << angleBetweenDirections(cible1.orientation.block<3, 3>(0, 0), cible2.orientation.block<3, 3>(0, 0)) << std::endl;

    Target entreeCerveau({0, 0, 0.025}, {0, 0, 0});
    //vecDisplay.push_back(sigleTarget(entreeCerveau, {0, 0.5, 0}, 0.001, true)); // true = oriented
    vecDisplay.push_back(sigleTarget(cible1, {0.9, 0.9, 0.5}, 0.001));
    vecDisplay.push_back(sigleTarget(cible2, {0.9, 0.9, 0}, 0.001));

    vecDisplay.push_back(sigleTarget(cibleTest, {0.9, 0.9, 0}, 0.001));

    std::cout << "Cible du bon côte = " << isPointOnSameSide(entreeCerveau, cible1.position) << std::endl;

    //----------------------------- maillage crâne ------------------------------------

    //Vector3d skullTranslation = {-1,-24,-22}; Vector3d skullRotation = {-90, 0, 0};
    Vector3d skullTranslation = {0,0,0}; Vector3d skullRotation = {-90, 0, 0};
    //std::string meshSkullFilePath = "../anat/minipig/demi_cranecoupe_droit_cranio_centered_resimplified.stl";
    std::string meshSkullFilePathShow = "../anat/minipig/demi_cranecoupe_droit_cranio_centered.stl";

    std::shared_ptr<geometry::TriangleMesh> meshToShow = MeshFromFile(meshSkullFilePathShow, {0.7, 0.7, 0.7}, brainTranslation+skullTranslation , brainRotation, skullRotation, centre_cerveau);
    
    std::string meshSkullFilePath = "../anat/minipig/demi_cranecoupe_droit_cranio_centered_resimplified.stl";
    std::shared_ptr<geometry::TriangleMesh> mesh = MeshFromFile(meshSkullFilePath, {0.7, 0.7, 0.7}, brainTranslation+skullTranslation , brainRotation, skullRotation, centre_cerveau);
    
    Kdtree::KdTree KdTreeSkull = KdTreeFromMesh(*mesh);

    Mesh meshSkull;
    if (!PMP::IO::read_polygon_mesh(meshSkullFilePath, meshSkull) || CGAL::is_empty(meshSkull) || !CGAL::is_triangle_mesh(meshSkull)) {
        throw std::runtime_error("Erreur lors du chargement du maillage STL ou maillage invalide.");
    }


    //PSO(num_dimensions, limite_basse, limite_haute);

    //--------------------------- Collision ----------------------

    /*Eigen::Vector3d NormalClosestTriangle;
    MatrixXd robotPos = ctr.GetYTot().topLeftCorner(3, YtotLastId[0]);
    for (int i = 0 ; i < robotPos.cols() ; i++)
    {
        Target cible(robotPos.col(i), {0, 0, 0});
        vecDisplay.push_back(sigleTarget(cible, {1, 0, 0.5}, 0.001));
    }
    //bool oui = isInMesh(mesh, &KdTreeSkull, boutDuRobot);
    //std::cout << "le bout du robot est dans le crâne ? : " << isInMesh(mesh, &KdTreeSkull, boutDuRobot) << std::endl;
    std::chrono::steady_clock::time_point TimerBegin = std::chrono::steady_clock::now();
    std::cout << "Le robot est en collision avec le crâne ? : " << isRobotContact(robotPos, mesh, &KdTreeSkull) << std::endl;
    std::chrono::steady_clock::time_point TimerEnd = std::chrono::steady_clock::now();
    std::cout << "Temps pour collision = " << std::chrono::duration_cast<std::chrono::microseconds>(TimerEnd - TimerBegin).count()*1e-3 << std::endl;
    */

    /*std::chrono::steady_clock::time_point TimerBegin2 = std::chrono::steady_clock::now();
    std::vector<bool> oui = arePointsInsideMesh(robotPos, &meshSkull);
    //std::cout << " points in mesh : " << oui << std::endl;
    for (size_t i = 0; i < oui.size(); ++i) {
        std::cout << "Point " << i << (oui[i] ? " est à l'intérieur" : " est à l'extérieur") << std::endl;
    }
    std::chrono::steady_clock::time_point TimerEnd2 = std::chrono::steady_clock::now();
    std::cout << "Temps pour collision = " << std::chrono::duration_cast<std::chrono::microseconds>(TimerEnd2 - TimerBegin2).count()*1e-3 << std::endl;*/

    //--------------------------- Dist Dura ----------------------

    /*vecDisplay.push_back(AfficheDistDura(posCtr, &KdTreeBrain, YtotLastId, entreeCerveau)); // distance dura
    std::cout << "Distance max du robot a la dura : " << DistMaxCtrDura(posCtr, &KdTreeBrain, entreeCerveau) << std::endl;
    surfaceCible = UpdateReachedObjective(&KdTreeBrain, surfaceCible, boutDuRobot, 0.005);
    PaintObjectif(points_int, surfaceCible, {1, 0.6, 0.6}, {0,1,0});*/

    //--------------------------- Workspace ----------------------

    /*WorkspaceData wrkSpace = buildWorkspace(ctr, 8, 8, pNominal.offset, deltaTubes, 0.1, true, true, 16);//CtrModel ctr, int discretAngle, int discretLength, VectorXd offset, double deltaTubes, double distUtile, int freqRevolution = 1, bool smart = false
    //WorkspaceData wrkSpace = buildWorkspace(ctr, 8, 8, pNominal.offset, deltaTubes, 0.1, true);
    //WorkspaceData wrkSpace = buildWorkspaceObstacle(ctr, 4, 4, pNominal.offset, deltaTubes, 0.1, true, mesh, &KdTreeSkull, true); 
    //wrkSpace = cleanWorkspace(wrkSpace, mesh, &KdTreeSkull);
    std::cout << "taille wrksp : " << wrkSpace.data.size() << std::endl;
    std::cout << "orientation 100 : " << wrkSpace.orientations[100] << std::endl;
    //vecDisplay.push_back(plotWorkspace(wrkSpace.data, false, 0.005));*/

    //--------------------------- RRT ----------------------
    
    // Initialisation de l'arbre RRT avec q_max = 10, max_step_size = 1
    /*Vector<double, 2*NB_TUBES> incrementsRRT;
    incrementsRRT << 0.002, 0.002, 0.002, 0.5, 0.5, 0.5;
    std::vector<MatrixXd> wrkspaceRRT = createRRT(10000, incrementsRRT, ctr, deltaTubes, *mesh, &KdTreeSkull, &KdTreeBrain);
    vecDisplay.push_back(plotWorkspace(wrkspaceRRT, false, 0.03));*/

    // Construction de l'arbre RRT avec 1000 itérations
    /*rrt.buildRRT(1000); */

    //TEST PSO-----------------------------------------------------------

    PsoParam PSOparam;
    PSOparam.nbDim = 6;
    PSOparam.nbParticules = 5;
    PSOparam.nbIteration = 10;

    PSOparam.lowLimit = {0, 0, 0, 0.0001, 0.0001, 0.0001};
    PSOparam.highLimit = {400, 400, 400, 0.05, 0.05, 0.05};

    PSOparam.nbObjFunction = 2;
    PSOparam.closeTargets = 10;
    PSOparam.discretLength = 10;
    PSOparam.discretAngle = 2;
    PSOparam.nbObjective_values = 2;
    
    //PSOparam.goodDistToTarget = 0.001;
    //PSOparam.goodDistToDura = 0.002;
    PSOparam.goodDist = {0.001, 0.002};

    /*std::chrono::steady_clock::time_point TimerBegin = std::chrono::steady_clock::now();
    vector<Particle> PSOout = PSO(PSOparam, targets, &KdTreeBrain, entreeCerveau, mesh, &KdTreeSkull);
    std::chrono::steady_clock::time_point TimerEnd = std::chrono::steady_clock::now();
    std::cout << "Temps pour PSO = " << std::chrono::duration_cast<std::chrono::microseconds>(TimerEnd - TimerBegin).count()*1e-6 << std::endl;
    std::cout << "resultats : " << endl;
    int cpt = 1;
    for (const auto& p : PSOout) {
        std::cout << "particle n" << cpt << std::endl;
        std::cout << " Ux1 = " << p.position[0] << " Ux2 = " << p.position[1] << " Ux3 = " << p.position[2] << " lk1 = " << p.position[3] << " lk2 = " << p.position[4] << " lk3 = " << p.position[5] << std::endl;

        std::cout << "obj 1: " << p.objective_values[0] << " | obj 2: " << p.objective_values[1] << " | obj 3: " << p.objective_values[2] << " | obj2 4: " << p.objective_values[3] << std::endl;
        cpt++;
    }
    cout << "-----------------" << endl;*/

    /*WorkspaceData wrkSpace = buildWorkspace(ctr, PSOparam.discretAngle, PSOparam.discretLength, pNominal.offset, deltaTubes, 0.1,false, true);//CtrModel ctr, int discretAngle, int discretLength, VectorXd offset, double deltaTubes, double distUtile, int freqRevolution = 1, bool smart = false
    //wrkSpace = cleanWorkspace(wrkSpace, mesh, &KdTreeSkull);
    std::cout << "----------wouhouuuu-------------" << std::endl;
    std::cout << "taille wrksp : " << wrkSpace.data.size() << std::endl;
    vecDisplay.push_back(plotWorkspace(wrkSpace.data, false, 0.005));*/
    //vecDisplay.push_back(PlotValidRobot(PSOparam, wrkSpace.data, wrkSpace.YtotLastIdList, targets, &KdTreeBrain, entreeCerveau, mesh, &KdTreeSkull));

    //std::vector<int> robotToPlot = {456507, 376536, 130397};
    //vecDisplay.push_back(plotTubes(wrkSpace.YtotLastIdList[robotToPlot[0]], wrkSpace.data[robotToPlot[0]]));
    //vecDisplay.push_back(plotTubes(wrkSpace.YtotLastIdList[robotToPlot[1]], wrkSpace.data[robotToPlot[1]]));


    visualizer.CreateVisualizerWindow("Open3D", 1600, 900);
   /*
    //q cible 1 :
q << -0.0929262,-0.068305,-0.0378125, 5.02, 40.45, 20;

    //q cible 2
q << -0.112421,-0.068305,-0.0378125, 3.53429, 32.5013, 32.2013;*/

    //q cible 
q << -0.102421,-0.078305,-0.0478125, pi/2, -pi/2, -pi/2;
    

    
 
    std::cout << q << std::endl;
    ctr.Compute(q, opt_LOAD);
    posCtr = ctr.GetYTot().topRows(3);
    for (int i = 0; i<NB_TUBES; i++){YtotLastId[i]=getYtotIndexFromIend(ctr.segmented.iEnd(i)) + 1;}
    std::cout << posCtr.col(YtotLastId[0]-1) << std::endl;
    //std::cout << YtotLastId << std::endl;
    auto afficheRobot = plotTubes(YtotLastId, posCtr);
    visualizer.AddGeometry(afficheRobot);

    /*double delta2 = q_final[1] - q_init[1];
    for (int i = 0 ; i < discretAnim ; i++){
        q[0] += delta2/discretAnim;
        q[1] += delta2/discretAnim;

        ctr.Compute(q, opt_LOAD);
        posCtr = ctr.GetYTot().topRows(3);
        for (int i = 0; i<NB_TUBES; i++){YtotLastId[i]=getYtotIndexFromIend(ctr.segmented.iEnd(i)) + 1;}
        //afficheRobot = plotTubes(YtotLastId, posCtr);
        maxdist.push_back(DistMaxCtrDura(posCtr, &KdTreeBrain, entreeCerveau));
        if (i == 5)
        {
            posCtr = ctr.GetYTot().topRows(3);
            for (int i = 0; i<NB_TUBES; i++){YtotLastId[i]=getYtotIndexFromIend(ctr.segmented.iEnd(i)) + 1;}
            std::cout << posCtr.col(YtotLastId[0]-1) << std::endl;
            std::cout << "uiiiiiiii" << std::endl;
            auto afficheRobot = plotTubes(YtotLastId, posCtr);
            visualizer.AddGeometry(afficheRobot);
        }

    }   */ 

    plotVec(&visualizer, vecDisplay);
    visualizer.AddGeometry(points_int);//20.4   28.1    26.1    0.15    0.1 0.06    0.011   0.001   0.026
    visualizer.AddGeometry(pts_objectif);
    visualizer.AddGeometry(meshToShow);
    //visualization::DrawGeometries({points_int});
    visualizer.Run();

    /*int discretAnim = 10;
    std::vector<double> maxdist;
    Vector<double,NB_Q> q_init;
    Vector<double,NB_Q> q_final;
    q_final << -0.0929262,-0.068305,-0.0378125, 5.02, 40.45, 20;
    q_init << q_final[2] + l3 - l1 +0.0002, q_final[2] + l3 - l2 + 0.0001, q_final[2], q_final[3], q_final[4], q_final[5];
    q = q_init;
    double delta1 = q_final[0] - q[0];
    for (int i = 0 ; i < discretAnim ; i++){
        q[0] += delta1/discretAnim;

        ctr.Compute(q, opt_LOAD);
        posCtr = ctr.GetYTot().topRows(3);
        for (int i = 0; i<NB_TUBES; i++){YtotLastId[i]=getYtotIndexFromIend(ctr.segmented.iEnd(i)) + 1;}
        //afficheRobot = plotTubes(YtotLastId, posCtr);
        maxdist.push_back(DistMaxCtrDura(posCtr, &KdTreeBrain, entreeCerveau));
        
        posCtr = ctr.GetYTot().topRows(3);
        for (int i = 0; i<NB_TUBES; i++){YtotLastId[i]=getYtotIndexFromIend(ctr.segmented.iEnd(i)) + 1;}
        std::cout << posCtr.col(YtotLastId[0]-1) << std::endl;
        auto afficheRobot = plotTubes(YtotLastId, posCtr);
        if (i == 1)
            visualizer.AddGeometry(afficheRobot);
        else
            visualizer.UpdateGeometry(afficheRobot);
        visualizer.PollEvents();
        visualizer.UpdateRender();
        std::this_thread::sleep_for(std::chrono::milliseconds(1000));

    }
    //std::cout << maxdist << std::endl;
    for (size_t i = 0; i < maxdist.size(); ++i) {
        std::cout << maxdist[i] << " ";
    }*/

    visualizer.DestroyVisualizerWindow();

    /*visualization::Visualizer vis;
    vis.CreateVisualizerWindow("Open3D Animation Example", 800, 600);

    auto sphere = geometry::TriangleMesh::CreateSphere(0.5);
    sphere->PaintUniformColor(Eigen::Vector3d(0.1, 0.2, 0.9));
    vis.AddGeometry(sphere);

    std::this_thread::sleep_for(std::chrono::milliseconds(3000));
    for (int i = 0; i < 100; ++i) {
        sphere->Translate(Eigen::Vector3d(0.01, 0, 0), true);
        vis.UpdateGeometry(sphere);
        vis.PollEvents();
        vis.UpdateRender();
        std::this_thread::sleep_for(std::chrono::milliseconds(30));
    }

    vis.DestroyVisualizerWindow();*/

    return 0;
}