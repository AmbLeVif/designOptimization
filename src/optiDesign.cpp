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
#include "kdtree.hpp"
#include "tools.hpp"
#include "Workspace.h"
//#include "Dessin.hpp"
#include "Target.h"
//#include "wrkspace.hpp"

//#include "optimisateurs.hpp"
//#include "RRT.hpp"

//#include "plotScene.cpp"

using namespace open3d;
using namespace CtrLib;
using namespace Eigen;

std::shared_ptr<open3d::geometry::TriangleMesh> meshFromLines(open3d::geometry::LineSet Lines, double D, Eigen::Vector3d color)
{
// Définir l'épaisseur de la ligne

    // Ajouter des cylindres à chaque ligne pour simuler l'épaisseur
    std::shared_ptr<open3d::geometry::TriangleMesh> thick_lines = std::make_shared<open3d::geometry::TriangleMesh>();

    for (const auto& line : Lines.lines_) {
        // Récupérer les points de la ligne
        Eigen::Vector3d start_point = Lines.points_[line(0)];
        Eigen::Vector3d end_point = Lines.points_[line(1)];

        // Créer un cylindre entre les points de la ligne
        if ((start_point - end_point).norm() > 0)
        {
          //std::cout << "hauteur cylindre : " << (start_point - end_point).norm() << std::endl;
          auto cylinder = geometry::TriangleMesh::CreateCylinder(D, (start_point - end_point).norm());
          // Placer le cylindre au bon endroit et l'orienter
          Eigen::Vector3d direction = (end_point - start_point);
          direction.normalize();
          Eigen::Vector3d axis(0, 0, 1);  // L'axe de base du cylindre
          Eigen::Vector3d rotation_axis = axis.cross(direction);
          double rotation_angle = acos(axis.dot(direction));
          // Si l'angle est très proche de 0 ou 180 degrés, ne pas faire de rotation
          if (rotation_angle > 1e-6) {
              Eigen::Matrix3d rotation_matrix;
              // Créer une matrice de rotation à partir de l'axe et de l'angle
              rotation_matrix = Eigen::AngleAxisd(rotation_angle, rotation_axis.normalized()).toRotationMatrix();

              // Appliquer la rotation sur le cylindre
              cylinder->Rotate(rotation_matrix, Vector3d(0,0,0));
          }

          //Eigen::Matrix3d rotation = open3d::geometry(axis, direction.normalized());
          //cylinder->Rotate(rotation, start_point);
          cylinder->Translate((start_point + end_point) / 2.0);  // Positionner au centre de la ligne

          // Ajouter le cylindre au maillage des lignes épaisses
          *thick_lines += *cylinder;
        }
    }
    thick_lines->PaintUniformColor(color);  

    return thick_lines;
}

geometry::LineSet* LineFromCtr(VectorXd YtotLastId, MatrixXd PosCtr, int i)
{
  int sz_end = YtotLastId[i];

  MatrixXd rbis = PosCtr.topLeftCorner(NB_TUBES, sz_end);

    MatrixXd rT = rbis.transpose();

    std::vector<Eigen::Vector3d> vecT;
    std::vector<Eigen::Vector2i> idVec;

    // Parcourir chaque ligne de la matrice et l'ajouter au vecteur
    for (int i = 0; i < rT.rows(); ++i) {
        Vector3d v;
        v << rT(i, 0), rT(i, 1), rT(i, 2);
        vecT.push_back(v);
        if (i !=0)
          idVec.push_back(Vector2i(i-1, i));
    }
    return new geometry::LineSet(vecT, idVec);
}

//std::shared_ptr<geometry::TriangleMesh> plotTubes(std::vector<std::shared_ptr< geometry::Geometry3D> >* vecDisplay,VectorXd YtotLastId, MatrixXd posCtr)
std::shared_ptr<geometry::TriangleMesh> plotTubes(VectorXd YtotLastId, MatrixXd posCtr)
{
    //std::shared_ptr<geometry::TriangleMesh> robotPlot;
    std::shared_ptr<geometry::TriangleMesh> robotPlot = std::make_shared<open3d::geometry::TriangleMesh>();


    Matrix3d clrTubes {{1.0, 0.0, 0.0},
                    {0.0, 0.0, 1.0},
                    {0.5, 1.0, 0.0}};
    /*Matrix3d clrTubes {{0.0, 0.0, 0.1},
                    {0.0, 0.0, 1.0},
                    {0.5, 1.0, 0.0}};*/

  for (int i = 0; i<NB_TUBES; i++)
  {
    float r = 0.0003 * (i + 1);
    geometry::LineSet* LT = LineFromCtr(YtotLastId, posCtr, i);
    std::shared_ptr<geometry::TriangleMesh> Tube = meshFromLines(*LT, r, clrTubes.row(i));
    //vecDisplay->push_back(Tube);
    *robotPlot += *Tube;

  }
  return robotPlot;
}

std::shared_ptr<geometry::TriangleMesh> sigleTarget(Target tar, Vector3d clr, double r, bool isOriented = false)
{

  Vector3d pos = tar.position;
  Matrix3d rot = tar.orientation.block<3, 3>(0, 0);
  Vector3d dir = tar.orientation.block<3, 1>(0, 2);
  double height = 4*r;
  std::shared_ptr<geometry::TriangleMesh> point = geometry::TriangleMesh::CreateSphere(r);
  point->Translate(pos);
  point->ComputeVertexNormals();
  
  if (isOriented)
  {
    auto cylinder = geometry::TriangleMesh::CreateCylinder(r/4, height);

    Eigen::Vector3d endPoint = pos + dir.normalized() * height;
    Eigen::Vector3d translation = (pos + endPoint) / 2.0;

    Eigen::Affine3d transformation = Eigen::Affine3d::Identity();
    transformation.translation() = translation;
    transformation.linear() = rot;

    cylinder->Transform(transformation.matrix());
    cylinder->ComputeVertexNormals();

    *point += *cylinder;
  }
  point->PaintUniformColor(clr);

  return point;
}

bool isPointOnSameSide(Target target, Eigen::Vector3d point) {
    // Récupérer la normale du plan (c'est la première colonne de la matrice orientation)
    Eigen::Vector3d normal = target.orientation.block<3,1>(0,2);
    
    // Vecteur du point à tester par rapport au point sur le plan
    Eigen::Vector3d vectorToPoint = point - target.position;
    
    // Produit scalaire entre la normale et le vecteur du point
    double dotProduct = normal.dot(vectorToPoint);
    
    // Si le produit scalaire est positif, le point est du même côté que le plan
    return dotProduct > 0;
}

std::shared_ptr<geometry::LineSet> AfficheDistDura(MatrixXd posCtr, Kdtree::KdTree* KdTreeBrain, VectorXd YtotLastId, Target entreeCerveau)
{    

    std::vector<Eigen::Vector3d> vecT = {};
    std::vector<Eigen::Vector2i> idVec = {};
    std::cout << "vect size " << vecT.size() << std::endl;

    Kdtree::KdNodeVector result;

    int cpt = 0;
    for (int i = 0; i < YtotLastId[0]; i++)
    {
        Vector3d ptCtr = posCtr.col(i).topRows(3);
        if (isPointOnSameSide(entreeCerveau, ptCtr))
        //if (true)
        {
            KdTreeBrain->k_nearest_neighbors(eigenToStdVector(ptCtr), 1, &result);

            vecT.push_back(ptCtr);
            vecT.push_back(stdVectorToEigen(result[0].point));
            idVec.push_back(Vector2i(2*cpt, 2*cpt+1));
            cpt++;
        }
    }
    auto lines = new geometry::LineSet(vecT, idVec);
    lines->PaintUniformColor({0.0, 0.0, 0.0});

    return std::shared_ptr<geometry::LineSet>(lines);
}

Kdtree::KdTree KdTreeFromPointCloud(std::vector<Eigen::Vector3d> points_int)
{
  Kdtree::KdNodeVector KdVec;

  for (int i = 0; i < points_int.size(); i++)
  {
    Vector3d coord = points_int[i];
    std::vector<double> pt = {coord[0], coord[1], coord[2]};
    Kdtree::KdNode KdPt(pt);
    KdPt.index = i;
    KdVec.push_back(KdPt);

  }
  return Kdtree::KdTree(&KdVec);
}

// Fonction pour calculer la couleur de chaque point en fonction du nombre de voisins dans un rayon donné
std::vector<Eigen::Vector3d> computePointCloudColors(const std::vector<Vector3d>& points, double radius) {
    std::vector<Eigen::Vector3d> colors;
    std::vector<int> VectNbNeighbor;
    int min = -1;
    int max = 0;

    Kdtree::KdTree KdTreeWrkSpace = KdTreeFromPointCloud(points);
    Kdtree::KdNodeVector result;
    int cpt = 0;

    // Parcourir chaque point du nuage de points
    for (const auto& point : points) {

        KdTreeWrkSpace.range_nearest_neighbors(eigenToStdVector(point), radius, &result);

        // Compter le nombre de voisins dans le rayon
        int neighborCount = result.size();
        VectNbNeighbor.push_back(neighborCount);

        if (min == -1 || neighborCount < min) min = neighborCount;
        if (neighborCount > max) max = neighborCount;
        std::cout << cpt << std::endl;
        cpt++;

    }
    cpt = 0;
    std::cout << "nb voisins : " << VectNbNeighbor[21940] << std::endl;

    for (int i = 0 ; i < VectNbNeighbor.size() ; i++) {

        int nbNeighbor = VectNbNeighbor[i];
        // Normaliser le nombre de voisins pour obtenir une couleur entre 0 et 1
        // On utilise la normalisation en fonction du nombre de points dans le nuage de points
        double shift = 0.2;

        double normalizedDensity = (1.0*nbNeighbor-min+shift*max) / (max - min+shift*max);

        // La couleur peut être par exemple définie sur une échelle de rouge à bleu
        Vector3d color = getRainbow(normalizedDensity);
        /*color[0] = normalizedDensity; // Rouge
        color[1] = 0.0;               // Vert (constant, peut être ajusté si besoin)
        color[2] = 1.0 - normalizedDensity; // Bleu (inverse du rouge)*/

        colors.push_back(color);
    }
    std::cout << "min de voisins = " << min << ", max de voisins = " << max << "taille vecteur couleur : " << colors.size() << std::endl;

    return colors;
}

std::shared_ptr< geometry::PointCloud > plotWorkspace(std::vector<MatrixXd> wrkSpace, bool colorPlot = false, double radius = 0.001)
{
    geometry::PointCloud pointcloud;
    std::shared_ptr<geometry::PointCloud> pointcloud_ptr(
            new geometry::PointCloud);
    int cpt =0;
    for(const MatrixXd& mat : wrkSpace) 
    {
        Vector3d pt = mat.col(mat.cols() - 1);
        pointcloud.points_.push_back(pt);
        cpt++;
    }

    *pointcloud_ptr = pointcloud;
    pointcloud_ptr->NormalizeNormals();
    if (colorPlot)
    {
        std::vector<Eigen::Vector3d> clrs = computePointCloudColors(pointcloud.points_, radius);
        pointcloud_ptr->colors_ = clrs;
    }
    else
        pointcloud_ptr->PaintUniformColor({0, 0, 0});

    return pointcloud_ptr;
}

void PaintObjectif(std::shared_ptr< geometry::PointCloud > points, VectorXi surfaceCible, Vector3d clrObj, Vector3d clrObjReached)
{
    for (int i = 0 ; i < points->colors_.size(); i++)
    {
        if (surfaceCible[i] == 1)
        {
            points->colors_[i] = clrObj;
            //std::cout << ":/";
        }
        else if(surfaceCible[i] == 2)
        {
            points->colors_[i] = clrObjReached;
        }
    }
}

Kdtree::KdTree KdTreeFromMesh(geometry::TriangleMesh mesh)
{
  Kdtree::KdNodeVector KdVec;

  for (int i = 0; i < mesh.triangles_.size(); i++)
  {
    Vector3i idVertTriangle = mesh.triangles_[i];
    Vector3d barycentre = (mesh.vertices_[idVertTriangle[0]] +
                          mesh.vertices_[idVertTriangle[1]] +
                          mesh.vertices_[idVertTriangle[2]]) / 3.0;
    Kdtree::KdNode KdPt(eigenToStdVector(barycentre));
    KdPt.index = i;
    KdVec.push_back(KdPt);

  }
  return Kdtree::KdTree(&KdVec);
}

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
/*
WorkspaceData buildWorkspaceObstacle(CtrModel ctr, int discretAngle, int discretLength, VectorXd offset, double deltaTubes, double distUtile, bool isOriented, std::shared_ptr<geometry::TriangleMesh> mesh, Kdtree::KdTree* KdTreeSkull, bool smart = false, int freqRevolution = -1)
{
    if(freqRevolution == -1)freqRevolution = discretAngle;
    Eigen::Vector<double, NB_TUBES> lTube = ctr.GetTubeParameters().l;
    double l1, l2, l3;
    l1 = lTube[0];
    l2 = lTube[1];
    l3 = lTube[2];

    double deltaP1 = offset[0] - offset[1] - 2*deltaTubes;
    double deltaP2 = offset[1] - offset[2] - deltaTubes;
    double deltaP3 = offset[2] - deltaTubes/2;
    double deltaD12 = 0.001;
    double deltaD23 = 0.001;
    double deltaD30 = 0.001;

    double eps = 0.00001;

    WorkspaceData workspace;
    std::vector<MatrixXd> result;
    std::vector<Vector3d> YtotLastIdList;
    std::vector<Matrix3d> orientations;
    std::vector<VectorXd> listQ;

    VectorXd YtotLastId(NB_TUBES);
    Vector<double,NB_Q> q;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    int cpt = 0;
    int cptError = 0;

    //std::cout << "l1 = " << l1 << "l2 = " << l2 << "l3 = " << l3 << " angle = " << discretAngle << " length = " << discretLength << std::endl;

    CtrModel ctrSave1 = ctr;
    CtrModel ctrSave2 = ctr;
    CtrModel ctrSave3 = ctr;
    CtrModel ctrSave4 = ctr;
    CtrModel ctrSave5 = ctr;
    int nThread = omp_get_max_threads();
    double alpha = 0;
    for (int ia = 0 ; ia < (smart ? 1 : discretAngle); ia++)
    {
        alpha += 2*pi/discretAngle;
        ctrSave1 = ctr;
        double beta = 0;
        //#pragma omp parallel for num_threads(nThread)
        for(int ib = 0 ; ib < discretAngle ; ib++)
        {
            beta += 2*pi/discretAngle;
            ctrSave2 = ctr;
            double gamma = 0;
            for(int ig = 0 ; ig < discretAngle ; ig++)
            {
                gamma += 2*pi/discretAngle;
                //std::cout << "----------------------------OOOOOH, alpha = " << alpha << ", beta = " << beta << " et gqmma = " << gamma << std::endl;
                ctrSave3 = ctr;

                double start3 = -l3 + deltaD30;
                double end3 = deltaP3;
                int cpt_q3 = 0;
                for (float q3 = start3; q3 <= end3+eps; q3 += (end3-start3)/discretLength)
                {
                    //std::cout << "----------------------------iq3 = " << cpt_q3 << " et q3 = " << q3 << " et " << q3 + (end3-start3)/discretLength << " < " << end3 << std::endl;
                    //std::cout << "----------------------------iq3 = " << cpt_q3 << std::endl;
                    cpt_q3++;
                    ctrSave4 = ctr;

                    double start2 = q3 + l3 + deltaD23 -l2 + eps; 
                    double end2 = q3 + deltaP2;
                    int cpt_q2 = 0;
                    for (float q2 = start2; q2 <= end2+eps; q2 += (end2-start2)/discretLength)
                    {
                        //std::cout << "----------------------------iq2 = " << cpt_q2 << std::endl;
                        //std::cout << "----------------------------iq2 = " << cpt_q2 << " et q2 = " << q2 << " et " << q2 + (end2-start2)/discretLength << " < " << end2 << std::endl;
                        cpt_q2++;
                        ctrSave5 = ctr;
                        double start1 = q2 + l2 + deltaD12 -l1 + eps; 
                        double end1 = q2 + deltaP1;
                        int cpt_q1 = 0;
                        for (float q1 = start1; q1 <= end1+eps; q1 += (end1-start1)/discretLength)
                        {
                            //std::cout << "----------------------------iq1 = " << cpt_q1 << std::endl;
                            cpt_q1++;
                            q << q1, q2, q3, alpha, beta, gamma; // coord. articulaires (angle en rad) pi/4
                            // Compute model
                            //std::cout << "ah" << std::endl;
                            int stableShape = (ctr.Compute(q, opt_LOAD) >=0);
                            for (int i = 0; i<NB_TUBES; i++){YtotLastId[i]=getYtotIndexFromIend(ctr.segmented.iEnd(i)) + 1;}
                            //std::cout << "wiiii cols = " << ctr.GetYTot().cols() << " rows = " << ctr.GetYTot().rows() << " enfin YtotLastId= " << YtotLastId << std::endl;
                            MatrixXd formeCtr = ctr.GetYTot().topLeftCorner(3, YtotLastId[0]);
                            if (!isRobotContact(formeCtr, mesh, KdTreeSkull))
                            {
                                workspace.isStable.push_back(stableShape);
                                if (!stableShape) cptError++;
                                result.push_back(formeCtr);
                                YtotLastIdList.push_back(YtotLastId);
                                listQ.push_back(q);
                                if (isOriented)
                                {
                                    Eigen::Matrix3d orientationMatrix;
                                    for (int i = 0; i < NB_TUBES; ++i) {
                                        orientationMatrix(0, i) = ctr.GetYTot()(i + 3, getYtotIndexFromIend(ctr.segmented.iEnd(0)));    // Ligne 3 à 5 (index 2 à 4)
                                        orientationMatrix(1, i) = ctr.GetYTot()(i + 6, getYtotIndexFromIend(ctr.segmented.iEnd(0)));    // Ligne 6 à 8 (index 5 à 7)
                                        orientationMatrix(2, i) = ctr.GetYTot()(i + 9, getYtotIndexFromIend(ctr.segmented.iEnd(0)));    // Ligne 9 à 11 (index 8 à 10)
                                        orientations.push_back(orientationMatrix);
                                    }
                                }
                                else
                                    orientations.push_back(Matrix3d::Identity());
                            }
                            cpt++;
                        }
                        ctr = ctrSave5;

                    }
                    ctr = ctrSave4;
                }
                ctr = ctrSave3;
            }
            ctr = ctrSave2;
        }
        ctr = ctrSave1;
    }
    result.erase(result.begin());
    YtotLastIdList.erase(YtotLastIdList.begin());


    std::vector<Vector3d> YtotLastIdListRevolution;
    std::vector<int> isStableRevolution;
    std::vector<VectorXd> listQRevolution;
    if (freqRevolution > 1 && smart)
    {
        result = genererRevolution(result, freqRevolution);
        orientations = genererRevolution(orientations, freqRevolution);
        
        for (auto value : YtotLastIdList) {for (int i = 0; i < freqRevolution; ++i) {YtotLastIdListRevolution.push_back(value);} }
        for (auto value : workspace.isStable) {for (int i = 0; i < freqRevolution; ++i) {isStableRevolution.push_back(value);} }
        for (auto value : listQ) {
            for (int i = 0; i < freqRevolution; ++i) {
                value[3] = 2*pi/freqRevolution*(i+1);
                listQRevolution.push_back(value);
            } 
        }
        
        YtotLastIdList = YtotLastIdListRevolution;
        workspace.isStable = isStableRevolution;
        listQ = listQRevolution;
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()*1e-6 << "[s] pour cpt = " << cpt << std::endl;
    std::cout << "size data " << result.size()  << "size stable " << isStableRevolution.size() << "lastId stable " << YtotLastIdList.size() << std::endl;

    workspace.data = result;
    workspace.orientations = orientations;
    workspace.YtotLastIdList = YtotLastIdList;
    return workspace;
}*/

/*WorkspaceData buildWorkspace(CtrModel ctr, int discretAngle, int discretLength, VectorXd offset, double deltaTubes, double distUtile, bool isOriented, bool smart = false, int freqRevolution = -1)
{
    if(freqRevolution == -1)freqRevolution = discretAngle;
    Eigen::Vector<double, NB_TUBES> lTube = ctr.GetTubeParameters().l;
    double l1, l2, l3;
    l1 = lTube[0];
    l2 = lTube[1];
    l3 = lTube[2];

    double deltaP1 = offset[0] - offset[1];
    double deltaP2 = offset[1] - offset[2];
    double deltaP3 = offset[2];
    /*double deltaP1 = offset[0] - offset[1] - 2*deltaTubes;
    double deltaP2 = offset[1] - offset[2] - deltaTubes;
    double deltaP3 = offset[2] - deltaTubes/2;*/
    /*double deltaD12 = 0.001;
    double deltaD23 = 0.001;
    double deltaD30 = 0.001;

    double eps = 0.00001;

    WorkspaceData workspace;
    std::vector<MatrixXd> result;
    std::vector<Vector3d> YtotLastIdList;
    std::vector<Matrix3d> orientations;
    std::vector<VectorXd> listQ;
    VectorXd YtotLastId(NB_TUBES);
    Vector<double,NB_Q> q;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    int cpt = 0;
    int cptError = 0;

    CtrModel ctrSave1 = ctr;
    CtrModel ctrSave2 = ctr;
    CtrModel ctrSave3 = ctr;
    CtrModel ctrSave4 = ctr;
    CtrModel ctrSave5 = ctr;
    double alpha = 0;
    for (int ia = 0 ; ia < (smart ? 1 : discretAngle); ia++)
    {
        alpha += 2*pi/discretAngle;
        ctrSave1 = ctr;
        double beta = 0;
        //#pragma omp parallel for num_threads(nThread)

        for(int ib = 0 ; ib < discretAngle ; ib++)
        {
            beta += 2*pi/discretAngle;
            ctrSave2 = ctr;
            double gamma = 0;
            for(int ig = 0 ; ig < discretAngle ; ig++)
            {
                gamma += 2*pi/discretAngle;
                //std::cout << "----------------------------OOOOOH, alpha = " << alpha << ", beta = " << beta << " et gqmma = " << gamma << std::endl;
                ctrSave3 = ctr;

                double start3 = -l3 + deltaD30;
                double end3 = deltaP3;
                int cpt_q3 = 0;
                for (float q3 = start3; q3 <= end3+eps; q3 += (end3-start3)/discretLength)
                {
                    //std::cout << "----------------------------iq3 = " << cpt_q3 << " et q3 = " << q3 << " et " << q3 + (end3-start3)/discretLength << " < " << end3 << std::endl;
                    //std::cout << "----------------------------iq3 = " << cpt_q3 << std::endl;
                    cpt_q3++;
                    ctrSave4 = ctr;

                    double start2 = q3 + l3 + deltaD23 -l2 + eps; 
                    double end2 = q3 + deltaP2;
                    int cpt_q2 = 0;
                    for (float q2 = start2; q2 <= end2+eps; q2 += (end2-start2)/discretLength)
                    {
                        //std::cout << "----------------------------iq2 = " << cpt_q2 << std::endl;
                        //std::cout << "----------------------------iq2 = " << cpt_q2 << " et q2 = " << q2 << " et " << q2 + (end2-start2)/discretLength << " < " << end2 << std::endl;
                        cpt_q2++;
                        ctrSave5 = ctr;
                        double start1 = q2 + l2 + deltaD12 -l1 + eps; 
                        double end1 = q2 + deltaP1;
                        int cpt_q1 = 0;
                        for (float q1 = start1; q1 <= end1+eps; q1 += (end1-start1)/discretLength)
                        {
                            //std::cout << "----------------------------iq1 = " << cpt_q1 << std::endl;
                            cpt_q1++;
                            q << q1, q2, q3, gamma, beta, alpha; // coord. articulaires (angle en rad) pi/4
                            // Compute model
                            //std::cout << "ah" << std::endl;
                            if (ctr.Compute(q, opt_LOAD) >=0)
                            {
                                workspace.isStable.push_back(1);
                                for (int i = 0; i<NB_TUBES; i++){YtotLastId[i]=getYtotIndexFromIend(ctr.segmented.iEnd(i)) + 1;}
                                //std::cout << "wiiii cols = " << ctr.GetYTot().cols() << " rows = " << ctr.GetYTot().rows() << " enfin YtotLastId= " << YtotLastId << std::endl;

                                MatrixXd formeCtr = ctr.GetYTot().topLeftCorner(3, YtotLastId[0]);
                                result.push_back(formeCtr);
                                YtotLastIdList.push_back(YtotLastId);
                                listQ.push_back(q);
                                if (isOriented)
                                {
                                    Eigen::Matrix3d orientationMatrix;
                                    for (int i = 0; i < NB_TUBES; ++i) {
                                        orientationMatrix(0, i) = ctr.GetYTot()(i + 3, getYtotIndexFromIend(ctr.segmented.iEnd(0)));    // Ligne 3 à 5 (index 2 à 4)
                                        orientationMatrix(1, i) = ctr.GetYTot()(i + 6, getYtotIndexFromIend(ctr.segmented.iEnd(0)));    // Ligne 6 à 8 (index 5 à 7)
                                        orientationMatrix(2, i) = ctr.GetYTot()(i + 9, getYtotIndexFromIend(ctr.segmented.iEnd(0)));    // Ligne 9 à 11 (index 8 à 10)
                                        orientations.push_back(orientationMatrix);
                                    }
                                }
                                else
                                    orientations.push_back(Matrix3d::Identity());
                            }
                            else
                            {
                                cptError++;
                                workspace.isStable.push_back(0);
                            }
                            cpt++;
                        }
                        ctr = ctrSave5;

                    }
                    ctr = ctrSave4;
                }
                ctr = ctrSave3;
            }
            ctr = ctrSave2;
        }
        ctr = ctrSave1;
    }
    result.erase(result.begin());
    YtotLastIdList.erase(YtotLastIdList.begin());
    listQ.erase(listQ.begin());

    std::vector<Vector3d> YtotLastIdListRevolution;
    std::vector<int> isStableRevolution;
    std::vector<VectorXd> listQRevolution;
    if (freqRevolution > 1 && smart)
    {
        result = genererRevolution(result, freqRevolution);
        orientations = genererRevolution(orientations, freqRevolution);
        
        for (auto value : YtotLastIdList) {for (int i = 0; i < freqRevolution; ++i) {YtotLastIdListRevolution.push_back(value);} }
        for (auto value : workspace.isStable) {for (int i = 0; i < freqRevolution; ++i) {isStableRevolution.push_back(value);} }
        for (auto value : listQ) {
            for (int i = 0; i < freqRevolution; ++i) {
                value[3] = 2*pi/freqRevolution*(i+1);
                value[4] = value[4]+pi;
                value[5] = value[5]+pi;
                
                listQRevolution.push_back(value);
            } 
        }
        
        YtotLastIdList = YtotLastIdListRevolution;
        workspace.isStable = isStableRevolution;
        listQ = listQRevolution;
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()*1e-6 << "[s] pour cpt = " << cpt << std::endl;
    std::cout << "size data " << result.size()  << "size stable " << isStableRevolution.size() << "lastId stable " << YtotLastIdList.size() << std::endl;

    workspace.q = listQ;
    workspace.data = result;
    workspace.orientations = orientations;
    workspace.YtotLastIdList = YtotLastIdList;
    return workspace;
}*/

/*std::vector<MatrixXd> buildWorkspaceOriented(CtrModel ctr, int discretAngle, int discretLength, VectorXd offset, double deltaTubes, double distUtile, bool smart = false, int freqRevolution = -1)
{
    if(freqRevolution == -1)freqRevolution = discretAngle;
    Eigen::Vector<double,NB_TUBES> lTube = ctr.GetTubeParameters().l;
    double l1, l2, l3;
    l1 = lTube[0];
    l2 = lTube[1];
    l3 = lTube[2];

    double deltaP1 = offset[0] - offset[1] - 2*deltaTubes;
    double deltaP2 = offset[1] - offset[2] - deltaTubes;
    double deltaP3 = offset[2] - deltaTubes/2;
    double deltaD12 = 0.001;
    double deltaD23 = 0.001;
    double deltaD30 = 0.001;

    double eps = 0.00001;


    std::vector<MatrixXd> result;
    Vector<double,NB_Q> q;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    int cpt = 0;

    CtrModel ctrSave1 = ctr;
    CtrModel ctrSave2 = ctr;
    CtrModel ctrSave3 = ctr;
    CtrModel ctrSave4 = ctr;
    CtrModel ctrSave5 = ctr;
    double alpha = 0;
    for (int ia = 0 ; ia < (smart ? 1 : discretAngle); ia++)
    {
        alpha += 2*pi/discretAngle;
        ctrSave1 = ctr;
        double beta = 0;
        for(int ib = 0 ; ib < discretAngle ; ib++)
        {
            beta += 2*pi/discretAngle;
            ctrSave2 = ctr;
            double gamma = 0;
            for(int ig = 0 ; ig < discretAngle ; ig++)
            {
                gamma += 2*pi/discretAngle;
                //std::cout << "----------------------------OOOOOH, alpha = " << alpha << ", beta = " << beta << " et gqmma = " << gamma << std::endl;
                ctrSave3 = ctr;

                double start3 = -l3 + deltaD30;
                double end3 = deltaP3;
                int cpt_q3 = 0;
                for (float q3 = start3; q3 <= end3+eps; q3 += (end3-start3)/discretLength)
                {
                    //std::cout << "----------------------------iq3 = " << cpt_q3 << " et q3 = " << q3 << " et " << q3 + (end3-start3)/discretLength << " < " << end3 << std::endl;
                    //std::cout << "----------------------------iq3 = " << cpt_q3 << std::endl;
                    cpt_q3++;
                    ctrSave4 = ctr;

                    double start2 = q3 + l3 + deltaD23 -l2 + eps; 
                    double end2 = q3 + deltaP2;
                    int cpt_q2 = 0;
                    for (float q2 = start2; q2 <= end2+eps; q2 += (end2-start2)/discretLength)
                    {
                        //std::cout << "----------------------------iq2 = " << cpt_q2 << std::endl;
                        //std::cout << "----------------------------iq2 = " << cpt_q2 << " et q2 = " << q2 << " et " << q2 + (end2-start2)/discretLength << " < " << end2 << std::endl;
                        cpt_q2++;
                        ctrSave5 = ctr;
                        double start1 = q2 + l2 + deltaD12 -l1 + eps; 
                        double end1 = q2 + deltaP1;
                        int cpt_q1 = 0;
                        for (float q1 = start1; q1 <= end1+eps; q1 += (end1-start1)/discretLength)
                        {
                            //std::cout << "----------------------------iq1 = " << cpt_q1 << std::endl;
                            cpt_q1++;
                            q << q1, q2, q3, alpha, beta, gamma; // coord. articulaires (angle en rad) pi/4
                            // Compute model
                            //std::cout << "ah" << std::endl;
                            ctr.Compute(q, opt_LOAD);
                            //add mutex
                            MatrixXd formeCtr = ctr.GetYTot().topLeftCorner(3, getYtotIndexFromIend(ctr.segmented.iEnd(0)) + 1);

                            Eigen::Matrix3d orientationMatrix;
                            for (int i = 0; i < NB_TUBES; ++i) {
                                orientationMatrix(0, i) = ctr.GetYTot()(i + 3, getYtotIndexFromIend(ctr.segmented.iEnd(0)));    // Ligne 3 à 5 (index 2 à 4)
                                orientationMatrix(1, i) = ctr.GetYTot()(i + 6, getYtotIndexFromIend(ctr.segmented.iEnd(0)));    // Ligne 6 à 8 (index 5 à 7)
                                orientationMatrix(2, i) = ctr.GetYTot()(i + 9, getYtotIndexFromIend(ctr.segmented.iEnd(0)));    // Ligne 9 à 11 (index 8 à 10)
                            }

                            MatrixXd FinalMatrix(formeCtr.rows(), formeCtr.cols()+orientationMatrix.cols());
                            FinalMatrix << formeCtr, orientationMatrix;

                            result.push_back(formeCtr);
                            cpt++;
                        }
                        ctr = ctrSave5;

                    }
                    ctr = ctrSave4;
                }
                ctr = ctrSave3;
            }
            ctr = ctrSave2;
        }
        ctr = ctrSave1;
    }
    result.erase(result.begin());
    if (freqRevolution > 1 && smart)
        result = genererRevolution(result, freqRevolution);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()*1e-6 << "[s] pour cpt = " << cpt << std::endl;

    return result;
}*/


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

bool isConfigOk(CtrModel ctr, Eigen::VectorXd actuators)
{
    double deltaTubes = 0.005;

    Vector<double,NB_TUBES> lTube = ctr.GetTubeParameters().l;
    VectorXd offset = ctr.GetOffset();
    double l1, l2, l3;
    l1 = lTube[0];
    l2 = lTube[1];
    l3 = lTube[2];

    double q1, q2, q3;
    q1 = actuators[0];
    q2 = actuators[1];
    q3 = actuators[2];

    double deltaP1 = offset[0] - offset[1] - 2*deltaTubes;
    double deltaP2 = offset[1] - offset[2] - deltaTubes;
    double deltaP3 = offset[2] - deltaTubes/2;
    double deltaD12 = 0.001;
    double deltaD23 = 0.001;
    double deltaD30 = 0.001;

    if (q3 >= -l3 + deltaD30 && q3 <= deltaP3){
        if (q2 >= q3 + l3 + deltaD23 - l2 && q2 <= q3 + deltaP2){
            if (q1 >= q2 + l2 + deltaD12 - l1 && q1 <= q2 + deltaP1)
            {
                //std::cout << "3 : " << -l3 + deltaD30 << " < " << q3 << " < " << deltaP3 << " | 2 : " << q3 + l3 + deltaD23 - l2 << " < " << q2 << " < " << q3 + deltaP2 << " | 1 : " << q2 + l2 + deltaD12 - l1 << " < " << q1 << " < " << q2 + deltaP1 << std::endl;
                return true;
            }
        }
    }
    //std::cout << "3 : " << -l3 + deltaD30 << " < " << q3 << " < " << deltaP3 << " | 2 : " << q3 + l3 + deltaD23 - l2 << " < " << q2 << " < " << q3 + deltaP2 << " | 1 : " << q2 + l2 + deltaD12 - l1 << " < " << q1 << " < " << q2 + deltaP1 << std::endl;
    return false;

}

// Fonction pour générer un nouvel échantillon d'actionneurs (random)
Node* generateRandomConfiguration(Node* parentNode, const Eigen::VectorXd& increments, geometry::TriangleMesh* mesh, Kdtree::KdTree* KdTreeSkull) {
    
    VectorXd oldActuators = parentNode->actuators;
    bool allOnes = false;

    while (!allOnes) {

        // Sélectionner un indice aléatoire parmi ceux qui sont encore 0 dans Vec2
        VectorXi exploredChildrenIndex = parentNode->exploredChildren;
        std::vector<int> availableIndices;
        for (int i = 0; i < exploredChildrenIndex.size(); ++i) {
            if (exploredChildrenIndex[i] == 0) {
                availableIndices.push_back(i);
            }
        }

        // Sélectionner un indice au hasard dans availableIndices
        size_t randomIndex = std::rand() % availableIndices.size();
        int i = availableIndices[randomIndex];

        // Tester la fonction isConfigOk
        VectorXd tempActuators = oldActuators; // !!!!!!!!!!!!!!!!!!!!!!!!!!!!! Attention constructeur par copie ???????????????????????
        tempActuators[i] += increments[i];

        bool goodConfig = true;
        parentNode->exploredChildren[i] = 1;

        if (isConfigOk(parentNode->ctr, tempActuators)) {
            //std::cout << "La configuration est correcte pour l'indice " << i << " avec la valeur " << Vec1[i] << "\n";
            if (availableIndices.size() <=1)
                parentNode->isFullyExplored = true;

            Node* newNode = new Node(tempActuators, parentNode);
            /*if (!isInMesh(mesh, KdTreeSkull, newNode->ctr.GetP()))
            {
                return newNode;  // Si la fonction renvoie true, on s'arrête
            }
            else
            {
                goodConfig = false;
            }*/

            
        }
        else
        {
            goodConfig = false;
        }
        

        if (goodConfig == false) {
            // Vérifier si tous les éléments de Vec2 sont maintenant égaux à 1
            if (availableIndices.size() <=1)
            {
                allOnes = true;
                parentNode->isFullyExplored = true;
            }
        }    
    }

    return nullptr;
}

// Fonction pour choisir un noeud basé sur son score
Node* selectNodeBasedOnScore(const std::vector<Node*>& tree) {
    // Calculer la somme totale des scores
    double totalScore = 0.0;
    for (const auto& node : tree) {
        if (!(node->isFullyExplored))
            totalScore += node->score;
    }

    // Choisir un noeud en fonction de la probabilité proportionnelle à son score
    double randomValue = static_cast<double>(rand()) / RAND_MAX * totalScore;
    double cumulativeScore = 0.0;
    
    for (auto& node : tree) {
        if (!(node->isFullyExplored))
        {
            cumulativeScore += node->score;
            //std::cout << "Selection : totalScore = " << totalScore << " randomValue = " << randomValue << " cumulativeScore = " << cumulativeScore << std::endl;
            if (cumulativeScore >= randomValue) {
                return node;
            }
        }
    }
    
    return nullptr;  // En cas de problème, retourner null
}

double score(Vector3d RobotTip, Kdtree::KdTree* KdTreeBrain){

    Kdtree::KdNodeVector result;

    KdTreeBrain->k_nearest_neighbors(eigenToStdVector(RobotTip), 1, &result);
    double dist = 1/((RobotTip - stdVectorToEigen(result[0].point)).norm());
    return dist;
}

// Fonction pour étendre un RRT (Randomly Exploring Random Tree)
Node* growRRT(std::vector<Node*>& tree, const Eigen::VectorXd& increments, geometry::TriangleMesh* mesh, Kdtree::KdTree* KdTreeSkull, Kdtree::KdTree* KdTreeBrain) {
    // Choisir un noeud existant basé sur son score

    Node* parentNode = selectNodeBasedOnScore(tree);
    if (parentNode == nullptr) {
        return nullptr;
    }

    // Générer une nouvelle configuration d'actionneurs
    Node* newNode = generateRandomConfiguration(parentNode, increments, mesh, KdTreeSkull);
    if (newNode == nullptr)
        return nullptr;


    // Calculer le score de cette nouvelle configuration
    double newScore = score(newNode->ctr.GetP(), KdTreeBrain);
    // Ajouter le nouveau noeud à l'arbre
    newNode->score = newScore;
    tree.push_back(newNode);  // Ajouter à la liste des noeuds de l'arbre

    return newNode;
}

// Fonction principale pour créer l'arbre RRT
std::vector<MatrixXd> createRRT(int iterations, const Eigen::VectorXd increments, CtrModel ctr, double deltaTubes, geometry::TriangleMesh mesh, Kdtree::KdTree* KdTreeSkull, Kdtree::KdTree* KdTreeBrain) {

    std::vector<MatrixXd> result;

    Eigen::Vector<double,NB_TUBES> lTube = ctr.GetTubeParameters().l;
    double l1, l2, l3;
    l1 = lTube[0];
    l2 = lTube[1];
    l3 = lTube[2];

    const VectorXd offset = ctr.GetOffset();

    double deltaP1 = offset[0] - offset[1] - 2*deltaTubes;
    double deltaP2 = offset[1] - offset[2] - deltaTubes;
    double deltaP3 = offset[2] - deltaTubes/2;
    double deltaD12 = 0.001;
    double deltaD23 = 0.001;
    double deltaD30 = 0.001;

    double eps = 0.00001;

    // Créer la racine de l'arbre (configuration initiale du robot)
    Vector<double, 2*NB_TUBES> initialConfig(increments.size());
    initialConfig << deltaD30 + deltaD23 + deltaD12 -l1 + 2*eps, deltaD30 + deltaD23 -l2 + eps, -l3 + deltaD30, 0, 0, 0;  // Initialisation avec des zéros (configuration de départ)
    ctr.Compute(initialConfig, opt_LOAD);
    result.push_back(ctr.GetYTot().topLeftCorner(3, getYtotIndexFromIend(ctr.segmented.iEnd(0)) + 1));

    std::cout << "initialConfig.size() = " << initialConfig.size() << " et score initial = " << score(ctr.GetP(), KdTreeBrain) << std::endl;
    Node* root = new Node(ctr, initialConfig, nullptr, score(ctr.GetP(), KdTreeBrain));

    // Initialiser l'arbre avec la racine
    std::vector<Node*> tree;
    tree.push_back(root);

    // Itérer pour construire l'arbre
    Node* currentNode = root;
    for (int i = 0; i < iterations; ++i) {
        Node* newNode = growRRT(tree, increments, &mesh, KdTreeSkull, KdTreeBrain);  // Créer un nouveau noeud

        if (newNode != nullptr) {
            currentNode = newNode;  // Si le noeud est ajouté, on met à jour la position du noeud courant
            result.push_back(currentNode->ctr.GetYTot().topLeftCorner(3, getYtotIndexFromIend(ctr.segmented.iEnd(0)) + 1));
        }
    }

    // Retourner la racine de l'arbre
    return result;
}

//-------------------------------------------------------------------------
//-----------------------------  Fin RRT  ---------------------------------
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
//------------------------  PSO (temporaire)  -----------------------------
//-------------------------------------------------------------------------

using namespace std;
#include <vector>

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

VectorXd ObjectiveFunctionTargetPosition(std::vector<MatrixXd> wrkSpace, Vector3d target, VectorXi* NclosestId)
{
    // Nombre de points à trouver
    const int numClosest = NclosestId->size();
    VectorXd NcloseDist(numClosest);

    std::vector<std::pair<double, int>> distances;

    for (int i = 0; i < wrkSpace.size(); ++i) {
        double dist = (target - wrkSpace[i].rightCols(1)).norm();;
        distances.push_back({dist, i});
    }

    // Trier les distances en ordre croissant (distance la plus faible en premier)
    std::sort(distances.begin(), distances.end(),
              [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                  return a.first < b.first;
              });
    for (int i = 0; i < numClosest; ++i) {
        (*NclosestId)[i] = distances[i].second;
        NcloseDist[i] = distances[i].first;
    }

    // Retourner la plus petite distance obtenue
    return NcloseDist;
}

VectorXd ObjectiveFunctionObstacle(std::vector<MatrixXd> wrkSpace, VectorXi* NclosestId, Kdtree::KdTree* KdTreeBrain, Target entreeCerveau)
{
    //double minOfMaxs = 99999;
    VectorXd distToDura;
    const int numClosest = NclosestId->size();
    VectorXd NcloseDist(numClosest);
    std::vector<std::pair<double, int>> distances;

    for (int i = 0; i < numClosest; ++i) {
        double distMax = DistMaxCtrDura(wrkSpace[(*NclosestId)[i]], KdTreeBrain, entreeCerveau);
        distances.push_back({distMax, (*NclosestId)[i]});
        //if (distMax < minOfMaxs) minOfMaxs = distMax;
    }

    // Trier les distances en ordre croissant (distance la plus faible en premier)
    std::sort(distances.begin(), distances.end(),
              [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                  return a.first < b.first;
              });

    for (int i = 0; i < numClosest; ++i) {
        (*NclosestId)[i] = distances[i].second;
        distToDura.conservativeResize(i+1);
        distToDura[i] = distances[i].first;
    } 

    return distToDura;
}

VectorXd ObjectiveFunction(PsoParam PsoParam, std::vector<MatrixXd> wrkSpace, Target target, Kdtree::KdTree* KdTreeBrain, Target entreeCerveau)
{
    VectorXd result(PsoParam.nbObjFunction);
    VectorXi NclosestId(PsoParam.closeTargets);

    VectorXd minDist = ObjectiveFunctionTargetPosition(wrkSpace, target.position, &NclosestId);
    //std::cout << "test ObjectiveFunctionTargetPosition : " << NclosestId[0] << std::endl;
    VectorXd minObst = ObjectiveFunctionObstacle(wrkSpace, &NclosestId, KdTreeBrain, entreeCerveau);

    result << minDist[0], minObst[0];
    return result;
}

std::shared_ptr< geometry::TriangleMesh > PlotValidRobot(PsoParam PSOparam, std::vector<MatrixXd> wrkSpace, std::vector<Vector3d> YtotLastId, std::vector<Target> targets, Kdtree::KdTree* KdTreeBrain, Target entreeCerveau, std::shared_ptr<geometry::TriangleMesh> mesh, Kdtree::KdTree* KdTreeSkull){

    std::shared_ptr<geometry::TriangleMesh> robotsPlot = std::make_shared<open3d::geometry::TriangleMesh>();

    for (int idTarget = 0 ; idTarget < targets.size() ; idTarget++)
    {
        VectorXi NclosestId(1000);
        VectorXd DistToTarget = ObjectiveFunctionTargetPosition(wrkSpace, targets[idTarget].position, &NclosestId);

        if (DistToTarget.size() > 0)
        {
            for (int i=0 ; i < NclosestId.size() ; i++)
            {
                if (DistToTarget[i] > PSOparam.goodDist[0])
                {
                    NclosestId.conservativeResize(i+1);
                    break;
                }
            }
            //std::cout << "NclosestId " << NclosestId << std::endl;

            VectorXd DistToDura = ObjectiveFunctionObstacle(wrkSpace, &NclosestId, KdTreeBrain, entreeCerveau);
            //std::cout << "Target " << idTarget << " dist target : " << DistToTarget[0] << " dist dura : " << DistToDura[0] << std::endl;
            int cpt = 0;
            for (int i=0 ; i < NclosestId.size() ; i++)
            {
                //std::cout << "okayy DistToDura[i] = " << DistToDura[i] << " & PSOparam.goodDistToDura = " << PSOparam.goodDistToDura << std::endl;
                if (DistToTarget.size() > 0 && DistToDura[i] < PSOparam.goodDist[1])
                {
                    //std::cout << "robot isok DistToTarget[i] = " << DistToTarget[i] << " et DistToDura[i] = " << DistToDura[i] << std::endl;
                    //std::cout << "YtotLastId : " << YtotLastId[NclosestId[i]] << std::endl;
                    if (!isRobotContact(wrkSpace[NclosestId[i]], mesh, KdTreeSkull))
                    {
                        std::cout << "Selected for target " << idTarget << " : no " << NclosestId[i] << std::endl;
                        *robotsPlot += *(plotTubes(YtotLastId[NclosestId[i]], wrkSpace[NclosestId[i]]));
                        cpt++;
                    }
                    /**robotsPlot += *(plotTubes(YtotLastId[NclosestId[i]], wrkSpace[NclosestId[i]]));
                    cpt++;*/
                }
            }
            //std::cout << "cible numero " << idTarget << " nb de robots trouvés : " << cpt+1 << std::endl;

            //*robotsPlot += *(plotTubes(VectorXd YtotLastId, MatrixXd posCtr));
        }
        else
            std::cout << "Pas de configuration proche de la cible " << idTarget << std::endl;
        
    }

    return robotsPlot;
}

vector<double> setObjectiveValue(PsoParam PsoParam, vector<double> ctrParam, std::vector<Target> targets, Kdtree::KdTree* KdTreeBrain, Target entreeCerveau,  std::shared_ptr<geometry::TriangleMesh> mesh, Kdtree::KdTree* KdTreeSkull)
{
    std::vector<parameters> vParameters;
    loadParameters("../parameters/parametersMaison.csv",vParameters);
    parameters &pNominal =  vParameters[0];

    std::vector<double> v(PsoParam.nbObjective_values*targets.size(), 99999);

    pNominal.Ux = {ctrParam[0], ctrParam[1], ctrParam[2]};
    pNominal.l_k = {ctrParam[3], ctrParam[4], ctrParam[5]};
    double deltaTubes = 0.005;

    CtrModel ctr(pNominal);

    //WorkspaceData wrkSpace = buildWorkspace(ctr, PsoParam.discretLength, PsoParam.discretAngle, pNominal.offset, deltaTubes, 0.1, false, true, 2*PsoParam.discretAngle);//CtrModel ctr, int discretAngle, int discretLength, VectorXd offset, double deltaTubes, double distUtile, int freqRevolution = 1, bool smart = false
    WorkspaceData wrkSpace = buildWorkspace(ctr, PsoParam.discretAngle, PsoParam.discretLength, pNominal.offset, false, true, 2*PsoParam.discretAngle);
    wrkSpace = cleanWorkspace(wrkSpace, mesh, KdTreeSkull);
    //PlotValidRobot(PsoParam, wrkSpace.data, wrkSpace.YtotLastIdList, targets, KdTreeBrain, entreeCerveau);

    for (int i = 0 ; i < targets.size() ; i++)
    {
        
        VectorXd result = ObjectiveFunction(PsoParam, wrkSpace.data, targets[i], KdTreeBrain, entreeCerveau);
        v[2*i] = result[0];
        v[2*i+1] = result[1];
    }

    return v;
}

// Implémentation de la fonction d'initialisation des particules
void initializeParticles(vector<Particle>& particles, PsoParam param, 
                         const vector<double>& vitesse_max, std::vector<Target> targets, Kdtree::KdTree* KdTreeBrain,
                         Target entreeCerveau, std::shared_ptr<geometry::TriangleMesh> mesh, Kdtree::KdTree* KdTreeSkull) {

    int num_dimensions = param.nbDim;
    std::vector<double> limite_basse = param.lowLimit;
    std::vector<double> limite_haute = param.highLimit;

    for (int i = 0; i < param.nbParticules; ++i) {
        Particle p;
        p.position.resize(num_dimensions);
        p.velocity.resize(num_dimensions);
        p.best_position.resize(num_dimensions);
        p.objective_values.resize(2*targets.size()); // Supposons 2 objectifs
        p.best_objective_values.resize(2*targets.size());

        // Initialisation des positions et vitesses aléatoires
        for (int j = 0; j < num_dimensions; ++j) {
            p.position[j] = limite_basse[j] + (rand() % 1000) / 1000.0 * (limite_haute[j] - limite_basse[j]); // Position entre limite_basse et limite_haute
            p.velocity[j] = (rand() % 200 - 100) / 100.0 * vitesse_max[j]; // Vitesse initiale aléatoire
        }

        std::cout << "ini particule n" << i << " Ux1 = " << p.position[0] << " Ux2 = " << p.position[1] << " Ux3 = " << p.position[2] << " lk1 = " << p.position[3] << " lk2 = " << p.position[4] << " lk3 = " << p.position[5] << std::endl;

        // Calcul des valeurs objectives initiales
        p.objective_values = setObjectiveValue(param, p.position, targets, KdTreeBrain, entreeCerveau, mesh, KdTreeSkull);
        p.best_position = p.position;
        
        particles.push_back(p);
    }
}

// Fonction pour évaluer si une solution est dominée par une autre (Dominance de Pareto)
bool dominates(const std::vector<double>& a, const std::vector<double>& b) {
    bool better = false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] < b[i]) return false;
        if (a[i] > b[i]) better = true;
    }
    return better;
}

// Recherche du meilleur ensemble global basé sur la dominance de Pareto
std::vector<Particle> getParetoFront(std::vector<Particle> oldPareto, std::vector<Particle> particles) {
    std::vector<Particle> candidates;
    candidates.reserve(particles.size() + oldPareto.size()); // optimisation mémoire
    candidates.insert(candidates.end(), particles.begin(), particles.end());
    candidates.insert(candidates.end(), oldPareto.begin(), oldPareto.end());

    //particles.insert(particles.end(), oldPareto.begin(), oldPareto.end());
    vector<Particle> pareto_front;

    for (size_t i = 0; i < candidates.size(); ++i) {
        bool dominated = false;
        for (size_t j = 0; j < candidates.size(); ++j) {
            if (i != j && dominates(candidates[j].objective_values, candidates[i].objective_values)) {
                dominated = true;
                break;
            }
        }
        if (!dominated) {
            pareto_front.push_back(candidates[i]);
        }
    }

    return pareto_front;
}

// Mise à jour des positions et vitesses des particules
void updateParticle(PsoParam param, Particle& p, const vector<double>& global_best_position, 
                    const vector<double>& vitesse_max, std::vector<Target> targets, 
                    Kdtree::KdTree* KdTreeBrain, Target entreeCerveau, 
                    std::shared_ptr<geometry::TriangleMesh> mesh, Kdtree::KdTree* KdTreeSkull) {
    for (size_t i = 0; i < p.position.size(); ++i) {
        // Mise à jour de la vitesse
        p.velocity[i] = 0.5 * p.velocity[i] 
                        + 1.5 * ((rand() % 100) / 100.0) * (p.best_position[i] - p.position[i]) 
                        + 1.5 * ((rand() % 100) / 100.0) * (global_best_position[i] - p.position[i]);
        
        // Limitation de la vitesse
        if (p.velocity[i] > vitesse_max[i]) p.velocity[i] = vitesse_max[i];
        if (p.velocity[i] < -vitesse_max[i]) p.velocity[i] = -vitesse_max[i];

        // Mise à jour de la position
        p.position[i] += p.velocity[i];

        // Limitation de la position dans les limites spécifiées
        if (p.position[i] > param.highLimit[i]) p.position[i] = param.highLimit[i];
        if (p.position[i] < param.lowLimit[i]) p.position[i] = param.lowLimit[i];
        
    }
    std::cout << " Ux1 = " << p.position[0] << " Ux2 = " << p.position[1] << " Ux3 = " << p.position[2] << " lk1 = " << p.position[3] << " lk2 = " << p.position[4] << " lk3 = " << p.position[5] << std::endl;

    // Mise à jour des valeurs objectives après le déplacement
    p.objective_values = setObjectiveValue(param, p.position, targets, KdTreeBrain, entreeCerveau, mesh, KdTreeSkull);

    // Mise à jour de la meilleure position personnelle si nécessaire
    if (dominates(p.objective_values, p.best_objective_values)) {
        p.best_position = p.position;
        p.best_objective_values = p.objective_values;
    }
}

bool isGoodParticle(PsoParam param, int nbTargets, Particle p)
{
    for (int nbFct = 0 ; nbFct < param.nbObjFunction ; nbFct++)
    {
        for (int nbT = 0 ; nbT < nbTargets ; nbT++)
        {
            if (p.objective_values[nbT*param.nbObjFunction + nbFct] > param.goodDist[nbFct])
                return false;
        }
    }
    return true;
}

void updateParticlesWithThreads(std::vector<Particle>& particles,
                                 std::vector<Particle>& result,
                                 const PsoParam& param,
                                 const std::vector<double>& global_best_position,
                                 const std::vector<double>& vitesse_max,
                                 const std::vector<Target>& targets,
                                 Kdtree::KdTree* KdTreeBrain,
                                 Target entreeCerveau,
                                 std::shared_ptr<geometry::TriangleMesh> mesh, 
                                 Kdtree::KdTree* KdTreeSkull) {
    std::mutex result_mutex;
    const size_t num_threads = std::thread::hardware_concurrency();
    const size_t num_particles = particles.size();
    const size_t block_size = (num_particles + num_threads - 1) / num_threads;

    std::vector<std::thread> threads;

    for (size_t t = 0; t < num_threads; ++t) {
        size_t start = t * block_size;
        size_t end = std::min(start + block_size, num_particles);

        threads.emplace_back([&, start, end]() {
            for (size_t i = start; i < end; ++i) {
                updateParticle(param, particles[i], global_best_position, vitesse_max, targets, KdTreeBrain, entreeCerveau, mesh, KdTreeSkull);

                if (isGoodParticle(param, targets.size(), particles[i])) {
                    std::lock_guard<std::mutex> lock(result_mutex);
                    result.push_back(particles[i]);
                }
            }
        });
    }

    for (auto& t : threads) {
        t.join();
    }
}

// Fonction principale d'optimisation PSO
vector<Particle> PSO(PsoParam param, std::vector<Target> targets, Kdtree::KdTree* KdTreeBrain, Target entreeCerveau, std::shared_ptr<geometry::TriangleMesh> mesh, Kdtree::KdTree* KdTreeSkull) {
    srand(time(0));

    std::vector<Particle> particles;
    int num_dimensions = param.nbDim;

    // Calcul des vitesses maximales
    vector<double> vitesse_max(num_dimensions);
    for (int i = 0; i < num_dimensions; ++i) {
        vitesse_max[i] = param.highLimit[i] - param.lowLimit[i];
    }

    initializeParticles(particles, param, vitesse_max, targets, KdTreeBrain, entreeCerveau, mesh, KdTreeSkull);

    std::vector<Particle> pareto_front = getParetoFront({}, particles);
    vector<Particle> result = {};


    for (int iter = 0; iter < param.nbIteration; ++iter) {
        
        // Mise à jour des meilleures positions globales
        vector<double> global_best_position(num_dimensions, 0.0);
        if (!pareto_front.empty()) {
            // On prend la meilleure position du front de Pareto (par exemple, la plus proche de l'origine)
            global_best_position = pareto_front[int(pareto_front.size()/2)].position;
        }

        updateParticlesWithThreads(particles, result, param, global_best_position, vitesse_max, targets, KdTreeBrain, entreeCerveau, mesh, KdTreeSkull);

        //-------------sans multithreading
        /*for (size_t i = 0; i < particles.size(); ++i) {
            // Mise à jour des positions et vitesses des particules
            std::cout << "particule n" << i << std::endl;
            updateParticle(param, particles[i], global_best_position, vitesse_max, targets, KdTreeBrain, entreeCerveau);

            if (isGoodParticle(param, targets.size(), particles[i]))
            {
                Particle newPart = particles[i];
                result.push_back(newPart);
            }

        }*/

        //-------------avec multithreading std::execution
        /*std::vector<std::future<std::optional<Particle>>> futures;

        for (size_t i = 0; i < particles.size(); ++i) {
            futures.push_back(std::async(std::launch::async, [&, i]() -> std::optional<Particle> {
                // Mise à jour de la particule
                updateParticle(param, particles[i], global_best_position, vitesse_max, targets, KdTreeBrain, entreeCerveau);

                // Vérification de la validité
                if (isGoodParticle(param, targets.size(), particles[i])) {
                    return particles[i]; // Retourne la particule valide
                } else {
                    return std::nullopt;
                }
            }));
        }

        // Collecte des résultats valides (séquentielle)
        for (auto& fut : futures) {
            auto optPart = fut.get();
            if (optPart.has_value()) {
                result.push_back(optPart.value());
            }
        }*/
        ////-------------fin multithread

        // Recherche du meilleur ensemble Pareto
        vector<Particle> newPareto = getParetoFront(pareto_front, particles);
        pareto_front = newPareto;

        // Affichage des résultats à chaque itération
        std::cout << "Itération " << iter + 1 << ": " << endl;
        int cpt = 1;
        for (const auto& p : pareto_front) {
            std::cout << "pareto front n" << cpt << std::endl;
            std::cout << " Ux1 = " << p.position[0] << " Ux2 = " << p.position[1] << " Ux3 = " << p.position[2] << " lk1 = " << p.position[3] << " lk2 = " << p.position[4] << " lk3 = " << p.position[5] << std::endl;

            std::cout << "obj 1: " << p.objective_values[0] << " | obj 2: " << p.objective_values[1] << " | obj 3: " << p.objective_values[2] << " | obj2 4: " << p.objective_values[3] << std::endl;
            cpt++;
        }
        cout << "-----------------" << endl;
    }
    return result;
}

//-------------------------------------------------------------------------
//-----------------------------  fin PSO  ---------------------------------
//-------------------------------------------------------------------------

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