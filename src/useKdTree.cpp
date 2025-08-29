#include "useKdTree.hpp"

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