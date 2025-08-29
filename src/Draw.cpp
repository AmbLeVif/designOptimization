#include "Draw.hpp"


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

std::shared_ptr<geometry::TriangleMesh> sigleTarget(Target tar, Vector3d clr, double r, bool isOriented)
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

std::shared_ptr< geometry::PointCloud > plotWorkspace(std::vector<MatrixXd> wrkSpace, bool colorPlot, double radius)
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
