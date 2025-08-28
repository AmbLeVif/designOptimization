#include "tools.hpp"

using namespace Eigen;

double piTool = 3.1415926535897932384626433832795028841971L;

// Fonction pour convertir les angles en degrés en radians
double degToRad(double degrees) {
    return degrees * piTool / 180.0;
}

// Fonction pour créer une matrice de transformation homogène à partir des angles d'Euler et d'un vecteur de translation
Matrix4d eulerToHomogeneousMatrix(const Vector3d& eulerAnglesDeg, const Vector3d& translation) {
    // Convertir les angles d'Euler en radians
    Vector3d eulerAnglesRad;
    eulerAnglesRad[0] = degToRad(eulerAnglesDeg[0]); // roll
    eulerAnglesRad[1] = degToRad(eulerAnglesDeg[1]); // pitch
    eulerAnglesRad[2] = degToRad(eulerAnglesDeg[2]); // yaw

    // Matrices de rotation autour des axes X, Y, Z
    Matrix3d R_x, R_y, R_z;

    // Matrice de rotation autour de l'axe X (roll)
    R_x << 1, 0, 0,
           0, cos(eulerAnglesRad[0]), -sin(eulerAnglesRad[0]),
           0, sin(eulerAnglesRad[0]), cos(eulerAnglesRad[0]);

    // Matrice de rotation autour de l'axe Y (pitch)
    R_y << cos(eulerAnglesRad[1]), 0, sin(eulerAnglesRad[1]),
           0, 1, 0,
           -sin(eulerAnglesRad[1]), 0, cos(eulerAnglesRad[1]);

    // Matrice de rotation autour de l'axe Z (yaw)
    R_z << cos(eulerAnglesRad[2]), -sin(eulerAnglesRad[2]), 0,
           sin(eulerAnglesRad[2]), cos(eulerAnglesRad[2]), 0,
           0, 0, 1;

    // La matrice de rotation totale est le produit des matrices de rotation
    Matrix3d R = R_z * R_y * R_x;

    // Créer la matrice de transformation homogène 4x4
    Matrix4d transformationMatrix = Eigen::Matrix4d::Identity();
    transformationMatrix.block<3, 3>(0, 0) = R;  // Rotation
    transformationMatrix.block<3, 1>(0, 3) = translation;  // Translation

    return transformationMatrix;
}

std::vector<double> eigenToStdVector(const Eigen::VectorXd& eigenVec) {
    std::vector<double> stdVec(eigenVec.data(), eigenVec.data() + eigenVec.size());
    return stdVec;
}

Eigen::VectorXd stdVectorToEigen(const std::vector<double>& stdVec) {
    Eigen::VectorXd eigenVec(stdVec.size());
    for (size_t i = 0; i < stdVec.size(); ++i) {
        eigenVec[i] = stdVec[i];
    }
    return eigenVec;
}

double angleBetweenDirections(const Eigen::Matrix3d& m1, const Eigen::Matrix3d& m2) {
    // Calcul du produit scalaire entre les deux vecteurs

    Vector3d v1 = m1 * Eigen::Vector3d::UnitZ();  // Direction après rotation R1
    Vector3d v2 = m2 * Eigen::Vector3d::UnitZ();  // Direction après rotation R2

    double dotProduct = v1.dot(v2);
    
    // Calcul des normes des vecteurs
    double normV1 = v1.norm();
    double normV2 = v2.norm();
    
    // Calcul de l'angle en radians
    double cosTheta = dotProduct / (normV1 * normV2);
    
    // Clamping pour éviter les erreurs numériques dues à des petites imprécisions
    cosTheta = std::max(-1.0, std::min(1.0, cosTheta));
    
    // Retourne l'angle en radians
    return std::acos(cosTheta);
}

Eigen::Vector3d getRainbow(double value) {
    // Assurez-vous que la valeur est bien dans l'intervalle [0, 1]
    if (value < 0.0) value = 0.0;
    if (value > 1.0) value = 1.0;

    double r, g, b;

    if (value < 0.33) {
        // Bleu -> Cyan
        r = 0.0;
        g = (value / 0.33);
        b = 1.0;
    } else if (value < 0.66) {
        // Cyan -> Vert
        r = 0.0;
        g = 1.0;
        b = 1.0 - (value - 0.33) / 0.33;
    } else {
        // Vert -> Rouge
        r = (value - 0.66) / 0.34;
        g = 1.0 - (value - 0.66) / 0.34;;
        b = 0.0;
    }

    // Retourner un vecteur Eigen contenant les valeurs RGB
    return Eigen::Vector3d(r, g, b);
}