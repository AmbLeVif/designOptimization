#ifndef TARGET_H
#define TARGET_H

#include <vector>
#include <cmath>
#include "Eigen/Dense"
#include "tools.hpp"

struct Target {
    Eigen::Vector3d position;      // Vecteur pour la position
    Eigen::Matrix4d orientation;   // Matrice pour l'orientation (matrice de rotation)

    // Constructeur avec un seul vecteur de position
    Target(Eigen::Vector3d pos) 
        : position(pos), orientation(Eigen::Matrix4d::Identity()) {}

    // Constructeur avec un vecteur de position et un vecteur d'orientation Euler
    Target(Eigen::Vector3d pos, Eigen::Vector3d orientationEuler) 
        : position(pos) {
        // Convertir l'orientation Euler en matrice de rotation
        orientation = eulerToHomogeneousMatrix(orientationEuler,{0,0,0});
    }

    // Méthode pour afficher la structure Point
    void display() const {
        std::cout << "Position: \n" << position.transpose() << std::endl;
        std::cout << "Orientation (Matrice de Rotation): \n" << orientation << std::endl;
    }
};

#endif  // TARGET_H