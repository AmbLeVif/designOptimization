cmake_minimum_required(VERSION 3.10)
project(designOptimization)

# The options need to be the same as Open3D's default
# If Open3D is configured and built with custom options, you'll also need to
# specify the same custom options.
option(STATIC_WINDOWS_RUNTIME "Use static (MT/MTd) Windows runtime" ON)
if(STATIC_WINDOWS_RUNTIME)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
else()
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Trouver la bibliothèque Open3D installée sur le système
find_package(Open3D REQUIRED)
find_package(CGAL REQUIRED)
find_package(OpenMP REQUIRED)

# Inclure les répertoires nécessaires pour ton projet
include_directories(${PROJECT_SOURCE_DIR}/library)  # Répertoire global des en-têtes
include_directories(${PROJECT_SOURCE_DIR}/library/CtrModel/src)
include_directories(${PROJECT_SOURCE_DIR}/library/CtrModel/Library/eigen-3.4.0)

# Ajouter le projet CtrModel en tant que sous-projet
add_subdirectory(library/CtrModel)

# Ajouter les fichiers sources de ton projet principal
file(GLOB SOURCES src/*.cpp library/*.cpp library/*.hpp)

# Créer l'exécutable de ton projet
add_executable(designOpti ${SOURCES})

# Lier Open3D et CtrModel à ton projet
target_link_libraries(designOpti PRIVATE Open3D::Open3D CtrModel CGAL::CGAL OpenMP::OpenMP_CXX)