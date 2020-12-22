cmake_minimum_required(VERSION 3.17)
project(DGTCM)
set(CMAKE_CXX_STANDARD 17)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
add_executable(DGTCM main.cpp Solver.cpp Solver.h PhysicalModel.h Matrix.cpp Matrix.h Thermo.h)
target_link_libraries (DGTCM Eigen3::Eigen)
