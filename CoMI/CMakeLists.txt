cmake_minimum_required(VERSION 3.16)
project(CoMI)

set(CMAKE_CXX_STANDARD 11)

find_package(Eigen3 REQUIRED)
include_directories($(EIGEN3_INCLUDE_DIRS))

find_package(pinocchio)
include_directories($(PINOCCHIO_INCLUDE_DIRS))



