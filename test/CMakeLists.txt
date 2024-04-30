cmake_minimum_required(VERSION 3.16)
project(CoMI_test)

#set(CMAKE_CXX_STANDARD 11)

#添加Eigen3和pincchio的库
find_package(PkgConfig REQUIRED)
pkg_check_modules(EIGEN3 REQUIRED eigen3)
include_directories(${eigen3_INCLUDE_DIRS})
find_package(pinocchio REQUIRED)
include_directories(${pinocchio_INCLUDE_DIRS})
link_directories(${pinocchio_LIBRARY_DIRS})
#添加CoMI的库
find_package(CoMI REQUIRED)
include_directories(${CoMI_INCLUDE_DIRS})
#link_directories(${CoMI_LIBRARY_DIRS})

include_directories(
    include
    ${catkin_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIRS}
    ${pinocchio_INCLUDE_DIRS}
)

add_executable(pinocchio_test src/pinocchio_test.cpp)
target_link_libraries(pinocchio_test ${pinocchio_LIBRARIES})