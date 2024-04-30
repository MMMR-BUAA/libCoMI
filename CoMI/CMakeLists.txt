#
# Copyright (c) 2015-2023 CNRS INRIA
# Copyright (c) 2015 Wandercraft, 86 rue de Paris 91400 Orsay, France.
#
CMAKE_MINIMUM_REQUIRED(VERSION 3.16)
project(CoMI)
#set(CoMI_VERSION_MAJOR 1)
#set(CoMI_VERSION_MINOR 0)

#set(CMAKE_CXX_STANDARD 11)

set(CMAKE_INSTALL_PREFIX /usr/local)

#添加Eigen3和pincchio的库
find_package(PkgConfig REQUIRED)
pkg_check_modules(EIGEN3 REQUIRED eigen3)
include_directories(${eigen3_INCLUDE_DIRS})
find_package(pinocchio REQUIRED)
include_directories(${pinocchio_INCLUDE_DIRS})
link_directories(${pinocchio_LIBRARY_DIRS})

include_directories(
    include
    ${catkin_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIRS}
    ${pinocchio_INCLUDE_DIRS}
)

aux_source_directory(${CMAKE_SOURCE_DIR}/include CoMI_H_LIST)
##########################构建库
add_library(CoMI INTERFACE)
target_include_directories(CoMI INTERFACE ${CoMI_H_LIST} $<INSTALL_INTERFACE:include>)
target_link_libraries(CoMI PUBLIC ${pinocchio_LIBRARY_DIRS})
target_sources(CoMI PUBLIC ${CoMI_H_LIST})
##########################安装
install(TARGETS CoMI
        DESTINATION lib
        EXPORT CoMITargets)
install(EXPORT CoMITargets
  FILE CoMITargets.cmake
  DESTINATION lib/cmake/CoMI
)
install(DIRECTORY "${CoMI_H_LIST}include/" DESTINATION include/CoMI)
##########################生成包
#include(InstallRequiredSystemLibraries)
#set(CPACK_PACKAGE_NAME "CoMI")
#set(CPACK_DEBIAN_PACKAGE_NAME "CoMI")
#set(CPACK_RESOURE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.txt")
#set(CPACK_PACKAGE_VERSION_MAJOR "${CoMI_VERSION_MAJOR}")
#set(CPACK_PACKAGE_VERSION_MINOR "${CoMI_VERSION_MINOR}")
#set(CPACK_GENERATOR "ZIP")
#set(CPACK_SET_DESTDIR ON)
#set(CPACK_OUTPUT_FILE_PRERIX "${pinocchio_INCLUDE_DIRS}/cmake/CoMI")
#include(CPack)
include(CMakePackageConfigHelpers)
# generate the config file that is includes the exports
configure_package_config_file(${PROJECT_SOURCE_DIR}/src/Config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/CoMIConfig.cmake"
  INSTALL_DESTINATION "./"
  NO_SET_AND_CHECK_MACRO
  NO_CHECK_REQUIRED_COMPONENTS_MACRO
  )
# generate the version file for the config file
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/CoMIConfigVersion.cmake"
  VERSION "1.0"
  COMPATIBILITY AnyNewerVersion
)
# install the configuration file
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/CoMIConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/CoMIConfigVersion.cmake
  DESTINATION lib/cmake/CoMI
  )


