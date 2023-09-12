# Copyright 2019, 2020, Collabora, Ltd.
# Copyright 2019, 2021, Visual Computing Lab, ISTI - Italian National Research Council
# SPDX-License-Identifier: BSL-1.0

cmake_minimum_required(VERSION 3.18)
project(MeshLab)
### Dev Settings
IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(CMAKE_INSTALL_PREFIX D:/application/MeshLab CACHE PATH comms FORCE)
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
set(Qt5_DIR D:/libraries/Qt/5.15.2/msvc2019_64/lib/cmake/Qt5)
set(EIGEN3_INCLUDE_DIR D:/libraries/Eigen3/include/eigen3)

### Build options
option(MESHLAB_BUILD_MINI "Build only common and meshlab - other plugin targets must be set manually" OFF)
option(MESHLAB_BUILD_STRICT "Strictly enforce resolution of all symbols" ON)
option(MESHLAB_BUILD_WITH_DOUBLE_SCALAR "Use double type instead of float type for scalars" OFF)
option(MESHLAB_ENABLE_DEBUG_LOG_FILE "If enabled, all the logs of MeshLab will be also saved into a log file" OFF)

option(MESHLAB_BUILD_ONLY_LIBRARIES "Build only meshlab-common and plugins, excluding executables" OFF)
option(MESHLAB_USE_DEFAULT_BUILD_AND_INSTALL_DIRS "If set to OFF, it expects that you set manually the binary and install directories" ON)

option(MESHLAB_IS_NIGHTLY_VERSION "Nightly version of meshlab will be used instead of ML_VERSION" OFF)

add_subdirectory(src)