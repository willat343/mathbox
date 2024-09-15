cmake_minimum_required(VERSION 3.5)
project(mathbox_catkin
    VERSION 1.1.0
    DESCRIPTION "IIO C++ catkin"
    LANGUAGES CXX
)

# Find catkin
find_package(catkin REQUIRED COMPONENTS)

# Install mathbox
include(ExternalProject)
file(MAKE_DIRECTORY ${CATKIN_DEVEL_PREFIX}/include)
ExternalProject_Add(mathbox
    PREFIX ${CMAKE_BINARY_DIR}
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../mathbox
    INSTALL_DIR ${CATKIN_DEVEL_PREFIX}
    CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=${CATKIN_DEVEL_PREFIX}
)

# Catkin package
catkin_package(
    INCLUDE_DIRS
        ${CATKIN_DEVEL_PREFIX}/include
)
