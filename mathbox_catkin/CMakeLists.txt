cmake_minimum_required(VERSION 3.5)
project(mathbox_catkin
    VERSION 0.3.0
    DESCRIPTION "mathbox catkin"
    LANGUAGES CXX
)

# Find catkin (which adds to CMAKE_PREFIX_PATH)
find_package(catkin REQUIRED)

# Setup
include(${CMAKE_CURRENT_SOURCE_DIR}/../mathbox/cmake/setup.cmake)

# Install mathbox
include(ExternalProject)
file(MAKE_DIRECTORY ${CATKIN_DEVEL_PREFIX}/include)
string(REPLACE ";" "|" CMAKE_PREFIX_PATH_EXT_PROJ "${CMAKE_PREFIX_PATH}")
ExternalProject_Add(mathbox
    PREFIX ${CMAKE_BINARY_DIR}
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../mathbox
    INSTALL_DIR ${CATKIN_DEVEL_PREFIX}
    BUILD_ALWAYS TRUE
    LIST_SEPARATOR | # Use alternate list separator to pass lists to CMAKE_ARGS
    CMAKE_ARGS
        -DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH_EXT_PROJ}
        -DCMAKE_INSTALL_PREFIX=${CATKIN_DEVEL_PREFIX}
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DBUILD_DOCUMENTATION=${BUILD_DOCUMENTATION}
        -DBUILD_TESTS=${BUILD_TESTS}
)

# Catkin package
catkin_package(
    INCLUDE_DIRS
        ${CATKIN_DEVEL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}
        ${mathbox_SYSTEM_INCLUDE_DIRS}
    LIBRARIES
        ${mathbox_SYSTEM_LIBRARIES}
)
