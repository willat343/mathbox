cmake_minimum_required(VERSION 3.10)
project(mathbox
    VERSION 0.5.0
    DESCRIPTION "Math Toolbox"
    LANGUAGES CXX
)

# Installation directories
include(GNUInstallDirs)

# CMake Build Type
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo" "MinSizeRel" )
endif()
message(STATUS "Build type set to ${CMAKE_BUILD_TYPE}")

# General Options
option(BUILD_DOCUMENTATION "Build Doxygen Documentation" OFF)
option(BUILD_TESTS "Build Tests" OFF)
option(DEVELOPER_MODE "Use local versions of cppbox" OFF)
option(USE_SYSTEM_Eigen3 "Use system install of Eigen3" ON)

# Get information about the CXX Compiler
message(STATUS "Detected ${CMAKE_CXX_COMPILER_ID} CXX compiler version ${CMAKE_CXX_COMPILER_VERSION}")

# Get information about the CXX Compiler
string(REPLACE "." ";" CMAKE_CXX_COMPILER_VERSION_LIST ${CMAKE_CXX_COMPILER_VERSION})
list(GET CMAKE_CXX_COMPILER_VERSION_LIST 0 CMAKE_CXX_COMPILER_VERSION_MAJOR)
list(GET CMAKE_CXX_COMPILER_VERSION_LIST 1 CMAKE_CXX_COMPILER_VERSION_MINOR)

# Dependency: Eigen3 (Target: Eigen3::Eigen, Note: Eigen3::Eigen sets INTERFACE_INCLUDE_DIRECTORIES)
if (NOT TARGET Eigen3::Eigen)
    set(EIGEN3_VERSION "3.3.7")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU") # Note: gcc/g++
        if (CMAKE_CXX_COMPILER_VERSION_MAJOR GREATER_EQUAL 10)
            # g++-10 with c++20 requires Eigen3 >=3.3.9. For further details see issues:
            #  - https://gitlab.com/libeigen/eigen/-/issues/2012
            #  - https://gitlab.com/libeigen/eigen/-/issues/1872
            set(EIGEN3_VERSION "3.3.9")
        endif()
    endif()
    if (${USE_SYSTEM_Eigen3})
        find_package(Eigen3 ${EIGEN3_VERSION} REQUIRED)
        message(STATUS "Found Eigen3 at ${Eigen3_DIR} with version ${Eigen3_VERSION}")
    else()
        FetchContent_Declare(
            eigen
            GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
            GIT_TAG ${EIGEN3_VERSION}
        )
        FetchContent_MakeAvailable(eigen)
        message(STATUS "Fetched eigen to ${eigen_SOURCE_DIR} with version ${EIGEN3_VERSION}")
    endif()
endif()

# Dependency: cppbox (Target: cppbox::cppbox)
if (NOT TARGET cppbox::cppbox)
    set(CPPBOX_VERSION "0.0.4")
    FetchContent_Declare(
        cppbox
        GIT_REPOSITORY git@github.com:willat343/cppbox.git
        GIT_TAG v${CPPBOX_VERSION}
    )
    FetchContent_MakeAvailable(cppbox)
    if (${DEVELOPER_MODE})
        message(STATUS "Found cppbox at ${cppbox_SOURCE_DIR} (DEVELOPER_MODE = ON)")
    else()
        message(STATUS "Fetched cppbox to ${cppbox_SOURCE_DIR} with version ${CPPBOX_VERSION}")
    endif()
endif()

# Interface Library
set(LIBRARY_NAME ${PROJECT_NAME})
add_library(${LIBRARY_NAME} INTERFACE)
add_library(${PROJECT_NAME}::${LIBRARY_NAME} ALIAS ${LIBRARY_NAME})

# Include Directories
target_include_directories(${LIBRARY_NAME} INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# Properties
set_target_properties(${LIBRARY_NAME} PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED ON
)

# Compile Options
if (CMAKE_CXX_COMPILER_VERSION_MAJOR EQUAL 9)
    # Add experimental concepts flag
    target_compile_options(${LIBRARY_NAME} INTERFACE -fconcepts)
endif()

# Link Libraries
target_link_libraries(${LIBRARY_NAME} INTERFACE
    Eigen3::Eigen
    cppbox::cppbox
)

# Configure package version and config files
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)

# Installation
install(
    TARGETS ${LIBRARY_NAME}
    EXPORT ${PROJECT_NAME}
    PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}"
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(
    EXPORT ${PROJECT_NAME}
    DESTINATION ${PROJECT_NAME}
    FILE ${PROJECT_NAME}-export.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)
install(
    FILES
        "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)
install(
    DIRECTORY include/${PROJECT_NAME}
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    FILES_MATCHING
        PATTERN "*.hpp"
        PATTERN "*.h"
)

# Uninstallation
if(NOT TARGET uninstall)
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
        @ONLY
    )

    # Add uninstall target
    add_custom_target(uninstall COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()

# Build Tests
if(BUILD_TESTS)
    add_subdirectory(test)
endif()

# Build Documentation
if (BUILD_DOCUMENTATION)
    find_package(Doxygen REQUIRED)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    configure_file(Doxyfile.in ${DOXYGEN_OUT})
    add_custom_target(documentation ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating documentation with doxygen."
        VERBATIM
    )
endif()
