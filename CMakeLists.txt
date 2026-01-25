cmake_minimum_required(VERSION 3.21)
project(mathbox
    VERSION 0.7.0
    DESCRIPTION "Math Toolbox"
    LANGUAGES CXX
)

## Include Standard CMake Packages
# GNUInstallDirs for Install directories
include(GNUInstallDirs)
# FetchContent for dependencies
include(FetchContent)
# CTest for BUILD_TESTING option (default ON)
include(CTest)
# Package version and config file generation
include(CMakePackageConfigHelpers)

# General Options
option(BUILD_DOCUMENTATION "Build Doxygen Documentation" OFF)
option(BUILD_SHARED_LIBS "Build shared library (if MATHBOX_HEADER_ONLY=OFF)" OFF)
option(DEVELOPER_MODE "Use local versions of cmakebox, cppbox" OFF)
option(MATHBOX_HEADER_ONLY "Build library as header-only" OFF)

# Configure FetchContent source directories for DEVELOPER_MODE
if (DEVELOPER_MODE)
    set(FETCHCONTENT_SOURCE_DIR_CMAKEBOX "${CMAKE_SOURCE_DIR}/dependencies/cmakebox")
    set(FETCHCONTENT_SOURCE_DIR_CPPBOX "${CMAKE_SOURCE_DIR}/dependencies/cppbox")
endif()

# Dependency: cmakebox
set(CMAKEBOX_VERSION "0.0.9")
FetchContent_Declare(
    cmakebox
    GIT_REPOSITORY git@github.com:willat343/cmakebox.git
    GIT_TAG v${CMAKEBOX_VERSION}
)
FetchContent_MakeAvailable(cmakebox)
list(APPEND CMAKE_MODULE_PATH "${cmakebox_SOURCE_DIR}/cmake")
include(CMakeBox)

# CMake Build Type
setup_build_type()

# Get CXX Compiler Information
get_cxx_compiler_info()

# Dependency: Eigen3 (Target: Eigen3::Eigen)
set(EIGEN3_VERSION "3.4.0")
import_Eigen3(
    VERSION ${EIGEN3_VERSION}
)

# Dependency: cppbox (Target: cppbox::cppbox)
set(CPPBOX_VERSION "0.1.5")
import_dependency(
    cppbox
    TARGET cppbox::cppbox
    METHOD FETCH_GIT
    FIND_PACKAGE_VERSION ${CPPBOX_VERSION}
    GIT_REPOSITORY git@github.com:willat343/cppbox.git
    GIT_TAG v${CPPBOX_VERSION}
)

## Build Library
set(LIBRARY_NAME ${PROJECT_NAME})
if (MATHBOX_HEADER_ONLY)
    # Interface Library
    add_library(${LIBRARY_NAME} INTERFACE)
    set(PUBLIC_OR_INTERFACE INTERFACE)
else()
    # Shared/Static Library (determined by BUILD_SHARED_LIBS)
    add_library(${LIBRARY_NAME}
        src/covariance.cpp
        src/decompose.cpp
        src/geometry.cpp
        src/lerp.cpp
        src/matrix_operations.cpp
        src/normal.cpp
        src/nsphere.cpp
        src/statistics.cpp
        src/stiffness.cpp
        src/time.cpp
        src/transform_components.cpp
        src/transform_manager.cpp
        src/vector_operations.cpp
    )
    set(PUBLIC_OR_INTERFACE PUBLIC)

    # Private Compile Options
    target_compile_options(${LIBRARY_NAME} PRIVATE
        -Wall -Wextra -Wpedantic
    )
endif()
add_library(${PROJECT_NAME}::${LIBRARY_NAME} ALIAS ${LIBRARY_NAME})

# Include Directories
target_include_directories(${LIBRARY_NAME} ${PUBLIC_OR_INTERFACE}
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# Properties
set_target_properties(${LIBRARY_NAME} PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED ON
    POSITION_INDEPENDENT_CODE ON
)

# Compile Definitions
target_compile_definitions(${LIBRARY_NAME} ${PUBLIC_OR_INTERFACE}
    MATHBOX_HEADER_ONLY=$<IF:$<BOOL:${MATHBOX_HEADER_ONLY}>,1,0>
)

# Compile Options
if (CMAKE_CXX_COMPILER_VERSION_MAJOR EQUAL 9)
    # Add experimental concepts flag
    target_compile_options(${LIBRARY_NAME} ${PUBLIC_OR_INTERFACE} -fconcepts)
endif()

# Link Libraries
target_link_libraries(${LIBRARY_NAME} ${PUBLIC_OR_INTERFACE}
    Eigen3::Eigen
    cppbox::cppbox
)

# Configure package version and config files
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
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
    NAMESPACE ${PROJECT_NAME}::
    FILE ${PROJECT_NAME}-export.cmake
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

# Build-Tree Export
export(
    EXPORT ${PROJECT_NAME}
    NAMESPACE ${PROJECT_NAME}::
    FILE "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-export.cmake"
)
export(PACKAGE ${PROJECT_NAME})

# Uninstallation
create_uninstall_target()

# Build Documentation
if (BUILD_DOCUMENTATION)
    set(DOXYGEN_VERSION "1.16.0")
    import_doxygen(VERSION ${DOXYGEN_VERSION})
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    configure_file(Doxyfile.in ${DOXYGEN_OUT})
    add_custom_target(${PROJECT_NAME}_documentation ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating documentation with doxygen."
        VERBATIM
    )
    add_dependencies(${PROJECT_NAME}_documentation doxygen)
endif()

# Build Tests
if(BUILD_TESTS)
    add_subdirectory(test)
endif()
