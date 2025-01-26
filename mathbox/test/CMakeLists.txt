# Enable ctest
enable_testing()

# Build Example
add_executable(${PROJECT_NAME}_tests
    src/covariance.cpp
    src/decompose.cpp
    src/geometry.cpp
    src/integrator.cpp
    src/matrix_operations.cpp
    src/normal.cpp
    src/nsphere.cpp
    src/stiffness.cpp
    src/vector_operations.cpp
)

# Link Libraries
target_link_libraries(${PROJECT_NAME}_tests PUBLIC
    mathbox
    ${mathbox_SYSTEM_LIBRARIES}
    gtest_main
    gtest
    pthread
)

# Properties
set_target_properties(${PROJECT_NAME}_tests PROPERTIES
    CXX_STANDARD 20
)

# Ensure the tests link to the built library, not installed library (this works by setting BUILD_RPATH, and then
# stopping the linker from setting RUNPATH instead of RPATH (RPATH is checked before LD_LIBRARY_PATH, while
# RUNPATH is checked after))
target_link_options(${PROJECT_NAME}_tests PUBLIC
    -Wl,--disable-new-dtags
)

# Compile Options
target_compile_options(${PROJECT_NAME}_tests PRIVATE
    -Wall -Wextra -Wpedantic
)

# Include Directories
target_include_directories(${PROJECT_NAME}_tests PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
)

# System Include Directories
target_include_directories(${PROJECT_NAME}_tests SYSTEM PUBLIC
    ${mathbox_SYSTEM_INCLUDE_DIRS}
)

# Add test
add_test(NAME ${PROJECT_NAME}_tests COMMAND ${PROJECT_NAME}_tests)
