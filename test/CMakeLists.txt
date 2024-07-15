# Enable ctest
enable_testing()

# Build Example
add_executable(${PROJECT_NAME}_tests
    src/decompose.cpp
    src/geometry.cpp
    src/integrator.cpp
    src/matrix_operations.cpp
    src/normal.cpp
    src/nsphere.cpp
    src/stiffness.cpp
)

# Properties
set_target_properties(${PROJECT_NAME}_tests PROPERTIES
    CXX_STANDARD 20
)

# Compile Options
target_compile_options(${PROJECT_NAME}_tests PRIVATE
    -Wall -Wextra -Wpedantic
)

# Include Directories
target_include_directories(${PROJECT_NAME}_tests PUBLIC
    ../include
    include
)

# System Include Directories
# target_include_directories(${PROJECT_NAME}_tests SYSTEM PUBLIC)

# Link Libraries
target_link_libraries(${PROJECT_NAME}_tests PUBLIC
    ${MATHBOX_SYSTEM_LIBRARIES}
    gtest_main
    gtest
    pthread
)

# Add test
add_test(NAME ${PROJECT_NAME}_tests COMMAND ${PROJECT_NAME}_tests)
