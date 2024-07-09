# Enable ctest
enable_testing()

# Build Example
add_executable(mathbox_tests
    src/decompose.cpp
    src/integrator.cpp
    src/normal.cpp
    src/nsphere.cpp
    src/stiffness.cpp
)

# Properties
set_target_properties(mathbox_tests PROPERTIES
    CXX_STANDARD 20
)

# Compile Options
target_compile_options(mathbox_tests PRIVATE
    -Wall -Wextra -Wpedantic
)

# Include Directories
target_include_directories(mathbox_tests PUBLIC
    ../include
    include
)

# System Include Directories
# target_include_directories(mathbox_tests SYSTEM PUBLIC)

# Link Libraries
target_link_libraries(mathbox_tests PUBLIC
    ${MATHBOX_SYSTEM_LIBRARIES}
    gtest_main
    gtest
    pthread
)

# Add test
add_test(NAME mathbox_tests COMMAND mathbox_tests)
