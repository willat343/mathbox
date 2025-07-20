# mathbox

A general purpose math library based on Eigen. It is currently a header-only library.

## Core C++ Library

| **Dependency** | **Version** | **Description** |
|----------------|-------------|-----------------|
| CMake | >= 3.21 | CMake Build Tool |
| [cmakebox](https://github.com/willat343/cmakebox) | >= 0.0.1 | CMake Functions and Utilities |
| [cppbox](https://github.com/willat343/cppbox) | >= 0.1.0 | C++ Package |
| Eigen3 | >= 3.4.0 | Linear Algebra Package |

There are several ways to include `cppbox` within your project:
- [Preferred] Via `FetchContent` allowing `cppbox` to be built as a submodule.
- Via `find_package`, requiring `cppbox` to be installed to the system, locally, or to a catkin workspace.

## Include via FetchContent

It is recommended to leverage the functionality of [cmakebox](https://github.com/willat343/cmakebox) by including the following lines in the `CMakeLists.txt` (replace `X.Y.Z` with version):
```CMake
set(CMAKEBOX_VERSION "0.0.1")
FetchContent_Declare(
    cmakebox
    GIT_REPOSITORY git@github.com:willat343/cmakebox.git
    GIT_TAG v${CMAKEBOX_VERSION}
)
FetchContent_MakeAvailable(cmakebox)
list(APPEND CMAKE_MODULE_PATH "${cmakebox_SOURCE_DIR}/cmake")
include(CMakeBox)

set(MATHBOX_VERSION "X.Y.Z")
import_dependency(
    mathbox
    TARGET mathbox::mathbox
    VERSION ${MATHBOX_VERSION}
    USE_SYSTEM_REQUIRED_VERSION ${MATHBOX_VERSION}
    GIT_REPOSITORY git@github.com:willat343/mathbox
    GIT_TAG v${MATHBOX_VERSION}
)
```

Without relying on [cmakebox](https://github.com/willat343/cmakebox), this can be achieved with (replace `X.Y.Z` with version):
```CMake
set(MATHBOX_VERSION "X.Y.Z")
FetchContent_Declare(
    mathbox
    GIT_REPOSITORY git@github.com:willat343/mathbox
    GIT_TAG        v${MATHBOX_VERSION}
)
FetchContent_MakeAvailable(mathbox)
```

## Include via Install

### Clone

```bash
git clone git@github.com:willat343/mathbox.git
cd mathbox
```

### Configure

For system install:
```bash
cmake -S . -B build
```

For local install:
```bash
cmake -S . -B build -DCMAKE_INSTALL_DIR=$HOME/.local
```

### Build

```bash
cmake --build build -j
```

### Install

```bash
sudo cmake --build build --target install
```

### Include

Include with the following lines in the `CMakeLists.txt`:
```CMake
find_package(mathbox REQUIRED)
target_link_libraries(<target> PUBLIC mathbox::mathbox)
```

### Uninstall

```bash
sudo cmake --build build --target uninstall
```

## Include in Catkin Workspace

A `package.xml` is supplied to facilitate an isolated installation within a catkin workspace (e.g. for ROS applications).

### Clone

```bash
cd /path/to/catkin_ws/src
git clone git@github.com:willat343/mathbox.git
```

### Build

```bash
cd /path/to/catkin_ws
catkin build mathbox
```

### Include

To use the package in a downstream project, one should add to their `package.xml`:
```xml
<depend>mathbox</depend>
```

One can then include `mathbox` package by includeing in the `CMakeLists.txt`:
```CMake
find_package(mathbox REQUIRED)
target_link_libraries(<target> PUBLIC mathbox::mathbox)
```

### Clean

```bash
cd /path/to/catkin_ws
catkin clean mathbox
```

## Documentation

Documentation must be turned on by setting the `-DBUILD_DOCUMENTATION=ON` cmake argument.

To view the HTML documentation, open the `build/docs/html/index.html` file.

To view the LaTeX documentation, build it with:
```bash
cd build/docs/latex
make
```
Then open the file `refman.pdf`.

## Tests

Tests must be turned on by setting the `-DBUILD_TESTS=ON` cmake argument.

They can then be run with `ctest`:
```bash
ctest --test-dir test
```

For more explicit output, the test executables can be run directly from the build directory.
