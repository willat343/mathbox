# mathbox

A general purpose math library based on Eigen. It is currently a header-only library.

## Core C++ Library

### Prerequisites

| **Dependency** | **Version** | **Description** |
|----------------|-------------|-----------------|
| Eigen3 | >= 3.3 | Linear Algebra Package |

### Installation

It is recommended that you configure with `ccmake` (`sudo apt install cmake-curses-gui`) to see the various options. Otherwise use `cmake` instead of `ccmake` and set flags manually.

```bash
cd mathbox
mkdir build && cd build
ccmake ..
make -j
sudo make install
```

### Uninstallation

```bash
cd build
sudo make uninstall
```

### Usage

Import the package into your project and add the dependency to your target `<target>` with:
```cmake
find_package(mathbox REQUIRED)
target_link_libraries(<target> <PUBLIC|INTERFACE|PRIVATE> ${mathbox_LIBRARIES})
target_include_directories(<target> SYSTEM <PUBLIC|INTERFACE|PRIVATE> ${mathbox_INCLUDE_DIRS})
```

For more information, see [mathbox/README.md](mathbox/README.md) and documentation.

### Documentation

Documentation must be turned on by setting the `-DBUILD_DOCUMENTATION=ON` cmake argument.

To view the HTML documentation, open the `build/docs/html/index.html` file.

To view the LaTeX documentation, build it with:
```bash
cd build/docs/latex
make
```
Then open the file `refman.pdf`.

### Tests

Tests must be turned on by setting the `-DBUILD_TESTS=ON` cmake argument.

```bash
cd build
cmake -DBUILD_TESTS=ON ..
make -j
```

They can then be run with `ctest`:
```bash
ctest --test-dir test
```

For more explicit output, the test executables can be run directly from the build directory.

## Catkin Support

A catkin wrapper is available to facilitate an isolated installation within a catkin workspace (e.g. for ROS applications).

### Prerequisites

| **Dependency** | **Version** | **Description** |
|----------------|-------------|-----------------|
| Eigen3 | >= 3.3 | Linear Algebra Package |

### Installation

Clone or symlink the repository to the workspace's `src` directory, for example:
```bash
ln -s /path/to/mathbox /path/to/catkin_ws/src
```

```bash
cd /path/to/catkin_ws
catkin build mathbox_catkin
```

### Uninstallation

```bash
cd /path/to/catkin_ws
catkin clean mathbox_catkin
```

### Usage

To use the package in a downstream project, one should add to their `package.xml`:
```xml
<exec_depend>mathbox_catkin</exec_depend>
```
One can then either use the workspace's isolated installation if the catkin wrapper exists in the workspace, or use the system installation otherwise.
Importing the dependency is then exactly the same as it would be in a non-catkin package as described above (do NOT rely on the `catkin` variables like `catkin_LIBRARIES` and `catkin_INCLUDE_DIRS`).

### Documentation

Documentation must be turned on by setting the `-DBUILD_DOCUMENTATION=ON` cmake argument. This can be done in catkin with:
```bash
catkin config --cmake-args -DBUILD_DOCUMENTATION=ON
```

### Tests

Tests must be turned on by setting the `-DBUILD_TESTS=ON` cmake argument. This can be done in catkin with:
```bash
catkin config --cmake-args -DBUILD_TESTS=ON
```
