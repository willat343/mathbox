# mathbox

A general purpose math library based on Eigen. It is currently a header-only library.

## Core C++ Library

### Prerequisites

| **Dependency** | **Version** | **Description** |
|----------------|-------------|-----------------|
| Eigen3 | >= 3.3.7 (< g++-10) or >= 3.3.9 (>= g++-10) | Linear Algebra Package |
| [cppbox](https://github.com/willat343/cppbox) | >= 0.0.2 | C++ Package |

#### Eigen3

If using g++-9, the default Ubuntu 20.04 compiler at time of writing, then Eigen3 can be installed with:
```bash
sudo apt install libeigen-dev
```

If using a more modern g++ compiler, then Eigen3 should be installed from source locally with:
```bash
git clone git@gitlab.com:libeigen/eigen.git
cd eigen
git checkout 3.3.9
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=$HOME/.local ..
make install
```

In this case `-DCMAKE_PREFIX_PATH=$HOME/.local` or `-DEigen3_DIR=$HOME/.local/lib/cmake/Eigen3` must be added to the cmake arguments.

### Installation

It is recommended that you configure with `ccmake` (`sudo apt install cmake-curses-gui`) to see the various options. Otherwise use `cmake` instead of `ccmake` and set flags manually.

```bash
cd mathbox
mkdir build && cd build
ccmake -DCMAKE_PREFIX_PATH=$HOME/.local ..
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

A `package.xml` is supplied to facilitate an isolated installation within a catkin workspace (e.g. for ROS applications).

### Prerequisites

Prerequisites of core C++ library plus the following:

| **Dependency** | **Version** | **Description** |
|----------------|-------------|-----------------|
| catkin | - | catkin build system |

### Installation

Clone or symlink the repository to the workspace's `src` directory, for example:
```bash
ln -s /path/to/mathbox /path/to/catkin_ws/src
cd /path/to/catkin_ws
```

For catkin to search for locally installed packages (e.g. Eigen3) it must be configured:
```bash
catkin config --append-args --cmake-args -DCMAKE_INSTALL_PREFIX=$HOME/.local
```

```bash
catkin build --summary mathbox
```

### Uninstallation

```bash
cd /path/to/catkin_ws
catkin clean mathbox
```

### Usage

To use the package in a downstream project, one should add to their `package.xml`:
```xml
<depend>mathbox</depend>
```
One can then either use the workspace's isolated installation or use the system installation otherwise.
Importing the dependency is then exactly the same as it would be in a non-catkin package as described above (do NOT rely on the `catkin` variables like `catkin_LIBRARIES` and `catkin_INCLUDE_DIRS`).

### Documentation

Documentation must be turned on by setting the `-DBUILD_DOCUMENTATION=ON` cmake argument. This can be done in catkin with:
```bash
catkin config --append-args --cmake-args -DBUILD_DOCUMENTATION=ON
```

### Tests

Tests must be turned on by setting the `-DBUILD_TESTS=ON` cmake argument. This can be done in catkin with:
```bash
catkin config --append-args --cmake-args -DBUILD_TESTS=ON
```
