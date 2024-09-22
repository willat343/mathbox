# mathbox

A general purpose math library based on Eigen. It is currently a header-only library.

## Core C++ Library

### Installation

It is recommended that you configure with `ccmake` (`sudo apt install cmake-curses-gui`) to see the various options. Otherwise use `cmake` instead of `ccmake` and set flags manually.

```bash
cd mathbox
mkdir build && cd build
ccmake ..
sudo make install
```

### Uninstallation

```bash
cd build
sudo make uninstall
```

## Catkin Support

A catkin wrapper is available to facilitate easy integration with the catkin build system (e.g. for ROS applications). To use mathbox with catkin, simply clone or symlink the mathbox repository to the catkin workspace `src` directory:
```bash
ln -s /path/to/mathbox /path/to/catkin_ws/src
```

Your project can then depend on `mathbox_catkin` in the `package.xml` and `mathbox_catkin` can be added to `COMPONENTS`, e.g. `find_package(catkin REQUIRED COMPONENTS mathbox_catkin)`.

## Changelog

* [09.04.2024] migrated functions from various sources
