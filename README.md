# mathbox

A general purpose math library based on Eigen. It is currently a header-only library.

## Installation

It is recommended that you configure with `ccmake` (`sudo apt install cmake-curses-gui`) to see the various options. Otherwise use `cmake` instead of `ccmake` and set flags manually.

```bash
mkdir build && cd build
ccmake ..
sudo make install
```

## Uninstallation

```bash
cd build
sudo make uninstall
```

## Changelog

* [09.04.2024] migrated functions from various sources
