#!/usr/bin/env bash
set -euxo pipefail

# Conda provides CPU_COUNT
JOBS=${CPU_COUNT}
export CMAKE_BUILD_PARALLEL_LEVEL="${JOBS}"

# Point to the specific Python being built against
export Python3_ROOT_DIR="${PREFIX}"

# Resolve NumPy include path from the host environment
NUMPY_INCLUDE="$(${PYTHON} -c 'import numpy; print(numpy.get_include())')"

# Use _conda_build to avoid collisions with a local dev build/ directory
cmake -S . -B _conda_build -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
    -DCMAKE_PREFIX_PATH="${PREFIX}" \
    -DRDK_BUILD_PYTHON_WRAPPERS=ON \
    -DRDK_BUILD_TESTS=OFF \
    -DRDK_BUILD_CPP_TESTS=OFF \
    -DRDK_BUILD_PYTHON_TESTS=OFF \
    -DRDK_BUILD_EXAMPLES=OFF \
    -DPYTHON_EXECUTABLE="${PYTHON}" \
    -DPython3_EXECUTABLE="${PYTHON}" \
    -DPython3_FIND_STRATEGY=LOCATION \
    -DPython3_INCLUDE_DIR="${PREFIX}/include/python${PY_VER}" \
    -DPython3_NumPy_INCLUDE_DIR="${NUMPY_INCLUDE}" \
    -DEIGEN3_INCLUDE_DIR="${PREFIX}/include/eigen3" \
    -DRDK_INSTALL_INTREE=OFF \
    -DRDK_INSTALL_DEV_COMPONENT=ON \
    -DBoost_NO_BOOST_CMAKE=ON

# Build and Install
cmake --build _conda_build -j"${JOBS}" -v
cmake --install _conda_build