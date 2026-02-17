#!/usr/bin/env bash
set -euxo pipefail

# Conda provides CPU_COUNT
JOBS=${CPU_COUNT}
export CMAKE_BUILD_PARALLEL_LEVEL="${JOBS}"

# Point to the specific Python being built against
export Python3_ROOT_DIR="${PREFIX}"

cmake -S . -B build -G Ninja \
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
    -DRDK_INSTALL_INTREE=OFF \
    -DRDK_INSTALL_DEV_COMPONENT=ON \
    -DBoost_NO_BOOST_CMAKE=ON

# Build and Install
cmake --build build -j"${JOBS}" -v
cmake --install build