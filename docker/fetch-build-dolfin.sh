#!/usr/bin/env bash
set -euo pipefail

mkdir -p ./build

cmake -S /work -B build -G Ninja \
    -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_INSTALL_PREFIX=/usr/local \
    -DCMAKE_CXX_COMPILER_LAUNCHER=ccache \
    -DDOLFIN_USE_CGAL=ON \
    -DDOLFIN_ENABLE_TESTS=ON \
    -DDOLFIN_ENABLE_DEMO=ON \
    -DDOLFIN_ENABLE_BENCHMARKS=OFF

cmake --build build -j$(nproc)

cmake --build build -j$(nproc) -- \
    test/unit/cpp/unittests \
    demo/undocumented/multimesh-poisson/cpp/demo_multimesh-poisson \
    demo/undocumented/multimesh-stokes/cpp/demo_multimesh-stokes \
    demo/undocumented/multimesh-3d/cpp/demo_multimesh-3d \
    demo/documented/nonmatching-interpolation/cpp/demo_nonmatching-interpolation

cmake --install build
ldconfig

cd /work/python
PYBIND11_DIR="$(python3 -c 'import pybind11; print(pybind11.get_cmake_dir(), end="")')" \
python3 -m pip install --no-cache-dir --no-deps --no-build-isolation .
