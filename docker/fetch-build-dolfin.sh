#!/usr/bin/env bash
set -euo pipefail

dockerfile="${DOLFIN_DOCKERFILE:-/work/docker/Dockerfile}"
src_dir="${DOLFIN_SRC_DIR:-/work}"
build_dir="${DOLFIN_BUILD_DIR:-${src_dir}/build}"
install_prefix="${DOLFIN_INSTALL_PREFIX:-/usr/local}"
build_type="${DOLFIN_BUILD_TYPE:-RelWithDebInfo}"
use_cgal="${DOLFIN_USE_CGAL:-OFF}"
jobs="${DOLFIN_JOBS:-$(nproc)}"

if [[ ! -f "${dockerfile}" ]]; then
    echo "Expected Dockerfile at ${dockerfile}" >&2
    exit 1
fi

if [[ ! -f "${src_dir}/CMakeLists.txt" ]]; then
    echo "Expected DOLFIN source tree at ${src_dir} (missing CMakeLists.txt)" >&2
    exit 1
fi

mkdir -p "${build_dir}"

if ! cmake -S "${src_dir}" -B "${build_dir}" -G Ninja \
    -DCMAKE_BUILD_TYPE="${build_type}" \
    -DCMAKE_INSTALL_PREFIX="${install_prefix}" \
    -DCMAKE_CXX_COMPILER_LAUNCHER=ccache \
    -DDOLFIN_USE_CGAL="${use_cgal}" \
    -DDOLFIN_ENABLE_TESTS=ON \
    -DDOLFIN_ENABLE_DEMO=ON \
    -DDOLFIN_ENABLE_BENCHMARKS=OFF
then
    echo "---- CMakeError.log ----"
    cat "${build_dir}/CMakeFiles/CMakeError.log" || true
    echo "---- CMakeOutput.log ----"
    cat "${build_dir}/CMakeFiles/CMakeOutput.log" || true
    exit 1
fi

cmake --build "${build_dir}" -j"${jobs}"

cmake --build "${build_dir}" -j"${jobs}" -- \
    test/unit/cpp/unittests \
    demo/undocumented/multimesh-poisson/cpp/demo_multimesh-poisson \
    demo/undocumented/multimesh-stokes/cpp/demo_multimesh-stokes \
    demo/undocumented/multimesh-3d/cpp/demo_multimesh-3d

cmake --install "${build_dir}"
ldconfig

cd "${src_dir}/python"
PYBIND11_DIR="$(python3 -c 'import pybind11; print(pybind11.get_cmake_dir(), end="")')" \
python3 -m pip install --no-cache-dir --no-deps --no-build-isolation .
