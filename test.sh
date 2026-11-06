#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: ./test.sh [-B BUILD_DIR] [-j JOBS]

Configure a Debug build, build all test targets, and run every CTest case,
including serial/MPI interfaces and I/O compatibility.

Options:
  -B DIR  Build directory (default: ./build-test)
  -j N    Parallel build/test jobs (default: 4)
  -h      Show this help

Environment:
  FC               Fortran compiler (CMake autodetects by default)
  MPIFC            MPI Fortran compiler wrapper (CMake autodetects by default)
  MPIEXEC          MPI launcher (CMake autodetects by default)

Requires CMake >= 3.16, a Fortran compiler, MPI, FFTW3 and Python >= 3.11.
Failure logs: BUILD_DIR/Testing/Temporary/.
EOF
}

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
build_dir="$script_dir/build-test"
jobs=4
while getopts ':B:j:h' option; do
    case "$option" in
        B) build_dir=$OPTARG ;;
        j) jobs=$OPTARG ;;
        h) usage; exit 0 ;;
        :) echo "test.sh: -$OPTARG requires an argument" >&2; exit 2 ;;
        \?) echo "test.sh: unknown option -$OPTARG" >&2; exit 2 ;;
    esac
done
shift "$((OPTIND - 1))"
if [[ $# -ne 0 || ! "$jobs" =~ ^[1-9][0-9]*$ || -z "$build_dir" ]]; then
    usage >&2
    exit 2
fi
# Resolve paths before entering the build tree; invocation is independent of cwd.
if [[ "$build_dir" != /* ]]; then
    build_dir="$PWD/$build_dir"
fi
cmake_args=(-S "$script_dir" -B "$build_dir"
    -DCMAKE_BUILD_TYPE=Debug -DENABLE_UNIT_TESTS=ON
)
if [[ -n ${MPIFC:-} ]]; then
    cmake_args+=("-DMPI_Fortran_COMPILER=$MPIFC")
fi
if [[ -n ${MPIEXEC:-} ]]; then
    cmake_args+=("-DMPIEXEC_EXECUTABLE=$MPIEXEC")
fi
cmake "${cmake_args[@]}"
cmake --build "$build_dir" --config Debug --parallel "$jobs"
# Running in the build directory works with the declared CMake/CTest 3.16 minimum.
cd -- "$build_dir"
ctest -C Debug --output-on-failure --parallel "$jobs"
