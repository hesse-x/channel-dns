#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: ./build.sh [-j JOBS]

Build the Fortran / pure MPI solver in Release mode (-O3) in ./build-channel-dns.
Executable: build-channel-dns/bin/channel_dns

Options:
  -j N  Parallel build jobs (default: 4)
  -h    Show this help

Environment:
  FC       Fortran compiler (CMake autodetects by default)
  MPIFC    MPI Fortran compiler wrapper (CMake autodetects by default)
  MPIEXEC  MPI launcher (CMake autodetects by default)

Tests are disabled. Run ./test.sh for all tests.
EOF
}

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
build_dir="$script_dir/build-channel-dns"
jobs=4
while getopts ':j:h' option; do
    case "$option" in
        j) jobs=$OPTARG ;;
        h) usage; exit 0 ;;
        :) echo "build.sh: -$OPTARG requires an argument" >&2; exit 2 ;;
        \?) echo "build.sh: unknown option -$OPTARG" >&2; exit 2 ;;
    esac
done
shift "$((OPTIND - 1))"
if [[ $# -ne 0 || ! "$jobs" =~ ^[1-9][0-9]*$ ]]; then
    usage >&2
    exit 2
fi
cmake_args=(-S "$script_dir" -B "$build_dir"
    -DCMAKE_BUILD_TYPE=Release -DENABLE_UNIT_TESTS=OFF
    "-DCMAKE_RUNTIME_OUTPUT_DIRECTORY=$build_dir/bin"
    "-DCMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE=$build_dir/bin")
if [[ -n ${MPIFC:-} ]]; then
    cmake_args+=("-DMPI_Fortran_COMPILER=$MPIFC")
fi
if [[ -n ${MPIEXEC:-} ]]; then
    cmake_args+=("-DMPIEXEC_EXECUTABLE=$MPIEXEC")
fi
cmake "${cmake_args[@]}"
cmake --build "$build_dir" --config Release --target channel_dns --parallel "$jobs"
printf '\nExecutable: %s/bin/channel_dns\n' "$build_dir"
