#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: ./format.sh

Format src/*.f90 and tests/fortran/*.f90 using .fprettify.rc.
Requires fprettify (0.3.7 recommended).

Options:
  -h, --help  Show this help

Environment:
  FPRETTIFY   Formatter executable (default: fprettify)
EOF
}

if [[ $# -eq 1 && ( "$1" == '-h' || "$1" == '--help' ) ]]; then
    usage
    exit 0
fi
if [[ $# -ne 0 ]]; then
    usage >&2
    exit 2
fi

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
formatter=${FPRETTIFY:-fprettify}
if ! command -v -- "$formatter" >/dev/null 2>&1; then
    echo "format.sh: formatter not found: $formatter" >&2
    echo 'Install fprettify with: pipx install fprettify==0.3.7' >&2
    echo 'Or set FPRETTIFY to its executable path.' >&2
    exit 1
fi

exec "$formatter" -c "$script_dir/.fprettify.rc" \
    "$script_dir"/src/*.f90 "$script_dir"/tests/fortran/*.f90
