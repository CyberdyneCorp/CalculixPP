#!/usr/bin/env bash
# Build and install the numerics dependencies for the CalculiX++ solver layer.
#
# NumPP ships a CMake config package (find_package(NumPP)); SciPP does not (it must be
# consumed via add_subdirectory and internally does find_package(NumPP CONFIG REQUIRED)).
# So we install NumPP into a local prefix and hand SciPP a source checkout.
#
# Usage:
#   scripts/bootstrap_deps.sh [--numpp DIR] [--scipp DIR] [--prefix DIR]
# Defaults assume sibling checkouts of NumPP and SciPP next to this repo.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
numpp="${here}/../NumPP"
scipp="${here}/../SciPP"
prefix="${here}/.deps/install"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --numpp)  numpp="$2"; shift 2;;
    --scipp)  scipp="$2"; shift 2;;
    --prefix) prefix="$2"; shift 2;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done

echo ">> NumPP  : ${numpp}"
echo ">> SciPP  : ${scipp}"
echo ">> prefix : ${prefix}"

[[ -f "${numpp}/CMakeLists.txt" ]] || { echo "NumPP not found at ${numpp}" >&2; exit 1; }
[[ -f "${scipp}/CMakeLists.txt" ]] || { echo "SciPP not found at ${scipp}" >&2; exit 1; }

# Build + install NumPP (portable CPU build; no GPU backends).
cmake -S "${numpp}" -B "${numpp}/build-cxpp" -G Ninja \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_INSTALL_PREFIX="${prefix}"
cmake --build "${numpp}/build-cxpp"
cmake --install "${numpp}/build-cxpp"

echo
echo "Done. Configure CalculiX++ with the solver layer via:"
echo "  cmake -S . -B build -G Ninja \\"
echo "     -DCALCULIXPP_WITH_SOLVER=ON \\"
echo "     -DCMAKE_PREFIX_PATH=${prefix} \\"
echo "     -DCALCULIXPP_SCIPP_DIR=${scipp}"
