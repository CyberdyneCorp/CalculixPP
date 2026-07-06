# CalculiX++ developer tasks — run `just` (or `just --list`) to see everything.
# Requires: cmake >= 3.24, a C++20 compiler, and a NumPP install + a SciPP checkout
# (see `just bootstrap`). No GPU toolkit is ever required — the solver runs on the
# portable CPU reference backend. `just gpu-detect` reports what your host has.
#
# Layout the recipes assume (override the variables below):
#   ../NumPP   — built + installed into .deps/install by `just bootstrap`
#   ../SciPP   — source checkout, consumed via add_subdirectory (has no install rules)

build_dir   := "build"
deps_prefix := ".deps/install"
numpp_src   := "../NumPP"
scipp_src   := "../SciPP"

# Show the available recipes (default).
default:
    @just --list

# Build + install NumPP into .deps/install and point at a SciPP checkout (override: just bootstrap ../NumPP ../SciPP).
bootstrap numpp=numpp_src scipp=scipp_src:
    ./scripts/bootstrap_deps.sh --numpp {{numpp}} --scipp {{scipp}}

# Configure a Release build with the solver layer; pass extra cmake flags (e.g. -DCALCULIXPP_BUILD_PYTHON=ON).
configure *FLAGS:
    cmake -S . -B {{build_dir}} -G Ninja \
      -DCMAKE_BUILD_TYPE=Release \
      -DCALCULIXPP_WITH_SOLVER=ON \
      -DCMAKE_PREFIX_PATH="{{justfile_directory()}}/{{deps_prefix}}" \
      -DCALCULIXPP_SCIPP_DIR={{scipp_src}} {{FLAGS}}

# Compile everything already configured into {{build_dir}}.
compile:
    cmake --build {{build_dir}} -j

# Configure (solver, CPU) + compile.
build: configure compile

# Configure + build + run the full CTest suite (CPU backend, no GPU toolkit).
test: build
    ctest --test-dir {{build_dir}} --output-on-failure
alias ctest := test

# Build the pybind11 module and run the Python regression suite (decks diffed vs CalculiX).
python:
    @just configure -DCALCULIXPP_BUILD_PYTHON=ON
    just compile
    ctest --test-dir {{build_dir}} -R python_regression --output-on-failure

# Build + test the dependency-free core alone (no NumPP/SciPP) — the mobile-toolchain path.
core:
    cmake -S . -B build-core -G Ninja -DCMAKE_BUILD_TYPE=Release -DCALCULIXPP_WITH_SOLVER=OFF
    cmake --build build-core -j
    ctest --test-dir build-core --output-on-failure

# Configure + build + test with debug symbols and assertions (separate dir).
debug:
    cmake -S . -B build-debug -G Ninja \
      -DCMAKE_BUILD_TYPE=Debug \
      -DCALCULIXPP_WITH_SOLVER=ON \
      -DCMAKE_PREFIX_PATH="{{justfile_directory()}}/{{deps_prefix}}" \
      -DCALCULIXPP_SCIPP_DIR={{scipp_src}}
    cmake --build build-debug -j
    ctest --test-dir build-debug --output-on-failure

# Detect usable GPU backends (CUDA / OpenCL / Metal) on this host — informational; the CPU backend is always used.
gpu-detect:
    #!/usr/bin/env bash
    set -uo pipefail
    os=$(uname -s); arch=$(uname -m)
    case "$os" in
      MINGW*|MSYS*|CYGWIN*) plat=Windows ;;
      Darwin)               plat=macOS   ;;
      Linux)                plat=Linux   ;;
      *)                    plat="$os"   ;;
    esac
    echo "Host: $plat ($os $arch)"
    echo
    cuda=no; opencl=no; metal=no

    # --- CUDA (NVIDIA) — driver at runtime, nvcc to build ----------------------
    if command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi -L >/dev/null 2>&1; then
      gpu=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)
      drv=$(nvidia-smi --query-gpu=driver_version --format=csv,noheader 2>/dev/null | head -1)
      if command -v nvcc >/dev/null 2>&1; then
        nvcc=$(nvcc --version | grep -oE 'release [0-9.]+' | head -1)
      else
        nvcc="nvcc not on PATH — install the CUDA toolkit to build"
      fi
      echo "✅ CUDA   : ${gpu:-NVIDIA GPU} (driver ${drv:-?}, ${nvcc})"
      cuda=yes
    else
      echo "❌ CUDA   : no NVIDIA driver (nvidia-smi) detected"
    fi

    # --- OpenCL — clinfo if present, else probe the ICD loader / framework -----
    if command -v clinfo >/dev/null 2>&1 && [ "$(clinfo -l 2>/dev/null | grep -ciE 'device|platform')" -gt 0 ]; then
      dev=$(clinfo -l 2>/dev/null | grep -iE 'device' | head -1 | sed 's/^[[:space:]]*//')
      echo "✅ OpenCL : ${dev:-device present} (clinfo)"; opencl=yes
    elif [ "$plat" = macOS ] && [ -d /System/Library/Frameworks/OpenCL.framework ]; then
      echo "✅ OpenCL : Apple OpenCL.framework present (install clinfo for details)"; opencl=yes
    elif [ "$plat" = Windows ] && { [ -f /c/Windows/System32/OpenCL.dll ] || reg query "HKLM\\SOFTWARE\\Khronos\\OpenCL\\Vendors" >/dev/null 2>&1; }; then
      echo "✅ OpenCL : Windows OpenCL ICD present (install clinfo to enumerate)"; opencl=yes
    elif ls /etc/OpenCL/vendors/*.icd >/dev/null 2>&1 || ldconfig -p 2>/dev/null | grep -q libOpenCL; then
      echo "✅ OpenCL : ICD loader present (install clinfo to enumerate devices)"; opencl=yes
    else
      echo "❌ OpenCL : no ICD loader / OpenCL runtime detected"
    fi

    # --- Metal (Apple GPU) — Apple platforms only ------------------------------
    if [ "$plat" = macOS ] && [ -d /System/Library/Frameworks/Metal.framework ]; then
      dev=$(system_profiler SPDisplaysDataType 2>/dev/null | grep -iE 'Chipset Model' | head -1 | sed 's/^[[:space:]]*//')
      echo "✅ Metal  : ${dev:-Metal framework present}"; metal=yes
    else
      echo "❌ Metal  : Apple platforms only"
    fi

    echo
    if [ "$cuda" = yes ] || [ "$metal" = yes ] || [ "$opencl" = yes ]; then
      echo "A GPU is present. CalculiX++ ships only the CPU reference backend today —"
      echo "GPU acceleration is an additive, later-phase capability. Run: just test"
    else
      echo "CPU only — the portable reference backend, always available. Run: just test"
    fi

# Validate all OpenSpec specs and changes.
spec:
    openspec validate --all --strict

# Full local CI: CPU tests + spec validation.
ci: test spec

# Remove all build directories and the local dependency prefix.
clean:
    rm -rf {{build_dir}} build-core build-debug .deps
