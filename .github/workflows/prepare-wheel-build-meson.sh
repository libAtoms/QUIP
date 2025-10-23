#!/bin/bash
set -e

echo "prepare-wheel-build-meson.sh: Building QUIP with Meson"
echo "Current directory: $(pwd)"
echo "Script location: $0"

# Find repository root relative to this script
# Script is in .github/workflows/, so repo root is ../..
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO_ROOT=$(cd "${SCRIPT_DIR}/../.." && pwd)

echo "Repository root: ${REPO_ROOT}"
cd "${REPO_ROOT}"

# Check if meson.build exists
if [ ! -f "meson.build" ]; then
    echo "ERROR: meson.build not found in ${REPO_ROOT}"
    ls -la
    exit 1
fi

# Build QUIP libraries if not already built
if [ ! -d "builddir" ]; then
    echo "Setting up Meson build..."
    meson setup builddir
    echo "Compiling QUIP..."
    meson compile -C builddir
else
    echo "Build directory already exists, skipping build"
fi

echo "QUIP build complete"
ls -la builddir/src/Programs/ || true
