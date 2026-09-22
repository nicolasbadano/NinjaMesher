#!/bin/sh
# ctest gate: generates the same icosphere fixture as
# test_closestpoint.sh, then runs the standalone smoothing checks
# (analytic cases + cost report -- see test_smoothing.cpp).
#
# Usage: test_smoothing.sh <test_smoothing-exe>
set -eu

EXE="$1"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TMP_STL="$(mktemp -u).stl"
trap 'rm -f "$TMP_STL"' EXIT

python3 "$SCRIPT_DIR/scripts/make_stls.py" sphere "$TMP_STL"

"$EXE" "$TMP_STL"
