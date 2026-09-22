#!/bin/sh
# ctest gate: generates the icosphere fixture
# then runs the closestPointOnSoup unit-style checks (hand-computed cases
# + seeded brute-force comparison).
#
# Usage: test_closestpoint.sh <test_closestpoint-exe>
set -eu

EXE="$1"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TMP_STL="$(mktemp -u).stl"
trap 'rm -f "$TMP_STL"' EXIT

python3 "$SCRIPT_DIR/scripts/make_stls.py" sphere "$TMP_STL"

"$EXE" "$TMP_STL"
