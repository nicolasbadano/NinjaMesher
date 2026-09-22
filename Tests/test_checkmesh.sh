#!/bin/sh
# Runs hex_10x8x6 through ninjaMesher, then OpenFOAM's checkMesh, and
# asserts the literal "Mesh OK." verdict plus the expected counts.
# Usage: test_checkmesh.sh <ninjaMesher-exe> <repoRoot>
set -eu

EXE="$1"
ROOT="$2"
CASE_DIR="$ROOT/Tests/cases/hex_10x8x6"
OUT_DIR="$CASE_DIR/constant/polyMesh"

rm -rf "$OUT_DIR"
mpirun --bind-to none -np "${NINJA_NP:-1}" "$EXE" "$CASE_DIR"

if ! command -v checkMesh >/dev/null 2>&1; then
    echo "FAIL: checkMesh not found on PATH (source OpenFOAM environment first)"
    exit 1
fi

LOG=$(mktemp)
trap 'rm -f "$LOG"' EXIT

if ! checkMesh -case "$CASE_DIR" >"$LOG" 2>&1; then
    echo "FAIL: checkMesh exited non-zero"
    cat "$LOG"
    exit 1
fi

status=0

if ! grep -q "Mesh OK\." "$LOG"; then
    echo "FAIL: 'Mesh OK.' verdict not found in checkMesh output"
    status=1
fi

# Patterns match checkMesh's "Mesh stats" block verbatim, e.g.:
#     points:           693
#     internal faces:   1252
#     cells:            480
if ! grep -Eq "^\s*cells:\s+480\s*$" "$LOG"; then
    echo "FAIL: expected cell count 480 not found"
    status=1
fi

if ! grep -Eq "^\s*points:\s+693\s*$" "$LOG"; then
    echo "FAIL: expected point count 693 not found"
    status=1
fi

if ! grep -Eq "^\s*internal faces:\s+1252\s*$" "$LOG"; then
    echo "FAIL: expected internal face count 1252 not found"
    status=1
fi

if [ "$status" -ne 0 ]; then
    cat "$LOG"
fi

exit $status
