#!/bin/sh
# refine_box gate: no STL, one level-1
# box region covering the (0,0,0)-(0.5,0.5,0.5) octant of a 4x4x4 domain.
# Asserts checkMesh "Mesh OK.", exactly 120 cells, total volume 1.0
# within 1e-9, and cellLevel counts (56 zeros, 64 ones).
# Usage: test_refinebox.sh <ninjaMesher-exe> <caseDir>
set -eu

EXE="$1"
CASE_DIR="$2"
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

if ! grep -Eq "^[[:space:]]*cells:[[:space:]]+120[[:space:]]*\$" "$LOG"; then
    echo "FAIL: expected cell count 120 not found"
    status=1
fi

VOL_LINE=$(grep -E "Total volume" "$LOG" || true)
if [ -z "$VOL_LINE" ]; then
    echo "FAIL: 'Total volume' line not found in checkMesh output"
    status=1
else
    ACTUAL_VOL=$(echo "$VOL_LINE" | sed -E 's/.*Total volume = ([0-9.eE+-]+).*/\1/')
    OK=$(awk -v a="$ACTUAL_VOL" -v e="1.0" -v t="1e-9" \
        'BEGIN { d = a - e; if (d < 0) d = -d; print (d <= t) ? "1" : "0" }')
    if [ "$OK" != "1" ]; then
        echo "FAIL: total volume $ACTUAL_VOL not within 1e-9 of expected 1.0"
        status=1
    fi
fi

ZEROS=$(grep -c '^0$' "$OUT_DIR/cellLevel" || true)
ONES=$(grep -c '^1$' "$OUT_DIR/cellLevel" || true)
if [ "$ZEROS" != "56" ] || [ "$ONES" != "64" ]; then
    echo "FAIL: expected cellLevel counts 56 zeros / 64 ones, got $ZEROS zeros / $ONES ones"
    status=1
fi

if [ "$status" -ne 0 ]; then
    cat "$LOG"
fi

exit $status
