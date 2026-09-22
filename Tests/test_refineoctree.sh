#!/bin/sh
# refine_octree gate: a single base
# cell, two nested box regions, forcing true per-octant refinement
# (sibling octants at mixed levels within ONE base cell). Hand-derivable:
# 7 level-1 leaves + 8 level-2 leaves = 15 cells, volume 1.0 exactly.
# Asserts checkMesh "Mesh OK.", exactly 15 cells, total volume 1.0
# within 1e-9, cellLevel counts (7 ones, 8 twos), and
# maxAdjacentLevelDiff==1 via --refine-stats.
# Usage: test_refineoctree.sh <ninjaMesher-exe> <caseDir>
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

if ! grep -Eq "^[[:space:]]*cells:[[:space:]]+15[[:space:]]*\$" "$LOG"; then
    echo "FAIL: expected cell count 15 not found"
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

ONES=$(grep -c '^1$' "$OUT_DIR/cellLevel" || true)
TWOS=$(grep -c '^2$' "$OUT_DIR/cellLevel" || true)
if [ "$ONES" != "7" ] || [ "$TWOS" != "8" ]; then
    echo "FAIL: expected cellLevel counts 7 ones / 8 twos, got $ONES ones / $TWOS twos"
    status=1
fi

STATS=$(mpirun --bind-to none -np "${NINJA_NP:-1}" "$EXE" --refine-stats "$CASE_DIR")
echo "--refine-stats output:"
echo "$STATS"
if ! echo "$STATS" | grep -q "maxAdjacentLevelDiff = 1"; then
    echo "FAIL: expected maxAdjacentLevelDiff = 1 (2:1 grading)"
    status=1
fi

if [ "$status" -ne 0 ]; then
    cat "$LOG"
fi

exit $status
