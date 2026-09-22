#!/bin/sh
# refine_grade gate: grading stressor
# — an 8x8x8 domain with one small level-2 box region at the centre, no
# STL. Asserts checkMesh "Mesh OK.", total volume 1.0 within 1e-9, and
# via --refine-stats: max face-adjacent cellLevel difference == 1.
# Also asserts determinism: --refine-stats output is byte-identical
# across two independent runs.
# Usage: test_refinegrade.sh <ninjaMesher-exe> <caseDir>
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
trap 'rm -f "$LOG" "$STATS1" "$STATS2"' EXIT

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

STATS1=$(mktemp)
STATS2=$(mktemp)
mpirun --bind-to none -np "${NINJA_NP:-1}" "$EXE" --refine-stats "$CASE_DIR" >"$STATS1"
mpirun --bind-to none -np "${NINJA_NP:-1}" "$EXE" --refine-stats "$CASE_DIR" >"$STATS2"

echo "--refine-stats output:"
cat "$STATS1"

if ! grep -q "maxAdjacentLevelDiff = 1" "$STATS1"; then
    echo "FAIL: expected maxAdjacentLevelDiff = 1 (2:1 grading)"
    status=1
fi

if ! diff -q "$STATS1" "$STATS2" >/dev/null; then
    echo "FAIL: --refine-stats output not deterministic across two runs"
    diff "$STATS1" "$STATS2" || true
    status=1
fi

if [ "$status" -ne 0 ]; then
    cat "$LOG"
fi

exit $status
