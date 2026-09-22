#!/bin/sh
# Runs a cut case through ninjaMesher, then checkMesh, and
# asserts "Mesh OK." plus an optional exact kept-cell count and an
# expected total-volume value (checked to `tol`).
#
# Usage: test_cutcase.sh <ninjaMesher-exe> <caseDir> <shape> <stlOut> \
#                         <expectCellsOrEmpty> <expectVolume> <tol>
set -eu

EXE="$1"
CASE_DIR="$2"
SHAPE="$3"
STL_OUT="$4"
EXPECT_CELLS="$5"
EXPECT_VOL="$6"
TOL="$7"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

python3 "$SCRIPT_DIR/scripts/make_stls.py" "$SHAPE" "$STL_OUT"

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

if [ -n "$EXPECT_CELLS" ]; then
    if ! grep -Eq "^[[:space:]]*cells:[[:space:]]+$EXPECT_CELLS[[:space:]]*\$" "$LOG"; then
        echo "FAIL: expected cell count $EXPECT_CELLS not found"
        status=1
    fi
fi

VOL_LINE=$(grep -E "Total volume" "$LOG" || true)
if [ -z "$VOL_LINE" ]; then
    echo "FAIL: 'Total volume' line not found in checkMesh output"
    status=1
else
    ACTUAL_VOL=$(echo "$VOL_LINE" | sed -E 's/.*Total volume = ([0-9.eE+-]+).*/\1/')
    OK=$(awk -v a="$ACTUAL_VOL" -v e="$EXPECT_VOL" -v t="$TOL" \
        'BEGIN { d = a - e; if (d < 0) d = -d; print (d <= t) ? "1" : "0" }')
    if [ "$OK" != "1" ]; then
        echo "FAIL: total volume $ACTUAL_VOL not within $TOL of expected $EXPECT_VOL"
        status=1
    fi
fi

if [ "$status" -ne 0 ]; then
    cat "$LOG"
fi

exit $status
