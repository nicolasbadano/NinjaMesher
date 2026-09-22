#!/bin/sh
# Gate: a full-cross-section plate splits the domain
# into two fluid pockets; `locationInMesh` picks one, and the flood-fill
# connectivity drop must remove the other. Asserts:
#   - stdout reports disconnectedCellsDropped > 0, discardedComponents == 1
#   - checkMesh "Mesh OK." + "Number of regions: 1"
#   - check_integrity.py clean
#   - total volume within [volMin, volMax] (the reachable region only)
#
# Usage: test_cutdisconnected.sh <ninjaMesher-exe> <caseDir> <volMin> <volMax>
set -eu

EXE="$1"; CASE_DIR="$2"; VOL_MIN="$3"; VOL_MAX="$4"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

python3 "$SCRIPT_DIR/scripts/make_stls.py" box_at 0.48 0 0 0.52 1 1 "$CASE_DIR/plate.stl"

OUT_DIR="$CASE_DIR/constant/polyMesh"
rm -rf "$OUT_DIR"
RUN_LOG=$(mktemp)
mpirun --bind-to none -np "${NINJA_NP:-1}" "$EXE" "$CASE_DIR" | tee "$RUN_LOG"

status=0

DROPPED=$(grep -E "^disconnectedCellsDropped" "$RUN_LOG" | sed -E 's/.*= ([0-9]+).*/\1/')
COMPONENTS=$(grep -E "^discardedComponents" "$RUN_LOG" | sed -E 's/.*= ([0-9]+).*/\1/')
echo "measured disconnectedCellsDropped = $DROPPED, discardedComponents = $COMPONENTS"
if [ -z "$DROPPED" ] || [ "$DROPPED" -le 0 ]; then
    echo "FAIL: expected disconnectedCellsDropped > 0"
    status=1
fi
if [ "$COMPONENTS" != "1" ]; then
    echo "FAIL: expected discardedComponents == 1, got $COMPONENTS"
    status=1
fi
rm -f "$RUN_LOG"

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

if ! grep -q "Mesh OK\." "$LOG"; then
    echo "FAIL: 'Mesh OK.' verdict not found in checkMesh output"
    status=1
fi

if ! grep -q "Number of regions: 1" "$LOG"; then
    echo "FAIL: expected 'Number of regions: 1' (the dropped pocket should no longer exist)"
    status=1
fi

if ! python3 "$SCRIPT_DIR/../Benchmarks/check_integrity.py" "$OUT_DIR"; then
    echo "FAIL: check_integrity.py reported a structural pathology"
    status=1
fi

VOL_LINE=$(grep -E "Total volume" "$LOG" || true)
if [ -z "$VOL_LINE" ]; then
    echo "FAIL: 'Total volume' line not found in checkMesh output"
    status=1
else
    ACTUAL_VOL=$(echo "$VOL_LINE" | sed -E 's/.*Total volume = ([0-9.eE+-]+).*/\1/')
    echo "measured volume = $ACTUAL_VOL"
    OK=$(awk -v a="$ACTUAL_VOL" -v lo="$VOL_MIN" -v hi="$VOL_MAX" \
        'BEGIN { print (a >= lo && a <= hi) ? "1" : "0" }')
    if [ "$OK" != "1" ]; then
        echo "FAIL: volume $ACTUAL_VOL not within [$VOL_MIN, $VOL_MAX]"
        status=1
    fi
fi

if [ "$status" -ne 0 ]; then
    cat "$LOG"
fi

exit $status
