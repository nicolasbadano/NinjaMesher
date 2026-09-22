#!/bin/sh
# Offset-cut ctest gate: generates the STL fixture,
# runs ninjaMesher (layers block present in the case dict), checkMesh,
# check_integrity.py, an optional volume window check, and prints
# --cut-stats (so multiRootEdges is captured in the ctest log).
#
# Usage: test_offsetcase.sh <ninjaMesher-exe> <caseDir> <stlgen-args...> \
#            -- <stlOut> <volMin> <volMax>
# <stlgen-args...> is passed verbatim to make_stls.py (shape + params);
# <stlOut> is where the STL lands (must match the case dict's geometry
# entry); <volMin>/<volMax> bound the checkMesh "Total volume" (use
# "" / "" to skip the volume check entirely, e.g. when the caller wants
# to inspect stdout instead).
set -eu

EXE="$1"; shift
CASE_DIR="$1"; shift

STLGEN_ARGS=""
while [ "$1" != "--" ]; do
    STLGEN_ARGS="$STLGEN_ARGS $1"
    shift
done
shift # consume --
STL_OUT="$1"; VOL_MIN="$2"; VOL_MAX="$3"; WALL_MIN="${4:-}"; WALL_MAX="${5:-}"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# shellcheck disable=SC2086
python3 "$SCRIPT_DIR/scripts/make_stls.py" $STLGEN_ARGS "$STL_OUT"

OUT_DIR="$CASE_DIR/constant/polyMesh"
rm -rf "$OUT_DIR"
RUN_LOG=$(mktemp)
mpirun --bind-to none -np "${NINJA_NP:-1}" "$EXE" "$CASE_DIR" | tee "$RUN_LOG"

echo "--- --cut-stats ---"
mpirun --bind-to none -np "${NINJA_NP:-1}" "$EXE" --cut-stats "$CASE_DIR"

status=0
if [ -n "$WALL_MIN" ] && [ -n "$WALL_MAX" ]; then
    WALL_LINE=$(grep -E "^wallArea" "$RUN_LOG" || true)
    if [ -z "$WALL_LINE" ]; then
        echo "FAIL: 'wallArea' line not found in ninjaMesher output"
        status=1
    else
        ACTUAL_WALL=$(echo "$WALL_LINE" | sed -E 's/.*= ([0-9.eE+-]+).*/\1/')
        echo "measured wallArea = $ACTUAL_WALL"
        OK=$(awk -v a="$ACTUAL_WALL" -v lo="$WALL_MIN" -v hi="$WALL_MAX" \
            'BEGIN { print (a >= lo && a <= hi) ? "1" : "0" }')
        if [ "$OK" != "1" ]; then
            echo "FAIL: wallArea $ACTUAL_WALL not within [$WALL_MIN, $WALL_MAX]"
            status=1
        fi
    fi
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
    if [ -n "$VOL_MIN" ] && [ -n "$VOL_MAX" ]; then
        OK=$(awk -v a="$ACTUAL_VOL" -v lo="$VOL_MIN" -v hi="$VOL_MAX" \
            'BEGIN { print (a >= lo && a <= hi) ? "1" : "0" }')
        if [ "$OK" != "1" ]; then
            echo "FAIL: volume $ACTUAL_VOL not within [$VOL_MIN, $VOL_MAX]"
            status=1
        fi
    fi
fi

if [ "$status" -ne 0 ]; then
    cat "$LOG"
fi

exit $status
