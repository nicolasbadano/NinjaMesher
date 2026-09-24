#!/bin/sh
# Gate for offset_forcedrop: the drop-locally seam-closure path,
# exercised deterministically via the layersDebug{forceDropSphere} test
# scaffold (a polar cap of the offset_sphere case's layer-top faces is
# force-reverted regardless of geometric validity).
#
# Asserts: droppedFaces > 0 (the scaffold actually fired), wallArea
# strictly between the landed wall (~0.755, all-faces-landed) and the
# offset surface (4*pi*0.3^2 = 1.131, nothing-landed) -- the dropped
# cap stays at the offset surface while the rest lands -- plus
# "Mesh OK." and structural integrity (seam quads close the reverted
# region watertight).
#
# Usage: test_offsetforcedrop.sh <ninjaMesher-exe> <caseDir>
set -eu

EXE="$1"
CASE_DIR="$2"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

python3 "$SCRIPT_DIR/scripts/make_stls.py" sphere_at 0.5 0.5 0.5 0.25 3 "$CASE_DIR/sphere.stl"

OUT_DIR="$CASE_DIR/constant/polyMesh"
rm -rf "$OUT_DIR"
RUN_LOG=$(mktemp)
mpirun --bind-to none -np "${NINJA_NP:-1}" "$EXE" "$CASE_DIR" | tee "$RUN_LOG"

echo "--- --cut-stats ---"
mpirun --bind-to none -np "${NINJA_NP:-1}" "$EXE" --cut-stats "$CASE_DIR"

status=0

DROPPED=$(grep -E "^droppedFaces" "$RUN_LOG" | sed -E 's/.*= ([0-9]+).*/\1/' || true)
echo "measured droppedFaces = ${DROPPED:-<missing>}"
if [ -z "${DROPPED:-}" ] || [ "$DROPPED" -le 0 ] 2>/dev/null; then
    echo "FAIL: droppedFaces not > 0 (forced drop did not fire)"
    status=1
fi

WALL=$(grep -E "^wallArea" "$RUN_LOG" | sed -E 's/.*= ([0-9.eE+-]+).*/\1/' || true)
echo "measured wallArea = ${WALL:-<missing>}"
if [ -z "${WALL:-}" ]; then
    echo "FAIL: 'wallArea' line not found in ninjaMesher output"
    status=1
else
    OK=$(awk -v a="$WALL" 'BEGIN { print (a > 0.765 && a < 1.10) ? "1" : "0" }')
    if [ "$OK" != "1" ]; then
        echo "FAIL: wallArea $WALL not strictly between landed (~0.755) and offset (1.131) bands (0.765, 1.10)"
        status=1
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

if [ "$status" -ne 0 ]; then
    cat "$LOG"
fi

exit $status
