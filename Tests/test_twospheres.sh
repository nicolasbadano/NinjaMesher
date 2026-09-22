#!/bin/sh
# Gate: two disjoint icospheres, each its own wall patch.
# Asserts "Mesh OK.", both wall patches (sphere1/sphere2) present with
# nonzero nFaces, and total kept volume within 2% of the analytic
# two-sphere-complement volume.
#
# Usage: test_twospheres.sh <ninjaMesher-exe> <caseDir>
set -eu

EXE="$1"
CASE_DIR="$2"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

python3 "$SCRIPT_DIR/scripts/make_stls.py" sphere_at 0.25 0.5 0.5 0.15 3 "$CASE_DIR/sphere1.stl"
python3 "$SCRIPT_DIR/scripts/make_stls.py" sphere_at 0.75 0.5 0.5 0.15 3 "$CASE_DIR/sphere2.stl"

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

BOUNDARY_FILE="$OUT_DIR/boundary"
for patch in sphere1 sphere2; do
    nfaces=$(awk -v p="$patch" '
        $0 ~ "^[[:space:]]*"p"$" { infound=1; next }
        infound && /nFaces/ { gsub(/[^0-9]/, "", $2); print $2; exit }
    ' "$BOUNDARY_FILE")
    if [ -z "${nfaces:-}" ] || [ "$nfaces" -le 0 ] 2>/dev/null; then
        echo "FAIL: patch '$patch' missing or has zero faces (nFaces='${nfaces:-}')"
        status=1
    fi
done

EXPECT_VOL="0.971726"
TOL="0.019435"
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
