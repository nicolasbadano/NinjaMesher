#!/bin/sh
# Gate for offset_overlap: two OVERLAPPING icospheres
# (r=0.2, centres 0.2 apart -- 2r=0.4 > 0.2, interiors intersect;
# lens volume analytic 0.010472), sphere1 layered (t=0.03), sphere2 a
# plain wall. Asserts:
#   - Mesh OK. + integrity clean
#   - total volume within a band of the direct-cut UNION complement
#     (1 - (V1+V2-Vlens) = 1 - 0.056549 = 0.943451, 5% band)
#   - both sphere1/sphere2 patches present with a nonzero, sane face
#     count (the winning-intercept attribution rule; a "buried
#     surface" bug would either crash the cutter or leave one patch
#     empty/oversized)
#
# Usage: test_offsetoverlap.sh <ninjaMesher-exe> <caseDir>
set -eu

EXE="$1"
CASE_DIR="$2"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

python3 "$SCRIPT_DIR/scripts/make_stls.py" sphere_at 0.4 0.5 0.5 0.2 3 "$CASE_DIR/sphere1.stl"
python3 "$SCRIPT_DIR/scripts/make_stls.py" sphere_at 0.6 0.5 0.5 0.2 3 "$CASE_DIR/sphere2.stl"

OUT_DIR="$CASE_DIR/constant/polyMesh"
rm -rf "$OUT_DIR"
mpirun --bind-to none -np "${NINJA_NP:-1}" "$EXE" "$CASE_DIR"

echo "--- --cut-stats ---"
mpirun --bind-to none -np "${NINJA_NP:-1}" "$EXE" --cut-stats "$CASE_DIR"

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

if ! python3 "$SCRIPT_DIR/../Benchmarks/check_integrity.py" "$OUT_DIR"; then
    echo "FAIL: check_integrity.py reported a structural pathology"
    status=1
fi

EXPECT_VOL="0.943451"
TOL="0.047173" # 5%
VOL_LINE=$(grep -E "Total volume" "$LOG" || true)
if [ -z "$VOL_LINE" ]; then
    echo "FAIL: 'Total volume' line not found in checkMesh output"
    status=1
else
    ACTUAL_VOL=$(echo "$VOL_LINE" | sed -E 's/.*Total volume = ([0-9.eE+-]+).*/\1/')
    echo "measured volume = $ACTUAL_VOL"
    OK=$(awk -v a="$ACTUAL_VOL" -v e="$EXPECT_VOL" -v t="$TOL" \
        'BEGIN { d = a - e; if (d < 0) d = -d; print (d <= t) ? "1" : "0" }')
    if [ "$OK" != "1" ]; then
        echo "FAIL: total volume $ACTUAL_VOL not within $TOL of expected $EXPECT_VOL"
        status=1
    fi
fi

AREAS=$(python3 "$SCRIPT_DIR/scripts/patch_areas.py" "$OUT_DIR" sphere1 sphere2) || {
    echo "FAIL: patch_areas.py could not find both sphere patches"
    exit 1
}
echo "$AREAS"
A1=$(echo "$AREAS" | awk '$1=="sphere1"{print $2}')
A2=$(echo "$AREAS" | awk '$1=="sphere2"{print $2}')
for pair in "sphere1 $A1" "sphere2 $A2"; do
    set -- $pair
    OK=$(awk -v a="$2" 'BEGIN { print (a > 0.0) ? "1" : "0" }')
    if [ "$OK" != "1" ]; then
        echo "FAIL: $1 patch has zero area (buried surface or attribution failure)"
        status=1
    fi
done

if [ "$status" -ne 0 ]; then
    cat "$LOG"
fi

exit $status
