#!/bin/sh
# Gate for offset_twospheres_layers: two congruent icospheres
# (r=0.15 at x=0.25 / x=0.75) with DIFFERENT layer thicknesses (t=0.03 /
# t=0.06) -- the multi-solid attribution probe for the offset cut +
# layer-by-layer march, and the regression gate for a two-sphere
# orientation bug (16 misoriented faces + sphere1's entire
# wall patch drop-reverted).
#
# Asserts: "Mesh OK.", integrity, total volume in a band around the
# direct-cut two-sphere complement, per-patch landed wall areas inside
# a SANITY band (bands are sanity checks, not
# fidelity assertions -- surface accuracy is the mission's documented
# trade), and the two patches' areas equal within 5% (congruent
# spheres land on the same true surface; the DIFFERENT t gives the two
# fronts legitimately different tessellations -- measured 52 vs 92
# wall faces, 2.3% area difference -- while a reverted-patch
# regression reads ~+50%, far outside).
#
# Usage: test_offsettwospheres.sh <ninjaMesher-exe> <caseDir>
set -eu

EXE="$1"
CASE_DIR="$2"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

python3 "$SCRIPT_DIR/scripts/make_stls.py" sphere_at 0.25 0.5 0.5 0.15 3 "$CASE_DIR/sphere1.stl"
python3 "$SCRIPT_DIR/scripts/make_stls.py" sphere_at 0.75 0.5 0.5 0.15 3 "$CASE_DIR/sphere2.stl"

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

# Total volume: landed walls approximate the direct two-sphere cut
# (complement 0.971726, same reference as the two_spheres gate), within
# 2%.
EXPECT_VOL="0.971726"
TOL="0.019435"
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

# Per-patch landed areas. Sanity band: [0.80, 1.05] x the smooth
# analytic 4*pi*0.15^2 = 0.282743 (the landed faceted wall
# under-samples chords at this resolution -- measured 0.2464, ~87% of
# smooth; a reverted-to-offset patch would read far outside the band,
# which is exactly the failure mode this band exists to
# catch).
AREAS=$(python3 "$SCRIPT_DIR/scripts/patch_areas.py" "$OUT_DIR" sphere1 sphere2) || {
    echo "FAIL: patch_areas.py could not find both sphere patches"
    exit 1
}
echo "$AREAS"
A1=$(echo "$AREAS" | awk '$1=="sphere1"{print $2}')
A2=$(echo "$AREAS" | awk '$1=="sphere2"{print $2}')
ANALYTIC="0.282743"
for pair in "sphere1 $A1" "sphere2 $A2"; do
    set -- $pair
    OK=$(awk -v a="$2" -v ref="$ANALYTIC" \
        'BEGIN { print (a >= 0.80 * ref && a <= 1.05 * ref) ? "1" : "0" }')
    if [ "$OK" != "1" ]; then
        echo "FAIL: $1 landed area $2 outside sanity band [0.80, 1.05] x $ANALYTIC"
        status=1
    fi
done
OK=$(awk -v a="$A1" -v b="$A2" 'BEGIN {
    m = (a > b) ? a : b; d = a - b; if (d < 0) d = -d;
    print (m > 0 && d / m <= 0.05) ? "1" : "0" }')
if [ "$OK" != "1" ]; then
    echo "FAIL: per-patch areas differ by more than 5% (sphere1=$A1 sphere2=$A2) -- attribution asymmetry"
    status=1
fi

if [ "$status" -ne 0 ]; then
    cat "$LOG"
fi

exit $status
