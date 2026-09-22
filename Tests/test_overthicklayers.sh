#!/bin/sh
# Gate for layers_overthick: a layer spec the requested cells cannot
# possibly carry must degrade to NO LAYERS, never to a bad cell.
#
# "La malla NUNCA tiene que salir con problemas de calidad. En todo caso
# tiene que simplemente quedar lejos de la pared." See the case dict for
# the spec (10 layers at ratio 5 inside 0.05 on h = 0.1 cells: a first
# layer of 2.0e-8 m under a ~0.1 m face) and for the measured
# before/after.
#
# Asserted:
#   * checkMesh -allGeometry reports NO high-aspect-ratio cells, no open
#     cells, no misoriented face pyramids, no bad face tets and no
#     highly skew faces. The high-aspect check is the point: before the
#     aspect-ratio ceiling in prismValid this case shipped
#     "High aspect ratio cells found, Max aspect ratio: 11375.3,
#     number of cells 480";
#   * max aspect ratio strictly below checkMesh's own threshold of 1000;
#   * check_integrity.py clean;
#   * the march actually DECLINED the wall rather than silently doing
#     nothing: infeasibleWallArea must be a large fraction (>= 90%) of
#     wallArea, i.e. the mesh really did stay at the offset cut. A
#     future change that made the guard reject everything by refusing to
#     run the march at all would still be wrong, so the cut's own wall
#     must still be there: wallArea > 0;
#   * determinism: a second np=1 run is byte-identical.
#
# Usage: test_overthicklayers.sh <ninjaMesher-exe> <caseDir> <stl-args...>
#   The trailing arguments are passed to Tests/scripts/make_stls.py to
#   regenerate the fixture STL, as the other synthetic gates do.
set -eu

EXE="$1"
CASE_DIR="$2"
shift 2

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

python3 "$SCRIPT_DIR/scripts/make_stls.py" "$@" "$CASE_DIR/sphere.stl"

OUT_DIR="$CASE_DIR/constant/polyMesh"
status=0

rm -rf "$OUT_DIR"
mpirun --bind-to none -np 1 "$EXE" "$CASE_DIR" 2>&1 | tee "$TMP/run1.log"
cp -r "$OUT_DIR" "$TMP/mesh_run1"

echo "--- check_integrity.py ---"
if python3 "$SCRIPT_DIR/../Benchmarks/check_integrity.py" "$OUT_DIR"; then
    echo "PASS: integrity"
else
    echo "FAIL: check_integrity.py reported defects"
    status=1
fi

echo "--- checkMesh -allGeometry ---"
( cd "$CASE_DIR" && checkMesh -allGeometry ) > "$TMP/checkmesh.log" 2>&1 || true
grep -E "Total volume|Open cells|face pyramids|face tets|skewness|aspect ratio|Failed .* mesh checks|Mesh OK" \
    "$TMP/checkmesh.log" || true
for pat in "Open cells found" "incorrectly oriented" "Error in face tets" "highly skew faces" "High aspect ratio"; do
    if grep -q "$pat" "$TMP/checkmesh.log"; then
        echo "FAIL: checkMesh reports '$pat'"
        status=1
    fi
done

AR=$(grep -E "Max aspect ratio" "$TMP/checkmesh.log" \
     | sed -E 's/.*Max aspect ratio(:| =) ([0-9.eE+-]+).*/\2/' | head -1 || true)
echo "measured Max aspect ratio = ${AR:-<missing>}"
if [ -z "${AR:-}" ]; then
    echo "FAIL: checkMesh reported no aspect ratio"
    status=1
elif ! python3 -c "import sys; sys.exit(0 if float(sys.argv[1]) < 1000.0 else 1)" "$AR"; then
    echo "FAIL: max aspect ratio $AR is at or above checkMesh's threshold of 1000"
    status=1
fi

WA=$(grep -E "^wallArea = " "$TMP/run1.log" | sed -E 's/.*= ([0-9.eE+-]+).*/\1/' | head -1 || true)
IWA=$(grep -E "^infeasibleWallArea = " "$TMP/run1.log" | sed -E 's/.*= ([0-9.eE+-]+).*/\1/' | head -1 || true)
echo "measured wallArea = ${WA:-<missing>}  infeasibleWallArea = ${IWA:-<missing>}"
if [ -z "${WA:-}" ] || [ -z "${IWA:-}" ]; then
    echo "FAIL: wallArea / infeasibleWallArea missing from the mesher output"
    status=1
elif ! python3 -c "import sys; w=float(sys.argv[1]); i=float(sys.argv[2]); sys.exit(0 if w > 0.0 and i >= 0.9*w else 1)" "$WA" "$IWA"; then
    echo "FAIL: the march did not decline the wall (expected infeasibleWallArea >= 90% of wallArea > 0)"
    status=1
fi

rm -rf "$OUT_DIR"
mpirun --bind-to none -np 1 "$EXE" "$CASE_DIR" >"$TMP/run2.log" 2>&1
cp -r "$OUT_DIR" "$TMP/mesh_run2"
if diff -r "$TMP/mesh_run1" "$TMP/mesh_run2" >/dev/null 2>&1; then
    echo "PASS: determinism (two np=1 runs byte-identical)"
else
    echo "FAIL: two np=1 runs produced DIFFERENT polyMesh output"
    status=1
fi

exit $status
