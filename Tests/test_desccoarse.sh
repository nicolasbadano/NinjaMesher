#!/bin/sh
# Gate for win_desc_coarse: a grid 8x too coarse for the layer stack the
# dict asks for must cost LAYERS, never validity.
#
# "La malla NUNCA tiene que salir con problemas de calidad. En todo caso
# tiene que simplemente quedar lejos de la pared." The window is a
# grid-aligned 4x2x3-cell sub-box of bm_layers_descargador at grid factor
# 8 (dx 0.5 -> 4.0 m) carrying the parent's layer block verbatim; see the
# case dict for the full derivation and the measured before/after.
#
# Asserted:
#   * checkMesh -allGeometry reports NO high-aspect-ratio cells, no open
#     cells, no misoriented face pyramids, no bad face tets, no highly
#     skew faces. The high-aspect check is the point: on the binary
#     before the aspect-ratio ceiling in prismValid this window shipped
#     "High aspect ratio cells found, Max aspect ratio: 1945.65,
#     number of cells 4";
#   * max aspect ratio strictly below checkMesh's own threshold of 1000;
#   * check_integrity.py clean;
#   * the march really did stay AWAY from the wall rather than silently
#     skipping the stage: the cut's wall is still there (wallArea > 0)
#     and infeasibleWallArea covers at least a quarter of it (measured
#     55% on the reconstructed STL -- see gates.toml; 34% on the raw CAD
#     counted on the fluid that is written, 45% while the wall faces of
#     discarded components were counted too, 50+% while the terrace-seam
#     cascade still dropped the healthy stacks beside every genuine
#     failure);
#   * total volume within 2% of the measured 393.974 m3, so the guard
#     cannot be "fixed" in future by deleting mesh (364.808 before the
#     cascade fix, 382.197 while a failing face also dropped its one-ring
#     neighbours: the layers it now keeps are volume it used to lose;
#     411.537 with layers off);
#   * determinism: a second np=1 run is byte-identical.
#
# Usage: test_desccoarse.sh <ninjaMesher-exe> <caseDir>
set -eu

EXE="$1"
CASE_DIR="$2"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

if [ ! -e "$CASE_DIR/descargadorG4bT.stl" ] || [ ! -e "$CASE_DIR/aireadorT.stl" ]; then
    echo "FAIL: $CASE_DIR STLs missing (relative symlinks into Benchmarks/cases/bm_layers_descargador)"
    exit 1
fi

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
elif ! python3 -c "import sys; w=float(sys.argv[1]); i=float(sys.argv[2]); sys.exit(0 if w > 0.0 and i >= 0.25*w else 1)" "$WA" "$IWA"; then
    echo "FAIL: the march did not stay away from the wall (expected infeasibleWallArea >= 25% of wallArea > 0)"
    status=1
fi

# checkMesh prints "Total volume = 393.974." -- the sentence's full stop
# is part of the match, so it is trimmed here.
VOL=$(grep -E "Total volume" "$TMP/checkmesh.log" | sed -E 's/.*Total volume = ([0-9.eE+-]+?)\.? .*/\1/' | head -1 || true)
VOL=${VOL%.}
echo "measured Total volume = ${VOL:-<missing>}"
if [ -z "${VOL:-}" ]; then
    echo "FAIL: checkMesh reported no Total volume"
    status=1
elif ! python3 -c "import sys; v=float(sys.argv[1]); sys.exit(0 if abs(v-393.974)/393.974 <= 0.02 else 1)" "$VOL"; then
    echo "FAIL: total volume $VOL outside 393.974 +- 2%"
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
