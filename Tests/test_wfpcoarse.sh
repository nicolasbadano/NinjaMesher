#!/bin/sh
# Gate for win_wfp_coarse_cut / win_wfp_coarse_layers: a grid 4x too coarse
# must never ship a highly skew BOUNDARY face.
#
# Both windows are grid-aligned sub-boxes of bm_layers_wfp at grid factor 4
# (dx 0.8 -> 3.2 m). checkMesh scores a boundary face against its owner's
# centre alone, and three stages create or re-centre boundary faces without
# having asked that question:
#   * cutMesh -- a clipped grid face ships as boundary once the across-cell
#     is removed (_cut: 1 face at skewness 15.5828 before the post-pass);
#   * mergeSlivers -- a merged group keeps every member's wall faces but
#     moves the centre (_layers: 1 face at 4.28628 before the group gate
#     scored boundary faces);
#   * planarize -- a fan triangle at the far end of a large wall polygon
#     (the parent case at factor 4, layers on: 4.23394).
#
# Asserted:
#   * checkMesh -allGeometry reports no highly skew faces, no open cells, no
#     misoriented face pyramids, no bad face tets;
#   * check_integrity.py clean;
#   * total volume within 2% of the measured value, so the guard cannot be
#     "fixed" in future by deleting mesh;
#   * determinism: a second np=1 run is byte-identical.
#
# Usage: test_wfpcoarse.sh <ninjaMesher-exe> <caseDir> <expectedVolume>
set -eu

EXE="$1"
CASE_DIR="$2"
EXPECTED_VOL="$3"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

if [ ! -e "$CASE_DIR/wfp_fw.stl" ] || [ ! -e "$CASE_DIR/wfp_ib.stl" ]; then
    echo "FAIL: $CASE_DIR STLs missing (relative symlinks into Benchmarks/cases/bm_layers_wfp)"
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
grep -E "Total volume|Open cells|face pyramids|face tets|skewness|Failed .* mesh checks|Mesh OK" \
    "$TMP/checkmesh.log" || true
if ! grep -q "Max skewness" "$TMP/checkmesh.log"; then
    echo "FAIL: checkMesh reported no skewness"
    status=1
fi
for pat in "Open cells found" "incorrectly oriented" "Error in face tets" "highly skew faces"; do
    if grep -q "$pat" "$TMP/checkmesh.log"; then
        echo "FAIL: checkMesh reports '$pat'"
        status=1
    fi
done

# checkMesh prints "Total volume = 185.422." -- the sentence's full stop
# is part of the match, so it is trimmed here.
VOL=$(grep -E "Total volume" "$TMP/checkmesh.log" | sed -E 's/.*Total volume = ([0-9.eE+-]+?)\.? .*/\1/' | head -1 || true)
VOL=${VOL%.}
echo "measured Total volume = ${VOL:-<missing>}"
if [ -z "${VOL:-}" ]; then
    echo "FAIL: checkMesh reported no Total volume"
    status=1
elif ! python3 -c "import sys; v=float(sys.argv[1]); e=float(sys.argv[2]); sys.exit(0 if abs(v-e)/e <= 0.02 else 1)" "$VOL" "$EXPECTED_VOL"; then
    echo "FAIL: total volume $VOL outside $EXPECTED_VOL +- 2%"
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
