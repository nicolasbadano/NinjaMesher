#!/bin/sh
# Gate for win_desc_slot: the layer march must never bury its own wall
# face inside the solid.
#
# The window is a grid-aligned sub-box of bm_layers_descargador. BEFORE the
# seal guard the parent shipped 5 sealed regions / 14.977 m2 of landed wall
# faces sitting INSIDE the CAD's thin plates; this window reproduced 90 such
# stacks. The STL is the reconstructed surface described in gates.toml; on
# the raw CAD this window's locationInMesh sat in the overlap of two bodies
# and the window meshed an inverted region (138.5 m3 of it).
#
# Asserted:
#   * smoothSealedRegionCount == 0 and smoothSealedRegionArea == 0 --
#     no landed wall face ends up inside the solid,
#   * the seal guard actually FIRED here (it is a live regression
#     tripwire, not a dormant one; measured 48 stacks),
#   * the layer-quality gate converged (no "stacks still bad" warning),
#   * check_integrity.py clean,
#   * checkMesh -allGeometry: no open cells, misoriented face pyramids,
#     bad face tets, highly skew faces or high aspect ratio cells,
#   * total volume within 1% of the measured 53.744 m3,
#   * determinism: a second np=1 run is byte-identical.
#
# Usage: test_descslot.sh <ninjaMesher-exe> <caseDir>
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

SC=$(grep -E "^smoothSealedRegionCount" "$TMP/run1.log" | sed -E 's/.*= ([0-9]+).*/\1/' || true)
echo "measured smoothSealedRegionCount = ${SC:-<missing>}"
if [ -z "${SC:-}" ]; then
    echo "FAIL: no smoothSealedRegionCount line in the mesher output"
    status=1
elif [ "$SC" -ne 0 ]; then
    echo "FAIL: $SC sealed region(s) -- landed wall faces inside the solid"
    status=1
fi

if ! grep -q "^layer seal guard: " "$TMP/run1.log"; then
    echo "FAIL: the seal guard did not fire (this window exists because it does)"
    status=1
fi
if grep -q "stacks still bad" "$TMP/run1.log"; then
    echo "FAIL: layer-quality gate did not converge"
    status=1
fi

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

# checkMesh prints "Total volume = 138.485." -- the sentence's full stop
# is part of the match, so it is trimmed here.
VOL=$(grep -E "Total volume" "$TMP/checkmesh.log" | sed -E 's/.*Total volume = ([0-9.eE+-]+?)\.? .*/\1/' | head -1 || true)
VOL=${VOL%.}
echo "measured Total volume = ${VOL:-<missing>}"
if [ -z "${VOL:-}" ]; then
    echo "FAIL: checkMesh reported no Total volume"
    status=1
elif ! python3 -c "import sys; v=float(sys.argv[1]); sys.exit(0 if abs(v-53.744)/53.744 <= 0.01 else 1)" "$VOL"; then
    echo "FAIL: total volume $VOL outside 53.744 +- 1%"
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
